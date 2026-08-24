// wax.js — the wax phase: a position-based fluid (Macklin & Muller 2013) with
// interfacial tension (Akinci et al. 2013), per-particle temperature, latent
// heat of fusion, and blob-scale coupling to the surrounding aqueous column.
//
// The model is a DISPERSED phase in a CONTINUOUS phase. The wax is resolved as
// particles; the water/glycol around it is the 1-D column in column.js. Three
// couplings run between them, and each is taken from the continuum result for a
// drop of the blob's actual size rather than invented per particle:
//
//   buoyancy  per particle, from its own temperature against the fluid
//             temperature at its own height. This is the engine of the lamp.
//   drag      per BLOB, quadratic, Cd ~ 1 for a wobbling deformable drop, using
//             the blob's equivalent radius. Applying drag per particle instead
//             would make a 3 cm blob rise at the same speed as a 2 mm droplet,
//             which is the one thing a lava lamp visibly does not do.
//   heat      per BLOB surface area, with the film coefficient from
//             Nu = 2 + 0.6 Re^0.5 Pr^(1/3), distributed over the blob surface
//             particles so the total flux equals the continuum flux exactly.
//
// Blobs are found by union-find over the neighbour graph every step, so "blob"
// means a connected lump of wax, not a label anyone assigned.
//
// The inner loops are written flat on purpose: local aliases for every typed
// array, squared distances until a square root is actually needed, kernels
// inlined, and no Math.hypot anywhere it runs per pair.

import {
  GEO, WAX, AQ, HEAT, IFACE, SOLVER, rhoWax, rhoAq, globeRadius,
} from './params.js';

const MAXN = 64;                 // neighbour cap per particle

export class Wax {
  constructor(n = SOLVER.particles, volume = WAX.volume) {
    this.alloc(n);
    this.setVolume(volume);
    this.blobCount = 0;
    this.biggestBlob = 0;
    this.stirTimer = 0;
    this.plateQ = 0;
    this.clamps = 0;
  }

  alloc(n) {
    this.n = n;
    const f = () => new Float32Array(n);
    this.x = f(); this.y = f(); this.z = f();
    this.vx = f(); this.vy = f(); this.vz = f();
    this.px = f(); this.py = f(); this.pz = f();
    this.T = f(); this.rho = f(); this.lam = f();
    this.nx = f(); this.ny = f(); this.nz = f();
    this.solid = f();
    this.blob = new Int32Array(n);        // union-find parent, then blob root
    this.nbr = new Int32Array(n * MAXN);
    this.nbrCount = new Int32Array(n);
    this.blobA = new Float32Array(n);     // equivalent radius of my blob
    this.blobVx = new Float32Array(n);
    this.blobVy = new Float32Array(n);
    this.blobVz = new Float32Array(n);
    this.expose = new Float32Array(n);    // surface-exposure weight, m^2
    this.id = new Int32Array(n);          // PERSISTENT blob identity
    this.nextId = 2;
    this.pairT = new Map();               // "a:b" -> seconds in contact
    this.merged = new Set();              // pairs whose film has drained
    this.contact = new Map();             // pairs touching this step
    this.lastDt = 1 / 30;
    this.dpx = f(); this.dpy = f(); this.dpz = f();
    this.dT = f();
    this._bc = new Int32Array(n);
    this._bvx = f(); this._bvy = f(); this._bvz = f(); this._bex = f();
    // globe radius lookup, so confine() is a table read not a profile walk
    this.profN = 512;
    this.profR = new Float32Array(this.profN + 1);
    for (let i = 0; i <= this.profN; i++) this.profR[i] = globeRadius(GEO.H * i / this.profN);
  }

  // Particle mass and smoothing length follow from the wax volume and count.
  setVolume(volume) {
    this.volume = volume;
    this.Vp = volume / this.n;
    this.rp = Math.cbrt(3 * this.Vp / (4 * Math.PI));
    this.dx = Math.cbrt(this.Vp);                  // rest spacing
    this.h = SOLVER.hFactor * this.dx;
    this.m = WAX.rho20 * this.Vp;
    const h = this.h;
    this.h2 = h * h;
    this.kPoly = 315 / (64 * Math.PI * h ** 9);
    this.kSpiky = -45 / (Math.PI * h ** 6);
    this.kVisc = 45 / (Math.PI * h ** 6);
    this.kCoh = 32 / (Math.PI * h ** 9);
    this.rho0 = this.latticeDensity();
    this.wDq = this.poly6(0.2 * h);                // tensile-instability ref
    // CFM regularisation must sit on the same scale as sum|grad C|^2, which
    // depends on h and the particle mass. A fixed epsilon would be a no-op for
    // one particle count and a straitjacket for another.
    const gref = (this.m / this.rho0) * Math.abs(this.kSpiky) * 0.25 * h * h;
    this.lamScale = 1 / (30 * gref * gref);        // natural magnitude of lambda
    this.eps = SOLVER.relaxEps / this.lamScale;
    // The artificial-pressure term of Macklin and Muller is quoted as a bare
    // constant because their lambda is O(0.1). Here lambda is O(1e-6): the same
    // bare constant is a thousand times the correction it is supposed to nudge,
    // and it detonates the fluid on the first projection. It has to be scaled to
    // lambda, not copied from the paper.
    this.scorr = SOLVER.scorr * this.lamScale;
    // Interfacial-tension coefficients, dimensionally closed so that changing
    // the particle count or sigma does not change the physics:
    //   cohesion   a = cohA * m * C(r),   [cohA] = L^4 / (M T^2)
    //   curvature  a = curvA * dn,        [curvA] = L / T^2
    this.cohScale = 1 / (this.rho0 * this.rho0 * this.dx * this.dx);
    // The curvature term is tied to the cohesion term through the cohesion
    // spline's own characteristic magnitude, C(0) = kCoh*h^6/64. Akinci et al.
    // use one gamma for both, but their two expressions do not carry the same
    // units, so a single constant cannot scale both: written literally the
    // curvature term came out some six hundred times the cohesion it is meant
    // to trim, and blew a free drop apart in a tenth of a second. Anchoring it
    // to C(0) leaves curvK as a pure shape parameter of order one.
    this.cohRef = this.kCoh * h ** 6 / 64;
    this.vMax = 0.35 * h;                          // per-second-of-dt CFL bound
    // Largest stable explicit step for conduction WITHIN the wax. This is a
    // separate limit from the mechanical one and it tightens as dx shrinks:
    // rescaling the globe cut the spacing from 4.5 mm to 3.5 mm, the limit fell
    // from 20 s to 12 s, and the 20-second relaxation step the warm start had
    // been using quietly went unstable and reported a base plate at -861 C.
    // The rate is the diagonal of the discrete operator, summed over a rest
    // lattice, so it follows the parameters instead of being guessed.
    let sw = 0;
    const kk = Math.ceil(h / this.dx);
    for (let i = -kk; i <= kk; i++)
      for (let j = -kk; j <= kk; j++)
        for (let l = -kk; l <= kk; l++) {
          const r = this.dx * Math.sqrt(i * i + j * j + l * l);
          if (r > 1e-12 && r < h) sw += h - r;
        }
    const rate = (WAX.k * this.Vp * this.kVisc) * (this.m / this.rho0) * sw /
                 (this.m * WAX.c);
    this.heatDtMax = 0.4 / Math.max(rate, 1e-9);
  }

  poly6(r) { const h = this.h; return r >= h ? 0 : this.kPoly * (h * h - r * r) ** 3; }

  // Rest density measured from a regular lattice at the rest spacing rather
  // than assumed equal to the bulk density. This is what makes the density
  // constraint mean "back to rest spacing".
  latticeDensity() {
    let s = 0;
    const k = Math.ceil(this.h / this.dx);
    for (let i = -k; i <= k; i++)
      for (let j = -k; j <= k; j++)
        for (let l = -k; l <= k; l++)
          s += this.poly6(this.dx * Math.sqrt(i * i + j * j + l * l));
    return this.m * s;
  }

  // Largest stable step for explicit surface tension: the capillary-wave CFL
  //
  //     dt <= 0.25 * sqrt( rho h^3 / (2 pi sigma) )
  //
  // The shortest capillary wave the kernel can carry has wavelength ~h, and an
  // explicit integrator has to resolve its period. Running past this does not
  // look like divergence -- it looks like a velocity clamp firing on two thirds
  // of the particles every step, which is what 1/30 s did here: 1.77x over the
  // limit, 660,000 clamp events in twelve seconds of lamp time, and a lamp that
  // otherwise seemed merely sluggish.
  //
  // Note what it costs. The limit scales as h^1.5, h goes as (V/n)^(1/3), so
  // dt ~ n^(-1/2) and the compute per second of lamp time goes as n^1.5. That is
  // the real ceiling on a particle lava lamp in a browser, and it is why the
  // default sits at 1800 particles rather than as many as memory allows.
  capillaryDt(sigma = IFACE.sigma) {
    return 0.25 * Math.sqrt(this.rho0 * this.h ** 3 / (2 * Math.PI * Math.max(sigma, 1e-5)));
  }

  radiusAt(y) {
    let t = y / GEO.H;
    if (t < 0) t = 0; else if (t > 1) t = 1;
    const f = t * this.profN;
    const i = f | 0;
    const j = i < this.profN ? i + 1 : this.profN;
    const a = this.profR[i];
    return a + (this.profR[j] - a) * (f - i);
  }

  // -------------------------------------------------------------------------
  // Cold start: a solid plug of wax resting on the base plate.
  //
  // Seeded on a lattice at the rest spacing, not scattered at random. A random
  // cloud starts at the wrong density everywhere, and the first density
  // projection then has to shove every particle several millimetres, which reads
  // out as metres-per-second of velocity. The lattice starts the solver at the
  // configuration it is trying to reach.
  // -------------------------------------------------------------------------
  reset(Tamb = HEAT.Tamb, rand = Math.random) {
    const dx = this.dx, n = this.n;
    let placed = 0, y = this.rp;
    const jit = 0.06 * dx;
    while (placed < n && y < GEO.H) {
      const R = this.radiusAt(y) - 1.3 * this.rp;
      const cols = Math.max(1, Math.floor(2 * R / dx));
      const off = -0.5 * (cols - 1) * dx;
      for (let a = 0; a < cols && placed < n; a++) {
        for (let b = 0; b < cols && placed < n; b++) {
          const cx = off + a * dx, cz = off + b * dx;
          if (cx * cx + cz * cz > R * R) continue;
          this.x[placed] = cx + (rand() * 2 - 1) * jit;
          this.z[placed] = cz + (rand() * 2 - 1) * jit;
          this.y[placed] = y + (rand() * 2 - 1) * jit;
          placed++;
        }
      }
      y += dx;
    }
    for (let i = 0; i < n; i++) {
      this.vx[i] = this.vy[i] = this.vz[i] = 0;
      this.T[i] = Tamb;
      this.solid[i] = 1;
      this.blobA[i] = this.rp;
      this.blobVx[i] = this.blobVy[i] = this.blobVz[i] = 0;
      this.rho[i] = this.rho0;
      this.px[i] = this.x[i]; this.py[i] = this.y[i]; this.pz[i] = this.z[i];
      this.confine(i);
      this.x[i] = this.px[i]; this.y[i] = this.py[i]; this.z[i] = this.pz[i];
    }
    this.fillHeight = y;
    this.stirTimer = 0;
    this.clamps = 0;
    this.id.fill(1);
    this.nextId = 2;
    this.pairT.clear(); this.merged.clear(); this.contact.clear();
    this.buildGrid();
    this.buildNeighbours();
    this.findBlobs();
  }

  // -------------------------------------------------------------------------
  // Neighbour search: counting sort into a uniform grid of cell size h.
  // -------------------------------------------------------------------------
  buildGrid() {
    const h = this.h, n = this.n, px = this.px, py = this.py, pz = this.pz;
    let x0 = Infinity, y0 = Infinity, z0 = Infinity;
    let x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
    for (let i = 0; i < n; i++) {
      const a = px[i], b = py[i], c = pz[i];
      if (a < x0) x0 = a; if (a > x1) x1 = a;
      if (b < y0) y0 = b; if (b > y1) y1 = b;
      if (c < z0) z0 = c; if (c > z1) z1 = c;
    }
    const nx = Math.max(1, Math.ceil((x1 - x0) / h) + 1);
    const ny = Math.max(1, Math.ceil((y1 - y0) / h) + 1);
    const nz = Math.max(1, Math.ceil((z1 - z0) / h) + 1);
    const nc = nx * ny * nz;
    if (!this.cellStart || this.cellStart.length < nc + 1) {
      this.cellStart = new Int32Array(nc + 1);
      this.cellCount = new Int32Array(nc + 1);
      this.order = new Int32Array(n);
      this.cellIdx = new Int32Array(n);
    }
    this.gnx = nx; this.gny = ny; this.gnz = nz;
    const cellCount = this.cellCount, cellStart = this.cellStart;
    const cellIdx = this.cellIdx, order = this.order;
    cellCount.fill(0, 0, nc);
    const ih = 1 / h;
    for (let i = 0; i < n; i++) {
      let cx = ((px[i] - x0) * ih) | 0; if (cx >= nx) cx = nx - 1;
      let cy = ((py[i] - y0) * ih) | 0; if (cy >= ny) cy = ny - 1;
      let cz = ((pz[i] - z0) * ih) | 0; if (cz >= nz) cz = nz - 1;
      const c = (cz * ny + cy) * nx + cx;
      cellIdx[i] = c;
      cellCount[c]++;
    }
    let acc = 0;
    for (let c = 0; c < nc; c++) { cellStart[c] = acc; acc += cellCount[c]; }
    cellStart[nc] = acc;
    for (let c = 0; c < nc; c++) cellCount[c] = cellStart[c];   // reuse as cursor
    for (let i = 0; i < n; i++) order[cellCount[cellIdx[i]]++] = i;
  }

  buildNeighbours() {
    const nx = this.gnx, ny = this.gny, nz = this.gnz, h2 = this.h2, n = this.n;
    const px = this.px, py = this.py, pz = this.pz;
    const nbr = this.nbr, nbrCount = this.nbrCount;
    const cellIdx = this.cellIdx, cellStart = this.cellStart, order = this.order;
    const nxy = nx * ny;
    for (let i = 0; i < n; i++) {
      const c = cellIdx[i];
      const cx = c % nx, cy = ((c / nx) | 0) % ny, cz = (c / nxy) | 0;
      let cnt = 0;
      const base = i * MAXN;
      const xi = px[i], yi = py[i], zi = pz[i];
      const z0 = cz > 0 ? cz - 1 : 0, z1 = cz < nz - 1 ? cz + 1 : nz - 1;
      const y0 = cy > 0 ? cy - 1 : 0, y1 = cy < ny - 1 ? cy + 1 : ny - 1;
      const x0 = cx > 0 ? cx - 1 : 0, x1 = cx < nx - 1 ? cx + 1 : nx - 1;
      for (let gz = z0; gz <= z1; gz++) {
        for (let gy = y0; gy <= y1; gy++) {
          const row = (gz * ny + gy) * nx;
          const s = cellStart[row + x0], e = cellStart[row + x1 + 1];
          for (let k = s; k < e; k++) {
            const j = order[k];
            if (j === i) continue;
            const dx = xi - px[j], dy = yi - py[j], dz = zi - pz[j];
            if (dx * dx + dy * dy + dz * dz < h2) {
              if (cnt >= MAXN) break;
              nbr[base + cnt++] = j;
            }
          }
        }
      }
      nbrCount[i] = cnt;
    }
  }

  // -------------------------------------------------------------------------
  // Blobs: connected components of the neighbour graph, then per-blob volume,
  // equivalent radius, mean velocity, and surface-exposure normalisation.
  // -------------------------------------------------------------------------
  findBlobs() {
    const n = this.n, p = this.blob;
    const px = this.px, py = this.py, pz = this.pz;
    const nbr = this.nbr, nbrCount = this.nbrCount;
    for (let i = 0; i < n; i++) p[i] = i;
    // "Touching" has to mean touching at SPH spacing, not at the ideal lattice
    // spacing. At 1.3 dx a settled drop of 1400 particles reported 230 blobs,
    // because surface particles sit a little further apart than the interior --
    // and every one of those phantom singletons then got particle-scale drag and
    // particle-scale heat transfer instead of its blob's.
    const bond2 = (1.75 * this.dx) ** 2;
    // Two particles are bonded when they are touching AND they belong to the
    // same blob, or to two blobs whose film has already drained. Touching alone
    // is not enough: a surfactant film keeps distinct blobs distinct, and if the
    // component search ignored that, blobs would merge geometrically before the
    // drainage clock had any say in it.
    const idA = this.id;
    const contact = this.contact;
    contact.clear();
    for (let i = 0; i < n; i++) {
      const base = i * MAXN, cnt = nbrCount[i];
      const xi = px[i], yi = py[i], zi = pz[i], idi = idA[i];
      for (let k = 0; k < cnt; k++) {
        const j = nbr[base + k];
        if (j < i) continue;
        const dx = xi - px[j], dy = yi - py[j], dz = zi - pz[j];
        if (dx * dx + dy * dy + dz * dz > bond2) continue;
        const idj = idA[j];
        if (idi !== idj) {
          const key = idi < idj ? idi + ':' + idj : idj + ':' + idi;
          contact.set(key, (contact.get(key) || 0) + 1);
          if (!this.merged.has(key)) continue;
        }
        let ra = i; while (p[ra] !== ra) { p[ra] = p[p[ra]]; ra = p[ra]; }
        let rb = j; while (p[rb] !== rb) { p[rb] = p[p[rb]]; rb = p[rb]; }
        if (ra !== rb) p[ra] = rb;
      }
    }
    for (let i = 0; i < n; i++) {
      let r = i; while (p[r] !== r) { p[r] = p[p[r]]; r = p[r]; }
      p[i] = r;
    }

    const cnt = this._bc, sVx = this._bvx, sVy = this._bvy, sVz = this._bvz, sEx = this._bex;
    cnt.fill(0); sVx.fill(0); sVy.fill(0); sVz.fill(0); sEx.fill(0);
    const irho0 = 1 / this.rho0;
    for (let i = 0; i < n; i++) {
      const r = p[i];
      cnt[r]++;
      sVx[r] += this.vx[i]; sVy[r] += this.vy[i]; sVz[r] += this.vz[i];
      let e = 1 - this.rho[i] * irho0;
      if (e < 0) e = 0; else if (e > 1) e = 1;
      this.expose[i] = e;
      sEx[r] += e;
    }
    let blobs = 0, biggest = 0;
    for (let i = 0; i < n; i++) if (p[i] === i) { blobs++; if (cnt[i] > biggest) biggest = cnt[i]; }
    this.blobCount = blobs;
    this.biggestBlob = biggest;
    this.relabel(p, cnt);
    const c43 = 3 / (4 * Math.PI), Vp = this.Vp;
    for (let i = 0; i < n; i++) {
      const r = p[i], c = cnt[r];
      const a = Math.cbrt(c43 * c * Vp);
      this.blobA[i] = a;
      const ic = 1 / c;
      this.blobVx[i] = sVx[r] * ic; this.blobVy[i] = sVy[r] * ic; this.blobVz[i] = sVz[r] * ic;
      const Ab = 4 * Math.PI * a * a;
      this.expose[i] = sEx[r] > 1e-9 ? this.expose[i] * Ab / sEx[r] : Ab * ic;
    }
  }

  // Carry blob identities across the step.
  //
  // Identity, not position in a list: a blob referenced by its index in this
  // frame's component array is a different blob next frame, and the drainage
  // clock would be timing nothing. Each component inherits the id most of its
  // members already carry; when one id turns up in two components the larger
  // keeps it and the rest are freshly named, which is how a pinch-off leaves a
  // blob whose film has not yet drained against its parent.
  relabel(p, cnt) {
    const n = this.n, id = this.id;
    const tally = new Map();                       // root -> Map(id -> count)
    for (let i = 0; i < n; i++) {
      const r = p[i];
      let t = tally.get(r);
      if (!t) tally.set(r, t = new Map());
      t.set(id[i], (t.get(id[i]) || 0) + 1);
    }
    const claim = new Map();                       // id -> {root, count}
    const major = new Map();                       // root -> id
    for (const [r, t] of tally) {
      let bid = 0, bc = -1;
      for (const [k, c] of t) if (c > bc) { bc = c; bid = k; }
      major.set(r, bid);
      const held = claim.get(bid);
      if (!held || bc > held.count) claim.set(bid, { root: r, count: bc });
    }
    const rootId = new Map();
    for (const [r, bid] of major) {
      rootId.set(r, claim.get(bid).root === r ? bid : this.nextId++);
    }
    for (let i = 0; i < n; i++) id[i] = rootId.get(p[i]);

    // Drainage clocks. A pair in contact accumulates time, a pair that separates
    // forgets it, and a pair past drainTime is allowed to merge from then on.
    const live = new Set();
    for (const [key, touching] of this.contact) {
      if (touching < 3) continue;                  // a grazing particle is not contact
      live.add(key);
      const t = (this.pairT.get(key) || 0) + this.lastDt;
      this.pairT.set(key, t);
      if (t >= SOLVER.drainTime) this.merged.add(key);
    }
    for (const key of [...this.pairT.keys()]) {
      if (!live.has(key)) { this.pairT.delete(key); this.merged.delete(key); }
    }
  }

  // Largest blobs, for the readout.
  blobStats(limit = 8) {
    const seen = new Map();
    for (let i = 0; i < this.n; i++) {
      const r = this.blob[i];
      let e = seen.get(r);
      if (!e) seen.set(r, e = { count: 0, radius: this.blobA[i], vy: this.blobVy[i], T: 0 });
      e.count++; e.T += this.T[i];
    }
    const out = [...seen.values()];
    for (const e of out) e.T /= e.count;
    return out.sort((a, b) => b.count - a.count).slice(0, limit);
  }

  // -------------------------------------------------------------------------
  // One mechanical + thermal step.
  // -------------------------------------------------------------------------
  step(dt, column, env = {}) {
    this.lastDt = dt;
    const n = this.n, h = this.h, h2 = this.h2, m = this.m, rho0 = this.rho0;
    const g = HEAT.g * (env.gravity ?? 1);
    const sigma = env.sigma ?? IFACE.sigma;
    const muAq = env.muAq ?? AQ.mu;
    const iters = env.iterations ?? SOLVER.iterations;
    const x = this.x, y = this.y, z = this.z;
    const vx = this.vx, vy = this.vy, vz = this.vz;
    const px = this.px, py = this.py, pz = this.pz;
    const rho = this.rho, lam = this.lam, solid = this.solid, T = this.T;
    const nbr = this.nbr, nbrCount = this.nbrCount;
    const dpx = this.dpx, dpy = this.dpy, dpz = this.dpz;
    const kPoly = this.kPoly, kSpiky = this.kSpiky, kVisc = this.kVisc;
    const drho = env.dRho ?? (WAX.rho20 - AQ.rho20);
    const rw20 = AQ.rho20 + drho;

    // The aqueous phase has to go somewhere. Mass conservation in a sealed
    // globe fixes its vertical velocity from the wax flux; see
    // Column.updateReturnFlow.
    const hasFlow = typeof column.updateReturnFlow === 'function';
    if (hasFlow) column.updateReturnFlow(this);

    // ---- 1. forces -> tentative velocity --------------------------------
    // Viscous resistance of the aqueous phase to motion relative to the blob.
    // The length scale is the BLOB, not the particle: the outside fluid shears
    // over the size of the lump it is flowing around. Written with the particle
    // radius it comes out 4.5 mu / (rho rp^2), which for rp = 1.3 mm is 27 per
    // second -- fourteen times Lamb's damping rate for a 15 mm drop, enough to
    // kill a shape oscillation inside a fifth of its period. With the blob
    // radius it lands at 0.20 per second against Lamb's 0.22. Same formula, and
    // the only difference is which length you put in it.
    const kSnum = 4.5 * muAq;
    const rg = env.reducedG;
    for (let i = 0; i < n; i++) {
      const Tf = column.at(y[i]);
      const rw = rw20 * (1 - WAX.beta * (T[i] - 20));
      const ra = rhoAq(Tf);
      // env.reducedG is a CALIBRATION HOOK: it replaces the buoyancy
      // calculation with a fixed downward reduced gravity, so a test can impose
      // a known Archimedes number instead of arranging temperatures to produce
      // one. The lamp itself never sets it.
      const ay = rg !== undefined ? -rg : g * (ra - rw) / rw;   // up positive
      const a = this.blobA[i] > this.rp ? this.blobA[i] : this.rp;
      const kD = 3 * SOLVER.Cd * ra / (8 * rw * a);   // blob-scale quadratic drag
      // Drag acts on motion RELATIVE to the fluid, and the fluid is moving --
      // pushed the other way by this very blob.
      const vaq = hasFlow ? column.vAq(y[i]) : 0;
      const bvx = this.blobVx[i], bvy = this.blobVy[i] - vaq, bvz = this.blobVz[i];
      const sp = Math.sqrt(bvx * bvx + bvy * bvy + bvz * bvz);
      const kS = kSnum / (rw * a * a);
      const mob = 1 - 0.98 * solid[i];                // solid wax barely moves
      vx[i] += dt * (-kD * sp * bvx - kS * (vx[i] - this.blobVx[i])) * mob;
      vy[i] += dt * (ay - kD * sp * bvy - kS * (vy[i] - this.blobVy[i])) * mob;
      vz[i] += dt * (-kD * sp * bvz - kS * (vz[i] - this.blobVz[i])) * mob;
    }

    // ---- 2. predict -----------------------------------------------------
    for (let i = 0; i < n; i++) {
      px[i] = x[i] + dt * vx[i];
      py[i] = y[i] + dt * vy[i];
      pz[i] = z[i] + dt * vz[i];
      this.confine(i);
    }
    this.buildGrid();
    this.buildNeighbours();

    // ---- 3. density projection (PBF) -----------------------------------
    const w0 = kPoly * h2 * h2 * h2;                  // W(0)
    const scale = m / rho0;
    for (let it = 0; it < iters; it++) {
      for (let i = 0; i < n; i++) {
        let d = w0;
        const base = i * MAXN, cnt = nbrCount[i];
        const xi = px[i], yi = py[i], zi = pz[i];
        let gx = 0, gy = 0, gz = 0, sumGrad = 0;
        for (let k = 0; k < cnt; k++) {
          const j = nbr[base + k];
          const ddx = xi - px[j], ddy = yi - py[j], ddz = zi - pz[j];
          const r2 = ddx * ddx + ddy * ddy + ddz * ddz;
          if (r2 >= h2 || r2 < 1e-24) continue;
          const q = h2 - r2;
          d += kPoly * q * q * q;
          const r = Math.sqrt(r2);
          const c = kSpiky * (h - r) * (h - r) / r * scale;
          const wx = c * ddx, wy = c * ddy, wz = c * ddz;
          gx += wx; gy += wy; gz += wz;
          sumGrad += wx * wx + wy * wy + wz * wz;
        }
        d *= m;
        rho[i] = d;
        sumGrad += gx * gx + gy * gy + gz * gz;
        // A liquid resists compression, not expansion. Free-surface particles
        // sit at rho < rho0 and must NOT be sucked inward by the density
        // constraint -- that is what interfacial tension is for. Without this
        // clamp a lone particle has sum|grad C|^2 -> 0 with C ~ -0.8, and the
        // correction goes to infinity: the sim detonates on the first step.
        const C = d / rho0 - 1;
        lam[i] = C > 0 ? -C / (sumGrad + this.eps) : 0;
      }
      const iwDq = 1 / this.wDq, scorr = this.scorr;
      for (let i = 0; i < n; i++) {
        let cx = 0, cy = 0, cz = 0;
        const base = i * MAXN, cnt = nbrCount[i];
        const xi = px[i], yi = py[i], zi = pz[i], li = lam[i];
        for (let k = 0; k < cnt; k++) {
          const j = nbr[base + k];
          const ddx = xi - px[j], ddy = yi - py[j], ddz = zi - pz[j];
          const r2 = ddx * ddx + ddy * ddy + ddz * ddz;
          if (r2 >= h2 || r2 < 1e-24) continue;
          const q = h2 - r2;
          const wp = kPoly * q * q * q * iwDq;
          const sc = -scorr * wp * wp * wp * wp;        // artificial pressure
          const r = Math.sqrt(r2);
          const c = kSpiky * (h - r) * (h - r) / r * scale;
          const w = (li + lam[j] + sc) * c;
          cx += w * ddx; cy += w * ddy; cz += w * ddz;
        }
        dpx[i] = cx; dpy[i] = cy; dpz[i] = cz;
      }
      for (let i = 0; i < n; i++) {
        const k = 1 - 0.92 * solid[i];      // unmelted wax holds its shape
        px[i] += dpx[i] * k; py[i] += dpy[i] * k; pz[i] += dpz[i] * k;
        this.confine(i);
      }
    }

    // ---- 4. velocity update -------------------------------------------
    const inv = 1 / dt;
    for (let i = 0; i < n; i++) {
      vx[i] = (px[i] - x[i]) * inv;
      vy[i] = (py[i] - y[i]) * inv;
      vz[i] = (pz[i] - z[i]) * inv;
      x[i] = px[i]; y[i] = py[i]; z[i] = pz[i];
      const f = 1 - solid[i];
      vx[i] *= f; vy[i] *= f; vz[i] *= f;
    }

    this.findBlobs();

    // ---- 5. interfacial tension (Akinci) ------------------------------
    const nxA = this.nx, nyA = this.ny, nzA = this.nz;
    for (let i = 0; i < n; i++) {
      let gx = 0, gy = 0, gz = 0;
      const base = i * MAXN, cnt = nbrCount[i];
      const xi = x[i], yi = y[i], zi = z[i];
      for (let k = 0; k < cnt; k++) {
        const j = nbr[base + k];
        const ddx = xi - x[j], ddy = yi - y[j], ddz = zi - z[j];
        const r2 = ddx * ddx + ddy * ddy + ddz * ddz;
        if (r2 >= h2 || r2 < 1e-24) continue;
        const r = Math.sqrt(r2);
        const c = kSpiky * (h - r) * (h - r) / r * (m / (rho[j] > 1e-6 ? rho[j] : 1e-6));
        gx += c * ddx; gy += c * ddy; gz += c * ddz;
      }
      nxA[i] = h * gx; nyA[i] = h * gy; nzA[i] = h * gz;
    }
    const coh = env.cohesion ?? 1;
    const cohA = SOLVER.cohK * sigma * coh * this.cohScale;
    const disjoinA = cohA * (env.disjoin ?? SOLVER.disjoinK);
    const idArr = this.id;
    const curvA = cohA * m * this.cohRef * SOLVER.curvK;
    const kCoh = this.kCoh, h6_64 = h ** 6 / 64;
    for (let i = 0; i < n; i++) {
      let fx = 0, fy = 0, fz = 0;
      const base = i * MAXN, cnt = nbrCount[i];
      const xi = x[i], yi = y[i], zi = z[i], ri = rho[i];
      for (let k = 0; k < cnt; k++) {
        const j = nbr[base + k];
        const ddx = xi - x[j], ddy = yi - y[j], ddz = zi - z[j];
        const r2 = ddx * ddx + ddy * ddy + ddz * ddz;
        if (r2 >= h2 || r2 < 1e-24) continue;
        const r = Math.sqrt(r2);
        const hr = h - r;
        const a3 = hr * hr * hr * r2 * r;
        const C = 2 * r > h ? kCoh * a3 : kCoh * (2 * a3 - h6_64);
        const Kij = 2 * rho0 / (ri + rho[j]);
        // Same spline, opposite sign, across a film that has not drained: the
        // surfactant layer holds the two interfaces apart instead of letting
        // them pull together.
        if (idArr[i] !== idArr[j]) {
          const w = disjoinA * m * Math.abs(C) / r;
          fx += Kij * w * ddx; fy += Kij * w * ddy; fz += Kij * w * ddz;
        } else {
          const w = -cohA * m * C / r;
          fx += Kij * (w * ddx - curvA * (nxA[i] - nxA[j]));
          fy += Kij * (w * ddy - curvA * (nyA[i] - nyA[j]));
          fz += Kij * (w * ddz - curvA * (nzA[i] - nzA[j]));
        }
      }
      const mob = (1 - 0.9 * solid[i]) * dt;
      vx[i] += fx * mob; vy[i] += fy * mob; vz[i] += fz * mob;
    }

    // ---- 6. wax viscosity, as XSPH velocity blending -----------------
    // An explicit Monaghan laplacian at this smoothing length is right at its
    // own stability limit for dt = 1/60, and a random velocity field walks it
    // straight over the edge. XSPH blending is the same diffusion written
    // unconditionally stable: nu_eff ~ c * h^2 / (kVisc2 * dt), so invert that
    // for c and clamp at the bound. If the clamp binds, the timestep cannot
    // resolve the requested viscosity and viscClamped says so rather than
    // silently simulating a thinner wax.
    const nu = WAX.mu / WAX.rho20;
    let cvis = nu * dt / (0.10 * h2);
    this.viscClamped = cvis > 0.5;
    if (cvis > 0.5) cvis = 0.5;
    // any extra smoothing has to be per second, not per step, or halving the
    // timestep doubles the damping
    const blend = cvis + SOLVER.xsph * dt;
    for (let i = 0; i < n; i++) {
      let sx = 0, sy = 0, sz = 0, wsum = 0;
      const base = i * MAXN, cnt = nbrCount[i];
      const xi = x[i], yi = y[i], zi = z[i];
      const vxi = vx[i], vyi = vy[i], vzi = vz[i];
      for (let k = 0; k < cnt; k++) {
        const j = nbr[base + k];
        const ddx = xi - x[j], ddy = yi - y[j], ddz = zi - z[j];
        const r2 = ddx * ddx + ddy * ddy + ddz * ddz;
        if (r2 >= h2) continue;
        const q = h2 - r2;
        const wq = kPoly * q * q * q;
        sx += (vx[j] - vxi) * wq; sy += (vy[j] - vyi) * wq; sz += (vz[j] - vzi) * wq;
        wsum += wq;
      }
      if (wsum > 1e-9) {
        const b = blend / wsum;
        dpx[i] = sx * b; dpy[i] = sy * b; dpz[i] = sz * b;
      } else { dpx[i] = 0; dpy[i] = 0; dpz[i] = 0; }
    }
    for (let i = 0; i < n; i++) { vx[i] += dpx[i]; vy[i] += dpy[i]; vz[i] += dpz[i]; }

    // ---- 6b. CFL safety net --------------------------------------------
    // A liquid blob in a lava lamp moves at centimetres per second. If anything
    // here is doing metres per second, the solver has failed, not the lamp. The
    // clamp keeps a bad step from destroying the state -- and it COUNTS, so the
    // diagnostics can say it fired instead of quietly papering over a bug.
    const vLim = this.vMax / dt;
    const vLim2 = vLim * vLim;
    for (let i = 0; i < n; i++) {
      const s2 = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
      if (s2 > vLim2) {
        const s = vLim / Math.sqrt(s2);
        vx[i] *= s; vy[i] *= s; vz[i] *= s;
        this.clamps++;
      }
    }

    // ---- 7. heat ----------------------------------------------------
    if (!env.noHeat) this.heat(dt, column, env);
    if (this.stirTimer > 0) this.stirTimer = Math.max(0, this.stirTimer - dt);
  }

  cohesionSpline(r) {
    const h = this.h;
    if (r > h) return 0;
    const a = (h - r) ** 3 * r ** 3;
    return 2 * r > h ? this.kCoh * a : this.kCoh * (2 * a - h ** 6 / 64);
  }

  // Push a particle back inside the globe.
  confine(i) {
    const rp = this.rp;
    let yy = this.py[i];
    if (yy < rp) { yy = rp; this.py[i] = rp; }
    else { const t = GEO.H - rp; if (yy > t) { yy = t; this.py[i] = t; } }
    const R = this.radiusAt(yy) - rp;
    const a = this.px[i], b = this.pz[i];
    const r2 = a * a + b * b;
    if (r2 > R * R && r2 > 1e-24) {
      const s = R / Math.sqrt(r2);
      this.px[i] = a * s; this.pz[i] = b * s;
    }
  }

  // -------------------------------------------------------------------------
  // Heat: plate conduction into the pool, film exchange with the fluid over the
  // blob surface, conduction within the blob, and latent heat across the melt
  // band via the apparent-heat-capacity method.
  // -------------------------------------------------------------------------
  heat(dt, column, env = {}) {
    const sub = Math.max(1, Math.ceil(dt / this.heatDtMax));
    // The deposits back into the column are POWERS, and the column applies them
    // over the whole outer dt. Substepping without weighting them by dt_sub/dt
    // hands the fluid `sub` times the heat the wax actually gave up -- energy
    // created out of a loop counter. tools/lamp-probe.mjs now audits the balance
    // for exactly this reason.
    const w = 1 / sub;
    this.plateQ = 0;
    for (let s = 0; s < sub; s++) this.heatStep(dt / sub, column, env, w);
    column.plateWaxQ = this.plateQ;
  }

  heatStep(dt, column, env = {}, weight = 1) {
    const n = this.n, h = this.h, h2 = this.h2, m = this.m;
    const x = this.x, y = this.y, z = this.z, T = this.T, rho = this.rho;
    const nbr = this.nbr, nbrCount = this.nbrCount, dT = this.dT;
    const Pr = (AQ.mu * AQ.c) / AQ.k;
    const cbrtPr = Math.cbrt(Pr);
    const plateBand = 2.5 * this.rp;
    dT.fill(0);
    let plateQ = 0;

    // The wax touching the plate cannot present more contact area than the plate
    // has. Charging every particle in the bottom film its own disc of pi*rp^2
    // over-counts the contact by the number of particle layers in the film, and
    // the plate then bleeds most of the bulb into the pool and never gets hot.
    // And the plate area the wax covers is area the FLUID does not: the two
    // films share one plate, so the column scales its own film by what is left.
    let onPlate = 0;
    for (let i = 0; i < n; i++) if (y[i] < plateBand) onPlate++;
    const plateArea = Math.PI * this.radiusAt(0) * this.radiusAt(0);
    const cover = Math.min(1, onPlate * Math.PI * this.rp * this.rp / plateArea);
    const plateK = onPlate > 0 ? HEAT.hPlate * plateArea * cover / onPlate : 0;
    column.plateWaxFrac = cover;
    // energy the wax is holding, for the audit
    let U = 0;
    for (let i = 0; i < n; i++) U += T[i];
    this.enthalpy = m * WAX.c * U + m * WAX.L * (n - this.solidCount());

    for (let i = 0; i < n; i++) {
      const Tf = column.at(y[i]);
      const a = this.blobA[i] > this.rp ? this.blobA[i] : this.rp;
      const bvx = this.blobVx[i], bvy = this.blobVy[i], bvz = this.blobVz[i];
      const spd = Math.sqrt(bvx * bvx + bvy * bvy + bvz * bvz);
      const Re = AQ.rho20 * spd * 2 * a / AQ.mu;
      const Nu = 2 + 0.6 * Math.sqrt(Re) * cbrtPr;
      const hf = Nu * AQ.k / (2 * a);
      const q = hf * this.expose[i] * (Tf - T[i]);        // W
      dT[i] += q;
      column.depositWaxHeat(y[i], -q * weight);
      if (y[i] < plateBand) {
        const qp = plateK * (column.Tplate - T[i]);
        dT[i] += qp;
        plateQ += qp * weight;
      }
    }
    // conduction inside the wax
    const kc = WAX.k * this.Vp * this.kVisc, mm = m;
    for (let i = 0; i < n; i++) {
      const base = i * MAXN, cnt = nbrCount[i];
      const xi = x[i], yi = y[i], zi = z[i], Ti = T[i];
      let s = 0;
      for (let k = 0; k < cnt; k++) {
        const j = nbr[base + k];
        const ddx = xi - x[j], ddy = yi - y[j], ddz = zi - z[j];
        const r2 = ddx * ddx + ddy * ddy + ddz * ddz;
        if (r2 >= h2) continue;
        const r = Math.sqrt(r2);
        s += (mm / (rho[j] > 1e-6 ? rho[j] : 1e-6)) * (T[j] - Ti) * (h - r);
      }
      dT[i] += kc * s;
    }
    // integrate with apparent heat capacity across the melt band
    const band = 2 * WAX.dTmelt;
    const cLatent = WAX.L / band;
    const lo = WAX.Tmelt - WAX.dTmelt, hi = WAX.Tmelt + WAX.dTmelt;
    for (let i = 0; i < n; i++) {
      const t = T[i];
      const c = (t > lo && t < hi) ? WAX.c + cLatent : WAX.c;
      const nt = t + dt * dT[i] / (m * c);
      T[i] = nt;
      let s = (hi - nt) / band;
      this.solid[i] = s < 0 ? 0 : (s > 1 ? 1 : s);
    }
    this.plateQ += plateQ;
  }

  // Shaking a lava lamp emulsifies it: the blobs shatter into droplets too
  // small to rise, and it takes a long time for them to find each other again.
  shake(strength = 0.05) {
    for (let i = 0; i < this.n; i++) {
      if (this.solid[i] > 0.5) continue;
      this.vx[i] += (Math.random() * 2 - 1) * strength;
      this.vy[i] += (Math.random() * 2 - 1) * strength;
      this.vz[i] += (Math.random() * 2 - 1) * strength;
    }
    this.stirTimer = 4;
  }

  // Diagnostics -------------------------------------------------------------
  solidCount() { let c = 0; for (let i = 0; i < this.n; i++) c += this.solid[i]; return c; }

  get moltenFraction() {
    let s = 0;
    for (let i = 0; i < this.n; i++) s += 1 - this.solid[i];
    return s / this.n;
  }
  get meanT() { let s = 0; for (let i = 0; i < this.n; i++) s += this.T[i]; return s / this.n; }
  get maxT() { let v = -Infinity; for (let i = 0; i < this.n; i++) if (this.T[i] > v) v = this.T[i]; return v; }
  get maxRise() { let v = 0; for (let i = 0; i < this.n; i++) if (this.blobVy[i] > v) v = this.blobVy[i]; return v; }

  // Volume-weighted mean upward speed over blobs that are actually moving up.
  // maxRise is the fastest particle in the tank at one instant; taken as a
  // running maximum over minutes it reports the worst transient the run ever
  // had and calls it the rise speed. This is the number to compare against
  // "a lava lamp blob rises at one to three centimetres a second".
  typicalRise() {
    const seen = new Map();
    for (let i = 0; i < this.n; i++) {
      const r = this.blob[i];
      if (!seen.has(r)) seen.set(r, { n: 0, vy: this.blobVy[i] });
      seen.get(r).n++;
    }
    let num = 0, den = 0;
    for (const b of seen.values()) {
      if (b.n < 8 || b.vy <= 0) continue;
      num += b.n * b.vy; den += b.n;
    }
    return den > 0 ? num / den : 0;
  }
  get centreOfMass() { let s = 0; for (let i = 0; i < this.n; i++) s += this.y[i]; return s / this.n; }

  // How much molten wax is currently lighter than the fluid around it, how
  // strongly, and how big a blob that buoyancy could tear off the pool.
  //
  // This is the gate for fast-forward, and the sign of the buoyancy alone is not
  // enough. One particle in a still-frozen plug goes marginally buoyant a full
  // hour before the lamp does anything, because a parcel that light cannot beat
  // the interfacial tension holding it to the pool. The physical threshold is
  // the capillary length: buoyancy only wins once the buoyant volume reaches
  //     a_c = sqrt(sigma / (drho * g))
  // which is the radius at which the buoyant pull on a blob overtakes the
  // surface tension trying to keep it flat. Below that, nothing lifts, so
  // advancing the heat equation alone loses no mechanics.
  buoyantSpan(column, drho = WAX.rho20 - AQ.rho20, sigma = IFACE.sigma) {
    const rw20 = AQ.rho20 + drho;
    let vol = 0, best = -Infinity, peakDrho = 0;
    for (let i = 0; i < this.n; i++) {
      if (this.solid[i] > 0.5) continue;
      const rw = rw20 * (1 - WAX.beta * (this.T[i] - 20));
      const ra = rhoAq(column.at(this.y[i]));
      const d = ra - rw;
      if (d > 0) { vol += this.Vp; if (d > peakDrho) peakDrho = d; }
      const a = HEAT.g * d / rw;
      if (a > best) best = a;
    }
    const ac = peakDrho > 1e-3
      ? Math.sqrt(sigma / (peakDrho * HEAT.g))
      : Infinity;
    const capVolume = ac === Infinity ? Infinity : (4 / 3) * Math.PI * ac * ac * ac;
    return {
      volume: vol,
      capillaryRadius: ac,
      capVolume,
      peakDrho,
      maxAccel: best === -Infinity ? -HEAT.g : best,
      lifts: vol >= capVolume,
    };
  }

  // Kept for the diagnostics: peak buoyant acceleration in the molten wax.
  buoyancyHeadroom(column, drho = WAX.rho20 - AQ.rho20) {
    return this.buoyantSpan(column, drho).maxAccel;
  }
}
