// volume.js — turn the wax particles into a 3-D field the renderer can march.
//
// Three channels, packed into one RGBA8 texture:
//   R  wax fraction (0 = clear fluid, 1 = wax at rest density)
//   G  temperature, mapped over TMIN..TMAX
//   A  how much bulb light reached this voxel (see the prefix sum below)
//
// There is no solid-fraction channel. Solid fraction is a pure function of
// temperature -- the same clamp wax.js applies -- so carrying it here would be
// storing a value twice and paying a third of the splat cost to do it. The
// shader derives it from G.
//
// The splat kernel is wider than a particle, because a lava lamp blob is a
// smooth-walled thing and sampling particles at their own radius renders the
// discretisation instead of the blob. But not much wider: the kernel radius is
// what sets the per-particle cell box, and the box is cubic in it. At 1.7 dx
// over a 2.4 mm grid this loop was 343 cells per particle and 7 ms a frame --
// the entire cost of drawing, while the GPU raymarch it feeds took 0.3 ms.

import { GEO, globeRadius } from './params.js';

export const TMIN = 18;
export const TMAX = 78;

export class Volume {
  constructor(nx = 28, ny = 82, nz = 28) {
    this.nx = nx; this.ny = ny; this.nz = nz;
    let rmax = 0;
    for (const [, r] of GEO.profile) rmax = Math.max(rmax, r);
    this.rmax = rmax;
    // bounds: a box that just contains the globe
    this.min = [-rmax, 0, -rmax];
    this.max = [rmax, GEO.H, rmax];
    this.dx = (this.max[0] - this.min[0]) / nx;
    this.dy = (this.max[1] - this.min[1]) / ny;
    this.dz = (this.max[2] - this.min[2]) / nz;
    const n = nx * ny * nz;
    this.w = new Float32Array(n);
    this.wT = new Float32Array(n);
    this.data = new Uint8Array(n * 4);
    // radius profile sampled for the shader
    this.profileN = 64;
    this.profile = new Float32Array(this.profileN);
    for (let i = 0; i < this.profileN; i++) {
      this.profile[i] = globeRadius(GEO.H * i / (this.profileN - 1));
    }
  }

  // Splat the wax into the field. `scale` converts accumulated kernel weight
  // into a wax fraction; it is set so a particle at rest spacing reads 1.0.
  build(wax) {
    const { nx, ny, nz } = this;
    this.w.fill(0); this.wT.fill(0);
    const R = 1.5 * wax.dx;                 // splat radius, metres
    const R2 = R * R;
    const kx = Math.ceil(R / this.dx), ky = Math.ceil(R / this.dy), kz = Math.ceil(R / this.dz);
    // Normalisation, analytic rather than fudged. A voxel buried in bulk wax
    // collects (number density) * integral of the kernel over its support:
    //   integral (1 - r^2/R^2)^2 dV = 4*pi*R^3 * 8/105 = 0.95760 R^3
    // and the number density of wax at rest is 1/Vp.
    const inv = wax.Vp / (0.957601 * R * R * R);

    for (let p = 0; p < wax.n; p++) {
      const px = wax.x[p], py = wax.y[p], pz = wax.z[p];
      const ci = ((px - this.min[0]) / this.dx) | 0;
      const cj = ((py - this.min[1]) / this.dy) | 0;
      const ck = ((pz - this.min[2]) / this.dz) | 0;
      const T = wax.T[p];
      for (let k = Math.max(0, ck - kz); k <= Math.min(nz - 1, ck + kz); k++) {
        const dz = (this.min[2] + (k + 0.5) * this.dz) - pz;
        const dz2 = dz * dz;
        if (dz2 > R2) continue;
        for (let j = Math.max(0, cj - ky); j <= Math.min(ny - 1, cj + ky); j++) {
          const dy = (this.min[1] + (j + 0.5) * this.dy) - py;
          const dyz = dz2 + dy * dy;
          if (dyz > R2) continue;
          const rowBase = (k * ny + j) * nx;
          for (let i = Math.max(0, ci - kx); i <= Math.min(nx - 1, ci + kx); i++) {
            const dx = (this.min[0] + (i + 0.5) * this.dx) - px;
            const d2 = dyz + dx * dx;
            if (d2 > R2) continue;
            const t = 1 - d2 / R2;
            const wgt = t * t;
            const idx = rowBase + i;
            this.w[idx] += wgt;
            this.wT[idx] += wgt * T;
          }
        }
      }
    }

    const d = this.data;
    const invT = 255 / (TMAX - TMIN);
    for (let idx = 0, n = this.w.length; idx < n; idx++) {
      const wv = this.w[idx];
      const o = idx * 4;
      if (wv <= 1e-8) { d[o] = 0; d[o + 1] = 128; d[o + 3] = 255; continue; }
      d[o] = Math.min(255, wv * inv * 255);
      d[o + 1] = Math.max(0, Math.min(255, (this.wT[idx] / wv - TMIN) * invT));
      d[o + 3] = 255;
    }

    // The bulb sits directly under the globe, so light travels straight up the
    // y axis and the optical depth a voxel sees is just the wax below it in its
    // own column. That makes shadowing a prefix sum, computed here once, packed
    // into the alpha channel, and free in the shader -- one texture fetch gives
    // density, temperature, solid fraction AND how much light reached it.
    // Extinction for the LIGHT path, not the view path. Dyed wax scatters
    // strongly forward, so transmitted light is attenuated far less than the
    // 110/m the view ray sees: the similarity relation gives an effective
    // sigma' = sigma (1 - g), and g ~ 0.77 for wax puts that near 25/m. This is
    // why a running lamp glows all the way up instead of going dark above the
    // bottom blob.
    const kAbs = 25;                        // 1/m, reduced extinction
    for (let k = 0; k < nz; k++) {
      for (let i = 0; i < nx; i++) {
        let tau = 0;
        for (let j = 0; j < ny; j++) {
          const o = ((k * ny + j) * nx + i) * 4;
          d[o + 3] = 255 * Math.exp(-kAbs * tau);
          tau += (d[o] / 255) * this.dy;
        }
      }
    }
    return d;
  }
}
