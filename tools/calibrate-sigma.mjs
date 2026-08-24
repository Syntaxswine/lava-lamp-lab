// calibrate-sigma.mjs — measure the interfacial tension the solver actually has,
// against a case with a closed-form answer, and report the coefficient that
// makes it match the tension we asked for.
//
// THE REFERENCE. A non-wetting liquid puddle on a flat floor cannot be made
// thicker by adding liquid: past a point it just spreads. The maximum thickness
// is set by the balance of hydrostatic pressure against surface tension,
//
//     e = 2 * sqrt( sigma / (drho * g) ) * sin(theta/2)
//
// (de Gennes, Brochard-Wyart and Quere, "Capillarity and Wetting Phenomena",
// section 2.3). The wax here has no adhesion to the glass, so theta = 180 deg
// and sin(theta/2) = 1: the puddle settles at exactly twice the capillary
// length. That makes the thickness of a settled puddle a direct reading of
// sigma, with no unknown constant in front of it.
//
// WHY NOT THE OSCILLATION TEST. The first version of this tool used Rayleigh's
// drop oscillation, omega^2 = 8 sigma / (rho a^3). It does not work here: with
// the wax viscosity and the aqueous drag both in play, a centimetre-scale drop
// at 3 mN/m is overdamped -- it overshoots once and creeps back without ringing,
// so there is no period to measure. Reading a frequency off that trace produced
// half-periods scattered over a factor of fifty and a "sigma" that moved by
// three orders of magnitude depending on which crossings you kept. A static
// measurement has no such problem.
//
// WHAT MAKES IT A MEASUREMENT AND NOT A NUMBER. e ~ (drho g)^(-1/2) is the
// signature of capillarity specifically. The tool runs several reduced gravities
// and fits the exponent. If the exponent is not near -1/2, then whatever is
// holding the puddle up is not surface tension, and reporting a sigma from it
// would be reporting the wrong quantity to three decimal places.
//
// AND A RESOLUTION FLOOR. A puddle only a couple of particle layers deep cannot
// be thinner than the particles, so at high reduced gravity the measurement
// stops reading capillarity and starts reading the discretisation -- which bends
// the fitted exponent toward zero and inflates sigma. Every point reports how
// many particle layers deep it is, and points below MIN_LAYERS are shown but
// excluded from the fit rather than quietly averaged in.
//
// usage: node tools/calibrate-sigma.mjs [particles] [cohK]

import { Wax } from '../js/wax.js';
import { WAX, IFACE, SOLVER, GEO, AQ } from '../js/params.js';

const n = Number(process.argv[2] ?? 2400);
if (process.argv[3]) SOLVER.cohK = Number(process.argv[3]);

const T0 = 45;
const bath = {
  Tplate: T0, plateWaxQ: 0, plateWaxFrac: 0,
  at: () => T0,
  depositWaxHeat: () => {},
};

// Settle a puddle of `volume` under a fixed reduced gravity and return its
// thickness, measured where it is flat rather than at the rim.
function puddle(volume, gRed, seconds = 10) {
  const wax = new Wax(n, volume);
  const dx = wax.dx;
  // start as a squat cylinder wider than it will end up, so it settles DOWN
  // into the equilibrium rather than spreading out to it
  const R0 = Math.sqrt(volume / (Math.PI * 0.010));
  let placed = 0, y = wax.rp;
  while (placed < n && y < GEO.H) {
    const R = Math.min(R0, wax.radiusAt(y) - 1.3 * wax.rp);
    const cols = Math.max(1, Math.floor(2 * R / dx));
    const off = -0.5 * (cols - 1) * dx;
    for (let a = 0; a < cols && placed < n; a++)
      for (let b = 0; b < cols && placed < n; b++) {
        const cx = off + a * dx, cz = off + b * dx;
        if (cx * cx + cz * cz > R * R) continue;
        wax.x[placed] = cx; wax.z[placed] = cz; wax.y[placed] = y;
        wax.vx[placed] = wax.vy[placed] = wax.vz[placed] = 0;
        wax.T[placed] = T0; wax.solid[placed] = 0; wax.rho[placed] = wax.rho0;
        wax.px[placed] = cx; wax.py[placed] = y; wax.pz[placed] = cz;
        placed++;
      }
    y += dx;
  }
  if (placed < n) throw new Error(`placed ${placed}/${n}`);
  wax.clamps = 0;
  wax.buildGrid(); wax.buildNeighbours(); wax.findBlobs();

  const env = {
    reducedG: gRed, noHeat: true, sigma: IFACE.sigma,
    muAq: AQ.mu, iterations: 4, cohesion: 1,
  };
  const dt = 1 / 480;
  const steps = Math.round(seconds / dt);
  const history = [];
  for (let s = 0; s < steps; s++) {
    wax.step(dt, bath, env);
    if (s > steps * 0.7 && s % 40 === 0) history.push(thickness(wax));
  }
  const e = history.reduce((a, b) => a + b, 0) / history.length;
  const drift = Math.max(...history) - Math.min(...history);
  return { e, drift, wax, layers: e / wax.dx };
}

// Thickness over the flat middle: the top surface height inside the inner half
// of the puddle radius, minus the floor. The rim is curved and would drag the
// average down.
function thickness(wax) {
  let rmax = 0;
  for (let i = 0; i < wax.n; i++) {
    const r = Math.hypot(wax.x[i], wax.z[i]);
    if (r > rmax) rmax = r;
  }
  const cut = 0.5 * rmax;
  let top = 0, cnt = 0;
  const bins = 24;
  const colTop = new Float64Array(bins * bins).fill(-1);
  for (let i = 0; i < wax.n; i++) {
    const x = wax.x[i], z = wax.z[i];
    if (x * x + z * z > cut * cut) continue;
    const bx = Math.min(bins - 1, Math.max(0, Math.floor((x / cut + 1) * bins / 2)));
    const bz = Math.min(bins - 1, Math.max(0, Math.floor((z / cut + 1) * bins / 2)));
    const k = bz * bins + bx;
    if (wax.y[i] > colTop[k]) colTop[k] = wax.y[i];
  }
  for (let k = 0; k < colTop.length; k++) if (colTop[k] >= 0) { top += colTop[k]; cnt++; }
  return cnt ? top / cnt + wax.rp : 0;      // + rp: surface is a radius above the last centre
}

console.log(`n = ${n}   cohK = ${SOLVER.cohK}   curvK = ${SOLVER.curvK}   ` +
  `target sigma = ${(IFACE.sigma * 1e3).toFixed(2)} mN/m`);
console.log('');
console.log('  g_red    predicted e   measured e   drift  layers   sigma_eff mN/m  clamps blobs  use');

const rho = WAX.rho20;
const V = 32e-6;                        // 32 mL: wide enough to be a puddle
const rows = [];
const MIN_LAYERS = 4.5;
for (const gRed of [0.05, 0.08, 0.13, 0.21]) {
  const ac = Math.sqrt(IFACE.sigma / (rho * gRed));
  let r;
  try { r = puddle(V, gRed); } catch (e) { console.log(`  ${gRed}  skipped: ${e.message}`); continue; }
  const sigmaEff = rho * gRed * (r.e / 2) ** 2;
  const ok = r.layers >= MIN_LAYERS;
  if (ok) rows.push({ gRed, e: r.e, sigmaEff });
  console.log(`  ${gRed.toFixed(3)}    ${(2 * ac * 1e3).toFixed(2).padStart(8)} mm   ` +
    `${(r.e * 1e3).toFixed(2).padStart(7)} mm  ${(r.drift * 1e3).toFixed(2).padStart(5)} mm  ` +
    `${r.layers.toFixed(1).padStart(5)}   ${(sigmaEff * 1e3).toFixed(3).padStart(12)}  ` +
    `${String(r.wax.clamps).padStart(6)} ${String(r.wax.blobCount).padStart(5)}  ` +
    `${ok ? 'yes' : 'NO (under-resolved)'}`);
}

console.log('');
if (rows.length < 3) {
  console.log('REFUSING TO REPORT: fewer than three usable gravities.');
  process.exit(2);
}
// fit log e = c + p * log(g); capillarity says p = -1/2
const lx = rows.map((r) => Math.log(r.gRed));
const ly = rows.map((r) => Math.log(r.e));
const mx = lx.reduce((a, b) => a + b) / lx.length;
const my = ly.reduce((a, b) => a + b) / ly.length;
let num = 0, den = 0;
for (let i = 0; i < lx.length; i++) { num += (lx[i] - mx) * (ly[i] - my); den += (lx[i] - mx) ** 2; }
const p = num / den;
const vals = rows.map((r) => r.sigmaEff);
const mean = vals.reduce((a, b) => a + b) / vals.length;
const spread = (Math.max(...vals) - Math.min(...vals)) / mean;

console.log(`thickness exponent  ${p.toFixed(3)}   (capillarity requires -0.500)`);
console.log(`sigma_eff spread    ${(spread * 100).toFixed(0)}%   mean ${(mean * 1e3).toFixed(3)} mN/m`);
if (Math.abs(p + 0.5) > 0.15) {
  console.log('');
  console.log('EXPONENT IS WRONG. The puddle is not being held up by surface');
  console.log('tension, so its thickness is not a reading of sigma. Not reporting');
  console.log('a calibration from it.');
  process.exit(3);
}
console.log('');
console.log(`suggested cohK = ${(SOLVER.cohK * IFACE.sigma / mean).toFixed(3)}`);
