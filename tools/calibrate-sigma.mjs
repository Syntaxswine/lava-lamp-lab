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
// usage: node tools/calibrate-sigma.mjs [particles] [cohK] [settle-seconds]

import { Wax } from '../js/wax.js';
import { WAX, IFACE, SOLVER, GEO, AQ } from '../js/params.js';

const n = Number(process.argv[2] ?? 3200);
if (process.argv[3]) SOLVER.cohK = Number(process.argv[3]);
const settleSeconds = Number(process.argv[4] ?? 10);

// CALIBRATE AT THE SPACING THE LAMP RUNS AT.
//
// The uncorrected pair acceleration scaled as dx^-2 while the Laplace-pressure
// argument asks for dx^-1. Wax now carries the missing reference length, but the
// calibration still runs at the shipping spacing: a transfer law is not an
// excuse to validate a different configuration. The puddle volume is therefore
// derived from the lamp's own volume-per-particle rather than fixed.
//
// The calibration vessel below is deliberately wider than the lamp.  That lets
// 3200 particles make a 107 mL puddle at the shipping spacing: broad enough for
// the infinite-puddle thickness law, without wall confinement supplying part of
// the answer.
const SPACING3 = WAX.volume / SOLVER.particles;    // m^3 per particle, as shipped

const T0 = 45;
const bath = {
  Tplate: T0, plateWaxQ: 0, plateWaxFrac: 0,
  at: () => T0,
  depositWaxHeat: () => {},
};

// Settle a puddle of `volume` under a fixed reduced gravity and return its
// thickness, measured where it is flat rather than at the rim.
function puddle(volume, gRed, seconds = settleSeconds) {
  const wax = new Wax(n, volume);
  // Calibration needs a flat, wide vessel rather than the lamp's tapered
  // globe.  At the shipping 3.22 mm spacing the old rig clipped the puddle at
  // the 27.5 mm foot, so wall confinement contaminated the thickness exactly
  // where the tool was supposed to isolate capillarity.  Wax caches its wall
  // profile, which makes a cylindrical calibration vessel a one-line override
  // without adding test-only branches to the solver.
  wax.profR.fill(0.060);
  // This instrument measures the cohesion law, so the surfactant-film state
  // machine must be absent, not merely assigned zero force.  `disjoin: 0` by
  // itself still lets a split component acquire a fresh persistent id; those
  // ids then refuse to bond again and the puddle turns into hundreds of pieces.
  // Keep one material identity while retaining the normal connected-component
  // pass for the blob and surface-area diagnostics.
  wax.relabel = function calibrationRelabel() { this.id.fill(1); };
  wax.id.fill(1);
  const dx = wax.dx;
  // Start above the analytic thickness so the measurement has to settle down
  // to the answer.  The old hard-coded 10 mm start was actually *below* every
  // predicted thickness in this sweep (15--24 mm), contradicting its own
  // comment and asking a viscous puddle to contract upward in only ten seconds.
  const predicted = 2 * Math.sqrt(IFACE.sigma / (WAX.rho20 * gRed));
  const startHeight = 1.35 * predicted;
  const R0 = Math.sqrt(volume / (Math.PI * startHeight));
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
    // Film OFF. This tool calibrates the COHESION coefficient against a
    // continuum result; the anti-coalescence film is a separate mechanism with
    // its own parameter. Leaving it on, the settling puddle sheds fragments, each
    // fragment takes a fresh blob id, and the film then holds them apart -- the
    // puddle came back as 71 pieces and the "thickness" was measured across a
    // shattered sheet. Adding a mechanism quietly changed what the instrument
    // measured, which is the kind of thing that only shows up if you re-run the
    // calibration after every model change.
    disjoin: 0,
    wetting: 0,             // non-wetting analytic puddle; plate adhesion is separate
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
console.log(`puddle ${(n * SPACING3 * 1e6).toFixed(1)} mL at the shipping spacing ` +
  `${(Math.cbrt(SPACING3) * 1e3).toFixed(2)} mm`);
console.log('');
console.log('  g_red    predicted e   measured e   drift  layers   sigma_eff mN/m  clamps blobs  use');

const rho = WAX.rho20;
const V = n * SPACING3;
const rows = [];
const MIN_LAYERS = 4.5;
// Keep the entire sweep comfortably above the 4.5-layer discretisation floor
// at the shipping spacing. Scale gravity with the requested tension so the
// analytic puddle thickness (and therefore its resolution) stays fixed when a
// surfactant-rich default is selected. The old 0.08 point sat on that floor and
// made the exponent partly a measurement of particle size.
const gravityScale = IFACE.sigma / 4.5e-3;
for (const gRed of [0.030, 0.050, 0.075].map((g) => g * gravityScale)) {
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
