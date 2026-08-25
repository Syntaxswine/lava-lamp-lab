// skin-probe.mjs — measure the two visible consequences of the interfacial
// model: hot-base wetting footprint and temperature-dependent film drainage.
//
// usage: node tools/skin-probe.mjs [seconds] [plate-wet-K]

import { Lamp } from '../js/lamp.js';
import { GEO, IFACE, SOLVER, globeRadius } from '../js/params.js';

const seconds = Number(process.argv[2] ?? 30);
if (process.argv[3]) IFACE.plateWetK = Number(process.argv[3]);
const lamp = new Lamp({ particles: SOLVER.particles });
lamp.warmStart();
const w = lamp.wax;

function baseFootprint() {
  const radii = [];
  const yCut = 0.065 * GEO.H;
  for (let i = 0; i < w.n; i++) {
    if (w.y[i] < yCut) radii.push(Math.hypot(w.x[i], w.z[i]));
  }
  radii.sort((a, b) => a - b);
  const r95 = radii.length ? radii[Math.floor(0.95 * (radii.length - 1))] : 0;
  const available = globeRadius(w.rp) - w.rp;
  return {
    particles: radii.length,
    r95_mm: Number((r95 * 1e3).toFixed(2)),
    coverage: Number((r95 / available).toFixed(3)),
  };
}

const skin = [40, 44, 48].map((T) => {
  const rate = w.filmDrainRate(T);
  return {
    T_C: T,
    sigma_mNm: Number((IFACE.sigma * w.surfaceScale(T) * 1e3).toFixed(3)),
    drainage_rate: Number(rate.toFixed(3)),
    effective_drain_s: Number((SOLVER.drainTime / rate).toFixed(1)),
    film_barrier: Number((w.surfaceScale(T) / Math.sqrt(rate)).toFixed(3)),
  };
});

const trace = [{ t_s: 0, ...baseFootprint() }];
const dt = lamp.physDt();
const sampleSteps = Math.max(1, Math.round(5 / dt));
let peakStretch = lamp.diagnostics().maxStretch;
for (let k = 1; k <= Math.ceil(seconds / dt); k++) {
  lamp.substep(dt);
  peakStretch = Math.max(peakStretch, lamp.diagnostics().maxStretch);
  if (k % sampleSteps === 0 || k === Math.ceil(seconds / dt)) {
    trace.push({ t_s: Number((k * dt).toFixed(1)), ...baseFootprint() });
  }
}

const end = lamp.diagnostics();
const result = {
  plate_wetting: trace,
  plate_wet_K: IFACE.plateWetK,
  skin,
  warm_to_cool_drain_time_ratio: Number((skin[2].effective_drain_s /
    skin[0].effective_drain_s).toFixed(3)),
  visible_blobs: end.blobs,
  coalescences: end.coalescences,
  peak_stretch: Number(peakStretch.toFixed(3)),
  clamps: end.clamps,
};
console.log(JSON.stringify(result, null, 2));

const earlyCoverage = Math.min(...trace.filter((p) => p.t_s <= 10).map((p) => p.coverage));
if (earlyCoverage < 0.85) {
  console.error('FAIL: developed hot pool does not sustain 85% of the plate radius');
  process.exitCode = 2;
}
if (result.warm_to_cool_drain_time_ratio >= 0.5) {
  console.error('FAIL: the warm film is not at least twice as quick to drain');
  process.exitCode = 3;
}
if (end.clamps !== 0) {
  console.error('FAIL: wetting/skin model leaned on the velocity clamp');
  process.exitCode = 4;
}
