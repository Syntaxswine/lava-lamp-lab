// Focused deformation/pinch-off probe. This is deliberately separate from the
// thermal lamp probe: it answers the visible question -- do bodies elongate and
// does one connected body become multiple macroscopic daughters?
//
// usage: node tools/shape-probe.mjs [sigma-mN/m] [wax-mPa.s] [wax-mL]
//                                    [seconds] [curvature-weight]

// Example:
//   node tools/shape-probe.mjs 2.5 15 75 45 0.4

import { Lamp } from '../js/lamp.js';
import { WAX, IFACE, SOLVER } from '../js/params.js';

function neckProfile(w, bid = 1) {
  const ys = [];
  for (let i = 0; i < w.n; i++) if (w.id[i] === bid) ys.push(w.y[i]);
  if (ys.length < 80) return null;
  const y0 = Math.min(...ys), y1 = Math.max(...ys), span = Math.max(y1 - y0, w.dx);
  const bins = new Int32Array(16);
  for (const y of ys) bins[Math.min(15, Math.floor(16 * (y - y0) / span))]++;
  let best = null;
  for (let b = 2; b <= 13; b++) {
    let lp = 0, rp = 0;
    for (let q = 0; q < b; q++) lp = Math.max(lp, bins[q]);
    for (let q = b + 1; q < 16; q++) rp = Math.max(rp, bins[q]);
    const shoulder = Math.min(lp, rp);
    if (!shoulder || !bins[b]) continue;
    const ratio = bins[b] / shoulder;
    const radiusDx = Math.sqrt(bins[b] * w.Vp / (Math.PI * span / 16)) / w.dx;
    if (!best || ratio < best.ratio) best = { bin: b, count: bins[b], ratio, radiusDx };
  }
  return best && {
    ...best,
    ratio: Number(best.ratio.toFixed(3)),
    radiusDx: Number(best.radiusDx.toFixed(3)),
    bins: [...bins],
  };
}

const sigma = Number(process.argv[2] ?? IFACE.sigma * 1e3) * 1e-3;
WAX.mu = Number(process.argv[3] ?? WAX.mu * 1e3) * 1e-3;
const volume = Number(process.argv[4] ?? WAX.volume * 1e6) * 1e-6;
const seconds = Number(process.argv[5] ?? 45);
SOLVER.curvK = Number(process.argv[6] ?? SOLVER.curvK);

const lamp = new Lamp({ particles: SOLVER.particles, waxVolume: volume });
lamp.env.sigma = sigma;
lamp.warmStart();

let minBlobs = Infinity;
let maxBlobs = 0;
let peak = { stretch: 0, time: 0, blobs: [] };
let postTenPeak = 0;
let splitEvents = 0;
let firstSplitAt = null;
const initial = lamp.diagnostics();
const initialBodies = lamp.wax.blobStats(8)
  .filter((b) => b.count >= 8)
  .map((b) => ({ count: b.count, stretch: Number(b.verticalStretch.toFixed(2)) }));
const initialNeck = neckProfile(lamp.wax);
let previous = initial.blobs;
const dt = lamp.physDt();
const steps = Math.ceil(seconds / dt);

for (let k = 0; k < steps; k++) {
  lamp.substep(dt);
  if (k % Math.max(1, Math.round(0.25 / dt)) !== 0 && k !== steps - 1) continue;
  const d = lamp.diagnostics();
  minBlobs = Math.min(minBlobs, d.blobs);
  maxBlobs = Math.max(maxBlobs, d.blobs);
  if (d.blobs > previous) {
    splitEvents += d.blobs - previous;
    if (firstSplitAt === null) firstSplitAt = lamp.lampTime - 12000;
  }
  previous = d.blobs;
  if (d.maxStretch > peak.stretch) {
    peak = {
      stretch: d.maxStretch,
      time: lamp.lampTime,
      blobs: lamp.wax.blobStats(8)
        .filter((b) => b.count >= 8)
        .map((b) => ({ count: b.count, stretch: Number(b.verticalStretch.toFixed(2)) })),
    };
  }
  if (lamp.lampTime - 12000 >= 10) postTenPeak = Math.max(postTenPeak, d.maxStretch);
}

const final = lamp.diagnostics();
const finalBodies = lamp.wax.blobStats(8)
  .filter((b) => b.count >= 8)
  .map((b) => ({ count: b.count, stretch: Number(b.verticalStretch.toFixed(2)) }));
console.log(JSON.stringify({
  sigma_mNm: sigma * 1e3,
  wax_mPas: WAX.mu * 1e3,
  wax_mL: volume * 1e6,
  curvK: SOLVER.curvK,
  simulated_s: seconds,
  particle_spacing_mm: lamp.wax.dx * 1e3,
  capillary_dt_ms: dt * 1e3,
  initial_visible_blobs: initial.blobs,
  initial_bodies: initialBodies,
  initial_neck: initialNeck,
  visible_blob_range: [minBlobs, maxBlobs],
  split_events: splitEvents,
  first_split_at_s: firstSplitAt === null ? null : Number(firstSplitAt.toFixed(2)),
  resolved_neck_pinches: final.pinches,
  peak_stretch: Number(peak.stretch.toFixed(3)),
  peak_stretch_after_10s: Number(postTenPeak.toFixed(3)),
  peak_at_s_after_handover: Number((peak.time - 12000).toFixed(2)),
  peak_bodies: peak.blobs,
  final_visible_blobs: final.blobs,
  final_bodies: finalBodies,
  final_droplets: final.droplets,
  clamps: final.clamps,
  visc_clamped: final.viscClamped,
}, null, 2));
