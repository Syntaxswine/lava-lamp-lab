// lamp-probe.mjs — run the whole lamp headless and ask whether it behaves like
// a lava lamp. Reports the warm-up, the thermal window, and then the mechanical
// signature: how many blobs, how big, how fast they rise, how long a cycle takes.
//
// usage: node tools/lamp-probe.mjs [minutes] [particles] [bulbW] [roomC] [warm]
//   pass "warm" as the fifth argument to start from the developed steady state
//   instead of a cold plug.

import { Lamp } from '../js/lamp.js';
import { WAX, GEO, IFACE, HEAT } from '../js/params.js';

const minutes = Number(process.argv[2] ?? 6);
const particles = Number(process.argv[3] ?? 3000);
const bulbW = Number(process.argv[4] ?? 40);
const roomC = Number(process.argv[5] ?? 22);
const warm = process.argv[6] === 'warm' || process.argv[2] === 'warm';

const lamp = new Lamp({ particles });
lamp.env.bulbW = bulbW;
lamp.env.Tamb = roomC;
lamp.reset();

const w = lamp.wax;
console.log(`particles ${w.n}   wax ${(w.volume * 1e6).toFixed(0)} mL   ` +
  `particle r ${(w.rp * 1000).toFixed(2)} mm   spacing ${(w.dx * 1000).toFixed(2)} mm   ` +
  `h ${(w.h * 1000).toFixed(2)} mm`);
console.log(`rest density ${w.rho0.toFixed(1)} kg/m3 (bulk ${WAX.rho20})   ` +
  `plug height ${(w.fillHeight * 100).toFixed(1)} cm of ${(GEO.H * 100).toFixed(0)} cm`);
console.log('');

// ---- start: cold plug plus a gated skip, or straight to the steady state ----
const t0 = Date.now();
let ff;
if (warm) {
  lamp.warmStart();
  ff = { skipped: lamp.lampTime, reason: 'warm start (initial condition, not a history)' };
  console.log(`warm start: column relaxed to steady state, wax molten at the ` +
    `local fluid temperature`);
} else {
  ff = lamp.fastForward(6 * 3600, 5);
  console.log(`fast-forward: skipped ${(ff.skipped / 60).toFixed(1)} min of lamp time, ` +
    `stopped because ${ff.reason}`);
}
const d0 = lamp.diagnostics();
console.log(`  at hand-over: plate ${d0.Tplate.toFixed(1)}  base ${d0.Tbase.toFixed(1)}  ` +
  `bulk ${d0.Tbulk.toFixed(1)}  top ${d0.Ttop.toFixed(1)}   wax ${d0.waxT.toFixed(1)}` +
  `  molten ${(d0.molten * 100).toFixed(0)}%`);
console.log(`  crossover ${d0.crossover.toFixed(1)} C   ` +
  `must reach ${d0.riseAtBase.toFixed(1)} C to lift at the base, ` +
  `fall below ${d0.sinkAtTop.toFixed(1)} C to sink at the top`);
console.log(`  window: ${d0.window}`);
console.log('');

// ---- coupled run ----------------------------------------------------------
const dt = lamp.physDt();
console.log(`solver step: 1/${(1 / dt).toFixed(0)} s from the capillary CFL`);
console.log('');
const steps = Math.round(minutes * 60 / dt);
const rows = [];
let peakRise = 0, peakRiseAt = 0;
const sampleEvery = Math.round(15 / dt);
for (let s = 1; s <= steps; s++) {
  lamp.substep(dt);
  if (w.maxRise > peakRise) { peakRise = w.maxRise; peakRiseAt = lamp.lampTime; }
  if (s % sampleEvery === 0) {
    const d = lamp.diagnostics();
    rows.push(d);
    console.log(`t+${((lamp.lampTime - ff.skipped) / 60).toFixed(1).padStart(5)} min  ` +
      `plate ${d.Tplate.toFixed(1)}  top ${d.Ttop.toFixed(1)}  ` +
      `wax ${d.waxT.toFixed(1)}/${d.waxTmax.toFixed(1)}  ` +
      `molten ${(d.molten * 100).toFixed(0).padStart(3)}%  ` +
      `blobs ${String(d.blobs).padStart(3)}  ` +
      `biggest ${String(d.biggest).padStart(4)}p  ` +
      `r ${(d.meanBlobRadius * 1000).toFixed(1).padStart(4)} mm  ` +
      `rise ${(d.riseSpeed * 100).toFixed(2).padStart(5)} cm/s  ` +
      `com ${(d.com * 100).toFixed(0).padStart(3)}%  ` +
      `clamp ${d.clamps}  ${d.window}`);
  }
}
const wall = (Date.now() - t0) / 1000;
console.log('');
const d = lamp.diagnostics();
console.log(`--- after ${minutes} min of coupled lamp time -------------------------`);
console.log(`blobs                 ${d.blobs}   largest ${d.biggest} particles ` +
  `(${(d.biggest * w.Vp * 1e6).toFixed(1)} mL)`);
console.log(`mean blob radius      ${(d.meanBlobRadius * 1000).toFixed(1)} mm`);
// A blob spends most of a cycle parked at one end. Taking the median over the
// whole run therefore measures how long it sits still, not how fast it travels;
// the comparable number is the median while it is actually in transit.
const moving = rows.map((r) => r.riseSpeed).filter((v) => v > 0.001).sort((a, b) => a - b);
const medRise = moving.length ? moving[moving.length >> 1] : 0;
console.log(`rise speed            ${(medRise * 100).toFixed(2)} cm/s median while in transit ` +
  `(${moving.length}/${rows.length} samples), peak ${(peakRise * 100).toFixed(2)} cm/s` +
  `   (real lamps: 1-3 cm/s)`);
const res = d.resolution;
if (!res.known) {
  console.log(`resolution            not measurable: no buoyant wax at the end of the run`);
} else {
  console.log(`resolution            capillary length ${(res.capillary * 1000).toFixed(1)} mm ` +
    `(from a measured drho of ${res.drho.toFixed(2)} kg/m3) vs floor ` +
    `${(res.floor * 1000).toFixed(1)} mm (3 particle spacings) -> ` +
    `${res.resolved ? 'RESOLVED' : 'UNDER-RESOLVED by ' + (1 / res.ratio).toFixed(1) + 'x'}`);
}
if (res.known && !res.resolved) {
  console.log(`                      blobs will come out larger, hotter and faster than`);
  console.log(`                      these parameters imply. Raise the particle count to`);
  console.log(`                      about ${Math.ceil(lamp.wax.volume / (res.capillary / 3) ** 3 / 100) * 100} to resolve them.`);
}
console.log(`cycle period          ${d.cyclePeriod ? d.cyclePeriod.toFixed(1) + ' s (spectral)' : 'not resolved -- refusing to guess'}`);
console.log(`first melt            ${d.firstMelt === null ? 'never' : (d.firstMelt / 60).toFixed(1) + ' min'}`);
console.log(`first lift            ${d.firstRise === null ? 'never'
  : ((d.firstRise - ff.skipped) / 60).toFixed(1) + ' min after hand-over'}`);
console.log(`window                ${d.window}`);
console.log(`energy in / out       ${d.inW.toFixed(1)} W / ${d.lossW.toFixed(1)} W`);
console.log('');
console.log(`compute: ${wall.toFixed(1)} s wall for ${(minutes * 60).toFixed(0)} s of lamp time ` +
  `=> ${(minutes * 60 / wall).toFixed(1)}x realtime at ${w.n} particles`);
console.log(`         ${(wall / steps * 1000).toFixed(2)} ms per substep`);
