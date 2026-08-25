// regression-probe.mjs — deterministic smoke test for the shipped lamp.
//
// It runs the same developed-state hand-off twice, proves the states are bitwise
// identical, and compares a state/configuration fingerprint with a committed
// baseline.  Deliberate constant changes therefore fail loudly instead of
// turning into a new anecdotal screenshot.

import fs from 'node:fs';
import assert from 'node:assert/strict';
import { Lamp } from '../js/lamp.js';
import { WAX, IFACE, SOLVER } from '../js/params.js';

const baseline = JSON.parse(fs.readFileSync(new URL('./lamp-baseline.json', import.meta.url), 'utf8'));
const mutation = process.argv.find((a) => a.startsWith('--mutate-coh='));
if (mutation) SOLVER.cohK += Number(mutation.split('=')[1]);

function fnv(bytes, h = 0x811c9dc5) {
  for (const b of bytes) h = Math.imul(h ^ b, 0x01000193) >>> 0;
  return h;
}

function fingerprint(lamp) {
  let h = 0x811c9dc5;
  const arrays = [lamp.wax.x, lamp.wax.y, lamp.wax.z, lamp.wax.T, lamp.column.T];
  for (const a of arrays) h = fnv(new Uint8Array(a.buffer, a.byteOffset, a.byteLength), h);
  const config = new TextEncoder().encode(JSON.stringify({
    cohK: SOLVER.cohK,
    cohRefDx: SOLVER.cohRefDx,
    curvK: SOLVER.curvK,
    waxMu: WAX.mu,
    drainTime: SOLVER.drainTime,
    disjoinK: SOLVER.disjoinK,
    pinchStretch: SOLVER.pinchStretch,
    pinchRatio: SOLVER.pinchRatio,
    pinchRadiusDx: SOLVER.pinchRadiusDx,
    pinchDelay: SOLVER.pinchDelay,
    sigma: IFACE.sigma,
    sigmaTempCoeff: IFACE.sigmaTempCoeff,
    filmRefT: IFACE.filmRefT,
    filmTempScale: IFACE.filmTempScale,
    plateWetK: IFACE.plateWetK,
    plateWetRangeH: IFACE.plateWetRangeH,
    dt: lamp.physDt(),
  }));
  return fnv(config, h).toString(16).padStart(8, '0');
}

function run() {
  const lamp = new Lamp({ particles: 1800, seed: baseline.seed });
  const start = lamp.warmStart();
  assert.equal(start.blobs, 4, 'developed warm start must contain four visible wax bodies');
  const dt = lamp.physDt();
  const steps = Math.round(baseline.seconds / dt);
  for (let i = 0; i < steps; i++) lamp.substep(dt);
  const end = lamp.diagnostics();
  assert.equal(end.clamps, 0, 'shipping step must not lean on the velocity clamp');
  assert.ok(end.blobs >= 4, 'surfactant film must preserve the developed blob population');
  return { fingerprint: fingerprint(lamp), end };
}

const a = run();
const b = run();
assert.equal(a.fingerprint, b.fingerprint, 'fixed seed produced two different runs');

function shakeFingerprint() {
  const lamp = new Lamp({ particles: 1800, seed: baseline.seed });
  lamp.warmStart();
  lamp.wax.shake();
  assert.ok(new Set(lamp.wax.id).size > 4, 'shake must give droplets fresh film identities');
  let h = fnv(new Uint8Array(lamp.wax.id.buffer));
  h = fnv(new Uint8Array(lamp.wax.vx.buffer), h);
  h = fnv(new Uint8Array(lamp.wax.vy.buffer), h);
  return h.toString(16).padStart(8, '0');
}
assert.equal(shakeFingerprint(), shakeFingerprint(), 'seeded shake produced two different states');

console.log(`fingerprint ${a.fingerprint}   seed ${baseline.seed}   ${baseline.seconds}s`);
console.log(`blobs ${a.end.blobs}   clamps ${a.end.clamps}   ` +
  `rise ${(a.end.riseSpeed * 100).toFixed(2)} cm/s   com ${(a.end.com * 100).toFixed(1)}%`);
if (baseline.fingerprint === 'PENDING') {
  console.log('baseline pending; record the fingerprint after reviewing this run');
  process.exit(2);
}
assert.equal(a.fingerprint, baseline.fingerprint,
  'lamp state/configuration differs from the reviewed baseline');
console.log('determinism baseline: PASS');
