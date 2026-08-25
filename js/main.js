// main.js — wire the lamp to the screen.
//
// The solver runs on a capillary-CFL step off an accumulator, not on whatever
// the frame time happened to be. Two reasons. A position-based fluid changes its own
// stiffness with dt, so a variable step makes the wax subtly firmer on fast
// frames; and a lava lamp moves at a centimetre a second, so thirty solver steps
// a second is far more than the motion needs while sixty would spend half the
// frame budget on nothing anyone can see.
//
// The step count per frame is bounded, so on a slow frame the lamp falls behind
// real time rather than taking a huge step and coming apart. The HUD clock is
// the lamp's own clock, so any gap between it and your wristwatch is visible
// rather than hidden.

import { Lamp } from './lamp.js';
import { Renderer, PALETTES } from './render.js';
import { WAX, AQ, IFACE, SOLVER } from './params.js';

const $ = (id) => document.getElementById(id);
const canvas = $('lamp-canvas');
const viewport = $('viewport');
const status = $('renderer-status');

let renderer;
try {
  renderer = new Renderer(canvas);
} catch (err) {
  viewport.classList.add('gl-unavailable');
  viewport.dataset.message =
    'This lamp needs WebGL2, which this browser did not provide. ' +
    'The physics runs headless too: see tools/lamp-probe.mjs in the repo.';
  canvas.style.display = 'none';
  status.textContent = 'WebGL2 unavailable';
  document.querySelector('.status-dot').classList.add('bad');
  throw err;
}

const lamp = new Lamp({ particles: SOLVER.particles });
lamp.env.iterations = SOLVER.iterations;
lamp.warmStart();
status.textContent = 'WebGL2 · developed warm start';

// ---------------------------------------------------------------------------
// camera
// ---------------------------------------------------------------------------
let drag = null;
canvas.addEventListener('pointerdown', (e) => {
  drag = { x: e.clientX, y: e.clientY };
  canvas.setPointerCapture(e.pointerId);
});
canvas.addEventListener('pointermove', (e) => {
  if (!drag) return;
  renderer.orbit.az -= (e.clientX - drag.x) * 0.006;
  renderer.orbit.el = Math.max(-0.5, Math.min(0.9,
    renderer.orbit.el + (e.clientY - drag.y) * 0.004));
  drag = { x: e.clientX, y: e.clientY };
});
const endDrag = () => { drag = null; };
canvas.addEventListener('pointerup', endDrag);
canvas.addEventListener('pointercancel', endDrag);
canvas.addEventListener('wheel', (e) => {
  e.preventDefault();
  renderer.orbit.dist = Math.max(0.24, Math.min(1.4,
    renderer.orbit.dist * (1 + Math.sign(e.deltaY) * 0.08)));
}, { passive: false });
$('reset-view').addEventListener('click', () => {
  renderer.orbit = { az: 0.6, el: 0.12, dist: 0.60 };
});

// ---------------------------------------------------------------------------
// controls
// ---------------------------------------------------------------------------
let running = true;
let timeScale = 1;
const adjustableControls = [...document.querySelectorAll('.control-panel input, .control-panel select')];
const defaultControlValues = new Map(adjustableControls.map((el) => [el, el.value]));

const bind = (id, out, fmt, apply) => {
  const el = $(id);
  const o = $(out);
  const push = () => {
    const v = Number(el.value);
    o.textContent = fmt(v);
    apply(v);
  };
  el.addEventListener('input', push);
  push();
};

bind('bulb-input', 'bulb-output', (v) => `${v} W`, (v) => { lamp.env.bulbW = v; });
bind('room-input', 'room-output', (v) => `${v.toFixed(1)} °C`, (v) => { lamp.env.Tamb = v; });
bind('drho-input', 'drho-output', (v) => `${v.toFixed(1)} kg/m³`, (v) => { lamp.env.dRho = v; });
bind('sigma-input', 'sigma-output', (v) => `${v.toFixed(1)} mN/m`,
  (v) => { lamp.env.sigma = v / 1000; });
bind('visc-input', 'visc-output', (v) => `${v} mPa·s`, (v) => { lamp.env.muAq = v / 1000; });

let pendingWax = null, pendingN = null;
bind('wax-input', 'wax-output', (v) => `${v} mL`, (v) => { pendingWax = v / 1e6; });
bind('particles-input', 'particles-output', (v) => v.toLocaleString(),
  (v) => { pendingN = v; });

$('solver-input').addEventListener('change', (e) => {
  lamp.env.iterations = Number(e.target.value);
});
$('time-input').addEventListener('change', (e) => { timeScale = Number(e.target.value); });
$('view-input').addEventListener('change', (e) => {
  renderer.view = Number(e.target.value);
  updateKey();
});
$('palette-input').addEventListener('change', (e) => {
  renderer.palette = PALETTES[e.target.value];
  updateKey();
});

function updateKey() {
  const bar = $('key-bar'), lo = $('key-lo'), hi = $('key-hi');
  if (renderer.view === 1) {
    bar.style.background = 'linear-gradient(90deg,#2f6bd8,#63d1c4,#f0d878,#ff4a24)';
    lo.textContent = '18 °C'; hi.textContent = '78 °C';
  } else if (renderer.view === 2) {
    bar.style.background = 'linear-gradient(90deg,#2669f2,#8f9ec4,#ff731e)';
    lo.textContent = 'SINKING'; hi.textContent = 'RISING';
  } else {
    const p = renderer.palette;
    const rgb = (c) => `rgb(${c.map((x) => Math.round(x * 255)).join(',')})`;
    bar.style.background = `linear-gradient(90deg,${rgb(p.deep)},${rgb(p.wax)},${rgb(p.bulb)})`;
    lo.textContent = 'SHADOW'; hi.textContent = 'LIT';
  }
}
updateKey();

$('toggle').addEventListener('click', (e) => {
  running = !running;
  e.target.textContent = running ? 'Pause' : 'Resume';
});

$('restore-defaults').addEventListener('click', (e) => {
  const button = e.currentTarget;
  button.disabled = true;
  status.textContent = 'restoring defaults…';

  for (const [el, value] of defaultControlValues) {
    el.value = value;
    el.dispatchEvent(new Event(el instanceof HTMLSelectElement ? 'change' : 'input'));
  }
  running = true;
  $('toggle').textContent = 'Pause';

  // Let the status paint before rebuilding the default developed state. Wax
  // charge and particle count are pending controls, so apply those first.
  setTimeout(() => {
    applyPending(true);
    lamp.warmStart();
    acc = 0;
    physDt = lamp.physDt();
    status.textContent = 'defaults restored · developed warm start';
    button.disabled = false;
  }, 0);
});

$('cold').addEventListener('click', () => {
  applyPending(true);
  lamp.startedWarm = false;
  lamp.reset();
  status.textContent = 'cold start · 22 °C solid plug';
});
$('warm').addEventListener('click', () => {
  applyPending(true);
  status.textContent = 'building steady state…';
  setTimeout(() => {
    lamp.warmStart();
    status.textContent = 'warm start · developed circulation';
  }, 0);
});
$('ff').addEventListener('click', () => {
  const r = lamp.fastForward();
  status.textContent = r.reason === 'liftoff'
    ? `skipped ${(r.skipped / 60).toFixed(0)} min — wax is buoyant`
    : `skipped ${(r.skipped / 60).toFixed(0)} min — still nothing lifting`;
});
$('shake').addEventListener('click', () => {
  lamp.wax.shake();
  status.textContent = 'shaken — that is how lamps get cloudy';
});

function applyPending(force = false) {
  if (pendingN !== null && (force || pendingN !== lamp.wax.n)) {
    lamp.wax = new (lamp.wax.constructor)(pendingN, pendingWax ?? lamp.wax.volume);
    lamp.wax.reset(lamp.env.Tamb);
  } else if (pendingWax !== null && (force || pendingWax !== lamp.wax.volume)) {
    lamp.setWaxVolume(pendingWax);
  }
  pendingN = null; pendingWax = null;
}

// A change of particle count or wax charge cannot be applied mid-flight without
// throwing away the state, so it waits for the next start instead of silently
// doing nothing.
for (const id of ['wax-input', 'particles-input']) {
  $(id).addEventListener('change', () => {
    status.textContent = 'press Warm start or Cold start to apply';
  });
}

// ---------------------------------------------------------------------------
// loop
// ---------------------------------------------------------------------------
// The solver step is not a chosen number: it comes from the capillary-wave CFL
// for the smoothing length and interfacial tension currently in force, so moving
// the tension slider moves the timestep with it. See Wax.capillaryDt.
let physDt = lamp.physDt();
let acc = 0;
let last = performance.now();
let fpsAcc = 0, fpsN = 0, fps = 60;
let hudAt = 0;
let frameCost = 16;
let drewOnce = false;

function fmtClock(s) {
  const h = Math.floor(s / 3600), m = Math.floor((s % 3600) / 60);
  return h > 0 ? `${h}:${String(m).padStart(2, '0')}` : `0:${String(m).padStart(2, '0')}`;
}

function frame(now) {
  const dtWall = Math.min(0.05, (now - last) / 1000);
  last = now;
  const t0 = performance.now();

  let steps = 0;
  if (running) {
    physDt = lamp.physDt();
    acc = Math.min(acc + dtWall * timeScale, 8 * physDt);
    const budget = timeScale > 1 ? 5 : 3;
    while (acc >= physDt && steps < budget) {
      lamp.substep(physDt);
      acc -= physDt;
      steps++;
    }
  }
  // nothing moved, nothing to redraw
  if (steps > 0 || !drewOnce) { renderer.draw(lamp.wax, lamp); drewOnce = true; }

  const cost = performance.now() - t0;
  frameCost += (cost - frameCost) * 0.1;
  // keep the frame inside a 60 Hz budget by moving pixels, not physics
  if (frameCost > 26 && renderer.scale > 0.45) renderer.scale -= 0.01;
  else if (frameCost < 16 && renderer.scale < 1.0) renderer.scale += 0.005;

  fpsAcc += dtWall; fpsN++;
  if (fpsAcc > 0.5) { fps = fpsN / fpsAcc; fpsAcc = 0; fpsN = 0; }

  if (now - hudAt > 220) {
    hudAt = now;
    const d = lamp.diagnostics();
    $('hud-clock').textContent = fmtClock(d.lampTime);
    $('hud-blobs').textContent = d.blobs;
    $('hud-rise').textContent = (d.riseSpeed * 100).toFixed(1);
    $('hud-fps').textContent = fps.toFixed(0);
    const win = $('hud-window');
    win.textContent = d.molten < 0.02 ? 'wax still solid' : d.window;
    win.classList.toggle('stalled', d.window !== 'cycling');
    $('r-plate').textContent = d.Tplate.toFixed(1);
    $('r-base').textContent = d.Tbase.toFixed(1);
    $('r-bulk').textContent = d.Tbulk.toFixed(1);
    $('r-top').textContent = d.Ttop.toFixed(1);
    $('r-cross').textContent = `${d.crossover.toFixed(1)} °C`;
    $('r-rise').textContent = `${d.riseAtBase.toFixed(1)} °C`;
    $('r-sink').textContent = `${d.sinkAtTop.toFixed(1)} °C`;
    $('r-blob').textContent = d.meanBlobRadius > 0
      ? `${(d.meanBlobRadius * 1000).toFixed(1)} mm` : '—';
    const res = d.resolution;
    $('r-cap').textContent = res.known ? `${(res.capillary * 1000).toFixed(1)} mm` : '—';
    $('r-molten').textContent = `${(d.molten * 100).toFixed(0)}%`;
    $('r-res').textContent = !res.known
      ? `blob scale not measurable yet — nothing in the tank is buoyant`
      : res.resolved
        ? `resolved: ${(res.capillary * 1000).toFixed(1)} mm capillary blobs vs a ${(res.floor * 1000).toFixed(1)} mm floor`
        : `UNDER-RESOLVED ${(1 / res.ratio).toFixed(1)}x — blobs come out bigger, hotter and faster than these settings imply`;
    $('r-res').classList.toggle('warn', res.known && !res.resolved);
    $('r-dt').textContent = `${(1 / d.physDt).toFixed(0)} Hz`;
    $('r-clamp').textContent = d.clamps > 0 ? d.clamps.toLocaleString() : 'none';
    $('r-clamp').classList.toggle('warn', d.clamps > 200);
  }
  requestAnimationFrame(frame);
}
requestAnimationFrame(frame);

// expose for the console and for browser-side checks
window.LAMP = lamp;
window.RENDERER = renderer;
