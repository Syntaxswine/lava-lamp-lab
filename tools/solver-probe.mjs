// solver-probe.mjs — look at the density projection itself, step by step, from
// the cold lattice. If the lamp explodes the instant the mechanics engage, the
// evidence is here: the density ratio the solver sees and the correction it
// applies, per particle, per iteration.
import { Lamp } from '../js/lamp.js';

const n = Number(process.argv[2] ?? 2000);
const steps = Number(process.argv[3] ?? 24);
const lamp = new Lamp({ particles: n });
lamp.reset();
const w = lamp.wax;

console.log(`n=${n}  dx=${(w.dx * 1e3).toFixed(2)}mm  h=${(w.h * 1e3).toFixed(2)}mm  ` +
  `rho0=${w.rho0.toFixed(1)}  eps=${w.eps.toExponential(2)}`);

// density of the untouched lattice
const stat = (arr) => {
  let lo = Infinity, hi = -Infinity, s = 0;
  for (let i = 0; i < arr.length; i++) { const v = arr[i]; if (v < lo) lo = v; if (v > hi) hi = v; s += v; }
  return { lo, hi, mean: s / arr.length };
};
let nb = 0;
for (let i = 0; i < w.n; i++) nb += w.nbrCount[i];
console.log(`lattice: mean neighbours ${(nb / w.n).toFixed(1)}`);

const dt = 1 / 60;
console.log('step   maxSpeed  maxDisp   rho/rho0 lo..hi   clamps  com');
for (let s = 1; s <= steps; s++) {
  const x0 = Float32Array.from(w.y);
  w.step(dt, lamp.column, lamp.env);
  let vmax = 0, dmax = 0;
  for (let i = 0; i < w.n; i++) {
    const sp = Math.hypot(w.vx[i], w.vy[i], w.vz[i]);
    if (sp > vmax) vmax = sp;
    const d = Math.abs(w.y[i] - x0[i]);
    if (d > dmax) dmax = d;
  }
  const r = stat(w.rho);
  console.log(`${String(s).padStart(4)}  ${(vmax * 100).toFixed(2).padStart(9)} cm/s ` +
    `${(dmax * 1000).toFixed(3).padStart(8)} mm  ` +
    `${(r.lo / w.rho0).toFixed(3)}..${(r.hi / w.rho0).toFixed(3)}  ` +
    `${String(w.clamps).padStart(6)}  ${(w.centreOfMass * 100).toFixed(2)} cm`);
}
