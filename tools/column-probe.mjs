// column-probe.mjs — warm-up curve and steady state of the aqueous column,
// with no wax present. Checks the thermal plumbing against measured lamps:
// base plate ~58 C, bulk ~45 C, top ~40 C for a 40 W bulb in a 22 C room.
import { Column } from '../js/column.js';
import { HEAT, T_CROSS, WAX, AQ } from '../js/params.js';

const bulb = Number(process.argv[2] ?? HEAT.bulbW);
const Tamb = Number(process.argv[3] ?? HEAT.Tamb);
const col = new Column(96);
col.reset(Tamb);
const P = bulb * HEAT.fCouple;

console.log(`bulb ${bulb} W  (${P.toFixed(1)} W into the globe)   room ${Tamb} C`);
console.log('  t/min   plate    base    bulk     top   loss W');
let t = 0;
const marks = new Set([0, 5, 10, 20, 30, 45, 60, 90, 120, 180, 240]);
for (let step = 0; step <= 240 * 60; step++) {
  if (marks.has(Math.round(t / 60)) && Math.abs(t % 60) < 0.5) {
    console.log(`  ${String(Math.round(t / 60)).padStart(5)}  ` +
      [col.Tplate, col.Tbase, col.Tbulk, col.Ttop].map(v => v.toFixed(1).padStart(6)).join('  ') +
      `   ${col.lossWatts(Tamb).toFixed(1).padStart(5)}`);
    marks.delete(Math.round(t / 60));
  }
  col.step(1, P, Tamb);
  t += 1;
}
console.log('');
console.log(`steady gradient base->top : ${(col.Tbase - col.Ttop).toFixed(1)} K`);
console.log(`crossover temperature     : ${T_CROSS.toFixed(1)} C`);
console.log(`wax melting point         : ${WAX.Tmelt.toFixed(1)} C`);
const window_ok = col.Ttop < T_CROSS && col.Tplate > T_CROSS && col.Ttop > WAX.Tmelt;
console.log(`cycling window            : ${window_ok ? 'OPEN' : 'CLOSED'}` +
  `  (need Tmelt ${WAX.Tmelt} < Ttop ${col.Ttop.toFixed(1)} < Tcross ${T_CROSS.toFixed(1)} < Tplate ${col.Tplate.toFixed(1)})`);
