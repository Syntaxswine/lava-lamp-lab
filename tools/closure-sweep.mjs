import { Column } from '../js/column.js';
import { HEAT } from '../js/params.js';
const P = HEAT.bulbW * HEAT.fCouple;
console.log(' nuC   plate   base   bulk    top   grad  t63/min');
for (const nuC of [0.55, 0.8, 1.2, 1.8, 2.6, 4.0]) {
  HEAT.nuC = nuC;
  const col = new Column(96); col.reset(HEAT.Tamb);
  let t63 = NaN;
  for (let i = 0; i < 600 * 60; i++) { col.step(1, P, HEAT.Tamb); }
  const fin = col.Tbulk;
  const col2 = new Column(96); col2.reset(HEAT.Tamb);
  for (let i = 0; i < 600 * 60; i++) {
    col2.step(1, P, HEAT.Tamb);
    if (isNaN(t63) && col2.Tbulk - HEAT.Tamb > 0.63 * (fin - HEAT.Tamb)) t63 = i / 60;
  }
  console.log(`${nuC.toFixed(2).padStart(4)}  ` +
    [col.Tplate, col.Tbase, col.Tbulk, col.Ttop, col.Tbase - col.Ttop]
      .map(v => v.toFixed(1).padStart(6)).join(' ') + `  ${t63.toFixed(0).padStart(6)}`);
}
