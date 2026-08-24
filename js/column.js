// column.js — the aqueous phase, as a 1-D vertical temperature column.
//
// WHY 1-D. The lamp's interesting fluid is the wax; the water/glycol around it
// matters only as (a) the density the wax is weighed against and (b) the heat
// bath the wax exchanges with. Resolving it in 3-D would triple the cost to
// render something invisible. So it is a stack of well-mixed horizontal slices
// with a turbulent-convection closure on the vertical heat flux, plus honest
// two-way coupling: whatever heat the wax takes out of a slice, that slice loses.
//
// The closure is the standard Rayleigh-Benard scaling Nu = C * Ra^(1/3), applied
// on a mixing length l = min(distance to base, distance to cap, globe diameter)
// and only where the stratification is unstable (warm fluid under cold). That is
// what produces the real structure of a running lamp: a hot film on the base
// plate, a near-isothermal convecting bulk, and a cooler cap.

import { GEO, AQ, HEAT, globeArea, globePerim, globeRadius } from './params.js';

const NU = AQ.mu / AQ.rho20;            // kinematic viscosity, m^2/s
const ALPHA = AQ.k / (AQ.rho20 * AQ.c); // thermal diffusivity, m^2/s

export class Column {
  constructor(n = 96) {
    this.n = n;
    this.dy = GEO.H / n;
    this.y = new Float64Array(n);       // cell centre heights
    this.T = new Float64Array(n);       // degC
    this.cap = new Float64Array(n);     // J/K per cell
    this.area = new Float64Array(n);    // m^2 cross-section
    this.perim = new Float64Array(n);   // m
    this.qWax = new Float64Array(n);    // W deposited into this cell by wax
    this.phi = new Float64Array(n);     // wax volume fraction in this cell
    this.vAqCell = new Float64Array(n); // aqueous vertical velocity, m/s
    this.kEff = new Float64Array(n);    // W/m/K at the upper face (diagnostic)
    for (let i = 0; i < n; i++) {
      const y = (i + 0.5) * this.dy;
      this.y[i] = y;
      this.area[i] = globeArea(y);
      this.perim[i] = globePerim(y);
      this.cap[i] = AQ.rho20 * AQ.c * this.area[i] * this.dy;
    }
    this.Tplate = HEAT.Tamb;            // base plate + coil
    this.plateWaxQ = 0;                 // W drawn off the plate by the wax pool
    this.plateWaxFrac = 0;              // fraction of the plate the wax covers
    this.plateArea = globeArea(0);
    this.coilCell = Math.min(n - 1, Math.floor(GEO.coilY / this.dy));
    this.capArea = globeArea(GEO.H);
    this.reset(HEAT.Tamb);
  }

  reset(Tamb) {
    this.T.fill(Tamb);
    this.Tplate = Tamb;
    this.qWax.fill(0);
  }

  // Temperature of the fluid at an arbitrary height (linear interpolation).
  at(y) {
    const f = y / this.dy - 0.5;
    const i = Math.floor(f);
    if (i < 0) return this.T[0];
    if (i >= this.n - 1) return this.T[this.n - 1];
    const t = f - i;
    return this.T[i] * (1 - t) + this.T[i + 1] * t;
  }

  cellOf(y) {
    return Math.min(this.n - 1, Math.max(0, Math.floor(y / this.dy)));
  }

  // Heat the wax pulls out of / dumps into the fluid at height y, in watts.
  depositWaxHeat(y, watts) {
    this.qWax[this.cellOf(y)] += watts;
  }

  // -------------------------------------------------------------------------
  // Return flow.
  //
  // The globe is sealed and both phases are effectively incompressible, so the
  // net volume crossing any horizontal plane is zero: whatever the wax carries
  // up, the aqueous phase carries down through the gap left beside it.
  //
  //     Q_wax(y) + Q_aq(y) = 0        =>    v_aq = -Q_wax / A_aq
  //
  // Leaving this out is what let a single coalesced slug fill the tube and sit
  // there: it was rising through a fluid that never had to get out of the way.
  // With the constraint in place a blob that spans most of the cross-section
  // faces a return jet in the narrow annulus beside it, the drag on it goes up
  // sharply, and it can no longer behave like a piston in a frictionless bore.
  //
  // A_aq is floored at a fraction of the tube area. As the wax approaches a
  // perfect seal the true return velocity diverges, and a divergence in a
  // forward Euler step is not physics, it is a crash. The floor is stated here
  // rather than buried: below it, the model stops resolving the film between
  // wax and glass, which is where a real lamp gets its slow slug motion.
  // -------------------------------------------------------------------------
  updateReturnFlow(wax) {
    const n = this.n, dy = this.dy;
    this.phi.fill(0);
    const flux = this.vAqCell;
    flux.fill(0);
    for (let i = 0; i < wax.n; i++) {
      let c = (wax.y[i] / dy) | 0;
      if (c < 0) c = 0; else if (c >= n) c = n - 1;
      this.phi[c] += wax.Vp;
      flux[c] += wax.Vp * wax.vy[i];
    }
    const FLOOR = 0.12;                 // minimum open fraction of the bore
    for (let c = 0; c < n; c++) {
      const cellVol = this.area[c] * dy;
      const f = Math.min(0.95, this.phi[c] / cellVol);
      this.phi[c] = f;
      const aOpen = Math.max(FLOOR, 1 - f) * this.area[c];
      flux[c] = -(flux[c] / dy) / aOpen;
    }
    // smooth once: the flux profile is a histogram of a few thousand particles
    // and its cell-to-cell noise is not a velocity field
    const tmp = this._vsm || (this._vsm = new Float64Array(n));
    for (let c = 0; c < n; c++) {
      const a = flux[c > 0 ? c - 1 : 0], b = flux[c], d = flux[c < n - 1 ? c + 1 : n - 1];
      tmp[c] = 0.25 * a + 0.5 * b + 0.25 * d;
    }
    for (let c = 0; c < n; c++) flux[c] = tmp[c];
  }

  // Aqueous vertical velocity at height y, m/s (positive up).
  vAq(y) {
    const f = y / this.dy - 0.5;
    const i = Math.floor(f);
    if (i < 0) return this.vAqCell[0];
    if (i >= this.n - 1) return this.vAqCell[this.n - 1];
    const t = f - i;
    return this.vAqCell[i] * (1 - t) + this.vAqCell[i + 1] * t;
  }

  // Effective conductivity at the face between cells i and i+1.
  faceConductivity(i) {
    const grad = (this.T[i] - this.T[i + 1]) / this.dy;   // >0 == unstable
    if (grad <= 0) return AQ.k;
    const yf = (i + 1) * this.dy;
    const l = Math.min(yf, GEO.H - yf, 2 * globeRadius(yf));
    const Ra = HEAT.g * AQ.beta * grad * l ** 4 / (NU * ALPHA);
    const Nu = Math.max(1, HEAT.nuC * Math.cbrt(Ra));
    return AQ.k * Nu;
  }

  // Largest explicit-diffusion timestep this state can take, in seconds.
  maxStep() {
    let kmax = AQ.k;
    for (let i = 0; i < this.n - 1; i++) kmax = Math.max(kmax, this.faceConductivity(i));
    const a = kmax / (AQ.rho20 * AQ.c);
    return 0.4 * this.dy * this.dy / a;
  }

  // Advance dt seconds. bulbWatts is what the bulb delivers into the globe.
  //
  // The substep is chosen from the CURRENT state each time round, not once from
  // the state on entry. A sharp transient -- the first seconds after the bulb
  // comes on, say -- can make the entry-state limit tiny, and dividing the whole
  // interval by it demanded tens of thousands of substeps and simply hung. The
  // limit relaxes again as the field smooths, so re-reading it each pass costs
  // one sweep of the faces and bounds the work honestly.
  //
  // If the guard binds, `stiffSteps` records it rather than the solver quietly
  // taking an unstable step.
  step(dt, bulbWatts, Tamb) {
    const n = this.n;
    const flux = this._flux || (this._flux = new Float64Array(n + 1));
    let done = 0, guard = 0;
    while (done < dt - 1e-12 && guard++ < 4096) {
      const h = Math.min(this.maxStep(), dt - done);
      done += h;
      // interior faces
      for (let i = 0; i < n - 1; i++) {
        const k = this.faceConductivity(i);
        this.kEff[i] = k;
        const a = 0.5 * (this.area[i] + this.area[i + 1]);
        flux[i + 1] = k * a * (this.T[i] - this.T[i + 1]) / this.dy;  // W upward
      }
      // base plate film, and the coil a centimetre higher
      flux[0] = HEAT.hPlate * this.plateArea * (1 - this.plateWaxFrac) *
                (this.Tplate - this.T[0]);
      const qCoil = HEAT.hPlate * HEAT.coilArea * (this.Tplate - this.T[this.coilCell]);
      // metal cap fin
      flux[n] = HEAT.Ucap * this.capArea * (this.T[n - 1] - Tamb);

      for (let i = 0; i < n; i++) {
        const wall = HEAT.Uwall * this.perim[i] * this.dy * (this.T[i] - Tamb);
        let net = flux[i] - flux[i + 1] - wall + this.qWax[i];
        if (i === this.coilCell) net += qCoil;
        this.T[i] += h * net / this.cap[i];
      }
      // plate: bulb in, film out, own loss to the room through the housing
      const plateLoss = 6 * this.plateArea * (this.Tplate - Tamb);
      this.Tplate += h * (bulbWatts - flux[0] - qCoil - plateLoss - this.plateWaxQ) / HEAT.plateC;
    }
    if (guard >= 4096) this.stiffSteps = (this.stiffSteps || 0) + 1;
    this.qWax.fill(0);
    this.plateWaxQ = 0;
  }

  // Diagnostics -------------------------------------------------------------
  get Tbase() { return this.T[0]; }
  get Ttop() { return this.T[this.n - 1]; }
  get Tbulk() {
    let s = 0, c = 0;
    for (let i = 0; i < this.n; i++) { s += this.T[i] * this.cap[i]; c += this.cap[i]; }
    return s / c;
  }
  // Watts currently leaving through the glass and cap.
  lossWatts(Tamb) {
    let q = HEAT.Ucap * this.capArea * (this.Ttop - Tamb);
    for (let i = 0; i < this.n; i++) q += HEAT.Uwall * this.perim[i] * this.dy * (this.T[i] - Tamb);
    return q;
  }
}
