// lamp.js — the whole lamp: bulb, base plate, aqueous column, wax.
//
// Owns the clock. Two clocks, actually, and the difference matters:
//   lampTime  seconds of simulated lamp operation
//   wallTime  seconds you have been watching
// timeScale is how many lamp-seconds pass per watched second. It is honest
// acceleration -- the same equations, integrated with more substeps -- not a
// separate fast thermal clock.
//
// fastForward() is the exception that earns its keep: while no wax anywhere has
// positive buoyancy, nothing can move, so the mechanics have nothing to lose and
// only the heat equation is doing work. It advances the thermal solution alone
// and STOPS THE MOMENT the first parcel of wax becomes buoyant. The gate is
// measured, not assumed: wax.buoyancyHeadroom() is checked every step.

import { Column } from './column.js';
import { Wax } from './wax.js';
import { HEAT, WAX, AQ, SOLVER, IFACE, T_CROSS, rhoWax, rhoAq, GEO } from './params.js';

const DT = 1 / 120;               // mechanical substep, seconds
const COM_SAMPLE = 0.5;           // seconds between centre-of-mass samples

function seededRandom(seed) {
  let s = seed >>> 0;
  return () => {
    s += 0x6D2B79F5;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export class Lamp {
  constructor(opts = {}) {
    this.seed = opts.seed ?? 0x1a2b3c4d;
    this.rand = seededRandom(this.seed);
    this.column = new Column(opts.cells ?? 96);
    this.wax = new Wax(opts.particles ?? SOLVER.particles, opts.waxVolume ?? WAX.volume);
    this.env = {
      bulbW: HEAT.bulbW,
      Tamb: HEAT.Tamb,
      sigma: IFACE.sigma,
      muAq: AQ.mu,
      iterations: SOLVER.iterations,
      gravity: 1,
      cohesion: 1,
      dRho: WAX.rho20 - AQ.rho20,     // exposed as a control: the whole trick
    };
    this.reset();
  }

  reset() {
    this.rand = seededRandom(this.seed);
    this.column.reset(this.env.Tamb);
    this.wax.reset(this.env.Tamb, this.rand);
    this.lampTime = 0;
    this.firstMelt = null;
    this.firstRise = null;
    this.comTrace = [];
    this._comAcc = 0; this._comN = 0; this._comAt = 0;
  }

  setParticles(n) {
    this.wax = new Wax(n, this.wax.volume);
    this.wax.reset(this.env.Tamb, this.rand);
  }

  setWaxVolume(v) {
    this.wax.setVolume(v);
    this.wax.reset(this.env.Tamb, this.rand);
  }

  get bulbIn() { return this.env.bulbW * HEAT.fCouple; }

  // -------------------------------------------------------------------------
  // Warm start. This 486 mL globe has a roughly 48-minute thermal time constant
  // plus the wax's latent heat, so an honest cold start takes about two hours to
  // lift. That wait is physical, and the cold-start path reproduces it.
  //
  // This builds a developed state directly instead: relax the column to the
  // steady state its own energy budget implies, then seed a pool drawing one
  // upper bulb plus three smaller bodies. It is an INITIAL CONDITION, not a simulated
  // history or an animation path. Everything after hand-over is fully coupled.
  // -------------------------------------------------------------------------
  warmStart() {
    this.reset();
    // dt = 10 keeps the wax conduction inside its own stability limit without
    // substepping, and 1000 steps is 2.8 hours -- about 3.5 time constants,
    // which lands within a per cent of the fixed point.
    const dt = 10;
    for (let i = 0; i < 1000; i++) {
      this.wax.heat(dt, this.column, this.env);
      this.column.step(dt, this.bulbIn, this.env.Tamb);
    }
    const w = this.wax;
    w.seedDeveloped(this.column, (Tf) => this.riseThreshold(Tf), this.rand);
    this.lampTime = 1200 * dt;
    this.startedWarm = true;
    this.firstMelt = 0;
    this.firstRise = null;
    this.comTrace = [];
    this._comAcc = 0; this._comN = 0; this._comAt = this.lampTime;
    return this.diagnostics();
  }

  // One display frame worth of simulation. dtWall is real seconds elapsed;
  // timeScale multiplies it. Returns the number of substeps actually taken.
  advance(dtWall, timeScale = 1, budgetSteps = 8) {
    const want = Math.min(dtWall * timeScale, 0.5);
    let steps = Math.min(budgetSteps, Math.max(1, Math.round(want / DT)));
    const dt = want / steps;
    for (let s = 0; s < steps; s++) this.substep(dt);
    return steps;
  }

  substep(dt) {
    this.wax.step(dt, this.column, this.env);
    this.column.step(dt, this.bulbIn, this.env.Tamb);
    this.lampTime += dt;
    this.track();
  }

  // Thermal-only fast-forward, valid only while nothing can rise. Returns the
  // lamp seconds actually skipped and why it stopped.
  fastForward(maxSeconds = 6 * 3600, dt = 5) {
    let skipped = 0;
    let reason = 'limit';
    while (skipped < maxSeconds) {
      const span = this.wax.buoyantSpan(this.column, this.env.dRho, this.env.sigma);
      if (span.lifts) { reason = 'liftoff'; this.span = span; break; }
      this.wax.heat(dt, this.column, this.env);
      this.column.step(dt, this.bulbIn, this.env.Tamb);
      this.lampTime += dt;
      skipped += dt;
      this.track();
    }
    return { skipped, reason };
  }

  // Event log, plus a record of the wax centre of mass. The cycle period comes
  // off a spectrum of that record rather than off turning points: counting
  // up-swings of a noisy signal with a minimum gap measures the gap, and this
  // one dutifully reported 13 s for a lamp whose blobs were plainly taking far
  // longer to make the round trip.
  track() {
    if (this.firstMelt === null && this.wax.moltenFraction > 0.02) this.firstMelt = this.lampTime;
    if (this.firstRise === null && this.wax.maxRise > 0.004) this.firstRise = this.lampTime;
    this._comAcc += this.wax.centreOfMass;
    this._comN++;
    if (this.lampTime - this._comAt >= COM_SAMPLE) {
      this.comTrace.push(this._comAcc / this._comN);
      this._comAcc = 0; this._comN = 0;
      this._comAt = this.lampTime;
      if (this.comTrace.length > 1024) this.comTrace.shift();
    }
  }

  // Dominant period of the centre-of-mass record, or null if the record is too
  // short or the peak is not prominent enough to call a period.
  get cyclePeriod() {
    const sig = this.comTrace;
    if (sig.length < 64) return null;
    const N = sig.length;
    const mean = sig.reduce((a, b) => a + b, 0) / N;
    // detrend: a lamp still warming up has a rising baseline that would swamp
    // any oscillation and put the spectral peak at the lowest bin
    let sx = 0, sxy = 0, sxx = 0;
    for (let i = 0; i < N; i++) { sx += i; sxy += i * (sig[i] - mean); sxx += i * i; }
    const slope = (sxy - sx * 0 / N) / (sxx - sx * sx / N);
    const w = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      const hann = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / (N - 1));
      w[i] = (sig[i] - mean - slope * (i - sx / N)) * hann;
    }
    let best = { f: 0, p: -1 };
    const spec = [];
    // periods from 6 s to a quarter of the record
    const pMin = 6, pMax = Math.max(12, N * COM_SAMPLE / 4);
    for (let b = 0; b < 240; b++) {
      const period = pMin + (pMax - pMin) * b / 239;
      const om = 2 * Math.PI * COM_SAMPLE / period;
      let re = 0, im = 0;
      for (let i = 0; i < N; i++) { re += w[i] * Math.cos(om * i); im += w[i] * Math.sin(om * i); }
      const pw = re * re + im * im;
      spec.push(pw);
      if (pw > best.p) best = { f: period, p: pw };
    }
    const sorted = spec.slice().sort((a, b) => a - b);
    const median = sorted[spec.length >> 1];
    if (best.p < 6 * Math.max(median, 1e-30)) return null;
    return best.f;
  }

  // -------------------------------------------------------------------------
  // The window rule. A lava lamp only runs when the crossover temperature --
  // where wax and fluid at the same temperature weigh the same -- falls between
  // the top of the globe and the base plate, AND the wax stays molten up top.
  //   too cold: the plate never reaches the crossover, nothing lifts
  //   too hot : the top never falls below it, everything parks at the cap
  // -------------------------------------------------------------------------
  get crossover() {
    const dr = this.env.dRho;
    const rw20 = AQ.rho20 + dr;
    return 20 + dr / (rw20 * WAX.beta - AQ.rho20 * AQ.beta);
  }

  // Three conditions, and all three have to hold at once:
  //   the plate can drive wax at the base past its local rise threshold
  //   the top is below the crossover, so a relaxed blob becomes heavy again
  //   the top is above the melting point, so the wax up there stays liquid
  // Note the first test is against riseThreshold(Tbase), NOT against the
  // crossover. The crossover is where the two phases match AT THE SAME
  // TEMPERATURE; at the base the fluid is hot and therefore light, so the wax
  // has to beat a lower bar than the crossover to lift. Testing the plate
  // against the crossover reports a stalled lamp that is in fact running.
  get windowState() {
    const c = this.column;
    if (c.Tplate < this.riseThreshold(c.Tbase)) return 'too cold to lift';
    if (c.Ttop > this.crossover) return 'too hot to fall';
    if (c.Ttop < WAX.Tmelt) return 'seizing at the top';
    return 'cycling';
  }

  // Temperature a wax parcel must reach to rise through fluid at Tf, and the
  // temperature it must fall to before it sinks again. These two numbers are
  // the lamp: the gap between them is how hard it has to work.
  riseThreshold(Tf) {
    const rw20 = AQ.rho20 + this.env.dRho;
    return 20 + (rw20 - rhoAq(Tf)) / (rw20 * WAX.beta);
  }

  // The step the physics allows: the smaller of the capillary-wave CFL and a
  // ceiling that keeps a paused-then-resumed lamp from taking one huge stride.
  physDt() {
    return Math.min(1 / 20, Math.max(1 / 240, this.wax.capillaryDt(this.env.sigma)));
  }

  // Can this particle count even represent the blobs the physics wants?
  //
  // The capillary length a_c = sqrt(sigma / (drho g)) is the radius at which
  // buoyancy tears a blob off the pool. SPH needs about three particle spacings
  // across a radius before it can hold a free surface at that scale, so when
  // a_c falls below 3 dx the lamp cannot make the blobs its own parameters call
  // for: it makes larger ones, which need more superheat to lift, so they rise
  // faster and cycle sooner. That is a resolution limit, not a lamp, and the
  // readout says which one you are looking at.
  resolution() {
    // The density difference that tears a blob off the pool is the THERMAL one
    // -- a few kg/m3 of superheat -- not the ~11 kg/m3 the wax carries at room
    // temperature. Using the room-temperature excess reports a capillary length
    // 40% too small and an under-resolution that is not there. Take the measured
    // peak from the wax itself, floored so a lamp with nothing buoyant in it
    // reports the cold-start value rather than infinity.
    const span = this.wax.buoyantSpan(this.column, this.env.dRho, this.env.sigma);
    const drho = span.peakDrho || 0;
    const floor = 3 * this.wax.dx;
    // With nothing buoyant there is no density difference to build a capillary
    // length from, and reporting the floored value comes out as a 30 mm blob and
    // a confident "resolved" -- a verdict about a quantity that was not measured.
    if (drho < 0.25) return { known: false, floor, drho, resolved: false, ratio: 0, capillary: 0 };
    const ac = Math.sqrt(this.env.sigma / (drho * HEAT.g));
    return { known: true, capillary: ac, floor, drho, resolved: ac >= floor, ratio: ac / floor };
  }

  diagnostics() {
    const c = this.column, w = this.wax;
    // A blob is a lump you can see. Counting every stray particle as a blob
    // makes the readout say 180 when the lamp plainly holds six, and hides the
    // one time the droplet count IS the story: right after a shake.
    const all = w.blobStats(4096);
    const stats = all.filter((b) => b.count >= 8);
    const droplets = all.length - stats.length;
    const big = stats.slice(0, 6);
    const mean = big.length ? big.reduce((a, b) => a + b.radius, 0) / big.length : 0;
    return {
      lampTime: this.lampTime,
      Tplate: c.Tplate, Tbase: c.Tbase, Tbulk: c.Tbulk, Ttop: c.Ttop,
      gradient: c.Tbase - c.Ttop,
      lossW: c.lossWatts(this.env.Tamb),
      inW: this.bulbIn,
      crossover: this.crossover,
      riseAtBase: this.riseThreshold(c.Tbase),
      sinkAtTop: this.riseThreshold(c.Ttop),
      window: this.windowState,
      molten: w.moltenFraction,
      waxT: w.meanT, waxTmax: w.maxT,
      blobs: stats.length,
      droplets,
      biggest: w.biggestBlob,
      meanBlobRadius: mean,
      maxStretch: big.reduce((a, b) => Math.max(a, b.verticalStretch || 0), 0),
      pinches: w.pinches,
      coalescences: w.coalescences,
      maxRise: w.maxRise,
      riseSpeed: w.typicalRise(),
      resolution: this.resolution(),
      physDt: this.physDt(),
      com: w.centreOfMass / GEO.H,
      cyclePeriod: this.cyclePeriod,
      firstMelt: this.firstMelt,
      firstRise: this.firstRise,
      startedWarm: !!this.startedWarm,
      clamps: w.clamps,
      viscClamped: !!w.viscClamped,
    };
  }
}

export { T_CROSS, rhoWax, rhoAq };
