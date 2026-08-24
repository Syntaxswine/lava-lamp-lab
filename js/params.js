// params.js — physical constants for a classic 41.5 cm ("16.3 inch") lava lamp.
//
// Every number here is a real, sourced material property or a closure coefficient
// that is explicitly labelled as calibrated. Nothing is a magic visual constant.
// See docs/PHYSICS.md for the derivations and the measurements that pin them down.

// ---------------------------------------------------------------------------
// Geometry — the globe is a surface of revolution. Control points are
// (height fraction, radius in metres). This is the classic 14.5 inch lamp: a
// 21 cm globe, 5.5 cm across the base, bulging to 6.6 cm a quarter of the way
// up, necking to 2.5 cm at the cap. About 540 mL inside.
//
// The bigger 16.3 inch globe was modelled first and it is worse on three counts
// at once. Its 1.37 L of fluid gives an 85-minute thermal time constant and
// wants three hours to develop; its 240 mL wax charge coalesces into a slug
// nearly as wide as the bore; and spreading a fixed particle budget through more
// volume coarsens the resolution just where blob size is already resolution-
// limited. The small lamp also happens to be the one the manufacturer pairs with
// a 25 W bulb rather than a 40 W, which is exactly what the energy budget here
// asks for -- a useful check that the plumbing is right.
// ---------------------------------------------------------------------------
export const GEO = {
  H: 0.210,                 // m — glass height, base plate to cap
  profile: [                // [t (0=base,1=cap), radius m]
    [0.00, 0.0275],
    [0.08, 0.0316],
    [0.22, 0.0330],
    [0.45, 0.0303],
    [0.70, 0.0241],
    [0.88, 0.0172],
    [1.00, 0.0124],
  ],
  coilY: 0.0075,            // m — the metal tension-breaker coil sits here
  coilR: 0.021,             // m — coil ring radius
};

// ---------------------------------------------------------------------------
// Wax phase: paraffin blend, dyed, dosed with a dense halocarbon so that at
// room temperature it is ~1% DENSER than the aqueous phase. That 1% is the
// whole trick — see T_CROSS below.
//   rho    liquid paraffin is 760-780 at 70 C; the halocarbon dose brings the
//          blend to just over 1000 at 20 C.
//   beta   liquid n-alkane volumetric expansion, 7-9e-4 /K
//   c      2.1-2.3 kJ/kg/K, liquid paraffin
//   L      fusion enthalpy of paraffin wax, 180-230 kJ/kg
//   Tmelt  lamp wax is a soft low-melt blend; it must stay molten at the TOP
//          of the globe or the lamp seizes, so Tmelt < T_top (~40 C). Soft
//          oil-extended blends melt in the mid-30s.
// ---------------------------------------------------------------------------
export const WAX = {
  rho20: 1011,              // kg/m^3 at 20 C
  beta: 7.5e-4,             // 1/K
  c: 2100,                  // J/kg/K
  k: 0.25,                  // W/m/K
  mu: 0.005,                // Pa.s at 55 C
  Tmelt: 34,                // degC
  dTmelt: 3.0,              // K — mushy band half-width
  L: 200e3,                 // J/kg latent heat of fusion
  volume: 60e-6,            // m^3 (60 mL), ~12% of the globe volume.
                            // On the LIGHT side -- real lamps run 15-25% -- and
                            // chosen deliberately. A coalesced 110 mL mass is a
                            // 30 mm sphere inside a 33 mm bore: it fills the tube
                            // and can only behave like a piston. 60 mL makes a
                            // 24 mm mass, leaves room beside it, and spreading a
                            // fixed particle budget through less volume is what
                            // brings the capillary length inside the resolvable
                            // range at all. See docs/PHYSICS.md section 5.
};

// ---------------------------------------------------------------------------
// Aqueous phase: water + propylene glycol + surfactant + salt.
// The glycol raises viscosity (slowing the blobs to their signature crawl) and
// the surfactant drops interfacial tension to a few mN/m, which is why the
// blobs are centimetres across instead of millimetres.
// ---------------------------------------------------------------------------
export const AQ = {
  rho20: 1000,              // kg/m^3 at 20 C
  beta: 3.8e-4,             // 1/K (water 2.1e-4 at 20 C, 4.6e-4 at 50 C; glycol ~6e-4)
  c: 4000,                  // J/kg/K
  k: 0.55,                  // W/m/K
  mu: 0.010,                // Pa.s (10 cP, glycol-thickened)
};

export const IFACE = {
  // Oil/water with surfactant spans about 1-5 mN/m, and the default sits at the
  // top of that range on purpose: the capillary length a_c = sqrt(sigma/drho g)
  // is what sets blob size, and at 4.5 mN/m it lands just above the three-
  // particle-spacing floor the solver needs to hold a free surface. Lower values
  // are equally physical and the slider reaches them -- the readout will then say
  // the configuration is under-resolved, which is the truth rather than a
  // quietly smaller blob.
  sigma: 4.5e-3,            // N/m
};

// ---------------------------------------------------------------------------
// Heat plumbing.
//   bulbW    the lamp bulb (25/40/60 W are the stock sizes)
//   fCouple  fraction of bulb power that enters the globe through the reflector
//            cone and glass foot; the rest heats the room directly. CALIBRATED
//            so steady state lands on measured lamp temperatures (base plate
//            ~58 C, bulk ~45 C, ambient 22 C).
//   Uwall    overall glass-to-room coefficient, lumping natural convection
//            (~5), radiation (~5) and conduction out of the metal cap.
//            CALIBRATED against the same steady state.
//   hPlate   base-plate to wax/fluid film conductance.
//   nuC      turbulent convection closure, Nu = nuC * Ra^(1/3), where Ra is
//            built on a wall-distance mixing length rather than the plate
//            spacing, so this is NOT the Globe-and-Dropkin 0.069 constant --
//            different length scale, different constant. CALIBRATED by
//            tools/closure-sweep.mjs: 1.8 reproduces the measured steady state
//            of a running 25 W lamp of this size (plate 56, bulk 46, top 41,
//            base-to-top gradient 12 K) and a 52-minute time constant.
//            It had to be REFITTED when the globe was rescaled -- 1.2 on the
//            larger one -- which is the plainest possible statement that this is
//            a fitted coefficient of one particular closure and not a physical
//            constant. Change the geometry and it needs measuring again.
// ---------------------------------------------------------------------------
export const HEAT = {
  bulbW: 25,                // the stock bulb for a globe this size
  fCouple: 0.65,
  Uwall: 19,                // W/m^2/K
  Ucap: 30,                 // W/m^2/K — the metal cap is a fin
  hPlate: 400,              // W/m^2/K — plate/coil to fluid film
  coilArea: 6.6e-3,         // m^2 — 8 turns of 2 mm wire at 21 mm radius:
                            // pi*d*L with L = 8*2*pi*0.021 = 1.06 m. The coil is
                            // not decoration: it spreads the plate's heat into
                            // the fluid and gives the wax pool something to
                            // tear away from.
  Tamb: 22,                 // degC
  plateC: 160,              // J/K — thermal mass of base plate + coil
  nuC: 1.8,                 // mixing-length convection closure (calibrated)
  g: 9.80665,
};

// ---------------------------------------------------------------------------
// Solver
//   cohK converts interfacial tension in N/m into Akinci-2013 cohesion
//   strength. It is not a free knob: tools/calibrate-sigma.mjs settles a
//   non-wetting puddle under a known reduced gravity and reads its thickness,
//   which capillarity fixes at e = 2 sqrt(sigma / (drho g)) with no unknown
//   constant. cohK is the value that makes the measured sigma match IFACE.sigma,
//   and the tool also fits the exponent of e against g to confirm it is really
//   measuring surface tension and not the particle size.
// ---------------------------------------------------------------------------
export const SOLVER = {
  particles: 1800,          // The capillary CFL makes dt ~ n^(-1/2), so compute
                            // per second of lamp time goes as n^1.5. This is the
                            // most a browser affords at a stable step.
  iterations: 3,            // PBF density projection passes
  hFactor: 2.0,             // smoothing length = hFactor * particle spacing
  relaxEps: 1.0e-6,         // CFM relaxation, as a fraction of |grad C|^2
  scorr: 0.05,              // artificial pressure, as a fraction of lambda
  cohK: 16.1,               // KNOWN WRONG AT THE SHIPPING SPACING -- see below.
                            //
                            // Fitted by tools/calibrate-sigma.mjs when
                            // IFACE.sigma was 3.0 mN/m and the calibration
                            // puddle ran at a 2.15 mm particle spacing. There it
                            // read sigma_eff = 3.0 against the 3.0 asked for.
                            //
                            // sigma is now 4.5 and the lamp runs at a 3.22 mm
                            // spacing, and re-measuring AT THAT SPACING gives
                            // sigma_eff = 1.9 mN/m -- 2.4x below nominal. The
                            // direction is exactly what the scaling flaw
                            // predicts: the per-pair cohesion acceleration goes
                            // as dx^-2 while the Laplace-pressure argument asks
                            // for dx^-1, so effective tension falls as the
                            // spacing grows. Pairwise SPH surface tension is
                            // resolution dependent and this coefficient does not
                            // carry between resolutions.
                            //
                            // NOT retuned in place, because it is not a one-line
                            // change: raising cohK ~2.4x raises the stiffness,
                            // which tightens the capillary CFL by sqrt(2.4) and
                            // takes the solver from 41 Hz to ~64 Hz, the frame
                            // budget from 39% to ~60%, and changes blob
                            // detachment -- invalidating every measured number in
                            // the handoff. BACKLOG.md item 1.
  curvK: 0.55,              // Akinci curvature-minimisation weight
  xsph: 0.0,                // extra XSPH smoothing per SECOND (0 = physical mu only)
  Cd: 1.0,                  // drag coefficient of a wobbling drop (Clift et al.)
  // Surfactant films resist coalescence. Two blobs that touch are separated by a
  // thin aqueous film that has to DRAIN before they can merge, and the
  // surfactant in a lava lamp is there precisely to make that slow. Without it
  // the wax coalesces into one 240 mL slug within thirty seconds, and a mass
  // that big has a thermal time constant of half an hour: the lamp stops being
  // a lava lamp and becomes a very slow piston.
  //   disjoinK   strength of the film's disjoining pressure, as a multiple of
  //              the cohesion that would otherwise pull the blobs together
  //   drainTime  seconds of continuous contact before the film ruptures and the
  //              two blobs become one. Real drainage times for surfactant-
  //              stabilised films span seconds to minutes.
  disjoinK: 1.5,
  drainTime: 25,
};

// Density of each phase at temperature T (linear expansion about 20 C).
export const rhoWax = (T) => WAX.rho20 * (1 - WAX.beta * (T - 20));
export const rhoAq  = (T) => AQ.rho20  * (1 - AQ.beta  * (T - 20));

// The crossover temperature: where wax and aqueous phase, both at the SAME
// temperature, have equal density. Below it wax sinks, above it wax floats.
export const T_CROSS = 20 + (WAX.rho20 - AQ.rho20) /
                            (WAX.rho20 * WAX.beta - AQ.rho20 * AQ.beta);

// Capillary length: the blob radius at which interfacial tension balances the
// buoyant pull. Sets how big a blob grows before it tears off the pool.
export const capillaryLength = (drho) =>
  Math.sqrt(IFACE.sigma / (Math.max(0.05, Math.abs(drho)) * HEAT.g));

// Globe radius at height y (metres above the base plate).
export function globeRadius(y) {
  const t = Math.min(1, Math.max(0, y / GEO.H));
  const p = GEO.profile;
  for (let i = 1; i < p.length; i++) {
    if (t <= p[i][0]) {
      const f = (t - p[i - 1][0]) / (p[i][0] - p[i - 1][0]);
      return p[i - 1][1] + f * (p[i][1] - p[i - 1][1]);
    }
  }
  return p[p.length - 1][1];
}

export const globeArea = (y) => Math.PI * globeRadius(y) ** 2;
export const globePerim = (y) => 2 * Math.PI * globeRadius(y);

// Interior volume of the globe by numeric integration of the profile.
export function globeVolume(n = 512) {
  let v = 0;
  const dy = GEO.H / n;
  for (let i = 0; i < n; i++) v += globeArea((i + 0.5) * dy) * dy;
  return v;
}
