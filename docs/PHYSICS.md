# The physics of a lava lamp, and how this one is built

## 1. The trick

A lava lamp is a heat engine with one moving part, and the part is a density
inversion.

The wax is dosed with a dense additive so that at room temperature it is about
**1% heavier** than the water-and-glycol around it. That is why a cold lamp is a
lump of wax on the bottom. But the wax expands with temperature about **twice as
fast** as the aqueous phase does — roughly 7.5 × 10⁻⁴ per kelvin against
3.8 × 10⁻⁴ — so heating closes the gap. Somewhere above room temperature the two
phases weigh the same, and above that the wax is the lighter one.

Set the two linear expansions equal and solve:

```
T_cross = 20 + (ρ_w20 − ρ_a20) / (ρ_w20 β_w − ρ_a20 β_a)
```

With ρ_w20 = 1011, ρ_a20 = 1000, β_w = 7.5e-4, β_a = 3.8e-4 kg/m³ and 1/K, that
is **49.1 °C**. Everything the lamp does is organised around that one number.

## 2. The window rule

The crossover alone does not make a lamp run. Buoyancy compares the wax against
the fluid *at its own height*, and the fluid is not at one temperature — a
running lamp holds a base plate near 56 °C, a convecting bulk near 45 °C, and a
cap near 39 °C. So the threshold a parcel of wax must beat depends on where it
is:

```
T_rise(T_fluid) = 20 + (ρ_w20 − ρ_a(T_fluid)) / (ρ_w20 β_w)
```

At the base, with fluid at 48 °C, the wax must reach **48.7 °C** to lift. At the
cap, with fluid at 39 °C, it must fall below **44.0 °C** to sink again. Three
conditions have to hold at once:

| condition | meaning | fails when |
|---|---|---|
| `T_plate > T_rise(T_base)` | the plate can make wax buoyant | the room is cold, or the bulb is weak |
| `T_top < T_cross` | a relaxed blob becomes heavy again | the room is hot, or the bulb too strong |
| `T_top > T_melt` | the wax stays liquid up top | the lamp is barely warm |

That is the whole behaviour of a lava lamp, including both of its famous failure
modes. **Too cold and nothing lifts. Too hot and everything parks at the cap.**
The lamp is not robust; it works inside a window a few kelvin wide, which is why
a draught kills it and why running one in a hot room ruins the effect.

Note the first test is against `T_rise(T_base)`, *not* against `T_cross`. The
crossover is where the phases match at the *same* temperature; at the base the
fluid is hot and therefore already light, so the wax has a lower bar to clear.
Testing the plate against the crossover reports a stalled lamp that is, in fact,
running — this simulator did exactly that until the readout was fixed.

## 3. Why it takes so long to start

The globe modelled here is the classic 14.5-inch lamp: a 21 cm glass holding
about **486 mL**, of which **60 mL** is wax. Its heat capacity is

```
C ≈ 0.43 kg × 4000 + 0.06 kg × 2100 + 160 (plate)  ≈  1.9 kJ/K
```

and it leaks to the room through `U·A ≈ 0.66 W/K`. A **25 W** bulb couples about
16.3 W into the globe, so the steady state sits ~25 K above ambient and the time
constant is `C / UA ≈ 48 minutes`. On top of that the wax has to *melt*:
0.06 kg × 200 kJ/kg = 12 kJ, another ~12 minutes of the bulb's whole output spent
at constant temperature.

The measured steady state is plate **55.2 °C**, base 52.0, bulk 45.4, cap 40.4 —
and the time to 63% of it is **52 minutes**, which is what the manufacturer says
about a lamp this size.

Note that 25 W is the bulb this globe actually ships with, and it is also the
power the energy budget here asks for: 40 W in a globe this small drives the cap
past the crossover and the lamp stops cycling. That pairing was not put in by
hand. It fell out.

The 16.3-inch globe was modelled first, at 1.37 L and 40 W, and behaved the same
way three times slower: an 85-minute time constant and about three hours to full
development. Also correct, and also what the makers of the larger lamps tell you
— but see §5 for why the smaller lamp is the one that got kept.

## 4. What is actually simulated

### 4.1 The aqueous phase: a 1-D column

The interesting fluid is the wax. The water/glycol matters only as the density
the wax is weighed against and the heat bath it exchanges with, so it is a stack
of 96 well-mixed horizontal slices carrying an energy equation: face conduction,
wall loss `U(T − T_amb)` per unit height, a metal cap fin, a base-plate film, and
two-way coupling with the wax.

Vertical transport is closed with the Rayleigh–Bénard scaling `Nu = C · Ra^(1/3)`
applied on a mixing length `l = min(y, H − y, 2R(y))` and only where the
stratification is unstable. That is what produces the real structure of a running
lamp: a hot film on the plate, a nearly isothermal convecting bulk, and a cooler
cap. Molecular conduction alone would take six days to warm the top.

`C = 1.8` is **calibrated**, by `tools/closure-sweep.mjs`, against the measured
steady state of a running lamp. It is *not* the Globe–Dropkin constant of 0.069
— that correlation is built on the plate spacing, this one on a wall-distance
mixing length, and a different length scale means a different constant.
Pretending otherwise would be borrowing authority from a paper that does not
cover this.

And it had to be **refitted** when the globe was rescaled: 1.2 on the larger one,
1.8 on this one. That is the plainest possible statement of what it is. A closure
coefficient that moves with the geometry is a fit, not a constant, and anyone
changing the globe has to measure it again.

### 4.2 The wax: a position-based fluid

Position-based fluids (Macklin & Müller 2013) with interfacial tension
(Akinci et al. 2013), per-particle temperature, and latent heat by the
apparent-heat-capacity method. Blobs are connected components of the neighbour
graph, found by union-find every step.

Three couplings run to the aqueous phase, and each is the continuum result for a
drop of the blob's *actual* size, not something invented per particle:

- **buoyancy**, per particle, from its own temperature against the fluid at its
  own height. This is the engine.
- **drag**, per *blob*, quadratic, with `Cd ≈ 1` for a wobbling deformable drop.
  Per-particle drag would make a 3 cm blob rise at the same speed as a 2 mm
  droplet, which is the one thing a lava lamp visibly does not do.
- **heat**, over the blob's surface area, with the film coefficient from
  `Nu = 2 + 0.6 Re^0.5 Pr^(1/3)`, distributed across the blob's surface particles
  so the total flux equals the continuum flux exactly.

### 4.3 Return flow

The globe is sealed and both phases are incompressible, so the net volume
crossing any horizontal plane is zero: whatever the wax carries up, the aqueous
phase carries down beside it.

```
v_aq(y) = − Q_wax(y) / A_open(y)
```

This was missing at first, and the consequence was unmistakable: the wax
coalesced into a single tube-filling slug and sat there, rising through a fluid
that never had to get out of the way. With the constraint in place a blob
spanning most of the bore faces a return jet in the annulus beside it, its drag
climbs sharply, and it stops behaving like a piston in a frictionless bore.

`A_open` is floored at 12% of the bore. As the wax approaches a perfect seal the
true return velocity diverges, and below the floor the model has stopped
resolving the film between wax and glass.

### 4.4 Why the blobs stay apart

Surfactant films resist coalescence. Two blobs that touch are separated by a thin
aqueous film that has to *drain* before they can merge, and the surfactant in a
lava lamp is there precisely to make that slow.

Without it the wax coalesced into a single slug within thirty seconds of the
hand-over, and a mass that big has a thermal time constant of

```
τ = ρ c a / 3h ≈ 1000 × 2100 × 0.0386 / (3 × 14)  ≈  30 minutes
```

so the lamp stopped being a lava lamp and became a very slow piston.

Blobs therefore carry a **persistent identity** across steps — not an index into
this frame's component list, which would be a different blob next frame and would
leave the drainage clock timing nothing — and particles belonging to different
blobs push each other apart, using the cohesion spline's own magnitude with the
opposite sign. A pair in continuous contact for longer than `drainTime` is
allowed to merge; a pair that separates forgets its clock. A component that
splits leaves the larger half holding the id and names the smaller half afresh,
which is how a pinch-off produces a blob whose film has not yet drained against
its parent.

That film is the temperature-sensitive “skin” visible in a real lamp. Elevated-
temperature oil-in-water experiments report much more frequent coalescence as
the aqueous-film viscosity falls and thermal capillary fluctuations rupture the
film; alkane/water tension itself falls only slightly. Surfactant gradients also
generate Marangoni stresses that oppose drainage. See the primary experiments on
[heated emulsion coalescence](https://www.nature.com/articles/s41598-021-89919-5)
and [surfactant stresses in thin films](https://www.nature.com/articles/s42005-025-02410-9).

The resolved surface therefore uses a small linear thermal correction about its
calibration temperature,

```
σ(T) = σ45 [1 − 0.006 (T − 45 °C)]
```

while the unresolved film drainage clock receives the stronger measured-sign
closure

```
r_film(T) = clamp( exp[(T − 44 °C) / 5.5 K], 0.25, 6 ) .
```

At 48 °C a continuously touching film reaches the 120 s rupture exposure in
about 58 s; at 40 °C it takes about 248 s. Its disjoining barrier scales as
`σ(T)/√r_film`, so a warm contact is also easier to maintain before rupture.
When one film ruptures, other clocks involving either parent reset: the violent
coalescence flow and newly redistributed surfactant make them new films rather
than an array of pre-aged clocks that all fail in an avalanche.

The hot plate is a separate three-phase boundary. A normal adhesion acceleration

```
a_wet = 3 σ(T)/(ρh²) · q² · smoothstep(38 °C, 47 °C, T)
```

acts through the bottom three smoothing lengths. Incompressibility converts that
normal attraction into lateral spreading; there is no prescribed radial target.
This is the particle-scale counterpart of a bridge growing radially after the
intervening film ruptures in substrate-wetting experiments
([Nature Communications](https://www.nature.com/articles/s41467-021-26015-2)).

The connected-component bond uses 1.95 particle spacings, essentially the full
compact-support neighbourhood. The old 1.75-spacing cutoff declared a drawn
surface disconnected while its particles were still exchanging pressure and
capillary force, shedding numerical spray. A second guard handles the opposite
resolution failure: an hourglass neck thinner than 2.3 spacings can remain
graph-connected after continuum capillarity would have pinched it. The solver
bins the cross-section, requires two substantial lobes and an interior population
below 55% of both shoulders, and requires that condition to persist for two
seconds before cutting only bonds that cross the unresolved waist. Ordinary
geometric splits still happen without this correction; fresh persistent ids and
the surfactant-film pressure then separate the daughters.

### 4.5 The timestep is not a choice

Explicit surface tension carries its own CFL condition. The shortest capillary
wave the kernel can hold has wavelength ~h, and the integrator has to resolve its
period:

```
dt ≤ 0.25 √( ρ h³ / (2π σ) )
```

The simulator computes this every frame from the smoothing length and whatever
the tension slider currently says. It uses a 0.18 safety factor, plus the
inverse-square-root of the cohesion resolution-transfer factor, rather than
riding the theoretical 0.25 limit. The coldest fully molten interface is about
5% stiffer than the 45 °C slider reference, and the step includes that worst
mobile surface. Move the tension from 3.1 to 8 mN/m and the shipping solver goes
from 59 Hz to about 94 Hz by itself.

Running past it does not look like divergence. It looks like a velocity clamp
firing on two-thirds of the particles every step while the lamp merely seems
sluggish — which is exactly what a hardcoded 1/30 s did here, at 1.77× over the
limit: **660,000 clamp events in twelve seconds of lamp time**. With the step
taken from the condition instead, the same twelve seconds produce **eight**.

This is also the ceiling on a particle lava lamp in a browser, and it is worth
writing down as a scaling law. With `h = 2(V/n)^(1/3)`,

```
dt ~ n^(-1/2)     so    compute per second of lamp time ~ n^1.5
```

Doubling the particle count nearly triples the cost, because you pay for the
particles *and* for the shorter step they force. That is why the default is 1800
particles and not as many as memory allows.

### 4.6 Interfacial tension, and how it was calibrated

Akinci's cohesion strength is quoted as a bare coefficient. Turning it into a
real surface tension in N/m takes a measurement against a case with a known
answer, and `tools/calibrate-sigma.mjs` uses the maximal thickness of a
non-wetting puddle:

```
e = 2 √(σ / (Δρ g)) · sin(θ/2)
```

(de Gennes, Brochard-Wyart & Quéré, *Capillarity and Wetting Phenomena*, §2.3.)
The wax has no adhesion to the glass here, so θ = 180° and the puddle settles at
exactly twice the capillary length — no unknown constant in front of it. Settle a
puddle under a known reduced gravity, read the thickness, invert for σ.

The tool refuses two things. It will not fit the exponent using puddles it cannot
resolve — a puddle two particle layers deep cannot be thinner than the particles,
and including those points bends the fitted exponent toward zero and inflates σ.
And it will not report a single-gravity answer: `e ~ (Δρ g)^(-1/2)` is the
signature of capillarity specifically, so the exponent is fitted and checked.

The first rig got this wrong in three different ways. It fitted at a 2.15 mm
spacing while the lamp ships at 3.22 mm; its 28 mL puddle was confined by the
lamp's tapered 27.5 mm foot; and setting the film force to zero did not disable
the film's persistent topology, so split pieces still refused to bond. Its
hard-coded 10 mm start was also below every 15–24 mm predicted thickness even
though the comment claimed the puddle would settle downward. The resulting
1.9 mN/m “measurement” was partly a reading of the wall, finite drop and initial
condition.

The corrected pair scale carries the missing resolution length:

```
cohScale = 1 / (rho0² dx cohRefDx),     cohRefDx = 2.15 mm
```

so the pair acceleration transfers as `dx^-1`, matching the Laplace-pressure
argument, instead of the old `dx^-2`. The calibration now uses 3200 particles in
a 60 mm-radius cylindrical vessel: a 107 mL puddle at the exact shipping spacing.
The film identity state machine is bypassed, the puddle starts 35% above the
analytic height, and only points deeper than 4.5 particle layers enter the fit.

At `cohK = 10.738`, target σ = 3.10 mN/m. The sweep gravity scales with the
requested tension so all target puddles retain the same resolvable thickness:

| reduced gravity | thickness | σ_eff |
|---|---:|---:|
| 0.021 m/s² | 24.14 mm | 3.044 mN/m |
| 0.034 m/s² | 17.93 mm | 2.798 mN/m |
| 0.052 m/s² | 15.25 mm | 3.036 mN/m |

The mean is **2.959 mN/m** (−4.5%), the thickness exponent is **−0.505**, the
spread is 8%, and all three runs record zero clamps. The calibration is explicitly
non-wetting; enabling the hot-plate boundary in this rig flattens the puddle and
would measure substrate adhesion rather than cohesion.

Two lessons worth keeping. **Re-run every calibration after every model change**:
adding the anti-coalescence film shattered the calibration puddle into 71 pieces
and silently changed what the tool was measuring (it now runs with the film off,
since it calibrates cohesion, not the film). And **calibrate at the configuration
you ship**, not at whatever the test rig found convenient.

An earlier version of this tool used Rayleigh's drop oscillation,
`ω² = 8σ/(ρa³)`. It does not work here: with the wax viscosity and the aqueous
drag both in play, a centimetre drop at 3 mN/m is *overdamped* — it overshoots
once and creeps back without ringing. Reading a frequency off that trace gave
half-periods scattered over a factor of fifty and a "σ" that moved three orders
of magnitude depending on which crossings you kept.

## 5. Known limits

**Blob size is resolution-limited, and the whole configuration was chosen around
that.** The capillary length `a_c = √(σ/Δρ g)` is the radius at which buoyancy
tears a blob off the pool, and SPH needs about three particle spacings across a
radius before it can hold a free surface at that scale. The default sits inside
that bound: `a_c ≈ 10.6 mm` against a `3·dx = 9.7 mm` floor, from a **measured**
Δρ near 2.8 kg/m³ at peak lift, and the readout says RESOLVED.

Getting there took three changes, all stated rather than hidden. The globe was
scaled down from 16.3 to 14.5 inches; the wax charge is 60 mL, on the light side
of a real lamp's 15–25%; and the default interfacial tension is 3.1 mN/m inside
the 1–5 mN/m surfactant-rich band. Each is a physically legitimate point
in the parameter space, and each was picked because it is where the solver can
express its own capillary length.

Push any of them the other way and the readout flips to UNDER-RESOLVED with the
factor, because the consequence chain is real and worth naming: the lamp cannot
shed a blob at the capillary scale, so the pool superheats further before *any*
blob is big enough to lift, so it lifts hotter, rises faster and cycles sooner.
At the original 2600 particles in 240 mL that factor was 2.6×, and the fix would
have been about 44,000 particles.

The resolution verdict also refuses to answer when it cannot. With nothing
buoyant in the tank there is no Δρ to build a capillary length from, and reporting
a floored value came out as a confident *resolved* about a quantity that was
never measured. It now says so instead.

**Warm start is constructed.** It seeds a developed circulation with one
connected pool/stem/upper bulb and three smaller parcels aloft. Their positions
are an initial condition, not paths: the normal PBF, buoyancy, drag, heat
exchange, return flow and surfactant film own the run immediately after
hand-over. Cold start still begins from the single solid plug and invents no
developed history.

The focused shape probe measures the connected feed body at 2.13× its equivalent
spherical diameter initially and 2.34× immediately before it divides into
macroscopic daughters. The hot footprint sustains about 89% of the foot radius
before returning wax spreads to the flared wall. A two-minute probe records five
isolated film coalescences while continuing pinch-off leaves five visible bodies;
peak speed is 1.80 cm/s and the velocity clamp never fires.

**A cycle is not necessarily one number.** Several asynchronous bodies do not have
to put one prominent peak in the global centre-of-mass spectrum. The detector
still refuses when that peak is absent; motion is reported separately as transit
fraction and visible-blob speed. Runs are now seeded and checked against a
committed bitwise fingerprint, so a refusal is reproducible rather than luck.

**Other simplifications, stated:**

- The return flow moves volume but does not advect heat; the aqueous energy
  equation sees only diffusion, wall loss and wax exchange.
- Paraffin expands ~10% on melting. That volume change is ignored; only thermal
  expansion is modelled. It matters at cold start, not to a running lamp.
- The aqueous phase is 1-D. There is no horizontal structure to its temperature,
  so a blob does not leave a warm wake beside it.
- The velocity clamp (`0.35 h` per step) is a safety net, not physics. It
  **counts**, and the count is in the diagnostics, so a run that leaned on it
  says so instead of quietly looking fine.

## 6. Numbers used

| quantity | value | source |
|---|---|---|
| globe | 486 mL, 21 cm | classic 14.5-inch lamp |
| wax charge | 60 mL (12%) | light; chosen for resolution, see §5 |
| wax density at 20 °C | 1011 kg/m³ | paraffin + halocarbon dose, set ~1% above the aqueous phase |
| wax expansion β | 7.5e-4 /K | liquid n-alkanes, 7–9e-4 |
| wax specific heat | 2100 J/kg/K | liquid paraffin, 2.1–2.3 kJ/kg/K |
| wax latent heat | 200 kJ/kg | paraffin fusion, 180–230 kJ/kg |
| wax melting point | 34 °C | soft oil-extended blend; must stay below T_top |
| wax viscosity | 15 mPa·s | viscous surfactant-rich blend; preserves a drawn neck |
| aqueous density at 20 °C | 1000 kg/m³ | water + glycol + salt |
| aqueous expansion β | 3.8e-4 /K | water 2.1e-4 at 20 °C, 4.6e-4 at 50 °C; glycol ~6e-4 |
| aqueous viscosity | 10 mPa·s | glycol-thickened |
| interfacial tension at 45 °C | 3.1 mN/m | oil/water with surfactant, 1–5 mN/m, §4.4–§5 |
| thermal tension slope | −0.6%/K | small compared with the film-viscosity effect, §4.4 |
| bulb | 25 W | the stock bulb for a globe this size |
| coupled fraction | 0.65 | **calibrated** to the measured steady state |
| glass–room coefficient | 19 W/m²/K | **calibrated**; lumps convection, radiation, the metal cap |
| convection closure C | 1.8 | **calibrated**, `tools/closure-sweep.mjs`; refitted per geometry |
| cohesion coefficient | 10.738 | **calibrated**, `tools/calibrate-sigma.mjs` |
| particles | 1800 | the most the capillary CFL affords in real time, §4.5 |
| solver step | 1/59 s | **derived**, including the colder molten-surface stiffness, §4.5 |

## 7. The instruments

| tool | measures | refuses to |
|---|---|---|
| `column-probe.mjs` | warm-up curve and steady state of the aqueous column | — |
| `closure-sweep.mjs` | steady state vs the convection constant | — |
| `calibrate-sigma.mjs` | σ from puddle thickness, and the exponent | report from under-resolved puddles, or from one gravity |
| `solver-probe.mjs` | the density projection, step by step from a cold plug | — |
| `shape-probe.mjs` | peak elongation, neck scale and macroscopic split events | — |
| `skin-probe.mjs` | hot-footprint radius, σ(T), film barrier and drainage time | accept a narrow hot pool or temperature-independent film |
| `lamp-probe.mjs` | the whole lamp: window, blobs, rise speed, cycle | name a cycle period whose spectral peak is not prominent |
| `regression-probe.mjs` | fixed-seed state and configuration fingerprint | accept a changed physical state without baseline review |

The cycle period comes off a detrended spectrum of the wax centre of mass. An
earlier version counted up-swings with a minimum gap between them and dutifully
reported "13 s" for a lamp whose blobs were plainly taking minutes — it was
measuring the minimum gap. Likewise "rise speed" is now the volume-weighted mean
over rising blobs; the running maximum of the fastest particle reports the worst
transient a run ever had and calls it a speed.
