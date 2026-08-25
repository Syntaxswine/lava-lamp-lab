# Lava Lamp Lab — 3D

A lava lamp simulated from the physics that actually makes one work: a wax phase
dosed ~1% denser than the fluid around it, expanding twice as fast with
temperature, crossing over into buoyancy at **49.1 °C**, and doing it inside a
temperature window a few kelvin wide.

Nothing here is animated. There is no blob path, no keyframe, no noise field
pretending to be convection. There is a 25 W bulb, a base plate, a column of
water-and-glycol solving an energy equation, and 1800 particles of paraffin with
temperatures.

Zero dependencies. One WebGL2 pass, plain ES modules, no build step.

```bash
node serve.mjs 8811
```

## Why it exists

It began as a request to turn a WebGPU Navier–Stokes particle demo into a lava
lamp. That demo's lab-instrument shell — orbit viewport, HUD, a panel of physical
controls — is a good frame for a lamp, so it was kept and the physics underneath
was replaced entirely.

## What you can do with it

| control | what it changes | what to watch |
|---|---|---|
| **Bulb** 10–50 W | power into the globe | at 10 W the plate never reaches the rise threshold and the lamp refuses to start; at 50 W the cap goes past the crossover and every blob parks at the top |
| **Room** 6–34 °C | ambient | the same two failure modes, reached from the other side |
| **Wax excess density** 3–22 kg/m³ | how much heavier the wax is at 20 °C | this *is* the crossover temperature — watch it move in the readout |
| **Interfacial tension** 0.6–9 mN/m | the surfactant dose | blob size through the capillary length, *and* the solver step, which is set by the capillary CFL |
| **Fluid viscosity** 2–34 mPa·s | the glycol | how fast blobs crawl |
| **Wax charge**, **Particles** | resolution and load | the readout says when the particle count cannot resolve the blobs the parameters call for |
| **View** | Lamp / Temperature / Rise–sink | the third one colours every parcel by whether it is currently heavier or lighter than the fluid beside it |

Buttons: **Warm start** seeds a developed four-body circulation directly: a hot
upper bulb drawing a neck from the pool, with three smaller parcels aloft. It is
an initial condition, not a prescribed animation; the connected feed body then
elongates and pinches under the normal solver. **Cold start**
puts a solid plug of wax on a cold plate and makes you wait, like the real thing.
**Fast-forward** skips the dead time — but only while nothing in the tank could
possibly rise, and it stops the instant that stops being true. **Shake it** does
what shaking a lava lamp does: emulsifies the wax into droplets too small to
lift, and takes a long time to clear.

## The five findings

**1. A lava lamp runs inside a window, not above a threshold.** Three conditions
have to hold at once: the plate must exceed the *local* rise threshold, the cap
must be below the crossover, and the cap must still be above the wax's melting
point. Miss any one and the lamp stops. This is why a draught kills it and why it
never quite works in a hot room.

**2. The bulb size is not a styling choice.** A 25 W bulb in this 486 mL globe
gives a steady state of plate 55 °C, bulk 45 °C, cap 40 °C — crossover neatly
between them. Put 40 W in and the cap climbs past 49 °C and the lamp stalls with
everything at the top. 25 W is what the manufacturer ships with a globe this
size, and the energy budget here asks for it independently.

**3. The hour of waiting is the wax's latent heat as much as the water's.**
1.9 kJ/K of thermal mass leaking 0.66 W/K gives a 48-minute time constant, and
then 12 kJ has to go into melting 60 g of paraffin at constant temperature. Cold
start reaches liftoff at **119 minutes** of lamp time, and the number was not
tuned to — it falls out of the energy budget.

**4. Surface tension, not particle count, is what caps a real-time particle lava
lamp.** Explicit interfacial tension carries a CFL condition
`dt ≤ 0.25 √(ρh³/2πσ)`. With `h ~ (V/n)^(1/3)` that makes `dt ~ n^(-1/2)`, so the
cost of a second of lamp time goes as **n^1.5** — you pay for the particles and
again for the shorter step they force. Doubling the count nearly triples the
bill. The solver here computes its step from that condition every frame, so
moving the tension slider from 3.1 to 8 mN/m visibly takes it from 59 Hz to
about 94 Hz.

**5. The wax “skin” is an interface and a trapped film, not a rubber shell.**
The oil/water tension falls slightly with temperature, while the thin aqueous
film between touching blobs drains much faster when warm. The hot plate also has
a wetting boundary energy. Together those effects make the bottom pool spread
across the foot and make base contacts more likely to coalesce, while cooler
blobs aloft retain their surfactant films for much longer.

## Verifying it

```bash
node tools/column-probe.mjs 25 22
node tools/closure-sweep.mjs
node tools/calibrate-sigma.mjs
node tools/solver-probe.mjs
node tools/shape-probe.mjs
node tools/skin-probe.mjs
node tools/lamp-probe.mjs 9 1800 25 22 warm
node tools/regression-probe.mjs
```

Three of these refuse to answer rather than answer badly.

`calibrate-sigma` will not fit a capillary exponent using puddles thinner than it
can resolve, and will not report from a single gravity — `e ~ (Δρg)^(-1/2)` is the
signature of capillarity, and checking the exponent is what makes the number a
measurement of surface tension rather than of the particle size.

`lamp-probe` will not name a cycle period whose spectral peak is not prominent.
An earlier version counted up-swings of the centre of mass with a minimum gap
between them and reported the gap.

The resolution check will not return a verdict when nothing in the tank is
buoyant, because with no Δρ there is no capillary length to compare against, and
the floored value came out as a confident *resolved* about a quantity that was
never measured.

Everything calibrated is labelled as calibrated in `js/params.js`, next to what
it was calibrated against. The cohesion law now carries the missing resolution
length, and the calibration runs in a wide 107 mL vessel at the lamp's shipping
3.22 mm spacing. Across three resolved gravities it reads **3.044, 2.798 and
3.036 mN/m** against the 3.100 control (mean 2.959), with a thickness exponent of
−0.505 and no clamp events. The non-wetting calibration explicitly disables the
separate hot-plate adhesion boundary. The rig and the earlier
false 1.9 mN/m result are written up in PHYSICS.md §4.6.

## Known limits

The default configuration **is** resolved — roughly 10.6 mm capillary radius
against a 9.7 mm floor at peak buoyancy — but only because the globe, the wax charge and the interfacial
tension were all chosen to put it there, and the readout flips to UNDER-RESOLVED
the moment you move them. The first configuration tried, a 16.3-inch globe with
240 mL of wax, was under-resolved by 2.6× and would have needed ~44,000
particles.

Warm start is now a deliberately constructed developed-state initial condition:
a connected pool/stem/upper bulb plus three smaller parcels. Their paths are not
prescribed; the coupled solver owns every step after hand-over. In the focused
probe the connected body grows from **2.13× to 2.34×** its equivalent diameter
and breaks into macroscopic daughters. The hot footprint sustains 88–89% of the
plate radius before returning blobs make it reach the flared wall. At 48 °C the
unresolved film drains in about 58 s versus 248 s at 40 °C; a two-minute run
records five isolated coalescences, ends with five visible bodies, peaks at
**1.80 cm/s**, remains resolved, and records **zero** velocity clamps. Cold start
remains the honest solid-plug experiment.

A four-blob asynchronous lamp need not have one globally coherent period, so the
spectral cycle instrument can still refuse even while blobs are plainly moving.
The fixed seed and committed state fingerprint make that refusal reproducible.
The largest remaining physical simplification is the 1-D aqueous temperature
field: it carries return volume but not a three-dimensional warm wake.

Full model, closures, calibrations and the rest:
[docs/PHYSICS.md](docs/PHYSICS.md).
