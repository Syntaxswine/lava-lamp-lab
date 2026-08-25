# Backlog — lava-lamp-lab

Reconciled against [docs/HANDOFF-2026-08-24.md](docs/HANDOFF-2026-08-24.md).
Ordered. Each item says what would count as done, because "improve the blobs" is
not a task.

Items 1 and 2 are open **defects with measurements attached**, not improvements.
Do them before anything that produces a new number, because every number in the
handoff was measured under both of them.

## 1. DONE — calibrate cohesion at the shipping spacing

**Done when** `node tools/calibrate-sigma.mjs` reports `σ_eff` within 10% of
`IFACE.sigma` **at the particle spacing the lamp ships with**, over at least
three resolvable gravities, with the exponent still near −0.5.

Completed in the hostile pass. The pair law now carries the missing reference
length, and a wide 107 mL rig runs at the exact 3.22 mm shipping spacing. With
`cohK = 10.738`, three resolved points read 2.143/2.011/2.154 mN/m against 2.100;
mean 2.103 (+0.1%), exponent −0.500, spread 7%, zero clamps.

`cohK = 16.1` was fitted when σ was 3.0 mN/m and the calibration puddle ran at a
2.15 mm spacing. Before the correction, the lamp ran at 3.22 mm while the label
said 4.5 mN/m. Re-measured there:

| spacing | σ target | σ_eff |
|---|---|---|
| 2.15 mm | 3.0 | 3.0 |
| 2.15 mm | 4.5 | 3.0–3.8 |
| **3.22 mm (as shipped)** | **4.5** | **1.9** |

The lamp's wax is held together by roughly 2.4× less tension than the panel says.
The direction matches the known scaling flaw: per-pair cohesion acceleration goes
as `dx^-2` where the Laplace-pressure argument asks for `dx^-1`.

**Two traps in the fix.** It is not a one-line retune — raising `cohK` 2.4×
raises the stiffness, tightens the capillary CFL by √2.4 (41 → ~64 Hz), takes the
frame budget from 39% to ~60%, and changes blob detachment. And at the shipping
spacing the calibration itself is barely possible: a 28 mL puddle is under five
particles deep, so only one of three gravities is resolvable and the tool
correctly refuses to fit.

So the real fix is probably not a new constant but a **corrected scaling** —
carry the extra power of `dx` in `cohScale` so the coefficient transfers between
resolutions, then re-fit once. A wider test container (the calibration is
currently confined by the globe profile) would let the puddle be deep enough to
measure at the shipping spacing.

## 2. PARTLY SUPERSEDED — get the thermal default off the edge of its window

**Done when** a warm-started lamp is in transit for ≥ 30% of samples over a
20-minute run, and `lamp-probe` resolves a cycle period in three runs out of
three.

Warm start now seeds four asynchronous developed bodies and a five-minute run is
in transit for 9/20 samples with zero clamps. A single global centre-of-mass
period is no longer the right acceptance test for four bodies, but the underlying
pool still sits close to its local rise threshold. Replace the old period
criterion with a 20-minute per-blob transit/dwell metric before closing this.

The shipping configuration sits on the **lower edge of the window rule**: the
pool equilibrates near 47.6 °C against a 47.9 °C rise threshold at the base. So
whether it lifts is decided by a residual heating rate of hundredths of a kelvin
per minute, the dwell dominates, and the cycle is erratic — three runs gave
128 s, *not resolved*, and *not resolved* (8 transits in 88 samples over 22 min).

This is a physically real state (a lamp running slightly too cool is exactly
this) and the window readout names it correctly. It is a poor **default**.

The margin is a handful of knobs and they interact, so measure rather than guess:
bulb power, `fCouple`, `Uwall`, and the wax excess density all move the plate
temperature relative to the threshold. Sweep the bulb over 25/30/35 W first,
since it is the one a user would reach for, and report transit fraction against
plate-minus-threshold margin. Note that item 1 will move this too — weaker
effective tension means a smaller blob can detach, which lowers the superheat
needed to lift.

## 3. DONE — a determinism baseline

**Done when** `tools/` has a check that runs a fixed seed to a fixed lamp time
and compares a fingerprint against a committed golden file, and it fails when any
physics constant is perturbed in the sixth decimal.

`Lamp` now owns a seeded generator used by reset, developed warm start and shake.
`tools/regression-probe.mjs` runs the same state twice, checks the committed
bitwise fingerprint, and has been mutation-tested: a `cohK` change of 0.000001
fails it.

## 4. DONE — make the blob population bigger than one

**Done when** a warm-started lamp holds ≥ 3 blobs of ≥ 8% of the wax each for a
full cycle, with the resolution readout still saying RESOLVED and the frame
budget still met.

Warm start now hands over a 55% connected pool/stem/bulb plus 18/15/12% round
parcels. The feed body reaches roughly 2.7× its equivalent diameter and pinches
into two macroscopic daughters after about 2.4 s, taking the visible population
from four to five. The 120 s surfactant-film drainage time keeps those daughters
separate; the focused probe records no velocity clamps and only about ten
microscopic strays.

At 60 mL a fully coalesced mass is a 24 mm sphere in a 33 mm bore, and the
capillary CFL makes brute force cost `n^1.5`. Do **not** simply raise the particle
count — measure first whether the barrier is coalescence or pinch-off resolution:
hold everything fixed and sweep `drainTime` (25 s) over 5/15/25/60, then
`disjoinK` (1.5) over 0.5/1.5/4. If blob count moves with the film parameters the
barrier is coalescence; if it does not move at all, the honest answer is to say so
rather than tune. Item 1 changes this too.

## 5. Mutation-test the instruments

**Done when** each of the five tools has been broken on purpose in the way most
likely to matter, and the ones that stayed green have been fixed or recorded.

| tool | mutation | what a pass would prove |
|---|---|---|
| `calibrate-sigma` | return a constant from `thickness()` | the exponent fit is doing work |
| `lamp-probe` | feed `cyclePeriod` pure noise | the prominence gate bites |
| `column-probe` | set `Uwall = 0` | the loss term is really in the balance |
| `solver-probe` | disable the λ clamp | it would have caught the original detonation |
| resolution check | force `peakDrho` constant | the "not measurable" branch is reachable |

Every "the instrument refuses" claim in the handoff is a claim about code that
has never been tried against a lie. This is the habit that found real defects in
the sibling project and has not been applied here at all.

## 6. Close the energy audit properly

**Done when** the probe reports `in − out − d(stored)/dt` and it is under 2% of
the input across a 10-minute run.

It currently prints in and out only, so the 16.3 / 14.7 W gap is *presumed* to be
the column and wax still warming and nobody has checked. `Wax` already computes
`enthalpy`; the column needs the matching sum. This is what would have caught the
substep-weighting bug at the moment it was introduced instead of via a plate at
−861 °C two hours later.

## 7. Advect heat with the return flow

**Done when** the aqueous energy equation carries a `v_aq ∂T/∂y` term, the
steady-state gradient is re-measured, and `nuC` is refitted.

The return flow moves volume but not heat, which matters most exactly when the
wax is moving fastest. Expect `nuC` to need refitting; it already had to move
from 1.2 to 1.8 when the globe was rescaled, which is what tells you it is a fit.

## 8. A preset for each failure mode

**Done when** three buttons put the lamp into *too cold to lift*, *too hot to
fall* and *seizing at the top* within a few seconds of lamp time.

The window rule is the best thing in this simulator and reaching its edges
currently means dragging a slider and waiting an hour. The states are already
detected and named by `Lamp.windowState`; what is missing is a way to *get*
there. Reuse `warmStart`, which builds an arbitrary steady state directly.

## 9. Smaller, once the above are settled

- `Volume.build` clears three full arrays per frame; only touched voxels need it.
  That clear is most of the ~2.4 ms draw.
- `blobStats` allocates a `Map` per call and runs on every HUD update.
- Mobile: the layout collapses correctly but 1800 particles at 47 Hz will not
  hold on a phone. Detect and drop the count, or say so.
