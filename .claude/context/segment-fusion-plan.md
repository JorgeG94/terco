# Weighted-segment fusion in the bin Fock launch: project plan

A follow-up to S15, and a candidate for the S16 list in
`.claude/context/symmetry-plan.md`. Written so that each stage can be handed
to an agent on its own: "implement stage W1 of
`.claude/context/segment-fusion-plan.md`". Read "Rules for every stage" first,
then your stage only. The working instructions in
`.claude/context/symmetry-agent-briefs.md` apply here too.

## Why this exists

S15 measured the skeleton Fock build on a V100 (cc70, nvfortran 26.5,
`-stdpar=gpu -acc -gpu=cc70,mem:separate`, `TERCO_LMAX=3`, cc-pVDZ,
`--econv 1e-10`, baseline `--symmetry off`). The quartet savings are real and
match `docs/SYMMETRY_FOCK.md` section 2.2's predictions for scheme 2a, but
below roughly twenty atoms they do not convert into wall time:

| system | atoms | bf | Fock/iter on | Fock/iter off | ratio | quartets on/off |
|---|---|---|---|---|---|---|
| SF6 | 7 | 102 | 0.0197 | 0.0070 | 0.36x | 2.44x |
| benzene | 12 | 114 | 0.164 (total) | 0.076 | 0.46x | 3.55x |
| benzene, cc-pVTZ | 12 | 264 | 2.377 | 1.147 | 0.48x | - |
| naphthalene | 18 | 180 | 0.0239 | 0.0234 | 0.98x | 3.81x |
| naphthalene, cc-pVTZ | 18 | 412 | 4.187 | 4.132 | 0.99x | - |
| pentacene | 36 | 378 | 0.0667 | 0.1429 | 2.14x | 3.87x |

Energies are unaffected: |dE| between the symmetric and the C1 run on the
same symmetrized geometry was 2.9e-13 to 1.4e-12 Eh against a 1e-10 gate,
with identical iteration counts on both sides. The two basis sets per
molecule show the ratio tracks the ATOM count, not the basis size:
naphthalene at 412 bf is still a wash while pentacene at 378 bf wins.

Where the loss is, for benzene (88 ms over 9 iterations), from nsys:

- **Not** in setup. `setup` is 0.018 s with symmetry against 0.015 s
  without, and its `pairs`, `schwarz` and `bins` sub-columns are the same to
  a millisecond. Detection, orbit tables and binning are not the cost.
- About a fifth is extra GPU kernel time. The ERI kernels (`trc_pc_*`) take
  90.36 ms over 1067 launches with symmetry against 69.65 ms over 495
  without: 2.16x the launches, but each only 1.66x cheaper (84.7 us against
  140.7 us), so the total RISES 1.30x. The device is left underfilled.
- About half is host-side per-segment cost: `cuStreamSynchronize` +26.2 ms
  (+1691 calls), `cuMemsetD32Async` +10.5 ms (+4576 calls),
  `cuMemcpyHtoDAsync` +8.2 ms (+2366 calls). `cudaLaunchKernel` barely
  moves, so raw launch latency is not the problem: it is the zeroing and
  synchronisation AROUND each segment.
- The remainder is host CPU work outside the CUDA API.

The same decomposition at 18 atoms comes out the other way: naphthalene's
ERI time is 177.0 ms with symmetry against 264.8 ms without, because the
per-launch time falls 3.7x while the launch count rises 2.44x. The fixed
per-segment costs are amortised once there is enough work per segment.

This is the risk `docs/SYMMETRY_FOCK.md` section 4.2(i) recorded and
deferred: "small segments cost warp occupancy at the tail. Measure on GPU
in S15." It has now been measured. The goal of this plan is to keep the
quartet savings and recover the occupancy, rather than to give up the
savings with a size rule.

## What is actually being changed

`binkernel_fock_bins` builds **one launch per (class, mask, weight)**
(`src/eri/trc_binkernel.F90`, the `do while (c0 <= nseg)` loop). The orbit
weight enters only as a scalar on the two operator factors:

```fortran
wt = plan%sW(c0)
jf = launch%jfac*wt
kf = launch%kfac*wt
do while (c1 < nseg)
   if (plan%class_key(c1 + 1, CLASS_RADIX) /= key) exit
   if (plan%sM(c1 + 1) /= mask) exit
   if (allocated(plan%sW)) then
      if (plan%sW(c1 + 1) /= wt) exit     ! this is what splits the launch
   end if
   ...
```

The realised weights are multiples of **one half**, not of one: `skeleton_walk`
and `bin_tiles_plan` both emit `half*real(orb%sec_w(r), dp)`, with
`half = 0.5` on a diagonal bin pair, and the sort key is `nint(2*sW)`
precisely because of it (`trc_bins.F90:1377`). W0 measured the realised set
as **{1/2, 1, 2, 4}**: weight 8 never occurs, because an orbit of 8 needs a
shell pair with trivial stabiliser and in all four test systems every shell
pair is fixed by at least one operation. So the ceiling is still 4 distinct
weights per class, for a different reason than "divides |G|".

Launch inflation per **Fock build** is 4.05x for benzene, not the 2.16x that
a whole-process count suggests: 306 of benzene's 1071 `trc_pc_*` launches
belong to the SAD atomic guess and are identical with symmetry on and off.
Per build it is 85 against 21 (naphthalene 4.19x, pentacene 4.29x, SF6
5.43x). Fusion does not recover all of that -- see W0's results.

**There is already a precedent in the tree for the fix.** The MD J engine
does not split its launches by weight: `weight_copies` in `src/eri/trc_mdj.F90`
multiplies the copies' primitive coefficient `mp_c` by the orbit size, then
pushes only that slice with `!$acc update device`, so one launch covers
every weight. The bin path is the outlier. Bringing it into line is the
substance of W2.

## Rules for every stage

- Read `AGENTS.md`, `.claude/skills/pr-review/SKILL.md`, this file's section
  for your stage, and `docs/SYMMETRY_FOCK.md` sections 4.2 and 6 before
  writing code.
- **`check_sym_fock` is the oracle and must not be weakened.** It pins
  `eri%nwork` (the quartet count), NOT `nlaunch`. Fusing changes how the
  work is launched and not which quartets run, so every pinned number in
  that test must be unchanged. If a change moves `nwork`, it is not this
  project.
- Kinds from `pic_types`. FORD `!!` comments only, never Doxygen. No bare
  `print`/`write(*,...)` in `src/`: use `pic_logger`'s `global_logger`.
- `src/generated/` is generator output. Change `scripts/`, rerun, commit
  both, and state the compile-time cost. CI regenerates and diffs.
- Every `!$acc enter data` has an `exit data` on every path out. Under
  `mem:separate` an unmapped array a kernel reads is a run-time illegal
  address that the host build will not show.
- A `do concurrent` body calls only intrinsics and same-module `pure`
  `!$acc routine seq` item routines.
- Host build first, then the OpenMP build at 4 threads (races only show
  there), then the device build.
- Performance claims state the system, basis, thresholds, GPU and baseline,
  and show unchanged energies. Per-class claims come from nsys, not wall
  time.
- Branch `sym/<stage>-<slug>` off the integration branch; merge commits, no
  squash; never merge without the maintainer.

## Stage map

| Stage | Title | Depends on | Status |
|---|---|---|---|
| W0 | Bound the win; no product code | - | done (GO for W2) |
| W1 | Host-side per-segment cost | W0 | todo |
| W2 | Per-segment weight in the kernels | W0, W1 | todo |
| W3 | Validation across builds | W1 or W2 | todo |
| W4 | Re-measure the S15 matrix | W3 | todo |

W1 and W2 are independent fixes to the two halves of the loss. W1 is cheap
and carries no register risk; W2 is the one that needs a generator change
and a full regeneration. Do W1 first and let W0's numbers decide whether W2
is worth its cost.

## W0 - Bound the win  [DONE: GO for W2]

**Measured** (cc-pVDZ, `--econv 1e-10`, `--symmetry auto`, one steady-state
Fock build, `build-gpu-l3`, V100 cc70). Launches are exact, not estimated:
the plan is sorted on the composite `(class, mask, weight)` key immediately
before the launch loop (`bins_plan_sort_by_class(CLASS_RADIX, by_mask=.true.)`,
`trc_binkernel.F90:305`), so a launch count IS a count of distinct triples
and a fused count IS a count of distinct `(class, mask)` pairs.

| system | atoms | segments | launches now | fused | ratio | off |
|---|---|---|---|---|---|---|
| benzene | 12 | 5712 | 85 | 33 | **2.58x** | 21 |
| naphthalene_d2h | 18 | 7871 | 88 | 36 | **2.44x** | 21 |
| pentacene_d2h | 36 | 8959 | 90 | 36 | **2.50x** | 21 |
| tops_sf6 | 7 | 4217 | 114 | 36 | **3.17x** | 21 |

Distinct weights per `(class, mask)` group: no group holds more than 4, and
the 4x ceiling is reached by only 3 of 33 groups (benzene) and 6 of 36
(SF6). SF6 is the worst case because every one of its groups splits 3 or 4
ways -- its small, heavily stabilised orbits put several weights into every
class.

Fusion lands at 33-36 launches and not at the symmetry-off 21, because
`bins_skel_cut_first` swaps bra and ket for some bin pairs and forms class
keys the full plan never does. That extra ~1.6x is the swap, not the weight,
and **W2 cannot remove it**.

**`launch_cap()` does not threaten W2: it is 0 by default and nothing in the
tree ever sets it.** `TRC_LAUNCH_CAP` appears only in `launch_cap()` itself
and two comments; no test, CMake file, script or CI job sets it, and the
loop guard is `if (lcap > 0 .and. ...)`. So zero fused groups exceed the cap
and it gives back none of the win. If it is ever enabled, the floor must sit
above the largest fused group (327k / 1.26M / 9.25M / 205k quartets for the
four systems) and **its setting and W2's benefit have to be decided
together**. It should not be scaled DOWN: the measured disease is launches
that are already too small. On 80 SMs at 128-thread blocks (~164k work items
to saturate), 68% of benzene's and 100% of SF6's quartets currently sit in
launches too small to fill the card; fusion takes both to ~46%.

**Throughput model** (each class's largest observed launch taken as the
throughput its fused launch would reach -- optimistic, and it charges nothing
for the per-segment weight load or any register cost): ERI time per build
falls benzene 10.00 -> 6.26 ms, naphthalene 16.55 -> 12.26, SF6 9.89 -> 7.29,
pentacene 70.21 -> 65.19. Benzene's symmetry-off ERI time is 7.74 ms per
build, so the predicted 6.26 ms is **below the C1 baseline** -- W2 alone
would turn benzene's ERI excess into a surplus, before W1 touches the host
side. Treat the launch counts as hard and the times as an upper bound.

Register pressure is untouched by W0 and remains W2's real gate.

**Original goal, for reference**: know, before touching product code, how
many launches fusion would remove and what that is worth.

- `TRC_SEG_STATS=1` already dumps the segment shape of each Fock build
  (`seg_stats_on`/`log_segment_shape` in `src/eri/trc_binkernel.F90`). Use
  it; do not add instrumentation.
- For benzene, naphthalene, pentacene and SF6 (cc-pVDZ, symmetry on):
  count segments and launches per (class, mask), and the number of DISTINCT
  weights in each group. The fused launch count is the number of
  (class, mask) groups; the current count is the number of
  (class, mask, weight) groups, further split by `launch_cap()`.
- Report current launches, fused launches, and the ratio, per system.
- **Check the `launch_cap()` interaction.** Fusing more segments into one
  launch reaches `lcap` sooner and re-splits. Report how many of the fused
  groups would exceed the current cap, and whether the cap wants rescaling.
  If the cap gives most of the win straight back, that is the finding and
  W2 should not start.
- **Done when** there is a table of current against fused launch counts for
  the four systems, with the cap's effect quantified, and a go/no-go for W2.

## W1 - Host-side per-segment cost

**Goal**: remove the synchronisation and zeroing that scale with the
segment count, without touching a kernel or a generator.

- Targets, measured on benzene (per-iteration figures are for the whole
  9-iteration run): +1691 `cuStreamSynchronize` calls, +4576
  `cuMemsetD32Async`, +2366 `cuMemcpyHtoDAsync`.
- Batch the per-segment buffer zeroing into one memset per launch group;
  hoist the host-to-device copies out of the segment loop; drop
  synchronisation between consecutive launches of the same class where
  nothing reads the intermediate result.
- No kernel signature changes, no regeneration, no new device arrays, so no
  register-pressure risk.
- **Done when** benzene's `cuStreamSynchronize`, `cuMemsetD32Async` and
  `cuMemcpyHtoDAsync` counts with symmetry on are within ~10% of the
  symmetry-off counts, `check_sym_fock` is unchanged and green, and the
  energies are bit-identical to before the change.

## W2 - Per-segment weight in the kernels

**Goal**: one launch per (class, mask). The weight becomes a per-segment
datum the kernel reads, instead of a scalar folded into `jfac`/`kfac`.

- Drop the `plan%sW(c1 + 1) /= wt` exit in `binkernel_fock_bins`.
  `bins_sort_by_weight` already keeps each weight contiguous and the
  representatives first, so no re-sorting is needed.
- Pass `plan%sW` into `pc_dispatch` and `rys_dispatch` alongside the
  per-segment arrays the kernels already index (`sA`, `sNB`, `sOA`, `sOB`,
  `sD`, `sOff`). Add it to the `enter data`/`exit data` lists.
- In `scripts/gen_perclass.py` and `scripts/gen_rys.py`, emit the new dummy
  and multiply the J and K contributions by the segment's weight at
  digestion. Hoist the weight into a scalar once per segment so the inner
  loops are unchanged.
- **Keep `jfac` and `kfac` as scalars.** The loop's own comment says the
  weight "cannot change a kernel's branch: they test jfac /= 0, kfac /= 0",
  and the generated kernels screen on `if (jfac /= 0.0_dp)`. A weight must
  never reach those tests. An orbit size is at least 1, so folding it in
  cannot zero them today, but the separation is the invariant that keeps
  J-only and K-only launches correct.
- 97 generated kernels change signature. Regenerate, commit generator and
  output together, and state the measured compile-time change.
- **Register pressure is the risk that can end this stage.** At
  `TERCO_LMAX=3` 648 of 1243 generated kernels already sit at the 255
  register ceiling, with spills currently at zero everywhere
  (`cuobjdump -res-usage` on the LINKED `libterco.so`, never the `.o`).
  Measure before and after; if the new argument introduces spills in the
  reachable canonical classes, stop and report rather than proceeding.
- **Done when** fused launch counts match W0's prediction, spills are still
  zero, and W3 passes.

## W3 - Validation

**Goal**: the numbers are unchanged on every build.

- `check_sym_fock` unchanged: symmetrized skeleton G equals full G for every
  group, every class up to `TERCO_LMAX`, Cartesian and spherical, J-only,
  K-only and J+K, exchange tiles on and off, incremental Fock on and off.
  The pinned `eri%nwork` counts must not move.
- The central SCF oracle for every group's test molecules, RHF, UHF and one
  DFT functional, skeleton path on.
- `symmetry=off` keeps the old launch plan exactly, with the goldens
  unchanged (`compare_goldens.py`).
- OpenMP build at 4 threads.
- **Device run comparing energies before and after the change.**
  `.claude/context/gpu-pitfalls.md` says a layout change to something
  kernels index "has triggered an nvfortran miscompile before" and must only
  be made "with a GPU run that compares energies before and after". W2 is
  exactly such a change. This check is not optional.
- `fortitude check` clean; `check_staleness.py` clean.

## W4 - Re-measure

**Goal**: say whether the loss is gone, in the same terms S15 used.

- Benzene, naphthalene, pentacene and SF6, cc-pVDZ, symmetry on against
  off: Fock time per iteration with iteration counts, quartets launched,
  launch counts, and the `cuStreamSynchronize`/`cuMemsetD32Async`/
  `cuMemcpyHtoDAsync` counts.
- nsys per-kernel breakdown for at least naphthalene and one other system.
- Targets: benzene from 0.46x to at least 1.0x; pentacene's 2.14x preserved
  or better. The ideal for benzene's ERI time is about 19.6 ms (69.65 ms
  divided by the measured 3.55x quartet reduction) against 90.36 ms today,
  so there is room for more than break-even.
- Report against the ideal |G| and against the |G|/2 shell-level ceiling
  that `docs/SYMMETRY_FOCK.md` section 6 open question 1 describes.

## If W2 cannot land

W1 plus a size rule still recovers most of the loss, and the size rule is
worth having regardless because it costs nothing. Key it on the atom count,
not the basis size: the measurements above show a molecule holds its own
ratio across a 2.3x basis increase (benzene 0.46 to 0.48, naphthalene 0.98
to 0.99) while 18 atoms break even and 36 atoms win. terco already picks the
J engine by a measured size rule (`check_sym_fock_mdrule` tests exactly that
pattern), so the mechanism exists.

## Aside: parity skipping appears inert on the CLI path

Found while measuring W0, and **not a W-stage item** -- it belongs to
whoever owns the parity work. On this GPU build, a CLI run engages **no**
z-masked variant kernel at all: `nsys` shows zero `trc_pc_*z_*` launches for
benzene, naphthalene, pentacene and SF6 under `--symmetry auto`, and
benzene shows zero with cc-pVDZ **and** with 6-31G*.

The same molecule and the same basis engage parity in the test:
`check_sym_parity.F90:424` drives `run_molecule("benzene", "D2h",
"basis_sets/6-31gs.json", .false., 1, 0.95_dp, ...)` -- expected dead axis
**1**, masked share at least **0.95** -- plus
`check(error, ctx%eri%nlaunch_variant > 0 .and. ctx%eri%nwork_variant > 0, ...)`
at line 344. That test passes. So the test asserts parity engages and the
CLI does not do it, for the same molecule and basis.

A plausible mechanism is `parity_choose_axis`
(`src/eri/trc_parity.F90:50-88`), whose gate is an **exact** floating-point
comparison -- `if (m(s) == 1 .and. sh_r(d, s) /= 0.0_dp) exact = .false.` --
and whose own comment says a plane carrying rounding noise "is no dead plane
for the kernels". A symmetrised benzene reports a 4.4e-16 Bohr displacement.
But that is a hypothesis, not a measurement: `nwork_masked`,
`nwork_variant`, `nlaunch_variant` and `parity_axis` are not surfaced through
the CLI or `--json`, so which step on the CLI path differs from the test's
has not been determined. Note also that zero *variant* launches does not by
itself prove `parity_axis == 0`: a masked class with no variant kernel falls
back to its base kernel.

**Suggested first step**: surface those four counters through `--json` (or a
`performance` log line) and run benzene through the CLI. That answers it in
one run, and it is the check both merge reviews independently asked for --
parity skipping has never been tested through the S17 default-auto path,
only through the Fortran API with `sym=` passed explicitly.

## Open questions for the maintainer

1. Is a size rule acceptable as a stopgap before W1/W2 land, given that
   S17 has already made symmetry the default and small molecules are
   currently slower than they were?
2. ~~`real(dp)` array or `integer` orbit size?~~ **Answered by W0: it must
   be `real(dp)`.** The realised weight set is {1/2, 1, 2, 4} -- the 1/2
   comes from a diagonal bin pair, not from an orbit size, so an integer
   orbit size cannot represent it. The alternative is the integer plus the
   existing `sD` diagonal flag with the kernel forming the half itself,
   which halves device memory at the cost of a branch in the per-segment
   prologue. Prefer `real(dp)` and spend the review effort on the register
   measurement.
3. SF6 reaches only 2.44x although it is nonplanar D2h, where
   `docs/SYMMETRY_FOCK.md` section 6 expects planarity to be the limiter.
   Its atoms sit on many symmetry elements at once (S is fixed by all eight
   operations, each F by four), so the orbits are small regardless. Should
   section 6 record local site symmetry as a second limiter beside
   planarity?
