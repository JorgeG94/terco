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
- The rest is host-side CUDA API cost that scales with the LAUNCH count:
  `cuStreamSynchronize` +26.2 ms (+1691 calls), `cuMemsetD32Async` +10.5 ms
  (+4576 calls), `cuMemcpyHtoDAsync` +8.2 ms (+2366 calls).
  `cudaLaunchKernel` barely moves, so raw launch latency is not the problem.
  **Do not add the sync time to the kernel time**: a `do concurrent` launch
  is synchronous, so the host blocks inside `cuStreamSynchronize` while that
  launch's own kernel runs, and nsys reports the host-side API duration and
  the device-side kernel duration separately -- they overlap (15.5 us mean
  sync against an 84.7 us mean kernel). The unambiguously additive host cost
  is memset + HtoD = **18.7 ms, about 21% of the 88 ms**, not half. W1
  established this from the published figures; it has not been confirmed by
  a fresh profile.
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
| W1 | Host-side per-segment cost | W0 | done (no target) |
| W2 | Per-segment weight in the kernels | W0 | done (landed; spills still zero) |
| W3 | Validation across builds | W1 or W2 | todo |
| W4 | Re-measure the S15 matrix | W3 | todo |

W1 and W2 were meant to be independent fixes to the two halves of the loss.
W1 found that there are not two halves: both the extra GPU kernel time and
the extra host-side API cost scale with the launch count, so W2 is the only
stage that addresses either. Go straight to W2.

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

## W1 - Host-side per-segment cost  [DONE: no target; do not attempt]

**W1 as specified cannot be done, because its premise was wrong.** The three
API costs are not per segment. They are a fixed multiple of the LAUNCH
count -- nvfortran's prologue for one `do concurrent` offload region, not
terco code:

| API | delta (calls) | per extra launch |
|---|---|---|
| `cuMemsetD32Async` | +4576 | 7.94 |
| `cuMemcpyHtoDAsync` | +2366 | 4.11 |
| `cuStreamSynchronize` | +1691 | 2.94 |

Benzene is 85 launches per Fock build with symmetry against 21 without, so
64 extra per build and 576 over 9 iterations. Taking the constants as 8/4/3
and applying them to the whole-process counts predicts 4608/2304/1728
against the measured 4576/2366/1691 -- within 1.4%, 2.7% and 2.1%, on both
sides of the comparison.

The per-segment reading is excluded by the same data: benzene's plan has
**5712 segments** against 85 launches (~67 segments per launch), so a call
issued per segment would make the deltas 5-14x larger than they are.

Audited and confirmed: the `do while (c0 <= nseg)` loop
(`trc_binkernel.F90:385-475`) contains **zero** `!$acc` directives, zero
allocations and zero zeroing -- only host integer arithmetic and one
dispatch call. The dispatchers are a `select case` to one generated driver,
one `do concurrent` per dispatch, with every per-thread workspace a
fixed-size automatic inside the item routine, so nvfortran emits no implicit
per-launch copy. All 30 array arguments the Fock kernels read are explicitly
mapped, and every module array in device code carries `!$acc declare
create`, so there is no implicit copy to hoist. terco issues no memset at
all: every zeroing on the path is a `do concurrent` kernel, once per build.

terco's own per-build device traffic is ~11 HtoD, ~2 DtoH and ~6 syncs, i.e.
~3-4% of the measured deltas and almost all of it common to both sides.
There is nothing to batch, hoist or drop that would register.

**Consequence: the entire measured host-side excess is proportional to the
launch count, so W2 is the only stage that can remove it.** W2 takes
benzene 85 -> 33 launches, which is ~61% off each of the three counts:

| per Fock build | off (21) | now (85) | after W2 (33) |
|---|---|---|---|
| `cuMemsetD32Async` | 168 | 680 | 264 |
| `cuMemcpyHtoDAsync` | 84 | 340 | 132 |
| `cuStreamSynchronize` | 63 | 255 | 99 |

**And W1's old acceptance criterion was unreachable even in principle.** It
asked for the three counts within ~10% of the symmetry-off counts, which
reduces to "launches within 10% of 21 per build". W0 measured full fusion
landing at 33, not 21, because `bins_skel_cut_first` forms class keys the
full plan never does. Restate any such criterion against **8/4/3 times the
fused launch count** (33 per build for benzene), or drop it for W4's
wall-time target.

Rejected on the brief's own rules, not on effort: fusing mask-split launches
would reorder segments, changing the order of accumulations into `jmat` and
so the last bits of the energy on the serial host build, which
`check_sym_parity.F90:504` relies on reproducing bit-for-bit. It is also
worth nothing on benzene, where W0 measured fused-by-`(class, mask)` = 33 =
the number of distinct classes, so the mask splits nothing.

## W2 - Per-segment weight in the kernels  [DONE]

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

### W2 results

**Measured** (cc-pVDZ, `--econv 1e-10`, `--symmetry auto`, V100 cc70,
nvfortran 26.5, `TERCO_LMAX=3`; per steady-state Fock build, from nsys as
the difference between two consecutive `--iterations` runs -- W0's method):

| system | launches before | after | W0 predicted | ERI ms before | after | sym off |
|---|---|---|---|---|---|---|
| benzene | 85 | **33** | 33 | 8.50 | **4.50** | 6.34 |
| naphthalene_d2h | 88 | **36** | 36 | 15.70 | **11.79** | 24.46 |
| pentacene_d2h | 90 | **36** | 36 | 71.48 | **66.61** | 184.92 |
| tops_sf6 | 218 | **72** | 36 | 14.89 | **7.29** | 3.23 |

Every system lands on exactly the number of distinct `(class, mask)` groups:
the per-kernel launch histograms after fusion are {1: 33}, {1: 36}, {1: 36}
and {1: 36, 2: 18}. Fusion is therefore complete; what is left over the
symmetry-off 21 is `bins_skel_cut_first`'s swap, as W0 said.

**SF6 does not match W0's 36, and its baseline is 218 and not 114, because
parity skipping is NOT inert on this branch.** Benzene and SF6 each engage
18 z-masked kernels under `--symmetry auto` with cc-pVDZ, against the zero
W0's aside reported. The masks split the `(class, mask)` groups further, so
SF6 has 72 groups rather than 36 -- a 3.03x reduction against W0's predicted
3.17x. **The "parity skipping appears inert on the CLI path" aside above is
stale for this branch** and whoever owns the parity work should re-check it.

Benzene's ERI time is now **below** its symmetry-off baseline (4.50 against
6.34 ms, and 3.86 against 4.44 on the next build), which is what W0's
throughput model predicted. naphthalene and pentacene come out 2.07x and
2.78x faster than symmetry-off on ERI time. SF6's ERI is still 2.3x above
its C1 baseline, down from 4.6x: at 7 atoms with heavily stabilised orbits
there is not enough work per class, and W2 does not fix that.

**Registers: no spills introduced.** `cuobjdump -res-usage` on the LINKED
`libterco.so`, `TERCO_LMAX=3`, kernel names normalised for the `_<line>_gpu`
suffix that moves with the source. Over 1225 device functions: REG at the
255 ceiling 587 -> 586, mean 183.6 -> 183.6, total STACK 2 684 272 ->
2 684 608 bytes (+336, +0.013%), max STACK unchanged at 169 600, LOCAL 0
throughout, and **every kernel whose STACK was 0 still has STACK 0**. For the
353 `trc_pc_*` Fock kernels alone, 4 REG up (largest +16, on a kernel at
142 -> 158) and 19 down; 51 STACK up and 46 down, largest +312 bytes on
`trc_pc_k1222` (+0.9%).

**Compile time: unchanged.** The generated kernels grow by comment lines and
one dummy only, no extra executable statement (`trc_pc_k2222.F90` 20 270 ->
20 313 lines). Interleaved nvfortran timings on a quiet box: (dd|dd) 562.4 ->
547.3 s, (dp|dp) 132.5 -> 132.8, (dd|ps) 64.6 -> 63.6, Rys (ff|ff) 11.53 ->
11.61, Rys (dd|dd) 4.49 -> 4.48, (pp|pp) 6.12 -> 6.11.

**The energies move in the last bits, and it is the FUSION, not the weight.**
Relocating the weight from `jfac`/`kfac` onto the digestion's `wq` is exact:
both factors are powers of two, and a host build carrying the new kernels
with the weight split still in place reproduces the base energies bit for
bit, skeleton and off, for water, B2H6 and benzene. Merging the runs is what
reorders: the chunked Rys and warp-role kernels map threads to work items in
blocks whose boundaries depend on the launch's own extent, so a fused launch
adds into G in a different order. On the serial host build benzene's
converged energy moves 9.3e-12 Eh against a 1e-10 gate with the iteration
count unchanged (water 2.8e-14, B2H6 6.5e-13). The symmetry-off plan never
had a weight to fuse on and stays bit-identical. On the device the
before/after spread (<= 4.6e-13 Eh total, <= 2.4e-12 on the two-electron
part) is inside the same binary's own run-to-run spread (<= 2.8e-13 and
<= 1.25e-12), i.e. the atomic accumulation order was never reproducible
there.

`eri%nwork` is unmoved: `check_sym_fock`'s pins are exact (water
14627/9968, B2H6 360415/107763, benzene 4741001/1342393, naphthalene
20671728/5451483).

**Open question 2 is answered and closed**: `real(dp)`, as W0 recommended.

**Note for W3/W4**: the tiled builder is the one in use and was measured, not
assumed. Benzene's device profile carries `trc_mdj_*` kernels, so `jfac == 0`
and the bin call is exchange-only and tiled (`bin_tiles_plan` ->
`walk_skeleton`). Both skeleton builders fill `sW` through the same
`bins_plan_emit`, so neither needed a change.

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
