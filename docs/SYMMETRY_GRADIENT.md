# Skeleton nuclear gradient: design (S16)

Status: **implemented** (see section 11 for what survived contact with the
code and what did not). Sections 0-10 are the design spike as it was signed
off; their measurements come from a CPU prototype (`test/spike_symgrad.F90`,
registered `DISABLED`) and from `nsys` on a GPU build of `main` at
**0f59b51**. Section 11 is the implementation: what was built, the realised
numbers, and the places the design was wrong.

S16 lists three optional follow-ups: DF-JK with symmetry, an XC grid
reduction, and skeleton gradients. This document recommends **skeleton
gradients**, says what it costs, and gives the numbers that rank it against
the other two (section 8).

The two-electron nuclear gradient is still a full C1 build under symmetry.
`eri%fock` takes `skeleton=`; `eri%gradient` takes no such argument, so
`grad_bins` builds `plan%build(..., use_density=.false.)` -- the full,
unweighted plan -- however symmetric the molecule is. Measured on a V100:
with `--symmetry auto` the Fock kernels of a pentacene/cc-pVDZ gradient step
fall from 1940 ms to 665 ms, and the gradient's own kernels stay at 441 ms
and 445 ms. The gradient goes from 17% to **35%** of all GPU kernel time in
the step, untouched.

## 0. The premise, checked against the code

| Claim | Verdict | Where |
|---|---|---|
| `eri_gradient` takes no `skeleton` argument | **true** | `src/eri/trc_eri.F90:1190`; the argument list is `(this, b, density, grad, k_scale, j_scale, omega)`. `eri_fock` (`:706`) and `eri_fock_resident` (`:1422`) both take `skeleton`. |
| so the gradient integral build runs the full quartet set | **true** | `grad_bins` calls `plan%build(b, launch%thresh, quartets%bin_pbase, .false.)` (`src/eri/trc_grad.F90:102`) -- the full plan, and with no bin-level density screen either. The only `build_skeleton` call site in `src/` is `trc_binkernel.F90:269`; the gradient has none. (`test/check_sym_fusion.F90:121` calls it directly to pin the launch count.) |
| and symmetry therefore does nothing for the gradient | **true, measured** | The gradient's GPU kernel time is 24.99 -> 24.16 ms (benzene), 70.69 -> 71.74 (naphthalene), 440.67 -> 445.30 (pentacene), 19.85 -> 19.48 (SF6) between `--symmetry off` and `auto`: the same to +-1%, i.e. run-to-run noise. Section 8. |
| CPHF/CPKS and the geometry optimiser are out of scope, "those paths run in C1" | **true** | `symmetry-plan.md`'s scope list and S17 decision 3; `trc_capi.F90:5134` (polarizability: `TRC_ERR_UNSUPPORTED` on an explicit request), `:4873` (Hessian), `:5356` (RI-MP2), `:4658` (the RI-MP2 gradient check). `--opt` calls `trc_symmetry_suppress_default`. A *default* context repeats the SCF in C1 instead of refusing; `check_capi_sym` pins both behaviours (`test/check_capi_sym.F90:481-495` for the refusals, `:779-828` for the C1 fallbacks, including that the symmetry result is gone afterwards). |
| DF is separately refused under symmetry | **true** | `trc_capi.F90:3838`, "symmetry is not supported with density fitting yet"; under the S17 default `--df` runs in C1. |
| XC matters only for DFT | **true**, and it is the larger item *for DFT* | Section 8: XC is 48% of a B3LYP naphthalene gradient step's GPU kernel time against the 4c gradient's 10%. It is still ranked second, for the reasons there. |
| so symmetry accelerates the part of a gradient step that is not the bottleneck | **true, and it is now the bottleneck** | After symmetry, the 4c gradient is 18% (benzene), 26% (naphthalene) and **35%** (pentacene) of a gradient step's GPU kernel time, up from 16-17% with symmetry off, because everything around it got faster and it did not. |

Two things the premise did not mention, both found by measurement and both
load-bearing for the design:

- the generated gradient kernels **cannot express a skeleton segment at all**
  (they take no `sSA`/`sSB`), so this is a generator change, not a plan change
  (section 2);
- the gradient's launch profile is **not** the one S15 was burned by: it is
  one screen and one compact launch per *class*, 42-46 for the whole gradient,
  at 0.08-15.7 ms each, and a weighted plan adds none (section 9).

## 1. The math

Notation as `docs/SYMMETRY_FOCK.md` section 1. `G` is the abelian point group,
`|G|` its order, and every operation is a diagonal sign matrix
`S_g = diag(s_x, s_y, s_z)` in the standard frame, so `g` is an involution
(`g = g^-1`). `perm(a, g)` is the atom `g` sends atom `a` to
(`trc_sym_geom_t%perm`); `R(g)` is the AO signed permutation
(`ao_map`, `ao_sign` of `trc_sym_ao`).

A shell quartet `q = {p, r}` (an unordered pair of shell pairs) contributes to
the forces on up to four centres. Write that contribution as a whole
`3 x natm` array `F_q(D)`: the derivative of the quartet's energy term with
its eight permutational orderings and its density weighting
(`Gamma = jfac D_ab D_cd - kfac/4 (D_ac D_bd + D_ad D_bc)`, `gen_rys_grad.py`),
supported on the at most four atoms of `q`. The full two-electron gradient is
`grad = sum_{q in S} F_q` over the screened quartet set `S`.

**Equivariance.** Rotating a quartet by `g` rotates its force contribution and
moves it to the image atoms:

    F_{g q}(:, perm(a, g)) = s_g * F_q(:, a)

(the force on the image atom is the sign-flipped force; `d/dR_{g A} = S_g d/dR_A`
because `S_g^-1 = S_g`). Define the **force symmetrizer**

    (P_g X)(:, a) := s_g * X(:, perm(a, g))

-- which is exactly `sym_grad_symmetrize` of `src/symmetry/trc_sym_scf.F90`
summed over one `g`. Then `F_{g q} = P_g F_q`, and `P` is a representation:
`P_g P_h = P_{gh}` (`s_g s_h = s_{gh}` componentwise and `perm(perm(a,g),h) =
perm(a, gh)` for an abelian group), so `(1/|G|) sum_g P_g` is an orthogonal
projector. **Condition**: `D` is invariant (`R D R^T = D`), which is what the
SO-blocked SCF produces and what the Fock skeleton already requires. `F_q` is
*quadratic* in `D`, not linear as `Dig_q` is, but the argument needs only the
invariance of `D` and the covariance of the integral derivative, so nothing
changes.

**The skeleton restriction and its weight.** For a `G`-closed set `S` with
representatives `Q0` and weights `w_q = |orbit(q)|`,

    grad = sum_{q in Q0} sum_{g in G/Stab(q)} P_g F_q
         = (1/|G|) sum_g P_g [ sum_{q in Q0} w_q F_q ]

(as in the Fock case: `P_g F_q = F_q` for `g in Stab(q)`, so the coset sum is
`1/|Stab|` of the full group sum, and `w_q = |G|/|Stab(q)|`). So

    grad_skel := sum_{q in Q0} w_q F_q,     grad = (1/|G|) sum_g P_g grad_skel.

**This is the same restriction and the same weight as the Fock build.** The
only difference is the symmetrizer: the Fock case averages a matrix,
`R G R^T`, at `O(nao^2 |G|)`; the gradient case averages a `3 x natm` array at
`O(3 natm |G|)` -- microseconds, and `sym_grad_symmetrize` already does it.

**At bin-pair granularity (scheme 2a).** The production skeleton plan cuts one
side of each admitted bin pair. For two different, orbit-closed bins `A`, `B`:

    sum_{p in A, r in B} F(p,r)  =  (1/|G|) sum_g P_g [ sum_{p in reps(A)} w_p sum_{r in B} F(p,r) ]

because for fixed `g` the map `r -> g r` is a bijection of `B`, so
`sum_r F(g p, g r)` depends only on `g p`, and summing over `g` hits each
member of `orbit(p)` exactly `|Stab(p)|` times. A bin against itself is the
rectangle `reps(A) x A` with an extra `1/2`, exactly as in the Fock case: the
kernels apply the `ab = cd` half only on a triangle
(`same_pair = sD(seg) .and. (iab == icd)` in every generated gradient kernel),
so the rectangle's `p = r` terms are not halved by the kernel and the
segment's half weight supplies the factor.

**What to test (the invariant).** For an invariant density,

    max | sym( grad_skel ) - grad_full |  <=  tol,     sym = (1/|G|) sum_g P_g

with `grad_skel` the weighted representative sum and `grad_full` what
`eri%gradient` returns today. Two corollaries worth asserting separately:

- `sum_a grad_skel(:, a) = 0`. Each `F_q` is translation invariant on its own
  (the kernels already exploit it: a quartet whose four shells sit on one atom
  returns immediately), so the *unsymmetrized* skeleton gradient is
  translation invariant too. Measured: 3e-15 to 9e-14 (section 4). This is a
  necessary but weak check -- it holds for a wrong weight as well.
- `sym(grad_skel) != grad_skel`. The raw skeleton gradient is **not** a
  physical gradient: measured, it differs from the full one by 0.58 to 1.8,
  against components of 0.27 to 2.1 (section 4). The symmetrization is doing
  all of the work, so a test that forgets it fails loudly rather than
  subtly.

**Where to symmetrize.** `trc_scf_gradient` already calls
`sym_grad_symmetrize` once on the finished gradient when `sym_geom` is present
(`src/scf/gradient/trc_scf_grad.F90:451`). Because the projector is linear and
every other term (one-electron, overlap, nuclear repulsion, XC on a symmetric
grid) is already equivariant, `P(skel + equivariant) = P(skel) + equivariant`:
that single existing call is in the right place and needs no change. **One
exception, and it is a decision** -- see section 7, item 3: `trc_capi.F90`
deliberately does *not* pass `sym_geom` for a **default DFT** run (the
input-axes grid breaks the symmetry at its own first-order error, so
projecting would move the gradient away from the `symmetry=off` one). A
skeleton 4c gradient cannot be left unprojected. The clean fix is to project
the **4c term alone**, inside `two_electron_4c`, where the skeleton is formed;
the end-of-routine projection then stays exactly as it is, including its
default-DFT suppression.

## 2. What the plan already provides, and what it does not

### Reusable unchanged

| What | Where | Note |
|---|---|---|
| Pair orbit tables | `trc_sym_pairs`, `trc_bin_orbits_t` (`src/eri/trc_eri_data.F90`) | Built by `eri%build(..., sym=geom)` into `eri%orb`; orbit-closed bins, representatives first in sections of equal orbit size. Already built on every symmetric engine, including the one the gradient runs on. |
| The weighted plan | `bins_plan_build_skeleton` / `skeleton_walk` / `bins_skel_cut_first` (`trc_bins.F90`) | Produces `sA/sNB/sOA/sOB/sPA/sPB/sSA/sSB/sD/sLA..sLD/sM` and `sW`. The gradient can call it instead of `plan%build`. |
| Per-segment weight | `trc_bin_plan_t%sW` | Already a mapped device array that the Fock kernels read per segment (W2 of `.claude/context/segment-fusion-plan.md`). Half-integer multiples; realised set `{1/2, 1, 2, 4, 8}`. |
| Class sort | `bins_plan_sort_by_class` | Folds `nint(2*sW)` into the key when `sW` is allocated, so one weight occupies a contiguous run. The gradient already calls it (`GCLASS_RADIX = 11`). |
| Invariant screen bound | `skeleton_dsh` (`trc_eri.F90`, private) | `dsh := max over each pair's orbit`, so the surviving quartet set is `G`-closed. The gradient's in-kernel screen reads `dsh`, so it needs this; the routine exists and is called by the Fock entries. |
| Force symmetrizer | `sym_grad_symmetrize` (`trc_sym_scf.F90`) | Already the exact `(1/|G|) sum_g P_g`, already called by `trc_scf_gradient`, already tested on a synthetic non-symmetric gradient by `check_sym_c2v`/`_c2`/`_c2h`/`_ci`/`_d2`/`_d2h`/`_tops`/`_scf_cs`. |
| Completed pair data | `complete_pairs` (`trc_eri.F90`) | `eri%gradient` already calls it, so the unique-pair build (U1) is a non-issue. `check_sym_unique` already pins `eri%gradient` on a unique-pair engine against a plain one. |
| `bins_rep_weights` | `trc_bins.F90` | Per-pair weights for a cut *by pair* rather than by bin segment. The MD engines use it; the bin-segment gradient does not need it. |

### Not reusable: the gradient kernels cannot express a skeleton segment

**This is the blocking finding of the spike, and it was measured, not reasoned
about.** The generated gradient kernels take **no `sSA`/`sSB`**. They index
the bin-interleaved primitive tables with the segment's own extent:

    kp = sPA(seg) + (k1 - 1)*sA(seg)  + iab      ! trc_pc_rysg<tag>.F90
    kq = sPB(seg) + (k2 - 1)*sNB(seg) + icd

where every energy kernel uses the bin stride instead:

    kp = sPA(seg) + (k1 - 1)*sSA(seg) + iab      ! trc_pc_rys<tag>.F90

That is correct today only because `plan%build` sets `sSA = sA` and
`sSB = sNB`: a gradient segment always covers whole bins. Every skeleton
segment is a *sub-range* of its bin (`sA = sec_cnt(r) < sSA = bin_cnt`), so the
kernels read the wrong primitives. The *screen* pass is unaffected -- it
indexes pairs, not primitives -- so the damage is confined to the compact
pass, which is where 81-89% of the gradient's GPU time is (section 8).

The prototype shows exactly this: feeding
the production skeleton plan to the existing kernels gives a symmetrized
gradient wrong by **16.6** (water), **3.4** (B2H6), **2.3e+05** (SF6) and
**1.8e+03** (benzene) against components of order one. Re-running the same
quartets with a stride-one primitive layout makes it exact to 5e-15 .. 2e-11
(section 4). No weight, no plan trick and no host-side repacking avoids this:
it is a kernel-signature change.

The same gap is why **exchange tiles are not available to the gradient**
either (`trc_bin_tiles` emits sub-range segments). Adding `sSA`/`sSB` unlocks
both.

### What a weighted gradient plan needs

1. `scripts/gen_rys_grad.py`: thread `sSA(nseg)`, `sSB(nseg)` and `sW(nseg)`
   through `*_compact` and its `crhf`/`cuhf` item routines; replace `sA`/`sNB`
   by `sSA`/`sSB` in the three primitive-slot expressions (`kp`, `kq` and the
   `qp_c(sPB(seg) + icd)` primitive bound); and one line in the item routine's
   prologue, next to the existing degeneracy factors:

   ```fortran
   wq = 1.0_dp
   if (same_ab) wq = wq*0.5_dp
   if (same_cd) wq = wq*0.5_dp
   if (same_pair) wq = wq*0.5_dp
   wq = wq*sW(seg)            ! the orbit weight, one load and one multiply
   ```

   `wq` already exists and already multiplies the digestion (`cj = 8 wq jfac`,
   `ck = 8 wq kfac_eff`), so the weight rides on the integral exactly as W2
   put it on the Fock side. `jfac`/`kfac` stay unweighted scalars, for the same
   reason as there.
2. The **screen** pass needs **nothing**. It decodes `(iab, icd)` from the
   segment's own extent `sA(seg)` -- correct for a sub-range -- and reaches
   the pairs only through `sOA`/`sOB`, never the interleaved primitive tables
   (it takes `sPA`/`sPB` as unused dummies and no `qp_*` array at all). Its
   bound is `sp_q * sp_q * dsh`, unweighted, which is exactly right: with a
   symmetrized `Q` and `dsh` a representative and all its images pass or fail
   together, so dropping the representative once accounts for its whole orbit.
   This matches what `count_kept` does on the Fock side.
3. `scripts/gen_rys_grad_dispatch.py`: `rysg_compact_dispatch`'s signature
   (`rysg_screen_dispatch` is unchanged).
4. `src/eri/trc_grad.F90`: `plan%build_skeleton(...)` when asked, map and pass
   `plan%sSA`, `plan%sSB`, `plan%sW` (filling ones when the plan carries no
   weights, as `binkernel_fock_bins` does), everything else unchanged.
5. `src/eri/trc_eri.F90`: a `skeleton` argument on `eri_gradient` with
   `eri_fock`'s semantics, a `skeleton_dsh` call, and the `orb%built` guard.
6. `src/scf/gradient/trc_scf_grad.F90`: pass it, and project the 4c term
   inside `two_electron_4c` (section 1).
7. Generators for the *energy* kernels, the MD engines, `trc_bins`,
   `trc_bin_tiles`, `trc_binkernel`, `trc_scf_fock`, the MPI split and the
   incremental Fock: **unchanged**.

## 3. Counts, measured

Host build, gfortran 15.1.0, `TERCO_LMAX=2`, commit **0f59b51**, cc-pVDZ,
Cartesian caller space, `tau = 1e-10`, the engine built with the detected
group, a group-projected pseudo-random density (`max|D - R D R^T| <= 1.1e-16`
on every system; the full gradient is already group symmetric to 2e-16 ..
5e-14, which is the prototype's own sanity gate). "Work items" is what a
launch evaluates before the in-kernel screens -- `plan%nwork` of the plan
`grad_bins` builds; "survivors" is what the gradient's own screen pass writes
to its survivor list, counted by running `rysg_screen_dispatch` and reading
its atomic counter.

**There is no `TRC_KEPT_STATS` for the gradient, and no counter at all.**
`eri%gradient` sets neither `eri%nwork` nor `eri%nkept` (only the `fock*`
entries do), and no test sets `TRC_KEPT_STATS`. The survivor counts below
come from the prototype; a `TRC_KEPT_STATS`-equivalent on the gradient would
be two lines (`plan%nwork`, and the sum of the screen pass's `ctr`) and is
recommended as part of S16 (section 7, item 6).

### 3.1 Work items the gradient launches

| system | group | natm | nao | full | skeleton | reduction |
|---|---|---|---|---|---|---|
| water | C2v | 3 | 25 | 14 627 | 9 968 | **1.467x** |
| b2h6_d2h (nonplanar) | D2h | 8 | 60 | 360 415 | 103 226 | **3.492x** |
| tops_sf6 | D2h | 7 | 109 | 2 883 269 | 1 192 250 | **2.418x** |
| benzene (D6h as D2h) | D2h | 12 | 120 | 4 725 138 | 1 337 126 | **3.534x** |
| naphthalene_d2h | D2h | 18 | 190 | 20 671 728 | 5 451 483 | **3.792x** |
| pentacene_d2h | D2h | 36 | 400 | 150 479 149 | 39 110 378 | **3.848x** |

Nothing about the gradient changes the quartet arithmetic -- it is the same
plan -- and the counts bear that out: the full counts reproduce
`docs/SYMMETRY_FOCK.md` section 2.2 **exactly** for water (14 627), B2H6
(360 415), naphthalene (20 671 728) and pentacene (150 479 149), and the
skeleton counts exactly for water (9 968), naphthalene (5 451 483) and
pentacene (39 110 378). Benzene's full count is 0.3% lower than S13's
(4 725 138 against 4 741 001) and B2H6's skeleton count 4.2% lower
(103 226 against 107 763); those two were not chased, and the candidates are a
different benzene geometry file and the difference between S13's CPU
prototype's "cut the smaller side" and production
`bins_skel_cut_first`'s `nrep(A) cnt(B) <= cnt(A) nrep(B)`. Note also that the
gradient's plan carries **no** bin-level density screen
(`plan%build(..., use_density=.false.)`: `Gamma` is not bounded by the Fock
density's bin bound, same reasoning as `eri_fock_nosym`), so its work-item
counts are not in general the Fock build's.

### 3.2 Survivors of the gradient's own screen

| system | full | skeleton | reduction | survivor fraction |
|---|---|---|---|---|
| water | 11 534 | 6 605 | 1.746x | 79% / 66% |
| b2h6_d2h | 350 878 | 98 871 | 3.549x | 97% / 96% |
| tops_sf6 | 2 761 476 | 1 130 616 | 2.442x | 96% / 95% |
| benzene | 4 609 837 | 1 303 367 | 3.537x | 98% / 97% |
| naphthalene_d2h | 20 113 684 | 5 305 368 | 3.791x | 97% / 97% |
| pentacene_d2h | 144 141 689 | 37 493 846 | 3.844x | 96% / 96% |

The reduction survives screening, as on the Fock side: a symmetrically
screened set is `G`-closed. Water's 1.75x against 1.47x is the one place the
skeleton does *better* after screening, and it is small-number noise (6.6k
quartets), not an effect to rely on.

### 3.3 Segments and the bra/ket swap

| system | full segs | skeleton segs | swapped segs | highest l swapped |
|---|---|---|---|---|
| water | 909 | 1 320 | 575 | 2 |
| b2h6_d2h | 4 474 | 6 275 | 2 525 | 2 |
| tops_sf6 | 7 304 | 11 787 | 5 116 | 2 |
| benzene | 3 106 | 5 090 | 1 706 | 2 |
| naphthalene_d2h | 5 186 | 6 553 | 2 315 | 2 |
| pentacene_d2h | 5 163 | 7 304 | 2 338 | 2 |

"Swapped" counts segments with `(sLC, sLD) > (sLA, sLB)`: `bins_skel_cut_first`
cut on the side with fewer representatives, which is the far side of the
admitted bin pair. A quarter to a half of all skeleton segments are swapped.
At `TERCO_LMAX=2` all 81 gradient classes are built, so none of them lacks a
kernel (the prototype checks `rysg_has_class` on every one: 0 missing). Above
d they would be the classes the gradient's filter drops -- section 5.

### 3.4 What forbidding the swap costs

Cutting always on the *first* bin of the admitted pair -- the full plan's own
bra -- keeps the full plan's class set exactly, so no class the gradient has no
kernel for can appear and no new class run is created. Measured:

| system | skeleton (swap allowed) | skeleton (swap forbidden) | reduction, forbidden | work penalty |
|---|---|---|---|---|
| water | 9 968 (1.467x) | 11 000 | **1.330x** | +10.4% |
| b2h6_d2h | 103 226 (3.492x) | 120 632 | **2.988x** | +16.9% |
| tops_sf6 | 1 192 250 (2.418x) | 1 365 994 | **2.111x** | +14.6% |
| benzene | 1 337 126 (3.534x) | 1 377 969 | **3.429x** | +3.1% |
| naphthalene_d2h | 5 451 483 (3.792x) | 5 565 270 | **3.714x** | +2.1% |
| pentacene_d2h | 39 110 378 (3.848x) | 39 546 981 | **3.805x** | +1.1% |

The penalty shrinks with size (it is paid on bin pairs whose two sides have
very different representative counts, a small-molecule effect). Section 5
shows it does not have to be paid at all.

## 4. Exactness, measured

The prototype has no weighted gradient kernel to call, so it applies the
weight by **repetition**: a segment of weight `w` is emitted `nint(2w)` times
(every weight is a multiple of one half) and the sum is halved. That is an
exact stand-in for a per-segment `sW` -- the kernels are linear in the
accumulation -- at about twice the full build's cost, and it changes nothing
in `src/`. To get past the missing `sSA`/`sSB` (section 2) the prototype
re-lays the primitive tables out with one pair's primitives contiguous and
emits **one segment per quartet**, where `sA = sNB = 1` *is* the stride, so the
production kernels are exact on it. Its full-plan run reproduces
`eri%gradient` to 1.3e-15 .. 6.2e-14, which is what makes the comparison
meaningful.

Host build, serial, commit 0f59b51, cc-pVDZ Cartesian, `tau = 1e-10`,
group-projected random density:

| system | group | max\|grad_full\| | raw skeleton err | **symmetrized err** | relative | trans. inv. of raw skeleton |
|---|---|---|---|---|---|---|
| water | C2v | 1.507 | 0.830 | **4.66e-15** | 3.1e-15 | 3.4e-15 |
| b2h6_d2h (nonplanar) | D2h | 0.272 | 0.993 | **1.73e-14** | 1.7e-14 | 3.9e-14 |
| tops_sf6 | D2h | 0.840 | 0.580 | **2.06e-12** | 2.1e-12 | 5.0e-14 |
| benzene | D2h | 2.139 | 1.835 | **2.17e-11** | 1.0e-11 | 8.8e-14 |

The raw (unsymmetrized) skeleton gradient is wrong by order one everywhere:
the projection is not cosmetic. Benzene's 1e-11 is the double-precision
summation noise of ~9.4M accumulated segment contributions (the prototype
accumulates about twice the full build's work items because of the repetition
trick; a real weighted kernel accumulates 3.5x *fewer* than the full build and
will be quieter). So the oracle's gate should be **relative to
`max(1, max|grad|)`**, as `check_sym_fock`'s `TOL_REL` is; see section 6.

The same table, run against the **production** skeleton plan with the
kernels as they are (whole-bin stride), for the record: 16.6, 3.4, 2.3e+05,
1.8e+03. That is the `sSA`/`sSB` finding.

## 5. The class filter: the blocker, and a better answer than the two on offer

`src/CMakeLists.txt`'s `_rys_grad_classes` filter builds, above d, only the
angular-momentum-canonical orderings `la >= lb`, `lc >= ld`,
`(la, lb) >= (lc, ld)`. The comment on `main` (0f59b51) records why the third
condition is not safe even today: `split_large_bins` and `arrange_orbits`
append a split bin's extra chunks past every original index while keeping the
parent's `(la, lb)`, so the full plan's ascending-bin-index walk can present a
low-l chunk as the bra against a high-l ket, and `grad_classes` `error stop`s
on a class it has no kernel for. Untested, because no gradient test sets
`TRC_BIN_CAP`.

Kernel counts, read from the configure step:

| build | gradient class kernels built | pair-canonical classes a plan can emit | missing |
|---|---|---|---|
| `TERCO_LMAX=2` | 81 | 36 | **0** |
| `TERCO_LMAX=3` | 115 | 100 | **30** |

(`TERCO_LMAX=4` is the same as 3: `_grad_lmax` is capped at f and `trc_capi`
refuses a g basis up front.) The 30 missing tags are exactly the
f-containing ones with `(lc, ld) > (la, lb)`:

    0030 0031 0032 0033 1030 1031 1032 1033 1130 1131 1132 1133
    2030 2031 2032 2033 2130 2131 2132 2133 2230 2231 2232 2233
    3031 3032 3033 3132 3133 3233

(the CMake comment's "~42" is the *energy* side's bra<ket count; the gradient
number is 30.)

A skeleton gradient makes this live rather than latent: `bins_skel_cut_first`
swaps bra and ket *by design* whenever the far side has fewer
representatives, up to `SKEL_SWAP_LMAX = 3`, and section 3.3 measures that as
a quarter to a half of all skeleton segments. At `TERCO_LMAX=3` with an f
basis those swaps would hit the 30 missing classes.

Three ways out, in increasing order of how much they are worth:

**(a) Build the 30 extra gradient kernels.** Keeps the cheapest-side cut and
the full work reduction. Costs 30 of the heaviest gradient kernels
(f-containing) of compile time at `TERCO_LMAX=3`, forever, and raises the
plan's distinct class count, which raises the launch count: the prototype
counts 21-22 class runs for the full plan against 32-36 for the swapping
skeleton plan, i.e. **+60%** launches per pass. That is the S15 pattern
(section 6 of `.claude/context/segment-fusion-plan.md`) reintroduced for no
necessary reason.

**(b) Forbid the swap** (always cut on the admitted pair's first bin). No new
kernels, class set identical to the full plan's, launch count unchanged (the
prototype's "swap forbidden" column is this scheme). Costs 1.1% to 16.9% of
the work reduction (section 3.4). It inherits the full plan's reliance on the
bin ordering, so it does not close the `split_large_bins` hole.

**(c) Recommended: cut the cheaper side, but emit the segment
angular-momentum-canonically.** The side that is *cut* and the side that is
written as the segment's *bra* are independent once the kernels take
`sSA`/`sSB`, because a sub-range is then expressible on either side -- that is
exactly what the exchange-tile mechanism does on the energy path -- and
because **`sW` is a scalar per segment, uniform over its quartets, so it does
not care which side the representatives are on**. So: pick the cheaper side to
cut (keeping the full 1.47x-3.85x), then order the two sides of the emitted
descriptor so that `(sLA, sLB) >= (sLC, sLD)`. Then

- no class outside the filter's own built set (`la >= lb`, `lc >= ld`,
  `(la, lb) >= (lc, ld)`) is ever emitted -- by construction, not by
  coincidence of the bin ordering, which is what the current argument relies
  on: no new kernels, no new class runs, no launch growth;
- `SKEL_SWAP_LMAX` becomes unnecessary on the gradient path;
- applying the same canonicalisation in `bins_plan_build` **closes the
  pre-existing `split_large_bins` hole** the `_rys_classes` comment records,
  which is worth doing whether or not the skeleton gradient is built;
- it costs nothing in work items (the quartet set is unchanged; only which
  side of the descriptor the data goes in moves).

The one thing (c) needs that (a) and (b) do not is for the kernels to accept a
sub-range on the **ket** side as well as the bra -- which `sSB` already
expresses and which the energy kernels already do.

## 6. The oracle

`check_sym_fock` is the model; `check_sym_scf_cs`'s `grad_case` and the
existing finite-difference checks supply the rest. Proposed
`test/check_sym_grad.F90` (test-drive, `run_testsuite(..., parallel=.false.)`),
with the axes `check_sym_fock` sweeps:

1. **The central identity.** `sym(grad_skel)` against `grad_full` from the
   *same* engine, for every group's test molecules x Cartesian/spherical x
   `{J+K, J only (k_scale=0), K only (j_scale=0), long-range
   (j_scale=0, k_scale, omega>0)}` x `{RHF (ndens=1), UHF (ndens=2)}`.
   Gate: `max|.| <= TOL_REL * max(1, max|grad_full|)` with
   **`TOL_REL = 1e-10`** -- `check_sym_fock` uses 5e-11 for `G`; benzene's
   measured 1.0e-11 (section 4, and that is the noisier repetition
   prototype) leaves a factor of ten. State the measurement in the
   declaration's own `!!` block, as `check_sym_fock.F90:87` does.
2. **The skeleton must actually be smaller.** `nwork` full against skeleton,
   `ns < nf`, so a silent fall back to the full plan fails
   (`check_sym_fock.F90:544`'s guard). This needs the counter of section 7
   item 6 to exist.
3. **The raw skeleton must be wrong.** `max|grad_skel - grad_full| > 0.1 *
   max|grad_full|` on at least one molecule, so a test that accidentally
   compares two full builds fails. `check_sym_d2h` already has the analogous
   `moved > 0.1` assertion for the projector.
4. **Translational invariance of the raw skeleton**, `<= 1e-10`, as
   `check_full_grad_fd:144` and `check_rys_grad_fd:263` assert for the full
   one.
5. **The forbidden components.** `check_sym_scf_cs:405`'s
   `forbidden <= TOL_GRAD = 1e-10` on the *raw* gradient of a symmetric run
   (an atom that is its own image under an operation has no force along the
   flipped axes) and `asym <= 1e-14` on the symmetrized one, now with the
   skeleton path on. The same three assertions already exist in
   `check_sym_c2v`, `_c2`, `_c2h`, `_ci`, `_d2`, `_d2h`, `_tops`; turning the
   skeleton on there is the cheapest coverage in the project.
6. **Against the energy path, end to end.** `check_full_grad_fd`'s and
   `check_rys_grad_fd`'s finite differences with the skeleton on: the
   analytic skeleton gradient against a central difference of terco's own
   reconverged SCF energy (`RELTOL = 1e-5`/`ABSTOL = 5e-7`) and of
   `Tr(D vhf(D))` at fixed `D` (`RELTOL = 1e-7`). These are the only checks
   that catch a *consistent* error in both the skeleton and the full path.
   `check_rys_grad_fd` runs its reference through HGP
   (`TRC_RYS_CLASSES=none`), which is an independent engine -- keep that.
7. **The C1 agreement**, `maxval(abs(gsym - g1)) <= 1e-9` on the same
   symmetrized geometry (`check_sym_scf_cs:409`), with the skeleton on.
8. **`TERCO_LMAX=3`**, an f basis (cc-pVTZ), **and `TRC_BIN_CAP` lowered** --
   the combination no gradient test has ever run. This is the only test that
   can fire the class-filter `error stop`, and after the fix it is the test
   that proves it cannot. `check_sym_fock_cap` (`TRC_BIN_CAP=40`) is the
   pattern; `check_binfock_split` (`TRC_BIN_CAP=3`) is the other.
9. **Exchange tiles on the gradient**, once `sSA`/`sSB` exist: an extended D2h
   molecule with a banded density (`check_sym_fock`'s `banded_density`,
   pentacene) and the tiled gradient against the untiled one.
10. The OpenMP build at 4 threads, and `symmetry=off` bit-for-bit (the
    goldens) -- the project's standing requirements.

The prototype `test/spike_symgrad.F90` stays on the branch as the
counts/exactness harness; it is not the oracle (it is one segment per quartet
and twice the work).

## 7. Decisions needed, and the compile cost

1. **Scheme.** Confirm (c) of section 5: cut the cheaper side, emit
   angular-momentum-canonically. The alternative (b) is simpler by one
   comparison and costs 1-17% of the gain.
2. **The kernel signature.** `sSA`, `sSB` and `sW` on every gradient kernel's
   compact pass and on `rysg_compact_dispatch`. There is no way to have a skeleton gradient
   without `sSA`/`sSB` (section 2). Confirm that regenerating all 81/115
   gradient kernels is acceptable -- it is a signature change, not a size
   change, so the compile cost is a re-run of what the build already pays.
3. **Default DFT.** `trc_capi.F90:4453` withholds `sym_geom` from the gradient
   for a default DFT run on purpose. Options: (i) project the 4c term inside
   `two_electron_4c` and leave the rest as today (recommended, exact for the
   4c term, no change to the DFT grid story); (ii) run the full 4c gradient
   for default DFT runs only. (i) is a few lines and strictly better.
4. **Also canonicalise `bins_plan_build`?** It closes the latent
   `split_large_bins` `error stop` in the gradient path today, independent of
   the skeleton. Recommended; it is three lines and one test
   (`TRC_BIN_CAP` + `TERCO_LMAX=3` + an f basis).
5. **Should the energy side drop its 42 bra<ket Rys kernels** once the plan is
   emitted canonically? That is a compile-time saving on the energy path
   (`src/CMakeLists.txt`'s own note says those orderings exist only because
   the skeleton plan swaps). Out of scope here; worth a separate look.
6. **A gradient work counter.** `eri%gradient` reports nothing. Two lines
   (`this%nwork = plan%nwork`, and `this%nkept` from the screen pass's
   counter, behind `TRC_KEPT_STATS` or a `count_survivors`-style option) make
   the oracle's "the skeleton is smaller" assertion possible and give the
   gradient the diagnostic the Fock path has had since S15.
7. **Range separation.** A range-separated hybrid runs `eri%gradient` twice
   (`j_scale=0, k_scale=exx_lr, omega`). Both passes take the same plan, so
   both are skeleton; the projector is applied once to the sum. No extra work,
   but it is an oracle axis (item 1).

### Compile cost

The change adds three array arguments and one multiply to each gradient
kernel; the emitted code does not grow. The cost is therefore **a re-run of
what the build already pays** for the gradient kernels, once -- not an
addition. Measured on this host (NVHPC 26.5,
`-stdpar=gpu -acc -gpu=mem:separate,cc70`, 40 cores, `nvidia-smi` idle),
recompiling all 81 `trc_pc_rysg*.F90` of the `TERCO_LMAX=2` GPU build with
their own compile commands, 20 at a time:

- **277.9 s wall at -j20**, 1473.4 s of summed single-file compile time
  (24.6 min serial), 0 failures;
- heavily skewed: `rysg2220` 240.0 s and `rysg2121` 177.9 s are 28% of the
  total on their own, then `rysg2210` 71.6 s, `rysg2120` 68.1 s,
  `rysg1111` 42.4 s; the mean is 18.2 s.

Relink is on top (NVHPC refuses to relink an existing `libterco.so`, so
delete it first -- `gpu-pitfalls.md`). A *full* `make terco` after deleting
those 81 objects took 8:04 wall at -j20 on this host, but that also rebuilt
the ~100 `src/` objects the deleted module files invalidated plus the link, so
4:38 is the number for the kernels themselves.

At `TERCO_LMAX=3` the same set is 115 kernels, and the 34 extra ones are the
f-containing ones, which are the expensive end of that distribution.
**Option (a) of section 5 would add 30 more f-containing gradient kernels,
permanently.** Extrapolating from the measured distribution -- `rysg2220`
(ltot 6) alone is 240 s and the new tags run to `rysg3233` (ltot 11) -- that
is somewhere between ~10 and ~100 minutes of added serial compile, i.e.
roughly +1 to +10 minutes of wall at -j20. That is an estimate from the
emitted sizes, not a measurement (no `TERCO_LMAX=3` GPU build was made), and
it is one of the reasons to prefer (c). **Options (b) and (c) add no
kernels at all.**

Host (gfortran 15.1.0, `-j16`) for comparison: a whole host build from
scratch, library plus every test program, took 3.5 minutes on this box
(03:04:35 to 03:08:05), of which `libterco.so` was 2.5 minutes. The gradient
kernels are not the host build's bottleneck.

## 8. Why this and not the other two S16 items

Measured on the same GPU build (V100-DGXS-32GB, cc70, nvfortran 26.5,
`TERCO_LMAX=2`, commit 0f59b51, cc-pVDZ, `--grad`, `--econv 1e-10`, GPU 0
otherwise idle, nsys 2025.1.3, `cuda_gpu_kern_sum` grouped by the class tag in
the kernel name). "Total" is all GPU kernel time of the whole run (SCF +
gradient, including linalg).

| system | symmetry | total | 4c gradient (`rysg*`) | gradient launches | Fock bins | Fock launches | gradient share |
|---|---|---|---|---|---|---|---|
| benzene | off | 157.1 ms | 24.99 ms | 42 | 68.14 ms | 495 | 15.9% |
| benzene | auto | 134.0 ms | **24.16 ms** | 42 | 49.65 ms | 603 | **18.0%** |
| naphthalene_d2h | off | 397.5 ms | 70.69 ms | 42 | 254.35 ms | 558 | 17.8% |
| naphthalene_d2h | auto | 276.0 ms | **71.74 ms** | 42 | 121.52 ms | 738 | **26.0%** |
| pentacene_d2h | off | 2533.4 ms | 440.67 ms | 46 | 1939.62 ms | 600 | 17.4% |
| pentacene_d2h | auto | 1285.6 ms | **445.30 ms** | 46 | 664.50 ms | 810 | **34.6%** |
| tops_sf6 | off | 191.8 ms | 19.85 ms | 42 | 65.26 ms | 924 | 10.4% |
| tops_sf6 | auto | 204.8 ms | **19.48 ms** | 42 | 85.54 ms | 1434 | **9.5%** |

The gradient's kernel time is the same to +-1% with symmetry on and off, on
every system -- the direct confirmation that `eri%gradient` ignores symmetry.
(SF6 is also the one system where symmetry still *loses* overall, 191.8 ->
204.8 ms: 7 atoms, 1434 Fock launches against 924. That is S15's small-molecule
launch problem, not a gradient question.)

Projecting the gradient's time by its measured quartet reduction (section 3,
and both passes scale with it -- the screen pass with work items, the compact
pass with survivors). The cost per survivor does not change: an orbit's
members share their shell types *and* their primitive count, because that is
what `bin_key` is keyed on and why an orbit never straddles two bins, so a
representative is statistically the same quartet as its images:

| system | step now (auto) | gradient now | gradient projected | step projected | step speedup | vs `off` |
|---|---|---|---|---|---|---|
| benzene | 134.0 ms | 24.16 | 6.8 | 116.6 ms | 1.15x | 1.35x |
| naphthalene_d2h | 276.0 ms | 71.74 | 18.9 | 223.2 ms | 1.24x | 1.78x |
| pentacene_d2h | 1285.6 ms | 445.30 | 115.7 | 956.0 ms | 1.34x | **2.65x** |
| tops_sf6 | 204.8 ms | 19.48 | 8.1 | 193.4 ms | 1.06x | 0.99x |

**DF-JK with symmetry is not an S16-sized item.** It is not an optimisation of
an existing path: symmetry and density fitting do not coexist at all.
`trc_capi.F90:3838` refuses the combination outright ("symmetry is not
supported with density fitting yet"), and under the S17 default `--df` runs in
C1. The DF-SCF forms J and K from the retained `B` tensor, not from quartets,
so "unique auxiliary-index blocks" means first giving the DF-SCF an SO-blocked
path and a symmetry-aware fit -- a project of S3-S14's size, with its own
gradient path (`trc_scf_grad_df`) to follow. Nothing in the bin/plan machinery
this document reuses applies to it.

**XC grid reduction is worth more than the gradient for DFT, and less
overall.** Same build, naphthalene_d2h/cc-pVDZ/B3LYP, `--grad --econv 1e-9`:

| component | symmetry off | symmetry auto |
|---|---|---|
| total GPU kernel time | 815.4 ms | 708.3 ms |
| XC (incl. grid build, 354 launches) | 345.0 ms (42.3%) | **340.7 ms (48.1%)** |
| Fock bins | 218.1 ms | 104.1 ms |
| 4c gradient (`rysg*`) | 70.7 ms (8.7%) | **73.3 ms (10.3%)** |
| MD J | 17.9 ms | 18.1 ms |
| other (linalg, 1e, transforms) | 163.7 ms | 172.1 ms |

So for a *DFT* step XC is 48% and the 4c gradient 10%; an |G|/2-fold XC
reduction would be the bigger win there. Three things rank it second anyway:

1. It does nothing for HF or for the four-centre part of a hybrid, which is
   where the 26-35% above lives.
2. It needs an **exactly** symmetric grid, and the S17 default deliberately
   does not build one: the grid sits on the user's axes
   (`dft_grid_options_t%angular_rotation`), and the plan records the residual
   as a 1e-7..1e-6 gradient error at `grid_level` 3 (2e-9 at 5). Pruning and
   the Becke partition have to be made exactly equivariant first -- the S16
   plan text already flags this ("the grid must then be symmetric, which
   Lebedev grids in the standard frame are, but pruning and Becke partitions
   must be checked"). Until that is done the reduction is available only to an
   *explicit* `--symmetry` request, i.e. not to the default path.
3. 166 of the 345 ms is a single `trc_xc_batch_batch_maxima` launch (grid
   construction, once per geometry). Some of the saving there is a grid-build
   optimisation rather than "integrate the unique atoms' points with weights",
   so the 3.6x the atom orbits suggest is an upper bound on part of the time,
   not on all of it.

The skeleton gradient, by contrast, is the only one of the three that changes
what a **geometry optimisation, a transition-state search or an MD step**
actually spends its time on in the default (HF or hybrid, four-centre)
configuration, and it is the one whose machinery -- orbit tables, weighted
plan, `sW`, the force projector -- is already built and already tested on the
Fock side.

## 9. Recommendation

Build the skeleton gradient: scheme 2a's plan (`bins_plan_build_skeleton`),
**emitted angular-momentum-canonically** (section 5c), with `sSA`, `sSB` and
`sW` threaded into the generated gradient kernels and the orbit weight folded
into the existing `wq`. Symmetrize the four-centre term alone, inside
`two_electron_4c`, with `sym_grad_symmetrize`.

Expected: **3.5x-3.85x fewer quartets on planar D2h acenes**, 3.49x on
nonplanar D2h (B2H6), 2.42x on SF6, 1.47x on C2v water, which on the measured
launch profile should give **1.15x (benzene), 1.24x (naphthalene), 1.34x
(pentacene)** on a whole HF gradient step, and take pentacene's gradient step
to 2.65x against `symmetry=off`. No gain on DFT-dominated steps beyond the
10% the 4c term holds there.

The GPU reality check says the launch profile has room, and says why this is
*not* the S15 trap:

- the gradient launches **one screen and one compact kernel per class**, not
  per segment and not per weight: nsys shows `Instances = 1` for every
  `trc_pc_rysg<tag>_{screen,compact}` kernel, 21 classes, 42-46 launches for
  the whole gradient, set by the class count and `GRAD_COMPACT_CAP`'s work-item
  batching. Segments and weights do not enter;
- those launches are large. Compact-pass mean per launch: **15.7 ms**
  (pentacene), 1.03 ms (benzene), 0.84 ms (SF6); screen pass 3.68 / 0.12 /
  0.08 ms. The smallest single class launch on the smallest system is 132 us.
  After a 3.5x cut the compact mean is ~300 us on benzene and ~4.5 ms on
  pentacene -- still far above the ~85 us per launch at which S15 measured the
  Fock build underfilling the card;
- the per-class distribution is flat (pentacene: largest class 13.8%, 16 of 21
  classes to reach 91%), so there is no single class whose shrinkage would
  leave the card idle;
- the launch **count** does not change under scheme (c). nvfortran's fixed
  per-offload-region cost (W1: ~33 us per launch of additive memset + HtoD) is
  ~1.5 ms over 46 launches, against 24-445 ms of gradient kernel time, and it
  stays there. Under scheme (a) -- allowing the swap and building the extra
  kernels -- the class count rises 21-22 to 32-36 per pass, which is the
  pattern W2 had to undo on the Fock side; that is the main reason to prefer
  (c).

## 10. Caveats and what the prototype does not show

- Every count and every exactness number is **host, CPU**, from the plan and
  from the production gradient kernels run serially. The projected step
  speedups in section 8 are the measured gradient time divided by the measured
  quartet reduction; they are **not** a measurement of a skeleton gradient,
  because no weighted gradient kernel exists. They assume the compact pass's
  time tracks its survivor count, which the flat per-class distribution and
  the large launch sizes make plausible but which only an implementation can
  confirm.
- The exactness table is **Cartesian caller space only**. The prototype's
  hand-rolled group projection of the density is not invariant in the
  spherical AO space (`max|D - R D R^T| = 0.27` there against 5.6e-17 in
  Cartesian), so the spherical runs were dropped rather than reported; the
  production symmetrizer (`trc_sym_scf_t%symmetrize`, which `check_sym_fock`
  uses) does not have this problem and the oracle must sweep both spaces. The
  force projector itself is basis-space independent, and the integrals are
  always formed in the Cartesian working basis, so nothing in the math is
  specific to one space -- but it is unverified here.
- The prototype applies the weight by repetition and runs one segment per
  quartet. Its accumulation order, and therefore its summation noise, is not
  the production one. Benzene's 1.0e-11 relative error is an upper bound on
  what a real weighted kernel would show.
- `TERCO_LMAX=2` only. The 30 missing f classes of section 5 are counted from
  the CMake filter (confirmed against the configure step's own
  "115 Rys gradient class kernels"), not observed firing: that needs an
  `TERCO_LMAX=3` build, an f basis and a lowered `TRC_BIN_CAP`, which is
  oracle item 8 and the reason it is in the list.
- No MPI, no multi-GPU. The gradient's rank split is by contiguous segment
  slice (`grad_classes`), whose balance is already unmeasured
  (`trc_grad.F90`'s own note); a weighted plan changes the per-segment work,
  so that balance would want a look before `TERCO_ENABLE_MPI=ON` runs with it.
  **Settled conservatively** by section 11 item 9: the driver declines the
  skeleton on more than one rank, so a multi-rank gradient runs the full
  plan and never reaches a weighted one. Lifting that restriction is what
  needs the balance measured.
- `compute-sanitizer` was not run: no new device array is proposed here
  beyond `plan%sSA`, `plan%sSB` and `plan%sW`, which `binkernel_fock_bins`
  already maps and unmaps the same way.
- The XC and DF numbers in section 8 are one run each, on one molecule, for
  ranking the three candidates. They are not a design study of either.
- Every nsys figure is a **single run**. The box is shared;
  `nvidia-smi --query-compute-apps` showed no other process and 1 MiB in use
  on all four cards before and during the measurements, and every profile was
  taken after the build that produced its binary had finished. Differences
  below a few percent (the gradient's +-1% between `off` and `auto`, for
  instance) are within run-to-run noise -- which is the point being made
  there, not a claim about the sign.
- Reproduce the counts and the exactness:
  `spike_symgrad <xyz> <basis.json> <sph|cart> <plan|count|verify|verify2> [thresh]`
  from `test/`, on `sym/s16-grad-spike`. `plan` is counts only, `count` adds
  the screen pass, `verify2` is the exact oracle of section 4 (`verify` is the
  same comparison with the production whole-bin stride, i.e. the wrong
  answer, kept because it is the evidence for section 2).

## 11. The implementation, and what the design got wrong

Built on `sym/s16-grad-impl`, based on this document's own branch.
`eri%gradient` takes `skeleton`, with `eri%fock`'s semantics;
`trc_scf_gradient`'s `sym_geom` turns it on and projects the four-centre
term inside `two_electron_4c`.

### 11.1 What was built

| Decision | What landed |
|---|---|
| 1, scheme (c) | `bins_plan_emit` takes `canon`: it writes the two sides so that `(sLA, sLB) >= (sLC, sLD)`, whichever side the walk cut. `skeleton_walk` threads it through and the gradient's own plan asks for it (`trc_grad.F90`). The cheaper-side cut is untouched, so the whole reduction is kept, and `SKEL_SWAP_LMAX` no longer applies to a canonical plan. |
| 2, the signature | `sSA`, `sSB`, `sW` on every gradient kernel's compact pass and on `rysg_compact_dispatch`. `_untile` is **gone**: the gradient kernels now carry what the Fock ones do, and the splice *asserts* `(k1 - 1)*sSA(seg)` instead of rewriting it to `sA(seg)`. The orbit weight is one line next to the existing degeneracy factors (`wq = wq*sW(seg)`). The screen pass is unchanged, as section 2 predicted. |
| 3, default DFT, (i) | `two_electron_4c` projects the four-centre term; `trc_scf_gradient` gained `sym_project` for the end-of-routine projection, which `trc_capi` sets false exactly where it used to withhold `sym_geom`. The group is now handed over for every symmetric run. |
| 4, `bins_plan_build` too | `bins_plan_build` is one `plain_walk` through `bins_plan_emit` with `canon`. It is a **strict no-op** on an unsplit plan -- a bin's index *is* its `bin_key`, which packs `(la, lb)` most significantly, and `live` is ascending -- so no golden moves. It changes the plan only where `split_large_bins`/`arrange_orbits` appended chunks past the original indices, which is the hole it closes. |
| 6, the counter | `eri%nwork` and `eri%nkept` after `eri%gradient`. `nkept` is free there: the ordered-survivor compaction already copies its screen pass's atomic counter to the host once per batch. |
| 7 item 7, range separation | Both passes take `skeleton`; `two_electron_4c` adds them and projects once. `check_sym_grad`'s `compare_two_pass` pins it. |
| 5 | Untouched, as scoped. The **Fock** skeleton plan does not ask for canonical emission (`canonical` defaults false), so its 42 bra<ket energy kernels stay live and its accumulation order is unchanged. |

### 11.2 Counts and exactness, as built

Host build, gfortran 15.1.0, `TERCO_LMAX=2`, `tau = 1e-10`, the engine built
with the detected group, the density made invariant by the production
symmetrizer (`trc_sym_scf_t%symmetrize`) in the caller's space and folded to
the Cartesian working basis. `check_sym_grad`'s own printout.

Work items the gradient launches, full -> skeleton:

| system | group | basis | full | skeleton | reduction |
|---|---|---|---|---|---|
| HOCl (planar) | Cs | 6-31G* | 9 142 | 9 355 | **0.98x** (see 11.3.5) |
| CHFCl2 | Cs | 6-31G* | 96 116 | 51 717 | 1.86x |
| H2O2 | C2 | cc-pVDZ | 85 784 | 44 874 | 1.91x |
| water | C2v | cc-pVDZ | 14 627 | 9 968 | **1.467x** |
| trans-N2H2 | C2h | 6-31G* | 9 312 | 4 875 | 1.91x |
| N2H2F2 | Ci | 6-31G* | 76 357 | 39 209 | 1.95x |
| twisted C2H4 | D2 | 6-31G* | 22 148 | 7 370 | 3.01x |
| B2H6 (nonplanar) | D2h | 6-31G* | 45 108 | 12 690 | 3.56x |
| SF6 | D2h | cc-pVDZ | 2 883 269 | 1 192 250 | **2.418x** |
| benzene (D6h as D2h) | D2h | 6-31G* | 649 008 | 180 146 | 3.60x |

The two systems section 3.1 measured at the same basis reproduce **exactly**:
water 14 627 -> 9 968 (1.467x) and SF6 2 883 269 -> 1 192 250 (2.418x). So do
the segment counts of section 3.3: water Cartesian, 909 full segments and
1 320 skeleton segments, of which 575 carry their sub-range on the ket side
-- which is section 3.3's "swapped" column, now reordered by the
canonicalisation instead of changing the class. Survivors are 11 541 -> 6 609
(water) and 2 770 356 -> 1 133 419 (SF6) against the prototype's 11 534 /
6 605 and 2 761 476 / 1 130 616; the few extra are the skeleton screen's
symmetrized `dsh`, which is a looser bound by construction.

Exactness: `max|sym(grad_skel) - grad_full|` relative to
`max(1, max|grad_full|)`, worst over {J+K, J, K, long range} x {RHF, UHF} x
{random, converged} density, and the same quantity for the **raw**
(unprojected) skeleton gradient:

| system | group | Cartesian | spherical | raw, relative |
|---|---|---|---|---|
| HOCl | Cs | 7.5e-15 | 2.3e-14 | 7.5e-15 (see 11.3.5) |
| CHFCl2 | Cs | 1.0e-12 | 1.9e-12 | 1.34 |
| H2O2 | C2 | 1.6e-14 | 1.6e-14 | 0.86 / 1.08 |
| water | C2v | 2.6e-15 | 2.7e-15 | 1.39 |
| trans-N2H2 | C2h | 7.2e-13 | 1.9e-12 | 0.58 |
| N2H2F2 | Ci | 1.8e-14 | 4.6e-14 | 0.61 / 0.68 |
| twisted C2H4 | D2 | 1.5e-12 | 1.5e-12 | 2.11 |
| B2H6 | D2h | 5.3e-14 | 6.6e-13 | 1.26 |
| SF6 | D2h | 3.8e-14 | 3.2e-13 | 1.56 |
| benzene | D2h | 3.1e-13 | 2.1e-12 | 1.88 |

Two things to read off. The **spherical** column is as good as the Cartesian
one, which settles section 10's open caveat: the spike's hand-rolled density
projection was not invariant in the spherical AO space, the production
symmetrizer is, and the identity holds there to the same precision. And the
**raw** skeleton gradient is wrong by 0.58 to 2.11 times the largest
component on every molecule whose orbits are not all singletons, so the
projection is doing all the work, exactly as section 4 said.

Translational invariance of the raw skeleton gradient: 1e-16 to 9e-12
absolute, on gradients up to 460 (the random densities are not physical).

The f-basis case (`TERCO_LMAX=3`, the uncontracted l = 0..3 set,
`TRC_BIN_CAP=40` so that `split_large_bins` fires -- the combination no
gradient test had ever run): water C2v, 82 621 -> 47 991 work items, the
projected skeleton gradient 6.2e-11 against `max|grad| = 1079`, i.e. 5.8e-14
relative. Of the 1 114 skeleton segments, 0 are non-canonical and 0 carry a
class the build has no gradient kernel for; 391 carry their sub-range on the
ket. With `canon` turned off: **391 non-canonical, of which 275 carry a class
with no gradient kernel**, so `grad_classes` would `error stop`. That is the
latent hole of section 5, fired and closed.

End to end, through `trc_scf_gradient` with `sym_geom` (the skeleton
four-centre term, its own projection, and the one-electron, overlap and
nuclear-repulsion terms on top) against the same molecule's C1 gradient:
3.8e-15 (water/6-31G*), 3.1e-15 (water/cc-pVDZ spherical), 9.5e-15
(B2H6/6-31G*), with the equivariance residual and every symmetry-forbidden
component exactly 0. Against a central difference of terco's own reconverged
SCF energy (water/6-31G*, h = 1e-3, 18 displaced C1 SCFs): worst absolute
1.4e-7, inside `check_full_grad_fd`'s own 5e-7 truncation floor.

### 11.3 What did not survive contact with the code

1. **The step speedups are real but smaller than projected, because the
   gradient's kernel time does not scale with its quartet count.** Measured
   on the same V100-DGXS-32GB (cc70, nvfortran 26.5, `TERCO_LMAX=2`, commit
   `c07a4c1`, cc-pVDZ, `--grad --econv 1e-10`, nsys 2025.1.3, GPU 0, no other
   compute process and 1 MiB in use on all four cards), two independent runs
   agreeing to within 2%:

   | system | `rysg*` off | `rysg*` skeleton | kernel reduction | quartet reduction | total off | total auto | step vs off |
   |---|---|---|---|---|---|---|---|
   | benzene | 25.71 ms | 16.47 ms | **1.56x** | 3.60x | 158.7 ms | 124.3 ms | 1.28x |
   | naphthalene_d2h | 71.46 | 27.53 | **2.60x** | 3.79x | 407.5 | 228.9 | 1.78x |
   | pentacene_d2h | 426.66 | 128.79 | **3.31x** | 3.85x | 2519.9 | 966.5 | 2.61x |
   | tops_sf6 | 20.17 | 14.64 | **1.38x** | 2.42x | 184.9 | 199.0 | 0.93x |

   Section 8 projected the gradient's time by its quartet reduction, on the
   argument that "the cost per survivor does not change". It does: a skeleton
   segment is a *sub-range* of its bin, so a class run is made of more and
   shorter segments (water: 1 320 skeleton segments against 909 full ones),
   the screen pass's RUNLEN=32 chunks share a ket less often, and the compact
   pass's survivor list is less ket-uniform. The shortfall is worst where the
   segments are smallest: 1.56x of a 3.60x quartet cut on benzene, against
   3.31x of 3.85x on pentacene.

   Isolating the skeleton's own contribution to the step -- the same `auto`
   run with the gradient's `off` time substituted back -- gives **1.07x
   (benzene), 1.19x (naphthalene), 1.31x (pentacene), 1.03x (SF6)** against
   section 8's projected 1.15x / 1.24x / 1.34x / 1.06x. Those substituted
   totals (272.8 ms naphthalene, 1264.4 ms pentacene, 204.5 ms SF6) reproduce
   section 8's own `auto` totals (276.0, 1285.6, 204.8) to within 1.6%, which
   is what makes the construction trustworthy; benzene's 133.5 against 134.0
   likewise. Against `symmetry=off` the whole step is 1.28x / 1.78x / 2.61x /
   0.93x, against 1.35x / 1.78x / 2.65x / 0.99x projected.

   The **launch count does not grow**, as scheme (c) promised: 42 gradient
   launches with symmetry off and 42 with the skeleton on, on every system but
   pentacene, where it *falls* from 46 to 42 (fewer `GRAD_COMPACT_CAP`
   batches). SF6 still loses overall, as it did before this work and for the
   same reason (S15's small-molecule Fock launch problem: 1 301 -> 1 811 Fock
   launches); its gradient term is 1.38x faster.
2. **The compile cost really is a re-run of what the build already pays, but
   the design's own figure for the old set does not reproduce.** Compiling
   all 81 `trc_pc_rysg*.F90` with their own compile commands on this host
   (NVHPC 26.5, `-stdpar=gpu -acc -gpu=cc70,mem:separate -fast -O3`), old
   against new, back to back: summed single-file time **1858.0 s new against
   1939.1 s old** (0.96x), mean 22.9 s against 23.9 s. The two files that set
   the critical path, recompiled four at a time on a quiet box, are
   **unchanged**: `rysg2220` 248.9 s new against 249.1 s old, `rysg2121`
   186.1 s new against 184.0 s old. At 20 concurrent compiles the *wall*
   figures were 395.8 s new against 302.5 s old, but that difference is
   scheduling under 20-way contention, not a property of the change -- the
   same two files measured 349.6 s and 268.9 s there, against 248.9 s and
   186.1 s at four. Section 7's 277.9 s / 1473.4 s for the old set does not
   reproduce on this box today (302.5 s / 1939.1 s for the same files), so the
   A/B above, not that figure, is the comparison to use.
3. **Oracle item 9, exchange tiles on the gradient, is still not available**,
   and `sSA`/`sSB` were not the only thing in the way. `bin_tiles_plan` is
   reached only from an **exchange-only** Fock launch
   (`tiled = ... .and. launch%jfac == 0.0_dp`), and the reach bound it prunes
   on, `bin_tiles_reach(tiles, dsh)`, bounds the **exchange** blocks alone.
   The gradient's Gamma contracts all six blocks, the two Coulomb ones
   included, so an exchange-tile prune would drop real contributions. The
   kernels can now express a tile on either side, so the mechanism is
   unblocked; using it needs a gradient-specific reach bound, which no
   decision in section 7 covered.
4. **Oracle item 6 cannot be done the way it is written.**
   "`check_full_grad_fd`'s and `check_rys_grad_fd`'s finite differences with
   the skeleton on" has nothing to turn on: both run in C1, with no
   `sym_geom` and no group, so there is no skeleton path through them. The
   finite difference is instead `check_sym_grad`'s own `test_fd`, which takes
   the symmetric molecule's skeleton gradient and differences the reconverged
   energy at **displaced, C1** geometries -- a displacement along one
   Cartesian component breaks the symmetry, so the difference cannot be taken
   inside the group.
5. **Oracle item 5 needed no test edits at all, and must not have had any.**
   `trc_scf_gradient`'s `sym_geom` now implies the skeleton build, so the
   eight existing `check_sym_*` gradient cases (`_c2v`, `_c2`, `_c2h`, `_ci`,
   `_d2`, `_d2h`, `_tops`, `_scf_cs`) exercise it, with their
   forbidden-component, equivariance and C1-agreement assertions, without a
   line changed. That is as well: all eight of those files are on a
   `test/pyscf_refs/*.ref` `# files:` line, so even a comment in them would
   have made eight anchors stale.
6. **The canonicalisation cannot be global.** Decision 4 asks for it on
   `bins_plan_build` and decision 1 on the gradient's skeleton plan, but
   decision 5 keeps the energy side's 42 bra<ket Rys kernels -- and those
   exist *only* because the skeleton Fock plan swaps. Canonicalising
   `bins_plan_build_skeleton` unconditionally would have made them dead
   kernels, which is decision 5's own territory. So `canon` is a per-plan
   flag, default off, and `bins_skel_cut_first`'s `SKEL_SWAP_LMAX` guard is
   lifted only for a canonical plan, where the cut side no longer decides
   which kernel runs.
7. **The skeleton is not always smaller, and nothing stops it.** On a
   molecule where every atom lies on every mirror -- a **planar Cs**
   molecule, whose only operation besides the identity is the molecular plane
   -- every pair orbit is a singleton. The restriction then drops nothing,
   while a diagonal bin pair still becomes the rectangle `reps(A) x A` at
   weight 1/2 instead of the triangle, so the plan *grows*: HOCl/6-31G*,
   9 355 work items against the full plan's 9 142 (+2.3%), with the raw
   skeleton gradient already equal to the full one to 7.5e-15 relative. The
   skeleton Fock build has exactly the same property, and *it* has a size
   rule (`sym_skeleton_size_rule`, ABI 1.25) — so the gradient takes that
   rule's verdict rather than a second one of its own. See item 9.
   `check_sym_grad`'s `test_cs_planar` pins the growth itself, so that it is
   a documented property rather than a surprise, and it is the one
   configuration where the oracle's "the plan must be smaller" and "a segment
   was reordered" assertions must *not* hold. Section 3's molecule list has
   no planar Cs case, which is why this is not in it.
8. **A symmetric gradient now requires an orbit-built engine.**
   `trc_scf_gradient` used to refuse to reuse one on purpose ("that is
   S16's"). It now requires one: a symmetric gradient on an engine without
   orbit tables would silently fall back to the full plan, so the reuse test
   admits exactly the orbit-built engine and `two_electron_4c` builds its own
   **with the group** when there is no cache. That is why the eight
   per-group gradient cases, which pass no cache, get the skeleton.
9. **One skeleton decision, not two.** Items 7 and 8 together are a bug when
   the skeleton Fock build has been *declined*. The size rule gates
   `skel_on` inside `trc_scf_driver`, which then calls `scf_context%setup`
   **without** `sym_geom`, so the SCF's engine has no orbit tables
   (`orbits_built` false) — while the C API hands the gradient the group
   whenever a group is active, from `cx%sym_active`, which the rule never
   touches. So a declined molecule got a skeleton gradient anyway: item 8's
   reuse test refused the SCF's engine, a **second** engine was built with
   the group, and the plan it ran was the larger one of item 7 on exactly
   the molecules the rule exists to protect.

   `trc_scf_gradient` therefore takes `sym_skeleton` (default true, so a
   direct Fortran caller that passes `sym_geom` still gets the skeleton),
   `two_electron_4c` takes it as `skeleton` rather than deriving it from
   `present(sym_geom)`, item 8's reuse test reads it instead of
   `present(sym_geom)`, and `trc_capi.F90` sets it from
   `cx%res%sym_skeleton` — the rule's verdict, or whatever
   `TRC_SYM_SKELETON` / `--sym-skeleton` / `[scf] symmetry_skeleton` /
   `trc_set_symmetry_skeleton` forced. One verdict governs both builds, and
   `trc_set_symmetry_skeleton`'s documented meaning widens to say so
   (ABI 1.26).

   The *projection* is deliberately not gated on it: it costs one pass over
   3N numbers and it is what zeroes the components the group forbids, so a
   declined molecule takes the full C1 four-centre term and
   `trc_scf_gradient`'s end-of-routine `sym_grad_symmetrize` — which is what
   this routine did before S16.

   Two consequences worth knowing. Under **MPI** the driver sets
   `skel_on = .false.` on more than one rank, so the gradient now declines
   the skeleton there too; that retires section 10's warning about
   `grad_classes`' rank balance under a weighted plan, conservatively, by
   not running one. And `--df` likewise declines it, which the `--df`
   gradient branch never reached anyway.

   `check_sym_grad`'s `test_declined` pins it: water (declined, 3 atoms in
   2 orbits), HOCl (declined, 3 in 3, item 7's molecule) and B2H6 (accepted,
   8 in 3) must each reuse the SCF's engine and match the C1 gradient, which
   they do at 2.2e-15, 8.9e-15 and 9.5e-15 with every forbidden component
   exactly zero. With item 8's reuse test reading `present(sym_geom)` again
   it fails on water with `engine reused F`. `check_capi_sym`'s
   `test_skeleton_switch` is the ABI-level half, and its gradient tolerance
   tightens from 1e-7 to 1e-9 now that the two sides are genuinely different
   builds rather than the same full gradient twice.
