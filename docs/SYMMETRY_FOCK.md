# Skeleton Fock build: design (S13)

Status: **design approved (section 6a) and implemented.** Scheme 2a is built
in `src/eri/trc_bins.F90`, `src/eri/trc_bin_tiles.F90`,
`src/eri/trc_binkernel.F90` and `src/symmetry/trc_sym_pairs.F90`, and pinned
by `check_sym_fock`, `check_sym_fusion` and `check_sym_skeleton_rule`. The
measurements in sections 1-7 come from a CPU prototype (`spike_symfock`,
fork-only history, commit 496bef9; not in this tree).

The Fock build is still full C1 under symmetry. This document designs how to
form G(D) from symmetry-unique integrals only and compares the candidate
schemes on measured quartet counts.

## 1. The common math

Notation. G is the abelian point group of order |G| with real orthogonal
representation R(g) in the AO space; for every operation of
`trc_sym_group` R(g) is a **signed permutation**: `ao_map(a,g)`,
`ao_sign(a,g)` of `trc_sym_ao`. A shell quartet q = (ab|cd) maps to g q, and
for the integrals

    (g a g b | g c g d) = s_a s_b s_c s_d (a b | c d)       (signs of the AOs)

Let Dig_q(D) be the contribution of the unique shell quartet q (with its
permutational degeneracy, the eight orderings) to G (J and K separately, or
J - K/2). It is linear in D and, because the AO signed permutation is an
orthogonal transformation that maps the integrals onto themselves,

    Dig_{gq}(D) = R Dig_q(Rᵀ D R) Rᵀ.

**Condition.** D is invariant: R D Rᵀ = D for all g (and symmetric, as the
folded digestion already assumes). Then Dig_{gq}(D) = R Dig_q(D) Rᵀ.

For any set S of quartets that is closed under G, with a representative set
Q0 and weights w_q = |orbit(q)|:

    G(D) = Σ_{q∈S} Dig_q = Σ_{q∈Q0} Σ_{g∈G/Stab(q)} R Dig_q Rᵀ
         = (1/|G|) Σ_g R [ Σ_{q∈Q0} w_q Dig_q ] Rᵀ  =  (1/|G|) Σ_g R G_skel Rᵀ

(the sum over G/Stab is (|orbit|/|G|) times the sum over G, because
R Dig_q Rᵀ is invariant under Stab(q) when q is). So **G_skel = Σ_{q∈Q0}
w_q Dig_q** and the symmetrizer is the one the code already has
(`sym_scf_symmetrize`, `trc_sym_scf`). The same derivation works if the sum
over q is over *ordered* pair-of-pairs, or over a coarser object (an atom
quartet, section 4.3) as long as the object's whole set of quartets is
evaluated on the representative.

**Why every ingredient of the existing build stays exact.**

- *Schwarz screening.* Q_ab is invariant, so a quartet and its images pass
  or fail together. The surviving set is G-closed, and G_skel built from
  it symmetrizes to the full screened G. One thing is not automatic: Q
  from two symmetry-equivalent pairs differs in the last bits, so a quartet
  at the threshold can pass for q and fail for gq, and bin buckets
  (`int(-log10 Q)`) can differ across an orbit. The prototype takes
  Q := max over the orbit and puts every pair in the bin of its orbit's
  representative (7 molecules: 0 pairs moved, so this costs nothing but it
  must be done by construction, not by luck). The production build must do
  the same.
- *Density-weighted screen.* dsh(I,J) = max |D| over a shell block is
  invariant when D is, so the same argument applies. The SO-blocked SCF
  produces a D that is invariant to roundoff (the symmetrizer is applied to
  it anyway); take dsh := max over the orbit as well, which costs O(nshell²).
- *Exchange tiles.* A tile pair is dropped only when every quartet in it
  fails the kernel's own screen, so which quartets are *enumerated* does not
  change the result; tiles are cut over whatever bra/ket lists the plan uses.
- *Incremental Fock.* G(D) = G(D_ref) + G(ΔD) is linear and ΔD = D - D_ref is
  invariant when both are, so the accumulated `dref_c`/`gacc` can hold
  *skeleton* matrices and the symmetrizer is applied once, after
  `scf_fock_fold_g`. The loosened increment screening and the "convergence
  only on a full build" rules are unchanged (the dropped quartet set is
  still G-closed). The prototype verifies G_skel(ΔD) symmetrized equals G(ΔD).
- *Idempotence.* R G Rᵀ = G for a G already invariant. So an engine left
  unrestricted (full) can be added to a skeleton one before the single
  symmetrization: mixed builds are exact. This is what makes a staged
  rollout possible (section 6).
- *Not covered.* Non-invariant trial densities (CPHF/CPKS) — those keep the
  existing `fock_nosym` path. Skeleton G requires an invariant D. UHF spin
  densities are each invariant (abelian groups, real 1D irreps).

**Where to symmetrize.** Everything commutes with the symmetrizer:
decontraction (the expand/fold maps use identical coefficients on equivalent
atoms), the Cartesian→spherical transform (the `sph_parities` check in
`trc_sym_ao` is exactly T R_cart = R_sph T). The prototype confirms it:
symmetrizing in the working Cartesian basis before the fold and in the
caller's spherical AO space after the fold give the same G (table in
section 3). The cheapest place is the caller's AO space, after
`scf_fock_fold_g`, where `scf_state%w1` and `sym%symmetrize` already exist.

## 2. What was built and measured

The prototype (`spike_symfock`, fork-only history, commit 496bef9; not in this tree) builds terco's own bins (`eri%bins`) for a
molecule, assigns every shell pair the bin of its orbit's representative,
and enumerates the quartets of five schemes with the production
bin-pair admission test (`bin_s(A)+bin_s(B) <= -log10 τ`). For the full
enumeration it reproduces `eri%nwork` exactly (14627, 360415, 20671728,
4741001 quartets for water, B2H6, naphthalene, benzene: checked against the
production count in each run). It counts

- **admitted**: quartets in admitted bin pairs (what a launch evaluates);
- **Schwarz**: after the exact Q_ab Q_cd > τ;
- **kept J+K**, **kept K**: after the density screen, with the J blocks
  counted or not (the exchange-only launch);
- the **weight sum** Σ w·(1, or 1/2 for a diagonal pair), which must be the
  same for every scheme and equal the full build's: all five schemes agree in every run
  (an integral-free check of the orbit bookkeeping);

and, for the verification runs, evaluates every enumerated quartet with
`hgp_batch` (the production HGP routine, Cartesian working basis after
decontraction), digests all eight orderings into J and K, symmetrizes, and
compares with the **full** enumeration under identical screening. Settings:
HF density from a converged terco SCF (symmetrized), cc-pVDZ, τ = 1e-10,
bin resolution 1 per decade, spherical and Cartesian AOs.

Prototype scheme names: `full` (production), `s1`, `s2a`, `s2`, `s3`
(sections 4.1-4.3 below).

### 2.1 Correctness: symmetrized G_skel against the full G

Max |difference| after symmetrization (elements of J up to 24, K up to 2.4
in naphthalene), tolerance asked for: 1e-12.

| system | AO | s1 | s2a | s2 | s3 |
|---|---|---|---|---|---|
| water C2v, SCF D | sph | 3.3e-14 | 2.6e-14 | 2.8e-14 | 2.0e-14 |
| water C2v, SCF D | cart | 3.3e-14 | 2.6e-14 | 2.8e-14 | 2.0e-14 |
| water, incremental ΔD | sph/cart | 1.9e-15 | 1.3e-15 | 1.6e-15 | 1.3e-15 |
| B2H6 D2h (nonplanar), SCF D | sph | 2.0e-13 | 1.5e-13 | 1.4e-13 | 1.2e-13 |
| B2H6, incremental ΔD | sph | 7e-15 | 6e-15 | 7e-15 | 6e-15 |
| benzene D2h (D6h via S12) | sph | 5.1e-13 | 3.6e-13 | 4.6e-13 | 4.4e-13 |
| benzene D2h | cart | 6.4e-13 | 4.4e-13 | 5.0e-13 | 5.9e-13 |
| naphthalene D2h | sph | 1.1e-12 | 1.5e-12 | 1.2e-12 | 1.1e-12 |
| naphthalene D2h | cart | 1.5e-12 | 2.0e-12 | 2.1e-12 | 2.1e-12 |

For water and B2H6 the tests were also run with J-only and K-only screens
and K-only with exchange tiles (J and K compared separately, errors <= 2.0e-13 (J) and 1.2e-14 (K) for the SCF density, <= 1e-15 for ΔD). The
raw, unsymmetrized skeleton differs from the full G by 3-50 (the
symmetrization is doing real work). The naphthalene errors are 5e-14
*relative* to the largest element, i.e. the summation noise of 2e7
accumulated quartets in double precision: they sit above the absolute 1e-12
only because the system is larger. S14's oracle should therefore state the
tolerance relative to max|G| as `check_nosym` does (1e-11·max(1,|G|) there).
The prototype's full G also agrees with `eri%fock` to 2.6e-13 (water),
1.5e-11 (B2H6), 2e-10 (benzene), 1.4e-9 (naphthalene): the differences are
the production kernel's different screening/primitive cut and are not part
of the symmetry question.

### 2.2 Quartet counts (cc-pVDZ, τ = 1e-10, density screening with the SCF density)

"ideal" is the number of *distinct shell permutations* the group
induces (printed by the prototype as the effective order). It is smaller
than |G| whenever an operation acts on every atom trivially: **in a planar
molecule the plane's mirror fixes every atom**, so D2h naphthalene, benzene
and pyrene have effective order 4, not 8; water (all atoms on the C2 axis
or swapped H) has 2 of 4. Only nonplanar B2H6 reaches 8 of 8. At the
*shell* level the mirror acts as a pure sign, so its gain would be an
AO-component zero-skipping optimisation (the odd-parity integrals vanish
exactly) and not a Fock-build one (section 7).

Quartets launched (admitted), and the reduction vs `full`:

| system (nop / ideal) | full | s1 (reps × all, ordered) | s2a (bin-ordered reps) | s2 (per-quartet orbit rep) | s3 (atom quartet) |
|---|---|---|---|---|---|
| water C2v (4 / 2) | 14 627 | 21 413 (0.68×) | 9 968 (1.47×) | 8 978 (1.63×) | 9 875 (1.48×) |
| B2H6 D2h (8 / 8) | 360 415 | 239 341 (1.51×) | 107 763 (3.34×) | 85 037 (4.24×) | 99 775 (3.61×) |
| naphthalene D2h (8 / 4) | 20 671 728 | 10 913 456 (1.89×) | 5 451 483 (3.79×) | 5 180 294 (3.99×) | 5 380 239 (3.84×) |
| benzene D2h (8 / 4) | 4 741 001 | 2 705 946 (1.75×) | 1 342 393 (3.53×) | 1 205 476 (3.93×) | 1 308 534 (3.62×) |
| pyrene D2h (8 / 4) | 76 742 559 | 41 382 127 (1.85×) | 20 718 357 (3.70×) | 19 287 735 (3.98×) | 19 702 169 (3.90×) |
| pentacene D2h (8 / 4) | 150 479 149 | 77 526 204 (1.94×) | 39 110 378 (3.85×) | 37 645 958 (4.00×) | 38 182 652 (3.94×) |

How much survives the triangular fold and screening. Reading naphthalene
(ideal 4): the ordered scheme s1 gets only 1.89× = ideal/2 (it gives up
the ab↔cd fold, the plan's |G|/2); keeping the fold at bin-pair level (s2a)
recovers 3.79×; the per-quartet orbit test (s2) reaches 3.99×, i.e. all of the
ideal. The remainder of s2a is the **diagonal bin pairs**, which s2a
treats as `reps × all` (a rectangle, gain ideal/2 on them only).

Screening does not change the picture: the survivor fraction is the same in
every scheme (symmetrically screened sets are G-closed), so the reductions
hold at every stage. Naphthalene, reduction vs full:

| stage | full count | s1 | s2a | s2 | s3 |
|---|---|---|---|---|---|
| admitted | 20 671 728 | 1.89× | 3.79× | 3.99× | 3.84× |
| + exact Schwarz | 20 255 420 | 1.89× | 3.79× | 3.99× | 3.84× |
| + density, J+K | 19 066 438 | 1.89× | 3.79× | 3.99× | 3.84× |
| + density, K only | 18 205 345 | 1.89× | 3.78× | 3.99× | 3.84× |

Survivor fractions (J+K) at τ = 1e-10: 92% (naphthalene), 95% (benzene), 90%
(pyrene), 88% (pentacene), 98% (B2H6), 100% (water): these compact cc-pVDZ
systems screen little, so the density-screen interaction is mostly a
non-event here; it matters for gly-class systems in S15. Tightening and loosening τ (naphthalene: 1e-8, 1e-10, 1e-12):
s2a reductions 3.77×, 3.79×, 3.81× and s2 3.99× throughout. The
**exchange tiles** (64-pair tiles, bound Q_A Q_B max dsh over the tiles'
shells) remove 0.02% of the admitted quartets on naphthalene
(20 671 728 → 20 667 313, full), 0.09% on pentacene (150 479 149 →
150 346 769) and the same fraction in every skeleton: no
interaction, but also no test of the tile code on a system where tiles
matter. That is an S14 oracle item (an extended D2h system with a density
that is small between its ends).

Launch granularity (naphthalene): bin pairs 5 186 (the same for full, s2a,
s2, s3); s1 10 291. A skeleton plan has one weight per bra pair
(`osz`, the orbit size); splitting the rep list of a bin by weight gives
6 553 weighted segments for s2a (+26%), 5 186 for s2, but s3 needs bins split
by atom pair: 3 513 bin×atom-pair groups instead of 118 bins (30×, 6 845 vs
121 for pyrene).

## 3. Mapping onto terco's architecture

What a "quartet" is here: no list exists. A launch walks admitted *bin
pairs*; each descriptor (segment) is (sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB,
sD, class) and the kernel decodes (segment, ab, cd) from the flat thread
index; `sD` is the triangle of a diagonal bin pair, otherwise a rectangle of
sA × sNB. Rank split is by work item (`rank`, `nranks` stride). Exchange tiles
already give a segment that covers a sub-range of each bin
(offset `sOA`/`sOB`, extent `sA`/`sNB`, bin stride `sSA`/`sSB`). Weights
today are per quartet and fixed: the factors 1/2 for i=j, k=l and ab=cd.
Scalars `jfac`, `kfac` are per launch and warp-uniform. The generated
kernels (`pc_dispatch`, `rys_dispatch`) take all these as arguments, so a
change to the *plan* needs no change to a kernel as long as it can be
expressed with those arguments.

## 4. The candidates

### 4.1 Scheme 1: unique bra pairs, ordered (the plan's candidate 1)

*Math.* Σ over ordered (p,q), p ∈ reps, q ∈ all pairs, weight |orbit(p)|;
each ordered pair carries half of the folded contribution:
G = ½ Σ_{ordered (p,q)} R(p,q), with R the eight-ordering digestion without
the same-pair factor (this makes p = q come out right with the same ½).
Exact for invariant D. Measured: 1.9× on planar D2h (ideal/2), and
**worse than full on water (0.68×)**, where most pairs have a nontrivial
stabilizer.

*Changes.* `trc_bins`: orbit tables, closed bins, reps first within a bin
(ordered by orbit size); `bins_plan_build`: for each *ordered* admitted bin
pair a rectangular segment reps(A) × B (so both (A,B) and (B,A) for A ≠ B)
with sD = false; launch scale ½·osz. `trc_binkernel`: launches grouped by
(class, weight), scale folded into `jfac`/`kfac`. Exchange tiles: cut the rep
sub-range. **No generator change**: rectangle segments, sub-range offsets and
launch scalars all exist. MD J: see 4.5. Compile-time cost: none (host
files only).

*Risks.* Poor for small molecules; leaves a factor 2 on the table for no
simplification over 4.2 (4.2 needs only the choice of which side of an
unordered bin pair to cut).

### 4.2 Scheme 2a: unique bra pairs with the fold kept at bin-pair level (recommended)

This is not in the plan's list; it falls out of the structure. Because the
bins are closed under G, an orbit never straddles two bins, and for an
unordered pair of *different* bins {A,B} the sum Σ_{p∈A, q∈B} f(p,q) can be
cut on **either** side:

    Σ_{p∈A,q∈B} f = Σ_{p∈reps(A)} |orbit(p)| Σ_{q∈B} f(p,q)     (or the mirror)

with no further fold left to break (A and B are different bins, the
rectangle appears once). Cut the side with fewer quartets. A *diagonal* bin
pair (A,A) is handled as the ordered rectangle reps(A) × A with scale ½, the
only place the fold is lost (naphthalene: this is the whole gap between 3.79×
and 3.99×; on a large molecule with many bins the diagonal pairs are a
shrinking fraction). Exactness: `wsum` identity, G_skel verified in 2.1.

*Weights.* w = osz(p) (× ½ on a diagonal bin pair): one value per bra pair.
Reps within a bin are stored contiguously, sorted by orbit size, so each
value is a contiguous sub-range = a tile-shaped segment; the weight is a
property of the segment.

*Files and routines that change* (all host except where noted):
- `src/eri/trc_bins.F90`: `pair_bins_t` gets `rep_cnt(bin)`, `rep_w(...)`
  (sub-range ends per orbit size); `bins_build_binned_pairs` makes the bins
  orbit-closed (the key uses the representative's size bucket; Q symmetrized
  over the orbit) and orders each bin reps-first, then by orbit size;
  `split_large_bins` cuts only at orbit boundaries (an orbit has at most
  |G| pairs, so a chunk can overshoot BIN_CAP by < |G|, and the segment size
  guard `nt > 2**31` keeps its margin); `bins_plan_build` produces the
  rep segments (choice of side, one segment per weight), `bins_plan_sort_by_class`
  adds the weight to the sort key.
- `src/eri/trc_binkernel.F90`: `launch_classes` passes the plan's `sW` to
  the kernels, one launch per `(class, mask)`, and leaves `jfac`/`kfac`
  unweighted. (As built for S14 it split a class run at weight changes and
  multiplied `jfac`, `kfac` by the weight for that launch; W2 (section 4.2b)
  replaced that -- see the compile-time note below.) `jfac`/`kfac` stay unweighted because the
  kernels test `jfac /= 0`, `kfac == 0` for warp-uniform branches, so a
  weight must not reach them. `count_kept` and `fock_all_nosym` are
  untouched (the nosym path stays full).
- `src/eri/trc_bin_tiles.F90`: `cut_bin` over the rep sub-range; `bin_tiles_plan`
  emits the weighted tile segments.
- `src/eri/trc_eri.F90` (`build_structures`, the `fock*` entries): take an
  orbit table (pair images, rep flag, orbit size) built from `trc_sym_ao`
  (a new module `trc_sym_pairs`, one more host table) and a flag; return the
  skeleton G. `mdg_pass`/`add_md_j` see 4.5.
- `src/scf/trc_scf_fock.F90`: after `scf_fock_fold_g`, `sym%symmetrize` on G
  (per spin) when the skeleton flag is on; nothing else. `trc_scf_incremental`
  is untouched.
- Generators: **none**. MPI: **none**.

*Compile-time cost.* As designed, no generated source changes, so the
per-class kernel build was unchanged (0 min of the 15-20) and the edited
host files are small modules (seconds). The per-segment weight array
`sW(seg)` was listed here as the alternative -- every `pc*`/`rys*` kernel
and dispatcher signature changes, a full regeneration and one full rebuild
of the kernels, with no growth in code (one load and multiply per quartet)
-- and judged unnecessary. **S15 and the W-stages made it necessary**: the
launch scalar costs one launch per `(class, mask, weight)`, which S15
measured as the reason the quartet savings do not convert into wall time
below ~20 atoms, and W1 showed both the extra kernel time and nvfortran's
fixed per-offload-region host cost scale with the launch count. The array
is what is built; see section 4.2b.

*Risks.* (i) More launches: `(class, weight)` instead of `class`; the
weights present in a D2h planar molecule are a few (1, 2, 4), most pairs are
at the maximum. Segments +26% on naphthalene, +67% on benzene, +63% on
pyrene, +41% on pentacene (weighted segments 6 553 / 5 105 / 8 356 / 7 304 vs
5 186 / 3 058 / 5 119 / 5 163);
small segments cost warp occupancy at the tail. Measure on GPU in S15.
(ii) A rep-first order destroys the "natural" shell order inside a bin that
tiles use for locality; reps are still in shell order among themselves and
the non-reps are never a bra, so tiles over reps stay compact, but the
sort-by-weight option `TRC_BIN_SORT` must be disabled or made orbit-aware.
(iii) Closure by construction (sec. 1).

### 4.2a When scheme 2a is used: the size rule

Scheme 2a is not used on every symmetric run. A measured size rule
(`sym_run_skeleton_size_rule` in `src/symmetry/trc_sym_run.F90`, applied by
`trc_scf_run`) decides, and it decides only *that*: the group, the standard
or input frame, the SO basis, the irrep labels, the orbital symmetries and
every reported quantity are the same whichever Fock build runs, and the two
agree to the project's gate. A declined skeleton build is the S4 path -- the
full build with F symmetrized after it -- and it is reported as
`trc_scf_result_t%sym_skeleton = .false.` and in the run's own log line
(`symmetry: full Fock build (size rule: 7 atoms in 4 orbits)`).

**The rule.** Use the skeleton build when the molecule has at least two atoms
per symmetry-inequivalent atom:

    natm >= 2 * norbit,    norbit = atom orbits of the group

`norbit` is `sym_geom_atom_orbits` of `trc_sym_geom`, counted from the atom
permutation `geom%perm` that detection has already built -- once per geometry,
on tables that exist either way. Its loop is O(natm^2 |G|), not O(natm |G|):
each of the `natm * |G|` merges relabels the whole array. That is ~10^4
integer compares at terco's largest symmetric molecule (36 atoms, |G| = 8) and
so is free in practice, but it is not linear, and a caller that wanted the
orbit count per SCF iteration rather than per geometry should be told the
right exponent.

**Provenance.** The four ERI-time numbers below were **not measured by the
change that added this rule**, which ran on a host build with no GPU. They are
stage **W2's** numbers, from section 4.2b's "W2
results" table, re-read here as the ratios the rule is fitted to; the rule
itself has not been re-measured on a GPU (item 5 of the open runs below).

System and settings, as W2 ran them: SF6, benzene, naphthalene and pentacene
in **cc-pVDZ**; integral screening threshold **`--thresh 1e-12`** (terco's
default, not overridden -- the `tau = 1e-10` on section 2.2's quartet table is
that pre-pass's own screen, not this one); convergence gate **`--econv
1e-10`**; the default **SAD** atomic guess; `--symmetry auto` against the same
build's `--symmetry off` as the baseline; Tesla V100-DGXS-32GB (cc70),
nvfortran 26.5, `-stdpar=gpu -acc -gpu=cc70,mem:separate`, `TERCO_LMAX=3`; ERI
time per steady-state Fock build, taken from nsys as the difference between
two consecutive `--iterations` runs. `*` marks what the rule chooses:

| system | atoms | orbits | atoms/orbit | quartets on/off | ERI on/off | rule |
|---|---|---|---|---|---|---|
| tops_sf6 | 7 | 4 | 1.75 | 2.44x | **0.44x** | full* |
| benzene | 12 | 4 | 3.00 | 3.55x | 1.41x | skeleton* |
| naphthalene_d2h | 18 | 5 | 3.60 | 3.81x | 2.07x | skeleton* |
| pentacene_d2h | 36 | 10 | 3.60 | 3.87x | 2.78x | skeleton* |

Energies are unaffected either way: |dE| between the symmetric and the C1 run
on the same symmetrized geometry was 2.9e-13 to 1.4e-12 Eh against a 1e-10
gate, with identical iteration counts (S15, unchanged by W2).

**Why orbit size and not size.** What scheme 2a saves is the orbit reduction
and what it costs is the extra launches the weights and
`bins_skel_cut_first`'s bra/ket swap create. The mean atom orbit size is a
cheap proxy for the saving, and it is available where the decision has to be
made -- before the engine is built, since the engine is built with the group
or without it and no plan exists yet. It tracks the measured scheme-2a
reduction across every system section 2.2 counted:

| system | group | atoms/orbit | s2a quartet reduction |
|---|---|---|---|
| water | C2v | 1.50 | 1.47x |
| tops_sf6 | D2h | 1.75 | 2.44x |
| b2h6_d2h | D2h | 2.67 | 3.34x |
| benzene | D2h | 3.00 | 3.53x |
| naphthalene_d2h | D2h | 3.60 | 3.79x |
| pentacene_d2h | D2h | 3.60 | 3.85x |

A molecule whose atoms sit on many symmetry elements at once has small orbits
however many atoms or functions it has -- SF6's sulfur is fixed by all eight
operations of D2h and each fluorine by four, and a planar molecule's own plane
fixes every atom, which is the same fact as section 2.2's "effective order 4,
not 8". That is also why the basis does not enter: S15 measured benzene at
0.46x (cc-pVDZ, 114 bf) and 0.48x (cc-pVTZ, 264 bf) and naphthalene at 0.98x
and 0.99x across the same 2.3x increase, while naphthalene at 412 bf was a
wash and pentacene at 378 bf a 2.14x win.

**Where the threshold comes from.** SF6 loses with a 2.44x reduction and
benzene wins with 3.53x, so break-even is a reduction in (2.44x, 3.53x), i.e.
a mean orbit size in (1.75, 3.00). The threshold is set at the **low** end of
that interval: the rule declines only what was measured to lose and keeps the
skeleton build wherever it might still pay. That asymmetry is deliberate --
post-W2 benzene at 12 atoms is a 1.41x win, so an over-eager rule would now
cost time, which is the mistake a rule fitted to the *pre*-W2 numbers would
make.

**Rejected predictors, on this same data.**
- *The atom count.* The four ERI rows allow a threshold anywhere in (7, 12],
  and every such threshold also declines 8-atom B2H6, whose 3.34x reduction
  is closer to benzene's than to SF6's.
- *The number of inequivalent atoms.* It does not separate the data at all:
  SF6 and benzene both have 4.
- *Predicted quartets per launch against the device's fill threshold* (~164k
  work items, 80 SMs at 128-thread blocks). This is the most principled
  candidate on paper, and W0's own measurement rules it out as the
  discriminator: after fusion benzene and SF6 sit at the *same* ~46% of
  quartets in launches too small to fill the card, yet benzene gains 1.41x
  and SF6 loses 2.3x. What separates them is the reduction (3.55x against
  2.44x) and the launch inflation over the C1 plan (1.57x against 3.43x), and
  both follow from the orbit sizes. It is also not available at the decision
  point.

**The small groups.** Cs, Ci, C2, C2v and C2h are declined wherever an atom
sits on a symmetry element, which is most molecules in those groups: water
reaches 1.47x, and a planar Cs molecule has no nontrivial shell-pair orbit at
all (its mirror is a pure sign operation, section 6 open question 1). This
answers section 6's open question 2 from measurement rather than by
assumption: a ~1.5x reduction cannot pay for the launches when SF6's 2.44x
does not.

**What the threshold rests on, and what would firm it up.** It rests on a
single measured loss -- SF6 -- with benzene the nearest measured win on the
other side, so the interval (1.75, 3.00) is bounded by one run at each end and
nothing inside it. These are the runs that would narrow it, in order of value
-- all on the GPU, none of them done here:

1. **B2H6 cc-pVDZ, D2h, skeleton forced on against off.** It is the only
   system inside the (1.75, 3.00) interval: 2.67 atoms per orbit, 3.34x
   reduction. If it wins the threshold belongs at 2 as set; if it loses the
   threshold belongs nearer 3 and benzene becomes the boundary case. One pair
   of runs decides it.
2. **C2H4 cc-pVDZ (3.00 atoms per orbit, 6 atoms)** against benzene (3.00,
   12 atoms). Same predictor, half the atoms: it says whether the atom count
   needs to appear in the rule at all, or whether the orbit size carries it.
3. **A large molecule with small orbits** -- SiF4 or SF6 at cc-pVTZ/cc-pVQZ,
   or a bigger octahedral cluster -- to check that size never rescues small
   orbits, and that the ratio is still basis-invariant where the rule
   declines (S15 showed basis-invariance for benzene and naphthalene only).
4. **A non-D2h pair either side of the threshold**: a nonplanar C2v or C2h
   molecule with every atom off the elements (mean orbit 2, accepted) against
   a planar one (declined). Every ERI measurement so far is D2h.
5. **A re-measurement of the four systems above with the rule in place**, to
   confirm the three accepted ones are bit-unchanged in time as well as in
   energy (the rule must be inert where it accepts).

**Forcing it.** `TRC_SYM_SKELETON=1` forces the skeleton build and `=0` forces
the full one, over any request, so both paths stay measurable from one binary
(the same shape as `TRC_MDJ` over the J engine's own size rule). From code:
`trc_sym_scf_t%skeleton = SYM_SKEL_ON / SYM_SKEL_OFF`. From the CLI:
`--sym-skeleton` / `--no-sym-skeleton`. From the input file:
`[scf] symmetry_skeleton = true/false`. From the C ABI:
`trc_set_symmetry_skeleton(h, 1/0)` (ABI 1.25). The tests are
`check_sym_skeleton_rule_{auto,on,off}`.

### 4.2b Segment fusion and launch cost (W0-W2)

S15 measured the skeleton Fock build on a Tesla V100-DGXS-32GB (cc70),
nvfortran 26.5, `-stdpar=gpu -acc -gpu=cc70,mem:separate`, `TERCO_LMAX=3`,
cc-pVDZ unless stated, `--econv 1e-10`, with `--symmetry off` of the same build
as the baseline. The quartet savings match section 2.2's predictions for
scheme 2a, but below roughly twenty atoms they did not convert into wall time
(before W2):

| system | atoms | bf | Fock/iter on | Fock/iter off | ratio | quartets on/off |
|---|---|---|---|---|---|---|
| SF6 | 7 | 102 | 0.0197 | 0.0070 | 0.36x | 2.44x |
| benzene | 12 | 114 | 0.164 (total) | 0.076 | 0.46x | 3.55x |
| benzene, cc-pVTZ | 12 | 264 | 2.377 | 1.147 | 0.48x | - |
| naphthalene | 18 | 180 | 0.0239 | 0.0234 | 0.98x | 3.81x |
| naphthalene, cc-pVTZ | 18 | 412 | 4.187 | 4.132 | 0.99x | - |
| pentacene | 36 | 378 | 0.0667 | 0.1429 | 2.14x | 3.87x |

The ratio tracks the atom count, not the basis size (naphthalene at 412 bf is
still a wash, pentacene at 378 bf wins). Energies are unaffected: |dE| between
the symmetric and the C1 run on the same symmetrized geometry was 2.9e-13 to
1.4e-12 Eh against a 1e-10 gate, with identical iteration counts.

**Where the loss was** (benzene, 88 ms over 9 iterations, from nsys):

- Not in setup: `setup` is 0.018 s with symmetry against 0.015 s without, and
  its `pairs`, `schwarz` and `bins` sub-columns agree to a millisecond.
- About a fifth is extra GPU kernel time. The ERI kernels (`trc_pc_*`) take
  90.36 ms over 1067 launches with symmetry against 69.65 ms over 495 without:
  2.16x the launches, each only 1.66x cheaper (84.7 us against 140.7 us), so
  the total rises 1.30x and the device is left underfilled.
- The rest is host-side CUDA API cost that scales with the launch count:
  `cuStreamSynchronize` +26.2 ms (+1691 calls), `cuMemsetD32Async` +10.5 ms
  (+4576 calls), `cuMemcpyHtoDAsync` +8.2 ms (+2366 calls). A `do concurrent`
  launch is synchronous, so the host's sync time overlaps the launch's own
  kernel time (15.5 us mean sync against an 84.7 us mean kernel) and the two
  must not be added. The unambiguously additive host cost is memset + HtoD,
  18.7 ms, about 21% of the 88 ms. This is derived from the published
  figures, not from a fresh profile.
- 306 of benzene's 1071 `trc_pc_*` launches belong to the SAD atomic guess and
  are identical with symmetry on and off. Per Fock build the count is 85
  against 21 (naphthalene 4.19x, pentacene 4.29x, SF6 5.43x), 4.05x for
  benzene.

At 18 atoms the same decomposition comes out the other way: naphthalene's ERI
time is 177.0 ms with symmetry against 264.8 ms without, because the
per-launch time falls 3.7x while the launch count rises 2.44x. The fixed
per-launch costs are amortised once there is enough work per launch. This is
the risk 4.2 (i) recorded: small segments cost occupancy at the tail.

**What was changed.** `binkernel_fock_bins` (`src/eri/trc_binkernel.F90`) used
to build one launch per `(class, mask, weight)`: the orbit weight was a scalar
folded into the launch's `jfac` and `kfac`, so a change of weight ended the
launch. The realised weights are multiples of one half (`skeleton_walk` and
`bin_tiles_plan` emit `half * sec_w`, `half = 0.5` on a diagonal bin pair; the
sort key is `nint(2*sW)`), measured as {1/2, 1, 2, 4} on benzene, naphthalene,
pentacene and SF6: weight 8 needs a shell pair with a trivial stabiliser, and
in all four every shell pair is fixed by at least one operation. The MD J
engine already did not split by weight (`weight_copies` in
`src/eri/trc_mdj.F90` scales the copies' primitive coefficient by the orbit
size); the bin path was the outlier.

- **W0 (bound the win, no product code).** A launch count is a count of
  distinct `(class, mask, weight)` triples, because the plan is sorted on that
  key immediately before the launch loop; a fused count is the count of
  distinct `(class, mask)` pairs. Predicted fusion (build-gpu-l3, one
  steady-state build, `--symmetry auto`):

  | system | atoms | segments | launches then | fused | ratio | off |
  |---|---|---|---|---|---|---|
  | benzene | 12 | 5712 | 85 | 33 | 2.58x | 21 |
  | naphthalene_d2h | 18 | 7871 | 88 | 36 | 2.44x | 21 |
  | pentacene_d2h | 36 | 8959 | 90 | 36 | 2.50x | 21 |
  | tops_sf6 | 7 | 4217 | 114 | 36 | 3.17x | 21 |

  No `(class, mask)` group holds more than 4 weights. Fusion lands at 33-36 and
  not at the symmetry-off 21 because `bins_skel_cut_first` swaps bra and ket for
  some bin pairs and forms class keys the full plan never does; weights do not
  cause that and fusion cannot remove it. `launch_cap()` (`TRC_LAUNCH_CAP`) is
  0 by default and nothing in the tree sets it, so it gives back none of the
  win; if it is ever enabled its floor must sit above the largest fused group
  (327k / 1.26M / 9.25M / 205k quartets for the four systems) and it should not
  be scaled down. On 80 SMs at 128-thread blocks (about 164k work items to
  saturate), 68% of benzene's and 100% of SF6's quartets sat in launches too
  small to fill the card; fusion takes both to about 46%.
- **W1 (host-side per-segment cost): no target.** The three API costs are not
  per segment; they are a fixed multiple of the launch count, nvfortran's
  prologue for one `do concurrent` offload region:

  | API | delta (calls) | per extra launch |
  |---|---|---|
  | `cuMemsetD32Async` | +4576 | 7.94 |
  | `cuMemcpyHtoDAsync` | +2366 | 4.11 |
  | `cuStreamSynchronize` | +1691 | 2.94 |

  Benzene has 64 extra launches per build, 576 over 9 iterations; constants of
  8/4/3 predict 4608/2304/1728 against the measured 4576/2366/1691 (within
  1.4%, 2.7%, 2.1%). Benzene's 5712 segments over 85 launches (about 67 per
  launch) exclude a per-segment cost, which would be 5-14x larger. The launch
  loop contains no `!$acc` directive, allocation or zeroing; all 30 array
  arguments of the Fock kernels are explicitly mapped; terco's own per-build
  device traffic is about 11 HtoD, 2 DtoH and 6 syncs, 3-4% of the deltas.
  Nothing on the host side can be batched or hoisted, so only a lower launch
  count removes the excess. Fusing mask-split launches was rejected: it would
  reorder segments and so the last bits of the serial host energy, which
  `check_sym_parity` relies on reproducing bit for bit.
- **W2 (per-segment weight in the kernels).** `plan%sW` is passed to
  `pc_dispatch` and `rys_dispatch` with the other per-segment arrays and
  multiplies the J and K contributions at digestion (hoisted to a scalar per
  segment); the `plan%sW(c1 + 1) /= wt` exit is gone, so there is one launch
  per `(class, mask)`. `jfac` and `kfac` stay scalars and unweighted, because
  the generated kernels branch warp-uniformly on `jfac /= 0` and `kfac == 0`
  and a weight must never reach those tests. 97 generated kernels change
  signature (`scripts/gen_perclass.py`, `scripts/gen_rys.py`). The weight array
  is `real(dp)`: the 1/2 of a diagonal bin pair is not an integer orbit size.

**W2 results.** cc-pVDZ, `--econv 1e-10`, `--symmetry auto`, V100 cc70,
nvfortran 26.5, `TERCO_LMAX=3`; per steady-state Fock build, from nsys as the
difference between two consecutive `--iterations` runs (W0's method), on a
build made after the parity-skip merge:

| system | launches before | after | W0 predicted | ERI ms before | after | sym off |
|---|---|---|---|---|---|---|
| benzene | 85 | 33 | 33 | 8.50 | 4.50 | 6.34 |
| naphthalene_d2h | 88 | 36 | 36 | 15.70 | 11.79 | 24.46 |
| pentacene_d2h | 90 | 36 | 36 | 71.48 | 66.61 | 184.92 |
| tops_sf6 | 218 | 72 | 36 | 14.89 | 7.29 | 3.23 |

- The per-kernel launch histograms after fusion are {1: 33}, {1: 36}, {1: 36}
  and {1: 36, 2: 18}: every system lands on exactly its number of distinct
  `(class, mask)` groups. The excess over the symmetry-off 21 is the
  `bins_skel_cut_first` swap.
- Host API calls fall in proportion to the launches. Benzene per build:
  `cuMemsetD32Async` 680 -> 264 (off: 168), `cuMemcpyHtoDAsync` 340 -> 132
  (84), `cuStreamSynchronize` 255 -> 99 (63).
- Parity skipping is engaged: benzene and SF6 each run 18 z-masked kernels
  under `--symmetry auto` at cc-pVDZ (`trc_pc_k*z*`). The masks split SF6's
  groups further, hence 72 and not W0's 36 (a 3.03x reduction against 3.17x
  predicted); W0's SF6 numbers came from a build that predates the merge. The
  other three systems' fused counts do not depend on the masks.
- Benzene's ERI time is below its symmetry-off baseline (4.50 against 6.34 ms,
  3.86 against 4.44 on the next build); naphthalene and pentacene are 2.07x and
  2.78x faster than off. SF6 is still 2.3x above its C1 baseline (down from
  4.6x): at 7 atoms with heavily stabilised orbits there is too little work per
  class, which section 4.2a's size rule handles.
- Registers: no spills introduced. `cuobjdump -res-usage` on the linked
  `libterco.so`, `TERCO_LMAX=3`, over 1225 device functions: functions at the
  255-register ceiling 587 -> 586, mean 183.6 -> 183.6, total STACK
  2 684 272 -> 2 684 608 bytes (+0.013%), max STACK unchanged at 169 600,
  LOCAL 0 throughout, and every kernel whose STACK was 0 still has STACK 0. For
  the 353 `trc_pc_*` Fock kernels, 4 gained registers (largest +16) and 19
  lost; largest STACK increase +312 bytes (`trc_pc_k1222`, +0.9%). (`LOCAL`
  is unused in this build and the frame is in STACK, which mixes declared
  scratch with spill slots; a zero `LOCAL` alone would not show an absence of
  spills.)
- Compile time is unchanged: the generated kernels grow by comment lines and
  one dummy, no executable statement (`trc_pc_k2222.F90` 20 270 -> 20 313
  lines). nvfortran, interleaved on a quiet box: (dd|dd) 562.4 -> 547.3 s,
  (dp|dp) 132.5 -> 132.8, (dd|ps) 64.6 -> 63.6, Rys (ff|ff) 11.53 -> 11.61,
  Rys (dd|dd) 4.49 -> 4.48, (pp|pp) 6.12 -> 6.11.
- Energies move in the last bits, from the fusion and not the weight. Moving
  the weight from `jfac`/`kfac` onto the digestion's `wq` is exact (both are
  powers of two; a host build with the new kernels and the weight split still
  in place reproduced the old energies bit for bit for water, B2H6 and
  benzene). Merging the launches reorders the accumulation into G, because
  the chunked Rys and warp-role kernels map threads to work items in blocks
  whose boundaries depend on the launch's extent. On the serial host build
  benzene moves 9.3e-12 Eh (water 2.8e-14, B2H6 6.5e-13) against a 1e-10
  gate, iteration counts unchanged. The symmetry-off plan has no weight to
  fuse on and stays bit-identical. On the device the before/after spread
  (<= 4.6e-13 Eh total, <= 2.4e-12 on the two-electron part) is inside the
  same binary's own run-to-run spread (<= 2.8e-13 and <= 1.25e-12).
- `eri%nwork` is unmoved; `check_sym_fock` pins it exactly (water 14627/9968,
  B2H6 360415/107763, benzene 4741001/1342393, naphthalene 20671728/5451483).
  Fusion changes how the quartets are launched, not which. `check_sym_fusion`
  pins one launch per `(class, mask)`.
- Benzene's device profile carries `trc_mdj_*` kernels, so `jfac == 0` and the
  bin call is exchange-only and tiled (`bin_tiles_plan`, `walk_skeleton`).
  Both skeleton builders fill `sW` through `bins_plan_emit`, so neither needed
  a separate change.

**Measuring what is active.** A GPU build configured before the parity-skip
merge held no masked kernels (zero `trc_pc_k*z*.F90` files against 30 after),
and a measurement on it concluded that parity skipping was inert; a second run
on the same binary "confirmed" it. Zero launches had meant the feature was
absent. Record the commit the measured binary was built from (not the
worktree's), and check that the expected generated files exist in the build;
a silent fall back to the full plan is otherwise invisible. `eri%parity_axis`,
`nwork_masked`, `nwork_variant` and `nlaunch_variant` exist on the engine but
are surfaced outside `check_sym_parity` nowhere.

**Not measured.** The tables above are the post-W2 matrix at cc-pVDZ on one
GPU; benzene was not re-run at cc-pVTZ after W2. The ideal for benzene's ERI
time is about 19.6 ms of the 90.36 ms (69.65 ms divided by the 3.55x quartet
reduction), against which the fused launch is not yet compared. SF6 reaches only 2.44x
although it is nonplanar D2h, where section 6 expects planarity to be the
limiter: its atoms sit on many symmetry elements at once (S is fixed by all
eight operations, each F by four), so the orbits are small regardless. Local
site symmetry is a second limiter beside planarity.

### 4.3 Scheme 2: unique pair-of-pairs with the fold (the plan's candidate 2)

*Math.* One quartet per orbit of the unordered {ab, cd}; weight
|G|/|Stab({ab,cd})|; G_skel = Σ w Dig. Measured ideal: 3.99× (planar), 4.24×
(B2H6, ideal 8; most gain lost to stabilizers, not to the scheme).
*Mechanism.* The orbit of a quartet is not a rectangle, so the per-quartet
condition "is {p,q} the least of its |G| images" has to be evaluated by the
kernel: |G| table lookups img(p,g), img(q,g) (a pair-image table of size
npair × |G|) and compares, then a weight in registers. The weight is
per quartet, so the *generated* kernels must change: the screen text
(`emit_screens` in `scripts/gen_perclass.py`, shared by rotaxis and
`gen_rys.py` which splices the same skeleton), the digestion's `wq`, and the
dispatch arguments (new device table), in every HGP and Rys class (the
chunked Rys kernels digest per chunk and need the weight in each chunk).
*Compile cost.* About 30 extra lines per kernel against kernels of hundreds
to thousands of lines, hence a few percent of the 15-20 min GPU compile
(**estimate**, not measured: no GPU here) plus one full regeneration and
rebuild; the real risk is register pressure in the heavy classes that are
already at 255 registers with local spill (the filter variables are dead
before the VRR, so likely fine, but not measured). The device table adds
2·|G| uncoalesced loads per quartet before the Schwarz test. Gain over 4.2:
at most the diagonal-bin-pair gap (0.2× on naphthalene, 0.9× on B2H6),
bought at the price of generator, register and per-quartet-latency risk.
Not recommended.

### 4.4 Scheme 3: unique atom quartets (the plan's candidate 3)

*Math.* Weight |orbit of the atom pair-of-pairs|, all shell quartets of a
representative atom quartet are evaluated (derivation in sec. 1). Exact, and
it needs no pair-orbit table, only the atom permutation (`geom%perm`).
Measured 3.84× (naphthalene), 3.62× (benzene), 3.90× (pyrene), 3.61×
(B2H6): between s2a and s2. It also loses the stabilizer's share inside an
atom quartet (all shell quartets of the rep atom quartet are done).
*Mechanism.* The condition is again per quartet unless the bins respect atom
pairs. Making them respect atom pairs means one bin per (type, size, atom
pair): measured 3 513 groups vs 118 bins (naphthalene), 6 845 vs 121
(pyrene), 532 vs 70 (B2H6): up to 30× more bin pairs (the square of that in
segments before admission) and tiny segments, the opposite of the
architecture's point (whole cross products admitted or rejected by an integer
test, big uniform launches). The alternative, the per-quartet test, has the
kernel and generator cost of 4.3 with less gain. Not recommended.

### 4.5 The engines around the bins

- **Rys vs HGP.** Both receive the same segment arguments from
  `launch_classes`; the plan change in 4.2 is invisible to them. The chunked
  Rys kernels digest per chunk with launch scalars `jfac`, `kfac`
  (weights multiply the same way). No difference between paths.
- **MD J (`trc_mdj`).** Its Hermite stage is a full ordered double loop over
  primitive pairs (bra P, ket Q; "every nranks-th ket chunk") with no
  ab↔cd fold. Restricting the **bra** list to representative pairs and
  scattering the result J_ab to the images (or weighting by osz and relying
  on the final symmetrization) gives the full effective |G| gain (ideal, not
  ideal/2) for the J that dominates f-containing bases. The ket stage stays
  replicated. Changes are in `order_pairs`/`mdj_coulomb` (bra range) and the
  `mdj_*` kernel generator is **not** touched (the bra stage already runs per
  bra shell pair, `scripts/gen_md_j.py`). Orbit weights: `jmat` of a rep bra
  pair times osz, folded into the existing `jfac` of `add_md_j`.
- **MD exchange (`trc_mdg`, g classes above `TERCO_RYS_G_LTOT`).** Windows
  of sorted pairs and a Schwarz staircase. The bra block restricted to
  reps with a per-window scale (the type already has `trc_mdg_scales_t`);
  by the idempotence argument it can simply stay full in the first PR.
- **Cartesian vs spherical, decontraction.** Skeleton in the Cartesian
  working basis; symmetrize once in the caller's AO space after the fold.
  Both verified. Orbits are taken over the *working* (decontracted) shells;
  `sym_ao_build` already supports the split basis (shell map positional per
  atom), and the prototype ran on decontracted cc-pVDZ throughout.
- **Several GPUs / MPI.** Ranks stride over the (weighted) work items of the
  same deterministic plan; every rank must build the same orbit table and rep
  order (deterministic: min pair id). Allreduce the skeleton G, then
  symmetrize once. Load balance is that of the item stride, unchanged.
- **Incremental Fock.** See sec. 1; `scf_incremental` needs no change; the
  accumulated reference is a skeleton matrix, never symmetrized until the
  fold. A restart that reads a symmetrized `gacc` is also fine (idempotence).
- **Gradients (out of scope).** The skeleton gradient is exact the same
  way: grad = (1/|G|) Σ_g P_g grad_skel with P_g the atom permutation and
  component signs (`sym_scf_grad_symmetrize` exists). What would break: the
  gradient plan (`grad_bins` calls `plan%build` directly and compacts
  survivors into a list, then 32 private gradient copies per thread) would
  need the same rep-first bins and weights; the weights are launch scalars
  there too. The response/CPHF paths use non-invariant densities and must
  stay on the full path.
- **DF-JK, XC grid.** Not touched (S16).

### 4.6 S14.4: the MD engines as built

Both McMurchie-Davidson engines restrict their bra to the orbit
representatives on the skeleton path (`eri%fock` and `eri%fock_resident` with
`skeleton=.true.`); with symmetry off, or the skeleton off, they run exactly
as before. No generator or kernel file changed: the restriction is in the
host code that sets up their launches (`trc_mdj`, `trc_mdg`, the dispatch in
`trc_eri`). The pair orbit weights come from `bins_rep_weights` (the bins'
representative sections, per pair).

- **MD J.** `mdj%build` given the weights appends a *copy* of each
  representative's primitive pairs after the ordinary ones (the exchange engine's
  view of the lists is unchanged), with `c` times the orbit size, in their
  own lp blocks of the Hermite arrays. The skeleton `coulomb` runs the
  Hermite stage with P over the copies and Q over all pairs (same kernels,
  other offsets and strides) and the bra stage over the copies' shell pairs;
  the ket stage is untouched. The result is the skeleton J, the sum over
  representatives of |orbit| J_p. The (P,Q) screen's bound I_P is formed from
  the plain `c`, before the weight. The Hermite stage has no ab-cd fold, so
  nothing is lost to a half weight: the gain is pairs / representatives,
  twice scheme 1's.
- **MD exchange.** The sorted groups (type, primitive-count bucket) are
  orbit-closed (an orbit shares its type, its primitive count and its
  symmetrized bound), so a group is split into representative groups, one
  per orbit size, in the group's own bound order (so the staircase and its
  early exit are intact), with their own E matrices. A skeleton window pairs a
  representative group (bra, weight = orbit size) with a whole group (ket);
  of two different groups of one type only the second-not-above-first pair is
  run as before, and a group against its own parent keeps the extra 1/2.
  Both J and K of the window carry the weight (omega > 0 forms J here as well).
  The fold of ab and cd is kept at group level (distinct groups, the cut is
  free), lost only on a group against itself; the choice of the cheaper
  side is not available, as swapping bra and ket changes the class.
- **Counters.** `eri%md%npq`: (P,Q) primitive pairs handed to the Hermite
  stage by the last `coulomb`. `eri%mg%nquart`: bra x ket shell pairs of
  the windows the last `fock` ran. A silent fall back to the full engine
  leaves them unchanged; `check_sym_fock_md/_mdrule/_mdg` require them nonzero
  in the full build and smaller in the skeleton build.
- **Not covered.** `fock_resident_many` (several densities) runs the engines
  per density through `fock_resident_lent`, which now honours `opts%skeleton`;
  the many-density kernel path itself has no skeleton option yet (S14.3).
  The `mg` object is rebuilt once, with the representative groups, the
  first time a skeleton build asks for them; the E matrices of the
  representative groups are the extra memory.
- **Device (S15).** New device arrays: `rsp_ao1/ao2/k0/nk` (the copies' shell
  pairs) in `mdj_build`'s and `mdj_release`'s data lists; `mp_*` and `mp_ip`/`mp_nq`
  are longer by the copies, and `mp_c`'s copy slice is updated after
  weighting. The `mg` side adds only host index lists (`gl` is longer). Not
  run on a GPU.

Measured (host, serial, `parallel=.false.`; symmetrized skeleton against the full G,
relative gate 5e-11·max(1,|G|) as in S14.2). Work counters, full -> skeleton:

| system, basis (engine) | full | skeleton | reduction | max err (abs), in brackets max G |
|---|---|---|---|---|
| water C2v, uncontracted l=0..3 (MD J, TRC_MDJ=1) | 164 836 | 99 876 | 1.65x | 5.8e-13 (196) |
| B2H6 D2h, uncontracted l=0..3 (MD J) | 5 683 456 | 1 535 296 | 3.70x | 2.6e-12 (233) |
| naphthalene D2h, d-aug-cc-pVDZ (MD J by the size rule) | 957 902 500 | 247 878 550 | 3.86x | 6.8e-11 (164) |
| water C2v, cc-pVQZ (MD J) | 1 306 449 | 833 247 | 1.57x | 3.2e-14 (50) |
| water C2v, cc-pVQZ (MD exchange quartets) | 3 303 | 2 548 | 1.30x | 3.2e-14 |
| B2H6 D2h, uncontracted l=0..4 (MD J) | 13 402 921 | 3 598 763 | 3.72x | 1.9e-12 (228) |
| B2H6 D2h, uncontracted l=0..4 (MD exchange quartets) | 1 778 423 | 502 742 | 3.54x | 1.9e-12 |

The naphthalene row's error is 4e-13 relative (the full J+K of 1e9 primitive
pairs). Without the orbit weight the same checks fail by 1.6 (g) and 48 (f);
that is, the weights, not the restriction alone, are what the checks test.

## 5. Recommendation

**Scheme 2a** (bin-ordered unique bra pairs, closed bins, reps first,
weights as launch scalars), for these reasons.

1. Gets 79-96% of the per-quartet optimum s2 (3.79× vs 3.99× naphthalene,
   3.53× vs 3.93× benzene, 3.70× vs 3.98× pyrene, 3.85× vs 4.00× pentacene, 3.34× vs 4.24× B2H6, 1.47×
   vs 1.63× water) and is 1-8% behind s3 (which costs 30x the bins, 4.4).
2. **Zero generator changes** and zero change of the kernels, the
   dispatchers, `rys_dispatch`/`pc_dispatch` or the GPU data layout beyond
   host-built, already-existing descriptors (rectangle segments with sub-range
   offsets, the tile mechanism). Compile-time cost: none. S14's risk is then
   host-side planning, which the check oracle covers completely.
3. Composable and staged by the idempotence property: PR boundaries can
   restrict one engine at a time and the others stay full.
4. The gain is what the physics allows at the shell level: on planar
   molecules ideal = |G|/2, and no scheme beats it without an AO-component
   sparsity change.

Expected benefit (quartets, K and J launches alike): ~3.5-3.85× fewer
quartets for planar D2h at cc-pVDZ (pentacene 3.85×), ~3.3× nonplanar D2h, ~1.5× C2v water-like
(small molecules with many atoms on the elements gain little). Real GPU
time follows the quartet count only to the extent that the per-segment
overhead is small (more, smaller segments); that is for S15 to measure.

Implementation plan for S14 (each PR keeps `symmetry=off` bit-for-bit, goldens
unchanged):

1. **Orbit tables** (`trc_sym_pairs`, host): from `trc_sym_ao` and
   `geom%perm`, the working-basis shell map (decontracted), pair images,
   representatives, orbit sizes, symmetrized Q and dsh. Unit check:
   closure, `Σ w` identity against the full enumeration (the prototype's
   integral-free check), orbit sizes for water/B2H6/naphthalene.
2. **Symmetry-aware bins and plan** (`trc_bins`, `trc_bin_tiles`,
   `trc_binkernel`/`launch_classes`): closed bins, reps first, weights as
   launch scalars, skeleton entry on `trc_eri_t`, `check_sym_fock`: for every
   group's test molecules × Cartesian/spherical × J/K/J+K × tiles on/off ×
   incremental on/off, symmetrized skeleton vs full to 1e-11·max(1,|G|),
   plus the l = 0..3 set (`TERCO_LMAX=3`). Host build, then OpenMP 4 threads.
3. **SCF wiring**: `scf_fock` symmetrizes after the fold when symmetry is on;
   central SCF oracle for all groups with the skeleton path on; UHF and
   range-separated (long-range exchange pass uses the same plan).
4. **MD J and MD exchange** bra restriction (required before symmetry becomes
   the default, decision 6; correct without it).
   Done (S14.4); what was built and measured: section 4.6.
5. Device mapping list (new arrays and their `enter data`/`exit data`) for
   S15: `rep` ranges and weights are host-only in design (segment
   descriptors), so the list is expected to be short.

## 6. Open questions for the maintainer

1. Accept the shell-level limit (planar molecules: effective |G|/2)? A
   further 2× on planar D2h needs AO-component parity skipping in the
   generated kernels (parity-zero integrals computed and discarded today),
   i.e. a generator project of its own.
2. Does default-on under `symmetry=on` apply to the C2v/Cs-size cases where the
   gain is ~1.5×, or only above a pair-count where it pays (the extra
   launches)? Suggest: always on, flag to disable, decide from S15.
3. Orbit-closed bin splitting past BIN_CAP can exceed 32768 pairs by < |G|;
   is changing the cap semantics to "at most BIN_CAP + |G|" acceptable (the
   segment-size guard has the margin: 32776² < 2³¹)?
4. Exchange tiles over rep ranges: the prototype's tile test says tiles are
   inert on compact molecules. An extended D2h test molecule (long acene,
   density screened) is needed to check tiles; OK to add pentacene/hexacene
   to the oracle?
5. Gradient and response are deferred; confirm the skeleton gradient
   (same plan, atom-permuted components) is the intended S16 route.
6. MD exchange bra restriction: stage it behind the bins (idempotence makes
   it safe), or require it before default-on for g bases?

## 6a. Maintainer decisions (2026-11-07)

The design (scheme 2a) is approved. The answers to section 6:

1. **Yes**: the shell-level limit is accepted (planar molecules gain about |G|/2).
   AO-component parity skipping is a separate generator project, not part
   of S14.
2. **Yes, always on**: under `symmetry=on` the skeleton build is used for
   every group, with a flag to disable it. Revisit from S15's measurements.
   *Revisited*: S15 measured a loss below roughly twenty atoms, W2 (section 4.2b)
   removed most of it, and what
   remained is now gated by the size rule of section 4.2a -- "always on"
   became "on where it was measured to pay, with the flag either way".
3. **Yes**: orbit-closed bin splitting may exceed BIN_CAP by fewer than |G|
   pairs. The cap semantics become "at most BIN_CAP + |G|".
4. **Yes**: an extended D2h molecule (pentacene, and hexacene if it runs in
   reasonable time) joins the S14 oracle to exercise the exchange tiles.
5. **Yes**: the skeleton gradient (same plan, atom-permuted components) is
   the intended S16 route. Response stays on `fock_nosym`.
6. **Required before default-on** (clarified by the user on 2026-11-07): the
   MD exchange bra restriction is staged behind the bins as S14 PR 4. It must
   be done before symmetry, and with it the skeleton build, becomes the
   default, for every basis and not only for bases with g functions. Until
   PR 4 lands, the MD exchange engine runs unrestricted, which is still exact
   by idempotence, and `symmetry` stays off by default.

## 7. Caveats and what the prototype does not show

- All counts are CPU, enumerated from the bins, with the production
  admission test. Wall times, occupancy and register effects are GPU
  questions (S15) and are not measured here.
- The tile numbers use a proxy for the tile test (exact shell sets of the
  64-pair tiles instead of the production block hulls), so production keeps at
  least as many tile pairs as counted.
- The compile-time numbers for 4.3 are estimates from the size of the
  emitted screen text, not a build.
- Pair-of-pairs weights beyond the abelian groups (non-abelian, orbits
  that are not G/Stab products) are out of scope; the project only needs
  abelian groups.
- Planar molecules: the z-mirror is a pure sign operation; the prototype's
  "ideal" column accounts for it by counting distinct shell permutations.
- Test molecules: pyrene and pentacene coordinates are idealised
  (CC 1.40, CH 1.08 Å on a hexagonal lattice), generated for this spike
  (`test/pentacene_d2h.xyz` is in the tree; the pyrene file is fork-only history, not in this tree).
- Run-to-run noise floor: the OpenMP and device builds accumulate the Fock matrix in a
  run-dependent order, so identical runs differ by about 1e-13 to 1e-12 Eh and, for ~1000 Eh
  systems, the commutator floors at a few 1e-10. Checks must not compare runs bitwise or
  converge below that floor (oracle SCFs use `conv_diis >= 1e-9`, `conv_energy >= 1e-10`; the
  energy is variational, so the 1e-10 gates still hold).
- Reproduce: the prototype and its command line are kept in the fork's history (AlbertThie/terco commit 496bef9).
