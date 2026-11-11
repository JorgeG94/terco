# Symmetry-unique shell-pair list: design (sub-project of S14)

Branch `sym/unique-pairs`, off `symmetry` at ee2cb3d (S0-S14.4). Status: the
design study is done, **stage U1 is implemented** (`src/symmetry/trc_sym_pairs.F90`,
checked by `check_sym_unique` and `check_sym_fock_nouniq`; section 7), stage U2 (the kernel reads the image) is designed and **not built**
(section 4, gated on a GPU measurement). The measurements are CPU; nothing
here was run on a GPU.

## 0. The idea and the answer in one page

The maintainer's idea (2026-11-07): "We dispatch batches of shell pairs that are
assembled as quartets. We assemble those pairs at startup after screening. Only
screened pairs are added to the list. Screened, symmetry-unique (point group)
pairs could be added only, and then in the Fock the necessary missing bits can
be done."

- **What is missing** (section 1). Evaluating quartets with *both* pairs
  orbit representatives is provably not enough: it reaches 89%, 47%, 32%, 28%,
  26% of the quartet orbits for water, B2H6, benzene, naphthalene, pentacene.
  Every other quartet orbit needs one *image pair*: the scheme-2a plan (bra =
  representatives, ket = the whole bin) is already the minimal form of that,
  and no choice of orbit representatives, transformed density or transformed
  Fock contribution removes the need. The quartet set, and so the quartet
  count, of scheme 2a does not change; only where the image's pair data comes
  from.
- **An image pair's data is its representative's, exactly** (section 2): same
  zeta, coefficient, exponent; the product centre P and the shell centres A, B
  with the axes flipped by the operation's mask. No component sign enters. The
  reconstruction is bit-identical to the directly built data (max deviation
  0.0 on every system).
- **What it buys** (section 3). Measured on the CPU: primitive pairs formed and
  stored 4.5-5.2x fewer (D2h), the pair object 2.6-3.9x smaller, the whole
  planning 2.4-3.8x faster (the Schwarz pass 4x), the host pair data 1.8-2.3x
  smaller with the device tables unchanged (U1), up to 3-4.5x smaller on the
  device with the kernel change (U2, modelled). Absolute numbers are small
  (0.09 s -> 0.04 s of planning at 792 shells): the gain is O(nshell^2)
  integer work and memory at the start, not Fock time.
- **Recommendation** (section 6): build U1 now (done, no generator change, no
  GPU risk). Build U2 only if a GPU measurement shows the bin-interleaved
  device tables are what limits system size: it is a generator change with the
  parity-skip project's files in play, and its gain is bounded by the
  O(nshell^2) tables that no pair-list scheme touches.

## 1. What is missing: the math

Notation as docs/SYMMETRY_FOCK.md section 1: G the abelian group, R(g) the
signed permutation of the AOs, Dig_q the contribution of the shell quartet q
(its eight orderings) to G, linear in the density D, with
Dig_{gq}(D) = R Dig_q(D) R^T for an invariant D. A quartet is an unordered
pair of shell pairs {p, q}; G acts on it by {gp, gq}.

**Claim A (what has to be evaluated).** For any G-closed quartet set S,
G(D) = (1/|G|) sum_g R G_skel R^T with G_skel = sum over *one member of each
G-orbit of S* of |orbit| Dig_q. The member can be chosen freely, and the number
of orbits is, by Burnside on the pair-of-pairs,

    N_orb = (1/2|G|) sum_g ( fix(g)^2 + npair ),   fix(g) = #{pairs p : gp = p}.

N_orb is the least number of quartets any exact scheme of this kind evaluates
(it is the S13 scheme s2, 3.99x on planar D2h).

**Claim B (a representative bra always exists, the ket then runs over images).**
Every quartet orbit has a member {p0, q} with p0 an orbit representative of the
pair-orbit of p (apply the group element that sends p to p0). Two members
{p0, q}, {p0, q'} lie in one orbit iff q' = h q with h in Stab(p0). So the kets
needed for a bra p0 are the Stab(p0)-orbit representatives of *all* pairs q,
and every pair q = h r (r its orbit's representative, h in G/Stab(r)) needs the
image pair h r, not r alone. This is scheme 2a: bra = representatives (weight
|orbit|), ket = the complete bin; it evaluates the Stab(p0)-orbits' union, which
is more than N_orb only by the diagonal bin pairs (3.79x against 3.99x).

**Claim C (refuted: representatives on both sides).** The quartets {p, q}
with p and q both representatives are nrep (nrep + 1) / 2 in number, at most one
per pair of pair-orbits, and distinct ones lie in distinct G-orbits (if g sends
{p, q} to {p', q'} with all four representatives, orbit(p) = orbit(p') forces
p' = p, then q' = q). A pair of pair-orbits (P, Q) holds |P||Q| / |G|
quartet orbits when the action is free, so all of them are reached only if
|P||Q| <= |G|; for any two large orbits (the generic case) most are not. Measured
(Burnside count against nrep (nrep + 1) / 2, cc-pVDZ, a spike program, fork-only history, commit 496bef9; not in this tree):

| system | pairs | orbits | quartet orbits N_orb | rep x rep | fraction reached |
|---|---|---|---|---|---|
| water C2v | 171 | 126 | 9 036 | 8 001 | 0.885 |
| B2H6 D2h | 866 | 288 | 88 842 | 41 616 | 0.468 |
| benzene D2h | 3 309 | 942 | 1 391 601 | 444 153 | 0.319 |
| naphthalene D2h | 7 279 | 1 911 | 6 637 551 | 1 826 916 | 0.275 |
| pentacene D2h | 21 433 | 5 495 | | 15 100 260 | 0.263 |
| hexacene D2h | 26 191 | 6 657 | | 22 161 153 | 0.258 |
| decacene D2h | 45 223 | 11 429 | | 65 316 735 | 0.255 |
| tetradecacene D2h | 64 255 | 16 201 | | 131 244 301 | 0.254 |

(the fraction tends to 1/|G_eff| = 1/4 on planar D2h).

**Claim D (transformed densities or Fock contributions do not help).** The
missing quartets are {p0, h r} with h not in Stab(p0) Stab(r). Applying h to
the whole quartet gives {h p0, r}: the same two pair-orbits in another relative
orientation, an *image* on one side again. The value (p0 | h r) is a different
number from (p0 | r) (the Gaussian-product centres lie at different
separations), and the symmetry relation that does exist,
(gp | gq) = s (p | q), links only quartets of one orbit, which is exactly what
the symmetrizer uses. A transformed density can only re-route where a computed
integral is digested, and it is already used: G_skel is symmetrized once. So
the integral of every missing orbit has to be computed from the pair data of
the representative on one side and an *image* on the other. "One side
complete" of scheme 2a is required (up to Stab(p0)); what is free is only how
the image's pair data reaches the kernel.

**What the image data is** (checked numerically, section 3). Group operations
are diagonal sign matrices g = diag(+-1, +-1, +-1) with g^2 = 1. The image of
the shell pair (a, b) is the pair (ga, gb) of two shells of the basis, each
with its own centre g A, g B, the same exponents and contraction. Its primitive
pairs are the representative's with zeta, the combined coefficient c
(K_ab = exp(-mu |AB|^2) is invariant) and the i-side exponent unchanged and

    P' = g P,   A' = g A,   B' = g B,   (A' - B') = g (A - B).

The Cartesian component signs `ao_sign` do **not** enter: the image shell is a
shell with the standard Cartesian components at its own centre; those signs
relate Fock *matrix elements* of equivalent shells and live in the symmetrizer.
The pair is stored as (ga, gb) in the bins' own order (higher l first, index
order breaking a tie), which can be the reverse order of the representative's:
the reverse order of the representative is then needed (`sym_pairs_wanted`).
The Q of an image is the representative's (Q is a function of the orbit), and
screening by Q := orbit maximum, dsh := orbit maximum is unchanged.

## 2. Reconstruction options

All of them leave the quartet set of scheme 2a unchanged: same bins, same
segments, same weights, same counts (`check_sym_fock` still asserts
14627/9968, 360415/107763, 4741001/1342393, 20671728/5451483).

### (a) In the kernel's pair read (U2)

The ket of a bin pair is the set of images of the representatives of bin B. The
plan emits it as segments (ket section = representatives of one stabilizer
type, operation g = one coset representative of G/Stab); inside a segment the
operation is **uniform**, so the sign is a segment constant, not per-thread
data. Kernel changes (HGP skeleton shared by `gen_perclass.py` and the Rys
splice of `gen_rys.py`, the dispatchers `pc_dispatch`/`rys_dispatch`):

- three signs sgx, sgy, sgz from a new segment array `sG(seg)` (one integer load);
- the ket reads `qp_px(kq)`, `qp_py(kq)`, `qp_pz(kq)`, `sq_ax/ay/az(slotcd)` and
  `sq_abx/aby/abz(slotcd)` become `sg*...`: in the (dd|dd) kernel these are
  the 24 lines 292 and 366-377 of `trc_pc_k2222.F90`; the bra stays;
- the ket shells `sk, sl` (digestion and the dsh screen) become
  `shmap(sp_i(slotcd), g)`: two integer loads from a table of
  nshell x nop;
- `sp_q` is the representative's (no change).

Estimated cost, not measured on a GPU: +3 live scalars (the signs; the ket
loop already holds dozens), +1 load per quartet for `sG`, +2 for the shell
map, +6 multiplications per ket primitive. Host compile-time measurement (the
only one available here): the (dd|dd) kernel with all 24 read sites patched
and the mask decoded from a segment integer, `gfortran-15 -O3`:
46.6 s -> 41.3 s, object 373 304 -> 373 400 bytes, i.e. no change beyond noise
(one run each). The nvfortran build of the heavy classes
(15-20 min today) and their register use (some already at 255 with a local
stack) is the unmeasured risk: the signs are dead after the ket reads, so the
expectation is no growth, but that is an expectation. A full regeneration of
~230 generated files and one full kernel rebuild are required (the dummy lists
of every kernel grow by `sG` and the shell map), so this is a one-time cost,
not a growth in code.

Plan consequences: the number of segments grows by about the effective group
order times the number of stabilizer types (<= ~30x of the 5 186 bin pairs of
naphthalene), each larger than a tile; the prefix-sum search of a thread gets
~5 steps longer. The weight and the diagonal-bin-pair rule are as in 2a.

| interaction | with (a) |
|---|---|
| screening (orbit-max Q, dsh) | unchanged: Q read from the representative, dsh read at the image shells |
| exchange tiles | tiles over representatives (as now) x operation; the tile reach uses image shell sets, built once from `shmap` |
| incremental Fock | unchanged (skeleton matrices accumulate) |
| MD engines (`trc_mdj`, `trc_mdg`) | read the basis, not `pairs`/bins tables; their ket needs all pairs, expanded host-side from the bins' representatives list (they build their own arrays already) |
| parity skipping (`sym/parity-skip`) | the dead-axis mask of a quartet (axes on which all four centres have zero coordinate) is invariant under every operation, so the mask computed for the representative quartet is valid for its image; the sign flips commute with it |
| DF-JK | not affected (3-centre pair lists, S16) |
| gradients | `gen_rys_grad.py` reads the same tables: either complete the pairs (as U1 does) or give the gradient kernels the same segment sign; S16 decides |
| risks | GPU registers and the segment count unmeasured; a second variant of every kernel dummy list; merge conflicts with the parity project in the generators |

### (b) Per-launch expansion into a device scratch buffer

No generator change. A launch covers a *class run* of segments of many ket
bins, so the buffer has to hold the images of every ket bin of that run: up to
all the bin tables of the ket type (the s-type pairs dominate the pair count),
i.e. no saving of memory unless launches are cut per ket bin (~50x more
launches). Traffic: one extra write and read of the whole ket data per launch
group, on top of the kernel's own reads. **Not recommended.**

### (c) Expand once at build time into the existing tables (U1, implemented)

`hgp_build_pairs` forms the primitive pairs of the wanted ordered shell pairs
only, Schwarz runs over them, and `bins_build_pair_quartets` fills the
bin-interleaved slots of an image from its representative with the signs.
Everything downstream (kernels, plan, tiles, MD engines, parity skipping, the
skeleton scheme) is untouched; the device holds the same tables as before. What
it saves is everything *before* the tables: the primitive-pair construction,
the Schwarz pass (the O(npair) device launch), the host `pairs` object and the
`keep` scratch of `hgp_build_pairs`. What it does not: the bin-interleaved
tables and the bins on the device, and every nshell^2 array.

The paths that read `eri%pairs` through the bins (the non-symmetric Fock build,
the gradient) call `complete_pairs` first, which appends the missing pairs'
data (host work, once, device copy replaced).

### (d) Something better?

Two refinements, neither built: (i) a device kernel that expands the
representative tables into the interleaved tables at build, so the host never
holds the full tables (only a memory-peak gain); (ii) the O(nshell^2) tables
(`pairs%off/n` 8 B x nshell^2, `dsh`, the `sym_pairs_build` pair index, `qs`)
dominate once the pair list is unique: 5 MB of 17 MB for tetradecacene already.
Compressing them (sparse by the screened pairs) is orthogonal to the symmetry
and larger than anything left in the pair lists at large size.

## 3. Measurements (CPU, cc-pVDZ, tau = 1e-10, decontracted working basis)

Programs: the prototype `spike_unique` (counts, memory model, reconstruction,
exactness; fork-only history, commit 496bef9; not in this tree) and
`test/check_sym_unique.F90` (the product path). The logs are measured
output, not in the tree.

### 3.1 Reconstruction and exactness

For every non-representative pair of the bins (e.g. naphthalene 7 455
primitive pairs, benzene 3 423): zeta, c, a bit-identical to the direct
build; max |P - g P_rep|, |A,B - g A,B| = 0.0 (bitwise, B2H6, benzene,
naphthalene); the primitive counts agree. G: the production kernels run on tables rebuilt from
reconstructed pairs, symmetrized skeleton against the full G, random invariant
density, gate 5e-11 max(1, |G|):

| system | max abs G | symmetrized error J+K | K only | raw skeleton error |
|---|---|---|---|---|
| water C2v | 23.6 | 2.5e-14 | 8.0e-15 | 11.5 |
| B2H6 D2h | 38.3 | 8.5e-14 | 1.9e-14 | 58.1 |
| benzene D2h | 51.3 | 1.8e-11 | 1.2e-12 | 81.8 |
| naphthalene D2h | 77.5 | 9.5e-12 | 1.6e-12 | 143.5 |

(gates 1.2e-9 ... 3.9e-9). In the product code the same is checked by
`check_sym_fock` (all groups, Cartesian and spherical, J/K/J+K, tiles, incremental,
converged and random densities) and `check_sym_unique`, section 7.

### 3.2 Pair list, primitive pairs, pair data

Primitive pairs: the merged build's `pairs` object vs the unique build
(`eri%upair_npp_stored`), and what the bins use (`upair_npp_bins`); shell pairs
formed of nshell^2; the `keep` scratch of `hgp_build_pairs` (maxnp^2 per pair
slot, here maxnp = 3) shrinks with the formed pairs.

| system | shells | pairs in bins | orbits | shell pairs formed | prim. pairs merged -> unique | ratio |
|---|---|---|---|---|---|---|
| water | 18 | 171 | 126 | 139 / 324 | 554 -> 246 | 2.3 |
| B2H6 | 42 | 866 | 288 | 434 / 1764 | 2 954 -> 692 | 4.3 |
| benzene | 90 | 3 309 | 942 | 1 820 / 8 100 | 9 414 -> 2 082 | 4.5 |
| naphthalene | 144 | 7 279 | 1 911 | 3 971 / 20 736 | 19 960 -> 3 862 | 5.2 |
| pentacene | 306 | 21 433 | 5 495 | 19 650 / 93 636 | 56 882 -> 11 547 | 4.9 |
| hexacene | 360 | 26 191 | 6 657 | 27 122 / 129 600 | 69 264 -> 13 965 | 5.0 |
| decacene | 576 | 45 223 | 11 429 | 69 442 / 331 776 | 118 792 -> 23 887 | 5.0 |
| tetradecacene | 792 | 64 255 | 16 201 | 131 346 / 627 264 | 168 320 -> 33 809 | 5.0 |
| tetradecacene 6-31G | 354 | 19 611 | 4 967 | 26 183 / 125 316 | 142 776 -> 29 788 | 4.8 |

Quartets evaluated: identical to scheme 2a by construction (asserted): 9 968
(1.47x), 107 763 (3.34x), 1 342 393 (3.53x), 5 451 483 (3.79x) fewer than
14 627 / 360 415 / 4 741 001 / 20 671 728; pentacene 39 110 378 (3.85x).

Pair data in MB (host + device; `pairs` object at 96 B per primitive pair plus
the 8 B nshell^2 off/n tables; bin-interleaved tables at 48 B per slot-primitive
plus 52 B per pair; bins 20 B per pair). "U1" is measured sizes of the
implemented build; "U2" is the model with representatives-only tables and bins
(prototype `spike_unique`, fork-only history, `[model]`):

| system | pairs object merged -> U1 | total merged -> U1 -> U2 [model] |
|---|---|---|
| benzene | 0.97 -> 0.26 | 1.44 -> 0.73 -> 0.34 |
| naphthalene | 2.08 -> 0.54 | 3.09 -> 1.55 -> 0.70 |
| pentacene | 6.21 -> 1.86 | 9.13 -> 4.78 -> 2.31 |
| hexacene | 7.69 -> 2.38 | 11.25 -> 5.94 -> 2.93 |
| decacene | 14.06 -> 4.95 | 20.18 -> 11.07 -> 5.89 |
| tetradecacene | 21.18 -> 8.26 | 29.87 -> 16.96 -> 9.60 |
| tetradecacene 6-31G | 14.71 -> 3.86 | 19.53 -> 8.68 -> 4.35 |

The reduction of the pair object falls with size (3.9x -> 2.6x) because the
8 B nshell^2 off/n tables (5 MB at 792 shells) are not reduced: they become the
floor. `dsh` (8 B nshell^2: 5.0 MB), the Q triangle, and the transient
`sym_pairs_build` index (3.3 MB) are likewise untouched. Device footprint of
U1 equals the merged scheme's except for `pairs` (7 of the 11 arrays of the
object are mapped).

### 3.3 Planning time (`eri%t_build`, serial host, CPU; ms, minimum of 3, identical to the ms)

| system | merged total | unique total | pairs stage | Schwarz stage | bins | quartet tables |
|---|---|---|---|---|---|---|
| naphthalene | 10 | 4 | 1 -> 1 | 8 -> 2 | 1 | 0 |
| pentacene | 28 | 11 | 3 -> 2 | 23 -> 6 | 2 | 1 |
| hexacene | 34 | 13 | 3 -> 3 | 28 -> 7 | 2 | 1 |
| decacene | 60 | 24 | 8 -> 6 | 47 -> 12 | 3 | 1-2 |
| tetradecacene | 88 | 36 | 13 -> 11 | 67 -> 17 | 5 | 2 |

The Schwarz stage falls 3.9-4.0x (it scales with the pairs formed). On this CPU
the Schwarz kernel is the host `do concurrent`; on the GPU it is the device
launch the pair list currently pays in full. The pairs stage falls less: it
loops over all nshell^2 ordered pairs for the cheap prune and (now) the
wanted-pair marking. All of planning is 0.1 s at 800 shells, against
seconds per Fock build: not a time that matters yet; at 10^4 shells the
O(nshell^2) prune and marking (10^8) and the Schwarz launch over the
pair list are what grow.

### 3.4 Where the extra work is

- U1: building the images' table slots (3 sign multiplications per primitive
  and a search over at most |G| operations per image pair) is inside the
  `quartet tables` stage: 1-2 ms at 800 shells, in the table above. No cost in
  the Fock build.
- U2: per-quartet cost as in section 2(a); no host-side cost; a longer segment
  array.

## 4. Candidate schemes, compared

| scheme | quartets | generator change | device pair data | planning | status |
|---|---|---|---|---|---|
| merged 2a | N_2a | none | full | full | merged (S14.2) |
| U1: unique build, tables expanded at build | N_2a | none | = merged (tables), `pairs` /3 | Schwarz /4 | implemented |
| U2: kernel reads the image | N_2a | `gen_perclass`, `gen_rys` (+ dispatchers) | /3..4.5 [model] | as U1 | designed, gated |
| (b) expand per launch | N_2a | none | scratch ~ bin tables of a run | as U1 | rejected |
| s2 (per-quartet orbit test) | N_orb (+5% over 2a) | per-quartet filter in every kernel | full | full | rejected in S13 |

## 5. Interactions common to all (U1 as built)

- *Screening*: the images' Q is zero in the Schwarz output; `sym_pairs_symmetrize_q`
  (orbit maximum) then yields the representative's Q for every member, which is
  exactly invariant (before: invariant by a maximum over rounding noise). The
  dsh maximum is unchanged. Bin keys use the representative's primitive count.
- *Exchange tiles, incremental Fock, MD engines, parity-skip*: unchanged
  (they see the same tables and bins).
- *Nonsymmetric Fock and gradient*: `complete_pairs`.
- *Symmetry off, or order-1 group*: `hgp_build_pairs` is called without `want`,
  the code takes the old branches; verified by the full suite and the goldens
  (section 7).

## 6. Recommendation and PR plan

1. **U1 (done, this branch)**: unique pair build, `complete_pairs`, counters,
   `check_sym_unique`, `check_sym_fock_nouniq`. Host files only.
2. **U2 (not built)**: only if a GPU run shows the bin-interleaved device tables
   limit the system size (the model says they are about two thirds of the
   pair data at 800 shells and the nshell^2 tables take over beyond). Needs:
   the generators and dispatchers (a), a plan with image segments
   (stabilizer section x coset operation), tile and MD changes, `shmap` on the
   device, and the parity project's files merged first. Compile-time: one full
   regeneration/rebuild; host measurement shows no change.
3. **U3 (independent)**: compress the O(nshell^2) tables.

## 7. Decisions for the fork implementation (2026-11-07)

Taken on the recommendation above; the coordinator treated the recommendation as
approved for the fork.

1. Build U1; default **on** whenever `eri%build` is given a point group of
   order > 1; `TRC_UNIQUE_PAIRS=0` restores the merged build (used by
   `check_sym_fock_nouniq`). *Jorge to confirm: default on.*
2. Do not build U2 now; keep this design. *Jorge to confirm the gate (a GPU
   memory measurement) and whether the generator/segment work is wanted at all.*
3. The images are rebuilt **once at build** into the existing tables, not per
   launch. 
4. `eri%pairs` is incomplete on a symmetric engine until `complete_pairs` runs
   (called by `fock_nosym` and `gradient`); no other reader of `eri%pairs` exists
   in `src` (checked); test programs that read it must build their engine without
   a group. *Jorge to confirm that this component semantic is acceptable for
   external users of `eri%pairs`.*
5. Wanted ordered pairs: the representative in the Schwarz order (i >= j) and in
   the order of every image's stored orientation (`sym_pairs_wanted`).
6. A rank builds the same list (deterministic canonical key); no MPI change.
7. Not changed: generators, kernels, plan, tiles, MD engines, ABI, API.

Shared files touched, for merging with the parity-skip work:
`src/eri/trc_hgp.F90` (`hgp_build_pairs` optional `want`),
`src/eri/trc_bins.F90` (`bins_build_pair_quartets` optional `orb`, new public
`bins_pair_image_source`), `src/eri/trc_binkernel.F90` (one call site: pass `orb` to
the rebuild after a resort), `src/eri/trc_eri.F90` (`build_structures`, new
`complete_pairs`, `unique_pairs_requested`, calls in `eri_gradient` and
`eri_fock_nosym`, four fields), `src/eri/trc_eri_data.F90` (`ops` in
`trc_bin_orbits_t`), `src/symmetry/trc_sym_pairs.F90` (`sym_pairs_wanted`).
Not touched: `src/api/`, `app/`, `trc_sym_run.F90`, `src/generated/`, the plan and
segment code.

## 8. Open questions for the maintainer

1. Is a ~2x smaller host pair data and a 4x shorter Schwarz stage enough to keep
   the default on, given that planning is ~0.1 s at 800 shells? (The gain is real
   only at much larger sizes; on the other hand it costs nothing in the Fock
   build.)
2. U2: GPU measurement of the device footprint of the bin-interleaved tables
   before deciding; who runs it?
3. Should the O(nshell^2) tables (`off/n`, `dsh`) be made sparse next (U3)?
4. With U2 the gradient either completes the pairs (as U1) or gets the same
   segment sign in `gen_rys_grad.py`: S16's choice.

## 9. Caveats

- All numbers are CPU/host; nothing ran on a GPU. `complete_pairs` replaces
  device arrays with `exit data`/`enter data` and has not been run on a device.
- The compile-time measurement is one gfortran run per variant on one kernel.
- U2 sizes are a model with the sizes of the production arrays.
- Acene geometries are idealised (CC 1.40, CH 1.08 Angstrom); decacene and
  tetradecacene were generated for this study (`acene.py`, not in the repository).
