# Skipping integrals that vanish by symmetry inside a shell (parity skipping)

Branch `sym/parity-skip`, off `symmetry` at 2291c9c (S14.1-S14.3). Separate
project from S14 by maintainer decision 1 in `docs/SYMMETRY_FOCK.md` §6a.
This is the **design study** (a spike in the style of S13): measurements, a
CPU prototype, candidate schemes and a recommendation. No product code was
written; `src/` and `src/generated/` are untouched. **No implementation
starts until the maintainer signs off on this design.**

## 1. Summary

- For a quartet whose four centres all have a coordinate equal to zero (the
  molecular plane of a planar molecule, or an axis), every component
  integral with odd combined parity in that coordinate is exactly zero. That
  is 34-52 % of the evaluated component integrals on the test molecules
  (45 % on naphthalene, benzene, pentacene).
- It is not only planar molecules: nonplanar B2H6 has 34 % (full build) or
  40 % (skeleton), because six of its eight atoms lie in one mirror plane
  (the expectation "small" was wrong). Mixed masks inside one bin are then
  the issue (§5.3).
- Counting zero components overstates the gain. Weighted by a per-class
  operation count, the work removed is **9-10 %** if only the digestion skips
  them, **10-11 %** if the HRR is also masked, and **30-32 %** if the
  VRR and the primitive accumulation are masked too (full per-class masked
  variants). That is a **1.43-1.46x** speed-up of the planar D2h skeleton
  on top of the orbit reduction (3.5-3.8x): about **5.0-5.6x of the ideal 8x**
  against 3.5-3.9x today. Sensitivity to the per-primitive overhead the mask
  cannot remove: 1.26x to 1.68x (§3.4).
- **Only the primitive-level skip is worth the trouble.** Skipping in the
  digestion alone buys 1.10-1.11x.
- Recommended: **scheme (a), masked per-class kernels, restricted to one
  dead axis, and only for the classes that carry the cost**, behind a
  launch-level variant key, with the plane normal canonicalised to z. Host
  compile cost per variant measured at 15-20 % of the unmasked kernel
  (smaller, not larger, objects). Details and the PR plan in §7.
- Open questions that need the maintainer, §8: the planar frame (D2h
  planar molecules come out with the normal on **x**, Cs ones on **z**),
  how many masks to build, mixed-mask bins.

## 2. The math

Let a centre have coordinate d (x, y or z) equal to 0 in the frame the basis
is built in. Let Z be the set of coordinates that are 0 for *all four*
centres of a quartet (the "dead" coordinates). Then

- PA_d = PB_d = QC_d = WP_d = WQ_d = AB_d = CD_d = 0 for d in Z, exactly;
- an Obara-Saika node [e0|f0]^(m) with odd e_d + f_d for some d in Z is
  exactly zero (the recurrence preserves the parity of the d-power because
  the only terms that change it carry a factor that is zero);
- a final component (ab|cd) with odd a_d+b_d+c_d+d_d for some d in Z is
  exactly zero.

The stabilizer view of the same fact: a symmetry operation g with sign
vector s (sx, sy, sz) that maps each of the four shells to itself multiplies
the component integral by prod_i s_i^(P_i), where P is the XOR of the four
components' power parities. So the integral vanishes unless s.P is even for
every g in the stabilizer. For a quartet of atoms that all lie in the plane
z = 0 the stabilizer contains sigma_z and the rule is P_z even. Spherical
harmonics have definite parity too (`trc_sym_ao` computes it), but the
engines are Cartesian, so Cartesian parities are what is counted.

The geometric rule (Z from the coordinates) is **at least as strong as** the
stabilizer rule. The prototype computes both; they agree on every molecule
measured (identical zero fractions in every table). They would differ only
for an atom that sits on a coordinate plane that is not a symmetry element
of the detected group.

Exactness requires the dead coordinates to be exactly 0. `sym_detect`
writes exact zeros to `r_std` for atoms on a symmetry element, and the
prototype confirms it: the largest |v| over the components it skipped is
**0.0** in every verify run.

## 3. Measurements

Tools (all CPU; nothing under `src/`):
`test/spike_parity.F90` (quartet enumeration on terco's own bins,
registered DISABLED), `scripts/spike_parity_gen.py` (operation counts per
class and mask from copies of the generator logic, and the report),
`test/spike_parity_kernel.F90` (masked kernels vs the shipped ones).

### 3.1 Method

`spike_parity` builds `eri%bins`, closes them under the group and
enumerates the quartets of two plans: **full** (the production enumeration)
and **s2a** (the skeleton plan of S14.2: bin-pair fold kept, one side cut to
orbit representatives). Quartets are counted after the exact Schwarz
product, tau = 1e-10, cc-pVDZ, working (decontracted) basis, Cartesian.
No density screening (it removes 2-12 % of quartets uniformly, S13 §2.2;
that moves no fraction here). For each quartet it records class
(la,lb,lc,ld) in canonical order, the dead-coordinate mask Z, the
stabilizer S, and the number of primitive quartets; the CSVs are in the
logs listed at the end. Component zero counts, and the costs, are computed
from the tables by `spike_parity_gen.py report`.

Cost model (operation counts, not time): per primitive quartet a fixed
overhead (prefactor, geometry, Boys: 90 + 4 L_tot, **assumed**), the VRR
statements of the sieved recurrence (counted from the text
`gen_perclass.py` really emits, asserted identical at mask 0), and the
accumulation into g1 (nv adds as shipped); per contracted quartet the HRR
statements (counted from `gen_hrr.py`) and the digestion (13 per component
plus the block loads/stores). Classes that production launches with Rys
(any f shell, and the seven heavy d classes) use an analytic Rys model
(roots, 2D tables, 3 flops per root per component). Four tiers:

| tier | what is skipped | needs |
|---|---|---|
| T1 | zero components in the digestion only | a runtime mask |
| T2 | + the zero HRR components and dead HRR terms | masked HRR/digestion text |
| T3 | + the odd VRR nodes, dead terms, accumulation (primitive level) | fully masked variants (scheme a) |
| TB | scheme (b): launch-uniform guards in the shipped kernel | zero-filled dead nodes |

### 3.2 Zero fractions and cost-weighted work (cc-pVDZ)

"Zero" = fraction of evaluated component integrals that are parity-zero.
"Removed" = share of the modelled flops. Speed-ups are 1/(1-removed).

| system (group, plane normal) | plan | quartets | zero | T1 | T2 | T3 (a) | TB (b) | x T3 | x TB |
|---|---|---|---|---|---|---|---|---|---|
| HOCl (Cs, z) | full | 152 776 | 0.459 | 0.140 | 0.162 | 0.322 | 0.278 | 1.47 | 1.39 |
| HOCl | s2a | 159 663 | 0.459 | 0.140 | 0.163 | 0.322 | 0.278 | 1.47 | 1.39 |
| C2H4 (D2h, x) | full | 192 728 | 0.489 | 0.084 | 0.094 | 0.287 | 0.233 | 1.40 | 1.30 |
| C2H4 | s2a | 64 530 | 0.519 | 0.092 | 0.104 | 0.317 | 0.258 | 1.46 | 1.35 |
| benzene (D2h, x) | full | 4 648 438 | 0.449 | 0.090 | 0.101 | 0.298 | 0.243 | 1.42 | 1.32 |
| benzene | s2a | 1 315 877 | 0.451 | 0.091 | 0.102 | 0.300 | 0.245 | 1.43 | 1.32 |
| naphthalene (D2h, x) | full | 20 255 420 | 0.450 | 0.095 | 0.107 | 0.307 | 0.252 | 1.44 | 1.34 |
| naphthalene | s2a | 5 344 135 | 0.450 | 0.096 | 0.107 | 0.308 | 0.253 | 1.45 | 1.34 |
| pentacene (D2h, x) | full | 145 127 813 | 0.453 | 0.101 | 0.113 | 0.316 | 0.261 | 1.46 | 1.35 |
| pentacene | s2a | 37 752 950 | 0.453 | 0.101 | 0.113 | 0.316 | 0.261 | 1.46 | 1.35 |
| B2H6 (D2h, mixed) | full | 357 274 | 0.338 | 0.055 | 0.062 | 0.189 | 0.153 | 1.23 | 1.18 |
| B2H6 | s2a | 106 829 | 0.398 | 0.068 | 0.077 | 0.231 | 0.188 | 1.30 | 1.23 |
| H2CO cc-pVTZ, f (C2v, x) | full | 426 234 | 0.553 | 0.158 | 0.166 | 0.437 | 0.385 | 1.78 | 1.63 |
| H2CO cc-pVTZ | s2a | 306 822 | 0.586 | 0.169 | 0.177 | 0.470 | 0.414 | 1.89 | 1.71 |

(The f case is formaldehyde, planar C2v, cc-pVTZ, built with `TERCO_LMAX=3`:
`test/basis_sets/cc-pvtz.json` has no chlorine, so HOCl was not possible.
Its f classes are Rys and use the analytic Rys model, so treat the f row as
the less certain one. Orbit reduction there is 1.39x: C2v has effective
order 2.)

Full-vs-s2a quartet ratios (orbit reduction, S14.2): naphthalene 3.79,
benzene 3.53, pentacene 3.84, C2H4 2.99, B2H6 3.34, HOCl 0.96. The
parity gain multiplies it: naphthalene 3.79 x 1.45 = **5.5x**, benzene
5.0x, pentacene 5.6x of an ideal 8 (3.5-3.9 today); HOCl 1.41x of an ideal
2.

Per-class breakdown, naphthalene s2a (the 14 most expensive classes; tag
(la lb|lc ld), share of modelled flops, zero fraction, work removed at T1,
T2, T3 and TB):

```
  class    engine  quartets   share   zero   rm_T1  rm_T2  rm_T3  rm_TB
  (ps|ss)  HGP     1513741   0.093   0.334  0.047  0.047  0.082  0.075
  (ds|ps)  HGP      433152   0.091   0.445  0.110  0.110  0.347  0.277
  (pp|ps)  HGP      202476   0.090   0.482  0.084  0.101  0.353  0.268
  (ps|ps)  HGP      634406   0.086   0.445  0.085  0.085  0.227  0.189
  (dp|ps)  HGP      141174   0.085   0.482  0.124  0.149  0.414  0.329
  (dp|pp)  Rys       23496   0.069   0.494  0.078  0.078  0.365  0.317
  (dp|ds)  Rys       52044   0.060   0.482  0.130  0.130  0.347  0.301
  (ds|pp)  HGP       74960   0.058   0.482  0.097  0.116  0.415  0.305
  (dp|dp)  Rys        9568   0.044   0.494  0.101  0.101  0.401  0.364
  (ds|ss)  HGP      551464   0.041   0.334  0.077  0.077  0.165  0.137
  (pp|ss)  HGP      264297   0.040   0.445  0.075  0.090  0.181  0.156
  (ss|ss)  HGP     1063842   0.037   0.000  0.000  0.000  0.000  0.000
  (dp|ss)  HGP      182180   0.035   0.445  0.120  0.143  0.272  0.237
  (pp|pp)  HGP       19061   0.034   0.494   0.065  0.114  0.420  0.316
```

(The tags (21|11) etc. in the raw tool output are the same classes written
la lb | lc ld.) No class is a free lunch: (ss|ss), 4 % of the cost, has no
odd component, and the bulk of the cost sits in classes with a p or d
shell at 33-49 % zero components.

### 3.3 Which zeros can be avoided where

| level | avoids | share of modelled work (naphthalene s2a) |
|---|---|---|
| digestion only | the 13 flops and 6 block updates of a zero component | 9.6 % |
| contracted level (HRR) | the HRR statements of the zero components and the dead terms of the live ones | + 1.1 % |
| primitive level (VRR, accumulation) | the odd nodes of the VRR (≈ half), the dead first-order terms of the even ones, the accumulation of nodes the HRR never reads | + 20 % |

The primitive level is where the saving is: a contracted-level or
digestion-only skip leaves the VRR, which dominates once a cc-pVDZ basis
has been decontracted to about 1.3 primitive quartets per contracted
quartet. This is also why the prototype's kernel-level speed-up (§4a)
is larger than the flop ratio: the dropped nodes are also the stores and
loads of `v`.

### 3.4 Sensitivity

The unremovable per-primitive overhead (Boys, geometry) is a guess. Naphthalene
s2a, T3 speed-up: overhead x0 (no overhead): 1.68x; x1 (default): 1.44x;
x3: 1.26x. Whatever the true value, T1 stays at 1.07-1.14x.

### 3.5 Exactness

`spike_parity verify` evaluates every quartet with `hgp_batch`, builds
G = J - K/2 from a symmetric density (not an SCF density: the identity does
not depend on it) three ways and compares with the S14.2 gate
5e-11 x max(1, max|G|):

| system | max abs G | skipped fraction (full / s2a) | largest abs v of a skipped component | max abs G(full, skip) - G(full) | max abs sym(G(s2a, skip)) - G(full) | s2a without skip |
|---|---|---|---|---|---|---|
| HOCl | 19.2 | 0.459 / 0.459 | 0 | 0 | 6.9e-14 | 6.9e-14 |
| C2H4 | 20.1 | 0.489 / 0.519 | 0 | 0 | 1.6e-13 | 1.6e-13 |
| B2H6 | 18.3 | 0.338 / 0.398 | 0 | 0 | 2.0e-13 | 2.0e-13 |
| benzene | 34.7 | 0.449 / 0.451 | 0 | 0 | 9.8e-13 | 9.8e-13 |
| H2CO cc-pVTZ (f) | 36.9 | 0.553 / 0.586 | 0 | 0 | 6.3e-13 | 6.3e-13 |

The gates are 9.6e-10 to 1.7e-9; every case is orders inside. The skip is
**bit-identical** to the unskipped G for the full plan (0 difference), and
identical to the S14.2 skeleton for s2a (the same number with and without
the skip). Both rules (stabilizer, dead coordinate) were run with the same
result.

### 3.6 The kernel prototype

`spike_parity_gen.py` emits the HGP VRR + HRR for a class and a mask with
the dead terms removed; at mask 0 the text is identical to `gen_perclass`
(`check0` asserts it). `test/spike_parity_kernel` compares the masked
kernel with the shipped one on random geometries that satisfy the mask
(6 primitive quartets per call, VRR + accumulation + HRR, no Boys or
digestion; gfortran-15 -O3, Apple Silicon, machine otherwise idle):

```
class  mask kind    max rel diff   ns/call  vs shipped
1010    4    reads   0.00E+00       11.     1.70x
1111    0    reads   0.00E+00      374.     1.15x   (shipped kernel, accumulation trimmed)
1111    4    reads   3.53E-16      150.     2.88x
1111    6    reads   3.06E-16       36.    11.91x   (centres on an axis)
2121    4    reads   4.36E-15      769.     3.28x
2222    0    reads   0.00E+00    12430.     1.10x
2222    4    reads   1.09E-14     3230.     4.25x
2222    6    reads   1.27E-15      926.    14.57x
```

All masked kernels agree with the unmasked ones to 1.1e-14 relative or
better, the three coordinates x, y, z give the same speed (1111: 2.88, 2.90,
2.93 x), and the x speed-up is larger than the flop count suggests (about
1.9x for (pp|pp)): fewer nodes means less traffic through `v`. Part of
the shipped accumulation (all nv nodes, 100 for (pp|pp)) is wasted whether
or not there is a mask: trimming it to the nodes the HRR reads is worth
1.10-1.15x by itself (the `reads` rows of mask 0), independent of this
project. Worth a separate look.

### 3.7 Host compile time of the variants

`scripts/spike_parity_compile_time.sh`: gfortran-15 -O3 -c of one emitted
kernel module, CPU seconds, while the coordinator's LMAX=4 build was running
(CPU time, not wall):

| class | mask | lines | seconds | object bytes |
|---|---|---|---|---|
| (pp\|pp) | 0 | 603 | 0.15 | 9 056 |
| (pp\|pp) | z (4) | 382 | 0.08 | 4 288 |
| (dp\|dp) | 0 | 2 368 | 1.27 | 42 296 |
| (dp\|dp) | z | 1 309 | 0.33 | 15 688 |
| (dp\|dp) | y,z (6) | 730 | 0.14 | 5 256 |
| (dd\|pp) | 0 | 3 228 | 1.74 | 53 784 |
| (dd\|pp) | z | 1 500 | 0.41 | 17 224 |
| (dd\|dd) | 0 | 19 553 | 36.3 | 329 856 |
| (dd\|dd) | z | 6 911 | 5.3 | 88 536 |
| (dd\|dd) | y,z | 2 884 | 1.2 | 21 384 |

For reference the shipped (dd|dd) file (VRR + HRR + digestion + prologue)
took 45 s in the full build's ninja log. One-axis variants cost **15-30 %**
of the unmasked kernel (the prototype omits the digestion and prologue,
which a variant repeats: ≈ 20 % of the base for (dd|dd), so ≈ 35 % in
total). Three axes, 3 x 35 %: roughly **double** the compile time of the
affected kernels. On the GPU, where the heavy classes take 15-20 minutes,
that is +5 minutes per class per axis. Register pressure falls (fewer
live nodes), not rises; this was not measured on a GPU.

## 4. Candidate schemes

### (a) Parity-masked kernels per (class, mask)

*Math*: §2. The generator drops the odd nodes, the first-order terms in the
dead directions (`PA_d*[..] + WP_d*[..]`), the matching sieve pathways, the
HRR terms with a power of AB_d or CD_d, and the odd components (`vbuf = 0`
or no digestion at all); the accumulation covers only the nodes the HRR
reads. Gill's sieve is rerun on the masked pathway set (the prototype does).

*Files*: `scripts/gen_perclass.py` (the VRR/HRR calls in `main` and
`_emit_block`), `gen_vrr.py` (masked `_pathways` and emission),
`gen_hrr.py` (term filter), the digestion skip in `emit_digestion`,
`gen_rys.py` for the Rys classes (component list and z-tables of the
chunk plan), `src/eri/trc_binkernel.F90` (dispatch key gets the mask:
`key + mask*CLASS_RADIX^4`; masks the build lacks fall back to 0).

*Gain*: §3.2 T3, 1.43-1.46x on planar D2h skeletons.

*Compile time*: §3.7, 15-30 % per axis per class.

*Skeleton*: weights unchanged. The mask is a property of the quartet's
four atoms; it does not enter the orbit weights, the plan or the
symmetrization. The two reductions are independent and multiply.

*Rys vs HGP*: both. Rys is the easier one: the component loop is already an
explicit list (chunk plan); a dead axis halves its 1D tables and drops the
odd components from the root sum. The model gives 0.36-0.40 removed on the
Rys classes, comparable to HGP. The seven heavy d classes and f classes
(where the compile cost is highest) are exactly the Rys ones. Rys chunk
plans are raced per class today; a masked variant has a different best
shape, so each variant needs its own race (a cost, not a risk).

*Risks*: kernel count x2-4; the chunk-plan races; mixed-mask launches
(§5.3).

### (b) Runtime mask inside the shipped kernels

*Math*: the same zeros, found by a launch-uniform mask word. The
generator groups the statements of each VRR level by parity vector (8
groups) and wraps each group in `if (live(g)) then ... else zero-fill`;
the HRR and digestion likewise per component.

*Gain*: TB in §3.2: 1.30-1.35x (0.25 removed vs 0.31 for (a)). It cannot
remove the dead terms *inside* live statements (the guard is per
statement), and the shipped sieve was chosen without the mask, so the
pathways it keeps are not the cheapest ones for a masked launch. Dead nodes
must hold zeros, not stale values, because the live statements still read
them (NaN-safe).

*Branch divergence*: none if the mask is uniform over the launch (a
planar molecule: every quartet is the same mask), because the guards are
on a kernel argument. For B2H6 four masks are mixed (§5.3) and a warp
would serialise up to four ways on the guarded groups; the segments
must be split by mask anyway, which is the same plumbing as (a).
Per guard cost: a uniform predicate branch, tens per primitive quartet
against hundreds to thousands of FMAs. The real cost lands on every
launch that is not masked (all non-planar molecules, and every run today): more code,
more live-range constraints through the if/else (a risk for the kernels
at 255 registers), no gain.

*Compile time*: the guard structure roughly keeps the unmasked size
(no new kernels) and adds a few percent; no multiplication by axes.

*Why consider it*: one kernel handles any axis, so no frame
canonicalisation and no variant explosion.

### (c) Basis-level reorganisation: parity blocks per shell

Split each shell into parity blocks (p: {px,py} | {pz}; d: {xx,yy,zz,xy} |
{xz,yz}) and skip zero block quartets at segment level. *This does not work
without kernel changes, and with them it is (a) in disguise*: the
kernels are generated per (la,lb,lc,ld) over whole shells, so a half-shell
has no kernel; run through the existing ones a sub-shell is a full shell
and nothing is skipped. If each of the 16 sub-quartets of a quartet (one
dead axis; 2^4) used the shipped kernel, 8 non-zero ones would each pay
the full VRR: **8x the VRR** of the unsplit class, far worse than the
baseline. Sharing one masked VRR across the blocks of a class is scheme
(a). Costs on top: the bins multiply (4 x pair-parity classes), the
Schwarz/density screening and the exchange tiles get 4x more granular,
contraction columns stay as they are. Rejected.

### (d) Better: tier it, and fix the shipped accumulation first

The measured structure suggests a staging that does not need the full
variant machinery at first:

1. *Trim the accumulation* (`g1 += v` over all nv nodes) to the nodes the HRR
   reads. Independent of symmetry, no new kernels. Kernel-level 1.10-1.15x.
2. *Masked digestion + HRR only* (T2): small compile cost, but only
   10-11 % of the work. Not worth a project.
3. *Full masked variants for the dominant classes only* (T3), the others
   unmasked: the ten most expensive classes carry ≈ 70 % of the cost
   (§3.2) and (ss|ss) 4 % has nothing to gain.

## 5. Interactions

### 5.1 The skeleton (S14.2/S14.3)

Weights unchanged (above). The skip is per evaluated quartet and the
skeleton evaluates the same quartets with fewer operations each. The
`orb` structures and the plan builder do not change. One caveat: the
representative of an orbit may have a different mask from its images only
if the orbit moves atoms off the plane, impossible for a plane that is a
group element; for a dead-coordinate rule from the geometry the masks are
those of the representative.

### 5.2 Rys, HGP, MD engines, DF, gradients

- **HGP and Rys four-centre**: §4a.
- **MD J engine** (`gen_md_j.py`, `trc_mdj`): the pair E-coefficients
  E^{ab}_t are non-zero only for t of the parity of the component when
  PA_d = 0, and the Hermite integrals R_tuv vanish for odd v when PQ_d = 0.
  With an invariant density (always so in the skeleton) the ket Hermite
  density has only even v, so the Hermite space and the R tree (the
  `md/rtree.py` sieve) halve in the dead axis. Same masks, applied to
  the Hermite index sets. Not measured. The mask has to compose with the
  S14.4 bra restriction (orthogonal: representatives of bra pairs vs a
  parity filter on Hermite indices); implement after S14.4 lands.
- **MD exchange** (`gen_md_k.py`, g classes above L=9): R per primitive
  quartet then S = PQ Hk^T; the same halving of R and of the Hk rows.
  Rarely relevant for a planar molecule at LMAX < 4.
- **DF three-centre** (`gen_df.py`): reuses the four-centre VRR with
  eta = gamma, so the masked VRR/HRR apply directly (all three centres,
  the auxiliary ones included, lie in the plane). The 2-centre metric
  (P|Q) is block-diagonal in parity; a smaller solve is possible but
  separate. Not prototyped.
- **Gradients**: derivative integrals flip the parity in the differentiated
  coordinate: dz(ab|cd) is non-zero only for *odd* combined z-parity.
  With an invariant density the contraction of the odd z-derivative with
  the (even) density blocks vanishes identically for in-plane quartets, so
  the dead-axis derivative can be dropped altogether, and the two live
  directions use only the even components. Estimated: of the three
  derivative directions one disappears and the others halve their
  components, so the saving is larger than for the energy. Not
  prototyped; needs its own measurement (`gen_rys_grad.py`).

### 5.3 Mixed masks in one bin (nonplanar molecules)

Quartets of B2H6 by dead-coordinate mask: Z=0 (no skip) 32.7 %, {x} 42.4 %,
{y} 15.5 %, {x,y} 9.2 %, {y,z} 0.1 %. In naphthalene 99.8 % are {x}, 0.2 %
{x,y} (atoms on the long axis). A launch has to be uniform in the mask, so
the plan builder splits segments by quartet mask (the bins hold pairs, the
quartet mask is the AND of two pair masks: add a per-pair mask to the
bin key, giving a bin-pair launch one mask). For planar molecules nothing
changes; for B2H6 the number of launches grows by up to the number of
masks and the gain is capped at 1.23-1.30x. Recommend: build variants for
**one** axis (the plane normal) and run everything else unmasked; B2H6
then gains only on the quartets that carry that axis (42 % of them for
x, 15 % for y), which should be an explicit non-goal.

## 6. What this does not claim

- All numbers are operation counts and CPU kernel timings, no GPU. The
  expected GPU gain follows the CPU kernel ratios only to the extent that
  the heavy classes stay compute-bound (they do: 255 registers, spills).
- The Rys cost model is analytic, not read from `gen_rys.py` text.
- The per-primitive overhead is assumed (sensitivity in §3.4).
- Gradient, DF and MD gains are analysis only.
- The density screen and exchange tiles were not applied to the counts.

## 7. Recommendation and plan

**Do (a) for one axis, in stages, after the maintainer's decision on the
frame (§8.1).** Not (b): it costs every launch, gives 4/5 of the gain, and
does not remove the need to split segments by mask. Not (c).

PRs, each stacked on the one before, each gated by "G equals unskipped G to
5e-11 x max(1, max|G|)", `symmetry=off` bit-identical, and the compile time
of the affected kernels stated:

1. **PR-A (no kernel change)**: per-atom dead-coordinate masks from
   `r_std` (exact zeros), per-pair mask in the bin key, quartet mask in the
   plan, the dispatch key extended with the mask (everything dispatches to
   mask 0 until a variant exists), stats in the segment log. Test: the
   mask table of a planar and a nonplanar molecule; planar plane
   normal canonicalised to z (if approved).
2. **PR-B (generator, HGP s/p/d)**: masked `gen_vrr`/`gen_hrr`/`gen_perclass`
   for the z axis, only for the classes that carry the cost (top ten by
   the §3.2 ranking, listed in the generator, the rest fall back). The
   trimmed accumulation as its own commit first (§4d.1), measured alone.
   Check: `check_sym_scf_skeleton` unchanged, plus a new check comparing a
   masked-launch Fock matrix with the unmasked one on benzene and a
   B2H6-like mixed case (gate above). Compile time per class stated from
   nvfortran, not gfortran.
3. **PR-C (Rys)**: masked chunk plans for the Rys classes in the same list
   (raced separately, `FOCK_CHUNK_PLAN` entries per mask).
4. **PR-D (DF 3c)** and **PR-E (gradients)**: after measuring each.
5. **PR-F (MD J/K)**: after S14.4.

Expected result on planar D2h at cc-pVDZ: 1.4x on the Fock build on top of
the skeleton, 5.0-5.6x of the ideal 8x overall.

## 8. Open questions for the maintainer

1. **Frame.** Planar D2h molecules in the standard frame come out with the
   plane normal on **x** (all atoms have x = 0: C2H4, benzene,
   naphthalene, pentacene), planar Cs ones on **z** (HOCl). One axis of
   variants needs a canonical normal. Options: rotate planar molecules so
   the normal is z (an axis permutation inside `sym_detect`; touches the
   pyscf ordering conventions of the irreps and the S12 tops), or build
   all three axes (3x the compile time, §3.7), or scheme (b).
2. How many classes and axes are we willing to build (GPU compile cost is
   not measured here; +30 % per axis per class on the host)?
3. Is a *geometric* dead-coordinate rule acceptable (it also fires with
   `symmetry=off` for a molecule whose atoms have exact zeros in the
   input frame; the integrals then differ from the unmasked ones only by
   the dropped exact zeros, i.e. not at all) or must it be tied to the
   detected group?
4. Mixed-mask bins (§5.3): split segments by mask, or only accelerate
   uniform launches?
5. The trimmed accumulation (§4d.1) is a free 10-15 % on the kernel
   independent of everything here: do it first?
6. GPU measurement: needs an A100/V100 run of one masked (dd|dd) and a
   Rys class before PR-B/C is merged (register spill, occupancy).

## 9. Reproducing

```sh
cmake --preset macos-host && cmake --build build --target spike_parity spike_parity_kernel
cd test && ../build/test/spike_parity naphthalene_d2h.xyz basis_sets/cc-pvdz.json count 1e-10 /tmp/naph.csv
../build/test/spike_parity benzene.xyz basis_sets/cc-pvdz.json both 1e-10 /tmp/benz.csv   # verify
python3 ../scripts/spike_parity_gen.py cost --lmax 3 -o /tmp/costs.json
python3 ../scripts/spike_parity_gen.py report /tmp/naph.csv /tmp/costs.json --plan s2a --top 14
../build/test/spike_parity_kernel 20000
../scripts/spike_parity_compile_time.sh ../build ~/.venvs/terco/bin/python /tmp/ct 1
```

The ctest entries `spike_parity` and `spike_parity_kernel` are
DISABLED. `git diff 2291c9c -- src` is empty.
