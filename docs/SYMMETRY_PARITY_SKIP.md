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

### 3.8 Axis-permutation round trip

Decision 1 of §8a (engine-internal axis permutation) is exact; this is the
measurement. `spike_parity ... roundtrip` (`test/spike_parity.F90`) builds J,
K and G = J - K/2 with the four-centre engine (same bins, `hgp_batch`
digestion, full and s2a plans as in §3.1) on the standard-frame coordinates.
It then rebuilds basis, pairs, bins, AO maps and orbit tables on coordinates
permuted cyclically so that the most common dead axis becomes z
((x,y,z)_new = (y,z,x)_old for a normal on x; the identity for z), maps the
symmetric density in with the induced AO permutation, runs the engine, maps
G back and, for the spherical case, folds both to spherical with the
production cart-to-sph transform. The AO permutation sends x^a y^b z^c to
the component of the same shell with the permuted powers, found with
`cart_perm_powers`, so it holds for any `TERCO_CART_ORDER`. The spherical
density is a Cartesian-symmetric density projected through T T^T. In the
permuted frame the dead-coordinate skip (rule 2 of §3.5) was also run.
cc-pVDZ, tau = 1e-10, working basis. Largest error over the J, K, G
matrices, absolute, against the gate 5e-11 max(1, max|G|) (gates 8e-10 to
2.5e-9):

| system (frame change) | space | full | s2a | full + skip | s2a + skip | gate (min) |
|---|---|---|---|---|---|---|
| HOCl (identity) | cart | 0 | 0 | 0 | 0 | 9.6e-10 |
| | sph | 0 | 0 | 0 | 0 | 8.4e-10 |
| benzene (x to z) | cart | 7.1e-14 | 7.1e-14 | 7.1e-14 | 7.1e-14 | 1.7e-09 |
| | sph | 3.2e-14 | 2.8e-14 | 3.2e-14 | 2.8e-14 | 1.5e-09 |
| naphthalene (x to z) | cart | 1.6e-13 | 6.4e-14 | 1.6e-13 | 6.4e-14 | 2.4e-09 |
| | sph | 7.1e-14 | 5.0e-14 | 7.1e-14 | 5.0e-14 | 2.2e-09 |
| B2H6 (x to z) | cart | 2.1e-14 | 1.4e-14 | (same) | (same) | 9.2e-10 |
| | sph | 8.9e-15 | 1.4e-14 | (same) | (same) | 8.5e-10 |
| H2CO cc-pVTZ, f, LMAX=3 (x to z) | cart | 1.1e-13 | 1.2e-13 | 1.1e-13 | 1.2e-13 | 1.9e-09 |
| | sph | 1.4e-13 | 1.1e-13 | 1.4e-13 | 1.1e-13 | 2.5e-09 |

J and K alone are never larger than G (K is 3e-15 to 1.2e-14, J equals G
within rounding). HOCl is bit-identical, as it must be. The other cases are
not bit-identical: the permuted frame orders bins and pairs differently, so
the accumulation order changes; the 1e-14 to 1e-13 differences are summation
noise, 4 to 5 orders inside the gate. The mapped density is invariant under
the permuted group to 2e-16 or better. With the skip in the permuted frame
the skipped fractions are those of §3.2 (0.449 benzene, 0.450 naphthalene,
0.338 / 0.398 B2H6 full / s2a, 0.553 / 0.586 H2CO) and the largest skipped
|v| is 0, so skip plus permutation equals the unskipped original frame to
the same noise. The same table was reproduced with `TERCO_CART_ORDER=gamess`
(HOCl, benzene, B2H6, all within 5e-14, log below).

**Mask counts.** The histogram of dead-coordinate masks (Z, quartets after
Schwarz) of the permuted frame is the histogram of the original with the
bits permuted, case by case, for both plans. For the planar molecules the
z-dead fraction after equals the x-dead fraction before (1.0000; the
stabilizer rule, sigma_z against sigma_x, gives 1.0000 as well; 98.5 % {z}
and 1.5 % {x,z} for benzene, 99.84 % / 0.16 % naphthalene, 73.2 % / 26.8 %
H2CO). B2H6, full plan (s2a in parentheses):

| mask | none | {x} | {y} | {x,y} | {y,z} |
|---|---|---|---|---|---|
| before | 32.74 % (24.05) | 42.45 % (41.89) | 15.54 % (17.24) | 9.21 % (16.69) | 0.06 % (0.13) |

| mask | none | {z} | {x} | {x,z} | {x,y} |
|---|---|---|---|---|---|
| after | 32.74 % (24.05) | 42.45 % (41.89) | 15.54 % (17.24) | 9.21 % (16.69) | 0.06 % (0.13) |

z now holds 42.45 % + 9.21 % = 51.7 % of the quartets (s2a: 58.6 %), the
same as x before; this is the share that runs masked under decision 1.

Logs: `/Users/albert/terco-parity-logs/lmax3.log` (H2CO f case, a helper
build at LMAX=3), `/Users/albert/terco-parity-logs/gamess.log` (gamess
order). Reproduce: `spike_parity <xyz> basis_sets/cc-pvdz.json roundtrip
1e-10 x [nskip]`.

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

## 8a. Decisions so far (2026-11-07)

1. **Frame: z-only variants, with an axis permutation inside the integral
   engine only.** The standard frame, irrep labels and every reported
   quantity stay pyscf's. In pyscf, `_adjust_planar_d2h` puts the plane
   normal on x on purpose (its comment cites the iOpenShell guide and pyscf
   issue #1201), so moving it would break label agreement.

   For the Fock build, the integral engine and its orbit tables are built on
   coordinates with the axes permuted, so that the most common dead axis is
   z. The Cartesian density goes in, and the Cartesian G comes back, through
   the induced AO permutation (x^a y^b z^c maps to another component of the
   same shell). That mapping is exact and costs O(N^2), and the fold to
   spherical happens after it as today.

   This covers every quartet of a planar molecule. On mixed-mask molecules
   it covers the most common axis (B2H6: x, so 42% + 9% of quartets); the
   rest runs unmasked. Variants for the other two axes, or scheme (b), are
   revisited only if S15 shows that mixed-mask molecules matter.

   The exactness of the permutation round trip (Cartesian and spherical) is
   shown by the prototype (§3.8): max|G_back - G| is 0 (HOCl, identity) to
   1.6e-13 (naphthalene) over full and s2a plans, J, K and G, Cartesian and
   spherical, including one f case, 4 to 5 orders inside the S14.2 gate; the
   residue is summation order, not the mapping.

Questions 2-6, decided for the fork's implementation (the user, acting as
maintainer on the fork, 2026-11-07). These are the design's own
recommendations:

2. **Scope:** z-only masked variants for the dominant classes only, i.e.
   the top ten by modelled cost (about 70% of it), with the list in the
   generator. Every other class runs unmasked. State the host
   compile-time cost; the GPU cost is measured in S15.
3. **Mask tied to the detected group.** The dead-axis mask comes from the
   stabilizer under the symmetry the SCF uses. With symmetry off or C1,
   nothing is masked, and runs stay bit-identical to today.
4. **Mixed-mask bins:** the plan builder splits segments by mask, so every
   launch is uniform in its variant.
5. **The trimmed accumulation (§4d.1) goes first,** as its own PR,
   independent of symmetry. It must be bit-identical, or justified exactly.
6. **GPU:** the masked kernels are validated on the host and the OpenMP
   build only. The GPU register/occupancy check is an S15 item, marked
   [GPU] not run.

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

## 10. Implementation (fork, decisions 1-6 of section 8a)

Commits on `sym/parity-skip`, one per PR boundary:

1. **Trimmed accumulation** (`gen_perclass.py`, `trim_accumulation`): the
   scalar kernels add only the final-buffer nodes the HRR reads. Bit-identical
   by construction (the HRR sees the same sums in the same order); the
   golden outputs agree to the last printed digit but one (scf_mgga, 3e-14
   absolute, the compiler vectorises the shorter loops differently).
2. **Plumbing** (`src/eri/trc_parity.F90`, `trc_bins.F90`, `trc_binkernel.F90`,
   `trc_eri.F90`): `parity_choose_axis` takes the axis `d` for which the group
   contains the mirror that flips `d` and every atom that mirror fixes has
   coordinate `d` exactly zero (most such shells, ties to z); a pair of two
   such shells bins apart (`BIN_G_MASKED` in the bin key), the plan carries
   `sM` per segment and sorts by (class, mask, weight), and `launch_classes`
   hands the kernels the coordinate arrays in the cyclic order that makes the
   axis their z and permutes D in and G out in `launch_in_frame` (2 nw^2
   doubles of scratch per build, only when the axis is not z). Counters on the
   engine: `parity_axis`, `nwork_masked`, `nwork_variant`, `nlaunch_variant`,
   `parity_ops_base/masked`, `parity_cost_masked/variant`. `eri%parity_mode`,
   `TRC_PARITY` and
   `parity_set_default` switch it (0 off, 1 auto, 2 force the frame).
3. **HGP variants**: 18 `pcs{tag}z` kernels (the seven HGP classes of the top
   ten, the Rys-routed ones in their bra/ket-swapped order, and the
   swaps of the HGP ones), dispatched at `key + CLASS_RADIX**4`.
4. **Rys variants**: `gen_rys.py --dead z` for (dp|pp), (dp|ds), (dp|dp),
   and (section 11) (dd|ds), (dd|pp), (dd|dp), (dd|dd).

Not done, deliberately (follow-ups):

- **DF three-centre** (`gen_df.py`) reuses the four-centre VRR and takes the
  same masks directly; the 2-centre metric is block-diagonal in parity.
- **Gradients** (`gen_rys_grad.py`): the dead-axis derivative drops out, the
  others halve their components (section 5.2); needs its own measurement.
- **MD J/K** (`gen_md_j.py`, `gen_md_k.py`): halve the Hermite index sets.
- The other 64 classes, the y and x axes without the permutation, and the
  GPU register/occupancy check of the variants ([GPU], S15): not run.
- `launch_in_frame` could use the caller's scratch instead of allocating.

## 11. The heavy classes (`sym/parity-heavy-classes`)

Section 10's variant set is the ten classes of decision 2. What it actually
covers was then measured on the shipped engine, and it is not what the work
counters suggested.

**Work items are the wrong unit.** `nwork_variant / nwork_masked` reads
67.7 % (benzene), 69.9 % (naphthalene_d2h) and 72.0 % (pentacene_d2h) --
6-31G*, `--symmetry auto`, skeleton plan, host build. Weighting each masked
launch by `pc_cost`, the operations one of its work items costs, the covered
share is **47.2 / 49.8 / 54.4 %**: a work item of (ss|ss) is 2 modelled
operations and one of (dd|dd) is 63 644. Half the masked work was reaching no
variant, and it was the expensive half. Shares of the masked modelled cost
with no variant, benzene / naphthalene / pentacene:

| class | % | class | % |
|---|---|---|---|
| (dp\|dd) | 9.3 / 11.5 / 9.9 | (dd\|dd) | 10.1 / 8.3 / 5.9 |
| (ds\|dd) | 6.2 / 5.5 / 4.3 | (ps\|dd) | 5.3 / 5.0 / 4.3 |
| (pp\|dd) | 4.6 / 4.6 / 4.2 | (ds\|ds) | 4.5 / 3.7 / 3.1 |
| (dd\|pp) | 3.9 / 2.4 / 2.4 | (pp\|pp) | 2.4 / 2.4 / 3.0 |

Every uncovered class at the top has a d shell on **both** sides, and the
largest single one is (dp|dd) -- the bra/ket-swapped ordering of (dd|dp).
`bins_skel_cut_first` cuts the ket side of those bin pairs, so production
launches them in the swapped order, which `rys_default_class` does not route
to Rys: they run HGP, where a masked variant is the cheaper kernel.

**What was added.** (dd|ds), (dd|pp), (dd|dp) and (dd|dd) joined
`PARITY_RYS_TOP`, which gives four `gen_rys.py --dead z` chunked kernels and,
through the swap rule, three HGP ones -- (ps|dd), (pp|dd), (dp|dd). The
cost-weighted covered share goes 0.472 -> 0.807 (benzene) and
0.498 -> 0.829 (naphthalene), 0.544 -> 0.842 (pentacene); the work-item share
only 0.677 -> 0.694, 0.699 -> 0.718 and 0.720 -> 0.739, which is why the
cost-weighted pair of counters exists and why `check_sym_parity` now gates on
it per molecule.

**What it is worth, measured, and why it is not what the coverage says.**
Host build, gfortran 15.1 `-O3` (`do concurrent` serial), `terco <mol>.xyz -b
basis_sets/6-31gs.json --symmetry auto`, the timing block's `fock` column,
best of 2-3 runs, `TRC_PARITY=0` over `TRC_PARITY=1`:

| system | before | after |
|---|---|---|
| benzene | 2.120 / 1.804 = **1.175x** | 2.097 / 1.748 = **1.200x** |
| naphthalene_d2h | 12.264 / 10.083 = **1.216x** | 11.885 / 9.591 = **1.239x** |
| pentacene_d2h | (not usable, below) | 97.806 / 79.734 = **1.227x** |

Two other agents were building and running on the machine throughout, so
pentacene's before pair (97.304 / 99.121 = 1.019x) saw load its after pair
did not -- it disagrees with the 1.16-1.19x this configuration is known to
give, so only its after value is quoted. Energies are unchanged to 3e-12 Eh
against a 1e-10 convergence (benzene `...128250` -> `...128247` with parity
on; parity off is bit-identical), and `check_sym_parity` puts masked against
unmasked G at 8.7e-16 to 6.2e-12 on a 5e-11 gate.

**+2 % of Fock time for 33 points of modelled coverage.** `pc_cost` badly
overstates what these classes cost in time: they have few quartets, and
production runs the (dd|..) ones through chunked Rys, whose cost is nothing
like HGP's statement count. The counter pair is the right instrument for
finding which classes have no variant and the wrong one for predicting what
adding them buys. It also says the shortfall against section 1's
1.43-1.46x is **not** mainly class coverage: on the host the ten mid classes
that already had variants were carrying nearly all of the available gain,
and the rest of the gap is the model's per-primitive overhead term
(section 3.4: x3 gives 1.26x) plus the host's own cost profile. The GPU does
agree that the heavy classes carry much more of the time -- 11.4 to 19.1 %
of it, against the ~2 % of host Fock time this stage buys -- and section 12
has the measurement.

**Compile cost** (AGENTS.md: a generator change that grows the heavy classes
states it). nvfortran 26.5, `-stdpar=gpu -acc -gpu=cc70,mem:separate -O3`,
one object, seconds, the kernel a variant shadows -> the variant:

| class | shadowed | variant | | class | shadowed | variant |
|---|---|---|---|---|---|---|
| (ps\|dd) HGP | 11.2 | 7.5 | | (dd\|ds) Rys | 6.5 | 8.8 |
| (pp\|dd) HGP | 42.5 | 19.6 | | (dd\|pp) Rys | 7.3 | 12.1 |
| (dp\|dd) HGP | 157.5 | 65.6 | | (dd\|dp) Rys | 4.7 | 10.3 |
| | | | | (dd\|dd) Rys | 6.6 | 10.8 |

The variants are extra objects, so the set grows by their own **135 s** --
against 547 s for the HGP (dd|dd) kernel alone. A masked HGP kernel is
smaller than the one it shadows (the sieve reruns on half the nodes:
(dp|dd) 7267 lines -> 3327), a masked chunked Rys kernel larger (the chunk
body exists once per parity of the outer components' z powers, each half the
size), which is why the four cheap Rys classes are the ones that cost more
than they replace. gfortran 15.1 `-O3`, the same seven objects inside the
host build: the three HGP variants 2.8 / 7.0 / 16.1 s and the four Rys ones
1.8-2.7 s, **34.9 s** in total, on a set whose per-object times sum to
2606 s and whose own HGP (dd|dd) kernel is 123 s.

**Registers.** cc70, read from a linked shared library with
`cuobjdump -res-usage` (never the `.o`), `STACK` not `LOCAL`. Every chunked
Rys kernel here already sits at 255 registers with a non-zero frame, before
and after; the variants do not change the kernels they shadow, and their own
frames are:

| class | unmasked | variant |
|---|---|---|
| (dd\|ds) | 192 B | 144 B |
| (dd\|pp) | 384 B | 464 B |
| (dd\|dp) | 288 B | 360 B |
| (dd\|dd) | 360 B | 568 B |

For comparison the variants section 10 already shipped: (dp|pp) 296 -> 376 B,
(dp|ds) 112 -> 72 B, (dp|dp) 216 -> 224 B. The new Rys ones are the same shape
of change, and the largest, (dd|dd) at 568 B, is inside the 56-1256 B the
chunked kernels already run with (docs/GENERATORS.md).

**The three HGP variants are the risk, and one of them is a regression.**
Same measurement, the scalar HGP kernels:

| class | unmasked | variant |
|---|---|---|
| (ps\|dd) | 2 192 B | 1 824 B |
| (pp\|dd) | 9 808 B | 5 040 B |
| (dp\|dd) | 15 904 B | **25 472 B** |

A masked HGP kernel does about a third of the arithmetic but declares the
SAME scratch -- `v(700, 0:1)`, `g1(700)`, `vbuf(648)`, 22 kB of doubles --
because the generator only stops writing the dead entries, it does not
compact the buffers the way the warp-role path does with its `vmap`. On
(ps|dd) and (pp|dd) the compiler still eliminates more of that scratch than
before; on (dp|dd) it eliminates less, and the frame grows by 9.6 kB per
thread on a kernel that already spilled 15.9 kB. docs/GENERATORS.md's
precedent is the array-form Rys kernels, which ran at 3.5 % of SM throughput
with 7.5 kB. **So (dp|dd) must not be believed until a GPU has run it**, even
though it is the largest single uncovered class; compacting the masked
kernel's buffers to the live node set is the generator work that would fix it
properly, and holding the (dp|dd) swap out of `parity_hgp_tags` is the
one-line alternative.

**The above was measured in a device link of these kernels and their module
dependencies only, not of `libterco.so`**, because this stage had no GPU
toolchain run over the whole tree. Section 12 closes that: a full device
build was made afterwards and the full-library registers, the occupancy and
the run-time effect are measured there.

**Not added, and why.** (ds|dd) 4.3-6.2 % and (ss|dd) 1.5 % would need
(dd|ds) and (dd|ss) in `PARITY_HGP_TOP`, and (dd|ds) is a `ROLE_SPLITS`
class: a masked launch would leave its warp-role kernel for a scalar one, and
the split is worth 1.1-4.5x exactly where the scalar working set spills.
Nothing on a host build can measure that, so a masked **role** kernel is the
next piece of generator work, not a list edit. (ds|ds) 3.1-4.5 % and
(pp|pp) 2.4-3.0 % are the next two uncovered classes without that
complication.

**Still [GPU].** The chunk shape of a masked Rys kernel is the unmasked
class's shape; the design (section 4a) says a masked variant has its own best
shape and needs its own race, which needs a device. So does the occupancy and
spill effect of the 568 B frame, and so does whether the HGP masked (dp|dd)
kernel beats the Rys route the canonical ordering takes.

## 12. On a GPU (measured after section 11)

Full device build of this branch: nvfortran 26.5, `-stdpar=gpu -acc
-gpu=cc70,mem:separate`, `TERCO_LMAX=3`, `TERCO_DF_LMAX=2`, Tesla
V100-DGXS-32GB (cc70), 33.8 min at `--parallel 12`. Baseline is the same
tree without these variants (`sym/w2-segment-weight`), identically
configured, so the W2 launch fusion is in both. ERI figures are summed over
`trc_pc_*` from `nsys cuda_gpu_kern_sum` for one steady-state Fock build (the
difference between two `--iterations` runs); 6-31G* is three repeats, per
class reproducible to under 0.5 %, cc-pVDZ single runs. `--symmetry auto`,
`--econv 1e-10`, `TRC_PARITY=1` against `TRC_PARITY=0`.

**The headline is not the increment. Before this stage, parity skipping was a
net loss on the GPU** in five of six configurations -- off/on of 0.964,
1.000, 1.005, 0.967, 0.964, 0.964 -- i.e. the shipped default cost time on a
device. With these variants it is 0.991 to 1.059. That is the case for the
stage, and it was invisible to section 11 because the host cost profile is
not the device's.

| system / basis | ERI gain | this stage's own effect on Fock |
|---|---|---|
| benzene 6-31G* | +3.93 % | +2.68 % |
| naphthalene 6-31G* | +5.46 % | +2.96 % |
| pentacene 6-31G* | +7.10 % | +4.88 % |
| benzene cc-pVDZ | +4.92 % | +3.30 % |
| naphthalene cc-pVDZ | +1.37 % | +2.59 % |
| pentacene cc-pVDZ | +2.20 % | +4.84 % |

All seven variants launch, 25 variant launches per build, and all 25 masked
kernels together are about 31 % of ERI time. Coverage reproduces section 11's
host numbers exactly (0.807 / 0.829 / 0.842).

Per class, masked against unmasked (6-31G*, pentacene):

| | gain |
|---|---|
| `pcs1122z` (pp\|dd) | **+58.5 %** |
| `pcs1022z` (ps\|dd) | +43.1 % |
| `rys2210z` (dd\|ds) | +17.5 % |
| `rys2221z` (dd\|dp) | +9.1 % |
| `rys2211z` (dd\|pp) | +5.6 % |
| `rys2222z` (dd\|dd) | -2.2 % |
| `pcs2122z` (dp\|dd) | **-18.2 % benzene, -9.7 % naphthalene, +4.1 % pentacene** |

The seven classes carry 11.4 / 14.0 / 19.1 % of ERI-off time, and the
realised saving *within* them is +3.9 / +13.2 / +24.8 % -- so 4 to 25 % of
the available ceiling, not the half a naive reading would give. The entire
shortfall is `(dp|dd)`.

**The `(dp|dd)` frame regression is real but costs no occupancy.**
`cuobjdump -res-usage` on both fully linked libraries: 0 of the common
functions changed `REG` or `STACK`, none went from `STACK` 0 to non-zero, the
maximum is unchanged at 169 600 B, and the total `STACK` delta is exactly the
sum of the seven new frames (33 872 B). `ncu` on `pcs2122z` against
`pcs2122`, identical grids: registers 255 both, the same four occupancy
limiters, theoretical occupancy 12.50 % both, achieved 9.1 against 9.2 %,
duration +0.1 %. So the 25 472 B is not an occupancy cost -- but it is not
merely declared either: the masked kernel executes **6.8 % fewer
instructions** and spends it on **15.5 % more local loads and 13.3 % more
local stores**, 39.5 % of all its instructions against 32.2 %. On a small
grid (benzene's `pcs2122` launches three blocks) there is no parallelism to
hide that latency, which is why it inverts to a loss there.

`STACK` mixes declared scratch with spill slots and cannot be decomposed from
a linked library. `nvfortran -gpu=...,ptxinfo` can, and it shows the frame
growing while the spill traffic falls:

| kernel | frame | spill stores | spill loads |
|---|---|---|---|
| `pcs2122` | 15 904 B | 23 344 B | 29 820 B |
| `pcs2122z` | 25 472 B | 15 532 B | 12 632 B |

So the claim that spills are zero anywhere in this library is false, and
`-res-usage`'s `LOCAL` field (0 throughout) does not say otherwise.

**The `(dp|dd)` retreat is not worth taking.** Dropping it changes the total
gain by +0.35 ms summed over all six configurations, and the two losses land
on pentacene, the most production-relevant case, at 0.19 to 0.29 points --
below the 0.5 to 0.7 % wall-clock trial spread. It would trade the coverage
gate for no measurable speed. What it would buy is 65.6 s of the stage's 135 s
of compile and 3 327 of its 6 457 generated lines, so it is defensible on
compile cost alone.

**The better follow-up** is the compaction section 10 already names: compact
the masked kernel's `v`/`g1`/`vbuf` to the live node set as `role_plan`'s
`vmap` does. `pcs2122z` already executes fewer instructions and loses it to
carrying the full uncompacted scratch; that would turn `(dp|dd)` from
break-even into a gain and likely lift `rys2222z` and `rys2211z` too.

Device suite: 276 of 277, the one failure being `terco/link_c`, a known
build-configuration defect fixed on `sym/capi-acc-gpu` and unrelated to this
stage. `check_sym_parity` and `check_sym_parity_f` both pass on the device.
Energies across all six systems, both builds, both parity settings: maximum
deviation 1.02e-12 Eh, under this device's own run-to-run noise (up to
1.4e-12, growing with system size), against a 1e-10 convergence. Iteration
counts identical everywhere.
