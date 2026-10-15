<!--
Kept verbatim from libfint (/home/jorge/dev/libfint/README_TERCO.md), written by
the libfint session on 2026-09-27 as a guide for porting its generators to
terco. It was taken out of libfint's git history on 2026-09-30, so this is
now the durable copy. Hashes and paths in it refer to libfint as it was then;
scripts/rotaxis/VENDORED.md maps the rewritten hashes. Some advice here was
later measured and overturned in terco: see the rotaxis commits (e2b401a,
ebfbb6a) for the GPU numbers, and note that GAMESS also needs a
per-component rescale, not only a reorder (src/core/trc_ao_norm.F90).
-->

# libfint's kernel generators — a guide for the terco session

You are adapting libfint's integral-kernel generators so that **terco** (the GPU,
`do concurrent` SCF engine) can evaluate ERIs by Rys quadrature and by rotated-axis
McMurchie–Davidson, matching what `$HOME/dev/gamess-libERI` does today and
eventually replacing it.

This document explains how libfint's generators are built, where the output
conventions are decided, and what will bite you. It assumes you have not read
libfint's source.

Cross-references: `doc/ROT_AXIS_MMD.md` (the derivation in prose), `doc/HGP.md`
(the second algebraic path, including §9/§9b on the crossover), `doc/PORTING.md`
(the Rys/Lisp side), `doc/PERFORMANCE.md`.

---

## 0. The single most important thing

**libfint's generators derive the algebra. gamess-libERI's does not.**

`$HOME/dev/gamess-libERI/rot_axis_generator/kernel_formulas.py` loads its formulas
by *parsing the existing hand-written `.F90` kernels* (`extract_sections.py`, a
24-section split: `header`, `declarations`, `axis_rotation`, `boys_function`,
`r_integrals`, `eri_assembly`, …) and re-assembling them. It is a scaffolding
generator. It cannot change a formula, reorder a component, or retarget a
language, because it never had the mathematics — only the text.

libfint's `scripts/rotaxis_mmd/` starts from the Hermite expansion and derives
every kernel symbolically, with exact rational arithmetic and no sympy
dependency. Nothing is hardcoded. That is the whole reason this is worth adopting:

> To emit GAMESS component ordering you change **one function**. To emit CUDA
> Fortran, OpenMP target, or `do concurrent` you change **one file**. Everything
> downstream re-derives.

Keep that property. It is easy to lose by hand-editing generated files.

---

## 1. The two generator families

libfint has two, and they are completely different in kind.

| | rotated-axis MMD (and HGP) | Rys / general path |
|---|---|---|
| Language | Python: 1,288 lines (rotaxis) + 1,165 (HGP) | Common Lisp: 1,732 lines of libfint's backend on 1,093 of upstream's | 
| Source of truth | symbolic derivation from scratch | upstream libcint's `parser.cl` + `derivator.cl` |
| Lives in | `scripts/rotaxis_mmd/`, `scripts/hgp/` | `scripts/gen_f90_all.cl`, `scripts/f90-backend.cl` |
| Emits | `src/rotaxis/`, `src/hgp/` + dispatchers | `src/cint_gen_*.f90` |
| Covers | s, p, L energies + gradients; HGP to d | every integral libcint declares, all l |
| Owned by | libfint | shared with upstream |

You almost certainly want the Python family for terco's hot kernels and the Lisp
family only as the correctness reference. Read §2 first, §3 when you need Rys.

---

## 2. The rotated-axis MMD generator

### 2.1 Pipeline

```
cartorder.py the Cartesian component order and its inverse, one definition,
   |         shared with scripts/hgp                        <-- ordering lives here
   v
poly.py      sparse multivariate polynomials, Fraction coefficients
   |         (136 lines; deliberately no sympy -- libfint has no build deps)
   v
derive.py    the frame, the Hermite E coefficients, the R^n_tuv recursion,
   |         and kind_components(kind, order)
   v
factor.py    group every term into three loop levels
   |
   v
emit.py      Fortran text: kernel bodies, the dispatcher, kernels.cmake
   |         (the ONLY file that knows what language we target)
   v
generate.py  CLI; picks the class list; writes files; deletes stale ones
```

`cartorder.py` + `poly.py` + `derive.py` + `factor.py` are **language- and
target-agnostic** (507 lines plus cartorder). `emit.py` is 689 lines and is entirely Fortran text. That
split is your port surface.

### 2.2 The frame (`derive.py` docstring)

```
origin at A, z along A->B, y along AB x CD, x = y x z

A = (0,0,0)        B = (0,0,R_AB)
C = (cx,cy,cz)     D = C + (Rs,0,Rc)      Rs = R_CD sin g, Rc = R_CD cos g
```

For a bra primitive pair, `p = a+b` and the product centre sits on z:
`P = (0,0,za)` with `za = b R_AB / p`, and `zb = za - R_AB`. So
`X_PA = X_PB = Y_PA = Y_PB = 0`, and for the ket `Y_QC = Y_QD = 0`.

That is the entire point of the rotation: the x and y Hermite expansions carry
only the `1/(2p)`, `1/(2q)` ladder terms, and **every geometry-dependent quantity
of the bra pair is a scalar** (`za`, `zb`) rather than a vector.

Symbols, grouped by what they depend on — this grouping *is* the factorisation:

| group | symbols |
|---|---|
| bra primitive | `ip = 1/(2p)`, `za`, `zb`, `zq = Z_PQ`, `ea` (gradients only), `B0..BL` |
| ket primitive | `iq = 1/(2q)`, `yc`, `yd`, `ax` |
| quartet only | `Rs`, `Rc`, `cy` |

with `B_n = (-2 alpha)^n F_n(T) / sqrt(p+q)`.

The expansion is Helgaker–Jørgensen–Olsen eq. 9.9.33 with recursions 9.9.18–20.
`derive.py:96 E()` is the Hermite coefficient recursion per Cartesian direction;
`derive.py:135 R()` is the `R^n_tuv` recursion; `derive.py:178 integral()`
multiplies the three direction sums and contracts with `R`.

### 2.3 Kinds, not angular momenta (`derive.py:59`)

A quartet slot is a **KIND**, not an `l`. The important case is `L`: libcint's
`KAPPA_SP_SHELL`, an s and a p on shared exponents with a separate coefficient
column for each. It is **one slot of four components**, not two slots.

```python
KINDS = {
  "s": [((0,0,0), 0)],
  "p": [(c, 0) for c in cart_components(1)],
  "d": [(c, 0) for c in cart_components(2)],
  "f": [(c, 0) for c in cart_components(3)],
  "L": [((0,0,0), 0)] + [(c, 1) for c in cart_components(1)],
}
KIND_LMAX  = {"s":0, "p":1, "d":2, "f":3, "L":1}
KIND_NTYPES= {"s":1, "p":1, "d":1, "f":1, "L":2}   # coefficient column blocks
KIND_RANK  = {"s":0, "p":1, "L":2, "d":3, "f":4}   # canonical sort order
```

The second element of each pair is the **coefficient type** — which column block
of the contraction coefficients that component draws from. An L shell's s
component takes column block 0 and its p components block 1. Getting this wrong
is the single most common bug; see §5.

### 2.4 The three loop levels (`factor.py`)

Every term of every component is a product of
`(bra monomial) x (ket monomial) x (quartet constant)`, and the kernel is two
nested contraction loops, ket outside, bra inside. So:

**Level 1, inside the bra loop** — one accumulator per distinct bra monomial:
```
S[g,n,k] += kab * g(ip,za,zb) * B_n * zq^k
```
Nothing here depends on which Cartesian component is being built. That is why
the inner loop of a `(pp|pp)` class is a few dozen multiply-adds rather than a
few hundred, and it is where rotated-axis wins.

**Level 2, once per ket primitive after the bra loop** — one accumulator per
distinct ket polynomial, shared across components that need the same combination:
```
r[j] += kcd * P_j(iq, yc, yd, ax; S)
```

**Level 3, once per quartet**:
```
comp = sum  scale * (Rs^a Rc^b cy^c) * r[j]
```

The factoriser is exact — it groups monomials, it never approximates. `scale` is
a `Fraction`, in practice always an integer.

`Factorised.summary()` prints the sizes; run with `--summary` to see them without
writing files. Those counts are the loop-level arithmetic used in `doc/HGP.md` §9
to predict the rotaxis/HGP crossover at ~1.5 primitive pairs, which was then
confirmed by measurement (§9b).

### 2.5 The emitted kernel ABI (`emit.py:103`)

```fortran
subroutine rotaxis_LLLL(nbra, ncb, bp, kab, nket, nck, kp, kcd, gc, cutoff, res, any)
   real(dp), intent(in)  :: bp(5, nbra), kab(NTT*ncb, nbra)
   real(dp), intent(in)  :: kp(7, nket), kcd(NTT*nck, nket), gc(3), cutoff
   real(dp), intent(out) :: res(NCOMP, ncb*nck)
   logical,  intent(out) :: any
```

| arg | contents |
|---|---|
| `bp(5, nbra)` | bra pair table: `p`, `1/(2p)`, `za`, `zb`, `mu R_AB^2`. **Six rows for gradients**, row 6 = the bra exponent `a` |
| `kab(NTT*ncb, nbra)` | bra pair weight `c_i c_j exp(-mu R_AB^2) / p`, at `tt + NTT*(cc-1) + 1`; `NTT` = product of `ntype` over the two bra slots |
| `kp(7, nket)` | ket pair table: `q`, `1/(2q)`, `c/q`, `d/q`, `ax`, `az`, `nu R_CD^2` |
| `kcd` | ket pair weight, laid out as `kab` |
| `gc(3)` | `Rs`, `Rc`, `cy` |
| `res(NCOMP, ncb*nck)` | local-frame Cartesian integrals, column `(ccb-1) + ncb*(cck-1) + 1`, **without** the `2 pi^(5/2)` and the per-l normalisation — the driver applies those |

The `1/sqrt(p+q)` is folded into `B_n` per primitive quartet; with the `1/p` in
`kab` and the `1/q` in `kcd` that is the full `1/(p q sqrt(p+q))`.

Two emission modes:

- **Unrolled** (`emit_kernel`) — straight-line statements. Default below
  `UNROLL_LIMIT = 3000` level-2 terms (`emit.py:270`).
- **Table-driven** (`emit_kernel_tabled`) — the same DAG as data arrays plus
  loops. A `(dd|dd)` class has ~96k level-2 terms and its unrolled form takes
  gfortran minutes.

Long procedures are chunked into contained procedures of ≤250 statements. One
8,584-statement procedure took gfortran 577 s; chunked it took 27 s, and AOCC
stopped reporting a bogus "killed by user".

**This is a decision you should probably reverse.** libfint backed off full
unrolling because it targets CPUs — GAMESS's loop-driven kernels are the right
shape there. terco is a GPU engine, and gamess-libERI's fully-unrolled kernels
are unrolled precisely because that is what GPUs want. Raise `--unroll-limit`
hard, or drop the table mode, and measure. The generator supports either; only
the flag changes.

### 2.6 `generate.py` — the CLI

```bash
python3 scripts/rotaxis_mmd/generate.py [--lmax N] [--no-L] [--grad] \
        [--summary] [--unroll-limit N] [--max-terms N] [--root DIR]
```

**Energy classes are canonicalised**: kinds ascend within each pair and pairs
ascend, by `KIND_RANK`. The driver permutes an incoming quartet into canonical
order, which is legal because the integral is symmetric under those swaps.

**Gradient classes are NOT canonicalised.** `d/dA` names one centre, so whichever
shell the caller put in slot 0 has to stay there. Every *ordered* combination of
kinds gets its own kernel — 81 rather than 21 for three kinds. This bit us once;
see §5.

The exact commands that reproduce the committed tree byte-for-byte:

```bash
python3 scripts/rotaxis_mmd/generate.py --lmax 1            # src/rotaxis/      21 kernels
python3 scripts/rotaxis_mmd/generate.py --lmax 1 --grad     # src/rotaxis_grad/ 81 kernels
python3 scripts/hgp/generate.py        --lmax 2             # src/hgp/          55 kernels
python3 scripts/hgp/generate.py        --lmax 1 --grad      # src/hgp_grad/     81 kernels
```

Note `--lmax 1` for rotaxis, not the default `--lmax 2`: d was deliberately cut
from the rotated-axis energy set because d is catastrophic for this factorisation
(a `dddd` gradient is 538,433 terms against HGP's 28,696) and Rys is simply
better there. `--lmax 2` will happily generate d kernels; they are correct and
slow.

Each run also writes the dispatcher (`src/cint_rotaxis_kernels.f90`, a `select
case` on `class_code`) and `src/rotaxis/kernels.cmake` (the file list).
`src/rotaxis/cint_rotaxis_boys.f90` is generated too -- it carries the Boys
function for these kernels, with the stability branches described in section 5.
Read its header comment before porting it; it explains the failure it works
around more precisely than section 5 does. Stale
class files are deleted, so a narrower `--lmax` cleans up after a wider one.

### 2.7 HGP, the sibling generator

`scripts/hgp/` is the same shape but emits a recurrence DAG rather than an
expanded polynomial: `recur.py` has `Graph.vrr`, `hrr_bra`, `hrr_ket`, `grad`,
and each intermediate is a node computed once. Read its module docstring — it
explains why expanding `[dd|dd]` into 96k-term tables lost to Rys, which is the
central lesson about when to expand and when to recur.

`scripts/hgp/check_numeric.py` is a pure-Python numerical check of the DAG,
independent of Fortran: `check()` compares it with McMurchie-Davidson,
`check_grad()` with a central finite difference. Run it as
`python3 -m hgp.check_numeric` from `scripts/`; it takes about a second and
exits non-zero on failure.

It was broken for a long time and the reason is worth knowing, because the
same trap is waiting in any evaluator you write. Two node kinds carry no
expression in the graph — `('v', e, f, m)`, a genuine Boys base whose fourth
element is the Boys order, and `('c', a, f, buf)`, a *contracted vertical
result* whose fourth element is a buffer index. `evaluate()` treated both as
bases and computed `pref * boys(k[3], t)`, so every `('c', …, 0)` node came
back as `boys(0, t)` — the same number for every component of a class. It
agreed with the reference at (ss|ss), where the buffer index coincides with
the right `m`, and nowhere else, which is exactly the shape that makes a bug
look like an angular-momentum problem. `check()` also ordered over the targets
alone, ignoring the `vroots` that `build()` returns, so the vertical nodes were
never evaluated at all. The emitter always did this correctly — it reads the
value out of the contracted buffer by `tindex` — so the kernels were never
affected.

Both are fixed, and the file now runs in CI (`.github/workflows/generators-ci.yml`).
Build the equivalent for anything you add, and wire yours into CI from the
start: nothing ran this file, which is the only reason it could sit broken.

---

## 3. The Rys generator

Rys is **not** libfint's own derivation. It is upstream libcint's autocode with a
Fortran backend bolted on.

```
scripts/gen_f90_all.cl     276 lines   libfint's: which integrals, in F90
scripts/f90-backend.cl   1,456 lines   libfint's: the Fortran emitter
   uses upstream's:
      parser.cl            684 lines   symbolic integral expression parsing   } 1,093 lines,
      derivator.cl         409 lines   symbolic differentiation               } language-agnostic
      auto_intor.cl        247 lines   the integral list, in a DSL
```

`parser.cl` + `derivator.cl` are **language-agnostic** — commutation rules, the
σ·p algebra, Gaunt and Breit gauge terms, GIAO factors. A port does not touch
them. Integrals are declared in a DSL:

```lisp
'("int1e_kin"      (.5 \| p dot p))
'("int2e_ip1"      (nabla \, \| r12 \| \,))
'("int1e_spsp"     (sigma dot p \| sigma dot p))
```

and `f90-backend.cl` emits `_cart`, `_sph` and `_spinor` from one description.
`scripts/fetch_libcint.sh` pulls upstream (unpinned, deliberately — the
bit-identity suite is a drift detector and a pin would suppress the news).

### What is NOT generated, and matters to you

Four Python scripts transcribe C that no generator produces:

| script | target | why |
|---|---|---|
| `c_rysroot_to_fortran.py` | `cint_rys_root_n.f90` | 1,376 lines of Chebyshev-like polynomial fits for 1–5 Rys roots/weights. No algorithm, only arithmetic |
| `c_g2e_unrolled_to_fortran.py` | `cint_g2e_unrolled.f90` | ~70 `_g0_2d4d_XXXX` routines, **hand-written in upstream libcint** — no script produces them |
| `c_tables_to_fortran.py` | static numeric tables | transforms the committed C rather than re-running the mpmath that made it |
| `c_cart2sph_to_fortran.py` | `cint_cart2sph.f90` | the unrolled l=0..4 transforms; they fix the summation order, which is what makes spherical integrals agree to the last bit |

The second row is the one to plan around. `CINTg0_2e_2d4d_unrolled` has **no
fallback** and is selected whenever the Rys order is ≤ 2 — most of chemistry. If
terco needs Rys on the GPU you must either port those ~70 routines or write the
generator upstream never wrote. Budget for it; do not discover it late.

The Rys algorithm proper lives in `src/cint_g2e.f90`, `cint_g2e_ops.f90`,
`cint_rys_roots.f90`, `cint_find_roots.f90`, `cint_tab_roots_x0.f90`.

---

## 4. Changing the component ordering to GAMESS

This is your main task and it is genuinely small.

### 4.1 What libfint emits

`scripts/cartorder.py` is now the single definition, shared by both generators:

```python
ORDERS = {"lexicographic": _lexicographic, "gamess": _gamess,
          "libcint": _lexicographic}        # alias, same function

@lru_cache(maxsize=None)
def components(l, order="lexicographic"): ...      # fixes the layout
def index(lx, ly, lz, order="lexicographic"): ...  # its inverse, BY LOOKUP
```

The default is called **`lexicographic`**, not `libcint`: it is the component
labels in alphabetical order with x < y < z (XX XY XZ YY YZ ZZ), equivalently
the exponent triples in *descending* lexicographic order, and libint, Psi4,
PySCF and the CCA standard all share it. Mind the direction — *ascending*
lexicographic on the triples is the exact reverse, ZZ YZ YY XZ XY XX. `libcint`
remains as an alias resolving to the same function.

`index` is derived from `components`, so the two cannot disagree. They used to
be written independently -- a comprehension in two copies and a closed-form
piece of arithmetic -- and the arithmetic is correct only for the lexicographic order:
reordering the list alone misplaces 5 of 6 d components and 9 of 10 f
components, with no error raised. `scripts/hgp/cart.py` is now a thin shim over
this module.

| l | lexicographic | GAMESS |
|---|---|---|
| 0 | S | S |
| 1 | X Y Z | X Y Z — **identical** |
| 2 | XX XY XZ YY YZ ZZ | XX YY ZZ XY XZ YZ |
| 3 | XXX XXY XXZ XYY XYZ XZZ YYY YYZ YZZ ZZZ | XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ |
| 4 | *(the same rule)* | XXXX YYYY ZZZZ XXXY XXXZ XYYY YYYZ XZZZ YZZZ XXYY XXZZ YYZZ XXYZ XYYZ XYZZ |

GAMESS's d order was verified empirically, not assumed: in
`gamess-libERI/rhf/rot_axis/int0022.F90` the 6×6 block's diagonal lands at
`eri_value(1, 8, 15)`; `eri_value(1)` carries `qx**4`, `eri_value(8)` carries no
geometry at all (the y direction, which has none in this frame), and
`eri_value(15)` carries `qz**4`. So 1,2,3 = XX,YY,ZZ.

The permutation lexicographic→GAMESS for d is `[0, 3, 5, 1, 2, 4]`.

f and g are now transcribed in `_GAMESS_TABLE`, from GAMESS's own
`int2a.src` (subroutine `SHELLS`, the `DATA LX/LY/LZ` tables, flat indices
11–20 for f and 21–35 for g) read across all three arrays — `LX` alone does
not determine the triple wherever it is 0. f is corroborated independently by
`gamess-libERI/rhf/rys/int3000_rysgen.F90` lines 69–116, whose `ix/iy/iz`
offsets minus one reproduce it exactly; that same table's indices 5–10 are a
second source for the d order above. g rests on `int2a.src` alone, which is
GAMESS itself rather than a derived copy, but has not been cross-checked
against a second implementation. l >= 5 still raises rather than guessing.

### 4.2 Where the order is decided

Three places, in order of preference for making the change:

1. **Pass `--order gamess`.** Both generators take it
   (`scripts/rotaxis_mmd/generate.py`, `scripts/hgp/generate.py`), threaded as
   an explicit parameter defaulting to `"lexicographic"`: `Derivation(..., order=)`,
   `build(..., order=)`, `ClassPlan(..., order=)`, `kind_components(kind, order)`.
   Never module state — a build-wide switch across several generators would
   otherwise cross-contaminate the first time one process emitted two orders,
   and the failure would look like a miscompiled kernel rather than a stale
   global. **This is the right place**, and `gamess` is already defined for
   l = 0..2.

   One distinction cost me a rewrite, and you will meet it in your own
   generators. `cart_components` serves two different purposes. The **output**
   layout — kind components, the transfer roots, `cart_index` — follows the
   requested order. The **internal** vertical targets must not: they are
   addressed through `tindex`, their sequence is invisible outside the kernel,
   and their l runs *above* the output l, because a gradient p-class needs f
   intermediates. An output order that stops at d cannot enumerate those, so
   parameterising the internal enumeration makes `--order gamess` die on an
   ordinary p gradient class.

2. **`derive.py all_components()`** (and `all_grad_components()`) — the
   loop nest that fixes the output index order across the four slots:

   ```python
   for cd, td in kd:            # l slowest
       for cc, tc in kc:
           for cb, tb in kb:
               for ca, ta in ka:    # i fastest
   ```

   That is libcint's column-major `(i,j,k,l)`. For gradients the three derivative
   directions are slowest, outside all four. Change here if GAMESS's *slot* order
   differs, not just its component order within a slot.

3. The driver's scatter (`src/cint_rotaxis_2e.f90:one_block`) — where `res` is
   written into the caller's buffer. Change here only if terco's output buffer
   layout differs from libcint's; do not use it to fix component order.

### 4.2b Order is not the only difference: GAMESS also rescales

**`--order gamess` gets you the permutation and nothing else.** GAMESS
additionally scales components individually *within* a shell, where libcint —
and therefore libfint and terco — normalise per l only.

`int2a.src`'s `GENRAL` carries one running factor `DUM1` across the component
loop. Its computed `GO TO` (near line 1251, indexed by the flat component
number 1..84) sends the first component of each l to a label that *sets*
`DUM1` from the contraction coefficient, some later components to labels that
*multiply* it, and everything else to label 220, which leaves it alone. The
factor is therefore cumulative and applies to a run of components:

| l | factor 1 | then | then | then |
|---|---|---|---|---|
| d | XX YY ZZ | √3 — XY XZ YZ | | |
| f | XXX YYY ZZZ | √5 — the six XXY-type | √15 — XYZ | |
| g | XXXX YYYY ZZZZ | √7 — the six XXXY-type | √(35/3) — XXYY XXZZ YYZZ | √35 — XXYZ XYYZ XYZZ |

Those factors have a closed form — the standard Cartesian normalisation ratio,

    factor(x^a y^b z^c)^2  =  (2l-1)!! / ((2a-1)!! (2b-1)!! (2c-1)!!)

which reproduces every entry decoded from the jump table exactly, and is
`gamess_norm_sq()` in `cartorder.py` (squared and returned as a `Fraction`, so
the comparisons stay exact). **The rescale therefore needs no transcription at
any l.**

The ordering does. The closed form says the GAMESS order is a grouping by
normalisation factor in nondecreasing order — which is exactly the block
structure of the table above, and what makes the g order look arbitrary until
you see it — but it does **not** determine the order, because the groups are
not singletons: the largest group at both f and g holds six components, and
the factor says nothing about their sequence within it. So `gamess_norm_sq`
gives a genuine import-time check (each table must be sorted by nondecreasing
factor, which catches a component moved across a factor boundary in a way a
permutation check would not) and it cannot stand in for reading `int2a.src`
when l >= 5 is wanted. Confirmed independently in
gamess-libERI, whose `rhf/eric/int0030_ericgen.F90` `angl` array is
`1, sqrt3, …` and `1, sqrt5, …, sqrt15, …` by component type in that code's own
slot order, with `shell_pair.F90` carrying the d √3.

So anyone diffing element by element against gamess-libERI will see sqrt3 /
sqrt5 / sqrt15 mismatches **on top of** any ordering mismatch, and will waste
a day deciding which one they are looking at. Expect both.

`cartorder.py` deliberately stays out of this. It says which component sits in
which slot and nothing more, because conflating a permutation with a diagonal
rescale produces a bug that presents as a permutation bug. The rescale belongs
at the GAMESS boundary as a per-AO diagonal factor — where terco intends to
put it — and an SCF energy is invariant to it anyway, so it only matters when
comparing integrals directly.

### 4.3 What not to do

Do **not** permute after the fact in the driver, and do not hand-edit generated
files. Both work once and then rot: the next regeneration silently reverts them,
and a permutation applied at the wrong layer interacts badly with the L-shell
coefficient-type bookkeeping (§5) and with the canonicalising quartet permutation
in the energy driver.

Change `cart_components`, regenerate, and let the validation suite tell you.

---

## 5. Landmines

These were all found the hard way in libfint. Every one costs a day if
rediscovered.

**L shells are one slot, not two.** Four components sharing exponents, with a
coefficient block **2 × nctr** wide and `NCTR_OF = 1`. Sizing a coefficient block
as `nprim × NCTR_OF` overruns `env` for every L shell. libfint has a shared
helper, `coeff_columns(bas, nbas, sh)`, precisely so this is computed in one
place; the bug was found in `env_extent`, `env_extent_all` *and*
`src/cint_c_abi.f90` independently.

**Gradients cannot permute the quartet.** The energy driver reorders into a
canonical class; `d/dA` names a centre, so it must use the identity permutation
and needs every *ordered* kind combination — 81, not 21.

**`supported()` must ask the dispatcher.** An early version tested an angular-
momentum bound, so a d quartet passed the predicate and then hit an `error stop`
inside. It must call `has_class(code)`.

**The Boys function cancels catastrophically.** `gamma_inc_like`'s
`turnover_point` is 0 for m ≤ 1, so the closed form loses all precision for F₁ at
T ~ 1e-33. libfint has a series branch below T = 1 and a far-field branch above
t = 50 that drops both `erf` and `exp`. Check this early — `gamess-libERI/boys.F90`
is a separate implementation and may or may not share the problem.

**Contract in stages, not flat.** Accumulating over every (bra column, ket column)
pair per primitive quartet costs the *product* of contraction counts. On cc-pVDZ
carbon that is 81 accumulations × 6,561 quartets. Stage it: fill a bra buffer,
then fold the ket. This was found by profiling in mqc, not by reading.

**Watch integer division in the scatter.** Eight divisions per output element made
`one_block` 29% of CPU. Walk component-then-contraction with carried strides.

**Chunk long procedures.** ≤250 statements per contained procedure. See §2.5.

**One module name per file.** A copied benchmark kept its source's module name and
ninja failed with "multiple rules generate test/fb_ws.mod" — an obscure error for
a trivial cause.

**`emit_kernel` had a latent `bp` stride bug, fixed 2026-09-27.** The gradient
work added the bra exponent as `bp` row 6 and changed the declaration in
`emit_kernel` *and* `emit_kernel_tabled` unconditionally. The committed energy
kernels declared `bp(5, nbra)`; regenerating produced `bp(6, nbra)` against a
driver that fills five rows — an explicit-shape stride mismatch that silently
misreads every bra primitive past the first. Fixed to
`bp({6 if grad else 5}, nbra)`. All four generated sets now reproduce
byte-identically, which is the invariant to hold: **after any generator change,
regenerate and diff against the committed tree.** If the diff is not empty and
you did not intend it, stop.

---

## 6. Validation

libfint's chain, which you should mirror:

```
libcint (C)  ==bit-identical==>  libfint Rys path  ==scaled-tolerance==>  rotaxis / HGP
```

- `test/rys_check.f90`, `rys_core_check.f90` + `rys_reference.c` hold the Rys path
  bit-identical to libcint. `scripts/fetch_libcint.sh` supplies the C; CI fetches
  rather than vendors.
- `test/rotaxis_check.F90` and `hgp_check.F90` compare the algebraic paths against
  libfint's own Rys path, so they need no C to run.
- `test/rotaxis_grad_molecule.F90` validates gradients by finite differences.
- `.github/workflows/generators-ci.yml` runs the Python layer: the recurrence
  check above, plus a regeneration of all four kernel sets that must leave the
  tree clean. That second one is the invariant this document keeps pressing —
  it is what catches a generator that has quietly stopped reproducing what it
  shipped, which is how the `bp` stride bug survived. Verified to fire: with
  the bug reintroduced it reports 21 changed files.

**The tolerance argument matters and is easy to get wrong.** rotaxis vs Rys
*cannot* be bit-identical — not because either is approximate, but because two
exact algorithms round differently. The bar is agreement scaled to **the largest
integral in the quartet**, not to the size of each answer, because both
algorithms' rounding is set by the largest term entering a sum. Comparing
element-wise relative error will produce false failures on near-cancelling
components and teach you to distrust a correct kernel.

Coverage in `rotaxis_check.F90`: every quartet of the reference bases up to the
generated l, Cartesian and spherical; then adversarial systems that push the Boys
argument to both ends; then a packed L-shell basis, because L shells are where
the column and component bookkeeping is most easily wrong.

**Suggested order for terco.** Generate with libcint ordering first and validate
against libfint. Only then flip `cart_components` to GAMESS order and validate the
permutation as a separate step against gamess-libERI. Two changes at once and you
will not know which one broke.

Note `rys_core_check` is red on Jorge's machine and green in CI — environmental,
known, do not chase it.

---

## 7. What you will have to write that libfint does not have

- **A GPU emitter.** `emit.py` is the only file that knows Fortran.
  `cartorder.py`, `poly.py`, `derive.py`, `factor.py` are target-agnostic. Write
  `emit_docon.py` beside `emit.py` rather than editing it, so both targets stay
  generatable from one derivation.
- **`do concurrent` loop structure.** libfint's kernels are called per shell
  quartet from a CPU driver. terco's launcher batches; the three loop levels
  (§2.4) map onto that differently, and level 1 (per bra primitive, component-
  independent) is the natural innermost.
- **Rys on the GPU.** See the `c_g2e_unrolled_to_fortran.py` warning in §3.
- **GAMESS's shell-pair data layout.** `gamess-libERI/rot_axis_generator/shell_config.py`
  has it: `PAIR_DATA` per `(la,lb)` with `loc`, `t_expon_ab`, `t_alpha`, `d_coeff`,
  `sq`, `ismlp`, and `XINT_NAMES` for the screening arrays. libfint's `bp`/`kp`
  tables (§2.5) carry the same information in a different shape; the mapping is
  mechanical but it is real work and it belongs in the emitter, not the driver.

## 8. Where things are

```
scripts/cartorder.py                                        the Cartesian order, once
scripts/rotaxis_mmd/{poly,derive,factor,emit,generate}.py   the rotaxis generator
scripts/hgp/{cart,kinds,recur,emit,generate,check_numeric}.py   the HGP generator
scripts/{gen_f90_all,f90-backend}.cl                        the Rys/Lisp backend
scripts/c_*_to_fortran.py                                   C transcribers
src/rotaxis/, src/hgp/, src/rotaxis_grad/, src/hgp_grad/     generated kernels
src/cint_{rotaxis,hgp}[_grad]_kernels.f90                    generated dispatchers
src/cint_{rotaxis,hgp}[_grad]_2e.f90                         hand-written drivers
include/libfint.f90                                          public ABI (libcint_fortran)
doc/{ROT_AXIS_MMD,HGP,PORTING,PERFORMANCE}.md                the prose
```

Dispatch policy deliberately does **not** live in libfint. libfint exposes
algorithms and `*_supported()` queries; the caller chooses, so a user can demand
all-Rys the way GAMESS offers "best vs your_algo". Keep that in terco: expose the
paths, let the engine pick.

---

*Written for the terco session by the libfint session, 2026-09-27. Questions back
to the libfint session — especially before hand-editing anything under `src/rotaxis*`
or `src/hgp*`.*
