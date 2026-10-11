# Rotated-axis MD for terco: the exact recipe

This is the numerically-verified recipe for a Fortran rotaxis emitter that
targets terco's per-class kernel ABI (`src/generated/trc_pc_kLLLL.F90`,
`src/trc_hgp.F90`'s `build_pairs_hgp`). Every claim below was checked by
`scripts/rotaxis/check_numeric.py` against `pyscf`'s `int2e_cart` on all 16
ordered s/p classes, including the degenerate geometries; see the report
handed back with this file for the run output. Nothing here is derived by
hand -- it is read off `scripts/rotaxis/{derive,factor}.py` (vendored
verbatim from libfint) and libfint's `scripts/rotaxis_mmd/emit.py` /
`src/cint_rotaxis_2e.f90`, translated onto terco's own `pp_*` arrays.

## 0. Ordered classes need no permutation

`derive.Derivation(ka, kb, kc, kd)` and `factor.Factorised` place **no**
canonicalisation requirement on their four kind arguments -- canonical
ordering (`KIND_RANK` ascending) is a choice `generate.py`'s CLI makes when
picking which classes to emit for libfint's *canonicalising* driver, not a
restriction inside the derivation or factoriser. Verified directly:

```
Derivation('p','s','s','s')  ->  Factorised: psss  L=1 comps=3 S=3 g=2 r=3 ...
Derivation('s','p','s','s')  ->  Factorised: spss  L=1 comps=3 S=3 g=2 r=3 ...
```

Both construct and evaluate correctly, and `check_numeric.py` runs all 16
ordered classes -- `Derivation(*kinds)` for every `kinds` in
`{'s','p'}^4` -- directly, with **no permutation trick**. This matters for
terco specifically because terco's kernels are one per *ordered* class
(`trc_pc_k1111.F90`'s vbuf digestion is written for shells `(si,sj,sk,sl)`
in the caller's own order, never reordered), so the emitter should derive
each of the 16 classes directly with `Derivation(ka, kb, kc, kd)` in the
literal order the kernel needs, exactly the way `trc_pc_kernels` already
has one module per ordered `(la,lb,lc,ld)` in `{0,1,2}^4`. There is no
"permute then un-permute" step to design or get wrong.

(Gradients would be different -- `README_TERCO.md` section 2.6 notes energy
classes are libfint-canonicalised while gradient classes are not, because
`d/dA` names a centre. For terco's energy kernels the point is moot: nothing
here is ever canonicalised in the first place.)

## 1. The frame, per quartet

Four lab-frame centres `A, B, C, D` (`pp_ra`/`pp_rb` of the bra and ket
pairs -- see section 2). Build:

```
ab = B - A;  rab = |ab|
zh = ab/rab                          if rab > 0
     [0, 0, 1]  (lab z)              if rab == 0   -- A==B, direction undefined

cd = D - C;  rcd = |cd|
w  = cd - (cd . zh) zh                -- component of CD perpendicular to zh
wn = |w|
xh = w/wn                             if wn > 1e-12 * max(rcd, 1.0)
     -- else: CD parallel to AB (or C==D so cd=0): pick the lab axis least
     -- aligned with zh, project out its zh component, normalise.
     kmin = argmin_i |zh_i|            (first index attaining the min)
     w  = e_kmin - (e_kmin . zh) zh
     xh = w / |w|

yh = zh x xh
```

`rot` is the 3x3 matrix whose **columns** are `xh, yh, zh` (local axes
expressed in lab coordinates): `rot(:,1)=xh; rot(:,2)=yh; rot(:,3)=zh`. It is
always a proper rotation, even in every degenerate branch above -- both
fallbacks pick an axis perpendicular to whatever `zh` turned out to be, so
`rot` never needs to be checked for validity afterward.

The threshold `1e-12 * max(rcd, 1.0)` is libfint's
(`cint_rotaxis_2e.f90:local_frame`); do not substitute a bare `1e-12`, the
`max(rcd,1)` scales it to the CD bond length so a long, nearly-parallel CD
doesn't fall through as "parallel" when it isn't (and a very short one
isn't falsely treated as non-degenerate).

`C` and `D` in local coordinates, relative to `A`:

```
cl = (xh.(C-A), yh.(C-A), zh.(C-A))
dl = (xh.(D-A), yh.(D-A), zh.(D-A))
Rs = dl(1) - cl(1)        ! R_CD sin(gamma)
Rc = dl(3) - cl(3)        ! R_CD cos(gamma)
cy = cl(2)                ! the y-offset of C (and of every Q on the ket)
```

`gc = (Rs, Rc, cy)` is the vector `factor.py`'s `cm_list` (`Rs`, `Rc`, `cy`
monomials) is built over, and is exactly libfint's `gc(3)`.

Degenerate cases checked, all pass at the check's tolerance (`1e-10`
scaled, actual worst measured `3.07e-14`, i.e. double-precision noise):

| case | what's undefined | what the frame does |
|---|---|---|
| `A==B` | z axis (`rab=0`) | falls back to lab z |
| `C==D` | x axis (`w=0` since `cd=0`) | falls back to lab-axis projection |
| `AB parallel CD` | y axis (`AB x CD = 0`) | same x fallback (w's non-parallel component vanishes) |
| all four coincide | both of the above at once, **and** `T=0` exactly | both fallbacks fire; `T=0` hits the Boys `t<1e-14` branch, not a division |
| tiny `T` (~1e-18, two nearly-superposed dense pairs) | none (frame is fine); stresses the Boys function only | closed-form Boys per order (section 3), no recursion to cancel |
| huge `T` (~1e7-1e8, distant tight pairs) | none | same closed form; no far-field branch needed because `scipy.special.gammainc`/`gammaln` don't lose precision here (checked) |

## 2. Bra primitive pair, in terms of terco's `pp_*`

terco's bra pair table comes from `build_pairs_hgp` (`src/trc_hgp.F90:83`)
called on the bra shell pair `(si, sj)`. For primitive pair `kp` in
`[pp_off(key)+1, pp_off(key)+pp_n(key)]`, `key = (si-1)*nbas + sj`:

```
p      = pp_p(kp)                                   ! a + b
A      = pp_ra(kp, :)                                ! shell si's centre
B      = pp_rb(kp, :)                                ! shell sj's centre
rab    = |A - B|                                      ! from the frame, section 1
ip     = 1/(2p)
za     = (P - A) . zh, but since P = pp_r(kp,:) lies on the bra-shell
         line by construction, this is simply  b*rab/p, i.e.
         za = (pp_p(kp) - a)*rab/pp_p(kp)     where a is shell si's own
         primitive exponent -- in practice just project:
             za = zh . (pp_r(kp,:) - A)
zb     = za - rab
```

(`za` can be had two ways -- projecting the already-built `pp_r(kp,:)`
onto `zh`, or from the bare exponents if the emitter kept them; both are
the same number. `check_numeric.py` uses the exponent form via
`bra_pairs()`, which is what an emitted kernel will actually have since
`pp_p` gives `p` but not `a`/`b` separately -- so **project `pp_r`,
`pp_ra` onto the frame's `zh`** rather than trying to recover `a`,`b` from
`pp_p` alone.)

The bra weight -- what `factor.py`'s `kab` (§2.5 of README_TERCO.md) wants,
`c_i c_j exp(-mu R_AB^2)/p` -- is **exactly** `pp_c(kp)/pp_p(kp)`:

```
kab(kp) = pp_c(kp) / pp_p(kp)
```

`pp_c(kp)` already **is** `c_i c_j exp(-mu R_AB^2)` (verified from
`src/trc_hgp.F90:179`: `pp_c(k) = sh_c(ki,i)*cfac(i)*sh_c(kj,j)*cfac(j)*
exp(-mu*ab2)`, with `cfac` always 1.0 on terco's per-class kernel path
because `trc_decontract` has already split every shell to one column
before this table is built -- see `src/trc_fock.F90:306-309`, which passes
`one` for `cfac`). `sh_c` is **not** the raw basis-file coefficient: see
section 4 for what it already carries.

## 3. Ket primitive pair

Symmetric, from the ket shell pair `(sk, sl)`, `key = (sk-1)*nbas + sl`:

```
q   = pp_p(kq)
C   = pp_ra(kq, :)
D   = pp_rb(kq, :)
rcd = |C - D|
iq  = 1/(2q)
yc  = c/q,  yd = d/q       ! c, d the two ket exponents; yc+yd=1
ax  = cl(1) + yd*Rs        ! cl, Rs from the frame (section 1)
az  = cl(3) + yd*Rc
```

`yc`, `yd` again come from the exponents; recoverable as
`yd = (Q - C).zh / Rc`-style projections if the emitter only kept `pp_r`,
but simplest is to keep the raw per-primitive exponents `c`, `d` (or
equivalently `yd = (pp_p(kq) - c_exp)/pp_p(kq)`) the same way the bra side
does.

The ket weight is `TWO_PI_2_5 * pp_c(kq) / pp_p(kq)`:

```
kcd(kq) = TWO_PI_2_5 * pp_c(kq) / pp_p(kq)
```

`TWO_PI_2_5` is terco's own named constant, already declared in every
per-class kernel module (`real(dp), parameter :: TWO_PI_2_5 =
34.986836655249725_dp`, `src/generated/trc_pc_k1111.F90:19`) -- reuse it,
don't redeclare.

## 4. Per primitive quartet: `T`, `B_n`, and where the prefactors land

```
zq    = za - az
alpha = p*q/(p+q)
T     = alpha*(ax*ax + cy*cy + zq*zq)
F_n   = boys(n, T),  n = 0..L               ! L = sum of the four kinds' l
B_n   = F_n(T) * (-2*alpha)^n / sqrt(p+q)
```

Putting `kab`, `kcd`, `B_n` together, one primitive quartet's contribution
to the `(ss|ss)` term is

```
kab(kp) * kcd(kq) * B_0
  = [pp_c(kp)/p] * [TWO_PI_2_5*pp_c(kq)/q] * [F_0(T)/sqrt(p+q)]
  = TWO_PI_2_5/(p*q*sqrt(p+q)) * pp_c(kp)*pp_c(kq) * F_0(T)
```

which is **bit-for-bit the same grouping** terco's own HGP kernel already
computes (`src/trc_hgp.F90:306`, `src/trc_binkernel.F90:900`,
`src/trc_batch.F90:386`):

```fortran
pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)
```

with `zeta=p`, `eta=q`, `zpe=p+q` -- so a rotaxis kernel built from
`kab`/`kcd`/`B_n` this way lands on **exactly** terco's existing
normalisation and prefactor convention with no extra factor to apply
anywhere else. This was verified, not assumed: `check_numeric.py` builds
`kab`/`kcd` this way and its output matches `pyscf`'s `int2e_cart` (which
is the same convention terco's shipped HGP kernel already targets) to
`3e-14` scaled, worst case, over all 16 classes and all 8 geometries per
class.

There is **no separate per-l normalisation step** to apply afterward the
way libfint's driver applies `cint_common_fac_sp` post-kernel
(`cint_rotaxis_2e.f90:one_block`, the `fc(...)` array) -- terco folds
`common_fac_sp` into `sh_c` **once, at basis-build time**
(`src/trc_api.F90:284`, `docs/GENERATORS.md` "Shared conventions"), so it
is already inside `pp_c` and therefore already inside `kab`/`kcd` above.
**Do not apply `common_fac_sp` a second time in the kernel or the driver.**
This is the one place terco's convention differs from libfint's own
rotaxis driver, and it is exactly the kind of thing `docs/GENERATORS.md`
warns has bitten this repository three times already (folding a smooth
per-shell factor in twice).

## 5. Rotating the local-frame result back to the lab frame

Only `p` (l=1) indices need rotating; an `s` index is 1-dimensional and the
identity. For a `p` index, the local-frame components are in
`cart_components(1) = [(1,0,0),(0,1,0),(0,0,1)]` order, i.e. **local x, y,
z**, and

```
lab[i] = sum_j rot(i,j) * local[j]           i,j in {x,y,z} = {1,2,3}
```

i.e. `lab = rot @ local`, using the **same** `rot` built in section 1 for
the whole quartet (one shared frame, not a per-slot rotation). This is
`cint_rotaxis_2e.f90:index_matrix`'s `case(1); tmat = rot` and
`transform_index`'s `out(io) = sum_ii tmat(io,ii)*loc(ii)`, confirmed
numerically here (every class with a `p` slot passes the pyscf check,
including `pppp` where all four indices rotate).

Applied per-index, independently, to whichever of the four slots is a `p`
kind -- the transform is separable across indices exactly because `rot` is
the same rigid rotation applied to each Cartesian index of the tensor, one
axis at a time (`check_numeric.py`'s `rotate_axis`/`evaluate_quartet`:
reshape the flat local-frame result to `(nd,nc,nb,na)` axes and apply
`np.tensordot`-style contraction with `rot` on each axis whose kind is
`'p'`, `s` axes untouched).

**Generalising to `d` (not implemented here, terco's rotaxis path is
`--lmax 1` only per `README_TERCO.md` section 2.6, and this task is s/p
only):** a local lab Cartesian `x^i y^j z^k` (l=2) is the product of two
factors each drawn from `rot`'s three rows dotted into the local axes, i.e.
the 6x6 matrix is built entry-by-entry as sums of products
`rot(dir1,a)*rot(dir2,b)` over which pair of local monomials `(a,b)`
combine to the lab monomial in question -- see
`cint_rotaxis_2e.f90:index_matrix`'s `case(2)` block for the exact
construction (it enumerates the two Cartesian "factors" of each lab
monomial and sums `rot(fa,a)*rot(fb,b)` over the local index pairs `(a,b)`
that produce it). Do not attempt this until the s/p emitter is done and
its own check is green; `doc/ROT_AXIS_MMD.md` section 9 records that d was
generated but not committed to libfint's rotaxis path because it is 2.5x
*slower* than Rys there, so it may not be worth porting to terco either
without first measuring.

## 6. Component order: derive's index vs terco's `vbuf` index

**Confirmed identical**, no permutation needed anywhere in the pipeline.

- `derive.Derivation.all_components()` iterates
  `for cd: for cc: for cb: for ca:` -- **`ca` (slot A) fastest**, then B,
  then C, then D slowest. This is libcint's column-major `(i,j,k,l)`.
- `src/generated/trc_pc_k1111.F90:970-976` (the digestion loop) reads:
  ```fortran
  idx = 0
  do id = 0, 2
  do ic = 0, 2
  do ib = 0, 2
  do ia = 0, 2
     idx = idx + 1
     sc = wq*vbuf(idx)
     ...
  ```
  i.e. `idx` increments with **`ia` fastest**, then `ib`, `ic`, `id`
  slowest -- the identical nesting order.

So a rotaxis kernel's local-frame `Factorised.components` list (or its
lab-frame-rotated counterpart from section 5) can be written straight into
`vbuf(1:ncomp)` in derivation order with **no index remapping**. This was
also confirmed empirically: `check_numeric.py`'s `evaluate_quartet` builds
its result in exactly this order and compares element-for-element (no
permutation applied) against `pyscf`'s Cartesian block, flattened with the
identical `for id: for ic: for ib: for ia` nesting
(`flatten_pyscf_order`).

## 7. Level sizes for all 16 ordered classes

From `Factorised.summary()`, `ka kb kc kd` in `{s,p}^4`, direct
(unpermuted) `Derivation`:

```
ssss: L=0 comps=1  S=1  g=1 r=1   (terms 1)   km=1  cm=1  assembly=1
sssp: L=1 comps=3  S=3  g=1 r=4   (terms 4)   km=3  cm=4  assembly=5
ssps: L=1 comps=3  S=3  g=1 r=4   (terms 4)   km=3  cm=4  assembly=5
sspp: L=2 comps=9  S=6  g=1 r=16  (terms 23)  km=9  cm=10 assembly=24
spss: L=1 comps=3  S=3  g=2 r=3   (terms 4)   km=2  cm=2  assembly=3
spsp: L=2 comps=9  S=8  g=2 r=12  (terms 18)  km=5  cm=7  assembly=16
spps: L=2 comps=9  S=8  g=2 r=12  (terms 18)  km=5  cm=7  assembly=16
sppp: L=3 comps=27 S=15 g=2 r=48  (terms 104) km=14 cm=16 assembly=80
psss: L=1 comps=3  S=3  g=2 r=3   (terms 4)   km=2  cm=2  assembly=3
pssp: L=2 comps=9  S=8  g=2 r=12  (terms 18)  km=5  cm=7  assembly=16
psps: L=2 comps=9  S=8  g=2 r=12  (terms 18)  km=5  cm=7  assembly=16
pspp: L=3 comps=27 S=15 g=2 r=48  (terms 104) km=14 cm=16 assembly=80
ppss: L=2 comps=9  S=10 g=5 r=9   (terms 21)  km=3  cm=3  assembly=10
ppsp: L=3 comps=27 S=24 g=5 r=36  (terms 95)  km=7  cm=10 assembly=54
ppps: L=3 comps=27 S=24 g=5 r=36  (terms 95)  km=7  cm=10 assembly=54
pppp: L=4 comps=81 S=43 g=5 r=143 (terms 549) km=19 cm=22 assembly=271
```

(`S`, `r`, `km`, `cm` are the accumulator/temporary counts of section 2.4
of README_TERCO.md; `assembly` is the total number of `(scale, cm, r)`
terms summed over all components at level 3.) The `pppp` row matches
`doc/ROT_AXIS_MMD.md`'s worked numbers exactly (43 accumulators over 5
distinct `g`, 143 ket-level accumulators), which cross-validates this
reading of `factor.py` against libfint's own prose. Note the asymmetry
between e.g. `sspp` (ket is `pp`, S=6) and `ppss` (bra is `pp`, S=10) --
the level-1 (bra) accumulator count depends on which pair carries the
angular momentum, because the bra loop is the one that varies fastest and
the ket polynomial substitution happens at level 2 regardless of which
side textually "has" the `p`.

## 8. What was NOT true and had to be checked rather than assumed

- Terco's `pp_c` is **not** "just" `exp(-mu R_AB^2)`; it already carries
  both contraction coefficients **and** `common_fac_sp(l)`. Confirmed
  empirically against `pyscf`'s raw (un-normalised-by-`common_fac_sp`)
  `_env` coefficients: multiplying by `common_fac_sp(l)` once per shell and
  nothing else reproduces `pyscf.intor('int2e_cart')` exactly for `ssss`.
- `common_fac_sp` is folded in **at basis-build time** in terco, not
  **after the kernel** the way libfint's own rotaxis driver applies it.
  A straight port of `cint_rotaxis_2e.f90`'s `fc(...)` scatter step would
  double-apply it for every quartet with an `s` or `p` index -- i.e. every
  quartet, since terco doesn't generate anything above `p` here. This is
  the single most important convention difference from libfint's own
  driver and is called out explicitly in section 4.
- `Derivation`/`Factorised` do **not** require canonical order -- section 0.
  This was checked, not assumed, because the alternative (permute-then-
  un-permute) would have been meaningfully more code and more places to
  get an index wrong.

## 9. Nothing here disagreed with pyscf

Every one of the 16 ordered classes, at 8 geometries each (2 generic random
plus the 6 adversarial cases of section 1), passed at a worst scaled error
of `3.07e-14` (`pppp`, a generic geometry) against a tolerance of `1e-10`.
That worst case is ordinary double-precision rounding, consistent with
libfint's own measured `5e-14` agreement between its rotaxis and Rys paths
(`doc/ROT_AXIS_MMD.md` section 7). There was no disagreement to paper over.

## 10. The role-swapped mapping (bra/ket-inner tuning, milestone-1 follow-up)

HGP's own kernels put the ket primitive loop innermost and the terco decode
that dispatches quartets across threads (`trc_pc_kernels`' segment decode)
makes the BRA pair index vary fastest across a warp, so 32 threads share one
ket pair and the innermost (ket) loop's loads broadcast. Sections 1-4 above
put rotaxis's bra loop innermost instead (level 1's accumulators are per
bra primitive), which is right for the factorisation but wrong for that
decode: the innermost loads now diverge across the warp instead of
broadcasting. Milestone-1's own review measured this and tried changing the
*decode* to match (making the ket index warp-uniform instead); it worked
for some classes and lost for others (mixed, moderately below the ~10-15%
noise floor for some, clearly outside it for the losses) and was not
shipped -- see the commit message and `emit_class_module`'s comment in
`emit_docon.py`.

This section is the other way to fix the same mismatch: change the
**algebra**, not the decode. Derive `Derivation(kc, kd, ka, kb)` instead of
`Derivation(ka, kb, kc, kd)`. Nothing in `derive.py`/`factor.py` cares which
physical pair is called "bra" and which "ket" -- section 0 already
established that no canonicalisation is required -- so this is legal for
the same reason ordered classes are: `(ab|cd) = (cd|ab)` exactly, index for
index, because it is the same four-centre Coulomb integral viewed from the
other pair's side, not merely the same integral after summing.

**The frame.** Build it from `(C, D, A, B)` instead of `(A, B, C, D)` --
i.e. call section 1's construction with its four inputs in that order. Every
formula is unchanged; only which physical centres feed "A", "B", "C", "D" in
it changes:

```
z along C -> D                    (was A -> B)
x spans AB perpendicular to CD     (was: spans CD perpendicular to AB)
origin, and cl/dl, at C            (was: at A)
Rs, Rc = the (A, B) pair's local x, z offset from C  (was: (C, D)'s from A)
```

**Level 1 (inner loop) is now over terco's KET primitives.** The
derivation's "bra" pair is physically (C, D): its `p`, `ip`, `za`, `zb` come
from `pp_p(kq)`, projected onto the NEW `zh` (along C->D), for `kq` in the
usual `[offcd+1, offcd+ncd]` range -- this is section 2's recipe, verbatim,
with the roles of "A","B" read as "C","D".

**Level 2 (outer loop) is now over terco's BRA primitives.** The
derivation's "ket" pair is physically (A, B): its `q`, `iq`, `yc`, `yd`,
`ax`, `az` come from `pp_p(kp)` for `kp` in `[offab+1, offab+nab]`, using
the NEW frame's `cl`/`Rs`/`Rc` -- section 3's recipe, with "C","D" read as
"A","B".

**Component order.** `Derivation(kc,kd,ka,kb).all_components()` iterates
its own four slots fastest-to-slowest as `(kc, kd, ka, kb)` (section 6's
rule applied to the swapped kinds tuple), i.e. **terco's `c` index is
fastest, not `a`**. A component at swapped-order position `(ic, id, ia,
ib)` (0-based within each KIND's own component list) belongs at terco's
`(ia, ib, ic, id)` -- this is a fixed relabelling, known completely at
generation time from the class alone, not a runtime permutation and not a
data-dependent computation. `scripts/rotaxis/check_numeric.py`'s
`swapped_component_map` builds this map directly (by construction, not by
reshaping or transposing a flat array) precisely so an axis-order mistake
cannot hide inside a `.reshape()`/`.transpose()` call.

**Rotation.** Unchanged in kind: rotate whichever of terco's four slots is
`p`, with the frame's `rot` -- the SAME shared frame built in this section,
used for every slot regardless of which pair (bra or ket, in either sense)
it belongs to. Only the frame itself differs from section 5's (built from
`(C,D,A,B)` here rather than `(A,B,C,D)`).

**Verified against pyscf**, not assumed: `scripts/rotaxis/check_numeric.py`
runs BOTH evaluations -- the original (`run()`) and this one
(`run_swapped()`) -- on the identical 16 classes x 8 geometries (including
every degenerate case of section 1; "C == D" is this evaluator's undefined-z
case, the way "A == B" is the original's, since `local_frame` does not know
or care which physical pair its first two arguments came from). Worst
scaled error `3.01e-14`, the same order as section 9's, no disagreement.
