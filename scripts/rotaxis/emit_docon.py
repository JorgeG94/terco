"""Emit terco's rotated-axis MD kernels as `do concurrent` Fortran.

This is the Fortran back end for the recipe in FRAME.md, one module per
ordered s/p class (`{s,p}^4`, 16 classes, see FRAME.md section 0 -- no
canonicalisation). It reuses the scaffolding gen_perclass.py already built
for the HGP per-class kernels (segment locate/decode, the exact Schwarz and
density-weighted screens, same_ab/same_cd/same_pair, the scalar-path AO
decode, the block-accumulated density-batched digestion and its
TRC_NO_DIGEST branch, the scalar launch driver, the item routine signature,
and the inline Chebyshev Boys evaluation) rather than re-implementing any of
it -- see the functions imported from `gen_perclass` below, and
docs/GENERATORS.md's "Shared conventions" for why an independent copy would
drift.

WHAT IS NEW HERE, PER FRAME.md
-------------------------------
Only the integral core: build the frame once per quartet (section 1), the
bra pair quantities inside the (now innermost) bra loop and the ket pair
quantities inside the (now outermost) ket loop (sections 2-3), accumulate
the three factorisation levels of `scripts/rotaxis/factor.py` (section 2.4
of libfint's README_TERCO.md: S per bra monomial, r per ket polynomial,
then the quartet assembly), and rotate every `p` slot of the result back to
the lab frame (section 5) before handing `vbuf` to the shared digestion.

The level generation itself -- which monomials become accumulators, how a
term's coefficient decomposes into a bra part, a ket part and a constant
part -- follows libfint's `scripts/rotaxis_mmd/emit.py` (`emit_kernel`,
non-tabled form; every one of the 16 s/p classes is far under its
UNROLL_LIMIT of 3000 ket-level terms, so table mode is never needed here).
That file's algorithm is translated onto terco's naming and single-column
ABI below; nothing is imported from libfint's tree -- terco has no
dependency on it -- and terco's own three-argument (bra pair, ket pair,
quartet) ABI is used throughout rather than libfint's general (ncb, nck,
NTT) one, since every terco per-class kernel already decontracts to one
column per shell (`trc_decontract`; see FRAME.md section 2's remark on
`cfac`).
"""

import os
import re
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.dirname(HERE)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)
if HERE not in sys.path:
    sys.path.insert(0, os.path.dirname(HERE))

import gen_perclass as gp  # noqa: E402  (the shared scaffolding)
from rotaxis import derive, factor  # noqa: E402

MAX_LINE = 100

#: The 16 ordered s/p classes, in the SAME nested order gen_perclass.py
#: walks (la,lb,lc,ld), so the two dispatchers list classes the same way.
ALL_CLASSES = [(a, b, c, d) for a in "sp" for b in "sp" for c in "sp" for d in "sp"]

_KIND_L = {"s": 0, "p": 1}


def tag_of(kinds):
    return "".join(str(_KIND_L[k]) for k in kinds)


def class_key(kinds):
    """The pc_dispatch/ra_dispatch key: kind ranks as digits of CLASS_RADIX.
    's' and 'p' rank the same way gen_perclass's angular momentum does (0
    and 1), since a class here IS just (la,lb,lc,ld) in {0,1}."""
    la, lb, lc, ld = (_KIND_L[k] for k in kinds)
    return ((la * gp.CLASS_RADIX + lb) * gp.CLASS_RADIX + lc) * gp.CLASS_RADIX + ld


# ---------------------------------------------------------------------
# Small Fortran-expression helpers, translated from libfint's
# scripts/rotaxis_mmd/emit.py (fnum/mono_expr/coef_term/wrap) onto terco's
# naming.  Pure functions, no state; see that file's docstring for the
# algebra these implement.
# ---------------------------------------------------------------------

def fnum(c):
    c = Fraction(c)
    if c.denominator == 1:
        v = c.numerator
        return f"{v}.0_dp" if v >= 0 else f"({v}.0_dp)"
    return f"({c.numerator}.0_dp/{c.denominator}.0_dp)"


def mono_expr(names, exps):
    """A monomial as a Fortran product: each nonzero exponent e becomes the
    name repeated e times (no `**`, matching gen_perclass's own style of
    keeping small integer powers as multiplies)."""
    parts = []
    for name, e in zip(names, exps):
        if e == 0:
            continue
        parts.extend([name] * e)
    return "*".join(parts) if parts else None


def coef_term(c, expr):
    c = Fraction(c)
    if expr is None:
        return fnum(c)
    if c == 1:
        return expr
    if c == -1:
        return f"-{expr}"
    return f"{fnum(c)}*{expr}"


def wrap(lhs, terms, indent="            "):
    """`lhs = t1 + t2 + ...`, continued onto further lines before MAX_LINE."""
    if not terms:
        terms = ["0.0_dp"]
    expr = terms[0]
    for t in terms[1:]:
        expr += (" - " + t[1:]) if t.startswith("-") else (" + " + t)
    out = []
    cur = f"{indent}{lhs} = "
    for tok in expr.split(" "):
        if len(cur) + len(tok) + 1 > MAX_LINE and not cur.endswith("= "):
            out.append(cur.rstrip() + " &")
            cur = indent + "   " + tok + " "
        else:
            cur += tok + " "
    out.append(cur.rstrip())
    return "\n".join(out)


# ---------------------------------------------------------------------
# The frame (FRAME.md section 1).  Fixed text: no class dependence, so it
# is built once and reused for every one of the 16 modules, the same way
# gen_perclass's own scaffolding functions are.
# ---------------------------------------------------------------------

def emit_frame():
    return """         ! === the frame (FRAME.md section 1) ===
         ! ab = B - A, the bra pair's own direction; z falls back to the lab
         ! axis when A == B (rab = 0), which happens whenever two shells on
         ! the same atom are paired -- routine, not rare.
         abx = pp_rb(offab + 1, 1) - pp_ra(offab + 1, 1)
         aby = pp_rb(offab + 1, 2) - pp_ra(offab + 1, 2)
         abz = pp_rb(offab + 1, 3) - pp_ra(offab + 1, 3)
         rab = sqrt(abx*abx + aby*aby + abz*abz)
         if (rab > 0.0_dp) then
            zhx = abx/rab; zhy = aby/rab; zhz = abz/rab
         else
            zhx = 0.0_dp; zhy = 0.0_dp; zhz = 1.0_dp
         end if
         cdx = pp_rb(offcd + 1, 1) - pp_ra(offcd + 1, 1)
         cdy = pp_rb(offcd + 1, 2) - pp_ra(offcd + 1, 2)
         cdz = pp_rb(offcd + 1, 3) - pp_ra(offcd + 1, 3)
         rcd = sqrt(cdx*cdx + cdy*cdy + cdz*cdz)
         ! x spans what of CD is perpendicular to z; falls back to the lab
         ! axis least aligned with z when CD is parallel to AB (or C == D,
         ! another same-atom case rather than a rare one) -- FRAME.md's
         ! 1e-12*max(rcd,1) threshold, taken from libfint verbatim.
         dotcz = cdx*zhx + cdy*zhy + cdz*zhz
         wx = cdx - dotcz*zhx; wy = cdy - dotcz*zhy; wz = cdz - dotcz*zhz
         wn = sqrt(wx*wx + wy*wy + wz*wz)
         if (wn > 1.0e-12_dp*max(rcd, 1.0_dp)) then
            xhx = wx/wn; xhy = wy/wn; xhz = wz/wn
         else
            ! z is a unit vector, so at least one lab axis has |component|
            ! <= 1/sqrt(3); projecting z out of THAT axis is never zero, so
            ! this branch's own division needs no further check (FRAME.md).
            if (abs(zhx) <= abs(zhy) .and. abs(zhx) <= abs(zhz)) then
               wx = 1.0_dp - zhx*zhx; wy = -zhx*zhy; wz = -zhx*zhz
            else if (abs(zhy) <= abs(zhz)) then
               wx = -zhy*zhx; wy = 1.0_dp - zhy*zhy; wz = -zhy*zhz
            else
               wx = -zhz*zhx; wy = -zhz*zhy; wz = 1.0_dp - zhz*zhz
            end if
            wn = sqrt(wx*wx + wy*wy + wz*wz)
            xhx = wx/wn; xhy = wy/wn; xhz = wz/wn
         end if
         yhx = zhy*xhz - zhz*xhy
         yhy = zhz*xhx - zhx*xhz
         yhz = zhx*xhy - zhy*xhx
         ! rot's COLUMNS are the local axes in lab coordinates (FRAME.md
         ! section 1); section 5 rotates each `p` slot back with it.
         rot(1,1) = xhx; rot(2,1) = xhy; rot(3,1) = xhz
         rot(1,2) = yhx; rot(2,2) = yhy; rot(3,2) = yhz
         rot(1,3) = zhx; rot(2,3) = zhy; rot(3,3) = zhz
         cl1 = xhx*(pp_ra(offcd + 1, 1) - pp_ra(offab + 1, 1)) &
             + xhy*(pp_ra(offcd + 1, 2) - pp_ra(offab + 1, 2)) &
             + xhz*(pp_ra(offcd + 1, 3) - pp_ra(offab + 1, 3))
         cl2 = yhx*(pp_ra(offcd + 1, 1) - pp_ra(offab + 1, 1)) &
             + yhy*(pp_ra(offcd + 1, 2) - pp_ra(offab + 1, 2)) &
             + yhz*(pp_ra(offcd + 1, 3) - pp_ra(offab + 1, 3))
         cl3 = zhx*(pp_ra(offcd + 1, 1) - pp_ra(offab + 1, 1)) &
             + zhy*(pp_ra(offcd + 1, 2) - pp_ra(offab + 1, 2)) &
             + zhz*(pp_ra(offcd + 1, 3) - pp_ra(offab + 1, 3))
         dl1 = xhx*(pp_rb(offcd + 1, 1) - pp_ra(offab + 1, 1)) &
             + xhy*(pp_rb(offcd + 1, 2) - pp_ra(offab + 1, 2)) &
             + xhz*(pp_rb(offcd + 1, 3) - pp_ra(offab + 1, 3))
         dl3 = zhx*(pp_rb(offcd + 1, 1) - pp_ra(offab + 1, 1)) &
             + zhy*(pp_rb(offcd + 1, 2) - pp_ra(offab + 1, 2)) &
             + zhz*(pp_rb(offcd + 1, 3) - pp_ra(offab + 1, 3))
         Rs = dl1 - cl1
         Rc = dl3 - cl3
         cy = cl2
"""


# ---------------------------------------------------------------------
# Rotating the p slots back to the lab frame (FRAME.md section 5).
#
# `vbuf` is laid out ia-fastest (FRAME.md section 6, matching the shared
# digestion), which is exactly Fortran's own column-major order for a
# (na,nb,nc,nd) array -- so rotating slot X in place is: for every fixed
# combination of the OTHER three slots' indices, replace the three (or one,
# for an s slot, a no-op) consecutive-with-stride entries by rot times
# themselves.  Unrolled fully at generation time (every size here is 1 or
# 3), so nothing runtime-loops over it.
# ---------------------------------------------------------------------

def emit_rotation(kinds, na, nb, nc, nd):
    sizes = [na, nb, nc, nd]
    strides = [1, na, na * nb, na * nb * nc]
    out = []
    for slot in range(4):
        if kinds[slot] != "p":
            continue
        n = sizes[slot]
        assert n == 3
        stride = strides[slot]
        others = [s for s in range(4) if s != slot]
        ranges = [range(sizes[s]) for s in others]
        seen = set()
        for i0 in ranges[0]:
            for i1 in ranges[1]:
                for i2 in ranges[2]:
                    base = 0
                    for s, i in zip(others, (i0, i1, i2)):
                        base += i * strides[s]
                    if base in seen:
                        continue
                    seen.add(base)
                    p0, p1, p2 = base + 1, base + 1 + stride, base + 1 + 2 * stride
                    out.append(f"         t1 = vbuf({p0}); t2 = vbuf({p1}); t3 = vbuf({p2})")
                    out.append(f"         vbuf({p0}) = rot(1,1)*t1 + rot(1,2)*t2 + rot(1,3)*t3")
                    out.append(f"         vbuf({p1}) = rot(2,1)*t1 + rot(2,2)*t2 + rot(2,3)*t3")
                    out.append(f"         vbuf({p2}) = rot(3,1)*t1 + rot(3,2)*t2 + rot(3,3)*t3")
    return "\n".join(out)


# ---------------------------------------------------------------------
# The integral core: frame, bra pair (inner loop), ket pair (outer loop),
# the three factorisation levels, the quartet assembly straight into
# `vbuf`, then the rotation above.  Translated from libfint's
# scripts/rotaxis_mmd/emit.py `emit_kernel` -- see that file for the
# reference algorithm this follows term for term.
# ---------------------------------------------------------------------

def V(base, i, scalarize):
    """A reference to element `i` (1-based) of one of the small per-quartet
    temporaries (sacc, racc, kmv, cmv, gb, bn, bz). Array syntax by default;
    `scalarize` switches every one of these to an individually NAMED scalar
    (`sacc_3` rather than `sacc(3)`), which is the only difference between
    the two kernel forms `--scalarize` compares (see gen_rotaxis.py). Every
    index used here is a Python-side literal, so nvfortran sees a compile-
    time constant subscript either way; the point of the comparison is
    whether ptxas still needs to be TOLD that to keep the array out of local
    memory, or already works it out on its own. `vbuf`, `f` (the Boys
    ladder) and the digestion's blocks are never touched by this -- they are
    either shared text with gen_perclass.py or, for `f`, the argument
    `emit_boys` itself expects -- so scalarizing is confined to the six
    temporaries the fold in emit_core actually shrank."""
    return f"{base}_{i}" if scalarize else f"{base}({i})"


def needs_frame(f):
    """False only for ssss among the 16 s/p classes (coordinator review item
    A): true unless EVERY g monomial is trivial (no za/zb dependence), every
    km monomial is trivial (no yc/yd/ax dependence), every cm monomial is
    trivial (no Rs/Rc/cy dependence) and no s_list entry carries a zq power.
    When all four hold, the only quantity the primitive loop ever needs is
    B_0 = F_0(T)/sqrt(zpe), and T = alpha*(ax^2 + cy^2 + zq^2) is exactly
    alpha*|P-Q|^2 -- a rotation is an isometry, so that sum of squares does
    not care what frame ax, cy, zq were measured in. Building the frame
    (two sqrt, a cross product, four dot products, ~15 live doubles) to
    recover a quantity that is already frame-independent is pure overhead
    for exactly this one class."""
    g_trivial = len(f.g_list) == 1 and not any(f.g_list[0])
    km_trivial = len(f.km_list) == 1 and not any(f.km_list[0])
    cm_trivial = len(f.cm_list) == 1 and not any(f.cm_list[0])
    return not (g_trivial and km_trivial and cm_trivial and f.max_zq_power == 0)


def emit_core_noframe(f, scalarize=False):
    """ssss's fast path: T computed directly from the lab-frame separation
    P - Q, no rotated frame at all. See `needs_frame` for why this is exact,
    not an approximation -- every quantity the frame would have supplied
    (za, zb, yc, yd, ax, cy, zq individually; Rs, Rc; the rotation itself)
    is provably unused by this class's algebra."""
    def sacc(i): return V("sacc", i, scalarize)
    def racc(i): return V("racc", i, scalarize)
    def kmv(i): return V("kmv", i, scalarize)

    assert f.d.L == 0 and len(f.s_list) == 1 and len(f.r_list) == 1

    o = [f"""
         ! === per-quartet integral core, NO-FRAME FAST PATH ===
         !
         ! ssss's only accumulator is B_0 = F_0(T)/sqrt(zpe); T needs only
         ! |P - Q|^2, so it is formed directly from pp_r in the lab frame
         ! (needs_frame's docstring proves this is exact for this class, not
         ! an approximation). No rot, no xh/yh/zh, no cl/dl/Rs/Rc/cy, no
         ! per-primitive za/zb/yc/yd/ax/az -- none of them would change the
         ! answer, and every one of the frame's sqrt/cross/dot calls would
         ! have been pure overhead for a class whose entire job is one Boys
         ! evaluation per primitive quartet.
         !
         {racc(1)} = 0.0_dp
         do kq = offcd + 1, offcd + ncd
            eta = pp_p(kq)
            bnd = TWO_PI_2_5*abs(pp_cs(kq))/(eta*sqrt(eta))
            if (bnd*abs(pp_cs(offab + 1))/pp_p(offab + 1) <= pcut) cycle
            kcdv = TWO_PI_2_5*pp_cs(kq)/eta
            {sacc(1)} = 0.0_dp
            do kp = offab + 1, offab + nab
               zeta = pp_p(kp)
               if (bnd*abs(pp_cs(kp))/zeta <= pcut) exit
               zpe = zeta + eta
               ! The exact (ss|ss)-type bound: bit-for-bit the HGP kernel's
               ! own `pref` (FRAME.md section 4).
               pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_cs(kp)*pp_cs(kq)
               if (abs(pref) <= pcut) cycle
               pqx = pp_r(kp, 1) - pp_r(kq, 1)
               pqy = pp_r(kp, 2) - pp_r(kq, 2)
               pqz = pp_r(kp, 3) - pp_r(kq, 3)
               alpha = zeta*eta/zpe
               tval = alpha*(pqx*pqx + pqy*pqy + pqz*pqz)

   {gp.emit_boys(0)}

               kabv = pp_cs(kp)/zeta
               w = kabv/sqrt(zpe)
               {sacc(1)} = {sacc(1)} + w*f(0)
            end do   ! kp (bra primitive)
            {kmv(1)} = kcdv
            {racc(1)} = {racc(1)} + {kmv(1)}*{sacc(1)}
         end do   ! kq (ket primitive)

         ! === quartet assembly, straight into vbuf (level 3) ===
         vbuf(1) = {racc(1)}

         ! === rotate every p slot back to the lab frame (FRAME.md section 5) ===
         ! every slot is s: the identity, nothing to rotate.
"""]
    return "\n".join(o), dict(NS=1, NR=1, NG=0, NKM=1, NCM=0, NCOMP=1, L=-1, KMAX=0,
                              NBZ=0, na=1, nb=1, nc=1, nd=1, noframe=True)


def emit_core(kinds, scalarize=False):
    d = derive.Derivation(*kinds)
    f = factor.Factorised(d)
    L = d.L
    NS, NR, NG = len(f.s_list), len(f.r_list), len(f.g_list)
    NKM, NCM = len(f.km_list), len(f.cm_list)
    NCOMP = len(f.components)
    KMAX = f.max_zq_power
    na, nb, nc, nd = (len(derive.KINDS[k]) for k in kinds)
    assert NCOMP == na * nb * nc * nd
    for g in f.g_list:
        assert g[3] == 0, "energy classes never touch the gradient's `ea` factor"
    if not needs_frame(f):
        return emit_core_noframe(f, scalarize=scalarize)

    def sacc(i): return V("sacc", i, scalarize)
    def racc(i): return V("racc", i, scalarize)
    def kmv(i): return V("kmv", i, scalarize)
    def cmv(i): return V("cmv", i, scalarize)
    def gb(i): return V("gb", i, scalarize)
    def bn(i): return V("bn", i, scalarize)
    def bz(i): return V("bz", i, scalarize)

    # The distinct (n, kz) combinations s_list actually asks for, sharing one
    # bn(n)*zq^kz product across every bra monomial that needs it (level 1,
    # coordinator review item 2). kz == 0 needs no separate slot: bn(n) IS
    # that product already (zq^0 = 1).
    nk_pairs = sorted({(n, kz) for _g, n, kz, _tt in f.s_list if kz})
    nk_index = {p: i for i, p in enumerate(nk_pairs)}
    NBZ = len(nk_pairs)

    o = []
    o.append(emit_frame())
    o.append(f"""
         ! === per-quartet integral core (FRAME.md sections 2-4) ===
         !
         ! Ket outside, bra inside (README_TERCO.md section 2.4): the level-1
         ! accumulators below are per BRA primitive and do not depend on which
         ! Cartesian component is being built, so putting the bra loop
         ! innermost is what keeps a (pp|pp) inner loop to a few dozen
         ! multiply-adds. This is the SAME order libfint's own rotaxis
         ! emitter uses (scripts/rotaxis_mmd/emit.py:emit_kernel).
         !
{chr(10).join(f'         {racc(i)} = 0.0_dp' for i in range(1, NR + 1))}
         do kq = offcd + 1, offcd + ncd
            eta = pp_p(kq)
            ! PRIMITIVE-QUARTET PRESCREEN, adapted from the HGP kernel's `bnd`
            ! bound to this loop's swapped order: the pair lists are sorted by
            ! |c|/p descending (build_pairs_hgp), so `bnd` bounds every
            ! remaining bra primitive using THIS ket's weight, and testing it
            ! against the best bra (offab + 1) before entering the bra loop
            ! at all, then EXITing that loop the first time a weaker bra fails
            ! it, is the same two-part test the HGP kernel runs -- just with
            ! bra and ket swapped, because here the ket is the outer loop.
            bnd = TWO_PI_2_5*abs(pp_cs(kq))/(eta*sqrt(eta))
            if (bnd*abs(pp_cs(offab + 1))/pp_p(offab + 1) <= pcut) cycle
            iq = 0.5_dp/eta
            ! yc, yd = c/q, d/q (FRAME.md section 3) are recovered by
            ! projecting the ket pair centre Q = pp_r(kq,:) onto CD, since
            ! terco's pp_p carries only q = c + d, not c and d separately.
            ! Q - C = yd*(D - C), so |Q - C|/rcd = yd exactly whenever
            ! rcd > 0. When C == D (rcd == 0, routine for two shells on one
            ! atom) that ratio is 0/0; yd is then set to an arbitrary finite
            ! value instead; it is provably safe because FRAME.md section 0's
            ! derivation shows every term carrying a yc or yd power also
            ! carries a matching Rs or Rc power from the SAME Hermite
            ! coefficient (they are only ever introduced together), and
            ! Rs = Rc = 0 exactly whenever C == D -- so the term vanishes via
            ! the constant factor regardless of what yd is, PROVIDED yd is
            ! finite: 0.0*NaN is NaN, and the whole point of this branch is to
            ! keep that finite so a degenerate frame does not poison the sum.
            if (rcd > 0.0_dp) then
               yd = sqrt((pp_r(kq,1) - pp_ra(kq,1))**2 + (pp_r(kq,2) - pp_ra(kq,2))**2 &
                       + (pp_r(kq,3) - pp_ra(kq,3))**2)/rcd
            else
               yd = 0.5_dp
            end if
            yc = 1.0_dp - yd
            ax = cl1 + yd*Rs
            az = cl3 + yd*Rc
            kcdv = TWO_PI_2_5*pp_cs(kq)/eta
{chr(10).join(f'            {sacc(i)} = 0.0_dp' for i in range(1, NS + 1))}
            do kp = offab + 1, offab + nab
               zeta = pp_p(kp)
               if (bnd*abs(pp_cs(kp))/zeta <= pcut) exit
               zpe = zeta + eta
               ! The exact (ss|ss)-type bound: bit-for-bit the HGP kernel's
               ! own `pref`, since kab*kcd*B_0 regroups to the identical
               ! product (FRAME.md section 4) -- same cutoff, same meaning.
               pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_cs(kp)*pp_cs(kq)
               if (abs(pref) <= pcut) cycle
               ip = 0.5_dp/zeta
               za = zhx*(pp_r(kp,1) - pp_ra(kp,1)) + zhy*(pp_r(kp,2) - pp_ra(kp,2)) &
                  + zhz*(pp_r(kp,3) - pp_ra(kp,3))
               zb = za - rab
               zq = za - az
               alpha = zeta*eta/zpe
               tval = alpha*(ax*ax + cy*cy + zq*zq)

   {gp.emit_boys(L)}

               ! kabv is folded into the Boys ladder's overall scale HERE,
               ! once per bra primitive, rather than multiplied into every
               ! one of the {NS} level-1 accumulators below (coordinator
               ! review item 2): bn(n) = kabv * (-2 alpha)^n F_n(T) /
               ! sqrt(zpe) is exactly kabv*B_n, so every accumulation that
               ! used to read "kabv*...*bn(n)..." now reads "...*bn(n)..."
               ! for the identical value.
               kabv = pp_cs(kp)/zeta
               w = kabv/sqrt(zpe)""")
    o.append(f"               {bn(0)} = w*f(0)")
    m2a_emitted = False
    for n in range(1, L + 1):
        if not m2a_emitted:
            o.append("               m2a = -2.0_dp*alpha")
            m2a_emitted = True
        o.append("               w = w*m2a")
        o.append(f"               {bn(n)} = w*f({n})")
    for k in range(2, KMAX + 1):
        prev = "zq" if k == 2 else f"zq{k - 1}"
        o.append(f"               zq{k} = {prev}*zq")
    # The (n, kz) products several bra monomials share, precomputed once each
    # (coordinator review item 2's second half) rather than recomputed -- and
    # multiplied -- inside every s_list entry that needs it.
    zqpow = {1: "zq"}
    for k in range(2, KMAX + 1):
        zqpow[k] = f"zq{k}"
    for (n, kz), idx in nk_index.items():
        o.append(f"               {bz(idx + 1)} = {bn(n)}*{zqpow[kz]}")
    for gi, g in enumerate(f.g_list):
        e = mono_expr(["ip", "za", "zb"], g[:3])
        o.append(f"               {gb(gi + 1)} = {e if e else '1.0_dp'}")
    for si, (g, n, kz, _tt) in enumerate(f.s_list):
        parts = []
        gi = f.g_index[g]
        if any(g):
            parts.append(gb(gi + 1))
        parts.append(bz(nk_index[(n, kz)] + 1) if kz else bn(n))
        o.append(f"               {sacc(si + 1)} = {sacc(si + 1)} + {'*'.join(parts)}")
    o.append("            end do   ! kp (bra primitive)")
    o.append("")
    # kcdv is folded into every ket monomial HERE, once per ket primitive
    # (coordinator review item 1), instead of scaling the {NR} level-2
    # accumulators by it after the fact -- which is what the removed `rtmp`
    # array existed to stage. kmv(i) is now kcdv itself when the monomial is
    # trivial (every ket exponent 0), so every r-entry below can multiply by
    # kmv(km_id) unconditionally; the old kme-is-None special case (skip the
    # multiply entirely for a trivial monomial) is gone because there is no
    # longer a bare, un-scaled monomial value to special-case around.
    for ki, km in enumerate(f.km_list):
        e = mono_expr(["iq", "yc", "yd", "ax"], km)
        rhs = f"kcdv*{e}" if e else "kcdv"
        o.append(f"            {kmv(ki + 1)} = {rhs}")
    for ri, (rterms, _tt) in enumerate(f.r_list):
        bykm = {}
        for (km_id, s_id), c in rterms:
            bykm.setdefault(km_id, []).append((c, s_id))
        terms = [racc(ri + 1)]
        for km_id in sorted(bykm):
            inner = [coef_term(c, sacc(s_id + 1)) for c, s_id in bykm[km_id]]
            kme = kmv(km_id + 1)
            if len(inner) == 1:
                c, s_id = bykm[km_id][0]
                terms.append(coef_term(c, f"{kme}*{sacc(s_id + 1)}"))
            else:
                body = " + ".join(inner).replace("+ -", "- ")
                terms.append(f"{kme}*({body})")
        o.append(wrap(racc(ri + 1), terms))
    o.append("         end do   ! kq (ket primitive)")
    o.append("")
    o.append("         ! === quartet assembly, straight into vbuf (level 3) ===")
    for ci, cm in enumerate(f.cm_list):
        e = mono_expr(["Rs", "Rc", "cy"], cm)
        o.append(f"         {cmv(ci + 1)} = {e if e else '1.0_dp'}")
    for comp_i, terms in enumerate(f.assembly):
        bycm = {}
        for scale, cm_id, r_id in terms:
            bycm.setdefault(cm_id, []).append((scale, r_id))
        parts = []
        for cm_id in sorted(bycm):
            inner = [coef_term(c, racc(r_id + 1)) for c, r_id in bycm[cm_id]]
            cme = None if not any(f.cm_list[cm_id]) else cmv(cm_id + 1)
            if len(inner) == 1:
                c, r_id = bycm[cm_id][0]
                base = racc(r_id + 1) if cme is None else f"{cme}*{racc(r_id + 1)}"
                parts.append(coef_term(c, base))
            else:
                body = " + ".join(inner).replace("+ -", "- ")
                parts.append(f"({body})" if cme is None else f"{cme}*({body})")
        o.append(wrap(f"vbuf({comp_i + 1})", parts))
    o.append("")
    o.append("         ! === rotate every p slot back to the lab frame (FRAME.md section 5) ===")
    rot_code = emit_rotation(kinds, na, nb, nc, nd)
    if rot_code:
        o.append(rot_code)
    else:
        o.append("         ! every slot is s: the identity, nothing to rotate.")
    return "\n".join(o), dict(NS=NS, NR=NR, NG=NG, NKM=NKM, NCM=NCM, NCOMP=NCOMP, L=L, KMAX=KMAX,
                              NBZ=NBZ, na=na, nb=nb, nc=nc, nd=nd)


def decl_temp(base, lo, hi, scalarize):
    """Declaration for one of the small per-quartet temporaries, `lo` to `hi`
    inclusive: one array (default) or `hi - lo + 1` individually named
    scalars (the register-promotion comparison of coordinator review item
    3). Returns None for an empty range (e.g. no class ever needs `bz` when
    KMAX < 2, so nothing should be declared)."""
    if hi < lo:
        return None
    if not scalarize:
        if lo == 1:
            return f"      real(dp) :: {base}({hi})"
        return f"      real(dp) :: {base}({lo}:{hi})"
    # One declaration statement, continued before MAX_LINE: a fully
    # scalarized pppp needs 143 names for racc alone, which is nowhere near
    # free-form Fortran's line-length limit on one line.
    names = [f"{base}_{i}" for i in range(lo, hi + 1)]
    lines = []
    cur, cur_has_name = "      real(dp) :: ", False
    for i, nm in enumerate(names):
        tok = nm + ("," if i < len(names) - 1 else "")
        if cur_has_name and len(cur) + len(tok) + 1 > MAX_LINE:
            lines.append(cur.rstrip() + " &")
            cur, cur_has_name = "         ", False
        cur += tok + " "
        cur_has_name = True
    lines.append(cur.rstrip())
    return "\n".join(lines)


def emit_decls(sizes, scalarize=False):
    NS, NR, NG = sizes["NS"], sizes["NR"], sizes["NG"]
    NKM, NCM, NCOMP = sizes["NKM"], sizes["NCM"], sizes["NCOMP"]
    NBZ = sizes["NBZ"]
    L, na, nb, nc, nd = sizes["L"], sizes["na"], sizes["nb"], sizes["nc"], sizes["nd"]
    noframe = sizes.get("noframe", False)
    o = []
    o.append("      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl")
    o.append("      integer(int64) :: nsa, u, kx")
    o.append("      real(dp) :: qcut")
    o.append("      integer :: keyab, keycd, offab, offcd, nab, ncd")
    o.append("      logical :: same_ab, same_cd, same_pair")
    o.append("      logical :: dij, dkl, dpq")
    o.append("      integer :: mui, nuj, lamk, sigl")
    o.append("      integer :: kp, kq, idens, ia, ib, ic, id, idx")
    o.append("      real(dp) :: sc, wq")
    if noframe:
        # No rotated frame at all (coordinator review item A / needs_frame):
        # ssss's T is formed directly from the lab-frame separation.
        o.append("      real(dp) :: pqx, pqy, pqz")
        o.append("      real(dp) :: zeta, eta, zpe, alpha, tval, bnd, pref")
        o.append("      real(dp) :: kabv, kcdv, w")
    else:
        o.append("      real(dp) :: abx, aby, abz, cdx, cdy, cdz, rab, rcd")
        o.append("      real(dp) :: zhx, zhy, zhz, xhx, xhy, xhz, yhx, yhy, yhz")
        o.append("      real(dp) :: wx, wy, wz, wn, dotcz, rot(3, 3)")
        o.append("      real(dp) :: cl1, cl2, cl3, dl1, dl3, Rs, Rc, cy")
        o.append("      real(dp) :: zeta, eta, zpe, alpha, tval, bnd, pref")
        o.append("      real(dp) :: ip, za, zb, zq, iq, yc, yd, ax, az, kabv, kcdv, w, m2a, t1, t2, t3")
    if not noframe and sizes["KMAX"] >= 2:
        o.append("      real(dp) :: " + ", ".join(f"zq{k}" for k in range(2, sizes["KMAX"] + 1)))
    o.append("      real(dp) :: f(0:BOYS_MMAX)")
    o.append("      integer  :: bi, bj, bbase")
    o.append("      real(dp) :: bx, bx2, b0, b1, b2, btt, bet")
    # sacc, racc, kmv, cmv, gb, bn, bz: the six small per-quartet temporaries
    # the fold in emit_core shrank. Array by default; `scalarize` names each
    # element instead, to compare against ptxas's own register allocation
    # (coordinator review item 3; gen_rotaxis.py --scalarize picks this up).
    for line in (decl_temp("gb", 1, NG, scalarize),
                 decl_temp("sacc", 1, NS, scalarize),
                 decl_temp("kmv", 1, NKM, scalarize),
                 decl_temp("racc", 1, NR, scalarize),
                 decl_temp("cmv", 1, NCM, scalarize),
                 decl_temp("bn", 0, L, scalarize),
                 decl_temp("bz", 1, NBZ, scalarize)):
        if line:
            o.append(line)
    o.append(f"      real(dp) :: vbuf({NCOMP})")
    o.append(f"      real(dp) :: jab({na*nb}), jcd({nc*nd}), kac({na*nc})")
    o.append(f"      real(dp) :: kad({na*nd}), kbc({nb*nc}), kbd({nb*nd})")
    o.append(f"      real(dp) :: dab({na*nb}), dcd({nc*nd}), dac({na*nc})")
    o.append(f"      real(dp) :: dad({na*nd}), dbc({nb*nc}), dbd({nb*nd})")
    return "\n".join(o)


def emit_class_module(kinds, cmd, scalarize=False, bra_fastest=True):
    tag = tag_of(kinds)
    driver, item = f"ras{tag}", f"rasi{tag}"
    mod = f"trc_ra_k{tag}"
    core, sizes = emit_core(kinds, scalarize=scalarize)
    na, nb, nc, nd = sizes["na"], sizes["nb"], sizes["nc"], sizes["nd"]

    body = []
    body.append(gp.emit_scalar_driver(driver, item))
    body.append("")
    body.append(gp.emit_item_signature(item))
    body.append(emit_decls(sizes, scalarize=scalarize))
    body.append("")
    # `bra_fastest=False` (bra-uniform, ket fastest across a warp) was tried
    # as coordinator review item B: it matches rotaxis's OWN ket-outside/
    # bra-inside primitive nesting instead of HGP's, so in principle the
    # innermost (bra-primitive) loop should broadcast across the warp the
    # way HGP's innermost (ket-primitive) loop does under its own decode.
    # MEASURED, not assumed, on a 20-water 6-31G cluster with per-class nsys
    # timing (both contracted and fully uncontracted): it is not a net win.
    # psps (1010) improved in both regimes (~1.0x, was ~1.04-1.10x), but
    # psss, ppss and ppps got WORSE, some substantially (ppss 1.05-1.11x ->
    # 1.23-1.50x). The likely mechanism: the "one warp = one bra pair"
    # assumption needs the ket-pair count per segment (the modulus, sNB) to
    # be at least ~32; classes whose ket side is trivial (an all-s shell
    # pair, as in psss's and ppss's kets here) tend to have FEW distinct ket
    # pairs per segment, so a warp spans several bra pairs after all and the
    # bra loads it was meant to keep uniform diverge again -- while the ket
    # loads that used to broadcast under HGP's own decode (still available
    # via `bra_fastest=True`) no longer do either. Shipped default is
    # therefore `bra_fastest=True`, HGP's own decode, for every rotaxis
    # class; `emit_segment_decode`'s (kept, correct, verified) mirrored form
    # is available for a future per-segment or per-class choice this
    # milestone did not have time to build and validate properly.
    body.append(gp.emit_segment_decode(bra_fastest=bra_fastest))
    body.append(gp.emit_screens())
    body.append(gp.emit_pair_offsets(False))
    body.append(core)
    body.append("")
    body.append(gp.emit_scalar_decode())
    body.append(gp.emit_digestion(na, nb, nc, nd))
    body.append(f"   end subroutine {item}")
    body_text = "\n".join(body)

    d = derive.Derivation(*kinds)
    f = factor.Factorised(d)
    header = f"""!
! Rotated-axis MD class ({kinds[0]}{kinds[1]}|{kinds[2]}{kinds[3]}) kernel.
!
! GENERATED by {cmd} -- do not edit.
!
! {f.summary()}
!
! One module per class, matching src/generated/trc_pc_k{tag}.F90: the item
! routine stays in the same file as its driver so it inlines (the whole
! reason the HGP kernels do this too, see docs/GENERATORS.md), and the
! dispatcher is host code free to cross modules.
!
module {mod}
   use, intrinsic :: iso_fortran_env, only: int64
   use trc_boys, only: dp, boys_eval, BOYS_MMAX, boys_table, &
                         BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, &
                         BOYS_DT, BOYS_DTINV
   use trc_tables, only: LMAX
   implicit none
   private
   public :: {driver}

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp

contains
{body_text}
end module {mod}
"""
    return header, driver, mod


def emit_dispatcher(cmd):
    """`trc_ra_kernels`: `ra_dispatch` (same arguments as `pc_dispatch`) and
    `ra_has_class`, which -- per the "supported() must ask the dispatcher"
    landmine in libfint's own README_TERCO.md section 5 -- answers from this
    module's own case list rather than from an angular-momentum bound.

    The case labels are CLASS_RADIX-symbolic expressions, not baked decimal
    literals: `use trc_pc_kernels, only: CLASS_RADIX` and build each label as
    `((la*CLASS_RADIX+lb)*CLASS_RADIX+lc)*CLASS_RADIX+ld`, the identical
    formula trc_binkernel.F90 uses to build the key it passes in. That means
    there is no baked constant here that generation time and compile time
    could ever disagree on -- a stronger guarantee than the `#if TRC_LMAX !=
    N` guard the other generators need, because rotaxis's class set does not
    depend on TRC_LMAX at all (every one of these 16 classes exists at every
    TRC_LMAX >= 1). The one thing that CAN go stale is compiling this file
    into a build that never defines TRC_PERCLASS -- there is then no
    trc_pc_kernels module to import CLASS_RADIX from -- so that is what the
    guard below catches, turning that mismatch into a compile error instead
    of a link error with a less obvious message."""
    mods = [(tag_of(k), k) for k in ALL_CLASSES]
    o = ["!",
         "! Rotated-axis kernel dispatcher, matching src/generated/trc_pc_kernels.F90.",
         "!",
         f"! GENERATED by {cmd} -- do not edit.",
         "!",
         "#ifndef TRC_PERCLASS",
         '#error "trc_ra_kernels.F90 needs TRC_PERCLASS (it imports CLASS_RADIX from trc_pc_kernels)."',
         "#endif",
         "module trc_ra_kernels",
         "   use, intrinsic :: iso_fortran_env, only: int64",
         "   use trc_boys, only: dp",
         "   use trc_pc_kernels, only: CLASS_RADIX"]
    for tag, _k in mods:
        o.append(f"   use trc_ra_k{tag}, only: ras{tag}")
    o.append("   implicit none")
    o.append("   private")
    o.append("   public :: ra_dispatch, ra_has_class")
    o.append("")
    o.append("contains")
    o.append("")
    o.append("   ! Answers from THIS module's own case list (the landmine libfint's own")
    o.append("   ! README_TERCO.md section 5 records: testing an angular-momentum bound")
    o.append("   ! instead let a class that has no kernel here reach `error stop` below).")
    o.append("   pure logical function ra_has_class(code) result(yes)")
    o.append("      integer, intent(in) :: code")
    o.append("      select case (code)")
    labels = []
    for tag, k in mods:
        la, lb, lc, ld = (_KIND_L[x] for x in k)
        labels.append(f"((({la}*CLASS_RADIX + {lb})*CLASS_RADIX + {lc})*CLASS_RADIX + {ld})")
    o.append("      case (" + ", &\n            ".join(", ".join(labels[i:i + 4]) for i in range(0, len(labels), 4)) + ")")
    o.append("         yes = .true.")
    o.append("      case default")
    o.append("         yes = .false.")
    o.append("      end select")
    o.append("   end function ra_has_class")
    o.append("")
    o.append("   subroutine ra_dispatch(key, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sD, &")
    o.append("                          npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &")
    o.append("                          pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, pp_cs, &")
    o.append("                          ndens, dmat, jmat, rank, nranks)")
    o.append("      integer, intent(in) :: key")
    o.append(gp.PROLOGUE)
    o.append("      integer, intent(in) :: rank, nranks")
    o.append("      select case (key)")
    for i, (tag, k) in enumerate(mods):
        la, lb, lc, ld = (_KIND_L[x] for x in k)
        label = f"((({la}*CLASS_RADIX + {lb})*CLASS_RADIX + {lc})*CLASS_RADIX + {ld})"
        o.append(f"""      case ({label}); call ras{tag}(lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sD, &
                          npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                          pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, pp_cs, ndens, dmat, jmat, rank, nranks)""")
    o.append("      case default")
    o.append("         error stop 'trc_ra_kernels: no rotaxis kernel for this class'")
    o.append("      end select")
    o.append("   end subroutine ra_dispatch")
    o.append("end module trc_ra_kernels")
    return "\n".join(o) + "\n"


def emit_files(cmd, scalarize=False, bra_fastest=True):
    """Return {filename: text} for every file this emitter writes, mirroring
    gen_perclass.py --split's own return shape."""
    files = {}
    for kinds in ALL_CLASSES:
        header, _driver, mod = emit_class_module(kinds, cmd, scalarize=scalarize,
                                                 bra_fastest=bra_fastest)
        files[f"{mod}.F90"] = header
    files["trc_ra_kernels.F90"] = emit_dispatcher(cmd)
    return files
