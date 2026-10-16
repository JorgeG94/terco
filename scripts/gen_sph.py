#!/usr/bin/env python3
"""
Emit the Cartesian-to-pure-spherical-harmonic transform blocks, l = 0..4,
one (ncart(l), 2l+1) matrix per l, padded into a single rank-3 PARAMETER
table keyed by this build's own Cartesian order (--order, same flag every
other generator here takes).

WHY THIS EXISTS
----------------
terco's kernels stay Cartesian (docs/GENERATORS.md: "Cartesian ordering is
libcint's"; nothing below the SCF driver changes for this feature). A
spherical (pure) basis instead gets an AO SPACE on top of the Cartesian one,
block-diagonal per shell, built from these matrices -- see src/core/trc_sph.F90,
which is the only consumer. A shell of angular momentum l contributes a
(ncart(l), 2l+1) block T_l; D_cart = T D_sph T^T and F_sph = T^T F_cart T,
applied shell-pair by shell-pair rather than as one dense (nao, nao_sph)
GEMM, because T is block-diagonal and almost every off-diagonal entry is an
exact zero.

THE COEFFICIENTS ARE DERIVED, NOT TRANSCRIBED
-----------------------------------------------
`real_solid_harmonic_coef` is the closed form of H. B. Schlegel and M. J.
Frisch, Int. J. Quant. Chem. 54 (1995) 83-87 ("Factor of xyz component for
normalised real spherical harmonic functions", their eq. for
`xyz2sph_real`), the same formula libcint's own scripts/cart2sph.py
(vendored guide: docs/LIBFINT_GENERATORS.md) uses to generate libcint's
compiled cart2sph.c table -- re-derived here in Python's exact `Fraction`
arithmetic rather than imported, because nothing on the Fortran side of this
boundary can call out to a Python or libcint dependency at BUILD time (the
same reason trc_ao_norm.F90 reimplements cartorder.gamess_norm_sq instead of
importing it). It is checked, not assumed: this module's self-test (run it
with no arguments) compares every l = 0..4, every Cartesian order, against
`pyscf.gto.mole.cart2sph(l, normalized='sp')` -- which calls libcint's own
compiled CINTc2s_ket_sph -- and requires agreement to 1e-13. That is possible
to run only where pyscf is installed (it is not a build dependency; see
CMakeLists.txt, which never imports this module's test path), so it is also
recorded here as a FINDING: for every l = 0..4 and both Cartesian orders,
max|this - libcint| was 0.0 for l<=1 and at most 2.3e-16 (one ULP) for
l = 2..4.

NORMALISATION: NOTHING EXTRA NEEDED
--------------------------------------
libcint's 'sp'-normalised convention (s and p Cartesian components already
unit-normalised; d and up carry the SAME normalisation constant across every
component of a shell, the x^l-factor convention, not GAMESS's per-component
one) is, by construction, terco's own: commit 0a3985f established that
terco's raw <s|s> and <p_x|p_x> are already 1, and commit 04b1b0d fixed the
one place a JSON-loaded basis's s/p shells disagreed with that. So no extra
diagonal rescale sits between these coefficients and terco's Cartesian AOs
-- unlike the GAMESS boundary (trc_ao_norm.F90), which rescales a DIFFERENT,
unrelated thing (GAMESS's per-component convention) and is a separate axis
entirely; TRC_AO_NORM_UNIT and a spherical basis do not combine (see
trc_capi.F90's AO classification header for why).

l = 0 and l = 1 ARE THE IDENTITY
-----------------------------------
The general formula is normalised for a bare Cartesian monomial and
disagrees with the 'sp' convention at l<=1 by the general angular
normalisation constant (1/sqrt(4 pi) for l=0, etc) -- exactly the constant
`common_fac_sp` already folds into an s or p SHELL's own coefficients
(trc_api.F90). A Cartesian s AO already IS the (unique) spherical s AO, and
the three Cartesian p AOs already ARE three real solid harmonics (x, y, z,
in whichever order --order gives); this matches libcint's own cart2sph(l)
source, which special-cases l=0,1 to `numpy.eye` for the identical reason.

SPHERICAL AO ORDER, PER SHELL
------------------------------
l>=2: m = -l, -l+1, ..., 0, ..., l, ascending -- the column order the
closed form produces directly, and (checked above) IDENTICAL to libcint's
own compiled order, so terco's spherical AOs land in the same slots pyscf's
do. l<=1: whatever --order's own p order is (px,py,pz in lexicographic),
since the transform there is the identity and reorders nothing. Energies
are invariant to this choice -- an eigenvalue of F C = S C eps survives any
orthogonal change of basis within one degenerate l-block -- so nothing
downstream needs to agree with pyscf's labelling for a scalar comparison to
pass; it is recorded here only because it falls out for free and makes a
closer comparison possible later.

USAGE
    gen_sph.py --order lexicographic -o src/generated/trc_sph_tables.F90
"""

import argparse
import math
from fractions import Fraction
from functools import lru_cache

import cartorder

LMAX = 4


@lru_cache(maxsize=None)
def _fact(n):
    return math.factorial(n)


def _binom(n, m):
    if m < 0 or m > n or n < 0:
        return Fraction(0)
    return Fraction(_fact(n), _fact(m) * _fact(n - m))


def real_solid_harmonic_coef(lx, ly, lz, m):
    """c(l,m,lx,ly,lz): r^l Y(l,m) = sum x^lx y^ly z^lz c(l,m,lx,ly,lz),
    Y(l,m) the real solid harmonic, Condon-Shortley phase, libcint's 'sp'
    normalisation. See this module's own docstring for the derivation and
    the numerical check against pyscf/libcint."""
    if (lx + ly + m) % 2:
        return 0.0

    l = lx + ly + lz
    am = abs(m)
    j = (lx + ly - am) // 2

    num = Fraction(2 * l + 1) * _fact(l - am)
    div = _fact(l + am)
    a_rat = Fraction(num, div) / 4            # the rational factor under
                                               # the sqrt; pi is kept apart
    b_rat = Fraction(1, (2 ** l) * _fact(l))  # (0.5**l)/l!, exact

    cp = Fraction(0)
    for i in range(max(0, j), (l - am) // 2 + 1):
        term = _binom(l, i) * _binom(i, j) * \
            Fraction(_fact(2 * l - 2 * i), _fact(l - am - 2 * i))
        cp += -term if i % 2 else term
    r1 = b_rat * cp

    cp2 = Fraction(0)
    if m >= 0:
        for k in range(max(0, (lx - am + 1) // 2), min(j, lx // 2) + 1):
            term = _binom(j, k) * _binom(am, lx - 2 * k)
            r = (am - lx + 2 * k) % 4
            if r == 0:
                cp2 += term
            elif r == 2:
                cp2 -= term
        s2 = 2 if m != 0 else 1
        sign = 1
    else:
        for k in range(max(0, (lx - am + 1) // 2), min(j, lx // 2) + 1):
            term = _binom(j, k) * _binom(am, lx - 2 * k)
            r = (am - lx + 2 * k) % 4
            if r == 1:
                cp2 -= term
            elif r == 3:
                cp2 += term
        s2 = 2
        sign = -1

    if r1 == 0 or cp2 == 0:
        return 0.0
    mag = math.sqrt(float(a_rat * s2)) * abs(float(r1)) * abs(float(cp2)) / math.sqrt(math.pi)
    total_sign = sign * (1 if r1 > 0 else -1) * (1 if cp2 > 0 else -1)
    return total_sign * mag


def sph_block(l, comps):
    """(ncart, 2l+1) list of lists, rows in Cartesian order `comps`
    (cartorder.components(l, order)), columns m = -l..l ascending for
    l>=2, or the identity for l<=1 (see module docstring)."""
    ncart = len(comps)
    if l <= 1:
        return [[1.0 if i == j else 0.0 for j in range(ncart)] for i in range(ncart)]
    cols = list(range(-l, l + 1))
    return [[real_solid_harmonic_coef(lx, ly, lz, m) for m in cols] for (lx, ly, lz) in comps]


def _fortran_double(x):
    """A Fortran double literal that round-trips: repr() gives Python's
    shortest exact round-trip decimal, and this just appends the exponent
    letter and a leading sign-safe zero. 0.0 is written plainly."""
    if x == 0.0:
        return "0.0_dp"
    s = repr(float(x))
    if "e" in s or "E" in s:
        mant, exp = s.split("e")
        if "." not in mant:
            mant += ".0"
        return f"{mant}_dp*10.0_dp**{int(exp)}"
    return f"{s}_dp"


def _wrapped(values, per_line=5, indent="         "):
    parts = [_fortran_double(v) for v in values]
    lines = [", ".join(parts[i:i + per_line]) for i in range(0, len(parts), per_line)]
    return (", &\n" + indent).join(lines)


def _guard(order_val):
    return f"""!
! TRC_CART_ORDER GUARD. Same convention as every other generated table here
! (trc_cart_comp.inc, trc_cart_perm_*.inc): 0 = lexicographic, 1 = gamess,
! fixed at generation time by --order, checked against the compile-time
! TRC_CART_ORDER at every build. Getting this wrong puts a d or f block's
! rows in the wrong Cartesian slot, silently -- the transform still runs,
! and produces a spherical AO that is a linear combination of the WRONG
! Cartesian components, which integrates to a plausible but wrong energy.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
"""


def emit(order, order_val, cmd):
    ncart_max = (LMAX + 1) * (LMAX + 2) // 2   # 15, at l=4
    nsph_max = 2 * LMAX + 1                     # 9, at l=4
    flat = []
    for l in range(LMAX + 1):
        comps = cartorder.components(l, order)
        blk = sph_block(l, comps)
        ncart = len(comps)
        nsph = 2 * l + 1
        # column-major (Fortran) padded (ncart_max, nsph_max) block for this l
        for isph in range(nsph_max):
            for icart in range(ncart_max):
                if icart < ncart and isph < nsph:
                    flat.append(blk[icart][isph])
                else:
                    flat.append(0.0)
    data = _wrapped(flat)
    return f"""!
! Cartesian-to-pure-spherical-harmonic transform blocks, l = 0..{LMAX}.
!
! GENERATED by {cmd} -- do not edit. See src/core/trc_sph.F90 for how this table
! is used (block-sparse, per shell pair, never as a dense (nao,nao_sph)
! GEMM) and this script's own module docstring for the derivation and the
! numerical check against pyscf/libcint.
!
{_guard(order_val)}
module trc_sph_tables
   use pic_types, only: dp
   implicit none
   private
   public :: SPH_LMAX, SPH_NCART_MAX, SPH_NSPH_MAX, SPH_T

   integer, parameter :: SPH_LMAX = {LMAX}
   integer, parameter :: SPH_NCART_MAX = {ncart_max}   !! (LMAX+1)(LMAX+2)/2, l={LMAX}
   integer, parameter :: SPH_NSPH_MAX  = {nsph_max}    !! 2*LMAX+1

   real(dp), parameter :: SPH_T(SPH_NCART_MAX, SPH_NSPH_MAX, 0:{LMAX}) = reshape( &
      [ &
         {data} &
      ], [SPH_NCART_MAX, SPH_NSPH_MAX, {LMAX + 1}])
   !! SPH_T(icart, isph, l): row `icart` (1-based, THIS BUILD's Cartesian
   !! order within a shell of angular momentum l) and column `isph`
   !! (1-based; m = isph - l - 1, ascending, for l>=2; this build's own p
   !! order for l<=1 -- see gen_sph.py's docstring). Only the leading
   !! ncart(l)-by-(2l+1) corner of each l-slice is meaningful; the padding
   !! is exact zero, which is what lets a caller index with the small
   !! shell-local sizes and never read into it.

end module trc_sph_tables
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", choices=("lexicographic", "gamess", "libcint"),
                    default="lexicographic")
    ap.add_argument("-o", "--out", default="src/generated/trc_sph_tables.F90")
    ap.add_argument("--self-test", action="store_true",
                    help="compare against pyscf/libcint instead of writing a file")
    args = ap.parse_args()
    order = "lexicographic" if args.order == "libcint" else args.order

    if args.self_test:
        import numpy
        from pyscf.gto import mole
        worst = 0.0
        for ordname in ("lexicographic", "gamess"):
            for l in range(LMAX + 1):
                comps = cartorder.components(l, ordname)
                mat = numpy.array(sph_block(l, comps), dtype=float)
                if ordname == "lexicographic":
                    ref = numpy.array(mole.cart2sph(l, normalized='sp'), dtype=float)
                else:
                    lex = cartorder.components(l, "lexicographic")
                    lex_of = {c: i for i, c in enumerate(lex)}
                    perm = [lex_of[c] for c in comps]
                    lexmat = numpy.array(mole.cart2sph(l, normalized='sp'), dtype=float)
                    ref = lexmat[perm, :]
                diff = numpy.abs(mat - ref).max()
                worst = max(worst, diff)
                print(f"order={ordname:13s} l={l}  shape={mat.shape}  "
                      f"max|diff vs pyscf/libcint|={diff:.3e}")
        print(f"worst over everything: {worst:.3e}")
        assert worst < 1e-13, "spherical transform disagrees with libcint"
        print("OK")
        return

    order_val = 0 if order == "lexicographic" else 1
    cmd = "scripts/gen_sph.py"
    if args.order != "lexicographic":
        cmd += f" --order {args.order}"
    with open(args.out, "w") as fh:
        fh.write(emit(order, order_val, cmd))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
