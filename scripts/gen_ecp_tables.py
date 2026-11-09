#!/usr/bin/env python3
"""
Emit the Cartesian-to-pure-spherical-harmonic coefficient blocks the ECP
angular integrals need, l = 0..ECP_ANG_LMAX.

WHY A SEPARATE TABLE FROM gen_sph.py
-------------------------------------
gen_sph.py's SPH_T stops at l = 4, which is exactly what the orbital AO space
needs (TERCO_LMAX never goes past f). The ECP angular code
(src/eri/trc_ecp_angular.F90, ported from PySCF's nr_ecp.c via libfint) needs
the SAME coefficients -- the "sp"-normalised real solid harmonic, written in
Cartesian monomials, libcint's own convention -- but at l = li + lc, an
ORBITAL angular momentum plus an ECP CHANNEL's, which reaches past 4 for a
d or f orbital next to a heavy element's l=3 local channel. Reusing SPH_T's
table would mean widening it for a use SPH_T's own consumer (trc_sph.F90)
never needs; a second, purpose-built table keeps the orbital AO path's table
untouched and the ECP path's own ceiling (`ECP_ANG_LMAX`) independent of it.

THE COEFFICIENTS ARE THE SAME CLOSED FORM
-------------------------------------------
`real_solid_harmonic_coef` (imported from gen_sph.py, not re-derived: a second
copy of a hand-transcribed closed form is exactly the failure mode
cartorder.py's own docstring warns about) is checked there against
`pyscf.gto.mole.cart2sph` to 1e-13 for l = 0..4. PySCF's own ECP code
(nr_ecp.c's `ECPcart2sph`/`ECPsph2cart`, which cint_ecp_ang.f90's
`cint_c2s_bra_sph`/`cint_s2c_bra_sph` port) is generated from the identical
convention, so the SAME table, read forward (Cartesian -> spherical, 1 for
every l it exports) and then backward with its own transpose, reproduces
exactly the c2s-then-s2c round trip `ecp_ang_nuc_in_cart` needs for l >= 2:
both directions read one (ncart, nsph) block, so a round trip through it is
mathematically the projector onto the l-pure harmonic subspace regardless of
which orthonormal basis of that subspace the table's columns happen to be --
and this table's columns ARE that basis, verified against pyscf at l <= 4
where a check is possible. l = 0 and l = 1 are not built from this table at
all: `trc_ecp_angular.F90` keeps PySCF's own hand-written constants for
those two (1/(4 pi), r*3/(4 pi)), exactly as nr_ecp.c does, so the table
starts contributing only from l = 2.

l = 0..1 ARE STILL EMITTED (as the identity), for the single OTHER caller
that reads this table directly rather than through the round trip:
`ecp_type2_facs_ang`'s own one-sided forward transform of the ECP channel's
own angular momentum `lc`, which libfint special-cases at lc = 0, 1 with the
same two constants and falls through to the general c2s machinery
otherwise. `trc_ecp_angular.F90` keeps that same special case, so this
table's l = 0, 1 slices are never actually read; they are emitted anyway so
every l up to ECP_ANG_LMAX is uniformly indexable and a future caller does
not have to special-case the table's own shape.

USAGE
    gen_ecp_tables.py -o src/generated/trc_ecp_tables.F90
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cartorder
from gen_sph import real_solid_harmonic_coef

ECP_ANG_LMAX = 8
# li (<= 3, f) + lc (<= 4, one below libfint's ECP_LMAX = 5 local-channel
# ceiling: the local channel is never itself a *projected* channel) reaches
# 7; type 1's li + lj reaches 6. 8 leaves one full unit of headroom past
# both, matching ECP_DESIGN.md's "extend past l=4 to 8".


def _c2s_block(l):
    """(ncart, 2l+1) list of lists, Cartesian rows in lexicographic order,
    columns m = -l..l ascending for l >= 2, identity for l <= 1 (never read
    there -- see module docstring)."""
    comps = cartorder.components(l, "lexicographic")
    ncart = len(comps)
    if l <= 1:
        return [[1.0 if i == j else 0.0 for j in range(ncart)] for i in range(ncart)]
    cols = list(range(-l, l + 1))
    return [[real_solid_harmonic_coef(lx, ly, lz, m) for m in cols] for (lx, ly, lz) in comps]


def _fortran_double(x):
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


def emit(cmd):
    ncart_max = (ECP_ANG_LMAX + 1) * (ECP_ANG_LMAX + 2) // 2
    nsph_max = 2 * ECP_ANG_LMAX + 1
    flat = []
    for l in range(ECP_ANG_LMAX + 1):
        blk = _c2s_block(l)
        ncart = len(blk)
        nsph = 2 * l + 1
        for isph in range(nsph_max):
            for icart in range(ncart_max):
                if icart < ncart and isph < nsph:
                    flat.append(blk[icart][isph])
                else:
                    flat.append(0.0)
    data = _wrapped(flat)
    return f"""!
! Cartesian-to-pure-spherical-harmonic coefficient blocks, l = 0..{ECP_ANG_LMAX},
! for the ECP angular integrals.
!
! GENERATED by {cmd} -- do not edit. See src/eri/trc_ecp_angular.F90 for how
! this table is used (the forward-and-back round trip that projects a
! Cartesian monomial onto its l-pure harmonic part) and this script's own
! module docstring for why it exists separately from trc_sph_tables.F90.
!
module trc_ecp_tables
   use pic_types, only: dp
   implicit none
   private
   public :: ECP_ANG_LMAX, ECP_ANG_NCART_MAX, ECP_ANG_NSPH_MAX, ECP_ANG_C2S

   integer, parameter :: ECP_ANG_LMAX = {ECP_ANG_LMAX}
   integer, parameter :: ECP_ANG_NCART_MAX = {ncart_max}   !! (ECP_ANG_LMAX+1)(ECP_ANG_LMAX+2)/2
   integer, parameter :: ECP_ANG_NSPH_MAX  = {nsph_max}    !! 2*ECP_ANG_LMAX+1

   real(dp), parameter :: ECP_ANG_C2S(ECP_ANG_NCART_MAX, ECP_ANG_NSPH_MAX, 0:{ECP_ANG_LMAX}) = reshape( &
      [ &
         {data} &
      ], [ECP_ANG_NCART_MAX, ECP_ANG_NSPH_MAX, {ECP_ANG_LMAX + 1}])
   !! ECP_ANG_C2S(icart, isph, l): row `icart` (1-based, lexicographic
   !! Cartesian order, independent of TRC_CART_ORDER -- the ECP integrals
   !! are computed entirely in the internal, always-lexicographic AO space,
   !! same as every other kernel; see AGENTS.md/pr-review's note that
   !! TRC_CART_ORDER only ever touches an output layout) and column `isph`
   !! (1-based; m = isph - l - 1, ascending). Only the leading
   !! ncart(l)-by-(2l+1) corner of each l-slice is meaningful; the padding
   !! is exact zero.

end module trc_ecp_tables
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="src/generated/trc_ecp_tables.F90")
    args = ap.parse_args()
    cmd = "scripts/gen_ecp_tables.py"
    with open(args.out, "w") as fh:
        fh.write(emit(cmd))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
