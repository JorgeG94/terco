#!/usr/bin/env python3
"""
Emit one-electron kernels -- overlap, kinetic, nuclear attraction -- one
straight-line routine per (la, lb) class.

STAGE ONE OF THE OPERATOR GENERALISATION
----------------------------------------
The four-centre generator is specialised to the ERI operator in exactly two
places: the base case ([00|00]^(m) = pref * F_m(T)) and the recurrence
coefficients. Everything else -- the per-class straight-line emission, the
sieve, the include-not-module discipline, the batched `do concurrent` driver --
is operator-independent. This file is the second operator, written to find out
what actually generalises before committing to an abstraction.

What it shows: the SHAPE differs more than the coefficients do.

  Overlap and kinetic FACTORISE. The operator is multiplicative and separable,
  so S = Sx*Sy*Sz with three independent one-dimensional recurrences, and there
  is no auxiliary index and no Boys function at all. Kinetic rides on the same
  1-D tables through
      T_i(a,b) = beta(2b+1) S_i(a,b) - 2 beta^2 S_i(a,b+2)
                 - (1/2) b(b-1) S_i(a,b-2)
  and T = Tx*Sy*Sz + Sx*Ty*Sz + Sx*Sy*Tz. Cheap, and nothing like the ERI.

  Nuclear attraction does NOT factorise -- 1/|r-C| couples the directions --
  so it has the familiar auxiliary index m and Boys function, and looks like a
  two-centre version of the ERI VRR:
      V[a+1_i]^(m) = PA_i V[a]^(m) - PC_i V[a]^(m+1)
                     + (N_i(a)/2zeta)(V[a-1_i]^(m) - V[a-1_i]^(m+1))
  built at b=0 and then transferred by the SAME horizontal recurrence the ERI
  uses, since the HRR is a property of the two Gaussian centres and not of the
  operator between them.

So the generalisation worth building is not "one recurrence with pluggable
constants". It is a small set of recurrence SHAPES -- separable-no-auxiliary,
two-centre-with-auxiliary, four-centre-with-auxiliary -- sharing the emission
machinery and the HRR.

USAGE
    gen_1e.py --lmax 2 -o src/trc_1e_kernels.F90
"""

import argparse

# FIXED, and not lmax + 1.
#
# The dispatch key was la*(ONE_E_LMAX + 1) + lb, which encodes the ceiling it
# was generated at. At ONE_E_LMAX = 2 the keys 0..8 happen to be unique, so
# nothing is wrong today -- but the moment the ceiling rises the encoding
# COLLIDES rather than merely missing: at a radix of 3, (s|f) keys to 3 and so
# does (p|s), so an f shell would not fall through to nothing. It would run the
# (p|s) kernel and return a plausible, wrong overlap.
#
# That is worse than the four-centre version of this bug, which dispatched a
# wrong kernel only across a rebuild. This one is wrong within a single
# consistent build.
#
# 11 matches CLASS_RADIX in gen_perclass.py, which had the same bug.
ONE_E_RADIX = 11


def cart(l):
    """(nx,ny,nz) for one shell, libcint order."""
    return [(lx, ly, l - lx - ly)
            for lx in range(l, -1, -1)
            for ly in range(l - lx, -1, -1)]


def cart_cum(lmax):
    out = []
    for l in range(lmax + 1):
        out.extend(cart(l))
    return {p: i + 1 for i, p in enumerate(out)}


def ncum(l):
    return (l + 1)*(l + 2)*(l + 3)//6


def binom(n, k):
    r = 1
    for i in range(1, k + 1):
        r = r*(n - k + i)//i
    return r


def powexpr(var, n):
    return "*".join([var]*n) if n else ""


# --------------------------------------------------------------------------
# separable part: 1-D overlap tables, and kinetic on top of them
# --------------------------------------------------------------------------

def emit_1d(la, lb, ax, extra=0):
    """1-D overlap table s{ax}(0:la, 0:lb+2), by the standard recurrence.

    Built to lb+2 because the kinetic operator reaches two units higher in b.
    Straight-line with literal indices: every trip count here is known.
    """
    o, w = [], None
    w = o.append
    hb = lb + 2
    la = la + extra          # multipoles need a raised by the operator degree
    w(f"         s{ax}(0,0) = e0{ax}")
    # first column: raise a at b = 0
    for a in range(1, la + 1):
        t = f"pa{ax}*s{ax}({a - 1},0)"
        if a - 1 >= 1:
            t += f" + {float(a - 1)}_dp*oo2z*s{ax}({a - 2},0)"
        w(f"         s{ax}({a},0) = {t}")
    # then raise b for every a
    for b in range(1, hb + 1):
        for a in range(0, la + 1):
            t = f"pb{ax}*s{ax}({a},{b - 1})"
            if a >= 1:
                t += f" + {float(a)}_dp*oo2z*s{ax}({a - 1},{b - 1})"
            if b - 1 >= 1:
                t += f" + {float(b - 1)}_dp*oo2z*s{ax}({a},{b - 2})"
            w(f"         s{ax}({a},{b}) = {t}")
    return "\n".join(o)


def emit_kin_1d(la, lb, ax):
    """t{ax}(0:la,0:lb) = beta(2b+1)S - 2 beta^2 S(a,b+2) - b(b-1)/2 S(a,b-2)."""
    o, w = [], None
    w = o.append
    for b in range(lb + 1):
        for a in range(la + 1):
            t = f"{float(2*b + 1)}_dp*bet*s{ax}({a},{b})" \
                f" - tbb*s{ax}({a},{b + 2})"
            if b >= 2:
                t += f" - {0.5*b*(b - 1)}_dp*s{ax}({a},{b - 2})"
            w(f"         t{ax}({a},{b}) = {t}")
    return "\n".join(o)


# --------------------------------------------------------------------------
# nuclear attraction: two-centre VRR with an auxiliary index, then the HRR
# --------------------------------------------------------------------------

def emit_nuc_vrr(lab, idx):
    """V[a]^(m) at b = 0, for every a with |a| <= lab and every needed m.

    Triangular in the usual way: level m only needs degree up to lab - m.
    """
    o, w = [], None
    w = o.append
    lst = []
    for l in range(lab + 1):
        for p in cart(l):
            lst.append(p)
    info = {}
    for p in lst:
        x, y, z = p
        if x > 0:
            d, dn1 = 0, (x - 1, y, z)
            dn2 = (x - 2, y, z) if x > 1 else None
            n = x
        elif y > 0:
            d, dn1 = 1, (x, y - 1, z)
            dn2 = (x, y - 2, z) if y > 1 else None
            n = y
        elif z > 0:
            d, dn1 = 2, (x, y, z - 1)
            dn2 = (x, y, z - 2) if z > 1 else None
            n = z
        else:
            continue
        info[p] = (d, dn1, dn2, n)

    PA = ["pax", "pay", "paz"]
    PC = ["pcx", "pcy", "pcz"]
    for m in range(lab, -1, -1):
        w(f"            ! --- m = {m} ---")
        w(f"            vn({idx[(0,0,0)]},{m}) = vpref*f({m})")
        for l in range(1, lab - m + 1):
            for p in cart(l):
                d, dn1, dn2, n = info[p]
                t = (f"{PA[d]}*vn({idx[dn1]},{m})"
                     f" - {PC[d]}*vn({idx[dn1]},{m + 1})")
                if dn2 is not None:
                    t += (f" &\n               + {float(n - 1)}_dp*oo2z*"
                          f"(vn({idx[dn2]},{m}) - vn({idx[dn2]},{m + 1}))")
                w(f"            vn({idx[p]},{m}) = {t}")
    return "\n".join(o)


def emit_hrr(la, lb, idx, src, dst):
    """(a|b) from [a+b|0] -- the same horizontal transfer the ERI uses.

    It is a statement about the two Gaussian centres, not about the operator
    between them, so it is reused verbatim.
    """
    o, w = [], None
    w = o.append
    na = len(cart(la))
    for ib, (bx, by, bz) in enumerate(cart(lb)):
        for ia, (axp, ayp, azp) in enumerate(cart(la)):
            terms = []
            for px in range(bx + 1):
                for py in range(by + 1):
                    for pz in range(bz + 1):
                        c = binom(bx, px)*binom(by, py)*binom(bz, pz)
                        j = idx[(axp + px, ayp + py, azp + pz)]
                        fac = [f for f in (powexpr("abx", bx - px),
                                           powexpr("aby", by - py),
                                           powexpr("abz", bz - pz)) if f]
                        head = "" if c == 1 else f"{float(c)}_dp*"
                        body = "*".join(fac)
                        if body:
                            body += "*"
                        terms.append(f"{head}{body}{src}({j})")
            expr = terms[0]
            for t in terms[1:]:
                expr += f" &\n               + {t}"
            w(f"            {dst}({1 + ia + na*ib}) = {expr}")
    return "\n".join(o)


# --------------------------------------------------------------------------

MULT_MAX = 3   # octupole


def emit_multipole(la, lb):
    """Cartesian multipole moments through the octupole, one class.

    SEPARABLE, like overlap and for the same reason: (r-O)^n is a product of
    powers of the three coordinates, so the integral factorises and the whole
    thing rides on the 1-D overlap tables already built here.

    Shifting the origin is a binomial expansion,

        (x - Ox)^p = sum_t C(p,t) (Ax - Ox)^(p-t) (x - Ax)^t

    and (x - Ax)^t times the bra is just the bra with t more units of angular
    momentum, so

        S^(p)(a,b) = sum_t C(p,t) (Ax-Ox)^(p-t) S(a+t, b)

    with S the ordinary 1-D overlap table extended by p in its first index.
    No new recurrence, no auxiliary index, no Boys function.

    mqc wants the FULL Cartesian tensor -- 3, 9 and 27 components, not the 4,
    6 and 10 unique ones -- in libcint's order with the last index fastest.
    The repeats are written out rather than mapped, because a component table
    is another thing to get out of step with the caller.
    """
    o, w = [], None
    w = o.append
    na, nb = len(cart(la)), len(cart(lb))
    axes = "xyz"

    # S^(p) for each direction and each power up to MULT_MAX
    for d, ax in enumerate(axes):
        for pw in range(1, MULT_MAX + 1):
            for b in range(lb + 1):
                for a in range(la + 1):
                    terms = []
                    for t in range(pw + 1):
                        c = binom(pw, t)
                        shift = powexpr(f"o{ax}", pw - t)
                        head = "" if c == 1 else f"{float(c)}_dp*"
                        body = (shift + "*") if shift else ""
                        terms.append(f"{head}{body}s{ax}({a + t},{b})")
                    joiner = " &" + chr(10) + "            + "
                    w(f"         m{ax}{pw}({a},{b}) = " + joiner.join(terms))

    # the components, in libcint order: last index fastest
    slot = 0
    for order in range(1, MULT_MAX + 1):
        for idx in range(3**order):
            # decode the index tuple, last fastest
            t, powers = idx, [0, 0, 0]
            for _ in range(order):
                powers[t % 3] += 1
                t //= 3
            slot += 1
            for ib, (bx, by, bz) in enumerate(cart(lb)):
                for ia, (ax_, ay_, az_) in enumerate(cart(la)):
                    k = 1 + ia + na*ib
                    fac = []
                    for d, ch in enumerate(axes):
                        pw = powers[d]
                        i0 = (ax_, ay_, az_)[d]
                        j0 = (bx, by, bz)[d]
                        fac.append(f"m{ch}{pw}({i0},{j0})" if pw
                                   else f"s{ch}({i0},{j0})")
                    w(f"         mout({k},{slot}) = mout({k},{slot}) + cc*"
                      + "*".join(fac))
    return chr(10).join(o), slot


def emit_class(la, lb, LMAX):
    lab = la + lb
    na, nb = len(cart(la)), len(cart(lb))
    idx = cart_cum(lab)
    nvn = ncum(lab)
    tag = f"{la}{lb}"

    sxyz = "\n".join(emit_1d(la, lb, a) for a in "xyz")
    txyz = "\n".join(emit_kin_1d(la, lb, a) for a in "xyz")

    # overlap and kinetic components, straight from the 1-D tables
    comps = []
    for ib, (bx, by, bz) in enumerate(cart(lb)):
        for ia, (ax, ay, az) in enumerate(cart(la)):
            k = 1 + ia + na*ib
            comps.append(
                f"         sout({k}) = sout({k}) + "
                f"cc*sx({ax},{bx})*sy({ay},{by})*sz({az},{bz})")
            comps.append(
                f"         tout({k}) = tout({k}) + cc*("
                f"tx({ax},{bx})*sy({ay},{by})*sz({az},{bz})"
                f" &\n            + sx({ax},{bx})*ty({ay},{by})*sz({az},{bz})"
                f" &\n            + sx({ax},{bx})*sy({ay},{by})*tz({az},{bz}))")
    comps = "\n".join(comps)

    vrr = emit_nuc_vrr(lab, idx)
    hrr = emit_hrr(la, lb, idx, "vn0", "vtmp")

    return f"""
   !> ({la}|{lb}) one-electron block: overlap, kinetic, nuclear attraction.
   !>
   !> One routine for all three because they share the shell-pair quantities
   !> (zeta, P, the Gaussian product prefactor) and mqc needs all three to form
   !> the core Hamiltonian. Computing them together reads the primitive data
   !> once instead of three times.
   pure subroutine one_e_{tag}(nprim_a, nprim_b, ea, ca, eb, cb, ra, rb, &
                               natm, zatm, ratm, sout, tout, vout)
      !$acc routine seq
      integer,  intent(in)  :: nprim_a, nprim_b, natm
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a)
      real(dp), intent(in)  :: eb(nprim_b), cb(nprim_b)
      real(dp), intent(in)  :: ra(3), rb(3)
      real(dp), intent(in)  :: zatm(natm), ratm(3, natm)
      real(dp), intent(out) :: sout({na*nb}), tout({na*nb}), vout({na*nb})

      integer  :: ip, jp, ic, k
      real(dp) :: al, bet, zeta, oo2z, xi, cc, tbb
      real(dp) :: px, py, pz, pax, pay, paz, pbx, pby, pbz
      real(dp) :: pcx, pcy, pcz, tval, vpref, rab2
      real(dp) :: abx, aby, abz, e0x, e0y, e0z
      real(dp) :: sx(0:{la},0:{lb + 2}), sy(0:{la},0:{lb + 2}), sz(0:{la},0:{lb + 2})
      real(dp) :: tx(0:{la},0:{lb}), ty(0:{la},0:{lb}), tz(0:{la},0:{lb})
      real(dp) :: vn({nvn},0:{lab}), vn0({nvn}), vtmp({na*nb})
      real(dp) :: f(0:BOYS_MMAX)

      abx = ra(1) - rb(1); aby = ra(2) - rb(2); abz = ra(3) - rb(3)
      rab2 = abx*abx + aby*aby + abz*abz

      do k = 1, {na*nb}
         sout(k) = 0.0_dp; tout(k) = 0.0_dp; vout(k) = 0.0_dp
      end do

      do ip = 1, nprim_a
         al = ea(ip)
      do jp = 1, nprim_b
         bet = eb(jp)
         zeta = al + bet
         oo2z = 0.5_dp/zeta
         xi = al*bet/zeta
         cc = ca(ip)*cb(jp)
         tbb = 2.0_dp*bet*bet

         px = (al*ra(1) + bet*rb(1))/zeta
         py = (al*ra(2) + bet*rb(2))/zeta
         pz = (al*ra(3) + bet*rb(3))/zeta
         pax = px - ra(1); pay = py - ra(2); paz = pz - ra(3)
         pbx = px - rb(1); pby = py - rb(2); pbz = pz - rb(3)

         ! 1-D Gaussian-product prefactors; their product is
         ! (pi/zeta)^(3/2) exp(-xi |AB|^2), the 3-D overlap of two s functions.
         e0x = sqrt(PI/zeta)*exp(-xi*abx*abx)
         e0y = sqrt(PI/zeta)*exp(-xi*aby*aby)
         e0z = sqrt(PI/zeta)*exp(-xi*abz*abz)

{sxyz}

{txyz}

{comps}

         ! --- nuclear attraction: one VRR per nucleus, accumulated ---
         vpref = 2.0_dp*PI/zeta*exp(-xi*rab2)
         do ic = 1, natm
            pcx = px - ratm(1, ic)
            pcy = py - ratm(2, ic)
            pcz = pz - ratm(3, ic)
            tval = zeta*(pcx*pcx + pcy*pcy + pcz*pcz)
            call boys_eval({lab}, tval, f)

{vrr}

            do k = 1, {nvn}
               vn0(k) = vn(k,0)
            end do

{hrr}

            do k = 1, {na*nb}
               vout(k) = vout(k) - cc*zatm(ic)*vtmp(k)
            end do
         end do
      end do
      end do
   end subroutine one_e_{tag}
"""


def emit_multipole_module(L, path):
    """Multipole kernels, one routine per class, in their own module."""
    parts = [f"""!
! Cartesian multipole moments through the octupole.
!
! GENERATED by scripts/gen_1e.py --lmax {L} --multipole -- do not edit.
!
! Separate from the S/T/V kernels ON PURPOSE. Fusing them shares the 1-D
! overlap tables, which sounds free and is not: 39 components per class in the
! same `!$acc routine seq` body took NVVM over thirty minutes at LMAX=2, for
! two quantities that are each computed ONCE per geometry. Rebuilding the
! tables here costs microseconds and the compile finishes.
!
! Full Cartesian tensors -- 3, 9 and 27 components, repeats included -- in
! libcint's order with the last index fastest, because that is what mqc's
! `multipole_matrices` returns and a component map is one more thing to drift.
!
module trc_mult_kernels
   use trc_boys, only: dp
   implicit none
   private
   public :: multipole_dispatch, mult_driver, TRC_NMULT, MULT_LMAX

   integer, parameter :: MULT_LMAX = {L}
   !! Largest Cartesian block a shell pair produces here; the driver's
   !! per-item buffer is this size because `do concurrent` locals are fixed.
   integer, parameter :: NC_MAX = (MULT_LMAX + 1)*(MULT_LMAX + 2)/2
   integer, parameter :: NBLK_MAX = NC_MAX*NC_MAX
   !! Dispatch radix, FIXED and independent of MULT_LMAX -- see gen_1e.py.
   integer, parameter :: ONE_E_RADIX = __RADIX__
   !! 3 dipole + 9 quadrupole + 27 octupole
   integer, parameter :: TRC_NMULT = 39

contains

#include "inc/trc_mult_driver.inc"
"""]
    names = []
    for la in range(L + 1):
        for lb in range(L + 1):
            body, nmult = emit_multipole(la, lb)
            na, nb = len(cart(la)), len(cart(lb))
            sxyz = "\n".join(emit_1d(la, lb, a, extra=MULT_MAX) for a in "xyz")
            names.append((la, lb))
            parts.append(f"""
   pure subroutine mult_{la}{lb}(nprim_a, nprim_b, ea, ca, eb, cb, ra, rb, &
                                 orig, mout)
      !$acc routine seq
      integer,  intent(in)  :: nprim_a, nprim_b
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a)
      real(dp), intent(in)  :: eb(nprim_b), cb(nprim_b)
      real(dp), intent(in)  :: ra(3), rb(3), orig(3)
      real(dp), intent(out) :: mout({na*nb},{nmult})

      integer  :: ip, jp, k, ic
      real(dp) :: al, bet, zeta, oo2z, xi, cc
      real(dp) :: px, py, pz, pax, pay, paz, pbx, pby, pbz
      real(dp) :: abx, aby, abz, e0x, e0y, e0z, ox, oy, oz
      real(dp) :: sx(0:{la + MULT_MAX},0:{lb + 2})
      real(dp) :: sy(0:{la + MULT_MAX},0:{lb + 2})
      real(dp) :: sz(0:{la + MULT_MAX},0:{lb + 2})
      real(dp) :: mx1(0:{la},0:{lb}), my1(0:{la},0:{lb}), mz1(0:{la},0:{lb})
      real(dp) :: mx2(0:{la},0:{lb}), my2(0:{la},0:{lb}), mz2(0:{la},0:{lb})
      real(dp) :: mx3(0:{la},0:{lb}), my3(0:{la},0:{lb}), mz3(0:{la},0:{lb})

      abx = ra(1) - rb(1); aby = ra(2) - rb(2); abz = ra(3) - rb(3)
      ox = ra(1) - orig(1); oy = ra(2) - orig(2); oz = ra(3) - orig(3)

      do ic = 1, {nmult}
         do k = 1, {na*nb}
            mout(k, ic) = 0.0_dp
         end do
      end do

      do ip = 1, nprim_a
         al = ea(ip)
      do jp = 1, nprim_b
         bet = eb(jp)
         zeta = al + bet
         oo2z = 0.5_dp/zeta
         xi = al*bet/zeta
         cc = ca(ip)*cb(jp)
         px = (al*ra(1) + bet*rb(1))/zeta
         py = (al*ra(2) + bet*rb(2))/zeta
         pz = (al*ra(3) + bet*rb(3))/zeta
         pax = px - ra(1); pay = py - ra(2); paz = pz - ra(3)
         pbx = px - rb(1); pby = py - rb(2); pbz = pz - rb(3)
         e0x = sqrt(PI_/zeta)*exp(-xi*abx*abx)
         e0y = sqrt(PI_/zeta)*exp(-xi*aby*aby)
         e0z = sqrt(PI_/zeta)*exp(-xi*abz*abz)

{sxyz}

{body}
      end do
      end do
   end subroutine mult_{la}{lb}
""")
    parts.append(f"""
   pure subroutine multipole_dispatch(la, lb, nprim_a, nprim_b, ea, ca, eb, cb, &
                                      ra, rb, orig, mout)
      !$acc routine seq
      integer,  intent(in)  :: la, lb, nprim_a, nprim_b
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a)
      real(dp), intent(in)  :: eb(nprim_b), cb(nprim_b)
      real(dp), intent(in)  :: ra(3), rb(3), orig(3)
      real(dp), intent(out) :: mout(*)
      select case (la*ONE_E_RADIX + lb)""")
    for la, lb in names:
        parts.append(f"""      case ({la*ONE_E_RADIX + lb}); call mult_{la}{lb}(nprim_a, nprim_b, ea, ca, eb, cb, &
                          ra, rb, orig, mout)""")
    parts.append("""      case default
         !! No kernel for this class. A sentinel rather than silence: `mout`
         !! is intent(out) on an assumed-size dummy, so doing nothing returns
         !! whatever the caller's buffer held. huge() and not a NaN because
         !! this is device code under -fast, where a NaN is laundered into a
         !! bound by the first min/max it meets.
         mout(1) = huge(1.0_dp)
      end select
   end subroutine multipole_dispatch

end module trc_mult_kernels
""")
    txt = "\n".join(parts).replace("PI_", "3.14159265358979323846_dp")
    txt = txt.replace("__RADIX__", str(ONE_E_RADIX))
    with open(path, "w") as fh:
        fh.write(txt)
    print(f"wrote {path}: {len(names)} classes, {txt.count(chr(10))} lines")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=2)
    ap.add_argument("-o", "--output", default="src/trc_1e_kernels.F90")
    ap.add_argument("--multipole", metavar="FILE", default=None,
                    help="emit the multipole kernels into their own module "
                         "instead. They are NOT fused with S/T/V: doing that "
                         "put 39 components per class into the same device "
                         "routine and NVVM ground for over thirty minutes at "
                         "LMAX=2. Both run once per geometry, so sharing the "
                         "1-D tables saved nothing worth that.")
    args = ap.parse_args()
    L = args.lmax

    if args.multipole:
        emit_multipole_module(L, args.multipole)
        return

    parts = [f"""!
! One-electron kernels: overlap, kinetic, nuclear attraction.
!
! GENERATED by scripts/gen_1e.py --lmax {L} -- do not edit.
!
! One straight-line routine per (la, lb) class, in the same style the
! four-centre kernels established: literal indices, no runtime loop bounds over
! angular momentum, and the routine in the same module as its caller so
! nvfortran will inline it.
!
! Overlap and kinetic are separable and need no auxiliary index; nuclear
! attraction is not separable and carries the Boys function and an m index.
! Both use the SAME horizontal transfer as the ERI, because the HRR is a
! statement about the two Gaussian centres rather than the operator between
! them. See the module docstring of scripts/gen_1e.py.
!
module trc_1e_kernels
   use trc_boys, only: dp, BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV
   implicit none
   private
   public :: one_e_dispatch, one_e_driver, ONE_E_LMAX, TRC_NMULT

   !! Dipole + quadrupole + octupole, full Cartesian tensors.
   integer, parameter :: TRC_NMULT = 3 + 9 + 27

   integer, parameter :: ONE_E_LMAX = {L}
   !! Dispatch radix, FIXED and independent of ONE_E_LMAX -- see gen_1e.py.
   integer, parameter :: ONE_E_RADIX = __RADIX__
   real(dp), parameter :: PI = 3.14159265358979323846_dp
   !! Largest Cartesian block a shell pair produces here; the driver's
   !! per-item buffers are this size because `do concurrent` locals are fixed.
   integer, parameter :: NC_MAX = (ONE_E_LMAX + 1)*(ONE_E_LMAX + 2)/2
   integer, parameter :: NBLK_MAX = NC_MAX*NC_MAX

contains

#include "inc/trc_boys_eval.inc"
#include "inc/trc_1e_driver.inc"
"""]

    names = []
    for la in range(L + 1):
        for lb in range(L + 1):
            parts.append(emit_class(la, lb, L))
            names.append((la, lb))

    parts.append(f"""
   !> Dispatch on the class. Host-side call, so the module boundary the
   !> kernels care about is not crossed here.
   pure subroutine one_e_dispatch(la, lb, nprim_a, nprim_b, ea, ca, eb, cb, &
                                  ra, rb, natm, zatm, ratm, sout, tout, vout)
      !$acc routine seq
      integer,  intent(in)  :: la, lb, nprim_a, nprim_b, natm
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a)
      real(dp), intent(in)  :: eb(nprim_b), cb(nprim_b)
      real(dp), intent(in)  :: ra(3), rb(3)
      real(dp), intent(in)  :: zatm(natm), ratm(3, natm)
      real(dp), intent(out) :: sout(*), tout(*), vout(*)
      select case (la*ONE_E_RADIX + lb)""")
    for la, lb in names:
        parts.append(
            f"""      case ({la*ONE_E_RADIX + lb}); call one_e_{la}{lb}(nprim_a, nprim_b, ea, ca, eb, cb, &
                          ra, rb, natm, zatm, ratm, sout, tout, vout)""")
    parts.append("""      case default
         !! No kernel for this class -- ONE_E_LMAX is below the basis. A
         !! sentinel rather than silence: these are intent(out) on assumed-size
         !! dummies, so doing nothing hands back the caller's stale buffer as
         !! an overlap matrix, and an SCF converges happily on it. huge() and
         !! not a NaN because this is device code under -fast, where the first
         !! min/max launders a NaN into a bound.
         !!
         !! trc_basis refuses such a basis on the host before any of this runs;
         !! this arm is the backstop for a path that skipped that check.
         sout(1) = huge(1.0_dp)
         tout(1) = huge(1.0_dp)
         vout(1) = huge(1.0_dp)
      end select
   end subroutine one_e_dispatch

end module trc_1e_kernels
""")

    txt = "\n".join(parts).replace("__RADIX__", str(ONE_E_RADIX))
    with open(args.output, "w") as fh:
        fh.write(txt)
    print(f"wrote {args.output}: {len(names)} classes, {txt.count(chr(10))} lines")


if __name__ == "__main__":
    main()
