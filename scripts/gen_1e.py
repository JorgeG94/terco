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
import os

import cartorder

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


def cart(l, order="lexicographic"):
    """(nx,ny,nz) for one shell, in `order`.

    Thin wrapper: the order itself is the single definition in cartorder.py.
    Called with NO explicit order (the default) at every INTERNAL use below
    -- `cart_cum`'s index into the VRR's `vn`/`g` arrays, and
    `emit_nuc_vrr`'s own level-by-level bookkeeping -- so those stay pinned
    to lexicographic regardless of what `--order` asks for. Only
    `emit_hrr`'s and `emit_class`'s/`emit_multipole`'s OUTPUT enumerations
    (the AO slot order) take `order` explicitly.
    """
    return cartorder.components(l, order)


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


def emit_hrr(la, lb, idx, src, dst, order="lexicographic"):
    """(a|b) from [a+b|0] -- the same horizontal transfer the ERI uses.

    It is a statement about the two Gaussian centres, not about the operator
    between them, so it is reused verbatim.

    `order` governs only the two outer loops (`dst`'s own AO slot); `idx`
    (built by `cart_cum` with no order argument by every caller) addresses
    `src`, an internal VRR array, and stays lexicographic regardless.
    """
    o, w = [], None
    w = o.append
    na = len(cart(la))
    for ib, (bx, by, bz) in enumerate(cart(lb, order)):
        for ia, (axp, ayp, azp) in enumerate(cart(la, order)):
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


def emit_multipole(la, lb, order="lexicographic"):
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
    6 and 10 unique ones -- with the last index fastest. That MOMENT-RANK
    axis (dipole/quadrupole/octupole component slot) is a fixed tensor
    convention unrelated to `--order`: it is not a within-shell Cartesian
    layout, so it is untouched here regardless of which order the AO axis
    (the `k` index below, this class's `(la|lb)` slot) follows. The repeats
    are written out rather than mapped, because a component table is another
    thing to get out of step with the caller.
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

    # the moment-rank components, last index fastest (see the docstring: this
    # is not the Cartesian shell order and does not take `order`)
    slot = 0
    for mrank in range(1, MULT_MAX + 1):
        for idx in range(3**mrank):
            # decode the index tuple, last fastest
            t, powers = idx, [0, 0, 0]
            for _ in range(mrank):
                powers[t % 3] += 1
                t //= 3
            slot += 1
            for ib, (bx, by, bz) in enumerate(cart(lb, order)):
                for ia, (ax_, ay_, az_) in enumerate(cart(la, order)):
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


def emit_class(la, lb, LMAX, order="lexicographic"):
    lab = la + lb
    na, nb = len(cart(la)), len(cart(lb))
    idx = cart_cum(lab)
    nvn = ncum(lab)
    tag = f"{la}{lb}"

    sxyz = "\n".join(emit_1d(la, lb, a) for a in "xyz")
    txyz = "\n".join(emit_kin_1d(la, lb, a) for a in "xyz")

    # overlap and kinetic components, straight from the 1-D tables. OUTPUT:
    # sout/tout are AO-indexed, so this enumeration follows `order`.
    comps = []
    for ib, (bx, by, bz) in enumerate(cart(lb, order)):
        for ia, (ax, ay, az) in enumerate(cart(la, order)):
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
    hrr = emit_hrr(la, lb, idx, "vn0", "vtmp", order)

    return f"""
   pure subroutine one_e_{tag}(nprim_a, nprim_b, ea, ca, eb, cb, ra, rb, &
                               natm, zatm, ratm, sout, tout, vout)
      !! ({la}|{lb}) one-electron block: overlap, kinetic, nuclear attraction.
      !!
      !! One routine for all three because they share the shell-pair quantities
      !! (zeta, P, the Gaussian product prefactor) and mqc needs all three to form
      !! the core Hamiltonian. Computing them together reads the primitive data
      !! once instead of three times.
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


def multipole_class_routine(la, lb, order="lexicographic"):
    """The pure mult_<la><lb> routine: one class of multipoles, own 1-D tables."""
    body, nmult = emit_multipole(la, lb, order)
    na, nb = len(cart(la)), len(cart(lb))
    sxyz = "\n".join(emit_1d(la, lb, a, extra=MULT_MAX) for a in "xyz")
    return f"""
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
"""


def emit_multipole_module(L, path, order="lexicographic"):
    """Multipole kernels, one routine per class, in their own module."""
    order_val = 0 if order == "lexicographic" else 1
    order_arg = "" if order == "lexicographic" else f" --order {order}"
    parts = [f"""!
! Cartesian multipole moments through the octupole.
!
! GENERATED by scripts/gen_1e.py --lmax {L} --multipole{order_arg} -- do not edit.
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
! TRC_CART_ORDER GUARD. An undefined TRC_CART_ORDER reads as 0 in a
! preprocessor #if, so the default (lexicographic) build needs no new -D
! flag; see docs/GENERATORS.md's TRC_LMAX guard for the failure mode this
! prevents.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
module trc_mult_kernels
   use pic_types, only: dp
   implicit none
   private
   public :: multipole_dispatch, mult_driver, TRC_NMULT, MULT_LMAX

   integer, parameter :: MULT_LMAX = {L}
   integer, parameter :: NC_MAX = (MULT_LMAX + 1)*(MULT_LMAX + 2)/2
   !! Largest Cartesian block a shell pair produces here; the driver's
   !! per-item buffer is this size because `do concurrent` locals are fixed.
   integer, parameter :: NBLK_MAX = NC_MAX*NC_MAX
   integer, parameter :: ONE_E_RADIX = __RADIX__
   !! Dispatch radix, FIXED and independent of MULT_LMAX -- see gen_1e.py.
   integer, parameter :: TRC_NMULT = 39
   !! 3 dipole + 9 quadrupole + 27 octupole

contains

#include "inc/trc_mult_driver.inc"
"""]
    names = []
    for la in range(L + 1):
        for lb in range(L + 1):
            names.append((la, lb))
            parts.append(multipole_class_routine(la, lb, order))
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


def emit_pcm_class(la, lb, order="lexicographic"):
    """One (la|lb) block of the smeared-charge potential, summed over points.

    THE OPERATOR. A PCM surface point k is a normalised s Gaussian of
    exponent zeta_k^2, whose potential is erf(zeta_k r)/r. So the integral is
    the nuclear attraction one with 1/|r - C| attenuated, and the same
    substitution the four-centre kernels make for range separation applies,
    with the pair exponent p in the role rho plays there and omega = zeta_k:

        erf(w r)/r = (2/sqrt(pi)) int_0^w exp(-t^2 r^2) dt

    Through the Gaussian product theorem, the t integrand for a primitive
    pair is an overlap-like integral centred on Q = (pP + t^2 C)/(p + t^2).
    With u^2 = t^2/(p + t^2) its Obara-Saika recurrence has QA = PA - u^2 PC
    and 1/(2(p + t^2)) = (1 - u^2)/(2p), which is exactly the point-charge
    recurrence above (emit_nuc_vrr) once each power u^(2m) is folded into an
    auxiliary index. Only the range of u changes: t in [0, inf) is u in
    [0, 1) for 1/r, while t in [0, w] is u in [0, sqrt(theta)],

        theta = w^2/(w^2 + p).

    So the auxiliary integrals are

        V^(m) = (2 pi/p) K_ab int_0^sqrt(theta) u^(2m) exp(-T u^2) du
              = (2 pi/p) K_ab sqrt(theta) theta^m F_m(theta T),   T = p |PC|^2

    and the recurrence that consumes them is emit_nuc_vrr's, UNCHANGED: the
    same PA, PC and 1/(2p). That is why this kernel scales F_m by theta^m and
    the prefactor by sqrt(theta) and nothing else. The four-centre kernels
    write the same result the other way round, W^(m) = V^(m)/theta^m, which
    puts theta on every coefficient that multiplies an (m+1) term instead
    (theta*W, theta*rho/zeta, ...); the two are equal, and this one leaves
    the shared recurrence emitter alone. theta -> 1 (a point charge) is the
    nuclear attraction exactly.

    THE SUM. HRR depends on A - B alone, not on the operator or on C, so the
    b = 0 VRR results are summed over primitive pairs and over points (with
    their weights) first, and transferred once.
    """
    lab = la + lb
    na, nb = len(cart(la)), len(cart(lb))
    idx = cart_cum(lab)
    nvn = ncum(lab)
    tag = f"{la}{lb}"
    vrr = emit_nuc_vrr(lab, idx)
    hrr = emit_hrr(la, lb, idx, "vacc", "blk", order)
    scale = ""
    im = ", m" if lab > 0 else ""
    rthm = ", thm" if lab > 0 else ""
    if lab > 0:
        scale = f"""            thm = 1.0_dp
            do m = 1, {lab}
               thm = thm*theta
               f(m) = f(m)*thm
            end do
"""
    return f"""
   pure subroutine pcm_{tag}(npp, ppz, ppk, ppp, ppa, abx, aby, abz, npt, k0, kstep, cpt, zk, wq, blk)
      !! ({la}|{lb}) block of sum_k wq_k erf(zeta_k|r - C_k|)/|r - C_k| over the
      !! npt points k = k0, k0 + kstep, ...
      !!
      !! Strided rather than a contiguous run so that the drivers can deal a
      !! pair's points out round-robin: neighbouring threads then read
      !! neighbouring points, coalesced, where contiguous chunks put each
      !! thread's reads a chunk apart (measured 6x slower on the Fock term).
      !$acc routine seq
      integer,  intent(in)  :: npp, npt, k0, kstep
      real(dp), intent(in)  :: ppz(npp), ppk(npp), ppp(3, npp), ppa(3, npp)
      real(dp), intent(in)  :: abx, aby, abz
      real(dp), intent(in)  :: cpt(3, *)     !! point centres
      real(dp), intent(in)  :: zk(*)         !! zeta_k; the Gaussian exponent is zeta_k^2
      real(dp), intent(in)  :: wq(*)         !! weight per point (a charge, or 1)
      real(dp), intent(out) :: blk({na*nb})

      integer  :: ip, ic, kk, k{im}
      real(dp) :: zeta, oo2z, pax, pay, paz, pcx, pcy, pcz
      real(dp) :: w2, theta, tval, vpref{rthm}
      real(dp) :: vn({nvn},0:{lab}), vacc({nvn})
      real(dp) :: f(0:BOYS_MMAX)

      do k = 1, {nvn}
         vacc(k) = 0.0_dp
      end do

      do ip = 1, npp
         zeta = ppz(ip)
         oo2z = 0.5_dp/zeta
         pax = ppa(1, ip); pay = ppa(2, ip); paz = ppa(3, ip)
         do ic = 1, npt
            kk = k0 + (ic - 1)*kstep
            pcx = ppp(1, ip) - cpt(1, kk)
            pcy = ppp(2, ip) - cpt(2, kk)
            pcz = ppp(3, ip) - cpt(3, kk)
            w2 = zk(kk)*zk(kk)
            theta = w2/(w2 + zeta)
            tval = theta*zeta*(pcx*pcx + pcy*pcy + pcz*pcz)
            call boys_eval({lab}, tval, f)
{scale}            vpref = wq(kk)*ppk(ip)*(2.0_dp*PI/zeta)*sqrt(theta)

{vrr}

            do k = 1, {nvn}
               vacc(k) = vacc(k) + vn(k,0)
            end do
         end do
      end do

{hrr}
   end subroutine pcm_{tag}
"""


def emit_pcm_module(L, path, order="lexicographic"):
    """The PCM surface-integral kernels, one routine per class, own module."""
    order_val = 0 if order == "lexicographic" else 1
    order_arg = "" if order == "lexicographic" else f" --order {order}"
    parts = [f"""!
! Smeared-charge surface integrals for the polarisable continuum.
!
! GENERATED by scripts/gen_1e.py --lmax {L} --pcm{order_arg} -- do not edit.
!
! <a| erf(zeta_k |r - C_k|)/|r - C_k| |b> for the Gaussian surface charges
! of trc_pcm: the nuclear-attraction recurrence with the Boys argument
! scaled by theta = zeta^2/(zeta^2 + p), F_m by theta^m and the prefactor by
! sqrt(theta). See emit_pcm_class in scripts/gen_1e.py for the derivation.
! One routine per (la, lb) class, summed over a list of points with a
! weight each, so a single call is either one point (the potential the
! density puts on it) or every point (the operator the charges add to the
! Fock matrix). The drivers are in inc/trc_pcm_driver.inc.
!
! A module of its own rather than more routines in trc_1e_kernels: the 1e
! kernels run once per geometry and the PCM ones twice per SCF iteration,
! and keeping them apart keeps each module's device compile small.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
module trc_pcm_kernels
   use pic_types, only: dp
   use trc_boys, only: BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV
   implicit none
   private
   public :: pcm_dispatch, pcm_potential_driver, pcm_fock_driver, PCM_LMAX

   integer, parameter :: PCM_LMAX = {L}
   integer, parameter :: ONE_E_RADIX = __RADIX__
   !! Dispatch radix, FIXED and independent of PCM_LMAX -- see gen_1e.py.
   real(dp), parameter :: PI = 3.14159265358979323846_dp
   integer, parameter :: NC_MAX = (PCM_LMAX + 1)*(PCM_LMAX + 2)/2
   integer, parameter :: NBLK_MAX = NC_MAX*NC_MAX
   !! Largest Cartesian block a shell pair produces; the drivers' per-item
   !! buffer is this size because `do concurrent` locals are fixed.
   integer, parameter :: PCM_PT_CHUNK = 128
   !! Points per Fock-term work item, roughly (inc/trc_pcm_driver.inc)
   integer, parameter :: PCM_PAIR_CHUNK = 256
   !! Shell pairs per potential work item

contains

#include "inc/trc_boys_eval.inc"
#include "inc/trc_pcm_driver.inc"
"""]
    names = []
    for la in range(L + 1):
        for lb in range(L + 1):
            parts.append(emit_pcm_class(la, lb, order))
            names.append((la, lb))
    parts.append("""
   pure subroutine pcm_dispatch(la, lb, npp, ppz, ppk, ppp, ppa, abx, aby, abz, &
                                npt, k0, kstep, cpt, zk, wq, blk)
      !$acc routine seq
      integer,  intent(in)  :: la, lb, npp, npt, k0, kstep
      real(dp), intent(in)  :: ppz(npp), ppk(npp), ppp(3, npp), ppa(3, npp)
      real(dp), intent(in)  :: abx, aby, abz
      real(dp), intent(in)  :: cpt(3, *), zk(*), wq(*)
      real(dp), intent(out) :: blk(*)
      select case (la*ONE_E_RADIX + lb)""")
    for la, lb in names:
        parts.append(f"""      case ({la*ONE_E_RADIX + lb}); call pcm_{la}{lb}(npp, ppz, ppk, ppp, ppa, abx, aby, abz, &
                          npt, k0, kstep, cpt, zk, wq, blk)""")
    parts.append("""      case default
         !! No kernel for this class: the basis is above PCM_LMAX. A
         !! sentinel rather than silence, as in one_e_dispatch; trc_pcm
         !! refuses such a basis on the host before any of this runs.
         blk(1) = huge(1.0_dp)
      end select
   end subroutine pcm_dispatch

end module trc_pcm_kernels
""")
    txt = "\n".join(parts).replace("__RADIX__", str(ONE_E_RADIX))
    with open(path, "w") as fh:
        fh.write(txt)
    print(f"wrote {path}: {len(names)} classes, {txt.count(chr(10))} lines")


def emit_dc_vrr(lab, idx):
    """Hellmann-Feynman derivative wrt the nucleus C: wx/wy/wz(:, 0:lab),
    built from the ALREADY-COMPUTED vn(:, 0:lab+1) -- one Boys order higher
    than this function's own `lab`, supplied by the caller's VRR call at
    lab+1 (see emit_grad_class).

    DERIVATION. C enters V[a]^(m) only through PC_i = P_i - C_i: once
    explicitly, in the lowering recurrence's own coefficient, and once more
    through every V[...]^(m+1) nested inside it. Differentiate

        V[a+1_i]^(m) = PA_i V[a]^(m) - PC_i V[a]^(m+1)
                       + (N_i(a)/2zeta)(V[a-1_i]^(m) - V[a-1_i]^(m+1))

    term by term wrt C_x. PA_i does not depend on C, so the first term's
    derivative is PA_i W_x[a]^(m) outright. PC_i does (d PC_i/dC_x =
    -delta_{ix}), so its term differentiates by the product rule into TWO
    pieces: -PC_i W_x[a]^(m+1) from V[a]^(m+1) itself varying, and
    +delta_{ix} V[a]^(m+1) from PC_i varying. The lowering term carries no
    explicit PC_i, so it differentiates straight through to its own W. So:

        W_x[a+1_i]^(m) = PA_i W_x[a]^(m) - PC_i W_x[a]^(m+1)
                         + delta_{ix} V[a]^(m+1)
                         + (N_i(a)/2zeta)(W_x[a-1_i]^(m) - W_x[a-1_i]^(m+1))

    the SAME recurrence as V, plus one inhomogeneous "source" term that
    fires only when the direction being raised (i) is the direction being
    differentiated (x): raising momentum in a direction and differentiating
    wrt that same direction both want the one explicit factor of PC_i, and
    whichever one is applied first leaves the other as a plain V, not a W.

    The base case (a = 0, no predecessor to inherit a recurrence from) is
    derived directly from the Boys substitution: V[0]^(m) = (2pi/zeta) K_ab
    F_m(T), T = zeta|PC|^2, dF_m/dT = -F_{m+1}(T), d T/dC_x = -2 zeta PC_x,
    so

        W_x[0]^(m) = (2pi/zeta) K_ab * (-F_{m+1}) * (-2 zeta PC_x)
                    = 2 zeta PC_x * (2pi/zeta) K_ab F_{m+1}
                    = 2 zeta PC_x V[0]^(m+1)

    which is exactly the general formula's pattern with the inhomogeneous
    term supplying the whole thing (there is no recurrence to add to at
    a = 0).

    Checked three ways in test/check_1e_grad_vs_pyscf.F90 and
    test/check_1e_grad_fd.F90: against pyscf's own Hellmann-Feynman term
    (grad.rhf's hcore_generator), against a central finite difference of
    the nuclear-attraction energy wrt a nuclear coordinate at fixed
    density and fixed basis-function positions, and via translational
    invariance (dV/dA + dV/dB + dV/dC = 0 for the one (A,B,C) triple) --
    never trusted on the strength of this derivation alone.
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
    WN = ["wx", "wy", "wz"]
    idx0 = idx[(0, 0, 0)]
    for m in range(lab, -1, -1):
        w(f"            ! --- dV/dC, m = {m} ---")
        for dd in range(3):
            w(f"            {WN[dd]}({idx0},{m}) = "
              f"2.0_dp*zeta*{PC[dd]}*vn({idx0},{m + 1})")
        for l in range(1, lab - m + 1):
            for p in cart(l):
                d, dn1, dn2, n = info[p]
                for dd in range(3):
                    t = (f"{PA[d]}*{WN[dd]}({idx[dn1]},{m})"
                         f" - {PC[d]}*{WN[dd]}({idx[dn1]},{m + 1})")
                    if dd == d:
                        t += f" &\n               + vn({idx[dn1]},{m + 1})"
                    if dn2 is not None:
                        t += (f" &\n               + {float(n - 1)}_dp*oo2z*"
                              f"({WN[dd]}({idx[dn2]},{m}) - {WN[dd]}({idx[dn2]},{m + 1}))")
                    w(f"            {WN[dd]}({idx[p]},{m}) = {t}")
    return "\n".join(o)


def emit_grad_class(la, lb, order="lexicographic"):
    """(la|lb) one-electron GRADIENT block: -W:dS/dR, +D:dT/dR and
    +D:dV/dR, contracted directly into ga (bra atom), gb (ket atom) and gc
    (every atom -- V's Hellmann-Feynman term), with no derivative matrix
    ever materialised.

    BRA/KET BASIS DERIVATIVE. A Cartesian Gaussian's derivative wrt its own
    centre shifts its angular momentum by one unit each way,

        d/dA_x [(x-Ax)^l exp(-a(x-Ax)^2)] = 2a (x-Ax)^(l+1) exp(...)
                                            - l (x-Ax)^(l-1) exp(...)

    independent of the ket or of whatever operator sits between bra and
    ket. S and T have no third centre, so dO/dA applied via that identity
    plus translational invariance (dO/dA + dO/dB = 0, exact -- O depends on
    nothing C could move) gives dO/dB = -dO/dA for free; that is why only
    one bra-raised set of 1-D tables (sx/sy/sz to la+1, tx/ty/tz likewise)
    is built below, and `gb`'s S/T lines are the negative of `ga`'s.

    V also depends on the nucleus C, so dV/dA + dV/dB = -dV/dC is no
    longer zero and the shortcut does not apply to B: dV/dB is built the
    same identity applied to the KET centre directly, which needs the
    (la,lb+1) and (la,lb-1) classes. Both of those and the bra's
    (la+1,lb)/(la-1,lb) come from ONE extra unit of VRR -- build vn to
    degree la+lb+1 (one more than the plain V kernel ever needs) and HRR-
    transfer the SAME vn0 at four (la',lb') pairs instead of one.

    OPERATOR DERIVATIVE (Hellmann-Feynman). emit_dc_vrr's own docstring has
    the derivation; it rides on the same vn(:, 0:la+lb+1) already built for
    the bra/ket shift above, so no further VRR is needed, only the HRR
    transfer of wx/wy/wz at the plain (la,lb) class.

    Translational invariance is used as a CHECK on the assembled gradient
    (test/check_1e_grad_fd.F90, test/check_1e_grad_vs_pyscf.F90), not as a
    shortcut here for anything beyond the dO/dB = -dO/dA already described.

    WEIGHT. `wt` folds the triangular shell-pair list's double count
    (ii != jj feeds both S(mu,nu) and S(nu,mu), equal numbers, both
    needed in the sum) into one scale on `cc`; every term below is linear
    in `cc`, so nothing else has to know about it.
    """
    lab = la + lb
    lab1 = lab + 1
    na, nb = len(cart(la)), len(cart(lb))
    na1, nb1 = len(cart(la + 1)), len(cart(lb + 1))
    nam1 = len(cart(la - 1)) if la >= 1 else 0
    nbm1 = len(cart(lb - 1)) if lb >= 1 else 0
    idx = cart_cum(lab1)
    nvn = ncum(lab1)
    nvnw = ncum(lab)
    tag = f"{la}{lb}"

    sxyz = "\n".join(emit_1d(la + 1, lb, a) for a in "xyz")
    txyz = "\n".join(emit_kin_1d(la + 1, lb, a) for a in "xyz")

    cla, clb = cart(la, order), cart(lb, order)
    cla1, clb1 = cart(la + 1, order), cart(lb + 1, order)
    clam1 = cart(la - 1, order) if la >= 1 else None
    clbm1 = cart(lb - 1, order) if lb >= 1 else None
    AX = {"x": 0, "y": 1, "z": 2}

    # --- S/T bra derivative; gb is -ga here (translational invariance) ---
    st_lines = []
    wst = st_lines.append
    for ib, bcomp in enumerate(clb):
        for ia, acomp in enumerate(cla):
            k = 1 + ia + na*ib
            comp = {"x": acomp[0], "y": acomp[1], "z": acomp[2]}
            bco = {"x": bcomp[0], "y": bcomp[1], "z": bcomp[2]}
            for d in "xyz":
                cd = comp[d]
                others = [o for o in "xyz" if o != d]
                o1, o2 = others

                def raised(tab, d=d, cd=cd, bco=bco):
                    hi = f"{tab}{d}({cd + 1},{bco[d]})"
                    if cd > 0:
                        return (f"(2.0_dp*al*{hi} - {float(cd)}_dp*"
                                f"{tab}{d}({cd - 1},{bco[d]}))")
                    return f"(2.0_dp*al*{hi})"

                ds_d = raised("s")
                dt_d = raised("t")
                sprod = "*".join(f"s{o}({comp[o]},{bco[o]})" for o in others)
                dS_total = f"{ds_d}*{sprod}"
                term2 = f"t{o1}({comp[o1]},{bco[o1]})*s{o2}({comp[o2]},{bco[o2]})"
                term3 = f"s{o1}({comp[o1]},{bco[o1]})*t{o2}({comp[o2]},{bco[o2]})"
                dT_total = f"{dt_d}*{sprod} + {ds_d}*({term2} + {term3})"
                di = AX[d] + 1
                wst(f"         ga({di}) = ga({di}) - cc*wblk({k})*({dS_total}) &")
                wst(f"                         + cc*dblk({k})*({dT_total})")
                wst(f"         gb({di}) = gb({di}) + cc*wblk({k})*({dS_total}) &")
                wst(f"                         - cc*dblk({k})*({dT_total})")
    st_block = "\n".join(st_lines)

    # --- V: VRR to lab+1, Hellmann-Feynman VRR, five HRR transfers ---
    vrr = emit_nuc_vrr(lab1, idx)
    dcvrr = emit_dc_vrr(lab, idx)
    hrr_ap = emit_hrr(la + 1, lb, idx, "vn0", "vAp", order)
    hrr_am = emit_hrr(la - 1, lb, idx, "vn0", "vAm", order) if la >= 1 else \
        "            ! la = 0: no (la-1,lb) class"
    hrr_bp = emit_hrr(la, lb + 1, idx, "vn0", "vBp", order)
    hrr_bm = emit_hrr(la, lb - 1, idx, "vn0", "vBm", order) if lb >= 1 else \
        "            ! lb = 0: no (la,lb-1) class"
    hrr_wx = emit_hrr(la, lb, idx, "wx0", "wxt", order)
    hrr_wy = emit_hrr(la, lb, idx, "wy0", "wyt", order)
    hrr_wz = emit_hrr(la, lb, idx, "wz0", "wzt", order)

    v_lines = []
    wv = v_lines.append
    for ib, bcomp in enumerate(clb):
        for ia, acomp in enumerate(cla):
            k = 1 + ia + na*ib
            for d in "xyz":
                di = AX[d] + 1
                j = AX[d]
                ap = tuple(c + (1 if jj == j else 0) for jj, c in enumerate(acomp))
                k1 = 1 + cla1.index(ap) + na1*ib
                bra_expr = f"2.0_dp*al*vAp({k1})"
                if acomp[j] > 0:
                    am = tuple(c - (1 if jj == j else 0) for jj, c in enumerate(acomp))
                    k2 = 1 + clam1.index(am) + nam1*ib
                    bra_expr = f"({bra_expr} - {float(acomp[j])}_dp*vAm({k2}))"
                else:
                    bra_expr = f"({bra_expr})"
                bp = tuple(c + (1 if jj == j else 0) for jj, c in enumerate(bcomp))
                k3 = 1 + ia + na*clb1.index(bp)
                ket_expr = f"2.0_dp*bet*vBp({k3})"
                if bcomp[j] > 0:
                    bm = tuple(c - (1 if jj == j else 0) for jj, c in enumerate(bcomp))
                    k4 = 1 + ia + na*clbm1.index(bm)
                    ket_expr = f"({ket_expr} - {float(bcomp[j])}_dp*vBm({k4}))"
                else:
                    ket_expr = f"({ket_expr})"
                wtab = {"x": "wxt", "y": "wyt", "z": "wzt"}[d]
                wv(f"            ga({di}) = ga({di}) - cc*zatm(ic)*dblk({k})*{bra_expr}")
                wv(f"            gb({di}) = gb({di}) - cc*zatm(ic)*dblk({k})*{ket_expr}")
                # `gc` is the ONLY one of the three that is natm-wide: ga/gb
                # each touch one atom (the bra/ket shell's own), same as
                # the validated eri%gradient ga/gb/gc/gd pattern, but gc
                # touches EVERY atom from EVERY shell pair -- the width
                # that turns "lane sharing happens to not collide" into
                # "collides often enough to be visibly wrong". Measured
                # non-deterministic, non-atomic: two runs of the same
                # input differed in their (wrong) answer, a clean race
                # signature neither compute-sanitizer memcheck nor
                # racecheck caught (gc lives in global memory; racecheck
                # only instruments __shared__). An atomic update is the fix
                # that does not depend on how many lanes happen to collide.
                wv("            !$acc atomic update")
                wv(f"            gc({di},ic) = gc({di},ic) - cc*zatm(ic)*dblk({k})*{wtab}({k})")
    v_block = "\n".join(v_lines)

    na1nb = na1*nb
    nanb1 = na*nb1
    decl_am = f"      real(dp) :: vAm({nam1*nb})\n" if la >= 1 else ""
    decl_bm = f"      real(dp) :: vBm({na*nbm1})\n" if lb >= 1 else ""

    return f"""
   pure subroutine one_e_grad_{tag}(nprim_a, nprim_b, ea, ca, eb, cb, ra, rb, &
                                    natm, zatm, ratm, wt, wblk, dblk, ga, gb, gc)
      !! ({la}|{lb}) one-electron gradient block: -W:dS/dR, +D:dT/dR, and
      !! +D:dV/dR (basis derivative on A and B, Hellmann-Feynman on every
      !! nucleus), contracted straight into ga (bra atom), gb (ket atom)
      !! and gc (every atom) -- no derivative matrix is ever formed. See
      !! emit_grad_class/emit_dc_vrr in scripts/gen_1e.py.
      !!
      !! gc is intent(inout) and ADDED into, never zeroed here, on
      !! purpose: the caller passes a slice of its own lane-private
      !! accumulator directly (trc_1e_grad_driver.inc's one_e_grad_item),
      !! not a per-thread local -- an automatic array sized by the RUNTIME
      !! value `natm` inside a `!$acc routine seq` body needs a per-thread
      !! dynamic stack allocation that a large `natm` (a 60-atom cluster,
      !! say) under thousands of concurrently resident threads overruns,
      !! which is exactly a CUDA_ERROR_LAUNCH_FAILED "invalid pointer
      !! dereference" and not a Fortran-level bug. Measured: that is what
      !! one_e_grad_item's own gc(3,natm) local used to be, and it is why
      !! this kernel does not get one either.
      !$acc routine seq
      integer,  intent(in)  :: nprim_a, nprim_b, natm
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a)
      real(dp), intent(in)  :: eb(nprim_b), cb(nprim_b)
      real(dp), intent(in)  :: ra(3), rb(3)
      real(dp), intent(in)  :: zatm(natm), ratm(3, natm)
      real(dp), intent(in)  :: wt, wblk({na*nb}), dblk({na*nb})
      real(dp), intent(out)   :: ga(3), gb(3)
      real(dp), intent(inout) :: gc(3, natm)

      integer  :: ip, jp, ic, k
      real(dp) :: al, bet, zeta, oo2z, xi, cc, tbb
      real(dp) :: px, py, pz, pax, pay, paz, pbx, pby, pbz
      real(dp) :: pcx, pcy, pcz, tval, vpref, rab2
      real(dp) :: abx, aby, abz, e0x, e0y, e0z
      real(dp) :: sx(0:{la + 1},0:{lb + 2}), sy(0:{la + 1},0:{lb + 2}), sz(0:{la + 1},0:{lb + 2})
      real(dp) :: tx(0:{la + 1},0:{lb}), ty(0:{la + 1},0:{lb}), tz(0:{la + 1},0:{lb})
      real(dp) :: vn({nvn},0:{lab1}), vn0({nvn})
      real(dp) :: wx({nvnw},0:{lab}), wy({nvnw},0:{lab}), wz({nvnw},0:{lab})
      real(dp) :: wx0({nvnw}), wy0({nvnw}), wz0({nvnw})
      real(dp) :: vAp({na1nb})
{decl_am}      real(dp) :: vBp({nanb1})
{decl_bm}      real(dp) :: wxt({na*nb}), wyt({na*nb}), wzt({na*nb})
      real(dp) :: f(0:BOYS_MMAX)

      abx = ra(1) - rb(1); aby = ra(2) - rb(2); abz = ra(3) - rb(3)
      rab2 = abx*abx + aby*aby + abz*abz

      ga = 0.0_dp; gb = 0.0_dp

      do ip = 1, nprim_a
         al = ea(ip)
      do jp = 1, nprim_b
         bet = eb(jp)
         zeta = al + bet
         oo2z = 0.5_dp/zeta
         xi = al*bet/zeta
         cc = ca(ip)*cb(jp)*wt
         tbb = 2.0_dp*bet*bet

         px = (al*ra(1) + bet*rb(1))/zeta
         py = (al*ra(2) + bet*rb(2))/zeta
         pz = (al*ra(3) + bet*rb(3))/zeta
         pax = px - ra(1); pay = py - ra(2); paz = pz - ra(3)
         pbx = px - rb(1); pby = py - rb(2); pbz = pz - rb(3)

         e0x = sqrt(PI/zeta)*exp(-xi*abx*abx)
         e0y = sqrt(PI/zeta)*exp(-xi*aby*aby)
         e0z = sqrt(PI/zeta)*exp(-xi*abz*abz)

{sxyz}

{txyz}

{st_block}

         ! --- nuclear attraction: basis derivative + Hellmann-Feynman, per nucleus ---
         vpref = 2.0_dp*PI/zeta*exp(-xi*rab2)
         do ic = 1, natm
            pcx = px - ratm(1, ic)
            pcy = py - ratm(2, ic)
            pcz = pz - ratm(3, ic)
            tval = zeta*(pcx*pcx + pcy*pcy + pcz*pcz)
            call boys_eval({lab1}, tval, f)

{vrr}

{dcvrr}

            do k = 1, {nvn}
               vn0(k) = vn(k,0)
            end do
            do k = 1, {nvnw}
               wx0(k) = wx(k,0); wy0(k) = wy(k,0); wz0(k) = wz(k,0)
            end do

{hrr_ap}

{hrr_am}

{hrr_bp}

{hrr_bm}

{hrr_wx}

{hrr_wy}

{hrr_wz}

{v_block}
         end do
      end do
      end do
   end subroutine one_e_grad_{tag}
"""


def emit_grad_module(L, path, order="lexicographic"):
    """The one-electron GRADIENT kernels, one routine per class, own module."""
    order_val = 0 if order == "lexicographic" else 1
    order_arg = "" if order == "lexicographic" else f" --order {order}"
    parts = [f"""!
! One-electron GRADIENT kernels: d/dR of overlap, kinetic and nuclear
! attraction, contracted directly with the energy-weighted density (S) and
! the density (T, V) into a (3, natm) gradient -- no derivative matrix is
! ever materialised. See scripts/gen_1e.py's emit_grad_class and
! emit_dc_vrr for the derivation (bra/ket shifted-angular-momentum
! relations; the Hellmann-Feynman VRR derivative wrt the nucleus).
!
! GENERATED by scripts/gen_1e.py --lmax {L} --deriv{order_arg} -- do not edit.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
module trc_1e_grad_kernels
   use pic_types, only: dp
   use trc_boys, only: BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV
   implicit none
   private
   public :: one_e_grad_dispatch, one_e_grad_driver, ONE_E_GRAD_LMAX

   integer, parameter :: ONE_E_GRAD_LMAX = {L}
   integer, parameter :: ONE_E_RADIX = __RADIX__
   !! Same fixed radix as gen_1e.py's other dispatchers -- see ONE_E_RADIX's
   !! own comment at the top of this file.
   real(dp), parameter :: PI = 3.14159265358979323846_dp
   integer, parameter :: NC_MAX = (ONE_E_GRAD_LMAX + 1)*(ONE_E_GRAD_LMAX + 2)/2
   integer, parameter :: NBLK_MAX = NC_MAX*NC_MAX
   integer, parameter :: ONE_E_GRAD_NPRIV = 8
   !! Lane-private gradient copies, same technique as trc_grad.F90's
   !! GRAD_NPRIV; fewer of them because this driver's work item is an
   !! entire shell pair (every primitive, every nucleus already summed
   !! inside one_e_grad_item) rather than one quartet, so there are far
   !! fewer concurrent writers to spread across lanes.

contains

#include "inc/trc_boys_eval.inc"
#include "inc/trc_1e_grad_driver.inc"
"""]
    names = []
    for la in range(L + 1):
        for lb in range(L + 1):
            parts.append(emit_grad_class(la, lb, order))
            names.append((la, lb))

    parts.append(f"""
   pure subroutine one_e_grad_dispatch(la, lb, nprim_a, nprim_b, ea, ca, eb, cb, &
                                       ra, rb, natm, zatm, ratm, wt, wblk, dblk, ga, gb, gc)
      !! Dispatch on the class. Host-side call, so the module boundary the
      !! kernels care about is not crossed here.
      !$acc routine seq
      integer,  intent(in)  :: la, lb, nprim_a, nprim_b, natm
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a)
      real(dp), intent(in)  :: eb(nprim_b), cb(nprim_b)
      real(dp), intent(in)  :: ra(3), rb(3)
      real(dp), intent(in)  :: zatm(natm), ratm(3, natm)
      real(dp), intent(in)  :: wt, wblk(*), dblk(*)
      real(dp), intent(out)   :: ga(3), gb(3)
      real(dp), intent(inout) :: gc(3, natm)
      select case (la*ONE_E_RADIX + lb)""")
    for la, lb in names:
        parts.append(
            f"""      case ({la*ONE_E_RADIX + lb}); call one_e_grad_{la}{lb}(nprim_a, nprim_b, ea, ca, eb, cb, &
                          ra, rb, natm, zatm, ratm, wt, wblk, dblk, ga, gb, gc)""")
    parts.append("""      case default
         !! No kernel for this class -- ONE_E_GRAD_LMAX is below the basis.
         !! A sentinel, as one_e_dispatch's: huge() rather than silence or a
         !! NaN (device code under -fast launders a NaN through the first
         !! min/max it meets). trc_basis refuses such a basis before this
         !! runs; this arm is the backstop for a path that skipped that.
         ga(1) = huge(1.0_dp)
      end select
   end subroutine one_e_grad_dispatch

end module trc_1e_grad_kernels
""")

    txt = "\n".join(parts).replace("__RADIX__", str(ONE_E_RADIX))
    with open(path, "w") as fh:
        fh.write(txt)
    print(f"wrote {path}: {len(names)} classes, {txt.count(chr(10))} lines")


# ---------------------------------------------------------------------------
# SPLIT EMITTERS: one module per (la, lb) class, each with its own launch.
#
# The single-file gradient and multipole modules put every class in one
# translation unit behind one device-side dispatcher, so one `do concurrent`
# could call any class and everything inlined. Above f that unit stopped
# compiling in reasonable time: at TERCO_LMAX=4 the gradient module was 114k
# lines and nvfortran's back end ran for over 90 minutes on it, alone, while
# the rest of the build waited (at LMAX=3: 413 s for it and 1165 s for the
# multipoles). Split as gen_df.py and gen_perclass.py already are: each class
# module launches its own `do concurrent` over its own pairs, so inlining
# stays inside the module and the classes compile in parallel. The driver
# module keeps the public interface trc_api uses and groups the pair list by
# class on the host.
# ---------------------------------------------------------------------------

def _ncart(l):
    return (l + 1) * (l + 2) // 2


def _guard(order):
    v = 0 if order == "lexicographic" else 1
    return (f"#if TRC_CART_ORDER != {v}\n"
            '#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."\n'
            "#endif\n")


GRAD_ARGS = ("nshell, nao, natm, maxnp, sh_l, sh_np, sh_ao, sh_at, sh_e, sh_c, sh_r, at_z, at_r, "
             "npair, p_i, p_j, wmat, dmat, nlane, gp")
GRAD_DECL = """      integer,  intent(in)  :: nshell, nao, natm, maxnp, npair, nlane
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell), sh_at(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      real(dp), intent(in)  :: at_z(natm), at_r(3, natm)
      real(dp), intent(in)  :: wmat(nao, nao), dmat(nao, nao)
      real(dp), intent(inout) :: gp(3, natm, nlane)"""


def emit_grad_split(L, outdir, order="lexicographic"):
    order_arg = "" if order == "lexicographic" else f" --order {order}"
    head = f"! GENERATED by scripts/gen_1e.py --lmax {L} --deriv-dir{order_arg} -- do not edit.\n"
    names = []
    for la in range(L + 1):
        for lb in range(L + 1):
            tag = f"{la}{lb}"
            nb = _ncart(la) * _ncart(lb)
            txt = f"""!
! One-electron gradient, class ({la}|{lb}): kernel, item and its own launch.
! Why one module per class: the split emitters in gen_1e.py.
!
{head}!
{_guard(order)}module trc_1e_grad_k{tag}
   use pic_types, only: dp
   use trc_boys, only: BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV
   implicit none
   private
   public :: one_e_grad_launch_{tag}
   real(dp), parameter :: PI = 3.14159265358979323846_dp

contains

#include "inc/trc_boys_eval.inc"

   subroutine one_e_grad_launch_{tag}(n, idx, {GRAD_ARGS})
      !! The pairs idx(1:n) of the pair list, all of class ({la}|{lb}).
      integer,  intent(in)  :: n, idx(n)
{GRAD_DECL}
      integer :: k
      do concurrent(k = 1:n)
         call one_e_grad_item_{tag}(idx(k), {GRAD_ARGS})
      end do
   end subroutine one_e_grad_launch_{tag}

   pure subroutine one_e_grad_item_{tag}(ij, {GRAD_ARGS})
      !$acc routine seq
      integer,  intent(in)  :: ij
{GRAD_DECL}
      integer  :: ii, jj, a, c, mu, nu, lane
      real(dp) :: wt, ga(3), gb(3)
      real(dp) :: wblk({nb}), dblk({nb})
      ii = p_i(ij); jj = p_j(ij)
      wt = merge(1.0_dp, 2.0_dp, ii == jj)
      lane = mod(ij - 1, nlane) + 1
      ! Affine index, not a running counter: nvfortran 26.5 fails "Missing
      ! increment value" on a counter used as a subscript in this loop nest
      ! (the same defect gen_df.py works around).
      do c = 0, {_ncart(lb) - 1}
         nu = sh_ao(jj) + c
         do a = 0, {_ncart(la) - 1}
            mu = sh_ao(ii) + a
            wblk(1 + a + {_ncart(la)}*c) = wmat(mu, nu)
            dblk(1 + a + {_ncart(la)}*c) = dmat(mu, nu)
         end do
      end do
      call one_e_grad_{tag}(sh_np(ii), sh_np(jj), sh_e(:, ii), sh_c(:, ii), sh_e(:, jj), sh_c(:, jj), &
                            sh_r(:, ii), sh_r(:, jj), natm, at_z, at_r, wt, wblk, dblk, ga, gb, &
                            gp(1:3, 1:natm, lane))
      ! ATOMIC, same reason as the gc write inside one_e_grad_{tag}: two
      ! shell pairs sharing BOTH a lane and an atom (sh_at(ii) == sh_at
      ! (jj) of some OTHER pair) run this read-modify-write concurrently.
      ! Measured on the GPU: water/6-31G* (few shell pairs, so few same-lane,
      ! same-atom collisions) passed with this non-atomic; water/cc-pVDZ
      ! and formamide/cc-pVDZ (more shells, more collisions) did not --
      ! wrong by 0.1 to 11 Eh/Bohr, the same race signature (small, not
      ! huge, because most pairs still land cleanly) as the gc case before
      ! ITS atomic fix. The OpenMP port at 4 threads shows it on water/6-31G*
      ! too. Scalar statements because OpenACC atomic update does not take a
      ! vector/array-section target.
      !$acc atomic update
      gp(1, sh_at(ii), lane) = gp(1, sh_at(ii), lane) + ga(1)
      !$acc atomic update
      gp(2, sh_at(ii), lane) = gp(2, sh_at(ii), lane) + ga(2)
      !$acc atomic update
      gp(3, sh_at(ii), lane) = gp(3, sh_at(ii), lane) + ga(3)
      !$acc atomic update
      gp(1, sh_at(jj), lane) = gp(1, sh_at(jj), lane) + gb(1)
      !$acc atomic update
      gp(2, sh_at(jj), lane) = gp(2, sh_at(jj), lane) + gb(2)
      !$acc atomic update
      gp(3, sh_at(jj), lane) = gp(3, sh_at(jj), lane) + gb(3)
   end subroutine one_e_grad_item_{tag}
{emit_grad_class(la, lb, order)}
end module trc_1e_grad_k{tag}
"""
            with open(os.path.join(outdir, f"trc_1e_grad_k{tag}.F90"), "w") as fh:
                fh.write(txt)
            names.append(tag)
    uses = "\n".join(f"   use trc_1e_grad_k{t}, only: one_e_grad_launch_{t}" for t in names)
    cases = "\n".join(f"""            case ({int(t[0]) * (L + 1) + int(t[1])})
               call one_e_grad_launch_{t}(cnt(c), idx(off(c) + 1), nshell, nao, natm, maxnp, sh_l, sh_np, &
                                        sh_ao, sh_at, sh_e, sh_c, sh_r, at_z, at_r, npair, p_i, p_j, &
                                        wmat, dmat, ONE_E_GRAD_NPRIV, gp)""" for t in names)
    drv = f"""!
! One-electron gradient driver over the per-class modules trc_1e_grad_k<la><lb>.
!
{head}!
module trc_1e_grad_kernels
   use pic_types, only: dp
{uses}
   implicit none
   private
   public :: one_e_grad_driver, ONE_E_GRAD_LMAX

   integer, parameter :: ONE_E_GRAD_LMAX = {L}
   integer, parameter :: ONE_E_GRAD_NPRIV = 8
   !! Lane-private gradient copies, as trc_grad.F90's GRAD_NPRIV.

contains

   subroutine one_e_grad_driver(nshell, nao, natm, maxnp, sh_l, sh_np, sh_ao, sh_at, &
                                sh_e, sh_c, sh_r, at_z, at_r, &
                                npair, p_i, p_j, wmat, dmat, grad)
      !! d/dR of S, T and V contracted with W and D into grad(3, natm). The
      !! pair list is grouped by class on the host (p_i, p_j and sh_l are read
      !! there) and each class runs its own launch.
      integer,  intent(in)  :: nshell, nao, natm, maxnp, npair
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell), sh_at(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      real(dp), intent(in)  :: at_z(natm), at_r(3, natm)
      real(dp), intent(in)  :: wmat(nao, nao), dmat(nao, nao)
      real(dp), intent(inout) :: grad(3, natm)
      integer, parameter :: NCLS = {(L + 1) ** 2}
      integer :: ij, c, ga_, gc_, gs, cnt(NCLS), off(NCLS), pos(NCLS)
      integer, allocatable :: idx(:)
      real(dp) :: gsum
      real(dp), allocatable :: gp(:, :, :)

      cnt = 0
      do ij = 1, npair
         c = sh_l(p_i(ij))*{L + 1} + sh_l(p_j(ij)) + 1
         cnt(c) = cnt(c) + 1
      end do
      off(1) = 0
      do c = 2, NCLS
         off(c) = off(c - 1) + cnt(c - 1)
      end do
      allocate (idx(max(npair, 1)))
      pos = off
      do ij = 1, npair
         c = sh_l(p_i(ij))*{L + 1} + sh_l(p_j(ij)) + 1
         pos(c) = pos(c) + 1
         idx(pos(c)) = ij
      end do
      allocate (gp(3, natm, ONE_E_GRAD_NPRIV))
      !$acc enter data copyin(idx) create(gp)
      do concurrent(gs=1:ONE_E_GRAD_NPRIV, ga_=1:natm, gc_=1:3)
         gp(gc_, ga_, gs) = 0.0_dp
      end do
      do c = 1, NCLS
         if (cnt(c) == 0) cycle
         select case (c - 1)
{cases}
         end select
      end do
      do concurrent(ga_=1:natm, gc_=1:3) local(gsum, gs)
         gsum = 0.0_dp
         do gs = 1, ONE_E_GRAD_NPRIV
            gsum = gsum + gp(gc_, ga_, gs)
         end do
         grad(gc_, ga_) = grad(gc_, ga_) + gsum
      end do
      !$acc exit data delete(idx, gp)
      deallocate (gp, idx)
   end subroutine one_e_grad_driver

end module trc_1e_grad_kernels
"""
    with open(os.path.join(outdir, "trc_1e_grad_kernels.F90"), "w") as fh:
        fh.write(drv)
    print(f"wrote {outdir}: {len(names)} gradient class modules + driver")


MULT_ARGS = "nshell, nao, maxnp, sh_l, sh_np, sh_ao, sh_e, sh_c, sh_r, npair, p_i, p_j, orig, mmat"
MULT_DECL = """      integer,  intent(in)  :: nshell, nao, maxnp, npair
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell), orig(3)
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      real(dp), intent(inout) :: mmat(nao, nao, 39)"""


def emit_multipole_split(L, outdir, order="lexicographic"):
    order_arg = "" if order == "lexicographic" else f" --order {order}"
    head = f"! GENERATED by scripts/gen_1e.py --lmax {L} --multipole-dir{order_arg} -- do not edit.\n"
    names = []
    for la in range(L + 1):
        for lb in range(L + 1):
            tag = f"{la}{lb}"
            nb = _ncart(la) * _ncart(lb)
            txt = f"""!
! Cartesian multipoles through the octupole, class ({la}|{lb}): kernel, item and
! its own launch. Why one module per class: the split emitters in gen_1e.py.
!
{head}!
{_guard(order)}module trc_mult_k{tag}
   use pic_types, only: dp
   implicit none
   private
   public :: mult_launch_{tag}

contains

   subroutine mult_launch_{tag}(n, idx, {MULT_ARGS})
      integer,  intent(in)  :: n, idx(n)
{MULT_DECL}
      integer :: k
      do concurrent(k = 1:n)
         call mult_item_{tag}(idx(k), {MULT_ARGS})
      end do
   end subroutine mult_launch_{tag}

   pure subroutine mult_item_{tag}(ij, {MULT_ARGS})
      !$acc routine seq
      integer,  intent(in)  :: ij
{MULT_DECL}
      integer  :: ii, jj, a, c, mu, nu, ic
      real(dp) :: mblk({nb}, 39)
      ii = p_i(ij); jj = p_j(ij)
      call mult_{tag}(sh_np(ii), sh_np(jj), sh_e(:, ii), sh_c(:, ii), sh_e(:, jj), sh_c(:, jj), &
                      sh_r(:, ii), sh_r(:, jj), orig, mblk)
      ! Affine index, not a running counter (nvfortran 26.5 "Missing increment
      ! value", as in the gradient items and gen_df.py).
      do c = 0, {_ncart(lb) - 1}
         nu = sh_ao(jj) + c
         do a = 0, {_ncart(la) - 1}
            mu = sh_ao(ii) + a
            do ic = 1, 39
               mmat(mu, nu, ic) = mblk(1 + a + {_ncart(la)}*c, ic)
               mmat(nu, mu, ic) = mblk(1 + a + {_ncart(la)}*c, ic)
            end do
         end do
      end do
   end subroutine mult_item_{tag}
{multipole_class_routine(la, lb, order)}
end module trc_mult_k{tag}
"""
            with open(os.path.join(outdir, f"trc_mult_k{tag}.F90"), "w") as fh:
                fh.write(txt.replace("PI_", "3.14159265358979323846_dp"))
            names.append(tag)
    uses = "\n".join(f"   use trc_mult_k{t}, only: mult_launch_{t}" for t in names)
    cases = "\n".join(f"""            case ({int(t[0]) * (L + 1) + int(t[1])})
               call mult_launch_{t}(cnt(c), idx(off(c) + 1), nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                                    sh_e, sh_c, sh_r, npair, p_i, p_j, orig, mmat)""" for t in names)
    drv = f"""!
! Multipole driver over the per-class modules trc_mult_k<la><lb>.
!
{head}!
module trc_mult_kernels
   use pic_types, only: dp
{uses}
   implicit none
   private
   public :: mult_driver, TRC_NMULT, MULT_LMAX

   integer, parameter :: MULT_LMAX = {L}
   integer, parameter :: TRC_NMULT = 39
   !! 3 dipole + 9 quadrupole + 27 octupole Cartesian components.

contains

   subroutine mult_driver(nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                          sh_e, sh_c, sh_r, npair, p_i, p_j, orig, mmat)
      !! Every multipole matrix over the pair list, grouped by class on the host.
      integer,  intent(in)  :: nshell, nao, maxnp, npair
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell), orig(3)
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      real(dp), intent(out) :: mmat(nao, nao, TRC_NMULT)
      integer, parameter :: NCLS = {(L + 1) ** 2}
      integer :: ij, c, i, j, k, cnt(NCLS), off(NCLS), pos(NCLS)
      integer, allocatable :: idx(:)

      do concurrent(k=1:TRC_NMULT, j=1:nao, i=1:nao)
         mmat(i, j, k) = 0.0_dp
      end do
      cnt = 0
      do ij = 1, npair
         c = sh_l(p_i(ij))*{L + 1} + sh_l(p_j(ij)) + 1
         cnt(c) = cnt(c) + 1
      end do
      off(1) = 0
      do c = 2, NCLS
         off(c) = off(c - 1) + cnt(c - 1)
      end do
      allocate (idx(max(npair, 1)))
      pos = off
      do ij = 1, npair
         c = sh_l(p_i(ij))*{L + 1} + sh_l(p_j(ij)) + 1
         pos(c) = pos(c) + 1
         idx(pos(c)) = ij
      end do
      !$acc enter data copyin(idx)
      do c = 1, NCLS
         if (cnt(c) == 0) cycle
         select case (c - 1)
{cases}
         end select
      end do
      !$acc exit data delete(idx)
      deallocate (idx)
   end subroutine mult_driver

end module trc_mult_kernels
"""
    with open(os.path.join(outdir, "trc_mult_kernels.F90"), "w") as fh:
        fh.write(drv)
    print(f"wrote {outdir}: {len(names)} multipole class modules + driver")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=2)
    ap.add_argument("-o", "--output", default="src/trc_1e_kernels.F90")
    ap.add_argument("--order", choices=("lexicographic", "gamess"),
                    default="lexicographic",
                    help="Cartesian order of the AO slots sout/tout/vout (and "
                         "the multipole AO axis, under --multipole) are "
                         "written in. INTERNAL VRR indices (cart_cum, "
                         "emit_nuc_vrr) stay lexicographic regardless.")
    ap.add_argument("--multipole", metavar="FILE", default=None,
                    help="emit the multipole kernels into their own module "
                         "instead. They are NOT fused with S/T/V: doing that "
                         "put 39 components per class into the same device "
                         "routine and NVVM ground for over thirty minutes at "
                         "LMAX=2. Both run once per geometry, so sharing the "
                         "1-D tables saved nothing worth that.")
    ap.add_argument("--pcm", metavar="FILE", default=None,
                    help="emit the PCM smeared-charge surface-integral "
                         "kernels (erf-attenuated nuclear attraction, one "
                         "routine per class) into their own module instead.")
    ap.add_argument("--deriv", metavar="FILE", default=None,
                    help="emit the one-electron GRADIENT kernels (dS, dT, "
                         "dV -- bra/ket basis derivative and the nuclear "
                         "Hellmann-Feynman term -- contracted directly "
                         "into a (3,natm) gradient) into their own module "
                         "instead.")
    ap.add_argument("--deriv-dir", metavar="DIR", default=None,
                    help="as --deriv, one module per class plus the driver module "
                         "(trc_1e_grad_kernels.F90) in DIR")
    ap.add_argument("--multipole-dir", metavar="DIR", default=None,
                    help="as --multipole, one module per class plus the driver "
                         "module (trc_mult_kernels.F90) in DIR")
    args = ap.parse_args()
    L = args.lmax
    order = args.order
    order_val = 0 if order == "lexicographic" else 1
    order_arg = "" if order == "lexicographic" else f" --order {order}"

    if args.deriv_dir:
        os.makedirs(args.deriv_dir, exist_ok=True)
        emit_grad_split(L, args.deriv_dir, order)
        return
    if args.multipole_dir:
        os.makedirs(args.multipole_dir, exist_ok=True)
        emit_multipole_split(L, args.multipole_dir, order)
        return
    if args.multipole:
        emit_multipole_module(L, args.multipole, order)
        return
    if args.pcm:
        emit_pcm_module(L, args.pcm, order)
        return
    if args.deriv:
        emit_grad_module(L, args.deriv, order)
        return

    parts = [f"""!
! One-electron kernels: overlap, kinetic, nuclear attraction.
!
! GENERATED by scripts/gen_1e.py --lmax {L}{order_arg} -- do not edit.
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
! TRC_CART_ORDER GUARD. sout/tout/vout below are AO-indexed and follow the
! order --order was given at generation time (0 = lexicographic, the
! default; 1 = gamess); TRC_CART_ORDER selects it at compile time, and an
! undefined TRC_CART_ORDER reads as 0 in a preprocessor #if, so the default
! build needs no new -D flag. Same idea as docs/GENERATORS.md's TRC_LMAX
! guard, for a different compile-time constant.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
module trc_1e_kernels
   use pic_types, only: dp
   use trc_boys, only: BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV
   implicit none
   private
   public :: one_e_dispatch, one_e_driver, ONE_E_LMAX, TRC_NMULT

   integer, parameter :: TRC_NMULT = 3 + 9 + 27
   !! Dipole + quadrupole + octupole, full Cartesian tensors.

   integer, parameter :: ONE_E_LMAX = {L}
   integer, parameter :: ONE_E_RADIX = __RADIX__
   !! Dispatch radix, FIXED and independent of ONE_E_LMAX -- see gen_1e.py.
   real(dp), parameter :: PI = 3.14159265358979323846_dp
   integer, parameter :: NC_MAX = (ONE_E_LMAX + 1)*(ONE_E_LMAX + 2)/2
   !! Largest Cartesian block a shell pair produces here; the driver's
   !! per-item buffers are this size because `do concurrent` locals are fixed.
   integer, parameter :: NBLK_MAX = NC_MAX*NC_MAX

contains

#include "inc/trc_boys_eval.inc"
#include "inc/trc_1e_driver.inc"
"""]

    names = []
    for la in range(L + 1):
        for lb in range(L + 1):
            parts.append(emit_class(la, lb, L, order))
            names.append((la, lb))

    parts.append(f"""
   pure subroutine one_e_dispatch(la, lb, nprim_a, nprim_b, ea, ca, eb, cb, &
                                  ra, rb, natm, zatm, ratm, sout, tout, vout)
      !! Dispatch on the class. Host-side call, so the module boundary the
      !! kernels care about is not crossed here.
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
