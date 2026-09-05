#!/usr/bin/env python3
"""
Emit two- and three-centre Coulomb kernels for density fitting.

    (a|c)   = int int  chi_a(r1) (1/r12) chi_c(r2)
    (ab|c)  = int int  phi_a(r1) phi_b(r1) (1/r12) chi_c(r2)

WHY THIS IS THE CHEAP ONE
-------------------------
Nothing new is needed. The four-centre VRR already builds [a0|c0]^(m), which
IS the two-centre quantity; the four-centre form only differs by having a
Gaussian PRODUCT on each side instead of a single function. A single shell is
formally a product with a unit function of zero exponent, so:

    2-centre:  zeta = alpha,      P = A,  K_ab = 1
               eta  = gamma,      Q = C,  K_cd = 1
    3-centre:  zeta = alpha+beta, P as usual, K_ab = exp(-xi |AB|^2)
               eta  = gamma,      Q = C,  K_cd = 1

and the base case is the same one the ERI uses,

    [00|00]^(m) = 2 pi^(5/2) / (zeta eta sqrt(zeta+eta)) K_ab K_cd F_m(rho R_PQ^2)

with rho = zeta*eta/(zeta+eta). The recurrence coefficients are untouched. The
only structural difference is that the ket carries no angular momentum to
transfer, so the ket HRR disappears and only the bra HRR runs (3-centre), or
neither does (2-centre).

So this reuses gen_vrr's sieved emitter verbatim -- including Gill's sieve --
and supplies its own target set: the ket index is pinned to |c| = lc rather
than ranging over [lc, lc+ld].

These produce INTEGRALS, not a Fock matrix. Density fitting wants the (mn|P)
tensor itself, so there is no digestion here.

USAGE
    gen_df.py --lmax 2 --lmax-aux 3 -o src/trc_df_kernels.F90
"""

import argparse
import importlib.util
import os


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HERE = os.path.dirname(os.path.abspath(__file__))
gv = _load("gv", os.path.join(HERE, "gen_vrr.py"))
gh = _load("gh", os.path.join(HERE, "gen_hrr.py"))


def cart(l):
    return [(lx, ly, l - lx - ly)
            for lx in range(l, -1, -1)
            for ly in range(l - lx, -1, -1)]


def binom(n, k):
    r = 1
    for i in range(1, k + 1):
        r = r*(n - k + i)//i
    return r


def powexpr(var, n):
    return "*".join([var]*n) if n else ""


def emit_bra_hrr(la, lb, lc, I, nca):
    """(ab|c) from [e0|c0], e in [la, la+lb]. Ket untouched: ld = 0."""
    o, w = [], None
    w = o.append
    na, nb = len(cart(la)), len(cart(lb))
    for ic, cp in enumerate(cart(lc)):
        jc = I[cp]
        for ib, (bx, by, bz) in enumerate(cart(lb)):
            for ia, (ax, ay, az) in enumerate(cart(la)):
                terms = []
                for px in range(bx + 1):
                    for py in range(by + 1):
                        for pz in range(bz + 1):
                            cbin = binom(bx, px)*binom(by, py)*binom(bz, pz)
                            ja = I[(ax + px, ay + py, az + pz)]
                            fac = [f for f in (powexpr("abx", bx - px),
                                               powexpr("aby", by - py),
                                               powexpr("abz", bz - pz)) if f]
                            head = "" if cbin == 1 else f"{float(cbin)}_dp*"
                            body = "*".join(fac)
                            if body:
                                body += "*"
                            terms.append(f"{head}{body}g({ja + (jc - 1)*nca})")
                expr = terms[0]
                for t in terms[1:]:
                    expr += f" &\n                  + {t}"
                slot = 1 + ia + na*(ib + nb*ic)
                w(f"               out({slot}) = out({slot}) + {expr}")
    return "\n".join(o)


def emit_class(la, lb, lc, three):
    """One class. `three` selects (ab|c); otherwise (a|c)."""
    lab = la + lb
    lt = lab + lc
    nca, ncc = gv.ncum(lab), gv.ncum(lc)
    na, nb, nc = len(cart(la)), len(cart(lb)), len(cart(lc))
    nout = na*nb*nc if three else na*nc
    nv = nca*ncc

    P = [None] + gv.cart_list(max(lab, lc))
    I = {p: i for i, p in enumerate(P) if p is not None}

    # Targets: the ket index is PINNED to |c| = lc, since there is no ket
    # transfer. That is the whole difference from the four-centre target set,
    # and it lets the sieve prune harder.
    kets = [I[p] for p in cart(lc)]
    if three:
        bras = [I[p] for l in range(la, lab + 1) for p in cart(l)]
    else:
        bras = [I[p] for p in cart(la)]
    targets = {(b, c, 0) for b in bras for c in kets}

    vrr, fin, nlive = gv.emit_class_sieved(lab, lc, targets, P, I, nca)
    vrr = "\n".join("         " + ln if ln.strip() else ln
                    for ln in vrr.split("\n"))

    if three:
        tag = f"3c{la}{lb}{lc}"
        body = emit_bra_hrr(la, lb, lc, I, nca)
        #
        # The bra arrives as PRECOMPUTED primitive-pair data, not as two sets
        # of exponents. This kernel is called once per (shell pair, auxiliary
        # shell); rebuilding zeta, P and K_ab from the exponents here would
        # redo the bra product once per auxiliary SHELL, which for a real
        # fitting basis is hundreds of times over. The pair list computes it
        # once per shell pair and this reads it.
        #
        sig = ("npp, ppz, ppk, ppp, ppa, abx, aby, abz, "
               "nprim_c, ec, cc_, rc, out")
        decls = """      integer,  intent(in)  :: npp, nprim_c
      real(dp), intent(in)  :: ppz(npp), ppk(npp), ppp(3, npp), ppa(3, npp)
      real(dp), intent(in)  :: abx, aby, abz
      real(dp), intent(in)  :: ec(nprim_c), cc_(nprim_c), rc(3)"""
        bra_loop_open = """      do ip = 1, npp
         zeta = ppz(ip)
         kab  = ppk(ip)
         px = ppp(1, ip); py = ppp(2, ip); pz = ppp(3, ip)
         pax = ppa(1, ip); pay = ppa(2, ip); paz = ppa(3, ip)"""
        bra_loop_close = ""
        abset = ""
    else:
        tag = f"2c{la}{lc}"
        # No bra transfer either: the surviving [a0|c0] IS the answer.
        o = []
        for ic, cp in enumerate(cart(lc)):
            jc = I[cp]
            for ia, ap in enumerate(cart(la)):
                ja = I[ap]
                slot = 1 + ia + na*ic
                o.append(f"               out({slot}) = out({slot})"
                         f" + g({ja + (jc - 1)*nca})")
        body = "\n".join(o)
        sig_extra = ""
        decl_extra = ""
        # Two-centre keeps the plain exponent interface: there is no pair to
        # precompute, and it is called once per (shell, shell) rather than
        # once per (pair, shell), so nothing is recomputed.
        sig = "nprim_a, ea, ca, ra, nprim_c, ec, cc_, rc, out"
        decls = """      integer,  intent(in)  :: nprim_a, nprim_c
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a), ra(3)
      real(dp), intent(in)  :: ec(nprim_c), cc_(nprim_c), rc(3)"""
        bra_loop_open = """      do ip = 1, nprim_a
         al = ea(ip)
         zeta = al
         kab = ca(ip)
         px = ra(1); py = ra(2); pz = ra(3)
         pax = 0.0_dp; pay = 0.0_dp; paz = 0.0_dp"""
        bra_loop_close = ""
        abset = ""   # no HRR in the two-centre case, so no AB vector

    return f"""
   !> {'(ab|c)' if three else '(a|c)'} class ({la}{lb}|{lc}) -- {nlive} VRR statements after the sieve.
   pure subroutine df_{tag}({sig})
      !$acc routine seq
{decls}
      real(dp), intent(out) :: out({nout})

      integer  :: ip, kp, k, cur
      real(dp) :: al, gam, zeta, eta, zpe, rho, kab, pref
      real(dp) :: px, py, pz, pqx, pqy, pqz, tval, wc
      real(dp) :: pax, pay, paz, qcx, qcy, qcz
      real(dp) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re
      real(dp) :: v({nv},0:1), g({nv})
      real(dp) :: f(0:BOYS_MMAX)
{abset}

      do k = 1, {nout}
         out(k) = 0.0_dp
      end do

{bra_loop_open}

         do k = 1, {nv}
            g(k) = 0.0_dp
         end do

         do kp = 1, nprim_c
            gam = ec(kp)
            eta = gam
            zpe = zeta + eta
            rho = zeta*eta/zpe
            ! The ket is a single shell, so Q = C and QC = 0. Kept as
            ! variables rather than folded away because the VRR is emitted
            ! by the four-centre generator and expects them.
            qcx = 0.0_dp; qcy = 0.0_dp; qcz = 0.0_dp
            pqx = px - rc(1); pqy = py - rc(2); pqz = pz - rc(3)
            wc = (zeta*px + eta*rc(1))/zpe; wpx = wc - px; wqx = wc - rc(1)
            wc = (zeta*py + eta*rc(2))/zpe; wpy = wc - py; wqy = wc - rc(2)
            wc = (zeta*pz + eta*rc(3))/zpe; wpz = wc - pz; wqz = wc - rc(3)
            oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = 0.5_dp/zpe
            rz = rho/zeta; re = rho/eta
            tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)
            ! K_cd = 1: a single Gaussian has no product exponential.
            pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*kab*cc_(kp)
            call boys_eval({lt}, tval, f)

{vrr}

            do k = 1, {nv}
               g(k) = g(k) + v(k,{fin})
            end do
         end do

{body}
      end do
{bra_loop_close}
   end subroutine df_{tag}
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=2,
                    help="max l on the orbital (bra) shells")
    # Default one above the orbital ceiling: a fitting basis normally
    # outranks the orbital basis by one, and the shipped set (--lmax 2) was
    # built with f on the auxiliary shell. CI and CMake regenerate with
    # --lmax alone, so the default IS the shipped configuration; a default
    # of 2 here once made CI's regeneration disagree with the committed
    # kernels by nine thousand lines while nothing was actually wrong.
    ap.add_argument("--lmax-aux", type=int, default=None,
                    help="max l on the auxiliary shell (default: lmax + 1)")
    ap.add_argument("-o", "--output", default="src/trc_df_kernels.F90")
    args = ap.parse_args()
    L = args.lmax
    LA = args.lmax_aux if args.lmax_aux is not None else L + 1

    parts = [f"""!
! Two- and three-centre Coulomb integrals for density fitting.
!
! GENERATED by scripts/gen_df.py --lmax {L} --lmax-aux {LA} -- do not edit.
!
! These reuse the four-centre VRR unchanged, including Gill's sieve: a single
! shell is a Gaussian product with a unit function of zero exponent, so the
! only differences are eta = gamma, Q = C, and K_cd = 1. The ket carries no
! momentum to transfer, so the ket HRR disappears; the bra HRR is the usual
! one for (ab|c) and absent for (a|c). See scripts/gen_df.py.
!
! Output is the integral block itself, not a Fock contribution -- density
! fitting wants the tensor.
!
module trc_df_kernels
   use trc_boys, only: dp, boys_eval, BOYS_MMAX
   implicit none
   private
   public :: df_2c_dispatch, df_3c_dispatch, DF_LMAX, DF_LMAX_AUX

   integer, parameter :: DF_LMAX = {L}, DF_LMAX_AUX = {LA}
   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp

contains
"""]

    two, three = [], []
    # Both sides of (a|c) are auxiliary shells: the metric (P|Q) needs the
    # bra up to LA as well, and a missing (f|f) case left the block unset and
    # the metric indefinite the first time an RI set with f functions was
    # tried.
    for la in range(LA + 1):
        for lc in range(LA + 1):
            parts.append(emit_class(la, 0, lc, three=False))
            two.append((la, lc))
    for la in range(L + 1):
        for lb in range(L + 1):
            for lc in range(LA + 1):
                parts.append(emit_class(la, lb, lc, three=True))
                three.append((la, lb, lc))

    parts.append(f"""
   pure subroutine df_2c_dispatch(la, lc, nprim_a, ea, ca, ra, &
                                  nprim_c, ec, cc_, rc, out)
      !$acc routine seq
      integer,  intent(in)  :: la, lc, nprim_a, nprim_c
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a), ra(3)
      real(dp), intent(in)  :: ec(nprim_c), cc_(nprim_c), rc(3)
      real(dp), intent(out) :: out(*)
      select case (la*({LA} + 1) + lc)""")
    for la, lc in two:
        parts.append(f"""      case ({la*(LA + 1) + lc}); call df_2c{la}{lc}(nprim_a, ea, ca, ra, &
                          nprim_c, ec, cc_, rc, out)""")
    parts.append("""      end select
   end subroutine df_2c_dispatch
""")

    parts.append(f"""
   pure subroutine df_3c_dispatch(la, lb, lc, npp, ppz, ppk, ppp, ppa, &
                                  abx, aby, abz, nprim_c, ec, cc_, rc, out)
      !$acc routine seq
      integer,  intent(in)  :: la, lb, lc, npp, nprim_c
      real(dp), intent(in)  :: ppz(npp), ppk(npp), ppp(3, npp), ppa(3, npp)
      real(dp), intent(in)  :: abx, aby, abz
      real(dp), intent(in)  :: ec(nprim_c), cc_(nprim_c), rc(3)
      real(dp), intent(out) :: out(*)
      select case ((la*({L} + 1) + lb)*({LA} + 1) + lc)""")
    for la, lb, lc in three:
        parts.append(
            f"""      case ({(la*(L + 1) + lb)*(LA + 1) + lc}); call df_3c{la}{lb}{lc}(npp, ppz, ppk, ppp, ppa, &
                          abx, aby, abz, nprim_c, ec, cc_, rc, out)""")
    parts.append("""      end select
   end subroutine df_3c_dispatch

end module trc_df_kernels
""")

    txt = "\n".join(parts)
    with open(args.output, "w") as fh:
        fh.write(txt)
    print(f"wrote {args.output}: {len(two)} two-centre + {len(three)} "
          f"three-centre classes, {txt.count(chr(10))} lines")


if __name__ == "__main__":
    main()
