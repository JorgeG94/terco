#!/usr/bin/env python3
"""Generate the McMurchie-Davidson Coulomb GRADIENT kernels.

d/dX Tr(D J[D]) = 2 sum_ab D_ab d_X (ab | rho) needs the bra side only
(derivation and numpy check: scripts/md/reference_grad.py). Beside the energy
engine's ket stage (unchanged, scripts/gen_md_j.py) it adds

  1. the Hermite stage with the BRA order raised by one,
         Jt'[P]_t = sum_Q sum_tau R_{t + tau}(P, Q) Dt[Q]_tau,  |t| <= lp + 1:
     the energy engine's gather kernel [lp + 1 | lq], emitted by
     gen_md_j.gen_r_module_gather(..., grad=True) as trc_mdjg_h<lp+1><lq>, one
     launch per R batch, a longer Jt' array of its own;
  2. a per-shell-pair gradient bra (trc_mdjg_e<la><lb>, mdjg_bra_<la><lb>),
     one thread per shell pair, which contracts D_ab with the derivative
     distributions against Jt' and adds three components to each atom of the
     pair, by atomics:
         gA_x = c_P ( 2 alpha sum_t E^{a+1_x,b}_t Jt'_t - a_x sum_t E^{a-1_x,b}_t Jt'_t )
         gT_x = c_P   sum_t E^{ab}_t Jt'_{t + 1_x}                     (d/dA + d/dB)
         grad[atom A] += 2 D' gA,  grad[atom B] += 2 D' (gT - gA);
  3. the screen's bound for the DERIVATIVE pair distribution (mdjg_bound_*):
     I'_P = |c_P| max over component pairs, directions and the three
     distributions an atom can accumulate, d/dA, d/dB and d/dA + d/dB (the
     last is what an atom carrying BOTH shells of the pair gets), of the
     Coulomb norm. By Cauchy-Schwarz in the Coulomb metric, one (P, Q)'s
     contribution to one component of one atom's gradient is at most
     2 |D'_ab| I'_P ||rho_Q|| per component pair: the factor 2 of w = 2 D'
     and the density |D'| are handled as in the energy engine (the screen is
     density-free, and the pair's cut carries the 1e-3 margin), and the max
     over the three distributions makes it hold per atom, whichever atom
     collects it. It needs the exponent alpha of each primitive pair, which
     the energy engine does not keep.

The E coefficients are recomputed in every item (e1d_loop of gen_md_k.py with
la + 1 on the first index: the a - 1, a and a + 1 tables come out of one
recursion).

Usage:
    gen_md_j_grad.py --lmax 2 [--order lexicographic] [--batches 7:6] -o DIR

The kernels are TERCO_LMAX-independent beyond --lmax (the highest shell l the
gradient covers; CMake passes min(TERCO_LMAX, 3), the Rys gradient's own
validated ceiling), so one shell above it makes the engine fall back to the
fused Rys gradient.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import gen_md_j as G  # noqa: E402
from gen_md_k import e1d_loop, int_table  # noqa: E402
from md.hermite import herm  # noqa: E402

BOYS_MMAX = 20   # src/core/trc_boys.F90


def norm_tables(hs):
    """(h, h') pairs of the Coulomb self-norm, h <= h', all-even sums: PI, PJ, PS, PW."""
    nh = len(hs)
    pi, pj, ps, pw = [], [], [], []
    for i, h in enumerate(hs):
        for j in range(i, nh):
            h2 = hs[j]
            r = tuple(h[d] + h2[d] for d in range(3))
            if any(x % 2 for x in r):
                continue
            w = (1 if i == j else 2) * G.dfact(r[0]) * G.dfact(r[1]) * G.dfact(r[2])
            if sum(h) % 2:
                w = -w
            pi.append(i + 1)
            pj.append(j + 1)
            ps.append(sum(r) // 2)
            pw.append(w)
    return pi, pj, ps, pw


def real_tab(name, vals):
    toks = [G.lit(float(v)) for v in vals]
    lines, cur = [], f"   real(dp), parameter :: {name}({len(vals)}) = ["
    for tok in toks:
        if len(cur) + len(tok) > 100:
            lines.append(cur + " &")
            cur = "      "
        cur += tok + ", "
    lines.append(cur.rstrip(", ") + "]")
    return lines


BRA_ARGS = ("n, s0, kl0, nkc, hb, nhp, mp_z, mp_ea, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, "
            "nsp, sp_ao1, sp_ao2, sp_k0, sp_nk, sp_at1, sp_at2, ngt, gjt, nao, dmat, natm, grad")
BND_ARGS = "n, k0, nhp, mp_z, mp_ea, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, mp_gip"

COMMON_DECL = """      real(dp), intent(in) :: mp_z(nhp), mp_ea(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)
      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)"""
BRA_DECL = f"""      integer,  intent(in) :: n, s0, kl0, nkc, hb, nhp, nsp, ngt, nao, natm
{COMMON_DECL}
      integer,  intent(in) :: sp_ao1(nsp), sp_ao2(nsp), sp_k0(nsp), sp_nk(nsp), sp_at1(nsp), sp_at2(nsp)
      real(dp), intent(in) :: gjt(ngt), dmat(nao, nao)
      real(dp), intent(inout) :: grad(3, natm)"""
BND_DECL = f"""      integer,  intent(in) :: n, k0, nhp
{COMMON_DECL}
      real(dp), intent(inout) :: mp_gip(nhp)"""


EE = "ee = e(ia(1), ib(1), HT(hh), 1)*e(ia(2), ib(2), HU(hh), 2)*e(ia(3), ib(3), HV(hh), 3)"


def gen_e_module(la, lb, order, header):
    """Gradient bra and derivative-bound stages of one (la, lb) type, as loops over tables."""
    lp = la + lb
    lp1 = lp + 1
    tag = f"{la}{lb}"
    h0, h1 = herm(lp), herm(lp1)
    nh0, nh1 = len(h0), len(h1)
    idx1 = {h: i + 1 for i, h in enumerate(h1)}
    ca, cb = G.comps(la, order), G.comps(lb, order)
    na, nb = len(ca), len(cb)
    pi, pj, ps, pw = norm_tables(h1)
    o = [header.format(what=f"Class ({la}{lb}) gradient bra and derivative bound of the MD J engine")]
    o.append('#include "inc/trc_dc_index.inc"')
    o.append(f"module trc_mdjg_e{tag}")
    o.append("   use pic_types, only: dp")
    o.append("   implicit none")
    o.append("   private")
    o.append(f"   public :: mdjg_bra_{tag}, mdjg_bound_{tag}")
    o.append(f"   real(dp), parameter :: TWO_PI_2_5 = {G.TWO_PI_2_5}")
    o.extend(int_table("HT", [h[0] for h in h1]))
    o.extend(int_table("HU", [h[1] for h in h1]))
    o.extend(int_table("HV", [h[2] for h in h1]))
    # Slot of h + 1_d in herm(lp + 1), for h in herm(lp) (a prefix of herm(lp + 1)).
    for d, ax in enumerate("XYZ"):
        o.extend(int_table(f"HS{ax}", [idx1[tuple(h[k] + (k == d) for k in range(3))] for h in h0]))
    for nm, cs in (("A", ca), ("B", cb)):
        for d, ax in enumerate("XYZ"):
            o.extend(int_table(f"C{nm}{ax}", [c[d] for c in cs]))
    o.extend(int_table("PI", pi))
    o.extend(int_table("PJ", pj))
    o.extend(int_table("PS", ps))
    o.extend(real_tab("PW", pw))
    o.extend(real_tab("RO", [1.0 / (2 * s_ + 1) for s_ in range(lp1 + 1)]))
    o.append("")
    o.append("contains")
    o.append("")
    for kind, args, dcl in (("bra", BRA_ARGS, BRA_DECL), ("bound", BND_ARGS, BND_DECL)):
        o.extend(G.wrap_args(f"subroutine mdjg_{kind}_{tag}", args, "   "))
        if kind == "bra":
            o.append(f"      !! Gradient of the {la}{lb} shell pairs s0+1..s0+n: d/dX Tr(D J[D]) into")
            o.append("      !! grad(3, natm) by atomics, from the primitive pairs' Jt' (Hermite order")
            o.append(f"      !! up to {lp1}). `kl0` is the first primitive of the pairs' lp block, less one.")
        else:
            o.append(f"      !! I'_P = |c_P| max ||.|| over the {la}{lb} primitive pairs k0+1..k0+n of the")
            o.append("      !! Coulomb norm of d/dA, d/dB and d/dA + d/dB of the pair distribution.")
        o.append(dcl)
        o.append("      integer :: t")
        o.append("      do concurrent(t = 1:n)")
        o.extend(G.wrap_args(f"call mdjg_{kind}_item_{tag}", "TRC_DC_INDEX(t), " + args, " " * 9))
        o.append("      end do")
        o.append(f"   end subroutine mdjg_{kind}_{tag}")
        o.append("")
    def e_setup(ind):
        return [f"{ind}oo2p = 0.5_dp/p",
                f"{ind}xa(1) = mp_ax(k); xa(2) = mp_ay(k); xa(3) = mp_az(k)",
                f"{ind}xb(1) = mp_bx(k); xb(2) = mp_by(k); xb(3) = mp_bz(k)"]

    def comp(ind):
        return [f"{ind}ia(1) = CAX(x1); ia(2) = CAY(x1); ia(3) = CAZ(x1)",
                f"{ind}ib(1) = CBX(x2); ib(2) = CBY(x2); ib(3) = CBZ(x2)"]
    locals_ = [f"      real(dp) :: e(0:{la + 1}, 0:{lb}, 0:{lp1}, 3), xa(3), xb(3), p, oo2p, v, ea",
               "      integer :: ia(3), ib(3), x1, x2, hh, d, ii, jj, tt, k"]

    # ---- bra item ----
    o.extend(G.wrap_args(f"pure subroutine mdjg_bra_item_{tag}", "i, " + BRA_ARGS, "   "))
    o.append("      !$acc routine seq")
    o.append("      integer,  intent(in) :: i")
    o.append(BRA_DECL)
    o.extend(locals_)
    o.append(f"      real(dp) :: jh({nh1}), ga(3), gt(3), su, sd, st, dab, cs, offd, ad, ee")
    o.append("      integer :: s, mu0, nu0, kk, pos, at1, at2, hs")
    o.append("      s = s0 + i")
    o.append("      mu0 = sp_ao1(s) - 1")
    o.append("      nu0 = sp_ao2(s) - 1")
    o.append("      ! One shell pair stands for both orders off the diagonal: D_ab + D_ba.")
    o.append("      offd = 1.0_dp")
    o.append("      if (mu0 == nu0) offd = 0.0_dp")
    o.append("      at1 = sp_at1(s)")
    o.append("      at2 = sp_at2(s)")
    o.append("      ga = 0.0_dp")
    o.append("      gt = 0.0_dp")
    o.append("      do kk = 1, sp_nk(s)")
    o.append("         k = sp_k0(s) + kk")
    o.append("         p = mp_z(k)")
    o.append("         ea = mp_ea(k)")
    o.extend(e_setup("         "))
    o.extend(e1d_loop(la + 1, lb, "         "))
    o.append("         pos = k - kl0")
    o.append("         cs = mp_c(k)")
    o.append(f"         do hh = 1, {nh1}")
    o.append("            jh(hh) = cs*gjt(hb + (hh - 1)*nkc + pos)")
    o.append("         end do")
    o.append(f"         do x2 = 1, {nb}")
    o.append(f"            do x1 = 1, {na}")
    o.extend(comp("               "))
    o.append("               dab = dmat(mu0 + x1, nu0 + x2) + offd*dmat(nu0 + x2, mu0 + x1)")
    o.append("               do d = 1, 3")
    o.append("                  ad = real(ia(d), dp)")
    o.append("                  ! d/dA of the Cartesian factor: 2 alpha (a + 1_d) - a_d (a - 1_d).")
    o.append("                  ia(d) = ia(d) + 1")
    o.append("                  su = 0.0_dp")
    o.append(f"                  do hh = 1, {nh1}")
    o.append("                     " + EE)
    o.append("                     su = su + ee*jh(hh)")
    o.append("                  end do")
    o.append("                  ia(d) = ia(d) - 2")
    o.append("                  sd = 0.0_dp")
    o.append("                  if (ia(d) >= 0) then")
    o.append(f"                     do hh = 1, {nh1}")
    o.append("                        " + EE)
    o.append("                        sd = sd + ee*jh(hh)")
    o.append("                     end do")
    o.append("                  end if")
    o.append("                  ia(d) = ia(d) + 1")
    o.append("                  ga(d) = ga(d) + dab*(2.0_dp*ea*su - ad*sd)")
    o.append("                  ! d/dA + d/dB: the Hermite index of the pair distribution up by one.")
    o.append("                  st = 0.0_dp")
    o.append(f"                  do hh = 1, {nh0}")
    o.append("                     if (d == 1) then")
    o.append("                        hs = HSX(hh)")
    o.append("                     else if (d == 2) then")
    o.append("                        hs = HSY(hh)")
    o.append("                     else")
    o.append("                        hs = HSZ(hh)")
    o.append("                     end if")
    o.append("                     " + EE)
    o.append("                     st = st + ee*jh(hs)")
    o.append("                  end do")
    o.append("                  gt(d) = gt(d) + dab*st")
    o.append("               end do")
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    o.append("      ! The factor 2 of d/dX Tr(D J[D]) = 2 sum D_ab d_X (ab|rho). Atoms are shared")
    o.append("      ! between shell pairs: atomics (a pair on one atom adds gA + gT - gA there).")
    o.append("      do d = 1, 3")
    o.append("         !$acc atomic update")
    o.append("         grad(d, at1) = grad(d, at1) + 2.0_dp*ga(d)")
    o.append("         !$acc atomic update")
    o.append("         grad(d, at2) = grad(d, at2) + 2.0_dp*(gt(d) - ga(d))")
    o.append("      end do")
    o.append(f"   end subroutine mdjg_bra_item_{tag}")
    o.append("")
    # ---- bound item ----
    o.extend(G.wrap_args(f"pure subroutine mdjg_bound_item_{tag}", "i, " + BND_ARGS, "   "))
    o.append("      !$acc routine seq")
    o.append("      integer,  intent(in) :: i")
    o.append(BND_DECL)
    o.extend(locals_)
    o.append(f"      real(dp) :: g({nh1}), acc(0:{lp1}), best, nq2, ws, ad, ee, gs({nh1})")
    o.append("      integer :: q, s, th(3), pass")
    o.append("      k = k0 + i")
    o.append("      p = mp_z(k)")
    o.append("      ea = mp_ea(k)")
    o.extend(e_setup("      "))
    o.extend(e1d_loop(la + 1, lb, "      "))
    o.append("      best = 0.0_dp")
    o.append(f"      do x2 = 1, {nb}")
    o.append(f"         do x1 = 1, {na}")
    o.extend(comp("            "))
    o.append("            do d = 1, 3")
    o.append("               ad = real(ia(d), dp)")
    o.append("               ia(d) = ia(d) + 1")
    o.append(f"               do hh = 1, {nh1}")
    o.append("                  " + EE)
    o.append("                  g(hh) = 2.0_dp*ea*ee")
    o.append("               end do")
    o.append("               ia(d) = ia(d) - 2")
    o.append("               if (ia(d) >= 0) then")
    o.append(f"                  do hh = 1, {nh1}")
    o.append("                     " + EE)
    o.append("                     g(hh) = g(hh) - ad*ee")
    o.append("                  end do")
    o.append("               end if")
    o.append("               ia(d) = ia(d) + 1")
    o.append(f"               do hh = 1, {nh1}")
    o.append("                  th(1) = HT(hh); th(2) = HU(hh); th(3) = HV(hh)")
    o.append("                  th(d) = th(d) - 1")
    o.append("                  gs(hh) = 0.0_dp")
    o.append("                  if (th(d) >= 0) then")
    o.append("                     gs(hh) = e(ia(1), ib(1), th(1), 1)*e(ia(2), ib(2), th(2), 2)*e(ia(3), ib(3), th(3), 3)")
    o.append("                  end if")
    o.append("               end do")
    o.append("               do pass = 1, 3")
    o.append("                  ! pass 1: d/dA. pass 2: d/dB = (d/dA + d/dB) - d/dA. pass 3: d/dA + d/dB, what an atom")
    o.append("                  ! carrying both shells of the pair accumulates (gA + (gT - gA) = gT).")
    o.append("                  if (pass == 2) then")
    o.append(f"                     do hh = 1, {nh1}")
    o.append("                        g(hh) = gs(hh) - g(hh)")
    o.append("                     end do")
    o.append("                  else if (pass == 3) then")
    o.append(f"                     do hh = 1, {nh1}")
    o.append("                        g(hh) = gs(hh)")
    o.append("                     end do")
    o.append("                  end if")
    o.append(f"                  do s = 0, {lp1}")
    o.append("                     acc(s) = 0.0_dp")
    o.append("                  end do")
    o.append(f"                  do q = 1, {len(pi)}")
    o.append("                     acc(PS(q)) = acc(PS(q)) + PW(q)*g(PI(q))*g(PJ(q))")
    o.append("                  end do")
    o.append("                  nq2 = acc(0)")
    o.append("                  ws = 1.0_dp")
    o.append(f"                  do s = 1, {lp1}")
    o.append("                     ws = -p*ws")
    o.append("                     nq2 = nq2 + ws*acc(s)*RO(s + 1)")
    o.append("                  end do")
    o.append("                  nq2 = nq2*TWO_PI_2_5/(p*p*sqrt(2.0_dp*p))")
    o.append("                  best = max(best, nq2)")
    o.append("               end do")
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    o.append("      mp_gip(k) = abs(mp_c(k))*sqrt(best)")
    o.append(f"   end subroutine mdjg_bound_item_{tag}")
    o.append("")
    o.append(f"end module trc_mdjg_e{tag}")
    return "\n".join(o) + "\n"


def gen_dispatch(lmax, header):
    lpmax = 2 * lmax
    o = [header.format(what="Host dispatch of the MD J gradient kernels")]
    o.append("module trc_mdjg_kernels")
    o.append("   use pic_types, only: dp")
    for la in range(lmax + 1):
        for lb in range(la + 1):
            o.append(f"   use trc_mdjg_e{la}{lb}, only: mdjg_bra_{la}{lb}, mdjg_bound_{la}{lb}")
    for lpb in range(1, lpmax + 2):
        for lq in range(lpmax + 1):
            o.append(f"   use trc_mdjg_h{lpb}{lq}, only: mdjg_h_{lpb}{lq}")
    o.append("   implicit none")
    o.append("   private")
    o.append("   public :: mdjg_bra, mdjg_bound, mdjg_hermite, MDJG_LMAX, MDJG_LPMAX")
    o.append("")
    o.append(f"   integer, parameter :: MDJG_LMAX = {lmax}")
    o.append("   !! Highest shell l the gradient kernels cover.")
    o.append(f"   integer, parameter :: MDJG_LPMAX = {lpmax}")
    o.append("   !! Largest Hermite degree of one shell pair they cover (the raised bra order is one more).")
    o.append("")
    o.append("contains")
    o.append("")

    def wrap(call, width=100, ind="         "):
        words = call.split(", ")
        lines, cur = [], ind
        for k, w in enumerate(words):
            tok = w + (", " if k < len(words) - 1 else "")
            if len(cur) + len(tok) > width:
                lines.append(cur.rstrip() + " &")
                cur = ind + "   "
            cur += tok
        lines.append(cur)
        return lines

    o.extend(G.wrap_args("subroutine mdjg_bra", "la, lb, " + BRA_ARGS, "   "))
    o.append("      integer,  intent(in) :: la, lb")
    o.append(BRA_DECL)
    o.append(f"      select case (la*{lmax + 1} + lb)")
    for la in range(lmax + 1):
        for lb in range(la + 1):
            o.append(f"      case ({la * (lmax + 1) + lb})")
            o.extend(wrap(f"call mdjg_bra_{la}{lb}({BRA_ARGS})"))
    o.append("      case default")
    o.append('         error stop "trc_mdjg_kernels: no gradient bra kernel for this shell-pair class"')
    o.append("      end select")
    o.append("   end subroutine mdjg_bra")
    o.append("")
    o.extend(G.wrap_args("subroutine mdjg_bound", "la, lb, " + BND_ARGS, "   "))
    o.append("      integer,  intent(in) :: la, lb")
    o.append(BND_DECL)
    o.append(f"      select case (la*{lmax + 1} + lb)")
    for la in range(lmax + 1):
        for lb in range(la + 1):
            o.append(f"      case ({la * (lmax + 1) + lb})")
            o.extend(wrap(f"call mdjg_bound_{la}{lb}({BND_ARGS})"))
    o.append("      case default")
    o.append('         error stop "trc_mdjg_kernels: no gradient bound kernel for this shell-pair class"')
    o.append("      end select")
    o.append("   end subroutine mdjg_bound")
    o.append("")
    args, idecl, rdecl = G.gather_signature(True)
    o.extend(G.wrap_args("subroutine mdjg_hermite", "lpb, lq, " + args, "   "))
    o.append(f"      !! Jt'[P] += [P|Q] Dt[Q] with the bra's Hermite order raised to lpb = lp + 1.")
    o.append(f"      integer,  intent(in) :: lpb, lq, {idecl}")
    o.append("      real(dp), intent(in) :: mp_z(nhp), mp_px(nhp), mp_py(nhp), mp_pz(nhp)")
    o.append(f"      real(dp), intent(in) :: mp_ip(nhp), mp_nq(nhp), {rdecl.split(', ')[0]}, cut")
    o.append(f"      real(dp), intent(inout) :: {rdecl.split(', ')[1]}")
    n = lpmax + 1
    o.append(f"      select case (lpb*{n} + lq)")
    for lpb in range(1, lpmax + 2):
        for lq in range(lpmax + 1):
            o.append(f"      case ({lpb * n + lq})")
            o.extend(wrap(f"call mdjg_h_{lpb}{lq}({args})"))
    o.append("      case default")
    o.append('         error stop "trc_mdjg_kernels: no raised-bra Hermite kernel for this class"')
    o.append("      end select")
    o.append("   end subroutine mdjg_hermite")
    o.append("")
    o.append("end module trc_mdjg_kernels")
    return "\n".join(o) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--lmax", type=int, required=True,
                    help="highest shell l the gradient covers (at most 3)")
    ap.add_argument("--order", default="lexicographic", choices=sorted(G.ORDER_CODE))
    ap.add_argument("--batches", default="", help="comma list L:n overriding the batch count per total L")
    ap.add_argument("--no-sieve", action="store_true")
    ap.add_argument("-o", "--outdir", required=True)
    a = ap.parse_args()
    if a.lmax > 3:
        raise SystemExit("gen_md_j_grad: the gradient is validated through f (--lmax 3); g would need "
                         "the split-module Hermite stage for the raised class [9|8]")
    lpmax = 2 * a.lmax
    if lpmax + 1 + lpmax > BOYS_MMAX:
        raise SystemExit(f"gen_md_j_grad: Hermite order {2 * lpmax + 1} exceeds the Boys table (BOYS_MMAX {BOYS_MMAX})")
    table = {}
    for item in filter(None, a.batches.split(",")):
        L, n = item.split(":")
        table[int(L)] = int(n)
    cmd = f"scripts/gen_md_j_grad.py --lmax {a.lmax}"
    if a.order not in ("lexicographic", "libcint"):
        cmd += f" --order {a.order}"
    if a.batches:
        cmd += f" --batches {a.batches}"
    if a.no_sieve:
        cmd += " --no-sieve"
    code = G.ORDER_CODE[a.order]
    header = ("!\n! {what}.\n!\n"
              f"! GENERATED by {cmd} -- do not edit.\n"
              "! See scripts/gen_md_j_grad.py for the algorithm and md/reference_grad.py for its check.\n!\n")
    e_header = header + (f"#if TRC_CART_ORDER != {code}\n"
                         '#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."\n'
                         "#endif\n")
    os.makedirs(a.outdir, exist_ok=True)

    def write(name, text):
        with open(os.path.join(a.outdir, name), "w") as fh:
            fh.write(text)

    for la in range(a.lmax + 1):
        for lb in range(la + 1):
            write(f"trc_mdjg_e{la}{lb}.F90", gen_e_module(la, lb, a.order, e_header))
    for lpb in range(1, lpmax + 2):
        for lq in range(lpmax + 1):
            write(f"trc_mdjg_h{lpb}{lq}.F90",
                  G.gen_r_module_gather(lpb, lq, G.nbatch_for(lpb + lq, table), header,
                                        sieve=not a.no_sieve, grad=True))
    write("trc_mdjg_kernels.F90", gen_dispatch(a.lmax, header))
    return 0


if __name__ == "__main__":
    sys.exit(main())
