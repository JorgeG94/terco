#!/usr/bin/env python3
"""Generate McMurchie-Davidson Fock kernels that contract through batched GEMMs.

The host driver (src/eri/trc_mdg.F90) walks, per class (la lb|lc ld), the
shell pairs of each (type, primitive count) group sorted by Schwarz bound,
in staircase windows: a block of ni bra pairs against the prefix of the
ket list that can still pass the screen. For each window it calls the
routines generated here:

  mdg_bra_<tag>   Hb[ab, pbar, i] = c_P E^ab_p, pbar = (p, k_ab)
  mdg_ket_<tag>   Hk likewise, then
                  R kernel      one thread per primitive quartet; stores the
                                n_r values R^(0)_r only (zero when screened)
                  expansion     PQ[pbar, i, qbar, j] = (-1)^|q| R_{p+q}, from
                                a literal (p, q) -> r table
                  GEMM 1        batched over j: S_j = PQ_j Hk_j^T
                  GEMM 2        batched over i: Eri_i = Hb_i S(i, :), S read in
                                place with leading dimension npb*ni
                  digestion     one thread per (quartet, ab row): K, and J when
                                jfac /= 0, with the HGP kernels' weights

Only the R recurrence is written out. The E matrices, the expansion and the
digestion are loops over literal tables: unrolled, they were 98% of an
(ff|ff) module's 26k lines and 751 s of its compile, for memory-bound work.

Usage: gen_md_kgemm.py --lmax 3 --classes 2121,3333 [--order ...] -o DIR
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cartorder  # noqa: E402
from md import rtree  # noqa: E402
from md.hermite import herm, herm_index  # noqa: E402
from md.schedule import dfs_order  # noqa: E402
from gen_md_j import ORDER_CODE, TWO_PI_2_5, decl, emit_sum, expr_of, wrap_args  # noqa: E402
from gen_md_k import class_key  # noqa: E402

DATA = ("nhp, mp_z, mp_px, mp_py, mp_pz, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, "
        "nsp, sp_ao1, sp_ao2, sp_k0, sp_nk, sp_sh1, sp_sh2, sp_q, ngl, gl")
DATA_DECL = """      integer,  intent(in) :: nhp, nsp, ngl
      real(dp), intent(in) :: mp_z(nhp), mp_px(nhp), mp_py(nhp), mp_pz(nhp)
      real(dp), intent(in) :: mp_ax(nhp), mp_ay(nhp), mp_az(nhp), mp_bx(nhp), mp_by(nhp), mp_bz(nhp)
      real(dp), intent(in) :: mp_c(nhp), sp_q(nsp)
      integer,  intent(in) :: sp_ao1(nsp), sp_ao2(nsp), sp_k0(nsp), sp_nk(nsp), sp_sh1(nsp), sp_sh2(nsp)
      integer,  intent(in) :: gl(ngl)"""
WIN = "ib0, ni, ka, jb0, nj, kb, thresh, pcut, jfac, kfac, wdiag, nbas, dsh, nao, ndens, dmat, jmat"
WIN_DECL = """      integer,  intent(in) :: ib0, ni, ka, jb0, nj, kb, nbas, nao, ndens
      real(dp), intent(in) :: thresh, pcut, jfac, kfac, wdiag, dsh(nbas, nbas)
      real(dp), intent(in) :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao)"""


def int_table(name, vals):
    body = ", ".join(str(v) for v in vals)
    lines, cur = [], f"   integer, parameter :: {name}({len(vals)}) = ["
    for tok in body.split(", "):
        if len(cur) + len(tok) > 100:
            lines.append(cur + " &")
            cur = "      "
        cur += tok + ", "
    lines.append(cur.rstrip(", ") + "]")
    return lines


def real_table(name, vals):
    toks = [("-1.0_dp" if v < 0 else "1.0_dp") for v in vals]
    lines, cur = [], f"   real(dp), parameter :: {name}({len(vals)}) = ["
    for tok in toks:
        if len(cur) + len(tok) > 100:
            lines.append(cur + " &")
            cur = "      "
        cur += tok + ", "
    lines.append(cur.rstrip(", ") + "]")
    return lines


def e1d_loop(l1, l2, ind):
    """E^{ij}_t, i <= l1, j <= l2, all three directions, into e(0:l1, 0:l2, 0:l1+l2, 3).

    The recursion of md/hermite.py without K_AB (the pair coefficient has it):
    i up first at j = 0, then j up at every i. `tt`, not `t`: the item
    routines' thread index is t, and ii, jj because i is the pair.
    """
    L = l1 + l2
    return [f"{ind}do d = 1, 3",
            f"{ind}   do tt = 0, {L}",
            f"{ind}      do jj = 0, {l2}",
            f"{ind}         do ii = 0, {l1}",
            f"{ind}            e(ii, jj, tt, d) = 0.0_dp",
            f"{ind}         end do",
            f"{ind}      end do",
            f"{ind}   end do",
            f"{ind}   e(0, 0, 0, d) = 1.0_dp",
            f"{ind}   do ii = 0, {l1 - 1}",
            f"{ind}      do tt = 0, ii + 1",
            f"{ind}         v = xa(d)*e(ii, 0, tt, d)",
            f"{ind}         if (tt > 0) v = v + oo2p*e(ii, 0, tt - 1, d)",
            f"{ind}         if (tt < {L}) v = v + real(tt + 1, dp)*e(ii, 0, tt + 1, d)",
            f"{ind}         e(ii + 1, 0, tt, d) = v",
            f"{ind}      end do",
            f"{ind}   end do",
            f"{ind}   do ii = 0, {l1}",
            f"{ind}      do jj = 0, {l2 - 1}",
            f"{ind}         do tt = 0, ii + jj + 1",
            f"{ind}            v = xb(d)*e(ii, jj, tt, d)",
            f"{ind}            if (tt > 0) v = v + oo2p*e(ii, jj, tt - 1, d)",
            f"{ind}            if (tt < {L}) v = v + real(tt + 1, dp)*e(ii, jj, tt + 1, d)",
            f"{ind}            e(ii, jj + 1, tt, d) = v",
            f"{ind}         end do",
            f"{ind}      end do",
            f"{ind}   end do",
            f"{ind}end do"]


def h_item(tag, side, l1, l2):
    """One primitive pair's E-matrix column block: h(ab, (k-1)*np + p, i)."""
    np_ = len(herm(l1 + l2))
    n1, n2 = (l1 + 1) * (l1 + 2) // 2, (l2 + 1) * (l2 + 2) // 2
    s = side.upper()
    nn, kk, b0 = ("ni", "ka", "ib0") if side == "b" else ("nj", "kb", "jb0")
    o = []
    o.extend(wrap_args(f"pure subroutine mdg_h{side}_item_{tag}", f"t, {DATA}, {b0}, {nn}, {kk}, h", "   "))
    o.append("      !$acc routine seq")
    o.append("      integer, intent(in) :: t")
    o.append(DATA_DECL)
    o.append(f"      integer,  intent(in) :: {b0}, {nn}, {kk}")
    o.append(f"      real(dp), intent(inout) :: h({n1 * n2}, {np_}*{kk}, {nn})")
    o.append("      integer :: i, k, sp, kp, col, x1, x2, ab, hh, d, ii, jj, tt")
    o.append("      real(dp) :: oo2p, cf, v, xa(3), xb(3)")
    o.append(f"      real(dp) :: e(0:{l1}, 0:{l2}, 0:{l1 + l2}, 3)")
    o.append("")
    o.append(f"      i = mod(t - 1, {nn}) + 1")
    o.append(f"      k = (t - 1)/{nn} + 1")
    o.append(f"      sp = gl({b0} + i)")
    o.append(f"      col = (k - 1)*{np_}")
    o.append("      if (k > sp_nk(sp)) then")
    o.append(f"         do hh = 1, {np_}")
    o.append(f"            do ab = 1, {n1 * n2}")
    o.append("               h(ab, col + hh, i) = 0.0_dp")
    o.append("            end do")
    o.append("         end do")
    o.append("         return")
    o.append("      end if")
    o.append("      kp = sp_k0(sp) + k")
    o.append("      oo2p = 0.5_dp/mp_z(kp)")
    o.append("      cf = mp_c(kp)")
    o.append("      xa(1) = mp_ax(kp); xa(2) = mp_ay(kp); xa(3) = mp_az(kp)")
    o.append("      xb(1) = mp_bx(kp); xb(2) = mp_by(kp); xb(3) = mp_bz(kp)")
    o.extend(e1d_loop(l1, l2, "      "))
    o.append(f"      do x2 = 1, {n2}")
    o.append(f"         do x1 = 1, {n1}")
    o.append(f"            ab = x1 + {n1}*(x2 - 1)")
    o.append(f"            do hh = 1, {np_}")
    o.append(f"               if (HT{np_}(hh) > C{s}1X(x1) + C{s}2X(x2) .or. HU{np_}(hh) > C{s}1Y(x1) + C{s}2Y(x2) .or. &")
    o.append(f"                   HV{np_}(hh) > C{s}1Z(x1) + C{s}2Z(x2)) then")
    o.append("                  h(ab, col + hh, i) = 0.0_dp")
    o.append("               else")
    o.append(f"                  h(ab, col + hh, i) = cf*e(C{s}1X(x1), C{s}2X(x2), HT{np_}(hh), 1) &")
    o.append(f"                                     *e(C{s}1Y(x1), C{s}2Y(x2), HU{np_}(hh), 2) &")
    o.append(f"                                     *e(C{s}1Z(x1), C{s}2Z(x2), HV{np_}(hh), 3)")
    o.append("               end if")
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    o.append(f"   end subroutine mdg_h{side}_item_{tag}")
    return o


def r_item(tag, L):
    g = rtree.RDag(herm(L))
    g.sieve()
    leaves = g.leaves()
    il = herm_index(L)
    o = []
    o.extend(wrap_args(f"pure subroutine mdg_r_item_{tag}", f"t, {DATA}, {WIN}, rb", "   "))
    o.append("      !$acc routine seq")
    o.append("      integer, intent(in) :: t")
    o.append(DATA_DECL)
    o.append(WIN_DECL)
    o.append(f"      real(dp), intent(inout) :: rb({len(herm(L))}, ni*ka*nj*kb)")
    o.append("      integer :: u, i, j, k1, k2, spa, spc, si, sj, sk, sl, kp, kq, r")
    o.append("      real(dp) :: zeta, eta, zpe, alpha, pref, tval, m2a, ww, bnd, pqx, pqy, pqz, qcut")
    o.append("      real(dp) :: f(0:BOYS_MMAX)")
    o.append(decl("real(dp)", [rtree.leaf_name(m) for m in leaves]))
    o.append(decl("real(dp)", [st.lhs for st in g.statements()]))
    o.append("")
    o.append("      ! Primitive quartet t = (i, k1, j, k2), bra pair fastest.")
    o.append("      u = t - 1")
    o.append("      i = mod(u, ni) + 1; u = u/ni")
    o.append("      k1 = mod(u, ka) + 1; u = u/ka")
    o.append("      j = mod(u, nj) + 1")
    o.append("      k2 = u/nj + 1")
    o.append("      spa = gl(ib0 + i); spc = gl(jb0 + j)")
    o.append("      si = sp_sh1(spa); sj = sp_sh2(spa); sk = sp_sh1(spc); sl = sp_sh2(spc)")
    o.append("      qcut = sp_q(spa)*sp_q(spc)")
    o.append("      ! The HGP kernels' quartet and primitive screens; what fails is a zero block.")
    o.append("      if (k1 > sp_nk(spa) .or. k2 > sp_nk(spc)) then")
    o.append(f"         do r = 1, {len(herm(L))}")
    o.append("            rb(r, t) = 0.0_dp")
    o.append("         end do")
    o.append("         return")
    o.append("      end if")
    o.append("      kp = sp_k0(spa) + k1")
    o.append("      kq = sp_k0(spc) + k2")
    o.append("      zeta = mp_z(kp); eta = mp_z(kq); zpe = zeta + eta")
    o.append("      pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))")
    o.append("      bnd = TWO_PI_2_5*abs(mp_c(kp))/(zeta*sqrt(zeta))")
    o.append("      if (qcut <= thresh .or. qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), dsh(si, sk), &")
    o.append("          dsh(si, sl), dsh(sj, sk), dsh(sj, sl)) <= thresh .or. &")
    o.append("          bnd*abs(mp_c(kq))/eta <= pcut .or. abs(pref*mp_c(kp)*mp_c(kq)) <= pcut) then")
    o.append(f"         do r = 1, {len(herm(L))}")
    o.append("            rb(r, t) = 0.0_dp")
    o.append("         end do")
    o.append("         return")
    o.append("      end if")
    o.append("      alpha = zeta*eta/zpe")
    o.append("      pqx = mp_px(kp) - mp_px(kq); pqy = mp_py(kp) - mp_py(kq); pqz = mp_pz(kp) - mp_pz(kq)")
    o.append("      tval = alpha*(pqx*pqx + pqy*pqy + pqz*pqz)")
    o.append(f"      call boys_down({max(leaves)}, tval, f)")
    o.append("      m2a = -2.0_dp*alpha")
    o.append("      ww = pref")
    for m in range(max(leaves) + 1):
        if m > 0:
            o.append("      ww = ww*m2a")
        if m in leaves:
            o.append(f"      {rtree.leaf_name(m)} = ww*f({m})")
    if ((0, 0, 0), 0) in g.tset:
        o.append(f"      rb(1, t) = {rtree.leaf_name(0)}")
    # Depth-first, each target stored as soon as it is formed.
    for st in dfs_order(g, sorted(g.targets, key=lambda n: (sum(n[0]), n[0]))):
        o.extend(emit_sum(st.lhs, expr_of(st.terms), " " * 6))
        fl = st.lhs[1:].split("_")
        if fl[3] == "0":
            r = tuple(int(x) for x in fl[:3])
            o.append(f"      rb({il[r] + 1}, t) = {st.lhs}")
    o.append(f"   end subroutine mdg_r_item_{tag}")
    return o


def x_item(tag, lp, lq, L):
    np_, nq_, nr = len(herm(lp)), len(herm(lq)), len(herm(L))
    o = []
    o.extend(wrap_args(f"pure subroutine mdg_x_item_{tag}", "t, ni, ka, nj, kb, rb, pq", "   "))
    o.append("      !! PQ[pbar, i, qbar, j] = (-1)^|q| R_{p+q}: one thread per (pbar, i, k2, j).")
    o.append("      !$acc routine seq")
    o.append("      integer,  intent(in) :: t, ni, ka, nj, kb")
    o.append(f"      real(dp), intent(in) :: rb({nr}, ni, ka, nj, kb)")
    o.append(f"      real(dp), intent(inout) :: pq({np_}, ka, ni, {nq_}, kb, nj)")
    o.append("      integer :: u, x, k1, i, k2, j, y")
    o.append("      u = t - 1")
    o.append(f"      x = mod(u, {np_}) + 1; u = u/{np_}")
    o.append("      k1 = mod(u, ka) + 1; u = u/ka")
    o.append("      i = mod(u, ni) + 1; u = u/ni")
    o.append("      k2 = mod(u, kb) + 1")
    o.append("      j = u/kb + 1")
    o.append(f"      do y = 1, {nq_}")
    o.append(f"         pq(x, k1, i, y, k2, j) = SQ{tag}(y)*rb(RM{tag}(x + {np_}*(y - 1)), i, k1, j, k2)")
    o.append("      end do")
    o.append(f"   end subroutine mdg_x_item_{tag}")
    return o


def d_item(tag, na, nb, nc, nd):
    nab, ncd = na * nb, nc * nd
    o = []
    o.extend(wrap_args(f"pure subroutine mdg_d_item_{tag}", f"t, {DATA}, {WIN}, eri", "   "))
    o.append("      !$acc routine seq")
    o.append("      integer, intent(in) :: t")
    o.append(DATA_DECL)
    o.append(WIN_DECL)
    o.append(f"      real(dp), intent(in) :: eri({nab}, ni, {ncd}, nj)")
    o.append("      integer :: u, i, j, ab, ia, ib, ic, id, spa, spc, si, sj, sk, sl, mui, nuj, lamk, sigl, idens")
    o.append("      real(dp) :: qcut, w, sm, x")
    o.append(f"      real(dp) :: row({ncd})")
    o.append("")
    o.append("      ! Row ab of quartet (i, j); ab fastest, so a warp reads contiguous rows.")
    o.append(f"      ab = mod(t - 1, {nab}) + 1")
    o.append(f"      u = (t - 1)/{nab}")
    o.append("      i = mod(u, ni) + 1")
    o.append("      j = u/ni + 1")
    o.append("      spa = gl(ib0 + i); spc = gl(jb0 + j)")
    o.append("      si = sp_sh1(spa); sj = sp_sh2(spa); sk = sp_sh1(spc); sl = sp_sh2(spc)")
    o.append("      qcut = sp_q(spa)*sp_q(spc)")
    o.append("      if (qcut <= thresh) return")
    o.append("      if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), dsh(si, sk), dsh(si, sl), &")
    o.append("                   dsh(sj, sk), dsh(sj, sl)) <= thresh) return")
    o.append("      ! wdiag is 1/2 on a group against itself, done as a full square.")
    o.append("      w = wdiag")
    o.append("      if (si == sj) w = 0.5_dp*w")
    o.append("      if (sk == sl) w = 0.5_dp*w")
    o.append("      mui = sp_ao1(spa); nuj = sp_ao2(spa); lamk = sp_ao1(spc); sigl = sp_ao2(spc)")
    o.append(f"      ia = mod(ab - 1, {na})")
    o.append(f"      ib = (ab - 1)/{na}")
    o.append(f"      do ic = 1, {ncd}")
    o.append("         row(ic) = eri(ab, i, ic, j)")
    o.append("      end do")
    o.append("      do idens = 1, ndens")
    o.append("         if (kfac /= 0.0_dp) then")
    for (outer, ol, inner, il, lhs_a, lhs_b, dm) in [
            ("ic", nc, "id", nd, "mui + ia", "lamk + ic", "nuj + ib, sigl + id"),
            ("ic", nc, "id", nd, "nuj + ib", "lamk + ic", "mui + ia, sigl + id"),
            ("id", nd, "ic", nc, "mui + ia", "sigl + id", "nuj + ib, lamk + ic"),
            ("id", nd, "ic", nc, "nuj + ib", "sigl + id", "mui + ia, lamk + ic")]:
        o.append(f"            do {outer} = 0, {ol - 1}")
        o.append("               sm = 0.0_dp")
        o.append(f"               do {inner} = 0, {il - 1}")
        o.append(f"                  sm = sm + row(1 + ic + {nc}*id)*dmat(idens, {dm})")
        o.append("               end do")
        o.append("               !$acc atomic update")
        o.append(f"               jmat(idens, {lhs_a}, {lhs_b}) = jmat(idens, {lhs_a}, {lhs_b}) - kfac*w*sm")
        o.append("            end do")
    o.append("         end if")
    o.append("         if (jfac /= 0.0_dp) then")
    o.append("            sm = 0.0_dp")
    o.append(f"            do id = 0, {nd - 1}")
    o.append(f"               do ic = 0, {nc - 1}")
    o.append(f"                  sm = sm + row(1 + ic + {nc}*id)*dmat(idens, lamk + ic, sigl + id)")
    o.append("               end do")
    o.append("            end do")
    o.append("            !$acc atomic update")
    o.append("            jmat(idens, mui + ia, nuj + ib) = jmat(idens, mui + ia, nuj + ib) + 4.0_dp*jfac*w*sm")
    o.append("            x = 4.0_dp*jfac*w*dmat(idens, mui + ia, nuj + ib)")
    o.append(f"            do id = 0, {nd - 1}")
    o.append(f"               do ic = 0, {nc - 1}")
    o.append("                  !$acc atomic update")
    o.append("                  jmat(idens, lamk + ic, sigl + id) = jmat(idens, lamk + ic, sigl + id) &")
    o.append(f"                                                    + x*row(1 + ic + {nc}*id)")
    o.append("               end do")
    o.append("            end do")
    o.append("         end if")
    o.append("      end do")
    o.append(f"   end subroutine mdg_d_item_{tag}")
    return o


def gen_class(la, lb, lc, ld, order, header):
    tag = f"{la}{lb}{lc}{ld}"
    lp, lq = la + lb, lc + ld
    L = lp + lq
    hp, hq = herm(lp), herm(lq)
    np_, nq_, nr = len(hp), len(hq), len(herm(L))
    ca, cb = cartorder.components(la, order), cartorder.components(lb, order)
    cc, cd = cartorder.components(lc, order), cartorder.components(ld, order)
    na, nb, nc, nd = len(ca), len(cb), len(cc), len(cd)
    nab, ncd = na * nb, nc * nd
    il = herm_index(L)
    o = [header.format(what=f"Class ({la}{lb}|{lc}{ld}) MD Fock kernels through batched GEMMs")]
    o.append(f"module trc_mdg_k{tag}")
    o.append("   use pic_types, only: dp, int64")
    o.append("   use trc_boys, only: BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV")
    o.append("   use trc_linalg, only: trc_linalg_t")
    o.append("   implicit none")
    o.append("   private")
    o.append(f"   public :: mdg_bra_{tag}, mdg_ket_{tag}, MDG_SIZES_{tag}")
    o.append(f"   real(dp), parameter :: TWO_PI_2_5 = {TWO_PI_2_5}")
    o.append(f"   integer, parameter :: MDG_SIZES_{tag}(6) = [{np_}, {nq_}, {nr}, {nab}, {ncd}, {L}]")
    o.append("   !! np, nq, n_r, bra and ket component pairs, total L: the host sizes buffers from these.")
    o.append("   ! Literal tables the loops index: Hermite decode, component exponents, and")
    o.append("   ! the (p, q) -> p + q map and (-1)^|q| sign of the expansion.")
    for lab, hs in sorted({len(hp): hp, len(hq): hq}.items()):
        o.extend(int_table(f"HT{lab}", [h[0] for h in hs]))
        o.extend(int_table(f"HU{lab}", [h[1] for h in hs]))
        o.extend(int_table(f"HV{lab}", [h[2] for h in hs]))
    for nm, cs in (("CB1", ca), ("CB2", cb), ("CK1", cc), ("CK2", cd)):
        for d, ax in enumerate("XYZ"):
            o.extend(int_table(f"{nm}{ax}", [c[d] for c in cs]))
    o.extend(int_table(f"RM{tag}", [il[tuple(p[d] + q[d] for d in range(3))] + 1 for q in hq for p in hp]))
    o.extend(real_table(f"SQ{tag}", [(-1) ** sum(q) for q in hq]))
    o.append("")
    o.append("contains")
    o.append("")
    o.append('#include "inc/trc_boys_down.inc"')
    o.append("")
    # host: bra block
    o.extend(wrap_args(f"subroutine mdg_bra_{tag}", f"{DATA}, ib0, ni, ka, hb", "   "))
    o.append("      !! The bra block's E matrices, once per block of the staircase.")
    o.append(DATA_DECL)
    o.append("      integer,  intent(in) :: ib0, ni, ka")
    o.append(f"      real(dp), intent(inout) :: hb({nab}, {np_}*ka, ni)")
    o.append("      integer :: t")
    o.append("      do concurrent(t = 1:ni*ka)")
    o.extend(wrap_args(f"call mdg_hb_item_{tag}", f"t, {DATA}, ib0, ni, ka, hb", " " * 9))
    o.append("      end do")
    o.append(f"   end subroutine mdg_bra_{tag}")
    o.append("")
    o.extend(wrap_args(f"subroutine mdg_ket_{tag}", f"lp, {DATA}, {WIN}, hb, hk, rb, pq, sm, eri", "   "))
    o.append("      !! One window: ket E matrices, R, PQ, the two GEMMs and the digestion.")
    o.append("      type(trc_linalg_t), intent(inout) :: lp")
    o.append(DATA_DECL)
    o.append(WIN_DECL)
    o.append(f"      real(dp), intent(in)    :: hb({nab}, {np_}*ka, ni)")
    o.append(f"      real(dp), intent(inout) :: hk({ncd}, {nq_}*kb, nj), rb({nr}, ni*ka*nj*kb)")
    o.append(f"      real(dp), intent(inout) :: pq({np_}*ka*ni, {nq_}*kb, nj), sm({np_}*ka*ni, {ncd}, nj)")
    o.append(f"      real(dp), intent(inout) :: eri({nab}, ni, {ncd}, nj)")
    o.append("      integer :: t, npb, nqb")
    o.append(f"      npb = {np_}*ka")
    o.append(f"      nqb = {nq_}*kb")
    o.append("      do concurrent(t = 1:nj*kb)")
    o.extend(wrap_args(f"call mdg_hk_item_{tag}", f"t, {DATA}, jb0, nj, kb, hk", " " * 9))
    o.append("      end do")
    o.append("      do concurrent(t = 1:ni*ka*nj*kb)")
    o.extend(wrap_args(f"call mdg_r_item_{tag}", f"t, {DATA}, {WIN}, rb", " " * 9))
    o.append("      end do")
    o.append(f"      do concurrent(t = 1:{np_}*ka*ni*kb*nj)")
    o.append(f"         call mdg_x_item_{tag}(t, ni, ka, nj, kb, rb, pq)")
    o.append("      end do")
    o.append("      ! S_j = PQ_j Hk_j^T, then Eri_i = Hb_i S(i, :) with S read in place.")
    o.append(f"      call lp%gemm_strided(\"N\", \"T\", npb*ni, {ncd}, nqb, 1.0_dp, pq, npb*ni, &")
    o.append(f"                           int(npb, int64)*ni*nqb, hk, {ncd}, int({ncd}*nqb, int64), 0.0_dp, &")
    o.append(f"                           sm, npb*ni, int(npb, int64)*ni*{ncd}, nj)")
    o.append(f"      call lp%gemm_strided(\"N\", \"N\", {nab}, {ncd}*nj, npb, 1.0_dp, hb, {nab}, &")
    o.append(f"                           int({nab}*npb, int64), sm, npb*ni, int(npb, int64), 0.0_dp, eri, &")
    o.append(f"                           {nab}*ni, int({nab}, int64), ni)")
    o.append(f"      do concurrent(t = 1:{nab}*ni*nj)")
    o.extend(wrap_args(f"call mdg_d_item_{tag}", f"t, {DATA}, {WIN}, eri", " " * 9))
    o.append("      end do")
    o.append(f"   end subroutine mdg_ket_{tag}")
    o.append("")
    o.extend(h_item(tag, "b", la, lb))
    o.append("")
    o.extend(h_item(tag, "k", lc, ld))
    o.append("")
    o.extend(r_item(tag, L))
    o.append("")
    o.extend(x_item(tag, lp, lq, L))
    o.append("")
    o.extend(d_item(tag, na, nb, nc, nd))
    o.append("")
    o.append(f"end module trc_mdg_k{tag}")
    return "\n".join(o) + "\n"


def gen_dispatch(classes, header):
    o = [header.format(what="Host dispatch of the MD GEMM Fock kernels")]
    o.append("module trc_mdg_kernels")
    o.append("   use pic_types, only: dp")
    o.append("   use trc_linalg, only: trc_linalg_t")
    for c in classes:
        tag = "".join(map(str, c))
        o.append(f"   use trc_mdg_k{tag}, only: mdg_bra_{tag}, mdg_ket_{tag}, MDG_SIZES_{tag}")
    o.append("   implicit none")
    o.append("   private")
    o.append("   public :: mdg_class_count, mdg_class_key, mdg_sizes, mdg_bra, mdg_ket")
    o.append("")
    o.append(f"   integer, parameter :: NCLASS = {len(classes)}")
    if classes:
        o.append(f"   integer, parameter :: KEYS(NCLASS) = [{', '.join(str(class_key(*c)) for c in classes)}]")
    else:
        o.append("   integer, parameter :: KEYS(1) = [-1]")
    o.append("")
    o.append("contains")
    o.append("")
    o.append("   pure integer function mdg_class_count() result(n)")
    o.append("      n = NCLASS")
    o.append("   end function mdg_class_count")
    o.append("")
    o.append("   pure integer function mdg_class_key(k) result(key)")
    o.append("      integer, intent(in) :: k")
    o.append("      key = KEYS(k)")
    o.append("   end function mdg_class_key")
    o.append("")
    o.append("   pure function mdg_sizes(key) result(sz)")
    o.append("      !! np, nq, n_r, nab, ncd, L of a generated class; zeros if none.")
    o.append("      integer, intent(in) :: key")
    o.append("      integer :: sz(6)")
    o.append("      sz = 0")
    for c in classes:
        tag = "".join(map(str, c))
        o.append(f"      if (key == {class_key(*c)}) sz = MDG_SIZES_{tag}")
    o.append("   end function mdg_sizes")
    o.append("")
    o.extend(wrap_args("subroutine mdg_bra", f"key, {DATA}, ib0, ni, ka, hb", "   "))
    o.append("      integer, intent(in) :: key")
    o.append(DATA_DECL)
    o.append("      integer,  intent(in) :: ib0, ni, ka")
    o.append("      real(dp), intent(inout) :: hb(*)")
    o.append("      select case (key)")
    for c in classes:
        tag = "".join(map(str, c))
        o.append(f"      case ({class_key(*c)})")
        o.extend(wrap_args(f"call mdg_bra_{tag}", f"{DATA}, ib0, ni, ka, hb", " " * 9))
    o.append("      case default")
    o.append('         error stop "trc_mdg_kernels: no MD GEMM kernel for this class"')
    o.append("      end select")
    o.append("   end subroutine mdg_bra")
    o.append("")
    o.extend(wrap_args("subroutine mdg_ket", f"key, lp, {DATA}, {WIN}, hb, hk, rb, pq, sm, eri", "   "))
    o.append("      integer, intent(in) :: key")
    o.append("      type(trc_linalg_t), intent(inout) :: lp")
    o.append(DATA_DECL)
    o.append(WIN_DECL)
    o.append("      real(dp), intent(in)    :: hb(*)")
    o.append("      real(dp), intent(inout) :: hk(*), rb(*), pq(*), sm(*), eri(*)")
    o.append("      select case (key)")
    for c in classes:
        tag = "".join(map(str, c))
        o.append(f"      case ({class_key(*c)})")
        o.extend(wrap_args(f"call mdg_ket_{tag}", f"lp, {DATA}, {WIN}, hb, hk, rb, pq, sm, eri", " " * 9))
    o.append("      case default")
    o.append('         error stop "trc_mdg_kernels: no MD GEMM kernel for this class"')
    o.append("      end select")
    o.append("   end subroutine mdg_ket")
    o.append("")
    o.append("end module trc_mdg_kernels")
    return "\n".join(o) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--classes", required=True)
    ap.add_argument("--lmax", type=int, required=True)
    ap.add_argument("--order", default="lexicographic", choices=sorted(ORDER_CODE))
    ap.add_argument("-o", "--outdir", required=True)
    a = ap.parse_args()
    classes = [tuple(int(ch) for ch in c) for c in filter(None, a.classes.split(","))]
    classes = [c for c in classes if max(c) <= a.lmax]
    for c in classes:
        if c[0] < c[1] or c[2] < c[3] or (c[0], c[1]) < (c[2], c[3]):
            raise SystemExit(f"class {c}: la >= lb, lc >= ld and (la, lb) >= (lc, ld), as fock_bins forms them")
    cmd = f"scripts/gen_md_kgemm.py --lmax {a.lmax} --classes {a.classes}"
    if a.order not in ("lexicographic", "libcint"):
        cmd += f" --order {a.order}"
    header = ("!\n! {what}.\n!\n"
              f"! GENERATED by {cmd} -- do not edit.\n!\n"
              f"#if TRC_CART_ORDER != {ORDER_CODE[a.order]}\n"
              '#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."\n'
              "#endif\n")
    os.makedirs(a.outdir, exist_ok=True)
    for c in classes:
        with open(os.path.join(a.outdir, f"trc_mdg_k{''.join(map(str, c))}.F90"), "w") as fh:
            fh.write(gen_class(*c, a.order, header))
    with open(os.path.join(a.outdir, "trc_mdg_kernels.F90"), "w") as fh:
        fh.write(gen_dispatch(classes, header))
    return 0


if __name__ == "__main__":
    sys.exit(main())
