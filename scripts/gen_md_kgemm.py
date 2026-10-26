#!/usr/bin/env python3
"""Generate McMurchie-Davidson Fock kernels that contract through batched GEMMs.

For a segment of fock_bins (bra bin x ket bin: fixed la lb K_ab and lc ld
K_cd), one window of ni bra pairs and nj ket pairs at a time:

  H kernels   Hb[ab, pbar, i] = c_P E^ab_p, pbar = (p, k_ab), likewise Hk
  PQ kernel   one thread per primitive quartet: R, then
              PQ[pbar, i, qbar, j] = (-1)^|q| R_{p+q}, zero when screened
  GEMM 1      batched over j: S_j = PQ_j Hk_j^T        (npb*ni x ncd)
  GEMM 2      batched over i: E_i = Hb_i S(i, :)       (nab x ncd*nj),
              S read in place with leading dimension npb*ni
  digestion   one thread per (quartet, ab): the row into K, and J when
              jfac /= 0, with the HGP kernels' degeneracy weights

Diagonal segments are done as full squares, every quartet at half weight:
an off-diagonal quartet then comes twice, at 1/2 each, and a diagonal one
once, at the 1/2 the HGP kernels give it.

The GEMMs are trc_linalg's gemm_strided: cuBLAS DGEMM on the device (the
FP64 units on a V100; FP64 tensor cores where the math mode allows), a
pic_gemm loop on the host. Nothing else in the class is a quartet in
registers: R is per primitive quartet, every contraction is a matrix.

Usage: gen_md_kgemm.py --lmax 3 --classes 2121,3333 [--order ...] -o DIR
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cartorder  # noqa: E402
from md import rtree  # noqa: E402
from md.hermite import herm, herm_index, e1d_statements, e_operand  # noqa: E402
from gen_md_j import ORDER_CODE, TWO_PI_2_5, decl, emit_sum, expr_of, lit, wrap_args  # noqa: E402
from gen_md_k import class_key, decode_lines, e_names_for, e_stmts, e3  # noqa: E402

SEG = ("s, nseg, sA, sNB, sOA, sOB, sPA, sPB, sD, npair, sp_i, sp_j, sp_q, thresh, pcut, "
       "jfac, kfac, dsh, nbas, nao, ao_off, nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, "
       "sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, ndens, dmat, jmat")
SEG_DECL = """      integer,  intent(in)    :: s, nseg, npair, nbas, nao, nqp, ndens
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg), sPA(nseg), sPB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair), ao_off(nbas), sq_n(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh, pcut, jfac, kfac, dsh(nbas, nbas)
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp), qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(npair), sq_ay(npair), sq_az(npair)
      real(dp), intent(in)    :: sq_abx(npair), sq_aby(npair), sq_abz(npair)
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao)"""
WIN = "i0, ni, j0, nj, ka, kb"
WIN_DECL = "      integer,  intent(in)    :: i0, ni, j0, nj, ka, kb"


def h_item(tag, side, l1, l2, order):
    """E-matrix of one primitive pair of the bra ('b') or ket ('k') bin."""
    lp = l1 + l2
    hp = herm(lp)
    c1, c2 = cartorder.components(l1, order), cartorder.components(l2, order)
    n1 = len(c1)
    need = e_names_for(l1, l2, c1, c2)
    lines, names = e_stmts(l1, l2, need, " " * 6, "")
    base, cnt, stride = ("sPA(s)", "sA(s)", "sOA(s)") if side == "b" else ("sPB(s)", "sNB(s)", "sOB(s)")
    off = "i0" if side == "b" else "j0"
    nn = "ni" if side == "b" else "nj"
    kk = "ka" if side == "b" else "kb"
    o = []
    o.extend(wrap_args(f"pure subroutine mdg_h{side}_item_{tag}", "t, " + SEG + ", " + WIN + ", h", "   "))
    o.append("      !$acc routine seq")
    o.append("      integer, intent(in) :: t")
    o.append(SEG_DECL)
    o.append(WIN_DECL)
    o.append(f"      real(dp), intent(inout) :: h({n1 * len(c2)}, {len(hp)}*{kk}, {nn})")
    o.append("      integer :: i, k, slot, kp, col")
    o.append("      real(dp) :: p, oo2p, pxa, pya, pza, pxb, pyb, pzb, cf")
    o.append(decl("real(dp)", names))
    o.append("")
    o.append(f"      i = mod(t - 1, {nn}) + 1")
    o.append(f"      k = (t - 1)/{nn} + 1")
    o.append(f"      slot = {stride} + {off} + i")
    o.append(f"      col = (k - 1)*{len(hp)}")
    o.append("      if (k > sq_n(slot)) then")
    o.append(f"         h(:, col + 1:col + {len(hp)}, i) = 0.0_dp")
    o.append("         return")
    o.append("      end if")
    o.append(f"      kp = {base} + (k - 1)*{cnt} + {off} + i")
    o.append("      p = qp_z(kp)")
    o.append("      oo2p = 0.5_dp/p")
    o.append("      cf = qp_c(kp)")
    o.append("      pxa = qp_px(kp) - sq_ax(slot); pya = qp_py(kp) - sq_ay(slot); pza = qp_pz(kp) - sq_az(slot)")
    o.append("      pxb = pxa + sq_abx(slot); pyb = pya + sq_aby(slot); pzb = pza + sq_abz(slot)")
    o.extend(lines)
    for x1, a in enumerate(c1):
        for x2, b in enumerate(c2):
            ab = 1 + x1 + n1 * x2
            for hi, h in enumerate(hp):
                f = e3("", a, b, h)
                rhs = "0.0_dp" if f is None else "cf*" + "*".join(f) if f else "cf"
                o.append(f"      h({ab}, col + {hi + 1}, i) = {rhs}")
    o.append(f"   end subroutine mdg_h{side}_item_{tag}")
    return o


def pq_item(tag, lp, lq):
    L = lp + lq
    g = rtree.RDag(herm(L))
    g.sieve()
    leaves = g.leaves()
    hp, hq = herm(lp), herm(lq)
    np_, nq_ = len(hp), len(hq)
    o = []
    o.extend(wrap_args(f"pure subroutine mdg_pq_item_{tag}", "t, " + SEG + ", " + WIN + ", pq", "   "))
    o.append("      !$acc routine seq")
    o.append("      integer, intent(in) :: t")
    o.append(SEG_DECL)
    o.append(WIN_DECL)
    o.append(f"      real(dp), intent(inout) :: pq({np_}*ka, ni, {nq_}*kb, nj)")
    o.append("      integer :: u, i, j, k1, k2, slotab, slotcd, si, sj, sk, sl, kp, kq, cp, cq, x, y")
    o.append("      real(dp) :: qcut, zeta, eta, zpe, alpha, pref, tval, m2a, ww, bnd")
    o.append("      real(dp) :: pqx, pqy, pqz")
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
    o.append(f"      cp = (k1 - 1)*{np_}")
    o.append(f"      cq = (k2 - 1)*{nq_}")
    o.append("      slotab = sOA(s) + i0 + i; slotcd = sOB(s) + j0 + j")
    o.append("      ! Screened primitive quartets are zero blocks: the GEMMs multiply them anyway.")
    o.append("      if (k1 > sq_n(slotab) .or. k2 > sq_n(slotcd)) then")
    o.append(f"         do y = 1, {nq_}")
    o.append(f"            do x = 1, {np_}")
    o.append("               pq(cp + x, i, cq + y, j) = 0.0_dp")
    o.append("            end do")
    o.append("         end do")
    o.append("         return")
    o.append("      end if")
    o.append("      qcut = sp_q(slotab)*sp_q(slotcd)")
    o.append("      si = sp_i(slotab); sj = sp_j(slotab); sk = sp_i(slotcd); sl = sp_j(slotcd)")
    o.append("      kp = sPA(s) + (k1 - 1)*sA(s) + i0 + i")
    o.append("      kq = sPB(s) + (k2 - 1)*sNB(s) + j0 + j")
    o.append("      zeta = qp_z(kp); eta = qp_z(kq); zpe = zeta + eta")
    o.append("      pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))")
    o.append("      bnd = TWO_PI_2_5*abs(qp_c(kp))/(zeta*sqrt(zeta))")
    o.append("      ! The HGP kernels' quartet and primitive screens, so both drop the same work.")
    o.append("      if (qcut <= thresh .or. qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), dsh(si, sk), &")
    o.append("          dsh(si, sl), dsh(sj, sk), dsh(sj, sl)) <= thresh .or. &")
    o.append("          bnd*abs(qp_c(kq))/eta <= pcut .or. abs(pref*qp_c(kp)*qp_c(kq)) <= pcut) then")
    o.append(f"         do y = 1, {nq_}")
    o.append(f"            do x = 1, {np_}")
    o.append("               pq(cp + x, i, cq + y, j) = 0.0_dp")
    o.append("            end do")
    o.append("         end do")
    o.append("         return")
    o.append("      end if")
    o.append("      alpha = zeta*eta/zpe")
    o.append("      pqx = qp_px(kp) - qp_px(kq); pqy = qp_py(kp) - qp_py(kq); pqz = qp_pz(kp) - qp_pz(kq)")
    o.append("      tval = alpha*(pqx*pqx + pqy*pqy + pqz*pqz)")
    o.append(f"      call boys_down({max(leaves)}, tval, f)")
    o.append("      m2a = -2.0_dp*alpha")
    o.append("      ww = pref")
    for m in range(max(leaves) + 1):
        if m > 0:
            o.append("      ww = ww*m2a")
        if m in leaves:
            o.append(f"      {rtree.leaf_name(m)} = ww*f({m})")
    for st in g.statements():
        o.extend(emit_sum(st.lhs, expr_of(st.terms), " " * 6))
    for y, q in enumerate(hq):
        sg = "-" if sum(q) % 2 else ""
        for x, p in enumerate(hp):
            r = tuple(p[d] + q[d] for d in range(3))
            val = rtree.leaf_name(0) if r == (0, 0, 0) else rtree.r_name(r, 0)
            o.append(f"      pq(cp + {x + 1}, i, cq + {y + 1}, j) = {sg}{val}")
    o.append(f"   end subroutine mdg_pq_item_{tag}")
    return o


def dig_item(tag, na, nb, nc, nd):
    nab, ncd = na * nb, nc * nd
    o = []
    o.extend(wrap_args(f"pure subroutine mdg_dig_item_{tag}", "t, " + SEG + ", " + WIN + ", eri", "   "))
    o.append("      !$acc routine seq")
    o.append("      integer, intent(in) :: t")
    o.append(SEG_DECL)
    o.append(WIN_DECL)
    o.append(f"      real(dp), intent(in) :: eri({nab}, ni, {ncd}, nj)")
    o.append("      integer :: u, i, j, ab, ia, ib, ic, id, slotab, slotcd, si, sj, sk, sl")
    o.append("      integer :: mui, nuj, lamk, sigl, idens")
    o.append("      real(dp) :: qcut, w, sm, x")
    o.append("")
    o.append("      ! Row ab of quartet (i, j); ab fastest so a warp reads contiguous rows.")
    o.append(f"      ab = mod(t - 1, {nab}) + 1")
    o.append(f"      u = (t - 1)/{nab}")
    o.append("      i = mod(u, ni) + 1")
    o.append("      j = u/ni + 1")
    o.append("      slotab = sOA(s) + i0 + i; slotcd = sOB(s) + j0 + j")
    o.append("      qcut = sp_q(slotab)*sp_q(slotcd)")
    o.append("      si = sp_i(slotab); sj = sp_j(slotab); sk = sp_i(slotcd); sl = sp_j(slotcd)")
    o.append("      if (qcut <= thresh) return")
    o.append("      if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), dsh(si, sk), dsh(si, sl), &")
    o.append("                   dsh(sj, sk), dsh(sj, sl)) <= thresh) return")
    o.append("      w = 1.0_dp")
    o.append("      if (si == sj) w = 0.5_dp*w")
    o.append("      if (sk == sl) w = 0.5_dp*w")
    o.append("      ! Diagonal segments are full squares: every quartet at half weight.")
    o.append("      if (sD(s)) w = 0.5_dp*w")
    o.append("      mui = ao_off(si); nuj = ao_off(sj); lamk = ao_off(sk); sigl = ao_off(sl)")
    o.append(f"      ia = mod(ab - 1, {na})")
    o.append(f"      ib = (ab - 1)/{na}")
    o.append("      do idens = 1, ndens")
    o.append("         if (kfac /= 0.0_dp) then")
    o.append(f"            do ic = 0, {nc - 1}")
    o.append("               sm = 0.0_dp")
    o.append(f"               do id = 0, {nd - 1}")
    o.append(f"                  sm = sm + eri(ab, i, 1 + ic + {nc}*id, j)*dmat(idens, nuj + ib, sigl + id)")
    o.append("               end do")
    o.append("               !$acc atomic update")
    o.append("               jmat(idens, mui + ia, lamk + ic) = jmat(idens, mui + ia, lamk + ic) - kfac*w*sm")
    o.append("               sm = 0.0_dp")
    o.append(f"               do id = 0, {nd - 1}")
    o.append(f"                  sm = sm + eri(ab, i, 1 + ic + {nc}*id, j)*dmat(idens, mui + ia, sigl + id)")
    o.append("               end do")
    o.append("               !$acc atomic update")
    o.append("               jmat(idens, nuj + ib, lamk + ic) = jmat(idens, nuj + ib, lamk + ic) - kfac*w*sm")
    o.append("            end do")
    o.append(f"            do id = 0, {nd - 1}")
    o.append("               sm = 0.0_dp")
    o.append(f"               do ic = 0, {nc - 1}")
    o.append(f"                  sm = sm + eri(ab, i, 1 + ic + {nc}*id, j)*dmat(idens, nuj + ib, lamk + ic)")
    o.append("               end do")
    o.append("               !$acc atomic update")
    o.append("               jmat(idens, mui + ia, sigl + id) = jmat(idens, mui + ia, sigl + id) - kfac*w*sm")
    o.append("               sm = 0.0_dp")
    o.append(f"               do ic = 0, {nc - 1}")
    o.append(f"                  sm = sm + eri(ab, i, 1 + ic + {nc}*id, j)*dmat(idens, mui + ia, lamk + ic)")
    o.append("               end do")
    o.append("               !$acc atomic update")
    o.append("               jmat(idens, nuj + ib, sigl + id) = jmat(idens, nuj + ib, sigl + id) - kfac*w*sm")
    o.append("            end do")
    o.append("         end if")
    o.append("         if (jfac /= 0.0_dp) then")
    o.append("            sm = 0.0_dp")
    o.append(f"            do id = 0, {nd - 1}")
    o.append(f"               do ic = 0, {nc - 1}")
    o.append(f"                  sm = sm + eri(ab, i, 1 + ic + {nc}*id, j)*dmat(idens, lamk + ic, sigl + id)")
    o.append("               end do")
    o.append("            end do")
    o.append("            !$acc atomic update")
    o.append("            jmat(idens, mui + ia, nuj + ib) = jmat(idens, mui + ia, nuj + ib) + 4.0_dp*jfac*w*sm")
    o.append("            x = 4.0_dp*jfac*w*dmat(idens, mui + ia, nuj + ib)")
    o.append(f"            do id = 0, {nd - 1}")
    o.append(f"               do ic = 0, {nc - 1}")
    o.append("                  !$acc atomic update")
    o.append(f"                  jmat(idens, lamk + ic, sigl + id) = jmat(idens, lamk + ic, sigl + id) &")
    o.append(f"                                                    + x*eri(ab, i, 1 + ic + {nc}*id, j)")
    o.append("               end do")
    o.append("            end do")
    o.append("         end if")
    o.append("      end do")
    o.append(f"   end subroutine mdg_dig_item_{tag}")
    return o


def gen_class(la, lb, lc, ld, order, header):
    tag = f"{la}{lb}{lc}{ld}"
    lp, lq = la + lb, lc + ld
    np_, nq_ = len(herm(lp)), len(herm(lq))
    na, nb = (la + 1) * (la + 2) // 2, (lb + 1) * (lb + 2) // 2
    nc, nd = (lc + 1) * (lc + 2) // 2, (ld + 1) * (ld + 2) // 2
    nab, ncd = na * nb, nc * nd
    o = [header.format(what=f"Class ({la}{lb}|{lc}{ld}) MD Fock kernels through batched GEMMs")]
    o.append(f"module trc_mdg_k{tag}")
    o.append("   use pic_types, only: dp, int64")
    o.append("   use trc_boys, only: BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV")
    o.append("   use trc_linalg, only: trc_linalg_t")
    o.append("   implicit none")
    o.append("   private")
    o.append(f"   public :: mdg_k{tag}")
    o.append(f"   real(dp), parameter :: TWO_PI_2_5 = {TWO_PI_2_5}")
    o.append("   integer(int64), parameter :: MDG_BUDGET = 2_int64**26")
    o.append("   !! Doubles in one window's PQ, S and integral buffers together (512 MB).")
    o.append("")
    o.append("contains")
    o.append("")
    o.append('#include "inc/trc_boys_down.inc"')
    o.append("")
    o.extend(wrap_args(f"subroutine mdg_k{tag}", "lp, " + SEG, "   "))
    o.append(f"      !! Segment s of class ({la}{lb}|{lc}{ld}), window by window.")
    o.append("      type(trc_linalg_t), intent(inout) :: lp")
    o.append(SEG_DECL)
    o.append("      integer :: i0, ni, j0, nj, ka, kb, npb, nqb, t, nimax, njmax")
    o.append("      integer(int64) :: per")
    o.append("      real(dp), allocatable :: hb(:), hk(:), pq(:), sm(:), eri(:)")
    o.append("      ka = maxval(sq_n(sOA(s) + 1:sOA(s) + sA(s)))")
    o.append("      kb = maxval(sq_n(sOB(s) + 1:sOB(s) + sNB(s)))")
    o.append(f"      npb = {np_}*ka")
    o.append(f"      nqb = {nq_}*kb")
    o.append("      nimax = min(sA(s), 256)")
    o.append(f"      per = int(nimax, int64)*(int(npb, int64)*nqb + int(npb, int64)*{ncd} + {nab}*{ncd})")
    o.append("      njmax = int(max(1_int64, min(int(sNB(s), int64), MDG_BUDGET/per)))")
    o.append(f"      allocate (hb({nab}*npb*nimax), hk({ncd}*nqb*njmax), pq(npb*nimax*nqb*njmax), &")
    o.append(f"                sm(npb*nimax*{ncd}*njmax), eri({nab}*nimax*{ncd}*njmax))")
    o.append("      !$acc enter data create(hb, hk, pq, sm, eri)")
    o.append("      do i0 = 0, sA(s) - 1, nimax")
    o.append("         ni = min(nimax, sA(s) - i0)")
    o.append("         do concurrent(t = 1:ni*ka)")
    o.extend(wrap_args(f"call mdg_hb_item_{tag}", "t, " + SEG + ", i0, ni, 0, 1, ka, kb, hb", " " * 12))
    o.append("         end do")
    o.append("         do j0 = 0, sNB(s) - 1, njmax")
    o.append("            nj = min(njmax, sNB(s) - j0)")
    o.append("            do concurrent(t = 1:nj*kb)")
    o.extend(wrap_args(f"call mdg_hk_item_{tag}", "t, " + SEG + ", 0, 1, j0, nj, ka, kb, hk", " " * 15))
    o.append("            end do")
    o.append("            do concurrent(t = 1:ni*ka*nj*kb)")
    o.extend(wrap_args(f"call mdg_pq_item_{tag}", "t, " + SEG + ", " + WIN + ", pq", " " * 15))
    o.append("            end do")
    o.append("            ! S_j = PQ_j Hk_j^T, then E_i = Hb_i S(i, :) with S read in place.")
    o.append(f"            call lp%gemm_strided(\"N\", \"T\", npb*ni, {ncd}, nqb, 1.0_dp, pq, npb*ni, &")
    o.append("                                 int(npb, int64)*ni*nqb, hk, " + f"{ncd}, int({ncd}*nqb, int64), 0.0_dp, sm, npb*ni, &")
    o.append(f"                                 int(npb, int64)*ni*{ncd}, nj)")
    o.append(f"            call lp%gemm_strided(\"N\", \"N\", {nab}, {ncd}*nj, npb, 1.0_dp, hb, {nab}, &")
    o.append(f"                                 int({nab}*npb, int64), sm, npb*ni, int(npb, int64), 0.0_dp, eri, &")
    o.append(f"                                 {nab}*ni, int({nab}, int64), ni)")
    o.append(f"            do concurrent(t = 1:{nab}*ni*nj)")
    o.extend(wrap_args(f"call mdg_dig_item_{tag}", "t, " + SEG + ", " + WIN + ", eri", " " * 15))
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    o.append("      !$acc exit data delete(hb, hk, pq, sm, eri)")
    o.append("      deallocate (hb, hk, pq, sm, eri)")
    o.append(f"   end subroutine mdg_k{tag}")
    o.append("")
    o.extend(h_item(tag, "b", la, lb, order))
    o.append("")
    o.extend(h_item(tag, "k", lc, ld, order))
    o.append("")
    o.extend(pq_item(tag, lp, lq))
    o.append("")
    o.extend(dig_item(tag, na, nb, nc, nd))
    o.append("")
    o.append(f"end module trc_mdg_k{tag}")
    return "\n".join(o) + "\n"


SEGRUN = ("lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, npair, sp_i, sp_j, sp_q, thresh, pcut, "
          "jfac, kfac, dsh, nbas, nao, ao_off, nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, "
          "sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, ndens, dmat, jmat")


def gen_dispatch(classes, header):
    o = [header.format(what="Host dispatch of the MD GEMM Fock kernels")]
    o.append("module trc_mdg_kernels")
    o.append("   use pic_types, only: dp, int64")
    o.append("   use trc_linalg, only: trc_linalg_t")
    for c in classes:
        tag = "".join(map(str, c))
        o.append(f"   use trc_mdg_k{tag}, only: mdg_k{tag}")
    o.append("   implicit none")
    o.append("   private")
    o.append("   public :: mdg_has_class, mdg_class_run")
    o.append("")
    o.append("   type(trc_linalg_t), save :: lp")
    o.append("   !! cuBLAS for the whole run: made on the first class run, never released")
    o.append("   !! (a module variable never goes out of scope). Prototype.")
    o.append("")
    o.append("contains")
    o.append("")
    o.append("   pure logical function mdg_has_class(key) result(res)")
    o.append("      integer, intent(in) :: key")
    o.append("      res = .false.")
    for c in classes:
        o.append(f"      if (key == {class_key(*c)}) res = .true.")
    o.append("   end function mdg_has_class")
    o.append("")
    o.extend(wrap_args("subroutine mdg_class_run", "key, " + SEGRUN, "   "))
    o.append("      integer, intent(in) :: key")
    o.append("      integer,  intent(in)    :: lo, hi, nseg, npair, nbas, nao, nqp, ndens")
    o.append("      integer(int64), intent(in) :: sOff(nseg + 1)")
    o.append(SEG_DECL.split("\n", 1)[1])
    o.append("      integer :: s")
    o.append("      if (.not. lp%ready) call lp%init(1)")
    o.append("      do s = lo, hi")
    o.append("         select case (key)")
    for c in classes:
        tag = "".join(map(str, c))
        o.append(f"         case ({class_key(*c)})")
        o.extend(wrap_args(f"call mdg_k{tag}", "lp, " + SEG, " " * 12))
    o.append("         case default")
    o.append('            error stop "trc_mdg_kernels: no MD GEMM kernel for this class"')
    o.append("         end select")
    o.append("      end do")
    o.append("   end subroutine mdg_class_run")
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
        if c[0] < c[1] or c[2] < c[3]:
            raise SystemExit(f"class {c}: la >= lb and lc >= ld, as the bins store them")
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
