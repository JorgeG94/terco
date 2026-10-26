#!/usr/bin/env python3
"""Generate McMurchie-Davidson Fock (J and K) kernels for chosen classes.

A PROTOTYPE, one class at a time, raced against HGP/Rys on the same
segments: fock_bins hands a class run to `mdk_dispatch` instead of
pc_dispatch/rys_dispatch when TRC_MDK_CLASSES names it, with the same
segment descriptors, screens and folded digestion, so the Fock matrix is
directly comparable and so is the time.

Per window of quartets of one class (la lb|lc ld):

  kernel 1, one thread per quartet: the HGP kernels' own decode, Schwarz,
     density and primitive screens; then for every primitive quartet
     R^(0)_r, |r| <= L, with c_P c_Q 2 pi^2.5/(pq sqrt(p+q)) folded in, into
     a buffer whose fastest index is the quartet, so warps coalesce.
  kernel 2, one thread per (bra component pair ab, quartet), quartet
     fastest: a warp shares ab, so the select case on it is warp-uniform and
     every E and R index is a literal. Per primitive quartet
         v_q    = sum_p E^ab_p R_{p+q}                 (q <= lc + ld)
         eri_cd += sum_q (-1)^|q| E^cd_q v_q
     then the thread digests its row of (ab|cd): K always, J when jfac /= 0,
     with the HGP kernels' degeneracy weights and atomics.

Nothing here is a quartet of Cartesian integrals held by one thread: a
thread holds one row, v and the E tables.

Usage: gen_md_k.py --classes 2121[,...] [--order lexicographic] -o DIR
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cartorder  # noqa: E402
from md import rtree  # noqa: E402
from md.hermite import DIRS, herm, herm_index, e1d_statements, e_operand  # noqa: E402
from gen_md_j import (ORDER_CODE, TWO_PI_2_5, decl, emit_sum, expr_of, lit,  # noqa: E402
                      wrap_args, nbatch_for)

CLASS_RADIX = 11  # trc_binkernel's CLASS_RADIX


def class_key(la, lb, lc, ld):
    return ((la * CLASS_RADIX + lb) * CLASS_RADIX + lc) * CLASS_RADIX + ld


# Arguments every launch takes, the HGP dispatch interface plus the window.
SEG_ARGS = ("lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, "
            "npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, dsh, nbas, nao, ao_off, "
            "nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, "
            "ndens, dmat, jmat")

SEG_DECL = """      integer,  intent(in)    :: lo, hi, nseg, npair, nbas, nao, nqp, ndens
      integer(int64), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg), sPA(nseg), sPB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair), ao_off(nbas), sq_n(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh, pcut, jfac, kfac, dsh(nbas, nbas)
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp), qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(npair), sq_ay(npair), sq_az(npair)
      real(dp), intent(in)    :: sq_abx(npair), sq_aby(npair), sq_abz(npair)
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao)"""

WIN_ARGS = "g0, nw, kk, nr, rbuf, wseg, wab, wcd, wq"
WIN_DECL = """      integer(int64), intent(in) :: g0
      integer,  intent(in)    :: nw, kk, nr
      real(dp), intent(inout) :: rbuf(nw, nr, kk)
      integer,  intent(inout) :: wseg(nw), wab(nw), wcd(nw)
      real(dp), intent(inout) :: wq(nw)"""


def decode_lines(ind):
    """The HGP kernels' segment search and triangle inversion, verbatim logic."""
    return f"""{ind}p = lo; q = hi
{ind}do while (p < q)
{ind}   mid = (p + q + 1)/2
{ind}   if (sOff(mid) < gt) then
{ind}      p = mid
{ind}   else
{ind}      q = mid - 1
{ind}   end if
{ind}end do
{ind}seg = p
{ind}t = int(gt - sOff(seg))
{ind}if (sD(seg)) then
{ind}   nsa = int(sA(seg), int64)
{ind}   u = int(t - 1, int64)
{ind}   kx = int((real(2*nsa + 1, dp) - sqrt(real(2*nsa + 1, dp)**2 - 8.0_dp*real(u, dp)))/2.0_dp, int64)
{ind}   if (kx < 0) kx = 0
{ind}   do while (kx > 0)
{ind}      if ((kx*(2*nsa + 1) - kx*kx)/2 <= u) exit
{ind}      kx = kx - 1
{ind}   end do
{ind}   do while (((kx + 1)*(2*nsa + 1) - (kx + 1)*(kx + 1))/2 <= u)
{ind}      kx = kx + 1
{ind}   end do
{ind}   icd = int(kx) + 1
{ind}   iab = icd + int(u - (kx*(2*nsa + 1) - kx*kx)/2)
{ind}else
{ind}   icd = (t - 1)/sA(seg) + 1
{ind}   iab = t - (icd - 1)*sA(seg)
{ind}end if""".split("\n")


def gen_class(la, lb, lc, ld, order, header):
    tag = f"{la}{lb}{lc}{ld}"
    lp, lq = la + lb, lc + ld
    L = lp + lq
    hp, hq, hl = herm(lp), herm(lq), herm(L)
    il = herm_index(L)
    nr = len(hl)
    ca, cb = cartorder.components(la, order), cartorder.components(lb, order)
    cc, cd = cartorder.components(lc, order), cartorder.components(ld, order)
    nab = len(ca) * len(cb)
    o = [header.format(what=f"Class ({la}{lb}|{lc}{ld}) MD Fock kernels, prototype")]
    o.append(f"module trc_mdk_k{tag}")
    o.append("   use pic_types, only: dp, int64")
    o.append("   use trc_boys, only: BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV")
    o.append("   implicit none")
    o.append("   private")
    o.append(f"   public :: mdk_k{tag}, MDK_NR_{tag}, MDK_NAB_{tag}")
    o.append(f"   integer, parameter :: MDK_NR_{tag} = {nr}")
    o.append(f"   integer, parameter :: MDK_NAB_{tag} = {nab}")
    o.append(f"   real(dp), parameter :: TWO_PI_2_5 = {TWO_PI_2_5}")
    o.append("")
    o.append("contains")
    o.append("")
    o.append('#include "inc/trc_boys_eval.inc"')
    o.append("")
    # ---- window driver ----
    o.extend(wrap_args(f"subroutine mdk_k{tag}", SEG_ARGS + ", " + WIN_ARGS, "   "))
    o.append(f"      !! One window of quartets g0+1..g0+nw of class ({la}{lb}|{lc}{ld}); kk is the")
    o.append("      !! largest primitive-quartet count in it, rbuf is (nw, nr, kk) on the device.")
    o.append(SEG_DECL)
    o.append(WIN_DECL)
    o.append("      integer :: i")
    o.append("      do concurrent(i = 1:nw)")
    o.extend(wrap_args(f"call mdk_r_item_{tag}", "i, " + SEG_ARGS + ", " + WIN_ARGS, " " * 9))
    o.append("      end do")
    o.append(f"      do concurrent(i = 1:nw*{nab})")
    o.extend(wrap_args(f"call mdk_d_item_{tag}", "i, " + SEG_ARGS + ", " + WIN_ARGS, " " * 9))
    o.append("      end do")
    o.append(f"   end subroutine mdk_k{tag}")
    o.append("")
    o.extend(r_item(tag, L, nr, il))
    o.append("")
    o.extend(d_item(tag, la, lb, lc, ld, ca, cb, cc, cd, hp, hq, il, nr))
    o.append("")
    o.append(f"end module trc_mdk_k{tag}")
    return "\n".join(o) + "\n"


def r_item(tag, L, nr, il):
    g = rtree.RDag(herm(L))
    g.sieve()
    stmts = g.statements()
    leaves = g.leaves()
    o = []
    o.extend(wrap_args(f"pure subroutine mdk_r_item_{tag}", "i, " + SEG_ARGS + ", " + WIN_ARGS, "   "))
    o.append("      !$acc routine seq")
    o.append("      integer, intent(in) :: i")
    o.append(SEG_DECL)
    o.append(WIN_DECL)
    o.append("      integer(int64) :: gt, nsa, u, kx")
    o.append("      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl, slotab, slotcd")
    o.append("      integer :: k1, k2, kp, kq, kq0, nab, ncd, idx, r")
    o.append("      real(dp) :: qcut, w, zeta, eta, zpe, alpha, pref, tval, m2a, ww, bnd")
    o.append("      real(dp) :: pqx, pqy, pqz, px, py, pz")
    o.append("      real(dp) :: f(0:BOYS_MMAX)")
    o.append(decl("real(dp)", [rtree.leaf_name(m) for m in leaves]))
    o.append(decl("real(dp)", [s.lhs for s in stmts]))
    o.append("")
    o.append("      gt = g0 + i")
    o.extend(decode_lines("      "))
    o.append("      wseg(i) = seg")
    o.append("      slotab = sOA(seg) + iab; slotcd = sOB(seg) + icd")
    o.append("      wab(i) = slotab; wcd(i) = slotcd")
    o.append("      wq(i) = 0.0_dp")
    o.append("      nab = sq_n(slotab); ncd = sq_n(slotcd)")
    o.append("      qcut = sp_q(slotab)*sp_q(slotcd)")
    o.append("      si = sp_i(slotab); sj = sp_j(slotab); sk = sp_i(slotcd); sl = sp_j(slotcd)")
    o.append("      ! The HGP kernels' screens, so both engines drop the same quartets.")
    o.append("      if (qcut <= thresh) return")
    o.append("      if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), dsh(si, sk), dsh(si, sl), &")
    o.append("                   dsh(sj, sk), dsh(sj, sl)) <= thresh) return")
    o.append("      w = 1.0_dp")
    o.append("      if (si == sj) w = 0.5_dp*w")
    o.append("      if (sk == sl) w = 0.5_dp*w")
    o.append("      if (sD(seg) .and. iab == icd) w = 0.5_dp*w")
    o.append("      wq(i) = w")
    o.append("      kq0 = sPB(seg) + icd")
    o.append("      do k1 = 1, nab")
    o.append("         kp = sPA(seg) + (k1 - 1)*sA(seg) + iab")
    o.append("         zeta = qp_z(kp)")
    o.append("         px = qp_px(kp); py = qp_py(kp); pz = qp_pz(kp)")
    o.append("         bnd = TWO_PI_2_5*abs(qp_c(kp))/(zeta*sqrt(zeta))")
    o.append("         do k2 = 1, ncd")
    o.append("            idx = (k1 - 1)*ncd + k2")
    o.append("            kq = sPB(seg) + (k2 - 1)*sNB(seg) + icd")
    o.append("            eta = qp_z(kq)")
    o.append("            zpe = zeta + eta")
    o.append("            pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*qp_c(kp)*qp_c(kq)")
    o.append("            ! The HGP kernels' primitive-quartet prescreen; a skipped one is zero.")
    o.append("            if (bnd*abs(qp_c(kq))/eta <= pcut .or. abs(pref) <= pcut) then")
    o.append("               do r = 1, nr")
    o.append("                  rbuf(i, r, idx) = 0.0_dp")
    o.append("               end do")
    o.append("               cycle")
    o.append("            end if")
    o.append("            alpha = zeta*eta/zpe")
    o.append("            pqx = px - qp_px(kq); pqy = py - qp_py(kq); pqz = pz - qp_pz(kq)")
    o.append("            tval = alpha*(pqx*pqx + pqy*pqy + pqz*pqz)")
    o.append(f"            call boys_eval({max(leaves)}, tval, f)")
    o.append("            m2a = -2.0_dp*alpha")
    o.append("            ww = pref")
    for m in range(max(leaves) + 1):
        if m > 0:
            o.append("            ww = ww*m2a")
        if m in leaves:
            o.append(f"            {rtree.leaf_name(m)} = ww*f({m})")
    for s in stmts:
        o.extend(emit_sum(s.lhs, expr_of(s.terms), " " * 12))
    for r in herm(L):
        val = rtree.leaf_name(0) if r == (0, 0, 0) else rtree.r_name(r, 0)
        o.append(f"            rbuf(i, {il[r] + 1}, idx) = {val}")
    o.append("         end do")
    o.append("      end do")
    o.append(f"   end subroutine mdk_r_item_{tag}")
    return o


def e_names_for(la, lb, comps_a, comps_b):
    need = [set(), set(), set()]
    for a in comps_a:
        for b in comps_b:
            for d in range(3):
                for t in range(a[d] + b[d] + 1):
                    need[d].add((a[d], b[d], t))
    return need


def e_stmts(la, lb, need, ind, side):
    """1D E statements, renamed per side ('b' bra, 'k' ket)."""
    lines, names = [], []
    for d in range(3):
        for s in e1d_statements(la, lb, d, needed=need[d]):
            parts = []
            for tm in s.terms:
                f = []
                if tm.const != 1.0 or (tm.factor is None and tm.operand is None):
                    f.append(lit(tm.const))
                if tm.factor is not None:
                    f.append(side + tm.factor)
                if tm.operand is not None:
                    f.append(side + tm.operand)
                parts.append("*".join(f))
            lines.extend(emit_sum(side + s.lhs, parts, ind))
            names.append(side + s.lhs)
    return lines, names


def e3(side, a, b, h):
    f = []
    for d in range(3):
        if h[d] > a[d] + b[d]:
            return None
        n = e_operand(d, a[d], b[d], h[d])
        if n is not None:
            f.append(side + n)
    return f


def d_item(tag, la, lb, lc, ld, ca, cb, cc, cd, hp, hq, il, nr):
    nab = len(ca) * len(cb)
    ncd = len(cc) * len(cd)
    need_k = e_names_for(lc, ld, cc, cd)
    k_lines, k_names = e_stmts(lc, ld, need_k, " " * 12, "k")
    o = []
    o.extend(wrap_args(f"pure subroutine mdk_d_item_{tag}", "i, " + SEG_ARGS + ", " + WIN_ARGS, "   "))
    o.append("      !$acc routine seq")
    o.append("      integer, intent(in) :: i")
    o.append(SEG_DECL)
    o.append(WIN_DECL)
    o.append("      integer :: w, ab, seg, slotab, slotcd, si, sj, sk, sl, mui, nuj, lamk, sigl")
    o.append("      integer :: k1, k2, kp, kq, nab, ncd, idx, idens")
    o.append("      real(dp) :: wqq, zeta, eta, rp, rq")
    o.append("      real(dp) :: bpxa, bpya, bpza, bpxb, bpyb, bpzb, boo2p")
    o.append("      real(dp) :: kpxa, kpya, kpza, kpxb, kpyb, kpzb, koo2p")
    o.append("      real(dp) :: s, dab")
    o.append(decl("real(dp)", k_names))
    # bra E names are per branch but declare the union
    need_b = e_names_for(la, lb, ca, cb)
    _, b_names = e_stmts(la, lb, need_b, "", "b")
    o.append(decl("real(dp)", b_names))
    o.append(decl("real(dp)", [f"v{j + 1}" for j in range(len(hq))]))
    o.append(decl("real(dp)", [f"x{c + 1}_{d + 1}" for c in range(len(cc)) for d in range(len(cd))]))
    o.append("")
    o.append("      ! Quartet fastest: a warp shares one bra component pair, so the select")
    o.append("      ! case below does not diverge.")
    o.append("      w = mod(i - 1, nw) + 1")
    o.append("      ab = (i - 1)/nw + 1")
    o.append("      wqq = wq(w)")
    o.append("      if (wqq == 0.0_dp) return")
    o.append("      seg = wseg(w); slotab = wab(w); slotcd = wcd(w)")
    o.append("      si = sp_i(slotab); sj = sp_j(slotab); sk = sp_i(slotcd); sl = sp_j(slotcd)")
    o.append("      mui = ao_off(si); nuj = ao_off(sj); lamk = ao_off(sk); sigl = ao_off(sl)")
    o.append("      nab = sq_n(slotab); ncd = sq_n(slotcd)")
    o.append("      select case (ab)")
    k = 0
    for xa, a in enumerate(ca):
        for xb, b in enumerate(cb):
            k += 1
            o.append(f"      case ({k})")
            o.extend(branch(la, lb, lc, ld, a, b, xa, xb, cc, cd, hp, hq, il, k_lines))
    o.append("      end select")
    o.append(f"   end subroutine mdk_d_item_{tag}")
    return o


def branch(la, lb, lc, ld, a, b, xa, xb, cc, cd, hp, hq, il, k_lines):
    ind = " " * 9
    need_b = [set(), set(), set()]
    for d in range(3):
        for t in range(a[d] + b[d] + 1):
            need_b[d].add((a[d], b[d], t))
    b_lines, _ = e_stmts(la, lb, need_b, " " * 9, "b")
    box = [h for h in hp if all(h[d] <= a[d] + b[d] for d in range(3))]
    xs = [f"x{c + 1}_{d + 1}" for c in range(len(cc)) for d in range(len(cd))]
    o = []
    for x in xs:
        o.append(f"{ind}{x} = 0.0_dp")
    o.append(f"{ind}do k1 = 1, nab")
    o.append(f"{ind}   kp = sPA(seg) + (k1 - 1)*sA(seg) + (slotab - sOA(seg))")
    o.append(f"{ind}   zeta = qp_z(kp)")
    o.append(f"{ind}   rp = 1.0_dp/zeta")
    o.append(f"{ind}   boo2p = 0.5_dp*rp")
    o.append(f"{ind}   bpxa = qp_px(kp) - sq_ax(slotab); bpya = qp_py(kp) - sq_ay(slotab)")
    o.append(f"{ind}   bpza = qp_pz(kp) - sq_az(slotab)")
    o.append(f"{ind}   bpxb = bpxa + sq_abx(slotab); bpyb = bpya + sq_aby(slotab); bpzb = bpza + sq_abz(slotab)")
    o.extend(["   " + ln for ln in b_lines])
    o.append(f"{ind}   do k2 = 1, ncd")
    o.append(f"{ind}      idx = (k1 - 1)*ncd + k2")
    o.append(f"{ind}      kq = sPB(seg) + (k2 - 1)*sNB(seg) + (slotcd - sOB(seg))")
    o.append(f"{ind}      eta = qp_z(kq)")
    o.append(f"{ind}      rq = 1.0_dp/eta")
    o.append(f"{ind}      koo2p = 0.5_dp*rq")
    o.append(f"{ind}      kpxa = qp_px(kq) - sq_ax(slotcd); kpya = qp_py(kq) - sq_ay(slotcd)")
    o.append(f"{ind}      kpza = qp_pz(kq) - sq_az(slotcd)")
    o.append(f"{ind}      kpxb = kpxa + sq_abx(slotcd); kpyb = kpya + sq_aby(slotcd); kpzb = kpza + sq_abz(slotcd)")
    o.extend(k_lines)
    # v_q = sum_p E^ab_p R_{p+q}
    for j, q in enumerate(hq):
        parts = []
        for h in box:
            f = e3("b", a, b, h)
            r = tuple(h[d] + q[d] for d in range(3))
            parts.append("*".join(f + [f"rbuf(w, {il[r] + 1}, idx)"]))
        o.extend(emit_sum(f"v{j + 1}", parts, ind + "      "))
    for xc, c in enumerate(cc):
        for xd, d in enumerate(cd):
            parts = []
            for j, q in enumerate(hq):
                f = e3("k", c, d, q)
                if f is None:
                    continue
                sg = "-" if sum(q) % 2 else ""
                term = "*".join(f + [f"v{j + 1}"])
                parts.append(term if not sg else f"(-{term})")
            o.extend(emit_sum(f"x{xc + 1}_{xd + 1}", parts, ind + "      ", accumulate=True))
    o.append(f"{ind}   end do")
    o.append(f"{ind}end do")
    # digestion: this thread's row (a, b) of (ab|cd)
    o.append(f"{ind}do idens = 1, ndens")
    for xc in range(len(cc)):
        terms = [f"x{xc + 1}_{xd + 1}*dmat(idens, nuj + {xb}, sigl + {xd})" for xd in range(len(cd))]
        o.extend(emit_sum("s", terms, ind + "   "))
        o.append(f"{ind}   !$acc atomic update")
        o.append(f"{ind}   jmat(idens, mui + {xa}, lamk + {xc}) = jmat(idens, mui + {xa}, lamk + {xc}) - kfac*wqq*s")
        terms = [f"x{xc + 1}_{xd + 1}*dmat(idens, mui + {xa}, sigl + {xd})" for xd in range(len(cd))]
        o.extend(emit_sum("s", terms, ind + "   "))
        o.append(f"{ind}   !$acc atomic update")
        o.append(f"{ind}   jmat(idens, nuj + {xb}, lamk + {xc}) = jmat(idens, nuj + {xb}, lamk + {xc}) - kfac*wqq*s")
    for xd in range(len(cd)):
        terms = [f"x{xc + 1}_{xd + 1}*dmat(idens, nuj + {xb}, lamk + {xc})" for xc in range(len(cc))]
        o.extend(emit_sum("s", terms, ind + "   "))
        o.append(f"{ind}   !$acc atomic update")
        o.append(f"{ind}   jmat(idens, mui + {xa}, sigl + {xd}) = jmat(idens, mui + {xa}, sigl + {xd}) - kfac*wqq*s")
        terms = [f"x{xc + 1}_{xd + 1}*dmat(idens, mui + {xa}, lamk + {xc})" for xc in range(len(cc))]
        o.extend(emit_sum("s", terms, ind + "   "))
        o.append(f"{ind}   !$acc atomic update")
        o.append(f"{ind}   jmat(idens, nuj + {xb}, sigl + {xd}) = jmat(idens, nuj + {xb}, sigl + {xd}) - kfac*wqq*s")
    o.append(f"{ind}   if (jfac /= 0.0_dp) then")
    terms = [f"x{xc + 1}_{xd + 1}*dmat(idens, lamk + {xc}, sigl + {xd})"
             for xc in range(len(cc)) for xd in range(len(cd))]
    o.extend(emit_sum("s", terms, ind + "      "))
    o.append(f"{ind}      !$acc atomic update")
    o.append(f"{ind}      jmat(idens, mui + {xa}, nuj + {xb}) = jmat(idens, mui + {xa}, nuj + {xb}) "
             "+ 4.0_dp*jfac*wqq*s")
    o.append(f"{ind}      dab = 4.0_dp*jfac*wqq*dmat(idens, mui + {xa}, nuj + {xb})")
    for xc in range(len(cc)):
        for xd in range(len(cd)):
            o.append(f"{ind}      !$acc atomic update")
            o.append(f"{ind}      jmat(idens, lamk + {xc}, sigl + {xd}) = jmat(idens, lamk + {xc}, sigl + {xd}) "
                     f"+ dab*x{xc + 1}_{xd + 1}")
    o.append(f"{ind}   end if")
    o.append(f"{ind}end do")
    return o


def gen_dispatch(classes, header):
    o = [header.format(what="Host dispatch of the MD Fock prototype kernels")]
    o.append("module trc_mdk_kernels")
    o.append("   use pic_types, only: dp, int64")
    for (la, lb, lc, ld) in classes:
        tag = f"{la}{lb}{lc}{ld}"
        o.append(f"   use trc_mdk_k{tag}, only: mdk_k{tag}, MDK_NR_{tag}, MDK_NAB_{tag}")
    o.append("   implicit none")
    o.append("   private")
    o.append("   public :: mdk_has_class, mdk_class_nr, mdk_window, mdk_class_run")
    o.append("")
    o.append("   integer(int64), parameter :: MDK_RBUF = 2_int64**25")
    o.append("   !! Doubles in the Hermite buffer of one window (256 MB).")
    o.append("")
    o.append("contains")
    o.append("")
    o.append("   pure logical function mdk_has_class(key) result(res)")
    o.append("      integer, intent(in) :: key")
    o.append("      res = .false.")
    for c in classes:
        o.append(f"      if (key == {class_key(*c)}) res = .true.")
    o.append("   end function mdk_has_class")
    o.append("")
    o.append("   pure integer function mdk_class_nr(key) result(nr)")
    o.append("      !! Hermite integrals per primitive quartet of the class.")
    o.append("      integer, intent(in) :: key")
    o.append("      nr = 0")
    for c in classes:
        o.append(f"      if (key == {class_key(*c)}) nr = MDK_NR_{''.join(map(str, c))}")
    o.append("   end function mdk_class_nr")
    o.append("")
    o.extend(wrap_args("subroutine mdk_window", "key, " + SEG_ARGS + ", " + WIN_ARGS, "   "))
    o.append("      integer, intent(in) :: key")
    o.append(SEG_DECL)
    o.append(WIN_DECL)
    o.append("      select case (key)")
    for c in classes:
        tag = "".join(map(str, c))
        o.append(f"      case ({class_key(*c)})")
        o.extend(wrap_args(f"call mdk_k{tag}", SEG_ARGS + ", " + WIN_ARGS, " " * 9))
    o.append("      case default")
    o.append('         error stop "trc_mdk_kernels: no MD kernel for this class"')
    o.append("      end select")
    o.append("   end subroutine mdk_window")
    o.append("")
    o.extend(wrap_args("subroutine mdk_class_run", "key, " + SEG_ARGS, "   "))
    o.append("      !! Every quartet of segments lo..hi, one buffer-sized window at a time.")
    o.append("      integer, intent(in) :: key")
    o.append(SEG_DECL)
    o.append("      integer(int64) :: g, glast, n64")
    o.append("      integer :: s, kk, nr, nw, n, ka, kb, j")
    o.append("      real(dp), allocatable :: rbuf(:)")
    o.append("      integer,  allocatable :: wseg(:), wab(:), wcd(:)")
    o.append("      real(dp), allocatable :: wq(:)")
    o.append("      if (sOff(hi + 1) <= sOff(lo)) return")
    o.append("      kk = 1")
    o.append("      do s = lo, hi")
    o.append("         ka = 1")
    o.append("         kb = 1")
    o.append("         do j = sOA(s) + 1, sOA(s) + sA(s)")
    o.append("            ka = max(ka, sq_n(j))")
    o.append("         end do")
    o.append("         do j = sOB(s) + 1, sOB(s) + sNB(s)")
    o.append("            kb = max(kb, sq_n(j))")
    o.append("         end do")
    o.append("         kk = max(kk, ka*kb)")
    o.append("      end do")
    o.append("      nr = mdk_class_nr(key)")
    o.append("      n64 = max(1_int64, MDK_RBUF/int(nr*kk, int64))")
    o.append("      nw = int(min(n64, sOff(hi + 1) - sOff(lo)))")
    o.append("      allocate (rbuf(int(nw, int64)*nr*kk), wseg(nw), wab(nw), wcd(nw), wq(nw))")
    o.append("      !$acc enter data create(rbuf, wseg, wab, wcd, wq)")
    o.append("      g = sOff(lo)")
    o.append("      glast = sOff(hi + 1)")
    o.append("      do while (g < glast)")
    o.append("         n = int(min(int(nw, int64), glast - g))")
    o.extend(wrap_args("call mdk_window", "key, " + SEG_ARGS + ", g, n, kk, nr, rbuf, wseg, wab, wcd, wq", " " * 9))
    o.append("         g = g + n")
    o.append("      end do")
    o.append("      !$acc exit data delete(rbuf, wseg, wab, wcd, wq)")
    o.append("      deallocate (rbuf, wseg, wab, wcd, wq)")
    o.append("   end subroutine mdk_class_run")
    o.append("")
    o.append("end module trc_mdk_kernels")
    return "\n".join(o) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--classes", required=True, help="comma list like 2121,2120")
    ap.add_argument("--lmax", type=int, required=True)
    ap.add_argument("--order", default="lexicographic", choices=sorted(ORDER_CODE))
    ap.add_argument("-o", "--outdir", required=True)
    a = ap.parse_args()
    classes = [tuple(int(ch) for ch in c) for c in filter(None, a.classes.split(","))]
    classes = [c for c in classes if max(c) <= a.lmax]
    for c in classes:
        if c[0] < c[1] or c[2] < c[3]:
            raise SystemExit(f"class {c}: la >= lb and lc >= ld, as the bins store them")
    cmd = f"scripts/gen_md_k.py --lmax {a.lmax} --classes {a.classes}"
    if a.order not in ("lexicographic", "libcint"):
        cmd += f" --order {a.order}"
    code = ORDER_CODE[a.order]
    header = ("!\n! {what}.\n!\n"
              f"! GENERATED by {cmd} -- do not edit.\n!\n"
              f"#if TRC_CART_ORDER != {code}\n"
              '#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."\n'
              "#endif\n")
    os.makedirs(a.outdir, exist_ok=True)
    for c in classes:
        with open(os.path.join(a.outdir, f"trc_mdk_k{''.join(map(str, c))}.F90"), "w") as fh:
            fh.write(gen_class(*c, a.order, header))
    with open(os.path.join(a.outdir, "trc_mdk_kernels.F90"), "w") as fh:
        fh.write(gen_dispatch(classes, header))
    return 0


if __name__ == "__main__":
    sys.exit(main())
