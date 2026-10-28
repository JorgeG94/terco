#!/usr/bin/env python3
"""Generate the McMurchie-Davidson J engine kernels.

J_ab = sum_cd (ab|cd) D_cd is built in Hermite space and never forms a
quartet (Palethorpe and Barca, JCTC 21, 9388 (2025); Shao and Head-Gordon's
ODLC ordering):

  1. ket, per primitive pair Q of shells c, d (one kernel per (lc, ld)):
         Dt[Q]_tau = (-1)^|tau| c_Q sum_cd E^cd_tau D'_cd,
     D' = D_cd + D_dc off the shell diagonal. Also the exact Coulomb norm
     ||rho_Q|| of that primitive pair's density, for the screen.
  2. Hermite, per primitive pair P against every Q (one kernel per
     (lp, lq) = (la + lb, lc + ld), one launch per R batch):
         Jt[P]_t += sum_Q sum_tau R_{t + tau}(P, Q) Dt[Q]_tau,
     skipping (P, Q) when I_P ||rho_Q|| <= cut. That product bounds the
     pair's contribution to any J element by Cauchy-Schwarz in the Coulomb
     metric, so the screen is rigorous and density-weighted at once.
  3. bra, per shell pair (one kernel per (la, lb)):
         J_ab = sum_P c_P sum_t E^ab_t Jt[P]_t.

Stage 2 is the only O(N^2)-in-pairs work and the only kernel that sees
angular momentum on all four centres, and it sees only lp + lq: (dp|dp) and
(ds|dd) share the (3|... ) and (4|...) kernels. One thread owns one bra
primitive pair and loops over the ket, so it accumulates Jt[P] in registers
and writes it once per launch: no atomics, deterministic, and the batches
of a class are consecutive launches over the same arrays.

Usage:
    gen_md_j.py --lmax 2 [--order lexicographic] [--batches 7:6,8:6] -o DIR

Writes trc_mdj_e<la><lb>.F90, trc_mdj_r<lp><lq>.F90 and the host
dispatcher trc_mdj_kernels.F90. Every kernel is checked first by
md/reference.py, which runs the same statement lists in numpy against
pyscf.

REFERENCES
----------
Inspiration only; the formulation and the code are terco's own.

  L. E. McMurchie, E. R. Davidson, J. Comput. Phys. 26, 218 (1978).
  E. Palethorpe, G. M. J. Barca, J. Chem. Theory Comput. 21, 9388 (2025).
  I. S. Ufimtsev, T. J. Martinez, J. Chem. Theory Comput. 4, 222 (2008);
    5, 1004 (2009).
  A. Asadchev, E. F. Valeev, J. Phys. Chem. A 127, 10889 (2023);
    J. Chem. Phys. 160, 244109 (2024).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cartorder  # noqa: E402
from md import rtree  # noqa: E402
from md.hermite import DIRS, herm, herm_index, e1d_statements, e_operand  # noqa: E402

ORDER_CODE = {"lexicographic": 0, "libcint": 0, "gamess": 1}

# Batches per total L when --batches does not say. Palethorpe and Barca's
# fastest counts on an A100 (their Table 4 / Figure 8 bolding): one batch up
# to L = 4, three at 5-6, six at 7-8, nine from 9. To be re-measured here.
DEFAULT_BATCHES = {5: 3, 6: 3, 7: 6, 8: 6}
DEFAULT_BATCHES_HIGH = 9

TWO_PI_2_5 = "34.986836655249725_dp"


def nbatch_for(L, table):
    if L in table:
        return table[L]
    if L >= 9:
        return DEFAULT_BATCHES_HIGH
    return DEFAULT_BATCHES.get(L, 1)


def lit(c):
    """A Fortran real literal for a float that is exactly representable."""
    if c == int(c):
        return f"{int(c)}.0_dp"
    return f"{c!r}_dp"


def dfact(n):
    """(n - 1)!! for n >= 0 even: 1, 1, 3, 15, ..."""
    out = 1
    for k in range(n - 1, 0, -2):
        out *= k
    return out


def expr_of(terms):
    parts = []
    for t in terms:
        f = []
        if t.const != 1.0 or (t.factor is None and t.operand is None):
            f.append(lit(t.const))
        if t.factor is not None:
            f.append(t.factor)
        if t.operand is not None:
            f.append(t.operand)
        parts.append("*".join(f))
    return parts


def emit_sum(lhs, parts, indent, accumulate=False, chunk=12):
    """lhs = p1 + p2 + ..., chunked so no statement needs many continuations."""
    out = []
    if not parts:
        if not accumulate:
            out.append(f"{indent}{lhs} = 0.0_dp")
        return out
    first = True
    for k in range(0, len(parts), chunk):
        piece = parts[k:k + chunk]
        head = f"{lhs} = " if (first and not accumulate) else f"{lhs} = {lhs} + "
        first = False
        lines = []
        cur = indent + head
        for j, p in enumerate(piece):
            tok = p if j == 0 else " + " + p
            if len(cur) + len(tok) > 110:
                lines.append(cur + " &")
                cur = indent + "   " + tok.lstrip()
            else:
                cur += tok
        lines.append(cur)
        # A continuation line that begins with "+" is fine Fortran.
        out.extend(lines)
    return out


def wrap_args(head, args, ind, width=100):
    """`head(args)` over continuation lines, broken at commas."""
    words = args.split(", ")
    lines, cur = [], ind + head + "("
    for k, w in enumerate(words):
        tok = w + (", " if k < len(words) - 1 else ")")
        if len(cur) + len(tok) > width:
            lines.append(cur.rstrip() + " &")
            cur = ind + " " * (len(head) + 1)
        cur += tok
    lines.append(cur)
    return lines


def comps(l, order):
    return cartorder.components(l, order)


# ---------------------------------------------------------------- E stages --

def e_needed(la, lb, order):
    """(i, j, t) per direction a (la, lb) class reads."""
    need = [set(), set(), set()]
    for a in comps(la, order):
        for b in comps(lb, order):
            for d in range(3):
                for t in range(a[d] + b[d] + 1):
                    need[d].add((a[d], b[d], t))
    return need


def e_block(la, lb, order, indent):
    need = e_needed(la, lb, order)
    lines = []
    names = []
    for d in range(3):
        for s in e1d_statements(la, lb, d, needed=need[d]):
            lines.extend(emit_sum(s.lhs, expr_of(s.terms), indent))
            names.append(s.lhs)
    return lines, names


def e3_factors(a, b, h):
    """The E^{ab}_h product as a list of factor names (empty means 1), or None."""
    f = []
    for d in range(3):
        if h[d] > a[d] + b[d]:
            return None
        n = e_operand(d, a[d], b[d], h[d])
        if n is not None:
            f.append(n)
    return f


def gen_e_module_loops(la, lb, order, header):
    """The ket, bra and bound stages of one (la, lb) type as loops over literal tables.

    These stages are O(pairs) per build and a few per cent of the J time at
    most, but written out they were the slowest part of the generated set to
    compile: the (gg) module took 753 s and 2 GB with every component pair and
    every Hermite product unrolled. Here only the tables are generated.
    """
    from gen_md_k import e1d_loop, int_table  # noqa: E402 (shared emitters)
    lp = la + lb
    tag = f"{la}{lb}"
    hs = herm(lp)
    nh = len(hs)
    ca, cb = comps(la, order), comps(lb, order)
    na, nb = len(ca), len(cb)
    # (h, h') pairs of the Coulomb self-norm, h <= h', all-even sums.
    pi, pj, ps, pw = [], [], [], []
    for i, h in enumerate(hs):
        for j in range(i, nh):
            h2 = hs[j]
            r = tuple(h[d] + h2[d] for d in range(3))
            if any(x % 2 for x in r):
                continue
            w = (1 if i == j else 2) * dfact(r[0]) * dfact(r[1]) * dfact(r[2])
            if sum(h) % 2:
                w = -w
            pi.append(i + 1)
            pj.append(j + 1)
            ps.append(sum(r) // 2)
            pw.append(w)
    npr = len(pi)

    def real_tab(name, vals):
        toks = [lit(float(v)) for v in vals]
        lines, cur = [], f"   real(dp), parameter :: {name}({len(vals)}) = ["
        for tok in toks:
            if len(cur) + len(tok) > 100:
                lines.append(cur + " &")
                cur = "      "
            cur += tok + ", "
        lines.append(cur.rstrip(", ") + "]")
        return lines

    o = [header.format(what=f"Class ({la}{lb}) ket, bra and bound stages of the MD J engine")]
    o.append(f"module trc_mdj_e{tag}")
    o.append("   use pic_types, only: dp")
    o.append("   implicit none")
    o.append("   private")
    o.append(f"   public :: mdj_ket_{tag}, mdj_bra_{tag}, mdj_bound_{tag}")
    o.append(f"   real(dp), parameter :: TWO_PI_2_5 = {TWO_PI_2_5}")
    o.extend(int_table("HT", [h[0] for h in hs]))
    o.extend(int_table("HU", [h[1] for h in hs]))
    o.extend(int_table("HV", [h[2] for h in hs]))
    o.extend(real_tab("SG", [(-1.0) ** sum(h) for h in hs]))
    for nm, cs in (("A", ca), ("B", cb)):
        for d, ax in enumerate("XYZ"):
            o.extend(int_table(f"C{nm}{ax}", [c[d] for c in cs]))
    o.extend(int_table("PI", pi))
    o.extend(int_table("PJ", pj))
    o.extend(int_table("PS", ps))
    o.extend(real_tab("PW", pw))
    o.extend(real_tab("RO", [1.0 / (2 * s_ + 1) for s_ in range(lp + 1)]))
    o.append("")
    o.append("contains")
    o.append("")
    ket_args = ("n, k0, kc, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, mp_sp, "
                "nsp, sp_ao1, sp_ao2, nao, dmat, ndt, dt, mp_nq")
    bra_args = ("n, s0, kl0, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, "
                "nsp, sp_ao1, sp_ao2, sp_k0, sp_nk, ndt, jt, nao, jmat")
    bnd_args = "n, k0, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, mp_ip"
    common = """      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)
      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)"""
    ket_decl = f"""      integer,  intent(in) :: n, k0, kc, nkc, hb, nhp, nsp, nao, ndt
{common}
      integer,  intent(in) :: mp_sp(nhp), sp_ao1(nsp), sp_ao2(nsp)
      real(dp), intent(in) :: dmat(nao, nao)
      real(dp), intent(inout) :: dt(ndt), mp_nq(nhp)"""
    bra_decl = f"""      integer,  intent(in) :: n, s0, kl0, nkc, hb, nhp, nsp, ndt, nao
{common}
      integer,  intent(in) :: sp_ao1(nsp), sp_ao2(nsp), sp_k0(nsp), sp_nk(nsp)
      real(dp), intent(in) :: jt(ndt)
      real(dp), intent(inout) :: jmat(nao, nao)"""
    bnd_decl = f"""      integer,  intent(in) :: n, k0, nhp
{common}
      real(dp), intent(inout) :: mp_ip(nhp)"""
    for kind, args, dcl in (("ket", ket_args, ket_decl), ("bra", bra_args, bra_decl), ("bound", bnd_args, bnd_decl)):
        o.extend(wrap_args(f"subroutine mdj_{kind}_{tag}", args, "   "))
        o.append(dcl)
        o.append("      integer :: t")
        o.append("      do concurrent(t = 1:n)")
        o.extend(wrap_args(f"call mdj_{kind}_item_{tag}", "t, " + args, " " * 9))
        o.append("      end do")
        o.append(f"   end subroutine mdj_{kind}_{tag}")
        o.append("")
    box = (f"HT(hh) <= CAX(x1) + CBX(x2) .and. HU(hh) <= CAY(x1) + CBY(x2) .and. "
           f"HV(hh) <= CAZ(x1) + CBZ(x2)")
    eprod = "e(CAX(x1), CBX(x2), HT(hh), 1)*e(CAY(x1), CBY(x2), HU(hh), 2)*e(CAZ(x1), CBZ(x2), HV(hh), 3)"
    e_setup = ["      oo2p = 0.5_dp/p",
               "      xa(1) = mp_ax(k); xa(2) = mp_ay(k); xa(3) = mp_az(k)",
               "      xb(1) = mp_bx(k); xb(2) = mp_by(k); xb(3) = mp_bz(k)"]
    norm = ["      do s = 0, " + str(lp), "         acc(s) = 0.0_dp", "      end do",
            f"      do q = 1, {npr}",
            "         acc(PS(q)) = acc(PS(q)) + PW(q)*g(PI(q))*g(PJ(q))",
            "      end do",
            "      nq2 = acc(0)", "      ws = 1.0_dp",
            f"      do s = 1, {lp}", "         ws = -p*ws", "         nq2 = nq2 + ws*acc(s)*RO(s + 1)", "      end do",
            "      nq2 = nq2*TWO_PI_2_5/(p*p*sqrt(2.0_dp*p))"]
    locals_ = [f"      real(dp) :: e(0:{la}, 0:{lb}, 0:{lp}, 3), g({nh}), acc(0:{lp})",
               "      real(dp) :: p, oo2p, xa(3), xb(3), v, nq2, ws",
               "      integer :: x1, x2, hh, q, s, d, ii, jj, tt"]
    # ---- ket item ----
    o.extend(wrap_args(f"pure subroutine mdj_ket_item_{tag}", "i, " + ket_args, "   "))
    o.append("      !$acc routine seq")
    o.append("      integer,  intent(in) :: i")
    o.append(ket_decl)
    o.extend(locals_)
    o.append("      integer :: k, sp, mu0, nu0")
    o.append("      real(dp) :: offd, dab, cs")
    o.append("      k = k0 + i")
    o.append("      sp = mp_sp(k)")
    o.append("      mu0 = sp_ao1(sp) - 1")
    o.append("      nu0 = sp_ao2(sp) - 1")
    o.append("      ! One shell pair stands for both orders off the diagonal: D_cd + D_dc.")
    o.append("      offd = 1.0_dp")
    o.append("      if (mu0 == nu0) offd = 0.0_dp")
    o.append("      p = mp_z(k)")
    o.extend(e_setup)
    o.extend(e1d_loop(la, lb, "      "))
    o.append(f"      do hh = 1, {nh}")
    o.append("         g(hh) = 0.0_dp")
    o.append("      end do")
    o.append(f"      do x2 = 1, {nb}")
    o.append(f"         do x1 = 1, {na}")
    o.append("            dab = dmat(mu0 + x1, nu0 + x2) + offd*dmat(nu0 + x2, mu0 + x1)")
    o.append(f"            do hh = 1, {nh}")
    o.append(f"               if ({box}) g(hh) = g(hh) + {eprod}*dab")
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    o.append("      ! Sign (-1)^|tau| and the pair coefficient (with K_cd) folded in.")
    o.append("      cs = mp_c(k)")
    o.append(f"      do hh = 1, {nh}")
    o.append("         g(hh) = cs*SG(hh)*g(hh)")
    o.append("         dt(hb + (hh - 1)*nkc + kc + i) = g(hh)")
    o.append("      end do")
    o.extend(norm)
    o.append("      mp_nq(k) = sqrt(max(nq2, 0.0_dp))")
    o.append(f"   end subroutine mdj_ket_item_{tag}")
    o.append("")
    # ---- bra item ----
    o.extend(wrap_args(f"pure subroutine mdj_bra_item_{tag}", "i, " + bra_args, "   "))
    o.append("      !$acc routine seq")
    o.append("      integer,  intent(in) :: i")
    o.append(bra_decl)
    o.extend(locals_)
    o.append(f"      real(dp) :: jab({na}, {nb})")
    o.append("      integer :: s2, mu0, nu0, kk, k")
    o.append("      real(dp) :: cs")
    o.append("      s2 = s0 + i")
    o.append("      mu0 = sp_ao1(s2) - 1")
    o.append("      nu0 = sp_ao2(s2) - 1")
    o.append(f"      do x2 = 1, {nb}")
    o.append(f"         do x1 = 1, {na}")
    o.append("            jab(x1, x2) = 0.0_dp")
    o.append("         end do")
    o.append("      end do")
    o.append("      do kk = 1, sp_nk(s2)")
    o.append("         k = sp_k0(s2) + kk")
    o.append("         p = mp_z(k)")
    o.extend(["   " + ln for ln in e_setup])
    o.extend(e1d_loop(la, lb, "         "))
    o.append("         cs = mp_c(k)")
    o.append(f"         do hh = 1, {nh}")
    o.append("            g(hh) = cs*jt(hb + (hh - 1)*nkc + k - kl0)")
    o.append("         end do")
    o.append(f"         do x2 = 1, {nb}")
    o.append(f"            do x1 = 1, {na}")
    o.append(f"               do hh = 1, {nh}")
    o.append(f"                  if ({box}) jab(x1, x2) = jab(x1, x2) + {eprod}*g(hh)")
    o.append("               end do")
    o.append("            end do")
    o.append("         end do")
    o.append("      end do")
    o.append("      ! Shell pairs are canonical and distinct, so no other thread writes")
    o.append("      ! either block: no atomics.")
    o.append(f"      do x2 = 1, {nb}")
    o.append(f"         do x1 = 1, {na}")
    o.append("            jmat(mu0 + x1, nu0 + x2) = jab(x1, x2)")
    o.append("            if (mu0 /= nu0) jmat(nu0 + x2, mu0 + x1) = jab(x1, x2)")
    o.append("         end do")
    o.append("      end do")
    o.append(f"   end subroutine mdj_bra_item_{tag}")
    o.append("")
    # ---- bound item ----
    o.extend(wrap_args(f"pure subroutine mdj_bound_item_{tag}", "i, " + bnd_args, "   "))
    o.append("      !$acc routine seq")
    o.append("      integer,  intent(in) :: i")
    o.append(bnd_decl)
    o.extend(locals_)
    o.append("      integer :: k")
    o.append("      real(dp) :: best")
    o.append("      k = k0 + i")
    o.append("      p = mp_z(k)")
    o.extend(e_setup)
    o.extend(e1d_loop(la, lb, "      "))
    o.append("      best = 0.0_dp")
    o.append(f"      do x2 = 1, {nb}")
    o.append(f"         do x1 = 1, {na}")
    o.append(f"            do hh = 1, {nh}")
    o.append("               g(hh) = 0.0_dp")
    o.append(f"               if ({box}) g(hh) = {eprod}")
    o.append("            end do")
    o.extend(["      " + ln for ln in norm])
    o.append("            best = max(best, nq2)")
    o.append("         end do")
    o.append("      end do")
    o.append("      mp_ip(k) = abs(mp_c(k))*sqrt(best)")
    o.append(f"   end subroutine mdj_bound_item_{tag}")
    o.append("")
    o.append(f"end module trc_mdj_e{tag}")
    return "\n".join(o) + "\n"


def gen_e_module(la, lb, order, header):
    lp = la + lb
    tag = f"{la}{lb}"
    hs = herm(lp)
    nh = len(hs)
    ca, cb = comps(la, order), comps(lb, order)
    na, nb = len(ca), len(cb)
    ind = " " * 6

    e_lines, e_names = e_block(la, lb, order, ind)
    decl_e = e_names

    # ---- ket item: density block -> Dt and ||rho_Q|| ----
    ket = []
    ket.append(f"   pure subroutine mdj_ket_item_{tag}(i, k0, kc, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, &")
    ket.append("                                    mp_bx, mp_by, mp_bz, mp_c, mp_sp, nsp, sp_ao1, sp_ao2, &")
    ket.append("                                    nao, dmat, ndt, dt, mp_nq)")
    ket.append("      !$acc routine seq")
    ket.append("      integer,  intent(in) :: i, k0, kc, nkc, hb, nhp, nsp, nao, ndt")
    ket.append("      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)")
    ket.append("      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)")
    ket.append("      integer,  intent(in) :: mp_sp(nhp), sp_ao1(nsp), sp_ao2(nsp)")
    ket.append("      real(dp), intent(in) :: dmat(nao, nao)")
    ket.append("      real(dp), intent(inout) :: dt(ndt), mp_nq(nhp)")
    ket.append("      integer  :: k, s, mu0, nu0, pos")
    ket.append("      real(dp) :: p, oo2p, pxa, pya, pza, pxb, pyb, pzb, offd, cs, nq2, ws")
    ket.append(decl("real(dp)", decl_e))
    dab = [f"dab{x + 1}_{y + 1}" for x in range(na) for y in range(nb)]
    ket.append(decl("real(dp)", dab))
    ket.append(decl("real(dp)", [f"g{h + 1}" for h in range(nh)]))
    ket.append(decl("real(dp)", [f"acc{s}" for s in range(lp + 1)]))
    ket.append("")
    ket.append("      k = k0 + i")
    ket.append("      s = mp_sp(k)")
    ket.append("      mu0 = sp_ao1(s) - 1")
    ket.append("      nu0 = sp_ao2(s) - 1")
    ket.append("      ! One shell pair stands for both orders off the diagonal: D_cd + D_dc.")
    ket.append("      offd = 1.0_dp")
    ket.append("      if (mu0 == nu0) offd = 0.0_dp")
    ket.append("      p = mp_z(k)")
    ket.append("      oo2p = 0.5_dp/p")
    ket.append("      pxa = mp_ax(k); pya = mp_ay(k); pza = mp_az(k)")
    ket.append("      pxb = mp_bx(k); pyb = mp_by(k); pzb = mp_bz(k)")
    ket.extend(e_lines)
    for x in range(na):
        for y in range(nb):
            ket.append(f"      dab{x + 1}_{y + 1} = dmat(mu0 + {x + 1}, nu0 + {y + 1}) "
                       f"+ offd*dmat(nu0 + {y + 1}, mu0 + {x + 1})")
    for hi, h in enumerate(hs):
        parts = []
        for x, a in enumerate(ca):
            for y, b in enumerate(cb):
                f = e3_factors(a, b, h)
                if f is None:
                    continue
                parts.append("*".join(f + [f"dab{x + 1}_{y + 1}"]))
        ket.extend(emit_sum(f"g{hi + 1}", parts, ind))
    ket.append("      ! Sign (-1)^|tau| and the pair coefficient (with K_cd) folded in.")
    ket.append("      cs = mp_c(k)")
    ket.append("      pos = kc + i")
    for hi, h in enumerate(hs):
        sg = "-" if sum(h) % 2 else ""
        ket.append(f"      g{hi + 1} = {sg}cs*g{hi + 1}")
        ket.append(f"      dt(hb + {hi}*nkc + pos) = g{hi + 1}")
    ket.extend(norm_lines(hs, "g", "p", ind))
    ket.append("      mp_nq(k) = sqrt(max(nq2, 0.0_dp))")
    ket.append(f"   end subroutine mdj_ket_item_{tag}")

    # ---- bra item: Jt -> J block ----
    bra = []
    bra.append(f"   pure subroutine mdj_bra_item_{tag}(i, s0, kl0, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, &")
    bra.append("                                    mp_bx, mp_by, mp_bz, mp_c, nsp, sp_ao1, sp_ao2, sp_k0, sp_nk, &")
    bra.append("                                    ndt, jt, nao, jmat)")
    bra.append("      !$acc routine seq")
    bra.append("      integer,  intent(in) :: i, s0, kl0, nkc, hb, nhp, nsp, ndt, nao")
    bra.append("      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)")
    bra.append("      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)")
    bra.append("      integer,  intent(in) :: sp_ao1(nsp), sp_ao2(nsp), sp_k0(nsp), sp_nk(nsp)")
    bra.append("      real(dp), intent(in) :: jt(ndt)")
    bra.append("      real(dp), intent(inout) :: jmat(nao, nao)")
    bra.append("      integer  :: s, mu0, nu0, kk, k, pos")
    bra.append("      real(dp) :: p, oo2p, pxa, pya, pza, pxb, pyb, pzb, cs")
    bra.append(decl("real(dp)", decl_e))
    bra.append(decl("real(dp)", [f"jh{h + 1}" for h in range(nh)]))
    jab = [f"jab{x + 1}_{y + 1}" for x in range(na) for y in range(nb)]
    bra.append(decl("real(dp)", jab))
    bra.append("")
    bra.append("      s = s0 + i")
    bra.append("      mu0 = sp_ao1(s) - 1")
    bra.append("      nu0 = sp_ao2(s) - 1")
    for n in jab:
        bra.append(f"      {n} = 0.0_dp")
    bra.append("      do kk = 1, sp_nk(s)")
    bra.append("         k = sp_k0(s) + kk")
    bra.append("         pos = k - kl0")
    bra.append("         p = mp_z(k)")
    bra.append("         oo2p = 0.5_dp/p")
    bra.append("         pxa = mp_ax(k); pya = mp_ay(k); pza = mp_az(k)")
    bra.append("         pxb = mp_bx(k); pyb = mp_by(k); pzb = mp_bz(k)")
    bra.append("         cs = mp_c(k)")
    for hi in range(nh):
        bra.append(f"         jh{hi + 1} = cs*jt(hb + {hi}*nkc + pos)")
    e_lines_in, _ = e_block(la, lb, order, " " * 9)
    bra.extend(e_lines_in)
    for x, a in enumerate(ca):
        for y, b in enumerate(cb):
            parts = []
            for hi, h in enumerate(hs):
                f = e3_factors(a, b, h)
                if f is None:
                    continue
                parts.append("*".join(f + [f"jh{hi + 1}"]))
            bra.extend(emit_sum(f"jab{x + 1}_{y + 1}", parts, " " * 9, accumulate=True))
    bra.append("      end do")
    bra.append("      ! Shell pairs are canonical and distinct, so no other thread writes")
    bra.append("      ! either block: no atomics.")
    for x in range(na):
        for y in range(nb):
            bra.append(f"      jmat(mu0 + {x + 1}, nu0 + {y + 1}) = jab{x + 1}_{y + 1}")
    bra.append("      if (mu0 /= nu0) then")
    for x in range(na):
        for y in range(nb):
            bra.append(f"         jmat(nu0 + {y + 1}, mu0 + {x + 1}) = jab{x + 1}_{y + 1}")
    bra.append("      end if")
    bra.append(f"   end subroutine mdj_bra_item_{tag}")

    # ---- bound item: I_P = |c| max_ab ||phi_a phi_b|| ----
    bnd = []
    bnd.append(f"   pure subroutine mdj_bound_item_{tag}(i, k0, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, &")
    bnd.append("                                      mp_c, mp_ip)")
    bnd.append("      !$acc routine seq")
    bnd.append("      integer,  intent(in) :: i, k0, nhp")
    bnd.append("      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)")
    bnd.append("      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)")
    bnd.append("      real(dp), intent(inout) :: mp_ip(nhp)")
    bnd.append("      integer  :: k")
    bnd.append("      real(dp) :: p, oo2p, pxa, pya, pza, pxb, pyb, pzb, nq2, ws, best")
    bnd.append(decl("real(dp)", decl_e))
    maxbox = 0
    boxes = []
    for a in ca:
        for b in cb:
            box = [h for h in hs if all(h[d] <= a[d] + b[d] for d in range(3))]
            boxes.append((a, b, box))
            maxbox = max(maxbox, len(box))
    bnd.append(decl("real(dp)", [f"g{h + 1}" for h in range(maxbox)]))
    bnd.append(decl("real(dp)", [f"acc{s}" for s in range(lp + 1)]))
    bnd.append("")
    bnd.append("      k = k0 + i")
    bnd.append("      p = mp_z(k)")
    bnd.append("      oo2p = 0.5_dp/p")
    bnd.append("      pxa = mp_ax(k); pya = mp_ay(k); pza = mp_az(k)")
    bnd.append("      pxb = mp_bx(k); pyb = mp_by(k); pzb = mp_bz(k)")
    bnd.extend(e_lines)
    bnd.append("      best = 0.0_dp")
    for a, b, box in boxes:
        bnd.append(f"      ! component pair {a} {b}")
        for hi, h in enumerate(box):
            f = e3_factors(a, b, h)
            bnd.append(f"      g{hi + 1} = " + ("*".join(f) if f else "1.0_dp"))
        bnd.extend(norm_lines(box, "g", "p", ind))
        bnd.append("      best = max(best, nq2)")
    bnd.append("      mp_ip(k) = abs(mp_c(k))*sqrt(best)")
    bnd.append(f"   end subroutine mdj_bound_item_{tag}")

    # ---- launches ----
    out = [header.format(what=f"Class ({la}{lb}) ket and bra stages of the MD J engine")]
    out.append(f"module trc_mdj_e{tag}")
    out.append("   use pic_types, only: dp")
    out.append("   implicit none")
    out.append("   private")
    out.append(f"   public :: mdj_ket_{tag}, mdj_bra_{tag}, mdj_bound_{tag}")
    out.append("")
    out.append("contains")
    out.append("")
    out.append(f"   subroutine mdj_ket_{tag}(n, k0, kc, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, &")
    out.append("                          mp_bx, mp_by, mp_bz, mp_c, mp_sp, nsp, sp_ao1, sp_ao2, &")
    out.append("                          nao, dmat, ndt, dt, mp_nq)")
    out.append(f"      !! Ket Hermite density of the {la}{lb} primitive pairs k0+1..k0+n, and")
    out.append("      !! each one's Coulomb norm. `kc` is k0's offset inside its lp block,")
    out.append("      !! whose Hermite rows start at `hb` with stride `nkc`.")
    out.append("      integer,  intent(in) :: n, k0, kc, nkc, hb, nhp, nsp, nao, ndt")
    out.append("      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)")
    out.append("      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)")
    out.append("      integer,  intent(in) :: mp_sp(nhp), sp_ao1(nsp), sp_ao2(nsp)")
    out.append("      real(dp), intent(in) :: dmat(nao, nao)")
    out.append("      real(dp), intent(inout) :: dt(ndt), mp_nq(nhp)")
    out.append("      integer :: i")
    out.append("      do concurrent(i = 1:n)")
    out.append(f"         call mdj_ket_item_{tag}(i, k0, kc, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, &")
    out.append("                              mp_bx, mp_by, mp_bz, mp_c, mp_sp, nsp, sp_ao1, sp_ao2, &")
    out.append("                              nao, dmat, ndt, dt, mp_nq)")
    out.append("      end do")
    out.append(f"   end subroutine mdj_ket_{tag}")
    out.append("")
    out.append(f"   subroutine mdj_bra_{tag}(n, s0, kl0, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, &")
    out.append("                          mp_bx, mp_by, mp_bz, mp_c, nsp, sp_ao1, sp_ao2, sp_k0, sp_nk, &")
    out.append("                          ndt, jt, nao, jmat)")
    out.append(f"      !! J blocks of the {la}{lb} shell pairs s0+1..s0+n from their primitive")
    out.append("      !! pairs' Jt; `kl0` is the first primitive of their lp block, less one.")
    out.append("      integer,  intent(in) :: n, s0, kl0, nkc, hb, nhp, nsp, ndt, nao")
    out.append("      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)")
    out.append("      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)")
    out.append("      integer,  intent(in) :: sp_ao1(nsp), sp_ao2(nsp), sp_k0(nsp), sp_nk(nsp)")
    out.append("      real(dp), intent(in) :: jt(ndt)")
    out.append("      real(dp), intent(inout) :: jmat(nao, nao)")
    out.append("      integer :: i")
    out.append("      do concurrent(i = 1:n)")
    out.append(f"         call mdj_bra_item_{tag}(i, s0, kl0, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, &")
    out.append("                              mp_bx, mp_by, mp_bz, mp_c, nsp, sp_ao1, sp_ao2, sp_k0, sp_nk, &")
    out.append("                              ndt, jt, nao, jmat)")
    out.append("      end do")
    out.append(f"   end subroutine mdj_bra_{tag}")
    out.append("")
    out.append(f"   subroutine mdj_bound_{tag}(n, k0, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, mp_ip)")
    out.append(f"      !! I_P = |c_P| max_ab ||phi_a phi_b||, the Coulomb self-norm of every")
    out.append(f"      !! component pair of the {la}{lb} primitive pairs k0+1..k0+n. Geometry only.")
    out.append("      integer,  intent(in) :: n, k0, nhp")
    out.append("      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)")
    out.append("      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)")
    out.append("      real(dp), intent(inout) :: mp_ip(nhp)")
    out.append("      integer :: i")
    out.append("      do concurrent(i = 1:n)")
    out.append(f"         call mdj_bound_item_{tag}(i, k0, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, &")
    out.append("                                mp_c, mp_ip)")
    out.append("      end do")
    out.append(f"   end subroutine mdj_bound_{tag}")
    out.append("")
    out.extend(ket)
    out.append("")
    out.extend(bra)
    out.append("")
    out.extend(bnd)
    out.append("")
    out.append(f"end module trc_mdj_e{tag}")
    return "\n".join(out) + "\n"


def norm_lines(hs, g, pexp, ind):
    """||rho||^2 = sum_{h,h'} (-1)^|h| g_h g_h' R^(0)_{h+h'}(PQ = 0), alpha = p/2.

    At PQ = 0 only all-even r survive, R^(0)_r = pref (r_x-1)!! (r_y-1)!!
    (r_z-1)!! (-p)^s / (2s + 1), s = |r|/2, pref = 2 pi^2.5/(p^2 sqrt(2p)):
    the recurrence with X_PQ = 0 halves each index at a factor of (r_d - 1).
    The sum is symmetric in h, h' (|h| and |h'| have the same parity), so it
    runs over h <= h' with the off-diagonal terms doubled.
    """
    lmax = max(sum(h) for h in hs)
    byS = {s: [] for s in range(lmax + 1)}
    for i, h in enumerate(hs):
        for j in range(i, len(hs)):
            h2 = hs[j]
            r = tuple(h[d] + h2[d] for d in range(3))
            if any(x % 2 for x in r):
                continue
            s = sum(r) // 2
            w = (1 if i == j else 2) * dfact(r[0]) * dfact(r[1]) * dfact(r[2])
            if sum(h) % 2:
                w = -w
            byS[s].append(f"{lit(float(w))}*{g}{i + 1}*{g}{j + 1}")
    lines = []
    for s in range(lmax + 1):
        lines.extend(emit_sum(f"acc{s}", byS[s], ind))
    lines.append(f"{ind}nq2 = acc0")
    lines.append(f"{ind}ws = 1.0_dp")
    for s in range(1, lmax + 1):
        lines.append(f"{ind}ws = -{pexp}*ws")
        lines.append(f"{ind}nq2 = nq2 + ws*acc{s}*{lit(1.0 / (2 * s + 1))}")
    lines.append(f"{ind}nq2 = nq2*{TWO_PI_2_5}/({pexp}*{pexp}*sqrt(2.0_dp*{pexp}))")
    return lines


def decl(typ, names, per=8):
    if not names:
        return ""
    rows = []
    for k in range(0, len(names), per):
        rows.append(f"      {typ} :: " + ", ".join(names[k:k + per]))
    return "\n".join(rows)


# ----------------------------------------------------------- Hermite stage --

def gen_r_module(lp, lq, nbatch, header, sieve=True):
    L = lp + lq
    tag = f"{lp}{lq}"
    hp, hq = herm(lp), herm(lq)
    iq = herm_index(lq)
    dags = rtree.batch_dags(L, nbatch, sieve=sieve)
    out = [header.format(what=f"Class [{lp}|{lq}] Hermite stage of the MD J engine, "
                         f"{len(dags)} R batch(es)")]
    out.append(f"module trc_mdj_r{tag}")
    out.append("   use pic_types, only: dp")
    out.append("   use trc_boys, only: BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV")
    out.append("   implicit none")
    out.append("   private")
    out.append(f"   public :: mdj_r_{tag}")
    out.append("")
    out.append(f"   real(dp), parameter :: TWO_PI_2_5 = {TWO_PI_2_5}")
    out.append("")
    out.append("contains")
    out.append("")
    out.append('#include "inc/trc_boys_eval.inc"')
    out.append("")
    args = ("np, kp0, hbp, nq, kq0, hbq, qlen, c0, cstep, nhp, mp_z, mp_px, mp_py, mp_pz, mp_ip, mp_nq, "
            "ndt, dt, jt, cut")
    out.append(f"   subroutine mdj_r_{tag}(np, kp0, hbp, nq, kq0, hbq, qlen, c0, cstep, nhp, mp_z, mp_px, &")
    out.append("                        mp_py, mp_pz, mp_ip, mp_nq, ndt, dt, jt, cut)")
    out.append(f"      !! Jt[P] += [P|Q] Dt[Q] for the lp = {lp} pairs kp0+1..kp0+np against the")
    out.append(f"      !! lq = {lq} pairs kq0+1..kq0+nq, one launch per R batch. hbp and hbq are")
    out.append("      !! the two blocks' first Hermite slot, less one; each row has the block's")
    out.append("      !! own pair count as stride. One thread per (P, chunk of `qlen` kets):")
    out.append("      !! one per P alone left a heavy class with a few hundred threads, each")
    out.append("      !! walking the whole ket list (measured: 128 threads, 27 ms, on w20).")
    out.append("      !! Only chunks c0, c0 + cstep, ... (0-based): an MPI rank's share, its")
    out.append("      !! partial Jt summed over ranks by the caller; (0, 1) is every chunk.")
    out.append("      integer,  intent(in) :: np, kp0, hbp, nq, kq0, hbq, qlen, c0, cstep, nhp, ndt")
    out.append("      real(dp), intent(in) :: mp_z(nhp), mp_px(nhp), mp_py(nhp), mp_pz(nhp)")
    out.append("      real(dp), intent(in) :: mp_ip(nhp), mp_nq(nhp), dt(ndt), cut")
    out.append("      real(dp), intent(inout) :: jt(ndt)")
    out.append("      integer :: i, nchunk, nmine")
    out.append("      nchunk = (nq + qlen - 1)/qlen")
    out.append("      nmine = max(0, (nchunk - c0 + cstep - 1)/cstep)")
    for b in range(len(dags)):
        out.append("      do concurrent(i = 1:np*nmine)")
        out.extend(wrap_args(f"call mdj_r_item_{tag}_b{b + 1}", "i, " + args, " " * 9))
        out.append("      end do")
    out.append(f"   end subroutine mdj_r_{tag}")
    for b, g in enumerate(dags):
        out.append("")
        out.extend(r_item(lp, lq, b + 1, g, hp, hq, iq, args))
    out.append("")
    out.append(f"end module trc_mdj_r{tag}")
    return "\n".join(out) + "\n"


def r_item(lp, lq, b, g, hp, hq, iq, args):
    tag = f"{lp}{lq}"
    stmts = g.statements()
    leaves = g.leaves()
    mb = max(leaves)
    inter, flops = g.counts()
    # Which Jt slots and Dt slots this batch touches.
    contr = {}
    for (r, _) in g.targets:
        lst = []
        for i, t in enumerate(hp):
            tau = tuple(r[d] - t[d] for d in range(3))
            if min(tau) < 0 or sum(tau) > lq:
                continue
            lst.append((i, iq[tau]))
        contr[r] = lst
    jused = sorted({i for lst in contr.values() for i, _ in lst})
    dused = sorted({j for lst in contr.values() for _, j in lst})
    rnames = [s.lhs for s in stmts]
    ind = " " * 9
    o = []
    o.extend(wrap_args(f"pure subroutine mdj_r_item_{tag}_b{b}", "i, " + args, " " * 3))
    o.append(f"      !! Batch {b}: {len(g.targets)} targets, {inter} intermediates, "
             f"{flops} R flops per primitive quartet.")
    o.append("      !$acc routine seq")
    o.append("      integer,  intent(in) :: i, np, kp0, hbp, nq, kq0, hbq, qlen, c0, cstep, nhp, ndt")
    o.append("      real(dp), intent(in) :: mp_z(nhp), mp_px(nhp), mp_py(nhp), mp_pz(nhp)")
    o.append("      real(dp), intent(in) :: mp_ip(nhp), mp_nq(nhp), dt(ndt), cut")
    o.append("      real(dp), intent(inout) :: jt(ndt)")
    o.append("      integer  :: kp, kq, k, ip1, k1, k2")
    o.append("      real(dp) :: p, px, py, pz, ip, q, zpe, alpha, pref, tval, m2a, w")
    o.append("      real(dp) :: pqx, pqy, pqz")
    o.append(f"      real(dp) :: f(0:BOYS_MMAX)")
    o.append(decl("real(dp)", [rtree.leaf_name(m) for m in leaves]))
    o.append(decl("real(dp)", rnames))
    o.append(decl("real(dp)", [f"j{i + 1}" for i in jused]))
    o.append(decl("real(dp)", [f"d{j + 1}" for j in dused]))
    o.append("")
    o.append("      ! P fastest, so the threads of a warp share one ket chunk.")
    o.append("      ip1 = mod(i - 1, np) + 1")
    o.append("      k1 = (c0 + ((i - 1)/np)*cstep)*qlen + 1")
    o.append("      k2 = min(nq, k1 + qlen - 1)")
    o.append("      kp = kp0 + ip1")
    o.append("      p = mp_z(kp)")
    o.append("      px = mp_px(kp); py = mp_py(kp); pz = mp_pz(kp)")
    o.append("      ip = mp_ip(kp)")
    for i in jused:
        o.append(f"      j{i + 1} = 0.0_dp")
    o.append("      do k = k1, k2")
    o.append("         kq = kq0 + k")
    o.append("         ! I_P ||rho_Q|| bounds this pair's share of any J element.")
    o.append("         if (ip*mp_nq(kq) <= cut) cycle")
    o.append("         q = mp_z(kq)")
    o.append("         zpe = p + q")
    o.append("         alpha = p*q/zpe")
    o.append("         pref = TWO_PI_2_5/(p*q*sqrt(zpe))")
    o.append("         pqx = px - mp_px(kq); pqy = py - mp_py(kq); pqz = pz - mp_pz(kq)")
    o.append("         tval = alpha*(pqx*pqx + pqy*pqy + pqz*pqz)")
    o.append(f"         call boys_eval({mb}, tval, f)")
    o.append("         m2a = -2.0_dp*alpha")
    o.append("         w = pref")
    for m in range(mb + 1):
        if m > 0:
            o.append("         w = w*m2a")
        if m in leaves:
            o.append(f"         {rtree.leaf_name(m)} = w*f({m})")
    for j in dused:
        o.append(f"         d{j + 1} = dt(hbq + {j}*nq + k)")

    def contract(r):
        val = rtree.leaf_name(0) if r == (0, 0, 0) else rtree.r_name(r, 0)
        for i, j in contr[r]:
            o.append(f"{ind}j{i + 1} = j{i + 1} + {val}*d{j + 1}")

    if ((0, 0, 0), 0) in g.tset:
        contract((0, 0, 0))
    for s in stmts:
        o.extend(emit_sum(s.lhs, expr_of(s.terms), ind))
        name = s.lhs
        # level-0 statements are targets: contract and let them die at once.
        fields = name[1:].split("_")
        if fields[3] == "0":
            r = tuple(int(x) for x in fields[:3])
            if (r, 0) in g.tset:
                contract(r)
    o.append("      end do")
    o.append("      ! Several chunks share one P: atomics, once per Hermite slot per thread.")
    for i in jused:
        o.append("      !$acc atomic update")
        o.append(f"      jt(hbp + {i}*np + ip1) = jt(hbp + {i}*np + ip1) + j{i + 1}")
    o.append(f"   end subroutine mdj_r_item_{tag}_b{b}")
    return o


# ------------------------------------------------------------- dispatcher --

def gen_dispatch(lmax, header):
    o = [header.format(what="Host dispatch of the MD J engine kernels")]
    o.append(f"#if TRC_LMAX != {lmax}")
    o.append('#error "trc_mdj_kernels was generated for a different TRC_LMAX. Regenerate it."')
    o.append("#endif")
    o.append("module trc_mdj_kernels")
    o.append("   use pic_types, only: dp")
    for la in range(lmax + 1):
        for lb in range(la + 1):
            o.append(f"   use trc_mdj_e{la}{lb}, only: mdj_ket_{la}{lb}, mdj_bra_{la}{lb}, mdj_bound_{la}{lb}")
    for lp in range(2 * lmax + 1):
        for lq in range(2 * lmax + 1):
            o.append(f"   use trc_mdj_r{lp}{lq}, only: mdj_r_{lp}{lq}")
    o.append("   implicit none")
    o.append("   private")
    o.append("   public :: mdj_ket, mdj_bra, mdj_bound, mdj_hermite, MDJ_LPMAX")
    o.append("")
    o.append(f"   integer, parameter :: MDJ_LPMAX = {2 * lmax}")
    o.append("   !! Largest Hermite degree of one shell pair these kernels cover.")
    o.append("")
    o.append("contains")
    o.append("")
    ket_args = ("n, k0, kc, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, mp_sp, "
                "nsp, sp_ao1, sp_ao2, nao, dmat, ndt, dt, mp_nq")
    bra_args = ("n, s0, kl0, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, "
                "nsp, sp_ao1, sp_ao2, sp_k0, sp_nk, ndt, jt, nao, jmat")
    r_args = ("np, kp0, hbp, nq, kq0, hbq, qlen, c0, cstep, nhp, mp_z, mp_px, mp_py, mp_pz, mp_ip, mp_nq, "
              "ndt, dt, jt, cut")

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

    o.append("   subroutine mdj_ket(la, lb, " + "n, k0, kc, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, &")
    o.append("                      mp_bx, mp_by, mp_bz, mp_c, mp_sp, nsp, sp_ao1, sp_ao2, nao, dmat, ndt, dt, mp_nq)")
    o.append("      integer,  intent(in) :: la, lb, n, k0, kc, nkc, hb, nhp, nsp, nao, ndt")
    o.append("      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)")
    o.append("      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)")
    o.append("      integer,  intent(in) :: mp_sp(nhp), sp_ao1(nsp), sp_ao2(nsp)")
    o.append("      real(dp), intent(in) :: dmat(nao, nao)")
    o.append("      real(dp), intent(inout) :: dt(ndt), mp_nq(nhp)")
    o.append(f"      select case (la*{lmax + 1} + lb)")
    for la in range(lmax + 1):
        for lb in range(la + 1):
            o.append(f"      case ({la * (lmax + 1) + lb})")
            o.extend(wrap(f"call mdj_ket_{la}{lb}({ket_args})"))
    o.append("      case default")
    o.append('         error stop "trc_mdj_kernels: no ket kernel for this shell-pair class"')
    o.append("      end select")
    o.append("   end subroutine mdj_ket")
    o.append("")
    o.append("   subroutine mdj_bra(la, lb, n, s0, kl0, nkc, hb, nhp, mp_z, mp_ax, mp_ay, mp_az, &")
    o.append("                      mp_bx, mp_by, mp_bz, mp_c, nsp, sp_ao1, sp_ao2, sp_k0, sp_nk, ndt, jt, nao, jmat)")
    o.append("      integer,  intent(in) :: la, lb, n, s0, kl0, nkc, hb, nhp, nsp, ndt, nao")
    o.append("      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)")
    o.append("      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)")
    o.append("      integer,  intent(in) :: sp_ao1(nsp), sp_ao2(nsp), sp_k0(nsp), sp_nk(nsp)")
    o.append("      real(dp), intent(in) :: jt(ndt)")
    o.append("      real(dp), intent(inout) :: jmat(nao, nao)")
    o.append(f"      select case (la*{lmax + 1} + lb)")
    for la in range(lmax + 1):
        for lb in range(la + 1):
            o.append(f"      case ({la * (lmax + 1) + lb})")
            o.extend(wrap(f"call mdj_bra_{la}{lb}({bra_args})"))
    o.append("      case default")
    o.append('         error stop "trc_mdj_kernels: no bra kernel for this shell-pair class"')
    o.append("      end select")
    o.append("   end subroutine mdj_bra")
    o.append("")
    o.append("   subroutine mdj_bound(la, lb, n, k0, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, &")
    o.append("                        mp_c, mp_ip)")
    o.append("      integer,  intent(in) :: la, lb, n, k0, nhp")
    o.append("      real(dp), intent(in) :: mp_z(nhp), mp_ax(nhp), mp_ay(nhp), mp_az(nhp)")
    o.append("      real(dp), intent(in) :: mp_bx(nhp), mp_by(nhp), mp_bz(nhp), mp_c(nhp)")
    o.append("      real(dp), intent(inout) :: mp_ip(nhp)")
    o.append(f"      select case (la*{lmax + 1} + lb)")
    for la in range(lmax + 1):
        for lb in range(la + 1):
            o.append(f"      case ({la * (lmax + 1) + lb})")
            o.extend(wrap(f"call mdj_bound_{la}{lb}(n, k0, nhp, mp_z, mp_ax, mp_ay, mp_az, mp_bx, mp_by, mp_bz, mp_c, mp_ip)"))
    o.append("      case default")
    o.append('         error stop "trc_mdj_kernels: no bound kernel for this shell-pair class"')
    o.append("      end select")
    o.append("   end subroutine mdj_bound")
    o.append("")
    o.append("   subroutine mdj_hermite(lp, lq, np, kp0, hbp, nq, kq0, hbq, qlen, c0, cstep, nhp, mp_z, mp_px, &")
    o.append("                          mp_py, mp_pz, mp_ip, mp_nq, ndt, dt, jt, cut)")
    o.append("      integer,  intent(in) :: lp, lq, np, kp0, hbp, nq, kq0, hbq, qlen, c0, cstep, nhp, ndt")
    o.append("      real(dp), intent(in) :: mp_z(nhp), mp_px(nhp), mp_py(nhp), mp_pz(nhp)")
    o.append("      real(dp), intent(in) :: mp_ip(nhp), mp_nq(nhp), dt(ndt), cut")
    o.append("      real(dp), intent(inout) :: jt(ndt)")
    n = 2 * lmax + 1
    o.append(f"      select case (lp*{n} + lq)")
    for lp in range(n):
        for lq in range(n):
            o.append(f"      case ({lp * n + lq})")
            o.extend(wrap(f"call mdj_r_{lp}{lq}({r_args})"))
    o.append("      case default")
    o.append('         error stop "trc_mdj_kernels: no Hermite kernel for this class"')
    o.append("      end select")
    o.append("   end subroutine mdj_hermite")
    o.append("")
    o.append("end module trc_mdj_kernels")
    return "\n".join(o) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--lmax", type=int, required=True)
    ap.add_argument("--order", default="lexicographic", choices=sorted(ORDER_CODE))
    ap.add_argument("--batches", default="",
                    help="comma list L:n overriding the batch count per total L")
    ap.add_argument("--no-sieve", action="store_true")
    ap.add_argument("-o", "--outdir", required=True)
    a = ap.parse_args()

    table = {}
    for item in filter(None, a.batches.split(",")):
        L, n = item.split(":")
        table[int(L)] = int(n)

    cmd = "scripts/gen_md_j.py --lmax {}".format(a.lmax)
    if a.order not in ("lexicographic", "libcint"):
        cmd += f" --order {a.order}"
    if a.batches:
        cmd += f" --batches {a.batches}"
    if a.no_sieve:
        cmd += " --no-sieve"
    code = ORDER_CODE[a.order]
    header = ("!\n! {what}.\n!\n"
              f"! GENERATED by {cmd} -- do not edit.\n"
              "! See scripts/gen_md_j.py for the algorithm and md/reference.py for its check.\n!\n")
    e_header = header + (f"#if TRC_CART_ORDER != {code}\n"
                         '#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."\n'
                         "#endif\n")

    os.makedirs(a.outdir, exist_ok=True)

    def write(name, text):
        with open(os.path.join(a.outdir, name), "w") as fh:
            fh.write(text)

    for la in range(a.lmax + 1):
        for lb in range(la + 1):
            write(f"trc_mdj_e{la}{lb}.F90", gen_e_module_loops(la, lb, a.order, e_header))
    for lp in range(2 * a.lmax + 1):
        for lq in range(2 * a.lmax + 1):
            L = lp + lq
            write(f"trc_mdj_r{lp}{lq}.F90",
                  gen_r_module(lp, lq, nbatch_for(L, table), header, sieve=not a.no_sieve))
    write("trc_mdj_kernels.F90", gen_dispatch(a.lmax, header))
    return 0


if __name__ == "__main__":
    sys.exit(main())
