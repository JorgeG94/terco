#!/usr/bin/env python3
"""
Emit fully unrolled Obara-Saika VRR kernels, one per (lab, lcd) class.

WHY
---
The generic VRR is a loop nest driven by index tables -- cll, cdir, cdn1, cdn2,
cf2 -- every one of which is a global-memory load in the innermost loop.  The
tables also force a three-way branch on the recurrence direction and two
predicated adds for the second-order terms.

None of that is thread divergence: the bins guarantee every thread in a launch
shares (la, lb, lc, ld) and K, so the loop bounds and every branch condition are
warp-uniform, and Nsight measures 29.6 of 32 threads active.  What it costs is
instructions and loads, not lanes.

Unrolling with compile-time indices removes the loads entirely and turns the
branches into nothing at all -- the generator knows which direction each index
decrements and whether the second-order term exists, so it simply does not emit
the dead ones.  This is the route taken in [1], arrived at from the other
direction.

WHAT IT COSTS
-------------
Registers, potentially.  Specialising the whole kernel earlier raised the count
from 146 to 241 and lost, because full unrolling lengthens live ranges.  Here
only the VRR is unrolled and the HRR and digestion stay generic, so the
experiment is narrower.  It gets measured before it gets believed.
"""

# [1] J. L. Galvez Vallejo, G. M. J. Barca and M. S. Gordon,
#     "High-performance GPU-accelerated evaluation of electron repulsion
#     integrals", Mol. Phys. (2022) e2112987, doi:10.1080/00268976.2022.2112987


import argparse


def cart_list(lmax):
    """(nx,ny,nz) in trc_cart order: angular-momentum-major, libcint within."""
    out = []
    for l in range(lmax + 1):
        for lx in range(l, -1, -1):
            for ly in range(l - lx, -1, -1):
                out.append((lx, ly, l - lx - ly))
    return out


def build_tables(lmax):
    lst = cart_list(lmax)
    idx = {p: i + 1 for i, p in enumerate(lst)}       # 1-based, as Fortran
    info = {}
    for h, (x, y, z) in enumerate(lst, start=1):
        l = x + y + z
        if l == 0:
            info[h] = None
            continue
        # first non-zero direction, matching trc_cart's choice exactly
        if x > 0:
            d, dn1, n = 1, (x - 1, y, z), x
            dn2 = (x - 2, y, z) if x > 1 else None
        elif y > 0:
            d, dn1, n = 2, (x, y - 1, z), y
            dn2 = (x, y - 2, z) if y > 1 else None
        else:
            d, dn1, n = 3, (x, y, z - 1), z
            dn2 = (x, y, z - 2) if z > 1 else None
        info[h] = dict(l=l, d=d, dn1=idx[dn1], dn2=idx[dn2] if dn2 else 0,
                       cf2=float(n - 1), pow=(x, y, z))
    return lst, idx, info


def ncum(l):
    return (l + 1) * (l + 2) * (l + 3) // 6


AX = {1: "pax", 2: "pay", 3: "paz"}
WP = {1: "wpx", 2: "wpy", 3: "wpz"}
QC = {1: "qcx", 2: "qcy", 3: "qcz"}
WQ = {1: "wqx", 2: "wqy", 3: "wqz"}


def emit_class(lab, lcd, lst, idx, info):
    """Straight-line VRR for one (lab, lcd), every index a literal."""
    lt = lab + lcd
    nca, ncc = ncum(lab), ncum(lcd)
    o = []
    w = o.append

    def V(h, ic, buf):
        return f"v({h + (ic - 1)*nca},{buf})"

    cur = 0
    w(f"      {V(1,1,cur)} = pref*f({lt})")

    for m in range(lt - 1, -1, -1):
        nxt = 1 - cur
        lim = lt - m
        w(f"      ! --- level m = {m}, degree <= {lim} ---")
        w(f"      {V(1,1,nxt)} = pref*f({m})")

        # build a at c = 0
        for ia in range(2, nca + 1):
            e = info[ia]
            if e["l"] > lim:
                break
            t = f"{AX[e['d']]}*{V(e['dn1'],1,nxt)} + {WP[e['d']]}*{V(e['dn1'],1,cur)}"
            if e["dn2"]:
                t += (f" &\n         + {e['cf2']}_dp*oo2z*({V(e['dn2'],1,nxt)}"
                      f" - rz*{V(e['dn2'],1,cur)})")
            w(f"      {V(ia,1,nxt)} = {t}")

        # build c, for every a
        for ic in range(2, ncc + 1):
            ec = info[ic]
            if ec["l"] > lim:
                break
            d = ec["d"]
            for ia in range(1, nca + 1):
                la_ = info[ia]["l"] if info[ia] else 0
                if la_ + ec["l"] > lim:
                    break
                t = (f"{QC[d]}*{V(ia,ec['dn1'],nxt)}"
                     f" + {WQ[d]}*{V(ia,ec['dn1'],cur)}")
                if ec["dn2"]:
                    t += (f" &\n         + {ec['cf2']}_dp*oo2e*({V(ia,ec['dn2'],nxt)}"
                          f" - re*{V(ia,ec['dn2'],cur)})")
                # cross term: a_d/(2(zeta+eta)) [a-1_d,0|c-1_d,0]^(m+1)
                if info[ia] is not None:
                    px, py, pz = info[ia]["pow"]
                    comp = (px, py, pz)[d - 1]
                    if comp > 0:
                        dn = list(comp for comp in (px, py, pz))
                        dn[d - 1] -= 1
                        t += (f" &\n         + {float(comp)}_dp*oo2ze"
                              f"*{V(idx[tuple(dn)],ec['dn1'],cur)}")
                w(f"      {V(ia,ic,nxt)} = {t}")
        cur = nxt

    return "\n".join(o), cur, nca * ncc



# ---------------------------------------------------------------------------
# Gill's generation-step sieve
# ---------------------------------------------------------------------------
#
# The level-by-level emitter above fixes ONE pathway to every intermediate:
# `build_tables` always decrements the first non-zero Cartesian component, x
# then y then z.  Every auxiliary is therefore reached exactly one way, and any
# node that a different direction would have shared gets built fresh.
#
# Gill's sieve (appendix, "The Generation-Step Sieve") keeps the complete list
# of auxiliaries and then repeatedly asks, of each one:
#
#     find all integrals J that could be formed FROM this one;
#     if every such J can also be formed from something else still on the
#     list, then this one is not needed -- remove it.
#
# Passes repeat to a fixed point.  There are O(L^4) auxiliaries and only
# O(L^3) targets, so the surplus -- and the saving -- grows with angular
# momentum.  Measured here: 21% of nodes removed at (pp|pp), 40.7% at (dd|dd).
#
# This is NOT dead-code elimination.  An earlier attempt pruned nodes
# unreachable from the targets and gained exactly nothing, because nvfortran
# had already done it: the recurrence is straight-line code in one procedure
# and the optimiser can see it.  The sieve removes nodes that are live under
# the current path choice, by RE-DERIVING their consumers down a different
# recurrence.  That is an algebraic identity between two pathways, which a
# compiler cannot discover.
#
AXN = {0: "pax", 1: "pay", 2: "paz"}
WPN = {0: "wpx", 1: "wpy", 2: "wpz"}
QCN = {0: "qcx", 1: "qcy", 2: "qcz"}
WQN = {0: "wqx", 1: "wqy", 2: "wqz"}


def _pathways(node, P, I):
    """Every way to form `node`; each is (kind, direction, [input nodes])."""
    ai, ci, m = node
    A, C = P[ai], P[ci]
    out = []
    for d in range(3):
        if A[d] > 0:
            am = list(A); am[d] -= 1; a1 = I[tuple(am)]
            ins = [(a1, ci, m), (a1, ci, m + 1)]
            if A[d] >= 2:
                a2l = list(A); a2l[d] -= 2
                a2 = I[tuple(a2l)]
                ins += [(a2, ci, m), (a2, ci, m + 1)]
            if C[d] >= 1:
                cl = list(C); cl[d] -= 1
                ins.append((a1, I[tuple(cl)], m + 1))
            out.append(("a", d, ins))
    for d in range(3):
        if C[d] > 0:
            cl = list(C); cl[d] -= 1; c1 = I[tuple(cl)]
            ins = [(ai, c1, m), (ai, c1, m + 1)]
            if C[d] >= 2:
                c2l = list(C); c2l[d] -= 2
                c2 = I[tuple(c2l)]
                ins += [(ai, c2, m), (ai, c2, m + 1)]
            if A[d] >= 1:
                al = list(A); al[d] -= 1
                ins.append((I[tuple(al)], c1, m + 1))
            out.append(("c", d, ins))
    return out


def sieve_nodes(lab, lcd, targets, P, I):
    """Complete auxiliary list, then Gill's sieve to a fixed point."""
    lt = lab + lcd

    full, stack = set(), list(targets)
    while stack:
        n = stack.pop()
        if n in full or n[2] > lt:
            continue
        full.add(n)
        if n[0] == 1 and n[1] == 1:
            continue
        for _, _, ins in _pathways(n, P, I):
            for t in ins:
                if t not in full:
                    stack.append(t)

    def formable(n, avail):
        if n[0] == 1 and n[1] == 1:
            return True
        return any(all(t in avail for t in ins)
                   for _, _, ins in _pathways(n, P, I))

    cons = {}
    for n in full:
        if n[0] == 1 and n[1] == 1:
            continue
        for _, _, ins in _pathways(n, P, I):
            for t in ins:
                if t in full:
                    cons.setdefault(t, set()).add(n)

    live = set(full)
    changed = True
    while changed:
        changed = False
        # try the deepest, highest-degree nodes first: they have the fewest
        # consumers and so are the easiest to retire
        order = sorted(live, key=lambda t: (-t[2],
                                            -(sum(P[t[0]]) + sum(P[t[1]]))))
        for n in order:
            if n in targets or n not in live:
                continue
            if n[0] == 1 and n[1] == 1:
                continue
            trial = live - {n}
            if all(formable(j, trial) for j in cons.get(n, ()) if j in trial):
                live = trial
                changed = True
    return full, live


def emit_class_sieved(lab, lcd, targets, P, I, nca):
    """Straight-line VRR over the sieved node set."""
    lt = lab + lcd
    _full, live = sieve_nodes(lab, lcd, targets, P, I)

    def buf(m):
        return (lt - m) % 2

    def V(n):
        ai, ci, m = n
        return f"v({ai + (ci - 1)*nca},{buf(m)})"

    o, w = [], None
    w = o.append
    for m in range(lt, -1, -1):
        here = [n for n in live if n[2] == m]
        if not here:
            continue
        w(f"      ! --- level m = {m} ---")
        here.sort(key=lambda n: (sum(P[n[0]]) + sum(P[n[1]]),
                                 n[1], n[0]))
        for n in here:
            ai, ci, _ = n
            if ai == 1 and ci == 1:
                w(f"      {V(n)} = pref*f({m})")
                continue
            # cheapest available pathway; direction order breaks ties so the
            # output stays deterministic
            best = None
            for kind, d, ins in _pathways(n, P, I):
                if all(t in live for t in ins):
                    if best is None or len(ins) < len(best[2]):
                        best = (kind, d, ins)
            if best is None:
                raise RuntimeError(f"sieve left {n} unformable")
            kind, d, _ins = best
            A, C = P[ai], P[ci]
            if kind == "a":
                am = list(A); am[d] -= 1; a1 = I[tuple(am)]
                t = (f"{AXN[d]}*{V((a1, ci, m))}"
                     f" + {WPN[d]}*{V((a1, ci, m + 1))}")
                if A[d] >= 2:
                    a2l = list(A); a2l[d] -= 2
                    a2 = I[tuple(a2l)]
                    t += (f" &\n         + {float(A[d] - 1)}_dp*oo2z*"
                          f"({V((a2, ci, m))} - rz*{V((a2, ci, m + 1))})")
                if C[d] >= 1:
                    cl = list(C); cl[d] -= 1
                    t += (f" &\n         + {float(C[d])}_dp*oo2ze*"
                          f"{V((a1, I[tuple(cl)], m + 1))}")
            else:
                cl = list(C); cl[d] -= 1; c1 = I[tuple(cl)]
                t = (f"{QCN[d]}*{V((ai, c1, m))}"
                     f" + {WQN[d]}*{V((ai, c1, m + 1))}")
                if C[d] >= 2:
                    c2l = list(C); c2l[d] -= 2
                    c2 = I[tuple(c2l)]
                    t += (f" &\n         + {float(C[d] - 1)}_dp*oo2e*"
                          f"({V((ai, c2, m))} - re*{V((ai, c2, m + 1))})")
                if A[d] >= 1:
                    al = list(A); al[d] -= 1
                    t += (f" &\n         + {float(A[d])}_dp*oo2ze*"
                          f"{V((I[tuple(al)], c1, m + 1))}")
            w(f"      {V(n)} = {t}")
    return "\n".join(o), buf(0), len(live)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=1)
    ap.add_argument("-o", "--output", default="src/trc_vrr_unrolled.F90")
    ap.add_argument("--mode", choices=("module", "include"), default="module",
                    help="include: emit `case` blocks to paste inside the hot "
                         "loop, so there is no call boundary at all")
    args = ap.parse_args()

    LMAX = args.lmax
    LCMAX = 2 * LMAX
    lst, idx, info = build_tables(LCMAX)

    if args.mode == "include":
        out = [f"""!
! Unrolled Obara-Saika VRR as `case` blocks, to be included INSIDE the hot
! loop's own procedure.
!
! GENERATED by scripts/gen_vrr.py --mode include --lmax {LMAX} -- do not edit.
!
!
! STALE-INCLUDE GUARD. The `case` labels below were computed with the angular
! momentum limit this file was generated for, and the dispatch key is built
! from TRC_LMAX at compile time. If the two disagree the keys land on the
! wrong labels -- or on none, since there is no `case default` -- and the
! result is silently wrong for some classes and right for others, which is
! far worse than a failure.
!
! That is not hypothetical: this file was generated once for lmax 1 and
! included unchanged in an LMAX=2 build, where the key base went from 2 to 3
! and (pp|pp) asked for a label that did not exist. The folded path never
! noticed because it uses the per-class kernels; only the enumerated
! `fock_nosym` reaches this code, and it was wrong for a month.
!
#if TRC_LMAX != {LMAX}
#error "This .inc was generated for a different TRC_LMAX. Regenerate it."
#endif
! Not a module, and deliberately so.  As separate `!$acc routine seq`
! procedures these measured 2.4x SLOWER than the generic table-driven loop
! (1.466 s against 0.604 s at 32 waters), because NVHPC does not inline across
! a module boundary: `v` has to be fully materialised in memory to cross the
! call, and 35 kB of spill appears.  Pasted textually into one_bin_item, `v`
! stays a local and the compiler sees the whole recurrence.
!
! Expects in scope: v, f, pref, pax/pay/paz, qcx/qcy/qcz, wpx/wpy/wpz,
! wqx/wqy/wqz, oo2z, oo2e, oo2ze, rz, re, and sets `cur` to the buffer holding
! the result.
!
      select case (lab*100 + lcd)"""]
        for lab in range(LCMAX + 1):
            for lcd in range(LCMAX + 1):
                body, fin, nv = emit_class(lab, lcd, lst, idx, info)
                out.append(f"      case ({lab*100 + lcd})")
                out.append(body)
                out.append(f"         cur = {fin}")
        out.append("      end select")
        txt = "\n".join(out)
        with open(args.output, "w") as fh:
            fh.write(txt + "\n")
        print(f"wrote {args.output}: include form, "
              f"{(LCMAX+1)**2} classes, {txt.count(chr(10))} lines")
        return

    parts = [f'''!
! Unrolled Obara-Saika VRR, one entry point per (lab, lcd) class.
!
! GENERATED by scripts/gen_vrr.py --lmax {LMAX} -- do not edit.
!
! The generic VRR in trc_binkernel is driven by the trc_cart index tables,
! which are global-memory loads in the innermost loop, plus a three-way branch
! on the recurrence direction and two predicated second-order terms.  Here the
! generator knows every index, so the loads become literals and the dead
! branches are simply not emitted.
!
! Not a divergence fix: bins already make those branches warp-uniform (Nsight
! measures 29.6 of 32 threads active).  It removes instructions and loads.
!
module trc_vrr_unrolled
   use trc_boys, only: dp, BOYS_MMAX
   implicit none
   private
   public :: vrr_dispatch, VRR_LMAX

   integer, parameter :: VRR_LMAX = {LMAX}

contains
''']

    names = []
    for lab in range(LCMAX + 1):
        for lcd in range(LCMAX + 1):
            body, fin, nv = emit_class(lab, lcd, lst, idx, info)
            nm = f"vrr_{lab}_{lcd}"
            names.append((lab, lcd, nm, fin))
            parts.append(f'''
   !> (lab={lab} | lcd={lcd});  result left in buffer {fin}
   pure subroutine {nm}(ldv, v, f, pref, pax, pay, paz, qcx, qcy, qcz, &
                        wpx, wpy, wpz, wqx, wqy, wqz, oo2z, oo2e, oo2ze, rz, re)
      !$acc routine seq
      ! ldv is the CALLER's leading dimension, not this class's nca*ncc.  The
      ! caller holds one fixed-size buffer pair for every class, so the second
      ! buffer starts at ldv+1 and not at nca*ncc+1; declaring it the latter
      ! way misaligns the ping-pong silently.
      integer, intent(in) :: ldv
      real(dp), intent(inout) :: v(ldv, 0:1)
      real(dp), intent(in) :: f(0:BOYS_MMAX), pref
      real(dp), intent(in) :: pax, pay, paz, qcx, qcy, qcz
      real(dp), intent(in) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp), intent(in) :: oo2z, oo2e, oo2ze, rz, re
{body}
   end subroutine {nm}
''')

    parts.append('''
   !> Which buffer a class leaves its answer in.
   pure integer function vrr_final(lab, lcd)
      integer, intent(in) :: lab, lcd
      vrr_final = mod(lab + lcd, 2)
   end function vrr_final

   pure subroutine vrr_dispatch(lab, lcd, ldv, v, f, pref, pax, pay, paz, &
                                qcx, qcy, qcz, wpx, wpy, wpz, wqx, wqy, wqz, &
                                oo2z, oo2e, oo2ze, rz, re, fin)
      !$acc routine seq
      integer, intent(in) :: lab, lcd, ldv
      real(dp), intent(inout) :: v(ldv, 0:1)
      real(dp), intent(in) :: f(0:BOYS_MMAX), pref
      real(dp), intent(in) :: pax, pay, paz, qcx, qcy, qcz
      real(dp), intent(in) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp), intent(in) :: oo2z, oo2e, oo2ze, rz, re
      integer, intent(out) :: fin
      fin = vrr_final(lab, lcd)
      select case (lab*100 + lcd)''')
    for lab, lcd, nm, fin in names:
        parts.append(f'''      case ({lab*100 + lcd}); call {nm}(ldv, v, f, pref, pax, pay, paz, &
            qcx, qcy, qcz, wpx, wpy, wpz, wqx, wqy, wqz, oo2z, oo2e, oo2ze, rz, re)''')
    parts.append('''      end select
   end subroutine vrr_dispatch

end module trc_vrr_unrolled
''')

    txt = "\n".join(parts)
    with open(args.output, "w") as fh:
        fh.write(txt)
    print(f"wrote {args.output}: {len(names)} classes, {txt.count(chr(10))} lines")


if __name__ == "__main__":
    main()
