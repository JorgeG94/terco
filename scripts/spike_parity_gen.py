#!/usr/bin/env python3
"""PARITY-SKIP SPIKE (docs/SYMMETRY_PARITY_SKIP.md) -- NOT USED BY THE BUILD.

An experimental copy of the HGP VRR/HRR emission of gen_vrr.py / gen_hrr.py /
gen_perclass.py in which the coordinates in a mask Z are "dead": every centre
of the quartet has that coordinate equal to 0 (the molecular plane of a planar
molecule is Z = {z}; centres on an axis give two coordinates).  Then

    PA_d = PB_d = QC_d = WP_d = WQ_d = AB_d = CD_d = 0      for d in Z,

and every Obara-Saika node (e0|f0) with an odd total power in a dead
coordinate is exactly zero, as is every final component (ab|cd) of odd
combined parity.  This script

  cost    emits the masked VRR + HRR for a class and counts the arithmetic
          (statements, flops) so a cost model per class and mask exists:
              scripts/spike_parity_gen.py cost --lmax 2 -o costs.json
  kernel  emits a standalone Fortran routine for one (class, mask) whose
          VRR/HRR carry only the live nodes, for the prototype test
          (test/spike_parity_kernel.F90) and for host compile-time
          measurements:
              scripts/spike_parity_gen.py kernel 2222 4 -o k2222z4.F90
  report  combines the counts of test/spike_parity.F90 (a CSV) with costs.json.

Mask bits: 1 = x dead, 2 = y dead, 4 = z dead.  Mask 0 reproduces
gen_perclass.py's text exactly (asserted), which is what makes the
comparison fair.  Nothing here writes into src/ or src/generated/.
"""
import argparse
import csv
import importlib.util
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


sys.path.insert(0, HERE)
gv = _load("gv", os.path.join(HERE, "gen_vrr.py"))
gh = _load("gh", os.path.join(HERE, "gen_hrr.py"))

AXN, WPN, QCN, WQN = gv.AXN, gv.WPN, gv.QCN, gv.WQN


def ncart(l):
    return (l + 1)*(l + 2)//2


# ---------------------------------------------------------------- masked VRR
def pathways_z(node, P, I, Z):
    """gv._pathways with the vanishing inputs removed.

    In a dead direction d the first-order terms PA_d*(a-1_d) + WP_d*(...)
    (and the c-side QC_d, WQ_d) are multiplied by exactly zero; the
    second-order and the cross term survive.  A pathway left with no input
    forms an odd (zero) node and is not a pathway.
    """
    ai, ci, m = node
    A, C = P[ai], P[ci]
    out = []
    for d in range(3):
        if A[d] > 0:
            am = list(A); am[d] -= 1; a1 = I[tuple(am)]
            ins = [] if (Z >> d) & 1 else [(a1, ci, m), (a1, ci, m + 1)]
            if A[d] >= 2:
                a2l = list(A); a2l[d] -= 2
                a2 = I[tuple(a2l)]
                ins += [(a2, ci, m), (a2, ci, m + 1)]
            if C[d] >= 1:
                cl = list(C); cl[d] -= 1
                ins.append((a1, I[tuple(cl)], m + 1))
            if ins:
                out.append(("a", d, ins))
    for d in range(3):
        if C[d] > 0:
            cl = list(C); cl[d] -= 1; c1 = I[tuple(cl)]
            ins = [] if (Z >> d) & 1 else [(ai, c1, m), (ai, c1, m + 1)]
            if C[d] >= 2:
                c2l = list(C); c2l[d] -= 2
                c2 = I[tuple(c2l)]
                ins += [(ai, c2, m), (ai, c2, m + 1)]
            if A[d] >= 1:
                al = list(A); al[d] -= 1
                ins.append((I[tuple(al)], c1, m + 1))
            if ins:
                out.append(("c", d, ins))
    return out


def emit_vrr_z(lab, lcd, targets, P, I, nca, Z):
    """gv.emit_class_sieved with the dead-coordinate terms removed.
    Returns (text, final buffer, live node count, live node set)."""
    lt = lab + lcd
    orig = gv._pathways
    gv._pathways = lambda node, P_, I_: pathways_z(node, P_, I_, Z)
    try:
        _full, live = gv.sieve_nodes(lab, lcd, targets, P, I)
        buf = lambda m: (lt - m) % 2
        V = lambda n: f"v({n[0] + (n[1] - 1)*nca},{buf(n[2])})"
        o = []
        w = o.append
        for m in range(lt, -1, -1):
            here = [n for n in live if n[2] == m]
            if not here:
                continue
            w(f"      ! --- level m = {m} ---")
            here.sort(key=lambda n: (sum(P[n[0]]) + sum(P[n[1]]), n[1], n[0]))
            for n in here:
                ai, ci, _ = n
                if ai == 1 and ci == 1:
                    w(f"      {V(n)} = pref*f({m})")
                    continue
                best = None
                for kind, d, ins in gv._pathways(n, P, I):
                    if all(t in live for t in ins):
                        if best is None or len(ins) < len(best[2]):
                            best = (kind, d, ins)
                if best is None:
                    raise RuntimeError(f"sieve left {n} unformable")
                kind, d, _ = best
                A, C = P[ai], P[ci]
                terms = []
                dead = (Z >> d) & 1
                if kind == "a":
                    am = list(A); am[d] -= 1; a1 = I[tuple(am)]
                    if not dead:
                        terms.append(f"{AXN[d]}*{V((a1, ci, m))}"
                                     f" + {WPN[d]}*{V((a1, ci, m + 1))}")
                    if A[d] >= 2:
                        a2l = list(A); a2l[d] -= 2
                        a2 = I[tuple(a2l)]
                        terms.append(f"{float(A[d] - 1)}_dp*oo2z*"
                                     f"({V((a2, ci, m))} - rz*{V((a2, ci, m + 1))})")
                    if C[d] >= 1:
                        cl = list(C); cl[d] -= 1
                        terms.append(f"{float(C[d])}_dp*oo2ze*"
                                     f"{V((a1, I[tuple(cl)], m + 1))}")
                else:
                    cl = list(C); cl[d] -= 1; c1 = I[tuple(cl)]
                    if not dead:
                        terms.append(f"{QCN[d]}*{V((ai, c1, m))}"
                                     f" + {WQN[d]}*{V((ai, c1, m + 1))}")
                    if C[d] >= 2:
                        c2l = list(C); c2l[d] -= 2
                        c2 = I[tuple(c2l)]
                        terms.append(f"{float(C[d] - 1)}_dp*oo2e*"
                                     f"({V((ai, c2, m))} - re*{V((ai, c2, m + 1))})")
                    if A[d] >= 1:
                        al = list(A); al[d] -= 1
                        terms.append(f"{float(A[d])}_dp*oo2ze*"
                                     f"{V((I[tuple(al)], c1, m + 1))}")
                w(f"      {V(n)} = " + " &\n         + ".join(terms))
        return "\n".join(o), buf(0), len(live), live
    finally:
        gv._pathways = orig


# ---------------------------------------------------------------- masked HRR
def comp_parity_zero(Z, *powers):
    """True when the combined power in some dead coordinate is odd."""
    for d in range(3):
        if (Z >> d) & 1 and sum(p[d] for p in powers) % 2:
            return True
    return False


def emit_hrr_z(la, lb, lc, ld, idx, Z, order="lexicographic"):
    """gh.emit_class with the vanishing terms dropped and the odd-parity
    components set to zero. Returns (text, n_zero_components)."""
    lab = la + lb
    nca = gh.ncum(lab)
    o = []
    nzero = 0
    na_, nb_, nc_ = len(gh.cart(la)), len(gh.cart(lb)), len(gh.cart(lc))
    for iid, (dx, dy, dz) in enumerate(gh.cart(ld, order)):
        for iic, (cx, cy, cz) in enumerate(gh.cart(lc, order)):
            for iib, (bx, by, bz) in enumerate(gh.cart(lb, order)):
                for iia, (ax, ay, az) in enumerate(gh.cart(la, order)):
                    slot = 1 + iia + na_*(iib + nb_*(iic + nc_*iid))
                    if comp_parity_zero(Z, (ax, ay, az), (bx, by, bz), (cx, cy, cz), (dx, dy, dz)):
                        o.append(f"         vbuf({slot}) = 0.0_dp")
                        nzero += 1
                        continue
                    B, D = (bx, by, bz), (dx, dy, dz)
                    terms = []
                    for px in range(bx + 1):
                        for py in range(by + 1):
                            for pz in range(bz + 1):
                                if any((Z >> d) & 1 and (px, py, pz)[d] != B[d] for d in range(3)):
                                    continue
                                cb = gh.binom(bx, px)*gh.binom(by, py)*gh.binom(bz, pz)
                                ja = idx[(ax + px, ay + py, az + pz)]
                                fb = [gh.powexpr("abx", bx - px), gh.powexpr("aby", by - py),
                                      gh.powexpr("abz", bz - pz)]
                                for qx in range(dx + 1):
                                    for qy in range(dy + 1):
                                        for qz in range(dz + 1):
                                            if any((Z >> d) & 1 and (qx, qy, qz)[d] != D[d] for d in range(3)):
                                                continue
                                            cd = gh.binom(dx, qx)*gh.binom(dy, qy)*gh.binom(dz, qz)
                                            jc = idx[(cx + qx, cy + qy, cz + qz)]
                                            fd = [gh.powexpr("cdx", dx - qx), gh.powexpr("cdy", dy - qy),
                                                  gh.powexpr("cdz", dz - qz)]
                                            fac = [f for f in fb + fd if f]
                                            c = cb*cd
                                            head = "" if c == 1 else f"{float(c)}_dp*"
                                            body = "*".join(fac)
                                            if body:
                                                body += "*"
                                            terms.append(f"{head}{body}g({ja + (jc - 1)*nca})")
                    assert terms, (la, lb, lc, ld, Z)
                    expr = terms[0]
                    for t in terms[1:]:
                        expr += f" &\n            + {t}"
                    o.append(f"         vbuf({slot}) = {expr}")
    return "\n".join(o), nzero


# --------------------------------------------------------------- the class
def build_class(la, lb, lc, ld, Z, idx_h=None, Pl=None, Il=None):
    """VRR + HRR text for (class, Z); Z = 0 is exactly gen_perclass's."""
    lab, lcd = la + lb, lc + ld
    L = max(la, lb, lc, ld)
    if idx_h is None:
        idx_h = gh.cart_cum(4*L)
    if Pl is None:
        Pl = [None] + gv.cart_list(2*L)
        Il = {q: i for i, q in enumerate(Pl) if q is not None}
    nca = gv.ncum(lab)
    if Z == 0:
        hrr = gh.emit_class(la, lb, lc, ld, idx_h)
        nzero = 0
    else:
        hrr, nzero = emit_hrr_z(la, lb, lc, ld, idx_h, Z)
    reads = sorted({int(x) for x in re.findall(r"g\((\d+)\)", hrr)})
    targets = {((r - 1) % nca + 1, (r - 1)//nca + 1, 0) for r in reads}
    if Z == 0:
        vrr, fin, nlive = gv.emit_class_sieved(lab, lcd, targets, Pl, Il, nca)
    else:
        vrr, fin, nlive, _ = emit_vrr_z(lab, lcd, targets, Pl, Il, nca, Z)
    return dict(vrr=vrr, hrr=hrr, reads=reads, nlive=nlive, nzero=nzero,
                nv=nca*gv.ncum(lcd), nca=nca)


def flops(text):
    """Arithmetic operators of a Fortran fragment (multiplications and
    additions/subtractions; a rough, FMA-agnostic count)."""
    n = len(re.findall(r"\*", text))
    n += len(re.findall(r"(?:^|\s)[+-](?=\s)", text, flags=re.M))
    return n


def nstmt(text):
    return len(re.findall(r"^\s*(?:v|vbuf)\(.*?\)\s*=", text, flags=re.M))



def split_stmts(text, head):
    """Statements (possibly continued with '&') of a fragment, by their
    left-hand-side prefix `head` ('v(' or 'vbuf(')."""
    out, cur = [], []
    for line in text.split("\n"):
        if re.match(r"\s*" + re.escape(head), line) and "=" in line and not cur:
            cur = [line]
        elif cur:
            cur.append(line)
        else:
            continue
        if not line.rstrip().endswith("&"):
            out.append("\n".join(cur)); cur = []
    return out


def runtime_mask_cost(la, lb, lc, ld, Z):
    """Scheme (b): the SHIPPED kernel with launch-uniform guards around the
    statements whose result is parity-odd. A guarded-off statement is
    skipped; the statements that stay keep ALL their terms (a term fed by a
    zero node still runs, since only the guard knows the mask). Returns
    (vrr_flops, acc, hrr_flops) for the shipped text under the mask."""
    lab, lcd = la + lb, lc + ld
    L = max(la, lb, lc, ld)
    Pl = [None] + gv.cart_list(2*L)
    nca = gv.ncum(lab)
    c = build_class(la, lb, lc, ld, 0)
    vf = 0
    live_idx = set()
    for st in split_stmts(c["vrr"], "v("):
        m = re.match(r"\s*v\((\d+),", st)
        idx = int(m.group(1))
        ai, ci = (idx - 1) % nca + 1, (idx - 1)//nca + 1
        if comp_parity_zero(Z, Pl[ai], Pl[ci]):
            continue
        live_idx.add(idx)
        vf += flops(st)
    nv_even = len([1 for i in range(1, c["nv"] + 1)
                   if not comp_parity_zero(Z, Pl[(i - 1) % nca + 1], Pl[(i - 1)//nca + 1])])
    hf = 0
    cs = [gh.cart(l) for l in (la, lb, lc, ld)]
    for st in split_stmts(c["hrr"], "vbuf("):
        slot = int(re.match(r"\s*vbuf\((\d+)\)", st).group(1)) - 1
        na_, nb_, nc_ = len(cs[0]), len(cs[1]), len(cs[2])
        ia = slot % na_; ib = (slot//na_) % nb_; ic = (slot//(na_*nb_)) % nc_; id_ = slot//(na_*nb_*nc_)
        if comp_parity_zero(Z, cs[0][ia], cs[1][ib], cs[2][ic], cs[3][id_]):
            continue
        hf += flops(st)
    return vf, nv_even, hf

# ------------------------------------------------------------- cost model
# Rough per-primitive-quartet overhead that no parity mask can remove: the
# prefactor, W/P/Q geometry and Boys function. Stated, not measured; the
# report prints results with and without it.
def fixed_prim(lt):
    return 90 + 4*lt


DIG_PER_COMP = 13      # sc = wq*v, then six scaled FMAs (J, K blocks)


def dig_blocks(la, lb, lc, ld):
    na, nb, nc, nd = map(ncart, (la, lb, lc, ld))
    return 2*(na*nb + nc*nd + na*nc + na*nd + nb*nc + nb*nd)


def rys_engine(la, lb, lc, ld):
    """Classes the pair bins launch with a Rys kernel (gen_rys.py)."""
    tag = f"{la}{lb}{lc}{ld}"
    return max(la, lb, lc, ld) >= 3 or tag in (
        "2111", "2120", "2121", "2210", "2211", "2221", "2222")


def rys_prim_cost(la, lb, lc, ld, Z):
    """Per primitive quartet, Rys: roots + 2D tables + 1D HRR + the
    component sum. Dead axes: only the even 1D entries are built, and an
    odd-parity component is not summed at all."""
    lt = la + lb + lc + ld
    nr = lt//2 + 1
    lab, lcd = la + lb, lc + ld
    nz = 0
    comps = 0
    from itertools import product
    cs = [gh.cart(l) for l in (la, lb, lc, ld)]
    for a, b, c, d in product(*cs):
        comps += 1
        if Z and comp_parity_zero(Z, a, b, c, d):
            nz += 1
    live_axes = 3 - bin(Z).count("1")
    axis_scale = lambda dead: 0.5 if dead else 1.0
    t2d = nr*sum(((lab + 1)*(lcd + 1)*4 + (la + 1)*(lb + 1)*(lc + 1)*(ld + 1)*2)
                 * axis_scale((Z >> d) & 1) for d in range(3))
    roots = 25*nr
    comp = (comps - nz)*(3*nr)   # two products + the weighted sum per root
    return roots + t2d + comp, comps, nz


def cost_table(lmax):
    rows = {}
    for la in range(lmax + 1):
        for lb in range(la + 1):
            for lc in range(la + 1):
                for ld in range(lc + 1):
                    if (la, lb) < (lc, ld):
                        continue
                    tag = f"{la}{lb}{lc}{ld}"
                    ncomp = ncart(la)*ncart(lb)*ncart(lc)*ncart(ld)
                    lt = la + lb + lc + ld
                    ent = {"ncomp": ncomp, "lt": lt, "rys": rys_engine(la, lb, lc, ld)}
                    for Z in range(8):
                        e = {}
                        if ent["rys"]:
                            pc, comps, nz = rys_prim_cost(la, lb, lc, ld, Z)
                            e.update(prim=pc + fixed_prim(lt), prim_fixed=fixed_prim(lt),
                                     nzero=nz, hrr=0, dig=DIG_PER_COMP*ncomp,
                                     digb=dig_blocks(la, lb, lc, ld), vrr=0, acc=0,
                                     prim_var=pc, stmts=0)
                        else:
                            c = build_class(la, lb, lc, ld, Z)
                            bv, bacc, bh = runtime_mask_cost(la, lb, lc, ld, Z)
                            e.update(b_vrr=bv, b_acc=bacc, b_hrr=bh,
                                     vrr=flops(c["vrr"]), vrr_stmts=nstmt(c["vrr"]),
                                     acc=(len(c["reads"]) if Z else c["nv"]),
                                     hrr=flops(c["hrr"]), hrr_stmts=nstmt(c["hrr"]),
                                     nzero=c["nzero"], prim_fixed=fixed_prim(lt),
                                     dig=DIG_PER_COMP*ncomp, digb=dig_blocks(la, lb, lc, ld),
                                     nlive=c["nlive"], nv=c["nv"], nreads=len(c["reads"]))
                        ent[str(Z)] = e
                    rows[tag] = ent
                    print(tag, "done", file=sys.stderr, flush=True)
    return rows


# ---------------------------------------------------------------- the kernel
def emit_kernel(la, lb, lc, ld, Z, acc="auto"):
    """Standalone Fortran for one class and mask: np primitive quartets in
    (packed in prm(ldp, np): prefactor, F_0..F_lt, PA, QC, WP, WQ, 1/2zeta,
    1/2eta, 1/2(zeta+eta), rho/zeta, rho/eta), the (ab|cd) block out,
    contracted over the primitives like the HGP kernel's g1, then HRR.
    acc: 'all' (as shipped) or 'reads' (only the nodes the HRR reads)."""
    tag = f"{la}{lb}{lc}{ld}"
    c = build_class(la, lb, lc, ld, Z)
    lt = la + lb + lc + ld
    nv, nca = c["nv"], c["nca"]
    nout = ncart(la)*ncart(lb)*ncart(lc)*ncart(ld)
    if acc == "auto":
        acc = "reads" if Z else "all"
    if acc == "all":
        accum = f"         do x = 1, {nv}\n            g1(x) = g1(x) + v(x, {lt % 2})\n         end do\n"
    else:
        accum = "".join(f"         g1({r}) = g1({r}) + v({r}, {lt % 2})\n" for r in c["reads"])
    hrr = re.sub(r"\bg\(", "g1(", c["hrr"])
    name = f"pk{tag}z{Z}{'r' if acc == 'reads' and not Z else ''}"
    k = lt + 3          # 1-based index of PA_x in the packed column
    txt = f"""   subroutine {name}(np, ldp, prm, geo, vbuf)
      !! ({la}{lb}|{lc}{ld}), dead coordinates mask {Z} (1 x, 2 y, 4 z); accumulation: {acc}.
      !! VRR statements {nstmt(c['vrr'])}, HRR statements {nstmt(c['hrr'])}, live nodes {c['nlive']} of {nv}.
      integer, intent(in) :: np, ldp
      real(dp), intent(in) :: prm(ldp, np), geo(6)
      real(dp), intent(out) :: vbuf(*)
      real(dp) :: v({nv}, 0:1), g1({nv}), f(0:{lt})
      real(dp) :: pref, pax, pay, paz, qcx, qcy, qcz, wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re, abx, aby, abz, cdx, cdy, cdz
      integer :: ip, x
      abx = geo(1); aby = geo(2); abz = geo(3); cdx = geo(4); cdy = geo(5); cdz = geo(6)
      g1 = 0.0_dp
      do ip = 1, np
         pref = prm(1, ip); f = prm(2:{lt + 2}, ip)
         pax = prm({k}, ip); pay = prm({k + 1}, ip); paz = prm({k + 2}, ip)
         qcx = prm({k + 3}, ip); qcy = prm({k + 4}, ip); qcz = prm({k + 5}, ip)
         wpx = prm({k + 6}, ip); wpy = prm({k + 7}, ip); wpz = prm({k + 8}, ip)
         wqx = prm({k + 9}, ip); wqy = prm({k + 10}, ip); wqz = prm({k + 11}, ip)
         oo2z = prm({k + 12}, ip); oo2e = prm({k + 13}, ip); oo2ze = prm({k + 14}, ip)
         rz = prm({k + 15}, ip); re = prm({k + 16}, ip)
{c['vrr']}
{accum}      end do
{hrr}
   end subroutine {name}
"""
    return name, dict(tag=tag, Z=Z, acc=acc, lt=lt, nout=nout, text=txt, name=name,
                      vrr=nstmt(c["vrr"]), hrr=nstmt(c["hrr"]))


def kernel_module(specs, modname, acc="auto"):
    """specs: (tag, mask) or (tag, mask, acc) triples."""
    ents = []
    for sp in specs:
        tag, Z = sp[0], sp[1]
        n, e = emit_kernel(*(int(ch) for ch in tag), Z, sp[2] if len(sp) > 2 else acc)
        ents.append(e)
    names = ", ".join(e["name"] for e in ents)
    tbl = []
    for i, e in enumerate(ents, start=1):
        tbl.append(f"      cases({i})%name = '{e['tag']}'; cases({i})%zmask = {e['Z']}; "
                   f"cases({i})%kind = '{e['acc']}'; cases({i})%lt = {e['lt']}; "
                   f"cases({i})%nout = {e['nout']}; cases({i})%run => {e['name']}")
    return (f"! GENERATED by scripts/spike_parity_gen.py -- spike, not part of the build.\n"
            f"module {modname}\n   use pic_types, only: dp\n   implicit none\n   private\n"
            f"   public :: case_t, init_cases, NCASE\n"
            f"   integer, parameter :: NCASE = {len(ents)}\n"
            f"   abstract interface\n      subroutine run_i(np, ldp, prm, geo, vbuf)\n         import :: dp\n"
            f"         integer, intent(in) :: np, ldp\n         real(dp), intent(in) :: prm(ldp, np), geo(6)\n"
            f"         real(dp), intent(out) :: vbuf(*)\n      end subroutine run_i\n   end interface\n"
            f"   type :: case_t\n      character(len=4) :: name = ''\n      character(len=8) :: kind = ''\n"
            f"      integer :: zmask = 0, lt = 0, nout = 0\n      procedure(run_i), pointer, nopass :: run => null()\n"
            f"   end type case_t\ncontains\n   subroutine init_cases(cases)\n      type(case_t), intent(out) :: cases(NCASE)\n"
            + "\n".join(tbl) + "\n   end subroutine init_cases\n\n" + "\n".join(e["text"] for e in ents)
            + f"end module {modname}\n")


# ------------------------------------------------------------------- report
def parity_counts(l):
    """number of Cartesian components of shell l per parity vector (bit d =
    power of d odd)."""
    cnt = [0]*8
    for (x, y, z) in gh.cart(l):
        cnt[(x & 1) | ((y & 1) << 1) | ((z & 1) << 2)] += 1
    return cnt


def xor_conv(a, b):
    out = [0]*8
    for i in range(8):
        for j in range(8):
            out[i ^ j] += a[i]*b[j]
    return out


def zero_comps_flipset(ls, S):
    """Components of class ls whose parity vector is rejected by the stabilizer
    S (a set of flip vectors 0..7): nonzero iff popcount(g & P) even for all g in S."""
    cnt = [1] + [0]*7
    for l in ls:
        cnt = xor_conv(cnt, parity_counts(l))
    bad = 0
    for P in range(8):
        if any(bin(g & P).count("1") % 2 for g in S):
            bad += cnt[P]
    return bad


def zero_comps_Z(ls, Z):
    return zero_comps_flipset(ls, [1 << d for d in range(3) if (Z >> d) & 1])


FIXED_SCALE = 1.0   # multiplies the per-primitive-quartet overhead (Boys, geometry)


def tier_costs(ent, Z, nprim, nq):
    """Total flops of nq quartets / nprim primitive quartets under each tier:
    T0 shipped; T1 digestion skip only; T2 + HRR of the masked variant;
    T3 full masked variants (scheme a); TB runtime guards (scheme b), last."""
    e0, ez = ent["0"], ent[str(Z)]
    ncomp = ent["ncomp"]
    nzero = ez["nzero"]
    dig_live = DIG_PER_COMP*(ncomp - nzero)
    if ent["rys"]:
        P0 = e0["prim_var"] + e0["prim_fixed"]*FIXED_SCALE
        Pz = ez["prim_var"] + ez["prim_fixed"]*FIXED_SCALE
        T0 = nprim*P0 + nq*(e0["dig"] + e0["digb"])
        T1 = nprim*P0 + nq*(dig_live + e0["digb"])
        T2 = T1
        T3 = nprim*Pz + nq*(dig_live + e0["digb"])
        # runtime guard on the component loop: components skipped, 2D tables not trimmed
        comp_sav = (ez["nzero"])*(3*(ent["lt"]//2 + 1))
        TB = nprim*(P0 - comp_sav) + nq*(dig_live + e0["digb"])
        return T0, T1, T2, T3, nprim*P0, nprim*Pz, 0, 0, TB
    P0 = e0["prim_fixed"]*FIXED_SCALE + e0["vrr"] + e0["acc"]
    Pz = ez["prim_fixed"]*FIXED_SCALE + ez["vrr"] + ez["acc"]
    C0 = e0["hrr"] + e0["dig"] + e0["digb"]
    T0 = nprim*P0 + nq*C0
    T1 = nprim*P0 + nq*(e0["hrr"] + dig_live + e0["digb"])
    T2 = nprim*P0 + nq*(ez["hrr"] + dig_live + e0["digb"])
    T3 = nprim*Pz + nq*(ez["hrr"] + dig_live + e0["digb"])
    PB = e0["prim_fixed"]*FIXED_SCALE + ez["b_vrr"] + ez["b_acc"]
    TB = nprim*PB + nq*(ez["b_hrr"] + dig_live + e0["digb"])
    return T0, T1, T2, T3, nprim*P0, nprim*Pz, nq*e0["hrr"], nq*ez["hrr"], TB


def report(csvfile, costfile, plan, top=12, boys=True):
    """Per-class and total zero fractions and cost-weighted savings."""
    costs = json.load(open(costfile))
    byclass = {}
    with open(csvfile) as fh:
        for r in csv.DictReader(fh):
            if r["plan"] != plan:
                continue
            ls = [int(r[k]) for k in ("la", "lb", "lc", "ld")]
            tag = "".join(map(str, ls))
            Z = int(r["Z"])
            S = [g for g in range(8) if (int(r["S"]) >> g) & 1]
            nq, nprim = int(r["nq"]), int(r["nprim"])
            ent = costs[tag]
            T = tier_costs(ent, Z, nprim, nq)
            b = byclass.setdefault(tag, dict(nq=0, nprim=0, ncomp=0, zS=0, zZ=0, T=[0.0]*5, T0v=0.0,
                                             rys=ent["rys"], prim_only=0.0, nrows=0))
            b["nq"] += nq; b["nprim"] += nprim; b["ncomp"] += nq*ent["ncomp"]
            b["zS"] += nq*zero_comps_flipset(ls, S)
            b["zZ"] += nq*zero_comps_Z(ls, Z)
            for i in range(4):
                b["T"][i] += T[i]
            b["T"][4] += T[8]
    tot = dict(nq=0, ncomp=0, zS=0, zZ=0, T=[0.0]*5)
    for b in byclass.values():
        for k in ("nq", "ncomp", "zS", "zZ"):
            tot[k] += b[k]
        for i in range(5):
            tot["T"][i] += b["T"][i]
    T0 = tot["T"][0]
    print(f"plan {plan}: quartets {tot['nq']}, component integrals {tot['ncomp']:.4g}")
    print(f"  parity-zero components: stabilizer rule {tot['zS']/tot['ncomp']:.4f}, "
          f"dead-coordinate rule {tot['zZ']/tot['ncomp']:.4f}")
    print(f"  cost-weighted work removed (model flops {T0:.4g}):  "
          f"digestion only {1 - tot['T'][1]/T0:.4f}   + HRR/contracted level {1 - tot['T'][2]/T0:.4f}   "
          f"+ primitive level (full variants) {1 - tot['T'][3]/T0:.4f}")
    print(f"  runtime-mask kernel (scheme b) removes {1 - tot['T'][4]/T0:.4f}  ({T0/tot['T'][4]:.3f}x)")
    print(f"  speed-up if realised: digestion-only {T0/tot['T'][1]:.3f}x, +HRR {T0/tot['T'][2]:.3f}x, "
          f"+primitive {T0/tot['T'][3]:.3f}x")
    print("  class   engine   quartets   share_cost  zeroS   zeroZ   rm_dig  rm_HRR  rm_prim  rm_b")
    rows = sorted(byclass.items(), key=lambda kv: -kv[1]["T"][0])
    for tag, b in rows[:top]:
        t0 = b["T"][0]
        print(f"  ({tag[0]}{tag[1]}|{tag[2]}{tag[3]})  {'Rys' if b['rys'] else 'HGP'}   {b['nq']:>10d}   "
              f"{t0/T0:8.4f}  {b['zS']/b['ncomp']:6.3f}  {b['zZ']/b['ncomp']:6.3f}  "
              f"{1 - b['T'][1]/t0:6.3f}  {1 - b['T'][2]/t0:6.3f}  {1 - b['T'][3]/t0:6.3f}  {1 - b['T'][4]/t0:6.3f}")
    return tot, byclass


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("cost")
    p.add_argument("--lmax", type=int, default=2)
    p.add_argument("-o", default="costs.json")
    p = sub.add_parser("kernel")
    p.add_argument("tag")
    p.add_argument("mask", type=int)
    p.add_argument("-o", default="-")
    p = sub.add_parser("module")
    p.add_argument("specs", help="comma list tag:mask[:r], e.g. 1111:0,1111:4,1111:0:r (r = accumulate only the "
                   "nodes the HRR reads, the trimmed baseline)")
    p.add_argument("--name", default="spike_parity_kernels")
    p.add_argument("--acc", default="auto")
    p.add_argument("-o", default="-")
    p = sub.add_parser("report")
    p.add_argument("csv")
    p.add_argument("costs")
    p.add_argument("--plan", default="full")
    p.add_argument("--top", type=int, default=12)
    p.add_argument("--fixed-scale", type=float, default=1.0)
    p = sub.add_parser("check0")
    p.add_argument("tags", help="comma list of tags: mask 0 must equal gen_perclass")
    args = ap.parse_args()
    if args.cmd == "cost":
        json.dump(cost_table(args.lmax), open(args.o, "w"), indent=0)
    elif args.cmd == "kernel":
        name, e = emit_kernel(*(int(c) for c in args.tag), args.mask)
        (sys.stdout if args.o == "-" else open(args.o, "w")).write(e["text"])
    elif args.cmd == "module":
        specs = []
        for s_ in args.specs.split(","):
            f_ = s_.split(":")
            specs.append((f_[0], int(f_[1])) + (("reads",) if len(f_) > 2 else ()))
        txt = kernel_module(specs, args.name, args.acc)
        (sys.stdout if args.o == "-" else open(args.o, "w")).write(txt)
    elif args.cmd == "report":
        global FIXED_SCALE
        FIXED_SCALE = args.fixed_scale
        report(args.csv, args.costs, args.plan, args.top)
    elif args.cmd == "check0":
        # mask 0 through the masked emitters equals the shipped emitters
        for tag in args.tags.split(","):
            la, lb, lc, ld = (int(c) for c in tag)
            lab, lcd = la + lb, lc + ld
            L = max(la, lb, lc, ld)
            idx_h = gh.cart_cum(4*L)
            Pl = [None] + gv.cart_list(2*L)
            Il = {q: i for i, q in enumerate(Pl) if q is not None}
            hrr = gh.emit_class(la, lb, lc, ld, idx_h)
            hz, _ = emit_hrr_z(la, lb, lc, ld, idx_h, 0)
            nca = gv.ncum(lab)
            reads = {int(x) for x in re.findall(r"g\((\d+)\)", hrr)}
            targets = {((r - 1) % nca + 1, (r - 1)//nca + 1, 0) for r in reads}
            v0, f0, n0 = gv.emit_class_sieved(lab, lcd, targets, Pl, Il, nca)[:3]
            vz, fz, nz, _ = emit_vrr_z(lab, lcd, targets, Pl, Il, nca, 0)
            print(tag, "HRR identical:", hrr == hz, " VRR identical:", v0 == vz and f0 == fz and n0 == nz)


if __name__ == "__main__":
    main()
