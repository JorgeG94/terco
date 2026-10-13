#!/usr/bin/env python3
"""Numerical check: the HRR DAG (gen_hrr.py's `emit_class_dag`, the only
form the generators emit) against its own closed-form expansion
(`emit_class`, kept in gen_hrr.py purely as this check's oracle -- it is
not reachable from any generator's command line), on random values for the
VRR's contracted output buffer and random AB/CD vectors.

WHY A SEPARATE, PURE-PYTHON CHECK
----------------------------------
`scripts/rotaxis/check_numeric.py` validates a different algebra (the
rotated-axis MMD path) against pyscf, and needs pyscf and a geometry to do
it. The HRR DAG's risk is narrower and purely algebraic: does building the
HRR from the targets back to the VRR buffer, sieved and shared the way
gen_vrr.py already builds the VRR, change what gets computed relative to
summing each of the na*nb*nc*nd output components independently in closed
form? That question has nothing to do with a basis set, a molecule, or the
Boys function, and answering it does not need Fortran, GPU or pyscf -- only
the generator's own recursion, evaluated numerically instead of emitted as
text, against the closed form it is meant to agree with to machine
precision for every (la,lb,lc,ld) up to TRC_LMAX=2.

Every output component's formula is a FINITE sum (emit_class) or a DAG
walk (emit_class_dag) over the same contracted VRR buffer; neither touches
a primitive loop, a screen, or a density. Agreeing on random buffer values
is therefore not a weak test -- there is no cancellation this could hide
behind, because the comparison is against the exact closed form, not
against a stored reference number.

Run as:  python3 scripts/check_hrr_dag_numeric.py
Exits nonzero if any class's worst relative error exceeds the tolerance.
"""
import itertools
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_hrr as gh  # noqa: E402
import gen_perclass as gp  # noqa: E402
import gen_vrr as gv  # noqa: E402

LMAX = 2
# ~1e-13 is the per-operation rounding floor; the DAG and closed forms sum
# the SAME terms in a DIFFERENT association order (the DAG's whole point),
# so double rounding compounds across the recursion depth rather than
# cancelling. Measured worst case across six seeds, AB/CD in [-2,2] and the
# VRR buffer in [-1,1]: 1.5e-11 for (dd|dd)-depth classes. 1e-10 keeps a
# comfortable margin above that without hiding a real algebra error, which
# would show as 1e-2 to 1 (a missing or doubled term), not 1e-11.
TOL = 1e-10


def binom(n, k):
    return gh.binom(n, k)


def eval_closed(la, lb, lc, ld, idx, nca, g, ab, cd):
    """Numeric twin of gen_hrr.emit_class: same term construction, a float
    accumulator instead of a string. Walks the SAME (px,py,pz)/(qx,qy,qz)
    splits in the SAME order, so a mismatch here would mean the two
    functions have drifted, not that this check is independently wrong."""
    abx, aby, abz = ab
    cdx, cdy, cdz = cd
    out = {}
    for d in gh.cart(ld):
        for c in gh.cart(lc):
            for b in gh.cart(lb):
                for a in gh.cart(la):
                    ax, ay, az = a
                    bx, by, bz = b
                    cx, cy, cz = c
                    dx, dy, dz = d
                    total = 0.0
                    for px in range(bx + 1):
                        for py in range(by + 1):
                            for pz in range(bz + 1):
                                cb = binom(bx, px)*binom(by, py)*binom(bz, pz)
                                ja = idx[(ax + px, ay + py, az + pz)]
                                fb = (abx**(bx - px))*(aby**(by - py))*(abz**(bz - pz))
                                for qx in range(dx + 1):
                                    for qy in range(dy + 1):
                                        for qz in range(dz + 1):
                                            ccd = binom(dx, qx)*binom(dy, qy)*binom(dz, qz)
                                            jc = idx[(cx + qx, cy + qy, cz + qz)]
                                            fd = (cdx**(dx - qx))*(cdy**(dy - qy))*(cdz**(dz - qz))
                                            total += cb*ccd*fb*fd*g[ja + (jc - 1)*nca]
                    out[(a, b, c, d)] = total
    return out


def eval_dag(la, lb, lc, ld, idx, g, ab, cd):
    """Numeric evaluation of gen_hrr.build_class_dag: walk the SAME node
    list the Fortran emitter would walk, with floats instead of statements."""
    coefval = {"abx": ab[0], "aby": ab[1], "abz": ab[2],
               "cdx": cd[0], "cdy": cd[1], "cdz": cd[2]}
    nodes, rec = gh.build_class_dag(la, lb, lc, ld, idx)
    h = {}

    def val(ref):
        kind, v = ref
        return g[v] if kind == "g" else h[v]

    for nid, r1, coef, r2 in nodes:
        h[nid] = val(r1) + coefval[coef]*val(r2)

    out = {}
    for d in gh.cart(ld):
        for c in gh.cart(lc):
            for b in gh.cart(lb):
                for a in gh.cart(la):
                    out[(a, b, c, d)] = val(rec(a, b, c, d))
    return out


def eval_role(h_nodes, role, coefval, g):
    """Mirror what emit_role_kernel.py emits for ONE role: walk its own
    `hset` (already topological -- node ids are assigned in dependency
    order, so a role's subset sorted by id stays in order) through its
    OWN `hmap`, exactly as the generated Fortran's role-local `h(...)`
    does, then copy `vb` through the SAME map. A bug in the role split
    (wrong closure, wrong renumbering) would show up here as a role
    computing a different value than the class's own DAG for the slots
    that role owns -- see `check_role_split` below."""
    hmap = role["hmap"]
    hlocal = {}
    for hid in role["hset"]:
        text, _deps = h_nodes[hid]
        rhs = text.split("=", 1)[1]
        m = re.match(r"\s*(g|h)\((\d+)\)\s*\+\s*(\w+)\*(g|h)\((\d+)\)\s*$", rhs)
        k1, v1, coef, k2, v2 = m.groups()
        a = g[int(v1)] if k1 == "g" else hlocal[hmap[int(v1)]]
        b = g[int(v2)] if k2 == "g" else hlocal[hmap[int(v2)]]
        hlocal[hmap[hid]] = a + coefval[coef]*b
    return hlocal, hmap


def check_role_split(rng):
    """Every ROLE_SPLITS class: does each role's own compacted HRR closure
    (role_plan's hset/hmap/vb/bmap -- the data emit_role_kernel renders to
    Fortran) reproduce the SAME value the class's own (unsplit) DAG gives
    for exactly the vbuf slots that role owns? This is the numeric check
    above's counterpart for the warp-role split (docs/GENERATORS.md):
    emit_class_dag's correctness says nothing about whether role_plan's
    per-role closure and renumbering is also right, and that logic -- not
    the DAG recursion itself -- is where a role split actually risks
    going wrong (an incomplete closure, a stale id after renumbering)."""
    worst = 0.0
    worst_tag = None
    lst, idx_v, info = gv.build_tables(2*LMAX)
    idx_h = gh.cart_cum(4*LMAX)
    for tag, shell in gp.ROLE_SPLITS.items():
        la, lb, lc, ld = (int(c) for c in tag)
        nca_ = gv.ncum(la + lb)
        hrr_pre, _ = gh.emit_class_dag(la, lb, lc, ld, idx_h)
        reads = {int(x) for x in re.findall(r"g\((\d+)\)", hrr_pre)}
        targets = {((r - 1) % nca_ + 1, (r - 1)//nca_ + 1, 0) for r in reads}
        Pl = [None] + gv.cart_list(2*LMAX)
        Il = {q: i for i, q in enumerate(Pl) if q is not None}
        vrr, _fin, _nlive = gv.emit_class_sieved(la + lb, lc + ld, targets, Pl, Il, nca_)
        vrr = "\n".join("      " + ln if ln.strip() else ln for ln in vrr.split("\n"))
        vrr += f"\n               cur = {(la + lb + lc + ld) % 2}"
        hrr, _nnodes = gh.emit_class_dag(la, lb, lc, ld, idx_h)

        _stmts, roles, (h_nodes, vbuf_refs) = gp.role_plan(la, lb, lc, ld, vrr, hrr, shell)
        ncc = gv.ncum(lc + ld)
        nv = nca_*ncc
        for _trial in range(3):
            g = {i: rng.uniform(-1.0, 1.0) for i in range(1, nv + 1)}
            coefval = {"abx": rng.uniform(-2, 2), "aby": rng.uniform(-2, 2), "abz": rng.uniform(-2, 2),
                      "cdx": rng.uniform(-2, 2), "cdy": rng.uniform(-2, 2), "cdz": rng.uniform(-2, 2)}
            memo = {}

            def val(ref):
                kind, v = ref
                if kind == "g":
                    return g[v]
                if v in memo:
                    return memo[v]
                text, _d = h_nodes[v]
                rhs = text.split("=", 1)[1]
                m = re.match(r"\s*(g|h)\((\d+)\)\s*\+\s*(\w+)\*(g|h)\((\d+)\)\s*$", rhs)
                k1, v1, coef, k2, v2 = m.groups()
                a = val(("g", int(v1)) if k1 == "g" else ("h", int(v1)))
                b = val(("g", int(v2)) if k2 == "g" else ("h", int(v2)))
                r = a + coefval[coef]*b
                memo[v] = r
                return r

            for role in roles:
                hlocal, hmap = eval_role(h_nodes, role, coefval, g)
                for slot in role["vb"]:
                    kind, v = vbuf_refs[slot]
                    ref_val = val((kind, v))
                    got = g[v] if kind == "g" else hlocal[hmap[v]]
                    denom = max(abs(ref_val), abs(got), 1e-300)
                    rel = abs(ref_val - got)/denom
                    if rel > worst:
                        worst = rel
                        worst_tag = (tag, shell, slot)
    return worst, worst_tag


def main():
    rng = random.Random(20261008)
    idx = gh.cart_cum(4*LMAX)
    worst = 0.0
    worst_tag = None
    nclass = 0
    for la, lb, lc, ld in itertools.product(range(LMAX + 1), repeat=4):
        lab, lcd = la + lb, lc + ld
        nca, ncd = gh.ncum(lab), gh.ncum(lcd)
        g = {}
        for ja in range(1, nca + 1):
            for jc in range(1, ncd + 1):
                g[ja + (jc - 1)*nca] = rng.uniform(-1.0, 1.0)
        ab = (rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0))
        cd = (rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0))

        closed = eval_closed(la, lb, lc, ld, idx, nca, g, ab, cd)
        dag = eval_dag(la, lb, lc, ld, idx, g, ab, cd)
        nclass += 1

        assert closed.keys() == dag.keys()
        for key, cval in closed.items():
            dval = dag[key]
            denom = max(abs(cval), abs(dval), 1e-300)
            rel = abs(cval - dval)/denom
            if rel > worst:
                worst = rel
                worst_tag = (la, lb, lc, ld, key)

    print(f"checked {nclass} classes up to LMAX={LMAX}, "
          f"worst relative error {worst:.3e} at {worst_tag}")
    rc = 0
    if worst > TOL:
        print(f"FAIL: worst relative error {worst:.3e} exceeds tolerance {TOL:.1e}")
        rc = 1
    else:
        print(f"PASS: closed-form and DAG HRR agree to {TOL:.1e}")

    role_worst, role_tag = check_role_split(rng)
    print(f"checked the warp-role split for every ROLE_SPLITS class, "
          f"worst relative error {role_worst:.3e} at {role_tag}")
    if role_worst > TOL:
        print(f"FAIL: role split worst relative error {role_worst:.3e} exceeds tolerance {TOL:.1e}")
        rc = 1
    else:
        print(f"PASS: every role's compacted HRR closure agrees with the class's own DAG to {TOL:.1e}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
