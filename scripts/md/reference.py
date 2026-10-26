"""numpy oracle for the MD J engine: executes the generator's own statements.

Not imported by any generator (numpy, scipy and pyscf are offline-only).
`python3 -m md.reference` from scripts/ builds J for a few molecules the way
the generated kernels will -- ket Hermite density, two-centre Hermite
contraction per (P, Q) primitive pair, bra back-transform -- running the
exact `Stmt` lists `hermite` and `rtree` produce, every batch split, and
compares it against pyscf's int2e_cart contracted with the same density.

That checks the algebra, the signs, the sieve and the batching before a line
of Fortran exists; the Fortran is then checked against terco's own Fock.
"""

import sys

import numpy as np
from scipy.special import gammainc, gamma

from .hermite import DIRS, herm, herm_index, e1d_statements, e_operand
from . import rtree

TWO_PI_2_5 = 2.0 * np.pi ** 2.5


def boys(mmax, t):
    """F_0..F_mmax(t) on an array, stable at t = 0."""
    t = np.asarray(t, dtype=float)
    out = np.empty((mmax + 1,) + t.shape)
    small = t < 1e-12
    ts = np.where(small, 1.0, t)
    for m in range(mmax + 1):
        a = m + 0.5
        val = gamma(a) * gammainc(a, ts) / (2.0 * ts ** a)
        out[m] = np.where(small, 1.0 / (2 * m + 1), val)
    return out


def run(stmts, env):
    for s in stmts:
        acc = 0.0
        for tm in s.terms:
            v = tm.const
            if tm.factor is not None:
                v = v * env[tm.factor]
            if tm.operand is not None:
                v = v * env[tm.operand]
            acc = acc + v
        env[s.lhs] = acc
    return env


def cart(l):
    return [(lx, ly, l - lx - ly) for lx in range(l, -1, -1) for ly in range(l - lx, -1, -1)]


def common_fac_sp(l):
    return {0: 0.282094791773878143, 1: 0.488602511902919921}.get(l, 1.0)


def shells_of(mol):
    """(l, exps, coefs, centre, ao offset) per contracted column."""
    out = []
    ao = 0
    for ib in range(mol.nbas):
        l = mol.bas_angular(ib)
        e = mol.bas_exp(ib)
        cs = mol._libcint_ctr_coeff(ib) * common_fac_sp(l)
        r = mol.bas_coord(ib)
        for k in range(cs.shape[1]):
            out.append((l, e, cs[:, k].copy(), r, ao))
            ao += (l + 1) * (l + 2) // 2
    return out, ao


def prim_pairs(shells):
    """Every ordered shell pair's primitive pairs, grouped by l_a + l_b."""
    pairs = {}
    for ia, (la, ea, ca, ra, oa) in enumerate(shells):
        for ib, (lb, eb, cb, rb, ob) in enumerate(shells):
            ab2 = np.sum((ra - rb) ** 2)
            for a, c1 in zip(ea, ca):
                for b, c2 in zip(eb, cb):
                    p = a + b
                    P = (a * ra + b * rb) / p
                    c = c1 * c2 * np.exp(-a * b / p * ab2)
                    pairs.setdefault(la + lb, []).append(
                        dict(ia=ia, ib=ib, la=la, lb=lb, p=p, P=P, c=c,
                             PA=P - ra, PB=P - rb, oa=oa, ob=ob))
    return pairs


def e3d(pp):
    """E^{ab}_{tuv} for every component pair: array (ncart_a, ncart_b, nherm)."""
    la, lb = pp["la"], pp["lb"]
    env = {"oo2p": 0.5 / pp["p"]}
    for d in range(3):
        env[f"p{DIRS[d]}a"] = pp["PA"][d]
        env[f"p{DIRS[d]}b"] = pp["PB"][d]
        run(e1d_statements(la, lb, d), env)

    def e1(d, i, j, t):
        if t > i + j:
            return 0.0
        n = e_operand(d, i, j, t)
        return 1.0 if n is None else env[n]

    hs = herm(la + lb)
    ca, cb = cart(la), cart(lb)
    out = np.zeros((len(ca), len(cb), len(hs)))
    for x, a in enumerate(ca):
        for y, b in enumerate(cb):
            for k, (t, u, v) in enumerate(hs):
                out[x, y, k] = (e1(0, a[0], b[0], t) * e1(1, a[1], b[1], u)
                                * e1(2, a[2], b[2], v))
    return out


def md_j(mol, dm, nbatch=1, sieve=True):
    shells, nao = shells_of(mol)
    pairs = prim_pairs(shells)
    # Stage 1: ket Hermite density, sign (-1)^|tau| and coefficient folded in.
    for lq, lst in pairs.items():
        sg = np.array([(-1.0) ** sum(h) for h in herm(lq)])
        for pp in lst:
            E = e3d(pp)
            pp["E"] = E
            nca, ncb = E.shape[:2]
            blk = dm[pp["oa"]:pp["oa"] + nca, pp["ob"]:pp["ob"] + ncb]
            pp["Dt"] = pp["c"] * sg * np.einsum("xyk,xy->k", E, blk)
            pp["Jt"] = np.zeros(len(herm(lq)))
    # Stage 2: per (lp, lq) class, all (P, Q) at once, statement by statement.
    for lp, bra in pairs.items():
        for lq, ket in pairs.items():
            L = lp + lq
            p = np.array([x["p"] for x in bra])[:, None]
            q = np.array([x["p"] for x in ket])[None, :]
            PQ = (np.array([x["P"] for x in bra])[:, None, :]
                  - np.array([x["P"] for x in ket])[None, :, :])
            alpha = p * q / (p + q)
            pref = TWO_PI_2_5 / (p * q * np.sqrt(p + q))
            F = boys(L, alpha * np.sum(PQ ** 2, axis=2))
            Dt = np.array([x["Dt"] for x in ket])          # (nQ, nherm(lq))
            hp, hq = herm(lp), herm(lq)
            Jt = np.zeros((len(bra), len(hp)))
            for g in rtree.batch_dags(L, nbatch, sieve=sieve):
                env = {f"pq{DIRS[d]}": PQ[:, :, d] for d in range(3)}
                for m in g.leaves():
                    env[rtree.leaf_name(m)] = pref * (-2.0 * alpha) ** m * F[m]
                run(g.statements(), env)
                for (r, _) in g.targets:
                    val = env[rtree.leaf_name(0)] if r == (0, 0, 0) else env[rtree.r_name(r, 0)]
                    for i, t in enumerate(hp):
                        tau = tuple(r[k] - t[k] for k in range(3))
                        if min(tau) < 0 or sum(tau) > lq:
                            continue
                        j = herm_index(lq)[tau]
                        Jt[:, i] += val @ Dt[:, j]
            for i, x in enumerate(bra):
                x["Jt"] = x["Jt"] + Jt[i]
    # Stage 3: back to Cartesian, coefficient folded in.
    J = np.zeros((nao, nao))
    for lp, lst in pairs.items():
        for pp in lst:
            nca, ncb = pp["E"].shape[:2]
            J[pp["oa"]:pp["oa"] + nca, pp["ob"]:pp["ob"] + ncb] += \
                pp["c"] * np.einsum("xyk,k->xy", pp["E"], pp["Jt"])
    return J


def _mols():
    from pyscf import gto
    water = "O 0.0 0.1 -0.05; H 0.76 0.55 0.02; H -0.71 0.6 0.1"
    yield "water 6-31G* (d)", gto.M(atom=water, basis="6-31g*", cart=True)
    yield "water cc-pVDZ (general contraction)", gto.M(atom=water, basis="cc-pvdz", cart=True)
    yield "HF cc-pVTZ (f)", gto.M(atom="F 0.1 -0.05 0.02; H 0.2 0.3 0.88", basis="cc-pvtz", cart=True)


def main():
    rng = np.random.default_rng(7)
    worst = 0.0
    for name, mol in _mols():
        eri = mol.intor("int2e_cart")
        nao = eri.shape[0]
        a = rng.standard_normal((nao, nao))
        dm = a + a.T
        jref = np.einsum("abcd,cd->ab", eri, dm)
        for nb, sv in ((1, False), (1, True), (3, True), (6, True), (9, True)):
            j = md_j(mol, dm, nbatch=nb, sieve=sv)
            err = np.max(np.abs(j - jref)) / np.max(np.abs(jref))
            worst = max(worst, err)
            print(f"{name:40s} nbatch={nb} sieve={sv!s:5s} max rel err {err:.2e}")
    print(f"worst {worst:.2e}")
    return 0 if worst < 1e-12 else 1


if __name__ == "__main__":
    sys.exit(main())
