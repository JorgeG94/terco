"""numpy oracle for the MD Coulomb gradient: executes the generator's statements.

Offline only (numpy, scipy, pyscf). `python3 -m md.reference_grad` from
scripts/ builds d/dX Tr(D J[D]) the way the generated gradient kernels will,
and checks it against pyscf's `int2e_ip1`.

THE DERIVATION (what gen_md_j.py's gradient kernels run)

Tr(D J[D]) = sum D_ab (ab|cd) D_cd is symmetric under bra <-> ket, so with
rho = sum_cd D_cd phi_c phi_d the nuclear derivative needs the bra side only:

    d/dX Tr(D J[D]) = 2 sum_ab D_ab  d_X (ab | rho),   d_X acting on a and b.

A bra primitive pair is Omega_ab = K_AB sum_t E^ab_t Lambda_t(r; P). Its
centre derivatives follow from those of the two Cartesian factors,

    d/dA_x G_a = 2 alpha G_{a+1_x} - a_x G_{a-1_x},

which is a pair distribution of Hermite order lp + 1 (E^{a+1,b}) or lp - 1,
and from translation invariance, (d/dA + d/dB) Omega = -grad_r Omega, with
grad_r Lambda_t = -Lambda_{t+1_x}:

    (d/dA_x + d/dB_x) Omega_ab = sum_t E^ab_t Lambda_{t + 1_x}.

So the whole bra needs ONE new Hermite quantity, the potential of rho on P for
Hermite orders up to lp + 1,

    Jt'[P]_t = sum_Q sum_tau R_{t + tau}(P, Q) Dt[Q]_tau,   |t| <= lp + 1,

which is the existing Hermite stage with the bra order raised by one and the ket
stage unchanged, and then, per component pair (a, b) and direction x,

    gA_x = c_P ( 2 alpha sum_t E^{a+1_x,b}_t Jt'_t - a_x sum_t E^{a-1_x,b}_t Jt'_t )
    gT_x = c_P   sum_t E^{ab}_t Jt'_{t + 1_x}                   (A + B together)

with gradient[atom(A)] += w gA, gradient[atom(B)] += w (gT - gA), w = 2 D'_ab
(D' = D_ab + D_ba off the shell diagonal: one shell pair stands for both
orders, as in the ket stage). A bra pair on ONE atom is not special: its
atom collects gT, the ket's charge is elsewhere.
"""

import sys

import numpy as np

from .hermite import DIRS, herm, herm_index, e1d_statements, e_operand
from . import rtree
from .reference import boys, run, cart, TWO_PI_2_5, shells_of


def atoms_of(mol):
    """Atom index of every contracted column, in shells_of's order."""
    out = []
    for ib in range(mol.nbas):
        nctr = mol.bas_nctr(ib)
        out.extend([mol.bas_atom(ib)] * nctr)
    return out


def _e_env(la, lb, pp):
    env = {"oo2p": 0.5 / pp["p"]}
    for d in range(3):
        env[f"p{DIRS[d]}a"] = pp["PA"][d]
        env[f"p{DIRS[d]}b"] = pp["PB"][d]
        run(e1d_statements(la, lb, d), env)
    return env


def _e1(env, d, i, j, t):
    if i < 0 or t < 0 or t > i + j:
        return 0.0
    n = e_operand(d, i, j, t)
    return 1.0 if n is None else env[n]


def _e3(env, a, b, t):
    return (_e1(env, 0, a[0], b[0], t[0]) * _e1(env, 1, a[1], b[1], t[1])
            * _e1(env, 2, a[2], b[2], t[2]))


def unordered_pairs(shells, atoms):
    """Primitive pairs of every shell pair ia >= ib (one stands for both orders)."""
    pairs = {}
    for ia, (la, ea, ca, ra, oa) in enumerate(shells):
        for ib in range(ia + 1):
            lb, eb, cb, rb, ob = shells[ib]
            ab2 = np.sum((ra - rb) ** 2)
            for a, c1 in zip(ea, ca):
                for b, c2 in zip(eb, cb):
                    p = a + b
                    P = (a * ra + b * rb) / p
                    c = c1 * c2 * np.exp(-a * b / p * ab2)
                    pairs.setdefault(la + lb, []).append(
                        dict(ia=ia, ib=ib, la=la, lb=lb, p=p, P=P, c=c, ea=a, eb=b,
                             PA=P - ra, PB=P - rb, oa=oa, ob=ob,
                             atA=atoms[ia], atB=atoms[ib], offd=0.0 if ia == ib else 1.0))
    return pairs


def md_j_grad(mol, dm, nbatch=1, sieve=True):
    """d/dX Tr(D J[D]), (natm, 3), for a symmetric dm, by the MD route."""
    shells, nao = shells_of(mol)
    atoms = atoms_of(mol)
    pairs = unordered_pairs(shells, atoms)
    # Ket stage, the energy engine's own: Dt = (-1)^|tau| c E (D + offd D^T).
    for lq, lst in pairs.items():
        sg = np.array([(-1.0) ** sum(h) for h in herm(lq)])
        hs = herm(lq)
        for pp in lst:
            la, lb = pp["la"], pp["lb"]
            env = _e_env(la, lb, pp)
            ca, cb = cart(la), cart(lb)
            dab = np.zeros((len(ca), len(cb)))
            for x, a in enumerate(ca):
                for y, b in enumerate(cb):
                    dab[x, y] = dm[pp["oa"] + x, pp["ob"] + y] \
                        + pp["offd"] * dm[pp["ob"] + y, pp["oa"] + x]
            g = np.zeros(len(hs))
            for k, h in enumerate(hs):
                for x, a in enumerate(ca):
                    for y, b in enumerate(cb):
                        g[k] += _e3(env, a, b, h) * dab[x, y]
            pp["dab"] = dab
            pp["env"] = env
            pp["Dt"] = pp["c"] * sg * g
    # Hermite stage with the bra raised to lp + 1.
    for lp, bra in pairs.items():
        hp1 = herm(lp + 1)
        for pp in bra:
            pp["Jt"] = np.zeros(len(hp1))
        for lq, ket in pairs.items():
            L = lp + 1 + lq
            p = np.array([x["p"] for x in bra])[:, None]
            q = np.array([x["p"] for x in ket])[None, :]
            PQ = (np.array([x["P"] for x in bra])[:, None, :]
                  - np.array([x["P"] for x in ket])[None, :, :])
            alpha = p * q / (p + q)
            pref = TWO_PI_2_5 / (p * q * np.sqrt(p + q))
            F = boys(L, alpha * np.sum(PQ ** 2, axis=2))
            Dt = np.array([x["Dt"] for x in ket])
            hq = herm(lq)
            Jt = np.zeros((len(bra), len(hp1)))
            for g in rtree.batch_dags(L, nbatch, sieve=sieve):
                env = {f"pq{DIRS[d]}": PQ[:, :, d] for d in range(3)}
                for m in g.leaves():
                    env[rtree.leaf_name(m)] = pref * (-2.0 * alpha) ** m * F[m]
                run(g.statements(), env)
                for (r, _) in g.targets:
                    val = env[rtree.leaf_name(0)] if r == (0, 0, 0) else env[rtree.r_name(r, 0)]
                    for i, t in enumerate(hp1):
                        tau = tuple(r[k] - t[k] for k in range(3))
                        if min(tau) < 0 or sum(tau) > lq:
                            continue
                        Jt[:, i] += val @ Dt[:, herm_index(lq)[tau]]
            for i, x in enumerate(bra):
                x["Jt"] = x["Jt"] + Jt[i]
    # Gradient bra.
    grad = np.zeros((mol.natm, 3))
    for lp, lst in pairs.items():
        hp1 = herm(lp + 1)
        idx1 = herm_index(lp + 1)
        hp = herm(lp)
        for pp in lst:
            la, lb = pp["la"], pp["lb"]
            envp = _e_env(la + 1, lb, pp)     # the a + 1, a, a - 1 tables at once
            ca, cb = cart(la), cart(lb)
            J = pp["Jt"]
            gA = np.zeros(3)
            gT = np.zeros(3)
            for x, a in enumerate(ca):
                for y, b in enumerate(cb):
                    w = pp["dab"][x, y] * pp["c"]
                    for d in range(3):
                        up = list(a); up[d] += 1
                        s = 0.0
                        for k, t in enumerate(hp1):
                            s += _e3(envp, up, b, t) * J[k]
                        acc = 2.0 * pp["ea"] * s
                        if a[d] > 0:
                            dn = list(a); dn[d] -= 1
                            s2 = 0.0
                            for k, t in enumerate(hp1):
                                s2 += _e3(envp, dn, b, t) * J[k]
                            acc -= a[d] * s2
                        gA[d] += w * acc
                        sh = 0.0
                        for t in hp:
                            t1 = list(t); t1[d] += 1
                            sh += _e3(envp, a, b, t) * J[idx1[tuple(t1)]]
                        gT[d] += w * sh
            grad[pp["atA"]] += 2.0 * gA
            grad[pp["atB"]] += 2.0 * (gT - gA)
    return grad


def pyscf_j_grad(mol, dm):
    """d/dX Tr(D J[D]) from int2e_ip1: -4 sum_{a on X} D_ab D_cd (nabla_a b|cd).

    libcint's ip1 is the derivative with respect to the electron coordinate,
    minus the one with respect to the centre."""
    ip1 = mol.intor("int2e_ip1_cart")            # (3, a, b, c, d), derivative on a
    rho = np.einsum("cd,xabcd->xab", dm, ip1)
    contr = np.einsum("ab,xab->xa", dm, rho)      # per function a
    g = np.zeros((mol.natm, 3))
    ao_atom = []
    for ib in range(mol.nbas):
        ao_atom.extend([mol.bas_atom(ib)] * ((mol.bas_angular(ib) + 1) * (mol.bas_angular(ib) + 2) // 2
                                             * mol.bas_nctr(ib)))
    for a, at in enumerate(ao_atom):
        g[at] -= 4.0 * contr[:, a]
    return g


def _mols():
    from pyscf import gto
    water = "O 0.0 0.1 -0.05; H 0.76 0.55 0.02; H -0.71 0.6 0.1"
    yield "water 6-31G* (d)", gto.M(atom=water, basis="6-31g*", cart=True)
    yield "water cc-pVDZ (general contraction)", gto.M(atom=water, basis="cc-pvdz", cart=True)
    yield "HF cc-pVTZ (f)", gto.M(atom="F 0.1 -0.05 0.02; H 0.2 0.3 0.88", basis="cc-pvtz", cart=True)


def main():
    rng = np.random.default_rng(11)
    worst = 0.0
    for name, mol in _mols():
        nao = mol.nao
        a = rng.standard_normal((nao, nao))
        dm = a + a.T
        ref = pyscf_j_grad(mol, dm)
        for nb, sv in ((1, False), (1, True), (3, True), (6, True)):
            g = md_j_grad(mol, dm, nbatch=nb, sieve=sv)
            err = np.max(np.abs(g - ref)) / np.max(np.abs(ref))
            tr = np.max(np.abs(g.sum(axis=0))) / np.max(np.abs(ref))
            worst = max(worst, err)
            print(f"{name:40s} nbatch={nb} sieve={sv!s:5s} max rel err {err:.2e}  sum/max {tr:.1e}")
    print(f"worst {worst:.2e}")
    return 0 if worst < 1e-11 else 1


if __name__ == "__main__":
    sys.exit(main())
