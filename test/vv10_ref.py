#!/usr/bin/env python3
"""VV10 (trc_vv10) against pyscf's numint.nr_nlc_vxc, on the same grid,
weights and density check_vv10_vs_pyscf.F90 dumps to vv10_probe.bin.

b = 6.0, C = 0.01 (wB97X-V's own NLC parameters) are not passed as bare
numbers to pyscf -- they are fetched through the SAME lookup the task's
own confirmation used, `dft.libxc.nlc_coeff('WB97X_V')` -- so a drift in
either code's constant would show up as a mismatch here rather than two
copies of the same literal silently agreeing. `numint.nr_nlc_vxc(xc_code)`
only ever uses `xc_code` to call `nlc_coeff` -- it never evaluates a
semilocal functional -- so passing 'WB97X_V' through it runs the VV10
kernel ALONE, with exactly wB97X-V's constants, nothing else mixed in.

AO ORDER/NORMALISATION: reuses rys_grad_ref.py's `bse_json_basis` and
`terco_to_pyscf_perm` (pyscf's cart=True groups a given atom's shells by
ascending l; terco keeps the JSON's literal order) and the same per-AO
rescaling check_rys_grad_vs_pyscf/rys_fock_ref use for a Cartesian d's
differently-normalised off-diagonal components: n_i =
sqrt(S_terco(i,i)/S_pyscf(i,i)).

V_xc IS LIKE THE OVERLAP, NOT LIKE THE DENSITY, under that rescaling: it is
a gradient dE/dD, and D_pyscf = N D_terco N (rys_grad_ref.py's own
identity) gives, by the chain rule through that same linear map,
V_terco = N V_pyscf N -- multiply BY n_i n_j to go from pyscf's matrix to
terco's, the same direction D goes from terco to pyscf, not the inverse
direction rys_fock_ref.py's plain operator matrix G needs (G_pyscf =
G_terco/(n_i n_j), because G is paired with D by a trace and D already
carries one factor of N on each side). Checked here against the
basis-independent scalar E_nlc as well, which needs no N at all.

The ADD contract: `vmat` in the dump is SEEDED with a nonzero matrix
before trc_vv10_rks/uks is called, so it is `vmat_after - vmat_seed`
(not `vmat_after`) that must match pyscf's nlc matrix -- proving the
module adds its contribution rather than overwriting whatever was there.
"""
import os
import sys

import numpy as np
from pyscf import gto, dft
from pyscf.dft import numint, libxc

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from rys_grad_ref import bse_json_basis, terco_to_pyscf_perm, _SYMBOL  # noqa: E402

TOL = 1.0e-9
XC_NLC = 'WB97X_V'

B_REF, C_REF = libxc.nlc_coeff(XC_NLC)[0][0]
print(f'  nlc_coeff({XC_NLC!r}) = (b={B_REF}, C={C_REF})')
assert abs(B_REF - 6.0) < 1e-12 and abs(C_REF - 0.01) < 1e-12, \
    'wB97X-V NLC parameters drifted from (6.0, 0.01)'


def read_matrix(f, n, m=None):
    m = n if m is None else m
    return np.fromfile(f, dtype=np.float64, count=n*m).reshape(n, m, order='F')


def run_case(f, kind, basjson):
    nat, nao, npts = (int(v) for v in np.fromfile(f, dtype=np.int64, count=3))
    at_z = np.fromfile(f, dtype=np.int64, count=nat)
    r = np.fromfile(f, dtype=np.float64, count=3*nat).reshape(3, nat, order='F').T
    coords = np.fromfile(f, dtype=np.float64, count=3*npts).reshape(3, npts, order='F').T
    w = np.fromfile(f, dtype=np.float64, count=npts)
    s_terco = read_matrix(f, nao)
    if kind == 0:
        dmat_terco = read_matrix(f, nao)
    else:
        da_terco = read_matrix(f, nao)
        db_terco = read_matrix(f, nao)
        dmat_terco = da_terco + db_terco
    vmat0 = read_matrix(f, nao) if kind == 0 else (read_matrix(f, nao), read_matrix(f, nao))
    e_nlc_terco = np.fromfile(f, dtype=np.float64, count=1)[0]
    if kind == 0:
        vmat_terco = read_matrix(f, nao)
    else:
        va_terco = read_matrix(f, nao)
        vb_terco = read_matrix(f, nao)

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    bpath = os.path.join(_HERE, basjson)
    basis_arg = bse_json_basis(bpath, elements)
    atoms = [(_SYMBOL[int(z)], tuple(r[i])) for i, z in enumerate(at_z)]
    spin = 1 if kind == 1 else 0   # OH radical: Nalpha - Nbeta = 1
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0, spin=spin)

    label = f'{"RKS" if kind == 0 else "UKS"} nao={nao} pts={npts}'
    print(f'  {label}: pyscf nao {mol.nao_nr()}')
    if mol.nao_nr() != nao:
        print('  RESULT: FAIL (basis mismatch)')
        return False

    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)

    s_terco_p = s_terco[np.ix_(inv_perm, inv_perm)]
    dmat_terco_p = dmat_terco[np.ix_(inv_perm, inv_perm)]

    s_pyscf = mol.intor('int1e_ovlp')
    n = np.sqrt(np.diag(s_terco_p) / np.diag(s_pyscf))
    s_pred = n[:, None] * s_pyscf * n[None, :]
    s_err = np.abs(s_terco_p - s_pred).max()
    print(f'    overlap check: max|S_terco - N S_pyscf N| = {s_err:.3e}'
          f'  ({"OK" if s_err < 1e-9 else "FAIL"})')
    if s_err >= 1e-9:
        print('  RESULT: FAIL (basis genuinely differs, not just AO normalisation)')
        return False

    dmat_pyscf = n[:, None] * dmat_terco_p * n[None, :]

    grids = dft.gen_grid.Grids(mol)
    grids.coords = np.ascontiguousarray(coords)
    grids.weights = w
    grids.non0tab = None
    ni = numint.NumInt()
    nelec_p, exc_p, vmat_p = ni.nr_nlc_vxc(mol, grids, XC_NLC, dmat_pyscf)

    de = abs(e_nlc_terco - exc_p)
    print(f'    E_nlc: terco {e_nlc_terco:18.14f}  pyscf {exc_p:18.14f}  |dE| {de:.2e}'
          f'  ({"OK" if de < TOL else "FAIL"})')

    v_pred_p = n[:, None] * vmat_p * n[None, :]   # predicted terco V, pyscf AO order

    ok = de < TOL
    if kind == 0:
        dv = (vmat_terco - vmat0)[np.ix_(inv_perm, inv_perm)]
        err = np.abs(dv - v_pred_p).max()
        print(f'    max|dV_nlc| = {err:.3e}  ({"OK" if err < TOL else "FAIL"})')
        ok = ok and err < TOL
    else:
        dva = (va_terco - vmat0[0])[np.ix_(inv_perm, inv_perm)]
        dvb = (vb_terco - vmat0[1])[np.ix_(inv_perm, inv_perm)]
        erra = np.abs(dva - v_pred_p).max()
        errb = np.abs(dvb - v_pred_p).max()
        print(f'    max|dV_nlc,a| = {erra:.3e}  max|dV_nlc,b| = {errb:.3e}'
              f'  ({"OK" if erra < TOL and errb < TOL else "FAIL"})')
        ok = ok and erra < TOL and errb < TOL
    print(f'  RESULT: {"PASS" if ok else "FAIL"}')
    return ok


def main():
    f = open(os.path.join(_HERE, 'vv10_probe.bin'), 'rb')
    ncase = int(np.fromfile(f, dtype=np.int64, count=1)[0])
    bases = ['basis_sets/6-31gs.json', 'basis_sets/cc-pvdz.json', 'basis_sets/cc-pvdz.json']
    bad = 0
    for k in range(ncase):
        kind = int(np.fromfile(f, dtype=np.int64, count=1)[0])
        if not run_case(f, kind, bases[k]):
            bad += 1
    f.close()
    print(f'\n  cases failing: {bad} / {ncase}')
    if bad:
        sys.exit(1)
    print('  OVERALL: PASS')


if __name__ == '__main__':
    main()
