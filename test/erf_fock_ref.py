#!/usr/bin/env python3
"""eri%fock / eri%fock_nosym with the erf-attenuated operator
(erf(omega*r12)/r12, omega = 0.3 bohr^-1, see check_erf_fock_vs_pyscf.F90's
own header for how J and K are isolated from the single g = j*J - k*K/2
entry point) against pyscf's own range-separated J/K.

Reuses rys_grad_ref.py's `bse_json_basis`/`terco_to_pyscf_perm`/`_SYMBOL`
and rys_fock_ref.py's AO-order/AO-normalisation handling verbatim -- same
JSON convention, same pyscf cart=True AO-ordering quirk, same covariant
transform for a Fock-like matrix (G_pyscf = G_terco / (n_i n_j), the
opposite direction from the density: see rys_fock_ref.py's own header for
the measured bug that confusing the two produces).

Usage: erf_fock_ref.py <basis.json> <dump.bin>
"""
import os
import sys

import numpy as np
from pyscf import gto, scf

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from rys_grad_ref import bse_json_basis, terco_to_pyscf_perm, _SYMBOL  # noqa: E402

TOL = 1.0e-9
#: Same absolute gate as every other pyscf comparison in this tree; never
#: loosened.
OMEGA = 0.3
#: Must match check_erf_fock_vs_pyscf.F90's OMEGA_TEST.


def main():
    if len(sys.argv) < 3:
        print('usage: erf_fock_ref.py <basis.json> <dump.bin>'); sys.exit(1)
    bpath = sys.argv[1]
    if not os.path.isabs(bpath):
        bpath = os.path.join(_HERE, bpath)

    # check_erf_fock_vs_pyscf.F90 skips itself (stop 0, no dump written) on a
    # basis this build's TRC_LMAX cannot reach yet -- same convention as
    # rys_fock_ref.py.
    if not os.path.exists(sys.argv[2]):
        print(f'  SKIP: {sys.argv[2]} was not written (the build-side check '
              f'skipped itself -- see its own stdout for the LMAX reason).')
        sys.exit(0)

    f = open(sys.argv[2], 'rb')
    nat, nao = (int(v) for v in np.fromfile(f, dtype=np.int64, count=2))
    at_z = np.fromfile(f, dtype=np.int64, count=nat)
    r = np.fromfile(f, dtype=np.float64, count=3*nat).reshape(3, nat, order='F').T
    d_terco = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    j_terco = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    k_terco = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    s_terco = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    f.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(bpath, elements)
    atoms = [(_SYMBOL[int(z)], tuple(r[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)

    print(f'  nao   terco {nao}   pyscf {mol.nao_nr()}')
    if mol.nao_nr() != nao:
        print('  RESULT: FAIL (basis mismatch)'); sys.exit(1)

    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)
    s_terco = s_terco[np.ix_(inv_perm, inv_perm)]
    d_terco = d_terco[np.ix_(inv_perm, inv_perm)]
    j_terco = j_terco[np.ix_(inv_perm, inv_perm)]
    k_terco = k_terco[np.ix_(inv_perm, inv_perm)]

    s_pyscf = mol.intor('int1e_ovlp')
    n = np.sqrt(np.diag(s_terco) / np.diag(s_pyscf))
    s_pred = n[:, None] * s_pyscf * n[None, :]
    s_err = np.abs(s_terco - s_pred).max()
    print(f'  overlap check: max|S_terco - N S_pyscf N| = {s_err:.3e}'
          f'  ({"OK" if s_err < TOL else "FAIL"})')
    if s_err >= TOL:
        print('  RESULT: FAIL (basis genuinely differs, not just AO normalisation)')
        sys.exit(1)

    # D contravariant (D_pyscf = N D_terco N); J/K covariant, same direction
    # as S (G_pyscf = G_terco / (n_i n_j)) -- see rys_fock_ref.py's header.
    d_pyscf = n[:, None] * d_terco * n[None, :]
    j_terco_ao = j_terco / n[:, None] / n[None, :]
    k_terco_ao = k_terco / n[:, None] / n[None, :]

    with mol.with_range_coulomb(OMEGA):
        vj_ref, vk_ref = scf.hf.get_jk(mol, d_pyscf)

    ok = True
    for name, terco_mat, ref_mat in (('J', j_terco_ao, vj_ref), ('K', k_terco_ao, vk_ref)):
        diff = np.abs(terco_mat - ref_mat)
        worst = diff.max()
        print(f'  worst |{name}_terco - {name}_pyscf| (omega={OMEGA}) = {worst:.3e}'
              f'  ({"OK" if worst < TOL else "FAIL"})')
        if worst >= TOL:
            i, j = np.unravel_index(np.argmax(diff), diff.shape)
            print(f'    largest disagreement at ({i}, {j}): '
                  f'terco={terco_mat[i, j]!r}  pyscf={ref_mat[i, j]!r}')
            ok = False
    print('  RESULT: ' + ('PASS' if ok else 'FAIL'))
    if not ok:
        sys.exit(1)


if __name__ == '__main__':
    main()
