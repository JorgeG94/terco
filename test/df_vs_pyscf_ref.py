#!/usr/bin/env python3
"""trc_df_2c ((P|Q)) and trc_df_3c ((mu nu|P)) against pyscf's int2c2e and
int3c2e (pyscf.df.incore.aux_e2), same geometry, orbital and auxiliary
bases read from the same JSON files check_df_vs_pyscf.F90 reads -- the
density-fitting analogue of rys_fock_ref.py, and the first comparison
exercising a g shell (l=4) on the auxiliary side.

Reuses rys_grad_ref.py's bse_json_basis/terco_to_pyscf_perm (see its own
header for the AO-order story) ONCE PER BASIS: the orbital and auxiliary
sets are two independent pyscf Mole objects, each with its own shell
order and its own permutation.

NORMALISATION. (mu nu|P) and (P|Q) are ordinary integrals, each index an
AO-like slot of its own basis, and every index transforms COVARIANTLY --
the same direction as the overlap matrix S in rys_fock_ref.py's own
derivation (S_terco = N S_pyscf N), not the contravariant direction a
density matrix uses. n is a property of the represented basis function,
not of which operator the integral happens to be, so it does not need to
come from an overlap specifically: check_df_vs_pyscf.F90 derives n_i from
trc_df_2c applied to the ORBITAL basis against itself (the (i|j)
two-centre Coulomb metric, diagonal i=j) rather than trc_1e, whose own
ceiling (ONE_E_LMAX) is independent of DF_LMAX and need not have been
raised alongside it -- see that program's own header. The pyscf side
matches: `mol.intor('int2c2e')`, not `int1e_ovlp`. With n_i from that
metric and n_P from the auxiliary (P|Q) diagonal:

    J_pyscf_pred(P,Q) = n_P n_Q J_pyscf(P,Q)   compared to J_terco directly
    tens_terco_pyscf_ao(i,j,P) = tens_terco(i,j,P) / (n_i n_j n_P)

compared to pyscf's raw aux_e2. Dividing by n (rather than multiplying,
D's direction) is the same bug this project already paid for once -- see
rys_fock_ref.py's own comment on the wrong-direction bug it replaced.

Usage: df_vs_pyscf_ref.py <basis.json> <aux.json> <dump.bin>
"""
import os
import sys

import numpy as np
from pyscf import df, gto

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from rys_grad_ref import bse_json_basis, terco_to_pyscf_perm, _SYMBOL  # noqa: E402

TOL = 1.0e-9
#: Not looser than the existing pyscf comparisons (rys_grad_ref.py,
#: rys_fock_ref.py both gate at 1e-9 absolute).


def main():
    if len(sys.argv) < 4:
        print('usage: df_vs_pyscf_ref.py <basis.json> <aux.json> <dump.bin>')
        sys.exit(1)
    bpath, apath = sys.argv[1], sys.argv[2]
    if not os.path.isabs(bpath):
        bpath = os.path.join(_HERE, bpath)
    if not os.path.isabs(apath):
        apath = os.path.join(_HERE, apath)

    # check_df_vs_pyscf.F90 skips itself (stop 0, no dump written) when
    # either basis needs more l than this build's DF_LMAX/DF_LMAX_AUX
    # reach -- see its own header. ctest's FIXTURES_REQUIRED still runs
    # this half of the pair in that case, so a missing dump is the same
    # SKIP, not a missing-file error (same convention as rys_fock_ref.py).
    if not os.path.exists(sys.argv[3]):
        print(f'  SKIP: {sys.argv[3]} was not written (the build-side check '
              f'skipped itself -- see its own stdout for the DF_LMAX reason).')
        sys.exit(0)

    f = open(sys.argv[3], 'rb')
    nat, nao, naux = (int(v) for v in np.fromfile(f, dtype=np.int64, count=3))
    at_z = np.fromfile(f, dtype=np.int64, count=nat)
    r = np.fromfile(f, dtype=np.float64, count=3*nat).reshape(3, nat, order='F').T
    o_terco = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    j_terco = np.fromfile(f, dtype=np.float64, count=naux*naux).reshape(naux, naux, order='F')
    tens_terco = np.fromfile(f, dtype=np.float64,
                             count=nao*nao*naux).reshape(nao, nao, naux, order='F')
    f.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(bpath, elements)
    aux_arg = bse_json_basis(apath, elements)
    atoms = [(_SYMBOL[int(z)], tuple(r[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)
    auxmol = gto.M(atom=atoms, basis=aux_arg, unit='Bohr', cart=True, verbose=0)

    print(f'  nao    terco {nao}   pyscf {mol.nao_nr()}')
    print(f'  naux   terco {naux}  pyscf {auxmol.nao_nr()}')
    if mol.nao_nr() != nao or auxmol.nao_nr() != naux:
        print('  RESULT: FAIL (basis mismatch)'); sys.exit(1)

    # Two independent permutations: the orbital and auxiliary bases are two
    # different Mole objects, each stable-sorted by l per atom on its own
    # (see terco_to_pyscf_perm's own header).
    perm_orb = terco_to_pyscf_perm(atoms, basis_arg)
    perm_aux = terco_to_pyscf_perm(atoms, aux_arg)
    inv_orb = np.argsort(perm_orb)
    inv_aux = np.argsort(perm_aux)
    o_terco = o_terco[np.ix_(inv_orb, inv_orb)]
    j_terco = j_terco[np.ix_(inv_aux, inv_aux)]
    tens_terco = tens_terco[np.ix_(inv_orb, inv_orb, inv_aux)]

    # Per-AO normalisation on the orbital side, from the (i|j) two-centre
    # Coulomb metric over the orbital basis against itself -- see this
    # script's own header for why not the overlap.
    o_pyscf = mol.intor('int2c2e')
    n_orb = np.sqrt(np.diag(o_terco) / np.diag(o_pyscf))
    o_pred = n_orb[:, None] * o_pyscf * n_orb[None, :]
    o_err = np.abs(o_terco - o_pred).max()
    print(f'  orbital (i|j) metric check: max|O_terco - N O_pyscf N| = {o_err:.3e}'
          f'  ({"OK" if o_err < TOL else "FAIL"})')
    if o_err >= TOL:
        print('  RESULT: FAIL (orbital basis genuinely differs)'); sys.exit(1)

    # Per-AO normalisation on the auxiliary side, same construction but
    # from the (P|Q) metric's diagonal -- (P|P) plays the same role (i|i)
    # does above: both are a diagonal, strictly positive, covariant
    # quantity.
    j_pyscf = auxmol.intor('int2c2e')
    n_aux = np.sqrt(np.diag(j_terco) / np.diag(j_pyscf))
    j_pred = n_aux[:, None] * j_pyscf * n_aux[None, :]
    j_err = np.abs(j_terco - j_pred).max()
    print(f'  2c (P|Q) check: max|J_terco - N J_pyscf N| = {j_err:.3e}'
          f'  ({"OK" if j_err < TOL else "FAIL"})')
    if j_err >= TOL:
        i, j = np.unravel_index(np.argmax(np.abs(j_terco - j_pred)), j_terco.shape)
        print(f'    largest disagreement at ({i}, {j}): '
              f'terco={j_terco[i, j]!r}  N J_pyscf N={j_pred[i, j]!r}')
        print('  RESULT: FAIL'); sys.exit(1)

    # 3c (mu nu|P): covariant in every index, same direction as S and J
    # above -- divide the terco tensor by n_i n_j n_P to land in pyscf's
    # AO convention, then compare directly to the raw aux_e2 tensor.
    tens_pyscf = df.incore.aux_e2(mol, auxmol, 'int3c2e')
    tens_terco_pyscf_ao = tens_terco / n_orb[:, None, None] / n_orb[None, :, None] \
        / n_aux[None, None, :]
    diff = np.abs(tens_terco_pyscf_ao - tens_pyscf)
    worst = diff.max()
    print(f'  3c (mu nu|P) check: max|tens_terco - tens_pyscf| = {worst:.3e}'
          f'  ({"OK" if worst < TOL else "FAIL"})')
    if worst >= TOL:
        i, j, k = np.unravel_index(np.argmax(diff), diff.shape)
        print(f'    largest disagreement at ({i}, {j}, {k}): '
              f'terco={tens_terco_pyscf_ao[i, j, k]!r}  pyscf={tens_pyscf[i, j, k]!r}')
        print('  RESULT: FAIL'); sys.exit(1)
    print('  RESULT: PASS')


if __name__ == '__main__':
    main()
