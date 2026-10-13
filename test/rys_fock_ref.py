#!/usr/bin/env python3
"""eri%fock (vhf = J - K/2 for the total density, src/trc_fock.F90's own
header) against pyscf's raw int2e, on the SAME fixed, deterministic,
non-SCF density check_rys_fock_vs_pyscf.F90 builds and dumps -- the pyscf
side of the first comparison that can exercise an f shell (l=3).

Reuses rys_grad_ref.py's `bse_json_basis` (BSE JSON -> pyscf basis dict,
handling the combined-SP and general-contraction shapes the same way) and
`terco_to_pyscf_perm` (terco's literal shell order -> pyscf's cart=True
per-atom ascending-l order) rather than a second copy: both scripts read
the identical JSON convention and hit the identical pyscf AO-ordering
quirk, and a second hand-copy of either function is a second place for the
two to silently drift apart. See rys_grad_ref.py's own header for the
AO-order and AO-normalisation story in full; the short version repeated
here is just the formula, not the reasoning.

G/vhf is an ordinary <mu|operator|nu> matrix element, COVARIANT under a
change of AO normalisation exactly like S (S_terco = N S_pyscf N, the
overlap check below) -- NOT like D, which is CONTRAVARIANT (D_pyscf =
N D_terco N, same direction rys_grad_ref.py uses for the density). So
G_pyscf = G_terco / (n_i n_j), the OPPOSITE direction from D. See the
comment at this script's own transform for the measured bug that confusing
the two directions produces -- a clean, reproducible per-d-AO factor, not
noise, caught by comparing the invariant scalar 0.5*Tr(D G) before trusting
any individual element.

Usage: rys_fock_ref.py <basis.json> <dump.bin>
"""
import os
import sys

import numpy as np
from pyscf import gto

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from rys_grad_ref import bse_json_basis, terco_to_pyscf_perm, _SYMBOL  # noqa: E402

TOL = 1.0e-9
#: Same absolute gate as rys_grad_ref.py's overlap check and gradient
#: comparison (test/rys_grad_ref.py) -- the existing pyscf-comparison
#: tolerance in this tree, and this check must not loosen it.


def fock_rhf(eri, dm):
    """vhf = J - K/2 for the TOTAL density, raw (uncompressed) int2e --
    no 2x/4x factor to apply here, unlike rys_grad_ref.py's gradient
    formulas: J and K are each a plain contraction of the untouched
    4-index tensor, not a derivative recovered through a permutation
    symmetry."""
    j = np.einsum('ijkl,lk->ij', eri, dm)
    k = np.einsum('ijkl,jk->il', eri, dm)
    return j - 0.5*k


def main():
    if len(sys.argv) < 3:
        print('usage: rys_fock_ref.py <basis.json> <dump.bin>'); sys.exit(1)
    bpath = sys.argv[1]
    if not os.path.isabs(bpath):
        bpath = os.path.join(_HERE, bpath)

    # check_rys_fock_vs_pyscf.F90 itself skips (stop 0, no dump written)
    # when the basis needs more l than this build's TRC_LMAX reaches -- see
    # its own header. ctest's FIXTURES_REQUIRED still runs this half of the
    # pair in that case, so it has to recognise "no dump" as the same SKIP
    # rather than a missing-file error.
    if not os.path.exists(sys.argv[2]):
        print(f'  SKIP: {sys.argv[2]} was not written (the build-side check '
              f'skipped itself -- see its own stdout for the LMAX reason).')
        sys.exit(0)

    f = open(sys.argv[2], 'rb')
    nat, nao = (int(v) for v in np.fromfile(f, dtype=np.int64, count=2))
    at_z = np.fromfile(f, dtype=np.int64, count=nat)
    r = np.fromfile(f, dtype=np.float64, count=3*nat).reshape(3, nat, order='F').T
    d_terco = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    g_terco = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    s_terco = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    f.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(bpath, elements)
    atoms = [(_SYMBOL[int(z)], tuple(r[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)

    print(f'  nao   terco {nao}   pyscf {mol.nao_nr()}')
    if mol.nao_nr() != nao:
        print('  RESULT: FAIL (basis mismatch)'); sys.exit(1)

    # terco's literal shell order -> pyscf's cart=True per-atom ascending-l
    # order (rys_grad_ref.py's own header has the full story and the
    # 6-31G* oxygen example that first caught this).
    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)
    s_terco = s_terco[np.ix_(inv_perm, inv_perm)]
    d_terco = d_terco[np.ix_(inv_perm, inv_perm)]
    g_terco = g_terco[np.ix_(inv_perm, inv_perm)]

    # Per-AO normalisation: pyscf's cart=True d/f components are not all
    # unit-normalised the way this build's are (xy/xz/yz carry a different
    # factor from xx/yy/zz at l=2, and the ten Cartesian f components carry
    # three distinct factors -- xxx/yyy/zzz, xxy-type, xyz -- at l=3). The
    # overlap check below is what proves the BASIS agrees before trusting
    # any density- or Fock-matrix number that follows it.
    s_pyscf = mol.intor('int1e_ovlp')
    n = np.sqrt(np.diag(s_terco) / np.diag(s_pyscf))
    s_pred = n[:, None] * s_pyscf * n[None, :]
    s_err = np.abs(s_terco - s_pred).max()
    print(f'  overlap check: max|S_terco - N S_pyscf N| = {s_err:.3e}'
          f'  ({"OK" if s_err < TOL else "FAIL"})')
    if s_err >= TOL:
        i, j = np.unravel_index(np.argmax(np.abs(s_terco - s_pred)), s_terco.shape)
        print(f'    largest off-pattern element at ({i}, {j}): '
              f'S_terco={s_terco[i, j]!r}  N S_pyscf N={s_pred[i, j]!r}')
        print('  RESULT: FAIL (basis genuinely differs, not just AO normalisation)')
        sys.exit(1)

    # D is a CONTRAVARIANT object (density-matrix coefficients of a state
    # expansion, like C C^T) and transforms D_pyscf = N D_terco N -- same
    # direction rys_grad_ref.py uses. G/vhf is COVARIANT, an ordinary
    # <AO_i|operator|AO_j> matrix element exactly like S, and transforms the
    # SAME way S does (S_terco = N S_pyscf N, confirmed above) -- which runs
    # the OTHER direction from D: G_pyscf = G_terco / (n_i n_j). Multiplying
    # by N here instead of dividing by it reproduces the bug this comment
    # replaces: every element touching a d AO came out scaled by n_d**2 too
    # far (n_d**4 on a pure d-d element), measured at 0.397887/0.158314 on
    # water/cc-pVDZ -- a clean 5/(4*pi) = n_d**2, not noise, and resolved by
    # checking the SCALAR 0.5*Tr(D G) (basis-convention-invariant) agreed
    # between terco and pyscf to 1e-13 while the elementwise one did not.
    d_pyscf = n[:, None] * d_terco * n[None, :]
    g_terco_pyscf_ao = g_terco / n[:, None] / n[None, :]

    eri = mol.intor('int2e', aosym='s1').reshape(nao, nao, nao, nao)
    g_ref = fock_rhf(eri, d_pyscf)

    diff = np.abs(g_terco_pyscf_ao - g_ref)
    worst = diff.max()
    print(f'  worst |G_terco - G_pyscf| = {worst:.3e}'
          f'  ({"OK" if worst < TOL else "FAIL"})')
    if worst >= TOL:
        i, j = np.unravel_index(np.argmax(diff), diff.shape)
        print(f'    largest disagreement at ({i}, {j}): '
              f'terco={g_terco_pyscf_ao[i, j]!r}  pyscf={g_ref[i, j]!r}')
        print('  RESULT: FAIL'); sys.exit(1)
    print('  RESULT: PASS')


if __name__ == '__main__':
    main()
