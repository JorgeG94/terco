#!/usr/bin/env python3
"""eri%gradient against pyscf's int2e_ip1, same basis, same density.

check_rys_grad_vs_pyscf.F90 writes the dump named on both command lines: water, a fixed
synthetic density (not a converged SCF one -- this is a kernel check on an
arbitrary density, same reasoning as check_rys_grad_fd.F90), RHF (one
density) and UHF (Da != Db).

The reference STARTED from pyscf's documented formula (pyscf/grad/rhf.py's
get_jk/get_veff/grad_elec, pyscf/grad/uhf.py's get_veff) -- not pyscf's
Gradients class, which assumes a CONVERGED SCF density (it adds the
energy-weighted overlap term) and so does not apply to an arbitrary fixed
D -- but the "2x" those use is specific to the symmetry-COMPRESSED (s2kl)
integral storage `_vhf.direct_mapdm` contracts against; this uses the
RAW, uncompressed `int2e_ip1` tensor (every ijkl, no s2kl folding), which
needs the full "4x" from the 8-fold permutation symmetry instead (one
factor from the bra bi<->ket, already present in vj - 0.5 vk; a second
from i<->j, since only i is differentiated and j is recovered by the
bra's own exchange symmetry, the same way terco's kernel gets its B
derivative from an AB-shift instead of a second raised 2D table). Checked
against a finite difference computed independently in pyscf itself
(perturb mol.atom, rebuild, recompute the plain int2e contraction) before
trusting it against terco -- 2x silently passed a dimensional sniff test
and was wrong by exactly half, every component, every system:

  RHF:  vj_ij = -sum_kl (nabla_i i j|k l) D_lk
        vk_il = -sum_jk (nabla_i i j|k l) D_jk
        vhf = vj - 0.5 vk
        dE2/dR_A = 4 * sum_{i in A, j} vhf_ij D_ij

  UHF:  vj_tot from D_tot = Da + Db, vk_a from Da, vk_b from Db
        veff_a = vj_tot - vk_a ; veff_b = vj_tot - vk_b
        dE2/dR_A = 4 * sum_{i in A, j} (veff_a_ij Da_ij + veff_b_ij Db_ij)

Done with a raw `int2e_ip1` contraction (no _vhf screening driver -- water
is tiny) so the formula is visible and checkable line by line, not hidden
behind pyscf's optimised C path.

STATUS: PASSES, both bases, to ~1e-13 absolute against a 1e-9 gate --
after two independent bugs in THIS SCRIPT, neither in the gradient kernel:

  1. AO ORDER. terco keeps every shell in the literal order the JSON's
     `electron_shells` lists them. pyscf's cart=True `Mole.build()` does
     NOT preserve that order: per atom, it stable-sorts shells by
     ascending angular momentum first -- checked directly against
     `mol.ao_labels()`, not assumed: 6-31G* oxygen's literal shell list
     is s, s, p, s, p, d (the two "SP" shells split apart for pyscf's
     basis dict), but its AOs come back grouped 1s 2s 3s / 2p 3p / 3d.
     AO index i meant two different functions in the two matrices, so
     every comparison before `terco_to_pyscf_perm` was comparing the
     wrong pairs -- which is exactly why the mismatch was not a clean
     scale factor (an earlier version of this status section said so and
     stopped there; the permutation is what the "basis genuinely
     differs" branch below would have caught if it had been tried first,
     since an S that disagrees by a permutation, not a rescaling, fails
     the overlap check with the wrong SIGN on some elements, not just the
     wrong magnitude -- the tell that sent this back to `mol.ao_labels()`
     instead of deeper into trc_basis_json.F90).
  2. AO NORMALISATION (the coordinator's diagnosis, confirmed once (1)
     stopped hiding it): an SCF energy is invariant to rescaling the AOs;
     a fixed synthetic density is not. `n_i = sqrt(S_terco(i,i) /
     S_pyscf(i,i))`; the overlap check asserts `S_terco == N S_pyscf N`
     to 1e-9 before anything else runs, so a future convention change on
     either side fails loudly here instead of silently in the gradient
     numbers.

Independently confirmed before (1) and (2) were found, and still true:
  - terco's eri%gradient agrees with a finite difference of terco's own
    eri%fock to ~1e-7 relative, RHF and UHF, both bases
    (check_rys_grad_fd.F90) -- the kernel's own derivative algebra was
    right against terco's own integrals the whole time.
  - This script's grad2e_rhf/grad2e_uhf agree with a finite difference
    computed entirely inside pyscf (plain `int2e`, no terco involved) to
    ~1e-7 relative, same two systems -- the reference formula was right
    against pyscf's own integrals the whole time too.
"""
import json
import os
import sys

import numpy as np
from pyscf import gto

_HERE = os.path.dirname(os.path.abspath(__file__))
_SYMBOL = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
           14: 'Si', 15: 'P', 16: 'S', 17: 'Cl'}


def bse_json_basis(path, elements):
    #: BSE JSON's electron_shells -> pyscf's basis={...} dict. Two shapes
    #: share this field, and they are not the same thing:
    #:  - len(angular_momentum) > 1: a COMBINED shell (6-31G-family "SP",
    #:    sharing one exponent set between an s and a p function). One row
    #:    of `coefficients` per l, in the SAME order as `angular_momentum`
    #:    -- row 0 is l=angular_momentum[0], row 1 is angular_momentum[1].
    #:  - len(angular_momentum) == 1: a GENERAL CONTRACTION at that one l
    #:    (cc-pVDZ's oxygen s-block: 9 exponents, 3 rows, all l=0). Every
    #:    row is that SAME l, not matched up positionally against a
    #:    length-1 list.
    #: scf_robust_ref.py's own bse_json_basis takes angular_momentum[0]
    #: unconditionally for every row -- right for the second shape, wrong
    #: for the first (every row becomes an s-shell); a first attempt here
    #: zipped angular_momentum against coefficients instead -- right for
    #: the first shape, wrong for the second (drops every row past the
    #: length of angular_momentum, e.g. 2 of cc-pVDZ oxygen's 3 s-functions
    #: and one of its two H contractions). Neither was caught by nao alone
    #: matching for the OTHER shape's basis until both were tried.
    d = json.load(open(path))
    out = {}
    for sym, z in elements.items():
        el = d['elements'][str(z)]
        shells = []
        for sh in el['electron_shells']:
            ls = sh['angular_momentum']
            exps = [float(e) for e in sh['exponents']]
            rows = sh['coefficients']
            if len(ls) > 1:
                pairs = zip(ls, rows)
            else:
                pairs = ((ls[0], row) for row in rows)
            for l, row in pairs:
                co = [float(c) for c in row]
                shells.append([l] + [[e, c] for e, c in zip(exps, co)])
        out[sym] = shells
    return out


_NCART = {0: 1, 1: 3, 2: 6, 3: 10}


def terco_to_pyscf_perm(atoms, basis_arg):
    """terco keeps every basis function in the LITERAL order its shell list
    gives them (confirmed against trc_basis_json.F90's reader: one pass
    over `electron_shells`, in file order, no regrouping). pyscf's cart=True
    `Mole.build()` does NOT preserve that order -- per atom, it stable-sorts
    shells by ascending angular momentum first (checked directly against
    `mol.ao_labels()`: for 6-31G* oxygen, the literal shell list is
    s, s, p, s, p, d -- two SP shells split apart -- but pyscf's AOs come
    back grouped 1s 2s 3s / 2p 3p / 3d). Returns `perm` with
    `perm[i] = j` meaning terco's AO slot i is pyscf's AO slot j, so
    `pyscf_mat = terco_mat[np.ix_(perm, perm)]` (or `terco_vec[perm]` for
    one axis) puts a terco-ordered array into pyscf's order.
    """
    perm = []
    base = 0
    for sym, _ in atoms:
        shells = basis_arg[sym]  # [(l, [[e, c], ...]), ...], literal order
        ls = [row[0] for row in shells]
        n_at = sum(_NCART[l] for l in ls)
        # terco's own AO offset of each shell, literal order
        terco_off = []
        off = 0
        for l in ls:
            terco_off.append(off)
            off += _NCART[l]
        # pyscf's AO offset of each shell: stable sort by l
        order = sorted(range(len(ls)), key=lambda k: ls[k])
        pyscf_off = [0]*len(ls)
        off = 0
        for k in order:
            pyscf_off[k] = off
            off += _NCART[ls[k]]
        for k, l in enumerate(ls):
            for c in range(_NCART[l]):
                perm.append(base + pyscf_off[k] + c)
        base += n_at
    return np.array(perm, dtype=np.int64)


def grad2e_rhf(eri1, aoslice, natm, dm):
    vj = -np.einsum('xijkl,lk->xij', eri1, dm)
    vk = -np.einsum('xijkl,jk->xil', eri1, dm)
    vhf = vj - 0.5*vk
    de = np.zeros((natm, 3))
    for ia in range(natm):
        p0, p1 = aoslice[ia, 2:]
        de[ia] = np.einsum('xij,ij->x', vhf[:, p0:p1], dm[p0:p1])*4.0
    return de


def grad2e_uhf(eri1, aoslice, natm, da, db):
    dtot = da + db
    vj_tot = -np.einsum('xijkl,lk->xij', eri1, dtot)
    vk_a = -np.einsum('xijkl,jk->xil', eri1, da)
    vk_b = -np.einsum('xijkl,jk->xil', eri1, db)
    veff_a = vj_tot - vk_a
    veff_b = vj_tot - vk_b
    de = np.zeros((natm, 3))
    for ia in range(natm):
        p0, p1 = aoslice[ia, 2:]
        de[ia] = (np.einsum('xij,ij->x', veff_a[:, p0:p1], da[p0:p1]) +
                  np.einsum('xij,ij->x', veff_b[:, p0:p1], db[p0:p1]))*4.0
    return de


def main():
    if len(sys.argv) < 3:
        print('usage: rys_grad_ref.py <basis.json> <dump.bin>'); sys.exit(1)
    bpath = sys.argv[1]
    if not os.path.isabs(bpath):
        bpath = os.path.join(_HERE, bpath)

    f = open(sys.argv[2], 'rb')
    nat, nao = (int(v) for v in np.fromfile(f, dtype=np.int64, count=2))
    at_z = np.fromfile(f, dtype=np.int64, count=nat)
    r = np.fromfile(f, dtype=np.float64, count=3*nat).reshape(3, nat, order='F').T
    drhf = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    da = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    db = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    grad_rhf_terco = np.fromfile(f, dtype=np.float64, count=3*nat).reshape(3, nat, order='F').T
    grad_uhf_terco = np.fromfile(f, dtype=np.float64, count=3*nat).reshape(3, nat, order='F').T
    s_terco = np.fromfile(f, dtype=np.float64, count=nao*nao).reshape(nao, nao, order='F')
    f.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(bpath, elements)
    atoms = [(_SYMBOL[int(z)], tuple(r[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)

    print(f'  nao   terco {nao}   pyscf {mol.nao_nr()}')
    if mol.nao_nr() != nao:
        print('  RESULT: FAIL (basis mismatch)'); sys.exit(1)

    # terco keeps every shell in the JSON's own literal order; pyscf's
    # cart=True Mole.build() groups a given atom's shells by ascending l
    # instead (checked directly against mol.ao_labels() -- 6-31G*'s split
    # "SP" shells come back as 1s 2s 3s / 2p 3p / 3d, not the s,s,p,s,p,d
    # the JSON actually lists). Put terco's arrays into pyscf's AO order
    # before anything else, or AO index i means two different functions in
    # the two matrices and every check below is comparing the wrong pairs.
    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)
    s_terco = s_terco[np.ix_(inv_perm, inv_perm)]
    drhf = drhf[np.ix_(inv_perm, inv_perm)]
    da = da[np.ix_(inv_perm, inv_perm)]
    db = db[np.ix_(inv_perm, inv_perm)]

    # Coordinator's diagnosis: an SCF energy is invariant to rescaling the
    # AOs, but a FIXED synthetic density is not -- pyscf's cart=True d
    # functions are not all unit-normalised the way this build's are
    # (xy/xz/yz carry a different factor from xx/yy/zz), so the basis can
    # be the SAME physical functions while S_terco and S_pyscf still differ
    # by a per-AO diagonal rescaling. n_i = sqrt(S_terco(i,i)/S_pyscf(i,i));
    # if the two overlaps agree up to N = diag(n), the basis is fine and
    # only the density's normalisation convention needs transforming before
    # contracting with pyscf's integrals: D_pyscf = N D_terco N.
    s_pyscf = mol.intor('int1e_ovlp')
    n = np.sqrt(np.diag(s_terco) / np.diag(s_pyscf))
    s_pred = n[:, None] * s_pyscf * n[None, :]
    s_err = np.abs(s_terco - s_pred).max()
    print(f'  overlap check: max|S_terco - N S_pyscf N| = {s_err:.3e}'
          f'  ({"OK" if s_err < 1e-9 else "FAIL"})')
    if s_err >= 1e-9:
        i, j = np.unravel_index(np.argmax(np.abs(s_terco - s_pred)), s_terco.shape)
        print(f'    largest off-pattern element at ({i}, {j}): '
              f'S_terco={s_terco[i, j]!r}  N S_pyscf N={s_pred[i, j]!r}')
        print('  RESULT: FAIL (basis genuinely differs, not just AO normalisation)')
        sys.exit(1)

    drhf_p = n[:, None] * drhf * n[None, :]
    da_p = n[:, None] * da * n[None, :]
    db_p = n[:, None] * db * n[None, :]

    eri1 = mol.intor('int2e_ip1', comp=3).reshape(3, nao, nao, nao, nao)
    aoslice = mol.aoslice_by_atom()

    grad_rhf_py = grad2e_rhf(eri1, aoslice, nat, drhf_p)
    grad_uhf_py = grad2e_uhf(eri1, aoslice, nat, da_p, db_p)

    ok = True
    for name, terco, py in (('RHF', grad_rhf_terco, grad_rhf_py),
                            ('UHF', grad_uhf_terco, grad_uhf_py)):
        d = np.abs(terco - py)
        worst = d.max()
        print(f'  {name}: worst |terco - pyscf| = {worst:.3e}'
              f'  ({"OK" if worst < 1e-9 else "FAIL"})')
        if worst >= 1e-9:
            ok = False
            print('    terco :\n', terco)
            print('    pyscf :\n', py)

    if not ok:
        print('  RESULT: FAIL'); sys.exit(1)
    print('  RESULT: PASS')


if __name__ == '__main__':
    main()
