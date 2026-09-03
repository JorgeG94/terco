#!/usr/bin/env python3
"""SCF total energies against pyscf, case by case, on terco's grid.

scf_rks.F90 writes `rks_probe.bin`: for each case the shells exactly as
built (6-31G uncontracted to primitives), the electron counts, the grid
when there is a functional, the converged energy and density. The same
basis is rebuilt here in the same order, pyscf's RHF/UHF/RKS/UKS is run
on the same points from terco's density, and the energies are compared.
"""

import os
import numpy as np
from pyscf import gto, scf, dft

SYMBOL = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
PYSCF = {'svwn': 'LDA_X,LDA_C_VWN', 'pbe': 'PBE', 'blyp': 'BLYP', 'b3lyp': 'HYB_GGA_XC_B3LYP',
         'b3lyp5': 'HYB_GGA_XC_B3LYP5', 'lda_x': 'LDA_X', 'pbe0': 'PBE0'}

f = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rks_probe.bin'), 'rb')
ncase = int(np.fromfile(f, dtype=np.int64, count=1)[0])
bad = 0
for ic in range(ncase):
    n, natm, nsh, nspin, na, nb = (int(v) for v in np.fromfile(f, dtype=np.int64, count=6))
    z = np.fromfile(f, dtype=np.float64, count=natm)
    r = np.fromfile(f, dtype=np.float64, count=3*natm).reshape(3, natm, order='F').T
    shells = []
    for _ in range(nsh):
        a, l = np.fromfile(f, dtype=np.int64, count=2)
        e = np.fromfile(f, dtype=np.float64, count=1)[0]
        shells.append((int(a), int(l), float(e)))
    name = f.read(8).decode().strip()
    npts = int(np.fromfile(f, dtype=np.int64, count=1)[0])
    if npts:
        pts = np.fromfile(f, dtype=np.float64, count=3*npts).reshape(3, npts, order='F').T
        w = np.fromfile(f, dtype=np.float64, count=npts)
    e_t = np.fromfile(f, dtype=np.float64, count=1)[0]
    dm = np.fromfile(f, dtype=np.float64, count=n*n*nspin).reshape(n, n, nspin, order='F')

    atoms, bas = [], {}
    for a in range(natm):
        lab = f'{SYMBOL[int(round(z[a]))]}{a + 1}'
        atoms.append([lab, tuple(r[a])])
        bas[lab] = [[l, [e, 1.0]] for (aa, l, e) in shells if aa == a + 1]
    charge = int(round(z.sum())) - (na + nb)
    mol = gto.M(atom=atoms, basis=bas, cart=True, unit='Bohr', charge=charge, spin=na - nb)
    if mol.nao_nr() != n:
        print(f'  case {ic}: nao terco {n} pyscf {mol.nao_nr()}')
        raise SystemExit(1)

    # pyscf groups an atom's shells by angular momentum; terco keeps the order
    # given. Permute terco's AO order into pyscf's through (atom, l, exponent).
    off, k = {}, 0
    for (aa, l, e) in shells:
        off[(aa, l, round(e, 10))] = k
        k += (l + 1)*(l + 2)//2
    perm = []
    for ib in range(mol.nbas):
        key = (mol.bas_atom(ib) + 1, mol.bas_angular(ib), round(float(mol.bas_exp(ib)[0]), 10))
        perm.extend(range(off[key], off[key] + mol.bas_len_cart(ib)))
    perm = np.array(perm)
    dm = dm[perm][:, perm]

    if name:
        mf = dft.RKS(mol) if nspin == 1 else dft.UKS(mol)
        mf.xc = PYSCF[name]
        mf.grids.coords = np.ascontiguousarray(pts)
        mf.grids.weights = w
        mf.grids.non0tab = None
    else:
        mf = scf.RHF(mol) if nspin == 1 else scf.UHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 200
    mf.verbose = 0
    dm0 = dm[:, :, 0] if nspin == 1 else (dm[:, :, 0], dm[:, :, 1])
    # The control: pyscf's energy AT terco's density, before any iteration.
    # A reconverged energy can agree while the density was wrong; this cannot.
    e_at = mf.energy_tot(dm0)
    e_p = mf.kernel(dm0=dm0)
    de = abs(e_t - e_p)
    dat = abs(e_t - e_at)
    ok = mf.converged and de < 1e-8 and dat < 1e-8
    bad += not ok
    label = (name or 'hf') + ('/U' if nspin == 2 else '/R')
    print(f'  {label:9s} nel {na}+{nb}  E terco {e_t:18.12f}  pyscf {e_p:18.12f}  |dE| {de:.1e}  '
          f'at terco D {dat:.1e}  {"OK" if ok else "MISMATCH"}{"" if mf.converged else " (pyscf not converged)"}')
f.close()
print(f'\n  cases disagreeing: {bad} / {ncase}')
if bad:
    raise SystemExit(1)
print('  RESULT: PASS')
