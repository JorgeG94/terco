#!/usr/bin/env python3
"""Restricted Kohn-Sham total energies against pyscf, on terco's grid.

scf_rks.F90 writes `rks_probe.bin`: the grid it integrated over, the
shells exactly as it built them (water, 6-31G uncontracted to primitives),
and for each functional the converged energy and density. The same basis
is rebuilt here in the same order, pyscf's RKS is run on the same points
and weights from terco's converged density, and the energies are compared.
"""

import os
import numpy as np
from pyscf import gto, dft

SYMBOL = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
NAMES = ['svwn', 'pbe', 'blyp', 'b3lyp']
PYSCF = {'svwn': 'LDA_X,LDA_C_VWN', 'pbe': 'PBE', 'blyp': 'BLYP', 'b3lyp': 'HYB_GGA_XC_B3LYP'}

f = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rks_probe.bin'), 'rb')
npts, n, nfunc, natm, nsh = np.fromfile(f, dtype=np.int64, count=5)
pts = np.fromfile(f, dtype=np.float64, count=3*npts).reshape(3, npts, order='F').T
w = np.fromfile(f, dtype=np.float64, count=npts)
z = np.fromfile(f, dtype=np.float64, count=natm)
r = np.fromfile(f, dtype=np.float64, count=3*natm).reshape(3, natm, order='F').T
shells = []
for _ in range(nsh):
    a, l = np.fromfile(f, dtype=np.int64, count=2)
    e = np.fromfile(f, dtype=np.float64, count=1)[0]
    shells.append((int(a), int(l), float(e)))
res = []
for k in range(nfunc):
    e = np.fromfile(f, dtype=np.float64, count=1)[0]
    dm = np.fromfile(f, dtype=np.float64, count=n*n).reshape(n, n, order='F')
    res.append((e, dm))
f.close()

# One label per atom so each carries its own shell list, in terco's order.
atoms, bas = [], {}
for a in range(natm):
    lab = f'{SYMBOL[int(round(z[a]))]}{a + 1}'
    atoms.append([lab, tuple(r[a])])
    bas[lab] = [[l, [e, 1.0]] for (aa, l, e) in shells if aa == a + 1]
mol = gto.M(atom=atoms, basis=bas, cart=True, unit='Bohr', charge=0, spin=0)
nao = mol.nao_nr()
print(f'  nao   terco {n}   pyscf {nao}      points {npts}   shells {nsh}')
if n != nao:
    raise SystemExit(1)

bad = 0
for k, name in enumerate(NAMES[:nfunc]):
    e_t, dm_t = res[k]
    mf = dft.RKS(mol)
    mf.xc = PYSCF[name]
    mf.grids.coords = np.ascontiguousarray(pts)
    mf.grids.weights = w
    mf.grids.non0tab = None
    mf.conv_tol = 1e-12
    mf.max_cycle = 200
    mf.verbose = 0
    e_p = mf.kernel(dm0=dm_t)
    de = abs(e_t - e_p)
    ok = mf.converged and de < 1e-8
    bad += not ok
    print(f'  {name:6s} E terco {e_t:20.12f}  pyscf {e_p:20.12f}  |dE| {de:.2e}  '
          f'{"OK" if ok else "MISMATCH"}{"" if mf.converged else " (pyscf not converged)"}')
print(f'\n  functionals disagreeing: {bad} / {nfunc}')
if bad:
    raise SystemExit(1)
print('  RESULT: PASS')
