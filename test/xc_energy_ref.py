#!/usr/bin/env python3
"""E_xc and V_xc against libxc through pyscf, on terco's grid and density.

check_xc_energy.F90 writes `xc_probe.bin`. This feeds the same points,
weights and density matrix to pyscf's numint for each functional and
compares the energy and every element of the potential matrix.
"""

import os
import numpy as np
from pyscf import gto, dft
from pyscf.dft import numint

coords = [(0.0, 0.0, 0.0), (1.3, 0.0, 0.7), (-0.4, 1.1, 1.9)]
labels = ['Be', 'Li', 'He']
bas, atoms = {}, []
for i, (lab, c) in enumerate(zip(labels, coords), start=1):
    bas[lab] = [[0, [0.9 + 0.17*i, 1.0]],
                [1, [1.4 + 0.11*i, 1.0]],
                [2, [1.9 + 0.13*i, 1.0]]]
    atoms.append([lab, c])
mol = gto.M(atom=atoms, basis=bas, cart=True, unit='Bohr', spin=None)
nao = mol.nao_nr()

# terco name -> pyscf name. b3lyp is libxc's (VWN_RPA); b3lyp5 is VWN5.
NAMES = ['lda_x', 'svwn', 'pbe', 'blyp', 'b3lyp', 'b3lyp5']
PYSCF = {'lda_x': 'LDA_X', 'svwn': 'LDA_X,LDA_C_VWN', 'pbe': 'PBE', 'blyp': 'BLYP',
         'b3lyp': 'HYB_GGA_XC_B3LYP', 'b3lyp5': 'HYB_GGA_XC_B3LYP5'}

f = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'xc_probe.bin'), 'rb')
npts, n, nfunc = np.fromfile(f, dtype=np.int64, count=3)
pts = np.fromfile(f, dtype=np.float64, count=3*npts).reshape(3, npts, order='F').T
w = np.fromfile(f, dtype=np.float64, count=npts)
dm = np.fromfile(f, dtype=np.float64, count=n*n).reshape(n, n, order='F')
results = []
for k in range(nfunc):
    exc, nelec = np.fromfile(f, dtype=np.float64, count=2)
    vxc = np.fromfile(f, dtype=np.float64, count=n*n).reshape(n, n, order='F')
    results.append((exc, nelec, vxc))
# The spin-polarised section: two densities, then per functional.
dm2 = np.fromfile(f, dtype=np.float64, count=2*n*n).reshape(n, n, 2, order='F')
uks = []
for k in range(nfunc):
    exc, nelec = np.fromfile(f, dtype=np.float64, count=2)
    vxc = np.fromfile(f, dtype=np.float64, count=2*n*n).reshape(n, n, 2, order='F')
    uks.append((exc, nelec, vxc))
f.close()
print(f'  nao   terco {n}   pyscf {nao}      points {npts}')
if n != nao:
    raise SystemExit(1)

grids = dft.gen_grid.Grids(mol)
grids.coords = np.ascontiguousarray(pts)
grids.weights = w
grids.non0tab = None
ni = numint.NumInt()

bad = 0
for k, name in enumerate(NAMES[:nfunc]):
    exc_t, nelec_t, vxc_t = results[k]
    nelec_p, exc_p, vxc_p = ni.nr_rks(mol, grids, PYSCF[name], dm)
    de = abs(exc_t - exc_p)
    dn = abs(nelec_t - nelec_p)
    dv = np.abs(vxc_t - vxc_p).max()
    ok = de < 1e-9 and dn < 1e-9 and dv < 1e-9
    bad += not ok
    print(f'  {name:7s} E_xc terco {exc_t:18.12f} pyscf {exc_p:18.12f}  |dE| {de:.2e}  '
          f'|dN| {dn:.2e}  max|dV| {dv:.2e}  {"OK" if ok else "MISMATCH"}')
print('  spin-polarised:')
for k, name in enumerate(NAMES[:nfunc]):
    exc_t, nelec_t, vxc_t = uks[k]
    nelec_p, exc_p, vxc_p = ni.nr_uks(mol, grids, PYSCF[name], (dm2[:, :, 0], dm2[:, :, 1]))
    de = abs(exc_t - exc_p)
    dn = abs(nelec_t - nelec_p.sum())
    dv = max(np.abs(vxc_t[:, :, 0] - vxc_p[0]).max(), np.abs(vxc_t[:, :, 1] - vxc_p[1]).max())
    ok = de < 1e-9 and dn < 1e-9 and dv < 1e-9
    bad += not ok
    print(f'  {name:7s} E_xc terco {exc_t:18.12f} pyscf {exc_p:18.12f}  |dE| {de:.2e}  '
          f'|dN| {dn:.2e}  max|dV| {dv:.2e}  {"OK" if ok else "MISMATCH"}')
print(f'\n  functionals disagreeing: {bad} / {2*nfunc}')
if bad:
    raise SystemExit(1)
print('  RESULT: PASS')
