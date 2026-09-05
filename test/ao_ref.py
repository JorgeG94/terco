#!/usr/bin/env python3
"""Basis function values and gradients against pyscf's eval_ao.

check_ao.F90 writes `ao_probe.bin`: the system is the uncontracted one from
mult_ref.py, so pyscf renormalises nothing and the two codes are over the
same functions or one of them is wrong.
"""

import numpy as np
from pyscf import gto
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

import os
f = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ao_probe.bin'), 'rb')
npts, n = np.fromfile(f, dtype=np.int64, count=2)
pts = np.fromfile(f, dtype=np.float64, count=3*npts).reshape(3, npts, order='F').T
chi = np.fromfile(f, dtype=np.float64, count=npts*n).reshape(npts, n, order='F')
gchi = np.fromfile(f, dtype=np.float64, count=3*npts*n).reshape(3, npts, n, order='F')
f.close()
print(f'  nao   terco {n}   pyscf {nao}')
if n != nao:
    raise SystemExit(1)

ao = numint.eval_ao(mol, np.ascontiguousarray(pts), deriv=1)   # (4, npts, nao)
scale = np.abs(ao[0]).max()
d0 = np.abs(chi - ao[0]).max() / scale
print(f'  values     max|diff|/max {d0:.3e}    {"OK" if d0 < 1e-12 else "MISMATCH"}')
bad = d0 >= 1e-12
for k, name in enumerate('xyz'):
    s = max(np.abs(ao[1 + k]).max(), 1.0)
    d = np.abs(gchi[k] - ao[1 + k]).max() / s
    print(f'  d/d{name}       max|diff|/max {d:.3e}    {"OK" if d < 1e-12 else "MISMATCH"}')
    bad = bad or d >= 1e-12
if bad:
    # Localise: which AO is off, in values
    per = np.abs(chi - ao[0]).max(axis=0) / scale
    print('  per-AO value error:', np.array2string(per, precision=2))
    raise SystemExit(1)
print('  RESULT: PASS')
