import numpy as np
from pyscf import gto

# EXACTLY the system in test/check_mult.F90: three centres, uncontracted
# s/p/d each. Uncontracted is the whole point -- pyscf renormalises a
# CONTRACTED shell to unit self-overlap and terco does not, which is what
# made the first comparison (on reference system 2) disagree on the overlap
# itself. One primitive per shell leaves nothing to renormalise.
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

s = mol.intor('int1e_ovlp_cart')
r1 = mol.intor('int1e_r_cart').reshape(3, nao, nao)
r2 = mol.intor('int1e_rr_cart').reshape(9, nao, nao)
r3 = mol.intor('int1e_rrr_cart').reshape(27, nao, nao)
ref = np.concatenate([r1, r2, r3], axis=0)          # 39, nao, nao

import os
f = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mult_probe.bin'), 'rb')
n, nm = np.fromfile(f, dtype=np.int64, count=2)
ts = np.fromfile(f, dtype=np.float64, count=n*n).reshape(n, n, order='F')
tm = np.fromfile(f, dtype=np.float64, count=n*n*nm).reshape(n, n, nm, order='F')
f.close()

print(f'  nao   terco {n}   pyscf {nao}')
ds = np.abs(ts - s).max()
print(f'  CONTROL overlap   max|diff| {ds:.3e}    {"OK" if ds < 1e-12 else "MISMATCH"}')
if ds >= 1e-12:
    print('  -> systems differ; nothing below means anything')
    raise SystemExit(1)

names = ['x', 'y', 'z'] + [a+b for a in 'xyz' for b in 'xyz'] \
        + [a+b+c for a in 'xyz' for b in 'xyz' for c in 'xyz']
bad = 0
for k in range(nm):
    d = np.abs(tm[:, :, k] - ref[k]).max()
    scale = max(np.abs(ref[k]).max(), 1.0)
    if d > 1e-11*scale:
        bad += 1
        if bad <= 12:
            print(f'  {k:3d} {names[k]:<4s} max|diff| {d:.3e}   '
                  f'terco sum {tm[:,:,k].sum():18.10f}  pyscf sum {ref[k].sum():18.10f}')
print(f'\n  components disagreeing: {bad} / {nm}')
if bad == 0:
    print('  RESULT: PASS')
