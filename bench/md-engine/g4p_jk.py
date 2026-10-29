"""GPU4PySCF get_jk at the converged density, with its own per-class timer.

    python g4p_jk.py <xyz> <basis> [direct_scf_tol]      (default 1e-13; terco runs use 1e-12)

Converges RHF once (the density is cached as dm_<xyz>_<basis>.npy), times
two get_jk calls, then one more at DEBUG1 so GPU4PySCF prints its per-class
"(xx|xx) wall time" lines (microseconds) -- join_profile.py reads those.
Compare full builds: terco builds J and K from the whole density.
"""
import sys, os, time
import numpy as np, cupy as cp
import pyscf, gpu4pyscf
from pyscf.lib import logger
from gpu4pyscf import scf
xyz, basis = sys.argv[1], sys.argv[2]
tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-13
tag = os.path.basename(xyz).split('.')[0] + '_' + basis
mol = pyscf.M(atom=xyz, basis=basis, verbose=0)
print('gpu4pyscf', gpu4pyscf.__version__, 'nao', mol.nao, 'nbas', mol.nbas, 'direct_scf_tol', tol)
mf = scf.RHF(mol); mf.direct_scf_tol = tol; mf.conv_tol = 1e-9
dmf = f'dm_{tag}.npy'
if os.path.exists(dmf):
    dm = cp.asarray(np.load(dmf))
else:
    t = time.time(); e = mf.kernel(); print('E', e, 'scf', f'{time.time()-t:.1f} s', 'cycles', mf.cycles)
    dm = mf.make_rdm1(); np.save(dmf, dm.get())
vj, vk = mf.get_jk(mol, dm)  # warm-up
for k in range(2):
    cp.cuda.Device().synchronize(); t = time.time()
    vj, vk = mf.get_jk(mol, dm)
    cp.cuda.Device().synchronize(); print(f'get_jk {k}: {time.time()-t:.3f} s')
mol.verbose = logger.DEBUG1; mf.verbose = logger.DEBUG1
vj, vk = mf.get_jk(mol, dm)
cp.cuda.Device().synchronize()
