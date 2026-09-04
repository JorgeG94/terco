#!/usr/bin/env python3
"""RI-MP2 against pyscf's DF-MP2 on the same orbital and auxiliary sets.

scf_rimp2.F90 writes `rimp2_probe.bin`: the primitive shells of both bases as
built, the RHF energy and terco's E_os, E_ss. pyscf rebuilds both bases (cart,
same exponents), runs RHF to 1e-12 and DF-MP2 with that auxiliary basis, and
the three energies are compared. The SCF is converged independently on each
side, so agreement is at the SCF's tolerance, ~1e-10, not the solver's.
"""
import os
import numpy as np
from pyscf import gto, scf, mp, df
SYMBOL = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
here = os.path.dirname(os.path.abspath(__file__))
f = open(os.path.join(here, 'rimp2_probe.bin'), 'rb')
n, natm, nsh, nsha, nocc = (int(v) for v in np.fromfile(f, dtype=np.int64, count=5))
z = np.fromfile(f, dtype=np.float64, count=natm)
r = np.fromfile(f, dtype=np.float64, count=3*natm).reshape(3, natm, order='F').T
def shells(count):
    out = []
    for _ in range(count):
        a, l = np.fromfile(f, dtype=np.int64, count=2)
        e = np.fromfile(f, dtype=np.float64, count=1)[0]
        out.append((int(a), int(l), float(e)))
    return out
ao = shells(nsh)
ax = shells(nsha)
e_scf, e_os, e_ss = np.fromfile(f, dtype=np.float64, count=3)

atoms, bas, auxbas = [], {}, {}
for a in range(natm):
    lab = f'{SYMBOL[int(round(z[a]))]}{a + 1}'
    atoms.append([lab, tuple(r[a])])
    bas[lab] = [[l, [e, 1.0]] for (aa, l, e) in ao if aa == a + 1]
    auxbas[lab] = [[l, [e, 1.0]] for (aa, l, e) in ax if aa == a + 1]
mol = gto.M(atom=atoms, basis=bas, cart=True, unit='Bohr')
assert mol.nao_nr() == n, (mol.nao_nr(), n)
mf = scf.RHF(mol)
mf.conv_tol = 1e-12
mf.verbose = 0
mf.kernel()
pt = mp.dfmp2.DFMP2(mf)
pt.with_df = df.DF(mol, auxbasis=auxbas)
pt.verbose = 0
pt.kernel()
ref_os, ref_ss = pt.e_corr_os, pt.e_corr_ss
print(f'  E_RHF   terco {e_scf:20.12f}  pyscf {mf.e_tot:20.12f}  diff {abs(e_scf - mf.e_tot):9.2e}')
print(f'  E_os    terco {e_os:20.12f}  pyscf {ref_os:20.12f}  diff {abs(e_os - ref_os):9.2e}')
print(f'  E_ss    terco {e_ss:20.12f}  pyscf {ref_ss:20.12f}  diff {abs(e_ss - ref_ss):9.2e}')
print(f'  E_corr  terco {e_os + e_ss:20.12f}  pyscf {pt.e_corr:20.12f}  diff {abs(e_os + e_ss - pt.e_corr):9.2e}')
bad = abs(e_scf - mf.e_tot) > 1e-9 or abs(e_os - ref_os) > 1e-9 or abs(e_ss - ref_ss) > 1e-9
print('rimp2_ref: ' + ('FAIL' if bad else 'PASS'))
raise SystemExit(1 if bad else 0)
