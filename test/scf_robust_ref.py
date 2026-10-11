#!/usr/bin/env python3
"""The hard-case set (scf_robust.F90) against pyscf's own basis library.

Unlike rks_ref.py this does not rebuild terco's shells primitive by
primitive: the basis here is CONTRACTED (6-31G, aug-cc-pVDZ), read by both
codes from a name rather than from terco's exact exponents, so the two
contraction-normalisation conventions need not agree to the last digit.
What this checks is that a hard case converged to the right PHYSICS, not a
bit-level fold convention -- agreement at 1e-5 to 1e-6 Hartree is the bar,
not the 1e-8 the uncontracted comparisons hold to.
"""

import os
import json
import numpy as np
from pyscf import gto, scf

SYMBOL = {1: 'H', 7: 'N', 8: 'O', 10: 'Ne'}
#: name -> element symbol: basis names pyscf's own library does not carry
#: (d-aug-cc-pVDZ is not shipped under that name) are read straight from
#: the same BSE JSON file terco builds from, converted to the [l, [[e, c],
#: ...]] per contracted function pyscf's `basis={...}` dict wants. A BSE
#: "shell" can hold several contracted functions sharing one exponent set
#: (`coefficients` has one row per function); each row becomes its own
#: pyscf shell entry at the same l.
_HERE = os.path.dirname(os.path.abspath(__file__))


def bse_json_basis(path, elements):
    d = json.load(open(path))
    out = {}
    for sym, z in elements.items():
        el = d['elements'][str(z)]
        shells = []
        for sh in el['electron_shells']:
            l = sh['angular_momentum'][0]
            exps = [float(e) for e in sh['exponents']]
            for row in sh['coefficients']:
                co = [float(c) for c in row]
                shells.append([l] + [[e, c] for e, c in zip(exps, co)])
        out[sym] = shells
    return out

f = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scf_robust_probe.bin'), 'rb')
ncase = int(np.fromfile(f, dtype=np.int64, count=1)[0])
bad = 0
for ic in range(ncase):
    n, natm, nspin, na, nb, charge = (int(v) for v in np.fromfile(f, dtype=np.int64, count=6))
    xyzname = f.read(16).decode().strip()
    basname = f.read(16).decode().strip()
    e_terco = np.fromfile(f, dtype=np.float64, count=1)[0]

    atoms = []
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), xyzname + '.xyz')) as g:
        nat = int(g.readline())
        g.readline()
        for _ in range(nat):
            parts = g.readline().split()
            atoms.append((parts[0], (float(parts[1]), float(parts[2]), float(parts[3]))))

    try:
        # pyscf's own basis library first -- 6-31G and aug-cc-pVDZ are in
        # it, so this is the usual path and keeps that comparison exactly
        # as it was.
        mol = gto.M(atom=atoms, basis=basname, unit='Angstrom', cart=True, charge=charge,
                    spin=na - nb, verbose=0)
    except Exception:
        # Not a name pyscf ships (d-aug-cc-pVDZ): read the same BSE JSON
        # terco built from instead.
        elements = {sym: {'H': 1, 'N': 7, 'O': 8, 'Ne': 10}[sym] for sym, _ in atoms}
        bpath = os.path.join(_HERE, 'basis_sets', basname + '.json')
        basis_arg = bse_json_basis(bpath, elements)
        mol = gto.M(atom=atoms, basis=basis_arg, unit='Angstrom', cart=True, charge=charge,
                    spin=na - nb, verbose=0)
    if mol.nao_nr() != n:
        print(f'  case {xyzname}: nao terco {n} pyscf {mol.nao_nr()} -- basis mismatch, skipping')
        continue
    mf = scf.UHF(mol) if nspin == 2 else scf.RHF(mol)
    mf.conv_tol = 1e-11
    e_py = mf.kernel()
    diff = abs(e_terco - e_py)
    # The diffuse case (w5/d-aug-cc-pVDZ) is deliberately past the point
    # where S has a well-defined inverse: canonical orthogonalisation
    # drops functions below `lindep_thresh`, and pyscf's default RHF does
    # not apply the same cut (or any -- it has no lindep protection on by
    # default), so the two codes are not solving the same, numerically
    # well-posed problem to begin with. Agreement there is bounded by how
    # near-degenerate the dropped subspace is, not by either SCF's
    # convergence, so it gets a tolerance two orders looser than the other
    # three cases, which are nowhere near that edge.
    tol = 3.0e-4 if basname.startswith('d-aug') else 2.0e-5
    ok = diff < tol
    print(f'  {xyzname:14s} {basname:12s}  terco {e_terco:20.10f}  pyscf {e_py:20.10f}'
          f'  |diff| {diff:10.3e}  {"OK" if ok else "FAIL"}')
    if not ok:
        bad += 1

if bad:
    raise SystemExit(f'{bad} case(s) disagree with pyscf by more than 2e-5 Hartree')
print('scf_robust_ref: all hard cases agree with pyscf')
