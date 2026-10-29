#!/usr/bin/env python3
"""Reference Hirshfeld charges for test/check_hirshfeld.F90's two cases.

mqc solves its free atoms at run time and caches them in memory -- it has no
stored reference numbers for water/6-31G or O2/cc-pVDZ (its own
test/test_mqc_hirshfeld.f90 pins water/cc-pVDZ, a different basis). This
script runs mqc's OWN cross-check (tools/hirshfeld/hirshfeld_reference.py's
recipe, imported from tools/spherical_atom/spherical_atom_reference.py,
both in a metalquicha checkout) directly in pyscf for exactly the two cases
test/check_hirshfeld.F90 needs:

  * water, RHF/6-31G, at mqc's own validation/inputs/sample_inputs/w1.xyz
    (test/hirshfeld_water.xyz here -- NOT test/water.xyz, a different
    fragment);
  * O2 triplet, UHF/cc-pVDZ, at mqc's validation/inputs/sample_inputs/o2.xyz
    (identical to terco's own test/o2.xyz).

Both bases are read from mqc's own basis_sets/*.json (BSE layout) rather
than pyscf's internal tables, so the AO contraction matches what
trc_basis_json reads for terco. The molecular integral runs on pyscf's own
(150, 974) atom_grid, well beyond terco's grid level 5 default, and
evaluates every free atom exactly from its own orbitals rather than from an
interpolated table -- so agreement with terco checks the table AND terco's
own (coarser) quadrature, not just the partition formula.

    python3 tools/hirshfeld/gen_test_references.py --mqc-root /path/to/metalquicha

Needs pyscf. Prints the numbers test/check_hirshfeld.F90's Q_WATER/Q_O2
parameters are copied from (see that file's own header for the numbers this
produced on 2026-10-29, and the CHARGE_TOL comment for how tight the
agreement with terco actually measured).
"""
import argparse
import pathlib
import sys

import numpy as np
from pyscf import dft, gto, scf

MOLECULAR_GRID = (150, 974)
ANGSTROM_PER_BOHR = 0.529177210903  # CODATA 2018


def load_solver(mqc_root):
    sys.path.insert(0, str(pathlib.Path(mqc_root) / "tools" / "spherical_atom"))
    import spherical_atom_reference as sar  # noqa
    import json

    def bse_to_pyscf_any_order(name, z):
        # mqc's own bse_to_pyscf asserts its shells are l-sorted, a check
        # that only matters for comparing AO order against mqc's own
        # matrices; this script never does that, and some BSE sets (6-31G
        # for O, Z=8) are not l-sorted.
        data = json.loads((sar.ROOT / "basis_sets" / (name + ".json")).read_text())
        shells = []
        for sh in data["elements"][str(z)]["electron_shells"]:
            exps = [float(e) for e in sh["exponents"]]
            cols = [[float(c) for c in col] for col in sh["coefficients"]]
            moms = sh["angular_momentum"]
            if len(moms) == 1:
                shells.append([moms[0]] + [[e] + [c[i] for c in cols]
                                           for i, e in enumerate(exps)])
            else:
                for l, col in zip(moms, cols):
                    shells.append([l] + [[e, col[i]] for i, e in enumerate(exps)])
        return shells

    sar.bse_to_pyscf = bse_to_pyscf_any_order
    return sar


FREE_ATOM_BASIS = "def2-tzvp"
FREE_ATOM_FUNCTIONAL = "wb97x"

_free_atoms = {}


def free_atom(sar, symbol):
    if symbol not in _free_atoms:
        mol, mf, dm, _ = sar.solve(symbol, FREE_ATOM_BASIS, FREE_ATOM_FUNCTIONAL, True)
        _free_atoms[symbol] = (mol, dm[0] + dm[1])
    return _free_atoms[symbol]


def molecule(sar, atoms_bohr, basis, charge=0, spin=0):
    symbols = sorted({a[0] for a in atoms_bohr})
    return gto.M(atom=atoms_bohr, unit="Bohr", charge=charge, spin=spin,
                 basis={s: sar.bse_to_pyscf(basis, gto.charge(s)) for s in symbols},
                 cart=False, verbose=0)


def converge(mol):
    mf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.run()
    assert mf.converged
    dm = np.asarray(mf.make_rdm1())
    if dm.ndim == 2:
        return mf, dm, None
    return mf, dm[0] + dm[1], dm[0] - dm[1]


def populations(sar, mol, total, spin=None):
    grids = dft.Grids(mol)
    grids.atom_grid = MOLECULAR_GRID
    grids.build()
    ni = dft.numint.NumInt()
    ao = ni.eval_ao(mol, grids.coords)
    rho = ni.eval_rho(mol, ao, total)
    spin_rho = None if spin is None else ni.eval_rho(mol, ao, spin)

    free = np.empty((mol.natm, len(grids.weights)))
    for a in range(mol.natm):
        amol, adm = free_atom(sar, mol.atom_symbol(a))
        shifted = grids.coords - mol.atom_coord(a)
        free[a] = ni.eval_rho(amol, ni.eval_ao(amol, shifted), adm)
    total_free = free.sum(axis=0)
    share = np.divide(free, total_free, out=np.zeros_like(free), where=total_free > 0)

    electrons = share @ (grids.weights * rho)
    charges = mol.atom_charges() - electrons
    spins = None if spin_rho is None else share @ (grids.weights * spin_rho)
    return charges, spins


def read_xyz_bohr(path):
    lines = pathlib.Path(path).read_text().splitlines()
    n = int(lines[0].strip())
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        xyz_bohr = [float(v) / ANGSTROM_PER_BOHR for v in parts[1:4]]
        atoms.append((parts[0], tuple(xyz_bohr)))
    return atoms


def show(label, values):
    print("  %-10s" % label + " ".join("%17.12f" % v for v in values))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mqc-root", default=str(pathlib.Path.home() / "dev" / "metalquicha"))
    args = ap.parse_args()
    sar = load_solver(args.mqc_root)

    root = pathlib.Path(__file__).resolve().parents[2]
    water = read_xyz_bohr(root / "test" / "hirshfeld_water.xyz")
    o2 = read_xyz_bohr(root / "test" / "o2.xyz")

    mol = molecule(sar, water, "6-31g")
    mf, dm, _ = converge(mol)
    q, _ = populations(sar, mol, dm)
    print("water_6-31g RHF  E=%.12f" % mf.e_tot)
    show("charges", q)
    print("  sum        %.3e" % q.sum())

    mol2 = molecule(sar, o2, "cc-pvdz", charge=0, spin=2)
    mf2, dm2, sdm2 = converge(mol2)
    q2, s2 = populations(sar, mol2, dm2, sdm2)
    print("o2_cc-pvdz UHF triplet  E=%.12f" % mf2.e_tot)
    show("charges", q2)
    show("spins", s2)
    print("  sum        %.3e" % q2.sum())


if __name__ == "__main__":
    main()
