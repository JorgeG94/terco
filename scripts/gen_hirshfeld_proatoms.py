#!/usr/bin/env python3
"""Generate test/basis_sets/hirshfeld_proatoms.json: Hirshfeld free-atom tables.

terco's Hirshfeld charges (src/scf/trc_hirshfeld.F90) partition the molecular
density by the stockholder recipe: each atom's share of space is its own
free, spherically averaged density divided by the sum of every atom's. The
free atoms are ONE FIXED RECIPE regardless of the molecule: the spherically
averaged, spin-polarized (Hund's rule within each angular momentum) atom in
wB97X/def2-TZVP, following metalquicha (mqc)'s mqc_czt_hirshfeld.f90 and
mqc_czt_spherical_atom.f90 exactly -- this is a straight port of mqc's
recipe, same owner, mqc_ -> trc_ renamed, which is explicitly wanted.

WHY PYSCF AND NOT MQC'S OWN FORTRAN BUILD
------------------------------------------
mqc solves the free atom itself, at run time, through its CPU Hartree-Fock/
DFT backend (`-DMQC_ENABLE_CZT=ON -DMQC_ENABLE_LIBXC=ON`), and caches the
result in memory for the life of the process -- it ships no precomputed
table. Reproducing that exactly would mean fetching and building libxc
(mqc's own cmake/modules/Findlibxc.cmake measures ~970s for that alone) and
then mqc's whole CPU ab initio path, just to run 36 single-atom SCFs.

mqc's OWN cross-check is pyscf: tools/spherical_atom/spherical_atom_reference.py
implements the identical spherically-averaged, Hund's-rule, spin-polarized
solver directly in pyscf (AtomSphericalUHF/UKS, overriding `eig`/`get_occ` on
pyscf's `atom_hf.AtomSphericAverageRHF` machinery), and
mqc_czt_hirshfeld.f90's own header records that mqc's Fortran result at its
production grid level (7) agrees with this pyscf recipe on a 200x1202 grid to
2e-11 hartree. So this script calls that same pyscf solver directly rather
than building mqc's Fortran/libxc stack -- it is mqc's own validation
reference, used as the source instead of as a check on the side.

The basis is read from mqc's own basis_sets/*.json (Basis Set Exchange
layout) via a local copy of spherical_atom_reference.bse_to_pyscf, with its
l-sorted assertion dropped: that assertion only guards an AO-ORDER
comparison against mqc's own matrices (irrelevant here -- only the
spherically symmetric electron density is read back, which cannot depend on
AO order), and some BSE sets (def2-tzvp for Al, Z=13) are not l-sorted.

WHAT IS STORED
---------------
Each element's free-atom density, uniform in ln r from R_MIN to R_MAX Bohr
(N_RADIAL points), as ln(rho) and the log-log slope r (drho/dr)/rho -- the
pair terco's cubic Hermite interpolator (trc_hirshfeld.F90's
`hirshfeld_table_density`) needs, exactly as mqc's own
`mqc_czt_hirshfeld.f90::tabulate_free_atom` builds it. The slope is read off
the EXACT analytic AO gradient (pyscf's GTOval_sph_deriv1), not a finite
difference: every sample point lies on the z-axis and the atom is
spherically symmetric about the origin, so d rho/dz at that point IS d
rho/dr exactly.

USAGE
-----
    python3 scripts/gen_hirshfeld_proatoms.py --mqc-root /path/to/metalquicha

Needs pyscf. Takes several minutes (the heaviest transition metals, solved
spin-polarized in a fairly large basis, are the slow part -- about 10-30s
each; the light elements are under 5s). Writes
test/basis_sets/hirshfeld_proatoms.json, picked up automatically by the
existing `install(DIRECTORY test/basis_sets/ ...)` glob in app/CMakeLists.txt
(it matches "*.json") and by trc_hirshfeld_data's own loader
(test/basis_sets/hirshfeld_proatoms.json relative default,
TERCO_HIRSHFELD_DATA environment override), the same pattern trc_minao and
trc_sap already use for minao.json/sap_helfem_large.json.
"""
import argparse
import json
import pathlib
import subprocess
import sys
import time

import numpy as np

Z_MAX = 36
N_RADIAL = 2000
R_MIN = 1.0e-6
R_MAX = 60.0
RHO_FLOOR = 1.0e-30
BASIS = "def2-tzvp"
FUNCTIONAL = "wb97x"

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "test" / "basis_sets" / "hirshfeld_proatoms.json"


def mqc_commit(mqc_root):
    try:
        return subprocess.check_output(
            ["git", "-C", str(mqc_root), "rev-parse", "HEAD"],
            text=True).strip()
    except Exception:
        return "unknown"


def load_solver(mqc_root):
    """Import mqc's spherical_atom_reference, with its l-sorted assertion
    dropped (see this file's module docstring)."""
    sys.path.insert(0, str(pathlib.Path(mqc_root) / "tools" / "spherical_atom"))
    import spherical_atom_reference as sar  # noqa

    def bse_to_pyscf_any_order(name, z):
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


def tabulate(sar, elements, z):
    from pyscf.data import elements as pyscf_elements

    log_r0 = np.log(R_MIN)
    step = (np.log(R_MAX) - np.log(R_MIN)) / (N_RADIAL - 1)
    radii = np.exp(log_r0 + np.arange(N_RADIAL) * step)
    points = np.zeros((N_RADIAL, 3))
    points[:, 2] = radii

    symbol = pyscf_elements.ELEMENTS[z]
    mol, mf, dm, s = sar.solve(symbol, BASIS, FUNCTIONAL, True)
    dm_tot = dm[0] + dm[1]
    ao = mol.eval_gto("GTOval_sph_deriv1", points)
    rho = np.einsum("pi,ij,pj->p", ao[0], dm_tot, ao[0])
    drho_dz = 2.0 * np.einsum("pi,ij,pj->p", ao[3], dm_tot, ao[0])

    last = N_RADIAL - 1
    for k in range(N_RADIAL):
        if rho[k] > RHO_FLOOR:
            last = k
    r_cut = radii[last]

    log_rho = np.full(N_RADIAL, np.log(RHO_FLOOR))
    slope = np.zeros(N_RADIAL)
    log_rho[: last + 1] = np.log(rho[: last + 1])
    slope[: last + 1] = radii[: last + 1] * drho_dz[: last + 1] / rho[: last + 1]

    n_check = 4.0 * np.pi * np.trapezoid(radii**3 * rho, dx=step)
    return symbol, float(mf.e_tot), float(n_check), float(r_cut), log_rho, slope


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mqc-root", default=str(pathlib.Path.home() / "dev" / "metalquicha"),
                    help="path to a metalquicha checkout (default: ~/dev/metalquicha)")
    ap.add_argument("--out", default=str(DEFAULT_OUT),
                    help="output path (default: %s)" % DEFAULT_OUT)
    args = ap.parse_args()

    sar = load_solver(args.mqc_root)
    commit = mqc_commit(args.mqc_root)

    out = {
        "schema": {"name": "terco-hirshfeld-proatoms", "version": "1.0"},
        "provenance": {
            "recipe": "spherically averaged, spin-polarized free atom, Hund's "
                      "rule within each l (Lehtola configuration table); H-Kr",
            "basis": BASIS,
            "functional": FUNCTIONAL,
            "polarized": True,
            "method": "pyscf %s, mqc's own spherical_atom_reference.py solver "
                      "(see this script's module docstring)" % __import__("pyscf").__version__,
            "mqc_reference_commit": commit,
            "mqc_reference_files": [
                "backends/cenzontle/analysis/mqc_czt_hirshfeld.f90",
                "backends/cenzontle/scf/mqc_czt_spherical_atom.f90",
                "tools/spherical_atom/spherical_atom_reference.py",
            ],
            "generated_by": "scripts/gen_hirshfeld_proatoms.py",
        },
        "n_radial": N_RADIAL,
        "r_min": R_MIN,
        "r_max": R_MAX,
        "rho_floor": RHO_FLOOR,
        "elements": {},
    }

    for z in range(1, Z_MAX + 1):
        t0 = time.time()
        symbol, energy, n_check, r_cut, log_rho, slope = tabulate(sar, None, z)
        dt = time.time() - t0
        print("Z=%2d %-2s  E=%.10f  N_check=%.6f  r_cut=%.3f  (%.1fs)" %
              (z, symbol, energy, n_check, r_cut, dt), flush=True)
        if abs(n_check - z) > 1.0e-4:
            raise SystemExit("Z=%d: radial integral %.6f does not match the "
                             "nuclear charge; the table is wrong" % (z, n_check))
        out["elements"][str(z)] = {
            "symbol": symbol,
            "r_cut": r_cut,
            "log_rho": log_rho.tolist(),
            "slope": slope.tolist(),
        }

    dest = pathlib.Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w") as f:
        json.dump(out, f, separators=(",", ":"))
    print("wrote", dest, dest.stat().st_size, "bytes")


if __name__ == "__main__":
    main()
