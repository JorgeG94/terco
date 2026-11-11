# Output

## The log

Each SCF iteration is one line. This is `terco water.xyz` (RHF/def2-SVP,
run against a CPU build):

```text
   it        E(total)            dE          RMS(D)      |FDS-SDF|   fock(s) linalg(s)
    1      -75.845437163849   -7.5845E+01    5.5990E-02    3.3350E+00     0.003     0.000
    2      -75.926042499268   -8.0605E-02    1.6019E-02    5.8492E-01     0.004     0.000 i
    ...
    9      -75.960983987068   -1.2523E-10    3.4864E-07    5.0742E-06     0.004     0.000
   10      -75.960983987078   -9.8197E-12    6.0952E-08    8.6926E-07     0.004     0.000
```

| Column | Meaning |
|---|---|
| `it` | SCF iteration number |
| `E(total)` | total energy, Eh, this iteration |
| `dE` | change in total energy from the previous iteration, Eh (`--econv`'s own quantity) |
| `RMS(D)` | RMS change in the density matrix |
| \|FDS-SDF\| | norm of the orthogonalised DIIS commutator (`--gconv`'s own quantity) |
| `fock(s)` | wall time this iteration spent in the Fock build |
| `linalg(s)` | wall time this iteration spent on the generalised eigenproblem/orthogonalisation |
| trailing `i` | this iteration's Fock build was incremental (from the density difference), the default; absent on the first iteration and whenever a full rebuild ran |

A DFT run (`-f <functional>`) adds an `xc(s)` column between `fock(s)` and
`linalg(s)` for the XC quadrature's own wall time.

A DFT run adds a `dft:` timing breakdown (grid generation, XC
evaluation, VV10) and, with `--df`, a `dfjk:` breakdown of the fitted J/K
build (metric Cholesky, the three-centre integrals, the half-transform,
...), each with its own achieved FLOP rate or bandwidth where that is
meaningful. An unrestricted run prints HOMO/LUMO per spin; `--rimp2` adds a
`RI-MP2 CORRELATION` and a `TOTAL RI-MP2 ENERGY` line and the opposite-/same-spin
components (the MP2 here is always the density-fitted one, and is labelled so);
`--riccsd`/`--riccsd-t` add the RI-MP2 and RI-CCSD correlation energies, the
(T) correction and the RI-CCSD and RI-CCSD(T) totals to the energy box (`TOTAL
ENERGY` stays the SCF energy, as with `--rimp2`; the correlated totals are
`TOTAL RI-CCSD ENERGY` and `TOTAL RI-CCSD(T) ENERGY`); `--pcm`
adds the dielectric solvation energy; `--grad` prints the Cartesian
gradient per atom (Eh/Bohr) with its RMS and max component.

The SCF converges when `dE` is below `--econv` **and** \|FDS-SDF\| is below
`--gconv`; the run then prints a boxed total/SCF energy and HOMO/LUMO, a
`SUMMARY` block (nuclear repulsion, one- and two-electron energy, and any
XC/solvation/correlation term that applies), and the wall time. The exit
status is 0, or 2 if the SCF did not converge within `--iterations`.

## `--json`

`--json <file>` writes a summary alongside the log. Every run's own
settings appear under `"input"`, grouped exactly as in the TOML file
({doc}`input`) -- including values nobody gave, filled in with their
defaults -- so the file written back as TOML reruns the same calculation.
`terco water.xyz --json water.json` (RHF/def2-SVP) writes:

```json
{
  "program": "terco",
  "version": "0.1.0",
  "input": {
    "molecule": { "geometry": "water.xyz" },
    "basis": {
      "main": "def2-svp",
      "auxiliary": "cc-pvdz-rifit",
      "jk_auxiliary": "def2-universal-jkfit",
      "cartesian": false
    },
    "scf": {
      "functional": "hf",
      "max_iterations": 100,
      "energy_tol": 0.1E-9,
      "commutator_tol": 0.10000000000000001E-4,
      "screening": 0.99999999999999998E-12,
      "incremental": true,
      "density_fitted": false,
      "convergence": { "diis_size": 10, "adiis": false, "soscf": false }
    },
    "dft": { "grid": 3, "nlc_grid": 1 },
    "correlation": { "method": "none", "frozen_core": -1 },
    "driver": { "gradient": false, "hessian": false, "optimize": false }
    /* every other [section] in docs/INPUT.md follows, each at its
       default, since none was given on this run */
  },
  "geometry": "water.xyz",
  "method": "RHF",
  "basis": ".../test/basis_sets/def2-svp.json",
  "guess": "sad",
  "charge": 0,
  "unpaired electrons": 0,
  "number of basis functions": 24,
  "spherical": true,
  "iterations": 10,
  "converged": true,
  "nuclear repulsion energy": 9.1895337626396856,
  "one-electron energy": -123.03292347071759,
  "two-electron energy": 37.8824057209996,
  "exchange-correlation energy": 0.0,
  "scf energy": -75.960983987078308,
  "total energy": -75.960983987078308,
  "homo energy": -0.49806581316084186,
  "lumo energy": 0.17614645009814475
}
```

A `--grad` run adds a `"gradient"` array; `--rimp2` adds an
`"mp2 opposite spin"`/`"mp2 same spin"` pair and a correlated total;
`--riccsd`/`--riccsd-t` add `"ri-ccsd converged"`, `"ri-ccsd iterations"`,
`"ri-ccsd frozen core orbitals"`, `"ri-ccsd active occupied"`, `"ri-ccsd virtual
orbitals"`, `"ri-ccsd auxiliary functions"`, `"ri-mp2 correlation
energy"`, `"ri-ccsd correlation energy"`, `"ri-ccsd t1 diagnostic"` and, when the
RI-CCSD converged, `"ri-ccsd total energy"` (and with (T) `"ri-ccsd triples
correction"` and `"ri-ccsd triples total energy"`; `"total energy"` stays the
SCF energy); an RI-CCSD that did not converge has `"ri-ccsd converged": false`
and no totals; the SUMMARY rows carry these same names;
`--opt`/`--hessian` add `"optimization"`/`"vibrational"` blocks. Keys
whose default is decided from the system (the charge, the guess, ...)
appear in `"input"` only when the run was given one explicitly.
