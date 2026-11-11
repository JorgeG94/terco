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

## Step tables

At the `performance` log level (the default `--logger info` prints them) the
correlated methods end with a table of timed steps, one row per step: the
device time to completion, the floating-point operations the step executed
(counted from the dimensions of its matrix products, `2 m n k` each, batched
calls times their count) and the rate they give, each count with the SI prefix
that keeps it under 1000. RI-MP2 prints `rimp2:` rows (metric, Cholesky,
three-index integrals, transforms, the K matrix products, the pair sums) and a
total with the work over the whole wall time.

`--riccsd` and `--riccsd-t` print the same table with `rccsd:` rows after the
iteration table:

```
  rccsd: step                           time (s)           work            rate
  rccsd:   integrals                      0.004              -               -
  rccsd:     cholesky                     0.000   197.57 kFLOP    2.49 GFLOP/s
  ...
  rccsd:   tau, intermediates             0.003     3.08 MFLOP  902.91 MFLOP/s
  rccsd:   W intermediates                0.012    31.01 MFLOP    2.64 GFLOP/s
  rccsd:   ladder, assembled              0.026   286.71 MFLOP   10.83 GFLOP/s
  ...
  rccsd:   total                          0.079   396.93 MFLOP    5.06 GFLOP/s   rank 0 of 1
  rccsd: mean time of an iteration      0.007 s over 11 iterations
```

The `integrals` row is the whole build of the fitted blocks; the indented rows
under it (metric, Cholesky, three-index integrals, half transform, `B = X L^-T`,
the four blocks) are pieces it times itself, so they do not add up to it (the
orbital transforms and the copy of the blocks back to the host are in the first
only). The other rows are the stages of one iteration summed over all of them:
tau and the one-particle intermediates, the W intermediates, the particle
ladder (named for the form that ran, assembled or direct) and its dressing, the
T1 equation, the T2 equation, its three ring terms, DIIS and the energy. A
stage's work counts its matrix products only; the packs, scatters, denominators
and per-element kernels (and the spherical-basis conversions of the integral
rows) add time but no count, so a row that mixes them (the intermediates, W,
the T1 and T2 equations) shows a rate below the card's matrix-product rate by
the share of its time those take. The ladder and the ring terms are nearly all
matrix products. Each step is also an NVTX range (`rccsd: <name>`) in a build
with `TERCO_NVTX`. The (T) correction reports its own timings.

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
and no totals; the SUMMARY rows show the same quantities under shorter labels that fit the
block (`ri-ccsd (t) correction`, `ri-ccsd(t) total`, `ri-ccsd frozen core`, ...);
`--opt`/`--hessian` add `"optimization"`/`"vibrational"` blocks. Keys
whose default is decided from the system (the charge, the guess, ...)
appear in `"input"` only when the run was given one explicitly.
