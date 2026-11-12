# Command-line program

The CMake build also makes a `terco` program (`-DTERCO_BUILD_APP=OFF`
skips it). It reads an xyz file and runs an SCF, in the spirit of xtb; the
same run can also be written as a TOML file (see {doc}`input`).

```sh
terco --help               # every flag, grouped as below
terco --help opt           # the --opt-* geometry-optimisation flags
```

A basis name (`-b def2-svp`) is looked up as `<name>.json` in each
directory of `$TERCO_BASIS_PATH`, then in `share/terco/basis_sets` under
the install prefix. The exit status is 0 on success, 1 on an error or a
refused request, and 2 if the SCF, a `--riccsd` iteration or an `--opt` did not
converge.

Every flag below is `terco --help`'s own text, copied verbatim (the code
names are set in code font); the
examples were run against a CPU (gfortran) build with `test/`'s own
molecules and basis sets, and each one reproduced the output shown.

## Molecule and basis

| Flag | Meaning |
|---|---|
| `geometry` (positional) | molecule in xyz format, Angstrom (or a TOML input file) |
| `-b, --basis` | basis set: a BSE JSON file, or a name in `$TERCO_BASIS_PATH` then the installed basis sets (`def2-svp`) |
| `--ecp` | effective core potential, looked up like `--basis` |
| `-c, --chrg` | molecular charge (default: a `.CHRG` file here, else 0) |
| `-u, --uhf` | number of unpaired electrons (default: a `.UHF` file here, else 0 or 1 by electron count) |
| `--cart` | Cartesian AOs even when the basis is defined spherical |
| `--lindep` | linear-dependence threshold on the normalised overlap's eigenvalues (1e-8) |

```sh
terco water.xyz                              # RHF/def2-SVP from a SAD guess
```

## DFT

| Flag | Meaning |
|---|---|
| `-f, --func` | `hf` (default) or a functional name -- see `terco --help` for the full list, or {doc}`functionals` |
| `--grid` | DFT integration grid level (3) |
| `--nlc-grid` | grid level for VV10 nonlocal correlation, e.g. `wb97x-v` (1) |

```sh
terco water.xyz -b cc-pvdz -f pbe0           # RKS/PBE0
terco oh.xyz -u 1 -f b3lyp                   # UHF/B3LYP doublet
```

## SCF convergence

| Flag | Meaning |
|---|---|
| `--iterations` | maximum SCF iterations (100) |
| `--econv` | energy convergence, Eh (1e-10) |
| `--gconv` | convergence on the DIIS commutator norm \|FDS-SDF\| (1e-5) |
| `--diis` | DIIS subspace size (10) |
| `--adiis` | ADIIS before CDIIS, for an SCF stuck near 1e-2 |
| `--level-shift`, `--damp`, `--damp-steps`, `--diis-start` | the damped-start schedule before DIIS takes over |
| `--soscf`, `--soscf-start` | second-order SCF hand-off from DIIS (RHF only) |
| `--thresh` | integral screening threshold (1e-12) |
| `--incremental` / `--no-incremental` | incremental Fock build (on by default; refused with `--df`) |

## Density fitting

| Flag | Meaning |
|---|---|
| `--df` | density-fitted J and K (restricted, Cartesian only) |
| `--jkfit` | fitting basis for `--df` (`def2-universal-jkfit`) |

```sh
terco water.xyz -b def2-svp -f pbe0 --df --grad --cart   # DF gradient needs a Cartesian basis
```

## Gradients, Hessian, geometry optimisation

| Flag | Meaning |
|---|---|
| `--grad` | analytic nuclear gradient after the SCF |
| `--hessian` | finite-difference Cartesian Hessian, then harmonic frequencies and RRHO thermochemistry |
| `--step` | `--hessian`'s displacement, Bohr (0.005) |
| `--temperature`, `--pressure`, `--symmetry-number`, `--masses` | `--hessian`'s thermochemistry inputs |
| `--opt` | optimise the geometry with DL-FIND (needs a `-DTERCO_ENABLE_DLFIND=ON` build; see {doc}`install`) |
| `--opt-freq` | `--opt` followed by `--hessian` at the optimised geometry |

`--opt` and its `--opt-*` sub-options (coordinates, algorithm, constraints,
NEB, basin-hopping/genetic global search, ...) are listed by
`terco --help opt`; they are unchanged by this release and not duplicated
here.

```sh
terco water.xyz -b def2-svp -f wb97x-v --grad
```

## Fukui / charge analysis

| Flag | Meaning |
|---|---|
| `--charges` | atomic charge analysis: `hirshfeld`, `chelpg`, or both as `hirshfeld,chelpg` |
| `--fukui` | condensed Fukui indices (from `hirshfeld` or `chelpg` charges; needs a closed-shell neutral) |

```sh
terco water.xyz -b sto-3g --charges hirshfeld,chelpg
```

## PCM solvation

| Flag | Meaning |
|---|---|
| `--pcm` | `cpcm` or `iefpcm` |
| `--eps` | solvent dielectric constant, e.g. `78.3553` for water |

```sh
terco water.xyz -f b3lyp --pcm cpcm --eps 78.3553
```

## RI-MP2

| Flag | Meaning |
|---|---|
| `--rimp2` | RI-MP2 correlation energy after the SCF (restricted or unrestricted; unrestricted refuses more than one rank) |
| `--aux` | auxiliary basis for `--rimp2`, `--riccsd` and `--riccsd-t` (`cc-pvdz-rifit` whatever the orbital basis: give the RI fitting set that matches it) |
| `--frozen` | frozen core orbitals for `--rimp2`, `--riccsd` and `--riccsd-t`; -1, the default, counts them by element (1 per Li-Ne atom, 5 per Na-Ar, ...); 0 correlates all electrons |
| `--aux-block` | auxiliary functions per three-index block for `--rimp2`, `--riccsd` and `--riccsd-t` (default: as many as keep the block under 2 GB) |

```sh
terco water.xyz -b 6-31g --rimp2             # RHF + RI-MP2 (cc-pvdz-rifit)
```

## RI-CCSD and RI-CCSD(T)

The method is the density-fitted one throughout: every name of it, on the
command line, in the input file and in the output, says RI-CCSD.

| Flag | Meaning |
|---|---|
| `--riccsd` | closed-shell RI-CCSD correlation energy after the SCF (RHF, no ECP; over several ranks the (T) is shared; also reports the RI-MP2 correlation; cannot be combined with `--rimp2`, a functional, `--uhf` or `--grad`/`--hessian`/`--opt`) |
| `--riccsd-t` | `--riccsd` followed by the perturbative triples (T): RI-CCSD(T) |
| `--cc-econv` | `--riccsd` energy convergence: \|dE\| between iterations, Eh (1e-8) |
| `--cc-tconv` | `--riccsd` amplitude convergence: the largest amplitude change (1e-6) |
| `--cc-iterations` | maximum `--riccsd` iterations (100); a run that does not converge is reported as such and exits with status 2 |
| `--cc-diis` | `--riccsd` DIIS history length (8, at most 50); 1 switches DIIS off |
| `--cc-ladder` | how `--riccsd` forms the particle ladder: auto (by operation count), assembled or direct; the energy does not depend on it (auto) |
| `--cc-triples` | which routine computes the (T) of `--riccsd-t`: auto (the ijk device routine in a GPU build, the host reference one otherwise; the ijk routine whenever there is more than one rank), reference (host, from the host copies; cross-checks a GPU build) tiled (the device routine over virtual triples, also on the host, slower there) or ijk (the device routine over occupied triples, the default in a GPU build); the energy does not depend on it (auto) |
| `--cc-triples-budget` | memory for the working tile of the tiled (T) routine of `--riccsd-t`, GiB (1; any positive value, the repacked amplitudes come on top); a budget below one triple (48 o^3 bytes, o active occupied orbitals) makes the (T) step fail after the RI-CCSD has converged, with exit status 1; no effect on the reference routine |

`--aux`, `--frozen` and `--aux-block` (above) apply to `--riccsd` and
`--riccsd-t` as to `--rimp2`. The auxiliary basis defaults to `cc-pvdz-rifit`
whatever the orbital basis, so give the RI fitting set that matches it
(`-b cc-pvtz --aux cc-pvtz-rifit`).

```sh
terco hirshfeld_water.xyz -b cc-pvdz --riccsd-t                 # RHF + RI-CCSD(T), cc-pvdz-rifit, core counted by element
terco hirshfeld_water.xyz -b cc-pvdz --riccsd --frozen 0        # all electrons correlated, RI-CCSD only
terco hirshfeld_water.xyz -b cc-pvdz --riccsd-t --cc-econv 1e-10 --cc-tconv 1e-8
```

The same from an input file (see `docs/INPUT.md`): `[correlation]
method = "ri-ccsd(t)"`, and the solver settings in a `[riccsd]` table. The
plain names `ccsd` and `ccsd(t)` are refused.

Both criteria must hold for the iteration to stop: the change of the
correlation energy is below `--cc-econv` **and** the largest change of any
amplitude is below `--cc-tconv`. An RI-CCSD that does not converge in
`--cc-iterations` is not an error you can miss: the energy box prints
`RI-CCSD DID NOT CONVERGE` above the last iterate, the SUMMARY block says
`NOT CONVERGED`, `--json` has `"ri-ccsd converged": false` and carries no total
energies, (T) is not run, the library's message is printed, and the exit status
is 2. If the RI-CCSD converged and the (T) step then failed (for example for
device memory), the summary and the JSON keep the RI-CCSD energies, say that (T)
failed and why, and the exit status is 1.

`--riccsd` is refused, before the SCF, together with `--rimp2` (the MP2
correlation energy is part of its own report), a functional (`-f`), unpaired
electrons (`--uhf`, or an odd electron count from `--chrg`), `--grad`,
`--hessian`, `--opt`, a molecule with an effective core potential and a
`--frozen` that leaves no correlated orbital; the `--cc-*` options need
`--riccsd` or `--riccsd-t`, and `--cc-triples` and `--cc-triples-budget` need `--riccsd-t`. Over
several ranks `--cc-triples reference` and `tiled` are refused (see below). A basis with linearly dependent functions
dropped by `--lindep` (`[basis] linear_dependence`) is refused after the SCF,
since the method needs every orbital of the basis: lower the threshold (0 drops
none) or remove the near-duplicate functions. Over a `--pcm` reference the
correlation is computed on the orbitals of the solvated SCF and has no solvent
response of its own (as `--rimp2`); terco warns.

### Several GPUs

`mpiexec -n 4 terco water.xyz -b cc-pvdz --riccsd-t` runs on four ranks, one
rank per GPU (terco binds each rank to its own device; do not put two ranks on
one GPU). What is distributed today is the (T) correction only. Every rank
builds the density-fitted integrals (replicated, so each GPU needs the memory of
a one-GPU run), the RI-CCSD iterations run on the first rank while the others
wait, the converged amplitudes are broadcast, and the (T) work, the occupied
pairs, is dealt over the ranks. The (T) energy is the same on every rank and
equals the one-rank value to rounding. The RI-CCSD energy is computed by the
first rank alone, so given the same orbitals it is the one-rank number; the
multi-rank SCF that precedes it in the CLI differs from a one-rank SCF at about
1e-14 Eh, and so does everything after it. Output (tables, summary, JSON) comes
from the first rank, and with `--logger verbose` the (T) step lists each rank's
triples and seconds.

Plain `--riccsd` over several ranks gains nothing today: the iterations run on
the first rank and the others only matter for (T). Distributing the iterations
is the next step.

Only the ijk (T) routine shares the work. `--cc-triples auto` picks it whenever
there is more than one rank, in a host build too; an explicit `reference` or
`tiled` is refused over several ranks, before any work. A failure on one rank (an
allocation, a refused setting) stops all ranks with the reason of the
lowest-numbered rank that met it; it does not leave the others waiting.

### What it is, and what it is not yet

Closed-shell RI-CCSD(T): the equations are those of PySCF's `rccsd` and
`ccsd_t` with density-fitted integrals, and the energies are validated against
PySCF's DF-CCSD(T) (water and hydrogen peroxide, cc-pVDZ with
cc-pVDZ-RIFIT: CCSD(T) total energies within 5e-12 Eh when run with
`--econv 1e-13 --gconv 1e-11 --cc-econv 1e-12 --cc-tconv 1e-10`; at the default
thresholds, 1e-8 and 1e-6, the agreement is about 1e-8 Eh).

Today it covers: a restricted Hartree-Fock reference with canonical orbitals,
one GPU per rank (the (T) is shared over ranks, the iterations are not yet), no effective core potential, **no gradients** (so no
optimisation or Hessian), no restart of the RI-CCSD itself. The dense
`(ia|bc)` block is what limits the system size for now.

The shipped kernels stop at d functions in the orbital basis, so cc-pVTZ needs
`-DTERCO_LMAX=3` (the SCF itself needs it for f functions, and `TERCO_DF_LMAX`
follows it by default). The auxiliary ceiling `TERCO_DF_LMAX_AUX` defaults to
4, so cc-pVDZ with cc-pVDZ-RIFIT works in a default build. Otherwise the run is
refused with, for example, `the orbital basis has l = 3, above this build's
density-fitting ceiling 2 (TERCO_DF_LMAX)`.

Memory, with `o` the active (not frozen) occupied orbitals, `v` the virtual
orbitals and `naux` the auxiliary functions: the device holds about
`14.8 o^2 v^2 + o v^3 + (o^2 + o v + v^2) naux + 2 diis (o v + o^2 v^2)`
doubles (`diis` is `--cc-diis`) plus about 1.2 GB of ladder scratch; (T) adds
`2 o^2 v^2 + o^3 v` doubles and its tile (`--cc-triples-budget`, 1 GiB by
default, tiled routine). For 25 active occupied, 130 virtual and 630 auxiliary functions that
is about 4.4 GB for the RI-CCSD and about 1.3 GB more for (T) at the default
tile. In a GPU build the host keeps a copy of the same arrays; in a host (CPU)
build this is all host memory.

One measured data point (adenine, cc-pVDZ with cc-pVDZ-RIFIT, frozen core: 25
active occupied, 130 virtual, 630 auxiliary functions; one V100, untuned): 0.43
s per RI-CCSD iteration, 15 iterations, 13.3 s for (T).

The (T) step has two routines. The tiled one is the device routine; it also
runs on the host, more slowly there, and `--cc-triples-budget` sizes its tile.
The reference one runs on the host from the host copies and has no budget; a
GPU build can be forced to it with `--cc-triples reference` to cross-check.
By default (`--cc-triples auto`) a host (CPU) build uses the reference routine.

The coupled-cluster work is mostly matrix products, so link a
threaded BLAS and set `OMP_NUM_THREADS`. With NVHPC, `-DTERCO_ENABLE_GPU=OFF
-DTERCO_ENABLE_MULTICORE=ON` builds the host version with `do concurrent`
running on CPU threads (`OMP_NUM_THREADS` or `ACC_NUM_CORES` sets the count).

The `TOTAL ENERGY` line of the energy box and the `"total energy"` key of the
JSON stay the SCF energy, as with `--rimp2`. The correlated totals are the
`TOTAL RI-CCSD ENERGY` and `TOTAL RI-CCSD(T) ENERGY` lines, and
`"ri-ccsd total energy"` and `"ri-ccsd triples total energy"` in the JSON.

## Output and restart

| Flag | Meaning |
|---|---|
| `--json` | also write a JSON summary (see {doc}`output`) |
| `--save-scf` | write the SCF result for a later `--restart` |
| `--restart` | read the SCF result from a `--save-scf` file (molecule, basis and method must match) |
| `-s, --silent` | only the summary: no SCF iterations, timings or guess details |
| `--logger` | log verbosity: `knowledge`, `error`, `warning`, `performance`, `info`, `large_info`, `verbose`, `debug` (`info`) |
| `--version` | print the version and stop |
| `-h, --help` | show the help |

```sh
terco water.xyz --json water.json            # also write a JSON summary
```
