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
the install prefix. The exit status is 0 on success, 1 on an error and 2
if the SCF did not converge.

Every flag below is `terco --help`'s own text, copied verbatim; the
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
| `--mom`, `--level-shift-until` | maximum overlap occupation; keep the level shift on through DIIS |
| `--gdm`, `--gdm-start`, `--gdm-max-rotation` | geometric direct minimisation (L-BFGS over orbital rotations) in place of DIIS near convergence, for stalled SCFs; refused with `--soscf` |
| `--thresh` | integral screening threshold (1e-12) |
| `--incremental` / `--no-incremental` | incremental Fock build (on by default; refused with `--df`) |

## Point-group symmetry

Off by default. See {doc}`theory` for the method and {doc}`output` for what a
symmetric run prints.

| Flag | Meaning |
|---|---|
| `--symmetry` | `off` (the default); `on` or `auto`: the largest Abelian group of the molecule (D2h or a subgroup; tops, linear molecules and atoms use their Abelian subgroup), results reported in your input frame, and a path the symmetric SCF does not support runs in C1 with one log line saying why; or a group name (`C2v`, `D2h`, ...), which must be the detected group or a subgroup of it, reports the standard frame of the group and refuses those paths |
| `--symmetry-tol` | geometry tolerance for the detection, Bohr (1e-5); an error with `--symmetry off` |
| `--irrep-nelec` | electrons per irrep, e.g. `A'=6,A"=4` (restricted: even counts) or `A'=4/3,A"=3/3` (alpha/beta); unlisted irreps hold none. Makes an `on`/`auto` request strict. A gradient is refused for occupations that are not the aufbau ones. Default: aufbau |
| `--sym-skeleton` / `--no-sym-skeleton` | force the skeleton build (Fock matrix and four-centre gradient from the orbit representatives of the integral quartets) or the full build; the result is the same to round-off. Without either, a size rule chooses: the skeleton build when the molecule has at least two atoms per symmetry-inequivalent atom. Needs a group other than C1; density fitting and several ranks take the full build |

`--opt`, `--hessian`, `--rimp2`, `--df`, `--pcm`, `--soscf`, `--mom`,
`--gdm`, `--fukui`, `--restart` and more than one MPI rank are not supported
with symmetry: `--symmetry on` runs them in C1, a named group refuses them.

```sh
terco water.xyz -b 6-31gs --symmetry on --grad       # C2v, input frame
terco water.xyz -b 6-31gs --symmetry C2v             # standard frame of C2v
```

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
| `--rimp2` | RI-MP2 correlation energy after the SCF |
| `--aux` | auxiliary basis for `--rimp2` (`cc-pvdz-rifit`) |
| `--frozen` | frozen core orbitals (-1: count by element) |
| `--aux-block` | auxiliary functions per three-index block |

```sh
terco water.xyz -b 6-31g --rimp2             # RHF + RI-MP2 (cc-pvdz-rifit)
```

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
