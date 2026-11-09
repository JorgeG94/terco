# terco input files

`terco` takes its settings either as command-line flags or from a TOML file:

```sh
terco water.xyz -b def2-svp -f wb97x-v --df --grad   # flags
terco water.toml                                     # the same run from a file
terco water.toml --grad                              # file, plus a flag on top
```

The file is the same run written down: every key below is one command-line
option, parsed and checked by the same code. Flags given after the file
override its keys. A key that is not in this list, or a value of the wrong
type, is an error that names the key; nothing is silently ignored.

## Example

```toml
[molecule]
geometry = "water.xyz"    # relative to this file
charge   = 0
unpaired = 0              # unpaired electrons, not the multiplicity

[basis]
main = "def2-svp"

[scf]
functional     = "wb97x-v"
density_fitted = true

[solvation]
model   = "cpcm"
epsilon = 78.3553

[driver]
gradient = true

[output]
json = "water.json"
```

## Rules

- **Sections and names.** Keys are `snake_case` and grouped by what they
  configure. Each maps to exactly one flag (tables below).
- **Defaults.** A key left out takes the flag's default, including the
  defaults that depend on the system (a `.CHRG`/`.UHF` file, the `sadq`
  guess for an ion, the `--opt-zero-modes` rule). An empty section is
  allowed and changes nothing.
- **Paths** (`geometry`, the `[basis]` keys when they are files, `output.json`,
  `output.save_scf`, `scf.restart`, `opt.neb.product`, `opt.neb.via`) are
  relative to the TOML file's directory, not the current directory, so a
  file and its inputs can move together. A basis *name* (`def2-svp`) is
  looked up exactly as on the command line.
- **Flags override the file.** `terco run.toml --grid 4` runs with grid 4
  whatever `[dft] grid` says. A switch set to `true` in the file cannot be
  turned off from the command line (there is no `--no-grad`); leave it out
  of the file instead. `incremental` is on by default and has both forms:
  `false` is `--no-incremental`, and the file's `false` with a command
  line `--incremental` (or the reverse) is refused as a contradiction.
- **The driver is explicit.** `[driver]` holds the switches that choose what
  runs after the SCF; there is no "section present means on".
- **The resolved input is recorded.** The JSON summary (`[output] json` or
  `--json`) has an `"input"` object in these sections, with every key as
  the run used it: from the file, from a flag that overrode it, or the
  flag's declared default. Keys whose default is decided from the system
  (the charge, the guess, ...) appear only when given. Written back as TOML
  it reruns the same calculation.
- `--version` and `--help` stay command-line only.

## Keys

### `[molecule]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `geometry` | string | *(positional)* | `.xyz` file; required |
| `charge` | integer | `--chrg` | molecular charge |
| `unpaired` | integer | `--uhf` | number of unpaired electrons |

### `[basis]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `main` | string | `--basis` | orbital basis: BSE JSON file or a basis name |
| `ecp` | string | `--ecp` | effective core potential: BSE JSON file or a basis name, overriding auto-attach from `main`'s own `ecp_potentials` |
| `auxiliary` | string | `--aux` | RI-MP2 auxiliary basis (cc-pvdz-rifit) |
| `jk_auxiliary` | string | `--jkfit` | J/K fitting basis for `density_fitted` (def2-universal-jkfit) |
| `cartesian` | bool | `--cart` | Cartesian AOs even for a spherical basis |
| `linear_dependence` | float | `--lindep` | drop eigenvectors of the unit-diagonal overlap below this (1e-8; 0 drops none) |

### `[scf]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `functional` | string | `--func` | `hf` (default) or a functional name |
| `guess` | string | `--guess` | `sad`, `sadq`, `sap`, `gwh`, `core` |
| `max_iterations` | integer | `--iterations` | maximum SCF iterations (100) |
| `energy_tol` | float | `--econv` | energy convergence, Eh (1e-10) |
| `commutator_tol` | float | `--gconv` | norm of the orthogonalised commutator \|FDS − SDF\| (1e-5) |
| `screening` | float | `--thresh` | integral screening threshold (1e-12) |
| `incremental` | bool | `--incremental` / `--no-incremental` | incremental Fock build (true, the default; false is a full build every iteration) |
| `density_fitted` | bool | `--df` | density-fitted J and K (basis: `[basis] jk_auxiliary`) |
| `restart` | string | `--restart` | read the SCF result from a `save_scf` file |

The SCF is converged when the energy change is below `energy_tol` **and** the
commutator norm is below `commutator_tol`. The energy error goes as the square
of the orbital error, a gradient's error linearly, so a tight optimisation
(`[opt] gradient_tol` 1e-5) wants `commutator_tol` near 1e-7. A Hessian
tightens both itself for its displaced SCFs.

### `[scf.convergence]`

How the SCF gets there; the defaults suit most systems.

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `diis_size` | integer | `--diis` | Fock/error pairs kept (10) |
| `adiis` | bool | `--adiis` | ADIIS before CDIIS, for SCFs stuck near 1e-2 |
| `level_shift` | float | `--level-shift` | virtual level shift during the damped start, Eh (1.0) |
| `damping` | float | `--damp` | weight of the new density while damping (0.3) |
| `damping_steps` | integer | `--damp-steps` | damped steps before DIIS (0 closed shell, 4 open shell) |
| `diis_start` | float | `--diis-start` | commutator norm below which DIIS may start |
| `soscf` | bool | `--soscf` | second-order SCF after DIIS (RHF only) |
| `soscf_start` | float | `--soscf-start` | commutator norm at which SOSCF takes over (1e-2) |

### `[dft]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `grid` | integer | `--grid` | integration grid level (3) |
| `nlc_grid` | integer | `--nlc-grid` | VV10 grid level (1) |

### `[solvation]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `model` | string | `--pcm` | `cpcm` or `iefpcm` |
| `epsilon` | float | `--eps` | dielectric constant (78.3553 for water) |

### `[correlation]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `method` | string | `--rimp2` | `rimp2`; absent means none |
| `frozen_core` | integer | `--frozen` | frozen core orbitals; −1 counts by element |
| `aux_block` | integer | `--aux-block` | auxiliary functions per three-index block |

### `[driver]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `gradient` | bool | `--grad` | analytic nuclear gradient |
| `hessian` | bool | `--hessian` | FD Hessian, frequencies, thermochemistry |
| `optimize` | bool | `--opt` | geometry optimisation (needs DL-FIND) |

`optimize` and `hessian` together are `--opt-freq`: optimise, then the Hessian
at the optimised geometry. With RI-MP2 (`[correlation] method`) or a double
hybrid, every optimisation step and every Hessian displacement reruns the
RI-MP2, so both act on the correlated total energy (a double hybrid's
E_KS + c_os E_os + c_ss E_ss), which is also the energy the thermochemistry
starts from. What the correlated gradient refuses (an unrestricted
reference, PCM) is refused: by `optimize` before the first SCF, by a
`hessian` without `optimize` after the reference SCF, before any
displacement runs.

### `[analysis]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `charges` | string | `--charges` | atomic charges: `hirshfeld`, `chelpg`, or both as `"hirshfeld,chelpg"` |
| `fukui` | string | `--fukui` | condensed Fukui indices from `hirshfeld` or `chelpg` charges |

### `[hessian]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `step` | float | `--step` | displacement, Bohr (0.005) |

### `[thermo]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `temperature` | float | `--temperature` | K (298.15) |
| `pressure` | float | `--pressure` | atm (1.0) |
| `symmetry_number` | integer | `--symmetry-number` | rotational symmetry number (1) |
| `masses` | string | `--masses` | `isotope` or `average` |

### `[opt]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `coordinates` | string | `--opt-coords` | `cartesian`, `hdlc`, `hdlc-tc`, `dlc`, `dlc-tc` |
| `max_steps` | integer | `--opt-maxiter` | maximum steps (100) |
| `gradient_tol` | float | `--opt-gtol` | max gradient component, Eh/Bohr (4.5e-4) |
| `energy_tol` | float | `--opt-etol` | energy change, Eh |
| `scf_guess` | string | `--opt-guess` | warm start between steps: `density` or `projected` |
| `algorithm` | string | `--opt-algorithm` | `sd`, `cg`, `cg-auto`, `lbfgs`, `nr`, `damped`, `stochastic`, `genetic` (P-RFO is a saddle search: `saddle_method`) |
| `target` | string | `--opt-target` | `minimum` or `saddle` |
| `saddle_method` | string | `--opt-saddle-method` | `prfo` or `dimer` |
| `connect` | bool | `--opt-connect` | minimise both ways from the saddle |
| `initial_hessian` | string | `--opt-init-hessian` | `terco` or `dlfind` |
| `hessian_update` | string | `--opt-hessian-update` | `none`, `powell`, `bofill`, `bfgs` |
| `max_hessian_updates` | integer | `--opt-max-updates` | updates before a recompute |
| `zero_modes` | integer | `--opt-zero-modes` | zero-frequency modes |
| `soft_modes` | float | `--opt-soft-modes` | soft-mode threshold, Eh/Bohr² |
| `max_step` | float | `--opt-max-step` | trust-radius ceiling, Bohr |
| `line_search` | string | `--opt-linesearch` | `auto`, `scaling`, `energy-trust`, `gradient-trust` |
| `step_scale` | float | `--opt-scalestep` | step scale for `line_search = "scaling"` |
| `lbfgs_memory` | integer | `--opt-lbfgs-memory` | L-BFGS history length |
| `print_level` | integer | `--opt-print-level` | DL-FIND verbosity: 0, 2, 4, 6 |
| `frozen_atoms` | integer array | `--opt-frozen` | 1-based atoms to freeze, e.g. `[1, 2, 3, 7]` |
| `constraints` | string array | `--opt-constrain` | e.g. `["angle 2 1 3", "bond 1 2"]` |
| `seed` | integer | `--opt-seed` | seed for stochastic/genetic/dimer |

### `[opt.neb]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `product` | string | `--opt-neb` | product geometry (`.xyz`); turns NEB on |
| `images` | integer | `--opt-neb-images` | number of images |
| `spring` | float | `--opt-neb-spring` | spring constant, Eh/Bohr² |
| `climb` | bool | `--opt-neb-climb` | climbing image |
| `endpoints` | string | `--opt-neb-endpoints` | `free`, `perpendicular`, `frozen` |
| `via` | string | `--opt-neb-via` | intermediate geometries (multi-frame `.xyz`) |

### `[opt.global]`

| Key | Type | Flag |
|---|---|---|
| `population` | integer | `--opt-global-pop-size` |
| `initial_population` | integer | `--opt-global-init-pop-size` |
| `max_cycles` | integer | `--opt-global-maxcycle` |
| `keep_best` | integer | `--opt-global-nsave` |
| `distribution` | integer | `--opt-global-distribution` |
| `reset_after` | integer | `--opt-global-reset` |
| `radius` | float | `--opt-global-radius` |
| `contraction` | float | `--opt-global-contraction` |
| `radius_tol` | float | `--opt-global-tolerance-r` |
| `gradient_tol` | float | `--opt-global-tolerance-g` |
| `mutation_rate` | float | `--opt-global-mutation-rate` |
| `death_rate` | float | `--opt-global-death-rate` |
| `scale_factor` | float | `--opt-global-scalefac` |

Meanings and defaults as in `terco --help`.

### `[output]`

| Key | Type | Flag | Meaning |
|---|---|---|---|
| `json` | string | `--json` | JSON summary file |
| `save_scf` | string | `--save-scf` | write the SCF result for a later restart |
| `silent` | bool | `--silent` | summary only |
| `log_level` | string | `--logger` | `knowledge`, `error`, `warning`, `performance`, `info`, `large_info`, `verbose`, `debug` |
