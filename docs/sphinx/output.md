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
with `TERCO_NVTX`. The (T) correction reports its own timings. Over several ranks (`mpiexec -n N`)
the table is rank 0's, the total row ends `rank 0 of N`, the step table gains the three transfer rows
(`amplitudes to ranks`, `ladder wait (idle)`, `ladder to first rank`) and the energy box says
`ladder and (T) shared by N ranks` ({doc}`cli`).

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
SCF energy); a (T) that failed after a converged RI-CCSD adds `"ri-ccsd triples error"`
with the library's message, in place of the two (T) keys; an RI-CCSD that did not converge has `"ri-ccsd converged": false`
and no totals; the SUMMARY rows show the same quantities under shorter labels that fit the
block (`ri-ccsd (t) correction`, `ri-ccsd(t) total`, `ri-ccsd frozen core`, ...);
`--opt`/`--hessian` add `"optimization"`/`"vibrational"` blocks. Keys
whose default is decided from the system (the charge, the guess, ...)
appear in `"input"` only when the run was given one explicitly.

### `"build"` and `"timings"`

Every `--json` ends with two objects. `"build"` says what the binary is, so
a timing can be filed against the right build: `"device"` (a GPU build or
the host build), `"gpu_arch"` (the `-gpu=ccNN` flag; written for device
builds only, a host build has no such flag), `"lmax"`, `"df_lmax"` and
`"df_lmax_aux"` (`TERCO_LMAX`, `TERCO_DF_LMAX`, `TERCO_DF_LMAX_AUX`),
`"cart_order"` (`TERCO_CART_ORDER`), `"mpi"` and `"multicore"` (booleans),
`"compiler"`, `"compiler_version"`, `"build_type"` and `"commit"` (`git
describe --always --dirty` when CMake last configured, `"unknown"` outside a
git checkout; a later `git pull` that does not reconfigure leaves it stale).
`terco --version` prints the same facts on a `build:` line.

`"timings"` is in seconds (doubles), except the counts whose key starts with
`n_`, which are integers. Only the stages the run had are written. The
times are terco's own timers, the ones the verbose tables print. What stops
each clock: the Fock timer stops after the (synchronous) device kernels and
the copy of the result back have returned; the XC, DF-JK step and RI-MP2 step
tables wait for the device explicitly before they read the clock. The SCF
loop itself has no extra wait added for the timers.

What the numbers cover:

- **A "Fock build"** (`scf.fock_*`) is, per iteration: the density transform
  from the spherical to the Cartesian basis, the J/K build (with `--df`, the
  whole DF-JK build, which is why `dfjk.total` equals `scf.fock_total`), and
  the fold of the result back. It does **not** include the XC build, the PCM
  solve or the linear algebra (Fock assembly, DIIS, diagonalisation), which have
  their own rows (`dft.xc`, `pcm.*`) or none.
- **`scf.wall`** runs from the start of the SCF routine (after the input is
  validated) to the last iteration and the final transfers of the density and
  orbitals to the host. It includes `scf.setup`, the DFT grid, the PCM cavity,
  every Fock build, XC, linear algebra, and the second-order iterations; it
  excludes what the CLI does before it (reading the geometry and the basis, building the initial guess)
  and after it (the gradient, analysis, writing the file).
- **`scf.setup`** is the once-per-geometry part of that: shell pairs,
  one-electron integrals, the ERI context, the linear algebra set-up. The DFT
  grid (`dft.grid`) and the PCM cavity (`pcm.setup`) are timed separately and
  are *not* in `scf.setup`, though they are in `scf.wall`.
- **`--soscf`**: the Newton iterations are counted in `scf.n_iterations` (and
  the top-level `"iterations"`) but are not Fock builds of the main loop, so
  they are in neither `scf.n_fock_builds` nor `scf.fock_total`; only the DIIS
  iterations before the hand-over are. `scf.wall` covers both.
- **Several MPI ranks**: the file is written by rank 0 and holds rank 0's own
  clocks. `scf.*`, `dft.*`, `pcm.*`, `dfjk.*`, `rimp2.total` and
  `gradient.total` are therefore rank 0's wall times (they include any wait
  for the other ranks inside the timed region). `rimp2.step.*` is the maximum
  over the ranks of each step. `dfjk.step.*` is rank 0's own.

| Key | Written when | Meaning |
|---|---|---|
| `total_wall` | always | wall time from the start of `terco` to just before this file is written (the file write itself is not included) |
| `scf.wall` | an SCF ran here | the whole SCF, see above |
| `scf.setup` | same | once per geometry, see above |
| `scf.n_iterations` | same | SCF iterations, DIIS and second-order, equal to the top-level `"iterations"` |
| `scf.n_fock_builds` | same | Fock builds of the main (DIIS) loop, one per iteration there; the second-order solver's own are not counted |
| `scf.fock_first` | same | build 1 alone: it includes one-off device and kernel setup |
| `scf.fock_total` | same | all `scf.n_fock_builds` builds |
| `scf.fock_mean_after_first` | `scf.n_fock_builds` >= 2 | mean of builds 2..n, `(fock_total - fock_first)/(n - 1)`: the steady-state build time. Incremental builds are included, so the mean falls when the schedule skips more of the integrals; use `--no-incremental` for a like-for-like kernel time |
| `dft.grid` | a functional was given | the integration grid(s), VV10's included |
| `dft.xc` | same | semilocal XC over all iterations (VV10 excluded) |
| `dft.vv10` | the functional has VV10 | the nonlocal correlation over all iterations |
| `dfjk.total` | `--df` | the DF-JK Fock builds, equal to `scf.fock_total` (with `--df` the build is the DF-JK) |
| `dfjk.step.<name>` | `--df` | the DF-JK step table (metric, Cholesky, three-index integrals, ...) as printed, names lower-cased with other characters as `_` |
| `pcm.setup` | `--pcm` | cavity, K matrix and nuclear potential, once |
| `pcm.per_iteration` | `--pcm` | the charge solve and potential, mean over the iterations |
| `pcm.total` | `--pcm` | the same, summed |
| `rimp2.total` | `--rimp2` and the energy step ran | the RI-MP2 energy step, whole |
| `rimp2.step.<name>` | same | its step table (`metric_p_q`, `cholesky`, `mu_nu_p`, `transform_ia_p`, ...), the slowest rank's time over MPI. With `--grad` on one rank the pair energies come from the gradient's amplitudes and no energy step runs: the `rimp2` group is then absent and `gradient.total` carries the cost |
| `gradient.total` | `--grad` | `trc_run_gradient` plus the copy of the gradient to the host, so the device has finished |

`--opt`, `--hessian` and `--fukui` write only `total_wall`: they run SCFs at
other geometries or charges after the one the run is about, and the library
holds only the last. The groups are nested objects (`"scf"`, `"dft"`,
`"dfjk"`, `"pcm"`, `"rimp2"`, `"gradient"`); a later method adds its own
group (its stage table) without changing any key above. Example, a
`--rimp2` run:

```json
"timings": {
  "total_wall": 0.34,
  "scf": { "wall": 0.147, "setup": 0.011, "n_iterations": 10, "n_fock_builds": 10,
           "fock_first": 0.0116, "fock_total": 0.1335, "fock_mean_after_first": 0.0135 },
  "rimp2": { "total": 0.005, "step": { "metric_p_q": 0.0006, "cholesky": 0.0001, "mu_nu_p": 0.0038 } }
}
```

The library side is the C entry points `trc_timing_count` and `trc_timing`
(ABI 1.24), which hand out the same named values one at a time: loop `index`
from 0 to the count minus one with a name buffer of `TRC_TIMING_NAME_MAX`
characters. A buffer too short for a name is refused (`TRC_ERR_BADARG`, and
`trc_message` says the length needed), never truncated. The `rimp2.*` entries
belong to the latest `trc_run_rimp2` only.
