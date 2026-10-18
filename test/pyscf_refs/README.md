# pyscf references

Jorge's rule, from the `scf_rks_vs_pyscf`/UKS build that went red because
pyscf's own UKS did not reach `conv_tol=1e-12` in 200 cycles (terco agreed
with it to 3.1e-11 at terco's density, 4.1e-9 SCF to SCF -- the regression
was pyscf's convergence, not terco's): **pyscf never runs in CI.** It runs
here, once, offline, on a developer machine with pyscf installed, and the
numbers it produces are committed as plain text in this directory (plus
`test/reference_xc.dat`, in the repository's existing recorded-oracle
convention). Every `check_*_vs_pyscf` and `scf_*` program under `test/`
reads its own file back and compares in Fortran -- none of them import
pyscf, or numpy, ever.

## Regenerating

```
/home/jorge/dev/metalquicha/.venv/bin/python3.12 test/pyscf_refs/generate.py            # every case
/home/jorge/dev/metalquicha/.venv/bin/python3.12 test/pyscf_refs/generate.py rks vv10   # named cases only
```

That venv has **pyscf 2.14.0** (numpy 2.5.2) -- the version every committed
file's `# pyscf: ...` header line records. Any pyscf recent enough to carry
the same functionals and libxc kernels should reproduce the same numbers
to the gates below; the version is recorded for provenance, not pinned by
the generator itself.

`generate.py` runs each case's `test/check_*`/`scf_*` binary (built, CPU is
fine -- nothing here touches the GPU) to get that program's fixed,
deterministic geometry/basis/density or grid, never terco's own converged
SCF output (see "Fixed inputs" below), then evaluates pyscf on exactly that
input. Build a CPU config first:

```
FC=gfortran cmake -S . -B build-cpu -DCMAKE_BUILD_TYPE=Release -DTERCO_ENABLE_GPU=OFF -DTERCO_ENABLE_TESTING=ON
cmake --build build-cpu -j4
TERCO_PYSCF_REF_BUILD=$PWD/build-cpu /path/to/venv/python3.12 test/pyscf_refs/generate.py
```

`TERCO_PYSCF_REF_BUILD` defaults to `<repo>/build-cpu`. The f/g cases
(`rys_fock_water_ccpvtz`, `rys_fock_formamide_def2tzvp`, `df_vs_pyscf`) need
a build at `-DTERCO_LMAX=3 -DTERCO_DF_LMAX=3` -- point
`TERCO_PYSCF_REF_BUILD` at that build for those cases; run default-ceiling
cases against a default build. A case whose own build-side program
self-skips (the basis needs more angular momentum than this build's
`TRC_LMAX`/`DF_LMAX` reaches) prints a note and leaves any existing
committed reference alone, rather than erroring -- run `generate.py` again
against the right build to fill it in.

`generate.py` always forces `TERCO_CART_ORDER=lexicographic` on the binary
it runs, because every within-shell-ordering check here
(`check_mult`/`check_ao`/`check_xc_energy`'s pyscf halves) is committed for
that convention only -- the CI default, and CI never overrides it. A
`TERCO_CART_ORDER=gamess` build can still run every OTHER pyscf comparison
in this tree (they are shell-order and AO-normalisation checks, not
within-shell-order ones), but those three will show a shape/number
mismatch against a gamess-order build; that is a known gap, not a bug (see
"Left open" at the bottom).

## Local-only references

`df_vs_pyscf.ref` is the one reference NOT committed -- it is the (ff|g)
3-centre tensor on water/cc-pVTZ x cc-pVTZ-RI, ~18 MB as `%.17e` text, and
`.gitignore` excludes it (`test/pyscf_refs/df_vs_pyscf.ref`). Without it,
`check_df_vs_pyscf` prints a `SKIP` line (ctest's `SKIP_REGULAR_EXPRESSION`
reports it as Skipped, not Passed or Failed) in every configuration,
including the `-DTERCO_DF_LMAX=3` CI job -- the one deliberate exception to
"a missing reference fails loudly past the LMAX/DF_LMAX self-skip" that the
rest of this directory holds to. `check_staleness.py` treats it the same
way: absent is reported and skipped, not an error, but if it IS present its
hash is still checked like every other reference.

To regenerate it locally and run the comparison for real:

```
# 1. A DF_LMAX=3 build -- this case needs the g-shell auxiliary classes,
#    which the shipped ceiling (DF_LMAX=2) does not have:
FC=gfortran cmake -S . -B build-lmax3 -DCMAKE_BUILD_TYPE=Release \
    -DTERCO_ENABLE_GPU=OFF -DTERCO_ENABLE_TESTING=ON \
    -DTERCO_LMAX=3 -DTERCO_DF_LMAX=3
cmake --build build-lmax3 -j4

# 2. generate.py's positional arguments select cases by name -- df_vs_pyscf
#    alone, not every case:
TERCO_PYSCF_REF_BUILD=$PWD/build-lmax3 \
    /home/jorge/dev/metalquicha/.venv/bin/python3.12 test/pyscf_refs/generate.py df_vs_pyscf

# 3. Compare for real:
ctest --test-dir build-lmax3 -R check_df_vs_pyscf --output-on-failure
```

Measured agreement (gfortran 15, CPU): 2.5e-13 against the 1e-10 gate, same
as every run before this file stopped being committed.

## The convergence-refusal rule

Every SCF comparison (`rks`, `scf_robust`, `rimp2`) converges pyscf
tightly -- `conv_tol=1e-12`, extra cycles, and a `.newton()` (second-order)
retry when the first-order solver does not reach it -- and **refuses to
write a case whose SCF did not converge**: `generate.py` raises
`ConvergenceError` (a `SystemExit`, with the fit measured) rather than
commit a number next to an unconverged energy. If a case starts refusing,
the fix is in `generate.py`'s SCF setup for that case (tighter
`conv_tol`/`max_cycle`, a level shift, `.newton()`), never a looser gate on
the Fortran side.

## Fixed inputs, never terco's output

Every fixed-density/fixed-grid check (`rys_fock`, `erf_fock`, `df_vs_pyscf`,
`rys_grad`, `vv10`, `check_xc_energy`, `check_ao`) builds its density (and,
where relevant, its grid) from a seeded LCG or a closed-form recipe baked
into the Fortran `check_*`/`scf_*` program itself -- never from an SCF
iteration, and never read back from a previous run. Because the recipe is
deterministic and lives in the Fortran source, re-running the SAME binary
always reproduces the SAME density/grid bit for bit; `generate.py` need not
(and does not) reimplement that recipe in Python. This is also why the
staleness hash (below) includes the `check_*.F90`/`scf_*.F90` source itself:
if that recipe ever changes, the committed reference is correctly flagged
stale.

The SCF comparisons (`rks`, `scf_robust`, `rimp2`) are the one place a
*converged* density matters, and there it is pyscf's own, not terco's:
`rks.ref` commits pyscf's converged energy AND its converged density (in
terco's AO order); the Fortran side's "at-density control" calls
`trc_scf_run` with `max_iter=1` and that committed density as `dguess`,
which runs exactly one Roothaan iteration -- `src/scf/trc_scf_driver.F90`'s own
"energy at this density" comment -- so the returned energy is terco's
functional evaluated at pyscf's density, equivalent to pyscf's own
`energy_tot(dm0)`, never an iterated answer.

## File format

Plain text, read by `test/pyscf_ref_io.F90` (see that file's own header for
the exact grammar): comment lines start with `#`; everything else is
`TAG value...`, one logical record per tag -- `INT`/`REAL`/`STR` scalars,
`VECTOR`/`MATRIX`/`MATRIX3` blocks (shape on the tag line, then that many
`%.17e` floats, column-major/Fortran order, list-directed so line breaks
inside a block don't matter). Every matrix here is already in **this
build's own AO order and normalisation** -- `generate.py` does the
shell-order permutation (pyscf's `cart=True` groups a given atom's shells
by ascending angular momentum; terco keeps a JSON's literal order) and the
per-AO rescale (pyscf's Cartesian d/f off-diagonal components are not all
unit-normalised the way this build's are) once, offline, so the Fortran
reader diffs straight across with no permutation of its own.

Every file's leading comment block also carries:

```
# generator: test/pyscf_refs/generate.py CASES["<case>"] -> <gen_* function>
# hash-scheme: case-ast-v1
# hash: <sha256 of that case's code and every file on the "files:" line>
# files: <path> <path> ...
```

The hash covers only what wrote this reference (`refhash.py`, standard
library only):

- its `CASES` entry in `generate.py` (the dict value, so a lambda's
  arguments count) and the transitive closure of every module-level
  function, class and constant of `generate.py` that entry refers to,
  found by walking the parsed source with `ast`. Each definition is
  serialised without comments, docstrings or line positions, so those edits
  never restale anything, and the serialisation is the same on Python 3.9
  through 3.14;
- the bytes of the files on `files:` -- geometry, basis JSON, the Fortran
  source that builds a fixed density or grid.

So editing one `gen_*` function restales only the references it writes,
and editing a shared helper (`converge_tight`, `RefWriter`, `run_bin`)
restales exactly the references whose generators reach it.
`check_staleness.py` (standard library only -- no pyscf, no numpy) imports
the same `refhash.py`, recomputes the hash from the `CASES` key on the
`generator:` line and the paths on `files:`, and names the case to rerun
when it does not match. Registered as `terco/check_pyscf_refs_fresh`, it
runs in every CI config; `terco/check_refhash` (`test_refhash.py`) tests
which edits restale which references.

The references written before this scheme hashed the whole of
`generate.py`. `migrate_ref_headers.py` moved every one of them to it
without regenerating anything: it checked each was fresh under the old
hash, then rewrote only its generator/hash/files header lines.

## Two tiers: pyscf anchors and terco goldens

Every method checked here is checked twice.

- **Anchors** -- the pyscf references in this directory. They say terco is
  *right*: an independent code agrees with it, at the physics tolerance of
  each comparison (`1e-8` Eh on a closed-shell SCF energy, for example).
  They are frozen. Regenerate one only for a new method, or when that
  case's own generator code or inputs change (`check_staleness.py` names
  exactly which); never to absorb a terco change. No anchor tolerance is
  ever loosened.
- **Goldens** -- `test/goldens/<case>.gold`, terco's own numbers from a
  terco run. They say terco has *not moved*: the same check also records its
  quantities (`test/golden_io.F90`) and `terco/golden/<case>`
  (`test/goldens/compare_goldens.py`) compares them at a tolerance set from
  the measured spread between CI's compilers, far below the anchor's. A
  refactor that shifts an energy by `1e-9` Eh passes every anchor and fails
  its golden.

When a golden fails:

- **unintended** -- it is a regression; fix the code.
- **intended** (a new algorithm, a different but equally valid order of
  operations) -- first confirm the anchors still pass, since that is what
  says the new numbers are still right. Then regenerate the goldens from a
  terco build, which runs terco only and never pyscf:

  ```
  python3 test/goldens/regen_goldens.py --build build [--build build-lmax3]
  ```

  and say in the commit what moved and by how much. A case that only runs
  at `TERCO_LMAX=3` comes from the second build.

The goldens cover the SCF and energy checks: `scf_rks`, `scf_df`,
`scf_robust`, `scf_rimp2` (both bases), `scf_pcm` (both sets), and the
fixed-density `rys_fock` and `erf_fock` builds (every kernel path). For
the Fock builds, which have no single energy, the golden is three scalars
of each matrix: `tr(D G)`, its Frobenius norm, and a sum with fixed uneven
weights, which moves when any element does.

Tolerances (`TOLERANCES` in `compare_goldens.py`), absolute, by quantity:

| quantity | key prefix | tolerance |
|---|---|---|
| total energies, Eh (SCF, at-density, E_diel, RI-MP2 components) | `E_` | `1e-10` |
| `tr(D G)` of a fixed-density Fock build | `TR_D` | `1e-10` |
| Frobenius norm of a fixed-density Fock build | `NORM_` | `1e-10` |
| fixed-weight element sum of a fixed-density Fock build | `WSUM_` | `1e-10` |

These are provisional, pending the spread between CI's compilers (gcc-15,
intel, nvidia-hpc, the OpenACC-to-OpenMP translation). The goldens come
from gfortran 15. Two gfortran 15 builds, at `TERCO_LMAX=2` and `=3`,
agree on every shared quantity bit for bit.

## Mapping from the old `*_ref.py` scripts

| old script | committed reference(s) | Fortran comparison |
|---|---|---|
| `ao_ref.py` | `ao.ref` | `check_ao.F90` |
| `mult_ref.py` | `mult.ref` | `check_mult.F90` |
| `xc_energy_ref.py` | `xc_energy.ref` | `check_xc_energy.F90` |
| `vv10_ref.py` | `vv10.ref` | `check_vv10_vs_pyscf.F90` |
| `rys_grad_ref.py` | `rys_grad_631gs.ref`, `rys_grad_ccpvdz.ref` | `check_rys_grad_vs_pyscf.F90` |
| `rys_fock_ref.py` | `rys_fock_water_ccpvdz.ref`, `rys_fock_water_ccpvtz.ref`, `rys_fock_formamide_def2tzvp.ref` | `check_rys_fock_vs_pyscf.F90` |
| `erf_fock_ref.py` | `erf_fock_water_631gs.ref`, `erf_fock_water_ccpvdz.ref`, `erf_fock_formamide_ccpvdz.ref` (shared by the hgp/nosym/rys kernel-path binaries) | `check_erf_fock_vs_pyscf.F90` |
| `df_vs_pyscf_ref.py` | `df_vs_pyscf.ref` (local-only, see above) | `check_df_vs_pyscf.F90` |
| `rks_ref.py` | `rks.ref` | `scf_rks.F90` |
| `rimp2_ref.py` | `rimp2.ref`, `rimp2_ccpvdz.ref` | `scf_rimp2.F90` |
| `scf_robust_ref.py` | `scf_robust.ref` | `scf_robust.F90` |
| `xc_ref.py` | `../reference_xc.dat` (unchanged location/format) | `check_xc.F90` (unchanged) |
| (new, PCM) | `pcm_cavity_water.ref` | `check_pcm_cavity.F90` |
| (new, PCM) | `pcm_ints_water_631gs.ref`, `pcm_ints_water_ccpvtz.ref` (LMAX=3 build) | `check_pcm_ints.F90` |
| (new, PCM) | `pcm_scf_main.ref`, `pcm_scf_f.ref` (LMAX=3 build) | `scf_pcm.F90` |

Every `terco/*_vs_pyscf_ref`/`*_ref` ctest entry these scripts used to back
is gone; the single `terco/check_*`/`scf_*` entry now does both halves
(build/compute, then compare) in one Fortran process.

## Tolerances

Unchanged from the scripts above, carried into the Fortran comparisons as
literal constants (or, where a script varied it by case -- `rys_fock`'s
`tol` argument -- as a CLI argument test/CMakeLists.txt still passes per
case). No gate was loosened; a few were tightened where the measured
agreement left ample margin (see git history for the measurements, and this
project's own commit history on the pre-existing Python scripts for the
original tightenings, e.g. water/cc-pVTZ's `rys_fock` gate at 1e-10).

## Left open

- References are generated for `TERCO_CART_ORDER=lexicographic` only (see
  above) -- a `gamess`-order build's `check_mult`/`check_ao`/
  `check_xc_energy` pyscf comparisons are not validated against this
  commit's references. CI never builds `gamess` order, so this does not
  affect CI.
- `erf_fock_*.ref` is shared across the hgp/nosym/rys kernel-path binaries
  by construction (the fixed density and pyscf's J/K do not depend on which
  terco path is under test) -- documented in `generate.py`'s own
  `gen_erf_fock`, not re-derived per path.
