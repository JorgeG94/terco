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
which runs exactly one Roothaan iteration -- `src/trc_scf.F90`'s own
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
# hash: <sha256 of generate.py + every file on the "files:" line>
# files: <path> <path> ...
```

`check_staleness.py` (standard library only -- no pyscf, no numpy)
recomputes that hash from exactly the paths on `files:` and compares; it
never hardcodes per-case knowledge of what feeds into which reference, so
a case that starts reading a new input is automatically covered. Registered
as `terco/check_pyscf_refs_fresh`, it runs in every CI config.

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
| `df_vs_pyscf_ref.py` | `df_vs_pyscf.ref` | `check_df_vs_pyscf.F90` |
| `rks_ref.py` | `rks.ref` | `scf_rks.F90` |
| `rimp2_ref.py` | `rimp2.ref`, `rimp2_ccpvdz.ref` | `scf_rimp2.F90` |
| `scf_robust_ref.py` | `scf_robust.ref` | `scf_robust.F90` |
| `xc_ref.py` | `../reference_xc.dat` (unchanged location/format) | `check_xc.F90` (unchanged) |

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
