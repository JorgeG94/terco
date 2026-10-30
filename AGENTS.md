# AGENTS.md

Guidance for AI coding agents (Claude Code, Codex, Copilot, ...) working on
terco. Humans are welcome to read it too. Keep it under 200 lines; put
anything longer under `.claude/` and link it from here.

## What terco is

A GPU electron-repulsion-integral and SCF library in **standard Fortran**.
Compute kernels are `do concurrent`, compiled by nvfortran with
`-stdpar=gpu -acc -gpu=mem:separate`. OpenACC is used only to move data.
The project's claim is that standard Fortran is competitive with hand-written
CUDA, so every change must keep two things true: **the numbers are right**,
and **the code stays portable standard Fortran**.

It covers HF and DFT (R/U), Fock builds over binned shell pairs, density
fitting, RI-MP2, PCM, VV10, analytic gradients and response (CPHF/CPKS), a
Fortran API, a `bind(c)` ABI and a `terco` CLI. See `README.md`.

## Layout

| Path | What lives there |
|---|---|
| `src/core/` | kinds glue, Boys, Cartesian/spherical tables, linalg, errors, MPI device binding |
| `src/basis/`, `src/grid/` | basis JSON reader, decontraction; DFT grids and Becke/Stratmann partitions |
| `src/eri/` | screening, shell-pair bins, HGP/Rys drivers, the McMurchie–Davidson J and exchange engines (`trc_mdj`, `trc_mdg`), Fock and gradient launchers |
| `src/scf/` | SCF driver, settings, state, DIIS, guesses, DF-JK, gradients, CPHF |
| `src/xc/`, `src/post_hf/`, `src/solvation/` | XC integration and functionals, VV10; RI-MP2; PCM |
| `src/api/` | `trc_api` (Fortran) and `trc_capi` (C ABI); `include/` holds the C interfaces |
| `src/opt/` | the geometry optimiser (DL-FIND via libdlfind), built as `libterco_opt` only with `-DTERCO_ENABLE_DLFIND=ON`; its own C ABI is `include/terco_opt.h` |
| `src/generated/` | **generator output: never edit by hand** |
| `scripts/` | the kernel generators (`gen_*.py`), see `docs/GENERATORS.md` |
| `tools/` | OpenMP converters, `pic_api.py`, XC registry and translators, lint helpers |
| `app/` | the `terco` CLI (`main.F90`, `trc_app.F90`) and its MPI checks |
| `test/` | every check program, `pyscf_refs/` (frozen anchors), `goldens/` (terco's own numbers) |

## Build and test

Host build (what CI runs, no GPU needed; gfortran 15 is the floor):

```sh
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DTERCO_ENABLE_GPU=OFF -S . -B build
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

GPU build: `-DCMAKE_Fortran_COMPILER=nvfortran -DTERCO_GPU_ARCH=cc70` (or your
arch). `TERCO_LMAX=3` / `TERCO_DF_LMAX=3` build f functions and regenerate
kernels. The host build is a validation target and must give the same
numbers as the device build: code that is correct only on the GPU is a bug.
Debug on the CPU build first (faster, better diagnostics). Use the GPU build
to confirm and to measure. GPU-only failures and build gotchas:
`.claude/context/gpu-pitfalls.md`.

Before pushing, run what CI runs on the files you touched:

- `fortitude check` (fortitude-lint **0.9.2**; rules in `fpm.toml`; `test/`
  and `src/generated/` are excluded). Fix the code; don't add ignores.
- `python3 test/pyscf_refs/check_staleness.py` from `test/` if you touched
  anything a `.ref` lists (see below).
- The OpenMP port, if you touched `do concurrent` code or `!$omp`/`!$acc`
  directives. CI compiles it with `-fopenmp` and runs ctest at 4 threads.
  Reproduce it with `.claude/context/ci.md`.

## Non-negotiable rules

The full review checklist is `.claude/skills/pr-review/SKILL.md`. Read it
before any non-trivial change. The rules that most often break builds:

**Correctness**
- Energies and integrals change only with a reference that justifies the
  new number (pyscf, libcint or a stored value). "Converges" is not
  validation. Never loosen a tolerance or edit a reference literal to make
  a test pass.
- libcint's `common_fac_sp` normalisation is applied exactly once, in basis
  construction. Å→Bohr is CODATA 2018 (0.529177210903).
- A coupled-perturbed solve uses `density_screen=.false.`.
- Iterative solvers report non-convergence on a breakdown or NaN; never let
  a NaN fall through a `resid >= tol` loop test as "converged".

**GPU / `do concurrent`**
- No CUDA, no OpenACC compute constructs, no `async`. OpenACC is limited to
  `enter`/`exit data`, `update`, `host_data`, `routine seq` and `atomic`.
- Every `enter data` has an `exit data` on every path out, including
  error returns. Every array a kernel reads must be mapped: under
  `mem:separate` a missing mapping is a run-time illegal address, and the
  host build will not show it.
- A `do concurrent` body calls only intrinsics and same-module `pure`
  `!$acc routine seq` item routines, with explicit-shape dummies and
  per-thread scratch declared inside the item routine.
- A procedure or program, test programs included, that makes a collective
  MPI call must have no internal (`contains`) procedures (an NVHPC 26.5 +
  HPC-X defect). Put helpers in modules.
- Always `%release()` a `trc_linalg_t` before it goes out of scope (it has no
  finalizer). Prefer one flattened index over a multi-index `do concurrent`
  with a small run-time outer extent.
- A hand-written `!$omp parallel default(none)` must list every variable it
  reads. Only the OpenMP CI job compiles it, so check it yourself.

**Generated code**
- `src/generated/` comes from `scripts/`. Change the generator, rerun it and
  commit both. CI regenerates and diffs.
- Source the OpenMP converters don't scan (`walk()` in `tools/dc_to_omp.py`
  and `tools/acc_to_omp.py`) breaks the OpenMP job. Check any new file type.
- State the compile-time cost of a generator change that grows the heavy
  classes; the per-class kernels already take 15–20 minutes to build.
- A Rys class above d needs a measured chunk plan (`FOCK_CHUNK_PLAN` in
  `scripts/gen_rys.py`); the generator refuses one without. Race the shapes,
  don't guess. Read register spills from the linked `libterco.so`
  (`cuobjdump -res-usage`), never the `.o`: `-rdc` reports 0 there.

**Fortran conventions**
- Kinds come from `pic_types` (`dp`, `int32`, `int64`, ...). No
  `iso_fortran_env` kinds, no `real(8)`, no `1.0d0`, no `iso_c_binding` kinds
  away from the C boundary.
- Use pic instead of local copies: `pic_timer` (not `system_clock`),
  `pic_strings: to_string` (not `write(buf,"(i0)")`), `pic_sorting`,
  ... The list and the pinned tag are in `docs/pic_api.md`.
  Device-side `pure` helpers are exempt.
- No bare `print` or `write(*,...)` to the console: everything goes through
  `pic_logger`'s `global_logger` (`info`, `performance`, `warning`, ...).
  Format into a buffer with an internal `write` if you need a fixed layout.
- Docs are **FORD only, never Doxygen**. No `!>`, `!<` or `@param`. A `!!`
  block goes *after* its entity, after a declaration or indented under a
  procedure header. A `!!` block above a declaration documents the one
  before it. Executable lines get plain `!` comments, never `!!`.
- `trc_` module prefix, match the surrounding code, comments explain *why*
  (often with the measurement), not history.

## Tests and references

- New behaviour needs a check under `test/` that fails without the change,
  against an external reference where one exists. A test that only checks
  that a run completes is not a test.
- **New checks use [test-drive](https://github.com/fortran-lang/test-drive)**
  (`cmake/modules/Findtest-drive.cmake`): a module whose `collect_*` routine
  returns `new_unittest` entries, each a `check` with a message, run by
  `run_testsuite`; link it with `test-drive::test-drive`. `test/check_si.F90`
  is the pattern. Most existing checks are still plain programs; port one
  when you touch it.
- **Two tiers** (`test/pyscf_refs/README.md`):
  - *Anchors*: `test/pyscf_refs/*.ref`, written offline by
    `generate.py` (pyscf never runs in CI). Each `.ref` carries a
    per-case hash (`hash-scheme: case-ast-v1`) over its `CASES` entry, the
    `generate.py` code it reaches and every file on its `# files:` line.
    Editing any of those makes `terco/check_pyscf_refs_fresh` fail.
    Re-stamp the hash only if you can show that the probe the test writes
    is unchanged. Otherwise regenerate with `generate.py <case>`, never by
    hand.
  - *Goldens*: `test/goldens/*.gold`, terco's own numbers at tight,
    measured tolerances (`compare_goldens.py`). Regenerate with
    `regen_goldens.py` for an intended numerical change that still passes
    the anchors.
- Performance claims state the system, basis, thresholds, GPU and baseline,
  and show unchanged energies. Per-class claims come from nsys/ncu, not wall
  time alone.

## Git and PR conventions

- Branch off `main`. Merge with merge commits, not squash. Never rewrite
  history someone else has pulled.
- Many PRs are stacked: base each one on the branch below it and say so in
  the body. Before merging, check overlap with the other open PRs. A clean
  textual merge can still fail to compile across a refactor. Details:
  `.claude/context/workflow.md`.
- No AI-session links in commits or PR bodies. A `Co-Authored-By:` trailer
  is fine. Commit messages say what changed and why. See SKILL.md section 9
  for the rest.
- Never merge without the maintainer's say-so.

## More context

- `.claude/skills/pr-review/SKILL.md`: the full review checklist (use the
  `pr-review` skill).
- `.claude/context/ci.md`: each CI job and how to reproduce it locally.
- `.claude/context/workflow.md`: writing tests that can fail, refreshing
  pyscf references, stacked PRs, and running long jobs as an agent.
- `.claude/context/gpu-pitfalls.md`: nvfortran and `mem:separate` traps,
  MPI and build gotchas.
- `docs/GENERATORS.md`: what each generator emits and how to extend it.
- `docs/pic_api.md`: everything pic provides at the pinned tag.
- `docs/functionals.md`: every registered functional (generated).
- `test/pyscf_refs/README.md`: regenerating anchors and the per-case hash.
