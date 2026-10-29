# CI jobs and how to reproduce them

`.github/workflows/ci.yml` runs these on every PR. None has a GPU: every
job builds for the host, where `do concurrent` is a serial loop and OpenACC
directives are inert.

| Job | What it checks |
|---|---|
| `host gcc-15`, `host intel-2025.0`, `host nvidia-hpc-latest` | host build and the full ctest suite, then the "Golden spread" report |
| `host gcc-15 (TERCO_LMAX=3)`, `(TERCO_DF_LMAX=3)` | the f-function and f/g-auxiliary builds; regenerates kernels at configure |
| `generated sources are up to date` | reruns the generators and requires no diff in `src/generated/` |
| `Rys root tables are up to date` | regenerates the Rys root tables with mpmath 1.3.0 |
| `OpenACC to OpenMP conversion` | converts the source to OpenMP, builds with `-fopenmp`, runs ctest at 4 threads |
| `fortitude lint` | `fortitude check`, fortitude-lint 0.9.2 |
| `fpm (gcc 15)` | the fpm build (often skipped) |

## Host build

```sh
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DTERCO_ENABLE_GPU=OFF -S . -B build
cmake --build build --parallel
ctest --test-dir build --output-on-failure --parallel 2 --no-tests=error
python3 test/goldens/compare_goldens.py --report build/test/golden_out
```

Add `-DTERCO_LMAX=3` or `-DTERCO_DF_LMAX=3` for the f builds. They take far
longer to compile.

Pull requests and pushes run ctest with `-LE LONG`. Tests labelled `LONG`
in `test/CMakeLists.txt` (`check_xc_grad_fd`, over half the suite's time) run
only in the nightly scheduled run on main (13:00 UTC, midnight AEDT) and in a
run started by hand (`workflow_dispatch`). Run them locally with `ctest -L LONG`
when you touch the XC gradient.

## Generators

```sh
python3 scripts/gen_perclass.py --lmax 2 --split src/generated
python3 scripts/gen_1e.py --lmax 2 -o src/generated/trc_1e_kernels.F90
python3 scripts/gen_1e.py --lmax 2 --multipole src/generated/trc_mult_kernels.F90
python3 scripts/gen_1e.py --lmax 2 --pcm src/generated/trc_pcm_kernels.F90
python3 scripts/gen_df.py --lmax 2 --lmax-aux 4 --unroll-max-comp 81 --split src/generated
git diff --exit-code src/generated
```

The CI workflow is the source of truth if these drift.

## OpenMP port

This is the job most often red for reasons the host jobs don't see. Work on
a scratch copy or worktree, because the converters rewrite `src/` in place:

```sh
python3 tools/dc_locality_lint.py src
python3 -m pytest tools/test_acc_to_omp.py -q
python3 tools/dc_to_omp.py  --target cpu --write src
python3 tools/acc_to_omp.py --target cpu --translate-atomics --write src
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DTERCO_ENABLE_GPU=OFF \
      -DCMAKE_Fortran_FLAGS=-fopenmp -S . -B build-omp
# The Rys and MD J/K kernels are generated into build-omp at configure time,
# after the conversion of src/, so convert them too.
python3 tools/dc_to_omp.py  --target cpu --write build-omp/generated_*
python3 tools/acc_to_omp.py --target cpu --translate-atomics --write build-omp/generated_*
cmake --build build-omp --parallel
OMP_NUM_THREADS=4 ctest --test-dir build-omp --output-on-failure
git checkout -- src        # undo the conversion
```

The job uses the runner's stock gfortran, which may be older than the
toolchain used elsewhere. Things that have failed only here:

- Generated kernels in the build tree (`build-omp/generated_rys`,
  `generated_mdj`): the converters only see the paths they are given.
- A compiler with an nvptx offload plug-in (Gadi's gcc 15.1) fails at link
  time under `-fopenmp` ("referenced in offloaded code"); add
  `-foffload=disable`. The runner's stock gfortran does not need it.
- Source the converters didn't scan. `walk()` in both tools lists the file
  extensions it converts. `*.inc` is added by the gradient work (#83). Any
  other new extension holding `do concurrent` or `!$acc` code must be added
  to both tools.
- A hand-written `!$omp parallel default(none)` missing a variable, dummy
  arguments included.
- Races that are invisible on one thread. That is why the job runs at 4
  threads.

## fortitude

```sh
pip install fortitude-lint==0.9.2
fortitude check
```

Rules live in `fpm.toml` under `[extra.fortitude.check]`. `test/` and
`src/generated/` are excluded. The common hits are S082 (one statement per
line), S241 (double quotes), S262, S221 (`result()`) and C031 (named
constants for array sizes).

## pyscf references

```sh
cd test && python3 pyscf_refs/check_staleness.py
python3 pyscf_refs/test_refhash.py
```

Regenerating needs pyscf and a CPU build; see `test/pyscf_refs/README.md`.
