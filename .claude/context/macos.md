# Building and testing on macOS (Apple Silicon)

Verified on an Apple M5 Pro (6 performance + 12 efficiency cores, 18 logical),
macOS 26.5.1, with no NVIDIA GPU. Everything marked [GPU] elsewhere cannot run here.

## Install (no sudo needed)

```sh
brew install gcc@15 cmake ninja openblas python@3.12
/opt/homebrew/bin/python3.12 -m venv ~/.venvs/terco
~/.venvs/terco/bin/pip install fortitude-lint==0.9.2 'pyscf==2.14.*' mpmath==1.3.0 pytest
```

Use Homebrew's Python for the venv. Apple's `/usr/bin/python3` is 3.9, and the
generators (`scripts/md/hermite.py` uses `str | None`) need Python 3.10 or
later; with 3.9 the configure step fails in `gen_md_j.py`.

Versions installed: gfortran 15.3.0 (gcc@15), cmake 4.4.5, ninja 1.13.2,
openblas 0.3.34, python 3.12.16, fortitude-lint 0.9.2, pyscf 2.14.0,
mpmath 1.3.0, pytest 9.2.0, numpy 2.5.4, scipy 1.18.1. `brew install gcc@15`
also pulled gcc 16.2.0 as a dependency; use `gfortran-15` explicitly. All
pinned versions installed as asked.

## Build and test

The `macos-host` preset (gfortran-15, Ninja, Release, GPU and MPI off, `Python3_EXECUTABLE`
set to the venv):

```sh
export PATH=/opt/homebrew/bin:$PATH
cmake --preset macos-host
cmake --build build --parallel
ctest --test-dir build --output-on-failure -LE LONG --parallel 4
~/.venvs/terco/bin/python test/goldens/compare_goldens.py --report build/test/golden_out
~/.venvs/terco/bin/fortitude check
(cd test && ~/.venvs/terco/bin/python pyscf_refs/check_staleness.py)
```

BLAS and LAPACK: `find_package(BLAS)` picks Apple Accelerate (via the
Command Line Tools SDK) and the whole suite passes with it, so no CMake change
was needed. openblas is installed but not used; to try it, add
`-DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas -DBLA_VENDOR=OpenBLAS`
(not run in S0).

Nothing in CMake needed an `if(APPLE)` guard: the shared library, `link_c.c` test,
`check_static_terco_opt.cmake` (only registered with DL-FIND) and the MPI-off
`pic-mpi` path all worked unchanged.

## OpenMP build

Use a scratch worktree because the converters rewrite `src/` in place:

```sh
git worktree add ../terco-omp -b scratch-omp <base> && cd ../terco-omp
export PATH=$HOME/.venvs/terco/bin:/opt/homebrew/bin:$PATH
python3 tools/dc_locality_lint.py src
python3 -m pytest tools/test_acc_to_omp.py -q
python3 tools/dc_to_omp.py  --target cpu --write src
python3 tools/acc_to_omp.py --target cpu --translate-atomics --write src
cmake --preset macos-omp
python3 tools/dc_to_omp.py  --target cpu --write build-omp/generated_*
python3 tools/acc_to_omp.py --target cpu --translate-atomics --write build-omp/generated_*
cmake --build build-omp --parallel
OMP_NUM_THREADS=4 ctest --test-dir build-omp --output-on-failure -LE LONG
```

The `macos-omp` preset only sets `-fopenmp` and `build-omp`; it does not convert the
source, so it is meaningless on an unconverted tree. The gcc 15 OpenMP runtime (libgomp) is
used, and no `-foffload=disable` was needed.

## Results (S0, base 663d956)

| Step | Wall time |
|---|---|
| configure (generates Rys and MD J kernels) | 46 s |
| host build, `--parallel` | 65 s |
| host `ctest -LE LONG --parallel 4` | 130 s, 227/227 passed |
| OpenMP build | 66 s |
| OpenMP `ctest -LE LONG` at 4 threads, run serially | 356 s, 227/227 passed |

19 tests report Skipped in both builds. They are the `TERCO_LMAX=3` cases (ccpvtz,
def2-tzvp, ccpvqz, `check_nosym_f`, `scf_pcm_f`, ...) that skip themselves at
the default `TERCO_LMAX=2`, as on Linux CI. No test is skipped for a macOS reason.
Golden spread report: largest spread about 9e-13 (`TR_DK`), no change. The 11k
"missing terminating '" warnings come from apostrophes inside `!` comments of
generated sources under the preprocessor (not compared with a Linux build).

Scaling of the OpenMP build, B3LYP/def2-SVP on `test/gly10.xyz` (73 atoms), energy
-2154.9207372146 Eh at every thread count (differences 4e-12 Eh, thread
reduction order):

| `OMP_NUM_THREADS` | wall time |
|---|---|
| 1 | 325 s |
| 4 | 131 s (2.5x) |
| 16 | 73 s (4.5x) |

Sixteen threads exceed the 6 performance cores, so part of them land on
efficiency cores and the speedup is well short of 16. RHF/def2-SVP benzene is too
small to scale (2.4 s, 1.2 s, 0.85 s).

## Skipped on macOS

None.
