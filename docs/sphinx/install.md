# Building terco

terco builds with CMake. There is also an fpm manifest for the shipped
configuration (`TERCO_LMAX=2`); see the end of this page.

## CPU (host) build

Needs **gfortran 15 or newer** (the floor, not a preference: `do concurrent`
locality specifiers are Fortran 2018 and landed there), ifx, or nvfortran.
This is the configuration that validates terco without a GPU -- `do
concurrent` runs as an ordinary serial loop and every OpenACC directive is
an inert comment, so a passing host build means the integrals are right
independent of the offload model.

```sh
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DTERCO_ENABLE_GPU=OFF -S . -B build
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## GPU build

Needs **NVHPC**: `-stdpar=gpu` is what turns `do concurrent` into device
kernels, and `TERCO_ENABLE_GPU=ON` with any other compiler is refused
outright rather than silently building a slow CPU library.

```sh
cmake -S . -B build -DCMAKE_Fortran_COMPILER=nvfortran -DTERCO_GPU_ARCH=cc70
cmake --build build -j
ctest --test-dir build --output-on-failure
```

`TERCO_GPU_ARCH` is spelled as nvfortran spells it (`cc70`, `cc80`, ...).

## Build options

From the top-level `CMakeLists.txt`:

| Option | Default | Meaning |
|---|---|---|
| `TERCO_LMAX` | `2` | four-centre angular-momentum ceiling (0=s, 1=p, 2=d, 3=f, 4=g). Anything but 2 regenerates the kernels and needs Python 3; 3 and 4 add real build time (minutes to tens of minutes). |
| `TERCO_DF_LMAX` | empty (follows `TERCO_LMAX`) | density-fitting kernels' orbital angular-momentum ceiling |
| `TERCO_DF_LMAX_AUX` | `4` | density-fitting kernels' auxiliary-basis angular-momentum ceiling |
| `TERCO_CART_ORDER` | `lexicographic` | within-shell Cartesian component order: `lexicographic` (libcint's own order) or `gamess` |
| `TERCO_GPU_ARCH` | `cc70` | target GPU architecture, as nvfortran spells it |
| `TERCO_CUDA_VERSION` | empty (compiler default) | CUDA version to build the device code against, e.g. `12.9` |
| `TERCO_BUILD_SHARED` | `ON` | build `libterco.so` (a host built with another compiler links against this); `OFF` for a static library |
| `TERCO_REGENERATE` | `OFF` | regenerate the committed kernels in `src/generated/` in place |
| `TERCO_FORCE_REGENERATE` | `OFF` | regenerate the build-tree kernels once, ignoring the generation stamps |
| `TERCO_ENABLE_MPI` | `OFF` | split the Fock build over MPI ranks, one GPU each |
| `TERCO_MPI_DEVICE` | follows `TERCO_ENABLE_GPU AND TERCO_ENABLE_MPI` | hand MPI device buffers directly (a CUDA-aware MPI); overridable at run time with `TERCO_MPI_DEVICE=0/1` |
| `TERCO_ENABLE_DLFIND` | `OFF` | build the `--opt` geometry optimiser against libdlfind (LGPL-3.0-or-later; see below) |
| `TERCO_ENABLE_TESTING` | `ON` | build the correctness checks |
| `TERCO_LONG_TESTS` | `OFF` | register the checks labelled `LONG` instead of `DISABLED` |
| `TERCO_BUILD_APP` | `ON` | build and install the `terco` command-line program |
| `TERCO_ENABLE_MULTICORE` | `OFF` | host build with `-stdpar=multicore -acc=multicore` (needs NVHPC; mutually exclusive with `TERCO_ENABLE_GPU`) |
| `TERCO_NVTX` | `OFF` | NVTX ranges for Nsight Systems (NVHPC builds) |

## MPI

```sh
cmake -S . -B build -DCMAKE_Fortran_COMPILER=nvfortran -DTERCO_GPU_ARCH=cc70 \
      -DTERCO_ENABLE_MPI=ON
cmake --build build -j
mpiexec -np 4 terco water.xyz -b 6-31g*     # one GPU per rank
```

With `TERCO_ENABLE_MPI=OFF` (the default), the library still links against
pic-mpi, but through its single-rank serial backend, so the same code path
runs either way: without MPI you simply have rank 0 of 1.

## Geometry optimisation (DL-FIND)

The optimiser (`--opt`) is **off by default**, for a licensing reason, not a
size one: it is driven by DL-FIND through libdlfind, both
LGPL-3.0-or-later, and terco itself is MIT. With `TERCO_ENABLE_DLFIND=OFF`
(the default), neither libdlfind nor its header is built at all, and
`terco --opt` is refused with a clear message.

```sh
cmake -S . -B build -DCMAKE_Fortran_COMPILER=nvfortran -DTERCO_GPU_ARCH=cc70 \
      -DTERCO_ENABLE_DLFIND=ON
cmake --build build -j
```

With this on, the optimiser is built as its **own** shared library,
`libterco_opt` (`terco::terco_opt`, its own C ABI at `include/terco_opt.h`),
never linked into `libterco` itself -- so the LGPL relinking obligation
stays confined to that one library, and a consumer who does not need
geometry optimisation never pulls in DL-FIND at all.

## fpm

```sh
tools/fpm_generate.sh                     # the Rys and MD kernels, into src/generated/fpm/
fpm build --compiler mpifort --flag "-ffree-line-length-none"            # gfortran
fpm build --compiler mpifort --flag "-stdpar=gpu -acc -gpu=cc70,mem:separate"  # nvfortran
```

fpm builds the **shipped configuration and only that one** (`TRC_LMAX=2`):
it cannot run a code generator, so `src/generated/` holds the committed HGP
kernels and `tools/fpm_generate.sh` generates the rest (Rys energy and
gradient kernels, the McMurchie-Davidson J/K kernels) that CMake would
otherwise generate at configure time. Any other angular-momentum ceiling
needs CMake.

`--opt` under fpm needs two edits to `fpm.toml` (adding `TERCO_ENABLE_DLFIND`
to `preprocess.cpp.macros` and a `[dependencies] libdlfind` entry pinned to
terco's libdlfind fork) -- see the top-level
[README](https://github.com/JorgeG94/terco#opt-dl-find-under-fpm) for the
exact snippet, which is kept there next to the licensing note above.
