# Building terco

## What you need

- **An NVIDIA GPU** (developed and benchmarked on V100, A100, H100, H200 and
  GH200).
- **The NVIDIA HPC SDK** (free download from NVIDIA), which provides
  `nvfortran` and, under `comm_libs/mpi/bin`, an MPI built with it. terco is
  developed with NVHPC 26.5.
- **CMake 3.28 or newer**, or **fpm** (0.13).

No CUDA toolkit, no other libraries: the SDK has everything, and CMake or
fpm fetch terco's Fortran dependencies (pic, json-fortran, toml-f) from
GitHub at build time, so the build machine needs network access.

The quickest GPU build:

```sh
export PATH=$NVHPC_ROOT/compilers/bin:$PATH      # nvfortran
cmake -S . -B build -DCMAKE_Fortran_COMPILER=nvfortran -DTERCO_GPU_ARCH=cc90
cmake --build build --parallel 16
cmake --install build --prefix $HOME/install/terco
```

`TERCO_GPU_ARCH` is your GPU's compute capability as nvfortran spells it:
`cc70` V100, `cc80` A100, `cc90` H100/H200/GH200. Without `TERCO_LMAX` this is
the default build, up to d functions; `-DTERCO_LMAX=3` (f) or `4` (g) needs
Python 3 to generate the extra kernels and takes much longer to compile.

### MPI and compiler wrappers

The CMake build needs no MPI unless you turn on `TERCO_ENABLE_MPI` (below).
The **fpm** build always does: one of its dependencies (pic-mpi) compiles MPI
backends, so fpm must be given an MPI compiler wrapper, not the bare
compiler, and that wrapper must wrap the compiler you want:

- With the HPC SDK alone, use its own wrapper:
  `export PATH=$NVHPC_ROOT/comm_libs/mpi/bin:$PATH` and `--compiler mpifort`.
- On a cluster, load the MPI module built for the NVHPC compiler (not one
  built for GCC or Intel) and use its `mpifort`. `mpifort --version` must
  print nvfortran.
- On a Cray system, the programming environment's wrapper is `ftn`
  (`module load PrgEnv-nvidia`), so pass `--compiler ftn`.

If CMake's `TERCO_ENABLE_MPI=ON` finds the wrong MPI (one built for another
compiler), unset `MPI_ROOT`/`MPI_HOME` in the environment and point
`MPI_Fortran_COMPILER` at the right wrapper.

## Release tarballs for fpm

CMake (above) is the recommended build. For those who prefer fpm, each
release also ships `terco-<version>-lmax2.tar.gz`, `-lmax3` and `-lmax4`:
the source with every kernel already generated for that angular-momentum
ceiling, so fpm needs nothing else. Pick the highest ceiling your basis sets
need (d, f or g), then:

```sh
tar xf terco-0.1.0-lmax3.tar.gz && cd terco-0.1.0-lmax3
export PATH=$NVHPC_ROOT/comm_libs/mpi/bin:$NVHPC_ROOT/compilers/bin:$PATH
mkdir -p blas && ln -sf $NVHPC_ROOT/compilers/lib/libblas.so blas/libopenblas.so
fpm install --prefix $HOME/install/terco/0.1.0 --compiler mpifort --profile release \
  --flag "-stdpar=gpu -acc -gpu=cc90,mem:separate -DTRC_DEVICE" \
  --link-flag "-L$PWD/blas -llapack -cudalib=cublas,cusolver"
export TERCO_BASIS_PATH=$PWD/basis_sets          # fpm does not install data
```

The link flags are what CMake adds for you: a BLAS dependency asks for
`-lopenblas`, which the SDK ships as `libblas.so` (hence the symlink; skip it
if an OpenBLAS is already on the linker path), and cuBLAS/cuSOLVER come in
through `-cudalib`. At run time, the MPI's own libraries must come first on
`LD_LIBRARY_PATH`; loading the SDK's or the cluster's MPI module does that,
and a mismatch shows up as an undefined `mpi_types_...` symbol. Each
tarball's `FPM.md` has the details.

## Basis sets

terco ships a library of Basis Set Exchange basis sets (`basis_sets/`: the
Pople, def2, Dunning and fitting families). `cmake --install` puts them in
`share/terco/basis_sets`, where `terco -b def2-tzvpd` finds them by name;
with fpm, set `TERCO_BASIS_PATH` as above. Any other BSE JSON file can be
given by path. A basis with f functions needs a `TERCO_LMAX=3` build, and g
functions (QZ sets) need `TERCO_LMAX=4`; a lower build refuses them and says
which to rebuild with.

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

## fpm from a git checkout

The release tarballs above are the easy way. From a checkout:

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
