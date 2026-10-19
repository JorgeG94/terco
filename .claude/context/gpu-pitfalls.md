# GPU and nvfortran pitfalls

Each of these has produced silently wrong numbers or memory corruption in
terco. The rules in `.claude/skills/pr-review/SKILL.md` section 2 come from
the same history. This file adds the debugging side.

## Working method

- Keep two builds: a CPU build (gfortran 15+) to iterate and debug, and an
  NVHPC GPU build to confirm and to measure. Debug on the CPU first: it is
  faster and its diagnostics are better.
- When a kernel is wrong on the GPU but right on the CPU, and
  compute-sanitizer reports nothing, dump the kernel's arrays from both
  builds and diff them. An all-zero slice points at the cases below.

## Kernel shape

- **Multi-index `do concurrent` with a small outer extent.** A two-index
  `do concurrent (i=1:n1, j=1:n2)` whose outer extent is small and set at
  run time has skipped its last outer slice on the GPU (nvfortran 26.5,
  cc70). The CPU is fine and memcheck is clean. Flatten it to
  `do concurrent (t=1:n1*n2)` and decode `i`/`j` in the body.
- **No automatic arrays in `!$acc routine seq` procedures.** Use
  fixed-size, lane-private scratch, and atomics for shared accumulation.
- **Item routines are `pure`, carry `routine seq`, sit in the driver's
  module and take explicit-shape dummies.** Cross-module device calls and
  assumed-shape dummies are measured slowdowns, and cross-module calls do
  not inline.
- **No CUDA Fortran and no streams.** Fill the GPU from inside one kernel.

## Data lifetime

- **`trc_linalg_t` has no finalizer.** Always call `%release()` before one
  goes out of scope. An unreleased handle leaves its device workspace
  mapped after the host array is freed. A later allocation at the same host
  address is then "present" with the wrong size and the wrong device
  memory, and corrupts results far from the cause.
- **`mem:separate` needs every array mapped.** A missing `copyin` is an
  illegal address at run time, or worse, a silent implicit copy that works
  until the compiler changes. Map arrays explicitly, as the neighbouring
  pass does.
- **Derived types mapped whole.** Adding allocatable components to a type
  that is mapped with `enter data copyin(obj)` has triggered an nvfortran
  internal compiler error. Prefer flat arrays.
- **Local memory.** Large generated kernels need CUDA local memory. The DF
  path reserves it and gives it back after the guess
  (`trc_release_local_mem_reservation` in `src/core/trc_linalg.F90`). New
  big kernels should follow the same pattern.
- **Data layouts that kernels index.** Changing a layout such as the pair
  list has triggered an nvfortran miscompile before. Change one only with a
  GPU run that compares energies before and after.

## MPI

- A program or procedure that makes a collective MPI call must have no
  internal (`contains`) procedures. Under NVHPC 26.5 with HPC-X/OpenMPI,
  even an unused sibling deterministically corrupted later GPU work, and in
  one case made every `MPI_Allreduce` return 0 or NaN. This applies to test
  programs too: put helpers in a module.
- Make sure FindMPI picks the MPI that matches the Fortran compiler (NVHPC's
  bundled OpenMPI for nvfortran). A stray `MPI_ROOT` or similar in the
  environment can make it link a gfortran-built MPI.

## Build gotchas

- NVHPC's linker refuses to relink a shared library that already exists
  ("input file libterco.so is the same as output file"). Delete
  `build/libterco.so` before a relink. That includes any helper target that
  pulls the terco link rule in again.
- A CMake reconfigure regenerates the `CTestTestfile.cmake` files. Any
  environment you injected there (for example `CUDA_VISIBLE_DEVICES` to pin
  tests to one GPU) is lost and has to be redone.
