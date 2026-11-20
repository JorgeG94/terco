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

## Codegen

- **`pure elemental` automatic-length character functions, in a very large
  program unit.** nvfortran 26.5 and 26.9 (`-fast -O3`, host and GPU builds
  alike) miscompile pic's `to_lower` -- `pure elemental`, automatic-length
  (not allocatable) result -- at several call sites inside `trc_main`
  (`app/main.F90`): the result comes back the right length but blank
  (all-space), so e.g. `--func pbe` silently fell back to a plain SCF with
  no exchange-correlation term. Not a GPU/OpenACC issue: it reproduces
  identically on the minimal host build. A non-elemental function with an
  **allocatable** result, called the same way, is unaffected --
  `to_lower_alloc` in `trc_run` (`app/main.F90`) is that workaround, used
  at every CLI normalisation site in `trc_main`. Scale of the call site's
  own procedure seems to matter (an isolated synthetic reproducer of similar
  size did not trip it); if you hit something that looks like this
  elsewhere, suspect any `pure elemental` automatic-length character
  function called from a large procedure first.
- **ifx 2025.0 / 2026.0, -O2 and above: a `do concurrent` index passed by
  reference to a pure routine that is not inlined.** The compiler treats the
  call as loop-invariant, so every element gets the first iteration's value
  (found through the Intel CI job: RI-CCSD T1 wrong by 3e-2 relative, (T) by
  1e-2 Eh; -O1 and ifx 2022.1 are right). The same holds for a scalar
  computed in the body and passed on. A six-line reproducer is in the GitHub
  issue titled "submit to Intel". Cure in terco: pass the index through the
  `TRC_DC_INDEX` macro of `src/inc/trc_dc_index.inc` at the call site, which
  is `(i)` under ifx (an expression, not the variable) and the plain
  variable for every other compiler. Every such call in `src/` uses it, the
  generators (`scripts/gen_*.py`) emit it, and new code must too. An input
  that is already an expression (`t + 1`) needs nothing, an output argument
  cannot take it, and the few older item routines with a `value` index
  dummy are already safe. The macro exists only because of this ifx bug:
  when a fixed ifx is the CI compiler, remove it and its call sites and go
  back to the plain variable.
- **nvfortran 26.5 HOST build (`-fast -O3`; the device build is right):
  `integer, intent(in), value` dummies in item routines.** They gave wrong
  RI-CCSD amplitudes (bisected: 16 failing checks in `check_rccsd_run` with
  the `value` dummies, 0 without). There is no small reproducer yet. Do not
  use `value` dummies in item routines; `TRC_DC_INDEX` is the cure for the
  ifx defect above, not a `value` dummy.
- **nvfortran 26.5 (not 26.9), `-O3`, Volta: one generated Hermite item of an h build
  (`TERCO_LMAX=5`, `[10|8]` batch 3 of 18, `mdj_s_item_108_b3`).** Symptom:
  `CUDA_ERROR_ILLEGAL_ADDRESS` in that item as soon as a primitive quartet
  survives screening (compute-sanitizer: invalid global read of `Dt()`, 10 GB
  from any allocation). It is clean at `-O1` and `-O2`, clean under gfortran
  `-fcheck=all`, and reproduces on synthetic data with every launch argument
  in bounds, so it is not an index error. `-O3 -Mnovect`, `-Mnounroll`,
  `-Mnopre` and `-Mnolre` still fail. A later `-O2` on the command line does
  not override an earlier `-O3`; the level is set with the routine directive
  `!pgi$r opt=2`, which `scripts/gen_md_j.py` emits for the batch named in
  `O2_ITEMS`, inside `#ifdef TRC_NVHPC_O3_ITEM_BUG` (CMake defines it for
  nvfortran older than 26.9; 26.9 is clean at `-O3` and matches `-O2` to 2e-16;
  earlier versions are untested and treated as affected; the run-time cost of
  the `-O2` item is not measured) (keyed by the batch's content, and the generator fails if the
  key stops matching one batch). The directive also stops the item being
  inlined into its launcher, so the cure may be the de-inlining. A crash is
  the loud form of this defect; a silent one would be a wrong number, and a
  crash-only sweep cannot see it.
  **Mandatory GPU smoke test for every h build and every nvfortran upgrade**,
  run directly on a device (CI has no GPU and no `TERCO_LMAX=5` job):
  `check_h_vs_pyscf ints_hco` and `check_h_vs_pyscf rhf_water_min`, then the
  HCO SCF cases (`uhf_hco_min`, `b3lyp_hco_min`). Re-run the sweep that
  launches every `mdj_sb_<lp><lq>_<b>` once on random data after any change of
  `SPLIT_BATCHES` or `SPLIT_WEIGHT_EXP`.

## MPI

- A program or procedure that makes a collective MPI call must have no
  internal (`contains`) procedures. Under NVHPC 26.5 with HPC-X/OpenMPI,
  even an unused sibling deterministically corrupted later GPU work, and in
  one case made every `MPI_Allreduce` return 0 or NaN. This applies to test
  programs too: put helpers in a module.
- Never link static terco or pic into a shared library that is loaded next
  to an executable that also has them (for example a SHARED libterco_opt
  beside a static libterco). The `.so` gets its own copy of pic-mpi and of
  the legacy `use mpi` common blocks, and the executable calls the `.so`'s
  pic_mpi. With Cray MPICH, MPI_IN_PLACE lives in `mpifcmb5_`, and libmpifort
  compares against the executable's copy, so MPI_IN_PLACE is taken for a
  real buffer. Symptom: a multi-rank run hangs at the first Fock allreduce
  with the GPUs at 0% (a SIGSEGV in MPI's local copy, retried forever by
  nvfortran's handler). Check `ldd app/terco` and `nm` for pic_mpi symbols
  in more than one object. Other MPIs keep MPI_IN_PLACE in their shared
  Fortran library, so the duplicate goes unnoticed there.
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
