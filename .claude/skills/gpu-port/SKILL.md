---
name: gpu-port
description: Bring a method into terco so that it runs on the GPU at the speed the rest of the library does — a new integral, gradient term, post-HF step or analysis, whether it starts from a paper, from pyscf, libcint, libfint, metalquicha or from terco's own host code. Covers how to find the parallel item, one launch instead of a host loop of launches, who owns each output element, contracting on the device, device data and memory budgets, the validation ladder and the measurements a PR needs. Use before designing or writing any new compute path, and when reviewing one that is correct but slow.
---

# Porting a method to the GPU in terco

A port is of the **mathematics**, not of the source's control flow. Reference
codes are written for one CPU core walking shells, atoms or grid batches in
nested loops and filling matrices as they go. Copied as written, that
structure gives right numbers and a GPU that waits for the host. This skill
is the order of work that avoids that. The hard rules it relies on are in
`AGENTS.md` and `.claude/skills/pr-review/SKILL.md` section 2, whose "GPU
strategy" questions are what a reviewer will ask of the result; the compiler
traps are in `.claude/context/gpu-pitfalls.md`. Read those too.

## 0. Before any code

1. **Write the quantity down.** Its formula, the size of every index, and
   the size of the **result**. A gradient is `3 natm` numbers, an energy is
   one: nothing in the plan should be larger than the inputs it needs.
2. **Name the reference** the port will be validated against (pyscf,
   libcint, finite differences) and the tolerance you expect, before
   writing the kernel.
3. **Check the licence of what you read.** terco is MIT. Apache/BSD/MIT
   sources may be ported with attribution in the file header. Copyleft
   (GPL, LGPL, MPL) may not: implement from the paper.
4. **Find the nearest path terco already has** (the energy of the term you
   differentiate, the same integral class with another operator). Its task
   list, screen, caches and data mappings are usually the ones you need,
   and its run time is your yardstick (section 6).

## 1. Find the item, then make one launch

The **item** is the smallest unit of work whose result does not depend on
another item's: a shell quartet, a shell pair, a (shell pair, atom) triple,
a grid point, an (i, j) orbital pair.

- **One flat task list, built on the host.** Screening decides what is in
  it. The list holds indices only (pair, atom, class, offset).
- **One `do concurrent` over the list**, one flat index, decoded in the
  body. Not a two-index loop with a small run-time outer extent
  (gpu-pitfalls: it has dropped its last slice).
- **No host loop of launches.** A loop over atoms, shells, perturbations or
  grid batches that launches a kernel per iteration is the commonest way to
  lose the GPU: each pass pays its setup, mapping and transfer, and the
  device idles in between. If the source loops over something, that
  something becomes an **index of the task**. Translational invariance,
  "one atom at a time", "one perturbation at a time" are properties of the
  result, not reasons for a loop.
- **Chunk only for memory**, with a measured block size (section 3), and
  say so in a comment. One chunk is the normal case.
- **Sort or bin tasks of one kind together** (angular-momentum class,
  contraction length) when the work per item varies widely: a warp runs at
  the speed of its slowest lane.

## 2. Decide who owns each output element

Race conditions and most slowdowns come from here.

1. **Exclusive slot per task** (best): task `t` writes `out(:, t)` and
   nothing else. No atomics, the same arithmetic on host and device, and a
   cheap reduction afterwards. Use it whenever the per-task result is small.
2. **Fold on the device**: accumulate a block in registers or
   item-private scratch, then one `!$acc atomic update` per output element.
   This is how the Fock builds work. Results then differ between runs at
   round-off: gate tests on a tolerance, never on bytes.
3. **Never** a shared array written without atomics, and never scratch
   declared beside the loop (it becomes shared).

**Contract on the device.** If the result is a contraction with a density,
amplitudes or weights, do the contraction inside the item or in a device
reduction. Do not bring an `(nao, nao, k)` intermediate to the host to dot
it with a density there. A linear post-processing step (a spherical fold, a
derivative-shell fold, a symmetrisation) can almost always be moved onto
the *other* factor of the contraction through a small index or coefficient
array built once on the host.

**Dense algebra goes through `trc_linalg_t`** (GEMM, TRSM, batched GEMM on
the device), in as few large calls as memory allows, not as a host loop of
small BLAS calls. Always `%release()` the handle.

## 3. Data and memory

- **Map once, outside every loop.** `enter data` for everything the kernels
  read before the first launch, `exit data` on every path out, error
  returns included. An `update self` or a re-map inside a loop is a design
  smell: find what could stay on the device.
- **Everything a kernel reads is mapped explicitly.** Under
  `-gpu=mem:separate` a missing mapping is a run-time illegal address that
  the host build cannot show. Derived-type components need the parent
  mapped and the component attached.
- **Precompute what does not depend on the item.** Per-shell, per-atom or
  per-channel tables go into mapped caches built by their own small kernel,
  so the main item holds a handful of scalars. Large per-thread arrays are
  stack overruns and local-memory traffic (the header of
  `src/eri/trc_ecp.F90` is a worked example).
- **Budget device memory before allocating.** Ask what is free
  (`src/core/trc_device_memory.F90`), keep a reserve, size blocks from it,
  and stop with a message that states the sizes rather than at a failed
  allocation (`choose_grad_blocks` in `src/scf/gradient/trc_scf_grad_df.F90`).
  Allocate a large temporary where it is first needed, not at the top.
- **Give the forced-chunk path a switch** (an environment variable read by
  the block chooser) so a small test can run it.

## 4. The item routine

`pure`, `!$acc routine seq`, in the **same module** as its launcher,
explicit-shape dummies, fixed-size scratch declared inside it, no automatic
arrays, no `value` dummies, loop indices passed through `TRC_DC_INDEX`
(`src/inc/trc_dc_index.inc`). It calls only intrinsics and other item
routines of its module. The reasons, each a measured failure, are in
`.claude/context/gpu-pitfalls.md`.

Kernels are standard `do concurrent`: no CUDA, no OpenACC compute
constructs, no `async`. A class of integrals that needs many near-identical
routines is **generated** (`scripts/gen_*.py`, `docs/GENERATORS.md`), never
written by hand into `src/generated/`.

## 5. Validate in this order

Each rung catches what the next cannot.

1. **Host build first** (gfortran): it is the same code, compiles in
   seconds and has bounds checking. Debug there.
2. **The raw quantity against an independent reference**: integrals
   against libcint/pyscf, derivatives against finite differences of your
   own integrals (Richardson, so truncation is not mistaken for a bug).
3. **The contracted result against pyscf** at terco's own converged
   orbitals, as a committed anchor (`test/pyscf_refs/README.md`).
4. **Screened against unscreened**, and **chunked against whole**, on a
   case where the screen and the chunking really bite: several atoms far
   apart, a forced small block.
5. **Host against device** on the same input. They must agree to
   round-off; code that is right only on one of them is a bug.
6. **A case that can fail**: confirm the new test fails with the term
   removed or a sign flipped. Test the refusals you added.

Then the gates in `AGENTS.md` (fortitude, archcheck, reference freshness,
the OpenMP conversion if you touched `do concurrent` or directives).

## 6. Measure before the PR

Correct is half of the claim. A heavy path is not finished until it runs
near what the device can do.

- **Every stage has a row** in a step table (`src/core/trc_step_table.F90`)
  or timer. Work done outside the table is invisible: compare the table's
  total with the wall time of the whole run, and account for the gap.
- **Compare with the yardstick from section 0.** A first derivative should
  cost a small multiple of its energy term. A result that scales worse than
  its inputs, or a GPU that sits at low utilisation while the run
  continues, means the structure is wrong, not that the kernel needs
  tuning.
- **Profile, then change one thing.** `nsys` for where the time is, `ncu`
  for why a kernel is slow. Read register spills from the linked
  `libterco.so`, never the `.o`. Do not cap registers.
- **State the numbers**: system, basis, thresholds, GPU, compiler, the
  baseline, and unchanged energies or gradients. Add or re-record a
  `bench/regression` case when the path is performance-sensitive.
- **State the compile-time cost** of new generated kernels.

## 7. Shapes that look right and are wrong

| What the code does | What it should do |
|---|---|
| Loops over atoms (or shells, perturbations, batches) on the host and launches a kernel in each pass | One task list with that index in it, one launch |
| Builds an `(n, n, k)` array on the device, downloads it, folds or contracts it in host loops | Per-task contraction on the device, a small slot per task |
| Rebuilds caches, screens and mappings in every pass | Build once, keep mapped for the whole routine |
| Reuses an energy routine `N` times to get `N` derivatives | Reuse its **items** inside one derivative launch |
| Allocates the largest temporary at entry "so the budget sees it" | Allocate at first use, after what it replaces is freed |
| Two-index `do concurrent` with a small run-time extent | One flattened index |
| Assumed-shape dummy read inside a kernel | Explicit shape, sized from arguments |
| Passes one array as both an `intent(in)` and the output dummy | An in-place variant; aliased dummies are not conforming |
| Times only the kernel | A row per stage, and the total against wall time |
| "Converges", "runs on the GPU" offered as validation | The ladder in section 5 |

## Checklist for the PR description

- The item, the task list and what screens it; how many launches a run makes.
- Who owns each output element (slot or atomic) and where the contraction happens.
- Device memory: what is resident, the peak, what happens when it does not fit.
- The validation rungs passed, with the measured differences.
- The timing table, the yardstick it was compared with, and what is not optimised yet.
