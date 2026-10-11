# Plan: a Fortran-native, device-resident DFT XC integrator

Written 2026-09-03. Self-contained: assumes the reader has no prior context.

## Getting started: where everything is

**Do not commit this file.** It is a scratch plan; delete it when the work lands.

### Toolchains — the system compilers are NOT sufficient
`/usr/bin/gfortran` is **13.3**, and terco's host build needs **15+**
(`do concurrent` locality specifiers are Fortran 2018). There is no nvfortran on
`PATH`. Both live outside the system tree:

| Purpose | Path |
|---|---|
| Host build / CPU validation | `/shared/compilers/gcc/16.1.0/bin/gfortran` (15.1.0 also present) |
| GPU build (`-stdpar=gpu`) | `/home/jorge/install/nvhpc/26.5/Linux_x86_64/26.5/compilers/bin/nvfortran` |

Existing terco build dirs, already configured, showing exactly this:
- `~/dev/terco/build` -> gcc 16.1.0, `TERCO_ENABLE_GPU=OFF`
- `~/dev/terco/build-cpu` -> gcc 15.1.0, `TERCO_ENABLE_GPU=OFF`
- `~/dev/terco/build-gpu-l3` -> nvhpc 26.5 nvfortran, `TERCO_ENABLE_GPU=ON`

`TERCO_ENABLE_GPU=ON` with anything but NVHPC is **refused**, deliberately, rather
than silently building a slow CPU library. Other options: `TERCO_ENABLE_TESTING`,
`TERCO_REGENERATE` (re-runs the kernel generators), `TERCO_BUILD_SHARED`.

### Python environment
The venv is at **`/home/jorge/dev/metalquicha/.venv`**, and the way to get it plus
everything else on `PATH` is:

```bash
source /home/jorge/dev/mqc_worktrees/mqc_env.sh
```

**Always source it before any tooling.** Without it, `python3`, `fortitude`,
`pre-commit` and friends resolve to broken copies under `~/.local`. It carries
PySCF 2.14, which is the reference for anything needing an independent number.

### The repositories

| What | Path |
|---|---|
| terco (the work lands here) | `~/dev/terco` |
| metalquicha main checkout | `~/dev/metalquicha` |
| metalquicha worktrees | `~/dev/mqc_worktrees/` (agent workspaces; do not delete) |
| GauXC (design reference) | `~/dev/GauXC` |
| ExchCXX (translate these) | `~/dev/ExchCXX` |
| IntegratorXX (ignore) | `~/dev/IntegratorXX` |
| libxc (the oracle) | reached through metalquicha's build |

### Building metalquicha, if you need it for reference numbers
```bash
source /home/jorge/dev/mqc_worktrees/mqc_env.sh
cmake --preset default && cmake --build --preset default
ctest --preset mqc          # never pass -j; the tests are OpenMP-parallel already
```
Cap threads for test runs (`OMP_NUM_THREADS=4`); an uncapped run takes every core.

## The claim

GauXC (C++/CUDA) is the only library that runs the whole exchange-correlation
pipeline on the device. Everything else evaluates the functional on the host via
libxc. **Doing the whole pipeline in standard Fortran with `do concurrent`, and
no libxc on the hot path, is a claim nobody else can make.**

## Decisions already taken (do not re-litigate)

1. **`do concurrent` stays. Do NOT convert to CUDA Fortran.** The opening idea
   was to abandon `do concurrent` for CUDA Fortran. That was reversed on
   analysis: the XC quadrature is a *flat, regular* loop over grid points with
   uniform work per point. It has none of the shell-quartet irregularity that
   forced terco's 81 generated permutation-class kernels. This is the case
   `do concurrent` compiles well, and it is a better target for the thesis than
   the ERIs were.
2. **libxc cannot be the device functional, under either model.** You cannot
   call host C from an offloaded kernel. This is why GauXC ships ExchCXX. The
   libxc boundary is identical under `do concurrent` and CUDA Fortran, so it is
   not an argument for either.
3. **Do not port IntegratorXX.** mqc already has the equivalent, in pure Fortran.
4. **Do not port GauXC.** Take its *design*, not its code (see Scale, below).
5. **pic-device handles memory/streams.** Already vendored at `terco/pic-device`,
   with CUDA *and* HIP backends: `backend_malloc/free`, `memcpy` sync+async,
   pinned host alloc, events, elapsed time.

## Why this is worth doing at all (the number that justifies it)

Measured on metalquicha: **the XC quadrature is 89% of an LDA run and 98% of a
meta-GGA one. The Fock build is under 2%.** Accelerating integrals does almost
nothing for DFT. All the value is in the quadrature — the part terco does not
have yet.

## What exists, and where

### terco — `~/dev/terco`
ERI/Fock only. **No grid, no XC, no DFT anything.**
- 118 Fortran files. `do concurrent` in 86: **5 hand-written, 81 generated**.
- Generator: `scripts/gen_perclass.py` (also gen_vrr/gen_hrr/gen_1e/gen_df).
- **Precedent for mechanical source translation already exists**:
  `tools/dc_to_omp.py`, `tools/acc_to_omp.py`. The ExchCXX translator below is
  the same shape of tool.
- Premise (README): "can we evaluate ERIs on the GPU just using do concurrent".
  DFT extends this thesis rather than abandoning it.

### metalquicha — `~/dev/metalquicha`, worktrees in `~/dev/mqc_worktrees`
Grid machinery is done, pure Fortran, and portable. In `src/methods/`:
- `mqc_dft_grid.f90` (301 lines), `mqc_dft_partition.f90` (503),
  `mqc_lebedev.f90` (235), plus `mqc_dft_radial`, `mqc_dft_prune`.
- Dependencies are only `pic_types`, `mqc_error`, `mqc_physical_constants`,
  and each other. No JSON, no MPI, no C. **Lift as-is.**

XC evaluation: `backends/cenzontle/mqc_czt_xc.F90` (was `backends/libcint/` —
renamed on `chore/cenzontle-rename`; modules are `mqc_czt_*`).
- The single libxc file in the whole tree (~50 `xc_f03` references).
- **The loop is already blocked over points**: `do g0 = 1, npts, ctx%point_block`
  around line 459, and libxc is called per block with flat arrays in and out
  (`rho_blk` -> `exc_i`, `vrho_i`). The structure you want is already there.
- **No non-libxc path exists**: without `MQC_WITH_LIBXC` it hard-errors.

### The three reference libraries (all present locally)
| Repo | Path | License | Verdict |
|---|---|---|---|
| GauXC | `~/dev/GauXC` | BSD | Copy the *design* only |
| ExchCXX | `~/dev/ExchCXX` | **MIT** | Translate the kernels |
| IntegratorXX | `~/dev/IntegratorXX` | BSD | Ignore — mqc has it |

## Scale, measured — why "port GauXC" is the wrong instinct

Inside `GauXC/src/xc_integrator/`:
- **201,375 lines are generated kernels** (79 collocation files per angular
  momentum; Obara-Saika files up to 66,537 lines in one file).
- **30,052 lines are architecture** — and much of even that is the sn-LinK
  *exact exchange* path, a separate feature from the XC quadrature.
- The part actually worth reading is small:
  `src/xc_integrator/integrator_util/` (60K) and `shell_batched/` (80K).

## ExchCXX: the transcription job

`~/dev/ExchCXX/include/exchcxx/impl/builtin/kernels/` — **104 built-in
device-callable functionals**. Target four to six, not all of them.

Sizes (generated lines, translated by script, not written by hand):
`slater_exchange.hpp` 557, `vwn3.hpp` 2,897, `pbe_x.hpp` 924, `pbe_c.hpp` 2,914.
LDA + PBE is roughly 7,300 lines. B3LYP additionally needs B88 and LYP.

**The bodies are straight-line scalar arithmetic** — Maple-generated
`const double tN = ...;` sequences. No classes, no templates, no pointers.
Only three constructs need handling:

| C++ | Fortran |
|---|---|
| `const double tN = expr;` | `real(dp) :: tN` + assignment |
| `safe_math::cbrt(x)` etc. | an `elemental pure` function |
| `piecewise_functor_N(c,a,b)` | `merge(a, b, c)` — **see risk below** |
| `constexpr` constants | `real(dp), parameter` |

Result: a `pure` scalar-in/scalar-out function, callable inside `do concurrent`
with no host round-trip.

### RISK, check this first
`merge` evaluates **both** arms. If any `piecewise_functor_*` exists to dodge a
NaN or a division by zero as rho -> 0, a naive `merge` translation reintroduces
exactly the failure the guard prevented. Confirm the semantics in
`include/exchcxx/impl/builtin/util.hpp` (or wherever `piecewise_functor_*` is
defined) **before** writing the translator. Symptom if missed: mysterious
grid-point blow-ups much later, on some functionals only.

## Sequence

**Phase 1 — translator + functionals.** Write `tools/exchcxx_to_fortran.py`
(same shape as the existing `dc_to_omp.py`). Emit `pure` Fortran for Slater,
VWN, PBE_x, PBE_c. Resolve the `piecewise_functor` risk first.

**Phase 2 — validate on the host, before any GPU.** `do concurrent` is standard
Fortran and runs on the CPU under gfortran, so the entire pipeline is
developable and testable with no device at all. Check each translated functional
against libxc pointwise via mqc's existing `mqc_czt_xc.F90` path. **This is the
big de-risking step: get numbers right on CPU, then offload.**

**Phase 3 — grid: move it into terco.** The five grid modules currently live in
metalquicha (`~/dev/metalquicha/src/methods/`) and **need to move to
`~/dev/terco`**. They are pure Fortran and depend only on `pic_types`,
`mqc_error`, `mqc_physical_constants` and each other, so this is a copy plus a
rename of the module prefix (`mqc_` -> `trc_`), not a port.

Decide early whether metalquicha then *depends* on terco for its grid, or whether
the two carry copies for a while. Copies diverge; a dependency means mqc gains a
terco build requirement. Worth asking Jorge rather than choosing silently.

**Phase 4 — batching and screening.** The one genuinely new architectural piece,
and the thing worth taking from GauXC: group grid points into tasks, screen each
task down to the shells with non-negligible support on it, and collocate only
that local basis. Without this the device just does more arithmetic faster.

**Phase 5 — collocation kernels.** chi, grad-chi, and (for meta-GGA) tau on
grid points, generated per angular momentum. terco's `scripts/` already
establishes the generator pattern.

**Phase 6 — rho and V_xc.** Both GEMM-shaped; `do concurrent` over points for
the pointwise parts, BLAS for the contractions.

**Phase 7 — offload and wire back.** pic-device for buffers and streams. Expose
as an alternative XC backend behind mqc's existing `xc_context_t` interface so
libxc remains selectable as the reference oracle.

## Verification notes

- Local box: sm_70 cards, no libcuest. cuEST needs sm_80. **But this plan does
  not need cuEST at all** — it needs nvfortran (or gfortran for CPU validation).
- Real GPU verification is Perlmutter A100s. Build there with tblite OFF —
  toml-f, reached through tblite, has a backslash string literal nvfortran
  rejects. `cmake --preset perlmutter` encodes the flags.
- libxc stays the reference forever. Never delete that path.

## Attribution

ExchCXX is MIT and GauXC is BSD; both permit translation with attribution.
Credit ExchCXX in the header of every generated functional file, and GauXC in
whatever documents the batching design.

---

## Status, 2026-09-03 (end of the first session)

Phases 1 through 6 are done on the host and phase 7 runs on the local
V100s. Nothing is committed; see the list at the end.

### What landed

- `tools/exchcxx_to_fortran.py` translates any LDA or GGA kernel. Seven are
  translated into `src/xc/`: Slater, VWN5, VWN_RPA, PBE x/c, B88, LYP.
  `test/check_xc.F90` + `test/xc_ref.py` (recorded `reference_xc.dat`)
  compare every derivative through second order, both spin cases, against
  libxc 7.0.0 through pyscf: energies and potentials agree to 1.2e-12,
  second derivatives to 3.7e-12.
- The grid modules moved to `src/grid/` with a `trc_` prefix and a minimal
  `trc_error`; `test/check_grid.F90` carries the identity checks.
- `src/xc/trc_collocation.F90`, `trc_xc_batch.F90`, `trc_xc_functional.F90`,
  `trc_xc.F90`: values and gradients on points, GauXC-style batching and
  screening, libxc's compositions, and the two-kernel integrator.
  `check_ao` (values vs pyscf eval_ao to 1e-16; numerical overlap vs
  `trc_1e`), `check_xc_energy` (E_xc and V_xc vs pyscf on the same grid
  and density to 1e-14), `scf_rks` + `rks_ref.py` (a full RKS on water in
  6-31G uncontracted to primitives, SVWN/PBE/BLYP/B3LYP, vs pyscf RKS on
  the same grid: 3e-12 on the host, 1.5e-12 on the GPU).
- GPU: `build-xc-gpu` (nvfortran, cc70, LMAX=2). Every XC check and the
  RKS pass on device 0 with the host's numbers; NV_ACC_NOTIFY shows the
  point kernel launched over the whole grid, one thread per point.
- `src/xc/bit_repro.f90` + `bit_repro_helpers.f90`: Jorge's bit-reproducible
  transcendentals, vendored from `~/nci/learning_tools/ci_enabled/bitwise_adventures`
  with an `!$acc routine seq` beside each `!$omp declare target`. `xc_cbrt`
  is now `cuberoot` and `xc_erfcx` is `exp_reprod(x*x)*erfc_reprod(x)`.
  The translated kernels still call the intrinsic exp/log/sqrt/atan;
  switching exp and log to the `_reprod` versions is a one-line change in
  the translator's FUNC_MAP, not done pending Jorge's call. NOTE the flag
  contract: bit_repro wants no FMA and no FTZ, and terco's nvfortran build
  uses `-fast` (FMA on) -- CPU/GPU bit identity needs `-Mnofma -gpu=nofma`
  and the `-Mnoflushz -Mnodaz` pair before it can be claimed.
- `test/trc_test_basis.F90` gained `build_631g(..., uncontracted=.true.)`;
  `test/mult_ref.py` had a cwd-relative path that failed under ctest and
  now resolves the file beside itself (pre-existing, fixed in passing).

### Findings worth keeping

- `piecewise_functor_*` takes VALUES, so C++ evaluates both arms exactly as
  `merge` does. The plan's risk did not exist.
- libxc is the noisy one at low density: its polarised GGA potentials drift
  from the closed form below a spin density of ~1e-11 (1e-10 relative at
  5e-13), and its polarised fxc carries ~1e-9 cancellation residue at a
  total density of 1e-10. The translated kernels match the closed form to
  1e-14 there. The reference sample therefore starts at 1e-9.
- NWChem pruning starves d-on-d cross-centre overlaps on the light probe
  system: 6.4e-4 at level 3 pruned, from pyscf's own grid too. Unpruned
  level 5 gives 9e-8.
- The three-centre probe (Be/Li/He or Be/He/He in the toy spd basis) is
  fine for a fixed density and unusable for an LDA SCF: no core function
  on Be, HOMO-LUMO gap ~1e-3, and pyscf's DIIS fails on it as well. The
  SCF check uses water for that reason.
- nvfortran: `erfc_scaled` has no device version (spell it out);
  assumed-size dummies cannot appear inside a device `do concurrent`; and
  a loop body given inline gets its INNER loops parallelised across the
  block -- keep the body in a `!$acc routine seq` procedure, as terco does.
- ExchCXX's B3LYP uses VWN5; libxc's (and pyscf's `B3LYP`) uses VWN_RPA.
  Both are here, as `b3lyp` (libxc) and `b3lyp5`.

### The SCF driver (later the same day)

- `src/trc_scf.F90` (module `trc_scf_driver`): HF and Kohn-Sham, R and U,
  DIIS, every matrix resident on the device; `src/trc_linalg.F90` is
  cuSOLVER + cuBLAS under -acc (through NVHPC's `cublas`/`cusolverdn`
  modules and `host_data use_device`, `-cudalib=cublas,cusolver`) and
  LAPACK/BLAS on the host. The library now links LAPACK in every build.
- C entry `trc_scf` in the skin and `include/trc_c_interfaces.f90`;
  exercised by `test/link_gfortran.F90` alongside its own loop.
- `test/scf_rks.F90` + `rks_ref.py`: seven cases (water RHF/SVWN/PBE/BLYP/
  B3LYP, OH radical UHF/UPBE), all 1e-12 against pyscf, with the control
  that pyscf's energy AT terco's density matches before it iterates.
- FINDINGS: (1) a binding label may not equal a module name -- the
  standard says so, nvfortran tolerates it, gfortran folds every call
  into the module onto the label and the C entry recurses to death. The
  pre-existing `trc_fock` entry had this; modules are now `trc_eri` and
  `trc_scf_driver`, labels unchanged. (2) pyscf groups an atom's shells
  by angular momentum, so a density crossing to pyscf in a basis with
  interleaved s and p shells must be permuted through (atom, l, exponent);
  every earlier check used s-p-d-per-centre bases where this is invisible.
- The mqc side (`run_terco_scf` beside the cuEST bridge) is to be written
  fresh off origin/main; the feat/terco-device-scf worktree predates the
  cenzontle rename and is not a base.
- MPI (ad78017): the Fock build is split over ranks by a stride through
  the sorted work list (generator emits it), pic-mpi is the dependency
  (serial backend when TERCO_ENABLE_MPI=OFF), reduced G is broadcast from
  rank 0 and D every iteration (allreduce bits differ per rank, XC atomics
  too). scf_mpi runs 1/2/4 ranks under ctest (launcher from the compiler's
  dir). Cholesterol RHF 5.0 s -> 2.2 s on 4 V100s; PBE 22.7 -> 20.1 s
  because the XC integration is still whole on every rank -- NEXT: split
  batches over ranks in xc_drive and allreduce V/E.
- SCF driver: GWH guess default, damp 0.3 + level shift 1.0 until the
  commutator < 0.5, then DIIS; core guess never converged cholesterol.

### Later still (same day): setup on the device, install, ABI

- Grid partition and the batching maxima run on the device (4b2670e):
  cholesterol setup 22 s -> 1.6 s (partition itself 0.2 s), same weights.
  Remaining host setup: Lebedev/radial assembly and the one-electron
  integrals, all small.
- TERCO_BUILD_SHARED defaults ON; `make install` ships only terco
  (include/terco + lib); -cudalib is a PUBLIC link option on the target,
  not a global flag (pic's and pic-blas's executables broke on it).
- pic-blas is a dependency now (host linear algebra, explicit interfaces;
  CI's -Werror=implicit-interface). DIIS is terco's own solve with the
  error block scaled, as mqc_diis does.
- trc_scf has a `verbose` argument (prints the iteration table).
- The mqc bridge lives on metalquicha branch feat/terco-gpu-scf (gs2):
  run_terco_scf beside the cuEST bridge, MQC_ENABLE_TERCO + TERCO_ROOT,
  general contractions split into segmented rows, validated against
  cenzontle to 1e-12 (HF) / 5e-11 (PBE) / 3e-13 (cc-pVDZ HF).
- Still pre-existing on main, not this branch: the generated DF kernels
  differ from gen_df.py's output, and ifx rejects trc_boys.F90's
  51840-element initialiser (token limit).

### Open, for Jorge

- Does metalquicha now depend on terco for its grid, or do the two carry
  copies? Copies exist as of today.
- Spin-polarised integrator (UKS) and meta-GGA: kernels and wrapper exist
  for polar; the translator refuses mGGA on purpose.
- Performance, first numbers (`test/bench_xc`, PBE, one evaluation of
  E_xc + V_xc, fixed density, level-3 pruned grid, 6-31G):

  | system | points | AOs | host gfortran 16 | V100 | split on the V100 |
  |---|---|---|---|---|---|
  | water | 20k | 13 | 64 ms | 1.0 ms | points 0.6, pairs 0.4 |
  | benzene | 144k | 66 | 3.9 s | 17.8 ms | points 12.4, pairs 5.9 |

  The point kernel (collocation + density + functional) is two thirds of
  the device time; the pair kernel with its one atomic per (batch, u, v)
  is not the bottleneck. Nothing has been tuned. Larger systems and the
  chunk boundary (benzene is one chunk at the default budget) are the
  next measurements.
- The test programs' pseudo-random densities were a floating-point LCG
  that gfortran and nvfortran evaluated to different bits (FMA), which
  looked like a GPU bug on contracted shells for half an hour. They are
  integer LCGs now, and host and GPU agree to every digit.
- Nothing is committed. New: `src/xc/*`, `src/grid/*`, `tools/exchcxx_to_fortran.py`,
  `test/{check_xc,check_grid,check_ao,check_xc_energy,scf_rks}.F90`,
  `test/{xc_ref,ao_ref,xc_energy_ref,rks_ref}.py`, `test/reference_xc.dat`.
  Modified: `src/CMakeLists.txt`, `test/CMakeLists.txt`. Do not commit
  this file, nor `pic-device/`.

## Status 2026-09-03 (later): XC batches over ranks (9faf654)

- Batches cost-sorted (npts*nloc*(4+nloc)) and dealt round-robin; other ranks'
  batches inactive via the screening flag. V/E_xc/N summed then bcast from 0.
- Density screen moved to the device (one thread per batch); it was 2.3 s of
  host time per cholesterol SCF and read a stale host density.
- Cholesterol PBE: 22.6 s (1 V100) -> 7.3 s (4 V100); RHF 4.95 -> 2.14 s.
  Breakdown on 4: setup 0.4, grid 1.7, fock 1.2, xc 3.7 (points 2.5, pairs 1.1).
  Grid build (1.7 s) is now the largest serial piece; the batching is per
  rank and could be split the same way.
- ctest in build-mpi-gpu needs NVHPC's OpenMPI first on LD_LIBRARY_PATH
  (comm_libs/13.2/hpcx/hpcx-2.50/ompi5/lib): the shell's openmpi-5.0.5 is
  gfortran-built and lacks nvfortran's mpi_f08 type descriptors.

## Status 2026-09-03 (evening): partition over ranks (c14dbbe), fpm (d7ecf5c)

- Becke partition split by contiguous point blocks per rank, summed + bcast;
  build_dft_grid_block is a separate entry (nvfortran 26.5 ICE on callers of
  the widened interface). Grid 1.69 -> 0.84 s; cholesterol PBE 6.4 s on 4 V100.
- The nvhpc CI SIGILL on 9faf654 did not reproduce on 26.5 or 26.1 locally
  and the next CI run was green: runner flake. Locality lint fixed.
- fpm job still red: pic-mpi needs an MPI under fpm (its manifest compiles all
  backends). Fix belongs in pic-mpi (guard the MPI backends, .F90 + USE_SERIAL)
  or the fpm job needs an MPI built for gfortran 15.

## Notes: Stocks & Barca, JCTC 2025, 21, 10263 (GPU XC evaluation)

**What they compared.** Four algorithms for rho on the grid and V_xc from it,
all in EXESS, single A100, B3LYP/def2-SVP, 55 molecules (glycine chains,
BN sheets, diamond, water clusters):

1. *Dense*: full N x N_grid GEMMs, no sparsity. Cubic, only for tiny systems.
2. *Batch Dense D*: per batch of grid points, slice D to the batch's
   significant shells (nloc x nloc), X = chi.D_loc (GEMM), rho = row dots,
   Z = diag scaling of chi by (w eps_rho/2 + 2 w eps_gamma grad rho . grad chi),
   V_loc = chi^T Z (GEMM), symmetrize. This is GauXC's algorithm; the whole
   difference is batching and launch grouping.
3. *Batch Dense C*: same but X = chi.C_occ (nloc x Nocc), rho = sum X^2.
   Wins when Nocc << N (large basis, dense systems: diamond, def2-TZVP).
4. *Direct*: one thread per point, loop over significant shell PAIRS with
   the Gaussian product cutoff, accumulate rho directly, no intermediate.
   Register bound, uncoalesced D loads; 5-10x slower than 2 at scale.

**Results.** Batch Dense D best for large sparse systems, C for small dense;
1.4-5.2x faster than GauXC and GPU4PySCF. 50-70% of A100 fp64 peak at large
N. Grouped GEMM is 2/3 of the time; the rest (collocation, diag scaling,
dot products) is memory bound and at the roofline.

**Design points that matter for us:**

- *Batching*: octree gives uneven batches; they run a Hilbert space-filling
  curve over the octree leaves and cut equal batches of exactly the target
  size. 10-25% more shell pairs per batch, but 20% faster overall for small
  systems. Optimal target batch 2048 points (GauXC 512, GPU4PySCF 4096).
- *Kernel grouping*: ~100 batches per launch; 4x over one batch per launch.
  cublasDgemmGroupedBatched cost 1 ms/call (cudaGetDeviceProperties inside);
  they used CUTLASS grouped GEMM instead.
- *Screening*: shell cutoff radius by binary search on the primitive bound
  c r^l exp(-z r^2) <= tau, tau = 1e-10; error is linear in tau down to
  1e-12, below that nothing. Bounding polyhedra (box + octahedron) per
  batch for the shell-batch test. Grid points with partition weight < 2e-15
  dropped before batching (15-20% of the grid).
- *Collocation*: one thread per (shell, point), templated on l and on
  contraction degree; exponent evaluated only if z|dr|^2 < -ln(tau)*1e-3;
  value and gradient in one kernel; chi evaluated ONCE per batch and reused
  for density, potential, and integration.
- *Two CUDA streams* alternating per batch so the copies overlap the GEMMs;
  pinned-memory pool of 100 buffers for the group descriptors.
- *Multi-GPU*: static split of the grid; scaling 2.7-3.1x on 4 GPUs because
  GauXC splits by atom centre and the batches lose locality. Splitting
  after forming the full grid is what they recommend (= what terco does).
- *Not yet done by them*: mixed precision (TF32 for the small-|D| block of
  the slice, with functions sorted by their batch bound); fusing the four
  diagonal-scaling kernels (5-10%); linear-scaling shell-batch assignment.

**Where terco stands.** terco's point kernel is the paper's *Direct*
algorithm (one thread per point, D slice read per point, then a pair kernel
with atomics for V). Cholesterol 6-31G (N=344, ~846k points, level 3) costs
0.60 s per iteration of XC on one V100. The paper's Batch Dense D at N~350
on a much finer grid (30k pruned points/atom vs our ~11k) is ~0.02 s on an
A100, and even their Direct is ~0.07 s. Allowing 2x for V100 vs A100 we are
10-15x off, per point worse than that. The screening work is done (batch
bounds, density screen, rank split); the kernel structure is the gap.

**Plan: Batch Dense D in terco.**

1. Per chunk of batches, one grouped GEMM pass rather than per-point loops:
   - collocation kernel writes chi(npts_b, nloc_b) and gchi per batch (we
     have shell_collocate; the tile layout is what changes),
   - slice kernel gathers D_loc(nloc_b, nloc_b) from D via b_ao,
   - X = chi . D_loc: batched GEMM (cublasDgemmBatched with pointer arrays
     built once per grid, since the batch structure is fixed across the SCF;
     or sort batches by nloc and use strided-batched over equal-nloc groups),
   - rho, grad rho, eps and derivatives: one thread per point (existing
     xc_eval_point), writes Z(npts_b, nloc_b) = diag scaling of chi and gchi,
   - V_loc = chi^T Z: batched GEMM, then scatter-add V_loc into V. The
     scatter is the only atomic, one per (batch, u, v) as now, or avoided by
     accumulating per batch and reducing by sorted (u,v) keys.
2. Batch size: raise grid_max_pts to 2048 and cut equal batches; the
   bisection we have gives ~half the target on average. A Hilbert order over
   the leaves is cheap to add on the host.
3. Drop points with partition weight < 2e-15 before batching.
4. Group launches: chunks already hold ~100+ batches, so one batched GEMM
   per chunk is the grouping.
5. Everything is cuBLAS from OpenACC host_data, no CUDA; the host build gets
   pic-blas dgemm in a loop, same numbers.
6. Batch Dense C later for large basis sets; needs C_occ from the driver,
   which trc_scf has.

Expected: XC per iteration from 0.6 s to <0.1 s on cholesterol, more on
larger systems, and the multi-GPU split carries over unchanged.

## Status 2026-09-03 (night): Batch Dense D landed (a5c825f, + weight drop)

- xc_drive is now collocate tiles -> gather D_b -> strided-batched GEMM X ->
  functional per point into Z -> strided-batched GEMM V_b -> atomic scatter.
  Tiles padded to (128 pts, 32 fns), batches stored per chunk in shape order
  so each run of equal shape is one cublasDgemmStridedBatched.
- Cholesterol PBE, 24 iters, 1 V100: XC 14.3 s -> 2.2 s (kernels 1.2 s:
  collocation 0.6, functional 0.36, screen 0.1; GEMM groups ~1.0 s).
  SCF 20.5 -> 8.3 s; on 4 V100s 6.4 -> 3.5 s. Points below 2e-15 dropped.
- Batch target 2048 was SLOWER (9.2 s) with the bisection batcher: bigger
  boxes reach more functions. Equal-size Hilbert batches are what makes
  2048 pay in the paper; not done.
- nvfortran 26.5 bugs met: ICE on callers of an interface widened by two
  optional dummies (build_dft_grid -> separate block entry), and a front-end
  segfault on a unit with eight calls to trc_xc_rks/uks (test split into
  two contained procedures).
- Next: collocation kernel as one thread per (shell, point) (paper's Alg 2),
  Hilbert equal batches, Batch Dense C for big basis sets, fewer chunks
  (bigger budget) to cut launches.
- Algorithm 2 collocation (one thread per shell,point): 0.61 -> 0.45 s;
  XC 2.10 s per cholesterol SCF on 1 V100. The kernel writes ~4 GB of tiles
  per iteration; ideal ~0.1 s, so 4x off, probably the exp evaluations and
  index arithmetic; the paper's exponent-skip test is the next cheap thing.

## Plans: GPU RI-MP2 and range-separated hybrids (2026-09-03)

**RI-MP2** (recommended first, ~2-3 days). Have: trc_df_2c (P|Q), trc_df_3c
(mu nu|P) as (nao,nao,naux), gemm_strided, cuSOLVER handles, C and eps from
trc_scf. Steps: (1) cuSOLVER potrf of (P|Q) in trc_linalg; B = (mu nu|P) L^-T
by trsm. (2) (ia|Q) = C_occ^T (mu nu|Q) C_virt, strided-batched over Q, chunk
over aux shells for memory. (3) per i: (ia|jb) for all j as one strided
batched GEMM, then a kernel over (j,a,b) with (2(ia|jb)-(ib|ja))(ia|jb)/de;
SCS/SOS at that line. (4) ranks: split i, allreduce the energy. (5) PySCF
dfmp2 with the same BSE aux basis.

**Range-separated** (~1 week). (1) erf-attenuated ERI: s = w/sqrt(rho+w^2);
[00|00]^(m) *= s^(2m+1), Boys at s^2 T; recurrences unchanged -> gen_vrr base
case + omega argument; validate vs libcint (int2e with omega). (2) second K
build per iteration for K_lr; F = H + J + a K + b K_lr + Vxc. (3) XC:
ExchCXX wb97x_xc.hpp, wb97x_d, lcwpbe/hse wPBEh SR exchange (GGA) through
the translator (expect a day on the attenuation helpers), functional
records with alpha, beta, omega; PySCF wB97X / CAM-B3LYP.

## Status 2026-09-04: RI-MP2 landed (src/trc_rimp2.F90)

- trc_rimp2_run(b, aux, pl, nocc, cmo, eps, res, nfrozen, aux_block, la, comm):
  (P|Q) potrf on device, (mu nu|P) in aux-shell blocks via trc_basis_t%subset,
  (ia|P) by two GEMMs per P, B = X L^-T by one trsm from the right, per i one
  GEMM K_i = B_i B^T and a (j,b) kernel with the denominators; i split over
  ranks. E_os/E_ss separate (SCS/SOS free). Water vs PySCF DF-MP2 on the same
  sets: 3e-11 (os), 8e-13 (ss).
- Bug found on the way: rebuilding an aux block through basis%build folded
  common_fac_sp a second time (s, p aux functions off by 1/sqrt(4pi),
  sqrt(3/4pi)) -> E_corr of -2.8e5. Third time this factor bit; hence subset.
- Next: C entry (trc_rimp2), aux basis from mqc (real RI sets, not the
  synthetic ladder), memory: B is nocc*nvir*naux per rank -- distribute over
  ranks like Stocks/Palethorpe/Barca for big systems; timings on cholesterol
  need a real RI basis reader in the tests.
- Deferred (Jorge, 2026-09-05): make the BSE JSON basis reader a small
  library of its own (json-fortran behind it, thin adapters in mqc and
  terco) instead of the port in src/trc_basis_json.F90. Not now.

## Status 2026-09-05 (late): what is slow and what is not

- RI-MP2 on cholesterol 6-31G with cc-pVDZ-RI (2538 aux): 0.85 s on one
  V100 (metric+Cholesky 0.01, (mu nu|P) 0.17, transform 0.10, trsm 0.01,
  pair energies 0.55). Not the problem.
- cc-pVDZ SCF on cholesterol: > 45 min on one V100 and still running. The
  Fock build over SEGMENTED general contractions: the three s columns of a
  9-primitive C/O contraction are three shells over the same primitives, so
  a core (ss|ss) quartet costs 3^4 shell quartets of 9^4 primitive quartets
  where a general-contraction-aware code pays one of 9^4. Fix belongs in
  the ERI path: contract primitive quartets once per primitive-shell
  quartet and scatter into all coefficient columns (libcint's way), i.e.
  the pair list keyed by primitive shells with a coefficient matrix per
  column. Substantial; not started.
- nvfortran host vectoriser bug on (sd|..) kernels: -Mnovect on the host
  target (747acbd). Reduced reproducer for NVIDIA still owed.
- SAD from mqc: convention rescale by common_fac_sp is terco's; the 0.013
  was an mqc averaging bug (fixed by gs2); a 0.0046 residual in O's s/d
  blocks of cc-pVDZ is with gs2 (per-block breakdown sent).

## 2026-09-05: general contraction in the ERI container (feat/general-contraction)

- ps_view_t (trc_bins): shells sharing centre, l, exponents -> one primitive
  shell with an (np x ncol) coefficient matrix, PS_NCOL_MAX = 4 columns
  (longer contractions split). Built in eri_build; primitive pairs over
  primitive shells with unit coefficients (K_ab only); Schwarz per ps pair =
  max over column pairs; bins over ps pairs; dshp = dsh folded (max over
  columns), folded inside fock_bins_ps where dsh is device-current.
- Kernels (gen_perclass): g(nv, ncmax) accumulated per column combination
  (cab x ccd weights), chunked at ncmax = min(16, 64/nv); HRR + digestion per
  combination with the contracted quartets' canonical order and weights.
  The nosym/enumerated kernel and the shared kernel keep the contracted view;
  a call without a view gets a trivial one-column view.
- Host: 45/45, water cc-pVDZ RHF vs PySCF 4e-14 through the new path.
