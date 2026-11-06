# Molecular point-group symmetry: project plan

Integration branch: `symmetry` (off `main` at a1fdef6, the merge of #199).
This plan is written so that each stage can be handed to an agent on its own:
"implement stage S5 of `.claude/context/symmetry-plan.md`". Every stage
lists its goal, the files it touches, what it must deliver, the tests that
prove it, and what it depends on. Read the "Rules for every stage" section
first, then your stage only. Agents follow the working instructions in
`.claude/context/symmetry-agent-briefs.md` (setup, hard limits, definition of
done, report format, per-stage brief).

## Goal and scope

- **Abelian point groups only**: C1, Cs, Ci, C2, C2v, C2h, D2, D2h. A
  molecule of higher symmetry (benzene, NH3, CH4, N2) runs in its largest
  Abelian subgroup, as in most quantum-chemistry codes.
- **Full use**: detect the group, orient and symmetrize the geometry, build a
  symmetry-adapted (SO) basis, block the SCF by irrep, label orbitals, and
  finally build the Fock matrix from symmetry-unique integrals only (a
  skeleton Fock build) for up to |G|-fold fewer quartets.
- **Off by default** until the maintainer decides otherwise. With
  `symmetry=off` (or group C1) every number must be bitwise what `main`
  produces.

Out of scope for this plan: non-Abelian groups and degenerate irreps,
symmetry in CPHF/CPKS, RI-MP2, PCM, the geometry optimiser, and vibrational
mode labels. Those paths run in C1. A stage that would need one of them
refuses with a clear `error_t` message or falls back to C1 with a
`global_logger%warning`, and never silently gives a wrong number.

## Why Abelian groups are tractable here

In a standard orientation every D2h operation is a diagonal sign matrix
diag(sx, sy, sz):

| Op | E | C2z | C2y | C2x | i | σxy | σxz | σyz |
|---|---|---|---|---|---|---|---|---|
| signs | +++ | −−+ | −+− | +−− | −−− | ++− | +−+ | −++ |

So an operation is a 3-bit mask, each subgroup is a subset of these eight
masks, and every irrep is one-dimensional with characters ±1. A Cartesian
Gaussian x^lx y^ly z^lz on atom A maps to the same function on atom g(A),
times sx^lx · sy^ly · sz^lz. Real spherical harmonics are likewise even or
odd under each axis flip. The AO representation of every operation is
therefore a signed permutation: cheap, exact, and the same code for all
eight groups. "Adding a group" mostly means detecting it, fixing its axis
convention, and proving it end to end, which is why stages S6–S11 are small.

## Hardware reality (read before S0)

The development machine is an Apple M5 Pro running macOS. terco's GPU path is
nvfortran `-stdpar=gpu`, which has no macOS release and no Apple GPU
target. Apple GPUs (Metal) also have **no fp64 arithmetic**, and every terco
kernel is `real(dp)`. Stage S0b checks this properly, but the working
assumption for the plan is:

- **Local development and testing**: host build with gfortran 15, plus the
  OpenMP CPU build on the 20 cores. AGENTS.md already says the host build
  is the reference and must match the GPU numbers.
- **GPU confirmation**: an NVIDIA machine (remote, or the maintainer's). Each
  stage marks its GPU-only checks as **[GPU]**. When no NVIDIA host is
  available, they are listed in the PR body as "not run on GPU" and
  collected in stage S15.

## Stage map

| Stage | Title | Depends on | Status |
|---|---|---|---|
| S0 | macOS host + OpenMP build | — | done (branch sym/s0-macos-build) |
| S0b | Apple GPU feasibility spike (time-boxed) | S0 | todo |
| S1 | Point-group tables (`trc_sym_group`) | S0 | todo |
| S2 | Detection, orientation, symmetrization (nondegenerate tops) | S1 | todo |
| S3 | AO representation and SO basis | S2 | todo |
| S4 | SO-blocked SCF (pilot group: Cs) | S3 | todo |
| S5 | Input, CLI, Fortran API, C ABI, reporting | S4 | todo |
| S6 | Group C2 | S5 | todo |
| S7 | Group Ci | S5 | todo |
| S8 | Group C2v | S5 | todo |
| S9 | Group C2h | S5 | todo |
| S10 | Group D2 | S5 | todo |
| S11 | Group D2h | S5 | todo |
| S12 | Symmetric and spherical tops, linear molecules | S11 | todo |
| S13 | Skeleton Fock design spike | S11 | todo |
| S14 | Skeleton Fock: host implementation | S13 + maintainer sign-off | todo |
| S15 | GPU validation and performance | S14, NVIDIA host | todo |
| S16 | Optional: DF-JK, XC grid and gradient reductions | S15 | todo |

S6–S11 are independent of each other once S5 is merged and may run in
parallel (one agent each). Do C2v (S8) early, because water is the main test
molecule.

## Rules for every stage

- Read `AGENTS.md`, `.claude/skills/pr-review/SKILL.md` and this file's
  section for your stage before writing code. Those rules (kinds from
  `pic_types`, `pic_logger`, FORD docs, `trc_` prefix, no edits in
  `src/generated/`, `enter`/`exit data` pairing) all apply here.
- New code lives in `src/symmetry/` with module prefix `trc_sym_`. Add it to
  the source list in `src/CMakeLists.txt` (or wherever `src/scf/` sources
  are registered). Symmetry setup is host code; only the parts on the SCF hot
  path (SO transforms, the skeleton Fock) need device data, and those
  follow the `do concurrent` rules in AGENTS.md.
- **The central oracle**: for a symmetric molecule, a run with symmetry and
  a C1 run *on the same symmetrized geometry* give the same energy (|ΔE| ≤
  1e-10 Eh for HF, ≤ 1e-9 for DFT unless measured otherwise), the same
  orbital energies, and the same density. Every stage from S4 on includes
  this check for RHF, UHF and one DFT functional.
- **External anchors**: point-group names, irrep labels, the number of
  electrons per irrep, and orbital energies per irrep against pyscf
  (`pyscf.symm`, `mf.get_irrep_nelec()`, `mf.mo_energy` with
  `mol.symmetry=True`). Generate them offline with
  `test/pyscf_refs/generate.py` as described in `test/pyscf_refs/README.md`.
  Never write a `.ref` by hand. pyscf must be installed in a local venv
  (S0 sets this up).
- New checks use test-drive (`test/check_si.F90` is the pattern). Name them
  `check_sym_*`.
- Branch: `sym/<stage>-<slug>` (for example `sym/s4-so-scf`), based on
  `symmetry`, or on the previous stage's branch if that is not merged yet.
  Say so in the PR body. PR target: `symmetry`. Merge commits, no squash,
  and never merge without the maintainer.
- Definition of done, for every PR:
  1. Host build green: `ctest --test-dir build --output-on-failure`, including
     the new checks.
  2. OpenMP build green at 4 threads if `do concurrent` or directive code
     was touched (recipe in `.claude/context/ci.md`).
  3. `fortitude check` clean; `check_staleness.py` clean if a `.ref` was
     added.
  4. With symmetry off, the existing goldens are unchanged
     (`compare_goldens.py`).
  5. The **[GPU]** items either run on an NVIDIA host or are listed as
     "not run" in the PR body.
  6. Run the `pr-review` skill on your own diff and fix every must-fix.
  7. Update the Status column above in the same PR.

## S0 — macOS host and OpenMP build

**Goal**: the full ctest suite passes on this Mac (Apple Silicon, gfortran
15), in both the host build and the OpenMP CPU build.

- Toolchain: Homebrew (the user installs it; it needs sudo), then
  `gcc@15`, `cmake`, `ninja`, `openblas` (or Apple Accelerate, whichever
  `pic-blas` accepts), and `pip install fortitude-lint==0.9.2` in a venv. The
  same venv gets `pyscf` (2.14.x to match the committed anchors), `mpmath`
  1.3.0 and `pytest`.
- Configure as in AGENTS.md with `-DCMAKE_Fortran_COMPILER=gfortran-15
  -DTERCO_ENABLE_GPU=OFF -DTERCO_ENABLE_MPI=OFF`. Expect macOS-specific
  breakage: `.so` versus `.dylib` assumptions in tests (`link_c.c`,
  `check_static_terco_opt.cmake`, `TERCO_CHECK_DEVICE_LINK`), linker flags,
  BLAS/LAPACK discovery in `pic-blas`, and `pic-mpi` without MPI. Fix these in
  CMake, guarded by `APPLE`, without changing Linux behaviour.
- Add a `macos-host` (and `macos-omp`) preset to `CMakePresets.json`, and a
  short `.claude/context/macos.md` with the exact commands and timings.
- Optional, and the maintainer's call: a `macos-latest` CI job.
- **Done when** the host build passes ctest `-LE LONG` and the OpenMP build
  passes ctest at `OMP_NUM_THREADS=4`, with no test skipped that Linux CI
  runs. Any test that must be skipped on macOS is listed with its reason.

## S0b — Apple GPU feasibility spike (time-boxed to one day, no product code)

**Goal**: a written answer to "can terco's `do concurrent` kernels run on the
M5 Pro GPU?", so the question is closed with evidence.

Check and record:
1. fp64 on Apple GPUs and in Metal (expected: none, which ends it for
   `real(dp)` integrals; emulated double-single would break the "numbers
   are right" rule).
2. LLVM flang `-fdo-concurrent-to-openmp=device` and its OpenMP offload
   targets (expected: NVIDIA/AMD only, no Metal).
3. Any other route (MLIR→SPIR-V→Metal, hand-written Metal kernels). Each of
   these would abandon the project's standard-Fortran premise.

Deliverable: `.claude/context/apple-gpu.md` with the findings, sources and a
go/no-go. If no-go, the plan continues as written (local CPU plus a remote
NVIDIA host). If, unexpectedly, go: stop and raise a separate proposal with
the maintainer before any code is written.

## S1 — Point-group tables

**Goal**: pure data and logic for the eight Abelian groups, unit-tested,
used by nothing yet.

- `src/symmetry/trc_sym_group.F90`: the eight operations as sign masks
  (table above); each group as its list of operations in a **fixed order**;
  the character tables; irrep names in Cotton order **matching pyscf's
  `pyscf.symm.param`** (for example C2v: A1, A2, B1, B2), so labels compare
  directly with anchors; the product table (irrep × irrep) for later use.
- `type :: trc_sym_group_t` with `name`, `order`, `ops(:)`, `chars(:,:)`,
  `irrep_names(:)`, and a constructor `trc_sym_group(name)` that returns an
  `error_t` for an unknown name.
- Tests (`check_sym_group`): each group is closed under composition (XOR of
  masks); characters are orthogonal (Σ_g χ_i χ_j = |G| δ_ij); the
  number of irreps equals |G|; names and order match a small table copied
  from pyscf (a `.ref` is overkill for static data; cite the pyscf file
  and version in a comment).

## S2 — Detection, orientation and symmetrization (nondegenerate tops)

**Goal**: given atoms (Z, coordinates, plus ghost and basis identity), find
the largest Abelian group within a tolerance, rotate to the standard frame,
and symmetrize the coordinates.

- `src/symmetry/trc_sym_geom.F90`. Algorithm: shift to the centre of nuclear
  charge (pyscf's choice, so frames agree); diagonalize the inertia tensor;
  when all three moments are distinct, the candidate symmetry elements are
  the principal axes, so test all 7 non-identity masks in that frame. Two
  atoms are equivalent only if Z, mass and **basis set** match (ghosts and
  per-atom bases count). The output is the group, the rotation matrix, the
  atom permutation table `perm(natm, nop)`, and the symmetrized coordinates
  (the average of g·x over the group).
- The axis convention for each group (which axis is C2, which plane is σv)
  must match pyscf, or the B1/B2 and B2g/B3g labels swap. Each group's
  stage (S6–S11) fixes and tests its own convention. S2 provides the hook
  and tests only C1 and Cs.
- Degenerate moments (symmetric or spherical tops, linear molecules) return
  C1 with a `warning` here. S12 replaces that.
- Tolerance: default 1e-5 Bohr, configurable, with the residual reported.
- Tests: water reported as at least Cs here (full C2v is S8); a distorted
  molecule as C1; symmetrized coordinates within tolerance of the input,
  and exactly invariant under each operation (to 1e-14) after
  symmetrization; the permutation table is a bijection per operation.

## S3 — AO representation and SO basis

**Goal**: for a `trc_basis_t` on a symmetrized geometry, the signed
permutation of every operation and the SO transformation U, as AO →
irrep blocks.

- `src/symmetry/trc_sym_ao.F90`:
  - `ao_map(nao, nop)` and `ao_sign(nao, nop)`. Cartesian: shell on atom A maps
    to the same shell index on `perm(A, g)`, and component (lx, ly, lz) takes
    sign sx^lx sy^ly sz^lz. Derive component order from the generated
    Cartesian tables so `TERCO_CART_ORDER` both ways works.
  - Spherical: the parity of each real solid harmonic under each axis flip,
    derived from the `trc_sph` coefficients. Do not hand-code a table;
    compute it and assert that the result is a pure ±1.
  - The SALC projector P_Γ = (1/|G|) Σ_g χ_Γ(g) R(g) applied to each AO, with
    duplicates removed and the result normalized. U is nao × nao,
    orthogonal, and stored as per-irrep column blocks (`irrep_off`,
    `irrep_dim`).
- Decontraction (`trc_decontract`) and pruned primitives must not break
  the shell-to-shell map. Check that equivalent atoms get identical shell
  lists, and refuse otherwise.
- Tests (`check_sym_ao`): Uᵀ U = I; Uᵀ S U and Uᵀ H_core U are
  block-diagonal (off-block elements ≤ 1e-12) using terco's own 1e
  integrals; the sum of irrep dimensions equals nao, and the per-irrep
  dimensions match pyscf's `mol.symm_orb` (anchor); run in Cartesian and
  spherical bases, 6-31G* and cc-pVDZ, and the uncontracted l = 0..3 set
  when `TERCO_LMAX ≥ 3` (skip pattern as in `check_nosym.F90`).

## S4 — SO-blocked SCF (pilot group: Cs)

**Goal**: the SCF diagonalizes per irrep and labels orbitals. The Fock
*build* is unchanged (full C1), so there is no speedup yet: this stage is
about correctness and plumbing.

- `src/scf/trc_scf_driver.F90` (`orbitals_and_density`, the orthogonaliser)
  and a new `src/symmetry/trc_sym_scf.F90`:
  - Transform F and S to the SO basis, orthogonalize and diagonalize each
    irrep block with `trc_linalg_t` (one handle per block; `%release()`
    each one), merge eigenpairs sorted by energy, and keep the irrep label
    per MO.
  - Occupations: aufbau across irreps by default; optional fixed per-irrep
    occupations (pyscf's `irrep_nelec`) for RHF and UHF.
  - Symmetrize the density each iteration (D ← (1/|G|) Σ R D Rᵀ), and the
    Fock/KS matrix as well. This hides grid noise in DFT.
  - Canonical orthogonalization with dropped functions is done per block.
  - DIIS, SOSCF and the incremental Fock must keep working. If SOSCF cannot
    respect the blocks, it refuses with symmetry on, or falls back to C1 with
    a warning. Document which you chose and why.
  - Gradients: symmetrize the final gradient over the group, with the atom
    permutation.
- Enable only C1 and Cs in this stage.
- Tests (`check_sym_scf_cs`): the central oracle on a Cs molecule
  (suggested: HOCl, or formamide only if it is exactly planar after
  symmetrization), RHF, UHF (an open-shell Cs radical such as HCO) and B3LYP;
  irrep occupations and orbital energies against a pyscf anchor; the
  gradient components that symmetry forces to zero are ≤ 1e-10. **[GPU]**
  the same checks on the device build, including `enter`/`exit data` for
  the SO blocks.

## S5 — Input, CLI, API, C ABI and reporting

**Goal**: users can ask for symmetry and see what it did.

- Input keyword `symmetry = off | auto | <group>` (default `off`) and
  `symmetry_tol`, in `app/trc_run_input.F90` and `app/trc_run_json.F90`.
  An explicit group must be a subgroup of the detected one, or the run
  refuses.
- Report: the detected and used group, the orientation, the symmetrization
  residual, a per-irrep occupation table, and an irrep label on every
  orbital in the orbital printout (`app/trc_run_report.F90`).
- `trc_api` / `trc_capi`: options in, group and labels out, and an ABI minor
  bump in `include/terco.h` following how earlier PRs bumped it (see the
  #199 merge, "ABI 1.20").
- Tests: CLI input parsing (`check_cli_input`-style), a C ABI round trip,
  and `symmetry=off` giving output byte-identical to before, apart from new
  report lines that only appear when symmetry is on.

## S6–S11 — one group per stage

Each of these stages is the same recipe for one group. Do not touch the
framework unless the group exposes a bug in it; if it does, fix it in a
separate small PR first.

1. Enable the group in `trc_sym_geom` and fix its axis convention to match
   pyscf. Verify this by running pyscf on the test molecules and comparing
   the orientation, not by reasoning about it.
2. Add test molecules under `test/` as `.xyz` files.
3. pyscf anchors (generated, never hand-written): group name, per-irrep SO
   dimensions, irrep occupations, and orbital energies per irrep for RHF;
   also UHF where an open-shell example exists.
4. The central oracle for RHF, UHF and B3LYP; gradient symmetrization; and
   SALC block-diagonality in Cartesian and spherical bases.
5. A golden for the symmetric run.

| Stage | Group | Order | Suggested molecules |
|---|---|---|---|
| S6 | C2 | 2 | H2O2 (gauche) |
| S7 | Ci | 2 | meso-CHFCl–CHFCl (staggered, anti) |
| S8 | C2v | 4 | water (`test/water.xyz`), formaldehyde, NH2 radical (`test/nh2.xyz`, UHF) |
| S9 | C2h | 4 | trans-diazene (N2H2), trans-glyoxal |
| S10 | D2 | 4 | twisted ethylene (e.g. 30°) |
| S11 | D2h | 8 | ethylene, naphthalene, O2 triplet (`test/o2.xyz`, UHF), once S12 lands |

## S12 — Symmetric and spherical tops, linear molecules

**Goal**: degenerate inertia moments no longer fall back to C1. Pick the
Abelian subgroup and axes that pyscf picks.

- With two equal moments (benzene, NH3), search for C2 axes and mirror
  planes among directions defined by atoms and atom pairs. With three equal
  moments (CH4, SF6), do the same, and expect D2 or D2h for Td and Oh. For
  linear molecules use C2v or D2h, as pyscf does.
- Tests: benzene → D2h, NH3 → Cs, CH4 → D2, N2 → D2h, CO → C2v, SF6 → D2h,
  each matching pyscf's chosen subgroup and orientation, with the central
  oracle on each.

## S13 — Skeleton Fock design spike

**Goal**: a reviewed design for building G(D) from symmetry-unique integrals
in terco's bin/segment architecture, before any kernel or generator changes.
The deliverable is a design document and a CPU prototype; it is not merged
as product code.

The difficulty: terco never forms a quartet list. Quartets come from
admitted bin pairs of shell pairs, the kernels fold the 8-fold permutational
symmetry, and screening is density-weighted. The classic petite-list
(Dupuis–King) works on unique quartets with weights; terco needs a variant
that keeps launches branch-free. Candidates to evaluate:

1. **Unique bra pairs**: restrict the ab pairs to orbit representatives
   under G with weight |orbit(ab)|, keep cd over all pairs, then
   G = (1/|G|) Σ R G_skel Rᵀ. This changes only the bin planner plus one
   weight per pair. Cost: it breaks the ab↔cd fold, so the net gain is
   about |G|/2.
2. **Unique pair-of-pairs**: orbit representatives of (ab, cd) with the
   triangular fold kept. This gains the full |G|, but the weights are per
   quartet, which means the generated kernels must change.
3. **Unique atom-quartet pre-pass**: weights per bin pair from atom orbits.
   Coarser, and only exact if bins respect atom orbits.

For each candidate: the change to `trc_bins`, `trc_binkernel`, exchange
tiles, the MD J/K engines (`trc_mdj`, `trc_mdg`), `scripts/gen_*.py`, the
incremental Fock, MPI segment splitting and screening bounds; and a CPU
prototype measuring quartet counts and verifying that the symmetrized
G_skel equals G_full to 1e-12 for C2v water, D2h naphthalene and D2h benzene
in cc-pVDZ. Write it up as `docs/SYMMETRY_FOCK.md` (or a section of
`docs/ARCHITECTURE.md`), with the recommended option and its compile-time
cost.

**Stop here for maintainer sign-off.**

## S14 — Skeleton Fock: host implementation

**Goal**: implement the approved design, with the host build first.

- Any kernel change goes through the generator (`scripts/`). Regenerate,
  commit both, and state the compile-time cost.
- The oracle extends `check_nosym`: symmetrized skeleton G equals full G to
  1e-12 for every group, every class up to `TERCO_LMAX`, Cartesian and
  spherical, J-only, K-only and J+K launches, with exchange tiles on and off,
  and the incremental Fock on and off.
- The central SCF oracle for every group's test molecules, now with the
  skeleton path on.
- The OpenMP build at 4 threads (races only show there).
- `symmetry=off` keeps the old launch plan exactly, with the goldens
  unchanged.

## S15 — GPU validation and performance [GPU]

**Goal**: everything marked [GPU] in S4–S14 runs on an NVIDIA device, plus a
performance report.

- Run the full ctest suite on the device build, with energies equal to the
  host build's.
- Measure (the rules from AGENTS.md apply): system, basis, thresholds, GPU,
  and baseline. Quartets launched and Fock time per iteration with symmetry
  on and off for benzene, naphthalene, a D2h acene and an SF6-class
  molecule, with nsys per-class breakdowns. Report speedup against the
  ideal |G|.
- MPI split (`TERCO_ENABLE_MPI`) with symmetry on, if a multi-GPU host is
  available.

## S16 — Optional follow-ups

Each of these is its own design-then-implement pair, decided after S15:

- DF-JK with symmetry (only unique auxiliary-index blocks).
- XC: integrate only over the grid points of unique atoms with weights. The
  grid must then be symmetric, which Lebedev grids in the standard frame
  are, but pruning and Becke partitions must be checked.
- Skeleton gradients.
- Turning symmetry on by default (maintainer decision).

## Open questions for the maintainer

1. Is a `symmetry` integration branch acceptable, or should each stage PR
   target `main` directly?
2. Should there be a macOS CI job after S0?
3. Where should NVIDIA [GPU] checks run, and who runs them?
4. Should orientation be forced to pyscf's convention, or is a terco
   convention acceptable, with a label mapping in the anchors?
