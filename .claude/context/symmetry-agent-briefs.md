# Symmetry work: briefs for Sonnet agents

These are the instructions for an agent (Claude Sonnet) implementing one
stage of `.claude/context/symmetry-plan.md`. The plan says *what* each stage
builds and why. This file says *how to work*: what to read, what to do in
what order, what not to touch, how to verify, when to stop, and what to
report.

A coordinator dispatches an agent with a prompt like:

> You are implementing stage **S8** of the terco symmetry project. Follow
> `.claude/context/symmetry-agent-briefs.md`: the "Common instructions"
> section, then the "S8" brief. Work in the worktree `../terco-s8`.

---

## Common instructions (every stage)

### 1. Read before writing code

In this order:
1. `AGENTS.md`, all of it. Its rules apply to every line you write.
2. `.claude/context/symmetry-plan.md`: the "Goal and scope", "Why Abelian
   groups are tractable", "Rules for every stage" sections, and your stage's
   section.
3. Your stage brief below.
4. `.claude/skills/pr-review/SKILL.md`. Your PR will be reviewed against it.
5. The files your brief lists under "Read". Read the code you will call or
   change; don't guess a signature.

### 2. Environment (this Mac)

- Apple Silicon, macOS. **No NVIDIA GPU.** You cannot run anything marked
  [GPU]. List it as "not run" in your report; never claim it passed.
- Toolchain from Homebrew in `/opt/homebrew/bin`: `gfortran-15`, `cmake`,
  `ninja`. Python tools are in the venv `~/.venvs/terco`; use
  `~/.venvs/terco/bin/python` for `fortitude`, `pyscf`, `pytest` and the
  anchor generator. If any of these is missing, stop and report; don't
  install system software yourself (Homebrew needs the user's password).
- Host build (after S0, the `macos-host` preset is preferred):
  ```sh
  cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran-15 \
        -DTERCO_ENABLE_GPU=OFF -DTERCO_ENABLE_MPI=OFF -S . -B build
  cmake --build build --parallel
  ctest --test-dir build --output-on-failure -LE LONG
  ```
- Run builds and full test suites **in the background** and wait for the
  completion notification. Don't write polling loops (see
  `.claude/context/workflow.md`). To iterate quickly, run your own checks with
  `ctest -R check_sym`.
- The OpenMP build rewrites `src/` in place. Do it only in a scratch
  worktree, never in your working tree with uncommitted changes. The recipe
  is in `.claude/context/ci.md`.

### 3. Git

- Work in your own worktree:
  `git worktree add ../terco-<stage> -b sym/<stage>-<slug> <base>`.
  `<base>` is `symmetry`, or the branch of the stage you depend on if that is
  not merged yet. Parallel agents must never share a worktree or a build
  directory.
- Make small commits whose messages say what changed and why. End each with
  `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`. No session links.
- **Do not push, open PRs or merge** unless your dispatch prompt explicitly
  says to. Never commit to `main` or `symmetry` directly. Never rewrite
  pushed history.
- **Everything targets the fork, `AlbertThie/terco` (remote `fork`).** That
  covers pushes, branches, pull requests and issues. `origin`
  (`JorgeG94/terco`) is fetch-only, and its push URL is disabled on purpose.
  Never re-enable it, never open a PR or issue against `JorgeG94/terco`, and
  with `gh` always pass `--repo AlbertThie/terco`. PRs target the fork's
  `symmetry` branch, not `main`, unless the dispatch prompt says otherwise.

### 4. Hard limits (stop and report instead of working around them)

- A test you didn't write starts failing. Don't loosen its tolerance, edit a
  reference literal, or mark it skipped.
- A `.ref` anchor needs to change: regenerate it with
  `test/pyscf_refs/generate.py <case>`, never by hand. If pyscf disagrees with
  terco beyond the gate, stop: it is a finding, not something to tune away.
- You need to edit `src/generated/`: change the generator in `scripts/`
  instead. If your brief doesn't say you may touch generators, stop and ask.
- With `symmetry=off`, any number changes (goldens, anchors, existing checks):
  stop. That breaks the stage's most basic guarantee.
- The design in the plan turns out wrong or impossible (for example, the
  pyscf convention can't be matched, or a framework type lacks something
  you need): stop and describe what you found. Don't redesign silently.
- You have spent more than 3 failed attempts on the same error: stop and
  report what you tried.

### 5. Definition of done

All of the following, run by you, with the output summarized in your
report:
1. Host build compiles with no new warnings in files you touched.
2. `ctest -LE LONG` is fully green, including your new `check_sym_*` checks.
3. Your new checks fail without your change. Show this once: revert the key
   line, run, and restore.
4. `~/.venvs/terco/bin/fortitude check` is clean on the files you touched.
5. `cd test && ~/.venvs/terco/bin/python pyscf_refs/check_staleness.py` is
   clean if you added or touched anything a `.ref` lists.
6. `python3 test/goldens/compare_goldens.py --report build/test/golden_out`
   shows no change with symmetry off.
7. The OpenMP build passes at `OMP_NUM_THREADS=4`, if you touched
   `do concurrent`, `!$acc` or `!$omp` code.
8. You ran the `pr-review` skill on your own diff and fixed every must-fix.
9. You did not edit the Status column in `symmetry-plan.md`; the
   coordinator updates it on merge (adjacent-row edits conflict).

### 6. Report format (your final message)

```
## Stage <id>: <title> — DONE | BLOCKED | PARTIAL
Branch: sym/<stage>-<slug> (base: <base>), commits: <n>, head: <sha>
### What changed
- <file>: <one line>
### Verification
- host ctest: <passed>/<total>, new checks: <names>
- new checks fail without the change: <how shown>
- fortitude / staleness / goldens / OpenMP: <result or "n/a: reason">
- [GPU] items not run: <list>
### Numbers
- <oracle results: ΔE, max off-block element, etc., with tolerances>
### Decisions and deviations from the plan
- <decision>: <why>
### Blockers / questions for the coordinator
- <or "none">
```

Report numbers you measured, not numbers you expect. If something wasn't
run, say so.

### 7. Delegating mechanical work to Haiku (`sym-helper`): the default

Subagents can't spawn subagents, so you can't start a helper yourself.
The coordinator (the main session) dispatches the `sym-helper` agent
(`.claude/agents/sym-helper.md`, model Haiku) on your behalf. **Delegate
by default**: for the task types below, send a request instead of doing
the work yourself, unless it takes under a minute or you are blocked on
the answer right now. Your report must list what you delegated and what
you kept, with a one-line reason for each mechanical task you kept.

- **Delegate these** (each is long-running or rote, and checkable):
  - **The final OpenMP run.** The helper creates a scratch worktree at your
    final commit, converts it, builds with the `macos-omp` preset, runs
    `OMP_NUM_THREADS=4 ctest -LE LONG --timeout 900`, removes the worktree,
    and returns the pass count and any failure output. Send it as soon as
    your final commit exists, then do your report and review while it runs.
  - **Full host suites** after an edit, when you can keep working meanwhile:
    `ctest -LE LONG` in your build dir, returning the pass count and
    failures.
  - Writing `.xyz` files from coordinates or internal coordinates you
    give.
  - Running a pyscf snippet you write, and returning its printed numbers.
  - Finding every call site or use of a routine, type or keyword.
  - `fortitude check` on named files, and applying its fixes to those files.
  - Regenerating one golden with the per-case method, or one `.ref` with
    `generate.py <case>`, and reporting which lines changed.
- **Never delegate**: design, correctness or numerical judgement,
  tolerances, choosing what a reference should be, editing source logic,
  generators, or git. Do those yourself, or report them as a blocker.
- **How to ask**: send the coordinator a message with `SendMessage` (to
  `main`) and keep working. Use one request per task, in this form:
  ```
  sym-helper request (stage <id>)
  Worktree: <absolute path>
  Task: <one self-contained instruction; the helper has none of your context>
  Files it may edit: <paths, or "none">
  Return: <exactly what you want back>
  ```
  If messaging isn't available, list the requests under "Blockers /
  questions" in your report instead.
- Pause work on the files a helper may edit until its result arrives, so
  you and the helper never edit the same file at once. Don't run a build
  in a directory a helper is building in.
- Check a helper's result before using it, for example by spot-running one
  test or diffing a file. You are responsible for it in your PR.

**Coordinator side.** The coordinator also uses Haiku helpers for its own
mechanical work: post-merge full suites (host, OpenMP, LMAX=3), reference
staleness and goldens reports, and conflict-free rebuilds of shared files
followed by a "reproduced exactly" diff check. Merge decisions and
conflict resolution that needs judgement stay with the coordinator.

### 8. Lessons from S1–S5 (each one cost a failed run)

- **test-drive suites run serially**: call `run_testsuite(..., parallel=.false.)`.
  The default runs cases in an OpenMP parallel do, and the basis reader is
  not thread-safe.
- **Never compare two runs bitwise**, whether energies, matrices or printed
  text. Under OpenMP (and on the GPU) the Fock build accumulates in a
  run-dependent order: identical runs differ by ~1e-13..1e-12 Eh and in the
  last printed digit. Compare with a tolerance, and justify it from a
  measurement.
- **Don't converge below the noise floor**: for ~1000 Eh systems the
  commutator floors at a few 1e-10 under OpenMP. Use `conv_diis >= 1e-9`
  and `conv_energy >= 1e-10` in oracle SCFs. The energy is variational, so
  the 1e-10 oracle gates still hold.
- **Reference hashes**: editing any file on a `.ref`'s `# files:` line makes
  it stale. Regenerate with `generate.py <case>` and confirm that only the
  hash line moved.
- **OpenMP builds** use the `macos-omp` preset, which links libomp (see
  `macos.md`), with `ctest --timeout 600`. Run your final OpenMP check on
  your final commit.
- **Test planar and nonplanar molecules**: in a planar molecule every atom
  lies on the mirror, so many symmetry errors vanish identically. Include
  a molecule whose operations move atoms.

---

## S0: macOS host and OpenMP build

**Objective**: the full suite is green on this Mac in the host build and the
OpenMP build, without changing Linux behaviour.

**Read**: `CMakeLists.txt`, `CMakePresets.json`, `cmake/`, `test/CMakeLists.txt`,
`.claude/context/ci.md`, `.claude/context/gpu-pitfalls.md` (build gotchas).

**Prerequisite (the user does this)**: Homebrew installed. If `brew` is
missing, stop immediately and report.

**Do**:
1. `brew install gcc@15 cmake ninja openblas`. Create the venv:
   `python3 -m venv ~/.venvs/terco`, then
   `~/.venvs/terco/bin/pip install fortitude-lint==0.9.2 pyscf==2.14.* mpmath==1.3.0 pytest`.
   Record every version in `.claude/context/macos.md`.
2. Configure the host build as in Common §2. Fix each failure in CMake,
   guarded with `if(APPLE)`, one commit per distinct problem. Likely areas:
   `.so` versus `.dylib` assumptions in tests (`link_c.c`,
   `check_static_terco_opt.cmake`, any `TERCO_CHECK_DEVICE_LINK` path), BLAS
   and LAPACK discovery for `pic-blas` (prefer Homebrew openblas; note what
   happens with Accelerate), `pic-mpi` with MPI off, rpath and `install_name`
   for the shared library.
3. Run `ctest -LE LONG`. For each failing test, decide whether it is a
   platform issue (fix the build) or a numerical difference. For a numerical
   difference, **stop and report** the test, the numbers and the tolerance.
4. OpenMP build in a scratch worktree, per `.claude/context/ci.md`, at
   `OMP_NUM_THREADS=4`. Also time it at 16 threads on one SCF example for the
   macos doc.
5. Add `macos-host` and `macos-omp` presets to `CMakePresets.json`.
6. Write `.claude/context/macos.md`: install steps, versions, the preset
   commands, build and test wall times, and any test skipped on macOS with
   its reason.

**Don't**: change any Fortran source to make a test pass; touch generators;
add a CI job (that's the maintainer's call; propose it in your report).

**Done when**: Common §5 holds, and `macos.md` lets a fresh agent build from
scratch.

## S0b: Apple GPU feasibility spike (research only, at most one working day)

**Objective**: a sourced go/no-go on running terco's `do concurrent` kernels
on the M5 Pro GPU.

**Do**:
1. Establish with primary sources (Apple Metal Shading Language
   specification, Metal feature-set tables) whether Apple GPUs support fp64
   arithmetic.
2. Establish which offload targets LLVM flang's
   `-fdo-concurrent-to-openmp=device` and LLVM OpenMP offload support. If
   flang is installable via Homebrew, try a 20-line `do concurrent` test
   program and record what happens.
3. Note any other route (SPIR-V to Metal, hand-written Metal) and why it
   would or would not keep the project's standard-Fortran premise.
4. Write `.claude/context/apple-gpu.md`: the question, each finding with its
   source link, and the conclusion.

**Don't**: change any product code or CMake. If the answer is
unexpectedly "go", stop and report; don't start a port.

## S1: point-group tables

**Objective**: `src/symmetry/trc_sym_group.F90` with data and logic for the
eight Abelian groups, plus `check_sym_group`.

**Read**: `test/check_si.F90` (test-drive pattern), `src/core/trc_error.F90`,
`src/CMakeLists.txt` (how sources are registered), `docs/pic_api.md`.
Look at pyscf's `pyscf/symm/param.py` in the venv for irrep names and the
operation order.

**Do**:
1. Represent operations as 3-bit sign masks, with composition by XOR. Fix the
   operation order of each group to pyscf's (`param.OPERATOR_TABLE`), and the
   irrep order and names to pyscf's (`param.IRREP_ID_TABLE` /
   `CHARACTER_TABLE`). Cite file and version in a comment.
2. `type :: trc_sym_group_t` (name, order, ops, chars, irrep_names, product
   table) and a constructor returning `error_t` for an unknown name.
3. Register the new directory in CMake. Make sure `tools/dc_to_omp.py` and
   `tools/acc_to_omp.py` `walk()` reach `src/symmetry/`, which they should
   for `.F90`. Confirm it.
4. `test/check_sym_group.F90`, test-drive, covering closure, character
   orthogonality, irrep count, names and order, and the product table
   consistent with the characters.

**Don't**: wire anything into the SCF.

## S2: detection, orientation and symmetrization

**Objective**: `src/symmetry/trc_sym_geom.F90` finds the largest Abelian
group for molecules with nondegenerate moments of inertia, orients them,
symmetrizes them, and builds the atom permutation table. In this stage only
C1 and Cs are reported; other groups are found internally but capped, with
a hook that S6–S11 enable.

**Read**: `src/core/trc_basis_types.F90` (atom data: `at_z`, `at_r`,
`at_elem`, ghost charge 0), `test/trc_test_basis.F90` (`read_xyz`),
`src/core/trc_constants.F90`, and pyscf's `pyscf/symm/geom.py`
(`detect_symm`, `as_subgroup`, and the charge-centre convention).

**Do**:
1. Centre on the centre of nuclear charge, diagonalize the inertia tensor
   (use the existing linalg helpers or a small `pure` 3×3 Jacobi), and test
   each of the 7 masks against atom equivalence. Equivalence means the same
   element, the same ghost status and the same basis assignment.
2. Return a result type: group name, rotation, origin shift, `perm(natm, nop)`,
   symmetrized coordinates and the residual. Symmetrize by averaging g·x
   over the group, and also rotate the molecule into the frame.
3. Degenerate moments (relative difference below a tolerance) return C1 with
   a `global_logger%warning` naming S12.
4. Tests in `check_sym_geom`: water and HOCl detected (water capped to Cs in
   this stage); a randomly distorted molecule gives C1; symmetrized
   coordinates are invariant to 1e-14; each `perm(:, g)` is a bijection; the
   rotation is orthogonal with determinant +1.

**Don't**: touch the basis or SCF; choose axis conventions for groups other
than Cs (those belong to their own stages).

## S3: AO representation and SO basis

**Objective**: `src/symmetry/trc_sym_ao.F90` with AO signed permutations
per operation (Cartesian and spherical), and the SO matrix U in per-irrep
column blocks.

**Read**: `src/core/trc_cart.F90`, `src/core/trc_cart_perm.F90`,
`scripts/cartorder.py`, `src/core/trc_sph.F90`, `scripts/gen_sph.py`,
`src/basis/trc_decontract.F90`, `src/core/trc_ao_norm.F90`, and the 1e
integral driver in `src/eri/trc_integrals.F90` (to get S and H for the
tests).

**Do**:
1. Cartesian signs from (lx, ly, lz) via the generated component tables, so
   both `TERCO_CART_ORDER` settings work. Build and run the test under both
   orders.
2. Spherical parities computed from the cart→sph coefficients: transform
   each Cartesian-signed component and assert the result is ±1 times the
   same harmonic. If it isn't, stop and report.
3. Shell mapping: equivalent atoms must have identical shell sequences
   (same l, exponents, coefficients). Refuse with `error_t` otherwise.
4. SALC projection, deduplication and normalization into U, with
   `irrep_off` and `irrep_dim`.
5. Anchor: per-irrep SO counts from pyscf (`mol.symm_orb`) for water (in
   Cs here), HOCl, Cartesian and spherical, 6-31G* and cc-pVDZ. Add a case to
   `test/pyscf_refs/generate.py` and generate the `.ref`.
6. `check_sym_ao`: Uᵀ U = I, Uᵀ S U and Uᵀ H U block-diagonal (≤ 1e-12),
   irrep dimensions summing to nao and matching the anchor; the l = 0..3 set
   gated on `TERCO_LMAX`, with the SKIP pattern from `test/check_nosym.F90`.

**Don't**: change basis construction or normalization (the
`common_fac_sp` rule).

## S4: SO-blocked SCF (pilot group Cs)

**Objective**: with symmetry on, the SCF diagonalizes per irrep, keeps
orbital irrep labels, and symmetrizes D and F. The energy is identical to
C1, and C1 runs are byte-identical to before.

**Read**: `src/scf/trc_scf_driver.F90` (all of it, especially
`orbitals_and_density`, the orthogonaliser, `scf_driver_null_shift`, the
gap check and failure reporting), `trc_scf_context.F90`,
`trc_scf_state.F90`, `trc_scf_settings.F90`, `trc_scf_diis.F90`,
`trc_scf_soscf.F90`, `trc_scf_incremental.F90`, `src/core/trc_linalg.F90`, and
the gradient driver in `src/scf/gradient/`.

**Do**:
1. Add an internal settings flag (the user-facing input comes in S5). When
   it is false or the group is C1, take exactly the old code path. Prove this
   with an unchanged goldens report.
2. New `src/symmetry/trc_sym_scf.F90`: transforms between AO and SO, per-block
   orthogonalization and diagonalization with one `trc_linalg_t` per block
   (release every one on every path, error paths included), eigenpair merge
   with labels, and per-irrep fixed occupations as an option.
3. Symmetrize D and F/KS each iteration. Keep DIIS working on the
   symmetrized matrices. SOSCF: try it; if it doesn't respect the blocks,
   refuse SOSCF with symmetry on through `error_t`, and explain the choice in
   the report.
4. Symmetrize the final gradient with `perm`.
5. Device data: any new array a kernel reads gets `enter data` and
   `exit data` on every path. You can't run the GPU build, so list the arrays
   and the mapping sites in the report for [GPU] review.
6. Tests in `check_sym_scf_cs`: the central oracle (symmetrized geometry, run
   with and without symmetry) for RHF HOCl, UHF HCO and B3LYP HOCl, with
   ΔE ≤ 1e-10 (HF) and ≤ 1e-9 (DFT); orbital energies equal; irrep
   occupations and per-irrep orbital energies against a pyscf anchor;
   symmetry-forbidden gradient components ≤ 1e-10.

**Don't**: touch the Fock build or `src/eri/` (that is S13–S14), or the CLI
(S5).

## S5: input, CLI, API, C ABI and reporting

**Objective**: users control symmetry and see its results; the ABI is
extended.

**Read**: `app/trc_run_input.F90`, `app/trc_run_json.F90`,
`app/trc_run_report.F90`, `app/trc_run.F90`, `src/api/trc_api.F90`,
`src/api/trc_capi.F90`, `include/terco.h`, and `git show a1fdef6 --stat` plus
the commit that bumped "ABI 1.20" (follow the same pattern).

**Do**:
1. Add the keywords `symmetry = off | auto | <group>` (default `off`) and
   `symmetry_tol`. An explicit group must be a subgroup of the detected one,
   or the run refuses.
2. Reporting: the group (detected and used), the residual, the orientation,
   a per-irrep occupation table, and an irrep label on each orbital. Nothing
   new is printed when symmetry is off.
3. Fortran API options and results; the C ABI getters and setters; the
   header; the ABI minor bump.
4. Tests: input parsing, the C round trip (pattern: `test/check_capi_*.F90`
   with `link_c.c`), and output with `symmetry=off` identical to before
   (diff a CLI run).

## S6–S11: one point group per stage

Each of these stages enables **one** group. Your dispatch prompt names it.
Use this table:

The "Ops" column is a guide only. The order that counts is pyscf's
`param.OPERATOR_TABLE`, which S1 encodes.

| Stage | Group | Ops | Test molecules | Open-shell case |
|---|---|---|---|---|
| S6 | C2 | E, C2z | H2O2 (gauche, dihedral ≈ 112°) | — |
| S7 | Ci | E, i | meso-CHFCl–CHFCl (staggered, anti) | — |
| S8 | C2v | E, C2z, σxz, σyz | water, formaldehyde | NH2 (UHF, `test/nh2.xyz`) |
| S9 | C2h | E, C2z, i, σxy | trans-N2H2, trans-glyoxal | — |
| S10 | D2 | E, C2z, C2y, C2x | ethylene twisted by 30° | — |
| S11 | D2h | all 8 | ethylene, naphthalene | O2 triplet (`test/o2.xyz`; needs S12 for linear, else skip it) |

**Do**, in order:
1. Generate the test geometries from a reliable source or an optimized
   structure (record where each comes from in the `.xyz` comment line), and
   put them under `test/`.
2. Run pyscf on each molecule with `symmetry=True`. Record the group it
   detects, its orientation of the molecule, and the irrep names. **The axis
   convention you implement must reproduce pyscf's labels.** For example,
   C2v water must give the same B1/B2 assignment pyscf gives. Verify this by
   running it.
3. Lift the cap for this group in `trc_sym_geom` and implement its axis
   choice.
4. Add the anchor case(s) to `generate.py` (group, SO dimensions per irrep,
   electrons per irrep, orbital energies per irrep, RHF and, where listed,
   UHF), and generate the `.ref`.
5. `check_sym_<group>`: detection, SALC block-diagonality in Cartesian and
   spherical bases, the central oracle for RHF, UHF and B3LYP, the anchors,
   and gradient symmetrization.
6. Add a golden for one symmetric run via `regen_goldens.py` (only for the
   new case; check that other goldens did not move).

**Don't**: modify the framework (S1–S5 code) beyond the group's own entries.
If the group exposes a framework bug, stop and report it with a minimal
reproduction; the coordinator will schedule a separate fix PR.

**Parallel note**: other agents may be doing other groups at the same time.
Touch only your group's lines in shared files (`generate.py`,
`test/CMakeLists.txt`) to keep merges trivial. Never edit the status table.

## S12: symmetric and spherical tops, linear molecules

**Objective**: molecules with degenerate moments of inertia get the same
Abelian subgroup and orientation that pyscf picks.

**Read**: `trc_sym_geom.F90` as it stands after S11, and pyscf's
`symm/geom.py` (`detect_symm`, its search for C2 axes and mirror planes, and
`as_subgroup`).

**Do**:
1. With two or three equal moments, search candidate axes from atoms, atom
   pairs and their cross products. Test each candidate frame against all 7
   masks; pick the largest group, and among frames that tie, pick pyscf's
   (match its tie-break rule, and document it).
2. Linear molecules: D2h with an inversion centre, C2v without, with the
   axis along z.
3. Tests: benzene → D2h, NH3 → Cs, CH4 → D2, N2 → D2h, CO → C2v, SF6 → D2h,
   each checking pyscf's group, orientation and labels plus the central
   oracle (RHF).

## S13: skeleton Fock design spike (no product code merged)

**Objective**: `docs/SYMMETRY_FOCK.md`, a measured comparison of the
candidate skeleton-Fock schemes in the plan and a recommendation for the
maintainer.

**Read**: `docs/ARCHITECTURE.md` sections 4 onwards in full; `src/eri/trc_bins.F90`,
`trc_binkernel.F90`, `trc_bin_tiles.F90`, `trc_screen.F90`, `trc_mdj.F90`,
`trc_mdg.F90`, `trc_eri.F90`; `src/scf/trc_scf_incremental.F90`;
`docs/GENERATORS.md`; `test/check_nosym.F90`.

**Do**:
1. Write a CPU-only prototype on a branch (`sym/s13-spike`) that is not
   meant to merge. It is plain Fortran, may be slow, and uses the existing
   integral routines to get (ab|cd). It implements each candidate scheme,
   verifies that the symmetrized G_skel equals G_full to 1e-12, and counts
   the quartets evaluated.
2. Measure on water (C2v), naphthalene (D2h) and benzene (D2h after S12) in
   cc-pVDZ: quartet counts per scheme against full, and the effect on
   density screening and exchange tiles.
3. For each scheme, list exactly which files and generators change, the
   expected compile-time cost, and the interaction with the incremental
   Fock, the MPI segment split and the MD J/K engines.
4. Recommend one scheme, with its risks.

**Then stop.** S14 must not start until the maintainer approves the design.
Say so explicitly in your report.

## S14: skeleton Fock, host implementation

**Objective**: implement the approved design from `docs/SYMMETRY_FOCK.md`,
host build first, with symmetry off unchanged.

**Read**: the approved design doc and the files it names. Read
`docs/GENERATORS.md` and the relevant `scripts/gen_*.py` in full before
changing any of them.

**Do**:
1. Follow the design. If the design is wrong somewhere, stop and report;
   don't improvise a different scheme.
2. Kernel changes go through generators only. Regenerate, commit both, and
   run the "generated sources are up to date" recipe from `ci.md`. Report the
   compile-time difference.
3. Oracle `check_sym_fock`, modelled on `check_nosym.F90`: the symmetrized
   skeleton G against the full G, ≤ 1e-12, for every group's test molecules
   × Cartesian/spherical × J/K/J+K × tiles on/off × incremental on/off, plus
   the l = 0..3 set gated on `TERCO_LMAX` (also build with `-DTERCO_LMAX=3`
   and run it).
4. The central SCF oracle for all groups with the skeleton path on.
5. The OpenMP build at 4 threads is mandatory here.
6. Device mapping: list every new array, with its `enter data` and
   `exit data` sites, for [GPU] review.

## S15: GPU validation and performance

**Prerequisite**: access to an NVIDIA host (the dispatch prompt gives the
host, the nvfortran version and the GPU). Without one, don't start.

**Do**:
1. Build per `AGENTS.md` (GPU build, the correct `TERCO_GPU_ARCH`). Run the
   full ctest suite. Every [GPU] item listed in the S4–S14 reports must pass,
   with energies equal to the host build's.
2. Performance: for benzene, naphthalene, a larger D2h acene and an
   SF6-class molecule, measure Fock time per iteration and quartets launched
   with symmetry on and off. Take nsys per-class breakdowns. State the
   system, basis, thresholds, GPU and baseline (the AGENTS.md rule), and show
   the energies are unchanged.
3. `compute-sanitizer` or the `gpu-pitfalls.md` checks on one symmetric run,
   for mapping errors.
4. Write a results section in `docs/SYMMETRY_FOCK.md`.

## S16: optional follow-ups

Not dispatchable as written. Each item (DF-JK, XC grid reduction, skeleton
gradients, default-on) gets its own design brief after S15, in the style of
S13 then S14.
