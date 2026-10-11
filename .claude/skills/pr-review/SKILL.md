---
name: pr-review
description: Review a terco pull request — correctness of integrals and SCF first, then the do concurrent / OpenACC rules, generated code, tests and performance claims. Use when asked to review a PR or a diff in this repository.
---

# Reviewing a terco pull request

terco is a GPU electron-repulsion-integral and SCF library written in standard
Fortran. Its compute kernels are `do concurrent`, compiled by nvfortran with
`-stdpar=gpu -acc -gpu=mem:separate`. The project's claim is that standard
Fortran reaches throughput competitive with hand-written CUDA, so a review
guards two things: the numbers are right, and the code stays portable
standard Fortran.

## How to review

1. Read the PR description, then the diff (`gh pr diff`). Read the files around
   each change, not only the hunks: most terco bugs are a mismatch between a
   hunk and code it doesn't show (a missing `exit data`, a generator changed
   without its output, a caller not updated).
2. Rank findings by severity: correctness, then GPU rules, then generated code,
   then tests, then performance claims, then style. Report problems, not
   praise. If there are none, say so in one line.
3. Post line-specific findings as inline comments and one summary comment.
   Every finding states the concrete failure: which input, what goes wrong.
   Don't post a finding you can't tie to a line and a consequence.
4. Never suggest a hand edit to `src/generated/`.

## 1. Correctness (most severe)

- **Energies move only with a reason.** A change that alters an SCF energy or
  an integral needs a reference that justifies the new number: pyscf, libcint
  or a stored value in `test/`. "Converges" is not validation.
- **Tolerances.** A loosened tolerance in `test/` is a defect unless the PR
  argues for it with a measurement. Watch for an assertion weakened to make a
  failing test pass, and for reference literals edited to match new output.
- **Normalisation.** libcint's `common_fac_sp` (s: 0.2820947917738781,
  p: 0.4886025119029199) is folded into the coefficients exactly once, in basis
  construction. Applying it twice, or compensating elsewhere, breaks the
  absolute Schwarz screen silently (this has happened). At the C API, unit
  normalised AOs follow `trc_set_ao_norm`.
- **Cartesian component order.** The output order is a build option
  (`TERCO_CART_ORDER`: lexicographic or GAMESS). Internal enumerations stay
  lexicographic; only output layouts follow the option. GAMESS order also
  carries per-component norms (d: √3; f: √5, √15). Any new layout must go
  through `trc_cart_perm`, and generated files must carry the
  `TRC_CART_ORDER` guard.
- **Units and data.** Å→Bohr uses CODATA 2018 (0.529177210903 Å). Basis sets
  come from Basis Set Exchange JSON. Comparisons against another code are
  meaningful only with the same basis file and the same constant.
- **Screening.** Schwarz plus density-weighted screening may only drop what
  can't matter at the threshold. Check:
  - degeneracy factors (½ for i=j, k=l, ij=kl);
  - the factor 4 on J in the density test;
  - that a coupled-perturbed solve passes `density_screen=.false.`, because a
    density-keyed screen on a trial vector changes the operator between
    matvecs.
- **SCF.** Convergence tests, DIIS/ADIIS history handling, the linear-
  dependence threshold in canonical orthogonalisation, and incremental Fock
  builds. An incremental build at loose screening drifts, so it needs periodic
  full rebuilds.

## 2. GPU rules

- **Kernels are standard `do concurrent`.** No CUDA Fortran, no CUDA C, no
  streams, no OpenACC `async`, no OpenACC compute constructs
  (`parallel`/`kernels`/`loop`) in the integral or XC path. OpenACC is used
  only for data movement (`enter`/`exit data`, `update`, `host_data`),
  `routine seq`, and `atomic`.
- **Data lifetime.**
  - Every `!$acc enter data` has a matching `exit data` on every path out,
    including early returns and error paths.
  - Every array a kernel reads is present on the device. With
    `mem:separate`, a missing mapping is an illegal address at run time, not
    a compile error.
  - Arrays updated on the host after they were mapped need an
    `!$acc update device`.
  - Reallocating a mapped array invalidates its mapping.
- **Item routines** called from a `do concurrent` body must be `pure`, carry
  `!$acc routine seq`, and live in the same module as their driver so they
  inline. A separate device procedure measured 2.4× slower.
- **No assumed-shape dummies** in anything called from `routine seq`. The
  descriptor walk per call is a measured slowdown. Use explicit-shape arrays.
- **Per-thread workspaces** (`v`, `g`, `vbuf`, the Boys `f`) must be declared
  inside the item routine. Declared alongside the loop, they become shared,
  the compiler copies them per launch, and threads race.
- **Atomics.** Updates to J/K/Fock from a kernel use `!$acc atomic update`.
  Prefer the folded block digestion: accumulate a block in registers, then one
  atomic per element.
- **Don't cap registers** (`-gpu=maxregcount`) to buy occupancy. Measured: it
  trades registers for local-memory spills and is slower.
- **Derived types and `enter data`.** Adding allocatable components to a
  derived type that is mapped whole by `!$acc enter data copyin(obj, ...)`
  has triggered an nvfortran 26.5 internal compiler error. Prefer flat
  arrays.
- **The host build must stay correct.** CI has no GPU. It builds for the host,
  where `do concurrent` is a serial loop and OpenACC directives are inert, and
  runs the whole suite. Code that is correct only on the device is a bug.

## 3. Generated code

- `src/generated/` is produced by the scripts in `scripts/`. These are
  defects:
  - a hand edit to a generated file;
  - a generator change committed without the regenerated output;
  - generated output that doesn't match its generator.

  CI's `generators` job regenerates and diffs.
- Generator changes should keep the emitted code readable: consistent
  indentation, no dead arguments, comments that match what is emitted.
- Watch compile time. The per-class kernels already take 15–20 minutes to
  build, and (dd|dd) alone is the longest. A generator change that grows the
  heavy classes' code or compile time needs that cost stated in the PR.
- A generator may use mpmath or sympy only if the outputs are committed and
  CI is told in the same change.

## 4. Tests

- New behaviour needs a check in `test/` that would fail without the change,
  against an external reference where one exists.
- Tests must pass on the host build.
- Don't accept a test that only checks a run completes.

## 5. Performance claims

terco's changes are justified by measurement. When a PR claims a speedup:
- it should state the system, basis, screening threshold, guess, GPU and the
  baseline, and show that the energies are unchanged;
- per-class claims should come from nsys or ncu, not wall time alone;
- be sceptical of a speedup that coincides with fewer quartets surviving
  screening: that's an accuracy change, not a speedup.

## 6. Style

- Match the surrounding code. terco's comments explain why, often with the
  measurement that decided it. Flag comments that restate the code, describe
  history instead of the code, or contradict it.
- Naming follows the file: `trc_` module prefix, short physics names
  (`zeta`, `eta`, `pqx`) inside kernels.
- Don't raise nits a formatter would settle.

## Output

1. Inline comments, one per finding, each with the failure scenario.
2. One summary comment: findings by severity, one line each, linking the
   inline comments. Then a single line on what you checked and found clean.
   No praise, no restating the PR.
