---
name: pr-review
description: Review a terco pull request — correctness of integrals and SCF first, then the do concurrent / OpenACC rules, generated code, tests (including stale pyscf reference hashes), performance claims, kinds from pic_types, pic reuse over local copies, overlapping PRs and merges across refactors, fortitude lint, and FORD-only doc comments (no Doxygen). Use when asked to review a PR or a diff in this repository.
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
   praise. If there are none, say so in one line. Always also run the
   documentation check in section 7, the kinds check in section 8 and the
   pic-reuse check in section 10: all three are must-fix whatever else the
   PR contains.
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
- **A `do concurrent` body calls only intrinsics and same-module `routine
  seq` item routines**, never a procedure from another module.
- **Derived-type components in a body are fine** (`st%fock(i, j)`), but with
  `mem:separate` the parent object must be on the device and the component
  attached: `!$acc enter data copyin(st)` and then `copyin(st%fock)` (the
  pattern `trc_linalg` uses for `this%dwork`). Check that every component a
  body reads or writes is mapped this way, and released in reverse order on
  every path out. A missing attach is an illegal address at run time, not a
  compile error, and the host build won't show it.
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
- **MPI and internal procedures.** A program or procedure that makes a
  collective MPI call must have no internal (`contains`) procedures at all:
  under NVHPC 26.5 with HPC-X, even an unused sibling internal procedure
  zeroed the two-electron part of the Fock build at more than one rank (see
  app/main.F90's `terco_run` header). New MPI-path helpers go in modules.
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
- **Committed pyscf references go stale on ANY edit to their inputs.** Each
  `test/pyscf_refs/*.ref` records a hash of `generate.py` plus every file on
  its `# files:` line, which usually includes the test program that writes
  its probe (e.g. `test/scf_rimp2.F90` for `rimp2.ref` and
  `rimp2_ccpvdz.ref`). `terco/check_pyscf_refs_fresh` fails CI on every host
  job until they are refreshed, even when the edit cannot change the numbers
  (this turned #65 red). So:
  - list every `.ref` whose `files:` names a file the PR touches
    (`grep -l '<path>' test/pyscf_refs/*.ref`); the PR must update them,
    or `python3 test/pyscf_refs/check_staleness.py` (run from `test/`) must
    still pass;
  - if only the hash line changes, the PR should say why the values still
    hold (the probe the test writes is unchanged), ideally checked by
    regenerating with `generate.py <case>`;
  - if the values change, they must come from `generate.py`, never a hand
    edit, and the PR must say what moved them.
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

## 7. Extras: documentation is FORD, never Doxygen

terco documents its Fortran with FORD. FORD attaches a `!!` block to the
statement **before** it, so a doc comment goes on the line(s) after the entity:
after a declaration (or inline, `real(dp) :: x  !! doc`), and indented one level
after a procedure, module, submodule, program, type or interface header. The
repository was converted from a Doxygen/FORD mix to FORD only, so don't let
Doxygen back in. These rules are not style nits: report every violation as
must-fix, even when the PR has nothing else wrong.

- **No Doxygen markers or tags.** A `!>` or `!<` doc comment, or a Doxygen tag
  (`@brief`, `@param`, `@return`, `@details`, `\param`, `\brief`, ...) in any
  added line is a defect, in `src/`, `include/`, `test/` and in the Fortran
  templates inside `scripts/gen_*.py`. Directives (`!$acc`, `!$omp`) are not
  comments and aren't affected.
- **No `!!` block above its entity.** That is Doxygen placement with FORD's
  marker. FORD then attaches the text to the previous declaration, so it
  silently documents the wrong variable. Check every added or moved `!!` block
  against the line above it, not only the line below:
  - a block between two declarations belongs to the one *above* it;
  - a block above `module`, `program` or a procedure header attaches to
    nothing useful and belongs under that header;
  - a heading for a group of constants ("Cutoff profile") is a plain `!`
    comment above the group, not a `!!` block.
- **Generated code.** A Doxygen marker or a misplaced `!!` in `src/generated/`
  is fixed in its generator, never in the output (see section 3).

To find them, run these on the diff's added lines:

```sh
gh pr diff <N> | grep -nE '^\+.*!(>|<)'
gh pr diff <N> | grep -nE '^\+.*[@\\](brief|param|returns?|details|note|see)\b'
```

Then read each added `!!` block with its surrounding lines to check its
placement; grep cannot tell which declaration a block is about.

## 8. Kinds from pic_types, and the lint CI enforces

- **Kinds come from `pic_types`.** `dp`, `sp`, `int32`, `int64` (and
  `default_int`, `int_index`) are `use pic_types, only: ...`. These are
  must-fix:
  - `dp` taken from another terco module (`use trc_boys, only: dp` was the
    old way; trc_boys keeps its `dp` private now, so it also fails to build);
  - `iso_fortran_env`'s kinds (`real64`, `int64`) or a literal kind
    (`integer(kind=8)`, `real(8)`, `1.0d0`) in new code. Its units
    (`output_unit`, `error_unit`) are fine;
  - `iso_c_binding` kinds (`c_int`, `c_double`, `c_ptr`) away from the C
    boundary. They belong in `bind(c)` interfaces and `trc_capi`, not in
    computational code.

  To find candidates on the added lines (read each hit: an
  `iso_fortran_env, only: output_unit` line is fine):

  ```sh
  gh pr diff <N> | grep -nE '^\+.*(use trc_boys, only: *dp|iso_fortran_env|kind *= *[0-9]|real\([0-9]\)|[0-9]d[-+]?[0-9])'
  ```
- **fortitude must pass.** CI's `fortitude lint` job runs
  `fortitude check` (fortitude-lint 0.9.2, rules in `fpm.toml`'s
  `[extra.fortitude.check]`; `test/` and `src/generated/` are excluded) and
  fails the PR. Run it on the PR's tree. A PR must fix the code, not add
  ignores to `fpm.toml` or per-line suppressions, unless it argues for the
  rule change. The common ones in new code: one statement per line (S082),
  double-quoted strings (S241), one entity kind per declaration (S262),
  `result()` on functions (S221), named constants for array sizes (C031).

## 9. Commits

- No claude.ai session links (`Claude-Session:` lines or
  `claude.ai/code/session_...` URLs) in commit messages or PR bodies. A
  `Co-Authored-By:` trailer is fine.
- Never mention "kgen" in commits, comments, docs or file headers.

## 10. Use pic, don't re-roll it

terco builds on pic, and `docs/pic_api.md` lists everything it provides at
the tag terco pins (generated by `tools/pic_api.py`; CI fails when it goes
stale). New code uses pic rather than a local copy. These are must-fix
unless the PR says why pic cannot do it:

| Local code | Use instead |
|---|---|
| `system_clock`, `cpu_time`, `date_and_time`, `omp_get_wtime`, a hand-rolled `wall()` | `pic_timer` (`timer_type`), `pic_clock` (`monotonic_ms`, `monotonic_us`) |
| an `itoa` / `int2str` / `to_str` helper, `write (buf, "(i0)")` into a scratch string | `pic_strings: to_string` |
| trimming, padding, prefix/suffix tests, joining, splitting on hand-written loops | `pic_strings` (`strip`, `padl`/`padr`, `starts_with`/`ends_with`, `join`, `replace_all`), `pic_tokenizer` |
| a sort or argsort (insertion sort, `sort_desc`, `sort_by_*`) | `pic_sorting` (`sort`, `sort_index`, `ord_sort`, `radix_sort`; `reverse=` for descending) |
| mean / variance / min-max loops | `pic_statistics` |
| a tolerance compare (`abs(a-b) < tol`) repeated as a helper | `pic_math: is_close` |
| fixed-size "grow by copy" arrays, linear-search lookup tables | `pic_vector`, `pic_hash_map`, `pic_heap` |
| `print` progress or diagnostic chatter | `pic_logger` |
| per-step timing or flop counting | `pic_profiler`, `pic_flop_recorder`, `pic_flop_rate` |

- **Exempt: anything that runs on the device.** pic is host code, so a `pure` helper called from a `do concurrent` body (a sort of a few register values, a device-side formatter) stays local. Such a line carries `! pic-reuse: <reason>`; check the reason is true.
- **Exempt: physical constants.** A local `pi` or conversion factor is fine.
- **CI:** the `pic reuse` job runs `tools/pic_reuse_check.py` on the PR's
  added lines and fails on the mechanical cases (clock calls, new
  `itoa`/sort/statistics routines). It cannot see the rest of the table;
  the reviewer does. Before flagging, check `docs/pic_api.md` that pic
  really has it at the pinned tag.
- To find candidates on the added lines:

  ```sh
  gh pr diff <N> | grep -nE '^\+.*(system_clock|cpu_time|date_and_time|omp_get_wtime|function +(itoa|int2str|to_str)|(subroutine|function) +(sort\w*|\w+_sort)\b|"\(i0\)")'
  ```

## 11. Overlapping PRs and merges across a refactor

A clean textual merge is not a clean merge. When #71 (a print in the SCF
loop) and #72 (the SCF driver split into modules) were both open, git
auto-merged #71's banner into the refactored driver, where the helpers it
called were no longer imported: no conflict marker, and code that does not
compile. So:

- **List the overlap.** For the PR under review, list the other open PRs
  that touch the same files:

  ```sh
  gh pr list --state open --json number,title,files \
    --jq '.[] | select(.number != <N>) | "\(.number) \(.title): \([.files[].path] | join(" "))"'
  ```

  For each one that shares a file, say which should merge first and what
  the second will have to redo. Prefer landing the larger refactor first and
  rebuilding the small change on top. This is a summary item, not a must-fix.
- **A refactor's removed names.** When the PR moves, renames or stops
  importing procedures or variables, list the old names and grep the other
  open PRs' diffs for them (`gh pr diff <M> | grep -nE '^\+.*\b(name1|name2)\b'`).
  Any hit is a merge that will auto-resolve into a build failure; name the
  PR and the line.
- **A PR that merges `main` in** (a merge commit in its history): read the
  merge commit's own diff (`git show <merge>`), not only the PR's net diff.
  For every auto-merged hunk in a file the other side restructured, check
  that each name it uses still exists and is still imported there. Then
  rely on CI to build it; never approve a merge resolution that CI has not
  built.

## Output

1. Inline comments, one per finding, each with the failure scenario.
2. One summary comment: findings by severity, one line each, linking the
   inline comments. Then a single line on what you checked and found clean.
   No praise, no restating the PR.
