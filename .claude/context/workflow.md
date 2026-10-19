# Workflow: tests, references, PRs

## Writing a test that can fail

- Compare against an external reference where one exists (pyscf, libcint,
  finite differences), not only against terco's previous output.
- **Open-shell and near-degenerate cases** (for example OH's degenerate pi
  pair) need a fixed initial density in both terco and the reference.
  Otherwise each run converges to a different state. For a tensor whose
  orientation depends on that choice, compare invariants (eigenvalues,
  trace).
- **Gradients against pyscf** need SCF convergence tighter than the
  gradient tolerance (for example a DIIS error of 1e-11). Otherwise CPU and
  GPU converge to slightly different densities, which looks like a
  gradient bug.
- **Finite-difference gradients.** A plain central difference is often
  limited by truncation error. Richardson extrapolation,
  `(4*D(h/2) - D(h))/3`, separates truncation from a real error.
  Translational invariance alone does not catch a wrong factor or a sign
  error that is symmetric across atoms.
- **Ill-conditioned DF metrics.** A metric with a condition number around
  1e14 (for example a Cartesian universal-jkfit on mid-size molecules)
  amplifies ~1e-14 relative integral differences into ~1e-10 to 1e-9
  differences against another code. Plain Cholesky is the stable solve; the
  sensitivity is in the fit itself, so set the gate from a measurement.
- **Tolerances.** Never loosen one to make a test pass. If a gate sits at a
  genuine numerical floor, measure the floor, explain it in the test, and
  get it agreed in review.
- **Test the refusals.** When an API refuses a case (VV10, multi-rank,
  meta-GGA, ...), a cheap test that it still refuses keeps the path from
  regressing.

## pyscf references after an edit

1. Run `python3 pyscf_refs/check_staleness.py` from `test/` to list the
   stale cases.
2. If only comments or code that cannot change the probe moved, re-stamp
   the hash, and say in the PR why the values still hold.
3. Otherwise regenerate with `generate.py <case>` (pyscf, offline), check
   that the data moved only at round-off or for the reason you intended,
   and say what moved it.
4. A hash computed against an uncommitted `generate.py` will not match.
   Commit the generator first, then stamp.

## Stacked PRs

- Base each PR on the branch it builds on, and name that PR in the body.
  When the base merges, retarget the next PR to `main` (or delete the
  merged branch, which retargets it).
- Merge from the bottom of the stack up. After a base merges, the next PR
  merges `main` in. Use a merge commit, and don't rewrite a branch someone
  else has pulled.
- Before merging, list the other open PRs that touch the same files. Two
  PRs that each add a line at the same spot in a CMake source list will
  conflict trivially: keep both lines. A refactor that removes or renames
  names can auto-merge into code that no longer compiles, with no conflict
  marker. Grep the other PRs' diffs for the old names, and let CI build any
  merge before trusting it.

## Running long jobs as an agent

- Run long builds and test suites in the background and wait for the
  completion notification.
- Don't write polling loops like `until ! pgrep -f X` or `ps | grep X`:
  they match their own command line and never exit.
