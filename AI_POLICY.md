# Policy on AI-assisted contributions to terco

Adapted, with thanks, from the MOM6 consortium's policy on AI-generated code.

## Guiding principle

AI coding assistants are welcome in terco the way any other tool is: a human
contributor is responsible for everything they submit. AI use is fine when
it supports the contributor's understanding of the code. It is not a
substitute for that understanding.

## Contributor responsibility

Every contribution has a human author who is accountable for its
correctness, its numbers and its conformance to terco's conventions
(`AGENTS.md`, `.claude/skills/pr-review/SKILL.md`). Every contribution has a
clear purpose, stated in its commit messages and pull request description.

A reviewer may always ask "Do you understand why this works?" and "Why are
you suggesting this?", and a satisfactory answer is always expected. Code
the contributor cannot explain should not be submitted.

## Disclosure

Say so when AI tooling materially shaped a commit or pull request. Either
is enough:

- a `Co-Authored-By:` trailer naming the tool on the commit, or
- a note in the pull request description, e.g. "An AI tool was used in
  making this contribution."

This is for the reviewer's benefit; the human contributor remains the
author and the one accountable. Links to AI chat sessions do not belong in
commits or pull requests.

## Risks of AI-assisted code

Code that compiles and passes the existing tests is not thereby correct.
The failures that matter in terco are quiet ones: an energy that is wrong in
the eighth digit looks exactly like a right one. AI tools write fluent
Fortran; they do not understand the quantum chemistry or the GPU memory
model. Contributors and reviewers should watch for:

- **Numerical correctness.** Normalisation applied twice, the wrong Bohr
  constant, a Cartesian component order or sign convention that differs
  from the reference, a screening bound that is too tight and drops real
  contributions. An energy, gradient or integral changes only with an
  external reference that justifies the new number (pyscf, libcint, or a
  stored value); "the SCF converges" is not validation.
- **Correct on one target only.** The host build must give the device
  build's numbers. A missing data mapping or an unsupported construct inside
  `do concurrent` can pass on the CPU and fail, or silently miscompute, on
  the GPU, and the reverse happens with compiler bugs.
- **Vacuous tests.** A test that only checks that a run completes is not a
  test. A new check must fail without the change it covers; show that it
  does (break the code, watch the test fail).
- **Comprehension debt.** Code its submitter does not understand gets
  harder to change with every later edit, especially when those edits are
  also AI-assisted.
- **Churn.** Changes without a compelling purpose, including performance
  changes that were not measured first, cost review time and stability.
- **Unsafe code.** AI tools can be steered by prompt injection or reproduce
  unsafe patterns. Review for logic that does not follow from the stated
  purpose: unexpected network or file-system access, changes to CI,
  generators or pinned dependencies that the change does not need. Treat
  anything that cannot be explained by its purpose as suspect.

## Roles for AI in the workflow

Encouraged, with human oversight:

- **Development:** exploring the code, planning changes, drafting
  implementations, explaining unfamiliar code, suggesting refactors.
- **Debugging and profiling:** narrowing down failures, reading compiler,
  sanitizer and profiler output.
- **Documentation:** drafting FORD comments, pages and descriptions.
- **CI triage:** summarising logs, proposing fixes for straightforward
  failures.
- **Review:** pointing out possible bugs or departures from terco's
  conventions. Only a human maintainer approves and merges a pull request,
  after understanding what it does.

Needing particular care:

- **Autonomous contribution** with minimal human oversight is not accepted:
  agents may draft, test and push branches, but every change is reviewed
  and merged by a human maintainer.
- **AI-written tests** must be shown to test something real.

## Equity

Nothing in terco's workflow requires paid AI access. `AGENTS.md` and the
files under `.claude/` are written to be read by people as well as by
tools, every CI check runs without an AI service, and the optional
`@claude` assistant on issues and pull requests is never a required check.

## Licensing

terco is MIT licensed (`LICENSE`; the parts named in `NOTICE` are
Apache-2.0). By submitting a contribution you confirm that you have the
right to license it under those terms, including any part an AI tool helped
write, and that your use of the tool complies with your employer's or
institution's policies. The legal status of AI-generated code is still
unsettled in many jurisdictions; if in doubt, ask before submitting. terco
does not endorse any particular AI vendor or tool.
