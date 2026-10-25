# Code structure

How the source tree is organised and the few rules that keep it that way.
`docs/ARCHITECTURE.md` explains the algorithms; this file explains where code
lives, who may depend on whom, and what names and sizes to aim for.

The rules are guidance with a tool behind them, not a gate on every line.
`python3 tools/modmap.py` draws the map, prints the layer matrix and lints
names. When a rule gets in the way, say why in the PR and move on.

## 1. Layers

A module may use modules in its own layer or any layer above it in this table
(lower rows depend on higher rows). Dependencies point towards `core`.

| Layer | Directory | Role |
|---|---|---|
| `core` | `src/core/` | Numerics and plumbing with no chemistry: Boys, tables, linalg, errors, MPI device binding |
| `basis` | `src/basis/` | Basis sets: JSON reader, decontraction |
| `grid` | `src/grid/` | DFT grids and partitions |
| `eri` | `src/eri/` | Screening, shell-pair bins, integral drivers, Fock and gradient launchers |
| `xc` | `src/xc/` | XC integration, functionals, VV10 |
| `solvation` | `src/solvation/` | PCM |
| `post_hf` | `src/post_hf/` | RI-MP2 |
| `scf` | `src/scf/` | SCF driver, DIIS, guesses, DF-JK, gradients, CPHF/CPKS |
| `api` | `src/api/` | Fortran API and C ABI |
| `app` | `app/` | The `terco` CLI |

`src/generated/` is output, not a layer: it belongs to the layer of the
generator that wrote it (`eri`, `xc`).

Two things follow from the direction of the arrows:

- **Entry points are leaves.** Nothing uses `app/`. `trc_api` and `trc_capi`
  are used only by `app/` and by tests.
- **`scf` orchestrates.** It may use everything below it. Nothing below it
  may call back into it.

### Known exceptions

The current tree has upward edges (see the "Upward dependencies" list in the
module map). They are tolerated, not endorsed, and each has a target:

| Edge | Why it exists | Target |
|---|---|---|
| Most modules in `basis`, `eri`, `xc`, `scf`, `solvation` use `trc_api` | `trc_basis_t` and friends are defined there | Move the shared types (`trc_basis_t`, `trc_pairlist_t`) to a `core`-level module so `trc_api` becomes a thin facade |
| `core/trc_ao_norm` uses `eri/trc_batch` | shared shell-count helpers | Move the helper to `core` |
| `post_hf/trc_rimp2_driver` uses `scf/trc_scf_df` | RI-MP2 reuses the DF factors | Move the DF factor build below `scf`, or pass the factors in |

New code must not add an upward edge. Removing one is always welcome.

## 2. Modules

- One module per file, named for the file: `trc_foo.F90` holds `trc_foo`.
- A module owns one idea (a type and its operations, or one stage of the
  pipeline). If you cannot say what it is in one sentence, split it.
- `use m, only: ...` always. `private` by default, then list what is public.
- Keep the public surface small. A module exporting more than about 40 names
  is usually two modules.
- Leaf data modules with no procedures (settings, constants) are fine and
  encouraged: they are what lets layers stay apart.

## 3. Naming

Names should read as a sentence fragment and tell you where to look.

- Module: `trc_<area>[_<part>]`, for example `trc_scf_diis`.
- Public procedure: `<module-area>_<noun>_<verb>`, noun before verb, so
  related procedures sort together: `scf_fock_build`, `scf_fock_assemble`,
  `scf_density_project`. Type-bound procedures drop the prefix because the
  type supplies it: `fock%build`, `context%setup`.
- Predicates start with `is_`/`has_`. Functions that return a count start
  with `n`.
- Private helpers may be short, but still say what they do.
- Types end in `_t`; keep the `trc_` prefix on the public ones.
- **Spell names out.** Readability beats brevity: `context`, not `ctx`;
  `density`, not `dens`; `iteration`, not `iter`; `threshold`, not `thr`.
  The same goes for variables, dummy arguments and types. Short names are
  fine where the maths is short and the meaning is local: loop indices
  (`i`, `j`, `ishell`), and symbols copied from a formula in the same
  procedure, with the formula named in a comment. Abbreviations from the
  domain's own vocabulary stay (`scf`, `eri`, `dft`, `diis`, `pcm`, `xc`);
  invented ones (`asm`, `blk`, `calc`, `tmp`, `ctx`, `cfg`) do not.

Existing names are renamed in batches, one module at a time, with the rename
and nothing else in the commit. `modmap.py lint` reports the names that
break the convention; the backlog is long, so new and touched code goes
first.

## 4. Size

Aim for procedures that fit on a screen or two. These are hand-written code
numbers, measured on `main`: 744 procedures, 66 over 100 lines, 18 over 200,
7 over 300, 1 over 500 (`trc_scf_run`, 560 lines).

| Size | Meaning |
|---|---|
| up to about 150 lines | Normal |
| 150 to 300 | Fine if it is one linear story (a driver, a setup). Split if it has independent phases |
| over 300 | Needs a reason in the review. Split by phase into private procedures or a small type |
| over 500 | Do not add. Existing ones get split when you next change them |

Signatures: more than 8 dummy arguments means a derived type is waiting to be
born. Bundle the arguments that travel together (settings, a work context,
the grid and its weights). `modmap.py sig` reports the count.

**Exempt:** everything under `src/generated/`, the translated functionals in
`src/xc/functionals/` (up to 2,700 lines each), and the per-class integral
kernels. Their size is a property of the generator, and they are changed by
changing the generator.

`do concurrent` item routines are exempt from signature length: the
explicit-shape dummy rules in `AGENTS.md` force long argument lists there.

## 5. Made for people and for LLMs

The same habits help both:

- A module's purpose goes in its FORD header (`!!` block) so a reader, human
  or model, can decide whether to open it from the first lines.
- Names carry the area, so a grep for `scf_fock_` finds the whole family.
- Small modules and short procedures fit in a context window; a 500-line
  procedure does not leave room for anything else.
- Layering means a change in `eri` can be understood without reading `scf`.
- `modmap.py json` is machine-readable. Point a model at it, or at
  `modmap.py calls <procedure>`, before it greps the tree.

## 6. Checking

`tools/archcheck.py` enforces sections 1 to 4 deterministically: same tree,
same output. It is regex based like `modmap.py`, so it checks terco's regular
style and does not parse all of Fortran.

| Rule | Flags |
|---|---|
| `layer-up` | a module using a module in a higher layer |
| `module-prefix` | a module without the `trc_` prefix |
| `proc-prefix` | a public procedure not starting with its module's name |
| `short-name` | a procedure, argument, local or module variable containing a word from the `SHORT_NAMES` table (`ctx`, `cfg`, `dens`, ...) |
| `proc-length` | a procedure over 300 lines |
| `arg-count` | a procedure with more than 8 dummy arguments (`routine seq` item routines excepted) |

```sh
python3 tools/archcheck.py                    # fail on new violations only
python3 tools/archcheck.py --summary          # counts per rule
python3 tools/archcheck.py --all --rule short-name
python3 tools/archcheck.py --update-baseline  # after paying debt down
```

Existing violations live in `tools/arch_baseline.json`. CI fails on a
violation that is not in it, and on a length or argument count that grew past
its baseline value. When you fix some, run `--update-baseline` and commit the
smaller file, so the debt cannot return. `SHORT_NAMES`, the limits and the
exempt paths are constants at the top of the script: change them there and
update the table in section 3 in the same PR.

`tools/modmap.py` is the exploration side:

```sh
python3 tools/modmap.py layers     # layer matrix and upward edges
python3 tools/modmap.py calls trc_scf_run --depth 2
python3 tools/modmap.py graph --dir scf > scf.mmd
```
