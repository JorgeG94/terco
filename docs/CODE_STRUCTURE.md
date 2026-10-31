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
| `scf` | `src/scf/` | SCF driver, DIIS, guesses, DF-JK, gradients, CPHF/CPKS |
| `post_hf` | `src/post_hf/` | RI-MP2 energy and gradient (built on the SCF's results, Z-vector and gradient pieces) |
| `api` | `src/api/` | Fortran API and C ABI |
| `opt` | `src/opt/` | The geometry optimiser (DL-FIND via libdlfind), built as `libterco_opt` with `trc_opt_capi.F90`'s own C ABI (`include/terco_opt.h`) -- only with `-DTERCO_ENABLE_DLFIND=ON`, see that option's own comment in the root `CMakeLists.txt` |
| `app` | `app/` | The `terco` CLI |

`src/generated/` is output, not a layer: it belongs to the layer of the
generator that wrote it (`eri`, `xc`).

Two things follow from the direction of the arrows:

- **Entry points are leaves.** Nothing uses `app/`. `trc_api` and `trc_capi`
  are used only by `app/`, `opt` and tests; `opt` is used only by `app/` and
  by an out-of-tree caller linking `terco::terco_opt` directly (never by
  anything else in `src/`) -- it is built as its OWN shared library
  (`libterco_opt`, never part of `libterco`) for a licensing reason, not a
  layering one: see `src/opt/CMakeLists.txt`.
- **`scf` orchestrates.** It may use everything below it. Nothing below it
  may call back into it.

### Known exceptions

None. The last upward edges were removed by moving `trc_basis_t` and
`trc_pairlist_t` to `core/trc_basis_types`, the device binding to
`core/trc_device_binding`, the 1e/multipole/PCM/DF drivers to
`eri/trc_integrals` (`trc_api` re-exports all of them and is now a thin
facade), the Cartesian-power helper `trc_ao_norm` needed to
`core/trc_cart_perm`, and by passing RI-MP2 the fitted DF factors as an array
instead of the `scf` object that owns them. `python3 tools/archcheck.py --rule
layer-up` reports none.

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
- **No POINTERs**, except where the language forces one: C ABI handles
  (`c_f_pointer` in `src/api/` and `src/opt/trc_opt_capi.F90`) and
  json-fortran's `type(json_value), pointer`. Pointer components and
  optional pointer dummies have been miscompiled by nvfortran, defeat
  no-alias optimisation and hide who owns what. Instead, an owner holds its
  parts by value (or as `allocatable` components when they are optional:
  `allocated()` replaces `associated()`) and passes them down as
  `intent(in)`/`intent(inout)` arguments; a bundle type holds plain values
  only. A procedure whose resource is optional for its callers makes a
  local one when it is absent and calls itself with it (`recursive`, see
  `trc_scf_run`), so its body always works on an argument.

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

`tools/archcheck.py` enforces sections 1 to 4 (and the no-pointer rule of
section 2) deterministically: same tree,
same output. It is regex based like `modmap.py`, so it checks terco's regular
style and does not parse all of Fortran.

| Rule | Flags |
|---|---|
| `layer-up` | a module using a module in a higher layer |
| `module-prefix` | a module without the `trc_` prefix |
| `proc-prefix` | a public procedure not starting with its module's name |
| `short-name` | a procedure, argument, local or module variable containing a word from the `SHORT_NAMES` table (`ctx`, `cfg`, `dens`, ...) |
| `proc-length` | a procedure over 300 lines |
| `no-pointer` | a declaration with the `pointer` attribute, or a `pointer` statement, outside `src/api/`, `src/opt/trc_opt_capi.F90` and `type(json_value)` declarations (section 2) |
| `arg-count` | a procedure with more than 8 dummy arguments (`routine seq` item routines excepted, and thin `do concurrent` kernel drivers: a procedure with a `do concurrent` and at most `KERNEL_HOST_STMTS` (10) executable statements outside its `do concurrent` bodies, not `trc_`-prefixed. Both must take explicit-shape arrays, since a `do concurrent` body may not read derived-type components) |

```sh
python3 tools/archcheck.py                    # fail on new violations only
python3 tools/archcheck.py --summary          # counts per rule
python3 tools/archcheck.py --all --rule short-name
python3 tools/archcheck.py --update-baseline  # after paying debt down
```

Existing violations live in `tools/arch_baseline.json`. CI fails on a
violation that is not in it, and on a length or argument count that grew past
its baseline value. When you fix some, run `--update-baseline` and commit the
smaller file, so the debt cannot return. Files under `src/generated/`, the translated functionals, and any file
whose first 40 lines carry a `Generated by ...` or `do not edit` comment are
exempt from every rule except `layer-up` and `no-pointer`. `SHORT_NAMES`, the limits and the
exempt paths are constants at the top of the script: change them there and
update the table in section 3 in the same PR.

## 7. The code score

`python3 tools/archcheck.py --score` adds up debt points. Lower is better and
zero is clean. It shows the total, points per rule and the worst files, so
you can pick a target and watch the number fall.

| Violation | Points |
|---|---|
| `layer-up` | 10 each |
| `module-prefix` | 5 each |
| `proc-length` | 5, plus 1 per 10 lines over 300 |
| `arg-count` | 2 per argument over 8, at most 30 per procedure |
| `no-pointer` | 5 each |
| `proc-prefix`, `short-name` | 1 each |

Weights are the `WEIGHTS` table and `points()` in the script. A PR that
lowers the score says so in its description (`score: 2693 -> 2540`).

`tools/modmap.py` is the exploration side:

```sh
python3 tools/modmap.py layers     # layer matrix and upward edges
python3 tools/modmap.py calls trc_scf_run --depth 2
python3 tools/modmap.py graph --dir scf > scf.mmd
```
