# The generators

Most of terco's compiled lines are emitted by Python in `scripts/`. This
document is what you need to change them.

---

## Why generate at all

The recurrences behind Gaussian integrals — Obara–Saika vertically, then
Head-Gordon–Pople horizontally — are loops over small, *compile-time known*
index sets. Written as loops they cost table lookups, bounds arithmetic and
branches on every trip. Written out, every index is a literal and the compiler
sees the whole dependency chain.

On a GPU that difference is larger than it is on a CPU, for three reasons that
were each measured rather than assumed:

1. **No branches in the inner loop.** One kernel per `(la,lb,lc,ld)` class
   means no `select case` on angular momentum at all.
2. **Register budgeting per class.** `ptxas` allocates for the code it is
   given. A single kernel covering every class is budgeted for `(dd|dd)`, and
   an `(ss|ss)` quartet then runs at that occupancy. Splitting the classes was
   worth more than any arithmetic change.
3. **Inlining does not cross a module boundary.** NVHPC will not inline an
   `!$acc routine seq` procedure from another module. The unrolled recurrence
   as separate procedures measured **2.4× slower** than the generic loop —
   `v` has to be materialised in memory to cross the call and 35 kB of spill
   appears. Pasted textually into the calling procedure, `v` stays a local.

That last point is why `trc_binkernel.F90` `#include`s a `.inc` of `case`
blocks instead of calling anything. It looks like a mistake and is the opposite.

---

## The five generators

| Script | Emits | Used by |
|---|---|---|
| `gen_vrr.py` | Obara–Saika vertical recurrence | `gen_perclass`, and `--mode include` for the shared kernel |
| `gen_hrr.py` | horizontal transfer plus the digestion, as `case` blocks | the shared kernel |
| `gen_perclass.py` | one four-centre Fock kernel per class, `--split` one module each | the four-centre path |
| `gen_1e.py` | overlap, kinetic, nuclear attraction; `--multipole` for the Cartesian tensors | `trc_1e`, `trc_multipoles` |
| `gen_df.py` | two- and three-centre Coulomb | `trc_df_2c`, `trc_df_3c` |

`gen_df.py` reuses `gen_vrr.py`'s four-centre recurrence with `eta = gamma`,
`Q = C` and `K_cd = 1` — a three-centre integral is a four-centre one with a
degenerate ket, and writing a second recurrence for it would be a second thing
to get wrong.

---

## Shared conventions

Everything below is assumed by more than one generator. Breaking one of them
in isolation produces working code that disagrees with the rest of the library.

**Cartesian ordering is libcint's.** `cart(l)` yields `(nx,ny,nz)` in the order
libcint uses, so a block terco writes drops straight into a caller's matrix
with no permutation. Do not "tidy" it.

**Angular momentum indices are per class, not global.** Inside a class the VRR
workspace is addressed `x = ia + (ic-1)*nca`, where `nca = ncum(lab)` is the
*cumulative Cartesian count for that class* — never the global `NCUM` that
dimensions the array. The two coincide for some classes, which is exactly what
makes getting it wrong survive testing.

**The scratch array is two-buffered.** `v(:, 0:1)` with `cur`/`nxt` swapping as
the auxiliary index `m` walks down. A generator that emits into the wrong
buffer produces a result that is right for `lt = 0` and wrong above.

**Contraction coefficients arrive pre-multiplied** by the primitive norm, and
terco's `basis%build` additionally folds in libcint's per-shell
`common_fac_sp`. A generator never applies either.

---

## Adding a new integral type

Work outward from something that already exists rather than from the algebra.

1. **Find the closest existing generator.** Separable and multiplicative, like
   a multipole? Extend `gen_1e.py`. Coulomb-like with an auxiliary index?
   Start from `gen_df.py`, which is already a reduced `gen_vrr.py`.
2. **Emit one class first**, the simplest non-trivial one — `(p|s)` or
   `(p|p)`. Read the emitted Fortran. If you cannot check it by eye at that
   size, the generator is too clever.
3. **Give it its own module** unless you have a reason not to. See the
   compile-time trap below.
4. **Write the dispatcher last**, and key it on a formula that cannot drift —
   see the stale-include trap.
5. **Wire it into `trc_api.F90`**, not into a caller. The API owns the
   device-resident basis and pair list; a kernel that reaches around it will
   work on the host and fail on the device.
6. **Validate against something outside terco** before believing anything. See
   *Oracle discipline*.

### Raising the angular momentum ceiling

For the **four-centre** kernels, `TERCO_LMAX` is the switch, and the cost is
`L⁴`. Before raising it past 2, check `BOYS_MMAX` in `trc_boys.F90`: it caps
`la+lb+lc+ld` and `(ff|ff)` needs 12 where the table currently stops at 8. The
generated volume is the real obstacle, not the mathematics.

For **one-electron, multipole and density-fitting** kernels, momentum sits on
two or three centres, so they grow as `L²` or `L³` and f or g is affordable.
They are generated at `--lmax 2` only because that is what the shipped set was
built with. The two ceilings are deliberately independent — a fitting basis
normally outranks the orbital basis, which is why `gen_df.py` defaults the
auxiliary ceiling to one above the orbital one (`--lmax 2` gives f on the
auxiliary shell, the shipped set) unless `--lmax-aux` says otherwise.

---

## Traps

Each of these cost real time. They are listed in the order they are likely to
bite.

### A stale generated file is a wrong number, not a failure

`gen_hrr.py` emits a dispatcher keyed on

```fortran
select case (((la*(LMAX + 1) + lb)*(LMAX + 1) + lc)*(LMAX + 1) + ld)
```

The key is built from `TRC_LMAX` **at compile time**; the `case` labels are
fixed **when the file is generated**. A file generated for `lmax 1` and
included in an `LMAX=2` build sends `(pp|pp)` to label 40 in a file whose
largest label is 15. There is no `case default`, so it falls through and
contributes nothing — while the classes where base 2 and base 3 coincide keep
working. Partly right is far worse than zero; a zero would have been noticed
immediately.

That happened, and it survived for a month because under `TRC_PERCLASS` only
the non-symmetric Fock build reaches the shared kernel: every other check
passes without executing a line of it.

**Both generators now emit a guard:**

```fortran
#if TRC_LMAX != 1
#error "This .inc was generated for a different TRC_LMAX. Regenerate it."
#endif
```

**If you write a new generator whose output depends on a compile-time
constant, emit the same guard.** It converts an entire class of silent wrongness
into a compile error, and it is four lines.

### nvfortran treats a rank mismatch as a warning

```fortran
real(dp) :: jmat(ndens, nao, nao)
jmat(mu, nu) = ...          ! two subscripts. Compiles. Runs. Wrong.
```

This is `NVFORTRAN-W-0155`, a **warning**, in a build log with a hundred files
in it. It has bitten twice, on different arrays, and the second time only
surfaced because an unrelated change caused that block to be compiled for the
first time. When adding an index to an array — as batching did with `ndens` —
grep for every reference rather than trusting the compiler to find them.

### The compile-time wall is real and arrives suddenly

Two instances, both fixed by splitting rather than by optimising:

- **39 multipole components fused into the S/T/V kernel.** One
  `!$acc routine seq` body per class; NVVM ground for over **thirty minutes** at
  LMAX=2 without finishing. Split into its own module: **2m27s**. Both
  quantities are computed once per geometry, so sharing the 1-D tables saved
  nothing and cost everything.
- **The unrolled recurrence in the shared kernel at LMAX=2.** That kernel holds
  two copies, so unrolling it puts twice the entire per-class kernel set into
  one translation unit. Still compiling after twenty-five minutes. The unroll
  is now emitted only for `LMAX<=1`.

The rule that falls out: **one module per class, and never fuse two quantities
into one device routine because they share intermediates.** Sharing arithmetic
that is not hot is not an optimisation.

### Gill's generation-step sieve

`gen_vrr.py` implements the sieve, and it is what makes d functions compile at
all — `(dd|dd)` drops from 3683 VRR statements to 2179.

It is **not** dead-code elimination, and implementing it as such gives exactly
zero, because the compiler already does that. The sieve keeps the complete
auxiliary list and removes any node **every one of whose consumers can be
formed from something else still on the list**, then iterates to a fixed point.
The removed nodes are ones that were *reachable and useful* — they are simply
not the cheapest route to anything.

If you extend the sieve, the fixed-point iteration is the part to be careful
with: removing a node can make another removable, and stopping after one pass
leaves most of the win on the table.

### The HRR is a sieved DAG too

`gen_vrr.py`'s sieve builds the VRR as a DAG of shared intermediates reachable
from the targets a class actually needs, kept to the cheapest route to each
one. `gen_hrr.py`'s HRR did not get the same treatment until now: `emit_class`
still summed each of a class's na*nb*nc*nd output components as its own
independent closed-form expression over the contracted VRR buffer `g`, so two
outputs that share a partial bra or ket transfer recomputed it from scratch.

`emit_class_dag` builds the HRR the same way the VRR already is built: from
the targets back to the VRR buffer. Bra transfer, `(a,b| = (a+1_k, b-1_k|
+ AB_k*(a, b-1_k|)`, and ket transfer, the same on `(c,d)` with `CD` once `b`
is exhausted, each recurse one elementary step at a time; a node is keyed on
the literal `(a,b,c,d)` tuple reached and memoised on that key, so whichever
output asks for `(a, b-1_k|c,d)` first builds it, and the output that would
otherwise recompute it hits the memo instead. The recursion direction (which
axis of `b` or `d` to reduce) is a fixed, deterministic rule — the first
nonzero axis in `x,y,z` order — rather than anything cost-modelled, which is
enough to let memoisation actually find the sharing: a value-dependent choice
would route two paths to equal content by unequal names and never collide in
the memo. Every output's own node is reachable this way, so nothing needs a
separate elimination pass — the DAG is sieved by construction, not pruned
afterward.

**Emission order matters more than it looks.** Building every node before
writing any output leaves a node live from its creation to the LAST output
that reads it — for `(dd|dd)` that meant 1406 values simultaneously live
tracked by hand. `emit_class_dag` instead walks the targets in their existing
output order (which already groups every `(a,b)` pair under a fixed `(c,d)`,
the DAG's own sharing axis) and emits each output's still-missing nodes right
before that output, so most nodes die within the block that created them.
Measured peak live count, same two schedules, same classes:

| class | nodes | peak live, all-nodes-first | peak live, interleaved by target |
|---|---|---|---|
| (pp\|pp) | 162 | 83 | 6 |
| (dp\|dp) | 612 | 326 | 12 |
| (dd\|dp) | 1980 | 703 | 31 |
| (dd\|dd) | 5293 | 1406 | 637 |

Interleaved is what ships. `(dd|dd)`'s 637 is still large by hand-count, and
is exactly the number this rewrite worried would turn into GPU spill — see
below for what `ptxas` actually did with it.

**Operation counts** (the `+ - * / ` characters on the right-hand side of
every assignment, closed form against the DAG, per contracted quartet):

| class | closed-form ops | DAG ops | ratio |
|---|---|---|---|
| (pp\|pp) | 567 | 324 | 0.57 |
| (dp\|dp) | 2268 | 1224 | 0.54 |
| (dd\|dp) | 11340 | 3960 | 0.35 |
| (dd\|dd) | 50544 | 10586 | 0.21 |

Fewer than a fifth of the arithmetic for `(dd|dd)`, and the regenerated
`src/generated` is smaller to match: 81 files, 49152 lines inserted against
88020 deleted net across the set, because each output is now a short `h(k) =
ref + coef*ref` chain instead of a long wrapped sum.

**Measured: spill went down, the stack frame went up, and compile time got
worse.** `ptxinfo` has to be read at the right entry. For the large classes
the `do concurrent` kernel (`trc_pc_kABCD_pcsABCD_79_gpu`, 64 registers for
`(dd|dd)`) is only a wrapper: it CALLS the item routine `pcsiABCD`, which is
not inlined at that size, and that routine is where the work and the local
memory are. `(dd|dd)`, closed form against DAG, same flags:

|                                   | closed form | DAG     |
|---|---|---|
| `pcsi2222` stack frame            | 59 KB       | 95 KB   |
| `pcsi2222` spill stores / loads   | 52 / 75 KB  | 38 / 58 KB |
| compile, alone on an idle machine | 485 s       | 639 s   |

`(pp|pp)` is inlined, and its kernel sits at 255 registers with about 2 KB
of spill in both forms. So the DAG cuts spill traffic by about a quarter on
`(dd|dd)` at the price of a larger frame, and nvfortran spends 32% LONGER
compiling it, despite a fifth of the arithmetic: it now has thousands of
named intermediates to schedule rather than independent expressions. (An
earlier reading of 58 and 64 registers with no spill came from the wrapper
entries, and does not describe the work.) It also shows up, modestly, in the Fock build itself: on
172 waters, RHF/6-31G\*, SAD guess, `eri_thresh` 1e-10 (`CUDA_VISIBLE_DEVICES`
pinned to one idle V100):

| basis | Fock s/iteration, old closed form | Fock s/iteration, DAG |
|---|---|---|
| 6-31G\* (has d) | 2.26 | 2.236 |
| 6-31G (no d) | 0.968 | 0.970 |

A small win where there are d functions to share, a tie within noise where
there are not — consistent with identical occupancy and a real but modest cut
in instructions actually executed. `E = -11399.3977137361` on 6-31G\*, bit
for bit the same value the closed form gives.

**Numeric correctness is checked in pure Python**, independent of Fortran,
GPU or a basis set: `scripts/check_hrr_dag_numeric.py` evaluates
`emit_class_dag`'s recursion and `emit_class`'s closed form (kept only as
this check's oracle; no generator's command line reaches it any more) on
random VRR-buffer and AB/CD values for every class up to `TRC_LMAX=2`, and
separately re-derives every warp-role split's own compacted closure
(`role_plan`'s `hset`/`hmap`) and checks it against the class's own
unsplit DAG — the role split's closure and renumbering logic is exactly
where a bug would actually live, since `emit_class_dag`'s own recursion says
nothing about whether a role's SLICE of it was taken correctly. Both run as
`ctest`'s `terco/check_hrr_dag_numeric`.

On these numbers the DAG is not yet the right default for nvfortran: a 1%
Fock-time gain does not pay for a 32% longer compile of the class that
already dominates the build. On this branch the closed form survives only as
the numeric check's oracle. Before this merges, it returns as the nvfortran
default, and the DAG stays as a generator option, chiefly for a CUDA backend,
where nvcc's compile cost on large unrolled code is far lower.

---

## Porting the directives to OpenMP

There are two scripts and they run in a fixed order:

```sh
python3 tools/dc_to_omp.py  --target cpu --write src
python3 tools/acc_to_omp.py --target cpu --translate-atomics --write src
```

1. **`dc_to_omp.py`** rewrites `do concurrent` as OpenMP worksharing, maps
   `local(...)` to `private(...)`, and **strips `pure`** from every procedure
   that ends up holding a worksharing directive or an atomic — cascading to
   their pure callers, since a pure procedure may only call pure ones.
2. **`acc_to_omp.py`** then translates the data directives and, with
   `--translate-atomics`, the atomics. That flag is only legal *after* step 1,
   because OpenMP forbids an atomic in a `PURE` procedure and OpenACC does not.

Reversing the two produces source that does not compile. CI runs both, builds
the result, and validates it **on four threads**.

### Three traps this port walked into, all of them silent

**A `do concurrent` that does not name what it assigns becomes a race.** The
standard privatises a scalar assigned before it is read, so omitting
`local(...)` is correct as written — and the moment the construct becomes
`!$omp parallel do`, whose default is *shared*, that scalar is one variable
torn between every thread. `fock_all_nosym` had exactly one such loop and it
was wrong after conversion while everything else passed.
`tools/dc_locality_lint.py` now refuses one, and CI runs it.

**A directive behind a preprocessor macro is invisible to a line-based scan.**

```fortran
#define FG_ATOMIC !$acc atomic update
```

does not begin with `!$acc`, so the converter never saw it, and every use of
the macro expanded to OpenACC in a tree that was otherwise OpenMP. Under a
compiler without `-acc` that is an inert comment — for an *atomic*, a silent
race. `acc_to_omp.py` now transforms the replacement text of a `#define` too.

**A race is invisible on one thread.** Both of the above passed at
`OMP_NUM_THREADS=1` and failed at 8, with the enumerated digestion disagreeing
with the folded one by a factor of 0.94 on a symmetric density — where the two
are the same operator by construction. That is why the CI job sets a thread
count rather than trusting the runner's default, and why `check_nosym` exists
at all.

### What this port is, and is not

`--target cpu` gives host multicore. A genuine **offload** port would use
`--target gpu`, which emits `!$omp target teams distribute parallel do` — the
same scripts, and untested here because validating it needs a machine with
both a GPU and an OpenMP-offload toolchain. The pieces are in place: the
purity that blocked the atomics is removed by step 1, and the device routines
carry `!$omp declare target` from `!$acc routine seq`.

## Oracle discipline

A new comparison harness is untested code like any other.

**Validate the oracle on a quantity you already trust before letting it convict
anything.** The multipole checker dumps the *overlap* alongside the multipoles
and refuses to compare a single component until that agrees. That control was
written last and should have been written first: two days went into a
disagreement that was in the reference system, not the kernel — it had a
contracted shell, and pyscf renormalises contracted shells to unit self-overlap
where terco's reference did not. Two different basis sets, being compared.

Practical consequences for a new check:

- **Use uncontracted shells** when comparing against another program, unless
  you have verified both normalise contractions identically. One primitive has
  nothing to renormalise.
- **Put the centres off-axis and make them unequal.** A symmetric test system
  makes half the possible component-ordering errors invisible: with everything
  on `z`, the `x` and `y` moments coincide.
- **Prefer a control that must hold by construction** over a recorded number.
  The folded and enumerated digestions must agree on a symmetric density; that
  is a stronger statement than any tolerance against a stored file, and it is
  what found the stale-include bug.

When two outputs that are *generated from the same expression* disagree — `xy`
and `yx` of a Cartesian tensor, say — stop reading the kernel. Identical code
cannot produce different numbers. The fault is in the addressing, and that
realisation was available on day one of a hunt that took considerably longer.

---

## Regenerating

```sh
# in-tree, refreshing the committed LMAX=2 set
cmake -S . -B build -DCMAKE_Fortran_COMPILER=nvfortran -DTERCO_REGENERATE=ON

# a different ceiling, into the build tree only
cmake -S . -B build-l1 -DCMAKE_Fortran_COMPILER=nvfortran -DTERCO_LMAX=1
```

Or call them directly, which is what CMake does:

```sh
python3 scripts/gen_perclass.py --lmax 2 --split src/generated
python3 scripts/gen_1e.py --lmax 2 -o src/generated/trc_1e_kernels.F90
python3 scripts/gen_1e.py --lmax 2 --multipole src/generated/trc_mult_kernels.F90
python3 scripts/gen_df.py --lmax 2 -o src/generated/trc_df_kernels.F90
```

`[1]` in the sources refers to the reference in the README.

**After regenerating, run `ctest`.** The generators have no tests of their own;
the emitted code is what is tested, and it is tested against libfint and pyscf
rather than against a stored copy of its own output.
