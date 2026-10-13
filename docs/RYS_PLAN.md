# Implementing Rys quadrature in `do concurrent`

Read `RYS.md` first for why.

STATUS: implemented, behind `TERCO_RYS` (default OFF) and runtime-selected
per class by `TRC_RYS_CLASSES` (`src/trc_rys_select.F90`), never the default.
The rest of this document is the plan as written before that work started;
it is kept because the reasoning in it is still the reasoning behind the
result, not because every phase below ended up matching the final shape --
`scripts/gen_rys.py`'s own docstring records where the splice points moved
once `--split` and `trc_decontract` existed. The measurements the plan asked
for (gates A-F and the acceptance criteria in the commit that reports them)
are what settled it, not this document.

## What we are copying, and what we are not

**gamess-libERI already does this** (`$GAMESS_LIBERI/rhf/rys/int*_rysgen.F90`):
per-class generated kernels, `do concurrent(iquart=...)` with `local(...)`,
`!$acc atomic update` on the digestion. The architecture is proven; the
algebra in those files is the reference for the recurrences.

What it does that we should NOT, and why terco can do better -- its per-quartet
state at (fd|dd) is ~4,100 doubles (~32 kB), barely under HGP:

| its choice | cost | ours |
|---|---|---|
| `xin/yin/zin(nroots*g_size)` | 5x the 2D array | one root at a time |
| `eri_value(2160)` materialised | the whole output block | fold into J/K as produced |
| `rgrid,wgrid,p0,p1,p2(35)` Stieltjes scratch | 175 doubles + divergent iteration | tabulated roots |

The middle row is the batched-shell-quartet advantage and it is already built:
terco's six-update folded digestion consumes a quartet's contribution without
ever forming the block.

## Target working set, (dd\|dd)

| | doubles |
|---|---|
| 2D integrals, ONE root at a time: `3 * g_size` | 243 |
| roots + weights | 10 |
| digestion accumulators (unchanged from today) | ~234 |
| **total** | **~490 (~4 kB)** |
| today, HGP | 4,971 (~39 kB) |

**~8x, not the 19x in RYS.md.** That figure compared intermediates only; the
accumulators do not shrink, and honesty about that is the difference between a
plan and a pitch.

## Phases, each with a gate that can stop the project

### 1. Roots on the device (the only real risk)

**Tabulate, do not construct.** GAMESS builds the roots per quartet: a fixed
35-point grid, a three-term recurrence for `alpha`/`beta`, then a 5x5 symmetric
tridiagonal QL with implicit shift. Correct, dp-only, and proven -- but it
carries `rgrid,wgrid,p0,p1,p2(35)` = 175 doubles of private scratch and the QL
iterates to a tolerance, which is warp divergence. gpu4pyscf instead
interpolates from a table (`rys_roots_dat.cu`, degree 13 per interval).

Take the table. It removes the scratch, the iteration and the divergence in one
move, and leaves only the n roots and weights in private state.

`scripts/gen_rys_roots.py` computes roots/weights at high precision, fits
Chebyshev coefficients per interval per `nroots` = 1..5, and emits Fortran.
Same shape as the Boys table already in `trc_boys`.

TWO CONSTRAINTS THAT DECIDE THE IMPLEMENTATION:

- **A generator may use mpmath/sympy, but CI must be told.** The existing
  generators import only the standard library and the `generators` job runs
  them on a bare `actions/setup-python` with no `pip install`, then requires
  no diff. The moment a generator imports mpmath that job goes red, so the
  install step lands in the SAME commit as the generator, never after.
  Pin the versions: this job's whole purpose is that regenerating reproduces
  the committed files byte for byte, and an unpinned library that changes its
  last digit breaks that silently.

  With mpmath available the roots are easy -- high-precision quadrature and
  root finding are built in, and the ill-conditioned arithmetic happens ONCE
  at generation time, on the host, which is where it belongs.
- **Not a `parameter` array.** `trc_boys` already records why: nvfortran's
  fort2 segfaults on parameter arrays of this size. A plain module array with
  an initialiser.

Compile time will be unpleasant -- gpu4pyscf's equivalent is ~50k lines,
though we need only `nroots <= 5` for LMAX=2 where they cover higher. That is
a known cost of this project's approach and `docs/GENERATORS.md` already names
the compile-time walls.

- **Gate A:** roots match a trusted reference to 1e-14 for n = 1..5 over the T
  range the LMAX=2 classes reach. libfint's validated implementation is the
  oracle; check the table against it before anything uses it.
- **Gate B:** table size and compile time are tolerable. If the table has to
  grow so far that a rebuild becomes impractical, reconsider constructing the
  roots after all.

### 2. One kernel by hand, before any generator

Write `(dd|dd)` by hand against the gamess-libERI reference. One class, no
codegen, wired into `pc_dispatch` behind a flag so both paths coexist.

- **Gate C:** agrees with the existing HGP kernel to ~1e-13 on the check_selftest
  systems. Class-by-class comparison against a validated implementation is the
  strongest harness this project has -- use it before generalising.
- **Gate D:** register count and occupancy from `cuobjdump -res-usage`. If it
  does not come off the 255 ceiling, the working-set argument was wrong and
  everything after this is wasted.
- **Gate E:** it must be FASTER on gly30. Rys does more arithmetic; if the
  kernels were compute-bound rather than capacity-bound, this is where that
  shows up as a regression.

### 3. `scripts/gen_rys.py`

Only after gate E. Emits the 2D recurrence and the root loop; reuses the
existing PROLOGUE (primitive-pair loops, screening) and digestion verbatim.
Simpler than `gen_vrr.py` + `gen_hrr.py` because the recurrence is uniform in
one index.

- **Gate F:** the generators job in CI regenerates with no diff, same as the
  HGP ones.

### 4. Dispatch by class

Rys above the crossover, HGP below. `RYS.md` puts it near (pp\|pp) -- Rys is
0.4x at (ss\|ss) and 6.3x at (pp\|pp) -- but MEASURE it per class rather than
trusting the table. `pc_dispatch` already selects per class; this is a keyed
choice, not new machinery.

## What this does not fix

`do concurrent` flattens the iteration space (`collapse(2)`, no `threadidx%y`,
`SHARED:0` across all 299 kernels). Threads still cannot cooperate on one
quartet, so if gpu4pyscf's remaining advantage is the decomposition rather than
the algebra, Rys will narrow the gap and not close it. That is the honest
expectation to set before starting.

## Order of work

1 -> 2 -> stop and look at the numbers. Phases 3 and 4 are mechanical once
gate E passes, and worthless if it does not.
