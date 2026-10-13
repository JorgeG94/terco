# Rys quadrature for the four-centre ERIs: what it would buy

Written as exploration, before any of it was implemented; see
`docs/RYS_PLAN.md` for where the work that followed ended up and what it
measured. Kept as the working-set argument for WHY Rys was worth trying --
that argument did not change, only the kernel template it had to be spliced
into.

## The number that matters

HGP materialises an intermediate whose size is the product of the Cartesian
component counts. Rys materialises 2D integrals indexed by ANGULAR MOMENTUM,
so its intermediate is `(la+1)(lb+1)(lc+1)(ld+1)` per axis. At (dd|dd) that is
81 against 1225.

Per-quartet working set, terco's committed kernels against Rys sized the way
gpu4pyscf sizes it (`nroots*2 + g_size*3 + 6`):

| class    | terco HGP | Rys | ratio |
|----------|-----------|-----|-------|
| (ss\|ss) |         4 |  11 | 0.4x  |
| (ps\|ss) |        15 |  14 | 1.1x  |
| (pp\|pp) |       381 |  60 | 6.3x  |
| (dp\|pd) |     1,524 | 122 | 12.5x |
| (dd\|pd) |     2,748 | 176 | 15.6x |
| (dd\|dd) |     4,971 | 259 | 19.2x |

At (dd|dd): **39 kB per quartet becomes 2 kB**. A thread has 255 registers,
which is 1 kB. That is why every non-trivial per-class kernel sits at the 255
register ceiling with ~12% occupancy.

Rys is WORSE below about (pp|pp), which is why gpu4pyscf tries
`rys_jk_unrolled()` first and falls back to general Rys above. Any terco
version wants the same split, and `pc_dispatch` is already the place to put it.

## What is already in hand

libfint has the Rys roots in Fortran, validated against the C bit for bit:
`cint_rys_roots.f90`, `cint_rys_root_n.f90`, `cint_find_roots.f90`,
`cint_tab_roots_x0.f90`. No allocate, no I/O, no stop, no save, no goto. That
is normally the hardest part of a Rys implementation and it is done.

terco's scaffolding is scheme-agnostic: the primitive-pair loops, the
screening, and the folded six-atomic digestion all take a Cartesian block and
do not care how it was produced. Rys replaces the VRR+HRR middle only.

## What is not in hand

1. **Quad precision cannot go to the device.** `cint_const.f90` has
   `qp = real128`; nvfortran has no 128-bit real and using it as a kind is a
   compile error, not a fallback. So `cint_schmidt_qp` and `cint_wheeler_qp`
   are unavailable and only the `_DP` branches can be taken. At LMAX=2
   (order <= 8, nroots <= 5) double is expected to suffice -- that needs
   checking against the qp path on the host, not assuming.

2. **Thread divergence in root finding.** `cint_find_roots.f90` iterates
   `do while (abs(x1 - x0) > x1*accrt)` with early exits. Data-dependent trip
   counts mean a warp runs at its slowest lane. gpu4pyscf sidesteps this with
   tabulated roots (`rys_roots_dat.cu`); that table is the fallback if the
   iterative version diverges badly.

3. **Assumed-shape dummies** (`u(0:)`, `w(0:)`) want explicit shape before
   they go anywhere near a device kernel.

## How far `do concurrent` goes, measured

nvfortran 26.5, `-stdpar=gpu -gpu=cc70`, read from `-Minfo=accel` and
`cuobjdump -res-usage`:

| construct                          | what the compiler does                     |
|------------------------------------|--------------------------------------------|
| rank-1 `do concurrent`             | `blockidx%x threadidx%x`, 128 threads      |
| rank-2 `do concurrent(j=..,i=..)`  | **`collapse(2)` -- auto-collapsed to 1D**  |
| rank-2 with `local()`              | same, collapsed                            |
| outer `dc` + inner ordinary `do`   | only the outer is parallelised             |

The iteration space is ALWAYS flattened. There is no `threadidx%y`, so there is
no way to make several threads cooperate on one quartet, and no way to give
them shared state. Confirmed on the shipped library: 299 kernels in
`libterco.so`, `SHARED:0` in every one. Not one byte of shared memory is
reachable from standard Fortran.

**This is the constraint, and it is what makes Rys the right answer rather than
a workaround.** If every iteration must own its working set privately, the only
lever is to make that working set smaller -- which is exactly what the table
above measures. Rys is the ERI scheme with the smallest private state, so it is
the scheme that fits the model this library is committed to.

## What gpu4pyscf does that terco cannot say in standard Fortran

Its kernel is a 2D block: `threadIdx.x` indexes quartets, `threadIdx.y` indexes
output components, so `gout_stride` threads COOPERATE on one quartet and the
intermediates live in dynamic shared memory -- over 48 kB per block, opted into
with `cudaFuncSetAttribute`. Block shape and shared-memory size are chosen at
runtime per class by `threads_scheme_for_k`.

`do concurrent` has no way to say "this array is shared across the block", so
that decomposition is out of reach without directives.

**But the table above does not depend on it.** Rys shrinks the per-quartet
working set by 6-19x whatever the thread mapping is, and 2 kB has a chance of
fitting in registers where 39 kB never did. The block-cooperative question is
separate and can stay unanswered.

## Caveat on the comparison

`root_g_cache_size` is what gpu4pyscf reserves in shared memory; the gout
accumulators live in registers per thread and are not counted in it. terco's
figure includes its `vbuf` HRR buffer. So the ratios above flatter Rys
somewhat -- the direction is not in doubt, the exact factor is.

## The measurement to make before building any of this

The kernels are compute bound with zero spill and a saturated GPU (see the
profiles in `terco-decks/`). Rys shrinks the working set but does MORE
arithmetic. If terco is already near peak FP64 then Rys is slower, and the
19x means nothing. One `ncu --section SpeedOfLight` on `pc2112` decides it.
