# Architecture: what the library hands a caller

The kernels work and are validated. What is missing is the thing that makes
them a *library*: a caller should hand over a basis and get work done, without
knowing that per-class kernels, bins, sieves or Boys tables exist.

This document fixes the object model before the refactor, so the pieces that
already exist can be moved into it rather than re-invented.

## The shape, and where it comes from

cuEST (`/home/jorge/dev/CUDALibrarySamples/cuEST`) solves exactly this problem
and its C examples read as a specification. Its one-electron flow is:

    cuestCreate(handle)                      context
    cuestAOShellCreate(...)      per shell   raw exponents, coefficients, centre
    cuestAOBasisCreate(...)                  basis from shells
    cuestAOPairListCreate(..., 1.0e-14, ...) SHELL PAIRS, with the threshold
    cuestOEIntPlanCreate(basis, pair_list)   a plan, like a cuFFT plan
    cuestOverlapCompute(...)                 and Kinetic, Nuclear, ...

Two things in that are worth copying outright:

1. **The pair list is a first-class object, built once, shared by every
   operator, and it takes the screening threshold at construction.** Screening
   is a property of the pair list, not of each compute call.

2. **There is a separate plan stage** between the pair list and the compute.
   That is where per-operator scheduling lives.

What we will NOT copy is the workspace-query protocol
(`*CreateWorkspaceQuery` returning descriptors, caller allocates persistent
and temporary buffers, temporary freed after construction). It exists because
cuEST is a C library that refuses to own device memory. We are Fortran with
allocatable components and a `release` type-bound; owning our own memory is
simpler for the caller and costs nothing. If a C ABI later needs the
query protocol it can be added at that boundary, not through the core.

## Objects

| object | holds | built from | lifetime |
|---|---|---|---|
| `trc_context_t` | device selection, stream, defaults | once per process | process |
| `trc_basis_t` | shells: l, exponents, coefficients, centres, AO offsets | RAW ARRAYS from the caller | per geometry |
| `trc_pairlist_t` | screened, binned shell pairs + primitive-pair data | basis + threshold | per geometry |
| `trc_plan_t` | per-operator launch schedule | basis + pairlist + operator | per operator |

Compute entry points take `(context, plan, pairlist, inputs) -> outputs` and
own no state.

## The caller supplies raw arrays, not a basis set

> "the consumer will just pass exponents, coefficients and whatever we need for
> a basis (we won't provide a basis)"

So the primary constructor takes flat arrays, one entry per shell:

    call basis%build(nshell, sh_l, sh_nprim, sh_exp, sh_coef, sh_centre, &
                     natm, at_z, at_r)

and NOT libcint's `atm`/`bas`/`env`. The current `from_libcint` was expedient
for testing against libfint and stays as a convenience adapter, but it is the
wrong primary interface: it forces every caller to build libcint's packed
integer layout, which only a caller that already links libcint can do.

Normalisation is the library's problem, not the caller's. `common_fac_sp` is
folded into the coefficients at import; whether the caller's coefficients
already carry the primitive norm is a documented flag on the constructor, not
a silent assumption. Getting this wrong is the single most common way to be
wrong by a smooth factor and look plausible -- it has already happened twice
in this repo.

## The pair list generalises what the four-centre path already has

`pair_bins_t` in `src/trc_bins.F90` IS an AO pair list: shell pairs binned
by `(la, lb, K, size class)` with a Schwarz pre-screen. It was written for the
Fock build and knows nothing about one-electron or density-fitting work. The
refactor is to promote it, not rewrite it:

  - keep the binning by `(type, size class)`, which is what keeps the kernel
    free of conditionals;
  - keep the Schwarz bound and the `qq*qmax` pre-screen;
  - add the primitive-pair data (`zeta`, `P`, `PA`, `PB`, `K_ab`) that
    `build_pairs_hgp` currently builds separately, so it is computed ONCE per
    shell pair and reused by every operator;
  - let the screening bound depend on the operator, since Schwarz
    `sqrt((ab|ab))` is the right bound for two-electron work and the wrong one
    for overlap.

That last point is the substantive design question. Overlap decays as
`exp(-xi |AB|^2)` and needs no two-electron quantity at all; nuclear
attraction needs the same overlap-scale bound; three-centre Coulomb needs
`Q_ab` on the bra and a one-centre bound on the auxiliary. So the pair list
should carry BOTH a Gaussian-product bound and a Schwarz bound, and each plan
picks the one its operator needs.

## Screening the one-electron and DF paths

Currently absent, and the cost is structural rather than constant:

  - **one-electron**: the drivers walk the full `nshell^2` square. Pairs whose
    Gaussian product is negligible produce a zero block and are computed
    anyway. For an extended molecule that is most of them.
  - **three-centre**: one thread per `(i, j, k)` with the bra pair product
    recomputed for EVERY auxiliary shell -- `nshell_aux` times redundant,
    which for a real auxiliary basis is a factor of several hundred. This is
    the single worst inefficiency currently in the library, and the pair list
    fixes it by construction.

## What "ready" means

> "the library won't be ready until we can basically just hand over work"

Concretely, a caller should be able to write an SCF iteration without touching
anything below the API:

    call ctx%init()
    call basis%build(...raw arrays...)
    call pairs%build(basis, thresh)
    call plan_1e%build(basis, pairs, TRC_OP_CORE)
    call plan_jk%build(basis, pairs, TRC_OP_ERI_FOCK)

    call trc_core_hamiltonian(ctx, plan_1e, pairs, hcore, smat)
    do while (.not. converged)
       call trc_fock(ctx, plan_jk, pairs, dmat, fock)     ! G = 2J - K
       ! diagonalisation, DIIS, density -- the caller's problem
    end do

DIIS, diagonalisation, convergence and the density update stay with the
caller. The library's contract is: given a density, return a Fock matrix; given
a basis, return the one-electron matrices; given a basis and an auxiliary
basis, return the density-fitting tensors.

## Order of work

1. ~~Raw-array basis constructor; demote `from_libcint` to an adapter.~~ DONE.
2. ~~Generalise `pair_bins_t` into `trc_pairlist_t`, carrying primitive-pair
   data and both bounds.~~ DONE.
3. ~~Rewrite the one-electron and DF drivers over the pair list.~~ DONE.
4. ~~Promote the segment/launch construction into an explicit
   `trc_plan_t`.~~ DONE, and it buys nothing measurable -- see below.
5. ~~Timing at realistic size for the 1e and DF paths.~~ DONE.
6. ~~A C ABI.~~ DONE, without the workspace-query protocol.

## What the measurements said

First timings for these paths, w64 (192 atoms, 832 orbital / 4480 auxiliary
functions), one V100, pair threshold 1e-12:

| | |
|---|---|
| pair list build (host) | 0.050 s |
| pairs kept | 39.7% of 166176 |
| one-electron S, T and V | 0.026 s |
| three-centre (mn\|P) | 0.924 s |

Three corrections to what was written above before anything was measured.

**Screening is the win, not the recomputation.** The old three-centre driver
walked `nshell^2 x nsx = 4.5e8` kernel calls; the pair-list version walks
`npair x nsx = 8.9e7`, about 5x fewer, from screening plus the triangle.
Precomputing the bra product removes per-call SETUP -- notably one `exp` per
primitive pair per call -- not redundant integral evaluation, since the VRR
count is unchanged. Calling it "the worst inefficiency in the library" was
overstated.

**The plan object buys nothing.** Rebuilding the launch plan costs 0.14 ms
against a 3.02 s Fock build at w128 -- 0.005%. It is worth having because the
API needs an object with that lifetime and because it is where per-operator
scheduling will go, but not for speed. The double loop it contains is over
live BIN pairs, of which there are thousands, not over shell pairs.

**One-electron work is already cheap.** 26 ms at 832 functions, against 924 ms
for the three-index tensor over the same pair list. Effort spent on the 1e
path from here is effort spent on 3% of the problem.

## What is still missing

- **Density fitting needs more than the tensor.** `(P|Q)` and `(mn|P)` are
  computed; forming and inverting the fitting metric, and contracting to J and
  K, are not here. That is where the makeEFP path actually spends its time.
- **The four-centre path and the new API are separate worlds.** `fock_bins`
  takes `pair_bins_t` and the `hp_*` arrays; the API takes `trc_basis_t` and
  `trc_pairlist_t`. A caller wanting both today builds two containers over
  the same geometry. They should be one.
- **No d functions through the C ABI by default.** The kernels support d; the
  API's `TRC_MAXL` is a compile-time parameter and the shipped build has to
  pick one.
- **No streams, no multi-GPU.** cuEST's handle carries a stream; ours carries
  nothing, because there is no context object yet.


## Density fitting: the metric is free, the tensor is the problem

### How cuEST does it

The metric never crosses cuEST's API. `cuestDFCoulombCompute(handle, plan,
params, workspace, d_D, d_J)` takes a density and returns J; everything
metric-related lives inside `cuestDFIntPlanCreate(handle, primaryBasis,
auxiliaryBasis, pair_list, ...)`. That is the same lifetime rule as the rest
of that API: the metric depends on the two bases and the geometry, not on the
density, so it is built once with the plan and never mentioned again.

Worth copying. A caller should not be handed `(P|Q)` and told to factor it.

### Never invert it

`(P|Q)` is symmetric positive definite -- it is a Coulomb metric over
auxiliary functions -- so the explicit inverse is both slower and less stable
than the factorisation. On device:

  - **Cholesky**, `cusolverDnDpotrf`, once per geometry.
  - **J needs only solves.** Contract `d_Q = sum_mn (mn|Q) D_mn`, solve
    `(P|Q) c = d` with `cusolverDnDpotrs`, then `J_mn = sum_P (mn|P) c_P`.
    The fitted tensor B is never formed at all.
  - **K needs B**, `B^P_mn = sum_Q (mn|Q) [L^-T]_QP`, one `cublasDtrsm`
    against the Cholesky factor. This is where the cost is.
  - **Fallback.** If `potrf` fails the auxiliary set is near
    linearly dependent; eigendecompose (`cusolverDnDsyevd`, already in
    `mqc_scf_device`) and build `(P|Q)^-1/2` dropping eigenvalues below a
    threshold. Rare, but it is the difference between a diagnostic and a
    crash.

### The sizes that actually decide the design

With naux ~ 3 nao, which is typical for a JK fitting set:

| system | nao | naux | (P\|Q) | potrf | full (mn\|P) | packed + screened |
|---|---|---|---|---|---|---|
| w64 | 832 | 2496 | 50 MB | 5 Gflop | 13.8 GB | 2.8 GB |
| Gly30 | 1273 | 3819 | 117 MB | 19 Gflop | 49.5 GB | 8.7 GB |
| w128 | 1664 | 4992 | 199 MB | 41 Gflop | 110.6 GB | 16.6 GB |

The metric is a rounding error in both memory and flops. The three-index
tensor is 14 to 110 GB, and 3 to 17 GB even packed over the lower triangle
and screened. A 32 GB V100 holds w128 in core only just, and nothing larger.

That is why cuEST's examples are named `core_df_jk`: in-core is a variant, not
the algorithm.

### What follows for terco

**Block over the auxiliary index.** Both J and K accumulate over P, so
processing the auxiliary basis in blocks bounds memory at
`npair x block x 8` bytes regardless of system size:

    for each block of auxiliary functions P:
        build (mn|P) for that block            <- trc_df_3c, already exists
        J:  accumulate d_P, later J += (mn|P) c_P
        K:  B_blk = (mn|P) L^-T ; K += B_blk B_blk^T

Neither J nor K ever needs all of B resident, and the block size becomes a
memory knob rather than a hard limit. This is the piece to build next; the
per-block integrals are done and validated.

Storage within a block should be the packed lower triangle over SCREENED
pairs -- the pair list already has both -- which is the difference between
2.8 GB and 13.8 GB at w64.

## What makeEFP needs from terco, traced through mqc

Read off `backends/libcint/mqc_efp_potential.f90` and what it calls, rather
than guessed. MakeFP is: one RHF, a Boys localization, a set of coupled-
perturbed solves, and the multipole/polarizability analysis on top.

    make_efp_potential
      -> run_libcint_rhf            SCF          -> build_fock_direct
      -> multipole_matrices         1-electron, through octupole
      -> boys_localize              linear algebra, no integrals
      -> distributed_polarizability -> mqc_libcint_cphf
      -> distributed_dynamic_cross  -> mqc_libcint_cphf

### The one that decides everything: batched densities

mqc's own note on `build_fock_direct_many` is the whole story:

> the coupled-perturbed equations for the dynamic polarizabilities need
> roughly a hundred right-hand sides -- nine perturbations times twelve
> imaginary frequencies -- and the matvec for each is a Fock build on a
> different response density, so without this the frequency loop pays for the
> integrals a hundred times over

So makeFP is not one Fock build, it is O(100). Their measured win from
batching is 3.9x at 6 densities, 4.3x at 12, 3.8x at 24 -- saturating around
fourfold.

The same ceiling applies to terco and can be predicted from a number already
measured here. Cost of a batched build is `eval + N*digest`, and terco's
digestion is 24% (Gly30) to 37% (RNA3) of the kernel, so with 0.30:

| N | speedup vs N separate builds |
|---|---|
| 4 | 2.1x |
| 12 | 2.8x |
| 24 | 3.0x |
| 100 | 3.3x |

The ceiling is `1/digestion fraction` = 3.3x, which is why mqc sees saturation
at about four: their digestion is a slightly smaller share. Batching cannot do
better than removing the integral evaluation from every density but the first.

**The GPU constraint that mqc does not have.** terco's digestion accumulates
six output blocks in registers. Per extra density that is 54 more doubles at
(pp|pp) and 216 at (dd|dd), on a kernel already at 255 registers with 2.5 kB
of spill. The batch size is therefore bounded by register pressure, not by the
algorithm, and the sweet spot has to be measured -- N = 4 already gets 2.1x of
the available 3.3x, so a small batch may be the whole win.

### Ranked, for four-centre work through d

1. **Batched multi-density Fock.** ~100 right-hand sides in CPHF, up to 3.3x,
   and nothing else on this list is worth as much. Bounded by registers, so
   measure N rather than assuming large is better.
2. **`k_scale` and `j_scale`.** 30 call sites in mqc's libcint backend, and
   CPHF passes `k_scale` explicitly. EASY here despite the folded digestion:
   they are per-call scalars, so they multiply the J and K terms inside the
   six atomic updates and no unfolding into separate J and K matrices is
   needed. `FGINT_FOCK6`'s 22% stays.
3. **Non-symmetric density.** Two call sites, both CPHF
   (`build_fock_direct_nosym`). The folded eight-fold digestion assumes
   `D_mn = D_nm`; a general density needs the permutations written out
   without that assumption. A separate kernel variant, not a flag.
4. **Multipole integrals through octupole**, for `multipole_matrices`. One-
   electron, and the separable machinery in `gen_1e.py` already covers the
   shape -- `r^n` is multiplicative, so it factorises exactly as overlap does,
   with the 1-D tables extended by the operator's own degree. Cheap to add,
   and out of scope while the focus is four-centre.

What terco already has and makeFP needs unchanged: the four-centre Fock at
LMAX=2 (d functions validated to 4.3e-13 against libfint and 1e-8 against
gpu4pyscf at 2773 functions), Schwarz and density screening, S, T and V.

What makeFP does NOT need: gradients. Which removes the whole derivative
generator from scope.
