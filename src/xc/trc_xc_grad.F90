!
! The exchange-correlation contribution to the nuclear gradient, LDA and
! GGA, restricted and unrestricted. Added into an existing `grad(3,natm)`,
! the same convention trc_1e_gradient and scf_grad_nuclear_repulsion use
! (trc_scf_grad.F90).
!
! THE FORMULA. For every atom A, the only thing that moves is the basis
! functions centred on it (the grid and its Becke weights are held fixed
! -- see GRID RESPONSE below), so differentiating
!
!     E_xc = sum_g w_g f(rho_g, sigma_g)
!
! with rho_g = sum_uv D_uv chi_u(g) chi_v(g) and sigma_g = |grad rho_g|^2
! gives, restricted,
!
!     dE_xc/dR_{A,k} = -2 sum_{u in A} sum_v D_uv sum_g w_g
!                      [ vrho_g dchi_u/dx_k(g) chi_v(g)
!                        + 2 vsigma_g ( chi_v(g) (grad_rho_g . H_u(g))_k
!                                       + dchi_u/dx_k(g) (grad_rho_g . grad_chi_v(g)) ) ]
!
! where H_u is chi_u's Hessian. LDA is the same formula with vsigma = 0 (a
! pure-LDA kernel composition already returns vsigma = 0 from
! xc_eval_point, so nothing here is skipped for correctness -- only for
! cost: LDA never collocates a Hessian). Unrestricted replaces the single
! vrho/vsigma pair with the five spin-resolved derivatives
! (va, vb, v_aa, v_ab, v_bb) and a density matrix per channel; see
! xc_grad_term for the channel-by-channel expansion. Checked,
! term for term, against metalquicha's independently-validated
! accumulate_channel/accumulate_gga_channel (mqc_czt_gradient.f90).
!
! GRID RESPONSE -- IMPLEMENTED, ON BY DEFAULT (`grid_response`, default
! .true.; .false. recovers the Pulay-only gradient above, e.g. to match
! pyscf's own `grid_response=False` default -- check_xc_grad_vs_pyscf.F90).
!
! The quadrature grid moves with the nuclei exactly as the basis does, and
! terco's own E_xc = sum_g w_g rho_g eps_g is a function of EVERY nuclear
! position through THREE channels, not one:
!
!   (1) the basis functions chi_u(r - R_u) centred on each atom -- THE
!       FORMULA above, the Pulay term;
!   (2) the grid point r_g = R_{owner(g)} + Delta_g, a FIXED offset from
!       whichever atom's atomic grid produced it, so it travels rigidly
!       with that atom;
!   (3) the Becke partition weight w_g = qw_g * P_g({R_A}), which depends
!       on every nuclear position through the confocal coordinates
!       mu_AB = (|r_g-R_A| - |r_g-R_B|)/R_AB (trc_dft_partition.F90) --
!       and, through r_g itself, on channel (2) as well.
!
! Channel (1) is unchanged by this round. Channels (2) and (3) are the
! "grid response", and they split into three pieces once the chain rule is
! carried through properly (not two, which is the trap: forgetting that
! w_g ALSO depends on the point's own position, through (3), double-counts
! nothing but drops a term -- see POINT C below):
!
!   A. POINT MOTION OF rho/eps. d(rho_g eps_g)/dR_B, at FIXED weight, for
!      g owned by B: the chain rule through r_g's dependence on R_B gives
!      + the SAME per-(point, function) force density the Pulay term
!      sums over points (xc_grad_term), but summed over EVERY local
!      function of the batch
!      (not just the ones on one atom) and landing on the POINT's owner
!      (not the function's), with the OPPOSITE sign:
!
!          dE_xc/dR_{B,k} += sum_{g: owner(g)=B} w_g [ vrho_g (drho_g/dx_k)
!                             + vsigma_g (dsigma_g/dx_k) ]
!
!      dsigma_g/dx_k = 2 (H_rho(g) . grad_rho_g)_k needs the density
!      HESSIAN at the point, which trc_xc_grad already collocates
!      (`hchi`) for the Pulay term's own vsigma piece. Its vsigma part
!      needs R D, R_u = grad rho . grad chi_u contracted with the density
!      the way X = chi D is chi's, which is a GEMM over the whole batch, so
!      the point-motion sum runs in its own kernel (xc_grad_pointmotion)
!      after it, one thread per point, keyed by `xg%atom` (the point's
!      owner), not `ao_atom`.
!
!   B. PARTITION WEIGHT, explicit nuclear dependence, AT FIXED r_g.
!      `trc_dft_partition.F90`'s `dft_partition_becke_derivatives` (already in
!      the tree, written for this and unused until this branch) returns
!      exactly dP_g/dR_A for every (atom, point) pair, carrying the VALUE
!      and the DERIVATIVE through the Becke product together rather than
!      as a logarithmic derivative (that routine's own header: a cutoff
!      factor legitimately underflows to zero far from the molecule, where
!      a log-derivative would divide by it). Multiplying by the quadrature
!      weight alone (`qw`, NOT `w` -- the partition's own factor must not
!      be double-counted) and by the point's energy density:
!
!          dE_xc/dR_{A,k} += sum_g qw_g (rho_g eps_g) dP_g/dR_{A,k}
!
!      for EVERY atom A, not only g's owner -- a distant atom's weight
!      still responds to a point it does not own.
!
!   C. PARTITION WEIGHT, motion of THE POINT ITSELF. P_g is a function of
!      inter-particle differences only (every `atom_r(i)` in
!      `partition_point` is |r_g - R_i|, and `shift`/`inv_distance` depend
!      only on the R_A themselves), so translating r_g AND every R_A by
!      the same vector leaves P_g unchanged:
!
!          grad_{r_g} P_g + sum_A dP_g/dR_A = 0
!          =>  grad_{r_g} P_g = - sum_A dP_g/dR_A
!
!      and since r_g only moves with its OWNER, this contributes
!
!          dE_xc/dR_{B,k} -= sum_{g: owner(g)=B} qw_g (rho_g eps_g)
!                             sum_A dP_g/dR_{A,k}
!
!      i.e. the negative of B's own sum, landed entirely on the point's
!      owner, for FREE out of the same dP/dR B already computes -- no new
!      derivative. trc_dft_partition's dft_partition_becke_grad_contract does
!      B and C together, on the device, one point per thread, with no
!      (3, natm, npts) array: it takes qw_g rho_g eps_g per point and adds
!      straight into a (3, natm) gradient. This is the piece a
!      two-term ("basis moves, grid moves") account of the gradient misses
!      silently: on water it is the difference between a ~1e-7-level
!      residual and a ~1 Hartree/Bohr one (metalquicha's own
!      mqc_czt_gradient.f90%xc_gradient derives and measures the identical
!      split; B and C below are a direct, GPU-first port of its algebra,
!      not an independent re-derivation).
!
! WHOSE BECKE PARTITION. Before porting B/C, checked against
! `trc_dft_partition.F90` itself, not assumed: terco's grid defaults to
! PARTITION_BECKE (three iterations of Becke's p(x)=x(3-x^2)/2) with
! ADJUST_TREUTLER (size adjustment on sqrt(Bragg radius), clamped to
! +-0.5) -- `trc_dft_grid.F90`'s own `used_scheme`/`used_adjust` defaults.
! This is EXACTLY pyscf's own default (`pyscf/dft/gen_grid.py`:
! `atomic_radii=radi.BRAGG_RADII`, `radii_adjust=radi.treutler_atomic_radii_adjust`,
! `becke_scheme=original_becke`, `radi_method=radi.treutler`,
! `prune=nwchem_prune`, `level=3`), term for term: pyscf's size-adjustment
! `a[i,j] = 0.25*(rad[j]/rad[i] - rad[i]/rad[j])` (`radi.py`'s
! `treutler_atomic_radii_adjust`) is terco's `shift(i,j)` in
! `size_adjustment` by direct substitution, and pyscf's cutoff
! `g + a[i,j]*(1-g**2)` is terco's `nu = mu + shift(i,j)*(1-mu*mu)`. That
! both codes land on IDENTICAL defaults is not a coincidence this header
! assumes -- check_xc_grad_vs_pyscf's `_gridresp_` cases below measure the
! consequence directly.
!
! VALIDATION, measured (test/check_xc_grad_fd.F90 and
! test/check_xc_grad_vs_pyscf.F90's `_gridresp_` cases; exact numbers in
! those tests' own output, not duplicated here to rot):
!
!   - FD, grid MOVING with the atoms (check_xc_grad_fd.F90: a central
!     difference of terco's own reconverged KS-DFT energy, grid rebuilt at
!     every displaced geometry). OFF (Pulay only) reproduces the OLD
!     numbers unchanged (water/SVWN ~8e-6 at level 3 / ~7e-7 at level 5,
!     water/PBE ~7e-6 / ~1e-6, NH2/PBE unrestricted ~4e-6 / ~9e-7) --
!     shrinking with the grid, the signature of a quadrature artefact. ON
!     collapses every one of those six cases to ~1-2e-7, THE SAME NUMBER
!     at level 3 AND level 5: the grid-LEVEL dependence is gone, which is
!     what "the physical gradient" means here -- what remains is the FD
!     stencil's own floor (H=1e-3 central difference against a
!     reconverged SCF energy), not a residual quadrature artefact.
!   - Translational invariance (sum over atoms): OFF is NOT required to
!     vanish and measured does not (~1e-5 to ~1e-7, shrinking with the
!     grid, the same artefact as above). ON is ~1e-13 to ~1e-14 on every
!     case -- six to seven orders of magnitude tighter, and flat across
!     grid levels.
!   - Against pyscf's OWN `grid_response=True` (`grad.RKS`/`grad.UKS`,
!     `get_vxc_full_response`/`grids_response_cc`), at terco's own
!     converged orbitals: pyscf's grid cannot be forced to terco's literal
!     point LIST (its `grids_response_cc` walks `grids.atom_grids`, built
!     only by `Grids.build()`, not by overriding `coords`/`weights` after
!     the fact -- test/pyscf_refs/generate.py's `gen_xc_grad` docstring),
!     so this builds a plain pyscf `Grids(mol)` at the SAME level instead
!     and lets it run ITS OWN defaults -- which the WHOSE BECKE PARTITION
!     paragraph above establishes are terco's defaults, term for term.
!     worst|terco - pyscf|, RKS water/SVWN and water/PBE: ~5-6e-13 --
!     essentially the two codes' shared double-precision floor. UKS
!     NH2/PBE (the cross-spin Q2/rmatx terms, and the only genuinely
!     open-shell case here): ~2.5e-11 -- still far below any grid or SCF
!     artefact, with more floating-point operations per point to
!     accumulate it over. Two independently written implementations of
!     the same physics agreeing this closely is a stronger check than
!     either FD comparison alone.
!
! Omitting the term (`grid_response=.false.`) stays available for a
! caller matching a reference that itself omits it (pyscf's own default),
! or isolating the Pulay term the way check_xc_grad_fixedgrid_fd.F90 does.
!
! BATCHING -- the energy path's (trc_xc.F90's xc_drive): tiles padded to
! (XC_NP_PAD, XC_NL_PAD) multiples, stored by padded shape so each run of
! equal shape is one strided-batched GEMM, every kernel one launch over a
! chunk of batches, and the workspace allocated once, at the largest
! chunk. The Pulay contraction needs only X = chi D and, for a GGA, R D
! (xc_grad_drive's own header): the (3 nloc x nloc) GEMMs the per-batch
! version formed, and their tiles, are gone. Measured, def2-SVP PBE on a
! V100: w20 1.5 s -> 0.18 s for everything but the grid build and the
! partition response; see the commit that made the change for the rest.
!
! CARTESIAN ONLY. XC is Cartesian throughout terco (the SCF folds to
! spherical); `b` and `dens` here are exactly what trc_scf_grad.F90 hands
! in, with no spherical folding of its own.
!
module trc_xc_grad
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_basis_types, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_eval_point, xc_eval_point_polar, &
                                 XC_FAMILY_LDA, XC_FAMILY_GGA, XC_MAX_KERNELS
   use trc_collocation, only: collocation_shell, collocation_shell_hessian, NCART_MAX, NHESS
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_xc, only: XC_NP_PAD, XC_NL_PAD, xc_sort_by_shape
   use trc_dft_grid, only: dft_grid_t, dft_grid_build, DEFAULT_GRID_LEVEL
   use trc_dft_partition, only: dft_partition_becke_grad_contract
   use trc_error, only: error_t
   use trc_step_table, only: step_table_t
   implicit none
   private

   public :: xc_grad_add, XC_GRAD_RHO_TOL, xc_grad_drive
   public :: xc_grad_options_t, xc_grad_options, xc_grad_point_io_t
   public :: XC_GRAD_NSTEP, xc_grad_step_names
   public :: ST_GRID, ST_POINTS, ST_PDERIV, ST_RESP
   integer, parameter :: XC_GRAD_NSTEP = 10
   !! Steps of the timing table `xc_grad_add` fills (see
   !! `xc_grad_step_names`); a caller's step_table_t is initialised with them.
   integer, parameter :: ST_GRID = 1, ST_SCREEN = 2, ST_ALLOC = 3, ST_COLLOC = 4, ST_GATHER = 5, &
                         ST_GEMM = 6, ST_POINTS = 7, ST_SCATTER = 8, ST_PDERIV = 9, ST_RESP = 10

   integer, parameter :: XC_GRAD_NPC = 12
   !! Per-point coefficients the force kernels read, at most (unrestricted):
   !! w, v_a, v_b, v_aa, v_ab, v_bb, grad rho_a, grad rho_b. Restricted uses 6.
#ifdef TRC_DEVICE
   integer, parameter :: XC_GRAD_LANES = 32
#else
   integer, parameter :: XC_GRAD_LANES = 1
#endif
   !! Threads that share one function's Pulay sum on the device, adjacent
   !! threads adjacent points (xc_grad_forces); one on the host, which
   !! walks a column in one iteration.
#ifdef TRC_DEVICE
   integer(int64), parameter :: XC_GRAD_TILE_CAP = 8388608_int64
#else
   integer(int64), parameter :: XC_GRAD_TILE_CAP = 262144_int64
#endif
   !! Doubles in one tile of a chunk (chi, say); a GGA chunk holds about
   !! thirteen such tiles, so 64 MB here is about 0.8 GB of device memory.
   !! On the host, 2 MB, so a chunk's tiles stay in cache from one kernel
   !! to the next: each kernel is a pass over every tile it reads, and at
   !! 64 MB every pass went to memory (w20/def2-SVP PBE, one core:
   !! collocate 3.9 -> 3.0 s, points 1.6 -> 1.2 s).

   real(dp), parameter :: XC_GRAD_RHO_TOL = 1.0e-12_dp
   !! Same role and value as trc_xc.F90's XC_RHO_TOL: a batch whose bound on
   !! rho is below this contributes nothing and is skipped before any
   !! device work for it.
   integer, parameter :: XC_GRAD_DEFAULT_MAX_PTS = 512
   !! Points per batch of a grid this module builds itself; trc_scf_context
   !! uses the same for the SCF's grid.

   type :: xc_grad_options_t
      !! The settings of one gradient term (xc_grad_add, and vv10_grad_add,
      !! which reads all but `rho_tol`). Build it with `xc_grad_options`,
      !! whose absent arguments leave these defaults.
      integer :: grid_level = DEFAULT_GRID_LEVEL
      !! Level of the DFT grid built here when `grid`/`xg` are not given.
      integer :: grid_max_pts = XC_GRAD_DEFAULT_MAX_PTS
      !! Points per batch of that grid.
      real(dp) :: rho_tol = XC_GRAD_RHO_TOL
      !! Batch screen on the density bound. VV10 always screens at its own
      !! VV10_RHO_TOL (trc_vv10_grad's header) and ignores this.
      logical :: grid_response = .true.
      !! See this module's header, GRID RESPONSE. OFF recovers the
      !! Pulay-only gradient (what check_xc_grad_fixedgrid_fd isolates, and
      !! what pyscf's own `grid_response=False` default compares against).
      type(step_table_t), pointer :: timing => null()
      !! Device-synced time per step, ADDED into; initialised by the caller
      !! with `xc_grad_step_names()`. Unassociated, a private table is kept
      !! and dropped.
      type(dft_grid_t), pointer :: grid => null()
      type(trc_xc_grid_t), pointer :: xg => null()
      !! The SCF's own grid, already batched and device-resident (`xg`).
      !! Used only when both are associated; see xc_grad_add's header.
   end type xc_grad_options_t

   type :: xc_grad_point_io_t
      !! What xc_grad_drive's point kernel exchanges with its caller, per
      !! point (per batch for `active`). Either `fxc` is associated, and the
      !! functional is evaluated and its rho*eps (no weight) written there
      !! for the partition response, left zero on a batch the screen skips;
      !! or `vrho` is, and a potential already in hand (VV10's, restricted
      !! GGA) is read from `grho`/`vrho`/`vsigma` on the device instead of
      !! evaluating the functional. All are the caller's arrays, mapped by
      !! the caller.
      real(dp), pointer, contiguous :: fxc(:) => null()
      !! (npts)
      integer, pointer, contiguous :: active(:) => null()
      !! (nbatch): the caller's own batch screen (1 = keep) in place of
      !! the density bound; optional with either form.
      real(dp), pointer, contiguous :: grho(:, :) => null()
      !! (3, npts)
      real(dp), pointer, contiguous :: vrho(:) => null()
      !! (npts)
      real(dp), pointer, contiguous :: vsigma(:) => null()
      !! (npts)
   end type xc_grad_point_io_t

   type :: xc_grad_work_t
      !! xc_grad_drive's workspace: the batch layout, the current chunk's
      !! storage order and offsets, and the tiles, allocated once at the
      !! largest chunk. Its arrays are mapped one by one, never the type
      !! (.claude/context/gpu-pitfalls.md), and reach the kernels as
      !! explicit-shape arguments.
      integer :: nspin = 1
      logical :: is_gga = .false.
      integer :: npc = 6      !! point coefficients per point (pc's columns)
      integer :: nhs = 1      !! Hessian components collocated: NHESS for a GGA
      integer(int64) :: ldt = 1, ldh = 1, ldr = 1, ldx = 1
      !! Leading dimensions of chi/gchi/xt, hchi, rt/rd, rdx.
      integer :: ldd = 1, ldp = 1, ldf = 1, ldq = 1
      !! Leading dimensions of dloc, pc, fp, qp/eb.
      integer :: nb = 0, p0 = 1, np = 0, npairs = 0, nwork = 0, nmu = 0
      !! The current chunk: stored batches, first point, points, (u,v)
      !! pairs, collocation items, functions of its active batches.
      integer(int64) :: nchi = 1
      !! The current chunk's tile length.
      integer, allocatable :: ao_atom(:)   !! (nao): the atom each function sits on
      integer, allocatable :: active(:), npp(:), nlp(:), c_off(:), pr_off(:)
      !! Per batch: kept, padded points and functions, tile and pair offsets.
      integer, allocatable :: ord(:), pb_off(:), wk_off(:), mw_off(:)
      !! Per stored batch of the current chunk: storage order, and the pair,
      !! collocation and force work offsets.
      real(dp), allocatable :: chi(:), gchi(:, :), hchi(:, :), xt(:, :), rt(:, :), rd(:, :), rdx(:, :)
      real(dp), allocatable :: dloc(:, :), pc(:, :), fp(:, :), qp(:, :), eb(:, :)
   end type xc_grad_work_t

contains

   !
   ! An xc_grad_options_t from the caller's own optional settings, passed
   ! straight through: an absent one keeps its default. `timing`, `grid`
   ! and `xg` are pointed at, so the caller's actuals need TARGET and must
   ! outlive the options. `grid` and `xg` are taken only together.
   !
   function xc_grad_options(grid_level, grid_max_pts, rho_tol, grid_response, timing, grid, xg) result(options)
      integer, intent(in), optional :: grid_level, grid_max_pts
      real(dp), intent(in), optional :: rho_tol
      logical, intent(in), optional :: grid_response
      type(step_table_t), intent(inout), optional, target :: timing
      type(dft_grid_t), intent(in), optional, target :: grid
      type(trc_xc_grid_t), intent(in), optional, target :: xg
      type(xc_grad_options_t) :: options

      if (present(grid_level)) options%grid_level = grid_level
      if (present(grid_max_pts)) options%grid_max_pts = grid_max_pts
      if (present(rho_tol)) options%rho_tol = rho_tol
      if (present(grid_response)) options%grid_response = grid_response
      if (present(timing)) options%timing => timing
      if (present(grid) .and. present(xg)) then
         options%grid => grid
         options%xg => xg
      end if
   end function xc_grad_options

   pure function xc_grad_step_names() result(names)
      !! The step names of `xc_grad_add`'s timing table, in order.
      character(len=16) :: names(XC_GRAD_NSTEP)
      names = [character(len=16) :: "grid", "screen", "alloc", "collocate", "gather", &
               "gemm", "points", "scatter", "partition deriv", "response"]
   end function xc_grad_step_names

   !
   ! The one entry point. Builds its own DFT quadrature at `b`'s geometry
   ! by default (trc_scf_run keeps neither its grid nor its functional
   ! alive past convergence, the same reason trc_scf_gradient rebuilds the
   ! pair list and the ERI engine) -- UNLESS the caller already has one at
   ! this exact geometry and settings and hands it in through `grid`/`xg`
   ! (trc_capi.F90's `trc_run_gradient`, reusing the SCF's own grid cache:
   ! see trc_scf_context.F90's `trc_grid_cache_t`). Present, both together,
   ! no rebuild and no upload happen at all, and neither is released here
   ! -- the caller owns them, exactly as it owns `dfjk`/`aux` in
   ! trc_scf_grad.F90's own `--df` branch. Absent (the default, and every
   ! existing caller), behaviour is unchanged: refuses what it cannot do,
   ! and ADDS the result into `grad`.
   !
   subroutine xc_grad_add(b, lp, func, nspin, density, grad, options)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_xc_functional_t), intent(in) :: func
      integer, intent(in) :: nspin
      real(dp), intent(in) :: density(b%nao, b%nao, nspin)
      !! The converged density, per spin channel. Present on both host and
      !! device: the GEMM pipeline below reads the device copy, the
      !! per-batch screen reads the host one (trc_scf_gradient's own
      !! convention -- see its header -- keeps both valid at this point).
      real(dp), intent(inout) :: grad(3, b%natm)
      !! ADDED into; device-resident for the whole call, exactly as
      !! trc_1e_gradient and scf_grad_nuclear_repulsion require of their
      !! own `grad` argument.
      type(xc_grad_options_t), intent(in), optional :: options
      !! Grid, screen, grid response, timing and a grid to reuse (REUSE in
      !! this subroutine's header); absent, every default.

      type(xc_grad_options_t) :: opts
      type(dft_grid_t), target :: grid_own
      type(trc_xc_grid_t), target :: xg_own
      type(dft_grid_t), pointer :: gp
      type(trc_xc_grid_t), pointer :: xgp
      logical :: reuse_grid
      type(step_table_t), target :: tbl_own
      type(step_table_t), pointer :: tbl
      type(xc_grad_point_io_t) :: pio
      type(timer_type) :: tc
      type(error_t) :: err
      integer, allocatable :: at_elem(:)
      real(dp), allocatable :: grad_extra(:, :), pcoef(:)
      real(dp), allocatable, target :: fxc(:)

      ! --- refuse what this module cannot do, with a reason -------------
      !
      ! Meta-GGA and double hybrids are not merely unimplemented here: on
      ! this branch `trc_xc_functional_t` has no meta-GGA family and no
      ! MP2-fraction field at all (see its own header), so neither can
      ! reach a converged SCF result in the first place -- there is no
      ! flag left to check for a double hybrid. The check below is for
      ! what CAN reach here that USED to be missing a kernel: VV10, whose
      ! nonlocal correlation gradient is now trc_vv10_grad%vv10_grad_add's
      ! job, called by trc_scf_gradient next to this routine rather than
      ! from inside it. `func%nlc_b`/`nlc_c` play NO role in this routine at
      ! all -- they gate the SEPARATE nonlocal double-sum, not the semilocal
      ! `kid`/`coef` kernel list this routine reads, so a VV10 functional
      ! (wB97X-V) reaches this point and gets exactly its correct semilocal
      ! GGA gradient term, same as any other GGA.
      !
      ! RANGE-SEPARATED HYBRIDS (func%omega/exx_lr /= 0, wb97x and
      ! friends) are NOT refused: `func%kid`/`func%coef` (read by
      ! xc_grad_drive below) are the functional's SEMILOCAL LDA/GGA
      ! kernels only -- exactly the same field an ordinary GGA hybrid
      ! (b3lyp) uses -- entirely separate from `exx`/`exx_lr`/`omega`,
      ! which only ever reach eri%gradient/eri%fock_resident's exact- and
      ! long-range-exchange passes (trc_scf_grad.F90, trc_scf_fock.F90).
      ! So what this routine computes for an RSH functional IS its
      ! correct semilocal (GGA) exchange-correlation gradient term; the
      ! erf(omega r)/r long-range exchange term is trc_scf_grad.F90's own
      ! job (its own header, RANGE-SEPARATED HYBRIDS), added next to this
      ! one rather than inside it.
      !
      if (func%family /= XC_FAMILY_LDA .and. func%family /= XC_FAMILY_GGA) then
         error stop "trc_xc_grad: only LDA and GGA exchange-correlation gradients "// &
            "are implemented (meta-GGA is not representable in trc_xc_functional_t "// &
            "on this branch, so this should be unreachable)"
      end if

      if (present(options)) opts = options
      if (associated(opts%timing)) then
         tbl => opts%timing
      else
         call tbl_own%init("xc grad", xc_grad_step_names())
         tbl => tbl_own
         opts%timing => tbl_own
      end if

      reuse_grid = associated(opts%grid) .and. associated(opts%xg)
      call tbl%step_on(lp, ST_GRID, tc)
      allocate (at_elem(b%natm), source=b%at_elem)
      if (reuse_grid) then
         ! The SCF's own grid, already batched and on the device -- no
         ! rebuild, no re-upload. See this subroutine's own header, REUSE.
         gp => opts%grid
         xgp => opts%xg
      else
         call dft_grid_build(b%at_r, at_elem, grid_own, err, level=opts%grid_level)
         if (err%has_error()) error stop "trc_xc_grad: "//trim(err%get_message())
         ! `atom`/`qw` are only ever read by the grid-response term below
         ! (`build` leaves both unallocated; `build_from_grid` fills
         ! them); built unconditionally anyway --
         ! they cost one more integer and one more real(dp) per point, next to
         ! the GEMMs a surviving batch pays for.
         call xg_own%build_from_grid(grid_own, b, max_pts=opts%grid_max_pts)
         call xg_own%to_device()
         gp => grid_own
         xgp => xg_own
      end if
      call tbl%step_off(lp, ST_GRID, tc)

      if (xgp%npts > 0) then
         ! The weight-response and point-motion-of-the-weight terms (this
         ! module's header, GRID RESPONSE, B and C) go into their own
         ! array, merged into the device-resident `grad` once at the end.
         ! They need every point's rho*eps (`fxc`, filled by the drive's
         ! point kernel and left zero where the screen skips a batch).
         allocate (grad_extra(3, b%natm), fxc(xgp%npts))
         grad_extra = 0.0_dp
         fxc = 0.0_dp
         !$acc enter data copyin(fxc)

         pio%fxc => fxc
         call xc_grad_drive(b, lp, func, density, xgp, grad, opts, pio)

         if (opts%grid_response) then
            ! The partition's own response, both pieces at once and with
            ! no (3, natm, npts) array: dft_partition_becke_grad_contract takes
            ! the per-point factor qw*rho*eps (the quadrature weight alone,
            ! not w: the partition's own factor is what is differentiated).
            call tbl%step_on(lp, ST_PDERIV, tc)
            allocate (pcoef(xgp%npts), source=xgp%qw)
            !$acc enter data copyin(pcoef)
            call xc_grad_weight_coef(xgp%npts, fxc, pcoef)
            call dft_partition_becke_grad_contract(xgp%r, b%at_r, at_elem, xgp%atom, pcoef, gp%partition, &
                                                   grad_extra, err)
            if (err%has_error()) error stop "trc_xc_grad: "//trim(err%get_message())
            !$acc exit data delete(pcoef)
            deallocate (pcoef)
            call tbl%step_off(lp, ST_PDERIV, tc)
            call tbl%step_on(lp, ST_RESP, tc)
            !$acc update self(grad)
            grad = grad + grad_extra
            !$acc update device(grad)
            call tbl%step_off(lp, ST_RESP, tc)
         end if
         !$acc exit data delete(fxc)
         deallocate (grad_extra, fxc)
      end if

      if (.not. reuse_grid) then
         call xg_own%release()
         call grid_own%destroy()
      end if
      deallocate (at_elem)
   end subroutine xc_grad_add

   !
   ! Every surviving batch, in chunks, the energy path's way (trc_xc.F90's
   ! xc_drive): each batch's tiles are padded to (XC_NP_PAD, XC_NL_PAD)
   ! multiples and stored in order of padded shape, so each run of equal
   ! shape is ONE strided-batched GEMM, and every kernel is one launch over
   ! the whole chunk. The workspace is sized once, for the largest chunk.
   !
   ! THE FORCE WITHOUT THE BIG GEMMS. The Pulay term contracts, for every
   ! local function mu and direction k,
   !
   !     F_k(mu) = sum_nu D(mu,nu) sum_g [P_k(g,mu) chi_nu(g) + Q_k(g,mu) R_nu(g)]
   !             = sum_g [P_k(g,mu) X(g,mu) + Q_k(g,mu) RD(g,mu)]
   !
   ! with X = chi D and RD = R D (D symmetric), so the (3 nloc x nloc)
   ! GEMMs pall^T chi and q^T R and their D-contraction collapse to a
   ! per-(point, function) product summed over the points. RD is the GEMM
   ! the point-motion term already needed; X is the one the density needs.
   ! P and Q are rebuilt from per-point coefficients (w, v_rho, v_sigma,
   ! grad rho) instead of being stored as tiles three columns per function.
   ! The point-motion term is the SAME per-(point, function) product summed
   ! over the functions instead, landing on the point's owner with the
   ! opposite sign -- this module's header, GRID RESPONSE, piece A.
   !
   subroutine xc_grad_drive(b, lp, func, density, xg, grad, options, points)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_xc_functional_t), intent(in) :: func
      !! LDA or GGA; its kernels are evaluated unless `points` brings a
      !! potential (then only its family is read).
      real(dp), intent(in), contiguous :: density(:, :, :)
      !! (b%nao, b%nao, nspin), on host and device.
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(inout) :: grad(3, b%natm)
      type(xc_grad_options_t), intent(in) :: options
      !! `rho_tol` and `grid_response` are read, and `timing`, which must be
      !! associated; the grid settings are the caller's business.
      type(xc_grad_point_io_t), intent(in) :: points

      type(xc_grad_work_t) :: work
      integer, allocatable :: ch_lo(:), ch_hi(:)
      integer(int64) :: need, nchi_max
      integer :: ib, npb, nloc, nch, ic, b0, b1, npairs, npairs_max, np_max, nb_max, nmu_c, nmu_max, nspin
      integer :: kid(XC_MAX_KERNELS)
      real(dp) :: coef(XC_MAX_KERNELS)
      logical :: is_gga, gresp
      type(step_table_t), pointer :: tbl
      type(timer_type) :: tc

      nspin = size(density, 3)
      is_gga = func%family == XC_FAMILY_GGA
      kid = func%kid
      coef = func%coef
      gresp = options%grid_response
      tbl => options%timing

      ! --- screen: which batches carry density at all -------------------
      call tbl%step_on(lp, ST_SCREEN, tc)
      allocate (work%active(xg%nbatch), work%npp(xg%nbatch), work%nlp(xg%nbatch))
      allocate (work%c_off(xg%nbatch + 1), work%pr_off(xg%nbatch + 1))
      do ib = 1, xg%nbatch
         npb = xg%b_off(ib + 1) - xg%b_off(ib)
         nloc = xg%b_aooff(ib + 1) - xg%b_aooff(ib)
         work%active(ib) = 0
         if (npb > 0 .and. nloc > 0) then
            if (associated(points%active)) then
               work%active(ib) = points%active(ib)
            else if (xc_grad_batch_active(ib, xg, b%nao, nspin, density, options%rho_tol)) then
               work%active(ib) = 1
            end if
         end if
         work%npp(ib) = ((npb + XC_NP_PAD - 1)/XC_NP_PAD)*XC_NP_PAD
         work%nlp(ib) = ((nloc*work%active(ib) + XC_NL_PAD - 1)/XC_NL_PAD)*XC_NL_PAD
      end do
      work%c_off = 0
      work%pr_off = 0

      ! The chunks, up front, so the workspace is allocated once at the
      ! largest one: as many batches as XC_GRAD_TILE_CAP doubles of one
      ! tile hold, and at least one.
      allocate (ch_lo(xg%nbatch), ch_hi(xg%nbatch))
      nch = 0
      nchi_max = 1
      npairs_max = 1
      np_max = 1
      nb_max = 1
      nmu_max = 1
      b0 = 1
      do while (b0 <= xg%nbatch)
         need = 0
         npairs = 0
         nmu_c = 0
         b1 = b0 - 1
         do ib = b0, xg%nbatch
            if (need + int(work%npp(ib), int64)*int(work%nlp(ib), int64) > XC_GRAD_TILE_CAP .and. ib > b0) exit
            need = need + int(work%npp(ib), int64)*int(work%nlp(ib), int64)
            npairs = npairs + work%nlp(ib)*work%nlp(ib)
            nmu_c = nmu_c + (xg%b_aooff(ib + 1) - xg%b_aooff(ib))*work%active(ib)
            b1 = ib
         end do
         nmu_max = max(nmu_max, nmu_c)
         nch = nch + 1
         ch_lo(nch) = b0
         ch_hi(nch) = b1
         nchi_max = max(nchi_max, need)
         npairs_max = max(npairs_max, npairs)
         np_max = max(np_max, xg%b_off(b1 + 1) - xg%b_off(b0))
         nb_max = max(nb_max, b1 - b0 + 1)
         b0 = b1 + 1
      end do
      call tbl%step_off(lp, ST_SCREEN, tc)

      ! --- the workspace, once ------------------------------------------
      call tbl%step_on(lp, ST_ALLOC, tc)
      work%nspin = nspin
      work%is_gga = is_gga
      call xc_grad_work_alloc(work, b, nchi_max, npairs_max, np_max, nb_max, nmu_max)
      call tbl%step_off(lp, ST_ALLOC, tc)

      do ic = 1, nch
         call xc_grad_work_layout(work, xg, ch_lo(ic), ch_hi(ic))
         if (work%nmu == 0) cycle

         ! Padding rows and columns of the GEMM inputs must be zero, not
         ! whatever the workspace last held: chi feeds X = chi D, R feeds
         ! RD = R D. D's own padding is written as zero by the gather.
         call tbl%step_on(lp, ST_COLLOC, tc)
         call xc_grad_zero(work%nchi, work%ldr, nspin, is_gga, work%chi, work%rt)
         call xc_grad_colloc_tiles(work%nwork, work%nb, work%ord, work%wk_off, b%maxnp, b%nshell, b%sh_l, &
                                   b%sh_np, b%sh_e, b%sh_c, b%sh_r, xg%npts, xg%r, xg%b_off, xg%b_shoff, &
                                   size(xg%b_sh), xg%b_sh, xg%b_shao, work%c_off, work%npp, xg%nbatch, work%ldt, &
                                   is_gga, work%ldh, work%nhs, work%chi, work%gchi, work%hchi)
         call tbl%step_off(lp, ST_COLLOC, tc)

         call tbl%step_on(lp, ST_GATHER, tc)
         call xc_grad_gather_tiles(work%npairs, work%nb, work%ord, work%pb_off, nspin, b%nao, density, &
                                   xg%b_aooff, size(xg%b_ao), xg%b_ao, work%nlp, xg%nbatch, work%ldd, work%dloc)
         call tbl%step_off(lp, ST_GATHER, tc)

         ! X = chi D, per spin: one strided-batched GEMM per run of shape.
         call tbl%step_on(lp, ST_GEMM, tc)
         call xc_grad_gemm_x(lp, work)
         call tbl%step_off(lp, ST_GEMM, tc)

         call tbl%step_on(lp, ST_POINTS, tc)
         if (associated(points%vrho)) then
            call xc_grad_points_ext(work%np, work%p0, xg%npts, xg%w, xg%batch_of, xg%b_off, xg%b_aooff, &
                                    work%c_off, work%npp, work%active, xg%nbatch, work%ldt, work%ldr, work%ldp, &
                                    work%npc, work%gchi, points%grho, points%vrho, points%vsigma, work%rt, work%pc)
         else
            call xc_grad_points(work%np, work%p0, nspin, is_gga, xg%npts, xg%w, xg%batch_of, xg%b_off, &
                                xg%b_aooff, work%c_off, work%npp, work%active, xg%nbatch, func%nk, kid, coef, &
                                work%ldt, work%ldr, work%ldp, work%npc, work%chi, work%gchi, work%xt, work%rt, &
                                work%pc, points%fxc)
         end if
         call tbl%step_off(lp, ST_POINTS, tc)

         ! RD = R D and, unrestricted, the cross-spin R_b D_a / R_a D_b.
         if (is_gga) then
            call tbl%step_on(lp, ST_GEMM, tc)
            call xc_grad_gemm_rd(lp, work)
            call tbl%step_off(lp, ST_GEMM, tc)
         end if

         ! The forces: Pulay onto each function's atom, point motion onto
         ! each point's owner.
         call tbl%step_on(lp, ST_SCATTER, tc)
         call xc_grad_forces(work, b, xg, gresp, grad)
         call tbl%step_off(lp, ST_SCATTER, tc)
      end do

      call tbl%step_on(lp, ST_ALLOC, tc)
      call xc_grad_work_release(work)
      deallocate (ch_lo, ch_hi)
      call tbl%step_off(lp, ST_ALLOC, tc)
   end subroutine xc_grad_drive

   !
   ! The workspace at the largest chunk (xc_grad_drive's chunk plan), and
   ! the AO -> atom map, built from the shell -> atom map `basis_build`
   ! already produced (sh_at), the same map the two-electron and
   ! one-electron gradients use. `work%nspin`/`is_gga` are set and
   ! `work%active`/`npp`/`nlp`/`c_off`/`pr_off` filled; everything goes to
   ! the device here, array by array.
   !
   subroutine xc_grad_work_alloc(work, b, nchi_max, npairs_max, np_max, nb_max, nmu_max)
      type(xc_grad_work_t), intent(inout) :: work
      type(trc_basis_t), intent(in) :: b
      integer(int64), intent(in) :: nchi_max
      integer, intent(in) :: npairs_max, np_max
      integer, intent(in) :: nb_max   !! stored batches in one chunk, at most
      integer, intent(in) :: nmu_max  !! local functions of active batches in one chunk, at most
      integer :: ish, i, nspin
      logical :: is_gga

      allocate (work%ao_atom(b%nao))
      do ish = 1, b%nshell
         do i = 0, (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2 - 1
            work%ao_atom(b%sh_ao(ish) + i) = b%sh_at(ish)
         end do
      end do

      nspin = work%nspin
      is_gga = work%is_gga
      work%npc = merge(6, XC_GRAD_NPC, nspin == 1)
      work%nhs = merge(NHESS, 1, is_gga)
      allocate (work%chi(nchi_max), work%gchi(nchi_max, 3), work%hchi(merge(nchi_max, 1_int64, is_gga), work%nhs))
      allocate (work%xt(nchi_max, nspin))
      if (is_gga) then
         allocate (work%rt(nchi_max, nspin), work%rd(nchi_max, nspin))
         allocate (work%rdx(merge(nchi_max, 1_int64, nspin == 2), nspin))
      else
         allocate (work%rt(1, nspin), work%rd(1, nspin), work%rdx(1, nspin))
      end if
      allocate (work%dloc(npairs_max, nspin), work%pc(np_max, work%npc))
      ! The force partials (xc_grad_forces): a row per (function, lane),
      ! and, on the host only, two rows per point.
      work%ldf = nmu_max*XC_GRAD_LANES
      work%ldq = merge(np_max, 1, XC_GRAD_LANES == 1)
      allocate (work%fp(work%ldf, 3), work%qp(work%ldq, 3), work%eb(work%ldq, 3))
      allocate (work%ord(nb_max), work%pb_off(nb_max + 1), work%wk_off(nb_max + 1), work%mw_off(nb_max + 1))
      work%ldt = nchi_max
      work%ldh = size(work%hchi, 1, int64)
      work%ldr = size(work%rt, 1, int64)
      work%ldx = size(work%rdx, 1, int64)
      work%ldd = npairs_max
      work%ldp = np_max
      !$acc enter data copyin(work%ao_atom, work%active, work%npp, work%nlp, work%c_off, work%pr_off) &
      !$acc      create(work%chi, work%gchi, work%hchi, work%xt, work%rt, work%rd, work%rdx, work%dloc, work%pc) &
      !$acc      create(work%fp, work%qp, work%eb, work%ord, work%pb_off, work%wk_off, work%mw_off)
   end subroutine xc_grad_work_alloc

   subroutine xc_grad_work_release(work)
      type(xc_grad_work_t), intent(inout) :: work
      !$acc exit data delete(work%ao_atom, work%active, work%npp, work%nlp, work%c_off, work%pr_off) &
      !$acc      delete(work%chi, work%gchi, work%hchi, work%xt, work%rt, work%rd, work%rdx, work%dloc, work%pc) &
      !$acc      delete(work%fp, work%qp, work%eb, work%ord, work%pb_off, work%wk_off, work%mw_off)
      deallocate (work%ao_atom, work%active, work%npp, work%nlp, work%c_off, work%pr_off)
      deallocate (work%chi, work%gchi, work%hchi, work%xt, work%rt, work%rd, work%rdx, work%dloc, work%pc)
      deallocate (work%fp, work%qp, work%eb, work%ord, work%pb_off, work%wk_off, work%mw_off)
   end subroutine xc_grad_work_release

   !
   ! Chunk b0..b1: its storage order (by padded shape) and the offsets in
   ! that order -- tiles (point fastest), (u,v) pairs, collocation work
   ! (shell, point), and force work (function) -- on host and device.
   !
   subroutine xc_grad_work_layout(work, xg, b0, b1)
      type(xc_grad_work_t), intent(inout) :: work
      type(trc_xc_grid_t), intent(in) :: xg
      integer, intent(in) :: b0, b1
      integer(int64) :: need
      integer :: k, ib, nb

      nb = b1 - b0 + 1
      work%nb = nb
      work%p0 = xg%b_off(b0)
      work%np = xg%b_off(b1 + 1) - work%p0
      do k = 1, nb
         work%ord(k) = b0 + k - 1
      end do
      call xc_sort_by_shape(nb, work%npp, work%nlp, work%ord)
      need = 0
      work%npairs = 0
      work%nwork = 0
      work%nmu = 0
      work%pb_off(1) = 0
      work%wk_off(1) = 0
      work%mw_off(1) = 0
      do k = 1, nb
         ib = work%ord(k)
         work%c_off(ib) = int(need)
         work%pr_off(ib) = work%npairs
         need = need + int(work%npp(ib), int64)*int(work%nlp(ib), int64)
         work%npairs = work%npairs + work%nlp(ib)*work%nlp(ib)
         work%pb_off(k + 1) = work%npairs
         work%nwork = work%nwork + (xg%b_shoff(ib + 1) - xg%b_shoff(ib))*(xg%b_off(ib + 1) - xg%b_off(ib)) &
                      *work%active(ib)
         work%wk_off(k + 1) = work%nwork
         work%nmu = work%nmu + (xg%b_aooff(ib + 1) - xg%b_aooff(ib))*work%active(ib)
         work%mw_off(k + 1) = work%nmu
      end do
      work%nchi = max(need, 1_int64)
      !$acc update device(work%c_off, work%pr_off, work%ord(1:nb), work%pb_off(1:nb+1), work%wk_off(1:nb+1)) &
      !$acc        device(work%mw_off(1:nb+1))
   end subroutine xc_grad_work_layout

   !
   ! X = chi D, per spin, for the current chunk: one strided-batched GEMM
   ! per run of equal padded shape.
   !
   subroutine xc_grad_gemm_x(lp, work)
      type(trc_linalg_t), intent(inout) :: lp
      type(xc_grad_work_t), intent(inout) :: work
      type(trc_gemm_batch_t) :: gemm_shape
      integer :: k0, k1, ib, cnt, mp, ml, s

      k0 = 1
      do while (k0 <= work%nb)
         call shape_run(k0, work%nb, work%ord, work%npp, work%nlp, k1, mp, ml)
         ib = work%ord(k0)
         cnt = k1 - k0 + 1
         if (ml > 0) then
            gemm_shape = trc_gemm_batch(mp, ml, ml, cnt)
            do s = 1, work%nspin
               call lp%gemm_batched(gemm_shape, work%chi(work%c_off(ib) + 1), work%dloc(work%pr_off(ib) + 1, s), &
                                    work%xt(work%c_off(ib) + 1, s))
            end do
         end if
         k0 = k1 + 1
      end do
   end subroutine xc_grad_gemm_x

   !
   ! RD = R D and, unrestricted, the cross-spin R_b D_a / R_a D_b, for the
   ! current chunk (a GGA only).
   !
   subroutine xc_grad_gemm_rd(lp, work)
      type(trc_linalg_t), intent(inout) :: lp
      type(xc_grad_work_t), intent(inout) :: work
      type(trc_gemm_batch_t) :: gemm_shape
      integer :: k0, k1, ib, cnt, mp, ml, s, co, po

      k0 = 1
      do while (k0 <= work%nb)
         call shape_run(k0, work%nb, work%ord, work%npp, work%nlp, k1, mp, ml)
         ib = work%ord(k0)
         cnt = k1 - k0 + 1
         if (ml > 0) then
            gemm_shape = trc_gemm_batch(mp, ml, ml, cnt)
            co = work%c_off(ib) + 1
            po = work%pr_off(ib) + 1
            do s = 1, work%nspin
               call lp%gemm_batched(gemm_shape, work%rt(co, s), work%dloc(po, s), work%rd(co, s))
            end do
            if (work%nspin == 2) then
               call lp%gemm_batched(gemm_shape, work%rt(co, 2), work%dloc(po, 1), work%rdx(co, 1))
               call lp%gemm_batched(gemm_shape, work%rt(co, 1), work%dloc(po, 2), work%rdx(co, 2))
            end if
         end if
         k0 = k1 + 1
      end do
   end subroutine xc_grad_gemm_rd

   !
   ! Host-side density screen, bounding rho on a batch through the same
   ! per-function maxima (`xg%b_amax`) trc_xc.F90's xc_screen_body uses --
   ! plain Fortran, not a device kernel, because `dens`'s host copy already
   ! holds the right numbers (trc_scf_gradient's `!$acc update self(dens)`,
   ! before this is ever called) and nloc is small enough that this costs
   ! nothing next to the GEMMs a surviving batch pays for.
   !
   logical function xc_grad_batch_active(ib, xg, nao, nspin, density, rtol) result(active)
      integer, intent(in) :: ib, nao, nspin
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: density(nao, nao, nspin)
      real(dp), intent(in) :: rtol
      real(dp) :: bound, su
      integer :: s, u, v, au, av
      bound = 0.0_dp
      do s = 1, nspin
         do u = xg%b_aooff(ib), xg%b_aooff(ib + 1) - 1
            au = xg%b_ao(u)
            su = 0.0_dp
            do v = xg%b_aooff(ib), xg%b_aooff(ib + 1) - 1
               av = xg%b_ao(v)
               su = su + abs(density(av, au, s))*xg%b_amax(v)
            end do
            bound = bound + su*xg%b_amax(u)
         end do
      end do
      active = bound > rtol
   end function xc_grad_batch_active

   ! The run of stored batches k0..k1 that share batch ord(k0)'s padded
   ! shape (mp, ml).
   subroutine shape_run(k0, nb, ord, npp, nlp, k1, mp, ml)
      integer, intent(in) :: k0, nb
      integer, intent(in) :: ord(:), npp(:), nlp(:)
      integer, intent(out) :: k1, mp, ml
      mp = npp(ord(k0))
      ml = nlp(ord(k0))
      k1 = k0
      do while (k1 < nb)
         if (npp(ord(k1 + 1)) /= mp .or. nlp(ord(k1 + 1)) /= ml) exit
         k1 = k1 + 1
      end do
   end subroutine shape_run

   !
   ! coef(g) = qw(g)*fxc(g): `coef` arrives holding qw, on the device.
   !
   subroutine xc_grad_weight_coef(npts, fxc, coef)
      integer, intent(in) :: npts
      real(dp), intent(in) :: fxc(npts)
      real(dp), intent(inout) :: coef(npts)
      integer :: g
      do concurrent(g=1:npts)
         coef(g) = coef(g)*fxc(g)
      end do
   end subroutine xc_grad_weight_coef

   ! Zero the tiles whose padding a GEMM reads: chi, and R for a GGA.
   subroutine xc_grad_zero(nchi, ldr, nspin, is_gga, chi, rt)
      integer(int64), intent(in) :: nchi, ldr
      integer, intent(in) :: nspin
      logical, intent(in) :: is_gga
      real(dp), intent(inout) :: chi(nchi), rt(ldr, nspin)
      integer(int64) :: i
      integer :: s
      do concurrent(i=1:nchi)
         chi(i) = 0.0_dp
      end do
      if (is_gga) then
         do s = 1, nspin
            do concurrent(i=1:nchi)
               rt(i, s) = 0.0_dp
            end do
         end do
      end if
   end subroutine xc_grad_zero

   !
   ! Which stored batch (1..nb, storage order) owns work item t: the last k
   ! with off(k) < t, and t's 0-based position inside it. trc_xc.F90's
   ! pair_batch, here because a `do concurrent` body calls only this
   ! module's own routines.
   !
   pure subroutine grad_item_batch(t, nb, off, k, tt)
      !$acc routine seq
      integer, intent(in) :: t, nb
      integer, intent(in) :: off(nb + 1)
      integer, intent(out) :: k, tt
      integer :: lo, hi, mid
      lo = 1
      hi = nb
      do while (lo < hi)
         mid = (lo + hi + 1)/2
         if (off(mid) < t) then
            lo = mid
         else
            hi = mid - 1
         end if
      end do
      k = lo
      tt = t - off(k) - 1
   end subroutine grad_item_batch

   !
   ! Collocation into the chunk's padded tiles: one thread per (shell,
   ! point) of every stored batch, point fastest, as trc_xc.F90's
   ! xc_collocate, with the Hessian for a GGA.
   !
   subroutine xc_grad_colloc_tiles(nwork, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                   npts, r, b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, ldt, &
                                   is_gga, ldh, nhs, chi, gchi, hchi)
      integer, intent(in) :: nwork, nb, maxnp, nshell, npts, nbsh, nbatch, nhs
      integer, intent(in) :: ord(nb), wk_off(nb + 1)
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: c_off(nbatch + 1), npp(nbatch)
      integer(int64), intent(in) :: ldt, ldh
      logical, intent(in) :: is_gga
      real(dp), intent(inout) :: chi(ldt), gchi(ldt, 3), hchi(ldh, nhs)
      integer :: t

      if (is_gga) then
         do concurrent(t=1:nwork)
            call xc_grad_colloc_tile2_body(t, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                           npts, r, b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, &
                                           ldt, chi, gchi, hchi)
         end do
      else
         do concurrent(t=1:nwork)
            call xc_grad_colloc_tile1_body(t, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                           npts, r, b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, &
                                           ldt, chi, gchi)
         end do
      end if
   end subroutine xc_grad_colloc_tiles

   pure subroutine xc_grad_colloc_tile1_body(t, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                             npts, r, b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, &
                                             ldt, chi, gchi)
      !$acc routine seq
      integer, intent(in) :: t, nb, maxnp, nshell, npts, nbsh, nbatch
      integer, intent(in) :: ord(nb), wk_off(nb + 1)
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: c_off(nbatch + 1), npp(nbatch)
      integer(int64), intent(in) :: ldt
      real(dp), intent(inout) :: chi(ldt), gchi(ldt, 3)
      integer :: k, ib, tt, npb, ks, gl, gg, ish, l, nc, m, ld, off, iao, i
      real(dp) :: d(3), cs(NCART_MAX), gs(3, NCART_MAX)

      call grad_item_batch(t, nb, wk_off, k, tt)
      ib = ord(k)
      npb = b_off(ib + 1) - b_off(ib)
      ks = tt/npb
      gl = tt - ks*npb
      gg = b_off(ib) + gl
      ish = b_sh(b_shoff(ib) + ks)
      iao = b_shao(b_shoff(ib) + ks)
      l = sh_l(ish)
      d(1) = r(1, gg) - sh_r(1, ish)
      d(2) = r(2, gg) - sh_r(2, ish)
      d(3) = r(3, gg) - sh_r(3, ish)
      call collocation_shell(l, sh_np(ish), sh_e(1:sh_np(ish), ish), sh_c(1:sh_np(ish), ish), d, cs, gs)
      nc = (l + 1)*(l + 2)/2
      ld = npp(ib)
      off = c_off(ib) + gl + 1
      do m = 1, nc
         i = off + (iao + m - 1)*ld
         chi(i) = cs(m)
         gchi(i, 1) = gs(1, m)
         gchi(i, 2) = gs(2, m)
         gchi(i, 3) = gs(3, m)
      end do
   end subroutine xc_grad_colloc_tile1_body

   pure subroutine xc_grad_colloc_tile2_body(t, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                             npts, r, b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, &
                                             ldt, chi, gchi, hchi)
      !$acc routine seq
      integer, intent(in) :: t, nb, maxnp, nshell, npts, nbsh, nbatch
      integer, intent(in) :: ord(nb), wk_off(nb + 1)
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: c_off(nbatch + 1), npp(nbatch)
      integer(int64), intent(in) :: ldt
      real(dp), intent(inout) :: chi(ldt), gchi(ldt, 3), hchi(ldt, NHESS)
      integer :: k, ib, tt, npb, ks, gl, gg, ish, l, nc, m, ld, off, iao, i
      real(dp) :: d(3), cs(NCART_MAX), gs(3, NCART_MAX), hs(NHESS, NCART_MAX)

      call grad_item_batch(t, nb, wk_off, k, tt)
      ib = ord(k)
      npb = b_off(ib + 1) - b_off(ib)
      ks = tt/npb
      gl = tt - ks*npb
      gg = b_off(ib) + gl
      ish = b_sh(b_shoff(ib) + ks)
      iao = b_shao(b_shoff(ib) + ks)
      l = sh_l(ish)
      d(1) = r(1, gg) - sh_r(1, ish)
      d(2) = r(2, gg) - sh_r(2, ish)
      d(3) = r(3, gg) - sh_r(3, ish)
      call collocation_shell_hessian(l, sh_np(ish), sh_e(1:sh_np(ish), ish), sh_c(1:sh_np(ish), ish), d, cs, gs, hs)
      nc = (l + 1)*(l + 2)/2
      ld = npp(ib)
      off = c_off(ib) + gl + 1
      do m = 1, nc
         i = off + (iao + m - 1)*ld
         chi(i) = cs(m)
         gchi(i, 1) = gs(1, m)
         gchi(i, 2) = gs(2, m)
         gchi(i, 3) = gs(3, m)
         hchi(i, 1) = hs(1, m)
         hchi(i, 2) = hs(2, m)
         hchi(i, 3) = hs(3, m)
         hchi(i, 4) = hs(4, m)
         hchi(i, 5) = hs(5, m)
         hchi(i, 6) = hs(6, m)
      end do
   end subroutine xc_grad_colloc_tile2_body

   !
   ! The density slice of every stored batch, one thread per (b, u, v),
   ! padding written as zero -- trc_xc.F90's xc_gather, with the leading
   ! dimension of the (once-allocated) workspace passed in.
   !
   subroutine xc_grad_gather_tiles(npairs, nb, ord, pb_off, nspin, nao, dmat, b_aooff, nbao, b_ao, nlp, &
                                   nbatch, ldd, dloc)
      integer, intent(in) :: npairs, nb, nspin, nao, nbao, nbatch, ldd
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(in) :: dmat(nao, nao, nspin)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(inout) :: dloc(ldd, nspin)
#ifdef TRC_DEVICE
      integer :: t
      do concurrent(t=1:npairs)
         call xc_grad_gather_tile_body(t, nb, ord, pb_off, nspin, nao, dmat, b_aooff, nbao, b_ao, nlp, nbatch, &
                                       ldd, dloc)
      end do
#else
      ! On the host, a batch per iteration, its tile column by column: the
      ! flat form's search for each element's batch, and its mod and
      ! divide, cost several times the copy itself.
      integer :: k
      do concurrent(k=1:nb)
         call xc_grad_gather_batch(k, nb, ord, pb_off, nspin, nao, dmat, b_aooff, nbao, b_ao, nlp, nbatch, ldd, &
                                   dloc)
      end do
#endif
   end subroutine xc_grad_gather_tiles

   ! Stored batch k's density tiles, as xc_grad_gather_tile_body writes
   ! them element by element.
   pure subroutine xc_grad_gather_batch(k, nb, ord, pb_off, nspin, nao, dmat, b_aooff, nbao, b_ao, nlp, nbatch, &
                                        ldd, dloc)
      !$acc routine seq
      integer, intent(in) :: k, nb, nspin, nao, nbao, nbatch, ldd
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(in) :: dmat(nao, nao, nspin)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(inout) :: dloc(ldd, nspin)
      integer :: ib, nloc, ld, a0, u, v, s, t

      ib = ord(k)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = nlp(ib)
      a0 = b_aooff(ib) - 1
      do s = 1, nspin
         do v = 1, ld
            t = pb_off(k) + (v - 1)*ld
            if (v > nloc) then
               do u = 1, ld
                  dloc(t + u, s) = 0.0_dp
               end do
            else
               do u = 1, nloc
                  dloc(t + u, s) = dmat(b_ao(a0 + u), b_ao(a0 + v), s)
               end do
               do u = nloc + 1, ld
                  dloc(t + u, s) = 0.0_dp
               end do
            end if
         end do
      end do
   end subroutine xc_grad_gather_batch

   pure subroutine xc_grad_gather_tile_body(t, nb, ord, pb_off, nspin, nao, dmat, b_aooff, nbao, b_ao, nlp, &
                                            nbatch, ldd, dloc)
      !$acc routine seq
      integer, intent(in) :: t, nb, nspin, nao, nbao, nbatch, ldd
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(in) :: dmat(nao, nao, nspin)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(inout) :: dloc(ldd, nspin)
      integer :: k, ib, tt, nloc, ld, u, v, s

      call grad_item_batch(t, nb, pb_off, k, tt)
      ib = ord(k)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = nlp(ib)
      u = mod(tt, ld) + 1
      v = tt/ld + 1
      if (u > nloc .or. v > nloc) then
         do s = 1, nspin
            dloc(t, s) = 0.0_dp
         end do
      else
         do s = 1, nspin
            dloc(t, s) = dmat(b_ao(b_aooff(ib) + u - 1), b_ao(b_aooff(ib) + v - 1), s)
         end do
      end if
   end subroutine xc_grad_gather_tile_body

   !
   ! One thread per point of the chunk: rho and grad rho from X = chi D,
   ! the functional, the point's coefficients for the force kernels (`pc`:
   ! w, the potential's derivatives, grad rho per spin), its energy density
   ! for the partition response, and, for a GGA, its row of
   ! R_u = grad rho . grad chi_u, the tile R D is formed from.
   !
   subroutine xc_grad_points(np, p0, nspin, is_gga, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, &
                             nbatch, nk, kid, coef, ldt, ldr, ldp, npc, chi, gchi, xt, rt, pc, fxc)
      integer, intent(in) :: np, p0, nspin, npts, nbatch, nk, ldp, npc
      logical, intent(in) :: is_gga
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), active(nbatch), kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(int64), intent(in) :: ldt, ldr
      real(dp), intent(in) :: chi(ldt), gchi(ldt, 3), xt(ldt, nspin)
      real(dp), intent(inout) :: rt(ldr, nspin), pc(ldp, npc), fxc(npts)
      integer :: g
      do concurrent(g=1:np)
         call xc_grad_point_body(g, p0, nspin, is_gga, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, &
                                 nbatch, nk, kid, coef, ldt, ldr, ldp, npc, chi, gchi, xt, rt, pc, fxc)
      end do
   end subroutine xc_grad_points

   pure subroutine xc_grad_point_body(g, p0, nspin, is_gga, npts, w, batch_of, b_off, b_aooff, c_off, npp, &
                                      active, nbatch, nk, kid, coef, ldt, ldr, ldp, npc, chi, gchi, xt, rt, &
                                      pc, fxc)
      !$acc routine seq
      integer, intent(in) :: g, p0, nspin, npts, nbatch, nk, ldp, npc
      logical, intent(in) :: is_gga
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), active(nbatch), kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(int64), intent(in) :: ldt, ldr
      real(dp), intent(in) :: chi(ldt), gchi(ldt, 3), xt(ldt, nspin)
      real(dp), intent(inout) :: rt(ldr, nspin), pc(ldp, npc), fxc(npts)
      integer :: gg, ib, nloc, ld, off, s, u, i
      real(dp) :: rho(2), grho(3, 2)
      real(dp) :: sigma, eps, vrho, vsigma
      real(dp) :: va, vb, s_aa, s_ab, s_bb, v_aa, v_ab, v_bb, x

      gg = p0 + g - 1
      ib = batch_of(gg)
      if (active(ib) == 0) return
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = npp(ib)
      off = c_off(ib) + (gg - b_off(ib)) + 1
      rho = 0.0_dp
      grho = 0.0_dp
      do s = 1, nspin
         do u = 1, nloc
            i = off + (u - 1)*ld
            x = xt(i, s)
            rho(s) = rho(s) + x*chi(i)
            grho(1, s) = grho(1, s) + 2.0_dp*x*gchi(i, 1)
            grho(2, s) = grho(2, s) + 2.0_dp*x*gchi(i, 2)
            grho(3, s) = grho(3, s) + 2.0_dp*x*gchi(i, 3)
         end do
      end do
      pc(g, 1) = w(gg)
      if (nspin == 1) then
         sigma = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         call xc_eval_point(nk, kid, coef, rho(1), sigma, eps, vrho, vsigma)
         pc(g, 2) = vrho
         pc(g, 3) = vsigma
         pc(g, 4) = grho(1, 1)
         pc(g, 5) = grho(2, 1)
         pc(g, 6) = grho(3, 1)
         fxc(gg) = rho(1)*eps
      else
         s_aa = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         s_bb = grho(1, 2)**2 + grho(2, 2)**2 + grho(3, 2)**2
         s_ab = grho(1, 1)*grho(1, 2) + grho(2, 1)*grho(2, 2) + grho(3, 1)*grho(3, 2)
         call xc_eval_point_polar(nk, kid, coef, rho(1), rho(2), s_aa, s_ab, s_bb, &
                                  eps, va, vb, v_aa, v_ab, v_bb)
         pc(g, 2) = va
         pc(g, 3) = vb
         pc(g, 4) = v_aa
         pc(g, 5) = v_ab
         pc(g, 6) = v_bb
         pc(g, 7) = grho(1, 1)
         pc(g, 8) = grho(2, 1)
         pc(g, 9) = grho(3, 1)
         pc(g, 10) = grho(1, 2)
         pc(g, 11) = grho(2, 2)
         pc(g, 12) = grho(3, 2)
         fxc(gg) = (rho(1) + rho(2))*eps
      end if
      if (is_gga) then
         do s = 1, nspin
            do u = 1, nloc
               i = off + (u - 1)*ld
               rt(i, s) = grho(1, s)*gchi(i, 1) + grho(2, s)*gchi(i, 2) + grho(3, s)*gchi(i, 3)
            end do
         end do
      end if
   end subroutine xc_grad_point_body

   !
   ! As xc_grad_points for a potential already evaluated on the grid
   ! (restricted GGA): the point's coefficients copied from it, and its
   ! row of R_u = grad rho . grad chi_u.
   !
   subroutine xc_grad_points_ext(np, p0, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, nbatch, &
                                 ldt, ldr, ldp, npc, gchi, grho, vrho, vsigma, rt, pc)
      integer, intent(in) :: np, p0, npts, nbatch, ldp, npc
      real(dp), intent(in) :: w(npts), grho(3, npts), vrho(npts), vsigma(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), active(nbatch)
      integer(int64), intent(in) :: ldt, ldr
      real(dp), intent(in) :: gchi(ldt, 3)
      real(dp), intent(inout) :: rt(ldr, 1), pc(ldp, npc)
      integer :: g
      do concurrent(g=1:np)
         call xc_grad_point_ext_body(g, p0, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, nbatch, &
                                     ldt, ldr, ldp, npc, gchi, grho, vrho, vsigma, rt, pc)
      end do
   end subroutine xc_grad_points_ext

   pure subroutine xc_grad_point_ext_body(g, p0, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, nbatch, &
                                          ldt, ldr, ldp, npc, gchi, grho, vrho, vsigma, rt, pc)
      !$acc routine seq
      integer, intent(in) :: g, p0, npts, nbatch, ldp, npc
      real(dp), intent(in) :: w(npts), grho(3, npts), vrho(npts), vsigma(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), active(nbatch)
      integer(int64), intent(in) :: ldt, ldr
      real(dp), intent(in) :: gchi(ldt, 3)
      real(dp), intent(inout) :: rt(ldr, 1), pc(ldp, npc)
      integer :: gg, ib, nloc, ld, off, u, i
      real(dp) :: grx, gry, grz

      gg = p0 + g - 1
      ib = batch_of(gg)
      if (active(ib) == 0) return
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = npp(ib)
      off = c_off(ib) + (gg - b_off(ib)) + 1
      grx = grho(1, gg)
      gry = grho(2, gg)
      grz = grho(3, gg)
      pc(g, 1) = w(gg)
      pc(g, 2) = vrho(gg)
      pc(g, 3) = vsigma(gg)
      pc(g, 4) = grx
      pc(g, 5) = gry
      pc(g, 6) = grz
      do u = 1, nloc
         i = off + (u - 1)*ld
         rt(i, 1) = grx*gchi(i, 1) + gry*gchi(i, 2) + grz*gchi(i, 3)
      end do
   end subroutine xc_grad_point_ext_body

   !
   ! The force density of one (point, function) pair, e_k = sum_s
   ! [P_sk X_s + Q1_sk RD_s + Q2_sk RDx_s], P and Q rebuilt from the point's
   ! coefficients exactly as the per-point tiles used to hold them:
   !
   !   restricted   P_k = w (v_rho dchi_k + 2 v_sigma (grad rho . H)_k)
   !                Q1_k = w 2 v_sigma dchi_k
   !   unrestricted P_ak = w (v_a dchi_k + 2 v_aa (grad rho_a . H)_k + v_ab (grad rho_b . H)_k)
   !                Q1_ak = w 2 v_aa dchi_k,  Q2_ak = w v_ab dchi_k   (b likewise)
   !
   ! with H packed xx, xy, xz, yy, yz, zz; an LDA keeps only the v_rho piece.
   !
   pure subroutine xc_grad_term(nspin, is_gga, pcv, gx, hx, x, rm, rx, e)
      !$acc routine seq
      integer, intent(in) :: nspin
      logical, intent(in) :: is_gga
      real(dp), intent(in) :: pcv(XC_GRAD_NPC), gx(3), hx(NHESS), x(2), rm(2), rx(2)
      real(dp), intent(out) :: e(3)
      real(dp) :: hda(3), hdb(3)
      real(dp) :: wg, pa, pb
      integer :: k

      wg = pcv(1)
      if (nspin == 1) then
         if (is_gga) then
            hda(1) = pcv(4)*hx(1) + pcv(5)*hx(2) + pcv(6)*hx(3)
            hda(2) = pcv(4)*hx(2) + pcv(5)*hx(4) + pcv(6)*hx(5)
            hda(3) = pcv(4)*hx(3) + pcv(5)*hx(5) + pcv(6)*hx(6)
            do k = 1, 3
               pa = wg*(pcv(2)*gx(k) + 2.0_dp*pcv(3)*hda(k))
               e(k) = pa*x(1) + wg*2.0_dp*pcv(3)*gx(k)*rm(1)
            end do
         else
            do k = 1, 3
               e(k) = wg*pcv(2)*gx(k)*x(1)
            end do
         end if
      else
         if (is_gga) then
            hda(1) = pcv(7)*hx(1) + pcv(8)*hx(2) + pcv(9)*hx(3)
            hda(2) = pcv(7)*hx(2) + pcv(8)*hx(4) + pcv(9)*hx(5)
            hda(3) = pcv(7)*hx(3) + pcv(8)*hx(5) + pcv(9)*hx(6)
            hdb(1) = pcv(10)*hx(1) + pcv(11)*hx(2) + pcv(XC_GRAD_NPC)*hx(3)
            hdb(2) = pcv(10)*hx(2) + pcv(11)*hx(4) + pcv(XC_GRAD_NPC)*hx(5)
            hdb(3) = pcv(10)*hx(3) + pcv(11)*hx(5) + pcv(XC_GRAD_NPC)*hx(6)
            do k = 1, 3
               pa = wg*(pcv(2)*gx(k) + 2.0_dp*pcv(4)*hda(k) + pcv(5)*hdb(k))
               pb = wg*(pcv(3)*gx(k) + 2.0_dp*pcv(6)*hdb(k) + pcv(5)*hda(k))
               e(k) = pa*x(1) + pb*x(2) &
                      + wg*2.0_dp*pcv(4)*gx(k)*rm(1) + wg*2.0_dp*pcv(6)*gx(k)*rm(2) &
                      + wg*pcv(5)*gx(k)*rx(1) + wg*pcv(5)*gx(k)*rx(2)
            end do
         else
            do k = 1, 3
               e(k) = wg*pcv(2)*gx(k)*x(1) + wg*pcv(3)*gx(k)*x(2)
            end do
         end if
      end if
   end subroutine xc_grad_term

   !
   ! The forces of one chunk: Pulay onto each function's atom and, with
   ! `gresp`, point motion onto each point's owner. The same two sums of
   ! xc_grad_term's force density, laid out for the device or the host.
   !
   ! DEVICE. The Pulay sum is one thread per (function, lane): the
   ! XC_GRAD_LANES threads of a function walk its tile column, contiguous,
   ! point fastest, at a stride of XC_GRAD_LANES points, so adjacent
   ! threads read adjacent points and every load coalesces. One thread per
   ! function put adjacent threads a padded column (npp) apart: w60/def2-SVP
   ! PBE on a V100, 470 ms -> 39 ms over the gradient's 36 launches, 16 and
   ! 64 lanes a little slower than 32. A second kernel sums each
   ! function's lanes and scatters, atomic: functions of different batches
   ! share atoms. The point-motion sum stays one thread per point
   ! (xc_grad_pointmotion), whose adjacent threads read adjacent points.
   !
   ! HOST. One iteration per stored batch walks each function's column once
   ! and adds every (point, function) density into both sums, the
   ! function's (Pulay) and the point's (point motion), so the density is
   ! formed once and from contiguous loads. The device's point-motion form,
   ! a point at a time across its functions a padded column apart, is what
   ! made the host's force step 3.8 s at w20/def2-SVP PBE (0.5 s now).
   ! Iterations write disjoint slices of fp, qp and eb (a batch's own
   ! functions and points), so only the scatter onto atoms is a plain loop,
   ! as xc_scatter's is. Every sum runs in the device's one-lane order.
   !
   subroutine xc_grad_forces(work, b, xg, gresp, grad)
      type(xc_grad_work_t), intent(inout) :: work
      !! The current chunk, its X (xt), RD (rd, rdx) and point
      !! coefficients (pc) formed; fp/qp/eb are its scratch.
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      logical, intent(in) :: gresp
      real(dp), intent(inout) :: grad(3, b%natm)
#ifdef TRC_DEVICE

      call xc_grad_pulay_lanes(work%nmu, work%nb, work%ord, work%mw_off, work%nspin, work%is_gga, work%p0, &
                               xg%b_off, work%c_off, work%npp, xg%nbatch, work%ldt, work%ldh, work%nhs, work%ldr, &
                               work%ldx, work%npc, work%gchi, work%hchi, work%xt, work%rd, work%rdx, work%ldp, &
                               work%pc, work%ldf, work%fp)
      call xc_grad_pulay_scatter(work%nmu, work%nb, work%ord, work%mw_off, xg%b_aooff, size(xg%b_ao), xg%b_ao, &
                                 b%nao, work%ao_atom, xg%nbatch, work%ldf, work%fp, b%natm, grad)
      if (gresp) then
         call xc_grad_pointmotion(work%np, work%p0, work%nspin, work%is_gga, xg%npts, xg%atom, xg%batch_of, &
                                  xg%b_off, xg%b_aooff, work%c_off, work%npp, work%active, xg%nbatch, work%ldt, &
                                  work%ldh, work%nhs, work%ldr, work%ldx, work%npc, work%gchi, work%hchi, work%xt, &
                                  work%rd, work%rdx, work%pc, work%ldp, b%natm, grad)
      end if
#else
      integer :: k, ib, mu, t, at, g, gg

      call xc_grad_force_batches(work%nb, work%ord, work%mw_off, work%nspin, work%is_gga, gresp, work%p0, &
                                 xg%b_off, work%c_off, work%npp, xg%nbatch, work%ldt, work%ldh, work%nhs, &
                                 work%ldr, work%ldx, work%npc, work%gchi, work%hchi, work%xt, work%rd, work%rdx, &
                                 work%ldp, work%pc, work%ldf, work%fp, work%ldq, work%qp, work%eb)
      do k = 1, work%nb
         ib = work%ord(k)
         do mu = 1, work%mw_off(k + 1) - work%mw_off(k)
            t = work%mw_off(k) + mu
            at = work%ao_atom(xg%b_ao(xg%b_aooff(ib) + mu - 1))
            grad(1, at) = grad(1, at) + work%fp(t, 1)
            grad(2, at) = grad(2, at) + work%fp(t, 2)
            grad(3, at) = grad(3, at) + work%fp(t, 3)
         end do
      end do
      if (gresp) then
         do g = 1, work%np
            gg = work%p0 + g - 1
            if (work%active(xg%batch_of(gg)) == 0) cycle
            at = xg%atom(gg)
            if (at > 0) then
               grad(1, at) = grad(1, at) + work%qp(g, 1)
               grad(2, at) = grad(2, at) + work%qp(g, 2)
               grad(3, at) = grad(3, at) + work%qp(g, 3)
            end if
         end do
      end if
#endif
   end subroutine xc_grad_forces

#ifdef TRC_DEVICE
   !
   ! DEVICE, the Pulay sum's first kernel: one thread per (function, lane)
   ! into fp's row t.
   !
   subroutine xc_grad_pulay_lanes(nmu, nb, ord, mw_off, nspin, is_gga, p0, b_off, c_off, npp, nbatch, ldt, ldh, &
                                  nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, ldp, pc, ldf, fp)
      integer, intent(in) :: nmu, nb, nspin, p0, nbatch, nhs, npc, ldp, ldf
      logical, intent(in) :: is_gga
      integer, intent(in) :: ord(nb), mw_off(nb + 1), b_off(nbatch + 1), c_off(nbatch + 1), npp(nbatch)
      integer(int64), intent(in) :: ldt, ldh, ldr, ldx
      real(dp), intent(in) :: gchi(ldt, 3), hchi(ldh, nhs), xt(ldt, nspin), rd(ldr, nspin), rdx(ldx, nspin)
      real(dp), intent(in) :: pc(ldp, npc)
      real(dp), intent(inout) :: fp(ldf, 3)
      real(dp) :: f(3)
      integer :: t

      do concurrent(t=1:nmu*XC_GRAD_LANES) local(f)
         call xc_grad_pulay_lane(t, nb, ord, mw_off, nspin, is_gga, p0, b_off, c_off, npp, nbatch, ldt, ldh, &
                                 nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, ldp, pc, f)
         fp(t, 1) = f(1)
         fp(t, 2) = f(2)
         fp(t, 3) = f(3)
      end do
   end subroutine xc_grad_pulay_lanes

   !
   ! DEVICE, the Pulay sum's second kernel: each function's lanes summed
   ! and added onto its atom, atomic -- functions of different batches
   ! share atoms.
   !
   subroutine xc_grad_pulay_scatter(nmu, nb, ord, mw_off, b_aooff, nbao, b_ao, nao, ao_atom, nbatch, ldf, fp, &
                                    natm, grad)
      integer, intent(in) :: nmu, nb, nbao, nao, nbatch, ldf, natm
      integer, intent(in) :: ord(nb), mw_off(nb + 1), b_aooff(nbatch + 1), b_ao(nbao), ao_atom(nao)
      real(dp), intent(in) :: fp(ldf, 3)
      real(dp), intent(inout) :: grad(3, natm)
      real(dp) :: f(3)
      integer :: m, at

      do concurrent(m=1:nmu) local(f, at)
         call xc_grad_pulay_sum(m, nb, ord, mw_off, b_aooff, nbao, b_ao, nao, ao_atom, nbatch, ldf, fp, f, at)
         !$acc atomic update
         grad(1, at) = grad(1, at) + f(1)
         !$acc atomic update
         grad(2, at) = grad(2, at) + f(2)
         !$acc atomic update
         grad(3, at) = grad(3, at) + f(3)
      end do
   end subroutine xc_grad_pulay_scatter

#else
   !
   ! HOST, every stored batch's force partials (xc_grad_force_batch): the
   ! iterations write disjoint rows of fp, qp and eb.
   !
   subroutine xc_grad_force_batches(nb, ord, mw_off, nspin, is_gga, gresp, p0, b_off, c_off, npp, nbatch, ldt, &
                                    ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, ldp, pc, ldf, fp, ldq, qp, eb)
      integer, intent(in) :: nb, nspin, p0, nbatch, nhs, npc, ldp, ldf, ldq
      logical, intent(in) :: is_gga, gresp
      integer, intent(in) :: ord(nb), mw_off(nb + 1), b_off(nbatch + 1), c_off(nbatch + 1), npp(nbatch)
      integer(int64), intent(in) :: ldt, ldh, ldr, ldx
      real(dp), intent(in) :: gchi(ldt, 3), hchi(ldh, nhs), xt(ldt, nspin), rd(ldr, nspin), rdx(ldx, nspin)
      real(dp), intent(in) :: pc(ldp, npc)
      real(dp), intent(inout) :: fp(ldf, 3), qp(ldq, 3), eb(ldq, 3)
      integer :: k

      do concurrent(k=1:nb)
         call xc_grad_force_batch(k, nb, ord, mw_off, nspin, is_gga, gresp, p0, b_off, c_off, npp, nbatch, &
                                  ldt, ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, ldp, pc, ldf, fp, ldq, qp, eb)
      end do
   end subroutine xc_grad_force_batches
#endif

   !
   ! One (function, lane) of the device's Pulay sum: lane l of function m
   ! (t = (m-1) XC_GRAD_LANES + l) sums the force density over the points
   ! l, l + XC_GRAD_LANES, ... of m's batch.
   !
   pure subroutine xc_grad_pulay_lane(t, nb, ord, mw_off, nspin, is_gga, p0, b_off, c_off, npp, nbatch, ldt, &
                                      ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, ldp, pc, f)
      !$acc routine seq
      integer, intent(in) :: t, nb, nspin, p0, nbatch, nhs, npc, ldp
      logical, intent(in) :: is_gga
      integer, intent(in) :: ord(nb), mw_off(nb + 1), b_off(nbatch + 1), c_off(nbatch + 1), npp(nbatch)
      integer(int64), intent(in) :: ldt, ldh, ldr, ldx
      real(dp), intent(in) :: gchi(ldt, 3), hchi(ldh, nhs), xt(ldt, nspin), rd(ldr, nspin), rdx(ldx, nspin)
      real(dp), intent(in) :: pc(ldp, npc)
      real(dp), intent(out) :: f(3)
      integer :: k, ib, mu, m, lane, npb, lg0
      integer(int64) :: col

      m = (t - 1)/XC_GRAD_LANES + 1
      lane = t - (m - 1)*XC_GRAD_LANES
      call grad_item_batch(m, nb, mw_off, k, mu)
      ib = ord(k)
      npb = b_off(ib + 1) - b_off(ib)
      col = int(c_off(ib), int64) + int(mu, int64)*int(npp(ib), int64)
      lg0 = b_off(ib) - p0
      call xc_grad_column_sum(lane, XC_GRAD_LANES, npb, col, lg0, nspin, is_gga, ldt, ldh, nhs, ldr, ldx, &
                              npc, gchi, hchi, xt, rd, rdx, ldp, pc, f)
   end subroutine xc_grad_pulay_lane

   !
   ! Function m's lanes summed, -2 x that, and the atom m sits on.
   !
   pure subroutine xc_grad_pulay_sum(m, nb, ord, mw_off, b_aooff, nbao, b_ao, nao, ao_atom, nbatch, ldf, fp, &
                                     f, at)
      !$acc routine seq
      integer, intent(in) :: m, nb, nbao, nao, nbatch, ldf
      integer, intent(in) :: ord(nb), mw_off(nb + 1), b_aooff(nbatch + 1), b_ao(nbao), ao_atom(nao)
      real(dp), intent(in) :: fp(ldf, 3)
      real(dp), intent(out) :: f(3)
      integer, intent(out) :: at
      integer :: k, mu, l, t0

      call grad_item_batch(m, nb, mw_off, k, mu)
      at = ao_atom(b_ao(b_aooff(ord(k)) + mu))
      t0 = (m - 1)*XC_GRAD_LANES
      f = 0.0_dp
      do l = 1, XC_GRAD_LANES
         f(1) = f(1) + fp(t0 + l, 1)
         f(2) = f(2) + fp(t0 + l, 2)
         f(3) = f(3) + fp(t0 + l, 3)
      end do
      f(1) = -2.0_dp*f(1)
      f(2) = -2.0_dp*f(2)
      f(3) = -2.0_dp*f(3)
   end subroutine xc_grad_pulay_sum

   !
   ! The force density of one function, summed over the points
   ! first, first + stride, ... (1-based, at most npb) of its tile column,
   ! which starts after element `col`; point gl's coefficients are pc's row
   ! lg0 + gl.
   !
   pure subroutine xc_grad_column_sum(first, stride, npb, col, lg0, nspin, is_gga, ldt, ldh, nhs, ldr, ldx, &
                                      npc, gchi, hchi, xt, rd, rdx, ldp, pc, f)
      !$acc routine seq
      integer, intent(in) :: first, stride, npb, lg0, nspin, nhs, npc, ldp
      integer(int64), intent(in) :: col, ldt, ldh, ldr, ldx
      logical, intent(in) :: is_gga
      real(dp), intent(in) :: gchi(ldt, 3), hchi(ldh, nhs), xt(ldt, nspin), rd(ldr, nspin), rdx(ldx, nspin)
      real(dp), intent(in) :: pc(ldp, npc)
      real(dp), intent(out) :: f(3)
      real(dp) :: pcv(XC_GRAD_NPC), gx(3), hx(NHESS), x(2), rm(2), rx(2), e(3)
      integer(int64) :: i
      integer :: gl

      pcv = 0.0_dp
      hx = 0.0_dp
      x = 0.0_dp
      rm = 0.0_dp
      rx = 0.0_dp
      f = 0.0_dp
      do gl = first, npb, stride
         i = col + gl
         call xc_grad_load(i, lg0 + gl, nspin, is_gga, ldt, ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, &
                           ldp, pc, pcv, gx, hx, x, rm, rx)
         call xc_grad_term(nspin, is_gga, pcv, gx, hx, x, rm, rx, e)
         f(1) = f(1) + e(1)
         f(2) = f(2) + e(2)
         f(3) = f(3) + e(3)
      end do
   end subroutine xc_grad_column_sum

   !
   ! The host's per-batch form (xc_grad_forces, HOST): batch k's functions
   ! into fp's rows mw_off(k) + 1.., -2 x their sums, and, with `gresp`,
   ! its points into qp's rows b_off - p0 + 1.., 2 x theirs. A column's
   ! densities go through eb's rows for those points first, so the loop
   ! that forms them carries no sum and vectorises.
   !
   pure subroutine xc_grad_force_batch(k, nb, ord, mw_off, nspin, is_gga, gresp, p0, b_off, c_off, npp, nbatch, &
                                       ldt, ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, ldp, pc, ldf, fp, &
                                       ldq, qp, eb)
      !$acc routine seq
      integer, intent(in) :: k, nb, nspin, p0, nbatch, nhs, npc, ldp, ldf, ldq
      logical, intent(in) :: is_gga, gresp
      integer, intent(in) :: ord(nb), mw_off(nb + 1), b_off(nbatch + 1), c_off(nbatch + 1), npp(nbatch)
      integer(int64), intent(in) :: ldt, ldh, ldr, ldx
      real(dp), intent(in) :: gchi(ldt, 3), hchi(ldh, nhs), xt(ldt, nspin), rd(ldr, nspin), rdx(ldx, nspin)
      real(dp), intent(in) :: pc(ldp, npc)
      real(dp), intent(inout) :: fp(ldf, 3), qp(ldq, 3), eb(ldq, 3)
      real(dp) :: f(3)
      integer(int64) :: col
      integer :: ib, mu, nloc, npb, lg0, gl, lg, t

      ib = ord(k)
      nloc = mw_off(k + 1) - mw_off(k)
      if (nloc == 0) return
      npb = b_off(ib + 1) - b_off(ib)
      lg0 = b_off(ib) - p0
      if (gresp) then
         do gl = 1, npb
            qp(lg0 + gl, 1) = 0.0_dp
            qp(lg0 + gl, 2) = 0.0_dp
            qp(lg0 + gl, 3) = 0.0_dp
         end do
      end if
      do mu = 1, nloc
         col = int(c_off(ib), int64) + int(mu - 1, int64)*int(npp(ib), int64)
         ! Literal nspin, is_gga and npc, so the inlined density has no
         ! branches left in it.
         if (nspin == 1 .and. is_gga) then
            call xc_grad_column(1, .true., 6, npb, col, lg0, ldt, ldh, nhs, ldr, ldx, gchi, hchi, xt, rd, rdx, &
                                ldp, pc, ldq, eb)
         else if (nspin == 1) then
            call xc_grad_column(1, .false., 6, npb, col, lg0, ldt, ldh, nhs, ldr, ldx, gchi, hchi, xt, rd, rdx, &
                                ldp, pc, ldq, eb)
         else if (is_gga) then
            call xc_grad_column(2, .true., XC_GRAD_NPC, npb, col, lg0, ldt, ldh, nhs, ldr, ldx, gchi, hchi, xt, &
                                rd, rdx, ldp, pc, ldq, eb)
         else
            call xc_grad_column(2, .false., XC_GRAD_NPC, npb, col, lg0, ldt, ldh, nhs, ldr, ldx, gchi, hchi, xt, &
                                rd, rdx, ldp, pc, ldq, eb)
         end if
         f = 0.0_dp
         do gl = 1, npb
            lg = lg0 + gl
            f(1) = f(1) + eb(lg, 1)
            f(2) = f(2) + eb(lg, 2)
            f(3) = f(3) + eb(lg, 3)
         end do
         t = mw_off(k) + mu
         fp(t, 1) = -2.0_dp*f(1)
         fp(t, 2) = -2.0_dp*f(2)
         fp(t, 3) = -2.0_dp*f(3)
         if (gresp) then
            do gl = 1, npb
               lg = lg0 + gl
               qp(lg, 1) = qp(lg, 1) + eb(lg, 1)
               qp(lg, 2) = qp(lg, 2) + eb(lg, 2)
               qp(lg, 3) = qp(lg, 3) + eb(lg, 3)
            end do
         end if
      end do
      if (gresp) then
         do gl = 1, npb
            qp(lg0 + gl, 1) = 2.0_dp*qp(lg0 + gl, 1)
            qp(lg0 + gl, 2) = 2.0_dp*qp(lg0 + gl, 2)
            qp(lg0 + gl, 3) = 2.0_dp*qp(lg0 + gl, 3)
         end do
      end if
   end subroutine xc_grad_force_batch

   ! One function's force density at each of its batch's npb points, into
   ! eb's rows lg0 + 1...
   pure subroutine xc_grad_column(nspin, is_gga, npc, npb, col, lg0, ldt, ldh, nhs, ldr, ldx, gchi, hchi, xt, &
                                  rd, rdx, ldp, pc, ldq, eb)
      !$acc routine seq
      integer, intent(in) :: nspin, npc, npb, lg0, nhs, ldp, ldq
      logical, intent(in) :: is_gga
      integer(int64), intent(in) :: col, ldt, ldh, ldr, ldx
      real(dp), intent(in) :: gchi(ldt, 3), hchi(ldh, nhs), xt(ldt, nspin), rd(ldr, nspin), rdx(ldx, nspin)
      real(dp), intent(in) :: pc(ldp, npc)
      real(dp), intent(inout) :: eb(ldq, 3)
      real(dp) :: pcv(XC_GRAD_NPC), gx(3), hx(NHESS), x(2), rm(2), rx(2), e(3)
      real(dp) :: wg, h1, h2, h3, gk
      integer(int64) :: i
      integer :: gl, lg, kk

      if (nspin == 1) then
         ! xc_grad_term's restricted lines, term for term in its order.
         if (is_gga) then
            do gl = 1, npb
               i = col + gl
               lg = lg0 + gl
               wg = pc(lg, 1)
               h1 = pc(lg, 4)*hchi(i, 1) + pc(lg, 5)*hchi(i, 2) + pc(lg, 6)*hchi(i, 3)
               h2 = pc(lg, 4)*hchi(i, 2) + pc(lg, 5)*hchi(i, 4) + pc(lg, 6)*hchi(i, 5)
               h3 = pc(lg, 4)*hchi(i, 3) + pc(lg, 5)*hchi(i, 5) + pc(lg, 6)*hchi(i, 6)
               eb(lg, 1) = wg*(pc(lg, 2)*gchi(i, 1) + 2.0_dp*pc(lg, 3)*h1)*xt(i, 1) &
                           + wg*2.0_dp*pc(lg, 3)*gchi(i, 1)*rd(i, 1)
               eb(lg, 2) = wg*(pc(lg, 2)*gchi(i, 2) + 2.0_dp*pc(lg, 3)*h2)*xt(i, 1) &
                           + wg*2.0_dp*pc(lg, 3)*gchi(i, 2)*rd(i, 1)
               eb(lg, 3) = wg*(pc(lg, 2)*gchi(i, 3) + 2.0_dp*pc(lg, 3)*h3)*xt(i, 1) &
                           + wg*2.0_dp*pc(lg, 3)*gchi(i, 3)*rd(i, 1)
            end do
         else
            do gl = 1, npb
               i = col + gl
               lg = lg0 + gl
               do kk = 1, 3
                  gk = gchi(i, kk)
                  eb(lg, kk) = pc(lg, 1)*pc(lg, 2)*gk*xt(i, 1)
               end do
            end do
         end if
         return
      end if
      pcv = 0.0_dp
      hx = 0.0_dp
      x = 0.0_dp
      rm = 0.0_dp
      rx = 0.0_dp
      do gl = 1, npb
         call xc_grad_load(col + gl, lg0 + gl, nspin, is_gga, ldt, ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, &
                           rdx, ldp, pc, pcv, gx, hx, x, rm, rx)
         call xc_grad_term(nspin, is_gga, pcv, gx, hx, x, rm, rx, e)
         eb(lg0 + gl, 1) = e(1)
         eb(lg0 + gl, 2) = e(2)
         eb(lg0 + gl, 3) = e(3)
      end do
   end subroutine xc_grad_column

   !
   ! What xc_grad_term reads for tile element i, point row lg of pc.
   !
   pure subroutine xc_grad_load(i, lg, nspin, is_gga, ldt, ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, &
                                ldp, pc, pcv, gx, hx, x, rm, rx)
      !$acc routine seq
      integer(int64), intent(in) :: i, ldt, ldh, ldr, ldx
      integer, intent(in) :: lg, nspin, nhs, npc, ldp
      logical, intent(in) :: is_gga
      real(dp), intent(in) :: gchi(ldt, 3), hchi(ldh, nhs), xt(ldt, nspin), rd(ldr, nspin), rdx(ldx, nspin)
      real(dp), intent(in) :: pc(ldp, npc)
      real(dp), intent(inout) :: pcv(XC_GRAD_NPC), gx(3), hx(NHESS), x(2), rm(2), rx(2)
      integer :: j, s

      do j = 1, npc
         pcv(j) = pc(lg, j)
      end do
      gx(1) = gchi(i, 1)
      gx(2) = gchi(i, 2)
      gx(3) = gchi(i, 3)
      if (is_gga) then
         do j = 1, NHESS
            hx(j) = hchi(i, j)
         end do
      end if
      do s = 1, nspin
         x(s) = xt(i, s)
         if (is_gga) rm(s) = rd(i, s)
         if (is_gga .and. nspin == 2) rx(s) = rdx(i, s)
      end do
   end subroutine xc_grad_load

   !
   ! GRID RESPONSE, point-motion piece (this module's header, A): the same
   ! force density summed over the point's local functions instead, +2 x
   ! that onto the atom that OWNS the point. One thread per point, atomic.
   !
   subroutine xc_grad_pointmotion(np, p0, nspin, is_gga, npts, g_atom, batch_of, b_off, b_aooff, c_off, npp, &
                                  active, nbatch, ldt, ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, pc, ldp, &
                                  natm, grad)
      integer, intent(in) :: np, p0, nspin, npts, nbatch, npc, ldp, natm, nhs
      logical, intent(in) :: is_gga
      integer, intent(in) :: g_atom(npts), batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1)
      integer, intent(in) :: c_off(nbatch + 1), npp(nbatch), active(nbatch)
      integer(int64), intent(in) :: ldt, ldh, ldr, ldx
      real(dp), intent(in) :: gchi(ldt, 3), hchi(ldh, nhs), xt(ldt, nspin), rd(ldr, nspin), rdx(ldx, nspin)
      real(dp), intent(in) :: pc(ldp, npc)
      real(dp), intent(inout) :: grad(3, natm)
      integer :: g, at
      real(dp) :: q(3)
#ifdef TRC_DEVICE
      do concurrent(g=1:np) local(q, at)
         call xc_grad_pointmotion_body(g, p0, nspin, is_gga, npts, g_atom, batch_of, b_off, b_aooff, c_off, npp, &
                                       active, nbatch, ldt, ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, pc, ldp, &
                                       q, at)
         if (at > 0) then
            !$acc atomic update
            grad(1, at) = grad(1, at) + q(1)
            !$acc atomic update
            grad(2, at) = grad(2, at) + q(2)
            !$acc atomic update
            grad(3, at) = grad(3, at) + q(3)
         end if
      end do
#else
      do g = 1, np
         call xc_grad_pointmotion_body(g, p0, nspin, is_gga, npts, g_atom, batch_of, b_off, b_aooff, c_off, npp, &
                                       active, nbatch, ldt, ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, rdx, pc, ldp, &
                                       q, at)
         if (at > 0) then
            grad(1, at) = grad(1, at) + q(1)
            grad(2, at) = grad(2, at) + q(2)
            grad(3, at) = grad(3, at) + q(3)
         end if
      end do
#endif
   end subroutine xc_grad_pointmotion

   pure subroutine xc_grad_pointmotion_body(g, p0, nspin, is_gga, npts, g_atom, batch_of, b_off, b_aooff, c_off, &
                                            npp, active, nbatch, ldt, ldh, nhs, ldr, ldx, npc, gchi, hchi, xt, rd, &
                                            rdx, pc, ldp, q, at)
      !$acc routine seq
      integer, intent(in) :: g, p0, nspin, npts, nbatch, npc, ldp, nhs
      logical, intent(in) :: is_gga
      integer, intent(in) :: g_atom(npts), batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1)
      integer, intent(in) :: c_off(nbatch + 1), npp(nbatch), active(nbatch)
      integer(int64), intent(in) :: ldt, ldh, ldr, ldx
      real(dp), intent(in) :: gchi(ldt, 3), hchi(ldh, nhs), xt(ldt, nspin), rd(ldr, nspin), rdx(ldx, nspin)
      real(dp), intent(in) :: pc(ldp, npc)
      real(dp), intent(out) :: q(3)
      integer, intent(out) :: at
      integer :: gg, ib, nloc, ld, u, s, j
      integer(int64) :: off, i
      real(dp) :: pcv(XC_GRAD_NPC), gx(3), hx(NHESS), x(2), rm(2), rx(2), e(3)

      q = 0.0_dp
      at = 0
      gg = p0 + g - 1
      ib = batch_of(gg)
      if (active(ib) == 0) return
      at = g_atom(gg)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = npp(ib)
      off = int(c_off(ib), int64) + int(gg - b_off(ib), int64) + 1_int64
      pcv = 0.0_dp
      do j = 1, npc
         pcv(j) = pc(g, j)
      end do
      hx = 0.0_dp
      x = 0.0_dp
      rm = 0.0_dp
      rx = 0.0_dp
      do u = 1, nloc
         i = off + int(u - 1, int64)*int(ld, int64)
         gx(1) = gchi(i, 1)
         gx(2) = gchi(i, 2)
         gx(3) = gchi(i, 3)
         if (is_gga) then
            do j = 1, NHESS
               hx(j) = hchi(i, j)
            end do
         end if
         do s = 1, nspin
            x(s) = xt(i, s)
            if (is_gga) rm(s) = rd(i, s)
            if (is_gga .and. nspin == 2) rx(s) = rdx(i, s)
         end do
         call xc_grad_term(nspin, is_gga, pcv, gx, hx, x, rm, rx, e)
         q(1) = q(1) + e(1)
         q(2) = q(2) + e(2)
         q(3) = q(3) + e(3)
      end do
      q(1) = 2.0_dp*q(1)
      q(2) = 2.0_dp*q(2)
      q(3) = 2.0_dp*q(3)
   end subroutine xc_grad_pointmotion_body


end module trc_xc_grad
