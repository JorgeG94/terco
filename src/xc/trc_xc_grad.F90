!
! The exchange-correlation contribution to the nuclear gradient, LDA and
! GGA, restricted and unrestricted. Added into an existing `grad(3,natm)`,
! the same convention trc_1e_gradient and nuclear_repulsion_gradient use
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
! xc_grad_point_gga_body for the channel-by-channel expansion. Checked,
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
!      + the SAME contraction the Pulay force tiles (pall/q1all/q2all/rall)
!      already build, but summed over EVERY local function of the batch
!      (not just the ones on one atom) and landing on the POINT's owner
!      (not the function's), with the OPPOSITE sign:
!
!          dE_xc/dR_{B,k} += sum_{g: owner(g)=B} w_g [ vrho_g (drho_g/dx_k)
!                             + vsigma_g (dsigma_g/dx_k) ]
!
!      dsigma_g/dx_k = 2 (H_rho(g) . grad_rho_g)_k needs the density
!      HESSIAN at the point, which trc_xc_grad already collocates
!      (`hchi`) for the Pulay term's own vsigma piece -- xc_grad_point_lda/
!      gga_body compute this directly: the LDA-like piece (vrho) is cheap
!      enough to fold into the SAME point kernel that builds `pall`
!      (it only needs `pall` and `xmat`, both already there); the GGA
!      vsigma piece needs `rmat`/`rmatx` -- `rall` contracted with the
!      density matrix the same way `xmat` is chi's -- which is a
!      whole-batch GEMM, so it is added by a follow-up kernel
!      (xc_grad_pointmotion_gga) once xc_grad_point_gga has filled `rall`
!      for every point. Scattered by xc_grad_scatter_points, keyed by
!      `xg%atom` (the point's owner), not `ao_atom`.
!
!   B. PARTITION WEIGHT, explicit nuclear dependence, AT FIXED r_g.
!      `trc_dft_partition.F90`'s `becke_partition_derivatives` (already in
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
!      owner, for FREE out of the same `dweights` array B already
!      computed -- no new derivative, no new routine. This is the piece a
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
! BATCHING -- deliberately SIMPLER than the energy path's. trc_xc.F90's
! xc_drive pads every batch's tile to a common shape and groups batches of
! equal padded shape into one strided-batched cuBLAS call, because it runs
! every SCF iteration and the grouping is what makes that affordable. This
! module runs once, after convergence, so it processes one spatial batch
! at a time -- collocate, gather, a handful of exact-sized GEMMs, scatter
! -- with no padding and no cross-batch grouping. Still cuBLAS (`trc_linalg`),
! still batched (one pass per `trc_xc_grid_t` batch, same screening), just
! not padding-grouped; see this branch's own report for the measured cost
! against the two-electron gradient. Grouping batches the way the energy
! path does is the natural next step if this ever dominates a gradient's
! wall time.
!
! CARTESIAN ONLY. XC is Cartesian throughout terco (the SCF folds to
! spherical); `b` and `dens` here are exactly what trc_scf_grad.F90 hands
! in, with no spherical folding of its own.
!
! GPU: CONFIRMED WRONG, NOT MERELY UNTESTED. check_full_grad_fd (the
! existing one/two-electron gradient) and scf_rks (the existing XC
! ENERGY path, same atomic-scatter idiom as this module's own
! xc_grad_scatter) both pass on this exact GPU build, which rules out the
! basis/pairlist/eri/trc_linalg/XC-energy machinery this module reuses.
! What is NEW here -- and the only thing left to suspect -- is this
! module's own per-batch device lifetime: unlike xc_drive, which
! allocates its device tiles ONCE outside the batch loop and reuses them
! (see BATCHING above), xc_grad_one_batch allocates (`!$acc enter data
! create`) and frees (`!$acc exit data delete`) a fresh set of device
! buffers on EVERY batch. check_xc_grad_fixedgrid_fd and
! check_xc_grad_vs_pyscf both fail the same way on GPU 2, deterministically
! (repeated runs give the identical worst|diff|, not a flaky range), with
! the error growing with how much per-batch work there is (UKS GGA worse
! than RKS GGA worse than RKS LDA) -- consistent with a stale or
! overlapping device allocation across batches rather than a sign/index
! mistake, which CPU's bitwise-identical formula code rules out anyway.
! NOT fixed this round: restructuring the per-batch allocation to the
! energy path's hoist-and-reuse convention is a real change to this
! module's data flow, not a one-line fix, and is left for the GPU
! follow-up this header's BATCHING section already calls out as the next
! step. Use the CPU build for this module until then.
!
module trc_xc_grad
   use pic_types, only: dp
   use trc_api, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_eval_point, xc_eval_point_polar, &
                                 XC_FAMILY_LDA, XC_FAMILY_GGA, XC_MAX_KERNELS
   use trc_collocation, only: shell_collocate, shell_collocate2, NCART_MAX, NHESS
   use trc_linalg, only: trc_linalg_t
   use trc_dft_grid, only: dft_grid_t, build_dft_grid, DEFAULT_GRID_LEVEL
   use trc_dft_partition, only: becke_partition_derivatives
   use trc_error, only: error_t
   use trc_step_table, only: step_table_t
   implicit none
   private

   public :: xc_gradient_add, XC_GRAD_RHO_TOL
   public :: XC_GRAD_NSTEP, xc_grad_step_names
   public :: xc_grad_collocate1, xc_grad_collocate2, xc_grad_gather
   public :: xc_grad_scatter, xc_grad_scatter_points
   ! Exported for trc_vv10_grad: the AO-value/gradient/Hessian collocation,
   ! the density gather onto a batch's local basis, and the two atom
   ! scatters are generic per-batch plumbing with nothing semilocal-XC-
   ! specific about them (no `kid`/`coef`/`xc_eval_point` in sight) -- VV10's
   ! gradient reuses them rather than keeping a second, parallel copy.

   integer, parameter :: XC_GRAD_NSTEP = 10
   !! Steps of the timing table `xc_gradient_add` fills (see
   !! `xc_grad_step_names`); a caller's step_table_t is initialised with them.
   integer, parameter :: ST_GRID = 1, ST_SCREEN = 2, ST_ALLOC = 3, ST_COLLOC = 4, ST_GATHER = 5, &
                         ST_GEMM = 6, ST_POINTS = 7, ST_SCATTER = 8, ST_PDERIV = 9, ST_RESP = 10

   real(dp), parameter :: XC_GRAD_RHO_TOL = 1.0e-12_dp
   !! Same role and value as trc_xc.F90's XC_RHO_TOL: a batch whose bound on
   !! rho is below this contributes nothing and is skipped before any
   !! device work for it.

contains

   pure function xc_grad_step_names() result(names)
      !! The step names of `xc_gradient_add`'s timing table, in order.
      character(len=16) :: names(XC_GRAD_NSTEP)
      names = [character(len=16) :: "grid", "screen", "alloc", "collocate", "gather", &
               "gemm", "points", "scatter", "partition deriv", "response"]
   end function xc_grad_step_names

   !
   ! The one entry point. Builds its own DFT quadrature at `b`'s geometry
   ! (trc_scf_run keeps neither its grid nor its functional alive past
   ! convergence, the same reason trc_scf_gradient rebuilds the pair list
   ! and the ERI engine), refuses what it cannot do, and ADDS the result
   ! into `grad`.
   !
   subroutine xc_gradient_add(b, lp, func, nspin, dens, grad, grid_level, grid_max_pts, rho_tol, &
                              grid_response, timing)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_xc_functional_t), intent(in) :: func
      integer, intent(in) :: nspin
      real(dp), intent(in) :: dens(b%nao, b%nao, nspin)
      !! The converged density, per spin channel. Present on both host and
      !! device: the GEMM pipeline below reads the device copy, the
      !! per-batch screen reads the host one (trc_scf_gradient's own
      !! convention -- see its header -- keeps both valid at this point).
      real(dp), intent(inout) :: grad(3, b%natm)
      !! ADDED into; device-resident for the whole call, exactly as
      !! trc_1e_gradient and nuclear_repulsion_gradient require of their
      !! own `grad` argument.
      integer, intent(in), optional :: grid_level, grid_max_pts
      real(dp), intent(in), optional :: rho_tol
      logical, intent(in), optional :: grid_response
      !! ON by default. See this module's header, GRID RESPONSE. OFF
      !! recovers the Pulay-only gradient the base of this branch shipped
      !! (what check_xc_grad_fixedgrid_fd's fixed-grid finite difference
      !! isolates, and what pyscf's own `grid_response=False` default
      !! compares against in check_xc_grad_vs_pyscf).
      type(step_table_t), intent(inout), optional, target :: timing
      !! Device-synced time per step, ADDED into; initialised by the caller
      !! with `xc_grad_step_names()`. Absent, a private table is kept and
      !! dropped.

      type(dft_grid_t) :: grid
      type(step_table_t), target :: tbl_own
      type(step_table_t), pointer :: tbl
      real(dp) :: tc
      type(trc_xc_grid_t) :: xg
      type(error_t) :: err
      logical :: is_gga, gresp
      integer :: glevel, gmaxpts, ish, i
      real(dp) :: rtol
      integer, allocatable :: ao_atom(:), at_z(:)
      real(dp), allocatable :: grad_extra(:, :)

      ! --- refuse what this module cannot do, with a reason -------------
      !
      ! Meta-GGA and double hybrids are not merely unimplemented here: on
      ! this branch `trc_xc_functional_t` has no meta-GGA family and no
      ! MP2-fraction field at all (see its own header), so neither can
      ! reach a converged SCF result in the first place -- there is no
      ! flag left to check for a double hybrid. The check below is for
      ! what CAN reach here that USED to be missing a kernel: VV10, whose
      ! nonlocal correlation gradient is now trc_vv10_grad%vv10_gradient_add's
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
      is_gga = func%family == XC_FAMILY_GGA

      glevel = DEFAULT_GRID_LEVEL
      if (present(grid_level)) glevel = grid_level
      gmaxpts = 512
      if (present(grid_max_pts)) gmaxpts = grid_max_pts
      rtol = XC_GRAD_RHO_TOL
      if (present(rho_tol)) rtol = rho_tol
      gresp = .true.
      if (present(grid_response)) gresp = grid_response

      if (present(timing)) then
         tbl => timing
      else
         call tbl_own%init("xc grad", xc_grad_step_names())
         tbl => tbl_own
      end if

      call tbl%step_on(lp, ST_GRID, tc)
      allocate (at_z(b%natm), source=nint(b%at_z))
      call build_dft_grid(b%at_r, at_z, grid, err, level=glevel)
      if (err%has_error()) error stop "trc_xc_grad: "//trim(err%get_message())
      ! `atom`/`qw` are only ever read by the grid-response term below
      ! (xg%build leaves both unallocated, and every other caller of
      ! xcgrid_build never asks for them); built unconditionally anyway --
      ! they cost one more integer and one more real(dp) per point, next to
      ! the GEMMs a surviving batch pays for.
      call xg%build(grid%n_points, grid%coords, grid%weights, b, max_pts=gmaxpts, &
                    atom=grid%atom, qw=grid%quad_weights)
      call xg%to_device()
      call tbl%step_off(lp, ST_GRID, tc)

      if (xg%npts > 0) then
         ! AO -> atom, built once from the shell -> atom map `basis_build`
         ! already produced (sh_at), the same map the two-electron and
         ! one-electron gradients use.
         allocate (ao_atom(b%nao))
         do ish = 1, b%nshell
            do i = 0, (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2 - 1
               ao_atom(b%sh_ao(ish) + i) = b%sh_at(ish)
            end do
         end do
         !$acc enter data copyin(ao_atom)

         ! The weight-response and point-motion-of-the-weight terms (this
         ! module's header, GRID RESPONSE) are accumulated HOST-side, into
         ! their own array: trc_dft_partition's becke_partition_derivatives
         ! is OpenMP, not OpenACC (it was already in the tree, unused, this
         ! branch is its first caller), and running it per batch is cheap
         ! enough (O(active_points x natm^2), natm tiny on every case this
         ! round validates) that porting it to the device is not this
         ! round's problem. Zeroed and merged into the device-resident
         ! `grad` once, after the drive loop, rather than after every batch.
         allocate (grad_extra(3, b%natm))
         grad_extra = 0.0_dp

         call xc_grad_drive(b, lp, is_gga, nspin, func%nk, func%kid, func%coef, dens, xg, rtol, &
                            ao_atom, grad, gresp, b%at_r, at_z, grid%scheme, grid%adjust, grad_extra, tbl)

         if (gresp) then
            !$acc update self(grad)
            grad = grad + grad_extra
            !$acc update device(grad)
         end if
         deallocate (grad_extra)

         !$acc exit data delete(ao_atom)
         deallocate (ao_atom)
      end if

      call xg%release()
      call grid%destroy()
      deallocate (at_z)
   end subroutine xc_gradient_add

   !
   ! One pass over every batch of the grid: screen, then collocate + gather
   ! + GEMM + scatter for the ones that survive. See this module's header
   ! for why this is per-batch rather than the energy path's padded,
   ! shape-grouped chunks.
   !
   subroutine xc_grad_drive(b, lp, is_gga, nspin, nk, kid, coef, dens, xg, rtol, ao_atom, grad, &
                            gresp, at_r, at_z, scheme, adjust, grad_extra, tbl)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      logical, intent(in) :: is_gga
      integer, intent(in) :: nspin, nk
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: dens(b%nao, b%nao, nspin)
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: rtol
      integer, intent(in) :: ao_atom(b%nao)
      real(dp), intent(inout) :: grad(3, b%natm)
      logical, intent(in) :: gresp
      real(dp), intent(in) :: at_r(3, b%natm)
      integer, intent(in) :: at_z(b%natm)
      integer, intent(in) :: scheme, adjust
      real(dp), intent(inout) :: grad_extra(3, b%natm)
      type(step_table_t), intent(inout) :: tbl

      integer :: ib, npb, nloc
      logical :: act
      real(dp) :: tc

      do ib = 1, xg%nbatch
         npb = xg%b_off(ib + 1) - xg%b_off(ib)
         nloc = xg%b_aooff(ib + 1) - xg%b_aooff(ib)
         if (npb == 0 .or. nloc == 0) cycle
         call tbl%step_on(lp, ST_SCREEN, tc)
         act = xc_grad_batch_active(ib, xg, b%nao, nspin, dens, rtol)
         call tbl%step_off(lp, ST_SCREEN, tc)
         if (.not. act) cycle
         call xc_grad_one_batch(b, lp, is_gga, nspin, nk, kid, coef, dens, xg, ib, npb, nloc, &
                                ao_atom, grad, gresp, at_r, at_z, scheme, adjust, grad_extra, tbl)
      end do
   end subroutine xc_grad_drive

   !
   ! Host-side density screen, bounding rho on a batch through the same
   ! per-function maxima (`xg%b_amax`) trc_xc.F90's xc_screen_body uses --
   ! plain Fortran, not a device kernel, because `dens`'s host copy already
   ! holds the right numbers (trc_scf_gradient's `!$acc update self(dens)`,
   ! before this is ever called) and nloc is small enough that this costs
   ! nothing next to the GEMMs a surviving batch pays for.
   !
   logical function xc_grad_batch_active(ib, xg, nao, nspin, dens, rtol) result(active)
      integer, intent(in) :: ib, nao, nspin
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: dens(nao, nao, nspin)
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
               su = su + abs(dens(av, au, s))*xg%b_amax(v)
            end do
            bound = bound + su*xg%b_amax(u)
         end do
      end do
      active = bound > rtol
   end function xc_grad_batch_active

   !
   ! One batch: collocate (order 1, or order 1+2 for a GGA), gather D onto
   ! the local basis, X = chi D, the per-point functional and force
   ! tiles, the force GEMMs, and the atom scatter. Every array here is
   ! sized to this batch exactly -- see the module header on why there is
   ! no padding to match.
   !
   subroutine xc_grad_one_batch(b, lp, is_gga, nspin, nk, kid, coef, dens, xg, ib, npb, nloc, &
                                ao_atom, grad, gresp, at_r, at_z, scheme, adjust, grad_extra, tbl)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      logical, intent(in) :: is_gga
      integer, intent(in) :: nspin, nk
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: dens(b%nao, b%nao, nspin)
      type(trc_xc_grid_t), intent(in) :: xg
      integer, intent(in) :: ib, npb, nloc
      integer, intent(in) :: ao_atom(b%nao)
      real(dp), intent(inout) :: grad(3, b%natm)
      logical, intent(in) :: gresp
      real(dp), intent(in) :: at_r(3, b%natm)
      integer, intent(in) :: at_z(b%natm)
      integer, intent(in) :: scheme, adjust
      real(dp), intent(inout) :: grad_extra(3, b%natm)
      type(step_table_t), intent(inout) :: tbl

      real(dp), allocatable :: chi(:, :), gchi(:, :, :), hchi(:, :, :)
      real(dp), allocatable :: dloc(:, :, :), xmat(:, :, :)
      real(dp), allocatable :: pall(:, :, :), q1all(:, :, :), q2all(:, :, :), rall(:, :, :)
      real(dp), allocatable :: tall(:, :, :)
      real(dp), allocatable :: rmat(:, :, :), rmatx(:, :, :)
      real(dp), allocatable :: fxc(:), qmove(:, :)
      real(dp), allocatable :: dweights(:, :, :)
      type(error_t) :: err
      integer :: s, g, a, own, g0
      real(dp) :: contrib, tc

      call tbl%step_on(lp, ST_ALLOC, tc)
      allocate (chi(npb, nloc), gchi(npb, nloc, 3))
      allocate (dloc(nloc, nloc, nspin), xmat(npb, nloc, nspin))
      allocate (pall(npb, 3*nloc, nspin), tall(3*nloc, nloc, nspin))
      allocate (fxc(npb), qmove(npb, 3))
      !$acc enter data create(chi, gchi, dloc, xmat, pall, tall, fxc, qmove)

      if (is_gga) then
         allocate (hchi(npb, nloc, NHESS))
         allocate (q1all(npb, 3*nloc, nspin), q2all(npb, 3*nloc, nspin), rall(npb, nloc, nspin))
         allocate (rmat(npb, nloc, nspin), rmatx(npb, nloc, nspin))
         !$acc enter data create(hchi, q1all, q2all, rall, rmat, rmatx)
      end if
      call tbl%step_off(lp, ST_ALLOC, tc)

      call tbl%step_on(lp, ST_COLLOC, tc)
      if (is_gga) then
         call xc_grad_collocate2(b%nshell, b%maxnp, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                                 xg%npts, xg%r, xg%nbatch, xg%b_off, xg%b_shoff, size(xg%b_sh), &
                                 xg%b_sh, xg%b_shao, ib, npb, nloc, chi, gchi, hchi)
      else
         call xc_grad_collocate1(b%nshell, b%maxnp, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                                 xg%npts, xg%r, xg%nbatch, xg%b_off, xg%b_shoff, size(xg%b_sh), &
                                 xg%b_sh, xg%b_shao, ib, npb, nloc, chi, gchi)
      end if
      call tbl%step_off(lp, ST_COLLOC, tc)

      call tbl%step_on(lp, ST_GATHER, tc)
      call xc_grad_gather(b%nao, size(xg%b_ao), xg%b_ao, xg%b_aooff(ib), nloc, nspin, dens, dloc)
      call tbl%step_off(lp, ST_GATHER, tc)

      call tbl%step_on(lp, ST_GEMM, tc)
      do s = 1, nspin
         call lp%gemm_mat(chi, dloc(:, :, s), xmat(:, :, s))
      end do
      call tbl%step_off(lp, ST_GEMM, tc)

      if (is_gga) then
         call tbl%step_on(lp, ST_POINTS, tc)
         call xc_grad_point_gga(npb, xg%npts, xg%w, xg%b_off(ib), nloc, nspin, nk, kid, coef, &
                                chi, gchi, hchi, xmat, pall, q1all, q2all, rall, fxc, qmove)
         call tbl%step_off(lp, ST_POINTS, tc)
         ! GRID RESPONSE, point-motion piece -- see this module's header.
         ! `rmat`/`rmatx` need every point's `rall` (the GEMM below
         ! contracts across the whole batch), so they cannot be fused into
         ! the point kernel above the way the LDA-like piece of `qmove`
         ! is; `xc_grad_pointmotion_gga` ADDS its own contribution into the
         ! `qmove` that kernel already seeded.
         if (gresp) then
            call tbl%step_on(lp, ST_GEMM, tc)
            do s = 1, nspin
               call lp%gemm_mat(rall(:, :, s), dloc(:, :, s), rmat(:, :, s))
            end do
            if (nspin == 2) then
               call lp%gemm_mat(rall(:, :, 2), dloc(:, :, 1), rmatx(:, :, 1))
               call lp%gemm_mat(rall(:, :, 1), dloc(:, :, 2), rmatx(:, :, 2))
            end if
            call tbl%step_off(lp, ST_GEMM, tc)
            call tbl%step_on(lp, ST_POINTS, tc)
            call xc_grad_pointmotion_gga(npb, nloc, nspin, q1all, q2all, rmat, rmatx, qmove)
            call tbl%step_off(lp, ST_POINTS, tc)
         end if
      else
         call tbl%step_on(lp, ST_POINTS, tc)
         call xc_grad_point_lda(npb, xg%npts, xg%w, xg%b_off(ib), nloc, nspin, nk, kid, coef, &
                                chi, gchi, xmat, pall, fxc, qmove)
         call tbl%step_off(lp, ST_POINTS, tc)
      end if

      ! GRID RESPONSE, point-motion piece, scattered: every point's
      ! `qmove` lands on the atom that OWNS the point (xg%atom), not on
      ! the atom of any basis function -- see this module's header.
      call tbl%step_on(lp, ST_SCATTER, tc)
      if (gresp) call xc_grad_scatter_points(xg%npts, xg%atom, xg%b_off(ib), npb, qmove, b%natm, grad)
      call tbl%step_off(lp, ST_SCATTER, tc)

      call tbl%step_on(lp, ST_GEMM, tc)
      do s = 1, nspin
         call lp%gemm_mat(pall(:, :, s), chi, tall(:, :, s), transa="T", beta=0.0_dp)
      end do
      if (is_gga) then
         do s = 1, nspin
            call lp%gemm_mat(q1all(:, :, s), rall(:, :, s), tall(:, :, s), transa="T", beta=1.0_dp)
         end do
         if (nspin == 2) then
            call lp%gemm_mat(q2all(:, :, 1), rall(:, :, 2), tall(:, :, 1), transa="T", beta=1.0_dp)
            call lp%gemm_mat(q2all(:, :, 2), rall(:, :, 1), tall(:, :, 2), transa="T", beta=1.0_dp)
         end if
      end if
      call tbl%step_off(lp, ST_GEMM, tc)

      call tbl%step_on(lp, ST_SCATTER, tc)
      call xc_grad_scatter(ib, xg%nbatch, xg%b_aooff, size(xg%b_ao), xg%b_ao, ao_atom, b%nao, &
                           nloc, nspin, dloc, tall, b%natm, grad)
      call tbl%step_off(lp, ST_SCATTER, tc)

      ! GRID RESPONSE, weight-derivative piece -- HOST side (see this
      ! module's header and xc_gradient_add's own comment on why): needs
      ! only `fxc` (rho*eps or (rho_a+rho_b)*eps per point, no weight
      ! baked in) brought back from the device, and `xg%r`/`xg%atom`,
      ! which are valid on the host regardless of their device mirror
      ! (xg%build never modifies them again after this batch's collocation
      ! reads the device copy).
      if (gresp) then
         call tbl%step_on(lp, ST_PDERIV, tc)
         !$acc update self(fxc)
         g0 = xg%b_off(ib)
         allocate (dweights(3, b%natm, npb))
         call becke_partition_derivatives(xg%r(:, g0:g0 + npb - 1), at_r, at_z, &
                                          xg%atom(g0:g0 + npb - 1), scheme, adjust, dweights, err)
         if (err%has_error()) error stop "trc_xc_grad: "//trim(err%get_message())
         call tbl%step_off(lp, ST_PDERIV, tc)
         call tbl%step_on(lp, ST_RESP, tc)
         do g = 1, npb
            contrib = xg%qw(g0 + g - 1)*fxc(g)
            own = xg%atom(g0 + g - 1)
            do a = 1, b%natm
               grad_extra(:, a) = grad_extra(:, a) + contrib*dweights(:, a, g)
            end do
            ! The point travels with `own`, and the partition is a function
            ! of relative positions only, so its derivative with respect to
            ! the POINT is minus the sum of its derivatives with respect to
            ! every nucleus (translational invariance of a function of
            ! r - R_1, ..., r - R_natm) -- see this module's header.
            do a = 1, b%natm
               grad_extra(:, own) = grad_extra(:, own) - contrib*dweights(:, a, g)
            end do
         end do
         deallocate (dweights)
         call tbl%step_off(lp, ST_RESP, tc)
      end if

      call tbl%step_on(lp, ST_ALLOC, tc)
      !$acc exit data delete(chi, gchi, dloc, xmat, pall, tall, fxc, qmove)
      deallocate (chi, gchi, dloc, xmat, pall, tall, fxc, qmove)
      if (is_gga) then
         !$acc exit data delete(hchi, q1all, q2all, rall, rmat, rmatx)
         deallocate (hchi, q1all, q2all, rall, rmat, rmatx)
      end if
      call tbl%step_off(lp, ST_ALLOC, tc)
   end subroutine xc_grad_one_batch

   !
   ! Collocation, order 1: one thread per (local shell, local point) of
   ! this one batch. Same mapping trc_xc.F90's xc_collocate_body uses
   ! (`b_shao` is the shell's first LOCAL function, 0-based), just over one
   ! batch's own (ks, gl) pair instead of a flattened multi-batch work list.
   !
   subroutine xc_grad_collocate1(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                 nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, chi, gchi)
      integer, intent(in) :: nshell, maxnp, npts, nbatch, nbsh, ib, npb, nloc
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      real(dp), intent(inout) :: chi(npb, nloc), gchi(npb, nloc, 3)
      integer :: nsh_b, t

      ! A genuine two-index `do concurrent(ks=1:nsh_b, gl=1:npb)` is the
      ! bug this module shipped with: on this GPU build, nvfortran's
      ! codegen for a 2-D `do concurrent` with a small, batch-varying outer
      ! extent (nsh_b, here as small as 1-4) against a much larger inner
      ! extent (npb, up to 512) silently drops the LAST outer slice --
      ! confirmed by an element-by-element CPU/GPU dump of `chi`/`gchi` for
      ! one batch (nsh_b=3, npb=86): columns for ks=1,2 matched CPU to
      ! 1e-9, the ENTIRE ks=3 column (shell index 4, local AO 3) came back
      ! bitwise 0 on GPU, (nsh_b-1)*npb = 172 of 258 elements right, the
      ! rest exactly the zero a freshly `create`-mapped device buffer
      ! starts as. Never seen on CPU, where `do concurrent` is just nested
      ! host loops. The energy path (trc_xc.F90's xc_collocate) never
      ! triggers this: it flattens every batch's work into ONE `do
      ! concurrent(t=1:nwork)` before it ever reaches the device, and
      ! recovers (ks, gl) INSIDE the pure body (xc_collocate_body) -- the
      ! exact convention adopted here, not just the flattening.
      nsh_b = b_shoff(ib + 1) - b_shoff(ib)
      do concurrent(t=1:nsh_b*npb)
         call xc_grad_colloc1_body(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                   nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                   t, chi, gchi)
      end do
   end subroutine xc_grad_collocate1

   pure subroutine xc_grad_colloc1_body(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                        nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                        t, chi, gchi)
      !$acc routine seq
      integer, intent(in) :: nshell, maxnp, npts, nbatch, nbsh, ib, npb, nloc, t
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      real(dp), intent(inout) :: chi(npb, nloc), gchi(npb, nloc, 3)
      integer :: ks, gl, ish, iao, l, nc, m, gg
      real(dp) :: d(3), cs(NCART_MAX), gs(3, NCART_MAX)

      ks = (t - 1)/npb + 1
      gl = t - (ks - 1)*npb
      ish = b_sh(b_shoff(ib) + ks - 1)
      iao = b_shao(b_shoff(ib) + ks - 1)
      gg = b_off(ib) + gl - 1
      l = sh_l(ish)
      d(1) = r(1, gg) - sh_r(1, ish)
      d(2) = r(2, gg) - sh_r(2, ish)
      d(3) = r(3, gg) - sh_r(3, ish)
      call shell_collocate(l, sh_np(ish), sh_e(1:sh_np(ish), ish), sh_c(1:sh_np(ish), ish), d, cs, gs)
      nc = (l + 1)*(l + 2)/2
      do m = 1, nc
         chi(gl, iao + m) = cs(m)
         gchi(gl, iao + m, 1) = gs(1, m)
         gchi(gl, iao + m, 2) = gs(2, m)
         gchi(gl, iao + m, 3) = gs(3, m)
      end do
   end subroutine xc_grad_colloc1_body

   !
   ! As `xc_grad_collocate1`, order 1 and 2 (the GGA Hessian term).
   !
   subroutine xc_grad_collocate2(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                 nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                 chi, gchi, hchi)
      integer, intent(in) :: nshell, maxnp, npts, nbatch, nbsh, ib, npb, nloc
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      real(dp), intent(inout) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      integer :: nsh_b, t

      ! Flattened to one index, (ks, gl) recovered inside the pure body,
      ! for the same reason xc_grad_collocate1 is -- see that subroutine's
      ! comment.
      nsh_b = b_shoff(ib + 1) - b_shoff(ib)
      do concurrent(t=1:nsh_b*npb)
         call xc_grad_colloc2_body(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                   nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                   t, chi, gchi, hchi)
      end do
   end subroutine xc_grad_collocate2

   pure subroutine xc_grad_colloc2_body(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                        nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                        t, chi, gchi, hchi)
      !$acc routine seq
      integer, intent(in) :: nshell, maxnp, npts, nbatch, nbsh, ib, npb, nloc, t
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      real(dp), intent(inout) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      integer :: ks, gl, ish, iao, l, nc, m, gg
      real(dp) :: d(3), cs(NCART_MAX), gs(3, NCART_MAX), hs(NHESS, NCART_MAX)

      ks = (t - 1)/npb + 1
      gl = t - (ks - 1)*npb
      ish = b_sh(b_shoff(ib) + ks - 1)
      iao = b_shao(b_shoff(ib) + ks - 1)
      gg = b_off(ib) + gl - 1
      l = sh_l(ish)
      d(1) = r(1, gg) - sh_r(1, ish)
      d(2) = r(2, gg) - sh_r(2, ish)
      d(3) = r(3, gg) - sh_r(3, ish)
      call shell_collocate2(l, sh_np(ish), sh_e(1:sh_np(ish), ish), sh_c(1:sh_np(ish), ish), d, cs, gs, hs)
      nc = (l + 1)*(l + 2)/2
      do m = 1, nc
         chi(gl, iao + m) = cs(m)
         gchi(gl, iao + m, 1) = gs(1, m)
         gchi(gl, iao + m, 2) = gs(2, m)
         gchi(gl, iao + m, 3) = gs(3, m)
         hchi(gl, iao + m, 1) = hs(1, m)
         hchi(gl, iao + m, 2) = hs(2, m)
         hchi(gl, iao + m, 3) = hs(3, m)
         hchi(gl, iao + m, 4) = hs(4, m)
         hchi(gl, iao + m, 5) = hs(5, m)
         hchi(gl, iao + m, 6) = hs(6, m)
      end do
   end subroutine xc_grad_colloc2_body

   !
   ! D sliced onto this batch's local basis, every (mu, nu, spin) of it --
   ! no padding, no screening within the batch, unlike xc_gather_body: this
   ! module's tiles are exact-sized, so there is no padded element to zero.
   !
   subroutine xc_grad_gather(nao, nbao, b_ao, aooff0, nloc, nspin, dens, dloc)
      integer, intent(in) :: nao, nbao, aooff0, nloc, nspin
      integer, intent(in) :: b_ao(nbao)
      real(dp), intent(in) :: dens(nao, nao, nspin)
      real(dp), intent(inout) :: dloc(nloc, nloc, nspin)
      integer :: mu, nu, s
      do concurrent(mu=1:nloc, nu=1:nloc, s=1:nspin)
         dloc(mu, nu, s) = dens(b_ao(aooff0 + mu - 1), b_ao(aooff0 + nu - 1), s)
      end do
   end subroutine xc_grad_gather

   !
   ! LDA: one thread per point. Builds rho, grad rho (needed even at LDA,
   ! since the moving function's OWN gradient enters the Pulay term), the
   ! potential, and the "P" tile `pall(g, (k-1)*nloc+mu, s)` the force
   ! GEMM `pall^T @ chi` contracts into T_k(mu,nu). Restricted and
   ! unrestricted in the same body, as trc_xc.F90's own xc_point_body does.
   !
   subroutine xc_grad_point_lda(npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, xmat, pall, &
                                fxc, qmove)
      integer, intent(in) :: npb, npts, p0, nloc, nspin, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), xmat(npb, nloc, nspin)
      real(dp), intent(inout) :: pall(npb, 3*nloc, nspin)
      real(dp), intent(inout) :: fxc(npb), qmove(npb, 3)
      integer :: g
      do concurrent(g=1:npb)
         call xc_grad_point_lda_body(g, npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, xmat, &
                                     pall, fxc, qmove)
      end do
   end subroutine xc_grad_point_lda

   pure subroutine xc_grad_point_lda_body(g, npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, &
                                          xmat, pall, fxc, qmove)
      !$acc routine seq
      integer, intent(in) :: g, npb, npts, p0, nloc, nspin, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), xmat(npb, nloc, nspin)
      real(dp), intent(inout) :: pall(npb, 3*nloc, nspin)
      real(dp), intent(inout) :: fxc(npb), qmove(npb, 3)
      real(dp) :: rho(2), grho(3, 2)
      real(dp) :: sigma, eps, vrho, vsigma, wg
      real(dp) :: va, vb, s_aa, s_ab, s_bb, v_aa, v_ab, v_bb
      real(dp) :: qm(3)
      integer :: s, u, mu, k

      wg = w(p0 + g - 1)
      rho = 0.0_dp
      grho = 0.0_dp
      do s = 1, nspin
         do u = 1, nloc
            rho(s) = rho(s) + xmat(g, u, s)*chi(g, u)
            grho(1, s) = grho(1, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 1)
            grho(2, s) = grho(2, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 2)
            grho(3, s) = grho(3, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 3)
         end do
      end do

      if (nspin == 1) then
         sigma = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         call xc_eval_point(nk, kid, coef, rho(1), sigma, eps, vrho, vsigma)
         do mu = 1, nloc
            do k = 1, 3
               pall(g, (k - 1)*nloc + mu, 1) = wg*vrho*gchi(g, mu, k)
            end do
         end do
         fxc(g) = rho(1)*eps
      else
         s_aa = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         s_bb = grho(1, 2)**2 + grho(2, 2)**2 + grho(3, 2)**2
         s_ab = grho(1, 1)*grho(1, 2) + grho(2, 1)*grho(2, 2) + grho(3, 1)*grho(3, 2)
         call xc_eval_point_polar(nk, kid, coef, rho(1), rho(2), s_aa, s_ab, s_bb, &
                                  eps, va, vb, v_aa, v_ab, v_bb)
         do mu = 1, nloc
            do k = 1, 3
               pall(g, (k - 1)*nloc + mu, 1) = wg*va*gchi(g, mu, k)
               pall(g, (k - 1)*nloc + mu, 2) = wg*vb*gchi(g, mu, k)
            end do
         end do
         fxc(g) = (rho(1) + rho(2))*eps
      end if

      ! GRID RESPONSE, point-motion piece, LDA -- see this module's
      ! header. The same contraction the Pulay force GEMM performs
      ! (pall^T @ chi, then dotted with D), but summed over EVERY local
      ! function rather than only the ones on one atom, with the
      ! OPPOSITE sign, and kept per-point rather than pre-summed over the
      ! batch -- because the atom this lands on is the POINT's owner,
      ! which can differ from point to point within one spatial batch.
      ! `xmat(g, mu, s)` already equals `sum_nu D(mu,nu,s) chi(nu,g)` (D
      ! symmetric), so this is the full `D : pall` contraction for this
      ! point alone, times the "+2" the module header derives.
      qm = 0.0_dp
      do mu = 1, nloc
         do k = 1, 3
            qm(k) = qm(k) + pall(g, (k - 1)*nloc + mu, 1)*xmat(g, mu, 1)
         end do
      end do
      if (nspin == 2) then
         do mu = 1, nloc
            do k = 1, 3
               qm(k) = qm(k) + pall(g, (k - 1)*nloc + mu, 2)*xmat(g, mu, 2)
            end do
         end do
      end if
      qmove(g, 1) = 2.0_dp*qm(1)
      qmove(g, 2) = 2.0_dp*qm(2)
      qmove(g, 3) = 2.0_dp*qm(3)
   end subroutine xc_grad_point_lda_body

   !
   ! GGA: as the LDA point kernel, plus the Hessian-dependent "P" piece and
   ! the "Q"/"R" pair the second force GEMM (`q^T @ r`) contracts. R_nu(g)
   ! = grad_rho(g) . grad_chi_nu(g) is the ONE quantity shared by every mu
   ! on this point -- built once per nu, not recomputed per (mu, nu) pair.
   !
   ! Unrestricted carries two cross terms the restricted branch does not:
   ! channel a's P and Q1/Q2 read v_ab against channel b's Hessian/R and
   ! vice versa (Ca_k = Pa_k phi_nu + Qa1_k Ra_nu + Qa2_k Rb_nu, this
   ! module's header), so `q2all(:,:,1)` pairs with `rall(:,:,2)` and
   ! `q2all(:,:,2)` with `rall(:,:,1)` in the caller's GEMM step.
   !
   subroutine xc_grad_point_gga(npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, hchi, xmat, &
                                pall, q1all, q2all, rall, fxc, qmove)
      integer, intent(in) :: npb, npts, p0, nloc, nspin, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      real(dp), intent(in) :: xmat(npb, nloc, nspin)
      real(dp), intent(inout) :: pall(npb, 3*nloc, nspin), q1all(npb, 3*nloc, nspin)
      real(dp), intent(inout) :: q2all(npb, 3*nloc, nspin), rall(npb, nloc, nspin)
      real(dp), intent(inout) :: fxc(npb), qmove(npb, 3)
      integer :: g
      do concurrent(g=1:npb)
         call xc_grad_point_gga_body(g, npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, hchi, &
                                     xmat, pall, q1all, q2all, rall, fxc, qmove)
      end do
   end subroutine xc_grad_point_gga

   pure subroutine xc_grad_point_gga_body(g, npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, &
                                          hchi, xmat, pall, q1all, q2all, rall, fxc, qmove)
      !$acc routine seq
      integer, intent(in) :: g, npb, npts, p0, nloc, nspin, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      real(dp), intent(in) :: xmat(npb, nloc, nspin)
      real(dp), intent(inout) :: pall(npb, 3*nloc, nspin), q1all(npb, 3*nloc, nspin)
      real(dp), intent(inout) :: q2all(npb, 3*nloc, nspin), rall(npb, nloc, nspin)
      real(dp), intent(inout) :: fxc(npb), qmove(npb, 3)
      real(dp) :: qm(3)
      real(dp) :: rho(2), grho(3, 2)
      real(dp) :: sigma, eps, vrho, vsigma, wg
      real(dp) :: va, vb, s_aa, s_ab, s_bb, v_aa, v_ab, v_bb
      real(dp) :: hdota(3), hdotb(3)
      integer :: s, u, mu, k

      wg = w(p0 + g - 1)
      rho = 0.0_dp
      grho = 0.0_dp
      do s = 1, nspin
         do u = 1, nloc
            rho(s) = rho(s) + xmat(g, u, s)*chi(g, u)
            grho(1, s) = grho(1, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 1)
            grho(2, s) = grho(2, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 2)
            grho(3, s) = grho(3, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 3)
         end do
      end do

      if (nspin == 1) then
         sigma = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         call xc_eval_point(nk, kid, coef, rho(1), sigma, eps, vrho, vsigma)
         do mu = 1, nloc
            ! (grad_rho . H_mu)_k, hchi packed xx, xy, xz, yy, yz, zz.
            hdota(1) = grho(1, 1)*hchi(g, mu, 1) + grho(2, 1)*hchi(g, mu, 2) + grho(3, 1)*hchi(g, mu, 3)
            hdota(2) = grho(1, 1)*hchi(g, mu, 2) + grho(2, 1)*hchi(g, mu, 4) + grho(3, 1)*hchi(g, mu, 5)
            hdota(3) = grho(1, 1)*hchi(g, mu, 3) + grho(2, 1)*hchi(g, mu, 5) + grho(3, 1)*hchi(g, mu, 6)
            do k = 1, 3
               pall(g, (k - 1)*nloc + mu, 1) = wg*(vrho*gchi(g, mu, k) + 2.0_dp*vsigma*hdota(k))
               q1all(g, (k - 1)*nloc + mu, 1) = wg*2.0_dp*vsigma*gchi(g, mu, k)
            end do
         end do
         do u = 1, nloc
            rall(g, u, 1) = grho(1, 1)*gchi(g, u, 1) + grho(2, 1)*gchi(g, u, 2) + grho(3, 1)*gchi(g, u, 3)
         end do
         fxc(g) = rho(1)*eps
      else
         s_aa = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         s_bb = grho(1, 2)**2 + grho(2, 2)**2 + grho(3, 2)**2
         s_ab = grho(1, 1)*grho(1, 2) + grho(2, 1)*grho(2, 2) + grho(3, 1)*grho(3, 2)
         call xc_eval_point_polar(nk, kid, coef, rho(1), rho(2), s_aa, s_ab, s_bb, &
                                  eps, va, vb, v_aa, v_ab, v_bb)
         do mu = 1, nloc
            hdota(1) = grho(1, 1)*hchi(g, mu, 1) + grho(2, 1)*hchi(g, mu, 2) + grho(3, 1)*hchi(g, mu, 3)
            hdota(2) = grho(1, 1)*hchi(g, mu, 2) + grho(2, 1)*hchi(g, mu, 4) + grho(3, 1)*hchi(g, mu, 5)
            hdota(3) = grho(1, 1)*hchi(g, mu, 3) + grho(2, 1)*hchi(g, mu, 5) + grho(3, 1)*hchi(g, mu, 6)
            hdotb(1) = grho(1, 2)*hchi(g, mu, 1) + grho(2, 2)*hchi(g, mu, 2) + grho(3, 2)*hchi(g, mu, 3)
            hdotb(2) = grho(1, 2)*hchi(g, mu, 2) + grho(2, 2)*hchi(g, mu, 4) + grho(3, 2)*hchi(g, mu, 5)
            hdotb(3) = grho(1, 2)*hchi(g, mu, 3) + grho(2, 2)*hchi(g, mu, 5) + grho(3, 2)*hchi(g, mu, 6)
            do k = 1, 3
               pall(g, (k - 1)*nloc + mu, 1) = wg*(va*gchi(g, mu, k) + 2.0_dp*v_aa*hdota(k) + v_ab*hdotb(k))
               pall(g, (k - 1)*nloc + mu, 2) = wg*(vb*gchi(g, mu, k) + 2.0_dp*v_bb*hdotb(k) + v_ab*hdota(k))
               q1all(g, (k - 1)*nloc + mu, 1) = wg*2.0_dp*v_aa*gchi(g, mu, k)   ! x rall(:,:,1)
               q2all(g, (k - 1)*nloc + mu, 1) = wg*v_ab*gchi(g, mu, k)         ! x rall(:,:,2)
               q1all(g, (k - 1)*nloc + mu, 2) = wg*2.0_dp*v_bb*gchi(g, mu, k)   ! x rall(:,:,2)
               q2all(g, (k - 1)*nloc + mu, 2) = wg*v_ab*gchi(g, mu, k)         ! x rall(:,:,1)
            end do
         end do
         do u = 1, nloc
            rall(g, u, 1) = grho(1, 1)*gchi(g, u, 1) + grho(2, 1)*gchi(g, u, 2) + grho(3, 1)*gchi(g, u, 3)
            rall(g, u, 2) = grho(1, 2)*gchi(g, u, 1) + grho(2, 2)*gchi(g, u, 2) + grho(3, 2)*gchi(g, u, 3)
         end do
         fxc(g) = (rho(1) + rho(2))*eps
      end if

      ! GRID RESPONSE, point-motion piece, the LDA-like part only -- see
      ! this module's header and xc_grad_point_lda_body's own comment on
      ! the sign and the "+2". The GGA-specific addition (the `q1`/`q2`
      ! pieces, which need `rmat`/`rmatx` -- a contraction across every
      ! point of the batch at once, not available inside this per-point
      ! body) is ADDED on top by xc_grad_pointmotion_gga, called once this
      ! kernel has run over the whole batch.
      qm = 0.0_dp
      do mu = 1, nloc
         do k = 1, 3
            qm(k) = qm(k) + pall(g, (k - 1)*nloc + mu, 1)*xmat(g, mu, 1)
         end do
      end do
      if (nspin == 2) then
         do mu = 1, nloc
            do k = 1, 3
               qm(k) = qm(k) + pall(g, (k - 1)*nloc + mu, 2)*xmat(g, mu, 2)
            end do
         end do
      end if
      qmove(g, 1) = 2.0_dp*qm(1)
      qmove(g, 2) = 2.0_dp*qm(2)
      qmove(g, 3) = 2.0_dp*qm(3)
   end subroutine xc_grad_point_gga_body

   !
   ! GRID RESPONSE, point-motion piece: the GGA-specific addition to
   ! `qmove` that xc_grad_point_gga_body could not compute itself -- see
   ! this module's header and that body's own comment. `rmat(g, mu, s) =
   ! sum_nu D(mu,nu,s) rall(nu,g,s)` and `rmatx` its cross-spin twin
   ! (xc_grad_one_batch's own comment on the GEMMs that build them) are
   ! the exact counterparts of `xmat` for the "R" tile the way `xmat`
   ! itself is for `chi`, which is what lets this follow the SAME "+2 x
   ! the full-basis contraction, no atom restriction" pattern as the
   ! LDA-like piece.
   !
   subroutine xc_grad_pointmotion_gga(npb, nloc, nspin, q1all, q2all, rmat, rmatx, qmove)
      integer, intent(in) :: npb, nloc, nspin
      real(dp), intent(in) :: q1all(npb, 3*nloc, nspin), q2all(npb, 3*nloc, nspin)
      real(dp), intent(in) :: rmat(npb, nloc, nspin), rmatx(npb, nloc, nspin)
      real(dp), intent(inout) :: qmove(npb, 3)
      integer :: g
      do concurrent(g=1:npb)
         call xc_grad_pointmotion_gga_body(g, npb, nloc, nspin, q1all, q2all, rmat, rmatx, qmove)
      end do
   end subroutine xc_grad_pointmotion_gga

   pure subroutine xc_grad_pointmotion_gga_body(g, npb, nloc, nspin, q1all, q2all, rmat, rmatx, qmove)
      !$acc routine seq
      integer, intent(in) :: g, npb, nloc, nspin
      real(dp), intent(in) :: q1all(npb, 3*nloc, nspin), q2all(npb, 3*nloc, nspin)
      real(dp), intent(in) :: rmat(npb, nloc, nspin), rmatx(npb, nloc, nspin)
      real(dp), intent(inout) :: qmove(npb, 3)
      integer :: mu, k
      real(dp) :: acc(3)

      acc = 0.0_dp
      if (nspin == 1) then
         do mu = 1, nloc
            do k = 1, 3
               acc(k) = acc(k) + q1all(g, (k - 1)*nloc + mu, 1)*rmat(g, mu, 1)
            end do
         end do
      else
         do mu = 1, nloc
            do k = 1, 3
               acc(k) = acc(k) + q1all(g, (k - 1)*nloc + mu, 1)*rmat(g, mu, 1) &
                        + q1all(g, (k - 1)*nloc + mu, 2)*rmat(g, mu, 2) &
                        + q2all(g, (k - 1)*nloc + mu, 1)*rmatx(g, mu, 1) &
                        + q2all(g, (k - 1)*nloc + mu, 2)*rmatx(g, mu, 2)
            end do
         end do
      end if
      qmove(g, 1) = qmove(g, 1) + 2.0_dp*acc(1)
      qmove(g, 2) = qmove(g, 2) + 2.0_dp*acc(2)
      qmove(g, 3) = qmove(g, 3) + 2.0_dp*acc(3)
   end subroutine xc_grad_pointmotion_gga_body

   !
   ! GRID RESPONSE, point-motion piece, scattered: `qmove(g, :)` lands on
   ! the atom OWNING point `g` (`xg%atom`, passed here as `g_atom` plus
   ! the batch's global offset `g0` -- the whole array and an offset, not
   ! a pre-sliced section, the same convention `w`/`p0` already use for
   ! the point weight, because a device-resident derived-type component
   ! is passed around this module as the FULL array). Atomic for the same
   ! reason xc_grad_scatter is: two batches, or two points of this very
   ! loop, can own the same atom.
   !
   subroutine xc_grad_scatter_points(npts, g_atom, g0, npb, qmove, natm, grad)
      integer, intent(in) :: npts, g0, npb, natm
      integer, intent(in) :: g_atom(npts)
      real(dp), intent(in) :: qmove(npb, 3)
      real(dp), intent(inout) :: grad(3, natm)
      integer :: g, k, t, at
#ifdef _OPENACC
      ! Flattened to one index, (g, k) recovered inside the loop body, for
      ! the same reason xc_grad_collocate1/2 are (this module's header, and
      ! those subroutines' own comments): a genuine two-index
      ! `do concurrent(g=1:npb, k=1:3)` is the exact construct family that
      ! silently dropped a slice on this GPU build in a2e6768, and k=3 here
      ! is every bit as small and fixed as the outer extent that triggered
      ! it there.
      do concurrent(t=1:3*npb) local(g, k, at)
         g = (t - 1)/3 + 1
         k = mod(t - 1, 3) + 1
         at = g_atom(g0 + g - 1)
         !$acc atomic update
         grad(k, at) = grad(k, at) + qmove(g, k)
      end do
#else
      do g = 1, npb
         at = g_atom(g0 + g - 1)
         do k = 1, 3
            grad(k, at) = grad(k, at) + qmove(g, k)
         end do
      end do
#endif
   end subroutine xc_grad_scatter_points

   !
   ! T_k(mu, nu) = tall((k-1)*nloc+mu, nu, s), contracted with D and summed
   ! over spin, lands as -2 * (that sum) on the atom owning local function
   ! mu -- one thread per (mu, k) of this batch, atomic because two batches
   ! (or two points of the atom-scatter loop below) can own the same atom.
   ! Same #ifdef pattern as trc_xc.F90's xc_scatter, for the same reason:
   ! a host compiler that takes `do concurrent` at its word reorders the
   ! read-modify-writes where the device's atomic makes them correct.
   !
   subroutine xc_grad_scatter(ib, nbatch, b_aooff, nbao, b_ao, ao_atom, nao, nloc, nspin, &
                              dloc, tall, natm, grad)
      integer, intent(in) :: ib, nbatch, nbao, nao, nloc, nspin, natm
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), ao_atom(nao)
      real(dp), intent(in) :: dloc(nloc, nloc, nspin), tall(3*nloc, nloc, nspin)
      real(dp), intent(inout) :: grad(3, natm)
      integer :: mu, k, at
      real(dp) :: fval
#ifdef _OPENACC
      do concurrent(mu=1:nloc, k=1:3) local(fval, at)
         call xc_grad_scatter_body(mu, k, nloc, nspin, dloc, tall, fval)
         at = ao_atom(b_ao(b_aooff(ib) + mu - 1))
         !$acc atomic update
         grad(k, at) = grad(k, at) + fval
      end do
#else
      do mu = 1, nloc
         do k = 1, 3
            call xc_grad_scatter_body(mu, k, nloc, nspin, dloc, tall, fval)
            at = ao_atom(b_ao(b_aooff(ib) + mu - 1))
            grad(k, at) = grad(k, at) + fval
         end do
      end do
#endif
   end subroutine xc_grad_scatter

   pure subroutine xc_grad_scatter_body(mu, k, nloc, nspin, dloc, tall, fval)
      !$acc routine seq
      integer, intent(in) :: mu, k, nloc, nspin
      real(dp), intent(in) :: dloc(nloc, nloc, nspin), tall(3*nloc, nloc, nspin)
      real(dp), intent(out) :: fval
      integer :: nu, s
      real(dp) :: acc
      acc = 0.0_dp
      do s = 1, nspin
         do nu = 1, nloc
            acc = acc + dloc(mu, nu, s)*tall((k - 1)*nloc + mu, nu, s)
         end do
      end do
      fval = -2.0_dp*acc
   end subroutine xc_grad_scatter_body

end module trc_xc_grad
