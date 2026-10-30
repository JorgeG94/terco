!
! VV10 (Vydrov & Van Voorhis 2010) nonlocal correlation's contribution to the
! nuclear gradient, RKS and UKS, with grid response. Added into an existing
! `grad(3,natm)`, the same convention trc_xc_grad%xc_grad_add's semilocal
! term uses.
!
! THE FORMULA. The non-locality is already inside vrho/vsigma: E_nl is a
! functional of rho and sigma alone,
!
!     E_nl = int dr rho(r) [ beta + 1/2 int dr' rho(r') Phi(r,r') ]
!     dE_nl/dR = int dr [ vrho(r) drho(r)/dR + vsigma(r) dsigma(r)/dR ]
!
! and the O(N_grid^2) pair sum (trc_vv10's PASS 2) is spent entirely on
! producing the two arrays `vrho`/`vsigma` at every point, exactly the way
! trc_vv10's own energy path does. What is left is the ORDINARY GGA
! contraction trc_xc_grad.F90 already implements for a semilocal kernel's
! vrho/vsigma -- this module reuses that pipeline (trc_xc_grad's
! xc_grad_drive, restricted GGA, with this module's own batch screen),
! just reading `vrho(g)`/`vsigma(g)`/grad rho from VV10's own precomputed
! arrays instead of calling `xc_eval_point`.
!
! THE NONLOCAL EXTRAS. Two things have no counterpart in a semilocal
! functional's gradient, and neither is optional (measured at ~4.5e-4 on
! water, and it does not shrink with the grid):
!
!   - the PARTITION-WEIGHT term's coefficient is `dE/dw_i`, not `rho_i eps_i`.
!     A point's weight multiplies the pair sum TWICE over (once as the outer
!     point, once inside every other point's inner sum), and the kernel's
!     symmetry in its two position arguments makes those equal, so `dE/dw`
!     carries `beta + F` where the energy density `eps` carries `beta + F/2`
!     (trc_vv10.F90's own header: `eps = beta + 0.5 F`, `vrho` already needs
!     `beta + F` on its own account -- `dedw` reuses that same combination,
!     this module's own `vv10_grad_kernel_body`);
!
!   - the points carry an EXPLICIT force of their own, because the kernel is
!     a function of where they ARE, not only of what the density is there --
!     nothing in a semilocal functional's `eps(rho, sigma)` depends on `r`
!     directly the way VV10's `r^2` denominator does. One pass over outer
!     point `i` alone is enough to capture both roles symmetry gives it (as
!     an outer point in its own sum, and as an inner point inside every
!     other point's sum): `vv10_nlc`'s `fexp` is `-3 r_i` times the same
!     per-pair sum `vrho`'s own `dk/dr`/`dw0/dr` terms are built from, not
!     `-1.5 r_i` -- the doubled coefficient IS that symmetry, not a second
!     loop. Landed on the point's OWNER at the point's FULL (Becke-weighted)
!     weight `xg%w`, not the bare quadrature weight `xg%qw` the partition
!     term below uses.
!
! Both grid-response pieces close the SAME header gap xc_grad_add's own
! "GRID RESPONSE" section documents for the semilocal term (point motion of
! both grids, partition-weight motion, and their cross-cancellation via
! `dft_partition_becke_grad_contract`) -- this module's own weight/point-motion
! block is that section's VV10 analogue, with `dedw`/`fexp` in place of
! `fxc = rho eps`/nothing.
!
! ENERGY-SIDE REUSE. `vv10_density_tile` (PASS 1: rho/grad-rho/sigma at every
! outer point) and `vv10_gather` (the inner-grid survivor list `vv10_kernel`
! reads) are trc_vv10's own, exported unchanged -- the gradient's outer-point
! sweep is bit-identical O(N) work to the energy path's, so it is not
! re-derived here. PASS 2 itself (`vv10_kernel_body`) is NOT reused as-is:
! it has no `dedw`/`fexp` outputs, and `vv10_grad_kernel_body` below is a
! second, gradient-only pass-2 body instead, reading the exact same
! gathered arrays and the exact same `kvv`/`beta`/`w0`/`k` formulas (see
! that routine's own comment for the one place its sign convention
! differs). Folding the two bodies into one, with `dedw`/`fexp` as
! optional outputs the energy path's own call site simply omits, is a
! reasonable follow-up refactor, not attempted here.
!
! RESTRICTED AND UNRESTRICTED. VV10 is a functional of the TOTAL density
! alone (trc_vv10.F90's own header), so this module's single entry point
! takes `dens(nao,nao,nspin)` exactly as xc_grad_add does, forms
! `dtot = dens(:,:,1) [+ dens(:,:,2)]` once, and runs the entire gradient
! against that one density -- there is no separate alpha/beta contraction to
! carry, unlike the semilocal GGA gradient's cross-spin Q2/rmatx terms.
!
! GRID. This module builds its OWN quadrature at `b`'s geometry, the same
! reason xc_grad_add does (trc_scf_context keeps neither grid nor
! functional alive past convergence) -- at the CALLER's chosen level
! (`grid_level`, expected to be the SCF's own `nlc_grid_level`, not the
! semilocal grid's `grid_level`: VV10 runs on its own, coarser grid
! end-to-end, energy and gradient alike), with `atom`/`qw` attached so the
! weight-response block below can read them -- the same two optional
! arguments xc_grad_add passes its own grid build.
!
! SINGLE RANK. No communicator: PASS 2 here is a second, small O(N_grid^2)
! sum next to the energy path's (coarser grid, no potential-matrix pass 3 to
! pay for), and splitting it the way trc_vv10's own PASS 2 does would need an
! allreduce over `vrho`/`vsigma`/`dedw`/`fexp` before the batch loop that
! consumes them could run device-resident on every rank -- a real plumbing
! change (a `comm_t` threaded through this module and trc_scf_gradient,
! which takes none today) rather than something that falls out of the
! existing single-rank call sites. `trc_run_gradient` already refuses a
! multi-rank context outright (trc_capi.F90, "the gradient path is not
! rank-split"), so every caller of this module is single-rank already;
! left that way rather than half-plumbed.
!
! GPU. Measured correct and consistent with CPU to ~1e-11/1e-12 (the
! fixed-density FD check, check_vv10_grad_fixedd_fd.F90, three repeated
! GPU runs) -- built from the same idioms (`do concurrent`/`routine seq`,
! the `TRC_DEVICE` CPU/GPU pass-2 split trc_vv10.F90's own header explains)
! as the two modules it reuses. Not yet measured for its OWN cost next to
! the rest of the gradient at more than one case; see the branch report
! that added this module for the one case it does measure.
!
module trc_vv10_grad
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_constants, only: PI
   use trc_basis_types, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_grad, only: xc_grad_drive, xc_grad_step_names, ST_GRID, ST_POINTS, ST_PDERIV, ST_RESP
   use trc_xc_functional, only: XC_MAX_KERNELS
   use trc_step_table, only: step_table_t
   use trc_vv10, only: vv10_density_tile, vv10_gather, VV10_RHO_TOL, PI43
   use trc_linalg, only: trc_linalg_t
   use trc_dft_grid, only: dft_grid_t, dft_grid_build, DEFAULT_GRID_LEVEL
   use trc_dft_partition, only: dft_partition_becke_grad_contract
   use trc_error, only: error_t
   implicit none
   private

   public :: vv10_grad_add

   integer(int64), parameter :: VV10_GRAD_TILE_CAP = 8388608_int64
   !! The chunk-size budget (padded rows x padded columns) handed to
   !! `vv10_density_tile`'s own batching loop -- trc_vv10.F90's `vv10_drive`
   !! uses the identical literal for its own PASS 1 call, unnamed; named
   !! here rather than duplicated unnamed a second time.
   integer, parameter :: VV10_GRAD_DEFAULT_MAX_PTS = 512
   !! Default points-per-batch cap for this module's OWN grid build
   !! (`xg%build`'s `max_pts`), when the caller's `grid_max_pts` is
   !! absent -- the same default xc_grad_add and trc_scf_context use
   !! for the semilocal grid.

contains

   !
   ! The one entry point. `dens(b%nao,b%nao,nspin)` -- the converged density,
   ! per spin channel, same convention xc_grad_add uses. `bcoef`/`ccoef`
   ! are VV10's own `b`/`C` (6.0, 0.01 for wB97X-V). ADDS into `grad`.
   !
   subroutine vv10_grad_add(b, lp, nspin, density, bcoef, ccoef, grad, grid_level, grid_max_pts, &
                                grid_response, &
                                timing, grid, xg)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nspin
      real(dp), intent(in) :: density(b%nao, b%nao, nspin)
      real(dp), intent(in) :: bcoef, ccoef
      real(dp), intent(inout) :: grad(3, b%natm)
      integer, intent(in), optional :: grid_level, grid_max_pts
      !! NO `rho_tol`: unlike xc_grad_add's semilocal term, this
      !! module always screens at `VV10_RHO_TOL` (trc_vv10's own
      !! constant), everywhere -- the gather, the batch screen below, AND
      !! `vv10_grad_kernel_body`'s own `ri < VV10_RHO_TOL` cutoff, which
      !! has no caller-supplied tolerance threaded into it (it is a
      !! `routine seq` body called from a `do concurrent`/GPU kernel, not
      !! worth a plumbed-through argument for a value every OTHER screen
      !! in this module already hard-codes). A `rho_tol` argument that
      !! moved three of the four screens and not the fourth was strictly
      !! worse than none at all. Matches the energy path exactly:
      !! `trc_vv10_rks`/`trc_vv10_uks` take no `rho_tol` either.
      logical, intent(in), optional :: grid_response
      !! ON by default, same meaning and default as xc_grad_add's own
      !! argument (that module's header, GRID RESPONSE). OFF recovers the
      !! Pulay-only term -- the piece check_dft_grad_vs_pyscf.F90 compares
      !! against pyscf, whose own nlc gradient has no grid-response support
      !! to match against (this module's own header, GRID RESPONSE caveat
      !! in trc_scf_grad.F90's docs).
      type(step_table_t), intent(inout), optional, target :: timing
      !! As xc_grad_add's: device-synced time per step, added into,
      !! initialised by the caller with `xc_grad_step_names()`.
      type(dft_grid_t), intent(in), optional, target :: grid
      type(trc_xc_grid_t), intent(in), optional, target :: xg
      !! The SCF's own VV10 (nlc) grid, already batched and device-
      !! resident -- present TOGETHER, same convention and same reuse
      !! this module's header (GRID) and xc_grad_add's own `grid`/`xg`
      !! document; see trc_scf_context.F90's `trc_grid_cache_t`.

      type(dft_grid_t), target :: grid_own
      type(trc_xc_grid_t), target :: xg_own
      type(dft_grid_t), pointer :: gp
      type(trc_xc_grid_t), pointer :: xgp
      logical :: reuse_grid
      type(step_table_t), target :: tbl_own
      type(step_table_t), pointer :: tbl
      integer, allocatable :: bact(:)
      integer :: kid0(XC_MAX_KERNELS)
      real(dp) :: coef0(XC_MAX_KERNELS)
      type(timer_type) :: tc
      real(dp), allocatable :: fxc0(:), pcoef(:)
      type(error_t) :: err
      integer :: glevel, gmaxpts, ish, i, j, ib, npb, nloc, g, own
      real(dp) :: kvv, beta
      integer, allocatable :: ao_atom(:), at_elem(:)
      real(dp), allocatable :: dtot(:, :, :)
      real(dp), allocatable :: rho(:), grho(:, :), sigma(:)
      real(dp), allocatable :: vrho(:), vsigma(:), dedw(:), fexp(:, :)
      real(dp), allocatable :: gx(:), gy(:), gz(:), gkp(:), gw0p(:), grpw(:)
      integer, allocatable :: gidx(:)
      integer :: npass, npts
      real(dp), allocatable :: grad_extra(:, :)
      logical :: gresp

      glevel = DEFAULT_GRID_LEVEL
      if (present(grid_level)) glevel = grid_level
      gmaxpts = VV10_GRAD_DEFAULT_MAX_PTS
      if (present(grid_max_pts)) gmaxpts = grid_max_pts
      gresp = .true.
      if (present(grid_response)) gresp = grid_response

      if (present(timing)) then
         tbl => timing
      else
         call tbl_own%init("vv10 grad", xc_grad_step_names())
         tbl => tbl_own
      end if

      reuse_grid = present(grid) .and. present(xg)
      call tbl%step_on(lp, ST_GRID, tc)
      allocate (at_elem(b%natm), source=b%at_elem)
      if (reuse_grid) then
         ! The SCF's own nlc grid, already batched and on the device -- no
         ! rebuild, no re-upload. See this subroutine's own header, GRID.
         gp => grid
         xgp => xg
      else
         call dft_grid_build(b%at_r, at_elem, grid_own, err, level=glevel)
         if (err%has_error()) error stop "trc_vv10_grad: "//trim(err%get_message())
         call xg_own%build(grid_own%n_points, grid_own%coords, grid_own%weights, b, max_pts=gmaxpts, &
                           atom=grid_own%atom, qw=grid_own%quad_weights)
         call xg_own%to_device()
         gp => grid_own
         xgp => xg_own
      end if
      npts = xgp%npts
      call tbl%step_off(lp, ST_GRID, tc)

      if (npts > 0) then
         allocate (ao_atom(b%nao))
         do ish = 1, b%nshell
            do i = 0, (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2 - 1
               ao_atom(b%sh_ao(ish) + i) = b%sh_at(ish)
            end do
         end do
         !$acc enter data copyin(ao_atom)

         ! VV10 sees the TOTAL density only (this module's own header) --
         ! built once, at nspin = 1, so the semilocal gradient's pipeline it
         ! reuses (xc_grad_drive, written for a spin-indexed
         ! `density(nao,nao,nspin)`) can be handed this array directly with
         ! nspin = 1 rather than reshaped.
         allocate (dtot(b%nao, b%nao, 1))
         !$acc enter data create(dtot)
         !$acc data copyin(density)
         if (nspin == 1) then
            do concurrent(i=1:b%nao, j=1:b%nao)
               dtot(i, j, 1) = density(i, j, 1)
            end do
         else
            do concurrent(i=1:b%nao, j=1:b%nao)
               dtot(i, j, 1) = density(i, j, 1) + density(i, j, 2)
            end do
         end if
         !$acc end data

         kvv = bcoef*1.5_dp*PI*(9.0_dp*PI)**(-1.0_dp/6.0_dp)
         beta = (3.0_dp/(bcoef*bcoef))**0.75_dp/32.0_dp

         allocate (rho(npts), grho(3, npts), sigma(npts))
         allocate (vrho(npts), vsigma(npts), dedw(npts), fexp(3, npts))
         !$acc data copyin(dtot) create(rho, grho, sigma, vrho, vsigma, dedw, fexp)

         call tbl%step_on(lp, ST_POINTS, tc)
         ! PASS 1: rho/grad-rho/sigma at every outer point -- trc_vv10's own
         ! tiled driver, reused bit-identically (this module's header).
         call vv10_density_tile(b%nao, b%nshell, b%maxnp, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                                xgp%npts, xgp%nbatch, xgp%r, xgp%batch_of, xgp%b_off, xgp%b_shoff, &
                                size(xgp%b_sh), xgp%b_sh, xgp%b_shao, xgp%b_aooff, size(xgp%b_ao), xgp%b_ao, &
                                dtot(:, :, 1), VV10_GRAD_TILE_CAP, lp, rho, grho, sigma)

         ! The inner-grid survivor gather (host-side compaction, then an
         ! O(N) device gather) -- trc_vv10's own `vv10_gather`, see
         ! vv10_drive's own comment on why this runs on the host first.
         !$acc update self(rho)
         allocate (gidx(max(npts, 1)))
         npass = 0
         do i = 1, npts
            if (rho(i) >= VV10_RHO_TOL) then
               npass = npass + 1
               gidx(npass) = i
            end if
         end do
         allocate (gx(max(npass, 1)), gy(max(npass, 1)), gz(max(npass, 1)))
         allocate (gkp(max(npass, 1)), gw0p(max(npass, 1)), grpw(max(npass, 1)))
         !$acc data copyin(gidx(1:max(npass, 1))) create(gx, gy, gz, gkp, gw0p, grpw)
         call vv10_gather(npts, xgp%r, xgp%w, rho, sigma, ccoef, kvv, npass, gidx, gx, gy, gz, gkp, gw0p, grpw)

         ! PASS 2, gradient variant: vrho/vsigma (fed into the Pulay
         ! contraction below) plus dedw/fexp (the two nonlocal-only grid-
         ! response pieces, this module's own header) -- see
         ! vv10_grad_kernel_body's own comment for why this duplicates
         ! trc_vv10's PASS 2 rather than extending it in place.
         call vv10_grad_kernel(npts, xgp%r, xgp%w, rho, sigma, ccoef, kvv, beta, npass, &
                               gx, gy, gz, gkp, gw0p, grpw, vrho, vsigma, dedw, fexp)
         !$acc end data
         deallocate (gidx, gx, gy, gz, gkp, gw0p, grpw)
         call tbl%step_off(lp, ST_POINTS, tc)

         ! `dedw`/`fexp` are only ever read host-side (the weight-response
         ! block below); `rho`/`grho`/`sigma`/`vrho`/`vsigma` stay device-
         ! resident inside the open `!$acc data` region for the batch loop.
         !$acc update self(dedw, fexp)

         allocate (grad_extra(3, b%natm))
         grad_extra = 0.0_dp

         ! The Pulay and point-motion terms: the semilocal gradient's own
         ! chunked pipeline (trc_xc_grad's xc_grad_drive), restricted GGA,
         ! with this module's PASS 2 potential in place of a functional and
         ! its own batch screen (any point at or above VV10_RHO_TOL).
         allocate (bact(xgp%nbatch), fxc0(npts))
         fxc0 = 0.0_dp
         do ib = 1, xgp%nbatch
            npb = xgp%b_off(ib + 1) - xgp%b_off(ib)
            nloc = xgp%b_aooff(ib + 1) - xgp%b_aooff(ib)
            bact(ib) = 0
            if (npb == 0 .or. nloc == 0) cycle
            if (vv10_grad_batch_active(xgp%b_off(ib), npb, rho, npts)) bact(ib) = 1
         end do
         kid0 = 0
         coef0 = 0.0_dp
         call xc_grad_drive(b, lp, .true., 1, 0, kid0, coef0, dtot, xgp, 0.0_dp, ao_atom, grad, gresp, &
                            fxc0, tbl, act_in=bact, ext_grho=grho, ext_vrho=vrho, ext_vsigma=vsigma)

         !$acc end data
         if (gresp) then
            ! GRID RESPONSE, the weight-derivative and the point's-own-force
            ! pieces (this module's header, THE NONLOCAL EXTRAS), over the
            ! batches the screen kept: the point's own explicit force at its
            ! FULL weight onto its owner, here; the partition's response to
            ! qw*dedw, both its pieces, on the device.
            call tbl%step_on(lp, ST_PDERIV, tc)
            allocate (pcoef(npts))
            pcoef = 0.0_dp
            do ib = 1, xgp%nbatch
               if (bact(ib) == 0) cycle
               do g = xgp%b_off(ib), xgp%b_off(ib + 1) - 1
                  own = xgp%atom(g)
                  grad_extra(:, own) = grad_extra(:, own) + xgp%w(g)*fexp(:, g)
                  pcoef(g) = xgp%qw(g)*dedw(g)
               end do
            end do
            call dft_partition_becke_grad_contract(xgp%r, b%at_r, at_elem, xgp%atom, pcoef, gp%partition, &
                                                   grad_extra, err)
            if (err%has_error()) error stop "trc_vv10_grad: "//trim(err%get_message())
            deallocate (pcoef)
            call tbl%step_off(lp, ST_PDERIV, tc)
            call tbl%step_on(lp, ST_RESP, tc)
            !$acc update self(grad)
            grad = grad + grad_extra
            !$acc update device(grad)
            call tbl%step_off(lp, ST_RESP, tc)
         end if
         deallocate (grad_extra, bact, fxc0)

         !$acc exit data delete(ao_atom, dtot)
         deallocate (ao_atom, dtot, rho, grho, sigma, vrho, vsigma, dedw, fexp)
      end if

      if (.not. reuse_grid) then
         call xg_own%release()
         call grid_own%destroy()
      end if
      deallocate (at_elem)
   end subroutine vv10_grad_add

   !
   ! Host-side screen: a batch whose every point's density is below
   ! threshold contributes nothing to the Pulay/grid-response contraction --
   ! `rho` is already the host array the gather above just read, so this is
   ! a cheap O(npb) bound, not a new array.
   !
   logical function vv10_grad_batch_active(g0, npb, rho, npts) result(active)
      integer, intent(in) :: g0, npb, npts
      real(dp), intent(in) :: rho(npts)
      integer :: g
      active = .false.
      do g = g0, g0 + npb - 1
         if (rho(g) >= VV10_RHO_TOL) then
            active = .true.
            return
         end if
      end do
   end function vv10_grad_batch_active

   !
   ! PASS 2, gradient variant: `vrho(i)`/`vsigma(i)` (the ordinary VV10
   ! potential, fed into the SAME GGA Pulay contraction xc_grad_add's
   ! semilocal term uses) plus `dedw(i)`/`fexp(:,i)`, the two pieces with no
   ! semilocal counterpart -- see this module's header, THE NONLOCAL EXTRAS.
   ! A duplicate of trc_vv10's own `vv10_kernel_body`, not a call to it:
   ! that body has no `dedw`/`fexp` outputs (see this module's header,
   ! ENERGY-SIDE REUSE).
   !
#ifdef TRC_DEVICE
   subroutine vv10_grad_kernel(npts, r, w, rho, sigma, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                               vrho, vsigma, dedw, fexp)
      integer, intent(in) :: npts, npass
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: c, kvv, beta
      real(dp), intent(in) :: gx(npass), gy(npass), gz(npass), gkp(npass), gw0p(npass), grpw(npass)
      real(dp), intent(out) :: vrho(npts), vsigma(npts), dedw(npts), fexp(3, npts)
      integer :: i
      do concurrent(i=1:npts)
         call vv10_grad_kernel_body(i, npts, r, w, rho, sigma, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                                    vrho(i), vsigma(i), dedw(i), fexp(1, i), fexp(2, i), fexp(3, i))
      end do
   end subroutine vv10_grad_kernel
#else
   subroutine vv10_grad_kernel(npts, r, w, rho, sigma, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                               vrho, vsigma, dedw, fexp)
      integer, intent(in) :: npts, npass
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: c, kvv, beta
      real(dp), intent(in) :: gx(npass), gy(npass), gz(npass), gkp(npass), gw0p(npass), grpw(npass)
      real(dp), intent(out) :: vrho(npts), vsigma(npts), dedw(npts), fexp(3, npts)
      integer :: i
      ! Plain sequential `do`, not `do concurrent` -- the same ifx
      ! mis-optimisation trc_vv10.F90's own PASS 2 header documents (its
      ! `vv10_kernel`), for the identical pattern: an output built from this
      ! body's own `do j` reduction, inside a `do concurrent`. Found there by
      ! bisection; applied here pre-emptively rather than re-discovered.
      do i = 1, npts
         call vv10_grad_kernel_body(i, npts, r, w, rho, sigma, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                                    vrho(i), vsigma(i), dedw(i), fexp(1, i), fexp(2, i), fexp(3, i))
      end do
   end subroutine vv10_grad_kernel
#endif

   !
   ! `i`'s own pair sum, matched term for term against pyscf's `_vv10nlc`.
   ! `w(i)` plays no role here (the weight multiplies `vrho`/`vsigma`/
   ! `dedw`/`fexp` only once they are CONTRACTED, same as trc_vv10's own
   ! PASS 2) -- the convention that matters is `dx = inner - outer`
   ! (`gx(j) - r(1,i)`, NOT the sign trc_vv10's own `vv10_kernel_body`
   ! uses for its unsigned `r2`), since `fexp` is LINEAR in `dx`, not
   ! squared.
   !
   pure subroutine vv10_grad_kernel_body(i, npts, r, w, rho, sigma, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                                         vrho_i, vsigma_i, dedw_i, fx_i, fy_i, fz_i)
      !$acc routine seq
      integer, intent(in) :: i, npts, npass
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: c, kvv, beta
      real(dp), intent(in) :: gx(npass), gy(npass), gz(npass), gkp(npass), gw0p(npass), grpw(npass)
      real(dp), intent(out) :: vrho_i, vsigma_i, dedw_i, fx_i, fy_i, fz_i

      integer :: j
      real(dp) :: ri, gi, w0tmp, w0, dw0dr, dw0dg, kk, dkdr
      real(dp) :: dx, dy, dz, r2
      real(dp) :: gg, ggp, gt, tt, t2, fsum, usum, wsum, ff, qq
      real(dp) :: fsx, fsy, fsz

      vrho_i = 0.0_dp
      vsigma_i = 0.0_dp
      dedw_i = 0.0_dp
      fx_i = 0.0_dp
      fy_i = 0.0_dp
      fz_i = 0.0_dp
      ri = rho(i)
      if (ri < VV10_RHO_TOL) return
      gi = sigma(i)
      w0tmp = c*(gi/(ri*ri))**2
      w0 = sqrt(w0tmp + PI43*ri)
      dw0dr = (0.5_dp*PI43*ri - 2.0_dp*w0tmp)/w0
      dw0dg = c*gi/(ri**3*w0)
      kk = kvv*ri**(1.0_dp/6.0_dp)
      dkdr = kk/6.0_dp

      fsum = 0.0_dp
      usum = 0.0_dp
      wsum = 0.0_dp
      fsx = 0.0_dp
      fsy = 0.0_dp
      fsz = 0.0_dp
      do j = 1, npass
         dx = gx(j) - r(1, i)
         dy = gy(j) - r(2, i)
         dz = gz(j) - r(3, i)
         r2 = dx*dx + dy*dy + dz*dz
         ggp = r2*gw0p(j) + gkp(j)
         gg = r2*w0 + kk
         gt = gg + ggp
         tt = grpw(j)/(gg*ggp*gt)
         fsum = fsum + tt
         t2 = tt*(1.0_dp/gg + 1.0_dp/gt)
         usum = usum + t2
         wsum = wsum + t2*r2
         ! The kernel's own dependence on the separation -- differentiated
         ! once more than `t2` above, hence the extra `1/gg`/`1/ggp`/`1/gt`
         ! factor rather than a sum of two.
         qq = tt*(w0/gg + gw0p(j)/ggp + (w0 + gw0p(j))/gt)
         fsx = fsx + qq*dx
         fsy = fsy + qq*dy
         fsz = fsz + qq*dz
      end do
      ff = -1.5_dp*fsum
      vrho_i = beta + ff + 1.5_dp*(usum*dkdr + wsum*dw0dr)
      vsigma_i = 1.5_dp*wsum*dw0dg
      ! dE/dw_i: `beta + F`, not `beta + F/2` -- see this module's header.
      dedw_i = ri*(beta + ff)
      ! d/dr_i of the pair sum: the factor of 3 (not 1.5) is the symmetry
      ! this module's header derives -- one pass over `i` already carries
      ! both its own-point and every-other-point's-inner-point role.
      fx_i = -3.0_dp*ri*fsx
      fy_i = -3.0_dp*ri*fsy
      fz_i = -3.0_dp*ri*fsz
   end subroutine vv10_grad_kernel_body

end module trc_vv10_grad
