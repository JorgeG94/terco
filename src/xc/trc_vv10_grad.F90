!
! VV10 (Vydrov & Van Voorhis 2010) nonlocal correlation's contribution to the
! nuclear gradient, RKS and UKS, with grid response. Added into an existing
! `grad(3,natm)`, the same convention trc_xc_grad%xc_gradient_add's semilocal
! term uses.
!
! THE FORMULA -- metalquicha's own derivation (mqc_czt_gradient.f90's
! `vv10_gradient_core`, this module's reference implementation), ported term
! for term. The non-locality is already inside vrho/vsigma: E_nl is a
! functional of rho and sigma alone,
!
!     E_nl = int dr rho(r) [ beta + 1/2 int dr' rho(r') Phi(r,r') ]
!     dE_nl/dR = int dr [ vrho(r) drho(r)/dR + vsigma(r) dsigma(r)/dR ]
!
! and the O(N_grid^2) pair sum (trc_vv10's PASS 2) is spent entirely on
! producing the two arrays `vrho`/`vsigma` at every point, exactly the way
! trc_vv10's own energy path does. What is left is the ORDINARY GGA
! contraction trc_xc_grad.F90 already implements for a semilocal kernel's
! vrho/vsigma -- this module reuses that contraction machinery
! (xc_grad_collocate2, xc_grad_gather, xc_grad_scatter, xc_grad_scatter_points,
! all exported from trc_xc_grad for exactly this), just reading `vrho(g)`/
! `vsigma(g)` from VV10's own precomputed arrays instead of calling
! `xc_eval_point`.
!
! THE NONLOCAL EXTRAS. Two things have no counterpart in a semilocal
! functional's gradient, and neither is optional (worth 4.5e-4 on water,
! metalquicha's own measurement, and it does not shrink with the grid):
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
!     weight, same as metalquicha's own `ctx%nlc_grid%weights(gg)*fexp(:,gg)`.
!
! Both grid-response pieces close the SAME header gap xc_gradient_add's own
! "GRID RESPONSE" section documents for the semilocal term (point motion of
! both grids, partition-weight motion, and their cross-cancellation via
! `becke_partition_derivatives`) -- this module's own weight/point-motion
! block is that section's VV10 analogue, with `dedw`/`fexp` in place of
! `fxc = rho eps`/nothing.
!
! ENERGY-SIDE REUSE. `vv10_density_tile` (PASS 1: rho/grad-rho/sigma at every
! outer point) and `vv10_gather` (the inner-grid survivor list `vv10_kernel`
! reads) are trc_vv10's own, exported unchanged -- the gradient's outer-point
! sweep is bit-identical O(N) work to the energy path's, so it is not
! re-derived here. PASS 2 itself (`vv10_kernel_body`) is NOT reused as-is:
! it has no `dedw`/`fexp` outputs, and adding them to a validated,
! performance-sensitive energy kernel under time pressure was judged a worse
! trade than a second, gradient-only pass-2 body that reads the exact same
! gathered arrays and the exact same `kvv`/`beta`/`w0`/`k` formulas
! (`vv10_grad_kernel_body` below) -- the two are checked against each other
! implicitly, since both are checked against pyscf/metalquicha independently.
!
! RESTRICTED AND UNRESTRICTED. VV10 is a functional of the TOTAL density
! alone (trc_vv10.F90's own header), so this module's single entry point
! takes `dens(nao,nao,nspin)` exactly as xc_gradient_add does, forms
! `dtot = dens(:,:,1) [+ dens(:,:,2)]` once, and runs the entire gradient
! against that one density -- there is no separate alpha/beta contraction to
! carry, unlike the semilocal GGA gradient's cross-spin Q2/rmatx terms.
!
! GRID. This module builds its OWN quadrature at `b`'s geometry, the same
! reason xc_gradient_add does (trc_scf_context keeps neither grid nor
! functional alive past convergence) -- at the CALLER's chosen level
! (`grid_level`, expected to be the SCF's own `nlc_grid_level`, not the
! semilocal grid's `grid_level`: VV10 runs on its own, coarser grid
! end-to-end, energy and gradient alike), with `atom`/`qw` attached so the
! weight-response block below can read them -- the same two optional
! arguments xc_gradient_add passes its own grid build.
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
! GPU. Not yet measured: see this branch's own report. Built from the same
! idioms (`do concurrent`/`routine seq`, the `_OPENACC` CPU/GPU pass-2 split
! trc_vv10.F90's own header explains) as the two modules it reuses.
!
module trc_vv10_grad
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_grad, only: xc_grad_collocate2, xc_grad_gather, xc_grad_scatter, xc_grad_scatter_points
   use trc_vv10, only: vv10_density_tile, vv10_gather, VV10_RHO_TOL, PI43
   use trc_collocation, only: NHESS
   use trc_linalg, only: trc_linalg_t
   use trc_dft_grid, only: dft_grid_t, build_dft_grid, DEFAULT_GRID_LEVEL
   use trc_dft_partition, only: becke_partition_derivatives
   use trc_error, only: error_t
   implicit none
   private

   public :: vv10_gradient_add

   real(dp), parameter :: PI = 3.14159265358979323846_dp

contains

   !
   ! The one entry point. `dens(b%nao,b%nao,nspin)` -- the converged density,
   ! per spin channel, same convention xc_gradient_add uses. `bcoef`/`ccoef`
   ! are VV10's own `b`/`C` (6.0, 0.01 for wB97X-V). ADDS into `grad`.
   !
   subroutine vv10_gradient_add(b, lp, nspin, dens, bcoef, ccoef, grad, grid_level, grid_max_pts, rho_tol, &
                                grid_response)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nspin
      real(dp), intent(in) :: dens(b%nao, b%nao, nspin)
      real(dp), intent(in) :: bcoef, ccoef
      real(dp), intent(inout) :: grad(3, b%natm)
      integer, intent(in), optional :: grid_level, grid_max_pts
      real(dp), intent(in), optional :: rho_tol
      logical, intent(in), optional :: grid_response
      !! ON by default, same meaning and default as xc_gradient_add's own
      !! argument (that module's header, GRID RESPONSE). OFF recovers the
      !! Pulay-only term -- the piece check_dft_grad_vs_pyscf.F90 compares
      !! against pyscf, whose own nlc gradient has no grid-response support
      !! to match against (this module's own header, GRID RESPONSE caveat
      !! in trc_scf_grad.F90's docs).

      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(error_t) :: err
      integer :: glevel, gmaxpts, ish, i, j, ib, npb, nloc
      real(dp) :: rtol, kvv, beta
      integer, allocatable :: ao_atom(:), at_z(:)
      real(dp), allocatable :: dtot(:, :, :)
      real(dp), allocatable :: rho(:), grho(:, :), sigma(:)
      real(dp), allocatable :: vrho(:), vsigma(:), dedw(:), fexp(:, :)
      real(dp), allocatable :: gx(:), gy(:), gz(:), gkp(:), gw0p(:), grpw(:)
      integer, allocatable :: gidx(:)
      integer :: npass, npts
      real(dp), allocatable :: grad_extra(:, :)
      logical :: active, gresp

      glevel = DEFAULT_GRID_LEVEL
      if (present(grid_level)) glevel = grid_level
      gmaxpts = 512
      if (present(grid_max_pts)) gmaxpts = grid_max_pts
      rtol = VV10_RHO_TOL
      if (present(rho_tol)) rtol = rho_tol
      gresp = .true.
      if (present(grid_response)) gresp = grid_response

      allocate (at_z(b%natm), source=nint(b%at_z))
      call build_dft_grid(b%at_r, at_z, grid, err, level=glevel)
      if (err%has_error()) error stop "trc_vv10_grad: "//trim(err%get_message())
      call xg%build(grid%n_points, grid%coords, grid%weights, b, max_pts=gmaxpts, &
                    atom=grid%atom, qw=grid%quad_weights)
      call xg%to_device()
      npts = xg%npts

      if (npts > 0) then
         allocate (ao_atom(b%nao))
         do ish = 1, b%nshell
            do i = 0, (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2 - 1
               ao_atom(b%sh_ao(ish) + i) = b%sh_at(ish)
            end do
         end do
         !$acc enter data copyin(ao_atom)

         ! VV10 sees the TOTAL density only (this module's own header) --
         ! built once, at nspin = 1, so every downstream call (the generic
         ! xc_grad_gather/xc_grad_scatter it reuses, both written for a
         ! spin-indexed `dens(nao,nao,nspin)`) can be handed this array
         ! directly with nspin = 1 rather than reshaped.
         allocate (dtot(b%nao, b%nao, 1))
         !$acc enter data create(dtot)
         !$acc data copyin(dens)
         if (nspin == 1) then
            do concurrent(i=1:b%nao, j=1:b%nao)
               dtot(i, j, 1) = dens(i, j, 1)
            end do
         else
            do concurrent(i=1:b%nao, j=1:b%nao)
               dtot(i, j, 1) = dens(i, j, 1) + dens(i, j, 2)
            end do
         end if
         !$acc end data

         kvv = bcoef*1.5_dp*PI*(9.0_dp*PI)**(-1.0_dp/6.0_dp)
         beta = (3.0_dp/(bcoef*bcoef))**0.75_dp/32.0_dp

         allocate (rho(npts), grho(3, npts), sigma(npts))
         allocate (vrho(npts), vsigma(npts), dedw(npts), fexp(3, npts))
         !$acc data copyin(dtot) create(rho, grho, sigma, vrho, vsigma, dedw, fexp)

         ! PASS 1: rho/grad-rho/sigma at every outer point -- trc_vv10's own
         ! tiled driver, reused bit-identically (this module's header).
         call vv10_density_tile(b%nao, b%nshell, b%maxnp, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                                xg%npts, xg%nbatch, xg%r, xg%batch_of, xg%b_off, xg%b_shoff, &
                                size(xg%b_sh), xg%b_sh, xg%b_shao, xg%b_aooff, size(xg%b_ao), xg%b_ao, &
                                dtot(:, :, 1), 8388608_int64, lp, rho, grho, sigma)

         ! The inner-grid survivor gather (host-side compaction, then an
         ! O(N) device gather) -- trc_vv10's own `vv10_gather`, see
         ! vv10_drive's own comment on why this runs on the host first.
         !$acc update self(rho)
         allocate (gidx(max(npts, 1)))
         npass = 0
         do i = 1, npts
            if (rho(i) >= rtol) then
               npass = npass + 1
               gidx(npass) = i
            end if
         end do
         allocate (gx(max(npass, 1)), gy(max(npass, 1)), gz(max(npass, 1)))
         allocate (gkp(max(npass, 1)), gw0p(max(npass, 1)), grpw(max(npass, 1)))
         !$acc data copyin(gidx(1:max(npass, 1))) create(gx, gy, gz, gkp, gw0p, grpw)
         call vv10_gather(npts, xg%r, xg%w, rho, sigma, ccoef, kvv, npass, gidx, gx, gy, gz, gkp, gw0p, grpw)

         ! PASS 2, gradient variant: vrho/vsigma (fed into the Pulay
         ! contraction below) plus dedw/fexp (the two nonlocal-only grid-
         ! response pieces, this module's own header) -- see
         ! vv10_grad_kernel_body's own comment for why this duplicates
         ! trc_vv10's PASS 2 rather than extending it in place.
         call vv10_grad_kernel(npts, xg%r, xg%w, rho, sigma, ccoef, kvv, beta, npass, &
                               gx, gy, gz, gkp, gw0p, grpw, vrho, vsigma, dedw, fexp)
         !$acc end data
         deallocate (gidx, gx, gy, gz, gkp, gw0p, grpw)

         ! `dedw`/`fexp` are only ever read host-side (the weight-response
         ! block below); `rho`/`grho`/`sigma`/`vrho`/`vsigma` stay device-
         ! resident inside the open `!$acc data` region for the batch loop.
         !$acc update self(dedw, fexp)

         allocate (grad_extra(3, b%natm))
         grad_extra = 0.0_dp

         do ib = 1, xg%nbatch
            npb = xg%b_off(ib + 1) - xg%b_off(ib)
            nloc = xg%b_aooff(ib + 1) - xg%b_aooff(ib)
            if (npb == 0 .or. nloc == 0) cycle
            active = vv10_grad_batch_active(xg%b_off(ib), npb, rho, npts, rtol)
            if (.not. active) cycle
            call vv10_grad_one_batch(b, lp, dtot, xg, ib, npb, nloc, ao_atom, grad, &
                                     b%at_r, at_z, grid%scheme, grid%adjust, gresp, &
                                     grho, vrho, vsigma, dedw, fexp, grad_extra)
         end do

         !$acc end data
         if (gresp) then
            !$acc update self(grad)
            grad = grad + grad_extra
            !$acc update device(grad)
         end if
         deallocate (grad_extra)

         !$acc exit data delete(ao_atom, dtot)
         deallocate (ao_atom, dtot, rho, grho, sigma, vrho, vsigma, dedw, fexp)
      end if

      call xg%release()
      call grid%destroy()
      deallocate (at_z)
   end subroutine vv10_gradient_add

   !
   ! Host-side screen: a batch whose every point's density is below
   ! threshold contributes nothing to the Pulay/grid-response contraction --
   ! `rho` is already the host array the gather above just read, so this is
   ! a cheap O(npb) bound, not a new array.
   !
   logical function vv10_grad_batch_active(g0, npb, rho, npts, rtol) result(active)
      integer, intent(in) :: g0, npb, npts
      real(dp), intent(in) :: rho(npts)
      real(dp), intent(in) :: rtol
      integer :: g
      active = .false.
      do g = g0, g0 + npb - 1
         if (rho(g) >= rtol) then
            active = .true.
            return
         end if
      end do
   end function vv10_grad_batch_active

   !
   ! PASS 2, gradient variant: `vrho(i)`/`vsigma(i)` (the ordinary VV10
   ! potential, fed into the SAME GGA Pulay contraction xc_gradient_add's
   ! semilocal term uses) plus `dedw(i)`/`fexp(:,i)`, the two pieces with no
   ! semilocal counterpart -- see this module's header, THE NONLOCAL EXTRAS.
   !
#ifdef _OPENACC
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
   ! `i`'s own pair sum, matched term for term against metalquicha's
   ! `vv10_nlc` (want_grad branch) and pyscf's `_vv10nlc`: `w(i)` plays no
   ! role here (the weight multiplies `vrho`/`vsigma`/`dedw`/`fexp` only once
   ! they are CONTRACTED, same as trc_vv10's own PASS 2) -- the convention
   ! that matters is `dx = inner - outer` (`gx(j) - r(1,i)`, NOT the sign
   ! trc_vv10's own `vv10_kernel_body` uses for its unsigned `r2`), since
   ! `fexp` is LINEAR in `dx`, not squared.
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

   !
   ! One batch of the Pulay + grid-response-point-motion contraction: the
   ! same shape as trc_xc_grad.F90's own `xc_grad_one_batch` GGA branch
   ! (nspin = 1, since VV10 only ever sees the total density), with
   ! `vrho`/`vsigma` read from this module's own PASS 2 output instead of
   ! `xc_eval_point`, and the weight/point-motion grid-response block built
   ! from `dedw`/`fexp` instead of `fxc = rho eps`.
   !
   subroutine vv10_grad_one_batch(b, lp, dtot, xg, ib, npb, nloc, ao_atom, grad, at_r, at_z, scheme, adjust, &
                                  gresp, grho, vrho, vsigma, dedw, fexp, grad_extra)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: dtot(b%nao, b%nao, 1)
      type(trc_xc_grid_t), intent(in) :: xg
      integer, intent(in) :: ib, npb, nloc
      integer, intent(in) :: ao_atom(b%nao)
      real(dp), intent(inout) :: grad(3, b%natm)
      real(dp), intent(in) :: at_r(3, b%natm)
      integer, intent(in) :: at_z(b%natm)
      integer, intent(in) :: scheme, adjust
      logical, intent(in) :: gresp
      real(dp), intent(in) :: grho(3, xg%npts), vrho(xg%npts), vsigma(xg%npts)
      real(dp), intent(in) :: dedw(xg%npts), fexp(3, xg%npts)
      real(dp), intent(inout) :: grad_extra(3, b%natm)

      real(dp), allocatable :: chi(:, :), gchi(:, :, :), hchi(:, :, :)
      real(dp), allocatable :: dloc(:, :, :), xmat(:, :, :)
      real(dp), allocatable :: pall(:, :, :), q1all(:, :, :), rall(:, :, :), rmat(:, :, :)
      real(dp), allocatable :: tall(:, :, :), qmove(:, :)
      real(dp), allocatable :: dweights(:, :, :)
      type(error_t) :: err
      integer :: g, a, own, g0
      real(dp) :: contrib

      allocate (chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS))
      allocate (dloc(nloc, nloc, 1), xmat(npb, nloc, 1))
      allocate (pall(npb, 3*nloc, 1), q1all(npb, 3*nloc, 1), rall(npb, nloc, 1), rmat(npb, nloc, 1))
      allocate (tall(3*nloc, nloc, 1), qmove(npb, 3))
      !$acc enter data create(chi, gchi, hchi, dloc, xmat, pall, q1all, rall, rmat, tall, qmove)

      call xc_grad_collocate2(b%nshell, b%maxnp, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                              xg%npts, xg%r, xg%nbatch, xg%b_off, xg%b_shoff, size(xg%b_sh), &
                              xg%b_sh, xg%b_shao, ib, npb, nloc, chi, gchi, hchi)

      call xc_grad_gather(b%nao, size(xg%b_ao), xg%b_ao, xg%b_aooff(ib), nloc, 1, dtot, dloc)

      call lp%gemm_mat(chi, dloc(:, :, 1), xmat(:, :, 1))

      call vv10_grad_point(npb, xg%npts, xg%w, xg%b_off(ib), nloc, grho, vrho, vsigma, &
                           chi, gchi, hchi, xmat(:, :, 1), pall(:, :, 1), q1all(:, :, 1), rall(:, :, 1), qmove)

      ! GRID RESPONSE, point-motion piece, GGA addition (needs every
      ! point's `rall` density-contracted) -- xc_grad_pointmotion_gga's
      ! exact counterpart; see this module's header for why this is a
      ! second GEMM rather than foldable into vv10_grad_point's own kernel.
      ! Skipped, along with the scatter below, when `gresp` is OFF -- same
      ! as xc_grad_one_batch's own `if (gresp)` guard around the identical
      ! pair.
      if (gresp) then
         call lp%gemm_mat(rall(:, :, 1), dloc(:, :, 1), rmat(:, :, 1))
         call vv10_grad_pointmotion(npb, nloc, q1all(:, :, 1), rmat(:, :, 1), qmove)
         call xc_grad_scatter_points(xg%npts, xg%atom, xg%b_off(ib), npb, qmove, b%natm, grad)
      end if

      call lp%gemm_mat(pall(:, :, 1), chi, tall(:, :, 1), transa="T", beta=0.0_dp)
      call lp%gemm_mat(q1all(:, :, 1), rall(:, :, 1), tall(:, :, 1), transa="T", beta=1.0_dp)

      call xc_grad_scatter(ib, xg%nbatch, xg%b_aooff, size(xg%b_ao), xg%b_ao, ao_atom, b%nao, &
                           nloc, 1, dloc, tall, b%natm, grad)

      ! GRID RESPONSE, weight-derivative and point's-own-force pieces --
      ! HOST side (dedw/fexp are already host-resident, this module's own
      ! call site); see this module's header, THE NONLOCAL EXTRAS. Skipped
      ! entirely when `gresp` is OFF, same as xc_gradient_add's own
      ! weight-response block.
      if (gresp) then
         g0 = xg%b_off(ib)
         allocate (dweights(3, b%natm, npb))
         call becke_partition_derivatives(xg%r(:, g0:g0 + npb - 1), at_r, at_z, &
                                          xg%atom(g0:g0 + npb - 1), scheme, adjust, dweights, err)
         if (err%has_error()) error stop "trc_vv10_grad: "//trim(err%get_message())
         do g = 1, npb
            own = xg%atom(g0 + g - 1)
            ! The point's own explicit force, at the point's FULL weight.
            grad_extra(:, own) = grad_extra(:, own) + xg%w(g0 + g - 1)*fexp(:, g0 + g - 1)
            contrib = xg%qw(g0 + g - 1)*dedw(g0 + g - 1)
            do a = 1, b%natm
               grad_extra(:, a) = grad_extra(:, a) + contrib*dweights(:, a, g)
            end do
            do a = 1, b%natm
               grad_extra(:, own) = grad_extra(:, own) - contrib*dweights(:, a, g)
            end do
         end do
         deallocate (dweights)
      end if

      !$acc exit data delete(chi, gchi, hchi, dloc, xmat, pall, q1all, rall, rmat, tall, qmove)
      deallocate (chi, gchi, hchi, dloc, xmat, pall, q1all, rall, rmat, tall, qmove)
   end subroutine vv10_grad_one_batch

   !
   ! Pulay term, one point per thread: `pall`/`q1all`/`rall` built from the
   ! PRECOMPUTED `vrho(g)`/`vsigma(g)` and the batch's own local collocation,
   ! exactly as xc_grad_point_gga_body's nspin = 1 branch -- see that
   ! routine's own comments for the algebra; the only change is that
   ! `grad_rho` is read from the GLOBAL array PASS 1 already built (`grho`)
   ! rather than reconstructed from `xmat`/`gchi` (the two are numerically
   ! the same quantity: this batch's local AO list is the same one PASS 1's
   ! `vv10_density_tile` screened to, so the local reconstruction a
   ! semilocal kernel needs equals the global value already in hand here).
   !
   subroutine vv10_grad_point(npb, npts, w, p0, nloc, grho, vrho, vsigma, chi, gchi, hchi, xmat, &
                              pall, q1all, rall, qmove)
      integer, intent(in) :: npb, npts, p0, nloc
      real(dp), intent(in) :: w(npts), grho(3, npts), vrho(npts), vsigma(npts)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      real(dp), intent(in) :: xmat(npb, nloc)
      real(dp), intent(inout) :: pall(npb, 3*nloc), q1all(npb, 3*nloc), rall(npb, nloc)
      real(dp), intent(inout) :: qmove(npb, 3)
      integer :: g
      do concurrent(g=1:npb)
         call vv10_grad_point_body(g, npb, npts, w, p0, nloc, grho, vrho, vsigma, chi, gchi, hchi, xmat, &
                                   pall, q1all, rall, qmove)
      end do
   end subroutine vv10_grad_point

   pure subroutine vv10_grad_point_body(g, npb, npts, w, p0, nloc, grho, vrho, vsigma, chi, gchi, hchi, xmat, &
                                        pall, q1all, rall, qmove)
      !$acc routine seq
      integer, intent(in) :: g, npb, npts, p0, nloc
      real(dp), intent(in) :: w(npts), grho(3, npts), vrho(npts), vsigma(npts)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      real(dp), intent(in) :: xmat(npb, nloc)
      real(dp), intent(inout) :: pall(npb, 3*nloc), q1all(npb, 3*nloc), rall(npb, nloc)
      real(dp), intent(inout) :: qmove(npb, 3)
      integer :: mu, k, gg
      real(dp) :: wg, vr, vs, grx, gry, grz
      real(dp) :: hdota(3), qm(3)

      gg = p0 + g - 1
      wg = w(gg)
      vr = vrho(gg)
      vs = vsigma(gg)
      grx = grho(1, gg)
      gry = grho(2, gg)
      grz = grho(3, gg)

      do mu = 1, nloc
         ! (grad_rho . H_mu)_k, hchi packed xx, xy, xz, yy, yz, zz.
         hdota(1) = grx*hchi(g, mu, 1) + gry*hchi(g, mu, 2) + grz*hchi(g, mu, 3)
         hdota(2) = grx*hchi(g, mu, 2) + gry*hchi(g, mu, 4) + grz*hchi(g, mu, 5)
         hdota(3) = grx*hchi(g, mu, 3) + gry*hchi(g, mu, 5) + grz*hchi(g, mu, 6)
         do k = 1, 3
            pall(g, (k - 1)*nloc + mu) = wg*(vr*gchi(g, mu, k) + 2.0_dp*vs*hdota(k))
            q1all(g, (k - 1)*nloc + mu) = wg*2.0_dp*vs*gchi(g, mu, k)
         end do
      end do
      do mu = 1, nloc
         rall(g, mu) = grx*gchi(g, mu, 1) + gry*gchi(g, mu, 2) + grz*gchi(g, mu, 3)
      end do

      ! GRID RESPONSE, point-motion piece, LDA-like part -- see
      ! xc_grad_point_gga_body's own comment on the sign and the "+2"; the
      ! GGA-specific addition needing `rmat` is added on top by
      ! vv10_grad_pointmotion below, once this kernel has run over the
      ! whole batch.
      qm = 0.0_dp
      do mu = 1, nloc
         do k = 1, 3
            qm(k) = qm(k) + pall(g, (k - 1)*nloc + mu)*xmat(g, mu)
         end do
      end do
      qmove(g, 1) = 2.0_dp*qm(1)
      qmove(g, 2) = 2.0_dp*qm(2)
      qmove(g, 3) = 2.0_dp*qm(3)
   end subroutine vv10_grad_point_body

   !
   ! GRID RESPONSE, point-motion piece: the GGA-specific addition to
   ! `qmove` -- xc_grad_pointmotion_gga_body's nspin = 1 branch exactly
   ! (its own header explains `rmat`).
   !
   subroutine vv10_grad_pointmotion(npb, nloc, q1all, rmat, qmove)
      integer, intent(in) :: npb, nloc
      real(dp), intent(in) :: q1all(npb, 3*nloc), rmat(npb, nloc)
      real(dp), intent(inout) :: qmove(npb, 3)
      integer :: g
      do concurrent(g=1:npb)
         call vv10_grad_pointmotion_body(g, npb, nloc, q1all, rmat, qmove)
      end do
   end subroutine vv10_grad_pointmotion

   pure subroutine vv10_grad_pointmotion_body(g, npb, nloc, q1all, rmat, qmove)
      !$acc routine seq
      integer, intent(in) :: g, npb, nloc
      real(dp), intent(in) :: q1all(npb, 3*nloc), rmat(npb, nloc)
      real(dp), intent(inout) :: qmove(npb, 3)
      integer :: mu, k
      real(dp) :: acc(3)
      acc = 0.0_dp
      do mu = 1, nloc
         do k = 1, 3
            acc(k) = acc(k) + q1all(g, (k - 1)*nloc + mu)*rmat(g, mu)
         end do
      end do
      qmove(g, 1) = qmove(g, 1) + 2.0_dp*acc(1)
      qmove(g, 2) = qmove(g, 2) + 2.0_dp*acc(2)
      qmove(g, 3) = qmove(g, 3) + 2.0_dp*acc(3)
   end subroutine vv10_grad_pointmotion_body

end module trc_vv10_grad
