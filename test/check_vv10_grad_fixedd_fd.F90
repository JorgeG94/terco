!
! The VV10 gradient's own formula, isolated from the SCF entirely: a
! central finite difference of E_nlc[D] alone, with D the converged wB97X-V
! AO density matrix held FIXED (not reconverged) while the basis and the
! VV10 grid (points, weights, AND the Becke partition) are rebuilt at every
! displaced geometry.
!
! WHY THIS IS SHARPER THAN check_xc_grad_fd.F90's wB97X-V CASE. That
! program differences the FULLY RECONVERGED KS-DFT energy, so its finite
! difference is dE_SCF/dR = Pulay + grid response + the Hellmann-Feynman
! piece ALREADY folded in by self-consistency + whatever SCF/DIIS
! convergence noise a harder (RSH+NLC) SCF carries from one displaced
! geometry to the next. RESOLVED: that program's own case originally
! measured worst abs ON 4.27e-7, flat across grid levels 3 and 5, and
! Richardson only reached 2.90e-7, not the ~440x collapse a purely
! truncation-limited case gives (w3.xyz's own SVWN case there) -- which
! this program's own clean numbers below (built BEFORE that bug was
! found, by ruling out the SCF entirely rather than guessing at it) show
! was neither SCF noise nor a missing gradient term: `trc_scf_gradient`
! was never handed `nlc_grid_level`, so that comparison differenced an
! energy on the SCF's own (coarse) VV10 grid against an analytic
! gradient on a DIFFERENT (finer) one. Fixed there; its own floor is now
! 1.38e-7, matching the plain-GGA cases.
!
! Holding D fixed removes the SCF from the comparison ENTIRELY: with no
! reconvergence at any displaced geometry, there is no DIIS error, no
! orbital-rotation (energy-weighted-density/orthonormality) term, and no
! iteration-count sensitivity -- the only things left in dE_nlc[D]/dR are
! EXACTLY the two terms trc_vv10_grad%vv10_grad_add computes (Pulay +
! grid response), because E_nlc[D] = sum_g w_g(R) rho_D(g; R) [...] is a
! plain functional of the FIXED matrix D and the geometry R, evaluated
! the same way at every displaced point. If this comparison is NOT tight
! to truncation, the gradient formula itself is wrong, full stop -- see
! this program's own suspects list below, not a step-size or SCF excuse.
!
! SUSPECTS, if this comes back loose (the formulas trc_vv10_grad.F90's
! own header derives and cites):
!   - the weight-response coefficient dE/dw_i = beta + F (NOT beta + F/2,
!     which is what the ENERGY density eps_i carries -- a point's weight
!     enters the pair sum twice over, see vv10_grad_kernel_body's own
!     comment);
!   - the point's own explicit force fexp, in particular its sign and the
!     factor of 3 (not 1.5 -- the pair-sum symmetry that lets one pass
!     over the outer index capture both the point's own-point AND
!     every-other-point's-inner-point role at once);
!   - dft_partition_becke_derivatives's own dP/dR for the INNER grid (used
!     identically to the semilocal term's, but through a different
!     point set -- VV10's own, coarser grid).
!
! NO INTERNAL PROCEDURES under the program, same reason check_xc_grad_fd.F90
! gives (a `contains`-ed procedure corrupted a later collective Fock build
! on this exact NVHPC/OpenMPI toolchain, measured elsewhere in this tree):
! everything here is a module procedure instead.
!
module check_vv10_grad_fixedd_fd_cases
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_vv10, only: trc_vv10_rks, trc_vv10_uks
   use trc_vv10_grad, only: vv10_grad_add
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none
   private
   public :: run_case

   real(dp), parameter, public :: VV10_B = 6.0_dp   !! wB97X-V's own b, C
   real(dp), parameter, public :: VV10_C = 0.01_dp
   real(dp), parameter, public :: ABSTOL = 1.0e-10_dp
   !! Richardson gate -- see run_case's own header comment on why this is
   !! allowed to be far tighter than check_xc_grad_fd.F90's reconverged-SCF
   !! case: there is no SCF noise left to absorb here. MEASURED (CPU,
   !! gfortran): worst abs 5.50e-13 (L1), 5.08e-12 (L2) -- ~20x headroom
   !! over the looser (L2) of the two, not a loosened gate (this case did
   !! not exist before this branch).
   real(dp), parameter, public :: INVTOL = 1.0e-14_dp
   !! MEASURED: ~4e-17 at both levels -- ~250x headroom.

contains

   subroutine run_case(xyzfile, basis_json, nalpha, nbeta, nspin, grid_level, label, any_fail)
      character(*), intent(in) :: xyzfile, basis_json, label
      integer, intent(in) :: nalpha, nbeta, nspin, grid_level
      !! `nspin` = 1 (RKS, `dens(:,:,1)` the total density) or 2 (UKS,
      !! `dens(:,:,1)`/`dens(:,:,2)` the FIXED alpha/beta densities
      !! separately -- the nspin = 2 branch of `vv10_grad_add`'s own
      !! `dtot` sum, and `e_nlc_at`'s own `trc_vv10_uks` path, are
      !! otherwise never exercised by this program's RKS-only cases).
      logical, intent(inout) :: any_fail

      integer :: nat, dev, a, d, nao
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :), rtmp(:, :)
      real(dp), allocatable :: dens(:, :, :), grad_an(:, :)
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      real(dp) :: eplus, eminus, fd_h, fd_hh, fd_rich
      real(dp) :: worst_abs_h, worst_rel_h, worst_abs_rich, worst_rel_rich, denom
      real(dp) :: inv(3)
      real(dp), parameter :: H1 = 5.0e-4_dp, H2 = 2.5e-4_dp
      !! NOT the 1e-3/5e-4 pair check_xc_grad_fd.F90 uses -- measured (CPU,
      !! gfortran) on this program's own RKS L2 case: fd(H=1e-3) lands
      !! 7.37e-8 off the analytic gradient while fd(H=5e-4) and
      !! fd(H=2.5e-4) are both within ~2e-10, in a clean ~4x-per-halving
      !! ratio (textbook O(H^2) central-difference truncation, confirmed
      !! against the analytic value itself, not merely against each
      !! other). H=1e-3 at this denser grid is a genuine outlier, not the
      !! leading term of a smooth series. VERIFIED, not merely guessed:
      !! probing `vv10_density_tile`'s own survivor count (points with
      !! rho >= VV10_RHO_TOL) at the oxygen's z +/-1e-3 and +/-5e-4
      !! displacements with the SAME converged D measures 18970 at -1e-3
      !! against 18966 everywhere else (+1e-3, +/-5e-4, and undisplaced)
      !! -- a discrete active-set change specific to that one step, not
      !! present at the smaller one, consistent in sign and scale with an
      !! O(1e-8) jump in a quantity this size. Richardson from 5e-4/2.5e-4
      !! alone lands 5.08e-12 off the analytic value -- essentially the
      !! FD/double-precision floor.

      dev = trc_bind_device(0)
      call read_xyz(xyzfile, nat, at_z, r0)
      call trc_basis_from_json(basis_json, nat, at_z, r0, bas, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_vv10_grad_fixedd_fd: '//err%get_message(); stop 1
      end if
      call bas%to_device()
      nao = bas%nao

      ! The density this program holds fixed: wB97X-V's own converged D at
      ! the UNDISPLACED geometry -- `res%dmat(:,:,1:nspin)` directly
      ! (restricted holds the TOTAL density at nspin=1; unrestricted
      ! holds alpha/beta separately at nspin=2, trc_scf_driver.F90's own
      ! per-channel `guess_from` build), at the Cartesian `nao`
      ! `force_cart=.true.` put the basis in, so no spherical fold is
      ! needed either.
      opts%functional = 'wb97x-v'
      opts%grid_level = 3
      opts%nlc_grid_level = grid_level
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-10_dp; opts%max_iteration = 200
      call trc_scf_run(bas, nalpha, nbeta, opts, res)
      if (.not. res%converged) then
         print '(a)', 'check_vv10_grad_fixedd_fd: '//label//': reference SCF did not converge'
         any_fail = .true.; return
      end if

      allocate (dens(nao, nao, nspin))
      dens(:, :, 1:nspin) = res%dmat(:, :, 1:nspin)
      call lp%init(nao)

      allocate (grad_an(3, nat))
      grad_an = 0.0_dp
      !$acc enter data copyin(dens, grad_an)
      call vv10_grad_add(bas, lp, nspin, dens, VV10_B, VV10_C, grad_an, grid_level=grid_level, &
                             grid_response=.true.)
      !$acc update self(grad_an)
      !$acc exit data delete(dens, grad_an)
      call lp%release()
      call bas%release()

      inv = sum(grad_an, dim=2)

      allocate (rtmp(3, nat))
      worst_abs_h = 0.0_dp; worst_rel_h = 0.0_dp
      worst_abs_rich = 0.0_dp; worst_rel_rich = 0.0_dp
      do a = 1, nat
         do d = 1, 3
            rtmp = r0; rtmp(d, a) = r0(d, a) + H1
            eplus = e_nlc_at(basis_json, nat, at_z, rtmp, nao, nspin, dens, grid_level)
            rtmp = r0; rtmp(d, a) = r0(d, a) - H1
            eminus = e_nlc_at(basis_json, nat, at_z, rtmp, nao, nspin, dens, grid_level)
            fd_h = (eplus - eminus)/(2.0_dp*H1)

            rtmp = r0; rtmp(d, a) = r0(d, a) + H2
            eplus = e_nlc_at(basis_json, nat, at_z, rtmp, nao, nspin, dens, grid_level)
            rtmp = r0; rtmp(d, a) = r0(d, a) - H2
            eminus = e_nlc_at(basis_json, nat, at_z, rtmp, nao, nspin, dens, grid_level)
            fd_hh = (eplus - eminus)/(2.0_dp*H2)

            ! Richardson, same combination check_xc_grad_fd.F90 uses
            ! (fd_h at step H1, fd_hh at step H2 = H1/2):
            fd_rich = (4.0_dp*fd_hh - fd_h)/3.0_dp

            denom = max(abs(grad_an(d, a)), abs(fd_h), 1.0e-10_dp)
            worst_abs_h = max(worst_abs_h, abs(grad_an(d, a) - fd_h))
            worst_rel_h = max(worst_rel_h, abs(grad_an(d, a) - fd_h)/denom)
            denom = max(abs(grad_an(d, a)), abs(fd_rich), 1.0e-10_dp)
            worst_abs_rich = max(worst_abs_rich, abs(grad_an(d, a) - fd_rich))
            worst_rel_rich = max(worst_rel_rich, abs(grad_an(d, a) - fd_rich)/denom)
         end do
      end do

      print '(a,a)', '  ', label
      print '(a,es10.2,a,es10.2)', '    plain H=5e-4:      worst abs ', worst_abs_h, '  worst rel ', worst_rel_h
      print '(a,es10.2,a,es10.2)', '    Richardson(H,H/2): worst abs ', worst_abs_rich, '  worst rel ', worst_rel_rich
      print '(a,3es10.2)', '    translational-invariance residual: ', inv
      if (worst_abs_rich > ABSTOL) any_fail = .true.
      if (any(abs(inv) > INVTOL)) any_fail = .true.

      deallocate (at_z, r0, rtmp, dens, grad_an)
   end subroutine run_case

   real(dp) function e_nlc_at(basis_json, nat, at_z, r, nao, nspin, dens, grid_level) result(e)
      !! E_nlc[D] at geometry `r`, D held fixed (the caller's own `dens`,
      !! in the SAME AO basis/order every geometry shares): rebuilds the
      !! basis and the VV10 grid at `r` alone, then one `trc_vv10_rks`
      !! (nspin=1) or `trc_vv10_uks` (nspin=2) call. No SCF, no potential
      !! matrix kept (`vmat`/`vmat_a`/`vmat_b` are scratch, discarded).
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nao, nspin, grid_level
      real(dp), intent(in) :: r(3, nat)
      real(dp), intent(in) :: dens(nao, nao, nspin)
      type(trc_basis_t) :: b2
      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(error_t) :: err
      real(dp), allocatable :: vmat(:, :), vmat_b(:, :)

      call trc_basis_from_json(basis_json, nat, at_z, r, b2, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_vv10_grad_fixedd_fd: '//err%get_message(); stop 1
      end if
      if (b2%nao /= nao) then
         print '(a)', 'check_vv10_grad_fixedd_fd: nao changed across geometries'; stop 1
      end if
      call b2%to_device()

      call dft_grid_build(r, at_z, grid, err, level=grid_level)
      if (err%has_error()) then
         print '(a)', 'check_vv10_grad_fixedd_fd: '//err%get_message(); stop 1
      end if
      call xg%build(grid%n_points, grid%coords, grid%weights, b2, max_pts=512)
      call xg%to_device()

      allocate (vmat(nao, nao))
      vmat = 0.0_dp
      if (nspin == 1) then
         call trc_vv10_rks(b2, xg, dens(:, :, 1), VV10_B, VV10_C, e, vmat)
      else
         allocate (vmat_b(nao, nao))
         vmat_b = 0.0_dp
         call trc_vv10_uks(b2, xg, dens(:, :, 1), dens(:, :, 2), VV10_B, VV10_C, e, vmat, vmat_b)
         deallocate (vmat_b)
      end if
      deallocate (vmat)

      call xg%release()
      call grid%destroy()
      call b2%release()
   end function e_nlc_at

end module check_vv10_grad_fixedd_fd_cases

program check_vv10_grad_fixedd_fd
   use pic_types, only: dp
   use check_vv10_grad_fixedd_fd_cases, only: run_case
   implicit none
   logical :: any_fail
   any_fail = .false.

   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 1, 1, &
                 'RKS wB97X-V water/6-31G*, VV10 only, fixed D, L1', any_fail)
   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 1, 2, &
                 'RKS wB97X-V water/6-31G*, VV10 only, fixed D, L2', any_fail)

   ! UKS: NH2, alpha and beta densities held fixed SEPARATELY (not their
   ! sum) -- the only case in this program that exercises
   ! vv10_grad_add's nspin=2 `dtot` branch and e_nlc_at's
   ! `trc_vv10_uks` path, neither reached by the RKS cases above.
   call run_case('nh2.xyz', 'basis_sets/6-31gs.json', 5, 4, 2, 1, &
                 'UKS wB97X-V NH2/6-31G*, VV10 only, fixed D, L1', any_fail)
   call run_case('nh2.xyz', 'basis_sets/6-31gs.json', 5, 4, 2, 2, &
                 'UKS wB97X-V NH2/6-31G*, VV10 only, fixed D, L2', any_fail)

   if (any_fail) then
      print '(a)', '  RESULT: FAIL'; stop 1
   end if
   print '(a)', '  RESULT: PASS'
end program check_vv10_grad_fixedd_fd
