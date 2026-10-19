!
! The density-fitted (`--df`) analytic gradient (trc_scf_grad%trc_scf_gradient
! with `dfjk`/`aux` present: trc_scf_grad_df's fitted J and K + trc_1e_gradient
! + the XC term + nuclear repulsion) against a Richardson-extrapolated
! central finite difference of terco's own FULLY RECONVERGED --df energy:
! every displaced geometry rebuilds the orbital and jkfit bases, refits
! (metric, Cholesky, (mu nu|P), B) and reconverges the SCF from scratch.
!
!   g_FD = (4 D(h/2) - D(h))/3,  D(h) = (E(+h) - E(-h))/(2h)
!
! cancels the O(h^2) truncation term of the plain central difference, so
! what is left is the SCF's own convergence noise -- the same scheme as
! metalquicha's check_df_ref_gradient.f90.
!
! RHF (fitted J and K), RKS PBE (fitted J only: K is skipped outright for a
! pure functional) and RKS B3LYP (fitted J and K scaled by 0.2), water,
! 6-31G*, def2-universal-jkfit (g auxiliary shells). Also checks
! translational invariance of the analytic gradient alone (the XC term's
! grid response is on, trc_scf_gradient's default, so the DFT cases are as
! invariant as the RHF one).
!
! GATES, from the measured agreement (gfortran 15, CPU): |analytic - FD|
! 1.4e-10 (RHF), 8.3e-11 (PBE), 8.4e-11 (B3LYP), at step 2e-3 bohr with the
! SCF converged to 1e-12 Eh; invariance 1.8e-13 to 1.9e-13. The FD gate is
! 1e-8 (the SCF's own 1e-12 Eh convergence over a 1e-3 bohr step is ~1e-9
! of noise), the invariance gate 1e-10.
!
! `check_df_grad_fd formamide` runs RHF on formamide/cc-pVDZ instead (its own
! ctest entry), the case whose pyscf comparison is limited by the fit's
! conditioning (check_df_grad_vs_pyscf.F90's `tol`): this is the
! pyscf-independent check of the same gradient. Measured 7.1e-10 against the
! same 1e-8 gate, invariance 1.7e-12.
!
! NO INTERNAL PROCEDURES (see check_full_grad_fd.F90's header): the cases
! are module procedures.
!
module check_df_grad_fd_cases
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_scf_df, only: df_jk_t
   use trc_df_kernels, only: DF_LMAX_AUX
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none
   private
   public :: run_case

   real(dp), parameter :: H = 2.0e-3_dp
   real(dp), parameter :: THRESH = 1.0e-12_dp
   character(len=*), parameter :: AUX_JSON = 'basis_sets/def2-universal-jkfit.json'

contains

   subroutine run_case(xyzfile, basis_json, nalpha, nbeta, functional, label, abstol, invtol, any_fail)
      character(*), intent(in) :: xyzfile, basis_json, functional, label
      integer, intent(in) :: nalpha, nbeta
      real(dp), intent(in) :: abstol, invtol
      logical, intent(inout) :: any_fail

      integer :: nat, dev, a, d, k
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :), rtmp(:, :), grad(:, :)
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(df_jk_t) :: dfjk
      type(trc_xc_functional_t) :: func

      real(dp) :: e(-2:2), d1, d2, fd, worst_abs, inv(3), hh
      logical :: ok

      dev = trc_bind_device(0)
      call read_xyz(xyzfile, nat, at_z, r0)
      opts%df = .true.
      opts%functional = functional
      opts%grid_level = 3
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-10_dp; opts%max_iter = 200
      if (len_trim(functional) > 0) then
         ok = xc_functional_by_name(functional, func)
         if (.not. ok) then
            print '(a)', 'check_df_grad_fd: unknown functional '//functional; stop 1
         end if
      end if

      call setup(basis_json, nat, at_z, r0, nalpha, bas, aux, pl, lp, dfjk)
      if (maxval(aux%sh_l) > DF_LMAX_AUX) then
         print '(a)', '  SKIP: def2-universal-jkfit needs g auxiliary shells, above this build''s DF_LMAX_AUX.'
         stop 0
      end if
      call trc_scf_run(bas, nalpha, nbeta, opts, res, la=lp, dfjk=dfjk)
      if (.not. res%converged) then
         print '(a)', 'check_df_grad_fd: '//label//': reference SCF did not converge'
         any_fail = .true.; return
      end if
      allocate (grad(3, nat))
      if (len_trim(functional) > 0) then
         call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad, &
                               func=func, grid_level=opts%grid_level, grid_max_pts=opts%grid_max_pts, &
                               rho_tol=opts%rho_tol, dfjk=dfjk, aux=aux)
      else
         call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad, &
                               dfjk=dfjk, aux=aux)
      end if
      call teardown(bas, aux, pl, lp, dfjk)

      inv = sum(grad, dim=2)
      if (maxval(abs(inv)) > invtol) any_fail = .true.

      allocate (rtmp(3, nat))
      worst_abs = 0.0_dp
      do a = 1, nat
         do d = 1, 3
            do k = -2, 2
               if (k == 0) cycle
               hh = real(k, dp)*0.5_dp*H
               rtmp = r0; rtmp(d, a) = r0(d, a) + hh
               e(k) = energy_at(basis_json, nat, at_z, rtmp, nalpha, nbeta, opts)
            end do
            d1 = (e(2) - e(-2))/(2.0_dp*H)
            d2 = (e(1) - e(-1))/H
            fd = (4.0_dp*d2 - d1)/3.0_dp
            worst_abs = max(worst_abs, abs(grad(d, a) - fd))
            print '(a,i0,a,i0,a,f18.12,a,f18.12,a,es9.2)', '    atom ', a, ' xyz ', d, ' analytic ', grad(d, a), &
               '  FD ', fd, '  diff ', grad(d, a) - fd
         end do
      end do

      print '(a,a,a,es10.2,a,es10.2,a,3es10.2)', '  ', label, ': worst |analytic - Richardson FD| ', worst_abs, &
         ' (gate ', abstol, ')  translational-invariance residual: ', inv
      if (worst_abs > abstol) any_fail = .true.
      deallocate (at_z, r0, rtmp, grad)
   end subroutine run_case

   subroutine setup(basis_json, nat, at_z, r, nalpha, bas, aux, pl, lp, dfjk)
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nalpha
      real(dp), intent(in) :: r(3, nat)
      type(trc_basis_t), intent(inout) :: bas, aux
      type(trc_pairlist_t), intent(inout) :: pl
      type(trc_linalg_t), intent(inout) :: lp
      type(df_jk_t), intent(inout) :: dfjk
      type(error_t) :: err
      call trc_basis_from_json(basis_json, nat, at_z, r, bas, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_df_grad_fd: '//err%get_message(); stop 1
      end if
      call trc_basis_from_json(AUX_JSON, nat, at_z, r, aux, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_df_grad_fd: '//err%get_message(); stop 1
      end if
      if (maxval(aux%sh_l) > DF_LMAX_AUX) return
      call bas%to_device()
      call aux%to_device()
      call pl%build(bas, THRESH)
      call pl%to_device()
      call lp%init(bas%nao)
      call dfjk%setup(bas, aux, pl, lp, err, nocc_hint=max(nalpha, 1))
      if (err%has_error()) then
         print '(a)', 'check_df_grad_fd: '//err%get_message(); stop 1
      end if
      call pl%release()
   end subroutine setup

   subroutine teardown(bas, aux, pl, lp, dfjk)
      !! `lp` released before it goes out of scope: trc_linalg_t has no
      !! finalizer (see check_full_grad_fd.F90).
      type(trc_basis_t), intent(inout) :: bas, aux
      type(trc_pairlist_t), intent(inout) :: pl
      type(trc_linalg_t), intent(inout) :: lp
      type(df_jk_t), intent(inout) :: dfjk
      call dfjk%release()
      call pl%release()
      call lp%release()
      call aux%release()
      call bas%release()
   end subroutine teardown

   real(dp) function energy_at(basis_json, nat, at_z, r, nalpha, nbeta, opts) result(e)
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nalpha, nbeta
      real(dp), intent(in) :: r(3, nat)
      type(trc_scf_options_t), intent(in) :: opts
      type(trc_basis_t) :: b2, a2
      type(trc_pairlist_t) :: p2
      type(trc_linalg_t) :: l2
      type(df_jk_t) :: f2
      type(trc_scf_result_t) :: r2
      call setup(basis_json, nat, at_z, r, nalpha, b2, a2, p2, l2, f2)
      call trc_scf_run(b2, nalpha, nbeta, opts, r2, la=l2, dfjk=f2)
      if (.not. r2%converged) then
         print '(a)', 'check_df_grad_fd: a displaced SCF did not converge'; stop 1
      end if
      e = r2%energy
      call teardown(b2, a2, p2, l2, f2)
   end function energy_at

end module check_df_grad_fd_cases

program check_df_grad_fd
   use pic_types, only: dp
   use check_df_grad_fd_cases, only: run_case
   implicit none
   logical :: any_fail
   character(len=32) :: which
   any_fail = .false.
   which = 'water'
   if (command_argument_count() >= 1) call get_command_argument(1, which)

   if (trim(which) == 'formamide') then
      call run_case('formamide.xyz', 'basis_sets/cc-pvdz.json', 12, 12, '', 'RHF --df formamide/cc-pVDZ', &
                    1.0e-8_dp, 1.0e-10_dp, any_fail)
      if (any_fail) then
         print '(a)', '  RESULT: FAIL'; stop 1
      end if
      print '(a)', '  RESULT: PASS'
      stop
   end if

   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, '', 'RHF --df water/6-31G*', &
                 1.0e-8_dp, 1.0e-10_dp, any_fail)
   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'pbe', 'RKS PBE --df water/6-31G*', &
                 1.0e-8_dp, 1.0e-10_dp, any_fail)
   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'b3lyp', 'RKS B3LYP --df water/6-31G*', &
                 1.0e-8_dp, 1.0e-10_dp, any_fail)

   if (any_fail) then
      print '(a)', '  RESULT: FAIL'; stop 1
   end if
   print '(a)', '  RESULT: PASS'
end program check_df_grad_fd
