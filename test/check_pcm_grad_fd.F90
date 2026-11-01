!
! The analytic PCM (C-PCM / IEF-PCM) nuclear gradient against a central
! finite difference of terco's own FULLY RECONVERGED SCF energy (same
! technique as check_full_grad_fd.F90, which this file otherwise mirrors),
! plus translational invariance of the analytic gradient alone.
!
! `trc_scf_gradient`'s vacuum terms (eri%gradient, trc_1e_gradient, nuclear
! repulsion) are ALREADY validated against pyscf and FD elsewhere
! (check_full_grad_fd.F90, check_full_grad_vs_pyscf.F90); this program is
! the FIRST thing that fails if the PCM piece this feature adds on top
! (trc_pcm_grad%pcm_gradient_add: grad_qv + grad_nuc + grad_solver) is
! wrong, since `trc_scf_run` already includes E_diel in the reconverged
! energy FD differentiates.
!
! Cases: water (RHF, C-PCM and IEF-PCM) and formamide (RHF, IEF-PCM) --
! formamide's four atoms sit close enough together that several surface
! points are in another atom's SWIG switching region, which the plain
! water case barely exercises (its three spheres are closer to isolated).
! eps = 78.3553 (water), 6-31G* throughout.
!
! NO INTERNAL PROCEDURES -- same reason check_full_grad_fd.F90's own
! header gives (an internal procedure anywhere in the executable corrupted
! a later collective Fock build on this toolchain, measured there): the
! gradient and FD helpers below are plain module procedures.
!
module check_pcm_grad_fd_suite
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_linalg, only: trc_linalg_t
   use trc_pcm, only: trc_pcm_t, pcm_model_from_name
   use trc_pcm_grad, only: pcm_grad_scf_add
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_pcm_grad_fd_tests

   real(dp), parameter :: H = 1.0e-4_dp
   !! The task's own step; a reconverged-SCF FD, so the floor is truncation
   !! (step^2 terms), not SCF noise at these tolerances.
   real(dp), parameter :: RELTOL = 1.0e-4_dp
   real(dp), parameter :: ABSTOL = 1.0e-6_dp
   real(dp), parameter :: INVTOL = 1.0e-9_dp

contains

   subroutine collect_pcm_grad_fd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF water/6-31G* C-PCM", test_water_cpcm), &
                  new_unittest("RHF water/6-31G* IEF-PCM", test_water_iefpcm), &
                  new_unittest("RHF formamide/6-31G* IEF-PCM", test_formamide_iefpcm) &
                  ]
   end subroutine collect_pcm_grad_fd_tests

   subroutine test_water_cpcm(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'cpcm', 'RHF water/6-31G* C-PCM')
   end subroutine test_water_cpcm

   subroutine test_water_iefpcm(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'iefpcm', 'RHF water/6-31G* IEF-PCM')
   end subroutine test_water_iefpcm

   subroutine test_formamide_iefpcm(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'formamide.xyz', 'basis_sets/6-31gs.json', 12, 12, 'iefpcm', &
                    'RHF formamide/6-31G* IEF-PCM')
   end subroutine test_formamide_iefpcm

   subroutine run_case(error, xyzfile, basis_json, nalpha, nbeta, model, label)
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: xyzfile, basis_json, model, label
      integer, intent(in) :: nalpha, nbeta

      integer :: nat, dev, a, d
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :), rtmp(:, :), grad(:, :)
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(trc_pcm_t) :: pcm
      type(error_t) :: err
      real(dp) :: eplus, eminus, fd, worst_abs, worst_rel, denom
      integer :: nnan
      real(dp) :: inv(3)

      dev = trc_bind_device(0)
      call read_xyz(xyzfile, nat, at_z, r0)
      call trc_basis_from_json(basis_json, nat, at_z, r0, bas, err)
      call check(error,.not. err%has_error(), label//': basis build failed: '//err%get_message())
      if (allocated(error)) return
      call bas%to_device()

      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iteration = 200
      opts%pcm_model = model; opts%pcm_eps = 78.3553_dp
      call trc_scf_run(bas, nalpha, nbeta, opts, res)
      call check(error, res%converged, label//': reference SCF did not converge')
      if (allocated(error)) then
         call bas%release()
         return
      end if

      allocate (grad(3, nat))
      call lp%init(bas%nao)
      call pcm%init(pcm_model_from_name(model), opts%pcm_eps, err)
      call check(error,.not. err%has_error(), label//': pcm%init failed: '//err%get_message())
      if (allocated(error)) then
         call lp%release(); call bas%release(); return
      end if
      call pcm%build(lp, bas, opts%pcm_lebedev, opts%pcm_radii_scale, err)
      call check(error,.not. err%has_error(), label//': pcm%build failed: '//err%get_message())
      if (allocated(error)) then
         call pcm%release(); call lp%release(); call bas%release(); return
      end if
      call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, 1.0e-12_dp, grad)
      call pcm_grad_scf_add(bas, lp, [nalpha, nbeta], res%cmo, res%eps, 1.0e-12_dp, pcm, grad)
      call pcm%release()
      ! Released before `lp` goes out of scope -- trc_linalg_t has no
      ! finalizer (check_full_grad_fd.F90's own comment).
      call lp%release()
      call bas%release()

      inv = sum(grad, dim=2)
      call check(error, all(abs(inv) <= INVTOL), label// &
                 ': translational invariance residual = '//to_string(maxval(abs(inv)), "es10.2"))
      if (allocated(error)) return

      allocate (rtmp(3, nat))
      worst_abs = 0.0_dp; worst_rel = 0.0_dp; nnan = 0
      do a = 1, nat
         do d = 1, 3
            rtmp = r0; rtmp(d, a) = r0(d, a) + H
            eplus = scf_energy_at(basis_json, nat, at_z, rtmp, nalpha, nbeta, opts)
            rtmp = r0; rtmp(d, a) = r0(d, a) - H
            eminus = scf_energy_at(basis_json, nat, at_z, rtmp, nalpha, nbeta, opts)
            fd = (eplus - eminus)/(2.0_dp*H)
            denom = max(abs(grad(d, a)), abs(fd), 1.0e-10_dp)
            worst_abs = max(worst_abs, abs(grad(d, a) - fd))
            worst_rel = max(worst_rel, abs(grad(d, a) - fd)/denom)
            if (ieee_is_nan(grad(d, a) - fd)) nnan = nnan + 1
         end do
      end do

      print '(a,a,a,es10.2,a,es10.2,a,3es10.2)', '  ', label, ': worst abs ', worst_abs, &
         '  worst rel ', worst_rel, '  translational-invariance residual: ', inv
      call check(error, nnan == 0 .and. .not. (worst_abs > ABSTOL .and. worst_rel > RELTOL), label// &
                 ': analytic vs FD gradient, worst abs = '//to_string(worst_abs, "es10.2")// &
                 ', worst rel = '//to_string(worst_rel, "es10.2")//', NaN components '//to_string(nnan))

      deallocate (at_z, r0, rtmp, grad)
   end subroutine run_case

   real(dp) function scf_energy_at(basis_json, nat, at_z, r, nalpha, nbeta, opts) result(e)
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nalpha, nbeta
      real(dp), intent(in) :: r(3, nat)
      type(trc_scf_options_t), intent(in) :: opts
      type(trc_basis_t) :: b2
      type(trc_scf_result_t) :: r2
      type(error_t) :: err
      call trc_basis_from_json(basis_json, nat, at_z, r, b2, err)
      if (err%has_error()) then
         print '(a)', 'check_pcm_grad_fd: '//err%get_message(); stop 1
      end if
      call b2%to_device()
      call trc_scf_run(b2, nalpha, nbeta, opts, r2)
      e = r2%energy
      call b2%release()
   end function scf_energy_at

end module check_pcm_grad_fd_suite

program check_pcm_grad_fd
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_pcm_grad_fd_suite, only: collect_pcm_grad_fd_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat

   stat = 0
   testsuites = [new_testsuite("check_pcm_grad_fd", collect_pcm_grad_fd_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_pcm_grad_fd
