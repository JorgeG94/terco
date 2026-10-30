!
! The FULL analytic RHF/UHF nuclear gradient (trc_scf_grad%trc_scf_gradient:
! one-electron + eri%gradient's two-electron part + nuclear repulsion)
! against pyscf's own `Gradients` class AT A CONVERGED density, on
! SPHERICAL cc-pVDZ bases -- water/RHF, formamide/RHF, OH/UHF.
!
! WHY THIS PROGRAM, ON TOP OF check_full_grad_fd.F90 AND
! check_1e_grad_vs_pyscf.F90. Neither of those exercises trc_scf_gradient's
! own spherical fold (D and W folded nao_sph -> the Cartesian nao the
! kernels need, T D T^T -- src/scf/trc_scf_grad.F90's own header):
! check_full_grad_fd.F90 only runs 6-31G* (Cartesian, no fold at all), and
! check_1e_grad_vs_pyscf.F90 hands trc_1e_gradient a synthetic density
! directly at the Cartesian nao, bypassing the SCF (and the fold) entirely
! even on a basis BSE tags spherical. This program is the one path that
! builds a converged SPHERICAL SCF, then asks trc_scf_gradient for its
! gradient -- the fold runs for real here.
!
! FIXED INITIAL DENSITY. See test/pyscf_refs/generate.py's gen_full_grad
! for the full reasoning (its own docstring): water and formamide need
! none (a closed-shell ground state has no nearby competing stationary
! point); OH does, the same recipe rks.ref's own OH cases use, because its
! unpaired electron sits in a genuinely degenerate pair of pi orbitals.
! The committed DM0 (OH only) seeds `trc_scf_run`'s own `dguess`, so both
! sides converge toward the SAME orientation.
!
! The three cases share one reference file (full_grad_ccpvdz.ref), read
! strictly in write order (pyscf_ref_io.F90's own header: tags are read
! sequentially, not by seeking a name, and generate.py's gen_full_grad
! writes water, formamide, then OH). All three are read once, up front
! (load_reference, called from the program below before run_testsuite),
! into g_cases; each test-drive case then just compares against its own
! already-read record.
!
module check_full_grad_vs_pyscf_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix, ref_read_matrix3
   implicit none
   private
   public :: collect_full_grad_vs_pyscf_tests, load_reference

   real(dp), parameter :: ABSTOL = 1.0e-6_dp
   !! Analytic-vs-analytic at two independently converged (but, for water/
   !! formamide, not bit-identical) SCFs -- tighter than check_full_grad_fd's
   !! FD-truncation-limited 5e-7 is not expected; this is the same order of
   !! magnitude, on top of whatever residual the two SCFs' own, slightly
   !! different convergence paths leave in the density.

   type :: ref_case_t
      character(len=64) :: name = ''
      integer :: nspin = 0, nao = 0, nat = 0
      real(dp) :: e = 0.0_dp, gnorm = 0.0_dp
      real(dp), allocatable :: grad(:, :)
      real(dp), allocatable :: dm0(:, :, :)
   end type ref_case_t

   type(ref_case_t) :: g_cases(3)
      !! Filled by load_reference, in the file's own write order: water,
      !! formamide, then OH (generate.py's gen_full_grad).

contains

   subroutine load_reference()
      integer :: u
      call ref_require('full_grad_ccpvdz.ref', u)
      call read_one_case(u, 'water_ccpvdz', g_cases(1))
      call read_one_case(u, 'formamide_ccpvdz', g_cases(2))
      call read_one_case(u, 'oh_ccpvdz', g_cases(3))
      call ref_close(u)
   end subroutine load_reference

   subroutine read_one_case(u, name, c)
      integer, intent(in) :: u
      character(*), intent(in) :: name
      type(ref_case_t), intent(out) :: c
      character(len=64) :: name_ref
      call ref_read_str(u, 'NAME_'//trim(name), name_ref)
      c%name = name_ref
      call ref_read_int(u, 'NSPIN_'//trim(name), c%nspin)
      call ref_read_int(u, 'NAO_'//trim(name), c%nao)
      call ref_read_int(u, 'NAT_'//trim(name), c%nat)
      if (c%nspin == 2) then
         allocate (c%dm0(c%nao, c%nao, 2))
         call ref_read_matrix3(u, 'DM0_'//trim(name), c%nao, c%nao, 2, c%dm0)
      end if
      call ref_read_real(u, 'E_'//trim(name), c%e)
      call ref_read_real(u, 'GRAD_NORM_'//trim(name), c%gnorm)
      allocate (c%grad(3, c%nat))
      call ref_read_matrix(u, 'GRAD_'//trim(name), 3, c%nat, c%grad)
   end subroutine read_one_case

   subroutine collect_full_grad_vs_pyscf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water/cc-pVDZ RHF gradient vs pyscf", test_water), &
                  new_unittest("formamide/cc-pVDZ RHF gradient vs pyscf", test_formamide), &
                  new_unittest("OH/cc-pVDZ UHF gradient vs pyscf", test_oh) &
                  ]
   end subroutine collect_full_grad_vs_pyscf_tests

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, g_cases(1), 'water.xyz', 5, 5)
   end subroutine test_water

   subroutine test_formamide(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, g_cases(2), 'formamide.xyz', 12, 12)
   end subroutine test_formamide

   subroutine test_oh(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, g_cases(3), 'oh.xyz', 5, 4)
   end subroutine test_oh

   subroutine run_case(error, c, xyzfile, nalpha, nbeta)
      type(error_type), allocatable, intent(out) :: error
      type(ref_case_t), intent(in) :: c
      character(*), intent(in) :: xyzfile
      integer, intent(in) :: nalpha, nbeta

      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :), grad(:, :)
      integer :: nat, nao, dev
      real(dp) :: worst
      character(len=64) :: name

      name = trim(c%name)
      dev = trc_bind_device(0)
      call read_xyz(xyzfile, nat, at_z, r0)
      call trc_basis_from_json('basis_sets/cc-pvdz.json', nat, at_z, r0, bas, err)
      call check(error,.not. err%has_error(), trim(name)//': basis build failed: '//err%get_message())
      if (allocated(error)) return
      call check(error, bas%spherical, trim(name)// &
                 ': cc-pVDZ did not come back spherical -- nothing to fold here')
      if (allocated(error)) then
         call bas%release()
         return
      end if
      call bas%to_device()
      nao = bas%nao_sph

      call check(error, c%nspin == merge(2, 1, nalpha /= nbeta) .and. c%nao == nao .and. c%nat == nat, &
                 trim(name)//': reference shape does not match this build (nspin/nao/nat)')
      if (allocated(error)) then
         call bas%release()
         return
      end if

      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iteration = 200
      if (c%nspin == 2) then
         call trc_scf_run(bas, nalpha, nbeta, opts, res, dguess=c%dm0)
      else
         call trc_scf_run(bas, nalpha, nbeta, opts, res)
      end if

      call check(error, res%converged, trim(name)//': terco SCF did not converge')
      if (allocated(error)) then
         call bas%release()
         return
      end if

      allocate (grad(3, nat))
      call lp%init(nao)
      call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, 1.0e-12_dp, grad)
      ! Released before `lp` goes out of scope -- trc_linalg_t has no
      ! finalizer (see check_full_grad_fd.F90's own comment on this).
      call lp%release()
      call bas%release()

      worst = maxval(abs(grad - c%grad))
      print '(a,a,a,es12.4,a,es12.4,a,es10.2)', '  ', trim(name), ': |E(terco)-E(pyscf)| ', &
         abs(res%energy - c%e), '  pyscf |g| ', c%gnorm, '  worst|diff| ', worst
      call check(error, worst < ABSTOL, trim(name)//': gradient vs pyscf, worst |diff| = '// &
                 to_string(worst, "es10.2"))

      deallocate (at_z, r0, grad)
   end subroutine run_case

end module check_full_grad_vs_pyscf_suite

program check_full_grad_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_full_grad_vs_pyscf_suite, only: collect_full_grad_vs_pyscf_tests, load_reference
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat

   call load_reference()

   stat = 0
   testsuites = [new_testsuite("check_full_grad_vs_pyscf", collect_full_grad_vs_pyscf_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_full_grad_vs_pyscf
