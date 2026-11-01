!
! The analytic PCM (C-PCM and IEF-PCM) nuclear gradient (trc_pcm_grad's
! pcm_grad_add: grad_qv + grad_nuc + grad_solver, added on top of
! trc_scf_gradient's own vacuum RHF gradient via pcm_grad_scf_add) against
! pyscf's own solvent-wrapped `Gradients` class -- water/6-31G*,
! eps = 78.3553, pyscf's own SWIG defaults (lebedev_order 29 -> 302 points
! per sphere, vdw_scale 1.2), matching terco's.
!
! FIXED INITIAL DENSITY: NOT NEEDED (test/pyscf_refs/generate.py's
! gen_pcm_grad: a closed-shell ground state has no nearby competing
! stationary point, continuum or not), so terco and pyscf converge
! independently and are compared at their own, separately converged
! densities -- the same reasoning check_full_grad_vs_pyscf.F90's own
! water/formamide cases use.
!
! check_pcm_grad_fd.F90 already validates this feature against a central
! finite difference of terco's OWN energy (pyscf-independent); this
! program is the external cross-check against a different code's
! analytic PCM gradient.
!
module check_pcm_grad_vs_pyscf_suite
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
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_str, ref_read_matrix
   implicit none
   private
   public :: collect_pcm_grad_vs_pyscf_tests, load_reference

   real(dp), parameter :: ABSTOL = 1.0e-6_dp
   !! Analytic-vs-analytic at two independently converged (not bit-
   !! identical) SCFs -- check_full_grad_vs_pyscf.F90's own ABSTOL and
   !! reasoning: not expected tighter than check_pcm_grad_fd.F90's
   !! FD-truncation-limited gate, the same order of magnitude on top of
   !! the two SCFs' own slightly different convergence paths.

   type :: ref_case_t
      character(len=16) :: model = ''
      real(dp) :: e = 0.0_dp, e_diel = 0.0_dp, gnorm = 0.0_dp
      real(dp), allocatable :: grad(:, :)
   end type ref_case_t

   integer :: g_nao = 0
   type(ref_case_t) :: g_cases(2)
      !! Filled by load_reference, in the file's own write order: C-PCM,
      !! then IEF-PCM (generate.py's gen_pcm_grad).

contains

   subroutine load_reference()
      integer :: u
      call ref_require('pcm_grad.ref', u)
      call ref_read_int(u, 'NAO', g_nao)
      call read_one_case(u, 'cpcm', g_cases(1))
      call read_one_case(u, 'iefpcm', g_cases(2))
      call ref_close(u)
   end subroutine load_reference

   subroutine read_one_case(u, name, c)
      integer, intent(in) :: u
      character(*), intent(in) :: name
      type(ref_case_t), intent(out) :: c
      character(len=16) :: model_ref
      call ref_read_str(u, 'MODEL_'//trim(name), model_ref)
      c%model = model_ref
      call ref_read_real(u, 'E_'//trim(name), c%e)
      call ref_read_real(u, 'E_DIEL_'//trim(name), c%e_diel)
      call ref_read_real(u, 'GRAD_NORM_'//trim(name), c%gnorm)
      allocate (c%grad(3, 3))   ! water: 3 atoms
      call ref_read_matrix(u, 'GRAD_'//trim(name), 3, 3, c%grad)
   end subroutine read_one_case

   subroutine collect_pcm_grad_vs_pyscf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water/6-31G* C-PCM gradient vs pyscf", test_cpcm), &
                  new_unittest("water/6-31G* IEF-PCM gradient vs pyscf", test_iefpcm) &
                  ]
   end subroutine collect_pcm_grad_vs_pyscf_tests

   subroutine test_cpcm(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, g_cases(1))
   end subroutine test_cpcm

   subroutine test_iefpcm(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, g_cases(2))
   end subroutine test_iefpcm

   subroutine run_case(error, c)
      type(error_type), allocatable, intent(out) :: error
      type(ref_case_t), intent(in) :: c

      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(trc_pcm_t) :: pcm
      type(error_t) :: err
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :), grad(:, :)
      integer :: nat, dev
      real(dp) :: worst
      character(len=16) :: model

      model = trim(c%model)
      dev = trc_bind_device(0)
      call read_xyz('water.xyz', nat, at_z, r0)
      call trc_basis_from_json('basis_sets/6-31gs.json', nat, at_z, r0, bas, err)
      call check(error,.not. err%has_error(), trim(model)//': basis build failed: '//err%get_message())
      if (allocated(error)) return
      call check(error, bas%nao == g_nao, trim(model)//': reference nao does not match this build')
      if (allocated(error)) then
         call bas%release(); return
      end if
      call bas%to_device()

      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iteration = 200
      opts%pcm_model = model; opts%pcm_eps = 78.3553_dp
      call trc_scf_run(bas, 5, 5, opts, res)
      call check(error, res%converged, trim(model)//': terco SCF did not converge')
      if (allocated(error)) then
         call bas%release(); return
      end if

      allocate (grad(3, nat))
      call lp%init(bas%nao)
      call trc_scf_gradient(bas, lp, res%nspin, 5, 5, res%cmo, res%eps, 1.0e-12_dp, grad)
      call pcm%init(pcm_model_from_name(model), opts%pcm_eps, err)
      call check(error,.not. err%has_error(), trim(model)//': pcm%init failed: '//err%get_message())
      if (allocated(error)) then
         call lp%release(); call bas%release(); return
      end if
      call pcm%build(lp, bas, opts%pcm_lebedev, opts%pcm_radii_scale, err)
      call check(error,.not. err%has_error(), trim(model)//': pcm%build failed: '//err%get_message())
      if (allocated(error)) then
         call pcm%release(); call lp%release(); call bas%release(); return
      end if
      call pcm_grad_scf_add(bas, lp, [5, 5], res%cmo, res%eps, 1.0e-12_dp, pcm, grad)
      call pcm%release()
      call lp%release()
      call bas%release()

      worst = maxval(abs(grad - c%grad))
      print '(a,a,a,es12.4,a,es12.4,a,es10.2)', '  ', trim(model), ': |E(terco)-E(pyscf)| ', &
         abs(res%energy - c%e), '  pyscf |g| ', c%gnorm, '  worst|diff| ', worst
      call check(error, worst < ABSTOL, trim(model)//': gradient vs pyscf, worst |diff| = '// &
                 to_string(worst, "es10.2"))

      deallocate (at_z, r0, grad)
   end subroutine run_case

end module check_pcm_grad_vs_pyscf_suite

program check_pcm_grad_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_pcm_grad_vs_pyscf_suite, only: collect_pcm_grad_vs_pyscf_tests, load_reference
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat

   call load_reference()

   stat = 0
   testsuites = [new_testsuite("check_pcm_grad_vs_pyscf", collect_pcm_grad_vs_pyscf_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_pcm_grad_vs_pyscf
