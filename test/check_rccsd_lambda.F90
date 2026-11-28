module check_rccsd_lambda_suite
   !!
   !! The RI-CCSD Lambda solver (src/post_hf/ccsd/trc_rccsd_lambda.F90) against pyscf's DF-CCSD Lambda of the
   !! gradient references (test/pyscf_refs/ccsd_grad_water_*.ref, generate.py's gen_ccsd_grad): water distorted to
   !! C1, STO-3G (Cartesian throughout) and cc-pVDZ (spherical), cc-pVDZ-RIFIT in the orbital basis's convention,
   !! frozen core 0 and 1. The probes do not depend on orbital signs: ||l1||, ||l2||, sum l1 t1 and sum l2 t2.
   !! The CCSD is converged as the reference's (1e-13 / 1e-11) and Lambda to a step of 1e-11.
   !!
   !!   references    the four cases
   !!   DIIS          history on the device, forced onto the host, sent there because it does not fit (a cap on the
   !!                 free memory), and no DIIS: the same multipliers
   !!   non-converged an iteration limit and a NaN are reported, never as converged
   !!   refusals      each by its message, with nothing left allocated
   !!
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_linalg, only: trc_linalg_t
   use trc_rccsd_driver, only: trc_rccsd_options_t, trc_rccsd_state_t, trc_rccsd_solve
   use trc_rccsd_lambda, only: trc_rccsd_lambda_solve_options_t, trc_rccsd_lambda_state_t, trc_rccsd_lambda_solve
   use rccsd_test_support, only: run_scf
   use ccsd_grad_support, only: NCASE, CASE_REF, CASE_BASIS, CASE_NFROZEN, ccsd_grad_ref_t, read_ccsd_grad_ref, &
                                geometry, supported
   implicit none
   private
   public :: collect_check_rccsd_lambda, release_fixture

   integer, parameter :: NPROBE = 4
   character(len=8), parameter :: PROBE(NPROBE) = [character(len=8) :: "L1_NORM", "L2_NORM", "L1_T1", "L2_T2"]
   ! absolute differences; worst measured over the four cases 2.5e-12
   real(dp), parameter :: TOL_PROBE = 3.0e-11_dp
   real(dp), parameter :: TOL_SAME = 2.0e-10_dp  !! max |l| difference between two ways of converging; worst 1.2e-11
   integer, parameter :: DIIS_CASE = 4          !! cc-pVDZ, frozen core 1

   type :: fixture_t
      logical :: built = .false.
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_rccsd_state_t) :: cc
      type(ccsd_grad_ref_t) :: ref
   end type fixture_t

   logical :: bound = .false.
   type(trc_linalg_t) :: la
   type(fixture_t) :: fix(NCASE)

contains

   subroutine collect_check_rccsd_lambda(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("STO-3G, frozen core 0 and 1: the references", test_sto3g), &
                   new_unittest("cc-pVDZ, frozen core 0 and 1: the references", test_ccpvdz), &
                   new_unittest("DIIS on the device, on the host, and off", test_diis), &
                   new_unittest("non-convergence is reported", test_nonconvergence), &
                   new_unittest("refusals", test_refusals)]
   end subroutine collect_check_rccsd_lambda

   subroutine test_sto3g(error)
      type(error_type), allocatable, intent(out) :: error
      call reference_case(error, 1)
      if (allocated(error)) return
      call reference_case(error, 2)
   end subroutine test_sto3g

   subroutine test_ccpvdz(error)
      type(error_type), allocatable, intent(out) :: error
      call reference_case(error, 3)
      if (allocated(error)) return
      call reference_case(error, 4)
   end subroutine test_ccpvdz

   function tight() result(opts)
      type(trc_rccsd_lambda_solve_options_t) :: opts

      opts%conv_amplitude = 1.0e-11_dp
      opts%conv_norm = 1.0e-10_dp
   end function tight

   subroutine reference_case(error, ic)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ic

      type(trc_rccsd_lambda_state_t) :: lam
      real(dp) :: got(NPROBE), diff(NPROBE)
      integer :: k

      if (.not. supported()) return
      call ensure_fixture(ic, error)
      if (allocated(error)) return
      call trc_rccsd_lambda_solve(fix(ic)%cc, tight(), la, lam)
      call check(error, lam%result%converged .and. lam%result%finite .and. .not. lam%result%refused, &
                 trim(CASE_REF(ic))//": "//trim(lam%result%message))
      if (allocated(error)) return
      call check(error, lam%result%step_max < 1.0e-11_dp .and. lam%result%step_norm2 < 1.0e-10_dp, &
                 trim(CASE_REF(ic))//": converged with a step above the thresholds")
      if (allocated(error)) return
      associate (l => lam%l, t => fix(ic)%cc%amps)
         got = [norm2(l%t1), norm2(l%t2), sum(l%t1*t%t1), sum(l%t2*t%t2)]
      end associate
      diff = abs(got - fix(ic)%ref%probes(1:NPROBE, 1))
      print "(a,a,i3,a,4es10.2)", trim(CASE_REF(ic)), ": iterations", lam%result%iterations, &
         "  |terco - pyscf| of L1_NORM L2_NORM L1_T1 L2_T2 ", diff
      do k = 1, NPROBE
         ! written so that a NaN fails
         call check(error, diff(k) < TOL_PROBE, trim(CASE_REF(ic))//" "//trim(PROBE(k))//": "// &
                    to_string(got(k), "es22.14")//" against "//to_string(fix(ic)%ref%probes(k, 1), "es22.14"))
         if (allocated(error)) exit
      end do
      call lam%release()
   end subroutine reference_case

   subroutine test_diis(error)
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_lambda_solve_options_t) :: opts
      type(trc_rccsd_lambda_state_t) :: device, other
      integer :: k
      real(dp) :: d1, d2
      character(len=32) :: what
      logical :: want_host

      if (.not. supported()) return
      call ensure_fixture(DIIS_CASE, error)
      if (allocated(error)) return
      call trc_rccsd_lambda_solve(fix(DIIS_CASE)%cc, tight(), la, device)
      call check(error, device%result%converged .and. .not. device%result%diis_host, &
                 "device history: "//trim(device%result%message))
      if (allocated(error)) return
      do k = 1, 3
         opts = tight()
         want_host = .false.
         select case (k)
         case (1)
            opts%diis_host = .true.
            want_host = .true.
            what = "host history"
         case (2)
            ! 200 entries of cc-pVDZ water are 9.6 MB, above 0.7 of the 8 MB the cap leaves; everything else fits
            opts%diis_size = 200
            opts%device_free_cap = 8000000_int64
            want_host = .true.
            what = "history that does not fit"
         case (3)
            opts%diis_size = 0
            what = "no DIIS"
         end select
         call trc_rccsd_lambda_solve(fix(DIIS_CASE)%cc, opts, la, other)
         call check(error, other%result%converged, trim(what)//": "//trim(other%result%message))
         if (allocated(error)) return
         call check(error, other%result%diis_host .eqv. want_host, trim(what)//": the history is not where it should be")
         if (allocated(error)) return
         d1 = maxval(abs(other%l%t1 - device%l%t1))
         d2 = maxval(abs(other%l%t2 - device%l%t2))
         print "(a,a,a,i0,a,i0,a,2es10.2)", "  ", trim(what), ": iterations ", other%result%iterations, " (device ", &
            device%result%iterations, "), max |difference| l1, l2 ", d1, d2
         call check(error, d1 < TOL_SAME .and. d2 < TOL_SAME, trim(what)//" gives other multipliers: "// &
                    to_string(d1, "es9.2")//" "//to_string(d2, "es9.2"))
         if (allocated(error)) return
         if (k == 3) then
            call check(error, other%result%iterations > device%result%iterations, "DIIS saves no iteration")
            if (allocated(error)) return
         end if
      end do
      call device%release()
      call other%release()
   end subroutine test_diis

   subroutine test_nonconvergence(error)
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_lambda_solve_options_t) :: opts
      type(trc_rccsd_lambda_state_t) :: lam
      real(dp) :: kept

      if (.not. supported()) return
      call ensure_fixture(DIIS_CASE, error)
      if (allocated(error)) return
      opts%max_iterations = 2
      call trc_rccsd_lambda_solve(fix(DIIS_CASE)%cc, opts, la, lam)
      call check(error, .not. lam%result%converged .and. lam%result%finite .and. .not. lam%result%refused .and. &
                 lam%result%iterations == 2 .and. index(lam%result%message, "not converged in 2 iterations") > 0, &
                 "two iterations: "//trim(lam%result%message))
      if (allocated(error)) return
      ! one criterion alone must not converge it
      opts = trc_rccsd_lambda_solve_options_t()
      opts%max_iterations = 3
      opts%conv_amplitude = 1.0e10_dp
      call trc_rccsd_lambda_solve(fix(DIIS_CASE)%cc, opts, la, lam)
      call check(error, .not. lam%result%converged, "converged on the max step alone")
      if (allocated(error)) return
      opts%conv_amplitude = 1.0e-6_dp
      opts%conv_norm = 1.0e10_dp
      call trc_rccsd_lambda_solve(fix(DIIS_CASE)%cc, opts, la, lam)
      call check(error, .not. lam%result%converged, "converged on the norm of the step alone")
      if (allocated(error)) return

      associate (t => fix(DIIS_CASE)%cc%amps)
         kept = t%t1(1, 1)
         t%t1(1, 1) = ieee_value(kept, ieee_quiet_nan)
         call t%push()
         opts = trc_rccsd_lambda_solve_options_t()
         call trc_rccsd_lambda_solve(fix(DIIS_CASE)%cc, opts, la, lam)
         t%t1(1, 1) = kept
         call t%push()
      end associate
      call check(error, .not. lam%result%converged .and. .not. lam%result%finite .and. &
                 index(lam%result%message, "NaN or infinite") > 0, "a NaN: "//trim(lam%result%message))
      call lam%release()
   end subroutine test_nonconvergence

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_lambda_solve_options_t) :: opts
      type(trc_rccsd_state_t) :: empty
      real(dp) :: kept

      if (.not. supported()) return
      call ensure_fixture(DIIS_CASE, error)
      if (allocated(error)) return
      associate (cc => fix(DIIS_CASE)%cc)
         call refused(error, empty, opts, "did not converge")
         if (allocated(error)) return
         empty%result%converged = .true.
         call refused(error, empty, opts, "holds no integrals or amplitudes")
         if (allocated(error)) return
         empty%result%refused = .true.
         empty%result%message = "trc_rccsd: an open shell"
         call refused(error, empty, opts, "the CCSD was refused: trc_rccsd: an open shell")
         if (allocated(error)) return
         cc%result%converged = .false.
         call refused(error, cc, opts, "did not converge")
         cc%result%converged = .true.
         if (allocated(error)) return
         allocate (cc%rank_seconds(4, 2))
         call refused(error, cc, opts, "several ranks")
         deallocate (cc%rank_seconds)
         if (allocated(error)) return
         kept = cc%amps%t2(1, 2, 1, 1)
         cc%amps%t2(1, 2, 1, 1) = kept + 1.0e-3_dp
         call refused(error, cc, opts, "t2 is not symmetric")
         cc%amps%t2(1, 2, 1, 1) = kept
         if (allocated(error)) return
         opts%max_iterations = 0
         call refused(error, cc, opts, "max_iterations and diis_start must be >= 1")
         if (allocated(error)) return
         opts = trc_rccsd_lambda_solve_options_t()
         opts%conv_norm = 0.0_dp
         call refused(error, cc, opts, "thresholds must be positive")
         if (allocated(error)) return
         opts = trc_rccsd_lambda_solve_options_t()
         opts%device_free_cap = 1000_int64
         call refused(error, cc, opts, "the amplitudes (o 4, v 19) needs")
         if (allocated(error)) return
         ! room for the amplitude sets, not for what the update builds
         opts%device_free_cap = 60000_int64
         call refused(error, cc, opts, "GB are free")
         if (allocated(error)) return
         ! and the same state still solves
         opts = trc_rccsd_lambda_solve_options_t()
         call solves(error, cc, opts)
      end associate
   end subroutine test_refusals

   subroutine solves(error, cc, opts)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_state_t), intent(in) :: cc
      type(trc_rccsd_lambda_solve_options_t), intent(in) :: opts
      type(trc_rccsd_lambda_state_t) :: lam

      call trc_rccsd_lambda_solve(cc, opts, la, lam)
      call check(error, lam%result%converged, "after the refusals: "//trim(lam%result%message))
      call lam%release()
   end subroutine solves

   subroutine refused(error, cc, opts, fragment)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_state_t), intent(in) :: cc
      type(trc_rccsd_lambda_solve_options_t), intent(in) :: opts
      character(len=*), intent(in) :: fragment
      type(trc_rccsd_lambda_state_t) :: lam

      call trc_rccsd_lambda_solve(cc, opts, la, lam)
      call check(error, lam%result%refused .and. .not. lam%result%converged .and. &
                 index(lam%result%message, fragment) > 0, &
                 "expected a refusal containing '"//fragment//"', got '"//trim(lam%result%message)//"'")
      if (allocated(error)) return
      call check(error, .not. allocated(lam%l%t1), "a refused solve left multipliers allocated")
   end subroutine refused

   subroutine ensure_fixture(ic, error)
      !! The SCF and the RI-CCSD of case `ic`, converged as the reference's
      integer, intent(in) :: ic
      type(error_type), allocatable, intent(out) :: error

      type(trc_scf_result_t) :: res
      type(trc_rccsd_options_t) :: opts
      type(error_t) :: err
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm, nocc, dev

      if (.not. fix(ic)%built) then
         if (.not. bound) then
            bound = .true.
            dev = trc_bind_device(0)
            call la%init(2)
         end if
         fix(ic)%built = .true.
         call read_ccsd_grad_ref(trim(CASE_REF(ic)), fix(ic)%ref)
         call geometry(natm, zint, at_r)
         nocc = sum(zint)/2
         associate (f => fix(ic))
            call trc_basis_from_json("basis_sets/"//trim(CASE_BASIS(ic))//".json", natm, zint, at_r, f%bas, err)
            if (err%has_error()) error stop "check_rccsd_lambda: "//err%get_message()
            ! the auxiliary basis follows the orbital basis's convention, as through the C API
            call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, f%aux, err, &
                                     force_cart=.not. f%bas%spherical)
            if (err%has_error()) error stop "check_rccsd_lambda: "//err%get_message()
            call f%bas%to_device()
            call f%aux%to_device()
            call f%pl%build(f%bas, 1.0e-12_dp)
            call f%pl%to_device()
            call run_scf(f%bas, nocc, res)
            opts%nocc = nocc
            opts%nfrozen = CASE_NFROZEN(ic)
            opts%conv_energy = 1.0e-13_dp
            opts%conv_amplitude = 1.0e-11_dp
            call trc_rccsd_solve(f%bas, f%aux, f%pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, f%cc)
         end associate
      end if
      call check(error, fix(ic)%cc%result%converged, trim(CASE_REF(ic))//": CCSD: "//trim(fix(ic)%cc%result%message))
      if (allocated(error)) return
      call check(error, abs(fix(ic)%cc%result%e_ccsd - fix(ic)%ref%e_ccsd_corr) < 1.0e-10_dp .and. &
                 fix(ic)%cc%result%no == fix(ic)%ref%no .and. fix(ic)%cc%result%naux == fix(ic)%ref%naux, &
                 trim(CASE_REF(ic))//": not the reference's CCSD")
   end subroutine ensure_fixture

   subroutine release_fixture()
      integer :: n

      do n = 1, NCASE
         call fix(n)%cc%release()
         call fix(n)%pl%release()
         call fix(n)%aux%release()
         call fix(n)%bas%release()
      end do
      call la%release()
   end subroutine release_fixture

end module check_rccsd_lambda_suite

program check_rccsd_lambda
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_lambda_suite, only: collect_check_rccsd_lambda, release_fixture
   use ccsd_grad_support, only: supported
   implicit none
   integer :: stat
   stat = 0
   if (.not. supported()) then
      print "(a)", "SKIP check_rccsd_lambda: this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
      stop 0
   end if
   call run_testsuite(collect_check_rccsd_lambda, error_unit, stat, parallel=.false.)
   call release_fixture()
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_lambda
