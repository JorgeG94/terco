!
! The CLI's --grad --df: what `terco` prints to its JSON summary for a
! density-fitted run must be exactly what the library itself
! (trc_run_gradient/trc_gradient over the context's own fit, the same C API
! app/main.F90 calls) computes for the identical molecule, basis, jkfit
! basis, functional and guess -- check_cli_grad.F90's wiring check, for the
! fitted gradient (trc_scf_grad_df.F90; check_df_grad_fd.F90 and
! check_df_grad_vs_pyscf.F90 are its physics checks). B3LYP, so the fitted
! J, the fitted K (scaled by 0.2) and the XC term all run, on water/6-31G*
! with def2-universal-jkfit, Cartesian as --df requires. Same tolerance and
! reasoning as check_cli_grad.F90's TOL.
!
module check_cli_grad_df_suite
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                       trc_set_method, trc_set_guess, trc_set_convergence, trc_set_screening, &
                       trc_set_incremental, trc_run_scf, trc_run_gradient, trc_gradient, &
                       trc_set_verbose, trc_set_force_cartesian, trc_set_nlc_grid, &
                       trc_set_df, trc_set_jkfit_json, &
                       TRC_OK, TRC_GUESS_GWH
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_cli_grad_df_tests, set_cli_grad_df_args

   real(dp), parameter :: TOL = 1.0e-6_dp
      !! Same value and reasoning as check_cli_grad.F90's own TOL: two
      !! separate OS processes, threaded BLAS that is not associative, so
      !! their GWH guesses and DIIS iterations (now also the XC quadrature's
      !! own batched GEMMs) can land at very slightly different densities.

   character(len=:), allocatable :: g_cli_exe, g_probe_json
      !! Set once from the command line (test-drive's new_unittest takes no
      !! arguments) by set_cli_grad_df_args, before run_testsuite runs the
      !! single case below.

contains

   subroutine set_cli_grad_df_args(cli_exe, probe_json)
      character(*), intent(in) :: cli_exe, probe_json
      g_cli_exe = cli_exe
      g_probe_json = probe_json
   end subroutine set_cli_grad_df_args

   subroutine collect_cli_grad_df_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("CLI --df --grad b3lyp JSON gradient matches a direct trc_gradient call", &
                               test_cli_grad_df) &
                  ]
   end subroutine collect_cli_grad_df_tests

   subroutine test_cli_grad_df(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: probe_out, cmd
      integer :: nat, i, rc
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: grad_api(:, :), grad_cli(:, :)
      type(c_ptr) :: h
      integer(c_int) :: st
      logical :: found
      type(json_file) :: json
      real(dp), allocatable :: row(:)
      real(dp) :: worst

      probe_out = g_probe_json//'.out'
      call read_xyz('water.xyz', nat, at_z, r)

      ! --- "the library's" own gradient, straight through the C API -------
      st = trc_create(h)
      call check(error, st == TRC_OK, "trc_create failed")
      if (allocated(error)) return
      st = trc_set_verbose(h, 1_c_int)
      call check(error, st == TRC_OK, "trc_set_verbose failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_molecule(h, int(nat, c_int), real(at_z, c_double), r, 0_c_int, 1_c_int)
      call check(error, st == TRC_OK, "trc_set_molecule failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_force_cartesian(h, 1_c_int)
      call check(error, st == TRC_OK, "trc_set_force_cartesian failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_basis_json(h, c_string_local('basis_sets/6-31gs.json'))
      call check(error, st == TRC_OK, "trc_set_basis_json failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_jkfit_json(h, c_string_local('basis_sets/def2-universal-jkfit.json'))
      call check(error, st == TRC_OK, "trc_set_jkfit_json failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_df(h, 1_c_int)
      call check(error, st == TRC_OK, "trc_set_df failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_method(h, c_string_local('b3lyp'), 3_c_int)
      call check(error, st == TRC_OK, "trc_set_method failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_nlc_grid(h, 1_c_int)
      call check(error, st == TRC_OK, "trc_set_nlc_grid failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-5_c_double, 100_c_int, 10_c_int)
      call check(error, st == TRC_OK, "trc_set_convergence failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_screening(h, 1.0e-12_c_double)
      call check(error, st == TRC_OK, "trc_set_screening failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_incremental(h, 0_c_int)
      call check(error, st == TRC_OK, "trc_set_incremental failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
      call check(error, st == TRC_OK, "trc_set_guess failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_run_scf(h)
      call check(error, st == TRC_OK, "trc_run_scf failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_run_gradient(h)
      call check(error, st == TRC_OK, "trc_run_gradient failed, status "//to_string(int(st)))
      if (allocated(error)) return
      allocate (grad_api(3, nat))
      st = trc_gradient(h, grad_api)
      call check(error, st == TRC_OK, "trc_gradient failed, status "//to_string(int(st)))
      if (allocated(error)) then
         st = trc_destroy(h)
         return
      end if
      st = trc_destroy(h)

      ! --- the CLI's own printed/JSON gradient, same molecule/basis/func ---
      cmd = g_cli_exe//' water.xyz -b basis_sets/6-31gs.json --df --jkfit basis_sets/def2-universal-jkfit.json '// &
            '--func b3lyp --guess gwh --grad --json '// &
            g_probe_json//' > '//probe_out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call check(error, rc == 0, "the CLI run itself failed, exit code "//to_string(rc))
      if (allocated(error)) return

      allocate (grad_cli(3, nat))
      call json%initialize()
      call json%load(filename=g_probe_json)
      call check(error,.not. json%failed(), "cannot read "//g_probe_json)
      if (allocated(error)) then
         call json%destroy()
         return
      end if
      do i = 1, nat
         call json%get('gradient('//to_string(i)//')', row, found)
         if (.not. found .or. size(row) /= 3) then
            call json%destroy()
            call check(error, .false., "no gradient row "//to_string(i)//" in the CLI's JSON output")
            return
         end if
         grad_cli(:, i) = row
      end do
      call json%destroy()

      worst = maxval(abs(grad_api - grad_cli))
      call check(error, worst, 0.0_dp, thr=TOL, &
                 message="worst |CLI JSON - direct trc_gradient| = "//to_string(worst, "es10.2")// &
                 " exceeds TOL = "//to_string(TOL, "es10.2"))
   end subroutine test_cli_grad_df

   pure function c_string_local(s) result(cs)
      !! `s` NUL-terminated, for a `character(kind=c_char) :: x(*)` dummy --
      !! the same convention app/trc_app.F90's own `c_string` uses.
      character(*), intent(in) :: s
      character(len=len(s) + 1) :: cs
      cs = s//achar(0)
   end function c_string_local

end module check_cli_grad_df_suite

program check_cli_grad_df
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_cli_grad_df_suite, only: collect_cli_grad_df_tests, set_cli_grad_df_args
   implicit none
   character(len=1024) :: arg
   character(len=:), allocatable :: cli_exe, probe_json
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat

   call get_command_argument(1, arg)
   cli_exe = trim(arg)
   call get_command_argument(2, arg)
   probe_json = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0) then
      write (error_unit, '(a)') &
         'check_cli_grad_df: usage: check_cli_grad_df <terco CLI binary> <scratch .json path>'
      error stop 1
   end if
   call set_cli_grad_df_args(cli_exe, probe_json)

   stat = 0
   testsuites = [new_testsuite("check_cli_grad_df", collect_cli_grad_df_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_grad_df
