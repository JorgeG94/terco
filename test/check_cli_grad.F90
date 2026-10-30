!
! The CLI's --grad: what `terco` prints to its JSON summary must be
! exactly what the library itself (trc_run_gradient/trc_gradient, the
! same C API app/main.F90 calls) computes for the identical molecule,
! basis and guess -- app/main.F90 adds no arithmetic of its own between
! `trc_gradient` and `write_json_gradient`, only a reshape, so this is a
! wiring check, not a physics one (check_full_grad_fd.F90 and
! check_full_grad_vs_pyscf.F90 are that).
!
! Runs the installed `terco` binary (its path given as this program's
! first command-line argument -- see test/CMakeLists.txt's `$<TARGET_FILE:
! terco_cli>`; its second argument is a scratch path for the CLI's own
! --json output, in the BINARY dir rather than hardcoded here, the same
! reason check_1e_grad_vs_pyscf.F90's own `dump` argument is) with
! --guess gwh (deterministic, no SAD atomic sub-SCF to keep in sync) on
! water/6-31G, reads its --json gradient array back, and
! compares it to a direct trc_create/.../trc_run_gradient/trc_gradient
! call on the SAME molecule/basis/guess/convergence settings, in the SAME
! order app/main.F90 itself calls them. Both paths run the identical code
! underneath, but as two separate processes on a shared host, so the
! tolerance below is calibrated to a measured threaded-BLAS noise floor,
! not machine precision -- see TOL's own comment.
!
module check_cli_grad_suite
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                       trc_set_method, trc_set_guess, trc_set_convergence, trc_set_screening, &
                       trc_set_incremental, trc_run_scf, trc_run_gradient, trc_gradient, &
                       trc_set_verbose, trc_set_force_cartesian, trc_set_nlc_grid, &
                       TRC_OK, TRC_GUESS_GWH
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_cli_grad_tests, set_cli_grad_args

   real(dp), parameter :: TOL = 1.0e-6_dp
      !! NOT 1e-10, even though both paths call the identical setter sequence
      !! with identical values: they are two separate OS processes sharing
      !! this (4-core) host, and the threaded BLAS under `trc_linalg_t` is not
      !! associative, so their GWH guesses and nine DIIS iterations can land
      !! at very slightly different (but both gconv=1e-5-converged) densities.
      !! Measured: two such runs' total energies agreed to ~6e-12 (consistent
      !! with that convergence gate) while their GRADIENTS -- first
      !! derivatives, linearly sensitive to the same residual density error
      !! an energy comparison is only quadratically sensitive to -- differed
      !! by ~5e-10. 1e-6 stays far above that measured noise floor while
      !! still catching a real wiring bug (a reshape, an atom-order mixup, a
      !! sign flip) by several orders of magnitude.

   character(len=:), allocatable :: g_cli_exe, g_probe_json
      !! Set once from the command line (test-drive's new_unittest takes no
      !! arguments) by set_cli_grad_args, before run_testsuite runs the
      !! single case below.

contains

   subroutine set_cli_grad_args(cli_exe, probe_json)
      character(*), intent(in) :: cli_exe, probe_json
      g_cli_exe = cli_exe
      g_probe_json = probe_json
   end subroutine set_cli_grad_args

   subroutine collect_cli_grad_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("CLI --json gradient matches a direct trc_gradient call", test_cli_grad) &
                  ]
   end subroutine collect_cli_grad_tests

   subroutine test_cli_grad(error)
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
      ! Same order app/main.F90 itself calls these in, not just the same
      ! values -- found by measurement to matter for the last few digits
      ! (see this module's own header for the reason once it is pinned down).
      st = trc_set_verbose(h, 1_c_int)
      call check(error, st == TRC_OK, "trc_set_verbose failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_molecule(h, int(nat, c_int), real(at_z, c_double), r, 0_c_int, 1_c_int)
      call check(error, st == TRC_OK, "trc_set_molecule failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_force_cartesian(h, 0_c_int)
      call check(error, st == TRC_OK, "trc_set_force_cartesian failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_basis_json(h, c_string_local('basis_sets/6-31g.json'))
      call check(error, st == TRC_OK, "trc_set_basis_json failed, status "//to_string(int(st)))
      if (allocated(error)) return
      st = trc_set_method(h, c_string_local(''), 3_c_int)
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

      ! --- the CLI's own printed/JSON gradient, same molecule/basis/guess --
      cmd = g_cli_exe//' water.xyz -b basis_sets/6-31g.json --guess gwh --grad --json '// &
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
   end subroutine test_cli_grad

   pure function c_string_local(s) result(cs)
      !! `s` NUL-terminated, for a `character(kind=c_char) :: x(*)` dummy --
      !! the same convention app/trc_app.F90's own `c_string` uses, copied
      !! here rather than pulled in: that module is compiled straight into
      !! the `terco_cli` executable (app/CMakeLists.txt), not a library
      !! test/ links against.
      character(*), intent(in) :: s
      character(len=len(s) + 1) :: cs
      cs = s//achar(0)
   end function c_string_local

end module check_cli_grad_suite

program check_cli_grad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_cli_grad_suite, only: collect_cli_grad_tests, set_cli_grad_args
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
      write (error_unit, '(a)') 'check_cli_grad: usage: check_cli_grad <terco CLI binary> <scratch .json path>'
      error stop 1
   end if
   call set_cli_grad_args(cli_exe, probe_json)

   stat = 0
   testsuites = [new_testsuite("check_cli_grad", collect_cli_grad_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_grad
