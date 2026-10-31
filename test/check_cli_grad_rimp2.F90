!
! The CLI's `--rimp2 --grad`: the whole path a user runs (app/main.F90 ->
! trc_run_scf, trc_run_rimp2, trc_run_gradient -> trc_rimp2_grad), its
! --json gradient against metalquicha's pyscf-validated RI-MP2 gradient for
! the same water geometry, cc-pVDZ / cc-pVDZ-RIFIT, no frozen core
! (metalquicha validation/validation_tests_cpu.json, "RI-MP2 gradient H2O
! cc-pvdz/cc-pvdz-rifit frozen 0"; provenance in check_rimp2_grad.F90's
! header), and the refusal of a --df reference, which has no RI-MP2
! gradient yet. The CLI keeps its default (pruned) reading of the basis:
! pruning preserves the span, so the energy and gradient do not move
! (measured: 1e-10 against the unpruned reference).
!
! Arguments: the terco executable and a scratch path for its --json output.
!
module check_cli_grad_rimp2_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_cli_grad_rimp2_tests, set_cli_grad_rimp2_args

   real(dp), parameter :: TOL = 2.0e-8_dp
   !! The library-level comparison's own (check_rimp2_grad.F90), measured
   !! 1e-10 through the CLI at --econv 1e-12 --gconv 1e-9.
   real(dp), parameter :: G_MQC(3, 3) = reshape([ &
                                       -0.0_dp, 2.31e-10_dp, -0.015717146972_dp, &
                                       0.0_dp, 0.002913907061_dp, 0.007858573573_dp, &
                                       -0.0_dp, -0.002913907293_dp, 0.0078585734_dp], [3, 3])

   character(len=:), allocatable :: g_cli_exe, g_probe_json

contains

   subroutine set_cli_grad_rimp2_args(cli_exe, probe_json)
      character(*), intent(in) :: cli_exe, probe_json
      g_cli_exe = cli_exe
      g_probe_json = probe_json
   end subroutine set_cli_grad_rimp2_args

   subroutine collect_cli_grad_rimp2_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("CLI --rimp2 --grad JSON gradient matches metalquicha", test_cli_rimp2_grad), &
                  new_unittest("CLI --rimp2 --grad over --df is refused", test_cli_rimp2_grad_df_refused) &
                  ]
   end subroutine collect_cli_grad_rimp2_tests

   subroutine test_cli_rimp2_grad(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe_out, what
      type(json_file) :: json
      real(dp), allocatable :: row(:)
      real(dp) :: grad(3, 3), worst
      logical :: found
      integer :: rc, i

      probe_out = g_probe_json//'.out'
      cmd = g_cli_exe//' water_mqc.xyz -b basis_sets/cc-pvdz.json --rimp2 --aux basis_sets/cc-pvdz-rifit.json'// &
            ' --frozen 0 --econv 1e-12 --gconv 1e-9 --grad --json '//g_probe_json//' > '//probe_out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call check(error, rc == 0, "the CLI run failed, exit code "//to_string(rc)//" (see "//probe_out//")")
      if (allocated(error)) return
      call json%initialize()
      call json%load(filename=g_probe_json)
      call check(error, .not. json%failed(), "cannot read "//g_probe_json)
      if (allocated(error)) then
         call json%destroy()
         return
      end if
      call json%get('gradient of', what, found)
      if (.not. found) what = ""
      do i = 1, 3
         call json%get('gradient('//to_string(i)//')', row, found)
         if (.not. found .or. size(row) /= 3) then
            call json%destroy()
            call check(error, .false., "no gradient row "//to_string(i)//" in the CLI's JSON output")
            return
         end if
         grad(:, i) = row
      end do
      call json%destroy()
      call check(error, what == "scf energy + mp2 correlation", &
                 "the JSON does not say the gradient includes the correlation: '"//what//"'")
      if (allocated(error)) return
      worst = maxval(abs(grad - G_MQC))
      call check(error, worst, 0.0_dp, thr=TOL, &
                 message="worst |CLI - metalquicha| = "//to_string(worst, "es10.2")//" exceeds "//to_string(TOL, "es10.2"))
   end subroutine test_cli_rimp2_grad

   subroutine test_cli_rimp2_grad_df_refused(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe_out
      character(len=512) :: line
      integer :: rc, u, ios
      logical :: said

      probe_out = g_probe_json//'.df.out'
      cmd = g_cli_exe//' water.xyz -b basis_sets/cc-pvdz.json --rimp2 --aux basis_sets/cc-pvdz-rifit.json'// &
            ' --df --grad > '//probe_out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call check(error, rc /= 0, "the CLI ran --rimp2 --grad over a --df reference instead of refusing it")
      if (allocated(error)) return
      said = .false.
      open (newunit=u, file=probe_out, status='old', action='read', iostat=ios)
      do while (ios == 0)
         read (u, '(a)', iostat=ios) line
         if (ios == 0) said = said .or. index(line, "over a --df (density-fitted) SCF reference") > 0
      end do
      close (u)
      call check(error, said, "the refusal does not name the --df reference (see "//probe_out//")")
   end subroutine test_cli_rimp2_grad_df_refused

end module check_cli_grad_rimp2_suite

program check_cli_grad_rimp2
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_grad_rimp2_suite, only: collect_cli_grad_rimp2_tests, set_cli_grad_rimp2_args
   implicit none
   character(len=1024) :: cli_exe, probe
   integer :: stat
   call get_command_argument(1, cli_exe)
   call get_command_argument(2, probe)
   call set_cli_grad_rimp2_args(trim(cli_exe), trim(probe))
   stat = 0
   call run_testsuite(collect_cli_grad_rimp2_tests, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_grad_rimp2
