!
! The CLI's `--rimp2 --grad`: the whole path a user runs (app/main.F90 ->
! trc_run_scf, trc_run_rimp2, trc_run_gradient -> trc_rimp2_grad), its
! --json gradient against metalquicha's pyscf-validated RI-MP2 gradient for
! the same water geometry, cc-pVDZ / cc-pVDZ-RIFIT, no frozen core
! (metalquicha validation/validation_tests_cpu.json, "RI-MP2 gradient H2O
! cc-pvdz/cc-pvdz-rifit frozen 0"; provenance in check_rimp2_grad.F90's
! header), and `--df --rimp2 --grad` against metalquicha's "RI-MP2 gradient
! over a fitted reference NH3 def2-svp/def2-svp-rifit" (one set doing both
! fits, `--jkfit` = `--aux`). The CLI keeps its default (pruned) reading of the basis:
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
   real(dp), parameter :: G_MQC_DF(3, 4) = reshape([ &
                                          -4e-11_dp, -2.9e-11_dp, 0.008288657367_dp, &
                                          -0.002408829768_dp, -2.38e-10_dp, -0.002762885685_dp, &
                                          0.001204414513_dp, -0.00208610766_dp, -0.002762886289_dp, &
                                          0.001204414784_dp, 0.00208610774_dp, -0.002762885861_dp], [3, 4])
   !! check_rimp2_grad.F90's own --df reference (metalquicha, PySCF finite
   !! differences), held there to 2e-8.
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
                  new_unittest("CLI --df --rimp2 --grad JSON gradient matches metalquicha", test_cli_rimp2_grad_df) &
                  ]
   end subroutine collect_cli_grad_rimp2_tests

   subroutine test_cli_rimp2_grad(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: grad(3, 3)
      call cli_gradient(error, ' water_mqc.xyz -b basis_sets/cc-pvdz.json --rimp2 --aux basis_sets/cc-pvdz-rifit.json', &
                        g_probe_json, 3, grad)
      if (allocated(error)) return
      call check(error, maxval(abs(grad - G_MQC)), 0.0_dp, thr=TOL, &
                 message="worst |CLI - metalquicha| = "//to_string(maxval(abs(grad - G_MQC)), "es10.2"))
   end subroutine test_cli_rimp2_grad

   subroutine test_cli_rimp2_grad_df(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: grad(3, 4)
      call cli_gradient(error, ' nh3_mqc.xyz -b basis_sets/def2-svp.json --df --jkfit basis_sets/def2-svp-rifit.json'// &
                        ' --rimp2 --aux basis_sets/def2-svp-rifit.json', g_probe_json//'.df.json', 4, grad)
      if (allocated(error)) return
      call check(error, maxval(abs(grad - G_MQC_DF)), 0.0_dp, thr=TOL, &
                 message="worst |CLI --df - metalquicha| = "//to_string(maxval(abs(grad - G_MQC_DF)), "es10.2"))
   end subroutine test_cli_rimp2_grad_df

   subroutine cli_gradient(error, args, probe, nat, grad)
      !! Run the CLI with `args` plus --frozen 0, tight convergence, --grad
      !! and --json `probe`, and read its gradient back.
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: args, probe
      integer, intent(in) :: nat
      real(dp), intent(out) :: grad(3, nat)
      character(len=:), allocatable :: cmd, what
      type(json_file) :: json
      real(dp), allocatable :: row(:)
      logical :: found
      integer :: rc, i

      grad = 0.0_dp
      cmd = g_cli_exe//args//' --frozen 0 --econv 1e-12 --gconv 1e-9 --grad --json '//probe//' > '// &
            probe//'.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call check(error, rc == 0, "the CLI run failed, exit code "//to_string(rc)//" (see "//probe//".out)")
      if (allocated(error)) return
      call json%initialize()
      call json%load(filename=probe)
      call check(error, .not. json%failed(), "cannot read "//probe)
      if (allocated(error)) then
         call json%destroy()
         return
      end if
      call json%get('gradient of', what, found)
      if (.not. found) what = ""
      do i = 1, nat
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
   end subroutine cli_gradient

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
