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
! The double hybrid: `-f b2plyp --df --grad` (water, Cartesian cc-pVDZ, one
! cc-pVDZ-RIFIT doing both fits, the default frozen core) against pyscf's
! own finite-difference gradient (check_dh_grad.F90's G_PYSCF_WATER, its
! header has the provenance), and its refusals: an unrestricted (UKS)
! reference and PCM.
!
! And for both RI-MP2 references, that `--grad` runs no separate RI-MP2
! energy step -- the energy comes from the gradient's own amplitudes
! (trc_run_rimp2 defers it) -- yet reports the same "mp2 opposite-spin
! energy" and "mp2 same-spin energy" as the run without `--grad` to 1e-12;
! "no step" is read off the --logger verbose output, where the energy
! step's own table ("rimp2:" lines) appears only in the energy-only run.
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

   real(dp), parameter :: ENERGY_TOL = 1.0e-12_dp
   !! --grad's E_os/E_ss against the energy-only run's: the same B and
   !! pairs, round-off apart.
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

   real(dp), parameter :: G_PYSCF_DH(3, 3) = reshape([ &
                                          -1.894780628694e-11_dp, -1.515824502955e-10_dp, -1.332378636339e-02_dp, &
                                          3.552713678801e-10_dp, -3.896076487801e-03_dp, 6.661893121892e-03_dp, &
                                          1.421085471520e-11_dp, 3.896076588461e-03_dp, 6.661893099391e-03_dp], [3, 3])
   !! pyscf 2.14, b2plyp --df: check_dh_grad.F90's G_PYSCF_WATER.
   real(dp), parameter :: TOL_DH = 5.0e-9_dp
   !! check_dh_grad.F90's own gate against the same numbers.

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
                  new_unittest("CLI --df --rimp2 --grad JSON gradient matches metalquicha", test_cli_rimp2_grad_df), &
                  new_unittest("CLI --rimp2 --grad: no energy step, the energy-only run's E_os/E_ss", &
                               test_cli_grad_energy), &
                  new_unittest("CLI --df --rimp2 --grad: no energy step, the energy-only run's E_os/E_ss", &
                               test_cli_grad_energy_df), &
                  new_unittest("CLI --rimp2 --grad over a UHF reference is refused", test_cli_rimp2_grad_uhf_refused), &
                  new_unittest("CLI -f b2plyp --df --grad JSON gradient matches pyscf", test_cli_dh_grad_df), &
                  new_unittest("CLI -f b2plyp --grad over a UKS reference is refused", test_cli_dh_grad_uks_refused), &
                  new_unittest("CLI -f b2plyp --pcm --grad is refused", test_cli_dh_grad_pcm_refused) &
                  ]
   end subroutine collect_cli_grad_rimp2_tests

   subroutine test_cli_rimp2_grad(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: grad(3, 3)
      call cli_gradient(error, ' water_mqc.xyz -b basis_sets/cc-pvdz.json --rimp2 --aux basis_sets/cc-pvdz-rifit.json'// &
                        ' --frozen 0', g_probe_json, 3, grad, "scf energy + mp2 correlation")
      if (allocated(error)) return
      call check(error, maxval(abs(grad - G_MQC)), 0.0_dp, thr=TOL, &
                 message="worst |CLI - metalquicha| = "//to_string(maxval(abs(grad - G_MQC)), "es10.2"))
   end subroutine test_cli_rimp2_grad

   subroutine test_cli_rimp2_grad_df(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: grad(3, 4)
      call cli_gradient(error, ' nh3_mqc.xyz -b basis_sets/def2-svp.json --df --jkfit basis_sets/def2-svp-rifit.json'// &
                        ' --rimp2 --aux basis_sets/def2-svp-rifit.json --frozen 0', g_probe_json//'.df.json', 4, grad, &
                        "scf energy + mp2 correlation")
      if (allocated(error)) return
      call check(error, maxval(abs(grad - G_MQC_DF)), 0.0_dp, thr=TOL, &
                 message="worst |CLI --df - metalquicha| = "//to_string(maxval(abs(grad - G_MQC_DF)), "es10.2"))
   end subroutine test_cli_rimp2_grad_df

   subroutine test_cli_grad_energy(error)
      type(error_type), allocatable, intent(out) :: error
      call grad_energy_matches(error, ' water_mqc.xyz -b basis_sets/cc-pvdz.json --rimp2'// &
                               ' --aux basis_sets/cc-pvdz-rifit.json', g_probe_json//'.e')
   end subroutine test_cli_grad_energy

   subroutine test_cli_grad_energy_df(error)
      type(error_type), allocatable, intent(out) :: error
      call grad_energy_matches(error, ' nh3_mqc.xyz -b basis_sets/def2-svp.json --df --jkfit basis_sets/def2-svp-rifit.json'// &
                               ' --rimp2 --aux basis_sets/def2-svp-rifit.json', g_probe_json//'.edf')
   end subroutine test_cli_grad_energy_df

   subroutine grad_energy_matches(error, args, probe)
      !! `args` with and without --grad (both --logger verbose): the same
      !! E_os and E_ss, and the energy step's table only in the run without.
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: args, probe
      real(dp) :: os_ref, ss_ref, os_grad, ss_grad
      logical :: table_ref, table_grad
      call cli_energy(error, args, probe//'.energy.json', os_ref, ss_ref, table_ref)
      if (allocated(error)) return
      call check(error, table_ref, "the energy-only run logs no RI-MP2 energy table (the probe is blind)")
      if (allocated(error)) return
      call cli_energy(error, args//' --grad', probe//'.grad.json', os_grad, ss_grad, table_grad)
      if (allocated(error)) return
      call check(error, .not. table_grad, "--grad ran the separate RI-MP2 energy step")
      if (allocated(error)) return
      call check(error, max(abs(os_grad - os_ref), abs(ss_grad - ss_ref)), 0.0_dp, thr=ENERGY_TOL, &
                 message="--grad's E_os/E_ss differ from the energy-only run's by "// &
                 to_string(max(abs(os_grad - os_ref), abs(ss_grad - ss_ref)), "es10.2"))
   end subroutine grad_energy_matches

   subroutine cli_energy(error, args, probe, e_os, e_ss, energy_table)
      !! Run the CLI with `args` plus --frozen 0, tight convergence, --logger
      !! verbose and --json `probe`; read its MP2 energies back, and whether
      !! its output carries trc_rimp2_run's own step table.
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: args, probe
      real(dp), intent(out) :: e_os, e_ss
      logical, intent(out) :: energy_table
      character(len=:), allocatable :: cmd
      character(len=512) :: line
      type(json_file) :: json
      logical :: found_os, found_ss
      integer :: rc, u, ios

      e_os = 0.0_dp
      e_ss = 0.0_dp
      energy_table = .false.
      cmd = g_cli_exe//args//' --frozen 0 --econv 1e-12 --gconv 1e-9 --logger verbose --json '//probe//' > '// &
            probe//'.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call check(error, rc == 0, "the CLI run failed, exit code "//to_string(rc)//" (see "//probe//".out)")
      if (allocated(error)) return
      call json%initialize()
      call json%load(filename=probe)
      call json%get('mp2 opposite-spin energy', e_os, found_os)
      call json%get('mp2 same-spin energy', e_ss, found_ss)
      call json%destroy()
      call check(error, found_os .and. found_ss, "no MP2 energies in "//probe)
      if (allocated(error)) return
      open (newunit=u, file=probe//'.out', status="old", action="read", iostat=ios)
      call check(error, ios == 0, "cannot read "//probe//".out")
      if (allocated(error)) return
      do
         read (u, "(a)", iostat=ios) line
         if (ios /= 0) exit
         if (index(line, " rimp2:") > 0) energy_table = .true.
      end do
      close (u)
   end subroutine cli_energy
   subroutine test_cli_rimp2_grad_uhf_refused(error)
      !! An unrestricted reference has no RI-MP2 gradient: the run must stop
      !! with a non-zero exit code and say why, not print an SCF-only gradient.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out
      character(len=512) :: line
      integer :: rc, u, ios
      out = g_probe_json//'.uhf.out'
      cmd = g_cli_exe//' water.xyz -b basis_sets/6-31g.json --uhf 2 --rimp2 --grad > '//out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call check(error, rc /= 0, "a UHF reference with --rimp2 --grad was not refused (exit code 0, see "//out//")")
      if (allocated(error)) return
      rc = 0
      open (newunit=u, file=out, status="old", action="read", iostat=ios)
      do while (ios == 0)
         read (u, "(a)", iostat=ios) line
         if (ios == 0 .and. index(line, "restricted (closed-shell RHF reference) only") > 0) rc = 1
      end do
      close (u)
      call check(error, rc == 1, "the refusal does not say the RI-MP2 gradient is restricted only (see "//out//")")
   end subroutine test_cli_rimp2_grad_uhf_refused

   subroutine test_cli_dh_grad_df(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: grad(3, 3)
      call cli_gradient(error, ' water.xyz -b basis_sets/cc-pvdz.json --cart -f b2plyp --df --jkfit '// &
                        'basis_sets/cc-pvdz-rifit.json --aux basis_sets/cc-pvdz-rifit.json', &
                        g_probe_json//'.dh.json', 3, grad, &
                        "double hybrid total: scf energy + pt2 fraction x mp2 correlation")
      if (allocated(error)) return
      call check(error, maxval(abs(grad - G_PYSCF_DH)), 0.0_dp, thr=TOL_DH, &
                 message="worst |CLI b2plyp --df - pyscf| = "//to_string(maxval(abs(grad - G_PYSCF_DH)), "es10.2"))
   end subroutine test_cli_dh_grad_df

   subroutine test_cli_dh_grad_uks_refused(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_refusal(error, ' oh.xyz -b basis_sets/6-31g.json --uhf 1 -f b2plyp --grad', '.uks.out', &
                          "an unrestricted (UKS) reference has no correlation gradient")
   end subroutine test_cli_dh_grad_uks_refused

   subroutine test_cli_dh_grad_pcm_refused(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_refusal(error, ' water.xyz -b basis_sets/6-31g.json -f b2plyp --pcm cpcm --eps 78.3553 --grad', &
                          '.pcm.out', "PCM gradient not implemented")
   end subroutine test_cli_dh_grad_pcm_refused

   subroutine expect_refusal(error, args, suffix, reason)
      !! The CLI with `args` must exit non-zero and print `reason`.
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: args, suffix, reason
      character(len=:), allocatable :: out
      character(len=512) :: line
      integer :: rc, u, ios
      out = g_probe_json//suffix
      call execute_command_line(g_cli_exe//args//' > '//out//' 2>&1', exitstat=rc, wait=.true.)
      call check(error, rc /= 0, "not refused (exit code 0, see "//out//")")
      if (allocated(error)) return
      rc = 0
      open (newunit=u, file=out, status="old", action="read", iostat=ios)
      do while (ios == 0)
         read (u, "(a)", iostat=ios) line
         if (ios == 0 .and. index(line, reason) > 0) rc = 1
      end do
      close (u)
      call check(error, rc == 1, "the refusal does not say '"//reason//"' (see "//out//")")
   end subroutine expect_refusal

   subroutine cli_gradient(error, args, probe, nat, grad, expect_of)
      !! Run the CLI with `args` plus tight convergence, --grad and --json
      !! `probe`, read its gradient back, and check its "gradient of" says
      !! `expect_of`.
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: args, probe, expect_of
      integer, intent(in) :: nat
      real(dp), intent(out) :: grad(3, nat)
      character(len=:), allocatable :: cmd, what
      type(json_file) :: json
      real(dp), allocatable :: row(:)
      logical :: found
      integer :: rc, i

      grad = 0.0_dp
      cmd = g_cli_exe//args//' --econv 1e-12 --gconv 1e-9 --grad --json '//probe//' > '// &
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
      call check(error, what == expect_of, &
                 "the JSON does not say what the gradient differentiates ('"//expect_of//"'): '"//what//"'")
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
