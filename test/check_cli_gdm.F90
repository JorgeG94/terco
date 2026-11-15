module check_cli_gdm_suite
   !!
   !! Geometric direct minimisation through the terco CLI and its TOML input, as a
   !! test-drive suite. The output line asserted on is one only GDM produces (the
   !! handover note), so a flag or key that never reached the SCF cannot pass.
   !!
   !! 1. `--gdm` hands over from DIIS (the note is printed, and a run without the
   !!    flag does not print it: the control for the text reading itself).
   !! 2. `--gdm-start 1e-12`, far below where DIIS converges, never hands over;
   !!    `--gdm-max-rotation 1e-4` still converges but takes more iterations than
   !!    the default cap.
   !! 3. The same three through `[scf.convergence]` (`gdm`, `gdm_start`,
   !!    `gdm_max_rotation`) give the same behaviour as the flags.
   !! 4. Bad values and `--gdm --soscf` are refused, by name, before any SCF.
   !!
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   implicit none
   private
   public :: collect_cli_gdm, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
   !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
   !! A scratch directory this check owns.

   character(len=*), parameter :: HANDOVER = "geometric direct minimisation steps from the next one"
   !! The note the handover iteration logs; nothing but GDM prints it.
   character(len=*), parameter :: BASE = "water.xyz -b basis_sets/sto-3g.json"

contains

   subroutine collect_cli_gdm(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("--gdm hands over from DIIS; a run without it does not", test_flag), &
                  new_unittest("--gdm-start below convergence never hands over", test_flag_start), &
                  new_unittest("--gdm-max-rotation 1e-4 costs iterations", test_flag_rotation), &
                  new_unittest("[scf.convergence] gdm = true hands over", test_toml_on), &
                  new_unittest("[scf.convergence] gdm_start below convergence never hands over", test_toml_start), &
                  new_unittest("[scf.convergence] gdm_max_rotation 1e-4 costs iterations", test_toml_rotation), &
                  new_unittest("bad --gdm values and --gdm with --soscf are refused", test_refusals) &
                  ]
   end subroutine collect_cli_gdm

   subroutine run_cmd(cmd, text, rc)
      !! The merged stdout and stderr of a command, and its exit code.
      character(len=*), intent(in) :: cmd
      character(len=:), allocatable, intent(out) :: text
      integer, intent(out) :: rc
      character(len=1024) :: line
      character(len=:), allocatable :: out
      integer :: u, ios
      out = work_dir//"/cli_gdm.out"
      ! A stale file must never be read as this command's output.
      call execute_command_line("rm -f "//out, wait=.true.)
      call execute_command_line(cmd//" > "//out//" 2>&1", exitstat=rc, wait=.true.)
      text = ""
      open (newunit=u, file=out, status="old", action="read", iostat=ios)
      if (ios /= 0) return
      do
         read (u, "(a)", iostat=ios) line
         if (ios /= 0) exit
         text = text//trim(line)//new_line("a")
      end do
      close (u)
   end subroutine run_cmd

   function cli(args) result(text)
      character(len=*), intent(in) :: args
      character(len=:), allocatable :: text
      integer :: rc
      call run_cmd(cli_exe//" "//args, text, rc)
   end function cli

   function toml_run(lines) result(text)
      !! Run the CLI on a TOML file with the given [scf.convergence] lines.
      character(len=*), intent(in) :: lines(:)
      character(len=:), allocatable :: text
      integer :: u, i, rc
      call execute_command_line("mkdir -p "//work_dir//"/case/basis_sets", wait=.true.)
      call execute_command_line("cp -f water.xyz "//work_dir//"/case/water.xyz", wait=.true.)
      call execute_command_line("cp -f basis_sets/sto-3g.json "//work_dir//"/case/basis_sets/sto-3g.json", wait=.true.)
      open (newunit=u, file=work_dir//"/case/gdm.toml", status="replace", action="write")
      write (u, "(a)") "[molecule]", 'geometry = "water.xyz"', "[basis]", 'main = "basis_sets/sto-3g.json"'
      write (u, "(a)") "[scf.convergence]"
      do i = 1, size(lines)
         write (u, "(a)") trim(lines(i))
      end do
      close (u)
      call run_cmd("cd "//work_dir//" && "//cli_exe//" case/gdm.toml", text, rc)
   end function toml_run

   function iterations(text) result(n)
      !! The iteration count of the line "SCF : N iterations, converged"; -1 if absent.
      character(len=*), intent(in) :: text
      integer :: n
      integer :: k, j
      n = -1
      k = index(text, " iterations, converged")
      if (k == 0) return
      j = k - 1
      do while (j > 0)
         if (text(j:j) < "0" .or. text(j:j) > "9") exit
         j = j - 1
      end do
      read (text(j + 1:k - 1), *) n
   end function iterations

   subroutine test_flag(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: plain, gdm
      plain = cli(BASE)
      gdm = cli(BASE//" --gdm")
      call check(error, iterations(plain) > 0 .and. index(plain, HANDOVER) == 0, &
                 "the plain run must converge and print no GDM note:"//new_line("a")//plain)
      if (allocated(error)) return
      call check(error, iterations(gdm) > 0, "the --gdm run must converge:"//new_line("a")//gdm)
      if (allocated(error)) return
      call check(error, index(gdm, HANDOVER) > 0, "--gdm did not reach the SCF: no handover note:"//new_line("a")//gdm)
   end subroutine test_flag

   subroutine test_flag_start(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: text
      text = cli(BASE//" --gdm --gdm-start 1e-12")
      call check(error, iterations(text) > 0, "the run must converge:"//new_line("a")//text)
      if (allocated(error)) return
      call check(error, index(text, HANDOVER) == 0, "--gdm-start 1e-12 handed over anyway:"//new_line("a")//text)
   end subroutine test_flag_start

   subroutine test_flag_rotation(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: dflt, capped
      dflt = cli(BASE//" --gdm")
      capped = cli(BASE//" --gdm --gdm-max-rotation 1e-4")
      call check(error, iterations(dflt) > 0 .and. iterations(capped) > 0, "both runs must converge:"// &
                 new_line("a")//capped)
      if (allocated(error)) return
      call check(error, iterations(capped) > iterations(dflt), "--gdm-max-rotation 1e-4 took "// &
                 to_string(iterations(capped))//" iterations against "//to_string(iterations(dflt))// &
                 " with the default: the cap was not honoured")
   end subroutine test_flag_rotation

   subroutine test_toml_on(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: text
      text = toml_run([character(len=64) :: "gdm = true"])
      call check(error, iterations(text) > 0 .and. index(text, HANDOVER) > 0, &
                 "gdm = true did not reach the SCF:"//new_line("a")//text)
   end subroutine test_toml_on

   subroutine test_toml_start(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: text
      text = toml_run([character(len=64) :: "gdm = true", "gdm_start = 1e-12"])
      call check(error, iterations(text) > 0 .and. index(text, HANDOVER) == 0, &
                 "gdm_start = 1e-12 handed over anyway (or the run failed):"//new_line("a")//text)
   end subroutine test_toml_start

   subroutine test_toml_rotation(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: dflt, capped
      dflt = toml_run([character(len=64) :: "gdm = true"])
      capped = toml_run([character(len=64) :: "gdm = true", "gdm_max_rotation = 1e-4"])
      call check(error, iterations(dflt) > 0 .and. iterations(capped) > iterations(dflt), &
                 "gdm_max_rotation = 1e-4 took "//to_string(iterations(capped))//" iterations against "// &
                 to_string(iterations(dflt))//" with the default:"//new_line("a")//capped)
   end subroutine test_toml_rotation

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: text
      integer :: rc
      call run_cmd(cli_exe//" "//BASE//" --gdm --soscf", text, rc)
      call check(error, rc /= 0 .and. index(text, "refuses --soscf") > 0, &
                 "--gdm with --soscf must be refused by name:"//new_line("a")//text)
      if (allocated(error)) return
      call run_cmd(cli_exe//" "//BASE//" --gdm --gdm-start 0", text, rc)
      call check(error, rc /= 0 .and. index(text, "--gdm-start must be > 0") > 0, &
                 "--gdm-start 0 must be refused:"//new_line("a")//text)
      if (allocated(error)) return
      call run_cmd(cli_exe//" "//BASE//" --gdm --gdm-max-rotation -1", text, rc)
      call check(error, rc /= 0 .and. index(text, "--gdm-max-rotation must be > 0") > 0, &
                 "--gdm-max-rotation -1 must be refused:"//new_line("a")//text)
   end subroutine test_refusals

end module check_cli_gdm_suite

program check_cli_gdm
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_gdm_suite, only: collect_cli_gdm, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_cli_gdm: usage: check_cli_gdm <terco CLI> <scratch dir>"
      error stop 1
   end if
   call execute_command_line("mkdir -p "//work_dir, wait=.true.)
   stat = 0
   call run_testsuite(collect_cli_gdm, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_gdm
