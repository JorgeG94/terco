module check_cli_incremental_suite
   !!
   !! The CLI's switches for the incremental Fock build, which is on by
   !! default: the default run builds increments (iteration lines marked `i`),
   !! while `--no-incremental` and a TOML file's `incremental = false` each turn
   !! them off -- no `i` line, and the energy of the full-build SCF. The two
   !! off switches must agree to run-to-run noise (the same full builds), and
   !! with the default to the convergence test's scale. `--incremental` together with
   !! `--no-incremental`, from flags or from the file and a flag, is refused
   !! before any SCF runs.
   !!
   !! A test-drive suite, in check_cli_input.F90's style: the CLI binary and
   !! a scratch directory come from the command line, the run happens in test/.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   implicit none
   private
   public :: collect_cli_incremental, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! Scratch directory for the TOML file, the logs and the JSON.

   character(len=*), parameter :: SYSTEM = "w3.xyz -b basis_sets/6-31gs.json"
   real(dp), parameter :: SAME_TOL = 1.0e-11_dp
      !! The two off switches run the same full builds; the OpenMP port's
      !! atomics sum in a different order every run, 6e-14 Eh measured.
   real(dp), parameter :: ETOL = 1.0e-10_dp
      !! The default (incremental) run against full builds: the CLI's own
      !! --econv 1e-10; measured 1.1e-13 (gfortran 15), 5.7e-14 (ifx 2026.2
      !! and the OpenMP port).

contains

   subroutine collect_cli_incremental(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("default on, --no-incremental and TOML false off", test_off_switches), &
                  new_unittest("--incremental with --no-incremental is refused", test_contradiction_flags), &
                  new_unittest("TOML false with --incremental is refused", test_contradiction_toml) &
                  ]
   end subroutine collect_cli_incremental

   subroutine test_off_switches(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e_on, e_flag, e_toml
      integer :: n_on, n_flag, n_toml, rc
      call setup()
      call run_case("on", "", e_on, n_on, error)
      if (allocated(error)) return
      call run_case("flag", "--no-incremental", e_flag, n_flag, error)
      if (allocated(error)) return
      call write_toml(work_dir//"/off.toml", false_input())
      call run_case("toml", work_dir//"/off.toml", e_toml, n_toml, error)
      if (allocated(error)) return
      call check(error, n_on > 0, "the default run built no increments")
      if (allocated(error)) return
      call check(error, n_flag == 0, "--no-incremental still built "//to_string(n_flag)//" increments")
      if (allocated(error)) return
      call check(error, n_toml == 0, "incremental = false still built "//to_string(n_toml)//" increments")
      if (allocated(error)) return
      call check(error, abs(e_flag - e_toml) <= SAME_TOL, "--no-incremental and incremental = false differ: "// &
                 rtoa(e_flag)//" against "//rtoa(e_toml))
      if (allocated(error)) return
      call check(error, abs(e_on - e_flag) <= ETOL, "incremental "//rtoa(e_on)//" against full builds "// &
                 rtoa(e_flag))
   end subroutine test_off_switches

   subroutine test_contradiction_flags(error)
      type(error_type), allocatable, intent(out) :: error
      call refused("both_flags", trim(cli_exe)//" "//SYSTEM//" --incremental --no-incremental", error)
   end subroutine test_contradiction_flags

   subroutine test_contradiction_toml(error)
      type(error_type), allocatable, intent(out) :: error
      call setup()
      call write_toml(work_dir//"/off2.toml", false_input())
      ! The file's switch comes first and the command line after it, so the
      ! flag does not silently win: the two are refused together.
      call refused("toml_flag", trim(cli_exe)//" "//work_dir//"/off2.toml --incremental", error)
   end subroutine test_contradiction_toml

   subroutine run_case(tag, extra, energy, n_incr, error)
      !! One CLI run on the water trimer: its total energy (JSON) and the
      !! number of iteration lines marked `i` in its log. `extra` is a
      !! flag, or a TOML file, which goes first (a file must be the first
      !! argument).
      character(len=*), intent(in) :: tag, extra
      real(dp), intent(out) :: energy
      integer, intent(out) :: n_incr
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, log_path, json_path
      integer :: rc
      log_path = work_dir//"/"//tag//".log"
      json_path = work_dir//"/"//tag//".json"
      if (index(extra, ".toml") > 0) then
         cmd = trim(cli_exe)//" "//extra
      else
         cmd = trim(cli_exe)//" "//SYSTEM//" "//extra
      end if
      cmd = cmd//" --json "//json_path//" > "//log_path//" 2>&1"
      call run(cmd, rc)
      energy = 0.0_dp
      n_incr = 0
      if (rc /= 0) then
         call test_failed(error, tag//": the CLI run failed, exit code "//to_string(rc)//" ("//log_path//")")
         return
      end if
      call read_energy(json_path, energy, error)
      if (allocated(error)) return
      n_incr = count_increments(log_path)
   end subroutine run_case

   subroutine refused(tag, cmd, error)
      !! `cmd` must exit nonzero, saying why, without running an SCF.
      character(len=*), intent(in) :: tag, cmd
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: log_path, text
      integer :: rc
      call setup()
      log_path = work_dir//"/"//tag//".log"
      call run(cmd//" > "//log_path//" 2>&1", rc)
      text = slurp(log_path)
      call check(error, rc /= 0, tag//": contradictory incremental switches were accepted")
      if (allocated(error)) return
      call check(error, index(text, "contradict") > 0, tag//": the refusal does not say why")
      if (allocated(error)) return
      call check(error, index(text, "TOTAL ENERGY") == 0, tag//": an SCF ran before the refusal")
   end subroutine refused

   integer function count_increments(path) result(n)
      !! Iteration lines marked as built from the density difference.
      character(len=*), intent(in) :: path
      character(len=1024) :: line
      integer :: u, ios, l
      n = 0
      open (newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         l = len_trim(line)
         if (l > 2) then
            if (line(l - 1:l) == " i") n = n + 1
         end if
      end do
      close (u)
   end function count_increments

   subroutine read_energy(path, energy, error)
      character(len=*), intent(in) :: path
      real(dp), intent(out) :: energy
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      logical :: found
      energy = 0.0_dp
      call json%initialize()
      call json%load(filename=path)
      if (json%failed()) then
         call test_failed(error, "cannot read "//path); return
      end if
      call json%get('total energy', energy, found)
      call json%destroy()
      if (.not. found) call test_failed(error, "no total energy in "//path)
   end subroutine read_energy

   function false_input() result(lines)
      !! The same system as SYSTEM, as a TOML file next to its copied inputs.
      character(len=200), allocatable :: lines(:)
      lines = [character(len=200) :: "[molecule]", 'geometry = "w3.xyz"', "[basis]", 'main = "6-31gs.json"', &
               "[scf]", "incremental = false"]
   end function false_input

   subroutine setup()
      !! The scratch directory, with the geometry and basis the TOML file names.
      integer :: rc
      call run('mkdir -p '//work_dir, rc)
      call run('cp -f w3.xyz '//work_dir//'/w3.xyz', rc)
      call run('cp -f basis_sets/6-31gs.json '//work_dir//'/6-31gs.json', rc)
   end subroutine setup

   subroutine write_toml(path, lines)
      character(len=*), intent(in) :: path, lines(:)
      integer :: u, i
      open (newunit=u, file=path, status='replace', action='write')
      do i = 1, size(lines)
         write (u, '(a)') trim(lines(i))
      end do
      close (u)
   end subroutine write_toml

   subroutine run(cmd, rc)
      character(len=*), intent(in) :: cmd
      integer, intent(out) :: rc
      rc = 0
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
   end subroutine run

   function slurp(path) result(text)
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: text
      character(len=1024) :: line
      integer :: u, ios
      text = ''
      open (newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         text = text//trim(line)//new_line('a')
      end do
      close (u)
   end function slurp

   pure function rtoa(x) result(s)
      !! Enough digits to tell two energies apart in a failure message.
      real(dp), intent(in) :: x
      character(len=:), allocatable :: s
      s = trim(adjustl(to_string(x, "es24.16")))
   end function rtoa

end module check_cli_incremental_suite

program check_cli_incremental
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_incremental_suite, only: collect_cli_incremental, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_cli_incremental: usage: check_cli_incremental <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_cli_incremental, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_incremental
