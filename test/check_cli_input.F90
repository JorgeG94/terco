!
! The CLI's TOML input files (docs/INPUT.md, app/trc_run_input.F90), through
! the real `terco` binary.
!
! The reference for a file is the command line it stands for: the same run
! given as flags must give the same numbers, since both reach the same cli_t
! and the same code from there. On the host build they are bit-identical; the
! GPU build is not bitwise reproducible run to run (two identical flag runs on
! a V100, nvfortran 26.5: energy 4e-14 Eh, gradient 5e-15 Eh/Bohr apart), so
! the comparison allows ENERGY_TOL/GRADIENT_TOL -- a thousand times that noise,
! and far below what any setting this run changes would move. The file sits in a subdirectory and the
! run starts from its parent, so the geometry, the basis file and the JSON
! output only resolve if paths are taken relative to the file, as the format
! says. Then: the "input" echo in --json shows the file's values, a flag
! overriding one, and a declared default; and every way a file can be wrong
! is refused with a message naming the problem, before any SCF runs.
!
! A test-drive suite.
!
module check_cli_input_suite
   use pic_types, only: dp
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   implicit none
   private
   public :: collect_cli_input, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns; the input files go in its
      !! `case/` subdirectory.
   real(dp), parameter :: ENERGY_TOL = 1.0e-11_dp
   real(dp), parameter :: GRADIENT_TOL = 1.0e-10_dp

contains

   subroutine collect_cli_input(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("a TOML file gives the same numbers as its flags", test_same_as_flags), &
         new_unittest("--json echoes the resolved input, with flags overriding the file", test_echo), &
         new_unittest("an unknown key is refused, by name", test_unknown_key), &
         new_unittest("an unknown section is refused, by name", test_unknown_section), &
         new_unittest("a value of the wrong type is refused", test_wrong_type), &
         new_unittest("a key outside any section is refused", test_outside_section), &
         new_unittest("a file without [molecule] geometry is refused", test_no_geometry), &
         new_unittest("an unknown [correlation] method is refused", test_bad_method), &
         new_unittest("a scalar where an array belongs is refused", test_not_array), &
         new_unittest("a TOML syntax error is refused", test_syntax_error)]
   end subroutine collect_cli_input

   subroutine test_same_as_flags(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e_toml, e_flags
      real(dp), allocatable :: g_toml(:), g_flags(:)
      integer :: rc

      call setup()
      call write_toml(work_dir//'/case/run.toml', [character(len=64) :: &
                      '[molecule]', 'geometry = "water.xyz"', 'charge = 0', &
                      '[basis]', 'main = "basis_sets/sto-3g.json"', &
                      '[scf]', 'max_iterations = 60', 'energy_tol = 1e-11', 'commutator_tol = 1e-7', &
                      '[driver]', 'gradient = true', &
                      '[output]', 'json = "toml.json"', 'silent = true'])
      ! From the PARENT of case/: nothing resolves unless paths follow the file.
      call run('cd '//work_dir//' && '//cli_exe//' case/run.toml > toml.out 2>&1', rc)
      if (rc /= 0) then
         call test_failed(error, "the TOML run failed, exit code "//itoa(rc)); return
      end if
      call run('cd '//work_dir//'/case && '//cli_exe//' water.xyz -b basis_sets/sto-3g.json --iterations 60 '// &
               '--econv 1e-11 --gconv 1e-7 --grad --json flags.json --silent > ../flags.out 2>&1', rc)
      if (rc /= 0) then
         call test_failed(error, "the flag run failed, exit code "//itoa(rc)); return
      end if

      call read_result(work_dir//'/case/toml.json', e_toml, g_toml, error)
      if (allocated(error)) return
      call read_result(work_dir//'/case/flags.json', e_flags, g_flags, error)
      if (allocated(error)) return
      call check(error, abs(e_toml - e_flags) < ENERGY_TOL, &
                 "total energy differs: file "//rtoa(e_toml)//", flags "//rtoa(e_flags))
      if (allocated(error)) return
      call check(error, size(g_toml) == size(g_flags) .and. size(g_toml) == 9, "gradient sizes differ")
      if (allocated(error)) return
      call check(error, maxval(abs(g_toml - g_flags)) < GRADIENT_TOL, &
                 "gradient differs between the file and the flags by "//rtoa(maxval(abs(g_toml - g_flags))))
   end subroutine test_same_as_flags

   subroutine test_echo(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      integer :: rc, ival
      real(dp) :: rval
      logical :: found, lval
      character(len=:), allocatable :: sval

      call setup()
      call write_toml(work_dir//'/case/echo.toml', [character(len=64) :: &
                      '[molecule]', 'geometry = "water.xyz"', &
                      '[basis]', 'main = "basis_sets/sto-3g.json"', &
                      '[scf]', 'energy_tol = 1e-9', &
                      '[dft]', 'grid = 2', &
                      '[output]', 'json = "echo.json"', 'silent = true'])
      call run('cd '//work_dir//' && '//cli_exe//' case/echo.toml --grid 4 > echo.out 2>&1', rc)
      if (rc /= 0) then
         call test_failed(error, "the run failed, exit code "//itoa(rc)); return
      end if
      call json%initialize()
      call json%load(filename=work_dir//'/case/echo.json')
      call json%get('input.scf.energy_tol', rval, found)
      call check(error, found .and. rval == 1.0e-9_dp, "input.scf.energy_tol should be the file's 1e-9")
      if (allocated(error)) return
      call json%get('input.dft.grid', ival, found)
      call check(error, found .and. ival == 4, "input.dft.grid should be the flag's 4, not the file's 2")
      if (allocated(error)) return
      call json%get('input.scf.max_iterations', ival, found)
      call check(error, found .and. ival == 100, "input.scf.max_iterations should be the declared default 100")
      if (allocated(error)) return
      call json%get('input.driver.gradient', lval, found)
      call check(error, found .and. .not. lval, "input.driver.gradient should be false")
      if (allocated(error)) return
      call json%get('input.molecule.geometry', sval, found)
      call check(error, found .and. sval == 'case/water.xyz', "input.molecule.geometry should be case/water.xyz")
      if (allocated(error)) return
      call json%get('input.molecule.charge', sval, found)
      call check(error, .not. found, "input.molecule.charge was never given and has no declared default")
      call json%destroy()
   end subroutine test_echo

   subroutine test_unknown_key(error)
      type(error_type), allocatable, intent(out) :: error
      call refused([character(len=64) :: '[molecule]', 'geometry = "water.xyz"', &
                    '[scf]', 'energy_tolerance = 1e-8'], "unknown key 'energy_tolerance'", error)
   end subroutine test_unknown_key

   subroutine test_unknown_section(error)
      type(error_type), allocatable, intent(out) :: error
      call refused([character(len=64) :: '[molecule]', 'geometry = "water.xyz"', &
                    '[solvent]', 'model = "cpcm"'], "unknown section [solvent]", error)
   end subroutine test_unknown_section

   subroutine test_wrong_type(error)
      type(error_type), allocatable, intent(out) :: error
      call refused([character(len=64) :: '[molecule]', 'geometry = "water.xyz"', 'charge = "0"'], &
                   "'charge' must be an integer", error)
   end subroutine test_wrong_type

   subroutine test_outside_section(error)
      type(error_type), allocatable, intent(out) :: error
      call refused([character(len=64) :: 'functional = "hf"', '[molecule]', 'geometry = "water.xyz"'], &
                   "outside any section", error)
   end subroutine test_outside_section

   subroutine test_no_geometry(error)
      type(error_type), allocatable, intent(out) :: error
      call refused([character(len=64) :: '[scf]', 'functional = "hf"'], "[molecule] geometry is required", error)
   end subroutine test_no_geometry

   subroutine test_bad_method(error)
      type(error_type), allocatable, intent(out) :: error
      call refused([character(len=64) :: '[molecule]', 'geometry = "water.xyz"', &
                    '[correlation]', 'method = "ccsd"'], "method must be 'rimp2' or 'none'", error)
   end subroutine test_bad_method

   subroutine test_not_array(error)
      type(error_type), allocatable, intent(out) :: error
      call refused([character(len=64) :: '[molecule]', 'geometry = "water.xyz"', &
                    '[opt]', 'frozen_atoms = 3'], "'frozen_atoms' must be an array", error)
   end subroutine test_not_array

   subroutine test_syntax_error(error)
      type(error_type), allocatable, intent(out) :: error
      call refused([character(len=64) :: '[molecule', 'geometry = "water.xyz"'], "", error)
   end subroutine test_syntax_error

   subroutine refused(lines, expect, error)
      !! A file that must stop the run before any SCF: nonzero exit, `expect`
      !! in the output (any message if `expect` is blank), and no SCF line.
      character(len=*), intent(in) :: lines(:), expect
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: text
      integer :: rc

      call setup()
      call write_toml(work_dir//'/case/bad.toml', lines)
      call run('cd '//work_dir//' && '//cli_exe//' case/bad.toml > bad.out 2>&1', rc)
      call check(error, rc /= 0, "the file was accepted (exit code 0)")
      if (allocated(error)) return
      text = slurp(work_dir//'/bad.out')
      if (len_trim(expect) > 0) then
         call check(error, index(text, expect) > 0, "the message does not say '"//expect//"':"// &
                    new_line('a')//text)
         if (allocated(error)) return
      end if
      call check(error, index(text, "TOTAL ENERGY") == 0, "an SCF ran for a file that should have been refused")
   end subroutine refused

   subroutine setup()
      !! case/ with the geometry and the basis file the inputs name.
      integer :: rc
      call run('mkdir -p '//work_dir//'/case/basis_sets', rc)
      call run('cp -f water.xyz '//work_dir//'/case/water.xyz', rc)
      call run('cp -f basis_sets/sto-3g.json '//work_dir//'/case/basis_sets/sto-3g.json', rc)
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

   subroutine read_result(path, energy, gradient, error)
      character(len=*), intent(in) :: path
      real(dp), intent(out) :: energy
      real(dp), allocatable, intent(out) :: gradient(:)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      logical :: found
      real(dp), allocatable :: row(:)
      integer :: ia
      call json%initialize()
      call json%load(filename=path)
      if (json%failed()) then
         call test_failed(error, "cannot read "//path); return
      end if
      call json%get('total energy', energy, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no total energy in "//path); return
      end if
      ! One row per atom; water has three.
      allocate (gradient(0))
      do ia = 1, 3
         call json%get('gradient('//itoa(ia)//')', row, found)
         if (.not. found) then
            call json%destroy(); call test_failed(error, "no gradient row "//itoa(ia)//" in "//path); return
         end if
         gradient = [gradient, row]
      end do
      call json%destroy()
   end subroutine read_result

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

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=16) :: buf
      write (buf, '(i0)') i
      s = trim(buf)
   end function itoa

   pure function rtoa(x) result(s)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, '(es24.16)') x
      s = trim(adjustl(buf))
   end function rtoa

end module check_cli_input_suite

program check_cli_input
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_input_suite, only: collect_cli_input, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_cli_input: usage: check_cli_input <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_cli_input, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_input
