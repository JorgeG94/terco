module check_cli_ecp_suite
   !!
   !! Effective core potentials through the `terco` CLI.
   !!
   !!   - HI in def2-SVP attaches the basis's own ECP: the energy is the pyscf
   !!     anchor E_hi (ecp_scf.ref), and the electron count the run reports (and
   !!     takes the default --uhf from) is the valence one, 26, not 54.
   !!   - --ecp names another ECP (basis_sets/ecp-i-modified.json, a test
   !!     fixture: iodine's def2-SVP ECP with one coefficient changed), which
   !!     replaces the basis's own; a TOML file's [basis] ecp gives the same
   !!     numbers as the flag.
   !!   - RbI with that iodine-only --ecp is refused: rubidium's def2-SVP
   !!     shells are valence-only and nothing would replace its core.
   !!
   !! A test-drive suite.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_real
   implicit none
   private
   public :: collect_cli_ecp, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns; the runs happen in `case/`.
   real(dp), parameter :: ANCHOR_TOL = 1.0e-8_dp
   !! The CLI's own run of check_ecp_scf's HI, at the same convergence
   !! (1e-11 Eh, 1e-8 commutator).
   real(dp), parameter :: SAME_TOL = 1.0e-10_dp
   !! A TOML file and the flags it stands for reach the same code
   !! (check_cli_input's own reasoning and tolerance).
   character(len=*), parameter :: CONV = ' --econv 1e-11 --gconv 1e-8 --iterations 150'

contains

   subroutine collect_cli_ecp(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("HI/def2-SVP: the anchor energy, and 26 electrons reported", test_auto_attach), &
                  new_unittest("--ecp replaces the basis's own ECP", test_ecp_flag), &
                  new_unittest("a TOML [basis] ecp gives the --ecp numbers", test_ecp_toml), &
                  new_unittest("RbI with an iodine-only --ecp is refused", test_partial_ecp)]
   end subroutine collect_cli_ecp

   subroutine test_auto_attach(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e
      integer :: rc
      character(len=:), allocatable :: text

      call setup()
      call run('cd '//work_dir//'/case && '//cli_exe//' hi.xyz -b basis_sets/def2-svp.json'//CONV// &
               ' --json auto.json > auto.out 2>&1', rc)
      if (rc /= 0) then
         call test_failed(error, "the HI run failed, exit code "//to_string(rc)); return
      end if
      text = slurp(work_dir//'/case/auto.out')
      call check(error, index(text, "electrons: 26") > 0, &
                 "the run should report HI's 26 valence electrons (def2-SVP's ECP takes 28 of iodine's 53):"// &
                 new_line('a')//text)
      if (allocated(error)) return
      call read_energy(work_dir//'/case/auto.json', e, error)
      if (allocated(error)) return
      call check(error, abs(e - anchor_e_hi()) < ANCHOR_TOL, &
                 "|E - E_hi(pyscf)| = "//rtoa(abs(e - anchor_e_hi())))
   end subroutine test_auto_attach

   subroutine test_ecp_flag(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e

      call flag_run(e, error)
      if (allocated(error)) return
      ! The changed coefficient moves the energy by about 2e-2 Eh.
      call check(error, abs(e - anchor_e_hi()) > 1.0e-4_dp, &
                 "--ecp should replace the basis's own ECP: E = "//rtoa(e)//", the def2-SVP anchor is "// &
                 rtoa(anchor_e_hi()))
   end subroutine test_ecp_flag

   subroutine test_ecp_toml(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e_flag, e_toml
      integer :: rc

      call flag_run(e_flag, error)
      if (allocated(error)) return
      call write_lines(work_dir//'/case/ecp.toml', [character(len=64) :: &
                       '[molecule]', 'geometry = "hi.xyz"', &
                       '[basis]', 'main = "basis_sets/def2-svp.json"', 'ecp = "basis_sets/ecp-i-modified.json"', &
                       '[scf]', 'max_iterations = 150', 'energy_tol = 1e-11', 'commutator_tol = 1e-8', &
                       '[output]', 'json = "toml.json"', 'silent = true'])
      call run('cd '//work_dir//' && '//cli_exe//' case/ecp.toml > toml.out 2>&1', rc)
      if (rc /= 0) then
         call test_failed(error, "the TOML run failed, exit code "//to_string(rc)); return
      end if
      call read_energy(work_dir//'/case/toml.json', e_toml, error)
      if (allocated(error)) return
      call check(error, abs(e_toml - e_flag) < SAME_TOL, &
                 "[basis] ecp: E = "//rtoa(e_toml)//", --ecp: "//rtoa(e_flag))
   end subroutine test_ecp_toml

   subroutine test_partial_ecp(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: rc
      character(len=:), allocatable :: text

      call setup()
      call run('cd '//work_dir//'/case && '//cli_exe//' rbi.xyz -b basis_sets/def2-svp.json '// &
               '--ecp basis_sets/ecp-i-modified.json --silent > rbi.out 2>&1', rc)
      text = slurp(work_dir//'/case/rbi.out')
      call check(error, rc /= 0, "RbI with an iodine-only ECP ran (exit code 0)")
      if (allocated(error)) return
      call check(error, index(text, "valence-only") > 0, "the refusal should name the valence-only basis:"// &
                 new_line('a')//text)
      if (allocated(error)) return
      call check(error, index(text, "TOTAL ENERGY") == 0, "an SCF ran for a refused input")
   end subroutine test_partial_ecp

   subroutine flag_run(e, error)
      !! ---- helpers ---------------------------------------------------------
      !! HI/def2-SVP with --ecp ecp-i-modified.json.
      real(dp), intent(out) :: e
      type(error_type), allocatable, intent(out) :: error
      integer :: rc
      e = 0.0_dp
      call setup()
      call run('cd '//work_dir//'/case && '//cli_exe//' hi.xyz -b basis_sets/def2-svp.json '// &
               '--ecp basis_sets/ecp-i-modified.json'//CONV//' --json flag.json --silent > flag.out 2>&1', rc)
      if (rc /= 0) then
         call test_failed(error, "the --ecp run failed, exit code "//to_string(rc)); return
      end if
      call read_energy(work_dir//'/case/flag.json', e, error)
   end subroutine flag_run

   subroutine setup()
      !! case/ with the geometries and the basis/ECP files the runs name.
      integer :: rc
      call run('mkdir -p '//work_dir//'/case/basis_sets', rc)
      call run('cp -f basis_sets/def2-svp.json basis_sets/ecp-i-modified.json '//work_dir//'/case/basis_sets/', rc)
      ! check_ecp_scf's HI: H at the origin, I 1.609 Angstrom along z.
      call write_lines(work_dir//'/case/hi.xyz', [character(len=64) :: '2', 'HI', &
                       'H   0.0   0.0   0.0', 'I   0.0   0.0   1.609'])
      call write_lines(work_dir//'/case/rbi.xyz', [character(len=64) :: '2', 'RbI', &
                       'Rb  0.0   0.0   0.0', 'I   0.0   0.0   3.18'])
   end subroutine setup

   function anchor_e_hi() result(e)
      !! E_hi, RHF HI/def2-SVP from pyscf (ecp_scf.ref's first value).
      real(dp) :: e
      integer :: u
      logical :: ok
      call ref_open('ecp_scf.ref', u, ok)
      if (.not. ok) error stop 'check_cli_ecp: ecp_scf.ref missing -- run test/pyscf_refs/generate.py ecp_scf'
      call ref_read_real(u, 'E_hi', e)
      call ref_close(u)
   end function anchor_e_hi

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
         call json%destroy(); call test_failed(error, "cannot read "//path); return
      end if
      call json%get('total energy', energy, found)
      call json%destroy()
      if (.not. found) call test_failed(error, "no total energy in "//path)
   end subroutine read_energy

   subroutine write_lines(path, lines)
      character(len=*), intent(in) :: path, lines(:)
      integer :: u, i
      open (newunit=u, file=path, status='replace', action='write')
      do i = 1, size(lines)
         write (u, '(a)') trim(lines(i))
      end do
      close (u)
   end subroutine write_lines

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
      real(dp), intent(in) :: x
      character(len=:), allocatable :: s
      s = trim(adjustl(to_string(x, "es24.16")))
   end function rtoa

end module check_cli_ecp_suite

program check_cli_ecp
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_ecp_suite, only: collect_cli_ecp, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_cli_ecp: usage: check_cli_ecp <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_cli_ecp, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_ecp
