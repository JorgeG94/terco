!
! The CLI's end-of-SCF frontier orbital energies, read back from its
! --json: "homo energy"/"lumo energy" (restricted, or alpha) and, for an
! unrestricted run, "beta homo energy"/"beta lumo energy".
!
! References: pyscf 2.x, the same water geometry (test/water.xyz) in
! 6-31G with cart=True, conv_tol 1e-12. pyscf's built-in 6-31G and the
! Basis Set Exchange JSON terco reads agree on water's RHF energy to 1e-8
! (-75.98397447 both), so 1e-6 Eh on the orbital energies is a check of
! the reporting, not of the basis.
!
! A test-drive suite.
!
module check_cli_homo_lumo_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   implicit none
   private
   public :: collect_cli_homo_lumo, cli_exe, probe_json

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: probe_json
      !! Where the CLI writes its --json, from the command line.

   real(dp), parameter :: TOL = 1.0e-6_dp
      !! Eh, against pyscf; see the header.

contains

   subroutine collect_cli_homo_lumo(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water RHF/6-31G: HOMO and LUMO in the JSON", test_rhf), &
                  new_unittest("water+ UHF/6-31G: alpha and beta HOMO/LUMO in the JSON", test_uhf) &
                  ]
   end subroutine collect_cli_homo_lumo

   subroutine run_cli(error, extra, json)
      !! One water/6-31G run with `extra` options, its --json loaded.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: extra
      type(json_file), intent(inout) :: json
      integer :: rc
      call execute_command_line(trim(cli_exe)//' water.xyz -b basis_sets/6-31g.json --econv 1e-11 '// &
                                extra//' --silent --json '//probe_json//' > '//probe_json//'.out 2>&1', &
                                exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run failed, exit code "//to_string(rc)); return
      end if
      call json%initialize()
      call json%load(filename=probe_json)
      if (json%failed()) call test_failed(error, "cannot read "//probe_json)
   end subroutine run_cli

   subroutine expect(error, json, key, ref)
      type(error_type), allocatable, intent(out) :: error
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: key
      real(dp), intent(in) :: ref
      real(dp) :: val
      logical :: found
      call json%get(key, val, found)
      if (.not. found) then
         call test_failed(error, "no '"//key//"' in "//probe_json); return
      end if
      call check(error, abs(val - ref) <= TOL, key//" = "//to_string(val)//", pyscf "//to_string(ref))
   end subroutine expect

   subroutine test_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      logical :: found
      real(dp) :: val
      call run_cli(error, '', json)
      if (allocated(error)) return
      call expect(error, json, 'homo energy', -0.501368125579746_dp)
      if (allocated(error)) then
         call json%destroy(); return
      end if
      call expect(error, json, 'lumo energy', 0.20364089470362523_dp)
      if (allocated(error)) then
         call json%destroy(); return
      end if
      ! A restricted run has no separate beta set.
      call json%get('beta homo energy', val, found)
      call json%destroy()
      call check(error, .not. found, "a restricted run wrote 'beta homo energy'")
   end subroutine test_rhf

   subroutine test_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      call run_cli(error, '--chrg 1 --uhf 1', json)
      if (allocated(error)) return
      call expect(error, json, 'homo energy', -1.093552467485519_dp)
      if (.not. allocated(error)) call expect(error, json, 'lumo energy', -0.13849112838386746_dp)
      if (.not. allocated(error)) call expect(error, json, 'beta homo energy', -1.0314319118318325_dp)
      if (.not. allocated(error)) call expect(error, json, 'beta lumo energy', -0.3226565429822299_dp)
      call json%destroy()
   end subroutine test_uhf

end module check_cli_homo_lumo_suite

program check_cli_homo_lumo
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_homo_lumo_suite, only: collect_cli_homo_lumo, cli_exe, probe_json
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); probe_json = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0) then
      write (error_unit, "(a)") "check_cli_homo_lumo: usage: check_cli_homo_lumo <terco CLI> <scratch .json>"
      error stop 1
   end if
   stat = 0
   ! Serial: both cases write the same probe file.
   call run_testsuite(collect_cli_homo_lumo, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_homo_lumo
