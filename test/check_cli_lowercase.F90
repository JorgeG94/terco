!
! The CLI's own ASCII-lowercasing of its string options (--func, --guess,
! ...), through the installed `terco` binary. app/main.F90's trc_run
! module normalises these with a small, non-elemental, allocatable-result
! helper (to_lower_alloc) INSTEAD of pic's own to_lower (pic_ascii.f90,
! `pure elemental`, automatic-length result): nvfortran 26.5 and 26.9,
! `trc_main`'s real `-fast -O3` flags, host and GPU builds alike,
! miscompile pic's to_lower at these call sites -- the result comes back
! the right LENGTH but BLANK (all-space) content, so `--func PBE` silently
! ran a plain SCF with no exchange-correlation term, and `--guess SAD`
! left `r%guess` blank (then refused as an unknown guess, since the CLI
! validates the normalised string). See .claude/context/gpu-pitfalls.md's
! "Codegen" entry.
!
! This never reproduces on gfortran (the only compiler that runs this
! suite in CI), so it cannot assert a difference between the broken and
! fixed code paths by itself; it is the test that WOULD have caught the
! bug on an nvfortran build, and does catch a regression of the same
! shape (reverting to_lower_alloc back to a bare to_lower call, or a
! typo that skips the lowercasing) on any compiler. Deliberately passes
! mixed-case CLI values (PBE, SAD) rather than already-lowercase ones, so
! a build that merely forwards the raw string unchanged still fails this
! check even where it happens to equal its own lowercasing.
!
! A test-drive suite.
!
module check_cli_lowercase_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   implicit none
   private
   public :: collect_cli_lowercase, cli_exe, probe_json

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: probe_json
      !! Where the CLI writes its --json, from the command line.

contains

   subroutine collect_cli_lowercase(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--func PBE --guess SAD: resolved method/guess are lowercased, exc is nonzero", &
                      test_func_and_guess_lowercased) &
         ]
   end subroutine collect_cli_lowercase

   subroutine test_func_and_guess_lowercased(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe_out, method_str, guess_str
      integer :: rc
      logical :: found
      real(dp) :: exc
      type(json_file) :: json

      probe_out = probe_json//'.out'
      ! Mixed case on purpose -- see this module's own header.
      cmd = trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --func PBE --guess SAD '// &
            '--silent --json '//probe_json//' > '//probe_out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run itself failed, exit code "//to_string(rc)); return
      end if

      call json%initialize()
      call json%load(filename=probe_json)
      if (json%failed()) then
         call test_failed(error, "cannot read "//probe_json); return
      end if

      call json%get('method', method_str, found)
      call json%get('guess', guess_str, found)
      call json%get('exchange-correlation energy', exc, found)
      call json%destroy()

      call check(error, method_str, 'RKS/pbe', &
                 "--func PBE did not resolve to method 'RKS/pbe' (blank-lowercase miscompile): got '"// &
                 method_str//"'")
      if (allocated(error)) return
      call check(error, guess_str, 'sad', &
                 "--guess SAD did not resolve to guess 'sad' (blank-lowercase miscompile): got '"// &
                 guess_str//"'")
      if (allocated(error)) return
      ! A plain HF/no-XC fallback (the original symptom) reports exactly
      ! 0.0 here; PBE/water/sto-3g's own XC energy is a few Eh in magnitude.
      if (abs(exc) < 1.0e-6_dp) then
         call test_failed(error, "exchange-correlation energy is "//to_string(exc)// &
                           " -- PBE did not actually run (func normalisation blanked to HF)")
      end if
   end subroutine test_func_and_guess_lowercased

end module check_cli_lowercase_suite

program check_cli_lowercase
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_lowercase_suite, only: collect_cli_lowercase, cli_exe, probe_json
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); probe_json = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0) then
      write (error_unit, "(a)") "check_cli_lowercase: usage: check_cli_lowercase <terco CLI> <scratch .json>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_cli_lowercase, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_lowercase
