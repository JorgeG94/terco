!
! The CLI's "timings" and "build" objects in --json (docs/sphinx/output.md):
! which keys a run writes for which options, that every value is finite and
! positive where it must be, that the stage times nest the way the stages
! do (total wall >= SCF wall >= the Fock builds), and that the counts agree
! with the summary's own "iterations". The values are wall times, so
! nothing here pins a number: only relations that hold on any machine.
!
! The build timings are checked against an independent source: the per-
! iteration "fock(s)" column that a non-silent run prints (the same timer read
! through the logger's own formatting, 3 decimals). What stays untestable on
! the host: that each timer stops only after the DEVICE has finished the
! build (a host build is synchronous); that is for a GPU run to show.
!
! A test-drive suite.
!
module check_cli_timings_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   implicit none
   private
   public :: collect_cli_timings, cli_exe, probe_json

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: probe_json
      !! Where the CLI writes its --json, from the command line.

   real(dp), parameter :: CLOCK_SLACK = 1.0e-3_dp
      !! s. Nested timers read the clock at different instants, so an inner
      !! total can exceed its outer one by a clock tick; a millisecond is far
      !! above any tick and far below any real stage.

contains

   subroutine collect_cli_timings(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF: scf stages only, counts match the summary", test_rhf), &
                  new_unittest("RHF: fock_first and the build count match the printed table", test_printed_table), &
                  new_unittest("--soscf: iterations exceed the Fock builds", test_soscf), &
                  new_unittest("--hessian and --fukui write total_wall only", test_total_wall_only), &
                  new_unittest("DFT with VV10: grid, xc and vv10 keys", test_dft), &
                  new_unittest("semilocal DFT: no vv10 key", test_dft_no_vv10), &
                  new_unittest("--df: dfjk keys and its steps", test_df), &
                  new_unittest("--pcm --grad: pcm and gradient keys", test_pcm_grad), &
                  new_unittest("--rimp2: rimp2 total and steps", test_rimp2), &
                  new_unittest("build object", test_build) &
                  ]
   end subroutine collect_cli_timings

   subroutine run_cli(error, extra, json, silent)
      !! One water run with `extra` options, its --json loaded. The probe file
      !! is deleted first, so a run that writes none fails here instead of
      !! reading the previous run's. `silent` false keeps the iteration table
      !! in `probe_json`.out.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: extra
      type(json_file), intent(inout) :: json
      logical, intent(in), optional :: silent
      integer :: rc, u
      logical :: quiet
      quiet = .true.
      if (present(silent)) quiet = silent
      open (newunit=u, file=probe_json, status="unknown", iostat=rc)
      if (rc == 0) close (u, status="delete")
      call execute_command_line(trim(cli_exe)//' water.xyz '//extra//merge(' --silent', '         ', quiet)//' --json '// &
                                probe_json//' > '//probe_json//'.out 2>&1', exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run failed, exit code "//to_string(rc)); return
      end if
      call json%initialize()
      call json%load(filename=probe_json)
      if (json%failed()) call test_failed(error, "the run wrote no readable "//probe_json)
   end subroutine run_cli

   subroutine printed_fock_times(error, times, n)
      !! The per-iteration "fock(s)" column of the table in `probe_json`.out
      !! (an RHF run: it, E, dE, RMS(D), |FDS-SDF|, fock, linalg, [i]).
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable, intent(out) :: times(:)
      integer, intent(out) :: n
      character(len=512) :: line
      character(len=32) :: tok(12)
      integer :: u, ios, nt, ios2, it
      logical :: in_table
      real(dp) :: v
      allocate (times(0))
      n = 0
      in_table = .false.
      open (newunit=u, file=probe_json//'.out', status="old", action="read", iostat=ios)
      if (ios /= 0) then
         call test_failed(error, "cannot read "//probe_json//".out"); return
      end if
      do
         read (u, "(a)", iostat=ios) line
         if (ios /= 0) exit
         if (index(line, "fock(s)") > 0) then
            in_table = .true.
            cycle
         end if
         if (.not. in_table) cycle
         call split_tokens(line, tok, nt)
         if (nt < 7) exit
         read (tok(1), *, iostat=ios2) it
         if (ios2 /= 0) exit
         if (trim(tok(nt)) == "i") nt = nt - 1
         read (tok(nt - 1), *, iostat=ios2) v
         if (ios2 /= 0) exit
         times = [times, v]
         n = n + 1
      end do
      close (u)
      if (n == 0) call test_failed(error, "no iteration table in "//probe_json//".out")
   end subroutine printed_fock_times

   subroutine split_tokens(line, tok, nt)
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: tok(:)
      integer, intent(out) :: nt
      integer :: i, start
      logical :: inside
      nt = 0
      inside = .false.
      start = 1
      tok = ""
      do i = 1, len_trim(line) + 1
         if (i > len_trim(line) .or. line(i:i) == ' ') then
            if (inside .and. nt < size(tok)) then
               nt = nt + 1
               tok(nt) = line(start:i - 1)
            end if
            inside = .false.
         else if (.not. inside) then
            inside = .true.
            start = i
         end if
      end do
   end subroutine split_tokens

   subroutine get_time(error, json, key, val, positive)
      !! timings.<key> must exist and be finite, and > 0 when `positive`.
      type(error_type), allocatable, intent(out) :: error
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: key
      real(dp), intent(out) :: val
      logical, intent(in) :: positive
      logical :: found
      val = 0.0_dp
      call json%get('timings.'//key, val, found)
      if (.not. found) then
         call test_failed(error, "no 'timings."//key//"' in "//probe_json); return
      end if
      if (.not. ieee_is_finite(val)) then
         call test_failed(error, "'timings."//key//"' is not finite"); return
      end if
      if (positive) then
         call check(error, val > 0.0_dp, "'timings."//key//"' = "//to_string(val)//" is not positive")
      else
         call check(error, val >= 0.0_dp, "'timings."//key//"' = "//to_string(val)//" is negative")
      end if
   end subroutine get_time

   subroutine expect_absent(error, json, key)
      type(error_type), allocatable, intent(out) :: error
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: key
      logical :: found
      found = json%valid_path('timings.'//key)
      call check(error, .not. found, "'timings."//key//"' is written but the run does not have that stage")
   end subroutine expect_absent

   subroutine test_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: total, wall, ftot, ffirst, fmean, nfock, niter
      integer :: it, nb
      logical :: found
      call run_cli(error, '-b 6-31g --econv 1e-9', json)
      if (allocated(error)) return
      call get_time(error, json, 'total_wall', total, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.wall', wall, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.setup', ftot, .false.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.fock_total', ftot, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.fock_first', ffirst, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.fock_mean_after_first', fmean, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.n_iterations', niter, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.n_fock_builds', nfock, .true.)
      if (allocated(error)) then
         call json%destroy(); return
      end if
      call json%get('iterations', it, found)
      nb = nint(nfock)
      call check(error, found .and. nint(niter) == it, "timings.scf.n_iterations does not match 'iterations'")
      if (.not. allocated(error)) call check(error, nb == it, "timings.scf.n_fock_builds "//to_string(nb)// &
                                             " /= iterations "//to_string(it))
      if (.not. allocated(error)) call check(error, total + CLOCK_SLACK >= wall, "total_wall < scf.wall")
      if (.not. allocated(error)) call check(error, wall + CLOCK_SLACK >= ftot, "scf.wall < sum of the Fock builds")
      ! Not first + (n-1) x mean == total: that holds by construction. The builds
      ! after the first are not free, and together they fit in what is left of the total.
      if (.not. allocated(error)) call check(error, fmean > 0.0_dp .and. ffirst > 0.0_dp .and. &
                                             real(nb - 1, dp)*fmean <= ftot - ffirst + 1.0e-9_dp + 1.0e-6_dp*ftot, &
                                             "the builds after the first do not fit in fock_total - fock_first")
      if (.not. allocated(error)) call expect_absent(error, json, 'dft')
      if (.not. allocated(error)) call expect_absent(error, json, 'dfjk')
      if (.not. allocated(error)) call expect_absent(error, json, 'pcm')
      if (.not. allocated(error)) call expect_absent(error, json, 'rimp2')
      if (.not. allocated(error)) call expect_absent(error, json, 'gradient')
      call json%destroy()
   end subroutine test_rhf

   subroutine test_printed_table(error)
      !! fock_first is build 1's own time and n_fock_builds the number of builds,
      !! against the iteration table the same run printed.
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp), allocatable :: printed(:)
      real(dp) :: ffirst, ftot, nfock
      integer :: n
      real(dp), parameter :: ROUNDING = 6.0e-4_dp
         !! The table prints 3 decimals, so each entry is good to 5e-4 s.
      call run_cli(error, '-b cc-pvdz --no-incremental', json, silent=.false.)
      if (allocated(error)) return
      call get_time(error, json, 'scf.fock_first', ffirst, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.fock_total', ftot, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.n_fock_builds', nfock, .true.)
      if (allocated(error)) then
         call json%destroy(); return
      end if
      call printed_fock_times(error, printed, n)
      if (.not. allocated(error)) call check(error, n == nint(nfock), "the table has "//to_string(n)// &
                                             " rows, timings.scf.n_fock_builds is "//to_string(nint(nfock)))
      if (.not. allocated(error)) call check(error, abs(printed(1) - ffirst) <= ROUNDING, &
                                             "fock_first "//to_string(ffirst)//" /= the table's build 1, "//to_string(printed(1)))
      if (.not. allocated(error)) call check(error, abs(sum(printed) - ftot) <= ROUNDING*real(n, dp), &
                                             "fock_total "//to_string(ftot)//" /= the table's sum, "//to_string(sum(printed)))
      call json%destroy()
   end subroutine test_printed_table

   subroutine test_soscf(error)
      !! Newton iterations are counted in n_iterations but are not timed Fock builds.
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: niter, nfock
      integer :: it
      logical :: found
      call run_cli(error, '-b 6-31g --soscf --soscf-start 1e-1', json)
      if (allocated(error)) return
      call get_time(error, json, 'scf.n_iterations', niter, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.n_fock_builds', nfock, .true.)
      if (.not. allocated(error)) then
         call json%get('iterations', it, found)
         call check(error, found .and. nint(niter) == it, "n_iterations /= iterations")
      end if
      if (.not. allocated(error)) call check(error, nint(nfock) < nint(niter), "the SOSCF iterations were counted as Fock builds")
      call json%destroy()
   end subroutine test_soscf

   subroutine test_total_wall_only(error)
      !! --hessian and --fukui run more SCFs after the one the run is about: only
      !! total_wall is honest there. (--opt does the same but needs a DL-FIND
      !! build, which the host CI build is not.)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: total
      call run_cli(error, '-b sto-3g --hessian', json)
      if (allocated(error)) return
      call get_time(error, json, 'total_wall', total, .true.)
      if (.not. allocated(error)) call expect_absent(error, json, 'scf')
      call json%destroy()
      if (allocated(error)) return
      call run_cli(error, '-b sto-3g --fukui hirshfeld', json)
      if (allocated(error)) return
      call get_time(error, json, 'total_wall', total, .true.)
      if (.not. allocated(error)) call expect_absent(error, json, 'scf')
      call json%destroy()
   end subroutine test_total_wall_only

   subroutine test_dft(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: grid, xc, vv10, wall
      call run_cli(error, '-b sto-3g -f wb97x-v', json)
      if (allocated(error)) return
      call get_time(error, json, 'scf.wall', wall, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'dft.grid', grid, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'dft.xc', xc, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'dft.vv10', vv10, .true.)
      if (.not. allocated(error)) call check(error, xc + vv10 <= wall + CLOCK_SLACK, "xc + vv10 exceeds the SCF wall time")
      call json%destroy()
   end subroutine test_dft

   subroutine test_dft_no_vv10(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: grid, xc
      call run_cli(error, '-b sto-3g -f pbe', json)
      if (allocated(error)) return
      call get_time(error, json, 'dft.grid', grid, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'dft.xc', xc, .true.)
      if (.not. allocated(error)) call expect_absent(error, json, 'dft.vv10')
      call json%destroy()
   end subroutine test_dft_no_vv10

   subroutine test_df(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: tot, chol, fock
      call run_cli(error, '-b 6-31g --df', json)
      if (allocated(error)) return
      call get_time(error, json, 'dfjk.total', tot, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'dfjk.step.cholesky', chol, .false.)
      if (.not. allocated(error)) call get_time(error, json, 'scf.fock_total', fock, .true.)
      if (.not. allocated(error)) call check(error, abs(tot - fock) <= 1.0e-9_dp, "dfjk.total /= scf.fock_total")
      call json%destroy()
   end subroutine test_df

   subroutine test_pcm_grad(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: setup, per, tot, g, total
      call run_cli(error, '-b 6-31g --pcm iefpcm --eps 78.39 --grad', json)
      if (allocated(error)) return
      call get_time(error, json, 'pcm.setup', setup, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'pcm.per_iteration', per, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'pcm.total', tot, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'gradient.total', g, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'total_wall', total, .true.)
      if (.not. allocated(error)) call check(error, total + CLOCK_SLACK >= g, "total_wall < gradient.total")
      if (.not. allocated(error)) call expect_absent(error, json, 'rimp2')
      call json%destroy()
   end subroutine test_pcm_grad

   subroutine test_rimp2(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: tot, chol
      call run_cli(error, '-b cc-pvdz --rimp2 --aux cc-pvdz-rifit', json)
      if (allocated(error)) return
      call get_time(error, json, 'rimp2.total', tot, .true.)
      if (.not. allocated(error)) call get_time(error, json, 'rimp2.step.cholesky', chol, .true.)
      if (.not. allocated(error)) call expect_absent(error, json, 'gradient')
      call json%destroy()
   end subroutine test_rimp2

   subroutine test_build(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      character(len=:), allocatable :: s
      logical :: found
      call run_cli(error, '-b sto-3g', json)
      if (allocated(error)) return
      call json%get('build.lmax', s, found)
      call check(error, found .and. len(s) > 0, "no build.lmax")
      if (.not. allocated(error)) then
         call json%get('build.compiler', s, found)
         call check(error, found .and. len(s) > 0, "no build.compiler")
      end if
      if (.not. allocated(error)) call check(error, json%valid_path('build.device'), "no build.device")
      call json%destroy()
   end subroutine test_build

end module check_cli_timings_suite

program check_cli_timings
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_timings_suite, only: collect_cli_timings, cli_exe, probe_json
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); probe_json = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0) then
      write (error_unit, "(a)") "check_cli_timings: usage: check_cli_timings <terco CLI> <scratch .json>"
      error stop 1
   end if
   stat = 0
   ! Serial: every case writes the same probe file.
   call run_testsuite(collect_cli_timings, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_timings
