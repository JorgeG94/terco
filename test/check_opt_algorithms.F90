!
! Every --opt-algorithm terco exposes (sd, cg, cg-auto, lbfgs, damped,
! nr, prfo-as-a-minimiser) converges water/STO-3G to the existing
! opt_water anchor energy (test/pyscf_refs/opt_water.ref) within a
! tolerance scaled to the --opt-gtol used -- the harmonic estimate
! check_opt_vs_pyscf.F90 already uses (dE ~ g^2/(2k)), generously
! loosened here since this check's job is "does the algorithm converge
! at all", not a tight anchor comparison (that is check_opt_vs_pyscf's).
! sd gets a generous maxiter (it is the slowest of the seven by a wide
! margin); every other algorithm uses the CLI's own default (100).
!
! Every --opt-coords terco exposes (cartesian, hdlc, hdlc-tc, dlc,
! dlc-tc) converges the SAME water minimisation with --opt-algorithm
! lbfgs, to the same anchor.
!
! A test-drive suite.
!
module check_opt_algorithms_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_report
   implicit none
   private
   public :: collect_opt_algorithms, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the xyz/JSON outputs --
      !! shared by every combination below: they are sequential
      !! (test-drive's own parallel=.false.), never concurrent.
   real(dp), parameter :: E_TOL = 2.0e-6_dp
      !! loose: "does it converge", not a tight anchor match.

contains

   subroutine collect_opt_algorithms(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt-algorithm sd converges water to the opt_water anchor", test_sd), &
         new_unittest("--opt-algorithm cg converges water to the opt_water anchor", test_cg), &
         new_unittest("--opt-algorithm cg-auto converges water to the opt_water anchor", test_cg_auto), &
         new_unittest("--opt-algorithm lbfgs converges water to the opt_water anchor", test_lbfgs), &
         new_unittest("--opt-algorithm damped converges water to the opt_water anchor", test_damped), &
         new_unittest("--opt-algorithm nr converges water to the opt_water anchor", test_nr), &
         new_unittest("--opt-algorithm prfo (as a plain minimiser) converges water to the "// &
                      "opt_water anchor", test_prfo), &
         new_unittest("--opt-coords hdlc (lbfgs) converges water to the opt_water anchor", test_hdlc), &
         new_unittest("--opt-coords hdlc-tc (lbfgs) converges water to the opt_water anchor", test_hdlc_tc), &
         new_unittest("--opt-coords dlc (lbfgs) converges water to the opt_water anchor", test_dlc), &
         new_unittest("--opt-coords dlc-tc (lbfgs) converges water to the opt_water anchor", test_dlc_tc), &
         new_unittest("--opt-algorithm: no probe JSON or xyz output was left in test/ itself", test_no_stray) &
         ]
   end subroutine collect_opt_algorithms

   subroutine test_sd(error); type(error_type), allocatable, intent(out) :: error
      call run_one('sd', 'cartesian', 600, error)
   end subroutine test_sd
   subroutine test_cg(error); type(error_type), allocatable, intent(out) :: error
      call run_one('cg', 'cartesian', 100, error)
   end subroutine test_cg
   subroutine test_cg_auto(error); type(error_type), allocatable, intent(out) :: error
      call run_one('cg-auto', 'cartesian', 100, error)
   end subroutine test_cg_auto
   subroutine test_lbfgs(error); type(error_type), allocatable, intent(out) :: error
      call run_one('lbfgs', 'cartesian', 100, error)
   end subroutine test_lbfgs
   subroutine test_damped(error); type(error_type), allocatable, intent(out) :: error
      call run_one('damped', 'cartesian', 150, error)
   end subroutine test_damped
   subroutine test_nr(error); type(error_type), allocatable, intent(out) :: error
      call run_one('nr', 'cartesian', 100, error)
   end subroutine test_nr
   subroutine test_prfo(error); type(error_type), allocatable, intent(out) :: error
      ! P-RFO is valid as a plain minimiser too, just unusual.
      call run_one('prfo', 'cartesian', 400, error)
   end subroutine test_prfo
   subroutine test_hdlc(error); type(error_type), allocatable, intent(out) :: error
      call run_one('lbfgs', 'hdlc', 100, error)
   end subroutine test_hdlc
   subroutine test_hdlc_tc(error); type(error_type), allocatable, intent(out) :: error
      call run_one('lbfgs', 'hdlc-tc', 100, error)
   end subroutine test_hdlc_tc
   subroutine test_dlc(error); type(error_type), allocatable, intent(out) :: error
      call run_one('lbfgs', 'dlc', 100, error)
   end subroutine test_dlc
   subroutine test_dlc_tc(error); type(error_type), allocatable, intent(out) :: error
      call run_one('lbfgs', 'dlc-tc', 100, error)
   end subroutine test_dlc_tc

   !
   ! water_opt.xyz/water_traj.xyz and every alg_*.json probe must not have
   ! landed in test/ itself -- the whole point of each run_one's own
   ! `cd work_dir &&`. Listed last (collect_opt_algorithms' own order,
   ! test-drive's sequential, non-parallel run) so every other case has
   ! already run and could have leaked a file by the time this checks.
   !
   subroutine test_no_stray(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: algs(7) = &
         [character(len=8) :: 'sd', 'cg', 'cg-auto', 'lbfgs', 'damped', 'nr', 'prfo']
      character(len=*), parameter :: coordsets(4) = &
         [character(len=8) :: 'hdlc', 'hdlc-tc', 'dlc', 'dlc-tc']
      integer :: k, nfail
      logical :: exists

      nfail = 0
      inquire (file='water_opt.xyz', exist=exists)
      call ref_report(nfail, .not. exists, 'check_opt_algorithms: water_opt.xyz was not left in test/ itself', &
                      merge(1.0_dp, 0.0_dp, exists))
      inquire (file='water_traj.xyz', exist=exists)
      call ref_report(nfail, .not. exists, 'check_opt_algorithms: water_traj.xyz was not left in test/ itself', &
                      merge(1.0_dp, 0.0_dp, exists))
      do k = 1, size(algs)
         inquire (file='alg_'//trim(algs(k))//'_cartesian.json', exist=exists)
         call ref_report(nfail, .not. exists, 'check_opt_algorithms: alg_'//trim(algs(k))// &
                         '_cartesian.json was not left in test/ itself', merge(1.0_dp, 0.0_dp, exists))
      end do
      do k = 1, size(coordsets)
         inquire (file='alg_lbfgs_'//trim(coordsets(k))//'.json', exist=exists)
         call ref_report(nfail, .not. exists, 'check_opt_algorithms: alg_lbfgs_'//trim(coordsets(k))// &
                         '.json was not left in test/ itself', merge(1.0_dp, 0.0_dp, exists))
      end do
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_no_stray

   subroutine run_one(alg, coords, maxiter, error)
      character(len=*), intent(in) :: alg, coords
      integer, intent(in) :: maxiter
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, tag, probe
      type(json_file) :: json
      integer :: rc, steps, nfail
      real(dp) :: e
      logical :: ok

      nfail = 0
      ! Run with WORKING_DIRECTORY set to work_dir (not test/ itself, this
      ! program's own WORKING_DIRECTORY, ctest's convention) -- --opt
      ! writes <geometry>_opt.xyz/_traj.xyz AND the --json probe into the
      ! CLI's own CWD; check_opt_connect.F90's own `cd work_dir &&`
      ! pattern. mkdir/cp are idempotent, so every case redoes them rather
      ! than depend on run order.
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)

      tag = alg//'_'//coords
      probe = 'alg_'//tag//'.json'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt --opt-algorithm '// &
            alg//' --opt-coords '//coords//' --opt-maxiter '//to_string(maxiter)//' --opt-gtol 1e-5 --silent '// &
            '--json '//probe//' > '//tag//'.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run itself failed, exit code "//to_string(rc)); return
      end if

      call json%initialize()
      call json%load(filename=work_dir//'/'//probe)
      if (json%failed()) then
         call test_failed(error, "cannot read "//probe); return
      end if
      call json%get('optimization.final_energy', e, ok)
      call ref_report(nfail, ok, 'check_opt_algorithms ('//tag//'): optimization.final_energy present', 0.0_dp)
      call json%get('optimization.steps', steps, ok)
      call json%destroy()
      if (.not. ok) then
         call test_failed(error, "no optimization.steps"); return
      end if

      print '(a,a10,a,a9,a,i5,a,f20.12)', '  ', alg, '  ', coords, '  steps=', steps, '  E=', e
      call ref_report(nfail, abs(e - water_anchor_energy()) < E_TOL, &
                      'check_opt_algorithms ('//tag//'): final energy matches the opt_water anchor', &
                      abs(e - water_anchor_energy()))
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine run_one

   !
   ! opt_water.ref's own converged energy (test/pyscf_refs/README.md's
   ! two-tier scheme) -- read from the committed reference rather than
   ! restated as a literal, so this can never drift from it. Cached after
   ! the first read: every run_one call needs it, and the file never
   ! changes mid-process.
   !
   function water_anchor_energy() result(e)
      use pyscf_ref_io, only: ref_read_int, ref_read_vector, ref_read_matrix
      real(dp) :: e
      real(dp), save :: cached = 0.0_dp
      logical, save :: have_cached = .false.
      integer :: uu, nat_ref
      real(dp), allocatable :: z_ref(:), coords_ref(:, :)
      if (.not. have_cached) then
         ! opt_water.ref's tags are read strictly in the order
         ! generate.py wrote them (NAT, Z, COORDS, E_FINAL, GMAX_FINAL):
         ! ref_read_* has no random access, so NAT/Z/COORDS have to be
         ! read (and discarded) before E_FINAL, the same sequence
         ! check_opt_vs_pyscf.F90/check_opt_connect.F90's check_hcn_side
         ! use.
         call ref_require('opt_water.ref', uu)
         call ref_read_int(uu, 'NAT', nat_ref)
         allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
         call ref_read_vector(uu, 'Z', nat_ref, z_ref)
         call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
         call ref_read_real(uu, 'E_FINAL', cached)
         call ref_close(uu)
         have_cached = .true.
      end if
      e = cached
   end function water_anchor_energy

end module check_opt_algorithms_suite

program check_opt_algorithms
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_algorithms_suite, only: collect_opt_algorithms, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_algorithms: usage: check_opt_algorithms <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_algorithms, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_algorithms
