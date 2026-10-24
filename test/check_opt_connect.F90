!
! --opt-connect (DL-FIND task 1011): after a saddle search, displace along
! the imaginary mode both ways and minimise each to a downhill minimum.
! Works with --opt-saddle-method dimer; --opt-saddle-method prfo is
! refused at the CLI (app/main.F90) because it crashes INSIDE libdlfind
! itself -- confirmed with a debug (-fcheck=all) build of libdlfind: the
! fault is a bad array copy at dlf_formstep.f90:2362 (`array(:)=
! tsmode_r(:)` in dlf_formstep_get_ra, called from dlf_task.f90's
! task=1011 case right after the P-RFO TS sub-search converges). The size
! check one line above the fault passes, so `tsmode_r` itself is a
! corrupted allocation inside libdlfind's own TS-mode bookkeeping
! (dlf_formstep_set_tsmode's mode 10/11 path, which only P-RFO/NR ever
! exercise) by the time the downhill phase asks for it -- not a size
! mismatch or missing argument this bridge supplies or could catch. The
! dimer method's own TS-mode tracking never goes through that code path
! and does not reproduce the fault, which is what this check exercises
! end to end: --opt-target saddle --opt-saddle-method dimer --opt-connect
! from test/hcn_ts_guess.xyz finds the TS and both downhill minima.
!
! Runs the CLI with WORKING_DIRECTORY set to a scratch dir (not test/
! itself): DL-FIND writes path.xyz/path_active.xyz/path_tsmode.xyz/
! path_force.xyz unconditionally into its CWD once printf>=3 (dl-find.f90
! -- confirmed by reading it; the condition is `printf`, not `printl`),
! and a `connect` run needs printf=4 for the TS/minima callbacks to fire
! at all. This check also asserts those files do NOT land in test/
! itself.
!
! Validation: the two minima are HCN and HNC. Whichever is energetically
! lower (HCN; opt_hcn.ref, -91.67520899058249 Eh) is checked against that
! EXISTING anchor directly, since a connect run from this start should
! land on the SAME HCN minimum --opt alone already validates there. The
! higher one (HNC) is checked with --hessian for ZERO imaginary modes --
! a genuine minimum, not a second saddle a bad connect run slid into.
!
! The TS sub-search's OWN convergence/gradient (trc_dlfind_bridge's own
! fix: cb_put_coords snapshots these at switch 3, before the two downhill
! minimisations run, rather than letting the LAST accepted evaluation
! overall -- the second minimisation's own -- leak into
! optimization.converged/final_gradient_max). test_ts_gradient_is_real
! below cannot easily construct a run whose TS stage hits its own
! --opt-maxiter while both minimisations still converge (this suite's own
! HCN<->HNC case converges the TS in a handful of dimer steps, long
! before either downhill run could matter), so it checks the next best
! thing instead: recomputing the gradient at the REPORTED TS geometry
! with a single --grad run and comparing maxval(abs(.)) to
! optimization.final_gradient_max. Before the fix, that quantity was the
! LAST downhill minimisation's own gmax against the TS's own geometry --
! two unrelated numbers that would not have agreed to anywhere near this
! tolerance.
!
! A test-drive suite.
!
module check_opt_connect_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   implicit none
   private
   public :: collect_opt_connect, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the xyz/JSON outputs.
   real(dp), parameter :: ANCHOR_ENERGY_TOL = 1.0e-7_dp
      !! Eh -- same reasoning as check_opt_ts.F90's own ANCHOR_FREQ_TOL
      !! header: both compare a converged terco run to a pyscf/mqc
      !! anchor well inside the optimiser's own --opt-gtol slack.
   real(dp), parameter :: GRAD_RECOMPUTE_TOL = 1.0e-7_dp
      !! Eh/Bohr -- the same geometry, same basis, same SCF converged
      !! twice (once inside --opt's last step, once by a fresh --grad
      !! run): this is round-off/SCF-convergence slack, not --opt-gtol.

   ! --- state shared between the single CLI run and the checks on it ------
   real(dp) :: e_ts = 0.0_dp, e_minus = 0.0_dp, e_plus = 0.0_dp
   real(dp) :: ts_gmax_reported = 0.0_dp
   logical :: ts_converged_reported = .false.
   logical :: min_minus_found = .false., min_plus_found = .false.
   real(dp), allocatable :: ts_xyz_angstrom(:, :)   !! (3, natm) -- row k is atom k's [x,y,z]

contains

   subroutine collect_opt_connect(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt-connect: TS and both downhill minima, dimer method, HCN<->HNC", &
                      test_connect_run), &
         new_unittest("--opt-connect: the reported TS gmax is the gradient AT the TS geometry, "// &
                      "not a downhill minimum's", test_ts_gradient_is_real) &
         ]
   end subroutine collect_opt_connect

   subroutine test_connect_run(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd
      integer :: rc, nfail, ia, k
      logical :: converged, ts_found, ts_converged_nested, found
      real(dp), allocatable :: row(:)
      type(json_file) :: json

      nfail = 0
      ! Everything the subprocess needs is copied INTO work_dir, so `cd
      ! work_dir &&` below can use plain relative paths: DL-FIND's own
      ! path*.xyz (gated on `printf`, confirmed by reading dl-find.f90)
      ! then land there too, not in test/ itself (this program's own
      ! WORKING_DIRECTORY, ctest's convention) -- see this file's own
      ! header.
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn_ts_guess.xyz '//work_dir//'/hcn_ts_guess.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)

      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn_ts_guess.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-target saddle --opt-saddle-method dimer '// &
            '--opt-connect --opt-maxiter 150 --opt-gtol 1e-5 --silent --json connect.json '// &
            '> connect.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run itself failed, exit code "//to_string(rc)); return
      end if

      call json%initialize()
      call json%load(filename=work_dir//'/connect.json')
      if (json%failed()) then
         call test_failed(error, "cannot read connect.json"); return
      end if

      ! optimization.converged/final_gradient_max are, for a connect run,
      ! the TS SEARCH's own verdict and gradient (trc_dlfind_bridge's own
      ! fix) -- never the last downhill minimisation's, which is what
      ! test_connect's own assertion on "converged" actually tests here.
      call json%get('optimization.converged', converged, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no optimization.converged"); return
      end if
      call ref_report(nfail, converged, 'check_opt_connect: the TS sub-search converged', 0.0_dp)
      ts_converged_reported = converged
      call json%get('optimization.final_gradient_max', ts_gmax_reported, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no optimization.final_gradient_max"); return
      end if

      call json%get('optimization.ts_found', ts_found, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no optimization.ts_found"); return
      end if
      call ref_report(nfail, ts_found, 'check_opt_connect: a transition state was found', 0.0_dp)
      call json%get('optimization.ts.energy', e_ts, found)
      call ref_report(nfail, found, 'check_opt_connect: optimization.ts.energy present', 0.0_dp)
      call json%get('optimization.ts.converged', ts_converged_nested, found)
      call ref_report(nfail, found, 'check_opt_connect: optimization.ts.converged present', 0.0_dp)
      if (found) then
         call ref_report(nfail, ts_converged_nested .eqv. converged, &
                         'check_opt_connect: optimization.ts.converged matches the top-level verdict '// &
                         '(both are the TS search''s own)', 0.0_dp)
      end if
      call json%get('optimization.minimum_minus_found', min_minus_found, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no optimization.minimum_minus_found"); return
      end if
      call ref_report(nfail, min_minus_found, 'check_opt_connect: the -distort downhill minimum was found', 0.0_dp)
      call json%get('optimization.minimum_minus.energy', e_minus, found)
      call ref_report(nfail, found, 'check_opt_connect: optimization.minimum_minus.energy present', 0.0_dp)
      call json%get('optimization.minimum_plus_found', min_plus_found, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no optimization.minimum_plus_found"); return
      end if
      call ref_report(nfail, min_plus_found, 'check_opt_connect: the +distort downhill minimum was found', 0.0_dp)
      call json%get('optimization.minimum_plus.energy', e_plus, found)
      call ref_report(nfail, found, 'check_opt_connect: optimization.minimum_plus.energy present', 0.0_dp)

      ! The TS's own geometry, for test_ts_gradient_is_real below --
      ! trc_dlfind_bridge overrides final_coords with the TS's for a
      ! connect run, so this is the TS, never a downhill minimum's.
      ! Counted first, THEN allocated to that size: growing
      ! ts_xyz_angstrom one column at a time by re-allocating on every
      ! row would leave the earlier columns pointing past the new
      ! (smaller, until the next row) array's own bound.
      ia = 0
      do
         call json%get('optimization.final_geometry_angstrom('//to_string(ia + 1)//')', row, found)
         if (.not. found) exit
         ia = ia + 1
      end do
      if (ia > 0) then
         allocate (ts_xyz_angstrom(3, ia))
         do k = 1, ia
            call json%get('optimization.final_geometry_angstrom('//to_string(k)//')', row, found)
            ts_xyz_angstrom(:, k) = row
         end do
      end if
      call ref_report(nfail, allocated(ts_xyz_angstrom) .and. ia == 3, &
                      'check_opt_connect: optimization.final_geometry_angstrom has 3 rows (HCN)', real(ia, dp))
      call json%destroy()

      print '(a,f20.12)', '  transition state  E = ', e_ts
      print '(a,f20.12)', '  minimum (-distort) E = ', e_minus
      print '(a,f20.12)', '  minimum (+distort) E = ', e_plus

      call ref_report(nfail, e_ts > e_minus .and. e_ts > e_plus, &
                      'check_opt_connect: the transition state sits above BOTH downhill minima', 0.0_dp)

      if (min_minus_found .and. min_plus_found) then
         if (e_minus < e_plus) then
            call check_hcn_side(e_minus, nfail)
            call check_hnc_side(cli_exe, work_dir, 'minus', nfail)
         else
            call check_hcn_side(e_plus, nfail)
            call check_hnc_side(cli_exe, work_dir, 'plus', nfail)
         end if
      end if

      ! DL-FIND's own path*.xyz/fort.32, and --opt's own
      ! <geometry>_opt.xyz/_traj.xyz/_ts.xyz/_minus.xyz/_plus.xyz
      ! (app/main.F90), must not have landed in test/ itself (this
      ! program's own WORKING_DIRECTORY, ctest's convention) -- the whole
      ! point of the `cd work_dir &&` above.
      block
         character(len=*), parameter :: stray(9) = &
            [character(len=22) :: 'path.xyz', 'path_active.xyz', 'path_tsmode.xyz', 'path_force.xyz', &
                                   'hcn_ts_guess_opt.xyz', 'hcn_ts_guess_traj.xyz', 'hcn_ts_guess_ts.xyz', &
                                   'hcn_ts_guess_minus.xyz', 'hcn_ts_guess_plus.xyz']
         integer :: k
         logical :: exists
         do k = 1, size(stray)
            inquire (file=trim(stray(k)), exist=exists)
            call ref_report(nfail, .not. exists, &
                            'check_opt_connect: '//trim(stray(k))//' was not left in test/ itself', &
                            merge(1.0_dp, 0.0_dp, exists))
         end do
      end block

      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_connect_run

   !
   ! Recomputes the gradient at the TS geometry test_connect_run just
   ! read back (a fresh single-point --grad, same basis, same molecule,
   ! atom order H/C/N as hcn_ts_guess.xyz) and compares maxval(abs(.)) to
   ! optimization.final_gradient_max from that same run. Depends on
   ! test_connect_run having already populated ts_xyz_angstrom/
   ! ts_gmax_reported -- test-drive runs a module's tests in the order
   ! collect_opt_connect lists them, not in parallel (run_testsuite's own
   ! parallel=.false. in the program below).
   !
   subroutine test_ts_gradient_is_real(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd
      integer :: rc, u, ia
      real(dp), allocatable :: row(:)
      real(dp) :: gmax_recomputed
      logical :: found
      type(json_file) :: json
      character(len=1), parameter :: ELEM(3) = ["H", "C", "N"]

      if (.not. allocated(ts_xyz_angstrom)) then
         call test_failed(error, "test_connect_run did not leave a TS geometry to recompute at"); return
      end if

      open (newunit=u, file=work_dir//'/ts_recompute.xyz', status='replace', action='write')
      write (u, '(i0)') size(ts_xyz_angstrom, 2)
      write (u, '(a)') "recomputed gradient at --opt-connect's own reported TS geometry"
      do ia = 1, size(ts_xyz_angstrom, 2)
         write (u, '(a,3(1x,f20.12))') ELEM(ia), ts_xyz_angstrom(:, ia)
      end do
      close (u)

      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' ts_recompute.xyz -b basis_sets/sto-3g.json --grad '// &
            '--silent --json grad_recompute.json > grad_recompute.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the --grad recompute CLI run itself failed, exit code "//to_string(rc)); return
      end if

      call json%initialize()
      call json%load(filename=work_dir//'/grad_recompute.json')
      if (json%failed()) then
         call test_failed(error, "cannot read grad_recompute.json"); return
      end if
      gmax_recomputed = 0.0_dp
      ia = 0
      do
         call json%get('gradient('//to_string(ia + 1)//')', row, found)
         if (.not. found) exit
         ia = ia + 1
         gmax_recomputed = max(gmax_recomputed, maxval(abs(row)))
      end do
      call json%destroy()
      if (ia == 0) then
         call test_failed(error, "no top-level 'gradient' array in grad_recompute.json"); return
      end if

      if (abs(ts_gmax_reported - gmax_recomputed) > GRAD_RECOMPUTE_TOL) then
         call test_failed(error, "optimization.final_gradient_max from the connect run ("// &
                          to_string(ts_gmax_reported)//") does not match the gradient recomputed AT that "// &
                          "same geometry ("//to_string(gmax_recomputed)//") -- it is describing a different point")
      end if
   end subroutine test_ts_gradient_is_real

   !
   ! The lower-energy minimum must be HCN: opt_hcn.ref's own anchor
   ! (test/pyscf_refs/README.md's two-tier scheme -- the SAME anchor
   ! check_opt_vs_pyscf.F90 reads for a plain --opt run).
   !
   subroutine check_hcn_side(e_hcn, nf)
      real(dp), intent(in) :: e_hcn
      integer, intent(inout) :: nf
      integer :: uu, nat_ref
      real(dp), allocatable :: z_ref(:), coords_ref(:, :)
      real(dp) :: e_ref, gmax_ref
      call ref_require('opt_hcn.ref', uu)
      call ref_read_int(uu, 'NAT', nat_ref)
      allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
      call ref_read_vector(uu, 'Z', nat_ref, z_ref)
      call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
      call ref_read_real(uu, 'E_FINAL', e_ref)
      call ref_read_real(uu, 'GMAX_FINAL', gmax_ref)
      call ref_close(uu)
      call ref_report(nf, abs(e_hcn - e_ref) < ANCHOR_ENERGY_TOL, &
                      'check_opt_connect: the lower downhill minimum matches the opt_hcn anchor, Eh', &
                      abs(e_hcn - e_ref))
   end subroutine check_hcn_side

   !
   ! The higher-energy minimum (HNC) must be a genuine minimum: --hessian
   ! on its own geometry shows zero imaginary modes.
   !
   subroutine check_hnc_side(exe, wdir, side, nf)
      character(len=*), intent(in) :: exe, wdir, side
      integer, intent(inout) :: nf
      character(len=:), allocatable :: cmd2
      type(json_file) :: hj
      integer :: rc2, n_imag
      logical :: ok
      cmd2 = 'cd '//wdir//' && '//trim(exe)//' hcn_ts_guess_'//side//'.xyz -b basis_sets/sto-3g.json '// &
             '--hessian --silent --json hnc_hess.json > hnc_hess.out 2>&1'
      call execute_command_line(cmd2, exitstat=rc2, wait=.true.)
      if (rc2 /= 0) then
         call ref_report(nf, .false., 'check_opt_connect: the HNC --hessian re-check CLI run itself failed', 0.0_dp)
         return
      end if
      call hj%initialize()
      call hj%load(filename=wdir//'/hnc_hess.json')
      if (hj%failed()) then
         call ref_report(nf, .false., 'check_opt_connect: cannot read hnc_hess.json', 0.0_dp); return
      end if
      call hj%get('frequencies.n_imaginary', n_imag, ok)
      if (.not. ok) then
         call hj%destroy()
         call ref_report(nf, .false., 'check_opt_connect: no frequencies.n_imaginary', 0.0_dp); return
      end if
      call hj%destroy()
      call ref_report(nf, n_imag == 0, &
                      'check_opt_connect: the higher downhill minimum (HNC) has zero imaginary modes', &
                      real(n_imag, dp))
   end subroutine check_hnc_side

end module check_opt_connect_suite

program check_opt_connect
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_connect_suite, only: collect_opt_connect, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_connect: usage: check_opt_connect <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_connect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_connect
