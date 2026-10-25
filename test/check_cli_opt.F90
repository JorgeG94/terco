!
! `terco water.xyz -b basis_sets/sto-3g.json --opt --json <probe>`, end to end: the
! JSON "optimization" block, and that <geometry>_opt.xyz/_traj.xyz were
! actually written with the right shape, in the style of
! check_cli_grad/check_cli_hessian.
! Numerical agreement with an external reference is check_opt_vs_pyscf.F90
! (tight) and check_opt_vs_mqc.F90 (loose sanity); this one is the CLI
! plumbing: does --opt --json write what it says it writes, in a shape a
! caller could actually parse, AND does it do so correctly on the three
! outcomes DL-FIND itself can leave --opt in (converged, stopped before
! converging, and a step that fails even after terco's own retry) --
! see this suite's own test headers below for each one's invariant.
!
! A test-drive suite.
!
module check_cli_opt_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use pyscf_ref_io, only: ref_report
   implicit none
   private
   public :: collect_cli_opt, cli_exe, probe_json, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: probe_json
      !! Where the CLI writes its --json, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the xyz outputs.
   real(dp), parameter :: TOL_SAME_POINT = 1.0e-8_dp
      !! test_truncated_keeps_accepted's own tolerance for "two CLI runs
      !! reported the SAME accepted point" -- see that test's own header.

contains

   subroutine collect_cli_opt(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt --json: shape of the optimization block and the written xyz files", &
                      test_cli_opt_json), &
         new_unittest("--opt: final_gradient_max is maxval(abs(gradient)), not a per-atom norm", &
                      test_convergence_quantity), &
         new_unittest("--opt: a run with a rejected trial reports steps < energy_evaluations, "// &
                      "and its trajectory energies strictly decrease", test_rejected_trial), &
         new_unittest("--opt: stopping right after a rejected trial reports the PRECEDING "// &
                      "accepted geometry/energy, not the rejected one", test_truncated_keeps_accepted), &
         new_unittest("--opt: a run that does not converge exits non-zero and says so", &
                      test_nonconverged_exit_code), &
         new_unittest("--opt: a step that fails even after the retry still writes _opt.xyz/_traj.xyz/json, "// &
                      "marked failed, before the process stops", test_failed_scf_outputs) &
         ]
   end subroutine collect_cli_opt

   subroutine test_cli_opt_json(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe_out
      integer :: rc, nfail, natm_json, steps, evals, u, ios
      logical :: converged, found, exists_opt, exists_traj
      type(json_file) :: json
      real(dp) :: final_energy, gmax, grms
      real(dp), allocatable :: row(:), traj_e(:)
      character(len=512) :: line
      character(len=:), allocatable :: coord_str

      nfail = 0
      call execute_command_line('mkdir -p '//work_dir, exitstat=rc, wait=.true.)
      call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
      probe_out = probe_json//'.out'

      cmd = trim(cli_exe)//' '//work_dir//'/water.xyz -b basis_sets/sto-3g.json --opt --opt-maxiter 50 '// &
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
      call json%get('optimization.converged', converged, found)
      call ref_report(nfail, found, 'check_cli_opt: no optimization.converged', 0.0_dp)
      call ref_report(nfail, converged, 'check_cli_opt: optimization.converged is true', 0.0_dp)
      call json%get('optimization.steps', steps, found)
      call ref_report(nfail, found, 'check_cli_opt: no optimization.steps', 0.0_dp)
      call ref_report(nfail, steps > 0, 'check_cli_opt: optimization.steps', real(steps, dp))
      call json%get('optimization.energy_evaluations', evals, found)
      call ref_report(nfail, found, 'check_cli_opt: no optimization.energy_evaluations', 0.0_dp)
      call ref_report(nfail, evals >= steps, 'check_cli_opt: energy_evaluations - steps', real(evals - steps, dp))
      call json%get('optimization.final_energy', final_energy, found)
      call ref_report(nfail, found, 'check_cli_opt: no optimization.final_energy', 0.0_dp)
      call ref_report(nfail, final_energy < -70.0_dp .and. final_energy > -80.0_dp, &
                      'check_cli_opt: optimization.final_energy is a plausible water RHF/STO-3G number, Eh', &
                      final_energy)
      call json%get('optimization.final_gradient_max', gmax, found)
      call ref_report(nfail, found, 'check_cli_opt: no optimization.final_gradient_max', 0.0_dp)
      call json%get('optimization.final_gradient_rms', grms, found)
      call ref_report(nfail, found, 'check_cli_opt: no optimization.final_gradient_rms', 0.0_dp)
      call ref_report(nfail, gmax >= grms, 'check_cli_opt: final_gradient_max - final_gradient_rms, Eh/Bohr', &
                      gmax - grms)
      call json%get('optimization.coordinates', coord_str, found)
      call ref_report(nfail, found, 'check_cli_opt: no optimization.coordinates', 0.0_dp)
      call ref_report(nfail, coord_str == 'cartesian', &
                      'check_cli_opt: optimization.coordinates is the default (cartesian)', 0.0_dp)
      call json%get('optimization.trajectory_energies', traj_e, found)
      call ref_report(nfail, found, 'check_cli_opt: no optimization.trajectory_energies', 0.0_dp)
      if (.not. allocated(traj_e)) allocate (traj_e(0))
      call ref_report(nfail, size(traj_e) == steps, 'check_cli_opt: size(trajectory_energies) - steps', &
                      real(size(traj_e) - steps, dp))
      call json%get('optimization.final_geometry_angstrom(1)', row, found)
      call ref_report(nfail, found .and. size(row) == 3, 'check_cli_opt: no final_geometry_angstrom(1)', 0.0_dp)
      natm_json = 0
      do
         call json%get('optimization.final_geometry_angstrom('//to_string(natm_json + 1)//')', row, found)
         if (.not. found) exit
         natm_json = natm_json + 1
      end do
      call ref_report(nfail, natm_json == 3, 'check_cli_opt: final_geometry_angstrom rows (water)', &
                      real(natm_json, dp))
      call json%destroy()

      ! --- the xyz files themselves -------------------------------------------
      inquire (file=work_dir//'/water_opt.xyz', exist=exists_opt)
      call ref_report(nfail, exists_opt, 'check_cli_opt: water_opt.xyz was written', 0.0_dp)
      inquire (file=work_dir//'/water_traj.xyz', exist=exists_traj)
      call ref_report(nfail, exists_traj, 'check_cli_opt: water_traj.xyz was written', 0.0_dp)
      if (exists_opt) then
         open (newunit=u, file=work_dir//'/water_opt.xyz', status='old', action='read', iostat=ios)
         read (u, '(a)', iostat=ios) line
         call ref_report(nfail, ios == 0 .and. trim(line) == '3', &
                         'check_cli_opt: water_opt.xyz starts with atom count 3', 0.0_dp)
         read (u, '(a)', iostat=ios) line
         call ref_report(nfail, ios == 0 .and. index(line, 'converged') == 1, &
                         'check_cli_opt: water_opt.xyz''s comment line starts with converged/NOT converged', 0.0_dp)
         close (u)
      end if
      if (exists_traj) then
         open (newunit=u, file=work_dir//'/water_traj.xyz', status='old', action='read', iostat=ios)
         block
            integer :: nlines, frames
            nlines = 0
            do
               read (u, '(a)', iostat=ios) line
               if (ios /= 0) exit
               nlines = nlines + 1
            end do
            frames = nlines/5   ! each frame is natm(1) + comment(1) + 3 atom lines
            call ref_report(nfail, frames == steps, 'check_cli_opt: traj.xyz frames - steps', &
                            real(frames - steps, dp))
         end block
         close (u)
      end if

      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_cli_opt_json

   !
   ! --opt-gtol's own help/docs/trc_dlfind_bridge's CONVERGENCE AUTHORITY
   ! comment both say final_gradient_max is the max CARTESIAN COMPONENT
   ! (matching DL-FIND's own dlf_convergence.f90 test) -- not the largest
   ! per-atom gradient NORM, a different, larger number whenever an atom's
   ! gradient has more than one non-negligible component. Checked by
   ! recomputing maxval(abs(.)) from the "gradient" array --opt forces on
   ! (app/main.F90 sets r%grad = .true. for --opt) and comparing to the
   ! JSON's own optimization.final_gradient_max, rather than trusting a
   ! second implementation of the same formula.
   !
   ! --opt-maxiter 5 on water_stretched.xyz/6-31g with a widened trust
   ! radius (--opt-max-step, same as test_rejected_trial) is deliberately
   ! NOT run to convergence: measured 2026-10-23, its step-5 gradient has
   ! an atom whose max |component| (0.1024) is LESS than another atom's
   ! own norm (0.1107) -- a converged, symmetric water geometry never
   ! shows this (every atom's largest component equals its own norm
   ! there), which is exactly why the bug this test guards survived
   ! check_opt_vs_mqc/check_opt_vs_pyscf (both compare converged results).
   !
   subroutine test_convergence_quantity(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe_out
      integer :: rc, ia
      logical :: found
      type(json_file) :: json
      real(dp) :: gmax_json, gmax_recomputed
      real(dp), allocatable :: row(:)

      probe_out = probe_json//'.out'
      cmd = trim(cli_exe)//' water_stretched.xyz -b basis_sets/6-31g.json --opt --opt-max-step 10.0 '// &
            '--opt-maxiter 5 --silent --json '//probe_json//' > '//probe_out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)

      call json%initialize()
      call json%load(filename=probe_json)
      if (json%failed()) then
         call test_failed(error, "cannot read "//probe_json); return
      end if
      call json%get('optimization.final_gradient_max', gmax_json, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no optimization.final_gradient_max"); return
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
         call test_failed(error, "no top-level 'gradient' array in the JSON (--opt should force --grad)"); return
      end if
      if (abs(gmax_json - gmax_recomputed) > 1.0e-12_dp) then
         call test_failed(error, "optimization.final_gradient_max ("//to_string(gmax_json)// &
                          ") does not match maxval(abs(gradient)) recomputed from the JSON's own gradient array ("// &
                          to_string(gmax_recomputed)//")")
      end if
   end subroutine test_convergence_quantity

   !
   ! water_stretched.xyz is far enough from equilibrium, and --opt-max-step
   ! widens the trust radius's own ceiling enough (trc_opt_options_t%
   ! max_step's own header), that DL-FIND's test_acceptance rejects at
   ! least one trial on the way down -- measured stable across several
   ! --opt-max-step values (5, 10, 20, 50 Bohr all give the same step/
   ! evaluation counts), so this is not a borderline, compiler-sensitive
   ! single-iteration fluke. steps < energy_evaluations is the one
   ! invariant a bridge that counts every evaluation as a step can never
   ! produce: steps and energy_evaluations would be equal by construction,
   ! and no comparison between them could ever catch that.
   !
   subroutine test_rejected_trial(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: rc, steps, evals, k
      real(dp), allocatable :: traj_e(:)
      character(len=:), allocatable :: probe

      probe = probe_json//'.rej'
      call run_opt_cli("water_stretched.xyz", "--opt-max-step 10.0", 60, probe, rc, steps, evals, traj_e)
      if (rc /= 0) then
         call test_failed(error, "the reference (uncapped) run did not converge, exit code "//to_string(rc)); return
      end if
      if (steps >= evals) then
         call test_failed(error, "steps ("//to_string(steps)//") >= energy_evaluations ("//to_string(evals)// &
                          ") -- no rejected trial happened, so this case no longer exercises the fix"); return
      end if
      if (.not. allocated(traj_e)) then
         call test_failed(error, "no optimization.trajectory_energies"); return
      end if
      if (size(traj_e) /= steps) then
         call test_failed(error, "trajectory_energies has "//to_string(size(traj_e))//" entries, steps is "// &
                          to_string(steps)); return
      end if
      do k = 2, size(traj_e)
         if (traj_e(k) >= traj_e(k - 1)) then
            call test_failed(error, "trajectory_energies("//to_string(k)//") = "//to_string(traj_e(k))// &
                             " did not decrease from the previous accepted step's "//to_string(traj_e(k - 1))// &
                             " -- every ACCEPTED step must be strictly lower than the one before it "// &
                             "(test_acceptance's own rule, dlf_scalestep.f90:238)")
            return
         end if
      end do
   end subroutine test_rejected_trial

   !
   ! A run stopped by --opt-maxiter immediately after a REJECTED trial
   ! must still report the geometry/energy of the step before it, never
   ! the rejected trial's own (higher) energy. Found by scanning
   ! --opt-maxiter upward one evaluation at a time (cheap: each run is a
   ! few SCF cycles on a 3-atom STO-3G/6-31G molecule) until `steps` stops
   ! growing while `energy_evaluations` keeps growing -- that evaluation is
   ! a rejected trial by test_rejected_trial's own definition -- then
   ! checking that stopping exactly there reports the SAME final_energy as
   ! stopping one evaluation earlier (the last point that WAS accepted).
   ! Before the fix (every evaluation counted as a step) this would have
   ! compared two DIFFERENT energies, the second one higher.
   !
   ! TOL_SAME_POINT is not round-off: these are two SEPARATE CLI process
   ! invocations, and on the GPU build the Fock digestion's atomic updates
   ! settle in a different order run to run, so even the identical point
   ! does not come back bit-for-bit (measured ~2e-9 Eh on an nvfortran GPU
   ! build at this size; 1e-8 held over repeated GPU runs). A leaked
   ! rejected trial's energy is a genuinely
   ! different point on the PES -- a rejection only happens because the
   ! energy moved enough to overshoot, orders of magnitude above this.
   !
   subroutine test_truncated_keeps_accepted(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: rc, steps_full, evals_full
      integer :: m, rc_m, steps_m, evals_m, rc_prev, steps_prev, evals_prev
      real(dp), allocatable :: dummy(:)
      real(dp) :: e_m, e_prev
      logical :: found_rejection
      character(len=:), allocatable :: probe_a, probe_b

      probe_a = probe_json//'.trunc_a'
      call run_opt_cli("water_stretched.xyz", "--opt-max-step 10.0", 60, probe_a, rc, steps_full, evals_full, dummy)
      if (rc /= 0) then
         call test_failed(error, "the reference (uncapped) run did not converge, exit code "//to_string(rc)); return
      end if
      if (evals_full <= 2) then
         call test_failed(error, "the reference run only took "//to_string(evals_full)// &
                          " evaluation(s), too few to scan for a rejected one"); return
      end if

      found_rejection = .false.
      probe_a = probe_json//'.trunc_prev'
      probe_b = probe_json//'.trunc_cur'
      call run_opt_cli_e("water_stretched.xyz", "--opt-max-step 10.0", 1, probe_a, rc_prev, steps_prev, evals_prev, &
                         e_prev)
      do m = 2, evals_full - 1
         call run_opt_cli_e("water_stretched.xyz", "--opt-max-step 10.0", m, probe_b, rc_m, steps_m, evals_m, e_m)
         if (evals_m == m .and. steps_m == steps_prev) then
            ! Evaluation m added no new accepted step: it was a rejected
            ! trial. The run truncated right after it must report the
            ! SAME final_energy as the run truncated right before it.
            found_rejection = .true.
            if (abs(e_m - e_prev) > TOL_SAME_POINT) then
               call test_failed(error, "--opt-maxiter "//to_string(m)//" (stopped right after a rejected "// &
                                "trial) reports final_energy "//to_string(e_m)//", but --opt-maxiter "// &
                                to_string(m - 1)//" (the preceding accepted step) reports "//to_string(e_prev)// &
                                " -- a rejected trial's energy leaked into the reported result")
            end if
            exit
         end if
         rc_prev = rc_m; steps_prev = steps_m; evals_prev = evals_m; e_prev = e_m
      end do
      if (.not. found_rejection) then
         call test_failed(error, "scanned --opt-maxiter 2.."//to_string(evals_full - 1)// &
                          " and found no rejected trial -- test_rejected_trial already found steps < "// &
                          "evaluations for this same run, so this scan should have found where")
      end if
   end subroutine test_truncated_keeps_accepted

   !
   ! --opt-maxiter 1 on a geometry far from equilibrium cannot possibly
   ! reach --opt-gtol's default in one accepted step -- a CLI run that
   ! does not converge must not exit 0 (the SCF's own non-convergence
   ! already works this way; this is the same rule for --opt).
   !
   subroutine test_nonconverged_exit_code(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe_out
      integer :: rc
      logical :: converged, found
      type(json_file) :: json

      probe_out = probe_json//'.short.out'
      cmd = trim(cli_exe)//' water_stretched.xyz -b basis_sets/sto-3g.json --opt --opt-maxiter 1 --silent --json '// &
            probe_json//'.short > '//probe_out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc == 0) then
         call test_failed(error, "--opt-maxiter 1 converged by exit code 0, which makes this case useless; "// &
                          "pick a geometry/tolerance that genuinely needs more than one step"); return
      end if

      call json%initialize()
      call json%load(filename=probe_json//'.short')
      if (json%failed()) then
         call test_failed(error, "cannot read "//probe_json//'.short'); return
      end if
      call json%get('optimization.converged', converged, found)
      call json%destroy()
      if (.not. found) then
         call test_failed(error, "no optimization.converged in the JSON"); return
      end if
      if (converged) then
         call test_failed(error, "the JSON says optimization.converged is true, but the process exited non-zero")
      end if
   end subroutine test_nonconverged_exit_code

   !
   ! TRC_OPT_TEST_FAIL_AFTER_EVALS (trc_dlfind_bridge%cb_get_gradient, test-
   ! only, read only there) forces the SECOND evaluation to fail -- cheap
   ! and deterministic, unlike trying to provoke a real SCF non-convergence.
   ! DL-FIND has no recovery from a failed evaluation (api.f90's `dlf_fail`
   ! ends in `error stop`), so api_dl_find never returns on this path: the
   ! only way _opt.xyz/_traj.xyz/the JSON block get written at all is
   ! fail_eval writing them itself, before handing status=1 back. The
   ! first evaluation (water's own
   ! equilibrium-ish geometry) is left to succeed and be accepted, so the
   ! files below carry a REAL accepted geometry, not the empty starting
   ! point a failure on evaluation 1 would leave behind.
   !
   subroutine test_failed_scf_outputs(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe_out, fail_work
      integer :: rc, u, ios, nfail
      logical :: found, exists_opt, exists_traj, failed_json
      type(json_file) :: json
      character(len=512) :: line
      character(len=:), allocatable :: reason_str

      nfail = 0
      fail_work = work_dir//'/fail'
      call execute_command_line('mkdir -p '//fail_work, exitstat=rc, wait=.true.)
      call execute_command_line('cp -f water.xyz '//fail_work//'/water.xyz', exitstat=rc, wait=.true.)
      probe_out = probe_json//'.fail.out'

      cmd = 'TRC_OPT_TEST_FAIL_AFTER_EVALS=2 '//trim(cli_exe)//' '//fail_work//'/water.xyz -b basis_sets/sto-3g.json '// &
            '--opt --opt-maxiter 50 --silent --json '//probe_json//'.fail > '//probe_out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call ref_report(nfail, rc /= 0, 'check_cli_opt: a forced failed evaluation exits non-zero', real(rc, dp))

      call json%initialize()
      call json%load(filename=probe_json//'.fail')
      call ref_report(nfail, .not. json%failed(), 'check_cli_opt: the JSON was still written on a forced failure', &
                      0.0_dp)
      if (.not. json%failed()) then
         call json%get('optimization.failed', failed_json, found)
         call ref_report(nfail, found .and. failed_json, 'check_cli_opt: optimization.failed is true', 0.0_dp)
         call json%get('optimization.failure_reason', reason_str, found)
         call ref_report(nfail, found .and. index(reason_str, 'TRC_OPT_TEST_FAIL_AFTER_EVALS') > 0, &
                         'check_cli_opt: optimization.failure_reason names the forced failure', 0.0_dp)
      end if
      call json%destroy()

      inquire (file=fail_work//'/water_opt.xyz', exist=exists_opt)
      call ref_report(nfail, exists_opt, 'check_cli_opt: water_opt.xyz was written on a forced failure', 0.0_dp)
      if (exists_opt) then
         open (newunit=u, file=fail_work//'/water_opt.xyz', status='old', action='read', iostat=ios)
         read (u, '(a)', iostat=ios)
         read (u, '(a)', iostat=ios) line
         call ref_report(nfail, ios == 0 .and. index(line, 'FAILED') > 0, &
                         'check_cli_opt: water_opt.xyz''s comment line says FAILED', 0.0_dp)
         close (u)
      end if
      inquire (file=fail_work//'/water_traj.xyz', exist=exists_traj)
      call ref_report(nfail, exists_traj, 'check_cli_opt: water_traj.xyz was written on a forced failure', 0.0_dp)

      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_failed_scf_outputs

   !
   ! Runs the CLI's --opt on `xyz` with `extra_flags` and --opt-maxiter
   ! `maxiter`, and reads back steps/energy_evaluations/trajectory_energies.
   ! `rc` is the process exit code -- 0 only if --opt converged (a
   ! deliberately small maxiter is expected to exit non-zero; callers that
   ! need convergence check `rc` themselves).
   !
   subroutine run_opt_cli(xyz, extra_flags, maxiter, probe, rc, steps, evals, traj_e)
      character(len=*), intent(in) :: xyz, extra_flags, probe
      integer, intent(in) :: maxiter
      integer, intent(out) :: rc, steps, evals
      real(dp), allocatable, intent(out) :: traj_e(:)
      character(len=:), allocatable :: cmd
      logical :: found
      type(json_file) :: json

      cmd = trim(cli_exe)//' '//xyz//' -b basis_sets/6-31g.json --opt --opt-maxiter '//to_string(maxiter)//' '// &
            trim(extra_flags)//' --silent --json '//probe//' > '//probe//'.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      steps = -1; evals = -1
      call json%initialize()
      call json%load(filename=probe)
      if (json%failed()) then
         call json%destroy(); return
      end if
      call json%get('optimization.steps', steps, found)
      if (.not. found) steps = -1
      call json%get('optimization.energy_evaluations', evals, found)
      if (.not. found) evals = -1
      call json%get('optimization.trajectory_energies', traj_e, found)
      call json%destroy()
   end subroutine run_opt_cli

   !
   ! Same as run_opt_cli, but also returns final_energy -- the one extra
   ! field test_truncated_keeps_accepted compares between two consecutive
   ! --opt-maxiter values.
   !
   subroutine run_opt_cli_e(xyz, extra_flags, maxiter, probe, rc, steps, evals, final_energy)
      character(len=*), intent(in) :: xyz, extra_flags, probe
      integer, intent(in) :: maxiter
      integer, intent(out) :: rc, steps, evals
      real(dp), intent(out) :: final_energy
      character(len=:), allocatable :: cmd
      logical :: found
      type(json_file) :: json

      cmd = trim(cli_exe)//' '//xyz//' -b basis_sets/6-31g.json --opt --opt-maxiter '//to_string(maxiter)//' '// &
            trim(extra_flags)//' --silent --json '//probe//' > '//probe//'.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      steps = -1; evals = -1; final_energy = 0.0_dp
      call json%initialize()
      call json%load(filename=probe)
      if (json%failed()) then
         call json%destroy(); return
      end if
      call json%get('optimization.steps', steps, found)
      if (.not. found) steps = -1
      call json%get('optimization.energy_evaluations', evals, found)
      if (.not. found) evals = -1
      call json%get('optimization.final_energy', final_energy, found)
      call json%destroy()
   end subroutine run_opt_cli_e

end module check_cli_opt_suite

program check_cli_opt
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_opt_suite, only: collect_cli_opt, cli_exe, probe_json, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); probe_json = trim(arg)
   call get_command_argument(3, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_cli_opt: usage: check_cli_opt <terco CLI> <scratch .json> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_cli_opt, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_opt
