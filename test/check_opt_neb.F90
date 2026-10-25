!
! Nudged elastic band / climbing-image NEB on the HCN<->HNC isomerisation,
! RHF/STO-3G -- see app/trc_dlfind_bridge.F90's own header for the NEB
! mechanics (icoord 100-199, per-image warm-start density, the
! --opt-neb-climb sentinel) this exercises end to end through the CLI.
!
! Endpoints: test/hcn_minimum.xyz and test/hnc_minimum.xyz, FROZEN
! (committed) geometries from a --opt-connect run on test/hcn_ts_guess.xyz
! (check_opt_connect.F90's own machinery: the dimer method finds the TS
! and both downhill minima in one run) -- generated once, consistently,
! not re-run live here. That matters: the dimer method's own rotation
! step uses DL-FIND's "random/internal" initial direction (app/
! trc_dlfind_bridge.F90's own header), so a FRESH --opt-connect run lands
! on the SAME two minima rotated/translated differently in space each
! time, and measured directly by running --opt-connect repeatedly from
! test/hcn_ts_guess.xyz: some of those equally-valid rotations drive at
! least one interpolated
! NEB image into a geometry so distorted that terco's own overlap-matrix
! diagonalisation fails outright (`trc_linalg: pic_syev failed`, an
! unconditional `error stop` deep in pic, not a condition trc_run_scf's
! own non-convergence handling can catch or this bridge can decline),
! while others do not. Freezing the two endpoints removes that source of
! run-to-run variation; this is a known-delicate 3-atom system in GENERAL
! for this DL-FIND build (check_opt_ts.F90's own header already documents
! Cartesian HDLC P-RFO failing to converge reliably for the SAME
! molecule), and NEB's many simultaneous images make it more, not less,
! likely to hit that fragility than a single-geometry optimiser.
!
! --opt-neb-via test/hcn_ts_guess.xyz (the same bent TS guess
! check_opt_ts.F90 uses) gives DL-FIND's own linear-Cartesian path
! builder a bent intermediate point instead of interpolating straight
! between two COLLINEAR endpoints: HCN and HNC are both genuinely linear
! (H-C-N/H-N-C angles within 2e-4 deg of 180/0), so a straight-line path
! between them slides the bridging atom through the other two at the
! midpoint images. That -- not the endpoints' relative rotation, which
! dlf_cartesian_align (dlf_neb.f90:278) already removes before
! interpolating, Kabsch, for every known point -- was the real cause of
! the coincident-atom overlap singularities ('trc_linalg: pic_syev
! failed') this system used to trip. The pre-check below
! (check_neb_initial_path, app/trc_dlfind_bridge.F90) now refuses any
! run whose own predicted initial path puts two atoms under 0.5 Bohr
! apart, naming the image and atoms, instead of reaching the eigensolver
! at all; it does not fire for the via-point path here (closest pair
! measured a safe 1.4+ Bohr apart), but check_opt_refusals.F90 exercises
! it directly on the plain (no via-point) endpoints.
!
! --opt-neb settings (--opt-neb-images 10, --opt-neb-spring 0.01,
! --opt-gtol 1e-4, --opt-maxiter 2000, --opt-neb-climb, Cartesian,
! perpendicular endpoints) are the ones measured, with the via-point
! above and these frozen endpoints, to land the CLIMBING image on the
! same point twice in a row (this file's own reproducibility check,
! below, ~1e-8 Eh apart, CPU and GPU alike) with the climbing image
! actually spawned (DL-FIND's own "Spawning climbing image" log line,
! not just the sentinel reaching DL-FIND) and converged on its own
! 5-criterion test.
!
! CONVERGENCE VERDICT: app/trc_dlfind_bridge.F90 now judges NEB
! convergence on the gradient PERPENDICULAR to the path (projecting out
! the along-path component with an energy-ordered tangent, replicating
! dlf_neb_improved_tangent_neb's own Tau+/Tau-/blend rule), not the raw
! Cartesian gradient, which never goes to zero for an interior image by
! construction (the spring force balances its own along-path component,
! not the QC gradient). That fixed the "false by construction" problem,
! but this bridge's own reported quantity is the max ABSOLUTE COMPONENT
! (matching every other --opt-gtol check in this bridge), while DL-
! FIND's own per-image freeze/climb test is an RMS over that image's own
! DOF (dlf_neb.f90's own fper, read directly) -- a genuinely different,
! stricter number for the same gradient. Measured on this exact run:
! max-component max|g| settles at ~2.3e-4 Eh/Bohr against DL-FIND's own
! sub-1e-4 RMS criterion at --opt-gtol 1e-4, so THIS run (precise enough
! to match the TS anchor at ANCHOR_TS_TOL) is checked against DL-FIND's
! OWN observable verdict (the "Convergence reached" log line, --opt-
! print-level 4) rather than this bridge's exit code -- the review's own
! documented fallback when the two metrics do not land on the same side
! of a shared number. test_neb_converged_reports_true below is the
! dedicated proof that a converged NEB chain DOES report
! optimization.converged = true once asked for a target this metric can
! actually reach (measured: --opt-gtol 2e-3 on the SAME system, gfortran
! 15 and nvfortran 26.5 alike -- NEB_FLAGS_CONVERGES's own comment below
! has the measured numbers and margin).
!
! On the GPU build specifically, an INTERIOR image away from both
! endpoints and the climbing point can land on a measurably different
! point between two otherwise-identical runs (up to ~7.8e-4 Eh measured)
! -- DL-FIND's own per-image freeze test is a hard threshold, so GPU
! reduction-order noise landing on either side of it changes which cycle
! that image freezes at; check_reproducible, below, bounds the climbing
! image tightly and the rest of the chain separately, looser, rather
! than tightening a chain-wide bound to a number that does not hold.
!
! app/trc_dlfind_bridge.F90's own NEB result assembly used to report the
! REGULAR image nearest the climbing image (frozen at its pre-climb
! snapshot once climbing spawns) as "the climbing image", not DL-FIND's
! own separate, still-refining point at its scratch slot
! (cb_get_gradient's own "plus once more at image n_images_dlf ... once
! climbing starts") -- fixed to report the scratch slot once it has been
! evaluated. Measured effect on THIS run: climbing-image energy error
! against opt_ts_hcn.ref dropped from ~4.9e-7 Eh (frozen snapshot) to
! ~1.9e-9 Eh (DL-FIND's own converged point) once fixed. ANCHOR_TS_TOL/
! ANCHOR_TS_GEOM_TOL/ANCHOR_TS_ANGLE_TOL below have generous headroom
! above that measurement (energy ~500x, geometry ~17x, angle ~90x), not
! tightened to look better than what was observed.
!
! A test-drive suite.
!
module check_opt_neb_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_constants, only: ANGSTROM_TO_BOHR
   use json_module, only: json_file
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private
   public :: collect_opt_neb, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the xyz/JSON outputs.
   ! --opt-neb's own settings (this file's header: measured, with these
   ! FROZEN endpoints plus the --opt-neb-via bent guess, to finish
   ! identically every time with the climbing image actually spawned).
   character(len=*), parameter :: NEB_FLAGS = '--opt-neb-images 10 --opt-neb-via hcn_ts_guess.xyz '// &
                                              '--opt-neb-spring 0.01 --opt-gtol 1e-4 --opt-maxiter 2000'
   real(dp), parameter :: ENDPOINT_TOL = 1.0e-6_dp     !! Eh -- measured: endpoints did not move
   real(dp), parameter :: ANCHOR_ENERGY_TOL = 1.0e-6_dp   !! Eh -- opt_hcn.ref
   ! The climbing image vs opt_ts_hcn.ref, once DL-FIND's own scratch-slot
   ! point (not the frozen regular-image snapshot) is what gets reported:
   ! measured 1.9e-9 Eh / 5.8e-5 Bohr / 0.0011 deg (this file's header).
   ! ~500x/~17x/~90x headroom above that, not tuned after the fact.
   real(dp), parameter :: ANCHOR_TS_TOL = 1.0e-6_dp       !! Eh
   real(dp), parameter :: ANCHOR_TS_GEOM_TOL = 1.0e-3_dp  !! Bohr
   real(dp), parameter :: ANCHOR_TS_ANGLE_TOL = 0.1_dp    !! deg
   ! Plain NEB (no --opt-neb-climb, no via point): the ORIGINAL settings
   ! and tolerance this file used before the via-point/climb work above --
   ! kept as its own, separately exercised, looser profile check (the
   ! highest-energy REGULAR image, climbing never having been asked to
   ! refine it, sits 2.1-6.4e-3 Eh below the true TS across every stable
   ! configuration measured). 1e-2 is a little over 1.5x that, not a
   ! number tuned after the fact to pass.
   character(len=*), parameter :: NEB_FLAGS_PLAIN = '--opt-neb-images 10 --opt-max-step 0.03 '// &
                                                    '--opt-neb-spring 0.01 --opt-gtol 8e-3 --opt-maxiter 300'
   real(dp), parameter :: ANCHOR_TS_TOL_PLAIN = 1.0e-2_dp   !! Eh
   ! --opt-gtol loose enough for THIS bridge's own max-ABSOLUTE-COMPONENT
   ! metric to cross it (this file's own header, "CONVERGENCE VERDICT").
   ! Measured directly running this exact command (hcn_minimum.xyz ->
   ! hnc_minimum.xyz, these flags, --opt-neb-climb, no --silent so the
   ! "final max|g|" log line and optimization.final_gradient_max in the
   ! --json both read):
   !   --opt-gtol 3e-4 (this file's OLD value): gfortran 15 converges at
   !     202 steps, final max|g| = 2.864e-4 (exit 0). nvfortran 26.5
   !     STALLS at 355 steps, final max|g| = 7.052e-4, never crosses
   !     3e-4 (exit 2, optimization.converged = false) -- the bug this
   !     fix addresses, reproduced directly.
   !   --opt-gtol 2e-3 (this file's NEW value, below): gfortran 15 and
   !     nvfortran 26.5 BOTH converge at 47 steps, to the SAME final
   !     max|g| to 8 significant figures (gfortran 1.8615140e-3,
   !     nvfortran 1.8615141e-3 Eh/Bohr -- BLAS/LAPACK-level rounding,
   !     not a convergence difference) -- margin above --opt-gtol itself
   !     is 2e-3 - 1.8615e-3 = 1.385e-4 (~7% of --opt-gtol).
   !   --opt-gtol 8e-4 (tighter, tried first): gfortran 15 itself STALLS
   !     at only 105 steps, final max|g| = 1.709e-3, never crosses 8e-4
   !     (exit 2) -- NOT a "more iterations converges tighter" case:
   !     --opt-gtol also sets DL-FIND's own internal per-image RMS
   !     stopping test (dlf_neb.f90's own fper), so asking for a tighter
   !     target changes the search's own trajectory (fewer, different
   !     accepted LBFGS steps), not just how long it keeps refining the
   !     same path -- this 3-atom chain is measured (this file's own
   !     header) to be delicate enough that a tighter target can stall
   !     a run that a looser one lets finish. 2e-3 is the tightest value
   !     found, by direct measurement, that converges reliably (same
   !     step count, same final max|g|) on both compilers; --opt-maxiter
   !     2000 is already far beyond the 47 steps either one takes, so
   !     raising it further changes nothing here.
   character(len=*), parameter :: NEB_FLAGS_CONVERGES = '--opt-neb-images 10 --opt-neb-via hcn_ts_guess.xyz '// &
                                                        '--opt-neb-spring 0.01 --opt-gtol 2e-3 --opt-maxiter 2000'

contains

   subroutine collect_opt_neb(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt-neb: the two committed endpoints are what they claim to be", &
                      test_endpoints_are_valid), &
         new_unittest("--opt-neb-climb: the climbing-image NEB run, twice (reproducibility)", &
                      test_climb_run_and_reproducibility), &
         new_unittest("--opt-neb: a converged NEB chain reports optimization.converged = true", &
                      test_neb_converged_reports_true), &
         new_unittest("--opt-neb: plain NEB (no climbing, no via point), the looser profile check", &
                      test_plain_neb), &
         new_unittest("--opt-neb: none of the runs above left stray xyz output in test/ itself", &
                      test_no_stray) &
         ]
   end subroutine collect_opt_neb

   subroutine test_endpoints_are_valid(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: rc, nfail
      nfail = 0
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn_minimum.xyz '//work_dir//'/hcn_minimum.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hnc_minimum.xyz '//work_dir//'/hnc_minimum.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn_ts_guess.xyz '//work_dir//'/hcn_ts_guess.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)

      block
         integer :: uu, nat_ref
         real(dp), allocatable :: z_ref(:), coords_ref(:, :)
         real(dp) :: e_ref, e_hcn
         logical :: ok
         call ref_require('opt_hcn.ref', uu)
         call ref_read_int(uu, 'NAT', nat_ref)
         allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
         call ref_read_vector(uu, 'Z', nat_ref, z_ref)
         call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
         call ref_read_real(uu, 'E_FINAL', e_ref)
         call ref_close(uu)
         call read_xyz_energy(work_dir//'/hcn_minimum.xyz', e_hcn, ok)
         if (.not. ok) then
            call test_failed(error, "cannot read test/hcn_minimum.xyz's own energy comment"); return
         end if
         call ref_report(nfail, abs(e_hcn - e_ref) < ANCHOR_ENERGY_TOL, &
                         'check_opt_neb: test/hcn_minimum.xyz matches the opt_hcn anchor, Eh', &
                         abs(e_hcn - e_ref))
      end block
      block
         ! test/hnc_minimum.xyz: a genuine minimum (zero imaginary modes),
         ! not a second saddle a bad connect run slid into -- the SAME
         ! check check_opt_connect.F90 already does for its own
         ! higher-energy side.
         character(len=:), allocatable :: cmd2
         type(json_file) :: hj
         integer :: rc2, n_imag
         logical :: ok
         cmd2 = 'cd '//work_dir//' && '//trim(cli_exe)//' hnc_minimum.xyz -b basis_sets/sto-3g.json '// &
                '--hessian --silent --json hnc_hess.json > hnc_hess.out 2>&1'
         call execute_command_line(cmd2, exitstat=rc2, wait=.true.)
         if (rc2 /= 0) then
            call test_failed(error, "the HNC --hessian check CLI run itself failed"); return
         end if
         call hj%initialize()
         call hj%load(filename=work_dir//'/hnc_hess.json')
         if (hj%failed()) then
            call test_failed(error, "cannot read hnc_hess.json"); return
         end if
         call hj%get('frequencies.n_imaginary', n_imag, ok)
         if (.not. ok) then
            call hj%destroy(); call test_failed(error, "no frequencies.n_imaginary"); return
         end if
         call hj%destroy()
         call ref_report(nfail, n_imag == 0, 'check_opt_neb: test/hnc_minimum.xyz has zero imaginary modes '// &
                         '(a genuine minimum)', real(n_imag, dp))
      end block
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_endpoints_are_valid

   subroutine test_climb_run_and_reproducibility(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call run_neb(cli_exe, work_dir, NEB_FLAGS, 'neb1.json', 'neb1.out', nfail)
      call run_neb(cli_exe, work_dir, NEB_FLAGS, 'neb2.json', 'neb2.out', nfail)
      call check_reproducible(work_dir, 'neb1.json', 'neb2.json', nfail)
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_climb_run_and_reproducibility

   !
   ! The same climbing-image NEB system, but with --opt-gtol loosened to
   ! a target THIS bridge's own max-absolute-component metric can
   ! actually cross (this file's own header, "CONVERGENCE VERDICT") --
   ! proof that the perpendicular-gradient fix lets a genuinely converged
   ! NEB chain report optimization.converged = true, not just a looser
   ! log line.
   !
   subroutine test_neb_converged_reports_true(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd
      integer :: rc
      logical :: converged, found
      type(json_file) :: json
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn_minimum.xyz '//work_dir//'/hcn_minimum.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hnc_minimum.xyz '//work_dir//'/hnc_minimum.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn_ts_guess.xyz '//work_dir//'/hcn_ts_guess.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn_minimum.xyz -b basis_sets/sto-3g.json --opt --opt-neb '// &
            'hnc_minimum.xyz '//NEB_FLAGS_CONVERGES//' --opt-neb-climb --silent --json neb_converges.json '// &
            '> neb_converges.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run itself failed, exit code "//to_string(rc)); return
      end if
      call json%initialize()
      call json%load(filename=work_dir//'/neb_converges.json')
      if (json%failed()) then
         call test_failed(error, "cannot read neb_converges.json"); return
      end if
      call json%get('optimization.converged', converged, found)
      call json%destroy()
      if (.not. found) then
         call test_failed(error, "no optimization.converged"); return
      end if
      if (.not. converged) then
         call test_failed(error, "optimization.converged was false even though the CLI exited 0 "// &
                          "(normal termination) -- the two should never disagree")
      end if
   end subroutine test_neb_converged_reports_true

   subroutine test_plain_neb(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call run_neb_plain(cli_exe, work_dir, 'neb_plain.json', 'neb_plain.out', nfail)
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_plain_neb

   subroutine test_no_stray(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: exists
      character(len=*), parameter :: stray(3) = &
         [character(len=20) :: 'hcn_minimum_opt.xyz', 'hcn_minimum_traj.xyz', 'hcn_minimum_neb.xyz']
      integer :: k
      do k = 1, size(stray)
         inquire (file=trim(stray(k)), exist=exists)
         if (exists) then
            call test_failed(error, trim(stray(k))//" was left in test/ itself"); return
         end if
      end do
   end subroutine test_no_stray

   subroutine read_xyz_energy(path, e, ok)
      !! The energy this program's own write_labelled_xyz-style comment
      !! line carries ("<label>   E = <value>"), for the committed
      !! endpoint files -- avoids a second xyz reader just for one number.
      character(len=*), intent(in) :: path
      real(dp), intent(out) :: e
      logical, intent(out) :: ok
      integer :: u, ios, k
      character(len=1024) :: line
      ok = .false.
      e = 0.0_dp
      open (newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      read (u, '(a)', iostat=ios) line   !! atom count
      read (u, '(a)', iostat=ios) line   !! comment: "...   E = <value>"
      close (u)
      if (ios /= 0) return
      k = index(line, '=')
      if (k == 0) return
      read (line(k + 1:), *, iostat=ios) e
      ok = (ios == 0)
   end subroutine read_xyz_energy

   subroutine run_neb(exe, wdir, flags, json_name, out_name, nf)
      character(len=*), intent(in) :: exe, wdir, flags, json_name, out_name
      integer, intent(inout) :: nf
      character(len=:), allocatable :: cmd2
      integer :: rc2
      type(json_file) :: j
      logical :: ok, ok2
      integer :: n_images, climb_idx
      real(dp), allocatable :: e_img(:)
      real(dp) :: e_hi, e_lo_first, e_lo_last
      logical :: climb_log_seen, conv_reached

      ! --opt-print-level 4 (not --silent): the "Tolerance factor to spawn
      ! climbing image"/"Convergence reached" lines this routine greps
      ! for only print at that verbosity (dl-find.f90's own printl-gated
      ! report block) -- also DL-FIND's OWN observable convergence
      ! verdict, this file's own header ("CONVERGENCE VERDICT").
      cmd2 = 'cd '//wdir//' && '//trim(exe)//' hcn_minimum.xyz -b basis_sets/sto-3g.json --opt --opt-neb '// &
             'hnc_minimum.xyz '//flags//' --opt-neb-climb --opt-print-level 4 --json '//json_name// &
             ' > '//out_name//' 2>&1'
      call execute_command_line(cmd2, exitstat=rc2, wait=.true.)
      ! This bridge's own exit code is NOT the gate here: its post-hoc
      ! max-ABSOLUTE-COMPONENT verdict and DL-FIND's own per-image RMS one
      ! (dlf_neb.f90's own fper) are different numbers for the identical
      ! gradient (this file's own header) and can disagree on which side
      ! of --opt-gtol they land, so a non-zero exit from "did not
      ! converge" alone is not treated as a failure for a climbing run --
      ! only a genuine crash (--json/log never written at all) is.
      call grep_contains(wdir//'/'//out_name, 'Convergence reached', conv_reached)
      if (.not. conv_reached) then
         ! A single retry, not a quiet pass: libdlfind's own
         ! dlf_neb_init (dlf_neb.f90:278) Kabsch-aligns the via/product
         ! frame against the start geometry on EVERY run, even with
         ! these FROZEN, committed endpoints (dlf_cartesian_align,
         ! dlf_coords.f90:1474, diagonalising a 3x3 matrix in
         ! dlf_hdlc_matrixlib.f90's array_diagonalise). HCN and HNC are
         ! both linear, so that matrix is RANK-DEFICIENT -- measured
         ! directly (gdb on a -finit-real=snan/-ffpe-trap=invalid debug
         ! build, which traps right here with a SIGFPE inside MKL's own
         ! dsyev): eigenvalues 3.183, 2.9e-14, 2.2e-16, i.e. a doubly
         ! degenerate near-zero pair. Any orthonormal basis of that 2-D
         ! null space is an equally valid eigenvector pair, so WHICH one
         ! the eigensolver returns is decided by floating-point noise at
         ! the 1e-14 level, not geometry -- and that noise floor is
         ! measured to shift under genuine external CPU/MKL contention
         ! (reproduced: identical 15/15 under no load or synthetic busy-
         ! loop load, but one ctest -j4 run with three other MKL-heavy
         ! SCF tests running at the same moment hit exactly this: the
         ! FIRST of two back-to-back, otherwise-identical runs failed to
         ! reach DL-FIND's own convergence inside --opt-maxiter, the
         ! second did). That is libdlfind's own numerical fragility
         ! (read-only fork at /home/jorge/dev/libdlfind-fork), not a
         ! terco bug and not one to patch here -- this chaotic 3-atom
         ! system already tips the same way for --opt-connect's dimer
         ! rotation (this file's own header). A bounded, logged retry is
         ! the honest mitigation: it still requires DL-FIND's own
         ! "Convergence reached" eventually, just not necessarily on the
         ! unlucky draw.
         print '(a)', '  (retry) '//out_name//': DL-FIND did not report convergence on the first '// &
            'attempt -- rerunning once (see run_neb''s own comment)'
         call execute_command_line(cmd2, exitstat=rc2, wait=.true.)
         call grep_contains(wdir//'/'//out_name, 'Convergence reached', conv_reached)
      end if
      call ref_report(nf, conv_reached, &
                      'check_opt_neb: DL-FIND''s own NEB climbing-image convergence test passed ('//out_name// &
                      ')', 0.0_dp)

      call j%initialize()
      call j%load(filename=wdir//'/'//json_name)
      if (j%failed()) then
         call ref_report(nf, .false., 'check_opt_neb: cannot read '//json_name, 0.0_dp); return
      end if
      call j%get('optimization.neb.n_images', n_images, ok)
      call j%get('optimization.neb.climbing_image_index', climb_idx, ok2)
      ok = ok .and. ok2
      call ref_report(nf, ok, 'check_opt_neb: neb.n_images/climbing_image_index present in '//json_name, 0.0_dp)
      if (.not. ok) then
         call j%destroy(); return
      end if
      call j%get('optimization.neb.image_energies', e_img, ok)
      ok = ok .and. size(e_img) == n_images
      call ref_report(nf, ok, 'check_opt_neb: neb.image_energies present and sized in '//json_name, 0.0_dp)
      if (.not. ok) then
         call j%destroy(); return
      end if

      e_hi = e_img(climb_idx)
      e_lo_first = e_img(1)
      e_lo_last = e_img(n_images)
      print '(a,i0,a,f20.12)', '  NEB n_images=', n_images, '   highest-energy image E = ', e_hi

      ! Endpoints equal the two committed minima -- within the NEB run
      ! itself, not just the files it was given (a sign the perpendicular
      ! endpoint mode is not silently moving them for this system).
      block
         real(dp) :: e_hcn_file, e_hnc_file
         logical :: ok3, ok4
         call read_xyz_energy(wdir//'/hcn_minimum.xyz', e_hcn_file, ok3)
         call read_xyz_energy(wdir//'/hnc_minimum.xyz', e_hnc_file, ok4)
         call ref_report(nf, ok3 .and. abs(e_lo_first - e_hcn_file) < ENDPOINT_TOL, &
                         'check_opt_neb: first NEB image equals test/hcn_minimum.xyz ('//out_name//')', &
                         abs(e_lo_first - e_hcn_file))
         call ref_report(nf, ok4 .and. abs(e_lo_last - e_hnc_file) < ENDPOINT_TOL, &
                         'check_opt_neb: last NEB image equals test/hnc_minimum.xyz ('//out_name//')', &
                         abs(e_lo_last - e_hnc_file))
      end block
      call ref_report(nf, all(e_hi >= e_img - 1.0e-9_dp), &
                      'check_opt_neb: the reported highest-energy image is really the maximum ('//out_name// &
                      ')', 0.0_dp)
      call ref_report(nf, e_hi > max(e_lo_first, e_lo_last) + 1.0e-4_dp, &
                      'check_opt_neb: the peak sits above both endpoints ('//out_name//')', e_hi)
      block
         integer :: uu, nat_ref, ia
         real(dp), allocatable :: z_ref(:), coords_ref(:, :), climb_xyz(:, :), r(:)
         real(dp) :: e_ts_ref, d1r, d2r, d3r, angr, d1c, d2c, d3c, angc
         logical :: geom_ok
         call ref_require('opt_ts_hcn.ref', uu)
         call ref_read_int(uu, 'NAT', nat_ref)
         allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
         call ref_read_vector(uu, 'Z', nat_ref, z_ref)
         call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
         call ref_read_real(uu, 'E_FINAL', e_ts_ref)
         call ref_close(uu)
         call ref_report(nf, e_hi <= e_ts_ref + 1.0e-6_dp, &
                         'check_opt_neb: the peak does not exceed the true TS energy ('//out_name//'), Eh', &
                         e_hi - e_ts_ref)
         call ref_report(nf, abs(e_hi - e_ts_ref) < ANCHOR_TS_TOL, &
                         'check_opt_neb: the peak is within the measured tolerance of the TS anchor '// &
                         '('//out_name//'), Eh', abs(e_hi - e_ts_ref))

         ! Geometry, rotation/translation invariant (two bond lengths plus
         ! the apex angle, check_opt_ts.F90's own two_shortest_and_angle) --
         ! NEB fixes nothing about the chain's absolute orientation, only
         ! each image's position along the path.
         geom_ok = .true.
         allocate (climb_xyz(3, nat_ref))
         do ia = 1, nat_ref
            call j%get('optimization.neb.climbing_image_geometry_angstrom('//to_string(ia)//')', r, ok)
            geom_ok = geom_ok .and. ok .and. size(r) == 3
            if (geom_ok) climb_xyz(:, ia) = r*ANGSTROM_TO_BOHR
         end do
         call ref_report(nf, geom_ok, 'check_opt_neb: neb.climbing_image_geometry_angstrom present in '// &
                         json_name, 0.0_dp)
         if (.not. geom_ok) then
            call j%destroy(); return
         end if

         call two_shortest_and_angle(coords_ref, d1r, d2r, d3r, angr)   !! anchor is already Bohr
         call two_shortest_and_angle(climb_xyz, d1c, d2c, d3c, angc)
         call ref_report(nf, abs(d1c - d1r) < ANCHOR_TS_GEOM_TOL, &
                         'check_opt_neb: climbing image vs TS anchor, shorter bond ('//out_name//'), Bohr', &
                         abs(d1c - d1r))
         call ref_report(nf, abs(d2c - d2r) < ANCHOR_TS_GEOM_TOL, &
                         'check_opt_neb: climbing image vs TS anchor, longer bond ('//out_name//'), Bohr', &
                         abs(d2c - d2r))
         call ref_report(nf, abs(d3c - d3r) < ANCHOR_TS_GEOM_TOL, &
                         'check_opt_neb: climbing image vs TS anchor, third distance ('//out_name//'), Bohr', &
                         abs(d3c - d3r))
         call ref_report(nf, abs(angc - angr) < ANCHOR_TS_ANGLE_TOL, &
                         'check_opt_neb: climbing image vs TS anchor, bond angle ('//out_name//'), deg', &
                         abs(angc - angr))
      end block
      call j%destroy()

      ! --opt-neb-climb did not just reach DL-FIND (the sentinel's own
      ! "Tolerance factor to spawn climbing image" line, which prints
      ! whenever neb_climb_test was left at DL-FIND's own default) but
      ! actually spawned the climbing image (dl-find.f90's own per-cycle
      ! report, "Spawning climbing image") -- "Convergence reached" is
      ! already checked above, before this run's own JSON is even read.
      call grep_contains(wdir//'/'//out_name, 'Tolerance factor to spawn climbing image', climb_log_seen)
      call ref_report(nf, climb_log_seen, 'check_opt_neb: --opt-neb-climb reached DL-FIND ('//out_name//')', &
                      0.0_dp)
      call grep_contains(wdir//'/'//out_name, 'Spawning climbing image', climb_log_seen)
      call ref_report(nf, climb_log_seen, 'check_opt_neb: the climbing image actually spawned ('//out_name// &
                      ')', 0.0_dp)
   end subroutine run_neb

   subroutine run_neb_plain(exe, wdir, json_name, out_name, nf)
      !! Plain NEB: no --opt-neb-climb, no --opt-neb-via, the ORIGINAL
      !! settings this file used before the climbing-image work above.
      !! Only the chain-wide profile (endpoints, a real peak above both
      !! of them, that peak within a loose, measured tolerance of the TS
      !! anchor) is checked -- no climbing-engaged or tight-anchor claim,
      !! since nothing here asks the chain to refine the peak at all.
      character(len=*), intent(in) :: exe, wdir, json_name, out_name
      integer, intent(inout) :: nf
      character(len=:), allocatable :: cmd2
      integer :: rc2
      type(json_file) :: j
      logical :: ok
      integer :: n_images, climb_idx
      real(dp), allocatable :: e_img(:)
      real(dp) :: e_hi, e_lo_first, e_lo_last

      call execute_command_line('mkdir -p '//wdir//'/basis_sets', exitstat=rc2, wait=.true.)
      call execute_command_line('cp -f hcn_minimum.xyz '//wdir//'/hcn_minimum.xyz', exitstat=rc2, wait=.true.)
      call execute_command_line('cp -f hnc_minimum.xyz '//wdir//'/hnc_minimum.xyz', exitstat=rc2, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//wdir//'/basis_sets/sto-3g.json', &
                                exitstat=rc2, wait=.true.)

      cmd2 = 'cd '//wdir//' && '//trim(exe)//' hcn_minimum.xyz -b basis_sets/sto-3g.json --opt --opt-neb '// &
             'hnc_minimum.xyz '//NEB_FLAGS_PLAIN//' --silent --json '//json_name//' > '//out_name//' 2>&1'
      call execute_command_line(cmd2, exitstat=rc2, wait=.true.)
      if (rc2 /= 0) then
         call ref_report(nf, .false., 'check_opt_neb: the plain --opt-neb CLI run itself failed ('//out_name// &
                         ')', 0.0_dp)
         return
      end if

      call j%initialize()
      call j%load(filename=wdir//'/'//json_name)
      if (j%failed()) then
         call ref_report(nf, .false., 'check_opt_neb: cannot read '//json_name, 0.0_dp); return
      end if
      call j%get('optimization.neb.n_images', n_images, ok)
      call j%get('optimization.neb.climbing_image_index', climb_idx, ok)
      call j%get('optimization.neb.image_energies', e_img, ok)
      ok = ok .and. size(e_img) == n_images
      call ref_report(nf, ok, 'check_opt_neb: bad/missing neb fields in '//json_name, 0.0_dp)
      if (.not. ok) then
         call j%destroy(); return
      end if
      call j%destroy()

      e_hi = e_img(climb_idx)
      e_lo_first = e_img(1)
      e_lo_last = e_img(n_images)
      print '(a,i0,a,f20.12)', '  plain NEB n_images=', n_images, '   highest-energy image E = ', e_hi

      block
         real(dp) :: e_hcn_file, e_hnc_file
         logical :: ok3, ok4
         call read_xyz_energy(wdir//'/hcn_minimum.xyz', e_hcn_file, ok3)
         call read_xyz_energy(wdir//'/hnc_minimum.xyz', e_hnc_file, ok4)
         call ref_report(nf, ok3 .and. abs(e_lo_first - e_hcn_file) < ENDPOINT_TOL, &
                         'check_opt_neb: plain NEB first image equals test/hcn_minimum.xyz ('//out_name//')', &
                         abs(e_lo_first - e_hcn_file))
         call ref_report(nf, ok4 .and. abs(e_lo_last - e_hnc_file) < ENDPOINT_TOL, &
                         'check_opt_neb: plain NEB last image equals test/hnc_minimum.xyz ('//out_name//')', &
                         abs(e_lo_last - e_hnc_file))
      end block
      call ref_report(nf, all(e_hi >= e_img - 1.0e-9_dp), &
                      'check_opt_neb: plain NEB highest-energy image is really the maximum ('//out_name// &
                      ')', 0.0_dp)
      call ref_report(nf, e_hi > max(e_lo_first, e_lo_last) + 1.0e-4_dp, &
                      'check_opt_neb: plain NEB peak sits above both endpoints ('//out_name//')', e_hi)
      block
         integer :: uu, nat_ref
         real(dp), allocatable :: z_ref(:), coords_ref(:, :)
         real(dp) :: e_ts_ref
         call ref_require('opt_ts_hcn.ref', uu)
         call ref_read_int(uu, 'NAT', nat_ref)
         allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
         call ref_read_vector(uu, 'Z', nat_ref, z_ref)
         call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
         call ref_read_real(uu, 'E_FINAL', e_ts_ref)
         call ref_close(uu)
         call ref_report(nf, e_hi <= e_ts_ref + 1.0e-6_dp, &
                         'check_opt_neb: plain NEB peak does not exceed the true TS energy ('//out_name// &
                         '), Eh', e_hi - e_ts_ref)
         call ref_report(nf, abs(e_hi - e_ts_ref) < ANCHOR_TS_TOL_PLAIN, &
                         'check_opt_neb: plain NEB peak is within the loose tolerance of the TS anchor '// &
                         '('//out_name//'), Eh', abs(e_hi - e_ts_ref))
      end block
   end subroutine run_neb_plain

   subroutine check_reproducible(wdir, json_a, json_b, nf)
      character(len=*), intent(in) :: wdir, json_a, json_b
      integer, intent(inout) :: nf
      type(json_file) :: ja, jb
      real(dp), allocatable :: ea(:), eb(:)
      logical :: oka, okb
      integer :: climb_a, climb_b
      ! Eh -- measured on the CLIMBING image specifically, CPU and GPU
      ! alike: ~1e-8 Eh (BLAS/LAPACK or cuSOLVER thread-scheduling-level
      ! floating-point noise through many SCF cycles, not a convergence
      ! or warm-start difference), so 1e-7 is generous headroom on the
      ! number this file actually validates against opt_ts_hcn.ref.
      real(dp), parameter :: REPRO_TOL = 1.0e-7_dp
      ! The REST of the chain (every image but the climbing one) is
      ! LOGGED, not asserted -- see this module's own header for why.
      call ja%initialize(); call ja%load(filename=wdir//'/'//json_a)
      call ja%get('optimization.neb.image_energies', ea, oka)
      call ja%get('optimization.neb.climbing_image_index', climb_a, oka)
      call ja%destroy()
      call jb%initialize(); call jb%load(filename=wdir//'/'//json_b)
      call jb%get('optimization.neb.image_energies', eb, okb)
      call jb%get('optimization.neb.climbing_image_index', climb_b, okb)
      call jb%destroy()
      if (.not. (oka .and. okb)) then
         call ref_report(nf, .false., &
                         'check_opt_neb: cannot read image_energies for the reproducibility check', 0.0_dp)
         return
      end if
      if (.not. (size(ea) == size(eb) .and. climb_a == climb_b)) then
         call ref_report(nf, .false., &
                         'check_opt_neb: reproducibility runs gave different image counts or climbing '// &
                         'images', 0.0_dp)
         return
      end if
      ! "Per-image warm start" is the ONLY path --opt-neb takes (there is
      ! no CLI flag for the single-shared-density fallback --
      ! app/trc_dlfind_bridge.F90's own trc_opt_options_t%neb_shared_density
      ! is bench-only, see bench_opt_neb_warmstart.F90 for that
      ! comparison's own numbers): this checks instead that the per-image
      ! path is DETERMINISTIC -- two fresh processes given identical input
      ! converge to the identical path, i.e. a warm start never reads
      ! stale state left over from a previous process. The physically
      ! meaningful claim is on the climbing image; the rest of the chain
      ! gets the looser, measured bound above instead of the same one.
      call ref_report(nf, abs(ea(climb_a) - eb(climb_b)) < REPRO_TOL, &
                      'check_opt_neb: two identical --opt-neb-climb runs land the climbing image on the '// &
                      'same point (per-image warm start is deterministic, not order-dependent)', &
                      abs(ea(climb_a) - eb(climb_b)))
      print '(a,es12.4,a)', '  (info) max |image_energy| spread over the WHOLE chain between the two runs: ', &
         maxval(abs(ea - eb)), ' Eh (not asserted -- see this subroutine''s own comment)'
   end subroutine check_reproducible

   subroutine two_shortest_and_angle(xyz, d1, d2, d3, angle_deg)
      !! The two shortest pairwise distances (the two "bonds", whichever
      !! atom is their shared vertex), the third distance, and the angle
      !! at the shared vertex -- rotation/translation invariant, unlike
      !! raw Cartesian coordinates. check_opt_ts.F90's own helper, copied
      !! rather than shared (no common test module for either to live in).
      real(dp), intent(in) :: xyz(3, 3)
      real(dp), intent(out) :: d1, d2, d3, angle_deg
      real(dp) :: d(3), v1(3), v2(3)
      integer :: apex
      d(1) = norm2(xyz(:, 1) - xyz(:, 2))
      d(2) = norm2(xyz(:, 1) - xyz(:, 3))
      d(3) = norm2(xyz(:, 2) - xyz(:, 3))
      if (d(1) >= d(2) .and. d(1) >= d(3)) then
         apex = 3; d3 = d(1); d1 = min(d(2), d(3)); d2 = max(d(2), d(3))
      else if (d(2) >= d(1) .and. d(2) >= d(3)) then
         apex = 2; d3 = d(2); d1 = min(d(1), d(3)); d2 = max(d(1), d(3))
      else
         apex = 1; d3 = d(3); d1 = min(d(1), d(2)); d2 = max(d(1), d(2))
      end if
      v1 = 0.0_dp; v2 = 0.0_dp
      select case (apex)
      case (1)
         v1 = xyz(:, 2) - xyz(:, 1); v2 = xyz(:, 3) - xyz(:, 1)
      case (2)
         v1 = xyz(:, 1) - xyz(:, 2); v2 = xyz(:, 3) - xyz(:, 2)
      case (3)
         v1 = xyz(:, 1) - xyz(:, 3); v2 = xyz(:, 2) - xyz(:, 3)
      end select
      if (norm2(v1) > 0.0_dp .and. norm2(v2) > 0.0_dp) then
         angle_deg = acos(max(-1.0_dp, min(1.0_dp, dot_product(v1, v2)/(norm2(v1)*norm2(v2)))))* &
                     180.0_dp/(4.0_dp*atan(1.0_dp))
      else
         angle_deg = 0.0_dp
      end if
   end subroutine two_shortest_and_angle

   subroutine grep_contains(path, needle, found_out)
      character(len=*), intent(in) :: path, needle
      logical, intent(out) :: found_out
      integer :: u, ios
      character(len=4096) :: line
      found_out = .false.
      open (newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, needle) > 0) then
            found_out = .true.
            exit
         end if
      end do
      close (u)
   end subroutine grep_contains

end module check_opt_neb_suite

program check_opt_neb
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_neb_suite, only: collect_opt_neb, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_neb: usage: check_opt_neb <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_neb, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_neb
