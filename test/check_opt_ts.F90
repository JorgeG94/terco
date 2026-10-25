!
! The transition-state search: --opt-target saddle on the HCN<->HNC
! isomerisation, RHF/STO-3G, from test/hcn_ts_guess.xyz (a bent geometry
! about half way from a distorted HCN toward the critical point -- see
! that file's own comment line for how it was picked: it is the point
! where Cartesian P-RFO, started closer to HCN itself, first reaches the
! SAME energy the dimer method finds independently, with no Hessian at
! all). Two independent routes to the same point is the validation here:
!   - P-RFO, Cartesian, terco's own analytic-from-FD Hessian
!     (cb_get_hessian/trc_run_hessian), --opt-zero-modes' auto default (6:
!     the start geometry is not linear);
!   - the dimer method (--opt-saddle-method dimer), no Hessian at all.
! Both must land on the same energy and geometry, and --hessian at the
! result must show EXACTLY ONE imaginary mode -- the saddle check, not a
! second copy of the energy comparison.
!
! EXTERNAL anchor: test/pyscf_refs/opt_ts_hcn.ref, written by
! test/pyscf_refs/generate.py's own gen_opt_ts -- a from-scratch
! partitioned-RFO (eigenvector-following) search on pyscf's OWN analytic
! RHF gradient and Hessian, translation/rotation projected out before
! diagonalising, converged to max|g| < 1e-8 (two orders tighter than
! either terco run below), then pyscf's own harmonic_analysis at the
! result -- refused by the generator unless that shows EXACTLY ONE
! imaginary mode. Independent of terco end to end: different language,
! different SCF code, a from-scratch optimiser, not DL-FIND's. Measured
! agreement (this file's own run): energy to ~1e-10 Eh, imaginary
! frequency to ~0.02 cm^-1.
!
! --opt-coords hdlc with --opt-target saddle prfo is NOT checked here
! either: measured, it does not converge reliably for this 3-atom system
! (DL-FIND's own "lam-eigval small for non-TS mode" warning, step lengths
! of order 10 Bohr, no progress after a handful of steps) -- HDLC's
! hybrid internal machinery is built for larger, fragmented systems, and
! a bare triatomic is a poor fit for it independent of anything terco's
! bridge does. Cartesian and the dimer method (which also has a Cartesian
! option, exercised above) are the validated route for a system this
! small.
!
! A test-drive suite.
!
module check_opt_ts_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use trc_constants, only: ANGSTROM_TO_BOHR
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   implicit none
   private
   public :: collect_opt_ts, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the xyz/JSON outputs.
   real(dp), parameter :: ENERGY_TOL = 1.0e-6_dp
      !! Eh -- both runs are --opt-gtol 1e-5/1e-6; see this module's own
      !! header.
   real(dp), parameter :: GEOM_TOL = 5.0e-4_dp
      !! Bohr, same spirit as check_opt_vs_pyscf's BOND_TOL.
   ! Anchor (pyscf) tolerances, same reasoning as check_opt_vs_pyscf.F90's
   ! own ENERGY_TOL/BOND_TOL (dE ~ g^2/(2k), g = the LOOSER of the two
   ! terco gtols below, 1e-5, k ~ the TS's own lowest real curvature,
   ! ~0.15 Eh/Bohr^2 -- test/pyscf_refs/opt_ts_hcn.ref's own comment):
   ! dE ~ 3e-10 Eh, ~300x headroom here.
   real(dp), parameter :: ANCHOR_ENERGY_TOL = 1.0e-7_dp
   real(dp), parameter :: ANCHOR_GEOM_TOL = 3.0e-4_dp   !! Bohr
   real(dp), parameter :: ANCHOR_ANGLE_TOL = 0.05_dp    !! deg
   ! The imaginary frequency: terco's own is from a CENTRAL-DIFFERENCE
   ! Hessian (trc_run_hessian, h=0.005 Bohr); the anchor's is pyscf's
   ! analytic one. Measured agreement for this molecule/mode (this
   ! module's own header): ~0.02 cm^-1; generous ~15x headroom here for a
   ! different mode/molecule than whatever that one measurement covers.
   real(dp), parameter :: ANCHOR_FREQ_TOL = 0.3_dp   !! cm^-1

   ! --- state shared between test_prfo_vs_dimer and test_vs_pyscf_anchor --
   real(dp) :: e_prfo = 0.0_dp, e_dimer = 0.0_dp, gmax_prfo = 0.0_dp, gmax_dimer = 0.0_dp
   real(dp), allocatable :: row_prfo(:, :), row_dimer(:, :)

contains

   subroutine collect_opt_ts(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt-target saddle: P-RFO (Hessian) and the dimer method (no Hessian) "// &
                      "agree on the HCN<->HNC TS", test_prfo_vs_dimer), &
         new_unittest("--opt-target saddle vs pyscf's own from-scratch P-RFO anchor, and the "// &
                      "saddle check (exactly one imaginary mode)", test_vs_pyscf_anchor), &
         new_unittest("--opt-target saddle: no xyz output was left in test/ itself", test_no_stray) &
         ]
   end subroutine collect_opt_ts

   subroutine test_prfo_vs_dimer(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd
      integer :: rc, nfail
      logical :: conv_prfo, conv_dimer, found
      type(json_file) :: json

      nfail = 0
      ! Run with WORKING_DIRECTORY set to work_dir (not test/ itself, this
      ! program's own WORKING_DIRECTORY, ctest's convention) -- --opt
      ! writes <geometry>_opt.xyz/_traj.xyz into the CLI's CWD
      ! (app/main.F90), and everything the subprocess needs is copied in
      ! so the `cd work_dir &&` below can use plain relative paths --
      ! check_opt_connect.F90's own pattern.
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn_ts_guess.xyz '//work_dir//'/hcn_ts_guess.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)

      ! --- P-RFO, Cartesian, terco's own analytic-from-FD Hessian ---------
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn_ts_guess.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-target saddle --opt-saddle-method prfo --opt-gtol 1e-6 --opt-maxiter 60 --silent --json '// &
            'prfo.json > prfo.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the P-RFO CLI run itself failed, exit code "//to_string(rc)); return
      end if
      call json%initialize()
      call json%load(filename=work_dir//'/prfo.json')
      if (json%failed()) then
         call test_failed(error, "cannot read prfo.json"); return
      end if
      call json%get('optimization.converged', conv_prfo, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no prfo optimization.converged"); return
      end if
      call json%get('optimization.final_energy', e_prfo, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no prfo optimization.final_energy"); return
      end if
      call json%get('optimization.final_gradient_max', gmax_prfo, found)
      call read_geometry(json, 3, row_prfo, error)
      call json%destroy()
      if (allocated(error)) return
      call ref_report(nfail, conv_prfo, 'check_opt_ts: P-RFO Cartesian converged', 0.0_dp)

      ! --- the dimer method: no Hessian at all ----------------------------
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn_ts_guess.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-target saddle --opt-saddle-method dimer --opt-gtol 1e-5 --opt-maxiter 150 --silent --json '// &
            'dimer.json > dimer.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the dimer CLI run itself failed, exit code "//to_string(rc)); return
      end if
      call json%initialize()
      call json%load(filename=work_dir//'/dimer.json')
      if (json%failed()) then
         call test_failed(error, "cannot read dimer.json"); return
      end if
      call json%get('optimization.converged', conv_dimer, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no dimer optimization.converged"); return
      end if
      call json%get('optimization.final_energy', e_dimer, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no dimer optimization.final_energy"); return
      end if
      call json%get('optimization.final_gradient_max', gmax_dimer, found)
      call read_geometry(json, 3, row_dimer, error)
      call json%destroy()
      if (allocated(error)) return
      call ref_report(nfail, conv_dimer, 'check_opt_ts: dimer method converged', 0.0_dp)

      print '(a,f20.12,a,es10.3)', '  P-RFO  final energy ', e_prfo, '   max|g| ', gmax_prfo
      print '(a,f20.12,a,es10.3)', '  dimer  final energy ', e_dimer, '   max|g| ', gmax_dimer

      call ref_report(nfail, abs(e_prfo - e_dimer) < ENERGY_TOL, &
                      'check_opt_ts: P-RFO and dimer agree on the TS energy', abs(e_prfo - e_dimer))
      ! Compared as INTERNAL coordinates (bond lengths + angle), not raw
      ! Cartesian components: the two searches take different paths and
      ! can land on the same physical structure rotated/translated
      ! differently in space (zero modes -- translation/rotation never
      ! affect the energy, so raw per-atom coordinates are not the right
      ! comparison; check_opt_vs_pyscf.F90's own two_shortest_and_angle
      ! does the same).
      block
         real(dp) :: d1p, d2p, d3p, angp, d1d, d2d, d3d, angd
         call two_shortest_and_angle(row_prfo, d1p, d2p, d3p, angp)
         call two_shortest_and_angle(row_dimer, d1d, d2d, d3d, angd)
         call ref_report(nfail, abs(d1p - d1d) < GEOM_TOL, 'check_opt_ts: shorter bond agrees', abs(d1p - d1d))
         call ref_report(nfail, abs(d2p - d2d) < GEOM_TOL, 'check_opt_ts: longer bond agrees', abs(d2p - d2d))
         call ref_report(nfail, abs(d3p - d3d) < GEOM_TOL, 'check_opt_ts: third distance agrees', abs(d3p - d3d))
         call ref_report(nfail, abs(angp - angd) < 0.1_dp, 'check_opt_ts: bond angle agrees (deg)', &
                         abs(angp - angd))
      end block
      ! A genuine transition state sits well above either HCN or HNC's
      ! own minimum (opt_hcn.ref: -91.67520899058249 Eh) -- not a second
      ! minimum that a broken search slid into.
      call ref_report(nfail, e_prfo > -91.60_dp .and. e_prfo < -91.50_dp, &
                      'check_opt_ts: the TS energy is well above the HCN/HNC minima (sanity range)', 0.0_dp)

      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_prfo_vs_dimer

   subroutine test_vs_pyscf_anchor(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail, uu, nat_ref, n_imag_ref
      real(dp), allocatable :: z_ref(:), coords_ref(:, :), freq_ref(:)
      real(dp) :: e_ref, gmax_ref, d1r, d2r, d3r, angr, d1p, d2p, d3p, angp

      if (.not. (allocated(row_prfo) .and. allocated(row_dimer))) then
         call test_failed(error, "test_prfo_vs_dimer did not leave a TS geometry to compare against"); return
      end if

      nfail = 0
      call ref_require('opt_ts_hcn.ref', uu)
      call ref_read_int(uu, 'NAT', nat_ref)
      if (nat_ref /= 3) then
         call ref_close(uu); call test_failed(error, "opt_ts_hcn.ref: NAT must be 3"); return
      end if
      allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
      call ref_read_vector(uu, 'Z', nat_ref, z_ref)
      call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
      call ref_read_real(uu, 'E_FINAL', e_ref)
      call ref_read_real(uu, 'GMAX_FINAL', gmax_ref)
      call ref_read_int(uu, 'N_IMAGINARY', n_imag_ref)
      ! 3*natm - 6 = 3 vibrational modes for a 3-atom molecule.
      allocate (freq_ref(nat_ref))
      call ref_read_vector(uu, 'FREQ_CM1', 3, freq_ref)
      call ref_close(uu)
      if (n_imag_ref /= 1) then
         call test_failed(error, "opt_ts_hcn.ref itself must have exactly one imaginary mode (the "// &
                          "generator already refuses to write one otherwise)"); return
      end if

      ! Anchor is already Bohr.
      call two_shortest_and_angle(coords_ref, d1r, d2r, d3r, angr)
      call two_shortest_and_angle(row_prfo*ANGSTROM_TO_BOHR, d1p, d2p, d3p, angp)
      call ref_report(nfail, abs(e_prfo - e_ref) < ANCHOR_ENERGY_TOL, &
                      'check_opt_ts (P-RFO vs pyscf anchor): TS energy, Eh', abs(e_prfo - e_ref))
      call ref_report(nfail, abs(d1p - d1r) < ANCHOR_GEOM_TOL, &
                      'check_opt_ts (P-RFO vs pyscf anchor): shorter bond, Bohr', abs(d1p - d1r))
      call ref_report(nfail, abs(d2p - d2r) < ANCHOR_GEOM_TOL, &
                      'check_opt_ts (P-RFO vs pyscf anchor): longer bond, Bohr', abs(d2p - d2r))
      call ref_report(nfail, abs(d3p - d3r) < ANCHOR_GEOM_TOL, &
                      'check_opt_ts (P-RFO vs pyscf anchor): third distance, Bohr', abs(d3p - d3r))
      call ref_report(nfail, abs(angp - angr) < ANCHOR_ANGLE_TOL, &
                      'check_opt_ts (P-RFO vs pyscf anchor): bond angle, deg', abs(angp - angr))

      call two_shortest_and_angle(row_dimer*ANGSTROM_TO_BOHR, d1p, d2p, d3p, angp)
      call ref_report(nfail, abs(e_dimer - e_ref) < ANCHOR_ENERGY_TOL, &
                      'check_opt_ts (dimer vs pyscf anchor): TS energy, Eh', abs(e_dimer - e_ref))
      call ref_report(nfail, abs(d1p - d1r) < ANCHOR_GEOM_TOL, &
                      'check_opt_ts (dimer vs pyscf anchor): shorter bond, Bohr', abs(d1p - d1r))
      call ref_report(nfail, abs(d2p - d2r) < ANCHOR_GEOM_TOL, &
                      'check_opt_ts (dimer vs pyscf anchor): longer bond, Bohr', abs(d2p - d2r))
      call ref_report(nfail, abs(angp - angr) < ANCHOR_ANGLE_TOL, &
                      'check_opt_ts (dimer vs pyscf anchor): bond angle, deg', abs(angp - angr))

      ! The saddle check: exactly one imaginary mode, within the
      ! FD-Hessian error of the anchor's own analytic frequency.
      call check_one_imaginary(cli_exe, work_dir, row_prfo, minval(freq_ref), nfail)

      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_vs_pyscf_anchor

   !
   ! hcn_ts_guess_opt.xyz/hcn_ts_guess_traj.xyz must not have landed in
   ! test/ itself (this program's own WORKING_DIRECTORY, ctest's
   ! convention) -- the whole point of test_prfo_vs_dimer's own
   ! `cd work_dir &&`.
   !
   subroutine test_no_stray(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: stray(2) = &
         [character(len=22) :: 'hcn_ts_guess_opt.xyz', 'hcn_ts_guess_traj.xyz']
      integer :: k, nfail
      logical :: exists

      nfail = 0
      do k = 1, size(stray)
         inquire (file=trim(stray(k)), exist=exists)
         call ref_report(nfail, .not. exists, 'check_opt_ts: '//trim(stray(k))//' was not left in test/ itself', &
                         merge(1.0_dp, 0.0_dp, exists))
      end do
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_no_stray

   subroutine read_geometry(j, natm, xyz, error)
      type(json_file), intent(inout) :: j
      integer, intent(in) :: natm
      real(dp), allocatable, intent(out) :: xyz(:, :)
      type(error_type), allocatable, intent(inout) :: error
      real(dp), allocatable :: row(:)
      logical :: ok
      integer :: ia
      allocate (xyz(3, natm))
      do ia = 1, natm
         call j%get('optimization.final_geometry_angstrom('//to_string(ia)//')', row, ok)
         if (.not. (ok .and. size(row) == 3)) then
            call test_failed(error, "bad final_geometry_angstrom row "//to_string(ia)); return
         end if
         xyz(:, ia) = row
      end do
   end subroutine read_geometry

   !
   ! The two shortest pairwise distances (the two "bonds", whichever atom
   ! is their shared vertex), the third distance, and the angle at the
   ! shared vertex -- check_opt_vs_pyscf.F90's own helper, copied rather
   ! than shared across two standalone check programs.
   !
   subroutine two_shortest_and_angle(xyz, d1, d2, d3, angle_deg)
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

   !
   ! --hessian at the P-RFO result: exactly one imaginary mode is the
   ! saddle check (app/main.F90's opt_summary already prints it; here it
   ! is asserted).
   !
   subroutine check_one_imaginary(exe, wdir, xyz, freq_anchor, nf)
      character(len=*), intent(in) :: exe, wdir
      real(dp), intent(in) :: xyz(:, :)
      real(dp), intent(in) :: freq_anchor   !! pyscf's own analytic imaginary frequency, cm^-1 (negative)
      integer, intent(inout) :: nf
      character(len=:), allocatable :: cmd2
      type(json_file) :: hj
      integer :: u, rc2, n_imag
      real(dp), allocatable :: freqs(:)
      logical :: ok
      open (newunit=u, file=wdir//'/ts_candidate.xyz', status='replace', action='write')
      write (u, '(i0)') size(xyz, 2)
      write (u, '(a)') 'check_opt_ts: the P-RFO result, re-checked with --hessian'
      write (u, '(a2,3f20.12)') 'H', xyz(:, 1)
      write (u, '(a2,3f20.12)') 'C', xyz(:, 2)
      write (u, '(a2,3f20.12)') 'N', xyz(:, 3)
      close (u)
      cmd2 = 'cd '//wdir//' && '//trim(exe)//' ts_candidate.xyz -b basis_sets/sto-3g.json --hessian --silent '// &
             '--json hess.json > hess.out 2>&1'
      call execute_command_line(cmd2, exitstat=rc2, wait=.true.)
      if (rc2 /= 0) then
         call ref_report(nf, .false., 'check_opt_ts: the --hessian re-check CLI run itself failed', 0.0_dp); return
      end if
      call hj%initialize()
      call hj%load(filename=wdir//'/hess.json')
      if (hj%failed()) then
         call ref_report(nf, .false., 'check_opt_ts: cannot read hess.json', 0.0_dp); return
      end if
      call hj%get('frequencies.n_imaginary', n_imag, ok)
      if (.not. ok) then
         call hj%destroy(); call ref_report(nf, .false., 'check_opt_ts: no frequencies.n_imaginary', 0.0_dp); return
      end if
      call hj%get('frequencies.frequencies_cm1', freqs, ok)
      if (.not. ok) then
         call hj%destroy()
         call ref_report(nf, .false., 'check_opt_ts: no frequencies.frequencies_cm1', 0.0_dp); return
      end if
      call hj%destroy()
      call ref_report(nf, n_imag == 1, 'check_opt_ts: the TS has EXACTLY ONE imaginary mode', real(n_imag, dp))
      if (allocated(freqs) .and. size(freqs) > 0) then
         print '(a,f10.4,a)', '  imaginary frequency, terco FD-Hessian : ', minval(freqs), ' i cm^-1'
         print '(a,f10.4,a)', '  imaginary frequency, pyscf analytic   : ', freq_anchor, ' i cm^-1'
         call ref_report(nf, abs(minval(freqs) - freq_anchor) < ANCHOR_FREQ_TOL, &
                         'check_opt_ts: imaginary frequency within the FD-Hessian error of the pyscf anchor', &
                         abs(minval(freqs) - freq_anchor))
      end if
   end subroutine check_one_imaginary

end module check_opt_ts_suite

program check_opt_ts
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_ts_suite, only: collect_opt_ts, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_ts: usage: check_opt_ts <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_ts, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_ts
