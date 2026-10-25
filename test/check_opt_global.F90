!
! --opt-algorithm stochastic, a fast fixed-population/fixed-cycle search
! on water/STO-3G: with --opt-seed set, reaches the same basin
! opt_water.ref's own minimum sits in, and running the SAME command
! twice reproduces the SAME search (DL-FIND's own trajectory through
! the population is bit-for-bit identical: --opt-seed's whole point;
! dlf_parallel_opt.f90 otherwise reseeds from the wall clock on every
! call, so two runs never agree without it -- app/trc_dlfind_bridge.F90's
! own header). Genetic is not covered here: it needs more cycles to
! reach its own po_tolerance_g test before it reports anything at all
! (dlf_parallel_opt.f90's dlf_genetic_endcycle only calls dlf_put_coords
! once that passes), so it is not the "fast" case this check wants;
! check_opt_refusals.F90 covers its CLI validation the same way
! stochastic's is covered there.
!
! --opt-global-tolerance-g matches --opt-gtol here on purpose:
! dlf_stoch_endcycle (dlf_parallel_opt.f90) stops the population search
! on its OWN po_tolerance_g test -- maxval(abs(igradient_best)) against
! it, read directly -- which is a SEPARATE DL-FIND parameter from
! --opt-gtol (this bridge's own convergence verdict, accepted_gmax
! against opts%gtol). Left at DL-FIND's own default (1e-3, dl-
! find.f90's own inline default), the population search was declaring
! itself converged and reporting its best candidate at a gradient
! inside 1e-3 but outside this run's own --opt-gtol 1e-4, which this
! bridge then (correctly) reported as "NOT converged" -- not a case of
! loosening a validation tolerance, but of giving DL-FIND's own
! internal stopping rule the SAME target this check asks the CLI for.
!
! TOLERANCES.
!
! REPRO_TOL = 1e-10 Eh/Angstrom, not exact equality: the two runs follow
! the IDENTICAL seeded trajectory through the population, but each
! geometry's own energy/gradient still comes from a fresh terco_cli
! process's own SCF, and on the GPU build a Fock build's reduction
! order is not guaranteed bit-identical run to run (thread-scheduling-
! dependent partial sums, same source of non-reproducibility as any
! other `do concurrent` reduction here) -- measured on a GPU build:
! |diff| 2.8e-14 Eh, 1.2e-13 Angstrom over 20 cycles' worth of
! accumulated round-off, four orders of headroom under this gate. The
! gfortran host build reproduces bit-for-bit (|diff| exactly 0).
!
! ENERGY_TOL = 1e-3 Eh is loose on purpose: population 8, 20 cycles is
! a SMALL search, and the exact trajectory through that population
! depends on the Fortran compiler's own RANDOM_NUMBER stream for the
! SAME seed array (the standard guarantees a seed reproduces a
! sequence only within one implementation, not across compilers) --
! a GPU build driving the same DL-FIND (itself CPU-only, serial)
! through nvfortran instead of gfortran is not guaranteed the
! IDENTICAL sequence this seed gives here, only its OWN reproducible
! one. Measured: final energy -74.965899825171 Eh (gfortran 15, host)
! and -74.965899686593 Eh (nvfortran 26.5, GPU) against opt_water.ref's
! -74.965901192300 Eh, |diff| 1.4e-6 and 1.5e-6 Eh -- three orders of
! headroom under the gate either way, left there for a different
! compiler's own stream to land somewhere else in the same basin.
!
! A test-drive suite.
!
module check_opt_global_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private
   public :: collect_opt_global, cli_exe, probe_json, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: probe_json
      !! The scratch .json base path, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the xyz/JSON outputs.
   real(dp), parameter :: ENERGY_TOL = 1.0e-3_dp   !! Eh -- see this file's own header
   real(dp), parameter :: REPRO_TOL = 1.0e-10_dp   !! Eh/Angstrom -- see this file's own header
   character(len=*), parameter :: SEED_ARGS = '--opt-algorithm stochastic --opt-seed 42 '// &
      '--opt-global-pop-size 8 --opt-global-maxcycle 20 --opt-global-tolerance-g 1e-4 --opt-gtol 1e-4'

   ! --- state shared between the two CLI runs and the checks on them ------
   real(dp) :: e1 = 0.0_dp, e2 = 0.0_dp
   real(dp), allocatable :: row1(:), row2(:)

contains

   subroutine collect_opt_global(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt-algorithm stochastic: --opt-seed reproduces the same run twice", &
                      test_global_reproducible), &
         new_unittest("--opt-algorithm stochastic: final energy matches the opt_water anchor", &
                      test_global_vs_anchor), &
         new_unittest("--opt-algorithm stochastic: no stray xyz output in test/ itself", test_no_stray) &
         ]
   end subroutine collect_opt_global

   subroutine test_global_reproducible(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe1, probe2
      integer :: rc
      logical :: found1, found2
      type(json_file) :: j1, j2

      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)

      probe1 = probe_json//'.1'
      probe2 = probe_json//'.2'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
            SEED_ARGS//' --silent --json '//probe1//' > '//probe1//'.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "run 1 (the CLI itself) failed, exit code "//to_string(rc)); return
      end if

      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
            SEED_ARGS//' --silent --json '//probe2//' > '//probe2//'.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "run 2 (the CLI itself) failed, exit code "//to_string(rc)); return
      end if

      call j1%initialize(); call j1%load(filename=probe1)
      if (j1%failed()) then
         call test_failed(error, "cannot read run 1's JSON"); return
      end if
      call j1%get('optimization.final_energy', e1, found1)
      if (.not. found1) then
         call j1%destroy(); call test_failed(error, "run 1: no optimization.final_energy"); return
      end if
      call j1%get('optimization.final_geometry_angstrom(1)', row1, found1)
      if (.not. found1) then
         call j1%destroy(); call test_failed(error, "run 1: no final_geometry_angstrom(1)"); return
      end if
      call j1%destroy()

      call j2%initialize(); call j2%load(filename=probe2)
      if (j2%failed()) then
         call test_failed(error, "cannot read run 2's JSON"); return
      end if
      call j2%get('optimization.final_energy', e2, found2)
      if (.not. found2) then
         call j2%destroy(); call test_failed(error, "run 2: no optimization.final_energy"); return
      end if
      call j2%get('optimization.final_geometry_angstrom(1)', row2, found2)
      if (.not. found2) then
         call j2%destroy(); call test_failed(error, "run 2: no final_geometry_angstrom(1)"); return
      end if
      call j2%destroy()

      print '(a,f20.12)', '  run 1 final energy = ', e1
      print '(a,f20.12)', '  run 2 final energy = ', e2

      if (abs(e1 - e2) >= REPRO_TOL) then
         call test_failed(error, "run 1 and run 2 final_energy do not reproduce, |diff|="// &
                          to_string(abs(e1 - e2))//" Eh"); return
      end if
      if (maxval(abs(row1 - row2)) >= REPRO_TOL) then
         call test_failed(error, "run 1 and run 2 final_geometry_angstrom(1) do not reproduce, |diff|="// &
                          to_string(maxval(abs(row1 - row2)))//" Angstrom")
      end if
   end subroutine test_global_reproducible

   !
   ! Depends on test_global_reproducible having already run (e1, above) --
   ! test-drive runs a module's tests in the order collect_opt_global
   ! lists them, not in parallel (run_testsuite's own parallel=.false. in
   ! the program below).
   !
   subroutine test_global_vs_anchor(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: uu, nat_ref
      real(dp), allocatable :: z_ref(:), coords_ref(:, :)
      real(dp) :: e_ref, gmax_ref
      call ref_require('opt_water.ref', uu)
      call ref_read_int(uu, 'NAT', nat_ref)
      allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
      call ref_read_vector(uu, 'Z', nat_ref, z_ref)
      call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
      call ref_read_real(uu, 'E_FINAL', e_ref)
      call ref_read_real(uu, 'GMAX_FINAL', gmax_ref)
      call ref_close(uu)
      if (abs(e1 - e_ref) >= ENERGY_TOL) then
         call test_failed(error, "final energy vs opt_water anchor, |diff|="//to_string(abs(e1 - e_ref))//" Eh")
      end if
   end subroutine test_global_vs_anchor

   !
   ! best.xyz/best_active.xyz (dlf_parallel_opt.f90's own printf>=2 output,
   ! needed for this bridge to ever see the result at all -- app/
   ! trc_dlfind_bridge.F90's own header) must not have landed in test/
   ! itself (this program's own WORKING_DIRECTORY, ctest's convention).
   !
   subroutine test_no_stray(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: stray(4) = &
         [character(len=16) :: 'best.xyz', 'best_active.xyz', 'water_opt.xyz', 'water_traj.xyz']
      integer :: k
      logical :: exists
      do k = 1, size(stray)
         inquire (file=trim(stray(k)), exist=exists)
         if (exists) then
            call test_failed(error, trim(stray(k))//" was left in test/ itself"); return
         end if
      end do
   end subroutine test_no_stray

end module check_opt_global_suite

program check_opt_global
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_global_suite, only: collect_opt_global, cli_exe, probe_json, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); probe_json = trim(arg)
   call get_command_argument(3, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_global: usage: check_opt_global <terco CLI> <scratch .json> "// &
         "<scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_global, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_global
