!
! --opt-algorithm stochastic, a fast fixed-population/fixed-cycle search
! on water/STO-3G: with --opt-seed set, reaches the same basin
! opt_water.ref's own minimum sits in, and running the SAME command
! twice gives IDENTICAL results -- the whole point of --opt-seed
! (dlf_parallel_opt.f90 otherwise reseeds from the wall clock on every
! call, so two runs never agree without it; app/trc_dlfind_bridge.F90's
! own header). Genetic is not covered here: it needs more cycles to
! reach its own po_tolerance_g test before it reports anything at all
! (dlf_parallel_opt.f90's dlf_genetic_endcycle only calls dlf_put_coords
! once that passes), so it is not the "fast" case this check wants;
! check_opt_refusals.F90 covers its CLI validation the same way
! stochastic's is covered there.
!
! TOLERANCE. ENERGY_TOL = 1e-3 Eh is loose on purpose: population 8,
! 20 cycles is a SMALL search, and the exact trajectory through that
! population depends on the Fortran compiler's own RANDOM_NUMBER
! stream for the SAME seed array (the standard guarantees a seed
! reproduces a sequence only within one implementation, not across
! compilers) -- a GPU build driving the same DL-FIND (itself CPU-only,
! serial) through nvfortran instead of gfortran is not guaranteed the
! IDENTICAL sequence this seed gives here, only its OWN reproducible
! one. Measured with gfortran 15 (host): final energy -74.965899825171
! Eh against opt_water.ref's -74.965901192300 Eh, |diff| 1.4e-6 Eh --
! three orders of headroom under the gate, left there for a different
! compiler's own stream to land somewhere else in the same basin.
!
program check_opt_global
   use pic_types, only: dp
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   use json_module, only: json_file
   implicit none

   character(len=:), allocatable :: cli_exe, probe_json, work_dir, cmd, probe1, probe2
   character(len=1024) :: arg
   integer :: rc, nfail
   real(dp) :: e1, e2, e_ref, gmax_ref
   real(dp), allocatable :: row1(:), row2(:)
   logical :: found1, found2
   type(json_file) :: j1, j2
   real(dp), parameter :: ENERGY_TOL = 1.0e-3_dp   !! Eh -- see this file's own header
   character(len=*), parameter :: SEED_ARGS = '--opt-algorithm stochastic --opt-seed 42 '// &
      '--opt-global-pop-size 8 --opt-global-maxcycle 20 --opt-gtol 1e-4'

   nfail = 0
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); probe_json = trim(arg)
   call get_command_argument(3, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0 .or. len(work_dir) == 0) then
      print '(a)', 'check_opt_global: usage: check_opt_global <terco CLI> <scratch .json> <scratch dir>'
      stop 1
   end if
   call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                             exitstat=rc, wait=.true.)

   probe1 = probe_json//'.1'
   probe2 = probe_json//'.2'
   cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
         SEED_ARGS//' --silent --json '//probe1//' > '//probe1//'.out 2>&1'
   call execute_command_line(cmd, exitstat=rc, wait=.true.)
   call expect(rc == 0, 'check_opt_global: run 1 (the CLI itself) failed')

   cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
         SEED_ARGS//' --silent --json '//probe2//' > '//probe2//'.out 2>&1'
   call execute_command_line(cmd, exitstat=rc, wait=.true.)
   call expect(rc == 0, 'check_opt_global: run 2 (the CLI itself) failed')

   call j1%initialize(); call j1%load(filename=probe1)
   call expect(.not. j1%failed(), 'check_opt_global: cannot read run 1''s JSON')
   call j1%get('optimization.final_energy', e1, found1)
   call expect(found1, 'check_opt_global: run 1: no optimization.final_energy')
   call j1%get('optimization.final_geometry_angstrom(1)', row1, found1)
   call expect(found1, 'check_opt_global: run 1: no final_geometry_angstrom(1)')
   call j1%destroy()

   call j2%initialize(); call j2%load(filename=probe2)
   call expect(.not. j2%failed(), 'check_opt_global: cannot read run 2''s JSON')
   call j2%get('optimization.final_energy', e2, found2)
   call expect(found2, 'check_opt_global: run 2: no optimization.final_energy')
   call j2%get('optimization.final_geometry_angstrom(1)', row2, found2)
   call expect(found2, 'check_opt_global: run 2: no final_geometry_angstrom(1)')
   call j2%destroy()

   print '(a,f20.12)', '  run 1 final energy = ', e1
   print '(a,f20.12)', '  run 2 final energy = ', e2

   call report(nfail, e1 == e2, 'check_opt_global: --opt-seed: run 1 and run 2 final_energy are IDENTICAL', &
              abs(e1 - e2))
   call report(nfail, all(row1 == row2), &
              'check_opt_global: --opt-seed: run 1 and run 2 final_geometry_angstrom(1) are IDENTICAL', &
              maxval(abs(row1 - row2)))

   block
      integer :: uu, nat_ref
      real(dp), allocatable :: z_ref(:), coords_ref(:, :)
      call ref_require('opt_water.ref', uu)
      call ref_read_int(uu, 'NAT', nat_ref)
      allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
      call ref_read_vector(uu, 'Z', nat_ref, z_ref)
      call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
      call ref_read_real(uu, 'E_FINAL', e_ref)
      call ref_read_real(uu, 'GMAX_FINAL', gmax_ref)
      call ref_close(uu)
   end block
   call ref_report(nfail, abs(e1 - e_ref) < ENERGY_TOL, &
                  'check_opt_global: final energy vs opt_water anchor, Eh', abs(e1 - e_ref))

   ! best.xyz/best_active.xyz (dlf_parallel_opt.f90's own printf>=2 output,
   ! needed for this bridge to ever see the result at all -- app/
   ! trc_dlfind_bridge.F90's own header) must not have landed in test/
   ! itself (this program's own WORKING_DIRECTORY, ctest's convention).
   block
      character(len=*), parameter :: stray(4) = &
         [character(len=16) :: 'best.xyz', 'best_active.xyz', 'water_opt.xyz', 'water_traj.xyz']
      integer :: k
      logical :: exists
      do k = 1, size(stray)
         inquire (file=trim(stray(k)), exist=exists)
         call report(nfail, .not. exists, 'check_opt_global: '//trim(stray(k))//' was not left in test/ itself')
      end do
   end block

   if (nfail > 0) then
      print '(a,i0,a)', 'check_opt_global: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_opt_global: RESULT: PASS'

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', what; stop 1
      end if
   end subroutine expect

   subroutine report(nf, cond, what, val)
      integer, intent(inout) :: nf
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      real(dp), intent(in), optional :: val
      if (cond) then
         print '(a,a)', '  PASS  ', what
      else
         if (present(val)) then
            print '(a,a,a,es12.4)', '  FAIL  ', what, '   value=', val
         else
            print '(a,a)', '  FAIL  ', what
         end if
         nf = nf + 1
      end if
   end subroutine report

end program check_opt_global
