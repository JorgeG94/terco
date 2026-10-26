!
! SCF iterations per NEB image evaluation: a per-image warm-start density
! (trc_opt_options_t%neb -- the default) against ONE density shared across
! every image in turn (neb_shared_density, bench-only -- app/
! trc_dlfind_bridge.F90's own header) -- this PR's own warm-start
! measurement for --opt-neb, the per-image analogue of
! bench_opt_warmstart.F90's single-state one. Informational only: prints a
! table, asserts nothing, not registered as a ctest. Run by hand:
!
!   build-cpu/test/bench_opt_neb_warmstart
!
program bench_opt_neb_warmstart
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_null_ptr, c_null_char
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_method, &
                       trc_set_convergence, trc_set_guess, TRC_GUESS_SAD
   use trc_api, only: trc_bind_device
   use trc_dlfind_bridge, only: trc_opt_options_t, trc_opt_result_t, trc_dlfind_optimize, &
                                TRC_OPT_GUESS_DENSITY, TRC_OPT_NEB_PERP
   use trc_test_basis, only: read_xyz
   implicit none

   integer :: dev
   integer :: nat
   integer, allocatable :: z(:)
   real(dp), allocatable :: r(:, :), r2(:, :)

   dev = trc_bind_device(0)

   call read_xyz('hcn_minimum.xyz', nat, z, r)
   call read_xyz('hnc_minimum.xyz', nat, z, r2)
   ! The same two FROZEN endpoints check_opt_neb.F90 uses (that program's
   ! own header explains why they are committed files, not a fresh
   ! --opt-connect run) -- this bench expects to be run from test/ itself,
   ! or with those two files copied next to it.

   print '(a)', 'warm start          images  steps  evals  total SCF iters  wall (s)'
   call one_run('per-image (default)', .false.)
   call one_run('single shared density', .true.)

contains

   subroutine one_run(label, shared_density)
      character(len=*), intent(in) :: label
      logical, intent(in) :: shared_density
      type(c_ptr) :: h
      integer(c_int) :: st
      integer, allocatable :: residues(:)
      type(trc_opt_options_t) :: opts
      type(trc_opt_result_t) :: res
      type(timer_type) :: clk
      character(len=1024) :: line

      st = trc_create(h)
      st = trc_set_molecule(h, int(nat, c_int), real(z, c_double), reshape(real(r, c_double), [3*nat]), &
                            0_c_int, 1_c_int)
      st = trc_set_basis_json(h, c_path('basis_sets/sto-3g.json'))
      st = trc_set_method(h, c_path(''), 3_c_int)
      st = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
      st = trc_set_guess(h, TRC_GUESS_SAD, c_null_ptr, 1_c_int)

      allocate (residues(nat)); residues = 1
      opts%coords_kind = 0   !! cartesian
      opts%guess_kind = TRC_OPT_GUESS_DENSITY
      opts%fallback_guess_kind = TRC_GUESS_SAD
      opts%maxiter = 300
      opts%gtol = 8.0e-3_dp
      opts%max_step = 0.03_dp
      opts%neb_spring = 0.01_dp
      opts%neb = .true.
      opts%neb_images = 10
      ! --opt-neb-climb is deliberately NOT set here: with
      ! neb_shared_density, cb_get_gradient's scratch-slot (climbing)
      ! image is never actually evaluated before dlf_neb.f90 runs its own
      ! climbing-convergence test on the very first cycle, which then
      ! reads an un-written (zero) gradient there and reports "converged"
      ! after 9 evaluations -- a degenerate interaction specific to that
      ! configuration, not a difference this bench is trying to measure
      ! (check_opt_neb.F90's own ctest DOES use --opt-neb-climb, and does
      ! not hit this, since it never asks for the shared-density path).
      ! Comparing the SCF cost of the warm start is orthogonal to
      ! climbing, so plain NEB here keeps the comparison clean.
      opts%neb_endpoints = TRC_OPT_NEB_PERP
      opts%neb_product = r2
      opts%neb_shared_density = shared_density

      call clk%start()
      call trc_dlfind_optimize(h, nat, z, residues, r, opts, res)
      call clk%stop()

      if (res%failed) then
         print '(a,a,a)', '  ', trim(label), ': FAILED -- '//res%message
      else
         write (line, '(a20,2x,i6,2x,i5,2x,i5,2x,i14,2x,f9.2)') label, res%neb_n_images, &
            res%steps, res%energy_evaluations, res%total_scf_iterations, clk%get_elapsed_time()
         print '(a)', trim(line)
         if (allocated(res%neb_scf_iterations_per_image)) then
            print '(a,*(i6))', '    SCF iterations per image: ', res%neb_scf_iterations_per_image
         end if
      end if

      st = trc_destroy(h)
      deallocate (residues)
   end subroutine one_run

   pure function c_path(s) result(cs)
      character(len=*), intent(in) :: s
      character(kind=c_char, len=len_trim(s) + 1) :: cs
      cs = trim(s)//c_null_char
   end function c_path

end program bench_opt_neb_warmstart
