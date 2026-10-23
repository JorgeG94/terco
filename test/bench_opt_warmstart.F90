!
! SCF iterations per --opt step (summed over every accepted step,
! trc_opt_result_t%total_scf_iterations), for the mqc baseline
! (TRC_OPT_GUESS_NONE: re-run SAD at every step, DLFIND_INTEGRATION_
! NOTE.md section 5) against the two warm starts (density, projected),
! water/STO-3G and a 6-water cluster/def2-SVP, HF and B3LYP -- this PR's
! own warm-start measurement. Informational only, like bench_warmstart.F90
! (which measures the SAME two warm starts for trc_run_hessian's
! displacements, at a much smaller step): prints a table, asserts
! nothing, and is not registered as a ctest. Run by hand:
!
!   build-cpu/test/bench_opt_warmstart <6-water cluster xyz>
!
! (the cluster is the first 18 atoms of .claude/notes/bench/fig/inp/
! w20.xyz in the maintainer's own notes checkout, not shipped here --
! pass any 6-water-sized xyz, or omit to measure water alone).
!
program bench_opt_warmstart
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_null_ptr, c_null_char
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_method, &
                       trc_set_convergence, trc_set_guess, TRC_GUESS_SAD
   use trc_api, only: trc_bind_device
   use trc_dlfind_bridge, only: trc_opt_options_t, trc_opt_result_t, trc_dlfind_optimize, &
                                TRC_OPT_GUESS_DENSITY, TRC_OPT_GUESS_PROJECTED, TRC_OPT_GUESS_NONE
   use trc_test_basis, only: read_xyz
   implicit none

   character(len=1024) :: arg
   character(len=:), allocatable :: w6_path
   integer :: dev

   dev = trc_bind_device(0)

   print '(a)', 'system           func    guess                 steps  evals  total SCF iters  wall (s)'
   call run_case('water/sto-3g', 'water.xyz', 'basis_sets/sto-3g.json', 'hf')
   call run_case('water/sto-3g', 'water.xyz', 'basis_sets/sto-3g.json', 'b3lyp')

   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      w6_path = trim(arg)
      call run_case('(H2O)6/def2-svp', w6_path, 'basis_sets/def2-svp.json', 'hf')
      call run_case('(H2O)6/def2-svp', w6_path, 'basis_sets/def2-svp.json', 'b3lyp')
   else
      print '(a)', '(pass a 6-water cluster xyz as argv(1) to also measure it)'
   end if

contains

   subroutine run_case(label, xyzfile, basis_json, func)
      character(len=*), intent(in) :: label, xyzfile, basis_json, func
      integer :: nat
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)

      call read_xyz(xyzfile, nat, z, r)
      call one_guess(label, basis_json, func, nat, z, r, TRC_OPT_GUESS_NONE, 'none (mqc baseline)')
      call one_guess(label, basis_json, func, nat, z, r, TRC_OPT_GUESS_DENSITY, 'density')
      call one_guess(label, basis_json, func, nat, z, r, TRC_OPT_GUESS_PROJECTED, 'projected')
      deallocate (z, r)
   end subroutine run_case

   subroutine one_guess(label, basis_json, func, nat, z, r, guess_kind, guess_name)
      character(len=*), intent(in) :: label, basis_json, func, guess_name
      integer, intent(in) :: nat, z(nat), guess_kind
      real(dp), intent(in) :: r(3, nat)
      type(c_ptr) :: h
      integer(c_int) :: st
      integer, allocatable :: residues(:)
      type(trc_opt_options_t) :: opts
      type(trc_opt_result_t) :: res
      type(timer_type) :: clk
      character(len=1024) :: line
      character(len=:), allocatable :: func_arg

      func_arg = func
      if (func_arg == 'hf') func_arg = ''

      st = trc_create(h)
      st = trc_set_molecule(h, int(nat, c_int), real(z, c_double), reshape(real(r, c_double), [3*nat]), &
                            0_c_int, 1_c_int)
      st = trc_set_basis_json(h, c_path(basis_json))
      st = trc_set_method(h, c_path(func_arg), 3_c_int)
      st = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
      st = trc_set_guess(h, TRC_GUESS_SAD, c_null_ptr, 1_c_int)

      allocate (residues(nat)); residues = 1
      opts%coords_kind = 0   !! cartesian
      opts%guess_kind = guess_kind
      opts%fallback_guess_kind = TRC_GUESS_SAD
      opts%maxiter = 60
      opts%gtol = 4.5e-4_dp

      call clk%start()
      call trc_dlfind_optimize(h, nat, z, residues, r, opts, res)
      call clk%stop()

      if (res%failed) then
         print '(a,a,a,a,a,a)', '  ', trim(label), '/', trim(func), '/', trim(guess_name)// &
            ': FAILED -- '//res%message
      else
         write (line, '(a16,2x,a6,2x,a20,i5,2x,i5,2x,i14,2x,f9.2)') label, func, guess_name, &
            res%steps, res%energy_evaluations, res%total_scf_iterations, clk%get_elapsed_time()
         print '(a)', trim(line)
      end if

      st = trc_destroy(h)
      deallocate (residues)
   end subroutine one_guess

   pure function c_path(s) result(cs)
      character(len=*), intent(in) :: s
      character(kind=c_char, len=len_trim(s) + 1) :: cs
      cs = trim(s)//c_null_char
   end function c_path

end program bench_opt_warmstart
