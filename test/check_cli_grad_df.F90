!
! The CLI's --grad --df: what `terco` prints to its JSON summary for a
! density-fitted run must be exactly what the library itself
! (trc_run_gradient/trc_gradient over the context's own fit, the same C API
! app/main.F90 calls) computes for the identical molecule, basis, jkfit
! basis, functional and guess -- check_cli_grad.F90's wiring check, for the
! fitted gradient (trc_scf_grad_df.F90; check_df_grad_fd.F90 and
! check_df_grad_vs_pyscf.F90 are its physics checks). B3LYP, so the fitted
! J, the fitted K (scaled by 0.2) and the XC term all run, on water/6-31G*
! with def2-universal-jkfit, Cartesian as --df requires. Same tolerance and
! reasoning as check_cli_grad.F90's TOL.
!
program check_cli_grad_df
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                       trc_set_method, trc_set_guess, trc_set_convergence, trc_set_screening, &
                       trc_set_incremental, trc_run_scf, trc_run_gradient, trc_gradient, &
                       trc_set_verbose, trc_set_force_cartesian, trc_set_nlc_grid, &
                       trc_set_df, trc_set_jkfit_json, &
                       TRC_OK, TRC_GUESS_GWH
   use trc_test_basis, only: read_xyz
   implicit none

   character(len=:), allocatable :: cli_exe, probe_json, probe_out
   character(len=1024) :: arg
   character(len=:), allocatable :: cmd
   integer :: nat, i, rc
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: r(:, :)
   real(dp), allocatable :: grad_api(:, :), grad_cli(:, :)
   type(c_ptr) :: h
   integer(c_int) :: st
   logical :: found
   type(json_file) :: json
   real(dp), allocatable :: row(:)
   real(dp) :: worst
   real(dp), parameter :: TOL = 1.0e-6_dp
   !! Same value and reasoning as check_cli_grad.F90's own TOL: two
   !! separate OS processes, threaded BLAS that is not associative, so
   !! their GWH guesses and DIIS iterations (now also the XC quadrature's
   !! own batched GEMMs) can land at very slightly different densities.

   call get_command_argument(1, arg)
   cli_exe = trim(arg)
   call get_command_argument(2, arg)
   probe_json = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0) then
      print '(a)', 'check_cli_grad_df: usage: check_cli_grad_df <terco CLI binary> <scratch .json path>'
      stop 1
   end if
   probe_out = probe_json//'.out'

   call read_xyz('water.xyz', nat, at_z, r)

   ! --- "the library's" own gradient, straight through the C API -------
   st = trc_create(h)
   if (st /= TRC_OK) then
      print '(a)', 'check_cli_grad_df: trc_create failed'; stop 1
   end if
   st = trc_set_verbose(h, 1_c_int)
   call must_ok(st, 'trc_set_verbose')
   st = trc_set_molecule(h, int(nat, c_int), real(at_z, c_double), r, 0_c_int, 1_c_int)
   call must_ok(st, 'trc_set_molecule')
   st = trc_set_force_cartesian(h, 1_c_int)
   call must_ok(st, 'trc_set_force_cartesian')
   st = trc_set_basis_json(h, c_string_local('basis_sets/6-31gs.json'))
   call must_ok(st, 'trc_set_basis_json')
   st = trc_set_jkfit_json(h, c_string_local('basis_sets/def2-universal-jkfit.json'))
   call must_ok(st, 'trc_set_jkfit_json')
   st = trc_set_df(h, 1_c_int)
   call must_ok(st, 'trc_set_df')
   st = trc_set_method(h, c_string_local('b3lyp'), 3_c_int)
   call must_ok(st, 'trc_set_method')
   st = trc_set_nlc_grid(h, 1_c_int)
   call must_ok(st, 'trc_set_nlc_grid')
   st = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-5_c_double, 100_c_int, 10_c_int)
   call must_ok(st, 'trc_set_convergence')
   st = trc_set_screening(h, 1.0e-12_c_double)
   call must_ok(st, 'trc_set_screening')
   st = trc_set_incremental(h, 0_c_int)
   call must_ok(st, 'trc_set_incremental')
   st = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
   call must_ok(st, 'trc_set_guess')
   st = trc_run_scf(h)
   call must_ok(st, 'trc_run_scf')
   st = trc_run_gradient(h)
   call must_ok(st, 'trc_run_gradient')
   allocate (grad_api(3, nat))
   st = trc_gradient(h, grad_api)
   call must_ok(st, 'trc_gradient')
   st = trc_destroy(h)

   ! --- the CLI's own printed/JSON gradient, same molecule/basis/func ---
   cmd = trim(cli_exe)//' water.xyz -b basis_sets/6-31gs.json --df --jkfit basis_sets/def2-universal-jkfit.json '// &
         '--func b3lyp --guess gwh --grad --json '// &
         probe_json//' > '//probe_out//' 2>&1'
   call execute_command_line(cmd, exitstat=rc, wait=.true.)
   if (rc /= 0) then
      print '(a,i0)', 'check_cli_grad_df: the CLI run itself failed, exit code ', rc; stop 1
   end if

   allocate (grad_cli(3, nat))
   call json%initialize()
   call json%load(filename=probe_json)
   if (json%failed()) then
      print '(a)', 'check_cli_grad_df: cannot read '//probe_json; stop 1
   end if
   do i = 1, nat
      call json%get('gradient('//to_string(i)//')', row, found)
      if (.not. found .or. size(row) /= 3) then
         print '(a,i0,a)', 'check_cli_grad_df: no gradient row ', i, ' in the CLI''s JSON output'; stop 1
      end if
      grad_cli(:, i) = row
   end do
   call json%destroy()

   worst = maxval(abs(grad_api - grad_cli))
   print '(a,es10.2)', '  check_cli_grad_df: worst|CLI JSON - direct trc_gradient| ', worst
   if (worst > TOL) then
      print '(a)', 'check_cli_grad_df: RESULT: FAIL'; stop 1
   end if
   print '(a)', 'check_cli_grad_df: RESULT: PASS'

contains

   subroutine must_ok(status, what)
      integer(c_int), intent(in) :: status
      character(*), intent(in) :: what
      if (status /= TRC_OK) then
         print '(a,a,a,i0)', 'check_cli_grad_df: ', what, ' failed, status ', status; stop 1
      end if
   end subroutine must_ok

   pure function c_string_local(s) result(cs)
      !! `s` NUL-terminated, for a `character(kind=c_char) :: x(*)` dummy --
      !! the same convention app/trc_app.F90's own `c_string` uses.
      character(*), intent(in) :: s
      character(len=len(s) + 1) :: cs
      cs = s//achar(0)
   end function c_string_local

end program check_cli_grad_df
