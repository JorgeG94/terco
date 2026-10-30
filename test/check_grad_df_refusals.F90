!
! The --grad --df refusals, at the C API: every out-of-scope request must
! come back as an error STATUS with a message naming the reason, never reach
! the library's own `error stop` (which would take the host process down).
! The fitted gradient is restricted, Cartesian and not range-separated
! (trc_run_gradient's own guards, trc_capi.F90).
!
! Two routes into each refusal:
!
!   RUN. trc_run_scf with --df on an unrestricted (OH), spherical (water,
!        cc-pVDZ as tagged) or range-separated (wB97X) setup. Where the SCF
!        refuses (unrestricted, range-separated), trc_run_gradient after it
!        refuses too (no SCF). The DF-SCF runs spherical bases, so there the
!        gradient's own guard must refuse, TRC_ERR_UNSUPPORTED, naming it.
!   LOAD. The same setup's SCF, run WITHOUT --df and saved (trc_save_scf),
!        then loaded into a --df context (trc_set_df, trc_load_scf):
!        trc_run_gradient is now past the SCF's refusals and its own guard
!        is what must refuse, TRC_ERR_UNSUPPORTED, naming the reason. A
!        restricted Cartesian result loaded the same way has no fit, and
!        must be refused with TRC_ERR_STATE.
!
! Plus the positive control: RHF water/6-31G* with --df runs the gradient.
!
module check_grad_df_refusals_suite
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr, c_char, c_null_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                       trc_set_method, trc_set_guess, trc_set_convergence, trc_set_verbose, &
                       trc_set_force_cartesian, trc_set_df, trc_set_jkfit_json, trc_run_scf, &
                       trc_run_gradient, trc_save_scf, trc_load_scf, trc_message, &
                       TRC_OK, TRC_ERR_STATE, TRC_ERR_UNSUPPORTED, TRC_GUESS_GWH
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_grad_df_refusals_tests, set_scratch_path

   character(len=*), parameter :: JKFIT = 'basis_sets/def2-universal-jkfit.json'

   character(len=512) :: g_scratch = ''
      !! The scratch .scf path ctest hands this program as its one
      !! command-line argument -- set once by the program below, before
      !! run_testsuite, and read by every case here (test-drive's
      !! new_unittest callback takes only `error`, so there is nowhere else
      !! to thread a run-time argument through).

contains

   subroutine set_scratch_path(path)
      character(*), intent(in) :: path
      g_scratch = path
   end subroutine set_scratch_path

   subroutine collect_grad_df_refusals_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF water/6-31G* --df runs the gradient (control)", test_control), &
                  new_unittest("UHF OH/6-31G* --df is refused (unrestricted)", test_uhf), &
                  new_unittest("spherical water/cc-pVDZ --df is refused (not Cartesian)", test_spherical), &
                  new_unittest("wB97X water/6-31G* --df is refused (range-separated)", test_rsh), &
                  new_unittest("RHF water/6-31G*, no fit, loaded into --df is refused", test_nofit) &
                  ]
   end subroutine collect_grad_df_refusals_tests

   subroutine test_control(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'RHF water/6-31G* (control)', 'water.xyz', 1, 'basis_sets/6-31gs.json', .true., '', &
                    TRC_OK, '', trim(g_scratch), .true.)
   end subroutine test_control

   subroutine test_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'UHF OH/6-31G*', 'oh.xyz', 2, 'basis_sets/6-31gs.json', .true., '', &
                    TRC_ERR_UNSUPPORTED, 'restricted', trim(g_scratch), .true.)
   end subroutine test_uhf

   subroutine test_spherical(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'spherical water/cc-pVDZ', 'water.xyz', 1, 'basis_sets/cc-pvdz.json', .false., '', &
                    TRC_ERR_UNSUPPORTED, 'Cartesian', trim(g_scratch), .true.)
   end subroutine test_spherical

   subroutine test_rsh(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'wB97X water/6-31G*', 'water.xyz', 1, 'basis_sets/6-31gs.json', .true., 'wb97x', &
                    TRC_ERR_UNSUPPORTED, 'range-separated', trim(g_scratch), .true.)
   end subroutine test_rsh

   subroutine test_nofit(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'RHF water/6-31G*, no fit', 'water.xyz', 1, 'basis_sets/6-31gs.json', .true., '', &
                    TRC_ERR_STATE, 'density fit', trim(g_scratch), .false.)
   end subroutine test_nofit

   subroutine setup(h, xyz, mult, basis, cart, func)
      type(c_ptr), intent(out) :: h
      character(*), intent(in) :: xyz, basis, func
      integer, intent(in) :: mult
      logical, intent(in) :: cart
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      call read_xyz(xyz, nat, at_z, r)
      call must(trc_create(h), 'trc_create')
      call must(trc_set_verbose(h, 0_c_int), 'trc_set_verbose')
      call must(trc_set_molecule(h, int(nat, c_int), real(at_z, c_double), r, 0_c_int, int(mult, c_int)), &
                'trc_set_molecule')
      call must(trc_set_force_cartesian(h, merge(1_c_int, 0_c_int, cart)), 'trc_set_force_cartesian')
      call must(trc_set_basis_json(h, cstr(basis)), 'trc_set_basis_json')
      call must(trc_set_method(h, cstr(func), 3_c_int), 'trc_set_method')
      call must(trc_set_convergence(h, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int), &
                'trc_set_convergence')
      call must(trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int), 'trc_set_guess')
   end subroutine setup

   subroutine run_case(error, label, xyz, mult, basis, cart, func, want, key, save_path, do_run)
      !! `want` = TRC_OK for the positive control; otherwise the status the
      !! LOAD route's gradient must return, `key` a word its message must
      !! contain. `do_run` = .false. skips the RUN route (the no-fit case,
      !! whose RUN route is the positive control).
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: label, xyz, basis, func, key, save_path
      integer, intent(in) :: mult
      logical, intent(in) :: cart, do_run
      integer(c_int), intent(in) :: want
      type(c_ptr) :: h
      integer(c_int) :: st_scf, st_grad
      character(len=512) :: msg
      logical :: ok

      ! --- RUN: --df from the start ---
      if (do_run) then
         call setup(h, xyz, mult, basis, cart, func)
         call must(trc_set_jkfit_json(h, cstr(JKFIT)), 'trc_set_jkfit_json')
         call must(trc_set_df(h, 1_c_int), 'trc_set_df')
         st_scf = trc_run_scf(h)
         st_grad = trc_run_gradient(h)
         msg = message(h)
         if (want == TRC_OK) then
            ok = st_scf == TRC_OK .and. st_grad == TRC_OK
         else if (st_scf == TRC_OK) then
            ! The DF-SCF supports this setup (spherical --df does), so the
            ! gradient's own guard is what must refuse it.
            ok = st_grad == want .and. index(msg, key) > 0
         else
            ok = st_grad /= TRC_OK
         end if
         print '(a,a,a,i0,a,i0,a,a)', merge('ok   ', 'FAIL ', ok), label, ' RUN : scf ', st_scf, ' grad ', &
            st_grad, '  -- ', trim(msg)
         call check(error, ok, label//' RUN: scf status '//to_string(int(st_scf))//', grad status '// &
                    to_string(int(st_grad))//', message "'//trim(msg)//'"')
         call must(trc_destroy(h), 'trc_destroy')
         if (allocated(error)) return
      end if
      if (want == TRC_OK) return

      ! --- LOAD: the SCF without --df, saved, loaded into a --df context ---
      call setup(h, xyz, mult, basis, cart, func)
      call must(trc_run_scf(h), 'trc_run_scf (no --df)')
      call must(trc_save_scf(h, cstr(save_path)), 'trc_save_scf')
      call must(trc_set_jkfit_json(h, cstr(JKFIT)), 'trc_set_jkfit_json')
      call must(trc_set_df(h, 1_c_int), 'trc_set_df')
      call must(trc_load_scf(h, cstr(save_path)), 'trc_load_scf')
      st_grad = trc_run_gradient(h)
      msg = message(h)
      ok = st_grad == want .and. index(msg, key) > 0
      print '(a,a,a,i0,a,a)', merge('ok   ', 'FAIL ', ok), label, ' LOAD: grad ', st_grad, '  -- ', trim(msg)
      call check(error, ok, label//' LOAD: expected grad status '//to_string(int(want))// &
                 ' containing "'//key//'", got status '//to_string(int(st_grad))// &
                 ', message "'//trim(msg)//'"')
      call must(trc_destroy(h), 'trc_destroy')
   end subroutine run_case

   function message(h) result(m)
      type(c_ptr), intent(in) :: h
      character(len=512) :: m
      character(kind=c_char) :: buf(512)
      integer :: i
      call must(trc_message(h, buf, 512_c_int), 'trc_message')
      m = ''
      do i = 1, 511
         if (buf(i) == c_null_char) exit
         m(i:i) = buf(i)
      end do
   end function message

   subroutine must(st, what)
      integer(c_int), intent(in) :: st
      character(*), intent(in) :: what
      if (st /= TRC_OK) then
         print '(a,a,a,i0)', 'check_grad_df_refusals: ', what, ' failed, status ', st
         stop 1
      end if
   end subroutine must

   pure function cstr(s) result(cs)
      character(*), intent(in) :: s
      character(len=len(s) + 1) :: cs
      cs = s//achar(0)
   end function cstr

end module check_grad_df_refusals_suite

program check_grad_df_refusals
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_grad_df_refusals_suite, only: collect_grad_df_refusals_tests, set_scratch_path
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=512) :: scratch
   integer :: stat

   if (command_argument_count() < 1) then
      print '(a)', 'usage: check_grad_df_refusals <scratch .scf path>'; stop 1
   end if
   call get_command_argument(1, scratch)
   call set_scratch_path(trim(scratch))

   stat = 0
   testsuites = [new_testsuite("check_grad_df_refusals", collect_grad_df_refusals_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_grad_df_refusals
