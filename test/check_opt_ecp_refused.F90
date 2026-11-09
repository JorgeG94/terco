!
! A geometry optimisation with an effective core potential attached is
! refused before its first step. The analytic gradient has no ECP term yet,
! so trc_run_gradient refuses it, and a refused gradient inside DL-FIND is
! `error stop` -- the end of a C caller's process. HI in def2-SVP (the
! basis's own ECP) through terco_opt.h's trc_opt_run must RETURN, with
! TRC_ERR_UNSUPPORTED and the gradient's message, without running an SCF.
!
! The CLI's own up-front refusal is in check_opt_refusals.
!
! A test-drive suite.
!
module check_opt_ecp_refused_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char, c_null_ptr
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_constants, only: ANGSTROM_TO_BOHR
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_iterations, &
                               TRC_OK, TRC_ERR_UNSUPPORTED
   use trc_opt_c_interfaces, only: trc_opt_result_create, trc_opt_result_destroy, trc_opt_run, &
                                   trc_opt_result_message
   implicit none
   private
   public :: collect_opt_ecp_refused

contains

   subroutine collect_opt_ecp_refused(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("trc_opt_run refuses an ECP molecule before its first step", test_refused)]
   end subroutine collect_opt_ecp_refused

   subroutine test_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h, res
      integer(c_int) :: rc, st, niter
      real(c_double) :: xyz(6)
      character(kind=c_char) :: msg(1024)
      character(len=:), allocatable :: text
      integer :: i

      ! check_ecp_scf's HI: H at the origin, I 1.609 Angstrom along z.
      xyz = [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, &
             real(1.609_dp*ANGSTROM_TO_BOHR, c_double)]
      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, 2_c_int, [1.0_c_double, 53.0_c_double], xyz, 0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, "basis_sets/def2-svp.json"//c_null_char)
      if (rc == TRC_OK) rc = trc_opt_result_create(res)
      call check(error, rc == TRC_OK, "HI/def2-SVP context and result setup failed")
      if (allocated(error)) return
      rc = trc_opt_run(h, 2_c_int, [1_c_int, 53_c_int], xyz, c_null_ptr, res)
      msg = c_null_char
      st = trc_opt_result_message(res, msg, int(size(msg), c_int))
      text = ""
      do i = 1, size(msg)
         if (msg(i) == c_null_char) exit
         text = text//msg(i)
      end do
      ! No SCF ran: trc_iterations refuses before the first one.
      st = trc_iterations(h, niter)
      call check(error, st /= TRC_OK, "an SCF ran before the refusal")
      if (.not. allocated(error)) then
         call check(error, rc == TRC_ERR_UNSUPPORTED .and. index(text, "effective-core-potential") > 0, &
                    "trc_opt_run should refuse with TRC_ERR_UNSUPPORTED and the gradient's ECP message", &
                    more=text)
      end if
      st = trc_opt_result_destroy(res)
      st = trc_destroy(h)
   end subroutine test_refused

end module check_opt_ecp_refused_suite

program check_opt_ecp_refused
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_ecp_refused_suite, only: collect_opt_ecp_refused
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_opt_ecp_refused, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_ecp_refused
