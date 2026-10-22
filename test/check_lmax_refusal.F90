module check_lmax_refusal_suite
   !! A basis with a shell above this build's four-centre ceiling (an f
   !! basis on the default TERCO_LMAX=2 build) is refused up front, with a
   !! message that names the shell and the rebuild that would take it. It
   !! used to have no kernel at all: its integrals came out as NaN and the
   !! run died much later in the eigensolver ("pic_syev failed"), with
   !! nothing pointing at the basis. Every case builds its basis one shell
   !! above `LMAX`, so it means the same thing on an LMAX=3 build.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_tables, only: LMAX
   use trc_api, only: trc_basis_t
   use trc_eri, only: eri_lmax_refusal
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_arrays, &
                               trc_run_scf, trc_message, TRC_OK, TRC_ERR_UNSUPPORTED
   implicit none
   private
   public :: collect_lmax_refusal

contains

   subroutine collect_lmax_refusal(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("a basis within the ceiling passes", test_within), &
                  new_unittest("one shell above it is named", test_above), &
                  new_unittest("trc_run_scf refuses it with a status", test_capi) &
                  ]
   end subroutine collect_lmax_refusal

   subroutine two_shell_basis(b, ltop)
      !! An s shell on each of two hydrogens and one shell of angular
      !! momentum `ltop` on the first.
      type(trc_basis_t), intent(out) :: b
      integer, intent(in) :: ltop
      real(dp) :: sh_e(1, 3), sh_c(1, 3), sh_r(3, 3), at_r(3, 2)
      at_r = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.4_dp], [3, 2])
      sh_e = reshape([1.0_dp, 1.0_dp, 0.8_dp], [1, 3])
      sh_c = 1.0_dp
      sh_r(:, 1) = at_r(:, 1)
      sh_r(:, 2) = at_r(:, 2)
      sh_r(:, 3) = at_r(:, 1)
      call b%build(3, [0, 0, ltop], [1, 1, 1], sh_e, sh_c, sh_r, 2, [1.0_dp, 1.0_dp], at_r, 1)
   end subroutine two_shell_basis

   subroutine test_within(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      call two_shell_basis(b, LMAX)
      call check(error, eri_lmax_refusal(b), "", "a basis at the ceiling must not be refused")
      call b%release()
   end subroutine test_within

   subroutine test_above(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      character(len=:), allocatable :: msg
      character(len=8) :: want
      call two_shell_basis(b, LMAX + 1)
      msg = eri_lmax_refusal(b)
      call b%release()
      write (want, "(a,i0,a)") "(l = ", LMAX + 1, ")"
      call check(error, index(msg, trim(want)) > 0, "the refusal names the shell's l: '"//msg//"'")
   end subroutine test_above

   subroutine test_capi(error)
      !! The C API path a host (and the terco CLI) takes: a status and a
      !! message, before any guess or integral runs.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, rc_run
      character(kind=c_char) :: buf(512)
      character(len=512) :: msg
      integer :: i
      rc = trc_create(h)
      rc = trc_set_molecule(h, 2_c_int, [1.0_c_double, 1.0_c_double], &
                            [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 1.4_c_double], &
                            0_c_int, 1_c_int)
      rc = trc_set_basis_arrays(h, 3_c_int, 1_c_int, [0_c_int, 0_c_int, int(LMAX + 1, c_int)], &
                                [1_c_int, 1_c_int, 1_c_int], [1.0_c_double, 1.0_c_double, 0.8_c_double], &
                                [1.0_c_double, 1.0_c_double, 1.0_c_double], &
                                [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 1.4_c_double, &
                                 0.0_c_double, 0.0_c_double, 0.0_c_double])
      call check(error, rc == TRC_OK, "trc_set_basis_arrays accepts the basis; only the run is refused")
      if (allocated(error)) return
      rc_run = trc_run_scf(h)
      buf = c_null_char
      rc = trc_message(h, buf, 512_c_int)
      msg = ""
      do i = 1, size(buf)
         if (buf(i) == c_null_char) exit
         msg(i:i) = buf(i)
      end do
      rc = trc_destroy(h)
      call check(error, rc_run == TRC_ERR_UNSUPPORTED, "trc_run_scf returns TRC_ERR_UNSUPPORTED")
      if (allocated(error)) return
      call check(error, index(msg, "four-centre integrals stop at") > 0, &
                 "trc_run_scf's message says why: '"//trim(msg)//"'")
   end subroutine test_capi

end module check_lmax_refusal_suite

program check_lmax_refusal
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_lmax_refusal_suite, only: collect_lmax_refusal
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_lmax_refusal, error_unit, stat)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_lmax_refusal
