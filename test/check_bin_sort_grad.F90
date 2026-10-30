!
! The four-centre gradient after a TRC_BIN_SORT resort.
!
! TRC_BIN_SORT=1 reorders the pairs inside each bin the first time a Fock
! build runs (fock_bins), and the bin-interleaved quartet tables are keyed
! on that order, so they are rebuilt then. The gradient reads those tables
! afterwards on the same engine -- including the exponents
! (`quartets%a`) only the Rys gradient kernels use, which a resort once
! left in the old order. So: one Fock build and then the gradient on an
! engine with the resort, against the same on an engine without it. The
! two differ only by summation order.
!
module check_bin_sort_grad_suite
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_basis_types, only: trc_basis_t
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_bin_sort_grad

   real(dp), parameter :: TOL = 1.0e-10_dp
      !! The resort only reorders sums; measured |dG| is at round-off.

   interface
      function c_setenv(name, value, overwrite) bind(c, name="setenv") result(r)
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: r
      end function c_setenv
   end interface

contains

   subroutine collect_bin_sort_grad(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("gradient after a TRC_BIN_SORT resort matches the unsorted one", &
                                test_resort_gradient)]
   end subroutine collect_bin_sort_grad

   subroutine test_resort_gradient(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: g_plain(:, :), g_sorted(:, :)
      real(dp) :: worst

      call gradient_after_fock("0", g_plain, error)
      if (allocated(error)) return
      call gradient_after_fock("1", g_sorted, error)
      if (allocated(error)) return
      worst = maxval(abs(g_sorted - g_plain))
      call check(error, worst < TOL, "gradient after the resort differs from the unsorted one: max |dG| = "// &
                 trim(real_text(worst)))
   end subroutine test_resort_gradient

   subroutine gradient_after_fock(sort, grad, error)
      ! A fresh engine with TRC_BIN_SORT=`sort`: one Fock build, then the
      ! two-electron gradient of the same density.
      character(*), intent(in) :: sort
      real(dp), allocatable, intent(out) :: grad(:, :)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      real(dp), allocatable :: dmat(:, :), gmat(:, :), r(:, :)
      integer, allocatable :: z(:)
      integer :: nat, nao, i, j
      integer(c_int) :: rc

      rc = c_setenv("TRC_BIN_SORT"//c_null_char, sort//c_null_char, 1_c_int)
      call check(error, rc == 0, "setenv(TRC_BIN_SORT) failed")
      if (allocated(error)) return
      call read_xyz("water.xyz", nat, z, r)
      call trc_basis_from_json("basis_sets/cc-pvdz.json", nat, z, r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao
      call bas%to_device()
      call eri%build(bas, 1.0e-13_dp)
      allocate (dmat(nao, nao), gmat(nao, nao), grad(3, nat))
      ! Deterministic, dense and symmetric, as check_rys_vs_hgp's.
      do j = 1, nao
         do i = 1, nao
            dmat(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j, dp))
         end do
      end do
      dmat = 0.5_dp*(dmat + transpose(dmat))
      call eri%fock(bas, dmat, gmat)
      grad = 0.0_dp
      call eri%gradient(bas, reshape(dmat, [nao, nao, 1]), grad)
      call eri%release()
      call bas%release()
      rc = c_setenv("TRC_BIN_SORT"//c_null_char, "0"//c_null_char, 1_c_int)
   end subroutine gradient_after_fock

   function real_text(x) result(s)
      real(dp), intent(in) :: x
      character(len=24) :: s
      write (s, "(es12.4)") x
   end function real_text

end module check_bin_sort_grad_suite

program check_bin_sort_grad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_bin_sort_grad_suite, only: collect_bin_sort_grad
   implicit none
   integer :: stat
   stat = 0
   ! Serial: the environment variable it sets is process-wide.
   call run_testsuite(collect_bin_sort_grad, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_bin_sort_grad
