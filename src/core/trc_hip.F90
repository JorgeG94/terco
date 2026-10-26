!
! hipBLAS, hipSOLVER and the HIP runtime, as C interfaces.
!
! Only what trc_linalg calls. ROCm's own Fortran route is hipfort, which is
! slow to build and one more thing to install; hipBLAS ships a Fortran module
! but hipSOLVER does not. The C ABI is stable, so these declarations are the
! whole dependency.
!
! Everything is a plain `bind(c)` declaration: handles are `type(c_ptr)`, enums
! are `integer(c_int)` with the values from hipblas-common.h and
! hipsolver-types.h, and scalars alpha/beta/result are passed by reference
! (the handles' default pointer mode is host). Arrays are `real(c_double)`
! assumed-size so the caller can hand a device address obtained from
! `!$omp target data use_device_addr(...)` straight through.
!
! Only compiled into a TERCO_HIP build; otherwise the module is empty.
!
module trc_hip
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double
   implicit none
   private

#ifdef TERCO_HIP
   public :: hipblasCreate, hipblasDestroy
   public :: hipsolverCreate, hipsolverDestroy
   public :: hipblasDgemm, hipblasDgemmStridedBatched, hipblasDgeam
   public :: hipblasDgemv, hipblasDtrsm, hipblasDdot
   public :: hipsolverDsyevd_bufferSize, hipsolverDsyevd
   public :: hipsolverDpotrf_bufferSize, hipsolverDpotrf
   public :: hipsolverDgetrf_bufferSize, hipsolverDgetrf
   public :: hipsolverDgetrs_bufferSize, hipsolverDgetrs
   public :: hipDeviceSynchronize

   ! hipblasOperation_t / hipsolverOperation_t (the same enum).
   integer(c_int), parameter, public :: HIP_OP_N = 111
   integer(c_int), parameter, public :: HIP_OP_T = 112
   ! hipblasFillMode_t / hipsolverFillMode_t.
   integer(c_int), parameter, public :: HIP_FILL_UPPER = 121
   integer(c_int), parameter, public :: HIP_FILL_LOWER = 122
   ! hipblasDiagType_t.
   integer(c_int), parameter, public :: HIP_DIAG_NON_UNIT = 131
   ! hipblasSideMode_t.
   integer(c_int), parameter, public :: HIP_SIDE_LEFT = 141
   integer(c_int), parameter, public :: HIP_SIDE_RIGHT = 142
   ! hipsolverEigMode_t.
   integer(c_int), parameter, public :: HIP_EIG_VECTOR = 202

   interface
      function hipblasCreate(handle) result(st) bind(c, name="hipblasCreate")
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), intent(out) :: handle
         integer(c_int) :: st
      end function hipblasCreate

      function hipblasDestroy(handle) result(st) bind(c, name="hipblasDestroy")
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int) :: st
      end function hipblasDestroy

      function hipsolverCreate(handle) result(st) bind(c, name="hipsolverCreate")
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), intent(out) :: handle
         integer(c_int) :: st
      end function hipsolverCreate

      function hipsolverDestroy(handle) result(st) bind(c, name="hipsolverDestroy")
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int) :: st
      end function hipsolverDestroy

      function hipblasDgemm(handle, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) &
         result(st) bind(c, name="hipblasDgemm")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: ta, tb, m, n, k, lda, ldb, ldc
         real(c_double), intent(in) :: alpha, beta
         real(c_double), intent(in) :: a(*), b(*)
         real(c_double), intent(inout) :: c(*)
         integer(c_int) :: st
      end function hipblasDgemm

      function hipblasDgemmStridedBatched(handle, ta, tb, m, n, k, alpha, a, lda, sa, b, ldb, sb, &
                                          beta, c, ldc, sc, count) &
         result(st) bind(c, name="hipblasDgemmStridedBatched")
         import :: c_ptr, c_int, c_int64_t, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: ta, tb, m, n, k, lda, ldb, ldc, count
         integer(c_int64_t), value, intent(in) :: sa, sb, sc
         real(c_double), intent(in) :: alpha, beta
         real(c_double), intent(in) :: a(*), b(*)
         real(c_double), intent(inout) :: c(*)
         integer(c_int) :: st
      end function hipblasDgemmStridedBatched

      function hipblasDgeam(handle, ta, tb, m, n, alpha, a, lda, beta, b, ldb, c, ldc) &
         result(st) bind(c, name="hipblasDgeam")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: ta, tb, m, n, lda, ldb, ldc
         real(c_double), intent(in) :: alpha, beta
         real(c_double), intent(in) :: a(*), b(*)
         real(c_double), intent(inout) :: c(*)
         integer(c_int) :: st
      end function hipblasDgeam

      function hipblasDgemv(handle, trans, m, n, alpha, a, lda, x, incx, beta, y, incy) &
         result(st) bind(c, name="hipblasDgemv")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: trans, m, n, lda, incx, incy
         real(c_double), intent(in) :: alpha, beta
         real(c_double), intent(in) :: a(*), x(*)
         real(c_double), intent(inout) :: y(*)
         integer(c_int) :: st
      end function hipblasDgemv

      function hipblasDtrsm(handle, side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb) &
         result(st) bind(c, name="hipblasDtrsm")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: side, uplo, trans, diag, m, n, lda, ldb
         real(c_double), intent(in) :: alpha
         real(c_double), intent(in) :: a(*)
         real(c_double), intent(inout) :: b(*)
         integer(c_int) :: st
      end function hipblasDtrsm

      function hipblasDdot(handle, n, x, incx, y, incy, r) result(st) bind(c, name="hipblasDdot")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: n, incx, incy
         real(c_double), intent(in) :: x(*), y(*)
         real(c_double), intent(out) :: r
         integer(c_int) :: st
      end function hipblasDdot

      function hipsolverDsyevd_bufferSize(handle, jobz, uplo, n, a, lda, w, lwork) &
         result(st) bind(c, name="hipsolverDsyevd_bufferSize")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: jobz, uplo, n, lda
         real(c_double), intent(in) :: a(*), w(*)
         integer(c_int), intent(out) :: lwork
         integer(c_int) :: st
      end function hipsolverDsyevd_bufferSize

      function hipsolverDsyevd(handle, jobz, uplo, n, a, lda, w, work, lwork, info) &
         result(st) bind(c, name="hipsolverDsyevd")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: jobz, uplo, n, lda, lwork
         real(c_double), intent(inout) :: a(*), w(*), work(*)
         integer(c_int), intent(inout) :: info
         integer(c_int) :: st
      end function hipsolverDsyevd

      function hipsolverDpotrf_bufferSize(handle, uplo, n, a, lda, lwork) &
         result(st) bind(c, name="hipsolverDpotrf_bufferSize")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: uplo, n, lda
         real(c_double), intent(in) :: a(*)
         integer(c_int), intent(out) :: lwork
         integer(c_int) :: st
      end function hipsolverDpotrf_bufferSize

      function hipsolverDpotrf(handle, uplo, n, a, lda, work, lwork, info) &
         result(st) bind(c, name="hipsolverDpotrf")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: uplo, n, lda, lwork
         real(c_double), intent(inout) :: a(*), work(*)
         integer(c_int), intent(inout) :: info
         integer(c_int) :: st
      end function hipsolverDpotrf

      function hipsolverDgetrf_bufferSize(handle, m, n, a, lda, lwork) &
         result(st) bind(c, name="hipsolverDgetrf_bufferSize")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: m, n, lda
         real(c_double), intent(in) :: a(*)
         integer(c_int), intent(out) :: lwork
         integer(c_int) :: st
      end function hipsolverDgetrf_bufferSize

      function hipsolverDgetrf(handle, m, n, a, lda, work, lwork, ipiv, info) &
         result(st) bind(c, name="hipsolverDgetrf")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: m, n, lda, lwork
         real(c_double), intent(inout) :: a(*), work(*)
         integer(c_int), intent(inout) :: ipiv(*)
         integer(c_int), intent(inout) :: info
         integer(c_int) :: st
      end function hipsolverDgetrf

      function hipsolverDgetrs_bufferSize(handle, trans, n, nrhs, a, lda, ipiv, b, ldb, lwork) &
         result(st) bind(c, name="hipsolverDgetrs_bufferSize")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: trans, n, nrhs, lda, ldb
         real(c_double), intent(in) :: a(*), b(*)
         integer(c_int), intent(in) :: ipiv(*)
         integer(c_int), intent(out) :: lwork
         integer(c_int) :: st
      end function hipsolverDgetrs_bufferSize

      function hipsolverDgetrs(handle, trans, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info) &
         result(st) bind(c, name="hipsolverDgetrs")
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value, intent(in) :: handle
         integer(c_int), value, intent(in) :: trans, n, nrhs, lda, ldb, lwork
         real(c_double), intent(in) :: a(*)
         real(c_double), intent(inout) :: b(*), work(*)
         integer(c_int), intent(in) :: ipiv(*)
         integer(c_int), intent(inout) :: info
         integer(c_int) :: st
      end function hipsolverDgetrs

      function hipDeviceSynchronize() result(st) bind(c, name="hipDeviceSynchronize")
         import :: c_int
         implicit none
         integer(c_int) :: st
      end function hipDeviceSynchronize
   end interface
#endif

end module trc_hip
