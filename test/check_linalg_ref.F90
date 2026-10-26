module check_linalg_ref_suite
   !! trc_linalg's raw entry points against plain Fortran: matmul, transpose
   !! and residual identities (L L^T = A, A V = V W, A X = B). check_linalg_mat
   !! only compares one form of trc_linalg with another, so a wrong backend
   !! (cuBLAS, hipBLAS or the host BLAS) would pass it; this suite checks the
   !! numbers. The matrices are closed-form fills, no RNG, so the host and
   !! device builds see the same data. Device builds map every array first:
   !! trc_linalg hands the library a device address, and a host array there
   !! is a segfault.
   use pic_types, only: dp, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_linalg, only: trc_linalg_t
   implicit none
   private
   public :: collect_linalg_ref

   real(dp), parameter :: TOL = 1.0e-10_dp

contains

   subroutine collect_linalg_ref(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("gemm, all transposes", test_gemm), &
                  new_unittest("gemm strided batched", test_gemm_strided), &
                  new_unittest("transpose and S + S^T", test_transpose), &
                  new_unittest("gemv", test_gemv), &
                  new_unittest("dot", test_dot), &
                  new_unittest("trsm", test_trsm), &
                  new_unittest("potrf: L L^T = A", test_potrf), &
                  new_unittest("getrf + getrs: A X = B", test_getrs), &
                  new_unittest("syev: A V = V W", test_syev) &
                  ]
   end subroutine collect_linalg_ref

   pure function fill(off, i, j) result(v)
      integer, intent(in) :: off, i, j
      real(dp) :: v
      v = sin(real(3*i - 2*j + off, dp)) + 0.07_dp*real(i, dp) - 0.03_dp*real(j, dp)
   end function fill

   function matrix(off, m, n) result(a)
      integer, intent(in) :: off, m, n
      real(dp) :: a(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            a(i, j) = fill(off, i, j)
         end do
      end do
   end function matrix

   function vec(off, n) result(v)
      integer, intent(in) :: off, n
      real(dp) :: v(n)
      integer :: i
      do i = 1, n
         v(i) = fill(off, i, 1)
      end do
   end function vec

   !! Symmetric positive definite: M M^T + n I.
   function spd(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n, n), m(n, n)
      integer :: i
      m = matrix(7, n, n)
      a = matmul(m, transpose(m))
      do i = 1, n
         a(i, i) = a(i, i) + real(n, dp)
      end do
   end function spd

   subroutine test_gemm(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: MM = 7, NN = 5, KK = 6
      character, parameter :: OPS(2) = ["N", "T"]
      real(dp), allocatable :: a(:, :), b(:, :), c(:, :), ref(:, :), op_a(:, :), op_b(:, :)
      real(dp) :: alpha, beta
      type(trc_linalg_t) :: lp
      integer :: ia, ib
      logical :: ta, tb
      alpha = 1.3_dp
      beta = 0.4_dp
      call lp%init(8)
      do ia = 1, 2
         do ib = 1, 2
            ta = ia == 2
            tb = ib == 2
            ! if/else, not merge(): the two shapes differ, which merge() does
            ! not allow (nvfortran segfaults in its transpose).
            if (ta) then
               a = matrix(0, KK, MM)
               op_a = transpose(a)
            else
               a = matrix(0, MM, KK)
               op_a = a
            end if
            if (tb) then
               b = matrix(100, NN, KK)
               op_b = transpose(b)
            else
               b = matrix(100, KK, NN)
               op_b = b
            end if
            c = matrix(200, MM, NN)
            ref = alpha*matmul(op_a, op_b) + beta*c
            !$acc enter data copyin(a, b, c)
            call lp%gemm(OPS(ia), OPS(ib), MM, NN, KK, alpha, a, size(a, 1), b, size(b, 1), beta, c, MM)
            !$acc update self(c)
            !$acc exit data delete(a, b, c)
            call check(error, maxval(abs(c - ref)) < TOL, "gemm " // OPS(ia) // OPS(ib) // " differs from matmul")
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
      end do
      call lp%release()
   end subroutine test_gemm

   subroutine test_gemm_strided(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: MM = 4, NN = 3, KK = 5, NB = 3
      real(dp), allocatable :: a(:), b(:), c(:), ref(:)
      real(dp) :: ai(MM, KK), bi(KK, NN)
      type(trc_linalg_t) :: lp
      integer :: q
      allocate (a(MM*KK*NB), b(KK*NN*NB), c(MM*NN*NB), ref(MM*NN*NB))
      c = 0.0_dp
      do q = 1, NB
         ai = matrix(10*q, MM, KK)
         bi = matrix(50*q, KK, NN)
         a((q - 1)*MM*KK + 1:q*MM*KK) = reshape(ai, [MM*KK])
         b((q - 1)*KK*NN + 1:q*KK*NN) = reshape(bi, [KK*NN])
         ref((q - 1)*MM*NN + 1:q*MM*NN) = reshape(matmul(ai, bi), [MM*NN])
      end do
      call lp%init(8)
      !$acc enter data copyin(a, b, c)
      call lp%gemm_strided("N", "N", MM, NN, KK, 1.0_dp, a, MM, int(MM*KK, int64), b, KK, int(KK*NN, int64), &
                           0.0_dp, c, MM, int(MM*NN, int64), NB)
      !$acc update self(c)
      !$acc exit data delete(a, b, c)
      call lp%release()
      call check(error, maxval(abs(c - ref)) < TOL, "gemm_strided differs from per-matrix matmul")
   end subroutine test_gemm_strided

   subroutine test_transpose(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: MM = 6, NN = 4
      real(dp), allocatable :: a(:, :), c(:, :), s(:, :), sc(:, :), sref(:, :)
      type(trc_linalg_t) :: lp
      integer :: i, j
      a = matrix(3, MM, NN)
      allocate (c(NN, MM))
      c = 0.0_dp
      ! S: upper triangle (diagonal zero) of a square matrix; C = S + S^T.
      s = matrix(5, NN, NN)
      do j = 1, NN
         do i = j, NN
            s(i, j) = 0.0_dp
         end do
      end do
      sref = s + transpose(s)
      allocate (sc(NN, NN))
      sc = 0.0_dp
      call lp%init(8)
      !$acc enter data copyin(a, c, s, sc)
      call lp%transpose(MM, NN, a, MM, c, NN)
      call lp%sym_from_upper(NN, s, NN, sc, NN)
      !$acc update self(c, sc)
      !$acc exit data delete(a, c, s, sc)
      call lp%release()
      call check(error, maxval(abs(c - transpose(a))) < TOL, "transpose differs from intrinsic transpose")
      if (allocated(error)) return
      call check(error, maxval(abs(sc - sref)) < TOL, "sym_from_upper differs from S + S^T")
   end subroutine test_transpose

   subroutine test_gemv(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: MM = 6, NN = 4
      real(dp), allocatable :: a(:, :), x(:), y(:), ref(:)
      type(trc_linalg_t) :: lp
      a = matrix(1, MM, NN)
      call lp%init(8)
      ! y = 1.5 A x + 0.5 y
      x = vec(9, NN)
      y = vec(11, MM)
      ref = 1.5_dp*matmul(a, x) + 0.5_dp*y
      !$acc enter data copyin(a, x, y)
      call lp%gemv("N", MM, NN, 1.5_dp, a, MM, x, 0.5_dp, y)
      !$acc update self(y)
      !$acc exit data delete(a, x, y)
      call check(error, maxval(abs(y - ref)) < TOL, "gemv N differs from matmul")
      if (allocated(error)) then
         call lp%release()
         return
      end if
      ! y = A^T x
      x = vec(13, MM)
      y = vec(15, NN)
      ref = matmul(transpose(a), x)
      !$acc enter data copyin(a, x, y)
      call lp%gemv("T", MM, NN, 1.0_dp, a, MM, x, 0.0_dp, y)
      !$acc update self(y)
      !$acc exit data delete(a, x, y)
      call lp%release()
      call check(error, maxval(abs(y - ref)) < TOL, "gemv T differs from matmul")
   end subroutine test_gemv

   subroutine test_dot(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 37
      real(dp), allocatable :: x(:), y(:)
      real(dp) :: r
      type(trc_linalg_t) :: lp
      x = vec(2, N)
      y = vec(4, N)
      call lp%init(8)
      !$acc enter data copyin(x, y)
      r = lp%dot(N, x, y)
      !$acc exit data delete(x, y)
      call lp%release()
      call check(error, abs(r - dot_product(x, y)) < TOL*max(1.0_dp, abs(r)), "dot differs from dot_product")
   end subroutine test_dot

   subroutine test_trsm(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: MM = 5, NN = 3
      real(dp), allocatable :: l(:, :), u(:, :), b(:, :), x(:, :), b0(:, :)
      type(trc_linalg_t) :: lp
      integer :: i, j
      ! Well-conditioned triangles: unit-ish diagonal dominates.
      l = 0.2_dp*matrix(21, MM, MM)
      u = 0.2_dp*matrix(23, NN, NN)
      do j = 1, MM
         do i = 1, MM
            if (i < j) l(i, j) = 0.0_dp
         end do
         l(j, j) = 2.0_dp + abs(l(j, j))
      end do
      do j = 1, NN
         do i = 1, NN
            if (i > j) u(i, j) = 0.0_dp
         end do
         u(j, j) = 2.0_dp + abs(u(j, j))
      end do
      b0 = matrix(30, MM, NN)
      call lp%init(8)
      ! Left, lower, N: solve L X = 0.7 B.
      b = b0
      !$acc enter data copyin(l, b)
      call lp%trsm("L", "L", "N", MM, NN, 0.7_dp, l, MM, b, MM)
      !$acc update self(b)
      !$acc exit data delete(l, b)
      call check(error, maxval(abs(matmul(l, b) - 0.7_dp*b0)) < TOL, "trsm L,L,N: L X /= alpha B")
      if (allocated(error)) then
         call lp%release()
         return
      end if
      ! Right, upper, T: solve X U^T = B.
      b = b0
      !$acc enter data copyin(u, b)
      call lp%trsm("R", "U", "T", MM, NN, 1.0_dp, u, NN, b, MM)
      !$acc update self(b)
      !$acc exit data delete(u, b)
      call lp%release()
      x = matmul(b, transpose(u))
      call check(error, maxval(abs(x - b0)) < TOL, "trsm R,U,T: X U^T /= B")
   end subroutine test_trsm

   subroutine test_potrf(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 9
      real(dp), allocatable :: a(:, :), a0(:, :), l(:, :)
      type(trc_linalg_t) :: lp
      integer :: i, j
      a0 = spd(N)
      a = a0
      call lp%init(N)
      !$acc enter data copyin(a)
      call lp%potrf(N, a, N)
      !$acc update self(a)
      !$acc exit data delete(a)
      call lp%release()
      allocate (l(N, N))
      l = 0.0_dp
      do j = 1, N
         do i = j, N
            l(i, j) = a(i, j)
         end do
      end do
      call check(error, maxval(abs(matmul(l, transpose(l)) - a0)) < TOL*maxval(abs(a0)), "potrf: L L^T /= A")
   end subroutine test_potrf

   subroutine test_getrs(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 8, NRHS = 3
      real(dp), allocatable :: a(:, :), a0(:, :), b(:, :), b0(:, :)
      integer, allocatable :: ipiv(:)
      integer :: info, i
      type(trc_linalg_t) :: lp
      a0 = matrix(40, N, N)
      do i = 1, N
         a0(i, i) = a0(i, i) + 0.5_dp   ! keep it away from singular
      end do
      a = a0
      b0 = matrix(60, N, NRHS)
      b = b0
      allocate (ipiv(N))
      ipiv = 0
      call lp%init(N)
      !$acc enter data copyin(a, b, ipiv)
      call lp%getrf(N, a, N, ipiv, info)
      call lp%getrs("N", N, NRHS, a, N, ipiv, b, N)
      !$acc update self(b)
      !$acc exit data delete(a, b, ipiv)
      call lp%release()
      call check(error, info == 0, "getrf reported a zero pivot")
      if (allocated(error)) return
      call check(error, maxval(abs(matmul(a0, b) - b0)) < 1.0e-8_dp*maxval(abs(b0)), "getrs: A X /= B")
   end subroutine test_getrs

   subroutine test_syev(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 10
      real(dp), allocatable :: a(:, :), a0(:, :), w(:)
      type(trc_linalg_t) :: lp
      integer :: i
      a0 = matrix(70, N, N)
      a0 = 0.5_dp*(a0 + transpose(a0))
      a = a0
      allocate (w(N))
      w = 0.0_dp
      call lp%init(N)
      !$acc enter data copyin(a) create(w)
      call lp%syev(N, a, N, w)
      !$acc update self(a, w)
      !$acc exit data delete(a, w)
      call lp%release()
      call check(error, maxval(abs(matmul(a0, a) - a*spread(w, 1, N))) < TOL*max(1.0_dp, maxval(abs(w))), &
                 "syev: A V /= V W")
      if (allocated(error)) return
      call check(error, maxval(abs(matmul(transpose(a), a))) > 0.0_dp, "syev: eigenvectors are zero")
      if (allocated(error)) return
      do i = 1, N
         if (i > 1) call check(error, w(i) >= w(i - 1), "syev: eigenvalues not ascending")
         if (allocated(error)) return
      end do
      ! V^T V = I
      call check(error, maxval(abs(matmul(transpose(a), a) - identity(N))) < TOL, "syev: V^T V /= I")
   end subroutine test_syev

   function identity(n) result(e)
      integer, intent(in) :: n
      real(dp) :: e(n, n)
      integer :: i
      e = 0.0_dp
      do i = 1, n
         e(i, i) = 1.0_dp
      end do
   end function identity

end module check_linalg_ref_suite

program check_linalg_ref
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_linalg_ref_suite, only: collect_linalg_ref
   implicit none
   integer :: stat
   stat = 0
   ! Serial: the device libraries and their handles are not thread-safe
   ! across test-drive's OpenMP workers.
   call run_testsuite(collect_linalg_ref, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_linalg_ref
