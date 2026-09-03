!
! Dense linear algebra on matrices that live on the device.
!
! The SCF keeps its nao x nao matrices resident across iterations, and a
! Fock build on the device is worth much less if the diagonalisation drags
! them back every step: metalquicha measured 0.33 s of Fock against 0.92 s
! of host LAPACK plus copies at 832 functions, and the ratio only worsens
! with size since the eigenproblem is cubic. So the eigenproblem, the
! products and the dot products are cuSOLVER and cuBLAS here, reached from
! OpenACC through `host_data use_device` on arrays the caller has already
! put on the device.
!
! Under a host build the same calls go to pic-blas on the same arrays,
! which are then simply host arrays, and the directives are comments. One
! source, two builds, the same numbers to rounding -- which is how every
! terco driver is arranged, and the reason a CI runner without a GPU can
! still validate the SCF. pic-blas rather than raw BLAS because its
! wrappers carry explicit interfaces: nothing in terco is called through
! an implicit one.
!
! Every array argument must be present on the device in the OpenACC build.
! A host array passed here under -acc is a segfault, not a wrong number,
! and that is by design: it cannot pass silently.
!
module trc_linalg
   use trc_boys, only: dp
   use pic_types, only: default_int
   use pic_blas_interfaces, only: pic_gemm, pic_dot
   use pic_lapack_interfaces, only: pic_syev
#ifdef _OPENACC
   use cublas
   use cusolverdn
   use openacc, only: acc_get_cuda_stream, acc_async_sync
#endif
   implicit none
   private

   public :: trc_linalg_t

   type :: trc_linalg_t
      integer :: n = 0
      logical :: ready = .false.
#ifdef _OPENACC
      type(cublasHandle) :: hb
      type(cusolverDnHandle) :: hs
      real(dp), allocatable :: dwork(:)
      integer :: lwork = 0
      integer :: devinfo = 0
#endif
   contains
      procedure :: init => linalg_init
      procedure :: release => linalg_release
      procedure :: gemm => linalg_gemm
      procedure :: syev => linalg_syev
      procedure :: dot => linalg_dot
   end type trc_linalg_t

contains

   !
   ! Handles, the solver's workspace for an n x n eigenproblem, and the
   ! library streams set to the OpenACC synchronous queue so a kernel that
   ! wrote a matrix has finished before cuBLAS reads it.
   !
   subroutine linalg_init(this, n)
      class(trc_linalg_t), intent(inout) :: this
      integer, intent(in) :: n
#ifdef _OPENACC
      real(dp), allocatable :: a(:, :), w(:)
      integer :: st
#endif
      call this%release()
      this%n = n
#ifdef _OPENACC
      st = cublasCreate(this%hb)
      if (st /= 0) error stop "trc_linalg: cublasCreate failed"
      st = cusolverDnCreate(this%hs)
      if (st /= 0) error stop "trc_linalg: cusolverDnCreate failed"
      st = cublasSetStream(this%hb, acc_get_cuda_stream(acc_async_sync))
      st = cusolverDnSetStream(this%hs, acc_get_cuda_stream(acc_async_sync))
      allocate (a(n, n), w(n))
      !$acc enter data create(a, w)
      !$acc host_data use_device(a, w)
      st = cusolverDnDsyevd_bufferSize(this%hs, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, &
                                       n, a, n, w, this%lwork)
      !$acc end host_data
      !$acc exit data delete(a, w)
      if (st /= 0) error stop "trc_linalg: cusolverDnDsyevd_bufferSize failed"
      allocate (this%dwork(max(this%lwork, 1)))
      !$acc enter data create(this%dwork, this%devinfo)
#endif
      this%ready = .true.
   end subroutine linalg_init

   subroutine linalg_release(this)
      class(trc_linalg_t), intent(inout) :: this
#ifdef _OPENACC
      integer :: st
      if (this%ready) then
         !$acc exit data delete(this%dwork, this%devinfo)
         deallocate (this%dwork)
         st = cublasDestroy(this%hb)
         st = cusolverDnDestroy(this%hs)
      end if
#endif
      this%ready = .false.
      this%n = 0
   end subroutine linalg_release

   ! C = alpha op(A) op(B) + beta C, with 'N' or 'T' for each operand.
   subroutine linalg_gemm(this, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      class(trc_linalg_t), intent(in) :: this
      character, intent(in) :: ta, tb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: a(lda, *), b(ldb, *)
      real(dp), intent(inout) :: c(ldc, *)
      logical :: tra, trb
      tra = ta == 'T' .or. ta == 't'
      trb = tb == 'T' .or. tb == 't'
#ifdef _OPENACC
      block
         integer :: st, oa, ob
         oa = merge(CUBLAS_OP_T, CUBLAS_OP_N, tra)
         ob = merge(CUBLAS_OP_T, CUBLAS_OP_N, trb)
         !$acc host_data use_device(a, b, c)
         st = cublasDgemm_v2(this%hb, oa, ob, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
         !$acc end host_data
         if (st /= 0) error stop "trc_linalg: cublasDgemm failed"
      end block
#else
      ! op(A) is m x k and op(B) is k x n; the stored block of a transposed
      ! operand has those extents the other way round.
      if (.not. tra .and. .not. trb) then
         call pic_gemm(a(1:m, 1:k), b(1:k, 1:n), c(1:m, 1:n), 'N', 'N', alpha, beta)
      else if (tra .and. .not. trb) then
         call pic_gemm(a(1:k, 1:m), b(1:k, 1:n), c(1:m, 1:n), 'T', 'N', alpha, beta)
      else if (.not. tra .and. trb) then
         call pic_gemm(a(1:m, 1:k), b(1:n, 1:k), c(1:m, 1:n), 'N', 'T', alpha, beta)
      else
         call pic_gemm(a(1:k, 1:m), b(1:n, 1:k), c(1:m, 1:n), 'T', 'T', alpha, beta)
      end if
#endif
   end subroutine linalg_gemm

   ! Symmetric eigenproblem: A is overwritten with its eigenvectors, W gets
   ! the eigenvalues ascending. Upper triangle is read.
   subroutine linalg_syev(this, n, a, lda, w)
      class(trc_linalg_t), intent(inout) :: this
      integer, intent(in) :: n, lda
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*)
#ifdef _OPENACC
      integer :: st, info
      if (n > this%n) error stop "trc_linalg: eigenproblem larger than the workspace"
      !$acc host_data use_device(a, w, this%dwork, this%devinfo)
      st = cusolverDnDsyevd(this%hs, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, a, lda, w, &
                            this%dwork, this%lwork, this%devinfo)
      !$acc end host_data
      if (st /= 0) error stop "trc_linalg: cusolverDnDsyevd failed"
      ! Non-convergence is reported through the device int, not the status.
      !$acc update self(this%devinfo)
      info = this%devinfo
      if (info /= 0) error stop "trc_linalg: cusolverDnDsyevd did not converge"
#else
      integer(default_int) :: info
      call pic_syev(a(1:n, 1:n), w(1:n), 'V', 'U', info)
      if (info /= 0) error stop "trc_linalg: pic_syev failed"
#endif
   end subroutine linalg_syev

   ! x . y over n elements, returned to the host.
   function linalg_dot(this, n, x, y) result(r)
      class(trc_linalg_t), intent(in) :: this
      integer, intent(in) :: n
      real(dp), intent(in) :: x(*), y(*)
      real(dp) :: r
#ifdef _OPENACC
      integer :: st
      !$acc host_data use_device(x, y)
      st = cublasDdot_v2(this%hb, n, x, 1, y, 1, r)
      !$acc end host_data
      if (st /= 0) error stop "trc_linalg: cublasDdot failed"
#else
      r = pic_dot(x(1:n), y(1:n))
#endif
   end function linalg_dot

end module trc_linalg
