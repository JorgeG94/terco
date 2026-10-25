!!
!! Dense linear algebra on matrices that live on the device.
!!
!! The SCF keeps its nao x nao matrices resident across iterations, and a
!! Fock build on the device is worth much less if the diagonalisation drags
!! them back every step: metalquicha measured 0.33 s of Fock against 0.92 s
!! of host LAPACK plus copies at 832 functions, and the ratio only worsens
!! with size since the eigenproblem is cubic. So the eigenproblem, the
!! products and the dot products are cuSOLVER and cuBLAS here, reached from
!! OpenACC through `host_data use_device` on arrays the caller has already
!! put on the device.
!!
!! Under a host build the same calls go to pic-blas on the same arrays,
!! which are then simply host arrays, and the directives are comments. One
!! source, two builds, the same numbers to rounding -- which is how every
!! terco driver is arranged, and the reason a CI runner without a GPU can
!! still validate the SCF. pic-blas rather than raw BLAS because its
!! wrappers carry explicit interfaces: nothing in terco is called through
!! an implicit one.
!!
!! Every array argument must be present on the device in the OpenACC build.
!! A host array passed here under -acc is a segfault, not a wrong number,
!! and that is by design: it cannot pass silently.
!!
module trc_linalg
   use, intrinsic :: iso_c_binding, only: c_int, c_size_t
   use pic_types, only: dp, int64, default_int
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
   public :: trc_release_local_mem_reservation

   type :: trc_linalg_t
    !! long lived object that on device enabled builds holds a buffer space and cublas and cusolver handles
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
      ! gemm and dot CANNOT be made into true `generic ::` bindings: every
      ! existing call site that offsets into a larger array -- a scalar
      ! element standing for "the sub-matrix starting here"
      ! (trc_rimp2_driver's tens(1,1,q)), or a higher-rank section fed to
      ! dot's flat x(*) (estore(:,:,:,i), everywhere in the ADIIS code) --
      ! relies on sequence association, and Fortran resolves a GENERIC
      ! reference by exact type/kind/rank agreement only: sequence
      ! association is a property of a resolved reference to one specific
      ! procedure, not of overload resolution, so gfortran (and, by the
      ! same standard rule, nvfortran) rejects the call outright before it
      ! ever reaches the raw specific. Confirmed with a minimal reproducer
      ! (`call gemm('T', big(1,1,2), 3)` inside a two-specific generic
      ! fails to resolve even though only one candidate has the right
      ! first-argument type). So the raw forms keep their own plain names,
      ! and the shape-inferred forms get their OWN names alongside them --
      ! every existing `lp%gemm("N", ...)` / `lp%dot(n, x, y)` call is
      ! completely unaffected.
      procedure :: gemm => linalg_gemm_raw
      procedure :: gemm_mat => linalg_gemm_mat
      procedure :: gemm_strided => linalg_gemm_strided
      procedure :: transpose => linalg_transpose
      procedure :: sym_from_upper => linalg_sym_from_upper
      procedure :: potrf => linalg_potrf
      procedure :: getrf => linalg_getrf
      procedure :: getrs => linalg_getrs
      procedure :: gemv => linalg_gemv
      procedure :: trsm => linalg_trsm
      ! syev has no such caller -- every existing call passes a plain,
      ! exactly-ranked array (never a scalar offset or a rank mismatch) --
      ! so it is safe to make it a real generic.
      procedure, private :: linalg_syev_raw
      procedure, private :: linalg_syev_mat
      generic :: syev => linalg_syev_raw, linalg_syev_mat
      procedure :: dot => linalg_dot_raw
      procedure :: dot_2d => linalg_dot_2d
      procedure :: dot_3d => linalg_dot_3d
      procedure :: density => linalg_density
      procedure :: sync => linalg_sync
   end type trc_linalg_t

contains


   subroutine linalg_init(this, n)
    !! Handles, the solver's workspace for an n x n eigenproblem, and the
    !! library streams set to the OpenACC synchronous queue so a kernel that
    !! wrote a matrix has finished before cuBLAS reads it.
    !!
    !! CREATED ONCE. A second init on a ready object keeps its handles and only
    !! grows the eigensolver's workspace if n is larger than any seen so far,
    !! so one object can be made at the start of a calculation and carried
    !! through all of it -- the SAD atoms, the SCF, the XC quadrature, RI-MP2 --
    !! instead of each of them creating and destroying cuBLAS and cuSOLVER.
      class(trc_linalg_t), intent(inout) :: this
      integer, intent(in) :: n
#ifdef _OPENACC
      integer :: st
#endif
      if (.not. this%ready) then
#ifdef _OPENACC
         st = cublasCreate(this%hb)
         if (st /= 0) error stop "trc_linalg: cublasCreate failed"
         st = cusolverDnCreate(this%hs)
         if (st /= 0) error stop "trc_linalg: cusolverDnCreate failed"
         st = cublasSetStream(this%hb, acc_get_cuda_stream(acc_async_sync))
         st = cusolverDnSetStream(this%hs, acc_get_cuda_stream(acc_async_sync))
         !$acc enter data create(this%devinfo)
#endif
         this%n = 0
         this%ready = .true.
      end if
      if (n > this%n) call grow_workspace(this, n)
   end subroutine linalg_init

   subroutine grow_workspace(this, n)
   !! The dsyevd workspace for n x n, replacing a smaller one.
      class(trc_linalg_t), intent(inout) :: this
      integer, intent(in) :: n
#ifdef _OPENACC
      real(dp), allocatable :: a(:, :), w(:)
      integer :: st
      if (allocated(this%dwork)) then
         !$acc exit data delete(this%dwork)
         deallocate (this%dwork)
      end if
      allocate (a(n, n), w(n))
      !$acc enter data create(a, w)
      !$acc host_data use_device(a, w)
      st = cusolverDnDsyevd_bufferSize(this%hs, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, &
                                       n, a, n, w, this%lwork)
      !$acc end host_data
      !$acc exit data delete(a, w)
      if (st /= 0) error stop "trc_linalg: cusolverDnDsyevd_bufferSize failed"
      allocate (this%dwork(max(this%lwork, 1)))
      !$acc enter data create(this%dwork)
#endif
      this%n = n
   end subroutine grow_workspace

   !
   ! Wait for everything this object's library calls queued on the device.
   ! cuBLAS and cuSOLVER return to the host before their kernels finish, so a
   ! host timer around them measures the launch; this is what a timer calls
   ! before it reads the clock. A no-op on the host build.
   !
   subroutine linalg_sync(this)
    !! GPU synchoronize, TODO(chore): abstract to its own type of procedure this does not belong here
      class(trc_linalg_t), intent(in) :: this
#ifdef _OPENACC
      interface
         function cudaDeviceSynchronize() result(res) bind(c, name="cudaDeviceSynchronize")
            import :: c_int
            implicit none
            integer(c_int) :: res
         end function cudaDeviceSynchronize
      end interface
      integer(c_int) :: st
      if (.not. this%ready) return
      st = cudaDeviceSynchronize()
#else
      if (.not. this%ready) return
#endif
   end subroutine linalg_sync

   subroutine linalg_release(this)
    !! at the end relrease the handles
      class(trc_linalg_t), intent(inout) :: this
#ifdef _OPENACC
      integer :: st
      if (this%ready) then
         if (allocated(this%dwork)) then
            !$acc exit data delete(this%dwork)
            deallocate (this%dwork)
         end if
         !$acc exit data delete(this%devinfo)
         st = cublasDestroy(this%hb)
         st = cusolverDnDestroy(this%hs)
      end if
#endif
      this%ready = .false.
      this%n = 0
   end subroutine linalg_release

   subroutine linalg_gemm_raw(this, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
   !! C = alpha op(A) op(B) + beta C, with 'N' or 'T' for each operand.
   !! if not a GPU enabled build, it will run on the host
      class(trc_linalg_t), intent(in) :: this
      character, intent(in) :: ta, tb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: a(lda, *), b(ldb, *)
      real(dp), intent(inout) :: c(ldc, *)
      logical :: tra, trb
      tra = ta == "T" .or. ta == "t"
      trb = tb == "T" .or. tb == "t"
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
         call pic_gemm(a(1:m, 1:k), b(1:k, 1:n), c(1:m, 1:n), "N", "N", alpha, beta)
      else if (tra .and. .not. trb) then
         call pic_gemm(a(1:k, 1:m), b(1:k, 1:n), c(1:m, 1:n), "T", "N", alpha, beta)
      else if (.not. tra .and. trb) then
         call pic_gemm(a(1:m, 1:k), b(1:n, 1:k), c(1:m, 1:n), "N", "T", alpha, beta)
      else
         call pic_gemm(a(1:k, 1:m), b(1:n, 1:k), c(1:m, 1:n), "T", "T", alpha, beta)
      end if
#endif
   end subroutine linalg_gemm_raw

   ! `count` GEMMs of one shape at once: matrix i of A starts at a(1 + (i-1) sa),
   ! likewise B and C. cublasDgemmStridedBatched on the device; a loop of
   ! pic_gemm on the host. The flat dummies are what lets a caller hand in
   ! the start of a tile inside a larger array.
   !
   ! C = A^T, out of place: A is m x n, C is n x m. cuBLAS's geam does it at
   ! close to the device's bandwidth, which a `do concurrent` copy cannot:
   ! one side of a transpose is always strided across a warp, and only a
   ! tile staged through shared memory hides that.
   !
   subroutine linalg_transpose(this, m, n, a, lda, c, ldc)
    !! transpose using cublasDgeam can be done better, see the transpose in the tools/
    !! however that tranpose needs cuda and we are a do concurrent codebase
      class(trc_linalg_t), intent(in) :: this
      integer, intent(in) :: m, n, lda, ldc
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(inout) :: c(ldc, *)
#ifdef _OPENACC
      integer :: st
      ! beta = 0: the second operand is not read, but cuBLAS still checks
      ! its leading dimension against C's shape, so it is C itself -- the
      ! in-place form geam documents (B = C, ldb = ldc, transb = N).
      !$acc host_data use_device(a, c)
      st = cublasDgeam(this%hb, CUBLAS_OP_T, CUBLAS_OP_N, n, m, 1.0_dp, a, lda, 0.0_dp, c, ldc, c, ldc)
      !$acc end host_data
      if (st /= 0) error stop "trc_linalg: cublasDgeam failed"
#else
      integer :: i, j
      do j = 1, n
         do i = 1, m
            c(j, i) = a(i, j)
         end do
      end do
#endif
   end subroutine linalg_transpose

   !
   ! C = S + S^T, n x n, at the device's bandwidth (the same cuBLAS geam
   ! as `transpose`, alpha = beta = 1 instead of beta = 0): the unpack a
   ! packed-triangular object (trc_scf_df.F90's K block) needs, once the
   ! packed run for each column has already been copied into S's UPPER
   ! triangle (diagonal included, lower zeroed) by a plain, both-sides-
   ! contiguous copy -- see that caller's own comment for why the packed
   ! layout makes that copy trivial where the old (mu, nu) gather was not.
   ! The caller is responsible for S's diagonal: this geam reads it from
   ! BOTH operands (S and S^T agree there), so a non-zero diagonal comes
   ! out doubled -- zero it in S first and add it back after, or halve
   ! C's diagonal afterward; either is the caller's O(n) problem, not
   ! this routine's.
   !
   subroutine linalg_sym_from_upper(this, n, s, lds, c, ldc)
      class(trc_linalg_t), intent(in) :: this
      integer, intent(in) :: n, lds, ldc
      real(dp), intent(in) :: s(lds, *)
      real(dp), intent(inout) :: c(ldc, *)
#ifdef _OPENACC
      integer :: st
      !$acc host_data use_device(s, c)
      st = cublasDgeam(this%hb, CUBLAS_OP_T, CUBLAS_OP_N, n, n, 1.0_dp, s, lds, 1.0_dp, s, lds, c, ldc)
      !$acc end host_data
      if (st /= 0) error stop "trc_linalg: cublasDgeam failed"
#else
      integer :: i, j
      do concurrent(j=1:n, i=1:n)
         c(i, j) = s(i, j) + s(j, i)
      end do
#endif
   end subroutine linalg_sym_from_upper

   subroutine linalg_gemm_strided(this, ta, tb, m, n, k, alpha, a, lda, sa, b, ldb, sb, beta, c, ldc, sc, count)
      class(trc_linalg_t), intent(in) :: this
      character, intent(in) :: ta, tb
      integer, intent(in) :: m, n, k, lda, ldb, ldc, count
      integer(int64), intent(in) :: sa, sb, sc
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: a(*), b(*)
      real(dp), intent(inout) :: c(*)
      logical :: tra, trb
      tra = ta == "T" .or. ta == "t"
      trb = tb == "T" .or. tb == "t"
      if (count < 1) return
#ifdef _OPENACC
      block
         integer :: st, oa, ob
         oa = merge(CUBLAS_OP_T, CUBLAS_OP_N, tra)
         ob = merge(CUBLAS_OP_T, CUBLAS_OP_N, trb)
         !$acc host_data use_device(a, b, c)
         st = cublasDgemmStridedBatched(this%hb, oa, ob, m, n, k, alpha, a, lda, sa, b, ldb, sb, &
                                        beta, c, ldc, sc, count)
         !$acc end host_data
         if (st /= 0) error stop "trc_linalg: cublasDgemmStridedBatched failed"
      end block
#else
      block
         integer :: i
         integer(int64) :: oa, ob, oc
         do i = 1, count
            oa = int(i - 1, int64)*sa
            ob = int(i - 1, int64)*sb
            oc = int(i - 1, int64)*sc
            call this%gemm(ta, tb, m, n, k, alpha, a(oa + 1), lda, b(ob + 1), ldb, beta, c(oc + 1), ldc)
         end do
      end block
#endif
   end subroutine linalg_gemm_strided

   ! Cholesky factor of a symmetric positive definite A: the lower triangle
   ! is overwritten with L, A = L L^T. The upper triangle is left as it was.
   subroutine linalg_potrf(this, n, a, lda)
      class(trc_linalg_t), intent(inout) :: this
      integer, intent(in) :: n, lda
      real(dp), intent(inout) :: a(lda, *)
#ifdef _OPENACC
      integer :: st, lw, info
      real(dp), allocatable :: work(:)
      !$acc host_data use_device(a)
      st = cusolverDnDpotrf_bufferSize(this%hs, CUBLAS_FILL_MODE_LOWER, n, a, lda, lw)
      !$acc end host_data
      if (st /= 0) error stop "trc_linalg: cusolverDnDpotrf_bufferSize failed"
      allocate (work(max(lw, 1)))
      !$acc data create(work)
      !$acc host_data use_device(a, work, this%devinfo)
      st = cusolverDnDpotrf(this%hs, CUBLAS_FILL_MODE_LOWER, n, a, lda, work, lw, this%devinfo)
      !$acc end host_data
      !$acc end data
      if (st /= 0) error stop "trc_linalg: cusolverDnDpotrf failed"
      !$acc update self(this%devinfo)
      info = this%devinfo
      if (info /= 0) error stop "trc_linalg: cusolverDnDpotrf: matrix not positive definite"
#else
      integer :: info
      interface
         subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            implicit none
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
         end subroutine dpotrf
      end interface
      call dpotrf("L", n, a, lda, info)
      if (info /= 0) error stop "trc_linalg: dpotrf: matrix not positive definite"
#endif
   end subroutine linalg_potrf

   subroutine linalg_getrf(this, n, a, lda, ipiv, info)
   !! LU factorisation with partial pivoting, A = P L U, in place; `ipiv` is
   !! LAPACK's (1-based rows). The pivots live where A does: on the device
   !! under the OpenACC build, so `ipiv` must be present there too.
      class(trc_linalg_t), intent(inout) :: this
      integer, intent(in) :: n, lda
      real(dp), intent(inout) :: a(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info   !! 0, or the first zero pivot as LAPACK reports it
#ifdef _OPENACC
      integer :: st, lw
      real(dp), allocatable :: work(:)
      !$acc host_data use_device(a)
      st = cusolverDnDgetrf_bufferSize(this%hs, n, n, a, lda, lw)
      !$acc end host_data
      if (st /= 0) error stop "trc_linalg: cusolverDnDgetrf_bufferSize failed"
      allocate (work(max(lw, 1)))
      !$acc data create(work)
      !$acc host_data use_device(a, work, ipiv, this%devinfo)
      st = cusolverDnDgetrf(this%hs, n, n, a, lda, work, ipiv, this%devinfo)
      !$acc end host_data
      !$acc end data
      if (st /= 0) error stop "trc_linalg: cusolverDnDgetrf failed"
      !$acc update self(this%devinfo)
      info = this%devinfo
#else
      ! TODO(chore): should use pic_lapack_interfaces
      interface
         subroutine dgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            implicit none
            integer, intent(in) :: m, n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
         end subroutine dgetrf
      end interface
      call dgetrf(n, n, a, lda, ipiv, info)
#endif
   end subroutine linalg_getrf

   ! Solve op(A) X = B ('N' or 'T') from getrf's factors; B (n x nrhs) is
   ! overwritten with X.
   subroutine linalg_getrs(this, trans, n, nrhs, a, lda, ipiv, b, ldb)
      class(trc_linalg_t), intent(inout) :: this
      character, intent(in) :: trans
      integer, intent(in) :: n, nrhs, lda, ldb
      real(dp), intent(in) :: a(lda, *)
      integer, intent(in) :: ipiv(*)
      real(dp), intent(inout) :: b(ldb, *)
#ifdef _OPENACC
      integer :: st, ct
      ct = merge(CUBLAS_OP_T, CUBLAS_OP_N, trans == "T" .or. trans == "t")
      !$acc host_data use_device(a, ipiv, b, this%devinfo)
      st = cusolverDnDgetrs(this%hs, ct, n, nrhs, a, lda, ipiv, b, ldb, this%devinfo)
      !$acc end host_data
      if (st /= 0) error stop "trc_linalg: cusolverDnDgetrs failed"
#else
      integer :: info
      interface
         subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            implicit none
            character, intent(in) :: trans
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            integer, intent(in) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
         end subroutine dgetrs
      end interface
      call dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
      if (info /= 0) error stop "trc_linalg: dgetrs failed"
#endif
   end subroutine linalg_getrs

   ! y = alpha op(A) x + beta y, A m x n, op 'N' or 'T'.
   subroutine linalg_gemv(this, trans, m, n, alpha, a, lda, x, beta, y)
      class(trc_linalg_t), intent(in) :: this
      character, intent(in) :: trans
      integer, intent(in) :: m, n, lda
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: a(lda, *), x(*)
      real(dp), intent(inout) :: y(*)
#ifdef _OPENACC
      integer :: st, ct
      ct = merge(CUBLAS_OP_T, CUBLAS_OP_N, trans == "T" .or. trans == "t")
      !$acc host_data use_device(a, x, y)
      st = cublasDgemv_v2(this%hb, ct, m, n, alpha, a, lda, x, 1, beta, y, 1)
      !$acc end host_data
      if (st /= 0) error stop "trc_linalg: cublasDgemv failed"
#else
      interface
         subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: dp
            implicit none
            character, intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            real(dp), intent(in) :: alpha, beta
            real(dp), intent(in) :: a(lda, *), x(*)
            real(dp), intent(inout) :: y(*)
         end subroutine dgemv
      end interface
      call dgemv(trans, m, n, alpha, a, lda, x, 1, beta, y, 1)
#endif
   end subroutine linalg_gemv

   ! Triangular solve: op(A) X = alpha B ('L') or X op(A) = alpha B ('R'),
   ! A lower ('L') or upper ('U') triangular, op 'N' or 'T', unit diagonal
   ! never. B (m x n) is overwritten with X.
   subroutine linalg_trsm(this, side, uplo, trans, m, n, alpha, a, lda, b, ldb)
      class(trc_linalg_t), intent(in) :: this
      character, intent(in) :: side, uplo, trans
      integer, intent(in) :: m, n, lda, ldb
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(lda, *)
      real(dp), intent(inout) :: b(ldb, *)
#ifdef _OPENACC
      integer :: st, cs, cu, ct
      cs = merge(CUBLAS_SIDE_RIGHT, CUBLAS_SIDE_LEFT, side == "R" .or. side == "r")
      cu = merge(CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER, uplo == "U" .or. uplo == "u")
      ct = merge(CUBLAS_OP_T, CUBLAS_OP_N, trans == "T" .or. trans == "t")
      !$acc host_data use_device(a, b)
      st = cublasDtrsm_v2(this%hb, cs, cu, ct, CUBLAS_DIAG_NON_UNIT, m, n, alpha, a, lda, b, ldb)
      !$acc end host_data
      if (st /= 0) error stop "trc_linalg: cublasDtrsm failed"
#else
      interface
         subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            import :: dp
            implicit none
            character, intent(in) :: side, uplo, transa, diag
            integer, intent(in) :: m, n, lda, ldb
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
         end subroutine dtrsm
      end interface
      call dtrsm(side, uplo, trans, "N", m, n, alpha, a, lda, b, ldb)
#endif
   end subroutine linalg_trsm

   ! Symmetric eigenproblem: A is overwritten with its eigenvectors, W gets
   ! the eigenvalues ascending. Upper triangle is read.
   subroutine linalg_syev_raw(this, n, a, lda, w)
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
      call pic_syev(a(1:n, 1:n), w(1:n), "V", "U", info)
      if (info /= 0) error stop "trc_linalg: pic_syev failed"
#endif
   end subroutine linalg_syev_raw

   ! x . y over n elements, returned to the host.
   function linalg_dot_raw(this, n, x, y) result(r)
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
   end function linalg_dot_raw

   !!
   !! Shape-inferred matrix product: C = alpha op(A) op(B) + beta C, where
   !! transa/transb ("N" by default) say whether A/B enter transposed and
   !! alpha/beta (1 and 0 by default) match BLAS's. m, n and k are read off
   !! the shapes of A and B after the requested transpose, not taken from
   !! the caller, so a conformant call needs nothing beyond the matrices
   !! themselves; a shape mismatch is an `error stop`, not a wrong number
   !! silently read past the end of an array.
   !!
   !! A, B and C must be `contiguous`. A non-contiguous actual argument
   !! (a strided section, or a section that is not a leading run of whole
   !! columns) makes the compiler materialise a HOST temporary for the
   !! dummy, and the `host_data use_device` below then hands cuBLAS a host
   !! address -- a segfault, not a wrong answer. Pass whole matrices, or a
   !! contiguous section such as `fock(:,:,s)` or `estore(:,:,:,k)`.
   !!
   !! Reaches the exact same library entry point as the raw gemm --
   !! pic_gemm on the host, cublasDgemm_v2 on the device -- with the shapes
   !! computed here, so the result is bitwise identical to the equivalent
   !! raw call.
   !!
   subroutine linalg_gemm_mat(this, a, b, c, transa, transb, alpha, beta)
      class(trc_linalg_t), intent(in) :: this
      real(dp), intent(in), contiguous :: a(:, :), b(:, :)
      real(dp), intent(inout), contiguous :: c(:, :)
      character, intent(in), optional :: transa, transb
      real(dp), intent(in), optional :: alpha, beta
      character :: ta, tb
      logical :: tra, trb
      integer :: m, n, k, lda, ldb, ldc
      real(dp) :: al, be

      ta = "N"
      if (present(transa)) ta = transa
      tb = "N"
      if (present(transb)) tb = transb
      al = 1.0_dp
      if (present(alpha)) al = alpha
      be = 0.0_dp
      if (present(beta)) be = beta
      tra = ta == "T" .or. ta == "t"
      trb = tb == "T" .or. tb == "t"

      ! op(A) is m x k; A is stored k x m when transposed, m x k otherwise.
      if (tra) then
         k = size(a, 1)
         m = size(a, 2)
      else
         m = size(a, 1)
         k = size(a, 2)
      end if
      ! op(B) is k x n; B is stored n x k when transposed, k x n otherwise.
      if (trb) then
         n = size(b, 1)
         if (size(b, 2) /= k) error stop "trc_linalg: gemm_mat: A and B do not conform"
      else
         n = size(b, 2)
         if (size(b, 1) /= k) error stop "trc_linalg: gemm_mat: A and B do not conform"
      end if
      if (size(c, 1) /= m .or. size(c, 2) /= n) error stop "trc_linalg: gemm_mat: C does not conform to A, B"

      lda = size(a, 1)
      ldb = size(b, 1)
      ldc = size(c, 1)
      call this%gemm(ta, tb, m, n, k, al, a, lda, b, ldb, be, c, ldc)
   end subroutine linalg_gemm_mat

   !!
   !! x . y over the whole of two same-shaped, contiguous rank-2 arrays,
   !! n = size(x). Same device/host path as the raw dot: the raw form's
   !! assumed-size dummies accept a contiguous array of any rank by
   !! sequence association, so this simply counts the elements and hands
   !! both arrays to it unchanged -- bitwise identical to calling the raw
   !! dot directly with n = size(x).
   !!
   function linalg_dot_2d(this, x, y) result(r)
      class(trc_linalg_t), intent(in) :: this
      real(dp), intent(in), contiguous :: x(:, :), y(:, :)
      real(dp) :: r
      r = this%dot(size(x), x, y)
   end function linalg_dot_2d

   !! As linalg_dot_2d, for contiguous rank-3 arrays -- e.g. a section of a
   !! rank-4 array such as `estore(:,:,:,k)`.
   function linalg_dot_3d(this, x, y) result(r)
      class(trc_linalg_t), intent(in) :: this
      real(dp), intent(in), contiguous :: x(:, :, :), y(:, :, :)
      real(dp) :: r
      r = this%dot(size(x), x, y)
   end function linalg_dot_3d

   !!
   !! Symmetric eigenproblem, shape-inferred: A is square and `contiguous`,
   !! n = size(a,1). A is overwritten with its eigenvectors, W gets the
   !! eigenvalues ascending, exactly as the raw syev leaves them -- this
   !! only reads n and lda off A's shape before calling it, so the result
   !! is bitwise identical to the equivalent raw call.
   !!
   !! A non-contiguous actual argument would force a HOST temporary and
   !! then hand the GPU a host address under `host_data use_device` inside
   !! the raw call -- a segfault. Pass a whole matrix or a contiguous
   !! section.
   !!
   subroutine linalg_syev_mat(this, a, w)
      class(trc_linalg_t), intent(inout) :: this
      real(dp), intent(inout), contiguous :: a(:, :)
      real(dp), intent(out), contiguous :: w(:)
      integer :: n
      n = size(a, 1)
      call this%linalg_syev_raw(n, a, n, w)
   end subroutine linalg_syev_mat

   !!
   !! The density matrix built from the occupied columns of C: D = occ *
   !! C(:,1:k) C(:,1:k)^T, one GEMM with inner dimension k. For k <= 0
   !! (nothing occupied) D is zeroed by a `do concurrent` here instead --
   !! there is no GEMM that would do that, and the zeroing has to run on
   !! whatever device D already lives on, which `do concurrent` does on its
   !! own. C and D are contiguous; C need only have at least k columns.
   !!
   subroutine linalg_density(this, c, k, occ, d)
      class(trc_linalg_t), intent(in) :: this
      real(dp), intent(in), contiguous :: c(:, :)
      integer, intent(in) :: k
      real(dp), intent(in) :: occ
      real(dp), intent(inout), contiguous :: d(:, :)
      integer :: n, i, j
      n = size(d, 1)
      if (k <= 0) then
         do concurrent(j=1:size(d, 2), i=1:n)
            d(i, j) = 0.0_dp
         end do
      else
         call this%gemm("N", "T", n, n, k, occ, c(:, 1:k), n, c(:, 1:k), n, 0.0_dp, d, n)
      end if
   end subroutine linalg_density

   !
   ! The CUDA driver's per-thread "local" stack, reserved across the WHOLE
   ! device (every resident thread the GPU could ever run concurrently, not
   ! what the kernel that triggered it actually launched) the first time any
   ! kernel with a frame bigger than the current limit runs -- and never
   ! given back on its own. Measured here: `trc_pc_k2222_pcs2222` (the
   ! (dd|dd) four-centre class) has the largest stack frame in this whole
   ! library, 59568 bytes (`cuobjdump -res-usage`); on a V100 (80 SMs x 2048
   ! threads/SM = 163840 resident threads) that reservation is
   ! 163840*59568 = 9.09 GiB, matching the ~8.9 GiB this branch measured for
   ! the SAD guess's first non-DF Fock build (SAD always needs this class --
   ! cc-pVDZ oxygen's own d shell against itself -- regardless of `--df`).
   ! Nothing after the guess launches that four-centre path under `--df`
   ! (build_g never calls `eri%fock_resident` there), so MOST of the
   ! reservation -- everything above what the DF path's own largest kernel
   ! still needs, see `SMALL_STACK` below -- is pure waste for the rest of
   ! a `--df` run. `cudaDeviceSetLimit` on `cudaLimitStackSize` shrinks it
   ! back down immediately to that smaller, still-sufficient size; a later
   ! kernel that genuinely needs a bigger frame than THAT grows it again on
   ! its own next launch, at that kernel's own one-time cost, not this
   ! one's -- which is exactly the failure this function's own SMALL_STACK
   ! choice has to avoid: regrowing at the worst possible moment, with
   ! almost no free memory left, is a crash instead of a reservation.
   !
   ! A no-op on the host build (no device, nothing to reserve).
   !
   subroutine trc_release_local_mem_reservation()
#ifdef _OPENACC
      interface
         function cudadevicesetlimit(limit, lim_value) bind(c, name="cudaDeviceSetLimit") result(ierr)
            import :: c_int, c_size_t
            implicit none
            integer(c_int), value :: limit
            integer(c_size_t), value :: lim_value
            integer(c_int) :: ierr
         end function cudadevicesetlimit
      end interface
      integer(c_int), parameter :: CUDA_LIMIT_STACK_SIZE = 0
      integer(c_size_t), parameter :: SMALL_STACK = 32768_c_size_t
      !! NOT CUDA's 1 KiB default: the three-centre DF kernels
      !! (trc_df_3c_k223, measured 28784 B/thread -- the next generation
      !! above it, the def2-universal-jkfit g-function classes, could be
      !! larger still) run AFTER this reset, inside the very fit this
      !! reset is clearing room for, and dropping below their own need
      !! only moves the regrowth -- at a moment chosen to have the LEAST
      !! free memory of the whole run, right after the fit's own setup
      !! block has claimed everything choose_setup_block thought was safe
      !! -- to a cuLaunchKernel OOM instead of a cuMemAlloc one (measured:
      !! it did, at trc_df_3c_k223, with 1024 here). 32 KiB clears every
      !! three-centre and two-centre DF kernel class measured so far with
      !! margin, while still giving back the four-centre classes' own
      !! excess (the next ones down from k2222's 59568 B are the other
      !! four-centre p/d combinations, all still well above 32 KiB).
      integer(c_int) :: ierr
      ierr = cudadevicesetlimit(CUDA_LIMIT_STACK_SIZE, SMALL_STACK)
      if (ierr /= 0) error stop "trc_linalg: cudaDeviceSetLimit(cudaLimitStackSize) failed"
#endif
   end subroutine trc_release_local_mem_reservation

end module trc_linalg
