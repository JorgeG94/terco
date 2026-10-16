! allowed in the repo, do not delete, use in case of emergency performance of transpose
!
! A hand-written CUDA Fortran transpose that beats cublasDgeam.
!
! BREAK IN CASE OF PERFORMANCE. This is not part of terco: its kernels are
! standard `do concurrent` by rule, and nothing builds this file. It is kept
! because it was measured, and it is the answer if terco ever takes an
! optional CUDA backend.
!
! WHAT IT DOES
! ------------
! C = A^T out of place, double precision: A is m x n with leading dimension
! lda, C is n x m with leading dimension ldc -- the operation trc_linalg's
! `transpose` gets from cublasDgeam(OP_T, OP_N, ...), which RI-MP2's pair
! step calls once per occupied orbital on a tall, narrow block of K.
!
! THE KERNEL
! ----------
! The classic shared-memory transpose: a block reads a 32 x 32 tile of A
! with consecutive threads on consecutive rows (coalesced), and writes it
! out transposed from shared memory, again with consecutive threads on
! consecutive rows of C. 32 x 8 threads, four rows each. The tile is padded
! to 33 by habit; shared memory has ~16x the bandwidth the copy needs, so
! the 2-way bank conflict it avoids was measured not to matter.
!
! What made it fast is the GRID ORIENTATION: the long n direction of A
! (the columns, whose tiles become C's row bands) on blockIdx%x, the m
! direction on blockIdx%y. Blocks are scheduled x-fastest, so the ~640
! resident blocks then share ONE m-tile and cover ~640 consecutive n-tiles:
! 32 rows of C, each written in a contiguous ~160 KB run. The other way
! round they spread over every m-tile and few n-tiles, and C is written in
! ~6.6 KB pieces across all its rows. (The likely reason, from the access
! pattern; not confirmed with ncu.)
!
! MEASURED (V100-DGXS-32GB, GB/s = 2 m n 8 / time, median of 12, bitwise
! equal to cublasDgeam on every shape incl. a 1002 x 7778 edge case)
!
!   shape (m x n, lda)            geam   n on y (old)   n on x (this)   cudaMemcpy
!   800 x 128000, 3200 (w40)       700        775             794            818
!   1200 x 288000, 4800 (w60)      694        771             796            820
!   4096 x 4096, 4096              773        715             806            799
!   8192 x 8192, 8192              785        776             817            815
!
! (From the tuning harness. The driver below times with CUDA events and
! reproduces these on the large shapes; at 4096 a call is ~0.3 ms and it
! reads ~766 against geam's ~747.)
!
! On RI-MP2's shapes this is the bandwidth a `do concurrent` copy of the
! same bytes reaches; the best pure `do concurrent` TRANSPOSE manages
! 147-294 GB/s there (one side is always strided across a warp, and do
! concurrent has no shared memory). At w60 the pair step's ~666 GB of
! transposes would take 0.84 s instead of geam's 0.96.
!
! Tried and NOT worth it (all within noise or worse): 16-byte (double2)
! loads/stores, a bounds-check-free path for interior tiles, streaming
! stores, read-only loads (nvfortran already emits ld.global.nc for an
! intent(in) device dummy), 32x16/32x32/32x4 thread blocks, persistent and
! software-pipelined grid-stride kernels, grouped rasterisation.
!
! BUILD AND RUN (the driver below checks and times it against geam)
!
!   nvfortran -O3 -cuda -gpu=cc70 -cudalib=cublas fast_transpose.F90 -o fast_transpose
!   CUDA_VISIBLE_DEVICES=<gpu> ./fast_transpose
!
module fast_transpose_m
   use cudafor, only: dim3, syncthreads, threadIdx, blockIdx
   implicit none
   private

   public :: fast_transpose

   integer, parameter :: dp = selected_real_kind(15, 307)
   integer, parameter :: TILE = 32   !! tile edge
   integer, parameter :: ROWS = 8    !! thread rows per block; each does TILE/ROWS rows

contains

   attributes(global) subroutine k_transpose(a, lda, c, ldc, m, n)
      !! One TILE x TILE tile of C = A^T. The tile's n index comes from
      !! blockIdx%x and its m index from blockIdx%y -- see the header.
      integer, value :: lda, ldc, m, n
      real(dp), device, intent(in) :: a(lda, *)
      real(dp), device, intent(inout) :: c(ldc, *)
      real(dp), shared :: tl(TILE + 1, TILE)
      integer :: tx, ty, k, i, j, bi, bj
      tx = threadIdx%x
      ty = threadIdx%y
      bi = (blockIdx%y - 1)*TILE
      bj = (blockIdx%x - 1)*TILE
      i = bi + tx
      do k = 0, TILE - ROWS, ROWS
         j = bj + ty + k
         if (i <= m .and. j <= n) tl(tx, ty + k) = a(i, j)
      end do
      call syncthreads()
      j = bj + tx
      do k = 0, TILE - ROWS, ROWS
         i = bi + ty + k
         if (j <= n .and. i <= m) c(j, i) = tl(ty + k, tx)
      end do
   end subroutine k_transpose

   subroutine fast_transpose(a, lda, c, ldc, m, n)
      !! C = A^T for device arrays: A is m x n (leading dimension lda), C is
      !! n x m (leading dimension ldc). Asynchronous on the default stream.
      real(dp), device, intent(in) :: a(lda, *)
      real(dp), device, intent(inout) :: c(ldc, *)
      integer, intent(in) :: lda, ldc, m, n
      type(dim3) :: grid, blk
      grid = dim3((n + TILE - 1)/TILE, (m + TILE - 1)/TILE, 1)
      blk = dim3(TILE, ROWS, 1)
      call k_transpose<<<grid, blk>>>(a, lda, c, ldc, m, n)
   end subroutine fast_transpose

end module fast_transpose_m

program fast_transpose_check
   !! Checks fast_transpose bitwise against cublasDgeam and times both, on
   !! RI-MP2's shapes and two square ones, plus an edge shape that is not a
   !! multiple of the tile.
   use cudafor, only: cudaEvent, cudaEventCreate, cudaEventRecord, cudaEventSynchronize, &
                      cudaEventElapsedTime
   use cublas, only: cublasHandle, cublasCreate, cublasDestroy, cublasDgeam, CUBLAS_OP_T, CUBLAS_OP_N
   use fast_transpose_m, only: fast_transpose
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 307)
   integer, parameter :: sp = selected_real_kind(6, 37)   !! CUDA event times are single precision
   integer, parameter :: NREP = 12
   call run(800, 128000, 3200, 128000, "w40 RI-MP2 block")
   call run(1200, 288000, 4800, 288000, "w60 RI-MP2 block")
   call run(4096, 4096, 4096, 4096, "square 4096")
   call run(8192, 8192, 8192, 8192, "square 8192")
   call run(1002, 7778, 1010, 7780, "edge 1002 x 7778")

contains

   subroutine run(m, n, lda, ldc, label)
      integer, intent(in) :: m, n, lda, ldc
      character(len=*), intent(in) :: label
      real(dp), allocatable :: a(:, :), c1(:, :), c2(:, :)
      real(dp), device, allocatable :: ad(:, :), cd(:, :)
      type(cublasHandle) :: h
      type(cudaEvent) :: e0, e1
      real(sp) :: ms
      real(sp) :: tg(NREP), tf(NREP)
      integer :: i, j, r, st
      allocate (a(lda, n), c1(ldc, m), c2(ldc, m))
      do j = 1, n
         do i = 1, lda
            a(i, j) = real(i, dp) + 1.0e-6_dp*real(j, dp)
         end do
      end do
      allocate (ad(lda, n), cd(ldc, m))
      ad = a
      st = cublasCreate(h)
      st = cudaEventCreate(e0)
      st = cudaEventCreate(e1)
      do r = 0, NREP
         st = cudaEventRecord(e0, 0)
         st = cublasDgeam(h, CUBLAS_OP_T, CUBLAS_OP_N, n, m, 1.0_dp, ad, lda, 0.0_dp, cd, ldc, cd, ldc)
         st = cudaEventRecord(e1, 0)
         st = cudaEventSynchronize(e1)
         st = cudaEventElapsedTime(ms, e0, e1)
         if (r > 0) tg(r) = ms
      end do
      c1 = cd
      cd = 0.0_dp
      do r = 0, NREP
         st = cudaEventRecord(e0, 0)
         call fast_transpose(ad, lda, cd, ldc, m, n)
         st = cudaEventRecord(e1, 0)
         st = cudaEventSynchronize(e1)
         st = cudaEventElapsedTime(ms, e0, e1)
         if (r > 0) tf(r) = ms
      end do
      c2 = cd
      print "(a20,a,f8.1,a,f8.1,a,a)", label, "   geam ", rate(m, n, median(tg)), &
         "   fast ", rate(m, n, median(tf)), " GB/s   ", &
         merge("bitwise equal    ", "MISMATCH         ", all(c1(1:n, 1:m) == c2(1:n, 1:m)))
      st = cublasDestroy(h)
      deallocate (ad, cd)
   end subroutine run

   function rate(m, n, t_ms) result(r)
      !! GB/s for reading and writing m x n doubles in t_ms milliseconds.
      integer, intent(in) :: m, n
      real(sp), intent(in) :: t_ms
      real(sp) :: r
      r = real(2.0_dp*8.0_dp*real(m, dp)*real(n, dp)/(real(t_ms, dp)*1.0e6_dp), sp)
   end function rate

   function median(t) result(r)
      !! The median of t, by insertion sort (t is short).
      real(sp), intent(in) :: t(:)
      real(sp) :: r
      real(sp) :: s(size(t))
      real(sp) :: x
      integer :: i, j
      s = t
      do i = 2, size(s)
         x = s(i)
         j = i - 1
         do while (j >= 1)
            if (s(j) <= x) exit
            s(j + 1) = s(j)
            j = j - 1
         end do
         s(j + 1) = x
      end do
      r = s((size(s) + 1)/2)
   end function median

end program fast_transpose_check
