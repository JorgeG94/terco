!
! trc_linalg's shape-inferred forms against the raw forms they must stay
! bitwise identical to.
!
! trc_linalg_t grows a shape-inferred form alongside each raw one. syev
! becomes a true `generic`, since every existing caller passes it a plain,
! exactly-ranked array. gemm and dot cannot: existing callers reach their
! raw forms through sequence association (trc_rimp2_driver's tens(1,1,q)
! stands for "the sub-matrix starting here"; trc_scf_driver's ADIIS feeds
! rank-3 estore(:,:,:,i) sections to dot's flat x(*)), and Fortran will
! not resolve a GENERIC reference that way -- only a reference to one
! already-identified specific procedure gets that trick. So gemm_mat and
! dot_2d/dot_3d are their own names, reached alongside -- not through --
! the unchanged lp%gemm(...) / lp%dot(...) raw calls. Either way, the new
! form and the raw form reach the SAME library entry point with the SAME
! arguments, so this check does not compare with a tolerance -- it
! compares with ==. Any difference at all is a bug in the new form, not
! rounding.
!
! Every array here is deterministic (a closed-form fill, no RNG), so the
! host and device builds see the same bit pattern and this check is
! meaningful on both.
!
!   1. gemm_mat vs the raw gemm, all four transa/transb combinations,
!      non-square shapes, default and explicit alpha/beta. One combination
!      passes A as a(:,:,2), a contiguous section of a rank-3 array --
!      this is what proves host_data on a contiguous assumed-shape dummy
!      that received a SECTION (not a whole array) still resolves to a
!      real device address on nvfortran, which is the failure mode the new
!      forms exist to avoid (see trc_linalg.F90's own doc comments).
!   2. dot_2d and dot_3d vs the raw dot, including a section
!      estore(:,:,:,k) of a rank-4 array.
!   3. syev_mat vs the raw syev.
!   4. density vs one raw gemm done by hand; density with k = 0 against
!      zero.
!   5. gemm_batched: trc_gemm_batch's tight defaults for each transpose,
!      and the batch against a plain triple loop, for all four transposes,
!      tight tiles and a shared, non-tight operand (the trc_sph pattern).
!      The fills are small integers, so every summation order is exact
!      and == holds for cuBLAS's batched kernels too. count = 0 must leave
!      C alone.
!   6. A conformance mismatch is NOT exercised here: gemm_mat's `error
!      stop` on a shape mismatch would simply end this program, not report
!      a FAIL, so there is nothing a test could usefully assert.
!
module check_linalg_mat_suite
   use pic_types, only: dp, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   implicit none
   private
   public :: collect_linalg_mat, lp

   type(trc_linalg_t) :: lp
      !! Shared across every test below: initialised before run_testsuite
      !! and released after (program check_linalg_mat), safe because this
      !! suite always runs with parallel=.false.

contains

   subroutine collect_linalg_mat(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("1. gemm_mat matches raw gemm: every transpose/alpha-beta combination "// &
                               "and a 3D array section", test_gemm_mat), &
                  new_unittest("2. dot_2d/dot_3d match raw dot, including a rank-4 section", test_dot_mat), &
                  new_unittest("3. syev_mat matches raw syev", test_syev_mat), &
                  new_unittest("3b. syev_mat reports a signalling-NaN matrix through info, not error stop", &
                               test_syev_failure), &
                  new_unittest("4. density matches a raw gemm; density(k<=0) is zero", test_density), &
                  new_unittest("5a. trc_gemm_batch defaults for each transpose", test_gemm_batch_defaults), &
                  new_unittest("5b. gemm_batched vs triple loop, ta/tb=NN", test_gemm_batched_nn), &
                  new_unittest("5b. gemm_batched vs triple loop, ta/tb=NT", test_gemm_batched_nt), &
                  new_unittest("5b. gemm_batched vs triple loop, ta/tb=TN", test_gemm_batched_tn), &
                  new_unittest("5b. gemm_batched vs triple loop, ta/tb=TT", test_gemm_batched_tt), &
                  new_unittest("5c. gemm_batched, shared non-tight B and offset tiles; count = 0 leaves C alone", &
                               test_gemm_batched_shared) &
                  ]
   end subroutine collect_linalg_mat

   ! A closed-form fill: no two (offset, i, j) combinations collide, so
   ! A, B and C never accidentally share a value that would mask a bug.
   pure real(dp) function fill(offset, i, j) result(v)
      integer, intent(in) :: offset, i, j
      v = sin(real(3*i - 2*j + offset, dp)) + 0.07_dp*real(i, dp) - 0.03_dp*real(j, dp)
   end function fill

   ! =======================================================================
   ! 1. gemm_mat vs the raw gemm
   ! =======================================================================
   subroutine test_gemm_mat(error)
      type(error_type), allocatable, intent(out) :: error
      character, parameter :: OPS(2) = ["N", "T"]
      integer, parameter :: MM = 6, KK = 4, NN = 5
      real(dp), allocatable :: a_store(:, :), b_store(:, :), c_mat(:, :), c_raw(:, :)
      real(dp), allocatable :: a3(:, :, :), b(:, :)
      real(dp) :: alpha, beta
      integer :: ia, ib, iv, i, j, lda, ldb
      logical :: tra, trb, explicit
      character(len=64) :: what
      character :: ta, tb

      do ia = 1, 2
         do ib = 1, 2
            do iv = 1, 2
               ta = OPS(ia); tb = OPS(ib); explicit = (iv == 2)
               tra = ta == "T"; trb = tb == "T"
               if (explicit) then
                  alpha = 1.7_dp; beta = 0.3_dp
               else
                  alpha = 1.0_dp; beta = 0.0_dp
               end if

               if (tra) then
                  allocate (a_store(KK, MM))
               else
                  allocate (a_store(MM, KK))
               end if
               do j = 1, size(a_store, 2)
                  do i = 1, size(a_store, 1)
                     a_store(i, j) = fill(0, i, j)
                  end do
               end do

               if (trb) then
                  allocate (b_store(NN, KK))
               else
                  allocate (b_store(KK, NN))
               end if
               do j = 1, size(b_store, 2)
                  do i = 1, size(b_store, 1)
                     b_store(i, j) = fill(100, i, j)
                  end do
               end do

               allocate (c_mat(MM, NN), c_raw(MM, NN))
               do j = 1, NN
                  do i = 1, MM
                     ! Beta = 0 never reads C, but a non-trivial starting value still
                     ! has to be fed to both forms identically for the explicit-beta
                     ! case to mean anything.
                     c_mat(i, j) = fill(200, i, j)
                  end do
               end do
               c_raw = c_mat

               lda = size(a_store, 1); ldb = size(b_store, 1)
               !$acc enter data copyin(a_store, b_store, c_mat, c_raw)
               call lp%gemm_mat(a_store, b_store, c_mat, transa=ta, transb=tb, alpha=alpha, beta=beta)
               call lp%gemm(ta, tb, MM, NN, KK, alpha, a_store, lda, b_store, ldb, beta, c_raw, MM)
               !$acc update self(c_mat, c_raw)
               !$acc exit data delete(a_store, b_store, c_mat, c_raw)

               write (what, '(a,a1,a1,a,l1)') 'gemm_mat vs raw, ta/tb=', ta, tb, ' explicit=', explicit
               call check(error, all(c_mat == c_raw), trim(what))
               deallocate (a_store, b_store, c_mat, c_raw)
               if (allocated(error)) return
            end do
         end do
      end do

      ! A as a(:,:,2): a contiguous section of a rank-3 array, exactly the
      ! shape trc_scf_driver passes (fock(:,:,s), estore(:,:,:,k)).
      allocate (a3(MM, KK, 2), b(KK, KK), c_mat(MM, KK), c_raw(MM, KK))
      do j = 1, KK
         do i = 1, MM
            a3(i, j, 1) = fill(300, i, j)
            a3(i, j, 2) = fill(400, i, j)
         end do
      end do
      do j = 1, KK
         do i = 1, KK
            b(i, j) = fill(500, i, j)
         end do
      end do
      c_mat = 0.0_dp; c_raw = 0.0_dp
      !$acc enter data copyin(a3, b, c_mat, c_raw)
      call lp%gemm_mat(a3(:, :, 2), b, c_mat)
      call lp%gemm("N", "N", MM, KK, KK, 1.0_dp, a3(:, :, 2), MM, b, KK, 0.0_dp, c_raw, MM)
      !$acc update self(c_mat, c_raw)
      !$acc exit data delete(a3, b, c_mat, c_raw)
      call check(error, all(c_mat == c_raw), 'gemm_mat vs raw, A as a 3D section a(:,:,2)')
      deallocate (a3, b, c_mat, c_raw)
   end subroutine test_gemm_mat

   ! =======================================================================
   ! 2. dot_2d / dot_3d vs the raw dot
   ! =======================================================================
   subroutine test_dot_mat(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: PP = 5, QQ = 3
      real(dp), allocatable :: x2(:, :), y2(:, :)
      real(dp), allocatable :: estore(:, :, :, :)
      real(dp) :: r_mat, r_raw
      integer :: i, j

      allocate (x2(PP, QQ), y2(PP, QQ))
      do j = 1, QQ
         do i = 1, PP
            x2(i, j) = fill(10, i, j)
            y2(i, j) = fill(20, i, j)
         end do
      end do
      !$acc enter data copyin(x2, y2)
      r_mat = lp%dot_2d(x2, y2)
      r_raw = lp%dot(size(x2), x2, y2)
      !$acc exit data delete(x2, y2)
      call check(error, r_mat == r_raw, 'dot_2d vs raw dot')
      deallocate (x2, y2)
      if (allocated(error)) return

      ! A rank-4 array's section estore(:,:,:,k) fed to dot_3d, exactly as
      ! trc_scf_driver's ADIIS feeds estore(:,:,:,i) to the raw dot today.
      allocate (estore(PP, QQ, 2, 4))
      block
         integer :: a, b, c, d
         do d = 1, 4
            do c = 1, 2
               do b = 1, QQ
                  do a = 1, PP
                     estore(a, b, c, d) = fill(1000 + 7*d, a, b + 10*c)
                  end do
               end do
            end do
         end do
      end block
      !$acc enter data copyin(estore)
      r_mat = lp%dot_3d(estore(:, :, :, 2), estore(:, :, :, 3))
      r_raw = lp%dot(size(estore(:, :, :, 2)), estore(:, :, :, 2), estore(:, :, :, 3))
      !$acc exit data delete(estore)
      call check(error, r_mat == r_raw, 'dot_3d vs raw dot, sections estore(:,:,:,k)')
      deallocate (estore)
   end subroutine test_dot_mat

   ! =======================================================================
   ! 3. syev_mat vs the raw syev
   ! =======================================================================
   subroutine test_syev_mat(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 6
      real(dp), allocatable :: a_mat(:, :), a_raw(:, :), w_mat(:), w_raw(:)
      integer :: i, j

      allocate (a_mat(N, N), a_raw(N, N), w_mat(N), w_raw(N))
      do j = 1, N
         do i = 1, N
            a_mat(i, j) = fill(40, i, j) + merge(5.0_dp*real(i, dp), 0.0_dp, i == j)
         end do
      end do
      a_mat = 0.5_dp*(a_mat + transpose(a_mat))
      a_raw = a_mat

      !$acc enter data copyin(a_mat, a_raw) create(w_mat, w_raw)
      call lp%syev(a_mat, w_mat)
      call lp%syev(N, a_raw, N, w_raw)
      !$acc update self(a_mat, a_raw, w_mat, w_raw)
      !$acc exit data delete(a_mat, a_raw, w_mat, w_raw)

      call check(error, all(a_mat == a_raw), 'syev_mat eigenvectors vs raw syev')
      if (allocated(error)) then
         deallocate (a_mat, a_raw, w_mat, w_raw)
         return
      end if
      call check(error, all(w_mat == w_raw), 'syev_mat eigenvalues vs raw syev')
      deallocate (a_mat, a_raw, w_mat, w_raw)
   end subroutine test_syev_mat

   ! =======================================================================
   ! 3b. A diagonalisation `pic_syev`/cuSOLVER genuinely cannot do must
   ! come back through the OPTIONAL `info` argument, not `error stop`.
   ! dsyev returns a nonzero info ("failed to converge") for a NaN
   ! matrix whether it started quiet or signalling -- any ordinary
   ! arithmetic on a signalling NaN quiets it (the standard IEEE 754
   ! trap-then-return-qNaN behaviour), so by the time a NaN seeded
   ! earlier reaches an eigensolver it is already quiet either way;
   ! `ieee_signaling_nan` here is just the intrinsic this test uses to
   ! seed one deliberately, not a distinction dsyev's own result
   ! depends on. A density read past the end of an undersized buffer
   ! (an SCF guess, or app/trc_dlfind_bridge.F90's own per-image warm
   ! start) is exactly how a NaN like this reaches the diagonalisation
   ! one level up, in trc_scf_driver.F90; this is the cheap,
   ! deterministic version of the same thing, at the trc_linalg level.
   ! =======================================================================
   subroutine test_syev_failure(error)
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 6
      real(dp), allocatable :: a(:, :), w(:)
      integer :: info, i

      allocate (a(N, N), w(N))
      a = 0.0_dp
      do i = 1, N
         a(i, i) = 1.0_dp
      end do
      a(2, 5) = ieee_value(a(2, 5), ieee_signaling_nan)
      a(5, 2) = a(2, 5)

      !$acc enter data copyin(a) create(w)
      call lp%syev(a, w, info)
      !$acc exit data delete(a, w)

      call check(error, info /= 0, &
                 'syev_mat reports a signalling-NaN matrix through info, not error stop')
      deallocate (a, w)
   end subroutine test_syev_failure

   ! =======================================================================
   ! 4. density vs one raw gemm; density(k<=0)
   ! =======================================================================
   subroutine test_density(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 5, NC = 4, K = 3
      real(dp), allocatable :: c(:, :), d_mat(:, :), d_raw(:, :), d_zero(:, :)
      real(dp) :: occ
      integer :: i, j

      ! density: D = occ * C(:,1:k) C(:,1:k)^T, C wider than K columns so
      ! C(:,1:k) is a genuine section, not the whole array.
      occ = 2.0_dp
      allocate (c(N, NC), d_mat(N, N), d_raw(N, N), d_zero(N, N))
      do j = 1, NC
         do i = 1, N
            c(i, j) = fill(70, i, j)
         end do
      end do
      d_mat = 0.0_dp; d_raw = 0.0_dp
      d_zero = 7.0_dp ! a nonzero pattern the k<=0 path must stamp out

      !$acc enter data copyin(c, d_mat, d_raw, d_zero)
      call lp%density(c, K, occ, d_mat)
      call lp%gemm("N", "T", N, N, K, occ, c(:, 1:K), N, c(:, 1:K), N, 0.0_dp, d_raw, N)
      call lp%density(c, 0, occ, d_zero)
      !$acc update self(d_mat, d_raw, d_zero)
      !$acc exit data delete(c, d_mat, d_raw, d_zero)

      call check(error, all(d_mat == d_raw), 'density vs raw gemm')
      if (allocated(error)) then
         deallocate (c, d_mat, d_raw, d_zero)
         return
      end if
      call check(error, all(d_zero == 0.0_dp), 'density(k<=0) is zero')
      deallocate (c, d_mat, d_raw, d_zero)
   end subroutine test_density

   ! =======================================================================
   ! 5. gemm_batched
   ! =======================================================================

   ! A small integer, so products and sums are exact in any order.
   pure real(dp) function ifill(offset, i) result(v)
      integer, intent(in) :: offset, i
      v = real(mod(7*i + offset, 11) - 5, dp)
   end function ifill

   subroutine test_gemm_batch_defaults(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_gemm_batch_t) :: g
      ! op(A) 6 x 4, op(B) 4 x 5, C 6 x 5.
      g = trc_gemm_batch(6, 5, 4, 3)
      call check(error, g%lda == 6 .and. g%ldb == 4 .and. g%ldc == 6 .and. &
                 g%sa == 24 .and. g%sb == 20 .and. g%sc == 30 .and. g%count == 3, &
                 'trc_gemm_batch NN defaults')
      if (allocated(error)) return
      g = trc_gemm_batch(6, 5, 4, 3, transa="T")
      call check(error, g%lda == 4 .and. g%ldb == 4 .and. g%ldc == 6 .and. &
                 g%sa == 24 .and. g%sb == 20 .and. g%sc == 30, &
                 'trc_gemm_batch TN defaults')
      if (allocated(error)) return
      g = trc_gemm_batch(6, 5, 4, 3, transb="T")
      call check(error, g%lda == 6 .and. g%ldb == 5 .and. g%ldc == 6 .and. &
                 g%sa == 24 .and. g%sb == 20 .and. g%sc == 30, &
                 'trc_gemm_batch NT defaults')
      if (allocated(error)) return
      ! An explicit leading dimension sets the default stride with it.
      g = trc_gemm_batch(6, 5, 4, 3, ldb=9, sa=0_int64, transa="T", transb="T")
      call check(error, g%lda == 4 .and. g%ldb == 9 .and. &
                 g%ldc == 6 .and. g%sa == 0 .and. g%sb == 36 .and. g%sc == 30, &
                 'trc_gemm_batch TT with ldb, sa overrides')
   end subroutine test_gemm_batch_defaults

   subroutine test_gemm_batched_nn(error)
      type(error_type), allocatable, intent(out) :: error
      call run_gemm_batched(error, "N", "N")
   end subroutine test_gemm_batched_nn

   subroutine test_gemm_batched_nt(error)
      type(error_type), allocatable, intent(out) :: error
      call run_gemm_batched(error, "N", "T")
   end subroutine test_gemm_batched_nt

   subroutine test_gemm_batched_tn(error)
      type(error_type), allocatable, intent(out) :: error
      call run_gemm_batched(error, "T", "N")
   end subroutine test_gemm_batched_tn

   subroutine test_gemm_batched_tt(error)
      type(error_type), allocatable, intent(out) :: error
      call run_gemm_batched(error, "T", "T")
   end subroutine test_gemm_batched_tt

   ! C_i = alpha op(A_i) op(B_i) + beta C_i over flat, tightly packed tiles,
   ! against a triple loop on the host copies.
   subroutine run_gemm_batched(error, ta, tb)
      type(error_type), allocatable, intent(out) :: error
      character, intent(in) :: ta, tb
      integer, parameter :: MM = 6, KK = 4, NN = 5, CNT = 3
      real(dp), parameter :: ALPHA = 2.0_dp, BETA = 0.5_dp
      real(dp), allocatable :: a(:), b(:), c(:), c_ref(:), c0(:)
      type(trc_gemm_batch_t) :: g
      integer :: i, j, l, it
      integer(int64) :: off_a, off_b, off_c
      real(dp) :: av, bv, acc
      character(len=64) :: what

      g = trc_gemm_batch(MM, NN, KK, CNT, transa=ta, transb=tb)
      allocate (a(MM*KK*CNT), b(KK*NN*CNT), c(MM*NN*CNT), c_ref(MM*NN*CNT))
      do i = 1, size(a)
         a(i) = ifill(1, i)
      end do
      do i = 1, size(b)
         b(i) = ifill(4, i)
      end do
      do i = 1, size(c)
         c(i) = ifill(8, i)
      end do
      c0 = c
      !$acc enter data copyin(a, b, c)
      call lp%gemm_batched(g, a, b, c, alpha=ALPHA, beta=BETA)
      !$acc update self(c)
      !$acc exit data delete(a, b, c)

      do it = 0, CNT - 1
         off_a = it*g%sa
         off_b = it*g%sb
         off_c = it*g%sc
         do j = 1, NN
            do i = 1, MM
               acc = 0.0_dp
               do l = 1, KK
                  if (ta == "T") then
                     av = a(off_a + l + (i - 1)*g%lda)
                  else
                     av = a(off_a + i + (l - 1)*g%lda)
                  end if
                  if (tb == "T") then
                     bv = b(off_b + j + (l - 1)*g%ldb)
                  else
                     bv = b(off_b + l + (j - 1)*g%ldb)
                  end if
                  acc = acc + av*bv
               end do
               c_ref(off_c + i + (j - 1)*g%ldc) = ALPHA*acc + BETA*c0(off_c + i + (j - 1)*g%ldc)
            end do
         end do
      end do
      write (what, '(a,a1,a1)') 'gemm_batched vs triple loop, ta/tb=', ta, tb
      call check(error, all(c == c_ref), trim(what))
      deallocate (a, b, c, c_ref, c0)
   end subroutine run_gemm_batched

   ! The trc_sph column fold: B is one block of a taller table shared by
   ! every item (sb = 0, ldb > k), the A and C tiles start inside larger
   ! arrays; then count = 0 must not touch C.
   subroutine test_gemm_batched_shared(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NROW = 5, NCI = 3, NSI = 2, LDT = 7, CNT = 4, LEAD = 3
      real(dp), allocatable :: m(:), t(:, :), s(:), s_ref(:), s_keep(:)
      type(trc_gemm_batch_t) :: g
      integer :: i, j, l, it
      integer(int64) :: oa, oc
      real(dp) :: acc

      allocate (m(LEAD + NROW*NCI*CNT), t(LDT, NSI), s(LEAD + NROW*NSI*CNT))
      do i = 1, size(m)
         m(i) = ifill(2, i)
      end do
      do j = 1, NSI
         do i = 1, LDT
            t(i, j) = ifill(3, i + 13*j)
         end do
      end do
      s = -99.0_dp
      s_ref = s
      g = trc_gemm_batch(NROW, NSI, NCI, CNT, ldb=LDT, sb=0_int64)
      !$acc enter data copyin(m, t, s)
      call lp%gemm_batched(g, m(LEAD + 1), t, s(LEAD + 1))
      !$acc update self(s)
      s_keep = s
      g = trc_gemm_batch(NROW, NSI, NCI, 0, ldb=LDT, sb=0_int64)
      call lp%gemm_batched(g, m(LEAD + 1), t, s(LEAD + 1), beta=2.0_dp)
      !$acc update self(s)
      !$acc exit data delete(m, t, s)

      do it = 0, CNT - 1
         oa = LEAD + it*int(NROW*NCI, int64)
         oc = LEAD + it*int(NROW*NSI, int64)
         do j = 1, NSI
            do i = 1, NROW
               acc = 0.0_dp
               do l = 1, NCI
                  acc = acc + m(oa + i + (l - 1)*NROW)*t(l, j)
               end do
               s_ref(oc + i + (j - 1)*NROW) = acc
            end do
         end do
      end do
      call check(error, all(s_keep == s_ref), 'gemm_batched, shared non-tight B, offset tiles')
      if (allocated(error)) then
         deallocate (m, t, s, s_ref, s_keep)
         return
      end if
      call check(error, all(s == s_keep), 'gemm_batched, count = 0 leaves C alone')
      deallocate (m, t, s, s_ref, s_keep)
   end subroutine test_gemm_batched_shared

end module check_linalg_mat_suite

program check_linalg_mat
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_linalg_mat_suite, only: collect_linalg_mat, lp
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   call lp%init(20)
   testsuites = [new_testsuite("check_linalg_mat", collect_linalg_mat)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   call lp%release()

   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_linalg_mat
