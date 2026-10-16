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
!   4. sandwich vs two raw gemms done by hand; density vs one raw gemm done
!      by hand; density with k = 0 against zero.
!   5. A conformance mismatch is NOT exercised here: gemm_mat's `error
!      stop` on a shape mismatch would simply end this program, not report
!      a FAIL, so there is nothing a test could usefully assert.
!
program check_linalg_mat
   use pic_types, only: dp
   use trc_linalg, only: trc_linalg_t
   implicit none

   integer :: nbad
   type(trc_linalg_t) :: lp

   nbad = 0
   call lp%init(20)

   call test_gemm_mat()
   call test_dot_mat()
   call test_syev_mat()
   call test_sandwich_density()

   call lp%release()

   print '(a)', ''
   if (nbad == 0) then
      print '(a)', '  RESULT: PASS'
   else
      print '(a,i0,a)', '  RESULT: FAIL (', nbad, ' checks differ)'
      stop 1
   end if

contains

   subroutine report(what, ok)
      character(len=*), intent(in) :: what
      logical, intent(in) :: ok
      character(len=4) :: verdict
      verdict = 'ok  '
      if (.not. ok) then
         verdict = 'FAIL'
         nbad = nbad + 1
      end if
      print '(a,a56,a,a)', '  ', what, '   ', verdict
   end subroutine report

   ! A closed-form fill: no two (offset, i, j) combinations collide, so
   ! A, B and C never accidentally share a value that would mask a bug.
   pure real(dp) function fill(offset, i, j) result(v)
      integer, intent(in) :: offset, i, j
      v = sin(real(3*i - 2*j + offset, dp)) + 0.07_dp*real(i, dp) - 0.03_dp*real(j, dp)
   end function fill

   ! =======================================================================
   ! 1. gemm_mat vs the raw gemm
   ! =======================================================================
   subroutine run_gemm_combo(ta, tb, explicit)
      character, intent(in) :: ta, tb
      logical, intent(in) :: explicit
      integer, parameter :: MM = 6, KK = 4, NN = 5
      real(dp), allocatable :: a_store(:, :), b_store(:, :), c_mat(:, :), c_raw(:, :)
      real(dp) :: alpha, beta
      integer :: i, j, lda, ldb
      logical :: tra, trb
      character(len=64) :: what

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
      call report(trim(what), all(c_mat == c_raw))
      deallocate (a_store, b_store, c_mat, c_raw)
   end subroutine run_gemm_combo

   subroutine test_gemm_mat()
      character, parameter :: OPS(2) = ["N", "T"]
      integer, parameter :: MM = 6, KK = 4
      real(dp), allocatable :: a3(:, :, :), b(:, :), c_mat(:, :), c_raw(:, :)
      integer :: ia, ib, iv, i, j

      do ia = 1, 2
         do ib = 1, 2
            do iv = 1, 2
               call run_gemm_combo(OPS(ia), OPS(ib), iv == 2)
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
      call report('gemm_mat vs raw, A as a 3D section a(:,:,2)', all(c_mat == c_raw))
      deallocate (a3, b, c_mat, c_raw)
   end subroutine test_gemm_mat

   ! =======================================================================
   ! 2. dot_2d / dot_3d vs the raw dot
   ! =======================================================================
   subroutine test_dot_mat()
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
      call report('dot_2d vs raw dot', r_mat == r_raw)
      deallocate (x2, y2)

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
      call report('dot_3d vs raw dot, sections estore(:,:,:,k)', r_mat == r_raw)
      deallocate (estore)
   end subroutine test_dot_mat

   ! =======================================================================
   ! 3. syev_mat vs the raw syev
   ! =======================================================================
   subroutine test_syev_mat()
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

      call report('syev_mat eigenvectors vs raw syev', all(a_mat == a_raw))
      call report('syev_mat eigenvalues vs raw syev', all(w_mat == w_raw))
      deallocate (a_mat, a_raw, w_mat, w_raw)
   end subroutine test_syev_mat

   ! =======================================================================
   ! 4. sandwich vs two raw gemms; density vs one raw gemm; density(k<=0)
   ! =======================================================================
   subroutine test_sandwich_density()
      integer, parameter :: N = 5, NC = 4, K = 3
      real(dp), allocatable :: x(:, :), f(:, :), out_s(:, :), work_s(:, :)
      real(dp), allocatable :: out_r(:, :), work_r(:, :)
      real(dp), allocatable :: c(:, :), d_mat(:, :), d_raw(:, :), d_zero(:, :)
      real(dp) :: occ
      integer :: i, j

      allocate (x(N, N), f(N, N), out_s(N, N), work_s(N, N), out_r(N, N), work_r(N, N))
      do j = 1, N
         do i = 1, N
            x(i, j) = fill(50, i, j)
            f(i, j) = fill(60, i, j)
         end do
      end do
      f = 0.5_dp*(f + transpose(f))
      out_s = 0.0_dp; work_s = 0.0_dp; out_r = 0.0_dp; work_r = 0.0_dp

      !$acc enter data copyin(x, f, out_s, work_s, out_r, work_r)
      call lp%sandwich(x, f, out_s, work_s)
      call lp%gemm("N", "N", N, N, N, 1.0_dp, f, N, x, N, 0.0_dp, work_r, N)
      call lp%gemm("T", "N", N, N, N, 1.0_dp, x, N, work_r, N, 0.0_dp, out_r, N)
      !$acc update self(out_s, work_s, out_r, work_r)
      !$acc exit data delete(x, f, out_s, work_s, out_r, work_r)

      call report('sandwich vs two raw gemms', all(out_s == out_r))
      deallocate (x, f, out_s, work_s, out_r, work_r)

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

      call report('density vs raw gemm', all(d_mat == d_raw))
      call report('density(k<=0) is zero', all(d_zero == 0.0_dp))
      deallocate (c, d_mat, d_raw, d_zero)
   end subroutine test_sandwich_density

end program check_linalg_mat
