!
! Pulay's DIIS on the commutator (cdiis_t), and ADIIS (Hu & Yang, JCP 132,
! 054109) on the convex hull of the Fock/density history while damping has
! ended but the commutator is still too large for CDIIS to trust (adiis_t):
! see trc_scf_driver's own header for the two owning most of the actual
! theory. This module owns their histories, pushes a new iterate into
! them, and extrapolates/combines a step, plus the plain subroutines both
! share: the DIIS linear solve, the level shift and the damped step while
! neither extrapolation is active, and the one sequential history-shift
! loop both histories use.
!
! EXPLICIT-SHAPE DUMMIES. Every one of these takes `fock`/`dold`/`errv`/
! history arrays as plain arguments, never touching `st%anything` or a
! type's own OTHER component inside a `do concurrent` body -- see
! trc_scf_fock's own header for why that matters on the GPU build.
! `cdiis_t`/`adiis_t` get their OWN enter/exit data pair rather than
! folding into trc_scf_state's single one: the two histories are built
! only when DIIS/ADIIS actually run (ADIIS not at all unless `opts%adiis`),
! so trc_scf_run maps them separately, exactly as the original code did.
!
module trc_scf_diis
   use pic_types, only: dp
   use trc_linalg, only: trc_linalg_t
   implicit none
   private

   public :: cdiis_t, adiis_t, solve_diis, level_shift, damp

   type :: cdiis_t
      real(dp), allocatable :: fstore(:, :, :, :), estore(:, :, :, :)
      real(dp), allocatable :: bmat(:, :), rhs(:)
      integer :: n_used = 0
      integer :: ndiis = 0
   contains
      procedure :: alloc => cdiis_alloc
      procedure :: to_device => cdiis_to_device
      procedure :: release_device => cdiis_release_device
      procedure :: push => cdiis_push
      procedure :: extrapolate => cdiis_extrapolate
   end type cdiis_t

   type :: adiis_t
      real(dp), allocatable :: fstore_a(:, :, :, :), dstore_a(:, :, :, :)
      real(dp), allocatable :: cad(:)
      integer :: n_used = 0
      integer :: ndiis = 0
   contains
      procedure :: alloc => adiis_alloc
      procedure :: to_device => adiis_to_device
      procedure :: release_device => adiis_release_device
      procedure :: push => adiis_push
      procedure :: extrapolate => adiis_extrapolate
   end type adiis_t

contains

   subroutine cdiis_alloc(this, nao, nspin, ndiis)
      class(cdiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin, ndiis
      this%ndiis = ndiis
      this%n_used = 0
      allocate (this%fstore(nao, nao, nspin, ndiis), this%estore(nao, nao, nspin, ndiis))
      allocate (this%bmat(ndiis + 1, ndiis + 1), this%rhs(ndiis + 1))
      this%fstore = 0.0_dp
      this%estore = 0.0_dp
   end subroutine cdiis_alloc

   subroutine cdiis_to_device(this)
      class(cdiis_t), intent(inout) :: this
      !$acc enter data copyin(this%fstore, this%estore)
   end subroutine cdiis_to_device

   subroutine cdiis_release_device(this)
      class(cdiis_t), intent(inout) :: this
      !$acc exit data delete(this%fstore, this%estore)
   end subroutine cdiis_release_device

   subroutine cdiis_push(this, nao, nspin, fock, errv)
      !! Slot n_used+1 <- (fock, errv), dropping the oldest entry first
      !! once the history is full.
      class(cdiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: fock(nao, nao, nspin), errv(nao, nao, nspin)
      if (this%n_used < this%ndiis) then
         this%n_used = this%n_used + 1
      else
         call shift_history(nao, nspin, this%ndiis, this%fstore, this%estore)
      end if
      call diis_store(nao, nspin, this%n_used, fock, errv, this%fstore, this%estore)
   end subroutine cdiis_push

   subroutine cdiis_extrapolate(this, lp, nao, nspin, fock, ok)
      !!
      !! An ill-conditioned subspace loses its oldest entry, not the step.
      !! A failed solve followed by a plain Roothaan step is not neutral: on
      !! water/def2-SVP PBE0 the history holds errors of order 1 next to the
      !! newest at 1e-7, B's scaled pivots drop under solve_diis's 1e-14, and
      !! plain steps grow the commutator x1.28 per iteration along one
      !! direction, which keeps B singular: 65 iterations of runaway. Dropping
      !! the oldest entries until the system solves keeps the extrapolation;
      !! two distinct error vectors always solve, and that step cancels one
      !! growing mode.
      !!
      class(cdiis_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nspin
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      logical, intent(out) :: ok
      integer :: m, i, j
      do
         m = this%n_used + 1
         this%bmat = -1.0_dp
         this%bmat(m, m) = 0.0_dp
         do j = 1, this%n_used
            do i = 1, j
               this%bmat(i, j) = lp%dot_3d(this%estore(:, :, :, i), this%estore(:, :, :, j))
               this%bmat(j, i) = this%bmat(i, j)
            end do
         end do
         call solve_diis(m, this%bmat(1:m, 1:m), this%rhs(1:m), ok)
         if (ok .or. this%n_used <= 2) exit
         call shift_history(nao, nspin, this%n_used, this%fstore, this%estore)
         this%n_used = this%n_used - 1
      end do
      if (ok) call diis_combine(nao, nspin, this%n_used, this%rhs(1:this%n_used), this%fstore, fock)
   end subroutine cdiis_extrapolate

   subroutine adiis_alloc(this, nao, nspin, ndiis)
      class(adiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin, ndiis
      this%ndiis = ndiis
      this%n_used = 0
      allocate (this%fstore_a(nao, nao, nspin, ndiis), this%dstore_a(nao, nao, nspin, ndiis))
      allocate (this%cad(ndiis))
      this%fstore_a = 0.0_dp
      this%dstore_a = 0.0_dp
   end subroutine adiis_alloc

   subroutine adiis_to_device(this)
      class(adiis_t), intent(inout) :: this
      !$acc enter data copyin(this%fstore_a, this%dstore_a)
   end subroutine adiis_to_device

   subroutine adiis_release_device(this)
      class(adiis_t), intent(inout) :: this
      !$acc exit data delete(this%fstore_a, this%dstore_a)
   end subroutine adiis_release_device

   subroutine adiis_push(this, nao, nspin, fock, dold)
      !! Slot n_used+1 <- (fock, dold): the Fock/density pair ADIIS's convex
      !! combination runs over, kept apart from cdiis_t's own history (see
      !! trc_scf_options_t%adiis's own comment on why).
      class(adiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: fock(nao, nao, nspin), dold(nao, nao, nspin)
      if (this%n_used < this%ndiis) then
         this%n_used = this%n_used + 1
      else
         call shift_history(nao, nspin, this%ndiis, this%fstore_a, this%dstore_a)
      end if
      call diis_store(nao, nspin, this%n_used, fock, dold, this%fstore_a, this%dstore_a)
   end subroutine adiis_push

   subroutine adiis_extrapolate(this, lp, nao, nspin, fock)
      !! One entry: the simplex is a point, c = 1 on it, a no-op. Skip the
      !! solve and leave `fock` as it is, which is also what `adiis_coeffs`
      !! would return from "start at the newest entry" after finding no
      !! improving step.
      class(adiis_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nspin
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      if (this%n_used > 1) then
         call adiis_coeffs(lp, this%n_used, this%fstore_a, this%dstore_a, this%cad(1:this%n_used))
         call diis_combine(nao, nspin, this%n_used, this%cad(1:this%n_used), this%fstore_a, fock)
      end if
   end subroutine adiis_extrapolate

   subroutine diis_store(nao, nspin, k, fock, companion, fstore, cstore)
      !! One push, shared by CDIIS (companion = the commutator `errv`) and
      !! ADIIS (companion = the density `dold`): slot k <- the current
      !! iterate in both histories at once.
      integer, intent(in) :: nao, nspin, k
      real(dp), intent(in) :: fock(nao, nao, nspin), companion(nao, nao, nspin)
      real(dp), intent(inout) :: fstore(nao, nao, nspin, *), cstore(nao, nao, nspin, *)
      integer :: i, j, s
      do concurrent(i=1:nao, j=1:nao, s=1:nspin)
         fstore(i, j, s, k) = fock(i, j, s)
         cstore(i, j, s, k) = companion(i, j, s)
      end do
   end subroutine diis_store

   subroutine diis_combine(nao, nspin, m, coeff, fstore, fock)
      !! fock := sum_k coeff(k) * fstore(:,:,:,k), the extrapolated/combined
      !! step both CDIIS (coeff = the solved weights) and ADIIS (coeff = the
      !! simplex coefficients) take once they have more than one history slot.
      integer, intent(in) :: nao, nspin, m
      real(dp), intent(in) :: coeff(m)
      real(dp), intent(in) :: fstore(nao, nao, nspin, *)
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      integer :: i, j, s, k
      do concurrent(i=1:nao, j=1:nao, s=1:nspin)
         fock(i, j, s) = 0.0_dp
      end do
      do k = 1, m
         do concurrent(i=1:nao, j=1:nao, s=1:nspin)
            fock(i, j, s) = fock(i, j, s) + coeff(k)*fstore(i, j, s, k)
         end do
      end do
   end subroutine diis_combine

   subroutine level_shift(lp, nao, nspin, mu, occ, smat, dmat, w1, w2, fock)
      !! F += mu (S - S D S / occ) on the virtual space while DIIS is off,
      !! per spin: see trc_scf_options_t%level_shift's own comment.
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: mu, occ
      real(dp), intent(in) :: smat(nao, nao), dmat(nao, nao, nspin)
      real(dp), intent(inout) :: w1(nao, nao), w2(nao, nao)
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      integer :: i, j, s
      do s = 1, nspin
         call lp%gemm_mat(smat, dmat(:, :, s), w1)
         call lp%gemm_mat(w1, smat, w2)
         do concurrent(i=1:nao, j=1:nao)
            fock(i, j, s) = fock(i, j, s) + mu*(smat(i, j) - w2(i, j)/occ)
         end do
      end do
   end subroutine level_shift

   subroutine damp(nao, nspin, alpha, dold, dmat)
      !! The damped step while DIIS/A-DIIS are off: dmat := alpha*dmat +
      !! (1 - alpha)*dold, the weight `alpha` given to the new density.
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: dold(nao, nao, nspin)
      real(dp), intent(inout) :: dmat(nao, nao, nspin)
      integer :: i, j, s
      do concurrent(i=1:nao, j=1:nao, s=1:nspin)
         dmat(i, j, s) = alpha*dmat(i, j, s) + (1.0_dp - alpha)*dold(i, j, s)
      end do
   end subroutine damp

   subroutine shift_history(nao, nspin, n, a, bb)
      !! Slot k <- slot k+1 for k = 1..n-1 in two histories at once, which
      !! drops the oldest entry and frees slot n.
      integer, intent(in) :: nao, nspin, n
      real(dp), intent(inout) :: a(nao, nao, nspin, *), bb(nao, nao, nspin, *)
      integer :: i, j, s, k
      ! Slot by slot, in order: slot k reads slot k+1, so the slot loop
      ! carries a dependence. As a `do concurrent` over k, a device slot
      ! was overwritten before being read and the commutator climbed for
      ! eight iterations on cholesterol; the host never saw it.
      do k = 1, n - 1
         do concurrent(i=1:nao, j=1:nao, s=1:nspin)
            a(i, j, s, k) = a(i, j, s, k + 1)
            bb(i, j, s, k) = bb(i, j, s, k + 1)
         end do
      end do
   end subroutine shift_history

   !
   ! ADIIS coefficients over the first `m` slots of fstore_a/dstore_a. Slot
   ! `m` is the newest, so Hu & Yang's "entry n" is the last index.
   ! Minimises over the simplex {c >= 0, sum c = 1}
   !
   !     f(c) = 2 sum_i c_i <D_i - D_m, F_m>
   !            + sum_ij c_i c_j <D_i - D_m, F_j - F_m>
   !
   ! by projected gradient with backtracking (a few hundred host flops once the
   ! device dot products are in hand). Starting with all weight on the newest
   ! entry means no progress degrades to the plain Roothaan step.
   !
   subroutine adiis_coeffs(lp, m, fstore_a, dstore_a, c)
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: m
      !! CONTIGUOUS: dot_3d below takes a trailing-index section of each --
      !! fstore_a(:,:,:,jj) -- which is only provably contiguous (no copy,
      !! no device/host mismatch) when the whole array it is sliced from
      !! is declared contiguous; see trc_linalg's own header.
      real(dp), intent(in), contiguous :: fstore_a(:, :, :, :), dstore_a(:, :, :, :)
      real(dp), intent(out) :: c(m)
      real(dp) :: f, ft, step
      real(dp) :: dvec(m), amat(m, m), g(m), trial(m)
      integer :: ii, jj, itr, ls
      logical :: improved
      do ii = 1, m
         dvec(ii) = lp%dot_3d(dstore_a(:, :, :, ii), fstore_a(:, :, :, m)) &
                    - lp%dot_3d(dstore_a(:, :, :, m), fstore_a(:, :, :, m))
         do jj = 1, m
            amat(ii, jj) = lp%dot_3d(dstore_a(:, :, :, ii), fstore_a(:, :, :, jj)) &
                           - lp%dot_3d(dstore_a(:, :, :, ii), fstore_a(:, :, :, m)) &
                           - lp%dot_3d(dstore_a(:, :, :, m), fstore_a(:, :, :, jj)) &
                           + lp%dot_3d(dstore_a(:, :, :, m), fstore_a(:, :, :, m))
         end do
      end do
      c = 0.0_dp
      c(m) = 1.0_dp
      f = adiis_value(m, dvec, amat, c)
      step = 1.0_dp
      do itr = 1, 200
         do ii = 1, m
            g(ii) = 2.0_dp*dvec(ii)
            do jj = 1, m
               g(ii) = g(ii) + (amat(ii, jj) + amat(jj, ii))*c(jj)
            end do
         end do
         improved = .false.
         do ls = 1, 30
            trial = c - step*g
            call project_simplex(m, trial)
            ft = adiis_value(m, dvec, amat, trial)
            if (ft < f - 1.0e-14_dp) then
               c = trial
               f = ft
               improved = .true.
               exit
            end if
            step = step*0.5_dp
         end do
         if (.not. improved) exit
         step = step*1.5_dp
      end do
   end subroutine adiis_coeffs

   pure function adiis_value(m, dvec, amat, c) result(f)
      integer, intent(in) :: m
      real(dp), intent(in) :: dvec(m), amat(m, m), c(m)
      real(dp) :: f
      integer :: ii, jj
      f = 0.0_dp
      do ii = 1, m
         f = f + 2.0_dp*c(ii)*dvec(ii)
         do jj = 1, m
            f = f + c(ii)*c(jj)*amat(ii, jj)
         end do
      end do
   end function adiis_value

   ! Projection onto the simplex {c : sum c = 1, c >= 0}: the standard
   ! sort-and-threshold algorithm, which is what enforces ADIIS's
   ! non-negativity constraint, the thing that keeps it from wandering
   ! off far from convergence.
   pure subroutine project_simplex(m, c)
      integer, intent(in) :: m
      real(dp), intent(inout) :: c(m)
      real(dp) :: css, theta, t
      real(dp) :: u(m)
      integer :: ii, rho
      u = c
      call sort_desc(m, u)
      css = 0.0_dp
      theta = 0.0_dp
      rho = 0
      do ii = 1, m
         css = css + u(ii)
         t = (css - 1.0_dp)/real(ii, dp)
         if (u(ii) - t > 0.0_dp) then
            rho = ii
            theta = t
         end if
      end do
      if (rho == 0) then
         c = 1.0_dp/real(m, dp)
         return
      end if
      c = max(c - theta, 0.0_dp)
   end subroutine project_simplex

   ! Insertion sort, descending: m <= ndiis, so this is never worth more
   ! than a handful of compares.
   pure subroutine sort_desc(m, u)
      integer, intent(in) :: m
      real(dp), intent(inout) :: u(m)
      integer :: ii, jj
      real(dp) :: t
      do ii = 1, m - 1
         do jj = ii + 1, m
            if (u(jj) > u(ii)) then
               t = u(ii)
               u(ii) = u(jj)
               u(jj) = t
            end if
         end do
      end do
   end subroutine sort_desc

   !
   ! The DIIS system B c = (0, ..., 0, -1), pivoted elimination on the host,
   ! with the error block scaled to O(1) against the -1 border first. The
   ! weights are invariant under that scaling; only the Lagrange multiplier
   ! moves, and nothing reads it. `ok` is false on a pivot below 1e-14, in
   ! which case the caller keeps the unextrapolated Fock matrix.
   !
   subroutine solve_diis(m, b, c, ok)
      integer, intent(in) :: m
      real(dp), intent(in) :: b(m, m)
      real(dp), intent(out) :: c(m)
      logical, intent(out) :: ok
      real(dp), allocatable :: a(:, :)
      real(dp) :: scale, pivot, factor
      integer :: i, j, p
      allocate (a(m, m + 1))
      a(:, 1:m) = b
      a(:, m + 1) = 0.0_dp
      a(m, m + 1) = -1.0_dp
      scale = maxval(abs(b(1:m - 1, 1:m - 1)))
      if (scale > 0.0_dp) a(1:m - 1, 1:m - 1) = a(1:m - 1, 1:m - 1)/scale
      ok = .true.
      c = 0.0_dp
      do i = 1, m
         p = i
         do j = i + 1, m
            if (abs(a(j, i)) > abs(a(p, i))) p = j
         end do
         if (p /= i) a([i, p], :) = a([p, i], :)
         pivot = a(i, i)
         if (abs(pivot) < 1.0e-14_dp) then
            ok = .false.
            return
         end if
         do j = i + 1, m
            factor = a(j, i)/pivot
            a(j, :) = a(j, :) - factor*a(i, :)
         end do
      end do
      do i = m, 1, -1
         c(i) = (a(i, m + 1) - sum(a(i, i + 1:m)*c(i + 1:m)))/a(i, i)
      end do
   end subroutine solve_diis

end module trc_scf_diis
