!
! The density-fitted reference potential projected on the occupied
! orbitals, G[X] C_occ, for any symmetric X -- what a response over a --df
! RHF reference needs, and all it needs.
!
! A coupled-perturbed solve, a Z-vector's right-hand side and the
! occupied-occupied potential of a relaxed density never use the full
! n x n G[X] = J[X] - K[X]/2: they read only C_v^T G C_o or C_o^T G C_o.
! With B = T L^-T (df_jk_t's own fitted tensor) and b^P = B^P C_o, built
! once,
!
!     K[X] C_o = sum_P B^P X b^P            J[X] C_o = sum_P b^P (B^P . X)
!
! so one application is the GEMM W = X b (every P at once), an elementwise
! W' = -W/2 + C_o t^T (t_P = B^P . X, the J part riding in the same
! contraction), and G[X] C_o = sum_P B^P W'^P, one GEMM per block of
! unpacked B -- 4 n^2 nocc naux flops, against the 4 n^2 rank(X) naux plus
! an n^3 eigenproblem df_jk_t%build spends on a general density (its K
! needs a factor X = Y Y^T, and a trial density's rank is 2 nocc, a relaxed
! density's n). gly10/cc-pVDZ, def2-universal-jkfit, V100: see
! trc_rimp2_grad.F90's own measurements.
!
! b and W are stored (n, naux, nocc), the occupied index LAST, so W = X b
! is one GEMM whose output already has the (nu, P) rows the second
! contraction reads contiguously for any block of P. Memory: two
! (n nocc naux) tensors plus one block of unpacked B. Exact: no screening,
! no eigenvalue cut.
!
module trc_scf_df_response
   use pic_types, only: dp, int64
   use trc_scf_df, only: df_jk_t
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   implicit none
   private
   public :: df_response_t

   integer(int64), parameter :: BLOCK_DOUBLES = 268435456_int64
   !! Unpacked B per block, at most 2 GB.

   type :: df_response_t
      !! b^P = B^P C_o for one fit and one set of orbitals, and the
      !! orbitals themselves (device-resident), for `apply`.
      integer :: n = 0, nocc = 0, nvir = 0, naux = 0, npack = 0, ablk = 0
      logical :: ready = .false.
      real(dp), allocatable :: cmo(:, :)
      !! (n, nocc + nvir): occupied columns first; device-resident.
      real(dp), allocatable :: eps(:)
      !! (nocc + nvir), host.
      real(dp), allocatable :: bocc(:, :, :)
      !! (n, naux, nocc): b^P_mu i = sum_nu B^P_mu nu C_nu i; device-resident.
      real(dp), allocatable :: work(:, :, :)
      !! (n, naux, nocc): W, then W'; device-resident.
      real(dp), allocatable :: tens(:, :, :)
      !! (n, n, ablk): one block of unpacked B; device-resident.
      real(dp), allocatable :: xpk(:), tq(:)
      !! X packed (off-diagonal doubled) and t = B^T xpk; device-resident.
   contains
      procedure :: setup => df_response_setup
      procedure :: apply => df_response_apply
      procedure :: release => df_response_release
   end type df_response_t

contains

   subroutine df_response_setup(this, dfjk, lp, cmo, eps, nocc)
      !! Keep `cmo` (n, nmo) and `eps` (nmo) -- host arrays, the Cartesian
      !! AO space of the fit, occupied columns first -- and build b = B C_o.
      class(df_response_t), intent(inout) :: this
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: cmo(:, :), eps(:)
      integer, intent(in) :: nocc
      type(trc_gemm_batch_t) :: gemm_shape
      integer :: n, naux, npack, p0, np

      call this%release()
      call dfjk%dims(n, naux, npack)
      if (size(cmo, 1) /= n) error stop "trc_scf_df_response: the orbitals are not in the fit's AO space"
      this%n = n
      this%naux = naux
      this%npack = npack
      this%nocc = nocc
      this%nvir = size(cmo, 2) - nocc
      this%ablk = int(max(1_int64, min(int(naux, int64), BLOCK_DOUBLES/(int(n, int64)*int(n, int64)))))
      allocate (this%cmo(n, size(cmo, 2)), this%eps(size(eps)))
      this%cmo = cmo
      this%eps = eps
      allocate (this%bocc(n, naux, nocc), this%work(n, naux, nocc), this%tens(n, n, this%ablk))
      allocate (this%xpk(npack), this%tq(naux))
      !$acc enter data copyin(this%cmo) create(this%bocc, this%work, this%tens, this%xpk, this%tq)
      do p0 = 1, naux, this%ablk
         np = min(this%ablk, naux - p0 + 1)
         call unpack_columns(n, npack, naux, np, p0, dfjk%bpk, this%tens)
         ! b(:, P, i) = B^P C_o(:, i): items step one P (n elements) in C,
         ! whose columns i are n naux apart.
         gemm_shape = trc_gemm_batch(n, nocc, n, np, ldc=n*naux, sb=0_int64, sc=int(n, int64))
         call lp%gemm_batched(gemm_shape, this%tens, this%cmo, this%bocc(1, p0, 1))
      end do
      this%ready = .true.
   end subroutine df_response_setup

   subroutine df_response_apply(this, dfjk, lp, x, gco)
      !! gco := G[X] C_o = (J[X] - K[X]/2) C_o, X (n, n) symmetric and
      !! device-resident, gco (n, nocc) device-resident.
      class(df_response_t), intent(inout) :: this
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(inout) :: gco(:, :)
      integer :: n, naux, nocc, p0, np, ld

      if (.not. this%ready) error stop "trc_scf_df_response: apply before setup"
      n = this%n
      naux = this%naux
      nocc = this%nocc
      ld = n*naux
      ! t = B^T x (the J part's fitted density), W = X b for every P.
      call pack_doubled(n, this%npack, x, this%xpk)
      call lp%gemv("T", this%npack, naux, 1.0_dp, dfjk%bpk, this%npack, this%xpk, 0.0_dp, this%tq)
      call lp%gemm("N", "N", n, naux*nocc, n, 1.0_dp, x, n, this%bocc, n, 0.0_dp, this%work, n)
      call combine_jk(n, naux, nocc, this%tq, this%cmo, this%work)
      ! G C_o = sum_P B^P W'^P, a block of P per GEMM: B as (mu, (nu, P)),
      ! W' as ((nu, P), i) with leading dimension n naux.
      do p0 = 1, naux, this%ablk
         np = min(this%ablk, naux - p0 + 1)
         call unpack_columns(n, this%npack, naux, np, p0, dfjk%bpk, this%tens)
         call lp%gemm("N", "N", n, nocc, n*np, 1.0_dp, this%tens, n, this%work(1, p0, 1), ld, &
                      merge(0.0_dp, 1.0_dp, p0 == 1), gco, n)
      end do
   end subroutine df_response_apply

   subroutine df_response_release(this)
      class(df_response_t), intent(inout) :: this
      if (this%ready) then
         !$acc exit data delete(this%cmo, this%bocc, this%work, this%tens, this%xpk, this%tq)
      end if
      if (allocated(this%cmo)) deallocate (this%cmo)
      if (allocated(this%eps)) deallocate (this%eps)
      if (allocated(this%bocc)) deallocate (this%bocc)
      if (allocated(this%work)) deallocate (this%work)
      if (allocated(this%tens)) deallocate (this%tens)
      if (allocated(this%xpk)) deallocate (this%xpk)
      if (allocated(this%tq)) deallocate (this%tq)
      this%ready = .false.
   end subroutine df_response_release

   ! ----------------------------------------------------------------------
   ! do concurrent kernels: explicit-shape dummies, device-resident arrays.
   ! ----------------------------------------------------------------------

   subroutine unpack_columns(n, npack, naux, np, p0, bpk, tens)
      !! tens(:, :, q) = B's column p0 - 1 + q, unpacked (both triangles).
      integer, intent(in) :: n, npack, naux, np, p0
      real(dp), intent(in) :: bpk(npack, naux)
      real(dp), intent(out) :: tens(n, n, np)
      integer(int64) :: t, n8
      integer :: mu, nu, q
      n8 = int(n, int64)
      do concurrent(t=0:n8*n8*int(np, int64) - 1) local(mu, nu, q)
         mu = int(mod(t, n8)) + 1
         nu = int(mod(t/n8, n8)) + 1
         q = int(t/(n8*n8)) + 1
         tens(mu, nu, q) = bpk(pack_index(max(mu, nu), min(mu, nu)), p0 - 1 + q)
      end do
   end subroutine unpack_columns

   subroutine pack_doubled(n, npack, x, xpk)
      !! xpk = the lower triangle of x, off-diagonal doubled (B^T xpk = t).
      integer, intent(in) :: n, npack
      real(dp), intent(in) :: x(n, n)
      real(dp), intent(out) :: xpk(npack)
      integer(int64) :: t, n8
      integer :: mu, nu
      n8 = int(n, int64)
      do concurrent(t=0:n8*n8 - 1) local(mu, nu)
         mu = int(mod(t, n8)) + 1
         nu = int(t/n8) + 1
         if (nu <= mu) xpk(pack_index(mu, nu)) = merge(1.0_dp, 2.0_dp, mu == nu)*x(mu, nu)
      end do
   end subroutine pack_doubled

   subroutine combine_jk(n, naux, nocc, tq, cmo, work)
      !! work(nu, P, i) := C(nu, i) t_P - work(nu, P, i)/2.
      integer, intent(in) :: n, naux, nocc
      real(dp), intent(in) :: tq(naux), cmo(n, nocc)
      real(dp), intent(inout) :: work(n, naux, nocc)
      integer(int64) :: t, n8, na8
      integer :: nu, p, i
      n8 = int(n, int64)
      na8 = int(naux, int64)
      do concurrent(t=0:n8*na8*int(nocc, int64) - 1) local(nu, p, i)
         nu = int(mod(t, n8)) + 1
         p = int(mod(t/n8, na8)) + 1
         i = int(t/(n8*na8)) + 1
         work(nu, p, i) = cmo(nu, i)*tq(p) - 0.5_dp*work(nu, p, i)
      end do
   end subroutine combine_jk

   pure integer function pack_index(mu, nu) result(p)
      !! trc_scf_df's packed index, mu >= nu.
      !$acc routine seq
      integer, intent(in) :: mu, nu
      p = mu*(mu - 1)/2 + nu
   end function pack_index

end module trc_scf_df_response
