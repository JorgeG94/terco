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
! AN ORBITAL ROTATION, X = Y C_o^T + C_o Y^T (a CPHF trial density, Y =
! C_v U): `apply_rotation` never forms X. X b^P = Y (C_o^T b^P) + C_o
! (Y^T b^P), and C_o^T b^P is fixed, so W is three (n or nocc) x nocc x
! (nocc naux) GEMMs -- 6 n nocc^2 naux flops instead of 2 n^2 nocc naux --
! and t_P = 2 sum_i (Y^T b^P)_ii needs no pass over the packed B at all.
!
! UNPACKED B. Every application reads all of B once, unpacked block by
! block; the blocks that fit in the device memory left after setup (half
! a gigabyte kept back) are unpacked once, at setup, and kept, so an
! iterative solve pays the unpack only for the rest.
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
      real(dp), allocatable :: boo(:, :, :)
      !! (nocc, naux, nocc): C_o^T b^P; device-resident.
      real(dp), allocatable :: m2(:, :, :)
      !! (nocc, naux, nocc): Y^T b^P, apply_rotation's scratch; device-resident.
      real(dp), allocatable :: cache(:, :, :)
      !! (n, n, ncache) unpacked B for the first ncache auxiliary functions
      !! (whole blocks); device-resident.
      integer :: ncache = 0
   contains
      procedure :: setup => df_response_setup
      procedure :: apply => df_response_apply
      procedure :: apply_rotation => df_response_apply_rotation
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
      allocate (this%xpk(npack), this%tq(naux), this%boo(nocc, naux, nocc), this%m2(nocc, naux, nocc))
      !$acc enter data copyin(this%cmo) create(this%bocc, this%work, this%tens, this%xpk, this%tq)
      !$acc enter data create(this%boo, this%m2)
      ! The unpacked blocks to keep, whole ones only, from what is free now.
      this%ncache = this%ablk*int(min(int((naux + this%ablk - 1)/this%ablk, int64), &
                                      max(0_int64, cache_budget()/(8_int64*int(n, int64)**2*this%ablk))))
      this%ncache = min(this%ncache, naux)
      allocate (this%cache(n, n, max(this%ncache, 1)))
      !$acc enter data create(this%cache)
      this%ready = .true.
      do p0 = 1, naux, this%ablk
         np = min(this%ablk, naux - p0 + 1)
         if (p0 <= this%ncache) then
            call unpack_columns(n, npack, naux, np, p0, dfjk%bpk, this%cache(1, 1, p0))
            ! b(:, P, i) = B^P C_o(:, i): items step one P (n elements) in C,
            ! whose columns i are n naux apart.
            gemm_shape = trc_gemm_batch(n, nocc, n, np, ldc=n*naux, sb=0_int64, sc=int(n, int64))
            call lp%gemm_batched(gemm_shape, this%cache(1, 1, p0), this%cmo, this%bocc(1, p0, 1))
         else
            call unpack_columns(n, npack, naux, np, p0, dfjk%bpk, this%tens)
            gemm_shape = trc_gemm_batch(n, nocc, n, np, ldc=n*naux, sb=0_int64, sc=int(n, int64))
            call lp%gemm_batched(gemm_shape, this%tens, this%cmo, this%bocc(1, p0, 1))
         end if
      end do
      ! boo(j, P, i) = sum_mu C(mu, j) b(mu, P, i)
      call lp%gemm("T", "N", nocc, naux*nocc, n, 1.0_dp, this%cmo, n, this%bocc, n, 0.0_dp, this%boo, nocc)
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
      call b_pass(this, dfjk, lp, gco)
   end subroutine df_response_apply

   subroutine df_response_apply_rotation(this, dfjk, lp, y, gco)
      !! gco := G[X] C_o for X = Y C_o^T + C_o Y^T, without X: see this
      !! module's header. y (n, nocc) and gco (n, nocc) device-resident.
      class(df_response_t), intent(inout) :: this
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: y(:, :)
      real(dp), intent(inout) :: gco(:, :)
      integer :: n, naux, nocc

      if (.not. this%ready) error stop "trc_scf_df_response: apply before setup"
      n = this%n
      naux = this%naux
      nocc = this%nocc
      ! m2(j, P, i) = Y_j . b(:, P, i);  t_P = 2 sum_i m2(i, P, i)
      call lp%gemm("T", "N", nocc, naux*nocc, n, 1.0_dp, y, n, this%bocc, n, 0.0_dp, this%m2, nocc)
      call trace_diagonal(nocc, naux, this%m2, this%tq)
      ! W = Y boo + C_o m2
      call lp%gemm("N", "N", n, naux*nocc, nocc, 1.0_dp, y, n, this%boo, nocc, 0.0_dp, this%work, n)
      call lp%gemm("N", "N", n, naux*nocc, nocc, 1.0_dp, this%cmo, n, this%m2, nocc, 1.0_dp, this%work, n)
      call combine_jk(n, naux, nocc, this%tq, this%cmo, this%work)
      call b_pass(this, dfjk, lp, gco)
   end subroutine df_response_apply_rotation

   subroutine b_pass(this, dfjk, lp, gco)
      !! gco := sum_P B^P W'^P, a block of P per GEMM: B as (mu, (nu, P)),
      !! W' as ((nu, P), i) with leading dimension n naux; the kept blocks
      !! read in place, the rest unpacked again. The output is only
      !! n x nocc against an inner dimension of n naux: cuBLAS splits it
      !! over 64x64 tiles (gly10/cc-pVDZ: 5.1 TFLOP/s on a V100, nocc = 155
      !! filling 81% of three 64-wide tiles); a strided batch of eight
      !! partial sums per block measured the same.
      type(df_response_t), intent(inout) :: this
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(inout) :: gco(:, :)
      integer :: n, naux, nocc, p0, np, ld
      n = this%n
      naux = this%naux
      nocc = this%nocc
      ld = n*naux
      do p0 = 1, naux, this%ablk
         np = min(this%ablk, naux - p0 + 1)
         if (p0 <= this%ncache) then
            call lp%gemm("N", "N", n, nocc, n*np, 1.0_dp, this%cache(1, 1, p0), n, this%work(1, p0, 1), ld, &
                         merge(0.0_dp, 1.0_dp, p0 == 1), gco, n)
         else
            call unpack_columns(n, this%npack, naux, np, p0, dfjk%bpk, this%tens)
            call lp%gemm("N", "N", n, nocc, n*np, 1.0_dp, this%tens, n, this%work(1, p0, 1), ld, &
                         merge(0.0_dp, 1.0_dp, p0 == 1), gco, n)
         end if
      end do
   end subroutine b_pass

   subroutine df_response_release(this)
      class(df_response_t), intent(inout) :: this
      if (this%ready) then
         !$acc exit data delete(this%cmo, this%bocc, this%work, this%tens, this%xpk, this%tq)
         !$acc exit data delete(this%boo, this%m2, this%cache)
      end if
      if (allocated(this%boo)) deallocate (this%boo)
      if (allocated(this%m2)) deallocate (this%m2)
      if (allocated(this%cache)) deallocate (this%cache)
      this%ncache = 0
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

   integer(int64) function cache_budget() result(bytes)
      !! Device bytes free for kept blocks of unpacked B, half a gigabyte
      !! held back for the libraries. The host build keeps up to 2 GB.
#ifdef TRC_DEVICE
      use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, &
                         acc_property_free_memory
      bytes = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                   acc_property_free_memory), int64) - 512_int64*1024*1024
#else
      bytes = 2_int64*1024*1024*1024
#endif
   end function cache_budget

   subroutine trace_diagonal(nocc, naux, m2, tq)
      !! tq(P) = 2 sum_i m2(i, P, i).
      integer, intent(in) :: nocc, naux
      real(dp), intent(in) :: m2(nocc, naux, nocc)
      real(dp), intent(out) :: tq(naux)
      integer :: p, i
      real(dp) :: acc
      do concurrent(p=1:naux) local(acc, i)
         acc = 0.0_dp
         do i = 1, nocc
            acc = acc + m2(i, p, i)
         end do
         tq(p) = 2.0_dp*acc
      end do
   end subroutine trace_diagonal

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
