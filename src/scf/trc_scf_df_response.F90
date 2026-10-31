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
! THE ORBITAL HESSIAN IN THE MO BASIS: df_rotation_t, what the fitted
! CPHF solve applies. [A' U]_ai = 2 [C_v^T G[D_U] C_o]_ai for
! D_U = C_v U C_o^T + C_o U^T C_v^T is, with the fitted MO blocks
! B^Q_ab, B^Q_ai and B^Q_ij built once (one pass over the unpacked B),
!
!     2 sum_Q B^Q_ai t_Q - sum_Q (B^Q_vv U B^Q_oo + B^Q_vo U^T B^Q_vo)_ai,
!     t_Q = 2 sum_bj B^Q_bj U_bj
!
! two GEMVs and four GEMMs per application, 2 nocc naux (nvir^2 + 3 nocc
! nvir) flops, no AO quantity and no unpacking at all: against the AO
! route's 2 n^2 nocc naux pass over B (n the Cartesian AO count, 1.3 nvir
! at cc-pVDZ) this is about half the work and all of it GEMM-shaped. It
! needs B^Q_vv, nvir^2 naux doubles (9.9 GB at gly10/cc-pVDZ with
! def2-universal-jkfit), so setup measures the free device memory and
! falls back to df_response_t's apply_rotation when the MO blocks do not
! fit.
!
module trc_scf_df_response
   use pic_types, only: dp, int64
   use trc_scf_df, only: df_jk_t
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use pic_logger, only: global_logger
   implicit none
   private
   public :: df_response_t, df_rotation_t

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

   type :: df_rotation_t
      !! The fitted RHF orbital Hessian's two-electron part, [A' U] =
      !! 2 C_v^T G[D_U] C_o, for a CPHF solve: in the MO basis when its
      !! blocks fit, else through df_response_t (this module's header).
      integer :: n = 0, nocc = 0, nvir = 0, naux = 0
      logical :: ready = .false.
      logical :: mo = .false.
      !! The MO-basis path is set up (else `ao` is).
      logical :: verbose_free = .true.
      !! Log the memory decision at setup.
      real(dp), allocatable :: eps(:)
      !! (nocc + nvir), host.
      real(dp), allocatable :: cvir(:, :), cocc(:, :)
      !! (n, nvir) and (n, nocc), device-resident (the AO path's C_v U and
      !! C_o^T G C_o).
      real(dp), allocatable :: bvo(:, :, :)
      !! (nvir, nocc, naux): B^Q_ai; device-resident.
      real(dp), allocatable :: bvv(:, :, :)
      !! (nvir, nvir, naux): B^Q_ab; device-resident.
      real(dp), allocatable :: boo(:, :, :)
      !! (nocc, naux, nocc): B^Q_ji, as (j, Q, i); device-resident.
      real(dp), allocatable :: wscr(:, :, :)
      !! (nvir, naux, nocc): U B_oo; device-resident.
      real(dp), allocatable :: mscr(:, :, :)
      !! (nocc, naux, nocc): U^T B_vo; device-resident.
      real(dp), allocatable :: tq(:), yao(:, :), gco(:, :)
      !! t, and the AO path's C_v U and G C_o; device-resident.
      type(df_response_t) :: ao
      !! The AO path, set up only when the MO blocks do not fit.
   contains
      procedure :: setup => df_rotation_setup
      procedure :: apply => df_rotation_apply
      procedure :: apply_occ => df_rotation_apply_occ
      procedure :: apply_general => df_rotation_apply_general
      procedure :: release => df_rotation_release
   end type df_rotation_t

contains

   subroutine df_rotation_setup(this, dfjk, lp, cmo, eps, nocc, force_ao)
      !! `cmo` (n, nmo) and `eps` (nmo): host arrays, the fit's Cartesian AO
      !! space, occupied columns first.
      class(df_rotation_t), intent(inout) :: this
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: cmo(:, :), eps(:)
      integer, intent(in) :: nocc
      logical, intent(in) :: force_ao
      !! Take the AO path whatever the memory: the fallback, for testing.
      real(dp), allocatable :: cmo_d(:, :), tens(:, :, :), half(:, :)
      type(trc_gemm_batch_t) :: gemm_shape
      integer(int64) :: need
      integer :: n, naux, npack, nv, no, ablk, p0, np

      call this%release()
      call dfjk%dims(n, naux, npack)
      no = nocc
      nv = size(cmo, 2) - nocc
      this%n = n
      this%nocc = no
      this%nvir = nv
      this%naux = naux
      allocate (this%eps(size(eps)), this%cvir(n, nv), this%cocc(n, no), this%tq(naux), this%yao(n, no), &
                this%gco(n, no))
      this%eps = eps
      this%cvir = cmo(:, no + 1:no + nv)
      this%cocc = cmo(:, 1:no)
      !$acc enter data copyin(this%cvir, this%cocc) create(this%tq, this%yao, this%gco)
      this%ready = .true.
      ! A small block: the build is GEMM-bound either way, and the MO blocks
      ! need the room.
      ablk = int(max(1_int64, min(int(naux, int64), (BLOCK_DOUBLES/4)/(int(n, int64)*int(n + nv, int64)))))
      ! B_vv, B_vo, B_oo, the two scratch tensors, one unpacked block and
      ! its half-transform.
      need = 8_int64*int(naux, int64)*(int(nv, int64)**2 + 2_int64*int(nv, int64)*int(no, int64) + &
                                       2_int64*int(no, int64)**2) + &
             8_int64*int(ablk, int64)*int(n, int64)*int(n + nv, int64)
      this%mo = need < cache_budget() .and. .not. force_ao
      if (this%verbose_free) then
         block
            character(len=160) :: logline
            write (logline, "(a,f8.2,a,f8.2,a,l1)") "  df rotation: MO blocks need ", real(need, dp)/2.0_dp**30, &
               " GiB, free ", real(cache_budget(), dp)/2.0_dp**30, " GiB: MO path ", this%mo
            call global_logger%performance(trim(logline))
         end block
      end if
      if (.not. this%mo) then
         call this%ao%setup(dfjk, lp, cmo, eps, nocc, keep_blocks=.true.)
         return
      end if
      allocate (this%bvo(nv, no, naux), this%bvv(nv, nv, naux), this%boo(no, naux, no))
      allocate (this%wscr(nv, naux, no), this%mscr(no, naux, no))
      allocate (cmo_d(n, no + nv), tens(n, n, ablk), half(nv, n*ablk))
      cmo_d = cmo(:, 1:no + nv)
      !$acc enter data create(this%bvo, this%bvv, this%boo, this%wscr, this%mscr)
      !$acc enter data copyin(cmo_d) create(tens, half)
      do p0 = 1, naux, ablk
         np = min(ablk, naux - p0 + 1)
         call unpack_columns(n, npack, naux, np, p0, dfjk%bpk, tens)
         ! half(a, nu q) = sum_mu C_v(mu, a) B(mu, nu, q); then per q
         ! B_vv(:, :, q) = half_q C_v and B_vo(:, :, q) = half_q C_o.
         call lp%gemm("T", "N", nv, n*np, n, 1.0_dp, this%cvir, n, tens, n, 0.0_dp, half, nv)
         gemm_shape = trc_gemm_batch(nv, nv, n, np, sb=0_int64)
         call lp%gemm_batched(gemm_shape, half, this%cvir, this%bvv(1, 1, p0))
         gemm_shape = trc_gemm_batch(nv, no, n, np, sb=0_int64)
         call lp%gemm_batched(gemm_shape, half, cmo_d, this%bvo(1, 1, p0))
         ! B_oo(j, q, i) = C_o^T B_q C_o: half_o reuses `half`'s storage.
         call lp%gemm("T", "N", no, n*np, n, 1.0_dp, cmo_d, n, tens, n, 0.0_dp, half, no)
         gemm_shape = trc_gemm_batch(no, no, n, np, ldc=no*naux, sb=0_int64, sc=int(no, int64))
         call lp%gemm_batched(gemm_shape, half, cmo_d, this%boo(1, p0, 1))
      end do
      !$acc exit data delete(cmo_d, tens, half)
      deallocate (cmo_d, tens, half)
   end subroutine df_rotation_setup

   subroutine df_rotation_apply(this, dfjk, lp, u, au)
      !! au := 2 C_v^T G[D_U] C_o for D_U = C_v U C_o^T + transpose, u and
      !! au (nvir, nocc) device-resident.
      class(df_rotation_t), intent(inout) :: this
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: u(:, :)
      real(dp), intent(inout) :: au(:, :)
      type(trc_gemm_batch_t) :: gemm_shape
      integer :: nv, no, naux, nov

      if (.not. this%ready) error stop "trc_scf_df_response: apply before setup"
      nv = this%nvir
      no = this%nocc
      naux = this%naux
      if (.not. this%mo) then
         call lp%gemm("N", "N", this%n, no, nv, 1.0_dp, this%cvir, this%n, u, nv, 0.0_dp, this%yao, this%n)
         call this%ao%apply_rotation(dfjk, lp, this%yao, this%gco)
         call lp%gemm("T", "N", nv, no, this%n, 2.0_dp, this%cvir, this%n, this%gco, this%n, 0.0_dp, au, nv)
         return
      end if
      nov = nv*no
      ! J: t = 2 B_vo^T u, au = 2 B_vo t
      call lp%gemv("T", nov, naux, 2.0_dp, this%bvo, nov, u, 0.0_dp, this%tq)
      call lp%gemv("N", nov, naux, 2.0_dp, this%bvo, nov, this%tq, 0.0_dp, au)
      ! K, first term: W(b, Q, i) = sum_j U(b, j) B_oo(j, Q, i);
      ! au -= sum_(b, Q) B_vv(a, b, Q) W(b, Q, i)
      call lp%gemm("N", "N", nv, naux*no, no, 1.0_dp, u, nv, this%boo, no, 0.0_dp, this%wscr, nv)
      call lp%gemm("N", "N", nv, no, nv*naux, -1.0_dp, this%bvv, nv, this%wscr, nv*naux, 1.0_dp, au, nv)
      ! K, second term: M(j, Q, i) = sum_b U(b, j) B_vo(b, i, Q);
      ! au -= sum_(j, Q) B_vo(a, j, Q) M(j, Q, i)
      gemm_shape = trc_gemm_batch(no, no, nv, naux, ldc=no*naux, sa=0_int64, sc=int(no, int64), transa="T")
      call lp%gemm_batched(gemm_shape, u, this%bvo, this%mscr)
      call lp%gemm("N", "N", nv, no, no*naux, -1.0_dp, this%bvo, nv, this%mscr, no*naux, 1.0_dp, au, nv)
   end subroutine df_rotation_apply

   subroutine df_rotation_apply_occ(this, dfjk, lp, u, goo)
      !! goo := C_o^T G[D_U] C_o (nocc, nocc), the occupied-occupied block of
      !! the same potential `apply` projects on C_v^T: with the MO blocks,
      !!
      !!   sum_Q B^Q_ki t_Q - (1/2) sum_Q (B^Q_ov U B^Q_oo + B^Q_oo U^T B^Q_vo)_ki
      !!
      !! u (nvir, nocc) and goo device-resident. Clobbers the scratch.
      class(df_rotation_t), intent(inout) :: this
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: u(:, :)
      real(dp), intent(inout) :: goo(:, :)
      type(trc_gemm_batch_t) :: gemm_shape
      integer :: nv, no, naux, nov

      if (.not. this%ready) error stop "trc_scf_df_response: apply before setup"
      nv = this%nvir
      no = this%nocc
      naux = this%naux
      if (.not. this%mo) then
         call lp%gemm("N", "N", this%n, no, nv, 1.0_dp, this%cvir, this%n, u, nv, 0.0_dp, this%yao, this%n)
         call this%ao%apply_rotation(dfjk, lp, this%yao, this%gco)
         call lp%gemm("T", "N", no, no, this%n, 1.0_dp, this%cocc, this%n, this%gco, this%n, 0.0_dp, goo, no)
         return
      end if
      nov = nv*no
      ! J: goo(:, i) = sum_Q B_oo(:, Q, i) t_Q, one (nocc x naux) GEMV per i
      call lp%gemv("T", nov, naux, 2.0_dp, this%bvo, nov, u, 0.0_dp, this%tq)
      gemm_shape = trc_gemm_batch(no, 1, naux, no, sb=0_int64, sc=int(no, int64))
      call lp%gemm_batched(gemm_shape, this%boo, this%tq, goo)
      ! K, (B_oo U^T B_vo): M(j, Q, i) = sum_b U(b, j) B_vo(b, i, Q);
      ! goo -= 1/2 sum_(j, Q) B_oo(j, Q, k) M(j, Q, i)
      gemm_shape = trc_gemm_batch(no, no, nv, naux, ldc=no*naux, sa=0_int64, sc=int(no, int64), transa="T")
      call lp%gemm_batched(gemm_shape, u, this%bvo, this%mscr)
      call lp%gemm("T", "N", no, no, no*naux, -0.5_dp, this%boo, no*naux, this%mscr, no*naux, 1.0_dp, goo, no)
      ! K, (B_ov U B_oo): W(b, Q, i) = sum_j U(b, j) B_oo(j, Q, i), then per Q
      ! B_vo(:, :, Q)^T W(:, Q, :) into mscr's room (nocc, nocc, naux), summed.
      call lp%gemm("N", "N", nv, naux*no, no, 1.0_dp, u, nv, this%boo, no, 0.0_dp, this%wscr, nv)
      gemm_shape = trc_gemm_batch(no, no, nv, naux, ldb=nv*naux, sb=int(nv, int64), sc=int(no, int64)**2, &
                                  transa="T")
      call lp%gemm_batched(gemm_shape, this%bvo, this%wscr, this%mscr)
      call add_slices(no, naux, -0.5_dp, this%mscr, goo)
   end subroutine df_rotation_apply_occ

   subroutine df_rotation_apply_general(this, dfjk, lp, xmo, xao, gvo, goo)
      !! gvo := C_v^T G[X] C_o and goo := C_o^T G[X] C_o for a general
      !! symmetric X: `xmo` (nmo, nmo) in the MO basis, `xao` the same in the
      !! AO basis (read only by the AO fallback). With the MO blocks,
      !! t_Q = sum_rs B^Q_rs x_rs, V^Q_r i = sum_s x_rs B^Q_si for r in o and
      !! in v, and K_pi = sum_Q sum_r B^Q_pr V^Q_ri -- ten GEMMs and GEMVs,
      !! the largest 2 nvir^2 nocc naux. Every array device-resident.
      class(df_rotation_t), intent(inout) :: this
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: xmo(:, :), xao(:, :)
      real(dp), intent(inout) :: gvo(:, :), goo(:, :)
      real(dp), allocatable :: xvo(:, :), xvv(:, :), xoo(:, :)
      type(trc_gemm_batch_t) :: gemm_shape
      integer :: nv, no, naux, nov, nmo

      if (.not. this%ready) error stop "trc_scf_df_response: apply before setup"
      nv = this%nvir
      no = this%nocc
      naux = this%naux
      if (.not. this%mo) then
         call this%ao%apply(dfjk, lp, xao, this%gco)
         call lp%gemm("T", "N", nv, no, this%n, 1.0_dp, this%cvir, this%n, this%gco, this%n, 0.0_dp, gvo, nv)
         call lp%gemm("T", "N", no, no, this%n, 1.0_dp, this%cocc, this%n, this%gco, this%n, 0.0_dp, goo, no)
         return
      end if
      nov = nv*no
      nmo = no + nv
      allocate (xvo(nv, no), xvv(nv, nv), xoo(no, no))
      !$acc enter data create(xvo, xvv, xoo)
      call copy_block(nmo, no + 1, 1, nv, no, xmo, xvo)
      call copy_block(nmo, no + 1, no + 1, nv, nv, xmo, xvv)
      call copy_block(nmo, 1, 1, no, no, xmo, xoo)
      ! J: t = 2 B_vo^T x_vo + B_vv^T x_vv + B_oo^T x_oo
      call lp%gemv("T", nov, naux, 2.0_dp, this%bvo, nov, xvo, 0.0_dp, this%tq)
      call lp%gemv("T", nv*nv, naux, 1.0_dp, this%bvv, nv*nv, xvv, 1.0_dp, this%tq)
      call add_oo_trace(no, naux, this%boo, xoo, this%tq)
      call lp%gemv("N", nov, naux, 1.0_dp, this%bvo, nov, this%tq, 0.0_dp, gvo)
      gemm_shape = trc_gemm_batch(no, 1, naux, no, sb=0_int64, sc=int(no, int64))
      call lp%gemm_batched(gemm_shape, this%boo, this%tq, goo)
      ! K, occupied r: V_o(j, Q, i) = x_oo B_oo + x_ov B_vo, then
      ! gvo -= B_vo V_o / 2 and goo -= B_oo^T V_o / 2
      call lp%gemm("N", "N", no, naux*no, no, 1.0_dp, xoo, no, this%boo, no, 0.0_dp, this%mscr, no)
      gemm_shape = trc_gemm_batch(no, no, nv, naux, ldc=no*naux, sa=0_int64, sc=int(no, int64), transa="T")
      call lp%gemm_batched(gemm_shape, xvo, this%bvo, this%mscr, beta=1.0_dp)
      call lp%gemm("N", "N", nv, no, no*naux, -0.5_dp, this%bvo, nv, this%mscr, no*naux, 1.0_dp, gvo, nv)
      call lp%gemm("T", "N", no, no, no*naux, -0.5_dp, this%boo, no*naux, this%mscr, no*naux, 1.0_dp, goo, no)
      ! K, virtual r: V_v(a, Q, i) = x_vo B_oo + x_vv B_vo, then
      ! gvo -= B_vv V_v / 2 and goo -= sum_Q B_vo(:, :, Q)^T V_v(:, Q, :) / 2
      call lp%gemm("N", "N", nv, naux*no, no, 1.0_dp, xvo, nv, this%boo, no, 0.0_dp, this%wscr, nv)
      gemm_shape = trc_gemm_batch(nv, no, nv, naux, ldc=nv*naux, sa=0_int64, sc=int(nv, int64))
      call lp%gemm_batched(gemm_shape, xvv, this%bvo, this%wscr, beta=1.0_dp)
      call lp%gemm("N", "N", nv, no, nv*naux, -0.5_dp, this%bvv, nv, this%wscr, nv*naux, 1.0_dp, gvo, nv)
      gemm_shape = trc_gemm_batch(no, no, nv, naux, ldb=nv*naux, sb=int(nv, int64), sc=int(no, int64)**2, &
                                  transa="T")
      call lp%gemm_batched(gemm_shape, this%bvo, this%wscr, this%mscr)
      call add_slices(no, naux, -0.5_dp, this%mscr, goo)
      !$acc exit data delete(xvo, xvv, xoo)
      deallocate (xvo, xvv, xoo)
   end subroutine df_rotation_apply_general

   subroutine df_rotation_release(this)
      class(df_rotation_t), intent(inout) :: this
      if (this%ready) then
         !$acc exit data delete(this%cvir, this%cocc, this%tq, this%yao, this%gco)
         if (this%mo) then
            !$acc exit data delete(this%bvo, this%bvv, this%boo, this%wscr, this%mscr)
         end if
      end if
      call this%ao%release()
      if (allocated(this%eps)) deallocate (this%eps)
      if (allocated(this%cvir)) deallocate (this%cvir)
      if (allocated(this%cocc)) deallocate (this%cocc)
      if (allocated(this%tq)) deallocate (this%tq)
      if (allocated(this%yao)) deallocate (this%yao)
      if (allocated(this%gco)) deallocate (this%gco)
      if (allocated(this%bvo)) deallocate (this%bvo)
      if (allocated(this%bvv)) deallocate (this%bvv)
      if (allocated(this%boo)) deallocate (this%boo)
      if (allocated(this%wscr)) deallocate (this%wscr)
      if (allocated(this%mscr)) deallocate (this%mscr)
      this%ready = .false.
      this%mo = .false.
   end subroutine df_rotation_release

   subroutine df_response_setup(this, dfjk, lp, cmo, eps, nocc, keep_blocks)
      !! Keep `cmo` (n, nmo) and `eps` (nmo) -- host arrays, the Cartesian
      !! AO space of the fit, occupied columns first -- and build b = B C_o.
      !! `keep_blocks`: keep what fits of the unpacked B, for an iterative
      !! caller (a single application gains nothing from it).
      class(df_response_t), intent(inout) :: this
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: cmo(:, :), eps(:)
      integer, intent(in) :: nocc
      logical, intent(in) :: keep_blocks
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
      if (.not. keep_blocks) this%ncache = 0
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

   subroutine copy_block(ld, r0, c0, m, n, x, y)
      !! y := x(r0 : r0 + m - 1, c0 : c0 + n - 1), x with leading dimension ld.
      integer, intent(in) :: ld, r0, c0, m, n
      real(dp), intent(in) :: x(ld, *)
      real(dp), intent(out) :: y(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         y(i, j) = x(r0 - 1 + i, c0 - 1 + j)
      end do
   end subroutine copy_block

   subroutine add_oo_trace(no, naux, boo, xoo, tq)
      !! tq(Q) += sum_jk B_oo(j, Q, k) x_oo(j, k).
      integer, intent(in) :: no, naux
      real(dp), intent(in) :: boo(no, naux, no), xoo(no, no)
      real(dp), intent(inout) :: tq(naux)
      integer :: q, j, k
      real(dp) :: acc
      do concurrent(q=1:naux) local(acc, j, k)
         acc = 0.0_dp
         do k = 1, no
            do j = 1, no
               acc = acc + boo(j, q, k)*xoo(j, k)
            end do
         end do
         tq(q) = tq(q) + acc
      end do
   end subroutine add_oo_trace

   subroutine add_slices(no, naux, alpha, slices, goo)
      !! goo += alpha sum_Q slices(:, :, Q), slices (no, no, naux) flat.
      integer, intent(in) :: no, naux
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: slices(no, no, naux)
      real(dp), intent(inout) :: goo(no, no)
      integer :: k, i, q
      real(dp) :: acc
      do concurrent(i=1:no, k=1:no) local(acc, q)
         acc = 0.0_dp
         do q = 1, naux
            acc = acc + slices(k, i, q)
         end do
         goo(k, i) = goo(k, i) + alpha*acc
      end do
   end subroutine add_slices

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
