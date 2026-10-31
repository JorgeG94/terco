!
! The correlation side of the RI-MP2 nuclear gradient, on the device.
!
! Stocks, Palethorpe and Barca, J. Chem. Theory Comput. 2024, 20, 2505,
! section 3, as metalquicha implements it (backends/cenzontle/correlation/
! mqc_czt_ri_mp2_gradient.f90, whose numpy prototype tools/cpu_validation/
! ri_mp2_gradient.py settled every index convention used here). The SCF
! side -- the Z-vector, the relaxed and energy-weighted densities, and the
! one-electron, overlap and four-centre derivative terms -- is
! src/scf/trc_rimp2_grad.F90's; this module only builds what needs the
! fitted integrals:
!
!   T_ij^ab = (ia|jb) / (e_i + e_j - e_a - e_b),  (ia|jb) = sum_Q B_ia^Q B_jb^Q
!   X_ij^ab = 2 T_ij^ab - T_ij^ba
!   doo_ij  = - sum_kab T_ki^ab X_kj^ab          dvv_ba = sum_ijc T_ij^ca X_ij^cb
!   Gamma^P_ia = sum_jbQ X_ij^ab B_jb^Q (L^-1)_QP        (eq 8, Cholesky form)
!   gamma_PQ   = sum_ia Gamma^P_ia c_ia^Q,  c = B L^-1    (eq 10)
!   L'_iq  = 2 sum_aP Gamma^P_ia (q a|P)                  (eq 13)
!   L''_aq = 2 sum_iP Gamma^P_ia (i q|P)                  (eq 14)
!
! with B = (ia|P) L^-T and M = (P|Q) = L L^T, exactly the B trc_rimp2_run
! fits. Only the ACTIVE occupied orbitals (nfrozen + 1 .. nocc) enter the
! amplitudes, Gamma and doo; q runs over every molecular orbital, which is
! where the occupied-frozen rows of the Lagrangian come from. And the
! first two terms of eq 6, the only ones the fit changes:
!
!   dE/dR += 4 sum_(mu nu P) Gamma^P_mu nu d(mu nu|P)/dR - 2 sum_PQ gamma_PQ d(P|Q)/dR
!
! Gamma^P_mu nu = C_act Gamma^P C_vir^T, through trc_df_grad_kernels'
! three- and two-centre derivative kernels (the --df gradient's own, which
! take an arbitrary three-index Y and an arbitrary symmetric Z: no
! generator change).
!
! THE ORDER OF WORK
!   1. the metric and its Cholesky factor; B, one auxiliary block of
!      (mu nu|P) at a time (trc_rimp2_run's own single-rank steps 1-4).
!   2. the amplitudes, a group of occupied i at a time: K = B_i B^T for
!      every j (not only j >= i: Gamma needs each row), then one kernel for
!      T and X, then Gamma'_i = X_i B (a strided batch), dvv (one GEMM over
!      the whole group) and doo (one GEMM per i). E_corr = sum Gamma' . B
!      comes for free and is returned, so a caller can check the gradient
!      differentiates the energy it thinks it does.
!   3. c = B L^-1 and Gamma = Gamma' L^-1 (two triangular solves), gamma =
!      Gamma^T c, and the two-centre derivative with Z = -(gamma + gamma^T).
!   4. (mu nu|P) again, block by block, for the Lagrangian and the three-
!      centre derivative together: G^P = C_vir Gamma^P, M^P = 2 C_act G^P^T
!      (the derivative kernels' one-sided Y, Y = M + M^T), R_o += (mu nu|P)
!      G^P, and U += Gamma^P H^P^T with H^P = (mu nu|P) C_act. Then
!      L' = 2 R_o^T C and L'' = 2 U C.
!
! A spherical auxiliary basis: B, Gamma and gamma live in its spherical
! space, as the energy's do; Gamma is expanded onto the Cartesian
! auxiliary functions once before step 4 (trc_sph's sph_cols_sph_to_cart)
! and gamma's Z before the two-centre kernels, so the kernels and the raw
! (mu nu|P) of step 4 stay Cartesian. The orbital basis is always the
! Cartesian AO space here: a spherical caller passes C_cart = T C_sph.
!
! MEMORY. B and Gamma are nocc nvir naux doubles each, plus three
! group-sized (nvir, nvir, nocc, ng) amplitude buffers; the group size is
! chosen from the device's free memory (amplitude_group). One rank.
!
module trc_rimp2_relax
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_df_2c, trc_df_3c
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_df_grad_kernels, only: df_grad_2c_driver, df_grad_3c_driver, df_grad_3c_plan, DFG_LMAX, &
                                  DFG_LMAX_AUX, DFG_NPCLS
   use trc_rys_select, only: rys_select_tables_init
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_sph, only: sph_metric_cart_to_sph, sph_aux3c_to_sph, sph_cols_sph_to_cart
   use trc_step_table, only: step_table_t
   use trc_rimp2_driver, only: trc_rimp2_options_t
   implicit none
   private

   public :: trc_rimp2_relax_t

   integer, parameter :: ST_METRIC = 1, ST_3C = 2, ST_TRANS = 3, ST_TRSM = 4, ST_KGEMM = 5, ST_AMP = 6, &
                         ST_GAMMA = 7, ST_DENSITY = 8, ST_FIT = 9, ST_D2C = 10, ST_LAG = 11, ST_D3C = 12, &
                         NSTEP = 12
   character(len=22), parameter :: STEP_NAME(NSTEP) = [character(len=22) :: &
                                   "metric + Cholesky", "(mu nu|P)", "transform (ia|P)", "B = X L^-T", &
                                   "K_i GEMMs", "amplitudes T, X", "Gamma' = X B (gemm)", &
                                   "doo, dvv (gemm)", "Gamma, gamma (trsm)", "2c derivative kernels", &
                                   "Lagrangian (gemm)", "3c derivative kernels"]

   integer, parameter :: NG_BEST = 8
   !! Occupied orbitals per amplitude group when memory allows: the K GEMM
   !! is (ng nvir) x (nocc nvir) x naux, and past a few hundred rows cuBLAS
   !! gains nothing from more.

   type :: trc_rimp2_relax_t
      !! What the correlation contributes to the RI-MP2 gradient beyond the
      !! one-particle quantities the SCF side owns: the unrelaxed density
      !! blocks, the Lagrangian, and the fitted two-particle derivative.
      !! Sized by `init`, filled by `build`; host arrays.
      integer :: nmo = 0       !! molecular orbitals (columns of the C passed to build)
      integer :: nocc = 0      !! doubly occupied orbitals, frozen ones included
      integer :: nfrozen = 0   !! of which frozen, from build's options
      integer :: nact = 0      !! nocc - nfrozen
      integer :: nvir = 0      !! nmo - nocc
      real(dp) :: e_corr = 0.0_dp
      !! The RI-MP2 correlation energy these amplitudes give.
      real(dp), allocatable :: doo(:, :)
      !! (nact, nact): eq 17, active occupied block, before symmetrising.
      real(dp), allocatable :: dvv(:, :)
      !! (nvir, nvir): eq 18, before symmetrising.
      real(dp), allocatable :: lag_o(:, :)
      !! (nact, nmo): L', the occupied Lagrangian, its own index first.
      real(dp), allocatable :: lag_v(:, :)
      !! (nvir, nmo): L'', the virtual Lagrangian.
      real(dp), allocatable :: grad(:, :)
      !! (3, natm): the three- and two-centre derivative terms of eq 6.
      character(len=256) :: message = ""
      !! Empty on success; the refusal otherwise.
   contains
      procedure :: init => rimp2_relax_init
      procedure :: build => rimp2_relax_build
   end type trc_rimp2_relax_t

   type :: relax_work_t
      !! One build's dimensions, timer table and device-resident arrays,
      !! passed whole to the private steps below. Each array component is
      !! mapped on its own (`!$acc enter data create(w%xr)`), never the
      !! object, and reaches a kernel or a library call only as an actual
      !! argument -- trc_scf_df's `df_jk_t%bpk` pattern.
      integer :: n = 0, nmo = 0, no = 0, nv = 0, nov = 0, naux = 0, naux_c = 0, nblk = 0, natm = 0
      integer(int64) :: nov8 = 0
      type(step_table_t) :: tbl
      real(dp), allocatable :: cmo(:, :), cocc(:, :), cvir(:, :), eo(:), ev(:), jm(:, :)
      real(dp), allocatable :: xr(:, :)
      !! B, then c = B L^-1: (nocc nvir, naux), rows a + (i-1) nvir.
      real(dp), allocatable :: gam(:, :)
      !! Gamma', then Gamma: (nocc nvir, naux).
      real(dp), allocatable :: gamc(:, :)
      !! Gamma on the Cartesian auxiliary functions (a spherical aux only).
      real(dp), allocatable :: doo(:, :), dvv(:, :), lag_o(:, :), lag_v(:, :), grad(:, :)
   end type relax_work_t

contains

   subroutine rimp2_relax_init(this, nmo, nocc)
      !! Set the orbital counts `build` reads its dummies' shapes from.
      class(trc_rimp2_relax_t), intent(inout) :: this
      integer, intent(in) :: nmo, nocc
      this%nmo = nmo
      this%nocc = nocc
      this%nvir = nmo - nocc
      this%message = ""
   end subroutine rimp2_relax_init

   recursive subroutine rimp2_relax_build(this, b, aux, pl, cmo, eps, opts, la)
      !! Everything in this module's header, at the orbitals `cmo`/`eps`
      !! (occupied first, ascending energy). `b`, `aux` and `pl` must be
      !! device-resident; `cmo`/`eps` are host arrays. `opts%nfrozen` core
      !! orbitals are left out of the correlation; `opts%aux_block` caps the
      !! auxiliary functions per (mu nu|P) block; `opts%comm` is not read
      !! (one rank). `la` absent: made and released here.
      class(trc_rimp2_relax_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: cmo(b%nao, this%nmo), eps(this%nmo)
      type(trc_rimp2_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(inout), optional :: la

      type(relax_work_t) :: w
      type(timer_type) :: tall
      integer :: nfz
      logical :: dense(NSTEP)

      if (.not. present(la)) then
         block
            type(trc_linalg_t) :: la_local
            call rimp2_relax_build(this, b, aux, pl, cmo, eps, opts, la_local)
            call la_local%release()
         end block
         return
      end if

      nfz = max(opts%nfrozen, 0)
      w%n = b%nao
      w%nmo = this%nmo
      w%no = this%nocc - nfz
      w%nv = this%nmo - this%nocc
      w%nov = w%no*w%nv
      w%nov8 = int(w%no, int64)*int(w%nv, int64)
      w%naux_c = aux%nao
      w%naux = merge(aux%nao_sph, aux%nao, aux%spherical)
      w%natm = b%natm
      w%nblk = w%naux_c
      if (opts%aux_block > 0) w%nblk = opts%aux_block
      this%nfrozen = nfz
      this%nact = w%no
      this%nvir = w%nv
      this%message = ""
      if (w%no < 1 .or. w%nv < 1 .or. w%nmo > w%n) then
         this%message = "trc_rimp2_relax: nothing to correlate (no active occupied or no virtual orbital)"
         return
      end if
      if (maxval(b%sh_l) > min(DF_LMAX, DFG_LMAX) .or. maxval(aux%sh_l) > min(DF_LMAX_AUX, DFG_LMAX_AUX)) then
         this%message = "trc_rimp2_relax: a basis is above the density-fitting (derivative) kernels' "// &
                        "ceiling (gen_df.py, gen_df_grad.py)"
         return
      end if
      if (aux%natm /= b%natm) then
         this%message = "trc_rimp2_relax: the auxiliary basis is on a different molecule"
         return
      end if

      call la%init(w%n)
      call rys_select_tables_init()
      call w%tbl%init("rimp2 grad", STEP_NAME)
      call la%sync()
      call tall%start()

      allocate (w%cmo(w%n, w%nmo), w%cocc(w%n, w%no), w%cvir(w%n, w%nv), w%eo(w%no), w%ev(w%nv))
      w%cmo = cmo
      w%cocc = cmo(:, nfz + 1:this%nocc)
      w%cvir = cmo(:, this%nocc + 1:w%nmo)
      w%eo = eps(nfz + 1:this%nocc)
      w%ev = eps(this%nocc + 1:w%nmo)
      allocate (w%jm(w%naux, w%naux), w%xr(w%nov, w%naux), w%gam(w%nov, w%naux))
      allocate (w%doo(w%no, w%no), w%dvv(w%nv, w%nv), w%lag_o(w%no, w%nmo), w%lag_v(w%nv, w%nmo), &
                w%grad(3, w%natm))
      !$acc enter data copyin(w%cmo, w%cocc, w%cvir, w%eo, w%ev)
      !$acc enter data create(w%jm, w%xr, w%gam, w%doo, w%dvv, w%lag_o, w%lag_v, w%grad)
      call zero_vector(int(w%no, int64)**2, w%doo)
      call zero_vector(int(w%nv, int64)**2, w%dvv)
      call zero_vector(3_int64*int(w%natm, int64), w%grad)

      ! --- 1. M = L L^T and B = (ia|P) L^-T --------------------------------
      call fitted_b(w, b, aux, pl, la)
      ! --- 2. amplitudes, Gamma', doo, dvv ----------------------------------
      call amplitude_pass(w, la)
      this%e_corr = chunked_dot(la, w%nov8, w%naux, w%gam, w%xr)
      ! --- 3. Gamma, gamma and the two-centre derivative --------------------
      call fitted_two_index(w, aux, la)
      ! B is spent: everything after this reads Gamma and the raw integrals.
      !$acc exit data delete(w%xr)
      deallocate (w%xr)
      ! --- 4. the Lagrangian and the three-centre derivative ----------------
      if (aux%spherical) then
         allocate (w%gamc(w%nov, w%naux_c))
         !$acc enter data create(w%gamc)
         call sph_cols_sph_to_cart(aux, la, w%nov, w%gam, w%gamc)
         !$acc exit data delete(w%gam)
         deallocate (w%gam)
         call lagrangian_pass(w, b, aux, pl, la, w%gamc)
         !$acc exit data delete(w%gamc)
         deallocate (w%gamc)
      else
         call lagrangian_pass(w, b, aux, pl, la, w%gam)
         !$acc exit data delete(w%gam)
         deallocate (w%gam)
      end if

      !$acc update self(w%doo, w%dvv, w%lag_o, w%lag_v, w%grad)
      this%doo = w%doo
      this%dvv = w%dvv
      this%lag_o = w%lag_o
      this%lag_v = w%lag_v
      this%grad = w%grad
      !$acc exit data delete(w%cmo, w%cocc, w%cvir, w%eo, w%ev, w%jm, w%doo, w%dvv, w%lag_o, w%lag_v, w%grad)
      if (opts%verbose) then
         dense = .false.
         dense([ST_TRANS, ST_TRSM, ST_KGEMM, ST_GAMMA, ST_DENSITY, ST_FIT, ST_LAG]) = .true.
         call la%sync()
         call w%tbl%report(tall%get_elapsed_time(), 0, 1, subtotal_mask=dense, subtotal_label="all dense algebra")
      end if
   end subroutine rimp2_relax_build

   subroutine fitted_b(w, b, aux, pl, la)
      !! w%jm := L, the metric's Cholesky factor (spherical when `aux` is),
      !! and w%xr := B = (ia|P) L^-T, rows a + (i-1) nvir -- trc_rimp2_run's
      !! own steps 1-4 on one rank, so this is the B that energy fitted.
      type(relax_work_t), intent(inout) :: w
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_linalg_t), intent(inout) :: la
      type(trc_basis_t) :: aux_subset
      real(dp), allocatable :: jmc(:, :), tens(:, :, :), tcart3(:, :, :), half(:, :, :)
      type(trc_gemm_batch_t) :: gemm_shape
      type(timer_type) :: tc
      integer :: n, no, nv, s0, s1, p0, np, p0c, npc

      n = w%n
      no = w%no
      nv = w%nv
      call w%tbl%step_on(la, ST_METRIC, tc)
      if (aux%spherical) then
         allocate (jmc(w%naux_c, w%naux_c))
         !$acc enter data create(jmc)
         call trc_df_2c(aux, jmc)
         call sph_metric_cart_to_sph(aux, la, jmc, w%jm)
         !$acc exit data delete(jmc)
         deallocate (jmc)
      else
         call trc_df_2c(aux, w%jm)
      end if
      call la%potrf(w%naux, w%jm, w%naux)
      call w%tbl%step_off(la, ST_METRIC, tc)
      call w%tbl%add_flop(ST_METRIC, int(w%naux, int64)**3/3)

      s0 = 1
      do while (s0 <= aux%nshell)
         call next_aux_block(aux, s0, w%nblk, s1, p0c, npc)
         call aux%subset(s0, s1, aux_subset)
         call aux_subset%to_device()
         if (aux%spherical) then
            p0 = aux%sh_ao_sph(s0)
            np = aux_subset%nao_sph
         else
            p0 = p0c
            np = npc
         end if
         allocate (tens(n, n, np), half(no, n, np))
         !$acc enter data create(tens, half)
         call w%tbl%step_on(la, ST_3C, tc)
         if (aux%spherical) then
            allocate (tcart3(n, n, npc))
            !$acc enter data create(tcart3)
            call trc_df_3c(b, pl, aux_subset, tcart3)
            call sph_aux3c_to_sph(aux_subset, n, tcart3, tens, la)
            !$acc exit data delete(tcart3)
            deallocate (tcart3)
         else
            call trc_df_3c(b, pl, aux_subset, tens)
         end if
         call w%tbl%step_off(la, ST_3C, tc)
         ! half(i, nu q) = C_act^T (mu nu|q); X(a + (i-1) nv, q) = C_vir^T half_q^T
         call w%tbl%step_on(la, ST_TRANS, tc)
         call la%gemm("T", "N", no, n*np, n, 1.0_dp, w%cocc, n, tens, n, 0.0_dp, half, no)
         gemm_shape = trc_gemm_batch(nv, no, n, np, sa=0_int64, sc=w%nov8, transa="T", transb="T")
         call la%gemm_batched(gemm_shape, w%cvir, half, w%xr(1, p0))
         call w%tbl%step_off(la, ST_TRANS, tc)
         call w%tbl%add_flop(ST_TRANS, int(np, int64)*(2*int(no, int64)*int(n, int64)**2 + &
                                                       2*int(nv, int64)*int(no, int64)*int(n, int64)))
         !$acc exit data delete(tens, half)
         deallocate (tens, half)
         call aux_subset%release()
         s0 = s1 + 1
      end do
      call w%tbl%step_on(la, ST_TRSM, tc)
      call la%trsm("R", "L", "T", w%nov, w%naux, 1.0_dp, w%jm, w%naux, w%xr, w%nov)
      call w%tbl%step_off(la, ST_TRSM, tc)
      call w%tbl%add_flop(ST_TRSM, w%nov8*int(w%naux, int64)**2)
   end subroutine fitted_b

   subroutine amplitude_pass(w, la)
      !! Step 2 of the header: w%gam := Gamma' = X B (rows a + (i-1) nv, the
      !! fitted auxiliary index), doo -= sum_k T_k^T X_k, dvv += X T^T, one
      !! group of occupied orbitals at a time.
      type(relax_work_t), intent(inout) :: w
      type(trc_linalg_t), intent(inout) :: la
      real(dp), allocatable :: ki(:, :), tt(:), xx(:)
      type(trc_gemm_batch_t) :: gemm_shape
      type(timer_type) :: tc
      integer :: ng, nb, i0, il, no, nv, nov, ldk
      integer(int64) :: slab

      no = w%no
      nv = w%nv
      nov = w%nov
      slab = int(nv, int64)*w%nov8
      ng = amplitude_group(no, nv)
      ldk = ng*nv
      allocate (ki(ldk, nov), tt(slab*ng), xx(slab*ng))
      !$acc enter data create(ki, tt, xx)
      do i0 = 1, no, ng
         nb = min(ng, no - i0 + 1)
         ! K(a + (il-1) nv, b + (j-1) nv) = (ia|jb), i = i0 - 1 + il
         call w%tbl%step_on(la, ST_KGEMM, tc)
         call la%gemm("N", "T", nb*nv, nov, w%naux, 1.0_dp, w%xr(1 + (i0 - 1)*nv, 1), nov, w%xr, nov, &
                      0.0_dp, ki, ldk)
         call w%tbl%step_off(la, ST_KGEMM, tc)
         call w%tbl%add_flop(ST_KGEMM, 2*int(nb*nv, int64)*w%nov8*int(w%naux, int64))
         call w%tbl%step_on(la, ST_AMP, tc)
         call amplitude_kernel(i0, nb, no, nv, ldk, nov, w%eo, w%ev, ki, tt, xx)
         call w%tbl%step_off(la, ST_AMP, tc)
         ! Gamma'(a + (i-1) nv, Q) = sum_(b,j) X_i(a, b + (j-1) nv) B(b + (j-1) nv, Q)
         call w%tbl%step_on(la, ST_GAMMA, tc)
         gemm_shape = trc_gemm_batch(nv, w%naux, nov, nb, ldb=nov, ldc=nov, sb=0_int64, sc=int(nv, int64))
         call la%gemm_batched(gemm_shape, xx, w%xr, w%gam(1 + (i0 - 1)*nv, 1))
         call w%tbl%step_off(la, ST_GAMMA, tc)
         call w%tbl%add_flop(ST_GAMMA, 2*int(nb, int64)*slab*int(w%naux, int64))
         ! dvv(b, a) += sum_(c, j, i) X(b, c, j, i) T(a, c, j, i): the group is
         ! one contiguous inner index; doo(i, j) -= sum_(a, b) T_k(ab, i) X_k(ab, j).
         call w%tbl%step_on(la, ST_DENSITY, tc)
         call la%gemm("N", "T", nv, nv, nv*no*nb, 1.0_dp, xx, nv, tt, nv, 1.0_dp, w%dvv, nv)
         do il = 1, nb
            call la%gemm("T", "N", no, no, nv*nv, -1.0_dp, tt(1 + (il - 1)*slab), nv*nv, &
                         xx(1 + (il - 1)*slab), nv*nv, 1.0_dp, w%doo, no)
         end do
         call w%tbl%step_off(la, ST_DENSITY, tc)
         call w%tbl%add_flop(ST_DENSITY, 2*int(nb, int64)*slab*int(nv + no, int64))
      end do
      !$acc exit data delete(ki, tt, xx)
      deallocate (ki, tt, xx)
   end subroutine amplitude_pass

   subroutine fitted_two_index(w, aux, la)
      !! Step 3: w%xr := c = B L^-1, w%gam := Gamma = Gamma' L^-1, gamma =
      !! Gamma^T c, and w%grad += sum_PQ Z_PQ d(P|Q)/dR with Z = -(gamma +
      !! gamma^T), expanded to the Cartesian auxiliary functions when `aux`
      !! is spherical.
      type(relax_work_t), intent(inout) :: w
      type(trc_basis_t), intent(in) :: aux
      type(trc_linalg_t), intent(inout) :: la
      real(dp), allocatable :: gpq(:, :), gpq_t(:, :), zmat(:, :), zhalf(:, :), zhalf_t(:, :)
      type(timer_type) :: tc
      integer :: naux, naux_c, nov

      naux = w%naux
      naux_c = w%naux_c
      nov = w%nov
      allocate (gpq(naux, naux), gpq_t(naux, naux))
      !$acc enter data create(gpq, gpq_t)
      call w%tbl%step_on(la, ST_FIT, tc)
      call la%trsm("R", "L", "N", nov, naux, 1.0_dp, w%jm, naux, w%xr, nov)
      call la%trsm("R", "L", "N", nov, naux, 1.0_dp, w%jm, naux, w%gam, nov)
      call la%gemm("T", "N", naux, naux, nov, 1.0_dp, w%gam, nov, w%xr, nov, 0.0_dp, gpq, naux)
      call la%transpose(naux, naux, gpq, naux, gpq_t, naux)
      call negated_sum(int(naux, int64)**2, gpq, gpq_t)
      call w%tbl%step_off(la, ST_FIT, tc)
      call w%tbl%add_flop(ST_FIT, 4*w%nov8*int(naux, int64)**2)
      ! gpq_t now holds Z = -(gamma + gamma^T), in the fitted space.
      allocate (zmat(naux_c, naux_c))
      !$acc enter data create(zmat)
      if (aux%spherical) then
         ! T Z T^T as two column expansions around a transpose (Z symmetric).
         allocate (zhalf(naux, naux_c), zhalf_t(naux_c, naux))
         !$acc enter data create(zhalf, zhalf_t)
         call sph_cols_sph_to_cart(aux, la, naux, gpq_t, zhalf)
         call la%transpose(naux, naux_c, zhalf, naux, zhalf_t, naux_c)
         call sph_cols_sph_to_cart(aux, la, naux_c, zhalf_t, zmat)
         !$acc exit data delete(zhalf, zhalf_t)
         deallocate (zhalf, zhalf_t)
      else
         call copy_vector(int(naux, int64)**2, gpq_t, zmat)
      end if
      call w%tbl%step_on(la, ST_D2C, tc)
      call df_grad_2c_driver(aux%nshell, naux_c, aux%maxnp, aux%sh_l, aux%sh_np, aux%sh_ao, aux%sh_at, &
                             aux%sh_e, aux%sh_c, aux%sh_r, w%natm, zmat, w%grad)
      call w%tbl%step_off(la, ST_D2C, tc)
      !$acc exit data delete(gpq, gpq_t, zmat)
      deallocate (gpq, gpq_t, zmat)
   end subroutine fitted_two_index

   subroutine lagrangian_pass(w, b, aux, pl, la, gamc)
      !! Step 4: the Lagrangian L', L'' and the three-centre derivative, from
      !! Gamma over the CARTESIAN auxiliary functions (`gamc`, rows
      !! a + (i-1) nv) and (mu nu|P) rebuilt block by block.
      type(relax_work_t), intent(inout) :: w
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_linalg_t), intent(inout) :: la
      real(dp), intent(in) :: gamc(w%nov, w%naux_c)
      type(trc_basis_t) :: aux_subset
      real(dp), allocatable :: tens(:, :, :), gg(:, :, :), gp(:, :, :), mm(:, :, :), hh(:, :, :)
      real(dp), allocatable :: ro(:, :), uu(:, :), cq0(:), pd(:, :)
      integer, allocatable :: q_a(:), q_b(:), q_o(:), q_n(:)
      integer :: coff(0:DFG_NPCLS)
      type(trc_gemm_batch_t) :: gemm_shape
      type(timer_type) :: tc
      integer :: n, no, nv, s0, s1, p0c, npc

      n = w%n
      no = w%no
      nv = w%nv
      allocate (ro(n, no), uu(nv, n), cq0(w%naux_c))
      !$acc enter data create(ro, uu, cq0)
      call zero_vector(int(n, int64)*int(no, int64), ro)
      call zero_vector(int(nv, int64)*int(n, int64), uu)
      call zero_vector(int(w%naux_c, int64), cq0)
      ! The three-centre kernels' pair plan over the energy's own pair list.
      call df_grad_3c_plan(b%nshell, b%maxnp, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, pl%npair, pl%p_i, &
                           pl%p_j, coff, q_a, q_b, q_o, q_n, pd)
      !$acc enter data copyin(q_a, q_b, q_o, q_n, pd)

      s0 = 1
      do while (s0 <= aux%nshell)
         call next_aux_block(aux, s0, w%nblk, s1, p0c, npc)
         call aux%subset(s0, s1, aux_subset)
         call aux_subset%to_device()
         allocate (tens(n, n, npc), gg(n, no, npc), gp(n, npc, no), mm(n, n, npc), hh(n, no, npc))
         !$acc enter data create(tens, gg, gp, mm, hh)
         call w%tbl%step_on(la, ST_3C, tc)
         call trc_df_3c(b, pl, aux_subset, tens)
         call w%tbl%step_off(la, ST_3C, tc)
         call w%tbl%step_on(la, ST_LAG, tc)
         ! G(nu, i, q) = sum_a C_vir(nu, a) Gamma(a, i, q)
         call la%gemm("N", "N", n, no*npc, nv, 1.0_dp, w%cvir, n, gamc(1, p0c), nv, 0.0_dp, gg, n)
         ! M(mu, nu, q) = 2 sum_i C_act(mu, i) G(nu, i, q) = 2 Gamma_ao^q
         gemm_shape = trc_gemm_batch(n, n, no, npc, sa=0_int64, transb="T")
         call la%gemm_batched(gemm_shape, w%cocc, gg, mm, alpha=2.0_dp)
         ! R_o(mu, i) += sum_(nu, q) (mu nu|q) G(nu, i, q), G reordered to (nu, q, i)
         call reorder_last_two(n, no, npc, gg, gp)
         call la%gemm("N", "N", n, no, n*npc, 1.0_dp, tens, n, gp, n*npc, 1.0_dp, ro, n)
         ! H(nu, i, q) = sum_mu (mu nu|q) C_act(mu, i); U(a, nu) += sum_(i, q) Gamma(a, i, q) H(nu, i, q)
         gemm_shape = trc_gemm_batch(n, no, n, npc, sb=0_int64, transa="T")
         call la%gemm_batched(gemm_shape, tens, w%cocc, hh)
         call la%gemm("N", "T", nv, n, no*npc, 1.0_dp, gamc(1, p0c), nv, hh, n, 1.0_dp, uu, nv)
         call w%tbl%step_off(la, ST_LAG, tc)
         call w%tbl%add_flop(ST_LAG, 2*int(npc, int64)*int(no, int64)*int(n, int64)* &
                                     (2*int(nv, int64) + 3*int(n, int64)))
         ! grad += sum (M + M^T)(mu, nu, q) d(mu nu|q)/dR = 4 Gamma_ao . d(mu nu|P).
         ! jc = 0: the kernels' own J part reads `cq0` (zeros) and `mm` as its
         ! density and contributes nothing.
         call w%tbl%step_on(la, ST_D3C, tc)
         call df_grad_3c_driver(b%nshell, n, b%sh_ao, b%sh_at, b%sh_r, coff, pl%npair, q_a, q_b, q_o, q_n, &
                                size(pd, 2), pd, aux_subset%nshell, aux_subset%maxnp, aux_subset%sh_l, &
                                aux_subset%sh_np, aux_subset%sh_ao, aux_subset%sh_at, aux_subset%sh_e, &
                                aux_subset%sh_c, aux_subset%sh_r, w%natm, 2, n, npc, mm, 0.0_dp, &
                                w%naux_c, 0, cq0, mm, w%grad)
         call w%tbl%step_off(la, ST_D3C, tc)
         !$acc exit data delete(tens, gg, gp, mm, hh)
         deallocate (tens, gg, gp, mm, hh)
         call aux_subset%release()
         s0 = s1 + 1
      end do
      ! L'(i, q) = 2 sum_mu R_o(mu, i) C(mu, q);  L''(a, q) = 2 sum_nu U(a, nu) C(nu, q)
      call w%tbl%step_on(la, ST_LAG, tc)
      call la%gemm("T", "N", no, w%nmo, n, 2.0_dp, ro, n, w%cmo, n, 0.0_dp, w%lag_o, no)
      call la%gemm("N", "N", nv, w%nmo, n, 2.0_dp, uu, nv, w%cmo, n, 0.0_dp, w%lag_v, nv)
      call w%tbl%step_off(la, ST_LAG, tc)
      !$acc exit data delete(q_a, q_b, q_o, q_n, pd, ro, uu, cq0)
      deallocate (q_a, q_b, q_o, q_n, pd, ro, uu, cq0)
   end subroutine lagrangian_pass

   subroutine next_aux_block(aux, s0, nblk, s1, p0, np)
      !! The run of whole auxiliary shells from s0 holding at most `nblk`
      !! Cartesian functions (at least one shell): last shell s1, first
      !! Cartesian function p0, Cartesian width np.
      type(trc_basis_t), intent(in) :: aux
      integer, intent(in) :: s0, nblk
      integer, intent(out) :: s1, p0, np
      integer :: ish
      s1 = s0
      do ish = s0 + 1, aux%nshell
         if (aux%sh_ao(ish) + ncart(aux%sh_l(ish)) - aux%sh_ao(s0) > nblk) exit
         s1 = ish
      end do
      p0 = aux%sh_ao(s0)
      np = aux%sh_ao(s1) + ncart(aux%sh_l(s1)) - p0
   end subroutine next_aux_block

   pure integer function ncart(l) result(res)
      integer, intent(in) :: l
      res = (l + 1)*(l + 2)/2
   end function ncart

   integer function amplitude_group(no, nv) result(ng)
      !! Occupied orbitals per amplitude group: NG_BEST, or as many as three
      !! (nvir, nocc nvir) buffers each fit in the device's free memory with
      !! half a gigabyte kept back for the libraries. At least one.
#ifdef TRC_DEVICE
      use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, &
                         acc_property_free_memory
#endif
      integer, intent(in) :: no, nv
#ifdef TRC_DEVICE
      integer(int64) :: free, per_i
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
      per_i = 24_int64*int(nv, int64)**2*int(no, int64)
      ng = int(min(int(NG_BEST, int64), max(1_int64, (free - 512_int64*1024*1024)/per_i)))
#else
      ng = NG_BEST
#endif
      ng = max(1, min(ng, no))
   end function amplitude_group

   function chunked_dot(la, nrow, ncol, x, y) result(s)
      !! sum(x * y) over an (nrow, ncol) pair, one panel per library call so
      !! no single call's count passes 2^30.
      type(trc_linalg_t), intent(in) :: la
      integer(int64), intent(in) :: nrow
      integer, intent(in) :: ncol
      real(dp), intent(in) :: x(*), y(*)
      real(dp) :: s
      integer(int64), parameter :: PANEL = 2_int64**30
      integer(int64) :: total, start, cnt
      total = nrow*int(ncol, int64)
      s = 0.0_dp
      start = 1
      do while (start <= total)
         cnt = min(PANEL, total - start + 1)
         s = s + la%dot(int(cnt), x(start), y(start))
         start = start + cnt
      end do
   end function chunked_dot

   ! ----------------------------------------------------------------------
   ! do concurrent kernels: explicit-shape dummies, device-resident arrays.
   ! ----------------------------------------------------------------------

   subroutine amplitude_kernel(i0, nb, no, nv, ldk, nov, eo, ev, ki, tt, xx)
      !! tt(a, b, j, il) = T_ij^ab and xx(a, b, j, il) = 2 T_ij^ab - T_ij^ba,
      !! i = i0 - 1 + il, from K(a + (il-1) nv, b + (j-1) nv) = (ia|jb).
      integer, intent(in) :: i0, nb, no, nv, ldk, nov
      real(dp), intent(in) :: eo(no), ev(nv), ki(ldk, nov)
      real(dp), intent(out) :: tt(nv, nv, no, nb), xx(nv, nv, no, nb)
      integer(int64) :: t, nv8, no8
      integer :: a, bb, j, il, i
      real(dp) :: denom, tab, tba
      nv8 = int(nv, int64)
      no8 = int(no, int64)
      do concurrent(t=0:nv8*nv8*no8*int(nb, int64) - 1) local(a, bb, j, il, i, denom, tab, tba)
         a = int(mod(t, nv8)) + 1
         bb = int(mod(t/nv8, nv8)) + 1
         j = int(mod(t/(nv8*nv8), no8)) + 1
         il = int(t/(nv8*nv8*no8)) + 1
         i = i0 - 1 + il
         denom = eo(i) + eo(j) - ev(a) - ev(bb)
         tab = ki((il - 1)*nv + a, (j - 1)*nv + bb)/denom
         tba = ki((il - 1)*nv + bb, (j - 1)*nv + a)/denom
         tt(a, bb, j, il) = tab
         xx(a, bb, j, il) = 2.0_dp*tab - tba
      end do
   end subroutine amplitude_kernel

   subroutine reorder_last_two(n, no, np, gg, gp)
      !! gp(nu, q, i) = gg(nu, i, q).
      integer, intent(in) :: n, no, np
      real(dp), intent(in) :: gg(n, no, np)
      real(dp), intent(out) :: gp(n, np, no)
      integer(int64) :: t, n8, np8
      integer :: nu, q, i
      n8 = int(n, int64)
      np8 = int(np, int64)
      do concurrent(t=0:n8*np8*int(no, int64) - 1) local(nu, q, i)
         nu = int(mod(t, n8)) + 1
         q = int(mod(t/n8, np8)) + 1
         i = int(t/(n8*np8)) + 1
         gp(nu, q, i) = gg(nu, i, q)
      end do
   end subroutine reorder_last_two

   subroutine negated_sum(n, x, y)
      !! y := -(x + y), elementwise.
      integer(int64), intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(inout) :: y(n)
      integer(int64) :: t
      do concurrent(t=1:n)
         y(t) = -(x(t) + y(t))
      end do
   end subroutine negated_sum

   subroutine zero_vector(n, x)
      !! x := 0, on the device.
      integer(int64), intent(in) :: n
      real(dp), intent(out) :: x(n)
      integer(int64) :: t
      do concurrent(t=1:n)
         x(t) = 0.0_dp
      end do
   end subroutine zero_vector

   subroutine copy_vector(n, x, y)
      !! y := x, on the device.
      integer(int64), intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: y(n)
      integer(int64) :: t
      do concurrent(t=1:n)
         y(t) = x(t)
      end do
   end subroutine copy_vector

end module trc_rimp2_relax
