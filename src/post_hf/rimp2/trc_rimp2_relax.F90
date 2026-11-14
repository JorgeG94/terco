module trc_rimp2_relax
   !!
   !! The correlation side of the RI-MP2 nuclear gradient, on the device.
   !!
   !! Stocks, Palethorpe and Barca, J. Chem. Theory Comput. 2024, 20, 2505,
   !! section 3, as metalquicha implements it (backends/cenzontle/correlation/
   !! mqc_czt_ri_mp2_gradient.f90, whose numpy prototype tools/cpu_validation/
   !! ri_mp2_gradient.py settled every index convention used here). The SCF
   !! side -- the Z-vector, the relaxed and energy-weighted densities, and the
   !! one-electron, overlap and four-centre derivative terms -- is
   !! src/post_hf/rimp2/trc_rimp2_grad.F90's; this module only builds what needs the
   !! fitted integrals:
   !!
   !!   T_ij^ab = (ia|jb) / (e_i + e_j - e_a - e_b),  (ia|jb) = sum_Q B_ia^Q B_jb^Q
   !!   X_ij^ab = 2 T_ij^ab - T_ij^ba
   !!   doo_ij  = - sum_kab T_ki^ab X_kj^ab          dvv_ba = sum_ijc T_ij^ca X_ij^cb
   !!   Gamma^P_ia = sum_jbQ X_ij^ab B_jb^Q (L^-1)_QP        (eq 8, Cholesky form)
   !!   gamma_PQ   = sum_ia Gamma^P_ia c_ia^Q,  c = B L^-1    (eq 10)
   !!   L'_iq  = 2 sum_aP Gamma^P_ia (q a|P)                  (eq 13)
   !!   L''_aq = 2 sum_iP Gamma^P_ia (i q|P)                  (eq 14)
   !!
   !! with B = (ia|P) L^-T and M = (P|Q) = L L^T, exactly the B trc_rimp2_run
   !! fits. Only the ACTIVE occupied orbitals (nfrozen + 1 .. nocc) enter the
   !! amplitudes, Gamma and doo; q runs over every molecular orbital, which is
   !! where the occupied-frozen rows of the Lagrangian come from. And the
   !! first two terms of eq 6, the only ones the fit changes:
   !!
   !!   dE/dR += 4 sum_(mu nu P) Gamma^P_mu nu d(mu nu|P)/dR - 2 sum_PQ gamma_PQ d(P|Q)/dR
   !!
   !! Gamma^P_mu nu = C_act Gamma^P C_vir^T, through trc_df_grad_kernels'
   !! three- and two-centre derivative kernels (the --df gradient's own, which
   !! take an arbitrary three-index Y and an arbitrary symmetric Z: no
   !! generator change).
   !!
   !! THE ORDER OF WORK
   !!   1. the metric and its Cholesky factor; B, one auxiliary block of
   !!      (mu nu|P) at a time (trc_rimp2_run's own single-rank steps 1-4).
   !!   2. the amplitudes, a group of occupied i at a time (a row of pairs):
   !!      K = B_i B_j^T only for the groups j from the row's own on, one kernel
   !!      for T and X (X also written back over K), and the earlier groups'
   !!      columns from the pair symmetry, T_ij^ab = T_ji^ba, out of mirror
   !!      tiles those groups' rows kept (amplitude_pass, pair_tile_plan): about
   !!      half of K, the largest step. Then Gamma' = X B for the whole row
   !!      (one GEMM the shape of K's), dvv (one GEMM over the row) and doo (one
   !!      GEMM per i) -- doo pairs (k, i) with (k, j), which is why each row is
   !!      completed rather than each tile used twice. E_corr = sum Gamma' . B,
   !!      and E_os, E_ss from T and X, come for free and are returned: the
   !!      energy step's numbers, so a caller with a gradient needs no separate
   !!      energy run, and can check the gradient differentiates its energy.
   !!   3. Gamma = Gamma' L^-1 (a triangular solve), gamma = L^-T (Gamma'^T B)
   !!      L^-1, and the two-centre derivative with Z = -(gamma + gamma^T).
   !!   4. (mu nu|P) again, block by block, for the Lagrangian and the three-
   !!      centre derivative together: G^P = C_vir Gamma^P, M^P = 2 C_act G^P^T
   !!      (the derivative kernels' one-sided Y, Y = M + M^T), R_o += (mu nu|P)
   !!      G^P, and U += Gamma^P H^P^T with H^P = (mu nu|P) C_act. Then
   !!      L' = 2 R_o^T C and L'' = 2 U C.
   !!
   !! A spherical auxiliary basis: B, Gamma and gamma live in its spherical
   !! space, as the energy's do; Gamma is expanded onto the Cartesian
   !! auxiliary functions once before step 4 (trc_sph's sph_cols_sph_to_cart)
   !! and gamma's Z before the two-centre kernels, so the kernels and the raw
   !! (mu nu|P) of step 4 stay Cartesian. The orbital basis is always the
   !! Cartesian AO space here: a spherical caller passes C_cart = T C_sph.
   !!
   !! MEMORY. B and Gamma are nocc nvir naux doubles each, plus three
   !! group-sized (nvir, nvir, nocc, ng) amplitude buffers; the group size is
   !! chosen from the device's free memory (amplitude_group), and the mirror
   !! tiles, (ng nvir)^2 doubles each, take what is left: about (nocc nvir)^2/4
   !! doubles hold them all (gly10/cc-pVDZ: 9.6 GB; next to a --df SCF's own
   !! fit on a 32 GB V100, 6.6 GB were free, a band of 9 of the 14 later
   !! groups). With less, the band of groups kept narrows and the rest of K is
   !! computed again, as without the symmetry -- the numbers do not change.
   !! One rank.
   !!
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_cart, only: ncart => cart_ncart_of
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
                         ST_TILE = 7, ST_PAIR = 8, ST_GAMMA = 9, ST_DENSITY = 10, ST_FIT = 11, ST_D2C = 12, &
                         ST_LAG = 13, ST_D3C = 14, NSTEP = 14
   character(len=22), parameter :: STEP_NAME(NSTEP) = [character(len=22) :: &
                                   "metric + Cholesky", "(mu nu|P)", "transform (ia|P)", "B = X L^-T", &
                                   "K_i GEMMs", "amplitudes T, X", "mirror tiles", "pair energies", &
                                   "Gamma' = X B (gemm)", &
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
      !! The RI-MP2 correlation energy these amplitudes give, sum X . (ia|jb)
      !! (= Gamma' . B).
      real(dp) :: e_os = 0.0_dp
      !! Its opposite-spin part, sum (ia|jb)^2 / D, summed from the same
      !! amplitudes (trc_rimp2_run's E_os, so a caller with a gradient needs
      !! no separate energy step).
      real(dp) :: e_ss = 0.0_dp
      !! Its same-spin part, sum [(ia|jb)^2 - (ia|jb)(ib|ja)] / D.
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
      integer :: group_size = 0
      !! Occupied orbitals per amplitude group; 0 sizes it from the free
      !! memory (amplitude_group). For tests.
      integer :: tile_slots = -1
      !! Mirror tiles the pool may hold; negative sizes it from the free
      !! memory, 0 recomputes every column (the pass without the pair
      !! symmetry). For tests.
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
      integer :: group_size = 0, tile_slots = -1
      logical :: verbose = .false.
      integer(int64) :: nov8 = 0
      type(step_table_t) :: tbl
      real(dp), allocatable :: cmo(:, :), cocc(:, :), cvir(:, :), eo(:), ev(:), jm(:, :)
      real(dp), allocatable :: xr(:, :)
      !! B: (nocc nvir, naux), rows a + (i-1) nvir.
      real(dp), allocatable :: gam(:, :)
      !! Gamma', then Gamma: (nocc nvir, naux).
      real(dp), allocatable :: gamc(:, :)
      !! Gamma on the Cartesian auxiliary functions (a spherical aux only).
      real(dp), allocatable :: doo(:, :), dvv(:, :), lag_o(:, :), lag_v(:, :), grad(:, :)
      real(dp) :: e_os = 0.0_dp, e_ss = 0.0_dp
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
      w%group_size = this%group_size
      w%tile_slots = this%tile_slots
      w%verbose = opts%verbose
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
      this%e_os = w%e_os
      this%e_ss = w%e_ss
      ! --- 3. Gamma, gamma and the two-centre derivative --------------------
      call fitted_two_index(w, aux, la)
      ! B is spent: everything after this reads Gamma and the raw integrals.
      !$acc exit data delete(w%xr)
      deallocate (w%xr)
      ! --- 4. the Lagrangian and the three-centre derivative ----------------
      if (aux%spherical) then
         allocate (w%gamc(w%nov, w%naux_c))
         !$acc enter data create(w%gamc)
         call sph_cols_sph_to_cart(aux, w%nov, w%gam, w%gamc)
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
      !! fitted auxiliary index), doo -= sum_k T_k^T X_k, dvv += X T^T and
      !! E_os, E_ss, one group of occupied orbitals (a row of pairs) at a
      !! time. Each group's row needs T and X for EVERY j, but K only for
      !! j from the group on: the columns of earlier groups are mirror
      !! tiles, T_ij^ab = T_ji^ba, kept from those groups' own rows in a
      !! pool (pair_tile_plan) -- what the pool cannot hold is computed
      !! again, as every column was before.
      type(relax_work_t), intent(inout) :: w
      type(trc_linalg_t), intent(inout) :: la
      real(dp), allocatable :: ki(:, :), tt(:), xx(:), part(:, :), pool(:, :)
      integer, allocatable :: slot_of(:, :), free_slot(:)
      type(timer_type) :: tc
      integer :: ng, ngrp, nb, i0, p, q, band, nslot, nfree, jstore, c1, no, nv, nov, ldk, npart
      integer(int64) :: slab, tile

      no = w%no
      nv = w%nv
      nov = w%nov
      slab = int(nv, int64)*w%nov8
      call amplitude_group(no, nv, w%group_size, w%tile_slots, ng, nslot)
      ngrp = (no + ng - 1)/ng
      call pair_tile_plan(ngrp, nslot, band)
      if (w%verbose) then
         call global_logger%performance("  rimp2 grad: "//to_string(ng)//" occupied per group, "// &
                                        to_string(ngrp)//" groups; mirror tiles kept "//to_string(band)// &
                                        " groups ahead, "//to_string(nslot)//" at once ("// &
                                        to_string(8.0_dp*real(nslot, dp)*(real(ng, dp)*real(nv, dp))**2/ &
                                                  1024.0_dp**3, "f8.2")//" GiB)")
      end if
      ldk = ng*nv
      tile = int(nv, int64)**2*int(ng, int64)**2
      if (nslot == 0) tile = 1
      allocate (ki(ldk, nov), tt(slab*ng), xx(slab*ng), part(2, nov*ng), pool(tile, max(nslot, 1)))
      allocate (slot_of(ngrp, ngrp), free_slot(max(nslot, 1)))
      !$acc enter data create(ki, tt, xx, part, pool)
      slot_of = 0
      nfree = nslot
      free_slot = [(q, q=1, max(nslot, 1))]
      w%e_os = 0.0_dp
      w%e_ss = 0.0_dp
      do p = 1, ngrp
         i0 = (p - 1)*ng + 1
         nb = min(ng, no - i0 + 1)
         ! Groups jstore .. p - 1 come from the pool; columns 1 .. c1 (the
         ! groups beyond the band) and i0 .. no from K.
         jstore = max(1, p - band)
         c1 = (jstore - 1)*ng
         call w%tbl%step_on(la, ST_KGEMM, tc)
         if (c1 == i0 - 1) then
            call k_columns(w, la, i0, nb, 1, no, ldk, ki)
         else
            if (c1 > 0) call k_columns(w, la, i0, nb, 1, c1, ldk, ki)
            call k_columns(w, la, i0, nb, i0, no, ldk, ki)
         end if
         call w%tbl%step_off(la, ST_KGEMM, tc)
         call w%tbl%step_on(la, ST_AMP, tc)
         if (c1 == i0 - 1) then
            call amplitude_kernel(i0, nb, 1, no, no, nv, ldk, nov, w%eo, w%ev, ki, tt, xx)
         else
            if (c1 > 0) call amplitude_kernel(i0, nb, 1, c1, no, nv, ldk, nov, w%eo, w%ev, ki, tt, xx)
            call amplitude_kernel(i0, nb, i0, no - i0 + 1, no, nv, ldk, nov, w%eo, w%ev, ki, tt, xx)
         end if
         call w%tbl%step_off(la, ST_AMP, tc)
         ! The mirror tiles: T and X of (i in p, j in q) from T of (j, i).
         call w%tbl%step_on(la, ST_TILE, tc)
         do q = jstore, p - 1
            call mirror_fill((q - 1)*ng + 1, ng, nb, no, nv, ldk, nov, pool(1, slot_of(q, p)), ki, tt, xx)
            nfree = nfree + 1
            free_slot(nfree) = slot_of(q, p)
         end do
         call w%tbl%step_off(la, ST_TILE, tc)
         call w%tbl%add_flop(ST_TILE, 3*int(nv, int64)**2*int(nb, int64)*int(i0 - 1 - c1, int64))
         ! E_os and E_ss of the row's pairs, which the energy step would
         ! otherwise compute again from its own K (trc_rimp2_run).
         call w%tbl%step_on(la, ST_PAIR, tc)
         ! Only j >= i0: an earlier group's pairs were counted by its own row.
         npart = nv*(no - i0 + 1)*nb
         call pair_energy_kernel(i0, nb, no, nv, w%eo, w%ev, tt, xx, part)
         !$acc update self(part(:, 1:npart))
         w%e_os = w%e_os + sum(part(1, 1:npart))
         w%e_ss = w%e_ss + sum(part(2, 1:npart))
         call w%tbl%step_off(la, ST_PAIR, tc)
         call w%tbl%add_flop(ST_PAIR, 7*int(npart, int64)*int(nv, int64))
         call row_densities(w, la, i0, nb, ldk, ki, tt, xx)
         ! This row's T for the later groups in the band, into the pool.
         call w%tbl%step_on(la, ST_TILE, tc)
         do q = p + 1, min(p + band, ngrp)
            slot_of(p, q) = free_slot(nfree)
            nfree = nfree - 1
            call tile_store((q - 1)*ng + 1, min(ng, no - (q - 1)*ng), ng, nb, no, nv, tt, pool(1, slot_of(p, q)))
         end do
         call w%tbl%step_off(la, ST_TILE, tc)
      end do
      !$acc exit data delete(ki, tt, xx, part, pool)
      deallocate (ki, tt, xx, part, pool, slot_of, free_slot)
   end subroutine amplitude_pass

   subroutine k_columns(w, la, i0, nb, j0, j1, ldk, ki)
      !! K(a + (il-1) nv, b + (j-1) nv) = (ia|jb) for i = i0 - 1 + il and the
      !! occupied columns j = j0 .. j1, one GEMM.
      type(relax_work_t), intent(inout) :: w
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: i0, nb, j0, j1, ldk
      real(dp), intent(inout) :: ki(ldk, w%nov)
      integer :: nv
      nv = w%nv
      call la%gemm("N", "T", nb*nv, (j1 - j0 + 1)*nv, w%naux, 1.0_dp, w%xr(1 + (i0 - 1)*nv, 1), w%nov, &
                   w%xr(1 + (j0 - 1)*nv, 1), w%nov, 0.0_dp, ki(1, 1 + (j0 - 1)*nv), ldk)
      call w%tbl%add_flop(ST_KGEMM, 2*int(nb*nv, int64)*int((j1 - j0 + 1)*nv, int64)*int(w%naux, int64))
   end subroutine k_columns

   subroutine row_densities(w, la, i0, nb, ldk, ki, tt, xx)
      !! Gamma' for the row's i, and its doo and dvv terms, from X in K's
      !! layout (`ki`) and T, X as (a, b, j, il) (`tt`, `xx`).
      type(relax_work_t), intent(inout) :: w
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: i0, nb, ldk
      real(dp), intent(in) :: ki(ldk, w%nov), tt(*), xx(*)
      type(timer_type) :: tc
      integer :: il, no, nv, nov
      integer(int64) :: slab
      no = w%no
      nv = w%nv
      nov = w%nov
      slab = int(nv, int64)*w%nov8
      ! Gamma'(a + (i-1) nv, Q) = sum_(b,j) X(a + (il-1) nv, b + (j-1) nv) B(b + (j-1) nv, Q):
      ! the kernels left X in K's own layout, so the group's rows of Gamma'
      ! are one GEMM the shape of K's. A strided batch of nv-row items
      ! (one per i) ran at 4.9 TFLOP/s on a V100 (gly10/cc-pVDZ), cuBLAS
      ! picking its 64x64 tile; this one runs at K's 7.0.
      call w%tbl%step_on(la, ST_GAMMA, tc)
      call la%gemm("N", "N", nb*nv, w%naux, nov, 1.0_dp, ki, ldk, w%xr, nov, 0.0_dp, &
                   w%gam(1 + (i0 - 1)*nv, 1), nov)
      call w%tbl%step_off(la, ST_GAMMA, tc)
      call w%tbl%add_flop(ST_GAMMA, 2*int(nb, int64)*slab*int(w%naux, int64))
      ! dvv(b, a) += sum_(c, j, i) X(b, c, j, i) T(a, c, j, i): the group is
      ! one contiguous inner index; doo(i, j) -= sum_(a, b) T_k(ab, i) X_k(ab, j).
      ! doo couples two pairs sharing an index (k, i) and (k, j), which is
      ! why each row is completed from the mirror tiles before this.
      call w%tbl%step_on(la, ST_DENSITY, tc)
      call la%gemm("N", "T", nv, nv, nv*no*nb, 1.0_dp, xx, nv, tt, nv, 1.0_dp, w%dvv, nv)
      do il = 1, nb
         call la%gemm("T", "N", no, no, nv*nv, -1.0_dp, tt(1 + (il - 1)*slab), nv*nv, &
                      xx(1 + (il - 1)*slab), nv*nv, 1.0_dp, w%doo, no)
      end do
      call w%tbl%step_off(la, ST_DENSITY, tc)
      call w%tbl%add_flop(ST_DENSITY, 2*int(nb, int64)*slab*int(nv + no, int64))
   end subroutine row_densities

   subroutine pair_tile_plan(ngrp, nslot, band)
      !! The widest band of mirror tiles the pool can hold: group p's row
      !! keeps its T for the groups p + 1 .. p + band, each tile freed when
      !! that group's row reads it. `nslot` comes back as the most tiles
      !! ever held at once (band 0: none, every column recomputed).
      integer, intent(in) :: ngrp
      integer, intent(inout) :: nslot
      integer, intent(out) :: band
      integer :: trial, peak
      band = 0
      do trial = 1, ngrp - 1
         peak = band_peak(ngrp, trial)
         if (peak > nslot) exit
         band = trial
      end do
      nslot = band_peak(ngrp, band)
   end subroutine pair_tile_plan

   pure integer function band_peak(ngrp, band) result(peak)
      !! Tiles held right after group p's row stores its own -- (q, r) with
      !! q <= p < r <= q + band -- at the worst p.
      integer, intent(in) :: ngrp, band
      integer :: p, q, held
      peak = 0
      do p = 1, ngrp
         held = 0
         do q = max(1, p - band + 1), p
            held = held + max(0, min(q + band, ngrp) - p)
         end do
         peak = max(peak, held)
      end do
   end function band_peak

   subroutine fitted_two_index(w, aux, la)
      !! Step 3: w%gam := Gamma = Gamma' L^-1, gamma = Gamma^T c with
      !! c = B L^-1, and w%grad += sum_PQ Z_PQ d(P|Q)/dR with Z = -(gamma +
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
      ! gamma = (Gamma' L^-1)^T (B L^-1) = L^-T (Gamma'^T B) L^-1: the
      ! (nocc nvir) x naux^2 product first and the two solves on the
      ! naux^2 result, so c = B L^-1 is never formed -- one full-size
      ! triangular solve fewer (gly10/cc-pVDZ: 0.41 -> 0.27 s).
      call la%gemm("T", "N", naux, naux, nov, 1.0_dp, w%gam, nov, w%xr, nov, 0.0_dp, gpq, naux)
      call la%trsm("L", "L", "T", naux, naux, 1.0_dp, w%jm, naux, gpq, naux)
      call la%trsm("R", "L", "N", naux, naux, 1.0_dp, w%jm, naux, gpq, naux)
      call la%trsm("R", "L", "N", nov, naux, 1.0_dp, w%jm, naux, w%gam, nov)
      call la%transpose(naux, naux, gpq, naux, gpq_t, naux)
      call negated_sum(int(naux, int64)**2, gpq, gpq_t)
      call w%tbl%step_off(la, ST_FIT, tc)
      call w%tbl%add_flop(ST_FIT, 3*w%nov8*int(naux, int64)**2 + 2*int(naux, int64)**3)
      ! gpq_t now holds Z = -(gamma + gamma^T), in the fitted space.
      allocate (zmat(naux_c, naux_c))
      !$acc enter data create(zmat)
      if (aux%spherical) then
         ! T Z T^T as two column expansions around a transpose (Z symmetric).
         allocate (zhalf(naux, naux_c), zhalf_t(naux_c, naux))
         !$acc enter data create(zhalf, zhalf_t)
         call sph_cols_sph_to_cart(aux, naux, gpq_t, zhalf)
         call la%transpose(naux, naux_c, zhalf, naux, zhalf_t, naux_c)
         call sph_cols_sph_to_cart(aux, naux_c, zhalf_t, zmat)
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

   subroutine amplitude_group(no, nv, group_size, tile_slots, ng, nslot)
      !! Occupied orbitals per group, `ng`: `group_size` when positive, else
      !! NG_BEST, or as many as three (nvir, nocc nvir) buffers each fit in
      !! the device's free memory with half a gigabyte kept back for the
      !! libraries; at least one. And `nslot`, the mirror tiles of (ng nvir)^2
      !! the pool may hold: `tile_slots` when not negative, else what the
      !! rest of the free memory takes (on the host, no limit:
      !! pair_tile_plan caps it at the whole band). Shrinking ng so that
      !! every tile fits does not pay: gly10/cc-pVDZ --df on a V100 at ng = 6
      !! (all tiles) against 8 (a band of 9 of 14 groups) saved 0.20 s of K
      !! and lost 0.21 s of Gamma' to the smaller GEMMs.
#ifdef TRC_DEVICE
      use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, &
                         acc_property_free_memory
#endif
      integer, intent(in) :: no, nv, group_size, tile_slots
      integer, intent(out) :: ng, nslot
      integer(int64), parameter :: HOST_TILE_BUDGET = 2_int64*1024*1024*1024
      !! The host build's mirror-tile pool: the same 2 GB trc_scf_df_response's
      !! cache_budget keeps on the host. Every tile is ~(nocc nvir)^2/4 doubles
      !! in total (9.6 GB at gly10/cc-pVDZ), which a CPU/OpenMP run that fit
      !! before the tiles existed need not have.
      integer(int64) :: tile_bytes
#ifdef TRC_DEVICE
      integer(int64) :: free, per_i
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64) - 512_int64*1024*1024
      per_i = 24_int64*int(nv, int64)**2*int(no, int64)
      ng = int(min(int(NG_BEST, int64), max(1_int64, free/per_i)))
#else
      ng = NG_BEST
#endif
      if (group_size > 0) ng = group_size
      ng = max(1, min(ng, no))
      tile_bytes = 8_int64*(int(ng, int64)*int(nv, int64))**2
#ifdef TRC_DEVICE
      nslot = int(min(int(huge(nslot), int64), max(0_int64, (free - per_i*ng)/tile_bytes)))
#else
      nslot = int(min(int(huge(nslot), int64), HOST_TILE_BUDGET/tile_bytes))
#endif
      if (tile_slots >= 0) nslot = tile_slots
   end subroutine amplitude_group

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

   subroutine amplitude_kernel(i0, nb, j0, nj, no, nv, ldk, nov, eo, ev, ki, tt, xx)
      !! ----------------------------------------------------------------------
      !! do concurrent kernels: explicit-shape dummies, device-resident arrays.
      !! ----------------------------------------------------------------------
      !! tt(a, b, j, il) = T_ij^ab and xx(a, b, j, il) = 2 T_ij^ab - T_ij^ba,
      !! i = i0 - 1 + il and j = j0 .. j0 + nj - 1, from K(a + (il-1) nv,
      !! b + (j-1) nv) = (ia|jb), and X written back over K in K's own layout
      !! (Gamma's GEMM reads it there). One item per unordered pair a <= b,
      !! which reads both (ia|jb) and (ib|ja) before writing either, so the
      !! in-place update has no race.
      integer, intent(in) :: i0, nb, j0, nj, no, nv, ldk, nov
      real(dp), intent(in) :: eo(no), ev(nv)
      real(dp), intent(inout) :: ki(ldk, nov)
      real(dp), intent(inout) :: tt(nv, nv, no, nb), xx(nv, nv, no, nb)
      integer(int64) :: t, nv8, nj8, tlast
      integer :: a, bb, j, il, i, ra, rb, ca, cb
      real(dp) :: denom, tab, tba, xab, xba
      nv8 = int(nv, int64)
      nj8 = int(nj, int64)
      tlast = nv8*nv8*nj8*int(nb, int64) - 1
      do concurrent(t=0:tlast) local(a, bb, j, il, i, ra, rb, ca, cb, denom, tab, tba, xab, xba)
         a = int(mod(t, nv8)) + 1
         bb = int(mod(t/nv8, nv8)) + 1
         j = j0 + int(mod(t/(nv8*nv8), nj8))
         il = int(t/(nv8*nv8*nj8)) + 1
         if (a <= bb) then
            i = i0 - 1 + il
            ra = (il - 1)*nv + a
            rb = (il - 1)*nv + bb
            ca = (j - 1)*nv + a
            cb = (j - 1)*nv + bb
            denom = eo(i) + eo(j) - ev(a) - ev(bb)
            tab = ki(ra, cb)/denom
            tba = ki(rb, ca)/denom
            xab = 2.0_dp*tab - tba
            xba = 2.0_dp*tba - tab
            tt(a, bb, j, il) = tab
            tt(bb, a, j, il) = tba
            xx(a, bb, j, il) = xab
            xx(bb, a, j, il) = xba
            ki(ra, cb) = xab
            ki(rb, ca) = xba
         end if
      end do
   end subroutine amplitude_kernel

   subroutine mirror_fill(j0, ng, nb, no, nv, ldk, nov, tile, ki, tt, xx)
      !! The row's columns j = j0 .. j0 + ng - 1 (an earlier, whole group)
      !! from that group's own row: tile(a, b, il, jl) = T_ji^ab with
      !! i = the row's il-th orbital and j = j0 - 1 + jl, so T_ij^ab =
      !! tile(b, a, il, jl) and X_ij^ab = 2 tile(b, a, ..) - tile(a, b, ..),
      !! written where amplitude_kernel would have put them (tt, xx, and X
      !! over K's layout in ki). One item per (a, b), a fastest: every write
      !! and one of the two reads are unit stride.
      integer, intent(in) :: j0, ng, nb, no, nv, ldk, nov
      real(dp), intent(in) :: tile(nv, nv, ng, ng)
      real(dp), intent(inout) :: ki(ldk, nov)
      real(dp), intent(inout) :: tt(nv, nv, no, nb), xx(nv, nv, no, nb)
      integer(int64) :: t, nv8, ng8, tlast
      integer :: a, bb, jl, il, j
      real(dp) :: tab, tba
      nv8 = int(nv, int64)
      ng8 = int(ng, int64)
      tlast = nv8*nv8*ng8*int(nb, int64) - 1
      do concurrent(t=0:tlast) local(a, bb, jl, il, j, tab, tba)
         a = int(mod(t, nv8)) + 1
         bb = int(mod(t/nv8, nv8)) + 1
         jl = int(mod(t/(nv8*nv8), ng8)) + 1
         il = int(t/(nv8*nv8*ng8)) + 1
         j = j0 - 1 + jl
         tab = tile(bb, a, il, jl)
         tba = tile(a, bb, il, jl)
         tt(a, bb, j, il) = tab
         xx(a, bb, j, il) = 2.0_dp*tab - tba
         ki((il - 1)*nv + a, (j - 1)*nv + bb) = 2.0_dp*tab - tba
      end do
   end subroutine mirror_fill

   subroutine tile_store(c0, ncol, ng, nb, no, nv, tt, tile)
      !! tile(a, b, cl, il) = tt(a, b, c0 - 1 + cl, il) = T_ij^ab for the
      !! row's il-th i and the later group's cl-th j: what that group's row
      !! reads back as its mirror (mirror_fill). A copy, unit stride.
      integer, intent(in) :: c0, ncol, ng, nb, no, nv
      real(dp), intent(in) :: tt(nv*nv, no, nb)
      real(dp), intent(out) :: tile(nv*nv, ng, ng)
      integer(int64) :: t, nvv8, nc8
      integer :: ab, cl, il
      nvv8 = int(nv, int64)**2
      nc8 = int(ncol, int64)
      do concurrent(t=0:nvv8*nc8*int(nb, int64) - 1) local(ab, cl, il)
         ab = int(mod(t, nvv8)) + 1
         cl = int(mod(t/nvv8, nc8)) + 1
         il = int(t/(nvv8*nc8)) + 1
         tile(ab, cl, il) = tt(ab, c0 - 1 + cl, il)
      end do
   end subroutine tile_store

   subroutine pair_energy_kernel(i0, nb, no, nv, eo, ev, tt, xx, part)
      !! part(1, a + (j-i0) nv + (il-1) nv (no-i0+1)) = w sum_b T^2 D and
      !! part(2, ..) = w sum_b T (X - T) D for j = i0 .. no, D = e_i + e_j -
      !! e_a - e_b, i = i0 - 1 + il: with T = (ia|jb)/D and X - T = T -
      !! T_ij^ba these are the pair's (ia|jb)^2/D and [(ia|jb)^2 -
      !! (ia|jb)(ib|ja)]/D summed over b, the energy step's E_os and E_ss
      !! terms. w = 2 for j beyond the row's own group, whose (j, i) pairs
      !! are the same numbers and are never summed; 1 inside it, where both
      !! orders are. One item per (a, j, il), a fastest, so both reads are
      !! unit stride; the host adds the items.
      integer, intent(in) :: i0, nb, no, nv
      real(dp), intent(in) :: eo(no), ev(nv)
      real(dp), intent(in) :: tt(nv, nv, no, nb), xx(nv, nv, no, nb)
      real(dp), intent(out) :: part(2, nv*(no - i0 + 1)*nb)
      integer :: t, a, j, il, b, ncol
      real(dp) :: eija, d, tab, pos, pss, w
      ncol = no - i0 + 1
      do concurrent(t=1:nv*ncol*nb) local(a, j, il, b, eija, d, tab, pos, pss, w)
         a = mod(t - 1, nv) + 1
         j = i0 + mod((t - 1)/nv, ncol)
         il = (t - 1)/(nv*ncol) + 1
         eija = eo(i0 - 1 + il) + eo(j) - ev(a)
         w = merge(1.0_dp, 2.0_dp, j < i0 + nb)
         pos = 0.0_dp
         pss = 0.0_dp
         do b = 1, nv
            d = eija - ev(b)
            tab = tt(a, b, j, il)
            pos = pos + tab*tab*d
            pss = pss + tab*(xx(a, b, j, il) - tab)*d
         end do
         part(1, t) = w*pos
         part(2, t) = w*pss
      end do
   end subroutine pair_energy_kernel

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
