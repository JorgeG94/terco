module trc_rccsd_ints
   !!
   !! The density-fitted integral blocks of a closed-shell RI-CCSD(T) calculation.
   !!
   !! With the metric J = (P|Q) = L L^T and B^P_pq = sum_Q (pq|Q) (L^-T)_QP, every
   !! two-electron integral the amplitude equations need is
   !!
   !!     (pq|rs) = sum_P B^P_pq B^P_rs.
   !!
   !! Orbitals are the RHF canonical ones, the first `nfrozen` occupied orbitals
   !! are left out of the correlation, and the index ranges are
   !!
   !!     i, j, k, l   the ACTIVE occupied orbitals   (no = nocc - nfrozen)
   !!     a, b, c      the virtual orbitals           (nv)
   !!
   !! The blocks, as Fortran arrays (first index fastest), are
   !!
   !!     oooo(i,j,k,l) = (ij|kl)      ooov(i,j,k,a) = (ij|ka)
   !!     oovv(i,j,a,b) = (ij|ab)      ovov(i,a,j,b) = (ia|jb)
   !!
   !! and the three-index factors that make them keep their compound index with
   !! the FIRST orbital fastest:
   !!
   !!     b_oo(i + (j-1) no, P)    b_ov(i + (a-1) no, P)    b_vv(a + (b-1) nv, P)
   !!
   !! so that b_oo b_ov^T is ooov seen as a (no*no, no*nv) matrix, and so on.
   !! (ab|cd) is never formed: the ladder term builds it from `b_vv` a batch of
   !! (cd) columns at a time. Nor is (ia|bc) (o v^3): every contraction that needs it
   !! is done from b_ov and b_vv directly (trc_rccsd_inter's N1, the W stages, the
   !! T1 and T2 terms, the ladder dressing, the (T) slices). `fill_ovvv` builds the
   !! dense block into a caller's array for the two legacy (T) routines and the tests.
   !! B enters only as B B^T, so the Cholesky triangular
   !! solve is as good as any other factor (the eigen-decomposed J^-1/2 would
   !! also do).
   !!
   !! The 3c pipeline is that of the RI-MP2 driver (trc_rimp2_driver.F90): the same
   !! `trc_df_2c` metric and Cholesky factor, the same `trc_df_3c` kernels over
   !! blocks of auxiliary shells with the same spherical P-fold, the same
   !! orbital fold for a spherical basis. Only the half-transform differs:
   !! RI-MP2 keeps (ia|P) alone, here the occupied-occupied, occupied-virtual and
   !! virtual-virtual blocks come out of one half-transformed tensor, each by a
   !! strided batch of GEMMs.
   !!
   !! RESIDENCY. After `build` every array of the type exists on the host AND on
   !! the device (`enter data` in `build`, `exit data` in `release`, both on the
   !! same `ints%name` expression), and the two copies are equal. A consumer that
   !! rewrites one of them (none does yet) must `update` the other side before
   !! anyone reads it. Each component is mapped on its own and the parent `ints`
   !! is never mapped, so a kernel gets a block only as an explicit-shape dummy
   !! argument: `ints%...` named inside a `do concurrent` body is an illegal
   !! address under mem:separate, which the host build does not show (the rule
   !! `relax_work_t` in trc_rimp2_relax follows).
   !!
   !! MEMORY. The largest array is b_vv (nv^2 naux: 396 MB at nv 265, naux 705),
   !! then the o^2 v^2 blocks; the host mirror of every array is kept, and with
   !! `aux_block` = 0 the three-index transients hold the whole auxiliary basis at
   !! once.
   !!
   !! BUILD STAMP. Every successful `build` gives the type a new, nonzero `stamp`
   !! (`release` resets it to 0). A consumer that packs integral-only arrays once
   !! at its own `init` (the scratch types of trc_rccsd_inter and trc_rccsd_w)
   !! records the stamp it saw and refuses to run when it no longer matches: after
   !! `ints%build` is called again, re-`init` the scratch.
   !!
   !! RANKS. `build` makes no MPI call: every rank of a multi-rank run builds the whole of it from its own orbitals
   !! (replicated), see trc_rccsd_driver.
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use pic_profiler, only: profiler_start, profiler_stop
   use pic_blas_interfaces, only: pic_gemm_x
   use trc_rccsd_steps, only: RCCSD_ST_METRIC, RCCSD_ST_CHOLESKY, RCCSD_ST_3C, RCCSD_ST_HALF, RCCSD_ST_TRSM, &
                              RCCSD_ST_BLOCKS, RCCSD_STEP_NAME, rccsd_steps_gemm_flops, rccsd_steps_trsm_flops, &
                              rccsd_steps_potrf_flops
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_cart, only: ncart => cart_ncart_of
   use trc_integrals, only: trc_df_2c, trc_df_3c
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_sph, only: sph_metric_cart_to_sph, sph_mo_to_cart, sph_aux3c_to_sph
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   integer :: build_counter = 0  !! last stamp handed out; only ever increases, so a rebuilt type never repeats one

   public :: trc_rccsd_ints_t, trc_rccsd_ints_options_t

   type :: trc_rccsd_ints_options_t
      !! What `build` needs to know about the reference besides its orbitals. Plain values only.
      integer :: nalpha = 0     !! doubly occupied orbitals (alpha electrons)
      integer :: nbeta = 0      !! beta electrons; anything but `nalpha` is an open shell and is refused
      integer :: nfrozen = 0    !! the lowest `nfrozen` occupied orbitals are left out of the correlation
      integer :: aux_block = 0  !! auxiliary functions per (mu nu|P) block; 0 means all at once (no number changes)
      logical :: has_ecp = .false.  !! the calculation carries an effective core potential: refused
   end type trc_rccsd_ints_options_t

   type :: trc_rccsd_ints_t
      !! The fitted MO integral blocks of one RHF reference, see the module header
      type(error_t) :: status  !! set by `build`: a refusal leaves its message here and the type empty
      integer :: no = 0       !! active occupied orbitals
      integer :: nv = 0       !! virtual orbitals
      integer :: naux = 0     !! auxiliary functions (spherical count when the auxiliary basis is)
      integer :: nfrozen = 0  !! frozen occupied orbitals, the lowest by index
      integer :: stamp = 0    !! 0 when empty, else a new value per successful `build`; scratch types compare it
      real(dp), allocatable :: oooo(:, :, :, :)   !! (ij|kl)
      real(dp), allocatable :: ooov(:, :, :, :)   !! (ij|ka)
      real(dp), allocatable :: oovv(:, :, :, :)   !! (ij|ab)
      real(dp), allocatable :: ovov(:, :, :, :)   !! (ia|jb)
      real(dp), allocatable :: b_oo(:, :)         !! (ij, P), i fastest
      real(dp), allocatable :: b_ov(:, :)         !! (ia, P), i fastest
      real(dp), allocatable :: b_vv(:, :)         !! (ab, P), a fastest
      real(dp), allocatable :: eps_o(:)           !! active occupied orbital energies
      real(dp), allocatable :: eps_v(:)           !! virtual orbital energies
      real(dp) :: seconds(RCCSD_ST_METRIC:RCCSD_ST_BLOCKS) = 0.0_dp
      !! device seconds of the last `build` per piece, indexed by the RCCSD_ST_* row (trc_rccsd_steps) it fills
      integer(int64) :: flops(RCCSD_ST_METRIC:RCCSD_ST_BLOCKS) = 0_int64
      !! 2 m n k summed over the GEMMs (and the solves and the factorisation) that piece executed, as called
   contains
      procedure :: build => rccsd_ints_build
      procedure :: release => rccsd_ints_release
      procedure :: mp2_energy => rccsd_ints_mp2_energy
      procedure :: fill_ovvv => rccsd_ints_fill_ovvv
   end type trc_rccsd_ints_t

contains

   subroutine rccsd_ints_release(ints)
      !! Frees the host and device copies; safe on a type that holds nothing
      class(trc_rccsd_ints_t), intent(inout) :: ints

      if (allocated(ints%oooo)) then
         !$acc exit data delete(ints%oooo)
         deallocate (ints%oooo)
      end if
      if (allocated(ints%ooov)) then
         !$acc exit data delete(ints%ooov)
         deallocate (ints%ooov)
      end if
      if (allocated(ints%oovv)) then
         !$acc exit data delete(ints%oovv)
         deallocate (ints%oovv)
      end if
      if (allocated(ints%ovov)) then
         !$acc exit data delete(ints%ovov)
         deallocate (ints%ovov)
      end if
      if (allocated(ints%b_oo)) then
         !$acc exit data delete(ints%b_oo)
         deallocate (ints%b_oo)
      end if
      if (allocated(ints%b_ov)) then
         !$acc exit data delete(ints%b_ov)
         deallocate (ints%b_ov)
      end if
      if (allocated(ints%b_vv)) then
         !$acc exit data delete(ints%b_vv)
         deallocate (ints%b_vv)
      end if
      if (allocated(ints%eps_o)) then
         !$acc exit data delete(ints%eps_o)
         deallocate (ints%eps_o)
      end if
      if (allocated(ints%eps_v)) then
         !$acc exit data delete(ints%eps_v)
         deallocate (ints%eps_v)
      end if
      ints%no = 0
      ints%nv = 0
      ints%naux = 0
      ints%nfrozen = 0
      ints%stamp = 0
      ints%seconds = 0.0_dp
      ints%flops = 0_int64
   end subroutine rccsd_ints_release

   subroutine piece_on(la, k, tc)
      !! Opens piece `k` of the build: once the device is idle the clock starts and the NVTX range opens, as a step
      !! table's `step_on` does (the build cannot take the driver's table: it already has eight arguments)
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: k
      type(timer_type), intent(out) :: tc

      call la%sync()
      call tc%start()
      call profiler_start("rccsd: "//trim(adjustl(RCCSD_STEP_NAME(k))), nvtx_only=.true.)
   end subroutine piece_on

   subroutine piece_off(ints, la, k, tc, flops)
      !! Closes piece `k` after the same sync and adds its time and the flops it executed to `ints`
      class(trc_rccsd_ints_t), intent(inout) :: ints
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: k
      type(timer_type), intent(in) :: tc
      integer(int64), intent(in) :: flops

      call la%sync()
      ints%seconds(k) = ints%seconds(k) + tc%get_elapsed_time()
      ints%flops(k) = ints%flops(k) + flops
      call profiler_stop("rccsd: "//trim(adjustl(RCCSD_STEP_NAME(k))))
   end subroutine piece_off

   subroutine rccsd_ints_check(b, aux, cmo, eps, opts, error)
      !! Every refusal, before anything is allocated or mapped
      type(trc_basis_t), intent(in) :: b, aux
      real(dp), intent(in) :: cmo(:, :), eps(:)
      type(trc_rccsd_ints_options_t), intent(in) :: opts
      type(error_t), intent(inout) :: error

      integer :: nmo, nalpha, nbeta, nfrozen

      nmo = merge(b%nao_sph, b%nao, b%spherical)
      nalpha = opts%nalpha
      nbeta = opts%nbeta
      nfrozen = opts%nfrozen
      if (opts%has_ecp) then
         call error%set(ERROR_VALIDATION, "trc_rccsd_ints: an ECP is not supported")
         return
      end if
      if (nalpha /= nbeta) then
         call error%set(ERROR_VALIDATION, "trc_rccsd_ints: an open shell is not supported, nalpha /= nbeta")
         return
      end if
      if (nfrozen < 0 .or. nfrozen > nalpha) then
         call error%set(ERROR_VALIDATION, "trc_rccsd_ints: the frozen-core count must lie between 0 and nocc")
         return
      end if
      if (nalpha - nfrozen < 1) then
         call error%set(ERROR_VALIDATION, "trc_rccsd_ints: no active occupied orbital to correlate")
         return
      end if
      if (nmo - nalpha < 1) then
         call error%set(ERROR_VALIDATION, "trc_rccsd_ints: no virtual orbital to correlate")
         return
      end if
      if (size(cmo, 1) /= nmo .or. size(cmo, 2) /= nmo .or. size(eps) /= nmo) then
         call error%set(ERROR_VALIDATION, "trc_rccsd_ints: the orbitals are not sized by the basis's own space "// &
                        "(spherical count for a spherical basis)")
         return
      end if
      if (maxval(b%sh_l) > DF_LMAX) then
         call error%set(ERROR_VALIDATION, "trc_rccsd_ints: the orbital basis carries angular momentum above "// &
                        "what the density-fitting kernels were generated for (gen_df.py --lmax)")
         return
      end if
      if (maxval(aux%sh_l) > DF_LMAX_AUX) then
         call error%set(ERROR_VALIDATION, "trc_rccsd_ints: the auxiliary basis carries angular momentum above "// &
                        "what the density-fitting kernels were generated for (gen_df.py --lmax-aux)")
         return
      end if
   end subroutine rccsd_ints_check

   subroutine rccsd_ints_build(ints, b, aux, pl, cmo, eps, la, opts)
      !! Builds every block from the RHF orbitals. `cmo` (nmo, nmo) and `eps` (nmo) come from the SCF,
      !! occupied first, in the space of the basis (the SPHERICAL count for a spherical basis: no
      !! cartesian padding as `trc_rimp2_run` needs). The basis and pair list are device-resident
      !! as the SCF leaves them; `la` belongs to the caller, who releases it.
      !! Success or refusal is reported in `ints%status` (a refusal leaves the type empty).
      class(trc_rccsd_ints_t), intent(inout) :: ints
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: cmo(:, :), eps(:)
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_ints_options_t), intent(in) :: opts

      integer :: nalpha, nfrozen, nao, nmo, no, nv, nact, naux, naux_c, nblk, s0, s1, ish, p0, np, p0c, npc
      real(dp), allocatable :: jm(:, :), jmc(:, :), tens(:, :, :), tcart3(:, :, :), half(:, :, :)
      real(dp), allocatable :: csph(:, :), cact(:, :)
      type(trc_basis_t) :: aux_subset
      type(trc_gemm_batch_t) :: shape
      type(timer_type) :: tc

      call ints%release()
      call ints%status%clear()
      call rccsd_ints_check(b, aux, cmo, eps, opts, ints%status)
      if (ints%status%has_error()) return
      nalpha = opts%nalpha
      nfrozen = opts%nfrozen

      ! Nothing below can return early: every `enter data` has its `exit data`
      ! before the routine ends. The kernels and library calls stop the program
      ! on their own failures (a metric that is not positive definite, an
      ! integral class the build lacks).
      nao = b%nao
      nmo = merge(b%nao_sph, b%nao, b%spherical)
      no = nalpha - nfrozen
      nv = nmo - nalpha
      nact = no + nv
      naux_c = aux%nao
      naux = merge(aux%nao_sph, aux%nao, aux%spherical)
      nblk = naux_c
      if (opts%aux_block > 0) nblk = opts%aux_block
      ints%no = no
      ints%nv = nv
      ints%naux = naux
      ints%nfrozen = nfrozen

      allocate (ints%eps_o(no), ints%eps_v(nv))
      ints%eps_o = eps(nfrozen + 1:nalpha)
      ints%eps_v = eps(nalpha + 1:nmo)
      !$acc enter data copyin(ints%eps_o, ints%eps_v)

      ! Active orbitals in the Cartesian AO space the 3c kernels work in,
      ! occupied then virtual: columns 1:no and no+1:nact.
      allocate (cact(nao, nact))
      if (b%spherical) then
         allocate (csph(nmo, nact))
         csph = cmo(:, nfrozen + 1:nmo)
         !$acc enter data copyin(csph) create(cact)
         call sph_mo_to_cart(b, nact, csph, cact)
         !$acc exit data delete(csph)
         deallocate (csph)
      else
         cact = cmo(:, nfrozen + 1:nmo)
         !$acc enter data copyin(cact)
      end if

      allocate (ints%b_oo(no*no, naux), ints%b_ov(no*nv, naux), ints%b_vv(nv*nv, naux))
      ! Touch the host mirror before the closing `update self` writes it (why: at the rank-four blocks).
      ints%b_oo = 0.0_dp
      ints%b_ov = 0.0_dp
      ints%b_vv = 0.0_dp
      !$acc enter data create(ints%b_oo, ints%b_ov, ints%b_vv)

      ! The metric and its Cholesky factor, L in the lower triangle.
      call la%init(naux)
      call piece_on(la, RCCSD_ST_METRIC, tc)
      allocate (jm(naux, naux))
      !$acc enter data create(jm)
      if (aux%spherical) then
         allocate (jmc(naux_c, naux_c))
         !$acc enter data create(jmc)
         call trc_df_2c(aux, jmc)
         call sph_metric_cart_to_sph(aux, la, jmc, jm)
         !$acc exit data delete(jmc)
         deallocate (jmc)
      else
         call trc_df_2c(aux, jm)
      end if
      call piece_off(ints, la, RCCSD_ST_METRIC, tc, 0_int64)
      call piece_on(la, RCCSD_ST_CHOLESKY, tc)
      call la%potrf(naux, jm, naux)
      call piece_off(ints, la, RCCSD_ST_CHOLESKY, tc, rccsd_steps_potrf_flops(naux))

      ! X(pq, P) = sum_mu,nu C(mu,p) (mu nu|P) C(nu,q) for the three blocks,
      ! a block of auxiliary shells at a time.
      s0 = 1
      do while (s0 <= aux%nshell)
         s1 = s0
         do ish = s0 + 1, aux%nshell
            if (aux%sh_ao(ish) + ncart(aux%sh_l(ish)) - aux%sh_ao(s0) > nblk) exit
            s1 = ish
         end do
         p0c = aux%sh_ao(s0)
         npc = aux%sh_ao(s1) + ncart(aux%sh_l(s1)) - p0c
         call aux%subset(s0, s1, aux_subset)
         call aux_subset%to_device()
         if (aux%spherical) then
            p0 = aux%sh_ao_sph(s0)
            np = aux_subset%nao_sph
         else
            p0 = p0c
            np = npc
         end if
         allocate (tens(nao, nao, np), half(nact, nao, np))
         !$acc enter data create(tens, half)
         call piece_on(la, RCCSD_ST_3C, tc)
         if (aux%spherical) then
            allocate (tcart3(nao, nao, npc))
            !$acc enter data create(tcart3)
            call trc_df_3c(b, pl, aux_subset, tcart3)
            call sph_aux3c_to_sph(aux_subset, nao, tcart3, tens, la)
            !$acc exit data delete(tcart3)
            deallocate (tcart3)
         else
            call trc_df_3c(b, pl, aux_subset, tens)
         end if
         call piece_off(ints, la, RCCSD_ST_3C, tc, 0_int64)
         call piece_on(la, RCCSD_ST_HALF, tc)
         ! half(p, nu, P) = sum_mu C(mu, p) (mu nu|P), all active p at once
         call la%gemm("T", "N", nact, nao*np, nao, 1.0_dp, cact, nao, tens, nao, 0.0_dp, half, nact)
         ! then one strided batch over P per block: out(p, q, P) = sum_nu half(p, nu, P) C(nu, q)
         shape = trc_gemm_batch(no, no, nao, np, lda=nact, ldb=nao, ldc=no, &
                                sa=int(nact, int64)*nao, sb=0_int64, sc=int(no, int64)*no)
         call la%gemm_batched(shape, half, cact, ints%b_oo(1, p0))
         shape = trc_gemm_batch(no, nv, nao, np, lda=nact, ldb=nao, ldc=no, &
                                sa=int(nact, int64)*nao, sb=0_int64, sc=int(no, int64)*nv)
         call la%gemm_batched(shape, half, cact(1, no + 1), ints%b_ov(1, p0))
         shape = trc_gemm_batch(nv, nv, nao, np, lda=nact, ldb=nao, ldc=nv, &
                                sa=int(nact, int64)*nao, sb=0_int64, sc=int(nv, int64)*nv)
         call la%gemm_batched(shape, half(no + 1, 1, 1), cact(1, no + 1), ints%b_vv(1, p0))
         call piece_off(ints, la, RCCSD_ST_HALF, tc, rccsd_steps_gemm_flops(nact, nao*np, nao) + &
                        rccsd_steps_gemm_flops(no, no, nao, np) + rccsd_steps_gemm_flops(no, nv, nao, np) + &
                        rccsd_steps_gemm_flops(nv, nv, nao, np))
         !$acc exit data delete(tens, half)
         deallocate (tens, half)
         call aux_subset%release()
         s0 = s1 + 1
      end do
      !$acc exit data delete(cact)
      deallocate (cact)

      ! B = X L^-T: solve B L^T = X from the right.
      call piece_on(la, RCCSD_ST_TRSM, tc)
      call la%trsm("R", "L", "T", no*no, naux, 1.0_dp, jm, naux, ints%b_oo, no*no)
      call la%trsm("R", "L", "T", no*nv, naux, 1.0_dp, jm, naux, ints%b_ov, no*nv)
      call la%trsm("R", "L", "T", nv*nv, naux, 1.0_dp, jm, naux, ints%b_vv, nv*nv)
      call piece_off(ints, la, RCCSD_ST_TRSM, tc, rccsd_steps_trsm_flops(no*no, naux) + rccsd_steps_trsm_flops(no*nv, naux) + &
                     rccsd_steps_trsm_flops(nv*nv, naux))
      !$acc exit data delete(jm)
      deallocate (jm)

      ! (pq|rs) = sum_P B_pq,P B_rs,P. The rank-four destinations go to the
      ! GEMM by sequence association: no reshaping temporary.
      allocate (ints%oooo(no, no, no, no), ints%ooov(no, no, no, nv), ints%oovv(no, no, nv, nv))
      allocate (ints%ovov(no, nv, no, nv))
      ! The host pages must exist before the `update self` below. On a GH200 a
      ! pageable device-to-host copy into memory that was allocated but never
      ! touched faults the pages in from the GPU side: the 2.98 GB (ia|bc) of a
      ! (no 20, nv 265) run, when it was built here, took 10.6 s that way, 0.01 s
      ! once touched, while the touch itself costs 0.1 s on the host. An x86 host
      ! stages the copy and does not show it. Zero on the host, before the GEMMs:
      ! on the host build these arrays are the GEMM output itself.
      ints%oooo = 0.0_dp
      ints%ooov = 0.0_dp
      ints%oovv = 0.0_dp
      ints%ovov = 0.0_dp
      !$acc enter data create(ints%oooo, ints%ooov, ints%oovv, ints%ovov)
      call piece_on(la, RCCSD_ST_BLOCKS, tc)
      call la%gemm("N", "T", no*no, no*no, naux, 1.0_dp, ints%b_oo, no*no, ints%b_oo, no*no, 0.0_dp, ints%oooo, no*no)
      call la%gemm("N", "T", no*no, no*nv, naux, 1.0_dp, ints%b_oo, no*no, ints%b_ov, no*nv, 0.0_dp, ints%ooov, no*no)
      call la%gemm("N", "T", no*no, nv*nv, naux, 1.0_dp, ints%b_oo, no*no, ints%b_vv, nv*nv, 0.0_dp, ints%oovv, no*no)
      call la%gemm("N", "T", no*nv, no*nv, naux, 1.0_dp, ints%b_ov, no*nv, ints%b_ov, no*nv, 0.0_dp, ints%ovov, no*nv)
      call piece_off(ints, la, RCCSD_ST_BLOCKS, tc, rccsd_steps_gemm_flops(no*no, no*no, naux) + &
                     rccsd_steps_gemm_flops(no*no, no*nv, naux) + rccsd_steps_gemm_flops(no*no, nv*nv, naux) + &
                     rccsd_steps_gemm_flops(no*nv, no*nv, naux))
      ! Into host pages already touched at allocation (why: at the rank-four blocks).
      !$acc update self(ints%b_oo, ints%b_ov, ints%b_vv, ints%oooo, ints%ooov, ints%oovv, ints%ovov)
      build_counter = build_counter + 1
      ints%stamp = build_counter
   end subroutine rccsd_ints_build

   subroutine rccsd_ints_fill_ovvv(ints, ovvv, la)
      !! Fills the caller's array with ovvv(i,a,b,c) = (ia|bc) = sum_P b_ov(ia,P) b_vv(bc,P), one GEMM
      !! (m = o v, n = v^2, k = naux). With `la` the device copies are used and the DEVICE copy of `ovvv` is
      !! written (the caller has mapped it); without, the host copies and the host BLAS fill the HOST copy.
      !!
      !! For the consumers that still want the dense block (the (T) reference and the virtual-triple device
      !! routine, and the tests): it is o v^3 doubles, which is what the CCSD iterations no longer hold.
      class(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(inout) :: ovvv(ints%no*ints%nv, ints%nv*ints%nv)
      !! (no, nv, nv, nv) in the caller, by sequence association: the GEMM's (o v, v^2) matrix
      type(trc_linalg_t), intent(in), optional :: la

      integer :: nov, nv2

      nov = ints%no*ints%nv
      nv2 = ints%nv*ints%nv
      if (present(la)) then
         call la%gemm("N", "T", nov, nv2, ints%naux, 1.0_dp, ints%b_ov, nov, ints%b_vv, nv2, 0.0_dp, ovvv, nov)
      else
         call pic_gemm_x("N", "T", nov, nv2, ints%naux, 1.0_dp, ints%b_ov, nov, ints%b_vv, nv2, 0.0_dp, ovvv, nov)
      end if
   end subroutine rccsd_ints_fill_ovvv

   function rccsd_ints_mp2_energy(ints) result(energy)
      !! E_MP2 = sum_ijab (ia|jb) [2 (ia|jb) - (ib|ja)] / (e_i + e_j - e_a - e_b)
      !!
      !! From the host copy of `ovov`, summed in one fixed order (i fastest,
      !! then a, j, b), so the number does not depend on the thread count.
      class(trc_rccsd_ints_t), intent(in) :: ints
      real(dp) :: energy

      integer :: i, j, a, b
      real(dp) :: iajb

      energy = 0.0_dp
      do b = 1, ints%nv
         do j = 1, ints%no
            do a = 1, ints%nv
               do i = 1, ints%no
                  iajb = ints%ovov(i, a, j, b)
                  energy = energy + iajb*(2.0_dp*iajb - ints%ovov(i, b, j, a))/ &
                           (ints%eps_o(i) + ints%eps_o(j) - ints%eps_v(a) - ints%eps_v(b))
               end do
            end do
         end do
      end do
   end function rccsd_ints_mp2_energy

end module trc_rccsd_ints
