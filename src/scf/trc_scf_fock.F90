!
! The Fock build, one iteration's worth: the two-electron part
! (build_g: g_rks / g_rks_incremental / g_uks, dispatched on nspin/
! incr_now), the exchange-correlation potential (build_vxc), the
! polarisable continuum's operator (build_pcm), the energy at the current
! density (scf_energy), assembling H + G (+ V_xc, + the continuum) into
! the Fock matrix (assemble_fock), the orthogonalised commutator FDS -
! SDF (commutator_error), and the Cartesian <-> spherical folds every one
! of them crosses through (to_cart/to_sph). See trc_scf_driver's own
! header for what each piece means physically; this module is just where
! the arithmetic now lives.
!
! EXPLICIT-SHAPE DUMMIES THROUGHOUT. Every array here is a plain
! argument, sized by an explicit dimension dummy (`nu_c`, `nlr_c`,
! `nin_c`, `nxc`, `npcm`: 1 when the caller has no use for the real
! array, matching trc_scf_state%alloc's own sizing). A `do concurrent`
! body never names a derived-type component -- trc_scf_state's `st%...`
! or trc_scf_context's `ctx%...` -- only one of these dummies: that is
! what lets this compile correctly on the GPU build once the caller's
! arrays are struct components instead of plain locals.
!
module trc_scf_fock
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use pic_mpi_lib, only: comm_t
   use pic_profiler, only: profiler_start, profiler_stop
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_sph, only: sph_cart_to_sph, sph_sph_to_cart
   use trc_eri, only: trc_eri_t
   use trc_pcm, only: trc_pcm_t
   use trc_xc, only: trc_xc_rks, trc_xc_uks
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_vv10, only: trc_vv10_rks, trc_vv10_uks
   use trc_linalg, only: trc_linalg_t
   use trc_scf_df, only: df_jk_t
   implicit none
   private

   public :: build_g, build_vxc, build_pcm, scf_energy, assemble_fock, commutator_error
   public :: to_cart, to_sph

contains

   !
   ! The two-electron build, resident and entirely in Cartesian space:
   ! g_rks (nspin = 1, a full build), g_rks_incremental (nspin = 1, from the
   ! density DIFFERENCE, see trc_scf_driver's own header) and g_uks (nspin = 2),
   ! dispatched on `nspin`/`incr_now` rather than three separately-named
   ! entry points, since the caller already has exactly one of them to ask
   ! for.
   !
   subroutine build_g(b, eri, nao_c, nspin, nu_c, nlr_c, nin_c, incr_now, incr_on, &
                      exx, exx_lr, omega, incr_scale_run, &
                      d_c, gmat_c, dtot_c, jmat_c, w1_c, klr_c, gacc_c, dref_c, ddif_c, &
                      df_on, dfj, lp)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: nao_c, nspin, nu_c, nlr_c, nin_c
      logical, intent(in) :: incr_now, incr_on
      real(dp), intent(in) :: exx, exx_lr, omega, incr_scale_run
      real(dp), intent(in) :: d_c(nao_c, nao_c, nspin)
      real(dp), intent(inout) :: gmat_c(nao_c, nao_c, nspin)
      real(dp), intent(inout) :: dtot_c(nu_c, nu_c), jmat_c(nu_c, nu_c)
      real(dp), intent(inout) :: w1_c(nao_c, nao_c), klr_c(nlr_c, nlr_c)
      real(dp), intent(inout) :: gacc_c(nin_c, nin_c), dref_c(nin_c, nin_c), ddif_c(nin_c, nin_c)
      logical, intent(in) :: df_on
      !! `--df`: density-fitted J and K instead of the four-centre engine.
      !! Restricted only (`nspin == 1`), never combined with `incr_on`
      !! (trc_scf_settings.F90's own refusal) -- the whole incremental
      !! block below is unreachable under it, so this branches BEFORE any
      !! of that rather than threading a flag through it.
      type(df_jk_t), intent(inout), pointer :: dfj
      !! Null and unused when `df_on` is false.
      type(trc_linalg_t), intent(inout) :: lp
      integer :: i, j, s

      if (df_on) then
         call dfj%build(lp, d_c(:, :, 1), gmat_c(:, :, 1), j_scale=1.0_dp, k_scale=exx)
         return
      end if

      if (nspin == 1) then
         if (incr_now) then
            do concurrent(i=1:nao_c, j=1:nao_c)
               ddif_c(i, j) = d_c(i, j, 1) - dref_c(i, j)
            end do
            !
            ! A fixed, tighter threshold on the increment, not one tracking ddif.
            !
            ! The Schwarz/density screen drops a quartet whose bound qq*dd is below
            ! an absolute threshold: fine for a full build, where dd is O(1). On
            ! ddif, quartets sit at the edge of the test, and which ones changes
            ! each iteration because ddif is a difference of iterates with a
            ! changing shape. The dropped sets differ and the errors do not average
            ! down as the full build's (always the same negligible quartets) do. On
            ! w172/6-31G (450 atoms, 1950 AOs) RMS(D) reaches 1e-10 while the
            ! energy bounces by 1e-4 to 1e-6 Eh per iteration.
            !
            ! Scaling `thresh` by ||ddif||/||dref|| is backwards: it loosens as
            ! ddif shrinks and loses the speed (25.6 s against 24.4 s baseline).
            ! A fixed tighter factor, `incr_thresh_scale` (default 1e-4), keeps
            ! the error floor low enough that |dE| clears 1e-8 at every reset
            ! through iteration 20, while each increment stays a few times cheaper
            ! than a full build.
            !
            block
               real(dp) :: thresh0
               thresh0 = eri%thresh
               eri%thresh = thresh0*incr_scale_run
               call eri%fock_resident(b, ddif_c, gmat_c(:, :, 1), k_scale=exx)
               if (exx_lr /= 0.0_dp) then
                  call eri%fock_resident(b, ddif_c, klr_c, j_scale=0.0_dp, k_scale=exx_lr, omega=omega)
                  do concurrent(i=1:nao_c, j=1:nao_c)
                     gmat_c(i, j, 1) = gmat_c(i, j, 1) + klr_c(i, j)
                  end do
               end if
               eri%thresh = thresh0
            end block
            do concurrent(i=1:nao_c, j=1:nao_c)
               gacc_c(i, j) = gacc_c(i, j) + gmat_c(i, j, 1)
               gmat_c(i, j, 1) = gacc_c(i, j)
               dref_c(i, j) = d_c(i, j, 1)
            end do
         else
            call eri%fock_resident(b, d_c(:, :, 1), gmat_c(:, :, 1), k_scale=exx)
            if (exx_lr /= 0.0_dp) then
               call eri%fock_resident(b, d_c(:, :, 1), klr_c, j_scale=0.0_dp, k_scale=exx_lr, omega=omega)
               do concurrent(i=1:nao_c, j=1:nao_c)
                  gmat_c(i, j, 1) = gmat_c(i, j, 1) + klr_c(i, j)
               end do
            end if
            ! The reference for the next increment, kept only when there
            ! will be one (gacc_c and dref_c are 1 x 1 otherwise).
            if (incr_on) then
               do concurrent(i=1:nao_c, j=1:nao_c)
                  gacc_c(i, j) = gmat_c(i, j, 1)
                  dref_c(i, j) = d_c(i, j, 1)
               end do
            end if
         end if
      else
         do concurrent(i=1:nao_c, j=1:nao_c)
            dtot_c(i, j) = d_c(i, j, 1) + d_c(i, j, 2)
         end do
         call eri%fock_resident(b, dtot_c, jmat_c, k_scale=0.0_dp)
         do s = 1, 2
            if (exx /= 0.0_dp .or. exx_lr /= 0.0_dp) then
               do concurrent(i=1:nao_c, j=1:nao_c)
                  w1_c(i, j) = 2.0_dp*d_c(i, j, s)
               end do
               call eri%fock_resident(b, w1_c, gmat_c(:, :, s), j_scale=0.0_dp, k_scale=exx)
               if (exx_lr /= 0.0_dp) then
                  call eri%fock_resident(b, w1_c, klr_c, j_scale=0.0_dp, k_scale=exx_lr, omega=omega)
                  do concurrent(i=1:nao_c, j=1:nao_c)
                     gmat_c(i, j, s) = gmat_c(i, j, s) + klr_c(i, j)
                  end do
               end if
               do concurrent(i=1:nao_c, j=1:nao_c)
                  gmat_c(i, j, s) = gmat_c(i, j, s) + jmat_c(i, j)
               end do
            else
               do concurrent(i=1:nao_c, j=1:nao_c)
                  gmat_c(i, j, s) = jmat_c(i, j)
               end do
            end if
         end do
      end if
   end subroutine build_g

   !
   ! The exchange-correlation potential and energy, finding its (Cartesian)
   ! density already present: the semilocal part (trc_xc_rks/uks), then
   ! VV10's nonlocal correlation on its own coarser grid when the functional
   ! has it, folded back to the run's own space. Only called when `dft` is
   ! true (the caller's own guard), so `vxc_c`/`vxc` are always full size
   ! here -- no 1 x 1 sizing dummy needed, unlike `build_g`'s arrays.
   !
   subroutine build_vxc(b, xg, xg_nlc, func, nlc, nspin, nao_c, nao, sph, rho_tol, comm, lp, &
                        d_c, vxc_c, vxc, exc, ngrid, tx_pts, tx_prs, t_vv10)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(inout) :: xg, xg_nlc
      type(trc_xc_functional_t), intent(in) :: func
      logical, intent(in) :: nlc, sph
      integer, intent(in) :: nspin, nao_c, nao
      real(dp), intent(in) :: rho_tol
      type(comm_t), intent(in), optional :: comm
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: d_c(nao_c, nao_c, nspin)
      real(dp), intent(inout) :: vxc_c(nao_c, nao_c, nspin), vxc(nao, nao, nspin)
      real(dp), intent(inout) :: exc, ngrid, tx_pts, tx_prs
      real(dp), intent(inout), optional :: t_vv10
      !! VV10's own elapsed time, added in (not overwritten): the DFT
      !! step table (`trc_scf_state`'s `xcst`) reports it as its own
      !! step, separate from the semilocal `tx_pts`/`tx_prs` above --
      !! see NVTX_RANGES below for why this is a plain wall time rather
      !! than a `step_table_t%step_on`/`step_off` pair.
      real(dp) :: e_nlc, tp1, tq1
      type(timer_type) :: tm_vv10
      integer :: s
      !
      ! NVTX_RANGES. Both the semilocal call and VV10's are bracketed here,
      ! at the one call site that has both the functional's `nlc` flag and
      ! `lp`. This file does not own a `step_table_t` -- `trc_scf_state`'s
      ! `xcst` does -- so these ranges are drawn with `pic_profiler`
      ! directly, the same two calls `step_table_t%step_on`/`step_off` make
      ! internally.
      !
      ! Both branches call `lp%sync()` before their own `profiler_stop`:
      ! without it, an async device launch left over from `trc_xc_rks`/
      ! `trc_xc_uks` could still be running when the "scf: xc" NVTX range
      ! closes, so a profiler would show it finishing inside "scf: vv10"
      ! instead -- the xc row undercounting, the vv10 row overcounting. The
      ! sync measured free here: `trc_xc_rks`/`trc_xc_uks` already force a
      ! device-to-host transfer of `exc`/`nelec` (`xc_drive`'s own `!$acc
      ! update self(eg, ng)`/`update self(vxc)`) before returning, so the
      ! device is already idle by the time `lp%sync()` runs. VV10's own
      ! sync (`tm_vv10`, a plain wall time handed back through `t_vv10` for
      ! `xcst` to accumulate, rather than a `step_table_t` tick) stays for
      ! the same reason.
      !
      call profiler_start("scf: xc", nvtx_only=.true.)
      if (nspin == 1) then
         call trc_xc_rks(b, xg, func, d_c(:, :, 1), vxc_c(:, :, 1), exc, ngrid, rho_tol=rho_tol, &
                         comm=comm, t_points=tp1, t_pairs=tq1, la=lp)
      else
         call trc_xc_uks(b, xg, func, d_c, vxc_c, exc, ngrid, rho_tol=rho_tol, &
                         comm=comm, t_points=tp1, t_pairs=tq1, la=lp)
      end if
      call lp%sync()
      call profiler_stop("scf: xc")
      ! VV10 adds its energy to exc and its potential into vxc_c, on its own
      ! (coarser) grid, before the spherical transform below.
      if (nlc) then
         call profiler_start("scf: vv10", nvtx_only=.true.)
         call lp%sync()
         call tm_vv10%start()
         if (nspin == 1) then
            call trc_vv10_rks(b, xg_nlc, d_c(:, :, 1), func%nlc_b, func%nlc_c, e_nlc, vxc_c(:, :, 1), la=lp)
         else
            call trc_vv10_uks(b, xg_nlc, d_c(:, :, 1), d_c(:, :, 2), func%nlc_b, func%nlc_c, e_nlc, &
                              vxc_c(:, :, 1), vxc_c(:, :, 2), la=lp)
         end if
         call lp%sync()
         call profiler_stop("scf: vv10")
         if (present(t_vv10)) t_vv10 = t_vv10 + tm_vv10%get_elapsed_time()
         exc = exc + e_nlc
      end if
      tx_pts = tx_pts + tp1
      tx_prs = tx_prs + tq1
      do s = 1, nspin
         call to_sph(b, sph, nao, nao_c, vxc_c(:, :, s), vxc(:, :, s))
      end do
   end subroutine build_vxc

   !
   ! The continuum's operator at the TOTAL density (pyscf's choice for UKS
   ! too: one set of charges, one operator added to both spins), folded back
   ! to the run's own space. Only called when `pcm_on` is true. `nu_c` is
   ! `dtot_c`'s true size (1 when `nspin == 1`, where it is never read).
   !
   subroutine build_pcm(b, pcm, pl, nspin, nao_c, nu_c, nao, sph, d_c, dtot_c, vpcm_c, vpcm, epcm, comm)
      type(trc_basis_t), intent(in) :: b
      type(trc_pcm_t), intent(inout) :: pcm
      type(trc_pairlist_t), intent(in) :: pl
      integer, intent(in) :: nspin, nao_c, nu_c, nao
      logical, intent(in) :: sph
      real(dp), intent(in) :: d_c(nao_c, nao_c, nspin), dtot_c(nu_c, nu_c)
      real(dp), intent(inout) :: vpcm_c(nao_c, nao_c), vpcm(nao, nao)
      real(dp), intent(inout) :: epcm
      type(comm_t), intent(in), optional :: comm
      ! An absent `comm` passes through as absent: one rank, unsplit.
      if (nspin == 1) then
         call pcm%operator(b, pl, d_c(:, :, 1), vpcm_c, epcm, comm)
      else
         call pcm%operator(b, pl, dtot_c, vpcm_c, epcm, comm)
      end if
      call to_sph(b, sph, nao, nao_c, vpcm_c, vpcm)
   end subroutine build_pcm

   ! The one- and two-electron energy at the density `dmat`: Tr(D H) and
   ! Tr(D G)/2, per spin; `etot` adds the nuclear repulsion, the XC and PCM
   ! energies the caller already has.
   subroutine scf_energy(lp, nao, nspin, dmat, hcore, gmat, e_nuc, exc, epcm, e1, e2, etot)
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: dmat(nao, nao, nspin), hcore(nao, nao), gmat(nao, nao, nspin)
      real(dp), intent(in) :: e_nuc, exc, epcm
      real(dp), intent(out) :: e1, e2, etot
      integer :: s
      e1 = 0.0_dp
      e2 = 0.0_dp
      do s = 1, nspin
         e1 = e1 + lp%dot_2d(dmat(:, :, s), hcore)
         e2 = e2 + 0.5_dp*lp%dot_2d(dmat(:, :, s), gmat(:, :, s))
      end do
      etot = e1 + e2 + exc + e_nuc + epcm
   end subroutine scf_energy

   subroutine assemble_fock(nao, nxc, npcm, nspin, dft, pcm_on, hcore, gmat, vxc, vpcm, fock)
      !! F = H + G (+ V_xc, + the continuum's operator), per spin. `vxc` is
      !! 1 x 1 (unused) when `.not. dft`, `vpcm` 1 x 1 when `.not. pcm_on`;
      !! `nxc`/`npcm` are their true sizes, so the explicit-shape dummies
      !! below match whichever the caller actually allocated.
      integer, intent(in) :: nao, nxc, npcm, nspin
      logical, intent(in) :: dft, pcm_on
      real(dp), intent(in) :: hcore(nao, nao)
      real(dp), intent(in) :: gmat(nao, nao, nspin), vxc(nxc, nxc, nspin)
      real(dp), intent(in) :: vpcm(npcm, npcm)
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      integer :: i, j, s
      do s = 1, nspin
         if (dft) then
            do concurrent(i=1:nao, j=1:nao)
               fock(i, j, s) = hcore(i, j) + gmat(i, j, s) + vxc(i, j, s)
            end do
         else
            do concurrent(i=1:nao, j=1:nao)
               fock(i, j, s) = hcore(i, j) + gmat(i, j, s)
            end do
         end if
         if (pcm_on) then
            do concurrent(i=1:nao, j=1:nao)
               fock(i, j, s) = fock(i, j, s) + vpcm(i, j)
            end do
         end if
      end do
   end subroutine assemble_fock

   ! W = F D S - S D F, errv = X^T W X, per spin: the orthogonalised
   ! commutator DIIS extrapolates on, and whose norm gates DIIS/A-DIIS.
   subroutine commutator_error(lp, nao, nspin, fock, dmat, smat, x, w1, w2, errv)
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: fock(nao, nao, nspin), dmat(nao, nao, nspin)
      real(dp), intent(in) :: smat(nao, nao), x(nao, nao)
      real(dp), intent(inout) :: w1(nao, nao), w2(nao, nao)
      real(dp), intent(inout) :: errv(nao, nao, nspin)
      integer :: s
      do s = 1, nspin
         ! w2 = F D S - S D F
         call lp%gemm_mat(fock(:, :, s), dmat(:, :, s), w1)
         call lp%gemm_mat(w1, smat, w2)
         call lp%gemm_mat(smat, dmat(:, :, s), w1)
         call lp%gemm_mat(w1, fock(:, :, s), w2, alpha=-1.0_dp, beta=1.0_dp)
         ! errv = X^T w2 X
         call lp%gemm_mat(w2, x, w1)
         call lp%gemm_mat(x, w1, errv(:, :, s), transa="T")
      end do
   end subroutine commutator_error

   subroutine to_cart(b, sph, nao, nao_c, a, a_c)
      !! A matrix from the SCF's space to the Cartesian one the kernels run
      !! in: T A T^T when the basis is spherical, a plain copy when it is not
      !! (no transform call at all on the Cartesian path, see trc_scf_driver's
      !! own header). `a` and `a_c` are the same size then.
      type(trc_basis_t), intent(in) :: b
      logical, intent(in) :: sph
      integer, intent(in) :: nao, nao_c
      real(dp), intent(in) :: a(nao, nao)
      real(dp), intent(out) :: a_c(nao_c, nao_c)
      integer :: i, j
      if (sph) then
         call sph_sph_to_cart(b, a, a_c)
      else
         do concurrent(i=1:nao_c, j=1:nao_c)
            a_c(i, j) = a(i, j)
         end do
      end if
   end subroutine to_cart

   subroutine to_sph(b, sph, nao, nao_c, a_c, a)
      !! The reverse fold: T^T A_cart T into the SCF's space, or a plain copy.
      type(trc_basis_t), intent(in) :: b
      logical, intent(in) :: sph
      integer, intent(in) :: nao, nao_c
      real(dp), intent(in) :: a_c(nao_c, nao_c)
      real(dp), intent(out) :: a(nao, nao)
      integer :: i, j
      if (sph) then
         call sph_cart_to_sph(b, a_c, a)
      else
         do concurrent(i=1:nao, j=1:nao)
            a(i, j) = a_c(i, j)
         end do
      end if
   end subroutine to_sph

end module trc_scf_fock
