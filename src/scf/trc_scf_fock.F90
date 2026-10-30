!
! The Fock build, one iteration's worth: the two-electron part
! (scf_fock_build_g: g_rks / g_rks_incremental / g_uks, dispatched on nspin/
! incr_now), the exchange-correlation potential (scf_fock_build_vxc), the
! polarisable continuum's operator (scf_fock_build_pcm), the energy at the current
! density (scf_fock_energy), assembling H + G (+ V_xc, + the continuum) into
! the Fock matrix (scf_fock_assemble), the orthogonalised commutator FDS -
! SDF (scf_fock_commutator_error), and the Cartesian <-> spherical folds every one
! of them crosses through (scf_fock_to_cart/scf_fock_to_sph). See trc_scf_driver's own
! header for what each piece means physically; this module is just where
! the arithmetic now lives.
!
! THE CALLER PASSES OBJECTS, THE KERNELS SEE ARRAYS. Each public routine
! takes the run's basis, context (`scf_ctx_t`), settings (`scf_settings_t`)
! and per-iteration state (`scf_state_t`) and picks out what it needs, so
! the SCF loop reads as one line per step. What the run borrows rather
! than builds -- the linear-algebra handles, the ERI engine and pair list,
! the DFT grids, the density fit -- comes as its own argument. A `do concurrent` body must never
! name a derived-type component, though: under `-gpu=mem:separate` that is
! what breaks. So `scf_fock_build_g` and `scf_fock_assemble`, the two with loops of their
! own, hand the arrays to a private explicit-shape worker (`*_arrays`), whose
! sizes come from the arrays themselves (1 x 1 where the run has no use for
! one, matching trc_scf_state%alloc).
!
module trc_scf_fock
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use pic_mpi_lib, only: comm_t
   use pic_profiler, only: profiler_start, profiler_stop
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_sph, only: sph_cart_to_sph, sph_sph_to_cart
   use trc_eri, only: trc_eri_t
   use trc_pcm, only: trc_pcm_t
   use trc_xc, only: trc_xc_rks, trc_xc_uks, trc_xc_options, trc_xc_result_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_vv10, only: trc_vv10_rks, trc_vv10_uks
   use trc_linalg, only: trc_linalg_t
   use trc_scf_df, only: df_jk_t
   use trc_scf_settings, only: scf_settings_t, scf_energies_t, trc_scf_result_t
   use trc_scf_context, only: scf_ctx_t, trc_grid_cache_t
   use trc_scf_state, only: scf_state_t, scf_timers_t
   implicit none
   private

   public :: scf_fock_build_g, scf_fock_build_vxc, scf_fock_build_pcm, scf_fock_energy, scf_fock_assemble, scf_fock_commutator_error
   public :: scf_fock_to_cart, scf_fock_to_sph

contains

   !
   ! The two-electron build, resident and entirely in Cartesian space:
   ! g_rks (nspin = 1, a full build), g_rks_incremental (nspin = 1, from the
   ! density DIFFERENCE, see trc_scf_driver's own header) and g_uks (nspin = 2),
   ! dispatched on `nspin`/`incr_now` rather than three separately-named
   ! entry points, since the caller already has exactly one of them to ask
   ! for.
   !
   subroutine scf_fock_build_g(scf_basis, scf_context, scf_config, scf_state, incr_now, eri, la, dfjk)
      type(trc_basis_t), intent(in) :: scf_basis
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      logical, intent(in) :: incr_now
      !! This iteration builds G from the density difference (the caller's
      !! schedule: `incr_on`, the commutator below `incr_start`, not a reset).
      type(trc_eri_t), intent(inout) :: eri         !! the four-centre engine (unused under `--df`)
      type(trc_linalg_t), intent(inout) :: la
      type(df_jk_t), intent(inout), optional :: dfjk
      !! The fitted J/K object, present exactly when `scf_config%df`.
      call build_g_arrays(scf_basis, eri, scf_context%nao_c, scf_config%nspin, &
                          size(scf_state%dtot_c, 1), size(scf_state%klr_c, 1), size(scf_state%dref_c, 1), &
                          incr_now, scf_config%incr_on, &
                          scf_config%exx, scf_config%exx_lr, scf_config%omega, scf_config%incr_scale_run, &
                          scf_state%d_c, scf_state%gmat_c, scf_state%dtot_c, scf_state%jmat_c, scf_state%w1_c, &
                          scf_state%klr_c, scf_state%gacc_c, scf_state%dref_c, scf_state%ddif_c, &
                          scf_config%df, la, dfjk)
   end subroutine scf_fock_build_g

   subroutine build_g_arrays(b, eri, nao_c, nspin, nu_c, nlr_c, nin_c, incr_now, incr_on, &
                             exx, exx_lr, omega, incr_scale_run, &
                             d_c, gmat_c, dtot_c, jmat_c, w1_c, klr_c, gacc_c, dref_c, ddif_c, &
                             df_on, lp, dfj)
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
      type(trc_linalg_t), intent(inout) :: lp
      type(df_jk_t), intent(inout), optional :: dfj
      !! Present when `df_on`; absent and unused otherwise.
      integer :: i, j, s
      real(dp), allocatable :: w2_c(:, :, :), k2_c(:, :, :)
      !! Unrestricted only: both spins' 2 D_s, and their exchange, for
      !! fock_resident_many.

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
            call copy_into_spin(nao_c, jmat_c, gmat_c(:, :, s))
         end do
         if (exx /= 0.0_dp .or. exx_lr /= 0.0_dp) then
            ! Both spins' exchange in ONE pass over the integrals per operator
            ! (fock_resident_many): every quartet evaluated once and digested
            ! into K[D_alpha] and K[D_beta] together, not in two passes. J is
            ! J[D_total] above (its own pass, or MD J); the spins' passes are
            ! exchange only. A full-range pass only when the functional has a
            ! full-range exact-exchange part at all.
            allocate (w2_c(nao_c, nao_c, 2), k2_c(nao_c, nao_c, 2))
            !$acc enter data create(w2_c, k2_c)
            do s = 1, 2
               call scale_into_spin(nao_c, 2.0_dp, d_c(:, :, s), w2_c(:, :, s))
            end do
            if (exx /= 0.0_dp) then
               call eri%fock_resident_many(b, 2, w2_c, k2_c, j_scale=0.0_dp, k_scale=exx)
               do s = 1, 2
                  call add_into_spin(nao_c, k2_c(:, :, s), gmat_c(:, :, s))
               end do
            end if
            if (exx_lr /= 0.0_dp) then
               call eri%fock_resident_many(b, 2, w2_c, k2_c, j_scale=0.0_dp, k_scale=exx_lr, omega=omega)
               do s = 1, 2
                  call add_into_spin(nao_c, k2_c(:, :, s), gmat_c(:, :, s))
               end do
            end if
            !$acc exit data delete(w2_c, k2_c)
            deallocate (w2_c, k2_c)
         end if
      end if
   end subroutine build_g_arrays

   subroutine copy_into_spin(n, a, c)
      ! c := a, resident (n, n) matrices; explicit shape so a spin slice of
      ! gmat_c can be the target.
      integer, intent(in) :: n
      real(dp), intent(in) :: a(n, n)
      real(dp), intent(out) :: c(n, n)
      integer :: i, j
      do concurrent(i=1:n, j=1:n)
         c(i, j) = a(i, j)
      end do
   end subroutine copy_into_spin

   subroutine scale_into_spin(n, f, a, c)
      ! c := f a, resident.
      integer, intent(in) :: n
      real(dp), intent(in) :: f
      real(dp), intent(in) :: a(n, n)
      real(dp), intent(out) :: c(n, n)
      integer :: i, j
      do concurrent(i=1:n, j=1:n)
         c(i, j) = f*a(i, j)
      end do
   end subroutine scale_into_spin

   subroutine add_into_spin(n, a, c)
      ! c := c + a, resident.
      integer, intent(in) :: n
      real(dp), intent(in) :: a(n, n)
      real(dp), intent(inout) :: c(n, n)
      integer :: i, j
      do concurrent(i=1:n, j=1:n)
         c(i, j) = c(i, j) + a(i, j)
      end do
   end subroutine add_into_spin

   !
   ! The exchange-correlation potential and energy, finding its (Cartesian)
   ! density already present: the semilocal part (trc_xc_rks/uks), then
   ! VV10's nonlocal correlation on its own coarser grid when the functional
   ! has it, folded back to the run's own space. Only called when `dft` is
   ! true (the caller's own guard), so `vxc_c`/`vxc` are always full size
   ! here -- no 1 x 1 sizing dummy needed, unlike `scf_fock_build_g`'s arrays.
   ! The fold back to the run's space reads its sizes off those arrays.
   !
   subroutine scf_fock_build_vxc(scf_basis, grid_cache, la, scf_config, scf_state, scf_energies, scf_timers, comm)
      type(trc_basis_t), intent(in) :: scf_basis
      type(trc_grid_cache_t), intent(in) :: grid_cache   !! the run's grids, built and on the device
      type(trc_linalg_t), intent(inout) :: la
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      type(scf_energies_t), intent(inout) :: scf_energies
      !! `xc` and `nelec_grid` are set here, VV10's energy included.
      type(scf_timers_t), intent(inout) :: scf_timers
      !! `tx_pts`/`tx_prs` and `t_vv10` are added to (not overwritten). The
      !! DFT step table (`xcst`) reports VV10 as its own step, separate from
      !! the rest of the XC build -- see NVTX_RANGES below for why
      !! VV10's is a plain wall time rather than a `step_on`/`step_off` pair.
      type(comm_t), intent(in), optional :: comm
      real(dp) :: e_nlc
      type(trc_xc_result_t) :: xres
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
      ! `trc_xc_uks` could still be running when the "dft: xc" NVTX range
      ! closes, so a profiler would show it finishing inside "dft: vv10"
      ! instead -- the xc row undercounting, the vv10 row overcounting. The
      ! sync measured free here: `trc_xc_rks`/`trc_xc_uks` already force a
      ! device-to-host transfer of `exc`/`nelec` (`xc_drive`'s own `!$acc
      ! update self(eg, ng)`/`update self(vxc)`) before returning, so the
      ! device is already idle by the time `lp%sync()` runs. VV10's own
      ! sync (`tm_vv10`, a plain wall time handed back through `t_vv10` for
      ! `xcst` to accumulate, rather than a `step_table_t` tick) stays for
      ! the same reason.
      !
      call profiler_start("dft: xc", nvtx_only=.true.)
      if (scf_config%nspin == 1) then
         call trc_xc_rks(scf_basis, grid_cache%xg, scf_config%func, scf_state%d_c(:, :, 1), scf_state%vxc_c(:, :, 1), &
                         xres, trc_xc_options(rho_tol=scf_config%rho_tol, comm=comm), la)
      else
         call trc_xc_uks(scf_basis, grid_cache%xg, scf_config%func, scf_state%d_c, scf_state%vxc_c, &
                         xres, trc_xc_options(rho_tol=scf_config%rho_tol, comm=comm), la)
      end if
      scf_energies%xc = xres%exc
      scf_energies%nelec_grid = xres%nelec
      call la%sync()
      call profiler_stop("dft: xc")
      ! VV10 adds its energy to the XC energy and its potential into vxc_c, on
      ! its own (coarser) grid, before the spherical transform below.
      if (scf_config%nlc) then
         call profiler_start("dft: vv10", nvtx_only=.true.)
         call la%sync()
         call tm_vv10%start()
         if (scf_config%nspin == 1) then
            call trc_vv10_rks(scf_basis, grid_cache%xg_nlc, scf_config%func, scf_state%d_c(:, :, 1), e_nlc, &
                              scf_state%vxc_c(:, :, 1), la=la, comm=comm)
         else
            call trc_vv10_uks(scf_basis, grid_cache%xg_nlc, scf_config%func, scf_state%d_c, e_nlc, &
                              scf_state%vxc_c, la=la, comm=comm)
         end if
         call la%sync()
         call profiler_stop("dft: vv10")
         scf_timers%t_vv10 = scf_timers%t_vv10 + tm_vv10%get_elapsed_time()
         scf_energies%xc = scf_energies%xc + e_nlc
      end if
      scf_timers%tx_pts = scf_timers%tx_pts + xres%t_points
      scf_timers%tx_prs = scf_timers%tx_prs + xres%t_pairs
      do s = 1, scf_config%nspin
         call scf_fock_to_sph(scf_basis, scf_basis%spherical, size(scf_state%vxc, 1), size(scf_state%vxc_c, 1), &
                     scf_state%vxc_c(:, :, s), scf_state%vxc(:, :, s))
      end do
   end subroutine scf_fock_build_vxc

   !
   ! The continuum's operator at the TOTAL density (pyscf's choice for UKS
   ! too: one set of charges, one operator added to both spins), folded back
   ! to the run's own space. Only called when `pcm_on` is true. `nu_c` is
   ! `dtot_c`'s true size (1 when `nspin == 1`, where it is never read).
   !
   subroutine scf_fock_build_pcm(scf_basis, scf_context, scf_config, scf_state, scf_energies, la, pl, comm)
      type(trc_basis_t), intent(in) :: scf_basis
      type(scf_ctx_t), intent(inout) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      type(scf_energies_t), intent(inout) :: scf_energies
      !! `pcm` (the dielectric energy) is set here.
      type(trc_linalg_t), intent(inout) :: la
      type(trc_pairlist_t), intent(in) :: pl   !! the run's pair list
      type(comm_t), intent(in), optional :: comm
      ! An absent `comm` passes through as absent: one rank, unsplit.
      if (scf_config%nspin == 1) then
         call scf_context%pcm%operator(la, scf_basis, pl, scf_state%d_c(:, :, 1), scf_state%vpcm_c, &
                                       scf_energies%pcm, comm)
      else
         call scf_context%pcm%operator(la, scf_basis, pl, scf_state%dtot_c, scf_state%vpcm_c, &
                                       scf_energies%pcm, comm)
      end if
      call scf_fock_to_sph(scf_basis, scf_context%sph, scf_context%nao, scf_context%nao_c, scf_state%vpcm_c, scf_state%vpcm)
   end subroutine scf_fock_build_pcm

   ! The one- and two-electron energy at the density `scf_result%dmat`: Tr(D H) and
   ! Tr(D G)/2, per spin; the total adds the nuclear repulsion and the XC and
   ! PCM energies already in `scf_energies`.
   subroutine scf_fock_energy(scf_context, scf_config, scf_state, scf_result, scf_energies, la)
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(in) :: scf_state
      type(trc_scf_result_t), intent(in) :: scf_result
      type(scf_energies_t), intent(inout) :: scf_energies
      type(trc_linalg_t), intent(inout) :: la
      integer :: s
      scf_energies%one = 0.0_dp
      scf_energies%two = 0.0_dp
      do s = 1, scf_config%nspin
         scf_energies%one = scf_energies%one + la%dot_2d(scf_result%dmat(:, :, s), scf_context%hcore)
         scf_energies%two = scf_energies%two + 0.5_dp*la%dot_2d(scf_result%dmat(:, :, s), scf_state%gmat(:, :, s))
      end do
      scf_energies%total = scf_energies%one + scf_energies%two + scf_energies%xc + scf_context%e_nuc + scf_energies%pcm
   end subroutine scf_fock_energy

   subroutine scf_fock_assemble(scf_context, scf_config, scf_state)
      !! F = H + G (+ V_xc, + the continuum's operator), per spin, into
      !! `scf_state%fock`.
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      call assemble_fock_arrays(scf_context%nao, size(scf_state%vxc, 1), size(scf_state%vpcm, 1), scf_config%nspin, &
                                scf_config%dft, scf_config%pcm_on, scf_context%hcore, scf_state%gmat, scf_state%vxc, &
                                scf_state%vpcm, scf_state%fock)
   end subroutine scf_fock_assemble

   subroutine assemble_fock_arrays(nao, nxc, npcm, nspin, dft, pcm_on, hcore, gmat, vxc, vpcm, fock)
      !! `scf_fock_assemble`'s loops, on explicit-shape arrays. `vxc` is
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
   end subroutine assemble_fock_arrays

   ! W = F D S - S D F, errv = X^T W X, per spin: the orthogonalised
   ! commutator DIIS extrapolates on, and whose norm gates DIIS/A-DIIS.
   subroutine scf_fock_commutator_error(scf_context, scf_config, scf_state, scf_result, la)
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      !! Reads `fock`, writes `errv`; `w1`/`w2` are scratch.
      type(trc_scf_result_t), intent(in) :: scf_result
      type(trc_linalg_t), intent(inout) :: la
      integer :: s
      do s = 1, scf_config%nspin
         ! w2 = F D S - S D F
         call la%gemm_mat(scf_state%fock(:, :, s), scf_result%dmat(:, :, s), scf_state%w1)
         call la%gemm_mat(scf_state%w1, scf_context%smat, scf_state%w2)
         call la%gemm_mat(scf_context%smat, scf_result%dmat(:, :, s), scf_state%w1)
         call la%gemm_mat(scf_state%w1, scf_state%fock(:, :, s), scf_state%w2, alpha=-1.0_dp, beta=1.0_dp)
         ! errv = X^T w2 X
         call la%gemm_mat(scf_state%w2, scf_context%x, scf_state%w1)
         call la%gemm_mat(scf_context%x, scf_state%w1, scf_state%errv(:, :, s), transa="T")
      end do
   end subroutine scf_fock_commutator_error

   subroutine scf_fock_to_cart(b, sph, nao, nao_c, a, a_c)
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
   end subroutine scf_fock_to_cart

   subroutine scf_fock_to_sph(b, sph, nao, nao_c, a_c, a)
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
   end subroutine scf_fock_to_sph

end module trc_scf_fock
