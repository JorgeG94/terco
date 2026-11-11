!
! The analytic gradients are validated (pyscf, finite differences, CPU =
! GPU); single rank only. See README.
!
! The full analytic HF/UHF nuclear gradient: one-electron (trc_1e_gradient),
! two-electron (eri%gradient) and nuclear repulsion, assembled at a
! converged trc_scf_result_t. With `func` present (optional, see
! trc_scf_gradient's own header) this is the KS-DFT gradient instead:
! eri%gradient's exact-exchange piece scaled by `func%exx`, plus
! trc_xc_grad%xc_grad_add's LDA/GGA/meta-GGA exchange-correlation term.
!
! ORDER MATTERS. eri%gradient starts `grad` from zero (its own header: "the
! two-electron part ONLY"); trc_1e_gradient and scf_grad_nuclear_repulsion
! both ADD into whatever `grad` already holds. So the assembly below calls
! eri%gradient first and the other two after -- reversing that order would
! have eri%gradient discard the one-electron and nuclear-repulsion terms.
!
! THE FACTOR OF 1/2. eri%gradient's own convention (test/check_rys_grad_fd.F90's
! E2 = Tr(D, eri%fock(D)), no 1/2 -- that is what it is checked against, both
! there and vs pyscf) is the derivative of Tr(D, J(D) - K(D)/2) WITHOUT the
! 1/2 that turns it into the actual two-electron ENERGY. scf_fock_energy's own
! e_two is 0.5*Tr(D, gmat) (trc_scf_fock.F90): half of that same trace. So
! eri%gradient's raw output is exactly TWICE d(e_two)/dR, and has to be
! halved here before it is added to the one-electron and nuclear-repulsion
! terms, which are both already in the "actual energy" convention (verified
! against finite differences of the true SCF energy, not a doubled one).
! Found by a clean 2x, atom-by-atom mismatch against a fixed-density finite
! difference of the real e1 + e2 + e_nuc, not by inspection.
!
! trc_scf_run does not keep its pair list or ERI engine alive past
! convergence (trc_scf_context%release at the end of the run), so
! trc_scf_gradient rebuilds both at the converged geometry -- the same
! thing every eri%gradient caller in test/ already does, and the Schwarz
! screen it needs is the expensive part regardless of whether the gradient
! is wanted.
!
module trc_scf_grad
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_1e_gradient
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_sph, only: sph_sph_to_cart
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_xc_grad, only: xc_grad_add, xc_grad_add_on_grid, xc_grad_step_names, xc_grad_options, xc_grad_options_t
   use trc_step_table, only: step_table_t
   use trc_scf_context, only: trc_pair_eri_cache_t, trc_grid_cache_t
   use trc_decontract, only: trc_decontract_options_t, decontract_options_for_exchange, decontract_options_equal
   use trc_scf_df, only: df_jk_t
   use trc_scf_grad_df, only: scf_grad_df_two_electron
   use trc_vv10_grad, only: vv10_grad_add, vv10_grad_add_on_grid
   use trc_scf_settings, only: DEFAULT_NLC_GRID_LEVEL
   use trc_sym_geom, only: trc_sym_geom_t
   use trc_sym_scf, only: sym_scf_grad_symmetrize
   implicit none
   private
   public :: trc_scf_gradient, scf_grad_nuclear_repulsion, scf_grad_energy_weighted_density
   public :: trc_scf_grad_sym_t, trc_scf_grad_reuse_t

   type :: trc_scf_grad_sym_t
      !! The symmetry the molecule was run in, as `trc_scf_gradient` takes it
      !! (its `sym`). Present, the finished gradient is averaged over the group
      !! with the atom permutation (`sym_scf_grad_symmetrize` of trc_sym_scf),
      !! which zeroes the components the group forbids. The molecule's basis
      !! must sit on `geometry%r_std`, so the gradient is in the standard
      !! frame. The density must be invariant under the group, which is what an
      !! SO-blocked SCF produces. The four-centre term is the full C1 build
      !! either way. Absent, nothing is done.
      type(trc_sym_geom_t) :: geometry   !! the group and its atom permutations (trc_sym_geom)
   end type trc_scf_grad_sym_t

   type :: trc_scf_grad_reuse_t
      !! Which of `trc_scf_gradient`'s `grid_cache` grids the XC terms may run
      !! on without a rebuild (its `grid_reuse`). The caller checks that each
      !! matches the settings the gradient's quadrature has to match the SCF's
      !! on; `grid_cache` alone says only what is built, not what matches.
      logical :: semilocal = .false.   !! `grid_cache%grid`/`xg`, for the semilocal XC term
      logical :: nlc = .false.         !! `grid_cache%grid_nlc`/`xg_nlc`, for the VV10 term
   end type trc_scf_grad_reuse_t

   type :: grad_exchange_t
      !! The exact-exchange weights of the four-centre term, which travel
      !! together: the full-range `exx`, the long-range `exx_lr` and its
      !! range-separation `omega` (1, 0, 0 for Hartree-Fock).
      real(dp) :: exx = 1.0_dp, exx_lr = 0.0_dp, omega = 0.0_dp
   end type grad_exchange_t

contains

   subroutine scf_grad_nuclear_repulsion(b, grad)
      !! dE_nuc/dR, ADDED into grad(3,natm). E_nuc = sum_{A<C} Z_A Z_C /
      !! |R_A - R_C|, and d/dR_A |R_A - R_C|^-1 = -(R_A - R_C)/|R_A-R_C|^3
      !! (the usual inverse-distance derivative), so
      !!
      !!   dE_nuc/dR_A = -sum_{C /= A} Z_A Z_C (R_A - R_C)/|R_A - R_C|^3
      !!
      !! the MINUS sign matters: it is what makes this the gradient of the
      !! energy rather than the Coulomb force on A (force = -gradient).
      !! Host arithmetic on b%at_z/b%at_r, which stay valid plain Fortran
      !! arrays whether or not `b%on_device`.
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(inout) :: grad(3, b%natm)
      integer :: a, c
      real(dp) :: rr
      real(dp) :: dvec(3)
      do a = 1, b%natm
         do c = 1, b%natm
            if (a == c) cycle
            dvec = b%at_r(:, a) - b%at_r(:, c)
            rr = norm2(dvec)
            grad(:, a) = grad(:, a) - b%at_z(a)*b%at_z(c)*dvec/(rr*rr*rr)
         end do
      end do
   end subroutine scf_grad_nuclear_repulsion

   subroutine scf_grad_energy_weighted_density(lp, nao, nspin, cmo, eps, nocc, wmat)
      !! W = sum_s occ_s * sum_{i=1}^{nocc(s)} eps(i,s) c(:,i,s) c(:,i,s)^T,
      !! occ_s = 2 restricted (nspin=1, the one spin channel IS both
      !! electrons) or 1 per channel unrestricted -- the standard
      !! energy-weighted (Lagrangian) density the overlap's Pulay force is
      !! contracted against. Built the same way trc_scf_driver's own
      !! fractional-occupation branch builds a column-scaled density
      !! (orbitals_and_density's `t(i,j) = cmo(i,j)*sqrt(ofrac(j))`), except
      !! the two factors here (C and C*eps) differ, so it is a plain GEMM
      !! rather than linalg_density's single-scalar-times-C-C^T.
      !!
      !! DEVICE RESIDENCY. `lp%gemm_mat` is cuBLAS under `!$acc
      !! host_data use_device(...)` (trc_linalg.F90's linalg_gemm_raw),
      !! which only FETCHES a device pointer for data already present --
      !! it creates none. `cmo`/`eps`/`wmat` must already be present
      !! (the caller's job, same as every other trc_linalg_t call in this
      !! tree); the local scratch `cw` is this routine's own, entered and
      !! exited here.
      type(trc_linalg_t), intent(inout) :: lp
      integer,  intent(in)  :: nao, nspin
      integer,  intent(in)  :: nocc(2)
      real(dp), intent(in)  :: cmo(nao, nao, nspin), eps(nao, nspin)
      real(dp), intent(out) :: wmat(nao, nao)
      real(dp), allocatable :: cw(:, :)
      real(dp) :: occ_s
      integer :: s, i, j, nk

      allocate (cw(nao, nao))
      !$acc enter data create(cw)
      occ_s = merge(2.0_dp, 1.0_dp, nspin == 1)
      do s = 1, nspin
         nk = nocc(s)
         if (nk <= 0) cycle
         do j = 1, nk
            do i = 1, nao
               cw(i, j) = cmo(i, j, s)*eps(j, s)
            end do
         end do
         !$acc update device(cw(:, 1:nk))
         call lp%gemm_mat(cw(:, 1:nk), cmo(:, 1:nk, s), wmat, transb="T", &
                          alpha=occ_s, beta=merge(0.0_dp, 1.0_dp, s == 1))
      end do
      !$acc exit data delete(cw)
      deallocate (cw)
   end subroutine scf_grad_energy_weighted_density

   subroutine trc_scf_gradient(b, lp, nspin, nalpha, nbeta, cmo, eps, thresh, grad, &
                               func, grid_level, grid_max_pts, grid_response, rho_tol, nlc_grid_level, grid_rotation, &
                               dfjk, aux, verbose, grid_cache, grid_reuse, pair_eri_cache, reused, sym)
      !! The full RHF/UHF analytic gradient at the converged orbitals
      !! `cmo`/`eps`: eri%gradient (two-electron) + trc_1e_gradient
      !! (overlap Pulay term, kinetic, nuclear attraction) +
      !! scf_grad_nuclear_repulsion. Takes orbitals and electron counts
      !! rather than a trc_scf_result_t so a C API caller that only kept
      !! `cmo`/`eps` (trc_capi.F90's box handle, after trc_run_scf moved
      !! them out of its own local `res`) can call this directly, with no
      !! second copy of the density/aufbau-occupation bookkeeping: D is
      !! rebuilt here from cmo/eps/nocc exactly the way trc_scf_driver
      !! itself builds it (lp%density, occ=2 restricted or 1 per channel).
      !!
      !! SPHERICAL BASES. `cmo`/`eps` are sized at whatever space the SCF
      !! itself ran in -- `b%nao_sph` when `b%spherical`, `b%nao`
      !! (Cartesian) otherwise, same `merge` trc_capi.F90's own
      !! `scf_space_nao` uses -- but every kernel below (eri%gradient,
      !! trc_1e_gradient) is Cartesian throughout (trc_scf_driver.F90's own
      !! "SPHERICAL BASES" header: "Everything touching a kernel stays
      !! Cartesian"). D and W are therefore built in the SCF's own space
      !! first, then folded to Cartesian with `sph_sph_to_cart` -- T D T^T,
      !! the same transform trc_scf_driver folds the Fock-build density on
      !! each iteration -- before anything here touches an integral. This
      !! is EXACT, not an approximation: contracting the Cartesian
      !! derivative integrals with T D T^T equals contracting the
      !! spherical ones with D directly, because T is exactly the
      !! (per-shell) change of basis the SCF's own spherical AOs are
      !! defined by. When the basis is not spherical, `b%nao_sph` still
      !! exists structurally (every shell's would-be spherical count) but
      !! plays no role here: the fold is skipped and `density`/`wmat` are
      !! filled directly at `b%nao`, unchanged from before this feature.
      !!
      !! `func` ABSENT is the Hartree-Fock gradient, unchanged. PRESENT
      !! adds the exchange-correlation term (trc_xc_grad%xc_grad_add,
      !! LDA/GGA only -- see its own refusals) and scales eri%gradient's
      !! exact-exchange piece by the functional's `exx`, exactly as
      !! trc_scf_fock.F90 scales the Fock build's K by the same field.
      !!
      !! RANGE-SEPARATED HYBRIDS: `func%exx_lr`/`func%omega`, same meaning
      !! as trc_scf_fock's scf_fock_build_g. When `func%exx_lr` is nonzero, a SECOND
      !! eri%gradient pass builds the long-range exchange gradient
      !! (j_scale = 0, k_scale = exx_lr, omega) and adds it on top of the
      !! full-range one, next to it -- exactly how scf_fock_build_g assembles the
      !! energy (eri%fock_resident with j_scale=0, k_scale=exx_lr,
      !! omega, added into gmat_c). This routine computes that term for
      !! ANY func with exx_lr /= 0, and so does xc_grad_add's
      !! semilocal GGA term below (func%omega/exx_lr do not reach it --
      !! see that module's own header) -- an RSH functional like wB97X is
      !! NOT refused anywhere in this path any more.
      !!
      !! VV10 (func%nlc_b/nlc_c nonzero, e.g. wB97X-V) is NOT refused
      !! either any more: trc_vv10_grad%vv10_grad_add adds the
      !! nonlocal correlation term on top, on its own grid (`nlc_grid_level`,
      !! below) -- single rank only, see that module's own header.
      type(trc_basis_t),  intent(inout) :: b
      type(trc_linalg_t), intent(inout) :: lp
      integer,            intent(in)    :: nspin, nalpha, nbeta
      real(dp),           intent(in)    :: cmo(merge(b%nao_sph, b%nao, b%spherical), &
                                               merge(b%nao_sph, b%nao, b%spherical), nspin)
      real(dp),           intent(in)    :: eps(merge(b%nao_sph, b%nao, b%spherical), nspin)
      real(dp),           intent(in)    :: thresh
      real(dp),           intent(out)   :: grad(3, b%natm)
      type(trc_xc_functional_t), intent(in), optional :: func
      integer, intent(in), optional :: grid_level, grid_max_pts
      logical, intent(in), optional :: grid_response
      !! Forwarded to xc_grad_add verbatim (its own default: ON). See
      !! that module's header, GRID RESPONSE.
      real(dp), intent(in), optional :: rho_tol
      !! The SAME grid level/max-points-per-batch/rho screen the SCF itself
      !! converged with (trc_capi.F90's trc_run_gradient passes cx%opts's
      !! own grid_level/grid_max_pts/rho_tol) -- the gradient's quadrature
      !! has to match the energy's, not merely be "a reasonable DFT grid",
      !! or the two are not derivatives of the same number. Absent (a
      !! Fortran-side caller with no opinion) falls back to xc_grad_add's
      !! own defaults.
      type(trc_pair_eri_cache_t), intent(inout), optional :: pair_eri_cache
      !! The SCF's own pair list and (outside `--df`) four-centre engine,
      !! kept alive since `trc_scf_run` (trc_scf_context.F90's
      !! `trc_pair_eri_cache_t`) -- reused INSTEAD of a fresh build when
      !! `pair_eri_cache%matches(b, thresh)`: the build threshold EXACTLY
      !! (both sides read the same deterministic `scf_settings_eri_thresh_env`)
      !! and the cache's own geometry/basis fingerprint, the backstop a
      !! Fortran caller needs (the C API's `cx%scf_ok` covers geometry for
      !! that caller). The engine is reused only when, on top of that,
      !! `batch_res_built` is the default (1), the batching this routine's
      !! own fresh build uses, or "with cache" and "without" would differ by
      !! more than round-off. A mismatch, or no cache, rebuilds both: a pure
      !! performance path, never a correctness one, same as the grid's.
      logical, intent(out), optional :: reused
      !! Set to whether this call actually took the reuse path (the pair
      !! list, and -- outside `--df` -- the four-centre engine too) rather
      !! than rebuilding: the only way a caller or a test can tell the two
      !! paths apart, since they agree on the numbers by construction.
      integer, intent(in), optional :: nlc_grid_level
      real(dp), intent(in), optional :: grid_rotation(3, 3)
      !! The Lebedev rotation of the SCF's grids (`trc_scf_options_t%grid_rotation`),
      !! for a grid this routine has to build itself; absent: the identity.
      !! VV10's own grid level -- see this subroutine's own header, VV10.
      !! Absent means `trc_scf_settings%DEFAULT_NLC_GRID_LEVEL`, the SAME
      !! constant `trc_scf_options_t%nlc_grid_level` defaults to, not the
      !! semilocal grid's: the gradient has to be taken on the grid the
      !! energy used, and a caller's own default must not drift from the
      !! SCF's.
      type(df_jk_t), intent(in), optional :: dfjk
      !! DENSITY FITTING (`--df`): present (with `aux`), the two-electron
      !! term is the FITTED one, trc_scf_grad_df's scf_grad_df_two_electron
      !! over the SAME fit (`dfjk`'s retained B and metric Cholesky factor)
      !! the SCF converged with, instead of eri%gradient's four-centre one --
      !! no ERI engine is built at all. Restricted and Cartesian only, the
      !! DF-SCF's own scope, and no range separation (the DF-SCF fits no
      !! erf exchange): each is an error stop here, refused up front by
      !! trc_capi.F90's trc_run_gradient.
      type(trc_basis_t), intent(inout), optional :: aux
      !! The auxiliary (jkfit) basis `dfjk` was fitted with.
      logical, intent(in), optional :: verbose
      !! Log the gradient's timing table, its phases and the XC term's
      !! steps, through pic's logger at PERFORMANCE. Default off.
      type(trc_grid_cache_t), intent(in), optional :: grid_cache
      type(trc_scf_grad_reuse_t), intent(in), optional :: grid_reuse
      !! The SCF's own semilocal grid (`grid_cache%grid`/`xg`) and, for VV10, its
      !! nlc grid (`grid_nlc`/`xg_nlc`), and which of them the caller has
      !! checked to match (`grid_reuse`): each so marked is used by its XC term
      !! (`xc_grad_add_on_grid`/`vv10_grad_add_on_grid`) with no rebuild; any
      !! other, or no cache, is built by the term itself. Both arguments are
      !! needed for a reuse. trc_capi.F90's `trc_run_gradient` is the one caller.
      type(trc_scf_grad_sym_t), intent(in), optional :: sym   !! see its type; absent, the C1 gradient
      type(trc_pairlist_t) :: pl_own
      !! The pair list built here when the cache's is not reused.
      type(step_table_t) :: gst, xst, vst
      type(xc_grad_options_t) :: xopts
      type(timer_type) :: tc
      logical :: talk, reuse_pl, reuse_eri, reuse_semilocal, reuse_nlc
      integer :: nlc_level
      integer :: nao
      integer :: nocc(2)
      type(grad_exchange_t) :: xch
      real(dp), allocatable :: wmat(:, :), dtot(:, :), density(:, :, :)

      talk = .false.
      if (present(verbose)) talk = verbose
      reuse_semilocal = .false.
      reuse_nlc = .false.
      if (present(grid_cache) .and. present(grid_reuse)) then
         reuse_semilocal = grid_reuse%semilocal
         reuse_nlc = grid_reuse%nlc
      end if
      call gst%init("grad", [character(len=12) :: "density", "two-electron", "one-electron", "xc", "vv10"])
      call xst%init("xc grad", xc_grad_step_names())
      call vst%init("vv10 grad", xc_grad_step_names())
      call gst%step_on(lp, 1, tc)

      nao = b%nao
      nocc = [nalpha, nbeta]
      if (present(func)) then
         xch%exx = func%exx
         xch%exx_lr = func%exx_lr
         xch%omega = func%omega
      end if

      ! DEVICE RESIDENCY, for the whole routine: `cmo`/`eps` arrive as
      ! plain host arrays (trc_scf_run's own `release_device` -- see
      ! trc_scf_driver.F90 -- takes them OFF the device before returning,
      ! so neither a Fortran-side caller's `res%cmo`/`res%eps` nor the C
      ! API's box handle's copy of them is resident); `density`/`wmat` are
      ! this routine's own, at the Cartesian `nao` the kernels below need
      ! -- filled directly there when the basis is Cartesian, or via the
      ! spherical fold (its own block, right below) otherwise.
      allocate (density(nao, nao, nspin), wmat(nao, nao), dtot(nao, nao))
      !$acc enter data copyin(cmo, eps) create(density, wmat)

      ! D, W (in the SCF's own space, folded to Cartesian if spherical), and D_total.
      call scf_grad_density_weighted(b, lp, nocc, cmo, eps, density, wmat, dtot)

      call gst%step_off(lp, 1, tc)
      call gst%step_on(lp, 2, tc)
      call scf_grad_reuse_decide(b, thresh, xch, pair_eri_cache, reuse_pl, reuse_eri)
      if (present(reused)) reused = reuse_pl .and. (reuse_eri .or. present(dfjk))
      if (.not. reuse_pl) then
         call pl_own%build(b, thresh)
         call pl_own%to_device()
      end if

      if (present(dfjk)) then
         if (.not. present(aux)) error stop "trc_scf_gradient: dfjk given without its auxiliary basis"
         if (nspin /= 1 .or. b%spherical .or. xch%exx_lr /= 0.0_dp) then
            error stop "trc_scf_gradient: the --df gradient is restricted, Cartesian and not range-separated"
         end if
         ! --- density-fitted two-electron term, in the actual-energy
         ! convention (no halving: scf_grad_df_two_electron differentiates
         ! E_2 itself, not eri%gradient's doubled trace) ---
         call aux%to_device()
         grad = 0.0_dp
         !$acc enter data copyin(dtot, grad)
         if (reuse_pl) then
            call scf_grad_df_two_electron(b, aux, pair_eri_cache%pl, dfjk, lp, nalpha, cmo(:, :, 1), dtot, 1.0_dp, &
                                          xch%exx, grad)
         else
            call scf_grad_df_two_electron(b, aux, pl_own, dfjk, lp, nalpha, cmo(:, :, 1), dtot, 1.0_dp, xch%exx, grad)
         end if
      else
         if (reuse_eri) then
            call two_electron_4c(b, thresh, density, xch, grad, eri_cache=pair_eri_cache%eri)
         else
            call two_electron_fresh(b, thresh, density, xch, grad, pair_eri_cache)
         end if
         ! --- one-electron, added on top. `wmat` is already present (and
         ! already holds the right numbers, from scf_grad_energy_weighted_density
         ! above) -- only `dtot` and `grad` are freshly copied up. ---
         !$acc enter data copyin(dtot, grad)
      end if
      call gst%step_off(lp, 2, tc)
      call gst%step_on(lp, 3, tc)
      if (reuse_pl) then
         call trc_1e_gradient(b, pair_eri_cache%pl, wmat, dtot, grad)
      else
         call trc_1e_gradient(b, pl_own, wmat, dtot, grad)
      end if
      call gst%step_off(lp, 3, tc)

      ! --- exchange-correlation, added on top (LDA, GGA, meta-GGA). `density` is
      ! still present and correct on the device at this point (see this
      ! routine's own DEVICE RESIDENCY note above); `grad` likewise. ---
      if (present(func)) then
         call gst%step_on(lp, 4, tc)
         xopts = xc_grad_options(grid_level=grid_level, grid_max_pts=grid_max_pts, grid_response=grid_response, &
                                 rho_tol=rho_tol, grid_rotation=grid_rotation)
         if (reuse_semilocal) then
            xopts%partition = grid_cache%grid%partition
            call xc_grad_add_on_grid(b, lp, func, density, grad, xopts, xst, grid_cache%xg)
         else
            call xc_grad_add(b, lp, func, nspin, density, grad, xopts, xst)
         end if
         call gst%step_off(lp, 4, tc)
         ! --- VV10 nonlocal correlation, added on top. Its own grid, its
         ! own screen -- see trc_vv10_grad's own header. `density` is still
         ! present and correct on the device (same note as xc_grad_add
         ! above); `grad` likewise. ---
         if (func%nlc_b /= 0.0_dp .or. func%nlc_c /= 0.0_dp) then
            nlc_level = DEFAULT_NLC_GRID_LEVEL
            if (present(nlc_grid_level)) nlc_level = nlc_grid_level
            ! vv10_grad_add has no `rho_tol` of its own: it always
            ! screens at trc_vv10's own VV10_RHO_TOL, the same constant
            ! the energy path's `trc_vv10_rks/uks` is hard-coded to --
            ! see that module's own header.
            call gst%step_on(lp, 5, tc)
            xopts = xc_grad_options(grid_level=nlc_level, grid_max_pts=grid_max_pts, grid_response=grid_response, &
                                    grid_rotation=grid_rotation)
            if (reuse_nlc) then
               xopts%partition = grid_cache%grid_nlc%partition
               call vv10_grad_add_on_grid(b, lp, func, density, grad, xopts, vst, grid_cache%xg_nlc)
            else
               call vv10_grad_add(b, lp, func, nspin, density, grad, xopts, vst)
            end if
            call gst%step_off(lp, 5, tc)
         end if
      end if

      !$acc update self(grad)
      !$acc exit data delete(wmat, dtot, grad, cmo, eps, density)

      ! --- nuclear repulsion, host arithmetic, added on top ---
      call scf_grad_nuclear_repulsion(b, grad)
      if (present(sym)) call sym_scf_grad_symmetrize(sym%geometry%group, sym%geometry%perm, grad)

      ! A reused `pair_eri_cache%pl` is the caller's: only `pl_own` is ours.
      if (.not. reuse_pl) call pl_own%release()
      deallocate (wmat, dtot, density)

      if (talk) then
         call lp%sync()
         call gst%report(sum(gst%t), 0, 1)
         if (present(func)) call xst%report(gst%t(4), 0, 1)
         if (gst%t(5) > 0.0_dp) call vst%report(gst%t(5), 0, 1)
      end if
   end subroutine trc_scf_gradient

   subroutine scf_grad_density_weighted(b, lp, nocc, cmo, eps, density_cart, wmat_cart, dtot)
      !! The converged density D (per spin), the energy-weighted density W and
      !! D_total, all Cartesian on the host and, for `density_cart`/`wmat_cart`,
      !! resident on the device afterwards, from the orbitals `cmo`/`eps`
      !! (present on the device already). A spherical basis builds D and W in
      !! the SCF's own space first and folds them to Cartesian with
      !! `sph_sph_to_cart` (trc_scf_gradient's header, SPHERICAL BASES);
      !! `density_cart`/`wmat_cart` are created on the device by the caller.
      type(trc_basis_t), intent(inout) :: b
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nocc(2)
      real(dp), contiguous, intent(in) :: cmo(:, :, :), eps(:, :)
      real(dp), contiguous, intent(inout) :: density_cart(:, :, :)
      real(dp), contiguous, intent(inout) :: wmat_cart(:, :)
      real(dp), intent(out) :: dtot(:, :)
      integer :: s, nspin, nsp
      real(dp) :: occ_s
      nspin = size(cmo, 3)
      nsp = size(cmo, 1)
      occ_s = merge(2.0_dp, 1.0_dp, nspin == 1)
      if (.not. b%spherical) then
         do s = 1, nspin
            call lp%density(cmo(:, :, s), nocc(s), occ_s, density_cart(:, :, s))
         end do
         call scf_grad_energy_weighted_density(lp, nsp, nspin, cmo, eps, nocc, wmat_cart)
      else
         ! Build D and W in the SCF's own (spherical) space, nsp = nao_sph,
         ! then fold each one to Cartesian -- see this subroutine's own
         ! header. `density_s`/`wmat_s` are local to this block, entered and
         ! exited here; `density_cart`/`wmat_cart` (already created on the
         ! device above) receive the fold's output and are what every line
         ! after this block reads, exactly as if the basis had been Cartesian
         ! all along.
         block
            real(dp), allocatable :: density_s(:, :, :), wmat_s(:, :)
            allocate (density_s(nsp, nsp, nspin), wmat_s(nsp, nsp))
            !$acc enter data create(density_s, wmat_s)
            do s = 1, nspin
               call lp%density(cmo(:, :, s), nocc(s), occ_s, density_s(:, :, s))
            end do
            call scf_grad_energy_weighted_density(lp, nsp, nspin, cmo, eps, nocc, wmat_s)
            do s = 1, nspin
               call sph_sph_to_cart(b, density_s(:, :, s), density_cart(:, :, s))
            end do
            call sph_sph_to_cart(b, wmat_s, wmat_cart)
            !$acc exit data delete(density_s, wmat_s)
            deallocate (density_s, wmat_s)
         end block
      end if

      ! `density` only has valid numbers on the device at this point; bring
      ! it back before any HOST arithmetic reads it (dtot below) or hands
      ! it to eri%gradient, which re-copies it up itself as a plain host
      ! array (trc_eri.F90's eri_gradient, "density_cart(:, :, :), intent(in)").
      !$acc update self(density_cart)
      if (nspin == 1) then
         dtot = density_cart(:, :, 1)
      else
         dtot = density_cart(:, :, 1) + density_cart(:, :, 2)
      end if

   end subroutine scf_grad_density_weighted

   subroutine scf_grad_reuse_decide(b, thresh, xch, pair_eri_cache, reuse_pl, reuse_eri)
      !! Whether the SCF's own pair list (`reuse_pl`) and four-centre engine
      !! (`reuse_eri`) in `pair_eri_cache` can serve this gradient instead of
      !! a fresh build. A pure performance decision, never a correctness one.
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      type(grad_exchange_t), intent(in) :: xch
      type(trc_pair_eri_cache_t), intent(in), optional :: pair_eri_cache
      logical, intent(out) :: reuse_pl, reuse_eri
      ! REUSE, see `pair_eri_cache`'s own header above: `matches` is the
      ! threshold-and-geometry gate (bit-identical thresholds, exact
      ! fingerprint), then batch_res 1 and the split policy for the engine:
      ! a fresh build (two_electron_fresh) never batches and splits as
      ! `gradient_policy` says, and only that build is reproduced exactly.
      reuse_pl = present(pair_eri_cache)
      if (reuse_pl) reuse_pl = pair_eri_cache%matches(b, thresh)
      reuse_eri = .false.
      ! The engine must also have been split by the policy this gradient
      ! would build with.
      if (reuse_pl) reuse_eri = pair_eri_cache%has_eri .and. pair_eri_cache%batch_res_built == 1
      if (reuse_eri) reuse_eri = decontract_options_equal(pair_eri_cache%decontract_built, &
                                                          gradient_policy(pair_eri_cache, xch%exx, xch%exx_lr))
   end subroutine scf_grad_reuse_decide

   function gradient_policy(pair_eri_cache, exx, exx_lr) result(policy)
      !! The split policy a gradient's own engine is built with: the one the
      !! SCF was run with (`trc_scf_options_t%decontract`, carried by its
      !! cache whether or not the cache can be reused here), resolved for
      !! this method's exchange as the SCF resolved it
      !! (`decontract_options_for_exchange`). The defaults without a cache.
      type(trc_pair_eri_cache_t), intent(in), optional :: pair_eri_cache
      real(dp), intent(in) :: exx, exx_lr
      type(trc_decontract_options_t) :: policy
      policy = trc_decontract_options_t()
      if (present(pair_eri_cache)) policy = pair_eri_cache%decontract_base
      policy = decontract_options_for_exchange(policy, exx, exx_lr)
   end function gradient_policy

   subroutine two_electron_fresh(b, thresh, density, xch, grad, pair_eri_cache)
      !! two_electron_4c over an engine built here, with `gradient_policy`.
      type(trc_basis_t), intent(inout) :: b
      real(dp), intent(in) :: thresh
      real(dp), intent(in) :: density(:, :, :)
      type(grad_exchange_t), intent(in) :: xch
      real(dp), intent(out) :: grad(3, b%natm)
      type(trc_pair_eri_cache_t), intent(in), optional :: pair_eri_cache
      type(trc_eri_t) :: eri_own
      call eri_own%build(b, thresh, decontract=gradient_policy(pair_eri_cache, xch%exx, xch%exx_lr))
      call two_electron_4c(b, thresh, density, xch, grad, eri_own)
      call eri_own%release()
   end subroutine two_electron_fresh

   recursive subroutine two_electron_4c(b, thresh, density, xch, grad, eri_cache)
      !! The four-centre two-electron gradient (eri%gradient), halved into
      !! the actual-energy convention -- trc_scf_gradient's non-`--df`
      !! branch. `grad` is OVERWRITTEN (host).
      type(trc_basis_t), intent(inout) :: b
      real(dp), intent(in) :: thresh
      real(dp), intent(in) :: density(:, :, :)
      type(grad_exchange_t), intent(in) :: xch   !! the exchange weights of the full-range and long-range passes
      real(dp), intent(out) :: grad(3, b%natm)
      type(trc_eri_t), intent(inout), optional :: eri_cache
      !! The SCF's own four-centre engine at this SAME geometry and
      !! threshold (trc_scf_gradient's own `pair_eri_cache`, checked
      !! there before this is ever passed) -- REUSED here instead of a
      !! fresh build: the Schwarz screen and quartet bins the Fock build
      !! already paid for cost nothing extra to run the gradient kernels
      !! over, since `%gradient` takes `omega`/`j_scale`/`k_scale` per
      !! call and reads nothing from `build` that depends on them.
      !! Absent, this routine builds its own and releases it (the
      !! recursive call below, with one).
      real(dp), allocatable :: grad_lr(:, :)

      if (.not. present(eri_cache)) then
         block
            type(trc_eri_t) :: eri_own
            call eri_own%build(b, thresh, decontract=gradient_policy(exx=xch%exx, exx_lr=xch%exx_lr))
            call two_electron_4c(b, thresh, density, xch, grad, eri_own)
            call eri_own%release()
         end block
         return
      end if
      ! --- two-electron, first: this call starts `grad` from zero ---
      call eri_cache%gradient(b, density, grad, k_scale=xch%exx)

      ! --- the long-range exchange gradient, for a range-separated
      ! hybrid: a SECOND eri%gradient pass (j_scale=0, k_scale=exx_lr,
      ! omega), added into `grad` next to the full-range term above --
      ! exactly how scf_fock_build_g assembles the energy (trc_scf_fock.F90's own
      ! j_scale=0.0_dp, k_scale=exx_lr, omega=omega fock_resident call,
      ! added into gmat_c). eri%gradient always starts ITS OWN `grad_lr`
      ! from zero (this module's own header), so the two passes are added
      ! here rather than accumulated into one eri%gradient call. Skipped
      ! entirely when exx_lr is zero (the default): no second build, no
      ! second allocation, same cost as before this argument existed. ---
      if (xch%exx_lr /= 0.0_dp) then
         allocate (grad_lr(3, b%natm))
         call eri_cache%gradient(b, density, grad_lr, j_scale=0.0_dp, k_scale=xch%exx_lr, omega=xch%omega)
         grad = grad + grad_lr
         deallocate (grad_lr)
      end if
      grad = 0.5_dp*grad   ! see this module's own header: eri%gradient's
                           ! convention is twice d(e_two)/dR -- true of
                           ! both passes above, added before this halving.
   end subroutine two_electron_4c

end module trc_scf_grad
