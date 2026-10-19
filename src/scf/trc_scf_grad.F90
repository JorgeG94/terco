!
! The full analytic HF/UHF nuclear gradient: one-electron (trc_1e_gradient),
! two-electron (eri%gradient) and nuclear repulsion, assembled at a
! converged trc_scf_result_t. With `func` present (optional, see
! trc_scf_gradient's own header) this is the KS-DFT gradient instead:
! eri%gradient's exact-exchange piece scaled by `func%exx`, plus
! trc_xc_grad%xc_gradient_add's LDA/GGA exchange-correlation term.
!
! ORDER MATTERS. eri%gradient starts `grad` from zero (its own header: "the
! two-electron part ONLY"); trc_1e_gradient and nuclear_repulsion_gradient
! both ADD into whatever `grad` already holds. So the assembly below calls
! eri%gradient first and the other two after -- reversing that order would
! have eri%gradient discard the one-electron and nuclear-repulsion terms.
!
! THE FACTOR OF 1/2. eri%gradient's own convention (test/check_rys_grad_fd.F90's
! E2 = Tr(D, eri%fock(D)), no 1/2 -- that is what it is checked against, both
! there and vs pyscf) is the derivative of Tr(D, J(D) - K(D)/2) WITHOUT the
! 1/2 that turns it into the actual two-electron ENERGY. scf_energy's own
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
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e_gradient
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_sph, only: sph_sph_to_cart
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_xc_grad, only: xc_gradient_add
   implicit none
   private
   public :: trc_scf_gradient, nuclear_repulsion_gradient, energy_weighted_density

contains

   subroutine nuclear_repulsion_gradient(b, grad)
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
   end subroutine nuclear_repulsion_gradient

   subroutine energy_weighted_density(lp, nao, nspin, cmo, eps, nocc, wmat)
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
   end subroutine energy_weighted_density

   subroutine trc_scf_gradient(b, lp, nspin, nalpha, nbeta, cmo, eps, thresh, grad, &
                               func, grid_level, grid_max_pts, grid_response, rho_tol)
      !! The full RHF/UHF analytic gradient at the converged orbitals
      !! `cmo`/`eps`: eri%gradient (two-electron) + trc_1e_gradient
      !! (overlap Pulay term, kinetic, nuclear attraction) +
      !! nuclear_repulsion_gradient. Takes orbitals and electron counts
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
      !! plays no role here: the fold is skipped and `dens`/`wmat` are
      !! filled directly at `b%nao`, unchanged from before this feature.
      !!
      !! `func` ABSENT is the Hartree-Fock gradient, unchanged. PRESENT
      !! adds the exchange-correlation term (trc_xc_grad%xc_gradient_add,
      !! LDA/GGA only -- see its own refusals) and scales eri%gradient's
      !! exact-exchange piece by the functional's `exx`, exactly as
      !! trc_scf_fock.F90 scales the Fock build's K by the same field.
      !!
      !! RANGE-SEPARATED HYBRIDS: `func%exx_lr`/`func%omega`, same meaning
      !! as trc_scf_fock's build_g. When `func%exx_lr` is nonzero, a SECOND
      !! eri%gradient pass builds the long-range exchange gradient
      !! (j_scale = 0, k_scale = exx_lr, omega) and adds it on top of the
      !! full-range one, next to it -- exactly how build_g assembles the
      !! energy (eri%fock_resident with j_scale=0, k_scale=exx_lr,
      !! omega, added into gmat_c). This routine computes that term for
      !! ANY func with exx_lr /= 0, and so does xc_gradient_add's
      !! semilocal GGA term below (func%omega/exx_lr do not reach it --
      !! see that module's own header) -- an RSH functional like wB97X is
      !! NOT refused anywhere in this path any more. Only VV10
      !! (func%nlc_b/nlc_c) still is, in xc_gradient_add and in
      !! trc_capi.F90's trc_run_gradient.
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
      !! Forwarded to xc_gradient_add verbatim (its own default: ON). See
      !! that module's header, GRID RESPONSE.
      real(dp), intent(in), optional :: rho_tol
      !! The SAME grid level/max-points-per-batch/rho screen the SCF itself
      !! converged with (trc_capi.F90's trc_run_gradient passes cx%opts's
      !! own grid_level/grid_max_pts/rho_tol) -- the gradient's quadrature
      !! has to match the energy's, not merely be "a reasonable DFT grid",
      !! or the two are not derivatives of the same number. Absent (a
      !! Fortran-side caller with no opinion) falls back to xc_gradient_add's
      !! own defaults.

      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      integer :: nao, nsp, s
      integer :: nocc(2)
      real(dp) :: occ_s, exx_use, exx_lr_use, omega_use
      real(dp), allocatable :: wmat(:, :), dtot(:, :), dens(:, :, :), grad_lr(:, :)

      nao = b%nao
      nsp = merge(b%nao_sph, b%nao, b%spherical)
      nocc = [nalpha, nbeta]
      occ_s = merge(2.0_dp, 1.0_dp, nspin == 1)
      exx_use = 1.0_dp
      exx_lr_use = 0.0_dp
      omega_use = 0.0_dp
      if (present(func)) then
         exx_use = func%exx
         exx_lr_use = func%exx_lr
         omega_use = func%omega
      end if

      ! DEVICE RESIDENCY, for the whole routine: `cmo`/`eps` arrive as
      ! plain host arrays (trc_scf_run's own `release_device` -- see
      ! trc_scf_driver.F90 -- takes them OFF the device before returning,
      ! so neither a Fortran-side caller's `res%cmo`/`res%eps` nor the C
      ! API's box handle's copy of them is resident); `dens`/`wmat` are
      ! this routine's own, at the Cartesian `nao` the kernels below need
      ! -- filled directly there when the basis is Cartesian, or via the
      ! spherical fold (its own block, right below) otherwise.
      allocate (dens(nao, nao, nspin), wmat(nao, nao), dtot(nao, nao))
      !$acc enter data copyin(cmo, eps) create(dens, wmat)

      if (.not. b%spherical) then
         do s = 1, nspin
            call lp%density(cmo(:, :, s), nocc(s), occ_s, dens(:, :, s))
         end do
         call energy_weighted_density(lp, nao, nspin, cmo, eps, nocc, wmat)
      else
         ! Build D and W in the SCF's own (spherical) space, nsp = nao_sph,
         ! then fold each one to Cartesian -- see this subroutine's own
         ! header. `dens_s`/`wmat_s` are local to this block, entered and
         ! exited here; `dens`/`wmat` (already created on the device above)
         ! receive the fold's output and are what every line after this
         ! block reads, exactly as if the basis had been Cartesian all
         ! along.
         block
            real(dp), allocatable :: dens_s(:, :, :), wmat_s(:, :)
            allocate (dens_s(nsp, nsp, nspin), wmat_s(nsp, nsp))
            !$acc enter data create(dens_s, wmat_s)
            do s = 1, nspin
               call lp%density(cmo(:, :, s), nocc(s), occ_s, dens_s(:, :, s))
            end do
            call energy_weighted_density(lp, nsp, nspin, cmo, eps, nocc, wmat_s)
            do s = 1, nspin
               call sph_sph_to_cart(b, dens_s(:, :, s), dens(:, :, s))
            end do
            call sph_sph_to_cart(b, wmat_s, wmat)
            !$acc exit data delete(dens_s, wmat_s)
            deallocate (dens_s, wmat_s)
         end block
      end if

      ! `dens` only has valid numbers on the device at this point; bring
      ! it back before any HOST arithmetic reads it (dtot below) or hands
      ! it to eri%gradient, which re-copies it up itself as a plain host
      ! array (trc_eri.F90's eri_gradient, "dens(:, :, :), intent(in)").
      !$acc update self(dens)
      if (nspin == 1) then
         dtot = dens(:, :, 1)
      else
         dtot = dens(:, :, 1) + dens(:, :, 2)
      end if

      call pl%build(b, thresh)
      call pl%to_device()
      call eri%build(b, thresh)

      ! --- two-electron, first: this call starts `grad` from zero ---
      call eri%gradient(b, dens, grad, k_scale=exx_use)

      ! --- the long-range exchange gradient, for a range-separated
      ! hybrid: a SECOND eri%gradient pass (j_scale=0, k_scale=exx_lr,
      ! omega), added into `grad` next to the full-range term above --
      ! exactly how build_g assembles the energy (trc_scf_fock.F90's own
      ! j_scale=0.0_dp, k_scale=exx_lr, omega=omega fock_resident call,
      ! added into gmat_c). eri%gradient always starts ITS OWN `grad_lr`
      ! from zero (this module's own header), so the two passes are added
      ! here rather than accumulated into one eri%gradient call. Skipped
      ! entirely when exx_lr is zero (the default): no second build, no
      ! second allocation, same cost as before this argument existed. ---
      if (exx_lr_use /= 0.0_dp) then
         allocate (grad_lr(3, b%natm))
         call eri%gradient(b, dens, grad_lr, j_scale=0.0_dp, k_scale=exx_lr_use, omega=omega_use)
         grad = grad + grad_lr
         deallocate (grad_lr)
      end if
      grad = 0.5_dp*grad   ! see this module's own header: eri%gradient's
                           ! convention is twice d(e_two)/dR -- true of
                           ! both passes above, added before this halving.

      ! --- one-electron, added on top. `wmat` is already present (and
      ! already holds the right numbers, from energy_weighted_density
      ! above) -- only `dtot` and `grad` are freshly copied up. ---
      !$acc enter data copyin(dtot, grad)
      call trc_1e_gradient(b, pl, wmat, dtot, grad)

      ! --- exchange-correlation, added on top, LDA/GGA only. `dens` is
      ! still present and correct on the device at this point (see this
      ! routine's own DEVICE RESIDENCY note above); `grad` likewise. ---
      if (present(func)) call xc_gradient_add(b, lp, func, nspin, dens, grad, &
                                              grid_level=grid_level, grid_max_pts=grid_max_pts, &
                                              grid_response=grid_response, rho_tol=rho_tol)

      !$acc update self(grad)
      !$acc exit data delete(wmat, dtot, grad, cmo, eps, dens)

      ! --- nuclear repulsion, host arithmetic, added on top ---
      call nuclear_repulsion_gradient(b, grad)

      call eri%release()
      call pl%release()
      deallocate (wmat, dtot, dens)
   end subroutine trc_scf_gradient

end module trc_scf_grad
