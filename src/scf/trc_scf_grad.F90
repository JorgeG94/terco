!
! The full analytic HF/UHF nuclear gradient: one-electron (trc_1e_gradient),
! two-electron (eri%gradient) and nuclear repulsion, assembled at a
! converged trc_scf_result_t.
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
      type(trc_linalg_t), intent(inout) :: lp
      integer,  intent(in)  :: nao, nspin
      integer,  intent(in)  :: nocc(2)
      real(dp), intent(in)  :: cmo(nao, nao, nspin), eps(nao, nspin)
      real(dp), intent(out) :: wmat(nao, nao)
      real(dp), allocatable :: cw(:, :)
      real(dp) :: occ_s
      integer :: s, i, j, nk

      allocate (cw(nao, nao))
      occ_s = merge(2.0_dp, 1.0_dp, nspin == 1)
      do s = 1, nspin
         nk = nocc(s)
         if (nk <= 0) cycle
         do j = 1, nk
            do i = 1, nao
               cw(i, j) = cmo(i, j, s)*eps(j, s)
            end do
         end do
         call lp%gemm_mat(cw(:, 1:nk), cmo(:, 1:nk, s), wmat, transb="T", &
                          alpha=occ_s, beta=merge(0.0_dp, 1.0_dp, s == 1))
      end do
      deallocate (cw)
   end subroutine energy_weighted_density

   subroutine trc_scf_gradient(b, lp, nspin, nalpha, nbeta, cmo, eps, thresh, grad)
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
      type(trc_basis_t),  intent(inout) :: b
      type(trc_linalg_t), intent(inout) :: lp
      integer,            intent(in)    :: nspin, nalpha, nbeta
      real(dp),           intent(in)    :: cmo(b%nao, b%nao, nspin), eps(b%nao, nspin)
      real(dp),           intent(in)    :: thresh
      real(dp),           intent(out)   :: grad(3, b%natm)

      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      integer :: nao, s
      integer :: nocc(2)
      real(dp) :: occ_s
      real(dp), allocatable :: wmat(:, :), dtot(:, :), dens(:, :, :)

      nao = b%nao
      nocc = [nalpha, nbeta]
      occ_s = merge(2.0_dp, 1.0_dp, nspin == 1)

      allocate (dens(nao, nao, nspin))
      do s = 1, nspin
         call lp%density(cmo(:, :, s), nocc(s), occ_s, dens(:, :, s))
      end do

      call pl%build(b, thresh)
      call pl%to_device()
      call eri%build(b, thresh)

      ! --- two-electron, first: this call starts `grad` from zero ---
      call eri%gradient(b, dens, grad)
      grad = 0.5_dp*grad   !! see this module's own header: eri%gradient's
                           !! convention is twice d(e_two)/dR.

      ! --- one-electron, added on top ---
      allocate (wmat(nao, nao), dtot(nao, nao))
      call energy_weighted_density(lp, nao, nspin, cmo, eps, nocc, wmat)
      if (nspin == 1) then
         dtot = dens(:, :, 1)
      else
         dtot = dens(:, :, 1) + dens(:, :, 2)
      end if

      !$acc enter data copyin(wmat, dtot, grad)
      call trc_1e_gradient(b, pl, wmat, dtot, grad)
      !$acc update self(grad)
      !$acc exit data delete(wmat, dtot, grad)

      ! --- nuclear repulsion, host arithmetic, added on top ---
      call nuclear_repulsion_gradient(b, grad)

      call eri%release()
      call pl%release()
      deallocate (wmat, dtot, dens)
   end subroutine trc_scf_gradient

end module trc_scf_grad
