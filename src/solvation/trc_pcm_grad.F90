!
! Analytic nuclear gradient of the PCM (C-PCM / IEF-PCM) solvation energy.
!
! E = 0.5 q.v with K q = R v; q is VARIATIONAL in the sense that substitutes
! away, not in the sense of an envelope-theorem shortcut -- q = K^-1 R v is
! plugged straight in, E = 0.5 v^T (K^-1 R) v, and differentiating that
! product rule gives
!
!     dE/dR = sum_k q_k d(v_nuc_k - v_elec_k)/dR            (grad_qv + grad_nuc)
!           + 0.5 vK_1^T dR_mat/dR v - 0.5 vK_1^T dK/dR q    (grad_solver)
!
! with vK_1 solving K^T x = v -- a SECOND solve through the already-factored
! K (trc_pcm's own charge solve is the untransposed K q = R v). No dq/dR
! term anywhere: that is the substitution, not an approximation. Ported
! term-for-term from pyscf.solvent.grad.pcm's grad_qv/grad_nuc/grad_solver
! but re-derived from this one identity rather than transcribed, since the
! -1/2 contraction's sign and the IEF cross terms are easy to get wrong by
! eye and this way they are checked (test/check_pcm_grad_fd.F90) against
! this module's own derivation, not against a copy of someone else's.
!
! K = S - fac*D*A*S, R = -f_eps*I + fac*D*A, fac = f_eps/(2 pi) for IEF-PCM
! and 0 for C-PCM (whose K is exactly S and R exactly -f_eps*I with fac=0
! substituted -- no special-casing needed in the algebra, only in which
! matrices get built at all). Differentiating,
!
!     dK = dS - fac*(dD*A*S + D*dA*S + D*A*dS)
!     dR = fac*(dD*A + D*dA)
!
! substituted into the grad_solver identity and collected by WHICH
! derivative each term carries:
!
!     dS-bilinear:  bra_s(i) q(j) dS_ij/dr_i,    bra_s = 0.5*(fac*A.*Dtv - vK_1)
!     dD-bilinear:  bra_d(i) ket_d(j) dD_ij/dr_i, bra_d = 0.5*fac*vK_1,
!                   ket_d = A.*(v + S q)
!     dA/dF fanout: combined_scalar(i) = weight_a(i)*area(i)
!                   + bra_s(i)*q(i)*dSii_dF(i)*switch_f(i), contracted with
!                   d(ln switch)/dR -- the diagonal dS term (which rides on
!                   dF, the SAME switching derivative, scaled by F instead
!                   of A) folds into the same per-point scalar as the dA/dF
!                   piece, weight_a(i) = 0.5*fac*Dtv(i)*(v(i)+Sq(i))
!
! where Dtv = D^T vK_1, Sq = S q are plain matrix-vector products. S, D
! (NOT D*A -- that is build_operators' own and not what this needs) are
! rebuilt locally with pcm_s_fill/pcm_da_fill(...,ones,...): they are not kept past
! build_operators, and rebuilding them here is the SAME O(n^2) cost
! build_operators itself already pays once per geometry.
!
! OFF-DIAGONAL S_ij, D_ij depend on r_i, r_j only through their difference
! (D_ij's normal is point j's own FIXED direction cosine from its sphere
! centre, not a function of either position), so d/dr_j = -d/dr_i exactly:
! one derivative formula per kernel, its negative for the other point,
! bucketed by the TWO points' own parent atoms (trc_pcm_grad_pair_item).
!
! THE SWITCHING FUNCTION is different: every point's cavity depends on
! EVERY other atom's sphere, so dF/dA is an O(n*natm) fanout
! (trc_pcm_grad_fanout_item), not a second O(n^2) pass.
!
! GRAD_NUC (nucleus-vs-surface-charge) turns out to be the SAME radial
! function as the dS kernel, with the point's own zeta in place of a
! combined zeta_ij (a nucleus is a true point charge, not a second smeared
! Gaussian, so there is nothing to combine) -- see pcm_grad_g's own
! derivation note.
!
! GRAD_QV (the basis-derivative, electron-density piece) is the one actual
! integral-derivative kernel, generated: trc_integrals%trc_pcm_gradient,
! from scripts/gen_1e.py --pcm-grad. Everything else here is dense,
! analytic host/device arithmetic, no Boys function, no VRR.
!
module trc_pcm_grad
   use pic_types, only: dp
   use trc_constants, only: PI
   use trc_linalg, only: trc_linalg_t
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_pcm_gradient
   use trc_pcm, only: trc_pcm_t, PCM_MODEL_NONE, PCM_MODEL_IEFPCM, pcm_s_fill, pcm_da_fill
   use trc_sph, only: sph_sph_to_cart
   implicit none
   private
   public :: pcm_grad_add, pcm_grad_scf_add

contains

   subroutine pcm_grad_add(this, la, b, pl, dtot, grad)
      !! The full PCM gradient, ADDED into grad(3, natm) -- same convention
      !! as trc_1e_gradient. `this` must already be built AND solved at
      !! `dtot` (this%v, this%q current; the caller's `pcm%operator` call).
      class(trc_pcm_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: la
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: dtot(b%nao, b%nao)
      real(dp), intent(inout) :: grad(3, b%natm)

      real(dp), allocatable :: vk1(:), smat(:, :), sq(:), q_unsym(:)
      real(dp), allocatable :: bra_s(:), bra_d(:), ket_d(:), combined(:)
      logical :: has_d
      integer :: n

      n = this%surface%n_points
      has_d = (this%model == PCM_MODEL_IEFPCM)

      allocate (vk1(n), smat(n, n), sq(n), q_unsym(n), bra_s(n), bra_d(n), ket_d(n), combined(n))
      !$acc enter data create(vk1, smat, sq, q_unsym, bra_s, bra_d, ket_d, combined)

      call pcm_grad_dcopy(n, this%v, vk1)
      call la%getrs("T", n, 1, this%kmat, n, this%ipiv, vk1, n)
      call pcm_s_fill(n, this%surface%coords, this%surface%zeta, this%surface%switch_f, smat)

      ! q_unsym = K^-1 R v, the SAME un-symmetrised solve pcm_solve itself
      ! does before averaging it with R^T K^-T v into this%q (its own
      ! comment). grad_solver's dK contraction needs THIS one, not the
      ! symmetrised this%q: E = 0.5 v^T (K^-1 R) v symmetrises the OPERATOR
      ! on its own (v^T M v = v^T M^T v for any M, a scalar's own
      ! transpose), so differentiating it needs M v wherever "the charge"
      ! appears, not 0.5(Mv + M^T v). For C-PCM, K^-1 R is already
      ! symmetric (K=S symmetric, R=-f_eps I), so q_unsym = this%q exactly
      ! and no second solve is needed.
      if (has_d) then
         call la%gemv("N", n, n, 1.0_dp, this%rmat, n, this%v, 0.0_dp, q_unsym)
         call la%getrs("N", n, 1, this%kmat, n, this%ipiv, q_unsym, n)
      else
         call pcm_grad_dcopy(n, this%q, q_unsym)
      end if
      call la%gemv("N", n, n, 1.0_dp, smat, n, q_unsym, 0.0_dp, sq)

      if (has_d) then
         call pcm_grad_ief_vecs(this, la, n, vk1, sq, q_unsym, bra_s, bra_d, ket_d, combined)
      else
         call pcm_grad_cpcm_vecs(n, vk1, q_unsym, this%surface%switch_f, this%surface%zeta, &
                                 bra_s, bra_d, ket_d, combined)
      end if

      call pcm_grad_offdiag(n, b%natm, this%surface%coords, this%surface%normals, &
                            this%surface%zeta, this%surface%parent, bra_s, q_unsym, &
                            has_d, bra_d, ket_d, grad)
      call pcm_grad_fanout(n, b%natm, this%surface%coords, this%surface%parent, b%at_r, &
                           this%surface%r_switch, this%surface%r_inner, combined, grad)
      call pcm_grad_nuc(n, b%natm, this%surface%coords, this%surface%zeta, &
                        this%surface%parent, this%q, b%at_z, b%at_r, grad)

      !$acc exit data delete(vk1, smat, sq, q_unsym, bra_s, bra_d, ket_d, combined)
      deallocate (vk1, smat, sq, q_unsym, bra_s, bra_d, ket_d, combined)

      ! --- the one integral-derivative piece: basis-centre + point-centre
      ! derivative of sum_k q_k erf(zeta_k r)/r, contracted with the
      ! density -- ADDS into grad, same as everything above. ---
      call trc_pcm_gradient(b, pl, this%surface%coords, this%surface%zeta, &
                            this%surface%parent, this%q, dtot, grad)
   end subroutine pcm_grad_add

   subroutine pcm_grad_scf_add(b, lp, nocc, cmo, eps, thresh, pcm, grad)
      !! A convenience wrapper around `pcm_grad_add` for a caller that only
      !! has the converged orbitals (cmo/eps/nocc), not the SCF gradient's
      !! own already-built `dtot`/pair list -- trc_capi.F90's
      !! `trc_run_gradient` is the one caller. Deliberately NOT wired into
      !! trc_scf_grad%trc_scf_gradient's own argument list or body: doing
      !! that pushed its argument count (24 -> 25) and line count (316 ->
      !! 352) past tools/arch_baseline.json's recorded measures, which
      !! AGENTS.md's "restructure instead" asks to avoid rather than
      !! accept as new debt. The cost is one extra O(n^2) pair-list build
      !! the vacuum gradient's own `pl_own`/`pair_eri_cache%pl` already
      !! paid for; a later change can thread that reuse through without
      !! touching either routine's contract.
      !!
      !! Builds its own density from cmo/eps/nocc, with the SAME
      !! spherical-to-Cartesian fold trc_scf_gradient's own header
      !! documents (T D T^T, exact, not an approximation) when `b%spherical`.
      type(trc_basis_t), intent(inout) :: b
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nocc(2)
      real(dp), intent(in) :: cmo(:, :, :), eps(:, :)
      real(dp), intent(in) :: thresh
      class(trc_pcm_t), intent(inout) :: pcm
      real(dp), intent(inout) :: grad(3, b%natm)

      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: dtot(:, :), vmat_pcm(:, :)
      real(dp) :: e_diel_pcm, occ_s
      integer :: nao, nsp, nspin

      if (pcm%model == PCM_MODEL_NONE) return
      nao = b%nao
      nsp = merge(b%nao_sph, b%nao, b%spherical)
      nspin = size(cmo, 3)
      occ_s = merge(2.0_dp, 1.0_dp, nspin == 1)

      allocate (dtot(nao, nao), vmat_pcm(nao, nao))
      !$acc enter data copyin(cmo, eps) create(dtot, vmat_pcm)

      if (.not. b%spherical) then
         call pcm_grad_density(lp, nao, nspin, nocc, occ_s, cmo, eps, dtot)
      else
         ! T(Da+Db)T^T = T Da T^T + T Db T^T (sph_sph_to_cart is linear),
         ! so folding the ALREADY-SUMMED spherical-space total density is
         ! exact and needs no per-channel fold.
         block
            real(dp), allocatable :: dtot_s(:, :)
            allocate (dtot_s(nsp, nsp))
            !$acc enter data create(dtot_s)
            call pcm_grad_density(lp, nsp, nspin, nocc, occ_s, cmo, eps, dtot_s)
            call sph_sph_to_cart(b, dtot_s, dtot)
            !$acc exit data delete(dtot_s)
            deallocate (dtot_s)
         end block
      end if

      call pl%build(b, thresh)
      call pl%to_device()
      call pcm%operator(lp, b, pl, dtot, vmat_pcm, e_diel_pcm)
      call pcm_grad_add(pcm, lp, b, pl, dtot, grad)
      call pl%release()

      !$acc exit data delete(dtot, vmat_pcm, cmo, eps)
      deallocate (dtot, vmat_pcm)
   end subroutine pcm_grad_scf_add

   subroutine pcm_grad_density(lp, n, nspin, nocc, occ_s, cmo, eps, dtot)
      !! The total density (alpha + beta) at `n` AOs, built from `nspin`
      !! channels of cmo/eps the same way trc_scf_gradient's own header
      !! does (lp%density per channel, summed) -- split out so
      !! pcm_grad_scf_add's spherical and Cartesian branches share it.
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: n, nspin
      integer, intent(in) :: nocc(2)
      real(dp), intent(in) :: occ_s
      real(dp), intent(in) :: cmo(:, :, :), eps(:, :)
      real(dp), intent(out) :: dtot(n, n)
      real(dp), allocatable :: density(:, :, :)
      integer :: s

      allocate (density(n, n, nspin))
      !$acc enter data create(density)
      do s = 1, nspin
         call lp%density(cmo(:, :, s), nocc(s), occ_s, density(:, :, s))
      end do
      if (nspin == 1) then
         call pcm_grad_dcopy(n*n, density(:, :, 1), dtot)
      else
         call pcm_grad_daxpy2(n*n, density(:, :, 1), density(:, :, 2), dtot)
      end if
      !$acc exit data delete(density)
      deallocate (density)
   end subroutine pcm_grad_density

   subroutine pcm_grad_daxpy2(n, x, y, z)
      !! z := x + y, device-resident arrays all three
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n), y(n)
      real(dp), intent(out) :: z(n)
      integer :: i
      do concurrent(i = 1:n)
         z(i) = x(i) + y(i)
      end do
   end subroutine pcm_grad_daxpy2

   subroutine pcm_grad_dcopy(n, x, y)
      !! y := x, device-resident arrays both
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: y(n)
      integer :: i
      do concurrent(i = 1:n)
         y(i) = x(i)
      end do
   end subroutine pcm_grad_dcopy

   subroutine pcm_grad_cpcm_vecs(n, vk1, q, switch_f, zeta, bra_s, bra_d, ket_d, combined)
      !! C-PCM: fac = 0, so bra_s = -0.5 vK_1 and the D/A terms vanish
      !! outright (K = S exactly, R = -f_eps I exactly -- no D, no A was
      !! ever built). Only the diagonal dSii/dF piece survives into
      !! `combined`.
      integer, intent(in) :: n
      real(dp), intent(in) :: vk1(n), q(n), switch_f(n), zeta(n)
      real(dp), intent(out) :: bra_s(n), bra_d(n), ket_d(n), combined(n)
      integer :: i
      real(dp) :: dsii_df
      do concurrent(i = 1:n) local(dsii_df)
         bra_s(i) = -0.5_dp*vk1(i)
         bra_d(i) = 0.0_dp
         ket_d(i) = 0.0_dp
         dsii_df = -zeta(i)*sqrt(2.0_dp/PI)/switch_f(i)**2
         combined(i) = bra_s(i)*q(i)*dsii_df*switch_f(i)
      end do
   end subroutine pcm_grad_cpcm_vecs

   subroutine pcm_grad_ief_vecs(this, la, n, vk1, sq, q_unsym, bra_s, bra_d, ket_d, combined)
      !! IEF-PCM: builds raw D (pcm_da_fill with an all-ones area, the SAME
      !! trick pcm_sd_matrices uses to get D rather than D*A), forms
      !! Dtv = D^T vK_1, then bra_s/bra_d/ket_d/combined per this module's
      !! own header.
      class(trc_pcm_t), intent(in) :: this
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: n
      real(dp), intent(in) :: vk1(n), sq(n), q_unsym(n)
      real(dp), intent(out) :: bra_s(n), bra_d(n), ket_d(n), combined(n)

      real(dp), allocatable :: draw(:, :), ones(:), dtv(:)
      real(dp) :: fac
      integer :: i
      real(dp) :: dsii_df, weight_a

      allocate (draw(n, n), ones(n), dtv(n))
      !$acc enter data create(draw, dtv) copyin(ones)
      do concurrent(i = 1:n)
         ones(i) = 1.0_dp
      end do
      call pcm_da_fill(n, this%surface%coords, this%surface%normals, this%surface%zeta, &
                   this%surface%r_parent, ones, draw)
      call la%gemv("T", n, n, 1.0_dp, draw, n, vk1, 0.0_dp, dtv)

      fac = this%f_eps/(2.0_dp*PI)
      do concurrent(i = 1:n) local(dsii_df, weight_a)
         bra_s(i) = 0.5_dp*(fac*this%surface%area(i)*dtv(i) - vk1(i))
         bra_d(i) = 0.5_dp*fac*vk1(i)
         ket_d(i) = this%surface%area(i)*(this%v(i) + sq(i))
         weight_a = 0.5_dp*fac*dtv(i)*(this%v(i) + sq(i))
         dsii_df = -this%surface%zeta(i)*sqrt(2.0_dp/PI)/this%surface%switch_f(i)**2
         combined(i) = weight_a*this%surface%area(i) + bra_s(i)*q_unsym(i)*dsii_df*this%surface%switch_f(i)
      end do

      !$acc exit data delete(draw, dtv, ones)
      deallocate (draw, ones, dtv)
   end subroutine pcm_grad_ief_vecs

   pure real(dp) function pcm_grad_g(r, zeta_r) result(g)
      !! -(1/r) d/dr [erf(zeta_r r)/r], the radial function both the
      !! off-diagonal dS kernel and grad_nuc are built from (a nucleus is
      !! the SAME erf(zeta r)/r interaction with a single exponent, not a
      !! combined zeta_ij, since there is only one smeared Gaussian in that
      !! pair). Sympy-checked against trc_pcm's own d_elem formula (its
      !! off-diagonal D_ij IS g(r)*nr, by definition): see trc_pcm_grad's
      !! own module header.
      !$acc routine seq
      real(dp), intent(in) :: r, zeta_r
      g = erf(zeta_r*r)/r**3 - 2.0_dp*zeta_r*exp(-zeta_r*zeta_r*r*r)/(sqrt(PI)*r*r)
   end function pcm_grad_g

   pure real(dp) function pcm_grad_gp(r, zij) result(gp)
      !! d/dr of pcm_grad_g -- the extra piece dD_ij/dr_i needs on top of
      !! g(r)*n_j (the normal itself is fixed, so only g(r) varies with r).
      !$acc routine seq
      real(dp), intent(in) :: r, zij
      gp = exp(-zij*zij*r*r)*(4.0_dp*zij**3*r*r + 6.0_dp*zij)/(sqrt(PI)*r**3) - 3.0_dp*erf(zij*r)/r**4
   end function pcm_grad_gp

   subroutine pcm_grad_offdiag(n, natm, coords, normals, zeta, parent, bra_s, q, &
                               has_d, bra_d, ket_d, grad)
      !! The dS (always) and dD (IEF only) bilinear terms, O(n^2) over
      !! UNORDERED pairs i < j rather than every ordered pair: g(r)/gp(r)
      !! (pcm_grad_g/_gp, the erf/exp-bearing radial pieces) depend only on
      !! r = |coords(i)-coords(j)| and zij = zeta(i)*zeta(j)/sqrt(zeta(i)^2
      !! + zeta(j)^2), both symmetric in i and j, so the ordered pair (j,i)
      !! needs the SAME g/gp the pair (i,j) already paid the sqrt/erf/exp
      !! for. pcm_grad_pair_item below folds both ordered contributions
      !! into one call, halving the transcendental-heavy work (and, for
      !! the dS term, the atomic scatter too -- see its own header).
      integer, intent(in) :: n, natm
      real(dp), intent(in) :: coords(3, n), normals(3, n), zeta(n)
      integer, intent(in) :: parent(n)
      real(dp), intent(in) :: bra_s(n), q(n)
      logical, intent(in) :: has_d
      real(dp), intent(in) :: bra_d(n), ket_d(n)
      real(dp), intent(inout) :: grad(3, natm)
      integer :: i, j
      do concurrent(j = 1:n, i = 1:n)
         if (i < j) then
            call pcm_grad_pair_item(i, j, n, natm, coords, normals, zeta, parent, &
                                    bra_s, q, has_d, bra_d, ket_d, grad)
         end if
      end do
   end subroutine pcm_grad_offdiag

   pure subroutine pcm_grad_pair_item(i, j, n, natm, coords, normals, zeta, parent, &
                                      bra_s, q, has_d, bra_d, ket_d, grad)
      !! Called once per unordered pair i < j (pcm_grad_offdiag's `if`),
      !! folding in BOTH ordered contributions (i,j) and (j,i).
      !!
      !! dS: dS_ji/dr_j = -dS_ij/dr_i exactly (translation invariance, this
      !! module's header), so bra_s(i)*q(j)*dS_ij/dr_i (pair (i,j)) and
      !! bra_s(j)*q(i)*dS_ji/dr_j (pair (j,i)) share the SAME direction
      !! vector dsx/dsy/dsz up to its own sign and collapse into one
      !! combined coefficient `coefs`, scattered once -- 6 atomics for
      !! this unordered pair where the old ordered-pair loop used 12.
      !!
      !! dD: D_ij's ket normal is point j's own (this module's header), so
      !! (i,j) and (j,i) read DIFFERENT normals and do not collapse the
      !! same way; they still share gval/gpval/r, computed once.
      !$acc routine seq
      integer,  value       :: i, j, n, natm
      real(dp), intent(in)  :: coords(3, n), normals(3, n), zeta(n)
      integer,  intent(in)  :: parent(n)
      real(dp), intent(in)  :: bra_s(n), q(n)
      logical,  value       :: has_d
      real(dp), intent(in)  :: bra_d(n), ket_d(n)
      real(dp), intent(inout) :: grad(3, natm)

      real(dp) :: vx, vy, vz, r, zij, gval, dsx, dsy, dsz, coefs
      real(dp) :: gpval, nrij, nrji, ddx_ij, ddy_ij, ddz_ij, ddx_ji, ddy_ji, ddz_ji
      real(dp) :: coefd_ij, coefd_ji
      integer  :: ai, aj

      vx = coords(1, i) - coords(1, j)
      vy = coords(2, i) - coords(2, j)
      vz = coords(3, i) - coords(3, j)
      r = sqrt(vx*vx + vy*vy + vz*vz)
      zij = zeta(i)*zeta(j)/sqrt(zeta(i)**2 + zeta(j)**2)
      gval = pcm_grad_g(r, zij)
      dsx = -gval*vx
      dsy = -gval*vy
      dsz = -gval*vz
      coefs = bra_s(i)*q(j) + bra_s(j)*q(i)

      ai = parent(i)
      aj = parent(j)
      !$acc atomic update
      grad(1, ai) = grad(1, ai) + coefs*dsx
      !$acc atomic update
      grad(2, ai) = grad(2, ai) + coefs*dsy
      !$acc atomic update
      grad(3, ai) = grad(3, ai) + coefs*dsz
      !$acc atomic update
      grad(1, aj) = grad(1, aj) - coefs*dsx
      !$acc atomic update
      grad(2, aj) = grad(2, aj) - coefs*dsy
      !$acc atomic update
      grad(3, aj) = grad(3, aj) - coefs*dsz

      if (has_d) then
         gpval = pcm_grad_gp(r, zij)
         nrij = vx*normals(1, j) + vy*normals(2, j) + vz*normals(3, j)
         ddx_ij = gpval*(nrij/r)*vx + gval*normals(1, j)
         ddy_ij = gpval*(nrij/r)*vy + gval*normals(2, j)
         ddz_ij = gpval*(nrij/r)*vz + gval*normals(3, j)
         coefd_ij = bra_d(i)*ket_d(j)
         ! (j,i): ket point i, its own normal; direction vector is -v, not
         ! a sign-flip of ddx_ij (DIFFERENT normal, this routine's header).
         nrji = -vx*normals(1, i) - vy*normals(2, i) - vz*normals(3, i)
         ddx_ji = -gpval*(nrji/r)*vx + gval*normals(1, i)
         ddy_ji = -gpval*(nrji/r)*vy + gval*normals(2, i)
         ddz_ji = -gpval*(nrji/r)*vz + gval*normals(3, i)
         coefd_ji = bra_d(j)*ket_d(i)
         ! Two atomics per component, not one combined update: OpenACC's
         ! `atomic update` requires the single-operator form `x = x op
         ! expr` (nvfortran rejects `x = x + a - b` outright, "Invalid
         ! atomic update expression"). Still only ONE gval/gpval/r per
         ! unordered pair -- this routine's own header.
         !$acc atomic update
         grad(1, ai) = grad(1, ai) + coefd_ij*ddx_ij
         !$acc atomic update
         grad(1, ai) = grad(1, ai) - coefd_ji*ddx_ji
         !$acc atomic update
         grad(2, ai) = grad(2, ai) + coefd_ij*ddy_ij
         !$acc atomic update
         grad(2, ai) = grad(2, ai) - coefd_ji*ddy_ji
         !$acc atomic update
         grad(3, ai) = grad(3, ai) + coefd_ij*ddz_ij
         !$acc atomic update
         grad(3, ai) = grad(3, ai) - coefd_ji*ddz_ji
         !$acc atomic update
         grad(1, aj) = grad(1, aj) - coefd_ij*ddx_ij
         !$acc atomic update
         grad(1, aj) = grad(1, aj) + coefd_ji*ddx_ji
         !$acc atomic update
         grad(2, aj) = grad(2, aj) - coefd_ij*ddy_ij
         !$acc atomic update
         grad(2, aj) = grad(2, aj) + coefd_ji*ddy_ji
         !$acc atomic update
         grad(3, aj) = grad(3, aj) - coefd_ij*ddz_ij
         !$acc atomic update
         grad(3, aj) = grad(3, aj) + coefd_ji*ddz_ji
      end if
   end subroutine pcm_grad_pair_item

   subroutine pcm_grad_fanout(n, natm, coords, parent, atom_r, r_switch, r_inner, combined, grad)
      !! The switching-function (dF/dA) term, O(n*natm): point i's cavity
      !! depends on every OTHER atom's sphere, so this is a fanout, not a
      !! second O(n^2) pass. `combined` already carries the dA and the
      !! diagonal dSii pieces together -- see this module's own header.
      integer, intent(in) :: n, natm
      real(dp), intent(in) :: coords(3, n)
      integer, intent(in) :: parent(n)
      real(dp), intent(in) :: atom_r(3, natm), r_switch(natm), r_inner(natm)
      real(dp), intent(in) :: combined(n)
      real(dp), intent(inout) :: grad(3, natm)
      integer :: i, ja
      do concurrent(ja = 1:natm, i = 1:n)
         if (ja /= parent(i)) then
            call pcm_grad_fanout_item(i, ja, n, natm, coords, parent, atom_r, r_switch, r_inner, combined, grad)
         end if
      end do
   end subroutine pcm_grad_fanout

   pure subroutine pcm_grad_fanout_item(i, ja, n, natm, coords, parent, atom_r, r_switch, r_inner, combined, grad)
      !! d(ln switch_h(d_iJ))/dR, point i's own switching factor against
      !! atom J = ja: SCALAR = combined(i) * h'(d_iJ) / (h(d_iJ) * r_iJ *
      !! r_switch(J)), the log-derivative of SWIG's eq. 3.19 factor,
      !! times the direction (coords_i - atom_J). Atom J moving changes
      !! only ITS factor (scattered with a minus, pyscf's own convention
      !! -- this module's header); point i's own parent atom gets the SAME
      !! quantity with a plus, and summing it over every J one atomic
      !! update at a time reproduces "every factor changes together"
      !! exactly, with no separate pre-summed array needed.
      !$acc routine seq
      integer,  value       :: i, ja, n, natm
      real(dp), intent(in)  :: coords(3, n)
      integer,  intent(in)  :: parent(n)
      real(dp), intent(in)  :: atom_r(3, natm), r_switch(natm), r_inner(natm)
      real(dp), intent(in)  :: combined(n)
      real(dp), intent(inout) :: grad(3, natm)

      real(dp) :: ux, uy, uz, rij, x, hp, fij, scalar_part
      integer  :: ia

      ux = coords(1, i) - atom_r(1, ja)
      uy = coords(2, i) - atom_r(2, ja)
      uz = coords(3, i) - atom_r(3, ja)
      rij = sqrt(ux*ux + uy*uy + uz*uz)
      x = (rij - r_inner(ja))/r_switch(ja)
      if (x <= 0.0_dp .or. x >= 1.0_dp) return
      ! h(x) = x^3(10-15x+6x^2); h'(x) = 30 x^2 (1-x)^2 -- pcm_switch_h's
      ! own derivative (trc_pcm.F90's eq. 3.19), evaluated directly since
      ! 0 < x < 1 here (the clamp above is pcm_switch_h's own, no need to
      ! call it and re-test the same bounds).
      hp = 30.0_dp*x*x*(1.0_dp - x)*(1.0_dp - x)
      fij = x**3*(10.0_dp - 15.0_dp*x + 6.0_dp*x**2)
      scalar_part = combined(i)*hp/(fij*rij*r_switch(ja))
      ia = parent(i)
      !$acc atomic update
      grad(1, ja) = grad(1, ja) - scalar_part*ux
      !$acc atomic update
      grad(2, ja) = grad(2, ja) - scalar_part*uy
      !$acc atomic update
      grad(3, ja) = grad(3, ja) - scalar_part*uz
      !$acc atomic update
      grad(1, ia) = grad(1, ia) + scalar_part*ux
      !$acc atomic update
      grad(2, ia) = grad(2, ia) + scalar_part*uy
      !$acc atomic update
      grad(3, ia) = grad(3, ia) + scalar_part*uz
   end subroutine pcm_grad_fanout_item

   subroutine pcm_grad_nuc(n, natm, coords, zeta, parent, q, at_z, at_r, grad)
      !! Nucleus-vs-surface-charge: the SAME radial function as the dS
      !! kernel (pcm_grad_g), with the point's own zeta -- see this
      !! module's header, GRAD_NUC.
      integer, intent(in) :: n, natm
      real(dp), intent(in) :: coords(3, n), zeta(n)
      integer, intent(in) :: parent(n)
      real(dp), intent(in) :: q(n), at_z(natm), at_r(3, natm)
      real(dp), intent(inout) :: grad(3, natm)
      integer :: k, ia
      do concurrent(ia = 1:natm, k = 1:n)
         call pcm_grad_nuc_item(k, ia, n, natm, coords, zeta, parent, q, at_z, at_r, grad)
      end do
   end subroutine pcm_grad_nuc

   pure subroutine pcm_grad_nuc_item(k, ia, n, natm, coords, zeta, parent, q, at_z, at_r, grad)
      !$acc routine seq
      integer,  value       :: k, ia, n, natm
      real(dp), intent(in)  :: coords(3, n), zeta(n)
      integer,  intent(in)  :: parent(n)
      real(dp), intent(in)  :: q(n), at_z(natm), at_r(3, natm)
      real(dp), intent(inout) :: grad(3, natm)
      real(dp) :: ux, uy, uz, r, g, coef
      integer  :: pa

      ux = coords(1, k) - at_r(1, ia)
      uy = coords(2, k) - at_r(2, ia)
      uz = coords(3, k) - at_r(3, ia)
      r = sqrt(ux*ux + uy*uy + uz*uz)
      g = pcm_grad_g(r, zeta(k))
      coef = q(k)*at_z(ia)*g
      pa = parent(k)
      !$acc atomic update
      grad(1, pa) = grad(1, pa) - coef*ux
      !$acc atomic update
      grad(2, pa) = grad(2, pa) - coef*uy
      !$acc atomic update
      grad(3, pa) = grad(3, pa) - coef*uz
      !$acc atomic update
      grad(1, ia) = grad(1, ia) + coef*ux
      !$acc atomic update
      grad(2, ia) = grad(2, ia) + coef*uy
      !$acc atomic update
      grad(3, ia) = grad(3, ia) + coef*uz
   end subroutine pcm_grad_nuc_item

end module trc_pcm_grad
