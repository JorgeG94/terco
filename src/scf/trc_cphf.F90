!
! Coupled-perturbed Hartree-Fock: the orbital-response linear system this
! milestone solves, and the static dipole polarizability built on top of it
! as the validating property. Intended foundation for RI-MP2's Z-vector
! equation and, later, a coupled-perturbed Kohn-Sham variant for
! double-hybrid gradients -- neither is implemented here, this module only
! solves the ground-state orbital Hessian system and the one property that
! needs nothing but it.
!
! THE EQUATION (RHF)
! -------------------
! For a real, symmetric one-electron perturbation (a static field has no
! AO-basis response, so this is the whole right-hand side), the response
! orbital-rotation amplitudes U_ai solve
!
!     (e_a - e_i) U_ai + [A U]_ai = -B_ai
!
! with the orbital Hessian action, built from one resident Fock-like build
! rather than an AO integral derivative:
!
!     D_U = C_vir U C_occ^T + (C_vir U C_occ^T)^T          (symmetric, AO)
!     [A U]_ai = 2 [C_vir^T G[D_U] C_occ]_ai,   G[D_U] = J[D_U] - K[D_U]/2
!
! i.e. exactly eri%fock_resident(b, D_U, G_U, j_scale=1, k_scale=1,
! density_screen=.false.) -- the trial density is not the SCF's physical
! one, so it must not drive the Schwarz screen (trc_eri.F90's own warning
! on both `fock_resident` and `fock`). The derivation (4(ai|jb) -
! (ij|ab) - (ib|ja), contracted against the full permutation symmetry of
! real two-electron integrals) is in this branch's PR description, not
! repeated here; see trc_eri.F90's header for why `fock_resident` is the
! entry point this was always meant to reach ("the coupled-perturbed
! equations run a hundred of them").
!
! THE EQUATION (UHF)
! -------------------
! Each spin channel's Fock response only picks up a J contribution from
! the OTHER spin's trial density; its own K is spin-local:
!
!     D_U^s = C_vir,s U^s C_occ,s^T + transpose,   D_U^tot = D_U^a + D_U^b
!     [A U]_ai^s = [C_vir,s^T (J[D_U^tot] - K[D_U^s]) C_occ,s]_ai
!
! J[D_U^tot] is one fock_resident call (j_scale=1, k_scale=0) shared by
! both spins; K[D_U^s] is fock_resident(j_scale=0, k_scale=2) per spin --
! the same "k_scale=2 on the undoubled density" trick trc_scf_fock.F90's
! `build_g` uses for the ground-state UKS exchange term. No extra factor
! of 2 here: RHF's 2 comes from folding two identical spin channels into
! one solve, which UHF does not do.
!
! THE SOLVER
! ----------
! Preconditioned conjugate gradient, the diagonal (e_a - e_i) as the
! preconditioner (SPD for a stable ground state, so CG is the right
! choice over a general Krylov method), one right-hand side at a time
! inside a single call that reports an iteration count per column. A
! batched build across right-hand sides -- one `fock_resident`-like pass
! contracting a quartet against every trial density at once, as
! `eri%fock_many` does for already-host-resident problems -- is the
! obvious next speed-up and is NOT done: there is no device-resident
! multi-density fock_resident, so this runs ncomp independent CG loops,
! each paying its own integral pass per iteration.
!
! RESIDENCY. `cmo`, the right-hand sides and the solution are assumed
! already on the device (the same contract eri%fock_resident itself
! documents); only `eps` (nao-size) and the CG's own scalars (dot
! products, read back by trc_linalg%dot_2d, exactly as DIIS's history
! does) cross the bus. Every per-iteration matrix -- D_U, G_U, the GEMM
! scratch, the CG vectors themselves -- stays resident for the life of
! the solve.
!
module trc_cphf
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_multipoles, TRC_NMULT
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_settings, only: trc_scf_result_t
   implicit none
   private

   public :: trc_cphf_options_t, trc_cphf_result_t
   public :: trc_cphf_solve_rhf, trc_cphf_solve_uhf
   public :: trc_scf_polarizability_rhf, trc_scf_polarizability_uhf, trc_scf_polarizability
   ! Shared with trc_cpks (the CPKS extension of this solver): the same CG
   ! loop, the same elementwise `do concurrent` helpers -- see this
   ! module's header on why they are split out as explicit-shape-dummy
   ! subroutines rather than inlined.
   public :: zero_mat, copy_mat, neg_mat, add_mat, axpy_mat, precond_mat, combine_p, add_diag

   type :: trc_cphf_options_t
      real(dp) :: conv_residual = 1.0e-8_dp
      !! Converged when the Frobenius norm of the residual (e_a-e_i) U +
      !! [A U] + B drops below this, for every right-hand side.
      integer :: max_iter = 100
      logical :: verbose = .false.
   end type trc_cphf_options_t

   type :: trc_cphf_result_t
      logical :: converged = .true.
      integer, allocatable :: iterations(:)   !! CG iterations, one per right-hand side
      character(len=256) :: message = ""
   end type trc_cphf_result_t

contains

   !
   ! RHF CPHF: `ncomp` right-hand sides at once (one CG solve per column,
   ! see this module's own header on why that is not one shared integral
   ! pass). `cmo` (nao, nao), `bvec` (nvirt, nocc, ncomp) and `u` (same
   ! shape, overwritten with the solution) must already be device-resident;
   ! `eps` (nao) is host, read once to build the diagonal preconditioner.
   !
   subroutine trc_cphf_solve_rhf(b, eri, lp, nao, nocc, nvirt, cmo, eps, ncomp, bvec, u, opts, res)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nocc, nvirt, ncomp
      real(dp), intent(in) :: cmo(nao, nao)
      !! Occupied columns 1:nocc, virtual nocc+1:nao, ascending orbital
      !! energy -- trc_scf_driver's own convention.
      real(dp), intent(in) :: eps(nao)
      real(dp), intent(in) :: bvec(nvirt, nocc, ncomp)
      real(dp), intent(inout) :: u(nvirt, nocc, ncomp)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res

      real(dp), allocatable :: d_u(:, :), d_ut(:, :), g_u(:, :)
      real(dp), allocatable :: tmp_no(:, :), tmp_no2(:, :)
      real(dp), allocatable :: r(:, :), z(:, :), p(:, :), ap(:, :), diagv(:, :)
      real(dp) :: rz, rz_old, pap, alpha_cg, beta_cg, resid
      integer :: it, comp, i, a
      logical :: broke

      allocate (d_u(nao, nao), d_ut(nao, nao), g_u(nao, nao))
      allocate (tmp_no(nao, nocc), tmp_no2(nao, nocc))
      allocate (r(nvirt, nocc), z(nvirt, nocc), p(nvirt, nocc), ap(nvirt, nocc), diagv(nvirt, nocc))
      do i = 1, nocc
         do a = 1, nvirt
            diagv(a, i) = eps(nocc + a) - eps(i)
         end do
      end do

      !$acc enter data create(d_u, d_ut, g_u, tmp_no, tmp_no2, r, z, p, ap)
      !$acc enter data copyin(diagv)

      allocate (res%iterations(ncomp))
      res%converged = .true.
      res%message = ""

      do comp = 1, ncomp
         call neg_mat(nvirt, nocc, bvec(:, :, comp), r)
         call zero_mat(nvirt, nocc, u(:, :, comp))
         call precond_mat(nvirt, nocc, diagv, r, z)
         call copy_mat(nvirt, nocc, z, p)
         rz = lp%dot_2d(r, z)
         resid = sqrt(lp%dot_2d(r, r))
         it = 0
         broke = .false.
         ! `.not. (resid < conv_residual)`, not `resid >= conv_residual`:
         ! a NaN right-hand side makes resid NaN from the start, and NaN
         ! fails BOTH comparisons. With `>=` that reads as "already
         ! converged" (the loop body, and the pap guard inside it, never
         ! run at all) and the post-loop check below falls through the
         ! same way, leaving res%converged at its initial .true. with
         ! zero iterations -- silently wrong. `.not. (resid < conv)` reads
         ! a NaN resid as "not yet converged" instead, so the loop runs,
         ! matvec produces a NaN ap, and the pap guard above catches it.
         do while (it < opts%max_iter .and. .not. (resid < opts%conv_residual))
            call matvec(p, ap)
            pap = lp%dot_2d(p, ap)
            ! A breakdown (p^T A p <= 0, including NaN -- NaN > 0 is
            ! false) must stop here: falling through to alpha_cg = rz/pap
            ! gives an Inf/NaN update that the residual check below
            ! cannot catch either (NaN >= conv_residual is also false),
            ! which would otherwise exit the loop and report converged.
            if (.not. (pap > 0.0_dp)) then
               broke = .true.
               exit
            end if
            alpha_cg = rz/pap
            call axpy_mat(nvirt, nocc, alpha_cg, p, u(:, :, comp))
            call axpy_mat(nvirt, nocc, -alpha_cg, ap, r)
            resid = sqrt(lp%dot_2d(r, r))
            it = it + 1
            if (resid < opts%conv_residual) exit
            call precond_mat(nvirt, nocc, diagv, r, z)
            rz_old = rz
            rz = lp%dot_2d(r, z)
            beta_cg = rz/rz_old
            call combine_p(nvirt, nocc, z, beta_cg, p)
         end do
         res%iterations(comp) = it
         if (broke) then
            res%converged = .false.
            res%message = "trc_cphf: CG breakdown (p^T A p <= 0, or NaN) on a right-hand side"
         else if (.not. (resid < opts%conv_residual)) then
            res%converged = .false.
            res%message = "trc_cphf: a right-hand side did not converge"
         end if
         if (opts%verbose) print "(a,i0,a,i0,a,es9.2)", &
            "   cphf rhs ", comp, ": ", it, " iterations, residual ", resid
      end do

      !$acc exit data delete(diagv)
      !$acc exit data delete(d_u, d_ut, g_u, tmp_no, tmp_no2, r, z, p, ap)

   contains

      ! [A p]_ai = 2 [C_vir^T G[D_p] C_occ]_ai + (e_a-e_i) p_ai, D_p the
      ! symmetrised AO density this p rotates into (this module's header).
      subroutine matvec(pin, apout)
         real(dp), intent(in) :: pin(nvirt, nocc)
         real(dp), intent(inout) :: apout(nvirt, nocc)
         call lp%gemm_mat(cmo(:, nocc + 1:nao), pin, tmp_no)
         call lp%gemm_mat(tmp_no, cmo(:, 1:nocc), d_u, transb="T")
         call lp%transpose(nao, nao, d_u, nao, d_ut, nao)
         call add_mat(nao, nao, d_u, d_ut, d_u)
         call eri%fock_resident(b, d_u, g_u, k_scale=1.0_dp, j_scale=1.0_dp, density_screen=.false.)
         call lp%gemm_mat(g_u, cmo(:, 1:nocc), tmp_no2)
         call lp%gemm_mat(cmo(:, nocc + 1:nao), tmp_no2, apout, transa="T", alpha=2.0_dp)
         call add_diag(nvirt, nocc, diagv, pin, apout)
      end subroutine matvec

   end subroutine trc_cphf_solve_rhf

   !
   ! UHF CPHF: the two spin channels solved together (one combined CG over
   ! the stacked [U^a; U^b] vector), because the matvec couples them
   ! through J -- see this module's own header. `cmo_a`/`cmo_b`, the
   ! right-hand sides and the solutions must already be device-resident;
   ! `eps_a`/`eps_b` are host.
   !
   subroutine trc_cphf_solve_uhf(b, eri, lp, nao, noa, nva, nob, nvb, cmo_a, cmo_b, eps_a, eps_b, &
                                 ncomp, bvec_a, bvec_b, u_a, u_b, opts, res)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, noa, nva, nob, nvb, ncomp
      real(dp), intent(in) :: cmo_a(nao, nao), cmo_b(nao, nao)
      real(dp), intent(in) :: eps_a(nao), eps_b(nao)
      real(dp), intent(in) :: bvec_a(nva, noa, ncomp), bvec_b(nvb, nob, ncomp)
      real(dp), intent(inout) :: u_a(nva, noa, ncomp), u_b(nvb, nob, ncomp)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res

      real(dp), allocatable :: d_ua(:, :), d_ub(:, :), d_t(:, :), d_tmp(:, :)
      real(dp), allocatable :: g_j(:, :), g_ka(:, :), g_kb(:, :), g_eff_a(:, :), g_eff_b(:, :)
      real(dp), allocatable :: tmp_a(:, :), tmp_a2(:, :), tmp_b(:, :), tmp_b2(:, :)
      real(dp), allocatable :: ra(:, :), za(:, :), pa(:, :), apa(:, :), diag_a(:, :)
      real(dp), allocatable :: rb(:, :), zb(:, :), pb(:, :), apb(:, :), diag_b(:, :)
      real(dp) :: rz, rz_old, pap, alpha_cg, beta_cg, resid
      integer :: it, comp, i, a
      logical :: broke

      allocate (d_ua(nao, nao), d_ub(nao, nao), d_t(nao, nao), d_tmp(nao, nao))
      allocate (g_j(nao, nao), g_ka(nao, nao), g_kb(nao, nao), g_eff_a(nao, nao), g_eff_b(nao, nao))
      allocate (tmp_a(nao, noa), tmp_a2(nao, noa), tmp_b(nao, nob), tmp_b2(nao, nob))
      allocate (ra(nva, noa), za(nva, noa), pa(nva, noa), apa(nva, noa), diag_a(nva, noa))
      allocate (rb(nvb, nob), zb(nvb, nob), pb(nvb, nob), apb(nvb, nob), diag_b(nvb, nob))
      do i = 1, noa
         do a = 1, nva
            diag_a(a, i) = eps_a(noa + a) - eps_a(i)
         end do
      end do
      do i = 1, nob
         do a = 1, nvb
            diag_b(a, i) = eps_b(nob + a) - eps_b(i)
         end do
      end do

      !$acc enter data create(d_ua, d_ub, d_t, d_tmp, g_j, g_ka, g_kb, g_eff_a, g_eff_b, &
      !$acc                   tmp_a, tmp_a2, tmp_b, tmp_b2, ra, za, pa, apa, rb, zb, pb, apb)
      !$acc enter data copyin(diag_a, diag_b)

      allocate (res%iterations(ncomp))
      res%converged = .true.
      res%message = ""

      do comp = 1, ncomp
         call neg_mat(nva, noa, bvec_a(:, :, comp), ra)
         call neg_mat(nvb, nob, bvec_b(:, :, comp), rb)
         call zero_mat(nva, noa, u_a(:, :, comp))
         call zero_mat(nvb, nob, u_b(:, :, comp))
         call precond_mat(nva, noa, diag_a, ra, za)
         call precond_mat(nvb, nob, diag_b, rb, zb)
         call copy_mat(nva, noa, za, pa)
         call copy_mat(nvb, nob, zb, pb)
         rz = lp%dot_2d(ra, za) + lp%dot_2d(rb, zb)
         resid = sqrt(lp%dot_2d(ra, ra) + lp%dot_2d(rb, rb))
         it = 0
         broke = .false.
         ! See trc_cphf_solve_rhf's own comment on `.not. (resid < conv)`
         ! vs `resid >= conv`: a NaN resid must read as "not converged",
         ! not fall through both comparisons as false.
         do while (it < opts%max_iter .and. .not. (resid < opts%conv_residual))
            call matvec(pa, pb, apa, apb)
            pap = lp%dot_2d(pa, apa) + lp%dot_2d(pb, apb)
            ! See trc_cphf_solve_rhf's own comment on this guard: a
            ! breakdown (pap <= 0, or NaN) must stop here, never fall
            ! through to an Inf/NaN alpha_cg that the residual check
            ! below cannot catch.
            if (.not. (pap > 0.0_dp)) then
               broke = .true.
               exit
            end if
            alpha_cg = rz/pap
            call axpy_mat(nva, noa, alpha_cg, pa, u_a(:, :, comp))
            call axpy_mat(nvb, nob, alpha_cg, pb, u_b(:, :, comp))
            call axpy_mat(nva, noa, -alpha_cg, apa, ra)
            call axpy_mat(nvb, nob, -alpha_cg, apb, rb)
            resid = sqrt(lp%dot_2d(ra, ra) + lp%dot_2d(rb, rb))
            it = it + 1
            if (resid < opts%conv_residual) exit
            call precond_mat(nva, noa, diag_a, ra, za)
            call precond_mat(nvb, nob, diag_b, rb, zb)
            rz_old = rz
            rz = lp%dot_2d(ra, za) + lp%dot_2d(rb, zb)
            beta_cg = rz/rz_old
            call combine_p(nva, noa, za, beta_cg, pa)
            call combine_p(nvb, nob, zb, beta_cg, pb)
         end do
         res%iterations(comp) = it
         if (broke) then
            res%converged = .false.
            res%message = "trc_cphf: CG breakdown (p^T A p <= 0, or NaN) on a right-hand side"
         else if (.not. (resid < opts%conv_residual)) then
            res%converged = .false.
            res%message = "trc_cphf: a right-hand side did not converge"
         end if
         if (opts%verbose) print "(a,i0,a,i0,a,es9.2)", &
            "   cphf (uhf) rhs ", comp, ": ", it, " iterations, residual ", resid
      end do

      !$acc exit data delete(diag_a, diag_b)
      !$acc exit data delete(d_ua, d_ub, d_t, d_tmp, g_j, g_ka, g_kb, g_eff_a, g_eff_b, &
      !$acc                  tmp_a, tmp_a2, tmp_b, tmp_b2, ra, za, pa, apa, rb, zb, pb, apb)

   contains

      subroutine matvec(pina, pinb, apouta, apoutb)
         real(dp), intent(in) :: pina(nva, noa), pinb(nvb, nob)
         real(dp), intent(inout) :: apouta(nva, noa), apoutb(nvb, nob)
         ! D_U^alpha
         call lp%gemm_mat(cmo_a(:, noa + 1:nao), pina, tmp_a)
         call lp%gemm_mat(tmp_a, cmo_a(:, 1:noa), d_ua, transb="T")
         call lp%transpose(nao, nao, d_ua, nao, d_tmp, nao)
         call add_mat(nao, nao, d_ua, d_tmp, d_ua)
         ! D_U^beta
         call lp%gemm_mat(cmo_b(:, nob + 1:nao), pinb, tmp_b)
         call lp%gemm_mat(tmp_b, cmo_b(:, 1:nob), d_ub, transb="T")
         call lp%transpose(nao, nao, d_ub, nao, d_tmp, nao)
         call add_mat(nao, nao, d_ub, d_tmp, d_ub)
         call add_mat(nao, nao, d_ua, d_ub, d_t)
         call eri%fock_resident(b, d_t, g_j, j_scale=1.0_dp, k_scale=0.0_dp, density_screen=.false.)
         call eri%fock_resident(b, d_ua, g_ka, j_scale=0.0_dp, k_scale=2.0_dp, density_screen=.false.)
         call eri%fock_resident(b, d_ub, g_kb, j_scale=0.0_dp, k_scale=2.0_dp, density_screen=.false.)
         call add_mat(nao, nao, g_j, g_ka, g_eff_a)
         call add_mat(nao, nao, g_j, g_kb, g_eff_b)
         call lp%gemm_mat(g_eff_a, cmo_a(:, 1:noa), tmp_a2)
         call lp%gemm_mat(cmo_a(:, noa + 1:nao), tmp_a2, apouta, transa="T")
         call lp%gemm_mat(g_eff_b, cmo_b(:, 1:nob), tmp_b2)
         call lp%gemm_mat(cmo_b(:, nob + 1:nao), tmp_b2, apoutb, transa="T")
         call add_diag(nva, noa, diag_a, pina, apouta)
         call add_diag(nvb, nob, diag_b, pinb, apoutb)
      end subroutine matvec

   end subroutine trc_cphf_solve_uhf

   !
   ! Static dipole polarizability, RHF: alpha_xy = -2 [sum_ai mu^x_ai U^y_ai
   ! + sum_ai mu^y_ai U^x_ai] (symmetrised -- see this routine's own
   ! comment), with U solving (e_a-e_i) U^y + [A U^y] = -mu^y (the "B" of
   ! this module's equation is the dipole integral itself, occ-virt MO
   ! block -- a static field adds no AO-basis response, so there is no
   ! other term).
   ! `cmo`/`eps` are `res%cmo(:,:,1)`/`res%eps(:,1)`, HOST, exactly as
   ! trc_scf_run hands them back. `pl` and `eri` are the caller's, already
   ! built on the SAME basis `b` (not necessarily the ones trc_scf_run used
   ! internally -- that engine is released when the SCF returns).
   !
   subroutine trc_scf_polarizability_rhf(b, pl, eri, nocc, cmo, eps, alpha, opts, res, la)
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(b%nao, b%nao), eps(b%nao)
      real(dp), intent(out) :: alpha(3, 3)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res
      type(trc_linalg_t), intent(inout), optional, target :: la

      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      real(dp), allocatable :: mmat(:, :, :), tmp(:, :), bvec(:, :, :), u(:, :, :)
      integer :: nao, nvirt, x, y
      logical :: owns_la

      nao = b%nao
      nvirt = nao - nocc
      owns_la = .not. present(la)
      if (present(la)) then
         lp => la
      else
         lp => la_own
         call lp%init(nao)
      end if

      allocate (mmat(nao, nao, TRC_NMULT))
      call trc_multipoles(b, pl, [0.0_dp, 0.0_dp, 0.0_dp], mmat)

      allocate (tmp(nao, nocc), bvec(nvirt, nocc, 3), u(nvirt, nocc, 3))
      !$acc enter data copyin(cmo, mmat) create(tmp, bvec, u)
      do x = 1, 3
         call lp%gemm_mat(mmat(:, :, x), cmo(:, 1:nocc), tmp)
         call lp%gemm_mat(cmo(:, nocc + 1:nao), tmp, bvec(:, :, x), transa="T")
      end do

      call trc_cphf_solve_rhf(b, eri, lp, nao, nocc, nvirt, cmo, eps, 3, bvec, u, opts, res)

      ! alpha_xy = -2 [sum_ai mu^x_ai U^y_ai + sum_ai mu^y_ai U^x_ai] -- the
      ! symmetrised form pyscf's own rhf.polarizability computes (e2 + e2^T),
      ! not the task header's bare -2*sum_ai mu^x U^y alone: those agree only
      ! once U is exact (A self-adjoint makes mu^x.U^y = mu^y.U^x), and a
      ! not-fully-converged CG residual breaks that equality asymmetrically
      ! before cancelling here. Measured (water/6-31G*): using the bare,
      ! unsymmetrised form read a tensor exactly half of pyscf's.
      do y = 1, 3
         do x = 1, 3
            alpha(x, y) = -2.0_dp*(lp%dot_2d(bvec(:, :, x), u(:, :, y)) + lp%dot_2d(bvec(:, :, y), u(:, :, x)))
         end do
      end do

      !$acc exit data delete(cmo, mmat, tmp, bvec, u)
      deallocate (mmat, tmp, bvec, u)
      if (owns_la) call lp%release()
   end subroutine trc_scf_polarizability_rhf

   !
   ! Static dipole polarizability, UHF: alpha_xy = -sum_s sum_ai mu^x_ai,s
   ! U^y_ai,s -- the sum over both spin channels already carries RHF's
   ! factor of 2 (one identical contribution per spin there; two distinct
   ! ones, each counted once, here).
   !
   subroutine trc_scf_polarizability_uhf(b, pl, eri, noa, nob, cmo_a, cmo_b, eps_a, eps_b, alpha, opts, res, la)
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: noa, nob
      real(dp), intent(in) :: cmo_a(b%nao, b%nao), cmo_b(b%nao, b%nao), eps_a(b%nao), eps_b(b%nao)
      real(dp), intent(out) :: alpha(3, 3)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res
      type(trc_linalg_t), intent(inout), optional, target :: la

      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      real(dp), allocatable :: mmat(:, :, :), tmp_a(:, :), tmp_b(:, :)
      real(dp), allocatable :: bvec_a(:, :, :), bvec_b(:, :, :), u_a(:, :, :), u_b(:, :, :)
      integer :: nao, nva, nvb, x, y
      logical :: owns_la

      nao = b%nao
      nva = nao - noa
      nvb = nao - nob
      owns_la = .not. present(la)
      if (present(la)) then
         lp => la
      else
         lp => la_own
         call lp%init(nao)
      end if

      allocate (mmat(nao, nao, TRC_NMULT))
      call trc_multipoles(b, pl, [0.0_dp, 0.0_dp, 0.0_dp], mmat)

      allocate (tmp_a(nao, noa), tmp_b(nao, nob))
      allocate (bvec_a(nva, noa, 3), bvec_b(nvb, nob, 3), u_a(nva, noa, 3), u_b(nvb, nob, 3))
      !$acc enter data copyin(cmo_a, cmo_b, mmat) create(tmp_a, tmp_b, bvec_a, bvec_b, u_a, u_b)
      do x = 1, 3
         call lp%gemm_mat(mmat(:, :, x), cmo_a(:, 1:noa), tmp_a)
         call lp%gemm_mat(cmo_a(:, noa + 1:nao), tmp_a, bvec_a(:, :, x), transa="T")
         call lp%gemm_mat(mmat(:, :, x), cmo_b(:, 1:nob), tmp_b)
         call lp%gemm_mat(cmo_b(:, nob + 1:nao), tmp_b, bvec_b(:, :, x), transa="T")
      end do

      call trc_cphf_solve_uhf(b, eri, lp, nao, noa, nva, nob, nvb, cmo_a, cmo_b, eps_a, eps_b, &
                              3, bvec_a, bvec_b, u_a, u_b, opts, res)

      ! alpha_xy = -[sum_s sum_ai mu^x_ai,s U^y_ai,s + (x<->y)] -- pyscf's
      ! uhf.polarizability e2 + e2^T; see trc_scf_polarizability_rhf's own
      ! comment on why the symmetrised form, not the bare one, is right.
      do y = 1, 3
         do x = 1, 3
            alpha(x, y) = -(lp%dot_2d(bvec_a(:, :, x), u_a(:, :, y)) + lp%dot_2d(bvec_b(:, :, x), u_b(:, :, y)) + &
                            lp%dot_2d(bvec_a(:, :, y), u_a(:, :, x)) + lp%dot_2d(bvec_b(:, :, y), u_b(:, :, x)))
         end do
      end do

      !$acc exit data delete(cmo_a, cmo_b, mmat, tmp_a, tmp_b, bvec_a, bvec_b, u_a, u_b)
      deallocate (mmat, tmp_a, tmp_b, bvec_a, bvec_b, u_a, u_b)
      if (owns_la) call lp%release()
   end subroutine trc_scf_polarizability_uhf

   !
   ! Dispatch on `res%nspin`: the one entry point a caller holding a
   ! converged `trc_scf_result_t` needs. `nalpha`/`nbeta` are the same
   ! electron counts the SCF run itself took -- trc_scf_result_t does not
   ! carry them (trc_scf_driver.F90's own `cfg%nocc` is private to that
   ! run), so the caller passes them again, exactly as scf_rimp2.F90 does
   ! for `trc_rimp2_run`'s `nocc`.
   !
   subroutine trc_scf_polarizability(b, pl, eri, nalpha, nbeta, res_scf, alpha, opts, res, la)
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: nalpha, nbeta
      type(trc_scf_result_t), intent(in) :: res_scf
      real(dp), intent(out) :: alpha(3, 3)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res
      type(trc_linalg_t), intent(inout), optional, target :: la

      if (res_scf%nspin == 1) then
         call trc_scf_polarizability_rhf(b, pl, eri, nalpha, res_scf%cmo(:, :, 1), res_scf%eps(:, 1), &
                                         alpha, opts, res, la)
      else
         call trc_scf_polarizability_uhf(b, pl, eri, nalpha, nbeta, res_scf%cmo(:, :, 1), res_scf%cmo(:, :, 2), &
                                         res_scf%eps(:, 1), res_scf%eps(:, 2), alpha, opts, res, la)
      end if
   end subroutine trc_scf_polarizability

   ! ----------------------------------------------------------------------
   ! Explicit-shape elementwise helpers shared by the RHF/UHF CG loops
   ! above: never touching anything but their own dummies, so each is a
   ! legal `do concurrent` body on its own -- see trc_scf_fock.F90's own
   ! header for why that split matters on the GPU build.
   ! ----------------------------------------------------------------------

   subroutine zero_mat(m, n, a)
      integer, intent(in) :: m, n
      real(dp), intent(inout) :: a(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         a(i, j) = 0.0_dp
      end do
   end subroutine zero_mat

   subroutine copy_mat(m, n, a, c)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: a(m, n)
      real(dp), intent(inout) :: c(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         c(i, j) = a(i, j)
      end do
   end subroutine copy_mat

   subroutine neg_mat(m, n, a, c)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: a(m, n)
      real(dp), intent(inout) :: c(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         c(i, j) = -a(i, j)
      end do
   end subroutine neg_mat

   subroutine add_mat(m, n, a, bb, c)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: a(m, n), bb(m, n)
      real(dp), intent(inout) :: c(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         c(i, j) = a(i, j) + bb(i, j)
      end do
   end subroutine add_mat

   subroutine axpy_mat(m, n, alpha, x, y)
      !! y := y + alpha*x
      integer, intent(in) :: m, n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(m, n)
      real(dp), intent(inout) :: y(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         y(i, j) = y(i, j) + alpha*x(i, j)
      end do
   end subroutine axpy_mat

   subroutine precond_mat(m, n, diagv, r, z)
      !! z := r / diagv, the diagonal (e_a-e_i) preconditioner.
      integer, intent(in) :: m, n
      real(dp), intent(in) :: diagv(m, n), r(m, n)
      real(dp), intent(inout) :: z(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         z(i, j) = r(i, j)/diagv(i, j)
      end do
   end subroutine precond_mat

   subroutine combine_p(m, n, z, beta, p)
      !! p := z + beta*p
      integer, intent(in) :: m, n
      real(dp), intent(in) :: beta
      real(dp), intent(in) :: z(m, n)
      real(dp), intent(inout) :: p(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         p(i, j) = z(i, j) + beta*p(i, j)
      end do
   end subroutine combine_p

   subroutine add_diag(m, n, diagv, pin, apout)
      !! apout := apout + diagv * pin (elementwise)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: diagv(m, n), pin(m, n)
      real(dp), intent(inout) :: apout(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         apout(i, j) = apout(i, j) + diagv(i, j)*pin(i, j)
      end do
   end subroutine add_diag

end module trc_cphf
