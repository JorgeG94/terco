!
! Coupled-perturbed Kohn-Sham: trc_cphf's orbital-Hessian solve extended
! with the exchange-correlation kernel response, for LDA and GGA
! functionals (restricted only here -- see this module's own limitations
! below). Intended, with trc_cphf, as the foundation for double-hybrid
! gradients; this module only adds the one missing term to trc_cphf's
! equation, not a new solver.
!
! THE EQUATION (RKS)
! -------------------
! Same orbital-rotation amplitudes U_ai, same right-hand side convention,
! as trc_cphf.F90's own header, with the orbital Hessian action now
! carrying the XC kernel response alongside J and K:
!
!     D_U = C_vir U C_occ^T + (C_vir U C_occ^T)^T          (symmetric, AO)
!     [A U]_ai = 2 [C_vir^T (J[D_U] - (a_x/2) K[D_U] + f_xc[D_U]) C_occ]_ai
!
! a_x is the functional's own fraction of exact exchange (`func%exx`,
! zero for a pure functional); `eri%fock_resident(..., k_scale=a_x)`
! already folds the 1/2 into its own K, matching trc_scf_fock.F90's
! `build_g` (k_scale = exx there too) -- see trc_eri.F90's own
! `eri_fock` header ("j*J - k*K/2"). f_xc[D_U], the kernel response, is
! `trc_xc_fxc_rks` (src/xc/trc_xc.F90): the second-derivative contraction
! against the trial density D_U, on the SAME grid batches trc_xc_rks
! builds the ground-state potential from.
!
! WHAT THIS DOES NOT DO
! ----------------------
! Unrestricted (UKS): not implemented. The spin-resolved kernel response
! needs v2rho2_aa/ab/bb, v2rhosigma_{a,b}_{aa,ab,bb} and
! v2sigma2_{aa,ab,bb}_{aa,ab,bb} contracted against BOTH trial spin
! densities, none of which `trc_xc_fxc_rks` computes -- a restricted-only
! first version, per this branch's own instructions.
!
! Meta-GGA, range-separated hybrids, VV10: refused outright by
! `cpks_refuse_reason`, called before the CG loop starts. A meta-GGA
! kernel is not translated at all here (`func%missing` is set and this
! never reaches the refusal logic needed for the other two); a
! range-separated hybrid needs a second, short-range-screened K response
! this module does not build; VV10's nonlocal correlation has no second
! derivative here. All three stop with a message naming the gap rather
! than silently running the wrong physics.
!
! SINGLE-RANK ONLY. `trc_xc_fxc_rks` (src/xc/trc_xc.F90) takes no `comm`
! argument -- unlike `trc_xc_rks`/`xc_drive`, which split batches over a
! communicator for the energy path, the XC kernel response here always
! runs the whole grid on whichever rank calls it. Nothing in this module
! checks a communicator's size or refuses one; a caller under MPI with
! more than one rank gets every rank redundantly repeating the same
! solve rather than an error naming the gap. Not refused at runtime, so
! noted here instead -- the same limitation trc_cphf.F90's CPHF solver
! already has, silently.
!
! THE SOLVER
! ----------
! trc_cphf's own preconditioned CG, reused verbatim (its elementwise
! `do concurrent` helpers are exported for exactly this); the only change
! is one extra call inside `matvec`, folding `f_xc[D_U]` into the
! resident Fock-like build before the two back-transforms.
!
module trc_cpks
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_multipoles, TRC_NMULT
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_xc, only: trc_xc_fxc_rks
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, &
                       zero_mat, copy_mat, neg_mat, add_mat, axpy_mat, precond_mat, combine_p, add_diag
   implicit none
   private

   public :: trc_cpks_solve_rks, trc_scf_polarizability_rks, cpks_refuse_reason

contains

   !
   ! Why this functional cannot run through trc_cpks, or an empty string if
   ! it can: meta-GGA (never reaches `ok=.true.` with an empty `missing`
   ! today -- see trc_xc_functional.F90's own table, so this arm is mostly
   ! future-proofing), a range-separated hybrid (`omega`/`exx_lr` nonzero)
   ! or VV10 (`nlc_b`/`nlc_c` nonzero).
   !
   pure function cpks_refuse_reason(func) result(msg)
      type(trc_xc_functional_t), intent(in) :: func
      character(len=256) :: msg
      msg = ""
      if (len_trim(func%missing) > 0) then
         msg = "trc_cpks: functional '"//trim(func%name)//"' is not available: "//trim(func%missing)
      else if (func%omega /= 0.0_dp .or. func%exx_lr /= 0.0_dp) then
         msg = "trc_cpks: range-separated hybrids are refused -- the CPKS response needs a second, " &
               //"short-range-screened exchange term ('"//trim(func%name)//"') this module does not build"
      else if (func%nlc_b /= 0.0_dp .or. func%nlc_c /= 0.0_dp) then
         msg = "trc_cpks: VV10 nonlocal correlation ('"//trim(func%name)//"') has no CPKS response " &
               //"term here -- its second derivative is not implemented"
      end if
   end function cpks_refuse_reason

   !
   ! RKS CPKS: `ncomp` right-hand sides, trc_cphf_solve_rhf's own CG loop
   ! with the XC kernel response folded into `matvec`. `cmo`, `bvec` and
   ! `u` are device-resident (trc_cphf_solve_rhf's own contract); `eps` and
   ! `dmat` (the ground-state SCF density this functional was converged
   ! with) are HOST, `dmat` copied to the device once for the life of this
   ! solve. `xg` is the grid `dmat`'s own XC build used -- matching it is
   ! what makes the quadrature error the same linear operator from one CG
   ! iteration to the next, and the same one pyscf is compared against
   ! when `xg` is built from a converged `trc_scf_result_t`'s own
   ! `grid_r`/`grid_w`.
   !
   subroutine trc_cpks_solve_rks(b, eri, lp, xg, func, nao, nocc, nvirt, cmo, eps, dmat, ncomp, bvec, u, &
                                 opts, res, budget)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      integer, intent(in) :: nao, nocc, nvirt, ncomp
      real(dp), intent(in) :: cmo(nao, nao)
      real(dp), intent(in) :: eps(nao)
      real(dp), intent(in) :: dmat(nao, nao)
      real(dp), intent(in) :: bvec(nvirt, nocc, ncomp)
      real(dp), intent(inout) :: u(nvirt, nocc, ncomp)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res
      integer(int64), intent(in), optional :: budget

      real(dp), allocatable :: d_u(:, :), d_ut(:, :), g_u(:, :), v_resp(:, :)
      real(dp), allocatable :: tmp_no(:, :), tmp_no2(:, :)
      real(dp), allocatable :: r(:, :), z(:, :), p(:, :), ap(:, :), diagv(:, :)
      real(dp) :: rz, rz_old, pap, alpha_cg, beta_cg, resid, ax
      integer :: it, comp, i, a
      logical :: broke
      character(len=256) :: reason

      reason = cpks_refuse_reason(func)
      if (len_trim(reason) > 0) then
         res%converged = .false.
         res%message = reason
         allocate (res%iterations(ncomp))
         res%iterations = 0
         return
      end if
      ax = func%exx

      allocate (d_u(nao, nao), d_ut(nao, nao), g_u(nao, nao), v_resp(nao, nao))
      allocate (tmp_no(nao, nocc), tmp_no2(nao, nocc))
      allocate (r(nvirt, nocc), z(nvirt, nocc), p(nvirt, nocc), ap(nvirt, nocc), diagv(nvirt, nocc))
      do i = 1, nocc
         do a = 1, nvirt
            diagv(a, i) = eps(nocc + a) - eps(i)
         end do
      end do

      !$acc enter data create(d_u, d_ut, g_u, v_resp, tmp_no, tmp_no2, r, z, p, ap)
      !$acc enter data copyin(diagv, dmat)

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
         ! see trc_cphf.F90's trc_cphf_solve_rhf for why a NaN resid must
         ! read as "not yet converged" here, not fall through both
         ! comparisons as false and report a silent converged.
         do while (it < opts%max_iter .and. .not. (resid < opts%conv_residual))
            call matvec(p, ap)
            pap = lp%dot_2d(p, ap)
            ! See trc_cphf.F90's trc_cphf_solve_rhf for why this guard is
            ! needed: a breakdown (pap <= 0, or NaN -- NaN > 0 is false)
            ! must stop here, never fall through to an Inf/NaN alpha_cg
            ! that the residual check below cannot catch either.
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
            res%message = "trc_cpks: CG breakdown (p^T A p <= 0, or NaN) on a right-hand side"
         else if (.not. (resid < opts%conv_residual)) then
            res%converged = .false.
            res%message = "trc_cpks: a right-hand side did not converge"
         end if
         if (opts%verbose) print "(a,i0,a,i0,a,es9.2)", &
            "   cpks rhs ", comp, ": ", it, " iterations, residual ", resid
      end do

      !$acc exit data delete(diagv, dmat)
      !$acc exit data delete(d_u, d_ut, g_u, v_resp, tmp_no, tmp_no2, r, z, p, ap)

   contains

      ! [A p]_ai = 2 [C_vir^T (J[D_p] - (a_x/2) K[D_p] + f_xc[D_p]) C_occ]_ai
      ! + (e_a-e_i) p_ai -- this module's own header.
      subroutine matvec(pin, apout)
         real(dp), intent(in) :: pin(nvirt, nocc)
         real(dp), intent(inout) :: apout(nvirt, nocc)
         call lp%gemm_mat(cmo(:, nocc + 1:nao), pin, tmp_no)
         call lp%gemm_mat(tmp_no, cmo(:, 1:nocc), d_u, transb="T")
         call lp%transpose(nao, nao, d_u, nao, d_ut, nao)
         call add_mat(nao, nao, d_u, d_ut, d_u)
         call eri%fock_resident(b, d_u, g_u, k_scale=ax, j_scale=1.0_dp, density_screen=.false.)
         call trc_xc_fxc_rks(b, xg, func, dmat, d_u, v_resp, budget, lp)
         call add_mat(nao, nao, g_u, v_resp, g_u)
         call lp%gemm_mat(g_u, cmo(:, 1:nocc), tmp_no2)
         call lp%gemm_mat(cmo(:, nocc + 1:nao), tmp_no2, apout, transa="T", alpha=2.0_dp)
         call add_diag(nvirt, nocc, diagv, pin, apout)
      end subroutine matvec

   end subroutine trc_cpks_solve_rks

   !
   ! Static dipole polarizability through trc_cpks_solve_rks: same
   ! construction as trc_cphf's `trc_scf_polarizability_rhf` (its own
   ! comment on the symmetrised alpha_xy applies verbatim), with the
   ! ground-state density and functional threaded through to the solver.
   ! `cmo`/`eps`/`dmat` are `res_scf%cmo(:,:,1)`/`res_scf%eps(:,1)`/
   ! `res_scf%dmat(:,:,1)`, HOST, exactly as trc_scf_run hands them back;
   ! `xg` is built from `res_scf%grid_r`/`grid_w` by the caller, so the
   ! quadrature here is the one the SCF converged against.
   !
   subroutine trc_scf_polarizability_rks(b, pl, eri, xg, func, nocc, cmo, eps, dmat, alpha, opts, res, la, budget)
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(b%nao, b%nao), eps(b%nao), dmat(b%nao, b%nao)
      real(dp), intent(out) :: alpha(3, 3)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res
      type(trc_linalg_t), intent(inout), optional, target :: la
      integer(int64), intent(in), optional :: budget

      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      real(dp), allocatable :: mmat(:, :, :), tmp(:, :), bvec(:, :, :), u(:, :, :)
      integer :: nao, nvirt, x, y
      logical :: owns_la
      character(len=256) :: reason

      reason = cpks_refuse_reason(func)
      if (len_trim(reason) > 0) then
         res%converged = .false.
         res%message = reason
         allocate (res%iterations(3))
         res%iterations = 0
         alpha = 0.0_dp
         return
      end if

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

      call trc_cpks_solve_rks(b, eri, lp, xg, func, nao, nocc, nvirt, cmo, eps, dmat, 3, bvec, u, opts, res, budget)

      do y = 1, 3
         do x = 1, 3
            alpha(x, y) = -2.0_dp*(lp%dot_2d(bvec(:, :, x), u(:, :, y)) + lp%dot_2d(bvec(:, :, y), u(:, :, x)))
         end do
      end do

      !$acc exit data delete(cmo, mmat, tmp, bvec, u)
      deallocate (mmat, tmp, bvec, u)
      if (owns_la) call lp%release()
   end subroutine trc_scf_polarizability_rks

end module trc_cpks
