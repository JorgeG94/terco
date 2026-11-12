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
! `scf_fock_build_g` (k_scale = exx there too) -- see trc_eri.F90's own
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
! RANKS: correct under a distributed `eri`. `eri%fock_resident` reduces
! `gmat` with `sum_ranks` (an allreduce then a broadcast), so every rank
! holds the full G[D_U]; `trc_xc_fxc_rks` takes no communicator and
! computes the full `v_resp` on every rank. `matvec`'s g_u + v_resp is
! therefore full + full and identical across ranks. The only cost is the
! XC response, done redundantly on each rank. trc_cphf's matvec gets the
! same full G from the same call.
!
! THE SOLVER
! ----------
! trc_cphf's own preconditioned CG, copied (not called -- its
! elementwise `do concurrent` helpers are exported and reused directly,
! but the CG loop itself is this module's own separate copy, see
! trc_cpks_solve_rks below); the only change from trc_cphf's loop is one
! extra call inside `matvec`, folding `f_xc[D_U]` into the resident
! Fock-like build before the two back-transforms.
!
! A --df REFERENCE: trc_cpks_fitted_t
! -----------------------------------
! A --df KS SCF is stationary for its FITTED J and K, so its response
! operator is the fitted one: trc_scf_df_response's df_rotation_t (the
! MO-basis blocks of trc_cphf_solve_rhf_fitted) at k_scale = a_x, plus the
! XC kernel, which is a grid quantity and never fitted: trc_cpks_kernel_t
! forms D_U in the AO basis, calls trc_xc_fxc_rks_apply against its own
! `cache` (collocated once, not every CG iteration -- trc_xc_fxc_cache_t's
! own header) and projects the result on the orbitals. The same CG loop
! again; the double-hybrid gradient's
! Z-vector (trc_rimp2_grad.F90) is its caller. The kernel also serves the
! right-hand side and the energy-weighted density there (`add_general`,
! `add_rotation`'s `goo`), so every response of that gradient is the same
! operator.
!
module trc_cpks
   use pic_types, only: dp, int64
   use pic_logger, only: global_logger
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_multipoles, TRC_NMULT
   use trc_eri, only: trc_eri_t, trc_fock_options_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, XC_FAMILY_MGGA
   use trc_xc, only: trc_xc_fxc_rks, trc_xc_fxc_cache_t, trc_xc_fxc_rks_apply
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, &
                       cphf_zero_matrix, cphf_copy_matrix, cphf_neg_matrix, cphf_add_matrix, &
                       cphf_axpy_matrix, cphf_precond_matrix, cphf_combine_p, cphf_add_diag
   use trc_scf_df, only: df_jk_t
   use trc_scf_df_response, only: df_rotation_t
   implicit none
   private

   public :: trc_cpks_solve_rks, trc_scf_polarizability_rks, cpks_refuse_reason
   public :: trc_cpks_kernel_t, trc_cpks_fitted_t

   type :: trc_cpks_kernel_t
      !! The exchange-correlation kernel's share of a restricted Kohn-Sham
      !! orbital Hessian, f_xc[X] (trc_xc_fxc_rks_apply, evaluated against
      !! `cache`, built once from the reference density `dmat`) projected
      !! on the orbitals: what a fitted CPKS solve and a Z-vector's
      !! right-hand side add to the fitted J - (a_x/2) K of
      !! trc_scf_df_response, whose MO-basis blocks cannot carry a grid term.
      !! `cache` collocates the grid once and keeps as much of it resident
      !! as device memory allows (trc_xc.F90's own trc_xc_fxc_cache_t) --
      !! every X of this type's own CG loop reuses it, never recollocating
      !! what fits. Every array is device-resident, mapped on its own
      !! (never the type), and reaches a kernel only as an actual argument.
      type(trc_xc_functional_t) :: func
      integer :: n = 0, nocc = 0, nvir = 0
      integer(int64) :: budget = 0
      !! trc_xc_fxc_rks_apply's chunk cap in doubles; zero is its own
      !! default. Unrelated to how much of the grid `cache` keeps
      !! resident, which is sized from free device memory on its own.
      logical :: ready = .false.
      real(dp), allocatable :: cocc(:, :), cvir(:, :)
      !! (n, nocc), (n, nvir): the orbitals, Cartesian AO.
      real(dp), allocatable :: xden(:, :), xden_t(:, :), vresp(:, :), half(:, :)
      !! (n, n) trial density and its transpose, (n, n) f_xc[X], (n, nocc) f_xc[X] C_o.
      type(trc_xc_fxc_cache_t) :: cache
   contains
      procedure :: setup => cpks_kernel_setup
      procedure :: add_rotation => cpks_kernel_add_rotation
      procedure :: add_general => cpks_kernel_add_general
      procedure :: release => cpks_kernel_release
   end type trc_cpks_kernel_t

   type :: trc_cpks_fitted_t
      !! The restricted Kohn-Sham orbital Hessian of a --df reference, and its
      !! CG solve: the fitted J - (a_x/2) K in the MO basis (`rot`,
      !! trc_scf_df_response's df_rotation_t at k_scale = a_x) plus the
      !! unfitted XC kernel on the grid (`kernel`). The KS analogue of
      !! trc_cphf_solve_rhf_fitted, same CG, same convergence, breakdown and
      !! NaN handling. `opts` is read and `res` written by `solve`.
      type(df_rotation_t) :: rot
      type(trc_cpks_kernel_t) :: kernel
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res
   contains
      procedure :: apply => cpks_fitted_apply
      procedure :: solve => cpks_fitted_solve
      procedure :: release => cpks_fitted_release
   end type trc_cpks_fitted_t

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
      else if (func%family == XC_FAMILY_MGGA) then
         msg = "trc_cpks: meta-GGA ('"//trim(func%name)//"') has no CPKS response term here -- " &
               //"its second derivatives (tau included) are not translated"
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
      real(dp), allocatable :: temporary_no(:, :), temporary_no2(:, :)
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
      allocate (temporary_no(nao, nocc), temporary_no2(nao, nocc))
      allocate (r(nvirt, nocc), z(nvirt, nocc), p(nvirt, nocc), ap(nvirt, nocc), diagv(nvirt, nocc))
      do i = 1, nocc
         do a = 1, nvirt
            diagv(a, i) = eps(nocc + a) - eps(i)
         end do
      end do

      !$acc enter data create(d_u, d_ut, g_u, v_resp, temporary_no, temporary_no2, r, z, p, ap)
      !$acc enter data copyin(diagv, dmat)

      allocate (res%iterations(ncomp))
      res%converged = .true.
      res%message = ""

      do comp = 1, ncomp
         call cphf_neg_matrix(nvirt, nocc, bvec(:, :, comp), r)
         call cphf_zero_matrix(nvirt, nocc, u(:, :, comp))
         call cphf_precond_matrix(nvirt, nocc, diagv, r, z)
         call cphf_copy_matrix(nvirt, nocc, z, p)
         rz = lp%dot_2d(r, z)
         resid = sqrt(lp%dot_2d(r, r))
         it = 0
         broke = .false.
         ! `.not. (resid < conv_residual)`, not `resid >= conv_residual`:
         ! see trc_cphf.F90's trc_cphf_solve_rhf for why a NaN resid must
         ! read as "not yet converged" here, not fall through both
         ! comparisons as false and report a silent converged.
         do while (it < opts%max_iteration .and. .not. (resid < opts%conv_residual))
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
            call cphf_axpy_matrix(nvirt, nocc, alpha_cg, p, u(:, :, comp))
            call cphf_axpy_matrix(nvirt, nocc, -alpha_cg, ap, r)
            resid = sqrt(lp%dot_2d(r, r))
            it = it + 1
            if (resid < opts%conv_residual) exit
            call cphf_precond_matrix(nvirt, nocc, diagv, r, z)
            rz_old = rz
            rz = lp%dot_2d(r, z)
            beta_cg = rz/rz_old
            call cphf_combine_p(nvirt, nocc, z, beta_cg, p)
         end do
         res%iterations(comp) = it
         if (broke) then
            res%converged = .false.
            res%message = "trc_cpks: CG breakdown (p^T A p <= 0, or NaN) on a right-hand side"
         else if (.not. (resid < opts%conv_residual)) then
            res%converged = .false.
            res%message = "trc_cpks: a right-hand side did not converge"
         end if
         if (opts%verbose) then
            block
               character(len=1024) :: logline
               write (logline, "(a,i0,a,i0,a,es9.2)") &
                  "   cpks rhs ", comp, ": ", it, " iterations, residual ", resid
               call global_logger%info(trim(logline))
            end block
         end if
      end do

      !$acc exit data delete(diagv, dmat)
      !$acc exit data delete(d_u, d_ut, g_u, v_resp, temporary_no, temporary_no2, r, z, p, ap)

   contains

      ! [A p]_ai = 2 [C_vir^T (J[D_p] - (a_x/2) K[D_p] + f_xc[D_p]) C_occ]_ai
      ! + (e_a-e_i) p_ai -- this module's own header.
      subroutine matvec(pin, apout)
         real(dp), intent(in) :: pin(nvirt, nocc)
         real(dp), intent(inout) :: apout(nvirt, nocc)
         call lp%gemm_mat(cmo(:, nocc + 1:nao), pin, temporary_no)
         call lp%gemm_mat(temporary_no, cmo(:, 1:nocc), d_u, transb="T")
         call lp%transpose(nao, nao, d_u, nao, d_ut, nao)
         call cphf_add_matrix(nao, nao, d_u, d_ut, d_u)
         call eri%fock_resident(b, d_u, g_u, &
                                opts=trc_fock_options_t(k_scale=ax, j_scale=1.0_dp, density_screen=.false.))
         call trc_xc_fxc_rks(b, xg, func, dmat, d_u, v_resp, budget, lp)
         call cphf_add_matrix(nao, nao, g_u, v_resp, g_u)
         call lp%gemm_mat(g_u, cmo(:, 1:nocc), temporary_no2)
         call lp%gemm_mat(cmo(:, nocc + 1:nao), temporary_no2, apout, transa="T", alpha=2.0_dp)
         call cphf_add_diag(nvirt, nocc, diagv, pin, apout)
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
   recursive subroutine trc_scf_polarizability_rks(b, pl, eri, xg, func, nocc, cmo, eps, dmat, alpha, opts, res, &
                                                   la, budget)
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
      type(trc_linalg_t), intent(inout), optional :: la
      !! The caller's handles; absent, made and released for this call.
      integer(int64), intent(in), optional :: budget

      real(dp), allocatable :: mmat(:, :, :), temporary(:, :), bvec(:, :, :), u(:, :, :)
      integer :: nao, nvirt, x, y
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
      if (.not. present(la)) then
         block
            type(trc_linalg_t) :: la_local
            call la_local%init(nao)
            call trc_scf_polarizability_rks(b, pl, eri, xg, func, nocc, cmo, eps, dmat, alpha, opts, res, la_local, &
                                            budget)
            call la_local%release()
         end block
         return
      end if

      allocate (mmat(nao, nao, TRC_NMULT))
      call trc_multipoles(b, pl, [0.0_dp, 0.0_dp, 0.0_dp], mmat)

      allocate (temporary(nao, nocc), bvec(nvirt, nocc, 3), u(nvirt, nocc, 3))
      !$acc enter data copyin(cmo, mmat) create(temporary, bvec, u)
      do x = 1, 3
         call la%gemm_mat(mmat(:, :, x), cmo(:, 1:nocc), temporary)
         call la%gemm_mat(cmo(:, nocc + 1:nao), temporary, bvec(:, :, x), transa="T")
      end do

      call trc_cpks_solve_rks(b, eri, la, xg, func, nao, nocc, nvirt, cmo, eps, dmat, 3, bvec, u, opts, res, budget)

      do y = 1, 3
         do x = 1, 3
            alpha(x, y) = -2.0_dp*(la%dot_2d(bvec(:, :, x), u(:, :, y)) + la%dot_2d(bvec(:, :, y), u(:, :, x)))
         end do
      end do

      !$acc exit data delete(cmo, mmat, temporary, bvec, u)
      deallocate (mmat, temporary, bvec, u)
   end subroutine trc_scf_polarizability_rks

   subroutine cpks_kernel_setup(this, b, xg, func, dmat, cmo, nocc, la)
      !! Keep the functional and the orbitals `cmo` (n, nmo), occupied
      !! columns first -- host arrays, the Cartesian AO space the grid's
      !! collocation runs in -- on the device, and build `cache` from the
      !! reference density `dmat` (n, n), once, for every CPKS/Z-vector CG
      !! iteration this kernel runs through (trc_xc_fxc_cache_t's own
      !! header). `dmat` need not stay resident after this call: the cache
      !! keeps its own share of it (coefp and, where it fits, the
      !! collocation itself). `la` is the caller's GEMM handles, passed
      !! through to the cache's own build; absent, it makes and releases
      !! its own.
      class(trc_cpks_kernel_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      real(dp), intent(in) :: dmat(:, :), cmo(:, :)
      integer, intent(in) :: nocc
      type(trc_linalg_t), intent(inout), optional :: la
      integer :: n
      call this%release()
      n = size(cmo, 1)
      this%func = func
      this%n = n
      this%nocc = nocc
      this%nvir = size(cmo, 2) - nocc
      allocate (this%cocc(n, nocc), this%cvir(n, this%nvir))
      allocate (this%xden(n, n), this%xden_t(n, n), this%vresp(n, n), this%half(n, nocc))
      this%cocc = cmo(:, 1:nocc)
      this%cvir = cmo(:, nocc + 1:size(cmo, 2))
      !$acc enter data copyin(this%cocc, this%cvir) create(this%xden, this%xden_t, this%vresp, this%half)
      if (present(la)) then
         call this%cache%setup(b, xg, func, dmat, la)
      else
         call this%cache%setup(b, xg, func, dmat)
      end if
      this%ready = .true.
   end subroutine cpks_kernel_setup

   subroutine cpks_kernel_add_rotation(this, b, xg, lp, u, au, goo)
      !! For the rotation D_U = C_v U C_o^T + (C_v U C_o^T)^T: au += 2 C_v^T
      !! f_xc[D_U] C_o (the orbital Hessian's convention, trc_cpks_solve_rks's
      !! matvec) and goo += C_o^T f_xc[D_U] C_o, each when present. u (nvir,
      !! nocc), au and goo device-resident.
      class(trc_cpks_kernel_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: u(:, :)
      real(dp), intent(inout), optional :: au(:, :), goo(:, :)
      if (.not. this%ready) error stop "trc_cpks: kernel used before setup"
      call lp%gemm_mat(this%cvir, u, this%half)
      call lp%gemm_mat(this%half, this%cocc, this%xden, transb="T")
      call lp%transpose(this%n, this%n, this%xden, this%n, this%xden_t, this%n)
      call cphf_add_matrix(this%n, this%n, this%xden, this%xden_t, this%xden)
      call kernel_response(this, b, xg, lp)
      if (present(au)) call lp%gemm_mat(this%cvir, this%half, au, transa="T", alpha=2.0_dp, beta=1.0_dp)
      if (present(goo)) call lp%gemm_mat(this%cocc, this%half, goo, transa="T", beta=1.0_dp)
   end subroutine cpks_kernel_add_rotation

   subroutine cpks_kernel_add_general(this, b, xg, lp, xao, gvo, goo)
      !! gvo += C_v^T f_xc[X] C_o and goo += C_o^T f_xc[X] C_o for a general
      !! symmetric AO X (n, n): a relaxed density. Every array device-resident.
      class(trc_cpks_kernel_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: xao(:, :)
      real(dp), intent(inout) :: gvo(:, :), goo(:, :)
      if (.not. this%ready) error stop "trc_cpks: kernel used before setup"
      call cphf_copy_matrix(this%n, this%n, xao, this%xden)
      call kernel_response(this, b, xg, lp)
      call lp%gemm_mat(this%cvir, this%half, gvo, transa="T", beta=1.0_dp)
      call lp%gemm_mat(this%cocc, this%half, goo, transa="T", beta=1.0_dp)
   end subroutine cpks_kernel_add_general

   subroutine kernel_response(this, b, xg, lp)
      !! this%half := f_xc[this%xden] C_o, through this%vresp -- against
      !! `this%cache` (cpks_kernel_setup), not a fresh collocation: every
      !! CG iteration through this routine reuses it.
      type(trc_cpks_kernel_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_linalg_t), intent(inout) :: lp
      if (this%budget > 0) then
         call trc_xc_fxc_rks_apply(b, xg, this%cache, this%xden, this%vresp, this%budget, lp)
      else
         call trc_xc_fxc_rks_apply(b, xg, this%cache, this%xden, this%vresp, la=lp)
      end if
      call lp%gemm_mat(this%vresp, this%cocc, this%half)
   end subroutine kernel_response

   subroutine cpks_kernel_release(this)
      class(trc_cpks_kernel_t), intent(inout) :: this
      if (this%ready) then
         !$acc exit data delete(this%cocc, this%cvir, this%xden, this%xden_t, this%vresp, this%half)
      end if
      call this%cache%release()
      if (allocated(this%cocc)) deallocate (this%cocc)
      if (allocated(this%cvir)) deallocate (this%cvir)
      if (allocated(this%xden)) deallocate (this%xden)
      if (allocated(this%xden_t)) deallocate (this%xden_t)
      if (allocated(this%vresp)) deallocate (this%vresp)
      if (allocated(this%half)) deallocate (this%half)
      this%ready = .false.
   end subroutine cpks_kernel_release

   subroutine cpks_fitted_apply(this, b, xg, dfjk, lp, u, au)
      !! au := 2 C_v^T (J_fit - (a_x/2) K_fit + f_xc)[D_U] C_o, the orbital
      !! Hessian's two-electron and kernel part (no (e_a - e_i) diagonal).
      !! u and au (nvir, nocc) device-resident.
      class(trc_cpks_fitted_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: u(:, :)
      real(dp), intent(inout) :: au(:, :)
      call this%rot%apply(dfjk, lp, u, au)
      call this%kernel%add_rotation(b, xg, lp, u, au=au)
   end subroutine cpks_fitted_apply

   subroutine cpks_fitted_solve(this, b, xg, dfjk, lp, bvec, u)
      !! (e_a - e_i) u_ai + [A u]_ai = -b_ai for every right-hand side of
      !! `bvec` (nvir, nocc, ncomp), A this type's orbital Hessian, `xg` the
      !! grid the reference converged on: trc_cphf_solve_rhf_fitted's CG
      !! with `apply` as its matvec. bvec and u device-resident; this%res
      !! says how it went (non-convergence, breakdown and NaN reported, never
      !! passed as converged -- trc_cphf_solve_rhf's own comment).
      class(trc_cpks_fitted_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: bvec(:, :, :)
      real(dp), intent(inout) :: u(:, :, :)

      real(dp), allocatable :: r(:, :), z(:, :), p(:, :), ap(:, :), diagv(:, :)
      real(dp) :: rz, rz_old, pap, alpha_cg, beta_cg, resid
      integer :: nocc, nvirt, ncomp, it, comp, i, a
      logical :: broke
      character(len=256) :: reason

      ! Reset field by field: a structure constructor of a type with an
      ! allocatable component is one nvfortran has miscompiled.
      if (allocated(this%res%iterations)) deallocate (this%res%iterations)
      this%res%converged = .true.
      this%res%message = ""
      nvirt = this%rot%nvir
      nocc = this%rot%nocc
      ncomp = size(bvec, 3)
      allocate (this%res%iterations(ncomp))
      this%res%iterations = 0
      reason = cpks_refuse_reason(this%kernel%func)
      if (len_trim(reason) > 0) then
         this%res%converged = .false.
         this%res%message = reason
         return
      end if
      if (size(bvec, 1) /= nvirt .or. size(bvec, 2) /= nocc .or. this%kernel%nocc /= nocc .or. &
          this%kernel%nvir /= nvirt) error stop "trc_cpks: the fitted solve's right-hand side, rotation and "// &
         "kernel disagree on the orbitals"
      allocate (r(nvirt, nocc), z(nvirt, nocc), p(nvirt, nocc), ap(nvirt, nocc), diagv(nvirt, nocc))
      do i = 1, nocc
         do a = 1, nvirt
            diagv(a, i) = this%rot%eps(nocc + a) - this%rot%eps(i)
         end do
      end do
      !$acc enter data create(r, z, p, ap)
      !$acc enter data copyin(diagv)
      this%res%converged = .true.
      this%res%message = ""
      do comp = 1, ncomp
         call cphf_neg_matrix(nvirt, nocc, bvec(:, :, comp), r)
         call cphf_zero_matrix(nvirt, nocc, u(:, :, comp))
         call cphf_precond_matrix(nvirt, nocc, diagv, r, z)
         call cphf_copy_matrix(nvirt, nocc, z, p)
         rz = lp%dot_2d(r, z)
         resid = sqrt(lp%dot_2d(r, r))
         it = 0
         broke = .false.
         ! `.not. (resid < conv)`: a NaN residual must read as not converged
         ! (trc_cphf_solve_rhf's own comment).
         do while (it < this%opts%max_iteration .and. .not. (resid < this%opts%conv_residual))
            call this%apply(b, xg, dfjk, lp, p, ap)
            call cphf_add_diag(nvirt, nocc, diagv, p, ap)
            pap = lp%dot_2d(p, ap)
            if (.not. (pap > 0.0_dp)) then
               broke = .true.
               exit
            end if
            alpha_cg = rz/pap
            call cphf_axpy_matrix(nvirt, nocc, alpha_cg, p, u(:, :, comp))
            call cphf_axpy_matrix(nvirt, nocc, -alpha_cg, ap, r)
            resid = sqrt(lp%dot_2d(r, r))
            it = it + 1
            if (resid < this%opts%conv_residual) exit
            call cphf_precond_matrix(nvirt, nocc, diagv, r, z)
            rz_old = rz
            rz = lp%dot_2d(r, z)
            beta_cg = rz/rz_old
            call cphf_combine_p(nvirt, nocc, z, beta_cg, p)
         end do
         this%res%iterations(comp) = it
         if (broke) then
            this%res%converged = .false.
            this%res%message = "trc_cpks: CG breakdown (p^T A p <= 0, or NaN) on a right-hand side"
         else if (.not. (resid < this%opts%conv_residual)) then
            this%res%converged = .false.
            this%res%message = "trc_cpks: a right-hand side did not converge"
         end if
         if (this%opts%verbose) then
            block
               character(len=1024) :: logline
               write (logline, "(a,i0,a,i0,a,es9.2)") &
                  "   cpks (fitted) rhs ", comp, ": ", it, " iterations, residual ", resid
               call global_logger%info(trim(logline))
            end block
         end if
      end do
      !$acc exit data delete(diagv)
      !$acc exit data delete(r, z, p, ap)
   end subroutine cpks_fitted_solve

   subroutine cpks_fitted_release(this)
      class(trc_cpks_fitted_t), intent(inout) :: this
      call this%rot%release()
      call this%kernel%release()
   end subroutine cpks_fitted_release

end module trc_cpks
