!
! The analytic RI-MP2 nuclear gradient over a closed-shell RHF reference.
!
! Two references, E = E_RHF + E_RI-MP2 in both:
!
!   EXACT ERI (the default): what `terco --rimp2` computes -- the auxiliary
!     basis fits the correlation, not the SCF, so the Z-vector's operator,
!     the two reference potentials and the separable two-electron
!     derivative below are the exact four-centre ones;
!   DENSITY FITTED (`dfjk`/`jkaux` present, `terco --df --rimp2`): the SCF
!     fitted its own J and K, so every one of those three moves to the fit
!     the SCF converged with (`dfjk`, over the jkfit basis `jkaux`): the
!     Z-vector is solved against the fitted orbital Hessian in the MO
!     basis (trc_cphf_solve_rhf_fitted over df_rotation_t), the two
!     potentials of a general density are its occupied projections
!     G_fit[X] C_o (df_response_t, trc_scf_df_response both), and the
!     two-electron derivative is
!     trc_scf_grad_df's scf_grad_df_bilinear -- three- and two-centre, no
!     four-centre integral anywhere. metalquicha's `fitted_reference`.
!
! Stocks, Palethorpe and
! Barca, J. Chem. Theory Comput. 2024, 20, 2505, section 3, ported from
! metalquicha (backends/cenzontle/correlation/mqc_czt_ri_mp2_gradient.f90
! with `fit_response = .false.`, its exact-operator mode, and the numpy
! prototype tools/cpu_validation/ri_mp2_gradient.py whose conventions it
! follows: `imat` is the transpose of the paper's L', L''; the response
! operator is the conventional 2 C_v^T (J - K/2)[D] C_o, half the paper's
! A_pqrs).
!
! WHAT IS WHERE
!   trc_rimp2_relax (post_hf)  amplitudes, the unrelaxed doo/dvv, the
!                              Lagrangian L', L'' from (pq|P), and the
!                              fitted 3c/2c derivative terms;
!   this module                the occupied-frozen block, the Z-vector
!                              (trc_cphf's solver, exact operator,
!                              density_screen = .false. throughout), the
!                              relaxed and energy-weighted densities, and
!                              the assembly.
!
! THE ASSEMBLY (MO quantities over every orbital, AO in the Cartesian space
! the kernels run in; C_cart = T C_sph for a spherical basis, exact):
!
!   P   = C dm1mo C^T, the relaxed correlation density (dm1mo: doo + doo^T,
!         dvv + dvv^T, the occupied-frozen block imat/(e_p - e_i), and the
!         Z-vector in the virtual-occupied blocks)
!   dE  = dE_nuc + sum h^x (D_HF + P)                  (trc_1e_gradient)
!         - sum S^x W                                  (same call, W below)
!         + 1/2 d[Tr D G(D)] + d[Tr P G(D)]            (integrals only)
!         + 4 Gamma . (mu nu|P)^x - 2 gamma . (P|Q)^x  (trc_rimp2_relax)
!   W   = W_HF + C (zeta o dm1mo) C^T + 2 P_o G(P) P_o - sym(C imat C^T)
!
! with zeta_pq = (e_p + e_q)/2 except e_i on the occupied-virtual blocks,
! P_o = C_o C_o^T and G = J - K/2. Exact reference: the four-centre term
! comes from eri%gradient, whose single-density convention is
! d[Tr X G(X)]/dR, so the bilinear 1/2 Q(D) + B(D, P) is
! 1/2 [Q(D + P) - Q(P)], two passes. Fitted reference: it is
! d[Tr (D/2 + P) G_fit(D)] directly, one bilinear pass.
!
! SCOPE. Restricted, one rank, every orbital the SCF kept (no functions
! dropped by canonical orthogonalisation: `cmo` must be square and full
! rank), Cartesian or spherical orbital and auxiliary bases, frozen core.
! The caller refuses the rest (trc_capi.F90's gradient_capability_check:
! an unrestricted reference, a functional, PCM).
!
module trc_rimp2_grad
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_1e_gradient
   use trc_eri, only: trc_eri_t, eri_lmax_refusal
   use trc_linalg, only: trc_linalg_t
   use trc_sph, only: sph_mo_to_cart
   use trc_step_table, only: step_table_t
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_cphf_solve_rhf, trc_cphf_solve_rhf_fitted
   use trc_scf_df, only: df_jk_t
   use trc_scf_df_response, only: df_response_t, df_rotation_t
   use trc_scf_grad_df, only: scf_grad_df_bilinear
   use trc_scf_grad, only: scf_grad_nuclear_repulsion
   use trc_rimp2_driver, only: trc_rimp2_options_t
   use trc_rimp2_relax, only: trc_rimp2_relax_t
   implicit none
   private

   public :: trc_rimp2_gradient, trc_rimp2_grad_options_t, trc_rimp2_grad_result_t, trc_rimp2_grad_t

   real(dp), parameter :: PAD_EPS = 1.0e6_dp
   !! The orbital energy of a zero-coefficient pad column: the CPHF solver
   !! takes a square Cartesian (nao, nao) C, and a spherical basis has
   !! nao - nao_sph columns fewer. Their right-hand side is exactly zero, so
   !! their solution is too; the energy only keeps the preconditioner finite
   !! (trc_rimp2_driver's own rimp2_driver_mo_to_cart pads the same way).

   integer, parameter :: ST_ORB = 1, ST_CORR = 2, ST_ZVEC = 3, ST_DENSITY = 4, ST_GRAD4C = 5, ST_GRAD1E = 6, &
                         NSTEP = 6
   character(len=20), parameter :: STEP_NAME(NSTEP) = [character(len=20) :: &
                                   "orbitals, pairs, ERI", "correlation (fitted)", "Z-vector (CPHF)", &
                                   "relaxed densities", "2e gradient", "1e gradient"]

   type :: trc_rimp2_grad_options_t
      !! The RI-MP2 gradient's knobs. Plain values only.
      integer :: nocc = 0
      !! Doubly occupied orbitals of the reference.
      integer :: nfrozen = 0
      !! Core orbitals left out of the correlation (orbitals 1..nfrozen).
      integer :: aux_block = 0
      !! Auxiliary functions per (mu nu|P) block; 0 is the whole basis.
      real(dp) :: pair_threshold = 1.0e-12_dp
      !! The shell-pair screen of the pair list the fitted integrals run
      !! over: the energy's own (trc_capi.F90's trc_run_rimp2 builds it at
      !! the SCF's eri_thresh), so the gradient differentiates that energy.
      real(dp) :: eri_threshold = 1.0e-12_dp
      !! The four-centre engine's Schwarz threshold (Z-vector, potentials,
      !! the separable derivative).
      real(dp) :: cphf_residual = 1.0e-9_dp
      !! The Z-vector's residual norm: it enters the gradient linearly.
      !! Measured on gly10/cc-pVDZ against a 1e-12 solve: the gradient
      !! moves by 2e-11 at 1e-9, 2e-10 at 1e-8 and 2e-9 at 1e-7 (about
      !! 0.05 x the residual), and each decade costs two CG iterations (one
      !! four-centre Fock build each).
      integer :: cphf_max_iteration = 200
      logical :: df_response_ao = .false.
      !! A --df reference's Z-vector through the AO form of the fitted
      !! orbital Hessian even when its MO blocks fit (the memory fallback;
      !! for testing).
      logical :: verbose = .false.
      !! Log the step tables (this module's and the correlation's).
   end type trc_rimp2_grad_options_t

   type :: trc_rimp2_grad_result_t
      real(dp) :: e_corr = 0.0_dp
      !! The correlation energy the gradient's amplitudes give.
      integer :: cphf_iterations = 0
      logical :: ok = .false.
      character(len=256) :: message = ""
      !! Empty when `ok`; the refusal or the failure otherwise.
   end type trc_rimp2_grad_result_t

   type :: trc_rimp2_grad_t
      !! One RI-MP2 gradient: what to compute (`opts`), how it went (`res`)
      !! and, when `res%ok`, the gradient itself.
      type(trc_rimp2_grad_options_t) :: opts
      type(trc_rimp2_grad_result_t) :: res
      real(dp), allocatable :: grad(:, :)
      !! (3, natm), Eh/Bohr: d(E_RHF + E_RI-MP2)/dR.
   end type trc_rimp2_grad_t

   type :: grad_work_t
      !! One gradient's dimensions and its MO/AO arrays. Each array is mapped
      !! on its own, never the object (see trc_rimp2_relax's relax_work_t).
      integer :: n = 0, nmo = 0, nocc = 0, nfz = 0, nv = 0, natm = 0
      real(dp), allocatable :: ccart(:, :)   !! (n, nmo) C in the Cartesian AO space
      real(dp), allocatable :: eps(:)        !! (nmo)
      real(dp), allocatable :: dm1mo(:, :)   !! (nmo, nmo) correlation density, MO
      real(dp), allocatable :: imat(:, :)    !! (nmo, nmo) the Lagrangian in mqc's slot
      real(dp), allocatable :: pao(:, :)     !! (n, n) relaxed correlation density, AO
      real(dp), allocatable :: wao(:, :)     !! (n, n) total energy-weighted density
      real(dp), allocatable :: dhf(:, :)     !! (n, n) the reference density, 2 C_o C_o^T
      real(dp), allocatable :: goo(:, :)     !! (nocc, nocc) C_o^T G[P] C_o, the --df path's own
   end type grad_work_t

contains

   recursive subroutine trc_rimp2_gradient(b, aux, cmo, eps, job, la, dfjk, jkaux)
      !! dE(RHF + RI-MP2)/dR, Eh/Bohr, into `job%grad`, at the converged RHF
      !! orbitals `cmo`/`eps` (the SCF's own space: nao_sph square when
      !! `b%spherical`, nao otherwise; host arrays) with `job%opts`. `aux` is
      !! the RI-MP2 auxiliary basis, on the same atoms. `dfjk` and `jkaux`
      !! present (both or neither): the reference is the --df SCF that fit
      !! converged, see this module's header. `job%res%ok` false: no gradient,
      !! and `job%res%message` says why. `la` absent: made and released here.
      type(trc_basis_t), intent(inout) :: b, aux
      real(dp), intent(in) :: cmo(merge(b%nao_sph, b%nao, b%spherical), merge(b%nao_sph, b%nao, b%spherical))
      real(dp), intent(in) :: eps(merge(b%nao_sph, b%nao, b%spherical))
      type(trc_rimp2_grad_t), intent(inout) :: job
      type(trc_linalg_t), intent(inout), optional :: la
      type(df_jk_t), intent(inout), optional :: dfjk
      type(trc_basis_t), intent(inout), optional :: jkaux

      type(grad_work_t) :: w
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_rimp2_relax_t) :: relax
      type(df_response_t) :: resp
      type(trc_rimp2_options_t) :: ropts
      type(step_table_t) :: tbl
      type(timer_type) :: tc, tall
      character(len=:), allocatable :: why
      logical :: fitted

      if (.not. present(la)) then
         block
            type(trc_linalg_t) :: la_local
            call trc_rimp2_gradient(b, aux, cmo, eps, job, la_local, dfjk, jkaux)
            call la_local%release()
         end block
         return
      end if

      job%res = trc_rimp2_grad_result_t()
      if (allocated(job%grad)) deallocate (job%grad)
      fitted = present(dfjk)
      why = rimp2_grad_refusal(b, aux, job%opts)
      if (len(why) == 0 .and. fitted .neqv. present(jkaux)) why = "a fitted reference needs both its fit and its basis"
      if (len(why) == 0 .and. fitted) why = fitted_reference_refusal(b, jkaux, dfjk)
      if (len(why) > 0) then
         job%res%message = "trc_rimp2_gradient: "//why
         return
      end if
      call tbl%init("rimp2 gradient", STEP_NAME)
      call la%init(b%nao)
      call la%sync()
      call tall%start()

      ! --- orbitals in the Cartesian AO space, the pair list, the engine ---
      call tbl%step_on(la, ST_ORB, tc)
      call b%to_device()
      call aux%to_device()
      if (fitted) call jkaux%to_device()
      call cartesian_orbitals(b, job%opts, cmo, eps, w)
      call pl%build(b, job%opts%pair_threshold)
      call pl%to_device()
      if (.not. fitted) call eri%build(b, job%opts%eri_threshold)
      call tbl%step_off(la, ST_ORB, tc)

      ! --- the correlation: doo, dvv, L', L'', the fitted 3c/2c terms -----
      call tbl%step_on(la, ST_CORR, tc)
      ropts%nfrozen = w%nfz
      ropts%aux_block = job%opts%aux_block
      ropts%verbose = job%opts%verbose
      call relax%init(w%nmo, w%nocc)
      call relax%build(b, aux, pl, w%ccart, w%eps, ropts, la)
      call tbl%step_off(la, ST_CORR, tc)
      if (len_trim(relax%message) > 0) then
         job%res%message = "trc_rimp2_gradient: "//trim(relax%message)
         !$acc exit data delete(w%ccart)
         call release_all(pl, eri, fitted)
         return
      end if
      job%res%e_corr = relax%e_corr

      ! --- the Z-vector, then the relaxed and energy-weighted densities ----
      call tbl%step_on(la, ST_ZVEC, tc)
      call unrelaxed_density(relax, w)
      ! After the correlation, whose B and Gamma it would otherwise sit next to.
      if (fitted) call resp%setup(dfjk, la, w%ccart, w%eps, w%nocc, keep_blocks=.false.)
      if (fitted) then
         call z_vector(b, la, job, w, dfjk=dfjk, resp=resp)
      else
         call z_vector(b, la, job, w, eri=eri)
      end if
      call tbl%step_off(la, ST_ZVEC, tc)
      if (len_trim(job%res%message) > 0) then
         call resp%release()
         call release_all(pl, eri, fitted)
         return
      end if
      call tbl%step_on(la, ST_DENSITY, tc)
      if (fitted) then
         call relaxed_densities(b, la, w, dfjk=dfjk)
      else
         call relaxed_densities(b, la, w, eri=eri)
      end if
      call tbl%step_off(la, ST_DENSITY, tc)

      ! --- the assembly ------------------------------------------------------
      allocate (job%grad(3, b%natm))
      call tbl%step_on(la, ST_GRAD4C, tc)
      if (fitted) then
         call fitted_two_electron(b, jkaux, pl, dfjk, la, w, job%grad)
      else
         call separable_two_electron(b, eri, w, job%grad)
      end if
      call tbl%step_off(la, ST_GRAD4C, tc)
      call tbl%step_on(la, ST_GRAD1E, tc)
      !$acc enter data copyin(job%grad)
      call trc_1e_gradient(b, pl, w%wao, w%dhf, job%grad)
      !$acc update self(job%grad)
      !$acc exit data delete(job%grad, w%ccart, w%pao, w%wao, w%dhf)
      call tbl%step_off(la, ST_GRAD1E, tc)
      call scf_grad_nuclear_repulsion(b, job%grad)
      job%grad = job%grad + relax%grad
      call release_all(pl, eri, fitted)
      job%res%ok = .true.
      if (job%opts%verbose) then
         call la%sync()
         call tbl%report(tall%get_elapsed_time(), 0, 1)
      end if
   end subroutine trc_rimp2_gradient

   function fitted_reference_refusal(b, jkaux, dfjk) result(why)
      !! Why `dfjk` over `jkaux` cannot be this basis's --df reference, or "".
      type(trc_basis_t), intent(in) :: b, jkaux
      type(df_jk_t), intent(in) :: dfjk
      character(len=:), allocatable :: why
      integer :: n, naux, npack
      why = ""
      if (.not. dfjk%is_ready()) then
         why = "the --df reference's fit is not set up"
         return
      end if
      call dfjk%dims(n, naux, npack)
      if (n /= b%nao .or. naux /= merge(jkaux%nao_sph, jkaux%nao, jkaux%spherical)) then
         why = "the --df fit does not match the orbital and jkfit bases"
      else if (jkaux%natm /= b%natm) then
         why = "the jkfit basis is on a different number of atoms than the orbital basis"
      else if (maxval(abs(jkaux%at_r - b%at_r)) > 1.0e-10_dp) then
         why = "the jkfit basis is on a different geometry than the orbital basis"
      else if (any(jkaux%sh_at == 0)) then
         why = "a jkfit shell sits on no atom: it has nowhere to put its force"
      end if
   end function fitted_reference_refusal

   function rimp2_grad_refusal(b, aux, opts) result(why)
      !! Why this gradient cannot run on these inputs, or "".
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_rimp2_grad_options_t), intent(in) :: opts
      character(len=:), allocatable :: why
      integer :: nmo
      why = ""
      nmo = merge(b%nao_sph, b%nao, b%spherical)
      if (opts%nocc < 1 .or. opts%nocc >= nmo) then
         why = "needs at least one occupied and one virtual orbital"
      else if (opts%nfrozen < 0 .or. opts%nfrozen >= opts%nocc) then
         why = "the frozen core must leave at least one occupied orbital to correlate"
      else if (aux%natm /= b%natm) then
         why = "the auxiliary basis is on a different number of atoms than the orbital basis"
      else if (maxval(abs(aux%at_r - b%at_r)) > 1.0e-10_dp) then
         ! Equal counts are not enough: the same atoms in another order put
         ! every auxiliary derivative on the wrong nucleus.
         why = "the auxiliary basis is on a different geometry than the orbital basis"
      else if (any(b%sh_at == 0) .or. any(aux%sh_at == 0)) then
         why = "a shell sits on no atom (a ghost or probe function): it has nowhere to put its force"
      else
         why = eri_lmax_refusal(b)
      end if
   end function rimp2_grad_refusal

   subroutine cartesian_orbitals(b, opts, cmo, eps, w)
      !! w's dimensions, and C in the Cartesian AO space: T C_sph for a
      !! spherical basis (exact: the derivative kernels contract T X T^T),
      !! C itself otherwise. Leaves w%ccart device-resident.
      type(trc_basis_t), intent(in) :: b
      type(trc_rimp2_grad_options_t), intent(in) :: opts
      real(dp), intent(in) :: cmo(merge(b%nao_sph, b%nao, b%spherical), merge(b%nao_sph, b%nao, b%spherical))
      real(dp), intent(in) :: eps(merge(b%nao_sph, b%nao, b%spherical))
      type(grad_work_t), intent(inout) :: w
      real(dp), allocatable :: csph(:, :)
      w%n = b%nao
      w%nmo = merge(b%nao_sph, b%nao, b%spherical)
      w%nocc = opts%nocc
      w%nfz = opts%nfrozen
      w%nv = w%nmo - w%nocc
      w%natm = b%natm
      allocate (w%ccart(w%n, w%nmo), w%eps(w%nmo))
      w%eps = eps
      if (b%spherical) then
         allocate (csph(w%nmo, w%nmo))
         csph = cmo
         !$acc enter data copyin(csph) create(w%ccart)
         call sph_mo_to_cart(b, w%nmo, csph, w%ccart)
         !$acc update self(w%ccart)
         !$acc exit data delete(csph)
         deallocate (csph)
      else
         w%ccart = cmo
         !$acc enter data copyin(w%ccart)
      end if
   end subroutine cartesian_orbitals

   subroutine unrelaxed_density(relax, w)
      !! dm1mo and imat before the Z-vector, host (mqc's layout): the
      !! amplitude blocks symmetrised, imat = -L^T in the columns of the
      !! correlated and virtual orbitals, and the occupied-frozen block
      !! imat(p, i)/(e_p - e_i) -- that rotation's own first-order response,
      !! resolved directly since both orbitals are occupied. It lands before
      !! the Z-vector's right-hand side, which responds to all of dm1mo.
      type(trc_rimp2_relax_t), intent(in) :: relax
      type(grad_work_t), intent(inout) :: w
      integer :: nmo, nocc, nfz, i, p
      nmo = w%nmo
      nocc = w%nocc
      nfz = w%nfz
      allocate (w%dm1mo(nmo, nmo), w%imat(nmo, nmo))
      w%dm1mo = 0.0_dp
      w%dm1mo(nfz + 1:nocc, nfz + 1:nocc) = relax%doo + transpose(relax%doo)
      w%dm1mo(nocc + 1:nmo, nocc + 1:nmo) = relax%dvv + transpose(relax%dvv)
      w%imat = 0.0_dp
      w%imat(:, nfz + 1:nocc) = -transpose(relax%lag_o)
      w%imat(:, nocc + 1:nmo) = -transpose(relax%lag_v)
      do i = nfz + 1, nocc
         do p = 1, nfz
            w%dm1mo(p, i) = w%imat(p, i)/(w%eps(p) - w%eps(i))
            w%dm1mo(i, p) = w%dm1mo(p, i)
         end do
      end do
   end subroutine unrelaxed_density

   subroutine z_vector(b, la, job, w, eri, dfjk, resp)
      !! Solve (e_a - e_i) z_ai + [A z]_ai = -X_ai for the virtual-occupied
      !! block of the relaxed density, X = 2 C_v^T G[P_unrelaxed] C_o +
      !! imat(i, a)^T - imat(a, i) (mqc's right-hand side), and add z (and
      !! z^T) into w%dm1mo. G and A are the reference's own: the exact RHF
      !! operator with `eri` (trc_cphf_solve_rhf, density_screen off), the
      !! fitted one with `dfjk` and `resp` (trc_cphf_solve_rhf_fitted). Sets
      !! job%res%message on non-convergence.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rimp2_grad_t), intent(inout) :: job
      type(grad_work_t), intent(inout) :: w
      type(trc_eri_t), intent(inout), optional :: eri
      type(df_jk_t), intent(inout), optional :: dfjk
      type(df_response_t), intent(inout), optional :: resp
      real(dp), allocatable :: half(:, :), xvo(:, :), cpad(:, :), epad(:), bvec(:, :, :), u(:, :, :)
      type(df_rotation_t) :: rot
      type(trc_cphf_options_t) :: copts
      type(trc_cphf_result_t) :: cres
      integer :: n, nmo, nocc, nv, nvp, i, a

      n = w%n
      nmo = w%nmo
      nocc = w%nocc
      nv = w%nv
      nvp = n - nocc
      allocate (w%pao(n, n), half(n, nocc), xvo(nv, nocc))
      !$acc enter data copyin(w%dm1mo) create(w%pao, half, xvo)
      call ao_from_mo(la, n, nmo, w%ccart, w%dm1mo, w%pao)
      call reference_gco(b, la, w, w%pao, half, eri, dfjk, resp)
      call la%gemm_mat(w%ccart(:, nocc + 1:nmo), half, xvo, transa="T", alpha=2.0_dp)
      if (present(dfjk)) then
         ! C_o^T G[P] C_o, which the energy-weighted density needs: this
         ! unrelaxed part now, the Z-vector's after the solve, so the AO
         ! response is set up once.
         allocate (w%goo(nocc, nocc))
         !$acc enter data create(w%goo)
         call la%gemm_mat(w%ccart(:, 1:nocc), half, w%goo, transa="T")
      end if
      !$acc update self(xvo)
      !$acc exit data delete(w%dm1mo, half, xvo)
      do i = 1, nocc
         do a = 1, nv
            xvo(a, i) = xvo(a, i) + w%imat(i, nocc + a) - w%imat(nocc + a, i)
         end do
      end do

      copts%conv_residual = job%opts%cphf_residual
      copts%max_iteration = job%opts%cphf_max_iteration
      copts%verbose = job%opts%verbose
      if (present(dfjk)) then
         ! The right-hand side was the AO response's last use: its blocks
         ! make room for the fitted orbital Hessian's MO ones.
         call resp%release()
         call rot%setup(dfjk, la, w%ccart, w%eps, w%nocc, job%opts%df_response_ao)
         allocate (bvec(nv, nocc, 1), u(nv, nocc, 1))
         bvec(:, :, 1) = xvo
         !$acc enter data copyin(bvec) create(u)
         call trc_cphf_solve_rhf_fitted(rot, dfjk, la, bvec, u, copts, cres)
         block
            real(dp), allocatable :: goo_z(:, :)
            allocate (goo_z(nocc, nocc))
            !$acc enter data create(goo_z)
            call rot%apply_occ(dfjk, la, u(:, :, 1), goo_z)
            !$acc update self(goo_z, w%goo)
            w%goo = w%goo + goo_z
            !$acc update device(w%goo)
            !$acc exit data delete(goo_z)
         end block
         call rot%release()
      else
         ! The exact solver's square Cartesian C: the nmo real orbitals, then
         ! zero-coefficient pad virtuals (PAD_EPS's own comment).
         allocate (cpad(n, n), epad(n), bvec(nvp, nocc, 1), u(nvp, nocc, 1))
         cpad = 0.0_dp
         cpad(:, 1:nmo) = w%ccart
         epad = PAD_EPS
         epad(1:nmo) = w%eps
         bvec = 0.0_dp
         bvec(1:nv, :, 1) = xvo
         !$acc enter data copyin(cpad, bvec) create(u)
         call trc_cphf_solve_rhf(b, eri, la, n, nocc, nvp, cpad, epad, 1, bvec, u, copts, cres)
         !$acc exit data delete(cpad)
      end if
      !$acc update self(u)
      !$acc exit data delete(bvec, u)
      job%res%cphf_iterations = cres%iterations(1)
      if (.not. cres%converged) then
         job%res%message = "trc_rimp2_gradient: the Z-vector did not converge -- "//trim(cres%message)
         if (allocated(w%goo)) then
            !$acc exit data delete(w%goo)
         end if
         !$acc exit data delete(w%pao)
         !$acc exit data delete(w%ccart)
         return
      end if
      w%dm1mo(nocc + 1:nmo, 1:nocc) = w%dm1mo(nocc + 1:nmo, 1:nocc) + u(1:nv, :, 1)
      w%dm1mo(1:nocc, nocc + 1:nmo) = w%dm1mo(1:nocc, nocc + 1:nmo) + transpose(u(1:nv, :, 1))
   end subroutine z_vector

   subroutine relaxed_densities(b, la, w, eri, dfjk, resp)
      !! w%pao := P = C dm1mo C^T (relaxed), w%dhf := D_HF and w%wao := W,
      !! the total energy-weighted density of this module's header, G the
      !! reference's own operator (`eri`, or `dfjk` with `resp`). On the
      !! device; w%pao stays mapped.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(grad_work_t), intent(inout) :: w
      type(trc_eri_t), intent(inout), optional :: eri
      type(df_jk_t), intent(inout), optional :: dfjk
      type(df_response_t), intent(inout), optional :: resp
      real(dp), allocatable :: zw(:, :), im1(:, :), im1_t(:, :), gp(:, :), half(:, :), small(:, :), work(:, :)
      integer :: n, nmo, nocc, p, q, i

      n = w%n
      nmo = w%nmo
      nocc = w%nocc
      ! The Lagrangian's virtual-occupied block becomes its occupied-virtual
      ! block's transpose, and zeta_pq = (e_p + e_q)/2 except e_i where one
      ! index is occupied and the other virtual.
      w%imat(nocc + 1:nmo, 1:nocc) = transpose(w%imat(1:nocc, nocc + 1:nmo))
      allocate (zw(nmo, nmo), small(nocc, nocc))
      do q = 1, nmo
         do p = 1, nmo
            if (p > nocc .and. q <= nocc) then
               zw(p, q) = w%eps(q)*w%dm1mo(p, q)
            else if (p <= nocc .and. q > nocc) then
               zw(p, q) = w%eps(p)*w%dm1mo(p, q)
            else
               zw(p, q) = 0.5_dp*(w%eps(p) + w%eps(q))*w%dm1mo(p, q)
            end if
         end do
      end do
      ! W_HF = 2 sum_i e_i c_i c_i^T, as an (nocc, nocc) diagonal.
      small = 0.0_dp
      do i = 1, nocc
         small(i, i) = 2.0_dp*w%eps(i)
      end do
      allocate (w%wao(n, n), w%dhf(n, n), im1(n, n), im1_t(n, n), gp(n, n), half(n, nocc), work(n, n))
      !$acc enter data copyin(w%dm1mo, w%imat, zw, small) create(w%wao, w%dhf, im1, im1_t, gp, half, work)
      call ao_from_mo(la, n, nmo, w%ccart, w%dm1mo, w%pao)
      call ao_from_mo(la, n, nmo, w%ccart, w%imat, im1)
      call ao_from_mo(la, n, nmo, w%ccart, zw, w%wao)
      call ao_from_mo(la, n, nocc, w%ccart, small, work)
      ! 2 P_o G(P) P_o: G(P + P^T) is 2 G(P), P symmetric.
      if (present(dfjk)) then
         ! C_o^T G[P] C_o from the Z-vector step (z_vector's own comment).
         call scale_copy(nocc, 2.0_dp, w%goo, small)
         !$acc exit data delete(w%goo)
      else
         call reference_gco(b, la, w, w%pao, half, eri, dfjk, resp)
         call la%gemm_mat(w%ccart(:, 1:nocc), half, small, transa="T", alpha=2.0_dp)
      end if
      call ao_from_mo(la, n, nocc, w%ccart, small, gp)
      call la%transpose(n, n, im1, n, im1_t, n)
      call energy_weighted_sum(n, work, gp, im1, im1_t, w%wao)
      call la%density(w%ccart, nocc, 2.0_dp, w%dhf)
      !$acc exit data delete(w%dm1mo, w%imat, zw, small, im1, im1_t, gp, half, work)
      deallocate (zw, small, im1, im1_t, gp, half, work)
   end subroutine relaxed_densities

   subroutine separable_two_electron(b, eri, w, grad)
      !! grad := 1/2 d[Tr D G(D)] + d[Tr P G(D)], the four-centre integral
      !! derivatives at the reference density D and the relaxed P, as
      !! 1/2 [Q(D + P) - Q(P)] with Q(X) = eri%gradient(X) = d[Tr X G(X)].
      !! Then w%dhf := D + P, the density the one-electron term contracts.
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(grad_work_t), intent(inout) :: w
      real(dp), intent(out) :: grad(3, b%natm)
      real(dp), allocatable :: density(:, :, :), grad_p(:, :)
      allocate (density(w%n, w%n, 1), grad_p(3, b%natm))
      !$acc update self(w%pao, w%dhf)
      density(:, :, 1) = w%pao
      call eri%gradient(b, density, grad_p)
      w%dhf = w%dhf + w%pao
      density(:, :, 1) = w%dhf
      call eri%gradient(b, density, grad)
      grad = 0.5_dp*(grad - grad_p)
      !$acc update device(w%dhf)
      deallocate (density, grad_p)
   end subroutine separable_two_electron

   subroutine fitted_two_electron(b, jkaux, pl, dfjk, la, w, grad)
      !! grad := 1/2 d[Tr D G_fit(D)] + d[Tr P G_fit(D)] = d[Tr (D/2 + P)
      !! G_fit(D)], the --df reference's own fitted two-electron derivative
      !! and its cross term with the relaxed P, in one bilinear pass
      !! (scf_grad_df_bilinear). Then w%dhf := D + P, as
      !! separable_two_electron leaves it.
      type(trc_basis_t), intent(in) :: b, jkaux
      type(trc_pairlist_t), intent(in) :: pl
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: la
      type(grad_work_t), intent(inout) :: w
      real(dp), intent(out) :: grad(3, b%natm)
      real(dp), allocatable :: xmat(:, :)
      !$acc update self(w%pao, w%dhf)
      xmat = 0.5_dp*w%dhf + w%pao
      grad = 0.0_dp
      call scf_grad_df_bilinear(b, jkaux, pl, dfjk, la, w%ccart(:, 1:w%nocc), xmat, grad)
      w%dhf = w%dhf + w%pao
      !$acc update device(w%dhf)
   end subroutine fitted_two_electron

   subroutine reference_gco(b, la, w, d, gco, eri, dfjk, resp)
      !! gco := G[d] C_o, G = J - K/2 for a symmetric, possibly indefinite d,
      !! with the reference's own integrals: exact (`eri`, no density screen
      !! -- d is not the SCF's density -- then the projection) or fitted
      !! (`resp` over `dfjk`, which forms the projection directly). Every
      !! array device-resident.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(grad_work_t), intent(in) :: w
      real(dp), intent(in) :: d(:, :)
      real(dp), intent(inout) :: gco(:, :)
      type(trc_eri_t), intent(inout), optional :: eri
      type(df_jk_t), intent(inout), optional :: dfjk
      type(df_response_t), intent(inout), optional :: resp
      real(dp), allocatable :: g(:, :)
      if (present(resp)) then
         call resp%apply(dfjk, la, d, gco)
      else
         allocate (g(w%n, w%n))
         !$acc enter data create(g)
         call eri%fock_resident(b, d, g, k_scale=1.0_dp, j_scale=1.0_dp, density_screen=.false.)
         call la%gemm_mat(g, w%ccart(:, 1:w%nocc), gco)
         !$acc exit data delete(g)
         deallocate (g)
      end if
   end subroutine reference_gco

   subroutine release_all(pl, eri, fitted)
      !! The pair list, and the four-centre engine when it was built.
      type(trc_pairlist_t), intent(inout) :: pl
      type(trc_eri_t), intent(inout) :: eri
      logical, intent(in) :: fitted
      if (.not. fitted) call eri%release()
      call pl%release()
   end subroutine release_all

   subroutine ao_from_mo(la, n, k, c, m, out)
      !! out := C M C^T, C the first k columns of an (n, k)-leading array;
      !! every operand device-resident.
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: n, k
      real(dp), intent(in) :: c(n, k), m(k, k)
      real(dp), intent(inout) :: out(n, n)
      real(dp), allocatable :: half(:, :)
      allocate (half(n, k))
      !$acc enter data create(half)
      call la%gemm("N", "N", n, k, k, 1.0_dp, c, n, m, k, 0.0_dp, half, n)
      call la%gemm("N", "T", n, n, k, 1.0_dp, half, n, c, n, 0.0_dp, out, n)
      !$acc exit data delete(half)
      deallocate (half)
   end subroutine ao_from_mo

   subroutine scale_copy(n, alpha, x, y)
      !! y := alpha x, (n, n), on the device.
      integer, intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(n, n)
      real(dp), intent(out) :: y(n, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:n)
         y(i, j) = alpha*x(i, j)
      end do
   end subroutine scale_copy

   subroutine energy_weighted_sum(n, whf, vso, im1, im1_t, wao)
      !! wao := wao + whf + vso - (im1 + im1^T)/2, elementwise on the device.
      integer, intent(in) :: n
      real(dp), intent(in) :: whf(n, n), vso(n, n), im1(n, n), im1_t(n, n)
      real(dp), intent(inout) :: wao(n, n)
      integer :: mu, nu
      do concurrent(nu=1:n, mu=1:n)
         wao(mu, nu) = wao(mu, nu) + whf(mu, nu) + vso(mu, nu) - 0.5_dp*(im1(mu, nu) + im1_t(mu, nu))
      end do
   end subroutine energy_weighted_sum

end module trc_rimp2_grad
