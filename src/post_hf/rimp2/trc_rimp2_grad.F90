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
!     potentials of the relaxed density are projections through the same
!     MO blocks (trc_scf_df_response), and the
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
! A KOHN-SHAM REFERENCE: THE DOUBLE HYBRID (`opts%func` named). E =
! E_KS[D] + c E_RI-MP2, the RI-MP2 on the converged KS orbitals and c the
! functional's PT2 fraction (`opts%pt2_scale`; one for a plain KS + RI-MP2).
! The KS energy is variational and the PT2 term is not, so the assembly is
! the one above with the reference's operator swapped and one term added
! (metalquicha's czt_ri_mp2_gradient with `xc` present, whose four
! consequences are these):
!
!   1. every response is the KS one: G(X) = J[X] - (a_x/2) K[X] + f_xc[X],
!      a_x the functional's exact exchange -- the Z-vector's operator
!      (trc_cpks_solve_rks exact, trc_cpks_fitted_t over a --df fit), its
!      right-hand side and the 2 P_o G(P) P_o of W;
!   2. the two-electron derivative is the hybrid's, K at a_x, and the
!      reference's own XC gradient at D joins it (trc_xc_grad, grid
!      response included: the SCF's grid moved with the atoms);
!   3. d/dR Tr(P V_xc[D]) at fixed P and D appears (trc_xc_grad_potential):
!      P meets the derivative of the reference operator, and a KS operator
!      contains V_xc[D];
!   4. c scales every correlation quantity once, where they are made
!      (doo, dvv and the Lagrangian, and the fitted 3c/2c terms): the
!      Z-vector is linear in its right-hand side, so P, the correlation's
!      share of W and everything downstream carry c from there.
!
! W_HF becomes the KS orbitals' own 2 sum_i e_i c_i c_i^T, unchanged in
! form. The grid is the SCF's (dft_grid_build at `opts%grid_level`, the
! same points the converged energy used), built here; the XC parts run in
! the Cartesian AO space, as all of terco's XC does.
!
! SCOPE. Restricted, one rank, every orbital the SCF kept (no functions
! dropped by canonical orthogonalisation: `cmo` must be square and full
! rank), Cartesian or spherical orbital and auxiliary bases, frozen core.
! A KS reference: LDA or GGA, global hybrid or pure (trc_cpks's
! cpks_refuse_reason refuses meta-GGA, range separation and VV10). The
! caller refuses the rest (trc_capi.F90's gradient_capability_check: an
! unrestricted reference, PCM, a spin-component-scaled double hybrid).
!
module trc_rimp2_grad
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_1e_gradient
   use trc_eri, only: trc_eri_t, trc_fock_options_t, eri_lmax_refusal
   use trc_decontract, only: trc_decontract_options_t, decontract_options_for_exchange
   use trc_linalg, only: trc_linalg_t
   use trc_sph, only: sph_mo_to_cart
   use trc_step_table, only: step_table_t
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_cphf_solve_rhf, trc_cphf_solve_rhf_fitted
   use trc_cpks, only: trc_cpks_solve_rks, trc_cpks_fitted_t, cpks_refuse_reason
   use trc_scf_df, only: df_jk_t
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc, only: trc_xc_fxc_rks
   use trc_xc_grad, only: xc_grad_add_on_grid, xc_grad_options, xc_grad_options_t, xc_grad_step_names, &
                          XC_GRAD_RHO_TOL
   use trc_xc_grad_potential, only: xc_grad_potential_add
   use trc_dft_grid, only: dft_grid_t, dft_grid_build, DEFAULT_GRID_LEVEL
   use trc_error, only: error_t
   use trc_scf_grad_df, only: scf_grad_df_bilinear, scf_grad_df_bilinear_t
   use trc_scf_grad, only: scf_grad_nuclear_repulsion
   use trc_rimp2_driver, only: trc_rimp2_options_t
   use trc_rimp2_relax, only: trc_rimp2_relax_t
   implicit none
   private

   public :: trc_rimp2_gradient, trc_rimp2_grad_options_t, trc_rimp2_grad_result_t, trc_rimp2_grad_t
   public :: rimp2_grad_reference_refusal

   real(dp), parameter :: PAD_EPS = 1.0e6_dp
   !! The orbital energy of a zero-coefficient pad column: the CPHF solver
   !! takes a square Cartesian (nao, nao) C, and a spherical basis has
   !! nao - nao_sph columns fewer. Their right-hand side is exactly zero, so
   !! their solution is too; the energy only keeps the preconditioner finite
   !! (trc_rimp2_driver's own rimp2_driver_mo_to_cart pads the same way).

   integer, parameter :: ST_ORB = 1, ST_CORR = 2, ST_ZVEC = 3, ST_DENSITY = 4, ST_GRAD4C = 5, ST_GRAD1E = 6, &
                         ST_XC = 7, ST_KSGRID = 8, NSTEP = 8
   character(len=20), parameter :: STEP_NAME(NSTEP) = [character(len=20) :: &
                                   "orbitals, pairs, ERI", "correlation (fitted)", "Z-vector (CPHF)", &
                                   "relaxed densities", "2e gradient", "1e gradient", "XC (KS reference)", &
                                   "KS grid (rebuilt)"]

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
      type(trc_decontract_options_t) :: decontract
      !! The reference SCF's split policy (`trc_scf_options_t%decontract`),
      !! resolved here for the reference's exchange, so the four-centre
      !! engine runs on the SCF's working basis.
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
      type(trc_xc_functional_t) :: func
      !! The Kohn-Sham reference's functional (a double hybrid's own: its
      !! exchange fraction and semilocal kernels); an empty name is the
      !! Hartree-Fock reference. See this module's header.
      real(dp) :: pt2_scale = 1.0_dp
      !! c in E = E_ref + c E_RI-MP2: a double hybrid's PT2 fraction
      !! (c_os = c_ss), one otherwise.
      integer :: grid_level = DEFAULT_GRID_LEVEL
      !! The KS reference's grid: the level, points per batch and density
      !! screen its SCF converged with.
      integer :: grid_max_pts = 512
      real(dp) :: rho_tol = XC_GRAD_RHO_TOL
   end type trc_rimp2_grad_options_t

   type :: trc_rimp2_grad_result_t
      real(dp) :: e_corr = 0.0_dp
      !! The correlation energy the gradient's amplitudes give.
      real(dp) :: e_os = 0.0_dp
      !! Its opposite-spin part, as trc_rimp2_run reports it.
      real(dp) :: e_ss = 0.0_dp
      !! Its same-spin part, as trc_rimp2_run reports it.
      logical :: has_energy = .false.
      !! The amplitudes ran, so e_corr, e_os and e_ss are set -- true even
      !! when a later step (the Z-vector) fails: a caller wanting E_RI-MP2
      !! too needs no separate trc_rimp2_run.
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
      real(dp), allocatable :: scratch(:)    !! the Z-vector's B_vv storage, device-mapped, for the 2e term
      logical :: kohn_sham = .false.
      !! A KS reference (`opts%func` named): the fields below are set.
      real(dp) :: exx = 1.0_dp
      !! a_x of G = J - (a_x/2) K: one for Hartree-Fock.
      real(dp) :: scale = 1.0_dp
      !! opts%pt2_scale.
      type(trc_xc_functional_t) :: func
      type(dft_grid_t) :: grid               !! the KS reference's grid
      type(trc_xc_grid_t) :: xg              !! the same, batched, on the device
      real(dp), allocatable :: dref(:, :)    !! (n, n) the reference density D, host, kept past w%dhf's D + P
   end type grad_work_t

contains

   pure function rimp2_grad_reference_refusal(have_func, func, nspin, have_pcm) result(refusal)
      !! Why the RI-MP2/double-hybrid gradient cannot run over this
      !! reference, or "" when it can: the part of the scope above that is
      !! known before any SCF runs (the functional, the spin and PCM).
      !! trc_capi.F90's gradient_capability_check asks this first and then
      !! adds what only a finished SCF knows (the --df fit, the ranks, the
      !! dropped functions); the CLI asks it before --opt/--hessian starts,
      !! so a run DL-FIND cannot finish is refused with the SAME message
      !! the gradient itself would give at the first step.
      logical, intent(in) :: have_func
      type(trc_xc_functional_t), intent(in) :: func
         !! read only when `have_func`
      integer, intent(in) :: nspin
         !! 1 restricted (nalpha == nbeta), 2 unrestricted
      logical, intent(in) :: have_pcm
      character(len=:), allocatable :: refusal
      refusal = ""
      ! PCM first: the Z-vector and densities carry no PCM response, so a
      ! PCM reference is refused here even where the SCF-level PCM gradient
      ! is supported.
      if (have_pcm) then
         refusal = "the RI-MP2/double-hybrid gradient over a PCM reference is not implemented "// &
                   "(no PCM response in the Z-vector)"
      else if (have_func) then
         if (func%pt2_cos /= func%pt2_css) then
            refusal = "a spin-component-scaled double hybrid ('"//trim(func%name)//"', c_os /= "// &
                      "c_ss) has no gradient yet -- the RI-MP2 gradient scales one correlation energy"
         else if (len_trim(cpks_refuse_reason(func)) > 0) then
            refusal = "the RI-MP2 gradient over this Kohn-Sham reference needs its CPKS response: "// &
                      trim(cpks_refuse_reason(func))
         end if
      end if
      if (len(refusal) > 0) return
      if (nspin /= 1 .and. have_func) then
         refusal = "the double-hybrid (KS + RI-MP2) gradient is restricted only; an unrestricted (UKS) "// &
                   "reference has no correlation gradient -- trc_cpks has no spin-resolved kernel response"
      else if (nspin /= 1) then
         refusal = "the RI-MP2 gradient is restricted (closed-shell RHF reference) only; "// &
                   "an unrestricted reference has no correlation gradient yet"
      end if
   end function rimp2_grad_reference_refusal

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
      if (len(why) == 0 .and. (fitted .neqv. present(jkaux))) why = "a fitted reference needs both its fit and its basis"
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
      if (.not. fitted) then
         call eri%build(b, job%opts%eri_threshold, &
                        decontract=decontract_options_for_exchange(job%opts%decontract, w%exx, 0.0_dp))
      end if
      call tbl%step_off(la, ST_ORB, tc)
      if (w%kohn_sham) then
         call tbl%step_on(la, ST_KSGRID, tc)
         call kohn_sham_setup(b, la, job%opts, w)
         call tbl%step_off(la, ST_KSGRID, tc)
      end if

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
         call release_all(pl, eri, fitted, w)
         return
      end if
      job%res%e_corr = relax%e_corr
      job%res%e_os = relax%e_os
      job%res%e_ss = relax%e_ss
      job%res%has_energy = .true.

      ! --- the Z-vector, then the relaxed and energy-weighted densities ----
      call tbl%step_on(la, ST_ZVEC, tc)
      call unrelaxed_density(relax, w)
      if (fitted) then
         call z_vector(b, la, job, w, dfjk=dfjk)
      else
         call z_vector(b, la, job, w, eri=eri)
      end if
      call tbl%step_off(la, ST_ZVEC, tc)
      if (len_trim(job%res%message) > 0) then
         call release_all(pl, eri, fitted, w)
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
      call tbl%step_off(la, ST_GRAD1E, tc)
      if (w%kohn_sham) then
         call tbl%step_on(la, ST_XC, tc)
         call kohn_sham_xc_gradient(b, la, job%opts, w, job%grad)
         call tbl%step_off(la, ST_XC, tc)
      end if
      !$acc update self(job%grad)
      !$acc exit data delete(job%grad, w%ccart, w%pao, w%wao, w%dhf)
      call scf_grad_nuclear_repulsion(b, job%grad)
      job%grad = job%grad + w%scale*relax%grad
      call release_all(pl, eri, fitted, w)
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
      if (len(why) == 0 .and. len_trim(opts%func%name) > 0) why = trim(cpks_refuse_reason(opts%func))
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
      w%kohn_sham = len_trim(opts%func%name) > 0
      w%scale = opts%pt2_scale
      if (w%kohn_sham) then
         w%func = opts%func
         w%exx = opts%func%exx
      end if
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
      !! Everything is scaled by w%scale here, once (this module's header).
      type(trc_rimp2_relax_t), intent(in) :: relax
      type(grad_work_t), intent(inout) :: w
      integer :: nmo, nocc, nfz, i, p
      nmo = w%nmo
      nocc = w%nocc
      nfz = w%nfz
      allocate (w%dm1mo(nmo, nmo), w%imat(nmo, nmo))
      w%dm1mo = 0.0_dp
      w%dm1mo(nfz + 1:nocc, nfz + 1:nocc) = w%scale*(relax%doo + transpose(relax%doo))
      w%dm1mo(nocc + 1:nmo, nocc + 1:nmo) = w%scale*(relax%dvv + transpose(relax%dvv))
      w%imat = 0.0_dp
      w%imat(:, nfz + 1:nocc) = -w%scale*transpose(relax%lag_o)
      w%imat(:, nocc + 1:nmo) = -w%scale*transpose(relax%lag_v)
      do i = nfz + 1, nocc
         do p = 1, nfz
            w%dm1mo(p, i) = w%imat(p, i)/(w%eps(p) - w%eps(i))
            w%dm1mo(i, p) = w%dm1mo(p, i)
         end do
      end do
   end subroutine unrelaxed_density

   subroutine z_vector(b, la, job, w, eri, dfjk)
      !! Solve (e_a - e_i) z_ai + [A z]_ai = -X_ai for the virtual-occupied
      !! block of the relaxed density, X = 2 C_v^T G[P_unrelaxed] C_o +
      !! imat(i, a)^T - imat(a, i) (mqc's right-hand side), and add z (and
      !! z^T) into w%dm1mo. G and A are the reference's own: the exact RHF
      !! operator with `eri` (trc_cphf_solve_rhf, density_screen off), the
      !! fitted one with `dfjk` (trc_cphf_solve_rhf_fitted over its MO
      !! blocks, df_rotation_t). A KS reference's are its KS operator:
      !! trc_cpks_solve_rks exact, trc_cpks_fitted_t fitted (the same MO
      !! blocks at K scale a_x, plus the XC kernel on the grid). Sets
      !! job%res%message on non-convergence.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rimp2_grad_t), intent(inout) :: job
      type(grad_work_t), intent(inout) :: w
      type(trc_eri_t), intent(inout), optional :: eri
      type(df_jk_t), intent(inout), optional :: dfjk
      real(dp), allocatable :: half(:, :), xvo(:, :), cpad(:, :), epad(:), bvec(:, :, :), u(:, :, :)
      type(trc_cpks_fitted_t) :: cpks
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
      if (present(dfjk)) then
         ! The fitted orbital Hessian's MO blocks serve the right-hand side,
         ! the solve and C_o^T G[P] C_o (which the energy-weighted density
         ! needs: this unrelaxed part now, the Z-vector's after the solve).
         call cpks%rot%setup(dfjk, la, w%ccart, w%eps, w%nocc, job%opts%df_response_ao)
         cpks%rot%k_scale = w%exx
         allocate (w%goo(nocc, nocc))
         !$acc enter data create(w%goo)
         call cpks%rot%apply_general(dfjk, la, w%dm1mo, w%pao, xvo, w%goo)
         if (w%kohn_sham) then
            call cpks%kernel%setup(b, w%xg, w%func, w%dref, w%ccart, nocc, la)
            call cpks%kernel%add_general(b, w%xg, la, w%pao, xvo, w%goo)
         end if
         call scale_copy_rect(nv, nocc, 2.0_dp, xvo)
      else
         call reference_gco(b, la, w, w%pao, half, eri)
         call la%gemm_mat(w%ccart(:, nocc + 1:nmo), half, xvo, transa="T", alpha=2.0_dp)
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
         allocate (bvec(nv, nocc, 1), u(nv, nocc, 1))
         bvec(:, :, 1) = xvo
         !$acc enter data copyin(bvec) create(u)
         if (w%kohn_sham) then
            cpks%opts = copts
            call cpks%solve(b, w%xg, dfjk, la, bvec, u)
            cres%converged = cpks%res%converged
            cres%message = cpks%res%message
            cres%iterations = cpks%res%iterations
         else
            call trc_cphf_solve_rhf_fitted(cpks%rot, dfjk, la, bvec, u, copts, cres)
         end if
         block
            real(dp), allocatable :: goo_z(:, :)
            allocate (goo_z(nocc, nocc))
            !$acc enter data create(goo_z)
            call cpks%rot%apply_occ(dfjk, la, u(:, :, 1), goo_z)
            if (w%kohn_sham) call cpks%kernel%add_rotation(b, w%xg, la, u(:, :, 1), goo=goo_z)
            !$acc update self(goo_z, w%goo)
            w%goo = w%goo + goo_z
            !$acc update device(w%goo)
            !$acc exit data delete(goo_z)
         end block
         call cpks%rot%release_keeping(w%scratch)
         call cpks%kernel%release()
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
         if (w%kohn_sham) then
            call trc_cpks_solve_rks(b, eri, la, w%xg, w%func, n, nocc, nvp, cpad, epad, w%dref, 1, bvec, u, &
                                    copts, cres)
         else
            call trc_cphf_solve_rhf(b, eri, la, n, nocc, nvp, cpad, epad, 1, bvec, u, copts, cres)
         end if
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
         call cpks%release()
         if (allocated(w%scratch)) then
            !$acc exit data delete(w%scratch)
            deallocate (w%scratch)
         end if
         !$acc exit data delete(w%pao)
         !$acc exit data delete(w%ccart)
         return
      end if
      w%dm1mo(nocc + 1:nmo, 1:nocc) = w%dm1mo(nocc + 1:nmo, 1:nocc) + u(1:nv, :, 1)
      w%dm1mo(1:nocc, nocc + 1:nmo) = w%dm1mo(1:nocc, nocc + 1:nmo) + transpose(u(1:nv, :, 1))
   end subroutine z_vector

   subroutine relaxed_densities(b, la, w, eri, dfjk)
      !! w%pao := P = C dm1mo C^T (relaxed), w%dhf := D_HF and w%wao := W,
      !! the total energy-weighted density of this module's header, G the
      !! reference's own operator (`eri`, or with `dfjk` the C_o^T G[P] C_o
      !! z_vector left in w%goo). On the device; w%pao stays mapped.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(grad_work_t), intent(inout) :: w
      type(trc_eri_t), intent(inout), optional :: eri
      type(df_jk_t), intent(inout), optional :: dfjk
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
         call reference_gco(b, la, w, w%pao, half, eri)
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
      !! 1/2 [Q(D + P) - Q(P)] with Q(X) = eri%gradient(X) = d[Tr X G(X)], G at
      !! the reference's own K scale (a_x for a KS reference, one for HF).
      !! Then w%dhf := D + P, the density the one-electron term contracts.
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(grad_work_t), intent(inout) :: w
      real(dp), intent(out) :: grad(3, b%natm)
      real(dp), allocatable :: density(:, :, :), grad_p(:, :)
      allocate (density(w%n, w%n, 1), grad_p(3, b%natm))
      !$acc update self(w%pao, w%dhf)
      density(:, :, 1) = w%pao
      call eri%gradient(b, density, grad_p, k_scale=w%exx)
      w%dhf = w%dhf + w%pao
      density(:, :, 1) = w%dhf
      call eri%gradient(b, density, grad, k_scale=w%exx)
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
      type(scf_grad_df_bilinear_t) :: bil
      !$acc update self(w%pao, w%dhf)
      xmat = 0.5_dp*w%dhf + w%pao
      bil%k_scale = w%exx
      if (allocated(w%scratch)) call move_alloc(w%scratch, bil%buffer)
      call scf_grad_df_bilinear(b, jkaux, pl, dfjk, la, w%ccart(:, 1:w%nocc), xmat, bil)
      grad = bil%grad
      if (allocated(bil%buffer)) then
         !$acc exit data delete(bil%buffer)
         deallocate (bil%buffer)
      end if
      w%dhf = w%dhf + w%pao
      !$acc update device(w%dhf)
   end subroutine fitted_two_electron

   subroutine reference_gco(b, la, w, d, gco, eri)
      !! gco := G[d] C_o for a symmetric, possibly indefinite d, with the
      !! exact integrals (no density screen: d is not the SCF's density): G =
      !! J - K/2, or a KS reference's J - (a_x/2) K + f_xc (f_xc at D, on the
      !! reference's grid). Every array device-resident.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(grad_work_t), intent(in) :: w
      real(dp), intent(in) :: d(:, :)
      real(dp), intent(inout) :: gco(:, :)
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: g(:, :), vresp(:, :)
      allocate (g(w%n, w%n))
      !$acc enter data create(g)
      call eri%fock_resident(b, d, g, opts=trc_fock_options_t(k_scale=w%exx, j_scale=1.0_dp, density_screen=.false.))
      if (w%kohn_sham) then
         allocate (vresp(w%n, w%n))
         !$acc enter data create(vresp)
         call trc_xc_fxc_rks(b, w%xg, w%func, w%dref, d, vresp, la=la)
         call add_matrix(w%n, vresp, g)
         !$acc exit data delete(vresp)
         deallocate (vresp)
      end if
      call la%gemm_mat(g, w%ccart(:, 1:w%nocc), gco)
      !$acc exit data delete(g)
      deallocate (g)
   end subroutine reference_gco

   subroutine kohn_sham_setup(b, la, opts, w)
      !! The KS reference's grid -- the SCF's, rebuilt at its level from the
      !! same atoms, so the same points -- batched with its atoms and
      !! quadrature weights (the grid response needs both) and on the device,
      !! and the reference density D = 2 C_o C_o^T in the Cartesian AO space.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rimp2_grad_options_t), intent(in) :: opts
      type(grad_work_t), intent(inout) :: w
      type(error_t) :: err
      integer, allocatable :: at_elem(:)
      allocate (at_elem(b%natm), source=b%at_elem)
      call dft_grid_build(b%at_r, at_elem, w%grid, err, level=opts%grid_level)
      if (err%has_error()) error stop "trc_rimp2_grad: the KS reference's grid: "//trim(err%get_message())
      call w%xg%build_from_grid(w%grid, b, max_pts=opts%grid_max_pts)
      call w%xg%to_device()
      allocate (w%dref(w%n, w%n))
      !$acc enter data create(w%dref)
      call la%density(w%ccart, w%nocc, 2.0_dp, w%dref)
      !$acc update self(w%dref)
      !$acc exit data delete(w%dref)
   end subroutine kohn_sham_setup

   subroutine kohn_sham_xc_gradient(b, la, opts, w, grad)
      !! grad += dE_xc/dR at D (trc_xc_grad, grid response on) + d/dR Tr(P
      !! V_xc[D]) at fixed P and D (trc_xc_grad_potential), P = w%pao (its
      !! host copy is current: the two-electron step updated it). `grad`
      !! device-resident.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rimp2_grad_options_t), intent(in) :: opts
      type(grad_work_t), intent(inout) :: w
      real(dp), intent(inout) :: grad(3, b%natm)
      type(xc_grad_options_t) :: xopts
      type(step_table_t) :: timing
      type(timer_type) :: tm
      real(dp), allocatable :: dens3(:, :, :)
      real(dp) :: t_exc
      xopts = xc_grad_options(grid_level=opts%grid_level, grid_max_pts=opts%grid_max_pts, rho_tol=opts%rho_tol, &
                              grid_response=.true., partition=w%grid%partition)
      call timing%init("xc grad", xc_grad_step_names())
      call la%sync()
      call tm%start()
      allocate (dens3(w%n, w%n, 1))
      dens3(:, :, 1) = w%dref
      !$acc enter data copyin(dens3)
      call xc_grad_add_on_grid(b, la, w%func, dens3, grad, xopts, timing, w%xg)
      !$acc exit data delete(dens3)
      deallocate (dens3)
      call la%sync()
      t_exc = tm%get_elapsed_time()
      call tm%start()
      call xc_grad_potential_add(b, la, w%func, w%xg, w%dref, w%pao, grad, xopts)
      call la%sync()
      if (opts%verbose) then
         call timing%report(t_exc, 0, 1)
         block
            character(len=160) :: logline
            write (logline, "(a,f9.3,a)") "  rimp2 gradient: d/dR Tr(P V_xc[D]) ", tm%get_elapsed_time(), " s"
            call global_logger%performance(trim(logline))
         end block
      end if
   end subroutine kohn_sham_xc_gradient

   subroutine release_all(pl, eri, fitted, w)
      !! The pair list, the four-centre engine when it was built, and a KS
      !! reference's grid.
      type(trc_pairlist_t), intent(inout) :: pl
      type(trc_eri_t), intent(inout) :: eri
      logical, intent(in) :: fitted
      type(grad_work_t), intent(inout) :: w
      if (.not. fitted) call eri%release()
      call pl%release()
      if (w%kohn_sham) then
         call w%xg%release()
         call w%grid%destroy()
      end if
   end subroutine release_all

   subroutine add_matrix(n, x, y)
      !! y := y + x, (n, n), on the device.
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n, n)
      real(dp), intent(inout) :: y(n, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:n)
         y(i, j) = y(i, j) + x(i, j)
      end do
   end subroutine add_matrix

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

   subroutine scale_copy_rect(m, n, alpha, x)
      !! x := alpha x, (m, n), on the device.
      integer, intent(in) :: m, n
      real(dp), intent(in) :: alpha
      real(dp), intent(inout) :: x(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         x(i, j) = alpha*x(i, j)
      end do
   end subroutine scale_copy_rect

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
