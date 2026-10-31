!
! Second-order SCF (`--soscf`): once Pulay DIIS has brought a restricted
! Hartree-Fock reference close, converge the rest of the way by trust-
! region Newton on the orbital rotations, matrix-free. Ported from
! metalquicha's `mqc_czt_soscf`/`mqc_orbital_rotation` (the step control)
! and `mqc_czt_rhf`'s `run_soscf_phase` (the phase itself), onto terco's
! own device-resident SCF.
!
! SCOPE. Restricted Hartree-Fock only, exactly like the backend this was
! ported from: `mqc_czt_rhf.f90` is a closed-shell driver with no
! `mqc_czt_uks` equivalent, and its own test suite refuses a continuum
! solvent outright (`test_refuses_pcm`). `trc_scf_settings%resolve`
! refuses `--soscf` with an open shell, a functional, `--pcm`, `--df` or
! more than one rank before any of this ever runs, so none of those
! combinations reaches here silently as plain DIIS.
!
! `trc_scf_soscf_run` IS THE PUBLIC ENTRY POINT. `trc_scf_run` is its one
! caller. `scf_soscf_rotation_matrix` is ALSO public, but only so
! check_soscf.F90 can test the matrix exponential on its own, at sizes
! and random rotations the SCF path itself never happens to exercise --
! nothing else in terco calls it.
!
! BUNDLED ARGUMENTS. Three plain-value types below (`trc_soscf_dims_t`,
! `trc_soscf_step_io_t`, `trc_soscf_hvp_scratch_t`) exist only to keep
! every procedure here at the project's 8-argument limit, not to own or
! share anything across calls. `trc_soscf_dims_t` is scalars only, so a
! structure-constructor literal builds it safely; the other two have
! allocatable components and are never built that way (this branch's own
! no-pointer rule, and the nvfortran constructor defect it travels with)
! -- each is a plain `type(...) :: x` declaration, filled by ordinary
! `allocate` statements in the procedure that owns it, and read only
! through intent(in)/(inout) dummies elsewhere. None of the three is ever
! assigned FROM an existing `trc_basis_t`/`trc_eri_t`/`trc_linalg_t`/
! `scf_ctx_t`-family object: those stay separate arguments everywhere
! (the same `b, eri, lp` triple `trc_cphf.F90` passes around), because
! copying one into a bundle's component would deep-copy its allocatable
! device arrays -- for `trc_eri_t` for a large screened pair list, an
! expensive wrong move, and for `mem:separate` as such, a residency bug:
! the copy's arrays would carry a HOST address never entered into the
! device present table.
!
! THE ORBITAL HESSIAN ACTION (`scf_soscf_hessian_apply_rhf`)
! -------------------------------------------------------------
! `(e_a-e_i) p + [A p]`, `trc_cphf.F90`'s own equation (its header derives
! it, citing Pople/Krishnan/Schlegel/Binkley 1979) -- but NOT a call into
! `trc_cphf`, which is Cartesian-AO-only (its own tests force `cart=true`
! at basis load): the trial density here is built in the SCF's own `nao`
! space (spherical, when the basis is), and `eri%fock_resident` only ever
! runs in the Cartesian one. So this module folds the trial density down
! with `scf_fock_to_cart` before the Fock-like build and folds the result
! back up with `scf_fock_to_sph` after it -- the SAME two calls
! `trc_scf_driver`'s own iteration makes for the physical density every
! step; `sph = .false.` makes both a plain copy, so a Cartesian basis
! costs nothing extra.
!
! THE ALGORITHM, ONE NEWTON ITERATION
! ------------------------------------
! 1. SEMICANONICALISE (`scf_soscf_semicanonicalize`): two small diagonalisations,
!    occupied-occupied and virtual-virtual, that make the Fock matrix block
!    diagonal within each -- a change of representation, not a step, and the
!    one `scf_soscf_hessian_apply_rhf`'s Delta (e_a - e_i) needs to mean the
!    Hessian's diagonal rather than `delta_ij F_ab - delta_ab F_ij`.
! 2. DIFFERENTIATE: the gradient is `4 F_ai` (`scf_soscf_gradient`) and the
!    Hessian is `4` times `scf_soscf_hessian_apply_rhf`'s `(e_a-e_i) p +
!    [A p]`. The factor is the one `mqc_czt_soscf`'s own header derives
!    (`HESSIAN_SCALE`); it cancels out of the Newton step and only sets the
!    UNITS the curvature thresholds below are read in.
! 3. STEP (`scf_soscf_newton_step_rhf`): the level-shifted trust-region Newton
!    equations, solved in a small Krylov subspace built from
!    Hessian-vector products -- one Fock build (`scf_soscf_hessian_apply_rhf`)
!    per product, never written down in full.
! 4. BACKTRACK (`trc_scf_soscf_run`): the step is scaled to the trust
!    radius, the orbitals rotated by the exact matrix exponential
!    (`scf_soscf_rotation_matrix`), and the energy rebuilt in full. Accepted
!    only if it fell (or the model's own predicted gain is too small to
!    measure); otherwise the radius halves and the step is retried.
!
! RESIDENCY. Every nao-by-anything matrix here -- the working orbitals,
! the Fock matrices, the Krylov basis and its images, the rotation and its
! scratch -- stays device-resident for the life of `trc_scf_soscf_run`,
! the same contract `trc_cphf`/`trc_scf_driver` keep. The handful of
! numbers that cross the bus are exactly the ones DIIS already crosses:
! dot products (`la%dot_2d`), orbital energies after a semicanonicalisation
! (small, read once per Newton iteration) and the trust-radius scaling
! factor (`maxval(abs(step))`, read once per iteration too -- `step` is
! `nvirt`-by-`nocc`, not the full rotation).
!
module trc_scf_soscf
   use pic_types, only: dp, default_int
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_lapack_interfaces, only: pic_syev
   use pic_blas_interfaces, only: pic_gemv
   use trc_basis_types, only: trc_basis_t
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_cphf, only: cphf_zero_matrix, cphf_copy_matrix, cphf_add_matrix, cphf_axpy_matrix, cphf_add_diag
   use trc_scf_context, only: scf_ctx_t
   use trc_scf_settings, only: scf_settings_t, trc_scf_result_t, scf_energies_t
   use trc_scf_state, only: scf_state_t
   use trc_scf_fock, only: scf_fock_build_g, scf_fock_energy, scf_fock_assemble, scf_fock_commutator_error, &
                           scf_fock_to_cart, scf_fock_to_sph
   implicit none
   private

   public :: trc_scf_soscf_run
   ! Exported for check_soscf.F90's own standalone test of the matrix
   ! exponential's orthogonality (this module's own header, point 3):
   ! nothing else in terco calls it.
   public :: scf_soscf_rotation_matrix
   public :: soscf_trust_after_rejection

   real(dp), parameter :: HESSIAN_SCALE = 4.0_dp
   !! `scf_soscf_hessian_apply_rhf`'s `(e_a-e_i) p + [A p]` to the true
   !! second derivative of the energy; the gradient carries the same
   !! factor (`scf_soscf_gradient`). See this module's own header, point 2.
   real(dp), parameter :: MIN_CURVATURE = 1.0e-3_dp
   !! Smallest curvature the Newton equations are allowed to divide by, in
   !! hartree per radian squared -- metalquicha's own constant, same units.
   real(dp), parameter :: SADDLE_CURVATURE = -1.0e-5_dp
   !! How negative a curvature has to be before a point is called a saddle.
   real(dp), parameter :: MAX_ROTATION = 0.2_dp
   !! Largest rotation angle the trust radius grows back to, in radians.
   real(dp), parameter :: MIN_ROTATION = 1.0e-6_dp
   !! Below this, backtracking gives up: the direction is not a descent
   !! direction and halving the step again will not help.
   real(dp), parameter :: ENERGY_RESOLUTION = 1.0e-12_dp
   !! A predicted energy drop below this is too small to measure, not a
   !! step that failed -- taken without being tested against the actual
   !! energy (`mqc_orbital_rotation`'s own constant and reasoning).
   real(dp), parameter :: TRUST_GROWTH = 1.3_dp
   !! How fast the trust radius recovers after an accepted step.
   integer, parameter :: DEFAULT_SOSCF_SUBSPACE = 20
   !! Krylov vectors the Newton equations are solved in, at most -- each one
   !! a Fock build, the per-Newton-iteration cost this method adds.
   real(dp), parameter :: SUBSPACE_TOLERANCE = 1.0e-3_dp
   !! How far the Newton equations are solved, relative to the gradient
   !! norm: deliberately loose, since the step is checked against the
   !! actual energy by backtracking rather than trusted outright.
   real(dp), parameter :: PRECONDITION_FLOOR = 1.0e-6_dp
   real(dp), parameter :: LINEAR_DEPENDENCE = 1.0e-8_dp
   integer, parameter :: MAX_BACKTRACKS = 12

   type :: trc_soscf_dims_t
      !! Sizes threaded through the Hessian-vector product and the
      !! semicanonicalisation. Scalars only -- see this module's own
      !! header on why that makes a structure-constructor literal safe
      !! here and nowhere else in this module.
      logical :: sph = .false.
      integer :: nao = 0, nao_c = 0, n_mo = 0, nocc = 0, nvirt = 0
   end type trc_soscf_dims_t

   type :: trc_soscf_hvp_scratch_t
      !! `scf_soscf_hessian_apply_rhf`'s own device-resident scratch,
      !! owned by `scf_soscf_newton_step_rhf` for the life of one call.
      real(dp), allocatable :: diagv(:, :)
      real(dp), allocatable :: d_u(:, :), d_ut(:, :), g_u(:, :)
      real(dp), allocatable :: d_u_c(:, :), g_u_c(:, :)
      real(dp), allocatable :: temporary_no(:, :), temporary_no2(:, :)
   end type trc_soscf_hvp_scratch_t

   type :: trc_soscf_step_io_t
      !! The Newton step's own in/out parameters, owned by
      !! `trc_scf_soscf_run` for the life of the phase: `gradient`/`step`
      !! are reused, never reallocated, every Newton iteration.
      real(dp), allocatable :: gradient(:, :)
      real(dp), allocatable :: step(:, :)
      real(dp) :: escape = 0.0_dp, lowest = 0.0_dp, predicted = 0.0_dp
      integer :: products = 0, max_subspace = 0
      logical :: ok = .true.
      !! False on a breakdown inside the Krylov/level-shift solve (a NaN
      !! or non-positive-definite projected Hessian the eigensolver
      !! itself refused, or no direction at all to build a subspace
      !! from): `trc_scf_soscf_run` must stop on this, not fall through
      !! `max_iteration` iterations on a step that was never actually
      !! solved for (AGENTS.md: an iterative solver reports
      !! non-convergence on a breakdown or NaN, never lets it fall
      !! through as "converged").
   end type trc_soscf_step_io_t

contains

   !
   ! The second-order phase: `trc_scf_run`'s own DIIS loop calls this once
   ! the commutator has fallen below `scf_config%soscf_start`, in place of
   ! the rest of its own iterations. The DIIS iteration this took over
   ! after is `scf_result%iterations` on entry (every iteration already
   ! sets it, including the one that triggered the handover); the steps
   ! left are `scf_config%max_iteration` beyond it. Updates `scf_result`
   ! (density, orbitals, orbital energies, `converged`, `iterations`) in
   ! place; `scf_state`/`scf_energies` are the same per-iteration scratch
   ! and energy bundle the DIIS loop itself writes through, reused here
   ! rather than duplicated.
   !
   subroutine trc_scf_soscf_run(scf_basis, scf_context, scf_config, scf_state, scf_result, scf_energies, la, eri)
      type(trc_basis_t), intent(in) :: scf_basis
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(scf_energies_t), intent(inout) :: scf_energies
      type(trc_linalg_t), intent(inout) :: la
      type(trc_eri_t), intent(inout) :: eri

      type(trc_soscf_dims_t) :: dims
      type(trc_soscf_step_io_t) :: io
      integer :: switch_iteration, max_steps, iteration, trial, info, products_total, builds
      real(dp) :: e_cur, e_try, e_old, de, errmax, scaling, trust, one_norm, factor
      logical :: accepted, saddle, converged_now
      character(len=200) :: line

      real(dp), allocatable :: fock_cur(:, :), coeff_work(:, :), coeff_try(:, :), dmat_accepted(:, :)
      real(dp), allocatable :: fock_mo(:, :), energies_mo(:)
      real(dp), allocatable :: step_scaled(:, :)
      real(dp), allocatable :: kappa(:, :), rotation(:, :), scratch_a(:, :), scratch_b(:, :), scratch_c(:, :)
      real(dp), allocatable :: step_host(:, :)

      switch_iteration = scf_result%iterations
      max_steps = scf_config%max_iteration - switch_iteration
      dims = trc_soscf_dims_t(sph=scf_context%sph, nao=scf_context%nao, nao_c=scf_context%nao_c, &
                              n_mo=scf_context%nao - scf_context%n_lindep, nocc=scf_config%nocc(1), &
                              nvirt=scf_context%nao - scf_context%n_lindep - scf_config%nocc(1))
      products_total = 0
      if (dims%nocc < 1 .or. dims%nvirt < 1) then
         scf_result%message = "trc_scf: --soscf needs at least one occupied and one virtual orbital; "// &
                              "there are no rotations to take"
         return
      end if
      if (max_steps < 1) return

      allocate (fock_cur(dims%nao, dims%nao), coeff_work(dims%nao, dims%n_mo), coeff_try(dims%nao, dims%n_mo))
      allocate (dmat_accepted(dims%nao, dims%nao), fock_mo(dims%n_mo, dims%n_mo), energies_mo(dims%n_mo))
      allocate (io%gradient(dims%nvirt, dims%nocc), io%step(dims%nvirt, dims%nocc))
      allocate (step_scaled(dims%nvirt, dims%nocc), step_host(dims%nvirt, dims%nocc))
      allocate (kappa(dims%n_mo, dims%n_mo), rotation(dims%n_mo, dims%n_mo))
      allocate (scratch_a(dims%n_mo, dims%n_mo), scratch_b(dims%n_mo, dims%n_mo), scratch_c(dims%n_mo, dims%n_mo))
      !$acc enter data create(fock_cur, coeff_work, coeff_try, dmat_accepted, fock_mo) &
      !$acc            create(io%gradient, io%step, step_scaled, kappa, rotation, scratch_a, scratch_b, scratch_c)
      !$acc enter data copyin(energies_mo)

      call cphf_copy_matrix(dims%nao, dims%n_mo, scf_result%cmo(:, 1:dims%n_mo, 1), coeff_work)
      call cphf_copy_matrix(dims%nao, dims%nao, scf_result%dmat(:, :, 1), dmat_accepted)

      ! The point the Newton phase starts from, built in full: the Fock
      ! matrix DIIS left behind is extrapolated and possibly level
      ! shifted, and its energy belongs to the density before the last
      ! diagonalisation -- neither is a state this can differentiate
      ! (mqc_czt_rhf's own run_soscf_phase, same reasoning).
      call scf_soscf_fock_rebuild(scf_basis, scf_context, scf_config, scf_state, scf_result, scf_energies, la, eri)
      call cphf_copy_matrix(dims%nao, dims%nao, scf_state%fock(:, :, 1), fock_cur)
      e_cur = scf_energies%total
      e_old = e_cur
      trust = soscf_trust0_env()
      io%max_subspace = DEFAULT_SOSCF_SUBSPACE
      if (scf_config%talk) then
         call global_logger%info("   second-order phase: trust-region Newton on the orbital rotations, "// &
                                 "handing over from DIIS at iteration "//to_string(switch_iteration))
      end if

      do iteration = 1, max_steps
         builds = 0
         call scf_soscf_semicanonicalize(la, dims, fock_cur, coeff_work, fock_mo, energies_mo, info)
         if (info /= 0) then
            scf_result%message = "trc_scf: --soscf's semicanonicalisation did not converge "// &
                                 "(eigensolver info="//to_string(info)//")"
            exit
         end if
         call scf_soscf_gradient(dims%n_mo, dims%nocc, dims%nvirt, fock_mo, io%gradient)

         ! The commutator at the current (accepted) point, the SAME metric
         ! and the same test trc_scf_run's own DIIS loop converges on.
         call cphf_copy_matrix(dims%nao, dims%nao, fock_cur, scf_state%fock(:, :, 1))
         call scf_fock_commutator_error(scf_context, scf_config, scf_state, scf_result, la)
         errmax = sqrt(la%dot_3d(scf_state%errv, scf_state%errv))
         de = abs(e_cur - e_old)

         !$acc update self(energies_mo)
         io%escape = MAX_ROTATION
         call scf_soscf_newton_step_rhf(scf_basis, eri, la, dims, coeff_work, energies_mo, io)
         products_total = products_total + io%products
         builds = builds + io%products
         ! A breakdown inside the Krylov/level-shift solve (NaN, a
         ! non-positive-definite projected Hessian the eigensolver
         ! refused, or no direction at all): stop here, not after
         ! max_iteration iterations on a step this was never actually
         ! solved for -- AGENTS.md's own rule on an iterative solver and
         ! a breakdown or NaN.
         if (.not. io%ok) then
            scf_result%message = "trc_scf: --soscf's Newton step broke down "// &
                                 "(the projected orbital Hessian could not be solved); stopping"
            exit
         end if

         saddle = io%lowest < SADDLE_CURVATURE
         converged_now = de < scf_config%conv_energy .and. errmax < scf_config%conv_diis
         if (converged_now) then
            if (.not. saddle) then
               scf_result%converged = .true.
               scf_result%iterations = switch_iteration + iteration - 1
               exit
            end if
            ! Leaving a saddle is a fresh direction: the old trust radius
            ! records how the PREVIOUS direction behaved and is usually
            ! small by the time a run has settled this close.
            trust = MAX_ROTATION
            call global_logger%warning("   the SCF has stopped at a saddle point; following the "// &
                                       "negative curvature out")
         end if

         call soscf_pull_matrix(dims%nvirt, dims%nocc, io%step, step_host)
         scaling = maxval(abs(step_host))
         ! A valid 1-norm bound on the FULL n_mo x n_mo kappa `soscf_embed_kappa`
         ! builds from a step of this (ov-block) shape, computed from the
         ! (small) ov block alone -- see `soscf_ov_one_norm`'s own header.
         ! `scf_soscf_rotation_matrix` needs this, not `scaling` (the
         ! largest single element, which bounds nothing): the largest
         ! element can be arbitrarily smaller than the operator norm the
         ! scaling-and-squaring Taylor series needs small, which is how
         ! `exp(kappa)` lost orthogonality at larger nvirt/nocc.
         one_norm = soscf_ov_one_norm(dims%nvirt, dims%nocc, step_host)

         accepted = .false.
         do trial = 1, MAX_BACKTRACKS
            ! Always rescaled by the SAME factor that caps the step at the
            ! current trust radius (1 when the step already fits inside
            ! it): unconditional, so a REJECTED trial's step is never
            ! identical to the one just tried. The step/radius comparison
            ! used to gate scaling off entirely while trust still exceeded
            ! the step's own size, which is a step this module's own
            ! history showed could repeat, unchanged, for every one of
            ! MAX_BACKTRACKS trials.
            factor = min(1.0_dp, trust/scaling)
            call soscf_scale_matrix(dims%nvirt, dims%nocc, io%step, factor, step_scaled)
            call soscf_embed_kappa(dims%n_mo, dims%nocc, dims%nvirt, step_scaled, kappa)
            call scf_soscf_rotation_matrix(la, dims%n_mo, factor*one_norm, kappa, rotation, scratch_a, &
                                           scratch_b, scratch_c)
            call la%gemm_mat(coeff_work, rotation, coeff_try)
            call la%density(coeff_try, dims%nocc, scf_config%occ, scf_result%dmat(:, :, 1))
            call scf_soscf_fock_rebuild(scf_basis, scf_context, scf_config, scf_state, scf_result, scf_energies, &
                                       la, eri)
            e_try = scf_energies%total
            builds = builds + 1

            ! `predicted` below its own resolution is a step too small to
            ! measure, not a step that failed -- see ENERGY_RESOLUTION.
            if (e_try < e_cur .or. io%predicted < ENERGY_RESOLUTION) then
               call cphf_copy_matrix(dims%nao, dims%n_mo, coeff_try, coeff_work)
               call cphf_copy_matrix(dims%nao, dims%nao, scf_state%fock(:, :, 1), fock_cur)
               call cphf_copy_matrix(dims%nao, dims%nao, scf_result%dmat(:, :, 1), dmat_accepted)
               e_old = e_cur
               e_cur = e_try
               trust = min(MAX_ROTATION, trust*TRUST_GROWTH)
               accepted = .true.
               exit
            end if
            trust = soscf_trust_after_rejection(trust, scaling)
            if (trust < MIN_ROTATION) exit
         end do

         if (scf_config%talk) then
            write (line, "(a,i0,a,f20.12,a,es9.2,a,es9.2,a,f7.4,a,es9.2,a,i0)") &
               "   soscf ", iteration, ": E = ", e_cur, "  dE ", de, "  |FDS-SDF| ", errmax, &
               "  trust ", trust, "  curvature ", io%lowest, "  fock builds ", builds
            call global_logger%info(trim(line))
         end if
         scf_result%iterations = switch_iteration + iteration

         if (.not. accepted) then
            call global_logger%warning("   the second-order SCF found no step downhill; stopping. "// &
                                       "the reported energy is the last one that descended.")
            exit
         end if
      end do

      ! `scf_result%dmat`/`scf_energies` may hold a REJECTED trial
      ! (backtracking exhausted, or the loop above hit `max_steps`):
      ! restore the last ACCEPTED density and rebuild the Fock/energy
      ! from it, so what leaves is the point that actually descended,
      ! never a half-written trial's -- `scf_energies%total` (and
      ! `%gmat`, read by the energy) was otherwise left by whichever
      ! trial's `scf_soscf_fock_rebuild` ran last, accepted or not.
      call cphf_copy_matrix(dims%nao, dims%nao, dmat_accepted, scf_result%dmat(:, :, 1))
      call scf_soscf_fock_rebuild(scf_basis, scf_context, scf_config, scf_state, scf_result, scf_energies, la, eri)
      call scf_soscf_semicanonicalize(la, dims, fock_cur, coeff_work, fock_mo, energies_mo, info)
      if (info == 0) then
         call cphf_copy_matrix(dims%nao, dims%n_mo, coeff_work, scf_result%cmo(:, 1:dims%n_mo, 1))
         call soscf_copy_eps(dims%nao, dims%n_mo, energies_mo, scf_result%eps(:, 1))
      else
         ! The orbitals returned are still the last ACCEPTED ones (coeff_work
         ! is only ever written on acceptance): not semicanonical, but not
         ! stale either. eps is left at whatever the main loop's last
         ! successful semicanonicalisation wrote -- the best available,
         ! though no longer exactly consistent with this coeff.
         call cphf_copy_matrix(dims%nao, dims%n_mo, coeff_work, scf_result%cmo(:, 1:dims%n_mo, 1))
      end if
      ! `scf_state%dold` feeds `scf_result%update` and the analytic
      ! gradient's own energy-weighted density -- it has to be the SAME
      ! density as `scf_result%dmat` above, not whatever `trc_scf_run`'s
      ! own DIIS loop last wrote there before handing over.
      call cphf_copy_matrix(dims%nao, dims%nao, scf_result%dmat(:, :, 1), scf_state%dold(:, :, 1))
      if (scf_config%talk) then
         write (line, "(a,i0,a,i0,a)") "   second-order phase: ", products_total, &
            " Hessian-vector products (Fock builds) over ", max(scf_result%iterations - switch_iteration, 0), &
            " Newton iterations"
         call global_logger%info(trim(line))
      end if

      !$acc exit data delete(fock_cur, coeff_work, coeff_try, dmat_accepted, fock_mo) &
      !$acc           delete(io%gradient, io%step, step_scaled, kappa, rotation, scratch_a, scratch_b, scratch_c)
      !$acc exit data delete(energies_mo)
   end subroutine trc_scf_soscf_run

   !
   ! One full (non-incremental) Fock build and the energy at the density
   ! `scf_result%dmat(:,:,1)`, restricted Hartree-Fock only (no XC, no
   ! continuum -- `trc_scf_settings%resolve` has already refused both for
   ! `--soscf`): `scf_fock_build_g`/`scf_fock_energy`/`scf_fock_assemble`,
   ! the SAME sequence `trc_scf_run`'s own iteration runs, reused rather
   ! than duplicated. Leaves `scf_state%fock(:,:,1)` and
   ! `scf_energies%total` at that density.
   !
   subroutine scf_soscf_fock_rebuild(scf_basis, scf_context, scf_config, scf_state, scf_result, scf_energies, la, eri)
      type(trc_basis_t), intent(in) :: scf_basis
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(scf_energies_t), intent(inout) :: scf_energies
      type(trc_linalg_t), intent(inout) :: la
      type(trc_eri_t), intent(inout) :: eri
      integer :: s
      do s = 1, scf_config%nspin
         call scf_fock_to_cart(scf_basis, scf_context%sph, scf_context%nao, scf_context%nao_c, &
                               scf_result%dmat(:, :, s), scf_state%d_c(:, :, s))
      end do
      call scf_fock_build_g(scf_basis, scf_context, scf_config, scf_state, .false., eri, la)
      do s = 1, scf_config%nspin
         call scf_fock_to_sph(scf_basis, scf_context%sph, scf_context%nao, scf_context%nao_c, &
                              scf_state%gmat_c(:, :, s), scf_state%gmat(:, :, s))
      end do
      call scf_fock_energy(scf_context, scf_config, scf_state, scf_result, scf_energies, la)
      call scf_fock_assemble(scf_context, scf_config, scf_state)
   end subroutine scf_soscf_fock_rebuild

   !
   ! Diagonalise the Fock matrix within the occupied block and within the
   ! virtual block, so the Hessian expression `scf_soscf_hessian_apply_rhf`
   ! computes means the Hessian again -- see this module's own header.
   ! `coeff` (dims%nao, dims%n_mo) is rotated in place, within each block;
   ! what comes out is the Fock matrix in the new basis (`fock_mo`), whose
   ! occupied-virtual block is the gradient and whose diagonal
   ! (`energies_mo`) is the orbital energies -- canonical once the
   ! occupied-virtual block has vanished, semicanonical otherwise.
   !
   subroutine scf_soscf_semicanonicalize(la, dims, fock_ao, coeff, fock_mo, energies_mo, info)
      type(trc_linalg_t), intent(inout) :: la
      type(trc_soscf_dims_t), intent(in) :: dims
      real(dp), intent(in) :: fock_ao(dims%nao, dims%nao)
      real(dp), intent(inout) :: coeff(dims%nao, dims%n_mo)
      real(dp), intent(inout) :: fock_mo(dims%n_mo, dims%n_mo)
      real(dp), intent(inout) :: energies_mo(dims%n_mo)
      integer, intent(out) :: info
      integer :: nao, n_mo, nocc, n_vir
      real(dp), allocatable :: work(:, :), rotated(:, :), block_matrix(:, :), block_values(:)

      nao = dims%nao
      n_mo = dims%n_mo
      nocc = dims%nocc
      n_vir = n_mo - nocc
      allocate (work(nao, n_mo), rotated(nao, n_mo))
      !$acc enter data create(work, rotated)

      call la%gemm_mat(fock_ao, coeff, work)
      call la%gemm_mat(coeff, work, fock_mo, transa="T")

      allocate (block_matrix(nocc, nocc), block_values(nocc))
      !$acc enter data create(block_matrix, block_values)
      call soscf_copy_block(n_mo, n_mo, fock_mo, 1, 1, nocc, nocc, block_matrix)
      call la%syev(block_matrix, block_values, info)
      if (info /= 0) then
         !$acc exit data delete(block_matrix, block_values)
         deallocate (block_matrix, block_values)
         !$acc exit data delete(work, rotated)
         deallocate (work, rotated)
         return
      end if
      call la%gemm_mat(coeff(:, 1:nocc), block_matrix, rotated(:, 1:nocc))
      !$acc exit data delete(block_matrix, block_values)
      deallocate (block_matrix, block_values)

      allocate (block_matrix(n_vir, n_vir), block_values(n_vir))
      !$acc enter data create(block_matrix, block_values)
      call soscf_copy_block(n_mo, n_mo, fock_mo, nocc + 1, nocc + 1, n_vir, n_vir, block_matrix)
      call la%syev(block_matrix, block_values, info)
      if (info /= 0) then
         !$acc exit data delete(block_matrix, block_values)
         deallocate (block_matrix, block_values)
         !$acc exit data delete(work, rotated)
         deallocate (work, rotated)
         return
      end if
      call la%gemm_mat(coeff(:, nocc + 1:n_mo), block_matrix, rotated(:, nocc + 1:n_mo))
      !$acc exit data delete(block_matrix, block_values)
      deallocate (block_matrix, block_values)

      call cphf_copy_matrix(nao, n_mo, rotated, coeff)
      ! Rebuilt in the new basis rather than assembled from the two block
      ! rotations: the occupied-virtual block is the gradient and has to
      ! be the transformed one.
      call la%gemm_mat(fock_ao, coeff, work)
      call la%gemm_mat(coeff, work, fock_mo, transa="T")
      call soscf_diagonal(n_mo, fock_mo, energies_mo)

      !$acc exit data delete(work, rotated)
      deallocate (work, rotated)
   end subroutine scf_soscf_semicanonicalize

   subroutine soscf_copy_block(m, n, a, i0, j0, bm, bn, block)
      !! block(1:bm,1:bn) := a(i0:i0+bm-1, j0:j0+bn-1): a(m,n) is read by
      !! offset rather than handed a (non-contiguous) section, so this is a
      !! legal device kernel on the WHOLE resident array.
      integer, intent(in) :: m, n, i0, j0, bm, bn
      real(dp), intent(in) :: a(m, n)
      real(dp), intent(inout) :: block(bm, bn)
      integer :: i, j
      do concurrent(j=1:bn, i=1:bm)
         block(i, j) = a(i0 + i - 1, j0 + j - 1)
      end do
   end subroutine soscf_copy_block

   subroutine soscf_diagonal(n, a, d)
      integer, intent(in) :: n
      real(dp), intent(in) :: a(n, n)
      real(dp), intent(inout) :: d(n)
      integer :: k
      do concurrent(k=1:n)
         d(k) = a(k, k)
      end do
   end subroutine soscf_diagonal

   subroutine soscf_copy_eps(nao, n_mo, energies_mo, eps)
      !! eps(1:n_mo) := energies_mo, device-resident: the final commit of
      !! the semicanonical orbital energies into `scf_result%eps(:,1)`.
      !! Done as a kernel on the device-resident arrays, not a host
      !! assignment after `!$acc update self`: the latter leaves the
      !! DEVICE copy of `eps` stale, so `trc_scf_run`'s own trailing
      !! `!$acc update self(..., scf_result%eps)` pulls the OLD value
      !! straight back over it (a GPU-only staleness bug, invisible on
      !! the host build, where there is no separate device copy to go
      !! stale). `eps(n_mo+1:nao)` -- the null states canonical
      !! orthogonalisation dropped -- is left untouched: never written
      !! here, never read as real by anything downstream.
      integer, intent(in) :: nao, n_mo
      real(dp), intent(in) :: energies_mo(n_mo)
      real(dp), intent(inout) :: eps(nao)
      integer :: k
      do concurrent(k=1:n_mo)
         eps(k) = energies_mo(k)
      end do
   end subroutine soscf_copy_eps

   subroutine scf_soscf_gradient(n_mo, nocc, nvirt, fock_mo, gradient)
      !! `dE/d kappa_ai = 4 F_ai`, virtual fastest -- see this module's own
      !! header, point 2.
      integer, intent(in) :: n_mo, nocc, nvirt
      real(dp), intent(in) :: fock_mo(n_mo, n_mo)
      real(dp), intent(inout) :: gradient(nvirt, nocc)
      integer :: i, a
      do concurrent(i=1:nocc, a=1:nvirt)
         gradient(a, i) = HESSIAN_SCALE*fock_mo(nocc + a, i)
      end do
   end subroutine scf_soscf_gradient

   !
   ! The trust-region Newton step, from Hessian-vector products alone:
   ! builds a Krylov subspace off the preconditioned gradient, projects the
   ! Hessian and the gradient into it, and solves the small, dense,
   ! level-shifted subproblem there (`soscf_level_shifted_step`) -- the
   ! same step control CASSCF's own optimiser would apply to a dense
   ! Hessian, applied here to a projected one. `energies_mo` is HOST
   ! (pulled once by the caller); `cmo` and `io%gradient`/`io%step` are
   ! device resident. `io%lowest` is an UPPER BOUND on the true smallest
   ! curvature (a Rayleigh-Ritz value), not the curvature itself. The
   ! Hessian-vector product's own scratch (`trc_soscf_hvp_scratch_t`) is
   ! local to this call, never a dummy here -- see this module's own
   ! header on why the bundling stops at the arg-count limit, not at
   ! "everything this touches".
   !
   subroutine scf_soscf_newton_step_rhf(b, eri, lp, dims, cmo, energies_mo, io)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_soscf_dims_t), intent(in) :: dims
      real(dp), intent(in) :: cmo(dims%nao, dims%n_mo)
      real(dp), intent(in) :: energies_mo(dims%n_mo)
      type(trc_soscf_step_io_t), intent(inout) :: io

      type(trc_soscf_hvp_scratch_t) :: scratch
      real(dp), allocatable :: basis(:, :, :), image(:, :, :)
      real(dp), allocatable :: direction(:, :), residual(:, :)
      real(dp), allocatable :: small(:, :), small_gradient(:), amplitude(:)
      real(dp) :: gradient_norm, resid_norm
      integer :: nao, nao_c, n_mo, nocc, nvirt, nmax, nsub, i, j, a0, i0
      integer :: softest(2)
      logical :: added, ok

      nao = dims%nao
      nao_c = dims%nao_c
      n_mo = dims%n_mo
      nocc = dims%nocc
      nvirt = dims%nvirt
      io%lowest = 0.0_dp
      io%predicted = 0.0_dp
      io%products = 0
      io%ok = .true.
      call cphf_zero_matrix(nvirt, nocc, io%step)
      if (nvirt < 1 .or. nocc < 1) return
      if (soscf_force_breakdown_env()) then
         ! TRC_SOSCF_FORCE_BREAKDOWN=1: an injected failure for
         ! check_soscf.F90's own test that `trc_scf_soscf_run` stops
         ! immediately on a Krylov/level-shift breakdown rather than
         ! looping to max_iteration -- nothing in ordinary use sets this.
         io%ok = .false.
         return
      end if

      allocate (scratch%diagv(nvirt, nocc))
      do i = 1, nocc
         do j = 1, nvirt
            scratch%diagv(j, i) = energies_mo(nocc + j) - energies_mo(i)
         end do
      end do
      !$acc enter data copyin(scratch%diagv)

      nmax = max(1, min(io%max_subspace, nvirt*nocc))
      allocate (basis(nvirt, nocc, nmax), image(nvirt, nocc, nmax))
      allocate (direction(nvirt, nocc), residual(nvirt, nocc))
      allocate (scratch%d_u(nao, nao), scratch%d_ut(nao, nao), scratch%g_u(nao, nao))
      allocate (scratch%d_u_c(nao_c, nao_c), scratch%g_u_c(nao_c, nao_c))
      allocate (scratch%temporary_no(nao, nocc), scratch%temporary_no2(nao, nocc))
      !$acc enter data create(basis, image, direction, residual)
      !$acc enter data create(scratch%d_u, scratch%d_ut, scratch%g_u, scratch%d_u_c, scratch%g_u_c)
      !$acc enter data create(scratch%temporary_no, scratch%temporary_no2)

      gradient_norm = sqrt(lp%dot_2d(io%gradient, io%gradient))

      ! ---- the seeds: the preconditioned gradient, and the softest mode --
      nsub = 0
      call soscf_preconditioned_direction(nvirt, nocc, scratch%diagv, PRECONDITION_FLOOR, io%gradient, direction)
      call soscf_add_direction(lp, nvirt, nocc, basis, nsub, direction, added)
      softest = minloc(scratch%diagv)
      a0 = softest(1)
      i0 = softest(2)
      call soscf_unit_matrix(nvirt, nocc, a0, i0, direction)
      call soscf_add_direction(lp, nvirt, nocc, basis, nsub, direction, added)

      if (nsub > 0) then
         do i = 1, nsub
            call scf_soscf_hessian_apply_rhf(b, eri, lp, dims, cmo, basis(:, :, i), image(:, :, i), scratch)
            call soscf_scale_matrix(nvirt, nocc, image(:, :, i), HESSIAN_SCALE, image(:, :, i))
            io%products = io%products + 1
         end do

         ! ---- solve, expand, repeat --------------------------------------
         do
            if (allocated(small)) deallocate (small, small_gradient)
            allocate (small(nsub, nsub), small_gradient(nsub))
            do j = 1, nsub
               do i = 1, nsub
                  small(i, j) = lp%dot_2d(basis(:, :, i), image(:, :, j))
               end do
               small_gradient(j) = lp%dot_2d(basis(:, :, j), io%gradient)
            end do
            small = 0.5_dp*(small + transpose(small))

            if (allocated(amplitude)) deallocate (amplitude)
            allocate (amplitude(nsub))
            call soscf_level_shifted_step(small, small_gradient, io%escape, amplitude, io%lowest, io%predicted, ok)
            if (.not. ok) then
               call global_logger%warning("trc_scf_soscf: the projected orbital Hessian failed to diagonalise")
               io%ok = .false.
               exit
            end if

            call cphf_zero_matrix(nvirt, nocc, io%step)
            call cphf_copy_matrix(nvirt, nocc, io%gradient, residual)
            do i = 1, nsub
               call cphf_axpy_matrix(nvirt, nocc, amplitude(i), basis(:, :, i), io%step)
               call cphf_axpy_matrix(nvirt, nocc, amplitude(i), image(:, :, i), residual)
            end do

            if (nsub >= nmax) exit
            resid_norm = sqrt(lp%dot_2d(residual, residual))
            if (resid_norm <= SUBSPACE_TOLERANCE*max(gradient_norm, tiny(1.0_dp))) exit

            call soscf_preconditioned_direction(nvirt, nocc, scratch%diagv, PRECONDITION_FLOOR, residual, direction)
            call soscf_add_direction(lp, nvirt, nocc, basis, nsub, direction, added)
            if (.not. added) exit

            call scf_soscf_hessian_apply_rhf(b, eri, lp, dims, cmo, basis(:, :, nsub), image(:, :, nsub), scratch)
            call soscf_scale_matrix(nvirt, nocc, image(:, :, nsub), HESSIAN_SCALE, image(:, :, nsub))
            io%products = io%products + 1
         end do
      else
         ! Unreachable for nvirt, nocc >= 1 (the second seed is always a
         ! unit vector, and the first call to soscf_add_direction has
         ! nothing yet to project it against) -- defensive only.
         call global_logger%warning("trc_scf_soscf: no direction to build a Krylov subspace from")
         io%ok = .false.
      end if

      if (allocated(small)) deallocate (small, small_gradient)
      if (allocated(amplitude)) deallocate (amplitude)
      !$acc exit data delete(basis, image, direction, residual)
      deallocate (basis, image, direction, residual)
      !$acc exit data delete(scratch%d_u, scratch%d_ut, scratch%g_u, scratch%d_u_c, scratch%g_u_c)
      !$acc exit data delete(scratch%temporary_no, scratch%temporary_no2)
      deallocate (scratch%d_u, scratch%d_ut, scratch%g_u, scratch%d_u_c, scratch%g_u_c)
      deallocate (scratch%temporary_no, scratch%temporary_no2)
      !$acc exit data delete(scratch%diagv)
      deallocate (scratch%diagv)
   end subroutine scf_soscf_newton_step_rhf

   !
   ! `(e_a-e_i) p + [A p]`, the RHF orbital-Hessian action on a trial
   ! rotation `pin`, folded through the Cartesian space `eri%fock_resident`
   ! runs in -- see this module's own header. `scratch%d_u`/`%d_ut`/`%g_u`
   ! are `nao`-by-`nao` (the SCF's own space, spherical when the basis
   ! is); `%d_u_c`/`%g_u_c` are `nao_c`-by-`nao_c` (always Cartesian);
   ! `dims%sph = .false.` makes `scf_fock_to_cart`/`scf_fock_to_sph` a
   ! plain copy, so a Cartesian basis (`nao == nao_c`) pays nothing beyond
   ! it.
   !
   subroutine scf_soscf_hessian_apply_rhf(b, eri, lp, dims, cmo, pin, apout, scratch)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_soscf_dims_t), intent(in) :: dims
      real(dp), intent(in) :: cmo(dims%nao, dims%n_mo)
      real(dp), intent(in) :: pin(dims%nvirt, dims%nocc)
      real(dp), intent(inout) :: apout(dims%nvirt, dims%nocc)
      type(trc_soscf_hvp_scratch_t), intent(inout) :: scratch
      integer :: nocc, n_mo

      nocc = dims%nocc
      n_mo = dims%n_mo
      call lp%gemm_mat(cmo(:, nocc + 1:n_mo), pin, scratch%temporary_no)
      call lp%gemm_mat(scratch%temporary_no, cmo(:, 1:nocc), scratch%d_u, transb="T")
      call lp%transpose(dims%nao, dims%nao, scratch%d_u, dims%nao, scratch%d_ut, dims%nao)
      call cphf_add_matrix(dims%nao, dims%nao, scratch%d_u, scratch%d_ut, scratch%d_u)
      call scf_fock_to_cart(b, dims%sph, dims%nao, dims%nao_c, scratch%d_u, scratch%d_u_c)
      call eri%fock_resident(b, scratch%d_u_c, scratch%g_u_c, k_scale=1.0_dp, j_scale=1.0_dp, density_screen=.false.)
      call scf_fock_to_sph(b, dims%sph, dims%nao, dims%nao_c, scratch%g_u_c, scratch%g_u)
      call lp%gemm_mat(scratch%g_u, cmo(:, 1:nocc), scratch%temporary_no2)
      call lp%gemm_mat(cmo(:, nocc + 1:n_mo), scratch%temporary_no2, apout, transa="T", alpha=2.0_dp)
      call cphf_add_diag(dims%nvirt, nocc, scratch%diagv, pin, apout)
   end subroutine scf_soscf_hessian_apply_rhf

   subroutine soscf_add_direction(lp, m, n, basis, nsub, direction, added)
      !! Normalise `direction`, orthogonalise it against the `nsub`
      !! existing columns of `basis` (two Gram-Schmidt passes), and append
      !! it as column `nsub+1` if what survives is not negligible.
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: m, n
      integer, intent(inout) :: nsub
      real(dp), intent(inout) :: basis(m, n, *)
      real(dp), intent(inout) :: direction(m, n)
      logical, intent(out) :: added
      real(dp) :: norm, overlap
      integer :: j, pass
      added = .false.
      norm = sqrt(lp%dot_2d(direction, direction))
      if (.not. (norm > tiny(1.0_dp))) return
      call soscf_scale_matrix(m, n, direction, 1.0_dp/norm, direction)
      do pass = 1, 2
         do j = 1, nsub
            overlap = lp%dot_2d(basis(:, :, j), direction)
            call cphf_axpy_matrix(m, n, -overlap, basis(:, :, j), direction)
         end do
      end do
      norm = sqrt(lp%dot_2d(direction, direction))
      if (norm < LINEAR_DEPENDENCE) return
      nsub = nsub + 1
      call soscf_scale_matrix(m, n, direction, 1.0_dp/norm, basis(:, :, nsub))
      added = .true.
   end subroutine soscf_add_direction

   subroutine soscf_preconditioned_direction(m, n, gaps, floor_value, g, direction)
      !! direction := -g / max(gaps, floor_value), elementwise.
      integer, intent(in) :: m, n
      real(dp), intent(in) :: floor_value
      real(dp), intent(in) :: gaps(m, n)
      real(dp), intent(in) :: g(m, n)
      real(dp), intent(inout) :: direction(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         direction(i, j) = -g(i, j)/max(gaps(i, j), floor_value)
      end do
   end subroutine soscf_preconditioned_direction

   subroutine soscf_scale_matrix(m, n, a, alpha, c)
      !! c := alpha*a, elementwise (a and c may be the same actual array).
      integer, intent(in) :: m, n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: a(m, n)
      real(dp), intent(inout) :: c(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         c(i, j) = alpha*a(i, j)
      end do
   end subroutine soscf_scale_matrix

   subroutine soscf_unit_matrix(m, n, a0, i0, a)
      !! a := the (a0, i0) unit matrix.
      integer, intent(in) :: m, n, a0, i0
      real(dp), intent(inout) :: a(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         a(i, j) = merge(1.0_dp, 0.0_dp, i == a0 .and. j == i0)
      end do
   end subroutine soscf_unit_matrix

   subroutine soscf_identity_matrix(n, a)
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(n, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:n)
         a(i, j) = merge(1.0_dp, 0.0_dp, i == j)
      end do
   end subroutine soscf_identity_matrix

   subroutine soscf_embed_kappa(n_mo, nocc, nvirt, step_in, kappa)
      !! kappa := the antisymmetric n_mo x n_mo matrix whose occupied-
      !! virtual block is `step_in` -- zero elsewhere: rotating occupied
      !! orbitals among themselves, or virtuals among themselves, is
      !! redundant for a closed-shell determinant (this module's own
      !! header), and `trc_scf_soscf_run`'s own n_mo = nao - n_lindep
      !! keeps the dropped (null) directions out of this matrix entirely.
      integer, intent(in) :: n_mo, nocc, nvirt
      real(dp), intent(in) :: step_in(nvirt, nocc)
      real(dp), intent(inout) :: kappa(n_mo, n_mo)
      integer :: i, j, a
      do concurrent(j=1:n_mo, i=1:n_mo)
         kappa(i, j) = 0.0_dp
      end do
      do concurrent(i=1:nocc, a=1:nvirt)
         kappa(nocc + a, i) = step_in(a, i)
         kappa(i, nocc + a) = -step_in(a, i)
      end do
   end subroutine soscf_embed_kappa

   subroutine soscf_pull_matrix(m, n, a, a_host)
      !! a_host := a, device to host -- the one place this module reads a
      !! non-scalar device array back, for `maxval(abs(.))` (a host-only
      !! intrinsic): trust-radius scaling, once per Newton iteration, on
      !! the (small) nvirt x nocc step, never on the full nmo x nmo
      !! rotation.
      integer, intent(in) :: m, n
      real(dp), intent(in) :: a(m, n)
      real(dp), intent(inout) :: a_host(m, n)
      !$acc update self(a)
      a_host = a
   end subroutine soscf_pull_matrix

   pure function soscf_ov_one_norm(nvirt, nocc, step) result(norm)
      !! A 1-norm bound (the largest absolute COLUMN sum -- a true matrix
      !! norm, unlike the largest single element) on the antisymmetric
      !! `n_mo`-by-`n_mo` matrix `soscf_embed_kappa` builds from a step of
      !! this shape, computed from the (small) occupied-virtual block
      !! alone rather than the full matrix: column `i` (occupied) is
      !! `step(:,i)`, so its abs-sum is the 1-norm of that column; column
      !! `nocc+a` (virtual) is `-step(a,:)`, so its abs-sum is the 1-norm
      !! of that ROW of `step`. `scf_soscf_rotation_matrix`'s own
      !! scaling-and-squaring needs a bound like this on the matrix it is
      !! actually exponentiating to converge; the largest single element
      !! bounds nothing (a matrix with one huge entry and the rest zero
      !! has a tiny operator norm; one with every entry equal to that same
      !! "largest element" does not), which is how `exp(kappa)` lost
      !! orthogonality at larger nvirt/nocc.
      integer, intent(in) :: nvirt, nocc
      real(dp), intent(in) :: step(nvirt, nocc)
      real(dp) :: norm
      integer :: i, a
      norm = 0.0_dp
      do i = 1, nocc
         norm = max(norm, sum(abs(step(:, i))))
      end do
      do a = 1, nvirt
         norm = max(norm, sum(abs(step(a, :))))
      end do
   end function soscf_ov_one_norm

   real(dp) function soscf_trust0_env() result(res)
      !! TRC_SOSCF_TRUST0: override the trust radius a --soscf phase
      !! starts from -- MAX_ROTATION unless set. An experiment hook, like
      !! trc_scf_settings.F90's own TRC_INCR_START and friends: nothing in
      !! ordinary use sets it. check_soscf.F90's own backtracking test
      !! sets it deliberately larger than a real Newton step's own size,
      !! to force at least one rejection.
      character(len=16) :: v
      integer :: ios
      real(dp) :: t
      call get_environment_variable("TRC_SOSCF_TRUST0", v)
      res = MAX_ROTATION
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) t
         if (ios == 0 .and. t > 0.0_dp) res = t
      end if
   end function soscf_trust0_env

   pure function soscf_trust_after_rejection(trust, scaling) result(new_trust)
      !! The trust radius a rejection leaves behind, given the step's own
      !! (fixed, unscaled) size `scaling`: halving `trust` alone is a
      !! no-op on the NEXT trial's actual step for as long as `trust`
      !! stays above `scaling` (the cap is not binding, so the step tried
      !! is the same `scaling`-sized one every time) -- the bug this
      !! module's own history had, where a backtracking loop could retry
      !! one identical step for every one of MAX_BACKTRACKS trials.
      !! Snapping `trust` to at most `scaling` first makes the cap bind
      !! from the very next trial on, so every rejection actually
      !! produces a smaller step. A pure function, exported, so
      !! check_soscf.F90 can test this arithmetic directly rather than
      !! only through a full SCF run that may or may not happen to reject
      !! a step on the system it has in hand.
      real(dp), intent(in) :: trust, scaling
      real(dp) :: new_trust
      new_trust = 0.5_dp*min(trust, scaling)
   end function soscf_trust_after_rejection

   logical function soscf_force_breakdown_env() result(res)
      !! TRC_SOSCF_FORCE_BREAKDOWN=1: an injected Krylov/level-shift
      !! breakdown, for check_soscf.F90's own test of the stop-immediately
      !! path -- see soscf_trust0_env's own header on this module's
      !! experiment-hook convention. Nothing in ordinary use sets it.
      character(len=8) :: v
      call get_environment_variable("TRC_SOSCF_FORCE_BREAKDOWN", v)
      res = trim(v) == "1"
   end function soscf_force_breakdown_env

   !
   ! `exp(kappa)` for antisymmetric `kappa`, by scaling and squaring:
   ! orthogonal to machine precision for any step size (ported from
   ! metalquicha's `mqc_orbital_rotation%rotation_matrix`, same algorithm,
   ! device-resident GEMMs). `norm` is a 1-norm bound on `kappa`
   ! (`soscf_ov_one_norm`), passed in rather than recomputed here: the
   ! caller already knows it (from the small occupied-virtual step, not
   ! the full n x n matrix), and pulling the full n x n `kappa` back to
   ! the host just to recompute a bound on it would be the one
   ! unnecessary device-to-host transfer this module otherwise avoids.
   ! `scaled`/`term`/`next_term` are the caller's own n x n scratch.
   !
   subroutine scf_soscf_rotation_matrix(lp, n, norm, kappa, rotation, scaled, term, next_term)
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: n
      real(dp), intent(in) :: norm
      real(dp), intent(in) :: kappa(n, n)
      real(dp), intent(inout) :: rotation(n, n)
      real(dp), intent(inout) :: scaled(n, n), term(n, n), next_term(n, n)
      real(dp) :: nrm, term_norm
      integer :: squarings, k
      integer, parameter :: TAYLOR_TERMS = 18

      nrm = norm
      squarings = 0
      do while (nrm > 0.5_dp)
         nrm = 0.5_dp*nrm
         squarings = squarings + 1
      end do
      call soscf_scale_matrix(n, n, kappa, 1.0_dp/real(2**squarings, dp), scaled)
      call soscf_identity_matrix(n, rotation)
      call soscf_identity_matrix(n, term)
      do k = 1, TAYLOR_TERMS
         call lp%gemm_mat(term, scaled, next_term)
         call soscf_scale_matrix(n, n, next_term, 1.0_dp/real(k, dp), term)
         call cphf_add_matrix(n, n, rotation, term, rotation)
         term_norm = sqrt(lp%dot_2d(term, term))
         if (term_norm < 1.0e-16_dp) exit
      end do
      do k = 1, squarings
         call lp%gemm_mat(rotation, rotation, next_term)
         call cphf_copy_matrix(n, n, next_term, rotation)
      end do
   end subroutine scf_soscf_rotation_matrix

   !
   ! The level-shifted trust-region Newton step on a SMALL (at most
   ! `DEFAULT_SOSCF_SUBSPACE` square), host-resident projected Hessian --
   ! ported from metalquicha's `mqc_orbital_rotation%level_shifted_step`,
   ! same algorithm, host LAPACK (`pic_syev`/`pic_gemv`) rather than
   ! device: the same reasoning `trc_scf_diis`'s own small linear solves
   ! use (a handful of numbers, not worth a device round trip).
   !
   subroutine soscf_level_shifted_step(hessian, gradient, escape, step, lowest, predicted, ok)
      real(dp), intent(in) :: hessian(:, :)
      real(dp), intent(in) :: gradient(:)
      real(dp), intent(in) :: escape
      real(dp), intent(out) :: step(size(gradient))
      real(dp), intent(out) :: lowest, predicted
      logical, intent(out) :: ok
      real(dp), allocatable :: vectors(:, :), values(:), projected(:), amplitude(:)
      real(dp) :: shift
      integer :: n_param, k
      integer(default_int) :: info

      n_param = size(gradient)
      lowest = 0.0_dp
      predicted = 0.0_dp
      step = 0.0_dp
      ok = .true.
      if (n_param == 0) return
      allocate (vectors(n_param, n_param), values(n_param), projected(n_param), amplitude(n_param))
      vectors = hessian
      call pic_syev(vectors, values, "V", "U", info)
      if (info /= 0) then
         ok = .false.
         deallocate (vectors, values, projected, amplitude)
         return
      end if
      lowest = values(1)

      shift = 0.0_dp
      if (values(1) < MIN_CURVATURE) shift = MIN_CURVATURE - values(1)
      call pic_gemv(vectors, gradient, projected, trans_a="T")

      do k = 1, n_param
         amplitude(k) = -projected(k)/(values(k) + shift)
         if (values(k) < SADDLE_CURVATURE .and. abs(amplitude(k)) < escape) then
            amplitude(k) = sign(escape, amplitude(k))
         end if
      end do
      do k = 1, n_param
         predicted = predicted - projected(k)*amplitude(k) - 0.5_dp*values(k)*amplitude(k)**2
      end do
      call pic_gemv(vectors, amplitude, step)
      deallocate (vectors, values, projected, amplitude)
   end subroutine soscf_level_shifted_step

end module trc_scf_soscf
