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
!    UNITS the curvature thresholds below are read in, which is why it is
!    re-derived independently, by finite difference, in this port's own
!    test (`check_soscf.F90`) rather than carried over on trust alone.
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
   use trc_scf_settings, only: scf_settings_t, trc_scf_options_t, trc_scf_result_t, scf_energies_t
   use trc_scf_state, only: scf_state_t
   use trc_scf_fock, only: scf_fock_build_g, scf_fock_energy, scf_fock_assemble, scf_fock_commutator_error, &
                           scf_fock_to_cart, scf_fock_to_sph
   implicit none
   private

   public :: trc_scf_soscf_run
   ! Exported for check_soscf.F90's own finite-difference validation of the
   ! gradient and the Hessian scale -- see this module's own header, point 2.
   public :: scf_soscf_semicanonicalize, scf_soscf_gradient, scf_soscf_newton_step_rhf, scf_soscf_rotation_matrix
   public :: scf_soscf_fock_rebuild, scf_soscf_hessian_apply_rhf
   public :: HESSIAN_SCALE, MIN_CURVATURE, SADDLE_CURVATURE, MAX_ROTATION, MIN_ROTATION, ENERGY_RESOLUTION, &
             TRUST_GROWTH, DEFAULT_SOSCF_SUBSPACE, SUBSPACE_TOLERANCE, MAX_BACKTRACKS

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

contains

   !
   ! The second-order phase: `trc_scf_run`'s own DIIS loop calls this once
   ! the commutator has fallen below `scf_options%soscf_start`, in place of
   ! the rest of its own iterations. `switch_iteration` is the DIIS iteration
   ! this took over after (reported, not used); `max_steps` is whatever is
   ! left of `scf_options%max_iteration`. Updates `scf_result` (density,
   ! orbitals, orbital energies, `converged`, `iterations`) in place;
   ! `scf_state`/`scf_energies` are the same per-iteration scratch and
   ! energy bundle the DIIS loop itself writes through, reused here rather
   ! than duplicated.
   !
   subroutine trc_scf_soscf_run(scf_basis, scf_context, scf_config, scf_state, scf_options, scf_result, &
                                scf_energies, la, eri, switch_iteration, max_steps)
      type(trc_basis_t), intent(in) :: scf_basis
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      type(trc_scf_options_t), intent(in) :: scf_options
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(scf_energies_t), intent(inout) :: scf_energies
      type(trc_linalg_t), intent(inout) :: la
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: switch_iteration
      integer, intent(in) :: max_steps

      integer :: nao, n_mo, nocc, nvirt, iteration, trial, info, products, products_total, builds
      real(dp) :: e_cur, e_try, e_old, de, errmax, scaling, trust, predicted, lowest
      logical :: accepted, saddle, converged_now
      character(len=200) :: line

      real(dp), allocatable :: fock_cur(:, :), coeff_work(:, :), coeff_try(:, :), dmat_accepted(:, :)
      real(dp), allocatable :: fock_mo(:, :), energies_mo(:)
      real(dp), allocatable :: gradient(:, :), step(:, :), step_scaled(:, :)
      real(dp), allocatable :: kappa(:, :), rotation(:, :), scratch_a(:, :), scratch_b(:, :), scratch_c(:, :)
      real(dp), allocatable :: step_host(:, :)

      nao = scf_context%nao
      n_mo = nao - scf_context%n_lindep
      nocc = scf_config%nocc(1)
      nvirt = n_mo - nocc
      products_total = 0
      if (nocc < 1 .or. nvirt < 1) then
         scf_result%message = "trc_scf: --soscf needs at least one occupied and one virtual orbital; "// &
                              "there are no rotations to take"
         return
      end if
      if (max_steps < 1) return

      allocate (fock_cur(nao, nao), coeff_work(nao, n_mo), coeff_try(nao, n_mo), dmat_accepted(nao, nao))
      allocate (fock_mo(n_mo, n_mo), energies_mo(n_mo))
      allocate (gradient(nvirt, nocc), step(nvirt, nocc), step_scaled(nvirt, nocc))
      allocate (kappa(n_mo, n_mo), rotation(n_mo, n_mo), scratch_a(n_mo, n_mo), scratch_b(n_mo, n_mo), &
               scratch_c(n_mo, n_mo))
      allocate (step_host(nvirt, nocc))
      !$acc enter data create(fock_cur, coeff_work, coeff_try, dmat_accepted, fock_mo, gradient, step) &
      !$acc            create(step_scaled, kappa, rotation, scratch_a, scratch_b, scratch_c)
      !$acc enter data copyin(energies_mo)

      call cphf_copy_matrix(nao, n_mo, scf_result%cmo(:, 1:n_mo, 1), coeff_work)
      call cphf_copy_matrix(nao, nao, scf_result%dmat(:, :, 1), dmat_accepted)

      ! The point the Newton phase starts from, built in full: the Fock
      ! matrix DIIS left behind is extrapolated and possibly level
      ! shifted, and its energy belongs to the density before the last
      ! diagonalisation -- neither is a state this can differentiate
      ! (mqc_czt_rhf's own run_soscf_phase, same reasoning).
      call scf_soscf_fock_rebuild(scf_basis, scf_context, scf_config, scf_state, scf_result, scf_energies, la, eri)
      call cphf_copy_matrix(nao, nao, scf_state%fock(:, :, 1), fock_cur)
      e_cur = scf_energies%total
      e_old = e_cur
      trust = MAX_ROTATION
      if (scf_config%talk) then
         call global_logger%info("   second-order phase: trust-region Newton on the orbital rotations, "// &
                                 "handing over from DIIS at iteration "//to_string(switch_iteration))
      end if

      do iteration = 1, max_steps
         builds = 0
         call scf_soscf_semicanonicalize(la, nao, n_mo, nocc, fock_cur, coeff_work, fock_mo, energies_mo, info)
         if (info /= 0) then
            scf_result%message = "trc_scf: --soscf's semicanonicalisation did not converge "// &
                                 "(eigensolver info="//to_string(info)//")"
            exit
         end if
         call scf_soscf_gradient(n_mo, nocc, nvirt, fock_mo, gradient)

         ! The commutator at the current (accepted) point, the SAME metric
         ! and the same test trc_scf_run's own DIIS loop converges on.
         call cphf_copy_matrix(nao, nao, fock_cur, scf_state%fock(:, :, 1))
         call scf_fock_commutator_error(scf_context, scf_config, scf_state, scf_result, la)
         errmax = sqrt(la%dot_3d(scf_state%errv, scf_state%errv))
         de = abs(e_cur - e_old)

         !$acc update self(energies_mo)
         call scf_soscf_newton_step_rhf(scf_basis, scf_context%sph, scf_context%nao_c, eri, la, nao, n_mo, nocc, &
                                    nvirt, coeff_work, energies_mo, gradient, MAX_ROTATION, step, lowest, &
                                    predicted, products, DEFAULT_SOSCF_SUBSPACE)
         products_total = products_total + products
         builds = builds + products

         saddle = lowest < SADDLE_CURVATURE
         converged_now = de < scf_options%conv_energy .and. errmax < scf_options%conv_diis
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

         call soscf_pull_matrix(nvirt, nocc, step, step_host)
         scaling = maxval(abs(step_host))

         accepted = .false.
         do trial = 1, MAX_BACKTRACKS
            if (scaling > trust) then
               call soscf_scale_matrix(nvirt, nocc, step, trust/scaling, step_scaled)
            else
               call cphf_copy_matrix(nvirt, nocc, step, step_scaled)
            end if
            call soscf_embed_kappa(n_mo, nocc, nvirt, step_scaled, kappa)
            call scf_soscf_rotation_matrix(la, n_mo, min(scaling, trust), kappa, rotation, scratch_a, scratch_b, &
                                       scratch_c)
            call la%gemm_mat(coeff_work, rotation, coeff_try)
            call la%density(coeff_try, nocc, scf_config%occ, scf_result%dmat(:, :, 1))
            call scf_soscf_fock_rebuild(scf_basis, scf_context, scf_config, scf_state, scf_result, scf_energies, &
                                    la, eri)
            e_try = scf_energies%total
            builds = builds + 1

            ! `predicted` below its own resolution is a step too small to
            ! measure, not a step that failed -- see ENERGY_RESOLUTION.
            if (e_try < e_cur .or. predicted < ENERGY_RESOLUTION) then
               call cphf_copy_matrix(nao, n_mo, coeff_try, coeff_work)
               call cphf_copy_matrix(nao, nao, scf_state%fock(:, :, 1), fock_cur)
               call cphf_copy_matrix(nao, nao, scf_result%dmat(:, :, 1), dmat_accepted)
               e_old = e_cur
               e_cur = e_try
               trust = min(MAX_ROTATION, trust*TRUST_GROWTH)
               accepted = .true.
               exit
            end if
            trust = 0.5_dp*trust
            if (trust < MIN_ROTATION) exit
         end do

         if (scf_config%talk) then
            write (line, "(a,i0,a,f20.12,a,es9.2,a,es9.2,a,f7.4,a,es9.2,a,i0)") &
               "   soscf ", iteration, ": E = ", e_cur, "  dE ", de, "  |FDS-SDF| ", errmax, &
               "  trust ", trust, "  curvature ", lowest, "  fock builds ", builds
            call global_logger%info(trim(line))
         end if
         scf_result%iterations = switch_iteration + iteration

         if (.not. accepted) then
            call global_logger%warning("   the second-order SCF found no step downhill; stopping. "// &
                                       "the reported energy is the last one that descended.")
            exit
         end if
      end do

      ! `scf_result%dmat` may hold a REJECTED trial (backtracking exhausted,
      ! or the loop above hit `max_steps`): restore the last ACCEPTED
      ! density, and re-semicanonicalise so the orbital energies leaving
      ! are canonical (converged) or honestly semicanonical (stalled),
      ! never a half-written trial's.
      call cphf_copy_matrix(nao, nao, dmat_accepted, scf_result%dmat(:, :, 1))
      call scf_soscf_semicanonicalize(la, nao, n_mo, nocc, fock_cur, coeff_work, fock_mo, energies_mo, info)
      if (info == 0) then
         call cphf_copy_matrix(nao, n_mo, coeff_work, scf_result%cmo(:, 1:n_mo, 1))
         !$acc update self(energies_mo)
         scf_result%eps(1:n_mo, 1) = energies_mo
      end if
      if (scf_config%talk) then
         write (line, "(a,i0,a,i0,a)") "   second-order phase: ", products_total, &
            " Hessian-vector products (Fock builds) over ", max(scf_result%iterations - switch_iteration, 0), &
            " Newton iterations"
         call global_logger%info(trim(line))
      end if

      !$acc exit data delete(fock_cur, coeff_work, coeff_try, dmat_accepted, fock_mo, gradient, step) &
      !$acc           delete(step_scaled, kappa, rotation, scratch_a, scratch_b, scratch_c)
      !$acc exit data delete(energies_mo)
   end subroutine trc_scf_soscf_run

   !
   ! One full (non-incremental) Fock build and the energy at the density
   ! `scf_result%dmat(:,:,1)`, restricted Hartree-Fock only (no XC, no
   ! continuum -- `trc_scf_settings%resolve` has already refused both for
   ! `--soscf`): `scf_fock_build_g`/`scf_fock_energy`/`scf_fock_assemble`,
   ! the SAME sequence `trc_scf_run`'s own iteration runs, reused rather
   ! than duplicated. Leaves `scf_state%fock(:,:,1)` and
   ! `scf_energies%total` at that density; exported for
   ! `check_soscf.F90`'s own finite-difference validation.
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
   ! `coeff` (nao, n_mo) is rotated in place, within each block; what
   ! comes out is the Fock matrix in the new basis (`fock_mo`), whose
   ! occupied-virtual block is the gradient and whose diagonal
   ! (`energies_mo`) is the orbital energies -- canonical once the
   ! occupied-virtual block has vanished, semicanonical otherwise.
   !
   subroutine scf_soscf_semicanonicalize(la, nao, n_mo, nocc, fock_ao, coeff, fock_mo, energies_mo, info)
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: nao, n_mo, nocc
      real(dp), intent(in) :: fock_ao(nao, nao)
      real(dp), intent(inout) :: coeff(nao, n_mo)
      real(dp), intent(inout) :: fock_mo(n_mo, n_mo)
      real(dp), intent(inout) :: energies_mo(n_mo)
      integer, intent(out) :: info
      integer :: n_vir, k
      real(dp), allocatable :: work(:, :), rotated(:, :), block_matrix(:, :), block_values(:)

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
   ! (pulled once by the caller); `cmo`/`gradient`/`step` are device
   ! resident. `lowest` is an UPPER BOUND on the true smallest curvature
   ! (a Rayleigh-Ritz value), not the curvature itself.
   !
   subroutine scf_soscf_newton_step_rhf(b, sph, nao_c, eri, lp, nao, n_mo, nocc, nvirt, cmo, energies_mo, gradient, &
                                    escape, step, lowest, predicted, products, max_subspace)
      type(trc_basis_t), intent(in) :: b
      logical, intent(in) :: sph
      integer, intent(in) :: nao_c
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, n_mo, nocc, nvirt
      real(dp), intent(in) :: cmo(nao, n_mo)
      real(dp), intent(in) :: energies_mo(n_mo)
      real(dp), intent(in) :: gradient(nvirt, nocc)
      real(dp), intent(in) :: escape
      real(dp), intent(inout) :: step(nvirt, nocc)
      real(dp), intent(out) :: lowest, predicted
      integer, intent(out) :: products
      integer, intent(in) :: max_subspace

      real(dp), allocatable :: gaps(:, :)
      real(dp), allocatable :: basis(:, :, :), image(:, :, :)
      real(dp), allocatable :: direction(:, :), residual(:, :)
      real(dp), allocatable :: d_u(:, :), d_ut(:, :), g_u(:, :), d_u_c(:, :), g_u_c(:, :)
      real(dp), allocatable :: temporary_no(:, :), temporary_no2(:, :)
      real(dp), allocatable :: small(:, :), small_gradient(:), amplitude(:)
      real(dp) :: gradient_norm, resid_norm
      integer :: nmax, nsub, i, j, a0, i0
      integer :: softest(2)
      logical :: added, ok

      lowest = 0.0_dp
      predicted = 0.0_dp
      products = 0
      call cphf_zero_matrix(nvirt, nocc, step)
      if (nvirt < 1 .or. nocc < 1) return

      allocate (gaps(nvirt, nocc))
      do i = 1, nocc
         do j = 1, nvirt
            gaps(j, i) = energies_mo(nocc + j) - energies_mo(i)
         end do
      end do
      !$acc enter data copyin(gaps)

      nmax = max(1, min(max_subspace, nvirt*nocc))
      allocate (basis(nvirt, nocc, nmax), image(nvirt, nocc, nmax))
      allocate (direction(nvirt, nocc), residual(nvirt, nocc))
      allocate (d_u(nao, nao), d_ut(nao, nao), g_u(nao, nao), d_u_c(nao_c, nao_c), g_u_c(nao_c, nao_c))
      allocate (temporary_no(nao, nocc), temporary_no2(nao, nocc))
      !$acc enter data create(basis, image, direction, residual, d_u, d_ut, g_u, d_u_c, g_u_c) &
      !$acc            create(temporary_no, temporary_no2)

      gradient_norm = sqrt(lp%dot_2d(gradient, gradient))

      ! ---- the seeds: the preconditioned gradient, and the softest mode --
      nsub = 0
      call soscf_preconditioned_direction(nvirt, nocc, gaps, PRECONDITION_FLOOR, gradient, direction)
      call soscf_add_direction(lp, nvirt, nocc, basis, nsub, direction, added)
      softest = minloc(gaps)
      a0 = softest(1)
      i0 = softest(2)
      call soscf_unit_matrix(nvirt, nocc, a0, i0, direction)
      call soscf_add_direction(lp, nvirt, nocc, basis, nsub, direction, added)

      if (nsub > 0) then
         do i = 1, nsub
            call scf_soscf_hessian_apply_rhf(b, sph, eri, lp, nao, nao_c, n_mo, nocc, nvirt, cmo, gaps, &
                                             basis(:, :, i), image(:, :, i), d_u, d_ut, g_u, d_u_c, g_u_c, &
                                             temporary_no, temporary_no2)
            call soscf_scale_matrix(nvirt, nocc, image(:, :, i), HESSIAN_SCALE, image(:, :, i))
            products = products + 1
         end do

         ! ---- solve, expand, repeat --------------------------------------
         do
            if (allocated(small)) deallocate (small, small_gradient)
            allocate (small(nsub, nsub), small_gradient(nsub))
            do j = 1, nsub
               do i = 1, nsub
                  small(i, j) = lp%dot_2d(basis(:, :, i), image(:, :, j))
               end do
               small_gradient(j) = lp%dot_2d(basis(:, :, j), gradient)
            end do
            small = 0.5_dp*(small + transpose(small))

            if (allocated(amplitude)) deallocate (amplitude)
            allocate (amplitude(nsub))
            call soscf_level_shifted_step(small, small_gradient, escape, amplitude, lowest, predicted, ok)
            if (.not. ok) then
               call global_logger%warning("trc_scf_soscf: the projected orbital Hessian failed to diagonalise")
               exit
            end if

            call cphf_zero_matrix(nvirt, nocc, step)
            call cphf_copy_matrix(nvirt, nocc, gradient, residual)
            do i = 1, nsub
               call cphf_axpy_matrix(nvirt, nocc, amplitude(i), basis(:, :, i), step)
               call cphf_axpy_matrix(nvirt, nocc, amplitude(i), image(:, :, i), residual)
            end do

            if (nsub >= nmax) exit
            resid_norm = sqrt(lp%dot_2d(residual, residual))
            if (resid_norm <= SUBSPACE_TOLERANCE*max(gradient_norm, tiny(1.0_dp))) exit

            call soscf_preconditioned_direction(nvirt, nocc, gaps, PRECONDITION_FLOOR, residual, direction)
            call soscf_add_direction(lp, nvirt, nocc, basis, nsub, direction, added)
            if (.not. added) exit

            call scf_soscf_hessian_apply_rhf(b, sph, eri, lp, nao, nao_c, n_mo, nocc, nvirt, cmo, gaps, &
                                             basis(:, :, nsub), image(:, :, nsub), d_u, d_ut, g_u, d_u_c, g_u_c, &
                                             temporary_no, temporary_no2)
            call soscf_scale_matrix(nvirt, nocc, image(:, :, nsub), HESSIAN_SCALE, image(:, :, nsub))
            products = products + 1
         end do
      else
         ! Unreachable for nvirt, nocc >= 1 (the second seed is always a
         ! unit vector, and the first call to soscf_add_direction has
         ! nothing yet to project it against) -- defensive only.
         call global_logger%warning("trc_scf_soscf: no direction to build a Krylov subspace from")
      end if

      if (allocated(small)) deallocate (small, small_gradient)
      if (allocated(amplitude)) deallocate (amplitude)
      !$acc exit data delete(basis, image, direction, residual, d_u, d_ut, g_u, d_u_c, g_u_c) &
      !$acc           delete(temporary_no, temporary_no2)
      deallocate (basis, image, direction, residual, d_u, d_ut, g_u, d_u_c, g_u_c, temporary_no, temporary_no2)
      !$acc exit data delete(gaps)
      deallocate (gaps)
   end subroutine scf_soscf_newton_step_rhf

   !
   ! `(e_a-e_i) p + [A p]`, the RHF orbital-Hessian action on a trial
   ! rotation `pin`, folded through the Cartesian space `eri%fock_resident`
   ! runs in -- see this module's own header. `d_u`/`d_ut`/`g_u` are
   ! `nao`-by-`nao` (the SCF's own space, spherical when the basis is);
   ! `d_u_c`/`g_u_c` are `nao_c`-by-`nao_c` (always Cartesian); `sph =
   ! .false.` makes `scf_fock_to_cart`/`scf_fock_to_sph` a plain copy, so a
   ! Cartesian basis (`nao == nao_c`) pays nothing beyond it.
   !
   subroutine scf_soscf_hessian_apply_rhf(b, sph, eri, lp, nao, nao_c, n_mo, nocc, nvirt, cmo, diagv, pin, apout, &
                                          d_u, d_ut, g_u, d_u_c, g_u_c, temporary_no, temporary_no2)
      type(trc_basis_t), intent(in) :: b
      logical, intent(in) :: sph
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nao_c, n_mo, nocc, nvirt
      real(dp), intent(in) :: cmo(nao, n_mo)
      real(dp), intent(in) :: diagv(nvirt, nocc)
      real(dp), intent(in) :: pin(nvirt, nocc)
      real(dp), intent(inout) :: apout(nvirt, nocc)
      real(dp), intent(inout) :: d_u(nao, nao), d_ut(nao, nao), g_u(nao, nao)
      real(dp), intent(inout) :: d_u_c(nao_c, nao_c), g_u_c(nao_c, nao_c)
      real(dp), intent(inout) :: temporary_no(nao, nocc), temporary_no2(nao, nocc)
      call lp%gemm_mat(cmo(:, nocc + 1:n_mo), pin, temporary_no)
      call lp%gemm_mat(temporary_no, cmo(:, 1:nocc), d_u, transb="T")
      call lp%transpose(nao, nao, d_u, nao, d_ut, nao)
      call cphf_add_matrix(nao, nao, d_u, d_ut, d_u)
      call scf_fock_to_cart(b, sph, nao, nao_c, d_u, d_u_c)
      call eri%fock_resident(b, d_u_c, g_u_c, k_scale=1.0_dp, j_scale=1.0_dp, density_screen=.false.)
      call scf_fock_to_sph(b, sph, nao, nao_c, g_u_c, g_u)
      call lp%gemm_mat(g_u, cmo(:, 1:nocc), temporary_no2)
      call lp%gemm_mat(cmo(:, nocc + 1:n_mo), temporary_no2, apout, transa="T", alpha=2.0_dp)
      call cphf_add_diag(nvirt, nocc, diagv, pin, apout)
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

   !
   ! `exp(kappa)` for antisymmetric `kappa`, by scaling and squaring:
   ! orthogonal to machine precision for any step size (ported from
   ! metalquicha's `mqc_orbital_rotation%rotation_matrix`, same algorithm,
   ! device-resident GEMMs). `norm` is `maxval(abs(kappa))`, passed in
   ! rather than recomputed here: the caller already knows it as the
   ! (host) trust-radius scaling factor, and pulling the full n x n
   ! `kappa` back to the host just to recompute the same number would be
   ! the one unnecessary device-to-host transfer this module otherwise
   ! avoids. `scaled`/`term`/`next_term` are the caller's own n x n
   ! scratch.
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
