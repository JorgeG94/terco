!
!
! The SCF: a basis and an electron count in, a converged energy and
! density out, with every matrix resident on the device. Hartree-Fock or
! Kohn-Sham (any functional in trc_xc_functional), restricted or unrestricted.
!
! WHAT CROSSES THE BUS
! --------------------
! Per iteration: two energies, a norm, and a DIIS matrix of at most ten by
! ten. The density, Fock matrices, orbitals, DIIS history and orthogonaliser
! stay on the device; the Fock build is the resident driver, the XC integrator
! finds its density already present, and the products, commutator and
! eigenproblem are trc_linalg (cuBLAS/cuSOLVER under OpenACC). This follows
! metalquicha's device SCF, which measured host LAPACK plus copies at three
! times the cost of the GPU Fock build. A host build runs the same source with
! LAPACK, and the numbers agree to rounding, so the host suite tests this module.
!
! THE ITERATION
! -------------
! A few damped steps from the guess, then Pulay's DIIS on the commutator
! FDS - SDF in the orthogonal basis, per spin with one set of coefficients.
! Converged when the energy change and the RMS density change are both below
! threshold. The DIIS system (at most 11 by 11) is solved by pivoted
! elimination as in metalquicha's mqc_diis, with the error block scaled to
! order one against the constraint border first: near convergence the overlaps
! are 1e-16 and an unscaled elimination makes the weights noise, showing as a
! density that drifts from the fixed point while the energy sits flat.
!
! UNRESTRICTED
! ------------
! F_s = H + J[D_a + D_b] - a_x K[D_s] + V_xc,s. The Fock driver returns
! j J[D] - k K[D]/2 for a density it takes to be a total, so J is one call on
! D_a + D_b with k = 0 and -a_x K[D_s] one call on 2 D_s with j = 0, k = a_x:
! three builds per iteration against one restricted. A batched build over the
! three densities is the obvious improvement and is not done.
!
! SPHERICAL BASES
! ---------------
! `b%spherical` (trc_api.F90, from a JSON basis's function_type) changes the
! size of every matrix this module owns (S, H, D, F, C, eps, X, the DIIS
! history) from `b%nao` to `b%nao_sph`; everything operating on them is
! unchanged, just at a smaller `nao`. Everything touching a kernel stays
! Cartesian: `trc_1e`, `trc_eri_t%fock_resident`, and the XC integrator.
! Cartesian-sized scratch matrices (suffix `_c`: `s_c`, `h_c`, `d_c`,
! `gmat_c`, `vxc_c`, ...) sit between the kernels and the SCF's own matrices,
! built once per geometry (`s_c`/`h_c`, now `trc_scf_context`'s `smat`/`hcore`)
! or once per iteration (`d_c` from the current density; `gmat_c`/`vxc_c` from
! the Fock/XC calls, now `trc_scf_state`'s), and folded by `sph_cart_to_sph`.
! `sph` guards every call into trc_sph.F90, so the Cartesian path is the same
! code path as before spherical support existed.
!
! The XC integrator gets `d_c`, the Cartesian density rebuilt from `res%dmat`
! each iteration, not spherical AO values on the grid: the transform is shared
! with the Coulomb/exchange build, and `collocation_shell` evaluates Cartesian
! components straight off a flat table, so a second AO space there would
! duplicate trc_sph.F90's block structure.
!
! SAD and GWH reach this module already in the right space (SAD's atomic
! sub-SCFs run on a `subset` basis carrying `spherical`; GWH is built here from
! `hcore`). SAP stays Cartesian and is transformed once at the trc_capi.F90
! call site, as it has no basis handle to carry `spherical`.
!
! POLARISABLE CONTINUUM
! ---------------------
! With `pcm_model` set, trc_pcm builds the cavity, K and the nuclear potential
! once. Every iteration, at the density the energy is evaluated at, the surface
! potential of the total Cartesian density gives the charges, q.v/2 into the
! energy and the charges' full operator into each spin's Fock matrix (as
! pyscf.solvent.pcm's SCFWithSolvent, which also adds it before DIIS). It is
! rebuilt in full each iteration, independent of the incremental-ERI path.
!
! MODULES
! -------
! `trc_scf_settings` resolves `trc_scf_options_t` into what the loop reads
! (`scf_settings_t`, "scf_config" below); `trc_scf_context` builds everything once
! per geometry (`scf_ctx_t`, "scf_context": the pair list, S/H, the ERI engine, the
! linear-algebra handles, the orthogonaliser, and -- when asked for -- the
! PCM cavity and the DFT grid); `trc_scf_state` owns the per-iteration arrays
! (`scf_state_t`, "scf_state") and the one enter/exit pair that puts them, `scf_context`'s
! own matrices and `scf_result` on the device together. `scf_fock_build_g`/`scf_fock_assemble`
! and the DIIS push/combine helpers below stay in THIS file for now (they
! are not yet split into `trc_scf_fock`/`trc_scf_diis`): every one of them
! already takes its arrays as explicit-shape dummies rather than touching
! `scf_context`/`scf_state` inside a `do concurrent` body, which is what lets this compile
! correctly once those arrays are struct components instead of plain locals.
!
module trc_scf_driver
   use pic_types, only: dp, int64
   use trc_basis_types, only: trc_basis_t
   use trc_error, only: error_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_settings, only: trc_scf_options_t, trc_scf_result_t, scf_settings_t, scf_energies_t, &
                               scf_settings_homo_lumo_gap
   use trc_scf_context, only: scf_ctx_t, scf_context_bcast_flat, trc_grid_cache_t, trc_pair_eri_cache_t, scf_context_project_density
   use trc_scf_state, only: scf_state_t, scf_timers_t
   use trc_scf_fock, only: scf_fock_build_g, scf_fock_build_vxc, scf_fock_build_pcm, scf_fock_energy, scf_fock_assemble, &
                           scf_fock_commutator_error, scf_fock_to_cart, scf_fock_to_sph
   use trc_scf_df, only: df_jk_t
   use trc_scf_diis, only: cdiis_t, adiis_t, scf_diis_level_shift, scf_diis_damp
   use pic_types, only: default_int
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_mpi_lib, only: comm_t
   implicit none
   private

   public :: trc_scf_options_t, trc_scf_result_t, trc_scf_run, trc_grid_cache_t, trc_pair_eri_cache_t

contains

   subroutine dump_density(nao, d)
      !! The density at a chosen iteration, written where a benchmark can read
      !! it back. Off unless TRC_DUMP_DENSITY names a file. Stream format:
      !! an 8-byte nao, then the matrix, which is what the test readers expect.
      integer, intent(in) :: nao
      real(dp), intent(in) :: d(nao, nao)
      character(len=512) :: path
      integer :: u, ios
      integer(int64) :: prev
      call get_environment_variable("TRC_DUMP_DENSITY", path)
      if (len_trim(path) == 0) return
      ! The atomic SCFs of a SAD guess come through here too, and would
      ! overwrite the molecular density with a 5x5 one. Keep the largest.
      open (newunit=u, file=trim(path), access="stream", form="unformatted", &
            status="old", action="read", iostat=ios)
      if (ios == 0) then
         read (u, iostat=ios) prev
         close (u)
         if (ios == 0) then
            if (prev >= int(nao, int64)) return
         end if
      end if
      open (newunit=u, file=trim(path), access="stream", form="unformatted", &
            status="replace", action="write", iostat=ios)
      if (ios /= 0) return
      write (u) int(nao, int64)
      write (u) d
      close (u)
      call global_logger%info("density dumped to "//trim(path)//", nao "//to_string(nao))
   end subroutine dump_density

   !
   ! Run the SCF. `nalpha` and `nbeta` are the electron counts; the case is
   ! restricted when they are equal and `unrestricted` is not forced. The
   ! basis must already be on the device. `dguess` is (nao, nao, nspin) in
   ! the result's convention; without it the guess is the core Hamiltonian.
   !
   subroutine trc_scf_run(scf_basis, nalpha, nbeta, scf_options, scf_result, dguess, comm, hguess, la, scf_dfjk, &
                          grid_cache, pair_eri_cache, cmo_project)
      type(trc_basis_t), intent(in) :: scf_basis
      integer, intent(in) :: nalpha, nbeta
      type(trc_scf_options_t), intent(in) :: scf_options
      type(trc_scf_result_t), intent(out) :: scf_result
      real(dp), intent(in), optional :: dguess(:, :, :)
      real(dp), intent(in), optional :: hguess(:, :)
      !! A guess HAMILTONIAN to occupy instead of hcore or the GWH matrix --
      !! how SAP gets in (see trc_sap), which cannot come through `dguess`
      !! because what it produces is a one-electron operator and not a
      !! density. Ignored if `dguess` is also given: a density the caller
      !! already has is more specific than a matrix to make one from.
      real(dp), intent(in), optional :: cmo_project(:, :, :)
      !! A PREVIOUS geometry's converged orbitals (SAME molecule/basis, so
      !! the same nao), re-orthonormalised into a guess density at THIS
      !! geometry's own overlap once `scf_context%setup` below has built
      !! it (`trc_scf_context%scf_context_project_density`) -- trc_capi.F90's
      !! TRC_GUESS_PROJECTED, lowest priority of the three guess-density
      !! arguments (checked after `dguess` and `hguess`, same as every
      !! other guess kind the plain `scf_options%guess` string selects).
      type(comm_t), intent(in), optional :: comm
      !! Ranks. Every rank runs the same iteration on the same matrices;
      !! the quartets are split across them and the Fock matrix summed
      !! back, so the density stays identical everywhere. Bind devices
      !! first (trc_bind_device). The XC batches are split across ranks too.
      type(trc_linalg_t), intent(inout), optional, target :: la
      !! cuBLAS/cuSOLVER handles made once by the caller and carried through
      !! the calculation (trc_linalg's init is create-once and grows the
      !! workspace to nao). Absent, the run makes its own and releases them.
      type(df_jk_t), intent(inout), optional, target :: scf_dfjk
      !! The density-fitted J/K object for `scf_options%df` (`--df`), already fitted
      !! by the caller -- required when `scf_options%df` is true (trc_scf_context's
      !! `setup` refuses otherwise), ignored when it is false. Owned by the
      !! caller across more than one run (trc_capi.F90's context handle), the
      !! same reason `la` is a caller-owned, optional argument rather than
      !! something this routine builds for itself.
      type(trc_grid_cache_t), intent(inout), optional, target :: grid_cache
      !! A caller-owned DFT grid cache -- forwarded verbatim to
      !! `scf_context%setup` (trc_scf_context.F90's own header on the
      !! argument). Present, the converged grid stays alive (and device-
      !! resident) in `grid_cache` after this call returns, for a
      !! gradient run right after to reuse (trc_capi.F90's
      !! `trc_run_gradient`); absent, the run's own grid is released with
      !! `scf_context` as before this argument existed.
      type(trc_pair_eri_cache_t), intent(inout), optional, target :: pair_eri_cache
      !! A caller-owned pair-list/ERI cache -- forwarded verbatim to
      !! `scf_context%setup` (trc_scf_context.F90's `trc_pair_eri_cache_t`
      !! own header). Present, the pair list and (outside `--df`) the
      !! four-centre engine stay alive in it after this call returns, for a
      !! gradient run right after at the SAME geometry and threshold to
      !! reuse (trc_capi.F90's `trc_run_gradient`) instead of rebuilding
      !! both; absent, the run's own copies are released with
      !! `scf_context` as before this argument existed.

      integer :: nao, nao_c, s, it, i, j, n2
      logical :: ok, sph
      type(scf_settings_t) :: scf_config
      type(scf_ctx_t), target :: scf_context
      !! TARGET: scf_context%setup points scf_context%lp at scf_context%la_own when the caller gives
      !! no linalg handles of its own; that pointer must stay valid for the
      !! rest of this routine, which needs scf_context itself to carry TARGET.
      type(scf_state_t) :: scf_state
      type(scf_timers_t) :: scf_timers
      !! This iteration's share of the two totals (tf_it/t_fock, tx_it/t_xc,
      !! tr_it/t_rest), so the table says where the time went while it is
      !! going rather than only at the end. A Fock build is not the same
      !! cost twice: the SAD guess is atomic, so density screening throws
      !! nearly all of iteration 1 away, and the per-iteration cost then
      !! climbs to a plateau. Comparing two runs by total/iterations is
      !! therefore only valid at equal iteration counts.
      character(len=256) :: msg
      character(len=160) :: line
      logical :: incr_now
      type(cdiis_t) :: cdiis
      type(adiis_t) :: adiis
      !! ADIIS's history (fock, density) is kept apart from CDIIS's own (fock,
      !! commutator): ADIIS runs exactly in the window CDIIS has not opened
      !! yet, so the two never need the same entries, and keeping them
      !! separate means `scf_options%adiis = .false.` leaves the CDIIS path
      !! untouched, to the bit.
      logical :: adiis_active
      type(trc_linalg_t), pointer :: linalg_provider
      type(error_t) :: err
      type(scf_energies_t) :: scf_energies
      !! Every energy and integrated electron count the loop produces; all zero
      !! until the first iteration, so a run with max_iteration=0 still reports zeros.
      logical :: diis_on
      real(dp) :: ctl(2)
      real(dp) :: eold, drms, errmax
      integer :: nfrac
      logical :: diag_failed      !! this iteration's Fock diagonalisation (every rank's own copy)
      integer :: diag_info        !! the FIRST failing spin's own eigensolver info, kept across the spin loop
      integer :: info_s           !! this spin's own raw info from the call just made
      real(dp) :: diag_flag(2)    !! [0.0/1.0 failed, info], broadcast from root the same way `ctl` is

      ! Everything `scf_options` resolves to before the loop starts: the functional
      ! lookup and its rejections, the spin rule, the by-spin damping/DIIS
      ! defaults, and the env-hook overrides (trc_scf_settings). `scf_config` is
      ! immutable from here on, and the rest of this routine reads it
      ! directly.
      call scf_timers%run%start()
      ok = scf_config%resolve(scf_options, nalpha, nbeta, msg=msg, comm=comm)
      if (.not. ok) then
         scf_result%message = msg
         return
      end if
      scf_result%nspin = scf_config%nspin
      ! The cavity, K and its LU are built and solved whole on every rank
      ! (replicated: cheap next to the integrals); pcm%operator splits the
      ! two integral passes over the ranks itself, given `comm`.

      ! --- everything once per geometry: pairs, 1e, S/H fold, ERI, linalg,
      ! orthogonaliser, and (when asked for) the PCM cavity and DFT grid ---
      call scf_context%setup(scf_basis, scf_options, scf_config, err, comm=comm, la=la, dfjk=scf_dfjk, &
                             grid_cache=grid_cache, pair_eri_cache=pair_eri_cache)
      if (err%has_error()) then
         scf_result%message = "trc_scf: "//err%get_message()
         return
      end if
      nao = scf_context%nao
      nao_c = scf_context%nao_c
      sph = scf_context%sph
      n2 = nao*nao
      linalg_provider => scf_context%lp
      scf_result%e_nuc = scf_context%e_nuc
      scf_result%n_lindep = scf_context%n_lindep
      if (scf_config%pcm_on) scf_result%pcm_npoints = scf_context%pcm%surface%n_points
      if (scf_config%dft) then
         scf_result%grid_r = scf_context%gc%grid%coords
         scf_result%grid_w = scf_context%gc%grid%weights
         if (scf_config%nlc) then
            scf_result%grid_r_nlc = scf_context%gc%grid_nlc%coords
            scf_result%grid_w_nlc = scf_context%gc%grid_nlc%weights
         end if
         ! The DFT pipeline's own step table (grid build, semilocal XC,
         ! VV10), reported alongside the existing fock(s)/xc(s) columns
         ! below rather than instead of them: `report`'s per-step rate
         ! columns are new here (grid build has no flop count to show
         ! them against), but the three totals are not new numbers, only
         ! a second, uniform place they land -- see trc_scf_state's own
         ! header on `xcst`. "dft", not "scf": its total is grid + XC only.
         call scf_timers%xcst%init("dft", [character(len=8) :: "grid", "xc", "vv10"])
         call scf_timers%xcst%add_time(1, scf_context%t_grid)
      end if

      ! ONLY WHAT THIS RUN USES IS FULL SIZE (trc_scf_state%alloc owns the
      ! array sizes, 1 x 1 for an array the run never touches). One that a
      ! run never touches is nao^2 of device memory for nothing: on ionic05
      ! in aug-cc-pVTZ (4255 functions, 4950 Cartesian) an HF run carried
      ! ~2.4 GB of XC, unrestricted, range-separated, incremental and PCM
      ! matrices it had no use for.
      call scf_state%alloc(scf_context, scf_config)
      call cdiis%alloc(nao, scf_config%nspin, scf_options%ndiis)
      allocate (scf_result%dmat(nao, nao, scf_config%nspin), scf_result%cmo(nao, nao, scf_config%nspin), &
             scf_result%eps(nao, scf_config%nspin))
      ! Only when asked for: a second ndiis-deep history the size of
      ! cdiis's own, which matters on a system the size of w172, where
      ! one nao^2*ndiis block is already hundreds of megabytes.
      if (scf_options%adiis) call adiis%alloc(nao, scf_config%nspin, scf_options%ndiis)

      ! --- the guess, on the host, then everything goes up -------------------
      diag_failed = .false.
      diag_info = 0
      if (present(dguess)) then
         scf_result%dmat = dguess
      else if (present(hguess)) then
         do s = 1, scf_config%nspin
            call guess_from(hguess, scf_config%nocc(s), scf_result%dmat(:, :, s), info_s)
            if (info_s /= 0 .and. .not. diag_failed) diag_info = info_s
            if (info_s /= 0) diag_failed = .true.
         end do
      else if (present(cmo_project)) then
         call scf_context_project_density(linalg_provider, nao, scf_config%nspin, scf_config%nocc, cmo_project, &
                              scf_context%smat, scf_result%dmat)
      else if (scf_options%guess == "core") then
         do s = 1, scf_config%nspin
            call guess_from(scf_context%hcore, scf_config%nocc(s), scf_result%dmat(:, :, s), info_s)
            if (info_s /= 0 .and. .not. diag_failed) diag_info = info_s
            if (info_s /= 0) diag_failed = .true.
         end do
      else
         ! GWH: F_ij = 1.75/2 S_ij (H_ii + H_jj), diagonal H_ii.
         do j = 1, nao
            do i = 1, nao
               scf_state%w1(i, j) = 0.875_dp*scf_context%smat(i, j)*(scf_context%hcore(i, i) + scf_context%hcore(j, j))
            end do
            scf_state%w1(j, j) = scf_context%hcore(j, j)
         end do
         do s = 1, scf_config%nspin
            call guess_from(scf_state%w1, scf_config%nocc(s), scf_result%dmat(:, :, s), info_s)
            if (info_s /= 0 .and. .not. diag_failed) diag_info = info_s
            if (info_s /= 0) diag_failed = .true.
         end do
      end if
      ! A diagonalisation failure in the INITIAL guess (core/GWH/a supplied
      ! Hamiltonian -- SAD/SADQ build their own density without this
      ! routine, through `dguess`): report it the same way `scf_ctx_setup`
      ! itself does on its own refusals (iterations is still 0 here, so
      ! trc_capi.F90's trc_run_scf reports TRC_ERR_UNSUPPORTED, not
      ! TRC_ERR_NOCONV -- this never started, rather than failed to finish).
      if (diag_failed) then
         scf_result%message = "trc_scf: the initial guess's own diagonalisation did not converge "// &
                              "(eigensolver info="//to_string(diag_info)//")"
         call scf_context%release()
         return
      end if
      scf_result%cmo = 0.0_dp
      scf_result%eps = 0.0_dp

      call scf_state%to_device(scf_context, scf_result)
      call cdiis%to_device()
      if (scf_options%adiis) call adiis%to_device()

      eold = 0.0_dp
      ! Read by the first iteration's `incr_now` before any commutator exists
      ! (Fortran's .and. does not short-circuit), so give it a value.
      errmax = huge(1.0_dp)
      diis_on = .false.
      ! What the SCF is actually running on. Two callers can hand over the
      ! same molecule and basis name and different shells -- one had a
      ! reader that kept dead primitives -- and that only shows here.
      if (scf_config%talk) call global_logger%info("   basis: "//to_string(nao)//" functions, " &
         //to_string(scf_basis%nshell)//" shells, "//to_string(sum(scf_basis%sh_np))//" primitives, max " &
         //to_string(maxval(scf_basis%sh_np)))
      ! Say so when increments are on: they are opt-in and change the timings.
      ! An `i` after an iteration's line marks a Fock matrix built from the
      ! density difference.
      if (scf_config%talk .and. scf_config%incr_on) then
         write (line, "(a,es8.1,a,i0,a)") "   incremental Fock: from |FDS-SDF| < ", scf_config%incr_start_run, &
            ", full build every ", max(1, scf_config%incr_reset_run), " iterations (marked i)"
         call global_logger%info(trim(line))
      end if
      if (scf_config%talk) then
         if (scf_config%dft) then
            call global_logger%info("   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s)     xc(s) linalg(s)")
         else
            call global_logger%info("   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s) linalg(s)")
         end if
      end if
      do it = 1, scf_options%max_iteration
         call scf_timers%fock%start()
         incr_now = .false.
         ! --- D_cart := T D_sph T^T, the one transform every Cartesian
         ! quantity below (the Fock build, the XC density) is built from.
         ! Not taken at all when `sph` is false -- see this module's own
         ! header on why that keeps the Cartesian path byte-for-byte
         ! unchanged rather than merely numerically equivalent.
         do s = 1, scf_config%nspin
            call scf_fock_to_cart(scf_basis, sph, nao, nao_c, scf_result%dmat(:, :, s), scf_state%d_c(:, :, s))
         end do
         ! --- two-electron part, resident, entirely in Cartesian space -------
         !
         ! INCREMENTAL. G(D) = G(D_ref) + G(D - D_ref), exactly, because G is
         ! linear in D. The point is the screening: the kernel screens on the
         ! density it is handed, and the difference is small and getting
         ! smaller while the density is neither. See `scf_fock_build_g`'s own header.
         !
         incr_now = scf_config%incr_on .and. it > 1 &
                    .and. errmax < scf_config%incr_start_run &
                    .and. mod(it - 1, max(1, scf_config%incr_reset_run)) /= 0
         call scf_fock_build_g(scf_basis, scf_context, scf_config, scf_state, incr_now)
         ! --- G_sph := T^T G_cart T, back up to the space the rest of the
         ! iteration runs in.
         do s = 1, scf_config%nspin
            call scf_fock_to_sph(scf_basis, sph, nao, nao_c, scf_state%gmat_c(:, :, s), scf_state%gmat(:, :, s))
         end do
         scf_timers%tf_it = scf_timers%fock%get_elapsed_time()
         scf_timers%t_fock = scf_timers%t_fock + scf_timers%tf_it
         call scf_timers%xc%start()
         ! --- exchange-correlation: finds its (Cartesian) density present ----
         if (scf_config%dft) then
            call scf_fock_build_vxc(scf_basis, scf_context, scf_config, scf_state, scf_options%rho_tol, &
                                    scf_energies, scf_timers, comm)
         end if
         scf_timers%tx_it = scf_timers%xc%get_elapsed_time()
         scf_timers%t_xc = scf_timers%t_xc + scf_timers%tx_it
         call scf_timers%rest%start()
         ! --- the continuum, at the TOTAL density (pyscf's choice for UKS
         ! too: one set of charges, one operator added to both spins) -----
         if (scf_config%pcm_on) then
            call scf_timers%pcm%start()
            call scf_fock_build_pcm(scf_basis, scf_context, scf_config, scf_state, scf_energies, comm)
            scf_timers%t_pcm = scf_timers%t_pcm + scf_timers%pcm%get_elapsed_time()
            ! the PCM step is its own column, not part of the rest
            call scf_timers%rest%start()
         end if
         ! --- energy at this density -----------------------------------------
         call scf_fock_energy(scf_context, scf_config, scf_state, scf_result, scf_energies)
         ! --- Fock and its commutator, per spin ------------------------------
         call scf_fock_assemble(scf_context, scf_config, scf_state)
         call scf_fock_commutator_error(scf_context, scf_config, scf_state, scf_result)
         ! `dold` must hold `scf_result%dmat` from BEFORE DIIS/diagonalisation touch
         ! it -- see trc_scf_state's own header.
         call scf_state%save_dold(nao, scf_config%nspin, scf_result%dmat)
         ! --- DIIS, once the commutator is small enough to extrapolate from --
         errmax = sqrt(linalg_provider%dot_3d(scf_state%errv, scf_state%errv))
         ! Over ranks the Fock matrices agree only to the XC atomics' 1e-13,
         ! so the numbers every branch below and the convergence test read
         ! are rank 0's, on every rank: one rank leaving the loop while
         ! another iterates on would hang the next collective.
         if (present(comm)) then
            if (comm%size() > 1) then
               ctl = [scf_energies%total, errmax]
               call scf_context_bcast_flat(comm, ctl, 2)
               scf_energies%total = ctl(1)
               errmax = ctl(2)
            end if
         end if
         if (.not. diis_on .and. it > scf_config%ndamp) diis_on = errmax < scf_config%diis_start
         !
         ! A/EDIIS while damping has ended but the commutator is still too large for
         ! CDIIS. CDIIS has no constraint on its coefficients, hence the
         ! `diis_start` gate: far from convergence it can overshoot. ADIIS (Hu &
         ! Yang, JCP 132, 054109) minimises an energy functional over the convex
         ! hull of the history, so it only moves towards a visited density: reliable
         ! where CDIIS is not and slow where CDIIS is fast, hence the two switch
         ! rather than blend (gmshpc: a blended step turned a stretched-water case
         ! converging in 21 CDIIS iterations into one that never converged).
         !
         adiis_active = scf_options%adiis .and. .not. diis_on .and. it > scf_config%ndamp
         ! The new density is rank 0's alone: extrapolation, eigensolve and density
         ! run on one device and the density is broadcast, since dsyevd need not
         ! return the same vectors on two devices. Orbitals and energies come over
         ! once, after the loop.
         diag_failed = .false.
         diag_info = 0
         if (scf_config%root) then
            if (adiis_active) then
               call adiis%push(linalg_provider, nao, scf_config%nspin, scf_state%fock, scf_state%dold)
               call adiis%extrapolate(nao, scf_config%nspin, scf_state%fock)
            end if
            ! --- level shift while DIIS and A/EDIIS are off: F += mu (S - S D S / occ)
            if (.not. diis_on .and. .not. adiis_active .and. scf_options%level_shift > 0.0_dp) then
               call scf_diis_level_shift(linalg_provider, nao, scf_config%nspin, scf_options%level_shift, scf_config%occ, &
                      scf_context%smat, scf_result%dmat, scf_state%w1, scf_state%w2, scf_state%fock)
            end if
            if (diis_on) then
               call cdiis%push(linalg_provider, nao, scf_config%nspin, scf_state%fock, scf_state%errv)
               call cdiis%extrapolate(nao, scf_config%nspin, scf_state%fock, ok)
            end if
            ! --- diagonalise, per spin, and the new density ---------------------
            do s = 1, scf_config%nspin
               call orbitals_and_density(scf_state%fock(:, :, s), scf_config%nocc(s), scf_options%frac_occ, &
                                         scf_state%w1, scf_state%w2, scf_result%eps(:, s), scf_result%cmo(:, :, s), &
                                         scf_result%dmat(:, :, s), scf_state%ofrac, info_s)
               ! The FIRST failing spin's own info, not whichever spin ran
               ! last: an unrestricted run that fails on alpha and then
               ! succeeds on beta must still report alpha's own info, not
               ! overwrite it with beta's 0.
               if (info_s /= 0 .and. .not. diag_failed) diag_info = info_s
               if (info_s /= 0) diag_failed = .true.
            end do
            if (.not. diag_failed .and. .not. diis_on .and. .not. adiis_active) then
               call scf_diis_damp(nao, scf_config%nspin, scf_options%damp, scf_state%dold, scf_result%dmat)
            end if
         end if
         ! Every rank must learn of a diagonalisation failure the same
         ! cycle it happened (root's own `diag_failed` above is otherwise
         ! invisible to the others), or a non-root rank iterates on into
         ! the next collective while root has already exited below --
         ! the same hang `ctl`'s own broadcast, right above, exists to
         ! avoid for the energy/commutator-norm pair.
         if (present(comm)) then
            if (comm%size() > 1) then
               diag_flag = [merge(1.0_dp, 0.0_dp, diag_failed), real(diag_info, dp)]
               call scf_context_bcast_flat(comm, diag_flag, 2)
               diag_failed = diag_flag(1) > 0.5_dp
               diag_info = nint(diag_flag(2))
            end if
         end if
         if (diag_failed) then
            ! A malformed Fock matrix (NOT the ordinary "has not converged
            ! yet" case -- dsyev/cuSOLVER could not even diagonalise it):
            ! stop iterating rather than build the next Fock from
            ! eps_s/cmo_s/d_s that linalg_provider%syev never actually
            ! wrote. scf_result%converged stays .false.; trc_capi.F90's
            ! trc_run_scf already turns that into TRC_ERR_NOCONV, and
            ! scf_result%iterations > 0 keeps it off the "setup failed"
            ! (TRC_ERR_UNSUPPORTED) path -- the SAME two fields, and the
            ! SAME mapping, an ordinary non-converged run already uses.
            scf_result%message = "trc_scf: the Fock-matrix diagonalisation did not converge "// &
                                 "(eigensolver info="//to_string(diag_info)//")"
            scf_result%iterations = it
            exit
         end if
         ! Every rank must hold the same density to the bit. The Fock sum is
         ! broadcast already; the XC potential is accumulated with atomics
         ! whose order differs per rank, 1e-13 in the energy on water, and a
         ! difference that small compounds through a diverging iteration.
         if (present(comm)) then
            if (comm%size() > 1) then
               !$acc update self(scf_result%dmat)
               call scf_context_bcast_flat(comm, scf_result%dmat, n2*scf_config%nspin)
               !$acc update device(scf_result%dmat)
            end if
         end if
         call scf_state%update_errv(nao, scf_config%nspin, scf_result%dmat)
         drms = sqrt(linalg_provider%dot_3d(scf_state%errv, scf_state%errv))/real(nao, dp)
         ! The linear algebra between the XC build and the next Fock build:
         ! energy, Fock assembly, commutator, DIIS, eigensolve, density and
         ! its broadcast (PCM excluded, it has its own timer).
         scf_timers%tr_it = scf_timers%rest%get_elapsed_time()
         scf_timers%t_rest = scf_timers%t_rest + scf_timers%tr_it
         if (scf_config%talk) then
            if (scf_config%dft) then
               write (line, "(i5,f22.12,es14.4,es14.4,es14.4,3f10.3,a)") &
                  it, scf_energies%total, scf_energies%total - eold, drms, errmax, scf_timers%tf_it, scf_timers%tx_it, &
                  scf_timers%tr_it, merge(" i", "  ", incr_now)
            else
               write (line, "(i5,f22.12,es14.4,es14.4,es14.4,2f10.3,a)") &
                  it, scf_energies%total, scf_energies%total - eold, drms, errmax, scf_timers%tf_it, &
                  scf_timers%tr_it, merge(" i", "  ", incr_now)
            end if
            call global_logger%info(trim(line))
         end if
         scf_result%iterations = it
         ! A real density, on the way past, for the benchmarks. The model
         ! density bench_gc invents is small and banded, so it screens far
         ! harder than anything an SCF actually meets, and a Fock build
         ! timed on it flatters the code by about three times. Set
         ! TRC_DUMP_DENSITY to a path and the density leaving iteration 2 is
         ! written there, in the stream format read_density already reads.
         if (it == 2) then
            ! The density lives on the DEVICE -- every update to it since the
            ! guess was a gemm up there -- so the host copy is the SAD guess
            ! until it is pulled back. Dumping without this wrote a
            ! block-diagonal atomic density and called it iteration 2, which
            ! screens about three times harder than anything real and made
            ! every benchmark built on the file flattering and wrong.
            !$acc update self(scf_result%dmat)
            call dump_density(nao, scf_result%dmat(:, :, 1))
         end if
         if (it > 1 .and. abs(scf_energies%total - eold) < scf_options%conv_energy .and. errmax < scf_options%conv_diis) then
            scf_result%converged = .true.
            exit
         end if
         eold = scf_energies%total
      end do

      ! The density the energy was evaluated at, and the orbitals that made it.
      !$acc update self(scf_state%dold, scf_result%cmo, scf_result%eps)
      if (present(comm)) then
         if (comm%size() > 1) then
            call scf_context_bcast_flat(comm, scf_result%cmo, n2*scf_config%nspin)
            call scf_context_bcast_flat(comm, scf_result%eps, nao*scf_config%nspin)
         end if
      end if
      call scf_state%release_device(scf_context, scf_result)
      call cdiis%release_device()
      if (scf_options%adiis) call adiis%release_device()
      scf_result%dmat = scf_state%dold
      scf_energies%pcm_qtotal = scf_context%pcm%q_total
      call scf_result%update(scf_energies, scf_state%dold)
      ! Only plain non-convergence gets this; `diag_failed` set a sharper message.
      if (.not. scf_result%converged .and. len_trim(scf_result%message) == 0) then
         scf_result%message = "trc_scf: not converged"
      end if
      call scf_driver_gap_check(scf_config, scf_options%frac_occ .or. scf_result%iterations == 0, nao, &
                                nao - scf_context%n_lindep, scf_result%eps)
      if (scf_config%dft) then
         ! The run's own totals, read once here: grid build was added after
         ! setup, above. xc is all of scf_fock_build_vxc but VV10 (t_xc includes it),
         ! not trc_xc's points/pairs timers, which miss the functional and,
         ! over ranks, the collectives on V (3.3 of 10.3 s on 16 ranks), so
         ! the rows add up to the total. Points/pairs are in the summary.
         call scf_timers%xcst%add_time(2, scf_timers%t_xc - scf_timers%t_vv10)
         call scf_timers%xcst%add_time(3, scf_timers%t_vv10)
      end if
      if (scf_config%talk) then
         call report_timing()
         if (scf_config%dft) then
            if (present(comm)) then
               call scf_timers%xcst%report(scf_context%t_grid + scf_timers%t_xc, comm%rank(), comm%size())
            else
               call scf_timers%xcst%report(scf_context%t_grid + scf_timers%t_xc, 0, 1)
            end if
         end if
         if (scf_config%pcm_on) then
            write (line, "(a,f20.12,a,f10.6,a,f10.3,a,f10.4,a)") "   pcm: E_diel = ", scf_energies%pcm, &
               "  q_total = ", scf_context%pcm%q_total, "  setup ", scf_context%t_pcm_setup, " s, ", &
               scf_timers%t_pcm/real(max(1, scf_result%iterations), dp), " s per iteration"
            call global_logger%info(trim(line))
         end if
         if (scf_config%df) then
            if (present(comm)) then
               call scf_context%dfj%report(scf_timers%t_fock, comm%rank(), comm%size())
            else
               call scf_context%dfj%report(scf_timers%t_fock, 0, 1)
            end if
         end if
      end if
      call scf_context%release()

   contains

      subroutine report_timing()
         !! Where the run's time went, through pic's logger at PERFORMANCE,
         !! which its default INFO level shows: the phases, then setup and
         !! the ERI engine's build broken into their steps.
         character(len=160) :: line
         ! `total` is the run's own wall time, so it also holds what no column
         ! times (PCM, the guess density's upload, the final transfers).
         write (line, "(a,8f10.3)") "          ", scf_context%t_setup, scf_context%t_grid, scf_timers%t_fock, scf_timers%t_xc, &
            scf_timers%tx_pts, scf_timers%tx_prs, scf_timers%t_rest, scf_timers%run%get_elapsed_time()
         call global_logger%performance( &
            "   time (s)   setup      grid      fock        xc   (points     pairs)      linalg     total")
         call global_logger%performance(trim(line))
         call global_logger%performance( &
            "   setup (s)  pairs        1e   sph fold       eri    linalg    orthog")
         write (line, "(a,6f10.3)") "          ", scf_context%ts
         call global_logger%performance(trim(line))
         call global_logger%performance( &
            "   eri (s) hgp pairs   schwarz      bins  quartets    upload")
         write (line, "(a,5f10.3)") "          ", scf_context%pec%eri%t_build
         call global_logger%performance(trim(line))
         if (scf_config%dft) then
            call global_logger%performance( &
               "   grid (s)   points  partition  batching    upload")
            write (line, "(a,4f10.3)") "          ", scf_context%t_grid_sub
            call global_logger%performance(trim(line))
         end if
      end subroutine report_timing

      ! Guess for one spin from a Fock-like h, into `dout`: the eigenproblem
      ! and the products on the device, the result copied back for the
      ! resident copy made below. h and x are host arrays here, copied up for
      ! the region; the loop runs the same sequence on resident ones.
      subroutine guess_from(h, nocc_s, dout, info)
         real(dp), intent(in) :: h(nao, nao)
         integer, intent(in) :: nocc_s
         real(dp), intent(out) :: dout(nao, nao)
         integer, intent(out) :: info
         real(dp), allocatable :: fp(:, :), c(:, :), t(:, :), e(:), d(:, :)
         allocate (fp(nao, nao), c(nao, nao), t(nao, nao), e(nao), d(nao, nao))
         !$acc data copyin(h, scf_context%x) create(fp, c, t, e) copyout(d)
         call orbitals_and_density(h, nocc_s, .false., t, fp, e, c, d, scf_state%ofrac, info)
         !$acc end data
         dout = d
      end subroutine guess_from

      subroutine orbitals_and_density(f, nocc_s, frac, t, fp, eps_s, cmo_s, d_s, ofrac, info)
         !! One spin's step from a Fock-like matrix to its density: F' = X^T F X,
         !! its eigenvectors, C = X C', then D = occ C(:,1:nocc) C(:,1:nocc)^T
         !! (or, with `frac`, the fractional occupations of the free-atom runs).
         !! Everything is on the device already (the loop's resident arrays, or
         !! the data region of `guess_from`); `t` and `fp` are the caller's
         !! n x n workspaces, which hold nothing on return that matters. `ofrac`
         !! is an explicit argument (not host-associated from `scf_state%ofrac`
         !! directly) so the `frac` branch's `do concurrent` below never names
         !! a derived-type component.
         !!
         !! `info` comes back non-zero when F' failed to diagonalise (a
         !! malformed Fock matrix -- a density read past the end of an
         !! undersized buffer was one real cause, found and fixed on this
         !! branch): the caller's own `do s` loop stops using `cmo_s`/`d_s`
         !! for this spin and the main iteration reports non-convergence
         !! instead of continuing on a matrix linalg_provider%syev never
         !! actually filled in.
         real(dp), intent(in) :: f(nao, nao)
         integer, intent(in) :: nocc_s
         logical, intent(in) :: frac
         real(dp), intent(inout) :: t(nao, nao), fp(nao, nao)
         real(dp), intent(inout) :: eps_s(nao), cmo_s(nao, nao), d_s(nao, nao)
         real(dp), intent(inout) :: ofrac(nao)
         integer, intent(out) :: info
         integer :: i, j
         ! F' = X^T F X
         call linalg_provider%gemm_mat(f, scf_context%x, t)
         call linalg_provider%gemm_mat(scf_context%x, t, fp, transa="T")
         call scf_driver_null_shift(linalg_provider, scf_context, nao, fp)
         call linalg_provider%syev(fp, eps_s, info)
         if (info /= 0) return
         call linalg_provider%gemm_mat(scf_context%x, fp, cmo_s)
         if (frac) then
            ! D = sum_i o_i c_i c_i^T through a column-scaled copy of C.
            !$acc update self(eps_s)
            if (scf_config%nspin == 1) then
               call frac_occupations(nao, eps_s, scf_options%nelec_frac, 2.0_dp, ofrac, nfrac)
            else
               call frac_occupations(nao, eps_s, real(nocc_s, dp), 1.0_dp, ofrac, nfrac)
            end if
            !$acc update device(ofrac)
            do concurrent(i=1:nao, j=1:nao)
               t(i, j) = cmo_s(i, j)*sqrt(ofrac(j))
            end do
            ! D = occ * C(:,1:k) C(:,1:k)^T with occ = 1, k = nfrac, on the
            ! column-scaled copy -- the same GEMM linalg_provider%density does on t with
            ! those arguments (A(1:nao,1:nfrac) is the same block either way).
            call linalg_provider%density(t, nfrac, 1.0_dp, d_s)
         else
            ! k <= 0 is linalg_provider%density's own zeroing branch (do concurrent),
            ! replacing zero_dev.
            call linalg_provider%density(cmo_s, nocc_s, scf_config%occ, d_s)
         end if
      end subroutine orbitals_and_density

      pure subroutine frac_occupations(n, eps, nelec, omax, o, nocc_out)
         !! Aufbau over the eigenvalues, `omax` per orbital (2 restricted, 1 per
         !! spin channel), the frontier's degenerate set (within FRAC_TOL)
         !! sharing the remainder evenly. `nocc_out` is the last orbital with
         !! any occupation, for the GEMM.
         integer, intent(in) :: n
         real(dp), intent(in) :: nelec, omax
         real(dp), intent(in) :: eps(n)
         real(dp), intent(out) :: o(n)
         integer, intent(out) :: nocc_out
         real(dp), parameter :: FRAC_TOL = 1.0e-3_dp
         integer :: nfull, i, g0, g1, ninside
         real(dp) :: rem, ef, share
         nfull = int(nelec/omax)
         rem = nelec - omax*real(nfull, dp)
         o = 0.0_dp
         if (nfull > 0) o(1:min(nfull, n)) = omax
         nocc_out = min(nfull, n)
         if (nfull >= n) return
         ! the frontier: the first not-full orbital when electrons are left,
         ! else the last full one (an open shell of even count, carbon 2p^2)
         if (rem > 0.0_dp .or. nfull == 0) then
            ef = eps(nfull + 1)
         else
            ef = eps(nfull)
         end if
         g0 = n + 1
         g1 = 0
         do i = 1, n
            if (abs(eps(i) - ef) < FRAC_TOL) then
               g0 = min(g0, i)
               g1 = max(g1, i)
            end if
         end do
         if (g1 < g0) return
         ninside = count([(i <= nfull, i=g0, g1)])
         share = (omax*real(ninside, dp) + rem)/real(g1 - g0 + 1, dp)
         o(g0:g1) = share
         nocc_out = g1
      end subroutine frac_occupations

   end subroutine trc_scf_run

   subroutine scf_driver_null_shift(la, scf_context, n, fp)
      !! With functions dropped by canonical orthogonalisation, F' = X^T F X
      !! has a zero eigenvalue per dropped direction, a zero orbital that
      !! aufbau occupied below a positive HOMO. Add shift*(I - X^T S X) to
      !! lift them above every real level (2|F'|_F + 1 bounds the spectrum):
      !! never occupied, always last, at eps(n - n_lindep + 1:n). A no-op
      !! when nothing was dropped.
      type(trc_linalg_t), intent(in) :: la
      type(scf_ctx_t), intent(in) :: scf_context
      integer, intent(in) :: n
      real(dp), intent(inout) :: fp(n, n)
      if (scf_context%n_lindep == 0) return
      call null_shift_add(n, 2.0_dp*sqrt(la%dot_2d(fp, fp)) + 1.0_dp, scf_context%pnull, fp)
   end subroutine scf_driver_null_shift

   subroutine null_shift_add(n, shift, p, a)
      !! a = a + shift*p, on the device under the GPU build.
      integer, intent(in) :: n
      real(dp), intent(in) :: shift
      real(dp), intent(in) :: p(n, n)
      real(dp), intent(inout) :: a(n, n)
      integer :: i, j
      do concurrent(i=1:n, j=1:n)
         a(i, j) = a(i, j) + shift*p(i, j)
      end do
   end subroutine null_shift_add

   subroutine scf_driver_gap_check(scf_config, skip, nao, nmo, eps)
      !! A warning, never a refusal, when the final orbitals' HOMO-LUMO gap
      !! is tiny. A closed-shell molecule run at the wrong charge stalled for
      !! 60 iterations with |FDS - SDF| stuck near 1e-2 before anyone looked
      !! at the charge; a gap this small says to look there first. Rank 0
      !! only; `skip` for fractional occupations (no aufbau gap) or no
      !! iteration run.
      type(scf_settings_t), intent(in) :: scf_config
      logical, intent(in) :: skip
      integer, intent(in) :: nao, nmo
      real(dp), intent(in) :: eps(nao, scf_config%nspin)
      real(dp), parameter :: GAP_WARN = 1.0e-2_dp
      !! 0.27 eV. Closed-shell organics sit well above it even with a pure
      !! functional (singlet O2 at PBE/6-31G*, a degenerate pi* pair half
      !! filled, is 0.039); only a genuinely degenerate or wrongly charged
      !! system gets this close.
      real(dp) :: gap
      character(len=160) :: line
      if (skip .or. .not. scf_config%root) return
      gap = scf_settings_homo_lumo_gap(nao, scf_config%nspin, nmo, scf_config%nocc, eps)
      if (gap < GAP_WARN) then
         write (line, "(a,es10.3,a)") "trc_scf: HOMO-LUMO gap is only ", gap, &
            " Eh; check the charge and multiplicity"
         call global_logger%warning(trim(line))
      end if
   end subroutine scf_driver_gap_check

end module trc_scf_driver
