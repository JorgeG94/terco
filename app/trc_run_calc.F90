submodule (trc_run) trc_run_calc
   !! The calculation of the `terco` command: device binding, the handle, the
   !! SCF and everything after it (bodies of `trc_run`).
   !!
   !! None of these procedures has an internal (`contains`) procedure, and
   !! `run_calculate_setup`/`run_calculate_scf` make the collective MPI calls: see the
   !! header comment of `trc_run` for why that is load-bearing.
   use, intrinsic :: iso_c_binding, only: c_null_ptr
   use pic_types, only: dp
   use pic_knowledge, only: get_knowledge
   use pic_logger, only: global_logger, verbose_level
   use pic_strings, only: to_string
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                       trc_set_aux_json, trc_set_method, trc_set_nlc_grid, trc_set_convergence, trc_set_scf_accel, &
                       trc_set_soscf, &
                       trc_set_level_shift, trc_set_lindep, trc_set_screening, &
                       trc_set_incremental, trc_set_df, trc_set_jkfit_json, &
                       trc_set_guess, trc_set_verbose, trc_set_rimp2, trc_set_force_cartesian, trc_set_comm, &
                       trc_run_scf, trc_run_rimp2, trc_nao, trc_is_spherical, trc_nspin, trc_energy, &
                       trc_energy_parts, trc_converged, trc_iterations, trc_rimp2_energy, &
                       trc_set_pcm, trc_pcm_energy, trc_save_scf, trc_load_scf, &
                       trc_run_gradient, trc_gradient, trc_set_geometry, trc_set_hessian_step, trc_run_hessian, &
                       trc_hessian, trc_hirshfeld_charges, trc_fukui_hirshfeld, &
                       TRC_OK, TRC_ERR_NOCONV
   use trc_api, only: trc_bind_device
   use trc_vib, only: trc_vibrational_analysis
   use trc_thermo, only: trc_thermochemistry
   use trc_constants, only: ATM_TO_PA
   use trc_error, only: hess_err_t => error_t
   use trc_app, only: app_c_string
   use trc_dlfind_bridge, only: trc_opt_options_t, trc_opt_result_t, trc_dlfind_optimize, &
                                TRC_OPT_COORDS_HDLC, TRC_OPT_COORDS_HDLC_TC
   implicit none

contains

   module subroutine run_calculation(r)
      !! The calculation: handle setup, SCF, post-SCF steps, reports and the exit code.
      type(run_t), intent(inout) :: r
   ! ---- the calculation -------------------------------------------------
   call run_calculate_setup(r)
   call run_calculate_scf(r)
   call run_calculate_post(r)
   call run_calculate_report(r)
   call run_calculate_finish(r)
   end subroutine run_calculation

   module subroutine run_calculate_setup(r)
      !! Bind the device, create the handle and hand it the molecule and every setting.
      type(run_t), intent(inout) :: r
      integer :: dev, aux_block
      integer(c_int) :: st
   ! One process, one device: the first, or (more than one rank) this
   ! rank's own, by NODE-LOCAL rank -- the same split trc_build_comm does
   ! inside the library, done here too because device binding must happen
   ! before anything reaches the device, and the guess the library builds
   ! starts doing that before it gets to call trc_build_comm itself. A host
   ! build has nothing to bind either way. When trc_bind_device_env above
   ! already bound this device before MPI_Init, this picks the same one and
   ! changes nothing.
   if (r%nranks > 1) then
      block
         type(comm_t) :: node_comm
         node_comm = r%wcomm%split()
         dev = trc_bind_device(node_comm%rank())
         call node_comm%finalize()
      end block
   else
      dev = trc_bind_device(0)
   end if
   st = trc_create(r%h)
   ! Every rank must make this call, with the same handle, before trc_run_scf
   ! -- trc_set_comm's own contract. One rank keeps today's path exactly:
   ! the context's fcomm defaults to -1 and is never touched.
   if (r%nranks > 1) call run_api(r, trc_set_comm(r%h, int(r%wcomm%fortran_handle(), c_int)))
   ! The SCF driver's own per-iteration table, basis/lindep/pcm lines and
   ! timing report (trc_scf_settings%verbose, "talk" in src/scf) are all
   ! gated by this one flag. `base_verbose` is today's rule (on unless
   ! --silent) and applies to every SCF a plain caller would recognise as
   ! "the" SCF/gradient/Hessian run -- a single point, --grad, --hessian/
   ! --freq's own reference SCF, DF, DFT, ... never quieted by --logger.
   ! `loop_verbose` is what --opt/--opt-neb/DL-FIND's own evaluations and
   ! the FD Hessian's displaced points are switched to instead, only for
   ! the duration of those two loops (below) -- on again at verbose/debug,
   ! same as mqc forcing its own logger down around a per-step QC call
   ! unless already at large_info+.
   r%base_verbose = merge(0_c_int, 1_c_int, r%cli%is_set("silent"))
   r%loop_verbose = merge(0_c_int, r%base_verbose, global_logger%log_level < verbose_level)
   call run_api(r, trc_set_verbose(r%h, r%base_verbose))
   call run_api(r, trc_set_molecule(r%h, int(r%natm, c_int), real(r%z, c_double), r%xyz, int(r%chrg, c_int), &
                                int(r%mult, c_int)))
   call run_api(r, trc_set_force_cartesian(r%h, merge(1_c_int, 0_c_int, r%cart)))
   call run_api(r, trc_set_basis_json(r%h, app_c_string(r%basis_path)))
   if (r%rimp2) call run_api(r, trc_set_aux_json(r%h, app_c_string(r%aux_path)))
   if (r%df) then
      call run_api(r, trc_set_jkfit_json(r%h, app_c_string(r%jkfit_path)))
      call run_api(r, trc_set_df(r%h, 1_c_int))
   end if
   call run_api(r, trc_set_method(r%h, app_c_string(r%func), int(r%grid, c_int)))
   call run_api(r, trc_set_nlc_grid(r%h, int(r%nlc_grid, c_int)))
   if (len(r%pcm) > 0) call run_api(r, trc_set_pcm(r%h, app_c_string(r%pcm), real(r%eps, c_double), 0_c_int, 0.0_c_double))
   call run_api(r, trc_set_convergence(r%h, real(r%econv, c_double), real(r%gconv, c_double), int(r%maxiter, c_int), &
                                   int(r%ndiis, c_int)))
   ! The convergence aids, only when asked for: otherwise the library's own
   ! defaults (by-spin damping, no ADIIS, a 1 Eh shift before DIIS) stand.
   if (r%cli%is_set("level-shift")) call run_api(r, trc_set_level_shift(r%h, real(r%level_shift, c_double)))
   if (r%accel_set) call run_api(r, trc_set_scf_accel(r%h, merge(1_c_int, 0_c_int, r%adiis), &
                                                  int(r%ndamp, c_int), real(r%damp, c_double), &
                                                  real(r%diis_start, c_double)))
   if (r%soscf) call run_api(r, trc_set_soscf(r%h, 1_c_int, real(r%soscf_start, c_double)))
   call run_api(r, trc_set_screening(r%h, real(r%thresh, c_double)))
   if (r%lindep >= 0.0_dp) call run_api(r, trc_set_lindep(r%h, real(r%lindep, c_double)))
   ! On unless asked otherwise: the library's default is off, so a host
   ! code keeps the build it had, but a command-line run gets the speed.
   call run_api(r, trc_set_incremental(r%h, merge(1_c_int, 0_c_int, r%cli%is_set("incremental"))))
   call run_api(r, trc_set_guess(r%h, int(r%guess_kind, c_int), c_null_ptr, 1_c_int))
   if (r%rimp2) then
      ! The three-index tensor (mu nu|P) is nao^2 per auxiliary function and
      ! is held one block of them at a time; the library's own default is
      ! one block of all of them, which at a thousand functions and a
      ! cc-pVDZ-RIFIT basis is 30 GB. 2 GB a block leaves the device for the
      ! transformed tensor, which is the one RI-MP2 cannot split.
      if (r%cli%is_set("aux-block")) then
         call r%cli%get("aux-block", aux_block, r%err)
         call run_check(r, "--aux-block")
         if (aux_block < 1) call run_die(r, "--aux-block must be at least 1")
      else
         call run_api(r, trc_nao(r%h, r%nao))
         aux_block = int(max(16.0_dp, min(real(huge(1), dp), 2.0e9_dp/(8.0_dp*real(r%nao, dp)**2))))
      end if
      call run_api(r, trc_set_rimp2(r%h, int(r%nfrozen, c_int), int(aux_block, c_int)))
   end if
   end subroutine run_calculate_setup

   module subroutine run_calculate_scf(r)
      !! The SCF itself: `--opt`, a `--restart` load, or a plain run; then the energies and RI-MP2.
      type(run_t), intent(inout) :: r
      integer(c_int) :: st, nspin
   if (r%opt) then
      call run_calculate_opt(r, st)
   else if (len(r%restart) > 0) then
      st = trc_load_scf(r%h, app_c_string(r%restart))
   else
      st = trc_run_scf(r%h)
   end if
   if (st /= TRC_OK .and. st /= TRC_ERR_NOCONV) call run_api(r, st)
   ! Rank 0 alone: every rank holds the same result.
   if (len(r%save_scf) > 0 .and. r%rank0) call run_api(r, trc_save_scf(r%h, app_c_string(r%save_scf)))
   call run_api(r, trc_nao(r%h, r%nao))
   call run_api(r, trc_is_spherical(r%h, r%sph))
   call run_api(r, trc_nspin(r%h, nspin))
   call run_api(r, trc_energy(r%h, r%etot))
   call run_api(r, trc_energy_parts(r%h, r%enuc, r%eone, r%etwo, r%exc))
   call run_api(r, trc_converged(r%h, r%conv))
   call run_api(r, trc_iterations(r%h, r%niter))
   call run_api(r, trc_pcm_energy(r%h, r%epcm, r%qpcm, r%npcm))
   r%eos = 0
   r%ess = 0
   if (r%rimp2 .and. r%conv /= 0) then
      call run_api(r, trc_run_rimp2(r%h))
      call run_api(r, trc_rimp2_energy(r%h, r%eos, r%ess))
      ! The double hybrid's own weighted total comes from trc_energy
      ! itself (now that RI-MP2 has run, cx%rimp2_ok is set and its
      ! pt2_cos/pt2_css branch applies), not recomputed here from
      ! r%pt2_cos/r%pt2_css -- those stay around only to report the
      ! fractions themselves (write_json), never to duplicate the
      ! weighting trc_energy already does from the SAME table.
      if (r%dh) call run_api(r, trc_energy(r%h, r%etot_dh))
   end if
   end subroutine run_calculate_scf

   module subroutine run_calculate_opt(r, st)
      !! The DL-FIND run (`--opt`) and the SCF re-converged at its final geometry.
      type(run_t), intent(inout) :: r
      integer(c_int), intent(out) :: st
      integer, allocatable :: opt_residues(:)
      integer :: opt_base_slash
      character(len=:), allocatable :: ext_lower
      type(trc_opt_options_t) :: opt_opts
      type(trc_opt_result_t) :: opt_res
      type(hess_err_t) :: herr
   allocate (opt_residues(r%natm))
   if (r%opt_coords_kind == TRC_OPT_COORDS_HDLC .or. r%opt_coords_kind == TRC_OPT_COORDS_HDLC_TC) then
      call run_build_opt_residues(r%natm, r%z, r%xyz, opt_residues, r%err)
      if (r%err%has_error()) call run_die(r, r%err%get_message())
   else
      opt_residues = 1   ! cartesian/dlc/dlc-tc: trc_dlfind_bridge's own driver needs no real partition
   end if
   ! --opt-constrain: a constraint that crosses an HDLC residue
   ! boundary, or touches a --opt-frozen atom, reaches DL-FIND's own
   ! assign_cons and dies there with "crosses residue boundary" ->
   ! hdlc_errflag('stop') -> dlf_fail -> error stop -- caught here,
   ! against the SAME residue partition cb_get_params hands DL-FIND
   ! (opt_residues, computed just above, BEFORE the frozen overwrite
   ! below turns every frozen atom's own real residue number into -1
   ! and would make this check see every pair of frozen atoms as
   ! "different residues" even when they share one).
   if (allocated(r%opt_cons_kind)) call run_check_opt_constrain_residues(r, opt_residues)
   ! --opt-frozen: -1 in the SAME spec slot residues already occupy
   ! (app/trc_dlfind_bridge.F90's own cb_get_params copies spec(1:natm)
   ! straight from this array, no translation) -- validated against
   ! pure DLC/DLC-TC before this point (app/main.F90's own --opt-frozen
   ! block, above).
   if (allocated(r%opt_frozen_atoms)) opt_residues(r%opt_frozen_atoms) = -1
   r%opt_is_linear = run_opt_geometry_is_linear(r%natm, r%z, r%xyz, herr)
   if (herr%has_error()) call run_die(r, "--opt: "//herr%get_message())

   opt_opts%coords_kind = r%opt_coords_kind
   opt_opts%guess_kind = r%opt_guess_kind
   opt_opts%fallback_guess_kind = int(r%guess_kind, c_int)
   opt_opts%maxiter = r%opt_maxiter
   opt_opts%gtol = r%opt_gtol
   opt_opts%algorithm = r%opt_algorithm_kind
   opt_opts%target_kind = r%opt_target_kind
   opt_opts%saddle_method = r%opt_saddle_method_kind
   opt_opts%connect = r%opt_connect
   opt_opts%is_linear = r%opt_is_linear
   opt_opts%hessian_source = r%opt_hessian_source_kind
   opt_opts%hessian_update = r%opt_hessian_update_kind
   opt_opts%max_updates = r%opt_max_updates
   opt_opts%zero_modes = r%opt_zero_modes
   opt_opts%soft_modes = r%opt_soft_modes
   opt_opts%max_step = r%opt_max_step
   opt_opts%linesearch = r%opt_linesearch_kind
   opt_opts%scalestep = r%opt_scalestep
   opt_opts%lbfgs_memory = r%opt_lbfgs_memory
   opt_opts%etol = r%opt_etol
   opt_opts%print_level = r%opt_print_level

   ! --opt-neb
   opt_opts%neb = r%opt_neb
   if (r%opt_neb) then
      opt_opts%neb_images = r%opt_neb_images
      opt_opts%neb_spring = r%opt_neb_spring
      opt_opts%neb_climb = r%opt_neb_climb
      opt_opts%neb_endpoints = r%opt_neb_endpoints_kind
      opt_opts%neb_product = r%opt_neb_product_xyz
      if (allocated(r%opt_neb_via_xyz)) opt_opts%neb_via = r%opt_neb_via_xyz
   end if

   ! --opt-constrain: this module's own plain parallel arrays
   ! (run_t%opt_cons_kind's own comment for why) converted to
   ! trc_dlfind_bridge's trc_opt_constraint_t here, the one place in
   ! this executable that imports both.
   if (allocated(r%opt_cons_kind)) then
      block
         integer :: ic
         allocate (opt_opts%constraints(size(r%opt_cons_kind)))
         do ic = 1, size(r%opt_cons_kind)
            opt_opts%constraints(ic)%kind = r%opt_cons_kind(ic)
            opt_opts%constraints(ic)%a = r%opt_cons_a(ic)
            opt_opts%constraints(ic)%b = r%opt_cons_b(ic)
            opt_opts%constraints(ic)%c = r%opt_cons_c(ic)
            opt_opts%constraints(ic)%d = r%opt_cons_d(ic)
         end do
      end block
   end if

   ! --opt-algorithm stochastic/genetic (iopt 51/52): opt_opts%po_*
   ! stay at their sentinel defaults (DL-FIND's own inline values,
   ! trc_dlfind_bridge's own header) unless run_parse_opt_global_options
   ! validated one above, which it refuses together with any other
   ! --opt-algorithm, so these are always still at their r%-type
   ! sentinel (-1/-1.0/huge(0)) when opt_is_global is false.
   opt_opts%seed = r%opt_seed
   opt_opts%po_pop_size = r%opt_global_pop_size
   opt_opts%po_maxcycle = r%opt_global_maxcycle
   opt_opts%po_init_pop_size = r%opt_global_init_pop_size
   opt_opts%po_nsave = r%opt_global_nsave
   opt_opts%po_distribution = r%opt_global_distribution
   opt_opts%po_reset = r%opt_global_reset
   opt_opts%po_radius = r%opt_global_radius
   opt_opts%po_contraction = r%opt_global_contraction
   opt_opts%po_tolerance_r = r%opt_global_tolerance_r
   opt_opts%po_tolerance_g = r%opt_global_tolerance_g
   opt_opts%po_mutation_rate = r%opt_global_mutation_rate
   opt_opts%po_death_rate = r%opt_global_death_rate
   opt_opts%po_scalefac = r%opt_global_scalefac

   ! <geometry base>_opt.xyz/_traj.xyz's own path, and --json's (if
   ! given): trc_dlfind_optimize needs these BEFORE it runs, not after,
   ! because an unrecoverable mid-run SCF failure writes them from
   ! inside a DL-FIND callback (trc_dlfind_bridge%fail_eval) -- DL-FIND
   ! itself stops the process right after, so nothing after this call
   ! returns would ever get the chance.
   r%opt_base = r%geom
   if (len(r%opt_base) > 4) then
      call run_to_lower_into(r%opt_base(len(r%opt_base) - 3:), ext_lower)
      if (ext_lower == ".xyz") r%opt_base = r%opt_base(1:len(r%opt_base) - 4)
   end if
   ! Written into the CURRENT WORKING DIRECTORY, not next to the input
   ! geometry: a directory part in r%geom (an absolute path, or one
   ! outside the CWD) must not steer these writes somewhere a test's
   ! own WORKING_DIRECTORY (ctest's convention, test/CMakeLists.txt)
   ! does not expect -- DL-FIND's own path*.xyz already go to the CWD
   ! (dl-find.f90), so this matches.
   opt_base_slash = scan(r%opt_base, "/", back=.true.)
   if (opt_base_slash > 0) r%opt_base = r%opt_base(opt_base_slash + 1:)
   opt_opts%opt_xyz_path = r%opt_base//"_opt.xyz"
   opt_opts%traj_xyz_path = r%opt_base//"_traj.xyz"
   opt_opts%json_path = r%json
   opt_opts%geometry_path = r%geom
   opt_opts%coords_name = r%opt_coords_arg
   opt_opts%guess_name = r%opt_guess_arg

   ! Quiet for the whole DL-FIND run (every evaluation inside it,
   ! --opt-neb included) and the reconvergence SCF right below, which
   ! is the SAME geometry/density DL-FIND's own last accepted step
   ! left behind -- not a fresh single point a plain caller asked for.
   ! Restored to base_verbose before `die` on a failed run, and before
   ! the --grad/--hessian/RI-MP2 reporting below, which DO behave like
   ! an ordinary caller's own run.
   call run_api(r, trc_set_verbose(r%h, r%loop_verbose))
   call trc_dlfind_optimize(r%h, r%natm, r%z, opt_residues, r%xyz, opt_opts, opt_res)
   if (opt_res%failed) then
      call run_api(r, trc_set_verbose(r%h, r%base_verbose))
      call run_die(r, "--opt: "//opt_res%message)
   end if

   r%opt_converged = opt_res%converged
   r%opt_steps = opt_res%steps
   r%opt_evals = opt_res%energy_evaluations
   r%opt_hess_evals = opt_res%hessian_evaluations
   r%opt_gmax = opt_res%final_gmax
   r%opt_grms = opt_res%final_grms
   r%opt_n_traj = opt_res%n_traj
   if (allocated(opt_res%traj_energy)) r%opt_traj_energy = opt_res%traj_energy
   if (allocated(opt_res%traj_coords)) r%opt_traj_xyz = opt_res%traj_coords
   r%opt_have_ts = opt_res%have_ts
   if (opt_res%have_ts) then
      r%opt_ts_energy = opt_res%ts_energy
      r%opt_ts_xyz = opt_res%ts_coords
      r%opt_ts_gmax = opt_res%ts_gmax
      r%opt_ts_converged = opt_res%ts_converged
   end if
   r%opt_have_min_plus = opt_res%have_minimum_plus
   if (opt_res%have_minimum_plus) then
      r%opt_min_plus_energy = opt_res%minimum_plus_energy
      r%opt_min_plus_xyz = opt_res%minimum_plus_coords
      r%opt_min_plus_gmax = opt_res%minimum_plus_gmax
      r%opt_min_plus_converged = opt_res%minimum_plus_converged
   end if
   r%opt_have_min_minus = opt_res%have_minimum_minus
   if (opt_res%have_minimum_minus) then
      r%opt_min_minus_energy = opt_res%minimum_minus_energy
      r%opt_min_minus_xyz = opt_res%minimum_minus_coords
      r%opt_min_minus_gmax = opt_res%minimum_minus_gmax
      r%opt_min_minus_converged = opt_res%minimum_minus_converged
   end if
   r%opt_neb_have = opt_res%have_neb
   if (opt_res%have_neb) then
      r%opt_neb_n_images = opt_res%neb_n_images
      r%opt_neb_image_energy = opt_res%neb_image_energy
      r%opt_neb_image_xyz = opt_res%neb_image_coords
      r%opt_neb_climb_image = opt_res%neb_climb_image
      r%opt_neb_scf_iters = opt_res%neb_scf_iterations_per_image
   end if
   r%xyz = opt_res%final_coords
   r%opt_final_xyz = r%xyz
   r%grad = .true.   ! so the existing --grad path below reports the final gradient too

   ! Re-converge the context at the reported final (ACCEPTED) geometry:
   ! cb_get_gradient calls trc_set_geometry/trc_run_scf for EVERY
   ! evaluation, including a rejected trial that happens after the last
   ! accepted step (e.g. --opt-maxiter stops the run right after one),
   ! so by the time api_dl_find returns, the context may be sitting at
   ! that later trial's geometry/SCF, not opt_res%final_coords's own.
   ! Every getter below must report exactly the geometry --opt writes
   ! out.
   call run_api(r, trc_set_geometry(r%h, reshape(real(r%xyz, c_double), [3*r%natm])))
   st = trc_run_scf(r%h)
   ! Back to today's rule for everything from here on: --grad above
   ! forced true already reports the final gradient at base_verbose
   ! (an ordinary trc_run_gradient call, never run inside DL-FIND's own
   ! loop), and so does --hessian/--opt-freq's own FD Hessian below,
   ! which sets loop_verbose again itself, only around trc_run_hessian.
   call run_api(r, trc_set_verbose(r%h, r%base_verbose))
   end subroutine run_calculate_opt

   module subroutine run_calculate_post(r)
      !! The post-SCF steps: gradient, Hessian and charge analyses.
      type(run_t), intent(inout) :: r
   call run_calculate_gradient(r)
   call run_calculate_hessian(r)
   call run_calculate_charges(r)
   end subroutine run_calculate_post

   module subroutine run_calculate_gradient(r)
      !! The analytic nuclear gradient (`--grad`, and always after `--opt`).
      type(run_t), intent(inout) :: r
      real(c_double), allocatable :: gradbuf(:)
   if (r%grad .and. r%conv /= 0) then
      ! trc_run_gradient's own refusals (RI-MP2, PCM) come
      ! back through `api`, which prints the library's own message and
      ! stops -- the same clear-refusal path every other unsupported
      ! combination in this CLI already takes. An LDA/GGA --func (global
      ! and range-separated hybrids and VV10 included) is NOT refused:
      ! trc_run_gradient builds the KS-DFT gradient at the SCF's own
      ! functional and grids, and with --df the fitted one.
      call run_api(r, trc_run_gradient(r%h))
      allocate (gradbuf(3*r%natm), r%gradient(3, r%natm))
      call run_api(r, trc_gradient(r%h, gradbuf))
      r%gradient = reshape(real(gradbuf, dp), [3, r%natm])
      deallocate (gradbuf)
   end if
   end subroutine run_calculate_gradient

   module subroutine run_calculate_hessian(r)
      !! The finite-difference Hessian, the vibrational analysis and the thermochemistry.
      type(run_t), intent(inout) :: r
      real(c_double), allocatable :: hessbuf(:)
      real(dp), allocatable :: mass(:), displacement(:, :, :)
      type(hess_err_t) :: herr
      integer :: ndof
      integer(c_int) :: st
   if (r%hessian .and. r%conv /= 0) then
      ! trc_run_hessian's own refusals (meta-GGA, --rimp2, a double
      ! hybrid, --pcm, --df outside its own scope, more than one rank)
      ! come back through `api`, exactly as --grad's do -- checked
      ! before any displacement, so a refused run costs nothing beyond
      ! the SCF this CLI already ran for --grad/the energy anyway.
      call run_api(r, trc_set_hessian_step(r%h, real(r%hess_step, c_double)))
      ! Quiet for the 2*3*natm displaced SCFs only -- the reference SCF
      ! trc_run_hessian itself requires (this context's own, already run
      ! above at base_verbose, or at loop_verbose if --opt-freq's own
      ! reconvergence SCF is what left it converged) is never re-run here.
      ! Restored before `st` reaches `api`, which may die.
      call run_api(r, trc_set_verbose(r%h, r%loop_verbose))
      st = trc_run_hessian(r%h)
      call run_api(r, trc_set_verbose(r%h, r%base_verbose))
      call run_api(r, st)
      ndof = 3*r%natm
      allocate (hessbuf(ndof*ndof), r%hess(ndof, ndof))
      call run_api(r, trc_hessian(r%h, hessbuf))
      r%hess = reshape(real(hessbuf, dp), [ndof, ndof])
      deallocate (hessbuf)

      call run_build_masses(r%z, r%masses_kind, mass, herr)
      if (herr%has_error()) call run_die(r, "--masses: "//herr%get_message())
      call trc_vibrational_analysis(mass, r%xyz, r%hess, r%freq_cm1, r%redmass_amu, r%force_const, &
                                    displacement, r%is_linear, r%n_tr, herr)
      if (herr%has_error()) call run_die(r, "--hessian: vibrational analysis: "//herr%get_message())
      call trc_thermochemistry(mass, r%xyz, r%freq_cm1, real(r%etot, dp), r%sym_number, r%mult, r%thermo, herr, &
                               temperature_k=r%temperature_k, pressure_pa=r%pressure_atm*ATM_TO_PA)
      if (herr%has_error()) call run_die(r, "--hessian: thermochemistry: "//herr%get_message())
      r%hess_ok = .true.
   end if
   end subroutine run_calculate_hessian

   module subroutine run_calculate_charges(r)
      !! The Hirshfeld charges and Fukui indices.
      type(run_t), intent(inout) :: r
      real(c_double), allocatable :: chargebuf(:)
      real(c_double), allocatable :: fplusbuf(:), fminusbuf(:), fzerobuf(:), qnbuf(:), qnp1buf(:), qnm1buf(:)
   if (r%charges == "hirshfeld" .and. r%conv /= 0) then
      ! trc_hirshfeld_charges's own refusals (an element outside H-Kr) come
      ! back through `api`, same as --grad's/--hessian's.
      allocate (chargebuf(r%natm), r%hirshfeld_charges(r%natm))
      call run_api(r, trc_hirshfeld_charges(r%h, chargebuf))
      r%hirshfeld_charges = real(chargebuf, dp)
      deallocate (chargebuf)
   end if

   if (r%fukui == "hirshfeld" .and. r%conv /= 0) then
      ! trc_fukui_hirshfeld's own refusals (a non-singlet neutral, an
      ! unconverged ion, an element outside H-Kr) come back through
      ! `api`, same as --charges's/--grad's/--hessian's. It runs two more
      ! SCFs on r%h (N+1, N-1) and restores the handle's charge,
      ! multiplicity and converged state before returning -- its own
      ! header -- so nothing here needs to save or restore anything.
      allocate (fplusbuf(r%natm), fminusbuf(r%natm), fzerobuf(r%natm))
      allocate (qnbuf(r%natm), qnp1buf(r%natm), qnm1buf(r%natm))
      allocate (r%fukui_f_plus(r%natm), r%fukui_f_minus(r%natm), r%fukui_f_zero(r%natm))
      allocate (r%fukui_q_n(r%natm), r%fukui_q_np1(r%natm), r%fukui_q_nm1(r%natm))
      call run_api(r, trc_fukui_hirshfeld(r%h, fplusbuf, fminusbuf, fzerobuf, qnbuf, qnp1buf, qnm1buf))
      r%fukui_f_plus = real(fplusbuf, dp)
      r%fukui_f_minus = real(fminusbuf, dp)
      r%fukui_f_zero = real(fzerobuf, dp)
      r%fukui_q_n = real(qnbuf, dp)
      r%fukui_q_np1 = real(qnp1buf, dp)
      r%fukui_q_nm1 = real(qnm1buf, dp)
      deallocate (fplusbuf, fminusbuf, fzerobuf, qnbuf, qnp1buf, qnm1buf)
   end if
   end subroutine run_calculate_charges

   module subroutine run_calculate_report(r)
      !! The summaries, the JSON file and the `--opt` geometry files, on rank 0.
      type(run_t), intent(inout) :: r
   if (r%rank0) call run_summary(r)
   if (r%rank0 .and. r%grad .and. r%conv /= 0) call run_gradient_summary(r)
   if (r%rank0 .and. r%charges == "hirshfeld" .and. r%conv /= 0) call run_charges_summary(r)
   if (r%rank0 .and. r%fukui == "hirshfeld" .and. r%conv /= 0) call run_fukui_summary(r)
   if (r%rank0 .and. r%opt) call run_opt_summary(r)
   if (r%rank0 .and. r%hessian .and. r%hess_ok) call run_hessian_summary(r)
   if (r%rank0 .and. len(r%json) > 0) call run_write_json(r, r%json)
   if (r%rank0 .and. r%opt) then
      ! opt_base (and opt_opts%opt_xyz_path/traj_xyz_path built from it)
      ! was already computed before trc_dlfind_optimize ran, above.
      call run_write_opt_xyz(r, r%opt_base//"_opt.xyz", r%xyz)
      if (r%opt_n_traj > 0) call run_write_traj_xyz(r, r%opt_base//"_traj.xyz", r%opt_traj_xyz, r%opt_traj_energy)
      call global_logger%info("   wrote "//r%opt_base//"_opt.xyz and "//r%opt_base//"_traj.xyz")
      if (r%opt_connect) then
         if (r%opt_have_ts) then
            call run_write_labelled_xyz(r, r%opt_base//"_ts.xyz", r%opt_ts_xyz, r%opt_ts_energy, "transition state")
         end if
         if (r%opt_have_min_minus) then
            call run_write_labelled_xyz(r, r%opt_base//"_minus.xyz", r%opt_min_minus_xyz, r%opt_min_minus_energy, &
                                    "downhill minimum (-distort)")
         end if
         if (r%opt_have_min_plus) then
            call run_write_labelled_xyz(r, r%opt_base//"_plus.xyz", r%opt_min_plus_xyz, r%opt_min_plus_energy, &
                                    "downhill minimum (+distort)")
         end if
         call global_logger%info("   opt-connect: wrote "//r%opt_base//"_ts.xyz, "//r%opt_base//"_minus.xyz and "// &
                                 r%opt_base//"_plus.xyz")
      end if
      if (r%opt_neb_have) then
         call run_write_neb_path_xyz(r, r%opt_base//"_neb.xyz")
         call global_logger%info("   opt-neb: wrote "//r%opt_base//"_neb.xyz ("//to_string(r%opt_neb_n_images)// &
                                 " images)")
      end if
   end if
   end subroutine run_calculate_report

   module subroutine run_calculate_finish(r)
      !! Destroy the handle, turn non-convergence into the exit code, print the closing lines.
      type(run_t), intent(inout) :: r
      integer(c_int) :: st
   st = trc_destroy(r%h)

   if (r%conv == 0) then
      if (r%rank0) call global_logger%error("terco: the SCF did not converge in "//to_string(int(r%niter))//" iterations")
      call run_quit(r, 2)
   end if
   if (r%opt .and. .not. r%opt_converged) then
      ! Matches the SCF's own non-convergence above: a CLI run that did
      ! not reach the thing it was asked to do must not exit 0 and claim
      ! "normal termination" -- opt_summary already printed "NOT
      ! converged" above, this just makes the exit code agree with it.
      if (r%rank0) call global_logger%error("terco: --opt did not converge in "// &
                                             to_string(r%opt_steps)//" accepted step(s)")
      call run_quit(r, 2)
   end if
   if (r%rank0) then
      call global_logger%info("")
      call global_logger%info(" normal termination of terco")
      ! Re-read the clock: the Hessian and JSON output above come after the
      ! summary wall time, which would otherwise be buried.
      block
         character(len=1024) :: logline
         call r%wall%stop()
         write (logline, "(a,f12.3,a)") " total wall time    : ", r%wall%get_elapsed_time(), " s"
         call global_logger%info(trim(logline))
      end block
   end if
   ! A line of pic's collected lore, as metalquicha prints at the end of a
   ! successful run: rank 0 only, and not under --silent, whose promise is
   ! the summary and nothing else.
   if (r%rank0 .and. .not. r%cli%is_set("silent")) call get_knowledge()
   end subroutine run_calculate_finish

end submodule trc_run_calc
