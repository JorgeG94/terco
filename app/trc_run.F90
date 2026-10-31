module trc_run
   !! Everything `trc_main` (below) used to keep as internal (`contains`)
   !! procedures, moved out to genuine module procedures that take the
   !! run's state as an explicit argument (`run_t`) instead of reaching it
   !! by host association.
   !!
   !! This is not a style preference: an internal procedure ANYWHERE in this
   !! executable -- even one that is never called, and even one with no
   !! state in common with the SCF call at all -- corrupts the two-electron
   !! part of the Fock build the moment `trc_run_scf` runs collectively
   !! (more than one rank, this build's NVHPC 26.5 + this OpenMPI/HPCX):
   !! `conv` and `niter` still look fine, but the two-electron energy comes
   !! back as exactly 0 and the total is wrong, deterministically, on every
   !! rank alike. Measured by elimination: the program (or procedure) that
   !! makes the collective call must have NO sibling internal procedure at
   !! all -- not "none that touch the same state", none, full stop -- and a
   !! genuine module procedure (this file's old `contains` turned into
   !! `trc_run`'s) does not trip it. Reproduced down to an EMPTY, never
   !! invoked `contains subroutine noop(); end subroutine` in the same
   !! program; unaffected by optimisation level, by `-acc`/`-stdpar=gpu` on
   !! the caller, by MPI transport (ob1 vs UCX), or by UCX's memtype cache.
   !! Looks like an NVHPC/OpenMPI codegen interaction specific to this
   !! toolchain; worth a closer look (and a vendor report) on its own.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use pic_types, only: dp
   use pic_error, only: error_t
   use pic_cli, only: cli_t
   use pic_logger, only: global_logger
   use pic_timer, only: timer_type
   use pic_mpi_lib, only: comm_t, abort_comm, pic_mpi_finalize
   use trc_capi, only: trc_message, TRC_OK
   use trc_app, only: app_f_string
   use trc_thermo, only: trc_thermo_result_t
   use trc_error, only: mass_error_t => error_t
   implicit none
   private
   public :: run_t, VERSION, run_die, run_quit, run_check, run_api, run_int_option, run_banner, run_summary, run_write_json, &
             run_gradient_summary, run_hessian_summary, run_build_masses, run_charges_summary, &
             run_fukui_summary
   public :: run_opt_capability_check, run_build_opt_residues, run_opt_summary, run_write_opt_xyz, run_write_traj_xyz
   public :: run_parse_opt_global_options
   public :: run_opt_geometry_is_linear, run_write_labelled_xyz, run_to_lower_alloc, run_to_lower_into
   public :: run_parse_frozen_list, run_parse_constraints, run_check_opt_constrain_residues, &
             run_check_neb_endpoint, run_write_neb_path_xyz
   public :: run_read_xyz_frames
   public :: run_define_cli, run_parse, run_parse_inputs, run_parse_scf, &
             run_parse_features, run_parse_opt_flags, run_parse_opt_basic, run_parse_opt_target, &
             run_parse_opt_neb, run_parse_opt_hessian, run_parse_opt_steps, run_parse_opt_geometry, &
             run_parse_paths, run_calculation, run_calculate_setup, run_calculate_scf, &
             run_calculate_opt, run_calculate_post, run_calculate_gradient, run_calculate_hessian, &
             run_calculate_charges, run_calculate_report, run_calculate_finish

#ifndef TERCO_VERSION
#define TERCO_VERSION "unknown"
#endif
   character(len=*), parameter :: VERSION = TERCO_VERSION

   type :: run_t
      !! Every setting and result `trc_main`'s procedures need, passed
      !! explicitly -- see this module's own header comment for why.
      real(dp), allocatable :: xyz(:, :)   !! (3, natm) as app_read_xyz returns it; --opt replaces it with the final geometry
      integer :: maxiter = 0    !! --iterations, parsed in run_parse_scf, handed to the library in run_calculate_setup
      integer :: nfrozen = 0    !! --frozen, parsed in run_parse_scf, read by --rimp2 in run_calculate_setup
      integer :: guess_kind = 0   !! the TRC_GUESS_* value of `guess`; parsed early, read by the setup and --opt
      logical :: cart = .false.   !! --cart, parsed in run_parse_features, read in run_calculate_setup
      integer(c_int) :: base_verbose = 1, loop_verbose = 1
         !! trc_set_verbose's own argument: `base_verbose` is today's value
         !! (on unless --silent), restored around every SCF a plain caller
         !! would recognise as "the" SCF/gradient/Hessian run; `loop_verbose`
         !! is on only at --logger verbose/debug, set for exactly the SCFs
         !! evaluated INSIDE --opt/--opt-neb/DL-FIND's own loop and the FD
         !! Hessian's displaced points. Both are chosen in run_calculate_setup and
         !! used by run_calculate_opt and run_calculate_hessian.
      character(len=:), allocatable :: opt_base
         !! <geometry base name>_opt.xyz et al. without the suffix: built before
         !! DL-FIND runs (run_calculate_opt), used by the file writers in run_calculate_report.
      real(dp), allocatable :: opt_traj_xyz(:, :, :)   !! (3, natm, opt_n_traj), Bohr -- the --opt trajectory
      type(cli_t) :: cli
      type(error_t) :: err
      type(c_ptr) :: h
      character(len=:), allocatable :: geom, basis_arg, basis_path, aux_path, jkfit_path, func, guess, title, json
      character(len=:), allocatable :: pcm
      character(len=:), allocatable :: save_scf, restart   !! "" when not asked for
      integer, allocatable :: z(:)
      integer :: natm, chrg, uhf, nelec, mult, grid, nlc_grid
      integer(c_int) :: nao, conv, niter, sph, npcm
      real(dp) :: econv, gconv, thresh, eps
      real(dp) :: lindep = -1.0_dp  !! --lindep, negative = the library's default
      integer :: ndiis = 10              !! --diis
      real(dp) :: level_shift = 1.0_dp   !! --level-shift (only passed when given)
      logical :: adiis = .false.         !! --adiis
      logical :: accel_set = .false.     !! any of --adiis/--damp/--damp-steps/--diis-start given
      integer :: ndamp = -1              !! --damp-steps; -1 = by spin
      real(dp) :: damp = 0.3_dp          !! --damp
      real(dp) :: diis_start = -1.0_dp   !! --diis-start; -1 = by spin
      logical :: soscf = .false.         !! --soscf
      real(dp) :: soscf_start = -1.0_dp  !! --soscf-start; -1 = the library's default (1e-2)
      real(c_double) :: etot, enuc, eone, etwo, exc, eos, ess, epcm, qpcm
      logical :: rimp2, df
      logical :: dh = .false.   !! --func names a double hybrid: rimp2 runs even without --rimp2
      logical :: vv10 = .false.  !! --func has VV10 nonlocal correlation, on its own --nlc-grid
      real(c_double) :: pt2_cos = 0.0_c_double, pt2_css = 0.0_c_double, etot_dh = 0.0_c_double
      type(timer_type) :: wall
      type(comm_t) :: wcomm
      integer :: nranks
      logical :: rank0
      logical :: grad
      real(dp), allocatable :: gradient(:, :)   !! (3, natm), Eh/Bohr, after --grad
      character(len=:), allocatable :: charges  !! "" (none) or "hirshfeld", --charges
      real(dp), allocatable :: hirshfeld_charges(:)  !! (natm), e, after --charges hirshfeld
      character(len=:), allocatable :: fukui  !! "" (none) or "hirshfeld", --fukui
      ! (natm) each, after --fukui hirshfeld: f+/f-/f0 and the three
      ! states' own Hirshfeld charges (N, N+1, N-1) trc_fukui_hirshfeld
      ! condensed them from.
      real(dp), allocatable :: fukui_f_plus(:), fukui_f_minus(:), fukui_f_zero(:)
      real(dp), allocatable :: fukui_q_n(:), fukui_q_np1(:), fukui_q_nm1(:)
      ! --hessian: the FD Hessian (trc_run_hessian/trc_hessian) and the
      ! vibrational analysis/thermochemistry built on it (trc_vib/
      ! trc_thermo).
      logical :: hessian = .false.
      real(dp) :: hess_step = 0.005_dp, temperature_k = 298.15_dp, pressure_atm = 1.0_dp
      integer :: sym_number = 1
      character(len=:), allocatable :: masses_kind   !! "isotope" or "average"
      real(dp), allocatable :: hess(:, :)            !! (3*natm, 3*natm), Eh/Bohr^2
      real(dp), allocatable :: freq_cm1(:), redmass_amu(:), force_const(:)
      logical :: hess_ok = .false.
      logical :: is_linear = .false.
      integer :: n_tr = 0
      type(trc_thermo_result_t) :: thermo
      ! --opt: geometry optimisation through app/trc_dlfind_bridge.F90's
      ! trc_dlfind_optimize (TERCO_ENABLE_DLFIND=ON) or its stub (OFF, the
      ! default -- always refused with a message). opt_coords_kind/
      ! opt_guess_kind are trc_dlfind_bridge's own TRC_OPT_COORDS_*/
      ! TRC_OPT_GUESS_* enums, resolved from --opt-coords/--opt-guess here
      ! so main.F90 need not import that module's constants under every
      ! build configuration.
      logical :: opt = .false.
      character(len=:), allocatable :: opt_coords_arg, opt_guess_arg
      integer :: opt_coords_kind = 0, opt_guess_kind = 0
      integer :: opt_maxiter = 100
      real(dp) :: opt_gtol = 4.5e-4_dp
      real(dp) :: opt_max_step = -1.0_dp   !! negative = DL-FIND's own default (trc_opt_options_t%max_step)
      logical :: opt_converged = .false.
      integer :: opt_steps = 0, opt_evals = 0, opt_hess_evals = 0
      real(dp) :: opt_gmax = 0.0_dp, opt_grms = 0.0_dp
      integer :: opt_n_traj = 0
      real(dp), allocatable :: opt_traj_energy(:)
      real(dp), allocatable :: opt_final_xyz(:, :)   !! (3, natm), Bohr -- the optimised geometry
      ! Algorithm, target (minimum/saddle), the two saddle methods,
      ! `connect`, and the external-Hessian control -- every one of these
      ! is also one of
      ! trc_dlfind_bridge's own TRC_OPT_* enums (app/trc_dlfind_bridge.F90,
      ! app/trc_dlfind_stub.F90), stored here as the plain integer both
      ! modules define identically, same reason opt_coords_kind/
      ! opt_guess_kind above need no import of that module.
      character(len=:), allocatable :: opt_algorithm_arg, opt_target_arg, opt_saddle_method_arg, &
                                       opt_init_hessian_arg, opt_hessian_update_arg, opt_linesearch_arg
      integer :: opt_algorithm_kind = 3   !! TRC_OPT_ALG_LBFGS
      integer :: opt_target_kind = 0      !! TRC_OPT_TARGET_MINIMUM
      integer :: opt_saddle_method_kind = 0   !! TRC_OPT_SADDLE_PRFO
      logical :: opt_connect = .false.
      logical :: opt_is_linear = .false.
      integer :: opt_hessian_source_kind = 0  !! TRC_OPT_HESS_TERCO
      integer :: opt_hessian_update_kind = -1, opt_max_updates = -1, opt_zero_modes = -1
      real(dp) :: opt_soft_modes = -1.0_dp, opt_scalestep = -1.0_dp
      integer :: opt_linesearch_kind = -1   !! TRC_OPT_LINE_AUTO
      integer :: opt_lbfgs_memory = -1
      real(dp) :: opt_etol = -1.0_dp
      integer :: opt_print_level = -1
      ! --opt-connect's own three structures (the TS and both downhill
      ! minima); populated only when opt_connect is true and the run did
      ! not fail.
      logical :: opt_have_ts = .false., opt_have_min_plus = .false., opt_have_min_minus = .false.
      real(dp) :: opt_ts_energy = 0.0_dp, opt_min_plus_energy = 0.0_dp, opt_min_minus_energy = 0.0_dp
      real(dp) :: opt_ts_gmax = 0.0_dp, opt_min_plus_gmax = 0.0_dp, opt_min_minus_gmax = 0.0_dp
      logical :: opt_ts_converged = .false., opt_min_plus_converged = .false., opt_min_minus_converged = .false.
      real(dp), allocatable :: opt_ts_xyz(:, :), opt_min_plus_xyz(:, :), opt_min_minus_xyz(:, :)

      ! --- NEB / climbing-image NEB -------------------------------------------
      logical :: opt_neb = .false.
      character(len=:), allocatable :: opt_neb_endpoint_path, opt_neb_endpoints_arg
      integer :: opt_neb_images = -1
      real(dp) :: opt_neb_spring = -1.0_dp
      logical :: opt_neb_climb = .false.
      integer :: opt_neb_endpoints_kind = 1   !! TRC_OPT_NEB_PERP, trc_dlfind_bridge's own default
      real(dp), allocatable :: opt_neb_product_xyz(:, :)   !! (3, natm) Bohr -- the endpoint file, validated
      character(len=:), allocatable :: opt_neb_via_path
      real(dp), allocatable :: opt_neb_via_xyz(:, :, :)   !! (3, natm, n_via) Bohr -- --opt-neb-via, validated
      logical :: opt_neb_have = .false.
      integer :: opt_neb_n_images = 0
      real(dp), allocatable :: opt_neb_image_energy(:)
      real(dp), allocatable :: opt_neb_image_xyz(:, :, :)   !! (3, natm, opt_neb_n_images), Bohr
      integer :: opt_neb_climb_image = 0
      integer, allocatable :: opt_neb_scf_iters(:)

      ! --- frozen atoms --------------------------------------------------------
      character(len=:), allocatable :: opt_frozen_arg
      integer, allocatable :: opt_frozen_atoms(:)   !! 1-based, parsed from opt_frozen_arg

      ! --- internal-coordinate constraints -------------------------------------
      character(len=:), allocatable :: opt_constrain_arg
         !! ';'-separated ("bond i j"/"angle i j k"/"torsion i j k l"/
         !! "cartesian i"/"bond-difference i j k") -- pic_cli's own options
         !! are single-valued (the last occurrence wins), so repeats are
         !! spelled this way rather than by repeating the flag.
      integer, allocatable :: opt_cons_kind(:), opt_cons_a(:), opt_cons_b(:), opt_cons_c(:), opt_cons_d(:)
         !! Parsed parallel arrays: trc_dlfind_bridge's own TRC_OPT_CONS_*
         !! values by construction (1=bond,2=angle,3=torsion,4=cartesian,
         !! 5=bond-difference -- this module never imports that module
         !! under TERCO_ENABLE_DLFIND=OFF, same reason opt_coords_kind
         !! above is a plain integer).

      ! --- global optimisation (stochastic/genetic) ---------------------------
      logical :: opt_is_global = .false.   !! --opt-algorithm stochastic or genetic
      integer :: opt_seed = huge(0)        !! --opt-seed; huge(0) sentinel, trc_dlfind_bridge's own
      integer :: opt_global_pop_size = -1, opt_global_maxcycle = -1, opt_global_init_pop_size = -1, &
                opt_global_nsave = -1, opt_global_distribution = -1, opt_global_reset = -1
      real(dp) :: opt_global_radius = -1.0_dp, opt_global_contraction = -1.0_dp, opt_global_tolerance_r = -1.0_dp, &
                 opt_global_tolerance_g = -1.0_dp, opt_global_mutation_rate = -1.0_dp, &
                 opt_global_death_rate = -1.0_dp, opt_global_scalefac = -1.0_dp
         !! trc_dlfind_bridge's own trc_opt_options_t%po_*/seed, one for
         !! one -- parsed and range-checked in run_parse_opt_global_options,
         !! below, then copied across verbatim in the --opt dispatch
         !! (same "plain fields here, enums/values match the bridge's
         !! own module exactly" reason opt_coords_kind above gives).
   end type run_t

   interface
      module subroutine run_define_cli(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_define_cli

      module subroutine run_parse(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse

      module subroutine run_parse_inputs(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_inputs

      module subroutine run_parse_scf(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_scf

      module subroutine run_parse_features(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_features

      module subroutine run_parse_opt_flags(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_opt_flags

      module subroutine run_parse_opt_basic(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_opt_basic

      module subroutine run_parse_opt_target(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_opt_target

      module subroutine run_parse_opt_neb(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_opt_neb

      module subroutine run_parse_opt_hessian(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_opt_hessian

      module subroutine run_parse_opt_steps(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_opt_steps

      module subroutine run_parse_opt_geometry(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_opt_geometry

      module subroutine run_parse_paths(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_paths

      module subroutine run_calculation(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_calculation

      module subroutine run_calculate_setup(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_calculate_setup

      module subroutine run_calculate_scf(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_calculate_scf

      module subroutine run_calculate_opt(r, st)
         implicit none
         type(run_t), intent(inout) :: r
         integer(c_int), intent(out) :: st
      end subroutine run_calculate_opt

      module subroutine run_calculate_post(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_calculate_post

      module subroutine run_calculate_gradient(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_calculate_gradient

      module subroutine run_calculate_hessian(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_calculate_hessian

      module subroutine run_calculate_charges(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_calculate_charges

      module subroutine run_calculate_report(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_calculate_report

      module subroutine run_calculate_finish(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_calculate_finish

      module subroutine run_banner(r)
         implicit none
         type(run_t), intent(in) :: r
      end subroutine run_banner

      module subroutine run_summary(r)
         implicit none
         type(run_t), intent(in) :: r
      end subroutine run_summary

      module subroutine run_gradient_summary(r)
         implicit none
         type(run_t), intent(in) :: r
      end subroutine run_gradient_summary

      module subroutine run_charges_summary(r)
         implicit none
         type(run_t), intent(in) :: r
      end subroutine run_charges_summary

      module subroutine run_fukui_summary(r)
         implicit none
         type(run_t), intent(in) :: r
      end subroutine run_fukui_summary

      module subroutine run_hessian_summary(r)
         implicit none
         type(run_t), intent(in) :: r
      end subroutine run_hessian_summary

      module subroutine run_build_masses(z, kind, mass, err)
         implicit none
         integer, intent(in) :: z(:)
         character(len=*), intent(in) :: kind
         real(dp), allocatable, intent(out) :: mass(:)
         type(mass_error_t), intent(inout) :: err
      end subroutine run_build_masses

      module function run_opt_capability_check(r) result(refusal)
         implicit none
         type(run_t), intent(in) :: r
         character(len=:), allocatable :: refusal
      end function run_opt_capability_check

      module function run_opt_geometry_is_linear(natm, z, xyz, err) result(linear)
         implicit none
         integer, intent(in) :: natm
         integer, intent(in) :: z(natm)
         real(dp), intent(in) :: xyz(3, natm)
         type(mass_error_t), intent(inout) :: err
         logical :: linear
      end function run_opt_geometry_is_linear

      module subroutine run_build_opt_residues(natm, z, xyz, residues, err)
         implicit none
         integer, intent(in) :: natm
         integer, intent(in) :: z(natm)
         real(dp), intent(in) :: xyz(3, natm)
         integer, intent(out) :: residues(natm)
         type(error_t), intent(inout) :: err
      end subroutine run_build_opt_residues

      module subroutine run_parse_opt_global_options(r)
         implicit none
         type(run_t), intent(inout) :: r
      end subroutine run_parse_opt_global_options

      module subroutine run_parse_frozen_list(spec, natm, frozen, err)
         implicit none
         character(len=*), intent(in) :: spec
         integer, intent(in) :: natm
         integer, allocatable, intent(out) :: frozen(:)
         type(error_t), intent(inout) :: err
      end subroutine run_parse_frozen_list

      module subroutine run_parse_constraints(spec, natm, cons_kind, cons_a, cons_b, cons_c, cons_d, err)
         implicit none
         character(len=*), intent(in) :: spec
         integer, intent(in) :: natm
         integer, allocatable, intent(out) :: cons_kind(:), cons_a(:), cons_b(:), cons_c(:), cons_d(:)
         type(error_t), intent(inout) :: err
      end subroutine run_parse_constraints

      module subroutine run_check_opt_constrain_residues(r, opt_residues)
         implicit none
         type(run_t), intent(in) :: r
         integer, intent(in) :: opt_residues(:)
      end subroutine run_check_opt_constrain_residues

      module subroutine run_check_neb_endpoint(natm, z, z2, natm2, err, what)
         implicit none
         integer, intent(in) :: natm, natm2
         integer, intent(in) :: z(natm), z2(natm2)
         type(error_t), intent(inout) :: err
         character(len=*), intent(in), optional :: what
      end subroutine run_check_neb_endpoint

      module subroutine run_read_xyz_frames(path, natm, z, nframes, xyz, err)
         implicit none
         character(len=*), intent(in) :: path
         integer, intent(in) :: natm
         integer, intent(in) :: z(natm)
         integer, intent(out) :: nframes
         real(dp), allocatable, intent(out) :: xyz(:, :, :)   !! (3, natm, nframes), Bohr
         type(error_t), intent(inout) :: err
      end subroutine run_read_xyz_frames

      module subroutine run_write_neb_path_xyz(r, path)
         implicit none
         type(run_t), intent(in) :: r
         character(len=*), intent(in) :: path
      end subroutine run_write_neb_path_xyz

      module subroutine run_opt_summary(r)
         implicit none
         type(run_t), intent(in) :: r
      end subroutine run_opt_summary

      pure module function opt_coords_name(kind) result(name)
         implicit none
         integer, intent(in) :: kind
         character(len=:), allocatable :: name
      end function opt_coords_name

      pure module function opt_target_name(kind) result(name)
         implicit none
         integer, intent(in) :: kind
         character(len=:), allocatable :: name
      end function opt_target_name

      pure module function opt_saddle_method_name(kind) result(name)
         implicit none
         integer, intent(in) :: kind
         character(len=:), allocatable :: name
      end function opt_saddle_method_name

      pure module function opt_algorithm_name(kind) result(name)
         implicit none
         integer, intent(in) :: kind
         character(len=:), allocatable :: name
      end function opt_algorithm_name

      module subroutine run_write_opt_xyz(r, path, xyz_bohr)
         implicit none
         type(run_t), intent(in) :: r
         character(len=*), intent(in) :: path
         real(dp), intent(in) :: xyz_bohr(3, r%natm)
      end subroutine run_write_opt_xyz

      module subroutine run_write_labelled_xyz(r, path, xyz_bohr, e, label)
         implicit none
         type(run_t), intent(in) :: r
         character(len=*), intent(in) :: path
         real(dp), intent(in) :: xyz_bohr(3, r%natm)
         real(dp), intent(in) :: e
         character(len=*), intent(in) :: label
      end subroutine run_write_labelled_xyz

      module subroutine run_write_traj_xyz(r, path, traj_xyz_bohr, traj_energy)
         implicit none
         type(run_t), intent(in) :: r
         character(len=*), intent(in) :: path
         real(dp), intent(in) :: traj_xyz_bohr(3, r%natm, r%opt_n_traj)
         real(dp), intent(in) :: traj_energy(r%opt_n_traj)
      end subroutine run_write_traj_xyz

      module subroutine run_write_json(r, path)
         implicit none
         type(run_t), intent(in) :: r
         character(len=*), intent(in) :: path
      end subroutine run_write_json

      module function method(r) result(m)
         implicit none
         type(run_t), intent(in) :: r
         character(len=:), allocatable :: m
      end function method
   end interface

contains

   subroutine run_die(r, msg)
      !! Every caller of `die` in this file reaches it the same way on
      !! every rank -- the same argv, the same files, the same collective
      !! library status -- so printing on rank 0 alone says it once, not
      !! `nranks` times, and every other rank is about to call `quit` with
      !! the same `code` regardless.
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: msg
      if (r%rank0) call global_logger%error("terco: "//msg)
      call run_quit(r, 1)
   end subroutine run_die

   subroutine run_quit(r, code)
      !! End the process with exit status `code`, printing nothing more.
      !! `stop code, quiet=.true.` is the standard way to say this, but
      !! nvfortran rejects it, and a plain `stop code` prints the code.
      !!
      !! More than one rank and a non-zero code: `abort_comm` rather than a
      !! clean `pic_mpi_finalize`, because a caller of `die` cannot promise
      !! every OTHER rank reached the same conclusion the same way -- and
      !! `MPI_Abort` tears down the whole job either way, so there is
      !! nothing to wait for. One rank keeps today's path exactly: no pic-mpi
      !! call before it that was not already a no-op there.
      type(run_t), intent(in) :: r
      integer, intent(in) :: code
      interface
         subroutine c_exit(status) bind(c, name="exit")
            import :: c_int
            implicit none
            integer(c_int), value :: status
         end subroutine c_exit
      end interface
      flush (output_unit)
      flush (error_unit)
      if (code /= 0 .and. r%nranks > 1) call abort_comm(r%wcomm, code)
      call pic_mpi_finalize()
      call c_exit(int(code, c_int))
   end subroutine run_quit

   subroutine run_check(r, what)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: what
      if (r%err%has_error()) call run_die(r, what//": "//r%err%get_message())
   end subroutine run_check

   subroutine run_api(r, status)
      !! Stop on a failed library call, with the library's own explanation.
      type(run_t), intent(in) :: r
      integer(c_int), intent(in) :: status
      integer, parameter :: MSG_LEN = 1024
      character(kind=c_char) :: buffer(MSG_LEN)
      integer(c_int) :: s
      if (status == TRC_OK) return
      s = trc_message(r%h, buffer, int(size(buffer), c_int))
      call run_die(r, app_f_string(buffer))
   end subroutine run_api

   integer function run_int_option(r, name, file, default) result(v)
      !! An integer option, else the first integer in `file` in the working
      !! directory (xtb's .CHRG/.UHF convention), else `default`.
      type(run_t), intent(inout) :: r
      character(len=*), intent(in) :: name, file
      integer, intent(in) :: default
      integer :: u, ios
      logical :: there
      v = default
      if (r%cli%is_set(name)) then
         call r%cli%get(name, v, r%err)
         call run_check(r, "--"//name)
         return
      end if
      inquire (file=file, exist=there)
      if (.not. there) return
      open (newunit=u, file=file, status="old", action="read", iostat=ios)
      if (ios /= 0) return
      read (u, *, iostat=ios) v
      close (u)
      if (ios /= 0) call run_die(r, file//" does not start with an integer")
   end function run_int_option

   function run_to_lower_alloc(s) result(lower_s)
      !! A non-elemental, allocatable-result ASCII lowercaser for every CLI
      !! normalisation site in `trc_main` below (`--func`, `--guess`,
      !! `--masses`, `--pcm`, the `--opt-*` string options, and the geometry
      !! file-extension check). It replaces pic's `to_lower`
      !! (`pic_ascii.f90`): nvfortran 26.5 and 26.9, host and GPU builds
      !! alike, with `trc_main`'s real `-fast -O3` flags, miscompile
      !! pic's `to_lower` -- `pure elemental`, automatic-length (not
      !! allocatable) result -- at those call sites. The result comes back
      !! with the right length but blank (all-space) content, so e.g.
      !! `--func pbe` silently falls back to a plain SCF with no exchange-
      !! correlation term. A non-elemental function with an allocatable
      !! result, called the same way, is unaffected; see
      !! `.claude/context/gpu-pitfalls.md` for the trap. `trc_main` has no
      !! internal procedures (see this module's header comment), so this
      !! lives here instead.
      character(len=*), intent(in) :: s
      character(len=:), allocatable :: lower_s
      integer :: i, c
      allocate (character(len=len(s)) :: lower_s)
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) then
            lower_s(i:i) = achar(c + 32)
         else
            lower_s(i:i) = s(i:i)
         end if
      end do
   end function run_to_lower_alloc

   subroutine run_to_lower_into(s, lower_s)
      !! `run_to_lower_alloc` as a subroutine: ASCII-lowercase `s` into an
      !! `intent(out)` deferred-length allocatable dummy. Use this, not the
      !! function, wherever the result lands in a deferred-length allocatable
      !! (a `run_t` component or a local): nvfortran (CI's nvidia-hpc host job)
      !! returns corrupted contents (right length, garbage bytes) when a
      !! deferred-length allocatable CHARACTER function result is assigned to
      !! a deferred-length allocatable in a submodule procedure of this
      !! module, e.g. `r%func = run_to_lower_alloc(trim(sval))` in
      !! `run_parse_inputs`. The subroutine form allocates the target itself.
      character(len=*), intent(in) :: s
      character(len=:), allocatable, intent(out) :: lower_s
      integer :: i, c
      allocate (character(len=len(s)) :: lower_s)
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) then
            lower_s(i:i) = achar(c + 32)
         else
            lower_s(i:i) = s(i:i)
         end if
      end do
   end subroutine run_to_lower_into

end module trc_run
