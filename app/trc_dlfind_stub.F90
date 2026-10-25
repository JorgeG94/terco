!
! The --opt stub for a build with TERCO_ENABLE_DLFIND=OFF (the default --
! see that option's own comment in the root CMakeLists.txt). Same public
! interface as app/trc_dlfind_bridge.F90 (the TERCO_ENABLE_DLFIND=ON
! build), so app/main.F90 never branches on which of the two it was built
! against: it calls trc_dlfind_available() first and, finding it false,
! refuses --opt with a message instead of ever calling trc_dlfind_optimize.
!
module trc_dlfind_bridge
   use, intrinsic :: iso_c_binding, only: c_int, c_ptr
   use pic_types, only: dp
   implicit none
   private

   public :: trc_dlfind_available
   public :: TRC_OPT_COORDS_CARTESIAN, TRC_OPT_COORDS_HDLC, TRC_OPT_COORDS_HDLC_TC, &
             TRC_OPT_COORDS_DLC, TRC_OPT_COORDS_DLC_TC
   public :: TRC_OPT_GUESS_DENSITY, TRC_OPT_GUESS_PROJECTED, TRC_OPT_GUESS_NONE
   public :: TRC_OPT_ALG_SD, TRC_OPT_ALG_CG_AUTO, TRC_OPT_ALG_CG, TRC_OPT_ALG_LBFGS, &
             TRC_OPT_ALG_PRFO, TRC_OPT_ALG_NR, TRC_OPT_ALG_DAMPED
   public :: TRC_OPT_TARGET_MINIMUM, TRC_OPT_TARGET_SADDLE
   public :: TRC_OPT_SADDLE_PRFO, TRC_OPT_SADDLE_DIMER
   public :: TRC_OPT_HESS_TERCO, TRC_OPT_HESS_DLFIND
   public :: TRC_OPT_HUPD_NONE, TRC_OPT_HUPD_POWELL, TRC_OPT_HUPD_BOFILL, TRC_OPT_HUPD_BFGS
   public :: TRC_OPT_LINE_AUTO, TRC_OPT_LINE_SCALING, TRC_OPT_LINE_ENERGY_TRUST, TRC_OPT_LINE_GRADIENT_TRUST
   public :: trc_opt_options_t, trc_opt_result_t
   public :: trc_dlfind_optimize

   integer, parameter :: TRC_OPT_COORDS_CARTESIAN = 0, TRC_OPT_COORDS_HDLC = 1, TRC_OPT_COORDS_HDLC_TC = 2, &
                         TRC_OPT_COORDS_DLC = 3, TRC_OPT_COORDS_DLC_TC = 4
   integer, parameter :: TRC_OPT_GUESS_DENSITY = 0, TRC_OPT_GUESS_PROJECTED = 1, TRC_OPT_GUESS_NONE = 2
   integer, parameter :: TRC_OPT_ALG_SD = 0, TRC_OPT_ALG_CG_AUTO = 1, TRC_OPT_ALG_CG = 2, TRC_OPT_ALG_LBFGS = 3, &
                         TRC_OPT_ALG_PRFO = 10, TRC_OPT_ALG_NR = 20, TRC_OPT_ALG_DAMPED = 30
   integer, parameter :: TRC_OPT_TARGET_MINIMUM = 0, TRC_OPT_TARGET_SADDLE = 1
   integer, parameter :: TRC_OPT_SADDLE_PRFO = 0, TRC_OPT_SADDLE_DIMER = 1
   integer, parameter :: TRC_OPT_HESS_TERCO = 0, TRC_OPT_HESS_DLFIND = 1
   integer, parameter :: TRC_OPT_HUPD_NONE = 0, TRC_OPT_HUPD_POWELL = 1, TRC_OPT_HUPD_BOFILL = 2, &
                         TRC_OPT_HUPD_BFGS = 3
   integer, parameter :: TRC_OPT_LINE_AUTO = -1, TRC_OPT_LINE_SCALING = 0, TRC_OPT_LINE_ENERGY_TRUST = 1, &
                         TRC_OPT_LINE_GRADIENT_TRUST = 2

   type :: trc_opt_options_t
      integer :: coords_kind = TRC_OPT_COORDS_CARTESIAN
      integer :: guess_kind = TRC_OPT_GUESS_DENSITY
      integer(c_int) :: fallback_guess_kind = 2_c_int
      integer :: maxiter = 100
      real(dp) :: gtol = 4.5e-4_dp
      integer :: algorithm = TRC_OPT_ALG_LBFGS
      integer :: target_kind = TRC_OPT_TARGET_MINIMUM
      integer :: saddle_method = TRC_OPT_SADDLE_PRFO
      logical :: connect = .false.
      logical :: is_linear = .false.
      integer :: hessian_source = TRC_OPT_HESS_TERCO
      integer :: hessian_update = -1
      integer :: max_updates = -1
      integer :: zero_modes = -1
      real(dp) :: soft_modes = -1.0_dp
      real(dp) :: max_step = -1.0_dp
      integer :: linesearch = TRC_OPT_LINE_AUTO
      real(dp) :: scalestep = -1.0_dp
      integer :: lbfgs_memory = -1
      real(dp) :: etol = -1.0_dp
      integer :: print_level = -1
      character(len=:), allocatable :: opt_xyz_path
      character(len=:), allocatable :: traj_xyz_path
      character(len=:), allocatable :: json_path
      character(len=:), allocatable :: geometry_path
      character(len=:), allocatable :: coords_name
      character(len=:), allocatable :: guess_name
   end type trc_opt_options_t

   type :: trc_opt_result_t
      logical :: converged = .false.
      logical :: failed = .false.
      character(len=:), allocatable :: message
      integer :: steps = 0
      integer :: energy_evaluations = 0
      integer :: hessian_evaluations = 0
      real(dp) :: final_energy = 0.0_dp
      real(dp) :: final_gmax = 0.0_dp
      real(dp) :: final_grms = 0.0_dp
      real(dp), allocatable :: final_coords(:, :)
      integer :: n_traj = 0
      real(dp), allocatable :: traj_energy(:)
      real(dp), allocatable :: traj_coords(:, :, :)
      integer :: total_scf_iterations = 0
      logical :: have_ts = .false.
      real(dp) :: ts_energy = 0.0_dp
      real(dp) :: ts_gmax = 0.0_dp, ts_grms = 0.0_dp
      logical :: ts_converged = .false.
      real(dp), allocatable :: ts_coords(:, :)
      logical :: have_minimum_plus = .false., have_minimum_minus = .false.
      real(dp) :: minimum_plus_energy = 0.0_dp, minimum_minus_energy = 0.0_dp
      real(dp) :: minimum_plus_gmax = 0.0_dp, minimum_plus_grms = 0.0_dp
      real(dp) :: minimum_minus_gmax = 0.0_dp, minimum_minus_grms = 0.0_dp
      logical :: minimum_plus_converged = .false., minimum_minus_converged = .false.
      real(dp), allocatable :: minimum_plus_coords(:, :), minimum_minus_coords(:, :)
   end type trc_opt_result_t

contains

   pure function trc_dlfind_available() result(ok)
      !! Always false: this is the TERCO_ENABLE_DLFIND=OFF build.
      logical :: ok
      ok = .false.
   end function trc_dlfind_available

   subroutine trc_dlfind_optimize(handle, n, z, residues_in, coords, opts, res)
      !! Never actually called -- app/main.F90 checks trc_dlfind_available()
      !! before this point and refuses --opt itself. Present anyway (rather
      !! than left out of this stub) so the two bridge sources keep IDENTICAL
      !! public interfaces and main.F90 need not know which one it was built
      !! against; calling it regardless is a refusal, not a crash, in case a
      !! future caller skips the availability check.
      type(c_ptr), intent(in) :: handle
      integer, intent(in) :: n
      integer, intent(in) :: z(n)
      integer, intent(in) :: residues_in(n)
      real(dp), intent(in) :: coords(3, n)
      type(trc_opt_options_t), intent(in) :: opts
      type(trc_opt_result_t), intent(out) :: res
      res%failed = .true.
      res%message = "terco was built with -DTERCO_ENABLE_DLFIND=OFF (the default): "// &
                    "--opt needs libdlfind. Reconfigure with -DTERCO_ENABLE_DLFIND=ON."
   end subroutine trc_dlfind_optimize

end module trc_dlfind_bridge
