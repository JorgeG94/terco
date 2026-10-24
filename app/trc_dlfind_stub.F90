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
   public :: TRC_OPT_COORDS_CARTESIAN, TRC_OPT_COORDS_HDLC, TRC_OPT_COORDS_DLC
   public :: TRC_OPT_GUESS_DENSITY, TRC_OPT_GUESS_PROJECTED, TRC_OPT_GUESS_NONE
   public :: trc_opt_options_t, trc_opt_result_t
   public :: trc_dlfind_optimize

   integer, parameter :: TRC_OPT_COORDS_CARTESIAN = 0, TRC_OPT_COORDS_HDLC = 1, TRC_OPT_COORDS_DLC = 2
   integer, parameter :: TRC_OPT_GUESS_DENSITY = 0, TRC_OPT_GUESS_PROJECTED = 1, TRC_OPT_GUESS_NONE = 2

   type :: trc_opt_options_t
      integer :: coords_kind = TRC_OPT_COORDS_CARTESIAN
      integer :: guess_kind = TRC_OPT_GUESS_DENSITY
      integer(c_int) :: fallback_guess_kind = 2_c_int
      integer :: maxiter = 100
      real(dp) :: gtol = 4.5e-4_dp
      real(dp) :: max_step = -1.0_dp
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
      real(dp) :: final_energy = 0.0_dp
      real(dp) :: final_gmax = 0.0_dp
      real(dp) :: final_grms = 0.0_dp
      real(dp), allocatable :: final_coords(:, :)
      integer :: n_traj = 0
      real(dp), allocatable :: traj_energy(:)
      real(dp), allocatable :: traj_coords(:, :, :)
      integer :: total_scf_iterations = 0
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
