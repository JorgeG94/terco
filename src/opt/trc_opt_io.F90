!
! Small, self-contained helpers trc_dlfind_bridge.F90 needs: NUL-terminated
! C strings, POSIX setenv, and the three XYZ/JSON writers DL-FIND's own
! `error stop` forces to run from INSIDE a callback (trc_dlfind_bridge's own
! fail_eval -- there is no "after api_dl_find returns" on that path).
!
! WHY THIS DUPLICATES PART OF app/trc_app.F90 RATHER THAN MOVING IT THERE.
! opt_io_f_string/opt_io_set_env/opt_io_write_opt_xyz_frame/opt_io_write_traj_xyz_frames/
! opt_io_write_opt_failure_json used to live in trc_app (app/trc_app.F90),
! which trc_dlfind_bridge.F90 could reach because both were compiled
! straight into the terco_cli executable together. Now that the bridge is
! built into libterco_opt (TERCO_ENABLE_DLFIND=ON only), it can no longer
! reach anything in app/: CODE_STRUCTURE.md's "entry points are leaves"
! rule, and not just a style point here -- app/ is not installed or
! exported, so a caller outside terco's own build (cyclone, linking only
! libterco_opt.so and terco_opt.h) could not link against it even if the
! rule were bent. Moving these five routines OUT of trc_app instead would
! break the OPPOSITE build: app/main.F90's ordinary (non-crash) --opt finish
! path (write_opt_xyz/write_traj_xyz) calls them unconditionally, compiled
! into terco_cli in EVERY build, DLFIND=ON or OFF, so trc_app.F90 must
! keep its own copies regardless. Five small file-writers, stable since
! they were written, duplicated once rather than reshuffling module
! boundaries across both the CLI and the library for them.
!
module trc_opt_io
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use pic_types, only: dp
   use trc_constants, only: ANGSTROM_TO_BOHR
   use json_module, only: json_core, json_value
   implicit none
   private

   public :: opt_io_c_string, opt_io_f_string, opt_io_set_env
   public :: opt_io_write_opt_xyz_frame, opt_io_write_traj_xyz_frames, opt_io_write_opt_failure_json
   public :: trc_opt_failure_t

   type :: trc_opt_failure_t
    !! The scalar fields of a failed optimization's JSON block, bundled so
    !! opt_io_write_opt_failure_json stays under the dummy-argument limit.
    !! `path` and `traj_energy` (a run's own array) stay separate arguments.
      character(len=:), allocatable :: geometry, reason, coords_name, guess_name
      integer :: steps = 0
      integer :: energy_evaluations = 0
      real(dp) :: final_energy = 0.0_dp
      real(dp) :: final_gradient_max = 0.0_dp
      real(dp) :: final_gradient_rms = 0.0_dp
   end type trc_opt_failure_t

   real(dp), parameter :: BOHR_PER_ANGSTROM = ANGSTROM_TO_BOHR
      !! CODATA 2018 (trc_constants' ANGSTROM_TO_BOHR) -- the same constant
      !! app/trc_app.F90's own copy of this name uses.

   integer, parameter :: NELEM = 118
   character(len=2), parameter :: SYMBOLS(NELEM) = [character(len=2) :: &
      "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
      "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", &
      "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
      "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", &
      "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
      "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
      "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
      "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
      "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
      "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", &
      "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

   interface
      function c_setenv(name, value, overwrite) result(rc) bind(c, name="setenv")
         import :: c_char, c_int
         implicit none
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: rc
      end function c_setenv
   end interface

contains

   pure function element_symbol(z) result(sym)
      !! "??" for an out-of-range Z rather than an out-of-bounds read -- a
      !! malformed C caller's atomic number is not this routine's job to
      !! validate further than that.
      integer, intent(in) :: z
      character(len=2) :: sym
      if (z >= 1 .and. z <= NELEM) then
         sym = SYMBOLS(z)
      else
         sym = "??"
      end if
   end function element_symbol

   pure function opt_io_c_string(s) result(cs)
      !! `s` NUL-terminated, for a `character(kind=c_char) :: x(*)` dummy.
      character(len=*), intent(in) :: s
      character(kind=c_char, len=len_trim(s) + 1) :: cs
      cs = trim(s)//c_null_char
   end function opt_io_c_string

   pure function opt_io_f_string(buffer) result(s)
      !! The text of a NUL-terminated C buffer.
      character(kind=c_char), intent(in) :: buffer(:)
      character(len=:), allocatable :: s
      integer :: i, n
      n = size(buffer)
      do i = 1, size(buffer)
         if (buffer(i) == c_null_char) then
            n = i - 1
            exit
         end if
      end do
      allocate (character(len=n) :: s)
      do i = 1, n
         s(i:i) = buffer(i)
      end do
   end function opt_io_f_string

   subroutine opt_io_set_env(name, value)
      !! POSIX setenv -- trc_dlfind_bridge's own DLF_PO_SEED (--opt-seed).
      character(len=*), intent(in) :: name, value
      integer(c_int) :: rc
      rc = c_setenv(opt_io_c_string(name), opt_io_c_string(value), 1_c_int)
   end subroutine opt_io_set_env

   subroutine opt_io_write_opt_xyz_frame(path, natm, z, energy, xyz_bohr, label)
      !! `<geometry base>_opt.xyz`: one geometry in Angstrom, atom order
      !! unchanged, comment line carrying `label` and the energy. Shared by
      !! an ordinary --opt finish (app/main.F90's own write_opt_xyz, which
      !! calls trc_app's OWN identical copy) and trc_dlfind_bridge's
      !! unrecoverable-failure path, which has no other chance to write
      !! this file: DL-FIND's own api.f90 calls `error stop` right after a
      !! gradient evaluation reports status=1, so nothing runs after
      !! api_dl_find returns on that path.
      character(len=*), intent(in) :: path
      integer, intent(in) :: natm
      integer, intent(in) :: z(natm)
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: xyz_bohr(3, natm)
      character(len=*), intent(in) :: label
      integer :: u, ia
      open (newunit=u, file=path, status="replace", action="write")
      write (u, "(i0)") natm
      write (u, "(a,a,f0.12)") trim(label), "   E = ", energy
      do ia = 1, natm
         write (u, "(a2,3f20.12)") element_symbol(z(ia)), xyz_bohr(:, ia)/BOHR_PER_ANGSTROM
      end do
      close (u)
   end subroutine opt_io_write_opt_xyz_frame

   subroutine opt_io_write_traj_xyz_frames(path, natm, z, n_traj, traj_xyz_bohr, traj_energy)
      !! `<geometry base>_traj.xyz`: one frame per geometry handed in, in
      !! Angstrom, comment line carrying the step index and energy.
      character(len=*), intent(in) :: path
      integer, intent(in) :: natm, n_traj
      integer, intent(in) :: z(natm)
      real(dp), intent(in) :: traj_xyz_bohr(3, natm, n_traj)
      real(dp), intent(in) :: traj_energy(n_traj)
      integer :: u, ia, k
      open (newunit=u, file=path, status="replace", action="write")
      do k = 1, n_traj
         write (u, "(i0)") natm
         write (u, "(a,i0,a,f0.12)") "step ", k, "   E = ", traj_energy(k)
         do ia = 1, natm
            write (u, "(a2,3f20.12)") element_symbol(z(ia)), traj_xyz_bohr(:, ia, k)/BOHR_PER_ANGSTROM
         end do
      end do
      close (u)
   end subroutine opt_io_write_traj_xyz_frames

   subroutine opt_io_write_opt_failure_json(path, info, traj_energy)
      !! A minimal "optimization" JSON block for the one outcome
      !! app/main.F90's own add_optimization never gets a chance to write:
      !! DL-FIND's `error stop` right after a failed evaluation.
      character(len=*), intent(in) :: path
      type(trc_opt_failure_t), intent(in) :: info
      real(dp), intent(in) :: traj_energy(:)
      type(json_core) :: json
      type(json_value), pointer :: root, oobj
      call json%initialize()
      call json%create_object(root, "")
      call json%add(root, "program", "terco")
      call json%add(root, "geometry", info%geometry)
      call json%create_object(oobj, "optimization")
      call json%add(oobj, "converged", .false.)
      call json%add(oobj, "failed", .true.)
      call json%add(oobj, "failure_reason", info%reason)
      call json%add(oobj, "steps", info%steps)
      call json%add(oobj, "energy_evaluations", info%energy_evaluations)
      call json%add(oobj, "final_energy", info%final_energy)
      call json%add(oobj, "final_gradient_max", info%final_gradient_max)
      call json%add(oobj, "final_gradient_rms", info%final_gradient_rms)
      call json%add(oobj, "coordinates", info%coords_name)
      call json%add(oobj, "guess", info%guess_name)
      if (size(traj_energy) > 0) call json%add(oobj, "trajectory_energies", traj_energy)
      call json%add(root, oobj)
      call json%print(root, path)
      call json%destroy(root)
   end subroutine opt_io_write_opt_failure_json

end module trc_opt_io
