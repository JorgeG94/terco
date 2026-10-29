!
! The Hirshfeld free-atom (proatom) radial density tables.
!
! Ported from metalquicha's mqc_czt_hirshfeld (the table side only -- see
! trc_hirshfeld.F90's own header for the partition itself), mqc_ -> trc_
! renamed. mqc solves each free atom itself, at run time, and caches it in
! memory; terco instead ships the already-converged tables as a data file,
! test/basis_sets/hirshfeld_proatoms.json, generated once offline by
! scripts/gen_hirshfeld_proatoms.py and loaded here the first time any
! element is needed -- the same "basis_sets/<name>.json" relative default
! plus an environment override that trc_minao's MINAO_BASIS_DEFAULT/
! TRC_MINAO_BASIS and trc_sap's SAP_BASIS_DEFAULT/TRC_SAP_BASIS already use.
!
module trc_hirshfeld_data
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   use json_module, only: json_file
   implicit none
   private

   public :: trc_hirshfeld_table_t
   public :: trc_hirshfeld_ensure_loaded
   public :: HIRSHFELD_Z_MAX

   integer, parameter :: HIRSHFELD_Z_MAX = 36
      !! Krypton, the last element the free-atom recipe covers (H-Kr).

   character(len=*), parameter :: PROATOMS_DEFAULT = "basis_sets/hirshfeld_proatoms.json"

   type :: trc_hirshfeld_table_t
      !! Every element's radial density, uniform in ln r, ready for cubic
      !! Hermite interpolation in (ln r, ln rho) -- see
      !! trc_hirshfeld.F90's `hirshfeld_table_density`, which is where the
      !! interpolation itself lives (it has to be in the same module as the
      !! `do concurrent` kernel that calls it).
      logical :: loaded = .false.
      integer :: n_radial = 0
      real(dp) :: log_r0 = 0.0_dp
      real(dp) :: step = 0.0_dp
      logical :: valid(HIRSHFELD_Z_MAX) = .false.
      real(dp) :: r_cut(HIRSHFELD_Z_MAX) = 0.0_dp
      real(dp), allocatable :: log_rho(:, :)   !! (n_radial, Z_MAX)
      real(dp), allocatable :: slope(:, :)     !! (n_radial, Z_MAX)
   end type trc_hirshfeld_table_t

contains

   subroutine trc_hirshfeld_ensure_loaded(table, error, path)
      !! Loads the proatom data file into `table` unless it is already
      !! loaded. Every element present in the file is read at once: the
      !! file is one unit (2.7 MB for H-Kr), not a per-element fetch, so
      !! there is nothing to cache incrementally.
      type(trc_hirshfeld_table_t), intent(inout) :: table
      type(error_t), intent(inout) :: error
      character(len=*), intent(in), optional :: path
      !! Overrides TRC_HIRSHFELD_DATA/PROATOMS_DEFAULT, for a caller (a
      !! test) that wants a specific file regardless of environment.

      type(json_file) :: json
      character(len=:), allocatable :: data_path
      character(len=512) :: env
      integer :: z
      logical :: found
      real(dp) :: r_min, r_max

      if (error%has_error()) return
      if (table%loaded) return

      data_path = PROATOMS_DEFAULT
      env = " "
      call get_environment_variable("TRC_HIRSHFELD_DATA", env)
      if (len_trim(env) > 0) data_path = trim(env)
      if (present(path)) data_path = path

      call json%initialize()
      call json%load(filename=data_path)
      if (json%failed()) then
         call error%set(ERROR_VALIDATION, "Hirshfeld charges: cannot read "//data_path// &
                        " (set TRC_HIRSHFELD_DATA, or regenerate with "// &
                        "scripts/gen_hirshfeld_proatoms.py)")
         call json%destroy()
         return
      end if

      call json%get("n_radial", table%n_radial, found)
      call json%get("r_min", r_min, found)
      call json%get("r_max", r_max, found)
      if (table%n_radial <= 1 .or. r_min <= 0.0_dp .or. r_max <= r_min) then
         call error%set(ERROR_VALIDATION, "Hirshfeld charges: "//data_path//" has a "// &
                        "malformed radial grid")
         call json%destroy()
         return
      end if
      table%log_r0 = log(r_min)
      table%step = (log(r_max) - log(r_min))/real(table%n_radial - 1, dp)

      allocate (table%log_rho(table%n_radial, HIRSHFELD_Z_MAX))
      allocate (table%slope(table%n_radial, HIRSHFELD_Z_MAX))
      table%log_rho = 0.0_dp
      table%slope = 0.0_dp
      do z = 1, HIRSHFELD_Z_MAX
         call read_element(json, data_path, z, table, error)
         if (error%has_error()) then
            call json%destroy()
            return
         end if
      end do
      call json%destroy()
      table%loaded = .true.
   end subroutine trc_hirshfeld_ensure_loaded

   subroutine read_element(json, data_path, z, table, error)
      ! One element's table, left %valid(z) = .false. (not an error) when
      ! the file simply has no entry for it -- trc_hirshfeld's own refusal
      ! for "no free atom for this Z" reads that flag.
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: data_path
      integer, intent(in) :: z
      type(trc_hirshfeld_table_t), intent(inout) :: table
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: key
      real(dp), allocatable :: log_rho(:), slope(:)
      logical :: found
      integer :: n_children

      key = "elements."//to_string(z)
      call json%info(key//".log_rho", found=found, n_children=n_children)
      if (.not. found) return
      if (n_children /= table%n_radial) then
         call error%set(ERROR_VALIDATION, "Hirshfeld charges: "//data_path//" element "// &
                        to_string(z)//" has "//to_string(n_children)//" radial points, "// &
                        "not "//to_string(table%n_radial))
         return
      end if

      call json%get(key//".r_cut", table%r_cut(z), found)
      call json%get(key//".log_rho", log_rho, found)
      call json%get(key//".slope", slope, found)
      if (.not. found .or. size(log_rho) /= table%n_radial .or. size(slope) /= table%n_radial) then
         call error%set(ERROR_VALIDATION, "Hirshfeld charges: "//data_path//" element "// &
                        to_string(z)//" is missing r_cut/log_rho/slope")
         return
      end if
      table%log_rho(:, z) = log_rho
      table%slope(:, z) = slope
      table%valid(z) = .true.
   end subroutine read_element

end module trc_hirshfeld_data
