module check_capi_timing_suite
   !! trc_timing and trc_timing_count through the C ABI: the status before any
   !! SCF, an enumeration that ends by its count (not by an error), a name
   !! buffer too short to hold a name being refused rather than cut, and the
   !! RI-MP2 timings not outliving the run they belong to (a later trc_run_scf,
   !! trc_load_scf or trc_run_rimp2 must not report the previous run's rimp2.*).
   !! Values are wall times, so nothing here pins a number.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_convergence, trc_set_rimp2, trc_run_scf, trc_run_rimp2, trc_rimp2_energy, &
                               trc_save_scf, trc_load_scf, trc_timing, trc_timing_count, trc_message, &
                               TRC_OK, TRC_ERR_BADARG, TRC_ERR_STATE, TRC_TIMING_NAME_MAX
   implicit none
   private
   public :: collect_capi_timing

   character(len=*), parameter :: SAVE_FILE = "check_capi_timing_scf.tmp"

contains

   subroutine collect_capi_timing(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("before any SCF: TRC_ERR_STATE from both entry points", test_before_scf), &
                  new_unittest("enumeration ends at the count; bad index and short buffer are refused", test_enumeration), &
                  new_unittest("rimp2.* does not outlive its run (run_scf, load_scf, run_rimp2)", test_stale_rimp2) &
                  ]
   end subroutine collect_capi_timing

   subroutine setup(h, error)
      !! water / cc-pVDZ / cc-pVDZ-RIFIT, loose convergence, nothing run yet.
      type(c_ptr), intent(out) :: h
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      integer(c_int) :: rc
      call read_xyz('water.xyz', natm, zint, at_r)
      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                                              0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      if (rc == TRC_OK) rc = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      if (rc == TRC_OK) rc = trc_set_convergence(h, 1.0e-8_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
      if (rc == TRC_OK) rc = trc_set_rimp2(h, -1_c_int, 0_c_int)
      call check(error, rc == TRC_OK, "setup failed, status "//to_string(int(rc)))
   end subroutine setup

   function entry_name(h, i, status) result(name)
      !! Entry i's name via a TRC_TIMING_NAME_MAX buffer; `status` is trc_timing's.
      type(c_ptr), intent(in) :: h
      integer, intent(in) :: i
      integer(c_int), intent(out) :: status
      character(len=:), allocatable :: name
      character(kind=c_char) :: buf(TRC_TIMING_NAME_MAX)
      real(c_double) :: v
      integer :: k
      status = trc_timing(h, int(i, c_int), buf, TRC_TIMING_NAME_MAX, v)
      name = ""
      if (status /= TRC_OK) return
      do k = 1, TRC_TIMING_NAME_MAX
         if (buf(k) == c_null_char) exit
         name = name//buf(k)
      end do
   end function entry_name

   function has_prefix(h, prefix) result(found)
      !! Whether any listed name starts with `prefix`.
      type(c_ptr), intent(in) :: h
      character(len=*), intent(in) :: prefix
      logical :: found
      integer(c_int) :: n, rc
      integer :: i
      character(len=:), allocatable :: nm
      found = .false.
      rc = trc_timing_count(h, n)
      if (rc /= TRC_OK) return
      do i = 0, n - 1
         nm = entry_name(h, i, rc)
         if (rc == TRC_OK) then
            if (nm(1:min(len(nm), len(prefix))) == prefix) found = .true.
         end if
      end do
   end function has_prefix

   subroutine test_before_scf(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, n
      character(kind=c_char) :: buf(TRC_TIMING_NAME_MAX)
      real(c_double) :: v
      call setup(h, error)
      if (allocated(error)) return
      rc = trc_timing_count(h, n)
      call check(error, rc == TRC_ERR_STATE, "trc_timing_count before any SCF: status "//to_string(int(rc)))
      if (.not. allocated(error)) then
         rc = trc_timing(h, 0_c_int, buf, TRC_TIMING_NAME_MAX, v)
         call check(error, rc == TRC_ERR_STATE, "trc_timing before any SCF: status "//to_string(int(rc)))
      end if
      rc = trc_destroy(h)
   end subroutine test_before_scf

   subroutine test_enumeration(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, n, i
      character(len=:), allocatable :: nm
      character(kind=c_char) :: buf(TRC_TIMING_NAME_MAX), small(4)
      character(len=256) :: msg
      real(c_double) :: v
      call setup(h, error)
      if (allocated(error)) return
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the SCF failed")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_timing_count(h, n)
      call check(error, rc == TRC_OK .and. n >= 6, "trc_timing_count: status "//to_string(int(rc))// &
                 ", n = "//to_string(int(n)))
      if (.not. allocated(error)) then
         do i = 0, n - 1
            nm = entry_name(h, int(i), rc)
            call check(error, rc == TRC_OK, "entry "//to_string(int(i))//" refused, status "//to_string(int(rc)))
            if (allocated(error)) exit
            call check(error, len(nm) > 4, "entry "//to_string(int(i))//" has no name")
            if (allocated(error)) exit
         end do
      end if
      ! The end of the list is a bad argument, and it is the count that says where it is.
      if (.not. allocated(error)) then
         nm = entry_name(h, int(n), rc)
         call check(error, rc == TRC_ERR_BADARG, "index n: status "//to_string(int(rc)))
      end if
      if (.not. allocated(error)) then
         rc = trc_timing(h, -1_c_int, buf, TRC_TIMING_NAME_MAX, v)
         call check(error, rc == TRC_ERR_BADARG, "index -1: status "//to_string(int(rc)))
      end if
      ! A buffer that cannot hold the first name and its NUL is refused, with the length needed.
      if (.not. allocated(error)) then
         rc = trc_timing(h, 0_c_int, small, 4_c_int, v)
         call check(error, rc == TRC_ERR_BADARG, "a 4-char buffer: status "//to_string(int(rc)))
         if (.not. allocated(error)) then
            rc = trc_message(h, msg, int(len(msg), c_int))
            call check(error, index(msg, "needs") > 0, "the refusal does not say the length needed: "//trim(msg))
         end if
      end if
      if (.not. allocated(error)) call check(error, .not. has_prefix(h, "rimp2."), "rimp2.* listed with no RI-MP2 run")
      rc = trc_destroy(h)
   end subroutine test_enumeration

   subroutine test_stale_rimp2(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: eos, ess
      integer :: u
      call setup(h, error)
      if (allocated(error)) return
      rc = trc_run_scf(h)
      if (rc == TRC_OK) rc = trc_run_rimp2(h)
      if (rc == TRC_OK) rc = trc_rimp2_energy(h, eos, ess)   ! settles the deferred energy step
      call check(error, rc == TRC_OK, "SCF + RI-MP2 failed, status "//to_string(int(rc)))
      if (.not. allocated(error)) call check(error, has_prefix(h, "rimp2.total"), "no rimp2.total after the RI-MP2 energy")
      ! The same SCF, RI-MP2 asked for again: the first run's rimp2.* is not the second's
      ! (its energy step is deferred until something asks for the energy).
      if (.not. allocated(error)) then
         rc = trc_run_rimp2(h)
         call check(error, rc == TRC_OK, "the repeated trc_run_rimp2 failed")
      end if
      if (.not. allocated(error)) call check(error, .not. has_prefix(h, "rimp2."), &
                                             "rimp2.* of the previous trc_run_rimp2 survived a new one")
      ! A new SCF: the previous run's rimp2.* must be gone.
      if (.not. allocated(error)) then
         rc = trc_run_scf(h)
         call check(error, rc == TRC_OK, "the second SCF failed")
      end if
      if (.not. allocated(error)) call check(error, .not. has_prefix(h, "rimp2."), "rimp2.* survived trc_run_scf")
      ! trc_run_rimp2 defers its energy step: until something asks for the energy, nothing was timed.
      if (.not. allocated(error)) then
         rc = trc_run_rimp2(h)
         call check(error, rc == TRC_OK, "the second trc_run_rimp2 failed")
      end if
      if (.not. allocated(error)) call check(error, .not. has_prefix(h, "rimp2."), &
                                             "rimp2.* listed before the deferred energy step ran")
      if (.not. allocated(error)) then
         rc = trc_rimp2_energy(h, eos, ess)
         call check(error, rc == TRC_OK .and. has_prefix(h, "rimp2.total"), "no rimp2.total after the second energy step")
      end if
      ! trc_run_rimp2 again, then a load: the loaded result ran no SCF loop and no RI-MP2.
      if (.not. allocated(error)) then
         rc = trc_save_scf(h, SAVE_FILE//c_null_char)
         if (rc == TRC_OK) rc = trc_load_scf(h, SAVE_FILE//c_null_char)
         call check(error, rc == TRC_OK, "save/load of the SCF failed, status "//to_string(int(rc)))
      end if
      if (.not. allocated(error)) call check(error, .not. has_prefix(h, "rimp2."), "rimp2.* survived trc_load_scf")
      rc = trc_destroy(h)
      open (newunit=u, file=SAVE_FILE, status="old", iostat=rc)
      if (rc == 0) close (u, status="delete")
   end subroutine test_stale_rimp2

end module check_capi_timing_suite

program check_capi_timing
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_capi_timing_suite, only: collect_capi_timing
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_capi_timing, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_timing
