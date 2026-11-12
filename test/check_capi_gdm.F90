!
! trc_set_gdm, the C ABI switch for geometric direct minimisation, as a
! test-drive suite. Every case fails if the wiring it names is broken.
!
! 1. Bad arguments are refused (on other than 0/1, a negative or NaN start or
!    max_rotation) and the documented good ones accepted.
!
! 2. `on` reaches the SCF. Water/6-31G* from the default guess, with and
!    without GDM: the energies agree to 1e-8 (GDM is another route to the same
!    minimum) and only the GDM run logs the handover line and the step lines,
!    so a handle that silently ignored the option cannot pass. The plain run
!    is the negative control for the log reading itself.
!
! 3. `start` reaches the SCF: with a start far below where DIIS converges the
!    handover never happens, so no GDM line is logged; the default start does
!    hand over on the same system.
!
! 4. `max_rotation` reaches the SCF: a rotation cap of 1e-4 rad limits every
!    step, so the run needs many more iterations than with the default cap
!    (measured on the host: 27 against 11) and reaches the same energy.
!
! 5. GDM is a handle option, so it reaches every SCF a Fukui run does on that
!    handle (the N+1 and N-1 ions and the re-run neutral): the handover is
!    logged once per SCF.
!
module check_capi_gdm_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char, c_associated
   use pic_types, only: dp
   use pic_strings, only: to_string
   use pic_logger, only: global_logger, info_level
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_method, &
                               trc_set_convergence, trc_set_gdm, trc_set_verbose, trc_run_scf, trc_energy, &
                               trc_converged, trc_iterations, trc_fukui_hirshfeld, TRC_OK, TRC_ERR_BADARG
   implicit none
   private
   public :: collect_capi_gdm, set_capi_gdm_logdir

   character(len=:), allocatable :: g_logdir
   !! Where the scratch logs go (the build directory, from the command line).

   character(len=*), parameter :: HANDOVER_LINE = "geometric direct minimisation steps from the next one"
   !! Logged once, at the handover iteration, by nothing but GDM.
   character(len=*), parameter :: STEP_LINE = " gdm: "
   !! Starts every GDM step line.

contains

   subroutine set_capi_gdm_logdir(dir)
      character(len=*), intent(in) :: dir
      g_logdir = dir
   end subroutine set_capi_gdm_logdir

   subroutine collect_capi_gdm(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("trc_set_gdm: bad arguments are refused, good ones accepted", test_validation), &
                  new_unittest("trc_set_gdm: on reaches the SCF (same energy, GDM lines only with it)", test_on), &
                  new_unittest("trc_set_gdm: start reaches the SCF (a start below convergence never hands over)", &
                               test_start), &
                  new_unittest("trc_set_gdm: max_rotation reaches the SCF (a tiny cap costs iterations)", &
                               test_max_rotation), &
                  new_unittest("trc_set_gdm: reaches all three SCFs of a Fukui run", test_fukui) &
                  ]
   end subroutine collect_capi_gdm

   subroutine water_handle(h, maxit, error)
      !! A water/6-31G* handle on the default guess, verbose so the SCF logs.
      type(c_ptr), intent(out) :: h
      integer, intent(in) :: maxit
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm
      integer(c_int) :: rc
      h = c_null_ptr
      call read_xyz('water.xyz', natm, z, at_r)
      rc = trc_create(h)
      call check(error, rc == TRC_OK, "trc_create")
      if (allocated(error)) return
      rc = trc_set_molecule(h, int(natm, c_int), real(z, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      if (rc == TRC_OK) rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      if (rc == TRC_OK) rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-6_c_double, int(maxit, c_int), 10_c_int)
      if (rc == TRC_OK) rc = trc_set_verbose(h, 1_c_int)
      call check(error, rc == TRC_OK, "water handle setup (status "//to_string(int(rc))//")")
      if (allocated(error)) then
         rc = trc_destroy(h)
         h = c_null_ptr
      end if
   end subroutine water_handle

   subroutine run_logged(h, tag, status, converged, niter, energy, text)
      !! trc_run_scf with the logger writing to a scratch file; what it logged.
      type(c_ptr), intent(in) :: h
      character(len=*), intent(in) :: tag
      integer(c_int), intent(out) :: status
      logical, intent(out) :: converged
      integer, intent(out) :: niter
      real(dp), intent(out) :: energy
      character(len=:), allocatable, intent(out) :: text
      character(len=:), allocatable :: file
      character(len=512) :: line
      integer(c_int) :: rc, flag, n
      real(c_double) :: e
      integer :: u, ios
      file = g_logdir//"/check_capi_gdm_"//tag//".log"
      open (newunit=u, file=file, status="replace")
      close (u)
      call global_logger%configure_file_output(file, level=info_level)
      status = trc_run_scf(h)
      call global_logger%close_log_file()
      rc = trc_converged(h, flag)
      converged = flag == 1
      rc = trc_iterations(h, n)
      niter = n
      rc = trc_energy(h, e)
      energy = e
      text = ""
      open (newunit=u, file=file, status="old", action="read", iostat=ios)
      if (ios /= 0) return
      do
         read (u, "(a)", iostat=ios) line
         if (ios /= 0) exit
         text = text//trim(line)//new_line("a")
      end do
      close (u)
   end subroutine run_logged

   subroutine test_validation(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: nan
      nan = 0.0_c_double
      nan = nan/nan
      call water_handle(h, 100, error)
      if (allocated(error)) return
      rc = trc_set_gdm(h, 2_c_int, 0.0_c_double, 0.0_c_double)
      call check(error, rc == TRC_ERR_BADARG, "on = 2 must be refused")
      if (.not. allocated(error)) then
         rc = trc_set_gdm(h, 1_c_int, -1.0_c_double, 0.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "a negative start must be refused")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_gdm(h, 1_c_int, nan, 0.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "a NaN start must be refused")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_gdm(h, 1_c_int, 0.0_c_double, -0.1_c_double)
         call check(error, rc == TRC_ERR_BADARG, "a negative max_rotation must be refused")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_gdm(h, 1_c_int, 0.0_c_double, nan)
         call check(error, rc == TRC_ERR_BADARG, "a NaN max_rotation must be refused")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_gdm(h, 0_c_int, 0.0_c_double, 0.0_c_double)
         call check(error, rc == TRC_OK, "(0, 0, 0) must be accepted")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_gdm(h, 1_c_int, 0.1_c_double, 0.2_c_double)
         call check(error, rc == TRC_OK, "(1, 0.1, 0.2) must be accepted")
      end if
      rc = trc_destroy(h)
   end subroutine test_validation

   subroutine test_on(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, status
      logical :: conv_plain, conv_gdm
      integer :: n_plain, n_gdm
      real(dp) :: e_plain, e_gdm
      character(len=:), allocatable :: text_plain, text_gdm
      call water_handle(h, 100, error)
      if (allocated(error)) return
      call run_logged(h, "plain", status, conv_plain, n_plain, e_plain, text_plain)
      rc = trc_destroy(h)
      call check(error, status == TRC_OK .and. conv_plain, "the plain run must converge")
      if (allocated(error)) return
      call check(error, index(text_plain, HANDOVER_LINE) == 0 .and. index(text_plain, STEP_LINE) == 0, &
                 "a plain run logged GDM lines: the log reading is not a control")
      if (allocated(error)) return

      call water_handle(h, 100, error)
      if (allocated(error)) return
      rc = trc_set_gdm(h, 1_c_int, 0.0_c_double, 0.0_c_double)
      call check(error, rc == TRC_OK, "trc_set_gdm")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      call run_logged(h, "gdm", status, conv_gdm, n_gdm, e_gdm, text_gdm)
      rc = trc_destroy(h)
      call check(error, status == TRC_OK .and. conv_gdm, "the GDM run must converge (status "// &
                 to_string(int(status))//")")
      if (allocated(error)) return
      call check(error, index(text_gdm, HANDOVER_LINE) > 0, "trc_set_gdm(on) did not reach the SCF: no GDM handover logged")
      if (allocated(error)) return
      call check(error, index(text_gdm, STEP_LINE) > 0, "no GDM step logged")
      if (allocated(error)) return
      call check(error, abs(e_gdm - e_plain) < 1.0e-8_dp, "GDM energy differs from the plain run: "// &
                 to_string(e_gdm - e_plain, "es10.2"))
   end subroutine test_on

   subroutine test_start(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, status
      logical :: conv
      integer :: n
      real(dp) :: e
      character(len=:), allocatable :: text_low, text_default
      call water_handle(h, 100, error)
      if (allocated(error)) return
      ! DIIS converges at a commutator of 1e-6 (conv_diis of this handle), long
      ! before it reaches 1e-12.
      rc = trc_set_gdm(h, 1_c_int, 1.0e-12_c_double, 0.0_c_double)
      call run_logged(h, "start_low", status, conv, n, e, text_low)
      rc = trc_destroy(h)
      call check(error, status == TRC_OK .and. conv, "the run with a low start must converge (as DIIS)")
      if (allocated(error)) return
      call check(error, index(text_low, HANDOVER_LINE) == 0, "start = 1e-12 handed over anyway: start was not honoured")
      if (allocated(error)) return

      call water_handle(h, 100, error)
      if (allocated(error)) return
      rc = trc_set_gdm(h, 1_c_int, 0.0_c_double, 0.0_c_double)
      call run_logged(h, "start_default", status, conv, n, e, text_default)
      rc = trc_destroy(h)
      call check(error, index(text_default, HANDOVER_LINE) > 0, "the default start did not hand over on the same system")
   end subroutine test_start

   subroutine test_max_rotation(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, status
      logical :: conv_default, conv_tiny
      integer :: n_default, n_tiny
      real(dp) :: e_default, e_tiny
      character(len=:), allocatable :: text
      integer, parameter :: EXTRA = 5
      !! Iterations the tiny cap must cost at least; measured 16.
      call water_handle(h, 100, error)
      if (allocated(error)) return
      rc = trc_set_gdm(h, 1_c_int, 0.0_c_double, 0.0_c_double)
      call run_logged(h, "rot_default", status, conv_default, n_default, e_default, text)
      rc = trc_destroy(h)
      call check(error, conv_default, "the default cap must converge")
      if (allocated(error)) return

      call water_handle(h, 100, error)
      if (allocated(error)) return
      rc = trc_set_gdm(h, 1_c_int, 0.0_c_double, 1.0e-4_c_double)
      call run_logged(h, "rot_tiny", status, conv_tiny, n_tiny, e_tiny, text)
      rc = trc_destroy(h)
      call check(error, conv_tiny, "max_rotation = 1e-4 must still converge")
      if (allocated(error)) return
      call check(error, n_tiny >= n_default + EXTRA, "max_rotation = 1e-4 took "//to_string(n_tiny)// &
                 " iterations against "//to_string(n_default)//" with the default: the cap was not honoured")
      if (allocated(error)) return
      call check(error, abs(e_tiny - e_default) < 1.0e-8_dp, "the capped run reached another energy: "// &
                 to_string(e_tiny - e_default, "es10.2"))
   end subroutine test_max_rotation

   integer function count_of(text, pattern) result(n)
      !! How many times `pattern` occurs in `text`.
      character(len=*), intent(in) :: text, pattern
      integer :: k, pos
      n = 0
      pos = 1
      do
         k = index(text(pos:), pattern)
         if (k == 0) exit
         n = n + 1
         pos = pos + k - 1 + len(pattern)
      end do
   end function count_of

   subroutine test_fukui(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, status, st
      logical :: conv
      integer :: n
      real(dp) :: e
      real(c_double) :: fp(3), fm(3), fz(3), q0(3), qp(3), qm(3)
      character(len=:), allocatable :: text, file
      character(len=512) :: line
      integer :: u, ios
      call water_handle(h, 100, error)
      if (allocated(error)) return
      rc = trc_set_gdm(h, 1_c_int, 0.0_c_double, 0.0_c_double)
      call check(error, rc == TRC_OK, "trc_set_gdm")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      ! The neutral SCF first, as the entry point requires; its handover is
      ! logged apart so the count below is the Fukui SCFs' alone.
      call run_logged(h, "fukui_ref", status, conv, n, e, text)
      call check(error, status == TRC_OK .and. count_of(text, HANDOVER_LINE) == 1, &
                 "the reference SCF should hand over once")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      file = g_logdir//"/check_capi_gdm_fukui.log"
      open (newunit=u, file=file, status="replace")
      close (u)
      call global_logger%configure_file_output(file, level=info_level)
      st = trc_fukui_hirshfeld(h, fp, fm, fz, q0, qp, qm)
      call global_logger%close_log_file()
      rc = trc_destroy(h)
      call check(error, st == TRC_OK, "trc_fukui_hirshfeld (status "//to_string(int(st))//")")
      if (allocated(error)) return
      text = ""
      open (newunit=u, file=file, status="old", action="read", iostat=ios)
      if (ios == 0) then
         do
            read (u, "(a)", iostat=ios) line
            if (ios /= 0) exit
            text = text//trim(line)//new_line("a")
         end do
         close (u)
      end if
      call check(error, count_of(text, HANDOVER_LINE) == 3, "GDM should hand over in each of the three Fukui SCFs, "// &
                 "counted "//to_string(count_of(text, HANDOVER_LINE)))
   end subroutine test_fukui

end module check_capi_gdm_suite

program check_capi_gdm
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_capi_gdm_suite, only: collect_capi_gdm, set_capi_gdm_logdir
   implicit none
   character(len=1024) :: dir
   integer :: stat
   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) dir = "."
   call set_capi_gdm_logdir(trim(dir))
   stat = 0
   call run_testsuite(collect_capi_gdm, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_gdm
