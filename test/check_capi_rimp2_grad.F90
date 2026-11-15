module check_capi_rimp2_grad_suite
   !!
   !! The C ABI's RI-MP2 energy when a gradient follows: trc_run_rimp2 then
   !! trc_run_gradient must not run the separate energy step (trc_run_rimp2
   !! defers the pair energies on one rank, and the gradient fills them from
   !! its own amplitudes), and what trc_rimp2_energy and trc_energy then report
   !! must be what an energy-only handle reports -- trc_rimp2_run's E_os and
   !! E_ss -- to 1e-12. Three references: the exact four-centre SCF, a --df SCF
   !! over def2-universal-jkfit (the energy's own B fitted separately), and a
   !! --df SCF whose fit is RI-MP2's own auxiliary set (the energy chained off
   !! the SCF's B, trc_capi.F90's packed_b path), all water/cc-pVDZ/cc-pVDZ-
   !! RIFIT with the elements' frozen core.
   !!
   !! "Did not run" is read off the verbose log: the energy step's own table
   !! (trc_rimp2_run's "rimp2:" lines) is in an energy-only handle's log and
   !! absent from the gradient handle's, which has the gradient's own
   !! ("rimp2 grad:") instead.
   !!
   !! Argument: a scratch directory for the logs (default: the working one).
   !! A test-drive suite.
   !!
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char
   use pic_types, only: dp
   use pic_logger, only: global_logger, verbose_level
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_jkfit_json, trc_set_df, trc_set_convergence, trc_set_rimp2, &
                               trc_set_verbose, trc_run_scf, trc_run_rimp2, trc_run_gradient, trc_rimp2_energy, &
                               trc_energy, TRC_OK
   implicit none
   private
   public :: collect_capi_rimp2_grad, set_capi_rimp2_grad_logdir

   real(dp), parameter :: ENERGY_TOL = 1.0e-12_dp
   !! The gradient's amplitudes against trc_rimp2_run's pair energies: the
   !! same B and pairs, so round-off apart.

   character(len=:), allocatable :: g_logdir

contains

   subroutine set_capi_rimp2_grad_logdir(dir)
      character(*), intent(in) :: dir
      g_logdir = dir
   end subroutine set_capi_rimp2_grad_logdir

   subroutine collect_capi_rimp2_grad(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("exact reference: the gradient's energy, no energy step", test_exact), &
                  new_unittest("--df (universal-jkfit): the gradient's energy, no energy step", test_df), &
                  new_unittest("--df chained (jkfit = aux): the gradient's energy, no energy step", test_df_chained) &
                  ]
   end subroutine collect_capi_rimp2_grad

   subroutine test_exact(error)
      type(error_type), allocatable, intent(out) :: error
      call compare(error, "exact", "")
   end subroutine test_exact

   subroutine test_df(error)
      type(error_type), allocatable, intent(out) :: error
      call compare(error, "df", "basis_sets/def2-universal-jkfit.json")
   end subroutine test_df

   subroutine test_df_chained(error)
      type(error_type), allocatable, intent(out) :: error
      call compare(error, "df_chained", "basis_sets/cc-pvdz-rifit.json")
   end subroutine test_df_chained

   subroutine compare(error, tag, jkfit)
      !! Two handles on the same input: one reads the RI-MP2 energy alone,
      !! the other runs the gradient first and reads it after.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: tag, jkfit
      real(c_double) :: os_ref, ss_ref, e_ref, os_grad, ss_grad, e_grad
      logical :: energy_table_ref, energy_table_grad, grad_table

      call run_handle(error, tag//": energy only", jkfit, .false., &
                      g_logdir//"/check_capi_rimp2_grad_"//tag//"_energy.log", os_ref, ss_ref, e_ref, &
                      energy_table_ref, grad_table)
      if (allocated(error)) return
      call check(error, energy_table_ref, tag//": the energy-only log has no energy-step table (the probe is blind)")
      if (allocated(error)) return
      call run_handle(error, tag//": gradient first", jkfit, .true., &
                      g_logdir//"/check_capi_rimp2_grad_"//tag//"_grad.log", os_grad, ss_grad, e_grad, &
                      energy_table_grad, grad_table)
      if (allocated(error)) return
      print "(a,a,es10.2,a,es10.2,a,es10.2)", "  ", tag//": |dE_os| ", abs(os_grad - os_ref), "  |dE_ss| ", &
         abs(ss_grad - ss_ref), "  |dE_total| ", abs(e_grad - e_ref)
      call check(error, grad_table, tag//": the gradient's own step table is missing from its log")
      if (allocated(error)) return
      call check(error, .not. energy_table_grad, tag//": the separate RI-MP2 energy step ran before the gradient")
      if (allocated(error)) return
      call check(error, abs(os_grad - os_ref) <= ENERGY_TOL .and. abs(ss_grad - ss_ref) <= ENERGY_TOL, &
                 tag//": E_os/E_ss from the gradient differ from the energy step's")
      if (allocated(error)) return
      call check(error, abs(e_grad - e_ref) <= ENERGY_TOL, tag//": trc_energy differs after the gradient")
   end subroutine compare

   subroutine run_handle(error, what, jkfit, gradient_first, logname, e_os, e_ss, e_total, energy_table, grad_table)
      !! water/cc-pVDZ RHF (+ --df over `jkfit` when it is not empty), then
      !! trc_run_rimp2, [trc_run_gradient,] trc_rimp2_energy, trc_energy,
      !! with the verbose log captured in `logname`.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: what, jkfit, logname
      logical, intent(in) :: gradient_first
      real(c_double), intent(out) :: e_os, e_ss, e_total
      logical, intent(out) :: energy_table, grad_table
      type(c_ptr) :: h
      integer :: natm, u, ios
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      character(len=512) :: line
      integer(c_int) :: rc

      e_os = 0.0_c_double
      e_ss = 0.0_c_double
      e_total = 0.0_c_double
      energy_table = .false.
      grad_table = .false.
      call read_xyz('water_mqc.xyz', natm, zint, at_r)
      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                                              0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      if (rc == TRC_OK) rc = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      if (rc == TRC_OK .and. len(jkfit) > 0) rc = trc_set_jkfit_json(h, jkfit//c_null_char)
      if (rc == TRC_OK .and. len(jkfit) > 0) rc = trc_set_df(h, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 100_c_int, 10_c_int)
      if (rc == TRC_OK) rc = trc_set_rimp2(h, -1_c_int, 0_c_int)
      call check(error, rc == TRC_OK, what//": setup failed")
      if (allocated(error)) return
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, what//": SCF failed")
      if (allocated(error)) return
      ! Only the RI-MP2 steps are logged, and verbose so their tables print.
      call global_logger%configure_file_output(logname, level=verbose_level)
      rc = trc_set_verbose(h, 1_c_int)
      rc = trc_run_rimp2(h)
      if (rc == TRC_OK .and. gradient_first) rc = trc_run_gradient(h)
      if (rc == TRC_OK) rc = trc_rimp2_energy(h, e_os, e_ss)
      if (rc == TRC_OK) rc = trc_energy(h, e_total)
      call global_logger%close_log_file()
      call check(error, rc == TRC_OK, what//": an RI-MP2 entry failed")
      if (allocated(error)) return
      rc = trc_destroy(h)
      open (newunit=u, file=logname, status="old", action="read", iostat=ios)
      call check(error, ios == 0, what//": cannot read "//logname)
      if (allocated(error)) return
      do
         read (u, "(a)", iostat=ios) line
         if (ios /= 0) exit
         if (index(line, " rimp2:") > 0) energy_table = .true.
         if (index(line, " rimp2 grad:") > 0) grad_table = .true.
      end do
      close (u)
   end subroutine run_handle

end module check_capi_rimp2_grad_suite

program check_capi_rimp2_grad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_capi_rimp2_grad_suite, only: collect_capi_rimp2_grad, set_capi_rimp2_grad_logdir
   implicit none
   character(len=1024) :: dir
   integer :: stat
   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) dir = "."
   call set_capi_rimp2_grad_logdir(trim(dir))
   stat = 0
   call run_testsuite(collect_capi_rimp2_grad, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_rimp2_grad
