!
! --opt with a correlated method, through the C ABI only (terco.h for the
! context, terco_opt.h for the optimisation): RI-MP2 over RHF, and the
! double hybrid B2PLYP (E_KS + c_os E_os + c_ss E_ss). trc_opt_set_rimp2
! makes every step rerun trc_run_rimp2 after its SCF, so the energy DL-FIND
! sees is trc_energy's correlated total and the gradient its derivative.
! Before it, the CLI refused both outright, and a step's energy could only
! ever have been the SCF's -- 0.065 Eh (c E_PT2) above the B2PLYP total here,
! paired with a correlated gradient.
!
! For each method, water from water.xyz (|g| ~ 1e-2 there), def2-SVP with
! def2-SVP-RIFIT, frozen 1:
!   - the optimisation converges, max |g| below the tolerance;
!   - the energy it reports for the FIRST step (the start geometry, far
!     from stationary) equals a fresh context's single-point total there;
!   - its final energy equals a fresh context's single-point total at the
!     final geometry, and that context's analytic gradient is below the
!     tolerance too.
! The single points share nothing with the optimisation (a new handle, a
! SAD guess), so agreement to 1e-9 says the optimiser's energies are the
! correlated totals at those geometries. trc_energy's totals are anchored
! to pyscf by test/scf_dh.F90 and the RI-MP2 energy checks.
!
! A test-drive suite.
!
module check_opt_correlated_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char
   use, intrinsic :: iso_fortran_env, only: output_unit
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_method, trc_set_convergence, trc_set_rimp2, trc_run_scf, trc_run_rimp2, &
                               trc_run_gradient, trc_gradient, trc_energy, trc_message, trc_set_pcm, &
                               trc_iterations, TRC_OK, TRC_ERR_UNSUPPORTED, TRC_ERR_STATE
   use trc_opt_c_interfaces, only: trc_opt_options_create, trc_opt_options_destroy, trc_opt_set_gtol, &
                                   trc_opt_set_rimp2, trc_opt_result_create, trc_opt_result_destroy, trc_opt_run, &
                                   trc_opt_result_converged, trc_opt_result_energy, trc_opt_result_gradient_max, &
                                   trc_opt_result_geometry, trc_opt_result_trajectory_count, &
                                   trc_opt_result_trajectory, trc_opt_result_message
   implicit none
   private
   public :: collect_opt_correlated

   real(dp), parameter :: GTOL = 1.0e-5_dp
   !! Eh/Bohr, max Cartesian component: tighter than the CLI's default
   !! 4.5e-4, so the final-gradient check means something.
   real(dp), parameter :: ENERGY_TOL = 1.0e-9_dp
   !! Measured: see the "diff" lines each case prints (worst 3.3e-13, gfortran
   !! 15 host build) -- the SCFs converge to 1e-12 Eh / 1e-9 commutator.

contains

   subroutine collect_opt_correlated(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("--opt b2plyp/def2-SVP: the double-hybrid total at every step", test_dh), &
                  new_unittest("--opt RI-MP2/def2-SVP: the RHF + RI-MP2 total at every step", test_rimp2), &
                  new_unittest("trc_opt_run refuses an unrestricted RI-MP2 setup before its first step", &
                               test_refuses_unrestricted), &
                  new_unittest("trc_opt_run refuses a double hybrid over PCM before its first step", &
                               test_refuses_pcm), &
                  new_unittest("trc_opt_run refuses a correlated run with no auxiliary basis", test_refuses_no_aux) &
                  ]
   end subroutine collect_opt_correlated

   subroutine test_dh(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "b2plyp")
   end subroutine test_dh

   subroutine test_rimp2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "")
   end subroutine test_rimp2

   subroutine test_refuses_unrestricted(error)
      !! The water cation, a doublet: UHF + RI-MP2 has an energy but no
      !! gradient yet.
      type(error_type), allocatable, intent(out) :: error
      call expect_refused(error, "", 1, 2, .false., .true., TRC_ERR_UNSUPPORTED, &
                          "restricted (closed-shell RHF reference) only")
   end subroutine test_refuses_unrestricted

   subroutine test_refuses_pcm(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_refused(error, "b2plyp", 0, 1, .true., .true., TRC_ERR_UNSUPPORTED, &
                          "no PCM response in the Z-vector")
   end subroutine test_refuses_pcm

   subroutine test_refuses_no_aux(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_refused(error, "", 0, 1, .false., .false., TRC_ERR_STATE, "no auxiliary basis was set")
   end subroutine test_refuses_no_aux

   subroutine expect_refused(error, functional, charge, multiplicity, pcm, aux, expected_status, expected)
      !! trc_opt_run with trc_opt_set_rimp2 on a setup the RI-MP2 gradient
      !! refuses: it must RETURN (DL-FIND's `error stop` on the first step
      !! would end this test program instead), with the refusal's status
      !! and message, and without having run a single SCF.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: functional, expected
      integer, intent(in) :: charge, multiplicity
      logical, intent(in) :: pcm, aux
      integer(c_int), intent(in) :: expected_status
      type(c_ptr) :: h, opts, res
      integer :: natm
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      integer(c_int) :: rc, st, niter
      character(kind=c_char) :: msg(1024)
      character(len=:), allocatable :: text
      call read_xyz("water.xyz", natm, z, r)
      call water_context(error, h, functional, r, charge, multiplicity, pcm, aux)
      if (allocated(error)) return
      rc = trc_opt_options_create(opts)
      if (rc == TRC_OK) rc = trc_opt_set_rimp2(opts, 1_c_int)
      if (rc == TRC_OK) rc = trc_opt_result_create(res)
      call check(error, rc == TRC_OK, "trc_opt options/result setup failed")
      if (allocated(error)) return
      rc = trc_opt_run(h, 3_c_int, int(z, c_int), reshape(real(r, c_double), [9]), opts, res)
      st = trc_opt_result_message(res, msg, int(size(msg), c_int))
      text = f_message(msg)
      write (output_unit, "(a,i0,a,a)") "  trc_opt_run status ", rc, ": ", text
      ! No SCF ran: trc_iterations refuses before the first one.
      st = trc_iterations(h, niter)
      rc = merge(rc, -999_c_int, st /= TRC_OK)
      st = trc_opt_result_destroy(res)
      st = trc_opt_options_destroy(opts)
      st = trc_destroy(h)
      call check(error, rc == expected_status .and. index(text, expected) > 0, &
                 "trc_opt_run should refuse before its first SCF with '"//expected//"'", &
                 more="status "//real_text(real(rc, dp))//": "//text)
   end subroutine expect_refused

   subroutine water_context(error, h, functional, xyz, charge, multiplicity, pcm, aux)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(out) :: h
      character(len=*), intent(in) :: functional
      real(dp), intent(in) :: xyz(3, 3)
      integer, intent(in), optional :: charge, multiplicity
      logical, intent(in), optional :: pcm, aux
      integer(c_int) :: rc, chg, mult
      logical :: with_pcm, with_aux
      chg = 0_c_int
      mult = 1_c_int
      if (present(charge)) chg = int(charge, c_int)
      if (present(multiplicity)) mult = int(multiplicity, c_int)
      with_pcm = .false.
      if (present(pcm)) with_pcm = pcm
      with_aux = .true.
      if (present(aux)) with_aux = aux
      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, 3_c_int, [8.0_c_double, 1.0_c_double, 1.0_c_double], &
                                              reshape(real(xyz, c_double), [9]), chg, mult)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, "basis_sets/def2-svp.json"//c_null_char)
      if (rc == TRC_OK .and. with_aux) rc = trc_set_aux_json(h, "basis_sets/def2-svp-rifit.json"//c_null_char)
      if (rc == TRC_OK .and. with_pcm) rc = trc_set_pcm(h, "cpcm"//c_null_char, 78.3553_c_double, 0_c_int, &
                                                         0.0_c_double)
      if (rc == TRC_OK) rc = trc_set_method(h, functional//c_null_char, 3_c_int)
      if (rc == TRC_OK) rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 200_c_int, 8_c_int)
      if (rc == TRC_OK) rc = trc_set_rimp2(h, 1_c_int, 0_c_int)
      call check(error, rc == TRC_OK, "context setup failed: "//message_of(h))
   end subroutine water_context

   subroutine single_point(error, functional, xyz, energy, gmax)
      !! A fresh handle's SCF + RI-MP2 total and max |gradient| at `xyz`.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: functional
      real(dp), intent(in) :: xyz(3, 3)
      real(dp), intent(out) :: energy, gmax
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: e, gbuf(9)
      energy = 0.0_dp
      gmax = huge(1.0_dp)
      call water_context(error, h, functional, xyz)
      if (allocated(error)) return
      rc = trc_run_scf(h)
      if (rc == TRC_OK) rc = trc_run_rimp2(h)
      if (rc == TRC_OK) rc = trc_run_gradient(h)
      if (rc == TRC_OK) rc = trc_gradient(h, gbuf)
      if (rc == TRC_OK) rc = trc_energy(h, e)
      call check(error, rc == TRC_OK, "single point failed: "//message_of(h))
      if (.not. allocated(error)) then
         energy = real(e, dp)
         gmax = maxval(abs(real(gbuf, dp)))
      end if
      rc = trc_destroy(h)
   end subroutine single_point

   subroutine run_case(error, functional)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: functional
      type(c_ptr) :: h, opts, res
      integer :: natm
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: x0(3, 3), xfinal(3, 3), e_start, e_final, g_start, g_final
      integer(c_int) :: rc, converged, n_traj
      real(c_double) :: e_opt, gmax_opt, xyz_out(9)
      real(c_double), allocatable :: traj_xyz(:), traj_e(:)
      character(kind=c_char) :: msg(512)

      call read_xyz("water.xyz", natm, z, r)
      x0 = r
      call water_context(error, h, functional, x0)
      if (allocated(error)) return
      rc = trc_opt_options_create(opts)
      if (rc == TRC_OK) rc = trc_opt_set_gtol(opts, real(GTOL, c_double))
      if (rc == TRC_OK) rc = trc_opt_set_rimp2(opts, 1_c_int)
      if (rc == TRC_OK) rc = trc_opt_result_create(res)
      call check(error, rc == TRC_OK, "trc_opt options/result setup failed")
      if (allocated(error)) return
      rc = trc_opt_run(h, 3_c_int, int(z, c_int), reshape(real(x0, c_double), [9]), opts, res)
      if (rc /= TRC_OK) then
         rc = trc_opt_result_message(res, msg, int(size(msg), c_int))
         call check(error, .false., "trc_opt_run failed: "//f_message(msg))
         return
      end if
      rc = trc_opt_result_converged(res, converged)
      rc = trc_opt_result_energy(res, e_opt)
      rc = trc_opt_result_gradient_max(res, gmax_opt)
      rc = trc_opt_result_geometry(res, 3_c_int, xyz_out)
      rc = trc_opt_result_trajectory_count(res, n_traj)
      allocate (traj_xyz(9*max(n_traj, 1)), traj_e(max(n_traj, 1)))
      if (n_traj > 0) rc = trc_opt_result_trajectory(res, 3_c_int, n_traj, traj_xyz, traj_e)
      rc = trc_opt_result_destroy(res)
      rc = trc_opt_options_destroy(opts)
      rc = trc_destroy(h)
      xfinal = reshape(real(xyz_out, dp), [3, 3])

      call check(error, converged /= 0 .and. gmax_opt <= GTOL, "the optimisation did not converge below the "// &
                 "gradient tolerance", more="max |g| "//real_text(real(gmax_opt, dp)))
      if (allocated(error)) return
      call check(error, n_traj >= 2, "expected at least two accepted steps from water.xyz")
      if (allocated(error)) return

      ! The first accepted step is the start geometry itself: far from
      ! stationary, so an SCF energy paired with a correlated gradient
      ! cannot hide there.
      call check(error, maxval(abs(traj_xyz(1:9) - reshape(real(x0, c_double), [9]))) <= 1.0e-12_c_double, &
                 "the trajectory's first frame is not the start geometry")
      if (allocated(error)) return
      call single_point(error, functional, x0, e_start, g_start)
      if (allocated(error)) return
      write (output_unit, "(a,es22.14,a,es22.14,a,es10.3,a,es10.3)") "  start: optimiser ", traj_e(1), &
         "  single point ", e_start, "  diff ", abs(traj_e(1) - e_start), "  max|g| ", g_start
      call check(error, abs(real(traj_e(1), dp) - e_start) <= ENERGY_TOL, "the optimiser's energy at the start "// &
                 "geometry is not the correlated total there", more="optimiser "//real_text(real(traj_e(1), dp))// &
                 ", single point "//real_text(e_start))
      if (allocated(error)) return

      call single_point(error, functional, xfinal, e_final, g_final)
      if (allocated(error)) return
      write (output_unit, "(a,es22.14,a,es22.14,a,es10.3,a,es10.3)") "  final: optimiser ", e_opt, &
         "  single point ", e_final, "  diff ", abs(e_opt - e_final), "  max|g| ", g_final
      call check(error, abs(real(e_opt, dp) - e_final) <= ENERGY_TOL, "the optimiser's final energy is not the "// &
                 "correlated total at the final geometry", more="optimiser "//real_text(real(e_opt, dp))// &
                 ", single point "//real_text(e_final))
      if (allocated(error)) return
      call check(error, g_final <= GTOL, "the single-point gradient at the final geometry is above the tolerance", &
                 more=real_text(g_final))
      if (allocated(error)) return
      call check(error, e_final < e_start, "the optimisation did not lower the correlated energy")
   end subroutine run_case

   function message_of(h) result(text)
      type(c_ptr), intent(in) :: h
      character(len=:), allocatable :: text
      character(kind=c_char, len=1024) :: buffer
      integer(c_int) :: rc
      integer :: k
      buffer = ""
      rc = trc_message(h, buffer, 1024_c_int)
      k = index(buffer, c_null_char)
      if (k == 0) k = len_trim(buffer) + 1
      text = buffer(1:k - 1)
   end function message_of

   function f_message(buf) result(s)
      character(kind=c_char), intent(in) :: buf(:)
      character(len=:), allocatable :: s
      integer :: i, n
      n = 0
      do i = 1, size(buf)
         if (buf(i) == c_null_char) exit
         n = i
      end do
      allocate (character(len=n) :: s)
      do i = 1, n
         s(i:i) = buf(i)
      end do
   end function f_message

   function real_text(x) result(text)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: text
      character(len=32) :: buffer
      write (buffer, "(es22.14)") x
      text = trim(adjustl(buffer))
   end function real_text

end module check_opt_correlated_suite

program check_opt_correlated
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_correlated_suite, only: collect_opt_correlated
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_opt_correlated, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_correlated
