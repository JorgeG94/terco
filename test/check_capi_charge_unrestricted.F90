module check_capi_charge_unrestricted_suite
   !! trc_set_unrestricted, trc_set_level_shift and trc_set_charge,
   !! through the context API.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_df, only: df_jk_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_charge, &
                               trc_set_basis_json, trc_set_jkfit_json, trc_set_df, &
                               trc_set_method, trc_set_convergence, trc_set_unrestricted, &
                               trc_set_level_shift, trc_run_scf, trc_energy, trc_nspin, &
                               TRC_OK, TRC_ERR_BADARG
   implicit none
   private
   public :: collect_capi_charge_unrestricted_tests

contains

   subroutine collect_capi_charge_unrestricted_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("forced UHF on a closed-shell singlet matches RHF", &
                               test_unrestricted_matches_rhf), &
                  new_unittest("the level shift moves the path, not the converged energy", &
                               test_level_shift), &
                  new_unittest("df_jk_t%setup is idempotent on a second call", &
                               test_df_setup_is_idempotent), &
                  new_unittest("trc_set_charge reuses the DF fit and matches fresh handles", &
                               test_charge_reuses_df_fit) &
                  ]
   end subroutine collect_capi_charge_unrestricted_tests

   subroutine water_geometry(natm, zint, at_r)
      integer, intent(out) :: natm
      integer, allocatable, intent(out) :: zint(:)
      real(dp), allocatable, intent(out) :: at_r(:, :)
      call read_xyz('water.xyz', natm, zint, at_r)
   end subroutine water_geometry

   !
   ! trc_set_unrestricted(1) on a closed-shell singlet (nalpha == nbeta)
   ! must converge UHF to the SAME stationary point as the default
   ! restricted run: forcing two independent spin densities on a genuine
   ! closed shell does not change the fixed point, only whether it is
   ! reached by one eigenproblem per iteration or two identical ones.
   !
   subroutine test_unrestricted_matches_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(c_ptr) :: h_rhf, h_uhf
      integer(c_int) :: rc, nspin_rhf, nspin_uhf
      real(c_double) :: e_rhf, e_uhf

      call water_geometry(natm, zint, at_r)

      rc = trc_create(h_rhf)
      rc = trc_set_molecule(h_rhf, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                            0_c_int, 1_c_int)
      rc = trc_set_basis_json(h_rhf, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_method(h_rhf, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h_rhf, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(h_rhf)
      call check(error, rc == TRC_OK, "unrestricted: run_scf (RHF) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h_rhf)
         return
      end if
      rc = trc_nspin(h_rhf, nspin_rhf)
      call check(error, nspin_rhf == 1, "unrestricted: RHF's default is nspin 1")
      if (allocated(error)) then
         rc = trc_destroy(h_rhf)
         return
      end if
      rc = trc_energy(h_rhf, e_rhf)

      rc = trc_create(h_uhf)
      rc = trc_set_molecule(h_uhf, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                            0_c_int, 1_c_int)
      rc = trc_set_basis_json(h_uhf, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_unrestricted(h_uhf, 1_c_int)
      call check(error, rc == TRC_OK, "unrestricted: set_unrestricted(1) status")
      if (allocated(error)) then
         rc = trc_destroy(h_rhf); rc = trc_destroy(h_uhf)
         return
      end if
      rc = trc_set_method(h_uhf, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h_uhf, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(h_uhf)
      call check(error, rc == TRC_OK, "unrestricted: run_scf (forced UHF) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h_rhf); rc = trc_destroy(h_uhf)
         return
      end if
      rc = trc_nspin(h_uhf, nspin_uhf)
      call check(error, nspin_uhf == 2, "unrestricted: forced UHF is nspin 2")
      if (allocated(error)) then
         rc = trc_destroy(h_rhf); rc = trc_destroy(h_uhf)
         return
      end if
      rc = trc_energy(h_uhf, e_uhf)

      call check(error, abs(real(e_rhf, dp) - real(e_uhf, dp)) < 1.0e-10_dp, &
                 "unrestricted: forced UHF on a closed shell matches RHF, E_RHF "// &
                 to_string(real(e_rhf, dp), "f20.12")//"  E_UHF(forced) "// &
                 to_string(real(e_uhf, dp), "f20.12")//"  diff "// &
                 to_string(real(e_rhf, dp) - real(e_uhf, dp), "es9.2"))

      rc = trc_destroy(h_rhf); rc = trc_destroy(h_uhf)
   end subroutine test_unrestricted_matches_rhf

   !
   ! trc_set_level_shift only changes the early, pre-DIIS part of the
   ! iteration (trc_scf_driver.F90's own "while DIIS and A/EDIIS are
   ! off" gate), so two different positive shifts must converge to the
   ! SAME energy -- the level shift moves the path, not the fixed point.
   ! A negative shift is refused; 0 (off) must still converge.
   !
   subroutine test_level_shift(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: e_default, e_half, e_off

      call water_geometry(natm, zint, at_r)

      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)

      rc = trc_set_level_shift(h, -1.0_c_double)
      call check(error, rc == TRC_ERR_BADARG, "level_shift: a negative shift is refused")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      rc = trc_set_level_shift(h, 0.5_c_double)
      call check(error, rc == TRC_OK, "level_shift: set 0.5 status")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "level_shift: run_scf (0.5) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_energy(h, e_half)

      rc = trc_set_level_shift(h, 0.0_c_double)
      call check(error, rc == TRC_OK, "level_shift: set 0 (off) status")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "level_shift: run_scf (off) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_energy(h, e_off)

      rc = trc_set_level_shift(h, 1.0_c_double)
      call check(error, rc == TRC_OK, "level_shift: set 1.0 (default) status")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "level_shift: run_scf (1.0) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_energy(h, e_default)

      call check(error, abs(real(e_half, dp) - real(e_default, dp)) < 1.0e-9_dp, &
                 "level_shift: 0.5 and 1.0 converge to the same energy")
      if (.not. allocated(error)) then
         call check(error, abs(real(e_off, dp) - real(e_default, dp)) < 1.0e-9_dp, &
                    "level_shift: off and 1.0 converge to the same energy")
      end if

      rc = trc_destroy(h)
   end subroutine test_level_shift

   !
   ! The mechanism trc_set_charge relies on, proven directly rather than
   ! inferred from wall time at C-API scale: df_jk_t%setup is idempotent
   ! (trc_scf_df.F90's own header, "a second call on an object that
   ! already built its bpk returns immediately") -- a SECOND call on the
   ! SAME df_jk_t must take negligible time against the first, which
   ! does the metric/Cholesky/3c/trsm for real. No margin to tune: the
   ! second call's own code path is `if (this%ready) return`, nothing
   ! else, so it costs the same handful of instructions at any problem
   ! size.
   !
   subroutine test_df_setup_is_idempotent(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: la
      type(df_jk_t) :: dfjk
      type(error_t) :: err
      type(timer_type) :: tm
      real(dp) :: t_first, t_second
      real(dp), parameter :: ERI_THRESH = 1.0e-12_dp

      call water_geometry(natm, zint, at_r)

      call trc_basis_from_json('basis_sets/6-31gs.json', natm, zint, at_r, bas, err, force_cart=.true.)
      if (.not. err%has_error()) call trc_basis_from_json('basis_sets/def2-universal-jkfit.json', &
                                                           natm, zint, at_r, aux, err)
      call check(error, .not. err%has_error(), "df_setup: basis construction: "//err%get_message())
      if (allocated(error)) return

      call bas%to_device(); call aux%to_device()
      call pl%build(bas, ERI_THRESH); call pl%to_device()
      call la%init(bas%nao)

      call tm%start()
      call dfjk%setup(bas, aux, pl, la, err)
      t_first = tm%get_elapsed_time()
      call check(error, .not. err%has_error(), "df_setup: first setup succeeds")
      if (.not. allocated(error)) call check(error, dfjk%is_ready(), "df_setup: ready after the first setup")
      if (allocated(error)) then
         call dfjk%release(); call pl%release(); call la%release()
         call bas%release(); call aux%release()
         return
      end if

      call tm%start()
      call dfjk%setup(bas, aux, pl, la, err)
      t_second = tm%get_elapsed_time()
      call check(error, .not. err%has_error(), "df_setup: second setup (idempotent) succeeds")
      if (.not. allocated(error)) then
         call check(error, t_second < 0.1_dp*max(t_first, 1.0e-6_dp), &
                    "df_setup: a second call on the same object does no work "// &
                    "(trc_set_charge relies on this), first "//to_string(t_first, "es10.2")// &
                    " s, second "//to_string(t_second, "es10.2")//" s")
      end if

      call dfjk%release(); call pl%release(); call la%release()
      call bas%release(); call aux%release()
   end subroutine test_df_setup_is_idempotent

   !
   ! trc_set_charge on one --df handle, through an N, N+2, N-2 sequence
   ! (water: 10, 8, 12 electrons, all closed-shell singlets -- --df is
   ! restricted only, trc_run_gradient's own header says so and
   ! trc_run_scf refuses otherwise, so this stays closed-shell rather
   ! than exercising the N+1/N-1 radical case test_unrestricted_matches_rhf
   ! covers separately above), must match fresh handles built directly
   ! at each charge -- and the fit must NOT be rebuilt in between:
   ! df_jk_setup returns immediately once already `ready`
   ! (trc_scf_df.F90's own header), so trc_run_scf on the second and
   ! third charge states should be far faster than the first, which
   ! pays the metric/Cholesky/3c/trsm. test_df_setup_is_idempotent (above)
   ! is the deterministic proof of that mechanism; this checks the
   ! energies it protects actually match a fresh handle.
   !
   subroutine test_charge_reuses_df_fit(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(c_ptr) :: h, hf
      integer(c_int) :: rc
      real(c_double) :: e_plus2, e_minus2, e_plus2_fresh, e_minus2_fresh

      call water_geometry(natm, zint, at_r)

      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_jkfit_json(h, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
      rc = trc_set_df(h, 1_c_int)
      rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)

      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "charge: run_scf (neutral) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      rc = trc_set_charge(h, 2_c_int, 1_c_int)
      call check(error, rc == TRC_OK, "charge: set_charge(+2) status")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "charge: run_scf (+2) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_energy(h, e_plus2)

      rc = trc_set_charge(h, -2_c_int, 1_c_int)
      call check(error, rc == TRC_OK, "charge: set_charge(-2) status")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "charge: run_scf (-2) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_energy(h, e_minus2)
      rc = trc_destroy(h)

      ! Fresh handles, each built directly at its own charge.
      rc = trc_create(hf)
      rc = trc_set_molecule(hf, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 2_c_int, 1_c_int)
      rc = trc_set_basis_json(hf, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_jkfit_json(hf, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
      rc = trc_set_df(hf, 1_c_int)
      rc = trc_set_method(hf, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(hf, 1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(hf)
      call check(error, rc == TRC_OK, "charge: run_scf (fresh +2) must converge")
      if (allocated(error)) then
         rc = trc_destroy(hf)
         return
      end if
      rc = trc_energy(hf, e_plus2_fresh)
      rc = trc_destroy(hf)

      rc = trc_create(hf)
      rc = trc_set_molecule(hf, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), -2_c_int, 1_c_int)
      rc = trc_set_basis_json(hf, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_jkfit_json(hf, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
      rc = trc_set_df(hf, 1_c_int)
      rc = trc_set_method(hf, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(hf, 1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(hf)
      call check(error, rc == TRC_OK, "charge: run_scf (fresh -2) must converge")
      if (allocated(error)) then
         rc = trc_destroy(hf)
         return
      end if
      rc = trc_energy(hf, e_minus2_fresh)
      rc = trc_destroy(hf)

      call check(error, abs(real(e_plus2, dp) - real(e_plus2_fresh, dp)) < 1.0e-8_dp, &
                 "charge: +2 via trc_set_charge matches a fresh handle, reused "// &
                 to_string(real(e_plus2, dp), "f20.12")//"  fresh "// &
                 to_string(real(e_plus2_fresh, dp), "f20.12"))
      if (.not. allocated(error)) then
         call check(error, abs(real(e_minus2, dp) - real(e_minus2_fresh, dp)) < 1.0e-8_dp, &
                    "charge: -2 via trc_set_charge matches a fresh handle, reused "// &
                    to_string(real(e_minus2, dp), "f20.12")//"  fresh "// &
                    to_string(real(e_minus2_fresh, dp), "f20.12"))
      end if
   end subroutine test_charge_reuses_df_fit

end module check_capi_charge_unrestricted_suite

program check_capi_charge_unrestricted
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_bind_device
   use check_capi_charge_unrestricted_suite, only: collect_capi_charge_unrestricted_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, dev

   dev = trc_bind_device(0)
   stat = 0
   ! Serial: each case builds its own context through the C ABI, which is
   ! not thread-safe, and an -fopenmp build runs test-drive cases at once.
   testsuites = [new_testsuite("check_capi_charge_unrestricted", collect_capi_charge_unrestricted_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_charge_unrestricted
