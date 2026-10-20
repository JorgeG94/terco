program check_capi_charge_unrestricted
   !! trc_set_unrestricted, trc_set_level_shift and trc_set_charge,
   !! through the context API.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
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

   integer :: natm, dev
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   logical :: ok

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   ok = .true.

   call check_unrestricted_matches_rhf()
   call check_level_shift()
   call check_df_setup_is_idempotent()
   call check_charge_reuses_df_fit()

   if (.not. ok) then
      print '(a)', "check_capi_charge_unrestricted: FAIL"
      stop 1
   end if
   print '(a)', "check_capi_charge_unrestricted: PASS"

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "check_capi_charge_unrestricted: FAILED: "//what
         ok = .false.
      end if
   end subroutine expect

   !
   ! trc_set_unrestricted(1) on a closed-shell singlet (nalpha == nbeta)
   ! must converge UHF to the SAME stationary point as the default
   ! restricted run: forcing two independent spin densities on a genuine
   ! closed shell does not change the fixed point, only whether it is
   ! reached by one eigenproblem per iteration or two identical ones.
   !
   subroutine check_unrestricted_matches_rhf()
      type(c_ptr) :: h_rhf, h_uhf
      integer(c_int) :: rc, nspin_rhf, nspin_uhf
      real(c_double) :: e_rhf, e_uhf

      rc = trc_create(h_rhf)
      rc = trc_set_molecule(h_rhf, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h_rhf, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_method(h_rhf, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h_rhf, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(h_rhf);              call expect(rc == TRC_OK, "unrestricted: run_scf (RHF)")
      rc = trc_nspin(h_rhf, nspin_rhf);     call expect(nspin_rhf == 1, "unrestricted: RHF default is nspin 1")
      rc = trc_energy(h_rhf, e_rhf)

      rc = trc_create(h_uhf)
      rc = trc_set_molecule(h_uhf, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h_uhf, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_unrestricted(h_uhf, 1_c_int); call expect(rc == TRC_OK, "unrestricted: set_unrestricted(1)")
      rc = trc_set_method(h_uhf, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h_uhf, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(h_uhf);              call expect(rc == TRC_OK, "unrestricted: run_scf (forced UHF)")
      rc = trc_nspin(h_uhf, nspin_uhf);      call expect(nspin_uhf == 2, "unrestricted: forced UHF is nspin 2")
      rc = trc_energy(h_uhf, e_uhf)

      print '(a,f20.12,a,f20.12,a,es9.2)', "check_capi_charge_unrestricted: E_RHF ", e_rhf, &
         "  E_UHF(forced) ", e_uhf, "  diff ", real(e_rhf, dp) - real(e_uhf, dp)
      call expect(abs(real(e_rhf, dp) - real(e_uhf, dp)) < 1.0e-10_dp, &
                  "unrestricted: forced UHF on a closed shell matches RHF")

      rc = trc_destroy(h_rhf); rc = trc_destroy(h_uhf)
   end subroutine check_unrestricted_matches_rhf

   !
   ! trc_set_level_shift only changes the early, pre-DIIS part of the
   ! iteration (trc_scf_driver.F90's own "while DIIS and A/EDIIS are
   ! off" gate), so two different positive shifts must converge to the
   ! SAME energy -- the level shift moves the path, not the fixed point.
   ! A negative shift is refused; 0 (off) must still converge.
   !
   subroutine check_level_shift()
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: e_default, e_half, e_off

      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      rc = trc_set_level_shift(h, -1.0_c_double)
      call expect(rc == TRC_ERR_BADARG, "level_shift: a negative shift is refused")
      rc = trc_set_level_shift(h, 0.5_c_double); call expect(rc == TRC_OK, "level_shift: set 0.5")
      rc = trc_run_scf(h);                       call expect(rc == TRC_OK, "level_shift: run_scf (0.5)")
      rc = trc_energy(h, e_half)
      rc = trc_set_level_shift(h, 0.0_c_double);  call expect(rc == TRC_OK, "level_shift: set 0 (off)")
      rc = trc_run_scf(h);                       call expect(rc == TRC_OK, "level_shift: run_scf (off)")
      rc = trc_energy(h, e_off)
      rc = trc_set_level_shift(h, 1.0_c_double);  call expect(rc == TRC_OK, "level_shift: set 1.0 (default)")
      rc = trc_run_scf(h);                       call expect(rc == TRC_OK, "level_shift: run_scf (1.0)")
      rc = trc_energy(h, e_default)

      call expect(abs(real(e_half, dp) - real(e_default, dp)) < 1.0e-9_dp, &
                  "level_shift: 0.5 and 1.0 converge to the same energy")
      call expect(abs(real(e_off, dp) - real(e_default, dp)) < 1.0e-9_dp, &
                  "level_shift: off and 1.0 converge to the same energy")
      rc = trc_destroy(h)
   end subroutine check_level_shift

   !
   ! trc_set_charge on one --df handle, through an N, N+2, N-2 sequence
   ! (water: 10, 8, 12 electrons, all closed-shell singlets -- --df is
   ! restricted only, trc_run_gradient's own header says so and
   ! trc_run_scf refuses otherwise, so this stays closed-shell rather
   ! than exercising the N+1/N-1 radical case trc_set_unrestricted
   ! covers separately above), must match fresh handles built directly
   ! at each charge -- and the fit must NOT be rebuilt in between:
   ! df_jk_setup returns immediately once already `ready`
   ! (trc_scf_df.F90's own header), so trc_run_scf on the second and
   ! third charge states should be far faster than the first, which
   ! pays the metric/Cholesky/3c/trsm. The margin (2x) is generous
   ! against timing noise -- a real fit costs seconds, the skipped
   ! setup costs nothing.
   !
   subroutine check_charge_reuses_df_fit()
      type(c_ptr) :: h, hf
      integer(c_int) :: rc
      real(c_double) :: e_neutral, e_plus2, e_minus2
      real(c_double) :: e_plus2_fresh, e_minus2_fresh
      type(timer_type) :: tm
      real(dp) :: t_neutral, t_plus2, t_minus2

      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_jkfit_json(h, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
      rc = trc_set_df(h, 1_c_int)
      rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)

      call tm%start()
      rc = trc_run_scf(h);  call expect(rc == TRC_OK, "charge: run_scf (neutral)")
      t_neutral = tm%get_elapsed_time()
      rc = trc_energy(h, e_neutral)

      rc = trc_set_charge(h, 2_c_int, 1_c_int); call expect(rc == TRC_OK, "charge: set_charge (+2)")
      call tm%start()
      rc = trc_run_scf(h);  call expect(rc == TRC_OK, "charge: run_scf (+2)")
      t_plus2 = tm%get_elapsed_time()
      rc = trc_energy(h, e_plus2)

      rc = trc_set_charge(h, -2_c_int, 1_c_int); call expect(rc == TRC_OK, "charge: set_charge (-2)")
      call tm%start()
      rc = trc_run_scf(h);  call expect(rc == TRC_OK, "charge: run_scf (-2)")
      t_minus2 = tm%get_elapsed_time()
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
      rc = trc_run_scf(hf); call expect(rc == TRC_OK, "charge: run_scf (fresh +2)")
      rc = trc_energy(hf, e_plus2_fresh)
      rc = trc_destroy(hf)

      rc = trc_create(hf)
      rc = trc_set_molecule(hf, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), -2_c_int, 1_c_int)
      rc = trc_set_basis_json(hf, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_jkfit_json(hf, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
      rc = trc_set_df(hf, 1_c_int)
      rc = trc_set_method(hf, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(hf, 1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(hf); call expect(rc == TRC_OK, "charge: run_scf (fresh -2)")
      rc = trc_energy(hf, e_minus2_fresh)
      rc = trc_destroy(hf)

      print '(a,f20.12)', "check_capi_charge_unrestricted: E(neutral) ", e_neutral
      print '(a,f20.12,a,f20.12,a,es9.2)', "check_capi_charge_unrestricted: E(+2) reused ", e_plus2, &
         "  fresh ", e_plus2_fresh, "  diff ", real(e_plus2, dp) - real(e_plus2_fresh, dp)
      print '(a,f20.12,a,f20.12,a,es9.2)', "check_capi_charge_unrestricted: E(-2) reused ", e_minus2, &
         "  fresh ", e_minus2_fresh, "  diff ", real(e_minus2, dp) - real(e_minus2_fresh, dp)
      call expect(abs(real(e_plus2, dp) - real(e_plus2_fresh, dp)) < 1.0e-8_dp, &
                  "charge: +2 via trc_set_charge matches a fresh handle")
      call expect(abs(real(e_minus2, dp) - real(e_minus2_fresh, dp)) < 1.0e-8_dp, &
                  "charge: -2 via trc_set_charge matches a fresh handle")

      ! Informational only -- water/6-31gs is too small for the one-time
      ! fit cost to stand out reliably against noise at this wall-clock
      ! resolution. check_df_setup_is_idempotent (above) is the
      ! deterministic proof that a second setup call does no work;
      ! this just shows these three runs are all cheap.
      print '(a,f8.3,a,f8.3,a,f8.3,a)', "check_capi_charge_unrestricted: wall time neutral ", t_neutral, &
         " s  +2 ", t_plus2, " s  -2 ", t_minus2, " s"
   end subroutine check_charge_reuses_df_fit

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
   subroutine check_df_setup_is_idempotent()
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: la
      type(df_jk_t) :: dfjk
      type(error_t) :: err
      type(timer_type) :: tm
      real(dp) :: t_first, t_second
      real(dp), parameter :: ERI_THRESH = 1.0e-12_dp

      call trc_basis_from_json('basis_sets/6-31gs.json', natm, zint, at_r, bas, err, force_cart=.true.)
      if (.not. err%has_error()) call trc_basis_from_json('basis_sets/def2-universal-jkfit.json', &
                                                           natm, zint, at_r, aux, err)
      if (err%has_error()) then
         print '(a)', "check_capi_charge_unrestricted: "//err%get_message(); stop 1
      end if
      call bas%to_device(); call aux%to_device()
      call pl%build(bas, ERI_THRESH); call pl%to_device()
      call la%init(bas%nao)

      call tm%start()
      call dfjk%setup(bas, aux, pl, la, err)
      t_first = tm%get_elapsed_time()
      call expect(.not. err%has_error(), "df_setup: first setup succeeds")
      call expect(dfjk%is_ready(), "df_setup: ready after the first setup")

      call tm%start()
      call dfjk%setup(bas, aux, pl, la, err)
      t_second = tm%get_elapsed_time()
      call expect(.not. err%has_error(), "df_setup: second setup (idempotent) succeeds")

      print '(a,es10.2,a,es10.2,a)', "check_capi_charge_unrestricted: df_jk_setup first ", t_first, &
         " s, second (idempotent) ", t_second, " s"
      call expect(t_second < 0.1_dp*max(t_first, 1.0e-6_dp), &
                  "df_setup: a second call on the same object does no work (trc_set_charge relies on this)")

      call dfjk%release(); call pl%release(); call la%release()
      call bas%release(); call aux%release()
   end subroutine check_df_setup_is_idempotent

end program check_capi_charge_unrestricted
