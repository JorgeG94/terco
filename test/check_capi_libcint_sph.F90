program check_capi_libcint_sph
   !! trc_set_basis_libcint/trc_set_aux_libcint/trc_set_jkfit_libcint with
   !! cartesian = 0 (spherical), against the equivalent trc_set_basis_json/
   !! trc_set_aux_json/trc_set_jkfit_json path.
   !!
   !! The libcint path builds exactly the same Cartesian Gaussians either
   !! way (basis_from_libcint has no cartesian/spherical bit of its own to
   !! read -- see trc_capi.F90's own comment on `basis_libcint`); what
   !! `cartesian = 0` changes is only the `trc_basis_t%spherical` tag, the
   !! same one `trc_basis_from_json` sets from a JSON file's
   !! "gto_spherical" function_type. If that tag is all that matters, a
   !! molecule built spherical through libcint arrays and through JSON
   !! must give the SAME SCF energy, the SAME density, the SAME MO count
   !! and -- since terco's spherical AOs are built to land m = -l..l,
   !! p as x/y/z, matching pyscf/libcint position for position
   !! (scripts/gen_sph.py) -- the SAME MO coefficients, up to an orbital's
   !! overall sign. This program is that check.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char, c_char
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz, pack_libcint
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_basis_libcint, trc_set_aux_libcint, trc_set_aux_json, &
                               trc_set_jkfit_json, trc_set_jkfit_libcint, trc_set_df, &
                               trc_set_method, trc_set_convergence, trc_set_ao_norm, &
                               trc_run_scf, trc_energy, trc_nao, trc_density, trc_mo_coefficients, &
                               trc_message, TRC_OK, TRC_AO_NORM_UNIT
   implicit none


   integer :: natm, dev
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   logical :: ok

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   ok = .true.

   call check_orbital_basis()
   call check_jkfit()
   call check_mismatch_refusals()

   if (.not. ok) then
      print '(a)', "check_capi_libcint_sph: FAIL"
      stop 1
   end if
   print '(a)', "check_capi_libcint_sph: PASS"

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "check_capi_libcint_sph: FAILED: "//what
         ok = .false.
      end if
   end subroutine expect

   !
   ! water/cc-pVDZ, spherical, RHF: trc_set_basis_json vs
   ! trc_set_basis_libcint(cartesian=0). Energy, nao (the SCF's own
   ! nao_sph), density and MO coefficients must agree to round-off.
   !
   subroutine check_orbital_basis()
      type(trc_basis_t) :: bj
      type(error_t) :: err
      type(c_ptr) :: hj, hl
      integer(c_int) :: rc, n_json, n_lib
      real(c_double) :: e_json, e_lib
      real(c_double), allocatable :: d_json(:, :), d_lib(:, :), c_json(:, :), c_lib(:, :)
      integer :: i, j

      call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, bj, err)
      if (err%has_error()) then
         print '(a)', "check_capi_libcint_sph: "//err%get_message(); stop 1
      end if
      call expect(bj%spherical, "cc-pVDZ loads spherical from its own JSON tag")

      rc = trc_create(hj);                   call expect(rc == TRC_OK, "orbital: create (json)")
      rc = trc_set_molecule(hj, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(rc == TRC_OK, "orbital: set_molecule (json)")
      rc = trc_set_basis_json(hj, 'basis_sets/cc-pvdz.json'//c_null_char)
      call expect(rc == TRC_OK, "orbital: set_basis_json")
      rc = trc_set_method(hj, ''//c_null_char, 3_c_int); call expect(rc == TRC_OK, "orbital: set_method (json)")
      rc = trc_set_convergence(hj, 1.0e-11_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
      call expect(rc == TRC_OK, "orbital: set_convergence (json)")
      rc = trc_run_scf(hj);                   call expect(rc == TRC_OK, "orbital: run_scf (json)")
      rc = trc_energy(hj, e_json);            call expect(rc == TRC_OK, "orbital: energy (json)")
      rc = trc_nao(hj, n_json);               call expect(rc == TRC_OK, "orbital: nao (json)")

      rc = trc_create(hl);                   call expect(rc == TRC_OK, "orbital: create (libcint)")
      rc = trc_set_molecule(hl, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(rc == TRC_OK, "orbital: set_molecule (libcint)")
      call set_basis_libcint(hl, bj, aux=.false., cartesian=0_c_int)
      rc = trc_set_method(hl, ''//c_null_char, 3_c_int); call expect(rc == TRC_OK, "orbital: set_method (libcint)")
      rc = trc_set_convergence(hl, 1.0e-11_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
      call expect(rc == TRC_OK, "orbital: set_convergence (libcint)")
      rc = trc_run_scf(hl);                   call expect(rc == TRC_OK, "orbital: run_scf (libcint)")
      rc = trc_energy(hl, e_lib);             call expect(rc == TRC_OK, "orbital: energy (libcint)")
      rc = trc_nao(hl, n_lib);                call expect(rc == TRC_OK, "orbital: nao (libcint)")

      print '(a,f20.12,a,f20.12,a,es9.2)', "check_capi_libcint_sph: E_scf json ", e_json, &
         "  libcint ", e_lib, "  diff ", e_json - e_lib
      call expect(n_lib == n_json, "orbital: nao (MO count) matches between json and libcint")
      call expect(abs(e_json - e_lib) < 1.0e-10_dp, "orbital: SCF energy matches between json and libcint")

      allocate (d_json(n_json, n_json), d_lib(n_lib, n_lib))
      allocate (c_json(n_json, n_json), c_lib(n_lib, n_lib))
      rc = trc_density(hj, d_json);          call expect(rc == TRC_OK, "orbital: density (json)")
      rc = trc_density(hl, d_lib);           call expect(rc == TRC_OK, "orbital: density (libcint)")
      call expect(maxval(abs(d_json - d_lib)) < 1.0e-9_dp, "orbital: density matches between json and libcint")

      rc = trc_mo_coefficients(hj, c_json);  call expect(rc == TRC_OK, "orbital: mo_coefficients (json)")
      rc = trc_mo_coefficients(hl, c_lib);   call expect(rc == TRC_OK, "orbital: mo_coefficients (libcint)")
      ! Each orbital's overall sign is arbitrary (the eigensolver's own
      ! convention), so compare |C| column by column rather than C itself.
      block
         real(dp) :: worst
         worst = 0.0_dp
         do j = 1, n_json
            do i = 1, n_json
               worst = max(worst, abs(abs(c_json(i, j)) - abs(c_lib(i, j))))
            end do
         end do
         print '(a,es10.2)', "check_capi_libcint_sph: |C| agreement json vs libcint ", worst
         call expect(worst < 1.0e-7_dp, "orbital: |MO coefficients| match between json and libcint")
      end block

      rc = trc_destroy(hj); rc = trc_destroy(hl)
      call bj%release()
   end subroutine check_orbital_basis

   !
   ! trc_set_jkfit_libcint against trc_set_jkfit_json: the --df energy on
   ! a spherical orbital basis (cc-pVDZ) with def2-universal-jkfit (also
   ! spherical) must agree to round-off, whichever way the jkfit basis
   ! itself was built.
   !
   subroutine check_jkfit()
      type(trc_basis_t) :: bj, bk
      type(error_t) :: err
      type(c_ptr) :: hj, hl
      integer(c_int) :: rc
      real(c_double) :: e_json, e_lib

      call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, bj, err)
      if (.not. err%has_error()) call trc_basis_from_json('basis_sets/def2-universal-jkfit.json', &
                                                           natm, zint, at_r, bk, err)
      if (err%has_error()) then
         print '(a)', "check_capi_libcint_sph: "//err%get_message(); stop 1
      end if

      rc = trc_create(hj);                   call expect(rc == TRC_OK, "jkfit: create (json)")
      rc = trc_set_molecule(hj, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(hj, 'basis_sets/cc-pvdz.json'//c_null_char)
      rc = trc_set_jkfit_json(hj, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
      call expect(rc == TRC_OK, "jkfit: set_jkfit_json")
      rc = trc_set_df(hj, 1_c_int)
      rc = trc_set_method(hj, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(hj, 1.0e-11_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(hj);                   call expect(rc == TRC_OK, "jkfit: run_scf (json)")
      rc = trc_energy(hj, e_json);            call expect(rc == TRC_OK, "jkfit: energy (json)")

      rc = trc_create(hl);                   call expect(rc == TRC_OK, "jkfit: create (libcint)")
      rc = trc_set_molecule(hl, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call set_basis_libcint(hl, bj, aux=.false., cartesian=0_c_int)
      call set_jkfit_libcint(hl, bk, cartesian=0_c_int)
      rc = trc_set_df(hl, 1_c_int)
      rc = trc_set_method(hl, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(hl, 1.0e-11_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(hl);                   call expect(rc == TRC_OK, "jkfit: run_scf (libcint)")
      rc = trc_energy(hl, e_lib);             call expect(rc == TRC_OK, "jkfit: energy (libcint)")

      print '(a,f20.12,a,f20.12,a,es9.2)', "check_capi_libcint_sph: E_df json ", e_json, &
         "  libcint ", e_lib, "  diff ", e_json - e_lib
      call expect(abs(e_json - e_lib) < 1.0e-10_dp, "jkfit: --df energy matches between json and libcint")

      rc = trc_destroy(hj); rc = trc_destroy(hl)
      call bj%release(); call bk%release()
   end subroutine check_jkfit

   !
   ! Convention mismatches must still be refused over the libcint path,
   ! exactly as they are over JSON (trc_capi.F90's own "CARTESIAN
   ! COMPONENT ORDER"/"AO NORMALISATION" refusals): an orbital basis
   ! Cartesian while its aux/jkfit basis is spherical (or the reverse),
   ! and TRC_AO_NORM_UNIT together with a spherical basis.
   !
   subroutine check_mismatch_refusals()
      type(trc_basis_t) :: b_cart, b_sph
      type(error_t) :: err
      type(c_ptr) :: h
      integer(c_int) :: rc

      call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, b_cart, err, force_cart=.true.)
      if (.not. err%has_error()) call trc_basis_from_json('basis_sets/cc-pvdz-rifit.json', natm, zint, at_r, &
                                                           b_sph, err)
      if (err%has_error()) then
         print '(a)', "check_capi_libcint_sph: "//err%get_message(); stop 1
      end if

      ! orbital Cartesian, aux spherical -> refused
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call set_basis_libcint(h, b_cart, aux=.false., cartesian=1_c_int)
      call set_basis_libcint(h, b_sph, aux=.true., cartesian=0_c_int, expect_ok=.false.)
      rc = trc_destroy(h)

      ! orbital spherical, aux Cartesian -> refused
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call set_basis_libcint(h, b_sph, aux=.false., cartesian=0_c_int)
      call set_basis_libcint(h, b_cart, aux=.true., cartesian=1_c_int, expect_ok=.false.)
      rc = trc_destroy(h)

      ! TRC_AO_NORM_UNIT then a spherical libcint basis -> refused
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_ao_norm(h, TRC_AO_NORM_UNIT); call expect(rc == TRC_OK, "mismatch: set_ao_norm(UNIT)")
      call set_basis_libcint(h, b_sph, aux=.false., cartesian=0_c_int, expect_ok=.false.)
      rc = trc_destroy(h)

      call b_cart%release(); call b_sph%release()
   end subroutine check_mismatch_refusals

   !
   ! Pack `b`'s shells into libcint's atm/bas/env and call
   ! trc_set_basis_libcint/trc_set_aux_libcint on `h`. `b` must come
   ! straight from trc_basis_from_json (common_fac_sp folded in once by
   ! `build`); `unpack_basis` takes that fold back out before packing, so
   ! `from_libcint`'s own fold (inside `build`, again) lands on the right
   ! number -- the same round-trip link_gfortran.F90's own
   ! build_libcint_arrays uses.
   !
   subroutine set_basis_libcint(h, b, aux, cartesian, expect_ok)
      type(c_ptr), intent(in) :: h
      type(trc_basis_t), intent(in) :: b
      logical, intent(in) :: aux
      integer(c_int), intent(in) :: cartesian
      logical, intent(in), optional :: expect_ok
      integer(c_int), allocatable :: atm(:, :), bas(:, :)
      real(c_double), allocatable :: env(:)
      integer(c_int) :: rc
      logical :: want_ok

      want_ok = .true.
      if (present(expect_ok)) want_ok = expect_ok
      call pack_libcint(b, atm, bas, env)
      if (.not. aux) then
         rc = trc_set_basis_libcint(h, atm, int(b%natm, c_int), bas, int(b%nshell, c_int), env, &
                                    int(size(env), c_int), cartesian)
         call expect((rc == TRC_OK) .eqv. want_ok, "set_basis_libcint status as expected")
      else
         rc = trc_set_aux_libcint(h, atm, int(b%natm, c_int), bas, int(b%nshell, c_int), env, &
                                  int(size(env), c_int), cartesian)
         call expect((rc == TRC_OK) .eqv. want_ok, "set_aux_libcint status as expected")
      end if
   end subroutine set_basis_libcint

   subroutine set_jkfit_libcint(h, b, cartesian, expect_ok)
      type(c_ptr), intent(in) :: h
      type(trc_basis_t), intent(in) :: b
      integer(c_int), intent(in) :: cartesian
      logical, intent(in), optional :: expect_ok
      integer(c_int), allocatable :: atm(:, :), bas(:, :)
      real(c_double), allocatable :: env(:)
      integer(c_int) :: rc
      logical :: want_ok

      want_ok = .true.
      if (present(expect_ok)) want_ok = expect_ok
      call pack_libcint(b, atm, bas, env)
      rc = trc_set_jkfit_libcint(h, atm, int(b%natm, c_int), bas, int(b%nshell, c_int), env, &
                                 int(size(env), c_int), cartesian)
      call expect((rc == TRC_OK) .eqv. want_ok, "set_jkfit_libcint status as expected")
   end subroutine set_jkfit_libcint


end program check_capi_libcint_sph
