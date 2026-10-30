!
! The context surface against the Fortran drivers: water/cc-pVDZ RHF and
! RI-MP2 through one handle, set up in steps, must reproduce trc_scf_run
! and trc_rimp2_run to the last digit, and every misuse must be a status.
!
! FORCED CARTESIAN. This test is about the C API surface, not about
! cc-pVDZ's physics, and RI-MP2 (trc_run_rimp2/trc_rimp2_run) refuses a
! spherical orbital basis outright (trc_capi.F90's own guard) -- so every
! basis load here, Fortran-driver and handle alike, forces Cartesian,
! which is what keeps this test exercising the API it always has rather
! than failing on that refusal.
!
! Six independent test-drive cases below. Each accumulates its own
! `expect()` failures (continuing past the first, like the original
! program's `ok` flag) into one message, reported through a single
! `check` at the end -- test-drive's `check(error, ...)` has
! `intent(out)` on `error`, so calling it more than once per test would
! silently wipe an earlier failure rather than add to it.
!
module context_scf_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char, c_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t, trc_rimp2_options_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_force_cartesian, &
                               trc_set_method, trc_set_convergence, trc_set_rimp2, trc_run_scf, trc_run_rimp2, &
                               trc_nao, trc_nspin, trc_energy, trc_converged, trc_iterations, trc_density, &
                               trc_mo_energies, trc_rimp2_energy, trc_message, trc_set_guess, &
                               trc_run_gradient, trc_gradient, trc_set_df, trc_set_jkfit_json, trc_cart_order, &
                               trc_abi_version, trc_functional_info, trc_functional_list, &
                               TRC_OK, TRC_ERR_STATE, TRC_ERR_NULL, TRC_ERR_BADARG, TRC_GUESS_GWH, &
                               TRC_ABI_MAJOR, TRC_ABI_MINOR
   implicit none
   private
   public :: collect_context_scf

contains

   subroutine collect_context_scf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("misuse before setup is always a status, never a crash", test_misuse), &
                  new_unittest("the handle reproduces the Fortran drivers", test_handle_reproduces_driver), &
                  new_unittest("b2plyp double hybrid total", test_dh_b2plyp), &
                  new_unittest("b2gpplyp double hybrid total", test_dh_b2gpplyp), &
                  new_unittest("aux/orbital basis convention order", test_aux_order), &
                  new_unittest("ABI version and functional_info/functional_list", test_abi_and_functionals) &
                  ]
   end subroutine collect_context_scf

   subroutine expect(ok, msgs, cond, what)
      logical, intent(inout) :: ok
      character(len=:), allocatable, intent(inout) :: msgs
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "context_scf: FAILED: "//what
         ok = .false.
         if (.not. allocated(msgs)) msgs = ""
         msgs = msgs//new_line('a')//"  FAILED: "//what
      end if
   end subroutine expect

   subroutine test_misuse(error)
      !! Every misuse before (or instead of) correct setup is a status,
      !! never a crash.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      character(kind=c_char) :: msg(256)
      logical :: ok
      character(len=:), allocatable :: msgs

      ok = .true.
      rc = trc_destroy(c_null_ptr);            call expect(ok, msgs, rc == TRC_ERR_NULL, "destroy(null) is NULL")
      rc = trc_create(h);                      call expect(ok, msgs, rc == TRC_OK, "create")
      rc = trc_run_scf(h);                     call expect(ok, msgs, rc == TRC_ERR_STATE, "run before setup is STATE")
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      call expect(ok, msgs, rc == TRC_ERR_STATE, "basis before molecule is STATE")
      rc = trc_run_rimp2(h);                   call expect(ok, msgs, rc == TRC_ERR_STATE, "rimp2 before scf is STATE")
      rc = trc_message(h, msg, 256_c_int)
      call expect(ok, msgs, rc == TRC_OK .and. msg(1) /= c_null_char, "a refusal leaves a message")
      rc = trc_set_method(h, 'pbe99'//c_null_char, 3_c_int)
      call expect(ok, msgs, rc /= TRC_OK, "an unknown functional is refused at the setter")
      rc = trc_message(h, msg, 256_c_int)
      call expect(ok, msgs, index(transfer(msg, repeat(' ', 256)), 'pbe99') > 0, "and named")
      rc = trc_set_method(h, ''//c_null_char, -1_c_int)
      call expect(ok, msgs, rc /= TRC_OK, "a negative grid level is refused")
      rc = trc_destroy(h)

      call check(error, ok, msgs)
   end subroutine test_misuse

   subroutine test_handle_reproduces_driver(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm, nocc
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), dsad(:, :), dg3(:, :, :), dlib(:, :, :), elib(:, :)
      type(trc_basis_t) :: b, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: ref
      type(trc_rimp2_result_t) :: mp2
      type(error_t) :: err
      type(c_ptr) :: h
      integer(c_int) :: rc, n, ns, it, flag
      real(c_double) :: e, e_os, e_ss, e_df
      real(c_double), allocatable :: grad(:)
      character(kind=c_char) :: msg(256)
      logical :: ok
      character(len=:), allocatable :: msgs
      integer :: dev

      ok = .true.
      dev = trc_bind_device(0)
      call read_xyz('water.xyz', natm, zint, at_r)
      nocc = sum(zint)/2

      ! --- the reference, through the Fortran drivers, from terco's own SAD ----
      call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, b, err, force_cart=.true.)
      if (.not. err%has_error()) call trc_basis_from_json('basis_sets/cc-pvdz-rifit.json', natm, zint, at_r, aux, err, &
                                                           force_cart=.true.)
      if (.not. err%has_error()) call trc_sad_build(b, dsad, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call b%to_device(); call aux%to_device()
      opts%conv_energy = 1.0e-10_dp; opts%conv_diis = 1.0e-6_dp
      allocate (dg3(b%nao, b%nao, 1)); dg3(:, :, 1) = dsad
      call trc_scf_run(b, nocc, nocc, opts, ref, dguess=dg3)
      call pl%build(b, 1.0e-12_dp); call pl%to_device()
      block
         type(trc_rimp2_options_t) :: rimp2_opts
         rimp2_opts%nfrozen = 1
         rimp2_opts%aux_block = 20
         call trc_rimp2_run(b, aux, pl, nocc, ref%cmo(:, :, 1), ref%eps(:, 1), mp2, opts=rimp2_opts)
      end block

      ! --- the same calculation through the handle ----------------------------
      rc = trc_create(h);                      call expect(ok, msgs, rc == TRC_OK, "create")
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(ok, msgs, rc == TRC_OK, "set_molecule")
      rc = trc_set_force_cartesian(h, 1_c_int)
      call expect(ok, msgs, rc == TRC_OK, "set_force_cartesian")
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      call expect(ok, msgs, rc == TRC_OK, "set_basis_json")
      rc = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      call expect(ok, msgs, rc == TRC_OK, "set_aux_json")
      rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      call expect(ok, msgs, rc == TRC_OK, "set_method")
      rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
      call expect(ok, msgs, rc == TRC_OK, "set_convergence")
      rc = trc_set_rimp2(h, 1_c_int, 20_c_int)
      call expect(ok, msgs, rc == TRC_OK, "set_rimp2")
      rc = trc_nao(h, n);                      call expect(ok, msgs, rc == TRC_OK .and. n == b%nao, "nao")
      rc = trc_run_scf(h);                     call expect(ok, msgs, rc == TRC_OK, "run_scf")
      rc = trc_nspin(h, ns);                   call expect(ok, msgs, rc == TRC_OK .and. ns == 1, "nspin")
      rc = trc_converged(h, flag);             call expect(ok, msgs, flag == 1, "converged")
      rc = trc_iterations(h, it);              call expect(ok, msgs, it == ref%iterations, "same iteration count as the driver")
      rc = trc_energy(h, e)
      print '(a,f20.12,a,f20.12,a,es9.2)', "context_scf: E_scf handle ", e, "  driver ", ref%energy, "  diff ", e - ref%energy
      call expect(ok, msgs, abs(e - ref%energy) < 1.0e-10_dp, "SCF energy matches the driver")
      allocate (dlib(n, n, 1), elib(n, 1))
      rc = trc_density(h, dlib); rc = trc_mo_energies(h, elib)
      call expect(ok, msgs, maxval(abs(dlib(:, :, 1) - ref%dmat(:, :, 1))) < 1.0e-8_dp, "density matches")
      call expect(ok, msgs, maxval(abs(elib(:, 1) - ref%eps(:, 1))) < 1.0e-8_dp, "orbital energies match")
      ! the auxiliary basis set AFTER the SCF must leave the SCF standing: that
      ! is the order a driver naturally uses
      rc = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      call expect(ok, msgs, rc == TRC_OK, "aux after scf")
      rc = trc_converged(h, flag);             call expect(ok, msgs, rc == TRC_OK .and. flag == 1, "SCF survives a new aux basis")
      rc = trc_run_rimp2(h);                   call expect(ok, msgs, rc == TRC_OK, "run_rimp2")
      rc = trc_rimp2_energy(h, e_os, e_ss)
      print '(a,f20.12,a,f20.12,a,es9.2)', "context_scf: E_corr handle ", e_os + e_ss, "  driver ", mp2%e_os + mp2%e_ss, &
         "  diff ", e_os + e_ss - mp2%e_os - mp2%e_ss
      call expect(ok, msgs, abs(e_os - mp2%e_os) < 1.0e-9_dp .and. abs(e_ss - mp2%e_ss) < 1.0e-9_dp, &
                  "RI-MP2 matches the driver")
      rc = trc_energy(h, e)
      call expect(ok, msgs, abs(e - (ref%energy + mp2%e_os + mp2%e_ss)) < 1.0e-9_dp, "total includes the correlation")

      ! --- a setter after the run invalidates, and a second run from GWH agrees -
      rc = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
      call expect(ok, msgs, rc == TRC_OK, "set_guess gwh")
      rc = trc_rimp2_energy(h, e_os, e_ss);    call expect(ok, msgs, rc == TRC_ERR_STATE, "rimp2 result gone after a setter")
      rc = trc_run_scf(h);                     call expect(ok, msgs, rc == TRC_OK, "second run")
      rc = trc_energy(h, e);                   call expect(ok, msgs, abs(e - ref%energy) < 1.0e-9_dp, "GWH start, same energy")

      ! --- the entry points host codes reach only through trc_c_interfaces:
      ! declared there by hand, so a wrong declaration shows up here ---
      call expect(ok, msgs, trc_cart_order() == 0 .or. trc_cart_order() == 1, "cart_order is 0 or 1")
      allocate (grad(3*natm))
      rc = trc_run_gradient(h);                call expect(ok, msgs, rc == TRC_OK, "run_gradient")
      rc = trc_gradient(h, grad);              call expect(ok, msgs, rc == TRC_OK, "gradient")
      call expect(ok, msgs, maxval(abs(grad)) > 1.0e-4_dp, "gradient is not zero")
      call expect(ok, msgs, abs(sum(grad(1::3))) + abs(sum(grad(2::3))) + abs(sum(grad(3::3))) < 1.0e-8_dp, &
                  "gradient is translationally invariant")
      rc = trc_set_jkfit_json(h, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
      call expect(ok, msgs, rc == TRC_OK, "set_jkfit_json")
      rc = trc_set_df(h, 1_c_int);             call expect(ok, msgs, rc == TRC_OK, "set_df")
      rc = trc_run_scf(h);                     call expect(ok, msgs, rc == TRC_OK, "run_scf --df")
      rc = trc_energy(h, e_df)
      ! The fitting error itself: tens of microhartree, never zero.
      call expect(ok, msgs, abs(e_df - ref%energy) > 0.0_dp .and. abs(e_df - ref%energy) < 1.0e-3_dp, &
                  "--df energy near the 4c one")
      rc = trc_set_df(h, 0_c_int);             call expect(ok, msgs, rc == TRC_OK, "set_df off")
      rc = trc_message(h, msg, 256_c_int);     call expect(ok, msgs, rc == TRC_OK, "message")
      rc = trc_destroy(h);                     call expect(ok, msgs, rc == TRC_OK, "destroy")

      call pl%release(); call b%release(); call aux%release()
      call check(error, ok, msgs)
   end subroutine test_handle_reproduces_driver

   subroutine test_dh_b2plyp(error)
      type(error_type), allocatable, intent(out) :: error
      call run_double_hybrid(error, 'b2plyp', 0.27_c_double, 0.27_c_double)
   end subroutine test_dh_b2plyp

   subroutine test_dh_b2gpplyp(error)
      type(error_type), allocatable, intent(out) :: error
      call run_double_hybrid(error, 'b2gpplyp', 0.36_c_double, 0.36_c_double)
   end subroutine test_dh_b2gpplyp

   ! water/cc-pVDZ through a fresh handle, `name` a double hybrid: checks
   ! trc_energy's combined total against E_scf + pt2_cos_lit*E_os +
   ! pt2_css_lit*E_ss, the LITERAL coefficients the caller passes in
   ! (never read back from trc_xc_functional's own table), and that E_os
   ! and E_ss are not coincidentally equal.
   !
   ! Every double hybrid terco has today (b2plyp, b2gp-plyp, mpw2plyp) has
   ! pt2_cos == pt2_css, and nothing public (no functional, no C-API
   ! setter) can set them unequally -- trc_xc_functional.F90's own header
   ! comment on why the field is split in the first place is forward-
   ! looking, not something exercised yet. A cos/css coefficient swap is
   ! therefore NOT distinguishable by this test on EITHER case below (the
   ! weighted sum is symmetric in cos/css when they are equal); what IS
   ! checked is that E_os and E_ss enter as two genuinely different
   ! numbers (asserted separately, not just a coincidental E_os == E_ss),
   ! each multiplied by its own literal coefficient.
   subroutine run_double_hybrid(error, name, pt2_cos_lit, pt2_css_lit)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name
      real(c_double), intent(in) :: pt2_cos_lit, pt2_css_lit
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(c_ptr) :: h2
      real(c_double) :: e_scf2, e_dh, e_os2, e_ss2, e_dh_expected
      integer(c_int) :: rc2
      logical :: ok
      character(len=:), allocatable :: msgs

      ok = .true.
      call read_xyz('water.xyz', natm, zint, at_r)

      rc2 = trc_create(h2);                      call expect(ok, msgs, rc2 == TRC_OK, name//": create")
      rc2 = trc_set_molecule(h2, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(ok, msgs, rc2 == TRC_OK, name//": set_molecule")
      rc2 = trc_set_force_cartesian(h2, 1_c_int)
      call expect(ok, msgs, rc2 == TRC_OK, name//": set_force_cartesian")
      rc2 = trc_set_basis_json(h2, 'basis_sets/cc-pvdz.json'//c_null_char)
      call expect(ok, msgs, rc2 == TRC_OK, name//": set_basis_json")
      rc2 = trc_set_aux_json(h2, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      call expect(ok, msgs, rc2 == TRC_OK, name//": set_aux_json")
      rc2 = trc_set_method(h2, name//c_null_char, 3_c_int); call expect(ok, msgs, rc2 == TRC_OK, name//": set_method")
      rc2 = trc_set_convergence(h2, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
      call expect(ok, msgs, rc2 == TRC_OK, name//": set_convergence")
      rc2 = trc_set_rimp2(h2, 1_c_int, 20_c_int);    call expect(ok, msgs, rc2 == TRC_OK, name//": set_rimp2")
      rc2 = trc_run_scf(h2);                         call expect(ok, msgs, rc2 == TRC_OK, name//": run_scf")
      rc2 = trc_energy(h2, e_scf2);                  call expect(ok, msgs, rc2 == TRC_OK, name//": E_scf alone")
      rc2 = trc_run_rimp2(h2);                       call expect(ok, msgs, rc2 == TRC_OK, name//": run_rimp2")
      rc2 = trc_rimp2_energy(h2, e_os2, e_ss2);      call expect(ok, msgs, rc2 == TRC_OK, name//": raw E_os/E_ss")
      rc2 = trc_energy(h2, e_dh);                    call expect(ok, msgs, rc2 == TRC_OK, name//": weighted total")

      call expect(ok, msgs, abs(e_os2 - e_ss2) > 1.0e-3_dp, name//": E_os and E_ss are genuinely different, not E_os == E_ss")

      e_dh_expected = e_scf2 + pt2_cos_lit*e_os2 + pt2_css_lit*e_ss2
      print '(a,a,a,f20.12,a,f20.12,a,es9.2)', "context_scf: ", name, " E_dh handle ", e_dh, &
         "  E_scf+cos*os+css*ss (literals) ", e_dh_expected, "  diff ", e_dh - e_dh_expected
      call expect(ok, msgs, abs(e_dh - e_dh_expected) < 1.0e-12_dp, &
                  name//": double hybrid total = E_scf + pt2_cos*E_os + pt2_css*E_ss (literals) to round-off")
      rc2 = trc_destroy(h2);                         call expect(ok, msgs, rc2 == TRC_OK, name//": destroy")

      call check(error, ok, msgs)
   end subroutine run_double_hybrid

   subroutine test_aux_order(error)
      !! The auxiliary (or jkfit) basis must end up on the SAME Cartesian/
      !! spherical convention as the orbital basis regardless of which is
      !! set first (pyscf's own auxmol.cart = mol.cart rule) --
      !! trc_set_basis_json/trc_set_aux_json only inherit the convention
      !! FORWARD (aux from an already-loaded orbital basis), so an aux
      !! basis set FIRST had nothing to inherit from and could silently
      !! disagree once the orbital basis loaded after it. Two bases with
      !! genuinely different native tags exercise both orders: 6-31g.json
      !! has no spherical shell anywhere (always Cartesian, no
      !! force_cartesian needed), cc-pvdz-rifit.json is spherical
      !! throughout -- orbital-first inherits correctly; aux-first must
      !! now be refused rather than silently mismatched.
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(c_ptr) :: h3, h4
      integer(c_int) :: rc3
      logical :: ok
      character(len=:), allocatable :: msgs

      ok = .true.
      call read_xyz('water.xyz', natm, zint, at_r)

      ! Correct order: orbital first, aux second -- aux inherits the
      ! orbital's own Cartesian convention and succeeds.
      rc3 = trc_create(h3);                      call expect(ok, msgs, rc3 == TRC_OK, "create (orbital-first)")
      rc3 = trc_set_molecule(h3, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(ok, msgs, rc3 == TRC_OK, "set_molecule (orbital-first)")
      rc3 = trc_set_basis_json(h3, 'basis_sets/6-31g.json'//c_null_char)
      call expect(ok, msgs, rc3 == TRC_OK, "orbital first succeeds")
      rc3 = trc_set_aux_json(h3, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      call expect(ok, msgs, rc3 == TRC_OK, "aux after orbital inherits and succeeds")
      rc3 = trc_destroy(h3);                      call expect(ok, msgs, rc3 == TRC_OK, "destroy (orbital-first)")

      ! Wrong order: aux first -- nothing to inherit from yet, so it
      ! keeps its own native (spherical) tag; orbital second keeps ITS
      ! own native (Cartesian) tag -- now refused instead of silently
      ! running RI-MP2 with two AO spaces that disagree.
      rc3 = trc_create(h4);                      call expect(ok, msgs, rc3 == TRC_OK, "create (aux-first)")
      rc3 = trc_set_molecule(h4, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(ok, msgs, rc3 == TRC_OK, "set_molecule (aux-first)")
      rc3 = trc_set_aux_json(h4, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      call expect(ok, msgs, rc3 == TRC_OK, "aux loaded first, nothing to inherit from yet, succeeds")
      rc3 = trc_set_basis_json(h4, 'basis_sets/6-31g.json'//c_null_char)
      call expect(ok, msgs, rc3 /= TRC_OK, "orbital after a mismatched aux is refused, not silently wrong")
      rc3 = trc_destroy(h4);                      call expect(ok, msgs, rc3 == TRC_OK, "destroy (aux-first)")

      call check(error, ok, msgs)
   end subroutine test_aux_order

   subroutine test_abi_and_functionals(error)
      !! trc_abi_version and trc_functional_info/trc_functional_list --
      !! both handle-free, so nothing needs to run first.
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc4, major, minor, flags, reqlen
      real(c_double) :: exx, c_os2, c_ss2
      character(kind=c_char), allocatable :: buf(:)
      character(len=:), allocatable :: text
      logical :: ok
      character(len=:), allocatable :: msgs

      ok = .true.

      rc4 = trc_abi_version(major, minor)
      call expect(ok, msgs, rc4 == TRC_OK, "abi_version: status")
      call expect(ok, msgs, major == TRC_ABI_MAJOR .and. minor == TRC_ABI_MINOR, &
                  "abi_version: matches the header's own TRC_ABI_MAJOR/MINOR")

      ! a known, runnable global hybrid
      rc4 = trc_functional_info('b3lyp'//c_null_char, flags, exx, c_os2, c_ss2)
      call expect(ok, msgs, rc4 == TRC_OK, "functional_info(b3lyp): status")
      call expect(ok, msgs, iand(flags, 1_c_int) /= 0, "functional_info(b3lyp): known (bit 1)")
      call expect(ok, msgs, iand(flags, 2_c_int) /= 0, "functional_info(b3lyp): hybrid (bit 2)")
      call expect(ok, msgs, iand(flags, 32_c_int) == 0, "functional_info(b3lyp): not a double hybrid")
      call expect(ok, msgs, abs(exx - 0.2_c_double) < 1.0e-12_c_double, "functional_info(b3lyp): exx = 0.2")

      ! a known double hybrid: bit 32, and the SAME pt2_cos/pt2_css
      ! the double-hybrid tests already asserted through trc_energy
      rc4 = trc_functional_info('b2plyp'//c_null_char, flags, exx, c_os2, c_ss2)
      call expect(ok, msgs, rc4 == TRC_OK, "functional_info(b2plyp): status")
      call expect(ok, msgs, iand(flags, 32_c_int) /= 0, "functional_info(b2plyp): double hybrid (bit 32)")
      call expect(ok, msgs, abs(c_os2 - 0.27_c_double) < 1.0e-12_c_double .and. abs(c_ss2 - 0.27_c_double) < 1.0e-12_c_double, &
                  "functional_info(b2plyp): c_os = c_ss = 0.27")

      ! a known range-separated hybrid with VV10
      rc4 = trc_functional_info('wb97x-v'//c_null_char, flags, exx, c_os2, c_ss2)
      call expect(ok, msgs, rc4 == TRC_OK, "functional_info(wb97x-v): status")
      call expect(ok, msgs, iand(flags, 1_c_int) /= 0, "functional_info(wb97x-v): known (bit 1)")
      call expect(ok, msgs, iand(flags, 2_c_int) /= 0, "functional_info(wb97x-v): hybrid (bit 2)")
      call expect(ok, msgs, iand(flags, 4_c_int) /= 0, "functional_info(wb97x-v): range-separated (bit 4)")
      call expect(ok, msgs, iand(flags, 8_c_int) /= 0, "functional_info(wb97x-v): VV10 (bit 8)")

      ! a runnable meta-GGA: bit 16 with bit 1, not 64
      rc4 = trc_functional_info('tpss'//c_null_char, flags, exx, c_os2, c_ss2)
      call expect(ok, msgs, rc4 == TRC_OK, "functional_info(tpss): status")
      call expect(ok, msgs, iand(flags, 1_c_int) /= 0, "functional_info(tpss): known (bit 1)")
      call expect(ok, msgs, iand(flags, 16_c_int) /= 0, "functional_info(tpss): meta-GGA (bit 16)")
      call expect(ok, msgs, iand(flags, 64_c_int) == 0, "functional_info(tpss): runnable, so not bit 64")

      ! listed but not runnable: bit 64, bit 1 (known) must then be clear;
      ! r2SCAN-L is a meta-GGA that needs the density Laplacian
      rc4 = trc_functional_info('r2scanl'//c_null_char, flags, exx, c_os2, c_ss2)
      call expect(ok, msgs, rc4 == TRC_OK, "functional_info(r2scanl): status (listed)")
      call expect(ok, msgs, iand(flags, 64_c_int) /= 0, "functional_info(r2scanl): listed but not runnable (bit 64)")
      call expect(ok, msgs, iand(flags, 1_c_int) == 0, "functional_info(r2scanl): NOT known, since it cannot run")
      call expect(ok, msgs, iand(flags, 16_c_int) /= 0, "functional_info(r2scanl): meta-GGA (bit 16)")

      ! an unknown name: non-OK, flags = 0
      rc4 = trc_functional_info('not-a-functional'//c_null_char, flags, exx, c_os2, c_ss2)
      call expect(ok, msgs, rc4 /= TRC_OK, "functional_info(unknown): non-OK status")
      call expect(ok, msgs, flags == 0, "functional_info(unknown): flags = 0")

      ! trc_functional_list: query the length with buflen = 0 (buf itself
      ! still needs SOME storage behind the dummy's assumed-size `buf(*)`,
      ! even though the library never reads or writes it at buflen = 0).
      allocate (buf(1))
      rc4 = trc_functional_list(buf, 0_c_int, reqlen)
      call expect(ok, msgs, rc4 == TRC_OK .and. reqlen > 0, "functional_list: buflen=0 query returns reqlen > 0")
      deallocate (buf)
      allocate (buf(reqlen))
      rc4 = trc_functional_list(buf, reqlen, reqlen)
      call expect(ok, msgs, rc4 == TRC_OK, "functional_list: filling call with the queried reqlen")
      allocate (character(len=reqlen) :: text)
      text = transfer(buf, text)
      call expect(ok, msgs, index(text, 'b3lyp') > 0, "functional_list: names b3lyp")
      call expect(ok, msgs, index(text, char(0)) == reqlen, "functional_list: NUL-terminated at reqlen")
      deallocate (buf, text)

      ! a buffer that is one byte too small is refused, untouched
      allocate (buf(reqlen - 1))
      rc4 = trc_functional_list(buf, reqlen - 1, reqlen)
      call expect(ok, msgs, rc4 == TRC_ERR_BADARG, "functional_list: a too-small buflen is refused")
      deallocate (buf)

      call check(error, ok, msgs)
   end subroutine test_abi_and_functionals

end module context_scf_suite

program context_scf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use context_scf_suite, only: collect_context_scf
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_context_scf, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program context_scf
