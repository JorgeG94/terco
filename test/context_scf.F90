program context_scf
   !! The context surface against the Fortran drivers: water/cc-pVDZ RHF and
   !! RI-MP2 through one handle, set up in steps, must reproduce trc_scf_run
   !! and trc_rimp2_run to the last digit, and every misuse must be a status.
   !!
   !! FORCED CARTESIAN. This test is about the C API surface, not about
   !! cc-pVDZ's physics, and RI-MP2 (trc_run_rimp2/trc_rimp2_run) refuses a
   !! spherical orbital basis outright (trc_capi.F90's own guard) -- so every
   !! basis load here, Fortran-driver and handle alike, forces Cartesian,
   !! which is what keeps this test exercising the API it always has rather
   !! than failing on that refusal.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char, c_char
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t
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
                               TRC_OK, TRC_ERR_STATE, TRC_ERR_NULL, TRC_GUESS_GWH
   implicit none
   integer :: natm, dev, nocc
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

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   nocc = sum(zint)/2
   ok = .true.

   ! --- the reference, through the Fortran drivers, from terco's own SAD ----
   call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, b, err, force_cart=.true.)
   if (.not. err%has_error()) call trc_basis_from_json('basis_sets/cc-pvdz-rifit.json', natm, zint, at_r, aux, err, &
                                                        force_cart=.true.)
   if (.not. err%has_error()) call trc_sad_build(b, dsad, err)
   if (err%has_error()) then
      print '(a)', "context_scf: "//err%get_message()
      stop 1
   end if
   call b%to_device(); call aux%to_device()
   opts%conv_energy = 1.0e-10_dp; opts%conv_diis = 1.0e-6_dp
   allocate (dg3(b%nao, b%nao, 1)); dg3(:, :, 1) = dsad
   call trc_scf_run(b, nocc, nocc, opts, ref, dguess=dg3)
   call pl%build(b, 1.0e-12_dp); call pl%to_device()
   call trc_rimp2_run(b, aux, pl, nocc, ref%cmo(:, :, 1), ref%eps(:, 1), mp2, nfrozen=1, aux_block=20)

   ! --- misuse first: every one a status, none a crash ---------------------
   rc = trc_destroy(c_null_ptr);            call expect(rc == TRC_ERR_NULL, "destroy(null) is NULL")
   rc = trc_create(h);                      call expect(rc == TRC_OK, "create")
   rc = trc_run_scf(h);                     call expect(rc == TRC_ERR_STATE, "run before setup is STATE")
   rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
   call expect(rc == TRC_ERR_STATE, "basis before molecule is STATE")
   rc = trc_run_rimp2(h);                   call expect(rc == TRC_ERR_STATE, "rimp2 before scf is STATE")
   rc = trc_message(h, msg, 256_c_int);     call expect(rc == TRC_OK .and. msg(1) /= c_null_char, "a refusal leaves a message")
   rc = trc_set_method(h, 'pbe99'//c_null_char, 3_c_int)
   call expect(rc /= TRC_OK, "an unknown functional is refused at the setter")
   rc = trc_message(h, msg, 256_c_int);     call expect(index(transfer(msg, repeat(' ', 256)), 'pbe99') > 0, "and named")
   rc = trc_set_method(h, ''//c_null_char, -1_c_int)
   call expect(rc /= TRC_OK, "a negative grid level is refused")

   ! --- the same calculation through the handle ----------------------------
   rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
   call expect(rc == TRC_OK, "set_molecule")
   rc = trc_set_force_cartesian(h, 1_c_int);                                  call expect(rc == TRC_OK, "set_force_cartesian")
   rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char);       call expect(rc == TRC_OK, "set_basis_json")
   rc = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char);   call expect(rc == TRC_OK, "set_aux_json")
   rc = trc_set_method(h, ''//c_null_char, 3_c_int);                          call expect(rc == TRC_OK, "set_method")
   rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
   call expect(rc == TRC_OK, "set_convergence")
   rc = trc_set_rimp2(h, 1_c_int, 20_c_int);                                  call expect(rc == TRC_OK, "set_rimp2")
   rc = trc_nao(h, n);                      call expect(rc == TRC_OK .and. n == b%nao, "nao")
   rc = trc_run_scf(h);                     call expect(rc == TRC_OK, "run_scf")
   rc = trc_nspin(h, ns);                   call expect(rc == TRC_OK .and. ns == 1, "nspin")
   rc = trc_converged(h, flag);             call expect(flag == 1, "converged")
   rc = trc_iterations(h, it);              call expect(it == ref%iterations, "same iteration count as the driver")
   rc = trc_energy(h, e)
   print '(a,f20.12,a,f20.12,a,es9.2)', "context_scf: E_scf handle ", e, "  driver ", ref%energy, "  diff ", e - ref%energy
   call expect(abs(e - ref%energy) < 1.0e-10_dp, "SCF energy matches the driver")
   allocate (dlib(n, n, 1), elib(n, 1))
   rc = trc_density(h, dlib); rc = trc_mo_energies(h, elib)
   call expect(maxval(abs(dlib(:, :, 1) - ref%dmat(:, :, 1))) < 1.0e-8_dp, "density matches")
   call expect(maxval(abs(elib(:, 1) - ref%eps(:, 1))) < 1.0e-8_dp, "orbital energies match")
   ! the auxiliary basis set AFTER the SCF must leave the SCF standing: that
   ! is the order a driver naturally uses
   rc = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char);   call expect(rc == TRC_OK, "aux after scf")
   rc = trc_converged(h, flag);             call expect(rc == TRC_OK .and. flag == 1, "SCF survives a new aux basis")
   rc = trc_run_rimp2(h);                   call expect(rc == TRC_OK, "run_rimp2")
   rc = trc_rimp2_energy(h, e_os, e_ss)
   print '(a,f20.12,a,f20.12,a,es9.2)', "context_scf: E_corr handle ", e_os + e_ss, "  driver ", mp2%e_os + mp2%e_ss, &
      "  diff ", e_os + e_ss - mp2%e_os - mp2%e_ss
   call expect(abs(e_os - mp2%e_os) < 1.0e-9_dp .and. abs(e_ss - mp2%e_ss) < 1.0e-9_dp, "RI-MP2 matches the driver")
   rc = trc_energy(h, e)
   call expect(abs(e - (ref%energy + mp2%e_os + mp2%e_ss)) < 1.0e-9_dp, "total includes the correlation")

   ! --- a setter after the run invalidates, and a second run from GWH agrees -
   rc = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int); call expect(rc == TRC_OK, "set_guess gwh")
   rc = trc_rimp2_energy(h, e_os, e_ss);    call expect(rc == TRC_ERR_STATE, "rimp2 result gone after a setter")
   rc = trc_run_scf(h);                     call expect(rc == TRC_OK, "second run")
   rc = trc_energy(h, e);                   call expect(abs(e - ref%energy) < 1.0e-9_dp, "GWH start, same energy")

   ! --- the entry points host codes reach only through trc_c_interfaces:
   ! declared there by hand, so a wrong declaration shows up here ---
   call expect(trc_cart_order() == 0 .or. trc_cart_order() == 1, "cart_order is 0 or 1")
   allocate (grad(3*natm))
   rc = trc_run_gradient(h);                call expect(rc == TRC_OK, "run_gradient")
   rc = trc_gradient(h, grad);              call expect(rc == TRC_OK, "gradient")
   call expect(maxval(abs(grad)) > 1.0e-4_dp, "gradient is not zero")
   call expect(abs(sum(grad(1::3))) + abs(sum(grad(2::3))) + abs(sum(grad(3::3))) < 1.0e-8_dp, &
               "gradient is translationally invariant")
   rc = trc_set_jkfit_json(h, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
   call expect(rc == TRC_OK, "set_jkfit_json")
   rc = trc_set_df(h, 1_c_int);             call expect(rc == TRC_OK, "set_df")
   rc = trc_run_scf(h);                     call expect(rc == TRC_OK, "run_scf --df")
   rc = trc_energy(h, e_df)
   ! The fitting error itself: tens of microhartree, never zero.
   call expect(abs(e_df - ref%energy) > 0.0_dp .and. abs(e_df - ref%energy) < 1.0e-3_dp, "--df energy near the 4c one")
   rc = trc_set_df(h, 0_c_int);             call expect(rc == TRC_OK, "set_df off")
   rc = trc_message(h, msg, 256_c_int);     call expect(rc == TRC_OK, "message")
   rc = trc_destroy(h);                     call expect(rc == TRC_OK, "destroy")

   ! --- a double hybrid's weighted total (trc_energy's pt2_cos/pt2_css
   ! branch): checked against LITERAL coefficients, not trc_xc_functional's
   ! own table, so a wrong table entry (or a coefficient accidentally
   ! copied from another functional) fails this test even though it would
   ! not fail a check that re-read the same table trc_energy itself reads.
   ! pyscf has no packaged functional that applies spin-component PT2
   ! weights itself -- it is the same manual E_scf + c_os*E_os + c_ss*E_ss
   ! composition terco uses (test/pyscf_refs/generate.py's gen_dh) -- so
   ! this checks self-consistency instead: the handle's combined
   ! trc_energy must equal E_scf (read back BEFORE trc_run_rimp2 sets
   ! cx%rimp2_ok) plus cos*E_os + css*E_ss, from the SAME raw E_os/E_ss
   ! trc_rimp2_energy reports, to round-off.
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
   call check_double_hybrid('b2plyp', 0.27_c_double, 0.27_c_double)
   call check_double_hybrid('b2gpplyp', 0.36_c_double, 0.36_c_double)
   call check_aux_order()

   call pl%release(); call b%release(); call aux%release()
   if (.not. ok) then
      print '(a)', "context_scf: FAIL"
      stop 1
   end if
   print '(a)', "context_scf: PASS"

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "context_scf: FAILED: "//what
         ok = .false.
      end if
   end subroutine expect

   ! water/cc-pVDZ through a fresh handle, `name` a double hybrid: checks
   ! trc_energy's combined total against E_scf + pt2_cos_lit*E_os +
   ! pt2_css_lit*E_ss, the LITERAL coefficients the caller passes in
   ! (never read back from trc_xc_functional's own table), and that E_os
   ! and E_ss are not coincidentally equal -- see this program's own
   ! comment above the two calls on why a cos/css swap is otherwise
   ! undetectable for every double hybrid terco has today.
   subroutine check_double_hybrid(name, pt2_cos_lit, pt2_css_lit)
      character(len=*), intent(in) :: name
      real(c_double), intent(in) :: pt2_cos_lit, pt2_css_lit
      type(c_ptr) :: h2
      real(c_double) :: e_scf2, e_dh, e_os2, e_ss2, e_dh_expected
      integer(c_int) :: rc2

      rc2 = trc_create(h2);                      call expect(rc2 == TRC_OK, name//": create")
      rc2 = trc_set_molecule(h2, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(rc2 == TRC_OK, name//": set_molecule")
      rc2 = trc_set_force_cartesian(h2, 1_c_int);    call expect(rc2 == TRC_OK, name//": set_force_cartesian")
      rc2 = trc_set_basis_json(h2, 'basis_sets/cc-pvdz.json'//c_null_char)
      call expect(rc2 == TRC_OK, name//": set_basis_json")
      rc2 = trc_set_aux_json(h2, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      call expect(rc2 == TRC_OK, name//": set_aux_json")
      rc2 = trc_set_method(h2, name//c_null_char, 3_c_int); call expect(rc2 == TRC_OK, name//": set_method")
      rc2 = trc_set_convergence(h2, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
      call expect(rc2 == TRC_OK, name//": set_convergence")
      rc2 = trc_set_rimp2(h2, 1_c_int, 20_c_int);    call expect(rc2 == TRC_OK, name//": set_rimp2")
      rc2 = trc_run_scf(h2);                         call expect(rc2 == TRC_OK, name//": run_scf")
      rc2 = trc_energy(h2, e_scf2);                  call expect(rc2 == TRC_OK, name//": E_scf alone")
      rc2 = trc_run_rimp2(h2);                       call expect(rc2 == TRC_OK, name//": run_rimp2")
      rc2 = trc_rimp2_energy(h2, e_os2, e_ss2);      call expect(rc2 == TRC_OK, name//": raw E_os/E_ss")
      rc2 = trc_energy(h2, e_dh);                    call expect(rc2 == TRC_OK, name//": weighted total")

      call expect(abs(e_os2 - e_ss2) > 1.0e-3_dp, name//": E_os and E_ss are genuinely different, not E_os == E_ss")

      e_dh_expected = e_scf2 + pt2_cos_lit*e_os2 + pt2_css_lit*e_ss2
      print '(a,a,a,f20.12,a,f20.12,a,es9.2)', "context_scf: ", name, " E_dh handle ", e_dh, &
         "  E_scf+cos*os+css*ss (literals) ", e_dh_expected, "  diff ", e_dh - e_dh_expected
      call expect(abs(e_dh - e_dh_expected) < 1.0e-12_dp, &
                  name//": double hybrid total = E_scf + pt2_cos*E_os + pt2_css*E_ss (literals) to round-off")
      rc2 = trc_destroy(h2);                         call expect(rc2 == TRC_OK, name//": destroy")
   end subroutine check_double_hybrid

   ! The auxiliary (or jkfit) basis must end up on the SAME Cartesian/
   ! spherical convention as the orbital basis regardless of which is
   ! set first (pyscf's own auxmol.cart = mol.cart rule) --
   ! trc_set_basis_json/trc_set_aux_json only inherit the convention
   ! FORWARD (aux from an already-loaded orbital basis), so an aux
   ! basis set FIRST had nothing to inherit from and could silently
   ! disagree once the orbital basis loaded after it. Two bases with
   ! genuinely different native tags exercise both orders: 6-31g.json
   ! has no spherical shell anywhere (always Cartesian, no
   ! force_cartesian needed), cc-pvdz-rifit.json is spherical
   ! throughout -- orbital-first inherits correctly (already exercised
   ! elsewhere in this program, confirmed again here); aux-first must
   ! now be refused rather than silently mismatched.
   subroutine check_aux_order()
      type(c_ptr) :: h3, h4
      integer(c_int) :: rc3

      ! Correct order: orbital first, aux second -- aux inherits the
      ! orbital's own Cartesian convention and succeeds.
      rc3 = trc_create(h3);                      call expect(rc3 == TRC_OK, "aux order: create (orbital-first)")
      rc3 = trc_set_molecule(h3, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(rc3 == TRC_OK, "aux order: set_molecule (orbital-first)")
      rc3 = trc_set_basis_json(h3, 'basis_sets/6-31g.json'//c_null_char)
      call expect(rc3 == TRC_OK, "aux order: orbital first succeeds")
      rc3 = trc_set_aux_json(h3, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      call expect(rc3 == TRC_OK, "aux order: aux after orbital inherits and succeeds")
      rc3 = trc_destroy(h3);                      call expect(rc3 == TRC_OK, "aux order: destroy (orbital-first)")

      ! Wrong order: aux first -- nothing to inherit from yet, so it
      ! keeps its own native (spherical) tag; orbital second keeps ITS
      ! own native (Cartesian) tag -- now refused instead of silently
      ! running RI-MP2 with two AO spaces that disagree.
      rc3 = trc_create(h4);                      call expect(rc3 == TRC_OK, "aux order: create (aux-first)")
      rc3 = trc_set_molecule(h4, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(rc3 == TRC_OK, "aux order: set_molecule (aux-first)")
      rc3 = trc_set_aux_json(h4, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      call expect(rc3 == TRC_OK, "aux order: aux loaded first, nothing to inherit from yet, succeeds")
      rc3 = trc_set_basis_json(h4, 'basis_sets/6-31g.json'//c_null_char)
      call expect(rc3 /= TRC_OK, "aux order: orbital after a mismatched aux is refused, not silently wrong")
      rc3 = trc_destroy(h4);                      call expect(rc3 == TRC_OK, "aux order: destroy (aux-first)")
   end subroutine check_aux_order

end program context_scf
