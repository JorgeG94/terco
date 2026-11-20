module check_capi_tddft_suite
   !!
   !! The C ABI's TDDFT entry points, called through the generated interface module (include/trc_c_interfaces.f90),
   !! so what a C caller sees is what runs: trc_set_tddft, trc_set_tddft_convergence, trc_check_tddft, trc_run_tddft,
   !! trc_tddft_nroots, trc_tddft_energies, trc_tddft_character and trc_tddft_energy, then trc_run_gradient.
   !!  1. Hartree-Fock TDA (CIS) and RPA (TDHF): the gradient of a chosen root of water 6-31G against pyscf's
   !!     grad.tdrhf (tddft_grad_water_6-31g_hf.ref, tddft_rpa_grad_water_6-31g_hf.ref), and the energies the
   !!     handle reports (E_SCF + omega of the target root);
   !!  2. water PBE: the gradient with trc_set_tddft's grid_response = 0 against pyscf's grad.tdrks
   !!     (tddft_grad_water_6-31g_pbe.ref, Pulay only), and the default, grid response, differing from it by the
   !!     grid's derivative (about 1e-5);
   !!  3. state: the roots are not there before trc_run_tddft, trc_run_scf drops them and the gradient is the SCF's
   !!     again, trc_run_hessian refuses after a TDDFT run;
   !!  4. what is refused, with TRC_ERR_UNSUPPORTED and a message: a triplet (trc_set_tddft), a meta-GGA, a
   !!     range-separated hybrid, a double hybrid, an open shell (trc_check_tddft), and a bad root (TRC_ERR_BADARG).
   !!  5. state rules: SCF, TDDFT, SCF, TDDFT gives identical roots and no roots in between; RI-MP2 after TDDFT drops the
   !!     roots and the gradient is RI-MP2's; an effective core potential (HI/def2-SVP) refuses the excited-state
   !!     gradient as the ground-state one does; TDDFT after trc_load_scf (no grid in the loaded result) gives the roots
   !!     of the SCF that was saved.
   !! Run from test/. A test-drive suite. The argument `pbe` runs only the water PBE gradient (2.), the rest only without
   !! it: two runs of one program, each around ten seconds.
   !!
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char, c_char, c_null_ptr
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix3
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_method, &
                               trc_set_convergence, trc_set_force_cartesian, trc_run_scf, trc_run_gradient, trc_gradient, &
                               trc_energy, trc_run_hessian, trc_message, trc_set_tddft, trc_set_tddft_convergence, &
                               trc_set_tddft_missed_root_iterations, trc_check_tddft, trc_run_tddft, trc_tddft_nroots, &
                               trc_tddft_energies, &
                               trc_tddft_character, trc_tddft_energy, trc_tddft_target_root, trc_abi_version, &
                               trc_set_aux_json, trc_set_rimp2, trc_run_rimp2, trc_save_scf, trc_load_scf, &
                               TRC_OK, TRC_ERR_UNSUPPORTED, TRC_ERR_BADARG, TRC_ERR_STATE, TRC_ERR_NULL, TRC_ABI_MINOR
   implicit none
   private
   public :: collect_capi_tddft, mode_arg

   character(len=16) :: mode_arg = ""

   real(dp), parameter :: GATE = 2.0e-8_dp
   !! Against pyscf at the C ABI's own SCF convergence (1e-12 energy, 1e-9 commutator).

contains

   subroutine collect_capi_tddft(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (trim(mode_arg) == "pbe") then
         testsuite = [new_unittest("water PBE: Pulay-only gradient vs pyscf, grid response is the default", test_pbe)]
         return
      end if
      testsuite = [ &
                  new_unittest("Hartree-Fock: TDA and RPA gradients of a root vs pyscf, energies", test_hf), &
                  new_unittest("state: before the run, after a new SCF, the Hessian", test_state), &
                  new_unittest("refusals", test_refusals), &
                  new_unittest("missed-root iterations: 24 is the default, 0 reaches the solver", &
                               test_missed_root), &
                  new_unittest("a second SCF and TDDFT give the same roots, none in between", test_rerun), &
                  new_unittest("RI-MP2 after TDDFT takes over the energy and the gradient", test_rimp2_after), &
                  new_unittest("an effective core potential: the excited-state gradient is refused", test_ecp), &
                  new_unittest("TDDFT on a loaded SCF (the grid is rebuilt)", test_loaded_scf) &
                  ]
   end subroutine collect_capi_tddft

   function message_of(handle) result(text)
      type(c_ptr), intent(in) :: handle
      character(len=:), allocatable :: text
      character(kind=c_char) :: buffer(512)
      integer(c_int) :: rc
      integer :: i
      text = ""
      rc = trc_message(handle, buffer, int(size(buffer), c_int))
      do i = 1, size(buffer)
         if (buffer(i) == c_null_char) exit
         text = text//buffer(i)
      end do
   end function message_of

   subroutine open_handle(error, handle, functional, basis, ok, xyz)
      !! A water (or `xyz`) handle with the basis, the method and the SCF settings of the pyscf references.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(out) :: handle
      character(len=*), intent(in) :: functional, basis
      logical, intent(out) :: ok
      character(len=*), intent(in), optional :: xyz
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      integer(c_int) :: rc

      if (present(xyz)) then
         call read_xyz(xyz, natm, zint, at_r)
      else
         call read_xyz("water.xyz", natm, zint, at_r)
      end if
      rc = trc_create(handle)
      if (rc == TRC_OK) rc = trc_set_molecule(handle, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                                              0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_force_cartesian(handle, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(handle, basis//c_null_char)
      if (rc == TRC_OK) rc = trc_set_method(handle, functional//c_null_char, 3_c_int)
      if (rc == TRC_OK) rc = trc_set_convergence(handle, 1.0e-12_c_double, 1.0e-9_c_double, 200_c_int, 10_c_int)
      ok = rc == TRC_OK
      call check(error, ok, "setting up the handle failed, status "//to_string(int(rc)))
   end subroutine open_handle

   subroutine gradient_of(error, handle, grad)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(in) :: handle
      real(c_double), intent(out) :: grad(9)
      integer(c_int) :: rc
      rc = trc_run_gradient(handle)
      call check(error, rc == TRC_OK, "trc_run_gradient: "//message_of(handle))
      if (allocated(error)) return
      rc = trc_gradient(handle, grad)
      call check(error, rc == TRC_OK, "trc_gradient failed")
   end subroutine gradient_of

   subroutine ref_gradient(name, root, gref, e_ref)
      character(len=*), intent(in) :: name
      integer, intent(in) :: root
      real(dp), intent(out) :: gref(3, 3), e_ref
      real(dp) :: all_roots(3, 3, 3)
      integer :: u
      call ref_require(name, u)
      call ref_read_real(u, "E_SCF", e_ref)
      call ref_read_matrix3(u, "GRAD", 3, 3, 3, all_roots)
      call ref_close(u)
      gref = all_roots(:, :, root)
   end subroutine ref_gradient

   subroutine hf_root(error, ref_name, tda, root)
      !! Hartree-Fock water/6-31G, three roots, TDA or RPA: the gradient of `root` against the reference.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: ref_name
      logical, intent(in) :: tda
      integer, intent(in) :: root
      type(c_ptr) :: handle
      real(c_double) :: grad(9), omega(3), osc(3), ynorm(3), dipole(9), e_total, e_scf, weight(3)
      real(dp) :: gref(3, 3), e_ref
      integer(c_int) :: rc, nroots, hole(3), particle(3)
      logical :: ok

      call open_handle(error, handle, "", "basis_sets/6-31g.json", ok)
      if (allocated(error)) return
      rc = trc_set_tddft(handle, 3_c_int, merge(1_c_int, 0_c_int, tda), 1_c_int, int(root, c_int), 1_c_int)
      if (rc == TRC_OK) rc = trc_set_tddft_convergence(handle, 1.0e-10_c_double, 0_c_int, 1.0e-11_c_double, 0_c_int, 0_c_int)
      call check(error, rc == TRC_OK, "trc_set_tddft: "//message_of(handle))
      if (.not. allocated(error)) then
         rc = trc_run_scf(handle)
         call check(error, rc == TRC_OK, "the SCF failed")
      end if
      if (.not. allocated(error)) then
         rc = trc_run_tddft(handle)
         call check(error, rc == TRC_OK, "trc_run_tddft: "//message_of(handle))
      end if
      if (.not. allocated(error)) then
         rc = trc_tddft_nroots(handle, nroots)
         call check(error, rc == TRC_OK .and. nroots == 3, "three roots were solved")
      end if
      if (.not. allocated(error)) then
         rc = trc_tddft_energies(handle, 3_c_int, omega, osc, ynorm, dipole)
         if (rc == TRC_OK) rc = trc_tddft_character(handle, 3_c_int, hole, particle, weight)
         if (rc == TRC_OK) rc = trc_tddft_energy(handle, e_total)
         if (rc == TRC_OK) rc = trc_energy(handle, e_scf)
         call check(error, rc == TRC_OK, "reading the roots failed")
      end if
      if (.not. allocated(error)) then
         call check(error, all(omega(2:3) > omega(1:2)), "the roots are ascending")
         if (.not. allocated(error)) call check(error, abs(e_total - (e_scf + omega(root))) < 1.0e-13_dp, &
                                                "trc_tddft_energy is E_SCF + omega of the target root")
         if (.not. allocated(error)) call check(error, merge(maxval(ynorm) == 0.0_dp, minval(ynorm) > 0.0_dp, tda), &
                                                "||Y||^2 is zero for the TDA and not for the RPA")
         if (.not. allocated(error)) call check(error, all(hole >= 1) .and. all(particle >= 1), "dominant pairs are 1-based")
      end if
      if (.not. allocated(error)) call gradient_of(error, handle, grad)
      if (.not. allocated(error)) then
         call ref_gradient(ref_name, root, gref, e_ref)
         call check(error, maxval(abs(reshape(grad, [3, 3]) - gref)) < GATE, "gradient of root "//to_string(root)// &
                    " vs pyscf, "//to_string(maxval(abs(reshape(grad, [3, 3]) - gref)), "es10.2"))
      end if
      rc = trc_destroy(handle)
   end subroutine hf_root

   subroutine test_hf(error)
      type(error_type), allocatable, intent(out) :: error
      call hf_root(error, "tddft_grad_water_6-31g_hf.ref", .true., 2)
      if (allocated(error)) return
      call hf_root(error, "tddft_rpa_grad_water_6-31g_hf.ref", .false., 3)
   end subroutine test_hf

   subroutine test_pbe(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      real(c_double) :: pulay(9), exact(9)
      real(dp) :: gref(3, 3), e_ref, diff
      integer(c_int) :: rc
      logical :: ok

      call open_handle(error, handle, "pbe", "basis_sets/6-31g.json", ok)
      if (allocated(error)) return
      rc = trc_set_tddft(handle, 3_c_int, 1_c_int, 1_c_int, 2_c_int, 0_c_int)
      if (rc == TRC_OK) rc = trc_set_tddft_convergence(handle, 1.0e-10_c_double, 0_c_int, 1.0e-11_c_double, 0_c_int, 0_c_int)
      if (rc == TRC_OK) rc = trc_run_scf(handle)
      if (rc == TRC_OK) rc = trc_run_tddft(handle)
      call check(error, rc == TRC_OK, "PBE setup, SCF or run: "//message_of(handle))
      if (.not. allocated(error)) call gradient_of(error, handle, pulay)
      if (.not. allocated(error)) then
         call ref_gradient("tddft_grad_water_6-31g_pbe.ref", 2, gref, e_ref)
         call check(error, maxval(abs(reshape(pulay, [3, 3]) - gref)) < GATE, "Pulay-only gradient vs pyscf, "// &
                    to_string(maxval(abs(reshape(pulay, [3, 3]) - gref)), "es10.2"))
      end if
      ! The default of trc_set_tddft's grid_response argument is its argument: ask for it.
      if (.not. allocated(error)) then
         rc = trc_set_tddft(handle, 3_c_int, 1_c_int, 1_c_int, 2_c_int, 1_c_int)
         if (rc == TRC_OK) rc = trc_set_tddft_convergence(handle, 1.0e-10_c_double, 0_c_int, 1.0e-11_c_double, 0_c_int, 0_c_int)
         if (rc == TRC_OK) rc = trc_run_tddft(handle)
         call check(error, rc == TRC_OK, "second run: "//message_of(handle))
      end if
      if (.not. allocated(error)) call gradient_of(error, handle, exact)
      if (.not. allocated(error)) then
         diff = maxval(abs(exact - pulay))
         call check(error, diff > 1.0e-6_dp .and. diff < 1.0e-4_dp, "grid response should move the gradient by about "// &
                    "1e-5, got "//to_string(diff, "es10.2"))
      end if
      rc = trc_destroy(handle)
   end subroutine test_pbe

   subroutine test_state(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      real(c_double) :: omega(2), osc(2), ynorm(2), dipole(6), g_excited(9), g_ground(9), e_total
      integer(c_int) :: rc, nroots, major, minor
      logical :: ok

      rc = trc_abi_version(major, minor)
      call check(error, rc == TRC_OK .and. major == 1 .and. minor >= 31, "the ABI is at least 1.31")
      if (allocated(error)) return
      call open_handle(error, handle, "", "basis_sets/sto-3g.json", ok)
      if (allocated(error)) return
      rc = trc_set_tddft(handle, 2_c_int, 1_c_int, 1_c_int, 1_c_int, 1_c_int)
      rc = trc_tddft_nroots(handle, nroots)
      call check(error, rc == TRC_ERR_STATE, "no roots before a run: status "//to_string(int(rc)))
      if (allocated(error)) return
      rc = trc_run_tddft(handle)
      call check(error, rc == TRC_ERR_STATE, "trc_run_tddft before the SCF is a state error")
      if (allocated(error)) return
      rc = trc_run_scf(handle)
      rc = trc_run_tddft(handle)
      call check(error, rc == TRC_OK, "trc_run_tddft: "//message_of(handle))
      if (.not. allocated(error)) call gradient_of(error, handle, g_excited)
      if (.not. allocated(error)) then
         rc = trc_run_hessian(handle)
         call check(error, rc == TRC_ERR_UNSUPPORTED .and. index(message_of(handle), "excited-state Hessian") > 0, &
                    "the Hessian is refused after a TDDFT run")
      end if
      ! A new SCF drops the roots and the gradient is the ground state's again.
      if (.not. allocated(error)) then
         rc = trc_run_scf(handle)
         rc = trc_tddft_energies(handle, 2_c_int, omega, osc, ynorm, dipole)
         call check(error, rc == TRC_ERR_STATE, "a new SCF drops the roots")
      end if
      if (.not. allocated(error)) call gradient_of(error, handle, g_ground)
      if (.not. allocated(error)) then
         call check(error, maxval(abs(g_excited - g_ground)) > 1.0e-3_dp, "the ground-state gradient is not the excited one")
      end if
      rc = trc_tddft_energy(handle, e_total)
      if (.not. allocated(error)) call check(error, rc == TRC_ERR_STATE, "no excited-state energy after a new SCF")
      rc = trc_destroy(handle)
   end subroutine test_state

   subroutine expect_unsupported(error, functional, tddft_spin, word)
      !! trc_check_tddft on a handle with `functional` (and the spin of trc_set_tddft) refuses with `word`.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: functional, word
      integer(c_int), intent(in) :: tddft_spin
      type(c_ptr) :: handle
      integer(c_int) :: rc
      logical :: ok

      call open_handle(error, handle, functional, "basis_sets/sto-3g.json", ok)
      if (allocated(error)) return
      rc = trc_set_tddft(handle, 2_c_int, 1_c_int, tddft_spin, 1_c_int, 1_c_int)
      if (tddft_spin == 0_c_int) then
         call check(error, rc == TRC_ERR_UNSUPPORTED .and. index(message_of(handle), word) > 0, &
                    "a triplet is refused by trc_set_tddft: "//message_of(handle))
      else
         rc = trc_check_tddft(handle, 1_c_int)
         call check(error, rc == TRC_ERR_UNSUPPORTED, "'"//functional//"' must be refused, status "//to_string(int(rc)))
         if (.not. allocated(error)) call check(error, index(message_of(handle), word) > 0, &
                                                "the message for '"//functional//"' says: "//message_of(handle))
      end if
      rc = trc_destroy(handle)
   end subroutine expect_unsupported

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: rc
      logical :: ok

      call expect_unsupported(error, "", 0_c_int, "triplet")
      if (allocated(error)) return
      call expect_unsupported(error, "tpss", 1_c_int, "meta-GGA")
      if (allocated(error)) return
      call expect_unsupported(error, "wb97x", 1_c_int, "range-separated")
      if (allocated(error)) return
      call expect_unsupported(error, "b2plyp", 1_c_int, "double hybrid")
      if (allocated(error)) return
      call open_handle(error, handle, "", "basis_sets/sto-3g.json", ok)
      if (allocated(error)) return
      rc = trc_set_tddft(handle, 2_c_int, 1_c_int, 1_c_int, 3_c_int, 1_c_int)
      call check(error, rc == TRC_ERR_BADARG, "a target root past nroots is a bad argument")
      if (.not. allocated(error)) then
         rc = trc_set_tddft(handle, 0_c_int, 1_c_int, 1_c_int, 1_c_int, 1_c_int)
         call check(error, rc == TRC_ERR_BADARG, "zero roots is a bad argument")
      end if
      ! A tolerance must be finite and not negative: a negative, NaN or infinite one is refused, and a refusal keeps the
      ! earlier settings.
      if (.not. allocated(error)) then
         rc = trc_set_tddft_convergence(handle, -1.0e-8_c_double, 0_c_int, 0.0_c_double, 0_c_int, 0_c_int)
         call check(error, rc == TRC_ERR_BADARG, "a negative conv_tol is a bad argument")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_tddft_convergence(handle, ieee_value(1.0_c_double, ieee_quiet_nan), 0_c_int, 0.0_c_double, 0_c_int, 0_c_int)
         call check(error, rc == TRC_ERR_BADARG, "a NaN conv_tol is a bad argument")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_tddft_convergence(handle, ieee_value(1.0_c_double, ieee_positive_inf), 0_c_int, 0.0_c_double, 0_c_int, &
                                        0_c_int)
         call check(error, rc == TRC_ERR_BADARG, "an infinite conv_tol is a bad argument")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_tddft_convergence(handle, 0.0_c_double, 0_c_int, ieee_value(1.0_c_double, ieee_positive_inf), 0_c_int, &
                                        0_c_int)
         call check(error, rc == TRC_ERR_BADARG, "an infinite cpks_tol is a bad argument")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_tddft_convergence(handle, 1.0e-9_c_double, 0_c_int, 1.0e-9_c_double, 0_c_int, 0_c_int)
         call check(error, rc == TRC_OK, "an ordinary pair of tolerances is accepted")
      end if
      ! The missed-root budget: a negative count is refused with a message, a null handle is TRC_ERR_NULL.
      if (.not. allocated(error)) then
         rc = trc_set_tddft_missed_root_iterations(handle, -1_c_int)
         call check(error, rc == TRC_ERR_BADARG .and. index(message_of(handle), "negative") > 0, &
                    "a negative missed-root iteration count is a bad argument with a message")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_tddft_missed_root_iterations(handle, 0_c_int)
         call check(error, rc == TRC_OK, "a missed-root iteration count of zero is accepted")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_tddft_missed_root_iterations(c_null_ptr, 0_c_int)
         call check(error, rc == TRC_ERR_NULL, "a null handle is TRC_ERR_NULL")
      end if
      rc = trc_destroy(handle)
   end subroutine test_refusals

   subroutine third_root(error, handle, tda, omega3)
      !! Runs the TDDFT of the handle's SCF (three roots, 1e-7) and returns the third root's energy.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(in) :: handle
      logical, intent(in) :: tda
      real(c_double), intent(out) :: omega3
      real(c_double) :: omega(3)
      integer(c_int) :: rc
      logical :: ok

      omega3 = 0.0_c_double
      rc = trc_set_tddft(handle, 3_c_int, merge(1_c_int, 0_c_int, tda), 1_c_int, 1_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_tddft_convergence(handle, 1.0e-7_c_double, 0_c_int, 0.0_c_double, 0_c_int, 0_c_int)
      if (rc == TRC_OK) rc = trc_run_tddft(handle)
      call check(error, rc == TRC_OK, "trc_run_tddft: "//message_of(handle))
      if (allocated(error)) return
      call roots_of(handle, 3, omega, ok)
      call check(error, ok, "the roots could not be read")
      omega3 = omega(3)
   end subroutine third_root

   subroutine test_missed_root(error)
      !! Formaldehyde 6-31+G Hartree-Fock, TDA and RPA: the third root is reached only by iterating the candidate
      !! pairs after the probe. Through the C ABI a count of 24 is the default (same root), and 0 loses that root.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      real(c_double) :: w_24, w_0
      integer(c_int) :: rc
      integer :: pass
      logical :: ok, tda

      call open_handle(error, handle, "", "../basis_sets/6-31+g.json", ok, xyz="formaldehyde.xyz")
      if (allocated(error)) return
      rc = trc_run_scf(handle)
      call check(error, rc == TRC_OK, "trc_run_scf: "//message_of(handle))
      do pass = 1, 2
         if (allocated(error)) exit
         tda = pass == 1
         rc = trc_set_tddft_missed_root_iterations(handle, 24_c_int)
         call check(error, rc == TRC_OK, "trc_set_tddft_missed_root_iterations(24)")
         if (allocated(error)) exit
         ! Set before the default run, so this is the explicit 24; a fresh handle's default is the same number.
         call third_root(error, handle, tda, w_24)
         if (allocated(error)) exit
         rc = trc_set_tddft_missed_root_iterations(handle, 0_c_int)
         call check(error, rc == TRC_OK, "trc_set_tddft_missed_root_iterations(0)")
         if (allocated(error)) exit
         call third_root(error, handle, tda, w_0)
         if (allocated(error)) exit
         call check(error, abs(w_0 - w_24) > 1.0e-3_dp, merge("TDA", "RPA", tda)// &
                    ": a count of 0 left the third root alone, so the setting reached no solver")
         if (allocated(error)) exit
         ! A second handle never given the setting must equal the explicit 24.
         call check_default_handle(error, tda, w_24)
      end do
      rc = trc_destroy(handle)
   end subroutine test_missed_root

   subroutine check_default_handle(error, tda, w_24)
      !! A handle that never sets the count has the explicit 24's third root.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: tda
      real(c_double), intent(in) :: w_24
      type(c_ptr) :: handle
      real(c_double) :: w
      integer(c_int) :: rc
      logical :: ok

      call open_handle(error, handle, "", "../basis_sets/6-31+g.json", ok, xyz="formaldehyde.xyz")
      if (allocated(error)) return
      rc = trc_run_scf(handle)
      call check(error, rc == TRC_OK, "trc_run_scf: "//message_of(handle))
      if (.not. allocated(error)) call third_root(error, handle, tda, w)
      if (.not. allocated(error)) call check(error, abs(w - w_24) < 1.0e-9_dp, &
                                             merge("TDA", "RPA", tda)//": the default is not a count of 24")
      rc = trc_destroy(handle)
   end subroutine check_default_handle

   subroutine roots_of(handle, nroots, omega, ok)
      type(c_ptr), intent(in) :: handle
      integer, intent(in) :: nroots
      real(c_double), intent(out) :: omega(nroots)
      logical, intent(out) :: ok
      real(c_double) :: osc(nroots), ynorm(nroots), dipole(3*nroots)
      ok = trc_tddft_energies(handle, int(nroots, c_int), omega, osc, ynorm, dipole) == TRC_OK
   end subroutine roots_of

   subroutine test_rerun(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      real(c_double) :: first(2), second(2), osc(2), ynorm(2), dipole(6), small(1)
      integer(c_int) :: rc, root, nroots
      logical :: ok

      call open_handle(error, handle, "", "basis_sets/sto-3g.json", ok)
      if (allocated(error)) return
      rc = trc_set_tddft(handle, 2_c_int, 1_c_int, 1_c_int, 2_c_int, 1_c_int)
      rc = trc_tddft_target_root(handle, root, nroots)
      call check(error, rc == TRC_OK .and. root == 2 .and. nroots == 2, "the target root is read from the context")
      if (.not. allocated(error)) rc = trc_run_scf(handle)
      if (.not. allocated(error)) rc = trc_run_tddft(handle)
      call check(error, rc == TRC_OK, "first run: "//message_of(handle))
      if (.not. allocated(error)) call roots_of(handle, 2, first, ok)
      if (.not. allocated(error)) then
         rc = trc_tddft_energies(handle, 1_c_int, small, osc, ynorm, dipole)
         call check(error, rc == TRC_ERR_BADARG, "an array that holds fewer roots than were solved is refused")
      end if
      if (.not. allocated(error)) then
         rc = trc_run_scf(handle)
         rc = trc_tddft_nroots(handle, nroots)
         call check(error, rc == TRC_ERR_STATE, "between the two runs the getters report none")
      end if
      if (.not. allocated(error)) then
         rc = trc_tddft_energies(handle, 2_c_int, second, osc, ynorm, dipole)
         call check(error, rc == TRC_ERR_STATE, "between the two runs the energies report none")
      end if
      if (.not. allocated(error)) rc = trc_run_tddft(handle)
      if (.not. allocated(error)) call roots_of(handle, 2, second, ok)
      if (.not. allocated(error)) then
         call check(error, ok .and. maxval(abs(first - second)) < 1.0e-9_dp, "the second run's roots equal the first's, "// &
                    to_string(maxval(abs(first - second)), "es10.2"))
      end if
      rc = trc_destroy(handle)
   end subroutine test_rerun

   subroutine test_rimp2_after(error)
      !! SCF -> TDDFT -> RI-MP2: the excited state is gone, whatever answers trc_run_gradient is RI-MP2's.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      real(c_double) :: g_excited(9), g_mp2(9), e_total
      integer(c_int) :: rc, nroots
      logical :: ok

      call open_handle(error, handle, "", "basis_sets/sto-3g.json", ok)
      if (allocated(error)) return
      rc = trc_set_aux_json(handle, "basis_sets/cc-pvdz-rifit.json"//c_null_char)
      if (rc == TRC_OK) rc = trc_set_rimp2(handle, -1_c_int, 0_c_int)
      if (rc == TRC_OK) rc = trc_set_tddft(handle, 2_c_int, 1_c_int, 1_c_int, 1_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_run_scf(handle)
      if (rc == TRC_OK) rc = trc_run_tddft(handle)
      call check(error, rc == TRC_OK, "setup: "//message_of(handle))
      if (.not. allocated(error)) call gradient_of(error, handle, g_excited)
      if (.not. allocated(error)) then
         rc = trc_run_rimp2(handle)
         call check(error, rc == TRC_OK, "trc_run_rimp2: "//message_of(handle))
      end if
      if (.not. allocated(error)) then
         rc = trc_tddft_nroots(handle, nroots)
         call check(error, rc == TRC_ERR_STATE, "RI-MP2 drops the TDDFT roots")
         rc = trc_tddft_energy(handle, e_total)
         if (.not. allocated(error)) call check(error, rc == TRC_ERR_STATE, "and the excited-state energy")
      end if
      if (.not. allocated(error)) call gradient_of(error, handle, g_mp2)
      if (.not. allocated(error)) then
         call check(error, maxval(abs(g_mp2 - g_excited)) > 1.0e-4_dp, "the gradient is RI-MP2's, not the excited state's")
      end if
      rc = trc_destroy(handle)
   end subroutine test_rimp2_after

   subroutine test_ecp(error)
      !! HI with the def2-SVP effective core potential: the SCF and the roots run, the gradient is refused (no U_ECP
      !! derivative), exactly as the ground-state gradient is.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      real(c_double), parameter :: Z(2) = [1.0_c_double, 53.0_c_double]
      real(c_double), parameter :: R(6) = [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, &
                                           3.0407_c_double]
      integer(c_int) :: rc

      rc = trc_create(handle)
      if (rc == TRC_OK) rc = trc_set_molecule(handle, 2_c_int, Z, R, 0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_force_cartesian(handle, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(handle, "basis_sets/def2-svp.json"//c_null_char)
      if (rc == TRC_OK) rc = trc_set_tddft(handle, 2_c_int, 1_c_int, 1_c_int, 1_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_run_scf(handle)
      if (rc == TRC_OK) rc = trc_run_tddft(handle)
      call check(error, rc == TRC_OK, "HI SCF and roots: "//message_of(handle))
      if (.not. allocated(error)) then
         rc = trc_run_gradient(handle)
         call check(error, rc == TRC_ERR_UNSUPPORTED, "the gradient over an ECP must be refused, status "//to_string(int(rc)))
      end if
      if (.not. allocated(error)) then
         call check(error, index(message_of(handle), "core-potential") > 0, "the refusal names the ECP: "//message_of(handle))
      end if
      rc = trc_destroy(handle)
   end subroutine test_ecp

   subroutine test_loaded_scf(error)
      !! A result read back from a file has no grid: trc_run_tddft rebuilds it, and the roots are those of the SCF saved.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: saved, loaded
      character(len=*), parameter :: PATH = "tddft_capi_restart.bin"
      real(c_double) :: before(2), after(2)
      integer(c_int) :: rc
      integer :: u
      logical :: ok

      call open_handle(error, saved, "pbe", "basis_sets/sto-3g.json", ok)
      if (allocated(error)) return
      rc = trc_set_tddft(saved, 2_c_int, 1_c_int, 1_c_int, 1_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_run_scf(saved)
      if (rc == TRC_OK) rc = trc_run_tddft(saved)
      if (rc == TRC_OK) rc = trc_save_scf(saved, PATH//c_null_char)
      call check(error, rc == TRC_OK, "the SCF to save: "//message_of(saved))
      if (.not. allocated(error)) call roots_of(saved, 2, before, ok)
      rc = trc_destroy(saved)
      if (allocated(error)) return
      call open_handle(error, loaded, "pbe", "basis_sets/sto-3g.json", ok)
      if (allocated(error)) return
      rc = trc_set_tddft(loaded, 2_c_int, 1_c_int, 1_c_int, 1_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_load_scf(loaded, PATH//c_null_char)
      call check(error, rc == TRC_OK, "trc_load_scf: "//message_of(loaded))
      if (.not. allocated(error)) then
         rc = trc_run_tddft(loaded)
         call check(error, rc == TRC_OK, "trc_run_tddft on the loaded SCF: "//message_of(loaded))
      end if
      if (.not. allocated(error)) then
         call roots_of(loaded, 2, after, ok)
         call check(error, ok .and. maxval(abs(after - before)) < 1.0e-7_dp, "roots of the loaded SCF vs the saved one, "// &
                    to_string(maxval(abs(after - before)), "es10.2"))
      end if
      rc = trc_destroy(loaded)
      open (newunit=u, file=PATH, status="old", iostat=rc)
      if (rc == 0) close (u, status="delete")
   end subroutine test_loaded_scf

end module check_capi_tddft_suite

program check_capi_tddft
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_capi_tddft_suite, only: collect_capi_tddft, mode_arg
   implicit none
   integer :: stat

   stat = 0
   if (command_argument_count() >= 1) call get_command_argument(1, mode_arg)
   call run_testsuite(collect_capi_tddft, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_tddft
