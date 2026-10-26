module check_capi_shells_suite
   !! trc_basis_shell_count/trc_basis_shell_max_primitives/
   !! trc_basis_shell_info/trc_basis_shell_primitives/
   !! trc_basis_is_spherical and trc_physical_constant, through
   !! include/trc_c_interfaces.f90 -- cyclone's own motivating gaps (it
   !! could previously only check a GAMESS-to-terco basis conversion
   !! through the AO count and the final energy, and kept a second copy
   !! of CODATA 2018 that could drift from terco's).
   !!
   !! water/6-31G* (test/basis_sets/6-31gs.json), through
   !! trc_set_basis_json, checks the shells terco actually built against
   !! that file's own literal exponents -- l, nprim and the atom map come
   !! straight from the file's shell list, and exponents are untouched by
   !! trc_basis_json's own renormalisation (trc_basis_shell_primitives'
   !! own doc), so they are the file's literal numbers. Coefficients on
   !! THAT path are NOT: a separate case builds a basis directly through
   !! trc_basis_create, where trc_basis_shell_primitives' own doc says
   !! the coefficient convention IS the caller's input, and checks that
   !! claim against basis_sets/6-31g.json's literal oxygen 1s shell.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_set_molecule, trc_set_basis_json, &
                               trc_context_basis, trc_basis_create, trc_basis_destroy, &
                               trc_basis_shell_count, trc_basis_shell_max_primitives, &
                               trc_basis_shell_info, trc_basis_shell_primitives, &
                               trc_basis_is_spherical, trc_physical_constant, &
                               TRC_OK, TRC_ERR_BADARG, &
                               TRC_CONST_BOHR_TO_ANGSTROM, TRC_CONST_HARTREE_TO_EV, &
                               TRC_CONST_HARTREE_TO_KCALMOL, TRC_CONST_AMU_TO_ME, &
                               TRC_CONST_KB_SI, TRC_CONST_NA_SI
   implicit none
   private
   public :: collect_capi_shells

contains

   subroutine collect_capi_shells(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water/6-31G* shell structure matches the JSON file's own shell list", &
                               test_json_shell_structure), &
                  new_unittest("water/6-31G* exponents match the JSON file's own literals", &
                               test_json_exponents), &
                  new_unittest("water/6-31G* is Cartesian (trc_basis_is_spherical)", test_json_cartesian), &
                  new_unittest("a basis built through trc_basis_create round-trips its own coefficients", &
                               test_direct_coefficients), &
                  new_unittest("trc_physical_constant matches trc_constants and CODATA 2018", &
                               test_physical_constants), &
                  new_unittest("trc_physical_constant refuses an unknown id", test_unknown_constant) &
                  ]
   end subroutine collect_capi_shells

   !
   ! A water/6-31G* context, through the ABI, with its basis handle
   ! borrowed back out -- what every case below needs.
   !
   subroutine water_631gs_basis(basis)
      type(c_ptr), intent(out) :: basis
      type(c_ptr) :: h
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      integer(c_int) :: rc
      call read_xyz('water.xyz', natm, zint, at_r)
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), &
                            reshape(real(at_r, c_double), [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_context_basis(h, basis)
      ! `h` is left open: `basis` is BORROWED (trc_context_basis' own
      ! doc) and stays valid only while its context lives. The context
      ! leaks for the life of this check program, same as every other
      ! capi check that never calls trc_destroy on every path (e.g.
      ! check_capi_libcint_sph's early-return cases).
   end subroutine water_631gs_basis

   subroutine test_json_shell_structure(error)
      !! basis_sets/6-31gs.json's own shell list, water atom order O/H/H
      !! (test/water.xyz): O is 1s(6prim), then an SP shell(3prim) split
      !! into s/p, then an SP shell(1prim) split into s/p, then a
      !! Cartesian d(1prim) ("gto_cartesian" in the file); each H is
      !! 1s(3prim), 1s(1prim), no polarisation function.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: basis
      integer(c_int) :: rc, nshell, maxnp
      integer(c_int) :: atom(10), l(10), nprim(10)
      integer(c_int), parameter :: WANT_ATOM(10) = [1, 1, 1, 1, 1, 1, 2, 2, 3, 3]
      integer(c_int), parameter :: WANT_L(10) = [0, 0, 1, 0, 1, 2, 0, 0, 0, 0]
      integer(c_int), parameter :: WANT_NP(10) = [6, 3, 3, 1, 1, 1, 3, 1, 3, 1]

      call water_631gs_basis(basis)

      rc = trc_basis_shell_count(basis, nshell)
      call check(error, rc == TRC_OK, "trc_basis_shell_count status")
      if (allocated(error)) return
      call check(error, nshell == 10_c_int, "water/6-31G* has 10 shells")
      if (allocated(error)) return

      rc = trc_basis_shell_max_primitives(basis, maxnp)
      call check(error, rc == TRC_OK, "trc_basis_shell_max_primitives status")
      if (allocated(error)) return
      call check(error, maxnp == 6_c_int, "the widest shell (oxygen's 1s) has 6 primitives")
      if (allocated(error)) return

      rc = trc_basis_shell_info(basis, atom, l, nprim)
      call check(error, rc == TRC_OK, "trc_basis_shell_info status")
      if (allocated(error)) return
      call check(error, all(atom == WANT_ATOM), "shell -> atom map matches the file's own atom order")
      if (allocated(error)) return
      call check(error, all(l == WANT_L), "angular momenta match the file's own shell list")
      if (allocated(error)) return
      call check(error, all(nprim == WANT_NP), "primitive counts match the file's own shell list")
   end subroutine test_json_shell_structure

   subroutine test_json_exponents(error)
      !! basis_sets/6-31gs.json's own literals for oxygen's 1s shell
      !! (shell 1) and its polarisation d shell (shell 6): exponents are
      !! untouched by trc_basis_json's renormalisation, so these must be
      !! the file's numbers exactly (to the precision BSE publishes
      !! them at).
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: basis
      integer(c_int) :: rc
      real(c_double) :: exps(6, 10), coefs(6, 10)

      call water_631gs_basis(basis)
      rc = trc_basis_shell_primitives(basis, exps, coefs)
      call check(error, rc == TRC_OK, "trc_basis_shell_primitives status")
      if (allocated(error)) return

      call check(error, real(exps(1, 1), dp), 5484.671660_dp, thr=1.0e-4_dp, rel=.true.)
      if (allocated(error)) return
      call check(error, real(exps(6, 1), dp), 5.799635340_dp, thr=1.0e-6_dp, rel=.true.)
      if (allocated(error)) return
      call check(error, real(exps(1, 6), dp), 0.8_dp, thr=1.0e-12_dp, rel=.true.)
   end subroutine test_json_exponents

   subroutine test_json_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: basis
      integer(c_int) :: rc, spherical
      call water_631gs_basis(basis)
      rc = trc_basis_is_spherical(basis, spherical)
      call check(error, rc == TRC_OK, "trc_basis_is_spherical status")
      if (allocated(error)) return
      call check(error, spherical == 0_c_int, "6-31G* has no gto_spherical tag, so this basis is Cartesian")
   end subroutine test_json_cartesian

   subroutine test_direct_coefficients(error)
      !! basis_sets/6-31g.json's own oxygen 1s shell, BSE's literal
      !! exponents and coefficients, built directly through
      !! trc_basis_create (no JSON, no renormalisation): what
      !! trc_basis_shell_primitives hands back must be EXACTLY the
      !! numbers this test handed in, since basis_build's own
      !! common_fac_sp fold is undone exactly (a multiply by a
      !! compile-time constant, inverted the same way
      !! test/trc_test_basis.F90's unpack_basis already does). This is
      !! the path cyclone's GAMESS-to-terco conversion check is in.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: basis
      integer(c_int) :: rc
      integer(c_int), parameter :: L(1) = [0_c_int], NPRIM(1) = [6_c_int]
      real(c_double), parameter :: EXPS_IN(6) = [5484.671660_c_double, 825.2349460_c_double, &
                                                  188.0469580_c_double, 52.96450000_c_double, &
                                                  16.89757040_c_double, 5.799635340_c_double]
      real(c_double), parameter :: COEFS_IN(6) = [0.0018310744_c_double, 0.0139501722_c_double, &
                                                   0.0684450781_c_double, 0.2327143360_c_double, &
                                                   0.4701928980_c_double, 0.3585208530_c_double]
      real(c_double) :: centres(3), zatm(1), ratm(3)
      real(c_double) :: exps_out(6), coefs_out(6)

      centres = 0.0_c_double
      zatm = 8.0_c_double
      ratm = 0.0_c_double
      rc = trc_basis_create(1_c_int, 6_c_int, L, NPRIM, EXPS_IN, COEFS_IN, centres, &
                            1_c_int, zatm, ratm, basis)
      call check(error, rc == TRC_OK, "trc_basis_create status")
      if (allocated(error)) return

      rc = trc_basis_shell_primitives(basis, exps_out, coefs_out)
      call check(error, rc == TRC_OK, "trc_basis_shell_primitives status")
      if (allocated(error)) return
      call check(error, real(maxval(abs(exps_out - EXPS_IN)), dp), 0.0_dp, thr=1.0e-10_dp)
      if (allocated(error)) return
      call check(error, real(maxval(abs(coefs_out - COEFS_IN)), dp), 0.0_dp, thr=1.0e-10_dp)
      if (allocated(error)) return

      rc = trc_basis_destroy(basis)
      call check(error, rc == TRC_OK, "trc_basis_destroy status")
   end subroutine test_direct_coefficients

   subroutine test_physical_constants(error)
      !! Every selector against trc_constants' own named parameter AND
      !! the CODATA 2018 literal (physics.nist.gov/cuu/Constants, the
      !! same source AGENTS.md and check_constants.F90's own
      !! test_base_literals cite) -- an identity against trc_constants
      !! alone cannot catch a typo in the TRC_CONST_* case that happens
      !! to read the wrong named parameter.
      use trc_constants, only: BOHR_TO_ANGSTROM, HARTREE_TO_EV, HARTREE_TO_KCALMOL, AMU_TO_ME, &
                               KB_SI, NA_SI
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc
      real(c_double) :: value

      rc = trc_physical_constant(TRC_CONST_BOHR_TO_ANGSTROM, value)
      call check(error, rc == TRC_OK, "trc_physical_constant(BOHR_TO_ANGSTROM) status")
      if (allocated(error)) return
      call check(error, real(value, dp), BOHR_TO_ANGSTROM, thr=1.0e-13_dp)
      if (allocated(error)) return
      call check(error, real(value, dp), 0.529177210903_dp, thr=1.0e-13_dp)
      if (allocated(error)) return

      rc = trc_physical_constant(TRC_CONST_HARTREE_TO_EV, value)
      call check(error, rc == TRC_OK, "trc_physical_constant(HARTREE_TO_EV) status")
      if (allocated(error)) return
      call check(error, real(value, dp), HARTREE_TO_EV, thr=1.0e-10_dp)
      if (allocated(error)) return
      call check(error, real(value, dp), 27.211386245988_dp, thr=1.0e-10_dp)
      if (allocated(error)) return

      rc = trc_physical_constant(TRC_CONST_HARTREE_TO_KCALMOL, value)
      call check(error, rc == TRC_OK, "trc_physical_constant(HARTREE_TO_KCALMOL) status")
      if (allocated(error)) return
      call check(error, real(value, dp), HARTREE_TO_KCALMOL, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, real(value, dp), 627.5094740631_dp, thr=1.0e-9_dp)
      if (allocated(error)) return

      rc = trc_physical_constant(TRC_CONST_AMU_TO_ME, value)
      call check(error, rc == TRC_OK, "trc_physical_constant(AMU_TO_ME) status")
      if (allocated(error)) return
      call check(error, real(value, dp), AMU_TO_ME, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, real(value, dp), 1822.888486209_dp, thr=1.0e-9_dp)
      if (allocated(error)) return

      rc = trc_physical_constant(TRC_CONST_KB_SI, value)
      call check(error, rc == TRC_OK, "trc_physical_constant(KB_SI) status")
      if (allocated(error)) return
      call check(error, real(value, dp), KB_SI, thr=1.0e-30_dp)
      if (allocated(error)) return
      call check(error, real(value, dp), 1.380649e-23_dp, thr=1.0e-30_dp)
      if (allocated(error)) return

      rc = trc_physical_constant(TRC_CONST_NA_SI, value)
      call check(error, rc == TRC_OK, "trc_physical_constant(NA_SI) status")
      if (allocated(error)) return
      call check(error, real(value, dp), NA_SI, thr=10.0_dp)
      if (allocated(error)) return
      call check(error, real(value, dp), 6.02214076e23_dp, thr=10.0_dp)
   end subroutine test_physical_constants

   subroutine test_unknown_constant(error)
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc
      real(c_double) :: value
      value = -1.0_c_double
      rc = trc_physical_constant(999_c_int, value)
      call check(error, rc == TRC_ERR_BADARG, "an unknown id is refused with TRC_ERR_BADARG")
      if (allocated(error)) return
      call check(error, value == 0.0_c_double, "and value is left at 0")
   end subroutine test_unknown_constant

end module check_capi_shells_suite

program check_capi_shells
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_api, only: trc_bind_device
   use check_capi_shells_suite, only: collect_capi_shells
   implicit none
   integer :: stat, dev
   stat = 0
   dev = trc_bind_device(0)
   call run_testsuite(collect_capi_shells, error_unit, stat)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_shells
