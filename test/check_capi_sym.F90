!
! The symmetry entries of the C ABI (ABI 1.23: trc_set_symmetry, trc_set_irrep_nelec,
! trc_set_symmetry_skeleton, the report-frame pair and the getters), called through the
! Fortran interfaces of include/trc_c_interfaces.f90. The C side of the same round trip
! (water, the header compiled by a C compiler) is in link_c.c.
!
! Symmetry is opt-in: a fresh context is symmetry off (C1, the standard path), and
! "default" is a name for off. "on" (alias "auto") is the LENIENT request; a named group
! (and any request with fixed irrep occupations) is STRICT.
!
! HOCl (planar, Cs) 6-31G*, STRICT requests (the group name Cs). The oracle is the project's
! central one: the symmetric run and a plain run on the standard-frame geometry the context
! reports (trc_geometry) give the same energy, to 1e-10 Eh. Then the frame (an
! orthogonal map from the input to that geometry), the labels, fixed irrep
! occupations (and the occupations trc_mo_occupations reports for them), and what
! is refused, each with its message: every entry that would otherwise run the
! C1 path silently, the cases the SCF cannot do with symmetry, and an explicit
! group that is not a subgroup or not supported yet.
!
! Lenient requests: test_fresh_is_off (a fresh context and "default" run C1 and equal
! "off"; "on" and "auto" are the same request). With "on", a rotated water reports in
! the caller's INPUT frame and runs the paths the symmetric SCF cannot do in C1.
! test_default_round_trip: a rotated water, an "on" context against one with
! symmetry off: energy 1e-10, gradient 1e-9, orbital energies, density and the
! density of the rotated-back orbitals 1e-7, the reported geometry (the input
! symmetrised in it), atomic charges; the report-frame setter and getter.
! test_default_borrowed: the borrowed basis, its overlap and the dipole computed
! from it. test_default_fallbacks:
! DF, PCM, SOSCF, a given guess, and every entry the symmetric SCF cannot serve
! (Hessian, polarizability, RI-MP2, Fukui indices, restarts) give the symmetry-off
! numbers. test_irrep_nelec: irrep_nelec needs a symmetry request (a fresh context
! and an explicit off refuse); with "on" it makes the request strict.
! Edge cases: atoms with the same element but other shells (trc_set_basis_arrays)
! are not equivalent; a projected guess after a symmetric run starts from
! input-frame orbitals; irrep_nelec is strict in either order with the request;
! after an entry reruns the SCF in C1, a DFT gradient uses the input axes' grid.
!
module check_capi_sym_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char, c_null_ptr, c_loc
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_method, &
                               trc_set_convergence, trc_run_scf, trc_run_gradient, trc_gradient, trc_energy, &
                               trc_set_symmetry, trc_set_irrep_nelec, trc_set_symmetry_skeleton, trc_symmetry_group, trc_symmetry_detected, &
                               trc_symmetry_residual, trc_symmetry_frame, trc_symmetry_nirrep, &
                               trc_symmetry_irrep_name, trc_symmetry_irrep_dim, trc_mo_irreps, trc_geometry, &
                               trc_mo_occupations, trc_nao, trc_message, trc_set_geometry, trc_set_df, trc_set_pcm, &
                               trc_set_soscf, trc_run_rimp2, trc_run_hessian, trc_run_polarizability, trc_load_scf, &
                               trc_fukui_hirshfeld, trc_check_rimp2_gradient, trc_set_guess, &
                               trc_symmetry_report_frame, trc_set_symmetry_report_frame, trc_density, &
                               trc_mo_coefficients, trc_mo_energies, trc_nmo, trc_hirshfeld_charges, trc_chelpg_charges, &
                               trc_context_basis, trc_pairs_create, trc_pairs_destroy, trc_compute_1e, &
                               trc_compute_multipoles, trc_set_force_cartesian, trc_set_aux_json, trc_hessian, &
                               trc_polarizability, trc_rimp2_energy, trc_set_jkfit_json, trc_save_scf, &
                               trc_set_basis_arrays, trc_iterations, TRC_GUESS_PROJECTED, &
                               TRC_OK, TRC_ERR_BADARG, TRC_ERR_STATE, TRC_ERR_UNSUPPORTED, TRC_GUESS_GIVEN
   implicit none
   private
   public :: collect_capi_sym

   real(dp), parameter :: TOL_E = 1.0e-10_dp

contains

   subroutine collect_capi_sym(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("getters have nothing before a symmetric run", test_state), &
                  new_unittest("round trip: energy vs C1 on the reported geometry, frame, labels", test_round_trip), &
                  new_unittest("a second run reuses the frame; trc_set_geometry starts over", test_rerun), &
                  new_unittest("trc_set_symmetry_skeleton: off gives the same SCF and gradient", test_skeleton_switch), &
                  new_unittest("irrep_nelec: aufbau counts, occupations, a gradient refused for others", test_irrep_nelec), &
                  new_unittest("explicit groups: not a subgroup, a subgroup of the detected group, C1", test_explicit), &
                  new_unittest("the unsupported SCF settings are refused, not run in C1", test_scf_refusals), &
                  new_unittest("entries that would run C1 silently are refused", test_entry_refusals), &
                  new_unittest("on: a rotated water equals symmetry off in the input frame", test_default_round_trip), &
                  new_unittest("on: the borrowed basis, the overlap and the dipole are in the input frame", &
                               test_default_borrowed), &
                  new_unittest("on: the unsupported paths run in C1 with the symmetry-off numbers", &
                               test_default_fallbacks), &
                  new_unittest("fresh context: symmetry is off, and default/off/on/auto are the documented requests", &
                               test_fresh_is_off), &
                  new_unittest("per-atom shells: atoms with other shells are not equivalent", test_per_atom_shells), &
                  new_unittest("a projected guess after on projects input-frame orbitals", test_projected_after_on), &
                  new_unittest("irrep_nelec makes on strict in either order; cleared, lenient again", &
                               test_irrep_nelec_order), &
                  new_unittest("after a C1 rerun the gradient's grid is the input axes'", test_grid_after_c1_rerun) &
                  ]
   end subroutine collect_capi_sym

   subroutine make_context(h, xyz_file, basis, geom)
      !! A context on `xyz_file` (test/) with `basis`; `geom` is the input geometry in Bohr.
      type(c_ptr), intent(out) :: h
      character(len=*), intent(in) :: xyz_file, basis
      real(dp), allocatable, intent(out) :: geom(:, :)
      integer, allocatable :: elem(:)
      integer :: nat, rc

      call read_xyz(xyz_file, nat, elem, geom)
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(nat, c_int), real(elem, c_double), reshape(geom, [3*nat]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, basis//c_null_char)
      rc = trc_set_method(h, c_null_char, 3_c_int)
      rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 100_c_int, 10_c_int)
   end subroutine make_context

   function message(h) result(text)
      type(c_ptr), intent(in) :: h
      character(len=:), allocatable :: text
      character(kind=c_char) :: buf(512)
      integer(c_int) :: rc
      integer :: i
      rc = trc_message(h, buf, 512_c_int)
      text = ''
      do i = 1, 512
         if (buf(i) == c_null_char) exit
         text = text//buf(i)
      end do
   end function message

   function name_of(h, what, irrep) result(text)
      !! A text getter: 1 group used, 2 group detected, 3 irrep name.
      type(c_ptr), intent(in) :: h
      integer, intent(in) :: what, irrep
      character(len=:), allocatable :: text
      character(kind=c_char) :: buf(16)
      integer(c_int) :: rc
      integer :: i
      buf = c_null_char
      select case (what)
      case (1)
         rc = trc_symmetry_group(h, buf, 16_c_int)
      case (2)
         rc = trc_symmetry_detected(h, buf, 16_c_int)
      case default
         rc = trc_symmetry_irrep_name(h, int(irrep, c_int), buf, 16_c_int)
      end select
      text = ''
      if (rc /= TRC_OK) then
         text = '<status '//to_string(int(rc))//'>'
         return
      end if
      do i = 1, 16
         if (buf(i) == c_null_char) exit
         text = text//buf(i)
      end do
   end function name_of

   subroutine test_state(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: geom(:, :)
      real(c_double) :: resid
      integer(c_int) :: rc, n

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      call check(error, rc, TRC_OK, 'trc_set_symmetry(Cs)')
      if (allocated(error)) return
      call check(error, name_of(h, 1, 0), '<status 5>', 'no group before the SCF (TRC_ERR_STATE)')
      if (allocated(error)) return
      rc = trc_symmetry_residual(h, resid)
      call check(error, rc, TRC_ERR_STATE, 'no residual before the SCF')
      rc = trc_symmetry_nirrep(h, n)
      call check(error, rc, TRC_ERR_STATE, 'no irreps before the SCF')
      rc = trc_set_symmetry(h, 'bogus'//c_null_char, 0.0_c_double)
      call check(error, rc, TRC_ERR_BADARG, 'an unknown group is refused')
      if (allocated(error)) return
      call check(error, index(message(h), 'off, on/auto or a point group') > 0, 'its message: '//message(h))
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, -1.0_c_double)
      call check(error, rc, TRC_ERR_BADARG, 'a negative tolerance is refused')
      rc = trc_destroy(h)
   end subroutine test_state

   subroutine test_round_trip(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h, h1
      real(dp), allocatable :: geom(:, :), occ(:)
      real(c_double) :: e_sym, e_c1, resid, rot(9), shift(3), std(9), x(3), dev, rr(3, 3)
      real(c_double), allocatable :: g_sym(:), g_c1(:)
      integer(c_int) :: rc, nao, nirrep, dims(2)
      integer(c_int), allocatable :: irreps(:)
      integer :: ia, i, j

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'symmetric SCF: '//message(h))
      if (allocated(error)) return
      rc = trc_energy(h, e_sym)
      call check(error, name_of(h, 1, 0), 'Cs', 'group used')
      call check(error, name_of(h, 2, 0), 'Cs', 'group detected')
      if (allocated(error)) return
      rc = trc_symmetry_residual(h, resid)
      call check(error, rc == TRC_OK .and. resid < 1.0e-9_dp, 'a planar molecule needs no symmetrisation')
      if (allocated(error)) return

      ! The reported geometry and frame: x_std = rot (x - shift) for the input.
      rc = trc_geometry(h, std)
      call check(error, rc, TRC_OK, 'trc_geometry')
      rc = trc_symmetry_frame(h, rot, shift)
      call check(error, rc, TRC_OK, 'trc_symmetry_frame')
      if (allocated(error)) return
      rr = transpose(reshape(rot, [3, 3]))
      dev = maxval(abs(matmul(transpose(rr), rr) - reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])))
      call check(error, dev < 1.0e-12_dp, 'the frame rotation is orthogonal')
      dev = 0.0_dp
      do ia = 1, 3
         x = matmul(rr, geom(:, ia) - shift)
         dev = max(dev, maxval(abs(x - std(3*ia - 2:3*ia))))
      end do
      call check(error, dev <= resid + 1.0e-10_dp, 'the geometry is the input mapped by the frame ('//to_string(dev, 'es10.2')//')')
      if (allocated(error)) return

      ! The central oracle: a plain run on the geometry the context reports.
      call make_context(h1, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      ! the C1 oracle: symmetry switched off (a fresh context has symmetry off)
      rc = trc_set_symmetry(h1, 'off'//c_null_char, 0.0_c_double)
      rc = trc_set_molecule(h1, 3_c_int, [1.0_c_double, 8.0_c_double, 17.0_c_double], std, 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h1, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_run_scf(h1)
      call check(error, rc, TRC_OK, 'C1 run on the standard frame')
      rc = trc_energy(h1, e_c1)
      call check(error, abs(e_sym - e_c1) <= TOL_E, '|E(sym) - E(C1)| = '//to_string(abs(e_sym - e_c1), 'es10.2'))
      if (allocated(error)) return
      call check(error, name_of(h1, 1, 0), '<status 5>', 'the C1 context has no symmetry result')
      rc = trc_run_gradient(h)
      call check(error, rc, TRC_OK, 'symmetric gradient: '//message(h))
      rc = trc_run_gradient(h1)
      allocate (g_sym(9), g_c1(9))
      rc = trc_gradient(h, g_sym)
      rc = trc_gradient(h1, g_c1)
      call check(error, maxval(abs(g_sym - g_c1)) < 1.0e-9_dp, 'gradient vs C1: '// &
                 to_string(maxval(abs(g_sym - g_c1)), 'es10.2'))
      if (allocated(error)) return
      call check(error, maxval(abs(g_sym([3, 6, 9]))) < 1.0e-10_dp, 'no force out of the mirror plane')

      ! Irreps and labels.
      rc = trc_symmetry_nirrep(h, nirrep)
      call check(error, rc == TRC_OK .and. nirrep == 2, 'two irreps')
      call check(error, name_of(h, 3, 1), "A'", 'irrep 1')
      call check(error, name_of(h, 3, 2), 'A"', 'irrep 2')
      call check(error, name_of(h, 3, 3), '<status 2>', 'irrep 3 does not exist (TRC_ERR_BADARG)')
      if (allocated(error)) return
      rc = trc_symmetry_irrep_dim(h, dims)
      rc = trc_nao(h, nao)
      call check(error, sum(dims) == nao, 'SO dimensions sum to nao')
      allocate (irreps(nao), occ(nao))
      rc = trc_mo_irreps(h, irreps)
      call check(error, rc, TRC_OK, 'trc_mo_irreps')
      call check(error, all(irreps >= 1 .and. irreps <= 2), 'every orbital has an irrep')
      if (allocated(error)) return
      rc = trc_mo_occupations(h, occ)
      i = count(occ(1:nao) > 0.0_dp .and. irreps == 1)
      j = count(occ(1:nao) > 0.0_dp .and. irreps == 2)
      call check(error, i == 10 .and. j == 3, 'HOCl: 10 A'' and 3 A" occupied orbitals (pyscf), got '// &
                 to_string(i)//' and '//to_string(j))
      rc = trc_destroy(h)
      rc = trc_destroy(h1)
   end subroutine test_round_trip

   subroutine test_rerun(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: geom(:, :)
      real(c_double) :: e1, e2, g1(9), g2(9), resid
      integer(c_int) :: rc

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(h)
      rc = trc_energy(h, e1)
      rc = trc_geometry(h, g1)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'second run')
      rc = trc_energy(h, e2)
      rc = trc_geometry(h, g2)
      call check(error, abs(e1 - e2) < 1.0e-10_dp, 'the second run gives the same energy')
      call check(error, maxval(abs(g1 - g2)) == 0.0_dp, 'and does not move the molecule again')
      if (allocated(error)) return
      ! A new geometry (the input frame again) is detected afresh.
      rc = trc_set_geometry(h, reshape(geom, [9]))
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'after trc_set_geometry')
      rc = trc_energy(h, e2)
      call check(error, abs(e1 - e2) < 1.0e-9_dp, 'the same molecule in the input frame, the same energy')
      rc = trc_symmetry_residual(h, resid)
      call check(error, rc, TRC_OK, 'residual after the new geometry')
      rc = trc_destroy(h)
   end subroutine test_rerun

   subroutine test_skeleton_switch(error)
      !! The skeleton Fock build forced on against the full build (then
      !! symmetrized): the energy and the gradient are the same, and the
      !! gradient that follows a skeleton SCF does not reuse its orbit-closed
      !! engine. Forced on explicitly because HOCl is a molecule the size rule
      !! declines (3 atoms in 3 orbits), and without the call this would
      !! compare the full build with itself.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: geom(:, :)
      real(c_double) :: e1, e2, g1(9), g2(9)
      integer(c_int) :: rc

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      rc = trc_set_symmetry_skeleton(h, 1_c_int)
      call check(error, rc, TRC_OK, 'trc_set_symmetry_skeleton(1)')
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'skeleton run: '//message(h))
      rc = trc_energy(h, e1)
      rc = trc_run_gradient(h)
      call check(error, rc, TRC_OK, 'gradient after the skeleton run')
      rc = trc_gradient(h, g1)
      if (allocated(error)) return
      rc = trc_set_symmetry_skeleton(h, 0_c_int)
      call check(error, rc, TRC_OK, 'trc_set_symmetry_skeleton(0)')
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'full-build run: '//message(h))
      rc = trc_energy(h, e2)
      call check(error, abs(e1 - e2) < 1.0e-9_dp, 'the same energy with the skeleton switched off')
      rc = trc_run_gradient(h)
      rc = trc_gradient(h, g2)
      call check(error, maxval(abs(g1 - g2)) < 1.0e-7_dp, 'and the same gradient')
      ! -1 hands the choice back to the size rule; anything else is refused
      rc = trc_set_symmetry_skeleton(h, -1_c_int)
      call check(error, rc, TRC_OK, 'trc_set_symmetry_skeleton(-1)')
      if (allocated(error)) return
      rc = trc_set_symmetry_skeleton(h, 2_c_int)
      call check(error, rc, TRC_ERR_BADARG, 'trc_set_symmetry_skeleton(2) is refused')
      if (allocated(error)) return
      rc = trc_destroy(h)
      rc = trc_set_symmetry_skeleton(c_null_ptr, 1_c_int)
      call check(error, rc /= TRC_OK, 'a null handle is refused')
   end subroutine test_skeleton_switch

   subroutine test_irrep_nelec(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: geom(:, :), occ(:)
      real(c_double) :: e_auto, e_fixed, e_exc
      integer(c_int), allocatable :: irreps(:)
      integer(c_int) :: rc, nao, fr

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      ! irrep_nelec with symmetry off is a mistake, not a no-op: on a fresh
      ! context (symmetry is opt-in) and on an explicit off alike.
      rc = trc_set_irrep_nelec(h, "A'=20,A"""//'=6'//c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_BADARG, 'irrep_nelec on a fresh (symmetry off) context')
      call check(error, index(message(h), 'needs a symmetry request') > 0, message(h))
      if (allocated(error)) return
      rc = trc_set_symmetry(h, 'off'//c_null_char, 0.0_c_double)
      rc = trc_set_irrep_nelec(h, "A'=20,A"""//'=6'//c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_BADARG, 'irrep_nelec needs symmetry')
      call check(error, index(message(h), 'needs a symmetry request') > 0, message(h))
      if (allocated(error)) return
      rc = trc_destroy(h)
      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'on'//c_null_char, 0.0_c_double)
      rc = trc_set_irrep_nelec(h, "A'=20,A"""//'=6'//c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'irrep_nelec with on is a strict request: '//message(h))
      call check(error, name_of(h, 1, 0), 'Cs', 'and it runs in Cs')
      rc = trc_symmetry_report_frame(h, fr)
      call check(error, fr, 0_c_int, 'with the standard frame reported (fixed occupations make the request strict)')
      if (allocated(error)) return
      rc = trc_set_irrep_nelec(h, c_null_char)
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(h)
      rc = trc_energy(h, e_auto)
      rc = trc_set_irrep_nelec(h, "A'=20,A"""//'=6'//c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'aufbau counts: '//message(h))
      rc = trc_energy(h, e_fixed)
      call check(error, abs(e_auto - e_fixed) < 1.0e-9_dp, 'the aufbau counts give the aufbau energy')
      rc = trc_run_gradient(h)
      call check(error, rc, TRC_OK, 'a gradient for the aufbau counts')
      if (allocated(error)) return

      rc = trc_set_irrep_nelec(h, "A'=18,A"""//'=8'//c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'excited counts: '//message(h))
      rc = trc_energy(h, e_exc)
      call check(error, e_exc > e_auto + 0.1_dp, 'an excited occupation is a higher state')
      if (allocated(error)) return
      rc = trc_nao(h, nao)
      allocate (irreps(nao), occ(nao))
      rc = trc_mo_irreps(h, irreps)
      rc = trc_mo_occupations(h, occ)
      call check(error, count(occ > 0.0_dp .and. irreps == 1) == 9 .and. count(occ > 0.0_dp .and. irreps == 2) == 4, &
                 'trc_mo_occupations follows the fixed counts (9 A'', 4 A")')
      if (allocated(error)) return
      rc = trc_run_gradient(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, 'a gradient for non-aufbau counts is refused')
      call check(error, index(message(h), 'not the aufbau occupation') > 0, message(h))
      if (allocated(error)) return

      rc = trc_set_irrep_nelec(h, "A'=20,A"""//'=8'//c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_BADARG, 'counts that do not add up')
      call check(error, index(message(h), 'does not add up') > 0, message(h))
      rc = trc_set_irrep_nelec(h, 'B1=26'//c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_BADARG, 'a label the group does not have')
      call check(error, index(message(h), 'is not an irrep of Cs') > 0, message(h))
      rc = trc_set_irrep_nelec(h, c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'an empty spec restores aufbau')
      rc = trc_destroy(h)
   end subroutine test_irrep_nelec

   subroutine test_explicit(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: geom(:, :)
      integer(c_int) :: rc

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'c2v'//c_null_char, 0.0_c_double)
      call check(error, rc, TRC_OK, 'the setter takes any group name')
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_BADARG, 'C2v on a Cs molecule is refused')
      call check(error, index(message(h), 'is not that group or one of its subgroups') > 0, message(h))
      if (allocated(error)) return
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'Cs on a Cs molecule')
      call check(error, name_of(h, 1, 0), 'Cs', 'Cs used')
      if (allocated(error)) return
      rc = trc_set_symmetry(h, 'c1'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'C1')
      call check(error, name_of(h, 1, 0), '<status 5>', 'C1 requested: the old path, no symmetry result')
      rc = trc_destroy(h)
      if (allocated(error)) return

      call make_context(h, 'water.xyz', 'basis_sets/6-31g.json', geom)
      rc = trc_set_symmetry(h, 'c2v'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'C2v on water')
      call check(error, name_of(h, 1, 0), 'C2v', 'C2v used')
      if (allocated(error)) return
      rc = trc_set_symmetry(h, 'cs'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'Cs on water (a subgroup of the detected C2v)')
      call check(error, name_of(h, 1, 0), 'Cs', 'Cs used')
      call check(error, name_of(h, 2, 0), 'C2v', 'C2v detected')
      if (allocated(error)) return
      rc = trc_set_symmetry(h, 'ci'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_BADARG, 'Ci on water is not a subgroup')
      call check(error, index(message(h), 'is not that group or one of its subgroups') > 0, message(h))
      rc = trc_destroy(h)
   end subroutine test_explicit

   subroutine test_scf_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: geom(:, :)
      real(c_double), allocatable, target :: dg(:)
      integer(c_int) :: rc, nao

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      rc = trc_set_df(h, 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, 'density fitting')
      call check(error, index(message(h), 'density fitting') > 0, message(h))
      rc = trc_set_df(h, 0_c_int)
      rc = trc_set_pcm(h, 'cpcm'//c_null_char, 78.4_c_double, 0_c_int, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, 'PCM')
      call check(error, index(message(h), 'PCM') > 0, message(h))
      if (allocated(error)) return
      rc = trc_destroy(h)

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      rc = trc_set_soscf(h, 1_c_int, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc /= TRC_OK, 'SOSCF is refused with symmetry')
      call check(error, index(message(h), 'soscf') > 0, message(h))
      rc = trc_destroy(h)
      if (allocated(error)) return

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      rc = trc_nao(h, nao)
      allocate (dg(nao*nao))
      dg = 0.0_c_double
      rc = trc_set_guess(h, TRC_GUESS_GIVEN, c_loc(dg), 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc /= TRC_OK, 'a given guess (input frame) is refused with symmetry')
      rc = trc_destroy(h)
   end subroutine test_scf_refusals

   subroutine test_entry_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: geom(:, :)
      real(c_double) :: f(3, 6)
      integer(c_int) :: rc

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'Cs'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, 'the SCF itself runs')
      if (allocated(error)) return
      rc = trc_run_rimp2(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, 'RI-MP2')
      call check(error, index(message(h), 'symmetry is not supported with RI-MP2') > 0, message(h))
      rc = trc_check_rimp2_gradient(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, 'the RI-MP2 gradient check')
      rc = trc_run_hessian(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, 'Hessian')
      call check(error, index(message(h), 'symmetry is not supported with the finite-difference Hessian') > 0, message(h))
      rc = trc_run_polarizability(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, 'CPHF/CPKS')
      call check(error, index(message(h), 'symmetry is not supported with CPHF/CPKS') > 0, message(h))
      rc = trc_load_scf(h, 'none.scf'//c_null_char)
      call check(error, rc, TRC_ERR_UNSUPPORTED, 'restart')
      rc = trc_fukui_hirshfeld(h, f(:, 1), f(:, 2), f(:, 3), f(:, 4), f(:, 5), f(:, 6))
      call check(error, rc, TRC_ERR_UNSUPPORTED, 'Fukui')
      rc = trc_destroy(h)
   end subroutine test_entry_refusals


   ! ---- lenient requests ("on"/"auto"): input frame, C1 fallback -------------

   subroutine make_pair(hd, ho, xyz_file, basis, cart)
      !! A lenient ("on") context and a symmetry-off one on the same molecule.
      type(c_ptr), intent(out) :: hd, ho
      character(len=*), intent(in) :: xyz_file, basis
      logical, intent(in), optional :: cart
      real(dp), allocatable :: geom(:, :)
      integer(c_int) :: rc
      call make_context(hd, xyz_file, basis, geom)
      call make_context(ho, xyz_file, basis, geom)
      rc = trc_set_symmetry(hd, 'on'//c_null_char, 0.0_c_double)
      rc = trc_set_symmetry(ho, 'off'//c_null_char, 0.0_c_double)
      if (present(cart)) then
         if (cart) then
            rc = trc_set_force_cartesian(hd, 1_c_int)
            rc = trc_set_force_cartesian(ho, 1_c_int)
            rc = trc_set_basis_json(hd, basis//c_null_char)
            rc = trc_set_basis_json(ho, basis//c_null_char)
         end if
      end if
   end subroutine make_pair

   subroutine test_default_round_trip(error)
      !! Water in a general orientation (C2v), cc-pVDZ (spherical, so the d
      !! shells use the spherical AO rotation), an "on" context against one with
      !! symmetry off.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hd, ho
      real(dp), allocatable :: geom(:, :)
      real(c_double), allocatable :: dd(:), do_(:), cd(:), eps_d(:), eps_o(:), p(:)
      real(c_double) :: ed, eo, gd(9), go(9), resid, xyz(9), q_d(3), q_o(3), std(9)
      integer(c_int) :: rc, nao, nmo, fr
      integer :: n, i, nocc, ia
      real(dp) :: dev

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/cc-pvdz.json')
      call read_geom('water_rot.xyz', geom)
      rc = trc_run_scf(hd)
      call check(error, rc, TRC_OK, 'on SCF: '//message(hd))
      rc = trc_run_scf(ho)
      call check(error, rc, TRC_OK, 'symmetry-off SCF')
      if (allocated(error)) return
      call check(error, name_of(hd, 1, 0), 'C2v', 'on used C2v on a rotated water')
      call check(error, name_of(ho, 1, 0), '<status 5>', 'the off context has no symmetry result')
      if (allocated(error)) return
      rc = trc_energy(hd, ed)
      rc = trc_energy(ho, eo)
      print '(a,es10.2)', '  on vs off, |dE| ', abs(ed - eo)
      call check(error, abs(ed - eo) <= 1.0e-10_dp, '|E(on) - E(off)| '//to_string(abs(ed - eo), 'es10.2'))
      if (allocated(error)) return

      ! orbital energies, density, rotated-back orbitals
      rc = trc_nao(hd, nao)
      rc = trc_nmo(hd, nmo)
      n = nao
      allocate (dd(n*n), do_(n*n), cd(n*n), eps_d(n), eps_o(n), p(n*n))
      rc = trc_mo_energies(hd, eps_d)
      rc = trc_mo_energies(ho, eps_o)
      print '(a,es10.2)', '  orbital energies, max diff ', maxval(abs(eps_d - eps_o))
      call check(error, maxval(abs(eps_d - eps_o)) < 1.0e-8_dp, 'orbital energies differ: '// &
                 to_string(maxval(abs(eps_d - eps_o)), 'es10.2'))
      rc = trc_density(hd, dd)
      rc = trc_density(ho, do_)
      print '(a,es10.2)', '  density (input frame), max diff ', maxval(abs(dd - do_))
      call check(error, maxval(abs(dd - do_)) < 1.0e-9_dp, 'density differs from off: '// &
                 to_string(maxval(abs(dd - do_)), 'es10.2'))
      if (allocated(error)) return
      rc = trc_mo_coefficients(hd, cd)
      nocc = 5
      p = 0.0_dp
      do i = 1, nocc
         p = p + 2.0_dp*reshape(spread(cd((i - 1)*n + 1:i*n), 2, n)*spread(cd((i - 1)*n + 1:i*n), 1, n), [n*n])
      end do
      print '(a,es10.2)', '  2 C_occ C_occ^T (rotated-back orbitals) vs the off density, max diff ', maxval(abs(p - do_))
      call check(error, maxval(abs(p - do_)) < 1.0e-9_dp, 'rotated-back orbitals do not reproduce the off density: '// &
                 to_string(maxval(abs(p - do_)), 'es10.2'))
      if (allocated(error)) return

      ! the geometry is the input, symmetrised in the input frame
      rc = trc_geometry(hd, xyz)
      rc = trc_symmetry_residual(hd, resid)
      dev = 0.0_dp
      do ia = 1, 3
         dev = max(dev, maxval(abs(xyz(3*ia - 2:3*ia) - geom(:, ia))))
      end do
      print '(a,es10.2,a,es10.2)', '  reported geometry vs input ', dev, '   residual ', resid
      call check(error, dev <= resid + 1.0e-12_dp, 'the reported geometry is not the input within the residual')
      rc = trc_symmetry_report_frame(hd, fr)
      call check(error, fr, 1_c_int, 'an on request reports the input frame')
      if (allocated(error)) return

      ! gradient
      rc = trc_run_gradient(hd)
      call check(error, rc, TRC_OK, 'on gradient: '//message(hd))
      rc = trc_run_gradient(ho)
      rc = trc_gradient(hd, gd)
      rc = trc_gradient(ho, go)
      print '(a,es10.2)', '  gradient (input frame), max diff ', maxval(abs(gd - go))
      call check(error, maxval(abs(gd - go)) < 1.0e-9_dp, 'gradient (input frame) differs from off: '// &
                 to_string(maxval(abs(gd - go)), 'es10.2'))
      if (allocated(error)) return

      ! atomic charges: atom-based, evaluated in the input frame, so equal
      rc = trc_hirshfeld_charges(hd, q_d)
      rc = trc_hirshfeld_charges(ho, q_o)
      print '(a,es10.2)', '  Hirshfeld charges, max diff ', maxval(abs(q_d - q_o))
      call check(error, maxval(abs(q_d - q_o)) < 1.0e-9_dp, 'Hirshfeld charges differ: '// &
                 to_string(maxval(abs(q_d - q_o)), 'es10.2'))
      rc = trc_chelpg_charges(hd, q_d)
      rc = trc_chelpg_charges(ho, q_o)
      print '(a,es10.2)', '  CHELPG charges, max diff ', maxval(abs(q_d - q_o))
      call check(error, maxval(abs(q_d - q_o)) < 1.0e-9_dp, 'CHELPG charges differ: '// &
                 to_string(maxval(abs(q_d - q_o)), 'es10.2'))
      if (allocated(error)) return

      ! the report frame can be switched: standard frame = the group's frame
      rc = trc_set_symmetry_report_frame(hd, 0_c_int)
      call check(error, rc, TRC_OK, 'report frame 0')
      rc = trc_symmetry_report_frame(hd, fr)
      call check(error, fr, 0_c_int, 'reads back as the standard frame')
      rc = trc_geometry(hd, std)
      call check(error, maxval(abs(std - xyz)) > 0.1_dp, 'the standard-frame geometry of a rotated molecule differs from the input')
      rc = trc_set_symmetry_report_frame(hd, 7_c_int)
      call check(error, rc, TRC_ERR_BADARG, 'a bad frame is refused')
      rc = trc_set_symmetry_report_frame(hd, -1_c_int)
      rc = trc_symmetry_report_frame(hd, fr)
      call check(error, fr, 1_c_int, 'back to the lenient rule: the input frame')
      rc = trc_destroy(hd)
      rc = trc_destroy(ho)
   end subroutine test_default_round_trip

   subroutine test_fresh_is_off(error)
      !! Symmetry is opt-in. A fresh context, and one asked for "default", run C1 (no
      !! symmetry result) with the numbers of an explicit "off"; "on" and "auto" are the
      !! same lenient request (a symmetry result, the input frame).
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hf, hd, ho, hn, ha
      real(dp), allocatable :: geom(:, :)
      real(c_double) :: ef, ed, eo, en, ea, gf(9), go(9), xf(9), xo(9)
      integer(c_int) :: rc, fr

      call make_context(hf, 'water_rot.xyz', 'basis_sets/cc-pvdz.json', geom)
      call make_context(hd, 'water_rot.xyz', 'basis_sets/cc-pvdz.json', geom)
      call make_context(ho, 'water_rot.xyz', 'basis_sets/cc-pvdz.json', geom)
      call make_context(hn, 'water_rot.xyz', 'basis_sets/cc-pvdz.json', geom)
      call make_context(ha, 'water_rot.xyz', 'basis_sets/cc-pvdz.json', geom)
      rc = trc_set_symmetry(hd, 'default'//c_null_char, 0.0_c_double)
      call check(error, rc, TRC_OK, '"default" is accepted')
      rc = trc_set_symmetry(ho, 'off'//c_null_char, 0.0_c_double)
      rc = trc_set_symmetry(hn, 'on'//c_null_char, 0.0_c_double)
      call check(error, rc, TRC_OK, '"on" is accepted')
      rc = trc_set_symmetry(ha, 'auto'//c_null_char, 0.0_c_double)
      call check(error, rc, TRC_OK, '"auto" is accepted')
      rc = trc_run_scf(hf)
      call check(error, rc, TRC_OK, 'fresh SCF: '//message(hf))
      rc = trc_run_scf(hd)
      rc = trc_run_scf(ho)
      rc = trc_run_scf(hn)
      rc = trc_run_scf(ha)
      if (allocated(error)) return
      call check(error, name_of(hf, 1, 0), '<status 5>', 'a fresh context has no symmetry result (C1)')
      call check(error, name_of(hd, 1, 0), '<status 5>', '"default" has no symmetry result (C1)')
      call check(error, name_of(ho, 1, 0), '<status 5>', '"off" has no symmetry result (C1)')
      call check(error, name_of(hn, 1, 0), 'C2v', '"on" used C2v on a rotated water')
      call check(error, name_of(ha, 1, 0), 'C2v', '"auto" is the same request as "on"')
      if (allocated(error)) return
      rc = trc_energy(hf, ef)
      rc = trc_energy(hd, ed)
      rc = trc_energy(ho, eo)
      rc = trc_energy(hn, en)
      rc = trc_energy(ha, ea)
      call check(error, abs(ef - eo) <= 1.0e-12_dp .and. abs(ed - eo) <= 1.0e-12_dp, 'fresh and "default" energies equal "off" exactly')
      call check(error, abs(en - eo) <= 1.0e-10_dp .and. abs(ea - en) <= 1.0e-12_dp, '"on" and "auto" agree with "off"')
      rc = trc_run_gradient(hf)
      rc = trc_run_gradient(ho)
      rc = trc_gradient(hf, gf)
      rc = trc_gradient(ho, go)
      rc = trc_geometry(hf, xf)
      rc = trc_geometry(ho, xo)
      call check(error, maxval(abs(gf - go)) <= 1.0e-12_dp, 'fresh gradient equals "off"')
      call check(error, maxval(abs(xf - xo)) == 0.0_dp, 'fresh geometry is the input, as with "off"')
      rc = trc_symmetry_report_frame(hn, fr)
      call check(error, fr, 1_c_int, '"on" reports the input frame')
      rc = trc_destroy(hf); rc = trc_destroy(hd); rc = trc_destroy(ho); rc = trc_destroy(hn); rc = trc_destroy(ha)
   end subroutine test_fresh_is_off

   subroutine read_geom(xyz_file, geom)
      character(len=*), intent(in) :: xyz_file
      real(dp), allocatable, intent(out) :: geom(:, :)
      integer, allocatable :: elem(:)
      integer :: nat
      call read_xyz(xyz_file, nat, elem, geom)
   end subroutine read_geom

   function dipole(basis, dens, geom, z, origin, error) result(d)
      !! Nuclear minus electronic dipole about `origin` from a (Cartesian) basis
      !! handle, its density and the first three multipole components (x, y, z).
      type(c_ptr), intent(in) :: basis
      real(c_double), intent(in) :: dens(:)
      real(dp), intent(in) :: geom(:, :), z(:)
      real(c_double), intent(in) :: origin(3)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: d(3)
      type(c_ptr) :: pairs
      real(c_double), allocatable :: mm(:)
      integer :: n, k, ia
      integer(c_int) :: rc
      n = nint(sqrt(real(size(dens), dp)))
      allocate (mm(n*n*39))
      rc = trc_pairs_create(basis, 1.0e-14_c_double, pairs)
      call check(error, rc, TRC_OK, 'trc_pairs_create')
      if (allocated(error)) return
      rc = trc_compute_multipoles(basis, pairs, origin, mm)
      call check(error, rc, TRC_OK, 'trc_compute_multipoles')
      rc = trc_pairs_destroy(pairs)
      if (allocated(error)) return
      do k = 1, 3
         d(k) = -sum(dens*mm((k - 1)*n*n + 1:k*n*n))
      end do
      do ia = 1, size(z)
         d = d + z(ia)*(geom(:, ia) - origin)
      end do
   end function dipole

   subroutine test_default_borrowed(error)
      !! The basis an "on" context lends (trc_context_basis) is the input-frame
      !! one: its overlap and the dipole it gives with trc_density equal those of
      !! a symmetry-off context. Cartesian AOs, as the multipole integrals are.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hd, ho, bd, bo, pd, po
      real(dp), allocatable :: geom(:, :)
      real(c_double), allocatable :: dd(:), do_(:), sd(:), so(:), td(:), vd(:)
      real(c_double) :: o1(3), o2(3), xyz(9)
      real(dp) :: dip_d(3), dip_o(3), dip_o2(3), z(3), g(3, 3)
      integer(c_int) :: rc, nao
      integer :: n

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/cc-pvdz.json', cart=.true.)
      call read_geom('water_rot.xyz', geom)
      rc = trc_run_scf(hd)
      rc = trc_run_scf(ho)
      call check(error, rc, TRC_OK, 'SCF')
      call check(error, name_of(hd, 1, 0), 'C2v', 'on used C2v')
      if (allocated(error)) return
      rc = trc_nao(hd, nao)
      n = nao
      allocate (dd(n*n), do_(n*n), sd(n*n), so(n*n), td(n*n), vd(n*n))
      rc = trc_density(hd, dd)
      rc = trc_density(ho, do_)
      print '(a,es10.2)', '  Cartesian density, max diff ', maxval(abs(dd - do_))
      call check(error, maxval(abs(dd - do_)) < 1.0e-9_dp, 'density (Cartesian) differs from off')
      if (allocated(error)) return
      rc = trc_context_basis(hd, bd)
      rc = trc_context_basis(ho, bo)
      rc = trc_pairs_create(bd, 1.0e-14_c_double, pd)
      rc = trc_pairs_create(bo, 1.0e-14_c_double, po)
      rc = trc_compute_1e(bd, pd, sd, td, vd)
      rc = trc_compute_1e(bo, po, so, td, vd)
      rc = trc_pairs_destroy(pd)
      rc = trc_pairs_destroy(po)
      print '(a,es10.2)', '  overlap from the borrowed basis, max diff ', maxval(abs(sd - so))
      call check(error, maxval(abs(sd - so)) < 1.0e-9_dp, 'the borrowed basis is not the input-frame basis: overlap differs '// &
                 to_string(maxval(abs(sd - so)), 'es10.2'))
      if (allocated(error)) return

      z = [8.0_dp, 1.0_dp, 1.0_dp]
      rc = trc_geometry(hd, xyz)
      g = reshape(real(xyz, dp), [3, 3])
      o1 = 0.0_c_double
      o2 = [0.7_c_double, -1.1_c_double, 2.3_c_double]
      dip_d = dipole(bd, dd, g, z, o1, error)
      if (allocated(error)) return
      dip_o = dipole(bo, do_, geom, z, o1, error)
      if (allocated(error)) return
      dip_o2 = dipole(bo, do_, geom, z, o2, error)
      if (allocated(error)) return
      print '(a,3f12.8,a,es10.2)', '  dipole (a.u.), on     ', dip_d, '   |on - off| ', maxval(abs(dip_d - dip_o))
      call check(error, maxval(abs(dip_o - dip_o2)) < 1.0e-6_dp, 'the neutral molecule dipole must not depend on the origin '// &
                 '(wrong sign convention in the test)')
      call check(error, maxval(abs(dip_d - dip_o)) < 1.0e-8_dp, 'the dipole vector is not in the input frame: '// &
                 to_string(maxval(abs(dip_d - dip_o)), 'es10.2'))
      if (allocated(error)) return
      rc = trc_destroy(hd)
      rc = trc_destroy(ho)
   end subroutine test_default_borrowed

   subroutine test_default_fallbacks(error)
      !! Each setting or entry the symmetric SCF cannot serve, on an "on"
      !! context, gives what a symmetry-off context gives, and (for the
      !! settings) leaves no symmetry result.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hd, ho
      real(c_double) :: ed, eo, hd9(81), ho9(81), ad(9), ao(9), ed2, eo2, f(3, 6), g(3, 6)
      real(c_double), allocatable, target :: dg(:)
      integer(c_int) :: rc, nao

      ! --- settings of the SCF ---
      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json', cart=.true.)
      rc = trc_set_jkfit_json(hd, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
      rc = trc_set_jkfit_json(ho, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
      rc = trc_set_df(hd, 1_c_int)
      rc = trc_set_df(ho, 1_c_int)
      call compare_runs('density fitting', hd, ho, error)
      if (allocated(error)) return
      rc = trc_destroy(hd); rc = trc_destroy(ho)

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json')
      rc = trc_set_pcm(hd, 'cpcm'//c_null_char, 78.4_c_double, 0_c_int, 0.0_c_double)
      rc = trc_set_pcm(ho, 'cpcm'//c_null_char, 78.4_c_double, 0_c_int, 0.0_c_double)
      call compare_runs('PCM', hd, ho, error)
      if (allocated(error)) return
      rc = trc_destroy(hd); rc = trc_destroy(ho)

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json')
      rc = trc_set_soscf(hd, 1_c_int, 0.0_c_double)
      rc = trc_set_soscf(ho, 1_c_int, 0.0_c_double)
      call compare_runs('SOSCF', hd, ho, error)
      if (allocated(error)) return
      rc = trc_destroy(hd); rc = trc_destroy(ho)

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json')
      rc = trc_nao(hd, nao)
      allocate (dg(nao*nao))
      dg = 0.0_c_double
      rc = trc_set_guess(hd, TRC_GUESS_GIVEN, c_loc(dg), 1_c_int)
      rc = trc_set_guess(ho, TRC_GUESS_GIVEN, c_loc(dg), 1_c_int)
      call compare_runs('a given guess', hd, ho, error)
      if (allocated(error)) return
      rc = trc_destroy(hd); rc = trc_destroy(ho)

      ! --- entries that follow a (symmetric) "on" SCF ---
      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json')
      rc = trc_run_scf(hd)
      rc = trc_run_scf(ho)
      call check(error, name_of(hd, 1, 0), 'C2v', 'the on SCF is symmetric before the Hessian')
      if (allocated(error)) return
      rc = trc_run_hessian(hd)
      call check(error, rc, TRC_OK, 'Hessian on an on context: '//message(hd))
      rc = trc_run_hessian(ho)
      if (allocated(error)) return
      call check(error, name_of(hd, 1, 0), '<status 5>', 'after the Hessian the context is C1')
      rc = trc_hessian(hd, hd9)
      rc = trc_hessian(ho, ho9)
      print '(a,es10.2)', '  Hessian (on context) vs off, max diff ', maxval(abs(hd9 - ho9))
      call check(error, maxval(abs(hd9 - ho9)) < 1.0e-9_dp, 'the Hessian differs from the symmetry-off Hessian: '// &
                 to_string(maxval(abs(hd9 - ho9)), 'es10.2'))
      if (allocated(error)) return
      rc = trc_destroy(hd); rc = trc_destroy(ho)

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json')
      rc = trc_run_scf(hd)
      rc = trc_run_scf(ho)
      rc = trc_run_polarizability(hd)
      call check(error, rc, TRC_OK, 'polarizability on an on context: '//message(hd))
      rc = trc_run_polarizability(ho)
      rc = trc_polarizability(hd, ad)
      rc = trc_polarizability(ho, ao)
      call check(error, maxval(abs(ad - ao)) < 1.0e-9_dp, 'the polarizability differs from the symmetry-off one: '// &
                 to_string(maxval(abs(ad - ao)), 'es10.2'))
      if (allocated(error)) return
      rc = trc_destroy(hd); rc = trc_destroy(ho)

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json')
      rc = trc_run_scf(hd)
      rc = trc_run_scf(ho)
      rc = trc_fukui_hirshfeld(hd, f(:, 1), f(:, 2), f(:, 3), g(:, 1), g(:, 2), g(:, 3))
      call check(error, rc, TRC_OK, 'Fukui on an on context: '//message(hd))
      rc = trc_fukui_hirshfeld(ho, f(:, 4), f(:, 5), f(:, 6), g(:, 4), g(:, 5), g(:, 6))
      call check(error, maxval(abs(f(:, 1:3) - f(:, 4:6))) < 1.0e-9_dp .and. maxval(abs(g(:, 1:3) - g(:, 4:6))) < 1.0e-9_dp, &
                 'the Fukui indices differ from the symmetry-off ones')
      if (allocated(error)) return
      rc = trc_destroy(hd); rc = trc_destroy(ho)

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/cc-pvdz.json')
      rc = trc_set_aux_json(hd, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      rc = trc_set_aux_json(ho, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      rc = trc_run_scf(hd)
      rc = trc_run_scf(ho)
      rc = trc_run_rimp2(hd)
      call check(error, rc, TRC_OK, 'RI-MP2 on an on context: '//message(hd))
      rc = trc_run_rimp2(ho)
      rc = trc_rimp2_energy(hd, ed, ed2)
      rc = trc_rimp2_energy(ho, eo, eo2)
      call check(error, abs(ed - eo) < 1.0e-10_dp .and. abs(ed2 - eo2) < 1.0e-10_dp, 'the RI-MP2 energy differs from the symmetry-off one')
      rc = trc_destroy(hd); rc = trc_destroy(ho)
   end subroutine test_default_fallbacks

   ! ---- review regressions: per-atom shells, frames of stored orbitals, ----
   ! ---- irrep_nelec order, the grid after a C1 rerun                     ----

   subroutine water_arrays(geom, h2_scale, nsh, l, np, ex, co, cen)
      !! STO-3G water on `geom` (O, H, H; Bohr) as trc_set_basis_arrays takes
      !! it; the second H's exponents are scaled by `h2_scale`.
      real(dp), intent(in) :: geom(:, :), h2_scale
      integer(c_int), intent(out) :: nsh
      integer(c_int), allocatable, intent(out) :: l(:), np(:)
      real(c_double), allocatable, intent(out) :: ex(:, :), co(:, :), cen(:, :)
      real(dp), parameter :: E1S(3) = [130.70932_dp, 23.808861_dp, 6.4436083_dp]
      real(dp), parameter :: ESP(3) = [5.0331513_dp, 1.1695961_dp, 0.380389_dp]
      real(dp), parameter :: EH(3) = [3.42525091_dp, 0.62391373_dp, 0.1688554_dp]
      real(dp), parameter :: C1S(3) = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
      real(dp), parameter :: C2S(3) = [-0.09996723_dp, 0.39951283_dp, 0.70011547_dp]
      real(dp), parameter :: C2P(3) = [0.15591627_dp, 0.60768372_dp, 0.39195739_dp]
      nsh = 5
      l = [0, 0, 1, 0, 0]
      np = [3, 3, 3, 3, 3]
      ex = reshape([E1S, ESP, ESP, EH, EH*h2_scale], [3, 5])
      co = reshape([C1S, C2S, C2P, C1S, C1S], [3, 5])
      cen = reshape([geom(:, 1), geom(:, 1), geom(:, 1), geom(:, 2), geom(:, 3)], [3, 5])
   end subroutine water_arrays

   subroutine test_per_atom_shells(error)
      !! Water whose two H carry different shells (trc_set_basis_arrays allows
      !! that): only the molecular plane is a symmetry of the basis. "on" runs
      !! Cs, with the symmetry-off energy; C2v is refused as not a subgroup.
      !! Matching atoms on the element alone called the H equivalent, and the
      !! SO construction then refused the run.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hd, ho
      real(dp), allocatable :: geom(:, :)
      integer(c_int), allocatable :: l(:), np(:)
      real(c_double), allocatable :: ex(:, :), co(:, :), cen(:, :)
      real(c_double) :: ed, eo
      integer(c_int) :: rc, nsh

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json')
      call read_geom('water_rot.xyz', geom)
      call water_arrays(geom, 1.1_dp, nsh, l, np, ex, co, cen)
      rc = trc_set_basis_arrays(hd, nsh, 3_c_int, l, np, ex, co, cen)
      call check(error, rc, TRC_OK, 'basis arrays: '//message(hd))
      rc = trc_set_basis_arrays(ho, nsh, 3_c_int, l, np, ex, co, cen)
      if (allocated(error)) return
      rc = trc_run_scf(hd)
      call check(error, rc, TRC_OK, 'on with different shells on the two H: '//message(hd))
      if (allocated(error)) return
      rc = trc_run_scf(ho)
      call check(error, name_of(hd, 1, 0), 'Cs', 'the H with other shells leaves only the plane')
      if (allocated(error)) return
      rc = trc_energy(hd, ed)
      rc = trc_energy(ho, eo)
      call check(error, abs(ed - eo) <= TOL_E, 'E(on) - E(off) '//to_string(abs(ed - eo), 'es10.2'))
      if (allocated(error)) return
      rc = trc_set_symmetry(hd, 'C2v'//c_null_char, 0.0_c_double)
      rc = trc_run_scf(hd)
      call check(error, rc /= TRC_OK, 'C2v is not a symmetry of this basis')
      if (allocated(error)) return
      ! equal shells on both H: C2v again
      call water_arrays(geom, 1.0_dp, nsh, l, np, ex, co, cen)
      rc = trc_set_symmetry(hd, 'on'//c_null_char, 0.0_c_double)
      rc = trc_set_basis_arrays(hd, nsh, 3_c_int, l, np, ex, co, cen)
      rc = trc_run_scf(hd)
      call check(error, rc, TRC_OK, 'on with equal shells: '//message(hd))
      if (allocated(error)) return
      call check(error, name_of(hd, 1, 0), 'C2v', 'equal shells keep C2v')
      rc = trc_destroy(hd); rc = trc_destroy(ho)
   end subroutine test_per_atom_shells

   subroutine test_projected_after_on(error)
      !! A projected guess after a symmetric "on" run (it falls back to C1)
      !! projects the stored orbitals onto the input geometry, so they must be
      !! the input frame's: the run converges as fast as the same sequence
      !! with symmetry off, to the same energy. Orbitals left in the standard
      !! frame of a rotated water are a poor start.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hd, ho
      real(c_double) :: ed, eo
      integer(c_int) :: rc, itd, ito

      call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json')
      rc = trc_run_scf(hd)
      rc = trc_run_scf(ho)
      call check(error, name_of(hd, 1, 0), 'C2v', 'the first on run is symmetric')
      if (allocated(error)) return
      rc = trc_set_guess(hd, TRC_GUESS_PROJECTED, c_null_ptr, 1_c_int)
      rc = trc_set_guess(ho, TRC_GUESS_PROJECTED, c_null_ptr, 1_c_int)
      rc = trc_run_scf(hd)
      call check(error, rc, TRC_OK, 'projected guess after on: '//message(hd))
      rc = trc_run_scf(ho)
      if (allocated(error)) return
      rc = trc_energy(hd, ed)
      rc = trc_energy(ho, eo)
      rc = trc_iterations(hd, itd)
      rc = trc_iterations(ho, ito)
      print '(a,i0,a,i0)', '  projected guess: iterations after on ', itd, ', after off ', ito
      call check(error, abs(ed - eo) <= TOL_E, 'E differs: '//to_string(abs(ed - eo), 'es10.2'))
      if (allocated(error)) return
      call check(error, itd <= ito + 1, 'the projected start after on is no worse than after off ('// &
                 to_string(itd)//' against '//to_string(ito)//' iterations)')
      rc = trc_destroy(hd); rc = trc_destroy(ho)
   end subroutine test_projected_after_on

   subroutine test_irrep_nelec_order(error)
      !! irrep_nelec makes "on" strict whether it is set before or after the
      !! request, and clearing it makes "on" lenient again: density fitting is
      !! refused, then runs in C1 with the symmetry-off energy.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hd, ho
      real(c_double) :: ed, eo
      integer(c_int) :: rc
      integer :: k

      do k = 1, 2
         call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json', cart=.true.)
         rc = trc_set_jkfit_json(hd, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
         rc = trc_set_jkfit_json(ho, 'basis_sets/def2-universal-jkfit.json'//c_null_char)
         rc = trc_set_df(hd, 1_c_int)
         rc = trc_set_df(ho, 1_c_int)
         if (k == 1) then  ! irrep_nelec first, then the request
            rc = trc_set_symmetry(hd, 'off'//c_null_char, 0.0_c_double)
            rc = trc_set_irrep_nelec(hd, 'A1=6,B1=2,B2=2'//c_null_char)
            rc = trc_set_symmetry(hd, 'on'//c_null_char, 0.0_c_double)
         else
            rc = trc_set_irrep_nelec(hd, 'A1=6,B1=2,B2=2'//c_null_char)
         end if
         rc = trc_run_scf(hd)
         call check(error, rc /= TRC_OK, 'irrep_nelec ('//merge('before', 'after ', k == 1)// &
                    ' on) with density fitting is refused, not run in C1')
         if (allocated(error)) return
         rc = trc_set_irrep_nelec(hd, ''//c_null_char)
         rc = trc_run_scf(hd)
         call check(error, rc, TRC_OK, 'cleared irrep_nelec: on is lenient again: '//message(hd))
         if (allocated(error)) return
         rc = trc_run_scf(ho)
         rc = trc_energy(hd, ed)
         rc = trc_energy(ho, eo)
         call check(error, abs(ed - eo) < 1.0e-11_dp, 'C1 fallback energy differs from off: '// &
                    to_string(abs(ed - eo), 'es10.2'))
         if (allocated(error)) return
         rc = trc_destroy(hd); rc = trc_destroy(ho)
      end do
   end subroutine test_irrep_nelec_order

   subroutine test_grid_after_c1_rerun(error)
      !! B3LYP on a rotated water: after an entry that reruns the SCF in C1
      !! (the polarizability, the Hessian), the gradient is the symmetry-off
      !! gradient. The C1 rerun builds its grid on the input axes, so a later
      !! grid for the gradient must not carry the symmetric run's rotation.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hd, ho
      real(c_double) :: gd(9), go(9), a(9), hs(81)
      integer(c_int) :: rc
      integer :: k

      do k = 1, 2
         call make_pair(hd, ho, 'water_rot.xyz', 'basis_sets/6-31gs.json')
         rc = trc_set_method(hd, 'b3lyp'//c_null_char, 3_c_int)
         rc = trc_set_method(ho, 'b3lyp'//c_null_char, 3_c_int)
         rc = trc_run_scf(hd)
         rc = trc_run_scf(ho)
         call check(error, name_of(hd, 1, 0), 'C2v', 'the B3LYP on run is symmetric')
         if (allocated(error)) return
         if (k == 1) then
            rc = trc_run_polarizability(hd)
            call check(error, rc, TRC_OK, 'polarizability: '//message(hd))
            rc = trc_run_polarizability(ho)
            rc = trc_polarizability(hd, a)
         else
            rc = trc_run_hessian(hd)
            call check(error, rc, TRC_OK, 'Hessian: '//message(hd))
            rc = trc_run_hessian(ho)
            rc = trc_hessian(hd, hs)
         end if
         if (allocated(error)) return
         rc = trc_run_gradient(hd)
         call check(error, rc, TRC_OK, 'gradient: '//message(hd))
         rc = trc_run_gradient(ho)
         if (allocated(error)) return
         rc = trc_gradient(hd, gd)
         rc = trc_gradient(ho, go)
         print '(a,a,es10.2)', '  gradient after the C1 rerun of the ', merge('polarizability', 'Hessian       ', k == 1), &
            maxval(abs(gd - go))
         call check(error, maxval(abs(gd - go)) < 1.0e-9_dp, 'gradient after '// &
                    merge('polarizability', 'Hessian       ', k == 1)//' differs from off: '// &
                    to_string(maxval(abs(gd - go)), 'es10.2'))
         if (allocated(error)) return
         rc = trc_destroy(hd); rc = trc_destroy(ho)
      end do
   end subroutine test_grid_after_c1_rerun

   subroutine compare_runs(label, hd, ho, error)
      !! SCF on both contexts: the same energy to 1e-11, no symmetry result on the on one.
      character(len=*), intent(in) :: label
      type(c_ptr), intent(in) :: hd, ho
      type(error_type), allocatable, intent(out) :: error
      real(c_double) :: ed, eo
      integer(c_int) :: rc
      rc = trc_run_scf(hd)
      call check(error, rc, TRC_OK, label//': the on context must run, not refuse: '//message(hd))
      rc = trc_run_scf(ho)
      if (allocated(error)) return
      rc = trc_energy(hd, ed)
      rc = trc_energy(ho, eo)
      print '(a,a,es10.2)', '  '//label//': |E(on) - E(off)| ', abs(ed - eo)
      call check(error, abs(ed - eo) < 1.0e-11_dp, label//': energy differs from symmetry off: '// &
                 to_string(abs(ed - eo), 'es10.2'))
      if (allocated(error)) return
      call check(error, name_of(hd, 1, 0), '<status 5>', label//': a C1 fallback leaves no symmetry result')
   end subroutine compare_runs

end module check_capi_sym_suite

program check_capi_sym
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_capi_sym_suite, only: collect_capi_sym
   implicit none
   integer :: stat
   stat = 0
   ! parallel=.false.: the basis reader is not thread-safe under -fopenmp.
   call run_testsuite(collect_capi_sym, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_sym
