!
! The symmetry entries of the C ABI (ABI 1.22; 1.23 adds trc_set_symmetry_skeleton), called through the Fortran
! interfaces of include/trc_c_interfaces.f90: trc_set_symmetry,
! trc_set_irrep_nelec and the getters. The C side of the same round trip (water,
! the header compiled by a C compiler) is in link_c.c.
!
! HOCl (planar, Cs) 6-31G*. The oracle is the project's central one: the symmetric
! run and a plain run on the standard-frame geometry the context reports
! (trc_geometry) give the same energy, to 1e-10 Eh. Then the frame (an
! orthogonal map from the input to that geometry), the labels, fixed irrep
! occupations (and the occupations trc_mo_occupations reports for them), and what
! is refused, each with its message: every entry that would otherwise run the
! C1 path silently, the cases the SCF cannot do with symmetry, and an explicit
! group that is not a subgroup or not supported yet.
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
                  new_unittest("entries that would run C1 silently are refused", test_entry_refusals) &
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
      rc = trc_set_symmetry(h, 'auto'//c_null_char, 0.0_c_double)
      call check(error, rc, TRC_OK, 'trc_set_symmetry(auto)')
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
      call check(error, index(message(h), 'off, auto or a point group') > 0, 'its message: '//message(h))
      rc = trc_set_symmetry(h, 'auto'//c_null_char, -1.0_c_double)
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
      rc = trc_set_symmetry(h, 'auto'//c_null_char, 0.0_c_double)
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
      rc = trc_set_symmetry(h, 'auto'//c_null_char, 0.0_c_double)
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
      !! The skeleton Fock build is the default under a symmetry request; with
      !! it switched off (the full build, then symmetrized) the energy and the
      !! gradient are the same, and the gradient that follows a skeleton SCF
      !! does not reuse its orbit-closed engine.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: geom(:, :)
      real(c_double) :: e1, e2, g1(9), g2(9)
      integer(c_int) :: rc

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'auto'//c_null_char, 0.0_c_double)
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
      integer(c_int) :: rc, nao

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      ! irrep_nelec without a symmetry request is a mistake, not a no-op.
      rc = trc_set_irrep_nelec(h, "A'=20,A"""//'=6'//c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_BADARG, 'irrep_nelec needs symmetry')
      call check(error, index(message(h), 'needs a symmetry request') > 0, message(h))
      if (allocated(error)) return
      rc = trc_set_irrep_nelec(h, c_null_char)
      rc = trc_set_symmetry(h, 'auto'//c_null_char, 0.0_c_double)
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
      rc = trc_set_symmetry(h, 'auto'//c_null_char, 0.0_c_double)
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
      rc = trc_set_symmetry(h, 'auto'//c_null_char, 0.0_c_double)
      rc = trc_set_soscf(h, 1_c_int, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc /= TRC_OK, 'SOSCF is refused with symmetry')
      call check(error, index(message(h), 'soscf') > 0, message(h))
      rc = trc_destroy(h)
      if (allocated(error)) return

      call make_context(h, 'hocl.xyz', 'basis_sets/6-31gs.json', geom)
      rc = trc_set_symmetry(h, 'auto'//c_null_char, 0.0_c_double)
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
      rc = trc_set_symmetry(h, 'auto'//c_null_char, 0.0_c_double)
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
