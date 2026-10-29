module check_capi_ghost_suite
   !! Ghost atoms through the C ABI (trc_set_elements, ABI 1.7): mqc's own
   !! counterpoise use. A GHOST is z == 0 (trc_set_molecule's own nuclear
   !! charge) with element >= 1 (trc_set_elements) -- basis functions sit
   !! there, but there is no nucleus and no electron.
   !!
   !! Geometry: a water dimer, cc-pVDZ, monomer A real (atoms 1-3) and
   !! monomer B placed as ghosts (atoms 4-6, element O/H/H but z = 0) --
   !! monomer A run in the DIMER's basis, the textbook counterpoise
   !! fragment. Coordinates are literal BOHR, computed from the Angstrom
   !! geometry with terco's own CODATA 2018 Bohr constant and duplicated
   !! by hand into test/pyscf_refs/generate.py's gen_ghost_water_dimer
   !! (same convention as gen_opt/gen_vib_thermo) -- so there is no
   !! Angstrom -> Bohr conversion constant for the two sides to disagree
   !! about at the last bit.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char
   use pic_types, only: dp
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: pack_libcint
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix, ref_read_int
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_elements, &
                               trc_set_basis_json, trc_set_basis_arrays, trc_set_basis_libcint, &
                               trc_set_method, &
                               trc_set_convergence, trc_set_pcm, &
                               trc_run_scf, trc_run_gradient, trc_gradient, &
                               trc_energy, TRC_OK, TRC_ERR_BADARG, TRC_ERR_STATE
   implicit none
   private
   public :: collect_capi_ghost

   integer, parameter :: NATM = 6
   !! water dimer: monomer A (real) then monomer B (ghost)
   real(dp), parameter :: Z6(NATM) = [8.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   integer, parameter :: ELEM6(NATM) = [8, 1, 1, 8, 1, 1]
   real(dp), parameter :: XYZ6(3, NATM) = reshape([ &
      0.0_dp, 0.0_dp, 0.221664874418603_dp, &
      0.0_dp, 1.430900621566633_dp, -0.886659497674411_dp, &
      0.0_dp, -1.430900621566633_dp, -0.886659497674411_dp, &
      7.558904498503081_dp, 0.0_dp, 0.221664874418603_dp, &
      7.558904498503081_dp, 1.430900621566633_dp, -0.886659497674411_dp, &
      7.558904498503081_dp, -1.430900621566633_dp, -0.886659497674411_dp], [3, NATM])

contains

   subroutine collect_capi_ghost(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("element 0 at run time is refused", test_element_zero), &
                  new_unittest("an all-ghost molecule is refused", test_all_ghost), &
                  new_unittest("trc_set_elements after a JSON basis is refused", test_after_basis), &
                  new_unittest("trc_set_elements natm mismatch is refused", test_natm_mismatch), &
                  new_unittest("RHF energy and gradient match pyscf's ghost atoms", test_rhf_vs_pyscf), &
                  new_unittest("B3LYP energy (grids on the ghosts) is close to pyscf", test_b3lyp_vs_pyscf), &
                  new_unittest("B3LYP through trc_set_basis_libcint keeps the ghosts' elements", &
                               test_b3lyp_libcint), &
                  new_unittest("PCM with ghost atoms is refused (cavity model undecided)", &
                               test_pcm_ghost_refused) &
                  ]
   end subroutine collect_capi_ghost

   !
   ! A context with the dimer molecule set (z above), optionally with
   ! trc_set_elements(ELEM6) already called -- the ONE place every test
   ! below builds the molecule, so the geometry/elements never drift
   ! between cases.
   !
   subroutine make_molecule(h, with_elements)
      type(c_ptr), intent(out) :: h
      logical, intent(in) :: with_elements
      integer(c_int) :: rc
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(NATM, c_int), real(Z6, c_double), reshape(XYZ6, [3*NATM]), &
                            0_c_int, 1_c_int)
      if (with_elements) rc = trc_set_elements(h, int(NATM, c_int), int(ELEM6, c_int))
   end subroutine make_molecule

   subroutine test_element_zero(error)
      !! No trc_set_elements call, and a caller passes z = 0 for the
      !! second atom by mistake (not through trc_set_elements, so no
      !! deliberate ghost): trc_set_molecule's own default
      !! (element := nint(z)) leaves it at element 0 -- the bug ghost
      !! atoms would have reopened had trc_run_scf not refused it. A
      !! trc_set_basis_arrays basis (placed explicitly, no JSON element
      !! lookup) reaches trc_run_scf rather than being refused earlier
      !! by the JSON reader's own "element not in <path>" error.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      rc = trc_create(h)
      rc = trc_set_molecule(h, 2_c_int, [1.0_c_double, 0.0_c_double], &
                            [0.0_c_double, 0.0_c_double, 0.0_c_double, &
                             0.0_c_double, 0.0_c_double, 1.4_c_double], 0_c_int, 1_c_int)
      rc = trc_set_basis_arrays(h, 2_c_int, 1_c_int, [0_c_int, 0_c_int], [1_c_int, 1_c_int], &
                                [1.0_c_double, 1.0_c_double], [1.0_c_double, 1.0_c_double], &
                                [0.0_c_double, 0.0_c_double, 0.0_c_double, &
                                 0.0_c_double, 0.0_c_double, 1.4_c_double])
      call check(error, rc, TRC_OK, "set_basis_arrays (no element lookup, so element 0 does not stop it)")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_BADARG, "run_scf refuses element 0")
      rc = trc_destroy(h)
   end subroutine test_element_zero

   subroutine test_all_ghost(error)
      !! Every atom's z is 0 (nothing left to compute), with valid
      !! elements (1..118) so this is NOT the element-0 refusal above.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: z0(NATM)
      z0 = 0.0_c_double
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(NATM, c_int), z0, reshape(XYZ6, [3*NATM]), 0_c_int, 1_c_int)
      rc = trc_set_elements(h, int(NATM, c_int), int(ELEM6, c_int))
      call check(error, rc, TRC_OK, "set_elements (valid 1..118 everywhere)")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      call check(error, rc, TRC_OK, "set_basis_json")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_BADARG, "run_scf refuses an all-ghost molecule")
      rc = trc_destroy(h)
   end subroutine test_all_ghost

   subroutine test_after_basis(error)
      !! trc_set_elements must run before the JSON basis is loaded --
      !! the basis is looked up by element, so changing it afterwards
      !! would leave shells on the wrong atoms. An ordinary molecule
      !! (no ghosts) here: this is purely the ORDERING refusal, not the
      !! element-0 one above.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      rc = trc_create(h)
      rc = trc_set_molecule(h, 2_c_int, [1.0_c_double, 1.0_c_double], &
                            [0.0_c_double, 0.0_c_double, 0.0_c_double, &
                             0.0_c_double, 0.0_c_double, 1.4_c_double], 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      call check(error, rc, TRC_OK, "set_basis_json")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_set_elements(h, 2_c_int, [1_c_int, 1_c_int])
      call check(error, rc, TRC_ERR_STATE, "set_elements after set_basis_json is refused")
      rc = trc_destroy(h)
   end subroutine test_after_basis

   subroutine test_natm_mismatch(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      call make_molecule(h, with_elements=.false.)
      rc = trc_set_elements(h, int(NATM - 1, c_int), int(ELEM6(1:NATM - 1), c_int))
      call check(error, rc, TRC_ERR_BADARG, "set_elements natm mismatch is refused")
      rc = trc_destroy(h)
   end subroutine test_natm_mismatch

   !
   ! RHF: analytic throughout (no quadrature grid on either side), so
   ! this is a tight, ~1e-10 match. The gradient has all 6 rows -- the
   ! ghost ones carry only basis-function (Pulay) terms, no nuclear
   ! repulsion/attraction term, since z = 0 there on both sides.
   !
   subroutine test_rhf_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      integer :: u
      real(c_double) :: e_terco
      real(c_double), allocatable :: grad(:)
      real(dp) :: e_ref
      real(dp) :: grad_ref(3, NATM)
      real(dp) :: worst
      integer :: i, natm_ref

      call make_molecule(h, with_elements=.true.)
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "run_scf (RHF, ghost dimer)")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_energy(h, e_terco)

      rc = trc_run_gradient(h)
      call check(error, rc, TRC_OK, "run_gradient")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      allocate (grad(3*NATM))
      rc = trc_gradient(h, grad)

      call ref_require('ghost_water_dimer.ref', u)
      call ref_read_int(u, 'NATM', natm_ref)
      call ref_read_real(u, 'E_RHF', e_ref)
      call ref_read_matrix(u, 'GRAD_RHF', 3, NATM, grad_ref)
      call ref_close(u)

      call check(error, natm_ref, NATM, "the reference's own atom count")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if

      print '(a,f20.14,a,f20.14,a,es9.2)', "check_capi_ghost: E_RHF terco ", real(e_terco, dp), &
         "  pyscf ", e_ref, "  diff ", real(e_terco, dp) - e_ref
      call check(error, abs(real(e_terco, dp) - e_ref) < 1.0e-9_dp, "RHF energy matches pyscf to 1e-9")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if

      worst = 0.0_dp
      do i = 1, 3*NATM
         worst = max(worst, abs(real(grad(i), dp) - grad_ref(mod(i - 1, 3) + 1, (i - 1)/3 + 1)))
      end do
      print '(a,es9.2)', "check_capi_ghost: RHF gradient worst component diff ", worst
      call check(error, worst < 1.0e-8_dp, "RHF gradient (all 6 rows, ghosts included) matches pyscf")

      rc = trc_destroy(h)
   end subroutine test_rhf_vs_pyscf

   !
   ! B3LYP: both sides build their OWN quadrature grid independently
   ! (same radial/angular tables, Becke partition -- not the same
   ! points, since trc's C ABI exposes no grid-injection entry). Measured
   ! 3.95e-11 (gfortran 15, CPU, level 6): the two grids are built from
   ! the SAME tables (trc_dft_grid.F90's own header) closely enough that
   ! this reaches the ~1e-10 class, not merely a loose quadrature-limited
   ! bound -- the gate below is set from that measurement, with headroom,
   ! not asserted blind. pyscf's own default grid DOES put points on a
   ! ghost atom (checked separately against its PCM cavity, which does
   ! NOT exclude ghosts the way terco's does -- see the report), so this
   ! is comparing two grids that both cover the ghost, not one that does
   ! against one that does not.
   !
   subroutine test_b3lyp_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      integer :: u, lvl, natm_ref
      real(c_double) :: e_terco
      real(dp) :: e_ref
      real(dp) :: grad_ref(3, NATM)

      call make_molecule(h, with_elements=.true.)
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      rc = trc_set_method(h, 'b3lyp'//c_null_char, 6_c_int)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "run_scf (B3LYP, ghost dimer)")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_energy(h, e_terco)

      ! Same file as test_rhf_vs_pyscf: read past its earlier tags (in
      ! the order generate.py's RefWriter wrote them) to reach E_B3LYP.
      call ref_require('ghost_water_dimer.ref', u)
      call ref_read_int(u, 'NATM', natm_ref)
      call ref_read_real(u, 'E_RHF', e_ref)
      call ref_read_matrix(u, 'GRAD_RHF', 3, NATM, grad_ref)
      call ref_read_int(u, 'B3LYP_GRID_LEVEL', lvl)
      call ref_read_real(u, 'E_B3LYP', e_ref)
      call ref_close(u)

      print '(a,f20.14,a,f20.14,a,es9.2)', "check_capi_ghost: E_B3LYP terco ", real(e_terco, dp), &
         "  pyscf ", e_ref, "  diff ", real(e_terco, dp) - e_ref
      ! Independent grids, same tables: measured 3.95e-11 on gfortran (CPU)
      ! and 3.96e-11 on nvfortran (V100) -- see the header above.
      call check(error, abs(real(e_terco, dp) - e_ref) < 1.0e-10_dp, &
                 "B3LYP energy (grids on the ghosts, independent quadratures) close to pyscf")
      rc = trc_destroy(h)
   end subroutine test_b3lyp_vs_pyscf

   !
   ! The same ghosted B3LYP, loaded the way mqc loads it: trc_set_molecule
   ! (ghosts at charge 0), trc_set_elements, then trc_set_basis_libcint,
   ! whose atm carries only the charge. That path built the basis without the
   ! elements, so every ghost was element 0 -- an element-0 DFT grid and Becke
   ! partition at a ghost oxygen, 1e-5 to 1e-4 Eh off in mqc's VMFC rows while
   ! RHF agreed. Same pyscf reference and gate as test_b3lyp_vs_pyscf.
   !
   subroutine test_b3lyp_libcint(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      type(trc_basis_t) :: b
      type(error_t) :: berr
      integer(c_int), allocatable :: atm(:, :), bas(:, :)
      real(c_double), allocatable :: env(:)
      integer(c_int) :: rc
      integer :: u, lvl, natm_ref
      real(c_double) :: e_terco
      real(dp) :: e_ref
      real(dp) :: grad_ref(3, NATM)

      ! The dimer basis from the JSON file, then packed into libcint's arrays
      ! with the ghosts' CHARGE_OF at 0 (Z6), as a libcint host passes it.
      call trc_basis_from_json('basis_sets/cc-pvdz.json', NATM, ELEM6, XYZ6, b, berr, charges=Z6)
      call check(error, .not. berr%has_error(), "dimer basis from JSON")
      if (allocated(error)) return
      call pack_libcint(b, atm, bas, env)

      call make_molecule(h, with_elements=.true.)
      rc = trc_set_basis_libcint(h, atm, int(NATM, c_int), bas, int(size(bas, 2), c_int), env, &
                                 int(size(env), c_int), merge(0_c_int, 1_c_int, b%spherical))
      call check(error, rc, TRC_OK, "set_basis_libcint (ghosts at charge 0)")
      call b%release()
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_set_method(h, 'b3lyp'//c_null_char, 6_c_int)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "run_scf (B3LYP, ghost dimer, libcint basis)")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_energy(h, e_terco)

      call ref_require('ghost_water_dimer.ref', u)
      call ref_read_int(u, 'NATM', natm_ref)
      call ref_read_real(u, 'E_RHF', e_ref)
      call ref_read_matrix(u, 'GRAD_RHF', 3, NATM, grad_ref)
      call ref_read_int(u, 'B3LYP_GRID_LEVEL', lvl)
      call ref_read_real(u, 'E_B3LYP', e_ref)
      call ref_close(u)

      print '(a,f20.14,a,f20.14,a,es9.2)', "check_capi_ghost: E_B3LYP libcint terco ", real(e_terco, dp), &
         "  pyscf ", e_ref, "  diff ", real(e_terco, dp) - e_ref
      call check(error, abs(real(e_terco, dp) - e_ref) < 1.0e-10_dp, &
                 "B3LYP through trc_set_basis_libcint matches pyscf's ghost atoms")
      rc = trc_destroy(h)
   end subroutine test_b3lyp_libcint

   !
   ! PCM with ghosts is refused until a cavity model is chosen. pyscf's own
   ! PCM (pyscf.solvent.pcm), and mqc's CPU backend after it, give a ghost
   ! the same sphere a real atom of its element gets; trc_pcm's ghost mask
   ! would give it none (measured on this geometry before the refusal: 556
   ! cavity points for monomer A in the dimer basis, 556 bare, 1112 for the
   ! real dimer). The two models differ, so neither is offered yet.
   !
   subroutine test_pcm_ghost_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc

      call make_molecule(h, with_elements=.true.)
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      rc = trc_set_pcm(h, 'cpcm'//c_null_char, 78.3553_c_double, 0_c_int, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc /= TRC_OK, "run_scf refuses PCM with ghost atoms")
      rc = trc_destroy(h)
   end subroutine test_pcm_ghost_refused

end module check_capi_ghost_suite

program check_capi_ghost
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_capi_ghost_suite, only: collect_capi_ghost
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_capi_ghost, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_ghost
