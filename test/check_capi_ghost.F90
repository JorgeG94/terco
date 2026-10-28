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
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix, ref_read_int
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_elements, &
                               trc_set_basis_json, trc_set_basis_arrays, trc_set_method, &
                               trc_set_convergence, trc_set_pcm, trc_pcm_energy, &
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
                  new_unittest("PCM gives a ghost no sphere: same cavity as the bare monomer", &
                               test_pcm_ghost_no_sphere) &
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
      ! Independent grids, same tables: measured 3.95e-11 -- see the
      ! header above.
      call check(error, abs(real(e_terco, dp) - e_ref) < 1.0e-8_dp, &
                 "B3LYP energy (grids on the ghosts, independent quadratures) close to pyscf")
      rc = trc_destroy(h)
   end subroutine test_b3lyp_vs_pyscf

   !
   ! PCM with ghosts. pyscf's own PCM (pyscf.solvent.pcm) does NOT give a
   ! ghost atom an empty cavity: checked separately (not here -- no
   ! pyscf import on this side of the boundary), a ghost gets the SAME
   ! sphere a real atom of its element would, built from its own
   ! radius -- confirmed by building the dimer-with-ghost cavity in
   ! pyscf directly and finding non-zero Lebedev points on the ghost's
   ! three atoms, with the real monomer's OWN points and weights
   ! unchanged by the ghost's presence only because the two monomers
   ! are far enough apart (4 Angstrom, this geometry) that neither
   ! sphere's switching function reaches the other -- not because
   ! pyscf excludes the ghost. That is a DIFFERENT model from the one
   ! this feature's design agreed on ("no sphere and no v_nuc for
   ! ghosts"), so pyscf's PCM energy is not a reference for terco's
   ! here.
   !
   ! What IS tested: terco's own monomer-in-the-dimer's-basis cavity
   ! has exactly as many points as the bare monomer's -- the ghost
   ! atoms (4-6) contribute none, and the real atoms' (1-3) spheres are
   ! built the same both times (same geometry, same radii) -- against a
   ! REAL (non-ghost) dimer, which must have MORE points than either,
   ! as a sanity check that this is not merely "PCM always returns the
   ! same count regardless of what is handed to it."
   !
   subroutine test_pcm_ghost_no_sphere(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, npts_bare, npts_ghost, npts_real
      real(c_double) :: ediel, qtot

      ! (a) the bare monomer: atoms 1-3 only, real.
      rc = trc_create(h)
      rc = trc_set_molecule(h, 3_c_int, real(Z6(1:3), c_double), reshape(XYZ6(:, 1:3), [9]), &
                            0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      rc = trc_set_pcm(h, 'cpcm'//c_null_char, 78.3553_c_double, 0_c_int, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "run_scf (bare monomer, PCM)")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_pcm_energy(h, ediel, qtot, npts_bare)
      rc = trc_destroy(h)

      ! (b) monomer A real, monomer B ghost: this feature's own case.
      call make_molecule(h, with_elements=.true.)
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      rc = trc_set_pcm(h, 'cpcm'//c_null_char, 78.3553_c_double, 0_c_int, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "run_scf (monomer-in-dimer-basis, PCM)")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_pcm_energy(h, ediel, qtot, npts_ghost)
      rc = trc_destroy(h)

      ! (c) the REAL dimer: both monomers real, no ghosts -- more points
      ! than (a)/(b), the sanity check that this is testing the cavity
      ! and not a constant.
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(NATM, c_int), &
                            [8.0_c_double, 1.0_c_double, 1.0_c_double, 8.0_c_double, 1.0_c_double, 1.0_c_double], &
                            reshape(XYZ6, [3*NATM]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      rc = trc_set_pcm(h, 'cpcm'//c_null_char, 78.3553_c_double, 0_c_int, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "run_scf (real dimer, PCM)")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_pcm_energy(h, ediel, qtot, npts_real)
      rc = trc_destroy(h)

      print '(a,i0,a,i0,a,i0)', "check_capi_ghost: PCM points bare monomer ", npts_bare, &
         "  monomer-in-dimer-basis(ghost) ", npts_ghost, "  real dimer ", npts_real
      call check(error, npts_ghost, npts_bare, &
                 "a ghost gets no sphere: same cavity point count as the bare monomer")
      if (allocated(error)) return
      call check(error, npts_real > npts_bare, &
                 "a REAL dimer has more cavity points than the bare monomer (sanity)")
   end subroutine test_pcm_ghost_no_sphere

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
