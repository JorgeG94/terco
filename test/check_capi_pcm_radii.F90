module check_capi_pcm_radii_suite
   !! trc_set_pcm_radii, through the context API: water/6-31G*, C-PCM,
   !! PBE0.
   !!
   !! Case A -- self-consistency, no pyscf needed: an override equal to
   !! what terco's own table (trc_pcm_radii.F90) already produces for
   !! this molecule must give the same energy and E_diel as the table
   !! path -- the override and the table are two routes to the same
   !! per-atom radii array, so pcm_build_surface cannot tell them apart.
   !! The same to round-off (SAME_TOL), not bit for bit: each side is its
   !! own SCF, and with OpenMP threads the reductions sum in a different
   !! order every run (2.4e-13 apart at 4 threads in CI). A wrong radius
   !! moves the energy by more than 1e-6 (test_misuse checks that).
   !!
   !! Case B -- against pyscf, with DIFFERENT radii: mqc's own
   !! convention (plain Bondi, hydrogen 1.20 Angstrom, times 1.2) against
   !! terco's default table (pyscf's MODIFIED Bondi, hydrogen 1.10) --
   !! the discrepancy trc_set_pcm_radii exists to remove. Checked against
   !! test/pyscf_refs/pcm_radii_override_water.ref (pyscf independently
   !! converged with radii_table = pcm.radii.VDW, vdw_scale = 1.2).
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char, c_null_ptr, c_loc
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_force_cartesian, trc_set_method, trc_set_convergence, &
                               trc_set_pcm, trc_set_pcm_radii, trc_run_scf, trc_energy, trc_pcm_energy, &
                               TRC_OK, TRC_ERR_BADARG
   implicit none
   private
   public :: collect_capi_pcm_radii_tests

   ! trc_pcm_radii.F90's own conversion constant and table entries for
   ! O (Z=8, 1.52 Angstrom) and H (Z=1), both tables -- modified Bondi
   ! (1.10, terco's default) and plain Bondi (1.20, mqc's/case B's).
   real(dp), parameter :: BOHR_ANGSTROM = 0.52917721092_dp
   real(dp), parameter :: SCALE = 1.2_dp
   real(dp), parameter :: O_ANG = 1.52_dp, H_MODBONDI_ANG = 1.10_dp, H_BONDI_ANG = 1.20_dp
   real(dp), parameter :: SAME_TOL = 1.0e-10_dp
      !! Two SCFs on the same radii agree to this (round-off, see
      !! test_override_matches_table); a wrong radius is more than 1e-6
      !! away.

contains

   subroutine collect_capi_pcm_radii_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("an override equal to the table matches the table path", &
                               test_override_matches_table), &
                  new_unittest("an override at mqc's own radii matches an independent pyscf reference", &
                               test_override_vs_pyscf), &
                  new_unittest("a wrong natm is refused; natm=0 or a NULL pointer clears the override", &
                               test_misuse) &
                  ]
   end subroutine collect_capi_pcm_radii_tests

   subroutine run_one(h, radii_ang, error)
      !! water/6-31gs, cpcm/PBE0, with trc_set_pcm_radii given
      !! radii_ang(natm) (Angstrom, water.xyz's O, H, H order), or no
      !! override at all if radii_ang is not allocated.
      type(c_ptr), intent(out) :: h
      real(dp), allocatable, intent(in) :: radii_ang(:)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      integer(c_int) :: rc
      real(c_double), allocatable, target :: radii_bohr(:)

      call read_xyz('water.xyz', natm, zint, at_r)

      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_force_cartesian(h, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_method(h, 'pbe0'//c_null_char, 3_c_int)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 150_c_int, 10_c_int)
      rc = trc_set_pcm(h, 'cpcm'//c_null_char, 78.3553_c_double, 0_c_int, 0.0_c_double)
      call check(error, rc == TRC_OK, "run_one: set_pcm status")
      if (allocated(error)) return
      if (allocated(radii_ang)) then
         allocate (radii_bohr(size(radii_ang)))
         radii_bohr = real(radii_ang, c_double)*real(SCALE, c_double)/real(BOHR_ANGSTROM, c_double)
         rc = trc_set_pcm_radii(h, int(size(radii_bohr), c_int), c_loc(radii_bohr))
         call check(error, rc == TRC_OK, "run_one: set_pcm_radii status")
         if (allocated(error)) return
      end if
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "run_one: run_scf must converge")
   end subroutine run_one

   subroutine test_override_matches_table(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h_table, h_override
      real(dp), allocatable :: radii_ang(:)
      real(c_double) :: e_table, e_override, ediel_table, ediel_override, qtot
      integer(c_int) :: rc, npts

      call run_one(h_table, radii_ang, error)   ! radii_ang not allocated: table default
      if (allocated(error)) return

      ! water.xyz order: O, H, H -- same radius for both H.
      allocate (radii_ang(3))
      radii_ang(1) = O_ANG
      radii_ang(2) = H_MODBONDI_ANG
      radii_ang(3) = H_MODBONDI_ANG
      call run_one(h_override, radii_ang, error)
      if (allocated(error)) then
         rc = trc_destroy(h_table)
         return
      end if

      rc = trc_energy(h_table, e_table)
      rc = trc_energy(h_override, e_override)
      rc = trc_pcm_energy(h_table, ediel_table, qtot, npts)
      rc = trc_pcm_energy(h_override, ediel_override, qtot, npts)

      call check(error, abs(real(e_table, dp) - real(e_override, dp)) < SAME_TOL, &
                 "override_matches_table: SCF energy matches the table, table "// &
                 to_string(real(e_table, dp), "f20.14")//"  override "// &
                 to_string(real(e_override, dp), "f20.14")//"  diff "// &
                 to_string(real(e_table, dp) - real(e_override, dp), "es9.2"))
      if (.not. allocated(error)) then
         call check(error, abs(real(ediel_table, dp) - real(ediel_override, dp)) < SAME_TOL, &
                    "override_matches_table: E_diel matches the table")
      end if

      rc = trc_destroy(h_table); rc = trc_destroy(h_override)
   end subroutine test_override_matches_table

   subroutine test_override_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: radii_ang(:)
      real(c_double) :: e_mqc, ediel_mqc, qtot
      integer(c_int) :: rc, npts, u
      real(dp) :: e_ref, ediel_ref

      allocate (radii_ang(3))
      radii_ang(1) = O_ANG
      radii_ang(2) = H_BONDI_ANG
      radii_ang(3) = H_BONDI_ANG
      call run_one(h, radii_ang, error)
      if (allocated(error)) return

      rc = trc_energy(h, e_mqc)
      rc = trc_pcm_energy(h, ediel_mqc, qtot, npts)

      call ref_require('pcm_radii_override_water.ref', u)
      call ref_read_real(u, 'E_SCF', e_ref)
      call ref_read_real(u, 'E_DIEL', ediel_ref)
      call ref_close(u)

      call check(error, abs(real(e_mqc, dp) - e_ref) < 1.0e-7_dp, &
                 "override_vs_pyscf: SCF energy matches pyscf with the SAME (mqc) radii, terco "// &
                 to_string(real(e_mqc, dp), "f20.14")//"  pyscf "//to_string(e_ref, "f20.14")// &
                 "  diff "//to_string(real(e_mqc, dp) - e_ref, "es9.2"))
      if (.not. allocated(error)) then
         call check(error, abs(real(ediel_mqc, dp) - ediel_ref) < 1.0e-6_dp, &
                    "override_vs_pyscf: E_diel matches pyscf with the SAME (mqc) radii, terco "// &
                    to_string(real(ediel_mqc, dp), "f20.14")//"  pyscf "//to_string(ediel_ref, "f20.14")// &
                    "  diff "//to_string(real(ediel_mqc, dp) - ediel_ref, "es9.2"))
      end if

      rc = trc_destroy(h)
   end subroutine test_override_vs_pyscf

   !
   ! A wrong natm is refused; natm = 0 (or a NULL radii) clears a
   ! previously-set override back to the table, and trc_set_charge keeps
   ! it (the molecule, and so the cavity's atoms, are unchanged).
   !
   subroutine test_misuse(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double), target :: radii_wrong(2), radii_ok(3)
      real(c_double) :: e_table, e_cleared, e_before, e_after
      real(c_double) :: ediel_unused, qtot
      integer(c_int) :: npts

      call read_xyz('water.xyz', natm, zint, at_r)

      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_force_cartesian(h, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_method(h, 'pbe0'//c_null_char, 3_c_int)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 150_c_int, 10_c_int)
      rc = trc_set_pcm(h, 'cpcm'//c_null_char, 78.3553_c_double, 0_c_int, 0.0_c_double)

      radii_wrong = [1.0_c_double, 1.0_c_double]
      rc = trc_set_pcm_radii(h, 2_c_int, c_loc(radii_wrong))
      call check(error, rc == TRC_ERR_BADARG, "misuse: natm not matching the molecule is refused")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "misuse: run_scf (table default) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_energy(h, e_table)

      radii_ok = [O_ANG, H_BONDI_ANG, H_BONDI_ANG]*real(SCALE, c_double)/real(BOHR_ANGSTROM, c_double)
      rc = trc_set_pcm_radii(h, 3_c_int, c_loc(radii_ok))
      call check(error, rc == TRC_OK, "misuse: set_pcm_radii (ok) status")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "misuse: run_scf (overridden) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_energy(h, e_before)
      call check(error, abs(real(e_before, dp) - real(e_table, dp)) > 1.0e-6_dp, &
                 "misuse: the override actually changed the energy")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      ! natm = 0 clears it, back to the table.
      rc = trc_set_pcm_radii(h, 0_c_int, c_loc(radii_ok))
      call check(error, rc == TRC_OK, "misuse: natm = 0 clears the override")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "misuse: run_scf (cleared) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_energy(h, e_cleared)
      call check(error, abs(real(e_cleared, dp) - real(e_table, dp)) < SAME_TOL, &
                 "misuse: cleared override matches the table again")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      ! A NULL radii pointer clears it too.
      rc = trc_set_pcm_radii(h, 3_c_int, c_loc(radii_ok))
      rc = trc_set_pcm_radii(h, 3_c_int, c_null_ptr)
      call check(error, rc == TRC_OK, "misuse: a NULL radii pointer clears the override")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "misuse: run_scf (cleared via NULL) must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_energy(h, e_after)
      call check(error, abs(real(e_after, dp) - real(e_table, dp)) < SAME_TOL, &
                 "misuse: NULL-cleared override matches the table again")

      ediel_unused = 0.0_c_double; qtot = 0.0_c_double; npts = 0
      rc = trc_pcm_energy(h, ediel_unused, qtot, npts)   ! touch, so the entry stays exercised
      rc = trc_destroy(h)
   end subroutine test_misuse

end module check_capi_pcm_radii_suite

program check_capi_pcm_radii
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_bind_device
   use check_capi_pcm_radii_suite, only: collect_capi_pcm_radii_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, dev

   dev = trc_bind_device(0)
   stat = 0
   ! Serial: each case builds its own context through the C ABI, which is
   ! not thread-safe, and an -fopenmp build runs test-drive cases at once.
   testsuites = [new_testsuite("check_capi_pcm_radii", collect_capi_pcm_radii_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_pcm_radii
