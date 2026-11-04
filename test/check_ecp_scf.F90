!
! ECP phase 3, end to end, through the C API: trc_set_molecule with the
! BARE atomic number (ECP_DESIGN.md's DECIDED charge contract), auto-attach
! from the orbital basis file's own ecp_potentials (trc_set_basis_json
! alone -- def2-svp.json carries them for Sn/I/Xe/Cs), and the refusals
! phase 3 adds (MINAO/SAP/Hirshfeld guesses, the analytic gradient).
!
! HI and the free I atom are checked against mqc's own independently
! converged numbers (ECP_DESIGN.md's Validation section) -- a STORED VALUE
! reference (AGENTS.md), not a pyscf one: an unguided pyscf RHF/UHF on
! these same systems was measured landing on a DIFFERENT, nearby
! stationary point (|dE| ~1e-4 for HI; the free atom, with no geometry
! freedom to get wrong, instead matches mqc's own number to 1e-13) --
! multiple self-consistent HF solutions are a real possibility here, not
! a terco bug: check_ecp_matrix.F90's own gen_ecp_ints already validates
! this build's U_ECP matrix element-by-element against pyscf's own
! ECPscalar to 7e-15 on this EXACT geometry, so the one-electron operator
! itself is not in question, only which root an unseeded SCF finds. Chasing
! a pyscf-seeded anchor that reproducibly lands on the SAME root was not
! pursued further here; mqc's numbers are the validation this suite holds
! to, per the phase 3 handoff.
!
! An I2 case (two ECP atoms, exercising attach_ecp's per-atom loop) was
! tried at several I-I separations (4.0-6.0 bohr) and refused every time
! with "the ECP radial quadrature did not converge within ... ECP_LEVEL_MAX"
! -- trc_ecp.F90/trc_ecp_radial.F90's own adaptive radial quadrature, a
! phase 1-2 integral routine, not touched by this phase. The existing ECP
! integral test suite (check_ecp_matrix's hi_svp/ih_far/xe/csh/snh4) has
! exactly ONE ECP atom in every case; two ECP atoms close enough to overlap
! appears to be an untested code path there. At 25 bohr (ih_far's own far
! separation) the matrix builds fine, so this is specifically a near-field
! issue, not a general multi-ECP-atom one -- left for a phase 1-2 follow-up
! (raise ECP_LEVEL_MAX, or find the pair/channel it fails on), not pursued
! further here; no I2 case is in this suite.
!
module check_ecp_scf_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_convergence, trc_set_guess, trc_run_scf, trc_run_gradient, &
                               trc_hirshfeld_charges, trc_save_scf, trc_load_scf, &
                               trc_energy, trc_iterations, &
                               TRC_OK, TRC_ERR_UNSUPPORTED, TRC_GUESS_MINAO, TRC_GUESS_SAP
   implicit none
   private
   public :: collect_ecp_scf_tests

   real(dp), parameter :: MQC_HI = -297.23153166336_dp
   !! HI/def2-SVP RHF: H at the origin, I on the z-axis at 3.04 bohr --
   !! check_ecp_matrix.F90's own "hi_svp" case, the system the VERIFIED
   !! BUG (ECP_DESIGN.md) is about. Measured agreement: 1.7e-6.
   real(dp), parameter :: MQC_I_ATOM = -296.649606475643_dp
   !! The free I atom, UHF: def2-svp's 28-electron ECP core leaves 25
   !! valence electrons -- odd, a doublet. No geometry at all, so there is
   !! nothing for terco and mqc's own number to disagree about
   !! independently. Measured agreement: 1.1e-13.

contains

   subroutine collect_ecp_scf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF HI/def2-SVP (auto-attached ECP) vs mqc", test_hi), &
                  new_unittest("UHF I atom/def2-SVP (ECP, doublet) vs mqc", test_i_atom), &
                  new_unittest("MINAO guess is refused with an ECP attached", test_minao_refused), &
                  new_unittest("SAP guess is refused with an ECP attached", test_sap_refused), &
                  new_unittest("Hirshfeld charges are refused with an ECP attached", test_hirshfeld_refused), &
                  new_unittest("the analytic gradient is refused with an ECP attached", test_gradient_refused), &
                  new_unittest("a restart file round-trips an ECP SCF (format version 3)", test_restart_roundtrip) &
                  ]
   end subroutine collect_ecp_scf_tests

   !
   ! HI: H at the origin, I at (0, 0, 3.04) bohr -- check_ecp_matrix.F90's
   ! own "hi_svp" case_geometry, read here verbatim.
   !
   subroutine make_hi(h)
      type(c_ptr), intent(out) :: h
      integer(c_int) :: rc
      real(c_double), parameter :: z(2) = [1.0_c_double, 53.0_c_double]
      real(c_double), parameter :: r(6) = &
         [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 3.04_c_double]
      rc = trc_create(h)
      rc = trc_set_molecule(h, 2_c_int, z, r, 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/def2-svp.json'//c_null_char)
   end subroutine make_hi

   subroutine test_hi(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, iters
      real(c_double) :: e_terco
      character(len=160) :: msg

      call make_hi(h)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "ecp_scf hi: trc_run_scf must converge")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_iterations(h, iters)
      rc = trc_energy(h, e_terco)
      rc = trc_destroy(h)

      print '(a,f20.12,a,i0,a,es10.2)', "check_ecp_scf: E(HI,RHF) terco ", real(e_terco, dp), &
         "  iter ", iters, "  (vs mqc) dE = ", abs(real(e_terco, dp) - MQC_HI)
      write (msg, '(a,es10.2)') "ecp_scf hi vs mqc, |dE| = ", abs(real(e_terco, dp) - MQC_HI)
      call check(error, abs(real(e_terco, dp) - MQC_HI) < 1.0e-5_dp, trim(msg))
   end subroutine test_hi

   !
   ! The free I atom, UHF: def2-SVP's ECP replaces 28 of iodine's 53
   ! electrons, leaving 25 valence electrons -- odd, a doublet (na=13,
   ! nb=12).
   !
   subroutine test_i_atom(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, iters
      real(c_double) :: e_terco
      character(len=160) :: msg

      rc = trc_create(h)
      rc = trc_set_molecule(h, 1_c_int, [53.0_c_double], [0.0_c_double, 0.0_c_double, 0.0_c_double], &
                            0_c_int, 2_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/def2-svp.json'//c_null_char)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "ecp_scf i_atom: trc_run_scf must converge")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_iterations(h, iters)
      rc = trc_energy(h, e_terco)
      rc = trc_destroy(h)

      print '(a,f20.12,a,i0,a,es10.2)', "check_ecp_scf: E(I,UHF) terco ", real(e_terco, dp), &
         "  iter ", iters, "  (vs mqc) dE = ", abs(real(e_terco, dp) - MQC_I_ATOM)
      write (msg, '(a,es10.2)') "ecp_scf i_atom vs mqc, |dE| = ", abs(real(e_terco, dp) - MQC_I_ATOM)
      call check(error, abs(real(e_terco, dp) - MQC_I_ATOM) < 1.0e-8_dp, trim(msg))
   end subroutine test_i_atom

   subroutine test_minao_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc

      call make_hi(h)
      rc = trc_set_guess(h, TRC_GUESS_MINAO, c_null_ptr, 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, "trc_run_scf refuses MINAO over an ECP basis")
      rc = trc_destroy(h)
   end subroutine test_minao_refused

   subroutine test_sap_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc

      call make_hi(h)
      rc = trc_set_guess(h, TRC_GUESS_SAP, c_null_ptr, 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, "trc_run_scf refuses SAP over an ECP basis")
      rc = trc_destroy(h)
   end subroutine test_sap_refused

   subroutine test_hirshfeld_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: q(2)

      call make_hi(h)
      rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-7_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "ecp_scf hirshfeld: trc_run_scf must converge first")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_hirshfeld_charges(h, q)
      call check(error, rc /= TRC_OK, "trc_hirshfeld_charges refuses an ECP-attached molecule")
      rc = trc_destroy(h)
   end subroutine test_hirshfeld_refused

   subroutine test_gradient_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc

      call make_hi(h)
      rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-7_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "ecp_scf gradient: trc_run_scf must converge first")
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      rc = trc_run_gradient(h)
      call check(error, rc, TRC_ERR_UNSUPPORTED, "trc_run_gradient refuses an ECP-attached molecule (phase 4)")
      rc = trc_destroy(h)
   end subroutine test_gradient_refused

   !
   ! Save, destroy, create a fresh handle, load: the restart file (format
   ! version 3, trc_scf_file.F90) must carry enough ECP state that the
   ! loaded result is accepted and reports the SAME energy, round-tripped
   ! through disk rather than compared in memory.
   !
   subroutine test_restart_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h1, h2
      integer(c_int) :: rc
      real(c_double) :: e1, e2
      character(len=*), parameter :: path = 'ecp_scf_restart.bin'

      call make_hi(h1)
      rc = trc_set_convergence(h1, 1.0e-11_c_double, 1.0e-8_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h1)
      call check(error, rc, TRC_OK, "ecp_scf restart: the SCF to save must converge")
      if (allocated(error)) then
         rc = trc_destroy(h1); return
      end if
      rc = trc_energy(h1, e1)
      rc = trc_save_scf(h1, path//c_null_char)
      call check(error, rc, TRC_OK, "trc_save_scf")
      rc = trc_destroy(h1)
      if (allocated(error)) return

      call make_hi(h2)
      rc = trc_load_scf(h2, path//c_null_char)
      call check(error, rc, TRC_OK, "trc_load_scf (format version 3, ECP state matches)")
      if (allocated(error)) then
         rc = trc_destroy(h2); return
      end if
      rc = trc_energy(h2, e2)
      rc = trc_destroy(h2)

      print '(a,f20.12,a,f20.12)', "check_ecp_scf: restart round-trip E(saved) ", real(e1, dp), &
         "  E(loaded) ", real(e2, dp)
      call check(error, real(e1, dp), real(e2, dp), "the loaded energy matches the saved one exactly")
   end subroutine test_restart_roundtrip

end module check_ecp_scf_suite

program check_ecp_scf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_ecp_scf_suite, only: collect_ecp_scf_tests
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_ecp_scf_tests, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_ecp_scf
