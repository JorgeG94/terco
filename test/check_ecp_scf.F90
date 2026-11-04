!
! ECP phase 3, end to end, through the C API for the refusal/restart
! checks and through the lower-level Fortran API (trc_basis_from_json,
! trc_ecp_from_json, trc_scf_run) for the four energy cases below, which
! need the converged DFT grid itself (trc_scf_result_t%grid_r/grid_w) to
! hand pyscf the identical quadrature -- the C API has no accessor for it.
!
! trc_set_molecule with the BARE atomic number (ECP_DESIGN.md's DECIDED
! charge contract), auto-attach from the orbital basis file's own
! ecp_potentials (Z - ncore, mirroring trc_capi.F90's attach_ecp), and the
! refusals phase 3 adds (MINAO/SAP/Hirshfeld guesses, the analytic
! gradient).
!
! GEOMETRY. HI's bond is mqc's own
! (metalquicha/validation/inputs/sample_inputs/hi.xyz): I at the origin, H
! at z = 1.609 Angstrom = 3.040569334522864 bohr (CODATA 2018, terco's own
! ANGSTROM_TO_BOHR). An earlier version of this file used the rounded
! literal 3.04 bohr (check_ecp_matrix.F90's own "hi_svp" case, which only
! checks the U_ECP matrix and so is insensitive to a 6e-4 bohr geometry
! error) and measured |dE| ~1e-6 to mqc for RHF/PBE/B3LYP -- a real but
! small geometry mismatch, not a method or integral problem.
!
! THE PYSCF ANCHOR (test/pyscf_refs/generate.py's gen_ecp_scf,
! ecp_scf.ref): cart=False (pyscf's own spherical AOs, matching this
! build's own default for a basis tagged gto_spherical -- def2-svp is),
! the SAME JSON basis+ECP (bse_json_basis/bse_json_ecp) and mqc's exact
! geometry above. For the two DFT cases, pyscf's grids.coords/weights are
! set to THIS run's own dumped grid (ecp_scf_probe.bin, this module's own
! run_all), the same "same-grid comparison" every other DFT anchor in this
! tree uses (gen_rks's own pattern), so a grid-quadrature difference is
! not a variable here either. Measured agreement against the anchor
! (2026-11-04, gfortran 15, LMAX=2): RHF HI 2.8e-13, UHF I atom 3.4e-13,
! RKS/PBE HI 7.4e-13, RKS/B3LYP HI 1.3e-12 -- all at the double-precision
! floor, run_all's own tolerance (1e-9) keeps 3+ orders of margin.
!
! AN EARLIER VERSION OF THIS HEADER CLAIMED A SECOND, GENUINE RHF
! STATIONARY POINT: that pyscf's own UNGUIDED RHF on this system, fed the
! identical JSON basis+ECP, converged to a different, LOWER energy than
! mqc/terco's (|dE| ~9.6e-5), and that this was confirmed (H_core
! eigenvalue identity, a warm-start experiment) rather than a mismatch.
! THAT WAS WRONG. The pyscf calls behind both checks used `cart=True`
! (33 AOs: Cartesian d has 6 components), not this build's own default
! `cart=False` (31 AOs: spherical d has 5) -- a genuinely LARGER,
! DIFFERENT basis, not a second root of the same one. At `cart=False`,
! the SAME geometry, pyscf's RHF reaches EXACTLY mqc/terco's own energy
! (measured |dE| 1.1e-12), stable under every init_guess tried (minao,
! atom, 1e, huckel) and pyscf's own `.stability()` check -- there is one
! RHF solution here, not two. The "H_core eigenvalues match" check from
! that version was itself done at `cart=True` on both sides, so it is
! still true (and still not interesting: it says nothing about which
! basis produces the lower variational energy, only that terco's and
! pyscf's own Cartesian integrals over the SAME larger basis agree); the
! "warm start converges to pyscf's energy" check also ran terco in
! `cart=True`, so it was a trivial same-basis, same-answer check, not
! evidence of two roots either. Measured 2026-11-04, see commit history
! for both the wrong claim and this correction.
!
! THE RESTART-TAKES-14-ITERATIONS QUESTION (a SEPARATE experiment: seeding
! terco's SPHERICAL run from pyscf's own converged SPHERICAL density via
! TRC_GUESS_GIVEN + TRC_AO_NORM_UNIT) is explained, not a bug, and not
! fixed here -- see this module's run_all for where that density would
! have to be rescaled and why TRC_AO_NORM_UNIT does not do it. Checked:
! terco's own spherical d-shell AOs are not unit-normalised (self-overlap
! diag(S) = 0.3979 for def2-SVP iodine's own d shell, against pyscf's 1 --
! generate.py's own terco_to_pyscf_perm_sph docstring already documents
! this exact number for a cc-pVDZ oxygen d shell and calls a uniform
! per-AO scale error in a SEEDED density harmless for convergence, which
! is exactly what is observed: every guess here, including that one,
! converges to the identical energy, just with more or fewer iterations).
! TRC_AO_NORM_UNIT (trc_set_ao_norm, trc_ao_norm.F90) is scoped to a
! CARTESIAN shell's own per-component convention (GAMESS's) by its own ABI
! doc comment; it was never meant to, and does not, rescale a spherical
! shell's cross-code normalisation at all, so passing it here changes
! nothing -- not a terco bug, and the seeding script's permutation itself
! is correct (checked: pyscf's and terco's s/p AOs here already agree
! bit-for-bit, and both d-shell orders also agree. Only the d-shell's
! OWN normalisation differs, a known, already-documented gap in
! TRC_GUESS_GIVEN for a spherical basis with l>=2, costing iterations,
! never correctness).
!
! RI-MP2, frozen core 4 (not a committed case here): measured E = mqc's
! -297.360082905691 to |dE| = 1.32e-5 at the corrected geometry, using
! def2-svp-rifit.json as the auxiliary basis -- an order of magnitude
! looser than RHF/DFT, consistent with a density-fitting/auxiliary-basis
! difference from whatever mqc's own RI-MP2 used (not confirmed to be the
! identical aux set), rather than a bug; left unadded rather than
! committing a loose, unexplained tolerance.
!
! An I2 case (two ECP atoms, exercising attach_ecp's per-atom loop) was
! tried at several I-I separations (4.0-10.0 bohr) and refused every time
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
! Writes ecp_scf_probe.bin for generate.py's gen_ecp_scf, same idea as
! scf_robust.F90's own probe/ref pair; the probe write and the sequential
! pyscf-reference read run once (run_all, below), for the same reason
! scf_robust's own suite stays parallel=.false. -- see that module's own
! header.
!
module check_ecp_scf_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char
   use pic_types, only: dp, int64
   use trc_constants, only: ANGSTROM_TO_BOHR
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_ecp_types, only: trc_ecp_t
   use trc_ecp_json, only: trc_ecp_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_real
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_convergence, trc_set_guess, trc_run_scf, trc_run_gradient, &
                               trc_hirshfeld_charges, trc_save_scf, trc_load_scf, &
                               trc_energy, trc_iterations, &
                               TRC_OK, TRC_ERR_UNSUPPORTED, TRC_GUESS_MINAO, TRC_GUESS_SAP
   implicit none
   private
   public :: collect_ecp_scf_tests

   real(dp), parameter :: HI_BOND_ANGSTROM = 1.609_dp
   !! mqc's own HI geometry (metalquicha/validation/inputs/sample_inputs/
   !! hi.xyz): I at the origin, H at z = 1.609 Angstrom.
   real(dp), parameter :: HI_BOND_BOHR = HI_BOND_ANGSTROM*ANGSTROM_TO_BOHR
   !! = 3.040569334522864 bohr (CODATA 2018, AGENTS.md).

   real(dp), parameter :: MQC_HI = -297.23153166336_dp
   !! RHF HI/def2-SVP at HI_BOND_BOHR -- a SECONDARY cross-check against
   !! the pyscf anchor (ecp_scf.ref), not the primary assertion.
   real(dp), parameter :: MQC_I_ATOM = -296.649606475643_dp
   !! The free I atom, UHF: def2-svp's 28-electron ECP core leaves 25
   !! valence electrons -- odd, a doublet. No geometry at all, so there is
   !! nothing for terco and mqc's own number to disagree about
   !! independently.

   integer, parameter :: NCASE = 4
   character(len=16), parameter :: CASE_NAME(NCASE) = [character(len=16) :: "hi", "i_atom", "hi_pbe", "hi_b3lyp"]
   !! generate.py's gen_ecp_scf reads this suite's own dump in this exact
   !! order and writes E_<name> back in the same order; keep both lists in
   !! step.

   logical :: computed = .false.
   logical :: ok_case(NCASE)
   character(len=256) :: msg_case(NCASE)
   real(dp) :: e_all(NCASE)

contains

   subroutine collect_ecp_scf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF HI/def2-SVP (auto-attached ECP) vs pyscf", case_hi), &
                  new_unittest("UHF I atom/def2-SVP (ECP, doublet) vs pyscf", case_i_atom), &
                  new_unittest("RKS/PBE HI/def2-SVP (ECP) vs pyscf (same grid)", case_hi_pbe), &
                  new_unittest("RKS/B3LYP HI/def2-SVP (ECP) vs pyscf (same grid)", case_hi_b3lyp), &
                  new_unittest("MINAO guess is refused with an ECP attached", test_minao_refused), &
                  new_unittest("SAP guess is refused with an ECP attached", test_sap_refused), &
                  new_unittest("Hirshfeld charges are refused with an ECP attached", test_hirshfeld_refused), &
                  new_unittest("the analytic gradient is refused with an ECP attached", test_gradient_refused), &
                  new_unittest("a restart file round-trips an ECP SCF (format version 3)", test_restart_roundtrip) &
                  ]
   end subroutine collect_ecp_scf_tests

   subroutine case_hi(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, ok_case(1), trim(msg_case(1)))
   end subroutine case_hi

   subroutine case_i_atom(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, ok_case(2), trim(msg_case(2)))
   end subroutine case_i_atom

   subroutine case_hi_pbe(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, ok_case(3), trim(msg_case(3)))
   end subroutine case_hi_pbe

   subroutine case_hi_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, ok_case(4), trim(msg_case(4)))
   end subroutine case_hi_b3lyp

   subroutine ensure_computed()
      if (computed) return
      call run_all()
      computed = .true.
   end subroutine ensure_computed

   !
   ! All four energy cases, run once: builds each molecule/basis with the
   ! ECP attached exactly as trc_capi.F90's attach_ecp does (bare Z minus
   ! at_ncore), runs trc_scf_run directly (not the C API) so the converged
   ! DFT grid (trc_scf_result_t%grid_r/grid_w) is available to dump, writes
   ! ecp_scf_probe.bin, then reads ecp_scf.ref back and compares.
   !
   subroutine run_all()
      integer :: ic, unit, natm, refu
      logical :: has_ref
      real(dp) :: e_ref, de, mqc_val
      logical :: has_mqc
      character(len=16) :: tag

      open (newunit=unit, file='ecp_scf_probe.bin', access='stream', form='unformatted', status='replace')
      write (unit) int(NCASE, int64)

      call ref_open('ecp_scf.ref', refu, has_ref)
      ok_case = .false.
      msg_case = "pyscf reference missing -- run test/pyscf_refs/generate.py ecp_scf"

      do ic = 1, NCASE
         natm = merge(1, 2, CASE_NAME(ic) == "i_atom")
         call run_one_case(ic, natm, unit)

         has_mqc = CASE_NAME(ic) == "hi" .or. CASE_NAME(ic) == "i_atom"
         mqc_val = merge(MQC_HI, MQC_I_ATOM, CASE_NAME(ic) == "hi")
         if (has_mqc) then
            print '(a,a,a,f20.12,a,es10.2)', "check_ecp_scf: E(", trim(CASE_NAME(ic)), ") terco ", &
               e_all(ic), "  (vs mqc) dE = ", abs(e_all(ic) - mqc_val)
         else
            print '(a,a,a,f20.12)', "check_ecp_scf: E(", trim(CASE_NAME(ic)), ") terco ", e_all(ic)
         end if

         if (.not. has_ref) cycle
         call ref_read_real(refu, 'E_'//trim(CASE_NAME(ic)), e_ref)
         de = abs(e_all(ic) - e_ref)
         write (tag, '(a)') trim(CASE_NAME(ic))
         write (msg_case(ic), '(a,a,a,es10.2)') "ecp_scf ", trim(tag), " vs pyscf, |dE| = ", de
         ! Measured (2026-11-04, this build, gfortran 15, LMAX=2): hi
         ! 2.8e-13, i_atom 3.4e-13, hi_pbe 7.4e-13, hi_b3lyp 1.3e-12 (same
         ! grid as pyscf for the two DFT cases) -- all at the double-
         ! precision floor. 1e-9 keeps 3+ orders of margin for a different
         ! compiler/platform's last-bit accumulation order.
         ok_case(ic) = de < 1.0e-9_dp
      end do
      if (has_ref) call ref_close(refu)
      close (unit)
   end subroutine run_all

   !
   ! One case's own build, SCF and probe dump; writes e_all(ic) for
   ! run_all's own two uses (the mqc print and the pyscf-ref compare).
   !
   subroutine run_one_case(ic, natm, unit)
      integer, intent(in) :: ic
      integer, intent(in) :: natm
      integer, intent(in) :: unit

      integer, allocatable :: elements(:)
      real(dp), allocatable :: at_r(:, :), zred(:)
      integer :: charge, mult, nelec, nunp, nalpha, nbeta, ia
      type(trc_basis_t) :: bas
      type(trc_ecp_t) :: ecp
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      logical :: is_dft
      integer(int64) :: npts

      allocate (elements(natm), at_r(3, natm))
      charge = 0
      if (natm == 2) then
         elements = [1, 53]
         at_r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
         at_r(:, 2) = [0.0_dp, 0.0_dp, HI_BOND_BOHR]
         mult = 1
      else
         elements = [53]
         at_r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
         mult = 2
      end if

      call trc_ecp_from_json('basis_sets/def2-svp.json', natm, elements, ecp, err)
      if (err%has_error()) error stop 'check_ecp_scf: '//err%get_message()
      allocate (zred(natm))
      do ia = 1, natm
         zred(ia) = real(elements(ia), dp) - real(ecp%at_ncore(ia), dp)
      end do
      call trc_basis_from_json('basis_sets/def2-svp.json', natm, elements, at_r, bas, err, &
                               overrides=trc_basis_json_overrides(charges=zred, allow_ecp=.true.))
      if (err%has_error()) error stop 'check_ecp_scf: '//err%get_message()
      call bas%to_device()

      nelec = nint(sum(zred)) - charge
      nunp = mult - 1
      nalpha = (nelec + nunp)/2
      nbeta = nelec - nalpha

      is_dft = CASE_NAME(ic) == "hi_pbe" .or. CASE_NAME(ic) == "hi_b3lyp"
      opts%functional = ""
      if (CASE_NAME(ic) == "hi_pbe") opts%functional = "pbe"
      if (CASE_NAME(ic) == "hi_b3lyp") opts%functional = "b3lyp"
      opts%grid_level = 3
      opts%max_iteration = 150
      opts%conv_energy = 1.0e-11_dp
      opts%conv_diis = 1.0e-8_dp

      call trc_scf_run(bas, nalpha, nbeta, opts, res, ecp=ecp)
      if (.not. res%converged) error stop 'check_ecp_scf: '//trim(CASE_NAME(ic))//' did not converge: '// &
                                          trim(res%message)
      e_all(ic) = res%energy

      write (unit) CASE_NAME(ic)
      write (unit) int(natm, int64), int(elements, int64)
      write (unit) at_r
      write (unit) int(charge, int64), int(mult, int64), int(nalpha, int64), int(nbeta, int64)
      npts = 0
      if (is_dft) npts = int(size(res%grid_w), int64)
      write (unit) npts
      if (npts > 0) then
         write (unit) res%grid_r
         write (unit) res%grid_w
      end if

      call bas%release()
      call ecp%release()
   end subroutine run_one_case

   !
   ! HI through the C API (make_hi), for the refusal/restart checks below
   ! -- none of those need the grid, so the plain context-handle path is
   ! simpler than run_one_case's own lower-level one.
   !
   subroutine make_hi(h)
      type(c_ptr), intent(out) :: h
      integer(c_int) :: rc
      real(c_double), parameter :: z(2) = [1.0_c_double, 53.0_c_double]
      real(c_double), parameter :: r(6) = &
         [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, real(HI_BOND_BOHR, c_double)]
      rc = trc_create(h)
      rc = trc_set_molecule(h, 2_c_int, z, r, 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/def2-svp.json'//c_null_char)
   end subroutine make_hi

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
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_bind_device
   use check_ecp_scf_suite, only: collect_ecp_scf_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, dev

   dev = trc_bind_device(0)
   stat = 0
   testsuites = [new_testsuite("check_ecp_scf", collect_ecp_scf_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_ecp_scf
