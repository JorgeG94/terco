!
! ECP phase 3, end to end: RHF HI and UHF the free I atom through the C
! API (bare Z into trc_set_molecule, trc_set_basis_json's own auto-attach,
! trc_run_scf -- the path a caller, mqc included, actually uses); RKS/PBE
! and RKS/B3LYP on HI through the lower-level Fortran API instead
! (trc_basis_from_json/trc_ecp_from_json/trc_scf_run), because only that
! path hands back trc_scf_result_t%grid_r/grid_w, which generate.py's
! gen_ecp_scf needs to feed pyscf terco's own grid; the C API has no
! accessor for it. Also checks the refusals phase 3 adds (MINAO/SAP/
! Hirshfeld guesses, the analytic gradient) and a restart-file round trip.
!
! GEOMETRY: mqc's own HI bond
! (metalquicha/validation/inputs/sample_inputs/hi.xyz), I at the origin, H
! at z = 1.609 Angstrom = HI_BOND_BOHR (CODATA 2018, terco's own
! ANGSTROM_TO_BOHR).
!
! THE PYSCF ANCHOR (test/pyscf_refs/generate.py's gen_ecp_scf,
! ecp_scf.ref): cart=False (pyscf's own spherical AOs, matching this
! build's own default for a basis tagged gto_spherical -- def2-svp is),
! the SAME JSON basis+ECP (bse_json_basis/bse_json_ecp) and the geometry
! above. The two DFT cases feed pyscf THIS run's own dumped grid
! (ecp_scf_probe.bin, this module's own run_all) via grids.coords/weights,
! the same same-grid comparison every other DFT anchor in this tree uses
! (gen_rks's own pattern). Measured agreement against the anchor
! (2026-11-04, gfortran 15, LMAX=2): RHF HI 2.8e-13, UHF I atom 3.4e-13,
! RKS/PBE HI 7.4e-13, RKS/B3LYP HI 1.3e-12 -- all at the double-precision
! floor, run_all's own tolerance (1e-9) keeps 3+ orders of margin. mqc's
! own literals (MQC_HI, MQC_I_ATOM) are a secondary print only.
!
! Seeding terco's spherical HI run from pyscf's own converged spherical
! density (TRC_GUESS_GIVEN + TRC_AO_NORM_UNIT) takes MORE iterations than
! a cold SAD guess, not fewer: terco's own spherical l>=2 AOs are not
! unit-normalised (diag(S) = 0.398 for iodine's def2-SVP d shell, against
! pyscf's 1), and TRC_AO_NORM_UNIT does not rescale them -- it is scoped
! to a Cartesian shell's own per-component convention by its own ABI doc
! comment, not a spherical one. The given density's d-shell block is
! mis-scaled; the SCF still converges to the identical energy regardless,
! just with more iterations. Not fixed here.
!
! RI-MP2, frozen core 4 (not a committed case here): measured E = mqc's
! -297.360082905691 to |dE| = 1.32e-5 at the corrected geometry, using
! def2-svp-rifit.json as the auxiliary basis -- an order of magnitude
! looser than RHF/DFT, consistent with a density-fitting/auxiliary-basis
! difference from whatever mqc's own RI-MP2 used (not confirmed to be the
! identical aux set), rather than a bug; left unadded rather than
! committing a loose, unexplained tolerance.
!
! I2 (two ECP atoms, attach_ecp's per-atom loop and the semi-local blocks
! between two heavy atoms) is a case here now. It used to be refused at
! every separation tried: the ECP radial quadrature failed instead of
! keeping its finest level as pyscf does (#194 review, fixed in trc_ecp).
!
! H2 next to a GHOST iodine (trc_set_elements, charge 0) checks attach_ecp's
! ghost rule: the ghost keeps iodine's basis and gets no ECP and no core
! count, as in pyscf (an atom labelled GHOST-, with no ecp entry). Before,
! attach_ecp gave the ghost iodine's 28-electron core and refused it for a
! negative charge.
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
                               trc_hirshfeld_charges, trc_save_scf, trc_load_scf, trc_set_elements, &
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
   real(dp), parameter :: I2_BOND_BOHR = 5.04_dp
   !! I2 near its 2.666 Angstrom bond, as check_ecp_matrix's i2 case.

   real(dp), parameter :: MQC_HI = -297.23153166336_dp
   !! RHF HI/def2-SVP at HI_BOND_BOHR -- a SECONDARY cross-check against
   !! the pyscf anchor (ecp_scf.ref), not the primary assertion.
   real(dp), parameter :: MQC_I_ATOM = -296.649606475643_dp
   !! The free I atom, UHF: def2-svp's 28-electron ECP core leaves 25
   !! valence electrons -- odd, a doublet. No geometry at all, so there is
   !! nothing for terco and mqc's own number to disagree about
   !! independently.

   integer, parameter :: NCASE = 6
   character(len=16), parameter :: CASE_NAME(NCASE) = [character(len=16) :: "hi", "i_atom", "hi_pbe", "hi_b3lyp", &
                                                       "i2", "h2_ghost_i"]
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
                  new_unittest("RHF I2/def2-SVP (two ECP atoms) vs pyscf", case_i2), &
                  new_unittest("RHF H2 + ghost I/def2-SVP (the ghost gets no ECP) vs pyscf", case_h2_ghost_i), &
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

   subroutine case_i2(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, ok_case(5), trim(msg_case(5)))
   end subroutine case_i2

   subroutine case_h2_ghost_i(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, ok_case(6), trim(msg_case(6)))
   end subroutine case_h2_ghost_i

   subroutine ensure_computed()
      if (computed) return
      call run_all()
      computed = .true.
   end subroutine ensure_computed

   !
   ! All the energy cases, run once: builds each molecule/basis with the
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
         natm = case_natm(CASE_NAME(ic))
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
   ! RHF HI and UHF the free I atom go through the C API -- bare Z into
   ! trc_set_molecule, trc_set_basis_json's own auto-attach, trc_run_scf --
   ! the exact path a caller (mqc included) actually uses, so this is the
   ! one place that path's own energy is checked against anything external
   ! at all; everything else in this suite either only checks a refusal or
   ! only checks that an SCF converges. The two DFT cases stay on the
   ! lower-level path (trc_basis_from_json/trc_ecp_from_json/trc_scf_run)
   ! because only that path hands back trc_scf_result_t%grid_r/grid_w, and
   ! generate.py's gen_ecp_scf needs this run's own grid to feed pyscf the
   ! identical quadrature; the C API has no accessor for it.
   !
   subroutine run_one_case(ic, natm, unit)
      integer, intent(in) :: ic
      integer, intent(in) :: natm
      integer, intent(in) :: unit

      integer, allocatable :: elements(:)
      logical, allocatable :: ghost(:)
      real(dp), allocatable :: at_r(:, :), zred(:)
      real(dp), allocatable :: grid_r(:, :), grid_w(:)
      integer :: charge, mult, nelec, nunp, nalpha, nbeta, ia
      type(trc_ecp_t) :: ecp
      type(error_t) :: err
      logical :: is_dft
      integer(int64) :: npts

      allocate (elements(natm), at_r(3, natm), ghost(natm))
      charge = 0
      mult = 1
      ghost = .false.
      at_r = 0.0_dp
      select case (trim(CASE_NAME(ic)))
      case ("i_atom")
         elements = [53]
         mult = 2
      case ("i2")
         elements = [53, 53]
         at_r(3, 2) = I2_BOND_BOHR
      case ("h2_ghost_i")
         ! H2 at 1.4 bohr, a ghost iodine 3.6 bohr past the second H.
         elements = [1, 1, 53]
         ghost = [.false., .false., .true.]
         at_r(3, 2) = 1.4_dp
         at_r(3, 3) = 5.0_dp
      case default
         elements = [1, 53]
         at_r(3, 2) = HI_BOND_BOHR
      end select

      ! The per-atom reduced charge, needed either way for the probe's own
      ! na/nb (generate.py's pyscf mol needs the spin): trc_set_basis_json
      ! does this same reduction internally for the C-API path below, so
      ! this is not a second, independent implementation of attach_ecp,
      ! just the part of its output this dump also needs.
      call trc_ecp_from_json('basis_sets/def2-svp.json', natm, elements, ecp, err)
      if (err%has_error()) error stop 'check_ecp_scf: '//err%get_message()
      allocate (zred(natm))
      do ia = 1, natm
         zred(ia) = real(elements(ia), dp) - real(ecp%at_ncore(ia), dp)
      end do
      where (ghost) zred = 0.0_dp
      call ecp%release()
      nelec = nint(sum(zred)) - charge
      nunp = mult - 1
      nalpha = (nelec + nunp)/2
      nbeta = nelec - nalpha

      is_dft = CASE_NAME(ic) == "hi_pbe" .or. CASE_NAME(ic) == "hi_b3lyp"
      npts = 0
      if (is_dft) then
         call run_case_grid(ic, natm, elements, at_r, nalpha, nbeta, grid_r, grid_w)
         npts = int(size(grid_w), int64)
      else
         call run_case_capi(ic, natm, elements, ghost, at_r, charge, mult)
      end if

      write (unit) CASE_NAME(ic)
      write (unit) int(natm, int64), int(elements, int64)
      write (unit) int(merge(1, 0, ghost), int64)
      write (unit) at_r
      write (unit) int(charge, int64), int(mult, int64), int(nalpha, int64), int(nbeta, int64)
      write (unit) npts
      if (npts > 0) then
         write (unit) grid_r
         write (unit) grid_w
      end if
   end subroutine run_one_case

   !
   ! RHF HI / UHF I atom through the C API: bare Z, trc_set_basis_json's
   ! own auto-attach, trc_run_scf -- see run_one_case's own header.
   !
   subroutine run_case_capi(ic, natm, elements, ghost, at_r, charge, mult)
      integer, intent(in) :: ic, natm, elements(natm), charge, mult
      logical, intent(in) :: ghost(natm)
      real(dp), intent(in) :: at_r(3, natm)
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: z(natm), xyz(3*natm), e

      ! A ghost is charge 0 in trc_set_molecule and keeps its element
      ! through trc_set_elements, before the basis -- that call's own rule.
      z = merge(0.0_c_double, real(elements, c_double), ghost)
      xyz = reshape(real(at_r, c_double), [3*natm])
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), z, xyz, int(charge, c_int), int(mult, c_int))
      if (any(ghost)) then
         rc = trc_set_elements(h, int(natm, c_int), int(elements, c_int))
         if (rc /= TRC_OK) error stop 'check_ecp_scf: trc_set_elements refused'
      end if
      rc = trc_set_basis_json(h, 'basis_sets/def2-svp.json'//c_null_char)
      if (rc /= TRC_OK) error stop "check_ecp_scf: "//trim(CASE_NAME(ic))//": trc_set_basis_json refused"
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h)
      if (rc /= TRC_OK) error stop 'check_ecp_scf: '//trim(CASE_NAME(ic))//' (C API) did not converge'
      rc = trc_energy(h, e)
      rc = trc_destroy(h)
      e_all(ic) = real(e, dp)
   end subroutine run_case_capi

   !
   ! RKS/PBE and RKS/B3LYP on HI through the lower-level path, for the
   ! grid -- see run_one_case's own header.
   !
   subroutine run_case_grid(ic, natm, elements, at_r, nalpha, nbeta, grid_r, grid_w)
      integer, intent(in) :: ic, natm, elements(natm), nalpha, nbeta
      real(dp), intent(in) :: at_r(3, natm)
      real(dp), allocatable, intent(out) :: grid_r(:, :), grid_w(:)
      type(trc_basis_t) :: bas
      type(trc_ecp_t) :: ecp
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      real(dp), allocatable :: zred(:)
      integer :: ia

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

      grid_r = res%grid_r
      grid_w = res%grid_w

      call bas%release()
      call ecp%release()
   end subroutine run_case_grid

   pure integer function case_natm(name) result(n)
      character(len=*), intent(in) :: name
      select case (trim(name))
      case ("i_atom"); n = 1
      case ("h2_ghost_i"); n = 3
      case default; n = 2
      end select
   end function case_natm

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
