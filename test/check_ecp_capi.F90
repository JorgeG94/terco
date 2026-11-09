!
! The C API's ECP state, beyond check_ecp_scf's energies: which ECP attaches
! and when, what is refused, and that a refusal leaves the context as it was.
!
!   - trc_set_ecp_json with the basis's own file gives the auto-attached HI
!     energy (the pyscf anchor E_hi, ecp_scf.ref); a different iodine ECP
!     with the same 28-electron core (basis_sets/ecp-i-modified.json, a test
!     fixture) gives a different one, so the explicit file is what attaches.
!   - RbI in def2-SVP with that iodine-only ECP file is refused: rubidium's
!     def2-SVP shells are valence-only and nothing replaces its core.
!   - trc_set_molecule after an ECP basis drops the basis (it was built for
!     the reduced charges); reloading it gives E_hi again.
!   - a basis load refused after the ECP was read (an orbital/auxiliary
!     convention mismatch; a core larger than an atom's charge) leaves the
!     old basis, charges and ECP: the SCF that follows matches a context
!     that never tried.
!   - an SCF whose basis was built for other charges than the molecule now
!     holds is refused (TRC_ERR_STATE).
!   - an orbital basis from arrays or libcint is refused while an ECP is set
!     or attached; trc_context_basis is refused while one is attached.
!   - RI-MP2's automatic frozen core is refused (TRC_ERR_UNSUPPORTED).
!   - the restart file: a file from another ECP with the same core count, a
!     different core count, and a version 2 or 3 file are refused by an ECP
!     run; a version 3 file still loads into a run with no ECP.
!
! A test-drive suite.
!
module check_ecp_capi_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char, c_null_ptr
   use pic_types, only: dp, int32
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use trc_constants, only: ANGSTROM_TO_BOHR
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_real
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_elements, trc_set_basis_json, &
                               trc_set_ecp_json, trc_set_aux_json, trc_set_force_cartesian, &
                               trc_set_basis_arrays, trc_set_basis_libcint, trc_context_basis, &
                               trc_set_convergence, trc_run_scf, trc_run_rimp2, trc_energy, trc_message, &
                               trc_save_scf, trc_load_scf, trc_set_guess, TRC_GUESS_CORE, &
                               TRC_OK, TRC_ERR_BADARG, TRC_ERR_STATE, TRC_ERR_UNSUPPORTED, TRC_ERR_NOCONV
   implicit none
   private
   public :: collect_ecp_capi_tests

   character(len=*), parameter :: DEF2 = 'basis_sets/def2-svp.json'
   character(len=*), parameter :: ECP_MOD = 'basis_sets/ecp-i-modified.json'
   !! def2-SVP's iodine ECP with one s-channel coefficient changed, and no
   !! other element (its own "description" says so).
   character(len=*), parameter :: RB_STO3G = 'basis_sets/sto-3g-rb-zr.json'
   character(len=*), parameter :: RIFIT = 'basis_sets/def2-svp-rifit.json'
   real(dp), parameter :: HI_BOND_BOHR = 1.609_dp*ANGSTROM_TO_BOHR
   !! check_ecp_scf's HI: H at the origin, I on z.
   real(dp), parameter :: RB2_BOND_BOHR = 8.0_dp
   real(dp), parameter :: ANCHOR_TOL = 1.0e-9_dp
   !! check_ecp_scf's own tolerance against the same anchor.
   real(dp), parameter :: SAME_TOL = 1.0e-10_dp
   !! Two contexts doing the same thing on the host give the same bits; this
   !! allows for a GPU build's atomics.

contains

   subroutine collect_ecp_capi_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("trc_set_ecp_json with the basis's own file gives the anchor", test_explicit_same), &
                  new_unittest("trc_set_ecp_json with another iodine ECP attaches that one", test_explicit_other), &
                  new_unittest("RbI with an iodine-only ECP file is refused (Rb is valence-only)", test_partial_ecp), &
                  new_unittest("trc_set_molecule after an ECP basis drops the basis", test_set_molecule_drops), &
                  new_unittest("a refused basis load leaves the old basis and charges (mismatch)", test_refused_load), &
                  new_unittest("a refused core count leaves the old charges", test_refused_core), &
                  new_unittest("an SCF on a basis built for other charges is refused", test_charge_guard), &
                  new_unittest("trc_set_basis_arrays is refused with an ECP set or attached", test_arrays_refused), &
                  new_unittest("trc_set_basis_libcint is refused with an ECP set", test_libcint_refused), &
                  new_unittest("trc_context_basis is refused with an ECP attached", test_context_basis_refused), &
                  new_unittest("RI-MP2's automatic frozen core is refused with an ECP", test_rimp2_frozen_refused), &
                  new_unittest("restart: a file from another ECP with the same core is refused", test_restart_other_ecp), &
                  new_unittest("restart: a file with another core count is refused", test_restart_core_count), &
                  new_unittest("restart: a version 2 file is refused by an ECP run", test_restart_v2), &
                  new_unittest("restart: a version 3 file is refused by an ECP run", test_restart_v3), &
                  new_unittest("restart: a version 3 file loads into a run with no ECP", test_restart_v3_no_ecp) &
                  ]
   end subroutine collect_ecp_capi_tests

   ! ---- molecules -------------------------------------------------------

   subroutine new_hi(h)
      type(c_ptr), intent(out) :: h
      integer(c_int) :: rc
      rc = trc_create(h)
      rc = trc_set_molecule(h, 2_c_int, [1.0_c_double, 53.0_c_double], &
                            [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, &
                             real(HI_BOND_BOHR, c_double)], 0_c_int, 1_c_int)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 150_c_int, 10_c_int)
   end subroutine new_hi

   subroutine set_hi(h, rc)
      !! trc_set_molecule for HI again, on an existing handle.
      type(c_ptr), intent(in) :: h
      integer(c_int), intent(out) :: rc
      rc = trc_set_molecule(h, 2_c_int, [1.0_c_double, 53.0_c_double], &
                            [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, &
                             real(HI_BOND_BOHR, c_double)], 0_c_int, 1_c_int)
   end subroutine set_hi

   subroutine new_rb2(h, z2)
      !! Rb2 with atom 2's charge `z2` (its element stays Rb), a few SCF
      !! iterations only: these cases compare two contexts doing the same
      !! thing, which agree whether or not the SCF has converged.
      type(c_ptr), intent(out) :: h
      real(c_double), intent(in) :: z2
      integer(c_int) :: rc
      rc = trc_create(h)
      rc = trc_set_molecule(h, 2_c_int, [37.0_c_double, z2], &
                            [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, &
                             real(RB2_BOND_BOHR, c_double)], 0_c_int, 1_c_int)
      rc = trc_set_elements(h, 2_c_int, [37_c_int, 37_c_int])
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 4_c_int, 10_c_int)
      ! The core guess: no atomic SCF to fail on atom 2's odd charge.
      rc = trc_set_guess(h, TRC_GUESS_CORE, c_null_ptr, 1_c_int)
   end subroutine new_rb2

   ! ---- trc_set_ecp_json ----------------------------------------------

   subroutine test_explicit_same(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(dp) :: e

      call new_hi(h)
      rc = trc_set_ecp_json(h, DEF2//c_null_char)
      call check(error, rc, TRC_OK, "trc_set_ecp_json before the basis")
      if (allocated(error)) return
      call run_energy(h, DEF2, e, error)
      rc = trc_destroy(h)
      if (allocated(error)) return
      call check(error, abs(e - anchor_e_hi()) < ANCHOR_TOL, &
                 "explicit def2-SVP ECP: |E - E_hi(pyscf)| = "//rtoa(abs(e - anchor_e_hi())))
   end subroutine test_explicit_same

   subroutine test_explicit_other(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(dp) :: e

      call new_hi(h)
      rc = trc_set_ecp_json(h, ECP_MOD//c_null_char)
      call run_energy(h, DEF2, e, error)
      rc = trc_destroy(h)
      if (allocated(error)) return
      ! The changed coefficient moves the energy by about 2e-2 Eh.
      call check(error, abs(e - anchor_e_hi()) > 1.0e-4_dp, &
                 "the explicit ECP should replace the basis's own: E = "//rtoa(e)// &
                 ", the def2-SVP anchor is "//rtoa(anchor_e_hi()))
   end subroutine test_explicit_other

   subroutine test_partial_ecp(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      character(len=:), allocatable :: msg

      rc = trc_create(h)
      rc = trc_set_molecule(h, 2_c_int, [37.0_c_double, 53.0_c_double], &
                            [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, &
                             6.0_c_double], 0_c_int, 1_c_int)
      rc = trc_set_ecp_json(h, ECP_MOD//c_null_char)
      rc = trc_set_basis_json(h, DEF2//c_null_char)
      msg = message(h)
      call check(error, rc, TRC_ERR_BADARG, "RbI, def2-SVP, iodine-only ECP: should be refused; got '"//msg//"'")
      if (.not. allocated(error)) then
         call check(error, index(msg, "atom 1 (element 37)") > 0 .and. index(msg, "valence-only") > 0, &
                    "the refusal should name rubidium's valence-only basis: "//msg)
      end if
      if (.not. allocated(error)) then
         ! Nothing was loaded.
         rc = trc_run_scf(h)
         call check(error, rc, TRC_ERR_STATE, "no basis after the refused load")
      end if
      rc = trc_destroy(h)
   end subroutine test_partial_ecp

   ! ---- context state ---------------------------------------------------

   subroutine test_set_molecule_drops(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(dp) :: e

      call new_hi(h)
      rc = trc_set_basis_json(h, DEF2//c_null_char)
      call set_hi(h, rc)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_STATE, "trc_run_scf after trc_set_molecule should need the basis again "// &
                 "(the ECP basis carries the reduced charges); message: "//message(h))
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      call run_energy(h, DEF2, e, error)
      rc = trc_destroy(h)
      if (allocated(error)) return
      call check(error, abs(e - anchor_e_hi()) < ANCHOR_TOL, &
                 "reloaded: |E - E_hi(pyscf)| = "//rtoa(abs(e - anchor_e_hi())))
   end subroutine test_set_molecule_drops

   subroutine test_refused_load(error)
      !! Rb2 in all-electron STO-3G with a Cartesian auxiliary basis; then
      !! def2-SVP (spherical) is refused for the convention mismatch, AFTER
      !! its ECP was read. The SCF must still be the STO-3G one.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h, href
      integer(c_int) :: rc
      real(dp) :: e, eref

      call new_rb2(href, 37.0_c_double)
      call rb2_sto3g_cartesian(href)
      call run_any(href, eref, error)
      rc = trc_destroy(href)
      if (allocated(error)) return

      call new_rb2(h, 37.0_c_double)
      call rb2_sto3g_cartesian(h)
      rc = trc_set_force_cartesian(h, 0_c_int)
      rc = trc_set_basis_json(h, DEF2//c_null_char)
      call check(error, rc /= TRC_OK .and. index(message(h), "auxiliary") > 0, &
                 "def2-SVP over a Cartesian auxiliary basis should be refused; message: "//message(h))
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      call run_any(h, e, error)
      rc = trc_destroy(h)
      if (allocated(error)) return
      call check(error, abs(e - eref) < SAME_TOL, "after the refused load the SCF is not the STO-3G one: E = "// &
                 rtoa(e)//", a context that never tried: "//rtoa(eref))
   end subroutine test_refused_load

   subroutine test_refused_core(error)
      !! Rb2 with atom 2's charge 27: def2-SVP's 28-electron core would leave
      !! it below one, so the load is refused -- after atom 1's own charge
      !! was reduced, in the loop that finds it. The SCF must still be the
      !! all-electron one.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h, href
      integer(c_int) :: rc
      real(dp) :: e, eref

      call new_rb2(href, 27.0_c_double)
      rc = trc_set_basis_json(href, RB_STO3G//c_null_char)
      call run_any(href, eref, error)
      rc = trc_destroy(href)
      if (allocated(error)) return

      call new_rb2(h, 27.0_c_double)
      rc = trc_set_basis_json(h, RB_STO3G//c_null_char)
      rc = trc_set_basis_json(h, DEF2//c_null_char)
      call check(error, rc, TRC_ERR_BADARG, "a core larger than atom 2's charge should be refused")
      if (.not. allocated(error)) then
         call check(error, index(message(h), "fewer than one nuclear charge") > 0, &
                    "unexpected message: "//message(h))
      end if
      if (allocated(error)) then
         rc = trc_destroy(h); return
      end if
      call run_any(h, e, error)
      rc = trc_destroy(h)
      if (allocated(error)) return
      call check(error, abs(e - eref) < SAME_TOL, "after the refused core the SCF is not the STO-3G one: E = "// &
                 rtoa(e)//", a context that never tried: "//rtoa(eref))
   end subroutine test_refused_core

   subroutine test_charge_guard(error)
      !! Water in STO-3G, then trc_set_molecule with oxygen's charge 9: the
      !! basis (kept, no ECP) was built for 8.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double), parameter :: r(9) = [0.0_c_double, 0.0_c_double, 0.0_c_double, &
                                           0.0_c_double, 1.43_c_double, 1.1_c_double, &
                                           0.0_c_double, -1.43_c_double, 1.1_c_double]
      rc = trc_create(h)
      rc = trc_set_molecule(h, 3_c_int, [8.0_c_double, 1.0_c_double, 1.0_c_double], r, 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/sto-3g.json'//c_null_char)
      rc = trc_set_molecule(h, 3_c_int, [9.0_c_double, 1.0_c_double, 1.0_c_double], r, 1_c_int, 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_ERR_STATE, "an SCF on a basis built for other charges should be refused; "// &
                 "message: "//message(h))
      if (.not. allocated(error)) then
         call check(error, index(message(h), "other nuclear charges") > 0, "unexpected message: "//message(h))
      end if
      rc = trc_destroy(h)
   end subroutine test_charge_guard

   subroutine test_arrays_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc

      ! An explicit ECP set, no basis yet.
      call new_hi(h)
      rc = trc_set_ecp_json(h, DEF2//c_null_char)
      rc = set_two_s_shells(h)
      call check(error, rc, TRC_ERR_STATE, "trc_set_basis_arrays with trc_set_ecp_json set; message: "//message(h))
      rc = trc_destroy(h)
      if (allocated(error)) return

      ! An ECP attached by an earlier JSON basis.
      call new_hi(h)
      rc = trc_set_basis_json(h, DEF2//c_null_char)
      rc = set_two_s_shells(h)
      call check(error, rc, TRC_ERR_STATE, "trc_set_basis_arrays after an ECP basis; message: "//message(h))
      if (.not. allocated(error)) then
         ! Still the JSON basis with its ECP: the anchor energy.
         block
            real(dp) :: e
            rc = trc_run_scf(h)
            call check(error, rc, TRC_OK, "the JSON basis should still run")
            if (.not. allocated(error)) then
               rc = trc_energy(h, e)
               call check(error, abs(e - anchor_e_hi()) < ANCHOR_TOL, "E after the refused arrays: "//rtoa(e))
            end if
         end block
      end if
      rc = trc_destroy(h)
   end subroutine test_arrays_refused

   subroutine test_libcint_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      integer(c_int) :: atm(12), bas(16)
      real(c_double) :: env(30)

      ! libcint's layout: env's first 20 slots are reserved, then the two
      ! centres, then one s exponent and coefficient per shell.
      env = 0.0_c_double
      env(21:26) = [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, &
                    real(HI_BOND_BOHR, c_double)]
      env(27:28) = [1.0_c_double, 1.0_c_double]
      atm = [1_c_int, 20_c_int, 1_c_int, 0_c_int, 0_c_int, 0_c_int, &
             53_c_int, 23_c_int, 1_c_int, 0_c_int, 0_c_int, 0_c_int]
      bas = [0_c_int, 0_c_int, 1_c_int, 1_c_int, 0_c_int, 26_c_int, 27_c_int, 0_c_int, &
             1_c_int, 0_c_int, 1_c_int, 1_c_int, 0_c_int, 26_c_int, 27_c_int, 0_c_int]
      call new_hi(h)
      rc = trc_set_ecp_json(h, DEF2//c_null_char)
      rc = trc_set_basis_libcint(h, atm, 2_c_int, bas, 2_c_int, env, 30_c_int, 1_c_int)
      call check(error, rc, TRC_ERR_STATE, "trc_set_basis_libcint with trc_set_ecp_json set; message: "//message(h))
      rc = trc_destroy(h)
   end subroutine test_libcint_refused

   subroutine test_context_basis_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h, basis
      integer(c_int) :: rc

      call new_hi(h)
      rc = trc_set_basis_json(h, DEF2//c_null_char)
      rc = trc_context_basis(h, basis)
      call check(error, rc, TRC_ERR_STATE, "trc_context_basis with an ECP attached; message: "//message(h))
      rc = trc_destroy(h)
   end subroutine test_context_basis_refused

   subroutine test_rimp2_frozen_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc

      call new_hi(h)
      rc = trc_set_basis_json(h, DEF2//c_null_char)
      rc = trc_set_aux_json(h, RIFIT//c_null_char)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "HI SCF")
      if (.not. allocated(error)) then
         rc = trc_run_rimp2(h)
         call check(error, rc, TRC_ERR_UNSUPPORTED, "RI-MP2 with the automatic frozen core over an ECP; "// &
                    "message: "//message(h))
      end if
      if (.not. allocated(error)) then
         call check(error, index(message(h), "frozen_core is automatic") > 0, "unexpected message: "//message(h))
      end if
      rc = trc_destroy(h)
   end subroutine test_rimp2_frozen_refused

   ! ---- restart files ---------------------------------------------------

   subroutine test_restart_other_ecp(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      character(len=*), parameter :: path = 'ecp_capi_restart_def2.bin'

      call save_hi(path, error)
      if (allocated(error)) return
      call new_hi(h)
      rc = trc_set_ecp_json(h, ECP_MOD//c_null_char)
      rc = trc_set_basis_json(h, DEF2//c_null_char)
      rc = trc_load_scf(h, path//c_null_char)
      call check(error, rc, TRC_ERR_BADARG, "a file from the def2-SVP ECP into a run with another iodine "// &
                 "ECP (same 28-electron core); message: "//message(h))
      if (.not. allocated(error)) then
         call check(error, index(message(h), "different effective core potential") > 0, &
                    "unexpected message: "//message(h))
      end if
      rc = trc_destroy(h)
   end subroutine test_restart_other_ecp

   subroutine test_restart_core_count(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: path = 'ecp_capi_restart_core.bin', bad = 'ecp_capi_restart_core46.bin'
      character(len=1), allocatable :: bytes(:)
      integer :: off

      call save_hi(path, error)
      if (allocated(error)) return
      call slurp_bytes(path, bytes)
      ! The core counts follow the molecule: iodine is atom 2.
      off = ncore_offset(2) + 4
      bytes(off + 1:off + 4) = transfer(46_int32, bytes(1:4))
      call spill_bytes(bad, bytes)
      call expect_hi_load_refused(bad, "core electron count per atom does not match", error)
   end subroutine test_restart_core_count

   subroutine test_restart_v2(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: path = 'ecp_capi_restart_v4.bin', old = 'ecp_capi_restart_v2.bin'
      character(len=1), allocatable :: bytes(:)

      call save_hi(path, error)
      if (allocated(error)) return
      call slurp_bytes(path, bytes)
      call spill_bytes(old, downgrade(bytes, 2, 2))
      call expect_hi_load_refused(old, "format version 2", error)
   end subroutine test_restart_v2

   subroutine test_restart_v3(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: path = 'ecp_capi_restart_v4b.bin', old = 'ecp_capi_restart_v3.bin'
      character(len=1), allocatable :: bytes(:)

      call save_hi(path, error)
      if (allocated(error)) return
      call slurp_bytes(path, bytes)
      call spill_bytes(old, downgrade(bytes, 2, 3))
      call expect_hi_load_refused(old, "format version 3", error)
   end subroutine test_restart_v3

   subroutine test_restart_v3_no_ecp(error)
      !! Water/STO-3G saved as version 4, rewritten as the version 3 file the
      !! same run would have written: it must still load, with the energy.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: path = 'ecp_capi_water_v4.bin', old = 'ecp_capi_water_v3.bin'
      real(c_double), parameter :: r(9) = [0.0_c_double, 0.0_c_double, 0.0_c_double, &
                                           0.0_c_double, 1.43_c_double, 1.1_c_double, &
                                           0.0_c_double, -1.43_c_double, 1.1_c_double]
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: e1, e2
      character(len=1), allocatable :: bytes(:)

      rc = trc_create(h)
      rc = trc_set_molecule(h, 3_c_int, [8.0_c_double, 1.0_c_double, 1.0_c_double], r, 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/sto-3g.json'//c_null_char)
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "water SCF")
      if (.not. allocated(error)) then
         rc = trc_energy(h, e1)
         rc = trc_save_scf(h, path//c_null_char)
         call check(error, rc, TRC_OK, "trc_save_scf")
      end if
      rc = trc_destroy(h)
      if (allocated(error)) return
      call slurp_bytes(path, bytes)
      call spill_bytes(old, downgrade(bytes, 3, 3))

      rc = trc_create(h)
      rc = trc_set_molecule(h, 3_c_int, [8.0_c_double, 1.0_c_double, 1.0_c_double], r, 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/sto-3g.json'//c_null_char)
      rc = trc_load_scf(h, old//c_null_char)
      call check(error, rc, TRC_OK, "a version 3 file into a run with no ECP; message: "//message(h))
      if (.not. allocated(error)) then
         rc = trc_energy(h, e2)
         call check(error, real(e1, dp), real(e2, dp), "the loaded energy")
      end if
      rc = trc_destroy(h)
   end subroutine test_restart_v3_no_ecp

   ! ---- helpers ---------------------------------------------------------

   function anchor_e_hi() result(e)
      !! E_hi, RHF HI/def2-SVP from pyscf (ecp_scf.ref's first value).
      real(dp) :: e
      integer :: u
      logical :: ok
      call ref_open('ecp_scf.ref', u, ok)
      if (.not. ok) error stop 'check_ecp_capi: ecp_scf.ref missing -- run test/pyscf_refs/generate.py ecp_scf'
      call ref_read_real(u, 'E_hi', e)
      call ref_close(u)
   end function anchor_e_hi

   subroutine run_energy(h, basis, e, error)
      !! Load `basis`, run the SCF to convergence, return the energy.
      type(c_ptr), intent(in) :: h
      character(len=*), intent(in) :: basis
      real(dp), intent(out) :: e
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc
      real(c_double) :: ec
      e = 0.0_dp
      rc = trc_set_basis_json(h, basis//c_null_char)
      call check(error, rc, TRC_OK, "trc_set_basis_json("//basis//"): "//message(h))
      if (allocated(error)) return
      rc = trc_run_scf(h)
      call check(error, rc, TRC_OK, "trc_run_scf: "//message(h))
      if (allocated(error)) return
      rc = trc_energy(h, ec)
      e = real(ec, dp)
   end subroutine run_energy

   subroutine run_any(h, e, error)
      !! trc_run_scf, converged or not (TRC_ERR_NOCONV): the energy it ends at.
      type(c_ptr), intent(in) :: h
      real(dp), intent(out) :: e
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc
      real(c_double) :: ec
      e = 0.0_dp
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK .or. rc == TRC_ERR_NOCONV, "trc_run_scf: "//message(h))
      if (allocated(error)) return
      rc = trc_energy(h, ec)
      e = real(ec, dp)
   end subroutine run_any

   subroutine rb2_sto3g_cartesian(h)
      !! All-electron STO-3G and the def2-SVP-RIFIT auxiliary basis, both
      !! Cartesian.
      type(c_ptr), intent(in) :: h
      integer(c_int) :: rc
      rc = trc_set_force_cartesian(h, 1_c_int)
      rc = trc_set_basis_json(h, RB_STO3G//c_null_char)
      if (rc /= TRC_OK) error stop 'check_ecp_capi: Rb2 STO-3G basis refused'
      rc = trc_set_aux_json(h, RIFIT//c_null_char)
      if (rc /= TRC_OK) error stop 'check_ecp_capi: Rb2 auxiliary basis refused'
   end subroutine rb2_sto3g_cartesian

   integer(c_int) function set_two_s_shells(h) result(rc)
      !! One s shell on each HI atom, through trc_set_basis_arrays.
      type(c_ptr), intent(in) :: h
      rc = trc_set_basis_arrays(h, 2_c_int, 1_c_int, [0_c_int, 0_c_int], [1_c_int, 1_c_int], &
                                [1.0_c_double, 1.0_c_double], [1.0_c_double, 1.0_c_double], &
                                [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, &
                                 real(HI_BOND_BOHR, c_double)])
   end function set_two_s_shells

   subroutine save_hi(path, error)
      !! HI/def2-SVP, auto-attached ECP, converged and saved to `path`.
      character(len=*), intent(in) :: path
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(dp) :: e
      call new_hi(h)
      call run_energy(h, DEF2, e, error)
      if (.not. allocated(error)) then
         rc = trc_save_scf(h, path//c_null_char)
         call check(error, rc, TRC_OK, "trc_save_scf: "//message(h))
      end if
      rc = trc_destroy(h)
   end subroutine save_hi

   subroutine expect_hi_load_refused(path, expect, error)
      !! `path` into HI/def2-SVP with its own ECP attached: refused, saying
      !! `expect`.
      character(len=*), intent(in) :: path, expect
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      call new_hi(h)
      rc = trc_set_basis_json(h, DEF2//c_null_char)
      rc = trc_load_scf(h, path//c_null_char)
      call check(error, rc, TRC_ERR_BADARG, path//" should be refused; message: "//message(h))
      if (.not. allocated(error)) then
         call check(error, index(message(h), expect) > 0, "the refusal should say '"//expect//"': "//message(h))
      end if
      rc = trc_destroy(h)
   end subroutine expect_hi_load_refused

   ! The restart file's layout (src/scf/trc_scf_file.F90's own header):
   ! "TERCOSCF", the version, 16 int32 sizes and flags, the functional (16
   ! characters) and the PCM model (8), 9 real64 scalars; then at_z, at_r and
   ! at_elem; then (version 3 on) at_ncore; then (version 4) nrow, maxnp and
   ! the ECP's rows.
   pure integer function header_bytes() result(n)
      n = 8 + 4 + 16*4 + 16 + 8 + 9*8
   end function header_bytes

   pure integer function ncore_offset(natm) result(n)
      !! Bytes before at_ncore.
      integer, intent(in) :: natm
      n = header_bytes() + natm*(8 + 24 + 4)
   end function ncore_offset

   function downgrade(bytes, natm, version) result(out)
      !! A version 4 file rewritten as the `version` (2 or 3) file the same
      !! run would have written: the ECP rows dropped, and for version 2
      !! at_ncore too.
      character(len=1), intent(in) :: bytes(:)
      integer, intent(in) :: natm, version
      character(len=1), allocatable :: out(:)
      integer :: p_ncore, p_ecp, p_rest
      integer(int32) :: nrow, maxnp
      p_ncore = ncore_offset(natm)
      p_ecp = p_ncore + 4*natm
      nrow = transfer(bytes(p_ecp + 1:p_ecp + 4), nrow)
      maxnp = transfer(bytes(p_ecp + 5:p_ecp + 8), maxnp)
      p_rest = p_ecp + 8
      if (nrow > 0) p_rest = p_rest + 16*nrow + 16*maxnp*nrow
      if (version == 2) then
         out = [bytes(1:p_ncore), bytes(p_rest + 1:)]
      else
         out = [bytes(1:p_ecp), bytes(p_rest + 1:)]
      end if
      out(9:12) = transfer(int(version, int32), out(1:4))
   end function downgrade

   subroutine slurp_bytes(path, bytes)
      character(len=*), intent(in) :: path
      character(len=1), allocatable, intent(out) :: bytes(:)
      integer :: u, n
      open (newunit=u, file=path, access='stream', form='unformatted', status='old', action='read')
      inquire (unit=u, size=n)
      allocate (bytes(n))
      read (u) bytes
      close (u)
   end subroutine slurp_bytes

   subroutine spill_bytes(path, bytes)
      character(len=*), intent(in) :: path
      character(len=1), intent(in) :: bytes(:)
      integer :: u
      open (newunit=u, file=path, access='stream', form='unformatted', status='replace', action='write')
      write (u) bytes
      close (u)
   end subroutine spill_bytes

   function message(h) result(text)
      type(c_ptr), intent(in) :: h
      character(len=:), allocatable :: text
      character(kind=c_char) :: buf(1024)
      integer(c_int) :: rc
      integer :: i
      buf = c_null_char
      rc = trc_message(h, buf, int(size(buf), c_int))
      text = ''
      do i = 1, size(buf)
         if (buf(i) == c_null_char) exit
         text = text//buf(i)
      end do
   end function message

   function rtoa(x) result(s)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: t
      write (t, '(es24.16)') x
      s = trim(adjustl(t))
   end function rtoa

end module check_ecp_capi_suite

program check_ecp_capi
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_api, only: trc_bind_device
   use check_ecp_capi_suite, only: collect_ecp_capi_tests
   implicit none
   integer :: stat, dev

   dev = trc_bind_device(0)
   stat = 0
   call run_testsuite(collect_ecp_capi_tests, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_ecp_capi
