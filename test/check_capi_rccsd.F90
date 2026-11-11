!
! Closed-shell RI-CCSD and RI-CCSD(T) through the C ABI (trc_set_rccsd, trc_set_rccsd_convergence,
! trc_set_rccsd_algorithm, trc_run_rccsd, trc_rccsd_energy, and trc_energy after them): the numbers, the setters'
! arguments, the refusals and the invalidation of a stored result.
!
! Reference. Water, cc-pVDZ with cc-pVDZ-RIFIT, frozen core 1, the geometry of test/hirshfeld_water.xyz: pyscf's
! DF-CCSD(T) in test/pyscf_refs/rccsd_ints_w1.ref (E_MP2, E_CCSD_CORR, E_T, E_CCSD_T, T1_DIAG), converged there to
! 1e-12 / 1e-10. The handle converges the SCF to 1e-13 Eh and a commutator of 1e-11 and the coupled-cluster
! iteration to 1e-12 / 1e-10, as check_rccsd_run and check_rccsd_triples do, with their tolerances:
!   E_CCSD correlation   1e-10        E(T)   3e-12        CCSD(T) total energy   3e-11
! (trc_energy is the CCSD(T) total, so the third also covers the SCF energy and the sum).
!
! What else is checked, each against something that can fail:
!   - nfrozen = -1 (the core counted by element) gives the energies of nfrozen = 1, and a count of 0 does not;
!   - without the triples flag E(T) is not reported and trc_energy is the CCSD total;
!   - every setter refuses its bad arguments with TRC_ERR_BADARG and a message, and a refusal leaves the earlier
!     settings alone (the run after it still reproduces pyscf);
!   - each request the method does not cover is refused with its status and message: no SCF, an SCF that did not
!     converge, no auxiliary basis, an unrestricted reference, a Kohn-Sham one, an effective core potential,
!     RI-MP2 already run, every occupied orbital frozen, an auxiliary basis above this build's ceiling. (More than
!     one rank is refused in trc_run_rccsd as well but needs an MPI build and a communicator, so it is not run here.)
!   - a continuum solvent is NOT refused (as for RI-MP2), the energy is finite;
!   - a second trc_run_scf, a new auxiliary basis and every RI-CCSD setter drop the stored result;
!   - max_iterations = 2 returns TRC_ERR_NOCONV with a message, the last energies stay readable with converged = 0,
!     (T) is not run, and trc_energy does not include the correlation;
!   - after a converged run trc_run_gradient refuses and says why;
!   - the (T) step refusing after a converged RI-CCSD: the tiled routine with a budget below one triple. The
!     call returns TRC_ERR_UNSUPPORTED with a message naming the budget, trc_rccsd_energy still gives the converged
!     RI-CCSD (has_triples = 0, e_t = 0) and trc_energy is its total; a run with enough budget, and one with the
!     reference routine (which has no budget), give (T). The two routines agree on E(T) to 1e-12. The routine is
!     always set explicitly, so the case means the same in a host build and in a device build.
!
! A build whose density-fitting ceilings are below the cc-pVDZ / cc-pVDZ-RIFIT shells prints SKIP.
!
module check_capi_rccsd_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char, c_null_ptr
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use rccsd_ref_energies, only: rccsd_refs_t, read_rccsd_refs
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_convergence, trc_set_rimp2, trc_set_rccsd, trc_set_lindep, &
                               trc_rccsd_dimensions, trc_set_rccsd_convergence, &
                               trc_set_rccsd_algorithm, trc_set_method, trc_set_pcm, trc_set_unrestricted, &
                               trc_run_scf, trc_run_rimp2, trc_run_rccsd, trc_rccsd_energy, trc_run_gradient, &
                               trc_energy, trc_message, TRC_OK, TRC_ERR_NULL, TRC_ERR_BADARG, TRC_ERR_UNSUPPORTED, &
                               TRC_ERR_NOCONV, TRC_ERR_STATE
   implicit none
   private
   public :: collect_capi_rccsd

   real(dp), parameter :: TOL_MP2 = 1.0e-10_dp
   !! E_MP2 against pyscf, check_rccsd_ints' TOL_MP2
   real(dp), parameter :: TOL_CORR = 1.0e-10_dp
   !! E_CCSD correlation against pyscf, check_rccsd_run's
   real(dp), parameter :: TOL_T = 3.0e-12_dp
   !! E(T) against pyscf, check_rccsd_triples' TOL_PYSCF_CONV
   real(dp), parameter :: TOL_TOTAL = 3.0e-11_dp
   !! the CCSD(T) total energy against pyscf, check_rccsd_triples' TOL_PYSCF_TOTAL
   real(dp), parameter :: TOL_T1 = 1.0e-9_dp
   !! the T1 diagnostic, first order in the amplitude error (1e-10 here)
   real(dp), parameter :: TOL_ROUTINES = 1.0e-12_dp
   !! the reference and the tiled (T) routines, and a tiled run of another budget, on the same amplitudes
   real(dp), parameter :: TOL_SUM = 1.0e-13_dp
   !! trc_energy against the SCF energy plus the reported parts: the same numbers added in another order

   real(c_double), parameter :: HI_BOND_BOHR = 1.609_c_double*1.8897261254578281_c_double
   !! HI along z, the bond of test/hi.xyz
   character(len=*), parameter :: WATER = "hirshfeld_water.xyz"
   character(len=*), parameter :: REF = "rccsd_ints_w1.ref"

contains

   subroutine collect_capi_rccsd(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water CCSD(T) matches pyscf: MP2, CCSD, (T), totals, T1", test_energies), &
                  new_unittest("nfrozen -1 counts the core by element", test_frozen_by_element), &
                  new_unittest("triples off: no (T), trc_energy is the CCSD total", test_no_triples), &
                  new_unittest("the setters refuse bad arguments and keep the old settings", test_setters), &
                  new_unittest("no SCF, an unconverged SCF and no auxiliary basis are refused", test_state_refusals), &
                  new_unittest("unrestricted, open-shell, Kohn-Sham and ECP references are refused", test_reference_refusals), &
                  new_unittest("every occupied orbital frozen, and an auxiliary basis above the ceiling, are refused", &
                               test_size_refusals), &
                  new_unittest("a continuum solvent is not refused", test_pcm), &
                  new_unittest("a second SCF, a new auxiliary basis and a setter drop the result", test_invalidation), &
                  new_unittest("each convergence threshold and the DIIS size reach the solver", test_convergence_settings), &
                  new_unittest("max_iterations 2: NOCONV, the last energies, no (T)", test_noconv), &
                  new_unittest("a gradient after a converged RI-CCSD is refused, after an unconverged one it is not", &
                               test_gradient_refused), &
                  new_unittest("RI-MP2 and RI-CCSD: the later run wins, a refusal changes nothing", test_methods_alternate), &
                  new_unittest("(T) refused after a converged RI-CCSD: the CCSD result stands, the call says so", &
                               test_triples_refused), &
                  new_unittest("the reference and the tiled (T) routines agree; the routine and the budget are validated", &
                               test_triples_routines), &
                  new_unittest("the dimensions report the frozen core, the spaces and the auxiliary functions", &
                               test_dimensions), &
                  new_unittest("dropped linearly dependent functions are refused and the message names the option", &
                               test_lindep) &
                  ]
   end subroutine collect_capi_rccsd

   function supported() result(ok)
      !! cc-pVDZ is d and cc-pVDZ-RIFIT has f on O: both within the build's density-fitting ceilings
      logical :: ok
      ok = DF_LMAX >= 2 .and. DF_LMAX_AUX >= 3
   end function supported

   function message(h) result(text)
      !! The handle's last message
      type(c_ptr), intent(in) :: h
      character(len=:), allocatable :: text
      character(kind=c_char) :: buffer(512)
      integer(c_int) :: rc
      integer :: i
      rc = trc_message(h, buffer, 512_c_int)
      text = ""
      do i = 1, size(buffer)
         if (buffer(i) == c_null_char) exit
         text = text//buffer(i)
      end do
   end function message

   function nan() result(x)
      real(c_double) :: x
      x = ieee_value(1.0_c_double, ieee_quiet_nan)
   end function nan

   function infinity() result(x)
      real(c_double) :: x
      x = ieee_value(1.0_c_double, ieee_positive_inf)
   end function infinity

   subroutine make_handle(h, xyz, basis, aux, charge, multiplicity, rc)
      !! A handle for `xyz` in `basis`, with `aux` when it is not empty and the SCF tightly converged
      type(c_ptr), intent(out) :: h
      character(len=*), intent(in) :: xyz, basis, aux
      integer, intent(in) :: charge, multiplicity
      integer(c_int), intent(out) :: rc
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)

      call read_xyz(xyz, natm, zint, at_r)
      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                                              int(charge, c_int), int(multiplicity, c_int))
      if (rc == TRC_OK) rc = trc_set_basis_json(h, 'basis_sets/'//basis//'.json'//c_null_char)
      if (rc == TRC_OK .and. len(aux) > 0) then
         ! a name under test/basis_sets, or a path from test/ (the full-size sets live under ../basis_sets)
         if (aux(1:3) == "../") then
            rc = trc_set_aux_json(h, aux//'.json'//c_null_char)
         else
            rc = trc_set_aux_json(h, 'basis_sets/'//aux//'.json'//c_null_char)
         end if
      end if
      if (rc == TRC_OK) rc = trc_set_convergence(h, 1.0e-13_c_double, 1.0e-11_c_double, 100_c_int, 10_c_int)
   end subroutine make_handle

   subroutine water_scf(error, h, triples, nfrozen, ok)
      !! Water, cc-pVDZ / cc-pVDZ-RIFIT, RI-CCSD set with the tight convergence of the references, SCF run
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(out) :: h
      integer, intent(in) :: triples, nfrozen
      logical, intent(out) :: ok
      integer(c_int) :: rc

      ok = supported()
      if (.not. ok) then
         print "(a)", "SKIP: this build's density-fitting ceilings are below cc-pVDZ-RIFIT"
         return
      end if
      call make_handle(h, WATER, "cc-pvdz", "cc-pvdz-rifit", 0, 1, rc)
      if (rc == TRC_OK) rc = trc_set_rccsd(h, int(triples, c_int), int(nfrozen, c_int), 0_c_int)
      if (rc == TRC_OK) rc = trc_set_rccsd_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 100_c_int, 0_c_int)
      call check(error, rc == TRC_OK, "setting up the water handle failed")
      if (allocated(error)) return
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the water SCF failed")
   end subroutine water_scf

   subroutine test_energies(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      type(rccsd_refs_t) :: refs
      real(c_double) :: e_mp2, e_ccsd, e_t, t1, e_total
      integer(c_int) :: rc, iterations, converged, has_triples
      logical :: ok

      call water_scf(error, h, 1, 1, ok)
      if (.not. ok .or. allocated(error)) return
      call read_rccsd_refs(REF, refs)
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "trc_run_rccsd returned "//to_string(int(rc))//": "//message(h))
      if (allocated(error)) return
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, rc == TRC_OK .and. converged == 1 .and. has_triples == 1, "trc_rccsd_energy: not a converged (T) result")
      if (allocated(error)) return
      call check(error, abs(e_mp2 - refs%e_mp2) < TOL_MP2, "E_MP2 "//to_string(e_mp2, "f20.12")//" vs pyscf "// &
                 to_string(refs%e_mp2, "f20.12")//", |difference| "//to_string(abs(e_mp2 - refs%e_mp2), "es9.2"))
      if (allocated(error)) return
      call check(error, abs(e_ccsd - refs%e_ccsd_corr) < TOL_CORR, "E_CCSD "//to_string(e_ccsd, "f20.12")//" vs pyscf "// &
                 to_string(refs%e_ccsd_corr, "f20.12")//", |difference| "//to_string(abs(e_ccsd - refs%e_ccsd_corr), "es9.2"))
      if (allocated(error)) return
      call check(error, abs(e_t - refs%e_t) < TOL_T, "E(T) "//to_string(e_t, "f20.14")//" vs pyscf "// &
                 to_string(refs%e_t, "f20.14")//", |difference| "//to_string(abs(e_t - refs%e_t), "es9.2"))
      if (allocated(error)) return
      rc = trc_energy(h, e_total)
      call check(error, rc == TRC_OK .and. abs(e_total - refs%e_ccsd_t) < TOL_TOTAL, &
                 "trc_energy (the CCSD(T) total) "//to_string(e_total, "f20.12")//" vs pyscf "// &
                 to_string(refs%e_ccsd_t, "f20.12")//", |difference| "//to_string(abs(e_total - refs%e_ccsd_t), "es9.2"))
      if (allocated(error)) return
      call check(error, abs(t1 - refs%t1_diag) < TOL_T1, "T1 diagnostic "//to_string(t1, "f12.8")//" vs pyscf "// &
                 to_string(refs%t1_diag, "f12.8"))
      if (allocated(error)) return
      call check(error, iterations >= 5 .and. iterations < 100, "iterations "//to_string(int(iterations)))
      if (allocated(error)) return
      rc = trc_destroy(h)
      print "(a)", "  measured: |E_MP2| "//to_string(abs(e_mp2 - refs%e_mp2), "es9.2")//", |E_CCSD| "// &
         to_string(abs(e_ccsd - refs%e_ccsd_corr), "es9.2")//", |E(T)| "//to_string(abs(e_t - refs%e_t), "es9.2")// &
         ", |total| "//to_string(abs(e_total - refs%e_ccsd_t), "es9.2")//", |T1| "//to_string(abs(t1 - refs%t1_diag), "es9.2")
   end subroutine test_energies

   subroutine test_frozen_by_element(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(c_double) :: e_mp2, e_ccsd, e_t, t1, e_ccsd_by_count, e_ccsd_none
      integer(c_int) :: rc, iterations, converged, has_triples
      logical :: ok

      call water_scf(error, h, 0, 1, ok)
      if (.not. ok .or. allocated(error)) return
      rc = trc_run_rccsd(h)
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd_by_count, e_t, t1, iterations, converged, has_triples)
      rc = trc_set_rccsd(h, 0_c_int, -1_c_int, 0_c_int)
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "nfrozen -1: "//message(h))
      if (allocated(error)) return
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, abs(e_ccsd - e_ccsd_by_count) < 1.0e-14_dp, "nfrozen -1 "//to_string(e_ccsd, "f20.14")// &
                 " is not nfrozen 1 "//to_string(e_ccsd_by_count, "f20.14"))
      if (allocated(error)) return
      rc = trc_set_rccsd(h, 0_c_int, 0_c_int, 0_c_int)
      rc = trc_run_rccsd(h)
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd_none, e_t, t1, iterations, converged, has_triples)
      call check(error, abs(e_ccsd_none - e_ccsd_by_count) > 1.0e-4_dp, "an all-electron run gives the frozen-core energy "// &
                 to_string(e_ccsd_none, "f20.14"))
      rc = trc_destroy(h)
   end subroutine test_frozen_by_element

   subroutine test_no_triples(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      type(rccsd_refs_t) :: refs
      real(c_double) :: e_mp2, e_ccsd, e_t, t1, e_total, e_scf_part
      integer(c_int) :: rc, iterations, converged, has_triples
      logical :: ok

      call water_scf(error, h, 0, 1, ok)
      if (.not. ok .or. allocated(error)) return
      call read_rccsd_refs(REF, refs)
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "trc_run_rccsd: "//message(h))
      if (allocated(error)) return
      e_t = 1.0_c_double
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, has_triples == 0 .and. e_t == 0.0_c_double, "(T) reported although it was not asked for")
      if (allocated(error)) return
      rc = trc_energy(h, e_total)
      e_scf_part = e_total - e_ccsd
      call check(error, abs(e_scf_part - refs%e_scf) < 1.0e-9_dp, "trc_energy minus E_CCSD is "// &
                 to_string(e_scf_part, "f20.12")//", not the SCF energy "//to_string(refs%e_scf, "f20.12"))
      if (allocated(error)) return
      call check(error, abs(e_ccsd - refs%e_ccsd_corr) < TOL_CORR, "E_CCSD "//to_string(e_ccsd, "f20.12"))
      rc = trc_destroy(h)
   end subroutine test_no_triples

   subroutine test_setters(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      type(rccsd_refs_t) :: refs
      real(c_double) :: e_mp2, e_ccsd, e_t, t1
      integer(c_int) :: rc, iterations, converged, has_triples, nf, na, nvv, nx
      logical :: ok

      call water_scf(error, h, 1, 1, ok)
      if (.not. ok .or. allocated(error)) return
      call read_rccsd_refs(REF, refs)
      call refused(error, "triples 2", trc_set_rccsd(h, 2_c_int, 1_c_int, 0_c_int), h, "triples")
      call refused(error, "triples -1", trc_set_rccsd(h, -1_c_int, 1_c_int, 0_c_int), h, "triples")
      call refused(error, "nfrozen -2", trc_set_rccsd(h, 1_c_int, -2_c_int, 0_c_int), h, "nfrozen")
      call refused(error, "aux_block -1", trc_set_rccsd(h, 1_c_int, 1_c_int, -1_c_int), h, "aux_block")
      call refused(error, "conv_energy -1", trc_set_rccsd_convergence(h, -1.0_c_double, 0.0_c_double, 0_c_int, 0_c_int), h, &
                   "conv_energy")
      call refused(error, "conv_energy NaN", trc_set_rccsd_convergence(h, nan(), 0.0_c_double, 0_c_int, 0_c_int), h, &
                   "conv_energy")
      call refused(error, "conv_energy Inf", trc_set_rccsd_convergence(h, infinity(), 0.0_c_double, 0_c_int, 0_c_int), h, &
                   "conv_energy")
      call refused(error, "conv_amplitude -1", trc_set_rccsd_convergence(h, 0.0_c_double, -1.0_c_double, 0_c_int, 0_c_int), &
                   h, "conv_amplitude")
      call refused(error, "conv_amplitude NaN", trc_set_rccsd_convergence(h, 0.0_c_double, nan(), 0_c_int, 0_c_int), h, &
                   "conv_amplitude")
      call refused(error, "conv_amplitude Inf", trc_set_rccsd_convergence(h, 0.0_c_double, infinity(), 0_c_int, 0_c_int), h, &
                   "conv_amplitude")
      call refused(error, "max_iterations -1", trc_set_rccsd_convergence(h, 0.0_c_double, 0.0_c_double, -1_c_int, 0_c_int), &
                   h, "max_iterations")
      call refused(error, "diis_size -1", trc_set_rccsd_convergence(h, 0.0_c_double, 0.0_c_double, 0_c_int, -1_c_int), h, &
                   "diis_size")
      call refused(error, "ladder 3", trc_set_rccsd_algorithm(h, 3_c_int, 0_c_int, 0.0_c_double), h, "ladder")
      call refused(error, "ladder -1", trc_set_rccsd_algorithm(h, -1_c_int, 0_c_int, 0.0_c_double), h, "ladder")
      call refused(error, "budget -1", trc_set_rccsd_algorithm(h, 0_c_int, 0_c_int, -1.0_c_double), h, &
                   "triples_budget_mib")
      if (allocated(error)) return
      ! One call per statement: the order in which a compiler evaluates the operands of one expression is its own.
      rc = trc_set_rccsd(c_null_handle(), 1_c_int, 1_c_int, 0_c_int)
      call check(error, rc == TRC_ERR_NULL, "trc_set_rccsd on a null handle is not TRC_ERR_NULL")
      rc = trc_set_rccsd_convergence(c_null_handle(), 0.0_c_double, 0.0_c_double, 0_c_int, 0_c_int)
      call check(error, rc == TRC_ERR_NULL, "trc_set_rccsd_convergence on a null handle is not TRC_ERR_NULL")
      rc = trc_set_rccsd_algorithm(c_null_handle(), 0_c_int, 0_c_int, 0.0_c_double)
      call check(error, rc == TRC_ERR_NULL, "trc_set_rccsd_algorithm on a null handle is not TRC_ERR_NULL")
      rc = trc_run_rccsd(c_null_handle())
      call check(error, rc == TRC_ERR_NULL, "trc_run_rccsd on a null handle is not TRC_ERR_NULL")
      rc = trc_rccsd_dimensions(c_null_handle(), nf, na, nvv, nx)
      call check(error, rc == TRC_ERR_NULL, "trc_rccsd_dimensions on a null handle is not TRC_ERR_NULL")
      if (allocated(error)) return
      ! none of the refusals above changed the settings water_scf made: the run still reproduces pyscf
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "trc_run_rccsd after the refusals: "//message(h))
      if (allocated(error)) return
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, has_triples == 1 .and. abs(e_t - refs%e_t) < TOL_T .and. abs(e_ccsd - refs%e_ccsd_corr) < TOL_CORR, &
                 "a refused setter changed the settings: E_CCSD "//to_string(e_ccsd, "f20.12")//", E(T) "// &
                 to_string(e_t, "f20.14"))
      ! the accepted extremes
      rc = trc_set_rccsd_algorithm(h, 1_c_int, 0_c_int, 2048.0_c_double)
      call check(error, rc == TRC_OK, "ladder 1 with a 2048 MiB tile was refused")
      rc = trc_set_rccsd_algorithm(h, 2_c_int, 0_c_int, 0.0_c_double)
      call check(error, rc == TRC_OK, "ladder 2 was refused")
      rc = trc_set_rccsd_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 100_c_int, 1_c_int)
      call check(error, rc == TRC_OK, "a DIIS size of 1 was refused")
      rc = trc_set_rccsd(h, 1_c_int, 1_c_int, 7_c_int)
      call check(error, rc == TRC_OK, "an aux_block of 7 was refused")
      rc = trc_destroy(h)
   end subroutine test_setters

   function c_null_handle() result(h)
      type(c_ptr) :: h
      h = c_null_ptr
   end function c_null_handle

   subroutine refused(error, what, rc, h, word)
      !! `rc` is TRC_ERR_BADARG and the handle's message names `word`
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: what, word
      integer(c_int), intent(in) :: rc
      type(c_ptr), intent(in) :: h
      if (allocated(error)) return
      call check(error, rc == TRC_ERR_BADARG, what//": status "//to_string(int(rc))//", not TRC_ERR_BADARG")
      if (allocated(error)) return
      call check(error, index(message(h), word) > 0, what//": the message '"//message(h)//"' does not name '"//word//"'")
   end subroutine refused

   subroutine expect(error, what, rc, want, h, word)
      !! `rc` is `want` and the handle's message contains `word`
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: what, word
      integer(c_int), intent(in) :: rc, want
      type(c_ptr), intent(in) :: h
      if (allocated(error)) return
      call check(error, rc == want, what//": status "//to_string(int(rc))//", wanted "//to_string(int(want))// &
                 " ("//message(h)//")")
      if (allocated(error)) return
      call check(error, index(message(h), word) > 0, what//": the message '"//message(h)//"' does not contain '"//word//"'")
   end subroutine expect

   subroutine test_state_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      logical :: ok

      ok = supported()
      if (.not. ok) then
         print "(a)", "SKIP: this build's density-fitting ceilings are below cc-pVDZ-RIFIT"
         return
      end if
      ! no SCF
      call make_handle(h, WATER, "cc-pvdz", "cc-pvdz-rifit", 0, 1, rc)
      call expect(error, "no SCF", trc_run_rccsd(h), TRC_ERR_STATE, h, "run the SCF first")
      if (allocated(error)) return
      ! an SCF that stopped at one iteration
      rc = trc_set_convergence(h, 1.0e-13_c_double, 1.0e-11_c_double, 1_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_ERR_NOCONV, "the one-iteration SCF returned "//to_string(int(rc)))
      call expect(error, "unconverged SCF", trc_run_rccsd(h), TRC_ERR_STATE, h, "did not converge")
      rc = trc_destroy(h)
      if (allocated(error)) return
      ! no auxiliary basis
      call make_handle(h, WATER, "cc-pvdz", "", 0, 1, rc)
      rc = trc_run_scf(h)
      call expect(error, "no auxiliary basis", trc_run_rccsd(h), TRC_ERR_STATE, h, "no auxiliary basis")
      rc = trc_destroy(h)
      if (allocated(error)) return
   end subroutine test_state_refusals

   subroutine test_reference_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      logical :: ok

      ok = supported()
      if (.not. ok) then
         print "(a)", "SKIP: this build's density-fitting ceilings are below cc-pVDZ-RIFIT"
         return
      end if
      ! an unrestricted reference (closed-shell water forced to UHF), and a genuine open shell
      call make_handle(h, WATER, "6-31g", "cc-pvdz-rifit", 0, 1, rc)
      rc = trc_set_unrestricted(h, 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the UHF water SCF failed")
      if (allocated(error)) return
      call expect(error, "UHF reference", trc_run_rccsd(h), TRC_ERR_UNSUPPORTED, h, "closed-shell only")
      rc = trc_destroy(h)
      if (allocated(error)) return
      call make_handle(h, "oh.xyz", "6-31g", "cc-pvdz-rifit", 0, 2, rc)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the OH SCF failed")
      if (allocated(error)) return
      call expect(error, "open shell", trc_run_rccsd(h), TRC_ERR_UNSUPPORTED, h, "closed-shell only")
      rc = trc_destroy(h)
      if (allocated(error)) return
      ! a Kohn-Sham reference
      call make_handle(h, WATER, "6-31g", "cc-pvdz-rifit", 0, 1, rc)
      rc = trc_set_method(h, "pbe"//c_null_char, 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the PBE water SCF failed")
      if (allocated(error)) return
      call expect(error, "Kohn-Sham reference", trc_run_rccsd(h), TRC_ERR_UNSUPPORTED, h, "Kohn-Sham")
      rc = trc_destroy(h)
      if (allocated(error)) return
      ! an effective core potential: HI in def2-SVP carries the basis's own
      rc = trc_create(h)
      rc = trc_set_molecule(h, 2_c_int, [53.0_c_double, 1.0_c_double], &
                            [0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, 0.0_c_double, HI_BOND_BOHR], &
                            0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/def2-svp.json'//c_null_char)
      rc = trc_set_aux_json(h, 'basis_sets/def2-svp-rifit.json'//c_null_char)
      rc = trc_set_convergence(h, 1.0e-8_c_double, 1.0e-6_c_double, 150_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the HI SCF failed: "//message(h))
      if (allocated(error)) return
      call expect(error, "ECP", trc_run_rccsd(h), TRC_ERR_UNSUPPORTED, h, "effective core potential")
      rc = trc_destroy(h)
      if (allocated(error)) return
   end subroutine test_reference_refusals

   function trc_rimp2_energy_of(h, e_os, e_ss) result(rc)
      !! trc_rimp2_energy settles the deferred pair energies, so the run is complete
      use trc_c_interfaces, only: trc_rimp2_energy
      type(c_ptr), intent(in) :: h
      real(c_double), intent(out) :: e_os, e_ss
      integer(c_int) :: rc
      rc = trc_rimp2_energy(h, e_os, e_ss)
   end function trc_rimp2_energy_of

   subroutine test_size_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      logical :: ok

      call water_scf(error, h, 0, 5, ok)
      if (.not. ok .or. allocated(error)) return
      call expect(error, "every orbital frozen", trc_run_rccsd(h), TRC_ERR_BADARG, h, "frozen")
      rc = trc_destroy(h)
      if (allocated(error)) return
      ! cc-pVQZ-RIFIT goes to l = 7, above the shipped auxiliary ceiling of 4
      if (DF_LMAX_AUX >= 7) then
         print "(a)", "SKIP: this build's auxiliary ceiling covers cc-pVQZ-RIFIT"
         return
      end if
      call make_handle(h, WATER, "6-31g", "../basis_sets/cc-pvqz-rifit", 0, 1, rc)
      rc = trc_run_scf(h)
      call expect(error, "auxiliary l above the ceiling", trc_run_rccsd(h), TRC_ERR_UNSUPPORTED, h, "TERCO_DF_LMAX_AUX")
      rc = trc_destroy(h)
      if (allocated(error)) return
   end subroutine test_size_refusals

   subroutine test_pcm(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(c_double) :: e_mp2, e_ccsd, e_t, t1
      integer(c_int) :: rc, iterations, converged, has_triples
      logical :: ok

      ok = supported()
      if (.not. ok) then
         print "(a)", "SKIP: this build's density-fitting ceilings are below cc-pVDZ-RIFIT"
         return
      end if
      call make_handle(h, WATER, "6-31g", "cc-pvdz-rifit", 0, 1, rc)
      rc = trc_set_pcm(h, "cpcm"//c_null_char, 78.3553_c_double, 0_c_int, 0.0_c_double)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the PCM SCF failed")
      if (allocated(error)) return
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "RI-CCSD over a PCM reference: "//message(h))
      if (allocated(error)) return
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, converged == 1 .and. e_ccsd < -0.05_dp .and. e_ccsd > e_mp2 - 0.05_dp, &
                 "E_CCSD over PCM "//to_string(e_ccsd, "f20.12"))
      rc = trc_destroy(h)
   end subroutine test_pcm

   subroutine test_invalidation(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(c_double) :: e_mp2, e_ccsd, e_t, t1, e_total, e_scf
      integer(c_int) :: rc, rc_set, rc_get, iterations, converged, has_triples
      logical :: ok

      call water_scf(error, h, 1, 1, ok)
      if (.not. ok .or. allocated(error)) return
      ! a second SCF
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "the first run: "//message(h))
      if (allocated(error)) return
      rc = trc_run_scf(h)
      call check(error, trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples) == TRC_ERR_STATE, &
                 "a second trc_run_scf left the RI-CCSD result readable")
      if (allocated(error)) return
      rc = trc_energy(h, e_scf)
      rc = trc_run_rccsd(h)
      rc = trc_energy(h, e_total)
      call check(error, e_total < e_scf - 0.1_dp, "the rerun is not a correlated energy")
      if (allocated(error)) return
      ! a new auxiliary basis
      rc = trc_set_aux_json(h, 'basis_sets/def2-svp-rifit.json'//c_null_char)
      call check(error, trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples) == TRC_ERR_STATE, &
                 "a new auxiliary basis left the RI-CCSD result readable")
      if (allocated(error)) return
      rc = trc_energy(h, e_total)
      call check(error, abs(e_total - e_scf) < 1.0e-12_dp, "trc_energy still adds a correlation energy after a new auxiliary basis")
      if (allocated(error)) return
      ! each setter
      rc = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      rc = trc_run_rccsd(h)
      ! Two statements: the order in which the operands of one expression are evaluated is the
      ! compiler's choice, and nvfortran read the result before the setter had run.
      rc_set = trc_set_rccsd(h, 1_c_int, 1_c_int, 0_c_int)
      rc_get = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, rc_set == TRC_OK .and. rc_get == TRC_ERR_STATE, &
                 "trc_set_rccsd left the result readable")
      if (allocated(error)) return
      rc = trc_run_rccsd(h)
      ! Two statements: the order in which the operands of one expression are evaluated is the
      ! compiler's choice, and nvfortran read the result before the setter had run.
      rc_set = trc_set_rccsd_convergence(h, 0.0_c_double, 0.0_c_double, 0_c_int, 0_c_int)
      rc_get = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, rc_set == TRC_OK .and. rc_get == TRC_ERR_STATE, &
                 "trc_set_rccsd_convergence left the result readable")
      if (allocated(error)) return
      rc = trc_run_rccsd(h)
      ! Two statements: the order in which the operands of one expression are evaluated is the
      ! compiler's choice, and nvfortran read the result before the setter had run.
      rc_set = trc_set_rccsd_algorithm(h, 0_c_int, 0_c_int, 0.0_c_double)
      rc_get = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, rc_set == TRC_OK .and. rc_get == TRC_ERR_STATE, &
                 "trc_set_rccsd_algorithm left the result readable")
      if (allocated(error)) return
      ! RI-MP2 afterwards replaces it, and the energy is then the RI-MP2 one
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "the run before RI-MP2: "//message(h))
      if (allocated(error)) return
      rc = trc_run_rimp2(h)
      call check(error, trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples) == TRC_ERR_STATE, &
                 "trc_run_rimp2 left the RI-CCSD result readable")
      rc = trc_destroy(h)
   end subroutine test_invalidation

   subroutine iterations_with(h, conv_energy, conv_amplitude, diis_size, n, e_ccsd)
      !! The iteration count and energy of one CCSD run with these settings on the SCF the handle holds
      type(c_ptr), intent(in) :: h
      real(c_double), intent(in) :: conv_energy, conv_amplitude
      integer(c_int), intent(in) :: diis_size
      integer, intent(out) :: n
      real(c_double), intent(out) :: e_ccsd
      real(c_double) :: e_mp2, e_t, t1
      integer(c_int) :: rc, iterations, converged, has_triples
      rc = trc_set_rccsd_convergence(h, conv_energy, conv_amplitude, 100_c_int, diis_size)
      rc = trc_run_rccsd(h)
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      n = int(iterations)
   end subroutine iterations_with

   subroutine test_convergence_settings(error)
      !! The stopping rule is energy AND amplitude, so tightening either alone must cost iterations, and a history
      !! of one (no DIIS) must cost more than the default of eight
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(c_double) :: e_loose, e_tight_energy, e_tight_amplitude, e_nodiis, e_default
      integer :: n_loose, n_tight_energy, n_tight_amplitude, n_nodiis, n_default
      integer(c_int) :: rc
      logical :: ok

      call water_scf(error, h, 0, 1, ok)
      if (.not. ok .or. allocated(error)) return
      call iterations_with(h, 1.0e-3_c_double, 1.0e-2_c_double, 0_c_int, n_loose, e_loose)
      call iterations_with(h, 1.0e-12_c_double, 1.0e-2_c_double, 0_c_int, n_tight_energy, e_tight_energy)
      call iterations_with(h, 1.0e-3_c_double, 1.0e-10_c_double, 0_c_int, n_tight_amplitude, e_tight_amplitude)
      call check(error, n_tight_energy > n_loose, "conv_energy 1e-12 took "//to_string(n_tight_energy)// &
                 " iterations, the loose run "//to_string(n_loose))
      if (allocated(error)) return
      call check(error, n_tight_amplitude > n_loose, "conv_amplitude 1e-10 took "//to_string(n_tight_amplitude)// &
                 " iterations, the loose run "//to_string(n_loose))
      if (allocated(error)) return
      call iterations_with(h, 1.0e-12_c_double, 1.0e-10_c_double, 0_c_int, n_default, e_default)
      call iterations_with(h, 1.0e-12_c_double, 1.0e-10_c_double, 1_c_int, n_nodiis, e_nodiis)
      call check(error, n_nodiis > n_default, "diis_size 1 took "//to_string(n_nodiis)//" iterations, the default "// &
                 to_string(n_default))
      if (allocated(error)) return
      call check(error, abs(e_nodiis - e_default) < 1.0e-9_dp, "the converged energy depends on the DIIS size: "// &
                 to_string(e_nodiis - e_default, "es9.2"))
      rc = trc_destroy(h)
   end subroutine test_convergence_settings

   subroutine test_noconv(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(c_double) :: e_mp2, e_ccsd, e_t, t1, e_total, e_scf
      integer(c_int) :: rc, iterations, converged, has_triples
      logical :: ok

      call water_scf(error, h, 1, 1, ok)
      if (.not. ok .or. allocated(error)) return
      rc = trc_energy(h, e_scf)
      rc = trc_set_rccsd_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 2_c_int, 0_c_int)
      rc = trc_run_rccsd(h)
      call expect(error, "two iterations", rc, TRC_ERR_NOCONV, h, "not converged")
      if (allocated(error)) return
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, rc == TRC_OK .and. converged == 0 .and. iterations == 2 .and. has_triples == 0 .and. e_t == 0.0_c_double, &
                 "the unconverged result: status "//to_string(int(rc))//", converged "//to_string(int(converged))// &
                 ", iterations "//to_string(int(iterations))//", has_triples "//to_string(int(has_triples)))
      if (allocated(error)) return
      call check(error, e_ccsd < e_mp2 .and. e_ccsd > e_mp2 - 0.05_dp, "the last iterate E_CCSD "//to_string(e_ccsd, "f20.12"))
      if (allocated(error)) return
      rc = trc_energy(h, e_total)
      call check(error, abs(e_total - e_scf) < 1.0e-12_dp, "trc_energy includes an unconverged correlation energy")
      rc = trc_destroy(h)
   end subroutine test_noconv

   subroutine test_gradient_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      logical :: ok

      call water_scf(error, h, 0, 1, ok)
      if (.not. ok .or. allocated(error)) return
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "trc_run_rccsd: "//message(h))
      if (allocated(error)) return
      rc = trc_run_gradient(h)
      call expect(error, "gradient after RI-CCSD", rc, TRC_ERR_UNSUPPORTED, h, "RI-CCSD")
      if (allocated(error)) return
      call check(error, index(message(h), "RI-CCSD(T)") == 0, "the CCSD-only refusal names (T): "//message(h))
      if (allocated(error)) return
      ! with (T) the message says so
      rc = trc_set_rccsd(h, 1_c_int, 1_c_int, 0_c_int)
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "trc_run_rccsd with (T): "//message(h))
      if (allocated(error)) return
      rc = trc_run_gradient(h)
      call expect(error, "gradient after RI-CCSD(T)", rc, TRC_ERR_UNSUPPORTED, h, "RI-CCSD(T)")
      rc = trc_destroy(h)
      if (allocated(error)) return
      ! an unconverged CCSD leaves trc_energy at the SCF energy, so the SCF gradient is legitimate
      call water_scf(error, h, 0, 1, ok)
      if (.not. ok .or. allocated(error)) return
      rc = trc_set_rccsd_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 2_c_int, 0_c_int)
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_ERR_NOCONV, "two iterations returned "//to_string(int(rc)))
      if (allocated(error)) return
      rc = trc_run_gradient(h)
      call check(error, rc == TRC_OK, "the SCF gradient after an unconverged CCSD: "//message(h))
      rc = trc_destroy(h)
   end subroutine test_gradient_refused

   subroutine test_methods_alternate(error)
      !! One correlated method owns trc_energy: RI-MP2 after RI-CCSD replaces it and the reverse, and a run that is
      !! refused leaves the other's result alone.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(c_double) :: e_mp2, e_ccsd, e_t, t1, e_scf, e_total, e_os, e_ss
      integer(c_int) :: rc, iterations, converged, has_triples
      logical :: ok

      call water_scf(error, h, 0, 1, ok)
      if (.not. ok .or. allocated(error)) return
      rc = trc_energy(h, e_scf)
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "trc_run_rccsd: "//message(h))
      if (allocated(error)) return
      rc = trc_energy(h, e_total)
      call check(error, e_total < e_scf - 0.1_dp, "the CCSD total is not below the SCF energy")
      if (allocated(error)) return
      ! RI-MP2 afterwards: not refused, and the energy is the RI-MP2 one
      rc = trc_set_rimp2(h, 1_c_int, 0_c_int)
      rc = trc_run_rimp2(h)
      call check(error, rc == TRC_OK, "trc_run_rimp2 after RI-CCSD: "//message(h))
      if (allocated(error)) return
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, rc == TRC_ERR_STATE, "RI-MP2 left the RI-CCSD result readable")
      if (allocated(error)) return
      rc = trc_rimp2_energy_of(h, e_os, e_ss)
      call check(error, rc == TRC_OK, "the RI-MP2 energy is not readable after its own run")
      if (allocated(error)) return
      rc = trc_energy(h, e_total)
      call check(error, abs(e_total - (e_scf + e_os + e_ss)) < 1.0e-12_dp, "trc_energy after RI-MP2 is not SCF + RI-MP2")
      if (allocated(error)) return
      ! RI-CCSD afterwards: not refused, and it drops the RI-MP2 result
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "trc_run_rccsd after RI-MP2: "//message(h))
      if (allocated(error)) return
      rc = trc_rimp2_energy_of(h, e_os, e_ss)
      call check(error, rc == TRC_ERR_STATE, "RI-CCSD left the RI-MP2 result readable")
      if (allocated(error)) return
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      rc = trc_energy(h, e_total)
      call check(error, abs(e_total - (e_scf + e_ccsd)) < 1.0e-12_dp, "trc_energy is not SCF + RI-CCSD after the rerun")
      if (allocated(error)) return
      ! a refused RI-MP2 (every occupied orbital frozen) keeps the RI-CCSD result
      rc = trc_set_rimp2(h, 5_c_int, 0_c_int)
      rc = trc_run_rimp2(h)
      call check(error, rc == TRC_ERR_BADARG, "trc_run_rimp2 with everything frozen returned "//to_string(int(rc)))
      if (allocated(error)) return
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
      call check(error, rc == TRC_OK, "a refused trc_run_rimp2 dropped the RI-CCSD result")
      if (allocated(error)) return
      ! and a refused RI-CCSD keeps the RI-MP2 result
      rc = trc_set_rimp2(h, 1_c_int, 0_c_int)
      rc = trc_run_rimp2(h)
      call check(error, rc == TRC_OK, "trc_run_rimp2: "//message(h))
      if (allocated(error)) return
      rc = trc_set_rccsd(h, 0_c_int, 5_c_int, 0_c_int)
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_ERR_BADARG, "trc_run_rccsd with everything frozen returned "//to_string(int(rc)))
      if (allocated(error)) return
      rc = trc_rimp2_energy_of(h, e_os, e_ss)
      call check(error, rc == TRC_OK, "a refused trc_run_rccsd dropped the RI-MP2 result")
      rc = trc_destroy(h)
   end subroutine test_methods_alternate

   subroutine ccsd_numbers(h, e_ccsd, e_t, converged, has_triples, rc)
      !! The correlation energies of the stored RI-CCSD result and its two flags
      type(c_ptr), intent(in) :: h
      real(c_double), intent(out) :: e_ccsd, e_t
      integer(c_int), intent(out) :: converged, has_triples, rc
      real(c_double) :: e_mp2, t1
      integer(c_int) :: iterations
      rc = trc_rccsd_energy(h, e_mp2, e_ccsd, e_t, t1, iterations, converged, has_triples)
   end subroutine ccsd_numbers

   subroutine set_triples(error, h, routine, budget_mib)
      !! The (T) routine and budget of the handle, ladder at its default
      type(error_type), allocatable, intent(inout) :: error
      type(c_ptr), intent(in) :: h
      integer(c_int), intent(in) :: routine
      real(c_double), intent(in) :: budget_mib
      integer(c_int) :: rc
      if (allocated(error)) return
      rc = trc_set_rccsd_algorithm(h, 0_c_int, routine, budget_mib)
      call check(error, rc == TRC_OK, "trc_set_rccsd_algorithm(routine "//to_string(int(routine))//", "// &
                 to_string(real(budget_mib, dp), "es9.2")//" MiB): "//message(h))
   end subroutine set_triples

   subroutine run_for_triples(error, h, routine, budget_mib, e_ccsd, e_t)
      !! A run with the given routine and budget that must succeed with (T): its correlation energies
      type(error_type), allocatable, intent(inout) :: error
      type(c_ptr), intent(in) :: h
      integer(c_int), intent(in) :: routine
      real(c_double), intent(in) :: budget_mib
      real(c_double), intent(out) :: e_ccsd, e_t
      integer(c_int) :: rc, converged, has_triples
      character(len=*), parameter :: ASKED = "routine "
      e_ccsd = 0.0_c_double
      e_t = 0.0_c_double
      call set_triples(error, h, routine, budget_mib)
      if (allocated(error)) return
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, ASKED//to_string(int(routine))//", "//to_string(real(budget_mib, dp), "es9.2")// &
                 " MiB: trc_run_rccsd returned "//to_string(int(rc))//": "//message(h))
      if (allocated(error)) return
      call ccsd_numbers(h, e_ccsd, e_t, converged, has_triples, rc)
      call check(error, rc == TRC_OK .and. converged == 1 .and. has_triples == 1, ASKED//to_string(int(routine))// &
                 ": the run did not leave a converged (T) result")
   end subroutine run_for_triples

   subroutine test_triples_refused(error)
      !! One water triple takes 48*4^3 = 3072 bytes of W; 0.001 MiB is 1048 bytes
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      type(rccsd_refs_t) :: refs
      real(c_double) :: e_ccsd_off, e_ccsd, e_t, e_total, e_total_off, e_t_default, e_t_tiled, e_t_ref
      integer(c_int) :: rc, rc_run, converged, has_triples, nf, na, nvv, nx
      logical :: ok

      call water_scf(error, h, 0, 1, ok)
      if (.not. ok .or. allocated(error)) return
      call read_rccsd_refs(REF, refs)
      ! the CCSD alone on this handle's settings, and its total
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "the triples-off run: "//message(h))
      if (allocated(error)) return
      call ccsd_numbers(h, e_ccsd_off, e_t, converged, has_triples, rc)
      rc = trc_energy(h, e_total_off)
      ! the same handle with (T) on, the tiled routine and a budget below one triple
      rc = trc_set_rccsd(h, 1_c_int, 1_c_int, 0_c_int)
      call set_triples(error, h, 2_c_int, 0.001_c_double)
      if (allocated(error)) return
      rc_run = trc_run_rccsd(h)
      call check(error, rc_run == TRC_ERR_UNSUPPORTED, "a (T) step with no room for one triple returned "// &
                 to_string(int(rc_run))//", not TRC_ERR_UNSUPPORTED")
      if (allocated(error)) return
      call check(error, index(message(h), "(T) failed") > 0 .and. index(message(h), "budget") > 0, &
                 "the message does not say that (T) failed for its budget: '"//message(h)//"'")
      if (allocated(error)) return
      ! One call per statement: the order in which a compiler evaluates the operands of one expression is its own.
      e_t = 1.0_c_double
      call ccsd_numbers(h, e_ccsd, e_t, converged, has_triples, rc)
      call check(error, rc == TRC_OK .and. converged == 1 .and. has_triples == 0, &
                 "after the refused (T): trc_rccsd_energy status "//to_string(int(rc))//", converged "// &
                 to_string(int(converged))//", has_triples "//to_string(int(has_triples)))
      if (allocated(error)) return
      call check(error, abs(e_ccsd - e_ccsd_off) < TOL_ROUTINES, "the RI-CCSD correlation after the refused (T) is "// &
                 to_string(e_ccsd, "f20.12")//", the triples-off run's "//to_string(e_ccsd_off, "f20.12"))
      if (allocated(error)) return
      call check(error, e_t == 0.0_c_double, "e_t is "//to_string(e_t, "es9.2")//" after the refused (T)")
      if (allocated(error)) return
      rc = trc_energy(h, e_total)
      call check(error, rc == TRC_OK .and. abs(e_total - (e_total_off - e_ccsd_off + e_ccsd)) < TOL_ROUTINES .and. &
                 abs(e_total - e_total_off) < TOL_ROUTINES .and. abs(e_total - e_ccsd - refs%e_scf) < 1.0e-9_dp, &
                 "trc_energy after the refused (T) is "//to_string(e_total, "f20.12")//", not E_SCF + E_CCSD "// &
                 to_string(e_total_off, "f20.12"))
      if (allocated(error)) return
      rc = trc_rccsd_dimensions(h, nf, na, nvv, nx)
      call check(error, rc == TRC_OK .and. nf == 1 .and. na == 4 .and. nvv == 19 .and. nx == 84, &
                 "trc_rccsd_dimensions after the refused (T): status "//to_string(int(rc))//", "//to_string(int(na))// &
                 " active occupied")
      if (allocated(error)) return
      ! enough budget: (T), and the default routine of the build gives the same number
      call run_for_triples(error, h, 2_c_int, 64.0_c_double, e_ccsd, e_t_tiled)
      call run_for_triples(error, h, 0_c_int, 0.0_c_double, e_ccsd, e_t_default)
      ! the reference routine has no budget, so the budget that stopped the tiled one does not stop it
      call run_for_triples(error, h, 1_c_int, 0.001_c_double, e_ccsd, e_t_ref)
      if (allocated(error)) return
      call check(error, abs(e_t_tiled - e_t_default) < TOL_ROUTINES .and. abs(e_t_ref - e_t_default) < TOL_ROUTINES, &
                 "E(T): tiled "//to_string(e_t_tiled, "f20.14")//", reference with the tiny budget "// &
                 to_string(e_t_ref, "f20.14")//", default "//to_string(e_t_default, "f20.14"))
      if (allocated(error)) return
      call check(error, abs(e_t_default - refs%e_t) < TOL_T, "E(T) "//to_string(e_t_default, "f20.14")//" vs pyscf "// &
                 to_string(refs%e_t, "f20.14"))
      rc = trc_destroy(h)
      print "(a)", "  measured: |tiled - default| "//to_string(abs(e_t_tiled - e_t_default), "es9.2")// &
         ", |reference (0.001 MiB) - default| "//to_string(abs(e_t_ref - e_t_default), "es9.2")
   end subroutine test_triples_refused

   subroutine test_triples_routines(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(c_double) :: e_ccsd, e_t_ref, e_t_tiled, e_t_frac
      integer(c_int) :: rc
      logical :: ok

      call water_scf(error, h, 1, 1, ok)
      if (.not. ok .or. allocated(error)) return
      call refused(error, "routine 3", trc_set_rccsd_algorithm(h, 0_c_int, 3_c_int, 0.0_c_double), h, "triples_routine")
      call refused(error, "routine -1", trc_set_rccsd_algorithm(h, 0_c_int, -1_c_int, 0.0_c_double), h, "triples_routine")
      call refused(error, "budget NaN", trc_set_rccsd_algorithm(h, 0_c_int, 0_c_int, nan()), h, "triples_budget_mib")
      call refused(error, "budget Inf", trc_set_rccsd_algorithm(h, 0_c_int, 0_c_int, infinity()), h, "triples_budget_mib")
      call refused(error, "budget -0.5", trc_set_rccsd_algorithm(h, 0_c_int, 0_c_int, -0.5_c_double), h, "triples_budget_mib")
      call refused(error, "budget 1e13", trc_set_rccsd_algorithm(h, 0_c_int, 0_c_int, 1.0e13_c_double), h, "triples_budget_mib")
      if (allocated(error)) return
      call run_for_triples(error, h, 1_c_int, 0.0_c_double, e_ccsd, e_t_ref)
      call run_for_triples(error, h, 2_c_int, 0.0_c_double, e_ccsd, e_t_tiled)
      if (allocated(error)) return
      call check(error, abs(e_t_ref - e_t_tiled) < TOL_ROUTINES, "E(T): reference "//to_string(e_t_ref, "f20.14")// &
                 ", tiled "//to_string(e_t_tiled, "f20.14")//", |difference| "//to_string(abs(e_t_ref - e_t_tiled), "es9.2"))
      if (allocated(error)) return
      ! a fractional budget that holds a few triples is accepted and changes nothing but the tile size
      call run_for_triples(error, h, 2_c_int, 0.01_c_double, e_ccsd, e_t_frac)
      if (allocated(error)) return
      call check(error, abs(e_t_frac - e_t_tiled) < TOL_ROUTINES, "a 0.01 MiB tile moved E(T) by "// &
                 to_string(e_t_frac - e_t_tiled, "es9.2"))
      rc = trc_destroy(h)
      print "(a)", "  measured: |reference - tiled| "//to_string(abs(e_t_ref - e_t_tiled), "es9.2")// &
         ", |0.01 MiB tile - tiled| "//to_string(abs(e_t_frac - e_t_tiled), "es9.2")
   end subroutine test_triples_routines

   subroutine test_dimensions(error)
      !! Water cc-pVDZ: 24 orbitals, 5 occupied, 1 frozen oxygen 1s; cc-pVDZ-RIFIT has 84 functions
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, nf, na, nvv, nx
      logical :: ok

      call water_scf(error, h, 0, -1, ok)
      if (.not. ok .or. allocated(error)) return
      rc = trc_rccsd_dimensions(h, nf, na, nvv, nx)
      call check(error, rc == TRC_ERR_STATE, "dimensions before any run: "//to_string(int(rc)))
      if (allocated(error)) return
      rc = trc_run_rccsd(h)
      call check(error, rc == TRC_OK, "trc_run_rccsd: "//message(h))
      if (allocated(error)) return
      rc = trc_rccsd_dimensions(h, nf, na, nvv, nx)
      call check(error, rc == TRC_OK .and. nf == 1 .and. na == 4 .and. nvv == 19 .and. nx == 84, &
                 "dimensions "//to_string(int(nf))//", "//to_string(int(na))//", "//to_string(int(nvv))//", "// &
                 to_string(int(nx))//": wanted 1 frozen, 4 active occupied, 19 virtual and 84 auxiliary functions")
      if (allocated(error)) return
      rc = trc_set_rccsd(h, 0_c_int, 0_c_int, 0_c_int)
      rc = trc_run_rccsd(h)
      rc = trc_rccsd_dimensions(h, nf, na, nvv, nx)
      call check(error, nf == 0 .and. na == 5 .and. nvv == 19, "an all-electron run reports "//to_string(int(nf))//" frozen, "// &
                 to_string(int(na))//" active occupied")
      rc = trc_destroy(h)
   end subroutine test_dimensions

   subroutine test_lindep(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      logical :: ok

      ok = supported()
      if (.not. ok) then
         print "(a)", "SKIP: this build's density-fitting ceilings are below cc-pVDZ-RIFIT"
         return
      end if
      call make_handle(h, WATER, "cc-pvdz", "cc-pvdz-rifit", 0, 1, rc)
      ! a threshold far above the basis's smallest overlap eigenvalue drops functions
      rc = trc_set_lindep(h, 0.3_c_double)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the SCF with a 0.3 linear-dependence threshold failed: "//message(h))
      if (allocated(error)) return
      rc = trc_run_rccsd(h)
      call expect(error, "dropped functions", rc, TRC_ERR_UNSUPPORTED, h, "trc_set_lindep")
      if (allocated(error)) return
      call check(error, index(message(h), "--lindep") > 0 .and. index(message(h), "linear_dependence") > 0, &
                 "the message does not name the command-line option and the input key: "//message(h))
      rc = trc_destroy(h)
   end subroutine test_lindep

end module check_capi_rccsd_suite

program check_capi_rccsd
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_capi_rccsd_suite, only: collect_capi_rccsd
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_capi_rccsd, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_rccsd
