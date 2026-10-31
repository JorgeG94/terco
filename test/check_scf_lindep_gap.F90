!
! Canonical orthogonalisation's dropped directions, and the HOMO-LUMO gap.
!
! X zeroes the directions it drops, so X^T F X used to have an eigenvalue 0
! for each of them: a zero orbital. Below a positive HOMO aufbau occupied
! it and the density lost two electrons. O2- (a dianion, HOMO +0.12 Eh) in
! d-aug-cc-pVDZ with 10 functions dropped (terco's OLD, un-normalised
! dropping decision) went from -74.1475 to an unconverged -67.68 Eh that
! way.
!
! `scf_context_orthogonalize` now decides dropping on the overlap
! normalised to a unit diagonal (Sn = D S D, D = diag(1/sqrt(S_ii))) --
! pyscf's `scf.addons.canonical_orth_`/ORCA's own convention -- rather
! than on S's own eigenvalues, so the same 0.2 threshold now drops 5
! functions, not 10 (terco's own AO norms put every spherical d at
! S_ii = 5/(4*pi), p and s at 1, so the un-normalised and normalised
! eigenvalues differ). Independently confirmed against pyscf 2.14.0,
! basis from generate.py's bse_json_basis (a plain literal, not a .ref --
! like the old reference, nothing here is checked by check_staleness.py):
!
!   S = mol.intor("int1e_ovlp"); d = 1/sqrt(diag(S)); Sn = outer(d,d)*S
!   eigh(Sn) has 1 eigenvalue below 0.1 and 5 below 0.2, exactly matching
!   terco's own n_lindep at those thresholds.
!
!   mf = scf.RHF(mol); mf.level_shift = 0.5
!   mf = scf.addons.remove_linear_dep_(mf, threshold=0.2, lindep=0.9)
!   mf.kernel()  # does not converge by plain DIIS at this gap
!   mf.newton().kernel()  -> -74.14746033894737, converged
!
! (0.1 is NOT used here: pyscf's own ground state there is
! -74.445149289825, but terco's DIIS lands on a different, higher
! stationary point at that threshold -- 0.2 is where terco's own
! converged state IS pyscf's ground state.) The literal basis
! (prune=.false.): pyscf reads the JSON's columns as they are, and
! pruning changes the individual functions, so the dropped space.
!
! A test-drive suite.
!
module check_scf_lindep_gap_suite
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_settings, only: scf_settings_homo_lumo_gap
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_scf_lindep_gap

   real(dp), parameter :: E_REF = -74.147460338947_dp
   !! pyscf, see the header.
   real(dp), parameter :: E_TOL = 1.0e-9_dp

contains

   subroutine collect_scf_lindep_gap(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("O2- with 5 functions dropped matches pyscf", test_dianion_dropped), &
                  new_unittest("gap: closed shell", test_gap_closed), &
                  new_unittest("gap: per spin, not across spins", test_gap_per_spin), &
                  new_unittest("gap: dropped null states are not a LUMO", test_gap_null_states), &
                  new_unittest("gap: undefined without a HOMO and a LUMO", test_gap_undefined) &
                  ]
   end subroutine collect_scf_lindep_gap

   subroutine test_dianion_dropped(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer :: z(1)
      real(dp) :: r(3, 1)

      z = 8
      r = 0.0_dp
      call trc_basis_from_json("basis_sets/d-aug-cc-pvdz.json", 1, z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      opts%lindep_thresh = 0.2_dp
      opts%conv_energy = 1.0e-12_dp
      opts%conv_density = 1.0e-10_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100
      call trc_scf_run(bas, 5, 5, opts, res)   ! O2-: 10 electrons, RHF
      call bas%release()
      call check(error, res%n_lindep == 5, "expected 5 functions dropped below 0.2 (unit-diagonal Sn)")
      if (allocated(error)) return
      call check(error, res%converged, "the SCF did not converge")
      if (allocated(error)) return
      call check(error, abs(res%energy - E_REF) < E_TOL, "energy differs from pyscf in the same subspace")
   end subroutine test_dianion_dropped

   subroutine test_gap_closed(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: eps(4, 1)
      eps(:, 1) = [-1.0_dp, -0.5_dp, 0.25_dp, 1.0_dp]
      call check(error, scf_settings_homo_lumo_gap(4, 1, 4, [2, 0], eps), 0.75_dp, thr=1.0e-15_dp)
   end subroutine test_gap_closed

   subroutine test_gap_per_spin(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: eps(4, 2)
      ! A radical: beta's LUMO (-0.6) is below alpha's HOMO (-0.4). Per spin
      ! the gaps are 0.5 and 0.3; across spins it would be -0.2.
      eps(:, 1) = [-1.0_dp, -0.4_dp, 0.1_dp, 1.0_dp]
      eps(:, 2) = [-0.9_dp, -0.6_dp, 0.2_dp, 1.0_dp]
      call check(error, scf_settings_homo_lumo_gap(4, 2, 4, [2, 1], eps), 0.3_dp, thr=1.0e-15_dp)
   end subroutine test_gap_per_spin

   subroutine test_gap_null_states(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: eps(4, 1)
      ! Two real levels, both occupied; the last two are lifted null states.
      eps(:, 1) = [-1.0_dp, -0.5_dp, 1.0e3_dp, 1.0e3_dp]
      call check(error, scf_settings_homo_lumo_gap(4, 1, 2, [2, 0], eps) == huge(1.0_dp), &
                 "a null state was taken as the LUMO")
   end subroutine test_gap_null_states

   subroutine test_gap_undefined(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: eps(2, 1)
      eps(:, 1) = [-1.0_dp, 0.5_dp]
      call check(error, scf_settings_homo_lumo_gap(2, 1, 2, [0, 0], eps) == huge(1.0_dp), &
                 "no occupied level, yet a gap")
   end subroutine test_gap_undefined

end module check_scf_lindep_gap_suite

program check_scf_lindep_gap
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_scf_lindep_gap_suite, only: collect_scf_lindep_gap
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_lindep_gap, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_scf_lindep_gap
