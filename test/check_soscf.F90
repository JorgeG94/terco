!
! Second-order SCF (`--soscf`, trc_scf_soscf.F90): a test-drive suite that
! fails without the feature.
!
! 1. Water/6-31G*, RHF: DIIS alone and --soscf (DIIS, then a Newton
!    handover) must land on the SAME energy, to 1e-9, which also matches
!    pyscf independently -- `test/pyscf_refs/soscf.ref`'s own E_WATER
!    (`test/pyscf_refs/generate.py CASES["soscf"]` -> `gen_soscf`,
!    read here with `pyscf_ref_io`, never a literal): pyscf 2.14.0,
!    RHF/6-31G* (cart=True), geometry test/water.xyz, basis
!    test/basis_sets/6-31gs.json (bse_json_basis), conv_tol 1e-12,
!    E = -76.0105049953126297 Eh.
!
! 2. O2- (an atomic dianion: a single O, 10 electrons restricted), aug2-
!    cc-pVDZ (`d-aug-cc-pvdz.json`, `prune=.false.` -- the literal basis,
!    matching check_scf_lindep_gap.F90's own reasoning), `lindep_thresh =
!    0.1`: 10 Eh above its HOMO the dropped-function null shift keeps
!    sound (check_scf_lindep_gap.F90's own test, at lindep 0.2, confirms
!    that part independently), but at the LOOSER 0.1 threshold plain DIIS
!    on this run lands on a HIGHER stationary point than the true minimum
!    -- it does not fail to converge, it converges somewhere wrong. pyscf
!    finds the true minimum from the same starting point by exactly the
!    method this module ports -- `soscf.ref`'s own E_O2_ANION, from the
!    SAME generator: `scf.addons.remove_linear_dep_(scf.RHF(mol),
!    threshold=0.1, lindep=0.9)` (pyscf's own, always unit-diagonal-
!    normalised dropping decision -- PR #172/fix/lindep-normalised-
!    overlap brought `scf_context_orthogonalize`'s own decision in line
!    with it, which is why this reference is now generated fresh against
!    that fix rather than carried as a pre-fix literal), `.newton()`:
!    E = -74.445149289824... Eh, converged (|g| 3.3e-7 at this subspace's
!    own AH-solver floor, independently confirmed against
!    check_scf_lindep_gap.F90's own threshold-0.2 case; the unit-diagonal
!    and raw-S dropping decisions happen to coincide in COUNT at
!    threshold 0.1 for this basis -- 1 mode either way -- which is why
!    this number did not move across the fix, only got re-verified
!    against it). Plain mf.kernel() (pyscf's own DIIS) on the SAME
!    subspace does not reach it either -- the stalled, higher stationary
!    point this module's own `--soscf` run has to escape from.
!
!    So: DIIS alone must NOT reach the pyscf minimum (it would not be
!    testing anything otherwise); `--soscf` must.
!
! 3. An open shell (unrestricted) with `--soscf` is refused before a
!    single iteration runs -- never run silently as plain DIIS. mqc's own
!    backend (`mqc_czt_rhf.f90`) has no unrestricted second-order SCF, so
!    this port does not invent one; `trc_scf_settings%resolve` refuses it.
!
! 4. The returned state after --soscf IS the converged state: Tr(D S) =
!    nelec, D = 2 C_occ C_occ^T from the returned cmo, a fresh Fock build
!    from the returned D reproduces res%energy and has a small commutator,
!    eps is diag(C^T F C) for the returned C, and the analytic gradient
!    after --soscf agrees with the one after plain DIIS on the same
!    minimum (water) -- all of this failed before the fix in trc_scf_run's
!    own tail (`scf_result%dmat = scf_state%dold` overwrote the SOSCF
!    density with the pre-handover one; `scf_result%eps` was written on
!    the host only, so a GPU build's own trailing `!$acc update
!    self(...eps)` pulled the stale device copy straight back over it).
!
! 5. Backtracking rescales the step on every rejection (TRC_SOSCF_TRUST0
!    forces an oversized starting trust radius, which used to retry the
!    SAME unscaled step for every one of MAX_BACKTRACKS tries until trust
!    happened to fall below the step's own size).
!
! 6. `scf_soscf_rotation_matrix`'s `exp(kappa)` stays orthogonal to
!    machine precision at sizes the SCF path itself does not happen to
!    exercise (a random antisymmetric kappa at 10 occupied x 200
!    virtual, and a couple of other shapes).
!
! 7. A Krylov/level-shift breakdown (TRC_SOSCF_FORCE_BREAKDOWN, an
!    injected failure -- a real one needs a pathological Hessian this
!    suite has no natural case for) stops the phase immediately, not
!    after `max_iteration` more iterations on a step that was never
!    solved for.
!
module check_soscf_suite
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_integrals, only: trc_1e
   use trc_scf_soscf, only: scf_soscf_rotation_matrix, scf_soscf_trust_after_rejection
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real
   implicit none
   private
   public :: collect_soscf

contains

   subroutine test_setenv(name, value)
      !! Standard bind(c) interop with the C library's own setenv(3),
      !! not a GNU/vendor extension subroutine: trc_scf_soscf.F90's
      !! experiment hooks are read with get_environment_variable, so
      !! setting one from a test has to reach the SAME process
      !! environment a C call would set, portably across gfortran and
      !! nvfortran.
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
      character(len=*), intent(in) :: name, value
      interface
         function c_setenv(name, value, overwrite) bind(c, name="setenv") result(stat)
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: name(*), value(*)
            integer(c_int), value, intent(in) :: overwrite
            integer(c_int) :: stat
         end function c_setenv
      end interface
      integer(c_int) :: stat
      stat = c_setenv(trim(name)//c_null_char, trim(value)//c_null_char, 1_c_int)
   end subroutine test_setenv

   subroutine test_unsetenv(name)
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
      character(len=*), intent(in) :: name
      interface
         function c_unsetenv(name) bind(c, name="unsetenv") result(stat)
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int) :: stat
         end function c_unsetenv
      end interface
      integer(c_int) :: stat
      stat = c_unsetenv(trim(name)//c_null_char)
   end subroutine test_unsetenv

   subroutine collect_soscf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water: DIIS and --soscf land on the same energy, and pyscf's", &
                               test_water_diis_and_soscf_agree), &
                  new_unittest("O2-: DIIS stalls above the true minimum, --soscf reaches it", &
                               test_o2_anion_escapes_the_diis_stationary_point), &
                  new_unittest("--soscf refuses an open shell", test_open_shell_refused), &
                  new_unittest("the returned state after --soscf is the converged one", &
                               test_soscf_returned_state_is_consistent), &
                  new_unittest("a rejection actually shrinks the next trial's step", &
                               test_soscf_trust_shrinks_on_rejection), &
                  new_unittest("--soscf converges from an oversized starting trust radius", &
                               test_soscf_recovers_from_oversized_trust), &
                  new_unittest("exp(kappa) stays orthogonal at larger sizes", &
                               test_rotation_orthogonal_at_various_sizes), &
                  new_unittest("exp(kappa) matches a closed form on a structured kappa", &
                               test_rotation_structured_kappa_matches_closed_form), &
                  new_unittest("a Krylov breakdown stops the phase immediately", &
                               test_soscf_breakdown_stops_cleanly) &
                  ]
   end subroutine collect_soscf

   subroutine test_water_diis_and_soscf_agree(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res_diis, res_soscf
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm, nocc, u
      real(dp) :: e_water_pyscf

      call ref_require("soscf.ref", u)
      call ref_read_real(u, "E_WATER", e_water_pyscf)
      call ref_close(u)

      call read_xyz("water.xyz", natm, z, at_r)
      call trc_basis_from_json("basis_sets/6-31gs.json", natm, z, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      nocc = nint(sum(real(z, dp)))/2

      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 100
      call trc_scf_run(bas, nocc, nocc, opts, res_diis)
      call check(error, res_diis%converged, "DIIS: "//trim(res_diis%message))
      if (allocated(error)) then
         call bas%release()
         return
      end if
      call check(error, abs(res_diis%energy - e_water_pyscf) < 1.0e-8_dp, &
                 "DIIS energy differs from pyscf")
      if (allocated(error)) then
         call bas%release()
         return
      end if

      opts%soscf = .true.
      call trc_scf_run(bas, nocc, nocc, opts, res_soscf)
      call bas%release()
      call check(error, res_soscf%converged, "--soscf: "//trim(res_soscf%message))
      if (allocated(error)) return
      call check(error, abs(res_soscf%energy - e_water_pyscf) < 1.0e-8_dp, &
                 "--soscf energy differs from pyscf")
      if (allocated(error)) return
      call check(error, abs(res_soscf%energy - res_diis%energy) < 1.0e-9_dp, &
                 "--soscf and plain DIIS did not land on the same energy")
   end subroutine test_water_diis_and_soscf_agree

   subroutine test_o2_anion_escapes_the_diis_stationary_point(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res_diis, res_soscf
      integer :: z(1), u
      real(dp) :: r(3, 1)
      real(dp) :: e_water_unused, e_o2_anion_pyscf_newton

      call ref_require("soscf.ref", u)
      call ref_read_real(u, "E_WATER", e_water_unused)   ! soscf.ref's tags are in file order
      call ref_read_real(u, "E_O2_ANION", e_o2_anion_pyscf_newton)
      call ref_close(u)

      z = 8
      r = 0.0_dp
      call trc_basis_from_json("basis_sets/d-aug-cc-pvdz.json", 1, z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()

      opts%lindep_thresh = 0.1_dp
      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100
      call trc_scf_run(bas, 5, 5, opts, res_diis)   ! O2-: 10 electrons, RHF
      ! Plain DIIS either does not converge, or converges to a stationary
      ! point above the true minimum -- either way it must NOT be the
      ! energy pyscf's Newton solver finds from the same subspace, or this
      ! case is not exercising what it is here to exercise.
      call check(error, res_diis%energy > e_o2_anion_pyscf_newton + 1.0e-3_dp, &
                 "plain DIIS unexpectedly reached the true minimum; pick a harder case")
      if (allocated(error)) then
         call bas%release()
         return
      end if

      opts%soscf = .true.
      call trc_scf_run(bas, 5, 5, opts, res_soscf)
      call bas%release()
      call check(error, res_soscf%converged, "--soscf: "//trim(res_soscf%message))
      if (allocated(error)) return
      call check(error, abs(res_soscf%energy - e_o2_anion_pyscf_newton) < 1.0e-7_dp, &
                 "--soscf did not reach pyscf's Newton minimum")
   end subroutine test_o2_anion_escapes_the_diis_stationary_point

   subroutine test_open_shell_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer :: z(1)
      real(dp) :: r(3, 1)

      z = 8
      r = 0.0_dp
      call trc_basis_from_json("basis_sets/sto-3g.json", 1, z, r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()

      opts%soscf = .true.
      opts%max_iteration = 50
      call trc_scf_run(bas, 5, 3, opts, res)   ! an open shell: nalpha /= nbeta
      call bas%release()
      call check(error, .not. res%converged, "an open shell with --soscf should be refused, not run")
      if (allocated(error)) return
      call check(error, res%iterations == 0, "--soscf's refusal should happen before any iteration runs")
      if (allocated(error)) return
      call check(error, index(res%message, "restricted") > 0, &
                 "the refusal message should say --soscf is restricted-only: "//trim(res%message))
   end subroutine test_open_shell_refused

   subroutine test_soscf_returned_state_is_consistent(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res_diis, res_soscf
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm, nao, nocc, k
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :), hcore(:, :)
      real(dp), allocatable :: gmat(:, :), fmat(:, :), dmat_check(:, :)
      real(dp), allocatable :: cf(:, :), cfc(:, :), tmp(:, :), comm_mat(:, :)
      real(dp), allocatable :: grad_diis(:, :), grad_soscf(:, :)
      real(dp) :: nelec, trds, e_fresh, comm_norm, eps_worst

      call read_xyz("water.xyz", natm, z, at_r)
      call trc_basis_from_json("basis_sets/6-31gs.json", natm, z, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      nocc = nint(sum(real(z, dp)))/2
      nao = bas%nao

      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 100
      call trc_scf_run(bas, nocc, nocc, opts, res_diis)
      opts%soscf = .true.
      call trc_scf_run(bas, nocc, nocc, opts, res_soscf)
      call check(error, res_diis%converged .and. res_soscf%converged, &
                 "setup: DIIS or --soscf did not converge")

      ! An independent S/H/eri on the SAME (Cartesian, force_cart) basis,
      ! not anything trc_scf_run built internally.
      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()
      call eri%build(bas, 1.0e-12_dp)
      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao), hcore(nao, nao))
      call trc_1e(bas, pl, smat, tmat, vmat)
      hcore = tmat + vmat

      ! (a) Tr(D S) == nelec
      nelec = real(2*nocc, dp)
      trds = sum(res_soscf%dmat(:, :, 1)*smat)
      if (.not. allocated(error)) call check(error, abs(trds - nelec) < 1.0e-8_dp, &
                                             "Tr(D S) != nelec after --soscf")

      ! (b) D == 2 C_occ C_occ^T from the returned cmo
      allocate (dmat_check(nao, nao))
      dmat_check = 2.0_dp*matmul(res_soscf%cmo(:, 1:nocc, 1), transpose(res_soscf%cmo(:, 1:nocc, 1)))
      if (.not. allocated(error)) call check(error, maxval(abs(dmat_check - res_soscf%dmat(:, :, 1))) < 1.0e-8_dp, &
                                             "returned dmat is not 2 C_occ C_occ^T after --soscf")

      ! (c) a fresh Fock build from the returned D reproduces res%energy,
      ! and its commutator FDS - SDF is small
      allocate (gmat(nao, nao), fmat(nao, nao))
      call eri%fock(bas, res_soscf%dmat(:, :, 1), gmat)
      fmat = hcore + gmat
      e_fresh = sum(res_soscf%dmat(:, :, 1)*hcore) + 0.5_dp*sum(res_soscf%dmat(:, :, 1)*gmat) + res_soscf%e_nuc
      if (.not. allocated(error)) call check(error, abs(e_fresh - res_soscf%energy) < 1.0e-10_dp, &
                                             "a fresh Fock build from the returned D does not reproduce res%energy")

      allocate (tmp(nao, nao), comm_mat(nao, nao))
      tmp = matmul(fmat, matmul(res_soscf%dmat(:, :, 1), smat))
      comm_mat = tmp - transpose(tmp)
      comm_norm = sqrt(sum(comm_mat**2))
      if (.not. allocated(error)) call check(error, comm_norm < 1.0e-6_dp, &
                                             "the commutator at the returned state is not small after --soscf")

      ! (d) eps == diag(C^T F C) for the returned C
      allocate (cf(nao, nao), cfc(nao, nao))
      cf = matmul(fmat, res_soscf%cmo(:, :, 1))
      cfc = matmul(transpose(res_soscf%cmo(:, :, 1)), cf)
      eps_worst = 0.0_dp
      do k = 1, nao
         eps_worst = max(eps_worst, abs(cfc(k, k) - res_soscf%eps(k, 1)))
      end do
      if (.not. allocated(error)) call check(error, eps_worst < 1.0e-6_dp, &
                                             "eps is not diag(C^T F C) for the returned C after --soscf")

      ! (e) the analytic gradient after --soscf matches the one after
      ! plain DIIS: both converge to the SAME minimum (checked above), so
      ! their gradients there must agree. A stale eps on the GPU build
      ! (the energy-weighted density's own orbital-energy weights) would
      ! show up here even where the energy already happened to be right.
      allocate (grad_diis(3, natm), grad_soscf(3, natm))
      call lp%init(nao)
      call trc_scf_gradient(bas, lp, res_diis%nspin, nocc, nocc, res_diis%cmo, res_diis%eps, 1.0e-12_dp, grad_diis)
      call trc_scf_gradient(bas, lp, res_soscf%nspin, nocc, nocc, res_soscf%cmo, res_soscf%eps, 1.0e-12_dp, &
                            grad_soscf)
      if (.not. allocated(error)) call check(error, maxval(abs(grad_diis - grad_soscf)) < 1.0e-8_dp, &
                                             "the gradient after --soscf differs from plain DIIS at the same minimum")
      call lp%release()

      call eri%release()
      call pl%release()
      call bas%release()
   end subroutine test_soscf_returned_state_is_consistent

   subroutine test_soscf_trust_shrinks_on_rejection(error)
      !! A deterministic test of `scf_soscf_trust_after_rejection`'s own
      !! arithmetic, not a full SCF run that may or may not happen to
      !! reject a step on whatever system it has in hand: a step whose
      !! own natural size (`scaling`) is far below the starting trust
      !! radius -- exactly the regime where the OLD rule (`trust =
      !! 0.5*trust`, with nothing capping the step until trust dropped
      !! BELOW scaling) left the next trial's ACTUAL step unchanged,
      !! since `min(1, trust/scaling)` stays 1 regardless of how many
      !! times a trust that is still above `scaling` gets halved.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: scaling, trust0, trust1, trust2, trust3
      real(dp) :: factor0, factor1_old, factor1_new

      scaling = 0.05_dp
      trust0 = 10.0_dp

      ! The bug, demonstrated directly: the OLD rule's own capping
      ! factor on trial 1 (trust = trust0) and what it would still be on
      ! trial 2 (trust = 0.5*trust0, the old rule's only update) are
      ! identical -- the step actually taken would not have changed at
      ! all.
      factor0 = min(1.0_dp, trust0/scaling)
      factor1_old = min(1.0_dp, (0.5_dp*trust0)/scaling)
      call check(error, abs(factor0 - factor1_old) < 1.0e-14_dp, &
                 "setup: this scaling/trust0 pair should reproduce the old rule's own stuck step")
      if (allocated(error)) return

      ! The fix: trust1 is capped to the step's own size on the FIRST
      ! rejection, so trial 2's capping factor already differs from
      ! trial 1's.
      trust1 = scf_soscf_trust_after_rejection(trust0, scaling)
      factor1_new = min(1.0_dp, trust1/scaling)
      call check(error, trust1 <= scaling, &
                 "trust should be capped to at most the step's own size after one rejection")
      if (allocated(error)) return
      call check(error, abs(factor0 - factor1_new) > 1.0e-6_dp, &
                 "the step-capping factor should differ from trial 1 to trial 2")
      if (allocated(error)) return

      ! And it keeps shrinking on further rejections, strictly.
      trust2 = scf_soscf_trust_after_rejection(trust1, scaling)
      trust3 = scf_soscf_trust_after_rejection(trust2, scaling)
      call check(error, trust2 < trust1 .and. trust3 < trust2, &
                 "trust should keep shrinking, strictly, on repeated rejections")
   end subroutine test_soscf_trust_shrinks_on_rejection

   subroutine test_soscf_recovers_from_oversized_trust(error)
      !! TRC_SOSCF_TRUST0 forces the phase to start from a trust radius
      !! deliberately larger than a real Newton step's own size. Before
      !! the fix, the step/radius comparison gated scaling off entirely
      !! while trust still exceeded the (fixed) step's size, so a
      !! rejection only halved an irrelevant number and the SAME unscaled
      !! step was retried, unchanged, for every one of MAX_BACKTRACKS
      !! tries -- the run then reported "found no step downhill" instead
      !! of converging. The fix snaps trust to at most the step's own
      !! size on the first rejection, so it binds immediately.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm, nocc, u
      real(dp) :: e_water_pyscf

      call ref_require("soscf.ref", u)
      call ref_read_real(u, "E_WATER", e_water_pyscf)
      call ref_close(u)

      call read_xyz("water.xyz", natm, z, at_r)
      call trc_basis_from_json("basis_sets/6-31gs.json", natm, z, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      nocc = nint(sum(real(z, dp)))/2

      call test_setenv("TRC_SOSCF_TRUST0", "10.0")
      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 100
      opts%soscf = .true.
      call trc_scf_run(bas, nocc, nocc, opts, res)
      call bas%release()
      call test_unsetenv("TRC_SOSCF_TRUST0")
      call check(error, res%converged, "--soscf with an oversized starting trust radius: "//trim(res%message))
      if (allocated(error)) return
      call check(error, abs(res%energy - e_water_pyscf) < 1.0e-8_dp, &
                 "--soscf with an oversized starting trust radius did not reach the same minimum")
   end subroutine test_soscf_recovers_from_oversized_trust

   subroutine test_rotation_orthogonal_at_various_sizes(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_linalg_t) :: lp
      integer :: sizes(3), isz, n
      character(len=64) :: msg

      sizes = [10 + 200, 5 + 5, 50 + 50]   ! (nocc, nvirt) pairs collapsed to n_mo
      do isz = 1, 3
         n = sizes(isz)
         call lp%init(n)
         write (msg, "(a,i0)") "exp(kappa) orthogonality at n_mo = ", n
         call check_rotation_orthogonal(error, lp, n, msg)
         call lp%release()
         if (allocated(error)) return
      end do
   end subroutine test_rotation_orthogonal_at_various_sizes

   subroutine check_rotation_orthogonal(error, lp, n, msg)
      !! A random antisymmetric kappa of realistic magnitude (max element
      !! ~0.15 rad, comparable to MAX_ROTATION): ||U^T U - I||_max < 1e-12,
      !! and, applied to an S-orthonormal C (S diagonal, C = S^-1/2,
      !! trivially S-orthonormal by construction), C'^T S C' = I and
      !! Tr(D' S) = nelec for D' = 2 C'(:,1:nocc) C'(:,1:nocc)^T.
      type(error_type), allocatable, intent(inout) :: error
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: n
      character(len=*), intent(in) :: msg
      real(dp), allocatable :: kappa(:, :), rotation(:, :), scaled(:, :), term(:, :), next_term(:, :)
      real(dp), allocatable :: utu(:, :), identity(:, :), s(:, :), c(:, :), cprime(:, :), sc(:, :)
      real(dp), allocatable :: ctsc(:, :), dprime(:, :)
      real(dp) :: norm, worst, trds
      integer :: i, j, nocc

      nocc = n/2
      allocate (kappa(n, n), rotation(n, n), scaled(n, n), term(n, n), next_term(n, n))
      allocate (utu(n, n), identity(n, n), s(n, n), c(n, n), cprime(n, n), sc(n, n))
      allocate (ctsc(n, n), dprime(n, n))

      call random_number(kappa)
      kappa = 0.3_dp*(kappa - 0.5_dp)   ! +/- 0.15 rad per element
      do j = 1, n
         do i = 1, j - 1
            kappa(j, i) = -kappa(i, j)
         end do
         kappa(j, j) = 0.0_dp
      end do
      ! The proper 1-norm bound (max absolute COLUMN sum) this test's own
      ! fix uses. The old caller used the largest single element instead,
      ! which is not a matrix norm; a random antisymmetric kappa does not
      ! reliably demonstrate the old convention failing here, because its
      ! spectral norm (what the Taylor series actually sees) grows only as
      ! sqrt(n) by random-matrix statistics, while the 1-norm bound below
      ! grows as n -- the two differ by about a factor of n for a random
      ! matrix, but the entrywise bound still happens to dominate the
      ! spectral norm at the sizes tested. The mathematical argument
      ! stands regardless (the largest element is provably not a bound on
      ! any operator norm, so scaling-and-squaring from it is unsound);
      ! this test only checks the FIXED routine's own correctness.
      norm = 0.0_dp
      do j = 1, n
         norm = max(norm, sum(abs(kappa(:, j))))
      end do

      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
      s = 0.0_dp
      c = 0.0_dp
      do i = 1, n
         s(i, i) = 0.5_dp + real(i, dp)/real(n, dp)   ! a positive diagonal "overlap"
         c(i, i) = 1.0_dp/sqrt(s(i, i))               ! C = S^-1/2: C^T S C = I by construction
      end do

      !$acc enter data copyin(kappa, s, c) create(rotation, scaled, term, next_term, utu, cprime, sc, ctsc, dprime)
      call scf_soscf_rotation_matrix(lp, n, norm, kappa, rotation, scaled, term, next_term)
      call lp%gemm_mat(rotation, rotation, utu, transa="T")
      call lp%gemm_mat(c, rotation, cprime)
      call lp%gemm_mat(s, cprime, sc)
      call lp%gemm_mat(cprime, sc, ctsc, transa="T")
      call lp%density(cprime, nocc, 2.0_dp, dprime)
      !$acc update self(utu, ctsc, dprime)

      worst = maxval(abs(utu - identity))
      call check(error, worst < 1.0e-12_dp, trim(msg)//": ||U^T U - I||_max too large")
      if (.not. allocated(error)) then
         worst = maxval(abs(ctsc - identity))
         call check(error, worst < 1.0e-10_dp, trim(msg)//": C'^T S C' is not the identity")
      end if
      if (.not. allocated(error)) then
         trds = sum(dprime*s)
         call check(error, abs(trds - real(2*nocc, dp)) < 1.0e-8_dp, trim(msg)//": Tr(D' S) != nelec")
      end if

      !$acc exit data delete(kappa, s, c, rotation, scaled, term, next_term, utu, cprime, sc, ctsc, dprime)
   end subroutine check_rotation_orthogonal

   subroutine test_rotation_structured_kappa_matches_closed_form(error)
      !! check_rotation_orthogonal's RANDOM kappa cannot demonstrate the
      !! largest-single-element norm actually breaking orthogonality (this
      !! module's own header note: a random antisymmetric matrix's
      !! spectral norm grows only as sqrt(n), not n, so the entrywise
      !! bound still dominates it at the sizes tested). A STRUCTURED kappa
      !! can: kappa(nocc+1:n, 1:nocc) = c, kappa(1:nocc, nocc+1:n) = -c
      !! (10 occ x 200 virt) is rank 2, c*sqrt(nocc*nvirt)*(v u^T - u v^T)
      !! for the normalised ones-vectors u (occ block) and v (virt block),
      !! so it has a CLOSED FORM: a rotation by theta = c*sqrt(nocc*nvirt)
      !! in the 2D span of u and v, identity elsewhere --
      !!   exp(kappa) = I + (cos(theta)-1)(u u^T + v v^T) + sin(theta)(v u^T - u v^T).
      !! Its 1-norm is nvirt*c (every occ column sums nvirt entries of c);
      !! its largest-element norm is c -- a factor of nvirt = 200 apart,
      !! big enough that the largest-element scale picks ZERO squarings
      !! for every c tried here (all < 0.5) and runs the Taylor series
      !! directly on the unscaled, large-angle matrix. Checked against the
      !! closed form, not just self-consistency, so this also catches a
      !! bug in the Taylor/squaring implementation itself, not only in
      !! which norm a caller passes it.
      type(error_type), allocatable, intent(out) :: error
      type(trc_linalg_t) :: lp
      integer, parameter :: nocc = 10, nvirt = 200, n = nocc + nvirt
      real(dp), allocatable :: kappa(:, :), rotation(:, :), rotation_bad(:, :)
      real(dp), allocatable :: scaled(:, :), term(:, :), next_term(:, :), u_exact(:, :)
      real(dp) :: c, norm_good, norm_bad, theta, worst, worst_bad, root_nocc_nvirt
      integer :: ic, i, j
      real(dp), parameter :: cvals(3) = [0.05_dp, 0.15_dp, 0.30_dp]
      character(len=64) :: msg

      allocate (kappa(n, n), rotation(n, n), rotation_bad(n, n))
      allocate (scaled(n, n), term(n, n), next_term(n, n), u_exact(n, n))
      call lp%init(n)
      !$acc enter data create(kappa, rotation, rotation_bad, scaled, term, next_term)

      do ic = 1, 3
         c = cvals(ic)
         write (msg, "(a,f4.2)") "exp(kappa), structured 10x200 occ-virt block, c = ", c

         kappa = 0.0_dp
         kappa(nocc + 1:n, 1:nocc) = c
         kappa(1:nocc, nocc + 1:n) = -c
         root_nocc_nvirt = sqrt(real(nocc, dp)*real(nvirt, dp))
         theta = c*root_nocc_nvirt
         norm_good = real(nvirt, dp)*c   ! max absolute column sum: nvirt entries of c per occ column
         norm_bad = c                    ! the OLD convention: the largest single element

         !$acc update device(kappa)
         call scf_soscf_rotation_matrix(lp, n, norm_good, kappa, rotation, scaled, term, next_term)
         call scf_soscf_rotation_matrix(lp, n, norm_bad, kappa, rotation_bad, scaled, term, next_term)
         !$acc update self(rotation, rotation_bad)

         u_exact = 0.0_dp
         do i = 1, n
            u_exact(i, i) = 1.0_dp
         end do
         do j = 1, nocc
            do i = 1, nocc
               u_exact(i, j) = u_exact(i, j) + (cos(theta) - 1.0_dp)/real(nocc, dp)
            end do
         end do
         do j = nocc + 1, n
            do i = nocc + 1, n
               u_exact(i, j) = u_exact(i, j) + (cos(theta) - 1.0_dp)/real(nvirt, dp)
            end do
         end do
         do j = 1, nocc
            do i = nocc + 1, n
               u_exact(i, j) = u_exact(i, j) + sin(theta)/root_nocc_nvirt
            end do
         end do
         do j = nocc + 1, n
            do i = 1, nocc
               u_exact(i, j) = u_exact(i, j) - sin(theta)/root_nocc_nvirt
            end do
         end do

         worst = maxval(abs(rotation - u_exact))
         call check(error, worst < 1.0e-12_dp, trim(msg)//": |U - U_exact|_max too large with the fixed norm")
         if (allocated(error)) exit

         ! At c = 0.05 (theta ~ 2.24 rad), this module's own TAYLOR_TERMS = 18
         ! is accurate enough even on the UNSCALED matrix that the largest-
         ! element norm's zero squarings do not visibly break orthogonality
         ! (measured: |U_bad^T U_bad - I| ~ 6e-12, below the 1e-6 gate) --
         ! 18 terms is more headroom than an order-8 Taylor series would
         ! have, so the bug needs a bigger angle to show through it. It
         ! does at c = 0.15 (theta ~ 6.71) and c = 0.30 (theta ~ 13.42, where
         ! it is not a small discrepancy but outright loss of orthogonality,
         ! |U_bad^T U_bad - I| > 1e7 measured): check those two, not 0.05.
         if (ic > 1) then
            worst_bad = maxval(abs(rotation_bad - u_exact))
            call check(error, worst_bad > 1.0e-6_dp, &
                       trim(msg)//": the largest-element norm should have broken this case but did not "// &
                       "-- this test's own setup is not exercising the bug")
            if (allocated(error)) exit
         end if
      end do

      !$acc exit data delete(kappa, rotation, rotation_bad, scaled, term, next_term)
      call lp%release()
   end subroutine test_rotation_structured_kappa_matches_closed_form

   subroutine test_soscf_breakdown_stops_cleanly(error)
      !! TRC_SOSCF_FORCE_BREAKDOWN injects a Krylov/level-shift breakdown
      !! on the very first Newton step. The phase must stop there --
      !! converged = .false., a clear message, iterations left at the
      !! DIIS handover point -- not keep iterating to max_iteration on a
      !! step that was never actually solved for.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm, nocc

      call read_xyz("water.xyz", natm, z, at_r)
      call trc_basis_from_json("basis_sets/6-31gs.json", natm, z, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      nocc = nint(sum(real(z, dp)))/2

      call test_setenv("TRC_SOSCF_FORCE_BREAKDOWN", "1")
      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 100
      opts%soscf = .true.
      call trc_scf_run(bas, nocc, nocc, opts, res)
      call bas%release()
      call test_unsetenv("TRC_SOSCF_FORCE_BREAKDOWN")
      call check(error, .not. res%converged, "an injected Krylov breakdown should not report convergence")
      if (allocated(error)) return
      call check(error, index(res%message, "broke down") > 0, &
                 "the breakdown should leave a clear message: "//trim(res%message))
      if (allocated(error)) return
      ! Stopped at the DIIS handover point (well under max_iteration =
      ! 100, which this is here to rule out), not one Newton iteration
      ! further: the breakdown is injected before the FIRST Newton step
      ! does anything, so no accepted step ever moves `iterations` past
      ! wherever DIIS handed over (water's own handover is well inside
      ! the first ~20 DIIS iterations).
      call check(error, res%iterations < 20, &
                 "an injected Krylov breakdown should stop at once, not iterate on: "//trim(res%message))
   end subroutine test_soscf_breakdown_stops_cleanly

end module check_soscf_suite

program check_soscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_soscf_suite, only: collect_soscf
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_soscf, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_soscf
