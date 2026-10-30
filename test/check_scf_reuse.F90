!
! Reuse across SCF/gradient calls, added for a sequence of nearby
! geometries at the same molecule/basis (the finite-difference Hessian,
! and later a geometry optimiser): `trc_pair_eri_cache_t`
! (trc_scf_context.F90) lets the SCF's own pair list and (outside --df)
! four-centre engine stay alive for a gradient right after it, instead of
! a second pl%build/eri%build at the SAME geometry and threshold
! (trc_scf_grad.F90's own `pair_eri_cache` argument, which also reports
! whether it actually took that path through its `reused` out argument);
! `scf_context_project_density` (same module) is the orbital-projection warm start,
! trc_capi.F90's TRC_GUESS_PROJECTED.
!
! A test-drive suite. Reuse is checked as a pure PERFORMANCE path, never a
! correctness one: every case below either compares a cached call against
! an independently fresh one (never against a stored literal), or checks
! that a cache that should NOT be trusted is provably not reused (`reused`
! comes back false) rather than merely "happens to give the right answer
! this time". Cases:
!
!   1. at one geometry, a gradient computed WITH the cache and WITHOUT it
!      agrees to round-off, AND the cached call's own `reused` flag comes
!      back true -- without that flag this would equally pass if the
!      reuse path never ran at all (comparing fresh against fresh).
!   2. the SAME cache, handed to a gradient at a DIFFERENT geometry with
!      NO SCF rerun in between (the geometry simply moved, so the
!      fingerprint `trc_pair_eri_cache_t` now stores -- nshell, natm, the
!      atom centres -- no longer matches what it was built with): `reused`
!      comes back false, and the result agrees with an independently
!      fresh build at the NEW geometry, not the old one.
!   3. a threshold mismatch (SAME geometry, cache built at a different
!      `eri_thresh`) is the same story: `reused` false, result matching a
!      fresh build at the NEW threshold.
!   4. `batch_res` mismatch -- the SCF's own engine built at a non-default
!      `opts%batch_res`, while this routine's own fresh build always uses
!      the default -- also forces a rebuild rather than risking a result
!      that differs from the no-cache path by more than round-off.
!   5. `trc_pair_eri_cache_t%release()` is safe called twice.
!   6. `scf_context_project_density` at an UNCHANGED geometry (S' = S) reproduces the
!      SCF's own converged density to round-off -- the identity case
!      every orbital-projection warm start must degenerate to.
!   7. `scf_context_project_density` at a MOVED geometry gives an S'-orthonormal C':
!      checked through D'S'D' = occ_s D' (idempotency under the NEW
!      overlap), the matrix identity that holds iff C'^T S' C' = I for a
!      full-rank C -- and an SCF started from it converges to the SAME
!      energy as one started from SAD.
!   8. the C API: `trc_run_scf` -> a geometry OR a setting change ->
!      `trc_run_gradient` refuses with TRC_ERR_STATE rather than reusing a
!      cache the change already invalidated (`invalidate_runs`).
!   9. `TRC_GUESS_PROJECTED` through `trc_set_guess`: both its own
!      refusals (no previous SCF orbitals; a size/spin mismatch against
!      the previous ones), each with TRC_ERR_STATE.
!   10. the full C-API path `trc_set_geometry` makes possible: trc_run_scf
!       at one geometry, trc_set_geometry to a moved one, TRC_GUESS_PROJECTED,
!       trc_run_scf again -- converges to the same energy as a fresh SAD
!       run at the moved geometry.
!
module check_scf_reuse_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char, c_null_ptr
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run, trc_pair_eri_cache_t
   use trc_scf_grad, only: trc_scf_gradient
   use trc_scf_context, only: scf_context_project_density
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_guess, trc_set_unrestricted, trc_set_screening, trc_run_scf, &
                               trc_run_gradient, trc_set_geometry, trc_energy, trc_converged, &
                               TRC_OK, TRC_ERR_STATE, TRC_GUESS_GWH, TRC_GUESS_SAD, TRC_GUESS_PROJECTED
   implicit none
   private
   public :: collect_scf_reuse

   real(dp), parameter :: THRESH = 1.0e-12_dp
   real(dp), parameter :: THRESH_OTHER = 1.0e-10_dp
   real(dp), parameter :: TOL_GRAD = 1.0e-12_dp
   !! Two builds of the SAME pairs/Schwarz/bins at the SAME threshold,
   !! contracted against the SAME density -- round-off, not truncation.
   real(dp), parameter :: TOL_PROJ = 1.0e-9_dp
   !! `scf_context_project_density` rebuilds S independently (its own pair list/trc_1e
   !! call, not the SCF's), and goes through an eigendecomposition
   !! (`sym_orthog`) the SCF's own density build never does -- looser than
   !! TOL_GRAD, still far tighter than any physical effect this test could
   !! be confusing with a bug.
   real(dp), parameter :: TOL_ORTHO = 1.0e-9_dp
   !! D'S'D' = occ_s D' idempotency, same eigendecomposition as TOL_PROJ.
   real(dp), parameter :: TOL_ENERGY = 1.0e-9_dp
   !! Two different warm starts (projected vs SAD) converging to the SAME
   !! stationary point, both at opts%conv_energy tighter than this.

contains

   subroutine collect_scf_reuse(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("gradient with vs without cache: bit-identical, reused flag true", &
                      test_reuse_bitidentical), &
         new_unittest("a cache at a MOVED geometry (no SCF rerun) is never reused", &
                      test_stale_geometry_rebuilds), &
         new_unittest("a threshold mismatch forces a rebuild", test_thresh_mismatch_rebuilds), &
         new_unittest("a batch_res mismatch blocks engine reuse but stays correct", &
                      test_batch_res_mismatch_blocks_reuse), &
         new_unittest("trc_pair_eri_cache_t%release is safe called twice", test_release_idempotent), &
         new_unittest("projected guess at S'=S reproduces the SCF density", test_projected_identity), &
         new_unittest("projected guess at a moved geometry is S'-orthonormal, same energy as SAD", &
                      test_projected_moved_geometry), &
         new_unittest("C API: a geometry change after trc_run_scf makes trc_run_gradient refuse", &
                      test_capi_geometry_change_refuses), &
         new_unittest("C API: a setting change after trc_run_scf makes trc_run_gradient refuse", &
                      test_capi_setting_change_refuses), &
         new_unittest("C API: TRC_GUESS_PROJECTED refuses with no previous SCF", &
                      test_capi_projected_no_scf), &
         new_unittest("C API: TRC_GUESS_PROJECTED refuses a size/spin mismatch", &
                      test_capi_projected_size_mismatch), &
         new_unittest("C API: trc_run_scf -> trc_set_geometry -> TRC_GUESS_PROJECTED -> trc_run_scf "// &
                      "converges like SAD", test_capi_projected_after_set_geometry) &
         ]
   end subroutine collect_scf_reuse

   subroutine water_geom(nat, z, r0, nalpha, nbeta)
      integer, intent(out) :: nat
      integer, allocatable, intent(out) :: z(:)
      real(dp), allocatable, intent(out) :: r0(:, :)
      integer, intent(out) :: nalpha, nbeta
      call read_xyz("water.xyz", nat, z, r0)
      nalpha = sum(z)/2
      nbeta = nalpha
   end subroutine water_geom

   subroutine test_reuse_bitidentical(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(trc_pair_eri_cache_t), target :: pec
      !! TARGET: see trc_scf_context.F90's own `pair_eri_cache` argument
      !! header -- `scf_ctx_setup` writes `built`/`thresh_built`/`pl`/`eri`
      !! into THIS variable's storage through a pointer, which needs the
      !! actual argument to carry TARGET too.
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      integer, allocatable :: z(:)
      real(dp), allocatable :: r0(:, :)
      real(dp), allocatable :: grad_cached(:, :), grad_fresh(:, :)
      real(dp) :: worst
      integer :: nat, nalpha, nbeta, nao
      logical :: reused_flag

      call water_geom(nat, z, r0, nalpha, nbeta)
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r0, b, err)
      if (err%has_error()) then
         call test_failed(error, "basis load: "//err%get_message()); return
      end if
      call b%to_device()

      opts%max_iteration = 30; opts%conv_energy = 1.0e-11_dp; opts%conv_diis = 1.0e-8_dp
      call trc_scf_run(b, nalpha, nbeta, opts, res, pair_eri_cache=pec)
      if (res%iterations <= 0) then
         call b%release(); call pec%release()
         call test_failed(error, "the SCF never ran: "//trim(res%message)); return
      end if

      nao = b%nao
      allocate (grad_cached(3, nat), grad_fresh(3, nat))
      call lp%init(nao)
      call trc_scf_gradient(b, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad_cached, &
                            pair_eri_cache=pec, reused=reused_flag)
      call trc_scf_gradient(b, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad_fresh)
      call lp%release()
      worst = maxval(abs(grad_cached - grad_fresh))
      call b%release()
      call pec%release()

      call check(error, reused_flag, "the cached call must report that it reused the pair list/engine")
      if (allocated(error)) return
      call check(error, worst, 0.0_dp, thr=TOL_GRAD, &
                 message="cached-vs-fresh gradient worst|diff| is not round-off")
   end subroutine test_reuse_bitidentical

   subroutine test_stale_geometry_rebuilds(error)
      !! One cache, filled by an SCF at the equilibrium geometry, handed
      !! STRAIGHT to a gradient at a DIFFERENT geometry -- one H moved 30
      !! Bohr along z, far enough that its shell pairs with the other two
      !! atoms fall below THRESH's Schwarz bound (water's primitives have
      !! no exponent small enough to keep that overlap alive at this
      !! range), so the surviving pair count actually differs: a naive
      !! reuse here would not merely be "a bit stale", it would index
      !! against the wrong pair list. No SCF reruns at the new geometry --
      !! orbitals from the FIRST one are reused as-is, since this check is
      !! only about whether the CACHE gets reused, not about what the
      !! resulting number means physically.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b0, b1
      type(trc_pair_eri_cache_t), target :: pec
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      integer, allocatable :: z(:)
      real(dp), allocatable :: r0(:, :), r1(:, :)
      real(dp), allocatable :: grad_stale(:, :), grad_fresh(:, :)
      real(dp) :: worst
      integer :: nat, nalpha, nbeta, nao
      logical :: reused_flag

      call water_geom(nat, z, r0, nalpha, nbeta)
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r0, b0, err)
      if (err%has_error()) then
         call test_failed(error, "basis load (equilibrium): "//err%get_message()); return
      end if
      call b0%to_device()
      opts%max_iteration = 30; opts%conv_energy = 1.0e-11_dp; opts%conv_diis = 1.0e-8_dp
      call trc_scf_run(b0, nalpha, nbeta, opts, res, pair_eri_cache=pec)
      if (res%iterations <= 0) then
         call b0%release(); call pec%release()
         call test_failed(error, "the reference SCF never ran: "//trim(res%message)); return
      end if
      call b0%release()

      allocate (r1, source=r0)
      r1(3, 3) = r1(3, 3) + 30.0_dp
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r1, b1, err)
      if (err%has_error()) then
         call pec%release()
         call test_failed(error, "basis load (moved): "//err%get_message()); return
      end if
      call b1%to_device()

      nao = b1%nao
      allocate (grad_stale(3, nat), grad_fresh(3, nat))
      call lp%init(nao)
      call trc_scf_gradient(b1, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad_stale, &
                            pair_eri_cache=pec, reused=reused_flag)
      call trc_scf_gradient(b1, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad_fresh)
      call lp%release()
      worst = maxval(abs(grad_stale - grad_fresh))
      call b1%release()
      call pec%release()

      call check(error, .not. reused_flag, "a cache fingerprinted to a different geometry must not be reused")
      if (allocated(error)) return
      call check(error, worst, 0.0_dp, thr=TOL_GRAD, &
                 message="the stale-cache call does not match an independent fresh build at the NEW geometry")
   end subroutine test_stale_geometry_rebuilds

   subroutine test_thresh_mismatch_rebuilds(error)
      !! SAME geometry, cache built at THRESH, gradient asked for at
      !! THRESH_OTHER: `matches` must gate on the threshold too, not only
      !! the geometry fingerprint.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(trc_pair_eri_cache_t), target :: pec
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      integer, allocatable :: z(:)
      real(dp), allocatable :: r0(:, :)
      real(dp), allocatable :: grad_stale(:, :), grad_fresh(:, :)
      real(dp) :: worst
      integer :: nat, nalpha, nbeta, nao
      logical :: reused_flag

      call water_geom(nat, z, r0, nalpha, nbeta)
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r0, b, err)
      if (err%has_error()) then
         call test_failed(error, "basis load: "//err%get_message()); return
      end if
      call b%to_device()
      opts%max_iteration = 30; opts%conv_energy = 1.0e-11_dp; opts%conv_diis = 1.0e-8_dp
      opts%eri_thresh = THRESH
      call trc_scf_run(b, nalpha, nbeta, opts, res, pair_eri_cache=pec)
      if (res%iterations <= 0) then
         call b%release(); call pec%release()
         call test_failed(error, "the SCF never ran: "//trim(res%message)); return
      end if

      nao = b%nao
      allocate (grad_stale(3, nat), grad_fresh(3, nat))
      call lp%init(nao)
      call trc_scf_gradient(b, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH_OTHER, grad_stale, &
                            pair_eri_cache=pec, reused=reused_flag)
      call trc_scf_gradient(b, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH_OTHER, grad_fresh)
      call lp%release()
      worst = maxval(abs(grad_stale - grad_fresh))
      call b%release()
      call pec%release()

      call check(error, .not. reused_flag, "a cache built at a different threshold must not be reused")
      if (allocated(error)) return
      call check(error, worst, 0.0_dp, thr=TOL_GRAD, &
                 message="the threshold-mismatch call does not match an independent fresh build at the NEW threshold")
   end subroutine test_thresh_mismatch_rebuilds

   subroutine test_batch_res_mismatch_blocks_reuse(error)
      !! The SCF's own engine, built with a non-default `opts%batch_res`
      !! (trc_scf_context.F90's own FINGERPRINT note on `batch_res_built`):
      !! `trc_scf_gradient`'s OWN fresh build always batches at the
      !! default (1), so reusing an engine built at a different value
      !! would make the cached and fresh paths differ by more than
      !! round-off. `reused` must come back false, and the numbers must
      !! still agree -- the correctness this gate exists for.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(trc_pair_eri_cache_t), target :: pec
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      integer, allocatable :: z(:)
      real(dp), allocatable :: r0(:, :)
      real(dp), allocatable :: grad_cached(:, :), grad_fresh(:, :)
      real(dp) :: worst
      integer :: nat, nalpha, nbeta, nao
      logical :: reused_flag

      call water_geom(nat, z, r0, nalpha, nbeta)
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r0, b, err)
      if (err%has_error()) then
         call test_failed(error, "basis load: "//err%get_message()); return
      end if
      call b%to_device()
      opts%max_iteration = 30; opts%conv_energy = 1.0e-11_dp; opts%conv_diis = 1.0e-8_dp
      opts%batch_res = 2   ! non-default: finer Schwarz-bound buckets than 1
      call trc_scf_run(b, nalpha, nbeta, opts, res, pair_eri_cache=pec)
      if (res%iterations <= 0) then
         call b%release(); call pec%release()
         call test_failed(error, "the SCF never ran: "//trim(res%message)); return
      end if

      nao = b%nao
      allocate (grad_cached(3, nat), grad_fresh(3, nat))
      call lp%init(nao)
      call trc_scf_gradient(b, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad_cached, &
                            pair_eri_cache=pec, reused=reused_flag)
      call trc_scf_gradient(b, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad_fresh)
      call lp%release()
      worst = maxval(abs(grad_cached - grad_fresh))
      call b%release()
      call pec%release()

      call check(error, .not. reused_flag, "a batch_res-2 engine must not be reused by the (batch_res-1) fresh path")
      if (allocated(error)) return
      call check(error, worst, 0.0_dp, thr=TOL_GRAD, &
                 message="falling back to a fresh build keeps the result correct at non-default batching")
   end subroutine test_batch_res_mismatch_blocks_reuse

   subroutine test_release_idempotent(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(trc_pair_eri_cache_t), target :: pec
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      integer, allocatable :: z(:)
      real(dp), allocatable :: r0(:, :)
      integer :: nat, nalpha, nbeta

      call water_geom(nat, z, r0, nalpha, nbeta)
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r0, b, err)
      if (err%has_error()) then
         call test_failed(error, "basis load: "//err%get_message()); return
      end if
      call b%to_device()
      opts%max_iteration = 5
      call trc_scf_run(b, nalpha, nbeta, opts, res, pair_eri_cache=pec)
      call b%release()
      call pec%release()
      call pec%release()   ! must not crash or double-free

      call check(error, .not. pec%built, "a released cache reports built = .false.")
   end subroutine test_release_idempotent

   subroutine test_projected_identity(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), h(:, :), v(:, :), dproj(:, :, :)
      real(dp) :: worst
      integer, allocatable :: z(:)
      real(dp), allocatable :: r0(:, :)
      integer :: nat, nalpha, nbeta, nao, nocc(2)

      call water_geom(nat, z, r0, nalpha, nbeta)
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r0, b, err)
      if (err%has_error()) then
         call test_failed(error, "basis load: "//err%get_message()); return
      end if
      call b%to_device()
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iteration = 60
      call trc_scf_run(b, nalpha, nbeta, opts, res)
      if (.not. res%converged) then
         call b%release()
         call test_failed(error, "the reference SCF did not converge"); return
      end if
      nao = b%nao
      allocate (s(nao, nao), h(nao, nao), v(nao, nao), dproj(nao, nao, res%nspin))
      call pl%build(b, THRESH)
      call pl%to_device()
      call trc_1e(b, pl, s, h, v)
      call pl%release()
      call lp%init(nao)
      nocc = [nalpha, nbeta]
      call scf_context_project_density(lp, nao, res%nspin, nocc, res%cmo, s, dproj)
      call lp%release()
      worst = maxval(abs(dproj - res%dmat))
      call b%release()
      deallocate (s, h, v, dproj)

      call check(error, worst, 0.0_dp, thr=TOL_PROJ, &
                 message="projecting at the SAME geometry must reproduce the SCF density")
   end subroutine test_projected_identity

   subroutine test_projected_moved_geometry(error)
      !! A converged SCF at geometry A projected into geometry B's overlap
      !! (a 0.3 Bohr displacement -- comfortably inside the step sizes
      !! test/bench_warmstart.F90 measured, and far from water's
      !! dissociation): C' is S'-orthonormal (checked as D'S'D' = occ_s D',
      !! which holds iff C'^T S' C' = I for a full-rank C -- see
      !! `scf_context_project_density`'s own header), and an SCF started from the
      !! projected density converges to the SAME energy as one started
      !! from SAD.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: ba, bb
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res_a, res_proj, res_sad
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), h(:, :), v(:, :), dproj(:, :, :)
      real(dp), allocatable :: sd(:, :), sds(:, :)
      real(dp) :: worst_ortho, worst_energy, occ_s
      integer, allocatable :: z(:)
      real(dp), allocatable :: r0(:, :), r1(:, :)
      integer :: nat, nalpha, nbeta, nao, nocc(2), s_spin

      call water_geom(nat, z, r0, nalpha, nbeta)
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r0, ba, err)
      if (err%has_error()) then
         call test_failed(error, "basis load (A): "//err%get_message()); return
      end if
      call ba%to_device()
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iteration = 60
      call trc_scf_run(ba, nalpha, nbeta, opts, res_a)
      if (.not. res_a%converged) then
         call ba%release()
         call test_failed(error, "the geometry-A SCF did not converge"); return
      end if
      call ba%release()

      allocate (r1, source=r0)
      r1(1, 1) = r1(1, 1) + 0.3_dp
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r1, bb, err)
      if (err%has_error()) then
         call test_failed(error, "basis load (B): "//err%get_message()); return
      end if
      call bb%to_device()

      nao = bb%nao
      allocate (s(nao, nao), h(nao, nao), v(nao, nao), dproj(nao, nao, res_a%nspin))
      call pl%build(bb, THRESH)
      call pl%to_device()
      call trc_1e(bb, pl, s, h, v)
      call pl%release()
      call lp%init(nao)
      nocc = [nalpha, nbeta]
      call scf_context_project_density(lp, nao, res_a%nspin, nocc, res_a%cmo, s, dproj)
      call lp%release()

      ! D' S' D' = occ_s D', the idempotency that holds iff C'^T S' C' = I
      ! for the (full-rank, by construction) projected C' -- see this
      ! test's own header and `scf_context_project_density`'s.
      allocate (sd(nao, nao), sds(nao, nao))
      occ_s = merge(2.0_dp, 1.0_dp, res_a%nspin == 1)
      worst_ortho = 0.0_dp
      do s_spin = 1, res_a%nspin
         sd = matmul(s, dproj(:, :, s_spin))
         sds = matmul(dproj(:, :, s_spin), sd)
         worst_ortho = max(worst_ortho, maxval(abs(sds - occ_s*dproj(:, :, s_spin))))
      end do
      deallocate (s, h, v, dproj, sd, sds)

      call trc_scf_run(bb, nalpha, nbeta, opts, res_proj, cmo_project=res_a%cmo)
      call trc_scf_run(bb, nalpha, nbeta, opts, res_sad)
      call bb%release()

      if (.not. res_proj%converged .or. .not. res_sad%converged) then
         call test_failed(error, "the geometry-B SCF (projected or SAD start) did not converge"); return
      end if
      worst_energy = abs(res_proj%energy - res_sad%energy)

      call check(error, worst_ortho, 0.0_dp, thr=TOL_ORTHO, &
                 message="the projected C' is not S'-orthonormal (D'S'D' /= occ_s D')")
      if (allocated(error)) return
      call check(error, worst_energy, 0.0_dp, thr=TOL_ENERGY, &
                 message="the projected and SAD warm starts converge to different energies at geometry B")
   end subroutine test_projected_moved_geometry

   subroutine test_capi_geometry_change_refuses(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, rc_grad
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)

      call read_xyz("water.xyz", natm, zint, at_r)
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, "basis_sets/6-31g.json"//c_null_char)
      rc = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the reference SCF through the C API must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      ! Move an atom -- a geometry change, through the SAME setter that
      ! drops `cx%peri` (invalidate_runs, trc_capi.F90).
      at_r(3, 1) = at_r(3, 1) + 0.1_dp
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc_grad = trc_run_gradient(h)
      rc = trc_destroy(h)
      call check(error, rc_grad == TRC_ERR_STATE, &
                 "trc_run_gradient after a geometry change with no rerun must refuse with TRC_ERR_STATE")
   end subroutine test_capi_geometry_change_refuses

   subroutine test_capi_setting_change_refuses(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, rc_grad
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)

      call read_xyz("water.xyz", natm, zint, at_r)
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, "basis_sets/6-31g.json"//c_null_char)
      rc = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the reference SCF through the C API must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      rc = trc_run_gradient(h)
      call check(error, rc == TRC_OK, "the gradient right after the SCF, with no setting change, must succeed")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      ! A setting change -- trc_set_screening calls invalidate_runs too,
      ! same as every other setter that changes something a cached result
      ! depends on.
      rc = trc_set_screening(h, 1.0e-10_c_double)
      rc_grad = trc_run_gradient(h)
      rc = trc_destroy(h)
      call check(error, rc_grad == TRC_ERR_STATE, &
                 "trc_run_gradient after a setting change with no rerun must refuse with TRC_ERR_STATE")
   end subroutine test_capi_setting_change_refuses

   subroutine test_capi_projected_no_scf(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, rc_run
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)

      call read_xyz("water.xyz", natm, zint, at_r)
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, "basis_sets/6-31g.json"//c_null_char)
      rc = trc_set_guess(h, TRC_GUESS_PROJECTED, c_null_ptr, 1_c_int)
      rc_run = trc_run_scf(h)
      rc = trc_destroy(h)
      call check(error, rc_run == TRC_ERR_STATE, &
                 "TRC_GUESS_PROJECTED with no previous SCF must refuse with TRC_ERR_STATE")
   end subroutine test_capi_projected_no_scf

   subroutine test_capi_projected_size_mismatch(error)
      !! A restricted SCF converges and stores `cmo` at nspin 1; forcing
      !! unrestricted (no basis change -- `cmo` is NOT dropped by that
      !! setter, only by the basis ones, trc_capi.F90's `take_basis`) asks
      !! the next run for nspin 2, a shape the stored orbitals do not
      !! have.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, rc_run
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)

      call read_xyz("water.xyz", natm, zint, at_r)
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_basis_json(h, "basis_sets/6-31g.json"//c_null_char)
      rc = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "the reference (restricted) SCF through the C API must converge")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      rc = trc_set_unrestricted(h, 1_c_int)
      rc = trc_set_guess(h, TRC_GUESS_PROJECTED, c_null_ptr, 1_c_int)
      rc_run = trc_run_scf(h)
      rc = trc_destroy(h)
      call check(error, rc_run == TRC_ERR_STATE, &
                 "TRC_GUESS_PROJECTED at a size/spin mismatch must refuse with TRC_ERR_STATE")
   end subroutine test_capi_projected_size_mismatch

   subroutine test_capi_projected_after_set_geometry(error)
      !! The path `trc_set_geometry` exists for: trc_run_scf at geometry A,
      !! trc_set_geometry to geometry B (no basis rebuild), TRC_GUESS_PROJECTED,
      !! trc_run_scf again. `cx%bb%cmo` (trc_run_scf's own last converged
      !! orbitals) survives trc_set_geometry's invalidate_runs -- that
      !! survival is the whole point of TRC_GUESS_PROJECTED's doc comment on
      !! trc_set_geometry -- so this must converge, and to the SAME energy a
      !! fresh SAD run at geometry B does.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h_proj, h_sad
      integer(c_int) :: rc, rc_run, conv_proj, conv_sad
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), at_r_moved(:, :)
      real(c_double) :: e_proj, e_sad

      call read_xyz("water.xyz", natm, zint, at_r)
      allocate (at_r_moved, source=at_r)
      at_r_moved(3, 1) = at_r_moved(3, 1) + 0.3_dp

      rc = trc_create(h_proj)
      rc = trc_set_molecule(h_proj, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                             0_c_int, 1_c_int)
      rc = trc_set_basis_json(h_proj, "basis_sets/6-31g.json"//c_null_char)
      rc = trc_set_guess(h_proj, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
      rc = trc_run_scf(h_proj)
      call check(error, rc == TRC_OK, "the geometry-A SCF through the C API must converge")
      if (allocated(error)) then
         rc = trc_destroy(h_proj)
         return
      end if

      rc = trc_set_geometry(h_proj, reshape(at_r_moved, [3*natm]))
      call check(error, rc == TRC_OK, "trc_set_geometry to the moved geometry must succeed")
      if (allocated(error)) then
         rc = trc_destroy(h_proj)
         return
      end if

      rc = trc_set_guess(h_proj, TRC_GUESS_PROJECTED, c_null_ptr, 1_c_int)
      rc_run = trc_run_scf(h_proj)
      call check(error, rc_run == TRC_OK, "TRC_GUESS_PROJECTED after trc_set_geometry must converge")
      if (allocated(error)) then
         rc = trc_destroy(h_proj)
         return
      end if
      rc = trc_converged(h_proj, conv_proj)
      rc = trc_energy(h_proj, e_proj)
      rc = trc_destroy(h_proj)

      ! An independent context started fresh at the SAME moved geometry
      ! with the default SAD guess -- the stationary point the projected
      ! warm start must agree with.
      rc = trc_create(h_sad)
      rc = trc_set_molecule(h_sad, int(natm, c_int), real(zint, c_double), reshape(at_r_moved, [3*natm]), &
                             0_c_int, 1_c_int)
      rc = trc_set_basis_json(h_sad, "basis_sets/6-31g.json"//c_null_char)
      rc = trc_set_guess(h_sad, TRC_GUESS_SAD, c_null_ptr, 1_c_int)
      rc = trc_run_scf(h_sad)
      call check(error, rc == TRC_OK, "the reference SAD SCF at the moved geometry must converge")
      if (allocated(error)) then
         rc = trc_destroy(h_sad)
         return
      end if
      rc = trc_converged(h_sad, conv_sad)
      rc = trc_energy(h_sad, e_sad)
      rc = trc_destroy(h_sad)

      call check(error, conv_proj /= 0_c_int, "trc_converged must report the projected run as converged")
      if (allocated(error)) return
      call check(error, conv_sad /= 0_c_int, "trc_converged must report the SAD reference run as converged")
      if (allocated(error)) return
      call check(error, real(e_proj, dp), real(e_sad, dp), thr=TOL_ENERGY, &
                 message="TRC_GUESS_PROJECTED after trc_set_geometry converges to a different "// &
                 "energy than a fresh SAD run")
   end subroutine test_capi_projected_after_set_geometry

end module check_scf_reuse_suite

program check_scf_reuse
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_api, only: trc_bind_device
   use check_scf_reuse_suite, only: collect_scf_reuse
   implicit none
   integer :: stat, dev
   stat = 0
   dev = trc_bind_device(0)
   call run_testsuite(collect_scf_reuse, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_scf_reuse
