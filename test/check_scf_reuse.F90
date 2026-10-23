!
! Reuse across SCF/gradient calls, added for a sequence of nearby
! geometries at the same molecule/basis (the finite-difference Hessian,
! and later a geometry optimiser): `trc_pair_eri_cache_t`
! (trc_scf_context.F90) lets the SCF's own pair list and (outside --df)
! four-centre engine stay alive for a gradient right after it, instead of
! a second pl%build/eri%build at the SAME geometry and threshold
! (trc_scf_grad.F90's own `pair_eri_cache` argument); `project_density`
! (same module) is the orbital-projection warm start,
! trc_capi.F90's TRC_GUESS_PROJECTED.
!
! This program FAILS WITHOUT THAT CHANGE: `pair_eri_cache`/`cmo_project`/
! `project_density` do not exist before it, so this file does not even
! compile against the old interfaces. With the change, it checks that
! reuse is a pure performance path -- never a correctness one:
!
!   1. at one geometry, a gradient computed WITH the cache and WITHOUT it
!      must agree to round-off (same pairs/Schwarz/bins, built once or
!      twice);
!   2. the SAME cache variable, reused across a SECOND, very different
!      geometry (one atom moved far enough away that a different number
!      of shell pairs survives Schwarz screening -- the "capacity" this
!      PR's design handles by always fully rebuilding into the cache,
!      never by growing a stale array), gives the same agreement: a
!      stale rebuild would show up here, at the geometry AFTER the one
!      that filled the cache, not at the first one;
!   3. `trc_pair_eri_cache_t%release()` is safe called twice;
!   4. `project_density` at an UNCHANGED geometry (S' = S) reproduces the
!      SCF's own converged density to round-off -- the identity case
!      every orbital-projection warm start must degenerate to.
!
program check_scf_reuse
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run, trc_pair_eri_cache_t
   use trc_scf_grad, only: trc_scf_gradient
   use trc_scf_context, only: project_density
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none

   real(dp), parameter :: THRESH = 1.0e-12_dp
   real(dp), parameter :: TOL_GRAD = 1.0e-12_dp
   !! Two builds of the SAME pairs/Schwarz/bins at the SAME threshold,
   !! contracted against the SAME density -- round-off, not truncation.
   real(dp), parameter :: TOL_PROJ = 1.0e-9_dp
   !! `project_density` rebuilds S independently (its own pair list/trc_1e
   !! call, not the SCF's), and goes through an eigendecomposition
   !! (`sym_orthog`) the SCF's own density build never does -- looser than
   !! TOL_GRAD, still far tighter than any physical effect this test could
   !! be confusing with a bug.

   type(trc_basis_t) :: bas
   type(trc_pair_eri_cache_t), target :: pec
   !! TARGET: `scf_ctx_setup`'s own `pair_eri_cache` dummy carries TARGET
   !! (trc_scf_context.F90's own header), and the pointer it sets up
   !! inside the call (`this%pec => pair_eri_cache`) writes `built`/
   !! `thresh_built`/`pl`/`eri` into THIS variable's storage -- the actual
   !! argument needs TARGET too so that write is standard-conforming, the
   !! same as trc_capi.F90's `context_box%peri` (reached through a POINTER
   !! handle there, which gives it TARGET implicitly).
   integer, allocatable :: z(:)
   real(dp), allocatable :: r0(:, :), r1(:, :)
   integer :: nat, dev, nalpha, nbeta
   type(error_t) :: err
   logical :: any_fail

   any_fail = .false.
   dev = trc_bind_device(0)
   call read_xyz('water.xyz', nat, z, r0)
   nalpha = sum(z)/2; nbeta = nalpha

   ! --- 1: the reference geometry, reusing a freshly-built cache ----------
   call trc_basis_from_json('basis_sets/6-31g.json', nat, z, r0, bas, err)
   if (err%has_error()) then
      print '(a)', 'check_scf_reuse: '//err%get_message(); stop 1
   end if
   call bas%to_device()
   call run_point('equilibrium', bas, nat, nalpha, nbeta, pec, any_fail)
   call bas%release()

   ! --- 2: a very different geometry, reusing the SAME cache variable -----
   ! One H moved 30 Bohr along z: its shell pairs with the other two atoms
   ! fall far below THRESH's Schwarz bound (water's primitives have no
   ! exponent small enough to keep an overlap alive at that range), so
   ! npair here is smaller than at 1 -- the "capacity" case. If reusing
   ! `pec` across this geometry change ever skipped a rebuild (it must
   ! not: `pairs_build`/`eri_build` both self-release and rebuild in
   ! full, unconditionally -- see `trc_pair_eri_cache_t`'s own header),
   ! `run_point`'s cached-vs-fresh comparison HERE is what would catch it,
   ! not at 1.
   allocate (r1, source=r0)
   r1(3, 3) = r1(3, 3) + 30.0_dp
   call trc_basis_from_json('basis_sets/6-31g.json', nat, z, r1, bas, err)
   if (err%has_error()) then
      print '(a)', 'check_scf_reuse: '//err%get_message(); stop 1
   end if
   call bas%to_device()
   call run_point('one atom far away', bas, nat, nalpha, nbeta, pec, any_fail)
   call bas%release()

   call pec%release()
   call pec%release()   ! idempotent, same invariant as scf_ctx_t%release

   ! --- 4: the orbital-projection warm start's identity case --------------
   call check_projected_guess_identity(any_fail)

   if (any_fail) then
      print '(a)', '  RESULT: FAIL'
      stop 1
   end if
   print '(a)', '  RESULT: PASS'

contains

   subroutine run_point(label, b, nat_, nalpha_, nbeta_, pec_, any_fail_)
      character(*), intent(in) :: label
      type(trc_basis_t), intent(inout) :: b
      integer, intent(in) :: nat_, nalpha_, nbeta_
      type(trc_pair_eri_cache_t), intent(inout) :: pec_
      logical, intent(inout) :: any_fail_

      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      real(dp), allocatable :: grad_cached(:, :), grad_fresh(:, :)
      real(dp) :: worst
      integer :: nao

      opts%max_iter = 30
      opts%conv_energy = 1.0e-11_dp
      opts%conv_diis = 1.0e-8_dp
      call trc_scf_run(b, nalpha_, nbeta_, opts, res, pair_eri_cache=pec_)
      if (res%iterations <= 0) then
         print '(a)', 'check_scf_reuse: '//trim(label)//': the SCF never ran: '//trim(res%message)
         any_fail_ = .true.; return
      end if
      nao = b%nao
      allocate (grad_cached(3, nat_), grad_fresh(3, nat_))
      call lp%init(nao)
      call trc_scf_gradient(b, lp, res%nspin, nalpha_, nbeta_, res%cmo, res%eps, THRESH, grad_cached, &
                            pair_eri_cache=pec_)
      call trc_scf_gradient(b, lp, res%nspin, nalpha_, nbeta_, res%cmo, res%eps, THRESH, grad_fresh)
      call lp%release()

      worst = maxval(abs(grad_cached - grad_fresh))
      print '(a,a,a,l1,a,es10.2)', '  ', trim(label), ': converged ', res%converged, &
         '  cached-vs-fresh worst|diff| ', worst
      if (worst >= TOL_GRAD) then
         print '(a,es10.2,a,es10.2)', '    FAIL: worst|diff| ', worst, ' >= tol ', TOL_GRAD
         any_fail_ = .true.
      end if
      deallocate (grad_cached, grad_fresh)
   end subroutine run_point

   subroutine check_projected_guess_identity(any_fail_)
      logical, intent(inout) :: any_fail_
      type(trc_basis_t) :: b
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), h(:, :), v(:, :), dproj(:, :, :)
      real(dp) :: worst
      integer :: nao, nocc(2)

      call trc_basis_from_json('basis_sets/6-31g.json', nat, z, r0, b, err)
      if (err%has_error()) then
         print '(a)', 'check_scf_reuse: '//err%get_message()
         any_fail_ = .true.; return
      end if
      call b%to_device()
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iter = 60
      call trc_scf_run(b, nalpha, nbeta, opts, res)
      if (.not. res%converged) then
         print '(a)', 'check_scf_reuse: projected guess: the reference SCF did not converge'
         any_fail_ = .true.; return
      end if
      nao = b%nao
      allocate (s(nao, nao), h(nao, nao), v(nao, nao), dproj(nao, nao, res%nspin))
      call pl%build(b, THRESH)
      call pl%to_device()
      call trc_1e(b, pl, s, h, v)
      call pl%release()
      call lp%init(nao)
      nocc = [nalpha, nbeta]
      call project_density(lp, nao, res%nspin, nocc, res%cmo, s, dproj)
      call lp%release()
      worst = maxval(abs(dproj - res%dmat))
      print '(a,es10.2)', '  projected guess at the SAME geometry, worst|D_proj - D_scf| ', worst
      if (worst >= TOL_PROJ) then
         print '(a,es10.2,a,es10.2)', '    FAIL: worst|diff| ', worst, ' >= tol ', TOL_PROJ
         any_fail_ = .true.
      end if
      call b%release()
      deallocate (s, h, v, dproj)
   end subroutine check_projected_guess_identity

end program check_scf_reuse
