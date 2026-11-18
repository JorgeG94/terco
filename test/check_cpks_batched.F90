module check_cpks_batched_suite
   !!
   !! The batched coupled-perturbed solvers (trc_cphf_solve_rhf and
   !! trc_cpks_solve_rks with opts%max_batch > 1, `trc_cphf_cg_t`) against the
   !! same right-hand sides solved ONE COLUMN AT A TIME (ncomp = 1 calls, which
   !! take the unbatched Fock and f_xc builds): the two must agree column by
   !! column, because the batch is one independent CG recurrence per column,
   !! only the matvecs share a pass over the integrals.
   !!
   !! The right-hand sides are chosen to converge at DIFFERENT iterations (a
   !! column scaled to below the tolerance, ones scaled by 1e-6, 1e-3, 1 and
   !! 1e2, a localised one): a batch that shares alpha, beta or a residual
   !! across columns, or mishandles a column that leaves the active set, moves
   !! the unequal columns by far more than rounding. Water, 6-31G*, closed
   !! shell: HF and B3LYP (a hybrid GGA: the shared Fock pass with exact exchange, and
   !! the cached many-density f_xc pass).
   !!
   use pic_types, only: dp
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_cphf_solve_rhf
   use trc_cpks, only: trc_cpks_solve_rks
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_cpks_batched

   real(dp), parameter :: ERI_THRESH = 1.0e-12_dp
   integer, parameter :: NRHS = 6
   real(dp), parameter :: RHS_SCALE(NRHS) = [1.0_dp, 1.0e-3_dp, 1.0e2_dp, 1.0e-12_dp, 1.0_dp, 1.0e-6_dp]
   !! Column 4 starts below the tolerance (zero iterations), column 5 is a
   !! single-element vector, the rest are smooth patterns of unlike size.

   type :: fixture_t
      !! One converged SCF (HF or B3LYP) with everything a CP solve needs.
      logical :: ready = .false.
      character(len=8) :: name = ""
      integer :: nao = 0, nocc = 0, nvirt = 0
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_linalg_t) :: la
      type(trc_xc_grid_t) :: xg
      type(trc_xc_functional_t) :: func
      real(dp), allocatable :: cmo(:, :), eps(:), dmat(:, :)
      logical :: ks = .false.
   end type fixture_t

   type(fixture_t), save, target :: fx_hf, fx_ks

contains

   subroutine collect_cpks_batched(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("HF: batched = one at a time, unequal convergence", test_hf_agree), &
                  new_unittest("B3LYP: batched = one at a time, unequal convergence", test_ks_agree), &
                  new_unittest("group size does not change the answer", test_group_sizes), &
                  new_unittest("iteration cap: every column reports its own count", test_max_iteration), &
                  new_unittest("breakdown on an indefinite operator", test_breakdown_indefinite), &
                  new_unittest("a NaN column does not poison its neighbours", test_nan_column), &
                  new_unittest("release the fixtures", test_release) &
                  ]
   end subroutine collect_cpks_batched

   subroutine setup(fx, error, name)
      !! Converge the SCF `name` ("hf" or "b3lyp") once, on water/6-31G*.
      type(fixture_t), intent(inout) :: fx
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name
      integer :: natm, dev
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      logical :: ok

      if (fx%ready) return
      dev = trc_bind_device(0)
      call read_xyz('water.xyz', natm, zint, at_r)
      call trc_basis_from_json('basis_sets/6-31gs.json', natm, zint, at_r, fx%bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call fx%bas%to_device()
      call fx%pl%build(fx%bas, ERI_THRESH)
      call fx%pl%to_device()
      fx%name = name
      fx%nao = fx%bas%nao
      fx%nocc = nint(sum(real(zint, dp)))/2
      fx%nvirt = fx%nao - fx%nocc
      fx%ks = name /= "hf"
      opts%conv_energy = 1.0e-12_dp
      opts%conv_density = 1.0e-10_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      if (fx%ks) then
         opts%functional = name
         opts%grid_level = 3
      end if
      call trc_scf_run(fx%bas, fx%nocc, fx%nocc, opts, res)
      call check(error, res%converged, "SCF did not converge for "//name)
      if (allocated(error)) return
      if (fx%ks) then
         ok = xc_functional_by_name(name, fx%func)
         call check(error, ok, "unrecognised functional "//name)
         if (allocated(error)) return
         call fx%xg%build(size(res%grid_w), res%grid_r, res%grid_w, fx%bas, max_pts=opts%grid_max_pts)
         call fx%xg%to_device()
      end if
      call fx%eri%build(fx%bas, ERI_THRESH)
      call fx%la%init(fx%nao)
      allocate (fx%cmo(fx%nao, fx%nao), fx%eps(fx%nao), fx%dmat(fx%nao, fx%nao))
      fx%cmo = res%cmo(:, :, 1)
      fx%eps = res%eps(:, 1)
      fx%dmat = res%dmat(:, :, 1)
      !$acc enter data copyin(fx%cmo)
      fx%ready = .true.
   end subroutine setup

   subroutine rhs(fx, ncol, bvec)
      !! Deterministic right-hand sides of unlike size and shape (see the
      !! header); column c is RHS_SCALE(c) times a pattern over (a, i).
      type(fixture_t), intent(in) :: fx
      integer, intent(in) :: ncol
      real(dp), intent(out) :: bvec(fx%nvirt, fx%nocc, ncol)
      integer :: a, i, c
      do c = 1, ncol
         do i = 1, fx%nocc
            do a = 1, fx%nvirt
               bvec(a, i, c) = sin(0.7_dp*real(a, dp) + 1.3_dp*real(i*c, dp)) &
                               + 0.3_dp*cos(0.11_dp*real(a*i*c, dp))
            end do
         end do
         if (c == 5) then
            bvec(:, :, c) = 0.0_dp
            bvec(2, 1, c) = 1.0_dp
         end if
         bvec(:, :, c) = RHS_SCALE(c)*bvec(:, :, c)
      end do
   end subroutine rhs

   subroutine solve(fx, ncol, bvec, u, opts, res, eps)
      !! The solver of `fx` for `ncol` columns; `eps` replaces the orbital
      !! energies (a deliberately wrong spectrum), absent the SCF's own.
      type(fixture_t), intent(inout) :: fx
      integer, intent(in) :: ncol
      real(dp), intent(in) :: bvec(fx%nvirt, fx%nocc, ncol)
      real(dp), intent(out) :: u(fx%nvirt, fx%nocc, ncol)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res
      real(dp), intent(in), optional :: eps(fx%nao)
      real(dp) :: e(fx%nao)

      e = fx%eps
      if (present(eps)) e = eps
      u = 0.0_dp
      !$acc enter data copyin(bvec, u)
      if (fx%ks) then
         call trc_cpks_solve_rks(fx%bas, fx%eri, fx%la, fx%xg, fx%func, fx%nao, fx%nocc, fx%nvirt, fx%cmo, e, &
                                 fx%dmat, ncol, bvec, u, opts, res)
      else
         call trc_cphf_solve_rhf(fx%bas, fx%eri, fx%la, fx%nao, fx%nocc, fx%nvirt, fx%cmo, e, ncol, bvec, u, &
                                 opts, res)
      end if
      !$acc update self(u)
      !$acc exit data delete(bvec, u)
   end subroutine solve

   subroutine one_at_a_time(fx, ncol, bvec, u, its, opts, all_converged)
      !! Each column through its own ncol = 1 call, the unbatched path.
      type(fixture_t), intent(inout) :: fx
      integer, intent(in) :: ncol
      real(dp), intent(in) :: bvec(fx%nvirt, fx%nocc, ncol)
      real(dp), intent(out) :: u(fx%nvirt, fx%nocc, ncol)
      integer, intent(out) :: its(ncol)
      type(trc_cphf_options_t), intent(in) :: opts
      logical, intent(out) :: all_converged
      real(dp), allocatable :: b1(:, :, :), u1(:, :, :)
      type(trc_cphf_result_t) :: r1
      type(trc_cphf_options_t) :: o1
      integer :: c

      o1 = opts
      o1%max_batch = 1
      all_converged = .true.
      allocate (b1(fx%nvirt, fx%nocc, 1), u1(fx%nvirt, fx%nocc, 1))
      do c = 1, ncol
         b1(:, :, 1) = bvec(:, :, c)
         call solve(fx, 1, b1, u1, o1, r1)
         u(:, :, c) = u1(:, :, 1)
         its(c) = r1%iterations(1)
         all_converged = all_converged .and. r1%converged
      end do
   end subroutine one_at_a_time

   function worst_diff(fx, ncol, u, v, c_worst) result(w)
      !! Largest column error, relative to max(1, the column's largest entry).
      type(fixture_t), intent(in) :: fx
      integer, intent(in) :: ncol
      real(dp), intent(in) :: u(fx%nvirt, fx%nocc, ncol), v(fx%nvirt, fx%nocc, ncol)
      integer, intent(out), optional :: c_worst
      real(dp) :: w, d
      integer :: c
      w = 0.0_dp
      if (present(c_worst)) c_worst = 1
      do c = 1, ncol
         d = maxval(abs(u(:, :, c) - v(:, :, c)))/max(1.0_dp, maxval(abs(v(:, :, c))))
         if (d > w) then
            w = d
            if (present(c_worst)) c_worst = c
         end if
      end do
   end function worst_diff

   subroutine agree(fx, error)
      !! Batched against one at a time on the NRHS columns.
      type(fixture_t), intent(inout) :: fx
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: bvec(:, :, :), ub(:, :, :), us(:, :, :)
      integer :: its_single(NRHS), c, cw
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res
      logical :: single_ok
      real(dp) :: w
      character(len=512) :: line

      call setup(fx, error, trim(fx%name))
      if (allocated(error)) return
      allocate (bvec(fx%nvirt, fx%nocc, NRHS), ub(fx%nvirt, fx%nocc, NRHS), us(fx%nvirt, fx%nocc, NRHS))
      call rhs(fx, NRHS, bvec)
      opts%conv_residual = 1.0e-9_dp
      opts%max_iteration = 100
      opts%max_batch = 8

      call solve(fx, NRHS, bvec, ub, opts, res)
      call one_at_a_time(fx, NRHS, bvec, us, its_single, opts, single_ok)
      call check(error, res%converged, "batched solve did not converge: "//trim(res%message))
      if (allocated(error)) return
      call check(error, single_ok, "the one-at-a-time solves did not converge")
      if (allocated(error)) return

      write (line, "(a,a,a,6(1x,i0),a,6(1x,i0))") "cpks_batched ", trim(fx%name), ": iterations batched", &
         res%iterations, "  single", its_single
      call global_logger%info(trim(line))
      w = worst_diff(fx, NRHS, ub, us, cw)
      call global_logger%info("cpks_batched "//trim(fx%name)//": worst column diff "//to_string(w, "es10.2"))
      call check(error, w < 1.0e-10_dp, "batched vs one-at-a-time, worst relative column diff "// &
                 to_string(w, "es10.2")//" at column "//to_string(cw))
      if (allocated(error)) return
      do c = 1, NRHS
         call check(error, abs(res%iterations(c) - its_single(c)) <= 1, &
                    "column "//to_string(c)//": "//to_string(res%iterations(c))//" iterations batched, "// &
                    to_string(its_single(c))//" alone")
         if (allocated(error)) return
      end do
      ! The columns really do converge at unlike iterations, including none.
      call check(error, res%iterations(4) == 0, "the below-tolerance column must take no iteration")
      if (allocated(error)) return
      call check(error, all(ub(:, :, 4) == 0.0_dp), "the below-tolerance column's solution must stay zero")
      if (allocated(error)) return
      call check(error, maxval(res%iterations) - minval(res%iterations) >= 4, &
                 "the right-hand sides should converge at unlike iterations")
      if (allocated(error)) return
      ! A solution actually solves something: it is not the zero vector.
      call check(error, maxval(abs(ub(:, :, 3))) > 1.0e-3_dp, "the large column's solution is empty")
   end subroutine agree

   subroutine test_hf_agree(error)
      type(error_type), allocatable, intent(out) :: error
      fx_hf%name = "hf"
      call agree(fx_hf, error)
   end subroutine test_hf_agree

   subroutine test_ks_agree(error)
      type(error_type), allocatable, intent(out) :: error
      fx_ks%name = "b3lyp"
      call agree(fx_ks, error)
   end subroutine test_ks_agree

   subroutine test_group_sizes(error)
      !! max_batch 2 and 3 split the NRHS columns into several groups, 4 and
      !! 8 into one or two: all give one answer (B3LYP, so the f_xc cache is
      !! exercised too).
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: bvec(:, :, :), uref(:, :, :), ug(:, :, :)
      integer, parameter :: GROUPS(4) = [2, 3, 4, 8]
      integer :: g, c
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res, res_ref
      real(dp) :: w

      fx_ks%name = "b3lyp"
      call setup(fx_ks, error, "b3lyp")
      if (allocated(error)) return
      allocate (bvec(fx_ks%nvirt, fx_ks%nocc, NRHS), uref(fx_ks%nvirt, fx_ks%nocc, NRHS), &
                ug(fx_ks%nvirt, fx_ks%nocc, NRHS))
      call rhs(fx_ks, NRHS, bvec)
      opts%conv_residual = 1.0e-9_dp
      opts%max_batch = 8
      call solve(fx_ks, NRHS, bvec, uref, opts, res_ref)
      do g = 1, size(GROUPS)
         opts%max_batch = GROUPS(g)
         call solve(fx_ks, NRHS, bvec, ug, opts, res)
         call check(error, res%converged, "max_batch "//to_string(GROUPS(g))//" did not converge")
         if (allocated(error)) return
         w = worst_diff(fx_ks, NRHS, ug, uref)
         call check(error, w < 1.0e-10_dp, "max_batch "//to_string(GROUPS(g))//" differs from max_batch 8 by "// &
                    to_string(w, "es10.2"))
         if (allocated(error)) return
         do c = 1, NRHS
            call check(error, abs(res%iterations(c) - res_ref%iterations(c)) <= 1, &
                       "max_batch "//to_string(GROUPS(g))//": column "//to_string(c)//" iteration count moved")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_group_sizes

   subroutine test_max_iteration(error)
      !! With the cap at 3 the hard columns stop at exactly 3 and report it,
      !! the easy one reports its own 0, and the iterates (capped) still
      !! agree with the one-at-a-time ones.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: bvec(:, :, :), ub(:, :, :), us(:, :, :)
      integer :: its_single(NRHS)
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res
      logical :: single_ok
      real(dp) :: w

      fx_hf%name = "hf"
      call setup(fx_hf, error, "hf")
      if (allocated(error)) return
      allocate (bvec(fx_hf%nvirt, fx_hf%nocc, NRHS), ub(fx_hf%nvirt, fx_hf%nocc, NRHS), &
                us(fx_hf%nvirt, fx_hf%nocc, NRHS))
      call rhs(fx_hf, NRHS, bvec)
      opts%conv_residual = 1.0e-9_dp
      opts%max_iteration = 3
      opts%max_batch = 8
      call solve(fx_hf, NRHS, bvec, ub, opts, res)
      call one_at_a_time(fx_hf, NRHS, bvec, us, its_single, opts, single_ok)
      call check(error, .not. res%converged, "a capped solve must not report converged")
      if (allocated(error)) return
      call check(error, index(res%message, "did not converge") > 0, "message: "//trim(res%message))
      if (allocated(error)) return
      call check(error, .not. single_ok, "the capped one-at-a-time solves should not converge either")
      if (allocated(error)) return
      call check(error, res%iterations(1) == 3 .and. res%iterations(3) == 3, &
                 "the hard columns must stop at the cap, got "//to_string(res%iterations(1))//" and "// &
                 to_string(res%iterations(3)))
      if (allocated(error)) return
      call check(error, all(res%iterations == its_single), "capped iteration counts differ from one at a time")
      if (allocated(error)) return
      call check(error, res%iterations(4) == 0, "the below-tolerance column must report 0")
      if (allocated(error)) return
      w = worst_diff(fx_hf, NRHS, ub, us)
      call check(error, w < 1.0e-10_dp, "capped batched vs one-at-a-time, worst diff "//to_string(w, "es10.2"))
   end subroutine test_max_iteration

   subroutine test_breakdown_indefinite(error)
      !! Orbital energies that put the virtuals far BELOW the occupied make
      !! the preconditioner (and the operator) negative: p^T A p <= 0 on the
      !! first step. Every column must stop there as a breakdown, with no
      !! iteration done, not fall through to a NaN alpha and "converge".
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: bvec(:, :, :), ub(:, :, :)
      real(dp) :: eps_bad(fx_hf%nao)
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res

      fx_hf%name = "hf"
      call setup(fx_hf, error, "hf")
      if (allocated(error)) return
      allocate (bvec(fx_hf%nvirt, fx_hf%nocc, 3), ub(fx_hf%nvirt, fx_hf%nocc, 3))
      call rhs(fx_hf, 3, bvec)
      eps_bad = fx_hf%eps
      eps_bad(fx_hf%nocc + 1:) = fx_hf%eps(1) - 20.0_dp
      opts%max_batch = 8
      call solve(fx_hf, 3, bvec, ub, opts, res, eps=eps_bad)
      call check(error, .not. res%converged, "an indefinite operator must not report converged")
      if (allocated(error)) return
      call check(error, index(res%message, "breakdown") > 0, "message: "//trim(res%message))
      if (allocated(error)) return
      call check(error, all(res%iterations == 0), "a breakdown on the first step leaves 0 iterations")
   end subroutine test_breakdown_indefinite

   subroutine test_nan_column(error)
      !! A NaN right-hand side in the middle of a batch (HF and B3LYP): that
      !! column breaks down, the others solve to what they give alone. The
      !! digestion and the f_xc pass are per density, so one NaN density must
      !! stay in its own slot.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: bvec(:, :, :), ub(:, :, :), us(:, :, :)
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res
      integer :: its_single(3), pass
      logical :: single_ok
      real(dp) :: w
      type(fixture_t), pointer :: fx
      character(len=5) :: names(2)

      names = ["hf   ", "b3lyp"]
      do pass = 1, 2
         if (pass == 1) then
            fx_hf%name = "hf"
            call setup(fx_hf, error, "hf")
            fx => fx_hf
         else
            fx_ks%name = "b3lyp"
            call setup(fx_ks, error, "b3lyp")
            fx => fx_ks
         end if
         if (allocated(error)) return
         allocate (bvec(fx%nvirt, fx%nocc, 3), ub(fx%nvirt, fx%nocc, 3), us(fx%nvirt, fx%nocc, 3))
         call rhs(fx, 3, bvec)
         bvec(2, 1, 2) = ieee_nan()
         opts%conv_residual = 1.0e-9_dp
         opts%max_batch = 8
         call solve(fx, 3, bvec, ub, opts, res)
         call check(error, .not. res%converged, trim(names(pass))//": a NaN column must not report converged")
         if (allocated(error)) return
         call check(error, index(res%message, "breakdown") > 0, trim(names(pass))//" message: "//trim(res%message))
         if (allocated(error)) return
         call check(error, res%iterations(2) == 0, trim(names(pass))//": the NaN column iterates 0 times")
         if (allocated(error)) return
         ! The neighbours, against their own one-at-a-time solves.
         bvec(:, :, 2) = bvec(:, :, 3)
         call one_at_a_time(fx, 3, bvec, us, its_single, opts, single_ok)
         call check(error, single_ok, "reference solves did not converge")
         if (allocated(error)) return
         call check(error, all(ub(:, :, 1) == ub(:, :, 1)) .and. all(ub(:, :, 3) == ub(:, :, 3)), &
                    trim(names(pass))//": a NaN leaked into a neighbouring column")
         if (allocated(error)) return
         w = max(maxval(abs(ub(:, :, 1) - us(:, :, 1))), maxval(abs(ub(:, :, 3) - us(:, :, 3))))
         w = w/max(1.0_dp, maxval(abs(us(:, :, 3))))
         call check(error, w < 1.0e-10_dp, trim(names(pass))//": neighbours of a NaN column differ from alone by "// &
                    to_string(w, "es10.2"))
         if (allocated(error)) return
         deallocate (bvec, ub, us)
      end do
   end subroutine test_nan_column

   function ieee_nan() result(x)
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      real(dp) :: x
      x = ieee_value(x, ieee_quiet_nan)
   end function ieee_nan

   subroutine test_release(error)
      type(error_type), allocatable, intent(out) :: error
      if (fx_hf%ready) call release_fixture(fx_hf)
      if (fx_ks%ready) call release_fixture(fx_ks)
      call check(error, .not. (fx_hf%ready .or. fx_ks%ready), "fixtures released")
   end subroutine test_release

   subroutine release_fixture(fx)
      type(fixture_t), intent(inout) :: fx
      !$acc exit data delete(fx%cmo)
      call fx%eri%release()
      call fx%la%release()
      if (fx%ks) call fx%xg%release()
      call fx%pl%release()
      call fx%bas%release()
      fx%ready = .false.
   end subroutine release_fixture

end module check_cpks_batched_suite

program check_cpks_batched
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cpks_batched_suite, only: collect_cpks_batched
   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(collect_cpks_batched, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cpks_batched
