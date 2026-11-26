module check_cpks_batched_suite
   !!
   !! The batched coupled-perturbed solvers (trc_cphf_solve_rhf and
   !! trc_cpks_solve_rks with opts%max_batch > 1, `trc_cphf_cg_t`) against a
   !! reference that shares nothing with them: this file's own textbook
   !! preconditioned CG, one column at a time, with its own matvec built from
   !! `eri%fock_resident` (and `trc_xc_fxc_rks` for KS). So a defect in the
   !! batched CG state machine (alpha, beta, the active set), which a
   !! comparison against max_batch = 1 would share, changes the iteration
   !! counts or the solution here. Each batched solution is also checked
   !! against the tolerance by recomputing its residual (A + D) u + b with
   !! that independent matvec.
   !!
   !! The right-hand sides are chosen to converge at DIFFERENT iterations (a
   !! column scaled to below the tolerance, ones scaled by 1e-6, 1e-3, 1 and
   !! 1e2, a localised one): a batch that shares alpha, beta or a residual
   !! across columns, or mishandles a column that leaves the active set, moves
   !! the unequal columns by far more than rounding. Water, 6-31G*, closed
   !! shell: HF and B3LYP (a hybrid GGA: the shared Fock pass with exact
   !! exchange and the f_xc response; coarse grid, since only consistency is
   !! tested).
   !!
   use pic_types, only: dp
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t, trc_fock_options_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc, only: trc_xc_fxc_rks
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_cphf_solve_rhf, &
                       cphf_zero_matrix, cphf_neg_matrix, cphf_axpy_matrix, cphf_precond_matrix, &
                       cphf_combine_p, cphf_add_diag, cphf_copy_matrix
   use trc_cpks, only: trc_cpks_solve_rks
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_cpks_batched

   real(dp), parameter :: ERI_THRESH = 1.0e-12_dp
   real(dp), parameter :: CONV = 1.0e-9_dp
   integer, parameter :: NRHS = 6
   real(dp), parameter :: RHS_SCALE(NRHS) = [1.0_dp, 1.0e-3_dp, 1.0e2_dp, 1.0e-12_dp, 1.0_dp, 1.0e-6_dp]
   !! Column 4 starts below the tolerance (zero iterations), column 5 is a
   !! single-element vector, the rest are smooth patterns of unlike size.

   type :: fixture_t
      !! One converged SCF (HF or B3LYP) with everything a CP solve needs,
      !! the device scratch of the reference matvec, and the reference
      !! solutions of the NRHS columns.
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
      real(dp), allocatable :: wd(:, :), wdt(:, :), wg(:, :), wv(:, :), wt1(:, :), wt2(:, :)
      real(dp), allocatable :: vb(:, :), vu(:, :), vr(:, :), vz(:, :), vp(:, :), vap(:, :), vdiag(:, :)
      real(dp), allocatable :: bref(:, :, :), uref(:, :, :)
      integer :: itref(NRHS) = 0
   end type fixture_t

   type(fixture_t), save, target :: fx_hf, fx_ks

contains

   subroutine collect_cpks_batched(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("HF: batched = independent CG, unequal convergence", test_hf_agree), &
                  new_unittest("B3LYP: batched = independent CG, unequal convergence", test_ks_agree), &
                  new_unittest("group size does not change the answer", test_group_sizes), &
                  new_unittest("iteration cap: every column reports its own count", test_max_iteration), &
                  new_unittest("breakdown on an indefinite operator", test_breakdown_indefinite), &
                  new_unittest("a NaN column does not poison its neighbours", test_nan_column), &
                  new_unittest("release the fixtures", test_release) &
                  ]
   end subroutine collect_cpks_batched

   subroutine setup(fx, error, name)
      !! Converge the SCF `name` ("hf" or "b3lyp") once, on water/6-31G*, and
      !! solve the NRHS columns with the independent CG.
      type(fixture_t), intent(inout) :: fx
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name
      integer :: natm, dev, c
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      logical :: ok
      integer :: nv, no, n

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
         ! Coarse grid: the CP operator is whatever the grid makes it, and
         ! only the batched against the independent solve is tested.
         opts%grid_level = 1
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
      n = fx%nao
      nv = fx%nvirt
      no = fx%nocc
      allocate (fx%cmo(n, n), fx%eps(n), fx%dmat(n, n))
      fx%cmo = res%cmo(:, :, 1)
      fx%eps = res%eps(:, 1)
      fx%dmat = res%dmat(:, :, 1)
      allocate (fx%wd(n, n), fx%wdt(n, n), fx%wg(n, n), fx%wv(n, n), fx%wt1(n, no), fx%wt2(n, no))
      allocate (fx%vb(nv, no), fx%vu(nv, no), fx%vr(nv, no), fx%vz(nv, no), fx%vp(nv, no), fx%vap(nv, no), &
                fx%vdiag(nv, no))
      do c = 1, no
         fx%vdiag(:, c) = fx%eps(no + 1:n) - fx%eps(c)
      end do
      !$acc enter data copyin(fx%cmo, fx%dmat, fx%vdiag) create(fx%wd, fx%wdt, fx%wg, fx%wv, fx%wt1, fx%wt2, &
      !$acc fx%vb, fx%vu, fx%vr, fx%vz, fx%vp, fx%vap)
      fx%ready = .true.

      allocate (fx%bref(nv, no, NRHS), fx%uref(nv, no, NRHS))
      call rhs(fx, NRHS, fx%bref)
      do c = 1, NRHS
         call reference_cg(fx, fx%bref(:, :, c), CONV, 100, fx%uref(:, :, c), fx%itref(c))
      end do
      call global_logger%info("cpks_batched "//trim(fx%name)//": independent CG iterations "// &
                              to_string(fx%itref(1))//" "//to_string(fx%itref(2))//" "//to_string(fx%itref(3))//" "// &
                              to_string(fx%itref(4))//" "//to_string(fx%itref(5))//" "//to_string(fx%itref(6)))
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

   subroutine apply_a(fx, pin, apout)
      !! apout = [A p] + (e_a - e_i) p on the device scratch of `fx`: the
      !! orbital Hessian (module header of trc_cphf/trc_cpks), built here
      !! from one eri%fock_resident (and f_xc for KS) per call.
      type(fixture_t), intent(inout) :: fx
      real(dp), intent(in) :: pin(fx%nvirt, fx%nocc)
      real(dp), intent(inout) :: apout(fx%nvirt, fx%nocc)
      integer :: n, no
      n = fx%nao
      no = fx%nocc
      call fx%la%gemm_mat(fx%cmo(:, no + 1:n), pin, fx%wt1)
      call fx%la%gemm_mat(fx%wt1, fx%cmo(:, 1:no), fx%wd, transb="T")
      call fx%la%transpose(n, n, fx%wd, n, fx%wdt, n)
      call cphf_axpy_matrix(n, n, 1.0_dp, fx%wdt, fx%wd)
      if (fx%ks) then
         call fx%eri%fock_resident(fx%bas, fx%wd, fx%wg, &
                                   opts=trc_fock_options_t(k_scale=fx%func%exx, j_scale=1.0_dp, density_screen=.false.))
         call trc_xc_fxc_rks(fx%bas, fx%xg, fx%func, fx%dmat, fx%wd, fx%wv, la=fx%la)
         call cphf_axpy_matrix(n, n, 1.0_dp, fx%wv, fx%wg)
      else
         call fx%eri%fock_resident(fx%bas, fx%wd, fx%wg, &
                                   opts=trc_fock_options_t(k_scale=1.0_dp, j_scale=1.0_dp, density_screen=.false.))
      end if
      call fx%la%gemm_mat(fx%wg, fx%cmo(:, 1:no), fx%wt2)
      call fx%la%gemm_mat(fx%cmo(:, no + 1:n), fx%wt2, apout, transa="T", alpha=2.0_dp)
      call cphf_add_diag(fx%nvirt, no, fx%vdiag, pin, apout)
   end subroutine apply_a

   subroutine reference_cg(fx, b, tol, maxit, u, its)
      !! Textbook diagonally preconditioned CG for (A + D) u = -b, one
      !! column, written against `apply_a` and nothing from trc_cphf_cg_t.
      type(fixture_t), intent(inout) :: fx
      real(dp), intent(in) :: b(fx%nvirt, fx%nocc)
      real(dp), intent(in) :: tol
      integer, intent(in) :: maxit
      real(dp), intent(out) :: u(fx%nvirt, fx%nocc)
      integer, intent(out) :: its
      real(dp) :: rz, rz_new, pap, alpha, beta, rnorm
      integer :: nv, no

      nv = fx%nvirt
      no = fx%nocc
      fx%vb = b
      !$acc update device(fx%vb)
      call cphf_neg_matrix(nv, no, fx%vb, fx%vr)
      call cphf_zero_matrix(nv, no, fx%vu)
      call cphf_precond_matrix(nv, no, fx%vdiag, fx%vr, fx%vz)
      call cphf_copy_matrix(nv, no, fx%vz, fx%vp)
      rz = fx%la%dot_2d(fx%vr, fx%vz)
      rnorm = sqrt(fx%la%dot_2d(fx%vr, fx%vr))
      its = 0
      do while (its < maxit .and. .not. (rnorm < tol))
         call apply_a(fx, fx%vp, fx%vap)
         pap = fx%la%dot_2d(fx%vp, fx%vap)
         if (.not. (pap > 0.0_dp)) exit
         alpha = rz/pap
         call cphf_axpy_matrix(nv, no, alpha, fx%vp, fx%vu)
         call cphf_axpy_matrix(nv, no, -alpha, fx%vap, fx%vr)
         rnorm = sqrt(fx%la%dot_2d(fx%vr, fx%vr))
         its = its + 1
         if (rnorm < tol) exit
         call cphf_precond_matrix(nv, no, fx%vdiag, fx%vr, fx%vz)
         rz_new = fx%la%dot_2d(fx%vr, fx%vz)
         beta = rz_new/rz
         rz = rz_new
         call cphf_combine_p(nv, no, fx%vz, beta, fx%vp)
      end do
      !$acc update self(fx%vu)
      u = fx%vu
   end subroutine reference_cg

   function true_residual(fx, u, b) result(r)
      !! || (A + D) u + b || recomputed with the independent matvec.
      type(fixture_t), intent(inout) :: fx
      real(dp), intent(in) :: u(fx%nvirt, fx%nocc), b(fx%nvirt, fx%nocc)
      real(dp) :: r
      fx%vp = u
      fx%vb = b
      !$acc update device(fx%vp, fx%vb)
      call apply_a(fx, fx%vp, fx%vap)
      call cphf_axpy_matrix(fx%nvirt, fx%nocc, 1.0_dp, fx%vb, fx%vap)
      r = sqrt(fx%la%dot_2d(fx%vap, fx%vap))
   end function true_residual

   subroutine solve(fx, ncol, bvec, u, opts, res, eps)
      !! The batched solver of `fx` for `ncol` columns; `eps` replaces the
      !! orbital energies (a deliberately wrong spectrum), absent the SCF's own.
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

   subroutine same_as_reference(fx, res, u, error, label)
      !! The batched columns' counts equal the independent CG's; the
      !! solutions agree to 1e-10; and the residual recomputed with the
      !! independent matvec is under the tolerance.
      !! The counts are compared exactly: on the host the two paths run the
      !! same arithmetic per column. (A device build may differ in the last
      !! bit of a reduction; if a count then moves by one at the threshold,
      !! that, and only that, is rounding.)
      type(fixture_t), intent(inout) :: fx
      type(trc_cphf_result_t), intent(in) :: res
      real(dp), intent(in) :: u(fx%nvirt, fx%nocc, NRHS)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer :: c, cw
      real(dp) :: w, r

      w = worst_diff(fx, NRHS, u, fx%uref, cw)
      call global_logger%info("cpks_batched "//label//": worst column diff "//to_string(w, "es10.2"))
      call check(error, w < 1.0e-10_dp, label//": batched vs independent CG, worst relative column diff "// &
                 to_string(w, "es10.2")//" at column "//to_string(cw))
      if (allocated(error)) return
      do c = 1, NRHS
         call check(error, res%iterations(c) == fx%itref(c), &
                    label//": column "//to_string(c)//": "//to_string(res%iterations(c))//" iterations batched, "// &
                    to_string(fx%itref(c))//" independent")
         if (allocated(error)) return
         r = true_residual(fx, u(:, :, c), fx%bref(:, :, c))
         call check(error, r < 2.0_dp*CONV, label//": column "//to_string(c)//": recomputed residual "// &
                    to_string(r, "es10.2")//" above the tolerance")
         if (allocated(error)) return
      end do
   end subroutine same_as_reference

   subroutine agree(fx, error)
      type(fixture_t), intent(inout) :: fx
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: ub(:, :, :)
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res

      call setup(fx, error, trim(fx%name))
      if (allocated(error)) return
      allocate (ub(fx%nvirt, fx%nocc, NRHS))
      opts%conv_residual = CONV
      opts%max_iteration = 100
      opts%max_batch = 8
      call solve(fx, NRHS, fx%bref, ub, opts, res)
      call check(error, res%converged, "batched solve did not converge: "//trim(res%message))
      if (allocated(error)) return
      call same_as_reference(fx, res, ub, error, trim(fx%name))
      if (allocated(error)) return
      ! The columns really do converge at unlike iterations, including none.
      call check(error, fx%itref(4) == 0, "the below-tolerance column must take no iteration")
      if (allocated(error)) return
      call check(error, all(ub(:, :, 4) == 0.0_dp), "the below-tolerance column's solution must stay zero")
      if (allocated(error)) return
      call check(error, maxval(fx%itref) - minval(fx%itref) >= 4, &
                 "the right-hand sides should converge at unlike iterations")
      if (allocated(error)) return
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
      !! max_batch 2 and 4 split the NRHS columns into three groups and two:
      !! both give the independent CG's answer (HF).
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: ug(:, :, :)
      integer, parameter :: GROUPS(2) = [2, 4]
      integer :: g
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res

      fx_hf%name = "hf"
      call setup(fx_hf, error, "hf")
      if (allocated(error)) return
      allocate (ug(fx_hf%nvirt, fx_hf%nocc, NRHS))
      opts%conv_residual = CONV
      do g = 1, size(GROUPS)
         opts%max_batch = GROUPS(g)
         call solve(fx_hf, NRHS, fx_hf%bref, ug, opts, res)
         call check(error, res%converged, "max_batch "//to_string(GROUPS(g))//" did not converge")
         if (allocated(error)) return
         call same_as_reference(fx_hf, res, ug, error, "max_batch "//to_string(GROUPS(g)))
         if (allocated(error)) return
      end do
   end subroutine test_group_sizes

   subroutine test_max_iteration(error)
      !! With the cap at 3 the hard columns stop at exactly 3 and report it,
      !! the easy one reports its own 0, and the capped iterates agree with
      !! the independent CG stopped at the same cap.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: ub(:, :, :), uc(:, :, :)
      integer :: its_ref(NRHS), c
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res
      real(dp) :: w

      fx_hf%name = "hf"
      call setup(fx_hf, error, "hf")
      if (allocated(error)) return
      allocate (ub(fx_hf%nvirt, fx_hf%nocc, NRHS), uc(fx_hf%nvirt, fx_hf%nocc, NRHS))
      opts%conv_residual = CONV
      opts%max_iteration = 3
      opts%max_batch = 8
      call solve(fx_hf, NRHS, fx_hf%bref, ub, opts, res)
      do c = 1, NRHS
         call reference_cg(fx_hf, fx_hf%bref(:, :, c), CONV, 3, uc(:, :, c), its_ref(c))
      end do
      call check(error, .not. res%converged, "a capped solve must not report converged")
      if (allocated(error)) return
      call check(error, index(res%message, "did not converge") > 0, "message: "//trim(res%message))
      if (allocated(error)) return
      call check(error, res%iterations(1) == 3 .and. res%iterations(3) == 3, &
                 "the hard columns must stop at the cap, got "//to_string(res%iterations(1))//" and "// &
                 to_string(res%iterations(3)))
      if (allocated(error)) return
      call check(error, all(res%iterations == its_ref), "capped iteration counts differ from the independent CG")
      if (allocated(error)) return
      call check(error, res%iterations(4) == 0, "the below-tolerance column must report 0")
      if (allocated(error)) return
      w = worst_diff(fx_hf, NRHS, ub, uc)
      call check(error, w < 1.0e-10_dp, "capped batched vs independent CG, worst diff "//to_string(w, "es10.2"))
   end subroutine test_max_iteration

   subroutine test_breakdown_indefinite(error)
      !! Orbital energies that put the virtuals far BELOW the occupied make
      !! the preconditioner (and the operator) negative: p^T A p <= 0 on the
      !! first step. Every column must stop there as a breakdown, with no
      !! iteration done, not fall through to a NaN alpha and "converge".
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: ub(:, :, :)
      real(dp), allocatable :: eps_bad(:)
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res

      fx_hf%name = "hf"
      call setup(fx_hf, error, "hf")
      if (allocated(error)) return
      allocate (ub(fx_hf%nvirt, fx_hf%nocc, 3), eps_bad(fx_hf%nao))
      eps_bad = fx_hf%eps
      eps_bad(fx_hf%nocc + 1:) = fx_hf%eps(1) - 20.0_dp
      opts%max_batch = 8
      call solve(fx_hf, 3, fx_hf%bref(:, :, 1:3), ub, opts, res, eps=eps_bad)
      call check(error, .not. res%converged, "an indefinite operator must not report converged")
      if (allocated(error)) return
      call check(error, index(res%message, "breakdown") > 0, "message: "//trim(res%message))
      if (allocated(error)) return
      call check(error, all(res%iterations == 0), "a breakdown on the first step leaves 0 iterations")
   end subroutine test_breakdown_indefinite

   subroutine test_nan_column(error)
      !! A NaN right-hand side in the middle of a batch (HF and B3LYP): that
      !! column breaks down, the others solve to what the independent CG gives.
      !! The digestion and the f_xc pass are per density, so one NaN density
      !! must stay in its own slot.
      type(error_type), allocatable, intent(out) :: error

      fx_hf%name = "hf"
      call setup(fx_hf, error, "hf")
      if (allocated(error)) return
      call nan_pass(fx_hf, "hf", error)
      if (allocated(error)) return
      fx_ks%name = "b3lyp"
      call setup(fx_ks, error, "b3lyp")
      if (allocated(error)) return
      call nan_pass(fx_ks, "b3lyp", error)
   end subroutine test_nan_column

   subroutine nan_pass(fx, label, error)
      type(fixture_t), intent(inout) :: fx
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: bvec(:, :, :), ub(:, :, :)
      type(trc_cphf_options_t) :: opts
      type(trc_cphf_result_t) :: res
      real(dp) :: w

      allocate (bvec(fx%nvirt, fx%nocc, 3), ub(fx%nvirt, fx%nocc, 3))
      bvec(:, :, 1) = fx%bref(:, :, 1)
      bvec(:, :, 2) = fx%bref(:, :, 2)
      bvec(:, :, 3) = fx%bref(:, :, 3)
      bvec(2, 1, 2) = ieee_nan()
      opts%conv_residual = CONV
      opts%max_batch = 8
      call solve(fx, 3, bvec, ub, opts, res)
      call check(error, .not. res%converged, label//": a NaN column must not report converged")
      if (allocated(error)) return
      call check(error, index(res%message, "breakdown") > 0, label//" message: "//trim(res%message))
      if (allocated(error)) return
      call check(error, res%iterations(2) == 0, label//": the NaN column iterates 0 times")
      if (allocated(error)) return
      call check(error, res%iterations(1) == fx%itref(1) .and. res%iterations(3) == fx%itref(3), &
                 label//": the neighbours' iteration counts moved")
      if (allocated(error)) return
      call check(error, all(ub(:, :, 1) == ub(:, :, 1)) .and. all(ub(:, :, 3) == ub(:, :, 3)), &
                 label//": a NaN leaked into a neighbouring column")
      if (allocated(error)) return
      w = max(maxval(abs(ub(:, :, 1) - fx%uref(:, :, 1))), maxval(abs(ub(:, :, 3) - fx%uref(:, :, 3))))
      w = w/max(1.0_dp, maxval(abs(fx%uref(:, :, 3))))
      call check(error, w < 1.0e-10_dp, label//": neighbours of a NaN column differ from the independent CG by "// &
                 to_string(w, "es10.2"))
   end subroutine nan_pass

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
      !$acc exit data delete(fx%cmo, fx%dmat, fx%vdiag, fx%wd, fx%wdt, fx%wg, fx%wv, fx%wt1, fx%wt2) &
      !$acc delete(fx%vb, fx%vu, fx%vr, fx%vz, fx%vp, fx%vap)
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
