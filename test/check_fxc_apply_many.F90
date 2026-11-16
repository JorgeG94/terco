module check_fxc_apply_many_suite
   !!
   !! trc_xc_fxc_rks_apply_many, the XC kernel response of many trial
   !! densities against one cached collocation (the TDDFT matvec's
   !! f_xc[D_k]), restricted, LDA and GGA.
   !!
   !! Three questions, one test-drive case each:
   !!  1. the batch is the loop: nvec random symmetric densities through
   !!     trc_xc_fxc_rks_apply_many equal one trc_xc_fxc_rks_apply call per
   !!     density, at every vector chunking `budget` can force (including a
   !!     remainder pass);
   !!  2. a NON-symmetric input gives what its symmetrised copy gives: f_xc
   !!     only sees (D + D^T)/2, and the point kernel's grad rho is wrong
   !!     for an unsymmetrised D;
   !!  3. pyscf's dft.numint.nr_rks_fxc (test/pyscf_refs/generate.py's
   !!     gen_fxc_apply_many) for a probe density, matrix element by element
   !!     and as the trace Tr(P v_resp).
   !!
   !! The reference density is D0 = 2 C C^T for a fixed pseudo-random C (so
   !! rho >= 0, as check_xc_energy's), the probes a fixed pseudo-random
   !! symmetric matrix and a fixed non-symmetric one: neither is an SCF
   !! output, so pyscf sees exactly these matrices. The grid is
   !! terco's own (level 3), written to the probe file for pyscf.
   !!
   !! Tolerances: batch vs loop and non-symmetric vs symmetrised 1e-13
   !! relative (the same arithmetic, so equal to roundoff); pyscf 1e-10
   !! absolute, against a measured 6e-13 on |v| up to 0.08 (the two codes'
   !! libxc and quadrature sums differ only in order).
   !!
   !! Arguments: xyz basis.json molname basname funcname probe.bin.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_cart_perm, only: cart_perm_build_ao
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_linalg, only: trc_linalg_t
   use trc_xc, only: trc_xc_fxc_cache_t, trc_xc_fxc_rks_apply, trc_xc_fxc_rks_apply_many, trc_xc_fxc_vec_chunk
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix
   implicit none
   private
   public :: collect_check_fxc_apply_many, xyz_arg, bas_json_arg, molname_arg, basname_arg, funcname_arg, probe_arg

   integer, parameter :: NVEC = 3
   real(dp), parameter :: ERI_THRESH = 1.0e-12_dp

   character(len=256) :: xyz_arg = 'water.xyz'
   character(len=256) :: bas_json_arg = 'basis_sets/6-31gs.json'
   character(len=16) :: molname_arg = 'water'
   character(len=16) :: basname_arg = '6-31gs'
   character(len=16) :: funcname_arg = 'pbe'
   character(len=256) :: probe_arg = 'fxc_apply_many_probe.bin'

   type :: fixture_t
      !! Everything the three cases share, built once per case.
      type(trc_basis_t) :: bas
      type(trc_xc_grid_t) :: xg
      type(trc_xc_functional_t) :: func
      type(trc_xc_fxc_cache_t) :: cache
      type(trc_linalg_t) :: la
      type(dft_grid_t) :: grid
      integer :: nao = 0, natm = 0
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), d0(:, :), psym(:, :), pnon(:, :), smat(:, :)
   end type fixture_t

contains

   subroutine collect_check_fxc_apply_many(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("batch equals the loop over trc_xc_fxc_rks_apply", test_batch_vs_loop), &
                  new_unittest("non-symmetric input equals its symmetrised copy", test_nonsymmetric), &
                  new_unittest("f_xc[D] vs pyscf nr_rks_fxc", test_vs_pyscf) &
                  ]
   end subroutine collect_check_fxc_apply_many

   function lcg_next(state) result(x)
      !! The next number of a fixed Lehmer generator, in (0, 1).
      integer(int64), intent(inout) :: state
      real(dp) :: x
      state = mod(48271_int64*state, 2147483647_int64)
      x = real(state, dp)/2147483647.0_dp
   end function lcg_next

   subroutine fixture_build(fx, error)
      type(fixture_t), intent(inout) :: fx
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: c_lex(:, :), tmat(:, :), vmat(:, :), plex(:, :), qlex(:, :), c(:, :)
      integer, allocatable :: perm_l2b(:), perm_b2l(:)
      integer(int64) :: lcg
      integer :: nao, nocc, i, j, dev
      logical :: ok

      dev = trc_bind_device(0)
      ok = xc_functional_by_name(trim(funcname_arg), fx%func)
      call check(error, ok, "unrecognised functional '"//trim(funcname_arg)//"'")
      if (allocated(error)) return
      call read_xyz(trim(xyz_arg), fx%natm, fx%zint, fx%at_r)
      call trc_basis_from_json(trim(bas_json_arg), fx%natm, fx%zint, fx%at_r, fx%bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call fx%bas%to_device()
      nao = fx%bas%nao
      fx%nao = nao
      nocc = nint(sum(real(fx%zint, dp)))/2

      call pl%build(fx%bas, ERI_THRESH)
      call pl%to_device()
      allocate (fx%smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      !$acc enter data create(fx%smat, tmat, vmat)
      call trc_1e(fx%bas, pl, fx%smat, tmat, vmat)
      !$acc update self(fx%smat)
      !$acc exit data delete(fx%smat, tmat, vmat)
      call pl%release()

      ! Drawn in a fixed lexicographic AO order, then permuted to this
      ! build's own (check_xc_energy's reason: a d shell's six components
      ! must draw the same numbers whatever order the build enumerates).
      allocate (perm_l2b(nao), perm_b2l(nao))
      call cart_perm_build_ao(fx%bas%nshell, fx%bas%sh_l, fx%bas%sh_ao, nao, perm_l2b, perm_b2l)
      allocate (c_lex(nao, nocc), c(nao, nocc), plex(nao, nao), qlex(nao, nao))
      allocate (fx%d0(nao, nao), fx%psym(nao, nao), fx%pnon(nao, nao))
      lcg = 97_int64
      do j = 1, nocc
         do i = 1, nao
            c_lex(i, j) = 0.3_dp*(lcg_next(lcg) - 0.5_dp)
         end do
      end do
      do j = 1, nocc
         do i = 1, nao
            c(i, j) = c_lex(perm_b2l(i), j)
         end do
      end do
      fx%d0 = 2.0_dp*matmul(c, transpose(c))
      do j = 1, nao
         do i = 1, nao
            plex(i, j) = 0.1_dp*(lcg_next(lcg) - 0.5_dp)
            qlex(i, j) = 0.1_dp*(lcg_next(lcg) - 0.5_dp)
         end do
      end do
      do j = 1, nao
         do i = 1, nao
            fx%psym(i, j) = 0.5_dp*(plex(perm_b2l(i), perm_b2l(j)) + plex(perm_b2l(j), perm_b2l(i)))
            fx%pnon(i, j) = qlex(perm_b2l(i), perm_b2l(j))
         end do
      end do

      call dft_grid_build(fx%bas%at_r, fx%bas%at_elem, fx%grid, err, level=3)
      call check(error, .not. err%has_error(), "grid: "//err%get_message())
      if (allocated(error)) return
      call fx%xg%build_from_grid(fx%grid, fx%bas)
      call fx%xg%to_device()
      call fx%la%init(nao)
      call fx%cache%setup(fx%bas, fx%xg, fx%func, fx%d0, fx%la)
      call check(error, fx%cache%ready, "fxc cache did not set up")
   end subroutine fixture_build

   subroutine fixture_release(fx)
      type(fixture_t), intent(inout) :: fx
      call fx%cache%release()
      call fx%la%release()
      call fx%xg%release()
      call fx%bas%release()
   end subroutine fixture_release

   subroutine apply_many_host(fx, nvec, dmats, v, budget)
      !! trc_xc_fxc_rks_apply_many on host arrays, mapped the way the
      !! device build requires (dmats and v resident).
      type(fixture_t), intent(inout) :: fx
      integer, intent(in) :: nvec
      real(dp), intent(in) :: dmats(:, :, :)
      real(dp), intent(out) :: v(:, :, :)
      integer(int64), intent(in), optional :: budget
      !$acc enter data copyin(dmats) create(v)
      call trc_xc_fxc_rks_apply_many(fx%bas, fx%xg, fx%cache, nvec, dmats, v, budget, fx%la)
      !$acc exit data copyout(v) delete(dmats)
   end subroutine apply_many_host

   subroutine apply_one_host(fx, dmat, v)
      type(fixture_t), intent(inout) :: fx
      real(dp), intent(in) :: dmat(:, :)
      real(dp), intent(out) :: v(:, :)
      !$acc enter data copyin(dmat) create(v)
      call trc_xc_fxc_rks_apply(fx%bas, fx%xg, fx%cache, dmat, v, la=fx%la)
      !$acc exit data copyout(v) delete(dmat)
   end subroutine apply_one_host

   subroutine random_symmetric_set(nao, nvec, dmats)
      !! nvec fixed pseudo-random symmetric matrices, indefinite.
      integer, intent(in) :: nao, nvec
      real(dp), intent(out) :: dmats(nao, nao, nvec)
      integer(int64) :: lcg
      integer :: i, j, k
      lcg = 31337_int64
      do k = 1, nvec
         do j = 1, nao
            do i = 1, j
               dmats(i, j, k) = 0.08_dp*(lcg_next(lcg) - 0.5_dp)
               dmats(j, i, k) = dmats(i, j, k)
            end do
         end do
      end do
   end subroutine random_symmetric_set

   subroutine test_batch_vs_loop(error)
      type(error_type), allocatable, intent(out) :: error
      type(fixture_t) :: fx
      real(dp), allocatable :: dmats(:, :, :), vb(:, :, :), v1(:, :), vb2(:, :, :)
      real(dp) :: scale, dev
      integer(int64) :: budgets(6)
      integer :: k, ib, nvc, seen(NVEC)

      call fixture_build(fx, error)
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      allocate (dmats(fx%nao, fx%nao, NVEC), vb(fx%nao, fx%nao, NVEC), vb2(fx%nao, fx%nao, NVEC), &
                v1(fx%nao, fx%nao))
      call random_symmetric_set(fx%nao, NVEC, dmats)
      call apply_many_host(fx, NVEC, dmats, vb)
      scale = maxval(abs(vb))
      call check(error, scale > 1.0e-6_dp, "f_xc response is trivially small, |v| = "//to_string(scale, "es10.2"))
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      dev = 0.0_dp
      do k = 1, NVEC
         call apply_one_host(fx, dmats(:, :, k), v1)
         dev = max(dev, maxval(abs(vb(:, :, k) - v1)))
      end do
      call check(error, dev < 1.0e-13_dp*max(scale, 1.0_dp), &
                 "batch vs loop, max|diff| = "//to_string(dev, "es10.2")//" for |v| "//to_string(scale, "es10.2"))
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if

      ! Every pass count the vector scratch cap can force: budgets from far
      ! below one density's scratch (one density per pass, three passes) up
      ! to the default. The same budget also shrinks the basis-value chunks,
      ! so the overflow batches' own chunking is hit too.
      budgets = [1_int64, 3000000_int64, 6000000_int64, 8388608_int64, 12500000_int64, 25000000_int64]
      seen = 0
      do ib = 1, size(budgets)
         nvc = trc_xc_fxc_vec_chunk(fx%bas, fx%xg, fx%cache, NVEC, budgets(ib))
         if (nvc >= 1 .and. nvc <= NVEC) seen(nvc) = seen(nvc) + 1
         call apply_many_host(fx, NVEC, dmats, vb2, budgets(ib))
         dev = maxval(abs(vb2 - vb))
         call check(error, dev < 1.0e-13_dp*max(scale, 1.0_dp), &
                    "budget "//to_string(budgets(ib))//" (vectors per pass "//to_string(nvc)// &
                    ") vs unchunked, max|diff| = "//to_string(dev, "es10.2"))
         if (allocated(error)) then
            call fixture_release(fx)
            return
         end if
      end do
      call check(error, seen(1) > 0 .and. seen(2) > 0 .and. seen(NVEC) > 0, &
                 "the budget sweep must hit one, two (a remainder pass) and all densities per pass, counts seen "// &
                 to_string(seen(1))//" "//to_string(seen(2))//" "//to_string(seen(3)))
      call fixture_release(fx)
   end subroutine test_batch_vs_loop

   subroutine test_nonsymmetric(error)
      type(error_type), allocatable, intent(out) :: error
      type(fixture_t) :: fx
      real(dp), allocatable :: dn(:, :, :), ds(:, :, :), vn(:, :, :), vs(:, :, :)
      real(dp) :: scale, dev, asym

      call fixture_build(fx, error)
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      allocate (dn(fx%nao, fx%nao, 2), ds(fx%nao, fx%nao, 2), vn(fx%nao, fx%nao, 2), vs(fx%nao, fx%nao, 2))
      dn(:, :, 1) = fx%pnon
      dn(:, :, 2) = fx%pnon + 0.05_dp*fx%psym
      ds(:, :, 1) = 0.5_dp*(dn(:, :, 1) + transpose(dn(:, :, 1)))
      ds(:, :, 2) = 0.5_dp*(dn(:, :, 2) + transpose(dn(:, :, 2)))
      asym = maxval(abs(dn(:, :, 1) - ds(:, :, 1)))
      call check(error, asym > 1.0e-3_dp, "the non-symmetric probe is not non-symmetric")
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      call apply_many_host(fx, 2, dn, vn)
      call apply_many_host(fx, 2, ds, vs)
      scale = maxval(abs(vs))
      dev = maxval(abs(vn - vs))
      call check(error, dev < 1.0e-13_dp*max(scale, 1.0_dp), &
                 "non-symmetric vs symmetrised, max|diff| = "//to_string(dev, "es10.2")// &
                 " for |v| "//to_string(scale, "es10.2"))
      call fixture_release(fx)
   end subroutine test_nonsymmetric

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(fixture_t) :: fx
      real(dp), allocatable :: dm(:, :, :), v(:, :, :), vref(:, :)
      real(dp) :: tr(2), tr_ref(2), dv, dtr, scale
      integer :: unit, u, k, ng
      character(len=256) :: refname

      call fixture_build(fx, error)
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      allocate (dm(fx%nao, fx%nao, 2), v(fx%nao, fx%nao, 2), vref(fx%nao, fx%nao))
      dm(:, :, 1) = fx%psym
      dm(:, :, 2) = fx%pnon
      call apply_many_host(fx, 2, dm, v)
      do k = 1, 2
         tr(k) = sum(dm(:, :, k)*v(:, :, k))
      end do

      ng = fx%grid%n_points
      open (newunit=unit, file=trim(probe_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(fx%nao, int64), int(fx%natm, int64), int(ng, int64)
      write (unit) molname_arg, basname_arg, funcname_arg
      write (unit) int(fx%zint, int64)
      write (unit) fx%at_r
      write (unit) fx%grid%coords
      write (unit) fx%grid%weights
      write (unit) fx%smat
      write (unit) fx%d0
      write (unit) dm
      write (unit) v
      close (unit)

      refname = 'fxc_apply_many_'//trim(molname_arg)//'_'//trim(basname_arg)//'_'//trim(funcname_arg)//'.ref'
      call ref_require(trim(refname), u)
      call ref_read_real(u, 'TR_P', tr_ref(1))
      call ref_read_real(u, 'TR_Q', tr_ref(2))
      scale = 0.0_dp
      dv = 0.0_dp
      do k = 1, 2
         if (k == 1) call ref_read_matrix(u, 'V_P', fx%nao, fx%nao, vref)
         if (k == 2) call ref_read_matrix(u, 'V_Q', fx%nao, fx%nao, vref)
         dv = max(dv, maxval(abs(v(:, :, k) - vref)))
         scale = max(scale, maxval(abs(vref)))
      end do
      call ref_close(u)
      dtr = maxval(abs(tr - tr_ref))
      call check(error, dv < 1.0e-10_dp, &
                 "v_resp vs pyscf nr_rks_fxc, max|diff| = "//to_string(dv, "es10.2")// &
                 " for |v| "//to_string(scale, "es10.2"))
      if (.not. allocated(error)) then
         call check(error, dtr < 1.0e-10_dp*max(maxval(abs(tr_ref)), 1.0_dp), &
                    "Tr(P v_resp) vs pyscf, max|diff| = "//to_string(dtr, "es10.2"))
      end if
      call fixture_release(fx)
   end subroutine test_vs_pyscf

end module check_fxc_apply_many_suite

program check_fxc_apply_many
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_fxc_apply_many_suite, only: collect_check_fxc_apply_many, xyz_arg, bas_json_arg, molname_arg, &
                                         basname_arg, funcname_arg, probe_arg
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 2) call get_command_argument(2, bas_json_arg)
   if (command_argument_count() >= 3) call get_command_argument(3, molname_arg)
   if (command_argument_count() >= 4) call get_command_argument(4, basname_arg)
   if (command_argument_count() >= 5) call get_command_argument(5, funcname_arg)
   if (command_argument_count() >= 6) call get_command_argument(6, probe_arg)

   stat = 0
   call run_testsuite(collect_check_fxc_apply_many, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_fxc_apply_many
