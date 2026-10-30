!
! VV10 (trc_vv10) against pyscf's own `numint.nr_nlc_vxc`, on the SAME grid
! points, weights and density -- same reasoning as check_xc_energy.F90, for
! the nonlocal piece that module does not cover.
!
! THREE CASES: water/6-31G* and water/cc-pVDZ through `trc_vv10_rks`
! (closed shell), the OH radical/cc-pVDZ through `trc_vv10_uks` (D_a != D_b,
! VV10 sees only D_a + D_b). Every basis is read `force_cart=.true.` --
! pyscf's cart=True Mole, same convention check_df_vs_pyscf.F90 and
! rys_fock_ref.py use -- and the density is a fixed pseudo-random PSD
! matrix (C C^T for a seeded C), not a converged SCF one: this is a kernel
! check on an arbitrary density, not a physical-energy check, the same
! reasoning check_xc_energy.F90 and check_rys_grad_vs_pyscf.F90 both give.
!
! `vmat` is SEEDED with an arbitrary nonzero matrix before the call and
! dumped alongside the result, so `vv10_ref.py` can check the ADD
! contract (trc_vv10_rks/uks must add onto an existing potential, never
! overwrite it) in the same comparison that checks the numbers: it is
! `vmat_after - vmat_seed` that must match pyscf's nlc matrix, not
! `vmat_after` on its own.
!
! AO ORDER/NORMALISATION: see rys_grad_ref.py's `terco_to_pyscf_perm` and
! its own header -- pyscf's cart=True groups a given atom's shells by
! ascending l where terco keeps the JSON's literal order, and a Cartesian
! d's off-diagonal components (xy/xz/yz) are not normalised the same way
! in the two codes. `vv10_ref.py` handles both exactly as that script
! does: permute, then rescale by n_i = sqrt(S_terco(i,i)/S_pyscf(i,i)).
!
! ONE shared pass writes vv10_probe.bin and reads vv10.ref (both files
! are sequential, so all three cases have to go through them in the SAME
! order generate.py wrote/expects); the three test-drive cases below each
! read back their OWN already-computed result -- see
! check_xc_energy_suite's own `next_ifn` comment for why the shared pass
! runs exactly once and this suite stays `parallel=.false.`.
!
module check_vv10_vs_pyscf_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_vv10, only: trc_vv10_rks, trc_vv10_uks
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_matrix
   implicit none
   private
   public :: collect_vv10_vs_pyscf

   real(dp), parameter :: VV10_B = 6.0_dp, VV10_C = 0.01_dp
   !! wB97X-V's own NLC parameters (pyscf's dft.libxc.nlc_coeff('WB97X_V')).
   type(trc_xc_functional_t), parameter :: VV10_FUNC = trc_xc_functional_t(nlc_b=VV10_B, nlc_c=VV10_C)
   !! Only its VV10 b and C are read.
   integer, parameter :: GRID_LEVEL = 1
   !! Deliberately coarse: this test checks the KERNEL, not the quadrature,
   !! and pass 2 is O(npts^2) on a build with no GPU to run it on.
   integer, parameter :: NCASE = 3

   logical :: computed = .false.
   logical :: ok_case(NCASE)
   character(len=256) :: msg_case(NCASE)

contains

   subroutine collect_vv10_vs_pyscf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water/6-31G* (RKS) vs pyscf", test_case_0), &
                  new_unittest("water/cc-pVDZ (RKS) vs pyscf", test_case_1), &
                  new_unittest("oh/cc-pVDZ (UKS) vs pyscf", test_case_2) &
                  ]
   end subroutine collect_vv10_vs_pyscf

   subroutine test_case_0(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, ok_case(1), trim(msg_case(1)))
   end subroutine test_case_0

   subroutine test_case_1(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, ok_case(2), trim(msg_case(2)))
   end subroutine test_case_1

   subroutine test_case_2(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, ok_case(3), trim(msg_case(3)))
   end subroutine test_case_2

   subroutine ensure_computed()
      if (computed) return
      call run_all()
      computed = .true.
   end subroutine ensure_computed

   subroutine run_all()
      integer :: unit, refu, ncase_ref
      logical :: has_ref

      ok_case = .false.
      msg_case = "pyscf reference missing"

      open (newunit=unit, file='vv10_probe.bin', access='stream', form='unformatted', status='replace')
      write (unit) int(NCASE, int64)   ! number of cases

      ! Opened (not required) BEFORE the dump is written: a missing reference
      ! must not stop this suite before vv10_probe.bin exists, or
      ! test/pyscf_refs/generate.py has nothing to bootstrap a fresh reference
      ! from. The comparisons below are what actually enforce the contract --
      ! FAIL, not silently skip, once the dump is safely written.
      call ref_open('vv10.ref', refu, has_ref)
      if (has_ref) then
         call ref_read_int(refu, 'NCASE', ncase_ref)
         if (ncase_ref /= NCASE) error stop "check_vv10_vs_pyscf: pyscf reference case count does not match this build"
      end if

      call run_rks(unit, refu, has_ref, 'water.xyz', 'basis_sets/6-31gs.json', 5, 0)
      call run_rks(unit, refu, has_ref, 'water.xyz', 'basis_sets/cc-pvdz.json', 5, 1)
      call run_uks(unit, refu, has_ref, 'oh.xyz', 'basis_sets/cc-pvdz.json', 5, 4, 2)

      close (unit)
      if (has_ref) call ref_close(refu)
   end subroutine run_all

   subroutine run_rks(unit, refu, has_ref, xyzfile, basjson, nocc, ic)
      integer, intent(in) :: unit, refu
      logical, intent(in) :: has_ref
      character(len=*), intent(in) :: xyzfile, basjson
      integer, intent(in) :: nocc, ic
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(error_t) :: err
      integer :: nat, nao, i, j
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: at_r(:, :), c(:, :), dmat(:, :), vmat0(:, :), vmat(:, :)
      real(dp), allocatable :: smat(:, :), tmat(:, :), vnmat(:, :)
      real(dp) :: e_nlc
      integer(int64) :: lcg
      real(dp) :: seed

      call read_xyz(xyzfile, nat, at_z, at_r)
      ! overrides%prune = .false.: the per-AO rescaling below is an identity against
      ! pyscf's own (literal-column) overlap, so the AO basis has to be
      ! the same one, not a span-equivalent rotation of it.
      call trc_basis_from_json(basjson, nat, at_z, at_r, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(prune=.false.))
      if (err%has_error()) error stop 'check_vv10_vs_pyscf: '//err%get_message()
      call bas%to_device()
      nao = bas%nao

      ! The overlap, for vv10_ref.py's per-AO rescaling (n_i =
      ! sqrt(S_terco(i,i)/S_pyscf(i,i))) -- same identity rys_grad_ref.py
      ! and rys_fock_ref.py use, needed here because a Cartesian d's
      ! off-diagonal components are not normalised the same way pyscf's
      ! cart=True basis normalises them.
      call pl%build(bas, 1.0e-30_dp)
      call pl%to_device()
      allocate (smat(nao, nao), tmat(nao, nao), vnmat(nao, nao))
      call trc_1e(bas, pl, smat, tmat, vnmat)
      call pl%release()

      call dft_grid_build(at_r, at_z, grid, err, level=GRID_LEVEL)
      if (err%has_error()) error stop 'check_vv10_vs_pyscf: '//err%get_message()
      call xg%build(grid%n_points, grid%coords, grid%weights, bas, max_pts=256)
      call xg%to_device()

      allocate (c(nao, nocc), dmat(nao, nao), vmat0(nao, nao), vmat(nao, nao))
      lcg = 97_int64
      do j = 1, nocc
         do i = 1, nao
            lcg = mod(48271_int64*lcg, 2147483647_int64)
            seed = real(lcg, dp)/2147483647.0_dp
            c(i, j) = 0.16_dp*(seed - 0.5_dp)
         end do
      end do
      dmat = 2.0_dp*matmul(c, transpose(c))
      do j = 1, nao
         do i = 1, nao
            vmat0(i, j) = 1.0e-3_dp*real(i*j, dp)
         end do
      end do
      vmat = vmat0
      !$acc enter data copyin(dmat, vmat)
      call trc_vv10_rks(bas, xg, VV10_FUNC, dmat, e_nlc, vmat)
      !$acc update self(vmat)
      !$acc exit data delete(dmat, vmat)

      print '(a,a,a,a,a,i0,a,i0,a,f20.14)', '  ', trim(xyzfile), ' / ', trim(basjson), &
         ': nao ', nao, ' points ', grid%n_points, '  E_nlc(terco) = ', e_nlc

      write (unit) int(0, int64)   ! 0 = RKS case
      write (unit) int(nat, int64), int(nao, int64), int(grid%n_points, int64)
      write (unit) int(at_z, int64)
      write (unit) at_r
      write (unit) grid%coords
      write (unit) grid%weights
      write (unit) smat
      write (unit) dmat
      write (unit) vmat0
      write (unit) e_nlc
      write (unit) vmat

      if (has_ref) then
         block
            integer :: kind_ref
            real(dp) :: e_nlc_ref, dv
            real(dp), allocatable :: vmat_ref(:, :)
            character(len=2) :: tag
            write (tag, '(i0)') ic
            call ref_read_int(refu, 'KIND_'//trim(tag), kind_ref)
            call ref_read_real(refu, 'E_NLC_'//trim(tag), e_nlc_ref)
            allocate (vmat_ref(nao, nao))
            call ref_read_matrix(refu, 'VMAT_'//trim(tag), nao, nao, vmat_ref)
            dv = maxval(abs((vmat - vmat0) - vmat_ref))
            ok_case(ic + 1) = kind_ref == 0 .and. abs(e_nlc - e_nlc_ref) < 1.0e-9_dp .and. dv < 1.0e-9_dp
            write (msg_case(ic + 1), '(a,a,es10.2)') trim(xyzfile)//'/'//trim(basjson), &
               ': E_nlc/V_nlc (RKS) vs pyscf, worst = ', max(abs(e_nlc - e_nlc_ref), dv)
         end block
      end if

      call xg%release(); call grid%destroy(); call bas%release()
      deallocate (c, dmat, vmat0, vmat, smat, tmat, vnmat, at_z, at_r)
   end subroutine run_rks

   subroutine run_uks(unit, refu, has_ref, xyzfile, basjson, nalpha, nbeta, ic)
      integer, intent(in) :: unit, refu
      logical, intent(in) :: has_ref
      character(len=*), intent(in) :: xyzfile, basjson
      integer, intent(in) :: nalpha, nbeta, ic
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(error_t) :: err
      integer :: nat, nao, i, j
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: at_r(:, :), ca(:, :), cb(:, :)
      real(dp), allocatable :: dmat_a(:, :), dmat_b(:, :), vmat0_a(:, :), vmat0_b(:, :), vmat_a(:, :), vmat_b(:, :)
      real(dp), allocatable :: dmat2(:, :, :), vmat2(:, :, :)
      real(dp), allocatable :: smat(:, :), tmat(:, :), vnmat(:, :)
      real(dp) :: e_nlc
      integer(int64) :: lcg
      real(dp) :: seed

      call read_xyz(xyzfile, nat, at_z, at_r)
      ! overrides%prune = .false.: the per-AO rescaling below is an identity against
      ! pyscf's own (literal-column) overlap, so the AO basis has to be
      ! the same one, not a span-equivalent rotation of it.
      call trc_basis_from_json(basjson, nat, at_z, at_r, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(prune=.false.))
      if (err%has_error()) error stop 'check_vv10_vs_pyscf: '//err%get_message()
      call bas%to_device()
      nao = bas%nao

      call pl%build(bas, 1.0e-30_dp)
      call pl%to_device()
      allocate (smat(nao, nao), tmat(nao, nao), vnmat(nao, nao))
      call trc_1e(bas, pl, smat, tmat, vnmat)
      call pl%release()

      call dft_grid_build(at_r, at_z, grid, err, level=GRID_LEVEL)
      if (err%has_error()) error stop 'check_vv10_vs_pyscf: '//err%get_message()
      call xg%build(grid%n_points, grid%coords, grid%weights, bas, max_pts=256)
      call xg%to_device()

      allocate (ca(nao, nalpha), cb(nao, nbeta))
      allocate (dmat_a(nao, nao), dmat_b(nao, nao))
      allocate (vmat0_a(nao, nao), vmat0_b(nao, nao), vmat_a(nao, nao), vmat_b(nao, nao))
      lcg = 131_int64
      do j = 1, nalpha
         do i = 1, nao
            lcg = mod(48271_int64*lcg, 2147483647_int64)
            seed = real(lcg, dp)/2147483647.0_dp
            ca(i, j) = 0.16_dp*(seed - 0.5_dp)
         end do
      end do
      do j = 1, nbeta
         do i = 1, nao
            lcg = mod(48271_int64*lcg, 2147483647_int64)
            seed = real(lcg, dp)/2147483647.0_dp
            cb(i, j) = 0.18_dp*(seed - 0.5_dp)
         end do
      end do
      dmat_a = matmul(ca, transpose(ca))
      dmat_b = matmul(cb, transpose(cb))
      do j = 1, nao
         do i = 1, nao
            vmat0_a(i, j) = 1.0e-3_dp*real(i*j, dp)
            vmat0_b(i, j) = -2.0e-3_dp*real(i*j, dp)
         end do
      end do
      vmat_a = vmat0_a; vmat_b = vmat0_b
      ! trc_vv10_uks takes both spins in one (nao, nao, 2) array, as trc_xc_uks.
      allocate (dmat2(nao, nao, 2), vmat2(nao, nao, 2))
      dmat2(:, :, 1) = dmat_a; dmat2(:, :, 2) = dmat_b
      vmat2(:, :, 1) = vmat_a; vmat2(:, :, 2) = vmat_b
      !$acc enter data copyin(dmat2, vmat2)
      call trc_vv10_uks(bas, xg, VV10_FUNC, dmat2, e_nlc, vmat2)
      !$acc update self(vmat2)
      !$acc exit data delete(dmat2, vmat2)
      vmat_a = vmat2(:, :, 1); vmat_b = vmat2(:, :, 2)
      deallocate (dmat2, vmat2)

      print '(a,a,a,a,a,i0,a,i0,a,f20.14)', '  ', trim(xyzfile), ' / ', trim(basjson), &
         ': nao ', nao, ' points ', grid%n_points, '  E_nlc(terco, UKS) = ', e_nlc

      write (unit) int(1, int64)   ! 1 = UKS case
      write (unit) int(nat, int64), int(nao, int64), int(grid%n_points, int64)
      write (unit) int(at_z, int64)
      write (unit) at_r
      write (unit) grid%coords
      write (unit) grid%weights
      write (unit) smat
      write (unit) dmat_a
      write (unit) dmat_b
      write (unit) vmat0_a
      write (unit) vmat0_b
      write (unit) e_nlc
      write (unit) vmat_a
      write (unit) vmat_b

      if (has_ref) then
         block
            integer :: kind_ref
            real(dp) :: e_nlc_ref, dva, dvb
            real(dp), allocatable :: vmat_ref(:, :)
            character(len=2) :: tag
            write (tag, '(i0)') ic
            call ref_read_int(refu, 'KIND_'//trim(tag), kind_ref)
            call ref_read_real(refu, 'E_NLC_'//trim(tag), e_nlc_ref)
            allocate (vmat_ref(nao, nao))
            call ref_read_matrix(refu, 'VMAT_'//trim(tag), nao, nao, vmat_ref)
            dva = maxval(abs((vmat_a - vmat0_a) - vmat_ref))
            dvb = maxval(abs((vmat_b - vmat0_b) - vmat_ref))
            ok_case(ic + 1) = kind_ref == 1 .and. abs(e_nlc - e_nlc_ref) < 1.0e-9_dp .and. &
                              dva < 1.0e-9_dp .and. dvb < 1.0e-9_dp
            write (msg_case(ic + 1), '(a,a,es10.2)') trim(xyzfile)//'/'//trim(basjson), &
               ': E_nlc/V_nlc (UKS) vs pyscf, worst = ', max(abs(e_nlc - e_nlc_ref), max(dva, dvb))
         end block
      end if

      call xg%release(); call grid%destroy(); call bas%release()
      deallocate (ca, cb, dmat_a, dmat_b, vmat0_a, vmat0_b, vmat_a, vmat_b, smat, tmat, vnmat, at_z, at_r)
   end subroutine run_uks

end module check_vv10_vs_pyscf_suite

program check_vv10_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_vv10_vs_pyscf_suite, only: collect_vv10_vs_pyscf
   implicit none
   integer :: stat
   stat = 0
   ! Serial: the shared setup (probe write, reference read) must run once,
   ! in declared order -- see check_vv10_vs_pyscf_suite's own header.
   call run_testsuite(collect_vv10_vs_pyscf, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_vv10_vs_pyscf
