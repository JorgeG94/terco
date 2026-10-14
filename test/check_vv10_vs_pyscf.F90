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
program check_vv10_vs_pyscf
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_basis_json, only: trc_basis_from_json
   use trc_dft_grid, only: dft_grid_t, build_dft_grid
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_vv10, only: trc_vv10_rks, trc_vv10_uks
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none

   real(dp), parameter :: VV10_B = 6.0_dp, VV10_C = 0.01_dp
   !! wB97X-V's own NLC parameters (pyscf's dft.libxc.nlc_coeff('WB97X_V')).
   integer, parameter :: GRID_LEVEL = 1
   !! Deliberately coarse: this test checks the KERNEL, not the quadrature,
   !! and pass 2 is O(npts^2) on a build with no GPU to run it on.

   character(len=*), parameter :: xyz_water = 'water.xyz', xyz_oh = 'oh.xyz'
   character(len=*), parameter :: bas_631gs = 'basis_sets/6-31gs.json', bas_ccpvdz = 'basis_sets/cc-pvdz.json'
   integer :: unit

   open (newunit=unit, file='vv10_probe.bin', access='stream', form='unformatted', status='replace')
   write (unit) int(3, kind=8)   ! number of cases

   call run_rks(unit, xyz_water, bas_631gs, 5)
   call run_rks(unit, xyz_water, bas_ccpvdz, 5)
   call run_uks(unit, xyz_oh, bas_ccpvdz, 5, 4)

   close (unit)
   print '(a)', 'check_vv10_vs_pyscf: wrote vv10_probe.bin; run vv10_ref.py for the pyscf comparison'

contains

   subroutine run_rks(unit, xyzfile, basjson, nocc)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: xyzfile, basjson
      integer, intent(in) :: nocc
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
      integer(kind=8) :: lcg
      real(dp) :: seed

      call read_xyz(xyzfile, nat, at_z, at_r)
      call trc_basis_from_json(basjson, nat, at_z, at_r, bas, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_vv10_vs_pyscf: '//err%get_message(); stop 1
      end if
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

      call build_dft_grid(at_r, at_z, grid, err, level=GRID_LEVEL)
      if (err%has_error()) then
         print '(a)', 'check_vv10_vs_pyscf: '//err%get_message(); stop 1
      end if
      call xg%build(grid%n_points, grid%coords, grid%weights, bas, max_pts=256)
      call xg%to_device()

      allocate (c(nao, nocc), dmat(nao, nao), vmat0(nao, nao), vmat(nao, nao))
      lcg = 97_8
      do j = 1, nocc
         do i = 1, nao
            lcg = mod(48271_8*lcg, 2147483647_8)
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
      call trc_vv10_rks(bas, xg, dmat, VV10_B, VV10_C, e_nlc, vmat)
      !$acc update self(vmat)
      !$acc exit data delete(dmat, vmat)

      print '(a,a,a,a,a,i0,a,i0,a,f20.14)', '  ', trim(xyzfile), ' / ', trim(basjson), &
         ': nao ', nao, ' points ', grid%n_points, '  E_nlc(terco) = ', e_nlc

      write (unit) int(0, kind=8)   ! 0 = RKS case
      write (unit) int(nat, kind=8), int(nao, kind=8), int(grid%n_points, kind=8)
      write (unit) int(at_z, kind=8)
      write (unit) at_r
      write (unit) grid%coords
      write (unit) grid%weights
      write (unit) smat
      write (unit) dmat
      write (unit) vmat0
      write (unit) e_nlc
      write (unit) vmat

      call xg%release(); call grid%destroy(); call bas%release()
      deallocate (c, dmat, vmat0, vmat, smat, tmat, vnmat, at_z, at_r)
   end subroutine run_rks

   subroutine run_uks(unit, xyzfile, basjson, nalpha, nbeta)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: xyzfile, basjson
      integer, intent(in) :: nalpha, nbeta
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(error_t) :: err
      integer :: nat, nao, i, j
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: at_r(:, :), ca(:, :), cb(:, :)
      real(dp), allocatable :: dmat_a(:, :), dmat_b(:, :), vmat0_a(:, :), vmat0_b(:, :), vmat_a(:, :), vmat_b(:, :)
      real(dp), allocatable :: smat(:, :), tmat(:, :), vnmat(:, :)
      real(dp) :: e_nlc
      integer(kind=8) :: lcg
      real(dp) :: seed

      call read_xyz(xyzfile, nat, at_z, at_r)
      call trc_basis_from_json(basjson, nat, at_z, at_r, bas, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_vv10_vs_pyscf: '//err%get_message(); stop 1
      end if
      call bas%to_device()
      nao = bas%nao

      call pl%build(bas, 1.0e-30_dp)
      call pl%to_device()
      allocate (smat(nao, nao), tmat(nao, nao), vnmat(nao, nao))
      call trc_1e(bas, pl, smat, tmat, vnmat)
      call pl%release()

      call build_dft_grid(at_r, at_z, grid, err, level=GRID_LEVEL)
      if (err%has_error()) then
         print '(a)', 'check_vv10_vs_pyscf: '//err%get_message(); stop 1
      end if
      call xg%build(grid%n_points, grid%coords, grid%weights, bas, max_pts=256)
      call xg%to_device()

      allocate (ca(nao, nalpha), cb(nao, nbeta))
      allocate (dmat_a(nao, nao), dmat_b(nao, nao))
      allocate (vmat0_a(nao, nao), vmat0_b(nao, nao), vmat_a(nao, nao), vmat_b(nao, nao))
      lcg = 131_8
      do j = 1, nalpha
         do i = 1, nao
            lcg = mod(48271_8*lcg, 2147483647_8)
            seed = real(lcg, dp)/2147483647.0_dp
            ca(i, j) = 0.16_dp*(seed - 0.5_dp)
         end do
      end do
      do j = 1, nbeta
         do i = 1, nao
            lcg = mod(48271_8*lcg, 2147483647_8)
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
      !$acc enter data copyin(dmat_a, dmat_b, vmat_a, vmat_b)
      call trc_vv10_uks(bas, xg, dmat_a, dmat_b, VV10_B, VV10_C, e_nlc, vmat_a, vmat_b)
      !$acc update self(vmat_a, vmat_b)
      !$acc exit data delete(dmat_a, dmat_b, vmat_a, vmat_b)

      print '(a,a,a,a,a,i0,a,i0,a,f20.14)', '  ', trim(xyzfile), ' / ', trim(basjson), &
         ': nao ', nao, ' points ', grid%n_points, '  E_nlc(terco, UKS) = ', e_nlc

      write (unit) int(1, kind=8)   ! 1 = UKS case
      write (unit) int(nat, kind=8), int(nao, kind=8), int(grid%n_points, kind=8)
      write (unit) int(at_z, kind=8)
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

      call xg%release(); call grid%destroy(); call bas%release()
      deallocate (ca, cb, dmat_a, dmat_b, vmat0_a, vmat0_b, vmat_a, vmat_b, smat, tmat, vnmat, at_z, at_r)
   end subroutine run_uks

end program check_vv10_vs_pyscf
