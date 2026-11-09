module check_mdj_suite
   !! The McMurchie-Davidson J engine (eri%coulomb, trc_mdj) against the
   !! four-centre Fock build's own J (eri%fock at k_scale = 0) on the same
   !! density and threshold. The two share nothing past the basis and the
   !! shell-pair bins: one forms Cartesian quartets by HGP or Rys, the other
   !! never forms a quartet, so agreement to 1e-11 of max|J| checks the E
   !! coefficients, the R batches, both stages' layouts and the screen at
   !! once. The algebra itself is checked against pyscf by
   !! scripts/md/reference.py, from the generator's own statement lists.
   !!
   !! The density is synthetic and symmetric. The first cases use a full one
   !! (every element order one), so no block is screened out by the density
   !! and every class runs. The decaying case scales D_ij by w_i w_j, w
   !! falling from 1 to 1e-12 along the long axis of five waters, so ||rho_P||
   !! of a pair at one end is orders above ||rho_Q|| at the other while I_P
   !! and I_Q are alike: the symmetric kernels' two screen directions
   !! (I_P ||rho_Q|| for Jt[P], I_Q ||rho_P|| for Jt[Q]) disagree. (A decay
   !! with |R_i - R_j| does not do this: it shrinks I_P with ||rho_P||.)
   !! Letting the gather direction's test stand for both fails this case at
   !! 2e-5 and passes every full-density one.
   !!
   !! Run with TRC_MDJ_SYM_LMAX=-1 in a TERCO_MDJ_GATHER_ALL build
   !! (terco/check_mdj_gather), every case checks the gather kernels instead.
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t
   use trc_tables, only: LMAX
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_mdj

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: TOL = 1.0e-11_dp
   !! relative to max|J|

contains

   subroutine collect_mdj(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water 6-31G* (d shells)", test_water_631gs), &
                  new_unittest("water cc-pVDZ (split general contraction)", test_water_ccpvdz), &
                  new_unittest("formamide def2-SVP", test_formamide_svp), &
                  new_unittest("five waters cc-pVDZ, density falling 1e-12 along the cluster", test_w5_decay), &
                  new_unittest("three waters cc-pVTZ (f, LMAX=3 only)", test_w3_ccpvtz) &
                  ]
   end subroutine collect_mdj

   subroutine test_water_631gs(error)
      type(error_type), allocatable, intent(out) :: error
      call compare(error, "water.xyz", "basis_sets/6-31gs.json")
   end subroutine test_water_631gs

   subroutine test_water_ccpvdz(error)
      type(error_type), allocatable, intent(out) :: error
      call compare(error, "water.xyz", "basis_sets/cc-pvdz.json")
   end subroutine test_water_ccpvdz

   subroutine test_formamide_svp(error)
      type(error_type), allocatable, intent(out) :: error
      call compare(error, "formamide.xyz", "basis_sets/def2-svp.json")
   end subroutine test_formamide_svp

   subroutine test_w5_decay(error)
      type(error_type), allocatable, intent(out) :: error
      call compare(error, "w5.xyz", "basis_sets/cc-pvdz.json", decay=12.0_dp)
   end subroutine test_w5_decay

   subroutine test_w3_ccpvtz(error)
      type(error_type), allocatable, intent(out) :: error
      if (LMAX < 3) return
      call compare(error, "w3.xyz", "basis_sets/cc-pvtz.json")
   end subroutine test_w3_ccpvtz

   subroutine compare(error, xyz, basis, decay)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis
      real(dp), intent(in), optional :: decay
      !! w_i = 10^(-decay t_i/2), t_i in [0, 1] the AO centre's place along the
      !! axis of largest extent, so D_ij falls by up to 10^-decay
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat, nao, i, j, ish, ncart, ax
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d(:, :), jref(:, :), jmd(:, :), w(:)
      real(dp) :: lo, hi
      real(dp) :: scale, worst
      character(len=128) :: msg

      call read_xyz(xyz, nat, at_z, r)
      call trc_basis_from_json(basis, nat, at_z, r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      if (maxval(bas%sh_l) > LMAX) return

      nao = bas%nao
      allocate (d(nao, nao), jref(nao, nao), jmd(nao, nao))
      do j = 1, nao
         do i = 1, nao
            d(i, j) = cos(0.37_dp*real(i + j, dp)) + 0.5_dp*exp(-0.1_dp*real(abs(i - j), dp))
         end do
      end do
      if (present(decay)) then
         ax = maxloc(maxval(bas%sh_r, dim=2) - minval(bas%sh_r, dim=2), dim=1)
         lo = minval(bas%sh_r(ax, :))
         hi = maxval(bas%sh_r(ax, :))
         allocate (w(nao))
         do ish = 1, bas%nshell
            ncart = (bas%sh_l(ish) + 1)*(bas%sh_l(ish) + 2)/2
            w(bas%sh_ao(ish):bas%sh_ao(ish) + ncart - 1) = &
               10.0_dp**(-0.5_dp*decay*(bas%sh_r(ax, ish) - lo)/(hi - lo))
         end do
         do j = 1, nao
            do i = 1, nao
               d(i, j) = d(i, j)*w(i)*w(j)
            end do
         end do
      end if

      call bas%to_device()
      call eri%build(bas, THRESH)
      call eri%fock(bas, d, jref, density_screen=.false., j_scale=1.0_dp, k_scale=0.0_dp)
      call eri%coulomb(bas, d, jmd)
      call eri%release()
      call bas%release()

      scale = maxval(abs(jref))
      worst = maxval(abs(jmd - jref))/scale
      write (msg, '(a,es9.2,a,es9.2)') "max |J_md - J_fock|/max|J| = ", worst, ", max|J| = ", scale
      print '(4x,a)', trim(msg)
      call check(error, worst <= TOL, trim(msg))
   end subroutine compare

end module check_mdj_suite

program check_mdj
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_mdj_suite, only: collect_mdj
   implicit none
   integer :: stat
   stat = 0
   ! Serial: every case builds device-resident engines.
   call run_testsuite(collect_mdj, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, '(i0,a)') stat, " test(s) failed"
      error stop 1
   end if
end program check_mdj
