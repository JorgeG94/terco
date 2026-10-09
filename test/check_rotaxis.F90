!
! Rotated-axis vs HGP: the same Fock build, twice, through fock_bins with
! TRC_ROTAXIS forced off then on, compared directly.
!
! Every other check in this suite compares against something outside terco
! (libfint, pyscf) or an identity that must hold by construction. This one is
! neither -- it is the FOCUSED check milestone 1 of rotaxis asks for: the two
! algorithms are independent code paths computing the same physical
! quantity, so if they disagree by more than double-precision noise, one of
! them is wrong. `TRC_ROTAXIS` is read once per Fock build (trc_binkernel's
! `rotaxis_on()`), so it can be forced with libc's `setenv` between two calls
! in the SAME process -- no second binary, no ctest-level ENVIRONMENT dance,
! and no risk of comparing two different builds by accident.
!
! Two geometries, both 6-31G uncontracted to primitives (s and p only, so
! EVERY class in both runs goes through rotaxis when it is forced on -- there
! is no d shell here to fall back to HGP for):
!
!   water     -- the ordinary case, plus the A == B frame every atom's own
!                shell pairs give for free (see FRAME.md's degenerate list).
!   a linear molecule (CO2) -- collinear AB/CD frames throughout, which is
!                where a rotated-axis kernel's frame construction is under
!                the most stress (FRAME.md section 1's "AB parallel to CD"
!                and lab-z fallbacks) and a GPU kernel most easily produces a
!                silent NaN if a fallback is missing.
!
program check_rotaxis
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use trc_boys, only: dp, boys_init
   use trc_tables, only: tables_init, LMAX
   use trc_batch, only: build_pairs, ncart, common_fac_sp
   use trc_hgp, only: build_pairs_hgp
   use trc_cart, only: cart_init
   use trc_bins, only: pair_bins_t, build_binned_pairs
   use trc_binkernel, only: fock_bins
   use trc_screen, only: schwarz_bounds
   use trc_test_basis, only: read_xyz, uncontract_basis
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none

   interface
      function c_setenv(name, val, overwrite) bind(c, name="setenv") result(rc)
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: name(*), val(*)
         integer(c_int), value, intent(in) :: overwrite
         integer(c_int) :: rc
      end function c_setenv
   end interface

   integer :: nbad_total, rc
   nbad_total = 0

   call tables_init()
   call cart_init()
   call boys_init()

   if (LMAX < 1) then
      print '(a)', '  SKIPPED: build has LMAX < 1, cannot represent p functions'
      print '(a)', '  RESULT: PASS'
      stop
   end if

   call run_case('water.xyz')
   call run_case('rotaxis_linear.xyz')

   if (nbad_total > 0) then
      print '(a,i0)', '  quartets disagreeing beyond tolerance: ', nbad_total
      print '(a)', '  RESULT: FAIL'
      stop 1
   end if
   print '(a)', '  RESULT: PASS'

contains

   subroutine run_case(fn)
      character(len=*), intent(in) :: fn
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      integer :: nat, nbas, mnp, i, j, nlaunch2
      integer, allocatable :: l(:), np(:)
      real(dp), allocatable :: e(:, :), c(:, :), sr(:, :), cfac(:)
      integer, allocatable :: pp_off(:), pp_n(:)
      real(dp), allocatable :: pp_p(:), pp_r(:, :), pp_c(:), pp_e(:, :)
      real(dp), allocatable :: hp_p(:), hp_r(:, :), hp_ra(:, :), hp_rb(:, :), hp_c(:)
      integer, allocatable :: hp_off(:), hp_n(:), hp_ki(:), hp_kj(:)
      integer, allocatable :: ao_off(:)
      real(dp), allocatable :: qs(:), dsh(:, :)
      real(dp), allocatable :: dmat(:, :), jmat_hgp(:, :), kmat_hgp(:, :)
      real(dp), allocatable :: jmat_ra(:, :), kmat_ra(:, :)
      type(pair_bins_t) :: bins
      integer :: npp, nhpp, nao
      integer(kind=8) :: nwork2
      real(dp) :: worstj, worstk, maxf
      type(trc_basis_t) :: jbas
      type(error_t) :: err

      call read_xyz(fn, nat, z, r)
      ! 6-31G, through the one basis-construction path, then uncontracted
      ! to primitives exactly as build_631g's old `uncontracted=.true.`
      ! argument used to (see uncontract_basis's own comment).
      call trc_basis_from_json('basis_sets/6-31g.json', nat, z, r, jbas, err)
      if (err%has_error()) then
         print '(a)', 'check_rotaxis: '//err%get_message(); stop 1
      end if
      call uncontract_basis(jbas, nbas, l, np, e, c, sr, mnp)
      if (maxval(l) > LMAX) then
         print '(a,a,a)', '  ', trim(fn), ': SKIPPED (needs p, build has LMAX < 1)'
         return
      end if

      allocate (ao_off(nbas), cfac(nbas))
      nao = 0
      do i = 1, nbas
         ao_off(i) = nao + 1
         nao = nao + ncart(l(i))
         cfac(i) = common_fac_sp(l(i))
      end do

      call build_pairs(nbas, l, np, e, c, sr, 1.0e-30_dp, &
                       pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, npp)
      allocate (qs(nbas*(nbas + 1)/2))
      call schwarz_bounds(nbas, npp, l, pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, qs)
      call build_pairs_hgp(nbas, l, np, e, c, sr, cfac, c, 1.0e-33_dp, &
                           hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, hp_ki, hp_kj, nhpp)
      call build_binned_pairs(nbas, l, np, sr, qs, 1.0e-30_dp, bins, pp_n=hp_n)

      ! No screening: the comparison should hold for every quartet either
      ! path would ever evaluate, not just the ones a shared screen admits
      ! on both sides identically (which would agree trivially at zero).
      allocate (dsh(nbas, nbas))
      dsh = huge(1.0_dp)*1.0e-30_dp

      allocate (dmat(nao, nao))
      allocate (jmat_hgp(nao, nao), kmat_hgp(nao, nao))
      allocate (jmat_ra(nao, nao), kmat_ra(nao, nao))
      do i = 1, nao
         do j = 1, nao
            dmat(j, i) = 0.4_dp + 0.1_dp*cos(real(3*j + 5*i, dp))
         end do
      end do
      dmat = 0.5_dp*(dmat + transpose(dmat))

      !$acc enter data copyin(l, ao_off, dmat, dsh, bins, bins%sp_i, bins%sp_j, bins%sp_q, &
      !$acc                   hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c)

      rc = c_setenv("TRC_ROTAXIS"//c_null_char, "0"//c_null_char, 1_c_int)
      jmat_hgp = 0.0_dp; kmat_hgp = 0.0_dp
      !$acc enter data copyin(jmat_hgp, kmat_hgp)
      call fock_bins(bins, nbas, nhpp, nao, l, ao_off, 1.0e-30_dp, .false., &
                     1.0_dp, 1.0_dp, .false., dsh, &
                     hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
                     1, dmat, jmat_hgp, kmat_hgp, 0, 1, nlaunch2, nwork2)
      !$acc update self(jmat_hgp, kmat_hgp)
      !$acc exit data delete(jmat_hgp, kmat_hgp)

      rc = c_setenv("TRC_ROTAXIS"//c_null_char, "1"//c_null_char, 1_c_int)
      jmat_ra = 0.0_dp; kmat_ra = 0.0_dp
      !$acc enter data copyin(jmat_ra, kmat_ra)
      call fock_bins(bins, nbas, nhpp, nao, l, ao_off, 1.0e-30_dp, .false., &
                     1.0_dp, 1.0_dp, .false., dsh, &
                     hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
                     1, dmat, jmat_ra, kmat_ra, 0, 1, nlaunch2, nwork2)
      !$acc update self(jmat_ra, kmat_ra)
      !$acc exit data delete(jmat_ra, kmat_ra)

      rc = c_setenv("TRC_ROTAXIS"//c_null_char, "0"//c_null_char, 1_c_int)

      !$acc exit data delete(l, ao_off, dmat, dsh, bins%sp_i, bins%sp_j, bins%sp_q, bins, &
      !$acc                  hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c)

      ! Under TRC_FOCK6 (this build's TERCO_TEST_DEFS), the per-class
      ! dispatchers fold BOTH the Coulomb and exchange contributions into
      ! `jmat` -- `kmat` is never written on this path (see trc_binkernel's
      ! pc_dispatch/ra_dispatch calls, which pass only `jmat`) -- so the K
      ! comparison below would always read 0.00E+00 and look like a second,
      ! independent pass when it is really the same J comparison twice. Only
      ! J is reported; the comment says why rather than printing a vacuous
      ! column.
      maxf = max(maxval(abs(jmat_hgp)), maxval(abs(kmat_hgp)), 1.0e-300_dp)
      call cmp(jmat_hgp, jmat_ra, nao, maxf, worstj)
      call cmp(kmat_hgp, kmat_ra, nao, maxf, worstk)
      print '(a,a,a,i0,a,es10.2,a,es10.2,a)', '  ', trim(fn), ': nao ', nao, &
         '  max|F| ', maxf, '  worst J+K diff (scaled) ', worstj/maxf, &
         '  [K folded into J under TRC_FOCK6, kmat unused on this path]'

      deallocate (l, np, e, c, sr, ao_off, cfac)
      deallocate (pp_off, pp_n, pp_p, pp_r, pp_c, pp_e)
      deallocate (hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, hp_ki, hp_kj)
      deallocate (dmat, jmat_hgp, kmat_hgp, jmat_ra, kmat_ra, dsh, qs)
   end subroutine run_case

   subroutine cmp(a, b, nao, maxf, worst)
      real(dp), intent(in) :: a(:, :), b(:, :), maxf
      integer, intent(in) :: nao
      real(dp), intent(inout) :: worst
      ! About 1e-11 relative to max|F|, per the task: rotaxis and HGP are two
      ! different floating-point paths to the same number, not the same
      ! arithmetic in a different order, so this is a genuine (if small)
      ! disagreement budget rather than atomic-update noise.
      real(dp), parameter :: RTOL = 1.0e-11_dp
      integer :: i, j
      real(dp) :: ad
      do j = 1, nao
         do i = 1, nao
            ad = abs(a(i, j) - b(i, j))
            if (ad > worst) worst = ad
            if (ad > RTOL*maxf) then
               nbad_total = nbad_total + 1
               if (nbad_total <= 8) print '(a,2(1x,i0),2es22.13)', '   MISS', i, j, a(i, j), b(i, j)
            end if
         end do
      end do
   end subroutine cmp

end program check_rotaxis
