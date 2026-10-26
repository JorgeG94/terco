!
! Rotated-axis vs HGP: the same Fock build, twice, through fock_bins with
! TRC_ROTAXIS forced off then on, compared directly.
!
! Every other check in this suite compares against something outside terco
! (libfint, pyscf) or an identity that must hold by construction. This one is
! neither -- it is the FOCUSED check the revived rotaxis generator asks for:
! the two algorithms are independent code paths computing the same physical
! quantity (every class here is s/p, so TRC_ROTAXIS=1 sends ALL of it through
! `ra_dispatch`), so if they disagree by more than double-precision noise,
! one of them is wrong. `TRC_ROTAXIS` is read once per Fock build
! (trc_binkernel's `rotaxis_on()`), so it can be forced with libc's `setenv`
! between two calls in the SAME process -- no second binary, no ctest-level
! ENVIRONMENT dance, and no risk of comparing two different builds by
! accident.
!
! Two systems, chosen for contraction degree rather than geometry (FRAME.md's
! degenerate-frame cases are exercised plenty by check_binfock's own
! systems):
!
!   SiH4 / cc-pVDZ   -- light contraction (Si s/p K ~ 4-6); cc-pVDZ's d
!                       shells on Si stay on HGP regardless of the knob
!                       (`ra_has_class` has no case above p), so this ALSO
!                       checks that a mixed s/p/d run routes correctly.
!   SiH4 / cc-pVTZ (H/O/Si only, cc-pvtz_sioh.json) -- the deep general
!                       contraction this generator was revived to measure
!                       (Si s K=13, p K=7); d/f shells on Si likewise stay
!                       on HGP/Rys. Kept to 5 atoms (not Si8O12H8) because
!                       this check turns off every screen (see dsh below),
!                       which an f-containing basis on a bigger system
!                       makes expensive fast; the deep contraction itself
!                       only needs one Si, not a slab of them.
!
program check_rotaxis
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use pic_types, only: dp, int64
   use trc_boys, only: boys_init
   use trc_tables, only: tables_init, LMAX
   use trc_batch, only: build_pairs, ncart, common_fac_sp
   use trc_hgp, only: build_pairs_hgp
   use trc_cart, only: cart_init
   use trc_bins, only: pair_bins_t, build_binned_pairs, build_pair_quartets
   use trc_binkernel, only: fock_bins
   use trc_screen, only: schwarz_bounds
   use trc_test_basis, only: read_xyz
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

   call run_case('sih4.xyz', 'basis_sets/cc-pvdz.json')
   call run_case('sih4.xyz', 'basis_sets/cc-pvtz_sioh.json')

   if (nbad_total > 0) then
      print '(a,i0)', '  AO pairs disagreeing beyond tolerance: ', nbad_total
      print '(a)', '  RESULT: FAIL'
      stop 1
   end if
   print '(a)', '  RESULT: PASS'

contains

   subroutine run_case(fn, basfn)
      character(len=*), intent(in) :: fn, basfn
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      integer :: nat, nbas, i, j, nlaunch
      integer, allocatable :: sh_l(:), sh_np(:), ao_off(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), cfac(:)
      integer, allocatable :: pp_off(:), pp_n(:)
      real(dp), allocatable :: pp_p(:), pp_r(:, :), pp_c(:), pp_e(:, :)
      real(dp), allocatable :: hp_p(:), hp_r(:, :), hp_ra(:, :), hp_rb(:, :), hp_c(:)
      integer, allocatable :: hp_off(:), hp_n(:), hp_ki(:), hp_kj(:)
      integer, allocatable :: qp_bin_pbase(:), qp_bin_k(:)
      integer :: nqp
      real(dp), allocatable :: qp_z(:), qp_c(:), qp_px(:), qp_py(:), qp_pz(:)
      real(dp), allocatable :: sq_ax(:), sq_ay(:), sq_az(:)
      real(dp), allocatable :: sq_abx(:), sq_aby(:), sq_abz(:)
      integer, allocatable :: sq_n(:)
      real(dp), allocatable :: qs(:), dsh(:, :)
      real(dp), allocatable :: dmat(:, :), jmat_hgp(:, :), kmat_hgp(:, :)
      real(dp), allocatable :: jmat_ra(:, :), kmat_ra(:, :)
      type(pair_bins_t) :: bins
      integer :: npp, nhpp, nao, maxnp
      integer(int64) :: nwork
      real(dp) :: worstj, worstk, maxf
      type(trc_basis_t) :: jbas
      type(error_t) :: err

      call read_xyz(fn, nat, z, r)
      call trc_basis_from_json(basfn, nat, z, r, jbas, err)
      if (err%has_error()) then
         print '(a)', 'check_rotaxis: '//err%get_message(); stop 1
      end if

      nbas = jbas%nshell
      maxnp = jbas%maxnp
      allocate (sh_l(nbas), sh_np(nbas), sh_r(3, nbas), ao_off(nbas), cfac(nbas))
      allocate (sh_e(maxnp, nbas), sh_c(maxnp, nbas))
      sh_l = jbas%sh_l; sh_np = jbas%sh_np; sh_r = jbas%sh_r
      sh_e = jbas%sh_e; sh_c = jbas%sh_c

      if (maxval(sh_l) > LMAX) then
         print '(a,a,a)', '  ', trim(fn), ': SKIPPED (basis needs angular momentum above this build''s LMAX)'
         return
      end if

      nao = 0
      do i = 1, nbas
         ao_off(i) = nao + 1
         nao = nao + ncart(sh_l(i))
         cfac(i) = common_fac_sp(sh_l(i))
      end do

      call build_pairs(nbas, sh_l, sh_np, sh_e, sh_c, sh_r, 1.0e-30_dp, &
                       pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, npp)
      allocate (qs(nbas*(nbas + 1)/2))
      call schwarz_bounds(nbas, npp, sh_l, pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, qs)
      call build_pairs_hgp(nbas, sh_l, sh_np, sh_e, sh_c, sh_r, cfac, sh_c, 1.0e-33_dp, &
                           hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, hp_ki, hp_kj, nhpp)
      call build_binned_pairs(nbas, sh_l, sh_np, sh_r, qs, 1.0e-30_dp, bins, pp_n=hp_n)
      call build_pair_quartets(bins, nbas, nhpp, hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
                               qp_bin_pbase, qp_bin_k, nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                               sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n)

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

      !$acc enter data copyin(sh_l, ao_off, dmat, dsh, bins, bins%sp_i, bins%sp_j, bins%sp_q, &
      !$acc                   hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
      !$acc                   qp_z, qp_c, qp_px, qp_py, qp_pz, &
      !$acc                   sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n)

      rc = c_setenv("TRC_ROTAXIS"//c_null_char, "0"//c_null_char, 1_c_int)
      jmat_hgp = 0.0_dp; kmat_hgp = 0.0_dp
      !$acc enter data copyin(jmat_hgp, kmat_hgp)
      call fock_bins(bins, nbas, nhpp, nao, sh_l, ao_off, 1.0e-30_dp, .false., &
                     1.0_dp, 1.0_dp, 0.0_dp, .false., dsh, &
                     hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
                     qp_bin_pbase, qp_bin_k, nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                     sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                     1, dmat, jmat_hgp, kmat_hgp, 0, 1, nlaunch, nwork)
      !$acc update self(jmat_hgp, kmat_hgp)
      !$acc exit data delete(jmat_hgp, kmat_hgp)

      rc = c_setenv("TRC_ROTAXIS"//c_null_char, "1"//c_null_char, 1_c_int)
      jmat_ra = 0.0_dp; kmat_ra = 0.0_dp
      !$acc enter data copyin(jmat_ra, kmat_ra)
      call fock_bins(bins, nbas, nhpp, nao, sh_l, ao_off, 1.0e-30_dp, .false., &
                     1.0_dp, 1.0_dp, 0.0_dp, .false., dsh, &
                     hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
                     qp_bin_pbase, qp_bin_k, nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                     sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                     1, dmat, jmat_ra, kmat_ra, 0, 1, nlaunch, nwork)
      !$acc update self(jmat_ra, kmat_ra)
      !$acc exit data delete(jmat_ra, kmat_ra)

      rc = c_setenv("TRC_ROTAXIS"//c_null_char, "0"//c_null_char, 1_c_int)

      !$acc exit data delete(sh_l, ao_off, dmat, dsh, bins%sp_i, bins%sp_j, bins%sp_q, bins, &
      !$acc                  hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
      !$acc                  qp_z, qp_c, qp_px, qp_py, qp_pz, &
      !$acc                  sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n)

      maxf = max(maxval(abs(jmat_hgp)), maxval(abs(kmat_hgp)), 1.0e-300_dp)
      call cmp(jmat_hgp, jmat_ra, nao, maxf, worstj)
      call cmp(kmat_hgp, kmat_ra, nao, maxf, worstk)
      print '(a,a,a,i0,a,i0,a,es10.2,a,es10.2,a,es10.2)', '  ', trim(fn), ': nbas ', nbas, &
         '  nao ', nao, '  max|F| ', maxf, '  worst J diff (scaled) ', worstj/maxf, &
         '  worst K diff (scaled) ', worstk/maxf

      deallocate (sh_l, sh_np, sh_e, sh_c, sh_r, ao_off, cfac)
      deallocate (pp_off, pp_n, pp_p, pp_r, pp_c, pp_e)
      deallocate (hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, hp_ki, hp_kj)
      deallocate (dmat, jmat_hgp, kmat_hgp, jmat_ra, kmat_ra, dsh, qs)
   end subroutine run_case

   subroutine cmp(a, b, nao, maxf, worst)
      real(dp), intent(in) :: a(:, :), b(:, :), maxf
      integer, intent(in) :: nao
      real(dp), intent(inout) :: worst
      ! Rotaxis and HGP are two different floating-point paths to the same
      ! number, not the same arithmetic in a different order, so this is a
      ! genuine (if small) disagreement budget rather than atomic-update
      ! noise.
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
