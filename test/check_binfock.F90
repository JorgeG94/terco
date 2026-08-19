!
! Fused Fock build against a host contraction of libfint's own integrals.
!
!   J(mu,nu) = sum_{la,si} (mu nu|la si) D(la,si)
!   K(mu,la) = sum_{nu,si} (mu nu|la si) D(nu,si)
!
! The reference side reads reference.dat and contracts it directly, so this
! tests the DIGESTION and nothing else -- the integrals themselves are already
! pinned by check_batch.  If the integrals were wrong both sides would move
! together and the test would pass, which is why it is a separate check run
! after that one rather than instead of it.
!
! The density is deterministic and dense (no zeros to hide an indexing slip),
! asymmetric-looking in its generator but symmetrised, since J and K assume a
! symmetric D.
!
! Tolerance is looser than the integral check on purpose: the device path
! accumulates through `!$acc atomic update`, so the FP addition order varies
! run to run and the result is not bit-reproducible.  That is the documented
! cost of atomics, not a defect.
!
program check_binfock
   use trc_boys, only: dp, boys_init
   use trc_tables, only: LTOT, NHERM_MAX, tables_init
   use trc_batch, only: build_pairs, ncart, common_fac_sp
   use trc_cart, only: cart_init
   use trc_hgp, only: build_pairs_hgp
   use trc_tables, only: LMAX
   use trc_bins, only: pair_bins_t, build_binned_pairs
   use trc_binkernel, only: fock_bins
   use trc_screen, only: schwarz_bounds, build_worklist, density_blockmax, pair_index
   implicit none

   integer, parameter :: ATM_SLOTS = 6, BAS_SLOTS = 8
   integer, parameter :: ATOM_OF = 0, ANG_OF = 1, NPRIM_OF = 2
   integer, parameter :: PTR_EXP = 5, PTR_COEFF = 6, PTR_COORD = 1

   integer, allocatable :: atm(:), bas(:)
   real(dp), allocatable :: env(:)
   integer :: natm, nbas, nenv, u, isys
   character(len=64) :: tag
   real(dp) :: worstj, worstk
   integer  :: nbad
   real(dp), parameter :: ATOL = 1.0e-11_dp, RTOL = 1.0e-11_dp

   worstj = 0.0_dp; worstk = 0.0_dp; nbad = 0

   call tables_init()
   call cart_init()
   call boys_init()

   open (newunit=u, file='reference.dat', status='old', action='read')
   do
      read (u, '(a)') tag
      if (trim(tag) == 'END') exit
      if (tag(1:6) /= 'SYSTEM') cycle
      read (tag(8:), *) isys
      call read_system(u, atm, bas, env, natm, nbas, nenv)
      call run_system(u, isys, atm, bas, env, nbas)
      deallocate (atm, bas, env)
   end do
   close (u)

   print '(a,es10.2)', '  worst J diff    : ', worstj
   print '(a,es10.2)', '  worst K diff    : ', worstk
   print '(a,i0)',     '  outside tol     : ', nbad
   if (nbad > 0) then
      print '(a)', '  RESULT: FAIL'
      stop 1
   end if
   print '(a)', '  RESULT: PASS'

contains

   subroutine read_system(u, atm, bas, env, natm, nbas, nenv)
      integer, intent(in) :: u
      integer, allocatable, intent(out) :: atm(:), bas(:)
      real(dp), allocatable, intent(out) :: env(:)
      integer, intent(out) :: natm, nbas, nenv
      character(len=64) :: t
      integer :: i
      read (u, '(a)') t; read (t(6:), *) natm
      read (u, '(a)') t; read (t(6:), *) nbas
      read (u, '(a)') t; read (t(6:), *) nenv
      allocate (atm(0:ATM_SLOTS*natm - 1), bas(0:BAS_SLOTS*nbas - 1), env(0:nenv - 1))
      read (u, '(a)') t; read (u, *) (atm(i), i=0, ATM_SLOTS*natm - 1)
      read (u, '(a)') t; read (u, *) (bas(i), i=0, BAS_SLOTS*nbas - 1)
      read (u, '(a)') t; read (u, *) (env(i), i=0, nenv - 1)
   end subroutine read_system

   subroutine skip_to_next(u)
      integer, intent(in) :: u
      character(len=64) :: t
      do
         read (u, '(a)') t
         if (trim(t) == 'END' .or. t(1:6) == 'SYSTEM') then
            backspace (u); return
         end if
      end do
   end subroutine skip_to_next

   subroutine run_system(u, isys, atm, bas, env, nbas)
      integer, intent(in) :: u, isys, nbas
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)

      integer, allocatable :: sh_l(:), sh_np(:), ao_off(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
      integer, allocatable :: pp_off(:), pp_n(:)
      real(dp), allocatable :: pp_p(:), pp_r(:, :), pp_c(:), pp_e(:, :)
      integer, allocatable :: q_i(:), q_j(:), q_k(:), q_l(:), q_off(:)
      real(dp), allocatable :: out(:), rscr(:, :)
      real(dp), allocatable :: dmat(:, :), jmat(:, :), kmat(:, :), dsh(:, :)
      real(dp), allocatable :: jref(:, :), kref(:, :), ref(:)
      real(dp), allocatable :: qs(:), dmax(:, :), cfac(:)
      real(dp), allocatable :: hp_p(:), hp_r(:,:), hp_ra(:,:), hp_rb(:,:), hp_c(:)
      integer, allocatable :: hp_off(:), hp_n(:)
      type(pair_bins_t) :: bins
      integer :: nhpp, nlaunch
      integer(kind=8) :: nwork
      integer :: npp, maxnp, nq, nout, nao, i, j, k, l, iat, iq, m
      integer :: sh(4), d(4), nc4, nqrt
      integer :: na, nb, nc, nd, ia, ib, ic, id, idx
      integer :: mu, nu, lam, sig
      character(len=64) :: t
      real(dp) :: dd, v

      allocate (sh_l(nbas), sh_np(nbas), sh_r(3, nbas), ao_off(nbas), cfac(nbas))
      maxnp = 0; nao = 0
      do i = 1, nbas
         sh_l(i) = bas(BAS_SLOTS*(i - 1) + ANG_OF)
         sh_np(i) = bas(BAS_SLOTS*(i - 1) + NPRIM_OF)
         maxnp = max(maxnp, sh_np(i))
         iat = bas(BAS_SLOTS*(i - 1) + ATOM_OF)
         sh_r(:, i) = env(atm(ATM_SLOTS*iat + PTR_COORD):atm(ATM_SLOTS*iat + PTR_COORD) + 2)
         ao_off(i) = nao + 1
         nao = nao + ncart(sh_l(i))
      end do
      do i = 1, nbas
         cfac(i) = common_fac_sp(sh_l(i))
      end do
      allocate (sh_e(maxnp, nbas), sh_c(maxnp, nbas))
      sh_e = 0.0_dp; sh_c = 0.0_dp
      do i = 1, nbas
         k = bas(BAS_SLOTS*(i - 1) + PTR_EXP)
         sh_e(1:sh_np(i), i) = env(k:k + sh_np(i) - 1)
         k = bas(BAS_SLOTS*(i - 1) + PTR_COEFF)
         sh_c(1:sh_np(i), i) = env(k:k + sh_np(i) - 1)
      end do

      ! dense, deterministic, symmetric density
      allocate (dmat(nao, nao), jmat(nao, nao), kmat(nao, nao))
      allocate (jref(nao, nao), kref(nao, nao))
      do i = 1, nao
         do j = 1, nao
            dd = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j, dp))
            dmat(i, j) = dd
         end do
      end do
      dmat = 0.5_dp*(dmat + transpose(dmat))
      jmat = 0.0_dp; kmat = 0.0_dp; jref = 0.0_dp; kref = 0.0_dp

      ! A build with a lower LMAX cannot represent every reference system.
      ! Skip rather than compute nonsense, and say so.
      if (maxval(sh_l) > LMAX) then
         print '(a,i0,a,i0,a,i0,a)', '  system ', isys, ': SKIPPED (needs l=', &
            maxval(sh_l), ', build has LMAX=', LMAX, ')'
         call skip_to_next(u)
         deallocate (sh_l, sh_np, sh_e, sh_c, sh_r, cfac)
         deallocate (dmat, jmat, kmat, jref, kref)
         return
      end if

      call build_pairs(nbas, sh_l, sh_np, sh_e, sh_c, sh_r, 1.0e-30_dp, &
                       pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, npp)

      ! Schwarz bounds, then the binned pair container.  No quartet list.
      allocate (qs(nbas*(nbas + 1)/2), dmax(nbas, nbas))
      call schwarz_bounds(nbas, npp, sh_l, pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, qs)
      call build_binned_pairs(nbas, sh_l, sh_np, sh_r, qs, 1.0e-30_dp, bins)
      call build_pairs_hgp(nbas, sh_l, sh_np, sh_e, sh_c, sh_r, cfac, &
                           hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, nhpp)
      print '(a,i0,a,i0,a,i0)', '  system ', isys, ': shell pairs ', bins%npair, &
         '   live bins ', bins%nlive

      ! The correctness oracle must screen NOTHING, so the density bound is
      ! set to a value no product can fall below.
      allocate (dsh(nbas, nbas))
      dsh = huge(1.0_dp)*1.0e-30_dp

      !$acc enter data copyin(sh_l, ao_off, dmat, dsh, bins, bins%sp_i, bins%sp_j, bins%sp_q, &
      !$acc                   hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c) &
      !$acc            copyin(jmat, kmat)
      call fock_bins(bins, nbas, nhpp, nao, sh_l, ao_off, 1.0e-30_dp, .false., &
                     1.0_dp, 1.0_dp, .false., dsh, &
                     hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
                     1, dmat, jmat, kmat, nlaunch, nwork)
      !$acc update self(jmat, kmat)
      !$acc exit data delete(sh_l, ao_off, dmat, dsh, bins%sp_i, bins%sp_j, bins%sp_q, bins, &
      !$acc                  hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, jmat, kmat)
      print '(a,i0,a,i0)', '           launches ', nlaunch, '   work items ', nwork

      ! --- reference: contract libfint's own integrals on the host ---
      do
         read (u, '(a)') t
         if (trim(t) == 'END' .or. t(1:6) == 'SYSTEM') then
            backspace (u); exit
         end if
         if (t(1:7) /= 'QUARTET') cycle
         read (t(8:), *) sh(1), sh(2), sh(3), sh(4), d(1), d(2), d(3), d(4)
         nc4 = d(1)*d(2)*d(3)*d(4)
         allocate (ref(nc4))
         read (u, *) (ref(m), m=1, nc4)
         na = d(1); nb = d(2); nc = d(3); nd = d(4)
         idx = 0
         do id = 0, nd - 1
            sig = ao_off(sh(4) + 1) + id
         do ic = 0, nc - 1
            lam = ao_off(sh(3) + 1) + ic
         do ib = 0, nb - 1
            nu = ao_off(sh(2) + 1) + ib
         do ia = 0, na - 1
            mu = ao_off(sh(1) + 1) + ia
            idx = idx + 1
            v = ref(idx)
            jref(mu, nu) = jref(mu, nu) + v*dmat(lam, sig)
            kref(mu, lam) = kref(mu, lam) + v*dmat(nu, sig)
         end do
         end do
         end do
         end do
         deallocate (ref)
      end do

#if defined(TRC_FOCK6)
      ! The folded path builds G = 2J - K in symmetric storage, so symmetrise
      ! and compare against the same combination of the SAME libfint integrals.
      ! The oracle is unchanged -- only which contraction of it we check.
      ! G was accumulated into one triangle-ish half: each unordered AO pair
      ! written once.  Symmetrising COPIES a write to its transpose, it does not
      ! double it -- so G_sym = G + G^T, and the diagonal, which was written in
      ! place, has to have its double-count removed.
      block
         real(dp), allocatable :: gs(:, :)
         allocate (gs(nao, nao))
         ! Plain G + G^T, diagonal included.  Correcting the diagonal for a
         ! double count was wrong and showed up as every diagonal element being
         ! exactly half the reference: the 0.5 degeneracy scales already put
         ! half the value there, and the transpose supplies the other half.
         gs = jmat + transpose(jmat)
         jref = 2.0_dp*jref - kref
         call cmp(gs, jref, nao, worstj)
         deallocate (gs)
      end block
#else
      call cmp(jmat, jref, nao, worstj)
      call cmp(kmat, kref, nao, worstk)
#endif

      deallocate (sh_l, sh_np, sh_e, sh_c, sh_r, ao_off)
      deallocate (pp_off, pp_n, pp_p, pp_r, pp_c, pp_e)
      deallocate (dmat, jmat, kmat, jref, kref, qs, dmax, cfac)
   end subroutine run_system

   subroutine cmp(got, ref, nao, worst)
      real(dp), intent(in) :: got(:, :), ref(:, :)
      integer, intent(in) :: nao
      real(dp), intent(inout) :: worst
      integer :: i, j
      real(dp) :: ad
      do j = 1, nao
         do i = 1, nao
            ad = abs(got(i, j) - ref(i, j))
            if (ad > worst) worst = ad
            if (ad > ATOL + RTOL*abs(ref(i, j))) then
               nbad = nbad + 1
               if (nbad <= 8) print '(a,2(1x,i0),a,l1,2es22.13)', '   MISS', i, j, &
                  '  diag=', i == j, ref(i, j), got(i, j)
            end if
         end do
      end do
   end subroutine cmp

end program check_binfock
