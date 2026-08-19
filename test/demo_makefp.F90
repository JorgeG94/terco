!
! makeEFP, driven end to end by terco.
!
! WHAT THIS IS
! ------------
! A standalone copy of mqc's `make_efp_potential` pipeline, with every integral
! coming from terco instead of libcint. It is a DEMO, not the backend: the
! backend (Phase 4 of PLAN_MAKEEFP.md) lives inside mqc and has to be built
! against mqc's own molecule type, its logger and its error handling. This
! program exists so the pathway can be run and shown before that wiring is
! possible, on a machine that has a GPU but not a full mqc build.
!
! It follows the real routine's stage order, which is forced by what depends on
! what:
!
!     1. basis, pair list, ERI object          terco
!     2. RHF SCF                               terco, four-centre
!     3. Boys localization                     host, on terco's dipole matrices
!     4. distributed multipoles                terco, through octupole
!     5. static polarizability by CPHF         terco, BATCHED over perturbations
!     6. write the potential
!
! The SCF gives the density the multipoles are taken from and the orbitals
! everything else needs; localization gives the centres the polarizabilities
! are expressed on; the multipoles have to exist before any screening could be
! fitted. That ordering is copied deliberately.
!
! WHAT IS FAITHFUL AND WHAT IS NOT
! ---------------------------------
! Faithful: the expansion points (atoms plus bond midpoints by covalent radii),
! the nearest-point partition of the primitive density with ties split evenly,
! the packed component orders a `.efp` carries, the Boys functional, and the
! CPHF equations.
!
! Not attempted: the charge-penetration screening fit, the dynamic (Casimir-
! Polder) blocks, the exchange-repulsion basis projection, and GAMESS's AO
! ordering for the orbital sections. Those are linear algebra and file format,
! not integrals -- they are exactly the parts that stay in mqc under terco's
! contract, and copying them here would prove nothing about terco.
!
! WHY STAGE 5 IS THE POINT
! ------------------------
! Everything above stage 5 is one pass. Stage 5 is the one that repeats: the
! real makeEFP solves about a hundred right-hand sides (nine perturbations at
! twelve imaginary frequencies), and each iteration of each is a Fock build.
! This demo runs the static block only -- three perturbations -- but runs them
! THROUGH `fock_many`, so the batching measured in Phase 2 is what is being
! exercised, and the printed per-iteration cost is the number that scales to
! the full potential.
!
!   demo_makefp <xyz> [charge] [thresh] [ncore]
!
program demo_makefp
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_multipoles, &
                      TRC_NMULT
   use trc_fock, only: trc_eri_t
   use trc_test_basis, only: read_xyz, build_631g
   implicit none

   real(dp), parameter :: BOHR_PER_ANG = 1.8897261254578281_dp
   integer,  parameter :: MAXEL = 36
   real(dp), parameter :: COVALENT(MAXEL) = [ &
      0.32_dp, 0.93_dp, &
      1.23_dp, 0.90_dp, 0.82_dp, 0.77_dp, 0.75_dp, 0.73_dp, 0.72_dp, 0.71_dp, &
      1.54_dp, 1.36_dp, 1.18_dp, 1.11_dp, 1.06_dp, 1.02_dp, 0.99_dp, 0.98_dp, &
      2.03_dp, 1.74_dp, 1.44_dp, 1.32_dp, 1.22_dp, 1.18_dp, 1.17_dp, 1.17_dp, &
      1.16_dp, 1.15_dp, 1.17_dp, 1.25_dp, 1.26_dp, 1.22_dp, 1.20_dp, 1.16_dp, &
      1.14_dp, 1.12_dp]
   !> Offsets into the full Cartesian tensors terco returns. The dipole is
   !> components 1:3, the quadrupole 4:12 (all nine) and the octupole 13:39
   !> (all twenty-seven); a `.efp` carries only the unique ones, in this order.
   integer, parameter :: QUAD_PACK(6)  = [1, 5, 9, 2, 3, 6]
   integer, parameter :: OCT_PACK(10)  = [1, 14, 27, 2, 3, 5, 15, 9, 18, 6]
   real(dp), parameter :: TIE_TOL = 1.0e-6_dp

   character(len=256) :: xyzfile, arg
   integer  :: charge, ncore_in
   real(dp) :: thresh

   ! --- the molecule and its contracted basis --------------------------------
   integer :: nat, nsh, maxnp, nao, nocc, nvirt, nval, ncore
   integer,  allocatable :: at_z(:), sh_l(:), sh_np(:)
   real(dp), allocatable :: at_r(:, :), zat(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
   type(trc_basis_t)    :: bas
   type(trc_pairlist_t) :: pl
   type(trc_eri_t)      :: eri

   ! --- one-electron and the SCF ---------------------------------------------
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :), hcore(:, :)
   real(dp), allocatable :: mm(:, :, :), dmat(:, :), gmat(:, :), fock(:, :)
   real(dp), allocatable :: x(:, :), cmo(:, :), eps(:)
   real(dp) :: enuc, escf

   ! --- localization ----------------------------------------------------------
   real(dp), allocatable :: loc(:, :), qrot(:, :), centroid(:, :)

   ! --- distributed multipoles ------------------------------------------------
   integer :: npoint
   character(len=8), allocatable :: plabel(:)
   real(dp), allocatable :: ppos(:, :), qnuc(:), qel(:)
   real(dp), allocatable :: pdip(:, :), pquad(:, :), poct(:, :)

   ! --- polarizability --------------------------------------------------------
   real(dp) :: alpha(3, 3)
   real(dp), allocatable :: alpha_lmo(:, :, :)

   real(dp) :: t0, t1, tstage
   integer  :: i, j, k

   ! ==========================================================================
   charge = 0; thresh = 1.0e-10_dp; ncore_in = -1
   if (command_argument_count() < 1) then
      print '(a)', 'usage: demo_makefp <xyz> [charge] [thresh] [ncore]'
      stop 1
   end if
   call get_command_argument(1, xyzfile)
   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg); read (arg, *) charge
   end if
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg); read (arg, *) thresh
   end if
   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg); read (arg, *) ncore_in
   end if

   print '(a)', '=================================================='
   print '(a)', ' makeEFP through terco'
   print '(a)', '=================================================='

   ! --- stage 1: the basis and terco's containers ----------------------------
   call tick(t0)
   call read_xyz(xyzfile, nat, at_z, at_r)
   call build_631g(nat, at_z, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   allocate (zat(nat)); zat = real(at_z, dp)
   call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, at_r, maxnp)
   call pl%build(bas, thresh)
   nao = bas%nao
   nocc = (sum(at_z) - charge)/2
   if (2*nocc /= sum(at_z) - charge) then
      print '(a)', '  odd electron count: this demo is closed-shell only'
      stop 1
   end if
   nvirt = nao - nocc
   ncore = frozen_core(nat, at_z)
   if (ncore_in >= 0) ncore = ncore_in
   nval = nocc - ncore

   call bas%to_device()
   call pl%to_device()
   call eri%build(bas, thresh)
   call tick(t1)
   print '(a,a)',        '  geometry   ', trim(xyzfile)
   print '(a,i0,a,i0,a,i0)', '  6-31G      ', nao, ' functions, ', nsh, &
      ' shells, atoms ', nat
   print '(a,i0,a,i0,a,i0)', '  occupied   ', nocc, '   core ', ncore, &
      '   valence (localized) ', nval
   print '(a,es9.2)',    '  threshold  ', thresh
   print '(a,f9.3,a)',   '  stage 1    ', t1 - t0, ' s   basis, pair list, Schwarz'

   ! --- stage 2: the SCF -----------------------------------------------------
   allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao), hcore(nao, nao))
   allocate (mm(nao, nao, TRC_NMULT))
   !$acc enter data create(smat, tmat, vmat, mm)
   call trc_1e(bas, pl, smat, tmat, vmat)
   !$acc update self(smat, tmat, vmat)
   hcore = tmat + vmat
   enuc = nuclear_repulsion(nat, zat, at_r)

   allocate (x(nao, nao), cmo(nao, nao), eps(nao))
   allocate (dmat(nao, nao), gmat(nao, nao), fock(nao, nao))
   call tick(t0)
   call run_scf()
   call tick(t1)
   print '(a,f9.3,a,f18.10)', '  stage 2    ', t1 - t0, ' s   RHF energy ', escf

   ! --- stage 3: Boys localization -------------------------------------------
   ! The dipole matrices about the origin. Boys is origin-independent in its
   ! maximizer -- shifting the origin adds only terms invariant under an
   ! occupied-occupied rotation -- so no centre of mass is needed here.
   call tick(t0)
   call trc_multipoles(bas, pl, [0.0_dp, 0.0_dp, 0.0_dp], mm)
   !$acc update self(mm)
   allocate (loc(nao, nval), qrot(nval, nval), centroid(3, nval))
   call boys_localize()
   call tick(t1)
   print '(a,f9.3,a,i0,a)', '  stage 3    ', t1 - t0, ' s   localized ', nval, &
      ' valence orbitals'

   ! --- stage 4: distributed multipoles --------------------------------------
   call tick(t0)
   call distributed_multipoles()
   call tick(t1)
   print '(a,f9.3,a,i0,a)', '  stage 4    ', t1 - t0, ' s   ', npoint, &
      ' expansion points, charge through octupole'

   ! --- stage 5: static polarizability by CPHF -------------------------------
   allocate (alpha_lmo(3, 3, nval))
   call tick(t0)
   call cphf_static()
   call tick(t1)
   print '(a,f9.3,a,f10.4,a)', '  stage 5    ', t1 - t0, ' s   isotropic alpha ', &
      (alpha(1, 1) + alpha(2, 2) + alpha(3, 3))/3.0_dp, ' a.u.'

   ! --- stage 6: write it ----------------------------------------------------
   call tick(t0)
   call write_efp('demo.efp')
   call tick(t1)
   print '(a,f9.3,a)', '  stage 6    ', t1 - t0, ' s   wrote demo.efp'
   print '(a)', ''

   call summary()

   !$acc exit data delete(smat, tmat, vmat, mm)
   call eri%release(); call pl%release(); call bas%release()

contains

   subroutine tick(t)
      real(dp), intent(out) :: t
      integer(kind=8) :: c, rate
      call system_clock(c, rate)
      t = real(c, dp)/real(rate, dp)
   end subroutine tick

   !> The standard frozen core, which is what MAKEFP uses: its polarizable
   !> points and its exchange-repulsion orbitals are valence only.
   pure integer function frozen_core(n, z) result(c)
      integer, intent(in) :: n, z(n)
      integer :: i
      c = 0
      do i = 1, n
         if (z(i) > 2) c = c + 1
         if (z(i) > 10) c = c + 4
         if (z(i) > 18) c = c + 4
      end do
   end function frozen_core

   pure real(dp) function nuclear_repulsion(n, zz, rr) result(e)
      integer,  intent(in) :: n
      real(dp), intent(in) :: zz(n), rr(3, n)
      integer :: i, j
      e = 0.0_dp
      do i = 2, n
         do j = 1, i - 1
            e = e + zz(i)*zz(j)/norm2(rr(:, i) - rr(:, j))
         end do
      end do
   end function nuclear_repulsion

   ! =======================================================================
   ! stage 2
   ! =======================================================================
   subroutine run_scf()
      real(dp) :: eold, de, ddm
      real(dp), allocatable :: dold(:, :)
      integer :: it
      allocate (dold(nao, nao))
      call sym_orthog(nao, smat, x)
      call diag_fock(nao, hcore, x, cmo, eps)
      call make_density(nao, nocc, cmo, dmat)
      eold = 0.0_dp
      do it = 1, 200
         call eri%fock(bas, dmat, gmat)
         fock = hcore + gmat
         escf = enuc + sum(dmat*(hcore + fock))*0.5_dp
         dold = dmat
         call diag_fock(nao, fock, x, cmo, eps)
         call make_density(nao, nocc, cmo, dmat)
         de = abs(escf - eold); ddm = maxval(abs(dmat - dold))
         eold = escf
         if (de < 1.0e-10_dp .and. ddm < 1.0e-8_dp) exit
      end do
      if (it > 200) then
         print '(a)', '  the SCF did not converge'
         stop 1
      end if
      deallocate (dold)
   end subroutine run_scf

   subroutine sym_orthog(n, s, xo)
      integer, intent(in) :: n
      real(dp), intent(in) :: s(n, n)
      real(dp), intent(out) :: xo(n, n)
      real(dp), allocatable :: u(:, :), w(:), work(:)
      integer :: info, lwork, ii, jj
      allocate (u(n, n), w(n), work(64*n))
      u = s; lwork = 64*n
      call dsyev('V', 'U', n, u, n, w, work, lwork, info)
      if (info /= 0) then
         print '(a,i0)', '  dsyev failed in sym_orthog, info = ', info; stop 1
      end if
      do jj = 1, n
         do ii = 1, n
            xo(ii, jj) = u(ii, jj)/sqrt(w(jj))
         end do
      end do
      xo = matmul(xo, transpose(u))
      deallocate (u, w, work)
   end subroutine sym_orthog

   subroutine diag_fock(n, f, xo, c, e)
      integer, intent(in) :: n
      real(dp), intent(in) :: f(n, n), xo(n, n)
      real(dp), intent(out) :: c(n, n), e(n)
      real(dp), allocatable :: fp(:, :), work(:)
      integer :: info, lwork
      allocate (fp(n, n), work(64*n)); lwork = 64*n
      fp = matmul(transpose(xo), matmul(f, xo))
      call dsyev('V', 'U', n, fp, n, e, work, lwork, info)
      if (info /= 0) then
         print '(a,i0)', '  dsyev failed in diag_fock, info = ', info; stop 1
      end if
      c = matmul(xo, fp)
      deallocate (fp, work)
   end subroutine diag_fock

   subroutine make_density(n, no, c, d)
      integer, intent(in) :: n, no
      real(dp), intent(in) :: c(n, n)
      real(dp), intent(out) :: d(n, n)
      d = 2.0_dp*matmul(c(:, 1:no), transpose(c(:, 1:no)))
   end subroutine make_density

   ! =======================================================================
   ! stage 3 -- Boys localization by Jacobi sweeps
   ! =======================================================================
   subroutine boys_localize()
      !! Maximize sum_i |<i|r|i>|^2 over the valence occupied block.
      !!
      !! Two-by-two rotations, the classical Edmiston-style sweep. For a pair
      !! (p,q) with d = <p|r|p> - <q|r|q> and b = <p|r|q>,
      !!
      !!     A = b.b - d.d/4        B = d.b
      !!
      !! and the rotation that maximizes the functional is one quarter of
      !! atan2(B, -A). Sweeps stop when the largest rotation in a whole sweep
      !! stops moving anything.
      real(dp), allocatable :: rmo(:, :, :), tmp(:, :)
      real(dp) :: d(3), b(3), aa, bb, gam, cg, sg, big, vp, vq
      integer :: it, p, q, c, m
      allocate (rmo(nval, nval, 3), tmp(nao, nval))

      loc = cmo(:, ncore + 1:nocc)
      qrot = 0.0_dp
      do p = 1, nval
         qrot(p, p) = 1.0_dp
      end do

      do c = 1, 3
         tmp = matmul(mm(:, :, c), loc)
         rmo(:, :, c) = matmul(transpose(loc), tmp)
      end do

      do it = 1, 200
         big = 0.0_dp
         do p = 1, nval - 1
            do q = p + 1, nval
               do c = 1, 3
                  d(c) = rmo(p, p, c) - rmo(q, q, c)
                  b(c) = rmo(p, q, c)
               end do
               aa = dot_product(b, b) - 0.25_dp*dot_product(d, d)
               bb = dot_product(d, b)
               if (abs(aa) < 1.0e-14_dp .and. abs(bb) < 1.0e-14_dp) cycle
               gam = 0.25_dp*atan2(bb, -aa)
               if (abs(gam) < 1.0e-10_dp) cycle
               big = max(big, abs(gam))
               cg = cos(gam); sg = sin(gam)
               ! rotate the orbitals and every dipole matrix with them
               do m = 1, nao
                  vp = loc(m, p); vq = loc(m, q)
                  loc(m, p) =  cg*vp + sg*vq
                  loc(m, q) = -sg*vp + cg*vq
               end do
               do m = 1, nval
                  vp = qrot(m, p); vq = qrot(m, q)
                  qrot(m, p) =  cg*vp + sg*vq
                  qrot(m, q) = -sg*vp + cg*vq
               end do
               do c = 1, 3
                  do m = 1, nval
                     vp = rmo(m, p, c); vq = rmo(m, q, c)
                     rmo(m, p, c) =  cg*vp + sg*vq
                     rmo(m, q, c) = -sg*vp + cg*vq
                  end do
                  do m = 1, nval
                     vp = rmo(p, m, c); vq = rmo(q, m, c)
                     rmo(p, m, c) =  cg*vp + sg*vq
                     rmo(q, m, c) = -sg*vp + cg*vq
                  end do
               end do
            end do
         end do
         if (big < 1.0e-9_dp) exit
      end do

      do p = 1, nval
         do c = 1, 3
            centroid(c, p) = rmo(p, p, c)
         end do
      end do
      deallocate (rmo, tmp)
   end subroutine boys_localize

   ! =======================================================================
   ! stage 4 -- distributed multipoles
   ! =======================================================================
   subroutine distributed_multipoles()
      !! Charge through octupole at every expansion point.
      !!
      !! The partition is done at the PRIMITIVE level, which is why an
      !! uncontracted copy of the basis is built here: a contracted function
      !! spans a range of exponents and so a range of effective centres, and
      !! assigning the whole of it to one point throws that away. Each
      !! primitive shell pair goes to the expansion point nearest its
      !! Gaussian-product centre, ties split evenly.
      type(trc_basis_t)    :: unc
      type(trc_pairlist_t) :: upl
      real(dp), allocatable :: tmat_u(:, :), dprim(:, :), work(:, :)
      real(dp), allocatable :: us(:, :), ut(:, :), uv(:, :), umm(:, :, :)
      integer,  allocatable :: owner(:, :), nowner(:, :)
      real(dp), allocatable :: uexp(:)
      integer,  allocatable :: ush_l(:), ush_np(:), ush_ao(:)
      real(dp), allocatable :: ue(:, :), uc(:, :), ur(:, :)
      integer :: nu, nao_u, ip, iq, kk, kpt, ii, jj, mu, nu_i, cc
      integer :: du, dq, ou, oq
      real(dp) :: centre(3), best, d2, weight, pop
      real(dp), allocatable :: mono(:)

      call expansion_points()
      call build_uncontracted(unc, upl, tmat_u, uexp, ush_l, ush_np, &
                              ush_ao, ue, uc, ur, nu, nao_u)

      ! D' = T^T D T, the density in the primitive basis
      allocate (work(nao, nao_u), dprim(nao_u, nao_u))
      work = matmul(dmat, tmat_u)
      dprim = matmul(transpose(tmat_u), work)

      ! assignment of each primitive shell pair to its nearest point(s)
      allocate (owner(nu*nu, npoint), nowner(nu, nu))
      owner = 0; nowner = 0
      do iq = 1, nu
         do ip = 1, nu
            centre = (uexp(ip)*ur(:, ip) + uexp(iq)*ur(:, iq)) &
                     /(uexp(ip) + uexp(iq))
            best = huge(1.0_dp)
            do kpt = 1, npoint
               d2 = sum((centre - ppos(:, kpt))**2)
               if (d2 < best) best = d2
            end do
            kk = 0
            do kpt = 1, npoint
               d2 = sum((centre - ppos(:, kpt))**2)
               if (d2 - best <= TIE_TOL) then
                  kk = kk + 1
                  owner((iq - 1)*nu + ip, kk) = kpt
               end if
            end do
            nowner(ip, iq) = kk
         end do
      end do

      allocate (us(nao_u, nao_u), ut(nao_u, nao_u), uv(nao_u, nao_u))
      allocate (umm(nao_u, nao_u, TRC_NMULT))
      !$acc enter data create(us, ut, uv, umm)
      call trc_1e(unc, upl, us, ut, uv)
      !$acc update self(us)

      allocate (mono(npoint))
      allocate (qel(npoint), pdip(3, npoint), pquad(6, npoint), poct(10, npoint))
      mono = 0.0_dp; pdip = 0.0_dp; pquad = 0.0_dp; poct = 0.0_dp

      ! The monopole needs no origin, so it is done once outside the point loop.
      do iq = 1, nu
         dq = ncart(ush_l(iq)); oq = ush_ao(iq)
         do ip = 1, nu
            du = ncart(ush_l(ip)); ou = ush_ao(ip)
            if (nowner(ip, iq) == 0) cycle
            weight = 1.0_dp/real(nowner(ip, iq), dp)
            pop = 0.0_dp
            do jj = 1, dq
               do ii = 1, du
                  pop = pop + dprim(ou + ii, oq + jj)*us(ou + ii, oq + jj)
               end do
            end do
            do kk = 1, nowner(ip, iq)
               kpt = owner((iq - 1)*nu + ip, kk)
               mono(kpt) = mono(kpt) + weight*pop
            end do
         end do
      end do
      qel = -mono          ! electrons are negative; the nuclear column is separate

      ! Dipole and above are expanded about each point in turn. The integrals
      ! are recomputed per point rather than translated -- npoint is a handful,
      ! and a translation is one more binomial to get wrong.
      do kpt = 1, npoint
         call trc_multipoles(unc, upl, ppos(:, kpt), umm)
         !$acc update self(umm)
         do iq = 1, nu
            dq = ncart(ush_l(iq)); oq = ush_ao(iq)
            do ip = 1, nu
               du = ncart(ush_l(ip)); ou = ush_ao(ip)
               if (nowner(ip, iq) == 0) cycle
               weight = 0.0_dp
               do kk = 1, nowner(ip, iq)
                  if (owner((iq - 1)*nu + ip, kk) == kpt) &
                     weight = 1.0_dp/real(nowner(ip, iq), dp)
               end do
               if (weight == 0.0_dp) cycle
               do jj = 1, dq
                  do ii = 1, du
                     mu = ou + ii; nu_i = oq + jj
                     do cc = 1, 3
                        pdip(cc, kpt) = pdip(cc, kpt) &
                           - weight*dprim(mu, nu_i)*umm(mu, nu_i, cc)
                     end do
                     do cc = 1, 6
                        pquad(cc, kpt) = pquad(cc, kpt) &
                           - weight*dprim(mu, nu_i)*umm(mu, nu_i, 3 + QUAD_PACK(cc))
                     end do
                     do cc = 1, 10
                        poct(cc, kpt) = poct(cc, kpt) &
                           - weight*dprim(mu, nu_i)*umm(mu, nu_i, 12 + OCT_PACK(cc))
                     end do
                  end do
               end do
            end do
         end do
      end do

      !$acc exit data delete(us, ut, uv, umm)
      call upl%release(); call unc%release()
      deallocate (work, dprim, owner, nowner, mono, us, ut, uv, umm)
      deallocate (tmat_u, uexp, ush_l, ush_np, ush_ao, ue, uc, ur)
   end subroutine distributed_multipoles

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   subroutine expansion_points()
      !! Every atom, then every bond midpoint.
      !!
      !! Labels follow GAMESS: `A<nn><symbol>` for atoms in input order and
      !! `BO<hi><lo>` for a midpoint. A midpoint carries no charge and no mass,
      !! which is how a `.efp` reader tells the two kinds of point apart.
      integer :: nb, ii, jj, kk
      integer, allocatable :: bi(:), bj(:)
      real(dp) :: r, lim
      allocate (bi(nat*(nat - 1)/2), bj(nat*(nat - 1)/2))
      nb = 0
      do ii = 2, nat
         do jj = 1, ii - 1
            r = norm2(at_r(:, ii) - at_r(:, jj))/BOHR_PER_ANG
            lim = COVALENT(at_z(ii)) + COVALENT(at_z(jj))
            if (r <= lim) then
               nb = nb + 1; bi(nb) = ii; bj(nb) = jj
            end if
         end do
      end do
      npoint = nat + nb
      allocate (ppos(3, npoint), plabel(npoint), qnuc(npoint))
      do ii = 1, nat
         ppos(:, ii) = at_r(:, ii)
         qnuc(ii) = real(at_z(ii), dp)
         write (plabel(ii), '(a,i2.2,a)') 'A', ii, trim(symbol(at_z(ii)))
      end do
      do kk = 1, nb
         ! the plain arithmetic mean, not a weighted or shifted midpoint
         ppos(:, nat + kk) = 0.5_dp*(at_r(:, bi(kk)) + at_r(:, bj(kk)))
         qnuc(nat + kk) = 0.0_dp
         write (plabel(nat + kk), '(a,i0,i0)') 'BO', bi(kk), bj(kk)
      end do
      deallocate (bi, bj)
   end subroutine expansion_points

   pure function symbol(z) result(s)
      integer, intent(in) :: z
      character(len=2) :: s
      character(len=2), parameter :: tab(MAXEL) = [ &
         'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
         'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
         'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', &
         'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']
      s = tab(z)
   end function symbol

   subroutine build_uncontracted(u, up, tr, sexp, ul, unp, uao, uee, ucc, &
                                 urr, nu, naou)
      !! One shell per primitive, and the transform back to the contracted set.
      !!
      !! terco's `build` multiplies every coefficient by `common_fac_sp`, and
      !! `build_631g` has already folded in the primitive norm. Giving the
      !! uncontracted shells the SAME stored coefficients the contracted ones
      !! carry makes each uncontracted function equal to its term in the
      !! contraction, so the transform is the identity within a shell: column
      !! mu of T picks out the primitives that build it, with coefficient one.
      !! Nothing here has to know what normalisation convention was used
      !! upstream, which is the point -- the alternative divides by a norm and
      !! is wrong the moment the caller's convention changes.
      type(trc_basis_t),    intent(out) :: u
      type(trc_pairlist_t), intent(out) :: up
      real(dp), allocatable, intent(out) :: tr(:, :), sexp(:), uee(:, :)
      real(dp), allocatable, intent(out) :: ucc(:, :), urr(:, :)
      integer,  allocatable, intent(out) :: ul(:), unp(:), uao(:)
      integer, intent(out) :: nu, naou
      integer :: is, ip, s, off, a, n

      nu = 0
      do is = 1, nsh
         nu = nu + sh_np(is)
      end do
      allocate (ul(nu), unp(nu), uao(nu), sexp(nu))
      allocate (uee(1, nu), ucc(1, nu), urr(3, nu))
      s = 0; naou = 0
      do is = 1, nsh
         do ip = 1, sh_np(is)
            s = s + 1
            ul(s) = sh_l(is); unp(s) = 1
            uee(1, s) = sh_e(ip, is)
            ucc(1, s) = sh_c(ip, is)
            urr(:, s) = sh_r(:, is)
            sexp(s) = sh_e(ip, is)
            uao(s) = naou
            naou = naou + ncart(sh_l(is))
         end do
      end do

      call u%build(nu, ul, unp, uee, ucc, urr, nat, zat, at_r, 1)
      call up%build(u, thresh)
      call u%to_device(); call up%to_device()

      ! T(primitive AO, contracted AO)
      allocate (tr(nao, naou))
      tr = 0.0_dp
      s = 0; off = 0
      do is = 1, nsh
         n = ncart(sh_l(is))
         do ip = 1, sh_np(is)
            s = s + 1
            do a = 1, n
               tr(off + a, uao(s) + a) = 1.0_dp
            end do
         end do
         off = off + n
      end do
   end subroutine build_uncontracted

   ! =======================================================================
   ! stage 5 -- static polarizability by CPHF, batched over perturbations
   ! =======================================================================
   subroutine cphf_static()
      !! alpha_ij = -Tr(D^(j) r_i), with D^(j) the density response to a field
      !! along j.
      !!
      !! The perturbation is the dipole matrix terco already produced for the
      !! localization, so nothing new is integrated here. The iteration is the
      !! textbook one:
      !!
      !!     F1 = h1 + G(D1)                  <- the Fock build, the whole cost
      !!     U_ai = F1_ai / (eps_i - eps_a)
      !!     D1 = 2 (Cv U Co^T + Co U^T Cv^T)
      !!
      !! ALL THREE PERTURBATIONS GO THROUGH `fock_many` IN ONE PASS. That is
      !! the demo: the real potential has about a hundred right-hand sides, and
      !! this is the routine they would go through. Three is well under the
      !! N = 4 optimum Phase 2 measured, so the batching is not being flattered
      !! here -- it is being used exactly as makeEFP would use it.
      !!
      !! The iteration itself is a damped fixed point, which is the one place
      !! this demo is deliberately cruder than mqc: it takes about sixty passes
      !! where a DIIS or conjugate-gradient solver takes ten. That is the
      !! CALLER's half of the contract and it stays in mqc, so making it good
      !! here would be work that the backend throws away. It does mean the
      !! stage 5 time below is roughly six times what a real solver would
      !! spend -- read it as a Fock-build count, not as a solver benchmark.
      real(dp), allocatable :: d1(:, :, :), g1(:, :, :), u(:, :, :)
      real(dp), allocatable :: f1(:, :), fmo(:, :), w(:, :), tmp(:, :)
      real(dp), allocatable :: co(:, :), cv(:, :)
      real(dp), allocatable :: wi(:, :, :), uv_(:, :, :), wv(:, :, :)
      real(dp) :: conv, damp
      integer :: it, c, a, ii, cj

      allocate (d1(3, nao, nao), g1(3, nao, nao), u(nvirt, nocc, 3))
      allocate (f1(nao, nao), fmo(nao, nao), w(nao, nao), tmp(nao, nao))
      allocate (co(nao, nocc), cv(nao, nvirt))
      co = cmo(:, 1:nocc)
      cv = cmo(:, nocc + 1:nao)

      d1 = 0.0_dp
      u = 0.0_dp
      damp = 0.5_dp
      do it = 1, 100
         call eri%fock_many(bas, 3, d1, g1)
         conv = 0.0_dp
         do c = 1, 3
            f1 = mm(:, :, c) + g1(c, :, :)
            tmp(:, 1:nocc) = matmul(f1, co)
            fmo(1:nvirt, 1:nocc) = matmul(transpose(cv), tmp(:, 1:nocc))
            do ii = 1, nocc
               do a = 1, nvirt
                  w(a, ii) = fmo(a, ii)/(eps(ii) - eps(nocc + a))
               end do
            end do
            do ii = 1, nocc
               do a = 1, nvirt
                  conv = max(conv, abs(w(a, ii) - u(a, ii, c)))
                  u(a, ii, c) = (1.0_dp - damp)*u(a, ii, c) + damp*w(a, ii)
               end do
            end do
            tmp = matmul(cv, matmul(u(:, :, c), transpose(co)))
            d1(c, :, :) = 2.0_dp*(tmp + transpose(tmp))
         end do
         if (conv < 1.0e-8_dp) exit
      end do
      if (it > 100) print '(a)', '  WARNING: CPHF did not converge'
      print '(a,i0,a)', '               CPHF converged in ', min(it, 100), &
         ' iterations (3 right-hand sides per Fock build)'

      do cj = 1, 3
         do c = 1, 3
            alpha(c, cj) = -sum(d1(cj, :, :)*mm(:, :, c))
         end do
      end do

      ! Per-LMO decomposition, which is what makes these polarizable POINTS.
      ! alpha = -4 Tr(U^T W) over the occupied block, and that trace is
      ! invariant under the occupied-occupied rotation the localization
      ! applied -- so rotating both factors gives a partition that sums back
      ! to the same total.
      allocate (wi(nvirt, nocc, 3), uv_(nvirt, nval, 3), wv(nvirt, nval, 3))
      do c = 1, 3
         tmp(:, 1:nocc) = matmul(mm(:, :, c), co)
         wi(:, :, c) = matmul(transpose(cv), tmp(:, 1:nocc))
         uv_(:, :, c) = matmul(u(:, ncore + 1:nocc, c), qrot)
         wv(:, :, c) = matmul(wi(:, ncore + 1:nocc, c), qrot)
      end do
      do ii = 1, nval
         do cj = 1, 3
            do c = 1, 3
               alpha_lmo(c, cj, ii) = -4.0_dp*sum(uv_(:, ii, cj)*wv(:, ii, c))
            end do
         end do
      end do
      deallocate (wi, uv_, wv)

      deallocate (d1, g1, u, f1, fmo, w, tmp, co, cv)
   end subroutine cphf_static

   ! =======================================================================
   ! stage 6
   ! =======================================================================
   subroutine write_efp(fn)
      character(len=*), intent(in) :: fn
      integer :: uu, p, c
      open (newunit=uu, file=fn, status='replace', action='write')
      write (uu, '(a)') ' $DEMO'
      write (uu, '(a)') '! Written by terco''s demo_makefp. The sections below are the'
      write (uu, '(a)') '! ones this demo computes; a real potential also carries the'
      write (uu, '(a)') '! dynamic polarizability blocks, the exchange-repulsion basis'
      write (uu, '(a)') '! and the charge-penetration screening.'
      write (uu, '(a,f20.10)') '! RHF energy ', escf

      write (uu, '(a)') ' COORDINATES (BOHR)'
      do p = 1, npoint
         write (uu, '(a8,3f16.10,2f10.4)') plabel(p), ppos(:, p), &
            merge(mass(at_z(min(p, nat))), 0.0_dp, p <= nat), qnuc(p)
      end do
      write (uu, '(a)') ' STOP'

      write (uu, '(a)') ' MONOPOLES'
      do p = 1, npoint
         write (uu, '(a8,2f20.10)') plabel(p), qel(p), qnuc(p)
      end do
      write (uu, '(a)') ' STOP'

      write (uu, '(a)') ' DIPOLES'
      do p = 1, npoint
         write (uu, '(a8,3f20.10)') plabel(p), pdip(:, p)
      end do
      write (uu, '(a)') ' STOP'

      write (uu, '(a)') ' QUADRUPOLES'
      do p = 1, npoint
         write (uu, '(a8,6f18.10)') plabel(p), pquad(:, p)
      end do
      write (uu, '(a)') ' STOP'

      write (uu, '(a)') ' OCTUPOLES'
      do p = 1, npoint
         write (uu, '(a8,10f18.10)') plabel(p), poct(:, p)
      end do
      write (uu, '(a)') ' STOP'

      write (uu, '(a)') ' POLARIZABLE POINTS'
      do p = 1, nval
         write (uu, '(a,i3,3f16.10)') ' CT', p, centroid(:, p)
         write (uu, '(9f14.8)') ((alpha_lmo(c, k, p), c=1, 3), k=1, 3)
      end do
      write (uu, '(a)') ' STOP'
      write (uu, '(a)') ' $END'
      close (uu)
   end subroutine write_efp

   pure real(dp) function mass(z)
      integer, intent(in) :: z
      real(dp), parameter :: tab(10) = [1.00794_dp, 4.0026_dp, 6.941_dp, &
         9.0122_dp, 10.811_dp, 12.0107_dp, 14.0067_dp, 15.9994_dp, &
         18.9984_dp, 20.1797_dp]
      mass = 0.0_dp
      if (z >= 1 .and. z <= 10) mass = tab(z)
   end function mass

   subroutine summary()
      integer :: p, c
      real(dp) :: qtot
      print '(a)', '  distributed multipoles'
      print '(a)', '    point        q(elec)     q(nuc)      |dipole|'
      qtot = 0.0_dp
      do p = 1, npoint
         qtot = qtot + qel(p) + qnuc(p)
         print '(a,a8,3f13.6)', '    ', plabel(p), qel(p), qnuc(p), &
            norm2(pdip(:, p))
      end do
      print '(a,f13.6)', '    net charge (should be the molecular charge) ', qtot
      print '(a)', ''
      print '(a)', '  localized orbital centroids and their polarizabilities'
      do p = 1, nval
         print '(a,i3,a,3f10.5,a,f10.5)', '    LMO', p, '  at', centroid(:, p), &
            '   alpha(iso)', (alpha_lmo(1,1,p) + alpha_lmo(2,2,p) &
                              + alpha_lmo(3,3,p))/3.0_dp
      end do
      print '(a)', ''
      print '(a)', '  static polarizability tensor (a.u.)'
      do c = 1, 3
         print '(a,3f14.6)', '    ', alpha(c, 1), alpha(c, 2), alpha(c, 3)
      end do
   end subroutine summary

end program demo_makefp
