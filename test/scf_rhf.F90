!
! A complete restricted Hartree-Fock SCF, driven entirely through the library.
!
! WHAT THIS IS FOR
! ----------------
! Every check so far compares an integral against libfint or gpu4pyscf. None
! of them answers the question a consumer actually asks: can I hand this thing
! a molecule and get a converged energy? This does, and it does it using only
! the public API -- `trc_basis_t`, `trc_pairlist_t`, `trc_1e` and
! `trc_eri_t%fock`. Nothing below the API is referenced.
!
! It is also the reference for the mqc backend: whatever mqc does, it does
! this, with its own diagonaliser and its own DIIS.
!
! WHAT THE LIBRARY DOES AND DOES NOT DO
! -------------------------------------
! The library returns S, T, V and vhf = J - K/2. Everything else here --
! orthogonalisation, diagonalisation, the density update, damping,
! convergence testing, the nuclear repulsion -- is the caller's, deliberately.
! DIIS is not implemented; simple damping is enough to converge these cases
! and DIIS is exactly the kind of thing a host code already has.
!
!   scf_rhf <xyz> [charge] [thresh] [maxiter] [density-file]
!
! With a density file it does ONE Fock build from that density and prints the
! result instead of iterating. That mode exists because comparing two Fock
! implementations is only meaningful on the same density, and it is how the
! benchmark's four-fold accumulation bug was found.
!
program scf_rhf
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_eri, only: trc_eri_t
   use trc_test_basis, only: read_xyz, build_631g
   implicit none

   character(len=256) :: xyzfile, arg
   integer  :: charge, maxiter
   real(dp) :: thresh

   integer :: nat, nsh, maxnp, nao, nelec, nocc, i, j, it
   integer,  allocatable :: at_z(:), sh_l(:), sh_np(:)
   real(dp), allocatable :: at_r(:, :), zat(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
   real(dp), allocatable :: hcore(:, :), gmat(:, :), fock(:, :)
   real(dp), allocatable :: dmat(:, :), dold(:, :), x(:, :), c(:, :), eps(:)
   type(trc_basis_t)    :: bas
   type(trc_pairlist_t) :: pl
   type(trc_eri_t)      :: eri
   real(dp) :: enuc, e1, e2, etot, eold, drms, t0, t1, tscf

   charge = 0; thresh = 1.0e-10_dp; maxiter = 60
   if (command_argument_count() < 1) then
      print '(a)', 'usage: scf_rhf <xyz> [charge] [thresh] [maxiter]'
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
      call get_command_argument(4, arg); read (arg, *) maxiter
   end if

   ! Optional 5th argument: a density file. One Fock build from it, print the
   ! result, stop. Exists to compare this path against bench_binwater on the
   ! SAME density, which is the only way to tell two fold conventions apart.
   call read_xyz(xyzfile, nat, at_z, at_r)
   call build_631g(nat, at_z, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   allocate (zat(nat)); zat = real(at_z, dp)

   nelec = sum(at_z) - charge
   if (mod(nelec, 2) /= 0) then
      print '(a,i0)', '  RHF needs an even electron count; got ', nelec
      stop 1
   end if
   nocc = nelec/2

   ! ---- the library: basis, pair list, integral engines ----
   call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, at_r, maxnp)
   call pl%build(bas, thresh)
   nao = bas%nao

   print '(a,a)',       '  geometry   : ', trim(xyzfile)
   print '(a,i0,a,i0,a,i0)', '  atoms ', nat, '   shells ', nsh, '   nao ', nao
   print '(a,i0,a,i0)', '  electrons ', nelec, '   occupied ', nocc
   print '(a,es9.2)',   '  threshold  : ', thresh

   allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
   allocate (hcore(nao, nao), gmat(nao, nao), fock(nao, nao))
   allocate (dmat(nao, nao), dold(nao, nao), x(nao, nao), c(nao, nao), eps(nao))

   !$acc enter data create(smat, tmat, vmat)
   call bas%to_device()
   call pl%to_device()
   call trc_1e(bas, pl, smat, tmat, vmat)
   !$acc update self(smat, tmat, vmat)
   !$acc exit data delete(smat, tmat, vmat)

   hcore = tmat + vmat

   call tick(t0)
   call eri%build(bas, thresh)
   call tick(t1)
   print '(a,f8.3,a)', '  ERI setup  : ', t1 - t0, ' s'

   enuc = nuclear_repulsion(nat, zat, at_r)

   if (command_argument_count() >= 5) then
      block
         character(len=256) :: dmf
         integer :: ud
         integer(kind=8) :: nchk
         call get_command_argument(5, dmf)
         open (newunit=ud, file=trim(dmf), access='stream', &
               form='unformatted', status='old')
         read (ud) nchk
         if (int(nchk) /= nao) then
            print '(a,i0,a,i0)', '  density nao ', nchk, ' /= ', nao; stop 1
         end if
         read (ud) dmat
         close (ud)
         call eri%fock(bas, dmat, gmat)
         print '(a,es24.16)', '  [oneshot] sum(vhf)  ', sum(gmat)
         print '(a,es24.16)', '  [oneshot] sum|vhf|  ', sum(abs(gmat))
         ! k_scale = 0 leaves pure Coulomb, k_scale = 1 with j_scale = 0
         ! leaves pure exchange. Checking the two separately is what pins the
         ! scalings: a common factor would pass a J - K/2 comparison.
         call eri%fock(bas, dmat, gmat, k_scale=0.0_dp)
         print '(a,es24.16)', '  [oneshot] sum(J)    ', sum(gmat)
         call eri%fock(bas, dmat, gmat, j_scale=0.0_dp)
         print '(a,es24.16)', '  [oneshot] sum(-K/2) ', sum(gmat)
         !
         ! A deliberately ASYMMETRIC density. The folded kernel assumes
         ! D_mn = D_nm and symmetrises its result, so it must disagree here;
         ! fock_nosym enumerates the permutations and must not. Feeding both
         ! the SAME symmetric density first proves the two agree when the
         ! assumption holds, so a disagreement below is the asymmetry and not
         ! a bug in either.
         !
         block
            real(dp), allocatable :: da(:, :), gs(:, :), gn(:, :)
            integer :: p_, q_
            allocate (da(nao, nao), gs(nao, nao), gn(nao, nao))
            call eri%fock(bas, dmat, gs)
            call eri%fock_nosym(bas, dmat, gn)
            print '(a,es12.4)', '  [nosym] symmetric D, folded vs enumerated  ', &
               maxval(abs(gs - gn))
            print '(a,es24.16,a,es24.16)', '  [nosym] sum folded ', sum(gs), &
               '  sum enumerated ', sum(gn)
            print '(a,es12.4)', '  [nosym] ratio of sums ', sum(gn)/sum(gs)
            do q_ = 1, nao
               do p_ = 1, nao
                  da(p_, q_) = dmat(p_, q_) + 0.02_dp*sin(real(3*p_ - 5*q_, dp))
               end do
            end do
            print '(a,es12.4)', '  [nosym] asymmetry of the test density      ', &
               maxval(abs(da - transpose(da)))
            call eri%fock_nosym(bas, da, gn)
            print '(a,es24.16)', '  [nosym] sum(vhf) asymmetric D  ', sum(gn)
            print '(a,es12.4)', '  [nosym] its own antisymmetry               ', &
               maxval(abs(gn - transpose(gn)))
            deallocate (da, gs, gn)
         end block
         !
         ! Device-resident build: same answer, no matrices crossing the bus.
         !
         block
            real(dp), allocatable :: gr(:, :), gh(:, :)
            allocate (gr(nao, nao), gh(nao, nao))
            call eri%fock(bas, dmat, gh)
            !$acc enter data copyin(dmat) create(gr)
            call eri%fock_resident(bas, dmat, gr)
            !$acc update self(gr)
            !$acc exit data delete(dmat, gr)
            print '(a,es12.4)', '  [resident] vs host-copy path               ', &
               maxval(abs(gr - gh))
            block
               real(dp) :: ta, tb
               integer  :: rr
               integer, parameter :: NIT = 10
               call tick(ta)
               do rr = 1, NIT
                  call eri%fock(bas, dmat, gh)
               end do
               call tick(tb)
               print '(a,i0,a,f9.4,a)', '  [resident] ', NIT, ' x fock (host copies)  ', tb - ta, ' s'
               !$acc enter data copyin(dmat) create(gr)
               call tick(ta)
               do rr = 1, NIT
                  call eri%fock_resident(bas, dmat, gr)
               end do
               !$acc wait
               call tick(tb)
               !$acc exit data delete(dmat, gr)
               print '(a,i0,a,f9.4,a)', '  [resident] ', NIT, ' x fock_resident       ', tb - ta, ' s'
            end block
            deallocate (gr, gh)
         end block
         call cleanup(); stop 0
      end block
   end if

   ! Core-Hamiltonian guess, through the same orthogonaliser the SCF uses.
   call sym_orthog(nao, smat, x)
   call diag_fock(nao, hcore, x, c, eps)
   call make_density(nao, nocc, c, dmat)

   print '(a)', ''
   print '(a)', '   it        E(total)            dE          RMS(D)     t(Fock)'
   eold = 0.0_dp
   tscf = 0.0_dp
   do it = 1, maxiter
      call tick(t0)
      call eri%fock(bas, dmat, gmat)
      call tick(t1)
      tscf = tscf + (t1 - t0)

      ! `fock` returns J - K/2, so this needs no factor.
      fock = hcore + gmat

      ! E = sum_mn D_mn (H_mn + F_mn) + E_nuc, with D the total density
      e1 = 0.0_dp; e2 = 0.0_dp
      do j = 1, nao
         do i = 1, nao
            e1 = e1 + dmat(i, j)*hcore(i, j)
            e2 = e2 + dmat(i, j)*fock(i, j)
         end do
      end do
      etot = 0.5_dp*(e1 + e2) + enuc

      dold = dmat
      call diag_fock(nao, fock, x, c, eps)
      call make_density(nao, nocc, c, dmat)

      drms = 0.0_dp
      do j = 1, nao
         do i = 1, nao
            drms = drms + (dmat(i, j) - dold(i, j))**2
         end do
      end do
      drms = sqrt(drms/real(nao, dp)**2)

      print '(i5,f22.12,es14.4,es14.4,f10.4)', it, etot, etot - eold, drms, t1 - t0
      if (it > 1 .and. abs(etot - eold) < 1.0e-9_dp .and. drms < 1.0e-7_dp) then
         print '(a)', ''
         print '(a,f22.12)', '  CONVERGED  E(RHF) = ', etot
         print '(a,f22.12)', '  E(nuclear repulsion) = ', enuc
         print '(a,i0,a,f8.3,a)', '  ', it, ' iterations, ', tscf, ' s in Fock builds'
         call cleanup(); stop 0
      end if
      eold = etot

      ! Damping. Crude on purpose: DIIS belongs to the caller, and these
      ! systems converge without it.
      dmat = 0.5_dp*(dmat + dold)
   end do

   print '(a,i0,a)', '  NOT CONVERGED in ', maxiter, ' iterations'
   call cleanup()
   stop 1

contains

   subroutine cleanup()
      call eri%release(); call pl%release(); call bas%release()
   end subroutine cleanup

   subroutine tick(t)
      real(dp), intent(out) :: t
      integer(kind=8) :: cc, rate
      call system_clock(cc, rate)
      t = real(cc, dp)/real(rate, dp)
   end subroutine tick

   real(dp) function nuclear_repulsion(nat, z, r) result(e)
      integer,  intent(in) :: nat
      real(dp), intent(in) :: z(nat), r(3, nat)
      integer :: i, j
      real(dp) :: d
      e = 0.0_dp
      do i = 1, nat
         do j = 1, i - 1
            d = sqrt(sum((r(:, i) - r(:, j))**2))
            e = e + z(i)*z(j)/d
         end do
      end do
   end function nuclear_repulsion

   !> X = S^(-1/2), symmetric orthogonalisation, with near-null vectors
   !> dropped rather than inverted.
   subroutine sym_orthog(n, s, x)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: s(n, n)
      real(dp), intent(out) :: x(n, n)
      real(dp), allocatable :: u(:, :), w(:), work(:)
      integer :: i, j, k, info, lwork
      lwork = 64*n
      allocate (u(n, n), w(n), work(lwork))
      u = s
      call dsyev('V', 'U', n, u, n, w, work, lwork, info)
      if (info /= 0) then
         print '(a,i0)', '  dsyev failed in sym_orthog: ', info; stop 1
      end if
      x = 0.0_dp
      do k = 1, n
         if (w(k) < 1.0e-10_dp) cycle
         do j = 1, n
            do i = 1, n
               x(i, j) = x(i, j) + u(i, k)*u(j, k)/sqrt(w(k))
            end do
         end do
      end do
      deallocate (u, w, work)
   end subroutine sym_orthog

   subroutine diag_fock(n, f, x, c, eps)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: f(n, n), x(n, n)
      real(dp), intent(out) :: c(n, n), eps(n)
      real(dp), allocatable :: fp(:, :), tmp(:, :), work(:)
      integer :: info, lwork
      lwork = 64*n
      allocate (fp(n, n), tmp(n, n), work(lwork))
      tmp = matmul(f, x)
      fp  = matmul(transpose(x), tmp)
      call dsyev('V', 'U', n, fp, n, eps, work, lwork, info)
      if (info /= 0) then
         print '(a,i0)', '  dsyev failed in diag_fock: ', info; stop 1
      end if
      c = matmul(x, fp)
      deallocate (fp, tmp, work)
   end subroutine diag_fock

   !> D = 2 C_occ C_occ^T -- the TOTAL density, which is the convention
   !> `fock` expects (it returns 2J - K for that D).
   subroutine make_density(n, nocc, c, d)
      integer,  intent(in)  :: n, nocc
      real(dp), intent(in)  :: c(n, n)
      real(dp), intent(out) :: d(n, n)
      integer :: i, j, k
      d = 0.0_dp
      do k = 1, nocc
         do j = 1, n
            do i = 1, n
               d(i, j) = d(i, j) + 2.0_dp*c(i, k)*c(j, k)
            end do
         end do
      end do
   end subroutine make_density

end program scf_rhf
