!
! A complete restricted Kohn-Sham SCF through the public API, and the
! end-to-end comparison against pyscf on the same grid.
!
! terco's contract stops at the matrices: given a density, a Fock matrix;
! given a density and a grid, E_xc and V_xc. The eigenproblem and the
! iteration are the caller's, which is why this lives in test/ beside
! scf_rhf.F90 and needs LAPACK when the library does not.
!
! The system is water in 6-31G UNCONTRACTED to its primitives: every
! primitive its own shell with coefficient one, 30 Cartesian functions.
! Uncontracted because pyscf renormalises a contracted shell to unit
! self-overlap and terco does not, and a total energy compared at 1e-8
! cannot carry the 1e-6 that leaves in the basis; a lone primitive is
! normalised identically by both. And water rather than the three-centre
! probe of the other checks, because that probe's toy basis has no core
! function on the beryllium and its LDA HOMO-LUMO gap is a millihartree:
! pyscf's own DIIS fails on it too. A fixed density does not care; an SCF
! does.
!
! Four functionals run in one go -- SVWN, PBE, BLYP and B3LYP -- so the
! LDA, GGA and hybrid paths through the whole SCF are all exercised. The
! grid, the shells exactly as built, the energies and the final densities
! go to `rks_probe.bin`; `rks_ref.py` rebuilds the same basis in the same
! order, runs pyscf's RKS on the same points, and compares.
!
! The SCF is four damped iterations from the core guess, then Pulay's
! DIIS on the commutator FDS - SDF in the orthogonal basis, ten vectors
! deep. A converged energy is what is compared, and the route to it is not
! the thing under test; the solver is only as elaborate as this cation
! needed -- plain damping took two hundred iterations on the hybrid, and
! DIIS from the core guess oscillated on the LDA.
!
program scf_rks
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_fock, only: trc_eri_t
   use trc_error, only: error_t
   use trc_dft_grid, only: dft_grid_t, build_dft_grid
   use trc_dft_prune, only: PRUNE_NONE
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_xc, only: trc_xc_rks
   use trc_test_basis, only: read_xyz, build_631g
   implicit none

   integer, parameter :: NFUNC = 4, NELEC = 10, NOCC = NELEC/2
   integer, parameter :: MAXITER = 100, NDIIS = 10, NDAMP = 4
   character(len=8), parameter :: names(NFUNC) = [character(len=8) :: "svwn", "pbe", "blyp", "b3lyp"]
   integer :: i, j, nao, ifn, it, unit, natm, nsh, maxnp, ish
   integer, allocatable :: sh_l(:), sh_np(:), zint(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), at_z(:), at_r(:, :)
   real(dp) :: enuc, e1, e2, exc, nelec_grid, etot, eold, drms, energies(NFUNC)
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :), hcore(:, :), gmat(:, :), fock(:, :)
   real(dp), allocatable :: vxc(:, :), dmat(:, :), dold(:, :), x(:, :), c(:, :), eps(:)
   real(dp), allocatable :: fstore(:, :, :), estore(:, :, :), errv(:, :)
   integer :: ndiis_used
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   type(trc_eri_t) :: eri
   type(dft_grid_t) :: grid
   type(trc_xc_grid_t) :: xg
   type(trc_xc_functional_t) :: func
   type(error_t) :: err
   logical :: converged, all_ok

   call read_xyz('water.xyz', natm, zint, at_r)
   allocate (at_z(natm))
   at_z = real(zint, dp)
   call build_631g(natm, zint, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp, uncontracted=.true.)
   call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, at_z, at_r, maxnp)
   nao = bas%nao
   call pl%build(bas, 1.0e-30_dp)
   allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao), hcore(nao, nao), gmat(nao, nao))
   allocate (fock(nao, nao), vxc(nao, nao), dmat(nao, nao), dold(nao, nao), x(nao, nao), c(nao, nao), eps(nao))
   allocate (fstore(nao, nao, NDIIS), estore(nao, nao, NDIIS), errv(nao, nao))
   call bas%to_device(); call pl%to_device()
   call trc_1e(bas, pl, smat, tmat, vmat)
   hcore = tmat + vmat
   call eri%build(bas, 1.0e-12_dp)
   enuc = nuclear_repulsion(NATM, at_z, at_r)

   call build_dft_grid(at_r, zint, grid, err, level=3, prune=PRUNE_NONE)
   if (err%has_error()) then
      print '(a)', "scf_rks: grid failed: "//err%get_message()
      stop 1
   end if
   call xg%build(grid%n_points, grid%coords, grid%weights, bas)
   call xg%to_device()
   print '(a,i0,a,i0,a,i0,a)', "scf_rks: water, 6-31G uncontracted, ", nsh, " shells, ", nao, " AOs, ", &
      grid%n_points, " grid points"

   open (newunit=unit, file='rks_probe.bin', access='stream', form='unformatted', status='replace')
   write (unit) int(grid%n_points, kind=8), int(nao, kind=8), int(NFUNC, kind=8), int(natm, kind=8), int(nsh, kind=8)
   write (unit) grid%coords
   write (unit) grid%weights
   write (unit) at_z, at_r
   ! Each shell: its centre as an atom index, its l, its one exponent.
   do ish = 1, nsh
      write (unit) int(centre_of(ish), kind=8), int(sh_l(ish), kind=8), sh_e(1, ish)
   end do

   call sym_orthog(nao, smat, x)
   all_ok = .true.
   do ifn = 1, NFUNC
      if (.not. xc_functional_by_name(names(ifn), func)) then
         print '(a)', "scf_rks: unknown functional "//trim(names(ifn))
         stop 1
      end if
      ! Core guess, DIIS iteration.
      call diag_fock(nao, hcore, x, c, eps)
      call make_density(nao, NOCC, c, dmat)
      eold = 0.0_dp
      converged = .false.
      ndiis_used = 0
      print '(a)', ""
      print '(a,a)', "  functional ", trim(names(ifn))
      print '(a)', "   it        E(total)            dE          RMS(D)         N(grid)"
      do it = 1, MAXITER
         call eri%fock(bas, dmat, gmat, k_scale=func%exx)
         call trc_xc_rks(bas, xg, func, dmat, vxc, exc, nelec_grid)
         fock = hcore + gmat + vxc
         e1 = 0.0_dp; e2 = 0.0_dp
         do j = 1, nao
            do i = 1, nao
               e1 = e1 + dmat(i, j)*hcore(i, j)
               e2 = e2 + dmat(i, j)*gmat(i, j)
            end do
         end do
         etot = e1 + 0.5_dp*e2 + exc + enuc
         dold = dmat
         ! DIIS: error in the orthogonal basis, then the extrapolated Fock.
         errv = matmul(transpose(x), matmul(matmul(fock, matmul(dmat, smat)) - matmul(smat, matmul(dmat, fock)), x))
         if (it > NDAMP) then
            call diis_push(nao, fock, errv, fstore, estore, ndiis_used)
            call diis_extrapolate(nao, fstore, estore, ndiis_used, fock)
         end if
         call diag_fock(nao, fock, x, c, eps)
         call make_density(nao, NOCC, c, dmat)
         if (it <= NDAMP) dmat = 0.5_dp*(dmat + dold)
         drms = sqrt(sum((dmat - dold)**2))/real(nao, dp)
         if (it <= 3 .or. mod(it, 10) == 0 .or. drms < 1.0e-8_dp) &
            print '(i5,f22.12,es14.4,es14.4,f16.10)', it, etot, etot - eold, drms, nelec_grid
         if (it > 1 .and. abs(etot - eold) < 1.0e-11_dp .and. drms < 1.0e-8_dp) then
            converged = .true.
            exit
         end if
         eold = etot
      end do
      if (.not. converged) then
         print '(a,i0,a)', "  NOT CONVERGED in ", MAXITER, " iterations"
         all_ok = .false.
      end if
      print '(a,a,a,f22.12)', "  E(", trim(names(ifn)), ") = ", etot
      energies(ifn) = etot
      write (unit) etot
      write (unit) dold
   end do
   close (unit)

   call eri%release(); call xg%release(); call grid%destroy(); call pl%release(); call bas%release()
   if (.not. all_ok) stop 1
   print '(a)', ""
   print '(a)', "scf_rks: converged; run rks_ref.py for the pyscf comparison"

contains

   integer function centre_of(ish)
      integer, intent(in) :: ish
      integer :: a
      centre_of = 0
      do a = 1, natm
         if (all(abs(sh_r(:, ish) - at_r(:, a)) < 1.0e-12_dp)) centre_of = a
      end do
      if (centre_of == 0) then
         print '(a)', "scf_rks: a shell is on no atom"; stop 1
      end if
   end function centre_of

   real(dp) function nuclear_repulsion(nat, z, r) result(e)
      integer, intent(in) :: nat
      real(dp), intent(in) :: z(nat), r(3, nat)
      integer :: a, b
      e = 0.0_dp
      do a = 1, nat
         do b = a + 1, nat
            e = e + z(a)*z(b)/norm2(r(:, a) - r(:, b))
         end do
      end do
   end function nuclear_repulsion

   ! X = S^(-1/2)
   subroutine sym_orthog(n, s, x)
      integer, intent(in) :: n
      real(dp), intent(in) :: s(n, n)
      real(dp), intent(out) :: x(n, n)
      real(dp), allocatable :: u(:, :), w(:), work(:)
      integer :: info, lwork, k
      allocate (u(n, n), w(n), work(1))
      u = s
      call dsyev('V', 'U', n, u, n, w, work, -1, info)
      lwork = int(work(1))
      deallocate (work); allocate (work(lwork))
      call dsyev('V', 'U', n, u, n, w, work, lwork, info)
      if (info /= 0) then
         print '(a,i0)', '  dsyev failed in sym_orthog: ', info; stop 1
      end if
      do k = 1, n
         u(:, k) = u(:, k)/sqrt(sqrt(w(k)))
      end do
      x = matmul(u, transpose(u))
      deallocate (u, w, work)
   end subroutine sym_orthog

   subroutine diag_fock(n, f, x, c, eps)
      integer, intent(in) :: n
      real(dp), intent(in) :: f(n, n), x(n, n)
      real(dp), intent(out) :: c(n, n), eps(n)
      real(dp), allocatable :: fp(:, :), work(:)
      integer :: info, lwork
      allocate (fp(n, n), work(1))
      fp = matmul(transpose(x), matmul(f, x))
      call dsyev('V', 'U', n, fp, n, eps, work, -1, info)
      lwork = int(work(1))
      deallocate (work); allocate (work(lwork))
      call dsyev('V', 'U', n, fp, n, eps, work, lwork, info)
      if (info /= 0) then
         print '(a,i0)', '  dsyev failed in diag_fock: ', info; stop 1
      end if
      c = matmul(x, fp)
      deallocate (fp, work)
   end subroutine diag_fock

   subroutine diis_push(n, f, e, fs, es, nused)
      integer, intent(in) :: n
      real(dp), intent(in) :: f(n, n), e(n, n)
      real(dp), intent(inout) :: fs(n, n, NDIIS), es(n, n, NDIIS)
      integer, intent(inout) :: nused
      integer :: k
      if (nused < NDIIS) then
         nused = nused + 1
      else
         do k = 1, NDIIS - 1
            fs(:, :, k) = fs(:, :, k + 1)
            es(:, :, k) = es(:, :, k + 1)
         end do
      end if
      fs(:, :, nused) = f
      es(:, :, nused) = e
   end subroutine diis_push

   subroutine diis_extrapolate(n, fs, es, nused, f)
      integer, intent(in) :: n, nused
      real(dp), intent(in) :: fs(n, n, NDIIS), es(n, n, NDIIS)
      real(dp), intent(out) :: f(n, n)
      real(dp), allocatable :: bmat(:, :), rhs(:)
      integer, allocatable :: ipiv(:)
      integer :: i, j, info, m
      m = nused + 1
      allocate (bmat(m, m), rhs(m), ipiv(m))
      bmat = -1.0_dp
      bmat(m, m) = 0.0_dp
      do j = 1, nused
         do i = 1, nused
            bmat(i, j) = sum(es(:, :, i)*es(:, :, j))
         end do
      end do
      rhs = 0.0_dp
      rhs(m) = -1.0_dp
      call dgesv(m, 1, bmat, m, ipiv, rhs, m, info)
      if (info /= 0) then
         f = fs(:, :, nused)
         return
      end if
      f = 0.0_dp
      do i = 1, nused
         f = f + rhs(i)*fs(:, :, i)
      end do
   end subroutine diis_extrapolate

   subroutine make_density(n, nocc, c, d)
      integer, intent(in) :: n, nocc
      real(dp), intent(in) :: c(n, n)
      real(dp), intent(out) :: d(n, n)
      d = 2.0_dp*matmul(c(:, 1:nocc), transpose(c(:, 1:nocc)))
   end subroutine make_density

end program scf_rks
