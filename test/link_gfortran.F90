!
! terco driven from GFORTRAN, through the C ABI only.
!
! WHY THIS TEST EXISTS
! --------------------
! mqc is built with gfortran and terco with nvfortran. Fortran modules cannot
! cross that boundary -- the `.mod` formats are unrelated -- so the only usable
! interface is `trc_capi.F90`'s `bind(c)` skin, and the only proof that it
! works is a gfortran program that links `libterco.so` and gets the right
! number out.
!
! It runs a complete RHF SCF: basis, pair list, Schwarz bounds, one-electron
! matrices and every Fock build come from terco across the C boundary, and
! only the eigenproblem is local. Water at 6-31G must give
!
!     -75.983974469890
!
! which is what `test/scf_rhf.F90` gets natively and what pyscf gets. Anything
! else means the ABI is misdescribing an array somewhere -- and a transposed
! or mis-strided matrix usually still converges, to the wrong energy, so the
! energy is the check and "it ran" is not.
!
! It also exercises EVERY entry point in the ABI, not just the ones an SCF
! needs, so that a signature change anywhere is caught here.
!
! Build (see the `link-test` target):
!   gfortran ... -L<build> -lterco
!
program link_gfortran
   use, intrinsic :: iso_c_binding
   !
   ! The interfaces come from terco, not from a copy written here. A copy is
   ! a contract kept in step by hand, and nothing catches it drifting: add an
   ! argument to `trc_capi.F90`, forget the copy, and it still compiles and
   ! still links, because a C symbol carries no signature. It fails at run
   ! time by reading whatever was in the register.
   !
   ! Compiling the shipped declaration here, with gfortran, and then asserting
   ! an energy is what makes that arrangement self-checking -- the drift shows
   ! up as a wrong number in terco's own test suite rather than in a host
   ! code months later.
   !
   use trc_c_interfaces
   use trc_test_basis, only: read_xyz, build_631g
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   character(len=256) :: xyzfile
   integer :: nat, nsh, maxnp, i, j, it, nocc, rc, nmult
   integer(c_int) :: nao_c
   integer,  allocatable :: at_z(:), sh_l(:), sh_np(:)
   real(dp), allocatable :: at_r(:, :), sh_e(:, :), sh_c(:, :), sh_r(:, :)
   real(c_double), allocatable :: zat(:)
   integer(c_int), allocatable :: lc(:), npc(:)
   type(c_ptr) :: hbas, hpair, heri, hbas2, hpair2
   integer :: nao
   real(dp), allocatable :: s(:, :), t(:, :), v(:, :), h(:, :)
   real(dp), allocatable :: d(:, :), g(:, :), f(:, :), x(:, :), c(:, :), e(:)
   real(dp), allocatable :: mm(:, :, :)
   real(dp) :: enuc, escf, eold, thresh

   thresh = 1.0e-10_dp
   if (command_argument_count() < 1) then
      print '(a)', 'usage: link_gfortran <xyz>'; stop 1
   end if
   call get_command_argument(1, xyzfile)

   call read_xyz(xyzfile, nat, at_z, at_r)
   call build_631g(nat, at_z, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   allocate (zat(nat), lc(nsh), npc(nsh))
   zat = real(at_z, c_double)
   lc = int(sh_l, c_int); npc = int(sh_np, c_int)

   rc = trc_basis_create(int(nsh, c_int), int(maxnp, c_int), lc, npc, &
                         sh_e, sh_c, sh_r, int(nat, c_int), zat, at_r, hbas)
   call must(rc, 'trc_basis_create')
   rc = trc_basis_nao(hbas, nao_c); call must(rc, 'trc_basis_nao')
   nao = int(nao_c)
   rc = trc_pairs_create(hbas, thresh, hpair); call must(rc, 'trc_pairs_create')
   rc = trc_eri_create(hbas, thresh, heri);    call must(rc, 'trc_eri_create')

   allocate (s(nao, nao), t(nao, nao), v(nao, nao), h(nao, nao))
   rc = trc_compute_1e(hbas, hpair, s, t, v); call must(rc, 'trc_compute_1e')
   h = t + v

   enuc = 0.0_dp
   do i = 2, nat
      do j = 1, i - 1
         enuc = enuc + zat(i)*zat(j)/norm2(at_r(:, i) - at_r(:, j))
      end do
   end do

   nocc = sum(at_z)/2
   allocate (d(nao, nao), g(nao, nao), f(nao, nao), x(nao, nao))
   allocate (c(nao, nao), e(nao))
   call sym_orthog(nao, s, x)
   call diag(nao, h, x, c, e)
   d = 2.0_dp*matmul(c(:, 1:nocc), transpose(c(:, 1:nocc)))
   eold = 0.0_dp
   do it = 1, 200
      rc = trc_fock(heri, hbas, d, g, 1.0_c_double, 1.0_c_double, 1_c_int)
      call must(rc, 'trc_fock')
      f = h + g
      escf = enuc + 0.5_dp*sum(d*(h + f))
      call diag(nao, f, x, c, e)
      d = 2.0_dp*matmul(c(:, 1:nocc), transpose(c(:, 1:nocc)))
      if (abs(escf - eold) < 1.0e-11_dp) exit
      eold = escf
   end do

   call cover_rest()

   nmult = int(trc_multipole_count())
   allocate (mm(nao, nao, nmult))
   rc = trc_compute_multipoles(hbas, hpair, [0.0_c_double, 0.0_c_double, &
                               0.0_c_double], mm)
   call must(rc, 'trc_compute_multipoles')

   print '(a)', '  terco through the C ABI, from gfortran'
   print '(a,i0,a,i0)', '    nao ', nao, '   iterations ', it
   print '(a,f18.10)',  '    RHF energy   ', escf
   print '(a,f18.10)',  '    sum(S)       ', sum(s)
   print '(a,3f14.8)',  '    dipole <r>   ', -sum(d*mm(:, :, 1)), &
                                             -sum(d*mm(:, :, 2)), &
                                             -sum(d*mm(:, :, 3))
   ! The same SCF through the library's own driver, across the same ABI.
   ! It must land on the energy the loop above found, and PBE must converge.
   block
      real(c_double) :: e_lib, e_xc
      real(c_double), allocatable :: d_lib(:, :), eps_lib(:)
      integer(c_int) :: niter
      allocate (d_lib(nao, nao), eps_lib(nao))
      rc = trc_scf(hbas, int(nocc, c_int), int(nocc, c_int), c_null_char, 3_c_int, &
                   1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, &
                   e_lib, e_xc, d_lib, eps_lib, niter)
      call must(rc, 'trc_scf (HF)')
      print '(a,f18.10,a,i0,a)', '    trc_scf RHF  ', e_lib, '   (', niter, ' iterations)'
      if (abs(e_lib - escf) > 1.0e-8_dp) then
         print '(a)', '    RESULT: FAIL (trc_scf disagrees with the loop above)'
         stop 1
      end if
      rc = trc_scf(hbas, int(nocc, c_int), int(nocc, c_int), 'pbe'//c_null_char, 3_c_int, &
                   1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, &
                   e_lib, e_xc, d_lib, eps_lib, niter)
      call must(rc, 'trc_scf (PBE)')
      print '(a,f18.10,a,f16.10,a,i0,a)', '    trc_scf PBE  ', e_lib, '   E_xc ', e_xc, '   (', niter, ' iterations)'
   end block
   if (abs(escf + 75.983974469890_dp) < 1.0e-8_dp) then
      print '(a)', '    RESULT: PASS (matches the native build and pyscf)'
   else
      print '(a)', '    RESULT: FAIL'
      stop 1
   end if

   rc = trc_eri_destroy(heri); rc = trc_pairs_destroy(hpair)
   rc = trc_basis_destroy(hbas)

contains

   !
   ! Every remaining entry point, each with an invariant that a wrong
   ! signature would break.
   !
   ! Coverage is the point. An SCF alone touches eight of the eighteen calls,
   ! so a changed argument list on any of the other ten would sail through
   ! this test and surface in a host code instead.
   !
   subroutine cover_rest()
      integer(c_int) :: maxl, npair
      real(dp), allocatable :: g1(:, :), gn(:, :), gm(:, :, :)
      real(dp), allocatable :: j2c(:, :), t3c(:, :, :)
      real(dp) :: worst
      integer :: i, j

      ! --- basis_maxl: 6-31G for H and O is s and p, nothing higher --------
      rc = trc_basis_maxl(hbas, maxl); call must(rc, 'trc_basis_maxl')
      if (maxl /= 1) then
         print '(a,i0)', '    basis_maxl should be 1 for 6-31G, got ', maxl
         stop 1
      end if

      ! --- pairs_count -----------------------------------------------------
      rc = trc_pairs_count(hpair, npair); call must(rc, 'trc_pairs_count')
      if (npair <= 0) then
         print '(a)', '    pairs_count returned nothing'; stop 1
      end if

      ! --- the libcint-layout constructor ----------------------------------
      !
      ! This is the path mqc uses, so it gets the strongest check available:
      ! build the SAME basis both ways and require the same overlap. A
      ! mistake in the packed layout -- a slot off by one, a pointer that is
      ! 1-based where libcint is 0-based -- changes the basis rather than
      ! failing, so comparing against the array constructor is what catches
      ! it.
      call build_libcint_arrays()
      rc = trc_basis_nao(hbas2, nao_c); call must(rc, 'nao of libcint basis')
      if (int(nao_c) /= nao) then
         print '(a,i0,a,i0)', '    libcint-layout basis has nao ', nao_c, &
            ' against ', nao
         stop 1
      end if
      block
         real(dp), allocatable :: s2(:, :), t2(:, :), v2(:, :)
         allocate (s2(nao, nao), t2(nao, nao), v2(nao, nao))
         rc = trc_compute_1e(hbas2, hpair2, s2, t2, v2)
         call must(rc, '1e over the libcint-layout basis')
         worst = maxval(abs(s2 - s))
         print '(a,es10.2)', '    libcint-layout basis: overlap agrees to ', worst
         if (worst > 1.0e-12_dp) then
            print '(a)', '    the two constructors describe different bases'
            stop 1
         end if
         deallocate (s2, t2, v2)
      end block

      ! --- fock_many against fock ------------------------------------------
      !
      ! Two copies of the converged density must give two copies of the
      ! answer `fock` already gave. That checks the batched entry point AND
      ! the (ndens, nao, nao) layout at once: with the density index in the
      ! wrong position this reads across the wrong stride and disagrees.
      allocate (g1(nao, nao), gn(nao, nao), gm(2, nao, nao))
      rc = trc_fock(heri, hbas, d, g1, 1.0_c_double, 1.0_c_double, 1_c_int)
      call must(rc, 'trc_fock')
      block
         real(dp), allocatable :: dm(:, :, :)
         allocate (dm(2, nao, nao))
         do j = 1, nao
            do i = 1, nao
               dm(1, i, j) = d(i, j)
               dm(2, i, j) = d(i, j)
            end do
         end do
         rc = trc_fock_many(heri, hbas, 2_c_int, dm, gm, 1.0_c_double, &
                            1.0_c_double, 1_c_int)
         call must(rc, 'trc_fock_many')
         worst = 0.0_dp
         do j = 1, nao
            do i = 1, nao
               worst = max(worst, abs(gm(1, i, j) - g1(i, j)))
               worst = max(worst, abs(gm(2, i, j) - g1(i, j)))
            end do
         end do
         deallocate (dm)
      end block
      print '(a,es10.2)', '    fock_many vs fock                    ', worst
      if (worst > 1.0e-10_dp) then
         print '(a)', '    the batched build disagrees with the single one'
         stop 1
      end if

      ! --- fock_nosym against fock -----------------------------------------
      !
      ! On a SYMMETRIC density the folded and enumerated digestions must
      ! agree. That is the control that makes a disagreement on an
      ! asymmetric density mean the asymmetry rather than a bug.
      rc = trc_fock_nosym(heri, hbas, d, gn, 1_c_int)
      call must(rc, 'trc_fock_nosym')
      worst = maxval(abs(gn - g1))
      print '(a,es10.2)', '    fock_nosym vs fock (symmetric D)     ', worst
      print '(a,es24.16,a,es24.16)', '    sum folded ', sum(g1), &
         '   sum enumerated ', sum(gn)
      print '(a,f18.10)', '    ratio ', sum(gn)/sum(g1)
      if (worst > 1.0e-10_dp) then
         print '(a)', '    the enumerated digestion disagrees on a symmetric density'
         stop 1
      end if

      ! --- density fitting, with the orbital basis as its own auxiliary ----
      allocate (j2c(nao, nao), t3c(nao, nao, nao))
      rc = trc_compute_df2c(hbas, j2c);              call must(rc, 'trc_compute_df2c')
      rc = trc_compute_df3c(hbas, hpair, hbas, t3c); call must(rc, 'trc_compute_df3c')
      ! (P|Q) is symmetric and positive definite, so its diagonal is positive
      if (.not. ieee_ok(j2c) .or. minval([(j2c(i, i), i=1, nao)]) <= 0.0_dp) then
         print '(a)', '    the two-centre metric is not positive on its diagonal'
         stop 1
      end if
      print '(a,f16.8,a,f16.8)', '    df 2c trace ', &
         sum([(j2c(i, i), i=1, nao)]), '   3c sum ', sum(t3c)

      rc = trc_pairs_destroy(hpair2); rc = trc_basis_destroy(hbas2)
      deallocate (g1, gn, gm, j2c, t3c)
   end subroutine cover_rest

   logical function ieee_ok(a)
      real(dp), intent(in) :: a(:, :)
      ieee_ok = all(a == a) .and. all(abs(a) < huge(1.0_dp))
   end function ieee_ok

   !
   ! Pack the same basis into libcint's atm/bas/env.
   !
   ! The slot numbers are libcint's, 0-based, so a Fortran array indexed from
   ! one carries them at index+1. `env` pointers are 0-based offsets into
   ! `env` itself, which is the part that is easy to get wrong from Fortran
   ! and the reason the overlap is compared afterwards rather than trusted.
   !
   subroutine build_libcint_arrays()
      integer, parameter :: PTR_ENV_START = 20
      integer(c_int), allocatable :: atm(:, :), bas(:, :)
      real(c_double), allocatable :: env(:)
      integer :: off, ia, is, k, ncoord

      ncoord = 3*nat
      allocate (atm(6, nat), bas(8, nsh))
      allocate (env(PTR_ENV_START + ncoord + 2*sum(sh_np)))
      env = 0.0_c_double
      atm = 0_c_int; bas = 0_c_int

      off = PTR_ENV_START
      do ia = 1, nat
         atm(1, ia) = int(at_z(ia), c_int)      ! CHARGE_OF
         atm(2, ia) = int(off, c_int)           ! PTR_COORD
         env(off + 1:off + 3) = real(at_r(:, ia), c_double)
         off = off + 3
      end do

      do is = 1, nsh
         bas(1, is) = int(atom_of(is) - 1, c_int)   ! ATOM_OF, 0-based
         bas(2, is) = int(sh_l(is), c_int)          ! ANG_OF
         bas(3, is) = int(sh_np(is), c_int)         ! NPRIM_OF
         bas(4, is) = 1_c_int                       ! NCTR_OF
         bas(5, is) = 0_c_int                       ! KAPPA_OF
         bas(6, is) = int(off, c_int)               ! PTR_EXP
         do k = 1, sh_np(is)
            env(off + k) = real(sh_e(k, is), c_double)
         end do
         off = off + sh_np(is)
         bas(7, is) = int(off, c_int)               ! PTR_COEFF
         do k = 1, sh_np(is)
            env(off + k) = real(sh_c(k, is), c_double)
         end do
         off = off + sh_np(is)
      end do

      rc = trc_basis_create_libcint(atm, int(nat, c_int), bas, &
                                    int(nsh, c_int), env, &
                                    int(size(env), c_int), 1_c_int, hbas2)
      call must(rc, 'trc_basis_create_libcint')
      rc = trc_pairs_create(hbas2, thresh, hpair2)
      call must(rc, 'pairs over the libcint-layout basis')
      deallocate (atm, bas, env)
   end subroutine build_libcint_arrays

   !> Which atom a shell sits on, by matching its centre.
   integer function atom_of(is)
      integer, intent(in) :: is
      integer :: ia
      atom_of = 1
      do ia = 1, nat
         if (maxval(abs(sh_r(:, is) - at_r(:, ia))) < 1.0e-12_dp) atom_of = ia
      end do
   end function atom_of

   subroutine must(code, what)
      integer, intent(in) :: code
      character(len=*), intent(in) :: what
      if (code /= 0) then
         print '(a,a,a,i0)', '  ', what, ' failed with status ', code
         stop 1
      end if
   end subroutine must

   subroutine sym_orthog(n, sm, xo)
      integer, intent(in) :: n
      real(dp), intent(in) :: sm(n, n)
      real(dp), intent(out) :: xo(n, n)
      real(dp), allocatable :: u(:, :), w(:), work(:)
      integer :: info, lw, ii, jj
      allocate (u(n, n), w(n), work(64*n)); u = sm; lw = 64*n
      call dsyev('V', 'U', n, u, n, w, work, lw, info)
      do jj = 1, n
         do ii = 1, n
            xo(ii, jj) = u(ii, jj)/sqrt(w(jj))
         end do
      end do
      xo = matmul(xo, transpose(u))
   end subroutine sym_orthog

   subroutine diag(n, fm, xo, cm, ev)
      integer, intent(in) :: n
      real(dp), intent(in) :: fm(n, n), xo(n, n)
      real(dp), intent(out) :: cm(n, n), ev(n)
      real(dp), allocatable :: fp(:, :), work(:)
      integer :: info, lw
      allocate (fp(n, n), work(64*n)); lw = 64*n
      fp = matmul(transpose(xo), matmul(fm, xo))
      call dsyev('V', 'U', n, fp, n, ev, work, lw, info)
      cm = matmul(xo, fp)
   end subroutine diag

end program link_gfortran
