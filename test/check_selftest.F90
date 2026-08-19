!
! Correctness with no external program and no recorded data.
!
! WHY THIS EXISTS ALONGSIDE THE ORACLES
! -------------------------------------
! The other checks compare against libfint's numbers recorded in a file, or
! against pyscf. Those are the strongest evidence available -- an independent
! implementation is the only thing that catches a shared misunderstanding --
! but they make the suite depend on a 6 MB blob and on somebody else's package
! being installed, and they cannot be run by someone who has neither.
!
! These can. Every check here is either a closed-form value or an identity the
! integrals must satisfy whatever they are, so a failure is a real defect and
! not a tolerance argument.
!
!   1. ANALYTIC OVERLAP. Two s primitives have an overlap in closed form.
!      This pins the normalisation convention as well as the arithmetic --
!      which is the part most likely to be quietly wrong, and the part an
!      identity cannot see.
!
!   2. TRANSLATIONAL INVARIANCE. Move every centre by the same vector and
!      every integral must be unchanged. This is a strong check for its
!      price: it exercises the whole path -- pair list, screening, bins,
!      kernels, digestion -- and it is sensitive to any place a coordinate is
!      used where a difference of coordinates was meant. Absolute positions
!      appear in a Gaussian product centre, so getting this right is not
!      automatic.
!
!   3. THE MULTIPOLE ORIGIN SHIFT. Moving the expansion origin by `a` must
!      take the dipole to `mu - a*S` exactly. This checks the multipole
!      kernel against the OVERLAP kernel -- two independently generated code
!      paths that must agree by an identity -- with nothing external at all.
!
!   4. SYMMETRY OF THE FOCK MATRIX. J and K built from a symmetric density
!      are symmetric. Cheap, and it catches a scatter that writes one of the
!      six folded updates to a transposed address.
!
!   5. THE BATCH AGAINST THE LOOP. N densities in one integral pass must give
!      what N separate passes give, set for set.
!
program check_selftest
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_multipoles, &
                      TRC_NMULT
   use trc_fock, only: trc_eri_t
   implicit none

   real(dp), parameter :: PI = 3.14159265358979323846_dp
   !> libcint's per-shell factor, which `basis%build` folds into the
   !> coefficients. The analytic comparison has to include it or it is
   !> measuring a different function.
   real(dp), parameter :: FAC_S = 0.282094791773878143_dp

   integer :: nbad
   nbad = 0

   call test_analytic_overlap()
   call test_translation()
   call test_multipole_shift()
   call test_fock_symmetry()
   call test_batch_matches_loop()

   print '(a)', ''
   if (nbad == 0) then
      print '(a)', '  RESULT: PASS'
   else
      print '(a,i0,a)', '  RESULT: FAIL (', nbad, ' checks outside tolerance)'
      stop 1
   end if

contains

   subroutine report(what, diff, tol)
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: diff, tol
      character(len=4) :: verdict
      verdict = 'ok  '
      if (.not. (diff <= tol)) then
         verdict = 'FAIL'
         nbad = nbad + 1
      end if
      print '(a,a40,a,es10.2,a,es8.1,a)', '  ', what, '  ', diff, &
         '   (tol ', tol, ')   '//verdict
   end subroutine report

   !> libcint's primitive normalisation, so a coefficient of one means what a
   !> basis-set file means by one.
   pure real(dp) function gto_norm(l, a)
      integer,  intent(in) :: l
      real(dp), intent(in) :: a
      real(dp) :: nn, gi
      nn = real(2*l + 2, dp)
      gi = 0.5_dp*(2.0_dp*a)**(-(nn + 1.0_dp)/2.0_dp)*gamma((nn + 1.0_dp)/2.0_dp)
      gto_norm = 1.0_dp/sqrt(gi)
   end function gto_norm

   !
   ! A small system, built from scratch so nothing is read from disk.
   ! `shift` displaces every centre, which is what test 2 needs.
   !
   subroutine build_system(b, p, shift, uncontracted_s_only)
      type(trc_basis_t), intent(out) :: b
      type(trc_pairlist_t), intent(out) :: p
      real(dp), intent(in) :: shift(3)
      logical, intent(in) :: uncontracted_s_only

      integer, parameter :: NATM = 3
      integer :: nsh, i, k
      integer,  allocatable :: sh_l(:), sh_np(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
      real(dp) :: at_r(3, NATM), at_z(NATM)

      at_r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp] + shift
      at_r(:, 2) = [1.3_dp, 0.0_dp, 0.7_dp] + shift
      at_r(:, 3) = [-0.4_dp, 1.1_dp, 1.9_dp] + shift
      at_z = [4.0_dp, 3.0_dp, 2.0_dp]

      if (uncontracted_s_only) then
         nsh = NATM
      else
         nsh = 2*NATM          ! one s and one p per centre
      end if
      allocate (sh_l(nsh), sh_np(nsh), sh_e(1, nsh), sh_c(1, nsh), sh_r(3, nsh))
      sh_np = 1
      k = 0
      do i = 1, NATM
         k = k + 1
         sh_l(k) = 0; sh_e(1, k) = 0.8_dp + 0.19_dp*real(i, dp)
         sh_r(:, k) = at_r(:, i)
         if (.not. uncontracted_s_only) then
            k = k + 1
            sh_l(k) = 1; sh_e(1, k) = 1.3_dp + 0.11_dp*real(i, dp)
            sh_r(:, k) = at_r(:, i)
         end if
      end do
      do i = 1, nsh
         sh_c(1, i) = gto_norm(sh_l(i), sh_e(1, i))
      end do

      call b%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, NATM, at_z, at_r, 1)
      call p%build(b, 1.0e-14_dp)
      call b%to_device(); call p%to_device()
   end subroutine build_system

   ! =======================================================================
   ! 1. the closed form
   ! =======================================================================
   subroutine test_analytic_overlap()
      type(trc_basis_t) :: b
      type(trc_pairlist_t) :: p
      real(dp), allocatable :: s(:, :), t(:, :), v(:, :)
      real(dp) :: ea, eb, r2, want, worst
      real(dp) :: ra(3), rb(3)
      integer :: n, i, j
      real(dp), allocatable :: ex(:), ce(:)

      call build_system(b, p, [0.0_dp, 0.0_dp, 0.0_dp], .true.)
      n = b%nao
      allocate (s(n, n), t(n, n), v(n, n))
      !$acc enter data create(s, t, v)
      call trc_1e(b, p, s, t, v)
      !$acc update self(s)
      !$acc exit data delete(s, t, v)

      ! S(a,b) = Na Nb (fac_s)^2 (pi/(ea+eb))^{3/2} exp(-ea eb R^2/(ea+eb))
      allocate (ex(n), ce(n))
      do i = 1, n
         ex(i) = b%sh_e(1, i)
         ce(i) = gto_norm(0, ex(i))
      end do
      worst = 0.0_dp
      do j = 1, n
         do i = 1, n
            ra = b%sh_r(:, i); rb = b%sh_r(:, j)
            ea = ex(i); eb = ex(j)
            r2 = sum((ra - rb)**2)
            want = ce(i)*ce(j)*FAC_S*FAC_S &
                   *(PI/(ea + eb))**1.5_dp*exp(-ea*eb*r2/(ea + eb))
            worst = max(worst, abs(s(i, j) - want))
         end do
      end do
      call report('overlap vs the closed form (s|s)', worst, 1.0e-14_dp)

      call p%release(); call b%release()
      deallocate (s, t, v, ex, ce)
   end subroutine test_analytic_overlap

   ! =======================================================================
   ! 2. move everything, change nothing
   ! =======================================================================
   subroutine test_translation()
      type(trc_basis_t) :: b1, b2
      type(trc_pairlist_t) :: p1, p2
      type(trc_eri_t) :: e1, e2
      real(dp), allocatable :: s1(:, :), t1(:, :), v1(:, :)
      real(dp), allocatable :: s2(:, :), t2(:, :), v2(:, :)
      real(dp), allocatable :: d(:, :), g1(:, :), g2(:, :)
      real(dp), parameter :: SHIFT(3) = [3.7_dp, -2.1_dp, 5.4_dp]
      integer :: n, i, j

      call build_system(b1, p1, [0.0_dp, 0.0_dp, 0.0_dp], .false.)
      call build_system(b2, p2, SHIFT, .false.)
      n = b1%nao
      allocate (s1(n, n), t1(n, n), v1(n, n), s2(n, n), t2(n, n), v2(n, n))
      !$acc enter data create(s1, t1, v1, s2, t2, v2)
      call trc_1e(b1, p1, s1, t1, v1)
      call trc_1e(b2, p2, s2, t2, v2)
      !$acc update self(s1, t1, v1, s2, t2, v2)
      !$acc exit data delete(s1, t1, v1, s2, t2, v2)

      call report('overlap under translation', maxval(abs(s1 - s2)), 1.0e-13_dp)
      call report('kinetic under translation', maxval(abs(t1 - t2)), 1.0e-13_dp)
      ! Nuclear attraction moves with the nuclei, so it is invariant too --
      ! and it is the one that would break if an absolute coordinate leaked in
      ! where a difference was meant.
      call report('nuclear attraction under translation', &
                  maxval(abs(v1 - v2)), 1.0e-12_dp)

      call e1%build(b1, 1.0e-14_dp)
      call e2%build(b2, 1.0e-14_dp)
      allocate (d(n, n), g1(n, n), g2(n, n))
      do j = 1, n
         do i = 1, n
            d(i, j) = 0.4_dp*cos(real(3*i + 5*j, dp)) + 0.2_dp*real(min(i, j), dp)
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
      call e1%fock(b1, d, g1)
      call e2%fock(b2, d, g2)
      call report('Fock matrix under translation', maxval(abs(g1 - g2)), 1.0e-12_dp)

      call e1%release(); call e2%release()
      call p1%release(); call p2%release(); call b1%release(); call b2%release()
      deallocate (s1, t1, v1, s2, t2, v2, d, g1, g2)
   end subroutine test_translation

   ! =======================================================================
   ! 3. the multipole kernel against the overlap kernel
   ! =======================================================================
   subroutine test_multipole_shift()
      type(trc_basis_t) :: b
      type(trc_pairlist_t) :: p
      real(dp), allocatable :: s(:, :), t(:, :), v(:, :)
      real(dp), allocatable :: m0(:, :, :), m1(:, :, :)
      real(dp), parameter :: A(3) = [0.9_dp, -1.4_dp, 2.2_dp]
      real(dp) :: worst
      integer :: n, c, i, j

      call build_system(b, p, [0.0_dp, 0.0_dp, 0.0_dp], .false.)
      n = b%nao
      allocate (s(n, n), t(n, n), v(n, n))
      allocate (m0(n, n, TRC_NMULT), m1(n, n, TRC_NMULT))
      !$acc enter data create(s, t, v, m0, m1)
      call trc_1e(b, p, s, t, v)
      call trc_multipoles(b, p, [0.0_dp, 0.0_dp, 0.0_dp], m0)
      call trc_multipoles(b, p, A, m1)
      !$acc update self(s, m0, m1)
      !$acc exit data delete(s, t, v, m0, m1)

      ! <i| (r - a)_c |j>  =  <i| r_c |j>  -  a_c <i|j>
      worst = 0.0_dp
      do c = 1, 3
         do j = 1, n
            do i = 1, n
               worst = max(worst, abs(m1(i, j, c) - (m0(i, j, c) - A(c)*s(i, j))))
            end do
         end do
      end do
      call report('dipole origin shift vs the overlap', worst, 1.0e-13_dp)

      call p%release(); call b%release()
      deallocate (s, t, v, m0, m1)
   end subroutine test_multipole_shift

   ! =======================================================================
   ! 4 and 5
   ! =======================================================================
   subroutine test_fock_symmetry()
      type(trc_basis_t) :: b
      type(trc_pairlist_t) :: p
      type(trc_eri_t) :: e
      real(dp), allocatable :: d(:, :), g(:, :)
      integer :: n, i, j

      call build_system(b, p, [0.0_dp, 0.0_dp, 0.0_dp], .false.)
      call e%build(b, 1.0e-14_dp)
      n = b%nao
      allocate (d(n, n), g(n, n))
      do j = 1, n
         do i = 1, n
            d(i, j) = 0.3_dp*sin(real(2*i - 7*j, dp)) + 0.15_dp*real(i + j, dp)
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
      call e%fock(b, d, g)
      call report('G is symmetric for a symmetric density', &
                  maxval(abs(g - transpose(g))), 1.0e-13_dp)

      call e%release(); call p%release(); call b%release()
      deallocate (d, g)
   end subroutine test_fock_symmetry

   subroutine test_batch_matches_loop()
      type(trc_basis_t) :: b
      type(trc_pairlist_t) :: p
      type(trc_eri_t) :: e
      integer, parameter :: ND = 3
      real(dp), allocatable :: dm(:, :, :), gm(:, :, :), d1(:, :), g1(:, :)
      real(dp) :: worst
      integer :: n, i, j, k

      call build_system(b, p, [0.0_dp, 0.0_dp, 0.0_dp], .false.)
      call e%build(b, 1.0e-14_dp)
      n = b%nao
      allocate (dm(ND, n, n), gm(ND, n, n), d1(n, n), g1(n, n))
      do k = 1, ND
         do j = 1, n
            do i = 1, n
               dm(k, i, j) = 0.2_dp*cos(real(i + 2*j + 11*k, dp))
            end do
         end do
         dm(k, :, :) = 0.5_dp*(dm(k, :, :) + transpose(dm(k, :, :)))
      end do
      call e%fock_many(b, ND, dm, gm)

      worst = 0.0_dp
      do k = 1, ND
         d1 = dm(k, :, :)
         call e%fock(b, d1, g1)
         do j = 1, n
            do i = 1, n
               worst = max(worst, abs(gm(k, i, j) - g1(i, j)))
            end do
         end do
      end do
      call report('batched build vs the same densities singly', worst, 1.0e-12_dp)

      call e%release(); call p%release(); call b%release()
      deallocate (dm, gm, d1, g1)
   end subroutine test_batch_matches_loop

end program check_selftest
