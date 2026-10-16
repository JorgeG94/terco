!
! The polarisable continuum's cavity, operators and charge solve, with no
! integrals: trc_pcm's half that does not touch the basis.
!
! Two kinds of evidence. The physics, ported from metalquicha's own PCM
! tests: an isolated sphere keeps every point and integrates to 4 pi R^2;
! overlapping spheres bury the points between them; a unit charge at the
! centre of a sphere is Born's problem, which C-PCM reproduces as
! -f/(2R) with its own f = (eps-1)/eps and IEF-PCM must reproduce as the
! TRUE Born energy -(eps-1)/(2 R eps) despite carrying (eps-1)/(eps+1)
! inside (only the D A S terms get it there); and the refusals -- an
! unknown model, a gas-phase dielectric, a Lebedev order with no fitted
! exponent, an element with no radius.
!
! And pyscf's own arrays, from test/pyscf_refs/pcm_cavity_water.ref
! (generate.py's gen_pcm_cavity, pyscf.solvent.pcm on water.xyz with the
! same radii, scale and Lebedev orders): at 302 points per sphere the
! point set, weights, exponents, switching values and areas, S x and D x
! for a fixed x, and the symmetrised charges and q.v/2 of both models for a
! fixed potential; at 110 points per sphere S and D element by element.
! pyscf's points come in its own Lebedev order, so each is matched to
! terco's by position (1e-10 bohr) before anything is compared; the fixed
! x and v are closed-form functions of the position, so they need no
! matching of their own.
!
program check_pcm_cavity
   use pic_types, only: dp
   use trc_error, only: error_t
   use trc_pcm, only: trc_pcm_t, pcm_surface_t, build_pcm_surface, pcm_sd_matrices, &
                      pcm_switch_h, pcm_model_from_name, PCM_MODEL_CPCM, PCM_MODEL_IEFPCM
   use trc_pcm_radii, only: pcm_cavity_radius
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, &
                           ref_read_vector, ref_read_matrix, ref_report
   implicit none

   real(dp), parameter :: PI = 3.14159265358979323846_dp
   real(dp), parameter :: EPS_WATER = 78.3553_dp
   real(dp), parameter :: TOL = 1.0e-12_dp
   !! Absolute, on everything but S. S's diagonal is zeta sqrt(2/pi)/F, and
   !! a point nearly switched off has F ~ 1e-6 and S_ii ~ 1e5, where the
   !! last bits of F's product of switches (1e-15 apart between the two
   !! codes) are 1e-8 absolute. So S, and S x which its diagonal dominates,
   !! are compared relative to max(1, |S|): the same 1e-12, scaled.
   real(dp), parameter :: ORIGIN(3, 1) = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])
   integer :: nfail, natm
   integer, allocatable :: zat(:)
   real(dp), allocatable :: rat(:, :)

   nfail = 0
   call physics()
   call refusals()
   call read_xyz('water.xyz', natm, zat, rat)
   call against_pyscf()
   if (nfail > 0) then
      print '(a,i0,a)', 'check_pcm_cavity: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_pcm_cavity: RESULT: PASS'

contains

   subroutine physics()
      type(pcm_surface_t) :: surf
      type(error_t) :: err
      real(dp) :: radius, e, qt, f, h

      ! The switching function's own shape.
      h = max(abs(pcm_switch_h(0.0_dp)), abs(pcm_switch_h(1.0_dp) - 1.0_dp), &
              abs(pcm_switch_h(0.5_dp) - 0.5_dp), abs(pcm_switch_h(0.3_dp) + pcm_switch_h(0.7_dp) - 1.0_dp), &
              abs(pcm_switch_h(-0.2_dp)), abs(pcm_switch_h(1.4_dp) - 1.0_dp))
      call ref_report(nfail, h < 1.0e-15_dp, 'switch_h: 0 and 1 at the ends, odd about 1/2', h)

      ! An isolated sphere keeps every point; its weights sum to one and
      ! nothing is switched, so the area is the sphere's to roundoff.
      call build_pcm_surface([8], ORIGIN, 110, 1.2_dp, surf, err)
      call pcm_cavity_radius(8, 1.2_dp, radius, err)
      call ref_report(nfail, .not. err%has_error() .and. surf%n_points == 110, &
                      'isolated sphere keeps all 110 points', real(110 - surf%n_points, dp))
      e = abs(sum(surf%area) - 4.0_dp*PI*radius**2)
      call ref_report(nfail, e < 1.0e-10_dp, 'isolated sphere area = 4 pi R^2', e)

      ! Two oxygens a bond apart lose the points between them.
      call build_pcm_surface([8, 8], reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.3_dp], [3, 2]), &
                             110, 1.2_dp, surf, err)
      call ref_report(nfail, .not. err%has_error() .and. surf%n_points < 220 .and. &
                      surf%n_points + surf%n_discarded == 220 .and. &
                      sum(surf%area) < 8.0_dp*PI*radius**2, &
                      'overlapping spheres bury points, kept + dropped = 220', real(surf%n_discarded, dp))

      ! Born: a unit charge at the centre of one oxygen cavity.
      f = (EPS_WATER - 1.0_dp)/EPS_WATER
      call born(PCM_MODEL_CPCM, e, qt)
      call ref_report(nfail, abs(e + 0.5_dp*f/radius) < 1.0e-5_dp*0.5_dp*f/radius, &
                      'Born ion, C-PCM: -f/(2R), relative', abs(e + 0.5_dp*f/radius)/(0.5_dp*f/radius))
      call ref_report(nfail, abs(qt + f) < 1.0e-4_dp, 'Born ion, C-PCM: Gauss, total charge -f', abs(qt + f))
      call born(PCM_MODEL_IEFPCM, e, qt)
      call ref_report(nfail, abs(e + 0.5_dp*f/radius) < 1.0e-4_dp*0.5_dp*f/radius, &
                      'Born ion, IEF-PCM: the exact -(eps-1)/(2 R eps), relative', &
                      abs(e + 0.5_dp*f/radius)/(0.5_dp*f/radius))
   end subroutine physics

   subroutine born(model, e, qt)
      integer, intent(in) :: model
      real(dp), intent(out) :: e, qt
      type(trc_pcm_t) :: ctx
      type(error_t) :: err
      real(dp), allocatable :: v(:), q(:)
      integer :: k
      call ctx%init(model, EPS_WATER, err)
      call build_pcm_surface([8], ORIGIN, 302, 1.2_dp, ctx%surface, err)
      call ctx%build_operators(err)
      if (err%has_error()) then
         print '(a)', 'check_pcm_cavity: '//err%get_message(); stop 1
      end if
      allocate (v(ctx%surface%n_points), q(ctx%surface%n_points))
      do k = 1, ctx%surface%n_points
         v(k) = erf(ctx%surface%zeta(k)*norm2(ctx%surface%coords(:, k)))/norm2(ctx%surface%coords(:, k))
      end do
      ! solve's vectors live where K does: on the device under the GPU build
      !$acc enter data copyin(v) create(q)
      call ctx%solve(v, q)
      !$acc update self(q)
      !$acc exit data delete(v, q)
      e = 0.5_dp*dot_product(q, v)
      qt = sum(q)
      call ctx%release()
   end subroutine born

   subroutine refusals()
      type(trc_pcm_t) :: ctx
      type(pcm_surface_t) :: surf
      type(error_t) :: e1, e2, e3, e4, e5
      call ctx%init(pcm_model_from_name('cosmo'), EPS_WATER, e1)
      call ctx%init(PCM_MODEL_CPCM, 1.0_dp, e2)
      call build_pcm_surface([8], ORIGIN, 74, 1.2_dp, surf, e3)
      call build_pcm_surface([26], ORIGIN, 302, 1.2_dp, surf, e4)
      call build_pcm_surface([8], ORIGIN, 303, 1.2_dp, surf, e5)
      call ref_report(nfail, e1%has_error(), "refuses 'cosmo' (a model this does not solve)", 0.0_dp)
      call ref_report(nfail, e2%has_error(), 'refuses eps = 1', 0.0_dp)
      call ref_report(nfail, e3%has_error(), 'refuses Lebedev 74 (no fitted exponent)', 0.0_dp)
      call ref_report(nfail, e4%has_error(), 'refuses iron (no cavity radius)', 0.0_dp)
      call ref_report(nfail, e5%has_error(), 'refuses 303 points (no such order)', 0.0_dp)
      call ref_report(nfail, pcm_model_from_name('IEF-PCM') == PCM_MODEL_IEFPCM .and. &
                      pcm_model_from_name('cpcm') == PCM_MODEL_CPCM, 'model names', 0.0_dp)
   end subroutine refusals

   pure real(dp) function vfix(s)
      !! The fixed potential: two inverse distances, a closed form of position
      real(dp), intent(in) :: s(3)
      vfix = 1.0_dp/norm2(s - [0.1_dp, 0.2_dp, 0.3_dp]) + 0.5_dp/norm2(s - [-0.5_dp, 0.4_dp, -0.2_dp])
   end function vfix

   pure real(dp) function xfix(s)
      !! The fixed vector S and D are applied to
      real(dp), intent(in) :: s(3)
      xfix = cos(1.3_dp*s(1)) + sin(0.7_dp*s(2)) + 0.2_dp*s(3)
   end function xfix

   subroutine match(n, mine, theirs, perm, worst)
      !! perm(k) = terco's point at pyscf's point k, by position
      integer, intent(in) :: n
      real(dp), intent(in) :: mine(3, n), theirs(3, n)
      integer, intent(out) :: perm(n)
      real(dp), intent(out) :: worst
      integer :: k, j, best
      real(dp) :: d, dbest
      logical :: used(n)
      used = .false.
      worst = 0.0_dp
      do k = 1, n
         best = 0; dbest = huge(1.0_dp)
         do j = 1, n
            d = norm2(mine(:, j) - theirs(:, k))
            if (d < dbest) then
               dbest = d; best = j
            end if
         end do
         if (used(best)) dbest = huge(1.0_dp)
         used(best) = .true.
         perm(k) = best
         worst = max(worst, dbest)
      end do
   end subroutine match

   subroutine against_pyscf()
      type(pcm_surface_t) :: surf
      type(trc_pcm_t) :: ctx
      type(error_t) :: err
      integer :: u, n, n110, k, j, model
      integer, allocatable :: perm(:)
      real(dp), allocatable :: rc(:, :), rv(:), smat(:, :), dmat(:, :), x(:), v(:), q(:), sx(:), dx(:)
      real(dp), allocatable :: sref(:, :), dref(:, :)
      real(dp) :: worst, eref, e

      call ref_require('pcm_cavity_water.ref', u)

      ! --- 302 points per sphere: the cavity itself ---
      call build_pcm_surface(zat, rat, 302, 1.2_dp, surf, err)
      if (err%has_error()) then
         print '(a)', 'check_pcm_cavity: '//err%get_message(); stop 1
      end if
      call ref_read_int(u, 'NPTS', n)
      call ref_report(nfail, n == surf%n_points, 'water, 302/sphere: point count vs pyscf', &
                      real(n - surf%n_points, dp))
      if (n /= surf%n_points) then
         print '(a,i0,a,i0)', '  terco ', surf%n_points, ' pyscf ', n; stop 1
      end if
      allocate (rc(3, n), rv(n), perm(n))
      call ref_read_matrix(u, 'COORDS', 3, n, rc)
      call match(n, surf%coords, rc, perm, worst)
      call ref_report(nfail, worst < TOL, 'every pyscf point is one of terco''s (bohr)', worst)
      call ref_read_vector(u, 'WEIGHTS', n, rv)
      call ref_report(nfail, maxval(abs(surf%weights(perm) - rv)) < TOL, 'weights', &
                      maxval(abs(surf%weights(perm) - rv)))
      call ref_read_vector(u, 'ZETA', n, rv)
      call ref_report(nfail, maxval(abs(surf%zeta(perm) - rv)) < TOL, 'exponents zeta', &
                      maxval(abs(surf%zeta(perm) - rv)))
      call ref_read_vector(u, 'SWITCH', n, rv)
      call ref_report(nfail, maxval(abs(surf%switch_f(perm) - rv)) < TOL, 'switching values', &
                      maxval(abs(surf%switch_f(perm) - rv)))
      call ref_read_vector(u, 'AREA', n, rv)
      call ref_report(nfail, maxval(abs(surf%area(perm) - rv)) < TOL, 'areas', &
                      maxval(abs(surf%area(perm) - rv)))
      call ref_read_matrix(u, 'NORMALS', 3, n, rc)
      call ref_report(nfail, maxval(abs(surf%normals(:, perm) - rc)) < TOL, 'normals', &
                      maxval(abs(surf%normals(:, perm) - rc)))

      ! --- S x and D x, all n^2 elements in aggregate ---
      allocate (smat(n, n), dmat(n, n), x(n), sx(n), dx(n))
      call pcm_sd_matrices(n, surf%coords, surf%normals, surf%zeta, surf%switch_f, surf%r_parent, smat, dmat)
      do k = 1, n
         x(k) = xfix(surf%coords(:, k))
      end do
      sx = matmul(smat, x)
      dx = matmul(dmat, x)
      call ref_read_vector(u, 'SX', n, rv)
      call ref_report(nfail, maxval(abs(sx(perm) - rv)/max(1.0_dp, abs(rv))) < TOL, &
                      'S x, relative to max(1, |S x|)', maxval(abs(sx(perm) - rv)/max(1.0_dp, abs(rv))))
      call ref_read_vector(u, 'DX', n, rv)
      call ref_report(nfail, maxval(abs(dx(perm) - rv)) < TOL, 'D x', maxval(abs(dx(perm) - rv)))
      deallocate (smat, dmat)

      ! --- the charge solve of both models for a fixed potential ---
      allocate (v(n), q(n))
      do model = PCM_MODEL_CPCM, PCM_MODEL_IEFPCM
         call ctx%init(model, EPS_WATER, err)
         call build_pcm_surface(zat, rat, 302, 1.2_dp, ctx%surface, err)
         call ctx%build_operators(err)
         if (err%has_error()) then
            print '(a)', 'check_pcm_cavity: '//err%get_message(); stop 1
         end if
         do k = 1, n
            v(k) = vfix(ctx%surface%coords(:, k))
         end do
         !$acc enter data copyin(v) create(q)
         call ctx%solve(v, q)
         !$acc update self(q)
         !$acc exit data delete(v, q)
         e = 0.5_dp*dot_product(q, v)
         if (model == PCM_MODEL_CPCM) then
            call ref_read_vector(u, 'Q_CPCM', n, rv)
            call ref_read_real(u, 'E_CPCM', eref)
            call ref_report(nfail, maxval(abs(q(perm) - rv)) < TOL, 'C-PCM charges for a fixed v', &
                            maxval(abs(q(perm) - rv)))
            call ref_report(nfail, abs(e - eref) < TOL, 'C-PCM q.v/2 for a fixed v', abs(e - eref))
         else
            call ref_read_vector(u, 'Q_IEFPCM', n, rv)
            call ref_read_real(u, 'E_IEFPCM', eref)
            call ref_report(nfail, maxval(abs(q(perm) - rv)) < TOL, 'IEF-PCM symmetrised charges for a fixed v', &
                            maxval(abs(q(perm) - rv)))
            call ref_report(nfail, abs(e - eref) < TOL, 'IEF-PCM q.v/2 for a fixed v', abs(e - eref))
         end if
         call ctx%release()
      end do
      deallocate (rc, rv, perm)

      ! --- 110 points per sphere: S and D element by element ---
      call build_pcm_surface(zat, rat, 110, 1.2_dp, surf, err)
      call ref_read_int(u, 'NPTS110', n110)
      call ref_report(nfail, n110 == surf%n_points, 'water, 110/sphere: point count vs pyscf', &
                      real(n110 - surf%n_points, dp))
      if (n110 /= surf%n_points) stop 1
      n = n110
      allocate (rc(3, n), perm(n), smat(n, n), dmat(n, n), sref(n, n), dref(n, n))
      call ref_read_matrix(u, 'COORDS110', 3, n, rc)
      call match(n, surf%coords, rc, perm, worst)
      call ref_report(nfail, worst < TOL, 'every pyscf point is one of terco''s (110)', worst)
      call pcm_sd_matrices(n, surf%coords, surf%normals, surf%zeta, surf%switch_f, surf%r_parent, smat, dmat)
      call ref_read_matrix(u, 'S110', n, n, sref)
      call ref_read_matrix(u, 'D110', n, n, dref)
      worst = 0.0_dp
      do j = 1, n
         do k = 1, n
            worst = max(worst, abs(smat(perm(k), perm(j)) - sref(k, j))/max(1.0_dp, abs(sref(k, j))))
         end do
      end do
      call ref_report(nfail, worst < TOL, 'S, every element, relative to max(1, |S|)', worst)

      worst = 0.0_dp
      do j = 1, n
         do k = 1, n
            worst = max(worst, abs(dmat(perm(k), perm(j)) - dref(k, j)))
         end do
      end do
      call ref_report(nfail, worst < TOL, 'D, every element', worst)
      worst = maxval(abs(smat - transpose(smat)))
      call ref_report(nfail, worst == 0.0_dp, 'S is symmetric', worst)
      call ref_close(u)
   end subroutine against_pyscf

end program check_pcm_cavity
