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
! The reference is read once, sequentially (pyscf_ref_io.F90's own
! convention), by load_reference below, into g_ref -- before run_testsuite
! runs any case -- so each test-drive case here just builds its own
! surface/operators and compares against the already-read arrays.
!
module check_pcm_cavity_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_pcm, only: trc_pcm_t, pcm_surface_t, pcm_build_surface, pcm_sd_matrices, &
                      pcm_switch_h, pcm_model_from_name, PCM_MODEL_CPCM, PCM_MODEL_IEFPCM
   use trc_pcm_radii, only: pcm_radii_cavity_radius
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, &
                           ref_read_vector, ref_read_matrix
   implicit none
   private
   public :: collect_pcm_cavity_tests, load_reference

   real(dp), parameter :: PI = 3.14159265358979323846_dp
   real(dp), parameter :: EPS_WATER = 78.3553_dp
   real(dp), parameter :: TOL = 1.0e-12_dp
   !! Absolute, on everything but S. S's diagonal is zeta sqrt(2/pi)/F, and
   !! a point nearly switched off has F ~ 1e-6 and S_ii ~ 1e5, where the
   !! last bits of F's product of switches (1e-15 apart between the two
   !! codes) are 1e-8 absolute. So S, and S x which its diagonal dominates,
   !! are compared relative to max(1, |S|): the same 1e-12, scaled.
   real(dp), parameter :: ORIGIN(3, 1) = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])

   type :: pcm_ref_t
      integer :: npts = 0, npts110 = 0
      real(dp), allocatable :: coords(:, :), weights(:), zeta(:), switch(:), area(:), normals(:, :)
      real(dp), allocatable :: sx(:), dx(:)
      real(dp), allocatable :: q_cpcm(:), q_iefpcm(:)
      real(dp) :: e_cpcm = 0.0_dp, e_iefpcm = 0.0_dp
      real(dp), allocatable :: coords110(:, :), s110(:, :), d110(:, :)
   end type pcm_ref_t

   type(pcm_ref_t) :: g_ref
   integer :: g_natm
   integer, allocatable :: g_zat(:)
   real(dp), allocatable :: g_rat(:, :)

contains

   subroutine load_reference()
      integer :: u, n, n110
      call read_xyz('water.xyz', g_natm, g_zat, g_rat)
      call ref_require('pcm_cavity_water.ref', u)
      call ref_read_int(u, 'NPTS', n)
      g_ref%npts = n
      allocate (g_ref%coords(3, n), g_ref%weights(n), g_ref%zeta(n), g_ref%switch(n), g_ref%area(n), &
                g_ref%normals(3, n), g_ref%sx(n), g_ref%dx(n), g_ref%q_cpcm(n), g_ref%q_iefpcm(n))
      call ref_read_matrix(u, 'COORDS', 3, n, g_ref%coords)
      call ref_read_vector(u, 'WEIGHTS', n, g_ref%weights)
      call ref_read_vector(u, 'ZETA', n, g_ref%zeta)
      call ref_read_vector(u, 'SWITCH', n, g_ref%switch)
      call ref_read_vector(u, 'AREA', n, g_ref%area)
      call ref_read_matrix(u, 'NORMALS', 3, n, g_ref%normals)
      call ref_read_vector(u, 'SX', n, g_ref%sx)
      call ref_read_vector(u, 'DX', n, g_ref%dx)
      call ref_read_vector(u, 'Q_CPCM', n, g_ref%q_cpcm)
      call ref_read_real(u, 'E_CPCM', g_ref%e_cpcm)
      call ref_read_vector(u, 'Q_IEFPCM', n, g_ref%q_iefpcm)
      call ref_read_real(u, 'E_IEFPCM', g_ref%e_iefpcm)
      call ref_read_int(u, 'NPTS110', n110)
      g_ref%npts110 = n110
      allocate (g_ref%coords110(3, n110), g_ref%s110(n110, n110), g_ref%d110(n110, n110))
      call ref_read_matrix(u, 'COORDS110', 3, n110, g_ref%coords110)
      call ref_read_matrix(u, 'S110', n110, n110, g_ref%s110)
      call ref_read_matrix(u, 'D110', n110, n110, g_ref%d110)
      call ref_close(u)
   end subroutine load_reference

   subroutine collect_pcm_cavity_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("switch_h: 0 and 1 at the ends, odd about 1/2", test_switch_h), &
                  new_unittest("isolated sphere keeps every point, area = 4 pi R^2", test_isolated_sphere), &
                  new_unittest("overlapping spheres bury points between them", test_overlapping_spheres), &
                  new_unittest("Born ion, C-PCM: -f/(2R) and total charge -f", test_born_cpcm), &
                  new_unittest("Born ion, IEF-PCM: the exact -(eps-1)/(2 R eps)", test_born_iefpcm), &
                  new_unittest("refusals: model, eps=1, Lebedev order, missing radius", test_refusals), &
                  new_unittest("water/302-per-sphere cavity vs pyscf", test_cavity_302), &
                  new_unittest("S x and D x for a fixed x vs pyscf", test_sx_dx), &
                  new_unittest("C-PCM/IEF-PCM charge solve for a fixed potential vs pyscf", test_charge_solve), &
                  new_unittest("water/110-per-sphere S and D, element by element, vs pyscf", test_cavity_110) &
                  ]
   end subroutine collect_pcm_cavity_tests

   subroutine test_switch_h(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: h
      h = max(abs(pcm_switch_h(0.0_dp)), abs(pcm_switch_h(1.0_dp) - 1.0_dp), &
              abs(pcm_switch_h(0.5_dp) - 0.5_dp), abs(pcm_switch_h(0.3_dp) + pcm_switch_h(0.7_dp) - 1.0_dp), &
              abs(pcm_switch_h(-0.2_dp)), abs(pcm_switch_h(1.4_dp) - 1.0_dp))
      call check(error, h < 1.0e-15_dp, 'switch_h shape residual = '//to_string(h, "es10.2"))
   end subroutine test_switch_h

   subroutine test_isolated_sphere(error)
      !! An isolated sphere keeps every point; its weights sum to one and
      !! nothing is switched, so the area is the sphere's to roundoff.
      type(error_type), allocatable, intent(out) :: error
      type(pcm_surface_t) :: surf
      type(error_t) :: err
      real(dp) :: radius, e
      call pcm_build_surface([8], ORIGIN, 110, 1.2_dp, surf, err)
      call pcm_radii_cavity_radius(8, 1.2_dp, radius, err)
      call check(error,.not. err%has_error() .and. surf%n_points == 110, &
                 'isolated sphere must keep all 110 points, got '//to_string(surf%n_points))
      if (allocated(error)) return
      e = abs(sum(surf%area) - 4.0_dp*PI*radius**2)
      call check(error, e < 1.0e-10_dp, 'isolated sphere area vs 4 pi R^2, |diff| = '//to_string(e, "es10.2"))
   end subroutine test_isolated_sphere

   subroutine test_overlapping_spheres(error)
      !! Two oxygens a bond apart lose the points between them.
      type(error_type), allocatable, intent(out) :: error
      type(pcm_surface_t) :: surf
      type(error_t) :: err
      real(dp) :: radius
      call pcm_build_surface([8, 8], reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.3_dp], [3, 2]), &
                             110, 1.2_dp, surf, err)
      call pcm_radii_cavity_radius(8, 1.2_dp, radius, err)
      call check(error,.not. err%has_error() .and. surf%n_points < 220 .and. &
                 surf%n_points + surf%n_discarded == 220 .and. &
                 sum(surf%area) < 8.0_dp*PI*radius**2, &
                 'overlapping spheres must bury points, kept + dropped = 220, got kept='// &
                 to_string(surf%n_points)//' dropped='//to_string(surf%n_discarded))
   end subroutine test_overlapping_spheres

   subroutine test_born_cpcm(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e, qt, f, radius
      type(error_t) :: err
      call pcm_radii_cavity_radius(8, 1.2_dp, radius, err)
      f = (EPS_WATER - 1.0_dp)/EPS_WATER
      call born(PCM_MODEL_CPCM, e, qt)
      call check(error, abs(e + 0.5_dp*f/radius) < 1.0e-5_dp*0.5_dp*f/radius, &
                 'Born ion, C-PCM: -f/(2R), relative residual = '// &
                 to_string(abs(e + 0.5_dp*f/radius)/(0.5_dp*f/radius), "es10.2"))
      if (allocated(error)) return
      call check(error, abs(qt + f) < 1.0e-4_dp, &
                 'Born ion, C-PCM: Gauss, total charge must be -f, |diff| = '//to_string(abs(qt + f), "es10.2"))
   end subroutine test_born_cpcm

   subroutine test_born_iefpcm(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e, qt, f, radius
      type(error_t) :: err
      call pcm_radii_cavity_radius(8, 1.2_dp, radius, err)
      f = (EPS_WATER - 1.0_dp)/EPS_WATER
      call born(PCM_MODEL_IEFPCM, e, qt)
      call check(error, abs(e + 0.5_dp*f/radius) < 1.0e-4_dp*0.5_dp*f/radius, &
                 'Born ion, IEF-PCM: the exact -(eps-1)/(2 R eps), relative residual = '// &
                 to_string(abs(e + 0.5_dp*f/radius)/(0.5_dp*f/radius), "es10.2"))
   end subroutine test_born_iefpcm

   subroutine born(model, e, qt)
      integer, intent(in) :: model
      real(dp), intent(out) :: e, qt
      type(trc_pcm_t) :: ctx
      type(trc_linalg_t) :: la
      type(error_t) :: err
      real(dp), allocatable :: v(:), q(:)
      integer :: k
      call la%init(1)
      call ctx%init(model, EPS_WATER, err)
      call pcm_build_surface([8], ORIGIN, 302, 1.2_dp, ctx%surface, err)
      call ctx%build_operators(la, err)
      if (err%has_error()) then
         print '(a)', 'check_pcm_cavity: '//err%get_message(); stop 1
      end if
      allocate (v(ctx%surface%n_points), q(ctx%surface%n_points))
      do k = 1, ctx%surface%n_points
         v(k) = erf(ctx%surface%zeta(k)*norm2(ctx%surface%coords(:, k)))/norm2(ctx%surface%coords(:, k))
      end do
      ! solve's vectors live where K does: on the device under the GPU build
      !$acc enter data copyin(v) create(q)
      call ctx%solve(la, v, q)
      !$acc update self(q)
      !$acc exit data delete(v, q)
      e = 0.5_dp*dot_product(q, v)
      qt = sum(q)
      call ctx%release()
      call la%release()
   end subroutine born

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_pcm_t) :: ctx
      type(pcm_surface_t) :: surf
      type(error_t) :: e1, e2, e3, e4, e5
      call ctx%init(pcm_model_from_name('cosmo'), EPS_WATER, e1)
      call ctx%init(PCM_MODEL_CPCM, 1.0_dp, e2)
      call pcm_build_surface([8], ORIGIN, 74, 1.2_dp, surf, e3)
      call pcm_build_surface([26], ORIGIN, 302, 1.2_dp, surf, e4)
      call pcm_build_surface([8], ORIGIN, 303, 1.2_dp, surf, e5)
      call check(error, e1%has_error(), "refuses 'cosmo' (a model this does not solve)")
      if (allocated(error)) return
      call check(error, e2%has_error(), 'refuses eps = 1')
      if (allocated(error)) return
      call check(error, e3%has_error(), 'refuses Lebedev 74 (no fitted exponent)')
      if (allocated(error)) return
      call check(error, e4%has_error(), 'refuses iron (no cavity radius)')
      if (allocated(error)) return
      call check(error, e5%has_error(), 'refuses 303 points (no such order)')
      if (allocated(error)) return
      call check(error, pcm_model_from_name('IEF-PCM') == PCM_MODEL_IEFPCM .and. &
                 pcm_model_from_name('cpcm') == PCM_MODEL_CPCM, 'model names')
   end subroutine test_refusals

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

   subroutine test_cavity_302(error)
      type(error_type), allocatable, intent(out) :: error
      type(pcm_surface_t) :: surf
      type(error_t) :: err
      integer :: n
      integer, allocatable :: perm(:)
      real(dp) :: worst

      n = g_ref%npts
      call pcm_build_surface(g_zat, g_rat, 302, 1.2_dp, surf, err)
      call check(error,.not. err%has_error(), 'check_pcm_cavity: '//err%get_message())
      if (allocated(error)) return
      call check(error, n == surf%n_points, &
                 'water, 302/sphere: point count vs pyscf, terco='//to_string(surf%n_points)// &
                 ' pyscf='//to_string(n))
      if (allocated(error)) return

      allocate (perm(n))
      call match(n, surf%coords, g_ref%coords, perm, worst)
      call check(error, worst < TOL, "every pyscf point is one of terco's (bohr), worst = "// &
                 to_string(worst, "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(abs(surf%weights(perm) - g_ref%weights)) < TOL, &
                 'weights, worst |diff| = '//to_string(maxval(abs(surf%weights(perm) - g_ref%weights)), "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(abs(surf%zeta(perm) - g_ref%zeta)) < TOL, &
                 'exponents zeta, worst |diff| = '//to_string(maxval(abs(surf%zeta(perm) - g_ref%zeta)), "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(abs(surf%switch_f(perm) - g_ref%switch)) < TOL, &
                 'switching values, worst |diff| = '// &
                 to_string(maxval(abs(surf%switch_f(perm) - g_ref%switch)), "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(abs(surf%area(perm) - g_ref%area)) < TOL, &
                 'areas, worst |diff| = '//to_string(maxval(abs(surf%area(perm) - g_ref%area)), "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(abs(surf%normals(:, perm) - g_ref%normals)) < TOL, &
                 'normals, worst |diff| = '//to_string(maxval(abs(surf%normals(:, perm) - g_ref%normals)), "es10.2"))
   end subroutine test_cavity_302

   subroutine test_sx_dx(error)
      !! S x and D x, all n^2 elements in aggregate, at 302 points/sphere.
      type(error_type), allocatable, intent(out) :: error
      type(pcm_surface_t) :: surf
      type(error_t) :: err
      integer :: n, k
      integer, allocatable :: perm(:)
      real(dp), allocatable :: smat(:, :), dmat(:, :), x(:), sx(:), dx(:)
      real(dp) :: worst

      n = g_ref%npts
      call pcm_build_surface(g_zat, g_rat, 302, 1.2_dp, surf, err)
      ! Guard the count before indexing: match and pcm_sd_matrices take
      ! surf%coords as (3, n) with pyscf's n.
      call check(error, .not. err%has_error(), 'check_pcm_cavity: '//err%get_message())
      if (allocated(error)) return
      call check(error, n == surf%n_points, &
                 'water, 302/sphere: point count vs pyscf, terco='//to_string(surf%n_points)// &
                 ' pyscf='//to_string(n))
      if (allocated(error)) return
      allocate (perm(n))
      call match(n, surf%coords, g_ref%coords, perm, worst)
      call check(error, worst < TOL, "point matching (bohr), worst = "//to_string(worst, "es10.2"))
      if (allocated(error)) return

      allocate (smat(n, n), dmat(n, n), x(n), sx(n), dx(n))
      call pcm_sd_matrices(n, surf%coords, surf%normals, surf%zeta, surf%switch_f, surf%r_parent, smat, dmat)
      do k = 1, n
         x(k) = xfix(surf%coords(:, k))
      end do
      sx = matmul(smat, x)
      dx = matmul(dmat, x)
      call check(error, maxval(abs(sx(perm) - g_ref%sx)/max(1.0_dp, abs(g_ref%sx))) < TOL, &
                 'S x, relative to max(1, |S x|), worst = '// &
                 to_string(maxval(abs(sx(perm) - g_ref%sx)/max(1.0_dp, abs(g_ref%sx))), "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(abs(dx(perm) - g_ref%dx)) < TOL, &
                 'D x, worst |diff| = '//to_string(maxval(abs(dx(perm) - g_ref%dx)), "es10.2"))
   end subroutine test_sx_dx

   subroutine test_charge_solve(error)
      !! The charge solve of both models for a fixed potential, at 302
      !! points/sphere.
      type(error_type), allocatable, intent(out) :: error
      type(pcm_surface_t) :: surf
      type(trc_pcm_t) :: ctx
      type(trc_linalg_t) :: la
      type(error_t) :: err
      integer :: n, k, model
      integer, allocatable :: perm(:)
      real(dp), allocatable :: v(:), q(:)
      real(dp) :: worst, e

      n = g_ref%npts
      call pcm_build_surface(g_zat, g_rat, 302, 1.2_dp, surf, err)
      ! Guard the count before indexing: match and pcm_sd_matrices take
      ! surf%coords as (3, n) with pyscf's n.
      call check(error, .not. err%has_error(), 'check_pcm_cavity: '//err%get_message())
      if (allocated(error)) return
      call check(error, n == surf%n_points, &
                 'water, 302/sphere: point count vs pyscf, terco='//to_string(surf%n_points)// &
                 ' pyscf='//to_string(n))
      if (allocated(error)) return
      allocate (perm(n))
      call match(n, surf%coords, g_ref%coords, perm, worst)
      call check(error, worst < TOL, "point matching (bohr), worst = "//to_string(worst, "es10.2"))
      if (allocated(error)) return

      allocate (v(n), q(n))
      call la%init(1)
      do model = PCM_MODEL_CPCM, PCM_MODEL_IEFPCM
         call ctx%init(model, EPS_WATER, err)
         call pcm_build_surface(g_zat, g_rat, 302, 1.2_dp, ctx%surface, err)
         call ctx%build_operators(la, err)
         if (err%has_error()) then
            print '(a)', 'check_pcm_cavity: '//err%get_message(); stop 1
         end if
         do k = 1, n
            v(k) = vfix(ctx%surface%coords(:, k))
         end do
         !$acc enter data copyin(v) create(q)
         call ctx%solve(la, v, q)
         !$acc update self(q)
         !$acc exit data delete(v, q)
         e = 0.5_dp*dot_product(q, v)
         if (model == PCM_MODEL_CPCM) then
            call check(error, maxval(abs(q(perm) - g_ref%q_cpcm)) < TOL, &
                       'C-PCM charges for a fixed v, worst |diff| = '// &
                       to_string(maxval(abs(q(perm) - g_ref%q_cpcm)), "es10.2"))
            if (.not. allocated(error)) &
               call check(error, abs(e - g_ref%e_cpcm) < TOL, &
                          'C-PCM q.v/2 for a fixed v, |diff| = '//to_string(abs(e - g_ref%e_cpcm), "es10.2"))
         else
            call check(error, maxval(abs(q(perm) - g_ref%q_iefpcm)) < TOL, &
                       'IEF-PCM symmetrised charges for a fixed v, worst |diff| = '// &
                       to_string(maxval(abs(q(perm) - g_ref%q_iefpcm)), "es10.2"))
            if (.not. allocated(error)) &
               call check(error, abs(e - g_ref%e_iefpcm) < TOL, &
                          'IEF-PCM q.v/2 for a fixed v, |diff| = '//to_string(abs(e - g_ref%e_iefpcm), "es10.2"))
         end if
         call ctx%release()
         if (allocated(error)) exit
      end do
      call la%release()
   end subroutine test_charge_solve

   subroutine test_cavity_110(error)
      !! S and D element by element, at 110 points/sphere.
      type(error_type), allocatable, intent(out) :: error
      type(pcm_surface_t) :: surf
      type(error_t) :: err
      integer :: n, j, k
      integer, allocatable :: perm(:)
      real(dp), allocatable :: smat(:, :), dmat(:, :)
      real(dp) :: worst

      n = g_ref%npts110
      call pcm_build_surface(g_zat, g_rat, 110, 1.2_dp, surf, err)
      call check(error, n == surf%n_points, &
                 'water, 110/sphere: point count vs pyscf, terco='//to_string(surf%n_points)// &
                 ' pyscf='//to_string(n))
      if (allocated(error)) return

      allocate (perm(n), smat(n, n), dmat(n, n))
      call match(n, surf%coords, g_ref%coords110, perm, worst)
      call check(error, worst < TOL, "every pyscf point is one of terco's (110), worst = "// &
                 to_string(worst, "es10.2"))
      if (allocated(error)) return

      call pcm_sd_matrices(n, surf%coords, surf%normals, surf%zeta, surf%switch_f, surf%r_parent, smat, dmat)
      worst = 0.0_dp
      do j = 1, n
         do k = 1, n
            worst = max(worst, abs(smat(perm(k), perm(j)) - g_ref%s110(k, j))/max(1.0_dp, abs(g_ref%s110(k, j))))
         end do
      end do
      call check(error, worst < TOL, 'S, every element, relative to max(1, |S|), worst = '// &
                 to_string(worst, "es10.2"))
      if (allocated(error)) return

      worst = 0.0_dp
      do j = 1, n
         do k = 1, n
            worst = max(worst, abs(dmat(perm(k), perm(j)) - g_ref%d110(k, j)))
         end do
      end do
      call check(error, worst < TOL, 'D, every element, worst |diff| = '//to_string(worst, "es10.2"))
      if (allocated(error)) return

      worst = maxval(abs(smat - transpose(smat)))
      call check(error, worst == 0.0_dp, 'S must be exactly symmetric, worst |S - S^T| = '// &
                 to_string(worst, "es10.2"))
   end subroutine test_cavity_110

end module check_pcm_cavity_suite

program check_pcm_cavity
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_pcm_cavity_suite, only: collect_pcm_cavity_tests, load_reference
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat

   call load_reference()

   stat = 0
   testsuites = [new_testsuite("check_pcm_cavity", collect_pcm_cavity_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_pcm_cavity
