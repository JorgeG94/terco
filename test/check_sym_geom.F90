module check_sym_geom_suite
   !! trc_sym_geom: detection, orientation and symmetrization of nondegenerate
   !! tops. Water and HOCl against pyscf 2.14.0, with the molecules rotated
   !! and shifted so that the frame is not trivially the input one.
   !!
   !! The expected frames (rows x, y, z, compared up to sign, and the origin
   !! in Bohr) were printed by pyscf 2.14.0 on the geometries below with
   !! symmetry=True (HOCl) and symmetry='Cs' (water and CH2Cl2, whose
   !! symmetry=True group is C2v, which the cap reports as Cs). Static
   !! data, so no .ref file; pyscf never runs in CI. The sign of an axis is
   !! LAPACK's in pyscf and does not affect any label, so it is not compared.
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_constants, only: ANGSTROM_TO_BOHR
   use trc_error, only: error_t
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect, SYM_GROUP_ENABLED
   implicit none
   private
   public :: collect_sym_geom

   integer, parameter :: NMOL = 3
   character(len=6), parameter :: MOL_NAME(NMOL) = [character(len=6) :: "HOCl", "water", "CH2Cl2"]
   integer, parameter :: NATM(NMOL) = [3, 3, 5]
   character(len=3), parameter :: DETECTED(NMOL) = [character(len=3) :: "Cs", "C2v", "C2v"]

contains

   subroutine collect_sym_geom(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("groups detected, capped to Cs", test_groups), &
                  new_unittest("frame matches pyscf", test_pyscf_frame), &
                  new_unittest("rotation is a proper rotation", test_rotation), &
                  new_unittest("distorted molecule is C1 with the identity frame", test_c1), &
                  new_unittest("symmetrized coordinates invariant to 1e-14", test_invariance), &
                  new_unittest("perm is a bijection per operation", test_perm), &
                  new_unittest("basis identity and ghost break symmetry", test_classes), &
                  new_unittest("degenerate moments give C1", test_degenerate), &
                  new_unittest("tolerance decides detection", test_tolerance), &
                  new_unittest("bad arguments are refused", test_bad_args) &
                  ]
   end subroutine collect_sym_geom

   subroutine geometry(imol, elem, r)
      !! Geometry in Bohr (input from pyscf's run, Angstrom), and elements.
      integer, intent(in) :: imol
      integer, allocatable, intent(out) :: elem(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      select case (imol)
      case (1)  ! HOCl, planar, rotated by an arbitrary rotation and shifted
         elem = [1, 8, 17]
         r = reshape([-0.4097355013_dp, 0.3301706256_dp, 0.8801192260_dp, &
                      0.3000000000_dp, -0.2000000000_dp, 0.5000000000_dp, &
                      1.4898435896_dp, -0.1547782806_dp, 1.6978932458_dp], [3, 3])
      case (2)  ! water
         elem = [8, 1, 1]
         r = reshape([0.2543212284_dp, -0.2962864483_dp, 0.5490067312_dp, &
                      0.0334200294_dp, 0.6171335452_dp, 0.7339409511_dp, &
                      0.9320101430_dp, -0.2468419586_dp, -0.1259948004_dp], [3, 3])
      case default  ! CH2Cl2, not planar
         elem = [6, 17, 17, 1, 1]
         r = reshape([0.3000000000_dp, -0.2000000000_dp, 0.5000000000_dp, &
                      -0.9303249111_dp, -0.1525786544_dp, 1.7202671139_dp, &
                      0.7904300607_dp, -1.8070483858_dp, 0.0735333056_dp, &
                      1.1684143839_dp, 0.3327600172_dp, 0.8721870790_dp, &
                      -0.0855356395_dp, 0.2851018405_dp, -0.3902463002_dp], [3, 5])
      end select
      r = r*ANGSTROM_TO_BOHR
   end subroutine geometry

   function detect(imol, err) result(res)
      integer, intent(in) :: imol
      type(error_t), intent(out) :: err
      type(trc_sym_geom_t) :: res
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      call geometry(imol, elem, r)
      res = sym_detect(r, elem, err)
   end function detect

   subroutine test_groups(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      integer :: im
      do im = 1, NMOL
         res = detect(im, err)
         call check(error, .not. err%has_error(), "no error for "//MOL_NAME(im))
         if (allocated(error)) return
         call check(error, res%detected, DETECTED(im), "detected group of "//MOL_NAME(im))
         if (allocated(error)) return
         call check(error, res%group%name, "Cs", "reported group of "//MOL_NAME(im))
         if (allocated(error)) return
         call check(error, res%group%order, 2, "two operations")
         if (allocated(error)) return
      end do
      call check(error, count(SYM_GROUP_ENABLED), 2, "only C1 and Cs are enabled")
   end subroutine test_groups

   subroutine test_pyscf_frame(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      real(dp) :: ax(3, 3, NMOL), orig(3, NMOL)
      integer :: im, k
      ! pyscf 2.14.0 _symm_axes (rows x, y, z) and _symm_orig (Bohr)
      ax(:, :, 1) = transpose(reshape([-0.53940512_dp, 0.5709197_dp, 0.61894492_dp, &
                                       0.74658923_dp, -0.0156706_dp, 0.66510071_dp, &
                                       0.38941834_dp, 0.82085634_dp, -0.41778969_dp], [3, 3], order=[2, 1]))
      ax(:, :, 2) = transpose(reshape([0.38941834_dp, 0.82085634_dp, -0.41778969_dp, &
                                       0.59336378_dp, -0.5705068_dp, -0.56783924_dp, &
                                       -0.70446631_dp, -0.02677426_dp, -0.70923224_dp], [3, 3], order=[2, 1]))
      ax(:, :, 3) = transpose(reshape([0.38941834_dp, 0.82085634_dp, -0.41778969_dp, &
                                       0.70446631_dp, 0.02677426_dp, 0.70923224_dp, &
                                       0.59336378_dp, -0.5705068_dp, -0.56783924_dp], [3, 3], order=[2, 1]))
      orig(:, 1) = [1.98549203_dp, -0.28353597_dp, 2.45259591_dp]
      orig(:, 2) = [0.56691784_dp, -0.37794522_dp, 0.94486306_dp]
      orig(:, 3) = [0.02270669_dp, -1.52508987_dp, 1.5287231_dp]
      do im = 1, NMOL
         res = detect(im, err)
         call check(error, maxval(abs(res%shift - orig(:, im))) < 1.0e-6_dp, "origin of "//MOL_NAME(im))
         if (allocated(error)) return
         do k = 1, 3
            ! rows are unit vectors, so parallel up to sign means |dot| = 1
            call check(error, abs(abs(dot_product(res%rot(k, :), ax(:, k, im))) - 1.0_dp) < 1.0e-7_dp, &
                       "axis of "//MOL_NAME(im))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_pyscf_frame

   subroutine test_rotation(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      integer :: im, i
      real(dp) :: m(3, 3), det
      do im = 1, NMOL
         res = detect(im, err)
         m = matmul(res%rot, transpose(res%rot))
         do i = 1, 3
            m(i, i) = m(i, i) - 1.0_dp
         end do
         call check(error, maxval(abs(m)) < 1.0e-14_dp, "orthogonal "//MOL_NAME(im))
         if (allocated(error)) return
         det = res%rot(1, 1)*(res%rot(2, 2)*res%rot(3, 3) - res%rot(2, 3)*res%rot(3, 2)) &
               - res%rot(1, 2)*(res%rot(2, 1)*res%rot(3, 3) - res%rot(2, 3)*res%rot(3, 1)) &
               + res%rot(1, 3)*(res%rot(2, 1)*res%rot(3, 2) - res%rot(2, 2)*res%rot(3, 1))
         call check(error, abs(det - 1.0_dp) < 1.0e-14_dp, "determinant +1 "//MOL_NAME(im))
         if (allocated(error)) return
      end do
   end subroutine test_rotation

   subroutine test_c1(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      call geometry(3, elem, r)  ! (three atoms are always planar, hence Cs)
      r(:, 4) = r(:, 4) + [0.01_dp, -0.02_dp, 0.015_dp]  ! no plane survives
      r(:, 2) = r(:, 2) + [-0.02_dp, 0.01_dp, 0.03_dp]
      res = sym_detect(r, elem, err)
      call check(error, res%detected, "C1", "detected")
      if (allocated(error)) return
      call check(error, res%group%name, "C1", "reported")
      if (allocated(error)) return
      call check(error, all(res%rot == reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])), "identity frame")
      if (allocated(error)) return
      call check(error, all(res%shift == 0.0_dp), "no shift")
      if (allocated(error)) return
      call check(error, all(res%r_sym == r) .and. res%residual == 0.0_dp, "coordinates untouched")
   end subroutine test_c1

   subroutine noisy(imol, elem, r)
      integer, intent(in) :: imol
      integer, allocatable, intent(out) :: elem(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      integer :: i, k
      call geometry(imol, elem, r)
      do i = 1, size(r, 2)
         do k = 1, 3
            r(k, i) = r(k, i) + 1.0e-7_dp*sin(7.0_dp*i + 3.0_dp*k)
         end do
      end do
   end subroutine noisy

   subroutine test_invariance(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res, again
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: y(3), worst
      integer :: im, i, k, d
      do im = 1, NMOL
         call noisy(im, elem, r)
         res = sym_detect(r, elem, err)
         call check(error, res%group%name, "Cs", "still Cs with 1e-7 noise, "//MOL_NAME(im))
         if (allocated(error)) return
         call check(error, res%residual > 0.0_dp .and. res%residual < 5.0e-7_dp, "residual is the noise")
         if (allocated(error)) return
         worst = 0.0_dp
         do k = 1, res%group%order
            do i = 1, size(r, 2)
               y = res%r_std(:, i)
               do d = 1, 3
                  if (btest(res%group%ops(k), d - 1)) y(d) = -y(d)
               end do
               worst = max(worst, maxval(abs(y - res%r_std(:, res%perm(i, k)))))
            end do
         end do
         call check(error, worst < 1.0e-14_dp, "g x = x[perm] for "//MOL_NAME(im))
         if (allocated(error)) return
         ! the input-frame coordinates are the standard-frame ones moved back
         again = sym_detect(res%r_sym, elem, err)
         call check(error, again%residual < 1.0e-14_dp, "symmetrized geometry is a fixed point")
         if (allocated(error)) return
      end do
   end subroutine test_invariance

   subroutine test_perm(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      integer :: im, k, i
      logical, allocatable :: seen(:)
      do im = 1, NMOL
         res = detect(im, err)
         call check(error, size(res%perm, 2), res%group%order, "one column per operation")
         if (allocated(error)) return
         call check(error, all(res%perm(:, 1) == [(i, i=1, NATM(im))]), "E is the identity")
         if (allocated(error)) return
         do k = 1, res%group%order
            allocate (seen(NATM(im)))
            seen = .false.
            do i = 1, NATM(im)
               seen(res%perm(i, k)) = .true.
            end do
            call check(error, all(seen), "bijection of "//MOL_NAME(im))
            deallocate (seen)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_perm

   subroutine test_classes(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      call geometry(3, elem, r)
      res = sym_detect(r, elem, err, basis_id=[0, 0, 0, 0, 1])
      call check(error, res%detected, "Cs", "a different basis on one H leaves one plane")
      if (allocated(error)) return
      res = sym_detect(r, elem, err, ghost=[.false., .false., .false., .false., .true.])
      call check(error, res%detected, "Cs", "a ghost H leaves one plane")
      if (allocated(error)) return
      res = sym_detect(r, elem, err, basis_id=[0, 0, 0, 1, 1])
      call check(error, res%detected, "C2v", "equal labels keep C2v")
      if (allocated(error)) return
      res = sym_detect(r, elem, err, basis_id=[0, 0, 1, 0, 0])
      call check(error, res%detected, "Cs", "a different basis on one Cl leaves one plane")
   end subroutine test_classes

   subroutine test_degenerate(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      ! CO2 is linear: two vanishing moments
      res = sym_detect(reshape([0.0_dp, 0.0_dp, -2.2_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.2_dp], &
                               [3, 3]), [8, 6, 8], err)
      call check(error, res%degenerate, "flagged")
      if (allocated(error)) return
      call check(error, res%group%name, "C1", "C1 until S12")
      if (allocated(error)) return
      ! a single atom
      res = sym_detect(reshape([1.0_dp, 2.0_dp, 3.0_dp], [3, 1]), [8], err)
      call check(error, res%degenerate .and. res%group%name == "C1", "single atom")
   end subroutine test_degenerate

   subroutine test_tolerance(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      call noisy(3, elem, r)
      res = sym_detect(r, elem, err, tol=1.0e-9_dp)
      call check(error, res%group%name, "C1", "1e-7 noise is not symmetric at 1e-9")
      if (allocated(error)) return
      res = sym_detect(r, elem, err, tol=1.0e-6_dp)
      call check(error, res%group%name, "Cs", "but is at 1e-6")
   end subroutine test_tolerance

   subroutine test_bad_args(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      res = sym_detect(reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), [1, 1], err)
      call check(error, err%has_error(), "size mismatch")
      if (allocated(error)) return
      res = sym_detect(reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), [1], err, tol=0.0_dp)
      call check(error, err%has_error(), "zero tolerance")
   end subroutine test_bad_args

end module check_sym_geom_suite

program check_sym_geom
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_sym_geom_suite, only: collect_sym_geom
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_sym_geom, error_unit, stat)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_geom
