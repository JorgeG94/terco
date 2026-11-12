module check_sym_geom_suite
   !! trc_sym_geom: detection, orientation and symmetrization. HOCl, water and
   !! CH2Cl2 against pyscf 2.14.0, with the molecules rotated and shifted so
   !! that the frame is not trivially the input one (test/sym_*_rot.xyz).
   !!
   !! The expected frames (rows x, y, z, compared up to sign, and the origin
   !! in Bohr) are pyscf's `_symm_axes` and `_symm_orig`, in
   !! pyscf_refs/sym_geom_frames.ref (generate.py `sym_geom_frames`). The sign
   !! of an axis is LAPACK's in pyscf and does not affect any label, so it is
   !! not compared.
   !!
   !! Tolerance edge cases: a requested subgroup on input axes tilted below
   !! pyscf's axis test, symmetric and spherical tops with noise below `tol`,
   !! a nearly linear molecule bent in different directions, and a `tol` so
   !! large that two atoms share a nearest image.
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect, SYM_GROUP_ENABLED
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_vector, ref_read_matrix
   implicit none
   private
   public :: collect_sym_geom

   integer, parameter :: NMOL = 3
   character(len=6), parameter :: MOL_NAME(NMOL) = [character(len=6) :: "HOCl", "water", "CH2Cl2"]
   character(len=14), parameter :: MOL_FILE(NMOL) = [character(len=14) :: "sym_hocl_rot", "sym_water_rot", &
                                                     "sym_ch2cl2_rot"]
   integer, parameter :: NATM(NMOL) = [3, 3, 5]
   character(len=3), parameter :: DETECTED(NMOL) = [character(len=3) :: "Cs", "C2v", "C2v"]
   integer, parameter :: ORDER(NMOL) = [2, 4, 4]  !! operations of the reported group

contains

   subroutine collect_sym_geom(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("groups detected and reported", test_groups), &
                  new_unittest("frame matches pyscf", test_pyscf_frame), &
                  new_unittest("rotation is a proper rotation", test_rotation), &
                  new_unittest("distorted molecule is C1 with the identity frame", test_c1), &
                  new_unittest("symmetrized coordinates invariant to 1e-14", test_invariance), &
                  new_unittest("perm is a bijection per operation", test_perm), &
                  new_unittest("basis identity and ghost break symmetry", test_classes), &
                  new_unittest("degenerate moments: linear molecule and atom in D2h/C2v", test_degenerate), &
                  new_unittest("tolerance decides detection", test_tolerance), &
                  new_unittest("requested C2 on water tilted by 1e-4 rad: pyscf's frame", test_requested_tilted), &
                  new_unittest("tops with noise below tol keep their group", test_noisy_tops), &
                  new_unittest("a bent near-linear molecule: same group in every direction", test_bent_linear), &
                  new_unittest("a large tol never gives a perm that is not a permutation", test_large_tol), &
                  new_unittest("bad arguments are refused", test_bad_args) &
                  ]
   end subroutine collect_sym_geom

   subroutine geometry(imol, elem, r)
      !! Geometry in Bohr (test/sym_*_rot.xyz, Angstrom, the input of pyscf's
      !! run), and elements.
      integer, intent(in) :: imol
      integer, allocatable, intent(out) :: elem(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      integer :: nat
      call read_xyz(trim(MOL_FILE(imol))//".xyz", nat, elem, r)
   end subroutine geometry

   function detect(imol, err) result(res)
      integer, intent(in) :: imol
      type(error_t), intent(out) :: err
      type(trc_sym_geom_t) :: res
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      call geometry(imol, elem, r)
      res = sym_geom_detect(r, elem, err)
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
         call check(error, res%group%name, DETECTED(im), "reported group of "//MOL_NAME(im))
         if (allocated(error)) return
         call check(error, res%group%order, ORDER(im), "operations of "//MOL_NAME(im))
         if (allocated(error)) return
      end do
      ! C1 and Cs stay enabled whatever a later stage adds (the order of the
      ! flags is that of sym_group_names: D2h, C2h, C2v, D2, Cs, Ci, C2, C1)
      call check(error, SYM_GROUP_ENABLED(5) .and. SYM_GROUP_ENABLED(8), "C1 and Cs are enabled")
   end subroutine test_groups

   subroutine pyscf_frames(ax, orig)
      !! pyscf 2.14.0 _symm_axes (rows x, y, z) and _symm_orig (Bohr) of the
      !! NMOL molecules, then of the tilted water, in the order of the file.
      real(dp), intent(out) :: ax(3, 3, NMOL + 1), orig(3, NMOL + 1)
      character(len=14), parameter :: TILT = "sym_water_tilt"
      integer :: im, u
      call ref_require("sym_geom_frames.ref", u)
      do im = 1, NMOL + 1
         call ref_read_matrix(u, "AXES_"//trim(merge(MOL_FILE(min(im, NMOL)), TILT, im <= NMOL)), 3, 3, ax(:, :, im))
         call ref_read_vector(u, "ORIG_"//trim(merge(MOL_FILE(min(im, NMOL)), TILT, im <= NMOL)), 3, orig(:, im))
      end do
      call ref_close(u)
   end subroutine pyscf_frames

   subroutine test_pyscf_frame(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      real(dp) :: ax(3, 3, NMOL + 1), orig(3, NMOL + 1)
      integer :: im
      call pyscf_frames(ax, orig)
      do im = 1, NMOL
         res = detect(im, err)
         call check_frame(error, res, ax(:, :, im), orig(:, im), MOL_NAME(im))
         if (allocated(error)) return
      end do
   end subroutine test_pyscf_frame

   subroutine check_frame(error, res, ax, orig, label)
      !! `res`'s frame is pyscf's: the same origin, each axis parallel up to sign.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t), intent(in) :: res
      real(dp), intent(in) :: ax(3, 3), orig(3)
      character(len=*), intent(in) :: label
      integer :: k
      call check(error, maxval(abs(res%shift - orig)) < 1.0e-6_dp, "origin of "//label)
      if (allocated(error)) return
      do k = 1, 3
         ! rows are unit vectors, so parallel up to sign means |dot| = 1
         call check(error, abs(abs(dot_product(res%rot(k, :), ax(k, :))) - 1.0_dp) < 1.0e-7_dp, "axis of "//label)
         if (allocated(error)) return
      end do
   end subroutine check_frame

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
      res = sym_geom_detect(r, elem, err)
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
         res = sym_geom_detect(r, elem, err)
         call check(error, res%group%name, DETECTED(im), "still the same group with 1e-7 noise, "//MOL_NAME(im))
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
         again = sym_geom_detect(res%r_sym, elem, err)
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
      res = sym_geom_detect(r, elem, err, basis_id=[0, 0, 0, 0, 1])
      call check(error, res%detected, "Cs", "a different basis on one H leaves one plane")
      if (allocated(error)) return
      res = sym_geom_detect(r, elem, err, ghost=[.false., .false., .false., .false., .true.])
      call check(error, res%detected, "Cs", "a ghost H leaves one plane")
      if (allocated(error)) return
      res = sym_geom_detect(r, elem, err, basis_id=[0, 0, 0, 1, 1])
      call check(error, res%detected, "C2v", "equal labels keep C2v")
      if (allocated(error)) return
      res = sym_geom_detect(r, elem, err, basis_id=[0, 0, 1, 0, 0])
      call check(error, res%detected, "Cs", "a different basis on one Cl leaves one plane")
   end subroutine test_classes

   subroutine test_degenerate(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      ! CO2 is linear: two vanishing moments; the group is D2h (inversion centre)
      res = sym_geom_detect(reshape([0.0_dp, 0.0_dp, -2.2_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.2_dp], &
                               [3, 3]), [8, 6, 8], err)
      call check(error, res%degenerate, "flagged")
      if (allocated(error)) return
      call check(error, res%group%name, "D2h", "linear with an inversion centre: D2h")
      if (allocated(error)) return
      ! the same without the centre (OCS-like): C2v
      res = sym_geom_detect(reshape([0.0_dp, 0.0_dp, -2.2_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.9_dp], &
                               [3, 3]), [8, 6, 16], err)
      call check(error, res%group%name, "C2v", "linear without an inversion centre: C2v")
      if (allocated(error)) return
      ! a single atom: D2h, frame = input axes, origin = the atom
      res = sym_geom_detect(reshape([1.0_dp, 2.0_dp, 3.0_dp], [3, 1]), [8], err)
      call check(error, res%degenerate .and. res%group%name == "D2h", "single atom")
      if (allocated(error)) return
      call check(error, maxval(abs(res%shift - [1.0_dp, 2.0_dp, 3.0_dp])) < 1.0e-14_dp, "atom: origin on the atom")
      if (allocated(error)) return
      call check(error, maxval(abs(res%rot - reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]))) < 1.0e-14_dp, "atom: input axes")
   end subroutine test_degenerate

   subroutine test_tolerance(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      call noisy(3, elem, r)
      res = sym_geom_detect(r, elem, err, tol=1.0e-9_dp)
      call check(error, res%group%name, "C1", "1e-7 noise is not symmetric at 1e-9")
      if (allocated(error)) return
      res = sym_geom_detect(r, elem, err, tol=1.0e-6_dp)
      call check(error, res%group%name, "C2v", "but is at 1e-6")
   end subroutine test_tolerance

   subroutine test_requested_tilted(error)
      !! Water turned by 1e-4 rad off the yz plane, with C2 requested. Its
      !! principal axes pass the axis test of pyscf (and of `input_frame_holds`),
      !! but its mirrors do not hold in the input frame, so pyscf keeps the C2v
      !! frame instead of the input axes. Keeping the input axes on the axis
      !! test alone, the C2 fails on the atoms ("no image of an atom").
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: ax(3, 3, NMOL + 1), orig(3, NMOL + 1)
      integer :: nat
      call read_xyz("sym_water_tilt.xyz", nat, elem, r)
      res = sym_geom_detect(r, elem, err, requested="C2")
      call check(error, .not. err%has_error(), "no error: "//err%get_message())
      if (allocated(error)) return
      call check(error, res%detected, "C2v", "detected")
      if (allocated(error)) return
      call check(error, res%group%name, "C2", "reported")
      if (allocated(error)) return
      call check(error, res%residual < 1.0e-9_dp, "the tilted geometry is already symmetric")
      if (allocated(error)) return
      call pyscf_frames(ax, orig)
      call check_frame(error, res, ax(:, :, NMOL + 1), orig(:, NMOL + 1), "tilted water")
   end subroutine test_requested_tilted

   pure function hash_noise(amp, i, k, trial) result(x)
      !! A reproducible value in [-amp, amp] (no RNG: the same on every compiler).
      real(dp), intent(in) :: amp
      integer, intent(in) :: i, k, trial
      real(dp) :: x, t
      t = sin(12.9898_dp*i + 78.233_dp*k + 37.719_dp*trial)*43758.5453_dp
      x = amp*(2.0_dp*(t - floor(t)) - 1.0_dp)
   end function hash_noise

   pure function turn(axis, theta) result(m)
      !! The rotation by `theta` about `axis` (Rodrigues).
      real(dp), intent(in) :: axis(3), theta
      real(dp) :: m(3, 3), a(3), kx(3, 3)
      integer :: i
      a = axis/norm2(axis)
      kx = reshape([0.0_dp, a(3), -a(2), -a(3), 0.0_dp, a(1), a(2), -a(1), 0.0_dp], [3, 3])
      m = sin(theta)*kx + (1.0_dp - cos(theta))*matmul(kx, kx)
      do i = 1, 3
         m(i, i) = m(i, i) + 1.0_dp
      end do
   end function turn

   subroutine top_geometry(im, elem, r)
      !! NH3 (C3v), CH4 (Td) or SF6 (Oh), exact, in Bohr, turned off the input axes.
      integer, intent(in) :: im
      integer, allocatable, intent(out) :: elem(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      real(dp) :: th
      integer :: i
      select case (im)
      case (1)  ! NH3: N above the plane of the three H
         elem = [7, 1, 1, 1]
         allocate (r(3, 4))
         r(:, 1) = [0.0_dp, 0.0_dp, 0.72_dp]
         do i = 1, 3
            th = 2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/3.0_dp
            r(:, i + 1) = [1.77_dp*cos(th), 1.77_dp*sin(th), 0.0_dp]
         end do
      case (2)  ! CH4: H on alternate corners of a cube
         elem = [6, 1, 1, 1, 1]
         r = 1.19_dp*real(reshape([0, 0, 0, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1], [3, 5]), dp)
      case default  ! SF6: F on the axes
         elem = [16, 9, 9, 9, 9, 9, 9]
         r = 2.96_dp*real(reshape([0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1], [3, 7]), dp)
      end select
      r = matmul(turn([0.3_dp, -0.5_dp, 0.8_dp], 0.7_dp), r)
   end subroutine top_geometry

   subroutine test_noisy_tops(error)
      !! NH3, CH4 and SF6 with each coordinate moved by up to 1e-6 Bohr, well
      !! inside the default tol of 1e-5. pyscf 2.14.0 runs them in Cs, D2 and
      !! D2h. A threshold on equal moments that is not tied to tol splits the
      !! degenerate moments and gives C1, Cs or Ci in most trials. (At 3e-6 an
      !! F image is off by twice the noise plus the error of a frame found on
      !! noisy atoms, close to tol: pyscf itself keeps D2h in 81 of 200 such
      !! trials, terco in 155, so the check stays at 1e-6, where both always do.)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      integer, parameter :: NTRIAL = 10
      character(len=4), parameter :: NAME(3) = [character(len=4) :: "NH3", "CH4", "SF6"]
      character(len=3), parameter :: WANT(3) = [character(len=3) :: "Cs", "D2", "D2h"]
      real(dp), parameter :: AMP(3) = [1.0e-6_dp, 1.0e-6_dp, 1.0e-6_dp]
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r0(:, :), r(:, :)
      integer :: im, trial, i, k
      character(len=2) :: tr
      do im = 1, 3
         call top_geometry(im, elem, r0)
         res = sym_geom_detect(r0, elem, err)
         call check(error, res%group%name, WANT(im), trim(NAME(im))//" without noise")
         if (allocated(error)) return
         do trial = 1, NTRIAL
            r = r0
            do i = 1, size(r, 2)
               do k = 1, 3
                  r(k, i) = r(k, i) + hash_noise(AMP(im), i, k, trial + 100*im)
               end do
            end do
            res = sym_geom_detect(r, elem, err)
            write (tr, "(i0)") trial
            call check(error, .not. err%has_error(), trim(NAME(im))//" trial "//trim(tr)//": "//err%get_message())
            if (allocated(error)) return
            call check(error, res%group%name, WANT(im), trim(NAME(im))//" with noise, trial "//trim(tr))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_noisy_tops

   subroutine test_bent_linear(error)
      !! H-C-N with the H 1e-3 Bohr off the axis, the bend turned to several
      !! azimuths. Any three atoms are planar, so this is Cs whatever the
      !! direction. Its second moment off the axis (~1e-7 Bohr^2) is far above
      !! tol^2 but below tol: a linearity test against tol runs it as a linear
      !! molecule, whose mirrors pass or fail with the azimuth (Cs or C1). It
      !! has no pair of equivalent atoms for the top search's mirror candidates,
      !! so its plane comes from the principal frame.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      real(dp) :: r(3, 3), phi
      integer :: ia
      character(len=2) :: tag
      do ia = 0, 6
         phi = real(ia, dp)*acos(-1.0_dp)/12.0_dp
         r(:, 1) = [1.0e-3_dp*cos(phi), 1.0e-3_dp*sin(phi), -2.01_dp]
         r(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp]
         r(:, 3) = [0.0_dp, 0.0_dp, 2.19_dp]
         res = sym_geom_detect(r, [1, 6, 7], err)
         write (tag, "(i0)") 15*ia
         call check(error, .not. err%has_error(), "no error at azimuth "//trim(tag))
         if (allocated(error)) return
         call check(error, res%group%name, "Cs", "bent H-C-N at azimuth "//trim(tag)//" deg")
         if (allocated(error)) return
      end do
   end subroutine test_bent_linear

   subroutine test_large_tol(error)
      !! Four H atoms (Bohr) and tol = 0.75: under one sign mask of the
      !! principal frame two atoms have the same nearest image. Such a mask is
      !! not a symmetry: accepted, its `perm` is not a permutation and the
      !! symmetrised geometry averages the wrong atoms. Which group remains at this tol is not
      !! the point; that every operation permutes the atoms, and maps the
      !! symmetrised geometry onto itself, is.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      real(dp), parameter :: R(3, 4) = reshape([-1.66_dp, -1.05_dp, 0.0_dp, 0.33_dp, -1.62_dp, 0.0_dp, &
                                                -0.08_dp, -1.36_dp, 0.0_dp, -1.55_dp, -0.44_dp, 0.0_dp], [3, 4])
      logical :: seen(4)
      real(dp) :: y(3), worst
      integer :: k, i, dd
      res = sym_geom_detect(R, [1, 1, 1, 1], err, tol=0.75_dp)
      call check(error, .not. err%has_error(), "no error: "//err%get_message())
      if (allocated(error)) return
      do k = 1, res%group%order
         seen = .false.
         do i = 1, 4
            seen(res%perm(i, k)) = .true.
         end do
         call check(error, all(seen), "perm is a permutation")
         if (allocated(error)) return
      end do
      worst = 0.0_dp
      do k = 1, res%group%order
         do i = 1, 4
            y = res%r_std(:, i)
            do dd = 1, 3
               if (btest(res%group%ops(k), dd - 1)) y(dd) = -y(dd)
            end do
            worst = max(worst, maxval(abs(y - res%r_std(:, res%perm(i, k)))))
         end do
      end do
      call check(error, worst < 1.0e-12_dp, "g x = x[perm] on the symmetrised geometry")
   end subroutine test_large_tol

   subroutine test_bad_args(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: res
      type(error_t) :: err
      res = sym_geom_detect(reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), [1, 1], err)
      call check(error, err%has_error(), "size mismatch")
      if (allocated(error)) return
      res = sym_geom_detect(reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), [1], err, tol=0.0_dp)
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
   call run_testsuite(collect_sym_geom, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_geom
