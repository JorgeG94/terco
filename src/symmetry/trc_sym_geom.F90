module trc_sym_geom
   !! Point-group detection, orientation and symmetrization of a geometry.
   !!
   !! Algorithm. Shift to the centre of nuclear charge (the pyscf choice: a
   !! ghost atom weighs the charge of its element + 0.3, as in the pyscf
   !! `SymmSys`), diagonalize the charge-weighted second-moment tensor
   !! sum_i (w_i / W) d_i d_i^T (the "first Cartesian tensor" of the pyscf
   !! `detect_symm`), and, when its three eigenvalues are distinct, test
   !! all 7 sign masks in the principal frame: a mask is a symmetry when
   !! every atom maps onto an equivalent atom within `tol`. Two atoms are
   !! equivalent when element, ghost status and basis identity all match.
   !! The passing masks form the largest Abelian group of the D2h family
   !! (`detected`); the group actually reported (`group`) is capped by
   !! `SYM_GROUP_ENABLED`, or is the subgroup the caller `requested`.
   !!
   !! REQUESTED SUBGROUPS. With `requested` (a group name), `group` is that
   !! group when it is a proper subgroup of `detected` that has a rule in
   !! `has_subgroup_rule` (today C2v -> Cs, C2v -> C2, D2 -> C2, C2h -> Cs,
   !! C2h -> C2, C2h -> Ci, and D2h -> C2v, C2h, Cs, C2, Ci; not D2h -> D2, which
   !! pyscf refuses: D2 is not in its SUBGROUP table for D2h); `detected` still
   !! names the full group. The elements and the frame are those of pyscf
   !! 2.14.0 for `gto.M(..., symmetry='H')` (`Mole._build_symmetry`):
   !!  1. if the detected axes are the input axes (each a signed unit
   !!     vector, within 1e-5) and the elements of H hold in the input frame,
   !!     pyscf keeps the input frame (axes = identity). terco does the same:
   !!     `rot` is the identity (`input_frame_holds`);
   !!  2. else `as_subgroup`: the frame of the full group, turned for H. For
   !!     C2v -> Cs the axes are turned by 90 degrees about y (x becomes the
   !!     C2 axis, z the old x, the normal of the mirror kept); for
   !!     C2v -> C2 the C2v frame is kept as it is; for D2 -> C2 the D2 frame
   !!     is kept (z is the D2 axis closest to the input z, which is the C2 axis).
   !!     For every subgroup of C2h (Cs, C2, Ci) pyscf keeps the C2h frame as
   !!     it is (`get_subgroup` returns the C2h axes, `as_subgroup` changes
   !!     them for C2v only).
   !!     For every subgroup of D2h (C2v, C2h, Cs, C2, Ci) pyscf keeps the D2h
   !!     frame as it is (`as_subgroup` turns the axes for a C2v group name only,
   !!     and a D2h top group is named D2h): z is the D2h z, so the C2v and C2
   !!     axis, the C2h axis and the Cs normal are the same axis of the molecule.
   !!  Step 1 matters: the two mirrors of a nonplanar C2v molecule give
   !!  different SO blocks, and which one is kept depends on the input
   !!  orientation. Planar molecules have one choice, the molecular plane.
   !!
   !! HOW TO ADD A PAIR (full group G, subgroup H), the local change for the
   !! stages that enable a group: (a) make H's name true in
   !! `has_subgroup_rule(G, H)`; (b) in `frame_order`, give H's frame when
   !! `top == G` (the `top` argument is the detected group; copy what
   !! pyscf `as_subgroup` does for that pair: a rotation of G's axes, or
   !! none), reusing G's own axis order; (c) add a pyscf anchor for the
   !! pair to `check_sym_subgroup`. Nothing else changes: `sym_choose_group`
   !! accepts any enabled subgroup that `sym_is_subgroup` allows.
   !!
   !! Result: the reported group, the frame (`rot`, `shift`), the atom
   !! permutation `perm(natm, nop)` in the operation order of the group, the
   !! symmetrized coordinates in the standard frame (`r_std`) and in the
   !! input frame (`r_sym`), and the residual |r_sym - r|.
   !!
   !! Frame conventions follow pyscf 2.14.0 `symm/geom.py`, `detect_symm`
   !! and `as_subgroup` (all eight Abelian groups):
   !!  - C1: no rotation, no shift (pyscf resets both).
   !!  - Cs: z is the mirror normal; x, y are the other principal axes in
   !!    cyclic order (mirror normal P1 -> axes P2, P3, P1; and so on). A Cs
   !!    reached from a C2v molecule follows the pyscf `symmetry='Cs'`
   !!    subgroup of C2v: x is the C2 axis, and z is the normal of the
   !!    molecular plane when the molecule is planar.
   !!  - C2: z is the rotation axis; x, y are the other principal axes in
   !!    cyclic order. C2 has only C1 below it, so a molecule whose largest
   !!    Abelian group is C2 is reported as C2 (and only such a molecule: a
   !!    C2h or D2 is reported once enabled). A C2 requested on a C2v
   !!    molecule has the C2v axes, on a D2 molecule the D2 axes, on a C2h
   !!    molecule the C2h axes.
   !!  - D2: pyscf `detect_symm` renames the principal axes with
   !!    `alias_axes(axes, eye)`: z is the principal axis whose component
   !!    along the input z is largest (ties within 1e-5: the last of the
   !!    principal axes, which are in ascending order of the moment), then x
   !!    the same on the input x among the remaining two, and y the last.
   !!    Then pyscf swaps x and y when the three axes, with the signs LAPACK
   !!    gave them, are left-handed. That sign is not defined (it depends on
   !!    the LAPACK build), so for a molecule in a general orientation the
   !!    B2/B3 labels of pyscf are those of terco or exchanged; terco takes
   !!    the right-handed choice (no swap, y flipped). In an input whose
   !!    axes are the molecular axes the swap never happens (axes = identity).
   !!  - C2v: z is the C2 axis; x, y are the other principal axes in cyclic
   !!    order (C2 axis P1 -> x, y = P2, P3; P2 -> P3, P1; P3 -> P1, P2). A
   !!    planar molecule is turned about z so that x is the normal of its
   !!    plane (pyscf `_adjust_planar_c2v`): water lies in the yz plane, its
   !!    in-plane O-H bonding orbital is B2 and its out-of-plane lone pair B1.
   !!  - C2h: z is the C2 axis, which is also the normal of the mirror sigma_h
   !!    (xy); x, y are the other principal axes in cyclic order, as for C2
   !!    (pyscf's `Ag, Bg, Au, Bu`: x and y span Bu, z is Au, Rz is Ag). A Cs, C2
   !!    or Ci requested on a C2h molecule has this frame.
   !!  - D2h: the principal axes in ascending order of the second moment
   !!    (x, y, z = P1, P2, P3), with no `alias_axes` (unlike D2), so the frame
   !!    does not depend on LAPACK's eigenvector signs and a sign change of an axis
   !!    changes no D2h label. A planar molecule is turned by pyscf
   !!    `_adjust_planar_d2h` (`d2h_order`): x is the plane normal, and z the in-plane
   !!    axis chosen by counting the atoms off each axis. Measured on ethylene
   !!    (in the xy plane, C=C along x): pyscf axes (z_in, y_in, x_in). pyscf's
   !!    frame there is left-handed; terco flips y, which changes no label.
   !!  - Ci: the inversion fixes no axis, so the orientation is free; the
   !!    principal axes in order (1, 2, 3). The labels (Ag, Au) do not depend on it.
   !!    A Ci requested on a C2h molecule has the C2h frame.
   !! pyscf takes the sign of each axis from LAPACK, so it is arbitrary
   !! there; here the largest component of each principal axis is positive
   !! and y is flipped when needed so that det(rot) = +1. A sign flip of an
   !! axis changes no label of the groups above.
   !!
   !! Moments that are not all distinct (symmetric and spherical tops, linear
   !! molecules, single atoms): `trc_sym_tops` ports pyscf's axis search and
   !! `get_subgroup`, which give the Abelian group and the axes (S12). The frame
   !! is then used as it is (pyscf's `_symm_axes`, signs aside); only a requested
   !! C2v -> Cs turns it. Linear molecules and atoms accept no requested subgroup.
   use pic_types, only: dp
   use pic_logger, only: global_logger
   use trc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use trc_sym_group, only: trc_sym_group_t, sym_group_from_name, SYM_NGROUPS
   use trc_sym_tops, only: sym_top_axes, closest_axes_eye
   implicit none
   private

   public :: trc_sym_geom_t, sym_detect, sym_geom_atom_orbits
   public :: SYM_GROUP_ENABLED, SYM_DEFAULT_TOL

   logical, parameter :: SYM_GROUP_ENABLED(SYM_NGROUPS) = &
                         [.true., .true., .true., .true., .true., .true., .true., .true.]
      !! THE CAP: which groups may be reported, in the order of
      !! `sym_group_names` (D2h, C2h, C2v, D2, Cs, Ci, C2, C1). All eight Abelian
      !! groups are enabled since S11, so the cap no longer changes what is
      !! reported; it stays as the switch `enabled` and `frame_order` consult
      !! (`frame_order` refuses a group it has no convention for).

   integer, parameter :: N_MASK = 7  !! non-identity sign masks
   real(dp), parameter :: SYM_DEFAULT_TOL = 1.0e-5_dp  !! default geometry tolerance, Bohr
   real(dp), parameter :: DEGENERACY_TOL = 1.0e-6_dp
      !! eigenvalues closer than this fraction of the largest count as equal

   type :: trc_sym_geom_t
      type(trc_sym_group_t) :: group  !! the reported group
      character(len=3) :: detected = "C1"  !! the largest group found, before the cap
      logical :: degenerate = .false.  !! moments not all distinct (top, linear molecule or atom): `trc_sym_tops` frame
      real(dp) :: rot(3, 3) = 0.0_dp  !! rows are the frame axes: x_std = rot (x - shift)
      real(dp) :: shift(3) = 0.0_dp  !! origin of the frame (centre of nuclear charge)
      real(dp) :: residual = 0.0_dp  !! max |r_sym - r|, Bohr
      integer, allocatable :: perm(:, :)  !! perm(i, k): atom that the operation k sends atom i to
      real(dp), allocatable :: r_std(:, :)  !! (3, natm) symmetrized coordinates, standard frame
      real(dp), allocatable :: r_sym(:, :)  !! (3, natm) symmetrized coordinates, input frame
   end type trc_sym_geom_t

contains

   pure integer function sym_geom_atom_orbits(perm) result(norbit)
      !! Atoms of the molecule that are inequivalent under the group: the
      !! number of orbits of `perm` (`trc_sym_geom_t%perm`), the atom
      !! permutation of each operation. Orbits are found by relabelling each
      !! atom with the lowest index it reaches, so the count does not depend
      !! on the order of the operations. O(natm * nop) on tables the detection
      !! has already built; the skeleton Fock build's size rule
      !! (`trc_sym_run`'s sym_skeleton_size_rule) is the only caller, and it
      !! runs once per geometry.
      integer, intent(in) :: perm(:, :)
      integer :: lab(size(perm, 1))
      integer :: natm, nop, ia, k, lo, hi
      natm = size(perm, 1)
      nop = size(perm, 2)
      norbit = 0
      if (natm < 1) return
      do ia = 1, natm
         lab(ia) = ia
      end do
      ! A union over every (atom, image) edge, relabelling the whole array at
      ! each merge, so the labels are the connected components of the edge set
      ! whatever order the edges come in. The edges of all `nop` operations
      ! generate the orbit relation, so one pass over them is enough.
      do k = 1, nop
         do ia = 1, natm
            lo = min(lab(ia), lab(perm(ia, k)))
            hi = max(lab(ia), lab(perm(ia, k)))
            if (lo == hi) cycle
            where (lab == hi) lab = lo
         end do
      end do
      ! Every member of a component carries the component's lowest atom index,
      ! so the fixed points of `lab` count the orbits.
      do ia = 1, natm
         if (lab(ia) == ia) norbit = norbit + 1
      end do
   end function sym_geom_atom_orbits

   pure subroutine eig3(a, w, v)
      !! Eigenvalues `w` (ascending) and eigenvectors (columns of `v`) of a
      !! symmetric 3x3 matrix, by cyclic Jacobi rotations.
      real(dp), intent(in) :: a(3, 3)
      real(dp), intent(out) :: w(3)
      real(dp), intent(out) :: v(3, 3)
      real(dp) :: b(3, 3), j(3, 3), vt(3)
      real(dp) :: theta, t, c, s, off, scale
      integer :: sweep, p, q, i, k
      integer, parameter :: PQ(2, 3) = reshape([1, 2, 1, 3, 2, 3], [2, 3])

      b = a
      v = 0.0_dp
      do i = 1, 3
         v(i, i) = 1.0_dp
      end do
      scale = sum(b*b)
      do sweep = 1, 60
         off = b(1, 2)**2 + b(1, 3)**2 + b(2, 3)**2
         if (off <= 1.0e-32_dp*scale) exit
         do k = 1, 3
            p = PQ(1, k)
            q = PQ(2, k)
            if (b(p, q) == 0.0_dp) cycle
            theta = (b(q, q) - b(p, p))/(2.0_dp*b(p, q))
            t = sign(1.0_dp, theta)/(abs(theta) + sqrt(theta*theta + 1.0_dp))
            c = 1.0_dp/sqrt(t*t + 1.0_dp)
            s = t*c
            j = 0.0_dp
            do i = 1, 3
               j(i, i) = 1.0_dp
            end do
            j(p, p) = c
            j(q, q) = c
            j(p, q) = s
            j(q, p) = -s
            b = matmul(transpose(j), matmul(b, j))
            v = matmul(v, j)
         end do
      end do
      do i = 1, 3
         w(i) = b(i, i)
      end do
      ! sort ascending
      do i = 1, 2
         do k = i + 1, 3
            if (w(k) < w(i)) then
               t = w(i)
               w(i) = w(k)
               w(k) = t
               vt = v(:, i)
               v(:, i) = v(:, k)
               v(:, k) = vt
            end if
         end do
      end do
   end subroutine eig3

   pure function apply_mask(mask, x) result(y)
      !! The image of x under the operation with sign mask `mask`.
      integer, intent(in) :: mask
      real(dp), intent(in) :: x(3)
      real(dp) :: y(3)
      integer :: d
      y = x
      do d = 1, 3
         if (btest(mask, d - 1)) y(d) = -y(d)
      end do
   end function apply_mask

   pure function image_of(mask, q, cls, i, tol) result(j)
      !! The atom of the same class as atom i that lies within `tol` of its
      !! image under `mask` (the nearest one), or 0.
      integer, intent(in) :: mask, i
      integer, intent(in) :: cls(:)
      real(dp), intent(in) :: q(:, :)
      real(dp), intent(in) :: tol
      integer :: j
      real(dp) :: y(3)
      real(dp) :: dist, best
      integer :: m
      y = apply_mask(mask, q(:, i))
      j = 0
      best = tol
      do m = 1, size(q, 2)
         if (cls(m) /= cls(i)) cycle
         dist = maxval(abs(q(:, m) - y))
         if (dist < best) then
            best = dist
            j = m
         end if
      end do
   end function image_of

   pure function is_symmetry(mask, q, cls, tol) result(ok)
      integer, intent(in) :: mask
      integer, intent(in) :: cls(:)
      real(dp), intent(in) :: q(:, :)
      real(dp), intent(in) :: tol
      logical :: ok
      integer :: i
      ok = .true.
      do i = 1, size(q, 2)
         if (image_of(mask, q, cls, i, tol) == 0) then
            ok = .false.
            return
         end if
      end do
   end function is_symmetry

   pure function popcnt(m) result(n)
      integer, intent(in) :: m
      integer :: n, d
      n = 0
      do d = 0, 2
         if (btest(m, d)) n = n + 1
      end do
   end function popcnt

   pure function detected_name(found) result(name)
      !! Name of the group formed by the masks in `found` (identity not listed).
      logical, intent(in) :: found(N_MASK)
      character(len=3) :: name
      integer :: m, nmir, nc2, ninv
      nmir = 0
      nc2 = 0
      ninv = 0
      do m = 1, N_MASK
         if (.not. found(m)) cycle
         select case (popcnt(m))
         case (1)
            nmir = nmir + 1
         case (2)
            nc2 = nc2 + 1
         case default
            ninv = ninv + 1
         end select
      end do
      if (nmir + nc2 + ninv == 7) then
         name = "D2h"
      else if (nc2 == 3) then
         name = "D2"
      else if (nc2 == 1 .and. nmir == 2) then
         name = "C2v"
      else if (nc2 == 1 .and. ninv == 1) then
         name = "C2h"
      else if (nc2 == 1) then
         name = "C2"
      else if (nmir == 1) then
         name = "Cs"
      else if (ninv == 1) then
         name = "Ci"
      else
         name = "C1"
      end if
   end function detected_name

   pure function enabled(name) result(on)
      character(len=*), intent(in) :: name
      logical :: on
      integer :: ig
      character(len=3), parameter :: names(SYM_NGROUPS) = &
                                     [character(len=3) :: "D2h", "C2h", "C2v", "D2", "Cs", "Ci", "C2", "C1"]
      on = .false.
      do ig = 1, SYM_NGROUPS
         if (trim(names(ig)) == trim(name)) on = SYM_GROUP_ENABLED(ig)
      end do
   end function enabled

   function reported_name(detected, nmirror) result(name)
      !! The group to report: the detected one if enabled, else the largest
      !! enabled subgroup. Since S11 every Abelian group is enabled, so this is
      !! the detected group; the fallback to Cs (a mirror exists) or C1 only
      !! matters if a group is switched off in `SYM_GROUP_ENABLED`. Its frame
      !! would then be a Cs frame of the detected group, which only the groups
      !! with a single mirror define (`frame_order`).
      character(len=3), intent(in) :: detected
      integer, intent(in) :: nmirror
      character(len=3) :: name
      if (enabled(detected)) then
         name = detected
      else if (enabled("Cs") .and. nmirror >= 1) then
         name = "Cs"
      else
         name = "C1"
      end if
   end function reported_name

   pure function has_subgroup_rule(top, sub) result(ok)
      !! Whether a frame is defined for the proper subgroup `sub` of the
      !! detected group `top` (module header, "How to add a pair").
      character(len=*), intent(in) :: top, sub
      logical :: ok
      select case (trim(top))
      case ("C2v")
         ok = trim(sub) == "Cs" .or. trim(sub) == "C2"
      case ("D2")
         ok = trim(sub) == "C2"
      case ("C2h")
         ok = trim(sub) == "Cs" .or. trim(sub) == "C2" .or. trim(sub) == "Ci"
      case ("D2h")
         ! pyscf's SUBGROUP['D2h'] has no D2: the request is refused there too.
         ok = any(trim(sub) == ["C2v", "C2h", "Cs ", "C2 ", "Ci "])
      case default
         ok = .false.
      end select
   end function has_subgroup_rule

   pure function input_frame_holds(ax, ops, found) result(ok)
      !! pyscf keeps the input axes for a requested group when the detected
      !! axes are the input axes and every operation of the group holds there.
      !! `ax` are the principal axes (rows), `ops` the sign masks of the group
      !! in its own frame, `found` the masks of the detected group in the
      !! principal frame.
      real(dp), intent(in) :: ax(3, 3)
      integer, intent(in) :: ops(:)
      logical, intent(in) :: found(N_MASK)
      logical :: ok
      integer :: k(3)
      integer :: j, i, mq
      ok = .false.
      do j = 1, 3
         k(j) = maxloc(abs(ax(:, j)), 1)
         if (abs(abs(ax(k(j), j)) - 1.0_dp) > 1.0e-5_dp) return
      end do
      do i = 1, size(ops)
         if (ops(i) == 0) cycle
         mq = 0
         do j = 1, 3
            if (btest(ops(i), j - 1)) mq = ior(mq, ishft(1, k(j) - 1))
         end do
         if (mq < 1 .or. mq > N_MASK) return
         if (.not. found(mq)) return
      end do
      ok = .true.
   end function input_frame_holds

   pure function c2v_order(found, q, tol) result(order)
      !! The axis order of the C2v standard frame (module header).
      logical, intent(in) :: found(N_MASK)
      real(dp), intent(in) :: q(:, :)
      real(dp), intent(in) :: tol
      integer :: order(3)
      integer :: m, kc2, a0, a1
      kc2 = 0
      do m = 1, N_MASK
         if (found(m) .and. popcnt(m) == 2) kc2 = merge(1, merge(2, 3, .not. btest(m, 1)), .not. btest(m, 0))
      end do
      a0 = modulo(kc2, 3) + 1
      a1 = modulo(kc2 + 1, 3) + 1
      if (all(abs(q(a1, :)) < tol)) then  ! every atom on the plane whose normal would be y
         order = [a1, a0, kc2]
      else
         order = [a0, a1, kc2]
      end if
   end function c2v_order

   pure function d2_order(ax) result(order)
      !! The axis order of the D2 standard frame: pyscf `alias_axes` /
      !! `closest_axes` against the input axes (module header), shared with the
      !! O and T families in `trc_sym_tops`. `ax` are the principal axes (rows,
      !! input coordinates).
      real(dp), intent(in) :: ax(3, 3)
      integer :: order(3)
      order = closest_axes_eye(ax)
   end function d2_order

   pure function d2h_order(q) result(order)
      !! The axis order of the D2h standard frame: pyscf `detect_symm` keeps the
      !! principal axes (ascending second moment) and, for a planar molecule,
      !! applies `_adjust_planar_d2h`, which counts the atoms off each axis.
      !! `q` are the coordinates in the principal frame. Signs are not tracked:
      !! they change no D2h label. The counts use pyscf's tolerance.
      real(dp), intent(in) :: q(:, :)
      integer :: order(3)
      real(dp) :: tl
      integer :: nx, ny, nz
      tl = 1.0e-5_dp/sqrt(1.0_dp + real(size(q, 2), dp))
      nx = count(abs(q(1, :)) > tl)
      ny = count(abs(q(2, :)) > tl)
      nz = count(abs(q(3, :)) > tl)
      order = [1, 2, 3]
      if (nz == 0) then  ! atoms on the plane normal to P3
         order = merge([3, 2, 1], [3, 1, 2], nx >= ny)
      else if (ny == 0) then  ! atoms on the plane normal to P2
         order = merge([2, 1, 3], [2, 3, 1], nx < nz)
      else if (nx == 0) then  ! atoms on the plane normal to P1
         if (ny > nz) order = [1, 3, 2]
      end if
   end function d2h_order

   pure function c2h_order(found) result(order)
      !! The axis order of the C2h standard frame (module header): z is the C2
      !! axis (the axis the only C2 mask does not flip), x and y follow in
      !! cyclic order.
      logical, intent(in) :: found(N_MASK)
      integer :: order(3)
      integer :: m, kc2
      kc2 = 3
      do m = 1, N_MASK
         if (found(m) .and. popcnt(m) == 2) kc2 = merge(1, merge(2, 3, .not. btest(m, 1)), .not. btest(m, 0))
      end do
      order = [modulo(kc2, 3) + 1, modulo(kc2 + 1, 3) + 1, kc2]
   end function c2h_order

   subroutine frame_order(rep, top, found, q, tol, ax, order, error)
      !! Which principal axes become x, y, z of the standard frame of group
      !! `rep`, following pyscf (module header). `found` marks the masks of
      !! the detected group, `q` the coordinates in the principal frame.
      character(len=3), intent(in) :: rep
      character(len=3), intent(in) :: top  !! the detected group (for a requested subgroup)
      logical, intent(in) :: found(N_MASK)
      real(dp), intent(in) :: q(:, :)
      real(dp), intent(in) :: tol
      real(dp), intent(in) :: ax(3, 3)  !! the principal axes (rows, input coordinates)
      integer, intent(out) :: order(3)
      type(error_t), intent(out) :: error
      integer :: k, kc2, a0, a1, m, nmir
      logical :: planar

      order = [1, 2, 3]
      if (top == "D2h" .and. rep /= "C1") then
         ! pyscf keeps the D2h frame for D2h and every subgroup of it.
         order = d2h_order(q)
         return
      end if
      select case (rep)
      case ("C1")
         return
      case ("Ci")
         ! The inversion fixes no axis: the orientation is arbitrary. pyscf
         ! keeps the input orientation; the labels (Ag, Au) do not depend on
         ! it, so the principal axes are used as they come. A Ci requested on a
         ! C2h molecule has the C2h frame (pyscf as_subgroup keeps it).
         order = [1, 2, 3]
         if (top == "C2h") order = c2h_order(found)
      case ("Cs")
         nmir = 0
         k = 0
         kc2 = 0
         do m = 1, N_MASK
            if (.not. found(m)) cycle
            if (popcnt(m) == 1) then
               nmir = nmir + 1
               k = merge(1, merge(2, 3, m == 2), m == 1)  ! masks 1, 2, 4 flip x, y, z
            else if (popcnt(m) == 2) then
               ! the C2 axis is the one the mask does not flip
               kc2 = merge(1, merge(2, 3, .not. btest(m, 1)), .not. btest(m, 0))
            end if
         end do
         if (nmir == 1) then
            ! Cs, or C2h (whose only mirror is perpendicular to its axis):
            ! the mirror normal is z, the rest follow in cyclic order.
            order = [modulo(k, 3) + 1, modulo(k + 1, 3) + 1, k]
         else if (top == "C2v" .and. nmir == 2 .and. kc2 > 0) then
            ! the pyscf Cs subgroup of C2v: x = the C2 axis; the plane
            ! normal is the second remaining axis for a planar molecule.
            a0 = modulo(kc2, 3) + 1
            a1 = modulo(kc2 + 1, 3) + 1
            planar = all(abs(q(a1, :)) < tol)
            if (planar) then
               order = [kc2, a0, a1]
            else
               order = [kc2, a1, a0]
            end if
         else
            call error%set(ERROR_GENERIC, "frame_order: no Cs frame for a "//top//" molecule with "// &
                           "several mirrors")
            return
         end if
      case ("C2")
         ! pyscf's C2 frame has the axis along z; x and y follow in cyclic
         ! order (which two principal axes they are changes no A/B label).
         ! The rotation axis is the one the (only) C2 mask does not flip.
         kc2 = 0
         do m = 1, N_MASK
            if (found(m) .and. popcnt(m) == 2) kc2 = merge(1, merge(2, 3, .not. btest(m, 1)), .not. btest(m, 0))
         end do
         if (kc2 == 0) then
            call error%set(ERROR_GENERIC, "frame_order: group C2 reported without a C2 axis")
            return
         end if
         if (top == "D2") then
            order = d2_order(ax)  ! pyscf as_subgroup keeps the D2 axes
         else if (top == "C2v") then
            order = c2v_order(found, q, tol)  ! pyscf as_subgroup keeps the C2v axes
         else
            order = [modulo(kc2, 3) + 1, modulo(kc2 + 1, 3) + 1, kc2]
         end if
      case ("C2v")
         order = c2v_order(found, q, tol)
      case ("D2")
         order = d2_order(ax)
      case ("C2h")
         order = c2h_order(found)
      case default
         call error%set(ERROR_GENERIC, "no axis convention yet for point group "//rep// &
                        " (enabled in trc_sym_geom, but its frame is missing)")
      end select
   end subroutine frame_order

   function sym_detect(r, elem, error, tol, ghost, basis_id, requested) result(res)
      !! Detect, orient and symmetrize the geometry `r` (Bohr, 3 x natm) of
      !! atoms with element numbers `elem`. `ghost` marks ghost atoms and
      !! `basis_id` labels the basis of each atom; atoms are equivalent only when
      !! element, ghost status and basis label all agree. `tol` is the
      !! geometry tolerance in Bohr (default `SYM_DEFAULT_TOL`). `requested`
      !! names a subgroup of the detected group to report instead (module
      !! header); it is ignored when it is not a proper subgroup with a rule,
      !! and the caller (`sym_choose_group`) then decides what that means.
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: elem(:)
      type(error_t), intent(out) :: error
      real(dp), intent(in), optional :: tol
      logical, intent(in), optional :: ghost(:)
      integer, intent(in), optional :: basis_id(:)
      character(len=*), intent(in), optional :: requested
      type(trc_sym_geom_t) :: res

      integer :: natm, i, j, k, m, n, nmir
      integer, allocatable :: cls(:)
      real(dp), allocatable :: w(:), d(:, :), q(:, :), qf(:, :)
      real(dp) :: t(3, 3), vec(3, 3), ax(3, 3)
      real(dp) :: ev(3), cm(3)
      real(dp) :: tl, wsum, det, big
      integer :: order(3)
      logical :: found(N_MASK)
      character(len=3) :: rep, abel
      character(len=4) :: topname
      logical :: subreq, top_frame, eq12, eq23, top
      integer :: kind
      real(dp) :: tolp

      natm = size(r, 2)
      tl = SYM_DEFAULT_TOL
      if (present(tol)) tl = tol
      if (size(r, 1) /= 3 .or. natm < 1 .or. size(elem) /= natm .or. tl <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "sym_detect: need r(3, natm), elem(natm) with natm >= 1 and tol > 0")
         return
      end if
      if (present(ghost)) then
         if (size(ghost) /= natm) call error%set(ERROR_VALIDATION, "sym_detect: ghost has the wrong size")
      end if
      if (present(basis_id)) then
         if (size(basis_id) /= natm) call error%set(ERROR_VALIDATION, "sym_detect: basis_id has the wrong size")
      end if
      if (error%has_error()) return

      allocate (cls(natm), w(natm))
      do i = 1, natm
         cls(i) = elem(i)
         w(i) = real(elem(i), dp)
         if (present(ghost)) then
            if (ghost(i)) then
               cls(i) = cls(i) + 1000
               w(i) = w(i) + 0.3_dp
            end if
         end if
         if (present(basis_id)) cls(i) = cls(i) + 100000*basis_id(i)
      end do
      wsum = sum(w)
      if (wsum <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "sym_detect: total nuclear charge must be positive")
         return
      end if

      allocate (res%perm(natm, 1), res%r_std(3, natm), res%r_sym(3, natm))
      res%group = sym_group_from_name("C1", error)
      res%perm(:, 1) = [(i, i = 1, natm)]
      res%rot = 0.0_dp
      do i = 1, 3
         res%rot(i, i) = 1.0_dp
      end do
      res%r_std = r
      res%r_sym = r

      cm = 0.0_dp
      do i = 1, natm
         cm = cm + w(i)*r(:, i)
      end do
      cm = cm/wsum
      allocate (d(3, natm), q(3, natm), qf(3, natm))
      do i = 1, natm
         d(:, i) = r(:, i) - cm
      end do
      t = 0.0_dp
      do i = 1, natm
         do j = 1, 3
            do k = 1, 3
               t(j, k) = t(j, k) + w(i)/wsum*d(j, i)*d(k, i)
            end do
         end do
      end do
      call eig3(t, ev, vec)

      big = ev(3)
      tolp = tl/sqrt(1.0_dp + real(natm, dp))  ! pyscf: TOLERANCE / sqrt(1 + natm)
      eq12 = ev(2) - ev(1) <= DEGENERACY_TOL*big
      eq23 = ev(3) - ev(2) <= DEGENERACY_TOL*big
      top = big <= tolp .or. ev(2) <= tolp .or. eq12 .or. eq23
      top_frame = .false.
      if (top) then
         ! symmetric or spherical top, linear molecule, atom: the pyscf search of trc_sym_tops
         res%degenerate = .true.
         call sym_top_axes(d, cls, tl, ev, vec, eq12, eq23, kind, topname, abel, ax)
         if (kind == 0) then
            call global_logger%warning("symmetry: the point group of this symmetric or spherical top is not "// &
                                       "one pyscf reduces to an Abelian group; using C1")
            return
         end if
         top_frame = kind == 1
      else
         ! Principal axes, sign fixed: the largest component is positive.
         do k = 1, 3
            ax(k, :) = vec(:, k)
            m = maxloc(abs(ax(k, :)), 1)
            do j = 1, 3
               if (abs(ax(k, j)) >= abs(ax(k, m)) - 1.0e-6_dp) then
                  m = j
                  exit
               end if
            end do
            if (ax(k, m) < 0.0_dp) ax(k, :) = -ax(k, :)
         end do
      end if
      do i = 1, natm
         q(:, i) = matmul(ax, d(:, i))
      end do

      found = .false.
      do m = 1, N_MASK
         found(m) = is_symmetry(m, q, cls, tl)
      end do
      res%detected = detected_name(found)
      if (top_frame) then
         if (res%detected /= abel) then
            ! the verified group differs from pyscf's (tolerance): use the generic frame rules
            call global_logger%warning("symmetry: the group "//trim(abel)//" that the pyscf search found does not "// &
                                       "hold to the tolerance; using "//trim(res%detected))
            top_frame = .false.
         end if
      end if
      nmir = count([(found(m) .and. popcnt(m) == 1, m=1, N_MASK)])
      rep = reported_name(res%detected, nmir)
      subreq = .false.
      if (present(requested)) then
         if (trim(requested) /= trim(rep) .and. trim(requested) /= "C1" .and. enabled(requested) .and. &
             has_subgroup_rule(res%detected, requested) .and. &
             .not. (top .and. (topname == "SO3" .or. topname == "Dooh" .or. topname == "Coov"))) then
            rep = trim(requested)
            subreq = .true.
         end if
      end if
      res%group = sym_group_from_name(rep, error)
      if (error%has_error()) return

      if (rep == "C1") return  ! identity frame, no shift: the old code path sees the input untouched

      if (top_frame) then
         ! `ax` is already pyscf's `_symm_axes`; only C2v -> Cs turns it (as_subgroup, l.416)
         order = [1, 2, 3]
         if (res%detected == "C2v" .and. rep == "Cs") order = [3, 2, 1]
      else
         call frame_order(rep, res%detected, found, q, tl, ax, order, error)
         if (error%has_error()) return
      end if
      res%rot = ax(order, :)
      det = res%rot(1, 1)*(res%rot(2, 2)*res%rot(3, 3) - res%rot(2, 3)*res%rot(3, 2)) &
              - res%rot(1, 2)*(res%rot(2, 1)*res%rot(3, 3) - res%rot(2, 3)*res%rot(3, 1)) &
              + res%rot(1, 3)*(res%rot(2, 1)*res%rot(3, 2) - res%rot(2, 2)*res%rot(3, 1))
      if (det < 0.0_dp) res%rot(2, :) = -res%rot(2, :)
      if (subreq) then
         if (input_frame_holds(ax, res%group%ops, found)) then  ! pyscf keeps the input axes
            res%rot = 0.0_dp
            do i = 1, 3
               res%rot(i, i) = 1.0_dp
            end do
         end if
      end if
      res%shift = cm
      do i = 1, natm
         qf(:, i) = matmul(res%rot, d(:, i))
      end do

      n = res%group%order
      deallocate (res%perm)
      allocate (res%perm(natm, n))
      do k = 1, n
         do i = 1, natm
            j = image_of(res%group%ops(k), qf, cls, i, tl)
            if (j == 0) then
               call error%set(ERROR_GENERIC, "sym_detect: an operation of the detected group has no image "// &
                              "of an atom in the standard frame")
               return
            end if
            res%perm(i, k) = j
         end do
      end do

      ! Symmetrize: x_i <- mean_k g_k x_perm(i,k). Every g_k is its own inverse.
      do i = 1, natm
         res%r_std(:, i) = 0.0_dp
         do k = 1, n
            res%r_std(:, i) = res%r_std(:, i) + apply_mask(res%group%ops(k), qf(:, res%perm(i, k)))
         end do
         res%r_std(:, i) = res%r_std(:, i)/real(n, dp)
      end do
      do i = 1, natm
         res%r_sym(:, i) = res%shift + matmul(transpose(res%rot), res%r_std(:, i))
      end do
      res%residual = maxval(abs(res%r_sym - r))
   end function sym_detect

end module trc_sym_geom
