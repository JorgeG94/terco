module trc_sym_tops
   !! Point groups of molecules whose principal moments are not all distinct
   !! (symmetric and spherical tops, linear molecules, single atoms): a port of
   !! the degenerate branch of pyscf 2.14.0 `pyscf/symm/geom.py`, so that terco
   !! runs such a molecule in the Abelian group and frame pyscf gives it.
   !!
   !! What is ported, with the pyscf functions (line numbers of 2.14.0):
   !!  - `detect_symm` (l.188-346): the dispatch SO3 / Dooh, Coov / T, O, I /
   !!    symmetric top, and the naming of the full group (Dn, Dnh, Dnd, Cnv, Cnh,
   !!    Sn, Cn) from the axis search;
   !!  - `SymmSys` (l.~600-812): the atoms grouped by class and rounded distance
   !!    (`group_atoms_by_distance`), `symmetric_for`, `has_icenter`,
   !!    `has_rotation`, `has_mirror`, `has_improper_rotation`,
   !!    `search_possible_rotations`, `search_c2x`, `search_mirrorx`,
   !!    `search_c_highest`, with the same candidate order and tolerances;
   !!  - `_search_ot_group`, `_search_i_group`, `_remove_dupvec`,
   !!    `_pseudo_vectors`, `_make_axes`, `_refine` (l.824-959);
   !!  - `get_subgroup` (l.347-415): the full group to the Abelian group pyscf
   !!    uses (`mol.groupname`), with its axes; `alias_axes` / `closest_axes`
   !!    (l.118-141) for the O, T families.
   !!
   !! NOT ported, and why:
   !!  - the 2nd and 3rd order Casimir tensors that `detect_symm` uses to gate
   !!    the T/O/I and symmetric-top branches (`w2`, `w3` degeneracies): the gate
   !!    only decides whether to try the search, and the search itself verifies
   !!    every axis with `has_rotation`. The difference is a molecule with
   !!    accidentally equal moments, where pyscf runs on LAPACK's arbitrary axes
   !!    of the degenerate plane and this port searches properly;
   !!  - the case `n = -1` of `detect_symm` (no equal moments): that is the
   !!    existing terco path in `trc_sym_geom`;
   !!  - requested subgroups of the cylindrical groups (see `trc_sym_geom`);
   !!  - groups that pyscf's `get_subgroup` maps to a name it cannot use (S6,
   !!    S10, ... give 'C3', 'C5'): `kind = 0` is returned and the caller keeps C1.
   !!
   !! LAPACK. pyscf takes the unique axis of a symmetric top from `eigh`, whose
   !! sign is arbitrary. `search_c2x` (candidates of nonzero height) and
   !! `search_mirrorx` depend on that sign, so for a top whose perpendicular C2
   !! axes or mirrors are not all equivalent by symmetry the choice between
   !! equivalent ones can differ from pyscf's on another LAPACK. The labels and
   !! SO dimensions do not depend on it. Here the sign is fixed (the largest
   !! component of the unique axis is positive).
   !!
   !! The linear molecule. pyscf takes z along the molecule and x, y from `eigh`
   !! (arbitrary in the plane). Here x is the projection of the input axis with
   !! the smallest component along z, so an input along z gives the identity.
   !! The group is D2h with an inversion centre, C2v without (`get_subgroup`
   !! keeps 'Dooh' / 'Coov' for the spherical basis, whose irreps descend to
   !! D2h / C2v by `irrep_id % 10`; with a Cartesian basis pyscf itself uses
   !! D2h / C2v). A single atom is D2h with the input axes ('SO3' descends to
   !! D2h for a Cartesian basis).
   use pic_types, only: dp
   use pic_sorting, only: sort
   implicit none
   private

   public :: sym_top_axes, closest_axes_eye

   real(dp), parameter :: PI = 3.14159265358979323846_dp

   type :: sys_t
      !! pyscf `SymmSys`: atoms about the charge centre, grouped by class and
      !! rounded distance.
      integer :: natm = 0
      real(dp) :: tol = 1.0e-5_dp  !! pyscf `TOLERANCE`
      real(dp), allocatable :: r(:, :)  !! (3, natm) about the charge centre
      integer :: ngrp = 0
      integer, allocatable :: gofs(:)  !! group g is gidx(gofs(g)+1 : gofs(g+1))
      integer, allocatable :: gidx(:)
   end type sys_t

contains

   pure function normalize(v) result(u)
      !! pyscf `_normalize` of one vector (l.813).
      real(dp), intent(in) :: v(3)
      real(dp) :: u(3)
      u = v/(sqrt(sum(v*v)) + 1.0e-200_dp)
   end function normalize

   pure function cross(a, b) result(c)
      real(dp), intent(in) :: a(3), b(3)
      real(dp) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   end function cross

   pure function rotation_mat(vec, theta) result(r)
      !! pyscf `rotation_mat` (l.99): new = R old.
      real(dp), intent(in) :: vec(3)
      real(dp), intent(in) :: theta
      real(dp) :: r(3, 3)
      real(dp) :: v(3), ux(3, 3)
      real(dp) :: c, s
      integer :: i
      v = normalize(vec)
      ux = reshape([0.0_dp, v(3), -v(2), -v(3), 0.0_dp, v(1), v(2), -v(1), 0.0_dp], [3, 3])
      c = cos(theta)
      s = sin(theta)
      r = s*ux
      do i = 1, 3
         r(:, i) = r(:, i) + (1.0_dp - c)*v*v(i)
         r(i, i) = r(i, i) + c
      end do
   end function rotation_mat

   pure function householder(vec) result(h)
      !! pyscf `householder` (l.114).
      real(dp), intent(in) :: vec(3)
      real(dp) :: h(3, 3)
      real(dp) :: v(3)
      integer :: i
      v = normalize(vec)
      h = 0.0_dp
      do i = 1, 3
         h(:, i) = -2.0_dp*v*v(i)
         h(i, i) = h(i, i) + 1.0_dp
      end do
   end function householder

   pure function parallel_vectors(v1, v2, tol) result(p)
      !! pyscf `parallel_vectors` (l.64).
      real(dp), intent(in) :: v1(3), v2(3)
      real(dp), intent(in) :: tol
      logical :: p
      real(dp) :: c
      if (all(abs(v1) <= tol) .or. all(abs(v2) <= tol)) then
         p = .true.
      else
         c = dot_product(normalize(v1), normalize(v2))
         p = abs(c - 1.0_dp) < tol .or. abs(c + 1.0_dp) < tol
      end if
   end function parallel_vectors

   pure function vec_in_vecs(vec, vecs, tol) result(found)
      !! pyscf `_vec_in_vecs` (l.820): the L1 distance to the nearest vector,
      !! divided by sqrt(number of vectors), is below the tolerance.
      real(dp), intent(in) :: vec(3), vecs(:, :)
      real(dp), intent(in) :: tol
      logical :: found
      real(dp) :: best
      integer :: i
      best = huge(1.0_dp)
      do i = 1, size(vecs, 2)
         best = min(best, sum(abs(vecs(:, i) - vec)))
      end do
      found = best/sqrt(real(size(vecs, 2), dp)) < tol
   end function vec_in_vecs

   subroutine build_sys(d, cls, tol, s)
      !! pyscf `SymmSys.__init__` (l.~600): the groups of equal class and
      !! rounded distance from the charge centre, classes in order of first
      !! appearance, distances ascending (`np.unique`).
      real(dp), intent(in) :: d(:, :)
      integer, intent(in) :: cls(:)
      real(dp), intent(in) :: tol
      type(sys_t), intent(out) :: s
      integer :: natm, i, j, ic, ncls, g, ng, nin
      integer, allocatable :: ucls(:), key(:), members(:), uk(:)
      real(dp) :: scale
      natm = size(d, 2)
      s%natm = natm
      s%tol = tol
      s%r = d
      scale = 10.0_dp**(int(-log10(tol)) - 1)  ! decimals = int(-log10(TOLERANCE)) - 1
      allocate (ucls(natm), key(natm), members(natm), s%gofs(natm + 1), s%gidx(natm), uk(natm))
      do i = 1, natm
         key(i) = nint(sqrt(sum(d(:, i)**2))*scale)
      end do
      ncls = 0
      do i = 1, natm
         if (.not. any(ucls(1:ncls) == cls(i))) then
            ncls = ncls + 1
            ucls(ncls) = cls(i)
         end if
      end do
      g = 0
      s%gofs(1) = 0
      do ic = 1, ncls
         nin = 0
         do i = 1, natm
            if (cls(i) == ucls(ic)) then
               nin = nin + 1
               members(nin) = i
            end if
         end do
         ! distinct rounded distances, ascending
         ng = 0
         do i = 1, nin
            if (.not. any(uk(1:ng) == key(members(i)))) then
               ng = ng + 1
               uk(ng) = key(members(i))
            end if
         end do
         call sort(uk(1:ng))
         do j = 1, ng
            g = g + 1
            s%gofs(g + 1) = s%gofs(g)
            do i = 1, nin
               if (key(members(i)) == uk(j)) then
                  s%gofs(g + 1) = s%gofs(g + 1) + 1
                  s%gidx(s%gofs(g + 1)) = members(i)
               end if
            end do
         end do
      end do
      s%ngrp = g
   end subroutine build_sys

   function group_coords(s, g) result(c)
      type(sys_t), intent(in) :: s
      integer, intent(in) :: g
      real(dp), allocatable :: c(:, :)
      c = s%r(:, s%gidx(s%gofs(g) + 1:s%gofs(g + 1)))
   end function group_coords

   function symmetric_for(s, op) result(ok)
      !! pyscf `symmetric_for` (l.~708): every image r.op of every atom of every
      !! group is an atom of the group.
      type(sys_t), intent(in) :: s
      real(dp), intent(in) :: op(3, 3)
      logical :: ok
      real(dp), allocatable :: r0(:, :)
      integer :: g, i
      ok = .true.
      do g = 1, s%ngrp
         r0 = group_coords(s, g)
         do i = 1, size(r0, 2)
            if (.not. vec_in_vecs(matmul(r0(:, i), op), r0, s%tol)) then
               ok = .false.
               return
            end if
         end do
      end do
   end function symmetric_for

   function has_icenter(s) result(ok)
      type(sys_t), intent(in) :: s
      logical :: ok
      real(dp) :: op(3, 3)
      integer :: i
      op = 0.0_dp
      do i = 1, 3
         op(i, i) = -1.0_dp
      end do
      ok = symmetric_for(s, op)
   end function has_icenter

   function has_rotation(s, axis, n) result(ok)
      type(sys_t), intent(in) :: s
      real(dp), intent(in) :: axis(3)
      integer, intent(in) :: n
      logical :: ok
      ok = symmetric_for(s, transpose(rotation_mat(axis, 2.0_dp*PI/real(n, dp))))
   end function has_rotation

   function has_mirror(s, perp) result(ok)
      type(sys_t), intent(in) :: s
      real(dp), intent(in) :: perp(3)
      logical :: ok
      ok = symmetric_for(s, transpose(householder(perp)))
   end function has_mirror

   function has_improper_rotation(s, axis, n) result(ok)
      type(sys_t), intent(in) :: s
      real(dp), intent(in) :: axis(3)
      integer, intent(in) :: n
      logical :: ok
      ok = symmetric_for(s, transpose(matmul(householder(axis), rotation_mat(axis, PI/real(n, dp)))))
   end function has_improper_rotation

   subroutine search_possible_rotations(s, zaxis, vk, nk, nv)
      !! pyscf `search_possible_rotations` (l.~716): candidate axes and orders
      !! from the atoms of equal distance (C2 from pairs, Cn from equal-distance
      !! triples), duplicates merged. `zaxis`, when present, keeps the parallel
      !! ones.
      type(sys_t), intent(in) :: s
      real(dp), intent(in), optional :: zaxis(3)
      real(dp), allocatable, intent(out) :: vk(:, :)
      integer, allocatable, intent(out) :: nk(:)
      integer, intent(out) :: nv
      real(dp), allocatable :: mv(:, :), vecs(:, :), r0(:, :), c(:, :)
      integer, allocatable :: mn(:), ns(:)
      real(dp), allocatable :: dist(:)
      logical, allocatable :: seen(:), w1(:), w2(:)
      real(dp) :: v(3), vs(3)
      real(dp) :: tol, cosv, ang, nfrac
      integer :: g, i, j, k, ng, n, cap, nm, nvec, nmax_n
      logical :: has_n

      tol = s%tol
      cap = 1
      do g = 1, s%ngrp
         ng = s%gofs(g + 1) - s%gofs(g)
         cap = cap + ng*ng + ng
      end do
      allocate (mv(3, cap), mn(cap))
      nm = 0
      do g = 1, s%ngrp
         ng = s%gofs(g + 1) - s%gofs(g)
         if (ng <= 1) cycle
         c = group_coords(s, g)
         do i = 2, ng
            nm = nm + 1
            mn(nm) = 2
            if (sum(abs(c(:, 1) + c(:, i))) > tol) then
               mv(:, nm) = c(:, 1) + c(:, i)
            else
               mv(:, nm) = c(:, 1) - c(:, i)
            end if
         end do
         allocate (r0(3, ng), dist(ng))
         do i = 1, ng
            r0(:, i) = c(:, i) - c(:, 1)
            dist(i) = sqrt(sum(r0(:, i)**2))
         end do
         do i = 3, ng
            do j = 1, i - 1
               if (abs(dist(i) - dist(j)) >= tol) cycle
               cosv = dot_product(r0(:, i), r0(:, j))/(dist(i)*dist(j))
               cosv = max(-1.0_dp, min(1.0_dp, cosv))
               ang = acos(cosv)
               nfrac = 2.0_dp*PI/(PI - ang)
               n = nint(nfrac)
               if (abs(nfrac - real(n, dp)) < tol) then
                  nm = nm + 1
                  mv(:, nm) = cross(r0(:, i), r0(:, j))
                  mn(nm) = n
               end if
            end do
         end do
         deallocate (r0, dist)
      end do

      ! remove zero vectors, normalise
      allocate (vecs(3, max(nm, 1)), ns(max(nm, 1)))
      nvec = 0
      do i = 1, nm
         if (sqrt(sum(mv(:, i)**2)) > tol) then
            nvec = nvec + 1
            vecs(:, nvec) = normalize(mv(:, i))
            ns(nvec) = mn(i)
         end if
      end do
      if (present(zaxis)) then
         v = normalize(zaxis)
         k = 0
         do i = 1, nvec
            cosv = dot_product(vecs(:, i), v)
            if (abs(cosv - 1.0_dp) < tol .or. abs(cosv + 1.0_dp) < tol) then
               k = k + 1
               vecs(:, k) = vecs(:, i)
               ns(k) = ns(i)
            end if
         end do
         nvec = k
      end if

      allocate (vk(3, max(nvec*4, 1)), nk(max(nvec*4, 1)), seen(max(nvec, 1)), w1(max(nvec, 1)), w2(max(nvec, 1)))
      nv = 0
      seen = .false.
      do k = 1, nvec
         if (seen(k)) cycle
         v = vecs(:, k)
         w1 = .false.
         w2 = .false.
         do i = k, nvec
            w1(i) = sum(abs(vecs(:, i) - v)) < tol
            w2(i) = sum(abs(vecs(:, i) + v)) < tol
         end do
         seen = seen .or. w1 .or. w2
         vs = 0.0_dp
         do i = k, nvec
            if (w1(i)) vs = vs + vecs(:, i)
            if (w2(i)) vs = vs - vecs(:, i)
         end do
         vs = normalize(vs)
         ! the orders, ascending (a python set of small integers)
         nmax_n = 0
         do i = k, nvec
            if (w1(i) .or. w2(i)) nmax_n = max(nmax_n, ns(i))
         end do
         do n = 1, nmax_n
            has_n = .false.
            do i = k, nvec
               if ((w1(i) .or. w2(i)) .and. ns(i) == n) has_n = .true.
            end do
            if (has_n) then
               nv = nv + 1
               if (nv > size(nk)) call grow(vk, nk)
               vk(:, nv) = vs
               nk(nv) = n
            end if
         end do
      end do
   end subroutine search_possible_rotations

   subroutine grow(vk, nk)
      real(dp), allocatable, intent(inout) :: vk(:, :)
      integer, allocatable, intent(inout) :: nk(:)
      real(dp), allocatable :: tv(:, :)
      integer, allocatable :: tn(:)
      integer :: n
      n = size(nk)
      allocate (tv(3, 2*n), tn(2*n))
      tv(:, 1:n) = vk
      tn(1:n) = nk
      call move_alloc(tv, vk)
      call move_alloc(tn, nk)
   end subroutine grow

   subroutine search_c_highest(s, cmax, nmax, zaxis)
      !! pyscf `search_c_highest` (l.~805).
      type(sys_t), intent(in) :: s
      real(dp), intent(out) :: cmax(3)
      integer, intent(out) :: nmax
      real(dp), intent(in), optional :: zaxis(3)
      real(dp), allocatable :: vk(:, :)
      integer, allocatable :: nk(:)
      integer :: nv, i
      nmax = 1
      cmax = [0.0_dp, 0.0_dp, 1.0_dp]
      if (all(s%gofs(2:s%ngrp + 1) - s%gofs(1:s%ngrp) <= 1)) return  ! no pair of equivalent atoms
      call search_possible_rotations(s, zaxis, vk, nk, nv)
      do i = 1, nv
         if (nk(i) > nmax) then
            if (has_rotation(s, vk(:, i), nk(i))) then
               nmax = nk(i)
               cmax = vk(:, i)
            end if
         end if
      end do
   end subroutine search_c_highest

   pure subroutine pseudo_vectors(vs)
      !! pyscf `_pseudo_vectors` (l.893): z > 0; else y > 0; else x > 0.
      real(dp), intent(inout) :: vs(:, :)
      real(dp), parameter :: TOL = 1.0e-5_dp
      logical, allocatable :: idy0(:), idz0(:)
      integer :: i
      allocate (idy0(size(vs, 2)), idz0(size(vs, 2)))
      idy0 = abs(vs(2, :)) < TOL
      idz0 = abs(vs(3, :)) < TOL
      do i = 1, size(vs, 2)
         if (vs(3, i) < 0.0_dp) vs(:, i) = -vs(:, i)
      end do
      do i = 1, size(vs, 2)
         if (vs(2, i) < 0.0_dp .and. idz0(i)) vs(:, i) = -vs(:, i)
      end do
      do i = 1, size(vs, 2)
         if (vs(1, i) < 0.0_dp .and. idy0(i) .and. idz0(i)) vs(:, i) = -vs(:, i)
      end do
   end subroutine pseudo_vectors

   subroutine remove_dupvec(vs, n, tol)
      !! pyscf `_remove_dupvec` (l.905): signs by `_pseudo_vectors`, then keep
      !! each vector that differs from every earlier kept one (L1 > tol).
      real(dp), intent(inout) :: vs(:, :)
      integer, intent(inout) :: n
      real(dp), intent(in) :: tol
      real(dp), allocatable :: kept(:, :)
      integer :: i, j, nk
      logical :: dup
      if (n <= 0) return
      call pseudo_vectors(vs(:, 1:n))
      allocate (kept(3, n))
      nk = 0
      do i = 1, n
         dup = .false.
         do j = 1, nk
            if (sum(abs(vs(:, i) - kept(:, j))) <= tol) dup = .true.
         end do
         if (.not. dup) then
            nk = nk + 1
            kept(:, nk) = vs(:, i)
         end if
      end do
      vs(:, 1:nk) = kept(:, 1:nk)
      n = nk
   end subroutine remove_dupvec

   subroutine search_c2x(s, zaxis, n, found, c2x)
      !! pyscf `search_c2x` (l.~745): a C2 axis perpendicular to z, from the
      !! height classes of each group; the first that passes `has_rotation`.
      type(sys_t), intent(in) :: s
      real(dp), intent(in) :: zaxis(3)
      integer, intent(in) :: n
      logical, intent(out) :: found
      real(dp), intent(out) :: c2x(3)
      real(dp), allocatable :: r0(:, :), mc(:, :)
      integer, allocatable :: zk(:), uk(:)
      real(dp) :: r1(3), r2(3), z(3)
      real(dp) :: scale, nrm
      integer :: g, ng, i, j, ia, ib, nu, nm, d, cnt_above, cnt_below, nfin, ncand
      logical :: want
      real(dp), allocatable :: cand(:, :)

      found = .false.
      c2x = 0.0_dp
      scale = 10.0_dp**(int(-log10(s%tol)) - 1)
      z = zaxis
      do g = 1, s%ngrp
         ng = s%gofs(g + 1) - s%gofs(g)
         if (ng <= 1) cycle
         r0 = group_coords(s, g)
         allocate (zk(ng), uk(ng), mc(3, ng + 2))
         do i = 1, ng
            zk(i) = nint(dot_product(r0(:, i), z)*scale)
         end do
         nu = 0
         do i = 1, ng
            if (.not. any(uk(1:nu) == zk(i))) then
               nu = nu + 1
               uk(nu) = zk(i)
            end if
         end do
         call sort(uk(1:nu))
         nm = 0
         do ia = 1, nu
            d = uk(ia)
            if (d > 0) then
               cnt_above = count(zk == d)
               cnt_below = count(zk == -d)
               if (cnt_below == cnt_above) then
                  ! above[0] + below[i]
                  do i = 1, ng
                     if (zk(i) == d) exit
                  end do
                  do j = 1, ng
                     if (zk(j) == -d) then
                        nm = nm + 1
                        mc(:, nm) = r0(:, i) + r0(:, j)
                     end if
                  end do
               end if
            else if (d == 0) then
               do i = 1, ng
                  if (zk(i) == 0) exit
               end do
               r1 = r0(:, i)
               nm = nm + 1
               mc(:, nm) = r1
               r2 = matmul(rotation_mat(z, 2.0_dp*PI/real(n, dp)), r1)
               nm = nm + 1
               if (sum(abs(r1 + r2)) > s%tol) then
                  mc(:, nm) = r1 + r2
               else
                  mc(:, nm) = r2 - r1
               end if
            end if
         end do
         if (nm > 0) then
            allocate (cand(3, nm))
            ncand = 0
            do i = 1, nm
               nrm = sqrt(sum(mc(:, i)**2))
               if (nrm > s%tol) then
                  ncand = ncand + 1
                  cand(:, ncand) = normalize(mc(:, i))
               end if
            end do
            nfin = ncand
            call remove_dupvec(cand, nfin, s%tol)
            do i = 1, nfin
               want = .not. parallel_vectors(cand(:, i), z, s%tol)
               if (want) want = has_rotation(s, cand(:, i), 2)
               if (want) then
                  found = .true.
                  c2x = cand(:, i)
                  deallocate (cand, zk, uk, mc)
                  return
               end if
            end do
            deallocate (cand)
         end if
         deallocate (zk, uk, mc)
      end do
   end subroutine search_c2x

   subroutine search_mirrorx(s, n, found, mirror, zaxis)
      !! pyscf `search_mirrorx` (l.~780): the normal of a mirror parallel to z
      !! (`n > 1`, `zaxis` present), or of any mirror (`n = 1`).
      type(sys_t), intent(in) :: s
      integer, intent(in) :: n
      logical, intent(out) :: found
      real(dp), intent(out) :: mirror(3)
      real(dp), intent(in), optional :: zaxis(3)
      real(dp), allocatable :: r0(:, :)
      real(dp) :: r00(3), r1(3), m(3)
      integer :: g, ng, i
      found = .false.
      mirror = 0.0_dp
      if (n > 1) then
         do g = 1, s%ngrp
            ng = s%gofs(g + 1) - s%gofs(g)
            r00 = s%r(:, s%gidx(s%gofs(g) + 1))
            if (ng > 1 .and. .not. parallel_vectors(r00, zaxis, s%tol)) then
               r1 = matmul(rotation_mat(zaxis, 2.0_dp*PI/real(n, dp)), r00)
               m = normalize(r1 - r00)
               if (has_mirror(s, m)) then
                  found = .true.
                  mirror = m
                  return
               end if
            end if
         end do
      else
         do g = 1, s%ngrp
            ng = s%gofs(g + 1) - s%gofs(g)
            if (ng <= 1) cycle
            r0 = group_coords(s, g)
            do i = 2, ng
               m = normalize(r0(:, i) - r0(:, 1))
               if (has_mirror(s, m)) then
                  found = .true.
                  mirror = m
                  return
               end if
            end do
         end do
      end if
   end subroutine search_mirrorx

   pure function make_axes(z, x) result(ax)
      !! pyscf `_make_axes` (l.915): rows x, y, z.
      real(dp), intent(in) :: z(3), x(3)
      real(dp) :: ax(3, 3)
      real(dp) :: y(3), xx(3)
      y = cross(z, x)
      xx = cross(y, z)
      ax(1, :) = normalize(xx)
      ax(2, :) = normalize(y)
      ax(3, :) = normalize(z)
   end function make_axes

   pure function refine(axes) result(ax)
      !! pyscf `_refine` (l.920): signs of z, x, and y for det = +1.
      real(dp), intent(in) :: axes(3, 3)
      real(dp) :: ax(3, 3)
      integer :: xid, yid
      real(dp) :: det
      ax = axes
      if (ax(3, 3) < 0.0_dp) ax(3, :) = -ax(3, :)
      if (abs(ax(1, 1)) > abs(ax(2, 1))) then
         xid = 1
         yid = 2
      else
         xid = 2
         yid = 1
      end if
      if (ax(xid, 1) < 0.0_dp) ax(xid, :) = -ax(xid, :)
      det = ax(1, 1)*(ax(2, 2)*ax(3, 3) - ax(2, 3)*ax(3, 2)) - ax(1, 2)*(ax(2, 1)*ax(3, 3) - ax(2, 3)*ax(3, 1)) &
            + ax(1, 3)*(ax(2, 1)*ax(3, 2) - ax(2, 2)*ax(3, 1))
      if (det < 0.0_dp) ax(yid, :) = -ax(yid, :)
   end function refine

   pure function closest_axes_eye(ax) result(order)
      !! pyscf `closest_axes(axes, eye)` (l.118): which rows of `ax` become x, y
      !! and z, as close as possible to the input axes. Ties within the
      !! tolerance take the last candidate for z and x.
      real(dp), intent(in) :: ax(3, 3)
      integer :: order(3)
      real(dp), parameter :: TIE = 1.0e-5_dp  ! pyscf symm.TOLERANCE
      real(dp) :: zc(3), xc(3), yc(3)
      integer :: i, x, y, z
      zc = abs(ax(:, 3))
      xc = abs(ax(:, 1))
      yc = abs(ax(:, 2))
      z = 1
      do i = 1, 3
         if (abs(zc(i) - maxval(zc)) < TIE) z = i
      end do
      xc(z) = 0.0_dp
      yc(z) = 0.0_dp
      x = 1
      do i = 1, 3
         if (abs(xc(i) - maxval(xc)) < TIE) x = i
      end do
      yc(x) = 0.0_dp
      y = maxloc(yc, 1)
      order = [x, y, z]
   end function closest_axes_eye

   pure function alias_axes_eye(ax) result(new)
      !! pyscf `alias_axes(axes, eye)` (l.133), including its exchange of x and y
      !! when the renamed axes are left-handed (deterministic here: the axes
      !! come from `_refine`, not from LAPACK).
      real(dp), intent(in) :: ax(3, 3)
      real(dp) :: new(3, 3)
      integer :: o(3)
      real(dp) :: det
      o = closest_axes_eye(ax)
      new = ax(o, :)
      det = new(1, 1)*(new(2, 2)*new(3, 3) - new(2, 3)*new(3, 2)) - new(1, 2)*(new(2, 1)*new(3, 3) - new(2, 3)*new(3, 1)) &
            + new(1, 3)*(new(2, 1)*new(3, 2) - new(2, 2)*new(3, 1))
      if (det < 0.0_dp) new = ax([o(2), o(1), o(3)], :)
   end function alias_axes_eye

   subroutine search_i_group(s, name, ax, found)
      !! pyscf `_search_i_group` (l.824).
      type(sys_t), intent(in) :: s
      character(len=*), intent(out) :: name
      real(dp), intent(out) :: ax(3, 3)
      logical, intent(out) :: found
      real(dp), allocatable :: vk(:, :), c5(:, :)
      integer, allocatable :: nk(:)
      integer :: nv, i, n5
      real(dp) :: za(3), c5b(3), c5a(3)
      found = .false.
      name = ""
      ax = 0.0_dp
      call search_possible_rotations(s, vk=vk, nk=nk, nv=nv)
      allocate (c5(3, max(nv, 1)))
      n5 = 0
      do i = 1, nv
         if (nk(i) == 5) then
            if (has_rotation(s, vk(:, i), 5)) then
               n5 = n5 + 1
               c5(:, n5) = vk(:, i)
            end if
         end if
      end do
      if (n5 <= 1) return
      za = c5(:, 1)
      name = merge("Ih", "I ", has_icenter(s))
      c5b = c5(:, 2)
      if (dot_product(c5b, za) < 0.0_dp) c5b = -c5b
      c5a = matmul(rotation_mat(za, PI*6.0_dp/5.0_dp), c5b)
      ax = make_axes(za, c5a + c5b)
      found = .true.
   end subroutine search_i_group

   subroutine search_ot_group(s, name, ax, found)
      !! pyscf `_search_ot_group` (l.848).
      type(sys_t), intent(in) :: s
      character(len=*), intent(out) :: name
      real(dp), intent(out) :: ax(3, 3)
      logical, intent(out) :: found
      real(dp), allocatable :: vk(:, :), c4(:, :), c3(:, :)
      integer, allocatable :: nk(:)
      integer :: nv, i, n4, n3
      real(dp) :: c3a(3), c3b(3), c3c(3)
      found = .false.
      name = ""
      ax = 0.0_dp
      call search_possible_rotations(s, vk=vk, nk=nk, nv=nv)
      allocate (c4(3, max(nv, 1)), c3(3, max(nv, 1)))
      n4 = 0
      n3 = 0
      do i = 1, nv
         if (nk(i) == 4) then
            if (has_rotation(s, vk(:, i), 4)) then
               n4 = n4 + 1
               c4(:, n4) = vk(:, i)
            end if
         end if
      end do
      if (n4 > 1) then
         name = merge("Oh", "O ", has_icenter(s))
         ax = make_axes(c4(:, 1), c4(:, 2))
         found = .true.
         return
      end if
      if (n4 == 1) then
         return  ! pyscf asserts len(c4_axes) > 1
      end if
      do i = 1, nv
         if (nk(i) == 3) then
            if (has_rotation(s, vk(:, i), 3)) then
               n3 = n3 + 1
               c3(:, n3) = vk(:, i)
            end if
         end if
      end do
      if (n3 <= 1) return
      if (has_icenter(s)) then
         name = "Th"
      else if (has_mirror(s, cross(c3(:, 1), c3(:, 2)))) then
         name = "Td"
      else
         name = "T "
      end if
      c3a = c3(:, 1)
      if (dot_product(c3a, c3(:, 2)) > 0.0_dp) c3a = -c3a
      c3b = matmul(rotation_mat(c3a, -2.0_dp*PI/3.0_dp), c3(:, 2))
      c3c = matmul(rotation_mat(c3a, 2.0_dp*PI/3.0_dp), c3(:, 2))
      ax = make_axes(c3a + c3b, c3a + c3c)
      found = .true.
   end subroutine search_ot_group

   subroutine get_subgroup(gp, n, ax, abel, axa)
      !! pyscf `get_subgroup` (l.347): the Abelian group and axes that
      !! `Mole.groupname` / `_symm_axes` use for the full group `gp` (with the
      !! order `n` of its main axis for the D, C, S families). `abel` is blank
      !! when pyscf's answer is not a group pyscf can use.
      character(len=*), intent(in) :: gp
      integer, intent(in) :: n
      real(dp), intent(in) :: ax(3, 3)
      character(len=3), intent(out) :: abel
      real(dp), intent(out) :: axa(3, 3)
      character(len=1) :: fam, suf
      integer :: lg

      axa = ax
      abel = ""
      select case (trim(gp))
      case ("Oh")
         abel = "D2h"
         axa = alias_axes_eye(ax)
      case ("O", "Td", "T", "Th")
         abel = "D2"
         axa = alias_axes_eye(ax)
      case ("Ih")
         abel = "Ci"
         axa = alias_axes_eye(ax)
      case ("I")
         abel = "C1"
      case default
         lg = len_trim(gp)
         fam = gp(1:1)
         suf = merge(gp(lg:lg), " ", gp(lg:lg) == "h" .or. gp(lg:lg) == "v" .or. gp(lg:lg) == "d")
         if (fam == "S") then
            if (mod(n, 2) == 0 .and. n/2 == 2) abel = "C2"  ! S4; S6 gives 'C3', which pyscf cannot use
            return
         end if
         if (mod(n, 2) == 0) then
            if (fam == "D" .and. suf == "d") then
               abel = "D2"
            else if (fam == "D" .and. suf == "h") then
               abel = "D2h"
            else if (fam == "D") then
               abel = "D2"
            else if (suf == "h") then
               abel = "C2h"
            else if (suf == "v") then
               abel = "C2v"
            else
               abel = "C2"
            end if
         else
            if (fam == "D" .and. suf == "h") then
               abel = "C2v"
               axa = ax([2, 3, 1], :)
            else if (fam == "D" .and. suf == "d") then
               abel = "C2h"
               axa = ax([2, 3, 1], :)
            else if (fam == "D") then
               abel = "C2"
               axa = ax([2, 3, 1], :)
            else if (suf == "h") then
               abel = "Cs"
            else if (suf == "v") then
               abel = "Cs"
               axa = ax([2, 3, 1], :)
            else
               abel = "C1"
            end if
         end if
      end select
   end subroutine get_subgroup

   function cn_name(fam, n, suf) result(name)
      character(len=1), intent(in) :: fam
      integer, intent(in) :: n
      character(len=1), intent(in) :: suf
      character(len=4) :: name
      character(len=2) :: ns
      write (ns, "(i0)") n
      name = trim(fam//trim(ns)//trim(suf))
   end function cn_name

   subroutine sym_top_axes(d, cls, tol, ev, vec, eq12, eq23, kind, topname, abel, ax)
      !! The degenerate branch of pyscf `detect_symm` (l.188-346) followed by
      !! `get_subgroup`. `d` are the atoms about the charge centre, `cls` the
      !! atom classes, `tol` pyscf's `TOLERANCE`, `ev`/`vec` the eigenvalues
      !! (ascending) and eigenvectors (columns) of the second-moment tensor,
      !! `eq12`/`eq23` whether the moments 1,2 and 2,3 are equal.
      !!
      !! `kind` = 1: the full group `topname` was found; `abel` and `ax` (rows
      !! x, y, z) are pyscf's `groupname` and `_symm_axes`. `kind` = 2: no axis
      !! of order 2 or more (accidental degeneracy): `ax` are the axes pyscf then
      !! tests for D2h and its subgroups (the generic path of `trc_sym_geom`
      !! continues with them). `kind` = 0: pyscf's group is not usable.
      real(dp), intent(in) :: d(:, :)
      integer, intent(in) :: cls(:)
      real(dp), intent(in) :: tol
      real(dp), intent(in) :: ev(3), vec(3, 3)
      logical, intent(in) :: eq12, eq23
      integer, intent(out) :: kind
      character(len=4), intent(out) :: topname
      character(len=3), intent(out) :: abel
      real(dp), intent(out) :: ax(3, 3)

      type(sys_t) :: s
      real(dp) :: zhat(3), c2x(3), mir(3), e(3), axes(3, 3), z(3), proj(3), xk(3)
      real(dp) :: tolp
      integer :: natm, n, i, k, m
      logical :: have_c2x, have_mir, found, is_top
      character(len=1) :: suf

      kind = 0
      topname = ""
      abel = ""
      ax = 0.0_dp
      natm = size(d, 2)
      tolp = tol/sqrt(1.0_dp + real(natm, dp))
      call build_sys(d, cls, tol, s)
      do i = 1, 3
         ax(i, i) = 1.0_dp
      end do

      if (ev(3) <= tolp) then  ! l.201: allclose(w1, 0): a single atom
         kind = 1
         topname = "SO3"
         abel = "D2h"
         return
      end if
      if (ev(2) <= tolp) then  ! l.206: linear
         ! z along the molecule, largest component positive
         z = vec(:, 3)
         m = maxloc(abs(z), 1)
         do i = 1, 3
            if (abs(z(i)) >= abs(z(m)) - 1.0e-6_dp) then
               m = i
               exit
            end if
         end do
         if (z(m) < 0.0_dp) z = -z
         k = minloc(abs(z), 1)
         e = 0.0_dp
         e(k) = 1.0_dp
         proj = e - z(k)*z
         xk = normalize(proj)
         ax(1, :) = xk
         ax(2, :) = cross(z, xk)
         ax(3, :) = z
         kind = 1
         if (has_icenter(s)) then
            topname = "Dooh"
            abel = "D2h"
         else
            topname = "Coov"
            abel = "C2v"
         end if
         return
      end if

      n = 0
      have_c2x = .false.
      have_mir = .false.
      is_top = .false.
      zhat = [0.0_dp, 0.0_dp, 1.0_dp]
      if (eq12 .and. eq23) then
         ! spherical top (l.224-246): I, then T/O
         call search_i_group(s, topname, axes, found)
         if (.not. found) call search_ot_group(s, topname, axes, found)
         if (found) then
            axes = refine(axes)
            call get_subgroup(topname, 0, axes, abel, ax)
            kind = merge(1, 0, len_trim(abel) > 0)
            return
         end if
         n = -2  ! fall to the generic search (n is None)
      else
         ! symmetric top (l.248-256): the unique axis, sign fixed
         if (eq12) then
            zhat = vec(:, 3)
         else
            zhat = vec(:, 1)
         end if
         m = maxloc(abs(zhat), 1)
         do i = 1, 3
            if (abs(zhat(i)) >= abs(zhat(m)) - 1.0e-6_dp) then
               m = i
               exit
            end if
         end do
         if (zhat(m) < 0.0_dp) zhat = -zhat
         call search_c_highest(s, e, n, zhat)
         if (n == 1) then
            n = -2
         else
            call search_c2x(s, zhat, n, have_c2x, c2x)
            call search_mirrorx(s, n, have_mir, mir, zhat)
         end if
      end if

      if (n == -2) then  ! n is None (l.273)
         call search_c_highest(s, zhat, n)
         if (n > 1) then
            call search_c2x(s, zhat, n, have_c2x, c2x)
            call search_mirrorx(s, n, have_mir, mir, zhat)
         else
            call search_mirrorx(s, 1, have_mir, mir)
            if (have_mir) then
               axes = make_axes(mir, [1.0_dp, 0.0_dp, 0.0_dp])
               if (parallel_vectors(mir, [1.0_dp, 0.0_dp, 0.0_dp], tol)) axes = make_axes(mir, [0.0_dp, 1.0_dp, 0.0_dp])
            else
               axes = 0.0_dp
               do i = 1, 3
                  axes(i, i) = 1.0_dp
               end do
            end if
            ax = axes
            kind = 2
            return
         end if
      end if

      ! n >= 2 (l.298-318)
      zhat = normalize(zhat)
      if (have_c2x) then
         axes = make_axes(zhat, c2x)
      else if (have_mir) then
         axes = make_axes(zhat, mir)
      else
         axes = 0.0_dp
         do i = 1, 3
            e = 0.0_dp
            e(i) = 1.0_dp
            if (.not. parallel_vectors(e, zhat, tol)) then
               axes = make_axes(zhat, e)
               exit
            end if
         end do
      end if
      suf = " "
      if (have_c2x) then
         if (has_mirror(s, axes(3, :))) then
            suf = "h"
         else if (has_improper_rotation(s, axes(3, :), n)) then
            suf = "d"
         end if
         topname = cn_name("D", n, suf)
      else if (have_mir) then
         topname = cn_name("C", n, "v")
      else if (has_mirror(s, axes(3, :))) then
         topname = cn_name("C", n, "h")
      else if (has_improper_rotation(s, axes(3, :), n)) then
         topname = cn_name("S", 2*n, " ")
      else
         topname = cn_name("C", n, " ")
      end if
      axes = refine(axes)
      if (topname(1:1) == "S") then
         call get_subgroup(topname, 2*n, axes, abel, ax)
      else
         call get_subgroup(topname, n, axes, abel, ax)
      end if
      kind = merge(1, 0, len_trim(abel) > 0)
   end subroutine sym_top_axes

end module trc_sym_tops
