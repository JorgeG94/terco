module trc_sym_ao
   !! AO representation of the point-group operations, and the SO basis.
   !!
   !! Every operation of the Abelian groups of `trc_sym_group` is a diagonal
   !! sign matrix diag(sx, sy, sz) in the standard frame of `trc_sym_geom`.
   !! Acting on an AO it sends the function on atom A to a function of the
   !! same shell and component on the atom `perm(A, g)`, times a sign:
   !!   R(g) |A, shell, comp> = ao_sign(.) |perm(A, g), shell, comp>
   !! and that is all `ao_map` and `ao_sign` record, per AO and operation.
   !!
   !!  - Cartesian AOs: the component with powers (px, py, pz) has sign
   !!    sx^px sy^py sz^pz. The powers come from `cart_perm_powers`, the
   !!    generated table of this build, so both `TERCO_CART_ORDER` settings
   !!    give the right signs with no table of their own here.
   !!  - Spherical AOs (`bas%spherical`): column m of the cart-to-sph block
   !!    T(:, m, l) of `trc_sph_tables` is a real solid harmonic. The signed
   !!    Cartesian components sign(c) T(c, m) must be +-1 times the same
   !!    column; that +-1 is the parity, computed here and never tabulated by
   !!    hand. Anything else (a non-+-1 parity, or a harmonic that is not
   !!    mapped to itself) sets `error` instead of guessing.
   !!
   !! The SO basis is the set of projections P_G |a> = (1/|G|) sum_g
   !! chi_G(g) R(g) |a> of the AOs, one orbit of AOs at a time: the
   !! projections of the members of one orbit are proportional, so only the
   !! first AO of each orbit is projected, and a projection that vanishes
   !! (the orbit carries no function of that irrep) adds no column. Columns
   !! are normalized in the Euclidean sense, so `U^T U = I`; they are NOT
   !! orthonormal under the AO overlap, which is why `U^T S U` is a
   !! block-diagonal matrix and not the identity. `u` is stored as per-irrep
   !! column blocks in the irrep order of the group (the pyscf order).
   !!
   !! Equivalent atoms must carry the same shells, in the same order, with
   !! identical angular momentum, exponents and coefficients: the shell map
   !! is positional within an atom. Anything else sets `error`.
   !!
   !! All host code, and the AO maps are in the AO space the SCF runs in:
   !! `nao_sph` functions when `bas%spherical`, else `nao`. The basis must be
   !! built on the standard-frame coordinates (`geom%r_std`); that is checked.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_basis_types, only: trc_basis_t
   use trc_cart_perm, only: cart_perm_powers
   use trc_sph_tables, only: SPH_T, SPH_LMAX, SPH_NCART_MAX
   use trc_sym_group, only: trc_sym_group_t
   use trc_sym_geom, only: trc_sym_geom_t
   implicit none
   private

   public :: trc_sym_ao_t, sym_ao_build

   real(dp), parameter :: PARITY_TOL = 1.0e-12_dp  !! how close a harmonic must be to +-itself
   real(dp), parameter :: SHELL_TOL = 1.0e-13_dp  !! relative tolerance of equal exponents/coefficients
   real(dp), parameter :: FRAME_TOL = 1.0e-8_dp  !! Bohr: basis centres against the standard frame

   type :: trc_sym_ao_t
      type(trc_sym_group_t) :: group  !! the group the maps belong to
      logical :: spherical = .false.  !! AO space is spherical (`nao_sph` functions)
      integer :: nao = 0  !! number of AOs in that space
      integer :: nop = 0  !! number of operations
      integer :: nirrep = 0  !! number of irreps (equal to `nop`)
      integer, allocatable :: ao_map(:, :)  !! ao_map(a, g): the AO that R(g) sends AO a to
      integer, allocatable :: ao_sign(:, :)  !! ao_sign(a, g): +1 or -1, the sign R(g) gives it
      real(dp), allocatable :: u(:, :)  !! (nao, nao) SO basis, columns grouped by irrep
      integer, allocatable :: irrep_off(:)  !! columns of irrep i are irrep_off(i)+1 .. irrep_off(i)+irrep_dim(i)
      integer, allocatable :: irrep_dim(:)  !! number of SOs of each irrep
   end type trc_sym_ao_t

contains

   function sym_ao_build(bas, geom, error) result(ao)
      !! AO maps, signs and the SO basis of `bas` under the group and atom
      !! permutations of `geom`. `bas` must sit on `geom%r_std`.
      type(trc_basis_t), intent(in) :: bas
      type(trc_sym_geom_t), intent(in) :: geom
      type(error_t), intent(out) :: error
      type(trc_sym_ao_t) :: ao

      integer, allocatable :: shmap(:, :)
      integer :: nop

      nop = geom%group%order
      if (nop < 1 .or. .not. allocated(geom%perm)) then
         call error%set(ERROR_VALIDATION, "sym_ao_build: geometry has no group or atom permutation table")
         return
      end if
      if (size(geom%perm, 1) /= bas%natm .or. size(geom%perm, 2) /= nop .or. bas%nshell < 1) then
         call error%set(ERROR_VALIDATION, "sym_ao_build: atom permutation table does not match the basis")
         return
      end if
      ao%group = geom%group
      ao%spherical = bas%spherical
      ao%nao = merge(bas%nao_sph, bas%nao, bas%spherical)
      ao%nop = nop
      ao%nirrep = geom%group%order

      call check_frame(bas, geom, error)
      if (error%has_error()) return
      call shell_map(bas, geom, shmap, error)
      if (error%has_error()) return
      call build_maps(bas, geom, shmap, ao, error)
      if (error%has_error()) return
      call build_so(ao)
   end function sym_ao_build

   subroutine check_frame(bas, geom, error)
      !! The basis centres must be mapped onto each other by the operations
      !! of the group, i.e. the basis was built in the standard frame.
      type(trc_basis_t), intent(in) :: bas
      type(trc_sym_geom_t), intent(in) :: geom
      type(error_t), intent(inout) :: error
      integer :: k, ia, ic
      real(dp) :: s(3), x(3)

      do k = 1, geom%group%order
         do ic = 1, 3
            s(ic) = merge(-1.0_dp, 1.0_dp, btest(geom%group%ops(k), ic - 1))
         end do
         do ia = 1, bas%natm
            x = s*bas%at_r(:, ia)
            if (maxval(abs(x - bas%at_r(:, geom%perm(ia, k)))) > FRAME_TOL) then
               call error%set(ERROR_VALIDATION, "sym_ao_build: the basis centres are not invariant under "// &
                              "the group operations; build the basis on the standard-frame coordinates")
               return
            end if
         end do
      end do
   end subroutine check_frame

   subroutine shell_map(bas, geom, shmap, error)
      !! shmap(s, g): the shell that operation g sends shell s to, the n-th
      !! shell of the atom it sits on going to the n-th shell of the image
      !! atom. Refuses when equivalent atoms do not have identical shells.
      type(trc_basis_t), intent(in) :: bas
      type(trc_sym_geom_t), intent(in) :: geom
      integer, allocatable, intent(out) :: shmap(:, :)
      type(error_t), intent(inout) :: error
      integer, allocatable :: cnt(:), first(:, :), nth(:)
      integer :: s, ia, ja, k, n, maxsh, t

      allocate (cnt(bas%natm), nth(bas%nshell))
      cnt = 0
      do s = 1, bas%nshell
         ia = bas%sh_at(s)
         if (ia < 1 .or. ia > bas%natm) then
            call error%set(ERROR_VALIDATION, "sym_ao_build: a shell sits on no atom, so it has no image")
            return
         end if
         cnt(ia) = cnt(ia) + 1
         nth(s) = cnt(ia)
      end do
      maxsh = maxval(cnt)
      allocate (first(maxsh, bas%natm), shmap(bas%nshell, geom%group%order))
      first = 0
      do s = 1, bas%nshell
         first(nth(s), bas%sh_at(s)) = s
      end do

      do k = 1, geom%group%order
         do s = 1, bas%nshell
            ia = bas%sh_at(s)
            ja = geom%perm(ia, k)
            if (cnt(ja) /= cnt(ia)) then
               call error%set(ERROR_VALIDATION, "sym_ao_build: equivalent atoms "//to_string(ia)//" and "//to_string(ja)// &
                              " have different numbers of shells")
               return
            end if
            t = first(nth(s), ja)
            n = bas%sh_np(s)
            if (bas%sh_l(t) /= bas%sh_l(s) .or. bas%sh_np(t) /= n) then
               call error%set(ERROR_VALIDATION, "sym_ao_build: equivalent atoms "//to_string(ia)//" and "//to_string(ja)// &
                              " have different shell sequences (angular momentum or primitive count)")
               return
            end if
            if (any(abs(bas%sh_e(1:n, t) - bas%sh_e(1:n, s)) > SHELL_TOL*abs(bas%sh_e(1:n, s))) .or. &
                any(abs(bas%sh_c(1:n, t) - bas%sh_c(1:n, s)) > SHELL_TOL*max(abs(bas%sh_c(1:n, s)), 1.0e-300_dp))) then
               call error%set(ERROR_VALIDATION, "sym_ao_build: equivalent atoms "//to_string(ia)//" and "//to_string(ja)// &
                              " have different shell exponents or coefficients")
               return
            end if
            shmap(s, k) = t
         end do
      end do
   end subroutine shell_map

   subroutine build_maps(bas, geom, shmap, ao, error)
      !! Fill ao%ao_map and ao%ao_sign from the shell map and the signs.
      type(trc_basis_t), intent(in) :: bas
      type(trc_sym_geom_t), intent(in) :: geom
      integer, intent(in) :: shmap(:, :)
      type(trc_sym_ao_t), intent(inout) :: ao
      type(error_t), intent(inout) :: error
      integer :: k, s, l, c, ncomp, a0, b0, mask, px, py, pz, sg
      integer :: par(2*SPH_LMAX + 1)

      allocate (ao%ao_map(ao%nao, ao%nop), ao%ao_sign(ao%nao, ao%nop))
      do k = 1, ao%nop
         mask = geom%group%ops(k)
         do s = 1, bas%nshell
            l = bas%sh_l(s)
            if (ao%spherical) then
               if (l > SPH_LMAX) then
                  call error%set(ERROR_VALIDATION, "sym_ao_build: shell beyond the spherical transform tables")
                  return
               end if
               ncomp = 2*l + 1
               a0 = bas%sh_ao_sph(s)
               b0 = bas%sh_ao_sph(shmap(s, k))
               call sph_parities(l, mask, par, error)
               if (error%has_error()) return
               do c = 1, ncomp
                  ao%ao_map(a0 + c - 1, k) = b0 + c - 1
                  ao%ao_sign(a0 + c - 1, k) = par(c)
               end do
            else
               ncomp = (l + 1)*(l + 2)/2
               a0 = bas%sh_ao(s)
               b0 = bas%sh_ao(shmap(s, k))
               do c = 1, ncomp
                  call cart_perm_powers(l, c - 1, px, py, pz)
                  sg = 1
                  if (btest(mask, 0) .and. mod(px, 2) == 1) sg = -sg
                  if (btest(mask, 1) .and. mod(py, 2) == 1) sg = -sg
                  if (btest(mask, 2) .and. mod(pz, 2) == 1) sg = -sg
                  ao%ao_map(a0 + c - 1, k) = b0 + c - 1
                  ao%ao_sign(a0 + c - 1, k) = sg
               end do
            end if
         end do
      end do
   end subroutine build_maps

   subroutine sph_parities(l, mask, par, error)
      !! Parity of each real solid harmonic of angular momentum `l` under
      !! the sign operation `mask`, from the cart-to-sph block: the signed
      !! Cartesian column must be par(m) times the same column.
      integer, intent(in) :: l, mask
      integer, intent(out) :: par(:)
      type(error_t), intent(inout) :: error
      integer :: m, c, nc, px, py, pz, sg
      real(dp) :: v(SPH_NCART_MAX)
      real(dp) :: norm, p, resid

      nc = (l + 1)*(l + 2)/2
      do m = 1, 2*l + 1
         norm = 0.0_dp
         p = 0.0_dp
         do c = 1, nc
            call cart_perm_powers(l, c - 1, px, py, pz)
            sg = 1
            if (btest(mask, 0) .and. mod(px, 2) == 1) sg = -sg
            if (btest(mask, 1) .and. mod(py, 2) == 1) sg = -sg
            if (btest(mask, 2) .and. mod(pz, 2) == 1) sg = -sg
            v(c) = real(sg, dp)*SPH_T(c, m, l)
            norm = norm + SPH_T(c, m, l)**2
            p = p + v(c)*SPH_T(c, m, l)
         end do
         p = p/norm
         resid = maxval(abs(v(1:nc) - p*SPH_T(1:nc, m, l)))
         if (resid > PARITY_TOL .or. abs(abs(p) - 1.0_dp) > PARITY_TOL) then
            call error%set(ERROR_VALIDATION, "sym_ao_build: a spherical harmonic is not mapped to +-itself "// &
                           "by a sign operation (parity is not a pure +-1)")
            return
         end if
         par(m) = nint(p)
      end do
   end subroutine sph_parities

   subroutine build_so(ao)
      !! Project the first AO of each orbit onto each irrep and keep the
      !! non-vanishing, normalized results as columns of `u`.
      type(trc_sym_ao_t), intent(inout) :: ao
      real(dp), allocatable :: v(:)
      integer :: ir, a, k, j, ncol, jmin
      real(dp) :: nrm

      allocate (ao%u(ao%nao, ao%nao), ao%irrep_off(ao%nirrep), ao%irrep_dim(ao%nirrep), v(ao%nao))
      ao%u = 0.0_dp
      ncol = 0
      do ir = 1, ao%nirrep
         ao%irrep_off(ir) = ncol
         do a = 1, ao%nao
            if (a /= minval(ao%ao_map(a, :))) cycle  ! first AO of its orbit only
            v = 0.0_dp
            do k = 1, ao%nop
               j = ao%ao_map(a, k)
               v(j) = v(j) + real(ao%group%chars(ir, k)*ao%ao_sign(a, k), dp)
            end do
            nrm = sqrt(sum(v*v))
            ! Entries are integer sums, so a vanishing projection is exactly 0.
            if (nrm < 0.5_dp) cycle
            jmin = ao%ao_map(a, 1)
            do k = 2, ao%nop
               jmin = min(jmin, ao%ao_map(a, k))
            end do
            ncol = ncol + 1
            ao%u(:, ncol) = v/nrm
            ! sign convention: the lowest-index nonzero entry is positive
            do j = jmin, ao%nao
               if (abs(ao%u(j, ncol)) > 0.0_dp) then
                  if (ao%u(j, ncol) < 0.0_dp) ao%u(:, ncol) = -ao%u(:, ncol)
                  exit
               end if
            end do
         end do
         ao%irrep_dim(ir) = ncol - ao%irrep_off(ir)
      end do
   end subroutine build_so

end module trc_sym_ao
