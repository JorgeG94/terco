module trc_sym_frame
   !! The exact AO rotation between the standard frame a symmetric run works
   !! in and the user's input frame it reports in.
   !!
   !! A symmetric run moves the molecule: `x_std = rot (x - shift)`
   !! (`trc_sym_geom_t`). An AO of the run is a function of `x_std`,
   !! centred on `A_std`, with the angular part `p(x_std - A_std)`. Seen from
   !! the input frame the same function is `p(rot (x - A_in))` around
   !! `A_in = shift + rot^T A_std`. A homogeneous polynomial of degree `l`
   !! rotated by `rot` is again one, so, per shell and exactly,
   !!
   !!     p_mu(rot y) = sum_nu  D(nu, mu)  p_nu(y)
   !!
   !! with `p_nu` the AOs of the same shell on the input-frame axes, and `D`
   !! the l-th symmetric tensor power of `rot` on the Cartesian monomials
   !! (every Cartesian component of a terco shell shares the shell's
   !! normalisation, so no normalisation factor enters; `check_sym_frame`
   !! confirms it against the overlap matrix). Orbitals follow as
   !! `C_in = D C_std`, a density as `D_in = D D_std D^T`, and the overlap
   !! satisfies `S_std = D^T S_in D`.
   !!
   !! Spherical shells: the real solid harmonics are `T`-combinations of the
   !! Cartesian components (`SPH_T`). Rotation keeps the harmonic subspace,
   !! so `D T = T X` is solvable exactly, `X = (T^T T)^-1 T^T D T`, and `X`
   !! is the block used in the spherical AO space.
   !!
   !! Only the 3x3 rotation matters, never a particular group; any `l` up to
   !! the build's `TERCO_LMAX` works for Cartesian shells, and up to
   !! `SPH_LMAX` (4) for spherical ones. The shell blocks are small dense
   !! matrices applied shell by shell, so the cost is O(nao^2) per call.
   use pic_types, only: dp
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_basis_types, only: trc_basis_t
   use trc_cart, only: cart_ncart_of
   use trc_cart_perm, only: cart_perm_powers
   use trc_sph_tables, only: SPH_T, SPH_LMAX
   implicit none
   private

   public :: sym_frame_rot_blocks, sym_frame_rotate_rows, sym_frame_rotate_cols, sym_frame_rotate_dmat
   public :: sym_frame_to_input, sym_frame_to_standard

   integer, parameter :: LMAX_BLOCK = 8  !! the dense blocks are sized for l <= 8 (45 Cartesian components)
   integer, parameter :: NCMAX = (LMAX_BLOCK + 1)*(LMAX_BLOCK + 2)/2

contains

   pure function sym_frame_to_input(rot, shift, x) result(y)
      !! Standard-frame point(s) `x(3, n)` in the input frame: `shift + rot^T x`.
      real(dp), intent(in) :: rot(3, 3), shift(3), x(:, :)
      real(dp) :: y(3, size(x, 2))
      integer :: i
      do i = 1, size(x, 2)
         y(:, i) = shift + matmul(transpose(rot), x(:, i))
      end do
   end function sym_frame_to_input

   pure function sym_frame_to_standard(rot, shift, x) result(y)
      !! Input-frame point(s) in the standard frame: `rot (x - shift)`.
      real(dp), intent(in) :: rot(3, 3), shift(3), x(:, :)
      real(dp) :: y(3, size(x, 2))
      integer :: i
      do i = 1, size(x, 2)
         y(:, i) = matmul(rot, x(:, i) - shift)
      end do
   end function sym_frame_to_standard

   subroutine cart_block(l, rot, d)
      !! D(nu, mu) of degree `l`: p_mu(rot y) = sum_nu D(nu, mu) p_nu(y), in
      !! this build's Cartesian component order.
      integer, intent(in) :: l
      real(dp), intent(in) :: rot(3, 3)
      real(dp), intent(out) :: d(:, :)
      real(dp), allocatable :: poly(:, :, :), nxt(:, :, :)
      integer, allocatable :: slot(:, :, :)
      integer :: nc, c, px, py, pz, i, k, qx, qy, qz
      integer :: q(3)

      nc = cart_ncart_of(l)
      allocate (slot(0:l, 0:l, 0:l), poly(0:l, 0:l, 0:l), nxt(0:l, 0:l, 0:l))
      slot = 0
      do c = 1, nc
         call cart_perm_powers(l, c - 1, px, py, pz)
         slot(px, py, pz) = c
      end do
      d(1:nc, 1:nc) = 0.0_dp
      do c = 1, nc
         call cart_perm_powers(l, c - 1, px, py, pz)
         q = [px, py, pz]
         poly = 0.0_dp
         poly(0, 0, 0) = 1.0_dp
         ! prod_i (sum_j rot(i, j) y_j)^q_i, one linear factor at a time
         do i = 1, 3
            do k = 1, q(i)
               nxt = 0.0_dp
               do qz = 0, l
                  do qy = 0, l
                     do qx = 0, l
                        if (poly(qx, qy, qz) == 0.0_dp) cycle
                        if (qx < l) nxt(qx + 1, qy, qz) = nxt(qx + 1, qy, qz) + rot(i, 1)*poly(qx, qy, qz)
                        if (qy < l) nxt(qx, qy + 1, qz) = nxt(qx, qy + 1, qz) + rot(i, 2)*poly(qx, qy, qz)
                        if (qz < l) nxt(qx, qy, qz + 1) = nxt(qx, qy, qz + 1) + rot(i, 3)*poly(qx, qy, qz)
                     end do
                  end do
               end do
               poly = nxt
            end do
         end do
         do qz = 0, l
            do qy = 0, l
               do qx = 0, l
                  if (slot(qx, qy, qz) > 0) d(slot(qx, qy, qz), c) = poly(qx, qy, qz)
               end do
            end do
         end do
      end do
   end subroutine cart_block

   subroutine sph_block(l, dc, x, error)
      !! The spherical block X = (T^T T)^-1 T^T D T from the Cartesian block.
      integer, intent(in) :: l
      real(dp), intent(in) :: dc(:, :)
      real(dp), intent(out) :: x(:, :)
      type(error_t), intent(inout) :: error
      real(dp) :: t(NCMAX, 2*SPH_LMAX + 1), g(2*SPH_LMAX + 1, 2*SPH_LMAX + 1), rhs(2*SPH_LMAX + 1, 2*SPH_LMAX + 1)
      real(dp) :: piv, f, res
      integer :: nc, ns, i, j, k, p

      if (l > SPH_LMAX) then
         call error%set(ERROR_VALIDATION, "sym_frame_rot_blocks: a spherical shell beyond the transform tables")
         return
      end if
      nc = cart_ncart_of(l)
      ns = 2*l + 1
      t(1:nc, 1:ns) = SPH_T(1:nc, 1:ns, l)
      g(1:ns, 1:ns) = matmul(transpose(t(1:nc, 1:ns)), t(1:nc, 1:ns))
      rhs(1:ns, 1:ns) = matmul(transpose(t(1:nc, 1:ns)), matmul(dc(1:nc, 1:nc), t(1:nc, 1:ns)))
      ! Gauss-Jordan with partial pivoting on the ns x ns Gram system
      do k = 1, ns
         p = k
         do i = k + 1, ns
            if (abs(g(i, k)) > abs(g(p, k))) p = i
         end do
         if (abs(g(p, k)) < 1.0e-12_dp) then
            call error%set(ERROR_VALIDATION, "sym_frame_rot_blocks: singular spherical transform block")
            return
         end if
         if (p /= k) then
            do j = 1, ns
               piv = g(k, j)
               g(k, j) = g(p, j)
               g(p, j) = piv
               piv = rhs(k, j)
               rhs(k, j) = rhs(p, j)
               rhs(p, j) = piv
            end do
         end if
         piv = g(k, k)
         g(k, 1:ns) = g(k, 1:ns)/piv
         rhs(k, 1:ns) = rhs(k, 1:ns)/piv
         do i = 1, ns
            if (i == k) cycle
            f = g(i, k)
            g(i, 1:ns) = g(i, 1:ns) - f*g(k, 1:ns)
            rhs(i, 1:ns) = rhs(i, 1:ns) - f*rhs(k, 1:ns)
         end do
      end do
      x(1:ns, 1:ns) = rhs(1:ns, 1:ns)
      ! D T must lie in the span of T: the residual says whether it does.
      res = maxval(abs(matmul(dc(1:nc, 1:nc), t(1:nc, 1:ns)) - matmul(t(1:nc, 1:ns), x(1:ns, 1:ns))))
      if (res > 1.0e-10_dp) then
         call error%set(ERROR_VALIDATION, "sym_frame_rot_blocks: a rotated spherical harmonic left its shell")
      end if
   end subroutine sph_block

   subroutine sym_frame_rot_blocks(rot, spherical, lmax, dl, error)
      !! The per-shell blocks `dl(:, :, l)` (leading `ncomp x ncomp`, with
      !! `ncomp` the Cartesian or the spherical count) for `l = 0..lmax`.
      real(dp), intent(in) :: rot(3, 3)
      logical, intent(in) :: spherical
      integer, intent(in) :: lmax
      real(dp), allocatable, intent(out) :: dl(:, :, :)
      type(error_t), intent(out) :: error
      real(dp) :: dc(NCMAX, NCMAX)
      integer :: l, nc, nb

      if (lmax > LMAX_BLOCK) then
         call error%set(ERROR_VALIDATION, "sym_frame_rot_blocks: angular momentum beyond "//char(48 + LMAX_BLOCK))
         return
      end if
      nb = merge(2*lmax + 1, cart_ncart_of(lmax), spherical)
      allocate (dl(nb, nb, 0:lmax))
      dl = 0.0_dp
      do l = 0, lmax
         nc = cart_ncart_of(l)
         call cart_block(l, rot, dc)
         if (spherical) then
            call sph_block(l, dc, dl(:, :, l), error)
            if (error%has_error()) return
         else
            dl(1:nc, 1:nc, l) = dc(1:nc, 1:nc)
         end if
      end do
   end subroutine sym_frame_rot_blocks

   pure integer function first_ao(bas, s) result(a0)
      type(trc_basis_t), intent(in) :: bas
      integer, intent(in) :: s
      a0 = merge(bas%sh_ao_sph(s), bas%sh_ao(s), bas%spherical)
   end function first_ao

   pure integer function ncomp_of(bas, s) result(nc)
      type(trc_basis_t), intent(in) :: bas
      integer, intent(in) :: s
      nc = merge(2*bas%sh_l(s) + 1, cart_ncart_of(bas%sh_l(s)), bas%spherical)
   end function ncomp_of

   subroutine sym_frame_rotate_rows(bas, dl, a)
      !! `a <- D a` for `a(nao, m)` in the AO space of `bas`: the orbitals
      !! `C_in = D C_std`, one shell block of rows at a time.
      type(trc_basis_t), intent(in) :: bas
      real(dp), intent(in) :: dl(:, :, 0:)
      real(dp), intent(inout) :: a(:, :)
      integer :: s, a0, nc
      do s = 1, bas%nshell
         a0 = first_ao(bas, s)
         nc = ncomp_of(bas, s)
         a(a0:a0 + nc - 1, :) = matmul(dl(1:nc, 1:nc, bas%sh_l(s)), a(a0:a0 + nc - 1, :))
      end do
   end subroutine sym_frame_rotate_rows

   subroutine sym_frame_rotate_cols(bas, dl, a)
      !! `a <- a D^T` for `a(m, nao)`.
      type(trc_basis_t), intent(in) :: bas
      real(dp), intent(in) :: dl(:, :, 0:)
      real(dp), intent(inout) :: a(:, :)
      integer :: s, a0, nc
      do s = 1, bas%nshell
         a0 = first_ao(bas, s)
         nc = ncomp_of(bas, s)
         a(:, a0:a0 + nc - 1) = matmul(a(:, a0:a0 + nc - 1), transpose(dl(1:nc, 1:nc, bas%sh_l(s))))
      end do
   end subroutine sym_frame_rotate_cols

   subroutine sym_frame_rotate_dmat(bas, dl, a)
      !! `a <- D a D^T` for a square `a(nao, nao)`, standard frame to input
      !! frame: right for a contravariant matrix (a density, which pairs with
      !! AO functions as D_uv chi_u chi_v). Not for a covariant one (Fock,
      !! overlap: integrals over AOs) in a Cartesian basis, where D is not
      !! orthogonal and such a matrix goes as `D^-T a D^-1` (module header:
      !! `S_std = D^T S_in D`).
      type(trc_basis_t), intent(in) :: bas
      real(dp), intent(in) :: dl(:, :, 0:)
      real(dp), intent(inout) :: a(:, :)
      call sym_frame_rotate_rows(bas, dl, a)
      call sym_frame_rotate_cols(bas, dl, a)
   end subroutine sym_frame_rotate_dmat

end module trc_sym_frame
