module trc_parity
   !! Parity skipping, host side (docs/SYMMETRY_PARITY_SKIP.md): which axis of
   !! the molecule is the dead one, which shells lie in its plane, and the
   !! axis permutation inside the integral engine that makes it the z axis of
   !! the kernels.
   !!
   !! A coordinate d is DEAD for a quartet when all four centres have it
   !! equal to zero. Every component integral of odd combined parity in a
   !! dead coordinate is then exactly zero, and a masked kernel variant leaves
   !! those out. The rule is tied to the detected group (decision 3 of the
   !! design): axis d is a candidate only if the group contains the mirror
   !! sigma_d (the operation that flips d alone), and a shell is in its plane
   !! when its centre has coordinate d EXACTLY zero, which is what `sym_detect`
   !! writes for an atom on a symmetry element. With the symmetry off, or the
   !! group C1, there is no group and nothing is masked.
   !!
   !! One axis only (decision 1): the one with the most shells in its plane,
   !! ties to z. The engine's frame then reads the molecule's axes in a
   !! cyclic order, kernel axis k being molecule axis axperm(k), with
   !! axperm(3) the dead axis.
   use pic_types, only: dp
   use trc_sym_group, only: trc_sym_group_t
   use trc_cart_perm, only: cart_perm_powers
   implicit none
   private

   public :: parity_choose_axis, parity_axperm, parity_ao_perm, parity_mode_env, parity_set_default
   public :: PARITY_OFF, PARITY_AUTO, PARITY_FORCE_FRAME

   integer, parameter :: PARITY_OFF = 0
   !! No dead axis is looked for; the engine is what it was without parity skipping.
   integer, parameter :: PARITY_AUTO = 1
   !! The default: masks and the frame are used when a masked kernel exists.
   integer, parameter :: PARITY_FORCE_FRAME = 2
   !! As AUTO, and the axis permutation is applied even when no masked kernel
   !! exists (the tests of the permutation round trip).

   integer :: default_mode = -1
   !! Process-wide default set by `parity_set_default`; -1: the environment.

contains

   subroutine parity_set_default(mode)
      !! Override the mode of every engine built after this that did not set
      !! its own (the SCF builds its engine itself): -1 returns to TRC_PARITY.
      integer, intent(in) :: mode
      default_mode = mode
   end subroutine parity_set_default

   subroutine parity_choose_axis(group, perm, sh_at, sh_r, axis, shmask)
      !! The dead axis (1 to 3, 0 for none) and `shmask(s)`, 1 for a shell
      !! whose centre the mirror sigma_axis fixes. The axis qualifies only if
      !! EVERY fixed centre has that coordinate exactly zero: a plane whose
      !! atoms carry rounding noise (a molecule rotated into its frame by a
      !! tops transformation) is no dead plane for the kernels, and a mask
      !! that took some of its atoms and not others would not be invariant
      !! under the group.
      type(trc_sym_group_t), intent(in) :: group
      integer, intent(in) :: perm(:, :)  !! (natm, nop) atom each operation sends an atom to (trc_sym_geom)
      integer, intent(in) :: sh_at(:)  !! (nshell) atom of each shell
      real(dp), intent(in) :: sh_r(:, :)  !! (3, nshell) shell centres, standard frame
      integer, intent(out) :: axis
      integer, allocatable, intent(out) :: shmask(:)
      integer, allocatable :: m(:)
      integer :: d, n, best, k, s
      integer :: nshell
      logical :: exact
      nshell = size(sh_r, 2)
      allocate (shmask(nshell), m(nshell))
      shmask = 0
      axis = 0
      best = 0
      do d = 3, 1, -1
         k = group%op_index(ishft(1, d - 1))
         if (k == 0) cycle
         exact = .true.
         do s = 1, nshell
            m(s) = merge(1, 0, perm(sh_at(s), k) == sh_at(s))
            if (m(s) == 1 .and. sh_r(d, s) /= 0.0_dp) exact = .false.
         end do
         if (.not. exact) cycle
         n = sum(m)
         if (n > best) then
            best = n
            axis = d
            shmask = m
         end if
      end do
   end subroutine parity_choose_axis

   pure function parity_axperm(axis) result(p)
      !! Kernel axis k is molecule axis p(k); p(3) = axis, the order cyclic
      !! ((x, y, z) for z, (y, z, x) for x, (z, x, y) for y).
      integer, intent(in) :: axis
      integer :: p(3)
      integer :: k
      do k = 1, 3
         p(k) = mod(axis + k - 1, 3) + 1
      end do
   end function parity_axperm

   subroutine parity_ao_perm(nshell, sh_l, sh_ao, nao, axperm, perm, inv)
      !! The AO map of the axis permutation, host code. `perm(i)` is the AO
      !! of the molecule's frame that AO i of the kernels' frame is: the
      !! component with powers (a, b, c) along the kernel axes has the powers
      !! along the molecule axes with a at axperm(1), b at axperm(2), c at
      !! axperm(3). Found with the build's own component order
      !! (cart_perm_powers), so it holds for any TRC_CART_ORDER. `inv` is the
      !! inverse permutation.
      integer, intent(in) :: nshell, nao
      integer, intent(in) :: sh_l(nshell), sh_ao(nshell), axperm(3)
      integer, intent(out) :: perm(nao), inv(nao)
      integer :: pk(3), po(3)
      integer :: s, ic, jc, ncart, base, qx, qy, qz, k
      do s = 1, nshell
         base = sh_ao(s) - 1
         ncart = (sh_l(s) + 1)*(sh_l(s) + 2)/2
         do ic = 0, ncart - 1
            call cart_perm_powers(sh_l(s), ic, pk(1), pk(2), pk(3))
            do k = 1, 3
               po(axperm(k)) = pk(k)
            end do
            do jc = 0, ncart - 1
               call cart_perm_powers(sh_l(s), jc, qx, qy, qz)
               if (qx == po(1) .and. qy == po(2) .and. qz == po(3)) exit
            end do
            perm(base + ic + 1) = base + jc + 1
            inv(base + jc + 1) = base + ic + 1
         end do
      end do
   end subroutine parity_ao_perm

   function parity_mode_env() result(mode)
      !! TRC_PARITY in the environment: 0 off, 1 auto (the default), 2 force
      !! the frame permutation.
      integer :: mode, ios, v
      character(len=8) :: e
      mode = PARITY_AUTO
      if (default_mode >= 0) then
         mode = default_mode
         return
      end if
      e = " "
      call get_environment_variable("TRC_PARITY", e)
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) v
         if (ios == 0 .and. v >= PARITY_OFF .and. v <= PARITY_FORCE_FRAME) mode = v
      end if
   end function parity_mode_env

end module trc_parity
