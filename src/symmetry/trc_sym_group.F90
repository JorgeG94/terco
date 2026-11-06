module trc_sym_group
   !! The eight Abelian point groups of the D2h family as data.
   !!
   !! In the standard orientation every operation of D2h is a diagonal sign
   !! matrix diag(sx, sy, sz), so an operation is a 3-bit mask: bit 0 set
   !! flips x, bit 1 flips y, bit 2 flips z, and composition is `ieor`.
   !! E = 0, C2x = 6 (flips y, z), C2y = 5, C2z = 3, i = 7, sx = 1 (the
   !! mirror perpendicular to x flips x), sy = 2, sz = 4.
   !!
   !! The order of operations in each group, and the order and names of its
   !! irreps, are pyscf's: `pyscf/symm/param.py` of pyscf 2.14.0
   !! (`OPERATOR_TABLE`, `IRREP_ID_TABLE`, `CHARACTER_TABLE`), so labels
   !! compare directly with pyscf anchors. Every irrep is one-dimensional
   !! and its characters are +1 or -1.
   !!
   !! Pure data and logic: nothing here touches the basis or the SCF.
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: trc_sym_group_t, trc_sym_group
   public :: SYM_MAX_ORDER, SYM_NGROUPS, sym_group_names

   integer, parameter :: SYM_MAX_ORDER = 8
   integer, parameter :: SYM_NGROUPS = 8

   !> Group names, in pyscf's `param.POINTGROUP` order.
   character(len=3), parameter :: sym_group_names(SYM_NGROUPS) = &
                                  [character(len=3) :: "D2h", "C2h", "C2v", "D2", "Cs", "Ci", "C2", "C1"]

   type :: trc_sym_group_t
      character(len=3) :: name = ""
      !> Number of operations, equal to the number of irreps.
      integer :: order = 0
      !> Sign masks of the operations (see the module header), pyscf's order.
      integer, allocatable :: ops(:)
      !> Operation names ("E", "C2z", "sx", ...), parallel to `ops`.
      character(len=3), allocatable :: op_names(:)
      !> chars(irrep, op): +1 or -1.
      integer, allocatable :: chars(:, :)
      character(len=3), allocatable :: irrep_names(:)
      !> prod(i, j): the irrep of the direct product of irreps i and j.
      integer, allocatable :: prod(:, :)
   contains
      !> Index in `ops` of the operation with this mask, or 0.
      procedure :: op_index => group_op_index
   end type trc_sym_group_t

contains

   pure function group_op_index(this, mask) result(k)
      class(trc_sym_group_t), intent(in) :: this
      integer, intent(in) :: mask
      integer :: k
      do k = 1, this%order
         if (this%ops(k) == mask) return
      end do
      k = 0
   end function group_op_index

   pure function lower(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i, c
      t = s
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) t(i:i) = achar(c + 32)
      end do
   end function lower

   function trc_sym_group(name, error) result(g)
      !! The group called `name` (case-insensitive: "C2v", "c2v"). An unknown
      !! name sets `error` and returns an empty group.
      character(len=*), intent(in) :: name
      type(error_t), intent(out) :: error
      type(trc_sym_group_t) :: g
      integer :: i, j, k, n, ig
      integer, allocatable :: ch(:, :)

      ig = 0
      do i = 1, SYM_NGROUPS
         if (lower(trim(adjustl(name))) == lower(trim(sym_group_names(i)))) ig = i
      end do
      if (ig == 0) then
         call error%set(ERROR_VALIDATION, "unknown point group '"//trim(name)// &
                        "' (expected D2h, C2h, C2v, D2, Cs, Ci, C2 or C1)")
         return
      end if

      ! pyscf CHARACTER_TABLE rows, irrep by irrep, columns in OPERATOR_TABLE order.
      select case (ig)
      case (1) ! D2h
         g%ops = [0, 6, 5, 3, 7, 1, 2, 4]
         g%op_names = [character(len=3) :: "E", "C2x", "C2y", "C2z", "i", "sx", "sy", "sz"]
         g%irrep_names = [character(len=3) :: "Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"]
         ch = reshape([ &
              1, 1, 1, 1, 1, 1, 1, 1, &
              1, -1, -1, 1, 1, -1, -1, 1, &
              1, -1, 1, -1, 1, -1, 1, -1, &
              1, 1, -1, -1, 1, 1, -1, -1, &
              1, 1, 1, 1, -1, -1, -1, -1, &
              1, -1, -1, 1, -1, 1, 1, -1, &
              1, -1, 1, -1, -1, 1, -1, 1, &
              1, 1, -1, -1, -1, -1, 1, 1], [8, 8], order=[2, 1])
      case (2) ! C2h
         g%ops = [0, 3, 7, 4]
         g%op_names = [character(len=3) :: "E", "C2z", "i", "sz"]
         g%irrep_names = [character(len=3) :: "Ag", "Bg", "Au", "Bu"]
         ch = reshape([ &
              1, 1, 1, 1, &
              1, -1, 1, -1, &
              1, 1, -1, -1, &
              1, -1, -1, 1], [4, 4], order=[2, 1])
      case (3) ! C2v
         g%ops = [0, 3, 1, 2]
         g%op_names = [character(len=3) :: "E", "C2z", "sx", "sy"]
         g%irrep_names = [character(len=3) :: "A1", "A2", "B1", "B2"]
         ch = reshape([ &
              1, 1, 1, 1, &
              1, 1, -1, -1, &
              1, -1, -1, 1, &
              1, -1, 1, -1], [4, 4], order=[2, 1])
      case (4) ! D2
         g%ops = [0, 6, 5, 3]
         g%op_names = [character(len=3) :: "E", "C2x", "C2y", "C2z"]
         g%irrep_names = [character(len=3) :: "A", "B1", "B2", "B3"]
         ch = reshape([ &
              1, 1, 1, 1, &
              1, -1, -1, 1, &
              1, -1, 1, -1, &
              1, 1, -1, -1], [4, 4], order=[2, 1])
      case (5) ! Cs
         g%ops = [0, 4]
         g%op_names = [character(len=3) :: "E", "sz"]
         g%irrep_names = [character(len=3) :: "A'", 'A"']
         ch = reshape([1, 1, 1, -1], [2, 2], order=[2, 1])
      case (6) ! Ci
         g%ops = [0, 7]
         g%op_names = [character(len=3) :: "E", "i"]
         g%irrep_names = [character(len=3) :: "Ag", "Au"]
         ch = reshape([1, 1, 1, -1], [2, 2], order=[2, 1])
      case (7) ! C2
         g%ops = [0, 3]
         g%op_names = [character(len=3) :: "E", "C2z"]
         g%irrep_names = [character(len=3) :: "A", "B"]
         ch = reshape([1, 1, 1, -1], [2, 2], order=[2, 1])
      case default ! C1
         g%ops = [0]
         g%op_names = [character(len=3) :: "E"]
         g%irrep_names = [character(len=3) :: "A"]
         ch = reshape([1], [1, 1])
      end select

      g%name = sym_group_names(ig)
      n = size(g%ops)
      g%order = n
      g%chars = ch

      ! Direct product: the irrep whose characters are the elementwise product.
      allocate (g%prod(n, n))
      g%prod = 0
      do i = 1, n
         do j = 1, n
            do k = 1, n
               if (all(g%chars(k, :) == g%chars(i, :)*g%chars(j, :))) g%prod(i, j) = k
            end do
         end do
      end do
   end function trc_sym_group

end module trc_sym_group
