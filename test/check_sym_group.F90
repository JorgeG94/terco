module check_sym_group_suite
   !! trc_sym_group, the Abelian point-group tables: every group is closed
   !! under mask composition, its characters are orthogonal and as many as
   !! its operations, the operation order and irrep names match pyscf, the
   !! characters are homomorphisms of the operation masks, x, y and z
   !! transform as the irreps pyscf gives them, and the product table
   !! agrees with the characters.
   !!
   !! The expected names and orders below are copied from pyscf 2.14.0,
   !! `pyscf/symm/param.py` (OPERATOR_TABLE, IRREP_ID_TABLE,
   !! CHARACTER_TABLE). Static data, so no .ref file.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_sym_group, only: trc_sym_group_t, trc_sym_group, sym_group_names, SYM_NGROUPS
   implicit none
   private
   public :: collect_sym_group

contains

   subroutine collect_sym_group(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("closed under composition", test_closure), &
                  new_unittest("characters orthogonal", test_orthogonality), &
                  new_unittest("irreps as many as operations", test_irrep_count), &
                  new_unittest("operation order and names match pyscf", test_ops_pyscf), &
                  new_unittest("irrep names and order match pyscf", test_irreps_pyscf), &
                  new_unittest("characters are homomorphisms of the masks", test_homomorphism), &
                  new_unittest("x, y, z irreps match pyscf", test_axes), &
                  new_unittest("product table matches characters", test_product), &
                  new_unittest("unknown name is an error", test_unknown), &
                  new_unittest("name is case-insensitive", test_case) &
                  ]
   end subroutine collect_sym_group

   function get(name) result(g)
      character(len=*), intent(in) :: name
      type(trc_sym_group_t) :: g
      type(error_t) :: err
      g = trc_sym_group(name, err)
      if (err%has_error()) error stop "check_sym_group: "//err%get_message()
   end function get

   subroutine test_closure(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: g
      integer :: ig, i, j
      do ig = 1, SYM_NGROUPS
         g = get(sym_group_names(ig))
         call check(error, g%ops(1), 0, "E first in "//g%name)
         if (allocated(error)) return
         do i = 1, g%order
            do j = 1, g%order
               call check(error, g%op_index(ieor(g%ops(i), g%ops(j))) > 0, &
                          "closure of "//g%name)
               if (allocated(error)) return
            end do
            ! every operation is its own inverse
            call check(error, ieor(g%ops(i), g%ops(i)), 0, "involution in "//g%name)
            if (allocated(error)) return
         end do
         ! all masks distinct
         do i = 1, g%order
            call check(error, count(g%ops == g%ops(i)), 1, "distinct ops in "//g%name)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_closure

   subroutine test_orthogonality(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: g
      integer :: ig, i, j
      do ig = 1, SYM_NGROUPS
         g = get(sym_group_names(ig))
         do i = 1, g%order
            do j = 1, g%order
               call check(error, sum(g%chars(i, :)*g%chars(j, :)), merge(g%order, 0, i == j), &
                          "orthogonality in "//g%name)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_orthogonality

   subroutine test_irrep_count(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: g
      integer :: ig
      integer, parameter :: expected(SYM_NGROUPS) = [8, 4, 4, 4, 2, 2, 2, 1]
      do ig = 1, SYM_NGROUPS
         g = get(sym_group_names(ig))
         call check(error, g%order, expected(ig), "order of "//g%name)
         if (allocated(error)) return
         call check(error, size(g%irrep_names), g%order, "irrep count of "//g%name)
         if (allocated(error)) return
         call check(error, all(shape(g%chars) == g%order), "chars shape of "//g%name)
         if (allocated(error)) return
         call check(error, all(abs(g%chars) == 1), "characters are +-1 in "//g%name)
         if (allocated(error)) return
         ! the totally symmetric irrep comes first
         call check(error, all(g%chars(1, :) == 1), "first irrep symmetric in "//g%name)
         if (allocated(error)) return
         ! E has character 1 in every irrep (one-dimensional)
         call check(error, all(g%chars(:, 1) == 1), "E character 1 in "//g%name)
         if (allocated(error)) return
      end do
   end subroutine test_irrep_count

   subroutine test_ops_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: g
      call expect_ops("D2h", "E C2x C2y C2z i sx sy sz", error)
      if (allocated(error)) return
      call expect_ops("C2h", "E C2z i sz", error)
      if (allocated(error)) return
      call expect_ops("C2v", "E C2z sx sy", error)
      if (allocated(error)) return
      call expect_ops("D2", "E C2x C2y C2z", error)
      if (allocated(error)) return
      call expect_ops("Cs", "E sz", error)
      if (allocated(error)) return
      call expect_ops("Ci", "E i", error)
      if (allocated(error)) return
      call expect_ops("C2", "E C2z", error)
      if (allocated(error)) return
      call expect_ops("C1", "E", error)
      if (allocated(error)) return
      g = get("D2h")
      ! the masks themselves: bit 0 flips x, bit 1 y, bit 2 z
      call check(error, all(g%ops == [0, 6, 5, 3, 7, 1, 2, 4]), "D2h masks")
   end subroutine test_ops_pyscf

   subroutine expect_ops(name, joined, error)
      character(len=*), intent(in) :: name, joined
      type(error_type), allocatable, intent(inout) :: error
      type(trc_sym_group_t) :: g
      character(len=:), allocatable :: got
      integer :: i
      g = get(name)
      got = trim(g%op_names(1))
      do i = 2, g%order
         got = got//" "//trim(g%op_names(i))
      end do
      call check(error, got, joined, "operation order of "//name)
   end subroutine expect_ops

   subroutine test_irreps_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_irreps("D2h", "Ag B1g B2g B3g Au B1u B2u B3u", error)
      if (allocated(error)) return
      call expect_irreps("C2h", "Ag Bg Au Bu", error)
      if (allocated(error)) return
      call expect_irreps("C2v", "A1 A2 B1 B2", error)
      if (allocated(error)) return
      call expect_irreps("D2", "A B1 B2 B3", error)
      if (allocated(error)) return
      call expect_irreps("Cs", "A' A""", error)
      if (allocated(error)) return
      call expect_irreps("Ci", "Ag Au", error)
      if (allocated(error)) return
      call expect_irreps("C2", "A B", error)
      if (allocated(error)) return
      call expect_irreps("C1", "A", error)
   end subroutine test_irreps_pyscf

   subroutine expect_irreps(name, joined, error)
      character(len=*), intent(in) :: name, joined
      type(error_type), allocatable, intent(inout) :: error
      type(trc_sym_group_t) :: g
      character(len=:), allocatable :: got
      integer :: i
      g = get(name)
      got = trim(g%irrep_names(1))
      do i = 2, g%order
         got = got//" "//trim(g%irrep_names(i))
      end do
      call check(error, got, joined, "irrep names of "//name)
   end subroutine expect_irreps

   subroutine test_homomorphism(error)
      !! chi(a)*chi(b) = chi(a*b) for every irrep, with a*b the XOR of masks.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: g
      integer :: ig, k, a, b, c
      do ig = 1, SYM_NGROUPS
         g = get(sym_group_names(ig))
         do k = 1, g%order
            do a = 1, g%order
               do b = 1, g%order
                  c = g%op_index(ieor(g%ops(a), g%ops(b)))
                  call check(error, g%chars(k, a)*g%chars(k, b), g%chars(k, c), &
                             "homomorphism in "//g%name)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_homomorphism

   subroutine test_axes(error)
      !! The irrep that x, y, z span, from the sign of each coordinate under
      !! each mask, against pyscf's: D2h x=B3u y=B2u z=B1u; C2v x=B1 y=B2
      !! z=A1; C2h z=Au (x, y span Bu); D2 x=B3 y=B2 z=B1; Cs z=A", x, y=A';
      !! C2 z=A; Ci all Au.
      type(error_type), allocatable, intent(out) :: error
      call expect_axes("D2h", "B3u", "B2u", "B1u", error)
      if (allocated(error)) return
      call expect_axes("C2v", "B1", "B2", "A1", error)
      if (allocated(error)) return
      call expect_axes("C2h", "Bu", "Bu", "Au", error)
      if (allocated(error)) return
      call expect_axes("D2", "B3", "B2", "B1", error)
      if (allocated(error)) return
      call expect_axes("Cs", "A'", "A'", 'A"', error)
      if (allocated(error)) return
      call expect_axes("C2", "B", "B", "A", error)
      if (allocated(error)) return
      call expect_axes("Ci", "Au", "Au", "Au", error)
      if (allocated(error)) return
      call expect_axes("C1", "A", "A", "A", error)
   end subroutine test_axes

   subroutine expect_axes(name, xn, yn, zn, error)
      character(len=*), intent(in) :: name, xn, yn, zn
      type(error_type), allocatable, intent(inout) :: error
      type(trc_sym_group_t) :: g
      character(len=3) :: want(3)
      integer :: ax, k, found, sgn(8)
      want = [character(len=3) :: xn, yn, zn]
      g = get(name)
      do ax = 1, 3
         do k = 1, g%order
            sgn(k) = 1 - 2*ibits(g%ops(k), ax - 1, 1)
         end do
         found = 0
         do k = 1, g%order
            if (all(g%chars(k, :) == sgn(1:g%order))) found = k
         end do
         call check(error, found > 0, "axis spans an irrep of "//name)
         if (allocated(error)) return
         call check(error, trim(g%irrep_names(found)), trim(want(ax)), &
                    "axis irrep in "//name)
         if (allocated(error)) return
      end do
   end subroutine expect_axes

   subroutine test_product(error)
      !! prod(i, j) has the characters chi_i * chi_j, is symmetric, has the
      !! totally symmetric irrep as identity, and squares to it.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: g
      integer :: ig, i, j
      do ig = 1, SYM_NGROUPS
         g = get(sym_group_names(ig))
         do i = 1, g%order
            do j = 1, g%order
               call check(error, g%prod(i, j) > 0, "product exists in "//g%name)
               if (allocated(error)) return
               call check(error, all(g%chars(g%prod(i, j), :) == g%chars(i, :)*g%chars(j, :)), &
                          "product characters in "//g%name)
               if (allocated(error)) return
               call check(error, g%prod(i, j), g%prod(j, i), "product symmetric in "//g%name)
               if (allocated(error)) return
            end do
            call check(error, g%prod(1, i), i, "identity irrep in "//g%name)
            if (allocated(error)) return
            call check(error, g%prod(i, i), 1, "irrep squares to A in "//g%name)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_product

   subroutine test_unknown(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: g
      type(error_t) :: err
      g = trc_sym_group("D3h", err)
      call check(error, err%has_error(), "D3h refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "D3h") > 0, "message names the group")
      if (allocated(error)) return
      call check(error, g%order, 0, "empty group on error")
      if (allocated(error)) return
      g = trc_sym_group("", err)
      call check(error, err%has_error(), "empty name refused")
   end subroutine test_unknown

   subroutine test_case(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: g
      type(error_t) :: err
      g = trc_sym_group("c2V", err)
      call check(error, .not. err%has_error(), "c2V accepted")
      if (allocated(error)) return
      call check(error, g%name, "C2v", "canonical name")
   end subroutine test_case

end module check_sym_group_suite

program check_sym_group
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_sym_group_suite, only: collect_sym_group
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_sym_group, error_unit, stat)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_group
