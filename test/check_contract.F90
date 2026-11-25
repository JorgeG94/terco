module check_contract_suite
   !!
   !! trc_contract against plain loops.
   !!
   !!   1. contract_permute: every permutation of two, three and four indices,
   !!      with beta = 0 (destination prefilled with NaN: it must not be read)
   !!      and accumulating.
   !!   2. Every pairwise contraction the RI-CCSD Lambda equations and
   !!      densities need (INVENTORY: their index patterns, relabelled in order
   !!      of appearance), each at two sets of extents, with beta = 0 into NaN
   !!      and with alpha, beta; the scratch beyond `work_size` must stay
   !!      untouched. The patterns together must reach every path of the
   !!      plan: each array reordered, both GEMM transposes, swapped operands,
   !!      the strided batch.
   !!      The inventory runs at extents that include 1, and again asking for the
   !!      least scratch; two further patterns are batches that reorder the shared
   !!      operand and C.
   !!   3. The plans that must need no scratch, among them the batched one that
   !!      leaves a b(c,d,P) factor in place, and one that must; both sides of
   !!      the speed rule and of the memory bound.
   !!   4. Empty extents, and a NaN beta on every path.
   !!   5. Refusals, among them extents whose int64 product wraps.
   !!
   !! A `work` shorter than `work_size` stops the program, so it is not tested
   !! here.
   !!
   !! Every index letter has its own extent, so a confusion of two of them
   !! changes a shape or a value. The fills are small integers: every
   !! summation order is exact and the comparison is ==, on the host and with
   !! cuBLAS.
   !!
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
   use pic_types, only: dp, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_contract, only: trc_contract_plan_t, contract_permute
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_contract, la

   type(trc_linalg_t) :: la
      !! shared by the tests: initialised and released by the program

   integer, parameter :: PAD = 8  !! guard elements past `work_size`
   real(dp), parameter :: SENTINEL = -12345.0_dp

   type :: seen_t
      !! Which paths of the plan the contractions run so far have reached
      logical :: reorder_a = .false., reorder_b = .false., reorder_c = .false.
      logical :: trans_a = .false., trans_b = .false., swapped = .false.
      integer :: batched = 0              !! plans with more than one GEMM
      integer :: reordered = 0            !! plans that reorder an array
      integer :: batched_reordered = 0    !! plans that do both
   end type seen_t

   ! the index patterns of the pairwise steps of the equations' einsum form, factors b(p,q,P) included
   character(len=15), parameter :: INVENTORY(118) = [ &
      "ab,abc->c      ", "ab,ac->bc      ", "ab,ac->cb      ", "ab,bac->c      ", "ab,bc->ac      ", &
      "ab,bc->ca      ", "ab,c->abc      ", "ab,c->bac      ", "ab,ca->cb      ", "ab,cb->ac      ", &
      "ab,cb->ca      ", "abc,ab->c      ", "abc,ba->c      ", "abc,c->ab      ", "ab,acbd->cd    ", &
      "ab,acd->bcd    ", "ab,acd->cbd    ", "ab,bcd->acd    ", "ab,cabd->cd    ", "ab,cad->bcd    ", &
      "ab,cad->cbd    ", "ab,cadb->cd    ", "ab,cbd->acd    ", "ab,cbd->cad    ", "ab,cd->abdc    ", &
      "ab,cd->acbd    ", "ab,cdba->cd    ", "abc,abd->cd    ", "abc,acd->bd    ", "abc,adc->db    ", &
      "abc,bd->acd    ", "abc,cdb->ad    ", "abc,dac->db    ", "abc,db->acd    ", "abc,dbc->da    ", &
      "abc,dc->abd    ", "abc,dca->db    ", "abcd,ab->cd    ", "abcd,ac->bd    ", "abcd,ad->bc    ", &
      "abcd,bc->ad    ", "abcd,bd->ac    ", "abcd,bd->ca    ", "abcd,cb->ad    ", "abcd,cd->ab    ", &
      "abcd,da->cb    ", "abcd,db->ac    ", "abcd,dc->ab    ", "ab,acde->cdbe  ", "ab,bcde->adce  ", &
      "ab,bcde->dace  ", "ab,cade->cdbe  ", "ab,cbde->cade  ", "ab,cdae->cdbe  ", "ab,cdbe->cdae  ", &
      "ab,cdea->cdeb  ", "ab,cdea->ecbd  ", "ab,cdeb->cdea  ", "abc,abde->cde  ", "abc,cde->abed  ", &
      "abc,dabe->cde  ", "abc,daeb->cde  ", "abc,dce->abde  ", "abc,dce->dbae  ", "abc,debc->ade  ", &
      "abc,dec->abde  ", "abc,dec->adbe  ", "abc,decb->ade  ", "abcd,abce->de  ", "abcd,abce->ed  ", &
      "abcd,abe->cde  ", "abcd,abed->ec  ", "abcd,ace->bde  ", "abcd,aced->eb  ", "abcd,ade->bce  ", &
      "abcd,ae->cbed  ", "abcd,aecd->be  ", "abcd,bae->cde  ", "abcd,bce->ade  ", "abcd,bde->ace  ", &
      "abcd,caeb->ed  ", "abcd,cde->abe  ", "abcd,ce->abed  ", "abcd,de->abce  ", "abcd,eb->acde  ", &
      "abcd,eb->adce  ", "abcd,eb->aecd  ", "abcd,eb->cdae  ", "abcd,ebca->ed  ", "abcd,ebcd->ae  ", &
      "abcd,ebcd->ea  ", "abcd,ec->abde  ", "abcd,ec->abed  ", "abcd,ecbd->ae  ", "abcd,ed->abce  ", &
      "abcd,edca->eb  ", "abcd,bedf->aecf", "abcd,cdef->abef", "abcd,cedf->abef", "abcd,cefb->aefd", &
      "abcd,eabf->cdfe", "abcd,eafb->cdfe", "abcd,eafb->cefd", "abcd,eafb->efcd", "abcd,ebaf->cedf", &
      "abcd,ebaf->cefd", "abcd,ebcf->aedf", "abcd,ebcf->daef", "abcd,ebfa->cedf", "abcd,ebfc->daef", &
      "abcd,ebfd->caef", "abcd,ecbf->adfe", "abcd,efac->efbd", "abcd,efbd->aecf", "abcd,efcb->aedf", &
      "abcd,efcd->abef", "abcd,efdb->aecf", "abcd,efdc->abef"]

contains

   subroutine collect_contract(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("1. every permutation of 2, 3 and 4 indices", test_permute), &
                  new_unittest("2. every contraction of the inventory against loops", test_inventory), &
                  new_unittest("3. which plans need scratch", test_plans), &
                  new_unittest("4. empty extents and NaN beta", test_edges), &
                  new_unittest("5. refusals", test_refusals) &
                  ]
   end subroutine collect_contract

   pure function value_at(seed, i) result(v)
      !! A small integer, so that every sum and product is exact
      integer, intent(in) :: seed, i
      real(dp) :: v
      v = real(modulo(7*i + 3*(i/11) + seed, 13) - 6, dp)
   end function value_at

   pure function offset_of(letters, all_letters, position, extent) result(at)
      !! The 1-based column-major offset of the tensor `letters` at the index values `position`
      character(len=*), intent(in) :: letters, all_letters
      integer, intent(in) :: position(:), extent(:)
      integer :: at, stride, i, k

      at = 1
      stride = 1
      do i = 1, len(letters)
         k = index(all_letters, letters(i:i))
         at = at + position(k)*stride
         stride = stride*extent(k)
      end do
   end function offset_of

   subroutine next_position(position, extent)
      !! Odometer over the index values
      integer, intent(inout) :: position(:)
      integer, intent(in) :: extent(:)
      integer :: k

      do k = 1, size(position)
         position(k) = position(k) + 1
         if (position(k) < extent(k)) return
         position(k) = 0
      end do
   end subroutine next_position

   subroutine check_one_permutation(from, to, extent, error)
      !! contract_permute "from->to" against a loop, copy then accumulate
      character(len=*), intent(in) :: from, to
      integer, intent(in) :: extent(:)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: status
      real(dp), allocatable :: source(:), destination(:), expected(:)
      integer :: position(4), total, t, i
      character(len=:), allocatable :: spec

      spec = from//"->"//to
      total = product(extent)
      allocate (source(total), destination(total), expected(total))
      do i = 1, total
         source(i) = value_at(1, i)
      end do
      destination = ieee_value(0.0_dp, ieee_quiet_nan)
      !$acc enter data copyin(source, destination)
      call contract_permute(spec, extent, source, destination, status)
      !$acc update self(destination)

      position = 0
      do t = 1, total
         expected(offset_of(to, from, position, extent)) = source(offset_of(from, from, position, extent))
         call next_position(position(1:len(from)), extent)
      end do
      call check(error, .not. status%has_error() .and. all(destination == expected), "permute "//spec)
      if (.not. allocated(error)) then
         call contract_permute(spec, extent, source, destination, status, alpha=2.0_dp, beta=-3.0_dp)
         !$acc update self(destination)
         expected = 2.0_dp*expected - 3.0_dp*expected
         call check(error, all(destination == expected), "permute and accumulate "//spec)
      end if
      !$acc exit data delete(source, destination)
   end subroutine check_one_permutation

   subroutine test_permute(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=4), parameter :: NAMES = "abcd"
      integer, parameter :: EXTENT(4) = [3, 5, 4, 7]
      character(len=4) :: to
      integer :: rank, code, i, k, digits(4), count

      do rank = 1, 4
         count = 0
         ! every word of `rank` letters over the first `rank` letters, keeping the permutations
         do code = 0, rank**rank - 1
            k = code
            do i = 1, rank
               digits(i) = modulo(k, rank) + 1
               k = k/rank
               to(i:i) = NAMES(digits(i):digits(i))
            end do
            if (any([(index(to(1:rank), NAMES(i:i)) == 0, i=1, rank)])) cycle
            count = count + 1
            call check_one_permutation(NAMES(1:rank), to(1:rank), EXTENT(1:rank), error)
            if (allocated(error)) return
         end do
         call check(error, count == product([(i, i=1, rank)]), "the permutation count")
         if (allocated(error)) return
      end do
   end subroutine test_permute

   subroutine check_one_contraction(spec, extent_of_letter, seen, error, extra_work_limit)
      !! "x,y->z" against a loop over every index, with the extent of each letter 'a', 'b', ... given
      character(len=*), intent(in) :: spec
      integer, intent(in) :: extent_of_letter(:)
      type(seen_t), intent(inout) :: seen
      type(error_type), allocatable, intent(out) :: error
      integer(int64), intent(in), optional :: extra_work_limit
      type(trc_contract_plan_t) :: plan
      type(error_t) :: status
      character(len=:), allocatable :: x, y, z, letters
      real(dp), allocatable :: a(:), b(:), c(:), work(:), expected(:), before(:)
      integer :: extent(8), position(8), comma, arrow, i, t, na, nb, nc, nw, combinations

      comma = index(spec, ",")
      arrow = index(spec, "->")
      x = spec(1:comma - 1)
      y = spec(comma + 1:arrow - 1)
      z = trim(spec(arrow + 2:))
      letters = x
      do i = 1, len(y)
         if (index(letters, y(i:i)) == 0) letters = letters//y(i:i)
      end do
      do i = 1, len(letters)
         extent(i) = extent_of_letter(iachar(letters(i:i)) - iachar("a") + 1)
      end do
      na = size_of_tensor(x, letters, extent)
      nb = size_of_tensor(y, letters, extent)
      nc = size_of_tensor(z, letters, extent)

      if (present(extra_work_limit)) then
         plan%extra_work_limit = extra_work_limit
      end if
      call plan%build(spec(1:arrow + 1)//z, [(extent(index(letters, x(i:i))), i=1, len(x))], &
                      [(extent(index(letters, y(i:i))), i=1, len(y))], status)
      call check(error, plan%ready .and. .not. status%has_error(), "plan "//spec//" "//status%get_message())
      if (allocated(error)) return
      seen%reorder_a = seen%reorder_a .or. plan%reorders_a
      seen%reorder_b = seen%reorder_b .or. plan%reorders_b
      seen%reorder_c = seen%reorder_c .or. plan%reorders_c
      seen%trans_a = seen%trans_a .or. plan%gemm%transa == "T"
      seen%trans_b = seen%trans_b .or. plan%gemm%transb == "T"
      seen%swapped = seen%swapped .or. plan%swapped
      if (plan%gemm%count > 1) seen%batched = seen%batched + 1
      if (plan%work_size > 0) seen%reordered = seen%reordered + 1
      if (plan%gemm%count > 1 .and. plan%work_size > 0) seen%batched_reordered = seen%batched_reordered + 1

      nw = int(plan%work_size)
      allocate (a(na), b(nb), c(nc), expected(nc), before(nc), work(nw + PAD))
      do i = 1, na
         a(i) = value_at(2, i)
      end do
      do i = 1, nb
         b(i) = value_at(5, i)
      end do
      do i = 1, nc
         before(i) = value_at(9, i)
      end do
      c = ieee_value(0.0_dp, ieee_quiet_nan)
      work = SENTINEL

      expected = 0.0_dp
      combinations = product(extent(1:len(letters)))
      position = 0
      do t = 1, combinations
         associate (ic => offset_of(z, letters, position, extent))
            expected(ic) = expected(ic) + a(offset_of(x, letters, position, extent))*b(offset_of(y, letters, position, extent))
         end associate
         call next_position(position(1:len(letters)), extent(1:len(letters)))
      end do

      !$acc enter data copyin(a, b, c, work)
      if (nw == 0) then
         call plan%apply(la, a, b, c, work(1:0))
      else
         call plan%apply(la, a, b, c, work)
      end if
      !$acc update self(c)
      call check(error, all(c == expected), "contraction "//spec)
      if (.not. allocated(error)) then
         c = before
         !$acc update device(c)
         call plan%apply(la, a, b, c, work, alpha=2.0_dp, beta=-3.0_dp)
         !$acc update self(c, work)
         call check(error, all(c == 2.0_dp*expected - 3.0_dp*before), "alpha/beta contraction "//spec)
      end if
      if (.not. allocated(error)) call check(error, all(work(nw + 1:) == SENTINEL), "scratch overrun "//spec)
      !$acc exit data delete(a, b, c, work)
   end subroutine check_one_contraction

   pure function size_of_tensor(tensor, letters, extent) result(n)
      !! Elements of the tensor with these index letters
      character(len=*), intent(in) :: tensor, letters
      integer, intent(in) :: extent(:)
      integer :: n, i

      n = 1
      do i = 1, len(tensor)
         n = n*extent(index(letters, tensor(i:i)))
      end do
   end function size_of_tensor

   subroutine test_inventory(error)
      type(error_type), allocatable, intent(out) :: error
      ! extents of the letters a, b, ...: all different in the first two sets; then with ones, then all ones
      integer, parameter :: EXTENTS(8, 5) = reshape([2, 3, 5, 4, 7, 6, 9, 8, &
                                                     7, 6, 4, 5, 2, 3, 8, 9, &
                                                     1, 3, 5, 1, 7, 1, 9, 8, &
                                                     4, 1, 1, 3, 1, 5, 9, 8, &
                                                     1, 1, 1, 1, 1, 1, 1, 1], [8, 5])
      type(seen_t) :: seen, least, deep
      integer :: i, k

      do k = 1, size(EXTENTS, 2)
         do i = 1, size(INVENTORY)
            call check_one_contraction(trim(INVENTORY(i)), EXTENTS(:, k), seen, error)
            if (allocated(error)) return
         end do
      end do
      call check(error, seen%reorder_a .and. seen%reorder_b .and. seen%reorder_c, "each array reordered somewhere")
      if (allocated(error)) return
      call check(error, seen%trans_a .and. seen%trans_b .and. seen%swapped, "transposes and swapped operands reached")
      if (allocated(error)) return
      call check(error, seen%batched > 0 .and. seen%reordered > 0, "the strided batch and the reorderings reached")
      if (allocated(error)) return

      do i = 1, size(INVENTORY)
         call check_one_contraction(trim(INVENTORY(i)), EXTENTS(:, 1), least, error, extra_work_limit=0_int64)
         if (allocated(error)) return
      end do
      call check(error, least%batched > 0 .and. least%batched_reordered > 0, "batches reached at the least scratch")
      if (allocated(error)) return

      ! a batch that reorders the shared operand and C, the batched one being A and then B
      call check_one_contraction("cbd,ebf->ecfd", EXTENTS(:, 1), deep, error)
      if (allocated(error)) return
      call check(error, deep%batched_reordered == 1 .and. deep%reorder_b .and. deep%reorder_c .and. &
                 .not. deep%reorder_a, "a batch over A that reorders B and C")
      if (allocated(error)) return
      deep = seen_t()
      call check_one_contraction("ebf,cbd->ecfd", EXTENTS(:, 1), deep, error)
      if (allocated(error)) return
      call check(error, deep%batched_reordered == 1 .and. deep%reorder_a .and. deep%reorder_c .and. &
                 .not. deep%reorder_b, "a batch over B that reorders A and C")
   end subroutine test_inventory

   subroutine expect_plan(spec, extent_a, extent_b, work_size, count, error, extra_work_limit)
      !! The scratch a plan asks for and the number of GEMMs in it
      character(len=*), intent(in) :: spec
      integer, intent(in) :: extent_a(:), extent_b(:), work_size, count
      type(error_type), allocatable, intent(out) :: error
      integer(int64), intent(in), optional :: extra_work_limit
      type(trc_contract_plan_t) :: plan
      type(error_t) :: status

      if (present(extra_work_limit)) then
         plan%extra_work_limit = extra_work_limit
      end if
      call plan%build(spec, extent_a, extent_b, status)
      call check(error, plan%ready, "plan "//spec)
      if (allocated(error)) return
      call check(error, plan%work_size == int(work_size, int64), "work_size of "//spec)
      if (allocated(error)) return
      call check(error, plan%gemm%count == count, "GEMM count of "//spec)
   end subroutine expect_plan

   subroutine test_plans(error)
      type(error_type), allocatable, intent(out) :: error

      ! the four plain GEMM layouts, and C with B's free index first
      call expect_plan("ik,kj->ij", [3, 4], [4, 5], 0, 1, error)
      if (allocated(error)) return
      call expect_plan("ki,jk->ij", [4, 3], [5, 4], 0, 1, error)
      if (allocated(error)) return
      call expect_plan("ik,kj->ji", [3, 4], [4, 5], 0, 1, error)
      if (allocated(error)) return
      ! compound indices that are already adjacent
      call expect_plan("ijab,abkl->ijkl", [2, 3, 4, 5], [4, 5, 6, 7], 0, 1, error)
      if (allocated(error)) return
      ! t(k,d) with b(c,d,P): one GEMM per P, neither array copied
      call expect_plan("kd,cdP->kcP", [3, 4], [5, 4, 6], 0, 6, error)
      if (allocated(error)) return
      ! a ring in a batch over b while the contracted pair (d,l) is narrow; one GEMM and a reordered C
      ! once it is wide enough that the items would cost more than the copy
      call expect_plan("dlib,jlda->ijab", [2, 2, 2, 3], [2, 2, 2, 3], 24, 3, error)
      if (allocated(error)) return
      call expect_plan("dlib,jlda->ijab", [20, 20, 2, 3], [2, 20, 20, 3], 2400 + 36, 1, error)
      if (allocated(error)) return
      ! an outer product that interleaves its operands needs C once
      call expect_plan("ia,jb->ijab", [2, 3], [4, 5], 120, 1, error)
      if (allocated(error)) return
      ! reorder the smaller operand, not the larger one
      call expect_plan("ikcd,cdk->i", [2, 3, 4, 5], [4, 5, 3], 60, 1, error)
      if (allocated(error)) return

      ! A(a,b) with a factor B(c,b,P), v = a = b = c and naux = 4 v: 2 v flops per element a batch saves.
      ! Batched up to v = 250; at 251 the factor is copied, while the copy is within 2^26 doubles
      ! (v = 256); above that, batched again. Plans only: nothing of this size is allocated.
      call expect_plan("ab,cbd->acd", [250, 250], [250, 250, 1000], 0, 1000, error)
      if (allocated(error)) return
      call expect_plan("ab,cbd->acd", [251, 251], [251, 251, 1004], 251*251*1004, 1, error)
      if (allocated(error)) return
      call expect_plan("ab,cbd->acd", [256, 256], [256, 256, 1024], 256*256*1024, 1, error)
      if (allocated(error)) return
      call expect_plan("ab,cbd->acd", [257, 257], [257, 257, 1028], 0, 1028, error)
      if (allocated(error)) return
      call expect_plan("ab,cbd->acd", [600, 600], [600, 600, 2400], 0, 2400, error)
      if (allocated(error)) return
      ! the caller's own bound
      call expect_plan("ab,cbd->acd", [251, 251], [251, 251, 1004], 0, 1004, error, extra_work_limit=0_int64)
      if (allocated(error)) return
      call expect_plan("ab,cbd->acd", [600, 600], [600, 600, 2400], 600*600*2400, 1, error, &
                       extra_work_limit=int(huge(0), int64))
   end subroutine test_plans

   subroutine check_beta(spec, extent_a, extent_b, nc, error)
      !! A NaN beta must reach every element of C, as it does through BLAS
      character(len=*), intent(in) :: spec
      integer, intent(in) :: extent_a(:), extent_b(:), nc
      type(error_type), allocatable, intent(out) :: error
      type(trc_contract_plan_t) :: plan
      type(error_t) :: status
      real(dp), allocatable :: a(:), b(:), c(:), work(:)

      call plan%build(spec, extent_a, extent_b, status)
      call check(error, plan%ready, "plan "//spec)
      if (allocated(error)) return
      allocate (a(max(product(extent_a), 1)), b(max(product(extent_b), 1)), c(nc), work(max(int(plan%work_size), 1)))
      a = 1.0_dp
      b = 1.0_dp
      c = 1.0_dp
      work = 0.0_dp
      !$acc enter data copyin(a, b, c, work)
      call plan%apply(la, a, b, c, work, beta=ieee_value(0.0_dp, ieee_quiet_nan))
      !$acc update self(c)
      !$acc exit data delete(a, b, c, work)
      call check(error, all(ieee_is_nan(c)), "NaN beta, "//spec)
   end subroutine check_beta

   subroutine check_empty_sum(error)
      !! A contracted extent of 0: C = beta C, whatever C holds when beta = 0
      type(error_type), allocatable, intent(out) :: error
      type(trc_contract_plan_t) :: plan
      type(error_t) :: status
      real(dp) :: none(0), c(12)

      call plan%build("ak,kb->ab", [3, 0], [0, 4], status)
      call check(error, plan%ready .and. plan%work_size == 0, "plan with an empty contracted index")
      if (allocated(error)) return
      c = 2.0_dp
      !$acc enter data copyin(c)
      call plan%apply(la, none, none, c, none, beta=-3.0_dp)
      !$acc update self(c)
      call check(error, all(c == -6.0_dp), "empty sum: C = beta C")
      if (.not. allocated(error)) then
         c = ieee_value(0.0_dp, ieee_quiet_nan)
         !$acc update device(c)
         call plan%apply(la, none, none, c, none)
         !$acc update self(c)
         call check(error, all(c == 0.0_dp), "empty sum, beta = 0: C = 0")
      end if
      !$acc exit data delete(c)
   end subroutine check_empty_sum

   subroutine test_edges(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_contract_plan_t) :: plan
      type(error_t) :: status
      real(dp) :: none(0), source(6), destination(6)

      call check_empty_sum(error)
      if (allocated(error)) return

      ! an empty C: nothing to do, with an empty or a non-empty contracted index
      call plan%build("ak,kb->ab", [0, 3], [3, 4], status)
      call check(error, plan%ready .and. plan%work_size == 0, "plan with an empty free index")
      if (allocated(error)) return
      call plan%apply(la, none, none, none, none)
      call contract_permute("ab->ba", [0, 3], none, none, status)
      call check(error, .not. status%has_error(), "an empty permutation is not an error")
      if (allocated(error)) return

      ! NaN beta: direct GEMM, batch, reordered C, empty sum, and the permutation
      call check_beta("ak,kb->ab", [3, 2], [2, 4], 12, error)
      if (allocated(error)) return
      call check_beta("kd,cdP->kcP", [3, 4], [5, 4, 6], 90, error)
      if (allocated(error)) return
      call check_beta("ia,jb->ijab", [2, 3], [4, 5], 120, error)
      if (allocated(error)) return
      call check_beta("ak,kb->ab", [3, 0], [0, 4], 12, error)
      if (allocated(error)) return
      source = 1.0_dp
      destination = 1.0_dp
      !$acc enter data copyin(source, destination)
      call contract_permute("ab->ba", [2, 3], source, destination, status, beta=ieee_value(0.0_dp, ieee_quiet_nan))
      !$acc update self(destination)
      !$acc exit data delete(source, destination)
      call check(error, all(ieee_is_nan(destination)), "NaN beta, permutation")
   end subroutine test_edges

   subroutine expect_refusal(spec, extent_a, extent_b, error)
      character(len=*), intent(in) :: spec
      integer, intent(in) :: extent_a(:), extent_b(:)
      type(error_type), allocatable, intent(out) :: error
      type(trc_contract_plan_t) :: plan
      type(error_t) :: status

      call plan%build(spec, extent_a, extent_b, status)
      call check(error, status%has_error() .and. .not. plan%ready, "refusal of "//spec)
   end subroutine expect_refusal

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: status
      real(dp) :: source(6), destination(6)
      integer :: i
      character(len=16), parameter :: BAD(8) = [character(len=16) :: "ab,ab->ab", "ab,bc->a", "aa,ab->b", "ab,ab->", &
                                                "abcde,e->abcd", "ab,bc->ac,d", "ab;bc->ac", "a1,1c->ac"]

      do i = 1, size(BAD)
         call expect_refusal(trim(BAD(i)), [2, 2], [2, 2], error)
         if (allocated(error)) return
      end do
      call expect_refusal("ab,bc->ac", [2, 3], [4, 5], error)   ! b has two extents
      if (allocated(error)) return
      call expect_refusal("ab,bc->ac", [2, 3, 1], [3, 5], error)  ! rank
      if (allocated(error)) return
      call expect_refusal("ab,bc->ac", [2, -1], [-1, 5], error)   ! negative
      if (allocated(error)) return
      ! 65536^4 = 2^64 wraps to 0 in int64
      call expect_refusal("abcd,d->abc", [65536, 65536, 65536, 65536], [65536], error)
      if (allocated(error)) return
      call expect_refusal("d,abcd->abc", [65536], [65536, 65536, 65536, 65536], error)
      if (allocated(error)) return
      ! an empty contracted index does not excuse a C of 2^64 elements
      call expect_refusal("abe,cde->abcd", [65536, 65536, 0], [65536, 65536, 0], error)
      if (allocated(error)) return
      call expect_refusal("ab,cd->abcd", [1000, 1000], [1000, 1000], error)  ! C above 2^31-1
      if (allocated(error)) return

      source = 1.0_dp
      destination = 7.0_dp
      !$acc enter data copyin(source, destination)
      call contract_permute("ab->aa", [2, 3], source, destination, status)
      call check(error, status%has_error(), "a repeated letter is no permutation")
      if (.not. allocated(error)) then
         call status%clear()
         call contract_permute("ab->bc", [2, 3], source, destination, status)
         call check(error, status%has_error(), "a new letter is no permutation")
      end if
      if (.not. allocated(error)) then
         call status%clear()
         call contract_permute("ab->ba", [2, 3, 1], source, destination, status)
         call check(error, status%has_error(), "the rank must be the extents'")
      end if
      if (.not. allocated(error)) then
         call status%clear()
         call contract_permute("abcd->dcba", [65536, 65536, 65536, 65536], source, destination, status)
         call check(error, status%has_error(), "2^64 elements are refused, not wrapped to none")
      end if
      !$acc update self(destination)
      !$acc exit data delete(source, destination)
      if (.not. allocated(error)) call check(error, all(destination == 7.0_dp), "a refusal writes nothing")
   end subroutine test_refusals

end module check_contract_suite

program check_contract
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_contract_suite, only: collect_contract, la
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   call la%init(1)
   testsuites = [new_testsuite("check_contract", collect_contract)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   call la%release()

   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_contract
