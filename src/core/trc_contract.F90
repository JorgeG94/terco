!
! Pairwise contraction of tensors of up to four indices: permute, then GEMM.
!
! WHY. The RI-CCSD Lambda equations and densities, written without any o v^3
! or v^4 array, are about 180 distinct pairwise contractions of tensors of at
! most four indices. Taken in the stored index orders, some 110 of them need an
! operand or the result reordered before a GEMM can read it, between them in 19
! of the 23 reorderings of four indices. A kernel per order and the m, n, k and
! leading dimensions of every call written by hand is more code, and more
! places to confuse two dimensions, than one plan that derives them from the
! index letters.
!
! WHAT. A contraction is spelled as an einsum, "abcd,ed->abce", arrays stored
! with the first index fastest. Every letter appears in exactly two of the
! three tensors: no Hadamard index, no index summed within one tensor, no
! repeated letter, no scalar result (the equations need none of these apart
! from two traces, which their callers loop over).
!
!     type(trc_contract_plan_t) :: plan
!     call plan%build("ijcd,kd->ijck", shape(l2), shape(t1), error)
!     ! work: contiguous, at least plan%work_size doubles (none when that is 0), present on the device
!     call plan%apply(la, l2, t1, out, work, alpha, beta)
!
!     C = alpha contraction(A, B) + beta C
!
! alpha and beta as in BLAS: C is not read when beta = 0, and a NaN beta gives NaN on every path.
! An extent may be 0 (an empty virtual space): with C empty `apply` does nothing, with only the
! contracted index empty it leaves C = beta C. A tensor above 2^31-1 elements is refused.
!
! HOW. `build` picks, among the layouts a GEMM can read, the one that reorders
! the fewest elements:
!
!   - the contracted letters in A's order or in B's, the free letters of each
!     operand in the order they have in the operand or in C;
!   - A's or B's free letters first in C (swapping the GEMM's operands);
!   - if C's last letter is also the last letter of one operand, a strided
!     batched GEMM over it, the other operand shared by every item. This is
!     what keeps the factors b(c,d,P) in place when d is contracted against
!     t(k,d): a reordered copy of v^2 naux doubles costs more than the GEMMs
!     it feeds when the other side is only o wide.
!
! A BATCH ONLY FOR CHEAP GEMMS. The items of a batch are smaller than the one
! GEMM they replace and cuBLAS runs them slower, so the batch is taken only
! while the GEMM does at most BATCH_FLOPS_PER_ELEMENT flops per element it
! saves from being reordered. Measured on one V100, o = 20, v = 200,
! naux = 700, batched against reordered, over the 28 contractions where the
! question arises: with 40 flops per saved element the batch took 0.66 ms
! against 2.8 ms (t(k,c) b(b,c,P)); with 8000 it took 78 ms against 20 ms
! (an o^3 v^3 ring in 200 items, 1.7 TFLOP/s against 7); every loss seen
! was above 1000.
!
! AND NEVER AT ANY MEMORY PRICE. The copy that speed asks for may be a whole factor, v^2 naux
! doubles. The reordered layout is taken only while it needs at most EXTRA_WORK_LIMIT (2^26
! doubles, 512 MiB) more scratch than the batch; above that the batch is taken whatever its speed.
! A plan's `extra_work_limit`, set before `build`, replaces the bound: 0 asks for the least scratch.
! (A component, not an optional argument of `build`: nvfortran 26.9 -O3 read the absent int64 dummy
! of "limit = DEFAULT; if (present(x)) limit = x" and segfaulted.)
!
! Each array that still does not fit is reordered by one kernel into the
! caller's work array (A, then B, then C's GEMM result, which is then
! reordered into C). Nothing is allocated here.
!
! RESIDENCY. A, B, C and the work array must be present on the device, as for
! trc_linalg. None of them may overlap C or the work array.
!
module trc_contract
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: trc_contract_plan_t
   public :: contract_permute

   integer, parameter :: MAX_RANK = 4
   integer(int64), parameter :: BATCH_FLOPS_PER_ELEMENT = 500  !! module header, A BATCH ONLY FOR CHEAP GEMMS
   integer(int64), parameter :: EXTRA_WORK_LIMIT = 2_int64**26  !! doubles; module header, AND NEVER AT ANY MEMORY PRICE

   type :: reorder_t
      !! One reordering as the kernel takes it: the extents of the output, padded with 1 to four, and
      !! the distance in the source between neighbours along each output index
      integer :: extent(MAX_RANK) = 1
      integer :: stride(MAX_RANK) = 0
      integer :: total = 0  !! elements
   end type reorder_t

   type :: letters_t
      !! The three tensors' index letters, blank padded
      character(len=MAX_RANK) :: a = ""
      character(len=MAX_RANK) :: b = ""
      character(len=MAX_RANK) :: c = ""
   end type letters_t

   type :: layout_t
      !! A candidate: the orders the GEMM reads the arrays in, and which arrays must be reordered to get them
      character(len=MAX_RANK) :: a = ""
      character(len=MAX_RANK) :: b = ""
      character(len=MAX_RANK) :: c = ""
      character(len=MAX_RANK) :: contracted = ""
      character(len=MAX_RANK) :: free_a = ""
      character(len=MAX_RANK) :: free_b = ""
      character :: batch = ""     !! the batched letter, blank for one GEMM
      logical :: swapped = .false.  !! C is (free of B)(free of A): B is the GEMM's left operand
   end type layout_t

   type :: trc_contract_plan_t
      !! One contraction at one set of extents: plain values, no device data. Built by `build`, run by `apply`.
      logical :: ready = .false.
      integer(int64) :: extra_work_limit = EXTRA_WORK_LIMIT
         !! set before `build`: the most scratch, in doubles, a faster layout may cost over a batched one
      integer(int64) :: work_size = 0  !! doubles of scratch `apply` needs; 0 when nothing is reordered
      logical :: reorders_a = .false.
      logical :: reorders_b = .false.
      logical :: reorders_c = .false.
      logical :: swapped = .false.     !! B is the GEMM's left operand
      logical :: empty_sum = .false.   !! an operand is empty: no GEMM, C = beta C over its `elements_c`
      integer :: elements_c = 0        !! elements of C
      type(trc_gemm_batch_t) :: gemm   !! the GEMM (count = 1) or the strided batch, in left/right order
      type(reorder_t), private :: order_a, order_b, order_c
   contains
      procedure :: build => plan_build
      procedure :: apply => plan_apply
   end type trc_contract_plan_t

contains

   subroutine plan_build(this, spec, extent_a, extent_b, error)
      !! Plan "x,y->z" for A of shape `extent_a` and B of shape `extent_b`; C's shape follows from `spec`.
      !! A spec outside the scope of the module header, a negative extent, extents that disagree or a tensor
      !! above 2^31-1 elements set `error` and leave the plan not ready. `this%extra_work_limit` is read,
      !! not reset.
      class(trc_contract_plan_t), intent(inout) :: this
      character(len=*), intent(in) :: spec
      integer, intent(in) :: extent_a(:), extent_b(:)
      type(error_t), intent(inout) :: error
      type(letters_t) :: given
      type(layout_t) :: best
      integer :: size_of(MAX_RANK*2)
      integer :: i
      character(len=MAX_RANK*2) :: names
      integer(int64) :: elements(3)

      this%ready = .false.
      this%work_size = 0
      this%reorders_a = .false.
      this%reorders_b = .false.
      this%reorders_c = .false.
      this%swapped = .false.
      this%empty_sum = .false.
      this%elements_c = 0
      this%gemm = trc_gemm_batch(0, 0, 0, 0)
      if (.not. spec_is_valid(spec, given)) then
         call error%set(ERROR_VALIDATION, "trc_contract: not a supported pairwise contraction: "//trim(spec))
         return
      end if
      if (size(extent_a) /= len_trim(given%a) .or. size(extent_b) /= len_trim(given%b)) then
         call error%set(ERROR_VALIDATION, "trc_contract: operand ranks do not match "//trim(spec))
         return
      end if
      if (any(extent_a < 0) .or. any(extent_b < 0)) then
         call error%set(ERROR_VALIDATION, "trc_contract: a negative extent in "//trim(spec))
         return
      end if

      ! one extent per distinct letter
      names = given%a(1:size(extent_a))
      size_of = 1
      size_of(1:size(extent_a)) = extent_a
      do i = 1, size(extent_b)
         associate (at => index(trim(names), given%b(i:i)))
            if (at == 0) then
               names = trim(names)//given%b(i:i)
               size_of(len_trim(names)) = extent_b(i)
            else if (size_of(at) /= extent_b(i)) then
               call error%set(ERROR_VALIDATION, "trc_contract: index "//given%b(i:i)//" has two extents in "//trim(spec))
               return
            end if
         end associate
      end do
      elements(1) = bounded_product(extent_a)
      elements(2) = bounded_product(extent_b)
      elements(3) = bounded_product(shape_of(trim(given%c), names, size_of))
      if (any(elements > int(huge(0), int64))) then
         call error%set(ERROR_VALIDATION, "trc_contract: a tensor of "//trim(spec)//" has more than 2^31-1 elements")
         return
      end if
      this%elements_c = int(elements(3))
      if (minval(elements) == 0) then
         ! nothing to multiply; `apply` only scales whatever C there is
         this%empty_sum = .true.
         this%ready = .true.
         return
      end if

      best = choose_layout(given, names, size_of, this%extra_work_limit)

      this%swapped = best%swapped
      this%reorders_a = best%a /= given%a
      this%reorders_b = best%b /= given%b
      this%reorders_c = best%c /= given%c
      if (this%reorders_a) then
         this%order_a = reorder_of(trim(given%a), trim(best%a), extent_a)
         this%work_size = this%work_size + elements(1)
      end if
      if (this%reorders_b) then
         this%order_b = reorder_of(trim(given%b), trim(best%b), extent_b)
         this%work_size = this%work_size + elements(2)
      end if
      if (this%reorders_c) then
         ! C is reordered back: the source is the GEMM's result
         this%order_c = reorder_of(trim(best%c), trim(given%c), shape_of(trim(best%c), names, size_of))
         this%work_size = this%work_size + elements(3)
      end if
      this%gemm = gemm_of(best, names, size_of)
      this%ready = .true.
   end subroutine plan_build

   subroutine plan_apply(this, la, a, b, c, work, alpha, beta)
      !! C = alpha contraction(A, B) + beta C; alpha and beta default to 1 and 0. `work` holds at least
      !! `this%work_size` doubles, and may be empty when that is 0; a shorter one stops the program.
      !! `la` is initialised.
      class(trc_contract_plan_t), intent(in) :: this
      class(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: a(*), b(*)
      real(dp), intent(inout) :: c(*)
      real(dp), intent(inout), contiguous :: work(:)
      real(dp), intent(in), optional :: alpha, beta
      real(dp) :: al, be

      if (.not. this%ready) error stop "trc_contract: apply of a plan that was not built"
      if (size(work, kind=int64) < this%work_size) error stop "trc_contract: apply: work is shorter than work_size"
      al = 1.0_dp
      if (present(alpha)) al = alpha
      be = 0.0_dp
      if (present(beta)) be = beta

      if (this%empty_sum .or. ieee_is_nan(be)) then
         ! A NaN beta makes every element of C NaN whatever the product is, and it is done here because
         ! cublasDgemmStridedBatched did not carry it into C (cuBLAS's plain GEMM did)
         if (this%elements_c > 0) call scale_run(this%elements_c, be, c)
      else if (this%work_size == 0) then
         ! `work` may have no element to name
         if (this%swapped) then
            call gemm_run(la, this%gemm, b, a, c, al, be)
         else
            call gemm_run(la, this%gemm, a, b, c, al, be)
         end if
      else
         call apply_reordered(this, la, a, b, c, work, al, be)
      end if
   end subroutine plan_apply

   subroutine apply_reordered(plan, la, a, b, c, work, alpha, beta)
      !! `apply` of a plan that reorders something: `work` holds A's copy, then B's, then C's GEMM result
      class(trc_contract_plan_t), intent(in) :: plan
      class(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: a(*), b(*)
      real(dp), intent(inout) :: c(*)
      real(dp), intent(inout) :: work(*)
      real(dp), intent(in) :: alpha, beta
      integer(int64) :: at_a, at_b, at_c

      at_a = 1
      at_b = at_a
      if (plan%reorders_a) at_b = at_a + plan%order_a%total
      at_c = 1  ! unused unless C is reordered
      if (plan%reorders_c) then
         if (plan%reorders_a) at_c = at_c + plan%order_a%total
         if (plan%reorders_b) at_c = at_c + plan%order_b%total
      end if

      if (plan%reorders_a) call reorder_run(plan%order_a, a, work(at_a), 1.0_dp, 0.0_dp)
      if (plan%reorders_b) call reorder_run(plan%order_b, b, work(at_b), 1.0_dp, 0.0_dp)

      if (plan%reorders_a .and. plan%reorders_b) then
         call gemm_into(plan, la, work(at_a), work(at_b), c, work(at_c), alpha, beta)
      else if (plan%reorders_a) then
         call gemm_into(plan, la, work(at_a), b, c, work(at_c), alpha, beta)
      else if (plan%reorders_b) then
         call gemm_into(plan, la, a, work(at_b), c, work(at_c), alpha, beta)
      else
         call gemm_into(plan, la, a, b, c, work(at_c), alpha, beta)
      end if

      if (plan%reorders_c) call reorder_run(plan%order_c, work(at_c), c, 1.0_dp, beta)
   end subroutine apply_reordered

   subroutine gemm_into(plan, la, a, b, c, result_c, alpha, beta)
      !! The GEMM of `plan` on A and B as it reads them: into C, or into `result_c` when C is reordered
      !! afterwards (then beta waits for that step)
      class(trc_contract_plan_t), intent(in) :: plan
      class(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: a(*), b(*)
      real(dp), intent(inout) :: c(*), result_c(*)
      real(dp), intent(in) :: alpha, beta

      if (plan%reorders_c) then
         if (plan%swapped) then
            call gemm_run(la, plan%gemm, b, a, result_c, alpha, 0.0_dp)
         else
            call gemm_run(la, plan%gemm, a, b, result_c, alpha, 0.0_dp)
         end if
      else
         if (plan%swapped) then
            call gemm_run(la, plan%gemm, b, a, c, alpha, beta)
         else
            call gemm_run(la, plan%gemm, a, b, c, alpha, beta)
         end if
      end if
   end subroutine gemm_into

   subroutine gemm_run(la, g, left, right, c, alpha, beta)
      !! One GEMM goes to the plain call, several to the strided batched one
      class(trc_linalg_t), intent(in) :: la
      type(trc_gemm_batch_t), intent(in) :: g
      real(dp), intent(in) :: left(*), right(*)
      real(dp), intent(inout) :: c(*)
      real(dp), intent(in) :: alpha, beta

      if (g%count == 1) then
         call la%gemm(g%transa, g%transb, g%m, g%n, g%k, alpha, left, g%lda, right, g%ldb, beta, c, g%ldc)
      else
         call la%gemm_batched(g, left, right, c, alpha, beta)
      end if
   end subroutine gemm_run

   subroutine contract_permute(spec, extent, source, destination, error, alpha, beta)
      !! destination = alpha P(source) + beta destination, `spec` "abcd->acbd" naming the indices of the
      !! source (shape `extent`) and of the destination. alpha and beta default to 1 and 0, and the
      !! destination is not read when beta = 0. Both arrays are present on the device and distinct.
      !! An extent of 0 is nothing to do. A bad spec, a negative extent or more than 2^31-1 elements set
      !! `error` and write nothing.
      character(len=*), intent(in) :: spec
      integer, intent(in) :: extent(:)
      real(dp), intent(in) :: source(*)
      real(dp), intent(inout) :: destination(*)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: alpha, beta
      real(dp) :: al, be
      integer :: arrow, rank, i
      logical :: valid

      al = 1.0_dp
      if (present(alpha)) al = alpha
      be = 0.0_dp
      if (present(beta)) be = beta

      arrow = index(spec, "->")
      rank = size(extent)
      valid = arrow == rank + 1 .and. rank >= 1 .and. rank <= MAX_RANK .and. len_trim(spec) == 2*rank + 2
      if (valid) then
         do i = 1, rank
            valid = valid .and. count_of(spec(1:rank), spec(i:i)) == 1 &
                    .and. count_of(spec(arrow + 2:arrow + 1 + rank), spec(i:i)) == 1
         end do
      end if
      if (valid) valid = all(extent >= 0)
      if (valid) valid = bounded_product(extent) <= int(huge(0), int64)
      if (.not. valid) then
         call error%set(ERROR_VALIDATION, "trc_contract: not a permutation of these extents: "//trim(spec))
         return
      end if
      if (any(extent == 0)) return

      call reorder_run(reorder_of(spec(1:rank), spec(arrow + 2:arrow + 1 + rank), extent), source, destination, al, be)
   end subroutine contract_permute

   pure function bounded_product(extent) result(n)
      !! The product of non-negative extents, or some value above 2^31-1 as soon as it passes that: a
      !! plain int64 product of four extents can wrap to a small number
      integer, intent(in) :: extent(:)
      integer(int64) :: n
      integer :: i

      n = 0
      if (any(extent == 0)) return
      n = 1
      do i = 1, size(extent)
         n = n*extent(i)
         if (n > int(huge(0), int64)) return
      end do
   end function bounded_product

   pure function reorder_of(from, to, extent) result(order)
      !! How to read a tensor stored as `from` (shape `extent`) in the order `to`
      character(len=*), intent(in) :: from, to
      integer, intent(in) :: extent(:)
      type(reorder_t) :: order
      integer :: i, at

      do i = 1, len(to)
         at = index(from, to(i:i))
         order%extent(i) = extent(at)
         order%stride(i) = product(extent(1:at - 1))
      end do
      order%total = product(extent)
   end function reorder_of

   subroutine reorder_run(order, source, destination, alpha, beta)
      !! Launch the kernel of `order`; beta = 0 takes the one that does not read the destination. The test
      !! is written so that a NaN beta accumulates, as BLAS does with it.
      type(reorder_t), intent(in) :: order
      real(dp), intent(in) :: source(*)
      real(dp), intent(inout) :: destination(*)
      real(dp), intent(in) :: alpha, beta

      if (abs(beta) <= 0.0_dp) then
         call kernel_reorder(order%total, order%extent(1), order%extent(2), order%extent(3), order%stride(1), &
                             order%stride(2), order%stride(3), order%stride(4), alpha, source, destination)
      else
         call kernel_reorder_add(order%total, order%extent(1), order%extent(2), order%extent(3), order%stride(1), &
                                 order%stride(2), order%stride(3), order%stride(4), alpha, source, beta, destination)
      end if
   end subroutine reorder_run

   subroutine scale_run(total, beta, c)
      !! C = beta C, C not read when beta = 0
      integer, intent(in) :: total
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: c(total)

      if (abs(beta) <= 0.0_dp) then
         call kernel_zero(total, c)
      else
         call kernel_scale(total, beta, c)
      end if
   end subroutine scale_run

   subroutine kernel_zero(total, c)
      !! C = 0
      integer, intent(in) :: total
      real(dp), intent(inout) :: c(total)
      integer :: t

      do concurrent(t=1:total)
         c(t) = 0.0_dp
      end do
   end subroutine kernel_zero

   subroutine kernel_scale(total, beta, c)
      !! C = beta C
      integer, intent(in) :: total
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: c(total)
      integer :: t

      do concurrent(t=1:total)
         c(t) = beta*c(t)
      end do
   end subroutine kernel_scale

   subroutine kernel_reorder(total, n1, n2, n3, s1, s2, s3, s4, alpha, source, destination)
      !! destination(i1,i2,i3,i4) = alpha source(1 + i1 s1 + i2 s2 + i3 s3 + i4 s4), indices from 0.
      !! One flattened index, in the destination's storage order: the writes are contiguous.
      integer, intent(in) :: total, n1, n2, n3, s1, s2, s3, s4
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: source(total)
      real(dp), intent(inout) :: destination(total)
      integer :: t, q, r, i1, i2, i3, i4

      do concurrent(t=0:total - 1) local(q, r, i1, i2, i3, i4)
         q = t/n1
         i1 = t - q*n1
         r = q/n2
         i2 = q - r*n2
         i4 = r/n3
         i3 = r - i4*n3
         destination(t + 1) = alpha*source(1 + i1*s1 + i2*s2 + i3*s3 + i4*s4)
      end do
   end subroutine kernel_reorder

   subroutine kernel_reorder_add(total, n1, n2, n3, s1, s2, s3, s4, alpha, source, beta, destination)
      !! `kernel_reorder` with beta destination added
      integer, intent(in) :: total, n1, n2, n3, s1, s2, s3, s4
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: source(total)
      real(dp), intent(inout) :: destination(total)
      integer :: t, q, r, i1, i2, i3, i4

      do concurrent(t=0:total - 1) local(q, r, i1, i2, i3, i4)
         q = t/n1
         i1 = t - q*n1
         r = q/n2
         i2 = q - r*n2
         i4 = r/n3
         i3 = r - i4*n3
         destination(t + 1) = alpha*source(1 + i1*s1 + i2*s2 + i3*s3 + i4*s4) + beta*destination(t + 1)
      end do
   end subroutine kernel_reorder_add

   function spec_is_valid(spec, given) result(valid)
      !! Split "x,y->z" and check the scope of the module header: ranks 1 to 4, no letter twice in one
      !! tensor, every letter in exactly two of the three
      character(len=*), intent(in) :: spec
      type(letters_t), intent(out) :: given
      logical :: valid
      integer :: comma, arrow, last, i
      character(len=3*MAX_RANK) :: every

      valid = .false.
      comma = index(spec, ",")
      arrow = index(spec, "->")
      last = len_trim(spec)
      if (comma < 2 .or. arrow < comma + 2 .or. last < arrow + 2) return
      if (comma - 1 > MAX_RANK .or. arrow - comma - 1 > MAX_RANK .or. last - arrow - 1 > MAX_RANK) return
      given%a = spec(1:comma - 1)
      given%b = spec(comma + 1:arrow - 1)
      given%c = spec(arrow + 2:last)
      every = trim(given%a)//trim(given%b)//trim(given%c)
      do i = 1, len_trim(every)
         associate (x => every(i:i))
            if (verify(x, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") /= 0) return
            if (count_of(given%a, x) > 1 .or. count_of(given%b, x) > 1 .or. count_of(given%c, x) > 1) return
            if (count_of(every, x) /= 2) return
         end associate
      end do
      valid = .true.
   end function spec_is_valid

   pure function count_of(string, letter) result(n)
      !! Occurrences of `letter` in `string`
      character(len=*), intent(in) :: string
      character, intent(in) :: letter
      integer :: n, i

      n = 0
      do i = 1, len(string)
         if (string(i:i) == letter) n = n + 1
      end do
   end function count_of

   pure function common_letters(order, other) result(found)
      !! The letters of `order` that are also in `other`, in `order`'s order
      character(len=*), intent(in) :: order, other
      character(len=MAX_RANK) :: found
      integer :: i, n

      found = ""
      n = 0
      do i = 1, len_trim(order)
         if (index(other, order(i:i)) > 0) then
            n = n + 1
            found(n:n) = order(i:i)
         end if
      end do
   end function common_letters

   pure function extent_of(letters, names, size_of) result(n)
      !! The product of the extents of `letters` (blank padded; 1 when empty)
      character(len=*), intent(in) :: letters, names
      integer, intent(in) :: size_of(:)
      integer(int64) :: n
      integer :: i

      n = 1
      do i = 1, len_trim(letters)
         n = n*size_of(index(names, letters(i:i)))
      end do
   end function extent_of

   pure function shape_of(letters, names, size_of) result(extent)
      !! The extents of a tensor with the index letters `letters`
      character(len=*), intent(in) :: letters, names
      integer, intent(in) :: size_of(:)
      integer :: extent(len(letters))
      integer :: i

      do i = 1, len(letters)
         extent(i) = size_of(index(names, letters(i:i)))
      end do
   end function shape_of

   function choose_layout(given, names, size_of, extra_work_limit) result(best)
      !! The layout that reorders the fewest elements, then the fewest arrays; the batched one instead when
      !! it reorders fewer still and either the GEMM is cheap enough or the other layout would cost more than
      !! `extra_work_limit` doubles over it (module header). The candidates are the 32 choices the bits of
      !! `choice` spell, the upper 16 of them batched.
      type(letters_t), intent(in) :: given
      character(len=*), intent(in) :: names
      integer, intent(in) :: size_of(:)
      integer(int64), intent(in) :: extra_work_limit
      type(layout_t) :: best
      type(layout_t) :: trial, batched
      integer(int64) :: elements, arrays, best_elements, best_arrays, batched_elements, flops
      integer :: choice, i

      best_elements = huge(best_elements)
      best_arrays = 0
      batched_elements = huge(batched_elements)
      do choice = 0, 31
         if (.not. candidate_layout(choice, given, trial)) cycle
         elements = 0
         arrays = 0
         if (trial%a /= given%a) then
            elements = elements + extent_of(given%a, names, size_of)
            arrays = arrays + 1
         end if
         if (trial%b /= given%b) then
            elements = elements + extent_of(given%b, names, size_of)
            arrays = arrays + 1
         end if
         if (trial%c /= given%c) then
            elements = elements + extent_of(given%c, names, size_of)
            arrays = arrays + 1
         end if
         if (trial%batch /= "") then
            if (elements < batched_elements) then
               batched_elements = elements
               batched = trial
            end if
         else if (elements < best_elements .or. (elements == best_elements .and. arrays < best_arrays)) then
            best_elements = elements
            best_arrays = arrays
            best = trial
         end if
      end do

      if (batched_elements >= best_elements) return
      flops = 2
      do i = 1, len_trim(names)
         flops = flops*size_of(i)
      end do
      if (flops <= BATCH_FLOPS_PER_ELEMENT*(best_elements - batched_elements) .or. &
          best_elements - batched_elements > extra_work_limit) best = batched
   end function choose_layout

   function candidate_layout(choice, given, trial) result(possible)
      !! Candidate number `choice` (bits: contracted in B's order, A's free letters in C's order, B's
      !! likewise, B's free letters first, batched). False when the batch it asks for does not exist.
      integer, intent(in) :: choice
      type(letters_t), intent(in) :: given
      type(layout_t), intent(out) :: trial
      logical :: possible
      character(len=MAX_RANK) :: a, b, c
      integer :: na, nb, nc
      character :: batch

      possible = .false.
      a = given%a
      b = given%b
      c = given%c
      na = len_trim(a)
      nb = len_trim(b)
      nc = len_trim(c)
      batch = ""
      if (btest(choice, 4)) then
         batch = c(nc:nc)
         if (a(na:na) == batch) then
            a(na:na) = ""
         else if (b(nb:nb) == batch) then
            b(nb:nb) = ""
         else
            return
         end if
         c(nc:nc) = ""
      end if

      if (btest(choice, 0)) then
         trial%contracted = common_letters(b, a)
      else
         trial%contracted = common_letters(a, b)
      end if
      if (btest(choice, 1)) then
         trial%free_a = common_letters(c, a)
      else
         trial%free_a = common_letters(a, c)
      end if
      if (btest(choice, 2)) then
         trial%free_b = common_letters(c, b)
      else
         trial%free_b = common_letters(b, c)
      end if
      ! a batch of dot products or of rank-1 updates is worth less than one reordering of a small array
      if (batch /= "" .and. (trial%free_a == "" .or. trial%free_b == "")) return
      trial%swapped = btest(choice, 3)
      trial%batch = batch

      ! an operand is read as it is if it already is (free)(contracted) or (contracted)(free)
      trial%a = gemm_order(a, trial%free_a, trial%contracted)
      trial%b = gemm_order(b, trial%contracted, trial%free_b)
      if (trial%swapped) then
         trial%c = trim(trial%free_b)//trim(trial%free_a)
      else
         trial%c = trim(trial%free_a)//trim(trial%free_b)
      end if
      ! put the batched letter back, last, on the array it came off
      if (batch /= "") then
         if (len_trim(a) < na) trial%a = trim(trial%a)//batch
         if (len_trim(b) < nb) trial%b = trim(trial%b)//batch
         trial%c = trim(trial%c)//batch
      end if
      possible = .true.
   end function candidate_layout

   pure function gemm_order(stored, first, second) result(order)
      !! `stored` if it is already `first` then `second` or the reverse (the GEMM's transpose covers
      !! that), else `first` then `second`
      character(len=*), intent(in) :: stored, first, second
      character(len=MAX_RANK) :: order

      order = trim(first)//trim(second)
      if (stored == trim(second)//trim(first)) order = stored
   end function gemm_order

   function gemm_of(layout, names, size_of) result(g)
      !! The GEMM, or the strided batch, that reads the arrays in `layout`'s orders
      type(layout_t), intent(in) :: layout
      character(len=*), intent(in) :: names
      integer, intent(in) :: size_of(:)
      type(trc_gemm_batch_t) :: g
      character(len=MAX_RANK) :: left, right, left_free, right_free
      character :: transa, transb
      integer :: m, n, k, count
      integer(int64) :: sa, sb

      if (layout%swapped) then
         left = layout%b
         right = layout%a
         left_free = layout%free_b
         right_free = layout%free_a
      else
         left = layout%a
         right = layout%b
         left_free = layout%free_a
         right_free = layout%free_b
      end if
      m = int(extent_of(left_free, names, size_of))
      n = int(extent_of(right_free, names, size_of))
      k = int(extent_of(layout%contracted, names, size_of))
      count = 1
      if (layout%batch /= "") count = size_of(index(names, layout%batch))

      transa = "N"
      if (index(left, trim(left_free)//trim(layout%contracted)) /= 1) transa = "T"
      transb = "N"
      if (index(right, trim(layout%contracted)//trim(right_free)) /= 1) transb = "T"

      ! only the operand that carries the batched letter advances from item to item
      sa = 0
      sb = 0
      if (layout%batch /= "") then
         if (index(left, layout%batch) > 0) sa = int(m, int64)*k
         if (index(right, layout%batch) > 0) sb = int(k, int64)*n
      end if
      g = trc_gemm_batch(m, n, k, count, sa=sa, sb=sb, transa=transa, transb=transb)
   end function gemm_of

end module trc_contract
