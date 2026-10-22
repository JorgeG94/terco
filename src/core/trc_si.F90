!
! A quantity in a base unit, printed with the SI prefix that keeps its
! number between 1 and 1000: 2.45e15 FLOP prints as "2.45 PFLOP", 5.3e13
! FLOP/s as "53.27 TFLOP/s". The step tables used to print every count in
! fixed giga units, which on an RI-MP2 run read "2453056.2 GFLOP".
!
module trc_si
   use pic_types, only: dp
   use pic_format, only: to_string_fixed
   implicit none
   private

   public :: si_value_t

   integer, parameter :: NPREFIX = 6
   !! Prefixes past none: kilo through exa.
   character(len=1), parameter :: PREFIXES(0:NPREFIX) = [" ", "k", "M", "G", "T", "P", "E"]

   type :: si_value_t
      real(dp) :: value = 0.0_dp
      !! In the base unit (FLOP, not GFLOP).
      character(len=:), allocatable :: unit
      !! The base unit, e.g. "FLOP", "FLOP/s", "B/s".
   contains
      procedure :: to_string => si_value_to_string
   end type si_value_t

contains

   function si_value_to_string(this, width) result(s)
      !! The value scaled to [1, 1000) with two decimals, then the prefixed
      !! unit; right-justified in `width` characters when given. Below 1 in
      !! the base unit, or past exa, it stays unprefixed or exa-scaled.
      class(si_value_t), intent(in) :: this
      integer, intent(in), optional :: width
      character(len=:), allocatable :: s
      real(dp) :: v
      integer :: p
      v = abs(this%value)
      p = 0
      do while (v >= 999.995_dp .and. p < NPREFIX)   ! 999.996 would round to "1000.00"
         v = v/1000.0_dp
         p = p + 1
      end do
      s = to_string_fixed(sign(v, this%value), 2)//" "//trim(PREFIXES(p))//this%unit
      if (present(width)) then
         if (len(s) < width) s = repeat(" ", width - len(s))//s
      end if
   end function si_value_to_string

end module trc_si
