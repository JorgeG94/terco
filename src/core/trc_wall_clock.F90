!
! Seconds since the first call, read from a pic timer. The step table, RI-MP2
! and the XC driver all subtract two readings, so any fixed origin will do and
! a timer_type (which only measures an interval) can supply it.
!
module trc_wall_clock
   use pic_types, only: dp
   use pic_timer, only: timer_type
   implicit none
   private
   public :: trc_wall

   type(timer_type), save :: epoch
   logical, save :: epoch_started = .false.

contains

   function trc_wall() result(t)
      real(dp) :: t
      if (.not. epoch_started) then
         call epoch%start()
         epoch_started = .true.
      end if
      t = epoch%get_elapsed_time()
   end function trc_wall

end module trc_wall_clock
