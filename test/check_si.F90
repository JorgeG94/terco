program check_si
   !! trc_si's prefix choice: each count printed with the prefix that keeps
   !! its number in [1, 1000), unprefixed below 1, exa-scaled past it, never
   !! rounding up to "1000.00", and right-justified when a width is given.
   !! The cases are the RI-MP2 table's own (2453056.2 GFLOP and friends).
   use pic_types, only: dp
   use trc_si, only: si_value_t
   implicit none
   logical :: ok

   ok = .true.
   call expect(si_value_t(2453056.2e9_dp, "FLOP"), "2.45 PFLOP")
   call expect(si_value_t(53270.3e9_dp, "FLOP/s"), "53.27 TFLOP/s")
   call expect(si_value_t(866.2e9_dp, "FLOP/s"), "866.20 GFLOP/s")
   call expect(si_value_t(124.7e9_dp, "FLOP"), "124.70 GFLOP")
   call expect(si_value_t(62.0e9_dp, "B/s"), "62.00 GB/s")
   call expect(si_value_t(1000.0_dp, "FLOP"), "1.00 kFLOP")
   call expect(si_value_t(999.996e12_dp, "FLOP"), "1.00 PFLOP")
   call expect(si_value_t(0.5_dp, "FLOP"), "0.50 FLOP")
   call expect(si_value_t(3.0e21_dp, "FLOP"), "3000.00 EFLOP")
   if (si_value_t(2453056.2e9_dp, "FLOP")%to_string(15) /= "     2.45 PFLOP") then
      print '(a)', "check_si: FAILED: width 15 gave '"//si_value_t(2453056.2e9_dp, "FLOP")%to_string(15)//"'"
      ok = .false.
   end if

   if (.not. ok) then
      print '(a)', "check_si: FAIL"
      stop 1
   end if
   print '(a)', "check_si: PASS"

contains

   subroutine expect(q, want)
      type(si_value_t), intent(in) :: q
      character(len=*), intent(in) :: want
      if (q%to_string() /= want) then
         print '(a)', "check_si: FAILED: got '"//q%to_string()//"', want '"//want//"'"
         ok = .false.
      end if
   end subroutine expect

end program check_si
