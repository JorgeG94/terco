!
! How the particle-particle ladder of an RI-CCSD iteration is cut over the ranks of a run (one rank per device).
!
! THE UNIT is the virtual index d of the half form (trc_rccsd_ladder): the passes of a rank run over a contiguous
! range d_lo..d_hi, and the ladder is linear in that range, so the partial arrays of the ranks add up to the whole
! one. A virtual d owns the columns (c,d) with c >= d, nv - d + 1 of them, and every column costs the same
! (`rccsd_ladder_column_flops`: the build of C and the contraction with tau), so a range costs that times its columns.
!
! THE SCHEDULE this balances (trc_rccsd_update): the first rank runs its share of the ladder AND all of its other work
! while the other ranks run theirs, and takes their partial ladders last, so the model is exact in what it counts:
! the first rank's time is (its ladder + O) / rate, and a follower's its ladder / rate.
!
! THE WEIGHTS. The first rank does everything that is not the ladder as well (tau, the one-particle intermediates, W,
! the T1 and T2 equations with the rings, the t1 t1 correction of the dressed ladder, DIIS), so it gets the smaller
! share. With L the flops of the whole ladder and O the flops of the first rank's other work
! (`rccsd_share_other_flops`: the counted flops of the stages, the same functions the step table reports), n ranks,
! the aim is that every rank has the same flops in all:
!
!     the first rank's ladder  s0 = max(0, (L + O) / n - O)         every other rank (L - s0) / (n - 1)
!
! and when O alone is above an equal share (L + O) / n, s0 is 0 and the first rank gets an EMPTY range. The cut of
! a boundary is then to the nearest whole d (the cumulative columns k nv - k (k - 1) / 2 after k virtuals).
! bt, the t1-dressed factor, is built by every rank and cancels. The rule counts flops, not seconds: the first
! rank's stages run at other rates than the ladder's GEMMs (the kernels in them have a time and no count), so the
! result is a starting point, not a fit; the log shows each rank's ladder seconds to judge it.
! It is a pure function of (no, nv, naux, ranks): the same on every rank, no run-time input, and independent of
! d_batch and of the other ladder options (the half form is the only one that can be cut).
!
! NUMBERS (counted flops, TFLOP; "ladder" is the flops of the half-form passes of the range, "in all" adds the first
! rank's other work O; ranges d = lo..hi, "-" is an empty range):
!
!     adenine cc-pVDZ, no 25, nv 130, naux 630: L 0.361, O 0.518, an equal share (L + O) / n = 0.440 at n = 2
!        2 ranks   rank 0  -          0 % of the ladder, 0.518 in all   rank 1  d 1..130  100 %, 0.361
!        4 ranks   rank 0  -          0 %, 0.518                        rank 1  d 1..24    33.4 %, 0.121
!                  rank 2  d 25..55   33.1 %, 0.120                     rank 3  d 56..130  33.5 %, 0.121
!     (H2O)9 cc-pVTZ, no 36, nv 477, naux 1269: L 133.07, O 66.91, an equal share 100.0 at n = 2
!        2 ranks   rank 0  d 1..64    25.0 % of the ladder, 100.19 in all   rank 1  d 65..477  75.0 %, 99.79
!        4 ranks   rank 0  -          0 %, 66.91                        rank 1  d 1..88    33.5 %, 44.53
!                  rank 2  d 89..202  33.3 %, 44.24                     rank 3  d 203..477 33.3 %, 44.30
!     At adenine the ladder is the smaller part of an iteration (0.36 of 0.88 TFLOP), so the first rank keeps no pass
!     at any rank count; at (H2O)9 it is the larger (133 of 200), and two ranks cut it 1 : 3.
!
module trc_rccsd_share
   use pic_types, only: dp, int64
   use trc_rccsd_inter, only: rccsd_inter_flops
   use trc_rccsd_w, only: rccsd_w_flops
   use trc_rccsd_ladder, only: rccsd_ladder_column_flops, rccsd_ladder_correction_flops
   use trc_rccsd_update, only: rccsd_update_t1_flops, rccsd_update_t2_flops, rccsd_update_ring_flops
   implicit none
   private

   public :: rccsd_share_ranges, rccsd_share_other_flops, rccsd_share_columns, rccsd_share_ladder_flops

contains

   pure function rccsd_share_columns(nv, k) result(n)
      !! Columns (c,d), c >= d, of the first k virtuals d: sum over d = 1..k of nv - d + 1
      integer, intent(in) :: nv, k
      integer(int64) :: n

      n = int(k, int64)*nv - (int(k, int64)*(k - 1))/2
   end function rccsd_share_columns

   pure function rccsd_share_ladder_flops(no, nv, naux, d_lo, d_hi) result(flops)
      !! Flops of the half-form passes over d = d_lo..d_hi (0 for an empty range)
      integer, intent(in) :: no, nv, naux, d_lo, d_hi
      real(dp) :: flops

      flops = 0.0_dp
      if (d_hi < d_lo) return
      flops = rccsd_ladder_column_flops(no, nv, naux)* &
              real(rccsd_share_columns(nv, d_hi) - rccsd_share_columns(nv, d_lo - 1), dp)
   end function rccsd_share_ladder_flops

   pure function rccsd_share_other_flops(no, nv, naux) result(flops)
      !! The counted flops of what only the first rank does in an iteration: tau and the one-particle intermediates, W,
      !! the T1 and T2 equations, the three rings and the t1 t1 correction of the dressed ladder
      integer, intent(in) :: no, nv, naux
      real(dp) :: flops

      flops = real(rccsd_inter_flops(no, nv, naux), dp) + real(rccsd_w_flops(no, nv, naux), dp) + &
              real(rccsd_update_t1_flops(no, nv, naux), dp) + real(rccsd_update_t2_flops(no, nv, naux), dp) + &
              real(rccsd_update_ring_flops(no, nv), dp) + real(rccsd_ladder_correction_flops(no, nv), dp)
   end function rccsd_share_other_flops

   pure subroutine rccsd_share_ranges(no, nv, naux, nranks, d_lo, d_hi)
      !! The d range of each rank, ranks 0..nranks-1 at index 1..nranks: contiguous, in order, disjoint, covering 1..nv.
      !! An empty range has d_hi = d_lo - 1. Module header for the rule.
      integer, intent(in) :: no, nv, naux, nranks
      integer, intent(out) :: d_lo(nranks), d_hi(nranks)

      real(dp) :: unit, total, other, first, target_cols, goal
      integer :: r, previous, k

      if (nranks <= 1) then
         d_lo = 1
         d_hi = nv
         return
      end if
      unit = rccsd_ladder_column_flops(no, nv, naux)
      total = unit*real(rccsd_share_columns(nv, nv), dp)
      other = rccsd_share_other_flops(no, nv, naux)
      first = max(0.0_dp, (total + other)/real(nranks, dp) - other)
      first = min(first, total)
      previous = 0
      do r = 1, nranks - 1
         ! the cumulative column count at the end of rank r-1 (r = 1 is the first rank)
         target_cols = first/unit + real(r - 1, dp)*(real(rccsd_share_columns(nv, nv), dp) - first/unit)/ &
                       real(nranks - 1, dp)
         k = previous
         do while (k < nv)
            if (real(rccsd_share_columns(nv, k), dp) >= target_cols) exit
            k = k + 1
         end do
         if (k > previous) then
            goal = real(rccsd_share_columns(nv, k), dp) - target_cols
            if (target_cols - real(rccsd_share_columns(nv, k - 1), dp) < goal) k = k - 1
         end if
         d_lo(r) = previous + 1
         d_hi(r) = k
         previous = k
      end do
      d_lo(nranks) = previous + 1
      d_hi(nranks) = nv
   end subroutine rccsd_share_ranges

end module trc_rccsd_share
