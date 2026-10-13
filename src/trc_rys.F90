!
! Rys roots and weights, by table lookup, for n = 1..5.
!
! WHAT THIS IS
! ------------
! The n-point Gauss rule for the weight exp(-X u) on u in [0,1] against the
! Jacobi measure u^(-1/2), which is what a Rys quadrature is. A root comes
! back as u_i; the kernel wants u/(1-u) and applies that itself.
!
! n = floor(L_total/2) + 1 where L_total = la+lb+lc+ld, so LMAX=2 (L_total up
! to 8) needs n up to 5.
!
! WHY A TABLE, AND WHY NO CONSTRUCTION IS KEPT HERE
! --------------------------------------------------
! GAMESS's own construction -- RYSDS/RYSGW: discretised Stieltjes on a fixed
! 35-point base rule, then Golub-Welsch by QL -- is correct in double
! precision (checked against an independent mpmath construction of the same
! rule to 1.3e-15 for X in [0.05, 50], while this file was being written), but
! it is unusable on the device this library targets: compiled for it, that
! routine alone reports REG:255, the entire register budget, before a single
! 2D integral exists. Add 35 exponentials and a QL whose trip count is data
! dependent -- warp divergence -- and tabulating wins outright. Only the
! table-lookup form is shipped; `scripts/gen_rys_roots.py` is what regenerates
! the tables these five routines read, and it is checked against mpmath
! directly, not against the construction this replaces.
!
! NO ASSUMED-SHAPE ARRAYS. Every array here is explicit-shape or a named
! constant from the table modules: an assumed-shape dummy made `boys_eval`
! 2.4x slower before it was removed, and the same cost applies to anything
! called from `!$acc routine seq`.
!
module trc_rys
   use trc_boys, only: dp
   use trc_rys_table1, only: RYS_C1, RYS_ASY_R1, RYS_ASY_W1, RYS_NROOTS1, &
                             RYS_DEGREE1, RYS_NINT1, RYS_WIDTH1, RYS_XMAX1
   use trc_rys_table2, only: RYS_C2, RYS_ASY_R2, RYS_ASY_W2, RYS_NROOTS2, &
                             RYS_DEGREE2, RYS_NINT2, RYS_WIDTH2, RYS_XMAX2
   use trc_rys_table3, only: RYS_C3, RYS_ASY_R3, RYS_ASY_W3, RYS_NROOTS3, &
                             RYS_DEGREE3, RYS_NINT3, RYS_WIDTH3, RYS_XMAX3
   use trc_rys_table4, only: RYS_C4, RYS_ASY_R4, RYS_ASY_W4, RYS_NROOTS4, &
                             RYS_DEGREE4, RYS_NINT4, RYS_WIDTH4, RYS_XMAX4
   use trc_rys_table5, only: RYS_C5, RYS_ASY_R5, RYS_ASY_W5, RYS_NROOTS5, &
                             RYS_DEGREE5, RYS_NINT5, RYS_WIDTH5, RYS_XMAX5
   implicit none
   private
   public :: trc_rys_roots1_tab, trc_rys_roots2_tab, trc_rys_roots3_tab
   public :: trc_rys_roots4_tab, trc_rys_roots5_tab

contains

   pure subroutine trc_rys_roots1_tab(x, u, w)
      !! One Rys root and weight by table lookup.
      !$acc routine seq
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: u(1), w(1)

      real(dp) :: t, t2, d0, d1, tmp, fr, fw
      integer  :: iv, s, j, base
      real(dp) :: xc

      if (x >= RYS_XMAX1) then
         fr = 1.0_dp/x
         fw = sqrt(fr)
         do s = 1, RYS_NROOTS1
            u(s) = RYS_ASY_R1(s)*fr
            w(s) = RYS_ASY_W1(s)*fw
         end do
         return
      end if

      xc = max(x, 0.0_dp)
      iv = min(int(xc/RYS_WIDTH1) + 1, RYS_NINT1)
      ! map [lo, hi] onto [-1, 1]; lo = (iv-1)*WIDTH, hi = lo + WIDTH
      t = (2.0_dp*xc - (2*iv - 1)*RYS_WIDTH1)/RYS_WIDTH1
      t2 = 2.0_dp*t

      do s = 1, 2*RYS_NROOTS1
         base = ((iv - 1)*2*RYS_NROOTS1 + (s - 1))*(RYS_DEGREE1 + 1)
         d0 = 0.0_dp
         d1 = 0.0_dp
         do j = RYS_DEGREE1, 1, -1
            tmp = d0
            d0 = t2*d0 - d1 + RYS_C1(base + j + 1)
            d1 = tmp
         end do
         tmp = t*d0 - d1 + RYS_C1(base + 1)
         if (s <= RYS_NROOTS1) then
            u(s) = tmp
         else
            w(s - RYS_NROOTS1) = tmp
         end if
      end do
   end subroutine trc_rys_roots1_tab

   pure subroutine trc_rys_roots2_tab(x, u, w)
      !! Two Rys roots and weights by table lookup.
      !$acc routine seq
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: u(2), w(2)

      real(dp) :: t, t2, d0, d1, tmp, fr, fw
      integer  :: iv, s, j, base
      real(dp) :: xc

      if (x >= RYS_XMAX2) then
         fr = 1.0_dp/x
         fw = sqrt(fr)
         do s = 1, RYS_NROOTS2
            u(s) = RYS_ASY_R2(s)*fr
            w(s) = RYS_ASY_W2(s)*fw
         end do
         return
      end if

      xc = max(x, 0.0_dp)
      iv = min(int(xc/RYS_WIDTH2) + 1, RYS_NINT2)
      t = (2.0_dp*xc - (2*iv - 1)*RYS_WIDTH2)/RYS_WIDTH2
      t2 = 2.0_dp*t

      do s = 1, 2*RYS_NROOTS2
         base = ((iv - 1)*2*RYS_NROOTS2 + (s - 1))*(RYS_DEGREE2 + 1)
         d0 = 0.0_dp
         d1 = 0.0_dp
         do j = RYS_DEGREE2, 1, -1
            tmp = d0
            d0 = t2*d0 - d1 + RYS_C2(base + j + 1)
            d1 = tmp
         end do
         tmp = t*d0 - d1 + RYS_C2(base + 1)
         if (s <= RYS_NROOTS2) then
            u(s) = tmp
         else
            w(s - RYS_NROOTS2) = tmp
         end if
      end do
   end subroutine trc_rys_roots2_tab

   pure subroutine trc_rys_roots3_tab(x, u, w)
      !! Three Rys roots and weights by table lookup. Covers (pp|pp).
      !$acc routine seq
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: u(3), w(3)

      real(dp) :: t, t2, d0, d1, tmp, fr, fw
      integer  :: iv, s, j, base
      real(dp) :: xc

      if (x >= RYS_XMAX3) then
         fr = 1.0_dp/x
         fw = sqrt(fr)
         do s = 1, RYS_NROOTS3
            u(s) = RYS_ASY_R3(s)*fr
            w(s) = RYS_ASY_W3(s)*fw
         end do
         return
      end if

      xc = max(x, 0.0_dp)
      iv = min(int(xc/RYS_WIDTH3) + 1, RYS_NINT3)
      t = (2.0_dp*xc - (2*iv - 1)*RYS_WIDTH3)/RYS_WIDTH3
      t2 = 2.0_dp*t

      do s = 1, 2*RYS_NROOTS3
         base = ((iv - 1)*2*RYS_NROOTS3 + (s - 1))*(RYS_DEGREE3 + 1)
         d0 = 0.0_dp
         d1 = 0.0_dp
         do j = RYS_DEGREE3, 1, -1
            tmp = d0
            d0 = t2*d0 - d1 + RYS_C3(base + j + 1)
            d1 = tmp
         end do
         tmp = t*d0 - d1 + RYS_C3(base + 1)
         if (s <= RYS_NROOTS3) then
            u(s) = tmp
         else
            w(s - RYS_NROOTS3) = tmp
         end if
      end do
   end subroutine trc_rys_roots3_tab

   pure subroutine trc_rys_roots4_tab(x, u, w)
      !! Four Rys roots and weights by table lookup.
      !$acc routine seq
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: u(4), w(4)

      real(dp) :: t, t2, d0, d1, tmp, fr, fw
      integer  :: iv, s, j, base
      real(dp) :: xc

      if (x >= RYS_XMAX4) then
         fr = 1.0_dp/x
         fw = sqrt(fr)
         do s = 1, RYS_NROOTS4
            u(s) = RYS_ASY_R4(s)*fr
            w(s) = RYS_ASY_W4(s)*fw
         end do
         return
      end if

      xc = max(x, 0.0_dp)
      iv = min(int(xc/RYS_WIDTH4) + 1, RYS_NINT4)
      t = (2.0_dp*xc - (2*iv - 1)*RYS_WIDTH4)/RYS_WIDTH4
      t2 = 2.0_dp*t

      do s = 1, 2*RYS_NROOTS4
         base = ((iv - 1)*2*RYS_NROOTS4 + (s - 1))*(RYS_DEGREE4 + 1)
         d0 = 0.0_dp
         d1 = 0.0_dp
         do j = RYS_DEGREE4, 1, -1
            tmp = d0
            d0 = t2*d0 - d1 + RYS_C4(base + j + 1)
            d1 = tmp
         end do
         tmp = t*d0 - d1 + RYS_C4(base + 1)
         if (s <= RYS_NROOTS4) then
            u(s) = tmp
         else
            w(s - RYS_NROOTS4) = tmp
         end if
      end do
   end subroutine trc_rys_roots4_tab

   pure subroutine trc_rys_roots5_tab(x, u, w)
      !! Five Rys roots and weights by table lookup. Covers (dd|dd).
      !$acc routine seq
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: u(5), w(5)

      real(dp) :: t, t2, d0, d1, tmp, fr, fw
      integer  :: iv, s, j, base
      real(dp) :: xc

      if (x >= RYS_XMAX5) then
         fr = 1.0_dp/x
         fw = sqrt(fr)
         do s = 1, RYS_NROOTS5
            u(s) = RYS_ASY_R5(s)*fr
            w(s) = RYS_ASY_W5(s)*fw
         end do
         return
      end if

      xc = max(x, 0.0_dp)
      iv = min(int(xc/RYS_WIDTH5) + 1, RYS_NINT5)
      t = (2.0_dp*xc - (2*iv - 1)*RYS_WIDTH5)/RYS_WIDTH5
      t2 = 2.0_dp*t

      do s = 1, 2*RYS_NROOTS5
         base = ((iv - 1)*2*RYS_NROOTS5 + (s - 1))*(RYS_DEGREE5 + 1)
         d0 = 0.0_dp
         d1 = 0.0_dp
         do j = RYS_DEGREE5, 1, -1
            tmp = d0
            d0 = t2*d0 - d1 + RYS_C5(base + j + 1)
            d1 = tmp
         end do
         tmp = t*d0 - d1 + RYS_C5(base + 1)
         if (s <= RYS_NROOTS5) then
            u(s) = tmp
         else
            w(s - RYS_NROOTS5) = tmp
         end if
      end do
   end subroutine trc_rys_roots5_tab

end module trc_rys
