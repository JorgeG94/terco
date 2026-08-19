!
! Hermite index tables.
!
! The R recurrence is naturally written over (t,u,v) with a branch on which
! index to decrement.  On a device that is three-deep nesting plus divergence,
! for an operation that is really just "walk a list".  Flattening (t,u,v) to a
! single Hermite index and precomputing the decrement targets turns the whole
! recurrence into
!
!     do h = 2, NHERM(deg)
!        r(h,n) = pq(dir(h)) * r(hm1(h), n+1) + cfac(h) * r(hm2(h), n+1)
!     end do
!
! -- one flat loop, no branches, no divergence.  That is the plan's claim that
! MMD is loop nests plus tables rather than generated straight-line code, made
! concrete.
!
! ORDERING (load-bearing)
! -----------------------
! Degree-major: all (t,u,v) with t+u+v = 0 first, then all with 1, and so on.
! Two things fall out and both are used:
!
!   * "every index of degree <= d" is exactly the first NHERM(d) entries, so
!     the recurrence's shrinking bounds are a single integer, not a test;
!   * a decrement always lands on a STRICTLY smaller index, so the level-n
!     sweep can read level n+1 with no ordering hazard.
!
module trc_tables
   use trc_boys, only: dp
   implicit none
   private

   ! Highest angular momentum per shell this build supports.
   !
   ! A BUILD-TIME knob, not a runtime one, and deliberately so: it sizes the
   ! contraction's working vectors, which the compiler keeps in registers, and
   ! registers are what caps occupancy.  An s,p-only build (6-31G and friends)
   ! wants LMAX=1 and gets much smaller vectors for it.  Shipping one binary
   ! per angular-momentum ceiling is what gamess-libERI does too, for the same
   ! reason.
#ifndef TRC_LMAX
#define TRC_LMAX 2
#endif
   integer, parameter, public :: LMAX = TRC_LMAX
   ! Highest total Hermite degree a 4-centre quartet can reach.
   integer, parameter, public :: LTOT = 4*LMAX
   ! Number of Hermite functions up to degree LTOT.
   integer, parameter, public :: NHERM_MAX = (LTOT + 1)*(LTOT + 2)*(LTOT + 3)/6
   ! ... and up to the degree ONE PAIR can reach, which is half of it.  The
   ! contraction's working vectors are indexed by a pair's Hermite degree, not
   ! the quartet's, and sizing them at NHERM_MAX instead cost 2.6 kB of stack
   ! frame per thread and held occupancy near 22%.
   integer, parameter, public :: LPAIR = 2*LMAX
   integer, parameter, public :: NHERM_PAIR = (LPAIR + 1)*(LPAIR + 2)*(LPAIR + 3)/6

   public :: nherm, tables_init, herm_index

   integer,  public :: hidx(0:LTOT, 0:LTOT, 0:LTOT)  !! (t,u,v) -> flat, 0 if degree > LTOT
   integer,  public :: hdir(NHERM_MAX)               !! 1/2/3 = which of x,y,z was decremented
   integer,  public :: hm1(NHERM_MAX)                !! flat index of the once-decremented parent
   integer,  public :: hm2(NHERM_MAX)                !! twice-decremented, or 0
   real(dp), public :: hcf(NHERM_MAX)                !! the (n-1) coefficient on the hm2 term
   integer,  public :: nherm_of(0:LTOT)              !! cumulative count by degree

   ! Decode: flat Hermite index -> its (t,u,v) and its (-1)^(t+u+v).
   integer,  public :: ht(NHERM_MAX), hu(NHERM_MAX), hv(NHERM_MAX)
   real(dp), public :: hsgn(NHERM_MAX)

   ! hshift(h1,h2) = hidx(t1+t2, u1+u2, v1+v2), or 0 when the sum overflows
   ! LTOT.  This is what makes the two-GEMM contraction a pure gather: the
   ! inner loop needs no index arithmetic at all, just a table load.
   integer,  public :: hshift(NHERM_MAX, NHERM_MAX)

   !$acc declare create(hidx, hdir, hm1, hm2, hcf, nherm_of, &
   !$acc                ht, hu, hv, hsgn, hshift)

contains

   pure integer function nherm(l)
      integer, intent(in) :: l
      nherm = (l + 1)*(l + 2)*(l + 3)/6
   end function nherm

   pure integer function herm_index(t, u, v)
      !$acc routine seq
      integer, intent(in) :: t, u, v
      herm_index = hidx(t, u, v)
   end function herm_index

   !
   ! Build the tables and push them to the device.  Call once, before any
   ! kernel.  Host-side loops -- this is setup, not work.
   !
   subroutine tables_init()
      integer :: n, t, u, v, h

      hidx = 0; hdir = 0; hm1 = 0; hm2 = 0; hcf = 0.0_dp
      ht = 0; hu = 0; hv = 0; hsgn = 0.0_dp; hshift = 0

      h = 0
      do n = 0, LTOT
         do t = n, 0, -1
            do u = n - t, 0, -1
               v = n - t - u
               h = h + 1
               hidx(t, u, v) = h
               ht(h) = t; hu(h) = u; hv(h) = v
               hsgn(h) = 1.0_dp
               if (mod(n, 2) == 1) hsgn(h) = -1.0_dp

               if (n == 0) cycle

               ! Decrement whichever index is non-zero, x first.  The choice
               ! only has to be consistent, not clever.
               if (t > 0) then
                  hdir(h) = 1
                  hm1(h) = hidx(t - 1, u, v)
                  if (t > 1) then
                     hm2(h) = hidx(t - 2, u, v)
                     hcf(h) = real(t - 1, dp)
                  end if
               else if (u > 0) then
                  hdir(h) = 2
                  hm1(h) = hidx(t, u - 1, v)
                  if (u > 1) then
                     hm2(h) = hidx(t, u - 2, v)
                     hcf(h) = real(u - 1, dp)
                  end if
               else
                  hdir(h) = 3
                  hm1(h) = hidx(t, u, v - 1)
                  if (v > 1) then
                     hm2(h) = hidx(t, u, v - 2)
                     hcf(h) = real(v - 1, dp)
                  end if
               end if
            end do
         end do
         nherm_of(n) = h
      end do

      ! Shift table, once every hidx entry exists.
      block
         integer :: h1, h2
         do h2 = 1, NHERM_MAX
            do h1 = 1, NHERM_MAX
               if (ht(h1) + ht(h2) <= LTOT .and. hu(h1) + hu(h2) <= LTOT .and. &
                   hv(h1) + hv(h2) <= LTOT) then
                  if (ht(h1) + ht(h2) + hu(h1) + hu(h2) + hv(h1) + hv(h2) <= LTOT) then
                     hshift(h1, h2) = hidx(ht(h1) + ht(h2), hu(h1) + hu(h2), hv(h1) + hv(h2))
                  end if
               end if
            end do
         end do
      end block

      !$acc update device(hidx, hdir, hm1, hm2, hcf, nherm_of, &
      !$acc               ht, hu, hv, hsgn, hshift)
   end subroutine tables_init

end module trc_tables
