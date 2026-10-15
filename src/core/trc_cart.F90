!
! Cartesian index tables for the Obara-Saika / Head-Gordon-Pople path.
!
! Same idea as trc_tables does for Hermite indices: flatten the (nx,ny,nz)
! triple to one index, precompute every neighbour the recurrences ask for, and
! leave the kernels as flat table-driven loops with no index arithmetic.
!
! ORDERING is angular-momentum-major, and within one L it is libcint's own
! order (x descending outer, y descending inner).  Two things follow:
!
!   * "every index with angular momentum <= L" is the first NCUM_OF(L) entries,
!     which is what the VRR's shrinking bounds want;
!   * the local component index within a shell is the global index minus
!     NCUM_OF(L-1), so mapping a finished block onto libcint's buffer layout is
!     a subtraction rather than a search.
!
module trc_cart
   use trc_boys, only: dp
   use trc_tables, only: LMAX
   implicit none
   private

   ! Highest angular momentum the VRR has to reach: a bra pair can carry
   ! la+lb, so 2*LMAX.
   integer, parameter, public :: LCMAX = 2*LMAX
   integer, parameter, public :: NCUM = (LCMAX + 1)*(LCMAX + 2)*(LCMAX + 3)/6

   public :: cart_init, ncart_of

   integer, public :: cidx(0:LCMAX, 0:LCMAX, 0:LCMAX)  !! (nx,ny,nz) -> flat, 0 if L > LCMAX
   integer, public :: cnx(NCUM), cny(NCUM), cnz(NCUM)  !! flat -> powers
   integer, public :: cll(NCUM)                        !! flat -> its angular momentum
   integer, public :: cdir(NCUM)                       !! VRR direction: first non-zero power
   integer, public :: cdn1(NCUM)                       !! index of a - 1_dir, 0 if none
   integer, public :: cdn2(NCUM)                       !! index of a - 2_dir, 0 if none
   real(dp), public :: cf2(NCUM)                       !! (a_dir - 1), the cdn2 coefficient
   integer, public :: cup(NCUM, 3)                     !! index of a + 1_d, 0 if beyond LCMAX
   integer, public :: ncum_of(0:LCMAX)                 !! cumulative count through L

   !$acc declare create(cidx, cnx, cny, cnz, cll, cdir, cdn1, cdn2, cf2, cup, ncum_of)

contains

   pure integer function ncart_of(l)
      !$acc routine seq
      integer, intent(in) :: l
      ncart_of = (l + 1)*(l + 2)/2
   end function ncart_of

   subroutine cart_init()
      integer :: l, lx, ly, lz, h, d, nx, ny, nz

      cidx = 0; cnx = 0; cny = 0; cnz = 0; cll = 0
      cdir = 0; cdn1 = 0; cdn2 = 0; cf2 = 0.0_dp; cup = 0

      h = 0
      do l = 0, LCMAX
         do lx = l, 0, -1
            do ly = l - lx, 0, -1
               lz = l - lx - ly
               h = h + 1
               cidx(lx, ly, lz) = h
               cnx(h) = lx; cny(h) = ly; cnz(h) = lz; cll(h) = l
            end do
         end do
         ncum_of(l) = h
      end do

      ! Neighbours, now that every index exists.
      do h = 1, NCUM
         lx = cnx(h); ly = cny(h); lz = cnz(h)

         ! VRR decrement direction: first non-zero, x then y then z.  The
         ! choice only has to be consistent.
         if (lx > 0) then
            cdir(h) = 1
            cdn1(h) = cidx(lx - 1, ly, lz)
            if (lx > 1) then
               cdn2(h) = cidx(lx - 2, ly, lz)
               cf2(h) = real(lx - 1, dp)
            end if
         else if (ly > 0) then
            cdir(h) = 2
            cdn1(h) = cidx(lx, ly - 1, lz)
            if (ly > 1) then
               cdn2(h) = cidx(lx, ly - 2, lz)
               cf2(h) = real(ly - 1, dp)
            end if
         else if (lz > 0) then
            cdir(h) = 3
            cdn1(h) = cidx(lx, ly, lz - 1)
            if (lz > 1) then
               cdn2(h) = cidx(lx, ly, lz - 2)
               cf2(h) = real(lz - 1, dp)
            end if
         end if

         ! Increment in each direction -- the HRR walks upward.
         do d = 1, 3
            nx = lx; ny = ly; nz = lz
            if (d == 1) nx = nx + 1
            if (d == 2) ny = ny + 1
            if (d == 3) nz = nz + 1
            if (nx + ny + nz <= LCMAX) cup(h, d) = cidx(nx, ny, nz)
         end do
      end do

      !$acc update device(cidx, cnx, cny, cnz, cll, cdir, cdn1, cdn2, cf2, cup, ncum_of)
   end subroutine cart_init

end module trc_cart
