!
! terco's own numbers for the golden tier ("Two tiers" in
! test/pyscf_refs/README.md).
!
! A check that is compared against pyscf (the anchor, test/pyscf_refs/) also
! records the same quantities here, from this build's own run. ctest sets
! TERCO_GOLDEN_OUT to a directory in the build tree (which must exist), and
! the paired terco/golden/<case> test (test/goldens/compare_goldens.py)
! diffs the output against the committed test/goldens/<case>.gold at a
! tolerance far tighter than the anchor's. With TERCO_GOLDEN_OUT unset (a
! check run by hand) every call here does nothing.
!
! Output, one record per line:
!   # terco golden output v1 <case>
!   <key> <value, es24.16e3>
!   END
! A check that skips itself after golden_begin leaves no values and no END;
! compare_goldens.py reports that as SKIP, so a skip never reads a previous
! run's numbers.
!
module golden_io
   use pic_types, only: dp
   implicit none
   private

   public :: golden_begin, golden_put, golden_fock, golden_end

   integer, save :: gu = 0
   logical, save :: active = .false.

contains

   !
   ! Start (and truncate) this run's output for `case`.
   !
   subroutine golden_begin(case)
      character(len=*), intent(in) :: case
      character(len=4096) :: dir
      integer :: n, st, ios
      call get_environment_variable("TERCO_GOLDEN_OUT", dir, length=n, status=st)
      if (st /= 0 .or. n == 0) return
      open (newunit=gu, file=dir(1:n)//"/"//trim(case)//".out", status="replace", &
            action="write", iostat=ios)
      if (ios /= 0) then
         print "(a)", "golden_io: cannot write "//dir(1:n)//"/"//trim(case)//".out"
         stop 1
      end if
      active = .true.
      write (gu, "(a)") "# terco golden output v1 "//trim(case)
      flush (gu)
   end subroutine golden_begin

   !
   ! One quantity. `key` must not contain spaces.
   !
   subroutine golden_put(key, value)
      character(len=*), intent(in) :: key
      real(dp), intent(in) :: value
      if (.not. active) return
      write (gu, "(a,1x,es24.16e3)") trim(key), value
      flush (gu)
   end subroutine golden_put

   !
   ! Three scalars of a Fock-like matrix `g` built from density `d`: its
   ! contraction with the density (energy-like), its Frobenius norm, and a
   ! contraction with fixed, uneven weights, which moves when any element
   ! does. Recorded as TR_D<tag>, NORM_<tag> and WSUM_<tag>.
   !
   subroutine golden_fock(tag, d, g)
      character(len=*), intent(in) :: tag
      real(dp), intent(in) :: d(:, :), g(:, :)
      integer :: i, j
      real(dp) :: wsum
      if (.not. active) return
      wsum = 0.0_dp
      do j = 1, size(g, 2)
         do i = 1, size(g, 1)
            wsum = wsum + sin(real(7*i + 11*j, dp))*g(i, j)
         end do
      end do
      call golden_put("TR_D"//tag, sum(d*g))
      call golden_put("NORM_"//tag, sqrt(sum(g*g)))
      call golden_put("WSUM_"//tag, wsum)
   end subroutine golden_fock

   !
   ! Mark the output complete.
   !
   subroutine golden_end()
      if (.not. active) return
      write (gu, "(a)") "END"
      close (gu)
      active = .false.
   end subroutine golden_end

end module golden_io
