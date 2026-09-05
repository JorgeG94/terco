!
! F_0..F_m(t) into f(0:m), for any m <= BOYS_MMAX, as an INCLUDE.
!
! Pure and side-effect free so it runs inside a device kernel. It is a
! file to include after `contains` rather than a module procedure to
! import, because stdpar only inlines within a translation unit: called
! across a module boundary it is a real device call per primitive
! quartet, with every argument passed by address. Each kernel module
! includes its own copy; the guard keeps a unit from getting two. The
! including module supplies dp and the table constants from trc_boys.
!
#ifndef TRC_BOYS_EVAL_INC
#define TRC_BOYS_EVAL_INC
   pure subroutine boys_eval(m, t, f)
      !$acc routine seq
      integer,  intent(in)  :: m
      real(dp), intent(in)  :: t
      real(dp), intent(out) :: f(0:)
      integer  :: i, k, j, jj, base
      real(dp) :: x, x2, b0, b1, b2, et, tt

      if (t >= BOYS_TMAX) then
         ! Large T.  Needs erf regardless, and is rare enough not to matter.
         tt = sqrt(t)
         f(0) = 0.88622692545275801365_dp*erf(tt)/tt
         et = exp(-t)
         do k = 1, m
            f(k) = (real(2*k - 1, dp)*f(k - 1) - et)*(0.5_dp/t)
         end do
         return
      end if

      !
      ! Degree-5 Chebyshev, evaluated per order by Clenshaw recurrence.
      !
      ! Direct per order -- no downward recurrence, so no division and no
      ! exp(-T).  This is the structure used in [1] -- one interpolator per
      ! order, at a degree that keeps the whole table inside L2.
      !
      i = int(t*BOYS_DTINV)
      if (i >= BOYS_NGRID) i = BOYS_NGRID - 1
      ! map the interval onto [-1, 1]
      x = 2.0_dp*(t - real(i, dp)*BOYS_DT)*BOYS_DTINV - 1.0_dp
      x2 = 2.0_dp*x
      base = i*(BOYS_MMAX + 1)*(BOYS_NCHEB + 1)

      do k = 0, m
         j = base + k*(BOYS_NCHEB + 1)
         b1 = 0.0_dp; b2 = 0.0_dp
         do jj = BOYS_NCHEB, 1, -1
            b0 = x2*b1 - b2 + boys_table(j + jj + 1)
            b2 = b1; b1 = b0
         end do
         f(k) = x*b1 - b2 + boys_table(j + 1)
      end do
   end subroutine boys_eval
#endif
