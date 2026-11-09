!
! F_0..F_m(t) into f(0:m) from ONE interpolated order and the downward
! recurrence, as an INCLUDE (same reason as trc_boys_eval.inc).
!
! boys_eval interpolates every order on its own: m + 1 Chebyshev series,
! six table loads each. In the MD Fock row kernels each thread rebuilds R
! per primitive quartet, so at L = 6 that was 42 table loads per quartet
! in every one of 18 threads, and ncu had the kernel stalled 43% of the
! time on a full load/store queue. Here only F_m is interpolated and
!
!     F_{k-1}(t) = (2t F_k(t) + exp(-t))/(2k - 1),
!
! which is stable downward for every t (each step adds positive terms), so
! the cost is six loads and one exponential whatever m is. Above BOYS_TMAX
! the asymptotic F_0 and the upward recurrence, as in boys_eval.
!
#ifndef TRC_BOYS_DOWN_INC
#define TRC_BOYS_DOWN_INC
   pure subroutine boys_down(m, t, f)
      !$acc routine seq
      integer,  intent(in)  :: m
      real(dp), intent(in)  :: t
      real(dp), intent(out) :: f(0:)
      integer  :: i, k, j, jj
      real(dp) :: x, x2, b0, b1, b2, et, tt
      ! The T from which the asymptotic F_0..F_m is exact to 1e-15 relative, by
      ! the highest order m (mpmath against gammainc, every order up to 16).
      real(dp), parameter :: BOYS_TASYM(0:16) = [ &
         32.5_dp, 36.5_dp, 40.0_dp, 43.0_dp, 46.0_dp, 48.5_dp, &
         51.0_dp, 53.0_dp, 55.5_dp, 57.5_dp, 60.0_dp, 62.0_dp, &
         64.0_dp, 66.0_dp, 68.0_dp, 70.0_dp, 72.0_dp]
      ! 1/(2k - 1): a division per step would cost more than the series.
      real(dp), parameter :: RODD(16) = [1.0_dp, 1.0_dp/3.0_dp, 1.0_dp/5.0_dp, 1.0_dp/7.0_dp, &
         1.0_dp/9.0_dp, 1.0_dp/11.0_dp, 1.0_dp/13.0_dp, 1.0_dp/15.0_dp, 1.0_dp/17.0_dp, &
         1.0_dp/19.0_dp, 1.0_dp/21.0_dp, 1.0_dp/23.0_dp, 1.0_dp/25.0_dp, 1.0_dp/27.0_dp, &
         1.0_dp/29.0_dp, 1.0_dp/31.0_dp]

      if (t >= BOYS_TMAX) then
         if (t >= BOYS_TASYM(m)) then
            ! The asymptotic F_m(T) = (2m-1)!!/(2T)^m sqrt(pi/(4T)): no erf,
            ! no exp(-T). Only from BOYS_TASYM(m) on, where it is exact to
            ! 1e-15 relative for every order up to m; at BOYS_TMAX itself it is
            ! not for high orders (1e-6 at m = 8, 1e-4 at m = 12 at T = 30),
            ! and check_nosym_f's f classes saw that. In a large molecule most
            ! primitive quartets are far apart, so this is the common path: the
            ! erf and exp below cost the MD J engine 0.9 s of a 5.4 s Fock
            ! build on gly50 cc-pVDZ (V100).
            tt = 1.0_dp/t
            f(0) = 0.88622692545275801365_dp*sqrt(tt)
            do k = 1, m
               f(k) = real(2*k - 1, dp)*f(k - 1)*(0.5_dp*tt)
            end do
            return
         end if
         tt = sqrt(t)
         f(0) = 0.88622692545275801365_dp*erf(tt)/tt
         et = exp(-t)
         do k = 1, m
            f(k) = (real(2*k - 1, dp)*f(k - 1) - et)*(0.5_dp/t)
         end do
         return
      end if

      i = int(t*BOYS_DTINV)
      if (i >= BOYS_NGRID) i = BOYS_NGRID - 1
      x = 2.0_dp*(t - real(i, dp)*BOYS_DT)*BOYS_DTINV - 1.0_dp
      x2 = 2.0_dp*x
      j = i*(BOYS_MMAX + 1)*(BOYS_NCHEB + 1) + m*(BOYS_NCHEB + 1)
      b1 = 0.0_dp; b2 = 0.0_dp
      do jj = BOYS_NCHEB, 1, -1
         b0 = x2*b1 - b2 + boys_table(j + jj + 1)
         b2 = b1; b1 = b0
      end do
      f(m) = x*b1 - b2 + boys_table(j + 1)
      if (m == 0) return
      et = exp(-t)
      tt = 2.0_dp*t
      do k = m, 1, -1
         f(k - 1) = (tt*f(k) + et)*RODD(k)
      end do
   end subroutine boys_down
#endif
