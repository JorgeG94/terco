module check_rccsd_diis_suite
   !!
   !! The packed DIIS kernels of the RI-CCSD solver (trc_rccsd_diis, `packed=.true.`) against plain loops over the
   !! (i,j,a,b) index, as a test-drive suite. The kernels pass values computed inside their `do concurrent` (the pair
   !! partner `r`, the packed row) to item routines through the TRC_DC_INDEX macro (src/inc/trc_dc_index.inc): ifx
   !! 2025/2026 at -O2 treats such a call as loop-invariant when it is handed the plain variable and every element gets
   !! the first one's result, which no energy check of a converged run need catch (the iteration still converges, to
   !! what the wrong DIIS lets it). Not reproduced at this size: the check passed with the unprotected call on ifx at -O2
   !! and -O3, so it pins the packed arithmetic, not the miscompilation. Checked on pair-symmetric amplitudes
   !! of no = 3, nv = 4:
   !!   - the stored error: (new1 - old1, then the symmetrised t2 step times sqrt(2) off the diagonal) at its packed place;
   !!   - the stored vector: (new1, (t2(s) + t2(r))/2) at its packed place;
   !!   - the step maxima and ||e||_2 over the full vector;
   !!   - after two pushes, `next`'s t2 is sum_m coef(m) vec(k,m) at BOTH (i,j,a,b) and (j,i,b,a) for every element.
   !! Tolerance 1e-14 (a sum of two pushed vectors with weights of O(1)); the first two are exact copies of the loops'
   !! arithmetic and are compared to 1e-15.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_linalg, only: trc_linalg_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_diis, only: trc_rccsd_diis_t
   implicit none
   private
   public :: collect_rccsd_diis

   integer, parameter :: NO = 3, NV = 4
   real(dp), parameter :: SQRT2 = 1.4142135623730951_dp

contains

   subroutine collect_rccsd_diis(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("packed push and next against plain loops", test_packed)]
   end subroutine collect_rccsd_diis

   subroutine make_amps(amps, seed)
      !! Pair-symmetric amplitudes, t2(i,j,a,b) = t2(j,i,b,a), from a smooth function of the indices
      type(trc_rccsd_amps_t), intent(inout) :: amps
      real(dp), intent(in) :: seed
      integer :: i, j, a, b

      call amps%init(NO, NV)
      do a = 1, NV
         do i = 1, NO
            amps%t1(i, a) = sin(seed + 0.7_dp*i + 1.3_dp*a)
         end do
      end do
      do b = 1, NV
         do a = 1, NV
            do j = 1, NO
               do i = 1, NO
                  amps%t2(i, j, a, b) = cos(seed + 0.31_dp*(i*a + j*b) + 0.11_dp*(i + j)*(a + b) + 0.2_dp*i*j*a*b)
               end do
            end do
         end do
      end do
      amps%t2 = 0.5_dp*(amps%t2 + reshape([(((( &
                                              amps%t2(j, i, b, a), i=1, NO), j=1, NO), a=1, NV), b=1, NV)], shape(amps%t2)))
      !$acc update device(amps%t1, amps%t2)
   end subroutine make_amps

   subroutine test_packed(error)
      !! Owns every resource, so each early return of `run_packed` still releases them
      type(error_type), allocatable, intent(out) :: error
      type(trc_linalg_t) :: la
      type(trc_rccsd_diis_t) :: diis
      type(trc_rccsd_amps_t) :: old, new, new2, out

      call make_amps(old, 0.1_dp)
      call make_amps(new, 0.9_dp)
      call make_amps(new2, 1.7_dp)
      call out%init(NO, NV)
      ! the handles of the device BLAS: the dot products of `push` need them (a host build ignores it)
      call la%init(2)
      call diis%init(NO, NV, 4, packed=.true.)
      call run_packed(error, la, diis, old, new, new2, out)
      call diis%release()
      call la%release()
      call old%release()
      call new%release()
      call new2%release()
      call out%release()
   end subroutine test_packed

   subroutine run_packed(error, la, diis, old, new, new2, out)
      type(error_type), allocatable, intent(out) :: error
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_diis_t), intent(inout) :: diis
      type(trc_rccsd_amps_t), intent(inout) :: old, new, new2, out
      real(dp), allocatable :: e(:), v(:)
      real(dp) :: dmax, d_e, d_v, d_max, d_norm, d_next, x, step, step_full
      integer :: i, j, a, b, p, q, k, n1, nov, ip, iq
      logical :: extrapolated

      call check(error, .not. diis%status%has_error(), "init: "//diis%status%get_message())
      if (allocated(error)) return
      call diis%push(la, old, new, .true.)
      call check(error, .not. diis%status%has_error() .and. diis%finite, "push: "//diis%status%get_message())
      if (allocated(error)) return
      n1 = NO*NV
      nov = n1
      allocate (e(diis%n), v(diis%n))
      !$acc update self(diis%err, diis%vec)
      e = diis%err(:, 1)
      v = diis%vec(:, 1)
      d_e = 0.0_dp
      d_v = 0.0_dp
      step_full = 0.0_dp
      step = 0.0_dp
      do a = 1, NV
         do i = 1, NO
            d_e = max(d_e, abs(e(i + NO*(a - 1)) - (new%t1(i, a) - old%t1(i, a))))
            d_v = max(d_v, abs(v(i + NO*(a - 1)) - new%t1(i, a)))
            step_full = step_full + (new%t1(i, a) - old%t1(i, a))**2
            step = max(step, abs(new%t1(i, a) - old%t1(i, a)))
         end do
      end do
      do b = 1, NV
         do a = 1, NV
            do j = 1, NO
               do i = 1, NO
                  x = new%t2(i, j, a, b) - old%t2(i, j, a, b)
                  step_full = step_full + x*x
                  step = max(step, abs(x))
                  ip = i + NO*(a - 1)
                  iq = j + NO*(b - 1)
                  if (ip >= iq) then
                     k = nov + ip + (iq - 1)*(2*nov - iq)/2
                     x = 0.5_dp*(x + (new%t2(j, i, b, a) - old%t2(j, i, b, a)))
                     if (ip /= iq) x = SQRT2*x
                     d_e = max(d_e, abs(e(k) - x))
                     x = 0.5_dp*(new%t2(i, j, a, b) + new%t2(j, i, b, a))
                     d_v = max(d_v, abs(v(k) - x))
                  end if
               end do
            end do
         end do
      end do
      call check(error, d_e < 1.0e-15_dp, "the packed error differs from the plain loop by "//to_string(d_e, "es9.2"))
      if (allocated(error)) return
      call check(error, d_v < 1.0e-15_dp, "the packed vector differs from the plain loop by "//to_string(d_v, "es9.2"))
      if (allocated(error)) return
      d_max = abs(diis%step_max - step)
      d_norm = abs(diis%step_norm2 - sqrt(step_full))
      call check(error, d_max < 1.0e-15_dp .and. d_norm < 1.0e-13_dp, "step_max off by "//to_string(d_max, "es9.2")// &
                 ", step_norm2 by "//to_string(d_norm, "es9.2"))
      if (allocated(error)) return

      ! a second entry, then the extrapolation
      call diis%push(la, new, new2, .true.)
      call check(error, .not. diis%status%has_error(), "second push: "//diis%status%get_message())
      if (allocated(error)) return
      call diis%next(la, new2, out, extrapolated)
      call check(error, extrapolated, "the extrapolation did not run: "//diis%status%get_message())
      if (allocated(error)) return
      !$acc update self(out%t1, out%t2, diis%vec, diis%coef)
      d_next = 0.0_dp
      do b = 1, NV
         do a = 1, NV
            do j = 1, NO
               do i = 1, NO
                  ip = i + NO*(a - 1)
                  iq = j + NO*(b - 1)
                  p = max(ip, iq)
                  q = min(ip, iq)
                  k = nov + p + (q - 1)*(2*nov - q)/2
                  x = 0.0_dp
                  do p = 1, diis%nstored
                     x = x + diis%coef(p)*diis%vec(k, p)
                  end do
                  d_next = max(d_next, abs(out%t2(i, j, a, b) - x))
               end do
            end do
         end do
      end do
      call check(error, d_next < 1.0e-14_dp, "the extrapolated t2 differs from the plain sum by "//to_string(d_next, "es9.2"))
      dmax = 0.0_dp
      do b = 1, NV
         do a = 1, NV
            do j = 1, NO
               do i = 1, NO
                  dmax = max(dmax, abs(out%t2(i, j, a, b) - out%t2(j, i, b, a)))
               end do
            end do
         end do
      end do
      if (.not. allocated(error)) call check(error, dmax == 0.0_dp, "the extrapolated t2 is not pair symmetric: "// &
                                             to_string(dmax, "es9.2"))
   end subroutine run_packed

end module check_rccsd_diis_suite

program check_rccsd_diis
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_diis_suite, only: collect_rccsd_diis
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_rccsd_diis, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_diis
