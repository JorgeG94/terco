!
! PARITY-SKIP SPIKE -- NOT FOR MERGING.  docs/SYMMETRY_PARITY_SKIP.md.
!
! Scheme (a) at kernel level: the VRR + HRR kernels that
! scripts/spike_parity_gen.py emits for a class with dead coordinates Z (all
! four centres have those coordinates equal to 0) against the kernel with no
! mask (the text gen_perclass.py ships), on random geometries that satisfy the
! mask. Compares the (ab|cd) block (relative 1e-13) and times both.
! usage: spike_parity_kernel [repeats]
!
module spike_pk_driver
   use pic_types, only: dp, int64
   use trc_boys, only: boys_eval
   implicit none
   real(dp), parameter :: PI = 3.14159265358979323846_dp
   integer, parameter :: NP = 6   !! primitive quartets per call
contains

   subroutine make_params(lt, zmask, seed, prm, geo)
      !! NP random primitive quartets over one set of four centres with the
      !! coordinates in zmask equal to 0, in the packed layout of the kernels.
      integer, intent(in) :: lt, zmask, seed
      real(dp), intent(out) :: prm(lt + 19, NP), geo(6)
      real(dp) :: ctr(3, 4), ex(4), pcen(3), qcen(3), wcen(3), zeta, eta, zpe, rho, t, ab2, cd2, f(0:lt)
      integer :: ip, d, k
      integer(int64) :: n
      real(dp) :: u
      n = int(seed, int64)
      do k = 1, 4
         do d = 1, 3
            n = mod(n*48271_int64, 2147483647_int64)
            u = real(n, dp)/2147483647.0_dp
            ctr(d, k) = merge(0.0_dp, 3.0_dp*u - 1.5_dp, btest(zmask, d - 1))
         end do
      end do
      geo = [ctr(:, 1) - ctr(:, 2), ctr(:, 3) - ctr(:, 4)]
      do ip = 1, NP
         do k = 1, 4
            n = mod(n*48271_int64, 2147483647_int64)
            ex(k) = 0.3_dp + 2.5_dp*real(n, dp)/2147483647.0_dp
         end do
         zeta = ex(1) + ex(2); eta = ex(3) + ex(4); zpe = zeta + eta
         pcen = (ex(1)*ctr(:, 1) + ex(2)*ctr(:, 2))/zeta
         qcen = (ex(3)*ctr(:, 3) + ex(4)*ctr(:, 4))/eta
         wcen = (zeta*pcen + eta*qcen)/zpe
         rho = zeta*eta/zpe
         t = rho*sum((pcen - qcen)**2)
         ab2 = sum((ctr(:, 1) - ctr(:, 2))**2); cd2 = sum((ctr(:, 3) - ctr(:, 4))**2)
         call boys_eval(lt, t, f)
         prm(1, ip) = 2.0_dp*PI**2.5_dp/(zeta*eta*sqrt(zpe))*exp(-ex(1)*ex(2)/zeta*ab2 - ex(3)*ex(4)/eta*cd2)
         prm(2:lt + 2, ip) = f(0:lt)
         k = lt + 3
         prm(k:k + 2, ip) = pcen - ctr(:, 1)                 ! PA
         prm(k + 3:k + 5, ip) = qcen - ctr(:, 3)             ! QC
         prm(k + 6:k + 8, ip) = wcen - pcen                  ! WP
         prm(k + 9:k + 11, ip) = wcen - qcen                 ! WQ
         prm(k + 12, ip) = 0.5_dp/zeta; prm(k + 13, ip) = 0.5_dp/eta; prm(k + 14, ip) = 0.5_dp/zpe
         prm(k + 15, ip) = rho/zeta; prm(k + 16, ip) = rho/eta
      end do
   end subroutine make_params
end module spike_pk_driver

program spike_parity_kernel
   use pic_types, only: dp
   use spike_pk_driver
   use spike_parity_kernels
   implicit none
   type(case_t) :: cases(NCASE)
   integer :: reps, nfail, i, j, r, ref, seed
   character(len=16) :: arg
   real(dp), allocatable :: prm(:, :), vref(:), vtest(:)
   real(dp) :: geo(6), t0, t1, tref, ttest, err, scl, sink

   reps = 20000
   if (command_argument_count() > 0) then
      call get_command_argument(1, arg); read (arg, *) reps
   end if
   call init_cases(cases)
   nfail = 0
   sink = 0.0_dp
   print '(a)', 'class  mask kind    max rel diff   ns/call  vs shipped'
   do i = 1, NCASE
      ref = 0
      do j = 1, NCASE
         if (cases(j)%name == cases(i)%name .and. cases(j)%zmask == 0 .and. cases(j)%kind == 'all') ref = j
      end do
      if (ref == 0) stop 'no unmasked reference kernel for a class'
      allocate (prm(cases(i)%lt + 19, NP), vref(cases(i)%nout), vtest(cases(i)%nout))
      do seed = 1, 5
         call make_params(cases(i)%lt, cases(i)%zmask, 7919*seed, prm, geo)
         call cases(ref)%run(NP, cases(i)%lt + 19, prm, geo, vref)
         call cases(i)%run(NP, cases(i)%lt + 19, prm, geo, vtest)
         scl = max(maxval(abs(vref)), tiny(1.0_dp))
         err = maxval(abs(vref - vtest))/scl
         if (err > 1.0e-13_dp) nfail = nfail + 1
      end do
      call cpu_time(t0)
      do r = 1, reps
         prm(1, 1) = prm(1, 1)*(1.0_dp + 1.0e-16_dp*r)
         call cases(ref)%run(NP, cases(i)%lt + 19, prm, geo, vref)
         sink = sink + vref(1)
      end do
      call cpu_time(t1); tref = (t1 - t0)/reps*1.0e9_dp
      call cpu_time(t0)
      do r = 1, reps
         prm(1, 1) = prm(1, 1)*(1.0_dp + 1.0e-16_dp*r)
         call cases(i)%run(NP, cases(i)%lt + 19, prm, geo, vtest)
         sink = sink + vtest(1)
      end do
      call cpu_time(t1); ttest = (t1 - t0)/reps*1.0e9_dp
      print '(a4,2x,i3,2x,a7,es11.2,f10.0,f9.2,a)', cases(i)%name, cases(i)%zmask, trim(cases(i)%kind), err, ttest, tref/ttest, 'x'
      deallocate (prm, vref, vtest)
   end do
   if (sink == 12345.678_dp) print *, sink
   if (nfail > 0) then
      print '(a,i0)', 'FAILED comparisons: ', nfail
      error stop 1
   end if
   print '(a)', 'all masked kernels agree with the unmasked ones'
end program spike_parity_kernel
