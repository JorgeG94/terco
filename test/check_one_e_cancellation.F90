module check_one_e_cancellation_suite
   !! The one-electron kernels where one primitive is diffuse and the other tight.
   !!
   !! The nuclear attraction builds its recurrence table about A or B per
   !! primitive pair (gen_1e.py, emit_class): about A the transfer of b across
   !! the bond cancels when A's primitive is diffuse and B's tight, the loss
   !! growing like (exponent ratio x |AB|)^lb. V relative error of that
   !! transfer at ratio 25, |AB| = 2.5: 3.0e-14 at l = 2, 1.9e-12 at l = 3,
   !! 2.3e-11 at l = 4 (4.5e-15 at l = 1, within the tolerance).
   !! Single primitives on two centres and a unit charge on a third, both ways
   !! round (the tight-first block is compared with the transpose of the
   !! diffuse-first one: the same integral, and libcint is accurate only in
   !! the diffuse-first order), at l = 1 .. ONE_E_LMAX (d at LMAX=2, f at 3,
   !! g at 4, h at 5), against pyscf's Cartesian int1e_ovlp, int1e_kin and
   !! -int1e_rinv (test/pyscf_refs/generate.py, gen_one_e_cancellation: libcint
   !! checked against mpmath there, 2.6e-15).
   !!
   !! The coefficient is the unit-normalised x^l, so the block is terco's own AO
   !! normalisation. Tolerances, relative to the block's largest element: 1e-14
   !! for V (2.1e-15 at most against pyscf). A kernel that always transfers
   !! about A fails it from l = 2 up, so the test discriminates at
   !! ONE_E_LMAX >= 2; a build with ONE_E_LMAX = 1 does not catch it. 1e-13
   !! for S and T (1.6e-15 and 1.1e-14).
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_1e_kernels, only: one_e_dispatch, ONE_E_LMAX
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_vector
   implicit none
   private
   public :: collect_one_e_cancellation

   real(dp), parameter :: TOL_V = 1.0e-14_dp, TOL_ST = 1.0e-13_dp
   !! relative to the largest element of the block

   integer, parameter :: LTOP = 5
   !! the highest l the reference holds
   integer, parameter :: L_BEFORE = 4, L_CTR = 3
   !! the reference holds the single-primitive blocks of l = 1 .. L_BEFORE, then the
   !! contracted blocks of l = 1 .. L_CTR, then the single-primitive block of l = LTOP
   real(dp), parameter :: E_DIFFUSE = 1.6_dp, E_TIGHT = 40.0_dp
   real(dp), parameter :: R_B(3) = 2.5_dp*[0.48_dp, 0.36_dp, 0.8_dp]
   real(dp), parameter :: R_B3(3) = 3.0_dp*[0.48_dp, 0.36_dp, 0.8_dp]
   real(dp), parameter :: R_C(3) = [1.1_dp, -0.7_dp, 0.4_dp]

contains

   subroutine collect_one_e_cancellation(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("S, T, V with diffuse and tight primitives, l up to this build's ceiling", test_cancellation), &
                   new_unittest("S, T, V of a contraction mixing diffuse and tight primitives", test_contracted)]
   end subroutine collect_one_e_cancellation

   pure real(dp) function gnorm(l, a)
      !! The coefficient that makes x^l exp(-a r^2) unit-normalised: <x^l|x^l>
      !! = (2l-1)!!/(4a)^l (pi/2a)^(3/2), the convention terco's shells have.
      integer, intent(in) :: l
      real(dp), intent(in) :: a
      real(dp) :: df
      integer :: i
      df = 1.0_dp
      do i = 2*l - 1, 1, -2
         df = df*real(i, dp)
      end do
      gnorm = 1.0_dp/sqrt(df/(4.0_dp*a)**l*(3.14159265358979323846_dp/(2.0_dp*a))**1.5_dp)
   end function gnorm

   subroutine test_cancellation(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: u, l, n, nc, cfg, k
      real(dp), allocatable :: ref(:), refs(:, :), s(:), t(:), v(:)
      real(dp) :: ra(3), rb(3), ea, eb, ca(1), cb(1), rat(3, 1), z(1), err(3), mx
      character(len=24) :: tag
      call ref_require("one_e_cancellation.ref", u)
      z = 1.0_dp
      rat(:, 1) = R_C
      do l = 1, LTOP
         if (l > L_BEFORE) call skip_contracted(u)
         n = ((l + 1)*(l + 2)/2)**2
         nc = (l + 1)*(l + 2)/2
         allocate (ref(n), refs(n, 3), s(n), t(n), v(n))
         do cfg = 1, 2
            ! cfg 1: A diffuse, B tight (the cancelling order); cfg 2: the mirror.
            ra = 0.0_dp
            rb = R_B
            ea = E_DIFFUSE
            eb = E_TIGHT
            if (cfg == 2) then
               ea = E_TIGHT
               eb = E_DIFFUSE
               ra = R_B
               rb = 0.0_dp
            end if
            ca = gnorm(l, ea)
            cb = gnorm(l, eb)
            if (l <= ONE_E_LMAX) call one_e_dispatch(l, l, 1, 1, [ea], ca, [eb], cb, ra, rb, 1, z, rat, s, t, v)
            do k = 1, 3
               if (cfg == 1) then
                  tag = "STV"(k:k)//"_l"//to_string(l)//"_c1"
                  call ref_read_vector(u, trim(tag), n, ref)
                  refs(:, k) = ref
               end if
               if (l > ONE_E_LMAX) cycle
               ref = refs(:, k)
               if (cfg == 2) ref = reshape(transpose(reshape(refs(:, k), [nc, nc])), [n])
               select case (k)
               case (1); err(k) = maxval(abs(s - ref))/maxval(abs(ref))
               case (2); err(k) = maxval(abs(t - ref))/maxval(abs(ref))
               case default; err(k) = maxval(abs(v - ref))/maxval(abs(ref))
               end select
            end do
            if (l <= ONE_E_LMAX) then
               print '(4x,a,i0,a,i0,a,3es10.2)', "l = ", l, ", config ", cfg, ": relative S, T, V error ", err
               call check(error, err(1) <= TOL_ST .and. err(2) <= TOL_ST .and. err(3) <= TOL_V, "l = "//to_string(l)//", config "//to_string(cfg)// &
                          ": error S "//to_string(err(1), "es9.2")//", T "//to_string(err(2), "es9.2")// &
                          ", V "//to_string(err(3), "es9.2")//" exceeds the tolerance (S, T "//to_string(TOL_ST, "es9.2")// &
                          ", V "//to_string(TOL_V, "es9.2")//")")
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
            end if
         end do
         deallocate (ref, refs, s, t, v)
      end do
      call ref_close(u)
   end subroutine test_cancellation

   subroutine skip_contracted(u)
      !! Step the sequential reader over the contracted blocks.
      integer, intent(in) :: u
      real(dp), allocatable :: ref(:)
      character(len=24) :: tag
      integer :: l, n, k
      do l = 1, L_CTR
         n = ((l + 1)*(l + 2)/2)**2
         allocate (ref(n))
         do k = 1, 3
            tag = "STV"(k:k)//"_ctr_l"//to_string(l)
            call ref_read_vector(u, trim(tag), n, ref)
         end do
         deallocate (ref)
      end do
   end subroutine skip_contracted

   subroutine test_contracted(error)
      !! One contraction on each centre, exponents (1.0, 60) on A and (100, 40)
      !! on B, |AB| = 3.0. Three of the four primitive pairs have the tight
      !! primitive on B and build the nuclear-attraction table about B: these
      !! are the cancelling pairs and they dominate the error of a kernel that
      !! always builds about A. The fourth (60 on A, 40 on B) builds about A but
      !! is tight against tight at 3 bohr and adds nothing measurable; the
      !! A-built branch is covered by configuration 2 of the single-primitive
      !! test. The reference is the weighted sum of the single-primitive blocks
      !! (S, T from pyscf, V from mpmath).
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: EXA(2) = [1.0_dp, 60.0_dp], EXB(2) = [100.0_dp, 40.0_dp]
      real(dp), parameter :: WA(2) = [0.6_dp, 0.5_dp], WB(2) = [0.45_dp, 0.55_dp]
      integer :: u, l, n, k, i
      real(dp), allocatable :: ref(:), s(:), t(:), v(:)
      real(dp) :: ca(2), cb(2), rat(3, 1), z(1), err(3)
      character(len=24) :: tag
      call ref_require("one_e_cancellation.ref", u)
      z = 1.0_dp
      rat(:, 1) = R_C
      ! The reference is read in file order: step over the single-primitive blocks
      ! that precede the contracted ones.
      do l = 1, L_BEFORE
         n = ((l + 1)*(l + 2)/2)**2
         allocate (ref(n))
         do k = 1, 3
            tag = "STV"(k:k)//"_l"//to_string(l)//"_c1"
            call ref_read_vector(u, trim(tag), n, ref)
         end do
         deallocate (ref)
      end do
      do l = 1, min(L_CTR, ONE_E_LMAX)
         n = ((l + 1)*(l + 2)/2)**2
         allocate (ref(n), s(n), t(n), v(n))
         do i = 1, 2
            ca(i) = WA(i)*gnorm(l, EXA(i))
            cb(i) = WB(i)*gnorm(l, EXB(i))
         end do
         call one_e_dispatch(l, l, 2, 2, EXA, ca, EXB, cb, [0.0_dp, 0.0_dp, 0.0_dp], R_B3, 1, z, rat, s, t, v)
         do k = 1, 3
            tag = "STV"(k:k)//"_ctr_l"//to_string(l)
            call ref_read_vector(u, trim(tag), n, ref)
            select case (k)
            case (1); err(k) = maxval(abs(s - ref))/maxval(abs(ref))
            case (2); err(k) = maxval(abs(t - ref))/maxval(abs(ref))
            case default; err(k) = maxval(abs(v - ref))/maxval(abs(ref))
            end select
         end do
         print '(4x,a,i0,a,3es10.2)', "contracted l = ", l, ": relative S, T, V error ", err
         call check(error, err(1) <= TOL_ST .and. err(2) <= TOL_ST .and. err(3) <= TOL_V, "contracted l = "//to_string(l)// &
                    ": error S "//to_string(err(1), "es9.2")//", T "//to_string(err(2), "es9.2")// &
                    ", V "//to_string(err(3), "es9.2")//" exceeds the tolerance")
         if (allocated(error)) then
            call ref_close(u)
            return
         end if
         deallocate (ref, s, t, v)
      end do
      call ref_close(u)
   end subroutine test_contracted

end module check_one_e_cancellation_suite

program check_one_e_cancellation
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_one_e_cancellation_suite, only: collect_one_e_cancellation
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_one_e_cancellation, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_one_e_cancellation
