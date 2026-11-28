module check_rccsd_lambda_update_suite
   !!
   !! One update of the RI-CCSD Lambda equations (src/post_hf/ccsd/trc_rccsd_lambda_update.F90) on water cc-pVDZ
   !! with cc-pVDZ-RIFIT, frozen core 1 and 0, at RANDOM multipliers (and, except for the adjoint, random
   !! amplitudes): t1 and l1 without symmetry, t2 and l2 with only x(j,i,b,a) = x(i,j,a,b), so that a swapped
   !! index or a factor shows (converged amplitudes and multipliers hide some). No reference file.
   !!
   !!   against the dense formula  pyscf's update_lambda with (ia|bc), (ab|cd) and wvvov in host loops
   !!                              (rccsd_lambda_host), elementwise
   !!   against the adjoint        terco's OWN CCSD residual r(t) = D (t' - t) of trc_rccsd_update: with
   !!                                  Lag(t) = E(t) + 2 l1 . r1(t) + (2 l2 - l2 swapped) . r2(t)
   !!                              the update's residual D (l' - l) must give dLag/dt along random directions,
   !!                              d Lag = 2 rho1 . d1 + (2 rho2 - rho2 swapped) . d2. Lag is a polynomial of
   !!                              degree 4 in t, so the five-point central difference is exact to round-off.
   !!                              This does not rest on pyscf's equations, or on the code this was ported from.
   !!                              At CONVERGED amplitudes only: the l2 terms of these (and of pyscf's)
   !!                              equations use the CCSD equations r(t) = 0, so away from convergence they are not
   !!                              the derivative (measured: relative differences of order one with l1 = 0, round-off
   !!                              with l2 = 0). The terms that vanish with r(t) are what the dense comparison covers.
   !!   forms and blocks           assembled, direct and half ladders, one virtual d per half-ladder pass, the
   !!                              intermediates' woovo term in blocks of 1 and 3 virtuals: the same update
   !!   pair symmetry              l2' (j,i,b,a) = l2' (i,j,a,b) bit for bit
   !!   refusals                   each by its message
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_scratch_t, rccsd_inter_energy
   use trc_rccsd_update, only: trc_rccsd_update_t
   use trc_rccsd_driver, only: trc_rccsd_options_t, trc_rccsd_state_t, trc_rccsd_solve
   use trc_rccsd_lambda_update, only: trc_rccsd_lambda_update_t, trc_rccsd_lambda_update_options_t
   use trc_rccsd_ladder, only: RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_DIRECT, RCCSD_LADDER_HALF
   use trc_device_memory, only: device_memory_set_cap
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: load_molecule, run_scf
   use rccsd_lambda_host, only: host_lambda_update
   implicit none
   private
   public :: collect_check_rccsd_lambda_update, release_fixture

   real(dp), parameter :: TOL_DENSE = 2.0e-14_dp    !! max abs difference over max(1, |reference|max); worst 1.1e-15
   real(dp), parameter :: TOL_ADJOINT = 2.0e-13_dp  !! |dLag - adjoint| over max(1, |dLag|); worst 1.1e-14
   real(dp), parameter :: TOL_FORMS = 2.0e-14_dp    !! another form or block size against the default; worst 8e-17
   integer, parameter :: NFIX = 2
   integer, parameter :: FIX_FROZEN(NFIX) = [1, 0]

   type :: fixture_t
      !! One frozen core: integrals, random amplitudes and multipliers, and terco's update of them
      logical :: built = .false.
      logical :: supported = .true.
      integer :: nocc = 0
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_amps_t) :: amps, l
      real(dp), allocatable :: new1(:, :), new2(:, :, :, :)
   end type fixture_t

   logical :: bound = .false.
   type(trc_linalg_t) :: la
   type(fixture_t) :: fix(NFIX)

contains

   subroutine collect_check_rccsd_lambda_update(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("frozen 1: against the dense formula", test_dense_frozen1), &
                   new_unittest("frozen 0: against the dense formula", test_dense_frozen0), &
                   new_unittest("frozen 1: adjoint of the CCSD residual", test_adjoint_frozen1), &
                   new_unittest("frozen 0: adjoint of the CCSD residual", test_adjoint_frozen0), &
                   new_unittest("ladder forms and block sizes", test_forms), &
                   new_unittest("refusals", test_refusals)]
   end subroutine collect_check_rccsd_lambda_update

   subroutine test_dense_frozen1(error)
      type(error_type), allocatable, intent(out) :: error
      call dense_case(error, 1)
   end subroutine test_dense_frozen1

   subroutine test_dense_frozen0(error)
      type(error_type), allocatable, intent(out) :: error
      call dense_case(error, 2)
   end subroutine test_dense_frozen0

   subroutine test_adjoint_frozen1(error)
      type(error_type), allocatable, intent(out) :: error
      call adjoint_case(error, 1)
   end subroutine test_adjoint_frozen1

   subroutine test_adjoint_frozen0(error)
      type(error_type), allocatable, intent(out) :: error
      call adjoint_case(error, 2)
   end subroutine test_adjoint_frozen0

   subroutine dense_case(error, which)
      !! terco's update against pyscf's dense one in host loops, element by element
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: which

      real(dp), allocatable :: ref1(:, :), ref2(:, :, :, :)
      real(dp) :: d1, d2
      integer :: no, nv

      call ensure_fixture(which, error)
      if (allocated(error) .or. .not. fix(which)%supported) return
      associate (f => fix(which))
         no = f%ints%no
         nv = f%ints%nv
         allocate (ref1(no, nv), ref2(no, no, nv, nv))
         call host_lambda_update(f%ints, f%amps%t1, f%amps%t2, f%l%t1, f%l%t2, ref1, ref2)
         d1 = maxval(abs(f%new1 - ref1))/max(1.0_dp, maxval(abs(ref1)))
         d2 = maxval(abs(f%new2 - ref2))/max(1.0_dp, maxval(abs(ref2)))
         print "(a,i0,a,2es10.2,a,2es10.2)", "  dense formula, case ", which, ": rel. diff l1, l2 ", d1, d2, &
            "   max |l1'|, |l2'| ", maxval(abs(ref1)), maxval(abs(ref2))
         ! written so that a NaN fails; the references must be of order one for the comparison to mean anything
         call check(error, maxval(abs(ref1)) > 1.0e-2_dp .and. maxval(abs(ref2)) > 1.0e-2_dp, "the reference vanishes")
         if (allocated(error)) return
         call check(error, d1 < TOL_DENSE, "l1 against the dense formula: "//to_string(d1, "es9.2"))
         if (allocated(error)) return
         call check(error, d2 < TOL_DENSE, "l2 against the dense formula: "//to_string(d2, "es9.2"))
         if (allocated(error)) return
         call check(error, pair_asymmetry(f%new2) == 0.0_dp, "l2' is not exactly pair-symmetric")
      end associate
   end subroutine dense_case

   function pair_asymmetry(x) result(gap)
      real(dp), intent(in) :: x(:, :, :, :)
      real(dp) :: gap
      integer :: i, j, a, b

      gap = 0.0_dp
      do b = 1, size(x, 4)
         do a = 1, size(x, 3)
            do j = 1, size(x, 2)
               do i = 1, size(x, 1)
                  gap = max(gap, abs(x(i, j, a, b) - x(j, i, b, a)))
               end do
            end do
         end do
      end do
   end function pair_asymmetry

   subroutine adjoint_case(error, which)
      !! dLag/dt along random directions against the update's residual (module header). Directions: t1 alone,
      !! t2 alone, both
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: which

      real(dp), parameter :: H = 0.25_dp, STENCIL(4) = [-2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp]
      real(dp), parameter :: WEIGHT(4) = [1.0_dp, -8.0_dp, 8.0_dp, -1.0_dp]
      type(trc_rccsd_update_t) :: upd
      type(trc_rccsd_scratch_t) :: scr
      type(trc_rccsd_amps_t) :: t, tnew, conv
      type(trc_rccsd_options_t) :: copts
      type(trc_rccsd_state_t) :: cc
      type(trc_rccsd_lambda_update_options_t) :: lopts
      real(dp), allocatable :: d1(:, :), d2(:, :, :, :), rho1(:, :), rho2(:, :, :, :), x(:, :, :, :)
      real(dp), allocatable :: new1(:, :), new2(:, :, :, :)
      real(dp) :: lag(4), fd, adjoint, diff
      integer :: no, nv, kind, k, i, j, a, b, form_used

      call ensure_fixture(which, error)
      if (allocated(error) .or. .not. fix(which)%supported) return
      associate (f => fix(which))
         no = f%ints%no
         nv = f%ints%nv
         ! the converged amplitudes (the fixture's integrals are those of this run), and the update of the
         ! fixture's random multipliers at them
         copts%nocc = f%nocc
         copts%nfrozen = FIX_FROZEN(which)
         copts%conv_energy = 1.0e-13_dp
         copts%conv_amplitude = 1.0e-11_dp
         call trc_rccsd_solve(f%bas, f%aux, f%pl, f%res%cmo(:, :, 1), f%res%eps(:, 1), copts, la, cc)
         call check(error, cc%result%converged, "CCSD: "//trim(cc%result%message))
         if (allocated(error)) return
         call conv%init(no, nv)
         conv%t1 = cc%amps%t1
         conv%t2 = cc%amps%t2
         call conv%push()
         call cc%release()
         call run_update(error, f, conv, lopts, new1, new2, form_used)
         if (allocated(error)) return

         call upd%init(f%ints)
         call check(error, .not. upd%status%has_error(), "CCSD update init: "//upd%status%get_message())
         if (allocated(error)) return
         call scr%init(f%ints)
         call t%init(no, nv)
         call tnew%init(no, nv)

         allocate (rho1(no, nv), rho2(no, no, nv, nv), d1(no, nv), d2(no, no, nv, nv), x(no, no, nv, nv))
         do a = 1, nv
            do i = 1, no
               rho1(i, a) = (f%ints%eps_o(i) - f%ints%eps_v(a))*(new1(i, a) - f%l%t1(i, a))
            end do
         end do
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     rho2(i, j, a, b) = (f%ints%eps_o(i) - f%ints%eps_v(a) + f%ints%eps_o(j) - f%ints%eps_v(b)) &
                                        *(new2(i, j, a, b) - f%l%t2(i, j, a, b))
                  end do
               end do
            end do
         end do

         do kind = 1, 3
            call random_number(d1)
            d1 = d1 - 0.5_dp
            call random_symmetric(x, d2)
            if (kind == 1) d2 = 0.0_dp
            if (kind == 2) d1 = 0.0_dp
            do k = 1, 4
               t%t1 = conv%t1 + STENCIL(k)*H*d1
               t%t2 = conv%t2 + STENCIL(k)*H*d2
               call t%push()
               call lagrangian(f, upd, scr, t, tnew, lag(k))
               call check(error, .not. upd%status%has_error(), "CCSD update: "//upd%status%get_message())
               if (allocated(error)) exit
            end do
            if (allocated(error)) exit
            fd = sum(WEIGHT*lag)/(12.0_dp*H)
            adjoint = 2.0_dp*sum(rho1*d1)
            do b = 1, nv
               do a = 1, nv
                  do j = 1, no
                     do i = 1, no
                        adjoint = adjoint + (2.0_dp*rho2(i, j, a, b) - rho2(i, j, b, a))*d2(i, j, a, b)
                     end do
                  end do
               end do
            end do
            diff = abs(fd - adjoint)/max(1.0_dp, abs(fd))
            print "(a,i0,a,i0,a,es22.14,a,es10.2)", "  adjoint, case ", which, " direction ", kind, ": dLag ", fd, &
               "  rel. diff ", diff
            call check(error, abs(fd) > 1.0e-2_dp, "the derivative vanishes: the direction tests nothing")
            if (allocated(error)) exit
            call check(error, diff < TOL_ADJOINT, "direction "//to_string(kind)//": finite difference "// &
                       to_string(fd, "es22.14")//" against the update's "//to_string(adjoint, "es22.14"))
            if (allocated(error)) exit
         end do
         call upd%release()
         call scr%release()
         call conv%release()
         call t%release()
         call tnew%release()
      end associate
   end subroutine adjoint_case

   subroutine lagrangian(f, upd, scr, t, tnew, lag)
      !! E(t) + 2 l1 . r1 + (2 l2 - l2 swapped) . r2, with r = D (update(t) - t) from terco's CCSD update
      type(fixture_t), intent(in) :: f
      type(trc_rccsd_update_t), intent(inout) :: upd
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_rccsd_amps_t), intent(in) :: t
      type(trc_rccsd_amps_t), intent(inout) :: tnew
      real(dp), intent(out) :: lag

      real(dp) :: es, ed, d
      integer :: i, j, a, b

      call upd%step(f%ints, t, tnew, la)
      call tnew%pull()
      call rccsd_inter_energy(t, scr, es, ed)
      lag = es + ed
      do a = 1, f%ints%nv
         do i = 1, f%ints%no
            lag = lag + 2.0_dp*f%l%t1(i, a)*(f%ints%eps_o(i) - f%ints%eps_v(a))*(tnew%t1(i, a) - t%t1(i, a))
         end do
      end do
      do b = 1, f%ints%nv
         do a = 1, f%ints%nv
            do j = 1, f%ints%no
               do i = 1, f%ints%no
                  d = f%ints%eps_o(i) - f%ints%eps_v(a) + f%ints%eps_o(j) - f%ints%eps_v(b)
                  lag = lag + (2.0_dp*f%l%t2(i, j, a, b) - f%l%t2(i, j, b, a))*d*(tnew%t2(i, j, a, b) - t%t2(i, j, a, b))
               end do
            end do
         end do
      end do
   end subroutine lagrangian

   subroutine test_forms(error)
      !! Every ladder form and forced small blocks give the update of the default
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_lambda_update_options_t) :: opts
      integer :: k, form_used
      character(len=32) :: what

      call ensure_fixture(1, error)
      if (allocated(error) .or. .not. fix(1)%supported) return
      do k = 1, 5
         opts = trc_rccsd_lambda_update_options_t()
         select case (k)
         case (1)
            opts%ladder%form = RCCSD_LADDER_ASSEMBLED
            opts%ladder%cd_batch = 7
            what = "assembled ladder"
         case (2)
            opts%ladder%form = RCCSD_LADDER_DIRECT
            opts%ladder%p_batch = 5
            what = "direct ladder"
         case (3)
            opts%ladder%form = RCCSD_LADDER_HALF
            opts%ladder%d_batch = 1
            what = "half ladder, one d a pass"
         case (4)
            opts%inter%woovo_batch = 1
            what = "woovo block of 1"
         case (5)
            opts%inter%woovo_batch = 3
            what = "woovo block of 3"
         end select
         call compare(error, 1, opts, trim(what), form_used)
         if (allocated(error)) return
         ! the forced form must be the one that ran
         if (opts%ladder%form /= 0) then
            call check(error, form_used == opts%ladder%form, trim(what)//": the form that ran is "//to_string(form_used))
            if (allocated(error)) return
         end if
      end do
   end subroutine test_forms

   subroutine compare(error, which, opts, what, form_used)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: which
      type(trc_rccsd_lambda_update_options_t), intent(in) :: opts
      character(len=*), intent(in) :: what
      integer, intent(out) :: form_used

      real(dp), allocatable :: new1(:, :), new2(:, :, :, :)
      real(dp) :: d1, d2

      call run_update(error, fix(which), fix(which)%amps, opts, new1, new2, form_used)
      if (allocated(error)) return
      d1 = maxval(abs(new1 - fix(which)%new1))/max(1.0_dp, maxval(abs(new1)))
      d2 = maxval(abs(new2 - fix(which)%new2))/max(1.0_dp, maxval(abs(new2)))
      print "(a,a,a,2es10.2)", "  ", what, ": rel. diff l1, l2 ", d1, d2
      call check(error, d1 < TOL_FORMS .and. d2 < TOL_FORMS, what//" differs from the default: "// &
                 to_string(d1, "es9.2")//" "//to_string(d2, "es9.2"))
   end subroutine compare

   subroutine run_update(error, f, amps, opts, new1, new2, form_used)
      !! terco's update of the fixture's multipliers at the amplitudes `amps`
      type(error_type), allocatable, intent(out) :: error
      type(fixture_t), intent(in) :: f
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_lambda_update_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: new1(:, :), new2(:, :, :, :)
      integer, intent(out) :: form_used

      type(trc_rccsd_lambda_update_t) :: upd
      type(trc_rccsd_amps_t) :: lnew

      form_used = 0
      call upd%init(f%ints, amps, la, opts)
      call check(error, .not. upd%status%has_error(), "init: "//upd%status%get_message())
      if (.not. allocated(error)) then
         form_used = upd%ls%form
         call lnew%init(f%ints%no, f%ints%nv)
         call upd%step(f%ints, amps, f%l, lnew, la)
         call check(error, .not. upd%status%has_error(), "step: "//upd%status%get_message())
      end if
      if (.not. allocated(error)) then
         call lnew%pull()
         new1 = lnew%t1
         new2 = lnew%t2
      end if
      call upd%release()
      call lnew%release()
   end subroutine run_update

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_lambda_update_t) :: upd
      type(trc_rccsd_amps_t) :: lnew, none, small
      type(trc_rccsd_ints_t) :: x, empty
      integer :: no, nv

      call ensure_fixture(1, error)
      if (allocated(error) .or. .not. fix(1)%supported) return
      associate (f => fix(1))
         no = f%ints%no
         nv = f%ints%nv
         call lnew%init(no, nv)
         call small%init(no, nv - 1)
         call upd%step(f%ints, f%amps, f%l, lnew, la)
         call expect(error, upd, "step before init")
         if (allocated(error)) return
         call upd%init(empty, f%amps, la)
         call expect(error, upd, "the integral blocks are not built")
         if (allocated(error)) return
         call upd%init(f%ints, none, la)
         call expect(error, upd, "the CCSD amplitudes are not initialised")
         if (allocated(error)) return
         call upd%init(f%ints, small, la)
         call expect(error, upd, "(no, nv) differ")
         if (allocated(error)) return
         call device_memory_set_cap(1000_int64)
         call upd%init(f%ints, f%amps, la)
         call device_memory_set_cap(0_int64)
         call expect(error, upd, "GB are free")
         if (allocated(error)) return
         call check(error, .not. allocated(upd%theta) .and. .not. allocated(upd%im%woooo), &
                    "a refused init must leave nothing allocated")
         if (allocated(error)) return

         ! integrals of its own, rebuilt after init
         call build_ints(f, x, FIX_FROZEN(1))
         call upd%init(x, f%amps, la)
         call check(error, .not. upd%status%has_error(), "init: "//upd%status%get_message())
         if (allocated(error)) return
         call upd%step(x, f%amps, small, lnew, la)
         call expect(error, upd, "(no, nv) of the amplitudes differ")
         if (allocated(error)) return
         call build_ints(f, x, FIX_FROZEN(1))
         call upd%step(x, f%amps, f%l, lnew, la)
         call expect(error, upd, "the integrals were rebuilt after init")
         call upd%release()
         call x%release()
         call lnew%release()
         call small%release()
      end associate
   end subroutine test_refusals

   subroutine expect(error, upd, fragment)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_lambda_update_t), intent(in) :: upd
      character(len=*), intent(in) :: fragment

      call check(error, upd%status%has_error() .and. index(upd%status%get_message(), fragment) > 0, &
                 "expected a refusal containing '"//fragment//"', got '"//upd%status%get_message()//"'")
   end subroutine expect

   subroutine build_ints(f, x, nfrozen)
      type(fixture_t), intent(in) :: f
      type(trc_rccsd_ints_t), intent(inout) :: x
      integer, intent(in) :: nfrozen
      type(trc_rccsd_ints_options_t) :: bopts

      bopts%nalpha = f%nocc
      bopts%nbeta = f%nocc
      bopts%nfrozen = nfrozen
      call x%build(f%bas, f%aux, f%pl, f%res%cmo(:, :, 1), f%res%eps(:, 1), la, bopts)
      if (x%status%has_error()) error stop "check_rccsd_lambda_update: build: "//x%status%get_message()
   end subroutine build_ints

   subroutine random_symmetric(x, s)
      !! s uniform in [-0.5, 0.5] with s(j,i,b,a) = s(i,j,a,b); x is its work array
      real(dp), intent(out) :: x(:, :, :, :), s(:, :, :, :)
      integer :: i, j, a, b

      call random_number(x)
      do b = 1, size(x, 4)
         do a = 1, size(x, 3)
            do j = 1, size(x, 2)
               do i = 1, size(x, 1)
                  s(i, j, a, b) = 0.5_dp*(x(i, j, a, b) + x(j, i, b, a)) - 0.5_dp
               end do
            end do
         end do
      end do
   end subroutine random_symmetric

   subroutine ensure_fixture(which, error)
      !! Water: the integrals, random amplitudes and multipliers, and the default update of them
      integer, intent(in) :: which
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_lambda_update_options_t) :: opts
      integer, allocatable :: zint(:), seed(:)
      real(dp), allocatable :: x(:, :, :, :)
      integer :: nseed, dev, no, nv, form_used

      if (fix(which)%built) return
      fix(which)%built = .true.
      if (.not. bound) then
         bound = .true.
         dev = trc_bind_device(0)
         call la%init(2)
      end if
      call load_molecule(1, fix(which)%bas, fix(which)%aux, fix(which)%pl, zint, fix(which)%nocc, fix(which)%supported)
      if (.not. fix(which)%supported) then
         print "(a)", "SKIP check_rccsd_lambda_update: this build's density-fitting kernels stop below the RIFIT shells"
         return
      end if
      call run_scf(fix(which)%bas, fix(which)%nocc, fix(which)%res)
      call build_ints(fix(which), fix(which)%ints, FIX_FROZEN(which))
      no = fix(which)%ints%no
      nv = fix(which)%ints%nv
      call random_seed(size=nseed)
      allocate (seed(nseed))
      seed = 20261127 + which
      call random_seed(put=seed)
      allocate (x(no, no, nv, nv))
      associate (amps => fix(which)%amps, l => fix(which)%l)
         call amps%init(no, nv)
         call l%init(no, nv)
         call random_number(amps%t1)
         amps%t1 = 0.4_dp*(amps%t1 - 0.5_dp)
         call random_symmetric(x, amps%t2)
         amps%t2 = 0.2_dp*amps%t2
         call random_number(l%t1)
         l%t1 = 0.4_dp*(l%t1 - 0.5_dp)
         call random_symmetric(x, l%t2)
         l%t2 = 0.2_dp*l%t2
         call amps%push()
         call l%push()
      end associate
      call run_update(error, fix(which), fix(which)%amps, opts, fix(which)%new1, fix(which)%new2, form_used)
   end subroutine ensure_fixture

   subroutine release_fixture()
      integer :: n

      ! every release is safe on what was never built
      do n = 1, NFIX
         call fix(n)%amps%release()
         call fix(n)%l%release()
         call fix(n)%ints%release()
         call fix(n)%pl%release()
         call fix(n)%aux%release()
         call fix(n)%bas%release()
      end do
      call la%release()
   end subroutine release_fixture

end module check_rccsd_lambda_update_suite

program check_rccsd_lambda_update
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_lambda_update_suite, only: collect_check_rccsd_lambda_update, release_fixture
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_lambda_update, error_unit, stat, parallel=.false.)
   call release_fixture()
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_lambda_update
