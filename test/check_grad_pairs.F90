module check_grad_pairs_suite
   !!
   !! eri%gradient_pairs: one two-electron gradient PER COLUMN from a pool of
   !! densities and a table of separable pairs, every pool density going
   !! through one integral-derivative sweep.
   !!
   !! The oracles are eri%gradient (the RHF/UHF contraction the pair kernel
   !! must reduce to, which pyscf already anchors), the identity that a
   !! density split across the pool must give the same gradient as the whole,
   !! single-column calls against a many-column one, and -- for an
   !! ANTISYMMETRIC density, which eri%gradient cannot take -- a central
   !! finite difference of Tr(A K[A]) with A fixed and the integrals moving,
   !! built from eri%fock_nosym (the enumerated kernel that keeps the
   !! antisymmetric part of K), a different code path from the gradient.
   !!
   !! The densities are deterministic synthetic matrices, as in
   !! check_rys_grad_fd: this checks d(ab|cd)/dR and its contraction, not an
   !! SCF. The basis comes from the command line.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t
   use trc_tables, only: LMAX
   use trc_eri, only: trc_eri_t, trc_grad_pair_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_grad_pairs_tests

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: AGREE = 1.0e-12_dp
   !! Two contractions of the same derivative integrals, summed in a different
   !! order; gradient components are O(1).
   real(dp), parameter :: H = 1.0e-4_dp
   real(dp), parameter :: RELTOL = 1.0e-7_dp
   !! The finite-difference tolerance check_rys_grad_fd uses.

contains

   subroutine collect_grad_pairs_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF pair equals eri%gradient", test_rhf), &
                  new_unittest("a density split over the pool counts cross pairs", test_split), &
                  new_unittest("UHF pairs equal eri%gradient with two densities", test_uhf), &
                  new_unittest("columns equal separate calls, across sweeps", test_columns), &
                  new_unittest("antisymmetric exchange matches finite differences", test_antisym_fd), &
                  new_unittest("bad input is refused", test_refusals), &
                  new_unittest("a 10% exchange error is caught", test_negative_control), &
                  new_unittest("range-separated pairs equal eri%gradient with omega", test_omega) &
                  ]
   end subroutine collect_grad_pairs_tests

   subroutine setup(bas, eri, nat, at_z, r, error)
      !! Water in the command-line basis, eri built; skips a basis above LMAX.
      type(trc_basis_t), intent(out) :: bas
      type(trc_eri_t), intent(out) :: eri
      integer, intent(out) :: nat
      integer, allocatable, intent(out) :: at_z(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      type(error_t) :: err

      basis_json = ''
      if (command_argument_count() >= 1) call get_command_argument(1, basis_json)
      call check(error, len_trim(basis_json) > 0, "usage: check_grad_pairs <basis.json>")
      if (allocated(error)) return
      call read_xyz('water.xyz', nat, at_z, r)
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "check_grad_pairs: "//err%get_message())
      if (allocated(error)) return
      if (maxval(bas%sh_l) > LMAX) then
         call skip_test(error, "SKIP: basis needs l = "//to_string(maxval(bas%sh_l)))
         call bas%release()
         return
      end if
      call bas%to_device()
      call eri%build(bas, THRESH)
   end subroutine setup

   subroutine synth_sym(nao, phase, d)
      integer, intent(in) :: nao, phase
      real(dp), intent(out) :: d(nao, nao)
      integer :: i, j
      do j = 1, nao
         do i = 1, nao
            d(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j + 7*phase, dp))
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
   end subroutine synth_sym

   subroutine synth_anti(nao, phase, a)
      integer, intent(in) :: nao, phase
      real(dp), intent(out) :: a(nao, nao)
      integer :: i, j
      do j = 1, nao
         do i = 1, nao
            a(i, j) = 0.3_dp*sin(real(2*i + 3*j + 11*phase, dp)) + 0.1_dp*cos(real(i*j + phase, dp))
         end do
      end do
      a = 0.5_dp*(a - transpose(a))
   end subroutine synth_anti

   function pair(l, r, col, c, x) result(p)
      integer, intent(in) :: l, r, col
      real(dp), intent(in) :: c, x
      type(trc_grad_pair_t) :: p
      p = trc_grad_pair_t(left=l, right=r, column=col, coulomb=c, exchange=x)
   end function pair

   subroutine report(tag, dev)
      character(*), intent(in) :: tag
      real(dp), intent(in) :: dev
      print '(a,a,a,es10.2)', '  ', tag, ': max |difference| ', dev
   end subroutine report

   subroutine test_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d(:, :, :), g0(:, :), g1(:, :, :)

      call setup(bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (d(bas%nao, bas%nao, 1), g0(3, nat), g1(3, nat, 1))
      call synth_sym(bas%nao, 1, d(:, :, 1))
      call eri%gradient(bas, d, g0)
      call eri%gradient_pairs(bas, d, [pair(1, 1, 1, 1.0_dp, 0.5_dp)], 1, g1, err)
      call check(error, .not. err%has_error(), "gradient_pairs refused a valid call")
      if (.not. allocated(error)) then
         call report("RHF", maxval(abs(g1(:, :, 1) - g0)))
         call check(error, maxval(abs(g1(:, :, 1) - g0)) <= AGREE, "RHF pair differs from eri%gradient by " &
                    //to_string(maxval(abs(g1(:, :, 1) - g0)), "es10.2"))
      end if
      call eri%release()
      call bas%release()
   end subroutine test_rhf

   subroutine test_split(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d(:, :, :), g0(:, :), g1(:, :, :), pool(:, :, :)

      call setup(bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (d(bas%nao, bas%nao, 1), pool(bas%nao, bas%nao, 2), g0(3, nat), g1(3, nat, 1))
      call synth_sym(bas%nao, 1, d(:, :, 1))
      pool(:, :, 1) = 0.3_dp*d(:, :, 1)
      pool(:, :, 2) = 0.7_dp*d(:, :, 1)
      call eri%gradient(bas, d, g0)
      ! One cross entry stands for both orderings.
      call eri%gradient_pairs(bas, pool, [pair(1, 1, 1, 1.0_dp, 0.5_dp), pair(2, 2, 1, 1.0_dp, 0.5_dp), &
                                          pair(1, 2, 1, 1.0_dp, 0.5_dp)], 1, g1, err)
      call check(error, .not. err%has_error(), "gradient_pairs refused a valid call")
      if (.not. allocated(error)) then
         call report("split", maxval(abs(g1(:, :, 1) - g0)))
         call check(error, maxval(abs(g1(:, :, 1) - g0)) <= AGREE, "split density differs from the whole by " &
                    //to_string(maxval(abs(g1(:, :, 1) - g0)), "es10.2"))
      end if
      call eri%release()
      call bas%release()
   end subroutine test_split

   subroutine test_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d2(:, :, :), g0(:, :), g1(:, :, :), pool(:, :, :)

      call setup(bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (d2(bas%nao, bas%nao, 2), pool(bas%nao, bas%nao, 3), g0(3, nat), g1(3, nat, 1))
      call synth_sym(bas%nao, 1, d2(:, :, 1))
      call synth_sym(bas%nao, 2, d2(:, :, 2))
      pool(:, :, 1) = d2(:, :, 1) + d2(:, :, 2)
      pool(:, :, 2:3) = d2
      call eri%gradient(bas, d2, g0)
      ! Coulomb on the total density only, exchange per spin only.
      call eri%gradient_pairs(bas, pool, [pair(1, 1, 1, 1.0_dp, 0.0_dp), pair(2, 2, 1, 0.0_dp, 1.0_dp), &
                                          pair(3, 3, 1, 0.0_dp, 1.0_dp)], 1, g1, err)
      call check(error, .not. err%has_error(), "gradient_pairs refused a valid call")
      if (.not. allocated(error)) then
         call report("UHF", maxval(abs(g1(:, :, 1) - g0)))
         call check(error, maxval(abs(g1(:, :, 1) - g0)) <= AGREE, "UHF pairs differ from eri%gradient by " &
                    //to_string(maxval(abs(g1(:, :, 1) - g0)), "es10.2"))
      end if
      call eri%release()
      call bas%release()
   end subroutine test_uhf

   subroutine test_columns(error)
      !! Six columns (two sweeps of RYSG_PAIRS_NCB) with different pair sets,
      !! against the six one-column calls.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat, c
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), pool(:, :, :), gall(:, :, :), g1(:, :, :)
      type(trc_grad_pair_t) :: ps(9)
      type(trc_grad_pair_t), allocatable :: one(:)
      real(dp) :: dev

      call setup(bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (pool(bas%nao, bas%nao, 3), gall(3, nat, 6), g1(3, nat, 1))
      do c = 1, 3
         call synth_sym(bas%nao, c, pool(:, :, c))
      end do
      ps = [pair(1, 1, 1, 1.0_dp, 0.5_dp), pair(2, 2, 2, 1.0_dp, 0.0_dp), pair(3, 3, 2, 0.0_dp, 1.0_dp), &
            pair(1, 2, 3, 1.0_dp, 0.3_dp), pair(2, 3, 4, 0.2_dp, 0.7_dp), pair(1, 3, 5, 0.4_dp, 0.1_dp), &
            pair(3, 3, 5, 0.1_dp, 0.2_dp), pair(1, 1, 6, 0.0_dp, 0.0_dp), pair(2, 1, 6, 1.0_dp, 1.0_dp)]
      call eri%gradient_pairs(bas, pool, ps, 6, gall, err)
      call check(error, .not. err%has_error(), "gradient_pairs refused a valid call")
      if (allocated(error)) then
         call eri%release()
         call bas%release()
         return
      end if
      dev = 0.0_dp
      do c = 1, 6
         one = pack(ps, ps%column == c)
         one%column = 1
         call eri%gradient_pairs(bas, pool, one, 1, g1, err)
         dev = max(dev, maxval(abs(g1(:, :, 1) - gall(:, :, c))))
      end do
      call report("columns", dev)
      call check(error, dev <= AGREE, "a column of the batch differs from its own call by "//to_string(dev, "es10.2"))
      call eri%release()
      call bas%release()
   end subroutine test_columns

   real(dp) function ek_at(basis_json, nat, at_z, r, nao, a1, a2, error)
      !! E(R) = Tr(A1 K[A2]) + Tr(A2 K[A1]) at the geometry in r, with the
      !! matrices fixed in the AO basis and the basis functions moving.
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nao
      real(dp), intent(in) :: r(3, nat), a1(nao, nao), a2(nao, nao)
      type(error_type), allocatable, intent(inout) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      real(dp) :: k(nao, nao)

      ek_at = 0.0_dp
      call trc_basis_from_json(basis_json, nat, at_z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "check_grad_pairs: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call eri%build(bas, THRESH)
      ! j_scale = 0 and k_scale = -2 turn eri%fock_nosym's J - K/2 into +K.
      call eri%fock_nosym(bas, a2, k, k_scale=-2.0_dp, j_scale=0.0_dp, density_screen=.false.)
      ek_at = sum(transpose(a1)*k)
      call eri%fock_nosym(bas, a1, k, k_scale=-2.0_dp, j_scale=0.0_dp, density_screen=.false.)
      ek_at = ek_at + sum(transpose(a2)*k)
      call eri%release()
      call bas%release()
   end function ek_at

   subroutine test_antisym_fd(error)
      !! d/dR [Tr(A1 K[A2]) + Tr(A2 K[A1])] for antisymmetric A1, A2 (the pair
      !! (1, 2, 0, -1) counts both orderings), and Tr(A1 K[A1]) alone for the
      !! diagonal pair, against central differences.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      character(len=256) :: basis_json
      integer :: nat, nao, ia, d
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), rp(:, :), rm(:, :), pool(:, :, :), g(:, :, :)
      real(dp) :: ep, em, fd, ref, worst

      call setup(bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      call get_command_argument(1, basis_json)
      nao = bas%nao
      allocate (pool(nao, nao, 2), g(3, nat, 2), rp(3, nat), rm(3, nat))
      call synth_anti(nao, 1, pool(:, :, 1))
      call synth_anti(nao, 2, pool(:, :, 2))
      ! Column 1: the cross pair. Column 2: the diagonal pair of A1 (so the
      ! functional is Tr(A1 K[A1]) = half of E(A1, A1)).
      call eri%gradient_pairs(bas, pool, [pair(1, 2, 1, 0.0_dp, -1.0_dp), pair(1, 1, 2, 0.0_dp, -1.0_dp)], 2, g, err)
      call check(error, .not. err%has_error(), "gradient_pairs refused an antisymmetric pool")
      call eri%release()
      call bas%release()
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, nat
         do d = 1, 3
            rp = r; rm = r
            rp(d, ia) = rp(d, ia) + H
            rm(d, ia) = rm(d, ia) - H
            ep = ek_at(trim(basis_json), nat, at_z, rp, nao, pool(:, :, 1), pool(:, :, 2), error)
            if (allocated(error)) return
            em = ek_at(trim(basis_json), nat, at_z, rm, nao, pool(:, :, 1), pool(:, :, 2), error)
            if (allocated(error)) return
            ! The functional is Tr(A1 K[A2]) + Tr(A2 K[A1]); a pair with
            ! exchange = -1 differentiates +Tr(Dl K[Dr]).
            fd = (ep - em)/(2.0_dp*H)
            ref = g(d, ia, 1)
            worst = max(worst, abs(fd - ref)/max(abs(ref), abs(fd), 1.0e-10_dp))
            ep = ek_at(trim(basis_json), nat, at_z, rp, nao, pool(:, :, 1), pool(:, :, 1), error)
            if (allocated(error)) return
            em = ek_at(trim(basis_json), nat, at_z, rm, nao, pool(:, :, 1), pool(:, :, 1), error)
            if (allocated(error)) return
            ! ek_at doubled the diagonal term: Tr(A1 K[A1]) twice.
            fd = 0.5_dp*(ep - em)/(2.0_dp*H)
            ref = g(d, ia, 2)
            worst = max(worst, abs(fd - ref)/max(abs(ref), abs(fd), 1.0e-10_dp))
         end do
      end do
      print '(a,es10.2)', '  antisymmetric K: worst relative FD deviation ', worst
      call check(error, worst <= RELTOL, "antisymmetric exchange gradient off finite differences by " &
                 //to_string(worst, "es10.2"))
   end subroutine test_antisym_fd

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), pool(:, :, :), g(:, :, :)

      call setup(bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (pool(bas%nao, bas%nao, 1), g(3, nat, 1))
      call synth_sym(bas%nao, 1, pool(:, :, 1))
      call eri%gradient_pairs(bas, pool, [pair(2, 1, 1, 0.5_dp, 0.25_dp)], 1, g, err)
      call check(error, err%has_error(), "a left index past the pool was accepted")
      if (.not. allocated(error)) then
         call eri%gradient_pairs(bas, pool, [pair(1, 0, 1, 0.5_dp, 0.25_dp)], 1, g, err)
         call check(error, err%has_error(), "a right index of zero was accepted")
      end if
      if (.not. allocated(error)) then
         call eri%gradient_pairs(bas, pool, [pair(1, 1, 2, 0.5_dp, 0.25_dp)], 1, g, err)
         call check(error, err%has_error(), "a column past ncol was accepted")
      end if
      if (.not. allocated(error)) then
         ! Neither symmetric nor antisymmetric.
         pool(1, 2, 1) = pool(1, 2, 1) + 0.5_dp
         call eri%gradient_pairs(bas, pool, [pair(1, 1, 1, 1.0_dp, 0.5_dp)], 1, g, err)
         call check(error, err%has_error(), "a general (non-symmetric, non-antisymmetric) density was accepted")
      end if
      call eri%release()
      call bas%release()
   end subroutine test_refusals

   subroutine test_negative_control(error)
      !! A 10% error in the exchange weight must fail the RHF and UHF
      !! comparisons of the tests above.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d2(:, :, :), g0(:, :), g1(:, :, :), pool(:, :, :)
      real(dp) :: dev_rhf, dev_uhf

      call setup(bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (d2(bas%nao, bas%nao, 2), pool(bas%nao, bas%nao, 3), g0(3, nat), g1(3, nat, 1))
      call synth_sym(bas%nao, 1, d2(:, :, 1))
      call synth_sym(bas%nao, 2, d2(:, :, 2))
      call eri%gradient(bas, d2(:, :, 1:1), g0)
      call eri%gradient_pairs(bas, d2(:, :, 1:1), [pair(1, 1, 1, 1.0_dp, 0.55_dp)], 1, g1, err)
      dev_rhf = maxval(abs(g1(:, :, 1) - g0))
      pool(:, :, 1) = d2(:, :, 1) + d2(:, :, 2)
      pool(:, :, 2:3) = d2
      call eri%gradient(bas, d2, g0)
      call eri%gradient_pairs(bas, pool, [pair(1, 1, 1, 1.0_dp, 0.0_dp), pair(2, 2, 1, 0.0_dp, 1.1_dp), &
                                          pair(3, 3, 1, 0.0_dp, 1.1_dp)], 1, g1, err)
      dev_uhf = maxval(abs(g1(:, :, 1) - g0))
      print '(a,es10.2,a,es10.2)', '  10% exchange error: RHF deviation ', dev_rhf, '  UHF deviation ', dev_uhf
      call check(error, dev_rhf > 1.0e-3_dp, "a 10% exchange error went unseen in the RHF comparison")
      if (.not. allocated(error)) call check(error, dev_uhf > 1.0e-3_dp, &
                                             "a 10% exchange error went unseen in the UHF comparison")
      call eri%release()
      call bas%release()
   end subroutine test_negative_control

   subroutine test_omega(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d(:, :, :), g0(:, :), g1(:, :, :)

      call setup(bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (d(bas%nao, bas%nao, 1), g0(3, nat), g1(3, nat, 1))
      call synth_sym(bas%nao, 1, d(:, :, 1))
      call eri%gradient(bas, d, g0, omega=0.4_dp)
      call eri%gradient_pairs(bas, d, [pair(1, 1, 1, 1.0_dp, 0.5_dp)], 1, g1, err, omega=0.4_dp)
      call check(error, .not. err%has_error(), "gradient_pairs refused omega > 0")
      if (.not. allocated(error)) then
         call report("omega", maxval(abs(g1(:, :, 1) - g0)))
         call check(error, maxval(abs(g1(:, :, 1) - g0)) <= AGREE, "omega pair differs from eri%gradient by " &
                    //to_string(maxval(abs(g1(:, :, 1) - g0)), "es10.2"))
      end if
      call eri%release()
      call bas%release()
   end subroutine test_omega

end module check_grad_pairs_suite

program check_grad_pairs
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_grad_pairs_suite, only: collect_grad_pairs_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("check_grad_pairs", collect_grad_pairs_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_grad_pairs
