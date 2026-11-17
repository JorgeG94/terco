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
   !! SCF. The basis comes from the command line. Each reduction is run on
   !! water and on h2o2: the pair kernels skip one-centre quartets, so a
   !! molecule with one heavy atom never runs the classes whose four shells
   !! all have l >= 1; two heavy atoms with d shells do.
   !! (TRC_RYS_CLASSES=none, set in CMake, only chooses the Fock path of the
   !! fock_nosym oracle; the gradient is always Rys.)
   !!
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
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
   real(dp), parameter :: AGREE = 1.0e-13_dp
   !! Two contractions of the same derivative integrals, summed in a different
   !! order, relative to max(1, max|component|) (see reldev).
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

   subroutine setup(xyz, bas, eri, nat, at_z, r, error)
      !! The molecule in `xyz` in the command-line basis, eri built; skips a basis above LMAX.
      character(*), intent(in) :: xyz
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
      call read_xyz(xyz, nat, at_z, r)
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

   real(dp) function reldev(g, g0)
      !! max |g - g0| per unit of the larger of one and max |g0|: the
      !! summation order differs, so the round-off follows the size of the
      !! components (a few hundred on h2o2), not an absolute scale.
      real(dp), intent(in) :: g(:, :), g0(:, :)
      reldev = maxval(abs(g - g0))/max(1.0_dp, maxval(abs(g0)))
   end function reldev

   subroutine report(tag, dev)
      character(*), intent(in) :: tag
      real(dp), intent(in) :: dev
      print '(a,a,a,es10.2)', '  ', tag, ': max |difference|/scale ', dev
   end subroutine report

   subroutine run_rhf(xyz, error)
      character(*), intent(in) :: xyz
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d(:, :, :), g0(:, :), g1(:, :, :)

      call setup(xyz, bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (d(bas%nao, bas%nao, 1), g0(3, nat), g1(3, nat, 1))
      call synth_sym(bas%nao, 1, d(:, :, 1))
      call eri%gradient(bas, d, g0)
      call eri%gradient_pairs(bas, d, [pair(1, 1, 1, 1.0_dp, 0.5_dp)], 1, g1, err)
      call check(error, .not. err%has_error(), "gradient_pairs refused a valid call")
      if (.not. allocated(error)) then
         call report("RHF "//xyz, reldev(g1(:, :, 1), g0))
         call check(error, reldev(g1(:, :, 1), g0) <= AGREE, "RHF pair differs from eri%gradient by " &
                    //to_string(reldev(g1(:, :, 1), g0), "es10.2"))
      end if
      call eri%release()
      call bas%release()
   end subroutine run_rhf

   subroutine test_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call run_rhf('water.xyz', error)
      if (allocated(error)) return
      call run_rhf('h2o2.xyz', error)
   end subroutine test_rhf

   subroutine run_split(xyz, error)
      character(*), intent(in) :: xyz
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d(:, :, :), g0(:, :), g1(:, :, :), pool(:, :, :)

      call setup(xyz, bas, eri, nat, at_z, r, error)
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
         call report("split "//xyz, reldev(g1(:, :, 1), g0))
         call check(error, reldev(g1(:, :, 1), g0) <= AGREE, "split density differs from the whole by " &
                    //to_string(reldev(g1(:, :, 1), g0), "es10.2"))
      end if
      call eri%release()
      call bas%release()
   end subroutine run_split

   subroutine test_split(error)
      type(error_type), allocatable, intent(out) :: error
      call run_split('water.xyz', error)
      if (allocated(error)) return
      call run_split('h2o2.xyz', error)
   end subroutine test_split

   subroutine run_uhf(xyz, error)
      character(*), intent(in) :: xyz
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d2(:, :, :), g0(:, :), g1(:, :, :), pool(:, :, :)

      call setup(xyz, bas, eri, nat, at_z, r, error)
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
         call report("UHF "//xyz, reldev(g1(:, :, 1), g0))
         call check(error, reldev(g1(:, :, 1), g0) <= AGREE, "UHF pairs differ from eri%gradient by " &
                    //to_string(reldev(g1(:, :, 1), g0), "es10.2"))
      end if
      call eri%release()
      call bas%release()
   end subroutine run_uhf

   subroutine test_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      call run_uhf('water.xyz', error)
      if (allocated(error)) return
      call run_uhf('h2o2.xyz', error)
   end subroutine test_uhf

   subroutine run_columns(xyz, error)
      !! Six columns (two sweeps of RYSG_PAIRS_NCB) with different pair sets,
      !! against the six one-column calls.
      character(*), intent(in) :: xyz
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

      call setup(xyz, bas, eri, nat, at_z, r, error)
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
         dev = max(dev, reldev(g1(:, :, 1), gall(:, :, c)))
      end do
      call report("columns "//xyz, dev)
      call check(error, dev <= AGREE, "a column of the batch differs from its own call by "//to_string(dev, "es10.2"))
      call eri%release()
      call bas%release()
   end subroutine run_columns

   subroutine test_columns(error)
      type(error_type), allocatable, intent(out) :: error
      call run_columns('water.xyz', error)
      if (allocated(error)) return
      call run_columns('h2o2.xyz', error)
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

   real(dp) function fd_deriv(basis_json, nat, at_z, r, nao, ia, d, hstep, a1, a2, error)
      !! Central difference of ek_at along coordinate (d, ia) with step hstep.
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nao, ia, d
      real(dp), intent(in) :: r(3, nat), hstep, a1(nao, nao), a2(nao, nao)
      type(error_type), allocatable, intent(inout) :: error
      real(dp) :: rp(3, nat), rm(3, nat), ep, em

      fd_deriv = 0.0_dp
      rp = r; rm = r
      rp(d, ia) = rp(d, ia) + hstep
      rm(d, ia) = rm(d, ia) - hstep
      ep = ek_at(basis_json, nat, at_z, rp, nao, a1, a2, error)
      if (allocated(error)) return
      em = ek_at(basis_json, nat, at_z, rm, nao, a1, a2, error)
      if (allocated(error)) return
      fd_deriv = (ep - em)/(2.0_dp*hstep)
   end function fd_deriv

   subroutine test_antisym_fd(error)
      !! d/dR [Tr(A1 K[A2]) + Tr(A2 K[A1])] for antisymmetric A1, A2 (the pair
      !! (1, 2, 1, -1) counts both orderings), and Tr(A1 K[A1]) alone for the
      !! diagonal pair, against Richardson-extrapolated central differences,
      !! (4 D(h/2) - D(h))/3, which removes the h^2 truncation term. The pairs
      !! carry coulomb = 1: the functional is exchange only, so a Coulomb
      !! contribution leaking from an antisymmetric member would show. A
      !! third column holds a mixed symmetric/antisymmetric pair, which must
      !! add exactly zero.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      character(len=256) :: basis_json
      integer :: nat, nao, ia, d
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), pool(:, :, :), g(:, :, :)
      real(dp) :: fd, fd1, ref, worst, worst_plain

      call setup('water.xyz', bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      call get_command_argument(1, basis_json)
      nao = bas%nao
      allocate (pool(nao, nao, 3), g(3, nat, 3))
      call synth_anti(nao, 1, pool(:, :, 1))
      call synth_anti(nao, 2, pool(:, :, 2))
      call synth_sym(nao, 3, pool(:, :, 3))
      ! Column 1: the cross pair. Column 2: the diagonal pair of A1 (so the
      ! functional is Tr(A1 K[A1]) = half of E(A1, A1)). Column 3: a mixed
      ! pair of the symmetric density with A1, whose J and K both vanish.
      call eri%gradient_pairs(bas, pool, [pair(1, 2, 1, 1.0_dp, -1.0_dp), pair(1, 1, 2, 1.0_dp, -1.0_dp), &
                                          pair(3, 1, 3, 1.0_dp, -1.0_dp)], 3, g, err)
      call check(error, .not. err%has_error(), "gradient_pairs refused an antisymmetric pool")
      call eri%release()
      call bas%release()
      if (allocated(error)) return
      call check(error, all(g(:, :, 3) == 0.0_dp), "a mixed symmetric/antisymmetric pair did not add exactly zero")
      if (allocated(error)) return

      worst = 0.0_dp
      worst_plain = 0.0_dp
      do ia = 1, nat
         do d = 1, 3
            ! A pair with exchange = -1 differentiates +Tr(Dl K[Dr]).
            fd1 = fd_deriv(trim(basis_json), nat, at_z, r, nao, ia, d, H, pool(:, :, 1), pool(:, :, 2), error)
            if (allocated(error)) return
            fd = fd_deriv(trim(basis_json), nat, at_z, r, nao, ia, d, 0.5_dp*H, pool(:, :, 1), pool(:, :, 2), error)
            if (allocated(error)) return
            ref = g(d, ia, 1)
            worst_plain = max(worst_plain, abs(fd1 - ref)/max(abs(ref), abs(fd1), 1.0e-10_dp))
            fd = (4.0_dp*fd - fd1)/3.0_dp
            worst = max(worst, abs(fd - ref)/max(abs(ref), abs(fd), 1.0e-10_dp))
            ! ek_at doubled the diagonal term: Tr(A1 K[A1]) twice.
            fd1 = 0.5_dp*fd_deriv(trim(basis_json), nat, at_z, r, nao, ia, d, H, pool(:, :, 1), pool(:, :, 1), error)
            if (allocated(error)) return
            fd = 0.5_dp*fd_deriv(trim(basis_json), nat, at_z, r, nao, ia, d, 0.5_dp*H, pool(:, :, 1), pool(:, :, 1), error)
            if (allocated(error)) return
            ref = g(d, ia, 2)
            worst_plain = max(worst_plain, abs(fd1 - ref)/max(abs(ref), abs(fd1), 1.0e-10_dp))
            fd = (4.0_dp*fd - fd1)/3.0_dp
            worst = max(worst, abs(fd - ref)/max(abs(ref), abs(fd), 1.0e-10_dp))
         end do
      end do
      print '(a,es10.2,a,es10.2)', '  antisymmetric K: worst relative FD deviation, Richardson ', worst, &
         '  plain central ', worst_plain
      call check(error, worst <= RELTOL, "antisymmetric exchange gradient off finite differences by " &
                 //to_string(worst, "es10.2"))
   end subroutine test_antisym_fd

   subroutine refused(error, bas, eri, pool, pairs, ncol, what, omega)
      !! One refusal: gradient_pairs must set err and leave grad exactly zero,
      !! whatever grad held before the call.
      type(error_type), allocatable, intent(inout) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), intent(in) :: pool(:, :, :)
      type(trc_grad_pair_t), intent(in) :: pairs(:)
      integer, intent(in) :: ncol
      character(*), intent(in) :: what
      real(dp), intent(in), optional :: omega
      type(error_t) :: err
      real(dp), allocatable :: g(:, :, :)

      if (allocated(error)) return
      allocate (g(3, bas%natm, max(1, ncol)))
      g = 1.0_dp
      call eri%gradient_pairs(bas, pool, pairs, ncol, g, err, omega=omega)
      call check(error, err%has_error(), what//" was accepted")
      if (.not. allocated(error)) call check(error, all(g(:, :, 1:max(0, ncol)) == 0.0_dp), what//" left a nonzero gradient")
   end subroutine refused

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), pool(:, :, :), bad(:, :, :)
      real(dp) :: nan, inf

      call setup('water.xyz', bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      nan = ieee_value(1.0_dp, ieee_quiet_nan)
      inf = ieee_value(1.0_dp, ieee_positive_inf)
      allocate (pool(bas%nao, bas%nao, 1))
      call synth_sym(bas%nao, 1, pool(:, :, 1))
      call refused(error, bas, eri, pool, [pair(2, 1, 1, 0.5_dp, 0.25_dp)], 1, "a left index past the pool")
      call refused(error, bas, eri, pool, [pair(1, 0, 1, 0.5_dp, 0.25_dp)], 1, "a right index of zero")
      call refused(error, bas, eri, pool, [pair(1, 1, 2, 0.5_dp, 0.25_dp)], 1, "a column past ncol")
      call refused(error, bas, eri, pool, [pair(1, 1, 1, nan, 0.25_dp)], 1, "a NaN Coulomb weight")
      call refused(error, bas, eri, pool, [pair(1, 1, 1, 0.5_dp, inf)], 1, "an infinite exchange weight")
      call refused(error, bas, eri, pool, [pair(1, 1, 1, 0.5_dp, 0.25_dp)], 1, "a negative omega", omega=-0.1_dp)
      call refused(error, bas, eri, pool, [pair(1, 1, 1, 0.5_dp, 0.25_dp)], 1, "a NaN omega", omega=nan)
      call refused(error, bas, eri, pool, [pair(1, 1, 1, 0.5_dp, 0.25_dp)], 0, "ncol = 0")
      call refused(error, bas, eri, pool, [pair(1, 1, 1, 0.5_dp, 0.25_dp)], -1, "ncol < 0")
      if (.not. allocated(error)) then
         ! Wrong pool shapes: nao - 1 rows, and a non-square matrix.
         allocate (bad(bas%nao - 1, bas%nao - 1, 1))
         bad = 0.0_dp
         call refused(error, bas, eri, bad, [pair(1, 1, 1, 0.5_dp, 0.25_dp)], 1, "a pool of the wrong shape")
         deallocate (bad)
         allocate (bad(bas%nao, bas%nao - 1, 1))
         bad = 0.0_dp
         call refused(error, bas, eri, bad, [pair(1, 1, 1, 0.5_dp, 0.25_dp)], 1, "a non-square pool")
         deallocate (bad)
      end if
      if (.not. allocated(error)) then
         ! A single NaN, off the diagonal and so not the maximum of anything
         ! a maxval would be asked for first; then an infinity; then a NaN on
         ! the diagonal.
         bad = pool
         bad(1, 2, 1) = nan
         call refused(error, bas, eri, bad, [pair(1, 1, 1, 0.5_dp, 0.25_dp)], 1, "a pool with one off-diagonal NaN")
         bad = pool
         bad(2, 1, 1) = inf
         call refused(error, bas, eri, bad, [pair(1, 1, 1, 0.5_dp, 0.25_dp)], 1, "a pool with an infinity")
         bad = pool
         bad(1, 1, 1) = nan
         call refused(error, bas, eri, bad, [pair(1, 1, 1, 0.5_dp, 0.25_dp)], 1, "a pool with a diagonal NaN")
      end if
      if (.not. allocated(error)) then
         ! Neither symmetric nor antisymmetric.
         pool(1, 2, 1) = pool(1, 2, 1) + 0.5_dp
         call refused(error, bas, eri, pool, [pair(1, 1, 1, 1.0_dp, 0.5_dp)], 1, &
                      "a general (non-symmetric, non-antisymmetric) density")
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

      call setup('water.xyz', bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (d2(bas%nao, bas%nao, 2), pool(bas%nao, bas%nao, 3), g0(3, nat), g1(3, nat, 1))
      call synth_sym(bas%nao, 1, d2(:, :, 1))
      call synth_sym(bas%nao, 2, d2(:, :, 2))
      call eri%gradient(bas, d2(:, :, 1:1), g0)
      call eri%gradient_pairs(bas, d2(:, :, 1:1), [pair(1, 1, 1, 1.0_dp, 0.55_dp)], 1, g1, err)
      dev_rhf = reldev(g1(:, :, 1), g0)
      pool(:, :, 1) = d2(:, :, 1) + d2(:, :, 2)
      pool(:, :, 2:3) = d2
      call eri%gradient(bas, d2, g0)
      call eri%gradient_pairs(bas, pool, [pair(1, 1, 1, 1.0_dp, 0.0_dp), pair(2, 2, 1, 0.0_dp, 1.1_dp), &
                                          pair(3, 3, 1, 0.0_dp, 1.1_dp)], 1, g1, err)
      dev_uhf = reldev(g1(:, :, 1), g0)
      print '(a,es10.2,a,es10.2)', '  10% exchange error: RHF deviation ', dev_rhf, '  UHF deviation ', dev_uhf
      call check(error, dev_rhf > 1.0e-3_dp, "a 10% exchange error went unseen in the RHF comparison")
      if (.not. allocated(error)) call check(error, dev_uhf > 1.0e-3_dp, &
                                             "a 10% exchange error went unseen in the UHF comparison")
      call eri%release()
      call bas%release()
   end subroutine test_negative_control

   subroutine run_omega(xyz, error)
      character(*), intent(in) :: xyz
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d(:, :, :), g0(:, :), g1(:, :, :)

      call setup(xyz, bas, eri, nat, at_z, r, error)
      if (allocated(error)) return
      allocate (d(bas%nao, bas%nao, 1), g0(3, nat), g1(3, nat, 1))
      call synth_sym(bas%nao, 1, d(:, :, 1))
      call eri%gradient(bas, d, g0, omega=0.4_dp)
      call eri%gradient_pairs(bas, d, [pair(1, 1, 1, 1.0_dp, 0.5_dp)], 1, g1, err, omega=0.4_dp)
      call check(error, .not. err%has_error(), "gradient_pairs refused omega > 0")
      if (.not. allocated(error)) then
         call report("omega "//xyz, reldev(g1(:, :, 1), g0))
         call check(error, reldev(g1(:, :, 1), g0) <= AGREE, "omega pair differs from eri%gradient by " &
                    //to_string(reldev(g1(:, :, 1), g0), "es10.2"))
      end if
      call eri%release()
      call bas%release()
   end subroutine run_omega

   subroutine test_omega(error)
      type(error_type), allocatable, intent(out) :: error
      call run_omega('water.xyz', error)
      if (allocated(error)) return
      call run_omega('h2o2.xyz', error)
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
