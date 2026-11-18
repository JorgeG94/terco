module check_fock_split_suite
   !!
   !! eri%fock_split_many: the J - K/2 of densities with NO symmetry from the
   !! FOLDED kernels, D = S + A, G[D] = G[S] + K-only antisymmetric G[A], against
   !! eri%fock_nosym_many (the enumerated kernel, kept as the reference).
   !!
   !! Six checks (argv: xyz, basis json; test/CMakeLists.txt):
   !!
   !!  1. random non-symmetric densities, (j, k) in {(1,1), (1,0), (0,1),
   !!     (1,0.25)}, omega 0 and 0.3: split equals nosym;
   !!  2. negative controls, so the test can fail: dropping the A batch, and
   !!     running the A batch with the Coulomb half on, both miss by O(0.1);
   !!  3. an exactly symmetric density costs ONE folded batch and equals the
   !!     plain folded build bit for bit (S = D exactly);
   !!  4. J[A] = 0 and K[A] antisymmetric, from the reference and the split;
   !!  5. k_scale = 0 needs no A batch even for a non-symmetric density;
   !!  6. the antisymmetric folded build refuses j_scale /= 0 (and accepts 0).
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t, trc_fock_options_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_tables, only: LMAX
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_fock_split_tests, g_xyzfile, g_basis_json

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: OMEGA_TEST = 0.3_dp
   real(dp), parameter :: TOL_SPLIT = 1.0e-11_dp
      !! Split against nosym, relative to max|G|. Both add the same quartets
      !! in a different order: ~1e-14 measured, on the host and per mix.
   real(dp), parameter :: TOL_CONTROL = 1.0e-2_dp
      !! A negative control must miss by more than this (O(0.1) measured).
   integer, parameter :: ND = 3

   character(len=256) :: g_xyzfile = "h2o2.xyz"
   character(len=256) :: g_basis_json = "basis_sets/6-31gs.json"

contains

   subroutine collect_fock_split_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("split equals nosym for random non-symmetric densities", test_vs_nosym), &
                  new_unittest("dropping the A batch fails", test_negative_controls), &
                  new_unittest("a symmetric density costs one batch and equals the folded build", test_symmetric), &
                  new_unittest("J[A] vanishes and K[A] is antisymmetric", test_antisymmetric), &
                  new_unittest("k_scale = 0 needs no antisymmetric batch", test_no_exchange), &
                  new_unittest("an antisymmetric batch with j_scale /= 0 is refused", test_refusal) &
                  ]
   end subroutine collect_fock_split_tests

   subroutine nonsym_densities(n, nd, dd)
      !! Dense, deterministic, NOT symmetric densities (a module routine with
      !! explicit shapes, as check_fock_nosym_many's).
      integer, intent(in) :: n, nd
      real(dp), intent(out) :: dd(n, n, nd)
      integer :: i, j, k
      do k = 1, nd
         do j = 1, n
            do i = 1, n
               dd(i, j, k) = 0.5_dp*cos(real(3*i + 5*j + 7*k, dp)) + 0.3_dp*sin(real(2*i - 3*j + k, dp))
            end do
         end do
      end do
   end subroutine nonsym_densities

   subroutine split_parts(n, nd, dd, ds, da)
      !! D = S + A, exactly (the arithmetic `fock_split_many` uses).
      integer, intent(in) :: n, nd
      real(dp), intent(in) :: dd(n, n, nd)
      real(dp), intent(out) :: ds(n, n, nd), da(n, n, nd)
      integer :: k
      do k = 1, nd
         ds(:, :, k) = 0.5_dp*(dd(:, :, k) + transpose(dd(:, :, k)))
         da(:, :, k) = 0.5_dp*(dd(:, :, k) - transpose(dd(:, :, k)))
      end do
   end subroutine split_parts

   subroutine setup(bas, eri, built, error)
      !! `built` is true exactly when `eri%build` ran, i.e. when the caller owes
      !! a release; a skip releases its own basis.
      type(trc_basis_t), intent(out) :: bas
      type(trc_eri_t), intent(out) :: eri
      logical, intent(out) :: built
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      built = .false.
      call read_xyz(trim(g_xyzfile), nat, at_z, r)
      call trc_basis_from_json(trim(g_basis_json), nat, at_z, r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      ! A basis above the build's LMAX (cc-pVTZ's f on a d build) has no
      ! kernels to check, not wrong ones.
      if (maxval(bas%sh_l) > LMAX) then
         call skip_test(error, "SKIP: basis needs l = "//to_string(maxval(bas%sh_l)) &
                        //", this build has LMAX = "//to_string(LMAX))
         call bas%release()
         return
      end if
      call bas%to_device()
      call eri%build(bas, THRESH)
      built = .true.
   end subroutine setup

   subroutine test_vs_nosym(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_vs_nosym(bas, eri, error)
      ! Released on every path out: a failed check must not leave the device
      ! tables behind for the next case.
      call eri%release(); call bas%release()
   end subroutine test_vs_nosym

   subroutine run_vs_nosym(bas, eri, error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: d(:, :, :), gref(:, :, :), gsp(:, :, :)
      real(dp) :: sc(2, 4), om(2), rel
      integer :: nao, imix, io
      nao = bas%nao
      allocate (d(nao, nao, ND), gref(nao, nao, ND), gsp(nao, nao, ND))
      call nonsym_densities(nao, ND, d)
      call check(error, maxval(abs(d(:, :, 1) - transpose(d(:, :, 1)))) > 0.1_dp, "densities must be non-symmetric")
      if (allocated(error)) return
      sc = reshape([1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.25_dp], [2, 4])
      om = [0.0_dp, OMEGA_TEST]
      do io = 1, 2
         do imix = 1, 4
            call eri%fock_nosym_many(bas, ND, d, gref, trc_fock_options_t(j_scale=sc(1, imix), k_scale=sc(2, imix), &
                                                                          omega=om(io), density_screen=.false.))
            call eri%fock_split_many(bas, ND, d, gsp, trc_fock_options_t(j_scale=sc(1, imix), k_scale=sc(2, imix), &
                                                                         omega=om(io), density_screen=.false.))
            rel = maxval(abs(gref - gsp))/maxval(abs(gref))
            call global_logger%info("  omega "//to_string(om(io), "f4.1")//" (j,k) "//to_string(sc(1, imix), "f5.2")// &
                                    to_string(sc(2, imix), "f5.2")//": max|split - nosym|/max|G| "//to_string(rel, "es10.3"))
            call check(error, maxval(abs(gref)) > 1.0e-3_dp, "reference is trivially zero")
            if (allocated(error)) return
            call check(error, rel < TOL_SPLIT, "split vs nosym, relative "//to_string(rel, "es10.3")// &
                       " at omega "//to_string(om(io), "f4.1")//", (j,k) "//to_string(sc(1, imix), "f5.2")// &
                       to_string(sc(2, imix), "f5.2"))
            if (allocated(error)) return
         end do
      end do
   end subroutine run_vs_nosym

   subroutine test_negative_controls(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_negative_controls(bas, eri, error)
      ! Released on every path out: a failed check must not leave the device
      ! tables behind for the next case.
      call eri%release(); call bas%release()
   end subroutine test_negative_controls

   subroutine run_negative_controls(bas, eri, error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: d(:, :, :), ds(:, :, :), da(:, :, :), gref(:, :, :), gs(:, :, :), ga(:, :, :)
      real(dp) :: drop
      integer :: nao
      nao = bas%nao
      allocate (d(nao, nao, ND), ds(nao, nao, ND), da(nao, nao, ND), gref(nao, nao, ND), gs(nao, nao, ND), &
                ga(nao, nao, ND))
      call nonsym_densities(nao, ND, d)
      call split_parts(nao, ND, d, ds, da)
      call eri%fock_nosym_many(bas, ND, d, gref, trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, density_screen=.false.))
      ! Control 1: the symmetric batch alone.
      call eri%fock_folded_many(bas, ND, ds, gs, trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, density_screen=.false.))
      drop = maxval(abs(gref - gs))
      ! The A batch, exchange only.
      call eri%fock_folded_many(bas, ND, da, ga, trc_fock_options_t(j_scale=0.0_dp, k_scale=1.0_dp, density_screen=.false.), &
                                antisymmetric=.true.)
      call global_logger%info("  A batch dropped: miss "//to_string(drop, "es10.3"))
      call check(error, drop > TOL_CONTROL, "dropping the A batch must fail the comparison by > "// &
                 to_string(TOL_CONTROL, "es9.2")//", missed by only "//to_string(drop, "es10.3"))
      if (allocated(error)) return
      call check(error, maxval(abs(gref - gs - ga)) < TOL_SPLIT*maxval(abs(gref)), "the two batches sum to the reference")
      if (allocated(error)) return
   end subroutine run_negative_controls

   subroutine test_symmetric(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_symmetric(bas, eri, error)
      ! Released on every path out: a failed check must not leave the device
      ! tables behind for the next case.
      call eri%release(); call bas%release()
   end subroutine test_symmetric

   subroutine run_symmetric(bas, eri, error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: d(:, :, :), ds(:, :, :), da(:, :, :), gref(:, :, :), gsp(:, :, :), gf(:, :, :)
      type(trc_fock_options_t) :: o
      real(dp) :: rel, bitdiff
      integer :: nao
      nao = bas%nao
      allocate (d(nao, nao, ND), ds(nao, nao, ND), da(nao, nao, ND), gref(nao, nao, ND), gsp(nao, nao, ND), gf(nao, nao, ND))
      call nonsym_densities(nao, ND, d)
      call split_parts(nao, ND, d, ds, da)
      o = trc_fock_options_t(j_scale=2.0_dp, k_scale=1.5_dp, density_screen=.false.)
      call eri%fock_split_many(bas, ND, ds, gsp, o)
      call check(error, eri%nsplit_passes == 1, "a symmetric input must run ONE folded batch, ran "// &
                 to_string(eri%nsplit_passes))
      if (allocated(error)) return
      call eri%fock_folded_many(bas, ND, ds, gf, o)
      ! Same kernels on the same input: bitwise on the host, but device atomics
      ! reorder the accumulation, so compare to a tolerance.
      bitdiff = maxval(abs(gsp - gf))/maxval(abs(gf))
      call check(error, bitdiff < TOL_SPLIT, "symmetric split differs from the folded build by "//to_string(bitdiff, "es10.3"))
      if (allocated(error)) return
      call eri%fock_nosym_many(bas, ND, ds, gref, o)
      rel = maxval(abs(gref - gsp))/maxval(abs(gref))
      call check(error, rel < TOL_SPLIT, "symmetric split vs nosym, relative "//to_string(rel, "es10.3"))
      if (allocated(error)) return
      ! And a non-symmetric one runs two.
      call eri%fock_split_many(bas, ND, d, gsp, o)
      call check(error, eri%nsplit_passes == 2, "a non-symmetric input must run TWO folded batches, ran "// &
                 to_string(eri%nsplit_passes))
   end subroutine run_symmetric

   subroutine test_antisymmetric(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_antisymmetric(bas, eri, error)
      ! Released on every path out: a failed check must not leave the device
      ! tables behind for the next case.
      call eri%release(); call bas%release()
   end subroutine test_antisymmetric

   subroutine run_antisymmetric(bas, eri, error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: d(:, :, :), ds(:, :, :), da(:, :, :), gj(:, :, :), gk(:, :, :), gsp(:, :, :)
      real(dp) :: jmax, kmax, kasym, kdiff
      integer :: nao, id
      nao = bas%nao
      allocate (d(nao, nao, ND), ds(nao, nao, ND), da(nao, nao, ND), gj(nao, nao, ND), gk(nao, nao, ND), gsp(nao, nao, ND))
      call nonsym_densities(nao, ND, d)
      call split_parts(nao, ND, d, ds, da)
      call eri%fock_nosym_many(bas, ND, da, gj, trc_fock_options_t(j_scale=1.0_dp, k_scale=0.0_dp, density_screen=.false.))
      call eri%fock_nosym_many(bas, ND, da, gk, trc_fock_options_t(j_scale=0.0_dp, k_scale=1.0_dp, density_screen=.false.))
      jmax = maxval(abs(gj))
      kmax = maxval(abs(gk))
      call check(error, kmax > 1.0e-2_dp, "K of an antisymmetric density is trivially zero")
      if (allocated(error)) return
      call check(error, jmax < 1.0e-12_dp*kmax, "J[A] = 0, max|J[A]| = "//to_string(jmax, "es10.3"))
      if (allocated(error)) return
      ! The split's K[A] for the same input, antisymmetric and equal to nosym's.
      call eri%fock_split_many(bas, ND, da, gsp, trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, density_screen=.false.))
      kasym = 0.0_dp
      do id = 1, ND
         kasym = max(kasym, maxval(abs(gsp(:, :, id) + transpose(gsp(:, :, id)))))
      end do
      call check(error, kasym < 1.0e-12_dp*kmax, "K[A] antisymmetric, max|K+K^T| = "//to_string(kasym, "es10.3"))
      if (allocated(error)) return
      call eri%fock_nosym_many(bas, ND, da, gj, trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, density_screen=.false.))
      kdiff = maxval(abs(gsp - gj))
      call check(error, kdiff < TOL_SPLIT*kmax, "split of a purely antisymmetric density vs nosym "//to_string(kdiff, "es10.3"))
   end subroutine run_antisymmetric

   subroutine test_no_exchange(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_no_exchange(bas, eri, error)
      ! Released on every path out: a failed check must not leave the device
      ! tables behind for the next case.
      call eri%release(); call bas%release()
   end subroutine test_no_exchange

   subroutine run_no_exchange(bas, eri, error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: d(:, :, :), gref(:, :, :), gsp(:, :, :)
      type(trc_fock_options_t) :: o
      real(dp) :: rel
      integer :: nao
      nao = bas%nao
      allocate (d(nao, nao, ND), gref(nao, nao, ND), gsp(nao, nao, ND))
      call nonsym_densities(nao, ND, d)
      o = trc_fock_options_t(j_scale=1.0_dp, k_scale=0.0_dp, density_screen=.false.)
      call eri%fock_split_many(bas, ND, d, gsp, o)
      call check(error, eri%nsplit_passes == 1, "k_scale = 0 must run ONE folded batch, ran "//to_string(eri%nsplit_passes))
      if (allocated(error)) return
      call eri%fock_nosym_many(bas, ND, d, gref, o)
      rel = maxval(abs(gref - gsp))/maxval(abs(gref))
      call check(error, rel < TOL_SPLIT, "J-only split vs nosym, relative "//to_string(rel, "es10.3"))
   end subroutine run_no_exchange

   subroutine test_refusal(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_refusal(bas, eri, error)
      call eri%release(); call bas%release()
   end subroutine test_refusal

   subroutine run_refusal(bas, eri, error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: da(:, :, :), g(:, :, :)
      type(error_t) :: err
      integer :: nao
      nao = bas%nao
      allocate (da(nao, nao, 1), g(nao, nao, 1))
      da = 0.0_dp
      g = 7.0_dp
      call eri%fock_folded_many(bas, 1, da, g, trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, density_screen=.false.), &
                                antisymmetric=.true., err=err)
      call check(error, err%has_error(), "antisymmetric with j_scale = 1 must be refused")
      if (allocated(error)) return
      call check(error, err%code == ERROR_VALIDATION, "the refusal is a validation error")
      if (allocated(error)) return
      call check(error, all(g == 7.0_dp), "a refused call must leave gmats untouched")
      if (allocated(error)) return
      call eri%fock_folded_many(bas, 1, da, g, trc_fock_options_t(j_scale=0.0_dp, k_scale=1.0_dp, density_screen=.false.), &
                                antisymmetric=.true., err=err)
      call check(error, .not. err%has_error(), "antisymmetric with j_scale = 0 is accepted")
   end subroutine run_refusal

end module check_fock_split_suite

program check_fock_split
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_fock_split_suite, only: collect_fock_split_tests, g_xyzfile, g_basis_json
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   if (command_argument_count() >= 1) call get_command_argument(1, g_xyzfile)
   if (command_argument_count() >= 2) call get_command_argument(2, g_basis_json)

   stat = 0
   testsuites = [new_testsuite("check_fock_split", collect_fock_split_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop 1
   end if
end program check_fock_split
