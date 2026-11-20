module check_fock_split_suite
   !!
   !! eri%fock_split_many: the J - K/2 of densities with NO symmetry from the
   !! FOLDED kernels, D = S + A, G[D] = G[S] + K-only antisymmetric G[A], against
   !! eri%fock_nosym_many (the enumerated kernel, kept as the reference).
   !!
   !! Fifteen checks (argv: xyz, basis json; test/CMakeLists.txt):
   !!
   !!  1. random non-symmetric densities, (j, k) in {(1,1), (1,0), (0,1),
   !!     (1,0.25)}, omega 0 and 0.3: split equals nosym;
   !!  2. negative controls, so the test can fail: dropping the A batch, and
   !!     running the A batch with the Coulomb half on, both miss by O(0.1);
   !!  3. an exactly symmetric density costs ONE folded batch and equals the
   !!     plain folded build bit for bit (S = D exactly);
   !!  4. J[A] = 0 and K[A] antisymmetric, from the reference and the split;
   !!  5. k_scale = 0 needs no A batch even for a non-symmetric density;
   !!  6. the antisymmetric folded build refuses j_scale /= 0 (and accepts 0);
   !!  7. `fock_split_many_lent` (the device-resident route the TDDFT matvec
   !!     takes) equals `fock_split_many` on several non-symmetric densities at
   !!     once, for (j, k) in {(1,1), (1,0), (0,1)}, an exactly symmetric batch
   !!     takes ONE pass, a mixed batch (one symmetric density, first or
   !!     last, among non-symmetric ones) runs the antisymmetric pass, omega > 0
   !!     and the density_screen fallback agree too; equal bit for bit on the
   !!     host build, within TOL_LENT on the device (_OPENACC) and in the threaded
   !!     OpenMP port (_OPENMP);
   !!  8. `fock_nosym_many_lent` equals `fock_nosym_many`.
   !!  9. the chunk rule (`eri_lent_chunk_size`): 3 densities per launch on the
   !!     two large systems it was measured on, the whole batch on a small one;
   !! 10. `fock_split_many_lent` in chunks (`eri%many_chunk`) equals the whole
   !!     batch and the host split, density for density, over a batch that
   !!     mixes exactly symmetric, round-off-asymmetric and non-symmetric
   !!     densities across the chunk boundaries (the antisymmetric pass is
   !!     decided over the batch, not per chunk), and the launch counter shows
   !!     the batch was split;
   !! 11. `fock_nosym_many_lent` in chunks equals the whole batch and the host;
   !! 12. density_screen on the device route, densities with exactly zero
   !!     blocks (including a density whose block (i,j) is large and (j,i)
   !!     zero): the screened result equals the unscreened one exactly (a
   !!     dropped quartet then reads only zero blocks), work is dropped, and it
   !!     equals the host's screened routines; split and nosym, whole and in
   !!     chunks;
   !! 13. density_screen with a real threshold: the screened device route
   !!     equals the host's screened route (one chunk), a chunked screened
   !!     build differs from the unscreened one by no more than the bound the
   !!     dropped quartets allow, and the bound is per launch (a tiny density
   !!     in a launch with a large one is not screened away);
   !! 14. density_screen off is the plain Schwarz bound whatever the chunking;
   !! 15. `fock_folded_many_lent` (the RPA matvec's route: one folded batch on
   !!     densities that are exactly symmetric, or exactly antisymmetric with
   !!     j_scale = 0) equals the enumerated `fock_nosym_many` within TOL_SPLIT
   !!     and the host `fock_folded_many` (bit for bit on the host build),
   !!     whole and in chunks, screen off and on, with k_scale = 0; a
   !!     non-symmetric density through it is wrong (negative control), and the
   !!     antisymmetric batch with j_scale /= 0 is refused.
   !! With a general-contraction basis (the cc-pVDZ variant) the working basis
   !! is decontracted, so 7 to 14 also cover the device expand and fold.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t, trc_fock_options_t, eri_lent_chunk_size
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
   real(dp), parameter :: TOL_LENT = 1.0e-13_dp
      !! Device-resident route against the host route: the same kernels and
      !! the same per-element arithmetic, so the difference is the order the
      !! accumulators were added in (atomics on the device).
   real(dp), parameter :: TOL_CONTROL = 1.0e-2_dp
      !! A negative control must miss by more than this (O(0.1) measured).
   integer, parameter :: ND = 3
   real(dp), parameter :: SCREEN_THRESH = 1.0e-8_dp
      !! The integral threshold of the real-threshold screening checks: loose enough that a
      !! density block of 1e-12 drops its quartets, tight enough to leave an O(1) density whole.

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
                  new_unittest("an antisymmetric batch with j_scale /= 0 is refused", test_refusal), &
                  new_unittest("the device-resident split equals the host split", test_lent_split), &
                  new_unittest("the device-resident nosym build equals the host one", test_lent_nosym), &
                  new_unittest("the chunk rule: 3 densities on the large systems, the batch on a small one", test_chunk_rule), &
                  new_unittest("the chunked split equals the whole batch and the host", test_chunked_split), &
                  new_unittest("the chunked nosym build equals the whole batch and the host", test_chunked_nosym), &
                  new_unittest("device screening, zero blocks: screened equals unscreened exactly", test_screen_zero_blocks), &
                  new_unittest("device screening, real threshold: host agreement, derived bound, per launch", &
                               test_screen_threshold), &
                  new_unittest("density_screen off ignores the chunking", test_screen_off_chunks), &
                  new_unittest("the device-resident folded build equals nosym and the host folded build", test_lent_folded) &
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

   subroutine setup(bas, eri, built, error, thr)
      !! `built` is true exactly when `eri%build` ran, i.e. when the caller owes
      !! a release; a skip releases its own basis. `thr` is the integral
      !! threshold (default THRESH, 1e-14: tight enough that the density screen
      !! of a dense density drops nothing).
      type(trc_basis_t), intent(out) :: bas
      type(trc_eri_t), intent(out) :: eri
      logical, intent(out) :: built
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in), optional :: thr
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
      if (present(thr)) then
         call eri%build(bas, thr)
      else
         call eri%build(bas, THRESH)
      end if
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

   subroutine test_lent_split(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_lent_split(bas, eri, error)
      call eri%release(); call bas%release()
   end subroutine test_lent_split

   subroutine run_lent_split(bas, eri, error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: d(:, :, :)
      real(dp) :: sc(2, 3)
      integer :: nao, imix, k
      nao = bas%nao
      allocate (d(nao, nao, ND))
      sc = reshape([1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
      ! k = 1 non-symmetric densities, k = 2 every density exactly symmetric
      ! (the antisymmetric pass is skipped), k = 3 and 4 a MIXED batch: the
      ! first (the last) density symmetric and the others not, so the decision to run the
      ! antisymmetric pass must be taken over the whole batch.
      do k = 1, 4
         call nonsym_densities(nao, ND, d)
         if (k == 2) call symmetrise_batch(nao, ND, d, 1, ND)
         if (k == 3) call symmetrise_batch(nao, ND, d, 1, 1)
         if (k == 4) call symmetrise_batch(nao, ND, d, ND, ND)
         do imix = 1, 3
            call compare_lent_split(bas, eri, d, trc_fock_options_t(j_scale=sc(1, imix), k_scale=sc(2, imix), &
                                                                    density_screen=.false.), &
                                    "batch "//to_string(k), error)
            if (allocated(error)) return
            if (k == 2) then
               call check(error, eri%nsplit_passes == 1, "an exactly symmetric batch must cost ONE pass")
               if (allocated(error)) return
            end if
            if (k >= 3 .and. sc(2, imix) /= 0.0_dp) then
               call check(error, eri%nsplit_passes == 2, "a mixed batch needs the antisymmetric pass")
               if (allocated(error)) return
            end if
         end do
      end do
      ! A range-separated operator.
      call nonsym_densities(nao, ND, d)
      call compare_lent_split(bas, eri, d, trc_fock_options_t(j_scale=0.0_dp, k_scale=1.0_dp, omega=OMEGA_TEST, &
                                                              density_screen=.false.), "omega", error)
      if (allocated(error)) return
      ! density_screen = .true. takes the host-route fallback (pull the lent
      ! densities down, push the result up): it must give what the host
      ! routine called the same way gives.
      call compare_lent_split(bas, eri, d, trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, density_screen=.true.), &
                              "density_screen", error)
   end subroutine run_lent_split

   subroutine compare_lent_split(bas, eri, d, o, label, error)
      !! `fock_split_many_lent` against `fock_split_many` for one batch and
      !! one set of options: equal bit for bit on the host build, within
      !! atomic-order noise on the device.
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), intent(in) :: d(bas%nao, bas%nao, ND)
      type(trc_fock_options_t), intent(in) :: o
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gref(:, :, :), glent(:, :, :)
      real(dp) :: rel
      integer :: nao, npass_ref
      nao = bas%nao
      allocate (gref(nao, nao, ND), glent(nao, nao, ND))
      call eri%fock_split_many(bas, ND, d, gref, o)
      npass_ref = eri%nsplit_passes
      glent = -99.0_dp
      !$acc enter data copyin(d) create(glent)
      call eri%fock_split_many_lent(bas, ND, d, glent, o)
      !$acc update self(glent)
      !$acc exit data delete(d, glent)
      call check(error, maxval(abs(gref)) > 1.0e-3_dp, "reference is trivially zero")
      if (allocated(error)) return
      rel = maxval(abs(gref - glent))/maxval(abs(gref))
      call global_logger%info("  lent split "//label//" (j,k,omega,screen) "//to_string(o%j_scale, "f5.2")// &
                              to_string(o%k_scale, "f5.2")//to_string(o%omega, "f5.2")//" "// &
                              merge("T", "F", o%density_screen)//": max|lent - host|/max|G| "//to_string(rel, "es10.3"))
      call check_same(error, gref, glent, "lent split vs host split, "//label)
      if (allocated(error)) return
      call check(error, eri%nsplit_passes == npass_ref, "lent split ran "//to_string(eri%nsplit_passes)// &
                 " passes, the host route "//to_string(npass_ref))
   end subroutine compare_lent_split


   subroutine lent_run(eri, bas, split, nd, d, o, chunk, g)
      !! One lent build of `d` (nd densities) at `chunk` densities per launch (0: the
      !! rule), through the device-resident route the matvec takes, the result on the host.
      type(trc_eri_t), intent(inout) :: eri
      type(trc_basis_t), intent(in) :: bas
      logical, intent(in) :: split
      integer, intent(in) :: nd, chunk
      real(dp), intent(in) :: d(bas%nao, bas%nao, nd)
      type(trc_fock_options_t), intent(in) :: o
      real(dp), intent(out) :: g(bas%nao, bas%nao, nd)
      g = -99.0_dp
      eri%many_chunk = chunk
      !$acc enter data copyin(d) create(g)
      if (split) then
         call eri%fock_split_many_lent(bas, nd, d, g, o)
      else
         call eri%fock_nosym_many_lent(bas, nd, d, g, o)
      end if
      !$acc update self(g)
      !$acc exit data delete(d, g)
      eri%many_chunk = 0
   end subroutine lent_run

   subroutine check_same(error, a, b, label)
      !! Equal bit for bit on the serial host build, within atomic-order noise on the device and in the threaded OpenMP
      !! port (its atomics add in a thread-dependent order: ~1e-15 measured at 4 threads).
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: a(:, :, :), b(:, :, :)
      character(len=*), intent(in) :: label
      real(dp) :: rel
      rel = maxval(abs(a - b))/max(maxval(abs(a)), tiny(1.0_dp))
#if defined(_OPENACC) || defined(_OPENMP)
      call check(error, rel < TOL_LENT, label//": relative "//to_string(rel, "es10.3"))
#else
      call check(error, all(a == b), label//": not bit-identical, relative "//to_string(rel, "es10.3"))
#endif
   end subroutine check_same

   subroutine test_chunk_rule(error)
      type(error_type), allocatable, intent(out) :: error
      ! The two systems the rule was measured on (1000 and 1525 functions, 20 densities) take 3, the
      ! hexamer (150 functions) the whole batch; the cap at the batch and the floor of 1.
      call check(error, eri_lent_chunk_size(1000, 20) == 3, "40 waters: "//to_string(eri_lent_chunk_size(1000, 20)))
      if (allocated(error)) return
      call check(error, eri_lent_chunk_size(1525, 20) == 3, "gly20: "//to_string(eri_lent_chunk_size(1525, 20)))
      if (allocated(error)) return
      call check(error, eri_lent_chunk_size(150, 20) == 20, "hexamer: "//to_string(eri_lent_chunk_size(150, 20)))
      if (allocated(error)) return
      call check(error, eri_lent_chunk_size(1000, 2) == 2, "a batch smaller than the chunk is the batch")
      if (allocated(error)) return
      call check(error, eri_lent_chunk_size(1000, 0) == 1, "an empty batch still has a positive chunk")
   end subroutine test_chunk_rule

   subroutine mixed_batch(n, nd, d)
      !! nd (>= 7) densities: 1 symmetric plus an antisymmetric part at 1e-16 of it (below
      !! SPLIT_ASYM_TOL, so only a batch-wide decision to run the antisymmetric pass adds it), 2 and 3
      !! exactly symmetric, the others non-symmetric.
      integer, intent(in) :: n, nd
      real(dp), intent(out) :: d(n, n, nd)
      real(dp), allocatable :: ds(:, :, :), da(:, :, :)
      allocate (ds(n, n, nd), da(n, n, nd))
      call nonsym_densities(n, nd, d)
      call split_parts(n, nd, d, ds, da)
      d(:, :, 1) = ds(:, :, 1) + 1.0e-16_dp*da(:, :, 1)
      d(:, :, 2) = ds(:, :, 2)
      d(:, :, 3) = ds(:, :, 3)
   end subroutine mixed_batch

   subroutine test_chunked_split(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_chunked(bas, eri, .true., error)
      call eri%release(); call bas%release()
   end subroutine test_chunked_split

   subroutine test_chunked_nosym(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_chunked(bas, eri, .false., error)
      call eri%release(); call bas%release()
   end subroutine test_chunked_nosym

   subroutine run_chunked(bas, eri, split, error)
      !! Chunks of 2 (7 densities: 2, 2, 2, 1, the last shorter) and 3 against the whole batch
      !! and the host routine, for (j, k) in {(1,1), (1,0), (0,1)}. The whole batch is one launch
      !! per pass, the chunked builds more, and each density's result is the same.
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      logical, intent(in) :: split
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NDC = 7
      real(dp), allocatable :: d(:, :, :), gw(:, :, :), gc(:, :, :), gh(:, :, :)
      real(dp) :: sc(2, 3)
      type(trc_fock_options_t) :: o
      integer :: nao, imix, ic, chunks(2), nl_whole, expect
      character(len=64) :: label
      nao = bas%nao
      allocate (d(nao, nao, NDC), gw(nao, nao, NDC), gc(nao, nao, NDC), gh(nao, nao, NDC))
      call mixed_batch(nao, NDC, d)
      sc = reshape([1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
      chunks = [2, 3]
      do imix = 1, 3
         o = trc_fock_options_t(j_scale=sc(1, imix), k_scale=sc(2, imix), density_screen=.false.)
         if (split) then
            call eri%fock_split_many(bas, NDC, d, gh, o)
         else
            call eri%fock_nosym_many(bas, NDC, d, gh, o)
         end if
         call lent_run(eri, bas, split, NDC, d, o, 0, gw)
         nl_whole = eri%nlent_launch
         call check(error, maxval(abs(gw)) > 1.0e-3_dp, "reference is trivially zero")
         if (allocated(error)) return
         call check_same(error, gh, gw, "whole batch vs host, (j,k) "//to_string(sc(1, imix), "f4.1")// &
                         to_string(sc(2, imix), "f4.1"))
         if (allocated(error)) return
         do ic = 1, 2
            call lent_run(eri, bas, split, NDC, d, o, chunks(ic), gc)
            label = "chunk "//to_string(chunks(ic))//", (j,k) "//to_string(sc(1, imix), "f4.1")//to_string(sc(2, imix), "f4.1")
            call check_same(error, gw, gc, trim(label)//" vs the whole batch")
            if (allocated(error)) return
            call check_same(error, gh, gc, trim(label)//" vs the host")
            if (allocated(error)) return
            ! chunks x passes launches: the batch really was split.
            expect = ((NDC + chunks(ic) - 1)/chunks(ic))
            if (split) expect = expect*eri%nsplit_passes
            call check(error, eri%nlent_launch == expect, trim(label)//": "//to_string(eri%nlent_launch)// &
                       " launches, expected "//to_string(expect))
            if (allocated(error)) return
            call check(error, eri%nlent_launch > nl_whole, &
                       trim(label)//": not more launches than the whole batch ("//to_string(nl_whole)//")")
            if (allocated(error)) return
         end do
         if (split .and. sc(2, imix) /= 0.0_dp) then
            ! The antisymmetric pass ran although the chunk of densities 1 and 2 (round-off
            ! asymmetric, exactly symmetric) would not have asked for it.
            call check(error, eri%nsplit_passes == 2, "the batch needs the antisymmetric pass")
            if (allocated(error)) return
         end if
      end do
   end subroutine run_chunked

   subroutine region_weight(n, kind, w)
      !! Per-function weights: kind 1 all ones, 2 one on the first half and tiny on the second, 3 the reverse.
      integer, intent(in) :: n, kind
      real(dp), intent(out) :: w(n)
      integer :: i
      do i = 1, n
         w(i) = 1.0_dp
         if (kind == 2 .and. i > n/2) w(i) = 1.0e-12_dp
         if (kind == 3 .and. i <= n/2) w(i) = 1.0e-12_dp
      end do
   end subroutine region_weight

   subroutine threshold_batch(n, d)
      !! Six non-symmetric densities for the real-threshold screen: 1 dense (10 O(1)), 2 tiny (1e-11 of a dense one), 3 and 5 large
      !! on the first half of the functions with a 1e-12 floor elsewhere, 4 and 6 the reverse.
      integer, intent(in) :: n
      real(dp), intent(out) :: d(n, n, 6)
      real(dp) :: w(n), v(n)
      integer :: k, i, j
      call nonsym_densities(n, 6, d)
      ! 10 times a dense density: its symmetric and antisymmetric parts are both O(10), so (but for the blocks that vanish
      ! identically, the antisymmetric part of an s shell's own block) a quartet that survives Schwarz survives the screen.
      d(:, :, 1) = 10.0_dp*d(:, :, 1)
      d(:, :, 2) = 1.0e-11_dp*d(:, :, 2)
      do k = 3, 6
         if (mod(k, 2) == 1) then
            call region_weight(n, 2, w)
         else
            call region_weight(n, 3, w)
         end if
         v = w
         do j = 1, n
            do i = 1, n
               d(i, j, k) = d(i, j, k)*v(i)*v(j)
            end do
         end do
      end do
   end subroutine threshold_batch

   subroutine zero_block_batch(n, d)
      !! Five densities whose blocks are exactly zero outside a pattern: 1 only the block of rows in the first third and columns in
      !! the last third (its transpose block is zero), 2 the transposed pattern, 3 symmetric on both of those blocks, 4 only the
      !! first-third diagonal block, non-symmetric, 5 exactly antisymmetric on the (first, last) third blocks, so its
      !! symmetric part vanishes (a bound kept from the symmetric pass would drop every quartet it needs). The nonzero elements are 10 + O(1): every nonzero block then has a maximum
      !! of at least 9, so a quartet that reads one survives the density screen at any threshold whenever it survives
      !! Schwarz (Q dscr > Q), and the screen drops only the quartets that read zero blocks alone.
      integer, intent(in) :: n
      real(dp), intent(out) :: d(n, n, 5)
      real(dp), allocatable :: b(:, :, :)
      integer :: i, j
      logical :: ra, rb, ca, cb
      allocate (b(n, n, 5))
      call nonsym_densities(n, 5, b)
      d = 0.0_dp
      do j = 1, n
         do i = 1, n
            ra = i <= n/3
            rb = i > 2*n/3
            ca = j <= n/3
            cb = j > 2*n/3
            if (ra .and. cb) d(i, j, 1) = 10.0_dp + b(i, j, 1)
            if (rb .and. ca) d(i, j, 2) = 10.0_dp + b(i, j, 2)
            if ((ra .and. cb) .or. (rb .and. ca)) d(i, j, 3) = 10.0_dp + cos(real(3*min(i, j) + 5*max(i, j), dp))
            if (ra .and. ca) d(i, j, 4) = 10.0_dp + b(i, j, 4)
            ! Exactly antisymmetric across the (A, B) blocks: D(i,j) = c, D(j,i) = -c, so the symmetric part is
            ! identically zero and only the antisymmetric (exchange-only) pass has anything to digest.
            if (ra .and. cb) d(i, j, 5) = 10.0_dp + b(i, j, 5)
            if (rb .and. ca) d(i, j, 5) = -(10.0_dp + b(j, i, 5))
         end do
      end do
   end subroutine zero_block_batch

   subroutine test_screen_zero_blocks(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_screen_zero_blocks(bas, eri, error)
      call eri%release(); call bas%release()
   end subroutine test_screen_zero_blocks

   subroutine run_screen_zero_blocks(bas, eri, error)
      !! A quartet the screen drops here has every density block it reads exactly zero, so it
      !! contributes exactly nothing: the screened result must equal the unscreened one to the bit
      !! (host build). A bound that is too small in ANY orientation drops a quartet whose
      !! transposed block is the nonzero one and breaks the equality: density 1 has block (A,B)
      !! nonzero and (B,A) zero, density 2 the reverse, which the unsymmetrised bound of the
      !! enumerated kernel gets wrong (checked by removing `symmetrise_bound_device`).
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: d(:, :, :), gs(:, :, :), gu(:, :, :), gh(:, :, :)
      type(trc_fock_options_t) :: os, ou
      integer :: nao, ic, chunks(3), ir
      integer(int64) :: nk_on, nk_off
      character(len=48) :: label
      nao = bas%nao
      allocate (d(nao, nao, 5), gs(nao, nao, 5), gu(nao, nao, 5), gh(nao, nao, 5))
      call zero_block_batch(nao, d)
      os = trc_fock_options_t(j_scale=1.0_dp, k_scale=0.7_dp, density_screen=.true., count_survivors=.true.)
      ou = trc_fock_options_t(j_scale=1.0_dp, k_scale=0.7_dp, density_screen=.false., count_survivors=.true.)
      ! Chunks of 1 launch density 1 (block (A,B) nonzero, (B,A) zero) alone: the bound is then its own, the one that is
      ! asymmetric unless the enumerated kernel's route symmetrises it. In larger chunks the transposed density 2 does.
      chunks = [0, 2, 1]
      do ir = 1, 2
         do ic = 1, 3
            label = merge("split", "nosym", ir == 1)//", chunk "//to_string(chunks(ic))
            call lent_run(eri, bas, ir == 1, 5, d, os, chunks(ic), gs)
            nk_on = eri%nkept
            call lent_run(eri, bas, ir == 1, 5, d, ou, chunks(ic), gu)
            nk_off = eri%nkept
            call check(error, maxval(abs(gu)) > 1.0e-3_dp, "reference is trivially zero")
            if (allocated(error)) return
            call check_same(error, gu, gs, trim(label)//": screened vs unscreened")
            if (allocated(error)) return
            if (ir == 1) then
               call global_logger%info("  zero blocks "//trim(label)//": survivors "//to_string(nk_on)//" of "//to_string(nk_off))
               call check(error, nk_on < nk_off, trim(label)//": the screen dropped no quartet ("//to_string(nk_on)//" of "// &
                          to_string(nk_off)//")")
               if (allocated(error)) return
            end if
            if (chunks(ic) == 0) then
               if (ir == 1) then
                  call eri%fock_split_many(bas, 5, d, gh, os)
               else
                  call eri%fock_nosym_many(bas, 5, d, gh, os)
               end if
               call check_same(error, gh, gs, trim(label)//": screened device route vs screened host route")
               if (allocated(error)) return
            end if
         end do
      end do
   end subroutine run_screen_zero_blocks

   subroutine test_screen_threshold(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error, thr=SCREEN_THRESH)
      if (.not. built) return
      call run_screen_threshold(bas, eri, error)
      call eri%release(); call bas%release()
   end subroutine test_screen_threshold

   subroutine run_screen_threshold(bas, eri, error)
      !! The screen at an integral threshold of 1e-8, on densities with blocks of every size down to
      !! 1e-24 (`threshold_batch`).
      !!
      !! The bound on what screening may change. A quartet is dropped only when Q_ab Q_cd times the
      !! kernel's density factor is at most the threshold t (dscr = 4 max(J blocks) when j /= 0, the
      !! largest of the four K blocks when k /= 0), so every integral of it times every density
      !! element it multiplies is at most t (J: t/4) in magnitude. One dropped shell quartet adds to
      !! one element of G at most 8 permutations x (at most ncart^2 functions of the quartet's other
      !! pair) such terms, each weighted j or k/2, and the folded kernels accumulate at twice the
      !! scale: |dG| <= 2 x 8 x ncart^2 x (|j| + |k|) x t per dropped quartet. The dropped quartets
      !! are those the survivor counts lose (`nkept` unscreened minus screened, summed over the
      !! launches, an upper bound for any one density's own). It is an upper bound with every
      !! multiplicity at its maximum, so the observed difference is far below it; it still excludes
      !! a bound that dropped quartets of size one, which is the failure to catch.
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: JF = 1.0_dp, KF = 0.7_dp
      real(dp), allocatable :: d(:, :, :), gs(:, :, :), gu(:, :, :), gh(:, :, :), g1(:, :, :)
      type(trc_fock_options_t) :: os, ou, osk, ouk
      integer :: nao, ncart, lmax
      integer(int64) :: nk_on, nk_off, nk_pair_on, nk_pair_off, nk_small, nk_big
      real(dp) :: diff, bound, small_alone, shared
      nao = bas%nao
      lmax = maxval(bas%sh_l)
      ncart = (lmax + 1)*(lmax + 2)/2
      allocate (d(nao, nao, 6), gs(nao, nao, 6), gu(nao, nao, 6), gh(nao, nao, 6), g1(nao, nao, 1))
      call threshold_batch(nao, d)
      os = trc_fock_options_t(j_scale=JF, k_scale=KF, density_screen=.true., count_survivors=.true.)
      ou = trc_fock_options_t(j_scale=JF, k_scale=KF, density_screen=.false., count_survivors=.true.)

      ! (a) One chunk: the device bound is the host routine's bound, so the results agree exactly.
      call lent_run(eri, bas, .true., 6, d, os, 0, gs)
      call eri%fock_split_many(bas, 6, d, gh, os)
      call check_same(error, gh, gs, "split, one chunk: screened device vs screened host")
      if (allocated(error)) return
      call lent_run(eri, bas, .false., 6, d, os, 0, gs)
      call eri%fock_nosym_many(bas, 6, d, gh, os)
      call check_same(error, gh, gs, "nosym, one chunk: screened device vs screened host")
      if (allocated(error)) return

      ! (b) Chunks of 2: per-launch bounds, never larger than the batch's. Against the unscreened build.
      call lent_run(eri, bas, .true., 6, d, os, 2, gs)
      nk_on = eri%nkept
      call lent_run(eri, bas, .true., 6, d, ou, 2, gu)
      nk_off = eri%nkept
      diff = maxval(abs(gs - gu))
      bound = 2.0_dp*8.0_dp*real(ncart*ncart, dp)*(abs(JF) + abs(KF))*SCREEN_THRESH*real(nk_off - nk_on, dp)
      call global_logger%info("  threshold screen, split, chunks of 2: survivors "//to_string(nk_on)//" of "//to_string(nk_off)// &
                              ", max|screened - unscreened| "//to_string(diff, "es10.3")//", bound "//to_string(bound, "es10.3"))
      call check(error, nk_on < nk_off, "the screen dropped no quartet")
      if (allocated(error)) return
      call check(error, diff > 0.0_dp, "work was dropped but the result did not change")
      if (allocated(error)) return
      call check(error, diff <= bound, "screened vs unscreened "//to_string(diff, "es10.3")//" exceeds the bound "// &
                 to_string(bound, "es10.3"))
      if (allocated(error)) return
      ! The enumerated kernel has no survivor counter, and with the Coulomb term on the (first, first) blocks of densities 3 and 5
      ! keep every quartet that the exchange-only antisymmetric pass of the split drops, so its screen is exercised exchange only
      ! (j = 0: dscr is the four exchange blocks): the unscreened difference must be nonzero (the screen engaged) and the
      ! one-chunk result is the host routine's.
      osk = trc_fock_options_t(j_scale=0.0_dp, k_scale=KF, density_screen=.true.)
      ouk = trc_fock_options_t(j_scale=0.0_dp, k_scale=KF, density_screen=.false.)
      call lent_run(eri, bas, .false., 6, d, osk, 0, gs)
      call eri%fock_nosym_many(bas, 6, d, gh, osk)
      call check_same(error, gh, gs, "nosym, exchange only, one chunk: screened device vs screened host")
      if (allocated(error)) return
      call lent_run(eri, bas, .false., 6, d, osk, 2, gs)
      call lent_run(eri, bas, .false., 6, d, ouk, 2, gu)
      diff = maxval(abs(gs - gu))
      call global_logger%info("  threshold screen, nosym K only, chunks of 2: max|screened - unscreened| "// &
                              to_string(diff, "es10.3"))
      call check(error, diff > 0.0_dp, "nosym: the screen changed nothing")
      if (allocated(error)) return

      ! (c) The bound belongs to the launch. Densities 1 (dense, O(10)) and 2 (1e-11 of a dense one) share one: the bound is
      ! the dense one's (the launch keeps exactly what the dense density alone keeps, the few quartets the dense one's
      ! own zero blocks drop); alone, the tiny density's own bound drops almost every quartet and its result is mostly gone.
      ! What it loses in the shared launch is far less than that (compared, not bounded: its whole Fock matrix is below
      ! the integral threshold in absolute terms, so an absolute bound says nothing about it).
      call lent_run(eri, bas, .true., 2, d(:, :, 1:2), os, 0, gs(:, :, 1:2))
      nk_pair_on = eri%nkept
      call lent_run(eri, bas, .true., 2, d(:, :, 1:2), ou, 0, gu(:, :, 1:2))
      nk_pair_off = eri%nkept
      g1(:, :, 1) = 0.0_dp
      call lent_run(eri, bas, .true., 1, d(:, :, 1:1), os, 0, g1)
      nk_big = eri%nkept
      call check(error, nk_pair_on == nk_big, "the tiny density changed the launch's bound: "//to_string(nk_pair_on)// &
                 " survivors with it, "//to_string(nk_big)//" without")
      if (allocated(error)) return
      shared = maxval(abs(gs(:, :, 2) - gu(:, :, 2)))/maxval(abs(gu(:, :, 2)))
      call lent_run(eri, bas, .true., 1, d(:, :, 2:2), os, 0, g1)
      nk_small = eri%nkept
      small_alone = maxval(abs(g1(:, :, 1) - gu(:, :, 2)))/maxval(abs(gu(:, :, 2)))
      call global_logger%info("  tiny density: survivors "//to_string(nk_small)//" alone, "//to_string(nk_pair_on)// &
                              " with the dense one, "//to_string(nk_pair_off)//" unscreened; relative change of its "// &
                              "result alone "// &
                              to_string(small_alone, "es10.3")//", in the shared launch "//to_string(shared, "es10.3"))
      call check(error, shared < small_alone, "the tiny density lost as much in the shared launch as alone")
      if (allocated(error)) return
      call check(error, nk_small < nk_pair_on, "the tiny density's own launch dropped no more than the shared one")
      if (allocated(error)) return
      call check(error, small_alone > 1.0e-3_dp, "the tiny density alone should lose its result to the screen")
   end subroutine run_screen_threshold

   subroutine test_screen_off_chunks(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error, thr=SCREEN_THRESH)
      if (.not. built) return
      call run_screen_off_chunks(bas, eri, error)
      call eri%release(); call bas%release()
   end subroutine test_screen_off_chunks

   subroutine run_screen_off_chunks(bas, eri, error)
      !! density_screen = .false. is the plain Schwarz bound, whatever the densities and the chunking: the sparse
      !! batch gives the same unscreened result whole and in chunks (the survivors are the same too).
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: d(:, :, :), gw(:, :, :), gc(:, :, :), gh(:, :, :)
      type(trc_fock_options_t) :: ou
      integer :: nao, ir
      integer(int64) :: nk_whole
      nao = bas%nao
      allocate (d(nao, nao, 6), gw(nao, nao, 6), gc(nao, nao, 6), gh(nao, nao, 6))
      call threshold_batch(nao, d)
      ou = trc_fock_options_t(j_scale=1.0_dp, k_scale=0.7_dp, density_screen=.false., count_survivors=.true.)
      do ir = 1, 2
         call lent_run(eri, bas, ir == 1, 6, d, ou, 0, gw)
         nk_whole = eri%nkept
         call lent_run(eri, bas, ir == 1, 6, d, ou, 2, gc)
         call check_same(error, gw, gc, merge("split", "nosym", ir == 1)//": unscreened, chunks of 2 vs whole")
         if (allocated(error)) return
         if (ir == 1) then
            ! The chunked count sums three launches of the same Schwarz survivors.
            call check(error, eri%nkept == 3*nk_whole, "chunked survivors "//to_string(eri%nkept)//", whole "//to_string(nk_whole))
            if (allocated(error)) return
         end if
         if (ir == 1) then
            call eri%fock_split_many(bas, 6, d, gh, ou)
         else
            call eri%fock_nosym_many(bas, 6, d, gh, ou)
         end if
         call check_same(error, gh, gc, merge("split", "nosym", ir == 1)//": unscreened chunks vs host")
         if (allocated(error)) return
      end do
   end subroutine run_screen_off_chunks

   subroutine symmetrise_batch(n, nd, d, k0, k1)
      !! Densities k0..k1 made exactly symmetric.
      integer, intent(in) :: n, nd, k0, k1
      real(dp), intent(inout) :: d(n, n, nd)
      integer :: k
      do k = k0, k1
         d(:, :, k) = 0.5_dp*(d(:, :, k) + transpose(d(:, :, k)))
      end do
   end subroutine symmetrise_batch

   subroutine test_lent_folded(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_lent_folded(bas, eri, error)
      call eri%release(); call bas%release()
   end subroutine test_lent_folded

   subroutine lent_folded(eri, bas, nd, d, o, anti, chunk, g)
      !! One `fock_folded_many_lent` build, the result on the host.
      type(trc_eri_t), intent(inout) :: eri
      type(trc_basis_t), intent(in) :: bas
      integer, intent(in) :: nd, chunk
      real(dp), intent(in) :: d(bas%nao, bas%nao, nd)
      type(trc_fock_options_t), intent(in) :: o
      logical, intent(in) :: anti
      real(dp), intent(out) :: g(bas%nao, bas%nao, nd)
      g = -99.0_dp
      eri%many_chunk = chunk
      !$acc enter data copyin(d) create(g)
      call eri%fock_folded_many_lent(bas, nd, d, g, o, antisymmetric=anti)
      !$acc update self(g)
      !$acc exit data delete(d, g)
      eri%many_chunk = 0
   end subroutine lent_folded

   subroutine run_lent_folded(bas, eri, error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: d(:, :, :), ds(:, :, :), da(:, :, :), gref(:, :, :), gh(:, :, :), gl(:, :, :)
      type(trc_fock_options_t) :: o
      type(error_t) :: err
      real(dp) :: sc(2, 5), rel
      integer :: nao, imix, ia, ichunk, nl
      logical :: anti
      character(len=:), allocatable :: label
      nao = bas%nao
      allocate (d(nao, nao, ND), ds(nao, nao, ND), da(nao, nao, ND), gref(nao, nao, ND), gh(nao, nao, ND), gl(nao, nao, ND))
      call nonsym_densities(nao, ND, d)
      call split_parts(nao, ND, d, ds, da)
      ! The exactness the route trusts, checked on the inputs themselves.
      call check(error, all(ds == transpose_batch(nao, ND, ds)), "symmetric input is not exactly symmetric")
      if (allocated(error)) return
      call check(error, all(da == -transpose_batch(nao, ND, da)), "antisymmetric input is not exactly antisymmetric")
      if (allocated(error)) return
      ! (j, k): the symmetric batch takes all four, the antisymmetric one only j = 0.
      sc = reshape([1.0_dp, 1.0_dp, 4.0_dp, 1.5_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], [2, 5])
      do ia = 1, 2
         anti = ia == 2
         do imix = 1, 5
            if (anti .and. sc(1, imix) /= 0.0_dp) cycle
            o = trc_fock_options_t(j_scale=sc(1, imix), k_scale=sc(2, imix), density_screen=.false.)
            label = merge("antisymmetric", "symmetric    ", anti)//" (j,k) "//to_string(sc(1, imix), "f5.2")// &
                    to_string(sc(2, imix), "f5.2")
            if (anti) then
               call eri%fock_nosym_many(bas, ND, da, gref, o)
               call eri%fock_folded_many(bas, ND, da, gh, o, antisymmetric=.true.)
               call lent_folded(eri, bas, ND, da, o, .true., 0, gl)
            else
               call eri%fock_nosym_many(bas, ND, ds, gref, o)
               call eri%fock_folded_many(bas, ND, ds, gh, o)
               call lent_folded(eri, bas, ND, ds, o, .false., 0, gl)
            end if
            nl = eri%nlent_launch
            if (sc(2, imix) /= 0.0_dp .or. sc(1, imix) /= 0.0_dp) then
               call check(error, maxval(abs(gref)) > 1.0e-3_dp, "reference is trivially zero")
               if (allocated(error)) return
            end if
            rel = maxval(abs(gref - gl))/max(maxval(abs(gref)), tiny(1.0_dp))
            call global_logger%info("  lent folded "//label//": max|lent - nosym|/max|G| "//to_string(rel, "es10.3"))
            call check(error, rel < TOL_SPLIT, "lent folded vs nosym, "//label//", relative "//to_string(rel, "es10.3"))
            if (allocated(error)) return
            call check_same(error, gh, gl, "lent folded vs host folded, "//label)
            if (allocated(error)) return
            call check(error, eri%nsplit_passes == 1 .and. nl >= 1, "one folded pass, "//label)
            if (allocated(error)) return
            ! Chunks of two (2 + 1): each density's arithmetic is the same.
            do ichunk = 2, 2
               if (anti) then
                  call lent_folded(eri, bas, ND, da, o, .true., ichunk, gl)
               else
                  call lent_folded(eri, bas, ND, ds, o, .false., ichunk, gl)
               end if
               call check(error, eri%nlent_launch == 2, "chunks of 2 over 3 densities are two launches, "//label)
               if (allocated(error)) return
               call check_same(error, gh, gl, "lent folded in chunks vs host folded, "//label)
               if (allocated(error)) return
            end do
         end do
      end do
      ! Screening on: one chunk is the host's bound over the batch, exactly.
      o = trc_fock_options_t(j_scale=2.0_dp, k_scale=1.0_dp, density_screen=.true.)
      call eri%fock_folded_many(bas, ND, ds, gh, o)
      call lent_folded(eri, bas, ND, ds, o, .false., 0, gl)
      call check_same(error, gh, gl, "lent folded, screened, vs host folded")
      if (allocated(error)) return
      o = trc_fock_options_t(j_scale=0.0_dp, k_scale=1.0_dp, density_screen=.true.)
      call eri%fock_folded_many(bas, ND, da, gh, o, antisymmetric=.true.)
      call lent_folded(eri, bas, ND, da, o, .true., 0, gl)
      call check_same(error, gh, gl, "lent folded, antisymmetric, screened, vs host folded")
      if (allocated(error)) return
      ! Negative control: a density with no symmetry is NOT what this route
      ! computes, so the test can fail.
      o = trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, density_screen=.false.)
      call eri%fock_nosym_many(bas, ND, d, gref, o)
      call lent_folded(eri, bas, ND, d, o, .false., 0, gl)
      rel = maxval(abs(gref - gl))/maxval(abs(gref))
      call check(error, rel > TOL_CONTROL, "a non-symmetric density through the folded route must be wrong, relative "// &
                 to_string(rel, "es10.3"))
      if (allocated(error)) return
      ! The refusal: the Coulomb half of an antisymmetric batch.
      gl = 7.0_dp
      !$acc enter data copyin(da, gl)
      call eri%fock_folded_many_lent(bas, ND, da, gl, trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, &
                                                                         density_screen=.false.), antisymmetric=.true., err=err)
      !$acc update self(gl)
      !$acc exit data delete(da, gl)
      call check(error, err%has_error(), "antisymmetric with j_scale = 1 must be refused")
      if (allocated(error)) return
      call check(error, err%code == ERROR_VALIDATION, "the refusal is a validation error")
      if (allocated(error)) return
      call check(error, all(gl == 7.0_dp), "a refused call must leave gmats untouched")
   end subroutine run_lent_folded

   function transpose_batch(n, nd, a) result(t)
      !! Each (n, n) slice of `a` transposed.
      integer, intent(in) :: n, nd
      real(dp), intent(in) :: a(n, n, nd)
      real(dp) :: t(n, n, nd)
      integer :: k
      do k = 1, nd
         t(:, :, k) = transpose(a(:, :, k))
      end do
   end function transpose_batch

   subroutine test_lent_nosym(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      logical :: built
      call setup(bas, eri, built, error)
      if (.not. built) return
      call run_lent_nosym(bas, eri, error)
      call eri%release(); call bas%release()
   end subroutine test_lent_nosym

   subroutine run_lent_nosym(bas, eri, error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      type(trc_eri_t), intent(inout) :: eri
      real(dp), allocatable :: d(:, :, :), gref(:, :, :), glent(:, :, :)
      type(trc_fock_options_t) :: o
      real(dp) :: rel
      integer :: nao
      nao = bas%nao
      allocate (d(nao, nao, ND), gref(nao, nao, ND), glent(nao, nao, ND))
      call nonsym_densities(nao, ND, d)
      o = trc_fock_options_t(j_scale=1.0_dp, k_scale=0.7_dp, density_screen=.false.)
      call eri%fock_nosym_many(bas, ND, d, gref, o)
      glent = -99.0_dp
      !$acc enter data copyin(d) create(glent)
      call eri%fock_nosym_many_lent(bas, ND, d, glent, o)
      !$acc update self(glent)
      !$acc exit data delete(d, glent)
      rel = maxval(abs(gref - glent))/maxval(abs(gref))
      call global_logger%info("  lent nosym: max|lent - host|/max|G| "//to_string(rel, "es10.3"))
      call check(error, maxval(abs(gref)) > 1.0e-3_dp, "reference is trivially zero")
      if (allocated(error)) return
      call check_same(error, gref, glent, "lent nosym vs host nosym")
   end subroutine run_lent_nosym

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
