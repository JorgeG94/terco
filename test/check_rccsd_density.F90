module check_rccsd_density_suite
   !!
   !! The unrelaxed RI-CCSD densities (src/post_hf/ccsd/trc_rccsd_density.F90: gamma and g3 = Gamma b), three ways.
   !!
   !!   synthetic     random factors, amplitudes and multipliers (no SCF, nothing converged): every element of gamma
   !!                 and g3 against the dense two-particle density of rccsd_density_host (pyscf's blocks, the v^4
   !!                 one included, in host loops)
   !!   references    water distorted to C1, STO-3G and cc-pVDZ, frozen core 0 and 1, at terco's converged t and l:
   !!                 GAMMA_SQ, GAMMA_TR_OCC, GAMMA_TR_VIR, G3_NORM and E_FROM_DENSITIES of pyscf's densities
   !!                 (test/pyscf_refs/ccsd_grad_water_*.ref), the energy of the densities against terco's own CCSD
   !!                 correlation energy, and the dense reference at the same amplitudes, whose ||Gamma|| is
   !!                 GAMMA2_NORM
   !!   slabs         one to four virtuals per slab and the least contraction scratch: the same densities; the bytes
   !!                 `rccsd_density_bytes` reports, a build that just fits them and one a byte short
   !!   refusals      by message, with nothing left allocated
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_linalg, only: trc_linalg_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_driver, only: trc_rccsd_options_t, trc_rccsd_state_t, trc_rccsd_solve
   use trc_rccsd_lambda, only: trc_rccsd_lambda_solve_options_t, trc_rccsd_lambda_state_t, trc_rccsd_lambda_solve
   use trc_rccsd_density, only: trc_rccsd_density_t, trc_rccsd_density_options_t, trc_rccsd_density_solve, &
                                rccsd_density_build, rccsd_density_bytes
   use rccsd_density_host, only: dense_densities
   use rccsd_test_support, only: run_scf
   use ccsd_grad_support, only: NCASE, CASE_REF, CASE_BASIS, CASE_NFROZEN, ccsd_grad_ref_t, read_ccsd_grad_ref, &
                                geometry, supported
   implicit none
   private
   public :: collect_check_rccsd_density, release_fixture

   ! Positions of the density probes in `ccsd_grad_ref_t%probes`
   integer, parameter :: I_GAMMA_SQ = 5, I_TR_OCC = 6, I_TR_VIR = 7, I_G3 = 8, I_GAMMA2 = 9, I_ENERGY = 10
   ! TOLERANCES, each at least ten times the worst difference measured (host builds with flang and nvfortran, V100)
   real(dp), parameter :: TOL_DENSE = 1.0e-12_dp
   !! max |gamma|, |g3| element difference from the dense reference; worst measured 9.8e-15
   real(dp), parameter :: TOL_PROBE = 3.0e-11_dp
   !! the pyscf probes, absolute; worst measured 2.9e-12 (GAMMA2_NORM, STO-3G)
   real(dp), parameter :: TOL_ENERGY = 1.0e-12_dp
   !! energy of the densities against terco's CCSD; worst measured 4.2e-14
   real(dp), parameter :: TOL_BATCH = 1.0e-12_dp
   !! between slab sizes and scratch budgets; worst measured 2.8e-14

   type :: fixture_t
      logical :: built = .false.
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_rccsd_state_t) :: cc
      type(trc_rccsd_lambda_state_t) :: lam
      type(ccsd_grad_ref_t) :: ref
   end type fixture_t

   type :: synthetic_t
      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_amps_t) :: t, l
   end type synthetic_t

   logical :: bound = .false.
   type(trc_linalg_t) :: la
   type(fixture_t) :: fix(NCASE)
   integer(int64) :: seed = 1_int64

contains

   subroutine collect_check_rccsd_density(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("random amplitudes: every element against the dense densities", test_synthetic), &
                   new_unittest("STO-3G, frozen core 0 and 1: the references", test_sto3g), &
                   new_unittest("cc-pVDZ, frozen core 0 and 1: the references", test_ccpvdz), &
                   new_unittest("slab sizes and scratch budgets give the same densities", test_batches), &
                   new_unittest("refusals", test_refusals)]
   end subroutine collect_check_rccsd_density

   subroutine setup()
      integer :: dev

      if (bound) return
      bound = .true.
      dev = trc_bind_device(0)
      call la%init(2)
   end subroutine setup

   function uniform() result(x)
      !! A reproducible number in (-1/2, 1/2): the same inputs from every compiler
      real(dp) :: x

      seed = mod(48271_int64*seed, 2147483647_int64)
      x = real(seed, dp)/2147483647.0_dp - 0.5_dp
   end function uniform

   subroutine make_synthetic(no, nv, naux, s)
      !! b symmetric in its orbital pair, t2 and l2 pair-symmetric, nothing converged
      integer, intent(in) :: no, nv, naux
      type(synthetic_t), intent(inout) :: s

      real(dp), allocatable :: b(:, :, :)
      integer :: n, p, q, x

      n = no + nv
      allocate (b(n, n, naux))
      do x = 1, naux
         do q = 1, n
            do p = 1, q
               b(p, q, x) = uniform()
               b(q, p, x) = b(p, q, x)
            end do
         end do
      end do
      s%ints%no = no
      s%ints%nv = nv
      s%ints%naux = naux
      allocate (s%ints%b_oo(no*no, naux), s%ints%b_ov(no*nv, naux), s%ints%b_vv(nv*nv, naux))
      s%ints%b_oo = reshape(b(1:no, 1:no, :), [no*no, naux])
      s%ints%b_ov = reshape(b(1:no, no + 1:n, :), [no*nv, naux])
      s%ints%b_vv = reshape(b(no + 1:n, no + 1:n, :), [nv*nv, naux])
      !$acc enter data copyin(s%ints%b_oo, s%ints%b_ov, s%ints%b_vv)
      call random_amps(no, nv, s%t)
      call random_amps(no, nv, s%l)
   end subroutine make_synthetic

   subroutine random_amps(no, nv, a)
      integer, intent(in) :: no, nv
      type(trc_rccsd_amps_t), intent(inout) :: a

      integer :: i, j, b, c
      real(dp) :: x

      call a%init(no, nv)
      do b = 1, nv
         do i = 1, no
            a%t1(i, b) = uniform()
         end do
      end do
      do c = 1, nv
         do b = 1, nv
            do j = 1, no
               do i = 1, no
                  ! each unordered pair of (i,b), (j,c) once
                  if (i + no*(b - 1) > j + no*(c - 1)) cycle
                  x = uniform()
                  a%t2(i, j, b, c) = x
                  a%t2(j, i, c, b) = x
               end do
            end do
         end do
      end do
      call a%push()
   end subroutine random_amps

   subroutine release_synthetic(s)
      type(synthetic_t), intent(inout) :: s

      call s%ints%release()
      call s%t%release()
      call s%l%release()
   end subroutine release_synthetic

   subroutine full_g3(den, g3)
      !! g3(p,q,P) over the active orbitals from the library's three blocks
      type(trc_rccsd_density_t), intent(in) :: den
      real(dp), allocatable, intent(out) :: g3(:, :, :)

      integer :: no, nv, x, p

      no = den%no
      nv = den%nv
      allocate (g3(no + nv, no + nv, den%naux))
      do p = 1, den%naux
         g3(1:no, 1:no, p) = reshape(den%g3_oo(:, p), [no, no])
         g3(1:no, no + 1:, p) = reshape(den%g3_ov(:, p), [no, nv])
         g3(no + 1:, no + 1:, p) = reshape(den%g3_vv(:, p), [nv, nv])
         do x = 1, nv
            g3(no + x, 1:no, p) = g3(1:no, no + x, p)
         end do
      end do
   end subroutine full_g3

   subroutine compare_dense(ints, t, l, den, d_gamma, d_g3, gamma2_norm)
      !! Max element differences of gamma and g3 from the dense reference at the same amplitudes (host copies), and
      !! ||Gamma|| of that reference
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: t, l
      type(trc_rccsd_density_t), intent(in) :: den
      real(dp), intent(out) :: d_gamma, d_g3, gamma2_norm

      real(dp), allocatable :: gam(:, :), big(:, :, :, :), b(:, :, :), want(:, :, :), got(:, :, :)
      integer :: no, nv, n, x, p, q

      no = ints%no
      nv = ints%nv
      n = no + nv
      call dense_densities(no, nv, t%t1, t%t2, l%t1, l%t2, gam, big)
      gamma2_norm = norm2(big)
      allocate (b(n, n, ints%naux), want(n, n, ints%naux))
      do p = 1, ints%naux
         b(1:no, 1:no, p) = reshape(ints%b_oo(:, p), [no, no])
         b(1:no, no + 1:, p) = reshape(ints%b_ov(:, p), [no, nv])
         b(no + 1:, no + 1:, p) = reshape(ints%b_vv(:, p), [nv, nv])
         do x = 1, nv
            b(no + x, 1:no, p) = b(1:no, no + x, p)
         end do
      end do
      want = reshape(matmul(reshape(big, [n*n, n*n]), reshape(b, [n*n, ints%naux])), [n, n, ints%naux])
      call full_g3(den, got)
      d_g3 = maxval(abs(got - want))
      d_gamma = max(maxval(abs(den%gamma_oo - gam(1:no, 1:no))), maxval(abs(den%gamma_ov - gam(1:no, no + 1:))), &
                    maxval(abs(den%gamma_vv - gam(no + 1:, no + 1:))))
      do q = 1, no
         do p = 1, nv
            d_gamma = max(d_gamma, abs(den%gamma_ov(q, p) - gam(no + p, q)))
         end do
      end do
   end subroutine compare_dense

   subroutine test_synthetic(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NSIZE = 4
      integer, parameter :: SIZES(3, NSIZE) = reshape([1, 1, 2, 3, 5, 11, 2, 9, 7, 4, 3, 5], [3, NSIZE])
      type(synthetic_t) :: s
      type(trc_rccsd_density_t) :: den
      type(trc_rccsd_density_options_t) :: opts
      real(dp) :: d_gamma, d_g3, norm
      integer :: k
      character(len=48) :: what

      call setup()
      do k = 1, NSIZE
         call make_synthetic(SIZES(1, k), SIZES(2, k), SIZES(3, k), s)
         call rccsd_density_build(s%ints, s%t, s%l, opts, la, den)
         what = "no "//to_string(SIZES(1, k))//" nv "//to_string(SIZES(2, k))//" naux "//to_string(SIZES(3, k))
         call check(error, den%built, trim(what)//": "//den%status%get_message())
         if (allocated(error)) return
         call den%pull()
         call compare_dense(s%ints, s%t, s%l, den, d_gamma, d_g3, norm)
         print "(a,a,a,2es10.2,a,es10.2)", "  ", trim(what), ": max |difference| gamma, g3 ", d_gamma, d_g3, &
            "   max |g3| ", maxval(abs(den%g3_vv))
         call check(error, d_gamma < TOL_DENSE, trim(what)//": gamma differs from the dense one by "// &
                    to_string(d_gamma, "es9.2"))
         if (allocated(error)) return
         call check(error, d_g3 < TOL_DENSE, trim(what)//": g3 differs from the dense one by "//to_string(d_g3, "es9.2"))
         if (allocated(error)) return
         call release_synthetic(s)
         call den%release()
      end do
   end subroutine test_synthetic

   subroutine test_sto3g(error)
      type(error_type), allocatable, intent(out) :: error
      call reference_case(error, 1)
      if (allocated(error)) return
      call reference_case(error, 2)
   end subroutine test_sto3g

   subroutine test_ccpvdz(error)
      type(error_type), allocatable, intent(out) :: error
      call reference_case(error, 3)
      if (allocated(error)) return
      call reference_case(error, 4)
   end subroutine test_ccpvdz

   subroutine reference_case(error, ic)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ic

      integer, parameter :: NP = 5
      integer, parameter :: AT(NP) = [I_GAMMA_SQ, I_TR_OCC, I_TR_VIR, I_G3, I_ENERGY]
      character(len=16), parameter :: PROBE(NP) = [character(len=16) :: "GAMMA_SQ", "GAMMA_TR_OCC", "GAMMA_TR_VIR", &
         "G3_NORM", "E_FROM_DENSITIES"]
      type(trc_rccsd_density_t) :: den
      type(trc_rccsd_density_options_t) :: opts
      real(dp) :: got(NP), diff(NP), e_one, e_two, tr(2), d_gamma, d_g3, gamma2_norm, d_e, d_2
      integer :: k, i

      if (.not. supported()) return
      call ensure_fixture(ic, error)
      if (allocated(error)) return
      call trc_rccsd_density_solve(fix(ic)%cc, fix(ic)%lam, opts, la, den)
      call check(error, den%built, trim(CASE_REF(ic))//": "//den%status%get_message())
      if (allocated(error)) return
      call den%pull()
      associate (ints => fix(ic)%cc%ints)
         tr = 0.0_dp
         e_one = 0.0_dp
         do i = 1, ints%no
            tr(1) = tr(1) + den%gamma_oo(i, i)
            e_one = e_one + den%gamma_oo(i, i)*ints%eps_o(i)
         end do
         do i = 1, ints%nv
            tr(2) = tr(2) + den%gamma_vv(i, i)
            e_one = e_one + den%gamma_vv(i, i)*ints%eps_v(i)
         end do
         e_two = 0.5_dp*(sum(den%g3_oo*ints%b_oo) + 2.0_dp*sum(den%g3_ov*ints%b_ov) + sum(den%g3_vv*ints%b_vv))
      end associate
      got(1) = sum(den%gamma_oo**2) + 2.0_dp*sum(den%gamma_ov**2) + sum(den%gamma_vv**2)
      got(2) = tr(1)
      got(3) = tr(2)
      got(4) = sqrt(sum(den%g3_oo**2) + 2.0_dp*sum(den%g3_ov**2) + sum(den%g3_vv**2))
      got(5) = e_one + e_two
      diff = abs(got - fix(ic)%ref%probes(AT, 1))
      d_e = abs(got(5) - fix(ic)%cc%result%e_ccsd)
      call compare_dense(fix(ic)%cc%ints, fix(ic)%cc%amps, fix(ic)%lam%l, den, d_gamma, d_g3, gamma2_norm)
      d_2 = abs(gamma2_norm - fix(ic)%ref%probes(I_GAMMA2, 1))
      print "(a,a,5es10.2)", trim(CASE_REF(ic)), ": |terco - pyscf| of GAMMA_SQ TR_OCC TR_VIR G3_NORM E_FROM_DENS ", diff
      print "(a,es10.2,a,2es10.2,a,es10.2)", "    |E from densities - terco CCSD| ", d_e, &
         "   against dense: gamma, g3 ", d_gamma, d_g3, "   |GAMMA2_NORM| of the dense one ", d_2
      do k = 1, NP
         ! written so that a NaN fails
         call check(error, diff(k) < TOL_PROBE, trim(CASE_REF(ic))//" "//trim(PROBE(k))//": "// &
                    to_string(got(k), "es22.14")//" against "//to_string(fix(ic)%ref%probes(AT(k), 1), "es22.14"))
         if (allocated(error)) return
      end do
      call check(error, d_e < TOL_ENERGY, trim(CASE_REF(ic))//": the densities give "//to_string(got(5), "es22.14")// &
                 " for the CCSD energy "//to_string(fix(ic)%cc%result%e_ccsd, "es22.14"))
      if (allocated(error)) return
      call check(error, d_gamma < TOL_DENSE .and. d_g3 < TOL_DENSE, trim(CASE_REF(ic))// &
                 ": differs from the dense densities: "//to_string(d_gamma, "es9.2")//" "//to_string(d_g3, "es9.2"))
      if (allocated(error)) return
      call check(error, d_2 < TOL_PROBE, trim(CASE_REF(ic))//": GAMMA2_NORM of the dense reference "// &
                 to_string(gamma2_norm, "es22.14"))
      call den%release()
   end subroutine reference_case

   subroutine test_batches(error)
      type(error_type), allocatable, intent(out) :: error

      type(synthetic_t) :: s
      type(trc_rccsd_density_t) :: den, other
      type(trc_rccsd_density_options_t) :: opts
      real(dp) :: d
      integer :: k, want

      call setup()
      ! nv 17: a slab holds at most 4 virtuals, so 1, 2, 3 and 4 give 17, 9, 6 and 5 slabs
      call make_synthetic(3, 17, 13, s)
      opts%slab_batch = 4
      call rccsd_density_build(s%ints, s%t, s%l, opts, la, den)
      call check(error, den%built, "four virtuals a slab: "//den%status%get_message())
      if (allocated(error)) return
      call den%pull()
      do k = 1, 4
         opts = trc_rccsd_density_options_t()
         if (k <= 3) opts%slab_batch = k
         if (k == 4) opts%extra_work_limit = 0_int64
         if (k == 4) opts%slab_batch = 17
         want = merge(k, 4, k <= 3)
         call rccsd_density_build(s%ints, s%t, s%l, opts, la, other)
         call check(error, other%built, "slab "//to_string(k)//": "//other%status%get_message())
         if (allocated(error)) return
         ! TERCO_CCDENS_D_CHUNK may cap it further
         call check(error, other%slab_batch <= want .and. other%slab_batch >= 1, &
                    "slab "//to_string(k)//": took "//to_string(other%slab_batch)//" virtuals a slab")
         if (allocated(error)) return
         call other%pull()
         d = max(maxval(abs(other%g3_oo - den%g3_oo)), maxval(abs(other%g3_ov - den%g3_ov)), &
                 maxval(abs(other%g3_vv - den%g3_vv)), maxval(abs(other%gamma_ov - den%gamma_ov)))
         print "(a,i0,a,i0,a,es10.2)", "  case ", k, ": slab ", other%slab_batch, ", max |difference| ", d
         call check(error, d < TOL_BATCH, "slab "//to_string(k)//" gives other densities: "//to_string(d, "es9.2"))
         if (allocated(error)) return
      end do
      call other%release()
      call memory_budget(error, s)
      call release_synthetic(s)
      call den%release()
   end subroutine test_batches

   subroutine memory_budget(error, s)
      !! With exactly the bytes of the result and one slab's temporaries free the build takes one virtual a slab;
      !! with a byte less it is refused
      type(error_type), allocatable, intent(out) :: error
      type(synthetic_t), intent(in) :: s

      type(trc_rccsd_density_t) :: den
      type(trc_rccsd_density_options_t) :: opts
      type(error_t) :: err
      integer(int64) :: held, temporary

      call rccsd_density_bytes(s%ints%no, s%ints%nv, s%ints%naux, 1, held, temporary, err)
      call check(error, .not. err%has_error() .and. held > 0 .and. temporary > 0, "rccsd_density_bytes")
      if (allocated(error)) return
      opts%device_free_cap = held + temporary
      call rccsd_density_build(s%ints, s%t, s%l, opts, la, den)
      call check(error, den%built, "with just enough memory: "//den%status%get_message())
      if (allocated(error)) return
      call check(error, den%slab_batch == 1 .and. den%temp_bytes == temporary, &
                 "with just enough memory: slab "//to_string(den%slab_batch)//", "//to_string(den%temp_bytes)// &
                 " temporary bytes for "//to_string(temporary))
      if (allocated(error)) return
      opts%device_free_cap = held + temporary - 1_int64
      call rccsd_density_build(s%ints, s%t, s%l, opts, la, den)
      call check(error, .not. den%built .and. .not. allocated(den%g3_vv) .and. &
                 index(den%status%get_message(), "GB are free") > 0, &
                 "a byte short was not refused: "//den%status%get_message())
      call den%release()
   end subroutine memory_budget

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_density_options_t) :: opts
      type(trc_rccsd_state_t) :: empty_cc
      type(trc_rccsd_lambda_state_t) :: empty_lam

      if (.not. supported()) return
      call ensure_fixture(2, error)
      if (allocated(error)) return
      associate (cc => fix(2)%cc, lam => fix(2)%lam)
         call refused(error, empty_cc, lam, opts, "the CCSD was refused or did not converge")
         if (allocated(error)) return
         empty_cc%result%converged = .true.
         call refused(error, empty_cc, lam, opts, "holds no integrals or amplitudes")
         if (allocated(error)) return
         allocate (cc%rank_seconds(4, 2))
         call refused(error, cc, lam, opts, "several ranks")
         deallocate (cc%rank_seconds)
         if (allocated(error)) return
         call refused(error, cc, empty_lam, opts, "Lambda equations were refused or did not converge")
         if (allocated(error)) return
         empty_lam%result%converged = .true.
         call refused(error, cc, empty_lam, opts, "holds no multipliers")
         if (allocated(error)) return
         opts%slab_batch = -1
         call refused(error, cc, lam, opts, "a negative slab_batch")
         if (allocated(error)) return
         opts = trc_rccsd_density_options_t()
         opts%device_free_cap = 1000_int64
         call refused(error, cc, lam, opts, "GB are free")
         if (allocated(error)) return
         ! and the same states still give densities
         opts = trc_rccsd_density_options_t()
         call solves(error, cc, lam, opts)
      end associate
   end subroutine test_refusals

   subroutine solves(error, cc, lam, opts)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_state_t), intent(in) :: cc
      type(trc_rccsd_lambda_state_t), intent(in) :: lam
      type(trc_rccsd_density_options_t), intent(in) :: opts
      type(trc_rccsd_density_t) :: den

      call trc_rccsd_density_solve(cc, lam, opts, la, den)
      call check(error, den%built, "after the refusals: "//den%status%get_message())
      call den%release()
   end subroutine solves

   subroutine refused(error, cc, lam, opts, fragment)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_state_t), intent(in) :: cc
      type(trc_rccsd_lambda_state_t), intent(in) :: lam
      type(trc_rccsd_density_options_t), intent(in) :: opts
      character(len=*), intent(in) :: fragment
      type(trc_rccsd_density_t) :: den

      call trc_rccsd_density_solve(cc, lam, opts, la, den)
      call check(error, .not. den%built .and. index(den%status%get_message(), fragment) > 0, &
                 "expected a refusal containing '"//fragment//"', got '"//den%status%get_message()//"'")
      if (allocated(error)) return
      call check(error, .not. allocated(den%gamma_oo) .and. .not. allocated(den%g3_vv), &
                 "a refused build left densities allocated")
   end subroutine refused

   subroutine ensure_fixture(ic, error)
      !! The SCF, the RI-CCSD and the Lambda solution of case `ic`, converged as the reference's
      integer, intent(in) :: ic
      type(error_type), allocatable, intent(out) :: error

      type(trc_scf_result_t) :: res
      type(trc_rccsd_options_t) :: opts
      type(trc_rccsd_lambda_solve_options_t) :: lopts
      type(error_t) :: err
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm, nocc

      call setup()
      if (.not. fix(ic)%built) then
         fix(ic)%built = .true.
         call read_ccsd_grad_ref(trim(CASE_REF(ic)), fix(ic)%ref)
         call geometry(natm, zint, at_r)
         nocc = sum(zint)/2
         associate (f => fix(ic))
            call trc_basis_from_json("basis_sets/"//trim(CASE_BASIS(ic))//".json", natm, zint, at_r, f%bas, err)
            if (err%has_error()) error stop "check_rccsd_density: "//err%get_message()
            ! the auxiliary basis follows the orbital basis's convention, as through the C API
            call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, f%aux, err, &
                                     force_cart=.not. f%bas%spherical)
            if (err%has_error()) error stop "check_rccsd_density: "//err%get_message()
            call f%bas%to_device()
            call f%aux%to_device()
            call f%pl%build(f%bas, 1.0e-12_dp)
            call f%pl%to_device()
            call run_scf(f%bas, nocc, res)
            opts%nocc = nocc
            opts%nfrozen = CASE_NFROZEN(ic)
            opts%conv_energy = 1.0e-13_dp
            opts%conv_amplitude = 1.0e-11_dp
            call trc_rccsd_solve(f%bas, f%aux, f%pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, f%cc)
            lopts%conv_amplitude = 1.0e-11_dp
            lopts%conv_norm = 1.0e-10_dp
            call trc_rccsd_lambda_solve(f%cc, lopts, la, f%lam)
         end associate
      end if
      call check(error, fix(ic)%cc%result%converged .and. fix(ic)%lam%result%converged, &
                 trim(CASE_REF(ic))//": CCSD or Lambda: "//trim(fix(ic)%cc%result%message)//" "// &
                 trim(fix(ic)%lam%result%message))
      if (allocated(error)) return
      call check(error, abs(fix(ic)%cc%result%e_ccsd - fix(ic)%ref%e_ccsd_corr) < 1.0e-10_dp .and. &
                 fix(ic)%cc%result%no == fix(ic)%ref%no .and. fix(ic)%cc%result%naux == fix(ic)%ref%naux, &
                 trim(CASE_REF(ic))//": not the reference's CCSD")
   end subroutine ensure_fixture

   subroutine release_fixture()
      integer :: n

      do n = 1, NCASE
         call fix(n)%lam%release()
         call fix(n)%cc%release()
         call fix(n)%pl%release()
         call fix(n)%aux%release()
         call fix(n)%bas%release()
      end do
      call la%release()
   end subroutine release_fixture

end module check_rccsd_density_suite

program check_rccsd_density
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_density_suite, only: collect_check_rccsd_density, release_fixture
   use ccsd_grad_support, only: supported
   implicit none
   integer :: stat
   stat = 0
   if (.not. supported()) then
      print "(a)", "SKIP check_rccsd_density: this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
      stop 0
   end if
   call run_testsuite(collect_check_rccsd_density, error_unit, stat, parallel=.false.)
   call release_fixture()
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_density
