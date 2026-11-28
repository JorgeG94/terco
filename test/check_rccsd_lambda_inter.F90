module check_rccsd_lambda_inter_suite
   !!
   !! The once-per-solve intermediates of the RI-CCSD Lambda equations (src/post_hf/ccsd/trc_rccsd_lambda_inter.F90)
   !! on water cc-pVDZ with cc-pVDZ-RIFIT, frozen core 1 and 0, and on a C1 hydrogen peroxide, at RANDOM amplitudes: t1 without symmetry and a
   !! t2 with only t2(j,i,b,a) = t2(i,j,a,b), so that a swapped index shows (converged or MP2 amplitudes hide some).
   !!
   !! The reference is host loops: every (pq|rs) of the active space is formed from b_oo, b_ov and b_vv, the
   !! (ia|bc) included, and each intermediate is written term by term as the einsums of pyscf's
   !! ccsd_lambda.make_intermediates in a closed-shell, canonical form, in THEIR index order (p_* arrays below); the
   !! stored arrays are compared through the index map of the module header. No routine, packed block or helper of
   !! the module or of the CCSD stages is used. v1, v2 and v4, which the module takes from the CCSD one-particle
   !! intermediates, are compared too.
   !!
   !!   frozen 1, frozen 0   every intermediate against the host loops
   !!   blocks               the woovo (id|cb) term with 1 and 3 virtuals per block
   !!   least scratch        the contraction plans with extra_work_limit = 0 on H2O2 (frozen 0), the smallest of these
   !!                        where that changes plans: more of them must run as batches of GEMMs than by default
   !!   refusals             each asserted by its message, with nothing written
   !!
   !! With TERCO_LAMBDA_WOOVO_CHUNK=<k> in the environment (a second ctest entry sets it) the block taken must be
   !! capped at k in every case.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t, trc_rccsd_scratch_t
   use trc_rccsd_lambda_inter, only: trc_rccsd_lambda_inter_t, trc_rccsd_lambda_options_t, rccsd_lambda_inter_build
   use trc_device_memory, only: device_memory_set_cap
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: load_molecule, run_scf
   implicit none
   private
   public :: collect_check_rccsd_lambda_inter, release_fixture

   real(dp), parameter :: TOL_HOST = 1.0e-12_dp  !! max abs difference over max(1, |reference|max)
   integer, parameter :: NNAME = 10
   character(len=8), parameter :: NAMES(NNAME) = [character(len=8) :: "v1", "v2", "v4", "woooo", "w_both", &
                                                  "w_diff", "wovvo", "wovvo_x", "woovo", "w3"]

   type :: ref_t
      !! The reference of one frozen core, in the stored index orders (filled from the p_* arrays)
      integer :: no = 0, nv = 0
      real(dp), allocatable :: v1(:, :), v2(:, :), v4(:, :), w3(:, :)
      real(dp), allocatable :: woooo(:, :, :, :), w_both(:, :, :, :), w_diff(:, :, :, :)
      real(dp), allocatable :: wovvo(:, :, :, :), wovvo_x(:, :, :, :), woovo(:, :, :, :)
   end type ref_t

   integer, parameter :: NFIX = 3
   integer, parameter :: FIX_MOL(NFIX) = [1, 1, 2]     !! molecule of rccsd_test_support
   integer, parameter :: FIX_FROZEN(NFIX) = [1, 0, 0]  !! frozen core
   integer, parameter :: WATER_F1 = 1, WATER_F0 = 2, H2O2_F0 = 3

   type :: fixture_t
      !! One molecule and frozen core: reference orbitals, integrals, amplitudes and the host-loop reference
      logical :: built = .false.
      logical :: supported = .true.
      integer :: nocc = 0
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_amps_t) :: amps
      type(ref_t) :: ref
   end type fixture_t

   logical :: bound = .false.
   type(trc_linalg_t) :: la
   type(fixture_t) :: fix(NFIX)

contains

   subroutine collect_check_rccsd_lambda_inter(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("frozen 1: against host loops", test_frozen1), &
                   new_unittest("frozen 0: against host loops", test_frozen0), &
                   new_unittest("woovo in blocks", test_blocks), &
                   new_unittest("least contraction scratch", test_least_scratch), &
                   new_unittest("refusals", test_refusals)]
   end subroutine collect_check_rccsd_lambda_inter

   subroutine test_frozen1(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_lambda_options_t) :: opts
      call run_case(error, WATER_F1, opts)
   end subroutine test_frozen1

   subroutine test_frozen0(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_lambda_options_t) :: opts
      call run_case(error, WATER_F0, opts)
   end subroutine test_frozen0

   subroutine test_blocks(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_lambda_options_t) :: opts

      opts%woovo_batch = 1
      call run_case(error, WATER_F1, opts)
      if (allocated(error)) return
      ! 19 virtuals in blocks of 3: six whole ones and one of 1
      opts%woovo_batch = 3
      call run_case(error, WATER_F1, opts)
   end subroutine test_blocks

   subroutine test_least_scratch(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_lambda_options_t) :: opts
      integer :: items_default, items_least

      call run_case(error, H2O2_F0, opts, items_default)
      if (allocated(error) .or. .not. fix(H2O2_F0)%supported) return
      opts%extra_work_limit = 0_int64
      call run_case(error, H2O2_F0, opts, items_least)
      if (allocated(error)) return
      call check(error, items_least > items_default, "extra_work_limit = 0 changed no plan: "// &
                 to_string(items_least)//" GEMMs against "//to_string(items_default))
   end subroutine test_least_scratch

   subroutine run_case(error, which, opts, gemm_items)
      !! One build of one fixture with one set of options, every intermediate against the host loops
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: which
      type(trc_rccsd_lambda_options_t), intent(in) :: opts
      integer, intent(out), optional :: gemm_items

      type(trc_rccsd_inter_t) :: inter
      type(trc_rccsd_scratch_t) :: scr
      type(trc_rccsd_lambda_inter_t) :: im
      real(dp) :: diff(NNAME)
      integer :: k, nv, want_batch, cap

      if (present(gemm_items)) gemm_items = 0
      call ensure_fixture(which)
      if (.not. fix(which)%supported) return
      nv = fix(which)%ints%nv
      call inter%init(fix(which)%ints%no, nv)
      call scr%init(fix(which)%ints)
      call im%init(fix(which)%ints%no, nv)
      call check(error, .not. im%status%has_error(), "init: "//im%status%get_message())
      if (.not. allocated(error)) then
         call rccsd_lambda_inter_build(fix(which)%ints, fix(which)%amps, inter, scr, im, la, opts)
         call check(error, .not. im%status%has_error() .and. im%built, "build: "//im%status%get_message())
      end if
      if (.not. allocated(error)) then
         ! the whole range of virtuals when nothing forces less: these are far below the block budget
         want_batch = nv
         if (opts%woovo_batch > 0) want_batch = min(opts%woovo_batch, nv)
         cap = env_cap()
         if (cap > 0) want_batch = min(want_batch, cap)
         call check(error, im%woovo_batch == want_batch, "woovo block taken "//to_string(im%woovo_batch)// &
                    ", expected "//to_string(want_batch))
      end if
      if (.not. allocated(error)) then
         call im%pull()
         call inter%pull()
         associate (r => fix(which)%ref)
            diff(1) = rel(maxval(abs(inter%lvv - r%v1)), maxval(abs(r%v1)))
            diff(2) = rel(maxval(abs(inter%loo - r%v2)), maxval(abs(r%v2)))
            diff(3) = rel(maxval(abs(inter%fov - r%v4)), maxval(abs(r%v4)))
            diff(4) = rel(maxval(abs(im%woooo - r%woooo)), maxval(abs(r%woooo)))
            diff(5) = rel(maxval(abs(im%w_both - r%w_both)), maxval(abs(r%w_both)))
            diff(6) = rel(maxval(abs(im%w_diff - r%w_diff)), maxval(abs(r%w_diff)))
            diff(7) = rel(maxval(abs(im%wovvo - r%wovvo)), maxval(abs(r%wovvo)))
            diff(8) = rel(maxval(abs(im%wovvo_x - r%wovvo_x)), maxval(abs(r%wovvo_x)))
            diff(9) = rel(maxval(abs(im%woovo - r%woovo)), maxval(abs(r%woovo)))
            diff(10) = rel(maxval(abs(im%w3 - r%w3)), maxval(abs(r%w3)))
         end associate
         do k = 1, NNAME
            ! written so that a NaN fails
            if (.not. diff(k) < TOL_HOST) then
               call check(error, .false., trim(NAMES(k))//" vs host loops (case "//to_string(which)//"): "// &
                          to_string(diff(k), "es9.2")//report(diff))
               exit
            end if
         end do
         if (present(gemm_items)) gemm_items = im%gemm_items
      end if
      call im%release()
      call scr%release()
      call inter%release()
   end subroutine run_case

   function report(diff) result(text)
      !! The names that differ, for a failure message
      real(dp), intent(in) :: diff(NNAME)
      character(len=:), allocatable :: text
      integer :: k

      text = "; all that differ:"
      do k = 1, NNAME
         if (.not. diff(k) < TOL_HOST) text = text//" "//trim(NAMES(k))
      end do
   end function report

   pure function rel(difference, scale) result(r)
      real(dp), intent(in) :: difference, scale
      real(dp) :: r

      r = difference/max(1.0_dp, scale)
   end function rel

   function env_cap() result(cap)
      !! TERCO_LAMBDA_WOOVO_CHUNK, 0 when unset
      integer :: cap
      character(len=16) :: text
      integer :: stat

      cap = 0
      call get_environment_variable("TERCO_LAMBDA_WOOVO_CHUNK", text, status=stat)
      if (stat /= 0) return
      read (text, *, iostat=stat) cap
      if (stat /= 0) cap = 0
   end function env_cap

   subroutine test_refusals(error)
      !! Every refusal, by its message; `im` and `inter` keep the sentinel
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_inter_t) :: inter
      type(trc_rccsd_scratch_t) :: scr
      type(trc_rccsd_lambda_inter_t) :: im, im_other
      type(trc_rccsd_amps_t) :: no_amps
      type(trc_rccsd_ints_t) :: rebuilt
      type(trc_rccsd_lambda_options_t) :: opts
      integer :: no, nv

      call ensure_fixture(WATER_F1)
      if (.not. fix(WATER_F1)%supported) return
      no = fix(WATER_F1)%ints%no
      nv = fix(WATER_F1)%ints%nv
      associate (ints => fix(WATER_F1)%ints, amps => fix(WATER_F1)%amps)
      call inter%init(no, nv)
      call scr%init(ints)

      call rccsd_lambda_inter_build(ints, amps, inter, scr, im, la)
      call expect(error, im, "the result is not initialised")
      if (allocated(error)) return

      call im%init(no, nv)
      call im_other%init(no, nv - 1)
      call rccsd_lambda_inter_build(ints, amps, inter, scr, im_other, la)
      call expect(error, im_other, "(no, nv) differ")
      if (allocated(error)) return
      call rccsd_lambda_inter_build(ints, no_amps, inter, scr, im, la)
      call expect(error, im, "the amplitudes are not initialised")
      if (allocated(error)) return
      call scr%release()
      call rccsd_lambda_inter_build(ints, amps, inter, scr, im, la)
      call expect(error, im, "scratch are not initialised")
      if (allocated(error)) return
      call scr%init(ints)
      opts%woovo_batch = -1
      call rccsd_lambda_inter_build(ints, amps, inter, scr, im, la, opts)
      call expect(error, im, "a negative woovo_batch")
      if (allocated(error)) return

      ! too little free memory: for the result, then for the build's work arrays
      call device_memory_set_cap(1000_int64)
      call im_other%init(no, nv)
      call device_memory_set_cap(0_int64)
      call expect(error, im_other, "the intermediates (o 4, v 19) needs")
      if (allocated(error)) return
      call check(error, .not. allocated(im_other%woooo), "a refused init must allocate nothing")
      if (allocated(error)) return
      call device_memory_set_cap(1000_int64)
      opts%woovo_batch = 0
      call rccsd_lambda_inter_build(ints, amps, inter, scr, im, la, opts)
      call device_memory_set_cap(0_int64)
      call expect(error, im, "the build's work arrays (o 4, v 19, naux")
      if (allocated(error)) return
      call check(error, index(im%status%get_message(), "woovo block 1") > 0 .and. &
                 index(im%status%get_message(), "GB are free") > 0, &
                 "the memory refusal states the block and the sizes: "//im%status%get_message())
      if (allocated(error)) return

      ! nothing was written by any of the refusals above
      call im%pull()
      call inter%pull()
      call check(error, all(im%woooo == 0.0_dp) .and. all(im%woovo == 0.0_dp) .and. all(im%w3 == 0.0_dp) .and. &
                 all(inter%tau == 0.0_dp) .and. all(inter%lvv == 0.0_dp) .and. .not. im%built, &
                 "a refused build must write neither the result nor the one-particle intermediates")
      if (allocated(error)) return

      ! a scratch of integrals since rebuilt
      call build_ints(fix(WATER_F1), rebuilt, 1)
      call rccsd_lambda_inter_build(rebuilt, amps, inter, scr, im, la)
      call expect(error, im, "re-init the scratch after ints%build")
      call rebuilt%release()
      call im%release()
      call im_other%release()
      call scr%release()
      call inter%release()
      end associate
   end subroutine test_refusals

   subroutine expect(error, im, fragment)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_lambda_inter_t), intent(in) :: im
      character(len=*), intent(in) :: fragment

      call check(error, im%status%has_error() .and. index(im%status%get_message(), fragment) > 0, &
                 "expected a refusal containing '"//fragment//"', got '"//im%status%get_message()//"'")
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
      if (x%status%has_error()) error stop "check_rccsd_lambda_inter: build: "//x%status%get_message()
   end subroutine build_ints

   subroutine ensure_fixture(which)
      integer, intent(in) :: which

      integer, allocatable :: zint(:), seed(:)
      real(dp), allocatable :: x(:, :, :, :)
      integer :: nseed, dev, no, nv, i, j, a, b

      if (fix(which)%built) return
      fix(which)%built = .true.
      if (.not. bound) then
         bound = .true.
         dev = trc_bind_device(0)
      end if
      call load_molecule(FIX_MOL(which), fix(which)%bas, fix(which)%aux, fix(which)%pl, zint, fix(which)%nocc, &
                         fix(which)%supported)
      if (.not. fix(which)%supported) then
         print "(a)", "SKIP check_rccsd_lambda_inter: this build's density-fitting kernels stop below the RIFIT shells"
         return
      end if
      call run_scf(fix(which)%bas, fix(which)%nocc, fix(which)%res)
      call build_ints(fix(which), fix(which)%ints, FIX_FROZEN(which))
      no = fix(which)%ints%no
      nv = fix(which)%ints%nv
      call random_seed(size=nseed)
      allocate (seed(nseed))
      seed = 20261126 + which
      call random_seed(put=seed)
      associate (amps => fix(which)%amps)
         call amps%init(no, nv)
         call random_number(amps%t1)
         amps%t1 = 0.4_dp*(amps%t1 - 0.5_dp)
         allocate (x(no, no, nv, nv))
         call random_number(x)
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     amps%t2(i, j, a, b) = 0.2_dp*(0.5_dp*(x(i, j, a, b) + x(j, i, b, a)) - 0.5_dp)
                  end do
               end do
            end do
         end do
         call amps%push()
         call host_reference(fix(which)%ints, amps%t1, amps%t2, fix(which)%ref)
      end associate
   end subroutine ensure_fixture

   subroutine host_integrals(x, g)
      !! g(p,q,r,s) = (pq|rs) over the active orbitals, occupied first, from the three-index factors
      type(trc_rccsd_ints_t), intent(in) :: x
      real(dp), allocatable, intent(out) :: g(:, :, :, :)

      real(dp), allocatable :: bf(:, :), gm(:, :)
      integer :: no, nv, n, p, q, r, s, i, j, a, b

      no = x%no
      nv = x%nv
      n = no + nv
      allocate (bf(n*n, x%naux))
      do j = 1, no
         do i = 1, no
            bf(i + n*(j - 1), :) = x%b_oo(i + no*(j - 1), :)
         end do
      end do
      do a = 1, nv
         do i = 1, no
            bf(i + n*(no + a - 1), :) = x%b_ov(i + no*(a - 1), :)
            bf(no + a + n*(i - 1), :) = x%b_ov(i + no*(a - 1), :)
         end do
      end do
      do b = 1, nv
         do a = 1, nv
            bf(no + a + n*(no + b - 1), :) = x%b_vv(a + nv*(b - 1), :)
         end do
      end do
      gm = matmul(bf, transpose(bf))
      allocate (g(n, n, n, n))
      do s = 1, n
         do r = 1, n
            do q = 1, n
               do p = 1, n
                  g(p, q, r, s) = gm(p + n*(q - 1), r + n*(s - 1))
               end do
            end do
         end do
      end do
   end subroutine host_integrals

   subroutine host_reference(x, t1, t2, r)
      !! The intermediates by loops (module header)
      type(trc_rccsd_ints_t), intent(in) :: x
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      type(ref_t), intent(inout) :: r

      real(dp), allocatable :: g(:, :, :, :), tau(:, :, :, :)
      ! the arrays of the einsums, in their index order: [jbck] for the rings, [ijck] for woovo, [ck] for w3
      real(dp), allocatable :: p_v(:, :, :, :), p_u(:, :, :, :), p_both(:, :, :, :), p_diff(:, :, :, :)
      real(dp), allocatable :: p_ov(:, :, :, :), p_x(:, :, :, :), p_woovo(:, :, :, :), p_woooo(:, :, :, :)
      real(dp), allocatable :: p_w3(:, :), p_v1(:, :), p_v2(:, :), a_term(:, :, :, :)
      integer :: no, nv, i, j, k, l, a, b, c, d
      real(dp) :: s

      no = x%no
      nv = x%nv
      r%no = no
      r%nv = nv
      call host_integrals(x, g)
      allocate (tau(no, no, nv, nv))
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  tau(i, j, a, b) = t2(i, j, a, b) + t1(i, a)*t1(j, b)
               end do
            end do
         end do
      end do

      ! v1[b,a], v2[i,j] without the orbital energies; v4[j,b]
      allocate (r%v1(nv, nv), r%v2(no, no), r%v4(no, nv))
      do a = 1, nv
         do b = 1, nv
            s = 0.0_dp
            do j = 1, no
               do c = 1, nv
                  do k = 1, no
                     s = s - ovov1(j, a, k, c)*tau(j, k, b, c)
                  end do
                  s = s + (2.0_dp*eri_ovvv(j, c, b, a) - eri_ovvv(j, a, b, c))*t1(j, c)
               end do
            end do
            r%v1(b, a) = s
         end do
      end do
      do j = 1, no
         do i = 1, no
            s = 0.0_dp
            do k = 1, no
               do b = 1, nv
                  do c = 1, nv
                     s = s + ovov1(i, b, k, c)*tau(j, k, b, c)
                  end do
                  s = s + ovoo1(k, b, i, j)*t1(k, b)
               end do
            end do
            r%v2(i, j) = s
         end do
      end do
      do b = 1, nv
         do j = 1, no
            s = 0.0_dp
            do c = 1, nv
               do k = 1, no
                  s = s + ovov1(j, b, k, c)*t1(k, c)
               end do
            end do
            r%v4(j, b) = s
         end do
      end do

      ! woooo[i,k,j,l]
      allocate (p_woooo(no, no, no, no))
      do l = 1, no
         do j = 1, no
            do k = 1, no
               do i = 1, no
                  p_woooo(i, k, j, l) = g(i, k, j, l)
               end do
            end do
         end do
      end do
      do l = 1, no
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  do c = 1, nv
                     p_woooo(i, k, j, l) = p_woooo(i, k, j, l) + eri_ovoo(i, c, j, l)*t1(k, c)
                     p_woooo(i, l, j, k) = p_woooo(i, l, j, k) + eri_ovoo(j, c, i, l)*t1(k, c)
                     do d = 1, nv
                        p_woooo(i, k, j, l) = p_woooo(i, k, j, l) + eri_ovov(i, c, j, d)*tau(k, l, c, d)
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! V and U [j,b,c,k], then the two combinations and the two dressed ones
      allocate (p_v(no, nv, nv, no), p_u(no, nv, nv, no), p_ov(no, nv, nv, no), p_x(no, nv, nv, no))
      do k = 1, no
         do c = 1, nv
            do b = 1, nv
               do j = 1, no
                  p_v(j, b, c, k) = eri_ovvo(j, b, c, k)
                  p_u(j, b, c, k) = -g(j, k, no + c, no + b)
                  do d = 1, nv
                     do l = 1, no
                        p_v(j, b, c, k) = p_v(j, b, c, k) + ovov1(l, d, j, b)*t2(k, l, c, d) &
                                          - eri_ovov(l, d, j, b)*t2(k, l, d, c)
                        p_u(j, b, c, k) = p_u(j, b, c, k) + eri_ovov(j, d, l, b)*t2(k, l, d, c)
                     end do
                  end do
                  p_ov(j, b, c, k) = p_v(j, b, c, k)
                  p_x(j, b, c, k) = p_u(j, b, c, k)
                  do l = 1, no
                     do d = 1, nv
                        p_ov(j, b, c, k) = p_ov(j, b, c, k) - eri_ovov(j, b, l, d)*t1(k, d)*t1(l, c)
                        p_x(j, b, c, k) = p_x(j, b, c, k) + eri_ovov(j, d, l, b)*t1(k, d)*t1(l, c)
                     end do
                     p_ov(j, b, c, k) = p_ov(j, b, c, k) - eri_ovoo(j, b, l, k)*t1(l, c)
                     p_x(j, b, c, k) = p_x(j, b, c, k) + eri_ovoo(l, b, j, k)*t1(l, c)
                  end do
                  do d = 1, nv
                     p_ov(j, b, c, k) = p_ov(j, b, c, k) + eri_ovvv(j, b, c, d)*t1(k, d)
                     p_x(j, b, c, k) = p_x(j, b, c, k) - eri_ovvv(j, d, c, b)*t1(k, d)
                  end do
               end do
            end do
         end do
      end do
      p_both = 2.0_dp*p_v + p_u
      p_diff = p_v - p_u

      ! woovo[i,j,c,k]
      allocate (p_woovo(no, no, nv, no), a_term(no, no, nv, no))
      do k = 1, no
         do c = 1, nv
            do j = 1, no
               do i = 1, no
                  s = 0.0_dp
                  do b = 1, nv
                     s = s + p_both(i, b, c, k)*t1(j, b)
                  end do
                  a_term(i, j, c, k) = s
               end do
            end do
         end do
      end do
      do k = 1, no
         do c = 1, nv
            do j = 1, no
               do i = 1, no
                  p_woovo(i, j, c, k) = a_term(i, j, c, k) - a_term(i, k, c, j) + ovoo1(k, c, j, i)
               end do
            end do
         end do
      end do
      do k = 1, no
         do c = 1, nv
            do j = 1, no
               do i = 1, no
                  do b = 1, nv
                     ! "ibck,jb->ikcj" of w_diff
                     p_woovo(i, k, c, j) = p_woovo(i, k, c, j) + p_diff(i, b, c, k)*t1(j, b)
                  end do
               end do
            end do
         end do
      end do
      do j = 1, no
         do b = 1, nv
            do k = 1, no
               do i = 1, no
                  do l = 1, no
                     do c = 1, nv
                        ! "lcik,jlbc->ikbj" (theta), "lcik,jlbc->ijbk" (t2), "iclk,ljbc->ijbk" (t2)
                        p_woovo(i, k, b, j) = p_woovo(i, k, b, j) + &
                                              ovoo1(l, c, i, k)*(2.0_dp*t2(j, l, b, c) - t2(j, l, c, b))
                        p_woovo(i, j, b, k) = p_woovo(i, j, b, k) - ovoo1(l, c, i, k)*t2(j, l, b, c) &
                                              - ovoo1(i, c, l, k)*t2(l, j, b, c)
                     end do
                  end do
               end do
            end do
         end do
      end do
      do k = 1, no
         do c = 1, nv
            do j = 1, no
               do i = 1, no
                  do b = 1, nv
                     do d = 1, nv
                        p_woovo(i, j, c, k) = p_woovo(i, j, c, k) + 2.0_dp*eri_ovvv(i, d, c, b)*tau(j, k, d, b) &
                                              - eri_ovvv(i, b, c, d)*tau(j, k, d, b)
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! w3[c,k], with the whole v1 and v2 (the orbital energies on their diagonals)
      p_v1 = r%v1
      do a = 1, nv
         p_v1(a, a) = p_v1(a, a) + x%eps_v(a)
      end do
      p_v2 = r%v2
      do i = 1, no
         p_v2(i, i) = p_v2(i, i) + x%eps_o(i)
      end do
      allocate (p_w3(nv, no))
      do k = 1, no
         do c = 1, nv
            s = 0.0_dp
            do b = 1, nv
               do j = 1, no
                  s = s + p_both(j, b, c, k)*t1(j, b)
               end do
            end do
            p_w3(c, k) = s
         end do
      end do
      do j = 1, no
         do b = 1, nv
            s = 0.0_dp
            do k = 1, no
               do c = 1, nv
                  s = s + r%v4(k, c)*t1(k, b)*t1(j, c)
                  do l = 1, no
                     s = s - ovoo1(l, c, k, j)*t2(k, l, b, c)
                  end do
                  do d = 1, nv
                     s = s + (2.0_dp*eri_ovvv(k, d, b, c) - eri_ovvv(k, c, b, d))*t2(j, k, c, d)
                  end do
               end do
            end do
            do c = 1, nv
               s = s + p_v1(b, c)*t1(j, c)
            end do
            do k = 1, no
               s = s - p_v2(k, j)*t1(k, b)
            end do
            p_w3(b, j) = p_w3(b, j) + s
         end do
      end do

      ! into the stored index orders (the module header of trc_rccsd_lambda_inter)
      allocate (r%woooo(no, no, no, no), r%woovo(no, no, no, nv), r%w3(no, nv))
      allocate (r%w_both(no, nv, no, nv), r%w_diff(no, nv, no, nv), r%wovvo(no, nv, no, nv), r%wovvo_x(no, nv, no, nv))
      do l = 1, no
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  r%woooo(i, j, k, l) = p_woooo(i, k, j, l)
               end do
            end do
         end do
      end do
      do c = 1, nv
         do k = 1, no
            do b = 1, nv
               do j = 1, no
                  r%w_both(j, b, k, c) = p_both(j, b, c, k)
                  r%w_diff(j, b, k, c) = p_diff(j, b, c, k)
                  r%wovvo(j, b, k, c) = 2.0_dp*p_ov(j, b, c, k) + p_x(j, b, c, k)
                  r%wovvo_x(j, b, k, c) = p_x(j, b, c, k)
               end do
            end do
         end do
      end do
      do c = 1, nv
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  r%woovo(i, j, k, c) = p_woovo(i, j, c, k)
               end do
            end do
         end do
      end do
      do b = 1, nv
         do j = 1, no
            r%w3(j, b) = p_w3(b, j)
         end do
      end do

   contains

      ! (pq|rs) by the spaces of its indices
      pure function eri_ovov(i, a, j, b) result(v)
         integer, intent(in) :: i, a, j, b
         real(dp) :: v
         v = g(i, x%no + a, j, x%no + b)
      end function eri_ovov

      pure function eri_ovvo(i, a, b, j) result(v)
         integer, intent(in) :: i, a, b, j
         real(dp) :: v
         v = g(i, x%no + a, x%no + b, j)
      end function eri_ovvo

      pure function eri_ovoo(i, a, j, k) result(v)
         integer, intent(in) :: i, a, j, k
         real(dp) :: v
         v = g(i, x%no + a, j, k)
      end function eri_ovoo

      pure function eri_ovvv(i, a, b, c) result(v)
         integer, intent(in) :: i, a, b, c
         real(dp) :: v
         v = g(i, x%no + a, x%no + b, x%no + c)
      end function eri_ovvv

      pure function ovov1(i, a, j, b) result(v)
         !! 2 (ia|jb) - (ib|ja)
         integer, intent(in) :: i, a, j, b
         real(dp) :: v
         v = 2.0_dp*eri_ovov(i, a, j, b) - eri_ovov(i, b, j, a)
      end function ovov1

      pure function ovoo1(k, a, i, j) result(v)
         !! 2 (ka|ij) - (ia|kj)
         integer, intent(in) :: k, a, i, j
         real(dp) :: v
         v = 2.0_dp*eri_ovoo(k, a, i, j) - eri_ovoo(i, a, k, j)
      end function ovoo1

   end subroutine host_reference

   subroutine release_fixture()
      integer :: n

      ! every release is safe on what was never built
      do n = 1, NFIX
         call fix(n)%amps%release()
         call fix(n)%ints%release()
         call fix(n)%pl%release()
         call fix(n)%aux%release()
         call fix(n)%bas%release()
      end do
      call la%release()
   end subroutine release_fixture

end module check_rccsd_lambda_inter_suite

program check_rccsd_lambda_inter
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_lambda_inter_suite, only: collect_check_rccsd_lambda_inter, release_fixture
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_lambda_inter, error_unit, stat, parallel=.false.)
   call release_fixture()
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_lambda_inter
