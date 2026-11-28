module check_rccsd_ladder_tensor_suite
   !!
   !! `rccsd_ladder_tensor` (src/post_hf/ccsd/trc_rccsd_ladder.F90), the ladder of an arbitrary o^2 v^2 tensor, on water
   !! cc-pVDZ with cc-pVDZ-RIFIT. Two tensors: a random t0 with no symmetry, and t = (t0(i,j,a,b) + t0(j,i,b,a))/2,
   !! which has the symmetry of the Lambda amplitudes; t1 is random too, so nothing of a CCSD iteration (tau, MP2
   !! amplitudes) is involved. The reference is host loops over b_vv and b_ov, no routine or helper of the ladder:
   !!    L(i,j,a,b) = sum_cd (ac|bd) t(i,j,c,d)
   !!    D(i,j,a,b) = - sum_cd sum_k [ (kd|ac) t1(k,b) + (kc|bd) t1(k,a) ] t(i,j,c,d)
   !! The general paths (assembled, direct, the half form's two contractions) run on both: a gather or transpose
   !! that takes t(j,i,d,c) for t(i,j,c,d) is invisible on t alone. The half form's one-contraction path runs on t,
   !! and on t0, where it must give the ladder of t; its GEMM tally is compared with the general path's, since on t
   !! the two give the same numbers. Each undressed against L and dressed against L + D, with the default batches and
   !! with small ones (which start GEMMs and gathers inside t).
   !!
   !! `rccsd_ladder_transposed`, the transposed dressing, on the same two tensors and the same paths against
   !!    R(i,j,a,b) = sum_cd sum_P bT(a,c;P) bT(b,d;P) t(i,j,c,d),   bT(a,c;P) = b_vv(a,c;P) - sum_k b_ov(k,a;P) t1(k,c)
   !! and the identity the Lambda equations use it by (module header of the ladder, TRANSPOSED DRESSING), whose other
   !! side is the (kb|ac) slab term written in host loops, at t and a random pair-symmetric t2.
   !!
   !! With TERCO_RCCSD_LADDER_D_CHUNK=<k> in the environment (a second ctest entry sets it) every half-form plan of
   !! the transposed cases must be capped at k virtuals d a pass.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_ladder, only: trc_rccsd_ladder_t, trc_rccsd_lscratch_t, trc_rccsd_ladder_options_t, rccsd_ladder_tensor, &
                               RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_DIRECT, RCCSD_LADDER_HALF, rccsd_ladder_transposed
   use trc_device_memory, only: device_memory_set_cap
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: MOL_NFROZEN, load_molecule, run_scf
   implicit none
   private
   public :: collect_check_rccsd_ladder_tensor, release_fixture

   real(dp), parameter :: TOL_HOST = 1.0e-12_dp  !! max abs difference over max(1, |reference|max)
   !! Largest measured over every case here: 1.7e-15 on the gfortran 15 host build, 1.3e-15 on the nvfortran 26.9 host
   !! build, 7e-16 on the ifx 2025.0 host build (one thread each) and 7.2e-16 on a V100.

   logical :: built = .false.
   logical :: supported = .true.
   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl
   type(trc_rccsd_ints_t) :: ints
   type(trc_rccsd_amps_t) :: amps
   type(trc_linalg_t) :: la
   real(dp), allocatable :: t(:, :, :, :)                          !! the pair-symmetric tensor, mapped to the device
   real(dp), allocatable :: t0(:, :, :, :)                         !! the tensor without symmetry, mapped to the device
   real(dp), allocatable :: ref_l(:, :, :, :), ref_w(:, :, :, :)   !! L and L + D of t
   real(dp), allocatable :: ref0_l(:, :, :, :), ref0_w(:, :, :, :)  !! L and L + D of t0
   real(dp), allocatable :: ref_t(:, :, :, :), ref0_t(:, :, :, :)   !! the transposed ladder of t and of t0

contains

   subroutine collect_check_rccsd_ladder_tensor(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("assembled", test_assembled), &
                   new_unittest("direct", test_direct), &
                   new_unittest("half, one contraction", test_half_symmetric), &
                   new_unittest("half, two contractions", test_half_general), &
                   new_unittest("refusals", test_refusals), &
                   new_unittest("transposed: assembled and direct", test_transposed_full), &
                   new_unittest("transposed: half", test_transposed_half), &
                   new_unittest("transposed: the Lambda identity", test_identity), &
                   new_unittest("transposed: refusals", test_transposed_refusals)]
   end subroutine collect_check_rccsd_ladder_tensor

   subroutine test_assembled(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t) :: opts
      call ensure_fixture()
      opts%form = RCCSD_LADDER_ASSEMBLED
      call run_general(error, opts)
      if (allocated(error)) return
      opts%cd_batch = 7
      call run_general(error, opts)
   end subroutine test_assembled

   subroutine test_direct(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t) :: opts
      call ensure_fixture()
      opts%form = RCCSD_LADDER_DIRECT
      call run_general(error, opts)
      if (allocated(error)) return
      opts%p_batch = 5
      call run_general(error, opts)
   end subroutine test_direct

   subroutine test_half_general(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t) :: opts
      call ensure_fixture()
      opts%form = RCCSD_LADDER_HALF
      call run_general(error, opts)
      if (allocated(error)) return
      opts%d_batch = 2
      call run_general(error, opts)
   end subroutine test_half_general

   subroutine test_half_symmetric(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t) :: opts
      call ensure_fixture()
      opts%form = RCCSD_LADDER_HALF
      call run_symmetric(error, opts)
      if (allocated(error)) return
      opts%d_batch = 2
      call run_symmetric(error, opts)
   end subroutine test_half_symmetric

   subroutine run_general(error, opts)
      !! One plan without `tau_symmetric`, on the symmetric tensor and on the one without symmetry
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t), intent(in) :: opts
      integer(int64) :: flops(2)

      call run_form(error, opts, t, .false., ref_l, ref_w, flops)
      if (allocated(error)) return
      call run_form(error, opts, t0, .false., ref0_l, ref0_w, flops)
   end subroutine run_general

   subroutine run_symmetric(error, opts)
      !! One half-form plan with `tau_symmetric`: right on t, and on t0 the ladder of t = its symmetrised part. It must
      !! call one tau GEMM per pass where the general path calls two: over all passes, whatever the batches, the
      !! difference is 2 o^2 v^2 flops on each of the v (v + 1) / 2 columns c >= d, dressed or not.
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t), intent(in) :: opts

      integer(int64) :: one(2), two(2), second, no, nv

      call run_form(error, opts, t, .true., ref_l, ref_w, one)
      if (allocated(error)) return
      call run_form(error, opts, t, .false., ref_l, ref_w, two)
      if (allocated(error) .or. .not. supported) return
      no = ints%no
      nv = ints%nv
      second = no*no*nv*nv*nv*(nv + 1)
      call check(error, all(two - one == second), "one contraction did not run: general minus symmetric GEMM flops, "// &
                 "undressed "//to_string(two(1) - one(1))//", dressed "//to_string(two(2) - one(2))// &
                 ", expected "//to_string(second))
      if (allocated(error)) return
      call run_form(error, opts, t0, .true., ref_l, ref_w, one)
   end subroutine run_symmetric

   subroutine run_form(error, opts, x, sym, want_l, want_w, flops)
      !! The tensor `x`, undressed against `want_l`, dressed against `want_w`, for one plan
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t), intent(in) :: opts
      real(dp), intent(in), contiguous :: x(:, :, :, :)
      logical, intent(in) :: sym
      real(dp), intent(in) :: want_l(:, :, :, :), want_w(:, :, :, :)
      integer(int64), intent(out) :: flops(2)  !! `ls%flops_pp` of the undressed and of the dressed call

      type(trc_rccsd_lscratch_t) :: ls
      type(trc_rccsd_ladder_t) :: lad
      real(dp) :: d_l, d_w

      flops = 0_int64
      if (.not. supported) return
      call ls%init(ints, opts)
      call check(error, .not. ls%status%has_error(), "scratch: "//ls%status%get_message())
      if (allocated(error)) return
      call lad%init(ints%no, ints%nv)

      call rccsd_ladder_tensor(ints, x, lad, ls, la, tau_symmetric=sym)
      call check(error, .not. ls%status%has_error(), "undressed: "//ls%status%get_message())
      flops(1) = ls%flops_pp
      if (.not. allocated(error)) then
         call lad%pull()
         d_l = maxval(abs(lad%h - want_l))/max(1.0_dp, maxval(abs(want_l)))
         call check(error, d_l < TOL_HOST, "undressed vs host loops: "//to_string(d_l, "es9.2"))
      end if

      if (.not. allocated(error)) then
         call rccsd_ladder_tensor(ints, x, lad, ls, la, amps=amps, tau_symmetric=sym)
         call check(error, .not. ls%status%has_error(), "dressed: "//ls%status%get_message())
         flops(2) = ls%flops_pp
      end if
      if (.not. allocated(error)) then
         call lad%pull()
         d_w = maxval(abs(lad%h - want_w))/max(1.0_dp, maxval(abs(want_w)))
         call check(error, d_w < TOL_HOST, "dressed vs host loops: "//to_string(d_w, "es9.2"))
      end if
      call lad%release()
      call ls%release()
   end subroutine run_form

   subroutine test_refusals(error)
      !! A tensor of another shape, a scratch cut to a share and `tau_symmetric` outside the half form are refused and
      !! write nothing
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_ladder_options_t) :: opts
      type(trc_rccsd_lscratch_t) :: ls
      type(trc_rccsd_ladder_t) :: lad
      logical :: shape_refused, share_refused

      call ensure_fixture()
      if (.not. supported) return
      opts%form = RCCSD_LADDER_HALF
      call ls%init(ints, opts)
      call check(error, .not. ls%status%has_error(), "scratch: "//ls%status%get_message())
      if (allocated(error)) return
      call lad%init(ints%no, ints%nv)
      call rccsd_ladder_tensor(ints, t(:, :, :, 1:ints%nv - 1), lad, ls, la)
      shape_refused = index(ls%status%get_message(), "does not have the shape") > 0
      ls%d_lo = 2
      call rccsd_ladder_tensor(ints, t, lad, ls, la)
      share_refused = index(ls%status%get_message(), "takes no share") > 0
      call lad%pull()
      call check(error, shape_refused .and. share_refused .and. all(lad%h == 0.0_dp), &
                 "a wrong shape and a share must be refused, each for its own reason, with nothing written")
      if (.not. allocated(error)) call refuse_symmetric(error, RCCSD_LADDER_ASSEMBLED, lad)
      if (.not. allocated(error)) call refuse_symmetric(error, RCCSD_LADDER_DIRECT, lad)
      call lad%release()
      call ls%release()
   end subroutine test_refusals

   subroutine refuse_symmetric(error, form, lad)
      !! `tau_symmetric` with a form that has no one-contraction path: refused for that reason, undressed and
      !! dressed, and `lad` is not written
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: form
      type(trc_rccsd_ladder_t), intent(inout) :: lad

      type(trc_rccsd_ladder_options_t) :: opts
      type(trc_rccsd_lscratch_t) :: ls
      logical :: refused, refused_dressed

      opts%form = form
      call ls%init(ints, opts)
      call check(error, .not. ls%status%has_error(), "scratch: "//ls%status%get_message())
      if (allocated(error)) return
      lad%h = 0.0_dp
      !$acc update device(lad%h)
      call rccsd_ladder_tensor(ints, t, lad, ls, la, tau_symmetric=.true.)
      refused = index(ls%status%get_message(), "tau_symmetric") > 0
      call rccsd_ladder_tensor(ints, t, lad, ls, la, amps=amps, tau_symmetric=.true.)
      refused_dressed = index(ls%status%get_message(), "tau_symmetric") > 0
      call lad%pull()
      call check(error, refused .and. refused_dressed .and. all(lad%h == 0.0_dp), &
                 "tau_symmetric outside the half form must be refused with nothing written, form "//to_string(form))
      call ls%release()
   end subroutine refuse_symmetric

   subroutine test_transposed_full(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t) :: opts

      call ensure_fixture()
      opts%form = RCCSD_LADDER_ASSEMBLED
      call run_transposed_general(error, opts)
      if (allocated(error)) return
      opts%cd_batch = 7
      call run_transposed_general(error, opts)
      if (allocated(error)) return
      opts = trc_rccsd_ladder_options_t(form=RCCSD_LADDER_DIRECT)
      call run_transposed_general(error, opts)
      if (allocated(error)) return
      opts%p_batch = 5
      call run_transposed_general(error, opts)
   end subroutine test_transposed_full

   subroutine run_transposed_general(error, opts)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t), intent(in) :: opts
      integer(int64) :: flops

      call run_transposed(error, opts, t, .false., ref_t, flops)
      if (allocated(error)) return
      call run_transposed(error, opts, t0, .false., ref0_t, flops)
   end subroutine run_transposed_general

   subroutine test_transposed_half(error)
      !! Two contractions on both tensors; one contraction on t, and on t0 where it must give the ladder of t, with
      !! the GEMM tally showing which ran
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_ladder_options_t) :: opts
      integer(int64) :: one, two, second, no, nv
      integer :: d_batch

      call ensure_fixture()
      do d_batch = 0, 2, 2
         opts = trc_rccsd_ladder_options_t(form=RCCSD_LADDER_HALF, d_batch=d_batch)
         call run_transposed_general(error, opts)
         if (allocated(error)) return
         call run_transposed(error, opts, t, .true., ref_t, one)
         if (allocated(error)) return
         call run_transposed(error, opts, t, .false., ref_t, two)
         if (allocated(error) .or. .not. supported) return
         no = ints%no
         nv = ints%nv
         second = no*no*nv*nv*nv*(nv + 1)
         call check(error, two - one == second, "transposed: one contraction did not run, general minus symmetric "// &
                    "GEMM flops "//to_string(two - one)//", expected "//to_string(second))
         if (allocated(error)) return
         call run_transposed(error, opts, t0, .true., ref_t, one)
         if (allocated(error)) return
      end do
   end subroutine test_transposed_half

   subroutine run_transposed(error, opts, x, sym, want, flops)
      !! `rccsd_ladder_transposed` of `x` against `want`, twice with one scratch (the assembled and direct plans
      !! allocate their dressed factor at the first)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t), intent(in) :: opts
      real(dp), intent(in), contiguous :: x(:, :, :, :)
      logical, intent(in) :: sym
      real(dp), intent(in) :: want(:, :, :, :)
      integer(int64), intent(out) :: flops  !! `ls%flops_pp`

      type(trc_rccsd_lscratch_t) :: ls
      type(trc_rccsd_ladder_t) :: lad
      integer(int64) :: bytes
      real(dp) :: d
      integer :: pass

      flops = 0_int64
      if (.not. supported) return
      call ls%init(ints, opts)
      call check(error, .not. ls%status%has_error(), "scratch: "//ls%status%get_message())
      if (allocated(error)) return
      if (ls%form == RCCSD_LADDER_HALF .and. env_d_chunk() > 0) then
         call check(error, ls%d_batch <= env_d_chunk(), "d_batch "//to_string(ls%d_batch)//" above the environment's "// &
                    "cap "//to_string(env_d_chunk()))
         if (allocated(error)) return
      end if
      call lad%init(ints%no, ints%nv)
      do pass = 1, 2
         bytes = ls%bytes
         call rccsd_ladder_transposed(ints, x, amps, lad, ls, la, tau_symmetric=sym)
         call check(error, .not. ls%status%has_error(), "transposed: "//ls%status%get_message())
         if (allocated(error)) exit
         flops = ls%flops_pp
         call lad%pull()
         d = maxval(abs(lad%h - want))/max(1.0_dp, maxval(abs(want)))
         call check(error, d < TOL_HOST, "transposed vs host loops, form "//to_string(ls%form)//": "// &
                    to_string(d, "es9.2"))
         if (allocated(error)) exit
         if (pass == 2) call check(error, ls%bytes == bytes, "the dressed factor was allocated twice")
         if (allocated(error)) exit
      end do
      call lad%release()
      call ls%release()
   end subroutine run_transposed

   function env_d_chunk() result(cap)
      !! TERCO_RCCSD_LADDER_D_CHUNK, 0 when unset
      integer :: cap
      character(len=16) :: text
      integer :: stat

      cap = 0
      call get_environment_variable("TERCO_RCCSD_LADDER_D_CHUNK", text, status=stat)
      if (stat /= 0) return
      read (text, *, iostat=stat) cap
      if (stat /= 0) cap = 0
   end function env_d_chunk

   subroutine test_identity(error)
      !! With l2 = t and a random pair-symmetric t2 (tau = t2 + t1 t1), after (i,j,a,b) + (j,i,b,a):
      !!   1/2 L(l2) + 1/2 sum_kl (ka|lb) [l2 tau](ij;kl) - sum_kc (kb|ac) sum_d l2(ij;cd) t1(k,d)
      !!     = 1/2 R(l2) + 1/2 sum_kl (ka|lb) [l2 t2](ij;kl)
      !! left side by host loops (L = sum_cd (ac|bd) l2), R from the library
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_lscratch_t) :: ls
      type(trc_rccsd_ladder_t) :: lad
      type(trc_rccsd_ladder_options_t) :: opts
      real(dp), allocatable :: x(:, :, :, :), t2(:, :, :, :), old(:, :, :, :), new(:, :, :, :), ovov(:, :), ovvv(:, :)
      real(dp), allocatable :: l2t1(:, :, :, :), l2tau(:, :, :, :), l2t2(:, :, :, :)
      integer :: no, nv, i, j, a, b, c, d, k, l
      real(dp) :: s, diff, scale, raw

      call ensure_fixture()
      if (.not. supported) return
      no = ints%no
      nv = ints%nv
      allocate (x(no, no, nv, nv), t2(no, no, nv, nv), old(no, no, nv, nv), new(no, no, nv, nv))
      call random_number(x)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  t2(i, j, a, b) = 0.2_dp*(0.5_dp*(x(i, j, a, b) + x(j, i, b, a)) - 0.5_dp)
               end do
            end do
         end do
      end do
      ovov = matmul(ints%b_ov, transpose(ints%b_ov))
      ovvv = matmul(ints%b_ov, transpose(ints%b_vv))
      allocate (l2t1(no, no, nv, no), l2tau(no, no, no, no), l2t2(no, no, no, no))
      do k = 1, no
         do c = 1, nv
            do j = 1, no
               do i = 1, no
                  l2t1(i, j, c, k) = sum(t(i, j, c, :)*amps%t1(k, :))
               end do
            end do
         end do
      end do
      do l = 1, no
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  l2t2(i, j, k, l) = sum(t(i, j, :, :)*t2(k, l, :, :))
                  s = l2t2(i, j, k, l)
                  do d = 1, nv
                     do c = 1, nv
                        s = s + t(i, j, c, d)*amps%t1(k, c)*amps%t1(l, d)
                     end do
                  end do
                  l2tau(i, j, k, l) = s
               end do
            end do
         end do
      end do

      opts%form = RCCSD_LADDER_HALF
      call ls%init(ints, opts)
      call lad%init(no, nv)
      call rccsd_ladder_transposed(ints, t, amps, lad, ls, la, tau_symmetric=.true.)
      call check(error, .not. ls%status%has_error(), "transposed: "//ls%status%get_message())
      if (allocated(error)) return
      call lad%pull()

      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  old(i, j, a, b) = 0.5_dp*ref_l(i, j, a, b)
                  new(i, j, a, b) = 0.5_dp*lad%h(i, j, a, b)
                  do l = 1, no
                     do k = 1, no
                        s = 0.5_dp*ovov(k + no*(a - 1), l + no*(b - 1))
                        old(i, j, a, b) = old(i, j, a, b) + s*l2tau(i, j, k, l)
                        new(i, j, a, b) = new(i, j, a, b) + s*l2t2(i, j, k, l)
                     end do
                  end do
                  do c = 1, nv
                     do k = 1, no
                        old(i, j, a, b) = old(i, j, a, b) - ovvv(k + no*(b - 1), a + nv*(c - 1))*l2t1(i, j, c, k)
                     end do
                  end do
               end do
            end do
         end do
      end do
      diff = 0.0_dp
      scale = 1.0_dp
      raw = maxval(abs(old - new))
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  diff = max(diff, abs(old(i, j, a, b) + old(j, i, b, a) - new(i, j, a, b) - new(j, i, b, a)))
                  scale = max(scale, abs(old(i, j, a, b) + old(j, i, b, a)))
               end do
            end do
         end do
      end do
      call check(error, diff/scale < TOL_HOST, "the symmetrised identity: "//to_string(diff/scale, "es9.2"))
      if (.not. allocated(error)) then
         ! the two sides must differ before symmetrisation, or the case would not need it and so not test it
         call check(error, raw > 1.0e-6_dp, "the unsymmetrised sides agree: the t1 terms vanish")
      end if
      call lad%release()
      call ls%release()
   end subroutine test_identity

   subroutine test_transposed_refusals(error)
      !! Shape, share, `tau_symmetric` outside the half form, missing amplitudes and no room for the dressed
      !! factor: each by its message, nothing written
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_ladder_options_t) :: opts
      type(trc_rccsd_lscratch_t) :: ls
      type(trc_rccsd_ladder_t) :: lad
      type(trc_rccsd_amps_t) :: empty
      logical :: ok(5)

      call ensure_fixture()
      if (.not. supported) return
      opts%form = RCCSD_LADDER_ASSEMBLED
      call ls%init(ints, opts)
      call check(error, .not. ls%status%has_error(), "scratch: "//ls%status%get_message())
      if (allocated(error)) return
      call lad%init(ints%no, ints%nv)
      call rccsd_ladder_transposed(ints, t(:, :, :, 1:ints%nv - 1), amps, lad, ls, la)
      ok(1) = index(ls%status%get_message(), "does not have the shape") > 0
      call rccsd_ladder_transposed(ints, t, amps, lad, ls, la, tau_symmetric=.true.)
      ok(2) = index(ls%status%get_message(), "tau_symmetric") > 0
      call rccsd_ladder_transposed(ints, t, empty, lad, ls, la)
      ok(3) = index(ls%status%get_message(), "the amplitudes are not initialised") > 0
      call device_memory_set_cap(1000_int64)
      call rccsd_ladder_transposed(ints, t, amps, lad, ls, la)
      call device_memory_set_cap(0_int64)
      ok(4) = index(ls%status%get_message(), "the dressed factor of the transposed ladder (v 19, naux") > 0 .and. &
              index(ls%status%get_message(), "GB are free") > 0
      ls%d_lo = 2
      call rccsd_ladder_transposed(ints, t, amps, lad, ls, la)
      ! a share is itself refused by the scratch check of a form that cannot be shared
      ok(5) = index(ls%status%get_message(), "share of the ladder") > 0
      call lad%pull()
      call check(error, all(ok) .and. all(lad%h == 0.0_dp), "transposed refusals: shape "//merge("T", "F", ok(1))// &
                 " symmetric "//merge("T", "F", ok(2))//" amplitudes "//merge("T", "F", ok(3))//" memory "// &
                 merge("T", "F", ok(4))//" share "//merge("T", "F", ok(5)))
      call lad%release()
      call ls%release()
      if (allocated(error)) return
      opts%form = RCCSD_LADDER_HALF
      call ls%init(ints, opts)
      call lad%init(ints%no, ints%nv)
      ls%d_lo = 2
      call rccsd_ladder_transposed(ints, t, amps, lad, ls, la)
      call lad%pull()
      call check(error, index(ls%status%get_message(), "takes no share") > 0 .and. all(lad%h == 0.0_dp), &
                 "a share of the half form must be refused: "//ls%status%get_message())
      call lad%release()
      call ls%release()
   end subroutine test_transposed_refusals

   subroutine host_transposed(x, want)
      !! sum_cd sum_P bT(a,c;P) bT(b,d;P) x(i,j,c,d) by loops
      real(dp), intent(in) :: x(:, :, :, :)
      real(dp), allocatable, intent(out) :: want(:, :, :, :)

      real(dp), allocatable :: bt(:, :), w4(:, :), w(:, :)
      integer :: no, nv, i, j, a, b, c, d, k

      no = ints%no
      nv = ints%nv
      allocate (bt(nv*nv, ints%naux))
      do c = 1, nv
         do a = 1, nv
            bt(a + nv*(c - 1), :) = ints%b_vv(a + nv*(c - 1), :)
            do k = 1, no
               bt(a + nv*(c - 1), :) = bt(a + nv*(c - 1), :) - ints%b_ov(k + no*(a - 1), :)*amps%t1(k, c)
            end do
         end do
      end do
      w4 = matmul(bt, transpose(bt))
      allocate (want(no, no, nv, nv), w(nv, nv))
      do b = 1, nv
         do a = 1, nv
            do d = 1, nv
               do c = 1, nv
                  w(c, d) = w4(a + nv*(c - 1), b + nv*(d - 1))
               end do
            end do
            do j = 1, no
               do i = 1, no
                  want(i, j, a, b) = sum(w*x(i, j, :, :))
               end do
            end do
         end do
      end do
   end subroutine host_transposed

   subroutine ensure_fixture()
      type(trc_scf_result_t) :: res
      type(trc_rccsd_ints_options_t) :: bopts
      integer, allocatable :: zint(:), seed(:)
      integer :: nocc, nseed, dev, no, nv, a, b, i, j

      if (built) return
      built = .true.
      dev = trc_bind_device(0)
      call load_molecule(1, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) then
         print "(a)", "SKIP check_rccsd_ladder_tensor: this build's density-fitting kernels stop below the RIFIT shells"
         return
      end if
      call run_scf(bas, nocc, res)
      bopts%nalpha = nocc
      bopts%nbeta = nocc
      bopts%nfrozen = MOL_NFROZEN(1)
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      if (ints%status%has_error()) error stop "check_rccsd_ladder_tensor: build: "//ints%status%get_message()
      no = ints%no
      nv = ints%nv

      call random_seed(size=nseed)
      allocate (seed(nseed))
      seed = 20261124
      call random_seed(put=seed)
      allocate (t0(no, no, nv, nv), t(no, no, nv, nv))
      call random_number(t0)
      t0 = t0 - 0.5_dp
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  t(i, j, a, b) = 0.5_dp*(t0(i, j, a, b) + t0(j, i, b, a))
               end do
            end do
         end do
      end do
      !$acc enter data copyin(t, t0)
      call amps%init(no, nv)
      call random_number(amps%t1)
      amps%t1 = 0.2_dp*(amps%t1 - 0.5_dp)
      call amps%push()
      call host_reference(t, ref_l, ref_w)
      call host_reference(t0, ref0_l, ref0_w)
      call host_transposed(t, ref_t)
      call host_transposed(t0, ref0_t)
      ! the transposed dressing must be neither no dressing nor the dressed Wvvvv, or its cases would not tell
      if (maxval(abs(ref_t - ref_l)) < 1.0e-6_dp .or. maxval(abs(ref_t - ref_w)) < 1.0e-6_dp) &
         error stop "check_rccsd_ladder_tensor: the transposed dressing is not distinct"
      ! the two tensors must give different ladders, or the t0 cases would only repeat those of t
      if (maxval(abs(ref0_l - ref_l)) < 1.0e-6_dp) error stop "check_rccsd_ladder_tensor: t0 is symmetric"
   end subroutine ensure_fixture

   subroutine host_reference(x, want_l, want_w)
      !! L and L + D of `x` from the fitted factors, with b_vv(a + nv (c-1), P) and b_ov(k + no (c-1), P)
      real(dp), intent(in) :: x(:, :, :, :)
      real(dp), allocatable, intent(out) :: want_l(:, :, :, :), want_w(:, :, :, :)

      real(dp), allocatable :: vvvv(:, :), ovvv(:, :), w(:, :)
      integer :: no, nv, i, j, a, b, c, d, k

      no = ints%no
      nv = ints%nv
      ! vvvv((a,c), (b,d)) = (ac|bd), ovvv((k,c), (a,d)) = (kc|ad)
      vvvv = matmul(ints%b_vv, transpose(ints%b_vv))
      ovvv = matmul(ints%b_ov, transpose(ints%b_vv))
      allocate (want_l(no, no, nv, nv), want_w(no, no, nv, nv), w(nv, nv))
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  do d = 1, nv
                     do c = 1, nv
                        w(c, d) = vvvv(a + nv*(c - 1), b + nv*(d - 1))
                     end do
                  end do
                  want_l(i, j, a, b) = sum(w*x(i, j, :, :))
                  do d = 1, nv
                     do c = 1, nv
                        w(c, d) = 0.0_dp
                        do k = 1, no
                           w(c, d) = w(c, d) - ovvv(k + no*(d - 1), a + nv*(c - 1))*amps%t1(k, b) &
                                     - ovvv(k + no*(c - 1), b + nv*(d - 1))*amps%t1(k, a)
                        end do
                     end do
                  end do
                  want_w(i, j, a, b) = want_l(i, j, a, b) + sum(w*x(i, j, :, :))
               end do
            end do
         end do
      end do
      ! the dressing must be visible, or the dressed cases would only repeat the undressed ones
      if (maxval(abs(want_w - want_l)) < 1.0e-6_dp) error stop "check_rccsd_ladder_tensor: the dressing vanishes"
   end subroutine host_reference

   subroutine release_fixture()
      ! every release below is safe on what was never built: an unsupported build has only loaded `bas` and `aux`
      if (allocated(t)) then
         !$acc exit data delete(t, t0)
         deallocate (t, t0)
      end if
      call amps%release()
      call ints%release()
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine release_fixture

end module check_rccsd_ladder_tensor_suite

program check_rccsd_ladder_tensor
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_ladder_tensor_suite, only: collect_check_rccsd_ladder_tensor, release_fixture
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_ladder_tensor, error_unit, stat, parallel=.false.)
   call release_fixture()
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_ladder_tensor
