module check_rccsd_ladder_tensor_suite
   !!
   !! `rccsd_ladder_tensor` (src/post_hf/ccsd/trc_rccsd_ladder.F90), the ladder of an arbitrary o^2 v^2 tensor, on water
   !! cc-pVDZ with cc-pVDZ-RIFIT. The tensor is random with t(i,j,a,b) = t(j,i,b,a), the symmetry the Lambda amplitudes
   !! have, and t1 is random too, so nothing of a CCSD iteration (tau, MP2 amplitudes) is involved. The reference is
   !! host loops over b_vv and b_ov, no routine or helper of the ladder:
   !!    L(i,j,a,b) = sum_cd (ac|bd) t(i,j,c,d)
   !!    D(i,j,a,b) = - sum_cd sum_k [ (kd|ac) t1(k,b) + (kc|bd) t1(k,a) ] t(i,j,c,d)
   !! One test per form (assembled, direct, the half form's one- and two-contraction paths), each undressed against L
   !! and dressed against L + D, with the default batches and with small ones (which start GEMMs and gathers inside t).
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_ladder, only: trc_rccsd_ladder_t, trc_rccsd_lscratch_t, trc_rccsd_ladder_options_t, rccsd_ladder_tensor, &
                               RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_DIRECT, RCCSD_LADDER_HALF
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: MOL_NFROZEN, load_molecule, run_scf
   implicit none
   private
   public :: collect_check_rccsd_ladder_tensor, release_fixture

   real(dp), parameter :: TOL_HOST = 1.0e-12_dp  !! max abs difference over max(1, |reference|max)

   logical :: built = .false.
   logical :: supported = .true.
   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl
   type(trc_rccsd_ints_t) :: ints
   type(trc_rccsd_amps_t) :: amps
   type(trc_linalg_t) :: la
   real(dp), allocatable :: t(:, :, :, :)                      !! the tensor, mapped to the device
   real(dp), allocatable :: ref_l(:, :, :, :), ref_w(:, :, :, :)  !! L and L + D

contains

   subroutine collect_check_rccsd_ladder_tensor(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("assembled", test_assembled), &
                   new_unittest("direct", test_direct), &
                   new_unittest("half, one contraction", test_half_symmetric), &
                   new_unittest("half, two contractions", test_half_general), &
                   new_unittest("refusals", test_refusals)]
   end subroutine collect_check_rccsd_ladder_tensor

   subroutine test_assembled(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t) :: opts
      call ensure_fixture()
      opts%form = RCCSD_LADDER_ASSEMBLED
      call run_form(error, opts, .false.)
      if (allocated(error)) return
      opts%cd_batch = 7
      call run_form(error, opts, .false.)
   end subroutine test_assembled

   subroutine test_direct(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t) :: opts
      call ensure_fixture()
      opts%form = RCCSD_LADDER_DIRECT
      call run_form(error, opts, .false.)
      if (allocated(error)) return
      opts%p_batch = 5
      call run_form(error, opts, .false.)
   end subroutine test_direct

   subroutine test_half_symmetric(error)
      type(error_type), allocatable, intent(out) :: error
      call half_form(error, .true.)
   end subroutine test_half_symmetric

   subroutine test_half_general(error)
      type(error_type), allocatable, intent(out) :: error
      call half_form(error, .false.)
   end subroutine test_half_general

   subroutine half_form(error, sym)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: sym
      type(trc_rccsd_ladder_options_t) :: opts
      call ensure_fixture()
      opts%form = RCCSD_LADDER_HALF
      call run_form(error, opts, sym)
      if (allocated(error)) return
      opts%d_batch = 2
      call run_form(error, opts, sym)
   end subroutine half_form

   subroutine run_form(error, opts, sym)
      !! Undressed against L, dressed against L + D, for one plan
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ladder_options_t), intent(in) :: opts
      logical, intent(in) :: sym

      type(trc_rccsd_lscratch_t) :: ls
      type(trc_rccsd_ladder_t) :: lad
      real(dp) :: d_l, d_w

      if (.not. supported) return
      call ls%init(ints, opts)
      call check(error, .not. ls%status%has_error(), "scratch: "//ls%status%get_message())
      if (allocated(error)) return
      call lad%init(ints%no, ints%nv)

      call rccsd_ladder_tensor(ints, t, lad, ls, la, tau_symmetric=sym)
      call check(error, .not. ls%status%has_error(), "undressed: "//ls%status%get_message())
      if (.not. allocated(error)) then
         call lad%pull()
         d_l = maxval(abs(lad%h - ref_l))/max(1.0_dp, maxval(abs(ref_l)))
         call check(error, d_l < TOL_HOST, "undressed vs host loops: "//to_string(d_l, "es9.2"))
      end if

      if (.not. allocated(error)) then
         call rccsd_ladder_tensor(ints, t, lad, ls, la, amps=amps, tau_symmetric=sym)
         call check(error, .not. ls%status%has_error(), "dressed: "//ls%status%get_message())
      end if
      if (.not. allocated(error)) then
         call lad%pull()
         d_w = maxval(abs(lad%h - ref_w))/max(1.0_dp, maxval(abs(ref_w)))
         call check(error, d_w < TOL_HOST, "dressed vs host loops: "//to_string(d_w, "es9.2"))
      end if
      call lad%release()
      call ls%release()
   end subroutine run_form

   subroutine test_refusals(error)
      !! A tensor of another shape and a scratch cut to a share are refused and write nothing
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_ladder_options_t) :: opts
      type(trc_rccsd_lscratch_t) :: ls
      type(trc_rccsd_ladder_t) :: lad
      logical :: shape_refused, share_refused

      call ensure_fixture()
      if (.not. supported) return
      opts%form = RCCSD_LADDER_HALF
      call ls%init(ints, opts)
      call lad%init(ints%no, ints%nv)
      call rccsd_ladder_tensor(ints, t(:, :, :, 1:ints%nv - 1), lad, ls, la)
      shape_refused = ls%status%has_error()
      ls%d_lo = 2
      call rccsd_ladder_tensor(ints, t, lad, ls, la)
      share_refused = ls%status%has_error()
      call lad%pull()
      call check(error, shape_refused .and. share_refused .and. all(lad%h == 0.0_dp), &
                 "a wrong shape and a share must be refused with nothing written")
      call lad%release()
      call ls%release()
   end subroutine test_refusals

   subroutine ensure_fixture()
      type(trc_scf_result_t) :: res
      type(trc_rccsd_ints_options_t) :: bopts
      integer, allocatable :: zint(:), seed(:)
      real(dp), allocatable :: t0(:, :, :, :)
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
      !$acc enter data copyin(t)
      call amps%init(no, nv)
      call random_number(amps%t1)
      amps%t1 = 0.2_dp*(amps%t1 - 0.5_dp)
      call amps%push()
      call host_reference()
   end subroutine ensure_fixture

   subroutine host_reference()
      !! `ref_l` and `ref_w` from the fitted factors, with b_vv(a + nv (c-1), P) and b_ov(k + no (c-1), P)
      real(dp), allocatable :: vvvv(:, :), ovvv(:, :), w(:, :)
      integer :: no, nv, i, j, a, b, c, d, k

      no = ints%no
      nv = ints%nv
      ! vvvv((a,c), (b,d)) = (ac|bd), ovvv((k,c), (a,d)) = (kc|ad)
      vvvv = matmul(ints%b_vv, transpose(ints%b_vv))
      ovvv = matmul(ints%b_ov, transpose(ints%b_vv))
      allocate (ref_l(no, no, nv, nv), ref_w(no, no, nv, nv), w(nv, nv))
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  do d = 1, nv
                     do c = 1, nv
                        w(c, d) = vvvv(a + nv*(c - 1), b + nv*(d - 1))
                     end do
                  end do
                  ref_l(i, j, a, b) = sum(w*t(i, j, :, :))
                  do d = 1, nv
                     do c = 1, nv
                        w(c, d) = 0.0_dp
                        do k = 1, no
                           w(c, d) = w(c, d) - ovvv(k + no*(d - 1), a + nv*(c - 1))*amps%t1(k, b) &
                                     - ovvv(k + no*(c - 1), b + nv*(d - 1))*amps%t1(k, a)
                        end do
                     end do
                  end do
                  ref_w(i, j, a, b) = ref_l(i, j, a, b) + sum(w*t(i, j, :, :))
               end do
            end do
         end do
      end do
      ! the dressing must be visible, or the dressed cases would only repeat the undressed ones
      if (maxval(abs(ref_w - ref_l)) < 1.0e-6_dp) error stop "check_rccsd_ladder_tensor: the dressing vanishes"
   end subroutine host_reference

   subroutine release_fixture()
      if (.not. allocated(t)) return
      !$acc exit data delete(t)
      deallocate (t)
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
