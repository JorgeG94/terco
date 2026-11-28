!
! The intermediates of the closed-shell RI-CCSD Lambda equations that are built once per solve, from the
! converged amplitudes: what pyscf/cc/ccsd_lambda.py `make_intermediates` forms, without any o v^3 or v^4 array.
! Chemists' notation, spatial orbitals, REAL orbitals, canonical (f_ov = 0); the integral blocks and the
! amplitudes are those of trc_rccsd_ints and trc_rccsd_inter. The Lambda update itself is not here.
!
! WHAT IS THE CCSD'S OWN. Four of pyscf's intermediates are terco's CCSD one-particle ones at the same
! amplitudes, built here by `rccsd_inter_one_particle` into the caller's `inter` (which also holds tau):
!
!     v1(b,a) = inter%lvv(b,a)    v2(i,j) = inter%loo(i,j)    v4(j,b) = inter%fov(j,b)
!
! v1 and v2 therefore do NOT hold the orbital energies (pyscf's w1 and w2 are v1 + diag(eps_v) and
! v2 + diag(eps_o); its update subtracts the diagonals again). L(i,a,j,b) = 2 (ia|jb) - (ib|ja) is the
! caller's `scr%lov` and 2 (ij|ka) - (kj|ia) its `scr%lam(i,j,k,a)`.
!
! WHAT IS STORED (`trc_rccsd_lambda_inter_t`), first index fastest, and the einsum of each. Names: theta(ijab) =
! 2 t2(ijab) - t2(ijba), tau = t2 + t1 t1, (pq|rs) = sum_P b(pq,P) b(rs,P).
!
!   woooo(i,j,k,l)   = (ik|jl) + sum_c [(ik|jc) t1(l,c) + (jl|ic) t1(k,c)] + sum_cd (ic|jd) tau(klcd)
!                      (pyscf's woooo; terco's CCSD Woooo with the same index order)
!   Two ring combinations, from V(jbkc) = (jb|kc) + sum_ld [L(ld,jb) t2(klcd) - (ld|jb) t2(kldc)] and
!   U(jbkc) = sum_ld (jd|lb) t2(kldc) - (jk|bc):
!   w_both(j,b,k,c)  = 2 V + U
!   w_diff(j,b,k,c)  = V - U
!   and two t1-dressed ones, from
!   wOV(jbkc) = V - sum_dl (jb|ld) t1(k,d) t1(l,c) - sum_l (jb|lk) t1(l,c) + sum_d (jb|cd) t1(k,d)
!   wX(jbkc)  = U + sum_dl (jd|lb) t1(k,d) t1(l,c) + sum_l (lb|jk) t1(l,c) - sum_d (jd|cb) t1(k,d)
!   wovvo(j,b,k,c)   = 2 wOV + wX        (pyscf's wvOOv + 2 wVOov)
!   wovvo_x(j,b,k,c) = wX                (pyscf's wvOOv)
!   woovo(i,j,k,c)   : pyscf's wvooo(c,k,i,j); the nine terms are listed at `build_woovo`
!   w3(j,b)          = sum_kc w_both(kcjb) t1(k,c) + v5(b,j) + sum_c v1(b,c) t1(j,c) - sum_k v2(k,j) t1(k,b)
!                      + (eps_v(b) - eps_o(j)) t1(j,b),
!       v5(b,j) = sum_kc v4(k,c) t1(k,b) t1(j,c) - sum_klc lam(k,j,l,c) t2(klbc) + sum_kcd (kd|bc) theta(jkcd)
!   (pyscf's w3(b,j): it is formed with the whole w1, w2, so w3 keeps the orbital energies v1, v2 lack).
!
! The consumer (the Lambda update) therefore writes, with l2(i,j,a,b):
!   sum_kl woooo(i,j,k,l) l2(k,l,a,b)       sum_jkc woovo(i,j,k,c) l2(k,j,c,a)       sum_jb l2(i,j,a,b) w3(j,b)
!
! HOW THE o v^3 TERMS ARE DONE. No (ia|bc) block exists. With the DRESSED three-index quantities
!   Mhat(j,k;P) = b_oo(jk,P) + sum_d b_ov(jd,P) t1(k,d)
!   Q(k,c;P)    = n1(k,c;P) - sum_l Mhat(l,k;P) t1(l,c)          (n1 of `rccsd_inter_lvv`)
! the two dressed rings are
!   wOV(jbkc) = V(jbkc) + sum_P b_ov(jb,P) Q(kc,P)
!   wX(jbkc)  = [U(jbkc) + (jk|bc)] - sum_P Mhat(jk;P) b_vv(bc,P) + sum_lP Mhat(jk;P) b_ov(lb,P) t1(l,c)
! The (kd|bc) theta term of w3 goes through Z(j,c;P) = sum_kd b_ov(kd,P) theta(j,k,c,d). The term of woovo,
! sum_bd (id|cb) [2 tau(jkdb) - tau(jkbd)], is cut into blocks of the free virtual c: (id|bc) for the block's c is ONE
! GEMM of the factors (o v^2 per c, as the Z packs of the CCSD ladder dressing) and a strided batch over c
! contracts it with the packed tau.
!
! CONTRACTIONS. Every other product is a `trc_contract` plan. All plans are built before anything is written, so
! the scratch they need is known and checked with the rest; `options%extra_work_limit` is theirs.
!
! SIZES (o active occupied, v virtual, doubles of 8 bytes, host and device):
!   resident here    o^4 + 4 o^2 v^2 + o^3 v + o v
!   caller's         tau and the one-particle intermediates (`inter`), and `scr` (3 o^2 v^2 + o^3 v + 3 o v naux)
!   during `build`   the plans' scratch (`work_doubles`, a few o^2 v^2) + o^2 naux + o v naux + o^3 v + o^4
!                    + o v^2 woovo_batch
! `woovo_batch` is `options%woovo_batch`, else as many virtuals as 512 MiB, or half of what is free after the
! rest if that is less, hold; the environment variable TERCO_LAMBDA_WOOVO_CHUNK=<k> caps it (TESTING: forces
! the blocked path on a small molecule). Not even one virtual: refused with the sizes.
!
! MEASURED on one V100, adenine cc-pVDZ (no 25, nv 130, naux 630): 0.095 s and 460 GFLOP against 0.183 s for one
! CCSD iteration of the same run; the ring 0.031 s, the woovo blocks 0.037 s, the rest of woovo 0.011 s, the other
! steps 0.015 s together. Peak of this stage 1.4 GB (0.76 resident with `inter` and `scr`, 0.65 temporaries).
! Small woovo blocks cost time, not memory alone: the whole build took 0.095 s at blocks of 130 or 32 virtuals,
! 0.120 s at 8, 0.170 s at 4 and 0.109 s at 1.
!
! RESIDENCY. As in trc_rccsd_w: every array of the types exists on host and device, `build` writes the device
! copies, `pull` brings them back; kernels take explicit-shape dummies and never name `x%name`. `build` reads
! the device copies of the amplitudes.
!
! Failures are reported in `im%status`; every check is made before the first write, and none stops. No MPI call.
#include "inc/trc_dc_index.inc"

module trc_rccsd_lambda_inter
   use pic_types, only: dp, int32, int64
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t, trc_rccsd_scratch_t, rccsd_inter_one_particle, &
                              rccsd_inter_flops
   use trc_rccsd_steps, only: rccsd_steps_gemm_flops, rccsd_steps_open, rccsd_steps_close
   use trc_contract, only: trc_contract_plan_t, contract_permute
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_step_table, only: step_table_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_device_memory, only: device_memory_fits, device_memory_free_bytes
   use pic_timer, only: timer_type
   use pic_strings, only: to_string
   implicit none
   private

   public :: trc_rccsd_lambda_inter_t, trc_rccsd_lambda_options_t, rccsd_lambda_inter_build

   integer, parameter, public :: RCCSD_LAMBDA_ST_ONE = 1      !! tau, v1, v2, v4 (`rccsd_inter_one_particle`)
   integer, parameter, public :: RCCSD_LAMBDA_ST_WOOOO = 2    !! woooo
   integer, parameter, public :: RCCSD_LAMBDA_ST_RING = 3     !! V, U, w_both, w_diff
   integer, parameter, public :: RCCSD_LAMBDA_ST_DRESS = 4    !! Mhat, Q, wovvo, wovvo_x
   integer, parameter, public :: RCCSD_LAMBDA_ST_W3 = 5       !! w3
   integer, parameter, public :: RCCSD_LAMBDA_ST_WOOVO = 6    !! woovo without its (id|cb) tau term
   integer, parameter, public :: RCCSD_LAMBDA_ST_WOOVO_V = 7  !! the (id|cb) tau term of woovo, in blocks of c
   integer, parameter, public :: RCCSD_LAMBDA_NSTEP = 7
   character(len=24), parameter, public :: RCCSD_LAMBDA_STEP_NAME(RCCSD_LAMBDA_NSTEP) = [character(len=24) :: &
      "tau, one-particle", "woooo", "ring", "ring, t1 dressed", "w3", "woovo", "woovo (id|cb) blocks"]
   !! rows of the step table `rccsd_lambda_inter_build` fills

   integer(int64), parameter :: BLOCK_BUDGET = 536870912_int64  !! bytes of the woovo (ov|vv) block when not forced

   ! the contraction plans, by what they form
   integer, parameter :: P_OOOO_X = 1, P_OOOO_TAU = 2, P_V_L = 3, P_V_X = 4, P_U = 5, P_MHAT = 6, P_Q = 7, &
                         P_WOV = 8, P_X_VV = 9, P_X_OV = 10, P_X_T1 = 11, P_W3_BOTH = 12, P_W3_X = 13, &
                         P_W3_XT = 14, P_W3_LAM = 15, P_W3_Z = 16, P_W3_VV = 17, P_W3_V1 = 18, P_W3_V2 = 19, &
                         P_VO_A = 20, P_VO_B = 21, P_VO_C1 = 22, P_VO_C2 = 23, P_VO_D = 24, P_VO_E = 25, NPLAN = 25

   type :: trc_rccsd_lambda_options_t
      !! What `build` may be told. Plain values only.
      integer :: woovo_batch = 0  !! virtuals c per block of the woovo (id|cb) term; 0 chooses
      integer(int64) :: extra_work_limit = -1
         !! the contraction plans' `extra_work_limit` (doubles); negative keeps their default
   end type trc_rccsd_lambda_options_t

   type :: trc_rccsd_lambda_inter_t
      !! The Lambda intermediates on host and device, valid after `build` (module header)
      type(error_t) :: status  !! set by `init` and `build`: a refusal leaves its message here
      integer :: no = 0
      integer :: nv = 0
      logical :: built = .false.       !! the last `build` ran through
      integer :: woovo_batch = 0       !! virtuals per block the last `build` took
      integer(int64) :: work_doubles = 0  !! doubles of contraction scratch the last `build` allocated
      integer :: gemm_items = 0        !! GEMMs of the last `build`'s plans, the items of batched ones counted
      integer(int64) :: temp_bytes = 0    !! bytes of every temporary the last `build` allocated
      integer(int64) :: flops(RCCSD_LAMBDA_NSTEP) = 0  !! 2 m n k of the GEMMs of the last `build`, per step
      real(dp), allocatable :: woooo(:, :, :, :)    !! woooo(i,j,k,l)
      real(dp), allocatable :: w_both(:, :, :, :)   !! w_both(j,b,k,c)
      real(dp), allocatable :: w_diff(:, :, :, :)   !! w_diff(j,b,k,c)
      real(dp), allocatable :: wovvo(:, :, :, :)    !! wovvo(j,b,k,c)
      real(dp), allocatable :: wovvo_x(:, :, :, :)  !! wovvo_x(j,b,k,c)
      real(dp), allocatable :: woovo(:, :, :, :)    !! woovo(i,j,k,c)
      real(dp), allocatable :: w3(:, :)             !! w3(j,b)
   contains
      procedure :: init => lambda_init
      procedure :: release => lambda_release
      procedure :: pull => lambda_pull
   end type trc_rccsd_lambda_inter_t

   type :: work_t
      !! The plans and the temporaries of one `build`
      type(trc_contract_plan_t) :: plan(NPLAN)
      type(error_t) :: failure  !! the first contraction that refused
      real(dp), allocatable :: work(:)   !! the plans' scratch
      real(dp), allocatable :: mhat(:)   !! o^2 naux, Mhat(j,k;P)
      real(dp), allocatable :: q(:)      !! o v naux, Q(k,c;P), then Z(j,c;P)
      real(dp), allocatable :: hmat(:)   !! o^3 v, sum_P Mhat(jk;P) b_ov(lb,P)
      real(dp), allocatable :: xmat(:)   !! o^4: sum_c (ik|jc) t1(l,c), then (o^2 of it) sum_c t1(j,c) v4(k,c)
      real(dp), allocatable :: block(:)  !! o v^2 woovo_batch, (id|bc) for the block's c
   end type work_t

contains

   subroutine lambda_init(this, no, nv)
      !! The result arrays on host and device, zeroed; a refusal leaves `status` and the type empty
      class(trc_rccsd_lambda_inter_t), intent(inout) :: this
      integer, intent(in) :: no, nv

      integer(int64) :: o, v

      call this%release()
      call this%status%clear()
      if (no < 1 .or. nv < 1) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: no active occupied or no virtual orbitals")
         return
      end if
      o = no
      v = nv
      if (o*o*v*v > int(huge(0_int32), int64) .or. o*o*o*max(o, v) > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: a block has more than 2^31-1 elements")
         return
      end if
      if (.not. device_memory_fits(this%status, "trc_rccsd_lambda_inter: the intermediates (o "//to_string(no)//", v "// &
                                   to_string(nv)//")", 8_int64*(o**4 + 4_int64*o*o*v*v + o*o*o*v + o*v))) return
      this%no = no
      this%nv = nv
      allocate (this%woooo(no, no, no, no), this%w_both(no, nv, no, nv), this%w_diff(no, nv, no, nv))
      allocate (this%wovvo(no, nv, no, nv), this%wovvo_x(no, nv, no, nv), this%woovo(no, no, no, nv), this%w3(no, nv))
      this%woooo = 0.0_dp
      this%w_both = 0.0_dp
      this%w_diff = 0.0_dp
      this%wovvo = 0.0_dp
      this%wovvo_x = 0.0_dp
      this%woovo = 0.0_dp
      this%w3 = 0.0_dp
      !$acc enter data copyin(this%woooo, this%w_both, this%w_diff, this%wovvo, this%wovvo_x, this%woovo, this%w3)
   end subroutine lambda_init

   subroutine lambda_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_lambda_inter_t), intent(inout) :: this

      if (allocated(this%woooo)) then
         !$acc exit data delete(this%woooo, this%w_both, this%w_diff, this%wovvo, this%wovvo_x, this%woovo, this%w3)
         deallocate (this%woooo, this%w_both, this%w_diff, this%wovvo, this%wovvo_x, this%woovo, this%w3)
      end if
      this%no = 0
      this%nv = 0
      this%built = .false.
      this%woovo_batch = 0
      this%work_doubles = 0
      this%gemm_items = 0
      this%temp_bytes = 0
      this%flops = 0_int64
   end subroutine lambda_release

   subroutine lambda_pull(this)
      !! Device copy to host
      class(trc_rccsd_lambda_inter_t), intent(inout) :: this

      if (.not. allocated(this%woooo)) return
      !$acc update self(this%woooo, this%w_both, this%w_diff, this%wovvo, this%wovvo_x, this%woovo, this%w3)
   end subroutine lambda_pull

   subroutine rccsd_lambda_inter_build(ints, amps, inter, scr, im, la, options, stab)
      !! tau, v1, v2, v4 into `inter` and the Lambda intermediates into `im` (module header), for the amplitudes
      !! `amps` as they are on the device. `inter` and `scr` are initialised for these integrals (`scr` after the
      !! last `ints%build`), `im` by `init(no, nv)`. A refusal sets `im%status` and writes neither `im`'s arrays
      !! nor `inter`. `stab`, if given, has the rows RCCSD_LAMBDA_STEP_NAME.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_rccsd_lambda_inter_t), intent(inout) :: im
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_lambda_options_t), intent(in), optional :: options
      type(step_table_t), intent(inout), optional :: stab

      type(trc_rccsd_lambda_options_t) :: opts
      type(work_t) :: wk
      type(timer_type) :: tc

      call im%status%clear()
      im%built = .false.
      if (.not. ready(ints, amps, inter, scr, im)) return
      opts = trc_rccsd_lambda_options_t()
      if (present(options)) opts = options
      if (opts%woovo_batch < 0) then
         call im%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: a negative woovo_batch")
         return
      end if
      call plan_all(ints, opts, wk, im%status)
      if (im%status%has_error()) return
      call allocate_temporaries(ints, opts, wk, im)
      if (im%status%has_error()) return

      im%flops = 0_int64
      call rccsd_steps_open(stab, la, RCCSD_LAMBDA_ST_ONE, tc)
      call rccsd_inter_one_particle(ints, amps, inter, scr, la)
      im%flops(RCCSD_LAMBDA_ST_ONE) = rccsd_inter_flops(ints%no, ints%nv, ints%naux)
      call rccsd_steps_close(stab, la, RCCSD_LAMBDA_ST_ONE, tc, im%flops(RCCSD_LAMBDA_ST_ONE))
      if (scr%status%has_error()) then
         call im%status%set(ERROR_VALIDATION, scr%status%get_message())
         call release_temporaries(wk)
         return
      end if

      call rccsd_steps_open(stab, la, RCCSD_LAMBDA_ST_WOOOO, tc)
      call build_woooo(ints, amps, inter, scr, wk, im, la)
      call rccsd_steps_close(stab, la, RCCSD_LAMBDA_ST_WOOOO, tc, im%flops(RCCSD_LAMBDA_ST_WOOOO))
      call rccsd_steps_open(stab, la, RCCSD_LAMBDA_ST_RING, tc)
      call build_ring(ints, amps, scr, wk, im, la)
      call rccsd_steps_close(stab, la, RCCSD_LAMBDA_ST_RING, tc, im%flops(RCCSD_LAMBDA_ST_RING))
      call rccsd_steps_open(stab, la, RCCSD_LAMBDA_ST_DRESS, tc)
      call build_dressed(ints, amps, scr, wk, im, la)
      call rccsd_steps_close(stab, la, RCCSD_LAMBDA_ST_DRESS, tc, im%flops(RCCSD_LAMBDA_ST_DRESS))
      call rccsd_steps_open(stab, la, RCCSD_LAMBDA_ST_W3, tc)
      call build_w3(ints, amps, inter, scr, wk, im, la)
      call rccsd_steps_close(stab, la, RCCSD_LAMBDA_ST_W3, tc, im%flops(RCCSD_LAMBDA_ST_W3))
      call rccsd_steps_open(stab, la, RCCSD_LAMBDA_ST_WOOVO, tc)
      call build_woovo(amps, scr, wk, im, la)
      call rccsd_steps_close(stab, la, RCCSD_LAMBDA_ST_WOOVO, tc, im%flops(RCCSD_LAMBDA_ST_WOOVO))
      call rccsd_steps_open(stab, la, RCCSD_LAMBDA_ST_WOOVO_V, tc)
      call build_woovo_ovvv(ints, inter, scr, wk, im, la)
      call rccsd_steps_close(stab, la, RCCSD_LAMBDA_ST_WOOVO_V, tc, im%flops(RCCSD_LAMBDA_ST_WOOVO_V))

      call release_temporaries(wk)
      if (wk%failure%has_error()) then
         call im%status%set(ERROR_VALIDATION, wk%failure%get_message())
         return
      end if
      im%built = .true.
   end subroutine rccsd_lambda_inter_build

   function ready(ints, amps, inter, scr, im) result(ok)
      !! True when everything `build` reads and writes is there with the sizes of these integrals; else
      !! `im%status` says what is not
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_scratch_t), intent(in) :: scr
      type(trc_rccsd_lambda_inter_t), intent(inout) :: im
      logical :: ok

      integer :: no, nv

      ok = .false.
      if (.not. allocated(ints%oooo) .or. .not. allocated(ints%ooov) .or. .not. allocated(ints%oovv) .or. &
          .not. allocated(ints%ovov) .or. .not. allocated(ints%b_oo) .or. .not. allocated(ints%b_ov) .or. &
          .not. allocated(ints%b_vv)) then
         call im%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: the integral blocks are not built")
         return
      end if
      if (.not. allocated(amps%t1) .or. .not. allocated(amps%t2)) then
         call im%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: the amplitudes are not initialised")
         return
      end if
      if (.not. allocated(inter%tau) .or. .not. allocated(scr%lov) .or. .not. allocated(scr%lam) .or. &
          .not. allocated(scr%n1) .or. .not. allocated(scr%pack1)) then
         call im%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: the one-particle intermediates or their "// &
                            "scratch are not initialised")
         return
      end if
      if (.not. allocated(im%woooo)) then
         call im%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: the result is not initialised")
         return
      end if
      no = ints%no
      nv = ints%nv
      if (amps%no /= no .or. amps%nv /= nv .or. inter%no /= no .or. inter%nv /= nv .or. scr%no /= no .or. &
          scr%nv /= nv .or. im%no /= no .or. im%nv /= nv) then
         call im%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: (no, nv) differ between integrals, "// &
                            "amplitudes, intermediates, scratch and result")
         return
      end if
      if (size(scr%n1, kind=int64) /= int(no, int64)*nv*ints%naux) then
         call im%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: the scratch was sized for another auxiliary count")
         return
      end if
      if (scr%stamp /= ints%stamp) then
         call im%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: the integrals were rebuilt after the scratch "// &
                            "was initialised; re-init the scratch after ints%build")
         return
      end if
      ok = .true.
   end function ready

   subroutine plan_one(wk, k, spec, extent_a, extent_b, error)
      !! Builds plan `k`, keeping the first error
      type(work_t), intent(inout) :: wk
      integer, intent(in) :: k
      character(len=*), intent(in) :: spec
      integer, intent(in) :: extent_a(:), extent_b(:)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      call wk%plan(k)%build(spec, extent_a, extent_b, error)
   end subroutine plan_one

   subroutine plan_all(ints, opts, wk, error)
      !! Every contraction of `build`, before anything is allocated or written
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_lambda_options_t), intent(in) :: opts
      type(work_t), intent(inout) :: wk
      type(error_t), intent(inout) :: error

      integer :: o, v, x, k
      integer :: t1(2), t2(4), ovov(4), lam(4), bov(3), bvv(3), mhat(3), ovx(3)

      o = ints%no
      v = ints%nv
      x = ints%naux
      if (opts%extra_work_limit >= 0) then
         do k = 1, NPLAN
            wk%plan(k)%extra_work_limit = opts%extra_work_limit
         end do
      end if
      t1 = [o, v]
      t2 = [o, o, v, v]
      ovov = [o, v, o, v]
      lam = [o, o, o, v]
      bov = [o, v, x]
      bvv = [v, v, x]
      mhat = [o, o, x]
      ovx = [o, v, x]
      call plan_one(wk, P_OOOO_X, "ikjc,lc->ikjl", lam, t1, error)
      call plan_one(wk, P_OOOO_TAU, "icjd,klcd->ijkl", ovov, t2, error)
      call plan_one(wk, P_V_L, "ldjb,klcd->jbkc", ovov, t2, error)
      call plan_one(wk, P_V_X, "ldjb,kldc->jbkc", ovov, t2, error)
      call plan_one(wk, P_U, "jdlb,kldc->jbkc", ovov, t2, error)
      call plan_one(wk, P_MHAT, "jdP,kd->jkP", bov, t1, error)
      call plan_one(wk, P_Q, "lkP,lc->kcP", mhat, t1, error)
      call plan_one(wk, P_WOV, "jbP,kcP->jbkc", bov, ovx, error)
      call plan_one(wk, P_X_VV, "jkP,bcP->jkbc", mhat, bvv, error)
      call plan_one(wk, P_X_OV, "jkP,lbP->jklb", mhat, bov, error)
      call plan_one(wk, P_X_T1, "jklb,lc->jkbc", lam, t1, error)
      call plan_one(wk, P_W3_BOTH, "kcjb,kc->jb", ovov, t1, error)
      call plan_one(wk, P_W3_X, "jc,kc->jk", t1, t1, error)
      call plan_one(wk, P_W3_XT, "jk,kb->jb", [o, o], t1, error)
      call plan_one(wk, P_W3_LAM, "kjlc,klbc->jb", lam, t2, error)
      call plan_one(wk, P_W3_Z, "jkcd,kdP->jcP", t2, bov, error)
      call plan_one(wk, P_W3_VV, "jcP,bcP->jb", ovx, bvv, error)
      call plan_one(wk, P_W3_V1, "jc,bc->jb", t1, [v, v], error)
      call plan_one(wk, P_W3_V2, "kj,kb->jb", [o, o], t1, error)
      call plan_one(wk, P_VO_A, "ibkc,jb->ijkc", ovov, t1, error)
      call plan_one(wk, P_VO_B, "ibjc,kb->ijkc", ovov, t1, error)
      call plan_one(wk, P_VO_C1, "ijld,klcd->ijkc", lam, t2, error)
      call plan_one(wk, P_VO_C2, "ijld,kldc->ijkc", lam, t2, error)
      call plan_one(wk, P_VO_D, "ikld,jlcd->ijkc", lam, t2, error)
      call plan_one(wk, P_VO_E, "lkid,ljcd->ijkc", lam, t2, error)
   end subroutine plan_all

   subroutine allocate_temporaries(ints, opts, wk, im)
      !! Chooses the woovo block, checks that the temporaries fit and maps them (module header, SIZES). Not
      !! fitting sets `im%status` and allocates nothing.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_lambda_options_t), intent(in) :: opts
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_lambda_inter_t), intent(inout) :: im

      integer(int64) :: o, v, x, fixed, per_c, free, budget, work
      integer :: k, nb

      o = ints%no
      v = ints%nv
      x = ints%naux
      work = 1
      do k = 1, NPLAN
         work = max(work, wk%plan(k)%work_size)
      end do
      fixed = 8_int64*(work + o*o*x + o*v*x + o*o*o*v + o*o*o*o)
      per_c = 8_int64*o*v*v
      if (opts%woovo_batch > 0) then
         nb = int(min(int(opts%woovo_batch, int64), v))
      else
         budget = BLOCK_BUDGET
         free = device_memory_free_bytes()
         if (free >= 0) budget = min(budget, max(0_int64, (free - fixed)/2))
         nb = int(min(v, max(1_int64, budget/per_c)))
      end if
      nb = forced_cap(nb)
      ! o v nv nb: the extent of the block's GEMM and of the array, held to a default integer
      nb = int(min(int(nb, int64), max(1_int64, int(huge(0_int32), int64)/(o*v*v))))
      if (o*v*v > int(huge(0_int32), int64)) then
         call im%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: one virtual of the woovo block has more than "// &
                            "2^31-1 elements")
         return
      end if
      if (.not. device_memory_fits(im%status, "trc_rccsd_lambda_inter: the build's work arrays (o "//to_string(ints%no)// &
                                   ", v "//to_string(ints%nv)//", naux "//to_string(ints%naux)//", woovo block "// &
                                   to_string(nb)//")", fixed + per_c*nb)) return
      im%woovo_batch = nb
      im%work_doubles = work
      im%gemm_items = 0
      do k = 1, NPLAN
         im%gemm_items = im%gemm_items + wk%plan(k)%gemm%count
      end do
      im%temp_bytes = fixed + per_c*nb
      allocate (wk%work(work), wk%mhat(o*o*x), wk%q(o*v*x), wk%hmat(o*o*o*v), wk%xmat(o*o*o*o), wk%block(o*v*v*nb))
      !$acc enter data create(wk%work, wk%mhat, wk%q, wk%hmat, wk%xmat, wk%block)
   end subroutine allocate_temporaries

   function forced_cap(nb) result(capped)
      !! `nb` capped by TERCO_LAMBDA_WOOVO_CHUNK when that is a positive integer
      integer, intent(in) :: nb
      integer :: capped

      character(len=16) :: text
      integer :: stat, cap

      capped = nb
      call get_environment_variable("TERCO_LAMBDA_WOOVO_CHUNK", text, status=stat)
      if (stat /= 0) return
      read (text, *, iostat=stat) cap
      if (stat == 0 .and. cap >= 1) capped = min(nb, cap)
   end function forced_cap

   subroutine release_temporaries(wk)
      !! Unmaps and frees what `allocate_temporaries` made
      type(work_t), intent(inout) :: wk

      if (.not. allocated(wk%work)) return
      !$acc exit data delete(wk%work, wk%mhat, wk%q, wk%hmat, wk%xmat, wk%block)
      deallocate (wk%work, wk%mhat, wk%q, wk%hmat, wk%xmat, wk%block)
   end subroutine release_temporaries

   subroutine run(wk, k, la, a, b, c, scale)
      !! C = scale(1) contraction_k(A, B) + scale(2) C; a refusal is kept in `wk%failure`
      type(work_t), intent(inout) :: wk
      integer, intent(in) :: k
      type(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: a(*), b(*)
      real(dp), intent(inout) :: c(*)
      real(dp), intent(in) :: scale(2)

      call wk%plan(k)%apply(la, a, b, c, wk%work, scale(1), scale(2))
      if (wk%plan(k)%status%has_error() .and. .not. wk%failure%has_error()) then
         call wk%failure%set(ERROR_VALIDATION, "trc_rccsd_lambda_inter: "//wk%plan(k)%status%get_message())
      end if
   end subroutine run

   pure function plan_flops(plan) result(flops)
      !! 2 m n k of the GEMMs of one apply of `plan`
      type(trc_contract_plan_t), intent(in) :: plan
      integer(int64) :: flops

      flops = rccsd_steps_gemm_flops(plan%gemm%m, plan%gemm%n, plan%gemm%k, plan%gemm%count)
   end function plan_flops

   subroutine build_woooo(ints, amps, inter, scr, wk, im, la)
      !! woooo = oooo permuted + its two t1 terms + sum_cd (ic|jd) tau(klcd)
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_lambda_inter_t), intent(inout) :: im
      type(trc_linalg_t), intent(in) :: la

      ! X(i,k,j,l) = sum_c (ik|jc) t1(l,c)
      call run(wk, P_OOOO_X, la, ints%ooov, amps%t1, wk%xmat, [1.0_dp, 0.0_dp])
      call kernel_woooo(im%no, ints%oooo, wk%xmat, im%woooo)
      call run(wk, P_OOOO_TAU, la, ints%ovov, inter%tau, im%woooo, [1.0_dp, 1.0_dp])
      im%flops(RCCSD_LAMBDA_ST_WOOOO) = plan_flops(wk%plan(P_OOOO_X)) + plan_flops(wk%plan(P_OOOO_TAU))
   end subroutine build_woooo

   subroutine build_ring(ints, amps, scr, wk, im, la)
      !! V into im%wovvo and U + (jk|bc) into im%wovvo_x, then w_both and w_diff
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_scratch_t), intent(in) :: scr
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_lambda_inter_t), intent(inout) :: im
      type(trc_linalg_t), intent(in) :: la

      integer :: n

      n = im%no*im%nv*im%no*im%nv
      call kernel_axpby(n, 1.0_dp, ints%ovov, 0.0_dp, im%wovvo)
      call run(wk, P_V_L, la, scr%lov, amps%t2, im%wovvo, [1.0_dp, 1.0_dp])
      call run(wk, P_V_X, la, ints%ovov, amps%t2, im%wovvo, [-1.0_dp, 1.0_dp])
      call run(wk, P_U, la, ints%ovov, amps%t2, im%wovvo_x, [1.0_dp, 0.0_dp])
      call kernel_ring(im%no, im%nv, ints%oovv, im%wovvo, im%wovvo_x, im%w_both, im%w_diff)
      im%flops(RCCSD_LAMBDA_ST_RING) = plan_flops(wk%plan(P_V_L)) + plan_flops(wk%plan(P_V_X)) + &
                                       plan_flops(wk%plan(P_U))
   end subroutine build_ring

   subroutine build_dressed(ints, amps, scr, wk, im, la)
      !! Mhat, Q, then wOV in im%wovvo, wX in im%wovvo_x and wovvo = 2 wOV + wX (module header). Takes over
      !! im%wovvo = V and im%wovvo_x = U + (jk|bc) of `build_ring`.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_lambda_inter_t), intent(inout) :: im
      type(trc_linalg_t), intent(in) :: la

      type(error_t) :: err
      integer :: no, nv, naux

      no = ints%no
      nv = ints%nv
      naux = ints%naux
      call run(wk, P_MHAT, la, ints%b_ov, amps%t1, wk%mhat, [1.0_dp, 0.0_dp])
      call kernel_add(no*no*naux, ints%b_oo, wk%mhat)
      call kernel_axpby(no*nv*naux, 1.0_dp, scr%n1, 0.0_dp, wk%q)
      call run(wk, P_Q, la, wk%mhat, amps%t1, wk%q, [-1.0_dp, 1.0_dp])
      call run(wk, P_WOV, la, ints%b_ov, wk%q, im%wovvo, [1.0_dp, 1.0_dp])
      ! T(j,k,b,c) = sum_l H(j,k,l,b) t1(l,c) - sum_P Mhat(jk;P) b_vv(bc,P), added to wX as (j,b,k,c)
      call run(wk, P_X_VV, la, wk%mhat, ints%b_vv, scr%pack1, [-1.0_dp, 0.0_dp])
      call run(wk, P_X_OV, la, wk%mhat, ints%b_ov, wk%hmat, [1.0_dp, 0.0_dp])
      call run(wk, P_X_T1, la, wk%hmat, amps%t1, scr%pack1, [1.0_dp, 1.0_dp])
      call contract_permute("jkbc->jbkc", [no, no, nv, nv], scr%pack1, im%wovvo_x, err, 1.0_dp, 1.0_dp)
      if (err%has_error() .and. .not. wk%failure%has_error()) call wk%failure%set(ERROR_VALIDATION, err%get_message())
      ! wovvo = 2 wOV + wX
      call kernel_axpby(no*nv*no*nv, 1.0_dp, im%wovvo_x, 2.0_dp, im%wovvo)
      im%flops(RCCSD_LAMBDA_ST_DRESS) = plan_flops(wk%plan(P_MHAT)) + plan_flops(wk%plan(P_Q)) + &
                                        plan_flops(wk%plan(P_WOV)) + plan_flops(wk%plan(P_X_VV)) + &
                                        plan_flops(wk%plan(P_X_OV)) + plan_flops(wk%plan(P_X_T1))
   end subroutine build_dressed

   subroutine build_w3(ints, amps, inter, scr, wk, im, la)
      !! w3(j,b), module header. Overwrites wk%q (with Z) and scr%pack1 (with theta).
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_lambda_inter_t), intent(inout) :: im
      type(trc_linalg_t), intent(in) :: la

      call run(wk, P_W3_BOTH, la, im%w_both, amps%t1, im%w3, [1.0_dp, 0.0_dp])
      call run(wk, P_W3_X, la, amps%t1, inter%fov, wk%xmat, [1.0_dp, 0.0_dp])
      call run(wk, P_W3_XT, la, wk%xmat, amps%t1, im%w3, [1.0_dp, 1.0_dp])
      call run(wk, P_W3_LAM, la, scr%lam, amps%t2, im%w3, [-1.0_dp, 1.0_dp])
      call kernel_theta(im%no, im%nv, amps%t2, scr%pack1)
      call run(wk, P_W3_Z, la, scr%pack1, ints%b_ov, wk%q, [1.0_dp, 0.0_dp])
      call run(wk, P_W3_VV, la, wk%q, ints%b_vv, im%w3, [1.0_dp, 1.0_dp])
      call run(wk, P_W3_V1, la, amps%t1, inter%lvv, im%w3, [1.0_dp, 1.0_dp])
      call run(wk, P_W3_V2, la, inter%loo, amps%t1, im%w3, [-1.0_dp, 1.0_dp])
      call kernel_w3_eps(im%no, im%nv, ints%eps_o, ints%eps_v, amps%t1, im%w3)
      im%flops(RCCSD_LAMBDA_ST_W3) = plan_flops(wk%plan(P_W3_BOTH)) + plan_flops(wk%plan(P_W3_X)) + &
                                     plan_flops(wk%plan(P_W3_XT)) + plan_flops(wk%plan(P_W3_LAM)) + &
                                     plan_flops(wk%plan(P_W3_Z)) + plan_flops(wk%plan(P_W3_VV)) + &
                                     plan_flops(wk%plan(P_W3_V1)) + plan_flops(wk%plan(P_W3_V2))
   end subroutine build_w3

   subroutine build_woovo(amps, scr, wk, im, la)
      !! woovo(i,j,k,c) without the (id|cb) tau term (`build_woovo_ovvv` adds it). With lam(i,j,k,a) =
      !! 2 (ij|ka) - (kj|ia) and theta = 2 t2 - t2 swapped, the terms are, in the order below:
      !!     lam(j,i,k,c)
      !!   + sum_b [w_both(ib,kc) t1(j,b) - w_both(ib,jc) t1(k,b) + w_diff(ib,jc) t1(k,b)]
      !!   + sum_ld lam(i,j,l,d) theta(k,l,c,d)
      !!   - sum_ld lam(i,k,l,d) t2(j,l,c,d)
      !!   - sum_ld lam(l,k,i,d) t2(l,j,c,d)
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_scratch_t), intent(in) :: scr
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_lambda_inter_t), intent(inout) :: im
      type(trc_linalg_t), intent(in) :: la

      type(error_t) :: err
      integer :: no, nv

      no = im%no
      nv = im%nv
      call contract_permute("jikc->ijkc", [no, no, no, nv], scr%lam, im%woovo, err)
      if (err%has_error() .and. .not. wk%failure%has_error()) call wk%failure%set(ERROR_VALIDATION, err%get_message())
      call run(wk, P_VO_A, la, im%w_both, amps%t1, im%woovo, [1.0_dp, 1.0_dp])
      call run(wk, P_VO_B, la, im%w_both, amps%t1, im%woovo, [-1.0_dp, 1.0_dp])
      call run(wk, P_VO_B, la, im%w_diff, amps%t1, im%woovo, [1.0_dp, 1.0_dp])
      call run(wk, P_VO_C1, la, scr%lam, amps%t2, im%woovo, [2.0_dp, 1.0_dp])
      call run(wk, P_VO_C2, la, scr%lam, amps%t2, im%woovo, [-1.0_dp, 1.0_dp])
      call run(wk, P_VO_D, la, scr%lam, amps%t2, im%woovo, [-1.0_dp, 1.0_dp])
      call run(wk, P_VO_E, la, scr%lam, amps%t2, im%woovo, [-1.0_dp, 1.0_dp])
      im%flops(RCCSD_LAMBDA_ST_WOOVO) = plan_flops(wk%plan(P_VO_A)) + 2_int64*plan_flops(wk%plan(P_VO_B)) + &
                                        plan_flops(wk%plan(P_VO_C1)) + plan_flops(wk%plan(P_VO_C2)) + &
                                        plan_flops(wk%plan(P_VO_D)) + plan_flops(wk%plan(P_VO_E))
   end subroutine build_woovo

   subroutine build_woovo_ovvv(ints, inter, scr, wk, im, la)
      !! woovo(i,j,k,c) += sum_bd (id|cb) [2 tau(jkdb) - tau(jkbd)], `im%woovo_batch` virtuals c at a time
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_lambda_inter_t), intent(inout) :: im
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: shape
      integer :: no, nv, nov, nv2, c0, nc

      no = ints%no
      nv = ints%nv
      nov = no*nv
      nv2 = nv*nv
      call kernel_ttau(no, nv, inter%tau, scr%pack1)
      do c0 = 1, nv, im%woovo_batch
         nc = min(im%woovo_batch, nv - c0 + 1)
         ! block(i,d,b,c) = sum_P b_ov(id,P) b_vv(bc,P): the rows (b,c) of b_vv for the block's c are contiguous
         call la%gemm("N", "T", nov, nv*nc, ints%naux, 1.0_dp, ints%b_ov, nov, ints%b_vv(1 + nv*(c0 - 1), 1), nv2, &
                      0.0_dp, wk%block, nov)
         ! woovo(i,jk,c) += sum_db block(i,db,c) ttau(db,jk), one item per c, ttau shared
         shape = trc_gemm_batch(no, no*no, nv2, nc, lda=no, ldb=nv2, ldc=no, sa=int(nov, int64)*nv, sb=0_int64, &
                                sc=int(no, int64)*no*no)
         call la%gemm_batched(shape, wk%block, scr%pack1, im%woovo(1, 1, 1, c0), 1.0_dp, 1.0_dp)
         im%flops(RCCSD_LAMBDA_ST_WOOVO_V) = im%flops(RCCSD_LAMBDA_ST_WOOVO_V) + &
                                             rccsd_steps_gemm_flops(nov, nv*nc, ints%naux) + &
                                             rccsd_steps_gemm_flops(shape%m, shape%n, shape%k, shape%count)
      end do
   end subroutine build_woovo_ovvv

   subroutine kernel_axpby(n, a, x, b, y)
      !! y = a x + b y; y is not read when b = 0
      integer, intent(in) :: n
      real(dp), intent(in) :: a, b
      real(dp), intent(in) :: x(n)
      real(dp), intent(inout) :: y(n)
      integer :: t

      if (abs(b) <= 0.0_dp) then
         do concurrent(t=1:n)
            y(t) = a*x(t)
         end do
      else
         do concurrent(t=1:n)
            y(t) = a*x(t) + b*y(t)
         end do
      end if
   end subroutine kernel_axpby

   subroutine kernel_add(n, x, y)
      !! y += x
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(inout) :: y(n)
      integer :: t

      do concurrent(t=1:n)
         y(t) = y(t) + x(t)
      end do
   end subroutine kernel_add

   subroutine kernel_woooo(no, oooo, x, w)
      !! w(i,j,k,l) = (ik|jl) + x(i,k,j,l) + x(j,l,i,k)
      integer, intent(in) :: no
      real(dp), intent(in) :: oooo(no, no, no, no), x(no, no, no, no)
      real(dp), intent(out) :: w(no*no*no*no)
      integer :: t

      do concurrent(t=1:no*no*no*no)
         w(t) = woooo_item(TRC_DC_INDEX(t), no, oooo, x)
      end do
   end subroutine kernel_woooo

   subroutine kernel_ring(no, nv, oovv, v, u, both, diff)
      !! With U = u - (jk|bc) at (j,b,k,c): both = 2 v + U, diff = v - U
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: oovv(no, no, nv, nv)
      real(dp), intent(in) :: v(no*nv*no*nv), u(no*nv*no*nv)
      real(dp), intent(out) :: both(no*nv*no*nv), diff(no*nv*no*nv)
      integer :: t

      do concurrent(t=1:no*nv*no*nv)
         both(t) = 2.0_dp*v(t) + (u(t) - oovv_item(TRC_DC_INDEX(t), no, nv, oovv))
         diff(t) = v(t) - (u(t) - oovv_item(TRC_DC_INDEX(t), no, nv, oovv))
      end do
   end subroutine kernel_ring

   subroutine kernel_theta(no, nv, t2, theta)
      !! theta(j,k,c,d) = 2 t2(j,k,c,d) - t2(j,k,d,c)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(no, no, nv, nv)
      real(dp), intent(out) :: theta(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         theta(t) = theta_item(TRC_DC_INDEX(t), no, nv, t2)
      end do
   end subroutine kernel_theta

   subroutine kernel_ttau(no, nv, tau, ttau)
      !! ttau(d,b,j,k) = 2 tau(j,k,d,b) - tau(j,k,b,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(out) :: ttau(nv*nv*no*no)
      integer :: t

      do concurrent(t=1:nv*nv*no*no)
         ttau(t) = ttau_item(TRC_DC_INDEX(t), no, nv, tau)
      end do
   end subroutine kernel_ttau

   subroutine kernel_w3_eps(no, nv, eps_o, eps_v, t1, w3)
      !! w3(j,b) += (eps_v(b) - eps_o(j)) t1(j,b)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: eps_o(no), eps_v(nv), t1(no, nv)
      real(dp), intent(inout) :: w3(no, nv)
      integer :: t

      do concurrent(t=1:no*nv)
         w3(mod(t - 1, no) + 1, (t - 1)/no + 1) = w3(mod(t - 1, no) + 1, (t - 1)/no + 1) + &
            (eps_v((t - 1)/no + 1) - eps_o(mod(t - 1, no) + 1))*t1(mod(t - 1, no) + 1, (t - 1)/no + 1)
      end do
   end subroutine kernel_w3_eps

   pure subroutine split4(t, n1, n2, n3, i1, i2, i3, i4)
      !! t - 1 = (i1-1) + n1 ((i2-1) + n2 ((i3-1) + n3 (i4-1)))
      !$acc routine seq
      integer, intent(in) :: t, n1, n2, n3
      integer, intent(out) :: i1, i2, i3, i4
      integer :: q

      q = t - 1
      i1 = mod(q, n1) + 1
      q = q/n1
      i2 = mod(q, n2) + 1
      q = q/n2
      i3 = mod(q, n3) + 1
      i4 = q/n3 + 1
   end subroutine split4

   pure function woooo_item(t, no, oooo, x) result(w)
      !! woooo(i,j,k,l) without its tau term, at position t of (no, no, no, no)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no
      real(dp), intent(in) :: oooo(no, no, no, no), x(no, no, no, no)
      real(dp) :: w
      integer :: i, j, k, l

      call split4(t, no, no, no, i, j, k, l)
      w = oooo(i, k, j, l) + x(i, k, j, l) + x(j, l, i, k)
   end function woooo_item

   pure function oovv_item(t, no, nv, oovv) result(w)
      !! (jk|bc) at position t of (j,b,k,c) in (no, nv, no, nv)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: oovv(no, no, nv, nv)
      real(dp) :: w
      integer :: j, b, k, c

      call split4(t, no, nv, no, j, b, k, c)
      w = oovv(j, k, b, c)
   end function oovv_item

   pure function theta_item(t, no, nv, t2) result(w)
      !! theta(j,k,c,d) at position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(no, no, nv, nv)
      real(dp) :: w
      integer :: j, k, c, d

      call split4(t, no, no, nv, j, k, c, d)
      w = 2.0_dp*t2(j, k, c, d) - t2(j, k, d, c)
   end function theta_item

   pure function ttau_item(t, no, nv, tau) result(w)
      !! ttau(d,b,j,k) at position t of (nv, nv, no, no)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp) :: w
      integer :: d, b, j, k

      call split4(t, nv, nv, no, d, b, j, k)
      w = 2.0_dp*tau(j, k, d, b) - tau(j, k, b, d)
   end function ttau_item

end module trc_rccsd_lambda_inter
