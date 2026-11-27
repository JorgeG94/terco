!
! One update of the closed-shell, density-fitted CCSD Lambda equations: from (l1, l2), the converged CCSD
! amplitudes and the intermediates of trc_rccsd_lambda_inter, the next (l1, l2), what pyscf/cc/ccsd_lambda.py
! `update_lambda` gives, without any o v^3 or v^4 array. Chemists' notation, spatial orbitals, REAL orbitals,
! canonical (f_ov = 0). The iteration is trc_rccsd_lambda.
!
! CONVENTION (pyscf's). With r(t) = D (t' - t) the CCSD residual of trc_rccsd_update (t' its update,
! D1 = e_i - e_a, D2 = e_i + e_j - e_a - e_b), (l1, l2) make
!     Lag(t) = E(t) + 2 sum_ia l1(i,a) r1(i,a) + sum_ijab [2 l2(i,j,a,b) - l2(i,j,b,a)] r2(i,j,a,b)
! stationary in t. The update is affine in (l1, l2): its residual D (l' - l) is half of dLag/dt1 for l1, and
! for l2 the R with 2 R(i,j,a,b) - R(i,j,b,a) = dLag/dt2 (pair-symmetric variations). These identities hold at
! CONVERGED CCSD amplitudes (r(t) = 0) only: as pyscf's, the l2-dependent terms use the amplitude equations
! (check_rccsd_lambda_update: the finite difference of Lag agrees to 1e-14 there, and differs at order one away
! from them once l2 is non-zero). The terms that vanish with r(t), w3 l2 among them, are tested at random
! amplitudes against the dense equations.
!
! MEASURED on one V100, adenine cc-pVDZ / cc-pVDZ-RIFIT (no 25, nv 130, naux 630), 0.138 s without the DIIS:
! the ladder 0.068 s, the two rings 0.021 s, N and Q 0.014 s, A and Y 0.013 s, the woooo and one-body terms
! 0.009 s, the rest 0.012 s; 6 TFLOP/s overall.
!
! l2 MUST be pair-symmetric, l2(j,i,b,a) = l2(i,j,a,b), as t2 is: the ladder takes its one-contraction path and
! the equations use the symmetry to swap index pairs. `step` does not check it. The result is pair-symmetric
! exactly (it is x + x^T over a symmetric denominator).
!
! THE TERMS. v1, v2, v4 are the CCSD's Lvv, Loo, Fov (no orbital energies), the stored intermediates are those of
! trc_rccsd_lambda_inter, theta(j,k,c,d) = 2 t2(j,k,c,d) - t2(j,k,d,c), L = 2 (ia|jb) - (ib|ja),
! lam(i,j,k,a) = 2 (ij|ka) - (kj|ia). The one-particle quantities of l are
!     mvv(b,a) = sum_klc l2(klca) theta(klcb)        mvv1(b,c) = mvv(b,c) + sum_j t1(j,b) l1(j,c)
!     moo(i,j) = sum_kcd l2(kicd) theta(kjcd)        moo1(i,k) = moo(i,k) + sum_c l1(i,c) t1(k,c)
! The l2 right-hand side is g(i,j,a,b) + g(j,i,b,a) over D2, with g summed from
!   - m3 = 1/2 R + 1/2 sum_kl woooo(ijkl) l2(klab) + 1/2 sum_kl (ka|lb) [sum_cd l2(ijcd) t2(klcd)],
!     R = rccsd_ladder_transposed(l2). After the symmetrisation this is pyscf's bare ladder, woooo term,
!     1/2 ovov (l2 tau) and (kb|ac) slab term together (trc_rccsd_ladder, TRANSPOSED DRESSING);
!   - the terms with v1, v2 and (jk|ia) l1;
!   - the two rings with wovvo and wovvo_x;
!   - 1/2 (ia|jb) + l1(ia) v4(jb) - (ia|jc) mvv1(c,b) - moo1(i,k) (ka|jb) + (ia|cb) l1(j,c), the last through
!     the three-index factors.
! The l1 right-hand side is over D1. Its terms with (ia|bc) all end in one contraction of an o v naux array with
! b_vv, plus one of an o^2 naux array with b_ov; wvvvo . l2 is never formed:
!     N(i,k,b,d) = - sum_jc [l2(ijbc) - 1/2 l2(ijcb)] theta(jkcd)
!     A(i,k,c,d) = 3/2 sum_jb l2(ijbc) t2(kjbd) - N(i,k,c,d)
!     Q(i,b;P)   = sum_c b_ov(ic,P) mvv1(c,b) + sum_kd [2 N(i,k,b,d) - 2 l2(ikbd) + l2(ikdb)] b_ov(kd,P)
!     Y(i,k;P)   = sum_cd A(i,k,c,d) b_vv(cd,P)
! entering as - sum_bP Q(i,b;P) b_vv(ab,P) - sum_kP Y(i,k;P) b_ov(ka,P).
!
! SIZES (o active occupied, v virtual, doubles of 8 bytes, host and device each). `init` holds, beyond the
! integrals and amplitudes: the CCSD one-particle intermediates and their scratch, the Lambda intermediates
! (trc_rccsd_lambda_inter, SIZES), the ladder and its scratch, theta, w2, w3 (o^2 v^2 each), s1, s2 (o^3 v),
! x1, x2 (o^4) and y (o^2 naux). BORROWED from the one-particle scratch, which nothing reads after the
! intermediates are built: scr%pack1 (o^2 v^2, here w1), scr%n1p (o v naux: Q, and the dressed factor of the
! l1 term of l2), scr%u (naux). `step` allocates nothing, except that the assembled and direct ladders take
! v^2 naux at their first call (rccsd_ladder_transposed). Every contraction above o^2 v^2 is a GEMM; o^3 v^3
! flops come from the ladder's correction-free path, the woooo term, the two rings and N, A.
!
! RESIDENCY. As in trc_rccsd_update: the arrays exist on host and device, the kernels read and write the device
! copies, `step` leaves `lnew` on the device.
!
! Failures are reported in `status`; none stops. No MPI call: one rank, the whole ladder.
#include "inc/trc_dc_index.inc"

module trc_rccsd_lambda_update
   use pic_types, only: dp, int32, int64
   use pic_timer, only: timer_type
   use pic_strings, only: to_string
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t, trc_rccsd_scratch_t
   use trc_rccsd_lambda_inter, only: trc_rccsd_lambda_inter_t, trc_rccsd_lambda_options_t, rccsd_lambda_inter_build
   use trc_rccsd_ladder, only: trc_rccsd_ladder_t, trc_rccsd_lscratch_t, trc_rccsd_ladder_options_t, &
                               rccsd_ladder_transposed, RCCSD_LADDER_HALF
   use trc_rccsd_steps, only: rccsd_steps_gemm_flops, rccsd_steps_open, rccsd_steps_close
   use trc_contract, only: contract_permute
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_step_table, only: step_table_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_device_memory, only: device_memory_fits
   implicit none
   private

   public :: trc_rccsd_lambda_update_t, trc_rccsd_lambda_update_options_t

   integer, parameter, public :: RCCSD_LUPD_ST_ONE = 1     !! mvv, moo, mvv1, moo1
   integer, parameter, public :: RCCSD_LUPD_ST_LADDER = 2  !! the ladder on l2
   integer, parameter, public :: RCCSD_LUPD_ST_OOVV = 3    !! m3, the v1, v2 and (jk|ia) l1 terms of l2
   integer, parameter, public :: RCCSD_LUPD_ST_RING = 4    !! the two rings of l2
   integer, parameter, public :: RCCSD_LUPD_ST_L2 = 5      !! the rest of l2, symmetrisation and denominators
   integer, parameter, public :: RCCSD_LUPD_ST_Q = 6       !! N and Q, and Q's contraction with b_vv
   integer, parameter, public :: RCCSD_LUPD_ST_A = 7       !! A and Y, and Y's contraction with b_ov
   integer, parameter, public :: RCCSD_LUPD_ST_L1 = 8      !! the rest of l1
   integer, parameter, public :: RCCSD_LUPD_NSTEP = 8
   character(len=24), parameter, public :: RCCSD_LUPD_STEP_NAME(RCCSD_LUPD_NSTEP) = [character(len=24) :: &
      "l one-particle", "ladder on l2", "l2 woooo, one-body", "l2 ring", "l2 rest, denominators", &
      "l1 N, Q", "l1 A, Y", "l1 rest"]
   !! rows of the step table `step` fills; a caller's table may have more rows after these

   type :: trc_rccsd_lambda_update_options_t
      !! What `init` may be told. Plain values only.
      type(trc_rccsd_ladder_options_t) :: ladder  !! form and batch sizes of the ladder on l2
      type(trc_rccsd_lambda_options_t) :: inter   !! passed to the build of the intermediates
   end type trc_rccsd_lambda_update_options_t

   type :: trc_rccsd_lambda_update_t
      !! Everything one Lambda update needs besides the integrals, the CCSD amplitudes and (l1, l2)
      type(error_t) :: status  !! set by `init` and `step`: a refusal leaves its message here
      integer :: no = 0
      integer :: nv = 0
      integer :: naux = 0
      integer :: stamp = 0            !! `ints%stamp` of the integrals `init` built from
      integer(int64) :: bytes = 0     !! device bytes of the arrays below that `init` allocated itself
      type(trc_rccsd_inter_t) :: inter       !! tau, v1 (lvv), v2 (loo), v4 (fov)
      type(trc_rccsd_scratch_t) :: scr       !! lov, lam, bovp; pack1, n1p and u are borrowed (module header)
      type(trc_rccsd_lambda_inter_t) :: im   !! the once-per-solve intermediates
      type(trc_rccsd_ladder_t) :: lad        !! the ladder on l2, then the running g
      type(trc_rccsd_lscratch_t) :: ls
      real(dp), allocatable :: theta(:)  !! o^2 v^2, theta(j,k,c,d)
      real(dp), allocatable :: w2(:)     !! o^2 v^2 work
      real(dp), allocatable :: w3(:)     !! o^2 v^2 work
      real(dp), allocatable :: s1(:)     !! o^3 v work
      real(dp), allocatable :: s2(:)     !! o^3 v work
      real(dp), allocatable :: x1(:)     !! o^4 work
      real(dp), allocatable :: x2(:)     !! o^4 work
      real(dp), allocatable :: y(:)      !! o^2 naux, Y(i,k;P)
      real(dp), allocatable :: mvv(:, :)   !! mvv(b,a)
      real(dp), allocatable :: mvv1(:, :)  !! mvv1(b,c)
      real(dp), allocatable :: moo(:, :)   !! moo(i,j)
      real(dp), allocatable :: moo1(:, :)  !! moo1(i,k)
      real(dp), allocatable :: moot(:, :)  !! moo1 transposed
      real(dp), allocatable :: lvec(:, :)   !! lvec(j,b), what L multiplies in l1
   contains
      procedure :: init => update_init
      procedure :: release => update_release
      procedure :: step => update_step
   end type trc_rccsd_lambda_update_t

contains

   subroutine update_init(this, ints, amps, la, options, stab)
      !! Builds the intermediates for the CCSD amplitudes `amps` (on the device; converged, with pair-symmetric t2)
      !! and allocates the work arrays. A refusal is in `status` and leaves the type empty. `stab`, if given, has
      !! the rows RCCSD_LAMBDA_STEP_NAME of the intermediates' build. The integrals and `amps` must stay as they
      !! are for the life of the object: `step` reads them and nothing but the stamp of `ints` is checked.
      class(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_lambda_update_options_t), intent(in), optional :: options
      type(step_table_t), intent(inout), optional :: stab

      type(trc_rccsd_lambda_update_options_t) :: opts
      type(error_t) :: err

      call this%release()
      call this%status%clear()
      opts = trc_rccsd_lambda_update_options_t()
      if (present(options)) opts = options
      if (.not. init_ready(this, ints, amps)) return
      call this%inter%init(ints%no, ints%nv)
      if (this%inter%status%has_error()) call fail(this, this%inter%status)
      if (this%status%has_error()) return
      call this%scr%init(ints)
      if (this%scr%status%has_error()) call fail(this, this%scr%status)
      if (this%status%has_error()) return
      call this%im%init(ints%no, ints%nv)
      if (this%im%status%has_error()) call fail(this, this%im%status)
      if (this%status%has_error()) return
      ! before the ladder's scratch exists, so the build's temporaries do not stand on top of it
      call rccsd_lambda_inter_build(ints, amps, this%inter, this%scr, this%im, la, opts%inter, stab)
      if (this%im%status%has_error()) call fail(this, this%im%status)
      if (this%status%has_error()) return
      call this%lad%init(ints%no, ints%nv)
      if (this%lad%status%has_error()) call fail(this, this%lad%status)
      if (this%status%has_error()) return
      ! rccsd_ladder_transposed never runs the dressing, whose packs the plan allocates all the same
      if (opts%ladder%virtual_batch == 0) opts%ladder%virtual_batch = 1
      call this%ls%init(ints, opts%ladder)
      if (this%ls%status%has_error()) call fail(this, this%ls%status)
      if (this%status%has_error()) return
      call allocate_work(this, ints)
      if (this%status%has_error()) return
      call contract_permute("jkdc->jkcd", [ints%no, ints%no, ints%nv, ints%nv], amps%t2, this%theta, err, -1.0_dp)
      call contract_permute("jkcd->jkcd", [ints%no, ints%no, ints%nv, ints%nv], amps%t2, this%theta, err, 2.0_dp, 1.0_dp)
      if (err%has_error()) then
         call fail(this, err)
         return
      end if
      this%stamp = ints%stamp
   end subroutine update_init

   function init_ready(this, ints, amps) result(ok)
      !! True when the integrals are built and the amplitudes have their sizes; else `this%status` says which not
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      logical :: ok

      ok = .false.
      if (.not. allocated(ints%b_ov) .or. .not. allocated(ints%b_vv) .or. .not. allocated(ints%ovov) .or. &
          ints%stamp == 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: the integral blocks are not built")
         return
      end if
      if (ints%no < 1 .or. ints%nv < 1) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: no active occupied or no virtual orbitals")
         return
      end if
      if (.not. allocated(amps%t1) .or. .not. allocated(amps%t2)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: the CCSD amplitudes are not initialised")
         return
      end if
      if (amps%no /= ints%no .or. amps%nv /= ints%nv) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: (no, nv) differ between integrals and amplitudes")
         return
      end if
      ok = .true.
   end function init_ready

   subroutine allocate_work(this, ints)
      !! The work arrays, after a check that they fit (module header, SIZES)
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints

      integer(int64) :: o, v, x, doubles
      type(error_t) :: refusal
      integer :: no, nv

      no = ints%no
      nv = ints%nv
      o = no
      v = nv
      x = ints%naux
      doubles = 3_int64*o*o*v*v + 2_int64*o*o*o*v + 2_int64*o*o*o*o + o*o*x + 2_int64*v*v + 3_int64*o*o + o*v
      ! every GEMM extent below is a default integer
      if (o*o*v*v > int(huge(0_int32), int64) .or. v*x > int(huge(0_int32), int64) .or. &
          o*x > int(huge(0_int32), int64) .or. o*o*o*max(o, v) > int(huge(0_int32), int64) .or. &
          o*o*x > int(huge(0_int32), int64)) then
         call refusal%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: a block has more than 2^31-1 elements (o "// &
                          to_string(no)//", v "//to_string(nv)//", naux "//to_string(ints%naux)//")")
         call fail(this, refusal)
         return
      end if
      if (.not. device_memory_fits(refusal, "trc_rccsd_lambda_update: the work arrays (o "//to_string(no)//", v "// &
                                   to_string(nv)//", naux "//to_string(ints%naux)//")", 8_int64*doubles)) then
         call fail(this, refusal)
         return
      end if
      this%no = no
      this%nv = nv
      this%naux = ints%naux
      this%bytes = 8_int64*doubles
      allocate (this%theta(no*no*nv*nv), this%w2(no*no*nv*nv), this%w3(no*no*nv*nv))
      allocate (this%s1(no*no*no*nv), this%s2(no*no*no*nv), this%x1(no*no*no*no), this%x2(no*no*no*no))
      allocate (this%y(no*no*ints%naux))
      allocate (this%mvv(nv, nv), this%mvv1(nv, nv), this%moo(no, no), this%moo1(no, no), this%moot(no, no))
      allocate (this%lvec(no, nv))
      !$acc enter data create(this%theta, this%w2, this%w3, this%s1, this%s2, this%x1, this%x2, this%y)
      !$acc enter data create(this%mvv, this%mvv1, this%moo, this%moo1, this%moot, this%lvec)
   end subroutine allocate_work

   subroutine fail(this, source)
      !! Copies a refusal into `this%status` and frees what `init` had allocated so far
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(error_t), intent(in) :: source

      character(len=:), allocatable :: text

      text = source%get_message()
      call this%release()
      call this%status%set(ERROR_VALIDATION, text)
   end subroutine fail

   subroutine update_release(this)
      !! Frees host and device copies of everything; safe on an empty type
      class(trc_rccsd_lambda_update_t), intent(inout) :: this

      if (allocated(this%theta)) then
         !$acc exit data delete(this%theta, this%w2, this%w3, this%s1, this%s2, this%x1, this%x2, this%y)
         !$acc exit data delete(this%mvv, this%mvv1, this%moo, this%moo1, this%moot, this%lvec)
         deallocate (this%theta, this%w2, this%w3, this%s1, this%s2, this%x1, this%x2, this%y)
         deallocate (this%mvv, this%mvv1, this%moo, this%moo1, this%moot, this%lvec)
      end if
      call this%inter%release()
      call this%scr%release()
      call this%im%release()
      call this%lad%release()
      call this%ls%release()
      this%no = 0
      this%nv = 0
      this%naux = 0
      this%stamp = 0
      this%bytes = 0
   end subroutine update_release

   function step_ready(this, ints, amps, l, lnew) result(ok)
      !! True when `init` ran for these integrals and every amplitude set has their sizes
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l, lnew
      logical :: ok

      ok = .false.
      if (this%stamp == 0 .or. .not. allocated(this%theta)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: step before init")
         return
      end if
      if (ints%stamp /= this%stamp .or. .not. allocated(ints%b_vv)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: the integrals were rebuilt after init; re-init")
         return
      end if
      if (.not. allocated(amps%t2) .or. .not. allocated(l%t2) .or. .not. allocated(lnew%t2)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: the amplitudes are not initialised")
         return
      end if
      if (amps%no /= this%no .or. amps%nv /= this%nv .or. l%no /= this%no .or. l%nv /= this%nv .or. &
          lnew%no /= this%no .or. lnew%nv /= this%nv) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: (no, nv) of the amplitudes differ from init's")
         return
      end if
      ok = .true.
   end function step_ready

   subroutine update_step(this, ints, amps, l, lnew, la, stab)
      !! (lnew%t1, lnew%t2) = update(l%t1, l%t2): `l` holds l1(i,a) in t1 and l2(i,j,a,b) in t2, on the device, with
      !! l2 pair-symmetric (module header). `amps` and `ints` are those of `init`. `lnew` must not be `l`; it is left
      !! on the device. A refusal is in `status` and `lnew` is not to be used. `stab` has the rows
      !! RCCSD_LUPD_STEP_NAME first.
      class(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(trc_rccsd_amps_t), intent(inout) :: lnew
      type(trc_linalg_t), intent(in) :: la
      type(step_table_t), intent(inout), optional :: stab

      type(timer_type) :: tc
      type(error_t) :: err
      integer(int64) :: flops

      call this%status%clear()
      if (.not. step_ready(this, ints, amps, l, lnew)) return

      call rccsd_steps_open(stab, la, RCCSD_LUPD_ST_ONE, tc)
      call one_particle(this, amps, l, la, flops)
      call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_ONE, tc, flops)

      call rccsd_steps_open(stab, la, RCCSD_LUPD_ST_LADDER, tc)
      call rccsd_ladder_transposed(ints, l%t2, amps, this%lad, this%ls, la, &
                                   tau_symmetric=this%ls%form == RCCSD_LADDER_HALF)
      if (this%ls%status%has_error()) then
         call this%status%set(ERROR_VALIDATION, this%ls%status%get_message())
         call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_LADDER, tc)  ! a refusal still closes the NVTX range
         return
      end if
      call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_LADDER, tc, this%ls%flops_pp + this%ls%flops_dressing)

      ! the l1 terms that read nothing of g, so that lnew%t1 is set before m3 is added
      call rccsd_steps_open(stab, la, RCCSD_LUPD_ST_L1, tc)
      call l1_first(this, ints, amps, l, lnew, la)
      call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_L1, tc, l1_first_flops(this%no, this%nv))

      call rccsd_steps_open(stab, la, RCCSD_LUPD_ST_OOVV, tc)
      call l2_oovv(this, ints, amps, l, lnew, la, flops, err)
      call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_OOVV, tc, flops)

      call rccsd_steps_open(stab, la, RCCSD_LUPD_ST_RING, tc)
      call l2_ring(this, l, la, flops, err)
      call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_RING, tc, flops)

      call rccsd_steps_open(stab, la, RCCSD_LUPD_ST_L2, tc)
      call l2_finish(this, ints, l, lnew, la, flops, err)
      call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_L2, tc, flops)

      call rccsd_steps_open(stab, la, RCCSD_LUPD_ST_Q, tc)
      call l1_q(this, ints, l, lnew, la, flops, err)
      call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_Q, tc, flops)

      call rccsd_steps_open(stab, la, RCCSD_LUPD_ST_A, tc)
      call l1_a(this, ints, amps, l, lnew, la, flops, err)
      call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_A, tc, flops)

      call rccsd_steps_open(stab, la, RCCSD_LUPD_ST_L1, tc)
      call l1_last(this, ints, amps, l, lnew, la, flops, err)
      call rccsd_steps_close(stab, la, RCCSD_LUPD_ST_L1, tc, flops)

      if (err%has_error()) call this%status%set(ERROR_VALIDATION, "trc_rccsd_lambda_update: "//err%get_message())
   end subroutine update_step

   subroutine one_particle(this, amps, l, la, flops)
      !! mvv, mvv1, moo, moo1 of l (module header), and mvv1's pieces of the factors
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(trc_linalg_t), intent(in) :: la
      integer(int64), intent(out) :: flops

      integer :: no, nv

      no = this%no
      nv = this%nv
      ! mvv(b,a): theta and l2 as (klc, .) matrices
      call la%gemm("T", "N", nv, nv, no*no*nv, 1.0_dp, this%theta, no*no*nv, l%t2, no*no*nv, 0.0_dp, this%mvv, nv)
      ! moo(i,j) = sum_kdc l2(ikdc) theta(jkdc), the pair-swapped form
      call la%gemm("N", "T", no, no, no*nv*nv, 1.0_dp, l%t2, no, this%theta, no, 0.0_dp, this%moo, no)
      call la%gemm("T", "N", nv, nv, no, 1.0_dp, amps%t1, no, l%t1, no, 0.0_dp, this%mvv1, nv)
      call kernel_axpby(nv*nv, 1.0_dp, this%mvv, 1.0_dp, this%mvv1)
      call la%gemm("N", "T", no, no, nv, 1.0_dp, l%t1, no, amps%t1, no, 0.0_dp, this%moo1, no)
      call kernel_axpby(no*no, 1.0_dp, this%moo, 1.0_dp, this%moo1)
      flops = rccsd_steps_gemm_flops(nv, nv, no*no*nv) + rccsd_steps_gemm_flops(no, no, no*nv*nv) + &
              rccsd_steps_gemm_flops(nv, nv, no) + rccsd_steps_gemm_flops(no, no, nv)
   end subroutine one_particle

   pure function l1_first_flops(no, nv) result(flops)
      !! The GEMMs of `l1_first` (its GEMVs and kernels are not counted)
      integer, intent(in) :: no, nv
      integer(int64) :: flops

      flops = rccsd_steps_gemm_flops(no, nv, nv, 3) + rccsd_steps_gemm_flops(no, nv, no, 3)
   end function l1_first_flops

   subroutine l1_first(this, ints, amps, l, lnew, la)
      !! Sets lnew%t1 to the l1 terms that need neither g nor an o^2 v^2 work array
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(trc_rccsd_amps_t), intent(inout) :: lnew
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, nov, no2

      no = this%no
      nv = this%nv
      nov = no*nv
      no2 = no*no
      ! lvec(j,b) = t1 + sum_kc l1(kc) theta(kjcb) - sum_d t1(jd) mvv1(b,d) - sum_k moo(k,j) t1(k,b)
      call kernel_axpby(nov, 1.0_dp, amps%t1, 0.0_dp, this%lvec)
      call kernel_l1_theta(no, nv, this%theta, l%t1, this%lvec)
      call la%gemm("N", "T", no, nv, nv, -1.0_dp, amps%t1, no, this%mvv1, nv, 1.0_dp, this%lvec, no)
      call la%gemm("T", "N", no, nv, no, -1.0_dp, this%moo, no, amps%t1, no, 1.0_dp, this%lvec, no)
      call la%gemv("N", nov, nov, 1.0_dp, this%scr%lov, nov, this%lvec, 0.0_dp, lnew%t1)
      ! - moo v4 - v4 mvv + l1 v1 - v2 l1
      call la%gemm("N", "N", no, nv, no, -1.0_dp, this%moo, no, this%inter%fov, no, 1.0_dp, lnew%t1, no)
      call la%gemm("N", "N", no, nv, nv, -1.0_dp, this%inter%fov, no, this%mvv, nv, 1.0_dp, lnew%t1, no)
      call la%gemm("N", "N", no, nv, nv, 1.0_dp, l%t1, no, this%inter%lvv, nv, 1.0_dp, lnew%t1, no)
      call la%gemm("N", "N", no, nv, no, -1.0_dp, this%inter%loo, no, l%t1, no, 1.0_dp, lnew%t1, no)
      ! 2 (ia|jb) l1(jb) - (ij|ba) l1(jb), and w3 (2 l2 - l2 swapped)
      call kernel_l1_two(no, nv, ints%ovov, ints%oovv, l%t1, l%t2, this%im%w3, lnew%t1)
      ! - sum_jk lam(j,k,i,a) moo1(k,j)
      call kernel_transpose(no, this%moo1, this%moot)
      call la%gemv("T", no2, nov, -1.0_dp, this%scr%lam, no2, this%moot, 1.0_dp, lnew%t1)
      ! 2 sum_P b_ov(ia,P) u(P), u(P) = sum_bc b_vv(bc,P) mvv1(b,c)
      call la%gemv("T", nv*nv, this%naux, 1.0_dp, ints%b_vv, nv*nv, this%mvv1, 0.0_dp, this%scr%u)
      call la%gemv("N", nov, this%naux, 2.0_dp, ints%b_ov, nov, this%scr%u, 1.0_dp, lnew%t1)
   end subroutine l1_first

   subroutine l2_oovv(this, ints, amps, l, lnew, la, flops, err)
      !! Turns the ladder in lad%h into m3, adds m3's l1 term to lnew%t1, then adds the v1, v2 and (jk|ia) l1
      !! terms of g to the same array
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(trc_rccsd_amps_t), intent(inout) :: lnew
      type(trc_linalg_t), intent(in) :: la
      integer(int64), intent(out) :: flops
      type(error_t), intent(inout) :: err

      integer :: no, nv, no2, nv2, o2v

      no = this%no
      nv = this%nv
      no2 = no*no
      nv2 = nv*nv
      o2v = no*no*nv
      ! m3 = (R + woooo l2)/2 + (ka|lb) X / 2, X(ij,kl) = sum_cd l2(ijcd) t2(klcd)
      call la%gemm("N", "N", no2, nv2, no2, 0.5_dp, this%im%woooo, no2, l%t2, no2, 0.5_dp, this%lad%h, no2)
      call la%gemm("N", "T", no2, no2, nv2, 1.0_dp, l%t2, no2, amps%t2, no2, 0.0_dp, this%x1, no2)
      call contract_permute("kalb->klab", [no, nv, no, nv], ints%ovov, this%scr%pack1, err)
      call la%gemm("N", "N", no2, nv2, no2, 0.5_dp, this%x1, no2, this%scr%pack1, no2, 1.0_dp, this%lad%h, no2)
      call kernel_l1_m3(no, nv, this%lad%h, amps%t1, lnew%t1)
      ! + sum_c l2(ijac) v1(c,b)
      call la%gemm("N", "N", o2v, nv, nv, 1.0_dp, l%t2, o2v, this%inter%lvv, nv, 1.0_dp, this%lad%h, o2v)
      ! - sum_k v2(i,k) l2(kjab), the pair-swapped image of -l2(ikab) v2(j,k)
      call la%gemm("N", "N", no, no*nv2, no, -1.0_dp, this%inter%loo, no, l%t2, no, 1.0_dp, this%lad%h, no)
      ! - sum_k (jk|ia) l1(kb), the pair-swapped image of -l1(ka) (jb|ik)
      call contract_permute("jkia->ijak", [no, no, no, nv], ints%ooov, this%s1, err)
      call la%gemm("N", "N", o2v, nv, no, -1.0_dp, this%s1, o2v, l%t1, no, 1.0_dp, this%lad%h, o2v)
      flops = rccsd_steps_gemm_flops(no2, nv2, no2, 2) + rccsd_steps_gemm_flops(no2, no2, nv2) + &
              rccsd_steps_gemm_flops(o2v, nv, nv) + rccsd_steps_gemm_flops(no, no*nv2, no) + &
              rccsd_steps_gemm_flops(o2v, nv, no)
   end subroutine l2_oovv

   subroutine l2_ring(this, l, la, flops, err)
      !! g(ijab) += R(jb,ia)/2 + X(jb,ia)/2 + X(ib,ja), R = wovvo (2 l2 - l2 swapped), X = wovvo_x l2
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_amps_t), intent(in) :: l
      type(trc_linalg_t), intent(in) :: la
      integer(int64), intent(out) :: flops
      type(error_t), intent(inout) :: err

      integer :: no, nv, nov
      integer :: oovv(4), ovov(4)

      no = this%no
      nv = this%nv
      nov = no*nv
      oovv = [no, no, nv, nv]
      ovov = [no, nv, no, nv]
      ! w1(k,c,i,a) = 2 l2(ikac) - l2(ikca)
      call contract_permute("ikca->kcia", oovv, l%t2, this%scr%pack1, err, -1.0_dp)
      call contract_permute("ikac->kcia", oovv, l%t2, this%scr%pack1, err, 2.0_dp, 1.0_dp)
      call la%gemm("N", "N", nov, nov, nov, 1.0_dp, this%im%wovvo, nov, this%scr%pack1, nov, 0.0_dp, this%w2, nov)
      ! w1(k,c,i,a) = l2(ikca)
      call contract_permute("ikca->kcia", oovv, l%t2, this%scr%pack1, err)
      call la%gemm("N", "N", nov, nov, nov, 1.0_dp, this%im%wovvo_x, nov, this%scr%pack1, nov, 0.0_dp, this%w3, nov)
      call contract_permute("jbia->ijab", ovov, this%w2, this%lad%h, err, 0.5_dp, 1.0_dp)
      call contract_permute("jbia->ijab", ovov, this%w3, this%lad%h, err, 0.5_dp, 1.0_dp)
      call contract_permute("ibja->ijab", ovov, this%w3, this%lad%h, err, 1.0_dp, 1.0_dp)
      flops = rccsd_steps_gemm_flops(nov, nov, nov, 2)
   end subroutine l2_ring

   subroutine l2_finish(this, ints, l, lnew, la, flops, err)
      !! The remaining terms of g, then lnew%t2(ijab) = (g(ijab) + g(jiba))/D2
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: l
      type(trc_rccsd_amps_t), intent(inout) :: lnew
      type(trc_linalg_t), intent(in) :: la
      integer(int64), intent(out) :: flops
      type(error_t), intent(inout) :: err

      integer :: no, nv, nov, ovo, naux

      no = this%no
      nv = this%nv
      naux = this%naux
      nov = no*nv
      ovo = no*nv*no
      ! U(i,a,j,b) = (ia|jb)/2 + l1(ia) v4(jb) - (ia|jc) mvv1(c,b) - moo1(i,k) (ka|jb) + (ia|cb) l1(jc)
      call kernel_uinit(no, nv, ints%ovov, l%t1, this%inter%fov, this%w2)
      call la%gemm("N", "N", ovo, nv, nv, -1.0_dp, ints%ovov, ovo, this%mvv1, nv, 1.0_dp, this%w2, ovo)
      call la%gemm("N", "N", no, nv*nov, no, -1.0_dp, this%moo1, no, ints%ovov, no, 1.0_dp, this%w2, no)
      ! scr%n1p(j,b;P) = sum_c l1(jc) b_vv(cb,P)
      call la%gemm("N", "N", no, nv*naux, nv, 1.0_dp, l%t1, no, ints%b_vv, nv, 0.0_dp, this%scr%n1p, no)
      call la%gemm("N", "T", nov, nov, naux, 1.0_dp, ints%b_ov, nov, this%scr%n1p, nov, 1.0_dp, this%w2, nov)
      call contract_permute("iajb->ijab", [no, nv, no, nv], this%w2, this%lad%h, err, 1.0_dp, 1.0_dp)
      call kernel_l2_final(no, nv, ints%eps_o, ints%eps_v, this%lad%h, lnew%t2)
      flops = rccsd_steps_gemm_flops(ovo, nv, nv) + rccsd_steps_gemm_flops(no, nv*nov, no) + &
              rccsd_steps_gemm_flops(no, nv*naux, nv) + rccsd_steps_gemm_flops(nov, nov, naux)
   end subroutine l2_finish

   subroutine l1_q(this, ints, l, lnew, la, flops, err)
      !! N into w3 as N'(i,b,k,d) = N(i,k,b,d), then Q(i,b;P) into scr%n1p and lnew%t1 -= sum_bP Q(i,b;P) b_vv(ab,P).
      !! Overwrites w1 = scr%pack1 and w2; `l1_a` reads w3.
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: l
      type(trc_rccsd_amps_t), intent(inout) :: lnew
      type(trc_linalg_t), intent(in) :: la
      integer(int64), intent(out) :: flops
      type(error_t), intent(inout) :: err

      type(trc_gemm_batch_t) :: shape
      integer :: no, nv, nov, naux
      integer :: oovv(4)

      no = this%no
      nv = this%nv
      naux = this%naux
      nov = no*nv
      oovv = [no, no, nv, nv]
      ! Q's mvv1 part: Q(i,c;P) = sum_b b_ov(ib,P) mvv1(b,c), one item per P, mvv1 shared
      shape = trc_gemm_batch(no, nv, nv, naux, lda=no, ldb=nv, ldc=no, sa=int(nov, int64), sb=0_int64, &
                             sc=int(nov, int64))
      call la%gemm_batched(shape, ints%b_ov, this%mvv1, this%scr%n1p)
      ! N'(ib,kd) = - sum_jc [l2(ijbc) - l2(ijcb)/2] theta(jkcd)
      call contract_permute("ijbc->ibjc", oovv, l%t2, this%scr%pack1, err)
      call contract_permute("ijcb->ibjc", oovv, l%t2, this%scr%pack1, err, -0.5_dp, 1.0_dp)
      call contract_permute("jkcd->jckd", oovv, this%theta, this%w2, err)
      call la%gemm("N", "N", nov, nov, nov, -1.0_dp, this%scr%pack1, nov, this%w2, nov, 0.0_dp, this%w3, nov)
      ! the amplitude part: sum_kd [2 N' (ib,kd) - 2 l2(ikbd) + l2(ikdb)] b_ov(kd,P)
      call kernel_axpby(nov*nov, 2.0_dp, this%w3, 0.0_dp, this%w2)
      call contract_permute("ikbd->ibkd", oovv, l%t2, this%w2, err, -2.0_dp, 1.0_dp)
      call contract_permute("ikdb->ibkd", oovv, l%t2, this%w2, err, 1.0_dp, 1.0_dp)
      call la%gemm("N", "N", nov, naux, nov, 1.0_dp, this%w2, nov, ints%b_ov, nov, 1.0_dp, this%scr%n1p, nov)
      ! b_vv(a, (c,P)) has leading dimension v
      call la%gemm("N", "T", no, nv, nv*naux, -1.0_dp, this%scr%n1p, no, ints%b_vv, nv, 1.0_dp, lnew%t1, no)
      flops = rccsd_steps_gemm_flops(shape%m, shape%n, shape%k, shape%count) + rccsd_steps_gemm_flops(nov, nov, nov) + &
              rccsd_steps_gemm_flops(nov, naux, nov) + rccsd_steps_gemm_flops(no, nv, nv*naux)
   end subroutine l1_q

   subroutine l1_a(this, ints, amps, l, lnew, la, flops, err)
      !! From N' in w3: A(i,k,c,d), Y(i,k;P) = sum_cd A(ikcd) b_vv(cd,P), lnew%t1 -= sum_kP Y(i,k;P) b_ov(ka,P)
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(trc_rccsd_amps_t), intent(inout) :: lnew
      type(trc_linalg_t), intent(in) :: la
      integer(int64), intent(out) :: flops
      type(error_t), intent(inout) :: err

      integer :: no, nv, nov, no2, nv2, naux

      no = this%no
      nv = this%nv
      naux = this%naux
      nov = no*nv
      no2 = no*no
      nv2 = nv*nv
      ! A'(i,c,k,d) = 3/2 sum_jb l2(ijbc) t2(kjbd) - N'(ic,kd)
      call contract_permute("ijbc->icjb", [no, no, nv, nv], l%t2, this%scr%pack1, err)
      call contract_permute("kjbd->jbkd", [no, no, nv, nv], amps%t2, this%w2, err)
      call la%gemm("N", "N", nov, nov, nov, 1.5_dp, this%scr%pack1, nov, this%w2, nov, -1.0_dp, this%w3, nov)
      call contract_permute("ickd->ikcd", [no, nv, no, nv], this%w3, this%scr%pack1, err)
      call la%gemm("N", "N", no2, naux, nv2, 1.0_dp, this%scr%pack1, no2, ints%b_vv, nv2, 0.0_dp, this%y, no2)
      ! Y(i,(k,P)) against scr%bovp, b_ov as ((k,P), a)
      call la%gemm("N", "N", no, nv, no*naux, -1.0_dp, this%y, no, this%scr%bovp, no*naux, 1.0_dp, lnew%t1, no)
      flops = rccsd_steps_gemm_flops(nov, nov, nov) + rccsd_steps_gemm_flops(no2, naux, nv2) + &
              rccsd_steps_gemm_flops(no, nv, no*naux)
   end subroutine l1_a

   subroutine l1_last(this, ints, amps, l, lnew, la, flops, err)
      !! The wvvvo terms without an (ia|bc), woovo, and the division by D1
      type(trc_rccsd_lambda_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(trc_rccsd_amps_t), intent(inout) :: lnew
      type(trc_linalg_t), intent(in) :: la
      integer(int64), intent(out) :: flops
      type(error_t), intent(inout) :: err

      integer :: no, nv, nov, no2, o2v, o3, nv2
      integer :: ooov(4)

      no = this%no
      nv = this%nv
      nov = no*nv
      no2 = no*no
      nv2 = nv*nv
      o2v = no*no*nv
      o3 = no*no*no
      ooov = [no, no, no, nv]
      ! s1(i,j,b,m) = sum_c l2(ijbc) t1(mc)
      call la%gemm("N", "T", o2v, no, nv, 1.0_dp, l%t2, o2v, amps%t1, no, 0.0_dp, this%s1, o2v)
      ! - sum_jmc s1(j,i,c,m) w_both(m,a,j,c)
      call contract_permute("jicm->ijmc", [no, no, nv, no], this%s1, this%s2, err)
      call contract_permute("majc->jmca", [no, nv, no, nv], this%im%w_both, this%w2, err)
      call la%gemm("N", "N", no, nv, o2v, -1.0_dp, this%s2, no, this%w2, o2v, 1.0_dp, lnew%t1, no)
      ! - sum_jmc s1(i,j,c,m) (w_diff - w_both)(m,a,j,c)
      call contract_permute("majc->jmca", [no, nv, no, nv], this%im%w_both, this%w2, err, -1.0_dp)
      call contract_permute("majc->jmca", [no, nv, no, nv], this%im%w_diff, this%w2, err, 1.0_dp, 1.0_dp)
      call contract_permute("ijcm->ijmc", [no, no, nv, no], this%s1, this%s2, err)
      call la%gemm("N", "N", no, nv, o2v, -1.0_dp, this%s2, no, this%w2, o2v, 1.0_dp, lnew%t1, no)
      ! + sum_pjl X(ij,lp) lam(p,j,l,a), X(ij,kl) = sum_cd l2(ijcd) tau(klcd)
      call la%gemm("N", "T", no2, no2, nv2, 1.0_dp, l%t2, no2, this%inter%tau, no2, 0.0_dp, this%x1, no2)
      call contract_permute("ijlp->ipjl", [no, no, no, no], this%x1, this%x2, err)
      call la%gemm("N", "N", no, nv, o3, 1.0_dp, this%x2, no, this%scr%lam, o3, 1.0_dp, lnew%t1, no)
      ! - sum_jkc woovo(i,j,k,c) l2(kjca)
      call contract_permute("ijkc->ikjc", ooov, this%im%woovo, this%s2, err)
      call la%gemm("N", "N", no, nv, o2v, -1.0_dp, this%s2, no, l%t2, o2v, 1.0_dp, lnew%t1, no)
      call kernel_denominator_l1(no, nv, ints%eps_o, ints%eps_v, lnew%t1)
      flops = rccsd_steps_gemm_flops(o2v, no, nv) + rccsd_steps_gemm_flops(no, nv, o2v, 3) + &
              rccsd_steps_gemm_flops(no2, no2, nv2) + rccsd_steps_gemm_flops(no, nv, o3)
   end subroutine l1_last

   subroutine kernel_axpby(n, a, x, b, y)
      !! ---------------------------------------------------------------------
      !! the kernels. Each takes its arrays as explicit-shape dummies and runs one
      !! `do concurrent` over ONE flattened index. The pure reorderings are
      !! trc_contract's `contract_permute`.
      !! ---------------------------------------------------------------------
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

   subroutine kernel_transpose(n, a, at)
      !! at(j,i) = a(i,j)
      integer, intent(in) :: n
      real(dp), intent(in) :: a(n, n)
      real(dp), intent(out) :: at(n*n)
      integer :: t

      do concurrent(t=1:n*n)
         at(t) = a((t - 1)/n + 1, mod(t - 1, n) + 1)
      end do
   end subroutine kernel_transpose

   subroutine kernel_l1_theta(no, nv, theta, l1, lvec)
      !! lvec(j,b) += sum_kc theta(k,j,c,b) l1(k,c)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: theta(no, no, nv, nv), l1(no, nv)
      real(dp), intent(inout) :: lvec(no*nv)
      integer :: t

      do concurrent(t=1:no*nv)
         lvec(t) = lvec(t) + l1_theta_item(TRC_DC_INDEX(t), no, nv, theta, l1)
      end do
   end subroutine kernel_l1_theta

   subroutine kernel_l1_two(no, nv, ovov, oovv, l1, l2, w3, r)
      !! r(i,a) += sum_jb [2 (ia|jb) - (ij|ba)] l1(jb) + [2 l2(ijab) - l2(ijba)] w3(jb)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv), oovv(no, no, nv, nv)
      real(dp), intent(in) :: l1(no, nv), l2(no, no, nv, nv), w3(no, nv)
      real(dp), intent(inout) :: r(no*nv)
      integer :: t

      do concurrent(t=1:no*nv)
         r(t) = r(t) + l1_two_item(TRC_DC_INDEX(t), no, nv, ovov, oovv, l1, l2, w3)
      end do
   end subroutine kernel_l1_two

   subroutine kernel_l1_m3(no, nv, m, t1, r)
      !! r(i,a) += sum_jb [2 (m(ijab) + m(jiba)) - (m(ijba) + m(jiab))] t1(jb)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: m(no, no, nv, nv), t1(no, nv)
      real(dp), intent(inout) :: r(no*nv)
      integer :: t

      do concurrent(t=1:no*nv)
         r(t) = r(t) + l1_m3_item(TRC_DC_INDEX(t), no, nv, m, t1)
      end do
   end subroutine kernel_l1_m3

   subroutine kernel_uinit(no, nv, ovov, l1, v4, u)
      !! u(i,a,j,b) = (ia|jb)/2 + l1(ia) v4(jb)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no*nv, no*nv), l1(no*nv), v4(no*nv)
      real(dp), intent(out) :: u(no*nv*no*nv)
      integer :: t, nov

      nov = no*nv
      do concurrent(t=1:nov*nov)
         u(t) = 0.5_dp*ovov(mod(t - 1, nov) + 1, (t - 1)/nov + 1) + l1(mod(t - 1, nov) + 1)*v4((t - 1)/nov + 1)
      end do
   end subroutine kernel_uinit

   subroutine kernel_l2_final(no, nv, eps_o, eps_v, g, l2)
      !! l2(ijab) = (g(ijab) + g(jiba)) / (e_i - e_a + e_j - e_b)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: eps_o(no), eps_v(nv), g(no, no, nv, nv)
      real(dp), intent(out) :: l2(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         l2(t) = l2_final_item(TRC_DC_INDEX(t), no, nv, eps_o, eps_v, g)
      end do
   end subroutine kernel_l2_final

   subroutine kernel_denominator_l1(no, nv, eps_o, eps_v, l1)
      !! l1(i,a) /= e_i - e_a
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: eps_o(no), eps_v(nv)
      real(dp), intent(inout) :: l1(no*nv)
      integer :: t

      do concurrent(t=1:no*nv)
         l1(t) = l1(t)/(eps_o(mod(t - 1, no) + 1) - eps_v((t - 1)/no + 1))
      end do
   end subroutine kernel_denominator_l1

   pure function l1_theta_item(t, no, nv, theta, l1) result(s)
      !! sum_kc theta(k,j,c,b) l1(k,c) for t the position of (j,b)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: theta(no, no, nv, nv), l1(no, nv)
      real(dp) :: s
      integer :: j, b, k, c

      j = mod(t - 1, no) + 1
      b = (t - 1)/no + 1
      s = 0.0_dp
      do c = 1, nv
         do k = 1, no
            s = s + theta(k, j, c, b)*l1(k, c)
         end do
      end do
   end function l1_theta_item

   pure function l1_two_item(t, no, nv, ovov, oovv, l1, l2, w3) result(s)
      !! The sum of `kernel_l1_two` for t the position of (i,a)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv), oovv(no, no, nv, nv)
      real(dp), intent(in) :: l1(no, nv), l2(no, no, nv, nv), w3(no, nv)
      real(dp) :: s
      integer :: i, a, j, b

      i = mod(t - 1, no) + 1
      a = (t - 1)/no + 1
      s = 0.0_dp
      do b = 1, nv
         do j = 1, no
            s = s + (2.0_dp*ovov(i, a, j, b) - oovv(i, j, b, a))*l1(j, b) &
                + (2.0_dp*l2(i, j, a, b) - l2(i, j, b, a))*w3(j, b)
         end do
      end do
   end function l1_two_item

   pure function l1_m3_item(t, no, nv, m, t1) result(s)
      !! The sum of `kernel_l1_m3` for t the position of (i,a)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: m(no, no, nv, nv), t1(no, nv)
      real(dp) :: s
      integer :: i, a, j, b

      i = mod(t - 1, no) + 1
      a = (t - 1)/no + 1
      s = 0.0_dp
      do b = 1, nv
         do j = 1, no
            s = s + (2.0_dp*(m(i, j, a, b) + m(j, i, b, a)) - (m(i, j, b, a) + m(j, i, a, b)))*t1(j, b)
         end do
      end do
   end function l1_m3_item

   pure function l2_final_item(t, no, nv, eps_o, eps_v, g) result(x)
      !! (g(ijab) + g(jiba))/D2 at position t of (no, no, nv, nv); the sum of the two denominator halves is written
      !! so that (ijab) and (jiba) get the same value bit for bit
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: eps_o(no), eps_v(nv), g(no, no, nv, nv)
      real(dp) :: x
      integer :: q, i, j, a, b

      q = t - 1
      i = mod(q, no) + 1
      q = q/no
      j = mod(q, no) + 1
      q = q/no
      a = mod(q, nv) + 1
      b = q/nv + 1
      x = (g(i, j, a, b) + g(j, i, b, a))/((eps_o(i) - eps_v(a)) + (eps_o(j) - eps_v(b)))
   end function l2_final_item

end module trc_rccsd_lambda_update
