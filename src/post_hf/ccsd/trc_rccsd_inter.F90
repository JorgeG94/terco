!
! The one-particle intermediates and the energy functional of a closed-shell
! RI-CCSD(T) calculation: tau, Foo, Fvv, Fov, Loo, Lvv, E(t1, t2) and the MP2
! starting amplitudes.
!
! Hirata, Podeszwa, Tobita and Bartlett, JCP 120, 2581 (2004), Eqs. (37)-(41),
! as pyscf/cc/rintermediates.py and rccsd.energy write them, ported from the
! same author's mqc_czt_rcc.f90 (build_tau, cc_foo, cc_fvv, cc_fov, cc_loo,
! cc_lvv, rccsd_correlation_energy, mp2_guess). Chemists' notation, spatial
! orbitals, REAL orbitals (so (ia|jb) = (jb|ia) and (ij|ka) = (ka|ij)).
!
! LAYOUT, as in mqc (first index fastest):
!
!     t1(i,a)   t2(i,j,a,b)   tau(i,j,a,b)   Foo(k,i)   Fvv(a,c)   Fov(k,c)
!     Loo(k,i)  Lvv(a,c)
!
! which is also PySCF's index order for each of them.
!
! CANONICAL ORBITALS are assumed: the Fock matrix is diagonal, f_ov = 0. The
! bare `foo`/`fvv` pyscf adds to Foo, Fvv (and so to Loo, Lvv) is therefore
! NOT in these arrays, exactly as in mqc: it is diag(eps_o), diag(eps_v), and
! update_amps subtracts the same diagonal straight back out. To compare with
! pyscf's cc_Foo, add eps_o on the diagonal. The `f_ov t1` terms of Loo, Lvv
! and the 2 f_ov t1 term of the energy are zero for the same reason. Leaving the
! diagonal out is exact only while the amplitude denominators are exactly the Fock
! diagonal (e_i - e_a, ...): a level shift would have to be added to the Fvv, Lvv
! diagonal (and the Foo, Loo one) here.
!
! The linear-algebra handle `la` is taken intent(in): the caller initialises it
! first (as `trc_rccsd_ints_t%build` does).
!
! How each becomes GEMMs (the work is O(o^3 v^2) for Foo, Fvv and O(o^2 v^2)
! or O(o v^2 naux) for the rest; only Foo and Fvv are O(N^5)):
!
!     L(k,c,l,d) = 2 (kc|ld) - (kd|lc)   built once (`init` of the scratch),
!                  it is the combination of Foo, Fvv, Fov and the energy
!     Foo  = L_(k;cld) . taup_(i;cld)^T         taup(i,c,l,d) = tau(i,l,c,d)
!     Fvv  = - taur_(kld;a)^T . lq_(kld;c)      taur(k,l,d,a) = tau(k,l,a,d)
!                                               lq(k,l,d,c)   = L(k,c,l,d)
!     Fov  = L_(kc;ld) . t1_(ld)                (a matrix-vector product)
!     Loo  = Foo + lam_(ki;lc) . t1_(lc)        lam = 2 (ki|lc) - (li|kc)
!     Lvv  = Fvv + sum_kd [2 (kd|ac) - (kc|ad)] t1(k,d), from the fitted factors
!            ((kd|ac) = sum_P b_ov(kd,P) b_vv(ac,P); no (ia|bc) block exists):
!              N1(k,a;P) = sum_d t1(k,d) b_vv(d,a;P)    ONE GEMM (o, v naux, v), kept
!                                                       in `n1` for the T2 equation
!              direct    2 b_vv(ac;P) u(P), u = b_ov^T t1   two matrix-vector products
!              exchange  - sum_(k,P) N1(k,a;P) b_ov(k,c;P)   ONE GEMM (v, v, o naux) of
!                        N1 and b_ov permuted to (k,P,.) (`n1p`, and `bovp` packed once);
!                        2 o v^2 naux flops, o times fewer than one GEMM into an
!                        o^2 v^2 R((k,c),(k',a)) followed by its trace over k = k'
!
! RESIDENCY. Every array of the three types here exists on the host AND on the
! device (`enter data` in `init`, `exit data` in `release`, on the same
! `x%name` expression). The device copy is the one the kernels use; `push`
! sends the host copy of the amplitudes, `pull` brings results back. Each
! component is mapped on its own and the parent type is never mapped, so a
! kernel gets an array only as an explicit-shape dummy argument of a
! `kernel_*` routine below and never names `x%name` in a `do concurrent`
! body (an illegal address under mem:separate that the host build hides).
! `trc_rccsd_ints_t` follows the same rule.
!
! DETERMINISM. The energy is summed from one partial per (j,b) in ascending
! order of that pair, on the host, so it does not depend on the thread count
! or on host against device. The GEMMs are the library's.
!
! SCRATCH, allocated once by `trc_rccsd_scratch_t%init` (o = active occupied,
! v = virtual, 8 bytes each, device and host):
!
!     lov      o^2 v^2    L(k,c,l,d), kept for the life of the integrals
!     pack1    o^2 v^2    taup (Foo), then taur (Fvv): the only per-call scratch
!     lq       o^2 v^2    lq(k,l,d,c), kept for the life of the integrals (Fvv); not general scratch
!     lam      o^3 v      Loo's combination, kept for the life of the integrals
!     part_s, part_d  o v   energy partials
!     n1       o v naux   N1(k,a;P), written by Lvv and read by the T2 equation
!     n1p      o v naux   N1 as (k,P,a) for Lvv; per call, the T1 equation borrows it
!     bovp     o v naux   b_ov as (k,P,c), kept for the life of the integrals
!     u        naux
!
! Every block is addressed with a default integer: a block above 2^31-1
! elements is refused by `init` (status), not wrapped.
!
! STALE PACKS. `init` packs lov, lq and lam from the integrals once and records
! `ints%stamp`. After `ints%build` is called again (new geometry or orbitals,
! even with the same no and nv) the stamp differs and every routine that takes
! `ints` refuses with a message: re-`init` the scratch. `rccsd_inter_energy`
! takes no `ints`, so it cannot check; its caller keeps the scratch and the
! integrals in step.
!
! Failures are reported in `scratch%status` (a routine whose `ints`, `amps`, `inter`
! or `scr` differ in (no, nv), or are not built, sets it and returns without touching
! anything); none stops.
!
! RANKS. No MPI call here: every rank that holds the integrals can run these routines on its own device; under
! several ranks only rank 0 iterates today (trc_rccsd_driver).
module trc_rccsd_inter
   use pic_types, only: dp, int32, int64
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_steps, only: rccsd_steps_gemm_flops
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   character(len=*), parameter :: STALE_MSG = "trc_rccsd_inter: the integrals were rebuilt after the scratch was "// &
                                              "initialised; re-init the scratch after ints%build"

   public :: trc_rccsd_amps_t, trc_rccsd_inter_t, trc_rccsd_scratch_t
   public :: rccsd_inter_build_tau, rccsd_inter_foo, rccsd_inter_fvv, rccsd_inter_fov, rccsd_inter_loo, rccsd_inter_lvv
   public :: rccsd_inter_one_particle, rccsd_inter_energy, rccsd_inter_mp2_guess, rccsd_inter_flops

   type :: trc_rccsd_amps_t
      !! The amplitudes t1(i,a), t2(i,j,a,b) on host and device
      integer :: no = 0  !! active occupied orbitals
      integer :: nv = 0  !! virtual orbitals
      real(dp), allocatable :: t1(:, :)        !! t1(i,a)
      real(dp), allocatable :: t2(:, :, :, :)  !! t2(i,j,a,b)
   contains
      procedure :: init => amps_init
      procedure :: release => amps_release
      procedure :: push => amps_push
      procedure :: pull => amps_pull
   end type trc_rccsd_amps_t

   type :: trc_rccsd_inter_t
      !! tau and the one-particle intermediates on host and device, valid after the matching routine
      integer :: no = 0
      integer :: nv = 0
      real(dp), allocatable :: tau(:, :, :, :)  !! t2(i,j,a,b) + t1(i,a) t1(j,b)
      real(dp), allocatable :: foo(:, :)        !! Foo(k,i), without the diagonal eps_o
      real(dp), allocatable :: fvv(:, :)        !! Fvv(a,c), without the diagonal eps_v
      real(dp), allocatable :: fov(:, :)        !! Fov(k,c)
      real(dp), allocatable :: loo(:, :)        !! Loo(k,i), without the diagonal eps_o
      real(dp), allocatable :: lvv(:, :)        !! Lvv(a,c), without the diagonal eps_v
   contains
      procedure :: init => inter_init
      procedure :: release => inter_release
      procedure :: pull => inter_pull
   end type trc_rccsd_inter_t

   type :: trc_rccsd_scratch_t
      !! The work arrays of the routines below, sizes in the module header
      type(error_t) :: status  !! set by `init` and by every routine that refuses its arguments
      integer :: no = 0
      integer :: nv = 0
      integer :: stamp = 0  !! `ints%stamp` of the integrals `init` packed from
      real(dp), allocatable :: lov(:, :, :, :)  !! L(k,c,l,d) = 2 (kc|ld) - (kd|lc)
      real(dp), allocatable :: pack1(:)         !! o^2 v^2
      real(dp), allocatable :: lq(:)            !! o^2 v^2, lq(k,l,d,c) = L(k,c,l,d), built by `init`
      real(dp), allocatable :: lam(:)           !! o^3 v, lam(k,i,l,c) = 2 (ki|lc) - (li|kc), built by `init`
      real(dp), allocatable :: part_s(:)        !! o v energy partials, singles part
      real(dp), allocatable :: part_d(:)        !! o v energy partials, doubles part
      real(dp), allocatable :: n1(:)            !! o v naux, N1(k,a;P) = sum_d t1(k,d) b_vv(d,a;P), from `rccsd_inter_lvv`
      real(dp), allocatable :: n1p(:)           !! o v naux, N1 permuted to (k,P,a) for Lvv: per-call scratch
      real(dp), allocatable :: bovp(:)          !! o v naux, b_ov permuted to (k,P,c), built by `init`
      real(dp), allocatable :: u(:)             !! naux, u(P) = sum_kd b_ov(kd,P) t1(k,d)
   contains
      procedure :: init => scratch_init
      procedure :: release => scratch_release
   end type trc_rccsd_scratch_t

contains

   ! ---------------------------------------------------------------------
   ! the types
   ! ---------------------------------------------------------------------

   subroutine amps_init(this, no, nv)
      !! Zeroed amplitudes on host and device
      class(trc_rccsd_amps_t), intent(inout) :: this
      integer, intent(in) :: no, nv

      call this%release()
      this%no = no
      this%nv = nv
      allocate (this%t1(no, nv), this%t2(no, no, nv, nv))
      this%t1 = 0.0_dp
      this%t2 = 0.0_dp
      !$acc enter data copyin(this%t1, this%t2)
   end subroutine amps_init

   subroutine amps_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_amps_t), intent(inout) :: this

      if (allocated(this%t1)) then
         !$acc exit data delete(this%t1)
         deallocate (this%t1)
      end if
      if (allocated(this%t2)) then
         !$acc exit data delete(this%t2)
         deallocate (this%t2)
      end if
      this%no = 0
      this%nv = 0
   end subroutine amps_release

   subroutine amps_push(this)
      !! Host copy to device
      class(trc_rccsd_amps_t), intent(inout) :: this

      if (.not. allocated(this%t1)) return
      !$acc update device(this%t1, this%t2)
   end subroutine amps_push

   subroutine amps_pull(this)
      !! Device copy to host
      class(trc_rccsd_amps_t), intent(inout) :: this

      if (.not. allocated(this%t1)) return
      !$acc update self(this%t1, this%t2)
   end subroutine amps_pull

   subroutine inter_init(this, no, nv)
      !! Intermediates on host and device, zeroed
      class(trc_rccsd_inter_t), intent(inout) :: this
      integer, intent(in) :: no, nv

      call this%release()
      this%no = no
      this%nv = nv
      allocate (this%tau(no, no, nv, nv), this%foo(no, no), this%fvv(nv, nv), this%fov(no, nv))
      allocate (this%loo(no, no), this%lvv(nv, nv))
      this%tau = 0.0_dp
      this%foo = 0.0_dp
      this%fvv = 0.0_dp
      this%fov = 0.0_dp
      this%loo = 0.0_dp
      this%lvv = 0.0_dp
      !$acc enter data copyin(this%tau, this%foo, this%fvv, this%fov, this%loo, this%lvv)
   end subroutine inter_init

   subroutine inter_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_inter_t), intent(inout) :: this

      if (allocated(this%tau)) then
         !$acc exit data delete(this%tau)
         deallocate (this%tau)
      end if
      if (allocated(this%foo)) then
         !$acc exit data delete(this%foo)
         deallocate (this%foo)
      end if
      if (allocated(this%fvv)) then
         !$acc exit data delete(this%fvv)
         deallocate (this%fvv)
      end if
      if (allocated(this%fov)) then
         !$acc exit data delete(this%fov)
         deallocate (this%fov)
      end if
      if (allocated(this%loo)) then
         !$acc exit data delete(this%loo)
         deallocate (this%loo)
      end if
      if (allocated(this%lvv)) then
         !$acc exit data delete(this%lvv)
         deallocate (this%lvv)
      end if
      this%no = 0
      this%nv = 0
   end subroutine inter_release

   subroutine inter_pull(this)
      !! Device copy to host
      class(trc_rccsd_inter_t), intent(inout) :: this

      if (.not. allocated(this%tau)) return
      !$acc update self(this%tau, this%foo, this%fvv, this%fov, this%loo, this%lvv)
   end subroutine inter_pull

   subroutine scratch_init(this, ints)
      !! Allocates the work arrays for these integrals and builds L(k,c,l,d). Refusals go to `status`.
      class(trc_rccsd_scratch_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints

      integer :: no, nv, naux

      call this%release()
      call this%status%clear()
      if (.not. allocated(ints%ovov) .or. .not. allocated(ints%ooov) .or. .not. allocated(ints%b_ov) .or. &
          .not. allocated(ints%b_vv)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_inter: the integral blocks are not built")
         return
      end if
      no = ints%no
      nv = ints%nv
      naux = ints%naux
      ! the kernels index o^2 v^2, o^3 v and o v naux arrays with default integers
      if (int(no, int64)*no*nv*nv > int(huge(0_int32), int64) .or. int(no, int64)*no*no*nv > int(huge(0_int32), int64) .or. &
          int(no, int64)*nv*naux > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_inter: an o^2 v^2, o^3 v or o v naux array has more than "// &
                              "2^31-1 elements (default-integer kernel index)")
         return
      end if
      this%no = no
      this%nv = nv
      this%stamp = ints%stamp
      allocate (this%lov(no, nv, no, nv), this%pack1(no*no*nv*nv), this%lq(no*no*nv*nv))
      allocate (this%lam(no*no*no*nv), this%part_s(no*nv), this%part_d(no*nv))
      allocate (this%n1(no*nv*naux), this%n1p(no*nv*naux), this%bovp(no*nv*naux), this%u(naux))
      !$acc enter data create(this%lov, this%pack1, this%lq, this%lam, this%part_s, this%part_d)
      !$acc enter data create(this%n1, this%n1p, this%bovp, this%u)
      call kernel_lov(no, nv, ints%ovov, this%lov)
      call kernel_pack_lq(no, nv, this%lov, this%lq)
      call kernel_pack_lam(no, nv, ints%ooov, this%lam)
      call kernel_pack_kpa(no, nv, naux, ints%b_ov, this%bovp)
   end subroutine scratch_init

   subroutine scratch_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_scratch_t), intent(inout) :: this

      if (allocated(this%lov)) then
         !$acc exit data delete(this%lov)
         deallocate (this%lov)
      end if
      if (allocated(this%pack1)) then
         !$acc exit data delete(this%pack1)
         deallocate (this%pack1)
      end if
      if (allocated(this%lq)) then
         !$acc exit data delete(this%lq)
         deallocate (this%lq)
      end if
      if (allocated(this%lam)) then
         !$acc exit data delete(this%lam)
         deallocate (this%lam)
      end if
      if (allocated(this%part_s)) then
         !$acc exit data delete(this%part_s)
         deallocate (this%part_s)
      end if
      if (allocated(this%part_d)) then
         !$acc exit data delete(this%part_d)
         deallocate (this%part_d)
      end if
      if (allocated(this%n1)) then
         !$acc exit data delete(this%n1)
         deallocate (this%n1)
      end if
      if (allocated(this%n1p)) then
         !$acc exit data delete(this%n1p)
         deallocate (this%n1p)
      end if
      if (allocated(this%bovp)) then
         !$acc exit data delete(this%bovp)
         deallocate (this%bovp)
      end if
      if (allocated(this%u)) then
         !$acc exit data delete(this%u)
         deallocate (this%u)
      end if
      this%no = 0
      this%nv = 0
      this%stamp = 0
   end subroutine scratch_release

   function ready_inter(ints, inter, scr) result(ok)
      !! True when ints, inter and scr are built and agree on (no, nv); else `scr%status` says so
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      logical :: ok

      ok = .false.
      if (.not. allocated(scr%lq) .or. .not. allocated(ints%ovov)) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: the scratch or the integrals are not initialised")
         return
      end if
      if (.not. allocated(inter%tau) .or. .not. allocated(inter%foo) .or. .not. allocated(inter%fvv)) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: the intermediates are not initialised")
         return
      end if
      if (inter%no /= ints%no .or. inter%nv /= ints%nv .or. scr%no /= ints%no .or. scr%nv /= ints%nv) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: (no, nv) differ between integrals, "// &
                             "intermediates and scratch")
         return
      end if
      if (scr%stamp /= ints%stamp) then
         call scr%status%set(ERROR_VALIDATION, STALE_MSG)
         return
      end if
      ok = .true.
   end function ready_inter

   function ready(ints, amps, inter, scr) result(ok)
      !! True when the four agree on (no, nv) and are built; else `scr%status` says which not
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      logical :: ok

      ok = .false.
      if (.not. allocated(scr%lov) .or. .not. allocated(scr%lq) .or. .not. allocated(scr%lam) .or. &
          .not. allocated(scr%n1) .or. .not. allocated(ints%ovov)) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: the scratch or the integrals are not initialised")
         return
      end if
      if (.not. allocated(amps%t1) .or. .not. allocated(inter%tau)) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: the amplitudes or the intermediates are not initialised")
         return
      end if
      if (amps%no /= ints%no .or. amps%nv /= ints%nv .or. inter%no /= ints%no .or. inter%nv /= ints%nv .or. &
          scr%no /= ints%no .or. scr%nv /= ints%nv) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: (no, nv) differ between integrals, amplitudes, "// &
                             "intermediates and scratch")
         return
      end if
      if (scr%stamp /= ints%stamp) then
         call scr%status%set(ERROR_VALIDATION, STALE_MSG)
         return
      end if
      if (size(scr%n1, kind=int64) /= int(ints%no, int64)*ints%nv*ints%naux) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: the scratch was sized for another auxiliary count")
         return
      end if
      ok = .true.
   end function ready

   ! ---------------------------------------------------------------------
   ! the intermediates
   ! ---------------------------------------------------------------------

   subroutine rccsd_inter_build_tau(ints, amps, inter, scr)
      !! tau(i,j,a,b) = t2(i,j,a,b) + t1(i,a) t1(j,b)       einsum 'ia,jb->ijab' + t2
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr

      call scr%status%clear()
      if (.not. ready(ints, amps, inter, scr)) return
      call kernel_tau(ints%no, ints%nv, amps%t1, amps%t2, inter%tau)
   end subroutine rccsd_inter_build_tau

   subroutine rccsd_inter_foo(ints, inter, scr, la)
      !! Foo(k,i) = sum_lcd [2 (kc|ld) - (kd|lc)] tau(i,l,c,d)
      !!
      !! pyscf cc_Foo: 2 einsum('kcld,ilcd->ki', ovov, t2) - einsum('kdlc,ilcd->ki', ovov, t2)
      !! + the same two with t1 t1, i.e. with tau, and without its `+ foo` (module header).
      !! One GEMM, K = nv o nv: L viewed as (k; c l d) against tau packed to (i; c l d).
      !! Needs `rccsd_inter_build_tau`.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv

      call scr%status%clear()
      if (.not. ready_inter(ints, inter, scr)) return
      no = ints%no
      nv = ints%nv
      call kernel_pack_taup(no, nv, inter%tau, scr%pack1)
      call la%gemm("N", "T", no, no, nv*no*nv, 1.0_dp, scr%lov, no, scr%pack1, no, 0.0_dp, inter%foo, no)
   end subroutine rccsd_inter_foo

   subroutine rccsd_inter_fvv(ints, inter, scr, la)
      !! Fvv(a,c) = - sum_kld [2 (kc|ld) - (kd|lc)] tau(k,l,a,d)
      !!
      !! pyscf cc_Fvv: -2 einsum('kcld,klad->ac', ovov, t2) + einsum('kdlc,klad->ac', ovov, t2)
      !! and the t1 t1 pair, without its `+ fvv`. One GEMM, K = o o v: tau packed to
      !! (k l d; a) against L packed to (k l d; c). Needs `rccsd_inter_build_tau`.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv

      call scr%status%clear()
      if (.not. ready_inter(ints, inter, scr)) return
      no = ints%no
      nv = ints%nv
      call kernel_pack_taur(no, nv, inter%tau, scr%pack1)
      call la%gemm("T", "N", nv, nv, no*no*nv, -1.0_dp, scr%pack1, no*no*nv, scr%lq, no*no*nv, 0.0_dp, &
                   inter%fvv, nv)
   end subroutine rccsd_inter_fvv

   subroutine rccsd_inter_fov(ints, amps, inter, scr, la)
      !! Fov(k,c) = sum_ld [2 (kc|ld) - (kd|lc)] t1(l,d)
      !!
      !! pyscf cc_Fov: 2 einsum('kcld,ld->kc', ovov, t1) - einsum('kdlc,ld->kc', ovov, t1),
      !! without its `+ fov` (zero for canonical orbitals). L as an (o v) x (o v) matrix times t1.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_linalg_t), intent(in) :: la

      integer :: ov

      call scr%status%clear()
      if (.not. ready(ints, amps, inter, scr)) return
      ov = ints%no*ints%nv
      call la%gemm("N", "N", ov, 1, ov, 1.0_dp, scr%lov, ov, amps%t1, ov, 0.0_dp, inter%fov, ov)
   end subroutine rccsd_inter_fov

   subroutine rccsd_inter_loo(ints, amps, inter, scr, la)
      !! Loo(k,i) = Foo(k,i) + sum_lc [2 (ki|lc) - (li|kc)] t1(l,c)
      !!
      !! pyscf Loo: Foo + 2 einsum('lcki,lc->ki', ovoo, t1) - einsum('kcli,lc->ki', ovoo, t1),
      !! with ovoo(l,c,k,i) = (lc|ki) = ooov(k,i,l,c) and ovoo(k,c,l,i) = ooov(l,i,k,c); the
      !! `fov t1` term is zero. Combination packed to (k i; l c), one matrix-vector product
      !! added onto a copy of Foo. Needs `rccsd_inter_foo`.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv

      call scr%status%clear()
      if (.not. ready(ints, amps, inter, scr)) return
      no = ints%no
      nv = ints%nv
      call kernel_copy(no*no, inter%foo, inter%loo)
      call la%gemm("N", "N", no*no, 1, no*nv, 1.0_dp, scr%lam, no*no, amps%t1, no*nv, 1.0_dp, inter%loo, no*no)
   end subroutine rccsd_inter_loo

   subroutine rccsd_inter_lvv(ints, amps, inter, scr, la)
      !! Lvv(a,c) = Fvv(a,c) + sum_kd [2 (kd|ac) - (kc|ad)] t1(k,d), and N1 in `scr%n1`
      !!
      !! pyscf Lvv: Fvv + 2 einsum('kdac,kd->ac', ovvv, t1) - einsum('kcad,kd->ac', ovvv, t1),
      !! the `fov t1` term zero. From the factors, (kd|ac) = sum_P b_ov(kd,P) b_vv(ac,P) (module
      !! header): N1(k,a;P) = sum_d t1(k,d) b_vv(d,a;P) is ONE GEMM, the direct term 2 b_vv u
      !! (u = b_ov^T t1) two matrix-vector products, the exchange term one GEMM over (k,P) of N1
      !! and b_ov, both permuted to (k,P,.). N1 stays in `scr%n1` for the T2 equation of the same
      !! amplitudes (trc_rccsd_update). Needs `rccsd_inter_fvv`.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, naux, nov, nv2

      call scr%status%clear()
      if (.not. ready(ints, amps, inter, scr)) return
      no = ints%no
      nv = ints%nv
      naux = ints%naux
      nov = no*nv
      nv2 = nv*nv
      ! N1(k,(a,P)) = sum_d t1(k,d) b_vv(d,(a,P)): b_vv is the matrix (d; (a,P)) with leading dimension nv
      call la%gemm("N", "N", no, nv*naux, nv, 1.0_dp, amps%t1, no, ints%b_vv, nv, 0.0_dp, scr%n1, no)
      call kernel_copy(nv2, inter%fvv, inter%lvv)
      ! + 2 sum_P b_vv(ac,P) u(P)
      call la%gemm("T", "N", naux, 1, nov, 1.0_dp, ints%b_ov, nov, amps%t1, nov, 0.0_dp, scr%u, naux)
      call la%gemm("N", "N", nv2, 1, naux, 2.0_dp, ints%b_vv, nv2, scr%u, naux, 1.0_dp, inter%lvv, nv2)
      ! - sum_(k,P) N1(k,a;P) b_ov(k,c;P), the (k,P) pair made contiguous in both
      call kernel_pack_kpa(no, nv, naux, scr%n1, scr%n1p)
      call la%gemm("T", "N", nv, nv, no*naux, -1.0_dp, scr%n1p, no*naux, scr%bovp, no*naux, 1.0_dp, inter%lvv, nv)
   end subroutine rccsd_inter_lvv

   pure function rccsd_inter_flops(no, nv, naux) result(flops)
      !! Flops of one `rccsd_inter_one_particle`, 2 m n k of each of its GEMMs (the packs and tau are kernels, not counted):
      !!   Foo  (o, o, v o v)      2 o^3 v^2         Fvv  (v, v, o o v)      2 o^2 v^3
      !!   Fov  (o v, 1, o v)      2 o^2 v^2         Loo  (o o, 1, o v)      2 o^3 v
      !!   N1   (o, v naux, v)     2 o v^2 naux      Lvv  u (naux, 1, o v)   2 o v naux
      !!   Lvv  direct (v^2, 1, naux) 2 v^2 naux     Lvv  exchange (v, v, o naux)  2 o v^2 naux
      integer, intent(in) :: no, nv, naux
      integer(int64) :: flops

      flops = rccsd_steps_gemm_flops(no, no, nv*no*nv) + rccsd_steps_gemm_flops(nv, nv, no*no*nv) + &
              rccsd_steps_gemm_flops(no*nv, 1, no*nv) + rccsd_steps_gemm_flops(no*no, 1, no*nv) + &
              rccsd_steps_gemm_flops(no, nv*naux, nv) + rccsd_steps_gemm_flops(naux, 1, no*nv) + &
              rccsd_steps_gemm_flops(nv*nv, 1, naux) + rccsd_steps_gemm_flops(nv, nv, no*naux)
   end function rccsd_inter_flops

   subroutine rccsd_inter_one_particle(ints, amps, inter, scr, la)
      !! tau, Foo, Fvv, Fov, Loo, Lvv for the current amplitudes (`amps` on the device: `push` after a host edit)
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_linalg_t), intent(in) :: la

      call rccsd_inter_build_tau(ints, amps, inter, scr)
      if (scr%status%has_error()) return
      call rccsd_inter_foo(ints, inter, scr, la)
      if (scr%status%has_error()) return
      call rccsd_inter_fvv(ints, inter, scr, la)
      if (scr%status%has_error()) return
      call rccsd_inter_fov(ints, amps, inter, scr, la)
      if (scr%status%has_error()) return
      call rccsd_inter_loo(ints, amps, inter, scr, la)
      if (scr%status%has_error()) return
      call rccsd_inter_lvv(ints, amps, inter, scr, la)
   end subroutine rccsd_inter_one_particle

   ! ---------------------------------------------------------------------
   ! the energy and the MP2 guess
   ! ---------------------------------------------------------------------

   subroutine rccsd_inter_energy(amps, scr, e_singles, e_doubles)
      !! E = sum_ijab [2 (ia|jb) - (ib|ja)] (t2(i,j,a,b) + t1(i,a) t1(j,b))
      !!
      !! rccsd.energy: 2 einsum('ijab,iajb', tau, ovov) - einsum('ijab,ibja', tau, ovov), its
      !! `2 f_ov t1` term zero. Split as mqc does: `e_doubles` is what t2 alone gives and
      !! `e_singles` the t1 t1 part. One partial per (j,b), added in ascending order on the host.
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      real(dp), intent(out) :: e_singles, e_doubles

      integer :: p

      e_singles = 0.0_dp
      e_doubles = 0.0_dp
      call scr%status%clear()
      if (.not. allocated(scr%lov) .or. .not. allocated(amps%t1)) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: rccsd_inter_energy before the scratch and amplitudes exist")
         return
      end if
      if (amps%no /= scr%no .or. amps%nv /= scr%nv) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: rccsd_inter_energy: amplitudes and scratch differ in (no, nv)")
         return
      end if
      call kernel_energy(scr%no, scr%nv, scr%lov, amps%t1, amps%t2, scr%part_s, scr%part_d)
      !$acc update self(scr%part_s, scr%part_d)
      do p = 1, scr%no*scr%nv
         e_singles = e_singles + scr%part_s(p)
         e_doubles = e_doubles + scr%part_d(p)
      end do
   end subroutine rccsd_inter_energy

   subroutine rccsd_inter_mp2_guess(ints, amps, scr, e_mp2)
      !! t1 = 0, t2(i,j,a,b) = (ia|jb) / (e_i + e_j - e_a - e_b), and the MP2 energy they give
      !!
      !! The energy is the doubles part of `rccsd_inter_energy` at these amplitudes. Leaves the host
      !! and device copies of `amps` equal.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(inout) :: amps
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      real(dp), intent(out) :: e_mp2

      real(dp) :: e_singles

      e_mp2 = 0.0_dp
      call scr%status%clear()
      if (.not. allocated(scr%lov) .or. .not. allocated(amps%t1) .or. .not. allocated(ints%ovov)) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: rccsd_inter_mp2_guess before the scratch, amplitudes "// &
                             "and integrals exist")
         return
      end if
      if (amps%no /= ints%no .or. amps%nv /= ints%nv .or. scr%no /= ints%no .or. scr%nv /= ints%nv) then
         call scr%status%set(ERROR_VALIDATION, "trc_rccsd_inter: rccsd_inter_mp2_guess: (no, nv) differ")
         return
      end if
      if (scr%stamp /= ints%stamp) then
         call scr%status%set(ERROR_VALIDATION, STALE_MSG)
         return
      end if
      call kernel_mp2(ints%no, ints%nv, ints%eps_o, ints%eps_v, ints%ovov, amps%t1, amps%t2)
      !$acc update self(amps%t1, amps%t2)
      call rccsd_inter_energy(amps, scr, e_singles, e_mp2)
   end subroutine rccsd_inter_mp2_guess

   ! ---------------------------------------------------------------------
   ! the kernels. Each takes its arrays as explicit-shape dummies and runs one
   ! `do concurrent` over ONE flattened index, whose body calls an item routine
   ! of this module. No `x%name` appears in any of them.
   ! ---------------------------------------------------------------------

   subroutine kernel_lov(no, nv, ovov, lov)
      !! lov(k,c,l,d) = 2 ovov(k,c,l,d) - ovov(k,d,l,c)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp), intent(out) :: lov(no*nv*no*nv)
      integer :: t

      do concurrent(t=1:no*nv*no*nv)
         lov(t) = lov_item(t, no, nv, ovov)
      end do
   end subroutine kernel_lov

   subroutine kernel_tau(no, nv, t1, t2, tau)
      !! tau(i,j,a,b) = t2(i,j,a,b) + t1(i,a) t1(j,b)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv)
      real(dp), intent(out) :: tau(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         tau(t) = tau_item(t, no, nv, t1, t2)
      end do
   end subroutine kernel_tau

   subroutine kernel_pack_taup(no, nv, tau, taup)
      !! taup(i,c,l,d) = tau(i,l,c,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(out) :: taup(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         taup(t) = taup_item(t, no, nv, tau)
      end do
   end subroutine kernel_pack_taup

   subroutine kernel_pack_taur(no, nv, tau, taur)
      !! taur(k,l,d,a) = tau(k,l,a,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(out) :: taur(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         taur(t) = taur_item(t, no, nv, tau)
      end do
   end subroutine kernel_pack_taur

   subroutine kernel_pack_lq(no, nv, lov, lq)
      !! lq(k,l,d,c) = lov(k,c,l,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: lov(no, nv, no, nv)
      real(dp), intent(out) :: lq(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         lq(t) = lq_item(t, no, nv, lov)
      end do
   end subroutine kernel_pack_lq

   subroutine kernel_pack_lam(no, nv, ooov, lam)
      !! lam(k,i,l,c) = 2 ooov(k,i,l,c) - ooov(l,i,k,c)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ooov(no, no, no, nv)
      real(dp), intent(out) :: lam(no*no*no*nv)
      integer :: t

      do concurrent(t=1:no*no*no*nv)
         lam(t) = lam_item(t, no, nv, ooov)
      end do
   end subroutine kernel_pack_lam

   subroutine kernel_copy(n, src, dst)
      !! dst = src
      integer, intent(in) :: n
      real(dp), intent(in) :: src(n)
      real(dp), intent(out) :: dst(n)
      integer :: t

      do concurrent(t=1:n)
         dst(t) = src(t)
      end do
   end subroutine kernel_copy

   subroutine kernel_pack_kpa(no, nv, naux, x, y)
      !! y(k,P,a) = x(k,a,P), for an (o v, naux) factor such as b_ov or N1
      integer, intent(in) :: no, nv, naux
      real(dp), intent(in) :: x(no, nv, naux)
      real(dp), intent(out) :: y(no*naux*nv)
      integer :: t

      do concurrent(t=1:no*naux*nv)
         y(t) = kpa_item(t, no, nv, naux, x)
      end do
   end subroutine kernel_pack_kpa

   subroutine kernel_energy(no, nv, lov, t1, t2, part_s, part_d)
      !! one (singles, doubles) partial per (j,b), j fastest
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: lov(no, nv, no, nv), t1(no, nv), t2(no, no, nv, nv)
      real(dp), intent(out) :: part_s(no*nv), part_d(no*nv)
      integer :: t

      do concurrent(t=1:no*nv)
         call energy_item(t, no, nv, lov, t1, t2, part_s(t), part_d(t))
      end do
   end subroutine kernel_energy

   subroutine kernel_mp2(no, nv, eps_o, eps_v, ovov, t1, t2)
      !! t1 = 0, t2(i,j,a,b) = ovov(i,a,j,b) / (e_i + e_j - e_a - e_b)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: eps_o(no), eps_v(nv), ovov(no, nv, no, nv)
      real(dp), intent(out) :: t1(no*nv), t2(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*nv)
         t1(t) = 0.0_dp
      end do
      do concurrent(t=1:no*no*nv*nv)
         t2(t) = mp2_item(t, no, nv, eps_o, eps_v, ovov)
      end do
   end subroutine kernel_mp2

   ! ---------------------------------------------------------------------
   ! the item routines: pure, `routine seq`, explicit-shape dummies, no
   ! arrays of their own. A flat 1-based index t is the position in the
   ! DESTINATION, first index fastest.
   ! The loop index is a VALUE dummy: see the item routines of trc_rccsd_update
   ! (an ifx miscompilation of an index passed by reference).
   ! ---------------------------------------------------------------------

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

   pure function lov_item(t, no, nv, ovov) result(v)
      !! L(k,c,l,d) at position t of (no, nv, no, nv)
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp) :: v
      integer :: k, c, l, d

      call split4(t, no, nv, no, k, c, l, d)
      v = 2.0_dp*ovov(k, c, l, d) - ovov(k, d, l, c)
   end function lov_item

   pure function tau_item(t, no, nv, t1, t2) result(v)
      !! tau(i,j,a,b) at position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv)
      real(dp) :: v
      integer :: i, j, a, b

      call split4(t, no, no, nv, i, j, a, b)
      v = t2(i, j, a, b) + t1(i, a)*t1(j, b)
   end function tau_item

   pure function taup_item(t, no, nv, tau) result(v)
      !! taup(i,c,l,d) = tau(i,l,c,d), position t of (no, nv, no, nv)
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp) :: v
      integer :: i, c, l, d

      call split4(t, no, nv, no, i, c, l, d)
      v = tau(i, l, c, d)
   end function taup_item

   pure function taur_item(t, no, nv, tau) result(v)
      !! taur(k,l,d,a) = tau(k,l,a,d), position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp) :: v
      integer :: k, l, d, a

      call split4(t, no, no, nv, k, l, d, a)
      v = tau(k, l, a, d)
   end function taur_item

   pure function lq_item(t, no, nv, lov) result(v)
      !! lq(k,l,d,c) = L(k,c,l,d), position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: lov(no, nv, no, nv)
      real(dp) :: v
      integer :: k, l, d, c

      call split4(t, no, no, nv, k, l, d, c)
      v = lov(k, c, l, d)
   end function lq_item

   pure function lam_item(t, no, nv, ooov) result(v)
      !! lam(k,i,l,c) = 2 ooov(k,i,l,c) - ooov(l,i,k,c), position t of (no, no, no, nv)
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ooov(no, no, no, nv)
      real(dp) :: v
      integer :: k, i, l, c

      call split4(t, no, no, no, k, i, l, c)
      v = 2.0_dp*ooov(k, i, l, c) - ooov(l, i, k, c)
   end function lam_item

   pure function kpa_item(t, no, nv, naux, x) result(v)
      !! x(k,a,P) at position t of (no, naux, nv)
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv, naux
      real(dp), intent(in) :: x(no, nv, naux)
      real(dp) :: v
      integer :: k, p, a, unused

      call split4(t, no, naux, nv, k, p, a, unused)
      v = x(k, a, p)
   end function kpa_item

   pure subroutine energy_item(t, no, nv, lov, t1, t2, es, ed)
      !! The (singles, doubles) partial of the pair t = j + no (b-1); a outer, i inner
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: lov(no, nv, no, nv), t1(no, nv), t2(no, no, nv, nv)
      real(dp), intent(out) :: es, ed
      integer :: i, j, a, b
      real(dp) :: l

      j = mod(t - 1, no) + 1
      b = (t - 1)/no + 1
      es = 0.0_dp
      ed = 0.0_dp
      do a = 1, nv
         do i = 1, no
            l = lov(i, a, j, b)
            ed = ed + l*t2(i, j, a, b)
            es = es + l*t1(i, a)*t1(j, b)
         end do
      end do
   end subroutine energy_item

   pure function mp2_item(t, no, nv, eps_o, eps_v, ovov) result(v)
      !! t2(i,j,a,b) at position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: eps_o(no), eps_v(nv), ovov(no, nv, no, nv)
      real(dp) :: v
      integer :: i, j, a, b

      call split4(t, no, no, nv, i, j, a, b)
      v = ovov(i, a, j, b)/(eps_o(i) + eps_o(j) - eps_v(a) - eps_v(b))
   end function mp2_item

end module trc_rccsd_inter
