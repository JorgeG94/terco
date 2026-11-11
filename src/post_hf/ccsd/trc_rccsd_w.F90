!
! The two-particle intermediates Woooo, Wvoov and Wvovo of a closed-shell
! RI-CCSD(T) calculation.
!
! Hirata, Podeszwa, Tobita and Bartlett, JCP 120, 2581 (2004), as
! pyscf/cc/rintermediates.py's cc_Woooo, cc_Wvoov and cc_Wvovo write them,
! ported from the same author's mqc_czt_rcc.f90 (cc_woooo, cc_wvoov, cc_wvovo,
! accumulate_wvoov, accumulate_wvovo). Wvvvv, the ladder and the amplitude
! update are not here. Chemists' notation, spatial orbitals, REAL orbitals
! (so (ia|jb) = (jb|ia) and (ij|ka) = (ka|ij)); the integral blocks and the
! amplitudes are those of trc_rccsd_ints and trc_rccsd_inter.
!
! LAYOUT, as in mqc (first index fastest), which is PySCF's index order for
! each of them, so no permutation is needed to compare:
!
!     Woooo(k,l,i,j)   Wklij
!     Wvoov(a,k,i,c)   Wakic
!     Wvovo(a,k,c,i)   Wakci
!
! THE TERMS (ovoo(l,c,k,i) = (lc|ki) = ooov(k,i,l,c); ovvo(k,c,a,i) = (kc|ai)
! = ovov(k,c,i,a); every integral below is written in the index order of the
! array that holds it here).
!
!   Woooo(k,l,i,j) = (ki|lj)                              oooo(k,i,l,j)
!     + sum_c ooov(k,i,l,c) t1(j,c)        einsum 'lcki,jc->klij'   (ovoo)
!     + sum_c ooov(l,j,k,c) t1(i,c)        einsum 'kclj,ic->klij'
!     + sum_cd ovov(k,c,l,d) tau(i,j,c,d)  einsum 'kcld,ijcd->klij' + 'kcld,ic,jd->klij'
!
!   Wvoov(a,k,i,c) =
!     + ovov(k,c,i,a)                                      eris.ovvo.transpose(2,0,3,1)
!     + sum_d ovvv(k,c,a,d) t1(i,d)                        'kcad,id->akic'
!     - sum_l ooov(l,i,k,c) t1(l,a)                        'kcli,la->akic'
!     - 1/2 sum_ld ovov(l,d,k,c) t2(i,l,d,a)               'ldkc,ilda->akic'
!     - 1/2 sum_ld ovov(l,c,k,d) t2(i,l,a,d)               'lckd,ilad->akic'
!     -     sum_ld ovov(l,d,k,c) t1(i,d) t1(l,a)           'ldkc,id,la->akic'
!     +     sum_ld ovov(l,d,k,c) t2(i,l,a,d)               'ldkc,ilad->akic'
!
!   Wvovo(a,k,c,i) =
!     + oovv(k,i,a,c)                                      eris.oovv.transpose(2,0,3,1)
!     + sum_d ovvv(k,d,a,c) t1(i,d)                        'kdac,id->akci'
!     - sum_l ooov(k,i,l,c) t1(l,a)                        'lcki,la->akci'
!     - 1/2 sum_ld ovov(l,c,k,d) t2(i,l,d,a)               'lckd,ilda->akci'
!     -     sum_ld ovov(l,c,k,d) t1(i,d) t1(l,a)           'lckd,id,la->akci'
!
! HOW EACH BECOMES GEMMS. The terms with the contracted index inside two
! compound pairs are the O(N^6) ones; the rest are O(N^5) or less and are one
! `do concurrent` thread per element of W, summing the contracted index in
! ascending order (like Lvv in trc_rccsd_inter), so they are deterministic.
!
!   Woooo   the tau term is one GEMM, M = o^2, N = o^2, K = v^2 (o^4 v^2):
!             W(kl;ij) += akl(kl;cd) . tau(ij;cd)^T         akl(k,l,c,d) = ovov(k,c,l,d)
!           tau is already (ij;cd) in memory. The copy and the two t1 terms
!           (o^4 v) are the base kernel that runs first.
!   Wvoov   the four amplitude terms group by the integral ordering they carry
!           (mqc's grouping) into two products over the compound (l,d), both
!           accumulated in R((k,c),(i,a)) (M = N = K = o v, o^3 v^3):
!             R  = a1(kc;ld) . b1(ld;ia) + a2(kc;ld) . b2(ld;ia)
!             a1(k,c,l,d) = ovov(l,d,k,c)   (the stored block as a (ld;kc) matrix, so a1 is
!                                                            its transpose: the GEMM reads it with "T")
!             a2(k,c,l,d) = ovov(l,c,k,d)                   (a pack, built once)
!             b1(l,d,i,a) = -1/2 t2(i,l,d,a) - t1(i,d) t1(l,a) + t2(i,l,a,d)
!             b2(l,d,i,a) = -1/2 t2(i,l,a,d)
!           then Wvoov(a,k,i,c) += R(k,c,i,a). The base kernel (ovov, ovvv t1, ooov t1)
!           runs first.
!   Wvovo   one product, a2(kc;ld) . b3(ld;ia) with b3 = -1/2 t2(i,l,d,a) - t1(i,d) t1(l,a)
!           (packed into the b1 buffer), then Wvovo(a,k,c,i) += R(k,c,i,a).
!
! DEVIATIONS FROM mqc, none of which changes a number beyond rounding:
!   - a1 is `ints%ovov` itself read with "T", not mqc's pack ovov(l,d,k,c) (no copy of
!     o^2 v^2 doubles; the stored block is symmetric only to rounding, which "T" respects);
!   - the packs that depend on the integrals only (akl, a2) are built once by `init`
!     of the scratch, not on every call; the a2 of Wvoov and of Wvovo is one array;
!   - the two Wvoov products are chained in one R (second GEMM with beta = 1) and
!     accumulated into W once, so W = base + (R1 + R2) where mqc has (base + R1) + R2;
!   - Woooo takes tau from `inter` (built by `rccsd_inter_build_tau`), as mqc takes
!     its `tau` argument.
!
! RESIDENCY. Every array of the types here exists on the host AND on the device
! (`enter data` in `init`, `exit data` in `release`, on the same `x%name`
! expression). The device copy is the one the kernels write; `pull` brings W
! back to the host, `push` sends the host copy up. Each component is mapped on
! its own and the parent type is never mapped, so a kernel gets an array only as
! an explicit-shape dummy of a `kernel_*` routine and never names `x%name` in a
! `do concurrent` body (an illegal address under mem:separate that the host
! build hides). The scratch's host copies of b1, b2, r are never read.
!
! SCRATCH, allocated once by `trc_rccsd_wscratch_t%init` (o = active occupied,
! v = virtual, 8 bytes each, host and device; no = 25, nv = 130 in the last column):
!
!     akl      o^2 v^2   ovov(k,c,l,d) as (kl;cd), for life of the integrals     84.5 MB
!     a2       o^2 v^2   ovov(l,c,k,d) as (kc;ld), for life of the integrals     84.5 MB
!     b1       o^2 v^2   b1 (Wvoov), then b3 (Wvovo): per call
!     b2       o^2 v^2   b2 (Wvoov): per call
!     r        o^2 v^2   the GEMM result R(kc;ia): per call
!
! and `trc_rccsd_w_t` holds Woooo (o^4), Wvoov, Wvovo (o^2 v^2 each). In all
! 7 o^2 v^2 + o^4 doubles: 594.6 MB a side at no = 25, nv = 130. A block above
! 2^31-1 elements is refused by either `init` (the scratch's in `ws%status`, W's in
! `w%status`), not wrapped.
!
! STALE PACKS. `init` of the scratch packs akl and a2 from the integrals once and
! records `ints%stamp`. After `ints%build` is called again (new geometry or
! orbitals, even with the same no and nv) every routine here refuses with a
! message and writes nothing: re-`init` the scratch after `ints%build`.
!
! Failures are reported in `ws%status` (a routine whose arguments differ in
! (no, nv), or are not built, sets it and returns without writing anything);
! none stops. Single rank only.
module trc_rccsd_w
   use pic_types, only: dp, int32, int64
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_steps, only: rccsd_steps_gemm_flops
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   character(len=*), parameter :: STALE_MSG = "trc_rccsd_w: the integrals were rebuilt after the scratch was "// &
                                              "initialised; re-init the scratch after ints%build"

   public :: trc_rccsd_w_t, trc_rccsd_wscratch_t
   public :: rccsd_w_woooo, rccsd_w_wvoov, rccsd_w_wvovo, rccsd_w_build, rccsd_w_flops

   type :: trc_rccsd_w_t
      !! Woooo, Wvoov and Wvovo on host and device, valid after the matching routine
      type(error_t) :: status  !! set by `init`: a refusal leaves its message here and the type empty
      integer :: no = 0
      integer :: nv = 0
      real(dp), allocatable :: woooo(:, :, :, :)  !! Woooo(k,l,i,j)
      real(dp), allocatable :: wvoov(:, :, :, :)  !! Wvoov(a,k,i,c)
      real(dp), allocatable :: wvovo(:, :, :, :)  !! Wvovo(a,k,c,i)
   contains
      procedure :: init => w_init
      procedure :: release => w_release
      procedure :: push => w_push
      procedure :: pull => w_pull
   end type trc_rccsd_w_t

   type :: trc_rccsd_wscratch_t
      !! The work arrays of the routines below, sizes in the module header
      type(error_t) :: status  !! set by `init` and by every routine that refuses its arguments
      integer :: no = 0
      integer :: nv = 0
      integer :: stamp = 0  !! `ints%stamp` of the integrals `init` packed akl and a2 from
      real(dp), allocatable :: akl(:)  !! o^2 v^2, akl(k,l,c,d) = ovov(k,c,l,d), built by `init`
      real(dp), allocatable :: a2(:)   !! o^2 v^2, a2(k,c,l,d) = ovov(l,c,k,d), built by `init`
      real(dp), allocatable :: b1(:)   !! o^2 v^2, b1(l,d,i,a) of Wvoov, then b3 of Wvovo
      real(dp), allocatable :: b2(:)   !! o^2 v^2, b2(l,d,i,a)
      real(dp), allocatable :: r(:)    !! o^2 v^2, R(k,c,i,a)
   contains
      procedure :: init => wscratch_init
      procedure :: release => wscratch_release
   end type trc_rccsd_wscratch_t

contains

   ! ---------------------------------------------------------------------
   ! the types
   ! ---------------------------------------------------------------------

   subroutine w_init(this, no, nv)
      !! Zeroed W on host and device. A block above 2^31-1 elements (o^4 or o^2 v^2) is refused in `status`.
      class(trc_rccsd_w_t), intent(inout) :: this
      integer, intent(in) :: no, nv

      call this%release()
      call this%status%clear()
      if (no < 0 .or. nv < 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_w: negative orbital count")
         return
      end if
      if (elements(no, no, nv, nv) > int(huge(0_int32), int64) .or. elements(no, no, no, no) > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_w: a block has more than 2^31-1 elements")
         return
      end if
      this%no = no
      this%nv = nv
      allocate (this%woooo(no, no, no, no), this%wvoov(nv, no, no, nv), this%wvovo(nv, no, nv, no))
      this%woooo = 0.0_dp
      this%wvoov = 0.0_dp
      this%wvovo = 0.0_dp
      !$acc enter data copyin(this%woooo, this%wvoov, this%wvovo)
   end subroutine w_init

   subroutine w_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_w_t), intent(inout) :: this

      if (allocated(this%woooo)) then
         !$acc exit data delete(this%woooo)
         deallocate (this%woooo)
      end if
      if (allocated(this%wvoov)) then
         !$acc exit data delete(this%wvoov)
         deallocate (this%wvoov)
      end if
      if (allocated(this%wvovo)) then
         !$acc exit data delete(this%wvovo)
         deallocate (this%wvovo)
      end if
      this%no = 0
      this%nv = 0
   end subroutine w_release

   subroutine w_push(this)
      !! Host copy to device
      class(trc_rccsd_w_t), intent(inout) :: this

      if (.not. allocated(this%woooo)) return
      !$acc update device(this%woooo, this%wvoov, this%wvovo)
   end subroutine w_push

   subroutine w_pull(this)
      !! Device copy to host
      class(trc_rccsd_w_t), intent(inout) :: this

      if (.not. allocated(this%woooo)) return
      !$acc update self(this%woooo, this%wvoov, this%wvovo)
   end subroutine w_pull

   subroutine wscratch_init(this, ints)
      !! Allocates the work arrays for these integrals and builds akl and a2. Refusals go to `status`.
      class(trc_rccsd_wscratch_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints

      integer :: no, nv

      call this%release()
      call this%status%clear()
      if (.not. allocated(ints%ovov)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_w: the integral blocks are not built")
         return
      end if
      no = ints%no
      nv = ints%nv
      if (elements(no, no, nv, nv) > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_w: a block has more than 2^31-1 elements")
         return
      end if
      this%no = no
      this%nv = nv
      this%stamp = ints%stamp
      allocate (this%akl(no*no*nv*nv), this%a2(no*no*nv*nv), this%b1(no*no*nv*nv), this%b2(no*no*nv*nv))
      allocate (this%r(no*no*nv*nv))
      !$acc enter data create(this%akl, this%a2, this%b1, this%b2, this%r)
      call kernel_pack_akl(no, nv, ints%ovov, this%akl)
      call kernel_pack_a2(no, nv, ints%ovov, this%a2)
   end subroutine wscratch_init

   subroutine wscratch_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_wscratch_t), intent(inout) :: this

      if (allocated(this%akl)) then
         !$acc exit data delete(this%akl)
         deallocate (this%akl)
      end if
      if (allocated(this%a2)) then
         !$acc exit data delete(this%a2)
         deallocate (this%a2)
      end if
      if (allocated(this%b1)) then
         !$acc exit data delete(this%b1)
         deallocate (this%b1)
      end if
      if (allocated(this%b2)) then
         !$acc exit data delete(this%b2)
         deallocate (this%b2)
      end if
      if (allocated(this%r)) then
         !$acc exit data delete(this%r)
         deallocate (this%r)
      end if
      this%no = 0
      this%nv = 0
      this%stamp = 0
   end subroutine wscratch_release

   ! ---------------------------------------------------------------------
   ! size checks. A refusal sets ws%status and the caller returns before it
   ! writes anything.
   ! ---------------------------------------------------------------------

   pure function elements(n1, n2, n3, n4) result(n)
      !! n1 n2 n3 n4 as a 64-bit count
      integer, intent(in) :: n1, n2, n3, n4
      integer(int64) :: n

      n = int(n1, int64)*n2*n3*n4
   end function elements

   function ready(ints, amps, w, ws) result(ok)
      !! True when the integrals, amplitudes, W and the scratch are built, agree on (no, nv) and hold
      !! arrays of those sizes; else `ws%status` says which not
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_w_t), intent(in) :: w
      type(trc_rccsd_wscratch_t), intent(inout) :: ws
      logical :: ok

      integer :: no, nv

      ok = .false.
      if (.not. allocated(ints%oooo) .or. .not. allocated(ints%ooov) .or. .not. allocated(ints%oovv) .or. &
          .not. allocated(ints%ovov) .or. .not. allocated(ints%ovvv)) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: the integral blocks are not built")
         return
      end if
      if (.not. allocated(ws%akl) .or. .not. allocated(ws%a2) .or. .not. allocated(ws%b1) .or. &
          .not. allocated(ws%b2) .or. .not. allocated(ws%r)) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: the scratch is not initialised")
         return
      end if
      if (.not. allocated(amps%t1) .or. .not. allocated(amps%t2)) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: the amplitudes are not initialised")
         return
      end if
      if (.not. allocated(w%woooo) .or. .not. allocated(w%wvoov) .or. .not. allocated(w%wvovo)) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: the W arrays are not initialised")
         return
      end if
      no = ints%no
      nv = ints%nv
      if (amps%no /= no .or. amps%nv /= nv .or. w%no /= no .or. w%nv /= nv .or. ws%no /= no .or. ws%nv /= nv) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: (no, nv) differ between integrals, amplitudes, "// &
                            "W and scratch")
         return
      end if
      if (size(ints%oooo, kind=int64) /= elements(no, no, no, no) .or. &
          size(ints%ooov, kind=int64) /= elements(no, no, no, nv) .or. &
          size(ints%oovv, kind=int64) /= elements(no, no, nv, nv) .or. &
          size(ints%ovov, kind=int64) /= elements(no, no, nv, nv) .or. &
          size(ints%ovvv, kind=int64) /= elements(no, nv, nv, nv)) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: an integral block does not have the size (no, nv) says")
         return
      end if
      if (size(amps%t1, kind=int64) /= elements(no, nv, 1, 1) .or. size(amps%t2, kind=int64) /= elements(no, no, nv, nv) .or. &
          size(w%woooo, kind=int64) /= elements(no, no, no, no) .or. size(w%wvoov, kind=int64) /= elements(no, no, nv, nv) .or. &
          size(w%wvovo, kind=int64) /= elements(no, no, nv, nv)) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: an amplitude or W array does not have the size (no, nv) says")
         return
      end if
      if (size(ws%akl, kind=int64) /= elements(no, no, nv, nv) .or. size(ws%a2, kind=int64) /= elements(no, no, nv, nv) .or. &
          size(ws%b1, kind=int64) /= elements(no, no, nv, nv) .or. size(ws%b2, kind=int64) /= elements(no, no, nv, nv) .or. &
          size(ws%r, kind=int64) /= elements(no, no, nv, nv)) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: a scratch array does not have the size (no, nv) says")
         return
      end if
      if (ws%stamp /= ints%stamp) then
         call ws%status%set(ERROR_VALIDATION, STALE_MSG)
         return
      end if
      ok = .true.
   end function ready

   function ready_tau(ints, inter, ws) result(ok)
      !! True when `inter` holds a tau of the integrals' (no, nv); else `ws%status` says so
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_wscratch_t), intent(inout) :: ws
      logical :: ok

      ok = .false.
      if (.not. allocated(inter%tau)) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: the intermediates (tau) are not initialised")
         return
      end if
      if (inter%no /= ints%no .or. inter%nv /= ints%nv .or. &
          size(inter%tau, kind=int64) /= elements(ints%no, ints%no, ints%nv, ints%nv)) then
         call ws%status%set(ERROR_VALIDATION, "trc_rccsd_w: tau does not have the size (no, nv) says")
         return
      end if
      ok = .true.
   end function ready_tau

   ! ---------------------------------------------------------------------
   ! the intermediates
   ! ---------------------------------------------------------------------

   subroutine rccsd_w_woooo(ints, amps, inter, w, ws, la)
      !! Woooo(k,l,i,j) = (ki|lj) + sum_c ooov(k,i,l,c) t1(j,c) + sum_c ooov(l,j,k,c) t1(i,c)
      !!                  + sum_cd ovov(k,c,l,d) tau(i,j,c,d)
      !!
      !! pyscf cc_Woooo: einsum('lcki,jc->klij', ovoo, t1) + einsum('kclj,ic->klij', ovoo, t1)
      !! + einsum('kcld,ijcd->klij', ovov, t2) + einsum('kcld,ic,jd->klij', ovov, t1, t1)
      !! + oooo.transpose(0,2,1,3); the last two ovov terms are the tau term. Needs
      !! `rccsd_inter_build_tau` at the current amplitudes. One base kernel, then one GEMM
      !! with K = nv nv (module header).
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_w_t), intent(inout) :: w
      type(trc_rccsd_wscratch_t), intent(inout) :: ws
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv

      call ws%status%clear()
      if (.not. ready(ints, amps, w, ws)) return
      if (.not. ready_tau(ints, inter, ws)) return
      no = ints%no
      nv = ints%nv
      call kernel_woooo_base(no, nv, ints%oooo, ints%ooov, amps%t1, w%woooo)
      call la%gemm("N", "T", no*no, no*no, nv*nv, 1.0_dp, ws%akl, no*no, inter%tau, no*no, 1.0_dp, w%woooo, no*no)
   end subroutine rccsd_w_woooo

   subroutine rccsd_w_wvoov(ints, amps, w, ws, la)
      !! Wvoov(a,k,i,c): the seven terms of the module header
      !!
      !! pyscf cc_Wvoov. The base kernel writes the first three, two GEMMs into one R the
      !! four amplitude terms, an accumulate kernel adds R. Needs the amplitudes on the device.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_w_t), intent(inout) :: w
      type(trc_rccsd_wscratch_t), intent(inout) :: ws
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, nov

      call ws%status%clear()
      if (.not. ready(ints, amps, w, ws)) return
      no = ints%no
      nv = ints%nv
      nov = no*nv
      call kernel_voov_base(no, nv, ints%ovov, ints%ovvv, ints%ooov, amps%t1, w%wvoov)
      call kernel_pack_b1(no, nv, amps%t1, amps%t2, ws%b1)
      call kernel_pack_b2(no, nv, amps%t2, ws%b2)
      call la%gemm("T", "N", nov, nov, nov, 1.0_dp, ints%ovov, nov, ws%b1, nov, 0.0_dp, ws%r, nov)
      call la%gemm("N", "N", nov, nov, nov, 1.0_dp, ws%a2, nov, ws%b2, nov, 1.0_dp, ws%r, nov)
      call kernel_acc_voov(no, nv, ws%r, w%wvoov)
   end subroutine rccsd_w_wvoov

   subroutine rccsd_w_wvovo(ints, amps, w, ws, la)
      !! Wvovo(a,k,c,i): the five terms of the module header
      !!
      !! pyscf cc_Wvovo. The base kernel writes the first three, one GEMM the two amplitude
      !! terms, an accumulate kernel adds R. Overwrites the scratch b1 (with b3). Needs the
      !! amplitudes on the device (the caller pushes after any host edit).
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_w_t), intent(inout) :: w
      type(trc_rccsd_wscratch_t), intent(inout) :: ws
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, nov

      call ws%status%clear()
      if (.not. ready(ints, amps, w, ws)) return
      no = ints%no
      nv = ints%nv
      nov = no*nv
      call kernel_vovo_base(no, nv, ints%oovv, ints%ovvv, ints%ooov, amps%t1, w%wvovo)
      call kernel_pack_b3(no, nv, amps%t1, amps%t2, ws%b1)
      call la%gemm("N", "N", nov, nov, nov, 1.0_dp, ws%a2, nov, ws%b1, nov, 0.0_dp, ws%r, nov)
      call kernel_acc_vovo(no, nv, ws%r, w%wvovo)
   end subroutine rccsd_w_wvovo

   pure function rccsd_w_flops(no, nv) result(flops)
      !! Flops of one `rccsd_w_build`, 2 m n k of each of its GEMMs (the base kernels, packs and accumulates are not counted):
      !!   Woooo  (o o, o o, v v)             2 o^4 v^2
      !!   Wvoov  two (o v, o v, o v)         4 (o v)^3
      !!   Wvovo  one (o v, o v, o v)         2 (o v)^3
      integer, intent(in) :: no, nv
      integer(int64) :: flops

      flops = rccsd_steps_gemm_flops(no*no, no*no, nv*nv) + 3_int64*rccsd_steps_gemm_flops(no*nv, no*nv, no*nv)
   end function rccsd_w_flops

   subroutine rccsd_w_build(ints, amps, inter, w, ws, la)
      !! Woooo, Wvoov and Wvovo for the current amplitudes (`amps` on the device: `push` after a host
      !! edit; tau from `rccsd_inter_build_tau`). The first refusal stops the rest.
      !!
      !! All three W are built from the SAME input amplitudes, and tau must have been built from
      !! them. The iteration driver must therefore not overwrite `amps%t1` or `amps%t2` before this
      !! returns, and must rebuild tau after any update of them.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_w_t), intent(inout) :: w
      type(trc_rccsd_wscratch_t), intent(inout) :: ws
      type(trc_linalg_t), intent(in) :: la

      call rccsd_w_woooo(ints, amps, inter, w, ws, la)
      if (ws%status%has_error()) return
      call rccsd_w_wvoov(ints, amps, w, ws, la)
      if (ws%status%has_error()) return
      call rccsd_w_wvovo(ints, amps, w, ws, la)
   end subroutine rccsd_w_build

   ! ---------------------------------------------------------------------
   ! the kernels. Each takes its arrays as explicit-shape dummies and runs one
   ! `do concurrent` over ONE flattened index, whose body calls an item routine
   ! of this module. No `x%name` appears in any of them.
   ! ---------------------------------------------------------------------

   subroutine kernel_pack_akl(no, nv, ovov, akl)
      !! akl(k,l,c,d) = ovov(k,c,l,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp), intent(out) :: akl(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         akl(t) = akl_item(t, no, nv, ovov)
      end do
   end subroutine kernel_pack_akl

   subroutine kernel_pack_a2(no, nv, ovov, a2)
      !! a2(k,c,l,d) = ovov(l,c,k,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp), intent(out) :: a2(no*nv*no*nv)
      integer :: t

      do concurrent(t=1:no*nv*no*nv)
         a2(t) = a2_item(t, no, nv, ovov)
      end do
   end subroutine kernel_pack_a2

   subroutine kernel_pack_b1(no, nv, t1, t2, b1)
      !! b1(l,d,i,a) = -1/2 t2(i,l,d,a) - t1(i,d) t1(l,a) + t2(i,l,a,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv)
      real(dp), intent(out) :: b1(no*nv*no*nv)
      integer :: t

      do concurrent(t=1:no*nv*no*nv)
         b1(t) = b1_item(t, no, nv, t1, t2)
      end do
   end subroutine kernel_pack_b1

   subroutine kernel_pack_b2(no, nv, t2, b2)
      !! b2(l,d,i,a) = -1/2 t2(i,l,a,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(no, no, nv, nv)
      real(dp), intent(out) :: b2(no*nv*no*nv)
      integer :: t

      do concurrent(t=1:no*nv*no*nv)
         b2(t) = b2_item(t, no, nv, t2)
      end do
   end subroutine kernel_pack_b2

   subroutine kernel_pack_b3(no, nv, t1, t2, b3)
      !! b3(l,d,i,a) = -1/2 t2(i,l,d,a) - t1(i,d) t1(l,a)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv)
      real(dp), intent(out) :: b3(no*nv*no*nv)
      integer :: t

      do concurrent(t=1:no*nv*no*nv)
         b3(t) = b3_item(t, no, nv, t1, t2)
      end do
   end subroutine kernel_pack_b3

   subroutine kernel_woooo_base(no, nv, oooo, ooov, t1, w)
      !! w(k,l,i,j) = oooo(k,i,l,j) + sum_c [ooov(k,i,l,c) t1(j,c) + ooov(l,j,k,c) t1(i,c)]
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: oooo(no, no, no, no), ooov(no, no, no, nv), t1(no, nv)
      real(dp), intent(out) :: w(no*no*no*no)
      integer :: t

      do concurrent(t=1:no*no*no*no)
         w(t) = woooo_item(t, no, nv, oooo, ooov, t1)
      end do
   end subroutine kernel_woooo_base

   subroutine kernel_voov_base(no, nv, ovov, ovvv, ooov, t1, w)
      !! w(a,k,i,c) = ovov(k,c,i,a) + sum_d ovvv(k,c,a,d) t1(i,d) - sum_l ooov(l,i,k,c) t1(l,a)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv), ovvv(no, nv, nv, nv), ooov(no, no, no, nv), t1(no, nv)
      real(dp), intent(out) :: w(nv*no*no*nv)
      integer :: t

      do concurrent(t=1:nv*no*no*nv)
         w(t) = voov_item(t, no, nv, ovov, ovvv, ooov, t1)
      end do
   end subroutine kernel_voov_base

   subroutine kernel_vovo_base(no, nv, oovv, ovvv, ooov, t1, w)
      !! w(a,k,c,i) = oovv(k,i,a,c) + sum_d ovvv(k,d,a,c) t1(i,d) - sum_l ooov(k,i,l,c) t1(l,a)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: oovv(no, no, nv, nv), ovvv(no, nv, nv, nv), ooov(no, no, no, nv), t1(no, nv)
      real(dp), intent(out) :: w(nv*no*nv*no)
      integer :: t

      do concurrent(t=1:nv*no*nv*no)
         w(t) = vovo_item(t, no, nv, oovv, ovvv, ooov, t1)
      end do
   end subroutine kernel_vovo_base

   subroutine kernel_acc_voov(no, nv, r, w)
      !! w(a,k,i,c) += R(k,c,i,a)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: r(no*nv*no*nv)
      real(dp), intent(inout) :: w(nv*no*no*nv)
      integer :: t

      do concurrent(t=1:nv*no*no*nv)
         w(t) = w(t) + acc_voov_item(t, no, nv, r)
      end do
   end subroutine kernel_acc_voov

   subroutine kernel_acc_vovo(no, nv, r, w)
      !! w(a,k,c,i) += R(k,c,i,a)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: r(no*nv*no*nv)
      real(dp), intent(inout) :: w(nv*no*nv*no)
      integer :: t

      do concurrent(t=1:nv*no*nv*no)
         w(t) = w(t) + acc_vovo_item(t, no, nv, r)
      end do
   end subroutine kernel_acc_vovo

   ! ---------------------------------------------------------------------
   ! the item routines: pure, `routine seq`, explicit-shape dummies, no
   ! arrays of their own. A flat 1-based index t is the position in the
   ! DESTINATION, first index fastest.
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

   pure function akl_item(t, no, nv, ovov) result(v)
      !! akl(k,l,c,d) = ovov(k,c,l,d), position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp) :: v
      integer :: k, l, c, d

      call split4(t, no, no, nv, k, l, c, d)
      v = ovov(k, c, l, d)
   end function akl_item

   pure function a2_item(t, no, nv, ovov) result(v)
      !! a2(k,c,l,d) = ovov(l,c,k,d), position t of (no, nv, no, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp) :: v
      integer :: k, c, l, d

      call split4(t, no, nv, no, k, c, l, d)
      v = ovov(l, c, k, d)
   end function a2_item

   pure function b1_item(t, no, nv, t1, t2) result(v)
      !! b1(l,d,i,a), position t of (no, nv, no, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv)
      real(dp) :: v
      integer :: l, d, i, a

      call split4(t, no, nv, no, l, d, i, a)
      v = -0.5_dp*t2(i, l, d, a) - t1(i, d)*t1(l, a) + t2(i, l, a, d)
   end function b1_item

   pure function b2_item(t, no, nv, t2) result(v)
      !! b2(l,d,i,a), position t of (no, nv, no, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: t2(no, no, nv, nv)
      real(dp) :: v
      integer :: l, d, i, a

      call split4(t, no, nv, no, l, d, i, a)
      v = -0.5_dp*t2(i, l, a, d)
   end function b2_item

   pure function b3_item(t, no, nv, t1, t2) result(v)
      !! b3(l,d,i,a), position t of (no, nv, no, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv)
      real(dp) :: v
      integer :: l, d, i, a

      call split4(t, no, nv, no, l, d, i, a)
      v = -0.5_dp*t2(i, l, d, a) - t1(i, d)*t1(l, a)
   end function b3_item

   pure function woooo_item(t, no, nv, oooo, ooov, t1) result(v)
      !! The tau-free part of Woooo(k,l,i,j) at position t of (no, no, no, no); c ascending
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: oooo(no, no, no, no), ooov(no, no, no, nv), t1(no, nv)
      real(dp) :: v
      integer :: k, l, i, j, c

      call split4(t, no, no, no, k, l, i, j)
      v = oooo(k, i, l, j)
      do c = 1, nv
         v = v + ooov(k, i, l, c)*t1(j, c) + ooov(l, j, k, c)*t1(i, c)
      end do
   end function woooo_item

   pure function voov_item(t, no, nv, ovov, ovvv, ooov, t1) result(v)
      !! The GEMM-free part of Wvoov(a,k,i,c) at position t of (nv, no, no, nv); d, then l, ascending
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv), ovvv(no, nv, nv, nv), ooov(no, no, no, nv), t1(no, nv)
      real(dp) :: v
      integer :: a, k, i, c, d, l

      call split4(t, nv, no, no, a, k, i, c)
      v = ovov(k, c, i, a)
      do d = 1, nv
         v = v + ovvv(k, c, a, d)*t1(i, d)
      end do
      do l = 1, no
         v = v - ooov(l, i, k, c)*t1(l, a)
      end do
   end function voov_item

   pure function vovo_item(t, no, nv, oovv, ovvv, ooov, t1) result(v)
      !! The GEMM-free part of Wvovo(a,k,c,i) at position t of (nv, no, nv, no); d, then l, ascending
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: oovv(no, no, nv, nv), ovvv(no, nv, nv, nv), ooov(no, no, no, nv), t1(no, nv)
      real(dp) :: v
      integer :: a, k, c, i, d, l

      call split4(t, nv, no, nv, a, k, c, i)
      v = oovv(k, i, a, c)
      do d = 1, nv
         v = v + ovvv(k, d, a, c)*t1(i, d)
      end do
      do l = 1, no
         v = v - ooov(k, i, l, c)*t1(l, a)
      end do
   end function vovo_item

   pure function acc_voov_item(t, no, nv, r) result(v)
      !! R(k,c,i,a) for the element Wvoov(a,k,i,c) at position t of (nv, no, no, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: r(no*nv, no*nv)
      real(dp) :: v
      integer :: a, k, i, c

      call split4(t, nv, no, no, a, k, i, c)
      v = r((c - 1)*no + k, (a - 1)*no + i)
   end function acc_voov_item

   pure function acc_vovo_item(t, no, nv, r) result(v)
      !! R(k,c,i,a) for the element Wvovo(a,k,c,i) at position t of (nv, no, nv, no)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: r(no*nv, no*nv)
      real(dp) :: v
      integer :: a, k, c, i

      call split4(t, nv, no, nv, a, k, c, i)
      v = r((c - 1)*no + k, (a - 1)*no + i)
   end function acc_vovo_item

end module trc_rccsd_w
