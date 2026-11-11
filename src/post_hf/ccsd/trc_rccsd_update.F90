!
! The amplitude update of a closed-shell RI-CCSD(T) calculation: one function of
! the amplitudes, (t1', t2') = update(t1, t2), for canonical RHF orbitals.
!
! Hirata, Podeszwa, Tobita and Bartlett, JCP 120, 2581 (2004), Eqs. (35) and (36),
! as pyscf/cc/rccsd.py `update_amps` writes them (the non-CC2 branch), ported from
! the same author's mqc_czt_rcc.f90 (`rccsd_iteration`, `add_symmetrised`,
! `scatter_ring`). Chemists' notation, spatial orbitals, REAL orbitals; the
! integral blocks, tau, the one- and two-particle intermediates and the ladder are
! those of trc_rccsd_ints, _inter, _w and _ladder, which this module drives.
!
! CANONICAL ORBITALS, f_ov = 0, the Fock diagonal is e_p. PySCF adds the bare
! diag(e) to Foo, Fvv, Loo, Lvv and subtracts it straight back out of the same
! arrays before using them; the intermediates here never carry it (trc_rccsd_inter),
! so the subtraction is already done. The denominators are the same diagonal, so a
! level shift would have to enter both.
!
! ALL INTERMEDIATES FROM THE SAME INPUT. `step` builds tau, the one-particle
! intermediates, the W and the ladder from `amps`, and only then forms the new
! amplitudes, into the SEPARATE arrays `newamps`. Nothing it computes reads
! `newamps`, and `amps` is never written.
!
! THE T1 EQUATION (Eq. 35; i,k,l occupied, a,c,d virtual; Foo, Fvv, Fov are
! the cc_F* of the inter module). Each line: PySCF einsum, then how it is computed.
!
!   'ac,ic->ia'  Fvv t1          t1n(i,a) +=  sum_c Fvv(a,c) t1(i,c)
!   '-ki,ka->ia' Foo t1          t1n(i,a) -=  sum_k Foo(k,i) t1(k,a)
!   '2 kc,kica->ia' - 'kc,ikca->ia' + 'kc,ic,ka->ia'          (fov = 0, so no -2 kc,ka,ic term)
!                                t1n(i,a) +=  sum_kc Fov(k,c) [2 t2(k,i,c,a) - t2(i,k,c,a) + t1(i,c) t1(k,a)]
!   '2 kcai,kc->ia' - 'kiac,kc->ia'
!                                t1n(i,a) +=  sum_kc [2 ovov(k,c,i,a) - oovv(k,i,a,c)] t1(k,c)
!     the four lines above are ONE kernel, one thread per (i,a), o v items of o v MACs
!   '2 kdac,ikcd->ia' - 'kcad,ikcd->ia' (+ the same pair with t2 -> t1 t1, together tau)
!                                t1n(i,a) +=  sum_kcd [2 ovvv(k,d,a,c) - ovvv(k,c,a,d)] tau(i,k,c,d)
!     for each virtual x one GEMM (m = o, n = v, K = o v, 2 o^2 v^2 MACs each, v of them), because for a fixed x
!     the block ovvv(:,:,:,x) IS the matrix (k d; a); the two ovvv terms share one
!     packed operand, y(k,x',i,x) = 2 tau(i,k,x,x') - tau(i,k,x',x), so it is one GEMM per x
!   '-2 lcki,klac->ia' + 'kcli,klac->ia' (+ the same pair with t1 t1, together tau)
!                                t1n(i,a) -=  sum_klc lam(k,i,l,c) tau(k,l,a,c),  lam = 2 ooov(k,i,l,c) - ooov(l,i,k,c)
!     one GEMM, m = o, n = v, K = o^2 v, from the packs lamp(k,l,c; i) and taur(k,l,c; a)
!   t1n(i,a) /= e_i - e_a
!
! THE T2 EQUATION (Eq. 36), in mqc's order of accumulation into t2n:
!
!   'ovov.transpose(0,2,1,3)'    t2n(i,j,a,b)  = ovov(i,a,j,b)                                     kernel
!   A: term = tmp2 . t1, t2n += term + term.T      tmp2(a,b,i,c) = ovvv(i,a,c,b) - sum_k oovv(k,i,b,c) t1(k,a)
!        term(i,j,a,b) = sum_c ovvv(i,a,c,b) t1(j,c) - sum_k t1(k,a) [sum_c oovv(k,i,b,c) t1(j,c)]
!        'abic,jc->ijab': R(i,a,j,b) = sum_c ovvv(i,a,c,b) t1(j,c) is a strided batched GEMM over b (A_b is
!        the matrix (i a; c) of ovvv itself, m = o v, n = o, K = v) and a permutation; the oovv half is
!        U(k,i,j,b) = sum_c oovv(k,i,b,c) t1(j,c), a batch over b (m = o^2, n = o, K = v), then
!        term(ij;a,b) -= U_b(k;ij)^T t1(k;a), a batch over b (m = o^2, n = v, K = o). The regrouping of
!        mqc's tmp2 (an o v^3 array) changes the order of the sums, not the sums
!   B: term = tmp2b . t1, t2n -= term + term.T     tmp2b(a,k,i,j) = sum_c ovov(k,c,i,a) t1(j,c) + ooov(j,k,i,a)
!        'kcai,jc->akij' + 'ovoo.transpose(1,3,0,2)', then 'akij,kb->ijab'.
!        Y(k,j,i,a) = sum_c ovov(k,c,i,a) t1(j,c) is a batch over (i,a) (m = n = o, K = v), packed with
!        ooov into x(k,i,j,a) = tmp2b(a,k,i,j), and term(ij a; b) = x(k; ija)^T t1(k;b) is ONE GEMM
!        (m = o^2 v, n = v, K = o)
!   'klij,klab->ijab'            t2n(ij;ab) += Woooo(kl;ij)^T tau(kl;ab)              ONE GEMM (m = o^2, n = v^2, K = o^2)
!   'abcd,ijcd->ijab'            t2n += lad%h        the ladder and its t1 dressing (cc_Wvvvv . tau), unsymmetrised,
!                                coefficient +1, as mqc adds it; trc_rccsd_ladder's GEMMs
!   'ac,ijcb->ijab'              term(ij,a,b) = sum_c Lvv(a,c) t2(ij,c,b); t2n += term + term.T
!                                a batch over b, m = o^2, n = v, K = v
!   '-ki,kjab->ijab'             term(i;jab) = sum_k Loo(k,i) t2(k;jab); t2n -= term + term.T
!                                ONE GEMM (m = o, n = o v^2, K = o)
!   the three ring terms, each ONE GEMM m = n = K = o v (O(o^3 v^3)), as mqc packs them:
!        '2 akic,kjcb->ijab' - 'akci,kjcb->ijab'   R1 = (2 Wvoov - Wvovo)(ai;kc) . t2(kc;jb)   y1(k,c,j,b) = t2(k,j,c,b)
!        '-akic,kjbc->ijab'                        R2 = Wvoov(ai;kc) . y2(kc;jb)              y2(k,c,j,b) = t2(k,j,b,c)
!        '-bkci,kjac->ijab'                        R3 = Wvovo(ai;kc) . y2(kc;jb), scattered with (b,i) against (j,a)
!        x(a,i,k,c): x2 = Wvoov(a,k,i,c), x3 = Wvovo(a,k,c,i), x1 = 2 x2 - x3. Each R is scattered to (i,j,a,b) and added
!        symmetrised in ONE kernel (the scatter and the symmetrisation of mqc fused: no intermediate array)
!   t2n(i,j,a,b) /= e_i + e_j - e_a - e_b
! with `tmp + tmp.T` meaning tmp(i,j,a,b) + tmp(j,i,b,a) (`term` in the code, `tmp` in PySCF and mqc), the one
! symmetrisation of every PySCF line, written once in `kernel_add_sym`.
!
! Every O(N^5) and O(N^6) contraction above is a GEMM; the kernels that are not (the T1 small terms, the
! packs and scatters, the denominators) are O(o^2 v^2) per call except the T1 small kernel, whose o^2 v^2 MACs are
! the cost of one (kc) sum per element.
!
! BUFFERS. Nothing the size of o^2 v^2 is allocated here: `step` borrows the scratch of the modules it drives,
! each free once the module that owns it has returned. o = active occupied, v = virtual; 8 bytes a double, host and
! device each; "owner" is the type that allocates and releases it.
!
!   ws%b1  o^2 v^2   owner trc_rccsd_wscratch_t (W's b1). Free after rccsd_w_build. Here: the ycomb pack of the T1 ovvv
!                    term, then R(i,a,j,b) of term A, then the x pack of the rings. Written before it is read each time.
!   ws%b2  o^2 v^2   owner trc_rccsd_wscratch_t (W's b2). Free after rccsd_w_build. Here: taur of the T1 ooov term,
!                    then the y pack of the rings.
!   ws%r   o^2 v^2   owner trc_rccsd_wscratch_t (W's r). Free after rccsd_w_build. Here: the ring GEMM result R.
!   scr%pack1 o^2 v^2 owner trc_rccsd_scratch_t. Per-call scratch of Foo/Fvv, free after rccsd_inter_one_particle.
!                    Here: `term` of the symmetrised terms.
!   sm1, sm2 o^3 v   owner this type. U (term A) then Y (term B); x of term B.
!   lamp     o^3 v   owner this type, built once by `init`: lam(k,i,l,c) as lamp(k,l,c,i); stamped like the packs of
!                    the modules above.
!   lad%h    o^2 v^2 owner trc_rccsd_ladder_t. The ladder, read once by `kernel_add_ladder`.
! so this module adds 3 o^3 v doubles. The ladder's own work arrays (`ls`) are not used after rccsd_ladder_build.
!
! RESIDENCY. Every array this module allocates exists on the host AND on the device (`enter data` in `init`,
! `exit data` in `release`, on the same `this%name` expression). The device copy is the one the kernels write: `step`
! leaves `newamps` on the device (pull it to read the host copy). A kernel gets an array only as an explicit-shape
! dummy argument of a `kernel_*` routine and no `x%name` appears in a `do concurrent` body.
!
! STALE PACKS. `init` records `ints%stamp`; `step` refuses (status, nothing written) after `ints%build` ran again,
! as the modules it drives do for their own packs. Re-`init`.
!
! `newamps` must not be `amps` (no check can see an alias; the result would be wrong). Single rank only.
module trc_rccsd_update
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t, trc_rccsd_scratch_t, rccsd_inter_one_particle, &
                              rccsd_inter_flops
   use trc_rccsd_w, only: trc_rccsd_w_t, trc_rccsd_wscratch_t, rccsd_w_build, rccsd_w_flops
   use trc_rccsd_ladder, only: trc_rccsd_ladder_t, trc_rccsd_lscratch_t, trc_rccsd_ladder_options_t, rccsd_ladder_build
   use trc_rccsd_steps, only: RCCSD_ST_INTER, RCCSD_ST_W, RCCSD_ST_T1, RCCSD_ST_T2, RCCSD_ST_RING, rccsd_steps_gemm_flops, &
                              rccsd_steps_open, rccsd_steps_close
   use trc_step_table, only: step_table_t
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   character(len=*), parameter :: STALE_MSG = "trc_rccsd_update: the integrals were rebuilt after init; re-init"

   public :: trc_rccsd_update_t, rccsd_update_t1_flops, rccsd_update_t2_flops, rccsd_update_ring_flops

   type :: trc_rccsd_update_t
      !! Everything one amplitude update needs besides the integrals and the amplitudes
      type(error_t) :: status  !! set by `init` and `step`: a refusal leaves its message here and writes nothing
      integer :: no = 0
      integer :: nv = 0
      integer :: stamp = 0  !! `ints%stamp` of the integrals `init` packed from
      logical :: tau_symmetric = .true.
      !! `step` tells the ladder that tau(j,i,d,c) = tau(i,j,c,d). True for the MP2 start and everything the
      !! iteration makes of it; the driver clears it for a run that starts from amplitudes without that symmetry
      !! (the half ladder would otherwise return the ladder of the symmetrised tau)
      type(trc_rccsd_inter_t) :: inter     !! tau, Foo, Fvv, Fov, Loo, Lvv of the last `step`
      type(trc_rccsd_scratch_t) :: scr
      type(trc_rccsd_w_t) :: w             !! Woooo, Wvoov, Wvovo of the last `step`
      type(trc_rccsd_wscratch_t) :: ws
      type(trc_rccsd_ladder_t) :: lad      !! the ladder of the last `step`
      type(trc_rccsd_lscratch_t) :: ls
      real(dp), allocatable :: lamp(:)  !! o^3 v, lamp(k,l,c,i) = 2 ooov(k,i,l,c) - ooov(l,i,k,c)
      real(dp), allocatable :: sm1(:)   !! o^3 v work
      real(dp), allocatable :: sm2(:)   !! o^3 v work
   contains
      procedure :: init => update_init
      procedure :: release => update_release
      procedure :: step => update_step
   end type trc_rccsd_update_t

contains

   ! ---------------------------------------------------------------------
   ! the type
   ! ---------------------------------------------------------------------

   subroutine update_init(this, ints, ladder_options)
      !! Allocates and packs everything for these integrals. A refusal of any module goes to `status`
      !! and leaves the type empty.
      class(trc_rccsd_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_ladder_options_t), intent(in), optional :: ladder_options

      integer :: no, nv

      call this%release()
      call this%status%clear()
      if (.not. allocated(ints%ovvv) .or. ints%stamp == 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_update: the integral blocks are not built")
         return
      end if
      no = ints%no
      nv = ints%nv
      call this%scr%init(ints)
      if (this%scr%status%has_error()) call fail(this, this%scr%status)
      if (this%status%has_error()) return
      call this%ws%init(ints)
      if (this%ws%status%has_error()) call fail(this, this%ws%status)
      if (this%status%has_error()) return
      call this%w%init(no, nv)
      if (this%w%status%has_error()) call fail(this, this%w%status)
      if (this%status%has_error()) return
      call this%lad%init(no, nv)
      if (this%lad%status%has_error()) call fail(this, this%lad%status)
      if (this%status%has_error()) return
      call this%ls%init(ints, ladder_options)
      if (this%ls%status%has_error()) call fail(this, this%ls%status)
      if (this%status%has_error()) return
      call this%inter%init(no, nv)
      this%no = no
      this%nv = nv
      this%stamp = ints%stamp
      allocate (this%lamp(no*no*no*nv), this%sm1(no*no*no*nv), this%sm2(no*no*no*nv))
      !$acc enter data create(this%lamp, this%sm1, this%sm2)
      call kernel_pack_lamp(no, nv, this%scr%lam, this%lamp)
   end subroutine update_init

   subroutine fail(this, source)
      !! Copies a refusal into `this%status` and frees what `init` had allocated so far
      class(trc_rccsd_update_t), intent(inout) :: this
      type(error_t), intent(in) :: source

      character(len=:), allocatable :: text

      text = source%get_message()
      call this%release()
      call this%status%set(ERROR_VALIDATION, text)
   end subroutine fail

   subroutine update_release(this)
      !! Frees host and device copies of everything; safe on an empty type
      class(trc_rccsd_update_t), intent(inout) :: this

      if (allocated(this%lamp)) then
         !$acc exit data delete(this%lamp)
         deallocate (this%lamp)
      end if
      if (allocated(this%sm1)) then
         !$acc exit data delete(this%sm1)
         deallocate (this%sm1)
      end if
      if (allocated(this%sm2)) then
         !$acc exit data delete(this%sm2)
         deallocate (this%sm2)
      end if
      call this%inter%release()
      call this%scr%release()
      call this%w%release()
      call this%ws%release()
      call this%lad%release()
      call this%ls%release()
      this%no = 0
      this%nv = 0
      this%stamp = 0
   end subroutine update_release

   function step_ready(this, ints, amps, newamps) result(ok)
      !! True when `init` ran for these integrals and the amplitudes have their sizes; else `this%status` says which not
      class(trc_rccsd_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, newamps
      logical :: ok

      ok = .false.
      if (this%stamp == 0 .or. .not. allocated(this%lamp)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_update: step before init")
         return
      end if
      if (.not. allocated(ints%ovvv) .or. ints%stamp == 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_update: the integral blocks are not built")
         return
      end if
      if (ints%stamp /= this%stamp) then
         call this%status%set(ERROR_VALIDATION, STALE_MSG)
         return
      end if
      if (.not. allocated(amps%t1) .or. .not. allocated(amps%t2) .or. .not. allocated(newamps%t1) .or. &
          .not. allocated(newamps%t2)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_update: the amplitudes are not initialised")
         return
      end if
      if (amps%no /= ints%no .or. amps%nv /= ints%nv .or. newamps%no /= ints%no .or. newamps%nv /= ints%nv .or. &
          size(amps%t2, kind=int64) /= int(ints%no, int64)*ints%no*ints%nv*ints%nv .or. &
          size(newamps%t2, kind=int64) /= int(ints%no, int64)*ints%no*ints%nv*ints%nv .or. &
          size(amps%t1, kind=int64) /= int(ints%no, int64)*ints%nv .or. &
          size(newamps%t1, kind=int64) /= int(ints%no, int64)*ints%nv) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_update: (no, nv) differ between integrals and amplitudes")
         return
      end if
      ok = .true.
   end function step_ready

   function sub_failed(this) result(failed)
      !! True, with `this%status` set, when a module `step` drives reported a refusal
      class(trc_rccsd_update_t), intent(inout) :: this
      logical :: failed

      failed = .true.
      if (this%scr%status%has_error()) then
         call this%status%set(ERROR_VALIDATION, this%scr%status%get_message())
      else if (this%ws%status%has_error()) then
         call this%status%set(ERROR_VALIDATION, this%ws%status%get_message())
      else if (this%ls%status%has_error()) then
         call this%status%set(ERROR_VALIDATION, this%ls%status%get_message())
      else
         failed = .false.
      end if
   end function sub_failed

   ! ---------------------------------------------------------------------
   ! the update
   ! ---------------------------------------------------------------------

   pure function rccsd_update_t1_flops(no, nv) result(flops)
      !! Flops of the GEMMs of `t1_equation` (the small T1 kernel, the packs and the denominator are not counted):
      !!   ovvv term, one (o, v, o v) per virtual   2 o^2 v^3         lam term (o, v, o o v)   2 o^3 v^2
      integer, intent(in) :: no, nv
      integer(int64) :: flops

      flops = rccsd_steps_gemm_flops(no, nv, no*nv, nv) + rccsd_steps_gemm_flops(no, nv, no*no*nv)
   end function rccsd_update_t1_flops

   pure function rccsd_update_t2_flops(no, nv) result(flops)
      !! Flops of the GEMMs of the T2 equation without its rings: `t2_base_terms`, `t2_two_body_terms`,
      !! `t2_one_body_terms` (the ladder is `rccsd_ladder_flops`'s). Shapes (m, n, k) x batch count:
      !!   A: ovvv t1 (o v, o, v) x v          2 o^2 v^3       U: oovv t1 (o o, o, v) x v        2 o^3 v^2
      !!   A: term (o o, v, o) x v             2 o^3 v^2       B: ovov t1 (o, o, v) x o v        2 o^3 v^2
      !!   B: (o o v, v, o)                    2 o^3 v^2       Woooo tau (o o, v v, o o)         2 o^4 v^2
      !!   Lvv: (o o, v, v) x v                2 o^2 v^3       Loo: (o, o v v, o)                2 o^3 v^2
      integer, intent(in) :: no, nv
      integer(int64) :: flops

      flops = rccsd_steps_gemm_flops(no*nv, no, nv, nv) + rccsd_steps_gemm_flops(no*no, no, nv, nv) + &
              rccsd_steps_gemm_flops(no*no, nv, no, nv) + rccsd_steps_gemm_flops(no, no, nv, no*nv) + &
              rccsd_steps_gemm_flops(no*no*nv, nv, no) + rccsd_steps_gemm_flops(no*no, nv*nv, no*no) + &
              rccsd_steps_gemm_flops(no*no, nv, nv, nv) + rccsd_steps_gemm_flops(no, no*nv*nv, no)
   end function rccsd_update_t2_flops

   pure function rccsd_update_ring_flops(no, nv) result(flops)
      !! Flops of the three ring GEMMs of `t2_ring_terms`, each (o v, o v, o v): 6 (o v)^3
      integer, intent(in) :: no, nv
      integer(int64) :: flops

      flops = 3_int64*rccsd_steps_gemm_flops(no*nv, no*nv, no*nv)
   end function rccsd_update_ring_flops

   subroutine update_step(this, ints, amps, newamps, la, stab)
      !! (newamps%t1, newamps%t2) = update(amps%t1, amps%t2), on the device. `amps` must be on the device
      !! (`push` after a host edit) and is not written. A refusal is in `this%status` and nothing is written.
      !! With `stab` every stage is a step of that table (trc_rccsd_steps), timed at the same device syncs the
      !! stages always had, and its flops are added; without it nothing is timed or counted.
      class(trc_rccsd_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_amps_t), intent(inout) :: newamps
      type(trc_linalg_t), intent(in) :: la
      type(step_table_t), intent(inout), optional :: stab

      type(timer_type) :: tc
      integer :: no, nv

      call this%status%clear()
      call this%scr%status%clear()
      call this%ws%status%clear()
      call this%ls%status%clear()
      if (.not. step_ready(this, ints, amps, newamps)) return
      no = this%no
      nv = this%nv

      ! The three builders read `amps` and write only their own arrays; none of them is told about `newamps`.
      call rccsd_steps_open(stab, la, RCCSD_ST_INTER, tc)
      call rccsd_inter_one_particle(ints, amps, this%inter, this%scr, la)
      if (sub_failed(this)) then
         call rccsd_steps_close(stab, la, RCCSD_ST_INTER, tc)  ! a refusal still closes the NVTX range
         return
      end if
      call rccsd_steps_close(stab, la, RCCSD_ST_INTER, tc, rccsd_inter_flops(no, nv))
      call rccsd_steps_open(stab, la, RCCSD_ST_W, tc)
      call rccsd_w_build(ints, amps, this%inter, this%w, this%ws, la)
      if (sub_failed(this)) then
         call rccsd_steps_close(stab, la, RCCSD_ST_W, tc)  ! a refusal still closes the NVTX range
         return
      end if
      call rccsd_steps_close(stab, la, RCCSD_ST_W, tc, rccsd_w_flops(no, nv))
      ! the ladder and its dressing are two steps of their own, opened and closed inside. tau = t2 + t1 t1 has
      ! tau(j,i,d,c) = tau(i,j,c,d): the starting t2 (MP2) has it, every term of the t2 equation keeps it and DIIS
      ! mixes such t2 linearly. Start amplitudes of the caller's may not: the driver checks them and clears the flag.
      call rccsd_ladder_build(ints, amps, this%inter, this%lad, this%ls, la, stab, tau_symmetric=this%tau_symmetric)
      if (sub_failed(this)) return

      ! From here the W scratch (b1, b2, r) and the inter scratch's pack1 are free.
      call rccsd_steps_open(stab, la, RCCSD_ST_T1, tc)
      call t1_equation(this, ints, amps, newamps, la)
      call rccsd_steps_close(stab, la, RCCSD_ST_T1, tc, rccsd_update_t1_flops(no, nv))
      call rccsd_steps_open(stab, la, RCCSD_ST_T2, tc)
      call t2_base_terms(this, ints, amps, newamps, la)
      call t2_two_body_terms(this, newamps, la)
      call t2_one_body_terms(this, amps, newamps, la)
      call rccsd_steps_close(stab, la, RCCSD_ST_T2, tc, rccsd_update_t2_flops(no, nv))
      call rccsd_steps_open(stab, la, RCCSD_ST_RING, tc)
      call t2_ring_terms(this, amps, newamps, la)
      call rccsd_steps_close(stab, la, RCCSD_ST_RING, tc, rccsd_update_ring_flops(no, nv))
      ! the denominators divide the finished sum, so they follow the rings; they belong to the T2 equation
      call rccsd_steps_open(stab, la, RCCSD_ST_T2, tc)
      call kernel_denominator_t2(no, nv, ints%eps_o, ints%eps_v, newamps%t2)
      call rccsd_steps_close(stab, la, RCCSD_ST_T2, tc)
   end subroutine update_step

   subroutine t1_equation(this, ints, amps, newamps, la)
      !! Eq. (35), the lines of the module header, into `newamps%t1`
      class(trc_rccsd_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_amps_t), intent(inout) :: newamps
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, nov, x

      no = this%no
      nv = this%nv
      nov = no*nv
      call kernel_t1_small(no, nv, amps%t1, amps%t2, this%inter%fvv, this%inter%foo, this%inter%fov, &
                           ints%ovov, ints%oovv, newamps%t1)
      ! 2 ovvv(k,d,a,c) tau(i,k,c,d) - ovvv(k,c,a,d) tau(i,k,c,d): per virtual x, C(i,a) += y_x(kd;i)^T ovvv_x(kd;a)
      call kernel_pack_ycomb(no, nv, this%inter%tau, this%ws%b1)
      do x = 1, nv
         call la%gemm("T", "N", no, nv, nov, 1.0_dp, this%ws%b1(1 + (x - 1)*nov*no), nov, &
                      ints%ovvv(1, 1, 1, x), nov, 1.0_dp, newamps%t1, no)
      end do
      ! - lam(k,i,l,c) tau(k,l,a,c)
      call kernel_pack_taur(no, nv, this%inter%tau, this%ws%b2)
      call la%gemm("T", "N", no, nv, no*no*nv, -1.0_dp, this%lamp, no*no*nv, this%ws%b2, no*no*nv, 1.0_dp, &
                   newamps%t1, no)
      call kernel_denominator_t1(no, nv, ints%eps_o, ints%eps_v, newamps%t1)
   end subroutine t1_equation

   subroutine t2_base_terms(this, ints, amps, newamps, la)
      !! ovov, then the t1-dressed ovvv/oovv term A and the ovov/ooov term B (each symmetrised)
      class(trc_rccsd_update_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_amps_t), intent(inout) :: newamps
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: shape
      integer :: no, nv, nov, no2

      no = this%no
      nv = this%nv
      nov = no*nv
      no2 = no*no
      ! t2new = ovov.transpose(0,2,1,3)
      call kernel_t2_base(no, nv, ints%ovov, newamps%t2)

      ! A. R(i,a,j,b) = sum_c ovvv(i,a,c,b) t1(j,c): over b, A_b = ovvv(:,:,:,b) as the matrix (i a; c)
      shape = trc_gemm_batch(nov, no, nv, nv, lda=nov, ldb=no, ldc=nov, &
                             sa=int(nov, int64)*nv, sb=0_int64, sc=int(nov, int64)*no, transa="N", transb="T")
      call la%gemm_batched(shape, ints%ovvv, amps%t1, this%ws%b1)
      call kernel_perm_r1(no, nv, this%ws%b1, this%scr%pack1)
      ! U(k,i,j,b) = sum_c oovv(k,i,b,c) t1(j,c), over b
      shape = trc_gemm_batch(no2, no, nv, nv, lda=no2*nv, ldb=no, ldc=no2, &
                             sa=int(no2, int64), sb=0_int64, sc=int(no2, int64)*no, transa="N", transb="T")
      call la%gemm_batched(shape, ints%oovv, amps%t1, this%sm1)
      ! term(ij,a,b) -= sum_k U(k,ij,b) t1(k,a), over b
      shape = trc_gemm_batch(no2, nv, no, nv, lda=no, ldb=no, ldc=no2, &
                             sa=int(no2, int64)*no, sb=0_int64, sc=int(no2, int64)*nv, transa="T", transb="N")
      call la%gemm_batched(shape, this%sm1, amps%t1, this%scr%pack1, -1.0_dp, 1.0_dp)
      call kernel_add_sym(no, nv, 1.0_dp, this%scr%pack1, newamps%t2)

      ! B. Y(k,j,i,a) = sum_c ovov(k,c,i,a) t1(j,c), over (i,a)
      shape = trc_gemm_batch(no, no, nv, nov, lda=no, ldb=no, ldc=no, &
                             sa=int(nov, int64), sb=0_int64, sc=int(no2, int64), transa="N", transb="T")
      call la%gemm_batched(shape, ints%ovov, amps%t1, this%sm1)
      call kernel_pack_x_b(no, nv, this%sm1, ints%ooov, this%sm2)
      call la%gemm("T", "N", no2*nv, nv, no, 1.0_dp, this%sm2, no, amps%t1, no, 0.0_dp, this%scr%pack1, no2*nv)
      call kernel_add_sym(no, nv, -1.0_dp, this%scr%pack1, newamps%t2)
   end subroutine t2_base_terms

   subroutine t2_two_body_terms(this, newamps, la)
      !! 'klij,klab->ijab' and the ladder, both added without symmetrisation
      class(trc_rccsd_update_t), intent(inout) :: this
      type(trc_rccsd_amps_t), intent(inout) :: newamps
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, no2

      no = this%no
      nv = this%nv
      no2 = no*no
      ! t2n(ij;ab) += sum_kl Woooo(kl;ij) tau(kl;ab)
      call la%gemm("T", "N", no2, nv*nv, no2, 1.0_dp, this%w%woooo, no2, this%inter%tau, no2, 1.0_dp, newamps%t2, no2)
      ! 'abcd,ijcd->ijab': the ladder, h(i,j,a,b), coefficient +1, after the Woooo term and before Lvv/Loo
      call kernel_add_ladder(no, nv, this%lad%h, newamps%t2)
   end subroutine t2_two_body_terms

   subroutine t2_one_body_terms(this, amps, newamps, la)
      !! 'ac,ijcb->ijab' and '-ki,kjab->ijab', each symmetrised, from the OLD t2
      class(trc_rccsd_update_t), intent(inout) :: this
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_amps_t), intent(inout) :: newamps
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: shape
      integer :: no, nv, no2

      no = this%no
      nv = this%nv
      no2 = no*no
      ! term(ij,a,b) = sum_c t2(ij,c,b) Lvv(a,c), over b
      shape = trc_gemm_batch(no2, nv, nv, nv, lda=no2, ldb=nv, ldc=no2, &
                             sa=int(no2, int64)*nv, sb=0_int64, sc=int(no2, int64)*nv, transa="N", transb="T")
      call la%gemm_batched(shape, amps%t2, this%inter%lvv, this%scr%pack1)
      call kernel_add_sym(no, nv, 1.0_dp, this%scr%pack1, newamps%t2)
      ! term(i;jab) = sum_k Loo(k,i) t2(k;jab)
      call la%gemm("T", "N", no, no*nv*nv, no, 1.0_dp, this%inter%loo, no, amps%t2, no, 0.0_dp, this%scr%pack1, no)
      call kernel_add_sym(no, nv, -1.0_dp, this%scr%pack1, newamps%t2)
   end subroutine t2_one_body_terms

   subroutine t2_ring_terms(this, amps, newamps, la)
      !! The three ring terms, from the OLD t2, each one GEMM and one fused scatter/symmetrise kernel
      class(trc_rccsd_update_t), intent(inout) :: this
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_amps_t), intent(inout) :: newamps
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, nov

      no = this%no
      nv = this%nv
      nov = no*nv
      ! '2 akic,kjcb->ijab' - 'akci,kjcb->ijab'
      call kernel_pack_ring_x(no, nv, 1, this%w%wvoov, this%w%wvovo, this%ws%b1)
      call kernel_pack_ring_y(no, nv, 1, amps%t2, this%ws%b2)
      call la%gemm("N", "N", nov, nov, nov, 1.0_dp, this%ws%b1, nov, this%ws%b2, nov, 0.0_dp, this%ws%r, nov)
      call kernel_add_ring(no, nv, 1, 1.0_dp, this%ws%r, newamps%t2)
      ! '-akic,kjbc->ijab'
      call kernel_pack_ring_x(no, nv, 2, this%w%wvoov, this%w%wvovo, this%ws%b1)
      call kernel_pack_ring_y(no, nv, 2, amps%t2, this%ws%b2)
      call la%gemm("N", "N", nov, nov, nov, 1.0_dp, this%ws%b1, nov, this%ws%b2, nov, 0.0_dp, this%ws%r, nov)
      call kernel_add_ring(no, nv, 1, -1.0_dp, this%ws%r, newamps%t2)
      ! '-bkci,kjac->ijab': the same y2, and (b,i) against (j,a)
      call kernel_pack_ring_x(no, nv, 3, this%w%wvoov, this%w%wvovo, this%ws%b1)
      call la%gemm("N", "N", nov, nov, nov, 1.0_dp, this%ws%b1, nov, this%ws%b2, nov, 0.0_dp, this%ws%r, nov)
      call kernel_add_ring(no, nv, 2, -1.0_dp, this%ws%r, newamps%t2)
   end subroutine t2_ring_terms

   ! ---------------------------------------------------------------------
   ! the kernels. Each takes its arrays as explicit-shape dummies and runs one
   ! `do concurrent` over ONE flattened index, whose body calls an item routine
   ! of this module. No `x%name` appears in any of them.
   ! ---------------------------------------------------------------------

   subroutine kernel_pack_lamp(no, nv, lam, lamp)
      !! lamp(k,l,c,i) = lam(k,i,l,c)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: lam(no, no, no, nv)
      real(dp), intent(out) :: lamp(no*no*nv*no)
      integer :: t

      do concurrent(t=1:no*no*nv*no)
         lamp(t) = lamp_item(t, no, nv, lam)
      end do
   end subroutine kernel_pack_lamp

   subroutine kernel_t1_small(no, nv, t1, t2, fvv, foo, fov, ovov, oovv, t1n)
      !! t1n(i,a) = the four small T1 lines of the module header, one thread per (i,a)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv), fvv(nv, nv), foo(no, no), fov(no, nv)
      real(dp), intent(in) :: ovov(no, nv, no, nv), oovv(no, no, nv, nv)
      real(dp), intent(out) :: t1n(no*nv)
      integer :: t

      do concurrent(t=1:no*nv)
         t1n(t) = t1_small_item(t, no, nv, t1, t2, fvv, foo, fov, ovov, oovv)
      end do
   end subroutine kernel_t1_small

   subroutine kernel_pack_ycomb(no, nv, tau, y)
      !! y(k,x',i,x) = 2 tau(i,k,x,x') - tau(i,k,x',x)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(out) :: y(no*nv*no*nv)
      integer :: t

      do concurrent(t=1:no*nv*no*nv)
         y(t) = ycomb_item(t, no, nv, tau)
      end do
   end subroutine kernel_pack_ycomb

   subroutine kernel_pack_taur(no, nv, tau, taur)
      !! taur(k,l,c,a) = tau(k,l,a,c)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(out) :: taur(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         taur(t) = taur_item(t, no, nv, tau)
      end do
   end subroutine kernel_pack_taur

   subroutine kernel_denominator_t1(no, nv, eps_o, eps_v, t1n)
      !! t1n(i,a) /= e_i - e_a
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: eps_o(no), eps_v(nv)
      real(dp), intent(inout) :: t1n(no*nv)
      integer :: t

      do concurrent(t=1:no*nv)
         t1n(t) = t1n(t)/(eps_o(mod(t - 1, no) + 1) - eps_v((t - 1)/no + 1))
      end do
   end subroutine kernel_denominator_t1

   subroutine kernel_t2_base(no, nv, ovov, t2n)
      !! t2n(i,j,a,b) = ovov(i,a,j,b)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp), intent(out) :: t2n(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         t2n(t) = base_item(t, no, nv, ovov)
      end do
   end subroutine kernel_t2_base

   subroutine kernel_perm_r1(no, nv, r1, term)
      !! term(i,j,a,b) = r1(i,a,j,b)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: r1(no, nv, no, nv)
      real(dp), intent(out) :: term(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         term(t) = perm_r1_item(t, no, nv, r1)
      end do
   end subroutine kernel_perm_r1

   subroutine kernel_pack_x_b(no, nv, y, ooov, x)
      !! x(k,i,j,a) = y(k,j,i,a) + ooov(j,k,i,a)     (tmp2b(a,k,i,j) of mqc)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: y(no, no, no, nv), ooov(no, no, no, nv)
      real(dp), intent(out) :: x(no*no*no*nv)
      integer :: t

      do concurrent(t=1:no*no*no*nv)
         x(t) = x_b_item(t, no, nv, y, ooov)
      end do
   end subroutine kernel_pack_x_b

   subroutine kernel_add_sym(no, nv, factor, term, t2n)
      !! t2n(i,j,a,b) += factor [term(i,j,a,b) + term(j,i,b,a)]: the `tmp + tmp.transpose(1,0,3,2)` of every PySCF line
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: factor
      real(dp), intent(in) :: term(no, no, nv, nv)
      real(dp), intent(inout) :: t2n(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         t2n(t) = t2n(t) + sym_item(t, no, nv, factor, term)
      end do
   end subroutine kernel_add_sym

   subroutine kernel_add_ladder(no, nv, h, t2n)
      !! t2n(i,j,a,b) += h(i,j,a,b)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: h(no, no, nv, nv)
      real(dp), intent(inout) :: t2n(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         t2n(t) = t2n(t) + ladder_item(t, no, nv, h)
      end do
   end subroutine kernel_add_ladder

   subroutine kernel_pack_ring_x(no, nv, which, wvoov, wvovo, x)
      !! x(a,i,k,c): which = 1: 2 Wvoov(a,k,i,c) - Wvovo(a,k,c,i); 2: Wvoov(a,k,i,c); 3: Wvovo(a,k,c,i)
      integer, intent(in) :: no, nv, which
      real(dp), intent(in) :: wvoov(nv, no, no, nv), wvovo(nv, no, nv, no)
      real(dp), intent(out) :: x(nv*no*no*nv)
      integer :: t

      do concurrent(t=1:nv*no*no*nv)
         x(t) = ring_x_item(t, no, nv, which, wvoov, wvovo)
      end do
   end subroutine kernel_pack_ring_x

   subroutine kernel_pack_ring_y(no, nv, which, t2, y)
      !! y(k,c,j,b): which = 1: t2(k,j,c,b); 2: t2(k,j,b,c)
      integer, intent(in) :: no, nv, which
      real(dp), intent(in) :: t2(no, no, nv, nv)
      real(dp), intent(out) :: y(no*nv*no*nv)
      integer :: t

      do concurrent(t=1:no*nv*no*nv)
         y(t) = ring_y_item(t, no, nv, which, t2)
      end do
   end subroutine kernel_pack_ring_y

   subroutine kernel_add_ring(no, nv, which, factor, r, t2n)
      !! t2n(i,j,a,b) += factor [S(i,j,a,b) + S(j,i,b,a)], S scattered from r(ai,jb): which = 1 reads
      !! r((i-1)v+a, (b-1)o+j), which = 2 (the swapped term) r((i-1)v+b, (a-1)o+j)
      integer, intent(in) :: no, nv, which
      real(dp), intent(in) :: factor
      real(dp), intent(in) :: r(no*nv, no*nv)
      real(dp), intent(inout) :: t2n(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         t2n(t) = t2n(t) + ring_item(t, no, nv, which, factor, r)
      end do
   end subroutine kernel_add_ring

   subroutine kernel_denominator_t2(no, nv, eps_o, eps_v, t2n)
      !! t2n(i,j,a,b) /= e_i + e_j - e_a - e_b
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: eps_o(no), eps_v(nv)
      real(dp), intent(inout) :: t2n(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         t2n(t) = t2n(t)/denominator_item(t, no, nv, eps_o, eps_v)
      end do
   end subroutine kernel_denominator_t2

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

   pure function lamp_item(t, no, nv, lam) result(v)
      !! lamp(k,l,c,i) = lam(k,i,l,c), position t of (no, no, nv, no)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: lam(no, no, no, nv)
      real(dp) :: v
      integer :: k, l, c, i

      call split4(t, no, no, nv, k, l, c, i)
      v = lam(k, i, l, c)
   end function lamp_item

   pure function t1_small_item(t, no, nv, t1, t2, fvv, foo, fov, ovov, oovv) result(v)
      !! The four small T1 lines at (i,a), position t of (no, nv); each in the order of the module header
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv), fvv(nv, nv), foo(no, no), fov(no, nv)
      real(dp), intent(in) :: ovov(no, nv, no, nv), oovv(no, no, nv, nv)
      real(dp) :: v
      integer :: i, a, c, k

      i = mod(t - 1, no) + 1
      a = (t - 1)/no + 1
      v = 0.0_dp
      ! 'ac,ic->ia' and '-ki,ka->ia'
      do c = 1, nv
         v = v + fvv(a, c)*t1(i, c)
      end do
      do k = 1, no
         v = v - foo(k, i)*t1(k, a)
      end do
      ! '2 kc,kica->ia', '-kc,ikca->ia', 'kc,ic,ka->ia'
      do c = 1, nv
         do k = 1, no
            v = v + fov(k, c)*(2.0_dp*t2(k, i, c, a) - t2(i, k, c, a)) + fov(k, c)*t1(i, c)*t1(k, a)
         end do
      end do
      ! '2 kcai,kc->ia' with (kc|ai) = ovov(k,c,i,a); '-kiac,kc->ia'
      do c = 1, nv
         do k = 1, no
            v = v + (2.0_dp*ovov(k, c, i, a) - oovv(k, i, a, c))*t1(k, c)
         end do
      end do
   end function t1_small_item

   pure function ycomb_item(t, no, nv, tau) result(v)
      !! y(k,x',i,x) = 2 tau(i,k,x,x') - tau(i,k,x',x), position t of (no, nv, no, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp) :: v
      integer :: k, xp, i, x

      call split4(t, no, nv, no, k, xp, i, x)
      v = 2.0_dp*tau(i, k, x, xp) - tau(i, k, xp, x)
   end function ycomb_item

   pure function taur_item(t, no, nv, tau) result(v)
      !! taur(k,l,c,a) = tau(k,l,a,c), position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp) :: v
      integer :: k, l, c, a

      call split4(t, no, no, nv, k, l, c, a)
      v = tau(k, l, a, c)
   end function taur_item

   pure function base_item(t, no, nv, ovov) result(v)
      !! ovov(i,a,j,b) at t2n(i,j,a,b), position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp) :: v
      integer :: i, j, a, b

      call split4(t, no, no, nv, i, j, a, b)
      v = ovov(i, a, j, b)
   end function base_item

   pure function perm_r1_item(t, no, nv, r1) result(v)
      !! r1(i,a,j,b) at term(i,j,a,b), position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: r1(no, nv, no, nv)
      real(dp) :: v
      integer :: i, j, a, b

      call split4(t, no, no, nv, i, j, a, b)
      v = r1(i, a, j, b)
   end function perm_r1_item

   pure function x_b_item(t, no, nv, y, ooov) result(v)
      !! x(k,i,j,a) = y(k,j,i,a) + ooov(j,k,i,a), position t of (no, no, no, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: y(no, no, no, nv), ooov(no, no, no, nv)
      real(dp) :: v
      integer :: k, i, j, a

      call split4(t, no, no, no, k, i, j, a)
      v = y(k, j, i, a) + ooov(j, k, i, a)
   end function x_b_item

   pure function sym_item(t, no, nv, factor, term) result(v)
      !! factor [term(i,j,a,b) + term(j,i,b,a)] at position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: factor
      real(dp), intent(in) :: term(no, no, nv, nv)
      real(dp) :: v
      integer :: i, j, a, b

      call split4(t, no, no, nv, i, j, a, b)
      v = factor*(term(i, j, a, b) + term(j, i, b, a))
   end function sym_item

   pure function ladder_item(t, no, nv, h) result(v)
      !! h(i,j,a,b) at position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: h(no, no, nv, nv)
      real(dp) :: v
      integer :: i, j, a, b

      call split4(t, no, no, nv, i, j, a, b)
      v = h(i, j, a, b)
   end function ladder_item

   pure function ring_x_item(t, no, nv, which, wvoov, wvovo) result(v)
      !! x(a,i,k,c) at position t of (nv, no, no, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv, which
      real(dp), intent(in) :: wvoov(nv, no, no, nv), wvovo(nv, no, nv, no)
      real(dp) :: v
      integer :: a, i, k, c

      call split4(t, nv, no, no, a, i, k, c)
      if (which == 1) then
         v = 2.0_dp*wvoov(a, k, i, c) - wvovo(a, k, c, i)
      else if (which == 2) then
         v = wvoov(a, k, i, c)
      else
         v = wvovo(a, k, c, i)
      end if
   end function ring_x_item

   pure function ring_y_item(t, no, nv, which, t2) result(v)
      !! y(k,c,j,b) at position t of (no, nv, no, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv, which
      real(dp), intent(in) :: t2(no, no, nv, nv)
      real(dp) :: v
      integer :: k, c, j, b

      call split4(t, no, nv, no, k, c, j, b)
      if (which == 1) then
         v = t2(k, j, c, b)
      else
         v = t2(k, j, b, c)
      end if
   end function ring_y_item

   pure function ring_item(t, no, nv, which, factor, r) result(v)
      !! factor [S(i,j,a,b) + S(j,i,b,a)] at position t of (no, no, nv, nv), S scattered from r as in `kernel_add_ring`
      !$acc routine seq
      integer, intent(in) :: t, no, nv, which
      real(dp), intent(in) :: factor
      real(dp), intent(in) :: r(no*nv, no*nv)
      real(dp) :: v
      integer :: i, j, a, b

      call split4(t, no, no, nv, i, j, a, b)
      if (which == 1) then
         v = factor*(r((i - 1)*nv + a, (b - 1)*no + j) + r((j - 1)*nv + b, (a - 1)*no + i))
      else
         v = factor*(r((i - 1)*nv + b, (a - 1)*no + j) + r((j - 1)*nv + a, (b - 1)*no + i))
      end if
   end function ring_item

   pure function denominator_item(t, no, nv, eps_o, eps_v) result(v)
      !! e_i + e_j - e_a - e_b at position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: eps_o(no), eps_v(nv)
      real(dp) :: v
      integer :: i, j, a, b

      call split4(t, no, no, nv, i, j, a, b)
      v = eps_o(i) + eps_o(j) - eps_v(a) - eps_v(b)
   end function denominator_item

end module trc_rccsd_update
