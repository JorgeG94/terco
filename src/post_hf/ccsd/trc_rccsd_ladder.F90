!
! The particle-particle ladder of the t2 equation of a closed-shell RI-CCSD(T)
! calculation, in its two forms, and the t1 dressing that goes with it.
!
! Ported from the same author's mqc_czt_rcc.f90 (ladder_accumulate,
! ri_ladder_pair, ri_ladder_direct, ladder_t1_dressing, the fitted branch of
! particle_ladder). Chemists' notation, spatial orbitals, REAL orbitals; the
! integral blocks, tau and the amplitudes are those of trc_rccsd_ints and
! trc_rccsd_inter.
!
! WHAT IS COMPUTED. With the fitted (ac|bd) = sum_P b_vv(a,c; P) b_vv(b,d; P)
! (b_vv(ab,P) has a fastest), mqc's ladder and its dressing add to the new t2,
! before any symmetrisation and before the denominators, exactly
!
!   ladder    L(i,j,a,b)  = sum_cd (ac|bd) tau(i,j,c,d)           einsum 'ijcd,acbd->ijab'
!   dressing  D(i,j,a,b)  = - sum_cd [ sum_k ovvv(k,d,a,c) t1(k,b)
!                                    + sum_k ovvv(k,c,b,d) t1(k,a) ] tau(i,j,c,d)
!   total     L + D       = sum_cd Wvvvv(a,b,c,d) tau(i,j,c,d)    einsum 'abcd,ijcd->ijab'
!
! with Wvvvv(a,b,c,d) = (ac|bd) - sum_k ovvv(k,d,a,c) t1(k,b) - sum_k ovvv(k,c,b,d) t1(k,a)
! (pyscf/cc/rintermediates.py cc_Wvvvv: einsum('kdac,kb->abcd', ovvv, -t1) and
! einsum('kcbd,ka->abcd', ovvv, t1) added to eris.vvvv.transpose(0,2,1,3)).
! So L + D is the `t2new += lib.einsum('abcd,ijcd->ijab', Wvvvv, tau)` line of
! pyscf/cc/rccsd.py update_amps: no factor, no symmetrisation. PySCF's
! `mycc._add_vvvv(t1, t2, eris)` (ccsd.py _add_vvvv_full, einsum
! 'ijcd,acbd->ijab' on tau = t2 + t1 t1) is only L, the UNDRESSED ladder: the
! ovvv terms are applied elsewhere in pyscf's update_amps. The test therefore
! compares L with the `ladder` probe, and L + D, D with two further probes
! (`ladder_w`, `dressing`) that generate.py writes from cc_Wvvvv.
!
! SYMMETRY. L, D and L + D satisfy R(i,j,a,b) = R(j,i,b,a) when tau does
! (tau(j,i,d,c) = tau(i,j,c,d), true for t2 = MP2 and every converged t2): for L
! swap the roles of c and d and use (bd|ac) = (ac|bd); the two terms of D are
! mapped into each other by the same swap.
!
! LAYOUT. Results are in mqc's (i,j,a,b) layout, first index fastest, in
! `trc_rccsd_ladder_t%h`. `rccsd_ladder_pp` OVERWRITES it with L,
! `rccsd_ladder_dressing` ADDS D to it (so on a freshly initialised, zeroed
! type it gives D alone), `rccsd_ladder_build` is the two in that order.
!
! THE TWO FORMS (o = active occupied, v = virtual, P = auxiliary functions; MACs
! are multiply-adds, flops are twice that). Both never hold (ab|cd) whole:
!
!   ASSEMBLED  for a batch of nb compound columns cd = c + v (d-1)
!                W(a,b,cd) = (ac|bd) = sum_P b_vv(a,c;P) b_vv(b,d;P)
!              is one strided batched GEMM per run of columns that share d
!              (m = n = v, k = naux, the c-th A item is the rows (c-1)v+1..c v of
!              b_vv, B is the rows of d for every item): the raw (ab|cd) = one GEMM
!              of b_vv against b_vv would give (ab|cd), not the (ac|bd) the ladder
!              needs, and re-sorting it would cost a pass over v^2 nb. Then
!                L(ij; ab) += tau(ij; cd-batch) . W(ab; cd-batch)^T
!              is ONE GEMM per batch (m = o^2, n = v^2, k = nb; tau and the
!              result are the matrices their memory already is).
!                MACs  naux v^4 + o^2 v^4     scratch  v^2 cd_batch
!              The default cd_batch is the largest the memory budget allows.
!   DIRECT     for each auxiliary function P: G_P(a,c) = b_vv(a,c;P),
!                L(ij; ab) += sum_P G_P tau_ij G_P^T,     tau_ij(c,d) = tau(i,j,c,d)
!              restructured for the GPU as three steps per batch of np auxiliary
!              functions instead of o^2 np small products:
!                1. taut(c, d ij) = tau(i,j,c,d) (once per call);
!                   Y(a, d ij; P) = sum_c G_P(a,c) taut(c, d ij): ONE strided
!                   batched GEMM over P with B shared (m = v, n = v o^2, k = v);
!                2. a permutation pack Y2(ij, a, d, P) = Y(a, d ij; P);
!                3. L(ij a; b) += sum_{d,P} Y2(ij a; d P) b_vv(b, d; P)
!                   -- ONE GEMM (m = o^2 v, n = v, k = v np), because b_vv's
!                   columns P..P+np-1 ARE the matrix (v; d P) with ld v.
!                MACs  2 naux o^2 v^3          scratch  taut + 2 o^2 v^2 np
!              (the pack moves o^2 v^2 np doubles twice per batch).
!   The choice (`rccsd_ladder_prefers_direct`) is mqc's: on counted MACs,
!   (naux + o^2) v^4 against 2 naux o^2 v^3, so direct wins when v > 2 naux o^2 /
!   (naux + o^2), about 2 o^2 for naux >> o^2. `options%form` forces either one.
!
!   The dressing is the same for both forms and O(o^3 v^3):
!     Z1(k,a; ij) = sum_cd ovvv(k,d,a,c) tau(i,j,c,d)     D1(ij,a,b) = - sum_k Z1(k,a; ij) t1(k,b)
!     Z2(k,b; ij) = sum_cd ovvv(k,c,b,d) tau(i,j,c,d)     D2(ij,a,b) = - sum_k Z2(k,b; ij) t1(k,a)
!   Per virtual a (b) the pack pk(k; cd) = ovvv(k,d,a,c) (ovvv(k,c,b,d)) is
!   o v^2, then Z = pk . taut (m = o, n = o^2, k = v^2) and cmat(ij; b) =
!   Z^T t1 (m = o^2, n = v, k = o); `virtual_batch` virtuals go in one strided
!   batched GEMM each (B shared), and a scatter kernel subtracts cmat. Packing
!   ovvv for all virtuals at once would be an o v^3 array, which is why it is
!   chunked.
!                MACs  2 o^3 v^3 + 2 o^3 v^2        (both halves)
!                scratch  taut + (o v^2 + o^3 + o^2 v) virtual_batch
!   Flops are twice the MACs: `rccsd_ladder_flops` for the undressed form that ran (its counter, `rccsd_ladder_macs`,
!   matches the GEMMs both forms call for any batch plan; `ls%flops_pp` tallies the calls to show it) and
!   `rccsd_ladder_dressing_flops` for the dressing, whose MACs `rccsd_ladder_macs` does not include.
!
! NUMBERS for no = 25, nv = 130, naux = 630 (default budget 512 MiB):
!     assembled   717 GFLOP (build 360 + contraction 357),  wbatch 3970 columns = 537 MB
!     direct     3460 GFLOP (4.8 x the assembled; the selector picks assembled),  np = 3: 2 * 254 MB
!     dressing    138 GFLOP (19 % of the assembled), virtual_batch 129 = 537 MB
!     taut        84.5 MB, the result 84.5 MB.
!
! DETERMINISM. The batches are visited in ascending order and each GEMM
! accumulates into the result with beta = 1, the scatter kernels write each
! element once: on the device the order of every sum is fixed. On the host a
! threaded BLAS regroups the sums inside a GEMM, which changes round-off only
! (assembled against direct: 1.7e-17 on one thread, 2.1e-17 on four). Changing a
! batch size regroups the sums over cd (assembled) or over (d,P) (direct) inside
! the GEMM: round-off only.
!
! RESIDENCY. Every array of the types here exists on the host AND on the device
! (`enter data` in `init`, `exit data` in `release`, on the same `x%name`
! expression). The device copy is the one the kernels write; `pull` brings the
! result back, `push` sends the host copy up. Each component is mapped on its own
! and the parent type is never mapped, so a kernel gets an array only as an
! explicit-shape dummy of a `kernel_*` routine and never names `x%name` in a
! `do concurrent` body (an illegal address under mem:separate that the host build
! hides). The host copies of the scratch are never read.
!
! STALE SCRATCH. The scratch holds no integral-derived data, but its sizes and
! batch plan are those of one `ints%build`: it records `ints%stamp` and every
! routine refuses (message in `ls%status`, nothing written) after the integrals
! were rebuilt. Re-`init` it.
!
! Failures are reported in `ls%status` (a routine whose arguments differ in
! (no, nv, naux), are not built, or hold arrays of other sizes sets it and
! returns without writing anything); none stops. A block of more than 2^31-1
! elements is refused by `init`, not wrapped. Single rank only.
module trc_rccsd_ladder
   use pic_types, only: dp, int32, int64
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t
   use trc_rccsd_steps, only: RCCSD_ST_LADDER, RCCSD_ST_DRESSING, rccsd_steps_gemm_flops, rccsd_steps_open, rccsd_steps_close
   use pic_timer, only: timer_type
   use trc_step_table, only: step_table_t
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   character(len=*), parameter :: STALE_MSG = "trc_rccsd_ladder: the integrals were rebuilt after the scratch was "// &
                                              "initialised; re-init the scratch after ints%build"

   integer, parameter, public :: RCCSD_LADDER_AUTO = 0       !! `options%form`: pick by counted operations
   integer, parameter, public :: RCCSD_LADDER_ASSEMBLED = 1  !! `options%form`: build (ac|bd) in column batches
   integer, parameter, public :: RCCSD_LADDER_DIRECT = 2     !! `options%form`: contract G_P tau G_P^T, never form (ac|bd)
   integer(int64), parameter :: DEFAULT_BUDGET = 536870912_int64  !! bytes per working set when `budget_bytes` is 0

   public :: trc_rccsd_ladder_t, trc_rccsd_lscratch_t, trc_rccsd_ladder_options_t
   public :: rccsd_ladder_pp, rccsd_ladder_dressing, rccsd_ladder_build
   public :: rccsd_ladder_prefers_direct, rccsd_ladder_macs, rccsd_ladder_flops
   public :: rccsd_ladder_dressing_macs, rccsd_ladder_dressing_flops

   type :: trc_rccsd_ladder_options_t
      !! What `init` of the scratch needs to know besides the integrals. Plain values only. A zero is "choose".
      integer :: form = RCCSD_LADDER_AUTO  !! RCCSD_LADDER_AUTO, _ASSEMBLED or _DIRECT
      integer :: cd_batch = 0              !! compound (cd) columns per pass of the assembled form
      integer :: p_batch = 0               !! auxiliary functions per pass of the direct form
      integer :: virtual_batch = 0         !! virtuals per pass of the dressing
      integer(int64) :: budget_bytes = 0   !! bytes per working set (the form's, the dressing's); 0 means 512 MiB
   end type trc_rccsd_ladder_options_t

   type :: trc_rccsd_ladder_t
      !! The ladder contribution on host and device, valid after the matching routine
      type(error_t) :: status  !! set by `init`: a refusal leaves its message here and the type empty
      integer :: no = 0
      integer :: nv = 0
      real(dp), allocatable :: h(:, :, :, :)  !! h(i,j,a,b), zeroed by `init`
   contains
      procedure :: init => ladder_init
      procedure :: release => ladder_release
      procedure :: push => ladder_push
      procedure :: pull => ladder_pull
   end type trc_rccsd_ladder_t

   type :: trc_rccsd_lscratch_t
      !! The work arrays of the routines below, sizes in the module header
      type(error_t) :: status  !! set by `init` and by every routine that refuses its arguments
      integer :: no = 0
      integer :: nv = 0
      integer :: naux = 0
      integer :: stamp = 0          !! `ints%stamp` of the integrals `init` planned for
      integer :: form = 0           !! the form in use: RCCSD_LADDER_ASSEMBLED or _DIRECT (never AUTO after `init`)
      integer :: cd_batch = 0       !! columns per pass, assembled form
      integer :: p_batch = 0        !! auxiliary functions per pass, direct form
      integer :: virtual_batch = 0  !! virtuals per pass, dressing
      integer(int64) :: bytes = 0   !! bytes of the work arrays, host and device each
      integer(int64) :: flops_pp = 0        !! 2 m n k summed over the GEMMs the last undressed ladder called, batches counted
      integer(int64) :: flops_dressing = 0  !! the same for the last dressing: the tallies the step table reports
      real(dp), allocatable :: taut(:)    !! o^2 v^2, taut(c,d,ij) = tau(i,j,c,d), per call
      real(dp), allocatable :: ovpack(:)  !! o v^2 virtual_batch, ovvv packed per virtual (dressing)
      real(dp), allocatable :: zmat(:)    !! o^3 virtual_batch, Z(k; ij) per virtual
      real(dp), allocatable :: cmat(:)    !! o^2 v virtual_batch, Z^T t1 per virtual
      real(dp), allocatable :: wbatch(:)  !! v^2 cd_batch, (ac|bd) per column (assembled only)
      real(dp), allocatable :: gtau(:)    !! o^2 v^2 p_batch, Y(a, d ij; P) (direct only)
      real(dp), allocatable :: gpack(:)   !! o^2 v^2 p_batch, Y2(ij, a, d, P) (direct only)
   contains
      procedure :: init => lscratch_init
      procedure :: release => lscratch_release
   end type trc_rccsd_lscratch_t

contains

   ! ---------------------------------------------------------------------
   ! the types
   ! ---------------------------------------------------------------------

   subroutine ladder_init(this, no, nv)
      !! A zeroed result on host and device. A block above 2^31-1 elements is refused in `status`.
      class(trc_rccsd_ladder_t), intent(inout) :: this
      integer, intent(in) :: no, nv

      call this%release()
      call this%status%clear()
      if (no < 0 .or. nv < 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: negative orbital count")
         return
      end if
      if (elements(no, no, nv, nv) > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: a block has more than 2^31-1 elements")
         return
      end if
      this%no = no
      this%nv = nv
      allocate (this%h(no, no, nv, nv))
      this%h = 0.0_dp
      !$acc enter data copyin(this%h)
   end subroutine ladder_init

   subroutine ladder_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_ladder_t), intent(inout) :: this

      if (allocated(this%h)) then
         !$acc exit data delete(this%h)
         deallocate (this%h)
      end if
      this%no = 0
      this%nv = 0
   end subroutine ladder_release

   subroutine ladder_push(this)
      !! Host copy to device
      class(trc_rccsd_ladder_t), intent(inout) :: this

      if (.not. allocated(this%h)) return
      !$acc update device(this%h)
   end subroutine ladder_push

   subroutine ladder_pull(this)
      !! Device copy to host
      class(trc_rccsd_ladder_t), intent(inout) :: this

      if (.not. allocated(this%h)) return
      !$acc update self(this%h)
   end subroutine ladder_pull

   pure function plan_batch(requested, budget, item_elements, upper) result(n)
      !! Items per pass: `requested` if positive, else as many as `budget` bytes hold; never above `upper`
      !! or above what a default-integer index addresses; 0 when not even one item fits
      integer, intent(in) :: requested, upper
      integer(int64), intent(in) :: budget, item_elements
      integer :: n

      integer(int64) :: cap

      cap = int(huge(0_int32), int64)/max(item_elements, 1_int64)
      if (requested > 0) then
         cap = min(cap, int(requested, int64))
      else
         cap = min(cap, max(1_int64, budget/(8_int64*max(item_elements, 1_int64))))
      end if
      cap = min(cap, int(upper, int64))
      n = int(cap)
      if (cap < 1) n = 0
   end function plan_batch

   subroutine lscratch_init(this, ints, options)
      !! Plans the batches and allocates the work arrays for these integrals. Refusals go to `status`.
      class(trc_rccsd_lscratch_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_ladder_options_t), intent(in), optional :: options

      type(trc_rccsd_ladder_options_t) :: opts
      integer :: no, nv, naux
      integer(int64) :: budget, no2, nv2

      call this%release()
      call this%status%clear()
      if (.not. allocated(ints%b_vv) .or. .not. allocated(ints%ovvv)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the integral blocks are not built")
         return
      end if
      opts = trc_rccsd_ladder_options_t()
      if (present(options)) opts = options
      if (opts%form < RCCSD_LADDER_AUTO .or. opts%form > RCCSD_LADDER_DIRECT .or. opts%cd_batch < 0 .or. &
          opts%p_batch < 0 .or. opts%virtual_batch < 0 .or. opts%budget_bytes < 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: an unknown form or a negative batch size or budget")
         return
      end if
      no = ints%no
      nv = ints%nv
      naux = ints%naux
      no2 = int(no, int64)*no
      nv2 = int(nv, int64)*nv
      if (no < 1 .or. nv < 1) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: no active occupied or no virtual orbitals")
         return
      end if
      ! (ia|bc) too: the dressing kernels index it with default integers, and
      ! o v^3 passes 2^31 before o^2 v^2 does once v > o.
      if (no2*nv2 > int(huge(0_int32), int64) .or. int(no, int64)*nv2*nv > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: a block has more than 2^31-1 elements")
         return
      end if
      budget = merge(opts%budget_bytes, DEFAULT_BUDGET, opts%budget_bytes > 0)
      this%no = no
      this%nv = nv
      this%naux = naux
      this%form = opts%form
      if (this%form == RCCSD_LADDER_AUTO) then
         this%form = merge(RCCSD_LADDER_DIRECT, RCCSD_LADDER_ASSEMBLED, rccsd_ladder_prefers_direct(no, nv, naux))
      end if
      this%virtual_batch = plan_batch(opts%virtual_batch, budget, int(no, int64)*nv2 + no2*no + no2*nv, nv)
      if (this%form == RCCSD_LADDER_ASSEMBLED) then
         this%cd_batch = plan_batch(opts%cd_batch, budget, nv2, int(nv2))
      else
         this%p_batch = plan_batch(opts%p_batch, budget, 2_int64*no2*nv2, naux)
      end if
      if (this%virtual_batch < 1 .or. (this%form == RCCSD_LADDER_ASSEMBLED .and. this%cd_batch < 1) .or. &
          (this%form == RCCSD_LADDER_DIRECT .and. this%p_batch < 1)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: not even one batch item fits in a default integer index")
         call this%release()
         return
      end if
      this%stamp = ints%stamp
      call allocate_work(this)
   end subroutine lscratch_init

   subroutine allocate_work(this)
      !! The arrays of the plan in `this`, mapped to the device
      type(trc_rccsd_lscratch_t), intent(inout) :: this

      integer(int64) :: no, nv, no2, nv2

      no = this%no
      nv = this%nv
      no2 = no*no
      nv2 = nv*nv
      allocate (this%taut(no2*nv2))
      allocate (this%ovpack(no*nv2*this%virtual_batch), this%zmat(no2*no*this%virtual_batch), &
                this%cmat(no2*nv*this%virtual_batch))
      !$acc enter data create(this%taut, this%ovpack, this%zmat, this%cmat)
      this%bytes = 8_int64*(no2*nv2 + (no*nv2 + no2*no + no2*nv)*this%virtual_batch)
      if (this%form == RCCSD_LADDER_ASSEMBLED) then
         allocate (this%wbatch(nv2*this%cd_batch))
         !$acc enter data create(this%wbatch)
         this%bytes = this%bytes + 8_int64*nv2*this%cd_batch
      else
         allocate (this%gtau(no2*nv2*this%p_batch), this%gpack(no2*nv2*this%p_batch))
         !$acc enter data create(this%gtau, this%gpack)
         this%bytes = this%bytes + 16_int64*no2*nv2*this%p_batch
      end if
   end subroutine allocate_work

   subroutine lscratch_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_lscratch_t), intent(inout) :: this

      if (allocated(this%taut)) then
         !$acc exit data delete(this%taut)
         deallocate (this%taut)
      end if
      if (allocated(this%ovpack)) then
         !$acc exit data delete(this%ovpack)
         deallocate (this%ovpack)
      end if
      if (allocated(this%zmat)) then
         !$acc exit data delete(this%zmat)
         deallocate (this%zmat)
      end if
      if (allocated(this%cmat)) then
         !$acc exit data delete(this%cmat)
         deallocate (this%cmat)
      end if
      if (allocated(this%wbatch)) then
         !$acc exit data delete(this%wbatch)
         deallocate (this%wbatch)
      end if
      if (allocated(this%gtau)) then
         !$acc exit data delete(this%gtau)
         deallocate (this%gtau)
      end if
      if (allocated(this%gpack)) then
         !$acc exit data delete(this%gpack)
         deallocate (this%gpack)
      end if
      this%no = 0
      this%nv = 0
      this%naux = 0
      this%stamp = 0
      this%form = 0
      this%cd_batch = 0
      this%p_batch = 0
      this%virtual_batch = 0
      this%bytes = 0
      this%flops_pp = 0_int64
      this%flops_dressing = 0_int64
   end subroutine lscratch_release

   ! ---------------------------------------------------------------------
   ! cost model and size checks
   ! ---------------------------------------------------------------------

   pure function rccsd_ladder_macs(form, no, nv, naux) result(macs)
      !! Multiply-adds of one undressed ladder: assembled (naux + o^2) v^4, direct 2 naux o^2 v^3; 0 for another form
      integer, intent(in) :: form, no, nv, naux
      real(dp) :: macs

      macs = 0.0_dp
      if (form == RCCSD_LADDER_ASSEMBLED) macs = (real(naux, dp) + real(no, dp)**2)*real(nv, dp)**4
      if (form == RCCSD_LADDER_DIRECT) macs = 2.0_dp*real(naux, dp)*real(nv, dp)**3*real(no, dp)**2
   end function rccsd_ladder_macs

   pure function rccsd_ladder_flops(form, no, nv, naux) result(flops)
      !! Flops of one undressed ladder, twice `rccsd_ladder_macs`: every GEMM of both forms is counted in it (the
      !! build of (ac|bd) and its contraction, or the two products around the pack), and `ls%flops_pp` tallies the
      !! calls actually made, which agree for any batch plan. 0 for another form.
      integer, intent(in) :: form, no, nv, naux
      integer(int64) :: flops

      flops = nint(2.0_dp*rccsd_ladder_macs(form, no, nv, naux), int64)
   end function rccsd_ladder_flops

   pure function rccsd_ladder_dressing_macs(no, nv) result(macs)
      !! Multiply-adds of the dressing, both halves: Z = pk . taut is o * o^2 * v^2 per virtual, cmat = Z^T t1 is
      !! o^2 * v * o per virtual, for v virtuals each: 2 o^3 v^3 + 2 o^3 v^2 (the macs the header quotes)
      integer, intent(in) :: no, nv
      real(dp) :: macs

      macs = 2.0_dp*real(no, dp)**3*real(nv, dp)**3 + 2.0_dp*real(no, dp)**3*real(nv, dp)**2
   end function rccsd_ladder_dressing_macs

   pure function rccsd_ladder_dressing_flops(no, nv) result(flops)
      !! Flops of one dressing, twice `rccsd_ladder_dressing_macs`; `ls%flops_dressing` tallies the calls actually made
      integer, intent(in) :: no, nv
      integer(int64) :: flops

      flops = nint(2.0_dp*rccsd_ladder_dressing_macs(no, nv), int64)
   end function rccsd_ladder_dressing_flops

   pure function rccsd_ladder_prefers_direct(no, nv, naux) result(direct)
      !! Whether the direct form needs fewer multiply-adds than the assembled one (mqc's rule: counts only,
      !! because both forms are GEMMs of comparable shape, so the counts predict the order but not the constant)
      integer, intent(in) :: no, nv, naux
      logical :: direct

      direct = rccsd_ladder_macs(RCCSD_LADDER_DIRECT, no, nv, naux) < rccsd_ladder_macs(RCCSD_LADDER_ASSEMBLED, no, nv, naux)
   end function rccsd_ladder_prefers_direct

   pure function elements(n1, n2, n3, n4) result(n)
      !! n1 n2 n3 n4 as a 64-bit count
      integer, intent(in) :: n1, n2, n3, n4
      integer(int64) :: n

      n = int(n1, int64)*n2*n3*n4
   end function elements

   function ready(ints, lad, ls) result(ok)
      !! True when the integrals, the result and the scratch are built, agree on (no, nv, naux), hold arrays of
      !! those sizes and the scratch is not stale; else `ls%status` says which not
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_ladder_t), intent(in) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      logical :: ok

      integer :: no, nv
      integer(int64) :: no2, nv2

      ok = .false.
      if (.not. allocated(ints%b_vv) .or. .not. allocated(ints%ovvv)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the integral blocks are not built")
         return
      end if
      if (.not. allocated(lad%h)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the result is not initialised")
         return
      end if
      if (.not. allocated(ls%taut) .or. .not. allocated(ls%ovpack) .or. .not. allocated(ls%zmat) .or. &
          .not. allocated(ls%cmat)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the scratch is not initialised")
         return
      end if
      no = ints%no
      nv = ints%nv
      if (lad%no /= no .or. lad%nv /= nv .or. ls%no /= no .or. ls%nv /= nv .or. ls%naux /= ints%naux) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: (no, nv, naux) differ between integrals, result and scratch")
         return
      end if
      no2 = int(no, int64)*no
      nv2 = int(nv, int64)*nv
      if (size(ints%b_vv, kind=int64) /= nv2*ints%naux .or. size(ints%ovvv, kind=int64) /= elements(no, nv, nv, nv)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: an integral block does not have the size (no, nv) says")
         return
      end if
      if (size(lad%h, kind=int64) /= no2*nv2) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the result does not have the size (no, nv) says")
         return
      end if
      if (.not. scratch_sized(ls, no2, nv2)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: a scratch array does not have the size its plan says")
         return
      end if
      if (ls%stamp /= ints%stamp) then
         call ls%status%set(ERROR_VALIDATION, STALE_MSG)
         return
      end if
      ok = .true.
   end function ready

   function scratch_sized(ls, no2, nv2) result(ok)
      !! True when every array of the plan is allocated with the size the plan gives
      type(trc_rccsd_lscratch_t), intent(in) :: ls
      integer(int64), intent(in) :: no2, nv2
      logical :: ok

      integer(int64) :: no, vb

      no = ls%no
      vb = ls%virtual_batch
      ok = size(ls%taut, kind=int64) == no2*nv2 .and. size(ls%ovpack, kind=int64) == no*nv2*vb .and. &
           size(ls%zmat, kind=int64) == no2*no*vb .and. size(ls%cmat, kind=int64) == no2*ls%nv*vb
      if (.not. ok) return
      if (ls%form == RCCSD_LADDER_ASSEMBLED) then
         ok = allocated(ls%wbatch)
         if (ok) ok = size(ls%wbatch, kind=int64) == nv2*ls%cd_batch
      else if (ls%form == RCCSD_LADDER_DIRECT) then
         ok = allocated(ls%gtau) .and. allocated(ls%gpack)
         if (ok) ok = size(ls%gtau, kind=int64) == no2*nv2*ls%p_batch .and. size(ls%gpack, kind=int64) == no2*nv2*ls%p_batch
      else
         ok = .false.
      end if
   end function scratch_sized

   function ready_tau(ints, inter, ls) result(ok)
      !! True when `inter` holds a tau of the integrals' (no, nv); else `ls%status` says so
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      logical :: ok

      ok = .false.
      if (.not. allocated(inter%tau)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the intermediates (tau) are not initialised")
         return
      end if
      if (inter%no /= ints%no .or. inter%nv /= ints%nv .or. &
          size(inter%tau, kind=int64) /= elements(ints%no, ints%no, ints%nv, ints%nv)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: tau does not have the size (no, nv) says")
         return
      end if
      ok = .true.
   end function ready_tau

   function ready_amps(ints, amps, ls) result(ok)
      !! True when `amps` holds a t1 of the integrals' (no, nv); else `ls%status` says so
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      logical :: ok

      ok = .false.
      if (.not. allocated(amps%t1)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the amplitudes are not initialised")
         return
      end if
      if (amps%no /= ints%no .or. amps%nv /= ints%nv .or. &
          size(amps%t1, kind=int64) /= elements(ints%no, ints%nv, 1, 1)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: t1 does not have the size (no, nv) says")
         return
      end if
      ok = .true.
   end function ready_amps

   ! ---------------------------------------------------------------------
   ! the ladder
   ! ---------------------------------------------------------------------

   subroutine rccsd_ladder_pp(ints, inter, lad, ls, la)
      !! L(i,j,a,b) = sum_cd (ac|bd) tau(i,j,c,d), overwriting `lad%h` (pyscf _add_vvvv, einsum 'ijcd,acbd->ijab')
      !!
      !! The form is `ls%form`. Needs tau at the current amplitudes (`rccsd_inter_build_tau`).
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      call ls%status%clear()
      if (.not. ready(ints, lad, ls)) return
      if (.not. ready_tau(ints, inter, ls)) return
      call ladder_pp_run(ints, inter, lad, ls, la)
   end subroutine rccsd_ladder_pp

   subroutine ladder_pp_run(ints, inter, lad, ls, la)
      !! `rccsd_ladder_pp` after its checks
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      ls%flops_pp = 0_int64
      call kernel_zero(ints%no*ints%no*ints%nv*ints%nv, lad%h)
      if (ls%form == RCCSD_LADDER_ASSEMBLED) then
         call ladder_assembled(ints, inter, lad, ls, la)
      else
         call ladder_direct(ints, inter, lad, ls, la)
      end if
   end subroutine ladder_pp_run

   subroutine ladder_assembled(ints, inter, lad, ls, la)
      !! The assembled form (module header): per batch of cd columns, W by batched GEMMs, then one GEMM with tau
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: shape
      integer :: no, nv, naux, nv2, cd0, cd1, cd, col, nb, c, d, nrun

      no = ints%no
      nv = ints%nv
      naux = ints%naux
      nv2 = nv*nv
      do cd0 = 1, nv2, ls%cd_batch
         cd1 = min(cd0 + ls%cd_batch - 1, nv2)
         nb = cd1 - cd0 + 1
         ! W(a,b,col) = (ac|bd) = sum_P b_vv(a,c;P) b_vv(b,d;P), one batched call per run of columns with the same d:
         ! item c reads the rows (c-1)v+1..c v of b_vv, B is the rows of d for all of them.
         cd = cd0
         col = 0
         do while (cd <= cd1)
            c = mod(cd - 1, nv) + 1
            d = (cd - 1)/nv + 1
            nrun = min(nv - c + 1, cd1 - cd + 1)
            shape = trc_gemm_batch(nv, nv, naux, nrun, lda=nv2, ldb=nv2, ldc=nv, &
                                   sa=int(nv, int64), sb=0_int64, sc=int(nv2, int64), transb="T")
            call la%gemm_batched(shape, ints%b_vv((c - 1)*nv + 1, 1), ints%b_vv((d - 1)*nv + 1, 1), &
                                 ls%wbatch(col*nv2 + 1))
            ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(shape%m, shape%n, shape%k, shape%count)
            col = col + nrun
            cd = cd + nrun
         end do
         ! L(ij; ab) += tau(ij; cd-batch) W(ab; cd-batch)^T. tau is (o^2, v^2) in memory and the batch is a
         ! contiguous column range of it, starting at the (c,d) of cd0.
         c = mod(cd0 - 1, nv) + 1
         d = (cd0 - 1)/nv + 1
         call la%gemm("N", "T", no*no, nv2, nb, 1.0_dp, inter%tau(1, 1, c, d), no*no, ls%wbatch, nv2, 1.0_dp, lad%h, no*no)
         ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(no*no, nv2, nb)
      end do
   end subroutine ladder_assembled

   subroutine ladder_direct(ints, inter, lad, ls, la)
      !! The direct form (module header): per batch of auxiliary functions, GEMM against tau, a pack, GEMM against b_vv
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: shape
      integer :: no, nv, naux, p0, np

      no = ints%no
      nv = ints%nv
      naux = ints%naux
      call kernel_pack_taut(no, nv, inter%tau, ls%taut)
      do p0 = 1, naux, ls%p_batch
         np = min(ls%p_batch, naux - p0 + 1)
         ! Y(a, d ij; P) = sum_c G_P(a,c) taut(c, d ij): G_P is the (v, v) column P of b_vv, taut is shared
         shape = trc_gemm_batch(nv, nv*no*no, nv, np, lda=nv, ldb=nv, ldc=nv, &
                                sa=int(nv, int64)*nv, sb=0_int64, sc=int(nv, int64)*nv*no*no)
         call la%gemm_batched(shape, ints%b_vv(1, p0), ls%taut, ls%gtau)
         ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(shape%m, shape%n, shape%k, shape%count)
         call kernel_pack_y2(no, nv, np, ls%gtau, ls%gpack)
         ! L(ij a; b) += sum_{d,P} Y2(ij a; d P) b_vv(b, d; P): columns p0.. of b_vv are the matrix (v; d P), ld v
         call la%gemm("N", "T", no*no*nv, nv, nv*np, 1.0_dp, ls%gpack, no*no*nv, ints%b_vv(1, p0), nv, 1.0_dp, &
                      lad%h, no*no*nv)
         ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(no*no*nv, nv, nv*np)
      end do
   end subroutine ladder_direct

   subroutine rccsd_ladder_dressing(ints, amps, inter, lad, ls, la)
      !! Adds D(i,j,a,b) = - sum_cd [sum_k ovvv(k,d,a,c) t1(k,b) + sum_k ovvv(k,c,b,d) t1(k,a)] tau(i,j,c,d) to `lad%h`
      !!
      !! pyscf cc_Wvvvv's einsum('kdac,kb->abcd', ovvv, -t1) and einsum('kcbd,ka->abcd', ovvv, t1), contracted with
      !! tau. Needs tau at the current amplitudes and `amps` on the device (the caller pushes after a host edit). It
      !! ADDS: call it on a zeroed result for D alone, or after `rccsd_ladder_pp` for the total.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      call ls%status%clear()
      if (.not. ready(ints, lad, ls)) return
      if (.not. ready_tau(ints, inter, ls)) return
      if (.not. ready_amps(ints, amps, ls)) return
      call ladder_dressing_run(ints, amps, inter, lad, ls, la)
   end subroutine rccsd_ladder_dressing

   subroutine ladder_dressing_run(ints, amps, inter, lad, ls, la)
      !! `rccsd_ladder_dressing` after its checks
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      integer :: half

      ls%flops_dressing = 0_int64
      call kernel_pack_taut(ints%no, ints%nv, inter%tau, ls%taut)
      do half = 1, 2
         call dressing_half(half, ints, amps, lad, ls, la)
      end do
   end subroutine ladder_dressing_run

   subroutine dressing_half(half, ints, amps, lad, ls, la)
      !! half = 1: D1 (free virtual a is the first of the pair), 2: D2 (b is the second); virtual_batch at a time
      integer, intent(in) :: half
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: shape_z, shape_c
      integer :: no, nv, no2, nv2, v0, nvb

      no = ints%no
      nv = ints%nv
      no2 = no*no
      nv2 = nv*nv
      do v0 = 1, nv, ls%virtual_batch
         nvb = min(ls%virtual_batch, nv - v0 + 1)
         if (half == 1) then
            call kernel_pack_z1(no, nv, v0, nvb, ints%ovvv, ls%ovpack)
         else
            call kernel_pack_z2(no, nv, v0, nvb, ints%ovvv, ls%ovpack)
         end if
         ! Z(k; ij) = sum_cd pk(k; cd) taut(cd; ij), one item per virtual, taut shared
         shape_z = trc_gemm_batch(no, no2, nv2, nvb, lda=no, ldb=nv2, ldc=no, &
                                  sa=int(no, int64)*nv2, sb=0_int64, sc=int(no, int64)*no2)
         call la%gemm_batched(shape_z, ls%ovpack, ls%taut, ls%zmat)
         ls%flops_dressing = ls%flops_dressing + rccsd_steps_gemm_flops(shape_z%m, shape_z%n, shape_z%k, shape_z%count)
         ! cmat(ij; x) = sum_k Z(k; ij) t1(k,x), t1 shared
         shape_c = trc_gemm_batch(no2, nv, no, nvb, transa="T", lda=no, ldb=no, ldc=no2, &
                                  sa=int(no, int64)*no2, sb=0_int64, sc=int(no2, int64)*nv)
         call la%gemm_batched(shape_c, ls%zmat, amps%t1, ls%cmat)
         ls%flops_dressing = ls%flops_dressing + rccsd_steps_gemm_flops(shape_c%m, shape_c%n, shape_c%k, shape_c%count)
         if (half == 1) then
            call kernel_sub_z1(no, nv, v0, nvb, ls%cmat, lad%h)
         else
            call kernel_sub_z2(no, nv, v0, nvb, ls%cmat, lad%h)
         end if
      end do
   end subroutine dressing_half

   subroutine rccsd_ladder_build(ints, amps, inter, lad, ls, la, stab)
      !! L + D = sum_cd Wvvvv(a,b,c,d) tau(i,j,c,d) in `lad%h`, overwriting it (pyscf rccsd.update_amps's
      !! einsum('abcd,ijcd->ijab', Wvvvv, tau)). Every check is made before the first write, so a refusal changes nothing.
      !! With `stab` the undressed ladder and the dressing are two steps of that table (RCCSD_ST_LADDER and
      !! RCCSD_ST_DRESSING), each with its counted flops.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      type(step_table_t), intent(inout), optional :: stab

      type(timer_type) :: tc

      call ls%status%clear()
      if (.not. ready(ints, lad, ls)) return
      if (.not. ready_tau(ints, inter, ls)) return
      if (.not. ready_amps(ints, amps, ls)) return
      call rccsd_steps_open(stab, la, RCCSD_ST_LADDER, tc)
      call ladder_pp_run(ints, inter, lad, ls, la)
      ! The tallies of the GEMMs actually called, not the formulas: the table then stays right when the
      ! GEMMs change (the check compares the tallies with the formulas).
      call rccsd_steps_close(stab, la, RCCSD_ST_LADDER, tc, ls%flops_pp)
      call rccsd_steps_open(stab, la, RCCSD_ST_DRESSING, tc)
      call ladder_dressing_run(ints, amps, inter, lad, ls, la)
      call rccsd_steps_close(stab, la, RCCSD_ST_DRESSING, tc, ls%flops_dressing)
   end subroutine rccsd_ladder_build

   ! ---------------------------------------------------------------------
   ! the kernels. Each takes its arrays as explicit-shape dummies and runs one
   ! `do concurrent` over ONE flattened index, whose body calls an item routine
   ! of this module. No `x%name` appears in any of them.
   ! ---------------------------------------------------------------------

   subroutine kernel_zero(n, h)
      !! h = 0
      integer, intent(in) :: n
      real(dp), intent(out) :: h(n)
      integer :: t

      do concurrent(t=1:n)
         h(t) = 0.0_dp
      end do
   end subroutine kernel_zero

   subroutine kernel_pack_taut(no, nv, tau, taut)
      !! taut(c,d,ij) = tau(i,j,c,d), ij = i + o (j-1)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(out) :: taut(nv*nv*no*no)
      integer :: t

      do concurrent(t=1:nv*nv*no*no)
         taut(t) = taut_item(t, no, nv, tau)
      end do
   end subroutine kernel_pack_taut

   subroutine kernel_pack_y2(no, nv, np, y, y2)
      !! y2(ij,a,d,P) = y(a,d,ij,P)
      integer, intent(in) :: no, nv, np
      real(dp), intent(in) :: y(nv*nv*no*no*np)
      real(dp), intent(out) :: y2(no*no*nv*nv*np)
      integer :: t

      do concurrent(t=1:no*no*nv*nv*np)
         y2(t) = y2_item(t, no, nv, np, y)
      end do
   end subroutine kernel_pack_y2

   subroutine kernel_pack_z1(no, nv, v0, nvb, ovvv, ovpack)
      !! ovpack(k, c+v(d-1), al) = ovvv(k, d, v0+al-1, c)
      integer, intent(in) :: no, nv, v0, nvb
      real(dp), intent(in) :: ovvv(no, nv, nv, nv)
      real(dp), intent(out) :: ovpack(no*nv*nv*nvb)
      integer :: t

      do concurrent(t=1:no*nv*nv*nvb)
         ovpack(t) = z1_item(t, no, nv, v0, ovvv)
      end do
   end subroutine kernel_pack_z1

   subroutine kernel_pack_z2(no, nv, v0, nvb, ovvv, ovpack)
      !! ovpack(k, c+v(d-1), bl) = ovvv(k, c, v0+bl-1, d)
      integer, intent(in) :: no, nv, v0, nvb
      real(dp), intent(in) :: ovvv(no, nv, nv, nv)
      real(dp), intent(out) :: ovpack(no*nv*nv*nvb)
      integer :: t

      do concurrent(t=1:no*nv*nv*nvb)
         ovpack(t) = z2_item(t, no, nv, v0, ovvv)
      end do
   end subroutine kernel_pack_z2

   subroutine kernel_sub_z1(no, nv, v0, nvb, cmat, h)
      !! h(i,j,a,b) -= cmat(ij; b, al), a = v0 + al - 1
      integer, intent(in) :: no, nv, v0, nvb
      real(dp), intent(in) :: cmat(no*no*nv*nvb)
      real(dp), intent(inout) :: h(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nvb)
         h(sub_z1_dest(t, no, nv, v0)) = h(sub_z1_dest(t, no, nv, v0)) - cmat(t)
      end do
   end subroutine kernel_sub_z1

   subroutine kernel_sub_z2(no, nv, v0, nvb, cmat, h)
      !! h(i,j,a,b) -= cmat(ij; a, bl), b = v0 + bl - 1
      integer, intent(in) :: no, nv, v0, nvb
      real(dp), intent(in) :: cmat(no*no*nv*nvb)
      real(dp), intent(inout) :: h(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nvb)
         h(sub_z2_dest(t, no, nv, v0)) = h(sub_z2_dest(t, no, nv, v0)) - cmat(t)
      end do
   end subroutine kernel_sub_z2

   ! ---------------------------------------------------------------------
   ! the item routines: pure, `routine seq`, explicit-shape dummies, no
   ! arrays of their own. A flat 1-based index t is the position in the
   ! array that is written (or, for the scatters, read), first index fastest.
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

   pure function taut_item(t, no, nv, tau) result(v)
      !! taut(c,d,ij) = tau(i,j,c,d), position t of (nv, nv, no, no)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp) :: v
      integer :: c, d, i, j

      call split4(t, nv, nv, no, c, d, i, j)
      v = tau(i, j, c, d)
   end function taut_item

   pure function y2_item(t, no, nv, np, y) result(v)
      !! y2(ij,a,d,P) = y(a,d,ij,P), position t of (no*no, nv, nv, np)
      !$acc routine seq
      integer, intent(in) :: t, no, nv, np
      real(dp), intent(in) :: y(nv*nv*no*no*np)
      real(dp) :: v
      integer :: ij, a, d, p

      call split4(t, no*no, nv, nv, ij, a, d, p)
      v = y(a + nv*(d - 1) + nv*nv*(ij - 1) + nv*nv*no*no*(p - 1))
   end function y2_item

   pure function z1_item(t, no, nv, v0, ovvv) result(v)
      !! ovvv(k,d,a,c) at position t of (no, nv[c], nv[d], nvb[a - v0 + 1])
      !$acc routine seq
      integer, intent(in) :: t, no, nv, v0
      real(dp), intent(in) :: ovvv(no, nv, nv, nv)
      real(dp) :: v
      integer :: k, c, d, al

      call split4(t, no, nv, nv, k, c, d, al)
      v = ovvv(k, d, v0 + al - 1, c)
   end function z1_item

   pure function z2_item(t, no, nv, v0, ovvv) result(v)
      !! ovvv(k,c,b,d) at position t of (no, nv[c], nv[d], nvb[b - v0 + 1])
      !$acc routine seq
      integer, intent(in) :: t, no, nv, v0
      real(dp), intent(in) :: ovvv(no, nv, nv, nv)
      real(dp) :: v
      integer :: k, c, d, bl

      call split4(t, no, nv, nv, k, c, d, bl)
      v = ovvv(k, c, v0 + bl - 1, d)
   end function z2_item

   pure function sub_z1_dest(t, no, nv, v0) result(dest)
      !! Position in h(i,j,a,b) of cmat position t of (no*no, nv[b], nvb[a - v0 + 1])
      !$acc routine seq
      integer, intent(in) :: t, no, nv, v0
      integer :: dest
      integer :: ij, b, al, unused

      call split4(t, no*no, nv, 1, ij, b, unused, al)
      dest = ij + no*no*(v0 + al - 2) + no*no*nv*(b - 1)
   end function sub_z1_dest

   pure function sub_z2_dest(t, no, nv, v0) result(dest)
      !! Position in h(i,j,a,b) of cmat position t of (no*no, nv[a], nvb[b - v0 + 1])
      !$acc routine seq
      integer, intent(in) :: t, no, nv, v0
      integer :: dest

      dest = t + no*no*nv*(v0 - 1)
   end function sub_z2_dest

end module trc_rccsd_ladder
