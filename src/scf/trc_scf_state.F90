!
! The per-iteration arrays an SCF run writes every pass through the loop:
! the Fock/Coulomb-exchange/XC matrices in both the run's own AO space and
! the Cartesian one the kernels build in, the incremental-build history,
! the continuum's operator, and the DIIS commutator's own workspace.
! `alloc` owns the "1 x 1 when this run never uses it" sizing rule (see
! `trc_scf_driver`'s own header); `to_device`/`release_device` are the ONE
! enter/exit pair for all of it, together with the context's `hcore`/
! `smat`/`x` and the result's `dmat`/`cmo`/`eps` -- one region, because
! that is what the original single `!$acc enter data` statement did, and
! splitting it by module would only have meant more regions mapping the
! same arrays, not fewer bugs.
!
module trc_scf_state
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use trc_step_table, only: step_table_t
   use trc_scf_context, only: scf_ctx_t
   use trc_scf_settings, only: scf_settings_t, trc_scf_result_t
   implicit none
   private

   public :: scf_state_t, scf_timers_t

   type :: scf_timers_t
      !! The loop's phases, one pic timer each, started at the phase's
      !! first statement and read at its last, and their running totals
      !! (`trc_scf_context`'s `scf_ctx_t` owns the once-per-geometry setup
      !! timings). A pic timer measures one interval, so the totals over
      !! iterations are kept beside it. `tf_it`/`tx_it`/`tr_it` are this
      !! iteration's own share, not a running total (see `trc_scf_driver`'s
      !! own header on why that differs between SAD's first iteration and
      !! its plateau).
      type(timer_type) :: fock, xc, pcm, rest, run
      !! `run` times the whole of `trc_scf_run`, setup to the last iteration.
      real(dp) :: t_fock = 0.0_dp, t_xc = 0.0_dp, t_rest = 0.0_dp
      real(dp) :: t_pcm = 0.0_dp
      real(dp) :: tf_it = 0.0_dp, tx_it = 0.0_dp, tr_it = 0.0_dp
      real(dp) :: tx_pts = 0.0_dp, tx_prs = 0.0_dp
      real(dp) :: t_vv10 = 0.0_dp
      !! Running total of the VV10 calls `scf_fock_build_vxc` wraps in their own
      !! NVTX range -- see `xcst` below, and `trc_scf_fock`'s own header
      !! on why VV10's time is not folded into `t_xc`/`tx_pts`/`tx_prs`
      !! the way the semilocal functional's points/pairs split is.
      type(step_table_t) :: xcst
      !! The DFT pipeline's steps (grid build, semilocal XC, VV10) in the
      !! same shared table RI-MP2/DF-JK report through -- NVTX ranges on
      !! the two steps whose own call sites can bracket them
      !! (`scf_fock_build_vxc`'s xc_rks/uks and vv10_rks/uks calls); grid build's
      !! total is added in after the fact (`add_time`, no range) since
      !! its own timer lives inside `trc_scf_context` and threading a
      !! second timing mechanism through that file was not worth it here.
   end type scf_timers_t

   type :: scf_state_t
      real(dp), allocatable :: fock(:, :, :), gmat(:, :, :), vxc(:, :, :)
      real(dp), allocatable :: d_c(:, :, :), gmat_c(:, :, :), vxc_c(:, :, :)
      real(dp), allocatable :: dtot_c(:, :), jmat_c(:, :), w1_c(:, :), klr_c(:, :)
      !! klr_c: the long-range exchange of a range-separated hybrid, built
      !! by its own pass (erf(omega r12)/r12) and added into G.
      real(dp), allocatable :: dref_c(:, :), gacc_c(:, :), ddif_c(:, :)
      !! The incremental build's reference density/G and their difference.
      !! Never zeroed here (see `alloc`): the full build on iteration 1
      !! always writes them before an increment ever reads them.
      real(dp), allocatable :: vpcm_c(:, :), vpcm(:, :)
      !! the continuum's operator, Cartesian as its kernels make it and in
      !! the run's own space; the same for both spins
      real(dp), allocatable :: dold(:, :, :), errv(:, :, :)
      real(dp), allocatable :: w1(:, :), w2(:, :)
      real(dp), allocatable :: ofrac(:)
   contains
      procedure :: alloc => scf_state_alloc
      procedure :: to_device => scf_state_to_device
      procedure :: release_device => scf_state_release_device
      procedure :: save_dold => scf_state_save_dold
      procedure :: update_errv => scf_state_update_errv
   end type scf_state_t

contains

   subroutine scf_state_alloc(this, context, config)
      !! Sizes follow `context`/`config` alone (never `opts` directly): every array
      !! here is nao/nao_c square, or 1 x 1 when this run has no use for it
      !! (unrestricted-only, range-separated-only, incremental-only, PCM-only).
      class(scf_state_t), intent(inout) :: this
      type(scf_ctx_t), intent(in) :: context
      type(scf_settings_t), intent(in) :: config
      integer :: nao, nao_c, nspin, nxc, nxc_c, nu_c, nlr_c, nin_c, npcm, npcm_c
      nao = context%nao
      nao_c = context%nao_c
      nspin = config%nspin
      nxc = merge(nao, 1, config%dft)
      nxc_c = merge(nao_c, 1, config%dft)
      nu_c = merge(nao_c, 1, nspin == 2)
      nlr_c = merge(nao_c, 1, config%exx_lr /= 0.0_dp)
      nin_c = merge(nao_c, 1, config%incr_on)
      npcm = merge(nao, 1, config%pcm_on)
      npcm_c = merge(nao_c, 1, config%pcm_on)
      allocate (this%fock(nao, nao, nspin), this%gmat(nao, nao, nspin), this%vxc(nxc, nxc, nspin))
      allocate (this%d_c(nao_c, nao_c, nspin), this%gmat_c(nao_c, nao_c, nspin), this%vxc_c(nxc_c, nxc_c, nspin))
      allocate (this%dtot_c(nu_c, nu_c), this%jmat_c(nu_c, nu_c), this%w1_c(nao_c, nao_c), this%klr_c(nlr_c, nlr_c))
      allocate (this%dref_c(nin_c, nin_c), this%gacc_c(nin_c, nin_c), this%ddif_c(nin_c, nin_c))
      allocate (this%vpcm_c(npcm_c, npcm_c), this%vpcm(npcm, npcm))
      allocate (this%dold(nao, nao, nspin), this%errv(nao, nao, nspin))
      allocate (this%w1(nao, nao), this%w2(nao, nao))
      allocate (this%ofrac(nao))
      this%fock = 0.0_dp
      this%gmat = 0.0_dp
      this%vxc = 0.0_dp
      this%dold = 0.0_dp
      this%errv = 0.0_dp
      this%w1 = 0.0_dp
      this%w2 = 0.0_dp
      this%d_c = 0.0_dp
      this%gmat_c = 0.0_dp
      this%vxc_c = 0.0_dp
      this%dtot_c = 0.0_dp
      this%jmat_c = 0.0_dp
      this%w1_c = 0.0_dp
      this%klr_c = 0.0_dp
      this%vpcm_c = 0.0_dp
      this%vpcm = 0.0_dp
      this%ofrac = 0.0_dp
   end subroutine scf_state_alloc

   subroutine scf_state_to_device(this, context, res)
      class(scf_state_t), intent(inout) :: this
      type(scf_ctx_t), intent(in) :: context
      type(trc_scf_result_t), intent(in) :: res
      !$acc enter data copyin(this%ofrac)
      !$acc enter data copyin(context%hcore, context%smat, context%x, res%dmat, res%cmo, res%eps) &
      !$acc            copyin(this%fock, this%gmat, this%vxc, this%dold, this%errv, this%w1, this%w2) &
      !$acc            copyin(this%d_c, this%gmat_c, this%vxc_c, this%dtot_c, this%jmat_c, this%w1_c, this%klr_c) &
      !$acc            create(this%dref_c, this%gacc_c, this%ddif_c) copyin(this%vpcm_c, this%vpcm)
   end subroutine scf_state_to_device

   subroutine scf_state_release_device(this, context, res)
      class(scf_state_t), intent(inout) :: this
      type(scf_ctx_t), intent(in) :: context
      type(trc_scf_result_t), intent(in) :: res
      !$acc exit data delete(this%ofrac)
      !$acc exit data delete(context%hcore, context%smat, context%x, res%dmat, res%cmo, res%eps) &
      !$acc           delete(this%fock, this%gmat, this%vxc, this%dold, this%errv, this%w1, this%w2) &
      !$acc           delete(this%d_c, this%gmat_c, this%vxc_c, this%dtot_c, this%jmat_c, this%w1_c, this%klr_c) &
      !$acc           delete(this%dref_c, this%gacc_c, this%ddif_c, this%vpcm_c, this%vpcm)
   end subroutine scf_state_release_device

   subroutine scf_state_save_dold(this, nao, nspin, dmat)
      !! `dold` must hold `res%dmat` from BEFORE DIIS/diagonalisation touch
      !! it -- the energy this iteration reports was evaluated at `dold`,
      !! and the damping step below blends the new density back towards it.
      class(scf_state_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: dmat(nao, nao, nspin)
      call copy3(nao, nspin, dmat, this%dold)
   end subroutine scf_state_save_dold

   subroutine scf_state_update_errv(this, nao, nspin, dmat)
      !! `errv` after the step: how far the new density moved from `dold`,
      !! the RMS(D) the convergence print and (over MPI) every rank's own
      !! broadcast-settled density are read from.
      class(scf_state_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: dmat(nao, nao, nspin)
      call diff3(nao, nspin, dmat, this%dold, this%errv)
   end subroutine scf_state_update_errv

   subroutine copy3(nao, nspin, src, dst)
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: src(nao, nao, nspin)
      real(dp), intent(out) :: dst(nao, nao, nspin)
      integer :: i, j, s
      do concurrent(i=1:nao, j=1:nao, s=1:nspin)
         dst(i, j, s) = src(i, j, s)
      end do
   end subroutine copy3

   subroutine diff3(nao, nspin, a, b, c)
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: a(nao, nao, nspin), b(nao, nao, nspin)
      real(dp), intent(out) :: c(nao, nao, nspin)
      integer :: i, j, s
      do concurrent(i=1:nao, j=1:nao, s=1:nspin)
         c(i, j, s) = a(i, j, s) - b(i, j, s)
      end do
   end subroutine diff3

end module trc_scf_state
