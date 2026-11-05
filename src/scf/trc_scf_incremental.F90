!
! When the SCF builds its two-electron matrix from the density DIFFERENCE
! and when it builds it in full: the schedule of the incremental Fock build
! (`trc_scf_fock`'s build_g_arrays does the arithmetic, this module only
! decides). One `scf_incremental_t` per run, asked once per iteration before
! the build (`next_build`) and once after the convergence test
! (`accept`). Everything it reads is the commutator norm and the
! convergence test's own outcome, which every rank already holds as rank
! 0's broadcast copy, so every rank takes the same decision without a
! collective of its own.
!
! THE ERROR ARGUMENT. G is linear in D, so G(D_n) = G(D_ref) + G(D_n -
! D_ref) exactly; what is not exact is the screening. A quartet the kernels
! drop changes a G element by less than the threshold tau times its
! degeneracy in a full build and in an increment alike (the density test
! is the same, on D or on D_n - D_ref), but the increments' errors ADD UP:
! each increment drops a different set, and the accumulated matrix carries
! the sum of all of them since the last full build. That sum is an offset
! eps in the Fock matrix the iteration sees, so the SCF heads for the fixed
! point of F + eps, not of F, and |FDS - SDF| stops falling where eps is.
! Four rules keep that from reaching the answer:
!
! 0. A looser screen only while eps cannot matter. Increments are screened
!    at `early_scale` times tau while |FDS - SDF| >= `early_until`, at
!    `late_scale` times tau after, and the first build past `early_until`
!    is full, which discards the loose increments' eps. Measured at 100
!    times tau 3.16e-11, the floor that eps sets is 2e-5 to 5e-5 (water20,
!    gly10, gly20 cc-pVTZ), two decades below the 1e-2 that ends the phase.
! 1. Convergence is declared only on a FULL build. When the test is met on
!    an increment, the next build is full and the test is taken again
!    there, against the true F: the commutator it sees is the commutator
!    of the exact (fully screened) Fock matrix, and the energy is the
!    full-build energy at that density. The energy is stationary, so the
!    density's own error from eps moves it only at second order.
! 2. The iteration that is expected to converge is built in full to begin
!    with, so rule 1 costs no extra iteration when the expectation holds:
!    DIIS shrinks the commutator geometrically, and when the last ratio
!    carries it below `conv_diis` the next build is full. When it falls
!    faster than that, the run takes one more iteration than full builds
!    would (w20 cc-pVQZ 12 against 11; water20 aug-cc-pVDZ 13 against 12,
!    where the jump onto the full build also exceeds conv_energy 1e-10).
!    Predicting to 3 x conv_diis (TRC_INCR_PREDICT=3) or building in full
!    below 10 x conv_diis (TRC_INCR_ENDGAME=10) removes those, but spends a
!    full build an iteration early elsewhere: over gly10, gly30, water20,
!    w20, w172 and the aug-cc-pVDZ case 441 s as is, 454 s with either.
! 3. A stall is a rebuild. Once the commutator stops falling over
!    STALL_ITERATIONS increments, it has reached the floor eps sets, and
!    only a full build lowers that floor; it also resets the sum.
!
! At `late_scale` 1 eps is about a hundred times below the early phase's,
! under conv_diis's 1e-5, and the full builds of rules 0 to 2 are the only
! ones a converging run takes: three or four per SCF, against one per
! iteration (trc_scf_options_t%incremental has the timings).
!
module trc_scf_incremental
   use pic_types, only: dp
   use trc_scf_settings, only: scf_settings_t
   implicit none
   private

   public :: scf_incremental_t

   integer, parameter :: STALL_ITERATIONS = 2
   !! Consecutive increments without a new lowest commutator that count as
   !! a stall. One is DIIS's ordinary non-monotone step; two in a row, and
   !! the increments' screening floor is the likelier reason.

   type :: scf_incremental_t
      logical :: on = .false.
      real(dp) :: start = huge(1.0_dp)
      !! Increment only once the previous commutator is below this.
      integer :: reset = 0
      !! At most this many increments in a row, then a full build; 0 is no limit.
      real(dp) :: late_scale = 1.0_dp
      real(dp) :: early_scale = 1.0_dp
      real(dp) :: early_until = 0.0_dp
      !! Increments are screened at `early_scale` times the threshold while
      !! the previous commutator is at least `early_until`, at `late_scale`
      !! times after (this module's own header).
      real(dp) :: scale = 0.0_dp
      !! This iteration's threshold scale, 0 for a full build: what
      !! `scf_fock_build_g` takes.
      logical :: early_since_full = .false.
      !! An early (loose) increment has gone into G since the last full build.
      real(dp) :: conv_diis = 0.0_dp
      real(dp) :: predict = 1.0_dp
      real(dp) :: endgame = 0.0_dp
      !! Full builds once the commutator is within this factor of `conv_diis`
      !! (an experiment hook, scf_settings_t%incr_endgame_run; 0 is off).
      logical :: full_due = .false.
      !! The next build must be full: a convergence to verify, a predicted
      !! convergence, or a stall (this module's own header).
      logical :: built_full = .true.
      !! This iteration's G came from a full build.
      integer :: since_full = 0
      !! Increments since the last full build.
      integer :: stalled = 0
      real(dp) :: best = huge(1.0_dp)
      !! The smallest commutator since the last full build.
      real(dp) :: err_prev = huge(1.0_dp)
      integer :: n_increments = 0
      integer :: n_full = 0
   contains
      procedure :: setup => scf_incremental_setup
      procedure :: next_build => scf_incremental_next_build
      procedure :: accept => scf_incremental_accept
   end type scf_incremental_t

contains

   subroutine scf_incremental_setup(this, config)
      class(scf_incremental_t), intent(out) :: this
      type(scf_settings_t), intent(in) :: config
      this%on = config%incr_on
      this%start = config%incr_start_run
      this%reset = config%incr_reset_run
      this%late_scale = config%incr_scale_run
      this%early_scale = config%incr_early_scale_run
      this%early_until = config%incr_early_until_run
      this%conv_diis = config%conv_diis
      this%endgame = config%incr_endgame_run
      this%predict = config%incr_predict_run
   end subroutine scf_incremental_setup

   logical function scf_incremental_next_build(this, iteration, errmax) result(incr_now)
      !! Whether iteration `iteration` builds G from the density difference.
      !! `errmax` is the previous iteration's commutator norm (huge before
      !! the first). The first build is always full: there is no reference.
      class(scf_incremental_t), intent(inout) :: this
      integer, intent(in) :: iteration
      real(dp), intent(in) :: errmax
      logical :: early
      early = errmax >= this%early_until
      incr_now = this%on .and. iteration > 1 .and. .not. this%full_due .and. errmax < this%start
      if (this%reset > 0) incr_now = incr_now .and. this%since_full < this%reset
      ! Leaving the early phase discards what its loose increments left in G.
      if (.not. early .and. this%early_since_full) incr_now = .false.
      this%full_due = .false.
      this%built_full = .not. incr_now
      if (incr_now) then
         this%since_full = this%since_full + 1
         this%n_increments = this%n_increments + 1
         this%scale = merge(this%early_scale, this%late_scale, early)
         if (early) this%early_since_full = .true.
      else
         this%scale = 0.0_dp
         this%since_full = 0
         this%stalled = 0
         this%best = huge(1.0_dp)
         this%n_full = this%n_full + 1
         this%early_since_full = .false.
      end if
   end function scf_incremental_next_build

   logical function scf_incremental_accept(this, met, errmax) result(converged)
      !! The convergence test's outcome `met` for this iteration, and this
      !! iteration's commutator norm: true only when `met` holds on a full
      !! build (rule 1 of this module's header), and the next build's kind
      !! decided from here on (rules 1 to 3). With increments off every
      !! build is full and this is `met`.
      class(scf_incremental_t), intent(inout) :: this
      logical, intent(in) :: met
      real(dp), intent(in) :: errmax
      real(dp) :: ratio
      converged = met .and. this%built_full
      if (.not. this%on .or. converged) return
      if (met) this%full_due = .true.
      ! A NaN commutator compares false everywhere below and leaves the
      ! decision alone; the loop's own test never calls NaN converged.
      ! No ratio before there are two commutators to take it from.
      ratio = 1.0_dp
      if (this%err_prev < huge(1.0_dp)) ratio = min(1.0_dp, errmax/max(this%err_prev, tiny(1.0_dp)))
      if (errmax*ratio < this%predict*this%conv_diis) this%full_due = .true.
      if (errmax < this%endgame*this%conv_diis) this%full_due = .true.
      if (.not. this%built_full) then
         if (errmax < this%best) then
            this%best = errmax
            this%stalled = 0
         else
            this%stalled = this%stalled + 1
            if (this%stalled >= STALL_ITERATIONS) this%full_due = .true.
         end if
      end if
      this%err_prev = errmax
   end function scf_incremental_accept

end module trc_scf_incremental
