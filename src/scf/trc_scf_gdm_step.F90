module trc_scf_gdm_step
   !! Geometric direct minimisation (GDM) as a step mode of `trc_scf_run`'s
   !! loop: Van Voorhis and Head-Gordon, Mol. Phys. 100, 1713 (2002).
   !!
   !! The loop's iteration is "Fock build at the current density, energy,
   !! commutator, new density". GDM replaces only the last part, the root-only
   !! "DIIS + level shift + diagonalise", with a step in the orbital rotation
   !! space, so everything the loop already supports (RHF/UHF, DFT with VV10,
   !! PCM, density fitting, incremental Fock builds, linear-dependence
   !! removal, several ranks) works unchanged: the step never builds a Fock
   !! matrix, it is handed the one the loop built. One Fock build per
   !! iteration, as with DIIS. The numerical pieces (L-BFGS, preconditioner,
   !! accept rule, step cap) are `trc_scf_gdm`; the orbital-rotation device
   !! kernels are `trc_scf_soscf`'s.
   !!
   !! THE ITERATION. The Fock matrix and energy handed in belong to the TRIAL
   !! orbitals C_try the previous call left in `scf_result%cmo`. Then
   !!
   !!  1. the first call accepts the trial unconditionally;
   !!  2. later calls judge it (`scf_gdm_step_judge`): accepted when
   !!     E_try <= E_acc + slack. The slack is the option's (`gdm_slack`,
   !!     1e-10) whenever both energies come from FULL Fock builds, which are
   !!     deterministic; it is the larger of that and `GDM_NOISE_BOUND` when
   !!     either E_acc or E_try came from an INCREMENTAL build, whose
   !!     accumulated screening error is a real noise source (the driver
   !!     passes `incr_build`, and the accepted point's kind is kept in
   !!     `acc_incr`). A negative `gdm_slack` is an explicit demand for a
   !!     decrease and is never widened. No other excuse exists: a step that
   !!     promises a tiny first-order gain is still judged by the same rule
   !!     (the earlier "predicted gain below 1e-9" excuse accepted real rises
   !!     of 8e-9 on a flat mode and let the energy random-walk);
   !!  3. rejected: the history is untouched, the step scale halves and
   !!     C_try = C_acc exp(scale kappa(d)) is rebuilt from the ACCEPTED
   !!     orbitals along the SAME direction. After eleven halvings the judge
   !!     reports breakdown: the history is dropped and the preconditioned
   !!     steepest-descent direction at the accepted point is tried with a
   !!     fresh scale; if that breaks down too the run stops, not converged.
   !!     When the direction that broke down WAS the steepest descent (empty
   !!     history, or the non-descent fallback of step 4: `restarted`) there is
   !!     nothing else to try and the run stops at once;
   !!  4. accepted: the gradient g_new = hessian_factor F_ai in the frame of
   !!     C_try is formed, which becomes C_acc; the pair s = scale d (the step
   !!     actually taken, d already capped), y = g_new - g_acc goes into the
   !!     L-BFGS history; d = -H g_new, replaced by preconditioned steepest
   !!     descent (and the history reset) if it is not a descent direction;
   !!     d is capped so no rotation exceeds the maximum EVERY time; and
   !!     C_try = C_acc exp(kappa(d));
   !!  5. the density is rebuilt from the occupied columns of C_try.
   !!
   !! THE MOVING FRAME. d, g and the history live in the orbital frame of the
   !! point they were formed at. s and y are compared in successive frames
   !! with no transport between them (GDM's own approximation, the paper's
   !! section 3): since each step is small the frames are close, and the
   !! curvature test of the history (`scf_gdm_curvature_tol`) drops a pair
   !! that the approximation spoils.
   !!
   !! VECTOR LAYOUT. One host vector for the whole problem: per spin block the
   !! occupied-virtual elements of C^T F C at index (i-1)*nvir + a, spins
   !! concatenated (alpha then beta); one history for both spins since there
   !! is one energy. The block of one spin is, in memory, the column-major
   !! array (nvir, nocc), which is how it goes to the device kernels.
   !! Everything lives in the retained orbital space n_mo = nao - n_lindep. A
   !! spin block with no occupied or no virtual orbital has nothing to rotate.
   !!
   !! THE ENERGY NOISE. The accept slack (`trc_scf_settings` default
   !! `DEFAULT_GDM_SLACK`, 1e-10) is the convergence threshold: a rise below
   !! it cannot matter, and with full builds (no noise source at all) it is
   !! the whole rule. The incremental Fock build's accumulated error has been
   !! measured at ~1e-7 Eh on a hard case in the loose early phase and about
   !! a hundred times less late (trc_scf_incremental, rule 0); `GDM_NOISE_BOUND`
   !! covers the late value, and only comparisons that involve an incremental
   !! energy get it. A REJECTED step requests a full Fock build for the next
   !! iteration (the driver sets `full_due`), so the retry is judged on a
   !! build free of accumulated increments' error; E_acc may still be an
   !! incremental energy, which is why the wide bound follows E_acc too.
   !!
   !! HANDOVER AND THE END. `scf_gdm_mode` latches: the first iteration past the
   !! damped start whose commutator is below `gdm_start` is a plain
   !! undamped, unshifted, unextrapolated diagonalisation (so the density it
   !! leaves is exactly the idempotent density of `cmo`); GDM steps run from
   !! the next one, never handing back. The occupation of that diagonalisation
   !! is the one GDM keeps (a rotation never changes the number of occupied
   !! orbitals): aufbau, or, when MOM had latched before the handover, maximum
   !! overlap (the driver keeps ORB_MOM for this one iteration).
   !!
   !! MOM AND GDM. A state MOM holds off the aufbau occupation (an excited
   !! state) is a saddle point of the energy in orbital-rotation space, which a
   !! minimiser leaves: measured on water cc-pVDZ with HOMO and LUMO swapped,
   !! GDM reaches the MOM state to 1e-11 Eh and, with a numerical seed along the
   !! one downhill direction, runs away from it to the ground state in ten more
   !! iterations, converged, no warning. So when the handover diagonalisation
   !! had MOM occupy a non-aufbau set (`mom_swaps` > 0 on either spin) GDM
   !! DECLINES (`scf_gdm_decline`, `declined`): it never takes a step and the
   !! run continues as a `mom` run. If MOM merely latched on the aufbau
   !! occupation (a stall at a minimum) GDM goes on.
   !!
   !! STARTED TOO EARLY. A minimiser started far from convergence keeps the
   !! occupation of its first diagonalisation and cannot see negative
   !! curvature, so it can end on an excited or saddle solution: after the
   !! final semicanonicalisation `scf_gdm_finish` counts, per spin, the virtual
   !! orbitals below the highest occupied one (`off_aufbau`, in
   !! `scf_result%gdm_off_aufbau`); non-zero is a WARNING. `scf_gdm_finish` leaves in
   !! `scf_result` the SEMICANONICAL orbitals (occupied-occupied and
   !! virtual-virtual blocks of the Fock matrix diagonalised) of the last
   !! ACCEPTED point, its density in `scf_state%dold` and its energies in
   !! `scf_energies` (the apparent surface charge included: the driver stores
   !! it in `scf_energies%pcm_qtotal` at every Fock build), never a rejected
   !! trial's.
   !!
   !! RANKS. Only the root rank holds the state and takes the step. Every rank
   !! calls `scf_gdm_mode`, `scf_gdm_report` and `scf_gdm_finish`, which
   !! need nothing else: `rejected` and `stopped` reach the others through the
   !! driver's existing root-to-all flag broadcast, and `scf_gdm_finish`
   !! broadcasts the accepted density and energies once at the end.
   !!
   !! DEVICE RESIDENCY. The accepted orbitals and Fock matrix, the work
   !! matrices of the rotation and of the gradient, and two small vectors stay
   !! on the device between iterations (`enter data` in `gdm_setup`, `exit
   !! data` in `scf_gdm_release`, which the driver calls on every path out of
   !! the loop). The only things that cross the bus per accepted step are the
   !! gradient and the diagonal Fock elements (device to host, nocc*nvir and
   !! n_mo numbers), and the rotation amplitudes (host to device).
   use pic_types, only: dp
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_blas_interfaces, only: pic_dot
   use pic_mpi_lib, only: comm_t
   use trc_error, only: error_t
   use trc_linalg, only: trc_linalg_t
   use trc_cphf, only: cphf_copy_matrix
   use trc_scf_context, only: scf_ctx_t, scf_context_bcast_flat
   use trc_scf_settings, only: scf_settings_t, trc_scf_result_t, scf_energies_t, scf_settings_count_off_aufbau
   use trc_scf_state, only: scf_state_t
   use trc_scf_gdm, only: trc_scf_gdm_history_t, scf_gdm_preconditioner_build, scf_gdm_step_judge, &
                          scf_gdm_step_limit, scf_gdm_is_descent
   use trc_scf_soscf, only: trc_soscf_dims_t, scf_soscf_rotation_matrix, scf_soscf_semicanonicalize, &
                            scf_soscf_diagonal, scf_soscf_copy_eps, scf_soscf_ov_one_norm
   implicit none
   private

   public :: trc_scf_gdm_step_t, GDM_OFF, GDM_HANDOVER, GDM_ITERATE, GDM_STOP_MESSAGE
   public :: scf_gdm_step_accept_slack, GDM_NOISE_BOUND

   integer, parameter :: GDM_OFF = 0, GDM_HANDOVER = 1, GDM_ITERATE = 2
   !! What `scf_gdm_mode` says about an iteration: not GDM's, the plain
   !! diagonalisation that hands over, a GDM step.

   character(len=*), parameter :: GDM_STOP_MESSAGE = "trc_scf: GDM found no downhill step (the energy did not "// &
                                                    "fall along the quasi-Newton or the steepest-descent direction, "// &
                                                    "or a gradient was not finite); stopping, not converged"

   real(dp), parameter :: GDM_NOISE_BOUND = 1.0e-8_dp
   !! The accept slack for a comparison in which E_acc or E_try came from an
   !! INCREMENTAL Fock build (never otherwise, never for a negative
   !! `gdm_slack`): the build's accumulated screening error late in a run
   !! (trc_scf_incremental: ~1e-7 early, about a hundred times less once the
   !! screen tightens) with margin. Full builds are deterministic and get the
   !! option's plain slack.
   !! A virtual orbital counts as below the highest occupied one only when it
   !! is lower by more than this, hartree: a degenerate pair split by
   !! round-off (an open-shell pi pair) is not a non-aufbau state. The same
   !! number the post-SCF gate uses, so the warning and the refusal agree.
   real(dp), parameter :: GDM_HESSIAN_FLOOR = 0.1_dp
   !! The smallest Hessian estimate hessian_factor*(e_a - e_i) the
   !! preconditioner divides by, Eh per radian squared (a gap of 25 mEh
   !! restricted, 50 mEh unrestricted). Below it a near-degenerate pair would
   !! make its own component dominate the direction and the rotation cap
   !! would then shrink every other component with it.
   integer, parameter :: N_ENERGIES = 7
   !! The numbers of `scf_energies_t`, as `scf_gdm_finish` broadcasts them,
   !! followed by the off-aufbau count.
   integer, parameter :: GDM_HISTORY = 12
   !! L-BFGS pairs kept; 2*n*12 doubles on the host, n = sum of nocc*nvir.

   type :: trc_scf_gdm_step_t
      !! The state of the GDM mode across iterations. Plain `type ::` filled by
      !! `gdm_setup`, never by a structure constructor (allocatable components,
      !! see the header of trc_scf_soscf). `active`, `rejected` and `stopped`
      !! are known on every rank; the rest is the root rank's alone.
      logical :: active = .false.
         !! The latch: GDM steps run from the iteration after the handover on.
      logical :: rejected = .false.
         !! This iteration's trial was rejected (it must not converge the run).
      logical :: stopped = .false.
         !! No downhill step (or a non-finite gradient): stop, not converged.
      logical :: ready = .false.
         !! `gdm_setup` has run: the arrays below are allocated and mapped.
      logical :: trivial = .false.
         !! Nothing to rotate in either spin: the plain orbitals are the answer.
      logical :: restarted = .false.
         !! The current direction is steepest descent (empty history, the
         !! non-descent fallback, or a breakdown's restart): a breakdown of it
         !! stops the run.
      logical :: acc_incr = .false.
         !! The accepted energy came from an incremental Fock build.
      logical :: noisy = .false.
         !! The comparison being judged involves an incremental build (the
         !! trial's or the accepted point's).
      logical :: declined = .false.
         !! GDM gave way to MOM at the handover (this module's header): latched,
         !! `scf_gdm_mode` says `GDM_OFF` from then on. Known on every rank.
      integer :: nao = 0, n_mo = 0, nspin = 1, n_total = 0, n_iterations = 0
      integer :: nocc(2) = 0, nvir(2) = 0, offset(2) = 0
      real(dp) :: hessian_factor = 4.0_dp
      real(dp) :: e_acc = 0.0_dp
         !! Energy at the accepted point.
      real(dp) :: max_rise = 0.0_dp
         !! Largest energy rise between consecutive accepted points.
      real(dp) :: max_rise_full = 0.0_dp
         !! The same over the pairs of accepted points that both came from
         !! full Fock builds (no noise source): at most the plain slack.
      real(dp) :: max_step = 0.0_dp
         !! Largest rotation element applied in any trial of the run, radians.
      integer :: off_aufbau = 0
         !! Virtual orbitals below the highest occupied one, over the spins, at the end.
      real(dp) :: step_scale = 1.0_dp
         !! Scale of the trial now being tried (the judge halves or doubles it).
      real(dp) :: scale_used = 1.0_dp
         !! Scale of the trial the judge is looking at.
      type(scf_energies_t) :: energies_acc
         !! The energy bundle at the accepted point.
      character(len=160) :: line = ""
      type(trc_scf_gdm_history_t) :: history
      real(dp), allocatable :: g_acc(:), g_new(:), h0_inv(:), direction(:)
      real(dp), allocatable :: vec_s(:), vec_y(:)
      real(dp), allocatable :: step(:)
         !! Rotation amplitudes of one spin block (host and device copies).
      real(dp), allocatable :: c_acc(:, :, :), fock_acc(:, :, :)
         !! Accepted orbitals (nao, n_mo, nspin) and Fock matrix (nao, nao, nspin), device.
      real(dp), allocatable :: c_work(:, :), work(:, :), fock_mo(:, :)
      real(dp), allocatable :: kappa(:, :), rotation(:, :), scratch_a(:, :), scratch_b(:, :), scratch_c(:, :)
      real(dp), allocatable :: energies_mo(:), grad(:)
   contains
      procedure :: mode => scf_gdm_mode
      procedure :: iterate => scf_gdm_iterate
      procedure :: finish => scf_gdm_finish
      procedure :: release => scf_gdm_release
      procedure :: report => scf_gdm_report
      procedure :: decline => scf_gdm_decline
   end type trc_scf_gdm_step_t

contains

   function scf_gdm_mode(this, scf_config, iteration, errmax) result(mode)
      !! What iteration `iteration` is: `GDM_OFF`, the plain diagonalisation
      !! that hands over (`GDM_HANDOVER`, which latches `active`), or a GDM
      !! step (`GDM_ITERATE`). Every rank asks, with the same broadcast
      !! `errmax`, so every rank agrees. Also clears last iteration's flags.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_settings_t), intent(in) :: scf_config
      integer, intent(in) :: iteration
      real(dp), intent(in) :: errmax
      integer :: mode
      mode = GDM_OFF
      this%rejected = .false.
      this%stopped = .false.
      if (.not. scf_config%gdm .or. this%declined) return
      if (this%active) then
         mode = GDM_ITERATE
      else if (iteration > scf_config%ndamp .and. errmax < scf_config%gdm_start) then
         mode = GDM_HANDOVER
         this%active = .true.
      end if
   end function scf_gdm_mode

   pure function scf_gdm_step_accept_slack(slack_option, trial_incremental, accepted_incremental) result(slack)
      !! The energy rise a trial may show and still be accepted. The option's
      !! slack when both energies come from full Fock builds (deterministic);
      !! at least `GDM_NOISE_BOUND` when either came from an incremental build.
      !! A negative option (an explicit demand for a decrease) is never widened.
      real(dp), intent(in) :: slack_option
      logical, intent(in) :: trial_incremental, accepted_incremental
      real(dp) :: slack
      slack = slack_option
      if ((trial_incremental .or. accepted_incremental) .and. slack_option >= 0.0_dp) slack = max(slack, GDM_NOISE_BOUND)
   end function scf_gdm_step_accept_slack

   subroutine scf_gdm_decline(this, mom_swaps)
      !! Called by every rank after the handover diagonalisation with the
      !! orbitals MOM occupied above the aufbau set, per spin (the driver's
      !! broadcast copy): any of them means MOM holds an excited state, a
      !! saddle point a minimiser would leave, so GDM stands down for the run.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      integer, intent(in) :: mom_swaps(2)
      if (any(mom_swaps > 0)) then
         this%declined = .true.
         this%active = .false.
      end if
   end subroutine scf_gdm_decline

   subroutine scf_gdm_iterate(this, scf_context, scf_config, scf_state, scf_result, scf_energies, la, incr_build)
      !! One GDM step (root rank): judge the trial the Fock matrix and energy
      !! in `scf_state`/`scf_energies` belong to, and leave the next trial's
      !! orbitals and density in `scf_result%cmo`/`%dmat`. Sets `rejected`
      !! and `stopped`. `incr_build`: that Fock matrix came from an incremental
      !! build. This module's header is the algorithm.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(scf_energies_t), intent(in) :: scf_energies
      type(trc_linalg_t), intent(inout) :: la
      logical, intent(in) :: incr_build
      logical :: first, accepted, broken
      real(dp) :: slack
      character(len=8) :: status

      this%rejected = .false.
      this%stopped = .false.
      first = .not. this%ready
      if (first) call gdm_setup(this, scf_context, scf_config)
      this%n_iterations = this%n_iterations + 1
      accepted = .true.
      broken = .false.
      if (.not. first) then
         this%scale_used = this%step_scale
         this%noisy = incr_build .or. this%acc_incr
         slack = scf_gdm_step_accept_slack(scf_config%gdm_slack, incr_build, this%acc_incr)
         call scf_gdm_step_judge(this%e_acc, scf_energies%total, slack, this%step_scale, accepted, broken)
      end if
      if (accepted) then
         call gdm_accept(this, scf_config, scf_state, scf_result, scf_energies, la, first, incr_build)
      else
         this%rejected = .true.
         call gdm_reject(this, scf_config, broken)
      end if
      status = merge("accepted", "rejected", accepted)
      if (this%stopped) status = "stopped "
      write (this%line, "(a,a,a,f8.5,a,i0,a,es10.2)") "         gdm: ", status, "  scale ", this%step_scale, &
         "  history ", this%history%npairs(), "  |g| ", sqrt(pic_dot(this%g_acc, this%g_acc))
      if (this%stopped .or. this%trivial) return
      call gdm_rotate(this, scf_config, scf_result, la)
   end subroutine scf_gdm_iterate

   subroutine gdm_setup(this, scf_context, scf_config)
      !! Allocate and map the state, on the first GDM iteration (root rank).
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(error_t) :: err
      integer :: s, nmax

      this%nao = scf_context%nao
      this%n_mo = scf_context%nao - scf_context%n_lindep
      this%nspin = scf_config%nspin
      this%hessian_factor = merge(4.0_dp, 2.0_dp, this%nspin == 1)
      this%nocc = 0
      this%nvir = 0
      this%offset = 0
      this%n_total = 0
      do s = 1, this%nspin
         this%nocc(s) = min(scf_config%nocc(s), this%n_mo)
         this%nvir(s) = this%n_mo - this%nocc(s)
         if (this%nocc(s) < 1) this%nvir(s) = 0
         this%offset(s) = this%n_total
         this%n_total = this%n_total + this%nocc(s)*this%nvir(s)
      end do
      this%trivial = this%n_total == 0
      nmax = max(this%n_total, 1)
      allocate (this%g_acc(nmax), this%g_new(nmax), this%h0_inv(nmax), this%direction(nmax))
      allocate (this%vec_s(nmax), this%vec_y(nmax))
      this%g_acc = 0.0_dp
      this%g_new = 0.0_dp
      this%h0_inv = 1.0_dp
      this%direction = 0.0_dp
      if (.not. this%trivial) call this%history%init(this%n_total, GDM_HISTORY, err)
      nmax = max(maxval(this%nocc(1:this%nspin)*this%nvir(1:this%nspin)), 1)
      allocate (this%step(nmax), this%grad(nmax), this%energies_mo(this%n_mo))
      allocate (this%c_acc(this%nao, this%n_mo, this%nspin), this%fock_acc(this%nao, this%nao, this%nspin))
      allocate (this%c_work(this%nao, this%n_mo), this%work(this%nao, this%n_mo), this%fock_mo(this%n_mo, this%n_mo))
      allocate (this%kappa(this%n_mo, this%n_mo), this%rotation(this%n_mo, this%n_mo))
      allocate (this%scratch_a(this%n_mo, this%n_mo), this%scratch_b(this%n_mo, this%n_mo), &
                this%scratch_c(this%n_mo, this%n_mo))
      this%step = 0.0_dp
      !$acc enter data create(this%c_acc, this%fock_acc, this%c_work, this%work, this%fock_mo) &
      !$acc            create(this%kappa, this%rotation, this%scratch_a, this%scratch_b, this%scratch_c) &
      !$acc            create(this%energies_mo, this%grad) copyin(this%step)
      this%ready = .true.
      this%step_scale = 1.0_dp
      this%restarted = .false.
   end subroutine gdm_setup

   subroutine scf_gdm_release(this)
      !! Unmap and free everything `gdm_setup` made; a no-op when it never ran.
      !! The driver calls it on every path out of the loop.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      if (.not. this%ready) return
      !$acc exit data delete(this%c_acc, this%fock_acc, this%c_work, this%work, this%fock_mo) &
      !$acc           delete(this%kappa, this%rotation, this%scratch_a, this%scratch_b, this%scratch_c) &
      !$acc           delete(this%energies_mo, this%grad, this%step)
      deallocate (this%c_acc, this%fock_acc, this%c_work, this%work, this%fock_mo)
      deallocate (this%kappa, this%rotation, this%scratch_a, this%scratch_b, this%scratch_c)
      deallocate (this%energies_mo, this%grad, this%step)
      deallocate (this%g_acc, this%g_new, this%h0_inv, this%direction, this%vec_s, this%vec_y)
      call this%history%release()
      this%ready = .false.
   end subroutine scf_gdm_release

   subroutine gdm_accept(this, scf_config, scf_state, scf_result, scf_energies, la, first, incr_build)
      !! The trial is the new accepted point: save it, form its gradient,
      !! push the (s, y) pair, and choose the next direction.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(in) :: scf_state
      type(trc_scf_result_t), intent(in) :: scf_result
      type(scf_energies_t), intent(in) :: scf_energies
      type(trc_linalg_t), intent(inout) :: la
      logical, intent(in) :: first, incr_build
      type(error_t) :: err
      logical :: pushed
      integer :: s

      if (.not. first) then
         this%max_rise = max(this%max_rise, scf_energies%total - this%e_acc)
         if (.not. this%noisy) then
            this%max_rise_full = max(this%max_rise_full, scf_energies%total - this%e_acc)
         end if
      end if
      this%acc_incr = incr_build
      this%e_acc = scf_energies%total
      this%energies_acc = scf_energies
      do s = 1, this%nspin
         call cphf_copy_matrix(this%nao, this%n_mo, scf_result%cmo(:, 1:this%n_mo, s), this%c_acc(:, :, s))
         call cphf_copy_matrix(this%nao, this%nao, scf_state%fock(:, :, s), this%fock_acc(:, :, s))
      end do
      this%restarted = .false.
      if (this%trivial) return
      call gdm_gradient(this, la)
      if (this%stopped .or. .not. gdm_is_finite(this%g_new(1:this%n_total))) then
         this%stopped = .true.
         return
      end if
      if (.not. first) then
         ! The step actually taken (scale * the capped direction) and the
         ! change of the gradient over it, each in its own frame: the moving
         ! frame approximation of this module's header.
         this%vec_s(1:this%n_total) = this%scale_used*this%direction(1:this%n_total)
         this%vec_y(1:this%n_total) = this%g_new(1:this%n_total) - this%g_acc(1:this%n_total)
         call this%history%push(this%vec_s(1:this%n_total), this%vec_y(1:this%n_total), pushed, err)
         if (err%has_error()) then
            this%stopped = .true.
            return
         end if
      end if
      this%g_acc(1:this%n_total) = this%g_new(1:this%n_total)
      call this%history%direction(this%g_acc(1:this%n_total), this%h0_inv(1:this%n_total), &
                                  this%direction(1:this%n_total), err)
      if (err%has_error()) then
         this%stopped = .true.
         return
      end if
      if (.not. scf_gdm_is_descent(this%g_acc(1:this%n_total), this%direction(1:this%n_total))) then
         call this%history%reset()
         call gdm_steepest_descent(this, scf_config)
      else
         call gdm_cap_direction(this, scf_config)
      end if
      ! With nothing in the history the direction is the steepest descent too.
      if (this%history%npairs() == 0) this%restarted = .true.
   end subroutine gdm_accept

   subroutine gdm_reject(this, scf_config, broken)
      !! The trial did not lower the energy: keep the history, and either
      !! retry the same direction at the halved scale the judge chose or,
      !! on breakdown, restart from steepest descent (stop if that was it).
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_settings_t), intent(in) :: scf_config
      logical, intent(in) :: broken
      if (.not. broken) return
      call this%history%reset()
      if (this%restarted) then
         this%stopped = .true.
         return
      end if
      call gdm_steepest_descent(this, scf_config)
   end subroutine gdm_reject

   subroutine gdm_steepest_descent(this, scf_config)
      !! d = -H0 g at the accepted point, capped; marks the direction as the
      !! steepest descent (`restarted`).
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_settings_t), intent(in) :: scf_config
      this%restarted = .true.
      this%direction(1:this%n_total) = -this%h0_inv(1:this%n_total)*this%g_acc(1:this%n_total)
      call gdm_cap_direction(this, scf_config)
   end subroutine gdm_steepest_descent

   subroutine gdm_cap_direction(this, scf_config)
      !! Scale the direction so no rotation exceeds the maximum, every time;
      !! a non-finite direction stops the run.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_settings_t), intent(in) :: scf_config
      real(dp) :: factor
      logical :: valid
      call scf_gdm_step_limit(this%direction(1:this%n_total), scf_config%gdm_max_rotation, factor, valid)
      if (.not. valid) then
         this%stopped = .true.
         return
      end if
      this%direction(1:this%n_total) = factor*this%direction(1:this%n_total)
   end subroutine gdm_cap_direction

   logical function gdm_is_finite(x) result(finite)
      !! False if any element is NaN or Inf (every comparison is false for NaN).
      real(dp), intent(in) :: x(:)
      finite = all(abs(x) <= huge(1.0_dp))
   end function gdm_is_finite

   subroutine gdm_gradient(this, la)
      !! The gradient hessian_factor*F_ai and the preconditioner at the
      !! accepted point, from C_acc and F_acc, into `g_new` and `h0_inv`:
      !! C^T F C per spin, its occupied-virtual block packed on the device
      !! (index (i-1)*nvir + a) and pulled with the diagonal.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: la
      type(error_t) :: err
      integer :: s, n, nocc, nvir, off

      do s = 1, this%nspin
         nocc = this%nocc(s)
         nvir = this%nvir(s)
         n = nocc*nvir
         if (n == 0) cycle
         off = this%offset(s)
         call la%gemm_mat(this%fock_acc(:, :, s), this%c_acc(:, :, s), this%work)
         call la%gemm_mat(this%c_acc(:, :, s), this%work, this%fock_mo, transa="T")
         call gdm_pack(this%n_mo, nocc, nvir, this%hessian_factor, this%fock_mo, size(this%grad), this%grad)
         call scf_soscf_diagonal(this%n_mo, this%fock_mo, this%energies_mo)
         !$acc update self(this%grad, this%energies_mo)
         this%g_new(off + 1:off + n) = this%grad(1:n)
         call scf_gdm_preconditioner_build(this%energies_mo(1:nocc), this%energies_mo(nocc + 1:this%n_mo), &
                                           this%hessian_factor, GDM_HESSIAN_FLOOR, this%h0_inv(off + 1:off + n), err)
         if (err%has_error()) this%stopped = .true.
      end do
   end subroutine gdm_gradient

   subroutine gdm_rotate(this, scf_config, scf_result, la)
      !! C_try = C_acc exp(step_scale kappa(d)) per spin, into `scf_result%cmo`,
      !! and the density from its occupied columns into `scf_result%dmat`.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_settings_t), intent(in) :: scf_config
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(trc_linalg_t), intent(inout) :: la
      real(dp) :: one_norm
      integer :: s, n, nocc, nvir, off

      do s = 1, this%nspin
         nocc = this%nocc(s)
         nvir = this%nvir(s)
         n = nocc*nvir
         off = this%offset(s)
         if (n == 0) then
            call cphf_copy_matrix(this%nao, this%n_mo, this%c_acc(:, :, s), this%c_work)
         else
            this%step(1:n) = this%step_scale*this%direction(off + 1:off + n)
            this%max_step = max(this%max_step, maxval(abs(this%step(1:n))))
            !$acc update device(this%step)
            one_norm = scf_soscf_ov_one_norm(nvir, nocc, reshape(this%step(1:n), [nvir, nocc]))
            call gdm_embed(this%n_mo, nocc, nvir, size(this%step), this%step, this%kappa)
            call scf_soscf_rotation_matrix(la, this%n_mo, one_norm, this%kappa, this%rotation, this%scratch_a, &
                                           this%scratch_b, this%scratch_c)
            call la%gemm_mat(this%c_acc(:, :, s), this%rotation, this%c_work)
         end if
         call cphf_copy_matrix(this%nao, this%n_mo, this%c_work, scf_result%cmo(:, 1:this%n_mo, s))
         call la%density(this%c_work, nocc, scf_config%occ, scf_result%dmat(:, :, s))
      end do
   end subroutine gdm_rotate

   subroutine scf_gdm_finish(this, scf_context, scf_config, scf_state, scf_result, scf_energies, la, comm)
      !! Leave the last ACCEPTED point in the result (every rank calls; a no-op
      !! when no GDM step ran): on the root, its density in `scf_state%dold`,
      !! its energies in `scf_energies` and its SEMICANONICAL orbitals and
      !! orbital energies in `scf_result%cmo`/`%eps`; then the density and
      !! energies to the other ranks (their last trial may have been a
      !! rejected one). Must run before the driver's device-to-host update.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(scf_energies_t), intent(inout) :: scf_energies
      type(trc_linalg_t), intent(inout) :: la
      type(comm_t), intent(in), optional :: comm
      real(dp) :: bundle(N_ENERGIES + 1)
      integer :: s

      if (scf_result%gdm_steps == 0) return
      if (scf_config%root .and. this%ready) then
         this%off_aufbau = 0
         do s = 1, this%nspin
            call la%density(this%c_acc(:, :, s), this%nocc(s), scf_config%occ, scf_state%dold(:, :, s))
            call gdm_semicanonical(this, scf_context, scf_result, la, s)
         end do
         scf_energies = this%energies_acc
         scf_result%gdm_max_rise = this%max_rise
         scf_result%gdm_max_rise_full = this%max_rise_full
         scf_result%gdm_max_step = this%max_step
         scf_result%gdm_off_aufbau = this%off_aufbau
         call gdm_say_off_aufbau(this)
      end if
      if (present(comm)) then
         if (comm%size() > 1) then
            !$acc update self(scf_state%dold)
            call scf_context_bcast_flat(comm, scf_state%dold, size(scf_state%dold))
            !$acc update device(scf_state%dold)
            bundle = [scf_energies%total, scf_energies%one, scf_energies%two, scf_energies%xc, scf_energies%nelec_grid, &
                      scf_energies%pcm, scf_energies%pcm_qtotal, real(scf_result%gdm_off_aufbau, dp)]
            call scf_context_bcast_flat(comm, bundle, N_ENERGIES + 1)
            scf_energies = scf_energies_t(bundle(1), bundle(2), bundle(3), bundle(4), bundle(5), bundle(6), bundle(7))
            scf_result%gdm_off_aufbau = nint(bundle(N_ENERGIES + 1))
         end if
      end if
   end subroutine scf_gdm_finish

   subroutine gdm_say_off_aufbau(this)
      !! The end of the run is not an aufbau state (root rank): the minimiser
      !! kept the occupation of its first diagonalisation and ended on an
      !! excited or saddle solution it cannot tell from the ground state.
      class(trc_scf_gdm_step_t), intent(in) :: this
      if (this%off_aufbau == 0) return
      call global_logger%warning("trc_scf: gdm converged to a state that is not the aufbau state ("// &
                                 to_string(this%off_aufbau)//" virtual orbital(s) below the highest occupied one); "// &
                                 "it is likely an excited or saddle solution -- GDM keeps the occupation of the "// &
                                 "diagonalisation it started from; start it later (a smaller gdm_start)")
   end subroutine gdm_say_off_aufbau

   subroutine gdm_semicanonical(this, scf_context, scf_result, la, s)
      !! Spin `s` of the accepted orbitals with their occupied-occupied and
      !! virtual-virtual Fock blocks diagonalised, into `scf_result%cmo`
      !! (columns 1:n_mo) and `%eps` (levels 1:n_mo; the lifted null states
      !! above are left as the handover diagonalisation set them). A spin with
      !! no occupied or no virtual orbital has one block: it is diagonalised
      !! whole. On an eigensolver failure the accepted orbitals go out as they
      !! are and `eps` is left stale, with a warning.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(scf_ctx_t), intent(in) :: scf_context
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: s
      type(trc_soscf_dims_t) :: dims
      integer :: info

      dims%sph = scf_context%sph
      dims%nao = this%nao
      dims%nao_c = scf_context%nao_c
      dims%n_mo = this%n_mo
      dims%nocc = this%nocc(s)
      dims%nvirt = this%nvir(s)
      call cphf_copy_matrix(this%nao, this%n_mo, this%c_acc(:, :, s), this%c_work)
      if (this%nocc(s) >= 1 .and. this%nvir(s) >= 1) then
         call scf_soscf_semicanonicalize(la, dims, this%fock_acc(:, :, s), this%c_work, this%fock_mo, &
                                         this%energies_mo, info)
      else
         call gdm_diagonalise_whole(this, la, s, info)
      end if
      call cphf_copy_matrix(this%nao, this%n_mo, this%c_work, scf_result%cmo(:, 1:this%n_mo, s))
      if (info == 0) then
         call scf_soscf_copy_eps(this%nao, this%n_mo, this%energies_mo, scf_result%eps(:, s))
         ! Per-block sorted: a virtual below the highest occupied is a non-aufbau state.
         ! Counted on the committed eps, by the same helper as the post-SCF refusal.
         this%off_aufbau = this%off_aufbau + scf_settings_count_off_aufbau(this%n_mo, this%nocc(s), &
                                                                           scf_result%eps(1:this%n_mo, s))
      else
         call global_logger%warning("trc_scf_gdm_step: the final semicanonicalisation did not converge "// &
                                    "(eigensolver info="//to_string(info)//"); the orbital energies are stale")
      end if
   end subroutine gdm_semicanonical

   subroutine gdm_diagonalise_whole(this, la, s, info)
      !! The single block of a spin with nothing to rotate: C <- C U with U
      !! the eigenvectors of C^T F C (`c_work` in and out), eigenvalues into
      !! `energies_mo`.
      class(trc_scf_gdm_step_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: s
      integer, intent(out) :: info
      call la%gemm_mat(this%fock_acc(:, :, s), this%c_work, this%work)
      call la%gemm_mat(this%c_work, this%work, this%fock_mo, transa="T")
      call la%syev(this%fock_mo, this%energies_mo, info)
      if (info /= 0) return
      call la%gemm_mat(this%c_work, this%fock_mo, this%work)
      call cphf_copy_matrix(this%nao, this%n_mo, this%work, this%c_work)
   end subroutine gdm_diagonalise_whole

   subroutine scf_gdm_report(this, scf_config, scf_result, mode)
      !! Count the iteration (every rank, so the counters agree) and, on the
      !! rank that talks, say what happened: one line at the handover, one
      !! after each GDM step's table line.
      class(trc_scf_gdm_step_t), intent(in) :: this
      type(scf_settings_t), intent(in) :: scf_config
      type(trc_scf_result_t), intent(inout) :: scf_result
      integer, intent(in) :: mode
      if (mode == GDM_OFF) return
      if (mode == GDM_HANDOVER) then
         if (scf_config%talk .and. this%declined) then
            call global_logger%info("   gdm: mom holds a non-aufbau occupation, a saddle point a minimiser would "// &
                                    "leave; gdm declines and the run continues as mom + DIIS")
         else if (scf_config%talk) then
            call global_logger%info("   gdm: this iteration diagonalises plainly; "// &
                                    "geometric direct minimisation steps from the next one")
         end if
         return
      end if
      scf_result%gdm_steps = scf_result%gdm_steps + 1
      if (this%rejected) scf_result%gdm_rejected = scf_result%gdm_rejected + 1
      if (scf_config%talk) call global_logger%info(trim(this%line))
   end subroutine scf_gdm_report

   subroutine gdm_pack(n_mo, nocc, nvirt, factor, fock_mo, ng, grad)
      !! grad((i-1)*nvirt + a) := factor*fock_mo(nocc + a, i): the
      !! occupied-virtual block of C^T F C scaled to dE/d kappa.
      integer, intent(in) :: n_mo, nocc, nvirt, ng
      real(dp), intent(in) :: factor
      real(dp), intent(in) :: fock_mo(n_mo, n_mo)
      real(dp), intent(inout) :: grad(ng)
      integer :: t, i, a
      do concurrent(t=1:nocc*nvirt) local(i, a)
         i = (t - 1)/nvirt + 1
         a = t - (i - 1)*nvirt
         grad(t) = factor*fock_mo(nocc + a, i)
      end do
   end subroutine gdm_pack

   subroutine gdm_embed(n_mo, nocc, nvirt, ng, step, kappa)
      !! kappa := the antisymmetric n_mo x n_mo matrix whose occupied-virtual
      !! block is `step` (index (i-1)*nvirt + a), zero elsewhere.
      integer, intent(in) :: n_mo, nocc, nvirt, ng
      real(dp), intent(in) :: step(ng)
      real(dp), intent(inout) :: kappa(n_mo, n_mo)
      integer :: i, j, t, a
      do concurrent(j=1:n_mo, i=1:n_mo)
         kappa(i, j) = 0.0_dp
      end do
      do concurrent(t=1:nocc*nvirt) local(i, a)
         i = (t - 1)/nvirt + 1
         a = t - (i - 1)*nvirt
         kappa(nocc + a, i) = step(t)
         kappa(i, nocc + a) = -step(t)
      end do
   end subroutine gdm_embed

end module trc_scf_gdm_step
