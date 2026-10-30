!
! Everything an SCF run builds ONCE PER GEOMETRY, before the iteration loop
! ever starts: the shell-pair list, the one-electron integrals folded into
! the run's own AO space, the four-centre engine, the linear-algebra
! handles and the canonical orthogonaliser, and -- when asked for -- the
! PCM cavity and the DFT quadrature. `setup` is the whole of what used to
! be `trc_scf_run`'s own preamble (its timing table's "setup" row and the
! PCM/grid blocks that followed it); `release` is every early-return's
! cleanup, consolidated into the one place that knows what it actually
! built, so a half-built context (an error building the VV10 grid, say)
! only releases what got that far.
!
module trc_scf_context
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   use pic_timer, only: timer_type
   use pic_profiler, only: profiler_start, profiler_stop
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_1e
   use trc_sph, only: sph_cart_to_sph
   use trc_eri, only: trc_eri_t, eri_lmax_refusal
   use trc_error, only: error_t, ERROR_GENERIC
   use trc_dft_grid, only: dft_grid_t, dft_grid_build, dft_grid_build_block
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_linalg, only: trc_linalg_t
   use trc_pcm, only: trc_pcm_t, pcm_model_from_name
   use trc_scf_settings, only: trc_scf_options_t, scf_settings_t
   use trc_scf_df, only: df_jk_t
   implicit none
   private

   public :: scf_ctx_t, scf_context_bcast_flat, trc_grid_cache_t, trc_pair_eri_cache_t, scf_context_project_density
   public :: scf_context_orthogonalize   !! the overlap pseudo-inverse builder, also trc_minao's own cross-basis projection
   integer, parameter, public :: NSETUP = 6
   !! setup's steps: pair list, 1e integrals, spherical fold, ERI engine,
   !! linear-algebra handles, orthogonaliser -- `ts`'s six slots, printed
   !! with the rest of the timing table by the caller.
   integer, parameter, public :: NGRID = 4
   !! grid's own sub-steps: point generation (radial x Lebedev x pruning),
   !! the Becke/Stratmann partition, `trc_xc_grid_t%build`'s spatial
   !! batching, and `%to_device`'s upload -- `t_grid_sub`'s four slots.

   type :: trc_grid_cache_t
      !! The semilocal and (for VV10) nonlocal DFT grids a run built, kept
      !! apart from `scf_ctx_t` so a caller that wants them ALIVE past
      !! `setup`'s own run (trc_capi.F90's context handle, so the gradient
      !! that follows can reuse exactly what the SCF converged with rather
      !! than rebuild it -- see trc_scf_driver's `trc_scf_run`'s own
      !! `grid_cache` argument) can own one across more than one call,
      !! the same pattern `la`/`dfjk` already use. Owned by `scf_ctx_t`'s
      !! own `gc_own` when no caller hands one in.
      type(dft_grid_t) :: grid, grid_nlc
      type(trc_xc_grid_t) :: xg, xg_nlc
      logical :: grid_built = .false., nlc_built = .false.
   contains
      procedure :: release => grid_cache_release
   end type trc_grid_cache_t

   type :: trc_pair_eri_cache_t
      !! The shell-pair list and (when built) the four-centre engine from
      !! ONE `scf_ctx_t%setup` call, kept apart from `scf_ctx_t` so a caller
      !! that wants them ALIVE past that call's own `release` (trc_capi.F90's
      !! context handle) can hand them straight to a gradient run at the
      !! SAME geometry and threshold right after -- the same pattern
      !! `trc_grid_cache_t` already uses for the DFT grid, and for the same
      !! reason: a gradient that follows the SCF immediately has nothing
      !! left to learn from rebuilding a pair list and a Schwarz screen the
      !! SCF already paid for at this geometry. `thresh_built` is what a
      !! reuser must match EXACTLY (bit for bit -- both sides read it from
      !! the same `scf_settings_eri_thresh_env` call) before trusting `pl`/`eri` at all;
      !! `has_eri` is false under `--df`, where `scf_context_setup` never calls
      !! `eri%build` (trc_scf_df's own J/K, not the four-centre engine).
      !!
      !! FINGERPRINT. `thresh_built` is the only gate a C API caller ever
      !! needs (a geometry or basis change already clears `built` through
      !! `invalidate_runs`, before a reuser could reach this far), but a
      !! Fortran caller that keeps one `trc_pair_eri_cache_t` across more
      !! than one `trc_scf_run`/`trc_scf_gradient` call, at a DIFFERENT
      !! geometry, with no such invalidation to rely on, would otherwise
      !! get a gradient from stale pair data, or out-of-bounds indexing if
      !! `nshell` changed. `nshell_built`/`natm_built`/`at_r_built` are
      !! that caller's own safety net: set once, at build time
      !! (`scf_context_setup`), and checked by `matches` before ANY reuse,
      !! Fortran or C API alike. `batch_res_built` is the quartet-bin
      !! batching the engine was actually built with (`trc_scf_settings`'s
      !! own `opts%batch_res`) -- `trc_scf_gradient`'s OWN fresh build
      !! always uses the default (1), so a cache built at any other value
      !! would give a reused-vs-fresh gradient that differs by more than
      !! round-off; `matches` is the pair list's own gate (thresh and
      !! geometry only), so the ERI engine's reuse checks `batch_res_built`
      !! on top, in trc_scf_grad.F90.
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      real(dp) :: thresh_built = -1.0_dp
      integer :: nshell_built = -1, natm_built = -1, batch_res_built = -1
      real(dp), allocatable :: at_r_built(:, :)
      logical :: built = .false., has_eri = .false.
   contains
      procedure :: release => pair_eri_cache_release
      procedure :: matches => pair_eri_cache_matches
   end type trc_pair_eri_cache_t

   type :: scf_ctx_t
      !! Geometry-lifetime state: everything above is built once and read
      !! by every iteration, never written by one. `nao`/`nao_c`/`sph` are
      !! `b`'s own sizes, read here once so every other new module can take
      !! them as plain arguments instead of a basis handle.
      integer :: nao = 0, nao_c = 0
      logical :: sph = .false.
      type(trc_pair_eri_cache_t) :: pec_own
      !! TARGET for `pec` to point to comes from `this` itself, the same
      !! reason `gc_own`/`la_own` are -- see `lp`'s own comment, right below.
      type(trc_pair_eri_cache_t), pointer :: pec => null()
      logical :: owns_pec = .false.
      type(trc_grid_cache_t) :: gc_own
      !! TARGET for `gc` to point to comes from `this` itself, the same
      !! reason `la_own` is -- see `lp`'s own comment, right below.
      type(trc_grid_cache_t), pointer :: gc => null()
      logical :: owns_grid = .false.
      type(trc_pcm_t) :: pcm
      type(trc_linalg_t) :: la_own
      !! TARGET for `lp` to point to comes from `this` itself (see `setup`'s
      !! dummy and the caller's own `ctx` declaration), not from a component
      !! attribute: Fortran does not allow TARGET on a type's component.
      type(trc_linalg_t), pointer :: lp => null()
      logical :: owns_la = .false.
      !! The density-fitted J/K object (`--df`), borrowed from the caller
      !! (trc_capi.F90's context handle owns one per geometry, so it
      !! outlives a single `trc_scf_run`; see trc_scf_df.F90's own header)
      !! and already fitted -- `setup` never calls `dfj%setup` itself, only
      !! points at what it was handed. Null when `cfg%df` is false.
      type(df_jk_t), pointer :: dfj => null()
      real(dp), allocatable :: smat(:, :), hcore(:, :), x(:, :)
      real(dp) :: e_nuc = 0.0_dp
      integer :: n_lindep = 0
      real(dp) :: ts(NSETUP) = 0.0_dp
      real(dp) :: t_grid_sub(NGRID) = 0.0_dp
      real(dp) :: t_setup = 0.0_dp, t_grid = 0.0_dp, t_pcm_setup = 0.0_dp
      logical :: released = .false.
   contains
      procedure :: setup => scf_context_setup
      procedure :: release => scf_context_release
   end type scf_ctx_t

contains

   subroutine scf_context_setup(this, b, opts, config, err, comm, la, dfjk, grid_cache, pair_eri_cache)
      class(scf_ctx_t), intent(inout), target :: this
      !! TARGET: `this%lp` may end up pointing at `this%la_own`, and that
      !! association must survive the return into `trc_scf_run`, which
      !! requires the actual argument (its own `ctx`) to carry TARGET too.
      type(trc_basis_t), intent(in) :: b
      type(trc_scf_options_t), intent(in) :: opts
      type(scf_settings_t), intent(in) :: config
      type(error_t), intent(out) :: err
      type(comm_t), intent(in), optional :: comm
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(df_jk_t), intent(inout), optional, target :: dfjk
      !! Already fitted by the caller (see `scf_ctx_t%dfj`'s own comment);
      !! required when `config%df` is true, checked here rather than left for
      !! `scf_fock_build_g` to find out at the first Fock build.
      type(trc_grid_cache_t), intent(inout), optional, target :: grid_cache
      !! A caller-owned grid cache (trc_capi.F90's context handle) that
      !! this run builds the DFT grid DIRECTLY INTO -- not into a local
      !! object that is then copied or moved -- so the grid (and, once
      !! `%to_device` runs below, its device mapping) is built once, at
      !! its final, long-lived address, and stays valid for a gradient
      !! that runs after this call returns (the same reason `take_basis`
      !! in trc_capi.F90 copies a basis up BEFORE uploading it, never
      !! after: the device present-table is keyed on host addresses).
      !! Absent, `this%gc_own` plays the same role but does not outlive
      !! this call -- unchanged from before this argument existed.
      type(trc_pair_eri_cache_t), intent(inout), optional, target :: pair_eri_cache
      !! A caller-owned pair-list/ERI cache (trc_capi.F90's context handle)
      !! that this run builds the pair list -- and, outside `--df`, the
      !! four-centre engine -- DIRECTLY INTO, same reasoning as
      !! `grid_cache` right above: built once, at its final address, so a
      !! gradient immediately after can reuse it without a second build at
      !! the same geometry and threshold (`trc_pair_eri_cache_t`'s own
      !! header). Absent, `this%pec_own` plays the same role but does not
      !! outlive this call.
      type(timer_type) :: tm
      type(timer_type) :: tphase   !! the whole setup, and the PCM and grid phases
      character(len=:), allocatable :: refusal

      ! Before anything is allocated or uploaded: a shell above this build's
      ! four-centre kernels has no kernel, and its integrals would come out
      ! as NaN. A --df run checks its own ceiling in df_jk_t%setup.
      if (.not. config%df) then
         refusal = eri_lmax_refusal(b)
         if (len(refusal) > 0) then
            call err%set(ERROR_GENERIC, "trc_scf: "//refusal)
            return
         end if
      end if

      this%nao_c = b%nao
      this%sph = b%spherical
      this%nao = merge(b%nao_sph, this%nao_c, this%sph)
      allocate (this%smat(this%nao, this%nao), this%hcore(this%nao, this%nao), this%x(this%nao, this%nao))
      this%owns_la = .not. present(la)
      if (present(la)) then
         this%lp => la
      else
         this%lp => this%la_own
      end if
      this%owns_grid = .not. present(grid_cache)
      if (present(grid_cache)) then
         this%gc => grid_cache
      else
         this%gc => this%gc_own
      end if
      this%owns_pec = .not. present(pair_eri_cache)
      if (present(pair_eri_cache)) then
         this%pec => pair_eri_cache
      else
         this%pec => this%pec_own
      end if
      if (config%df) then
         if (.not. present(dfjk)) then
            call err%set(ERROR_GENERIC, "trc_scf: --df needs a density-fitted object already set up "// &
                         "(trc_set_df/trc_set_jkfit_json)")
            return
         end if
         if (.not. dfjk%is_ready()) then
            call err%set(ERROR_GENERIC, "trc_scf: --df's density-fitted object was never fitted "// &
                         "(trc_set_jkfit_json before trc_run_scf)")
            return
         end if
         this%dfj => dfjk
      end if

      ! --- once per geometry, on the host ------------------------------------
      call tphase%start()
      call tm%start()
      call this%pec%pl%build(b, config%eri_thresh_run)
      call this%pec%pl%to_device()
      this%pec%thresh_built = config%eri_thresh_run
      this%pec%built = .true.
      this%pec%has_eri = .false.
      ! The fingerprint a Fortran caller's reuse is checked against --
      ! see `trc_pair_eri_cache_t`'s own header, FINGERPRINT.
      this%pec%nshell_built = b%nshell
      this%pec%natm_built = b%natm
      if (allocated(this%pec%at_r_built)) deallocate (this%pec%at_r_built)
      this%pec%at_r_built = b%at_r
      this%ts(1) = tm%get_elapsed_time()
      call tm%start()
      ! The one-electron integrals are always Cartesian (trc_1e never
      ! changes); `smat`/`hcore` are folded up to the run's own space once,
      ! here, for the whole run -- S and H never change after this.
      block
         real(dp), allocatable :: s_c(:, :), h_c(:, :), v_c(:, :)
         allocate (s_c(this%nao_c, this%nao_c), h_c(this%nao_c, this%nao_c), v_c(this%nao_c, this%nao_c))
         call trc_1e(b, this%pec%pl, s_c, h_c, v_c)
         h_c = h_c + v_c
         this%ts(2) = tm%get_elapsed_time()
         call tm%start()
         if (this%sph) then
            !$acc enter data copyin(s_c, h_c) create(this%smat, this%hcore)
            call sph_cart_to_sph(b, s_c, this%smat)
            call sph_cart_to_sph(b, h_c, this%hcore)
            !$acc update self(this%smat, this%hcore)
            !$acc exit data delete(s_c, h_c, this%smat, this%hcore)
         else
            this%smat = s_c
            this%hcore = h_c
         end if
      end block
      this%ts(3) = tm%get_elapsed_time()
      call tm%start()
      ! Under --df, scf_fock_build_g never calls eri%fock_resident at all (the
      ! four-centre engine is not called for G, by design -- see
      ! trc_scf_df.F90's own header) -- so its setup (pair-quartet bins,
      ! Schwarz bounds, the device-resident work lists) is pure device
      ! memory this run has no use for. `eri%release` and scf_fock_build_g's own
      ! `ctx%eri` argument are both safe either way: release() guards
      ! every deallocation on `allocated`/`on_device`, and scf_fock_build_g reads
      ! nothing off `eri` before checking `df_on`.
      if (.not. config%df) then
         this%pec%has_eri = .true.
         this%pec%batch_res_built = opts%batch_res
         if (present(comm)) then
            call this%pec%eri%build(b, config%eri_thresh_run, comm, batch_res=opts%batch_res)
            ! Not when `comm` is a real split (nranks > 1): each rank then
            ! holds only ITS share of the quartets, which is not what a
            ! reuser (a single-rank gradient right after, trc_capi.F90's
            ! trc_run_gradient) expects of "the SCF's own engine" --
            ! trc_run_hessian already refuses more than one rank, so this
            ! is a backstop, not the gate. `.and.` is not guaranteed to
            ! short-circuit, hence the nested `if` rather than one
            ! expression that would call `comm%size()` even when absent.
            ! Marked unusable for a reuser, NOT released: this SCF's own
            ! Fock builds (scf_fock_build_g -> fock_resident) read `pec%eri` every
            ! iteration, and on deallocated bins they segfault on the host,
            ! which nvfortran's handler turns into a silent spin.
            if (comm%size() > 1) this%pec%has_eri = .false.
         else
            call this%pec%eri%build(b, config%eri_thresh_run, batch_res=opts%batch_res)
         end if
      end if
      this%e_nuc = nuclear_repulsion(b)
      this%ts(4) = tm%get_elapsed_time()
      call tm%start()
      call this%lp%init(this%nao)
      this%ts(5) = tm%get_elapsed_time()
      call tm%start()
      call scf_context_orthogonalize(this%lp, this%nao, this%smat, this%x, err, opts%lindep_thresh, this%n_lindep)
      if (err%has_error()) then
         call this%release()
         return
      end if
      if (config%talk .and. this%n_lindep > 0) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a,es9.2,a)") &
               "   canonical orthogonalisation: ", this%n_lindep, &
               " function(s) dropped, S eigenvalue below ", opts%lindep_thresh, ""
            call global_logger%info(trim(logline))
         end block
      end if
      this%ts(6) = tm%get_elapsed_time()
      this%t_setup = tphase%get_elapsed_time()

      ! --- the continuum: cavity, K and the nuclear potential, once -------
      if (config%pcm_on) then
         call tphase%start()
         call this%pcm%init(pcm_model_from_name(opts%pcm_model), opts%pcm_eps, err, la=this%lp)
         if (allocated(opts%pcm_radii_override)) then
            call this%pcm%build(b, opts%pcm_lebedev, opts%pcm_radii_scale, err, &
                                radii_override=opts%pcm_radii_override)
         else
            call this%pcm%build(b, opts%pcm_lebedev, opts%pcm_radii_scale, err)
         end if
         if (err%has_error()) then
            call this%release()
            return
         end if
         this%t_pcm_setup = tphase%get_elapsed_time()
         if (config%talk) then
            block
               character(len=1024) :: logline
               write (logline, "(a,a,a,f0.4,a,i0,a,i0,a,f0.3,a)") "   pcm: ", trim(opts%pcm_model), &
                  ", eps ", opts%pcm_eps, ", ", opts%pcm_lebedev, " points per sphere, ", &
                  this%pcm%surface%n_points, " kept; setup ", this%t_pcm_setup, " s"
               call global_logger%info(trim(logline))
            end block
         end if
      end if

      if (config%dft) then
         call tphase%start()
         call profiler_start("scf: grid", nvtx_only=.true.)
         ! The partition is split over the ranks; the summed cell weights are
         ! then rank 0's everywhere, as every reduced quantity is.
         if (present(comm)) then
            call dft_grid_build_block(b%at_r, b%at_elem, this%gc%grid, err, comm%rank(), comm%size(), &
                                      level=opts%grid_level, &
                                      t_points=this%t_grid_sub(1), t_partition=this%t_grid_sub(2))
         else
            call dft_grid_build(b%at_r, b%at_elem, this%gc%grid, err, level=opts%grid_level, &
                                t_points=this%t_grid_sub(1), t_partition=this%t_grid_sub(2))
         end if
         if (err%has_error()) then
            call profiler_stop("scf: grid")
            call this%release()
            return
         end if
         if (present(comm)) then
            if (comm%size() > 1) call sum_flat(comm, this%gc%grid%weights, this%gc%grid%n_points)
         end if
         call tm%start()
         call this%gc%xg%build(this%gc%grid%n_points, this%gc%grid%coords, this%gc%grid%weights, b, &
                               max_pts=opts%grid_max_pts, atom=this%gc%grid%atom, qw=this%gc%grid%quad_weights)
         this%t_grid_sub(3) = tm%get_elapsed_time()
         call tm%start()
         call this%gc%xg%to_device()
         this%t_grid_sub(4) = tm%get_elapsed_time()
         this%gc%grid_built = .true.
         if (config%nlc) then
            call dft_grid_build(b%at_r, b%at_elem, this%gc%grid_nlc, err, level=opts%nlc_grid_level)
            if (err%has_error()) then
               call profiler_stop("scf: grid")
               call this%release()
               return
            end if
            call this%gc%xg_nlc%build(this%gc%grid_nlc%n_points, this%gc%grid_nlc%coords, this%gc%grid_nlc%weights, &
                                      b, max_pts=opts%grid_max_pts, &
                                      atom=this%gc%grid_nlc%atom, qw=this%gc%grid_nlc%quad_weights)
            call this%gc%xg_nlc%to_device()
            this%gc%nlc_built = .true.
         end if
         call profiler_stop("scf: grid")
         this%t_grid = tphase%get_elapsed_time()
      end if
   end subroutine scf_context_setup

   subroutine scf_context_release(this)
      !! Idempotent: every early return in `setup` calls this once, and
      !! `trc_scf_run` calls it again at the end of a successful run, so a
      !! second call -- or a call on a context that never finished `setup`
      !! -- must be safe and must not release twice.
      class(scf_ctx_t), intent(inout) :: this
      if (this%released) return
      call this%pcm%release()
      if (this%owns_la) call this%lp%release()
      ! A caller-owned cache (pair_eri_cache present at `setup`) outlives
      ! this context on purpose, same reasoning as `gc` right below -- only
      ! an OWNED cache (`this%pec_own`, the default) is released here.
      if (this%owns_pec) call this%pec%release()
      ! A caller-owned cache (grid_cache present at `setup`) outlives this
      ! context on purpose -- see `scf_context_setup`'s own header on
      ! `grid_cache` -- so only an OWNED cache (`this%gc_own`, the default)
      ! is released here.
      if (this%owns_grid) call this%gc%release()
      this%released = .true.
   end subroutine scf_context_release

   subroutine grid_cache_release(this)
      !! Idempotent, same reason as `scf_context_release`'s own: a cache that
      !! was never built (no DFT run reached it, or it has already been
      !! released) costs nothing to release again.
      class(trc_grid_cache_t), intent(inout) :: this
      if (this%grid_built) then
         call this%xg%release()
         call this%grid%destroy()
         if (this%nlc_built) then
            call this%xg_nlc%release()
            call this%grid_nlc%destroy()
         end if
         this%grid_built = .false.
         this%nlc_built = .false.
      end if
   end subroutine grid_cache_release

   subroutine pair_eri_cache_release(this)
      !! Idempotent, same reason as `scf_context_release`'s own: `pl%release`/
      !! `eri%release` already guard every deallocation on `allocated`/
      !! `on_device` (that module's own invariant), so releasing a cache
      !! that was never built, or has already been released, costs nothing.
      class(trc_pair_eri_cache_t), intent(inout) :: this
      call this%eri%release()
      call this%pl%release()
      this%built = .false.
      this%has_eri = .false.
      this%thresh_built = -1.0_dp
      this%nshell_built = -1
      this%natm_built = -1
      this%batch_res_built = -1
      if (allocated(this%at_r_built)) deallocate (this%at_r_built)
   end subroutine pair_eri_cache_release

   logical function pair_eri_cache_matches(this, b, thresh) result(ok)
      !! The reuse gate every caller (Fortran or, through
      !! `trc_scf_gradient`'s `pair_eri_cache`, the C API) must pass before
      !! trusting `pl`/`eri`: the build threshold EXACTLY (bit for bit, see
      !! `thresh_built`'s own comment above) and the geometry/basis
      !! fingerprint (`nshell`, `natm`, the atom centres) unchanged since
      !! `scf_context_setup` built it. A cache that was never built, or has
      !! since been released, fails on `built` alone.
      class(trc_pair_eri_cache_t), intent(in) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      ok = this%built .and. this%thresh_built == thresh .and. &
           this%nshell_built == b%nshell .and. this%natm_built == b%natm
      if (ok) ok = allocated(this%at_r_built)
      if (ok) ok = size(this%at_r_built, 2) == b%natm
      if (ok) ok = all(this%at_r_built == b%at_r)
   end function pair_eri_cache_matches

   real(dp) function nuclear_repulsion(b) result(e)
      type(trc_basis_t), intent(in) :: b
      integer :: a, c
      e = 0.0_dp
      do a = 1, b%natm
         do c = a + 1, b%natm
            e = e + b%at_z(a)*b%at_z(c)/norm2(b%at_r(:, a) - b%at_r(:, c))
         end do
      end do
   end function nuclear_repulsion

   ! X = S^(-1/2), symmetric orthogonalisation, or, once `lindep_thresh` is positive
   ! and some S eigenvalue falls below it, CANONICAL orthogonalisation: those
   ! columns of U are zeroed rather than divided by sqrt(sqrt(w)), which is what
   ! makes symmetric orthogonalisation diverge on a near-null mode. A zeroed column
   ! decouples that mode from X^T F X. Once per geometry, on the host side.
   !
   ! `err` instead of `error stop` on either failure (the overlap's own
   ! diagonalisation not converging, or -- when the caller left
   ! lindep_thresh at 0, so canonical orthogonalisation never engages --
   ! a non-positive eigenvalue): a malformed overlap must come back as a
   ! clean SCF failure, the same way `scf_context_setup`'s other refusals
   ! (the LMAX ceiling, a missing --df fit) already do, not take the
   ! whole process down. Structured so every path still reaches the
   ! matching `!$acc end data` -- `return`ing out of an OpenACC `data`
   ! region partway through is not something the directive allows.
   subroutine scf_context_orthogonalize(la, n, s, x, err, lindep_thresh, ndropped)
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: n
      real(dp), intent(in) :: s(n, n)
      real(dp), intent(out) :: x(n, n)
      type(error_t), intent(inout) :: err
      real(dp), intent(in), optional :: lindep_thresh
      integer, intent(out), optional :: ndropped
      real(dp), allocatable :: u(:, :), w(:)
      real(dp) :: threshold
      integer :: i, k, ndrop, info
      logical :: bad_info, bad_pd
      allocate (u(n, n), w(n))
      u = s
      x = 0.0_dp
      threshold = 0.0_dp
      if (present(lindep_thresh)) threshold = lindep_thresh
      info = 0
      bad_info = .false.
      bad_pd = .false.
      ndrop = 0
      !$acc data copy(u) create(w) copyout(x)
      call la%syev(u, w, info)
      if (info /= 0) then
         bad_info = .true.
      else
         !$acc update self(w)
         ! w is ascending (LAPACK/cuSOLVER convention), so the near-null
         ! modes -- the ones canonical orthogonalisation drops -- are the
         ! leading entries, the same ones the old unconditional check read
         ! as w(1).
         if (threshold > 0.0_dp) then
            do k = 1, n
               if (w(k) < threshold) ndrop = ndrop + 1
            end do
         else if (w(1) <= 0.0_dp) then
            bad_pd = .true.
         end if
         if (.not. bad_pd) then
            ! Columns scaled by w^-1/4 on the device, then X = U U^T; the
            ! first `ndrop` columns -- the near-null modes -- are zeroed
            ! instead, so they carry no weight into X^T F X and never has
            ! to divide by a near-zero w.
            !$acc update device(w)
            do concurrent(i=1:n, k=1:n)
               if (k <= ndrop) then
                  u(i, k) = 0.0_dp
               else
                  u(i, k) = u(i, k)/sqrt(sqrt(w(k)))
               end if
            end do
            call la%gemm_mat(u, u, x, transb="T")
         end if
      end if
      !$acc end data
      if (bad_info) then
         call err%set(ERROR_GENERIC, "trc_scf: the overlap diagonalisation did not converge "// &
                      "(pic_syev info="//to_string(info)//")")
         return
      end if
      if (bad_pd) then
         call err%set(ERROR_GENERIC, "trc_scf: the overlap is not positive definite")
         return
      end if
      if (present(ndropped)) ndropped = ndrop
   end subroutine scf_context_orthogonalize

   ! Sum a contiguous array over the ranks, then rank 0's copy everywhere.
   subroutine sum_flat(comm, g, n)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call allreduce(comm, g, op=MPI_SUM)
      call bcast(comm, g, n, 0)
   end subroutine sum_flat

   subroutine scf_context_project_density(lp, nao, nspin, nocc, cmo_old, s_new, dmat)
      !! A warm-start density at a NEW geometry (same molecule/basis, so
      !! `nao` is unchanged) from a PREVIOUS geometry's converged occupied
      !! orbitals, re-orthonormalised in the new overlap `s_new`: for each
      !! spin's occupied block C (nao x nocc), M = C^T S' C is its overlap
      !! under S'; `scf_context_orthogonalize`'s own symmetric/canonical orthogonalisation
      !! works on ANY symmetric positive matrix handed to it, not only the
      !! AO overlap, so calling it on M gives X = M^(-1/2); C' = C X is
      !! once again orthonormal under S' (C'^T S' C' = X^T M X = X M X = I
      !! by construction, since X = M^(-1/2) commutes with M). D' =
      !! occ_s C' C'^T. At S' = S this is the identity transform (M = I
      !! already, X = I); trc_capi.F90's TRC_GUESS_PROJECTED is the public
      !! way to this. test/bench_warmstart.F90 measured it against the
      !! plain density warm start at water/def2-SVP, steps 0.005-1.0 Bohr:
      !! both needed the same iteration count at every step (HF 10/13/15,
      !! B3LYP 8/10/12).
      !!
      !! A near-singular M (an occupied orbital that left the space, e.g a
      !! very large step) is not refused here -- `scf_context_orthogonalize`'s own
      !! `lindep_thresh=0` default only error-stops on a NON-positive-
      !! definite M, which this call never raises
      !! (`lindep_thresh` absent): nocc columns of a full-rank C projected
      !! through a nonsingular S' make M positive definite by
      !! construction (it is a Gram matrix of independent vectors), so a
      !! near-null mode is never actually reached by this routine's own
      !! use. Host arrays in and out, the same residency the guess block
      !! calling this already assumes; `scf_context_orthogonalize`'s own `!$acc data`
      !! moves the (small, nocc x nocc) eigenproblem to the device and
      !! back.
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nspin
      integer, intent(in) :: nocc(2)
      real(dp), intent(in) :: cmo_old(nao, nao, nspin)
      real(dp), intent(in) :: s_new(nao, nao)
      real(dp), intent(out) :: dmat(nao, nao, nspin)
      real(dp) :: occ_s
      integer :: s, n
      real(dp), allocatable :: cocc(:, :), m(:, :), x(:, :), cproj(:, :)
      type(error_t) :: err

      occ_s = merge(2.0_dp, 1.0_dp, nspin == 1)
      dmat = 0.0_dp
      do s = 1, nspin
         n = nocc(s)
         if (n <= 0) cycle
         allocate (cocc(nao, n), m(n, n), x(n, n), cproj(nao, n))
         cocc = cmo_old(:, 1:n, s)
         m = matmul(transpose(cocc), matmul(s_new, cocc))
         m = 0.5_dp*(m + transpose(m))   ! host round-off symmetrisation
         call scf_context_orthogonalize(lp, n, m, x, err)
         if (err%has_error()) then
            ! Argued unreachable by this routine's own header (M is a Gram
            ! matrix of independent vectors, positive definite by
            ! construction) -- defensive only: leave this spin's slice at
            ! the zero `dmat` already set above rather than use a
            ! half-written `x`, and let the caller fall back to its own
            ! ordinary guess rather than crash.
            call global_logger%warning("trc_scf: the projected-guess overlap failed to diagonalise "// &
                                       "for spin "//to_string(s)//" ("//err%get_message()// &
                                       "); falling back to a zero density for that spin")
            deallocate (cocc, m, x, cproj)
            cycle
         end if
         cproj = matmul(cocc, x)
         dmat(:, :, s) = occ_s*matmul(cproj, transpose(cproj))
         deallocate (cocc, m, x, cproj)
      end do
   end subroutine scf_context_project_density

   ! Broadcast rank 0's copy of a contiguous array of any rank. Public: the
   ! iteration loop (trc_scf_driver) broadcasts the density, orbitals and
   ! orbital energies the same way `setup` broadcasts the grid weights.
   subroutine scf_context_bcast_flat(comm, g, n)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call bcast(comm, g, n, 0)
   end subroutine scf_context_bcast_flat
end module trc_scf_context
