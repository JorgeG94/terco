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
   use pic_logger, only: global_logger
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   use pic_timer, only: timer_type
   use pic_profiler, only: profiler_start, profiler_stop
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_sph, only: sph_cart_to_sph
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t, ERROR_GENERIC
   use trc_dft_grid, only: dft_grid_t, build_dft_grid, build_dft_grid_block
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_linalg, only: trc_linalg_t
   use trc_pcm, only: trc_pcm_t, pcm_model_from_name
   use trc_scf_settings, only: trc_scf_options_t, scf_settings_t
   use trc_scf_df, only: df_jk_t
   implicit none
   private

   public :: scf_ctx_t, bcast_flat
   integer, parameter, public :: NSETUP = 6
   !! setup's steps: pair list, 1e integrals, spherical fold, ERI engine,
   !! linear-algebra handles, orthogonaliser -- `ts`'s six slots, printed
   !! with the rest of the timing table by the caller.
   integer, parameter, public :: NGRID = 4
   !! grid's own sub-steps: point generation (radial x Lebedev x pruning),
   !! the Becke/Stratmann partition, `trc_xc_grid_t%build`'s spatial
   !! batching, and `%to_device`'s upload -- `t_grid_sub`'s four slots.

   type :: scf_ctx_t
      !! Geometry-lifetime state: everything above is built once and read
      !! by every iteration, never written by one. `nao`/`nao_c`/`sph` are
      !! `b`'s own sizes, read here once so every other new module can take
      !! them as plain arguments instead of a basis handle.
      integer :: nao = 0, nao_c = 0
      logical :: sph = .false.
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(dft_grid_t) :: grid, grid_nlc
      type(trc_xc_grid_t) :: xg, xg_nlc
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
      logical :: grid_built = .false., nlc_built = .false., released = .false.
   contains
      procedure :: setup => scf_ctx_setup
      procedure :: release => scf_ctx_release
   end type scf_ctx_t

contains

   subroutine scf_ctx_setup(this, b, opts, cfg, err, comm, la, dfjk)
      class(scf_ctx_t), intent(inout), target :: this
      !! TARGET: `this%lp` may end up pointing at `this%la_own`, and that
      !! association must survive the return into `trc_scf_run`, which
      !! requires the actual argument (its own `ctx`) to carry TARGET too.
      type(trc_basis_t), intent(in) :: b
      type(trc_scf_options_t), intent(in) :: opts
      type(scf_settings_t), intent(in) :: cfg
      type(error_t), intent(out) :: err
      type(comm_t), intent(in), optional :: comm
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(df_jk_t), intent(inout), optional, target :: dfjk
      !! Already fitted by the caller (see `scf_ctx_t%dfj`'s own comment);
      !! required when `cfg%df` is true, checked here rather than left for
      !! `build_g` to find out at the first Fock build.
      type(timer_type) :: tm
      type(timer_type) :: tphase   !! the whole setup, and the PCM and grid phases

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
      if (cfg%df) then
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
      call this%pl%build(b, cfg%eri_thresh_run)
      call this%pl%to_device()
      this%ts(1) = tm%get_elapsed_time()
      call tm%start()
      ! The one-electron integrals are always Cartesian (trc_1e never
      ! changes); `smat`/`hcore` are folded up to the run's own space once,
      ! here, for the whole run -- S and H never change after this.
      block
         real(dp), allocatable :: s_c(:, :), h_c(:, :), v_c(:, :)
         allocate (s_c(this%nao_c, this%nao_c), h_c(this%nao_c, this%nao_c), v_c(this%nao_c, this%nao_c))
         call trc_1e(b, this%pl, s_c, h_c, v_c)
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
      ! Under --df, build_g never calls eri%fock_resident at all (the
      ! four-centre engine is not called for G, by design -- see
      ! trc_scf_df.F90's own header) -- so its setup (pair-quartet bins,
      ! Schwarz bounds, the device-resident work lists) is pure device
      ! memory this run has no use for. `eri%release` and build_g's own
      ! `ctx%eri` argument are both safe either way: release() guards
      ! every deallocation on `allocated`/`on_device`, and build_g reads
      ! nothing off `eri` before checking `df_on`.
      if (.not. cfg%df) then
         if (present(comm)) then
            call this%eri%build(b, cfg%eri_thresh_run, comm, batch_res=opts%batch_res)
         else
            call this%eri%build(b, cfg%eri_thresh_run, batch_res=opts%batch_res)
         end if
      end if
      this%e_nuc = nuclear_repulsion(b)
      this%ts(4) = tm%get_elapsed_time()
      call tm%start()
      call this%lp%init(this%nao)
      this%ts(5) = tm%get_elapsed_time()
      call tm%start()
      call sym_orthog(this%lp, this%nao, this%smat, this%x, opts%lindep_thresh, this%n_lindep)
      if (cfg%talk .and. this%n_lindep > 0) then
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
      if (cfg%pcm_on) then
         call tphase%start()
         call this%pcm%init(pcm_model_from_name(opts%pcm_model), opts%pcm_eps, err, la=this%lp)
         call this%pcm%build(b, opts%pcm_lebedev, opts%pcm_radii_scale, err)
         if (err%has_error()) then
            call this%release()
            return
         end if
         this%t_pcm_setup = tphase%get_elapsed_time()
         if (cfg%talk) then
            block
               character(len=1024) :: logline
               write (logline, "(a,a,a,f0.4,a,i0,a,i0,a,f0.3,a)") "   pcm: ", trim(opts%pcm_model), &
                  ", eps ", opts%pcm_eps, ", ", opts%pcm_lebedev, " points per sphere, ", &
                  this%pcm%surface%n_points, " kept; setup ", this%t_pcm_setup, " s"
               call global_logger%info(trim(logline))
            end block
         end if
      end if

      if (cfg%dft) then
         call tphase%start()
         call profiler_start("scf: grid", nvtx_only=.true.)
         ! The partition is split over the ranks; the summed cell weights are
         ! then rank 0's everywhere, as every reduced quantity is.
         if (present(comm)) then
            call build_dft_grid_block(b%at_r, nint(b%at_z), this%grid, err, comm%rank(), comm%size(), &
                                      level=opts%grid_level, &
                                      t_points=this%t_grid_sub(1), t_partition=this%t_grid_sub(2))
         else
            call build_dft_grid(b%at_r, nint(b%at_z), this%grid, err, level=opts%grid_level, &
                                t_points=this%t_grid_sub(1), t_partition=this%t_grid_sub(2))
         end if
         if (err%has_error()) then
            call profiler_stop("scf: grid")
            call this%release()
            return
         end if
         if (present(comm)) then
            if (comm%size() > 1) call sum_flat(comm, this%grid%weights, this%grid%n_points)
         end if
         call tm%start()
         call this%xg%build(this%grid%n_points, this%grid%coords, this%grid%weights, b, max_pts=opts%grid_max_pts)
         this%t_grid_sub(3) = tm%get_elapsed_time()
         call tm%start()
         call this%xg%to_device()
         this%t_grid_sub(4) = tm%get_elapsed_time()
         this%grid_built = .true.
         if (cfg%nlc) then
            call build_dft_grid(b%at_r, nint(b%at_z), this%grid_nlc, err, level=opts%nlc_grid_level)
            if (err%has_error()) then
               call profiler_stop("scf: grid")
               call this%release()
               return
            end if
            call this%xg_nlc%build(this%grid_nlc%n_points, this%grid_nlc%coords, this%grid_nlc%weights, b, &
                                   max_pts=opts%grid_max_pts)
            call this%xg_nlc%to_device()
            this%nlc_built = .true.
         end if
         call profiler_stop("scf: grid")
         this%t_grid = tphase%get_elapsed_time()
      end if
   end subroutine scf_ctx_setup

   subroutine scf_ctx_release(this)
      !! Idempotent: every early return in `setup` calls this once, and
      !! `trc_scf_run` calls it again at the end of a successful run, so a
      !! second call -- or a call on a context that never finished `setup`
      !! -- must be safe and must not release twice.
      class(scf_ctx_t), intent(inout) :: this
      if (this%released) return
      call this%pcm%release()
      if (this%owns_la) call this%lp%release()
      call this%eri%release()
      call this%pl%release()
      if (this%grid_built) then
         call this%xg%release()
         call this%grid%destroy()
         if (this%nlc_built) then
            call this%xg_nlc%release()
            call this%grid_nlc%destroy()
         end if
      end if
      this%released = .true.
   end subroutine scf_ctx_release

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
   subroutine sym_orthog(la, n, s, x, lindep_thresh, ndropped)
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: n
      real(dp), intent(in) :: s(n, n)
      real(dp), intent(out) :: x(n, n)
      real(dp), intent(in), optional :: lindep_thresh
      integer, intent(out), optional :: ndropped
      real(dp), allocatable :: u(:, :), w(:)
      real(dp) :: thr
      integer :: i, k, ndrop
      allocate (u(n, n), w(n))
      u = s
      thr = 0.0_dp
      if (present(lindep_thresh)) thr = lindep_thresh
      !$acc data copy(u) create(w) copyout(x)
      call la%syev(u, w)
      !$acc update self(w)
      ! w is ascending (LAPACK/cuSOLVER convention), so the near-null modes
      ! -- the ones canonical orthogonalisation drops -- are the leading
      ! entries, the same ones the old unconditional check read as w(1).
      ndrop = 0
      if (thr > 0.0_dp) then
         do k = 1, n
            if (w(k) < thr) ndrop = ndrop + 1
         end do
      else if (w(1) <= 0.0_dp) then
         error stop "trc_scf: the overlap is not positive definite"
      end if
      if (present(ndropped)) ndropped = ndrop
      ! Columns scaled by w^-1/4 on the device, then X = U U^T; the first
      ! `ndrop` columns -- the near-null modes -- are zeroed instead, so
      ! they carry no weight into X^T F X and never has to divide by a
      ! near-zero w.
      !$acc update device(w)
      do concurrent(i=1:n, k=1:n)
         if (k <= ndrop) then
            u(i, k) = 0.0_dp
         else
            u(i, k) = u(i, k)/sqrt(sqrt(w(k)))
         end if
      end do
      call la%gemm_mat(u, u, x, transb="T")
      !$acc end data
   end subroutine sym_orthog

   ! Sum a contiguous array over the ranks, then rank 0's copy everywhere.
   subroutine sum_flat(comm, g, n)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call allreduce(comm, g, op=MPI_SUM)
      call bcast(comm, g, n, 0)
   end subroutine sum_flat

   ! Broadcast rank 0's copy of a contiguous array of any rank. Public: the
   ! iteration loop (trc_scf_driver) broadcasts the density, orbitals and
   ! orbital energies the same way `setup` broadcasts the grid weights.
   subroutine bcast_flat(comm, g, n)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call bcast(comm, g, n, 0)
   end subroutine bcast_flat
end module trc_scf_context
