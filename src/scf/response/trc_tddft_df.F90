module trc_tddft_df
   !!
   !! The two-electron part of the TDA and RPA excitation operators over a
   !! density-fitted (--df) RHF/RKS reference: the products
   !!
   !!     J[x]_ia  = sum_jb (ia|jb) x_jb
   !!     K1[x]_ia = sum_jb (ij|ab) x_jb
   !!     K2[x]_ia = sum_jb (ib|ja) x_jb
   !!
   !! with every four-centre integral replaced by the fit it was converged
   !! with, (pq|rs) = sum_Q B^Q_pq B^Q_rs (B = T L^-T, `df_jk_t%bpk`), so the
   !! excitation energies are the response of the FITTED energy.
   !! `trc_tddft` combines them (its header has the derivation, term by term
   !! against the four-centre matvec): TDA 2 J - c_x K1, RPA (A+B) 4 J - c_x
   !! (K1 + K2) and (A-B) -c_x (K1 - K2). The XC kernel is not fitted and not here.
   !!
   !! THE PACKED B. The rebuilds read it, so `setup` takes it over from the
   !! `df_jk_t` (a move, no copy, the device mapping untouched) and `restore`
   !! gives it back: no application needs the fit passed down every call, and
   !! nothing else touches the fit in between.
   !!
   !! THE TENSORS, built once per solve by one pass over the packed B. For a
   !! chunk of np auxiliary functions, with C_o, C_v the occupied and virtual
   !! orbitals (the fit's Cartesian AO space; a spherical auxiliary set needs
   !! nothing here, `df_jk_t` has folded it into B's auxiliary index):
   !!
   !!     boo(i, j, Q)  = B^Q_ij        (nocc, nocc, np)   only with exchange
   !!     bov(i, a, Q)  = B^Q_ia        (nocc, nvir, np)
   !!     bvv(a, Q, b)  = B^Q_ab        (nvir, np, nvir)   only with exchange
   !!
   !! (a pure functional needs bov alone: nocc nvir naux doubles). They come
   !! from ONE half-transform GEMM of the unpacked block, C^T B, and three
   !! batched GEMMs against C_v and C_o -- the construction `df_rotation_t`
   !! uses for the CPHF Hessian -- and no AO quantity exists after that.
   !!
   !! ONE APPLICATION on g <= vmax trial vectors x (nocc, nvir, g), per chunk,
   !! all GEMMs on the device through `trc_linalg_t`:
   !!
   !!   J   t = bov^T x (nvec at once), J += bov t               2 GEMMs
   !!   K1  W(i, k, b, Q) = sum_j boo(j, i, Q) x(j, b, k)        1 batched GEMM
   !!       (x reordered once to (j, k, b), shared by every Q), then
   !!       K1(i, k, a) += W((i, k), (b, Q)) bvv((b, Q), a)       1 GEMM
   !!   K2  per vector: R(b, Q, a) = sum_j x(j, b) bov(j, a, Q)  1 batched GEMM
   !!       (x^T shared by every Q), K2(:, :, k) += bov((i), (b, Q)) R   1 GEMM
   !!
   !! so the contraction over the auxiliary index is always INSIDE one GEMM
   !! (the (b, Q) pair is one contiguous index of both operands) and there is
   !! no reduction kernel. Flops per chunk and vector 2 np nocc nvir (2 nocc +
   !! 2 nvir) for K1, the same for K2, and 4 np nocc nvir for J. All the
   !! vectors of a call share one pass over the blocks.
   !!
   !! MEMORY AND BLOCKING. bvv is nvir^2 naux doubles and does not always fit
   !! beside the packed B the SCF keeps. `setup` sizes, from the device's free
   !! memory (a half gigabyte kept back, and the transition densities the
   !! caller keeps), the chunk `qblk` (auxiliary functions per pass) and how
   !! many chunks keep their blocks resident for the whole solve (`nres`
   !! functions). The remaining chunks are REBUILT from the packed B inside
   !! every application, into one reused slot -- device work (an unpack, the
   !! half-transform and three batched GEMMs per chunk), never a host
   !! fallback -- so a solve whose blocks do not fit is slower, not refused.
   !! A system for which not even one auxiliary function fits stops with the
   !! sizes (`err`). TERCO_TDDFT_DF_BLOCK=<q> forces the chunk to q functions
   !! and TERCO_TDDFT_DF_RESIDENT=<m> keeps at most m functions resident, 0 none
   !! (TESTING: both force the blocked and the rebuilt paths on a small
   !! molecule; `tddft_df_set_blocks` is the same hook for a test program,
   !! which cannot set the environment).
   !!
   use pic_types, only: dp, int64
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_scf_df, only: df_jk_t
   use trc_scf_df_response, only: scf_df_response_unpack_columns
   use trc_device_memory, only: device_memory_free_bytes
   use trc_step_table, only: step_table_t
   use pic_timer, only: timer_type
   implicit none
   private

   public :: df_excitation_t, tddft_df_set_blocks
   public :: tddft_df_combine_tda, tddft_df_combine_tda_j, tddft_df_combine_rpa, tddft_df_combine_rpa_j
   public :: ST_DF_SETUP, ST_DF_STREAM, ST_DF_J, ST_DF_K1, ST_DF_K2, ST_DF_FXC, ST_DF_REST

   integer, parameter :: ST_DF_SETUP = 1, ST_DF_STREAM = 2, ST_DF_J = 3, ST_DF_K1 = 4, ST_DF_K2 = 5, &
                         ST_DF_FXC = 6, ST_DF_REST = 7, N_STEPS = 7
   !! The rows of `df_excitation_t%tbl`: the MO blocks at setup, the blocks rebuilt inside applications, the
   !! Coulomb, (ij|ab) and (ib|ja) products, the XC kernel (filled by trc_tddft) and everything else of an
   !! application (reorderings, combinations, the diagonal, the transfers of the trial vectors).
   character(len=24), parameter :: STEP_NAMES(N_STEPS) = [character(len=24) :: "MO blocks (setup)", &
                                                          "MO blocks (rebuilt)", "coulomb  (ia|jb)", &
                                                          "exchange (ij|ab)", "exchange (ib|ja)", &
                                                          "xc kernel", "rest of the matvec"]

   integer, parameter :: VMAX = 16
   !! Trial vectors per pass over the blocks: the Davidson's new vectors of an iteration are 5 to 30.
   integer(int64), parameter :: RESERVE = 536870912_int64
   !! Device bytes kept back for the libraries.
   integer(int64), parameter :: HOST_FREE = 8589934592_int64
   !! What the host build treats as free (it has no device to ask).

   integer, save :: force_qblk = 0
   integer, save :: force_resident = -1
   logical, save :: env_read = .false.

   type :: df_blocks_t
      !! The MO blocks of consecutive auxiliary functions, flat: for a chunk of np functions starting at p0,
      !! boo (nocc, nocc, np), bov (nocc, nvir, np) and bvv (nvir, np, nvir) (module header). Resident blocks hold
      !! every chunk, chunk c at the offset of its first function; the slot holds one chunk at the start.
      logical :: slot = .false.
      real(dp), allocatable :: boo(:), bov(:), bvv(:)
      !! Device-resident.
   end type df_blocks_t

   type :: df_plan_t
      !! The sizes `choose_blocks` plans from.
      integer :: n = 0, nocc = 0, nvir = 0, naux = 0
      logical :: exchange = .false.
   end type df_plan_t

   type :: df_excitation_t
      !! The fitted MO blocks of one fit and one set of orbitals, device-resident, and the scratch of an application.
      logical :: ready = .false.
      integer :: n = 0, nocc = 0, nvir = 0, naux = 0, npack = 0
      integer :: qblk = 0
      !! Auxiliary functions per chunk.
      integer :: nres = 0
      !! Auxiliary functions (a whole number of chunks, or all) whose blocks stay resident.
      logical :: exchange = .false.
      !! boo and bvv exist (a hybrid reference).
      real(dp), allocatable :: bpk(:, :)
      !! The fit's packed B, taken over from the `df_jk_t` by `setup` and handed back by `restore`; device-resident.
      real(dp), allocatable :: cmo(:, :)
      !! (n, nocc + nvir): the fit's AO orbitals, occupied first; device-resident.
      type(df_blocks_t) :: res, slot
      !! The resident chunks, and the slot a rebuilt chunk is made in (allocated only when some chunk is not resident).
      real(dp), allocatable :: tens(:, :, :), half(:)
      !! One chunk of unpacked B and its half-transform; device-resident, held only for the rebuilt chunks.
      real(dp), allocatable :: xk(:), w(:), rk(:), tq(:), k1t(:)
      !! An application's scratch; device-resident.
      real(dp), allocatable :: xg(:, :), jg(:, :), k1g(:, :), k2g(:, :)
      !! (nocc nvir, vmax): the pass's vectors and results; device-resident.
      type(step_table_t) :: tbl
   contains
      procedure :: setup => df_excitation_setup
      procedure :: apply => df_excitation_apply
      procedure :: release => df_excitation_release
      procedure :: restore => df_excitation_restore
      procedure :: report => df_excitation_report
   end type df_excitation_t

contains

   subroutine tddft_df_set_blocks(qblk, resident)
      !! TESTING: force the chunk to `qblk` auxiliary functions (0 lifts it) and keep at most `resident` functions'
      !! blocks resident (a negative value lifts it, 0 rebuilds every chunk in every application). Read by the next
      !! `setup`; the environment's TERCO_TDDFT_DF_BLOCK and TERCO_TDDFT_DF_RESIDENT are read the first time `setup`
      !! runs, this overrides them.
      integer, intent(in) :: qblk, resident
      force_qblk = qblk
      force_resident = resident
      env_read = .true.
   end subroutine tddft_df_set_blocks

   subroutine read_block_environment()
      !! The two testing switches of the module header, once.
      character(len=32) :: v
      integer :: st, q
      if (env_read) return
      env_read = .true.
      call get_environment_variable("TERCO_TDDFT_DF_BLOCK", v, status=st)
      if (st == 0) then
         read (v, *, iostat=st) q
         if (st == 0 .and. q >= 1) force_qblk = q
      end if
      call get_environment_variable("TERCO_TDDFT_DF_RESIDENT", v, status=st)
      if (st == 0) then
         read (v, *, iostat=st) q
         if (st == 0 .and. q >= 0) force_resident = q
      end if
   end subroutine read_block_environment

   subroutine choose_blocks(plan, extra, qblk, nres, err)
      !! The chunk and the resident functions of the module header, from the free device memory less `extra` (bytes the
      !! caller will hold beside this object). `err` carries the sizes when not even one auxiliary function fits.
      type(df_plan_t), intent(in) :: plan
      integer(int64), intent(in) :: extra
      integer, intent(out) :: qblk, nres
      type(error_t), intent(inout) :: err
      integer(int64) :: free, avail, per_res, per_apply, per_build, fixed, per_all, left
      integer(int64) :: n8, no8, nv8, nrow, v8
      integer :: nchunk

      n8 = plan%n
      no8 = plan%nocc
      nv8 = plan%nvir
      v8 = VMAX
      nrow = merge(no8 + nv8, no8, plan%exchange)
      per_res = 8_int64*(no8*nv8 + merge(no8*no8 + nv8*nv8, 0_int64, plan%exchange))
      per_apply = 8_int64*(no8*v8*nv8 + merge(nv8*nv8, 0_int64, plan%exchange) + v8)
      per_build = 8_int64*(n8*n8 + nrow*n8)
      fixed = 8_int64*(2_int64*no8*v8*nv8 + 4_int64*no8*nv8*v8)
      free = device_memory_free_bytes()
      if (free < 0) free = HOST_FREE
      avail = free - RESERVE - extra
      per_all = per_res + per_apply + per_build
      if (avail < fixed + per_all) then
         call err%set(ERROR_VALIDATION, "trc_tddft: the density-fitted excitation operator needs at least "// &
                      to_string(real(fixed + per_all, dp)/1.0e9_dp, "f8.3")//" GB of device memory (one auxiliary "// &
                      "function of "//to_string(plan%naux)//", "//to_string(plan%nocc)//" occupied and "// &
                      to_string(plan%nvir)//" virtual orbitals) and "// &
                      to_string(real(max(avail, 0_int64), dp)/1.0e9_dp, "f8.3")// &
                      " GB are free: use a smaller basis or a device with more memory")
         qblk = 0
         nres = 0
         return
      end if
      ! A third of what is left goes to the chunk's scratch, the rest to resident blocks.
      qblk = int(min(int(plan%naux, int64), max(1_int64, (avail/3 - fixed)/per_all)))
      call read_block_environment()
      if (force_qblk >= 1) qblk = min(qblk, force_qblk)
      nchunk = (plan%naux + qblk - 1)/qblk
      left = avail - fixed - int(qblk, int64)*per_all
      nres = int(min(int(nchunk, int64), left/(per_res*int(qblk, int64))))*qblk
      if (nres >= plan%naux) nres = plan%naux
      if (force_resident >= 0 .and. force_resident < plan%naux) nres = min(nres, (force_resident/qblk)*qblk)
   end subroutine choose_blocks

   subroutine allocate_blocks(blocks, plan, nfunc, slot)
      !! Storage for the blocks of `nfunc` auxiliary functions, on the device.
      type(df_blocks_t), intent(inout) :: blocks
      type(df_plan_t), intent(in) :: plan
      integer, intent(in) :: nfunc
      logical, intent(in) :: slot
      integer(int64) :: no8, nv8
      no8 = plan%nocc
      nv8 = plan%nvir
      blocks%slot = slot
      allocate (blocks%bov(max(1_int64, no8*nv8*nfunc)))
      !$acc enter data create(blocks%bov)
      if (plan%exchange) then
         allocate (blocks%boo(max(1_int64, no8*no8*nfunc)), blocks%bvv(max(1_int64, nv8*nv8*nfunc)))
         !$acc enter data create(blocks%boo, blocks%bvv)
      else
         ! Never read without exchange, but passed along: one element keeps the actual argument defined.
         allocate (blocks%boo(1), blocks%bvv(1))
      end if
   end subroutine allocate_blocks

   subroutine free_blocks(blocks, exchange)
      type(df_blocks_t), intent(inout) :: blocks
      logical, intent(in) :: exchange
      if (allocated(blocks%bov)) then
         !$acc exit data delete(blocks%bov)
         deallocate (blocks%bov)
      end if
      if (allocated(blocks%boo) .and. exchange) then
         !$acc exit data delete(blocks%boo, blocks%bvv)
      end if
      if (allocated(blocks%boo)) deallocate (blocks%boo)
      if (allocated(blocks%bvv)) deallocate (blocks%bvv)
   end subroutine free_blocks

   subroutine df_excitation_setup(this, dfjk, lp, cocc, cvir, exchange, extra_bytes, err)
      !! The MO blocks of `dfjk`'s fit for the orbitals `cocc` (n, nocc) and `cvir` (n, nvir) (host, the fit's AO
      !! space), exchange blocks only when `exchange`, `extra_bytes` of device memory left to the caller. The packed B
      !! moves from `dfjk` into this object until `restore`. `err` carries a memory refusal; nothing is taken then.
      class(df_excitation_t), intent(inout) :: this
      type(df_jk_t), intent(inout) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: cocc(:, :), cvir(:, :)
      logical, intent(in) :: exchange
      integer(int64), intent(in) :: extra_bytes
      type(error_t), intent(inout) :: err
      type(df_plan_t) :: plan
      integer(int64) :: no8, nv8, nrow
      integer :: n, naux, npack, no, nv, c, nchunk, p0, np
      type(timer_type) :: tc

      call this%release()
      if (allocated(this%bpk) .or. .not. allocated(dfjk%bpk)) then
         call err%set(ERROR_VALIDATION, "trc_tddft: the density-fitted excitation operator needs the fit's packed B "// &
                      "in place: it is lent to another operator, or the fit was released")
         return
      end if
      call dfjk%dims(n, naux, npack)
      no = size(cocc, 2)
      nv = size(cvir, 2)
      if (size(cocc, 1) /= n .or. size(cvir, 1) /= n) then
         call err%set(ERROR_VALIDATION, "trc_tddft: the orbitals are not in the fit's AO space (" &
                      //to_string(size(cocc, 1))//" functions against the fit's "//to_string(n)//")")
         return
      end if
      plan = df_plan_t(n, no, nv, naux, exchange)
      call choose_blocks(plan, extra_bytes, this%qblk, this%nres, err)
      if (err%has_error()) return
      this%exchange = exchange
      this%n = n
      this%nocc = no
      this%nvir = nv
      this%naux = naux
      this%npack = npack
      no8 = no
      nv8 = nv
      nrow = merge(no8 + nv8, no8, exchange)
      nchunk = (naux + this%qblk - 1)/this%qblk
      call global_logger%info("  tddft df: "//to_string(naux)//" auxiliary functions in "//to_string(nchunk)// &
                              " chunk(s) of "//to_string(this%qblk)//", "//to_string(this%nres)//" resident, exchange "// &
                              "blocks "//merge("on ", "off", exchange))
      call this%tbl%init("tddft df", STEP_NAMES)

      call move_alloc(dfjk%bpk, this%bpk)
      allocate (this%cmo(n, no + nv))
      this%cmo(:, 1:no) = cocc
      this%cmo(:, no + 1:no + nv) = cvir
      !$acc enter data copyin(this%cmo)
      call allocate_blocks(this%res, plan, this%nres, .false.)
      if (this%nres < naux) call allocate_blocks(this%slot, plan, this%qblk, .true.)
      allocate (this%tens(n, n, this%qblk), this%half(nrow*int(n, int64)*int(this%qblk, int64)))
      allocate (this%xg(no*nv, VMAX), this%jg(no*nv, VMAX), this%k1g(no*nv, VMAX), this%k2g(no*nv, VMAX))
      !$acc enter data create(this%tens, this%half, this%xg, this%jg, this%k1g, this%k2g)
      this%ready = .true.

      call this%tbl%step_on(lp, ST_DF_SETUP, tc)
      do c = 1, nchunk
         p0 = (c - 1)*this%qblk + 1
         np = min(this%qblk, naux - p0 + 1)
         if (p0 + np - 1 > this%nres) exit
         call build_chunk(this, lp, this%res, p0, np)
      end do
      call this%tbl%step_off(lp, ST_DF_SETUP, tc)
      ! The unpack and half-transform scratch is needed again only by chunks rebuilt inside an application.
      if (this%nres >= naux) then
         !$acc exit data delete(this%tens, this%half)
         deallocate (this%tens, this%half)
      end if
   end subroutine df_excitation_setup

   subroutine build_chunk(this, lp, blocks, p0, np)
      !! The blocks of auxiliary functions p0 .. p0 + np - 1 into `blocks` (at their own position when resident, at
      !! the start for the slot): one unpack, one half-transform GEMM, three batched GEMMs (module header).
      type(df_excitation_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: lp
      type(df_blocks_t), intent(inout) :: blocks
      integer, intent(in) :: p0, np
      type(trc_gemm_batch_t) :: gb
      integer :: n, no, nv, nrow
      integer(int64) :: n8, nrow8, first, ob, ov, vv

      n = this%n
      no = this%nocc
      nv = this%nvir
      n8 = n
      nrow = merge(no + nv, no, this%exchange)
      nrow8 = nrow
      first = merge(0, p0 - 1, blocks%slot)
      ob = first*no*no
      ov = first*no*nv
      vv = first*nv*nv
      call scf_df_response_unpack_columns(n, this%npack, this%naux, np, p0, this%bpk, this%tens)
      ! half(p, nu q) = sum_mu C(mu, p) B^q(mu, nu), p over the occupied (and virtual) orbitals.
      call lp%gemm("T", "N", nrow, n*np, n, 1.0_dp, this%cmo, n, this%tens, n, 0.0_dp, this%half, nrow)
      ! B^q_ia = half_o^q C_v
      gb = trc_gemm_batch(no, nv, n, np, lda=nrow, ldb=n, ldc=no, sa=nrow8*n8, sb=0_int64, sc=int(no, int64)*nv)
      call lp%gemm_batched(gb, this%half, this%cmo(1, no + 1), blocks%bov(ov + 1))
      if (.not. this%exchange) return
      ! B^q_ij = half_o^q C_o
      gb = trc_gemm_batch(no, no, n, np, lda=nrow, ldb=n, ldc=no, sa=nrow8*n8, sb=0_int64, sc=int(no, int64)**2)
      call lp%gemm_batched(gb, this%half, this%cmo, blocks%boo(ob + 1))
      ! B^q_ab = half_v^q C_v, written as (b, q, a): the (b, q) pair contiguous for the contractions
      gb = trc_gemm_batch(nv, nv, n, np, lda=nrow, ldb=n, ldc=nv*np, sa=nrow8*n8, sb=0_int64, sc=int(nv, int64))
      call lp%gemm_batched(gb, this%half(no + 1), this%cmo(1, no + 1), blocks%bvv(vv + 1))
   end subroutine build_chunk

   subroutine df_excitation_apply(this, lp, nvec, x, jx, k1, k2, want_k2)
      !! jx(:, k) = J[x_k] and, with exchange, k1(:, k) = K1[x_k] and (`want_k2`) k2(:, k) = K2[x_k] for the nvec
      !! vectors x(:, k) (nocc nvir flattened), all device-resident, through passes of at most VMAX vectors over the
      !! blocks (module header). `k1` is not touched without exchange, `k2` unless `want_k2`.
      class(df_excitation_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nvec
      real(dp), intent(in) :: x(this%nocc*this%nvir, nvec)
      real(dp), intent(inout) :: jx(this%nocc*this%nvir, nvec)
      real(dp), intent(inout), optional :: k1(this%nocc*this%nvir, nvec), k2(this%nocc*this%nvir, nvec)
      !! `k1` is needed with exchange, `k2` with `want_k2`; a pure functional's caller allocates neither.
      logical, intent(in) :: want_k2
      integer :: g0, g

      if (.not. this%ready) error stop "trc_tddft_df: apply before setup"
      do g0 = 0, nvec - 1, VMAX
         g = min(VMAX, nvec - g0)
         if (want_k2) then
            call apply_group(this, lp, g, x(1, g0 + 1), jx(1, g0 + 1), k1(1, g0 + 1), k2(1, g0 + 1), want_k2)
         else if (this%exchange) then
            call apply_group(this, lp, g, x(1, g0 + 1), jx(1, g0 + 1), k1(1, g0 + 1), want_k2=want_k2)
         else
            call apply_group(this, lp, g, x(1, g0 + 1), jx(1, g0 + 1), want_k2=want_k2)
         end if
      end do
   end subroutine df_excitation_apply

   subroutine apply_group(this, lp, g, x, jx, k1, k2, want_k2)
      !! One pass over the blocks for g <= VMAX vectors.
      type(df_excitation_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: g
      real(dp), intent(in) :: x(this%nocc*this%nvir, g)
      real(dp), intent(inout) :: jx(this%nocc*this%nvir, g)
      real(dp), intent(inout), optional :: k1(this%nocc*this%nvir, g)
      real(dp), intent(inout), optional :: k2(this%nocc*this%nvir, g)
      logical, intent(in) :: want_k2
      integer :: c, nchunk, p0, np, no, nv
      integer(int64) :: no8, nv8
      logical :: resident
      type(timer_type) :: tc

      no = this%nocc
      nv = this%nvir
      no8 = no
      nv8 = nv
      nchunk = (this%naux + this%qblk - 1)/this%qblk
      if (this%exchange .and. .not. allocated(this%xk)) then
         allocate (this%xk(no8*VMAX*nv8), this%k1t(no8*VMAX*nv8), this%w(no8*VMAX*nv8*this%qblk), &
                   this%rk(nv8*nv8*this%qblk))
         !$acc enter data create(this%xk, this%k1t, this%w, this%rk)
      end if
      if (.not. allocated(this%tq)) then
         allocate (this%tq(this%qblk*VMAX))
         !$acc enter data create(this%tq)
      end if
      call this%tbl%step_on(lp, ST_DF_REST, tc)
      call copy_vectors(no*nv, g, x, this%xg)
      if (this%exchange) call reorder_x(no, nv, g, this%xg, this%xk)
      call this%tbl%step_off(lp, ST_DF_REST, tc)
      do c = 1, nchunk
         p0 = (c - 1)*this%qblk + 1
         np = min(this%qblk, this%naux - p0 + 1)
         resident = p0 + np - 1 <= this%nres
         if (resident) then
            call chunk_products(this, lp, this%res, p0, np, g, want_k2)
         else
            call this%tbl%step_on(lp, ST_DF_STREAM, tc)
            call build_chunk(this, lp, this%slot, p0, np)
            call this%tbl%step_off(lp, ST_DF_STREAM, tc)
            call chunk_products(this, lp, this%slot, p0, np, g, want_k2)
         end if
      end do
      call this%tbl%step_on(lp, ST_DF_REST, tc)
      if (this%exchange) call restore_k1(no, nv, g, this%k1t, this%k1g)
      call copy_vectors(no*nv, g, this%jg, jx)
      if (this%exchange) call copy_vectors(no*nv, g, this%k1g, k1)
      if (this%exchange .and. want_k2) call copy_vectors(no*nv, g, this%k2g, k2)
      call this%tbl%step_off(lp, ST_DF_REST, tc)
   end subroutine apply_group

   subroutine chunk_products(this, lp, blocks, p0, np, g, want_k2)
      !! One chunk's contributions of np auxiliary functions to J, K1 and K2 of the g vectors in `xg` (module header).
      !! The first chunk starts the accumulations (beta = 0), the others add.
      type(df_excitation_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: lp
      type(df_blocks_t), intent(inout) :: blocks
      integer, intent(in) :: p0, np, g
      logical, intent(in) :: want_k2
      type(trc_gemm_batch_t) :: gb
      type(timer_type) :: tc
      integer :: no, nv, nov, k
      integer(int64) :: first, ob, ov, vv
      real(dp) :: beta

      no = this%nocc
      nv = this%nvir
      nov = no*nv
      first = merge(0, p0 - 1, blocks%slot)
      ob = first*no*no
      ov = first*no*nv
      vv = first*nv*nv
      beta = merge(0.0_dp, 1.0_dp, p0 == 1)
      ! J: t = B_ov^T x, J += B_ov t
      call this%tbl%step_on(lp, ST_DF_J, tc)
      call lp%gemm("T", "N", np, g, nov, 1.0_dp, blocks%bov(ov + 1), nov, this%xg, nov, 0.0_dp, this%tq, np)
      call lp%gemm("N", "N", nov, g, np, 1.0_dp, blocks%bov(ov + 1), nov, this%tq, np, beta, this%jg, nov)
      call this%tbl%step_off(lp, ST_DF_J, tc)
      if (.not. this%exchange) return
      ! K1: W(i, k, b, Q) = sum_j B_oo(j, i, Q) x(j, k, b), K1(i, k, a) += W((i, k), (b, Q)) B_vv((b, Q), a)
      call this%tbl%step_on(lp, ST_DF_K1, tc)
      gb = trc_gemm_batch(no, g*nv, no, np, lda=no, ldb=no, ldc=no, sa=int(no, int64)**2, sb=0_int64, &
                          sc=int(no, int64)*g*nv, transa="T")
      call lp%gemm_batched(gb, blocks%boo(ob + 1), this%xk, this%w)
      call lp%gemm("N", "N", no*g, nv, nv*np, 1.0_dp, this%w, no*g, blocks%bvv(vv + 1), nv*np, beta, this%k1t, no*g)
      call this%tbl%step_off(lp, ST_DF_K1, tc)
      if (.not. want_k2) return
      ! K2, vector by vector: R(b, Q, a) = sum_j x(j, b) B_ov(j, a, Q), K2(:, :, k) += B_ov(i, (b, Q)) R((b, Q), a)
      call this%tbl%step_on(lp, ST_DF_K2, tc)
      gb = trc_gemm_batch(nv, nv, no, np, lda=no, ldb=no, ldc=nv*np, sa=0_int64, sb=int(no, int64)*nv, &
                          sc=int(nv, int64), transa="T")
      do k = 1, g
         call lp%gemm_batched(gb, this%xg(1, k), blocks%bov(ov + 1), this%rk)
         call lp%gemm("N", "N", no, nv, nv*np, 1.0_dp, blocks%bov(ov + 1), no, this%rk, nv*np, beta, this%k2g(1, k), no)
      end do
      call this%tbl%step_off(lp, ST_DF_K2, tc)
   end subroutine chunk_products

   subroutine copy_vectors(nov, g, a, b)
      !! b(:, 1:g) = a(:, 1:g), on the device: the vectors of a pass into and out of the object's scratch.
      integer, intent(in) :: nov, g
      real(dp), intent(in) :: a(nov, g)
      real(dp), intent(out) :: b(nov, g)
      integer(int64) :: t
      integer :: p, k
      do concurrent(t=0:int(nov, int64)*g - 1) local(p, k)
         p = int(mod(t, int(nov, int64))) + 1
         k = int(t/int(nov, int64)) + 1
         b(p, k) = a(p, k)
      end do
   end subroutine copy_vectors

   subroutine reorder_x(no, nv, g, x, xk)
      !! xk(j, k, b) = x(j, b, k): the vectors beside each other along the second index, so one GEMM
      !! against a shared block serves every vector.
      integer, intent(in) :: no, nv, g
      real(dp), intent(in) :: x(no*nv*g)
      real(dp), intent(out) :: xk(no*g*nv)
      integer(int64) :: t, no8, nv8
      integer :: j, b, k
      no8 = no
      nv8 = nv
      do concurrent(t=0:no8*nv8*g - 1) local(j, b, k)
         j = int(mod(t, no8)) + 1
         b = int(mod(t/no8, nv8)) + 1
         k = int(t/(no8*nv8)) + 1
         xk(j + no*(k - 1) + no*g*(b - 1)) = x(t + 1)
      end do
   end subroutine reorder_x

   subroutine restore_k1(no, nv, g, k1t, k1)
      !! k1(i, a, k) = k1t(i, k, a), the inverse of `reorder_x`'s layout.
      integer, intent(in) :: no, nv, g
      real(dp), intent(in) :: k1t(no*g*nv)
      real(dp), intent(out) :: k1(no*nv*g)
      integer(int64) :: t, no8, nv8
      integer :: i, a, k
      no8 = no
      nv8 = nv
      do concurrent(t=0:no8*nv8*g - 1) local(i, a, k)
         i = int(mod(t, no8)) + 1
         a = int(mod(t/no8, nv8)) + 1
         k = int(t/(no8*nv8)) + 1
         k1(t + 1) = k1t(i + no*(k - 1) + no*g*(a - 1))
      end do
   end subroutine restore_k1

   subroutine tddft_df_combine_tda(nov, nvec, jscale, kscale, jx, k1, ax)
      !! ax := jscale J + (-kscale) K1 over the amplitudes flattened to one index, on the device.
      integer, intent(in) :: nov, nvec
      real(dp), intent(in) :: jscale, kscale
      real(dp), intent(in) :: jx(nov, nvec), k1(nov, nvec)
      real(dp), intent(out) :: ax(nov, nvec)
      integer(int64) :: t
      integer :: p, k
      do concurrent(t=0:int(nov, int64)*nvec - 1) local(p, k)
         p = int(mod(t, int(nov, int64))) + 1
         k = int(t/int(nov, int64)) + 1
         ax(p, k) = jscale*jx(p, k) - kscale*k1(p, k)
      end do
   end subroutine tddft_df_combine_tda

   subroutine tddft_df_combine_tda_j(nov, nvec, jscale, jx, ax)
      !! ax := jscale J (a pure functional: no exchange), on the device.
      integer, intent(in) :: nov, nvec
      real(dp), intent(in) :: jscale
      real(dp), intent(in) :: jx(nov, nvec)
      real(dp), intent(out) :: ax(nov, nvec)
      integer(int64) :: t
      integer :: p, k
      do concurrent(t=0:int(nov, int64)*nvec - 1) local(p, k)
         p = int(mod(t, int(nov, int64))) + 1
         k = int(t/int(nov, int64)) + 1
         ax(p, k) = jscale*jx(p, k)
      end do
   end subroutine tddft_df_combine_tda_j

   subroutine tddft_df_combine_rpa(nov, nvec, cx, jx, k1, k2, ax)
      !! ax(:, :, 1) := 4 J - c_x (K1 + K2) (A+B without its diagonal and XC kernel) and ax(:, :, 2) :=
      !! -c_x (K1 - K2) (A-B without its diagonal), on the device.
      integer, intent(in) :: nov, nvec
      real(dp), intent(in) :: cx
      real(dp), intent(in) :: jx(nov, nvec), k1(nov, nvec), k2(nov, nvec)
      real(dp), intent(out) :: ax(nov, nvec, 2)
      integer(int64) :: t
      integer :: p, k
      do concurrent(t=0:int(nov, int64)*nvec - 1) local(p, k)
         p = int(mod(t, int(nov, int64))) + 1
         k = int(t/int(nov, int64)) + 1
         ax(p, k, 1) = 4.0_dp*jx(p, k) - cx*(k1(p, k) + k2(p, k))
         ax(p, k, 2) = -cx*(k1(p, k) - k2(p, k))
      end do
   end subroutine tddft_df_combine_rpa

   subroutine tddft_df_combine_rpa_j(nov, nvec, jx, ax)
      !! ax(:, :, 1) := 4 J and ax(:, :, 2) := 0 (a pure functional: A-B is the diagonal), on the device.
      integer, intent(in) :: nov, nvec
      real(dp), intent(in) :: jx(nov, nvec)
      real(dp), intent(out) :: ax(nov, nvec, 2)
      integer(int64) :: t
      integer :: p, k
      do concurrent(t=0:int(nov, int64)*nvec - 1) local(p, k)
         p = int(mod(t, int(nov, int64))) + 1
         k = int(t/int(nov, int64)) + 1
         ax(p, k, 1) = 4.0_dp*jx(p, k)
         ax(p, k, 2) = 0.0_dp
      end do
   end subroutine tddft_df_combine_rpa_j

   subroutine df_excitation_report(this, t_total)
      !! The step table, one line per row.
      class(df_excitation_t), intent(in) :: this
      real(dp), intent(in) :: t_total
      if (allocated(this%tbl%names)) call this%tbl%report(t_total, 0, 1)
   end subroutine df_excitation_report

   subroutine df_excitation_restore(this, dfjk)
      !! Hand the packed B back to the `df_jk_t` it came from and free everything else. Safe on an object that was
      !! never set up. Call it once the solve is over, on every path.
      class(df_excitation_t), intent(inout) :: this
      type(df_jk_t), intent(inout) :: dfjk
      call this%release()
      if (allocated(this%bpk)) call move_alloc(this%bpk, dfjk%bpk)
   end subroutine df_excitation_restore

   subroutine df_excitation_release(this)
      !! Free the blocks and the scratch. The packed B stays with this object until `restore` (it is not ours to free).
      class(df_excitation_t), intent(inout) :: this
      if (this%ready) then
         !$acc exit data delete(this%cmo, this%xg, this%jg, this%k1g, this%k2g)
         if (allocated(this%tens)) then
            !$acc exit data delete(this%tens, this%half)
         end if
         if (allocated(this%xk)) then
            !$acc exit data delete(this%xk, this%k1t, this%w, this%rk)
         end if
         if (allocated(this%tq)) then
            !$acc exit data delete(this%tq)
         end if
      end if
      call free_blocks(this%res, this%exchange)
      call free_blocks(this%slot, this%exchange)
      if (allocated(this%cmo)) deallocate (this%cmo)
      if (allocated(this%tens)) deallocate (this%tens)
      if (allocated(this%half)) deallocate (this%half)
      if (allocated(this%xk)) deallocate (this%xk)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%rk)) deallocate (this%rk)
      if (allocated(this%tq)) deallocate (this%tq)
      if (allocated(this%k1t)) deallocate (this%k1t)
      if (allocated(this%xg)) deallocate (this%xg, this%jg, this%k1g, this%k2g)
      this%ready = .false.
      this%exchange = .false.
      this%nres = 0
      this%qblk = 0
   end subroutine df_excitation_release

end module trc_tddft_df
