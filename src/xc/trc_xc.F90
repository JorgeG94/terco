!
! The exchange-correlation integrator: a density in, an energy and a
! potential matrix out, everything on the grid done inside `do concurrent`.
!
! For a closed-shell density D and a functional f(rho, sigma):
!
!     E_xc  = sum_g w_g rho_g eps(rho_g, sigma_g)
!     V_uv  = sum_g w_g [ v_rho chi_u chi_v
!                       + 2 v_sigma grad rho . (grad chi_u chi_v + chi_u grad chi_v) ]
!
! with rho_g = sum_uv D_uv chi_u chi_v and sigma = |grad rho|^2. The factor
! of two on the sigma term and the symmetrisation are the usual place a GGA
! goes wrong by a few millihartree while converging perfectly; check_xc_energy
! compares V_uv element by element against libxc for that reason.
!
! For a spin-polarised pair (D_a, D_b) the functional takes rho_a, rho_b
! and the three invariants sigma_aa, sigma_ab, sigma_bb, and the potential
! for spin a carries
!
!     2 v_aa grad rho_a + v_ab grad rho_b
!
! in place of 2 v_sigma grad rho, with a and b exchanged for spin b. The
! same two kernels serve both; `nspin` is the only difference.
!
! A meta-GGA also depends on the kinetic-energy density of each spin,
!
!     tau = 1/2 sum_i |grad psi_i|^2 = 1/2 sum_uv D_uv grad chi_u . grad chi_v
!
! with D the total density for RKS and the spin's own for UKS -- libxc's
! and pyscf's convention, so tau of a closed shell is the total and not
! half of it. Its potential adds
!
!     V_uv += sum_g w_g v_tau/2 grad chi_u . grad chi_v
!
! from d tau / d D_uv. Both come out of the gradient tiles a GGA already
! collocates: see THE PIPELINE PER CHUNK.
!
! THE PIPELINE PER CHUNK
! ----------------------
! Batched dense linear algebra, the "Batch Dense D" algorithm of Stocks
! and Barca (JCTC 2025, 21, 10263), which is also GauXC's. Per batch b of
! points with local basis of nloc functions, four kernels and two GEMMs:
!
!   1. collocate:  chi_b (npts x nloc) and its three gradient tiles,
!                  one thread per point, evaluated ONCE and reused below
!   2. gather:     D_b (nloc x nloc) sliced out of D through b_ao
!      GEMM:       X_b = chi_b D_b
!   3. points:     rho = sum_u X_gu chi_gu, grad rho = 2 sum_u X_gu grad chi_gu,
!                  the functional, and the scaled tile
!                  Z_gu = w (v_rho/2 chi_gu + 2 v_sigma grad rho . grad chi_gu)
!      GEMM:       V_b = chi_b^T Z_b
!   4. scatter:    V_uv += V_b(u,v) + V_b(v,u), one thread per (b, u, v),
!                  atomic because batches share global functions
!
! A meta-GGA adds three GEMMs per spin on each side, nothing else: before
! the point kernel Y_x = (d chi/dx)_b D_b, one per Cartesian direction, so
! that tau = 1/2 sum_x sum_u Y_x,gu (d chi/dx)_gu; and after it, in Y_x's own
! storage (each point thread reads its row of Y before overwriting it with
! W_x,gu = w v_tau/4 (d chi/dx)_gu),
!
!     V_b += sum_x (d chi/dx)_b^T W_x
!
! accumulated onto the chi^T Z product. W_x is symmetric under the
! scatter's V_b + V_b^T, hence the quarter. LDA and GGA functionals take
! none of this and run exactly the kernels they always did.
!
! The earlier form contracted D per point and summed the potential per
! pair in loop nests, their "Direct" algorithm, and on cholesterol that was
! 0.6 s per iteration where the GEMM form is expected under 0.1. The
! GEMMs are cuBLAS from OpenACC through trc_linalg, pic-blas on the host;
! nothing else on the path is a library.
!
! One cuBLAS call per batch cost 37 us each, three quarters of the time.
! So the tiles are padded -- points up to a multiple of XC_NP_PAD,
! functions up to a multiple of XC_NL_PAD, padding zero -- and the batches
! of a chunk are STORED in order of their padded shape, so that every run
! of equal shape is one strided-batched GEMM. The padding costs flops at
! the rate the shapes vary; the grouping removes the launches. Storage
! order is not batch order: c_off(ib) says where batch ib's tile is, and
! the pair kernels walk storage order through `ord`.
!
! The loop bodies live in `routine seq` procedures and the loops only call
! them. That is terco's pattern throughout, and it is not cosmetic: given
! the body inline, nvfortran parallelises the inner loops over the local
! basis across the threads of a block and hands each POINT a whole block,
! which is the wrong mapping by two orders of magnitude.
!
! CHUNKS
! ------
! The stored tiles are n_points x n_local per batch, which for a large
! system exceeds any device. Batches are therefore processed in chunks
! whose value tile fits `budget` doubles (the whole set of tiles is about
! ten times that), and the pipeline runs once per chunk. Blocking is designed in rather than retrofitted: retrofitting
! it is the coupled-cluster lesson metalquicha already paid for.
!
! DENSITY SCREENING
! -----------------
! Before the chunk loop, every batch is bounded: with m_u the largest
! value of local function u on the batch (recorded at batching time),
!
!     rho on the batch  <=  sum_uv |D_uv| m_u m_v
!
! and a batch whose bound is below `rho_tol` is skipped by both kernels
! and takes no storage. This is what turns the cost from every batch
! against every nearby function into every batch against the density
! that reaches it, and on a large molecule it is the difference in
! scaling rather than in constant. It costs one small matrix-vector per
! batch on the host per call.
!
! RANKS
! -----
! With a communicator the batches are split across ranks: sorted by their
! cost -- points times local functions times (4 + local functions), the
! point kernel plus the pair kernel -- and dealt round-robin, so each rank
! carries the same work to within one batch. A batch another rank owns is
! simply inactive here, the same flag the density screening uses, so the
! kernels never see it and it takes no storage. Afterwards V, E_xc and the
! electron count are summed across ranks, and V is broadcast from rank 0
! so every rank holds the identical matrix: an allreduce does not promise
! the same bits everywhere and the pair kernel's atomics do not either.
!
! WHAT RUNS ON THE HOST
! ---------------------
! The sums for E_xc and the electron count. `do concurrent` has no portable
! reduction -- gfortran took `reduce` after locality specifiers, and terco
! builds on the host with what a CI runner has -- so the per-point terms
! come back as arrays and are summed here. They are one double per point;
! the traffic is nothing next to the values.
!
module trc_xc
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_eval_point, xc_eval_point_polar, &
      xc_eval_point_fxc, xc_eval_point_mgga, xc_eval_point_mgga_polar, XC_MAX_KERNELS, XC_FAMILY_MGGA
   use trc_collocation, only: shell_collocate, NCART_MAX
   use trc_linalg, only: trc_linalg_t
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   implicit none
   private

   public :: trc_xc_rks, trc_xc_uks, trc_xc_fxc_rks, XC_RHO_TOL
   public :: XC_NP_PAD, XC_NL_PAD, sort_by_shape

   real(dp), parameter :: XC_RHO_TOL = 1.0e-12_dp
   !! A batch whose bound on rho is below this is skipped. The functionals'
   !! own density thresholds are 1e-12 to 1e-32; check_xc_energy measures
   !! what this value costs against an unscreened evaluation.
   integer, parameter :: XC_NP_PAD = 128   !! tile rows are padded to a multiple of this
   integer, parameter :: XC_NL_PAD = 32    !! and tile columns to a multiple of this

contains

   !
   ! E_xc, V_xc and the integrated electron count for a closed-shell density.
   !
   ! `budget` caps the stored basis values per chunk, in doubles; the
   ! gradients take three times as much again. The default is 8M doubles,
   ! 256 MB with the gradients; a meta-GGA adds three more tiles per spin.
   !
   subroutine trc_xc_rks(b, xg, func, dmat, vxc, exc, nelec, budget, t_points, t_pairs, rho_tol, n_skipped, comm, la)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      real(dp), intent(in) :: dmat(b%nao, b%nao)
      real(dp), intent(out) :: vxc(b%nao, b%nao)
      real(dp), intent(out) :: exc, nelec
      integer(int64), intent(in), optional :: budget
      real(dp), intent(out), optional :: t_points, t_pairs
      !! Wall time spent in the two kernels, for the benchmark. A `do
      !! concurrent` under -stdpar returns when the device is done, so the
      !! clock around it is the kernel's.
      real(dp), intent(in), optional :: rho_tol
      !! Density screening threshold on the per-batch bound of rho; zero
      !! screens nothing. Default XC_RHO_TOL.
      integer, intent(out), optional :: n_skipped   !! batches the screening removed
      type(comm_t), intent(in), optional :: comm     !! split the batches over its ranks
      type(trc_linalg_t), intent(inout), optional :: la  !! GEMM handles; made and released here if absent
      real(dp) :: tp, tq, rtol
      integer(int64) :: cap
      integer :: nsk

      cap = 8388608_int64
      if (present(budget)) cap = max(budget, 1_int64)
      rtol = XC_RHO_TOL
      if (present(rho_tol)) rtol = rho_tol
      call xc_drive(b, xg, func, 1, dmat, vxc, exc, nelec, cap, rtol, tp, tq, nsk, comm, la)
      if (present(t_points)) t_points = tp
      if (present(t_pairs)) t_pairs = tq
      if (present(n_skipped)) n_skipped = nsk
   end subroutine trc_xc_rks

   !
   ! The same for a spin-polarised density: D(:,:,1) is alpha, D(:,:,2)
   ! beta, and the potentials come back in the same slots.
   !
   subroutine trc_xc_uks(b, xg, func, dmat, vxc, exc, nelec, budget, t_points, t_pairs, rho_tol, n_skipped, comm, la)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      real(dp), intent(in) :: dmat(b%nao, b%nao, 2)
      real(dp), intent(out) :: vxc(b%nao, b%nao, 2)
      real(dp), intent(out) :: exc, nelec
      integer(int64), intent(in), optional :: budget
      real(dp), intent(out), optional :: t_points, t_pairs
      real(dp), intent(in), optional :: rho_tol
      integer, intent(out), optional :: n_skipped
      type(comm_t), intent(in), optional :: comm
      type(trc_linalg_t), intent(inout), optional :: la
      real(dp) :: tp, tq, rtol
      integer(int64) :: cap
      integer :: nsk

      cap = 8388608_int64
      if (present(budget)) cap = max(budget, 1_int64)
      rtol = XC_RHO_TOL
      if (present(rho_tol)) rtol = rho_tol
      call xc_drive(b, xg, func, 2, dmat, vxc, exc, nelec, cap, rtol, tp, tq, nsk, comm, la)
      if (present(t_points)) t_points = tp
      if (present(t_pairs)) t_pairs = tq
      if (present(n_skipped)) n_skipped = nsk
   end subroutine trc_xc_uks

   subroutine xc_drive(b, xg, func, nspin, dmat, vxc, exc, nelec, cap, rho_tol, t_points, t_pairs, n_skipped, comm, la)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      integer, intent(in) :: nspin
      real(dp), intent(in) :: dmat(b%nao, b%nao, nspin)
      real(dp), intent(out) :: vxc(b%nao, b%nao, nspin)
      real(dp), intent(out) :: exc, nelec, t_points, t_pairs
      integer(int64), intent(in) :: cap
      real(dp), intent(in) :: rho_tol
      integer, intent(out) :: n_skipped
      type(comm_t), intent(in), optional :: comm
      type(trc_linalg_t), intent(inout), optional, target :: la

      integer(int64) :: need, nchi, sa, sb, sc
      integer :: b0, b1, ib, nloc, npb, np, p0, npairs, rank, nranks, k, s, nb, k0, k1, cnt, mp, ml
      integer, allocatable :: c_off(:), pr_off(:), active(:), order(:), npp(:), nlp(:), ord(:), pb_off(:), wk_off(:)
      integer :: nwork
      real(dp), allocatable :: cost(:)
      real(dp) :: es(2)
      logical :: split
      real(dp), allocatable :: chi(:), gchi(:, :), xt(:, :), zt(:, :), dloc(:, :), vloc(:, :), eg(:), ng(:)
      real(dp), allocatable :: yt(:, :, :)
      !! meta-GGA only: (d chi/dx) D per direction and spin, then the tau
      !! potential tile in the same storage (see the module header)
      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      integer :: kid(XC_MAX_KERNELS)
      real(dp) :: coef(XC_MAX_KERNELS)
      real(dp) :: t0, t1, t2
      logical :: mgga
      integer :: x

      mgga = func%family == XC_FAMILY_MGGA
      kid = func%kid
      coef = func%coef
      vxc = 0.0_dp
      exc = 0.0_dp
      nelec = 0.0_dp
      t_points = 0.0_dp
      t_pairs = 0.0_dp
      if (xg%npts == 0) return

      allocate (c_off(xg%nbatch + 1), pr_off(xg%nbatch + 1), active(xg%nbatch))
      allocate (npp(xg%nbatch), nlp(xg%nbatch))

      ! Ranks: cost-sorted round-robin ownership. A batch another rank owns
      ! starts inactive here and is never even bounded.
      rank = 0
      nranks = 1
      split = .false.
      if (present(comm)) then
         rank = comm%rank()
         nranks = comm%size()
         split = nranks > 1
      end if
      active = 1
      if (split) then
         allocate (cost(xg%nbatch), order(xg%nbatch))
         do ib = 1, xg%nbatch
            nloc = xg%b_aooff(ib + 1) - xg%b_aooff(ib)
            npb = xg%b_off(ib + 1) - xg%b_off(ib)
            cost(ib) = real(npb, dp)*real(nloc, dp)*real(4 + nloc, dp)
            order(ib) = ib
         end do
         call sort_desc(xg%nbatch, cost, order)
         do k = 1, xg%nbatch
            if (mod(k - 1, nranks) /= rank) active(order(k)) = 0
         end do
      end if

      !$acc data copyin(dmat, kid, coef) copy(active, vxc) create(npp, nlp)
      ! Density screening: bound rho on each owned batch, on the device where
      ! the density is, and drop what cannot reach the tolerance.
      call xc_screen(xg%nbatch, nspin, b%nao, dmat, xg%b_aooff, size(xg%b_ao), xg%b_ao, xg%b_amax, &
                     rho_tol, active)
      !$acc update self(active)
      ! What the screen removed, not what the other ranks own.
      n_skipped = count(active == 0) - (xg%nbatch - n_owned(xg%nbatch, rank, nranks))

      ! A caller that keeps vxc resident (the SCF) has its device copy
      ! untouched by the host zeroing above; the accumulation must start
      ! from zero there too. Inside the region vxc is present either way.
      !$acc update device(vxc)
      if (present(la)) then
         lp => la
      else
         call la_own%init(1)
         lp => la_own
      end if
      ! Padded tile shape of every batch; inactive ones have no functions.
      do ib = 1, xg%nbatch
         nloc = (xg%b_aooff(ib + 1) - xg%b_aooff(ib))*active(ib)
         npb = xg%b_off(ib + 1) - xg%b_off(ib)
         npp(ib) = ((npb + XC_NP_PAD - 1)/XC_NP_PAD)*XC_NP_PAD
         nlp(ib) = ((nloc + XC_NL_PAD - 1)/XC_NL_PAD)*XC_NL_PAD
      end do
      !$acc update device(npp, nlp)
      b0 = 1
      do while (b0 <= xg%nbatch)
         ! The chunk: as many batches as the budget holds, and at least one.
         need = 0
         b1 = b0 - 1
         do ib = b0, xg%nbatch
            if (need + int(npp(ib), int64)*int(nlp(ib), int64) > cap .and. ib > b0) exit
            need = need + int(npp(ib), int64)*int(nlp(ib), int64)
            b1 = ib
         end do
         nchi = max(need, 1_int64)
         p0 = xg%b_off(b0)
         np = xg%b_off(b1 + 1) - p0
         nb = b1 - b0 + 1
         ! Storage order: by padded shape, so equal shapes are contiguous.
         allocate (ord(nb), pb_off(nb + 1), wk_off(nb + 1))
         do k = 1, nb
            ord(k) = b0 + k - 1
         end do
         call sort_by_shape(nb, npp, nlp, ord)
         ! Offsets in storage order: tiles (point fastest), and (u,v) pairs.
         need = 0
         npairs = 0
         nwork = 0
         pb_off(1) = 0
         wk_off(1) = 0
         do k = 1, nb
            ib = ord(k)
            c_off(ib) = int(need)
            pr_off(ib) = npairs
            need = need + int(npp(ib), int64)*int(nlp(ib), int64)
            npairs = npairs + nlp(ib)*nlp(ib)
            pb_off(k + 1) = npairs
            ! Collocation work items: (shell, point) of every active batch,
            ! point fastest, so a warp holds one shell at consecutive points.
            nwork = nwork + (xg%b_shoff(ib + 1) - xg%b_shoff(ib))*(xg%b_off(ib + 1) - xg%b_off(ib))*active(ib)
            wk_off(k + 1) = nwork
         end do

         allocate (chi(nchi), gchi(nchi, 3), xt(nchi, nspin), zt(nchi, nspin))
         allocate (dloc(max(npairs, 1), nspin), vloc(max(npairs, 1), nspin), eg(np), ng(np))
         !$acc enter data create(chi, gchi, xt, zt, dloc, vloc, eg, ng) copyin(c_off, pr_off, ord, pb_off, wk_off)
         if (mgga) then
            allocate (yt(nchi, 3, nspin))
            !$acc enter data create(yt)
         end if

         t0 = wall()
         ! Padding rows and columns must be zero, not whatever the allocation
         ! held: a NaN there would poison the GEMM. chi feeds both GEMMs, Z
         ! the second; X, D_b and V_b are written in full.
         call xc_zero(nchi, nspin, chi, zt)
         ! A meta-GGA's GEMMs read grad chi too, so its padding must be zero
         ! as well; Y and W then come out zero there without further help.
         if (mgga) call xc_zero_grad(nchi, gchi)
         call xc_collocate(nwork, nb, ord, wk_off, b%maxnp, b%nshell, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                           xg%npts, xg%r, xg%b_off, xg%b_shoff, size(xg%b_sh), xg%b_sh, xg%b_shao, &
                           c_off, npp, xg%nbatch, nchi, chi, gchi)
         call xc_gather(npairs, nb, ord, pb_off, nspin, b%nao, dmat, xg%b_aooff, size(xg%b_ao), xg%b_ao, &
                        nlp, xg%nbatch, dloc)
         ! X_b = chi_b D_b: one strided-batched GEMM per run of equal shape.
         k0 = 1
         do while (k0 <= nb)
            ib = ord(k0)
            mp = npp(ib)
            ml = nlp(ib)
            k1 = k0
            do while (k1 < nb)
               if (npp(ord(k1 + 1)) /= mp .or. nlp(ord(k1 + 1)) /= ml) exit
               k1 = k1 + 1
            end do
            cnt = k1 - k0 + 1
            if (ml > 0) then
               sa = int(mp, int64)*int(ml, int64)
               sb = int(ml, int64)*int(ml, int64)
               do s = 1, nspin
                  call lp%gemm_strided("N", "N", mp, ml, ml, 1.0_dp, chi(c_off(ib) + 1), mp, sa, &
                                       dloc(pr_off(ib) + 1, s), ml, sb, 0.0_dp, xt(c_off(ib) + 1, s), mp, sa, cnt)
                  if (mgga) then
                     ! Y_x = (d chi/dx)_b D_b, for tau.
                     do x = 1, 3
                        call lp%gemm_strided("N", "N", mp, ml, ml, 1.0_dp, gchi(c_off(ib) + 1, x), mp, sa, &
                                             dloc(pr_off(ib) + 1, s), ml, sb, 0.0_dp, yt(c_off(ib) + 1, x, s), &
                                             mp, sa, cnt)
                     end do
                  end if
               end do
            end if
            k0 = k1 + 1
         end do
         if (mgga) then
            call xc_points_mgga(np, p0, nspin, xg%npts, xg%w, xg%batch_of, xg%b_off, xg%b_aooff, c_off, npp, &
                                active, xg%nbatch, func%nk, kid, coef, nchi, chi, gchi, xt, zt, yt, eg, ng)
         else
            call xc_points(np, p0, nspin, xg%npts, xg%w, xg%batch_of, xg%b_off, xg%b_aooff, c_off, npp, active, &
                           xg%nbatch, func%nk, kid, coef, nchi, chi, gchi, xt, zt, eg, ng)
         end if
         t1 = wall()

         !$acc update self(eg, ng)
         exc = exc + sum(eg)
         nelec = nelec + sum(ng)

         ! V_b = chi_b^T Z_b, same grouping, then into V with its transpose.
         k0 = 1
         do while (k0 <= nb)
            ib = ord(k0)
            mp = npp(ib)
            ml = nlp(ib)
            k1 = k0
            do while (k1 < nb)
               if (npp(ord(k1 + 1)) /= mp .or. nlp(ord(k1 + 1)) /= ml) exit
               k1 = k1 + 1
            end do
            cnt = k1 - k0 + 1
            if (ml > 0) then
               sa = int(mp, int64)*int(ml, int64)
               sc = int(ml, int64)*int(ml, int64)
               do s = 1, nspin
                  call lp%gemm_strided("T", "N", ml, ml, mp, 1.0_dp, chi(c_off(ib) + 1), mp, sa, &
                                       zt(c_off(ib) + 1, s), mp, sa, 0.0_dp, vloc(pr_off(ib) + 1, s), ml, sc, cnt)
                  if (mgga) then
                     ! V_b += sum_x (d chi/dx)_b^T W_x, W_x in Y_x's storage.
                     do x = 1, 3
                        call lp%gemm_strided("T", "N", ml, ml, mp, 1.0_dp, gchi(c_off(ib) + 1, x), mp, sa, &
                                             yt(c_off(ib) + 1, x, s), mp, sa, 1.0_dp, vloc(pr_off(ib) + 1, s), &
                                             ml, sc, cnt)
                     end do
                  end if
               end do
            end if
            k0 = k1 + 1
         end do
         call xc_scatter(npairs, nb, ord, pb_off, nspin, b%nao, vxc, xg%b_aooff, size(xg%b_ao), xg%b_ao, &
                         nlp, xg%nbatch, vloc)
         t2 = wall()
         t_points = t_points + (t1 - t0)
         t_pairs = t_pairs + (t2 - t1)

         !$acc exit data delete(chi, gchi, xt, zt, dloc, vloc, eg, ng, c_off, pr_off, ord, pb_off, wk_off)
         deallocate (chi, gchi, xt, zt, dloc, vloc, eg, ng, ord, pb_off, wk_off)
         if (mgga) then
            !$acc exit data delete(yt)
            deallocate (yt)
         end if
         b0 = b1 + 1
      end do
      if (.not. present(la)) call la_own%release()
      !$acc end data

      ! Across ranks: sum, then rank 0's copy everywhere. A resident caller's
      ! V lives on the device and the data region above did not bring it
      ! back, so fetch it, reduce on the host, and put the result back.
      if (split) then
         !$acc update self(vxc) if_present
         call sum_ranks(comm, vxc, b%nao*b%nao*nspin)
         !$acc update device(vxc) if_present
         es = [exc, nelec]
         call sum_ranks(comm, es, 2)
         exc = es(1)
         nelec = es(2)
      end if
   end subroutine xc_drive

   !
   ! The CPKS XC kernel response, restricted: f_xc[D_U]_uv = integral
   ! phi_u phi_v (v2rho2 rho_U + v2rhosigma dsigma_U) for LDA, plus the GGA
   ! sigma-response terms through grad rho (ground state) and grad rho_U
   ! (trial) -- xc_resp_point_body below has the full per-point formula
   ! (trc_cpks's own header only names f_xc[D_U] as a term, it does not
   ! derive it). Same grid
   ! batches, same collocate-then-GEMM pipeline as `trc_xc_rks` (xc_drive's
   ! own header): the one difference is gathering TWO densities per batch
   ! (`dmat`, fixed for a whole CPKS solve, and `dmat_u`, rebuilt every CG
   ! iteration) instead of one, and a point kernel that evaluates the
   ! second derivatives instead of the first. No density screening and no
   ! rank splitting here -- unlike the energy path, every CPKS iteration
   ! needs the SAME linear operator, and the trial density is a CG
   ! iterate the solve drives towards zero as it converges, so screening
   ! on ITS magnitude would tighten from one matvec to the next and break
   ! that assumption -- see `eri_fock`'s own `density_screen` doc
   ! (trc_eri.F90) for the identical reason on the ERI side, where the
   ! coupled-perturbed equations have the same requirement. A later
   ! optimisation can screen on the GROUND density alone (fixed across
   ! the CG loop, unlike `dmat_u`).
   !
   ! `dmat` and `dmat_u` are both already device-resident (same contract as
   ! trc_cphf's `cmo`/`D_U`): the caller copies `dmat` in once, before the
   ! CG loop starts, and keeps `dmat_u` resident from one matvec to the
   ! next. `v_resp` is written in full (not accumulated into); the caller
   ! adds it into its own Fock-response matrix.
   !
   subroutine trc_xc_fxc_rks(b, xg, func, dmat, dmat_u, v_resp, budget, la)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      real(dp), intent(in) :: dmat(b%nao, b%nao), dmat_u(b%nao, b%nao)
      real(dp), intent(out) :: v_resp(b%nao, b%nao)
      integer(int64), intent(in), optional :: budget
      type(trc_linalg_t), intent(inout), optional, target :: la

      integer :: nao
      integer(int64) :: need, nchi, sa, sb, sc, cap
      integer :: b0, b1, ib, nloc, npb, np, p0, npairs, nb, k0, k1, cnt, mp, ml, s
      integer, allocatable :: c_off(:), pr_off(:), npp(:), nlp(:), ord(:), pb_off(:), wk_off(:)
      integer :: nwork
      real(dp), allocatable :: chi(:), gchi(:, :), xt(:, :), zt(:, :), dloc(:, :), vloc(:), dpack(:, :, :)
      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      integer :: kid(XC_MAX_KERNELS)
      real(dp) :: coef(XC_MAX_KERNELS)
      logical :: owns_la

      nao = b%nao
      kid = func%kid
      coef = func%coef
      v_resp = 0.0_dp
      if (xg%npts == 0) return
      cap = 8388608_int64
      if (present(budget)) cap = max(budget, 1_int64)

      allocate (c_off(xg%nbatch + 1), pr_off(xg%nbatch + 1))
      allocate (npp(xg%nbatch), nlp(xg%nbatch))
      allocate (dpack(nao, nao, 2))

      owns_la = .not. present(la)
      if (present(la)) then
         lp => la
      else
         lp => la_own
         call la_own%init(1)
      end if

      !$acc data copyin(dmat, dmat_u, kid, coef) copy(v_resp) create(npp, nlp, dpack)
      !$acc update device(v_resp)
      call xc_pack2(nao, dmat, dmat_u, dpack)
      do ib = 1, xg%nbatch
         nloc = xg%b_aooff(ib + 1) - xg%b_aooff(ib)
         npb = xg%b_off(ib + 1) - xg%b_off(ib)
         npp(ib) = ((npb + XC_NP_PAD - 1)/XC_NP_PAD)*XC_NP_PAD
         nlp(ib) = ((nloc + XC_NL_PAD - 1)/XC_NL_PAD)*XC_NL_PAD
      end do
      !$acc update device(npp, nlp)

      b0 = 1
      do while (b0 <= xg%nbatch)
         need = 0
         b1 = b0 - 1
         do ib = b0, xg%nbatch
            if (need + int(npp(ib), int64)*int(nlp(ib), int64) > cap .and. ib > b0) exit
            need = need + int(npp(ib), int64)*int(nlp(ib), int64)
            b1 = ib
         end do
         nchi = max(need, 1_int64)
         p0 = xg%b_off(b0)
         np = xg%b_off(b1 + 1) - p0
         nb = b1 - b0 + 1
         allocate (ord(nb), pb_off(nb + 1), wk_off(nb + 1))
         do k0 = 1, nb
            ord(k0) = b0 + k0 - 1
         end do
         call sort_by_shape(nb, npp, nlp, ord)
         need = 0
         npairs = 0
         nwork = 0
         pb_off(1) = 0
         wk_off(1) = 0
         do k0 = 1, nb
            ib = ord(k0)
            c_off(ib) = int(need)
            pr_off(ib) = npairs
            need = need + int(npp(ib), int64)*int(nlp(ib), int64)
            npairs = npairs + nlp(ib)*nlp(ib)
            pb_off(k0 + 1) = npairs
            nwork = nwork + (xg%b_shoff(ib + 1) - xg%b_shoff(ib))*(xg%b_off(ib + 1) - xg%b_off(ib))
            wk_off(k0 + 1) = nwork
         end do

         allocate (chi(nchi), gchi(nchi, 3), xt(nchi, 2), zt(nchi, 1))
         allocate (dloc(max(npairs, 1), 2), vloc(max(npairs, 1)))
         !$acc enter data create(chi, gchi, xt, zt, dloc, vloc) copyin(c_off, pr_off, ord, pb_off, wk_off)

         call xc_zero(nchi, 1, chi, zt)
         call xc_collocate(nwork, nb, ord, wk_off, b%maxnp, b%nshell, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                           xg%npts, xg%r, xg%b_off, xg%b_shoff, size(xg%b_sh), xg%b_sh, xg%b_shao, &
                           c_off, npp, xg%nbatch, nchi, chi, gchi)
         call xc_gather(npairs, nb, ord, pb_off, 2, nao, dpack, xg%b_aooff, size(xg%b_ao), xg%b_ao, &
                        nlp, xg%nbatch, dloc)
         k0 = 1
         do while (k0 <= nb)
            ib = ord(k0)
            mp = npp(ib)
            ml = nlp(ib)
            k1 = k0
            do while (k1 < nb)
               if (npp(ord(k1 + 1)) /= mp .or. nlp(ord(k1 + 1)) /= ml) exit
               k1 = k1 + 1
            end do
            cnt = k1 - k0 + 1
            if (ml > 0) then
               sa = int(mp, int64)*int(ml, int64)
               sb = int(ml, int64)*int(ml, int64)
               do s = 1, 2
                  call lp%gemm_strided("N", "N", mp, ml, ml, 1.0_dp, chi(c_off(ib) + 1), mp, sa, &
                                       dloc(pr_off(ib) + 1, s), ml, sb, 0.0_dp, xt(c_off(ib) + 1, s), mp, sa, cnt)
               end do
            end if
            k0 = k1 + 1
         end do
         call xc_resp_points(np, p0, xg%npts, xg%w, xg%batch_of, xg%b_off, xg%b_aooff, c_off, npp, &
                             xg%nbatch, func%nk, kid, coef, nchi, chi, gchi, xt, zt)

         k0 = 1
         do while (k0 <= nb)
            ib = ord(k0)
            mp = npp(ib)
            ml = nlp(ib)
            k1 = k0
            do while (k1 < nb)
               if (npp(ord(k1 + 1)) /= mp .or. nlp(ord(k1 + 1)) /= ml) exit
               k1 = k1 + 1
            end do
            cnt = k1 - k0 + 1
            if (ml > 0) then
               sa = int(mp, int64)*int(ml, int64)
               sc = int(ml, int64)*int(ml, int64)
               call lp%gemm_strided("T", "N", ml, ml, mp, 1.0_dp, chi(c_off(ib) + 1), mp, sa, &
                                    zt(c_off(ib) + 1, 1), mp, sa, 0.0_dp, vloc(pr_off(ib) + 1), ml, sc, cnt)
            end if
            k0 = k1 + 1
         end do
         call xc_resp_scatter(npairs, nb, ord, pb_off, nao, v_resp, xg%b_aooff, size(xg%b_ao), xg%b_ao, &
                              nlp, xg%nbatch, vloc)

         !$acc exit data delete(chi, gchi, xt, zt, dloc, vloc, c_off, pr_off, ord, pb_off, wk_off)
         deallocate (chi, gchi, xt, zt, dloc, vloc, ord, pb_off, wk_off)
         b0 = b1 + 1
      end do
      if (owns_la) call la_own%release()
      !$acc end data
   end subroutine trc_xc_fxc_rks

   ! dpack(:,:,1) = a, dpack(:,:,2) = bb -- packs the ground and trial
   ! densities into one array so a single `xc_gather` call (nspin = 2)
   ! fetches both per-batch tiles, rather than two differently-shaped calls.
   subroutine xc_pack2(nao, a, bb, dpack)
      integer, intent(in) :: nao
      real(dp), intent(in) :: a(nao, nao), bb(nao, nao)
      real(dp), intent(inout) :: dpack(nao, nao, 2)
      integer :: i, j
      do concurrent(j=1:nao, i=1:nao)
         dpack(i, j, 1) = a(i, j)
         dpack(i, j, 2) = bb(i, j)
      end do
   end subroutine xc_pack2

   !
   ! Kernel 3 for the response pass: one thread per point, reading the
   ! ground state's row of X (slot 1) and the trial density's (slot 2),
   ! writing the single response tile `zt(:,1)`. See trc_xc_fxc_rks's own
   ! header for the formula.
   !
   subroutine xc_resp_points(np, p0, npts, w, batch_of, b_off, b_aooff, c_off, npp, nbatch, &
                             nk, kid, coef, nchi, chi, gchi, xt, zt)
      integer, intent(in) :: np, p0, npts, nbatch, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3), xt(nchi, 2)
      real(dp), intent(inout) :: zt(nchi, 1)
      integer :: g

      do concurrent(g=1:np)
         call xc_resp_point_body(p0 + g - 1, npts, w, batch_of, b_off, b_aooff, c_off, npp, nbatch, &
                                 nk, kid, coef, nchi, chi, gchi, xt, zt)
      end do
   end subroutine xc_resp_points

   pure subroutine xc_resp_point_body(gg, npts, w, batch_of, b_off, b_aooff, c_off, npp, nbatch, &
                                      nk, kid, coef, nchi, chi, gchi, xt, zt)
      !$acc routine seq
      integer, intent(in) :: gg, npts, nbatch, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3), xt(nchi, 2)
      real(dp), intent(inout) :: zt(nchi, 1)

      integer :: ib, nloc, ld, off, u, i
      real(dp) :: wg, x, xu, sigma, dsig_u, vsigma, v2rho2, v2rhosigma, v2sigma2
      real(dp) :: rho, rho_u, zg_r
      real(dp) :: grho(3), grho_u(3), zv_r(3)

      ib = batch_of(gg)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = npp(ib)
      off = c_off(ib) + (gg - b_off(ib)) + 1
      wg = w(gg)
      rho = 0.0_dp
      grho = 0.0_dp
      rho_u = 0.0_dp
      grho_u = 0.0_dp
      do u = 1, nloc
         i = off + (u - 1)*ld
         x = xt(i, 1)
         xu = xt(i, 2)
         rho = rho + x*chi(i)
         rho_u = rho_u + xu*chi(i)
         grho(1) = grho(1) + 2.0_dp*x*gchi(i, 1)
         grho(2) = grho(2) + 2.0_dp*x*gchi(i, 2)
         grho(3) = grho(3) + 2.0_dp*x*gchi(i, 3)
         grho_u(1) = grho_u(1) + 2.0_dp*xu*gchi(i, 1)
         grho_u(2) = grho_u(2) + 2.0_dp*xu*gchi(i, 2)
         grho_u(3) = grho_u(3) + 2.0_dp*xu*gchi(i, 3)
      end do
      sigma = grho(1)**2 + grho(2)**2 + grho(3)**2
      call xc_eval_point_fxc(nk, kid, coef, rho, sigma, vsigma, v2rho2, v2rhosigma, v2sigma2)
      ! dsigma_U = 2 grad(rho).grad(rho_U) -- sigma's own response to the
      ! trial density, through the UNPERTURBED grad(rho) (ground state).
      dsig_u = 2.0_dp*(grho(1)*grho_u(1) + grho(2)*grho_u(2) + grho(3)*grho_u(3))
      zg_r = wg*(v2rho2*rho_u + v2rhosigma*dsig_u)
      ! The GGA term has two pieces: the second-derivative response
      ! multiplying the UNPERTURBED grad(rho) (vsigma's own value changes),
      ! plus vsigma (unperturbed) multiplying the TRIAL gradient directly
      ! (the explicit grad(rho) factor in the energy path's own vxc becomes
      ! grad(rho_U) here) -- trc_xc_fxc_rks's own header.
      zv_r(1) = wg*2.0_dp*(v2rhosigma*rho_u + v2sigma2*dsig_u)*grho(1) + wg*2.0_dp*vsigma*grho_u(1)
      zv_r(2) = wg*2.0_dp*(v2rhosigma*rho_u + v2sigma2*dsig_u)*grho(2) + wg*2.0_dp*vsigma*grho_u(2)
      zv_r(3) = wg*2.0_dp*(v2rhosigma*rho_u + v2sigma2*dsig_u)*grho(3) + wg*2.0_dp*vsigma*grho_u(3)
      do u = 1, nloc
         i = off + (u - 1)*ld
         zt(i, 1) = 0.5_dp*zg_r*chi(i) + zv_r(1)*gchi(i, 1) + zv_r(2)*gchi(i, 2) + zv_r(3)*gchi(i, 3)
      end do
   end subroutine xc_resp_point_body

   !
   ! Kernel 4 for the response pass: V_resp(u,v) += V_b(u,v) + V_b(v,u),
   ! one output channel -- xc_scatter's own twin, duplicated rather than
   ! made generic over nspin because its `vxc` here is genuinely rank 2
   ! (trc_cpks has one response matrix, not a per-spin pair).
   !
   subroutine xc_resp_scatter(npairs, nb, ord, pb_off, nao, vxc, b_aooff, nbao, b_ao, nlp, nbatch, vloc)
      integer, intent(in) :: npairs, nb, nao, nbao, nbatch
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(inout) :: vxc(nao, nao)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(in) :: vloc(max(npairs, 1))
      integer :: t, au, av
      real(dp) :: acc_a

#ifdef _OPENACC
      do concurrent(t=1:npairs) local(au, av, acc_a)
         call xc_resp_scatter_body(t, nb, ord, pb_off, b_aooff, nbao, b_ao, nlp, nbatch, npairs, vloc, &
                                   au, av, acc_a)
         if (au > 0) then
            !$acc atomic update
            vxc(au, av) = vxc(au, av) + acc_a
         end if
      end do
#else
      do t = 1, npairs
         call xc_resp_scatter_body(t, nb, ord, pb_off, b_aooff, nbao, b_ao, nlp, nbatch, npairs, vloc, &
                                   au, av, acc_a)
         if (au > 0) vxc(au, av) = vxc(au, av) + acc_a
      end do
#endif
   end subroutine xc_resp_scatter

   pure subroutine xc_resp_scatter_body(t, nb, ord, pb_off, b_aooff, nbao, b_ao, nlp, nbatch, npairs, vloc, &
                                        au, av, acc_a)
      !$acc routine seq
      integer, intent(in) :: t, nb, nbao, nbatch, npairs
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(in) :: vloc(npairs)
      integer, intent(out) :: au, av
      real(dp), intent(out) :: acc_a
      integer :: k, ib, tt, nloc, ld, u, v, tvu

      call pair_batch(t, nb, pb_off, k, tt)
      ib = ord(k)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = nlp(ib)
      u = mod(tt, ld) + 1
      v = tt/ld + 1
      acc_a = 0.0_dp
      au = 0
      av = 0
      if (u > nloc .or. v > nloc) return
      tvu = pb_off(k) + (u - 1)*ld + v
      acc_a = vloc(t) + vloc(tvu)
      au = b_ao(b_aooff(ib) + u - 1)
      av = b_ao(b_aooff(ib) + v - 1)
   end subroutine xc_resp_scatter_body

   ! Sum a contiguous array over the ranks, then rank 0's copy everywhere.
   subroutine sum_ranks(comm, g, n)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call allreduce(comm, g, op=MPI_SUM)
      call bcast(comm, g, n, 0)
   end subroutine sum_ranks

   ! Batches a rank owns under the round-robin: the k with mod(k-1,nranks)=rank.
   pure integer function n_owned(nbatch, rank, nranks) result(n)
      integer, intent(in) :: nbatch, rank, nranks
      n = nbatch/nranks
      if (mod(nbatch, nranks) > rank) n = n + 1
   end function n_owned

   !
   ! Density screen: one thread per batch bounds |rho| on it through the
   ! stored per-function maxima and clears its flag if the bound is under
   ! the tolerance. A batch already inactive (another rank's) is skipped.
   !
   subroutine xc_screen(nbatch, nspin, nao, dmat, b_aooff, nbao, b_ao, b_amax, rho_tol, active)
      integer, intent(in) :: nbatch, nspin, nao, nbao
      real(dp), intent(in) :: rho_tol
      real(dp), intent(in) :: dmat(nao, nao, nspin), b_amax(nbao)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao)
      integer, intent(inout) :: active(nbatch)
      integer :: ib
      do concurrent(ib=1:nbatch)
         call xc_screen_body(ib, nbatch, nspin, nao, dmat, b_aooff, nbao, b_ao, b_amax, rho_tol, active(ib))
      end do
   end subroutine xc_screen

   pure subroutine xc_screen_body(ib, nbatch, nspin, nao, dmat, b_aooff, nbao, b_ao, b_amax, rho_tol, flag)
      !$acc routine seq
      integer, intent(in) :: ib, nbatch, nspin, nao, nbao
      real(dp), intent(in) :: rho_tol
      real(dp), intent(in) :: dmat(nao, nao, nspin), b_amax(nbao)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao)
      integer, intent(inout) :: flag
      integer :: s, u, v, au, av
      real(dp) :: bound, su
      if (flag == 0) return
      bound = 0.0_dp
      do s = 1, nspin
         do u = b_aooff(ib), b_aooff(ib + 1) - 1
            au = b_ao(u)
            su = 0.0_dp
            do v = b_aooff(ib), b_aooff(ib + 1) - 1
               av = b_ao(v)
               su = su + abs(dmat(av, au, s))*b_amax(v)
            end do
            bound = bound + su*b_amax(u)
         end do
      end do
      if (bound <= rho_tol) flag = 0
   end subroutine xc_screen_body

   ! Insertion sort of `ord` by (npp, nlp) ascending; a chunk's worth of batches.
   subroutine sort_by_shape(n, npp, nlp, ord)
      integer, intent(in) :: n
      integer, intent(in) :: npp(:), nlp(:)
      integer, intent(inout) :: ord(n)
      integer :: i, j, t
      do i = 2, n
         t = ord(i)
         j = i - 1
         do while (j >= 1)
            if (.not. later(ord(j), t)) exit
            ord(j + 1) = ord(j)
            j = j - 1
         end do
         ord(j + 1) = t
      end do
   contains
      pure logical function later(a, b) result(res)
         integer, intent(in) :: a, b
         res = npp(a) > npp(b) .or. (npp(a) == npp(b) .and. nlp(a) > nlp(b))
      end function later
   end subroutine sort_by_shape

   ! Insertion sort of `order` by descending `cost`; a few thousand batches.
   subroutine sort_desc(n, cost, order)
      integer, intent(in) :: n
      real(dp), intent(in) :: cost(n)
      integer, intent(inout) :: order(n)
      integer :: i, j, t
      do i = 2, n
         t = order(i)
         j = i - 1
         do while (j >= 1)
            if (cost(order(j)) >= cost(t)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = t
      end do
   end subroutine sort_desc

   function wall() result(t)
      real(dp) :: t
      integer(int64) :: cc, rate
      ! TODO(fix): system_clock -> pic_timer timer_type
      call system_clock(cc, rate)
      t = real(cc, dp)/real(rate, dp)
   end function wall

   !
   ! Which stored batch (1..nb, storage order) owns pair t: the last k with
   ! pb_off(k) < t.
   !
   pure subroutine pair_batch(t, nb, pb_off, k, tt)
      !$acc routine seq
      integer, intent(in) :: t, nb
      integer, intent(in) :: pb_off(nb + 1)
      integer, intent(out) :: k, tt
      integer :: lo, hi, mid
      lo = 1
      hi = nb
      do while (lo < hi)
         mid = (lo + hi + 1)/2
         if (pb_off(mid) < t) then
            lo = mid
         else
            hi = mid - 1
         end if
      end do
      k = lo
      tt = t - pb_off(k) - 1
   end subroutine pair_batch

   ! Zero the tiles whose padding the GEMMs read.
   subroutine xc_zero(nchi, nspin, chi, zt)
      integer(int64), intent(in) :: nchi
      integer, intent(in) :: nspin
      real(dp), intent(out) :: chi(nchi), zt(nchi, nspin)
      integer(int64) :: i
      integer :: s
      do concurrent(i=1:nchi)
         chi(i) = 0.0_dp
      end do
      do concurrent(i=1:nchi, s=1:nspin)
         zt(i, s) = 0.0_dp
      end do
   end subroutine xc_zero

   ! Zero the gradient tiles, whose padding only a meta-GGA's GEMMs read.
   ! One flattened index rather than (i, x): the outer extent is small.
   subroutine xc_zero_grad(nchi, gchi)
      integer(int64), intent(in) :: nchi
      real(dp), intent(out) :: gchi(nchi*3)
      integer(int64) :: i
      do concurrent(i=1:3*nchi)
         gchi(i) = 0.0_dp
      end do
   end subroutine xc_zero_grad

   !
   ! Kernel 1: collocation -- one thread per (shell, point), Algorithm 2 of
   ! Stocks and Barca. Work item t maps to a stored batch, then a shell of
   ! its local list and a point of its box, point fastest: the threads of a
   ! warp evaluate one shell -- one l, one contraction, no divergence -- at
   ! consecutive points, and write consecutive tile addresses. Tile layout
   ! is (points, functions), point fastest, the (npp x nlp) column-major
   ! matrix the GEMMs read.
   !
   subroutine xc_collocate(nwork, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                           b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, nchi, chi, gchi)
      integer, intent(in) :: nwork, nb, maxnp, nshell, npts, nbsh, nbatch
      integer, intent(in) :: ord(nb), wk_off(nb + 1)
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: c_off(nbatch + 1), npp(nbatch)
      integer(int64), intent(in) :: nchi
      real(dp), intent(inout) :: chi(nchi), gchi(nchi, 3)
      integer :: t

      do concurrent(t=1:nwork)
         call xc_collocate_body(t, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, nchi, chi, gchi)
      end do
   end subroutine xc_collocate

   pure subroutine xc_collocate_body(t, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                     b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, nchi, chi, gchi)
      !$acc routine seq
      integer, intent(in) :: t, nb, maxnp, nshell, npts, nbsh, nbatch
      integer, intent(in) :: ord(nb), wk_off(nb + 1)
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: c_off(nbatch + 1), npp(nbatch)
      integer(int64), intent(in) :: nchi
      real(dp), intent(inout) :: chi(nchi), gchi(nchi, 3)
      integer :: k, ib, tt, npb, ks, gl, gg, ish, l, nc, m, ld, off, iao
      real(dp) :: d(3), cs(NCART_MAX), gs(3, NCART_MAX)

      call pair_batch(t, nb, wk_off, k, tt)
      ib = ord(k)
      npb = b_off(ib + 1) - b_off(ib)
      ks = tt/npb
      gl = tt - ks*npb
      gg = b_off(ib) + gl
      ish = b_sh(b_shoff(ib) + ks)
      iao = b_shao(b_shoff(ib) + ks)
      l = sh_l(ish)
      d(1) = r(1, gg) - sh_r(1, ish)
      d(2) = r(2, gg) - sh_r(2, ish)
      d(3) = r(3, gg) - sh_r(3, ish)
      call shell_collocate(l, sh_np(ish), sh_e(1:sh_np(ish), ish), sh_c(1:sh_np(ish), ish), d, cs, gs)
      nc = (l + 1)*(l + 2)/2
      ld = npp(ib)
      off = c_off(ib) + gl + 1
      do m = 1, nc
         chi(off + (iao + m - 1)*ld) = cs(m)
         gchi(off + (iao + m - 1)*ld, 1) = gs(1, m)
         gchi(off + (iao + m - 1)*ld, 2) = gs(2, m)
         gchi(off + (iao + m - 1)*ld, 3) = gs(3, m)
      end do
   end subroutine xc_collocate_body

   !
   ! Kernel 2: the density slice of every batch -- one thread per stored
   ! (b, u, v), padding included and written as zero.
   !
   subroutine xc_gather(npairs, nb, ord, pb_off, nspin, nao, dmat, b_aooff, nbao, b_ao, nlp, nbatch, dloc)
      integer, intent(in) :: npairs, nb, nspin, nao, nbao, nbatch
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(in) :: dmat(nao, nao, nspin)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(out) :: dloc(max(npairs, 1), nspin)
      integer :: t

      do concurrent(t=1:npairs)
         call xc_gather_body(t, nb, ord, pb_off, nspin, nao, dmat, b_aooff, nbao, b_ao, nlp, nbatch, npairs, dloc)
      end do
   end subroutine xc_gather

   pure subroutine xc_gather_body(t, nb, ord, pb_off, nspin, nao, dmat, b_aooff, nbao, b_ao, nlp, nbatch, npairs, dloc)
      !$acc routine seq
      integer, intent(in) :: t, nb, nspin, nao, nbao, nbatch, npairs
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(in) :: dmat(nao, nao, nspin)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(inout) :: dloc(npairs, nspin)
      integer :: k, ib, tt, nloc, ld, u, v, s

      call pair_batch(t, nb, pb_off, k, tt)
      ib = ord(k)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = nlp(ib)
      u = mod(tt, ld) + 1
      v = tt/ld + 1
      if (u > nloc .or. v > nloc) then
         do s = 1, nspin
            dloc(t, s) = 0.0_dp
         end do
      else
         do s = 1, nspin
            dloc(t, s) = dmat(b_ao(b_aooff(ib) + u - 1), b_ao(b_aooff(ib) + v - 1), s)
         end do
      end if
   end subroutine xc_gather_body

   !
   ! Kernel 3: the functional -- one thread per point. Reads its row of
   ! chi, grad chi and X = chi D, writes its row of Z and its energy term.
   !
   subroutine xc_points(np, p0, nspin, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, nbatch, &
                        nk, kid, coef, nchi, chi, gchi, xt, zt, eg, ng)
      integer, intent(in) :: np, p0, nspin, npts, nbatch, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), active(nbatch), kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3), xt(nchi, nspin)
      real(dp), intent(inout) :: zt(nchi, nspin)
      real(dp), intent(out) :: eg(np), ng(np)
      integer :: g

      do concurrent(g=1:np)
         call xc_point_body(p0 + g - 1, nspin, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, nbatch, &
                            nk, kid, coef, nchi, chi, gchi, xt, zt, eg(g), ng(g))
      end do
   end subroutine xc_points

   pure subroutine xc_point_body(gg, nspin, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, nbatch, &
                                 nk, kid, coef, nchi, chi, gchi, xt, zt, eg, ng)
      !$acc routine seq
      integer, intent(in) :: gg, nspin, npts, nbatch, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), active(nbatch), kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3), xt(nchi, nspin)
      real(dp), intent(inout) :: zt(nchi, nspin)
      real(dp), intent(out) :: eg, ng

      integer :: ib, nloc, ld, off, u, s, i
      real(dp) :: wg, x, sigma, eps, vrho, vsigma
      real(dp) :: rho(2), grho(3, 2), zg(2), zv(3, 2)
      real(dp) :: s_aa, s_ab, s_bb, va, vb, v_aa, v_ab, v_bb

      ib = batch_of(gg)
      eg = 0.0_dp
      ng = 0.0_dp
      if (active(ib) == 0) return
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = npp(ib)
      off = c_off(ib) + (gg - b_off(ib)) + 1
      wg = w(gg)
      ! rho = sum_u X_u chi_u, grad rho = 2 sum_u X_u grad chi_u (D symmetric)
      rho = 0.0_dp
      grho = 0.0_dp
      do s = 1, nspin
         do u = 1, nloc
            i = off + (u - 1)*ld
            x = xt(i, s)
            rho(s) = rho(s) + x*chi(i)
            grho(1, s) = grho(1, s) + 2.0_dp*x*gchi(i, 1)
            grho(2, s) = grho(2, s) + 2.0_dp*x*gchi(i, 2)
            grho(3, s) = grho(3, s) + 2.0_dp*x*gchi(i, 3)
         end do
      end do
      if (nspin == 1) then
         sigma = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         call xc_eval_point(nk, kid, coef, rho(1), sigma, eps, vrho, vsigma)
         eg = wg*rho(1)*eps
         ng = wg*rho(1)
         zg(1) = wg*vrho
         zv(1, 1) = 2.0_dp*wg*vsigma*grho(1, 1)
         zv(2, 1) = 2.0_dp*wg*vsigma*grho(2, 1)
         zv(3, 1) = 2.0_dp*wg*vsigma*grho(3, 1)
      else
         s_aa = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         s_bb = grho(1, 2)**2 + grho(2, 2)**2 + grho(3, 2)**2
         s_ab = grho(1, 1)*grho(1, 2) + grho(2, 1)*grho(2, 2) + grho(3, 1)*grho(3, 2)
         call xc_eval_point_polar(nk, kid, coef, rho(1), rho(2), s_aa, s_ab, s_bb, &
                                  eps, va, vb, v_aa, v_ab, v_bb)
         eg = wg*(rho(1) + rho(2))*eps
         ng = wg*(rho(1) + rho(2))
         zg(1) = wg*va
         zg(2) = wg*vb
         zv(1, 1) = wg*(2.0_dp*v_aa*grho(1, 1) + v_ab*grho(1, 2))
         zv(2, 1) = wg*(2.0_dp*v_aa*grho(2, 1) + v_ab*grho(2, 2))
         zv(3, 1) = wg*(2.0_dp*v_aa*grho(3, 1) + v_ab*grho(3, 2))
         zv(1, 2) = wg*(2.0_dp*v_bb*grho(1, 2) + v_ab*grho(1, 1))
         zv(2, 2) = wg*(2.0_dp*v_bb*grho(2, 2) + v_ab*grho(2, 1))
         zv(3, 2) = wg*(2.0_dp*v_bb*grho(3, 2) + v_ab*grho(3, 1))
      end if
      ! Z_u = zg/2 chi_u + zv . grad chi_u: chi^T Z plus its transpose is
      ! sum_g zg chi_u chi_v + zv . (grad chi_u chi_v + chi_u grad chi_v).
      do s = 1, nspin
         do u = 1, nloc
            i = off + (u - 1)*ld
            zt(i, s) = 0.5_dp*zg(s)*chi(i) + zv(1, s)*gchi(i, 1) + zv(2, s)*gchi(i, 2) + zv(3, s)*gchi(i, 3)
         end do
      end do
   end subroutine xc_point_body

   !
   ! Kernel 3 for a meta-GGA: xc_point_body plus tau. Reads its row of Y_x
   ! = (d chi/dx) D as well, and overwrites it with W_x = w v_tau/4 d chi/dx
   ! (the module header). A separate kernel rather than a branch in
   ! xc_point_body, so the LDA/GGA one compiles to exactly what it did.
   !
   subroutine xc_points_mgga(np, p0, nspin, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, nbatch, &
                             nk, kid, coef, nchi, chi, gchi, xt, zt, yt, eg, ng)
      integer, intent(in) :: np, p0, nspin, npts, nbatch, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), active(nbatch), kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3), xt(nchi, nspin)
      real(dp), intent(inout) :: zt(nchi, nspin), yt(nchi, 3, nspin)
      real(dp), intent(out) :: eg(np), ng(np)
      integer :: g

      do concurrent(g=1:np)
         call xc_point_mgga_body(p0 + g - 1, nspin, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, nbatch, &
                                 nk, kid, coef, nchi, chi, gchi, xt, zt, yt, eg(g), ng(g))
      end do
   end subroutine xc_points_mgga

   pure subroutine xc_point_mgga_body(gg, nspin, npts, w, batch_of, b_off, b_aooff, c_off, npp, active, nbatch, &
                                      nk, kid, coef, nchi, chi, gchi, xt, zt, yt, eg, ng)
      !$acc routine seq
      integer, intent(in) :: gg, nspin, npts, nbatch, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_aooff(nbatch + 1), c_off(nbatch + 1)
      integer, intent(in) :: npp(nbatch), active(nbatch), kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3), xt(nchi, nspin)
      real(dp), intent(inout) :: zt(nchi, nspin), yt(nchi, 3, nspin)
      real(dp), intent(out) :: eg, ng

      integer :: ib, nloc, ld, off, u, s, i
      real(dp) :: wg, x, sigma, eps, vrho, vsigma
      real(dp) :: rho(2), grho(3, 2), tau(2), vtau(2), zg(2), zv(3, 2), zw(2)
      real(dp) :: s_aa, s_ab, s_bb, va, vb, v_aa, v_ab, v_bb

      ib = batch_of(gg)
      eg = 0.0_dp
      ng = 0.0_dp
      if (active(ib) == 0) return
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = npp(ib)
      off = c_off(ib) + (gg - b_off(ib)) + 1
      wg = w(gg)
      ! rho and grad rho as xc_point_body; tau = 1/2 sum_u Y_u . grad chi_u
      rho = 0.0_dp
      grho = 0.0_dp
      tau = 0.0_dp
      do s = 1, nspin
         do u = 1, nloc
            i = off + (u - 1)*ld
            x = xt(i, s)
            rho(s) = rho(s) + x*chi(i)
            grho(1, s) = grho(1, s) + 2.0_dp*x*gchi(i, 1)
            grho(2, s) = grho(2, s) + 2.0_dp*x*gchi(i, 2)
            grho(3, s) = grho(3, s) + 2.0_dp*x*gchi(i, 3)
            tau(s) = tau(s) + yt(i, 1, s)*gchi(i, 1) + yt(i, 2, s)*gchi(i, 2) + yt(i, 3, s)*gchi(i, 3)
         end do
         tau(s) = 0.5_dp*tau(s)
      end do
      vtau = 0.0_dp
      if (nspin == 1) then
         sigma = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         call xc_eval_point_mgga(nk, kid, coef, rho(1), sigma, tau(1), eps, vrho, vsigma, vtau(1))
         eg = wg*rho(1)*eps
         ng = wg*rho(1)
         zg(1) = wg*vrho
         zv(1, 1) = 2.0_dp*wg*vsigma*grho(1, 1)
         zv(2, 1) = 2.0_dp*wg*vsigma*grho(2, 1)
         zv(3, 1) = 2.0_dp*wg*vsigma*grho(3, 1)
      else
         s_aa = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         s_bb = grho(1, 2)**2 + grho(2, 2)**2 + grho(3, 2)**2
         s_ab = grho(1, 1)*grho(1, 2) + grho(2, 1)*grho(2, 2) + grho(3, 1)*grho(3, 2)
         call xc_eval_point_mgga_polar(nk, kid, coef, rho(1), rho(2), s_aa, s_ab, s_bb, tau(1), tau(2), &
                                       eps, va, vb, v_aa, v_ab, v_bb, vtau(1), vtau(2))
         eg = wg*(rho(1) + rho(2))*eps
         ng = wg*(rho(1) + rho(2))
         zg(1) = wg*va
         zg(2) = wg*vb
         zv(1, 1) = wg*(2.0_dp*v_aa*grho(1, 1) + v_ab*grho(1, 2))
         zv(2, 1) = wg*(2.0_dp*v_aa*grho(2, 1) + v_ab*grho(2, 2))
         zv(3, 1) = wg*(2.0_dp*v_aa*grho(3, 1) + v_ab*grho(3, 2))
         zv(1, 2) = wg*(2.0_dp*v_bb*grho(1, 2) + v_ab*grho(1, 1))
         zv(2, 2) = wg*(2.0_dp*v_bb*grho(2, 2) + v_ab*grho(2, 1))
         zv(3, 2) = wg*(2.0_dp*v_bb*grho(3, 2) + v_ab*grho(3, 1))
      end if
      ! W_x,u = w v_tau/4 d chi_u/dx: (d chi/dx)^T W plus its transpose is
      ! sum_g w v_tau/2 grad chi_u . grad chi_v.
      zw = 0.25_dp*wg*vtau
      do s = 1, nspin
         do u = 1, nloc
            i = off + (u - 1)*ld
            zt(i, s) = 0.5_dp*zg(s)*chi(i) + zv(1, s)*gchi(i, 1) + zv(2, s)*gchi(i, 2) + zv(3, s)*gchi(i, 3)
            yt(i, 1, s) = zw(s)*gchi(i, 1)
            yt(i, 2, s) = zw(s)*gchi(i, 2)
            yt(i, 3, s) = zw(s)*gchi(i, 3)
         end do
      end do
   end subroutine xc_point_mgga_body

   !
   ! Kernel 4: V_uv += V_b(u,v) + V_b(v,u) -- one thread per stored (b, u, v);
   ! a padding element has nothing to add and says so with au = 0.
   !
   subroutine xc_scatter(npairs, nb, ord, pb_off, nspin, nao, vxc, b_aooff, nbao, b_ao, nlp, nbatch, vloc)
      integer, intent(in) :: npairs, nb, nspin, nao, nbao, nbatch
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(inout) :: vxc(nao, nao, nspin)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(in) :: vloc(max(npairs, 1), nspin)
      integer :: t, au, av
      real(dp) :: acc_a, acc_b

      ! Two batches can share a global (u, v), so iterations write the same
      ! element: on the device the atomic makes that correct, but it is not
      ! what `do concurrent` promises, and a host compiler that takes the
      ! promise at its word -- ifx does -- reorders the read-modify-writes
      ! and the matrix comes out wrong while the energy, from the other
      ! kernel, is exact. So the concurrent form only under OpenACC, where
      ! the atomic exists; a plain loop otherwise. terco's Fock digestion
      ! is arranged the same way.
      ! Two complete loops rather than one with its header under #ifdef:
      ! the locality lint reads the source before the preprocessor does.
#ifdef _OPENACC
      do concurrent(t=1:npairs) local(au, av, acc_a, acc_b)
         call xc_scatter_body(t, nb, ord, pb_off, nspin, b_aooff, nbao, b_ao, nlp, nbatch, npairs, vloc, &
                              au, av, acc_a, acc_b)
         if (au > 0) then
            !$acc atomic update
            vxc(au, av, 1) = vxc(au, av, 1) + acc_a
            if (nspin == 2) then
               !$acc atomic update
               vxc(au, av, 2) = vxc(au, av, 2) + acc_b
            end if
         end if
      end do
#else
      do t = 1, npairs
         call xc_scatter_body(t, nb, ord, pb_off, nspin, b_aooff, nbao, b_ao, nlp, nbatch, npairs, vloc, &
                              au, av, acc_a, acc_b)
         if (au > 0) then
            vxc(au, av, 1) = vxc(au, av, 1) + acc_a
            if (nspin == 2) vxc(au, av, 2) = vxc(au, av, 2) + acc_b
         end if
      end do
#endif
   end subroutine xc_scatter

   pure subroutine xc_scatter_body(t, nb, ord, pb_off, nspin, b_aooff, nbao, b_ao, nlp, nbatch, npairs, vloc, &
                                   au, av, acc_a, acc_b)
      !$acc routine seq
      integer, intent(in) :: t, nb, nspin, nbao, nbatch, npairs
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(in) :: vloc(npairs, nspin)
      integer, intent(out) :: au, av
      real(dp), intent(out) :: acc_a, acc_b
      integer :: k, ib, tt, nloc, ld, u, v, tvu

      call pair_batch(t, nb, pb_off, k, tt)
      ib = ord(k)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = nlp(ib)
      u = mod(tt, ld) + 1
      v = tt/ld + 1
      acc_a = 0.0_dp
      acc_b = 0.0_dp
      au = 0
      av = 0
      if (u > nloc .or. v > nloc) return
      tvu = pb_off(k) + (u - 1)*ld + v   ! the (v, u) element
      acc_a = vloc(t, 1) + vloc(tvu, 1)
      if (nspin == 2) acc_b = vloc(t, 2) + vloc(tvu, 2)
      au = b_ao(b_aooff(ib) + u - 1)
      av = b_ao(b_aooff(ib) + v - 1)
   end subroutine xc_scatter_body

end module trc_xc
