!
! VV10 (Vydrov & Van Voorhis 2010) nonlocal correlation: an energy and a
! potential matrix, added onto whatever the semilocal functional already
! put there. wB97X-V uses it with b = 6.0, C = 0.01 (confirmed against
! pyscf's `dft.libxc.nlc_coeff('WB97X_V')`); this module takes b and C as
! arguments rather than hardcoding them, so a caller is free to run plain
! VV10 (b = 5.9, C = 0.0093) with the same code.
!
! THE FORMULAS, matched term for term against pyscf's `numint._vv10nlc`
! (the only place the kernel is written down outside the original paper):
!
!     Pi43 = 4 pi / 3
!     Kvv  = 1.5 b pi (9 pi)^(-1/6)
!     Beta = (3 / b^2)^0.75 / 32
!
!     W0(rho, sigma) = sqrt(C (sigma/rho^2)^2 + Pi43 rho)
!     K(rho)         = Kvv rho^(1/6)
!
! and, for every pair of points i (the one being evaluated) and j (summed
! over), with R2 = |r_i - r_j|^2,
!
!     g  = R2 W0(rho_i, sigma_i) + K(rho_i)
!     gp = R2 W0(rho_j, sigma_j) + K(rho_j)
!     gt = g + gp
!     F_i  = -1.5 sum_j w_j rho_j / (g gp gt)
!     U_i  =      sum_j w_j rho_j (1/g + 1/gt) / (g gp gt)
!     W_i  =      sum_j w_j rho_j (1/g + 1/gt) R2 / (g gp gt)
!
!     eps_i   = Beta + 0.5 F_i
!     vrho_i  = Beta + F_i + 1.5 (U_i dK/dR_i + W_i dW0/dR_i)
!     vsigma_i = 1.5 W_i dW0/dG_i
!
!     E_nlc = sum_i w_i rho_i eps_i
!
! A point with rho below `VV10_RHO_TOL` (1e-8, pyscf's own `thresh` in
! `_vv10nlc`) contributes nothing, neither as an i nor as a j -- dropped
! from the outer loop (its own eps/vrho/vsigma left zero) and skipped in
! every inner sum that would otherwise divide by its vanishing rho.
!
! The potential is then contracted exactly as trc_xc.F90's semilocal
! pipeline contracts vrho/vsigma into V_xc (same Z_u, same chi^T Z + its
! transpose): for singlet rho with sigma = |grad rho|^2,
!
!     Z_u = w (vrho/2 chi_u + 2 vsigma grad rho . grad chi_u)
!     V_uv = sum_g (chi_u(g) Z_v(g) + chi_v(g) Z_u(g))
!
! UKS CONVENTION (checked against pyscf's uks.py `get_veff`): VV10 is a
! functional of the TOTAL density rho_a + rho_b alone, and its potential
! -- a single nao x nao matrix -- is added UNCHANGED to both the alpha and
! beta Fock matrices (`vxc += vnlc` broadcasts the one matrix onto both
! spin channels there; `trc_vv10_uks` does the same).
!
! SCALING. Pass 2 below is the textbook O(N_grid^2) double sum over point
! pairs the paper itself is O(N^2) in -- no spatial cutoff, no screening,
! because VV10's kernel has no compact support the way a Coulomb operator's
! Schwarz bound has. This module is a capability, not a performance path:
! correct on a modest grid, and the caller chooses how coarse a grid VV10
! itself runs on (pyscf and this module both default to "coarser than the
! semilocal grid" for exactly that reason). Pass 1 and pass 3 are each
! O(N_grid x n_loc^2) -- the collocation-and-contraction cost the semilocal
! pipeline also pays -- so for any grid worth running VV10 on at all, pass
! 2 dominates total time by orders of magnitude.
!
! DEVICE RESIDENCY, matching trc_xc_rks/trc_xc_uks exactly: `bas` and
! `grid` are assumed ALREADY on the device (`bas%to_device()`,
! `grid%to_device()`, called by whoever built them -- this module reads
! their flat arrays without a single `copyin` of its own, the same
! assumption xc_drive makes of `b`/`xg`). `dmat`/`dmat_a`/`dmat_b` are
! `copyin` on entry (a no-op refcount bump if already resident); `vmat`
! and friends are `copy` -- present or not, the accumulation inside this
! call sees the value the caller already put there and the result comes
! back out whichever way it went in. Every other array (rho, grad rho,
! sigma, vrho, vsigma, the per-point energy) is local to this call: `enter
! data create`d here and deleted before returning.
!
! PASS 3's RACE. Two points scatter into the same (u, v) when they share a
! batch's local basis or lie in different batches whose local lists
! overlap -- exactly trc_xc's kernel-4 situation, and the same fix: the
! outer loop over points is `do concurrent` with an OpenACC atomic update
! only under `_OPENACC`; a plain sequential `do` otherwise, because `do
! concurrent`'s concurrency promise is not one every host compiler honours
! under a shared, non-atomic update (ifx reorders it and gets the matrix
! wrong while the energy, summed independently, stays exact). Passes 1 and
! 2 write only to the index they own (rho(g), grad rho(g); vrho(i),
! vsigma(i)) and need no such split.
!
! MULTI-GPU. With a communicator, pass 2's outer points are split into one
! CONTIGUOUS range [i0, i1] per rank -- npts/nranks points each, the first
! npts mod nranks ranks getting one extra -- not round-robin: every point
! at or above VV10_RHO_TOL costs the same O(npass) inner sum (one below it
! returns immediately -- vv10_kernel_body's own `if (ri < VV10_RHO_TOL)
! return` -- so its cost is O(1), not O(npass)), so a contiguous range is
! already load-balanced, and it keeps every GPU warp either entirely inside a
! rank's range or entirely outside it, which round-robin (mod(i-1,
! nranks)) does not -- round-robin interleaves an owned and a skipped
! point into the SAME warp throughout, so the warp's time is always the
! owned thread's time no matter how many ranks split the work, measured
! as no speedup at all between 1 and 2 ranks. A rank leaves its
! vrho(i)/vsigma(i)/eg(i) at the 0.0_dp vv10_kernel_body already sets
! before the range check, for every i outside its own [i0, i1], and one
! allreduce(MPI_SUM) afterward -- a single buffer packing all three
! arrays, as trc_xc.F90's sum_ranks does for vxc -- reconstructs the full
! arrays on every rank: each point was computed by exactly one rank and
! left 0.0_dp by every other, so SUM is exact, not an approximation. Passes
! 1 and 3 stay REPLICATED on every rank rather than split too -- they are
! the O(N) cost the SCALING section above says pass 2 dominates "by orders
! of magnitude", so redistributing them would buy nothing worth the extra
! communication, and running them identically from the identical
! (replicated) density is what lets pass 3's scatter into `vmat` come out
! identical on every rank with no reduction of its own. No comm given
! (i0 = 1, i1 = npts) makes every point everyone's, unchanged from the
! no-MPI path.
!
module trc_vv10
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_collocation, only: shell_collocate, NCART_MAX
   use trc_linalg, only: trc_linalg_t
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   implicit none
   private

   public :: trc_vv10_rks, trc_vv10_uks

   real(dp), parameter :: VV10_RHO_TOL = 1.0e-8_dp
   !! pyscf's `_vv10nlc` threshold, outer and inner grid alike: a density
   !! below this contributes nothing, neither evaluated nor summed over.
   real(dp), parameter :: PI = 3.14159265358979323846_dp
   real(dp), parameter :: PI43 = 4.0_dp*PI/3.0_dp
   integer, parameter :: VV10_NP_PAD = 128  !! tile rows padded to a multiple of this
   integer, parameter :: VV10_NL_PAD = 32   !! tile columns padded to a multiple of this
   ! Same padding constants trc_xc.F90 uses for its batch tiles (XC_NP_PAD/
   ! XC_NL_PAD) -- not shared by `use`, since they are two bare parameters
   ! and duplicating them costs nothing a cross-module dependency would be
   ! worth paying for.
   integer(int64), parameter :: VV10_TILE_CAP = 8388608_int64
   !! Bounds one tile chunk's chi/gchi/xt element count -- the running sum
   !! of npp(ib) * nlp(ib) over the batches folded into that chunk in
   !! passes 1 and 3 -- so the chunk's device scratch (chi/gchi/xt/dloc or
   !! zt/vloc) stays a bounded size; once a chunk would exceed this, it is
   !! closed and a new one starts (see vv10_density_tile/vv10_potential_tile).

contains

   !
   ! Closed shell: E_nlc and the VV10 potential, ADDED into `vmat` (not
   ! overwritten -- the caller has already put the semilocal V_xc there).
   ! `la`: GEMM/solver handles for passes 1 and 3's batched GEMMs, shared
   ! with the caller's own (trc_xc_rks's `la`) when given; made and
   ! released here otherwise, exactly as trc_xc's xc_drive does.
   !
   subroutine trc_vv10_rks(bas, grid, dmat, b, c, e_nlc, vmat, la, comm)
      type(trc_basis_t), intent(in) :: bas
      type(trc_xc_grid_t), intent(in) :: grid
      real(dp), intent(in) :: dmat(bas%nao, bas%nao)
      real(dp), intent(in) :: b, c   !! VV10's own b and C (6.0, 0.01 for wB97X-V)
      real(dp), intent(out) :: e_nlc
      real(dp), intent(inout) :: vmat(bas%nao, bas%nao)
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(comm_t), intent(in), optional :: comm
      !! Splits pass 2's outer points into a contiguous range per rank,
      !! one allreduce afterward -- see this module's own MULTI-GPU
      !! header for why passes 1 and 3 are left replicated on every rank
      !! rather than split too.

      call vv10_drive(bas%nao, bas%nshell, bas%maxnp, bas%sh_l, bas%sh_np, bas%sh_e, bas%sh_c, bas%sh_r, &
                      grid%npts, grid%nbatch, grid%r, grid%w, grid%batch_of, grid%b_off, grid%b_shoff, &
                      size(grid%b_sh), grid%b_sh, grid%b_shao, grid%b_aooff, size(grid%b_ao), grid%b_ao, &
                      dmat, b, c, e_nlc, vmat, la, comm)
   end subroutine trc_vv10_rks

   !
   ! Spin-polarised: VV10 sees only the total density D_a + D_b, and the
   ! SAME potential matrix is added to both `vmat_a` and `vmat_b` (pyscf's
   ! `uks.py` convention -- see this module's header).
   !
   subroutine trc_vv10_uks(bas, grid, dmat_a, dmat_b, b, c, e_nlc, vmat_a, vmat_b, la, comm)
      type(trc_basis_t), intent(in) :: bas
      type(trc_xc_grid_t), intent(in) :: grid
      real(dp), intent(in) :: dmat_a(bas%nao, bas%nao), dmat_b(bas%nao, bas%nao)
      real(dp), intent(in) :: b, c
      real(dp), intent(out) :: e_nlc
      real(dp), intent(inout) :: vmat_a(bas%nao, bas%nao), vmat_b(bas%nao, bas%nao)
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(comm_t), intent(in), optional :: comm
      real(dp), allocatable :: dtot(:, :), vtmp(:, :)
      integer :: nao, i, j

      nao = bas%nao
      allocate (dtot(nao, nao), vtmp(nao, nao))
      vtmp = 0.0_dp
      !$acc enter data create(dtot) copyin(vtmp)
      !$acc data copyin(dmat_a, dmat_b)
      do concurrent(i=1:nao, j=1:nao)
         dtot(i, j) = dmat_a(i, j) + dmat_b(i, j)
      end do
      !$acc end data

      call vv10_drive(bas%nao, bas%nshell, bas%maxnp, bas%sh_l, bas%sh_np, bas%sh_e, bas%sh_c, bas%sh_r, &
                      grid%npts, grid%nbatch, grid%r, grid%w, grid%batch_of, grid%b_off, grid%b_shoff, &
                      size(grid%b_sh), grid%b_sh, grid%b_shao, grid%b_aooff, size(grid%b_ao), grid%b_ao, &
                      dtot, b, c, e_nlc, vtmp, la, comm)

      !$acc data present(vtmp) copy(vmat_a, vmat_b)
      do concurrent(i=1:nao, j=1:nao)
         vmat_a(i, j) = vmat_a(i, j) + vtmp(i, j)
         vmat_b(i, j) = vmat_b(i, j) + vtmp(i, j)
      end do
      !$acc end data
      !$acc exit data delete(dtot, vtmp)
      deallocate (dtot, vtmp)
   end subroutine trc_vv10_uks

   !
   ! The three passes, on the flat arrays `bas`/`grid` carry -- never the
   ! derived types themselves, the same discipline xc_drive uses, since a
   ! derived-type component read from inside a `do concurrent`/`routine
   ! seq` is one of the two nvfortran traps this file has to avoid (the
   ! other is an assumed-shape dummy in a device routine; every array
   ! below is explicit-shape).
   !
   subroutine vv10_drive(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                         npts, nbatch, r, w, batch_of, b_off, b_shoff, nbsh, b_sh, b_shao, &
                         b_aooff, nbao, b_ao, dmat, b, c, e_nlc, vmat, la, comm)
      integer, intent(in) :: nao, nshell, maxnp, npts, nbatch, nbsh, nbao
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: r(3, npts), w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_shoff(nbatch + 1)
      integer, intent(in) :: b_sh(nbsh), b_shao(nbsh), b_aooff(nbatch + 1), b_ao(nbao)
      real(dp), intent(in) :: dmat(nao, nao)
      real(dp), intent(in) :: b, c
      real(dp), intent(out) :: e_nlc
      real(dp), intent(inout) :: vmat(nao, nao)
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(comm_t), intent(in), optional :: comm

      real(dp), allocatable :: rho(:), grho(:, :), sigma(:), vrho(:), vsigma(:), eg(:)
      real(dp), allocatable :: gx(:), gy(:), gz(:), gkp(:), gw0p(:), grpw(:)
      real(dp), allocatable :: buf(:)
      integer, allocatable :: gidx(:)
      real(dp) :: kvv, beta
      integer :: npass, g, rank, nranks, i0, i1
      integer(int64) :: cap
      logical :: split
      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp

      rank = 0
      nranks = 1
      split = .false.
      if (present(comm)) then
         rank = comm%rank()
         nranks = comm%size()
         split = nranks > 1
      end if

      e_nlc = 0.0_dp
      if (npts == 0) return

      ! A contiguous range [i0, i1] per rank -- see this module's own
      ! MULTI-GPU header for why.
      i0 = rank*(npts/nranks) + min(rank, mod(npts, nranks)) + 1
      i1 = (rank + 1)*(npts/nranks) + min(rank + 1, mod(npts, nranks))

      kvv = b*1.5_dp*PI*(9.0_dp*PI)**(-1.0_dp/6.0_dp)
      beta = (3.0_dp/(b*b))**0.75_dp/32.0_dp
      cap = VV10_TILE_CAP
      if (present(la)) then
         lp => la
      else
         call la_own%init(1)
         lp => la_own
      end if

      allocate (rho(npts), grho(3, npts), sigma(npts), vrho(npts), vsigma(npts), eg(npts))
      !$acc data copyin(dmat) create(rho, grho, sigma, vrho, vsigma, eg)

      call vv10_density_tile(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                             npts, nbatch, r, batch_of, b_off, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                             dmat, cap, lp, rho, grho, sigma)

      !
      ! GATHER before the O(N^2) pass: every quantity the inner sum needs
      ! from point j -- w0p(j), K(j), rho(j)*w(j), and r(j) itself -- is a
      ! function of j ALONE, yet the unguarded double loop recomputed a
      ! sqrt and a 1/6 power of it once per (i, j) PAIR, N times more often
      ! than the N values that exist. Below, one O(N) pass (`vv10_gather`)
      ! filters the points below VV10_RHO_TOL (the same ones the inner sum
      ! always `cycle`d past) and writes the survivors' precomputed
      ! quantities into contiguous, struct-of-arrays storage -- r split
      ! into gx/gy/gz rather than the strided r(3, npts) -- in ascending
      ! original order. Ascending order is what makes this bit-identical:
      ! the inner sum over j only ever adds terms in that order (skipping
      ! the ones below tolerance exactly as before), and inserting or
      ! removing a step that adds literally 0.0_dp -- which is what a
      ! below-tolerance j always contributed -- changes nothing at any
      ! position in a floating-point sum.
      !
      ! The filter itself runs on the host -- it is inherently sequential
      ! (building a compacted index list) and O(N), nothing a device launch
      ! would help -- so the density this call just built on the device is
      ! brought back first.
      !$acc update self(rho)
      allocate (gidx(max(npts, 1)))
      npass = 0
      do g = 1, npts
         if (rho(g) >= VV10_RHO_TOL) then
            npass = npass + 1
            gidx(npass) = g
         end if
      end do

      allocate (gx(max(npass, 1)), gy(max(npass, 1)), gz(max(npass, 1)))
      allocate (gkp(max(npass, 1)), gw0p(max(npass, 1)), grpw(max(npass, 1)))
      !$acc data copyin(gidx(1:max(npass, 1))) create(gx, gy, gz, gkp, gw0p, grpw)
      call vv10_gather(npts, r, w, rho, sigma, c, kvv, npass, gidx, gx, gy, gz, gkp, gw0p, grpw)

      ! MULTI-GPU: pass 2's outer points, a contiguous range per rank; see
      ! this module's own MULTI-GPU header for why passes 1 and 3 stay
      ! replicated instead.
      call vv10_kernel(npts, r, w, rho, sigma, b, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                       i0, i1, vrho, vsigma, eg)

      !$acc end data
      deallocate (gidx, gx, gy, gz, gkp, gw0p, grpw)

      !
      ! The reduction runs with no `!$acc data` construct open at all --
      ! the region that carried `vrho`/`vsigma`/`eg`/`grho` on the device
      ! for passes 1-2 is closed here first, the same thing trc_xc.F90's
      ! `sum_ranks` call site does (its `!$acc end data` closes xc_drive's
      ! chunk loop before `sum_ranks` is ever reached), and pass 3 reopens
      ! its own, narrower region below. A `block` placed inside the open
      ! region instead was rejected by nvfortran outright ("Unimplemented
      ! feature: BLOCK construct in the scope of a parallel directive"),
      ! which is why this is split out rather than folded into the data
      ! region above.
      !
      !$acc update self(vrho, vsigma, eg, grho)
      !$acc end data

      if (split) then
         allocate (buf(3*npts))
         buf(1:npts) = vrho
         buf(npts + 1:2*npts) = vsigma
         buf(2*npts + 1:3*npts) = eg
         call allreduce(comm, buf, op=MPI_SUM)
         ! Then rank 0's copy everywhere: an allreduce does not promise the
         ! same bits on every rank, same reason trc_xc.F90's sum_ranks
         ! bcasts after its own allreduce.
         call bcast(comm, buf, 3*npts, 0)
         vrho = buf(1:npts)
         vsigma = buf(npts + 1:2*npts)
         eg = buf(2*npts + 1:3*npts)
         deallocate (buf)
      end if

      e_nlc = sum(eg)

      !$acc data copy(vmat) copyin(grho, vrho, vsigma)
      call vv10_potential_tile(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                               npts, nbatch, r, w, batch_of, b_off, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                               vrho, vsigma, grho, cap, lp, vmat)
      !$acc end data

      deallocate (rho, grho, sigma, vrho, vsigma, eg)
      if (.not. present(la)) call la_own%release()
   end subroutine vv10_drive

   !
   ! Pass 1, TILED: rho(g) and grad rho(g) for every point, through the
   ! same cached-batch-tile + GEMM shape trc_xc.F90's semilocal pipeline
   ! uses (its own header: "Batched Dense D", Stocks and Barca JCTC 2025).
   ! Per chunk of batches: collocate chi/gchi ONCE per (local shell, point)
   ! -- not once per (shell PAIR, point), which is what the plain nested
   ! loop over local shells cost before this -- gather the density block
   ! D_b, GEMM X_b = chi_b D_b, then one thread per point reads its row of
   ! chi/gchi/X to form rho and grad rho. No atomics: every point writes
   ! only its own index, same as the loop form it replaces.
   !
   subroutine vv10_density_tile(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                npts, nbatch, r, batch_of, b_off, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                                dmat, cap, lp, rho, grho, sigma)
      integer, intent(in) :: nao, nshell, maxnp, npts, nbatch, nbsh, nbao
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao)
      real(dp), intent(in) :: dmat(nao, nao)
      integer(int64), intent(in) :: cap
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(out) :: rho(npts), grho(3, npts), sigma(npts)

      integer :: b0, b1, ib, nb, np, p0, nloc, npb, k, k0, k1, cnt, mp, ml
      integer, allocatable :: ord(:), npp(:), nlp(:), c_off(:), pr_off(:), pb_off(:), wk_off(:)
      integer(int64) :: need, nchi, sa, sb, nlsq
      integer :: nwork
      real(dp), allocatable :: chi(:), gchi(:, :), xt(:), dloc(:)

      allocate (npp(nbatch), nlp(nbatch), c_off(nbatch + 1), pr_off(nbatch + 1))
      do ib = 1, nbatch
         nloc = b_aooff(ib + 1) - b_aooff(ib)
         npb = b_off(ib + 1) - b_off(ib)
         npp(ib) = ((npb + VV10_NP_PAD - 1)/VV10_NP_PAD)*VV10_NP_PAD
         nlp(ib) = ((nloc + VV10_NL_PAD - 1)/VV10_NL_PAD)*VV10_NL_PAD
      end do
      !$acc data copyin(dmat) copy(rho, grho, sigma) copyin(npp, nlp)

      b0 = 1
      do while (b0 <= nbatch)
         need = 0
         b1 = b0 - 1
         do ib = b0, nbatch
            if (need + int(npp(ib), int64)*int(nlp(ib), int64) > cap .and. ib > b0) exit
            need = need + int(npp(ib), int64)*int(nlp(ib), int64)
            b1 = ib
         end do
         nchi = max(need, 1_int64)
         p0 = b_off(b0)
         np = b_off(b1 + 1) - p0
         nb = b1 - b0 + 1
         allocate (ord(nb), pb_off(nb + 1), wk_off(nb + 1))
         do k = 1, nb
            ord(k) = b0 + k - 1
         end do
         call vv10_sort_by_shape(nb, npp, nlp, ord)
         need = 0
         wk_off(1) = 0
         pb_off(1) = 0
         nwork = 0
         do k = 1, nb
            ib = ord(k)
            ! need and nlp(ib)*nlp(ib) are formed in int64 before the narrowing
            ! cast into c_off/pb_off (both plain `integer`, the kind every
            ! device kernel below indexes chi/gchi/xt/dloc with): guard each
            ! cast so a tile chunk that would overflow int32 dies with a clear
            ! message instead of wrapping into a corrupted offset.
            if (need > int(huge(0), int64)) then
               error stop "trc_vv10_density_tile: tile chunk offset overflowed int32; lower VV10_TILE_CAP"
            end if
            c_off(ib) = int(need)
            pr_off(ib) = pb_off(k)
            need = need + int(npp(ib), int64)*int(nlp(ib), int64)
            nlsq = int(nlp(ib), int64)*int(nlp(ib), int64)
            if (int(pb_off(k), int64) + nlsq > int(huge(0), int64)) then
               error stop "trc_vv10_density_tile: tile pair-block offset overflowed int32; lower VV10_TILE_CAP"
            end if
            pb_off(k + 1) = pb_off(k) + nlp(ib)*nlp(ib)
            nwork = nwork + (b_shoff(ib + 1) - b_shoff(ib))*(b_off(ib + 1) - b_off(ib))
            wk_off(k + 1) = nwork
         end do

         allocate (chi(nchi), gchi(nchi, 3), xt(nchi), dloc(max(pb_off(nb + 1), 1)))
         !$acc enter data create(chi, gchi, xt, dloc) copyin(c_off, pr_off, ord, pb_off, wk_off)

         call vv10_tile_zero(nchi, chi)
         call vv10_tile_zero(int(pb_off(nb + 1), int64), dloc)
         call vv10_tile_collocate(nwork, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                  npts, r, b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, nchi, chi, gchi)
         call vv10_tile_gather(pb_off(nb + 1), nb, ord, pb_off, nao, dmat, b_aooff, nbao, b_ao, &
                               nlp, nbatch, dloc)

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
               call lp%gemm_strided("N", "N", mp, ml, ml, 1.0_dp, chi(c_off(ib) + 1), mp, sa, &
                                    dloc(pr_off(ib) + 1), ml, sb, 0.0_dp, xt(c_off(ib) + 1), mp, sa, cnt)
            end if
            k0 = k1 + 1
         end do

         call vv10_tile_rho(np, p0, npts, batch_of, b_off, c_off, npp, b_aooff, nbatch, nchi, chi, gchi, xt, &
                            rho, grho, sigma)

         !$acc exit data delete(chi, gchi, xt, dloc, c_off, pr_off, ord, pb_off, wk_off)
         deallocate (chi, gchi, xt, dloc, ord, pb_off, wk_off)
         b0 = b1 + 1
      end do
      !$acc end data
      deallocate (npp, nlp, c_off, pr_off)
   end subroutine vv10_density_tile

   !
   ! Collocation: one thread per (local shell, point) work item, exactly
   ! trc_xc.F90's xc_collocate (see its own header) with no nspin/padding
   ! difference -- VV10 always runs on one density, so this is that
   ! kernel's shape kept independent rather than imported, the same
   ! module-boundary choice this file's header explains for shell_collocate.
   !
   subroutine vv10_tile_collocate(nwork, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
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
         call vv10_tile_collocate_body(t, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                       b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, nchi, chi, gchi)
      end do
   end subroutine vv10_tile_collocate

   pure subroutine vv10_tile_collocate_body(t, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
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

      call vv10_pair_batch(t, nb, wk_off, k, tt)
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
   end subroutine vv10_tile_collocate_body

   !
   ! Density gather: one thread per stored (b, u, v), padding included and
   ! written as zero -- trc_xc's xc_gather with nspin dropped (VV10 runs
   ! on one density, built by `trc_vv10_uks` from D_a + D_b beforehand).
   !
   subroutine vv10_tile_gather(npairs, nb, ord, pb_off, nao, dmat, b_aooff, nbao, b_ao, nlp, nbatch, dloc)
      integer, intent(in) :: npairs, nb, nao, nbao, nbatch
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(in) :: dmat(nao, nao)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(out) :: dloc(max(npairs, 1))
      integer :: t

      do concurrent(t=1:npairs)
         call vv10_tile_gather_body(t, nb, ord, pb_off, nao, dmat, b_aooff, nbao, b_ao, nlp, nbatch, npairs, dloc)
      end do
   end subroutine vv10_tile_gather

   pure subroutine vv10_tile_gather_body(t, nb, ord, pb_off, nao, dmat, b_aooff, nbao, b_ao, nlp, nbatch, npairs, dloc)
      !$acc routine seq
      integer, intent(in) :: t, nb, nao, nbao, nbatch, npairs
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(in) :: dmat(nao, nao)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(inout) :: dloc(npairs)
      integer :: k, ib, tt, nloc, ld, u, v

      call vv10_pair_batch(t, nb, pb_off, k, tt)
      ib = ord(k)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = nlp(ib)
      u = mod(tt, ld) + 1
      v = tt/ld + 1
      if (u > nloc .or. v > nloc) then
         dloc(t) = 0.0_dp
      else
         dloc(t) = dmat(b_ao(b_aooff(ib) + u - 1), b_ao(b_aooff(ib) + v - 1))
      end if
   end subroutine vv10_tile_gather_body

   !
   ! rho(g)/grad rho(g)/sigma(g) from the tile: one thread per point,
   ! reading its row of chi/gchi/X -- same arithmetic as the old
   ! vv10_density_body's inner accumulation, just against a tile instead
   ! of a freshly-collocated pair.
   !
   subroutine vv10_tile_rho(np, p0, npts, batch_of, b_off, c_off, npp, b_aooff, nbatch, nchi, chi, gchi, xt, &
                            rho, grho, sigma)
      integer, intent(in) :: np, p0, npts, nbatch
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), c_off(nbatch + 1), npp(nbatch), b_aooff(nbatch + 1)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3), xt(nchi)
      real(dp), intent(inout) :: rho(npts), grho(3, npts), sigma(npts)
      integer :: g

      do concurrent(g=1:np)
         call vv10_tile_rho_body(p0 + g - 1, npts, batch_of, b_off, c_off, npp, b_aooff, nbatch, nchi, chi, gchi, xt, &
                                 rho(p0 + g - 1), grho(1, p0 + g - 1), grho(2, p0 + g - 1), grho(3, p0 + g - 1), &
                                 sigma(p0 + g - 1))
      end do
   end subroutine vv10_tile_rho

   pure subroutine vv10_tile_rho_body(gg, npts, batch_of, b_off, c_off, npp, b_aooff, nbatch, nchi, chi, gchi, xt, &
                                      rho_g, gx, gy, gz, sig_g)
      !$acc routine seq
      integer, intent(in) :: gg, npts, nbatch
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), c_off(nbatch + 1), npp(nbatch), b_aooff(nbatch + 1)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3), xt(nchi)
      real(dp), intent(out) :: rho_g, gx, gy, gz, sig_g
      integer :: ib, nloc, ld, off, u, i
      real(dp) :: x

      ib = batch_of(gg)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = npp(ib)
      off = c_off(ib) + (gg - b_off(ib)) + 1
      rho_g = 0.0_dp
      gx = 0.0_dp
      gy = 0.0_dp
      gz = 0.0_dp
      do u = 1, nloc
         i = off + (u - 1)*ld
         x = xt(i)
         rho_g = rho_g + x*chi(i)
         gx = gx + 2.0_dp*x*gchi(i, 1)
         gy = gy + 2.0_dp*x*gchi(i, 2)
         gz = gz + 2.0_dp*x*gchi(i, 3)
      end do
      sig_g = gx*gx + gy*gy + gz*gz
   end subroutine vv10_tile_rho_body

   !
   ! Gather: one thread per SURVIVING j (`npass` of them, `gidx` holds their
   ! original point indices in ascending order), precomputing everything
   ! pass 2's inner sum reads from j alone -- its position split into
   ! gx/gy/gz (contiguous, unit-stride reads for the pair loop below,
   ! unlike the strided r(3, npts) this is drawn from), K(j) = kvv rho(j)^
   ! (1/6), W0(j) and rho(j) w(j). See vv10_drive's call site for why this
   ! makes pass 2 bit-identical rather than merely faster.
   !
   subroutine vv10_gather(npts, r, w, rho, sigma, c, kvv, npass, gidx, gx, gy, gz, gkp, gw0p, grpw)
      integer, intent(in) :: npts, npass
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: c, kvv
      integer, intent(in) :: gidx(npass)
      real(dp), intent(out) :: gx(npass), gy(npass), gz(npass), gkp(npass), gw0p(npass), grpw(npass)
      integer :: k

      do concurrent(k=1:npass)
         call vv10_gather_body(k, npts, r, w, rho, sigma, c, kvv, npass, gidx, &
                               gx(k), gy(k), gz(k), gkp(k), gw0p(k), grpw(k))
      end do
   end subroutine vv10_gather

   pure subroutine vv10_gather_body(k, npts, r, w, rho, sigma, c, kvv, npass, gidx, &
                                    gx_k, gy_k, gz_k, gkp_k, gw0p_k, grpw_k)
      !$acc routine seq
      integer, intent(in) :: k, npts, npass
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: c, kvv
      integer, intent(in) :: gidx(npass)
      real(dp), intent(out) :: gx_k, gy_k, gz_k, gkp_k, gw0p_k, grpw_k
      integer :: j
      real(dp) :: rp, gp, w0ptmp

      j = gidx(k)
      gx_k = r(1, j)
      gy_k = r(2, j)
      gz_k = r(3, j)
      rp = rho(j)
      gp = sigma(j)
      w0ptmp = c*(gp/(rp*rp))**2
      gw0p_k = sqrt(w0ptmp + PI43*rp)
      gkp_k = kvv*rp**(1.0_dp/6.0_dp)
      grpw_k = rp*w(j)
   end subroutine vv10_gather_body

   !
   ! Pass 2: the O(N^2) double sum -- one launch, the outer point i calling
   ! vv10_kernel_body once each. The O(N) reduction over j that builds
   ! vrho(i)/vsigma(i)/eg(i) is entirely inside that one call, owned by
   ! iteration i alone -- nothing here can legally alias across i, so this
   ! is as race-free as pass 1, and `do concurrent` promises the compiler
   ! nothing a plain `do` does not already guarantee.
   !
   ! GPU form: `do concurrent`, same as pass 1 and pass 3's GPU form.
   !
   ! CPU form: a plain sequential `do`, same split as pass 3 but for a
   ! different reason, found by bisection rather than read off the standard.
   ! ifx (2026.0, Release, the CI host build) makes vrho(i)/vsigma(i)/eg(i)
   ! come out exactly 0 -- not noisy, not NaN, silently dropped -- but ONLY
   ! when the value assigned to one of them is built from vv10_kernel_body's
   ! own `do j` reduction (fsum/usum/wsum) AND the call sits inside a `do
   ! concurrent`. The same call through a plain `do` is correct to pyscf's
   ! precision; an output that does not depend on that reduction (a
   ! constant, or a bare dummy argument) survives `do concurrent` fine; and
   ! the loop itself still runs either way (its local accumulators compute
   ! correctly, they just never reach the caller's array through `do
   ! concurrent`). That is a mis-optimisation of ifx's inlined-reduction
   ! code path for this one pattern, not a standards question -- gfortran
   ! and nvfortran both compute pass 2 correctly through `do concurrent`
   ! already, which is why only the CPU form changes here.
   !
#ifdef _OPENACC
   subroutine vv10_kernel(npts, r, w, rho, sigma, b, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                          i0, i1, vrho, vsigma, eg)
      integer, intent(in) :: npts, npass, i0, i1
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: b, c, kvv, beta
      real(dp), intent(in) :: gx(npass), gy(npass), gz(npass), gkp(npass), gw0p(npass), grpw(npass)
      real(dp), intent(out) :: vrho(npts), vsigma(npts), eg(npts)
      integer :: i

      do concurrent(i=1:npts)
         call vv10_kernel_body(i, npts, r, w, rho, sigma, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                               i0, i1, vrho(i), vsigma(i), eg(i))
      end do
   end subroutine vv10_kernel
#else
   subroutine vv10_kernel(npts, r, w, rho, sigma, b, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                          i0, i1, vrho, vsigma, eg)
      integer, intent(in) :: npts, npass, i0, i1
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: b, c, kvv, beta
      real(dp), intent(in) :: gx(npass), gy(npass), gz(npass), gkp(npass), gw0p(npass), grpw(npass)
      real(dp), intent(out) :: vrho(npts), vsigma(npts), eg(npts)
      integer :: i

      do i = 1, npts
         call vv10_kernel_body(i, npts, r, w, rho, sigma, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                               i0, i1, vrho(i), vsigma(i), eg(i))
      end do
   end subroutine vv10_kernel
#endif

   !
   ! `j`'s own work (w0p(j), K(j), rho(j) w(j), its gathered position) now
   ! comes in through gx/gy/gz/gkp/gw0p/grpw, precomputed once by
   ! vv10_gather instead of once per (i, j) pair here -- the sqrt and the
   ! 1/6 power that used to run inside this `do j` are gone from it, and
   ! the inner loop itself is `do j = 1, npass` (survivors only), not
   ! `do j = 1, npts` with a `cycle`. Everything else, including the sum
   ! order (gidx is ascending), is unchanged.
   !
   pure subroutine vv10_kernel_body(i, npts, r, w, rho, sigma, c, kvv, beta, npass, gx, gy, gz, gkp, gw0p, grpw, &
                                    i0, i1, vrho_i, vsigma_i, eg_i)
      !$acc routine seq
      integer, intent(in) :: i, npts, npass, i0, i1
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: c, kvv, beta
      real(dp), intent(in) :: gx(npass), gy(npass), gz(npass), gkp(npass), gw0p(npass), grpw(npass)
      real(dp), intent(out) :: vrho_i, vsigma_i, eg_i

      integer :: j
      real(dp) :: ri, gi, w0tmp, w0, dw0dr, dw0dg, kk, dkdr
      real(dp) :: dx, dy, dz, r2
      real(dp) :: gg, ggp, gt, tt, t2, fsum, usum, wsum, ff

      vrho_i = 0.0_dp
      vsigma_i = 0.0_dp
      eg_i = 0.0_dp
      ! MULTI-GPU: point i belongs to this rank only if it falls in
      ! [i0, i1] (vv10_drive's own contiguous split) -- every other rank
      ! leaves it at the 0.0_dp just set above; see this module's own
      ! MULTI-GPU header for why. i0 = 1, i1 = npts with no comm given.
      if (i < i0 .or. i > i1) return
      ri = rho(i)
      if (ri < VV10_RHO_TOL) return
      gi = sigma(i)
      w0tmp = c*(gi/(ri*ri))**2
      w0 = sqrt(w0tmp + PI43*ri)
      dw0dr = (0.5_dp*PI43*ri - 2.0_dp*w0tmp)/w0
      ! Algebraically w0tmp*ri/(gi*w0) = c*(gi/ri**2)**2*ri/(gi*w0)
      ! = c*gi/(ri**3*w0): the gi in the numerator and denominator cancel
      ! down to one power rather than none, which is what makes this form
      ! finite at sigma = 0 (gi = 0, ri >= VV10_RHO_TOL): the old form was
      ! literally 0/0 there, since w0tmp ~ gi**2 vanishes together with the
      ! gi it was divided by.
      dw0dg = c*gi/(ri**3*w0)
      kk = kvv*ri**(1.0_dp/6.0_dp)
      dkdr = kk/6.0_dp

      fsum = 0.0_dp
      usum = 0.0_dp
      wsum = 0.0_dp
      do j = 1, npass
         dx = r(1, i) - gx(j)
         dy = r(2, i) - gy(j)
         dz = r(3, i) - gz(j)
         r2 = dx*dx + dy*dy + dz*dz
         ggp = r2*gw0p(j) + gkp(j)
         gg = r2*w0 + kk
         gt = gg + ggp
         tt = grpw(j)/(gg*ggp*gt)
         fsum = fsum + tt
         t2 = tt*(1.0_dp/gg + 1.0_dp/gt)
         usum = usum + t2
         wsum = wsum + t2*r2
      end do
      ff = -1.5_dp*fsum
      eg_i = w(i)*ri*(beta + 0.5_dp*ff)
      vrho_i = beta + ff + 1.5_dp*(usum*dkdr + wsum*dw0dr)
      vsigma_i = 1.5_dp*wsum*dw0dg
   end subroutine vv10_kernel_body

   !
   ! Pass 3, TILED: the contraction, through the same tile shape pass 1
   ! above uses -- collocate chi/gchi per chunk, form the Z tile from
   ! vrho/vsigma/grad rho (one thread per point, looping its batch's local
   ! functions), GEMM V_b = chi_b^T Z_b, then scatter V_b + V_b^T into
   ! `vmat` exactly as trc_xc's kernel 4 does (its own header: the atomic
   ! form only under `_OPENACC`, a plain loop otherwise, because `do
   ! concurrent` promises a host compiler nothing about the shared,
   ! non-atomic update two batches' overlapping local lists create here).
   !
   subroutine vv10_potential_tile(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                  npts, nbatch, r, w, batch_of, b_off, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                                  vrho, vsigma, grho, cap, lp, vmat)
      integer, intent(in) :: nao, nshell, maxnp, npts, nbatch, nbsh, nbao
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts), w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao)
      real(dp), intent(in) :: vrho(npts), vsigma(npts), grho(3, npts)
      integer(int64), intent(in) :: cap
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(inout) :: vmat(nao, nao)

      integer :: b0, b1, ib, nb, np, p0, nloc, npb, k, k0, k1, cnt, mp, ml
      integer, allocatable :: ord(:), npp(:), nlp(:), c_off(:), pr_off(:), pb_off(:), wk_off(:)
      integer(int64) :: need, nchi, sa, sc, nlsq
      integer :: nwork
      real(dp), allocatable :: chi(:), gchi(:, :), zt(:), vloc(:)

      allocate (npp(nbatch), nlp(nbatch), c_off(nbatch + 1), pr_off(nbatch + 1))
      do ib = 1, nbatch
         nloc = b_aooff(ib + 1) - b_aooff(ib)
         npb = b_off(ib + 1) - b_off(ib)
         npp(ib) = ((npb + VV10_NP_PAD - 1)/VV10_NP_PAD)*VV10_NP_PAD
         nlp(ib) = ((nloc + VV10_NL_PAD - 1)/VV10_NL_PAD)*VV10_NL_PAD
      end do
      !$acc data copyin(npp, nlp)

      b0 = 1
      do while (b0 <= nbatch)
         need = 0
         b1 = b0 - 1
         do ib = b0, nbatch
            if (need + int(npp(ib), int64)*int(nlp(ib), int64) > cap .and. ib > b0) exit
            need = need + int(npp(ib), int64)*int(nlp(ib), int64)
            b1 = ib
         end do
         nchi = max(need, 1_int64)
         p0 = b_off(b0)
         np = b_off(b1 + 1) - p0
         nb = b1 - b0 + 1
         allocate (ord(nb), pb_off(nb + 1), wk_off(nb + 1))
         do k = 1, nb
            ord(k) = b0 + k - 1
         end do
         call vv10_sort_by_shape(nb, npp, nlp, ord)
         need = 0
         wk_off(1) = 0
         pb_off(1) = 0
         nwork = 0
         do k = 1, nb
            ib = ord(k)
            ! Same int32 overflow guard as vv10_density_tile: need and
            ! nlp(ib)*nlp(ib) are formed in int64 before narrowing into
            ! c_off/pb_off, the kind every device kernel below indexes
            ! chi/gchi/zt/vloc with.
            if (need > int(huge(0), int64)) then
               error stop "trc_vv10_potential_tile: tile chunk offset overflowed int32; lower VV10_TILE_CAP"
            end if
            c_off(ib) = int(need)
            pr_off(ib) = pb_off(k)
            need = need + int(npp(ib), int64)*int(nlp(ib), int64)
            nlsq = int(nlp(ib), int64)*int(nlp(ib), int64)
            if (int(pb_off(k), int64) + nlsq > int(huge(0), int64)) then
               error stop "trc_vv10_potential_tile: tile pair-block offset overflowed int32; lower VV10_TILE_CAP"
            end if
            pb_off(k + 1) = pb_off(k) + nlp(ib)*nlp(ib)
            nwork = nwork + (b_shoff(ib + 1) - b_shoff(ib))*(b_off(ib + 1) - b_off(ib))
            wk_off(k + 1) = nwork
         end do

         allocate (chi(nchi), gchi(nchi, 3), zt(nchi), vloc(max(pb_off(nb + 1), 1)))
         !$acc enter data create(chi, gchi, zt, vloc) copyin(c_off, pr_off, ord, pb_off, wk_off)

         call vv10_tile_zero(nchi, chi)
         call vv10_tile_zero(nchi, zt)
         call vv10_tile_collocate(nwork, nb, ord, wk_off, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                  npts, r, b_off, b_shoff, nbsh, b_sh, b_shao, c_off, npp, nbatch, nchi, chi, gchi)
         call vv10_tile_z(np, p0, npts, batch_of, b_off, c_off, npp, b_aooff, nbatch, w, vrho, vsigma, grho, &
                          nchi, chi, gchi, zt)

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
                                    zt(c_off(ib) + 1), mp, sa, 0.0_dp, vloc(pr_off(ib) + 1), ml, sc, cnt)
            end if
            k0 = k1 + 1
         end do

         call vv10_tile_scatter(pb_off(nb + 1), nb, ord, pb_off, nao, vmat, b_aooff, nbao, b_ao, nlp, nbatch, vloc)

         !$acc exit data delete(chi, gchi, zt, vloc, c_off, pr_off, ord, pb_off, wk_off)
         deallocate (chi, gchi, zt, vloc, ord, pb_off, wk_off)
         b0 = b1 + 1
      end do
      !$acc end data
      deallocate (npp, nlp, c_off, pr_off)
   end subroutine vv10_potential_tile

   !
   ! Z tile: one thread per point, writing Z_gu for every local function u
   ! of its batch -- the same Z_u trc_xc's kernel 3 forms for nspin = 1,
   ! just against chi/gchi read from the tile instead of a fresh
   ! `shell_collocate`.
   !
   subroutine vv10_tile_z(np, p0, npts, batch_of, b_off, c_off, npp, b_aooff, nbatch, w, vrho, vsigma, grho, &
                          nchi, chi, gchi, zt)
      integer, intent(in) :: np, p0, npts, nbatch
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), c_off(nbatch + 1), npp(nbatch), b_aooff(nbatch + 1)
      real(dp), intent(in) :: w(npts), vrho(npts), vsigma(npts), grho(3, npts)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3)
      real(dp), intent(inout) :: zt(nchi)
      integer :: g

      do concurrent(g=1:np)
         call vv10_tile_z_body(p0 + g - 1, npts, batch_of, b_off, c_off, npp, b_aooff, nbatch, w, vrho, vsigma, grho, &
                               nchi, chi, gchi, zt)
      end do
   end subroutine vv10_tile_z

   pure subroutine vv10_tile_z_body(gg, npts, batch_of, b_off, c_off, npp, b_aooff, nbatch, w, vrho, vsigma, grho, &
                                    nchi, chi, gchi, zt)
      !$acc routine seq
      integer, intent(in) :: gg, npts, nbatch
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), c_off(nbatch + 1), npp(nbatch), b_aooff(nbatch + 1)
      real(dp), intent(in) :: w(npts), vrho(npts), vsigma(npts), grho(3, npts)
      integer(int64), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(nchi, 3)
      real(dp), intent(inout) :: zt(nchi)
      integer :: ib, nloc, ld, off, u, i
      real(dp) :: zg, zvx, zvy, zvz

      ib = batch_of(gg)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = npp(ib)
      off = c_off(ib) + (gg - b_off(ib)) + 1
      zg = w(gg)*vrho(gg)
      zvx = 2.0_dp*w(gg)*vsigma(gg)*grho(1, gg)
      zvy = 2.0_dp*w(gg)*vsigma(gg)*grho(2, gg)
      zvz = 2.0_dp*w(gg)*vsigma(gg)*grho(3, gg)
      do u = 1, nloc
         i = off + (u - 1)*ld
         zt(i) = 0.5_dp*zg*chi(i) + zvx*gchi(i, 1) + zvy*gchi(i, 2) + zvz*gchi(i, 3)
      end do
   end subroutine vv10_tile_z_body

   !
   ! Scatter: V_uv += V_b(u,v) + V_b(v,u) -- one thread per stored (b, u, v);
   ! a padding element has nothing to add and says so with au = 0. Trc_xc's
   ! kernel 4 exactly, nspin dropped: see this module's header (PASS 3's
   ! RACE) for why the `do concurrent` form only runs under `_OPENACC`.
   !
   subroutine vv10_tile_scatter(npairs, nb, ord, pb_off, nao, vmat, b_aooff, nbao, b_ao, nlp, nbatch, vloc)
      integer, intent(in) :: npairs, nb, nao, nbao, nbatch
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      real(dp), intent(inout) :: vmat(nao, nao)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(in) :: vloc(max(npairs, 1))
      integer :: t, au, av
      real(dp) :: acc

#ifdef _OPENACC
      do concurrent(t=1:npairs) local(au, av, acc)
         call vv10_tile_scatter_body(t, nb, ord, pb_off, b_aooff, nbao, b_ao, nlp, nbatch, npairs, vloc, au, av, acc)
         if (au > 0) then
            !$acc atomic update
            vmat(au, av) = vmat(au, av) + acc
         end if
      end do
#else
      do t = 1, npairs
         call vv10_tile_scatter_body(t, nb, ord, pb_off, b_aooff, nbao, b_ao, nlp, nbatch, npairs, vloc, au, av, acc)
         if (au > 0) vmat(au, av) = vmat(au, av) + acc
      end do
#endif
   end subroutine vv10_tile_scatter

   pure subroutine vv10_tile_scatter_body(t, nb, ord, pb_off, b_aooff, nbao, b_ao, nlp, nbatch, npairs, vloc, &
                                          au, av, acc)
      !$acc routine seq
      integer, intent(in) :: t, nb, nbao, nbatch, npairs
      integer, intent(in) :: ord(nb), pb_off(nb + 1)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), nlp(nbatch)
      real(dp), intent(in) :: vloc(npairs)
      integer, intent(out) :: au, av
      real(dp), intent(out) :: acc
      integer :: k, ib, tt, nloc, ld, u, v, tvu

      call vv10_pair_batch(t, nb, pb_off, k, tt)
      ib = ord(k)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      ld = nlp(ib)
      u = mod(tt, ld) + 1
      v = tt/ld + 1
      acc = 0.0_dp
      au = 0
      av = 0
      if (u > nloc .or. v > nloc) return
      tvu = pb_off(k) + (u - 1)*ld + v   ! the (v, u) element
      acc = vloc(t) + vloc(tvu)
      au = b_ao(b_aooff(ib) + u - 1)
      av = b_ao(b_aooff(ib) + v - 1)
   end subroutine vv10_tile_scatter_body

   ! Zero a flat tile array -- the padding rows/columns the GEMMs read
   ! must be zero, not whatever the allocation held.
   subroutine vv10_tile_zero(n, a)
      integer(int64), intent(in) :: n
      real(dp), intent(out) :: a(n)
      integer(int64) :: i
      do concurrent(i=1:n)
         a(i) = 0.0_dp
      end do
   end subroutine vv10_tile_zero

   !
   ! Which stored batch (1..nb, storage order) owns pair/work-item t: the
   ! last k with off(k) < t. trc_xc's `pair_batch`, renamed so this
   ! module's own copy cannot be confused with trc_xc's across a future
   ! refactor that exports one of them.
   !
   pure subroutine vv10_pair_batch(t, nb, off, k, tt)
      !$acc routine seq
      integer, intent(in) :: t, nb
      integer, intent(in) :: off(nb + 1)
      integer, intent(out) :: k, tt
      integer :: lo, hi, mid
      lo = 1
      hi = nb
      do while (lo < hi)
         mid = (lo + hi + 1)/2
         if (off(mid) < t) then
            lo = mid
         else
            hi = mid - 1
         end if
      end do
      k = lo
      tt = t - off(k) - 1
   end subroutine vv10_pair_batch

   ! Insertion sort of `ord` by (npp, nlp) ascending -- trc_xc's
   ! `sort_by_shape`; a chunk's worth of batches, a few hundred at most.
   subroutine vv10_sort_by_shape(n, npp, nlp, ord)
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
   end subroutine vv10_sort_by_shape

end module trc_vv10
