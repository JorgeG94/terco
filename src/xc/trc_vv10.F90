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
module trc_vv10
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_collocation, only: shell_collocate, NCART_MAX
   implicit none
   private

   public :: trc_vv10_rks, trc_vv10_uks

   real(dp), parameter :: VV10_RHO_TOL = 1.0e-8_dp
   !! pyscf's `_vv10nlc` threshold, outer and inner grid alike: a density
   !! below this contributes nothing, neither evaluated nor summed over.
   real(dp), parameter :: PI = 3.14159265358979323846_dp
   real(dp), parameter :: PI43 = 4.0_dp*PI/3.0_dp

contains

   !
   ! Closed shell: E_nlc and the VV10 potential, ADDED into `vmat` (not
   ! overwritten -- the caller has already put the semilocal V_xc there).
   !
   subroutine trc_vv10_rks(bas, grid, dmat, b, c, e_nlc, vmat)
      type(trc_basis_t), intent(in) :: bas
      type(trc_xc_grid_t), intent(in) :: grid
      real(dp), intent(in) :: dmat(bas%nao, bas%nao)
      real(dp), intent(in) :: b, c   !! VV10's own b and C (6.0, 0.01 for wB97X-V)
      real(dp), intent(out) :: e_nlc
      real(dp), intent(inout) :: vmat(bas%nao, bas%nao)

      call vv10_drive(bas%nao, bas%nshell, bas%maxnp, bas%sh_l, bas%sh_np, bas%sh_e, bas%sh_c, bas%sh_r, &
                      grid%npts, grid%nbatch, grid%r, grid%w, grid%batch_of, grid%b_off, grid%b_shoff, &
                      size(grid%b_sh), grid%b_sh, grid%b_shao, grid%b_aooff, size(grid%b_ao), grid%b_ao, &
                      dmat, b, c, e_nlc, vmat)
   end subroutine trc_vv10_rks

   !
   ! Spin-polarised: VV10 sees only the total density D_a + D_b, and the
   ! SAME potential matrix is added to both `vmat_a` and `vmat_b` (pyscf's
   ! `uks.py` convention -- see this module's header).
   !
   subroutine trc_vv10_uks(bas, grid, dmat_a, dmat_b, b, c, e_nlc, vmat_a, vmat_b)
      type(trc_basis_t), intent(in) :: bas
      type(trc_xc_grid_t), intent(in) :: grid
      real(dp), intent(in) :: dmat_a(bas%nao, bas%nao), dmat_b(bas%nao, bas%nao)
      real(dp), intent(in) :: b, c
      real(dp), intent(out) :: e_nlc
      real(dp), intent(inout) :: vmat_a(bas%nao, bas%nao), vmat_b(bas%nao, bas%nao)
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
                      dtot, b, c, e_nlc, vtmp)

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
                         b_aooff, nbao, b_ao, dmat, b, c, e_nlc, vmat)
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

      real(dp), allocatable :: rho(:), grho(:, :), sigma(:), vrho(:), vsigma(:), eg(:)
      real(dp) :: kvv, beta

      e_nlc = 0.0_dp
      if (npts == 0) return

      kvv = b*1.5_dp*PI*(9.0_dp*PI)**(-1.0_dp/6.0_dp)
      beta = (3.0_dp/(b*b))**0.75_dp/32.0_dp

      allocate (rho(npts), grho(3, npts), sigma(npts), vrho(npts), vsigma(npts), eg(npts))
      !$acc data copyin(dmat) copy(vmat) create(rho, grho, sigma, vrho, vsigma, eg)

      call vv10_density(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                        npts, nbatch, r, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                        dmat, rho, grho, sigma)

      call vv10_kernel(npts, r, w, rho, sigma, b, c, kvv, beta, vrho, vsigma, eg)

      !$acc update self(eg)
      e_nlc = sum(eg)

#ifdef _OPENACC
      call vv10_potential_gpu(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                              npts, nbatch, r, w, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                              vrho, vsigma, grho, vmat)
#else
      call vv10_potential_cpu(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                              npts, nbatch, r, w, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                              vrho, vsigma, grho, vmat)
#endif

      !$acc end data
      deallocate (rho, grho, sigma, vrho, vsigma, eg)
   end subroutine vv10_drive

   !
   ! Pass 1: rho(g) and grad rho(g) for every point -- one launch, one
   ! thread per point, no race (each writes only its own index). The local
   ! basis of the point's batch is walked as a double loop over LOCAL
   ! SHELLS (not a cached per-point tile, which would need an array sized
   ! by the batch's runtime local-function count -- not a legal device-
   ! routine automatic): `shell_collocate` is called once per (u-shell,
   ! v-shell) pair and its result reused across their Cartesian components,
   ! so the redundant work past trc_xc's tiled-GEMM approach is a constant
   ! factor on this already-O(N) pass, not a new scaling.
   !
   subroutine vv10_density(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                           npts, nbatch, r, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                           dmat, rho, grho, sigma)
      integer, intent(in) :: nao, nshell, maxnp, npts, nbatch, nbsh, nbao
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: batch_of(npts), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao)
      real(dp), intent(in) :: dmat(nao, nao)
      real(dp), intent(out) :: rho(npts), grho(3, npts), sigma(npts)
      integer :: g

      do concurrent(g=1:npts)
         call vv10_density_body(g, nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                npts, nbatch, r, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                                dmat, rho(g), grho(1, g), grho(2, g), grho(3, g), sigma(g))
      end do
   end subroutine vv10_density

   pure subroutine vv10_density_body(g, nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                     npts, nbatch, r, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                                     dmat, rho_g, gx, gy, gz, sig_g)
      !$acc routine seq
      integer, intent(in) :: g, nao, nshell, maxnp, npts, nbatch, nbsh, nbao
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: batch_of(npts), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao)
      real(dp), intent(in) :: dmat(nao, nao)
      real(dp), intent(out) :: rho_g, gx, gy, gz, sig_g

      integer :: ib, ku, kv, ishu, ishv, iaou0, iaov0, ncu, ncv, m, n, au, av
      real(dp) :: du(3), dv(3), csu(NCART_MAX), gsu(3, NCART_MAX), csv(NCART_MAX), gsv(3, NCART_MAX)
      real(dp) :: dd

      ib = batch_of(g)
      rho_g = 0.0_dp; gx = 0.0_dp; gy = 0.0_dp; gz = 0.0_dp
      do ku = b_shoff(ib), b_shoff(ib + 1) - 1
         ishu = b_sh(ku)
         iaou0 = b_shao(ku)
         du(1) = r(1, g) - sh_r(1, ishu)
         du(2) = r(2, g) - sh_r(2, ishu)
         du(3) = r(3, g) - sh_r(3, ishu)
         call shell_collocate(sh_l(ishu), sh_np(ishu), sh_e(1:sh_np(ishu), ishu), sh_c(1:sh_np(ishu), ishu), &
                              du, csu, gsu)
         ncu = (sh_l(ishu) + 1)*(sh_l(ishu) + 2)/2
         do kv = b_shoff(ib), b_shoff(ib + 1) - 1
            ishv = b_sh(kv)
            iaov0 = b_shao(kv)
            dv(1) = r(1, g) - sh_r(1, ishv)
            dv(2) = r(2, g) - sh_r(2, ishv)
            dv(3) = r(3, g) - sh_r(3, ishv)
            call shell_collocate(sh_l(ishv), sh_np(ishv), sh_e(1:sh_np(ishv), ishv), sh_c(1:sh_np(ishv), ishv), &
                                 dv, csv, gsv)
            ncv = (sh_l(ishv) + 1)*(sh_l(ishv) + 2)/2
            do m = 1, ncu
               au = b_ao(b_aooff(ib) + iaou0 + m - 1)
               do n = 1, ncv
                  av = b_ao(b_aooff(ib) + iaov0 + n - 1)
                  dd = dmat(av, au)
                  rho_g = rho_g + dd*csu(m)*csv(n)
                  gx = gx + 2.0_dp*dd*gsu(1, m)*csv(n)
                  gy = gy + 2.0_dp*dd*gsu(2, m)*csv(n)
                  gz = gz + 2.0_dp*dd*gsu(3, m)*csv(n)
               end do
            end do
         end do
      end do
      sig_g = gx*gx + gy*gy + gz*gz
   end subroutine vv10_density_body

   !
   ! Pass 2: the O(N^2) double sum -- one launch, `do concurrent` over the
   ! outer point i, a plain sequential loop over j inside it. Each i owns
   ! vrho(i), vsigma(i), eg(i) outright, so no atomics anywhere here,
   ! unlike pass 3.
   !
   subroutine vv10_kernel(npts, r, w, rho, sigma, b, c, kvv, beta, vrho, vsigma, eg)
      integer, intent(in) :: npts
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: b, c, kvv, beta
      real(dp), intent(out) :: vrho(npts), vsigma(npts), eg(npts)
      integer :: i

      do concurrent(i=1:npts)
         call vv10_kernel_body(i, npts, r, w, rho, sigma, c, kvv, beta, vrho(i), vsigma(i), eg(i))
      end do
   end subroutine vv10_kernel

   pure subroutine vv10_kernel_body(i, npts, r, w, rho, sigma, c, kvv, beta, vrho_i, vsigma_i, eg_i)
      !$acc routine seq
      integer, intent(in) :: i, npts
      real(dp), intent(in) :: r(3, npts), w(npts), rho(npts), sigma(npts)
      real(dp), intent(in) :: c, kvv, beta
      real(dp), intent(out) :: vrho_i, vsigma_i, eg_i

      integer :: j
      real(dp) :: ri, gi, w0tmp, w0, dw0dr, dw0dg, kk, dkdr
      real(dp) :: rp, gp, w0ptmp, w0p, kp, rpw, dx, dy, dz, r2
      real(dp) :: gg, ggp, gt, tt, t2, fsum, usum, wsum, ff

      vrho_i = 0.0_dp; vsigma_i = 0.0_dp; eg_i = 0.0_dp
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

      fsum = 0.0_dp; usum = 0.0_dp; wsum = 0.0_dp
      do j = 1, npts
         rp = rho(j)
         if (rp < VV10_RHO_TOL) cycle
         gp = sigma(j)
         w0ptmp = c*(gp/(rp*rp))**2
         w0p = sqrt(w0ptmp + PI43*rp)
         kp = kvv*rp**(1.0_dp/6.0_dp)
         rpw = rp*w(j)
         dx = r(1, i) - r(1, j)
         dy = r(2, i) - r(2, j)
         dz = r(3, i) - r(3, j)
         r2 = dx*dx + dy*dy + dz*dz
         ggp = r2*w0p + kp
         gg = r2*w0 + kk
         gt = gg + ggp
         tt = rpw/(gg*ggp*gt)
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
   ! Pass 3: the contraction, GPU form -- `do concurrent` over points with
   ! an OpenACC atomic on every (u, v) scatter, legal inside a `routine
   ! seq` body. See this module's header for why the CPU form below is a
   ! separate, complete loop rather than this one under an `#ifdef`.
   !
   subroutine vv10_potential_gpu(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                 npts, nbatch, r, w, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                                 vrho, vsigma, grho, vmat)
      integer, intent(in) :: nao, nshell, maxnp, npts, nbatch, nbsh, nbao
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts), w(npts)
      integer, intent(in) :: batch_of(npts), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao)
      real(dp), intent(in) :: vrho(npts), vsigma(npts), grho(3, npts)
      real(dp), intent(inout) :: vmat(nao, nao)
      integer :: g

      do concurrent(g=1:npts)
         call vv10_potential_body(g, nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                  npts, nbatch, r, w, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                                  vrho, vsigma, grho, vmat)
      end do
   end subroutine vv10_potential_gpu

   !
   ! Pass 3, CPU form: a plain sequential loop, since `do concurrent` makes
   ! no atomicity promise a host compiler is bound to honour and this scan
   ! writes the same (u, v) from more than one point.
   !
   subroutine vv10_potential_cpu(nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                 npts, nbatch, r, w, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                                 vrho, vsigma, grho, vmat)
      integer, intent(in) :: nao, nshell, maxnp, npts, nbatch, nbsh, nbao
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts), w(npts)
      integer, intent(in) :: batch_of(npts), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao)
      real(dp), intent(in) :: vrho(npts), vsigma(npts), grho(3, npts)
      real(dp), intent(inout) :: vmat(nao, nao)
      integer :: g

      do g = 1, npts
         call vv10_potential_body(g, nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                  npts, nbatch, r, w, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, b_ao, &
                                  vrho, vsigma, grho, vmat)
      end do
   end subroutine vv10_potential_cpu

   !
   ! One point's contribution to every (u, v) its batch's local basis
   ! reaches: Z_u/Z_v formed from vrho(g)/vsigma(g)/grad rho(g) exactly as
   ! trc_xc's kernel 3 forms them for the semilocal case (nspin = 1), and
   ! scattered with `chi_u Z_v + chi_v Z_u` -- the full ordered double loop
   ! over local shells visits every (u, v) pair exactly once per point, so
   ! this already IS `V_b + V_b^T`; there is no separate transpose step.
   !
   pure subroutine vv10_potential_body(g, nao, nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, &
                                       npts, nbatch, r, w, batch_of, b_shoff, nbsh, b_sh, b_shao, b_aooff, nbao, &
                                       b_ao, vrho, vsigma, grho, vmat)
      !$acc routine seq
      integer, intent(in) :: g, nao, nshell, maxnp, npts, nbatch, nbsh, nbao
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts), w(npts)
      integer, intent(in) :: batch_of(npts), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao)
      real(dp), intent(in) :: vrho(npts), vsigma(npts), grho(3, npts)
      real(dp), intent(inout) :: vmat(nao, nao)

      integer :: ib, ku, kv, ishu, ishv, iaou0, iaov0, ncu, ncv, m, n, au, av
      real(dp) :: du(3), dv(3), csu(NCART_MAX), gsu(3, NCART_MAX), csv(NCART_MAX), gsv(3, NCART_MAX)
      real(dp) :: zg, zvx, zvy, zvz, zu_m, zv_n, contrib

      zg = w(g)*vrho(g)
      zvx = 2.0_dp*w(g)*vsigma(g)*grho(1, g)
      zvy = 2.0_dp*w(g)*vsigma(g)*grho(2, g)
      zvz = 2.0_dp*w(g)*vsigma(g)*grho(3, g)
      if (zg == 0.0_dp .and. zvx == 0.0_dp .and. zvy == 0.0_dp .and. zvz == 0.0_dp) return

      ib = batch_of(g)
      do ku = b_shoff(ib), b_shoff(ib + 1) - 1
         ishu = b_sh(ku)
         iaou0 = b_shao(ku)
         du(1) = r(1, g) - sh_r(1, ishu)
         du(2) = r(2, g) - sh_r(2, ishu)
         du(3) = r(3, g) - sh_r(3, ishu)
         call shell_collocate(sh_l(ishu), sh_np(ishu), sh_e(1:sh_np(ishu), ishu), sh_c(1:sh_np(ishu), ishu), &
                              du, csu, gsu)
         ncu = (sh_l(ishu) + 1)*(sh_l(ishu) + 2)/2
         do kv = b_shoff(ib), b_shoff(ib + 1) - 1
            ishv = b_sh(kv)
            iaov0 = b_shao(kv)
            dv(1) = r(1, g) - sh_r(1, ishv)
            dv(2) = r(2, g) - sh_r(2, ishv)
            dv(3) = r(3, g) - sh_r(3, ishv)
            call shell_collocate(sh_l(ishv), sh_np(ishv), sh_e(1:sh_np(ishv), ishv), sh_c(1:sh_np(ishv), ishv), &
                                 dv, csv, gsv)
            ncv = (sh_l(ishv) + 1)*(sh_l(ishv) + 2)/2
            do m = 1, ncu
               au = b_ao(b_aooff(ib) + iaou0 + m - 1)
               zu_m = 0.5_dp*zg*csu(m) + zvx*gsu(1, m) + zvy*gsu(2, m) + zvz*gsu(3, m)
               do n = 1, ncv
                  av = b_ao(b_aooff(ib) + iaov0 + n - 1)
                  zv_n = 0.5_dp*zg*csv(n) + zvx*gsv(1, n) + zvy*gsv(2, n) + zvz*gsv(3, n)
                  contrib = csu(m)*zv_n + csv(n)*zu_m
                  !$acc atomic update
                  vmat(au, av) = vmat(au, av) + contrib
               end do
            end do
         end do
      end do
   end subroutine vv10_potential_body

end module trc_vv10
