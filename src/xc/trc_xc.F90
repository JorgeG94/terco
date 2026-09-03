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
! TWO KERNELS PER CHUNK
! ---------------------
! Kernel 1 is one thread per grid point. It collocates the batch's local
! basis at its point, stores the values and gradients, contracts them with
! D into rho and grad rho, and evaluates the functional -- so the density
! and every per-point quantity leave the kernel already weighted:
! z = w v_rho, zv = 2 w v_sigma grad rho, e = w rho eps.
!
! Kernel 2 is one thread per (batch, u, v). It sums over the batch's points
! reading the stored values, and adds its one number into V_uv. That is the
! GEMM V = chi^T Z chi written as a loop nest, which is what the thesis
! asks for: no library on the path. Distinct batches share global AOs, so
! the add is atomic -- one atomic per (batch, u, v), not per point.
!
! The loop bodies live in `routine seq` procedures and the loops only call
! them. That is terco's pattern throughout, and it is not cosmetic: given
! the body inline, nvfortran parallelises the inner loops over the local
! basis across the threads of a block and hands each POINT a whole block,
! which is the wrong mapping by two orders of magnitude.
!
! CHUNKS
! ------
! The stored values are n_points x n_local per batch, which for a large
! system exceeds any device. Batches are therefore processed in chunks
! whose stored values fit `budget` doubles, and the two kernels run once
! per chunk. Blocking is designed in rather than retrofitted: retrofitting
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
! WHAT RUNS ON THE HOST
! ---------------------
! The sums for E_xc and the electron count. `do concurrent` has no portable
! reduction -- gfortran took `reduce` after locality specifiers, and terco
! builds on the host with what a CI runner has -- so the per-point terms
! come back as arrays and are summed here. They are one double per point;
! the traffic is nothing next to the values.
!
module trc_xc
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_eval_point, xc_eval_point_polar, XC_MAX_KERNELS
   use trc_collocation, only: shell_collocate, NCART_MAX
   implicit none
   private

   public :: trc_xc_rks, trc_xc_uks, XC_RHO_TOL

   !> A batch whose bound on rho is below this is skipped. The functionals'
   !> own density thresholds are 1e-12 to 1e-32; check_xc_energy measures
   !> what this value costs against an unscreened evaluation.
   real(dp), parameter :: XC_RHO_TOL = 1.0e-12_dp

contains

   !
   ! E_xc, V_xc and the integrated electron count for a closed-shell density.
   !
   ! `budget` caps the stored basis values per chunk, in doubles; the
   ! gradients take three times as much again. The default is 8M doubles,
   ! 256 MB with the gradients.
   !
   subroutine trc_xc_rks(b, xg, func, dmat, vxc, exc, nelec, budget, t_points, t_pairs, rho_tol, n_skipped)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      real(dp), intent(in) :: dmat(b%nao, b%nao)
      real(dp), intent(out) :: vxc(b%nao, b%nao)
      real(dp), intent(out) :: exc, nelec
      integer(kind=8), intent(in), optional :: budget
      !! Wall time spent in the two kernels, for the benchmark. A `do
      !! concurrent` under -stdpar returns when the device is done, so the
      !! clock around it is the kernel's.
      real(dp), intent(out), optional :: t_points, t_pairs
      !! Density screening threshold on the per-batch bound of rho; zero
      !! screens nothing. Default XC_RHO_TOL.
      real(dp), intent(in), optional :: rho_tol
      integer, intent(out), optional :: n_skipped   !! batches the screening removed
      real(dp) :: tp, tq, rtol
      integer(kind=8) :: cap
      integer :: nsk

      cap = 8388608_8
      if (present(budget)) cap = max(budget, 1_8)
      rtol = XC_RHO_TOL
      if (present(rho_tol)) rtol = rho_tol
      call xc_drive(b, xg, func, 1, dmat, vxc, exc, nelec, cap, rtol, tp, tq, nsk)
      if (present(t_points)) t_points = tp
      if (present(t_pairs)) t_pairs = tq
      if (present(n_skipped)) n_skipped = nsk
   end subroutine trc_xc_rks

   !
   ! The same for a spin-polarised density: D(:,:,1) is alpha, D(:,:,2)
   ! beta, and the potentials come back in the same slots.
   !
   subroutine trc_xc_uks(b, xg, func, dmat, vxc, exc, nelec, budget, t_points, t_pairs, rho_tol, n_skipped)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      real(dp), intent(in) :: dmat(b%nao, b%nao, 2)
      real(dp), intent(out) :: vxc(b%nao, b%nao, 2)
      real(dp), intent(out) :: exc, nelec
      integer(kind=8), intent(in), optional :: budget
      real(dp), intent(out), optional :: t_points, t_pairs
      real(dp), intent(in), optional :: rho_tol
      integer, intent(out), optional :: n_skipped
      real(dp) :: tp, tq, rtol
      integer(kind=8) :: cap
      integer :: nsk

      cap = 8388608_8
      if (present(budget)) cap = max(budget, 1_8)
      rtol = XC_RHO_TOL
      if (present(rho_tol)) rtol = rho_tol
      call xc_drive(b, xg, func, 2, dmat, vxc, exc, nelec, cap, rtol, tp, tq, nsk)
      if (present(t_points)) t_points = tp
      if (present(t_pairs)) t_pairs = tq
      if (present(n_skipped)) n_skipped = nsk
   end subroutine trc_xc_uks

   subroutine xc_drive(b, xg, func, nspin, dmat, vxc, exc, nelec, cap, rho_tol, t_points, t_pairs, n_skipped)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_functional_t), intent(in) :: func
      integer, intent(in) :: nspin
      real(dp), intent(in) :: dmat(b%nao, b%nao, nspin)
      real(dp), intent(out) :: vxc(b%nao, b%nao, nspin)
      real(dp), intent(out) :: exc, nelec, t_points, t_pairs
      integer(kind=8), intent(in) :: cap
      real(dp), intent(in) :: rho_tol
      integer, intent(out) :: n_skipped

      integer(kind=8) :: need, nchi
      integer :: b0, b1, ib, nloc, npb, np, p0, npairs, u, v, s, au, av
      integer, allocatable :: c_off(:), pr_off(:), active(:)
      real(dp) :: bound, su
      real(dp), allocatable :: chi(:), gchi(:, :), zg(:, :), zv(:, :, :), eg(:), ng(:)
      integer :: kid(XC_MAX_KERNELS)
      real(dp) :: coef(XC_MAX_KERNELS)
      real(dp) :: t0, t1, t2

      kid = func%kid
      coef = func%coef
      vxc = 0.0_dp
      exc = 0.0_dp
      nelec = 0.0_dp
      t_points = 0.0_dp
      t_pairs = 0.0_dp
      if (xg%npts == 0) return

      allocate (c_off(xg%nbatch + 1), pr_off(xg%nbatch + 1), active(xg%nbatch))

      ! Density screening: bound rho on each batch and drop what cannot reach.
      n_skipped = 0
      do ib = 1, xg%nbatch
         bound = 0.0_dp
         do s = 1, nspin
            do u = xg%b_aooff(ib), xg%b_aooff(ib + 1) - 1
               au = xg%b_ao(u)
               su = 0.0_dp
               do v = xg%b_aooff(ib), xg%b_aooff(ib + 1) - 1
                  av = xg%b_ao(v)
                  su = su + abs(dmat(av, au, s))*xg%b_amax(v)
               end do
               bound = bound + su*xg%b_amax(u)
            end do
         end do
         if (bound > rho_tol) then
            active(ib) = 1
         else
            active(ib) = 0
            n_skipped = n_skipped + 1
         end if
      end do

      !$acc data copyin(dmat, kid, coef, active) copy(vxc)
      b0 = 1
      do while (b0 <= xg%nbatch)
         ! The chunk: as many batches as the budget holds, and at least one.
         need = 0
         b1 = b0 - 1
         do ib = b0, xg%nbatch
            nloc = (xg%b_aooff(ib + 1) - xg%b_aooff(ib))*active(ib)
            npb = xg%b_off(ib + 1) - xg%b_off(ib)
            if (need + int(nloc, 8)*int(npb, 8) > cap .and. ib > b0) exit
            need = need + int(nloc, 8)*int(npb, 8)
            b1 = ib
         end do
         nchi = max(need, 1_8)
         p0 = xg%b_off(b0)
         np = xg%b_off(b1 + 1) - p0
         ! Offsets within the chunk: stored values, and (u,v) pairs.
         c_off(b0) = 0
         pr_off(b0) = 0
         do ib = b0, b1
            nloc = (xg%b_aooff(ib + 1) - xg%b_aooff(ib))*active(ib)
            npb = xg%b_off(ib + 1) - xg%b_off(ib)
            c_off(ib + 1) = c_off(ib) + nloc*npb
            pr_off(ib + 1) = pr_off(ib) + nloc*nloc
         end do
         npairs = pr_off(b1 + 1)

         allocate (chi(nchi), gchi(3, nchi), zg(nspin, np), zv(3, nspin, np), eg(np), ng(np))
         !$acc enter data create(chi, gchi, zg, zv, eg, ng) copyin(c_off, pr_off)

         t0 = wall()
         call xc_points(np, p0, nspin, b%maxnp, b%nshell, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                        b%nao, dmat, xg%npts, xg%r, xg%w, xg%batch_of, xg%b_off, xg%b_shoff, &
                        size(xg%b_sh), xg%b_sh, xg%b_aooff, size(xg%b_ao), xg%b_ao, c_off, active, &
                        xg%nbatch, func%nk, kid, coef, nchi, chi, gchi, zg, zv, eg, ng)
         t1 = wall()

         !$acc update self(eg, ng)
         exc = exc + sum(eg)
         nelec = nelec + sum(ng)

         call xc_potential(npairs, b0, b1, np, p0, nspin, b%nao, vxc, xg%b_off, xg%b_aooff, &
                           size(xg%b_ao), xg%b_ao, c_off, pr_off, xg%nbatch, nchi, chi, gchi, zg, zv)
         t2 = wall()
         t_points = t_points + (t1 - t0)
         t_pairs = t_pairs + (t2 - t1)

         !$acc exit data delete(chi, gchi, zg, zv, eg, ng, c_off, pr_off)
         deallocate (chi, gchi, zg, zv, eg, ng)
         b0 = b1 + 1
      end do
      !$acc end data
   end subroutine xc_drive

   function wall() result(t)
      real(dp) :: t
      integer(kind=8) :: cc, rate
      call system_clock(cc, rate)
      t = real(cc, dp)/real(rate, dp)
   end function wall

   !
   ! Kernel 1: collocation, density, functional -- one thread per point.
   !
   subroutine xc_points(np, p0, nspin, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, nao, dmat, &
                        npts, r, w, batch_of, b_off, b_shoff, nbsh, b_sh, b_aooff, nbao, b_ao, &
                        c_off, active, nbatch, nk, kid, coef, nchi, chi, gchi, zg, zv, eg, ng)
      integer, intent(in) :: np, p0, nspin, maxnp, nshell, nao, npts, nbsh, nbao, nbatch, nk
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: dmat(nao, nao, nspin), r(3, npts), w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), c_off(nbatch + 1), active(nbatch)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(kind=8), intent(in) :: nchi
      real(dp), intent(out) :: chi(nchi), gchi(3, nchi)
      real(dp), intent(out) :: zg(nspin, np), zv(3, nspin, np), eg(np), ng(np)
      integer :: g

      ! Nothing but the call in the body: a branch with array syntax here
      ! cost the point kernel a factor of ten on the device, because the
      ! compiler took the zeroing as an inner loop to parallelise.
      do concurrent(g=1:np)
         call xc_point_body(p0 + g - 1, nspin, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, nao, dmat, &
                            npts, r, w, batch_of, b_off, b_shoff, nbsh, b_sh, b_aooff, nbao, b_ao, &
                            c_off, active, nbatch, nk, kid, coef, nchi, chi, gchi, &
                            zg(1, g), zv(1, 1, g), eg(g), ng(g))
      end do
   end subroutine xc_points

   pure subroutine xc_point_body(gg, nspin, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, nao, dmat, &
                                 npts, r, w, batch_of, b_off, b_shoff, nbsh, b_sh, b_aooff, nbao, b_ao, &
                                 c_off, active, nbatch, nk, kid, coef, nchi, chi, gchi, zg, zv, eg, ng)
      !$acc routine seq
      integer, intent(in) :: gg, nspin, maxnp, nshell, nao, npts, nbsh, nbao, nbatch, nk
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: dmat(nao, nao, nspin), r(3, npts), w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), c_off(nbatch + 1), active(nbatch)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(kind=8), intent(in) :: nchi
      real(dp), intent(inout) :: chi(nchi), gchi(3, nchi)
      real(dp), intent(out) :: zg(nspin), zv(3, nspin), eg, ng

      integer :: ib, nloc, off, iao, k, ish, m, nc, u, v, au, av, l, s
      real(dp) :: d(3), cs(NCART_MAX), gs(3, NCART_MAX), rho(2), grho(3, 2), su, wg
      real(dp) :: sigma, eps, vrho, vsigma, s_aa, s_ab, s_bb, va, vb, v_aa, v_ab, v_bb

      ib = batch_of(gg)
      eg = 0.0_dp
      ng = 0.0_dp
      do s = 1, nspin
         zg(s) = 0.0_dp
         zv(1, s) = 0.0_dp; zv(2, s) = 0.0_dp; zv(3, s) = 0.0_dp
      end do
      if (active(ib) == 0) return
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      off = c_off(ib) + (gg - b_off(ib))*nloc
      wg = w(gg)
      ! Collocate the local basis at this point.
      iao = 0
      do k = b_shoff(ib), b_shoff(ib + 1) - 1
         ish = b_sh(k)
         l = sh_l(ish)
         d(1) = r(1, gg) - sh_r(1, ish)
         d(2) = r(2, gg) - sh_r(2, ish)
         d(3) = r(3, gg) - sh_r(3, ish)
         call shell_collocate(l, sh_np(ish), sh_e(1:sh_np(ish), ish), sh_c(1:sh_np(ish), ish), d, cs, gs)
         nc = (l + 1)*(l + 2)/2
         do m = 1, nc
            chi(off + iao + m) = cs(m)
            gchi(1, off + iao + m) = gs(1, m)
            gchi(2, off + iao + m) = gs(2, m)
            gchi(3, off + iao + m) = gs(3, m)
         end do
         iao = iao + nc
      end do
      ! Per spin: rho = chi . D chi, grad rho = 2 (D chi) . grad chi
      rho = 0.0_dp
      grho = 0.0_dp
      do s = 1, nspin
         do u = 1, nloc
            au = b_ao(b_aooff(ib) + u - 1)
            su = 0.0_dp
            do v = 1, nloc
               av = b_ao(b_aooff(ib) + v - 1)
               su = su + dmat(av, au, s)*chi(off + v)
            end do
            rho(s) = rho(s) + su*chi(off + u)
            grho(1, s) = grho(1, s) + 2.0_dp*su*gchi(1, off + u)
            grho(2, s) = grho(2, s) + 2.0_dp*su*gchi(2, off + u)
            grho(3, s) = grho(3, s) + 2.0_dp*su*gchi(3, off + u)
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
   end subroutine xc_point_body

   !
   ! Kernel 2: V_uv over the chunk -- one thread per (batch, u, v).
   !
   subroutine xc_potential(npairs, b0, b1, np, p0, nspin, nao, vxc, b_off, b_aooff, nbao, b_ao, &
                           c_off, pr_off, nbatch, nchi, chi, gchi, zg, zv)
      integer, intent(in) :: npairs, b0, b1, np, p0, nspin, nao, nbatch, nbao
      real(dp), intent(inout) :: vxc(nao, nao, nspin)
      integer, intent(in) :: b_off(nbatch + 1), b_aooff(nbatch + 1), b_ao(nbao)
      integer, intent(in) :: c_off(nbatch + 1), pr_off(nbatch + 1)
      integer(kind=8), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(3, nchi), zg(nspin, np), zv(3, nspin, np)
      integer :: t, au, av
      real(dp) :: acc_a, acc_b

      do concurrent(t=1:npairs) local(au, av, acc_a, acc_b)
         call xc_pair_body(t, b0, b1, np, p0, nspin, b_off, b_aooff, nbao, b_ao, c_off, pr_off, nbatch, &
                           nchi, chi, gchi, zg, zv, au, av, acc_a, acc_b)
         !$acc atomic update
         vxc(au, av, 1) = vxc(au, av, 1) + acc_a
         if (nspin == 2) then
            !$acc atomic update
            vxc(au, av, 2) = vxc(au, av, 2) + acc_b
         end if
      end do
   end subroutine xc_potential

   pure subroutine xc_pair_body(t, b0, b1, np, p0, nspin, b_off, b_aooff, nbao, b_ao, c_off, pr_off, nbatch, &
                                nchi, chi, gchi, zg, zv, au, av, acc_a, acc_b)
      !$acc routine seq
      integer, intent(in) :: t, b0, b1, np, p0, nspin, nbatch, nbao
      integer, intent(in) :: b_off(nbatch + 1), b_aooff(nbatch + 1), b_ao(nbao)
      integer, intent(in) :: c_off(nbatch + 1), pr_off(nbatch + 1)
      integer(kind=8), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(3, nchi), zg(nspin, np), zv(3, nspin, np)
      integer, intent(out) :: au, av
      real(dp), intent(out) :: acc_a, acc_b
      integer :: lo, hi, mid, ib, tt, nloc, u, v, g, gl, col
      real(dp) :: cu, cv, sx, sy, sz

      ! Which batch owns pair t: the last ib in [b0, b1] with pr_off(ib) < t.
      lo = b0; hi = b1
      do while (lo < hi)
         mid = (lo + hi + 1)/2
         if (pr_off(mid) < t) then
            lo = mid
         else
            hi = mid - 1
         end if
      end do
      ib = lo
      tt = t - pr_off(ib) - 1
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      u = mod(tt, nloc) + 1
      v = tt/nloc + 1
      acc_a = 0.0_dp
      acc_b = 0.0_dp
      do g = b_off(ib), b_off(ib + 1) - 1
         gl = g - p0 + 1
         col = c_off(ib) + (g - b_off(ib))*nloc
         cu = chi(col + u)
         cv = chi(col + v)
         sx = gchi(1, col + u)*cv + cu*gchi(1, col + v)
         sy = gchi(2, col + u)*cv + cu*gchi(2, col + v)
         sz = gchi(3, col + u)*cv + cu*gchi(3, col + v)
         acc_a = acc_a + zg(1, gl)*cu*cv + zv(1, 1, gl)*sx + zv(2, 1, gl)*sy + zv(3, 1, gl)*sz
         if (nspin == 2) acc_b = acc_b + zg(2, gl)*cu*cv + zv(1, 2, gl)*sx + zv(2, 2, gl)*sy + zv(3, 2, gl)*sz
      end do
      au = b_ao(b_aooff(ib) + u - 1)
      av = b_ao(b_aooff(ib) + v - 1)
   end subroutine xc_pair_body

end module trc_xc
