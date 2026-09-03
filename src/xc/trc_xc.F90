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
! CHUNKS
! ------
! The stored values are n_points x n_local per batch, which for a large
! system exceeds any device. Batches are therefore processed in chunks
! whose stored values fit `budget` doubles, and the two kernels run once
! per chunk. Blocking is designed in rather than retrofitted: retrofitting
! it is the coupled-cluster lesson metalquicha already paid for.
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
   use trc_xc_functional, only: trc_xc_functional_t, xc_eval_point, XC_MAX_KERNELS
   use trc_collocation, only: shell_collocate, NCART_MAX
   implicit none
   private

   public :: trc_xc_rks

contains

   !
   ! E_xc, V_xc and the integrated electron count for a closed-shell density.
   !
   ! `budget` caps the stored basis values per chunk, in doubles; the
   ! gradients take three times as much again. The default is 8M doubles,
   ! 256 MB with the gradients.
   !
   subroutine trc_xc_rks(b, xg, func, dmat, vxc, exc, nelec, budget, t_points, t_pairs)
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
      real(dp) :: t0, t1, t2

      integer(kind=8) :: cap, need, nchi
      integer :: b0, b1, ib, nloc, npb, np, p0, npairs
      integer, allocatable :: c_off(:), pr_off(:)
      real(dp), allocatable :: chi(:), gchi(:, :), zg(:), zv(:, :), eg(:), ng(:)
      integer :: kid(XC_MAX_KERNELS)
      real(dp) :: coef(XC_MAX_KERNELS)

      cap = 8388608_8
      if (present(budget)) cap = max(budget, 1_8)
      kid = func%kid
      coef = func%coef
      vxc = 0.0_dp
      exc = 0.0_dp
      nelec = 0.0_dp
      if (present(t_points)) t_points = 0.0_dp
      if (present(t_pairs)) t_pairs = 0.0_dp
      if (xg%npts == 0) return

      allocate (c_off(xg%nbatch + 1), pr_off(xg%nbatch + 1))

      !$acc data copyin(dmat, kid, coef) copy(vxc)
      b0 = 1
      do while (b0 <= xg%nbatch)
         ! The chunk: as many batches as the budget holds, and at least one.
         need = 0
         b1 = b0 - 1
         do ib = b0, xg%nbatch
            nloc = xg%b_aooff(ib + 1) - xg%b_aooff(ib)
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
            nloc = xg%b_aooff(ib + 1) - xg%b_aooff(ib)
            npb = xg%b_off(ib + 1) - xg%b_off(ib)
            c_off(ib + 1) = c_off(ib) + nloc*npb
            pr_off(ib + 1) = pr_off(ib) + nloc*nloc
         end do
         npairs = pr_off(b1 + 1)

         allocate (chi(nchi), gchi(3, nchi), zg(np), zv(3, np), eg(np), ng(np))
         !$acc enter data create(chi, gchi, zg, zv, eg, ng) copyin(c_off, pr_off)

         t0 = wall()
         call xc_points(np, p0, b%maxnp, b%nshell, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                        b%nao, dmat, xg%npts, xg%r, xg%w, xg%batch_of, xg%b_off, xg%b_shoff, &
                        size(xg%b_sh), xg%b_sh, xg%b_aooff, size(xg%b_ao), xg%b_ao, c_off, &
                        xg%nbatch, func%nk, kid, coef, nchi, chi, gchi, zg, zv, eg, ng)

         t1 = wall()
         !$acc update self(eg, ng)
         exc = exc + sum(eg)
         nelec = nelec + sum(ng)

         call xc_potential(npairs, b0, b1, np, p0, b%nao, vxc, xg%b_off, xg%b_aooff, &
                           size(xg%b_ao), xg%b_ao, c_off, pr_off, xg%nbatch, nchi, chi, gchi, zg, zv)

         t2 = wall()
         if (present(t_points)) t_points = t_points + (t1 - t0)
         if (present(t_pairs)) t_pairs = t_pairs + (t2 - t1)

         !$acc exit data delete(chi, gchi, zg, zv, eg, ng, c_off, pr_off)
         deallocate (chi, gchi, zg, zv, eg, ng)
         b0 = b1 + 1
      end do
      !$acc end data
   end subroutine trc_xc_rks

   function wall() result(t)
      real(dp) :: t
      integer(kind=8) :: cc, rate
      call system_clock(cc, rate)
      t = real(cc, dp)/real(rate, dp)
   end function wall

   !
   ! Kernel 1: collocation, density, functional -- one thread per point.
   !
   ! The body is a `routine seq` procedure and the loop only calls it. That
   ! is terco's pattern throughout, and it is not cosmetic: given the body
   ! inline, nvfortran parallelises the inner loops over the local basis
   ! across the threads of a block and hands each POINT a whole block, which
   ! is the wrong mapping by two orders of magnitude. A call is opaque to
   ! that analysis, so the loop maps one point to one thread.
   !
   subroutine xc_points(np, p0, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, nao, dmat, &
                        npts, r, w, batch_of, b_off, b_shoff, nbsh, b_sh, b_aooff, nbao, b_ao, &
                        c_off, nbatch, nk, kid, coef, nchi, chi, gchi, zg, zv, eg, ng)
      integer, intent(in) :: np, p0, maxnp, nshell, nao, npts, nbsh, nbao, nbatch, nk
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: dmat(nao, nao), r(3, npts), w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), c_off(nbatch + 1)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(kind=8), intent(in) :: nchi
      real(dp), intent(out) :: chi(nchi), gchi(3, nchi)
      real(dp), intent(out) :: zg(np), zv(3, np), eg(np), ng(np)
      integer :: g

      do concurrent(g=1:np)
         call xc_point_body(p0 + g - 1, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, nao, dmat, &
                            npts, r, w, batch_of, b_off, b_shoff, nbsh, b_sh, b_aooff, nbao, b_ao, &
                            c_off, nbatch, nk, kid, coef, nchi, chi, gchi, &
                            zg(g), zv(1, g), zv(2, g), zv(3, g), eg(g), ng(g))
      end do
   end subroutine xc_points

   pure subroutine xc_point_body(gg, maxnp, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, nao, dmat, &
                                 npts, r, w, batch_of, b_off, b_shoff, nbsh, b_sh, b_aooff, nbao, b_ao, &
                                 c_off, nbatch, nk, kid, coef, nchi, chi, gchi, zg, zvx, zvy, zvz, eg, ng)
      !$acc routine seq
      integer, intent(in) :: gg, maxnp, nshell, nao, npts, nbsh, nbao, nbatch, nk
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: dmat(nao, nao), r(3, npts), w(npts)
      integer, intent(in) :: batch_of(npts), b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh)
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), c_off(nbatch + 1)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      integer(kind=8), intent(in) :: nchi
      real(dp), intent(inout) :: chi(nchi), gchi(3, nchi)
      real(dp), intent(out) :: zg, zvx, zvy, zvz, eg, ng

      integer :: ib, nloc, off, iao, k, ish, m, nc, u, v, au, av, l
      real(dp) :: d(3), cs(NCART_MAX), gs(3, NCART_MAX), rho, grho(3), su, sigma
      real(dp) :: eps, vrho, vsigma

      ib = batch_of(gg)
      nloc = b_aooff(ib + 1) - b_aooff(ib)
      off = c_off(ib) + (gg - b_off(ib))*nloc
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
      ! rho = chi . D chi, grad rho = 2 (D chi) . grad chi
      rho = 0.0_dp
      grho = 0.0_dp
      do u = 1, nloc
         au = b_ao(b_aooff(ib) + u - 1)
         su = 0.0_dp
         do v = 1, nloc
            av = b_ao(b_aooff(ib) + v - 1)
            su = su + dmat(av, au)*chi(off + v)
         end do
         rho = rho + su*chi(off + u)
         grho(1) = grho(1) + 2.0_dp*su*gchi(1, off + u)
         grho(2) = grho(2) + 2.0_dp*su*gchi(2, off + u)
         grho(3) = grho(3) + 2.0_dp*su*gchi(3, off + u)
      end do
      sigma = grho(1)*grho(1) + grho(2)*grho(2) + grho(3)*grho(3)
      call xc_eval_point(nk, kid, coef, rho, sigma, eps, vrho, vsigma)
      eg = w(gg)*rho*eps
      ng = w(gg)*rho
      zg = w(gg)*vrho
      zvx = 2.0_dp*w(gg)*vsigma*grho(1)
      zvy = 2.0_dp*w(gg)*vsigma*grho(2)
      zvz = 2.0_dp*w(gg)*vsigma*grho(3)
   end subroutine xc_point_body

   !
   ! Kernel 2: V_uv over the chunk -- one thread per (batch, u, v).
   !
   subroutine xc_potential(npairs, b0, b1, np, p0, nao, vxc, b_off, b_aooff, nbao, b_ao, &
                           c_off, pr_off, nbatch, nchi, chi, gchi, zg, zv)
      integer, intent(in) :: npairs, b0, b1, np, p0, nao, nbatch, nbao
      real(dp), intent(inout) :: vxc(nao, nao)
      integer, intent(in) :: b_off(nbatch + 1), b_aooff(nbatch + 1), b_ao(nbao)
      integer, intent(in) :: c_off(nbatch + 1), pr_off(nbatch + 1)
      integer(kind=8), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(3, nchi), zg(np), zv(3, np)
      integer :: t, au, av
      real(dp) :: acc

      do concurrent(t=1:npairs) local(au, av, acc)
         call xc_pair_body(t, b0, b1, np, p0, b_off, b_aooff, nbao, b_ao, c_off, pr_off, nbatch, &
                           nchi, chi, gchi, zg, zv, au, av, acc)
         !$acc atomic update
         vxc(au, av) = vxc(au, av) + acc
      end do
   end subroutine xc_potential

   pure subroutine xc_pair_body(t, b0, b1, np, p0, b_off, b_aooff, nbao, b_ao, c_off, pr_off, nbatch, &
                                nchi, chi, gchi, zg, zv, au, av, acc)
      !$acc routine seq
      integer, intent(in) :: t, b0, b1, np, p0, nbatch, nbao
      integer, intent(in) :: b_off(nbatch + 1), b_aooff(nbatch + 1), b_ao(nbao)
      integer, intent(in) :: c_off(nbatch + 1), pr_off(nbatch + 1)
      integer(kind=8), intent(in) :: nchi
      real(dp), intent(in) :: chi(nchi), gchi(3, nchi), zg(np), zv(3, np)
      integer, intent(out) :: au, av
      real(dp), intent(out) :: acc
      integer :: lo, hi, mid, ib, tt, nloc, u, v, g, gl, col
      real(dp) :: cu, cv

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
      acc = 0.0_dp
      do g = b_off(ib), b_off(ib + 1) - 1
         gl = g - p0 + 1
         col = c_off(ib) + (g - b_off(ib))*nloc
         cu = chi(col + u)
         cv = chi(col + v)
         acc = acc + zg(gl)*cu*cv &
               + zv(1, gl)*(gchi(1, col + u)*cv + cu*gchi(1, col + v)) &
               + zv(2, gl)*(gchi(2, col + u)*cv + cu*gchi(2, col + v)) &
               + zv(3, gl)*(gchi(3, col + u)*cv + cu*gchi(3, col + v))
      end do
      au = b_ao(b_aooff(ib) + u - 1)
      av = b_ao(b_aooff(ib) + v - 1)
   end subroutine xc_pair_body

end module trc_xc
