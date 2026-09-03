!
! Basis functions on grid points: the two checks that fix the conventions.
!
! The collocation kernel is a few dozen lines, and the whole difficulty is
! that its output must be the SAME functions terco's integrals are over --
! same normalisation, same Cartesian order, same everything. A wrong
! convention does not fail: the density is assembled in one basis and the
! Fock matrix in another, the SCF converges, and the energy is wrong by an
! amount that looks like a grid or a functional. So two checks, both exact:
!
!   1. ELEMENTWISE AGAINST PYSCF. The values and gradients at a cloud of
!      points are written to `ao_probe.bin`; `ao_ref.py` evaluates the same
!      system with pyscf's eval_ao (cart=True, deriv=1) and compares to
!      1e-12. Direct comparison of the quantity itself, localised to one
!      component of one shell. The system is the uncontracted one from
!      check_mult, so pyscf has nothing to renormalise.
!
!   2. THE NUMERICAL OVERLAP AGAINST THE ANALYTIC ONE, entirely inside
!      terco. sum_g w_g chi_u chi_v over a real molecular grid must match
!      trc_1e's S to the grid's own accuracy, so it validates the values,
!      the grid weights and the Becke partition together, which the pyscf
!      comparison does not touch.
!
!      The grid is level 5 UNPRUNED, and the tolerance is measured, not
!      chosen: on this system pyscf's own grids give max|dS| of 6.4e-4 at
!      level 3 pruned (which is what the first version of this test used,
!      and it failed), 1.2e-5 at level 3 unpruned, 9.2e-8 at level 5
!      unpruned. NWChem pruning fixes the inner spheres at 50 and 86
!      points, which cannot integrate a d-on-d product across centres 1.3
!      Bohr apart. terco's grid at the same settings gives the same numbers
!      to the digit, which is the finding: the grid is pyscf's, and so is
!      its limitation on this probe.
!
! And one control: the overlap through the BATCHED, SCREENED grid must
! agree with the unscreened sum to the screening tolerance. That is the
! only test of the batching, and it is a sharp one, because a shell
! wrongly dropped from a batch changes S by the size of that shell.
!
program check_ao
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_error, only: error_t
   use trc_dft_grid, only: dft_grid_t, build_dft_grid
   use trc_dft_prune, only: PRUNE_NONE
   use trc_collocation, only: shell_collocate, NCART_MAX
   use trc_xc_batch, only: trc_xc_grid_t
   implicit none

   integer, parameter :: NSH = 9, NATM = 3, MAXNP = 1, NPROBE = 200
   integer :: sh_l(NSH), sh_np(NSH), i, k, ish, m, nc, u, v, nao, ib, g, iao, au, av
   real(dp) :: sh_e(MAXNP, NSH), sh_c(MAXNP, NSH), sh_r(3, NSH)
   real(dp) :: at_z(NATM), at_r(3, NATM)
   real(dp) :: pt(3), d(3), cs(NCART_MAX), gs(3, NCART_MAX), seed
   integer(kind=8) :: lcg
   real(dp), allocatable :: probe(:, :), chi(:, :), gchi(:, :, :)
   real(dp), allocatable :: s(:, :), t(:, :), vn(:, :), s_num(:, :), s_bat(:, :), col(:)
   type(trc_basis_t) :: b
   type(trc_pairlist_t) :: pl
   type(dft_grid_t) :: grid
   type(trc_xc_grid_t) :: xg
   type(error_t) :: err
   real(dp) :: err_grid, err_batch, err_trace
   integer :: unit, nfail

   ! --- the probe system, exactly check_mult's ------------------------------
   at_r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
   at_r(:, 2) = [1.3_dp, 0.0_dp, 0.7_dp]
   at_r(:, 3) = [-0.4_dp, 1.1_dp, 1.9_dp]
   at_z = [4.0_dp, 3.0_dp, 2.0_dp]
   do i = 1, NATM
      sh_l(3*i - 2) = 0; sh_e(1, 3*i - 2) = 0.9_dp + 0.17_dp*real(i, dp)
      sh_l(3*i - 1) = 1; sh_e(1, 3*i - 1) = 1.4_dp + 0.11_dp*real(i, dp)
      sh_l(3*i) = 2; sh_e(1, 3*i) = 1.9_dp + 0.13_dp*real(i, dp)
      sh_r(:, 3*i - 2) = at_r(:, i)
      sh_r(:, 3*i - 1) = at_r(:, i)
      sh_r(:, 3*i) = at_r(:, i)
   end do
   sh_np = 1
   do i = 1, NSH
      sh_c(1, i) = gto_norm(sh_l(i), sh_e(1, i))
   end do
   call b%build(NSH, sh_l, sh_np, sh_e, sh_c, sh_r, NATM, at_z, at_r, MAXNP)
   nao = b%nao
   nfail = 0

   ! --- 1. values and gradients at a cloud of points, for pyscf ------------
   allocate (probe(3, NPROBE), chi(NPROBE, nao), gchi(3, NPROBE, nao))
   lcg = 37_8
   do k = 1, NPROBE
      do i = 1, 3
         lcg = mod(25214903917_8*lcg + 11_8, 281474976710656_8)
         seed = real(lcg, dp)/281474976710656.0_dp
         probe(i, k) = -3.0_dp + 7.0_dp*seed
      end do
   end do
   do k = 1, NPROBE
      do ish = 1, NSH
         d = probe(:, k) - b%sh_r(:, ish)
         call shell_collocate(b%sh_l(ish), b%sh_np(ish), b%sh_e(1:1, ish), b%sh_c(1:1, ish), d, cs, gs)
         nc = (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2
         do m = 1, nc
            chi(k, b%sh_ao(ish) + m - 1) = cs(m)
            gchi(:, k, b%sh_ao(ish) + m - 1) = gs(:, m)
         end do
      end do
   end do
   open (newunit=unit, file='ao_probe.bin', access='stream', form='unformatted', status='replace')
   write (unit) int(NPROBE, kind=8), int(nao, kind=8)
   write (unit) probe
   write (unit) chi
   write (unit) gchi
   close (unit)
   print '(a,i0,a,i0,a)', "check_ao: wrote ao_probe.bin (", NPROBE, " points, ", nao, " AOs) for ao_ref.py"

   ! --- 2. the numerical overlap on a real grid -----------------------------
   call pl%build(b, 1.0e-30_dp)
   allocate (s(nao, nao), t(nao, nao), vn(nao, nao))
   call b%to_device(); call pl%to_device()
   call trc_1e(b, pl, s, t, vn)

   call build_dft_grid(at_r, nint(at_z), grid, err, level=5, prune=PRUNE_NONE)
   if (err%has_error()) then
      print '(a)', "check_ao: grid failed: "//err%get_message()
      stop 1
   end if

   ! Unscreened: every AO at every point.
   allocate (s_num(nao, nao), s_bat(nao, nao), col(nao))
   s_num = 0.0_dp
   do g = 1, grid%n_points
      do ish = 1, NSH
         d = grid%coords(:, g) - b%sh_r(:, ish)
         call shell_collocate(b%sh_l(ish), b%sh_np(ish), b%sh_e(1:1, ish), b%sh_c(1:1, ish), d, cs, gs)
         nc = (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2
         col(b%sh_ao(ish):b%sh_ao(ish) + nc - 1) = cs(1:nc)
      end do
      do v = 1, nao
         do u = 1, nao
            s_num(u, v) = s_num(u, v) + grid%weights(g)*col(u)*col(v)
         end do
      end do
   end do

   ! Through the batched, screened container: only each batch's local basis.
   call xg%build(grid%n_points, grid%coords, grid%weights, b, max_pts=256, tol=1.0e-10_dp)
   s_bat = 0.0_dp
   do ib = 1, xg%nbatch
      do g = xg%b_off(ib), xg%b_off(ib + 1) - 1
         iao = 0
         do k = xg%b_shoff(ib), xg%b_shoff(ib + 1) - 1
            ish = xg%b_sh(k)
            d = xg%r(:, g) - b%sh_r(:, ish)
            call shell_collocate(b%sh_l(ish), b%sh_np(ish), b%sh_e(1:1, ish), b%sh_c(1:1, ish), d, cs, gs)
            nc = (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2
            col(iao + 1:iao + nc) = cs(1:nc)
            iao = iao + nc
         end do
         do v = 1, iao
            av = xg%b_ao(xg%b_aooff(ib) + v - 1)
            do u = 1, iao
               au = xg%b_ao(xg%b_aooff(ib) + u - 1)
               s_bat(au, av) = s_bat(au, av) + xg%w(g)*col(u)*col(v)
            end do
         end do
      end do
   end do

   err_grid = maxval(abs(s_num - s))
   err_trace = 0.0_dp
   do u = 1, nao
      err_trace = err_trace + s_num(u, u) - s(u, u)
   end do
   err_batch = maxval(abs(s_bat - s_num))

   print '(a,i0,a,i0,a,i0)', "check_ao: grid ", grid%n_points, " points in ", xg%nbatch, &
      " batches, largest local basis ", xg%max_nloc
   ! pyscf's level-5 unpruned grid on this system: 9.2e-8 and 1.8e-7.
   call report(err_grid < 1.0e-6_dp, "numerical overlap vs analytic, max|dS|", err_grid)
   call report(abs(err_trace) < 1.0e-6_dp, "numerical overlap vs analytic, trace", err_trace)
   call report(err_batch < 1.0e-9_dp, "batched/screened vs unscreened, max|dS|", err_batch)

   call xg%release(); call grid%destroy(); call pl%release(); call b%release()
   if (nfail > 0) then
      print '(a,i0,a)', "check_ao: ", nfail, " check(s) FAILED"
      stop 1
   end if
   print '(a)', "check_ao: internal checks passed; run ao_ref.py for the pyscf comparison"

contains

   subroutine report(ok, what, value)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: value
      if (.not. ok) nfail = nfail + 1
      print '(a,1x,a,t56,es10.2)', merge("ok  ", "FAIL", ok), what, value
   end subroutine report

   pure real(dp) function gto_norm(l, a)
      integer, intent(in) :: l
      real(dp), intent(in) :: a
      real(dp) :: tt
      tt = 2.0_dp**(2*l + 3)*fact(l + 1)*(2.0_dp*a)**(l + 1.5_dp) &
           /(fact(2*l + 2)*sqrt(3.14159265358979323846_dp))
      gto_norm = sqrt(tt)
   end function gto_norm

   pure real(dp) function fact(n)
      integer, intent(in) :: n
      integer :: j
      fact = 1.0_dp
      do j = 2, n
         fact = fact*real(j, dp)
      end do
   end function fact

end program check_ao
