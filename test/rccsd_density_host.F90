module rccsd_density_host
   !!
   !! The CCSD correlation densities by host loops with every block of the two-particle density formed, the v^4 one
   !! included: pyscf/cc/ccsd_rdm.py `_gamma1_intermediates`, `_gamma2_outcore` and `_make_rdm2` (Apache-2.0,
   !! Qiming Sun), then the eightfold symmetrisation test/pyscf_refs/generate.py's gen_ccsd_grad applies. It is the
   !! dense side of check_rccsd_density, for a handful of orbitals only: it shares no formula grouping with
   !! trc_rccsd_density, which never forms the o v^3 and v^4 blocks.
   !!
   use pic_types, only: dp
   implicit none
   private
   public :: dense_densities

contains

   subroutine dense_densities(no, nv, t1, t2, l1, l2, gam, big)
      !! gam(p,q) and big(p,q,r,s), active orbitals with the occupied first, symmetrised
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv), l1(no, nv), l2(no, no, nv, nv)
      real(dp), allocatable, intent(out) :: gam(:, :), big(:, :, :, :)

      real(dp), allocatable :: pvoov_a(:, :, :, :), pvoov_b(:, :, :, :), moo(:, :), mvv(:, :)
      integer :: n

      n = no + nv
      allocate (gam(n, n), big(n, n, n, n))
      gam = 0.0_dp
      big = 0.0_dp
      call one_particle(no, nv, t1, t2, l1, l2, gam)
      call pv_arrays(no, nv, t2, l2, pvoov_a, pvoov_b, moo, mvv)
      call occupied_blocks(no, nv, t1, t2, l1, l2, pvoov_a, pvoov_b, moo, big)
      call ovov_block(no, nv, t1, t2, l1, l2, pvoov_a, pvoov_b, moo, mvv, big)
      call ovvo_blocks(no, nv, t1, l1, l2, pvoov_a, pvoov_b, big)
      call virtual_blocks(no, nv, t1, t2, l1, l2, pvoov_a, pvoov_b, mvv, big)
      call symmetrise(n, big)
   end subroutine dense_densities

   subroutine one_particle(no, nv, t1, t2, l1, l2, gam)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv), l1(no, nv), l2(no, no, nv, nv)
      real(dp), intent(inout) :: gam(no + nv, no + nv)

      real(dp) :: doo(no, no), dvv(nv, nv), dvo(nv, no), xtv(no, no), xt1(no, no), xt2(nv, nv)
      real(dp) :: theta(no, no, nv, nv)
      integer :: i, j, k, m, n, a, b, e, f, c

      do b = 1, nv
         do a = 1, nv
            theta(:, :, a, b) = 2.0_dp*t2(:, :, a, b) - t2(:, :, b, a)
         end do
      end do
      doo = 0.0_dp
      dvv = 0.0_dp
      xtv = 0.0_dp
      xt1 = 0.0_dp
      xt2 = 0.0_dp
      do j = 1, no
         do i = 1, no
            do a = 1, nv
               doo(i, j) = doo(i, j) - t1(j, a)*l1(i, a)
               do b = 1, nv
                  do k = 1, no
                     doo(i, j) = doo(i, j) - theta(j, k, a, b)*l2(i, k, a, b)
                  end do
               end do
            end do
         end do
      end do
      do b = 1, nv
         do a = 1, nv
            do i = 1, no
               dvv(a, b) = dvv(a, b) + t1(i, a)*l1(i, b)
               do j = 1, no
                  do c = 1, nv
                     dvv(a, b) = dvv(a, b) + theta(j, i, c, a)*l2(j, i, c, b)
                  end do
               end do
            end do
         end do
      end do
      do m = 1, no
         do i = 1, no
            do e = 1, nv
               xtv(i, m) = xtv(i, m) + t1(i, e)*l1(m, e)
            end do
            do f = 1, nv
               do e = 1, nv
                  do n = 1, no
                     xt1(m, i) = xt1(m, i) + l2(m, n, e, f)*theta(i, n, e, f)
                  end do
               end do
            end do
         end do
      end do
      do a = 1, nv
         do e = 1, nv
            do f = 1, nv
               do n = 1, no
                  do m = 1, no
                     xt2(e, a) = xt2(e, a) + l2(m, n, a, f)*theta(m, n, e, f)
                  end do
               end do
            end do
         end do
      end do
      do i = 1, no
         do a = 1, nv
            dvo(a, i) = t1(i, a)
            do m = 1, no
               dvo(a, i) = dvo(a, i) - xtv(i, m)*t1(m, a) - xt1(m, i)*t1(m, a)
               do e = 1, nv
                  dvo(a, i) = dvo(a, i) + theta(i, m, a, e)*l1(m, e)
               end do
            end do
            do e = 1, nv
               dvo(a, i) = dvo(a, i) - t1(i, e)*xt2(a, e)
            end do
         end do
      end do
      gam(1:no, 1:no) = doo + transpose(doo)
      gam(1:no, no + 1:) = l1 + transpose(dvo)
      gam(no + 1:, 1:no) = transpose(gam(1:no, no + 1:))
      gam(no + 1:, no + 1:) = dvv + transpose(dvv)
   end subroutine one_particle

   subroutine pv_arrays(no, nv, t2, l2, pa, pb, moo, mvv)
      !! pyscf pvOOv(a,i,j,b) and pvoOV(a,i,j,b), moo, mvv
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(no, no, nv, nv), l2(no, no, nv, nv)
      real(dp), allocatable, intent(out) :: pa(:, :, :, :), pb(:, :, :, :), moo(:, :), mvv(:, :)

      integer :: i, j, k, l, a, b, c, d

      allocate (pa(nv, no, no, nv), pb(nv, no, no, nv), moo(no, no), mvv(nv, nv))
      pa = 0.0_dp
      pb = 0.0_dp
      do b = 1, nv
         do j = 1, no
            do i = 1, no
               do a = 1, nv
                  do c = 1, nv
                     do k = 1, no
                        pa(a, i, j, b) = pa(a, i, j, b) + l2(i, k, c, a)*t2(j, k, c, b)
                        pb(a, i, j, b) = pb(a, i, j, b) - l2(i, k, c, a)*t2(j, k, b, c) &
                                         + l2(i, k, a, c)*(2.0_dp*t2(j, k, b, c) - t2(j, k, c, b))
                     end do
                  end do
               end do
            end do
         end do
      end do
      moo = 0.0_dp
      mvv = 0.0_dp
      do l = 1, no
         do j = 1, no
            do d = 1, nv
               moo(j, l) = moo(j, l) + 2.0_dp*pa(d, l, j, d) + pb(d, l, j, d)
            end do
         end do
      end do
      do b = 1, nv
         do d = 1, nv
            do l = 1, no
               mvv(d, b) = mvv(d, b) + 2.0_dp*pa(b, l, l, d) + pb(b, l, l, d)
            end do
         end do
      end do
   end subroutine pv_arrays

   subroutine goooo_of(no, nv, t1, t2, l2, g)
      !! goooo(i,j,k,l) = sum_ab tau(ijab) l2(klab) / 2
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv), l2(no, no, nv, nv)
      real(dp), intent(out) :: g(no, no, no, no)

      integer :: i, j, k, l, a, b

      g = 0.0_dp
      do b = 1, nv
         do a = 1, nv
            do l = 1, no
               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        g(i, j, k, l) = g(i, j, k, l) + 0.5_dp*(t2(i, j, a, b) + t1(i, a)*t1(j, b))*l2(k, l, a, b)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine goooo_of

   subroutine occupied_blocks(no, nv, t1, t2, l1, l2, pa, pb, moo, big)
      !! doooo and dooov, placed as _make_rdm2 does
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv), l1(no, nv), l2(no, no, nv, nv)
      real(dp), intent(in) :: pa(nv, no, no, nv), pb(nv, no, no, nv), moo(no, no)
      real(dp), intent(inout) :: big(no + nv, no + nv, no + nv, no + nv)

      real(dp) :: g(no, no, no, no), d4(no, no, no, no), gooov(no, no, no, nv), dooov(no, no, no, nv)
      integer :: i, j, k, l, a, b, c, p, q, r, s

      call goooo_of(no, nv, t1, t2, l2, g)
      do s = 1, no
         do r = 1, no
            do q = 1, no
               do p = 1, no
                  d4(p, q, r, s) = 2.0_dp*g(p, r, q, s) - g(p, r, s, q)
               end do
            end do
         end do
      end do
      do a = 1, nv
         do i = 1, no
            do k = 1, no
               do j = 1, no
                  gooov(j, k, i, a) = -0.5_dp*moo(j, i)*t1(k, a)
                  do c = 1, nv
                     gooov(j, k, i, a) = gooov(j, k, i, a) + t1(k, c)*pa(c, i, j, a) - t1(j, c)*pb(c, i, k, a)
                  end do
                  do l = 1, no
                     gooov(j, k, i, a) = gooov(j, k, i, a) + 2.0_dp*t1(l, a)*g(j, k, i, l)
                  end do
                  do b = 1, nv
                     gooov(j, k, i, a) = gooov(j, k, i, a) - l1(i, b)*(t2(j, k, b, a) + t1(j, b)*t1(k, a)) &
                                         - l2(j, k, b, a)*t1(i, b)
                  end do
               end do
            end do
         end do
      end do
      do a = 1, nv
         do r = 1, no
            do q = 1, no
               do p = 1, no
                  dooov(p, q, r, a) = 2.0_dp*gooov(p, r, q, a) - gooov(r, p, q, a)
               end do
            end do
         end do
      end do
      do s = 1, no
         do r = 1, no
            do q = 1, no
               do p = 1, no
                  big(p, q, r, s) = 2.0_dp*(d4(p, q, r, s) + d4(q, p, s, r))
               end do
            end do
         end do
      end do
      do a = 1, nv
         b = no + a
         do r = 1, no
            do q = 1, no
               do p = 1, no
                  big(p, q, r, b) = dooov(p, q, r, a)       ! ooov
                  big(r, b, p, q) = dooov(p, q, r, a)       ! ovoo
                  big(q, p, b, r) = dooov(p, q, r, a)       ! oovo
                  big(b, r, q, p) = dooov(p, q, r, a)       ! vooo
               end do
            end do
         end do
      end do
   end subroutine occupied_blocks

   subroutine ovov_block(no, nv, t1, t2, l1, l2, pa, pb, moo, mvv, big)
      !! dovov from pyscf's goovv, placed as _make_rdm2 does
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv), l1(no, nv), l2(no, no, nv, nv)
      real(dp), intent(in) :: pa(nv, no, no, nv), pb(nv, no, no, nv), moo(no, no), mvv(nv, nv)
      real(dp), intent(inout) :: big(no + nv, no + nv, no + nv, no + nv)

      real(dp) :: g(no, no, no, no), goovv(no, no, nv, nv), tau(no, no, nv, nv), mia(no, nv), mab(nv, nv)
      real(dp) :: mij(no, no), d(no, nv, no, nv), s
      integer :: i, j, k, l, a, b, c, dd

      call goooo_of(no, nv, t1, t2, l2, g)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  tau(i, j, a, b) = t2(i, j, a, b) + t1(i, a)*t1(j, b)
               end do
            end do
         end do
      end do
      mia = 0.0_dp
      do a = 1, nv
         do i = 1, no
            do c = 1, nv
               do k = 1, no
                  mia(i, a) = mia(i, a) + l1(k, c)*(2.0_dp*t2(i, k, a, c) - t2(i, k, c, a))
               end do
            end do
         end do
      end do
      do b = 1, nv
         do c = 1, nv
            mab(c, b) = sum(l1(:, c)*t1(:, b))
         end do
      end do
      do k = 1, no
         do j = 1, no
            mij(j, k) = sum(l1(k, :)*t1(j, :)) + 0.5_dp*moo(j, k)
         end do
      end do
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  s = mia(i, a)*t1(j, b) + 0.5_dp*tau(i, j, a, b) + 0.5_dp*l2(i, j, a, b)
                  do l = 1, no
                     do k = 1, no
                        s = s + g(i, j, k, l)*tau(k, l, a, b)
                     end do
                  end do
                  do k = 1, no
                     s = s - mij(j, k)*tau(i, k, a, b)
                  end do
                  do c = 1, nv
                     s = s - mab(c, b)*t2(i, j, a, c)
                  end do
                  do dd = 1, nv
                     s = s - 0.5_dp*mvv(b, dd)*tau(i, j, a, dd)
                  end do
                  do l = 1, no
                     do dd = 1, nv
                        s = s + pa(dd, l, i, b)*(tau(j, l, dd, a) - 0.5_dp*t2(j, l, dd, a)) &
                            - pb(dd, l, i, a)*(tau(j, l, dd, b) - 0.5_dp*t2(j, l, dd, b)) &
                            + (pb(dd, l, i, a) + 0.5_dp*pa(dd, l, i, a))*t2(j, l, b, dd)
                     end do
                  end do
                  goovv(i, j, a, b) = s
               end do
            end do
         end do
      end do
      do b = 1, nv
         do j = 1, no
            do a = 1, nv
               do i = 1, no
                  d(i, a, j, b) = 2.0_dp*goovv(i, j, a, b) - goovv(j, i, a, b)
               end do
            end do
         end do
      end do
      do b = 1, nv
         do j = 1, no
            do a = 1, nv
               do i = 1, no
                  big(i, no + a, j, no + b) = d(i, a, j, b) + d(j, b, i, a)
                  big(no + a, i, no + b, j) = d(i, a, j, b) + d(j, b, i, a)
               end do
            end do
         end do
      end do
   end subroutine ovov_block

   subroutine ovvo_blocks(no, nv, t1, l1, l2, pa, pb, big)
      !! dovvo and doovv, placed as _make_rdm2 does
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), l1(no, nv), l2(no, no, nv, nv)
      real(dp), intent(in) :: pa(nv, no, no, nv), pb(nv, no, no, nv)
      real(dp), intent(inout) :: big(no + nv, no + nv, no + nv, no + nv)

      real(dp) :: dovvo(no, nv, nv, no), doovv(no, no, nv, nv), govvo_big, govvo_small
      integer :: i, j, k, a, b, c

      do j = 1, no
         do b = 1, nv
            do a = 1, nv
               do i = 1, no
                  govvo_big = pa(a, i, j, b)                              ! gOvvO
                  govvo_small = l1(i, a)*t1(j, b) + pb(a, i, j, b)        ! govVO
                  do c = 1, nv
                     do k = 1, no
                        govvo_big = govvo_big + l2(k, i, a, c)*t1(j, c)*t1(k, b)
                        govvo_small = govvo_small - l2(i, k, a, c)*t1(j, c)*t1(k, b)
                     end do
                  end do
                  dovvo(i, a, b, j) = 2.0_dp*govvo_small + govvo_big
                  doovv(j, i, a, b) = -2.0_dp*govvo_big - govvo_small
               end do
            end do
         end do
      end do
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  big(i, j, no + a, no + b) = doovv(i, j, a, b) + doovv(j, i, b, a)
                  big(no + a, no + b, i, j) = doovv(i, j, a, b) + doovv(j, i, b, a)
                  big(i, no + a, no + b, j) = dovvo(i, a, b, j) + dovvo(j, b, a, i)
                  big(no + a, i, j, no + b) = dovvo(i, a, b, j) + dovvo(j, b, a, i)
               end do
            end do
         end do
      end do
   end subroutine ovvo_blocks

   subroutine virtual_blocks(no, nv, t1, t2, l1, l2, pa, pb, mvv, big)
      !! dvvvv and dovvv with the v^4 intermediate, placed as _make_rdm2 does
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv), l1(no, nv), l2(no, no, nv, nv)
      real(dp), intent(in) :: pa(nv, no, no, nv), pb(nv, no, no, nv), mvv(nv, nv)
      real(dp), intent(inout) :: big(no + nv, no + nv, no + nv, no + nv)

      real(dp), allocatable :: g(:, :, :, :), dv(:, :, :, :), gvovv(:, :, :, :), dovvv(:, :, :, :)
      integer :: i, j, k, a, b, c, d
      real(dp) :: s

      allocate (g(nv, nv, nv, nv), dv(nv, nv, nv, nv), gvovv(nv, no, nv, nv), dovvv(no, nv, nv, nv))
      g = 0.0_dp
      do d = 1, nv
         do c = 1, nv
            do b = 1, nv
               do a = 1, nv
                  do j = 1, no
                     do i = 1, no
                        g(a, b, c, d) = g(a, b, c, d) + l2(i, j, a, b)*(t2(i, j, c, d) + t1(i, c)*t1(j, d))
                     end do
                  end do
               end do
            end do
         end do
      end do
      do d = 1, nv
         do c = 1, nv
            do b = 1, nv
               do a = 1, nv
                  dv(a, c, b, d) = g(a, b, c, d) - 0.5_dp*g(a, b, d, c)
               end do
            end do
         end do
      end do
      do d = 1, nv
         do c = 1, nv
            do b = 1, nv
               do a = 1, nv
                  big(no + a, no + b, no + c, no + d) = 2.0_dp*(dv(a, b, c, d) + dv(b, a, d, c))
               end do
            end do
         end do
      end do
      do c = 1, nv
         do b = 1, nv
            do i = 1, no
               do a = 1, nv
                  s = 0.0_dp
                  do d = 1, nv
                     s = s - g(a, d, b, c)*t1(i, d)
                  end do
                  do k = 1, no
                     s = s + pb(a, k, i, c)*t1(k, b) - pa(a, k, i, b)*t1(k, c)
                  end do
                  do j = 1, no
                     s = s + l1(j, a)*t2(j, i, b, c) + l1(j, a)*t1(j, b)*t1(i, c) + t1(j, a)*l2(j, i, b, c)
                  end do
                  gvovv(a, i, b, c) = s + 0.5_dp*mvv(b, a)*t1(i, c)
               end do
            end do
         end do
      end do
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               do i = 1, no
                  dovvv(i, a, b, c) = 2.0_dp*gvovv(b, i, c, a) - gvovv(b, i, a, c)
               end do
            end do
         end do
      end do
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               do i = 1, no
                  big(i, no + a, no + b, no + c) = dovvv(i, a, b, c)       ! ovvv
                  big(no + b, no + c, i, no + a) = dovvv(i, a, b, c)       ! vvov
                  big(no + c, no + b, no + a, i) = dovvv(i, a, b, c)       ! vvvo
                  big(no + a, i, no + c, no + b) = dovvv(i, a, b, c)       ! vovv
               end do
            end do
         end do
      end do
   end subroutine virtual_blocks

   subroutine symmetrise(n, big)
      !! Average over p<->q and r<->s, then over (pq)<->(rs)
      integer, intent(in) :: n
      real(dp), intent(inout) :: big(n, n, n, n)

      real(dp), allocatable :: t(:, :, :, :)
      integer :: p, q, r, s

      allocate (t(n, n, n, n))
      do s = 1, n
         do r = 1, n
            do q = 1, n
               do p = 1, n
                  t(p, q, r, s) = 0.25_dp*(big(p, q, r, s) + big(q, p, r, s) + big(p, q, s, r) + big(q, p, s, r))
               end do
            end do
         end do
      end do
      do s = 1, n
         do r = 1, n
            do q = 1, n
               do p = 1, n
                  big(p, q, r, s) = 0.5_dp*(t(p, q, r, s) + t(r, s, p, q))
               end do
            end do
         end do
      end do
   end subroutine symmetrise

end module rccsd_density_host
