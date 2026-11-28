module rccsd_lambda_host
   !!
   !! The independent reference of check_rccsd_lambda_update: pyscf/cc/ccsd_lambda.py `make_intermediates` and
   !! `update_lambda` (Apache-2.0, The PySCF Developers), transcribed line by line in their DENSE form, with
   !! (ia|bc), (ab|cd) and the v^3 o intermediate wvvov, for canonical RHF orbitals (f_ov = 0, no level shift)
   !! and without the block loops. Nothing of the modules under test is used: the integrals are formed here from
   !! the three-index factors and every contraction is the plain loop of `es`.
   !!
   use pic_types, only: dp
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   implicit none
   private
   public :: host_lambda_update

   type :: ten_t
      !! A dense tensor of rank 2 or 4, first index fastest
      integer :: rank = 0
      integer :: d(4) = 1
      real(dp), allocatable :: x(:)
   end type ten_t

   interface operator(+)
      module procedure add
   end interface
   interface operator(-)
      module procedure subtract
   end interface
   interface operator(*)
      module procedure scaled
   end interface

   type :: ints_t
      !! The integral blocks, pyscf's names
      type(ten_t) :: oooo, ovoo, oovv, ovvo, ovvv, vvvv
   end type ints_t

   type :: imds_t
      !! pyscf's _IMDS
      type(ten_t) :: woooo, wvooo, wvoov_a, wvoov_b, wvvov, w1, w2, w3, w4
      !! wvoov_a is wVOov, wvoov_b is wvOOv
   end type imds_t

contains

   subroutine host_lambda_update(x, t1_in, t2_in, l1_in, l2_in, l1_out, l2_out)
      !! (l1_out, l2_out) = pyscf's update_lambda of (l1_in, l2_in) at the amplitudes (t1_in, t2_in)
      type(trc_rccsd_ints_t), intent(in) :: x
      real(dp), intent(in) :: t1_in(:, :), t2_in(:, :, :, :), l1_in(:, :), l2_in(:, :, :, :)
      real(dp), intent(out) :: l1_out(:, :), l2_out(:, :, :, :)

      type(ints_t) :: g
      type(imds_t) :: w
      type(ten_t) :: t1, t2, l1, l2, new1, new2
      integer :: no, nv, i, j, a, b
      real(dp) :: dia, djb

      no = x%no
      nv = x%nv
      call integrals(x, g)
      t1 = tensor([no, nv], reshape(t1_in, [no*nv]))
      t2 = tensor([no, no, nv, nv], reshape(t2_in, [no*no*nv*nv]))
      l1 = tensor([no, nv], reshape(l1_in, [no*nv]))
      l2 = tensor([no, no, nv, nv], reshape(l2_in, [no*no*nv*nv]))
      call make_intermediates(x, g, t1, t2, w)
      call update_lambda(x, g, w, [t1, t2, l1, l2], new1, new2)
      do a = 1, nv
         do i = 1, no
            l1_out(i, a) = new1%x(i + no*(a - 1))/(x%eps_o(i) - x%eps_v(a))
         end do
      end do
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  dia = x%eps_o(i) - x%eps_v(a)
                  djb = x%eps_o(j) - x%eps_v(b)
                  l2_out(i, j, a, b) = (new2%x(1 + (i - 1) + no*((j - 1) + no*((a - 1) + nv*(b - 1)))) &
                                        + new2%x(1 + (j - 1) + no*((i - 1) + no*((b - 1) + nv*(a - 1)))))/(dia + djb)
               end do
            end do
         end do
      end do
   end subroutine host_lambda_update

   subroutine make_intermediates(x, g, t1, t2, w)
      !! pyscf's make_intermediates, canonical orbitals (fov = 0, w4 starts at 0)
      type(trc_rccsd_ints_t), intent(in) :: x
      type(ints_t), intent(in) :: g
      type(ten_t), intent(in) :: t1, t2
      type(imds_t), intent(out) :: w

      type(ten_t) :: vvov, theta, tau, g2vovv, g2vvov, vackb, g2ovoo, vcjik, voov, g2voov, tmpw4, tmp, vo_a, vo_b
      type(ten_t) :: ov1, ov2

      w%w1 = diagonal(x%eps_v)
      w%w2 = diagonal(x%eps_o)
      vvov = tr("iacd->cdia", g%ovvv)
      w%w1 = w%w1 + 2.0_dp*es("jcba,jc->ba", g%ovvv, t1)
      w%w1 = w%w1 - es("jabc,jc->ba", g%ovvv, t1)
      theta = 2.0_dp*t2 - tr("ijab->jiab", t2)
      w%w3 = es("jkcd,kdcb->bj", 2.0_dp*t2 - tr("ijab->ijba", t2), g%ovvv)
      w%wvoov_a = es("jbcd,kd->bjkc", g%ovvv, t1)
      w%wvoov_b = es("cbjd,kd->cjkb", vvov, (-1.0_dp)*t1)
      g2vovv = 2.0_dp*tr("abcd->acbd", vvov) - tr("abcd->acdb", vvov)
      tau = t2 + es("ia,jb->ijab", t1, t1)
      w%wvooo = es("cibd,jkbd->ckij", g2vovv, tau)

      w%wvvov = (-1.5_dp)*es("jabd,jkcd->abkc", g%ovvv, t2)
      w%wvvov = w%wvvov + 2.0_dp*tr("abcd->adcb", vvov) - vvov
      g2vvov = 2.0_dp*vvov - tr("abcd->bcad", g%ovvv)
      vackb = es("acjd,kjbd->ackb", g2vvov, 2.0_dp*t2 - tr("ijab->ijba", t2))
      w%wvvov = w%wvvov + tr("abcd->adcb", vackb) - 0.5_dp*vackb

      w%w2 = w%w2 + 2.0_dp*es("kbij,kb->ij", g%ovoo, t1)
      w%w2 = w%w2 - es("ibkj,kb->ij", g%ovoo, t1)
      theta = 2.0_dp*tr("ijab->jiab", t2) - t2
      w%w3 = w%w3 - es("lckj,klcb->bj", g%ovoo, theta)

      tmp = es("lc,jcik->ijkl", t1, g%ovoo)
      w%woooo = tmp + tr("ijkl->jilk", tmp)
      w%wvoov_b = w%wvoov_b + es("lbjk,lc->bjkc", g%ovoo, t1)
      w%wvoov_a = w%wvoov_a - es("jbkl,lc->bjkc", g%ovoo, t1)
      w%wvooo = w%wvooo + 2.0_dp*tr("abcd->bdca", g%ovoo)
      w%wvooo = w%wvooo - tr("abcd->bacd", g%ovoo)
      w%wvooo = w%wvooo - es("klbc,iblj->ckij", t2, 1.5_dp*g%ovoo)

      g2ovoo = 2.0_dp*g%ovoo - tr("abcd->cbad", g%ovoo)
      theta = 2.0_dp*t2 - tr("ijab->jiab", t2)
      vcjik = es("jlcb,lbki->cjki", theta, g2ovoo)
      w%wvooo = w%wvooo + tr("abcd->adcb", vcjik) - 0.5_dp*vcjik

      voov = tr("abcd->badc", g%ovvo)
      tau = t2 + es("ia,jb->ijab", t1, t1)
      w%woooo = w%woooo + es("cijd,klcd->ijkl", voov, tau)

      g2voov = 2.0_dp*voov - tr("abcd->acbd", voov)
      tmpw4 = es("ckld,ld->kc", g2voov, t1)
      w%w1 = w%w1 - es("ckja,kjcb->ba", g2voov, t2)
      w%w1 = w%w1 - es("ja,jb->ba", tmpw4, t1)
      w%w2 = w%w2 + es("jkbc,bikc->ij", t2, g2voov)
      w%w2 = w%w2 + es("ib,jb->ij", tmpw4, t1)
      w%w3 = w%w3 + es("bc,jc->bj", es("kb,kc->bc", t1, tmpw4), t1)
      w%w4 = tmpw4

      w%wvoov_b = w%wvoov_b + es("bljk,lc->bjkc", es("bljd,kd->bljk", voov, t1), t1)
      w%wvoov_a = w%wvoov_a - es("bjlk,lc->bjkc", es("bjld,kd->bjlk", voov, t1), t1)

      vo_a = es("bjld,klcd->bjkc", g2voov, t2) - es("bjld,kldc->bjkc", voov, t2) + voov
      vo_b = es("bljd,kldc->bjkc", voov, t2) - tr("abcd->cbad", g%oovv)
      w%wvoov_a = w%wvoov_a + vo_a
      w%wvoov_b = w%wvoov_b + vo_b

      ov1 = 2.0_dp*vo_b + vo_a
      ov2 = 2.0_dp*vo_a + vo_b
      w%wvooo = w%wvooo - es("jb,bikc->ckij", t1, ov1)
      w%wvooo = w%wvooo + es("kb,bijc->ckij", t1, ov2)
      w%w3 = w%w3 + es("ckjb,kc->bj", ov2, t1)

      w%wvvov = w%wvvov + es("ajkc,jb->abkc", ov1, t1)
      w%wvvov = w%wvvov - es("ajkb,jc->abkc", ov2, t1)

      g2ovoo = 2.0_dp*g%ovoo - tr("abcd->cbad", g%ovoo)
      tau = t2 + es("ia,jb->ijab", t1, t1)
      w%wvvov = w%wvvov + es("laki,klbc->abic", g2ovoo, tau)

      w%woooo = w%woooo + tr("abcd->acbd", g%oooo)
      w%w3 = w%w3 + es("bc,jc->bj", w%w1, t1)
      w%w3 = w%w3 - es("kj,kb->bj", w%w2, t1)
   end subroutine make_intermediates

   subroutine update_lambda(x, g, w, amps, l1new, l2new)
      !! pyscf's update_lambda up to the denominators; amps = [t1, t2, l1, l2]
      type(trc_rccsd_ints_t), intent(in) :: x
      type(ints_t), intent(in) :: g
      type(imds_t), intent(in) :: w
      type(ten_t), intent(in) :: amps(4)
      type(ten_t), intent(out) :: l1new, l2new

      type(ten_t) :: t1, t2, l1, l2, theta, mba, mij, mba1, mij1, mia1, w1, w2, tau, l2tau, l2t1, m4, voov, tmp, m3

      t1 = amps(1)
      t2 = amps(2)
      l1 = amps(3)
      l2 = amps(4)
      theta = 2.0_dp*t2 - tr("ijab->ijba", t2)
      mba = es("klca,klcb->ba", l2, theta)
      mij = es("ikcd,jkcd->ij", l2, theta)
      mba1 = es("jc,jb->bc", l1, t1) + mba
      mij1 = es("kb,jb->kj", l1, t1) + mij
      mia1 = t1 + 2.0_dp*es("kc,jkbc->jb", l1, t2)
      mia1 = mia1 - es("kc,jkcb->jb", l1, t2)
      mia1 = mia1 - es("jk,kb->jb", es("jc,kc->jk", t1, l1), t1)
      mia1 = mia1 - es("bd,jd->jb", mba, t1)
      mia1 = mia1 - es("lj,lb->jb", mij, t1)

      ! the particle ladder of l2
      l2new = es("ijab,acbd->ijcd", l2, g%vvvv)
      l1new = 2.0_dp*es("ijab,jb->ia", l2new, t1)
      l1new = l1new - es("jiab,jb->ia", l2new, t1)
      l2new = 0.5_dp*l2new

      w1 = w%w1 - diagonal(x%eps_v)
      w2 = w%w2 - diagonal(x%eps_o)
      l1new = l1new + es("ib,ba->ia", l1, w1)
      l1new = l1new - es("ja,ij->ia", l1, w2)
      l1new = l1new - es("ik,ka->ia", mij, w%w4)
      l1new = l1new - es("ca,ic->ia", mba, w%w4)
      l1new = l1new + 2.0_dp*es("ijab,bj->ia", l2, w%w3)
      l1new = l1new - es("ijba,bj->ia", l2, w%w3)

      l2new = l2new + es("ia,jb->ijab", l1, w%w4)
      l2new = l2new + es("jibc,ca->jiba", l2, w1)
      l2new = l2new - es("jk,kiba->jiba", w2, l2)

      l1new = l1new - 2.0_dp*es("iajk,kj->ia", g%ovoo, mij1)
      l1new = l1new + es("jaik,kj->ia", g%ovoo, mij1)
      l2new = l2new - es("jbki,ka->jiba", g%ovoo, l1)

      tau = t2 + es("ia,jb->ijab", t1, t1)
      l2tau = es("ijcd,klcd->ijkl", l2, tau)
      l2t1 = es("jidc,kc->ijkd", l2, t1)

      l1new = l1new - es("jb,jiab->ia", l1, g%oovv)
      l1new = l1new + 2.0_dp*es("iabc,bc->ia", g%ovvv, mba1)
      l1new = l1new - es("ibca,bc->ia", g%ovvv, mba1)
      l2new = l2new + es("jbac,ic->jiba", g%ovvv, l1)
      m4 = es("ijkd,kadb->ijab", l2t1, g%ovvv)
      l2new = l2new - m4
      l1new = l1new - 2.0_dp*es("ijab,jb->ia", m4, t1)
      l1new = l1new - 2.0_dp*es("ijab,ia->jb", m4, t1)
      l1new = l1new + es("jiab,jb->ia", m4, t1)
      l1new = l1new + es("jiab,ia->jb", m4, t1)

      voov = tr("abcd->badc", g%ovvo)
      l1new = l1new + 2.0_dp*es("jb,aijb->ia", l1, voov)
      l2new = l2new + 0.5_dp*tr("abcd->bcad", voov)
      l2new = l2new - es("bjic,ca->jiba", voov, mba1)
      l2new = l2new - es("bjka,ik->jiba", voov, mij1)
      l1new = l1new + 2.0_dp*es("aijb,jb->ia", voov, mia1)
      l1new = l1new - es("bija,jb->ia", voov, mia1)
      m4 = es("ijkl,aklb->ijab", l2tau, voov)
      l2new = l2new + 0.5_dp*m4
      l1new = l1new + 2.0_dp*es("ijab,jb->ia", m4, t1)
      l1new = l1new - es("ijba,jb->ia", m4, t1)

      l1new = l1new - es("ckij,jkca->ia", w%wvooo, l2)
      l1new = l1new + es("abkc,kibc->ia", w%wvvov, l2)

      tmp = tr("abcd->acbd", l2) - 0.5_dp*tr("abcd->adbc", l2)
      l2new = l2new + es("iakc,bjkc->jiba", tmp, 2.0_dp*w%wvoov_a + w%wvoov_b)
      tmp = es("jkca,bikc->jiba", l2, w%wvoov_b)
      l2new = l2new + tmp + 0.5_dp*tr("abcd->bacd", tmp)

      m3 = es("ijkl,klab->ijab", w%woooo, l2)
      l2new = l2new + 0.5_dp*m3
      l1new = l1new + 2.0_dp*es("ijab,jb->ia", m3, t1)
      l1new = l1new - es("ijba,jb->ia", m3, t1)
   end subroutine update_lambda

   subroutine integrals(x, g)
      !! The blocks of (pq|rs) = sum_P b(pq,P) b(rs,P), occupied first
      type(trc_rccsd_ints_t), intent(in) :: x
      type(ints_t), intent(out) :: g

      type(ten_t) :: boo, bov, bvv
      real(dp), allocatable :: oo(:, :), ov(:, :), vv(:, :)
      integer :: no, nv

      no = x%no
      nv = x%nv
      oo = x%b_oo
      ov = x%b_ov
      vv = x%b_vv
      g%oooo = block4([no, no, no, no], oo, oo)
      g%ovoo = block4([no, nv, no, no], ov, oo)
      g%oovv = block4([no, no, nv, nv], oo, vv)
      g%ovvv = block4([no, nv, nv, nv], ov, vv)
      g%vvvv = block4([nv, nv, nv, nv], vv, vv)
      ! (ia|bj) = (ia|jb)
      g%ovvo = tr("iajb->iabj", block4([no, nv, no, nv], ov, ov))
   end subroutine integrals

   function block4(extent, left, right) result(c)
      !! sum_P left(pq,P) right(rs,P) as a tensor (p,q,r,s)
      integer, intent(in) :: extent(4)
      real(dp), intent(in) :: left(:, :), right(:, :)
      type(ten_t) :: c

      real(dp), allocatable :: product_(:, :)

      product_ = matmul(left, transpose(right))
      c = tensor(extent, reshape(product_, [size(product_)]))
   end function block4

   function tensor(extent, values) result(c)
      integer, intent(in) :: extent(:)
      real(dp), intent(in) :: values(:)
      type(ten_t) :: c

      c%rank = size(extent)
      c%d = 1
      c%d(1:c%rank) = extent
      c%x = values
   end function tensor

   function diagonal(e) result(c)
      real(dp), intent(in) :: e(:)
      type(ten_t) :: c

      integer :: i, n

      n = size(e)
      c%rank = 2
      c%d = [n, n, 1, 1]
      allocate (c%x(n*n))
      c%x = 0.0_dp
      do i = 1, n
         c%x(i + n*(i - 1)) = e(i)
      end do
   end function diagonal

   function add(a, b) result(c)
      type(ten_t), intent(in) :: a, b
      type(ten_t) :: c

      if (any(a%d /= b%d)) error stop "rccsd_lambda_host: shapes differ in a sum"
      c = a
      c%x = a%x + b%x
   end function add

   function subtract(a, b) result(c)
      type(ten_t), intent(in) :: a, b
      type(ten_t) :: c

      if (any(a%d /= b%d)) error stop "rccsd_lambda_host: shapes differ in a difference"
      c = a
      c%x = a%x - b%x
   end function subtract

   function scaled(s, a) result(c)
      real(dp), intent(in) :: s
      type(ten_t), intent(in) :: a
      type(ten_t) :: c

      c = a
      c%x = s*a%x
   end function scaled

   function tr(spec, a) result(c)
      !! numpy's transpose spelled as an einsum of one operand, "abcd->acbd"
      character(len=*), intent(in) :: spec
      type(ten_t), intent(in) :: a
      type(ten_t) :: c

      type(ten_t) :: one
      integer :: arrow

      arrow = index(spec, "->")
      one = tensor([1, 1], [1.0_dp])
      c = es(spec(1:arrow - 1)//",yz->"//spec(arrow + 2:), a, one)
   end function tr

   function es(spec, a, b) result(c)
      !! einsum of two operands: every index not in the output is summed. One loop over the whole index space,
      !! the three offsets advanced like an odometer.
      character(len=*), intent(in) :: spec
      type(ten_t), intent(in) :: a, b
      type(ten_t) :: c

      character(len=12) :: letters
      integer :: extent(12), sa(12), sb(12), sc(12), at(12)
      integer :: comma, arrow, n, i, k, oa, ob, oc, stride, total, step

      comma = index(spec, ",")
      arrow = index(spec, "->")
      letters = ""
      n = 0
      sa = 0
      sb = 0
      sc = 0
      extent = 1
      stride = 1
      do i = 1, comma - 1
         k = letter(spec(i:i), a%d(i))
         sa(k) = sa(k) + stride
         stride = stride*a%d(i)
      end do
      stride = 1
      do i = 1, arrow - comma - 1
         k = letter(spec(comma + i:comma + i), b%d(i))
         sb(k) = sb(k) + stride
         stride = stride*b%d(i)
      end do
      c%rank = len_trim(spec) - arrow - 1
      c%d = 1
      stride = 1
      do i = 1, c%rank
         k = index(letters(1:n), spec(arrow + 1 + i:arrow + 1 + i))
         if (k == 0) error stop "rccsd_lambda_host: an output index of no operand"
         c%d(i) = extent(k)
         sc(k) = stride
         stride = stride*extent(k)
      end do
      allocate (c%x(stride))
      c%x = 0.0_dp
      total = product(extent(1:n))
      at = 0
      oa = 1
      ob = 1
      oc = 1
      do step = 1, total
         c%x(oc) = c%x(oc) + a%x(oa)*b%x(ob)
         do k = 1, n
            at(k) = at(k) + 1
            oa = oa + sa(k)
            ob = ob + sb(k)
            oc = oc + sc(k)
            if (at(k) < extent(k)) exit
            at(k) = 0
            oa = oa - sa(k)*extent(k)
            ob = ob - sb(k)*extent(k)
            oc = oc - sc(k)*extent(k)
         end do
      end do

   contains

      function letter(name, length) result(k)
         character, intent(in) :: name
         integer, intent(in) :: length
         integer :: k

         k = index(letters(1:n), name)
         if (k == 0) then
            n = n + 1
            letters(n:n) = name
            extent(n) = length
            k = n
         else if (extent(k) /= length) then
            error stop "rccsd_lambda_host: one index with two extents in "//spec
         end if
      end function letter

   end function es

end module rccsd_lambda_host
