!
! The translated exchange-correlation kernels against libxc.
!
! `reference_xc.dat` holds, for every kernel in both spin cases, a set of
! densities and gradient invariants and what libxc 7 returned for them
! through second order: eps, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2.
! This program evaluates terco's kernels -- ExchCXX's Maple arithmetic,
! translated by tools/exchcxx_to_fortran.py -- at the same points and
! reports the worst relative disagreement per output.
!
! Two independent implementations of the same functional, from two
! different codes' Maple sessions, agreeing to 1e-10 across eleven decades
! of density is the evidence that the translation is right. It is a
! stronger check than any identity, because a wrong phase or a dropped
! factor of two inside a `merge` passes every identity there is.
!
! The four evaluation routines of each kernel are also compared against
! each other where they overlap: `exc` against the eps `exc_vxc` gives, and
! `fxc` against the second derivatives `vxc_fxc` gives. Those are different
! Maple expressions for the same quantity and must agree to rounding.
!
program check_xc
   use trc_boys, only: dp
   use trc_xc_slater_exchange, only: slater_exchange_exc_unpolar, slater_exchange_exc_vxc_unpolar, &
                                     slater_exchange_fxc_unpolar, slater_exchange_vxc_fxc_unpolar, &
                                     slater_exchange_exc_polar, slater_exchange_exc_vxc_polar, &
                                     slater_exchange_fxc_polar, slater_exchange_vxc_fxc_polar
   use trc_xc_vwn, only: vwn_exc_unpolar, vwn_exc_vxc_unpolar, vwn_fxc_unpolar, vwn_vxc_fxc_unpolar, &
                         vwn_exc_polar, vwn_exc_vxc_polar, vwn_fxc_polar, vwn_vxc_fxc_polar
   use trc_xc_vwn_rpa, only: vwn_rpa_exc_unpolar, vwn_rpa_exc_vxc_unpolar, vwn_rpa_fxc_unpolar, &
                             vwn_rpa_vxc_fxc_unpolar, vwn_rpa_exc_polar, vwn_rpa_exc_vxc_polar, &
                             vwn_rpa_fxc_polar, vwn_rpa_vxc_fxc_polar
   use trc_xc_pbe_x, only: pbe_x_exc_unpolar, pbe_x_exc_vxc_unpolar, pbe_x_fxc_unpolar, &
                           pbe_x_vxc_fxc_unpolar, pbe_x_exc_polar, pbe_x_exc_vxc_polar, &
                           pbe_x_fxc_polar, pbe_x_vxc_fxc_polar
   use trc_xc_pbe_c, only: pbe_c_exc_unpolar, pbe_c_exc_vxc_unpolar, pbe_c_fxc_unpolar, &
                           pbe_c_vxc_fxc_unpolar, pbe_c_exc_polar, pbe_c_exc_vxc_polar, &
                           pbe_c_fxc_polar, pbe_c_vxc_fxc_polar
   use trc_xc_b88, only: b88_exc_unpolar, b88_exc_vxc_unpolar, b88_fxc_unpolar, b88_vxc_fxc_unpolar, &
                         b88_exc_polar, b88_exc_vxc_polar, b88_fxc_polar, b88_vxc_fxc_polar
   use trc_xc_lyp, only: lyp_exc_unpolar, lyp_exc_vxc_unpolar, lyp_fxc_unpolar, lyp_vxc_fxc_unpolar, &
                         lyp_exc_polar, lyp_exc_vxc_polar, lyp_fxc_polar, lyp_vxc_fxc_polar
   implicit none

   ! Two independent Maple derivations agree to this on energies and first
   ! derivatives; measured 1.2e-12 at worst, so a factor of ten of slack.
   real(dp), parameter :: tol_v = 1.0e-11_dp
   ! Second derivatives, measured 3.7e-12 at worst.
   !
   ! Every output is measured against the largest value in its derivative
   ! group at that point, not against itself. The spin cross terms of
   ! exchange are identically zero, and we return exactly zero; libxc
   ! returns the residue of cancelling terms of order 1e26, and a relative
   ! error against that residue is meaningless.
   !
   ! The reference stops at a spin density of 5e-12 because that is where
   ! libxc's polarised GGA starts to drift: checked against the closed form
   ! at sigma = 0, where PBE exchange is Slater exchange, ours matches it to
   ! fourteen digits down to 5e-13 and libxc to nine at 5e-13, twelve at
   ! 5e-12, thirteen at 1e-11. A looser tolerance would have hidden that;
   ! a tighter sample records it.
   real(dp), parameter :: tol_f = 1.0e-10_dp
   ! The kernel's own routines against each other: same expressions
   ! reassembled, so only rounding.
   real(dp), parameter :: tol_self = 1.0e-12_dp
   real(dp), parameter :: floor = 1.0e-14_dp

   integer, parameter :: maxout = 21, maxin = 5
   character(len=32) :: kernel, family
   character(len=256) :: line
   integer :: polar, npts, nin, nout, ios, ip, io, unit, nblocks, nfail
   real(dp) :: xin(maxin), ref(maxout), got(maxout), self_err, worst(maxout), worst_self
   real(dp) :: err, worst_got(maxout), worst_ref(maxout), worst_x(maxin, maxout)
   integer :: worst_at(maxout), grp(maxout), order(maxout), ig
   real(dp) :: scale(maxout), tol(maxout), worst_order(0:2)
   logical :: ok

   open (newunit=unit, file="reference_xc.dat", status="old", action="read", iostat=ios)
   if (ios /= 0) then
      print '(a)', "check_xc: cannot open reference_xc.dat (run from test/)"
      stop 1
   end if

   nblocks = 0
   nfail = 0
   do
      read (unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == "#") cycle
      read (line, *) kernel, family, polar, npts, nin, nout
      nblocks = nblocks + 1
      call column_groups(trim(family), polar == 1, nout, grp, order)
      tol(1:nout) = merge(tol_f, tol_v, order(1:nout) == 2)
      worst = 0.0_dp
      worst_self = 0.0_dp
      worst_at = 0
      worst_got = 0.0_dp
      worst_ref = 0.0_dp
      worst_x = 0.0_dp
      do ip = 1, npts
         read (unit, *) xin(1:nin), ref(1:nout)
         call evaluate(trim(kernel), polar == 1, xin, got, self_err)
         worst_self = max(worst_self, self_err)
         do ig = 1, maxval(grp(1:nout))
            scale(ig) = max(floor, maxval(abs(ref(1:nout)), mask=grp(1:nout) == ig))
         end do
         do io = 1, nout
            err = abs(got(io) - ref(io))/scale(grp(io))
            if (err > worst(io)) then
               worst(io) = err
               worst_at(io) = ip
               worst_got(io) = got(io)
               worst_ref(io) = ref(io)
               worst_x(:, io) = xin
            end if
         end do
      end do
      ok = all(worst(1:nout) <= tol(1:nout)) .and. worst_self <= tol_self
      if (.not. ok) nfail = nfail + 1
      do ig = 0, 2
         worst_order(ig) = maxval(worst(1:nout), mask=order(1:nout) == ig)
      end do
      print '(a,1x,a16,1x,a5,1x,l1,1x,a,3es9.2,a,es9.2)', &
         merge("ok  ", "FAIL", ok), kernel, family, polar == 1, &
         "vs libxc: eps, vxc, fxc ", worst_order, "  self ", worst_self
      if (.not. ok) then
         do io = 1, nout
            if (worst(io) > tol(io)) then
               print '(a,i3,a,es10.3,a,i4,a,2es24.16)', "     output ", io, " off by ", worst(io), &
                  " at point ", worst_at(io), ": ours, libxc =", worst_got(io), worst_ref(io)
               print '(a,5es14.6)', "        inputs ", worst_x(1:nin, io)
            end if
         end do
      end if
   end do
   close (unit)

   print '(a,i0,a,i0,a)', "check_xc: ", nblocks - nfail, " of ", nblocks, " blocks agree with libxc"
   if (nfail > 0 .or. nblocks == 0) stop 1

contains

   !
   ! Which reference columns belong together -- eps; vrho; vsigma; v2rho2;
   ! v2rhosigma; v2sigma2 -- and the derivative order of each.
   !
   subroutine column_groups(family, polar, nout, grp, order)
      character(len=*), intent(in) :: family
      logical, intent(in) :: polar
      integer, intent(in) :: nout
      integer, intent(out) :: grp(:), order(:)
      integer :: sizes(6), orders(6), ng, g, k, io

      if (family == "lda") then
         ng = 3
         sizes(1:3) = merge([1, 2, 3], [1, 1, 1], polar)
         orders(1:3) = [0, 1, 2]
      else
         ng = 6
         sizes(1:6) = merge([1, 2, 3, 3, 6, 6], [1, 1, 1, 1, 1, 1], polar)
         orders(1:6) = [0, 1, 1, 2, 2, 2]
      end if
      io = 0
      do g = 1, ng
         do k = 1, sizes(g)
            io = io + 1
            grp(io) = g
            order(io) = orders(g)
         end do
      end do
      if (io /= nout) then
         print '(a)', "check_xc: column count does not match the family"
         stop 1
      end if
   end subroutine column_groups

   !
   ! Fill `out` in the reference file's column order -- eps, vrho..., vsigma...,
   ! v2rho2..., v2rhosigma..., v2sigma2... -- from `exc_vxc` and `fxc`, and
   ! return in `self` the worst disagreement of `exc` and `vxc_fxc` with them.
   !
   subroutine evaluate(name, polar, x, out, self)
      character(len=*), intent(in) :: name
      logical, intent(in) :: polar
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: out(:), self
      real(dp) :: e, v(5), f(15), v2(5), f2(15)
      integer :: nv, nf

      out = 0.0_dp
      v = 0.0_dp; f = 0.0_dp; v2 = 0.0_dp; f2 = 0.0_dp
      select case (name)
      case ("slater_exchange")
         if (polar) then
            call slater_exchange_exc_vxc_polar(x(1), x(2), out(1), v(1), v(2))
            call slater_exchange_fxc_polar(x(1), x(2), f(1), f(2), f(3))
            call slater_exchange_exc_polar(x(1), x(2), e)
            call slater_exchange_vxc_fxc_polar(x(1), x(2), v2(1), v2(2), f2(1), f2(2), f2(3))
            nv = 2; nf = 3
         else
            call slater_exchange_exc_vxc_unpolar(x(1), out(1), v(1))
            call slater_exchange_fxc_unpolar(x(1), f(1))
            call slater_exchange_exc_unpolar(x(1), e)
            call slater_exchange_vxc_fxc_unpolar(x(1), v2(1), f2(1))
            nv = 1; nf = 1
         end if
      case ("vwn")
         if (polar) then
            call vwn_exc_vxc_polar(x(1), x(2), out(1), v(1), v(2))
            call vwn_fxc_polar(x(1), x(2), f(1), f(2), f(3))
            call vwn_exc_polar(x(1), x(2), e)
            call vwn_vxc_fxc_polar(x(1), x(2), v2(1), v2(2), f2(1), f2(2), f2(3))
            nv = 2; nf = 3
         else
            call vwn_exc_vxc_unpolar(x(1), out(1), v(1))
            call vwn_fxc_unpolar(x(1), f(1))
            call vwn_exc_unpolar(x(1), e)
            call vwn_vxc_fxc_unpolar(x(1), v2(1), f2(1))
            nv = 1; nf = 1
         end if
      case ("vwn_rpa")
         if (polar) then
            call vwn_rpa_exc_vxc_polar(x(1), x(2), out(1), v(1), v(2))
            call vwn_rpa_fxc_polar(x(1), x(2), f(1), f(2), f(3))
            call vwn_rpa_exc_polar(x(1), x(2), e)
            call vwn_rpa_vxc_fxc_polar(x(1), x(2), v2(1), v2(2), f2(1), f2(2), f2(3))
            nv = 2; nf = 3
         else
            call vwn_rpa_exc_vxc_unpolar(x(1), out(1), v(1))
            call vwn_rpa_fxc_unpolar(x(1), f(1))
            call vwn_rpa_exc_unpolar(x(1), e)
            call vwn_rpa_vxc_fxc_unpolar(x(1), v2(1), f2(1))
            nv = 1; nf = 1
         end if
      case ("pbe_x")
         if (polar) then
            call pbe_x_exc_vxc_polar(x(1), x(2), x(3), x(4), x(5), out(1), v(1), v(2), v(3), v(4), v(5))
            call pbe_x_fxc_polar(x(1), x(2), x(3), x(4), x(5), f(1), f(2), f(3), f(4), f(5), f(6), &
                                 f(7), f(8), f(9), f(10), f(11), f(12), f(13), f(14), f(15))
            call pbe_x_exc_polar(x(1), x(2), x(3), x(4), x(5), e)
            call pbe_x_vxc_fxc_polar(x(1), x(2), x(3), x(4), x(5), v2(1), v2(2), v2(3), v2(4), v2(5), &
                                     f2(1), f2(2), f2(3), f2(4), f2(5), f2(6), f2(7), f2(8), f2(9), &
                                     f2(10), f2(11), f2(12), f2(13), f2(14), f2(15))
            nv = 5; nf = 15
         else
            call pbe_x_exc_vxc_unpolar(x(1), x(2), out(1), v(1), v(2))
            call pbe_x_fxc_unpolar(x(1), x(2), f(1), f(2), f(3))
            call pbe_x_exc_unpolar(x(1), x(2), e)
            call pbe_x_vxc_fxc_unpolar(x(1), x(2), v2(1), v2(2), f2(1), f2(2), f2(3))
            nv = 2; nf = 3
         end if
      case ("pbe_c")
         if (polar) then
            call pbe_c_exc_vxc_polar(x(1), x(2), x(3), x(4), x(5), out(1), v(1), v(2), v(3), v(4), v(5))
            call pbe_c_fxc_polar(x(1), x(2), x(3), x(4), x(5), f(1), f(2), f(3), f(4), f(5), f(6), &
                                 f(7), f(8), f(9), f(10), f(11), f(12), f(13), f(14), f(15))
            call pbe_c_exc_polar(x(1), x(2), x(3), x(4), x(5), e)
            call pbe_c_vxc_fxc_polar(x(1), x(2), x(3), x(4), x(5), v2(1), v2(2), v2(3), v2(4), v2(5), &
                                     f2(1), f2(2), f2(3), f2(4), f2(5), f2(6), f2(7), f2(8), f2(9), &
                                     f2(10), f2(11), f2(12), f2(13), f2(14), f2(15))
            nv = 5; nf = 15
         else
            call pbe_c_exc_vxc_unpolar(x(1), x(2), out(1), v(1), v(2))
            call pbe_c_fxc_unpolar(x(1), x(2), f(1), f(2), f(3))
            call pbe_c_exc_unpolar(x(1), x(2), e)
            call pbe_c_vxc_fxc_unpolar(x(1), x(2), v2(1), v2(2), f2(1), f2(2), f2(3))
            nv = 2; nf = 3
         end if
      case ("b88")
         if (polar) then
            call b88_exc_vxc_polar(x(1), x(2), x(3), x(4), x(5), out(1), v(1), v(2), v(3), v(4), v(5))
            call b88_fxc_polar(x(1), x(2), x(3), x(4), x(5), f(1), f(2), f(3), f(4), f(5), f(6), &
                               f(7), f(8), f(9), f(10), f(11), f(12), f(13), f(14), f(15))
            call b88_exc_polar(x(1), x(2), x(3), x(4), x(5), e)
            call b88_vxc_fxc_polar(x(1), x(2), x(3), x(4), x(5), v2(1), v2(2), v2(3), v2(4), v2(5), &
                                   f2(1), f2(2), f2(3), f2(4), f2(5), f2(6), f2(7), f2(8), f2(9), &
                                   f2(10), f2(11), f2(12), f2(13), f2(14), f2(15))
            nv = 5; nf = 15
         else
            call b88_exc_vxc_unpolar(x(1), x(2), out(1), v(1), v(2))
            call b88_fxc_unpolar(x(1), x(2), f(1), f(2), f(3))
            call b88_exc_unpolar(x(1), x(2), e)
            call b88_vxc_fxc_unpolar(x(1), x(2), v2(1), v2(2), f2(1), f2(2), f2(3))
            nv = 2; nf = 3
         end if
      case ("lyp")
         if (polar) then
            call lyp_exc_vxc_polar(x(1), x(2), x(3), x(4), x(5), out(1), v(1), v(2), v(3), v(4), v(5))
            call lyp_fxc_polar(x(1), x(2), x(3), x(4), x(5), f(1), f(2), f(3), f(4), f(5), f(6), &
                               f(7), f(8), f(9), f(10), f(11), f(12), f(13), f(14), f(15))
            call lyp_exc_polar(x(1), x(2), x(3), x(4), x(5), e)
            call lyp_vxc_fxc_polar(x(1), x(2), x(3), x(4), x(5), v2(1), v2(2), v2(3), v2(4), v2(5), &
                                   f2(1), f2(2), f2(3), f2(4), f2(5), f2(6), f2(7), f2(8), f2(9), &
                                   f2(10), f2(11), f2(12), f2(13), f2(14), f2(15))
            nv = 5; nf = 15
         else
            call lyp_exc_vxc_unpolar(x(1), x(2), out(1), v(1), v(2))
            call lyp_fxc_unpolar(x(1), x(2), f(1), f(2), f(3))
            call lyp_exc_unpolar(x(1), x(2), e)
            call lyp_vxc_fxc_unpolar(x(1), x(2), v2(1), v2(2), f2(1), f2(2), f2(3))
            nv = 2; nf = 3
         end if
      case default
         print '(a)', "check_xc: no kernel named "//name
         stop 1
      end select

      out(2:1 + nv) = v(1:nv)
      out(2 + nv:1 + nv + nf) = f(1:nf)
      self = abs(e - out(1))/max(abs(out(1)), floor)
      self = max(self, maxval(abs(v2(1:nv) - v(1:nv))/max(abs(v(1:nv)), floor)))
      self = max(self, maxval(abs(f2(1:nf) - f(1:nf))/max(abs(f(1:nf)), floor)))
   end subroutine evaluate

end program check_xc
