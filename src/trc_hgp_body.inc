!
! HGP quartet body, with LA/LB/LC/LD supplied as PARAMETERS by the including
! procedure.  Included once per angular-momentum class.
!
! THE POINT, and it is different from the MMD specialisation that failed:
! sizing.  The generic kernel declares its VRR workspace at NCUM*NCUM = 35x35,
! the worst case, so an (ss|ss) quartet carries 2450 doubles to hold 1.  With
! LA..LD known at compile time the workspace is exactly NCA*NCC -- 2 doubles
! for (ss|ss), 200 for (pp|pp), and only (dd|dd) pays the full 2450.
!
! That is a reduction in live state, not merely unrolled loops.  The earlier
! MMD attempt only got the latter and made things worse by raising registers.
!
! The VRR is inlined here rather than called, so its bounds are compile-time
! too.
!
      integer, parameter :: LAB = LA + LB, LCD = LC + LD, LT = LAB + LCD
      integer, parameter :: NCA = (LAB + 1)*(LAB + 2)*(LAB + 3)/6
      integer, parameter :: NCC = (LCD + 1)*(LCD + 2)*(LCD + 3)/6
      integer, parameter :: NA = (LA + 1)*(LA + 2)/2, NB = (LB + 1)*(LB + 2)/2
      integer, parameter :: NC = (LC + 1)*(LC + 2)/2, ND = (LD + 1)*(LD + 2)/2

      integer  :: keyab, keycd, offab, offcd, nab, ncd
      integer  :: kp, kq, m, d, x, cur, nxt, lim
      integer  :: ia, ic, xd1, xd2, xc1, xc2, xac, ai, dd
      real(dp) :: zeta, eta, zpe, rho, tval, pref, wc, acc
      real(dp) :: pa(3), qc(3), wp(3), wq(3), pq(3)
      real(dp) :: oo2z, oo2e, oo2ze, rz, re
      real(dp) :: f(0:BOYS_MMAX)
      real(dp) :: v(NCA*NCC, 0:1)
      real(dp) :: g(NCA*NCC)
      real(dp) :: abx, aby, abz, cdx, cdy, cdz
      integer  :: ib, id, idx, base
      integer  :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz
      integer  :: px, py, pz, qqx, qqy, qqz, ja, jc
      real(dp) :: cb, cdc

      keyab = (q_i(iq) - 1)*nbas + q_j(iq)
      keycd = (q_k(iq) - 1)*nbas + q_l(iq)
      offab = pp_off(keyab); nab = pp_n(keyab)
      offcd = pp_off(keycd); ncd = pp_n(keycd)

      do x = 1, NCA*NCC
         g(x) = 0.0_dp
      end do

      do kp = offab + 1, offab + nab
         zeta = pp_p(kp)
         do kq = offcd + 1, offcd + ncd
            eta = pp_p(kq)
            zpe = zeta + eta
            rho = zeta*eta/zpe
            do d = 1, 3
               pq(d) = pp_r(kp, d) - pp_r(kq, d)
               wc = (zeta*pp_r(kp, d) + eta*pp_r(kq, d))/zpe
               pa(d) = pp_r(kp, d) - pp_ra(kp, d)
               qc(d) = pp_r(kq, d) - pp_ra(kq, d)
               wp(d) = wc - pp_r(kp, d)
               wq(d) = wc - pp_r(kq, d)
            end do
            tval = rho*(pq(1)*pq(1) + pq(2)*pq(2) + pq(3)*pq(3))

            call boys_eval(LT, tval, f)

            oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = 0.5_dp/zpe
            rz = rho/zeta; re = rho/eta
            pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)

            cur = 0
            v(1, cur) = pref*f(LT)
            do m = LT - 1, 0, -1
               nxt = 1 - cur
               lim = LT - m
               v(1, nxt) = pref*f(m)

               do ia = 2, NCA
                  if (cll(ia) > lim) exit
                  dd = cdir(ia); xd1 = cdn1(ia); xd2 = cdn2(ia)
                  acc = pa(dd)*v(xd1, nxt) + wp(dd)*v(xd1, cur)
                  if (xd2 > 0) acc = acc + cf2(ia)*oo2z*(v(xd2, nxt) - rz*v(xd2, cur))
                  v(ia, nxt) = acc
               end do

               do ic = 2, NCC
                  if (cll(ic) > lim) exit
                  dd = cdir(ic); xc1 = cdn1(ic); xc2 = cdn2(ic)
                  do ia = 1, NCA
                     if (cll(ia) + cll(ic) > lim) exit
                     x = ia + (ic - 1)*NCA
                     acc = qc(dd)*v(ia + (xc1 - 1)*NCA, nxt) &
                           + wq(dd)*v(ia + (xc1 - 1)*NCA, cur)
                     if (xc2 > 0) acc = acc + cf2(ic)*oo2e* &
                        (v(ia + (xc2 - 1)*NCA, nxt) - re*v(ia + (xc2 - 1)*NCA, cur))
                     ai = 0; xac = 0
                     if (dd == 1) then
                        ai = cnx(ia)
                        if (ai > 0) xac = cidx(cnx(ia) - 1, cny(ia), cnz(ia))
                     else if (dd == 2) then
                        ai = cny(ia)
                        if (ai > 0) xac = cidx(cnx(ia), cny(ia) - 1, cnz(ia))
                     else
                        ai = cnz(ia)
                        if (ai > 0) xac = cidx(cnx(ia), cny(ia), cnz(ia) - 1)
                     end if
                     if (xac > 0) acc = acc + real(ai, dp)*oo2ze*v(xac + (xc1 - 1)*NCA, cur)
                     v(x, nxt) = acc
                  end do
               end do
               cur = nxt
            end do

            do x = 1, NCA*NCC
               g(x) = g(x) + v(x, cur)
            end do
         end do
      end do

      ! --- horizontal transfer, closed form, on contracted data ---
      abx = pp_ra(offab + 1, 1) - pp_rb(offab + 1, 1)
      aby = pp_ra(offab + 1, 2) - pp_rb(offab + 1, 2)
      abz = pp_ra(offab + 1, 3) - pp_rb(offab + 1, 3)
      cdx = pp_ra(offcd + 1, 1) - pp_rb(offcd + 1, 1)
      cdy = pp_ra(offcd + 1, 2) - pp_rb(offcd + 1, 2)
      cdz = pp_ra(offcd + 1, 3) - pp_rb(offcd + 1, 3)

      base = q_off(iq)
      idx = 0
      do id = 0, ND - 1
         call comp(LD, id, dx, dy, dz)
      do ic = 0, NC - 1
         call comp(LC, ic, cx, cy, cz)
      do ib = 0, NB - 1
         call comp(LB, ib, bx, by, bz)
      do ia = 0, NA - 1
         call comp(LA, ia, ax, ay, az)
         idx = idx + 1
         acc = 0.0_dp
         do pz = 0, bz
         do py = 0, by
         do px = 0, bx
            cb = binom(bx, px)*binom(by, py)*binom(bz, pz) &
                 *powi(abx, bx - px)*powi(aby, by - py)*powi(abz, bz - pz)
            if (cb == 0.0_dp) cycle
            ja = cidx(ax + px, ay + py, az + pz)
            do qqz = 0, dz
            do qqy = 0, dy
            do qqx = 0, dx
               cdc = binom(dx, qqx)*binom(dy, qqy)*binom(dz, qqz) &
                     *powi(cdx, dx - qqx)*powi(cdy, dy - qqy)*powi(cdz, dz - qqz)
               if (cdc == 0.0_dp) cycle
               jc = cidx(cx + qqx, cy + qqy, cz + qqz)
               acc = acc + cb*cdc*g(ja + (jc - 1)*NCA)
            end do
            end do
            end do
         end do
         end do
         end do
         out(base + idx) = out(base + idx) + acc
      end do
      end do
      end do
      end do
