!
! Body of one shell-quartet evaluation, with LA/LB/LC/LD supplied as PARAMETERS
! by the including procedure.
!
! Included once per angular-momentum class rather than written once and
! branched, because ptxas budgets registers for the worst path a kernel can
! take: a single generic routine costs every class the (dd|dd) register count.
! This is the same reason libERI ships int0000, int1111, int2222gen ... as
! separate entry points.
!
! Not generated code -- the body is written once and the compiler specialises
! it, which keeps the "no code generator" property of PLAN.md 1.0.2.
!
      lab = LA + LB; lcd = LC + LD; lt = lab + lcd
      na = ((LA+1)*(LA+2)/2); nb = ((LB+1)*(LB+2)/2); nc = ((LC+1)*(LC+2)/2); nd = ((LD+1)*(LD+2)/2)

      keyab = (q_i(iq) - 1)*nbas + q_j(iq)
      keycd = (q_k(iq) - 1)*nbas + q_l(iq)
      offab = pp_off(keyab); nab = pp_n(keyab)
      offcd = pp_off(keycd); ncd = pp_n(keycd)

      do kp = offab + 1, offab + nab
         p = pp_p(kp)
         do kq = offcd + 1, offcd + ncd
            q = pp_p(kq)
            alpha = p*q/(p + q)
            pqx = pp_r(kp, 1) - pp_r(kq, 1)
            pqy = pp_r(kp, 2) - pp_r(kq, 2)
            pqz = pp_r(kp, 3) - pp_r(kq, 3)
            pqd(1) = pqx; pqd(2) = pqy; pqd(3) = pqz
            tval = alpha*(pqx*pqx + pqy*pqy + pqz*pqz)

            call boys_eval(lt, tval, f)

            !
            ! R recurrence, PING-PONG over two buffers.
            !
            ! Level n is built only from level n+1, so the whole (LTOT+1)-level
            ! table was never needed -- two buffers are.  That cuts the scratch
            ! from NHERM_MAX*(LTOT+1) = 1485 doubles per work item to 2*165,
            ! i.e. 4.5x less footprint and, more to the point, 4.5x less DRAM
            ! traffic.  Measured to matter far more than the contraction did.
            !
            ! fac accumulates (-2 alpha)^n, which seeds R^n_000 at each level.
            !
            fac = 1.0_dp
            do n = 1, lt
               fac = fac*(-2.0_dp*alpha)
            end do
            cur = 0
            rscr(ir, cur*NHERM_MAX + 1) = fac*f(lt)

            do n = lt - 1, 0, -1
               nxt = 1 - cur
               fac = fac/(-2.0_dp*alpha)
               basen = nxt*NHERM_MAX
               base1 = cur*NHERM_MAX
               rscr(ir, basen + 1) = fac*f(n)
               deg = nherm_of(lt - n)
               do h = 2, deg
                  rscr(ir, basen + h) = pqd(hdir(h))*rscr(ir, base1 + hm1(h))
                  if (hm2(h) > 0) &
                     rscr(ir, basen + h) = rscr(ir, basen + h) &
                                           + hcf(h)*rscr(ir, base1 + hm2(h))
               end do
               cur = nxt
            end do
            rbase = cur*NHERM_MAX

            pref = TWO_PI_2_5/(p*q*sqrt(p + q))*pp_c(kp)*pp_c(kq)
            base = q_off(iq)
            nh_ab = nherm_of(lab)
            nh_cd = nherm_of(lcd)

            !
            ! TWO-PASS CONTRACTION.  The naive form evaluates
            !     sum_{tuv} sum_{TUV} E^ab E^cd R
            ! separately for every one of the ncart_ab * ncart_cd components,
            ! costing ncart_ab*ncart_cd*nherm_ab*nherm_cd.  Factoring it as
            !
            !     T(h_ab)      = sum_{h_cd} sgn E^cd(h_cd) R(h_ab + h_cd)
            !     (ab|cd)      = sum_{h_ab} E^ab(h_ab) T(h_ab)
            !
            ! costs nherm_ab*nherm_cd + ncart_ab*nherm_ab per ket component.
            ! For (dd|dd) that is 89.5k against 1.59M -- about 18x.
            !
            ! T is only nherm_ab long (<= 35), so it stays a small thread-local
            ! rather than needing scratch, which is the reason for doing one ket
            ! component at a time instead of forming the whole T matrix.
            !
            ! hshift() carries the index arithmetic, so the hot loop is a pure
            ! gather-multiply-add.
            !
            do icd = 0, nc*nd - 1
               id = icd/nc
               ic = icd - id*nc
               call cart_pow(LC, ic, cx, cy, cz)
               call cart_pow(LD, id, dx, dy, dz)

               ! E^cd over the flat Hermite index, sign folded in
               do h2 = 1, nh_cd
                  ecdv(h2) = 0.0_dp
               end do
               do h2 = 1, nh_cd
                  if (ht(h2) > cx + dx) cycle
                  if (hu(h2) > cy + dy) cycle
                  if (hv(h2) > cz + dz) cycle
                  ecdv(h2) = hsgn(h2) &
                             *pp_e(kq, e_off(cx, dx, ht(h2), 1)) &
                             *pp_e(kq, e_off(cy, dy, hu(h2), 2)) &
                             *pp_e(kq, e_off(cz, dz, hv(h2), 3))
               end do

               ! pass 1: T = R . E^cd
               do h1 = 1, nh_ab
                  tvec(h1) = 0.0_dp
               end do
               do h2 = 1, nh_cd
                  ecd = ecdv(h2)
                  if (ecd == 0.0_dp) cycle
                  do h1 = 1, nh_ab
                     tvec(h1) = tvec(h1) + ecd*rscr(ir, rbase + hshift(h1, h2))
                  end do
               end do

               ! pass 2: (ab|cd) = E^ab . T
               do iab = 0, na*nb - 1
                  ib = iab/na
                  ia = iab - ib*na
                  call cart_pow(LA, ia, ax, ay, az)
                  call cart_pow(LB, ib, bx, by, bz)
                  acc = 0.0_dp
                  do h1 = 1, nh_ab
                     if (ht(h1) > ax + bx) cycle
                     if (hu(h1) > ay + by) cycle
                     if (hv(h1) > az + bz) cycle
                     acc = acc + pp_e(kp, e_off(ax, bx, ht(h1), 1)) &
                                 *pp_e(kp, e_off(ay, by, hu(h1), 2)) &
                                 *pp_e(kp, e_off(az, bz, hv(h1), 3))*tvec(h1)
                  end do
                  out(base + 1 + iab + na*nb*icd) = &
                     out(base + 1 + iab + na*nb*icd) + pref*acc
               end do
            end do
         end do
      end do
