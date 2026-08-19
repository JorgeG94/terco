!
! Dummy-argument and local declarations shared by every class-specialised
! quartet kernel.  LA/LB/LC/LD are NOT here -- each including procedure
! declares them as parameters, which is the whole point.
!
      !$acc routine seq
      integer,  intent(in)    :: iq, ir, nchunk, nq, nbas, npp, nout
      integer,  intent(in)    :: q_i(nq), q_j(nq), q_k(nq), q_l(nq), q_off(nq)
      integer,  intent(in)    :: sh_l(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_c(npp), pp_e(npp, NE_PAIR)
      real(dp), intent(inout) :: rscr(nchunk, 2*NHERM_MAX)
      real(dp), intent(inout) :: out(nout)

      integer  :: lab, lcd, lt
      integer  :: na, nb, nc, nd
      integer  :: keyab, keycd, offab, offcd, nab, ncd
      integer  :: kp, kq, n, h, deg, base, basen, base1, cur, nxt, rbase
      integer  :: ia, ib, ic, id, idx
      integer  :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz
      integer  :: h1, h2, nh_ab, nh_cd, iab, icd
      real(dp) :: p, q, alpha, tval, pref, acc, ecd, fac
      ! Small thread-locals, fixed extent (never sized by a dummy -- nvfortran
      ! crashes on that in a device routine).  Sized by NHERM_PAIR, the degree
      ! ONE PAIR reaches (35 at LMAX=2), not the quartet's NHERM_MAX (165):
      ! these are indexed by a pair's Hermite index and the wrong bound cost
      ! 2.6 kB of stack frame per thread.
      real(dp) :: tvec(NHERM_PAIR), ecdv(NHERM_PAIR)
      real(dp) :: pqx, pqy, pqz, pqd(3)
      real(dp) :: f(0:BOYS_MMAX)

