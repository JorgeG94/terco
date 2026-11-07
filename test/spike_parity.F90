!
! PARITY-SKIP SPIKE -- NOT FOR MERGING.  docs/SYMMETRY_PARITY_SKIP.md.
!
! CPU prototype in the style of the S13 spike_symfock.F90, same bins.
!
! For a molecule/basis it builds the terco bins (eri%bins), closes them under
! the point group, enumerates the quartets of two plans
!   full   the production enumeration (unordered bin pairs, triangle on the
!          diagonal ones)
!   s2a    the skeleton plan of S14.2 (bin-pair level fold kept, one side cut
!          to the orbit representatives of its bins)
! and for every quartet that survives the exact Schwarz product tallies, per
! class (la,lb,lc,ld; canonical order), per set Z of "dead" coordinates (the
! coordinates on which all four centres vanish) and per stabilizer S of the
! quartet (the operations of the group that map each of the four shells to
! itself, as a set of sign-flip vectors), the number of quartets and of
! primitive quartets.  `count` writes that table as CSV; scripts/
! spike_parity_gen.py turns it into the zero fractions and the cost-weighted
! work with its operation counts.
!
! `verify` evaluates every quartet with hgp_batch (production routine,
! Cartesian working basis), checks that every component that the stabilizer
! rule calls zero IS zero, and builds G = J - K/2 from a symmetric density
! three ways: unskipped, with the odd-parity components left out of the
! digestion, and the s2a skeleton with the skip, symmetrized, against the
! unskipped full G (the S14.2 gate: 5e-11 * max(1, max|G|)).
!
! usage: spike_parity xyz basis.json mode(count|verify|both) thresh csv
!
module spike_parity_mod
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t
   use trc_hgp, only: hgp_batch, hgp_scratch_size
   use trc_cart_perm, only: cart_perm_powers
   implicit none

   type :: seg_t
      integer :: ia = 0, ib = 0, ka = 1, kb = 1
      logical :: tri = .false., wosz = .false.
      real(dp) :: scale = 1.0_dp
   end type seg_t

   type :: ctx_t
      integer :: nop = 1, nsh = 0, nw = 0, np = 0, nlb = 0, natm = 0
      real(dp) :: thresh = 1.0e-10_dp
      integer :: smax_keep = 0
      integer, allocatable :: pi(:), pj(:), sh_l(:), sh_ao(:), sh_at(:), opmask(:)
      real(dp), allocatable :: pq(:)
      integer, allocatable :: pmap(:, :), img(:, :), rep(:), osz(:), sm(:, :), bs(:)
      integer, allocatable :: aoff(:), alist(:), roff(:), rlist(:)
      integer, allocatable :: zatom(:)          !! dead-coordinate mask of each atom
      integer, allocatable :: npp_of(:)         !! primitive pairs of each shell pair
      type(trc_eri_t), pointer :: eri => null()
      ! tallies (la, lb, lc, ld, Z, S)
      integer(int64), allocatable :: tq(:, :, :, :, :, :), tp(:, :, :, :, :, :)
      integer(int64) :: nquart = 0
      ! verify
      logical :: dig = .false.
      integer :: skip = 0         !! 0 none, 1 stabilizer rule, 2 dead-coordinate rule
      real(dp), allocatable :: dw(:, :), jx(:, :), kx(:, :), out(:), work(:)
      real(dp) :: maxflag = 0.0_dp         !! largest |v| of a component the rule calls zero
      integer(int64) :: ncomp = 0, nflag = 0
   end type ctx_t

   type(ctx_t), save, target :: c
   integer, parameter :: LM = 4

contains

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   subroutine setup(wb, aomap, geomperm, ops, nop, rstd, eri, thresh)
      type(trc_basis_t), intent(in) :: wb
      integer, intent(in) :: aomap(:, :), geomperm(:, :), ops(:), nop
      real(dp), intent(in) :: rstd(:, :)
      type(trc_eri_t), target, intent(inout) :: eri
      real(dp), intent(in) :: thresh
      integer :: s, g, k, p, a, t, nlive, cnt, nb
      integer, allocatable :: shell_of_ao(:), binorig(:), tmp(:), lbin(:), pa(:), pr(:)
      logical, allocatable :: seen(:)
      real(dp) :: qm

      c%eri => eri
      c%thresh = thresh
      c%nop = nop
      c%opmask = ops(1:nop)
      c%nsh = wb%nshell; c%nw = wb%nao; c%np = eri%bins%npair; c%natm = wb%natm
      c%sh_l = wb%sh_l; c%sh_ao = wb%sh_ao; c%sh_at = wb%sh_at
      c%pi = eri%bins%sp_i(1:c%np); c%pj = eri%bins%sp_j(1:c%np)
      c%pq = eri%bins%sp_q(1:c%np)
      c%smax_keep = int(-log10(thresh))*eri%bins%sres
      nb = c%nsh
      allocate (c%npp_of(c%np))
      do p = 1, c%np
         c%npp_of(p) = eri%pairs%n((c%pi(p) - 1)*nb + c%pj(p))
      end do
      allocate (c%zatom(c%natm))
      do a = 1, c%natm
         c%zatom(a) = 0
         do k = 1, 3
            if (abs(rstd(k, a)) < 1.0e-6_dp) c%zatom(a) = ior(c%zatom(a), ishft(1, k - 1))
         end do
      end do
      allocate (shell_of_ao(c%nw), c%sm(c%nsh, nop))
      do s = 1, c%nsh
         do a = 0, ncart(c%sh_l(s)) - 1
            shell_of_ao(c%sh_ao(s) + a) = s
         end do
      end do
      do g = 1, nop
         do s = 1, c%nsh
            c%sm(s, g) = shell_of_ao(aomap(c%sh_ao(s), g))
         end do
      end do
      allocate (c%pmap(c%nsh, c%nsh), c%img(c%np, nop), c%rep(c%np), c%osz(c%np))
      c%pmap = 0
      do p = 1, c%np
         c%pmap(c%pi(p), c%pj(p)) = p
         c%pmap(c%pj(p), c%pi(p)) = p
      end do
      do p = 1, c%np
         do g = 1, nop
            k = c%pmap(c%sm(c%pi(p), g), c%sm(c%pj(p), g))
            if (k == 0) k = p
            c%img(p, g) = k
         end do
      end do
      do p = 1, c%np
         c%rep(p) = minval(c%img(p, :))
         cnt = 0
         do g = 1, nop
            if (.not. any(c%img(p, 1:g - 1) == c%img(p, g))) cnt = cnt + 1
         end do
         c%osz(p) = cnt
      end do
      do p = 1, c%np
         qm = 0.0_dp
         do g = 1, nop
            qm = max(qm, eri%bins%sp_q(c%img(p, g)))
         end do
         c%pq(p) = qm
      end do
      allocate (binorig(c%np))
      do k = 1, eri%bins%nbin
         do p = eri%bins%bin_off(k) + 1, eri%bins%bin_off(k) + eri%bins%bin_cnt(k)
            binorig(p) = k
         end do
      end do
      allocate (tmp(c%np))
      do p = 1, c%np
         tmp(p) = binorig(c%rep(p))
      end do
      allocate (seen(eri%bins%nbin)); seen = .false.
      do p = 1, c%np
         seen(tmp(p)) = .true.
      end do
      nlive = count(seen)
      c%nlb = nlive
      allocate (c%bs(nlive), c%aoff(nlive + 1), c%roff(nlive + 1), c%alist(c%np), c%rlist(c%np))
      allocate (lbin(eri%bins%nbin)); lbin = 0
      t = 0
      do k = 1, eri%bins%nbin
         if (seen(k)) then
            t = t + 1; lbin(k) = t; c%bs(t) = eri%bins%bin_s(k)
         end if
      end do
      c%aoff = 0; c%roff = 0
      do p = 1, c%np
         c%aoff(lbin(tmp(p)) + 1) = c%aoff(lbin(tmp(p)) + 1) + 1
         if (c%rep(p) == p) c%roff(lbin(tmp(p)) + 1) = c%roff(lbin(tmp(p)) + 1) + 1
      end do
      do t = 1, nlive
         c%aoff(t + 1) = c%aoff(t + 1) + c%aoff(t)
         c%roff(t + 1) = c%roff(t + 1) + c%roff(t)
      end do
      allocate (pa(nlive), pr(nlive))
      pa = c%aoff(1:nlive); pr = c%roff(1:nlive)
      do p = 1, c%np
         t = lbin(tmp(p))
         pa(t) = pa(t) + 1; c%alist(pa(t)) = p
         if (c%rep(p) == p) then
            pr(t) = pr(t) + 1; c%rlist(pr(t)) = p
         end if
      end do
      allocate (c%tq(0:LM, 0:LM, 0:LM, 0:LM, 0:7, 0:255), c%tp(0:LM, 0:LM, 0:LM, 0:LM, 0:7, 0:255))
      allocate (c%out(4096), c%work(4096))
   end subroutine setup

   pure integer function listp(kind, lb, i) result(p)
      integer, intent(in) :: kind, lb, i
      if (kind == 1) then
         p = c%alist(c%aoff(lb) + i)
      else
         p = c%rlist(c%roff(lb) + i)
      end if
   end function listp

   pure integer function listn(kind, lb) result(n)
      integer, intent(in) :: kind, lb
      if (kind == 1) then
         n = c%aoff(lb + 1) - c%aoff(lb)
      else
         n = c%roff(lb + 1) - c%roff(lb)
      end if
   end function listn

   subroutine plan(scheme, segs, nseg)
      !! scheme 1: full, scheme 3: s2a (S13 numbering)
      integer, intent(in) :: scheme
      type(seg_t), allocatable, intent(out) :: segs(:)
      integer, intent(out) :: nseg
      type(seg_t), allocatable :: tmp(:)
      integer :: ia, ib, na, nb, nra, nrb
      allocate (tmp(c%nlb*(c%nlb + 1)/2 + 2))
      nseg = 0
      do ia = 1, c%nlb
         do ib = 1, ia
            if (c%bs(ia) + c%bs(ib) > c%smax_keep) cycle
            na = listn(1, ia); nb = listn(1, ib)
            nra = listn(2, ia); nrb = listn(2, ib)
            if (scheme == 1) then
               nseg = nseg + 1; tmp(nseg) = seg_t(ia, ib, 1, 1, ia == ib, .false., 1.0_dp)
            else
               nseg = nseg + 1
               if (ia == ib) then
                  tmp(nseg) = seg_t(ia, ib, 2, 1, .false., .true., 0.5_dp)
               else if (int(nra, int64)*nb <= int(na, int64)*nrb) then
                  tmp(nseg) = seg_t(ia, ib, 2, 1, .false., .true., 1.0_dp)
               else
                  tmp(nseg) = seg_t(ib, ia, 2, 1, .false., .true., 1.0_dp)
               end if
            end if
         end do
      end do
      segs = tmp(1:nseg)
   end subroutine plan

   subroutine walk(scheme, seg)
      integer, intent(in) :: scheme
      type(seg_t), intent(in) :: seg
      integer :: i, j, p, q, jhi
      do i = 1, listn(seg%ka, seg%ia)
         p = listp(seg%ka, seg%ia, i)
         jhi = listn(seg%kb, seg%ib)
         if (seg%tri) jhi = i
         do j = 1, jhi
            q = listp(seg%kb, seg%ib, j)
            call visit(seg, p, q)
         end do
      end do
   end subroutine walk

   subroutine quartet_masks(p, q, zmask, smask)
      !! Z: coordinates dead for all four centres; S: bit v set when the
      !! operation with flip vector v maps all four shells to themselves.
      integer, intent(in) :: p, q
      integer, intent(out) :: zmask, smask
      integer :: g, s4(4)
      s4 = [c%pi(p), c%pj(p), c%pi(q), c%pj(q)]
      zmask = iand(iand(c%zatom(c%sh_at(s4(1))), c%zatom(c%sh_at(s4(2)))), &
                   iand(c%zatom(c%sh_at(s4(3))), c%zatom(c%sh_at(s4(4)))))
      smask = 0
      do g = 1, c%nop
         if (all(c%sm(s4, g) == s4)) smask = ior(smask, ishft(1, c%opmask(g)))
      end do
   end subroutine quartet_masks

   subroutine visit(seg, p, q)
      type(seg_t), intent(in) :: seg
      integer, intent(in) :: p, q
      integer :: l(4), pl(2, 2), z, s
      real(dp) :: w
      if (c%pq(p)*c%pq(q) <= c%thresh) return
      c%nquart = c%nquart + 1
      call quartet_masks(p, q, z, s)
      pl(:, 1) = [max(c%sh_l(c%pi(p)), c%sh_l(c%pj(p))), min(c%sh_l(c%pi(p)), c%sh_l(c%pj(p)))]
      pl(:, 2) = [max(c%sh_l(c%pi(q)), c%sh_l(c%pj(q))), min(c%sh_l(c%pi(q)), c%sh_l(c%pj(q)))]
      if (pl(1, 1) < pl(1, 2) .or. (pl(1, 1) == pl(1, 2) .and. pl(2, 1) < pl(2, 2))) then
         l = [pl(:, 2), pl(:, 1)]
      else
         l = [pl(:, 1), pl(:, 2)]
      end if
      c%tq(l(1), l(2), l(3), l(4), z, s) = c%tq(l(1), l(2), l(3), l(4), z, s) + 1
      c%tp(l(1), l(2), l(3), l(4), z, s) = c%tp(l(1), l(2), l(3), l(4), z, s) &
                                          + int(c%npp_of(p), int64)*c%npp_of(q)
      if (.not. c%dig) return
      w = 1.0_dp
      if (seg%wosz) w = real(c%osz(p), dp)
      w = w*seg%scale
      call digest(p, q, w, z, s, (p == q .and. .not. seg%wosz))
   end subroutine visit

   pure logical function rejected(smask, par)
      !! The stabilizer rule: the component of parity vector `par` (bit d =
      !! odd power of coordinate d) is odd under some operation of the set.
      integer, intent(in) :: smask, par
      integer :: v
      rejected = .false.
      do v = 0, 7
         if (btest(smask, v)) then
            if (mod(popcnt(iand(v, par)), 2) == 1) rejected = .true.
         end if
      end do
   end function rejected

   subroutine digest(p, q, w, zmask, smask, half)
      integer, intent(in) :: p, q, zmask, smask
      real(dp), intent(in) :: w
      logical, intent(in) :: half
      integer :: si, sj, sk, sl, la, lb, lc, ld, na, nb, nc, nd, nout, idx
      integer :: ia, ib, ic, id, o(4), x(4), r, aoff4(4), px, py, pz, zs
      integer :: pa(15), pb(15), pc(15), pd(15)
      integer(int64) :: nwk
      integer :: qi(1), qj(1), qk(1), ql(1), qo(1)
      integer(int64) :: wo(1)
      real(dp) :: f, v
      logical :: flag
      integer, parameter :: ORD(4, 8) = reshape([ &
         1, 2, 3, 4,  2, 1, 3, 4,  1, 2, 4, 3,  2, 1, 4, 3, &
         3, 4, 1, 2,  4, 3, 1, 2,  3, 4, 2, 1,  4, 3, 2, 1], [4, 8])

      si = c%pi(p); sj = c%pj(p); sk = c%pi(q); sl = c%pj(q)
      la = c%sh_l(si); lb = c%sh_l(sj); lc = c%sh_l(sk); ld = c%sh_l(sl)
      na = ncart(la); nb = ncart(lb); nc = ncart(lc); nd = ncart(ld)
      nout = na*nb*nc*nd
      nwk = hgp_scratch_size(la + lb, lc + ld)
      if (size(c%out) < nout) then
         deallocate (c%out); allocate (c%out(2*nout))
      end if
      if (size(c%work) < nwk) then
         deallocate (c%work); allocate (c%work(2*nwk))
      end if
      c%out(1:nout) = 0.0_dp
      qi = si; qj = sj; qk = sk; ql = sl; qo = 0; wo = 0
      call hgp_batch(1, 1, 1, c%nsh, c%eri%pairs%npp, nout, qi, qj, qk, ql, qo, c%sh_l, &
                     c%eri%pairs%off, c%eri%pairs%n, c%eri%pairs%p, c%eri%pairs%r, &
                     c%eri%pairs%ra, c%eri%pairs%rb, c%eri%pairs%c, c%out(1:nout), &
                     nwk, wo, c%work(1:nwk))
      f = w
      if (si == sj) f = 0.5_dp*f
      if (sk == sl) f = 0.5_dp*f
      if (half) f = 0.5_dp*f
      aoff4 = [c%sh_ao(si), c%sh_ao(sj), c%sh_ao(sk), c%sh_ao(sl)] - 1
      do ia = 1, na
         call cart_perm_powers(la, ia - 1, px, py, pz)
         pa(ia) = ior(iand(px, 1), ior(ishft(iand(py, 1), 1), ishft(iand(pz, 1), 2)))
      end do
      do ia = 1, nb
         call cart_perm_powers(lb, ia - 1, px, py, pz)
         pb(ia) = ior(iand(px, 1), ior(ishft(iand(py, 1), 1), ishft(iand(pz, 1), 2)))
      end do
      do ia = 1, nc
         call cart_perm_powers(lc, ia - 1, px, py, pz)
         pc(ia) = ior(iand(px, 1), ior(ishft(iand(py, 1), 1), ishft(iand(pz, 1), 2)))
      end do
      do ia = 1, nd
         call cart_perm_powers(ld, ia - 1, px, py, pz)
         pd(ia) = ior(iand(px, 1), ior(ishft(iand(py, 1), 1), ishft(iand(pz, 1), 2)))
      end do
      zs = 0
      do ia = 0, 2
         if (btest(zmask, ia)) zs = ior(zs, ishft(1, ishft(1, ia)))   ! flip vector 1<<ia
      end do
      idx = 0
      do id = 0, nd - 1
         do ic = 0, nc - 1
            do ib = 0, nb - 1
               do ia = 0, na - 1
                  idx = idx + 1
                  v = f*c%out(idx)
                  c%ncomp = c%ncomp + 1
                  flag = .false.
                  select case (c%skip)
                  case (1); flag = rejected(smask, ieor(ieor(pa(ia + 1), pb(ib + 1)), ieor(pc(ic + 1), pd(id + 1))))
                  case (2); flag = rejected(zs, ieor(ieor(pa(ia + 1), pb(ib + 1)), ieor(pc(ic + 1), pd(id + 1))))
                  end select
                  if (flag) then
                     c%nflag = c%nflag + 1
                     c%maxflag = max(c%maxflag, abs(c%out(idx)))
                     cycle
                  end if
                  if (v == 0.0_dp) cycle
                  o = [ia, ib, ic, id] + aoff4 + 1
                  do r = 1, 8
                     x = o(ORD(:, r))
                     c%jx(x(1), x(2)) = c%jx(x(1), x(2)) + v*c%dw(x(3), x(4))
                     c%kx(x(1), x(3)) = c%kx(x(1), x(3)) + v*c%dw(x(2), x(4))
                  end do
               end do
            end do
         end do
      end do
   end subroutine digest

   subroutine run_plan(scheme, dig, skip)
      integer, intent(in) :: scheme, skip
      logical, intent(in) :: dig
      type(seg_t), allocatable :: segs(:)
      integer :: nseg, s
      c%dig = dig; c%skip = skip
      c%nquart = 0; c%tq = 0; c%tp = 0; c%ncomp = 0; c%nflag = 0; c%maxflag = 0
      if (dig) then
         if (allocated(c%jx)) deallocate (c%jx, c%kx)
         allocate (c%jx(c%nw, c%nw), c%kx(c%nw, c%nw)); c%jx = 0; c%kx = 0
      end if
      call plan(scheme, segs, nseg)
      do s = 1, nseg
         call walk(scheme, segs(s))
      end do
   end subroutine run_plan

   subroutine symmetrize(x, aomap, aosign, y)
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: aomap(:, :), aosign(:, :)
      real(dp), intent(out) :: y(:, :)
      integer :: g, a, b, n
      n = size(x, 1)
      y = 0.0_dp
      do g = 1, size(aomap, 2)
         do b = 1, n
            do a = 1, n
               y(aomap(a, g), aomap(b, g)) = y(aomap(a, g), aomap(b, g)) &
                                             + real(aosign(a, g)*aosign(b, g), dp)*x(a, b)
            end do
         end do
      end do
      y = y/real(size(aomap, 2), dp)
   end subroutine symmetrize

   subroutine write_csv(unit, label)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      integer :: a, b, cc, d, z, s
      do s = 0, 255
         do z = 0, 7
            do d = 0, LM
               do cc = 0, LM
                  do b = 0, LM
                     do a = 0, LM
                        if (c%tq(a, b, cc, d, z, s) == 0) cycle
                        write (unit, '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0)') &
                           trim(label), a, b, cc, d, z, s, c%tq(a, b, cc, d, z, s), c%tp(a, b, cc, d, z, s)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine write_csv

end module spike_parity_mod

program spike_parity
   use pic_types, only: dp, int64
   use spike_parity_mod
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect
   use trc_test_basis, only: read_xyz
   implicit none
   character(len=256) :: xyz, bfile, mode, csv, arg
   real(dp) :: thresh
   integer :: nat, i, u, nw, g, h, ne, pl
   integer, allocatable :: elem(:)
   real(dp), allocatable :: r(:, :), dpw(:, :), t(:, :), jf(:, :), kf(:, :), gf(:, :), js(:, :), ks(:, :), gs(:, :), gy(:, :)
   type(trc_basis_t), target :: bas
   type(trc_basis_t) :: wbc
   type(trc_sym_geom_t) :: geom
   type(trc_sym_ao_t) :: aow
   type(trc_eri_t), target :: eri
   type(error_t) :: err
   type(trc_basis_t), pointer :: wb
   logical :: do_count, do_verify
   character(len=5), parameter :: PNAME(2) = [character(len=5) :: 'full', 's2a']
   integer, parameter :: PSCHEME(2) = [1, 3]
   integer(int64) :: nq(2)
   real(dp) :: gate, scl

   i = trc_bind_device(0)
   call get_command_argument(1, xyz)
   call get_command_argument(2, bfile)
   call get_command_argument(3, mode)
   call get_command_argument(4, arg); read (arg, *) thresh
   call get_command_argument(5, csv)
   do_count = index(mode, 'count') > 0 .or. index(mode, 'both') > 0
   do_verify = index(mode, 'verify') > 0 .or. index(mode, 'both') > 0

   call read_xyz(trim(xyz), nat, elem, r)
   geom = sym_detect(r, elem, err, basis_id=elem)
   if (err%has_error()) stop 'sym_detect failed'
   print '(a,a,a,i0,a,a)', 'molecule ', trim(xyz), '  natm ', nat, '  group ', trim(geom%group%name)
   call trc_basis_from_json(trim(bfile), nat, elem, geom%r_std, bas, err, force_cart=.true.)
   if (err%has_error()) stop 'basis failed'
   bas%spherical = .false.
   call bas%to_device()
   call eri%build(bas, thresh)
   if (eri%dec%active) then
      wb => eri%pb
   else
      wb => bas
   end if
   wbc = wb
   wbc%spherical = .false.
   aow = sym_ao_build(wbc, geom, err)
   if (err%has_error()) stop 'sym_ao_build failed'
   nw = wb%nao
   call setup(wb, aow%ao_map, geom%perm, geom%group%ops, geom%group%order, geom%r_std, eri, thresh)
   print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,l1)', 'nop ', c%nop, '  shells ', c%nsh, '  working nao ', nw, '  pairs ', c%np, &
      '  live bins ', c%nlb, '  decontracted ', eri%dec%active
   ne = 0
   do g = 1, c%nop
      if (.not. any([(all(c%sm(:, g) == c%sm(:, h)), h=1, g - 1)])) ne = ne + 1
   end do
   print '(a,i0,a,i0)', 'effective group order at shell level ', ne, ' of |G| = ', c%nop
   print '(a,*(i0,1x))', 'atoms with dead x/y/z coordinates (counts) ', count(iand(c%zatom, 1) /= 0), &
      count(iand(c%zatom, 2) /= 0), count(iand(c%zatom, 4) /= 0)

   if (do_count) then
      open (newunit=u, file=trim(csv), status='replace')
      write (u, '(a)') 'plan,la,lb,lc,ld,Z,S,nq,nprim'
      do pl = 1, 2
         call run_plan(PSCHEME(pl), .false., 0)
         nq(pl) = c%nquart
         call write_csv(u, trim(PNAME(pl)))
         print '(a,a,a,i0)', 'plan ', trim(PNAME(pl)), ': quartets (after Schwarz) ', c%nquart
      end do
      close (u)
      print '(a,f6.3)', 'full / s2a quartets ', real(nq(1), dp)/real(nq(2), dp)
   end if

   if (.not. do_verify) stop
   ! a symmetric density (deterministic), scaled like a density matrix
   allocate (t(nw, nw), dpw(nw, nw))
   do i = 1, nw
      do g = 1, nw
         t(i, g) = 0.3_dp*cos(real(2*i + 3*g, dp)) + 0.3_dp*cos(real(2*g + 3*i, dp)) + merge(1.0_dp, 0.0_dp, i == g)
      end do
   end do
   call symmetrize(t, aow%ao_map, aow%ao_sign, dpw)
   c%dw = dpw
   allocate (jf(nw, nw), kf(nw, nw), gf(nw, nw), js(nw, nw), ks(nw, nw), gs(nw, nw), gy(nw, nw))
   call run_plan(1, .true., 0)
   jf = c%jx; kf = c%kx; gf = jf - 0.5_dp*kf
   print '(a,i0,a,es10.3)', 'VERIFY full unskipped: quartets ', c%nquart, '   max|G| ', maxval(abs(gf))
   scl = max(1.0_dp, maxval(abs(gf)))
   gate = 5.0e-11_dp*scl
   do i = 1, 2
      call run_plan(1, .true., i)
      gs = c%jx - 0.5_dp*c%kx
      print '(a,i0,a,f7.3,a,es10.2,a,es10.2,a,es10.2)', '  full + skip rule ', i, ': skipped fraction ', &
         real(c%nflag, dp)/real(c%ncomp, dp), '  max|v| of a skipped component ', c%maxflag, &
         '  max|G-G_full| ', maxval(abs(gs - gf)), '  gate ', gate
      call run_plan(3, .true., i)
      gs = c%jx - 0.5_dp*c%kx
      call symmetrize(gs, aow%ao_map, aow%ao_sign, gy)
      print '(a,i0,a,f7.3,a,es10.2,a,l1)', '  s2a  + skip rule ', i, ': skipped fraction ', &
         real(c%nflag, dp)/real(c%ncomp, dp), '  max|sym(G_s2a)-G_full| ', maxval(abs(gy - gf)), &
         '  within gate ', maxval(abs(gy - gf)) <= gate
   end do
   call run_plan(3, .true., 0)
   gs = c%jx - 0.5_dp*c%kx
   call symmetrize(gs, aow%ao_map, aow%ao_sign, gy)
   print '(a,es10.2)', '  s2a unskipped (S14.2 reference): max|sym(G_s2a)-G_full| ', maxval(abs(gy - gf))
end program spike_parity
