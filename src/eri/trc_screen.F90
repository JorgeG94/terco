!
! Screening and the compacted, symmetry-reduced work list.
!
! Three reductions, in this order because each makes the next cheaper:
!
!   1. PERMUTATIONAL SYMMETRY, 8-fold: only canonical quartets are stored
!      (i >= j, k >= l, pair(i,j) >= pair(k,l)).
!
!   2. SCHWARZ: |(ij|kl)| <= Q(ij) Q(kl), Q(ij) = sqrt(max |(ij|ij)|).
!      Rigorous (Cauchy-Schwarz), so discarding below a threshold bounds the
!      error. Q is built with the full Coulomb operator and never adjusted for
!      the range-separation omega: in k-space erf(omega r)/r has transform
!      4*pi/k^2 * exp(-k^2/(4 omega^2)) <= 4*pi/k^2, so (ab|ab)_erf <=
!      (ab|ab)_Coulomb for any charge distribution. A tighter omega-specific
!      bound would only screen more, so this costs no accuracy. The prescreens
!      in trc_batch.F90 and the generated kernels' `pcut` rest on the same fact.
!
!   3. DENSITY: a quartet's Fock contribution, not the integral's size, is what
!      matters, so the bound carries the largest density block the quartet can
!      touch. This makes late SCF iterations cheap and is redone whenever D
!      changes.
!
! The list is built once on the host and only real quartets reach the device;
! screening inside a kernel would not help, since an early-exiting GPU thread
! still holds its warp for the busiest lane.
!
! Q comes from our own engine on the diagonal quartets (O(N^2)), not libfint,
! which is gfortran-only (see NOTES.md). libfint remains the oracle that
! checks the bounds, in check_screen.
!
module trc_screen
   use pic_types, only: dp
   use trc_tables, only: LTOT, NHERM_MAX
   use trc_batch, only: eri_batch, ncart
   use pic_types, only: int64
   use trc_hgp, only: hgp_batch, hgp_scratch_size
   implicit none
   private

   public :: schwarz_bounds, schwarz_bounds_hgp, pair_index, build_worklist, density_blockmax

contains

   !
   ! Canonical compound index of a shell pair, i >= j, 1-based.
   !
   pure integer function pair_index(i, j)
      !$acc routine seq
      integer, intent(in) :: i, j
      pair_index = i*(i - 1)/2 + j
   end function pair_index

   !
   ! Q(i,j) = sqrt( max_{ab} |(ij|ij)_{ab,ab}| ) for every i >= j.
   !
   ! Runs the ordinary engine over the diagonal quartets.  N(N+1)/2 of them, so
   ! this is negligible against the N^4 it then removes.
   !
   subroutine schwarz_bounds(nbas, npp, sh_l, pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, q)
      integer,  intent(in)  :: nbas, npp
      integer,  intent(in)  :: sh_l(nbas)
      integer,  intent(in)  :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)  :: pp_p(npp), pp_r(npp, 3), pp_c(npp)
      real(dp), intent(in)  :: pp_e(npp, *)
      real(dp), intent(out) :: q(nbas*(nbas + 1)/2)

      integer, allocatable :: d_i(:), d_j(:), d_k(:), d_l(:), d_off(:), d_q(:)
      real(dp), allocatable :: out(:), rscr(:, :), qc(:)
      integer :: npair, i, j, n, nout, nb, m, na, nlive, c0, c1, nc, key
      real(dp) :: vmax
      ! One thread per shell pair, so the chunk IS the launch width. At 8192
      ! a V100 ran a fraction of one wave and this routine cost more than the
      ! Fock build it exists to accelerate. The scratch is 2*NHERM_MAX doubles
      ! per pair, about 2.6 kB, so 65536 pairs is ~170 MB -- affordable, and
      ! four launches instead of twenty-five.
      integer, parameter :: CHUNK = 65536

      !
      ! ONLY THE PAIRS THAT HAVE PRIMITIVES, AND ONLY A CHUNK AT A TIME.
      !
      ! This used to enumerate every one of the nbas(nbas+1)/2 pairs and
      ! allocate `out` and a Hermite scratch of npair x 2*NHERM_MAX for all
      ! of them at once, zero both on the HOST and let them cross the bus:
      ! on a 123-atom slice in cc-pVDZ that is 200 thousand pairs and half
      ! a gigabyte of scratch, for a result that is one number per pair.
      ! Pairs whose primitives were all dropped contribute nothing and are
      ! skipped outright; the rest go in chunks, with the scratch created on
      ! the device and never copied, and the reduction to q done there too,
      ! so only the bound itself comes back.
      !
      npair = nbas*(nbas + 1)/2
      q = 0.0_dp
      allocate (d_i(npair), d_j(npair), d_k(npair), d_l(npair), d_off(npair), d_q(npair))
      nlive = 0
      do i = 1, nbas
         do j = 1, i
            key = (i - 1)*nbas + j
            if (pp_n(key) == 0) cycle
            nlive = nlive + 1
            d_i(nlive) = i; d_j(nlive) = j; d_k(nlive) = i; d_l(nlive) = j
            d_q(nlive) = pair_index(i, j)
         end do
      end do
      ! Longest first. The cost of a pair is the square of its primitive-pair
      ! count, which on a general contraction spans four orders of magnitude
      ! within one launch; leaving them in shell order puts the 20000-quartet
      ! silicon pairs at random points and every wave waits on its slowest
      ! thread. This is a partial sort by bucket, not a full one -- the order
      ! inside a bucket does not matter and n log n on 200000 pairs would show.
      call order_by_work(nlive, nbas, pp_n, d_i, d_j, d_k, d_l, d_q)
      if (nlive == 0) then
         deallocate (d_i, d_j, d_k, d_l, d_off, d_q)
         return
      end if
      ! Each pair takes only the room its own block needs. A uniform slot of
      ! the largest block in the molecule -- 36x36 for a d pair -- gives an
      ! (ss) pair 1296 doubles for the one it uses, and the zeroing alone
      ! then moves gigabytes: it took this routine from 10 s to 30 s on the
      ! silica slice before the sizes were made per-pair.
      nout = 0
      m = 0
      do c0 = 1, nlive, CHUNK
         c1 = min(c0 + CHUNK - 1, nlive)
         m = 0
         do n = c0, c1
            nb = ncart(sh_l(d_i(n)))*ncart(sh_l(d_j(n)))
            m = m + nb*nb
         end do
         nout = max(nout, m)
      end do
      allocate (out(nout), rscr(min(nlive, CHUNK), 2*NHERM_MAX), qc(min(nlive, CHUNK)))

      !$acc enter data copyin(d_i, d_j, d_k, d_l, sh_l) create(d_off, rscr, out, qc)
      do c0 = 1, nlive, CHUNK
         c1 = min(c0 + CHUNK - 1, nlive)
         nc = c1 - c0 + 1
         ! Offsets for this chunk, uniform blocks of nb*nb so the host does
         ! not have to walk the class list to find them again. Written at
         ! the GLOBAL pair index, because the kernel indexes q_off by the
         ! quartet's own index and not by its position in the chunk.
         m = 0
         do n = c0, c1
            d_off(n) = m
            nb = ncart(sh_l(d_i(n)))*ncart(sh_l(d_j(n)))
            m = m + nb*nb
         end do
         !$acc update device(d_off(c0:c1))
         ! The kernel accumulates into `out`, so it starts at zero -- on the
         ! device, where it lives, rather than on the host as it used to.
         do concurrent(na = 1:m)
            out(na) = 0.0_dp
         end do
         call eri_batch(c0, c1, nc, nlive, nbas, npp, m, &
                        d_i, d_j, d_k, d_l, d_off, sh_l, &
                        pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, rscr, out)
         ! The diagonal elements (ab|ab) are what Cauchy-Schwarz needs; the
         ! max over the block is a valid, looser bound and avoids the
         ! strided walk. Reduced here so `out` never leaves the device.
         do concurrent(n = 1:nc) local(na, nb, vmax)
            na = ncart(sh_l(d_i(c0 + n - 1)))*ncart(sh_l(d_j(c0 + n - 1)))
            vmax = 0.0_dp
            do nb = 1, na*na
               vmax = max(vmax, abs(out(d_off(c0 + n - 1) + nb)))
            end do
            qc(n) = sqrt(vmax)
         end do
         !$acc update self(qc(1:nc))
         do n = 1, nc
            q(d_q(c0 + n - 1)) = qc(n)
         end do
      end do
      !$acc exit data delete(d_i, d_j, d_k, d_l, d_off, sh_l, rscr, out, qc)

      deallocate (d_i, d_j, d_k, d_l, d_off, d_q, out, rscr, qc)
   end subroutine schwarz_bounds

   !
   ! Q(i,j) = sqrt( max_{ab} |(ij|ij)_{ab,ab}| ) for every i >= j, computed from
   ! the same contracted HGP primitive pairs (hp_*) the Fock kernels read, so no
   ! second primitive-pair set is built (unlike `schwarz_bounds` above, which is
   ! no longer called in production).
   !
   ! Scale: `build_pairs` (MMD) multiplied `common_fac_sp` into every primitive
   ! coefficient, but `basis_build` already folds it into `sh_c`, so the old
   ! route folded it twice and s/p bounds came out too small. That is not
   ! conservative: a too-small Cauchy-Schwarz bound can screen away a quartet
   ! that matters. Reading hp_c as the kernels do (build_structures passes
   ! cfac = 1) avoids it.
   !
   ! `hgp_quartet` is not class-specialised (its scratch is sized per class in
   ! `work`), so it is correct for every class including f, even though
   ! production routes f classes to Rys.
   !
   ! `hp_*` must already be device-resident (they stay resident for the Fock
   ! kernels); staging a copy here would be an implicit per-launch transfer.
   !
   subroutine schwarz_bounds_hgp(nbas, nhpp, sh_l, hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, q)
      integer,  intent(in)  :: nbas, nhpp
      integer,  intent(in)  :: sh_l(nbas)
      integer,  intent(in)  :: hp_off(nbas*nbas), hp_n(nbas*nbas)
      real(dp), intent(in)  :: hp_p(nhpp), hp_r(nhpp, 3), hp_ra(nhpp, 3), hp_rb(nhpp, 3), hp_c(nhpp)
      real(dp), intent(out) :: q(nbas*(nbas + 1)/2)

      integer, allocatable :: d_i(:), d_j(:), d_k(:), d_l(:), d_off(:), d_q(:), cst(:)
      integer(int64), allocatable :: w_off(:)
      real(dp), allocatable :: out(:), qc(:), work(:)
      integer :: npair, i, j, n, nout, nb, m, na, nlive, c0, c1, nc, key, lab, nchunk, ic
      integer(int64) :: nw, nwmax, need
      real(dp) :: vmax
      ! Same launch-width reasoning as schwarz_bounds above: one thread per
      ! shell pair, chunked. A chunk also ends before its VRR scratch (which
      ! hgp_batch takes from `work`, sized per quartet, rather than from each
      ! thread's stack) would pass WBUDGET doubles: 64 M, 512 MB. An (ss|ss)
      ! pair needs 3 doubles of it and an (ff|ff) one 21168, so a chunk of
      ! light pairs stays at CHUNK and one of f pairs is about 3000 wide.
      integer, parameter :: CHUNK = 65536
      integer(int64), parameter :: WBUDGET = 64_int64*1024*1024

      npair = nbas*(nbas + 1)/2
      q = 0.0_dp
      allocate (d_i(npair), d_j(npair), d_k(npair), d_l(npair), d_off(npair), d_q(npair))
      nlive = 0
      do i = 1, nbas
         do j = 1, i
            key = (i - 1)*nbas + j
            if (hp_n(key) == 0) cycle
            nlive = nlive + 1
            d_i(nlive) = i; d_j(nlive) = j; d_k(nlive) = i; d_l(nlive) = j
            d_q(nlive) = pair_index(i, j)
         end do
      end do
      ! Longest first, same reason as schwarz_bounds: a general contraction's
      ! primitive-pair count spans orders of magnitude across shell pairs.
      call order_by_work(nlive, nbas, hp_n, d_i, d_j, d_k, d_l, d_q)
      if (nlive == 0) then
         deallocate (d_i, d_j, d_k, d_l, d_off, d_q)
         return
      end if

      ! The chunks: starts in cst(1:nchunk), cst(nchunk + 1) = nlive + 1.
      allocate (cst(nlive + 1), w_off(nlive))
      nchunk = 0; nout = 0; nwmax = 0
      c0 = 1
      do while (c0 <= nlive)
         nchunk = nchunk + 1
         cst(nchunk) = c0
         m = 0; nw = 0
         c1 = c0 - 1
         do n = c0, min(c0 + CHUNK - 1, nlive)
            lab = sh_l(d_i(n)) + sh_l(d_j(n))
            need = hgp_scratch_size(lab, lab)
            if (n > c0 .and. nw + need > WBUDGET) exit
            nb = ncart(sh_l(d_i(n)))*ncart(sh_l(d_j(n)))
            m = m + nb*nb
            nw = nw + need
            c1 = n
         end do
         nout = max(nout, m); nwmax = max(nwmax, nw)
         c0 = c1 + 1
      end do
      cst(nchunk + 1) = nlive + 1
      allocate (out(nout), qc(min(nlive, CHUNK)), work(nwmax))

      ! d_i/d_j/d_k/d_l/d_off/w_off/sh_l/out/qc/work are this routine's own
      ! scratch -- entered and exited here, same as schwarz_bounds. hp_* are
      ! NOT entered here: see the header comment above.
      !$acc enter data copyin(d_i, d_j, d_k, d_l, sh_l) create(d_off, w_off, out, qc, work)
      do ic = 1, nchunk
         c0 = cst(ic); c1 = cst(ic + 1) - 1
         nc = c1 - c0 + 1
         m = 0; nw = 0
         do n = c0, c1
            d_off(n) = m
            w_off(n) = nw
            nb = ncart(sh_l(d_i(n)))*ncart(sh_l(d_j(n)))
            m = m + nb*nb
            lab = sh_l(d_i(n)) + sh_l(d_j(n))
            nw = nw + hgp_scratch_size(lab, lab)
         end do
         !$acc update device(d_off(c0:c1), w_off(c0:c1))
         do concurrent(na = 1:m)
            out(na) = 0.0_dp
         end do
         call hgp_batch(c0, c1, nlive, nbas, nhpp, m, &
                        d_i, d_j, d_k, d_l, d_off, sh_l, &
                        hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, out, &
                        nwmax, w_off, work)
         ! The diagonal elements (ab|ab) are what Cauchy-Schwarz needs; the
         ! max over the block is a valid, looser bound (same choice
         ! schwarz_bounds makes) and avoids the strided walk to the diagonal
         ! proper. Reduced here so `out` never leaves the device.
         do concurrent(n = 1:nc) local(na, nb, vmax)
            na = ncart(sh_l(d_i(c0 + n - 1)))*ncart(sh_l(d_j(c0 + n - 1)))
            vmax = 0.0_dp
            do nb = 1, na*na
               vmax = max(vmax, abs(out(d_off(c0 + n - 1) + nb)))
            end do
            qc(n) = sqrt(vmax)
         end do
         !$acc update self(qc(1:nc))
         do n = 1, nc
            q(d_q(c0 + n - 1)) = qc(n)
         end do
      end do
      !$acc exit data delete(d_i, d_j, d_k, d_l, d_off, w_off, sh_l, out, qc, work)

      deallocate (d_i, d_j, d_k, d_l, d_off, d_q, cst, w_off, out, qc, work)
   end subroutine schwarz_bounds_hgp

   !
   ! Reorder the live diagonal pairs so the expensive ones come first.
   !
   ! Bucketed by the power of two of the primitive-pair count, descending.
   ! That is enough: the point is to keep a launch from ending on one long
   ! thread, not to sort exactly.
   !
   subroutine order_by_work(nlive, nbas, pp_n, d_i, d_j, d_k, d_l, d_q)
      integer, intent(in)    :: nlive, nbas
      integer, intent(in)    :: pp_n(nbas*nbas)
      integer, intent(inout) :: d_i(*), d_j(*), d_k(*), d_l(*), d_q(*)

      integer, parameter :: NB = 32
      integer :: cnt(NB), start(NB), n, b, p, key
      integer, allocatable :: bkt(:), si(:), sj(:), sq(:)

      allocate (bkt(nlive), si(nlive), sj(nlive), sq(nlive))
      cnt = 0
      do n = 1, nlive
         key = (d_i(n) - 1)*nbas + d_j(n)
         p = pp_n(key)
         b = 1
         do while (p > 1 .and. b < NB)
            p = p/2; b = b + 1
         end do
         bkt(n) = NB + 1 - b          ! largest count first
         cnt(bkt(n)) = cnt(bkt(n)) + 1
      end do
      start(1) = 1
      do b = 2, NB
         start(b) = start(b - 1) + cnt(b - 1)
      end do
      do n = 1, nlive
         b = bkt(n)
         si(start(b)) = d_i(n); sj(start(b)) = d_j(n); sq(start(b)) = d_q(n)
         start(b) = start(b) + 1
      end do
      do n = 1, nlive
         d_i(n) = si(n); d_j(n) = sj(n); d_k(n) = si(n); d_l(n) = sj(n)
         d_q(n) = sq(n)
      end do
      deallocate (bkt, si, sj, sq)
   end subroutine order_by_work

   !
   ! Largest |D| over each shell pair's block.  Recomputed whenever D changes.
   !
   subroutine density_blockmax(nbas, nao, sh_l, ao_off, dmat, dmax)
      integer,  intent(in)  :: nbas, nao
      integer,  intent(in)  :: sh_l(nbas), ao_off(nbas)
      real(dp), intent(in)  :: dmat(nao, nao)
      real(dp), intent(out) :: dmax(nbas, nbas)
      integer :: i, j, a, b
      real(dp) :: m
      do i = 1, nbas
         do j = 1, nbas
            m = 0.0_dp
            do a = 0, ncart(sh_l(i)) - 1
               do b = 0, ncart(sh_l(j)) - 1
                  m = max(m, abs(dmat(ao_off(i) + a, ao_off(j) + b)))
               end do
            end do
            dmax(i, j) = m
         end do
      end do
   end subroutine density_blockmax

   !
   ! The compacted work list.
   !
   ! Outer loop over significant shell PAIRS rather than over all shells, so
   ! the enumeration is O(npair_significant^2) and never touches the quartets
   ! it is going to reject.  Enumerating all N^4 and testing would defeat the
   ! purpose.
   !
   ! `use_density` toggles screen 3.  With it off this is a pure Schwarz list,
   ! which is what the correctness tests want -- a density screen makes the
   ! answer depend on D, so it is validated separately by tightening thresh to
   ! zero and checking the list is complete.
   !
   subroutine build_worklist(nbas, sh_l, q, thresh, use_density, dmax, &
                             q_i, q_j, q_k, q_l, nq)
      integer,  intent(in)  :: nbas
      integer,  intent(in)  :: sh_l(nbas)
      real(dp), intent(in)  :: q(:)
      real(dp), intent(in)  :: thresh
      logical,  intent(in)  :: use_density
      real(dp), intent(in)  :: dmax(nbas, nbas)
      integer, allocatable, intent(out) :: q_i(:), q_j(:), q_k(:), q_l(:)
      integer,  intent(out) :: nq

      integer, allocatable :: pi(:), pj(:)
      real(dp), allocatable :: pq(:)
      integer :: npair, np, i, j, k, l, a, b, cap
      real(dp) :: qmax, bound, dfac

      ! --- significant pairs ---
      npair = nbas*(nbas + 1)/2
      allocate (pi(npair), pj(npair), pq(npair))
      qmax = 0.0_dp
      do i = 1, nbas
         do j = 1, i
            qmax = max(qmax, q(pair_index(i, j)))
         end do
      end do
      np = 0
      do i = 1, nbas
         do j = 1, i
            if (q(pair_index(i, j))*qmax > thresh) then
               np = np + 1
               pi(np) = i; pj(np) = j; pq(np) = q(pair_index(i, j))
            end if
         end do
      end do

      ! --- count, then fill.  Two passes so the arrays are exactly sized. ---
      cap = 0
      do a = 1, np
         do b = 1, a
            bound = pq(a)*pq(b)
            if (bound <= thresh) cycle
            if (use_density) then
               i = pi(a); j = pj(a); k = pi(b); l = pj(b)
               dfac = max(4.0_dp*dmax(i, j), 4.0_dp*dmax(k, l), &
                          dmax(i, k), dmax(i, l), dmax(j, k), dmax(j, l))
               if (bound*dfac <= thresh) cycle
            end if
            cap = cap + 1
         end do
      end do

      allocate (q_i(cap), q_j(cap), q_k(cap), q_l(cap))
      nq = 0
      do a = 1, np
         do b = 1, a
            bound = pq(a)*pq(b)
            if (bound <= thresh) cycle
            i = pi(a); j = pj(a); k = pi(b); l = pj(b)
            if (use_density) then
               dfac = max(4.0_dp*dmax(i, j), 4.0_dp*dmax(k, l), &
                          dmax(i, k), dmax(i, l), dmax(j, k), dmax(j, l))
               if (bound*dfac <= thresh) cycle
            end if
            nq = nq + 1
            q_i(nq) = i; q_j(nq) = j; q_k(nq) = k; q_l(nq) = l
         end do
      end do

      call sort_by_class(nbas, sh_l, nq, q_i, q_j, q_k, q_l)

      deallocate (pi, pj, pq)
   end subroutine build_worklist

   !
   ! Sort the work list by angular-momentum class.
   !
   ! Neighbouring work items then share (la,lb,lc,ld), so every thread in a warp
   ! runs identical trip counts and touches the same shaped data.  Without this
   ! a warp mixes (ss|ss) with (dd|dd) and runs at the pace of its worst lane:
   ! Nsight measured 12.4 of 32 threads active per warp, which is most of a 2.6x
   ! sitting on the floor.
   !
   ! Counting sort, O(nq + nclass), on the host once per geometry.
   !
   subroutine sort_by_class(nbas, sh_l, nq, q_i, q_j, q_k, q_l)
      integer, intent(in)    :: nbas, nq
      integer, intent(in)    :: sh_l(nbas)
      integer, intent(inout) :: q_i(nq), q_j(nq), q_k(nq), q_l(nq)

      integer :: lm, nclass, n, key, c
      integer, allocatable :: cnt(:), pos(:), keyv(:)
      integer, allocatable :: ti(:), tj(:), tk(:), tl(:)

      lm = maxval(sh_l) + 1
      nclass = lm*lm*lm*lm
      allocate (cnt(0:nclass - 1), pos(0:nclass - 1), keyv(nq))
      allocate (ti(nq), tj(nq), tk(nq), tl(nq))
      cnt = 0

      do n = 1, nq
         key = ((sh_l(q_i(n))*lm + sh_l(q_j(n)))*lm + sh_l(q_k(n)))*lm + sh_l(q_l(n))
         keyv(n) = key
         cnt(key) = cnt(key) + 1
      end do

      pos(0) = 0
      do c = 1, nclass - 1
         pos(c) = pos(c - 1) + cnt(c - 1)
      end do

      do n = 1, nq
         key = keyv(n)
         pos(key) = pos(key) + 1
         ti(pos(key)) = q_i(n); tj(pos(key)) = q_j(n)
         tk(pos(key)) = q_k(n); tl(pos(key)) = q_l(n)
      end do

      q_i = ti; q_j = tj; q_k = tk; q_l = tl
      deallocate (cnt, pos, keyv, ti, tj, tk, tl)
   end subroutine sort_by_class

end module trc_screen
