!
! Screening and the compacted, symmetry-reduced work list.
!
! Three reductions, applied in this order because each makes the next cheaper:
!
!   1. PERMUTATIONAL SYMMETRY, 8-fold.  Only canonical quartets are stored:
!      i >= j, k >= l, and pair(i,j) >= pair(k,l).  Not an optimisation -- an
!      8x factor is the difference between a calculation and a stunt.
!
!   2. SCHWARZ.  |(ij|kl)| <= Q(ij) Q(kl) with Q(ij) = sqrt(max |(ij|ij)|).
!      Rigorous (Cauchy-Schwarz on the Coulomb metric), so discarding below a
!      threshold bounds the error rather than hoping.
!
!   3. DENSITY.  What actually matters is a quartet's contribution to the Fock
!      matrix, not the integral's size, so the bound carries the largest
!      density block the quartet can touch.  This is the screen that makes late
!      SCF iterations cheap, and it has to be redone whenever D changes.
!
! The list is built ONCE on the host and only real quartets reach the device.
! Screening inside the kernel would not help: a GPU thread that exits early has
! already been scheduled, and its warp still waits for the busiest lane.  That
! is the whole point of section 0.2.
!
! WHY NOT libfint FOR THE BOUNDS
! ------------------------------
! It cannot be linked here -- libfint is gfortran-only (real128 is -1 under
! nvfortran, see NOTES.md), and in any case a library that needed a second
! integral package at runtime to screen would be a poor thing.  Q is computed
! with our own engine on the diagonal quartets, which is O(N^2) and reuses code
! that is already validated.  libfint's role stays what it has been: the oracle
! that proves the bounds are right, in check_screen.
!
module trc_screen
   use trc_boys, only: dp
   use trc_tables, only: LTOT, NHERM_MAX
   use trc_batch, only: eri_batch, ncart
   implicit none
   private

   public :: schwarz_bounds, pair_index, build_worklist, density_blockmax

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
      integer, parameter :: CHUNK = 8192

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
      if (nlive == 0) then
         deallocate (d_i, d_j, d_k, d_l, d_off, d_q)
         return
      end if
      nout = 0
      do n = 1, min(nlive, CHUNK)
         nout = max(nout, ncart(sh_l(d_i(n)))*ncart(sh_l(d_j(n))))
      end do
      nb = 0
      do n = 1, nlive
         nb = max(nb, ncart(sh_l(d_i(n)))*ncart(sh_l(d_j(n))))
      end do
      nout = min(nlive, CHUNK)*nb*nb
      allocate (out(nout), rscr(min(nlive, CHUNK), 2*NHERM_MAX), qc(min(nlive, CHUNK)))

      !$acc enter data copyin(d_i, d_j, d_k, d_l, sh_l) create(d_off, rscr, out, qc)
      do c0 = 1, nlive, CHUNK
         c1 = min(c0 + CHUNK - 1, nlive)
         nc = c1 - c0 + 1
         ! Offsets for this chunk, uniform blocks of nb*nb so the host does
         ! not have to walk the class list to find them again. Written at
         ! the GLOBAL pair index, because the kernel indexes q_off by the
         ! quartet's own index and not by its position in the chunk.
         do n = c0, c1
            d_off(n) = (n - c0)*nb*nb
         end do
         !$acc update device(d_off(c0:c1))
         ! The kernel accumulates into `out`, so it starts at zero -- on the
         ! device, where it lives, rather than on the host as it used to.
         do concurrent(m = 1:nc*nb*nb)
            out(m) = 0.0_dp
         end do
         call eri_batch(c0, c1, nc, nlive, nbas, npp, nc*nb*nb, &
                        d_i, d_j, d_k, d_l, d_off, sh_l, &
                        pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, rscr, out)
         ! The diagonal elements (ab|ab) are what Cauchy-Schwarz needs; the
         ! max over the block is a valid, looser bound and avoids the
         ! strided walk. Reduced here so `out` never leaves the device.
         do concurrent(n = 1:nc) local(na, m, vmax)
            na = ncart(sh_l(d_i(c0 + n - 1)))*ncart(sh_l(d_j(c0 + n - 1)))
            vmax = 0.0_dp
            do m = 1, na*na
               vmax = max(vmax, abs(out((n - 1)*nb*nb + m)))
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
