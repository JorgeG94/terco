!
! Pulay's DIIS on the commutator (cdiis_t), and ADIIS (Hu & Yang, JCP 132,
! 054109) on the convex hull of the Fock/density history while damping has
! ended but the commutator is still too large for CDIIS to trust (adiis_t):
! see trc_scf_driver's own header for the two owning most of the actual
! theory. This module owns their histories, pushes a new iterate into
! them, and extrapolates/combines a step, plus the plain subroutines both
! share: the DIIS linear solve, the level shift and the damped step while
! neither extrapolation is active.
!
! HOST-RESIDENT, PACKED HISTORY. Every matrix the histories hold is
! symmetric (F, D) or antisymmetric (the commutator), so each is kept as
! its lower triangle, n(n+1)/2 values, and on the HOST, not the device. A
! full-square device history was the largest single item of an SCF's
! device memory: 2 x ndiis x nao^2 doubles, 20 of the 50 nao^2 live at the
! peak of RHF/cc-pVTZ on gly20 (3298 functions), which is what ran
! gly100 (16,258 functions) out of a 141 GB H200. Packed on the host it is
! no device memory at all. A push packs the new pair on the device into
! the caller's nao x nao scratch (`staging`), brings it down, and copies
! it into its slot; the new B (or P) row is one host pass over the
! history; the extrapolated step is another, into `staging`, uploaded and
! unpacked into `fock`. GPU4PySCF's CDIIS packs the same way
! (pack_tril), and moves its history to the host when the device is short
! of memory. The cost is three triangles over the bus and two passes over
! the history per DIIS iteration: on a V100 box, RHF/cc-pVTZ gly20 (3298
! functions, the window full at ten) spends 0.12-0.18 s more per iteration
! than the device history did, against 24 s per Fock build; gly10 (1678)
! 0.04 s against 5.9 s. Iteration counts and energies did not move.
!
! What packing changes, numerically. F is stored as (F + F^T)/2 and the
! commutator as (E - E^T)/2, their exactly symmetric and antisymmetric
! parts: the full-square history carried each one's round-off-sized other
! part along. <E_i, E_j> is 2 x the strict lower triangle's dot product,
! and <D_i, F_j> the diagonal's plus 2 x the strict lower triangle's (the
! density's off-diagonal is stored doubled for exactly that, and is read
! by nothing but that dot product). The extrapolated F comes out exactly
! symmetric. The dot products are host BLAS instead of cuBLAS. All of it is
! round-off: iteration counts and energies are unchanged.
!
! RING BUFFER. Neither history is ever shifted. `head` is the physical
! slot the newest entry landed in; a push advances it one slot (wrapping
! at `ndiis`) and writes there, so the oldest entry's slot is simply the
! one about to be overwritten -- nothing else moves. `ring_phys(head,
! n_used, ndiis, l)` maps a LOGICAL position (1 = oldest, n_used = newest)
! to the physical slot holding it; every loop that used to assume slot k
! IS logical position k now goes through it instead.
!
! CACHED B / P. cdiis_t's `bcache(ndiis,ndiis)` and adiis_t's
! `pcache(ndiis,ndiis)` hold every dot product <e_i,e_j> / <D_i,F_j> this
! run has ever needed, indexed by PHYSICAL slot, kept across iterations.
! A push computes only the new row and column against every slot written
! so far (cdiis: one, symmetric; adiis: both, P is not), reusing every
! older pair untouched since nothing about an entry still in the window
! has changed. CDIIS dropping its oldest entry on a failed solve (see
! cdiis_extrapolate's own comment) is then just `n_used = n_used - 1` --
! the window shrinks, no data moves, and the freed slot's stale cache
! entries are silently overwritten the next time a push lands there.
!
! EXPLICIT-SHAPE DUMMIES. Every one of these takes `fock`/`dold`/`errv`/
! `staging` as plain arguments, never touching `st%anything` or a type's
! own OTHER component inside a `do concurrent` body -- see trc_scf_fock's
! own header for why that matters on the GPU build. The histories and
! caches are host-only, so nothing here needs `!$acc enter data`;
! `staging` must be device-resident (the SCF state's `w1`).
!
module trc_scf_diis
   use pic_types, only: dp, int64
   use pic_sorting, only: sort
   use trc_linalg, only: trc_linalg_t
   implicit none
   private

   public :: cdiis_t, adiis_t, scf_diis_level_shift, scf_diis_damp
   ! solve_diis stays private: nothing outside this module uses it.

   real(dp), parameter :: SYMMETRIC = 1.0_dp
   !! `pack_device`'s sign: store (A + A^T)/2
   real(dp), parameter :: ANTISYMMETRIC = -1.0_dp
   !! `pack_device`'s sign: store (A - A^T)/2
   integer, parameter :: HOST_BLOCK = 4096
   !! elements per block of the host passes over the history: 32 KB a
   !! vector, so a block of every slot stays in L2

   type :: cdiis_t
      real(dp), allocatable :: fhist(:, :, :), ehist(:, :, :)
      !! (npack, ndiis, nspin): the packed Fock matrices and commutators,
      !! host only. Slot-major inside a spin so one spin's window is a
      !! contiguous block for the host passes.
      real(dp), allocatable :: bmat(:, :), rhs(:)
      real(dp), allocatable :: bcache(:, :)   !! <e_i,e_j> by physical slot, see this module's own header
      real(dp), allocatable :: row(:), coeff(:)   !! a new B row and the extrapolation weights, by physical slot
      integer :: n_used = 0
      integer :: ndiis = 0
      integer :: head = 0   !! physical slot of the newest entry, 0 before the first push
      integer :: nfilled = 0   !! physical slots 1..nfilled have been written at least once
      integer :: npack = 0   !! nao(nao+1)/2
   contains
      procedure :: alloc => cdiis_alloc
      procedure :: push => cdiis_push
      procedure :: restart => cdiis_restart
      procedure :: extrapolate => cdiis_extrapolate
   end type cdiis_t

   type :: adiis_t
      real(dp), allocatable :: fhist(:, :, :), dhist(:, :, :)
      !! (npack, ndiis, nspin), as cdiis_t's; dhist holds the density with
      !! its off-diagonal doubled (this module's header).
      real(dp), allocatable :: cad(:)
      real(dp), allocatable :: pcache(:, :)   !! <D_i,F_j> by physical slot, see this module's own header
      real(dp), allocatable :: row(:), column(:), coeff(:)
      integer :: n_used = 0
      integer :: ndiis = 0
      integer :: head = 0
      integer :: nfilled = 0
      integer :: npack = 0
   contains
      procedure :: alloc => adiis_alloc
      procedure :: push => adiis_push
      procedure :: extrapolate => adiis_extrapolate
   end type adiis_t

contains

   subroutine cdiis_alloc(this, nao, nspin, ndiis)
      class(cdiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin, ndiis
      this%ndiis = ndiis
      this%n_used = 0
      this%head = 0
      this%nfilled = 0
      this%npack = npacked(nao)
      allocate (this%fhist(this%npack, ndiis, nspin), this%ehist(this%npack, ndiis, nspin))
      allocate (this%bmat(ndiis + 1, ndiis + 1), this%rhs(ndiis + 1))
      allocate (this%bcache(ndiis, ndiis), this%row(ndiis), this%coeff(ndiis))
      this%bcache = 0.0_dp
   end subroutine cdiis_alloc

   subroutine cdiis_restart(this)
      !! Forget the history; the next push starts a new one. Nothing is
      !! reallocated or zeroed: `push` rewrites every slot it then reads.
      class(cdiis_t), intent(inout) :: this
      this%n_used = 0
      this%head = 0
      this%nfilled = 0
   end subroutine cdiis_restart

   subroutine cdiis_push(this, nao, nspin, fock, errv, staging)
      !! The ring advances one slot (wrapping at ndiis, overwriting the
      !! oldest entry once the history is full) and (fock, errv) land
      !! there, packed; bcache's new row/column is the only part of it
      !! this push computes -- see this module's own header. `staging` is
      !! device-resident nao x nao scratch, clobbered.
      class(cdiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: fock(nao, nao, nspin), errv(nao, nao, nspin)
      real(dp), intent(inout) :: staging(nao, nao)
      integer :: l, p, s
      this%head = modulo(this%head, this%ndiis) + 1
      if (this%n_used < this%ndiis) this%n_used = this%n_used + 1
      this%nfilled = max(this%nfilled, this%head)
      do s = 1, nspin
         call pack_to_host(nao, this%npack, fock(:, :, s), SYMMETRIC, 1.0_dp, staging, this%fhist(:, this%head, s))
         call pack_to_host(nao, this%npack, errv(:, :, s), ANTISYMMETRIC, 1.0_dp, staging, this%ehist(:, this%head, s))
      end do
      ! <e_p, e_head> for every slot written so far; the diagonal of an
      ! antisymmetric triangle is zero, so the full dot product is twice it.
      call history_dots(this%npack, this%ndiis, nspin, this%nfilled, this%ehist, this%ehist, this%head, this%row)
      do l = 1, this%n_used
         p = ring_phys(this%head, this%n_used, this%ndiis, l)
         this%bcache(p, this%head) = 2.0_dp*this%row(p)
         this%bcache(this%head, p) = this%bcache(p, this%head)
      end do
   end subroutine cdiis_push

   subroutine cdiis_extrapolate(this, nao, nspin, fock, staging, ok)
      !!
      !! An ill-conditioned subspace loses its oldest entry, not the step.
      !! A failed solve followed by a plain Roothaan step is not neutral: on
      !! water/def2-SVP PBE0 the history holds errors of order 1 next to the
      !! newest at 1e-7, B's scaled pivots drop under solve_diis's 1e-14, and
      !! plain steps grow the commutator x1.28 per iteration along one
      !! direction, which keeps B singular: 65 iterations of runaway. Dropping
      !! the oldest entries until the system solves keeps the extrapolation;
      !! two distinct error vectors always solve, and that step cancels one
      !! growing mode.
      !!
      class(cdiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      real(dp), intent(inout) :: staging(nao, nao)
      logical, intent(out) :: ok
      integer :: m, i, j, pi, pj
      do
         m = this%n_used + 1
         this%bmat = -1.0_dp
         this%bmat(m, m) = 0.0_dp
         do j = 1, this%n_used
            pj = ring_phys(this%head, this%n_used, this%ndiis, j)
            do i = 1, j
               pi = ring_phys(this%head, this%n_used, this%ndiis, i)
               this%bmat(i, j) = this%bcache(pi, pj)
               this%bmat(j, i) = this%bmat(i, j)
            end do
         end do
         call solve_diis(m, this%bmat(1:m, 1:m), this%rhs(1:m), ok)
         if (ok .or. this%n_used <= 2) exit
         ! Drop the oldest: the window shrinks by one logical position: no
         ! data moves, since ring_phys reads head/n_used fresh every call.
         this%n_used = this%n_used - 1
      end do
      if (ok) then
         call weights_by_slot(this%n_used, this%rhs(1:this%n_used), this%head, this%ndiis, this%coeff)
         call combine_to_device(nao, nspin, this%ndiis, this%nfilled, this%fhist, this%coeff, &
                                staging, fock)
      end if
   end subroutine cdiis_extrapolate

   subroutine adiis_alloc(this, nao, nspin, ndiis)
      class(adiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin, ndiis
      this%ndiis = ndiis
      this%n_used = 0
      this%head = 0
      this%nfilled = 0
      this%npack = npacked(nao)
      allocate (this%fhist(this%npack, ndiis, nspin), this%dhist(this%npack, ndiis, nspin))
      allocate (this%cad(ndiis))
      allocate (this%pcache(ndiis, ndiis), this%row(ndiis), this%column(ndiis), this%coeff(ndiis))
      this%pcache = 0.0_dp
   end subroutine adiis_alloc

   subroutine adiis_push(this, nao, nspin, fock, dold, staging)
      !! The ring advances and stores (fock, dold) -- the Fock/density pair
      !! ADIIS's convex combination runs over, kept apart from cdiis_t's own
      !! history (see trc_scf_options_t%adiis's own comment on why) -- then
      !! pcache gets the new row AND column P(head,:)/P(:,head): P is not
      !! symmetric (D first, F second, always -- see this module's header),
      !! so, unlike bcache, both are needed.
      class(adiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: fock(nao, nao, nspin), dold(nao, nao, nspin)
      real(dp), intent(inout) :: staging(nao, nao)
      integer :: l, p, s
      this%head = modulo(this%head, this%ndiis) + 1
      if (this%n_used < this%ndiis) this%n_used = this%n_used + 1
      this%nfilled = max(this%nfilled, this%head)
      do s = 1, nspin
         call pack_to_host(nao, this%npack, fock(:, :, s), SYMMETRIC, 1.0_dp, staging, this%fhist(:, this%head, s))
         call pack_to_host(nao, this%npack, dold(:, :, s), SYMMETRIC, 2.0_dp, staging, this%dhist(:, this%head, s))
      end do
      ! column(p) = <D_p, F_head>, row(p) = <D_head, F_p>
      call history_dots(this%npack, this%ndiis, nspin, this%nfilled, this%dhist, this%fhist, this%head, this%column)
      call history_dots(this%npack, this%ndiis, nspin, this%nfilled, this%fhist, this%dhist, this%head, this%row)
      do l = 1, this%n_used
         p = ring_phys(this%head, this%n_used, this%ndiis, l)
         this%pcache(p, this%head) = this%column(p)
         this%pcache(this%head, p) = this%row(p)
      end do
   end subroutine adiis_push

   subroutine adiis_extrapolate(this, nao, nspin, fock, staging)
      !! One entry: the simplex is a point, c = 1 on it, a no-op. Skip the
      !! solve and leave `fock` as it is, which is also what `adiis_coeffs`
      !! would return from "start at the newest entry" after finding no
      !! improving step.
      class(adiis_t), intent(inout) :: this
      integer, intent(in) :: nao, nspin
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      real(dp), intent(inout) :: staging(nao, nao)
      if (this%n_used > 1) then
         call adiis_coeffs(this%n_used, this%pcache, this%head, this%ndiis, this%cad(1:this%n_used))
         call weights_by_slot(this%n_used, this%cad(1:this%n_used), this%head, this%ndiis, this%coeff)
         call combine_to_device(nao, nspin, this%ndiis, this%nfilled, this%fhist, this%coeff, &
                                staging, fock)
      end if
   end subroutine adiis_extrapolate

   pure function npacked(nao) result(npack)
      !! The lower triangle's length, diagonal included.
      integer, intent(in) :: nao
      integer :: npack
      npack = int(int(nao, int64)*int(nao + 1, int64)/2_int64)
   end function npacked

   subroutine pack_to_host(nao, npack, a, sign, offdiag, staging, slot)
      !! slot := the lower triangle of offdiag (A + sign A^T)/2 (diagonal
      !! unscaled), packed on the device into `staging` and brought down.
      integer, intent(in) :: nao, npack
      real(dp), intent(in) :: a(nao, nao)
      real(dp), intent(in) :: sign, offdiag
      real(dp), intent(inout) :: staging(npack)
      !! the caller's nao x nao scratch, seen as its first npack elements
      real(dp), intent(out) :: slot(npack)
      integer :: k
      call pack_device(nao, npack, a, sign, offdiag, staging)
      !$acc update self(staging(1:npack))
      ! A plain loop on the host: `staging` is the host copy here.
      do k = 1, npack
         slot(k) = staging(k)
      end do
   end subroutine pack_to_host

   subroutine pack_device(nao, npack, a, sign, offdiag, packed)
      !! packed(k) = w 0.5 (a(i,j) + sign a(j,i)) for i >= j, column by
      !! column (k from `packed_index`), w = offdiag below the diagonal and
      !! 1 on it. 0.5 (x + x) is x exactly, so a symmetric pack keeps the
      !! diagonal bit for bit and an antisymmetric one zeroes it.
      integer, intent(in) :: nao, npack
      real(dp), intent(in) :: a(nao, nao)
      real(dp), intent(in) :: sign, offdiag
      real(dp), intent(out) :: packed(npack)
      integer :: i, j
      do concurrent(j=1:nao, i=1:nao)
         if (i > j) then
            packed(packed_index(nao, i, j)) = offdiag*(0.5_dp*(a(i, j) + sign*a(j, i)))
         else if (i == j) then
            packed(packed_index(nao, i, j)) = 0.5_dp*(a(i, j) + sign*a(j, i))
         end if
      end do
   end subroutine pack_device

   subroutine unpack_device(nao, npack, packed, a)
      !! a := the symmetric matrix whose lower triangle is `packed`.
      integer, intent(in) :: nao, npack
      real(dp), intent(in) :: packed(npack)
      real(dp), intent(out) :: a(nao, nao)
      integer :: i, j
      do concurrent(j=1:nao, i=1:nao)
         a(i, j) = packed(packed_index(nao, max(i, j), min(i, j)))
      end do
   end subroutine unpack_device

   pure function packed_index(nao, i, j) result(k)
      !! Position of (i, j), i >= j, in the column-major lower triangle:
      !! the j-1 columns before it hold (j-1)(2 nao - j + 2)/2 entries (an
      !! even product), then i - j + 1 into column j. int64: 2 nao^2
      !! overflows a default integer past 32768 functions.
      !$acc routine seq
      integer, intent(in) :: nao, i, j
      integer(int64) :: k
      k = int(j - 1, int64)*int(2*nao - j + 2, int64)/2_int64 + int(i - j + 1, int64)
   end function packed_index

   subroutine history_dots(npack, ndiis, nspin, nfilled, hist, other, slot, dots)
      !! dots(p) = sum over spins of <hist(:, p, s), other(:, slot, s)>, p =
      !! 1..nfilled: every slot written so far against one, in one pass
      !! over the history, block by block so the block of `other` stays in
      !! cache. `other` may be `hist` itself (both are read only).
      !!
      !! Plain host loops, not the host BLAS: the GPU build's host BLAS is
      !! whatever the toolchain links (MKL's threaded GEMV, on the V100
      !! box), and its OpenMP runtime next to NVHPC's crashed at the first
      !! call. Serial, the pass streams the history at memory bandwidth.
      integer, intent(in) :: npack, ndiis, nspin, nfilled, slot
      real(dp), intent(in) :: hist(npack, ndiis, nspin), other(npack, ndiis, nspin)
      real(dp), intent(out) :: dots(ndiis)
      integer :: s, p, k, k0, k1
      real(dp) :: acc
      dots = 0.0_dp
      do s = 1, nspin
         do k0 = 1, npack, HOST_BLOCK
            k1 = min(npack, k0 + HOST_BLOCK - 1)
            do p = 1, nfilled
               acc = 0.0_dp
               do k = k0, k1
                  acc = acc + hist(k, p, s)*other(k, slot, s)
               end do
               dots(p) = dots(p) + acc
            end do
         end do
      end do
   end subroutine history_dots

   subroutine weights_by_slot(m, weights, head, ndiis, coeff)
      !! coeff(p) = the weight of whichever logical position slot p holds,
      !! zero for a slot outside the window -- what combine_to_device's
      !! pass over every written slot reads.
      integer, intent(in) :: m, head, ndiis
      real(dp), intent(in) :: weights(m)
      real(dp), intent(out) :: coeff(ndiis)
      integer :: k
      coeff = 0.0_dp
      do k = 1, m
         coeff(ring_phys(head, m, ndiis, k)) = weights(k)
      end do
   end subroutine weights_by_slot

   subroutine combine_to_device(nao, nspin, ndiis, nfilled, fhist, coeff, staging, fock)
      !! fock(:, :, s) := sum_p coeff(p) fhist(:, p, s), per spin: the
      !! extrapolated/combined step both CDIIS (coeff = the solved weights)
      !! and ADIIS (coeff = the simplex coefficients) take once they have
      !! more than one history slot. Formed on the host in `staging`, uploaded
      !! and unpacked on the device. Slots outside the window carry weight
      !! zero and are skipped.
      integer, intent(in) :: nao, nspin, ndiis, nfilled
      real(dp), intent(in) :: fhist(npacked(nao), ndiis, nspin), coeff(ndiis)
      real(dp), intent(inout) :: staging(npacked(nao))
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      integer :: s, npack, p, k, k0, k1
      npack = npacked(nao)
      do s = 1, nspin
         ! Host loops, blocked, for the reason history_dots gives.
         do k0 = 1, npack, HOST_BLOCK
            k1 = min(npack, k0 + HOST_BLOCK - 1)
            do k = k0, k1
               staging(k) = 0.0_dp
            end do
            do p = 1, nfilled
               if (coeff(p) == 0.0_dp) cycle
               do k = k0, k1
                  staging(k) = staging(k) + coeff(p)*fhist(k, p, s)
               end do
            end do
         end do
         !$acc update device(staging(1:npack))
         call unpack_device(nao, npack, staging, fock(:, :, s))
      end do
   end subroutine combine_to_device

   subroutine scf_diis_level_shift(lp, nao, nspin, mu, occ, smat, dmat, w1, w2, fock)
      !! F += mu (S - S D S / occ) on the virtual space while DIIS is off,
      !! per spin: see trc_scf_options_t%scf_diis_level_shift's own comment.
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: mu, occ
      real(dp), intent(in) :: smat(nao, nao), dmat(nao, nao, nspin)
      real(dp), intent(inout) :: w1(nao, nao), w2(nao, nao)
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      integer :: i, j, s
      do s = 1, nspin
         call lp%gemm_mat(smat, dmat(:, :, s), w1)
         call lp%gemm_mat(w1, smat, w2)
         do concurrent(i=1:nao, j=1:nao)
            fock(i, j, s) = fock(i, j, s) + mu*(smat(i, j) - w2(i, j)/occ)
         end do
      end do
   end subroutine scf_diis_level_shift

   subroutine scf_diis_damp(nao, nspin, alpha, dold, dmat)
      !! The damped step while DIIS/A-DIIS are off: dmat := alpha*dmat +
      !! (1 - alpha)*dold, the weight `alpha` given to the new density.
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: dold(nao, nao, nspin)
      real(dp), intent(inout) :: dmat(nao, nao, nspin)
      integer :: i, j, s
      do concurrent(i=1:nao, j=1:nao, s=1:nspin)
         dmat(i, j, s) = alpha*dmat(i, j, s) + (1.0_dp - alpha)*dold(i, j, s)
      end do
   end subroutine scf_diis_damp

   pure function ring_phys(head, n_used, ndiis, l) result(p)
      !! The physical slot holding logical position `l` (1 = oldest, n_used
      !! = newest) of a ring buffer of `ndiis` slots whose newest entry is
      !! at physical slot `head`. `modulo`, not `mod`: head - n_used + l - 1
      !! is negative whenever l is not the newest position and the ring has
      !! wrapped, and modulo is the one that stays in 0..ndiis-1 there.
      integer, intent(in) :: head, n_used, ndiis, l
      integer :: p
      p = modulo(head - n_used + l - 1, ndiis) + 1
   end function ring_phys

   !
   ! ADIIS coefficients over the `m` logical positions cached in `pcache`
   ! (physical slots, see this module's own header). Logical position `m`
   ! is the newest, so Hu & Yang's "entry n" is the last index. Minimises
   ! over the simplex {c >= 0, sum c = 1}
   !
   !     f(c) = 2 sum_i c_i <D_i - D_m, F_m>
   !            + sum_ij c_i c_j <D_i - D_m, F_j - F_m>
   !
   ! read off pcache(i,j) = <D_i,F_j> as pcache(i,j) - pcache(i,m) -
   ! pcache(m,j) + pcache(m,m) -- the same values combined in the same
   ! expression as the un-cached form -- by projected gradient with
   ! backtracking (a few hundred host flops once the dot products are in
   ! hand). Starting with all weight on the newest entry means no progress
   ! degrades to the plain Roothaan step.
   !
   subroutine adiis_coeffs(m, pcache, head, ndiis, c)
      integer, intent(in) :: m, head, ndiis
      real(dp), intent(in) :: pcache(ndiis, ndiis)
      real(dp), intent(out) :: c(m)
      real(dp) :: f, ft, step
      real(dp) :: dvec(m), amat(m, m), g(m), trial(m)
      integer :: ii, jj, itr, ls, pii, pjj, pm
      logical :: improved
      pm = head   !! logical position m (the newest) is always at head
      do ii = 1, m
         pii = ring_phys(head, m, ndiis, ii)
         dvec(ii) = pcache(pii, pm) - pcache(pm, pm)
         do jj = 1, m
            pjj = ring_phys(head, m, ndiis, jj)
            amat(ii, jj) = pcache(pii, pjj) - pcache(pii, pm) - pcache(pm, pjj) + pcache(pm, pm)
         end do
      end do
      c = 0.0_dp
      c(m) = 1.0_dp
      f = adiis_value(m, dvec, amat, c)
      step = 1.0_dp
      do itr = 1, 200
         do ii = 1, m
            g(ii) = 2.0_dp*dvec(ii)
            do jj = 1, m
               g(ii) = g(ii) + (amat(ii, jj) + amat(jj, ii))*c(jj)
            end do
         end do
         improved = .false.
         do ls = 1, 30
            trial = c - step*g
            call project_simplex(m, trial)
            ft = adiis_value(m, dvec, amat, trial)
            if (ft < f - 1.0e-14_dp) then
               c = trial
               f = ft
               improved = .true.
               exit
            end if
            step = step*0.5_dp
         end do
         if (.not. improved) exit
         step = step*1.5_dp
      end do
   end subroutine adiis_coeffs

   pure function adiis_value(m, dvec, amat, c) result(f)
      integer, intent(in) :: m
      real(dp), intent(in) :: dvec(m), amat(m, m), c(m)
      real(dp) :: f
      integer :: ii, jj
      f = 0.0_dp
      do ii = 1, m
         f = f + 2.0_dp*c(ii)*dvec(ii)
         do jj = 1, m
            f = f + c(ii)*c(jj)*amat(ii, jj)
         end do
      end do
   end function adiis_value

   ! Projection onto the simplex {c : sum c = 1, c >= 0}: the standard
   ! sort-and-threshold algorithm, which is what enforces ADIIS's
   ! non-negativity constraint, the thing that keeps it from wandering
   ! off far from convergence.
   pure subroutine project_simplex(m, c)
      integer, intent(in) :: m
      real(dp), intent(inout) :: c(m)
      real(dp) :: css, theta, t
      real(dp) :: u(m)
      integer :: ii, rho
      u = c
      call sort(u, reverse=.true.)
      css = 0.0_dp
      theta = 0.0_dp
      rho = 0
      do ii = 1, m
         css = css + u(ii)
         t = (css - 1.0_dp)/real(ii, dp)
         if (u(ii) - t > 0.0_dp) then
            rho = ii
            theta = t
         end if
      end do
      if (rho == 0) then
         c = 1.0_dp/real(m, dp)
         return
      end if
      c = max(c - theta, 0.0_dp)
   end subroutine project_simplex

   !
   ! The DIIS system B c = (0, ..., 0, -1), pivoted elimination on the host,
   ! with the error block scaled to O(1) against the -1 border first. The
   ! weights are invariant under that scaling; only the Lagrange multiplier
   ! moves, and nothing reads it. `ok` is false on a pivot below 1e-14, in
   ! which case the caller keeps the unextrapolated Fock matrix.
   !
   subroutine solve_diis(m, b, c, ok)
      integer, intent(in) :: m
      real(dp), intent(in) :: b(m, m)
      real(dp), intent(out) :: c(m)
      logical, intent(out) :: ok
      real(dp), allocatable :: a(:, :)
      real(dp) :: scale, pivot, factor
      integer :: i, j, p
      allocate (a(m, m + 1))
      a(:, 1:m) = b
      a(:, m + 1) = 0.0_dp
      a(m, m + 1) = -1.0_dp
      scale = maxval(abs(b(1:m - 1, 1:m - 1)))
      if (scale > 0.0_dp) a(1:m - 1, 1:m - 1) = a(1:m - 1, 1:m - 1)/scale
      ok = .true.
      c = 0.0_dp
      do i = 1, m
         p = i
         do j = i + 1, m
            if (abs(a(j, i)) > abs(a(p, i))) p = j
         end do
         if (p /= i) a([i, p], :) = a([p, i], :)
         pivot = a(i, i)
         if (abs(pivot) < 1.0e-14_dp) then
            ok = .false.
            return
         end if
         do j = i + 1, m
            factor = a(j, i)/pivot
            a(j, :) = a(j, :) - factor*a(i, :)
         end do
      end do
      do i = m, 1, -1
         c(i) = (a(i, m + 1) - sum(a(i, i + 1:m)*c(i + 1:m)))/a(i, i)
      end do
   end subroutine solve_diis

end module trc_scf_diis
