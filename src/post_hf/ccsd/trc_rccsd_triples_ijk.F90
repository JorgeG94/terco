!
! The perturbative triples (T) of closed-shell RI-CCSD(T) on the DEVICE, driven by OCCUPIED triples i >= j >= k: the
! same quantity as `rccsd_triples_reference` (trc_rccsd_triples) and `rccsd_triples_dev_run` (trc_rccsd_triples_dev,
! the virtual-triple-driven routine), restructured so that every library call is sized by the VIRTUAL space.
!
! WHY. The virtual-triple routine stores W(i,j,k) blocks of no^3 per triple a >= b >= c, so its gemms have m = no
! and k = nv or no: millions of skinny products. Here one occupied triple holds three arrays over ALL virtuals (nv^3
! each), filled by tall gemms of (nv^2 x (nv + no))((nv + no) x n), n = nv or nv times a run of triples, with only
! no (no + 1) (no + 2) / 6 triples. Measured on one GH200 (H10O5, cc-pVTZ / cc-pVTZ-RIFIT, frozen core, no 20,
! nv 265, naux 705, default budget, E(T) the same to all 12 printed digits), seconds: the virtual-triple routine
! 11.90 (gemms 6.82 + 2.14, energy 2.92); this one 3.42 (wide gemms 1.73 at 35 TFLOP/s, per-array gemms 1.04 at
! 27 TFLOP/s, energy 0.59, repacks and slice gathers 0.02).
!
! THE ALGEBRA (an exact rewrite of the reference, no approximation). For an ORDERED occupied triple (o1,o2,o3) and
! ordered virtuals (x,y,z), pairs (o1,x), (o2,y), (o3,z):
!
!     q[o1,o2,o3](x,y,z) = sum_f (o1 x|f y) t2(o3,o2,z,f) - sum_m (o1 x|o2 m) t2(m,o3,y,z)
!         (o1 x|f y) = ovvv(o1,x,f,y),  (o1 x|o2 m) = ooov(o2,m,o1,x)
!
! the reference's get_w(a,b,c)(i,j,k). For an occupied triple (i,j,k) and Q1 = q[i,j,k], Q2 = q[i,k,j],
! Q3 = q[j,i,k], Q4 = q[j,k,i], Q5 = q[k,i,j], Q6 = q[k,j,i], for any ordered virtual triple (p,q,r)
!
!     W(p,q,r) = Q1(p,q,r) + Q2(p,r,q) + Q3(q,p,r) + Q4(q,r,p) + Q5(r,p,q) + Q6(r,q,p)
!     V(p,q,r) = ovov(i,p,j,q) t1(k,r) + ovov(i,p,k,r) t1(j,q) + ovov(j,q,k,r) t1(i,p),   X = W + V
!     D(p,q,r) = e_i + e_j + e_k - e_p - e_q - e_r
!
! (the array of the occupied arrangement sigma(i,j,k) is read at the same arrangement sigma of the virtual arguments).
! Then, with R(P) = 3 W(P) + S_same - 2 S_opp, S_even = W(abc) + W(bca) + W(cab), S_odd = W(acb) + W(bac) + W(cba)
! and "same" the parity class of the ordering P,
!
!     E_ijk = sum_{a>=b>=c} 1/(s_v D) sum_{P in orderings of (a,b,c)} X(P) R(P),   s_v = 1, 2, 6 (distinct, two, three equal)
!     E(T)  = (1/3) sum_{i>=j>=k} m_o E_ijk,                                      m_o = 6, 3, 1 (distinct, two, three equal)
!
! The 1/3 is the reference's factor 2 times the 1/6 of symmetrising X over the six pair permutations. V in three
! terms uses (ia|jb) = (jb|ia), which the fitted ovov has; the products below read t2 in the index order of q and use
! no symmetry of it, and the integral slices are built from the factors with b_vv(a,c;P) = b_vv(c,a;P) (SLICE CACHE),
! which the fitted factor has to round-off. The a = b = c terms are kept. An i = j = k triple contributes zero (W is then
! symmetric in (a,b,c) and R vanishes): it is skipped, and `diagonal_triples` computes it (round-off zeros, tested).
!
! TERMS. Only the sum W matters, so each array may carry the particle part of one arrangement and the hole part of
! another. For rho = (r1,r2,r3) the term T_rho pairs the particle part of rho with the hole part of (r3,r2,r1) (a
! bijection of the six, so every hole part appears once), read in rho's argument order:
!
!   T_rho(u,v,w) = sum_f ovvv(r1,u,f,v) t2(r3,r2,w,f) - sum_m t2(m,r1,v,u) ooov(r2,m,r3,w)
!                = sum_{K = 1..nv+no} A'_{r1}((u,v), K) B'_{r2 r3}(K, w)
!   A'_o = [ovvv(o,x,f,y) | t2(m,o,y,x)]             (nv^2 x (nv + no)): depends on r1 only
!   B'_{o2 o3} = [t2(o3,o2,w,f) ; -ooov(o2,m,o3,w)]  ((nv + no) x nv): depends on (r2,r3) only
!
! The particle and hole parts share one gemm with nv + no contraction columns; the hole alone, (nv x no)(no x nv^2),
! was bound by writing Q twice (3.70 s of separate particle and hole gemms became 2.78 s fused). All B' are built
! once per call into tvb (nv + no, nv, no, no), and tvb2(:,:,o3,o2) = tvb(:,:,o2,o3) so that a run over o3 is also
! contiguous (242 MB each here). Every gemm operand is read in place.
!
! THREE FRAMES. With A'_o, T_rho comes out of F((u,v),w) = A'_{r1} B'_{r2 r3} in the frame rho (storage positions
! 1,2,3 carry the virtuals of r1,r2,r3). With the slot held swapped, A2'_o((v,u),K) = A'_o((u,v),K), the same tall
! gemm writes it in the frame (r2,r1,r3). Two terms that end in the same frame share an array (beta = 0, then 1):
!
!   U1 frame (i,j,k) = T(i,j,k) [A'_i, B'_jk] + T(j,i,k) [A2'_j, B'_ik]
!   U2 frame (i,k,j) = T(i,k,j) [A'_i, B'_kj] + T(k,i,j) [A2'_k, B'_ij]
!   U3 frame (j,k,i) = T(j,k,i) [A'_j, B'_ki] + T(k,j,i) [A2'_k, B'_ji]
!   W(p,q,r) = U1(p,q,r) + U2(p,r,q) + U3(q,r,p),   (p,q,r) the virtuals of (i,j,k);  i = j: U3 = U2;  j = k: U2 = U1
!
! 18 Q loads per point. Tall rather than cyclic: writing a term in another frame by the transposed product
! F(w,(u,v)) = B'^T A'^T (four arrays, 24 loads) has the short-m shape 265 x 70225, which measured 21 TFLOP/s on this
! card with or without pre-transposed operands, against 28.6 for the tall 70225 x 265 and 33-35 for tall and wide
! (70225 x 1060 or more). The four-array version took 3.08 s of gemms for the same flops and 4.11 s in all.
!
! SLICE CACHE. No (ia|bc) block exists (trc_rccsd_ints): the integral columns of A'_o and A2'_o are built from the
! factors into S slots keyed by (orbital, layout). The rows b_ov(o,x;P) are strided (o fastest), so one kernel gathers
! them into bo(x,P) (nv naux); then A'_o((x,(y,f))) = sum_P bo(x,P) b_vv((y,f),P) = (o x|y f) = (o x|f y) is ONE gemm
! (nv, nv^2, naux, 2 nv^3 naux flops) with b_vv read in place, straight into the slot. A2'_o(x,y,f) = A'_o(y,x,f): when
! the other layout of the same orbital is cached it is a permuting copy from that slot (2 nv^3 doubles moved instead of
! the gemm), else the gemm and an in-place exchange of x and y. The run needs A'_j and A2'_j together, so most second
! layouts are copies. A kernel appends the t2 columns. S = clamp(slice_bytes / (8 nv^2 (nv + no)), 4, 2 no); with
! `slice_bytes` = 0 the device build sizes the whole plan to 80% of the free device memory, and at least 2 GiB short
! of all of it (queried before the plan is mapped, so what the CCSD left mapped is counted); the slices get what the
! rest of the plan leaves of that. The host build gives them 16 GiB.
!
! Measured on GH200 with every slice cached (seconds of slice work, "repacks and slice builds", which also holds the
! tvb repacks): no 20, nv 265, naux 705, S = 40 = 2 no: 0.029 for 40 builds, against 0.022 gathering them from the
! dense block (the 20 gemms are 0.52 TFLOP); (H2O)6 cc-pVTZ, no 24, nv 318, naux 846: 0.187 for 272 builds with
! 8 GiB (S = 31 of 48), against 0.246 gathering, and 0.042 for 48 builds with S = 48 = 2 no, where the whole (T) went
! from 9.34 s to 9.06 s: hence 16 GiB (13.3 GB of slots there, 5.1 GB more than at 8 GiB; the dense block it
! replaces was 6.2 GB). Three role slots in the first version took 1520 gathers and 0.69 s.
!
! BLOCKED ORDER. A triple reads A'_i, A'_j, A2'_j and A2'_k. With S < 2 no the canonical order cycles A2'_k, k = 1..j,
! for every (i,j), and least-recently-used eviction then drops exactly the slice needed next. So the occupied orbitals
! are cut into blocks of B = S / 4 (B = no, one block and the canonical order, when S >= 2 no) and the triples are
! visited block triple by block triple, IB >= JB >= KB with KB fastest, inside one i, j, k-fastest with i in IB, j in
! JB, k in KB (and i >= j >= k). A block triple reads at most 4 B <= S slices (A' of IB and JB, A2' of JB and KB):
! `pin_block` stamps the cached ones, then builds the rest (A' before A2' of one orbital, so the second is a copy), so
! nothing is evicted inside a block triple and the eviction between them only drops slices the next one does not read.
! Builds become about B per IB, per (IB,JB) and per block triple. Measured on GH200, (H2O)9 cc-pVTZ / cc-pVTZ-RIFIT
! (no 36, nv 477, naux 1269), E(T) the same to all 12 printed digits, seconds of the (T) less the mapping in `init`
! (4 to 47 s for the same plan from run to run): S = 18 (16 GiB) in the canonical order 6208 builds, 53.7 s of slice
! work, 233 (the sum of the stages); the same S blocked (B = 4) 729 builds (559 GEMM, 170 copies), 5.2 s, 186-189;
! S sized from free memory, 72 = 2 no (74.9 GB of scratch, device peak 90.1 of 97.9 GB), 72 builds, 0.31 s, 184-185.
!
! WIDE GEMMS. A tile is a run of consecutive k at one (i,j), cut where (i,j) changes, at the end of the block KB or at
! capacity, stored as q(nv^3, triple, frame) so one frame's arrays are consecutive in k. Over a run k = k0..kb four of
! the six products are one gemm each with n = nv (kb - k0 + 1): U1 [A'_i, tvb2(:,:,k,j)] and [A2'_j, tvb2(:,:,k,i)],
! U2 [A'_i, tvb(:,:,k,j)], U3 [A'_j, tvb(:,:,k,i)]; the swapped parts of U2 and U3 have the slot A2'_k and stay one
! gemm per triple. At i = j U3 is not computed (it is U2); at the j = k end of a run the U2 products stop at k = j - 1
! (U2 is U1 there); a skipped i = j = k triple is left out of the run. No array is computed twice.
! Measured (energy and gemms per tile size, seconds): tiles of 1 triple 3.91 in all (gemms 2.02 + 1.05), 2 triples
! 3.61, 4 triples 3.42 (wide 1.73 at 35 TFLOP/s), 8 triples 3.40, 18 triples 3.40: the default budget, 2 GiB, holds
! 4 triples at nv = 265 (447 MB each).
!
! ENERGY. Per tile, V's blocks are gathered compactly (vm = the three nv x nv matrices ovov(i,:,j,:), ovov(i,:,k,:),
! ovov(j,:,k,:) and vt = the three t1 rows, 1.7 MB per triple) so that the item reads them from cache, with what
! depends on (b,c) only read outside the loop over a. Then a `do concurrent` over the flattened (c, b, triple) index;
! items with c > b return 0, the others loop a = b..nv, form the six W, X and R and return the row sum. Two more
! `do concurrent` add the rows of each triple in index order (first over c for each (b, triple), then over b). Which
! index the item loops over was measured with six arrays (energy kernels, s): over a 2.94, over c 3.42, over b 4.76;
! threads sweeping 2-D (c,b) blocks were slower (8x4 2.82, 16x2 2.82, 4x8 3.04 against 2.33). Energy kernels went
! 2.25 (six arrays) -> 2.05 (compact V) -> 0.97 (four arrays) -> 0.59 (three, tiles of 4). The per-triple energies
! come to the host, are stored at the triple's canonical position (i, then j, then k fastest) whatever block triple
! computed them, and are added with m_o in that order; then times 1/3. No atomics, no dependence on thread
! scheduling: the tile size and the block size only change how the work is cut, and E(T) and the per-triple
! energies are the same for any slice cache (tested): bitwise on the V100 and the GH200, to the last bit of a triple
! (3e-20) with a host BLAS, whose gemm may round a column differently as the tile's width changes.
!
! WORK SHARE. `parts` runs (one per MPI rank, one device each) can split the triples: the unit is the occupied PAIR
! (i,j), i >= j, whose cost is its number of computed triples (j, less the i = j = k one unless `diagonal_triples`).
! The pairs are visited by decreasing cost, ties in canonical order, and each goes to the part with the fewest triples
! so far (ties to the lowest part): every rank recomputes the same table (`rccsd_triples_ijk_pair_owners`), and the
! largest part exceeds the mean by at most (1 - 1/parts) times the costliest pair (no - 1, or no with the diagonal;
! equality at no = 1, diagonal on, 2 parts). Pairs, not triples, because a tile is a run of k at one (i,j) and so never
! spans two owners. The deal does NOT shrink the slice work: pairs of equal j go round-robin, so every part reads
! nearly every slice and builds about the slices of the one-part run (model, no 36, 18 slots, blocks of 4: 729 builds
! for one part, 729 for each part of 2 and of 4). That is accepted: slice building is a small share of (T) (0.3 s of
! 184 s at that size, 0.03 s of 3.4 s with every slice cached); smaller cases build fewer (h2o2 test, 4 slots: 94 builds
! for one part, at most 50 for a part of 2); a block-aware deal is the lever if it ever matters. A
! block triple with no owned triple is skipped before its slices are built. A part computes only its triples, leaves 0
! at the others' canonical positions of `triple_energies` and sums its own in canonical order, so the per-triple
! vectors of all parts add to the one-part run's (tested, to round-off of the gemm width) and their E(T) add to its
! E(T).
!
! TILES. Q is 3 nv^3 T doubles (24 nv^3 bytes per triple); T = budget / (24 nv^3), `tile_triples` forces it
! (tests). `budget_bytes` = 0 means 2 GiB, raised to one triple if that is more (a default must not refuse a large
! virtual space); an explicit budget below one triple is a refusal. The three slots are kept for a degenerate triple
! too, so the plan does not depend on which triples a tile holds.
!
! DEVICE MEMORY, doubles: tvb and tvb2 (nv + no) nv no^2 each, the slots S nv^2 (nv + no), Q 3 nv^3 T, vm 3 nv^2 T,
! vt 3 nv T, row sums nv^2 T, column sums nv T, bo nv naux, plus T energies and 7 T integers. At no = 20, nv = 265,
! naux = 705, S = 40, T = 4: 484 MB + 6401 MB + 1787 MB + 7 MB + 3 MB + 1.5 MB = 8685 MB (logged as the device scratch).
!
! RESIDENCY, as trc_rccsd_triples_dev: `ints` and `amps` must be current on the device; the scratch maps its own
! arrays (`enter data create` in `init`, `exit data delete` in `release`) and records `ints%stamp`.
!
! INDEX RANGES. Q, the slots and every flat kernel index are 64-bit (3 nv^3 T passes 2^31 already at nv = 265). The
! gemm dimensions and leading dimensions are 32-bit: nv^2 and no^2 are checked in `init`.
#include "inc/trc_dc_index.inc"

module trc_rccsd_triples_ijk
   use pic_types, only: dp, int32, int64
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_linalg, only: trc_linalg_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_triples_dev, only: trc_rccsd_triples_dev_options_t
   implicit none
   private

   public :: trc_rccsd_triples_ijk_scratch_t
   public :: rccsd_triples_ijk_run, rccsd_triples_ijk_run_scratch
   public :: rccsd_triples_ijk_pair_owners, rccsd_triples_ijk_pair_cost, rccsd_triples_ijk_pair_index
   public :: rccsd_triples_ijk_triple_bytes

   integer(int64), parameter :: DEFAULT_BUDGET = 2147483648_int64  !! bytes of the Q tile when `budget_bytes` is 0
   integer(int64), parameter :: DEFAULT_SLICE_BYTES = 17179869184_int64
   !! bytes of cached slices when `slice_bytes` is 0 on the host, or when the free device memory is unknown (header)
   real(dp), parameter :: SLICE_FREE_FRACTION = 0.8_dp
   !! share of the free device memory the whole plan (slices and the rest) may take when the cache is sized from it:
   !! the rest stays for cuBLAS workspaces, later allocations and other processes on the device
   integer(int64), parameter :: SLICE_RESERVE = 2147483648_int64
   !! the least device memory left free when the slice cache is sized from free memory (the floor of the 20% margin
   !! on small devices)

   type :: trc_rccsd_triples_ijk_scratch_t
      !! The device arrays of one occupied-driven (T) and the plan they were sized for, see the module header
      type(error_t) :: status   !! set by `init`: a refusal leaves its message here and the type empty
      integer :: no = 0
      integer :: nv = 0
      integer :: stamp = 0               !! `ints%stamp` of the integrals `init` planned for
      integer(int64) :: ncap = 0         !! occupied triples per tile
      integer(int64) :: ntotal = 0       !! no (no + 1) (no + 2) / 6
      integer(int64) :: bytes = 0        !! device bytes of the arrays below
      integer(int64) :: ntiles = 0       !! tiles of the last run
      integer(int64) :: ngathers = 0     !! integral slices built in the last run
      integer(int64) :: ngemm_builds = 0  !! of those, built by the GEMM from the factors
      integer(int64) :: ncopy_builds = 0  !! of those, permuting copies of the other layout of the same orbital
      integer :: nslot = 0               !! A' slices the cache holds, 4..2 no
      integer :: nblock = 0              !! occupied orbitals per block: no when every slice fits, else nslot / 4
      integer(int64) :: tick = 0         !! block triples begun: the clock of the cache
      integer :: parts = 1               !! work parts the pairs are dealt to (module header, WORK SHARE)
      integer :: my_part = 0             !! this run's part, 0-based
      integer(int64) :: ncomputed = 0    !! triples this part computed in the last run
      integer, allocatable :: owner(:)   !! part (0-based) of the pair (i,j), at `rccsd_triples_ijk_pair_index`; host only
      integer, allocatable :: orb_slot(:)
      !! slot holding the slice of key o + no (layout - 1), layout 1 A' and 2 A2', 0 for none (2 no)
      integer, allocatable :: slot_orb(:)            !! key in each slot, 0 for none (nslot)
      integer(int64), allocatable :: slot_used(:)    !! tick of the last triple that used each slot (nslot)
      real(dp) :: seconds_wide = 0.0_dp, seconds_per = 0.0_dp, seconds_energy = 0.0_dp  !! device seconds, last run
      real(dp) :: seconds_gather = 0.0_dp                              !! repacks and slice gathers, last run
      real(dp) :: flops_wide = 0.0_dp, flops_per = 0.0_dp  !! counted flops of the wide and the per-array gemms, last run
      real(dp), allocatable :: tvb(:)    !! B' of every (o2,o3), (nv + no, nv, no, no): t2(o3,o2,w,f), then -ooov(o2,m,o3,w)
      real(dp), allocatable :: tvb2(:)   !! tvb2(:,:,o3,o2) = tvb(:,:,o2,o3): runs over o3 contiguous
      real(dp), allocatable :: gs(:)
      !! slots, (nv, nv, nv + no, nslot): A'(x,y,:) = [(o x|f y) | t2(m,o,y,x)] or A2'(x,y,:) = A'(y,x,:)
      real(dp), allocatable :: bo(:)     !! (nv, naux): b_ov(o,x;P) of the orbital whose slice is being built
      real(dp), allocatable :: q(:)      !! the Q tile, (nv^3, 6, ncap)
      real(dp), allocatable :: part(:)   !! row sums of the energy, (nv[c], nv[b], ncap)
      real(dp), allocatable :: vm(:)     !! per triple of the tile ovov(i,p,j,q), ovov(i,p,k,q), ovov(j,p,k,q): (nv, nv, 3, ncap)
      real(dp), allocatable :: vt(:)     !! per triple of the tile t1(i,p), t1(j,p), t1(k,p): (nv, 3, ncap)
      real(dp), allocatable :: col(:)    !! column sums of the rows, (nv, ncap)
      real(dp), allocatable :: etrip(:)  !! E_ijk of each triple of the tile, ncap
      integer, allocatable :: trip(:, :)  !! (i,j,k) of each triple of the tile, 3 x ncap
      integer, allocatable :: amap(:, :)
      !! per triple of the tile: for U1..U3 the array that holds it, then 1 (0 for a skipped i = j = k triple), 4 x ncap
      logical :: diagonal = .false.      !! compute the i = j = k triples (`diagonal_triples`)
   contains
      procedure :: init => scratch_init
      procedure :: release => scratch_release
   end type trc_rccsd_triples_ijk_scratch_t

contains

   pure function rccsd_triples_ijk_triple_bytes(nv) result(bytes)
      !! Bytes of Q for one occupied triple, 24 nv^3: the least a budget must hold (`init` refuses below it, and the
      !! API warns of it before the CCSD iterations)
      integer, intent(in) :: nv
      integer(int64) :: bytes

      bytes = 24_int64*nv*nv*nv
   end function rccsd_triples_ijk_triple_bytes

   subroutine scratch_init(this, ints, opts)
      !! ---------------------------------------------------------------------
      !! the scratch
      !! ---------------------------------------------------------------------
      !! Plans the tile and maps the arrays; a refusal is in `this%status` and leaves the type empty
      class(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_triples_dev_options_t), intent(in) :: opts

      integer(int64) :: no, nv, nv3, budget, per_triple, slice, other
      integer :: ierr

      call this%release()
      call this%status%clear()
      if (ints%stamp == 0 .or. .not. allocated(ints%b_ov) .or. .not. allocated(ints%b_vv) .or. .not. allocated(ints%ooov) .or. &
          .not. allocated(ints%ovov) .or. .not. allocated(ints%eps_o) .or. .not. allocated(ints%eps_v)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: the integral blocks are empty (no successful build)")
         return
      end if
      if (ints%no < 1 .or. ints%nv < 1) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: no occupied or no virtual orbitals: there are no triples")
         return
      end if
      if (opts%budget_bytes < 0 .or. opts%tile_triples < 0 .or. opts%slice_bytes < 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: a negative budget or tile size")
         return
      end if
      if (opts%parts < 1 .or. opts%part < 0 .or. opts%part >= opts%parts) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: parts "//to_string(opts%parts)//", part "// &
                              to_string(opts%part)//": need parts >= 1 and 0 <= part < parts")
         return
      end if
      no = ints%no
      nv = ints%nv
      nv3 = nv*nv*nv
      if (nv*nv > int(huge(0_int32), int64) .or. no*no > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: a gemm dimension would pass 2^31-1 ("// &
                              to_string(ints%no)//" occupied, "//to_string(ints%nv)//" virtual)")
         return
      end if
      this%no = ints%no
      this%nv = ints%nv
      this%ntotal = no*(no + 1)*(no + 2)/6
      per_triple = rccsd_triples_ijk_triple_bytes(ints%nv)
      if (opts%budget_bytes > 0) then
         budget = opts%budget_bytes
      else
         budget = max(DEFAULT_BUDGET, per_triple)
      end if
      if (opts%tile_triples > 0) then
         this%ncap = min(opts%tile_triples, this%ntotal)
      else
         this%ncap = min(budget/per_triple, this%ntotal)
      end if
      if (this%ncap < 1) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: the budget of "//to_string(budget)// &
                              " bytes does not hold one occupied triple ("//to_string(per_triple)//" bytes of Q)")
         this%no = 0
         this%nv = 0
         return
      end if
      ! every device byte of the plan but the slices (module header, DEVICE MEMORY); the slice cache takes what is left
      slice = 8_int64*nv*nv*(nv + no)
      other = 8_int64*(2*(nv + no)*nv*no*no + 3_int64*nv3*this%ncap + nv*ints%naux + nv*nv*this%ncap + &
                       3*nv*nv*this%ncap + 4*nv*this%ncap + this%ncap) + 28_int64*this%ncap
      this%nslot = int(min(max(slice_cache_bytes(opts%slice_bytes, other)/slice, 4_int64), 2*no))
      if (this%nslot >= 2*no) then
         this%nblock = this%no
      else
         this%nblock = max(1, this%nslot/4)
      end if
      allocate (this%tvb((nv + no)*nv*no*no), stat=ierr)
      if (ierr == 0) allocate (this%tvb2((nv + no)*nv*no*no), stat=ierr)
      if (ierr == 0) allocate (this%gs(this%nslot*nv*nv*(nv + no)), stat=ierr)
      if (ierr == 0) allocate (this%bo(nv*ints%naux), stat=ierr)
      if (ierr == 0) allocate (this%orb_slot(2*no), stat=ierr)
      if (ierr == 0) allocate (this%slot_orb(this%nslot), stat=ierr)
      if (ierr == 0) allocate (this%slot_used(this%nslot), stat=ierr)
      if (ierr == 0) allocate (this%q(3_int64*nv3*this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%part(nv*nv*this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%vm(nv*nv*3*this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%vt(nv*3*this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%col(nv*this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%etrip(this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%trip(3, this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%amap(4, this%ncap), stat=ierr)
      this%diagonal = opts%diagonal_triples
      this%parts = opts%parts
      this%my_part = opts%part
      if (ierr == 0) then
         call rccsd_triples_ijk_pair_owners(ints%no, opts%parts, opts%diagonal_triples, this%owner)
         if (.not. allocated(this%owner)) ierr = 1
      end if
      if (ierr /= 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: the host arrays of the plan could not be allocated")
         call scratch_free(this)
         return
      end if
      !$acc enter data create(this%tvb, this%tvb2, this%gs, this%q, this%part, this%vm, this%vt, this%col, this%etrip, this%trip, &
      !$acc                   this%amap, this%bo)
      this%stamp = ints%stamp
      this%bytes = other + this%nslot*slice
   end subroutine scratch_init

   function slice_cache_bytes(requested, other) result(bytes)
      !! Bytes for the slice cache: `requested` when above 0; else, on a device build, SLICE_FREE_FRACTION of the free
      !! device memory (at most all of it less SLICE_RESERVE) less `other` (the rest of the plan), queried before the
      !! plan is mapped, so that whatever the CCSD left mapped is already counted; 16 GiB on the host or when the
      !! query gives nothing.
      !! The caller clamps the result to 4..2 no slices
#ifdef TRC_DEVICE
      use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, acc_property_free_memory
#endif
      integer(int64), intent(in) :: requested, other
      integer(int64) :: bytes
#ifdef TRC_DEVICE
      integer(int64) :: free
#endif

      if (requested > 0) then
         bytes = requested
         return
      end if
      bytes = DEFAULT_SLICE_BYTES
#ifdef TRC_DEVICE
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
      ! the whole plan within SLICE_FREE_FRACTION of what is free, and never closer than SLICE_RESERVE to all of it
      if (free > 0) bytes = max(min(int(SLICE_FREE_FRACTION*real(free, dp), int64), free - SLICE_RESERVE) - other, 0_int64)
#endif
   end function slice_cache_bytes

   subroutine scratch_free(this)
      !! Frees what is allocated; the device copies only for the arrays `init` mapped (all or none)
      class(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: this

      if (this%stamp /= 0) then
         !$acc exit data delete(this%tvb, this%tvb2, this%gs, this%q, this%part, this%vm, this%vt, this%col, &
         !$acc                  this%etrip, this%trip, &
         !$acc                  this%amap, this%bo)
      end if
      if (allocated(this%tvb)) deallocate (this%tvb)
      if (allocated(this%bo)) deallocate (this%bo)
      if (allocated(this%tvb2)) deallocate (this%tvb2)
      if (allocated(this%gs)) deallocate (this%gs)
      if (allocated(this%q)) deallocate (this%q)
      if (allocated(this%part)) deallocate (this%part)
      if (allocated(this%col)) deallocate (this%col)
      if (allocated(this%vm)) deallocate (this%vm)
      if (allocated(this%vt)) deallocate (this%vt)
      if (allocated(this%etrip)) deallocate (this%etrip)
      if (allocated(this%trip)) deallocate (this%trip)
      if (allocated(this%orb_slot)) deallocate (this%orb_slot)
      if (allocated(this%slot_orb)) deallocate (this%slot_orb)
      if (allocated(this%slot_used)) deallocate (this%slot_used)
      this%nslot = 0
      this%nblock = 0
      if (allocated(this%amap)) deallocate (this%amap)
      if (allocated(this%owner)) deallocate (this%owner)
      this%parts = 1
      this%my_part = 0
      this%no = 0
      this%nv = 0
      this%stamp = 0
      this%ncap = 0
      this%ntotal = 0
      this%bytes = 0
   end subroutine scratch_free

   subroutine scratch_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: this

      call scratch_free(this)
   end subroutine scratch_release

   subroutine rccsd_triples_ijk_run(ints, amps, e_t, status, opts, seconds, triple_energies)
      !! ---------------------------------------------------------------------
      !! the entry points
      !! ---------------------------------------------------------------------
      !! E(T) in Hartree on the device from the DEVICE copies of `ints` and `amps`, driven by occupied triples
      !!
      !! Refused, with `e_t` = 0 and the reason in `status`: empty integrals, amplitudes missing or of another size
      !! than the integrals, no occupied or no virtual orbitals, a budget that does not hold one occupied triple, a
      !! gemm dimension above 2^31-1, and a result that is NaN or infinite. `triple_energies`, when present, is
      !! allocated to no (no + 1) (no + 2) / 6 in the order of the occupied triples i >= j >= k (i = 1..no, j = 1..i,
      !! k = 1..j, k fastest); each entry is m_o E_ijk / 3, and E(T) is their sum taken in that order.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      real(dp), intent(out) :: e_t
      type(error_t), intent(out) :: status
      type(trc_rccsd_triples_dev_options_t), intent(in), optional :: opts
      real(dp), intent(out), optional :: seconds  !! device wall seconds of the whole call
      real(dp), intent(out), allocatable, optional :: triple_energies(:)

      type(trc_rccsd_triples_dev_options_t) :: o
      type(trc_rccsd_triples_ijk_scratch_t) :: ts
      type(timer_type) :: clock
      type(trc_linalg_t) :: la

      e_t = 0.0_dp
      if (present(seconds)) seconds = 0.0_dp
      if (present(opts)) o = opts
      call la%init(1)
      call clock%start()
      call ts%init(ints, o)
      if (ts%status%has_error()) then
         call status%set(ERROR_VALIDATION, ts%status%get_message())
      else
         call rccsd_triples_ijk_run_scratch(ints, amps, ts, la, e_t, status, triple_energies)
         call la%sync()
         if (present(seconds)) seconds = clock%get_elapsed_time()
         if (o%verbose .and. .not. status%has_error()) call log_summary(ints, ts, e_t, clock%get_elapsed_time())
      end if
      call ts%release()
      call la%release()
   end subroutine rccsd_triples_ijk_run

   subroutine rccsd_triples_ijk_run_scratch(ints, amps, ts, la, e_t, status, triple_energies)
      !! One (T) with a scratch that `ts%init` planned for these integrals; see `rccsd_triples_ijk_run`
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      real(dp), intent(out) :: e_t
      type(error_t), intent(out) :: status
      real(dp), intent(out), allocatable, optional :: triple_energies(:)

      character(len=:), allocatable :: why
      type(timer_type) :: clock
      real(dp), allocatable :: e_canon(:)
      integer(int64) :: n
      integer :: s, i, j, k, nt, ib, jb, kb, nblocks, k1
      real(dp) :: mo, e_sum
      logical :: cached

      e_t = 0.0_dp
      why = inputs_refusal(ints, amps, ts)
      if (len(why) > 0) then
         call status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: "//why)
         return
      end if
      if (present(triple_energies)) then
         if (allocated(triple_energies)) deallocate (triple_energies)
         allocate (triple_energies(ts%ntotal))
      end if
      ts%seconds_wide = 0.0_dp
      ts%seconds_per = 0.0_dp
      ts%seconds_energy = 0.0_dp
      ts%seconds_gather = 0.0_dp
      ts%flops_wide = 0.0_dp
      ts%flops_per = 0.0_dp
      ts%ntiles = 0
      ts%ngathers = 0
      ts%ngemm_builds = 0
      ts%ncopy_builds = 0
      ts%ncomputed = 0
      call la%sync()
      call clock%start()
      call kernel_pack_tvb(ts%no, ts%nv, amps%t2, ints%ooov, ts%tvb)
      call kernel_pack_tvb2(ts%no, ts%nv, ts%tvb, ts%tvb2)
      call la%sync()
      ts%seconds_gather = ts%seconds_gather + clock%get_elapsed_time()
      ! The cache starts empty on every call: the amplitudes or integrals may have changed since the last.
      ts%orb_slot = 0
      ts%slot_orb = 0
      ts%slot_used = 0
      ts%tick = 0
      ! E_ijk at the triple's canonical position, so the sum below is in canonical order whatever the visiting order
      allocate (e_canon(ts%ntotal))
      e_canon = 0.0_dp
      nblocks = (ts%no + ts%nblock - 1)/ts%nblock
      do ib = 1, nblocks
         do jb = 1, ib
            do kb = 1, jb
               if (.not. block_owned(ts, ib, jb, kb)) cycle   ! no slice is built for another part's block triple
               ! every slice the block triple reads is cached before its first tile, and none is evicted until its end
               call pin_block(ints, amps, ts, la, ib, jb, kb, cached)
               if (.not. cached) then
                  call status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: the slice cache cannot hold the "// &
                                  "slices of a block triple ("//to_string(ts%nslot)//" slots, blocks of "// &
                                  to_string(ts%nblock)//")")
                  return
               end if
               do i = block_first(ts, ib), block_last(ts, ib)
                  do j = block_first(ts, jb), min(block_last(ts, jb), i)
                     if (ts%owner(rccsd_triples_ijk_pair_index(i, j)) /= ts%my_part) cycle
                     ! A tile is a run of consecutive k at one (i,j): cut at the end of the block KB or at capacity.
                     k = block_first(ts, kb)
                     k1 = min(block_last(ts, kb), j)
                     do while (k <= k1)
                        nt = int(min(ts%ncap, int(k1 - k + 1, int64)))
                        do s = 1, nt
                           ts%trip(:, s) = [i, j, k + s - 1]
                        end do
                        call run_products(ts, la, nt)
                        !$acc update device(ts%trip(1:3, 1:nt), ts%amap(1:4, 1:nt))
                        call tile_energy(ints, amps, ts, la, nt)
                        !$acc update self(ts%etrip(1:nt))
                        do s = 1, nt
                           e_canon(canonical_index(i, j, k + s - 1)) = ts%etrip(s)
                        end do
                        k = k + nt
                        ts%ntiles = ts%ntiles + 1
                        ts%ncomputed = ts%ncomputed + count(ts%amap(4, 1:nt) /= 0)
                     end do
                  end do
               end do
            end do
         end do
      end do
      ! Added in triple order (i, then j, then k fastest) with the multiplicity of the ordered triples each one stands for.
      e_sum = 0.0_dp
      n = 0
      do i = 1, ts%no
         do j = 1, i
            do k = 1, j
               n = n + 1
               mo = 6.0_dp
               if (i == k) then
                  mo = 1.0_dp
               else if (i == j .or. j == k) then
                  mo = 3.0_dp
               end if
               e_sum = e_sum + mo*e_canon(n)
               if (present(triple_energies)) triple_energies(n) = mo*e_canon(n)/3.0_dp
            end do
         end do
      end do
      e_t = e_sum/3.0_dp
      if (.not. abs(e_t) <= huge(1.0_dp)) then
         e_t = 0.0_dp
         call status%set(ERROR_VALIDATION, "trc_rccsd_triples_ijk: E(T) is NaN or infinite (a zero denominator or a "// &
                         "non-finite amplitude or integral)")
      end if
   end subroutine rccsd_triples_ijk_run_scratch

   function inputs_refusal(ints, amps, ts) result(why)
      !! Empty when the inputs can be used with the scratch, else the reason
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_ijk_scratch_t), intent(in) :: ts
      character(len=:), allocatable :: why

      why = ""
      if (ts%stamp == 0 .or. .not. allocated(ts%q)) then
         why = "the scratch is not initialised"
      else if (ints%stamp == 0 .or. .not. allocated(ints%b_ov) .or. .not. allocated(ints%b_vv) .or. &
               .not. allocated(ints%ooov) .or. .not. allocated(ints%ovov) .or. .not. allocated(ints%eps_o) .or. &
               .not. allocated(ints%eps_v)) then
         why = "the integral blocks are empty (no successful build)"
      else if (ts%stamp /= ints%stamp) then
         why = "the scratch was planned for an earlier build of the integrals (stale stamp): init it again"
      else if (ts%no /= ints%no .or. ts%nv /= ints%nv) then
         why = "the scratch has the sizes of other integrals"
      else if (.not. allocated(amps%t1) .or. .not. allocated(amps%t2)) then
         why = "the amplitudes are not allocated"
      else if (amps%no /= ints%no .or. amps%nv /= ints%nv) then
         why = "the amplitudes are of another size than the integrals"
      else if (any(shape(amps%t1) /= [ints%no, ints%nv]) .or. any(shape(amps%t2) /= [ints%no, ints%no, ints%nv, ints%nv])) then
         why = "the amplitude arrays do not have the shape (no, nv) and (no, no, nv, nv)"
      else if (any(shape(ints%b_ov) /= [ints%no*ints%nv, ints%naux]) .or. &
               any(shape(ints%b_vv) /= [ints%nv*ints%nv, ints%naux]) .or. &
               any(shape(ints%ooov) /= [ints%no, ints%no, ints%no, ints%nv]) .or. &
               any(shape(ints%ovov) /= [ints%no, ints%nv, ints%no, ints%nv]) .or. &
               size(ints%eps_o) /= ints%no .or. size(ints%eps_v) /= ints%nv) then
         why = "the integral blocks do not have the shapes of (no, nv)"
      end if
      if (len(why) > 0) why = trim(why)//" ("//to_string(ints%no)//" occupied, "//to_string(ints%nv)//" virtual)"
   end function inputs_refusal

   subroutine run_products(ts, la, nt)
      !! ---------------------------------------------------------------------
      !! one triple's steps
      !! ---------------------------------------------------------------------
      !! The arrays of the `nt` triples (i,j,k0..k1) of the tile: the products whose operands do not depend on k are one
      !! wide gemm over the run, the two whose slot does (A2'(k)) one gemm per triple (module header). Every slice the
      !! tile reads is cached (`pin_block`)
      type(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: nt

      type(timer_type) :: clock
      integer :: i, j, k, k0, kb, k2, s, nv, no, kk, sa, sb, sc, sk
      real(dp) :: v

      nv = ts%nv
      no = ts%no
      kk = nv + no
      i = ts%trip(1, 1)
      j = ts%trip(2, 1)
      k0 = ts%trip(3, 1)
      ! The array map: U2 = U1 when k = j, U3 = U2 when i = j; an i = j = k triple (k = j = i, the end of an i = j
      ! run) is skipped unless asked for.
      do s = 1, nt
         k = ts%trip(3, s)
         ts%amap(1, s) = 1
         ts%amap(2, s) = merge(1, 2, k == j)
         ts%amap(3, s) = merge(ts%amap(2, s), 3, i == j)
         ts%amap(4, s) = merge(0, 1, i == k .and. .not. ts%diagonal)
      end do
      kb = ts%trip(3, nt)
      if (ts%amap(4, nt) == 0) kb = kb - 1
      if (kb < k0) return
      k2 = min(kb, j - 1)   ! U2 is its own array only for k < j
      v = real(nv, dp)
      ! the run's fixed slices, A'(i), A'(j), A2'(j)
      sa = ts%orb_slot(i)
      sb = ts%orb_slot(j)
      sc = ts%orb_slot(j + no)
      call la%sync()
      call clock%start()
      ! U1 = T(i,j,k) nat [A'(i), tvb(j,k)] + T(j,i,k) swp [A2'(j), tvb(i,k)], k = k0..kb: B' columns from tvb2
      call wide_gemm(ts, la, 1, 1, sa, kb - k0 + 1, 0.0_dp, ts%tvb2(tvb2_pos(ts, k0, j)))
      call wide_gemm(ts, la, 1, 1, sc, kb - k0 + 1, 1.0_dp, ts%tvb2(tvb2_pos(ts, k0, i)))
      ! U2 nat T(i,k,j) [A'(i), tvb(k,j)], k = k0..k2
      if (k2 >= k0) call wide_gemm(ts, la, 2, 1, sa, k2 - k0 + 1, 0.0_dp, ts%tvb(tvb2_pos(ts, k0, j)))
      ! U3 nat T(j,k,i) [A'(j), tvb(k,i)], k = k0..kb, unless U3 = U2 (i = j)
      if (i /= j) call wide_gemm(ts, la, 3, 1, sb, kb - k0 + 1, 0.0_dp, ts%tvb(tvb2_pos(ts, k0, i)))
      call la%sync()
      ts%seconds_wide = ts%seconds_wide + clock%get_elapsed_time()
      ts%flops_wide = ts%flops_wide + 2.0_dp*v*v*v*kk*(2*(kb - k0 + 1) + max(k2 - k0 + 1, 0) + &
                                                       merge(kb - k0 + 1, 0, i /= j))
      ! the swapped terms of U2 and U3 have the slot A2'(k): one gemm per triple
      call clock%start()
      do s = 1, kb - k0 + 1
         k = k0 + s - 1
         sk = ts%orb_slot(k + no)
         ! U2 += T(k,i,j) swp [A2'(k), tvb(i,j)]
         if (k <= k2) call wide_gemm(ts, la, 2, s, sk, 1, 1.0_dp, ts%tvb(tvb2_pos(ts, i, j)))
         ! U3 += T(k,j,i) swp [A2'(k), tvb(j,i)]
         if (i /= j) call wide_gemm(ts, la, 3, s, sk, 1, 1.0_dp, ts%tvb(tvb2_pos(ts, j, i)))
         ts%flops_per = ts%flops_per + 2.0_dp*v*v*v*kk*(merge(1, 0, k <= k2) + merge(1, 0, i /= j))
      end do
      call la%sync()
      ts%seconds_per = ts%seconds_per + clock%get_elapsed_time()
   end subroutine run_products

   pure function tvb2_pos(ts, a, b) result(pos)
      !! Position of column block (:, :, a, b) in tvb or tvb2, both (nv + no, nv, no, no)
      type(trc_rccsd_triples_ijk_scratch_t), intent(in) :: ts
      integer, intent(in) :: a, b
      integer(int64) :: pos

      pos = 1_int64 + int(ts%nv + ts%no, int64)*ts%nv*((a - 1) + int(ts%no, int64)*(b - 1))
   end function tvb2_pos

   subroutine wide_gemm(ts, la, u, s0, slot1, nk, beta, bt)
      !! Frame u of the tile's triples s0..s0+nk-1 (consecutive in storage) += A B, A the slot, B the nk column blocks
      !! of B' that start at `bt` (the caller passes the table from that position): gemm (nv^2 x (nv + no))((nv + no) x nv nk)
      type(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: u, s0, slot1, nk
      real(dp), intent(in) :: beta
      real(dp), intent(in) :: bt(*)
      integer :: nv, no, nv2, kk
      integer(int64) :: nv3, gbase, qbase

      nv = ts%nv
      no = ts%no
      nv2 = nv*nv
      kk = nv + no
      nv3 = int(nv, int64)**3
      gbase = 1_int64 + int(slot1 - 1, int64)*nv2*kk
      qbase = 1_int64 + (int(u - 1, int64)*ts%ncap + (s0 - 1))*nv3
      call la%gemm("N", "N", nv2, nv*nk, kk, 1.0_dp, ts%gs(gbase), nv2, bt, kk, beta, ts%q(qbase), nv2)
   end subroutine wide_gemm

   subroutine pin_block(ints, amps, ts, la, ib, jb, kb, cached)
      !! Caches every slice the block triple (ib,jb,kb) reads: A'(o) for o in blocks ib and jb, A2'(o) for o in blocks
      !! jb and kb, at most 4 nblock <= nslot. The cached ones are stamped first, so the builds that follow evict only
      !! slices this block triple does not read; layout 1 is built before layout 2 of the same orbital, so the second
      !! is a copy. `cached` is false when a slot could not be found (impossible while 4 nblock <= nslot)
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: ib, jb, kb
      logical, intent(out) :: cached

      type(timer_type) :: clock
      integer :: o, layout, ob, key

      cached = .true.
      ts%tick = ts%tick + 1
      do o = 1, ts%no
         ob = (o - 1)/ts%nblock + 1
         do layout = 1, 2
            key = o + ts%no*(layout - 1)
            if (block_reads(ob, layout, ib, jb, kb) .and. ts%orb_slot(key) /= 0) ts%slot_used(ts%orb_slot(key)) = ts%tick
         end do
      end do
      call la%sync()
      call clock%start()
      do o = 1, ts%no
         ob = (o - 1)/ts%nblock + 1
         do layout = 1, 2
            if (.not. block_reads(ob, layout, ib, jb, kb)) cycle
            if (ensure_slot(ints, amps, ts, o, layout, la) == 0) then
               cached = .false.
               return
            end if
         end do
      end do
      call la%sync()
      ts%seconds_gather = ts%seconds_gather + clock%get_elapsed_time()
   end subroutine pin_block

   pure function block_reads(ob, layout, ib, jb, kb) result(reads)
      !! Whether the block triple (ib,jb,kb) reads the slice of layout `layout` of an orbital of block `ob`: A'(i),
      !! A'(j) and A2'(j), A2'(k) (module header, THREE FRAMES)
      integer, intent(in) :: ob, layout, ib, jb, kb
      logical :: reads

      if (layout == 1) then
         reads = ob == ib .or. ob == jb
      else
         reads = ob == jb .or. ob == kb
      end if
   end function block_reads

   pure function block_first(ts, b) result(o)
      !! The first occupied orbital of block b
      type(trc_rccsd_triples_ijk_scratch_t), intent(in) :: ts
      integer, intent(in) :: b
      integer :: o

      o = (b - 1)*ts%nblock + 1
   end function block_first

   pure function block_last(ts, b) result(o)
      !! The last occupied orbital of block b (the last block may be shorter)
      type(trc_rccsd_triples_ijk_scratch_t), intent(in) :: ts
      integer, intent(in) :: b
      integer :: o

      o = min(b*ts%nblock, ts%no)
   end function block_last

   pure function canonical_index(i, j, k) result(n)
      !! The position of the triple i >= j >= k in the order i = 1..no, j = 1..i, k = 1..j, k fastest: the triples
      !! before i are sum_{i' < i} i' (i' + 1) / 2 = (i - 1) i (i + 1) / 6, those of i before j are (j - 1) j / 2
      integer, intent(in) :: i, j, k
      integer(int64) :: n

      n = int(i - 1, int64)*i*(i + 1)/6 + int(j - 1, int64)*j/2 + k
   end function canonical_index

   pure function rccsd_triples_ijk_pair_index(i, j) result(p)
      !! The position of the pair i >= j in the order i = 1..no, j = 1..i
      integer, intent(in) :: i, j
      integer :: p

      p = (i - 1)*i/2 + j
   end function rccsd_triples_ijk_pair_index

   pure function rccsd_triples_ijk_pair_cost(i, j, diagonal) result(c)
      !! The triples the run computes for the pair (i,j), i >= j: k = 1..j, less the i = j = k one unless `diagonal`
      integer, intent(in) :: i, j
      logical, intent(in) :: diagonal
      integer :: c

      c = j
      if (i == j .and. .not. diagonal) c = j - 1
   end function rccsd_triples_ijk_pair_cost

   pure subroutine rccsd_triples_ijk_pair_owners(no, parts, diagonal, owner)
      !! The part (0-based) that computes each pair, at `rccsd_triples_ijk_pair_index`, size no (no + 1) / 2: pairs by
      !! decreasing cost, ties in canonical order, each to the part with the fewest triples so far (ties to the
      !! lowest). Deterministic, so every rank recomputes the same table (module header, WORK SHARE). `owner` is left
      !! unallocated if the memory is not there
      integer, intent(in) :: no, parts
      logical, intent(in) :: diagonal
      integer, intent(out), allocatable :: owner(:)

      integer(int64), allocatable :: load(:)
      integer, allocatable :: pi(:), pj(:)
      integer :: c, i, n, m, p, ierr

      allocate (owner(no*(no + 1)/2), stat=ierr)
      if (ierr == 0) allocate (load(0:max(parts, 1) - 1), stat=ierr)
      if (ierr == 0) allocate (pi(no + 1), stat=ierr)
      if (ierr == 0) allocate (pj(no + 1), stat=ierr)
      if (ierr /= 0) then
         if (allocated(owner)) deallocate (owner)
         return
      end if
      owner = 0
      load = 0
      if (parts <= 1) return
      do c = no, 0, -1
         ! the pairs of cost c in canonical order: (i,c) for i >= c (diagonal) or i > c, with (c+1,c+1) after (c+1,c)
         ! (no diagonal); at c = 0 only (1,1) without it
         n = 0
         if (diagonal) then
            do i = c, no
               if (c < 1) exit
               n = n + 1
               pi(n) = i
               pj(n) = c
            end do
         else
            do i = c + 1, no
               if (c >= 1) then
                  n = n + 1
                  pi(n) = i
                  pj(n) = c
               end if
               if (i == c + 1) then
                  n = n + 1
                  pi(n) = c + 1
                  pj(n) = c + 1
               end if
            end do
         end if
         do m = 1, n
            p = minloc(load, 1) - 1
            owner(rccsd_triples_ijk_pair_index(pi(m), pj(m))) = p
            load(p) = load(p) + c
         end do
      end do
   end subroutine rccsd_triples_ijk_pair_owners

   function block_owned(ts, ib, jb, kb) result(owned)
      !! Whether the block triple (ib,jb,kb) holds a computed triple of a pair (i,j) of this part: i in IB, j in JB,
      !! j <= i, and a k in KB with k <= j (k <= j - 1 for i = j when the i = j = k triple is skipped). Always true for
      !! one part
      type(trc_rccsd_triples_ijk_scratch_t), intent(in) :: ts
      integer, intent(in) :: ib, jb, kb
      logical :: owned
      integer :: i, j, kmax

      owned = .true.
      if (ts%parts <= 1) return
      owned = .false.
      do i = block_first(ts, ib), block_last(ts, ib)
         do j = block_first(ts, jb), min(block_last(ts, jb), i)
            if (ts%owner(rccsd_triples_ijk_pair_index(i, j)) /= ts%my_part) cycle
            kmax = min(block_last(ts, kb), j)
            if (i == j .and. .not. ts%diagonal) kmax = min(kmax, j - 1)
            if (block_first(ts, kb) <= kmax) then
               owned = .true.
               return
            end if
         end do
      end do
   end function block_owned

   function ensure_slot(ints, amps, ts, o, layout, la) result(slot)
      !! The slot holding the slice (o, layout), built if missing and stamped with the current tick; 0 when every slot
      !! holds a slice stamped with it (the block triple's)
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: ts
      integer, intent(in) :: o, layout
      type(trc_linalg_t), intent(in) :: la
      integer :: slot, key

      key = o + ts%no*(layout - 1)
      slot = cache_slot(ts, key)
      if (slot == 0) return
      if (ts%slot_orb(slot) /= key) then
         if (ts%slot_orb(slot) /= 0) ts%orb_slot(ts%slot_orb(slot)) = 0
         call build_slice(ints, amps, ts, o, layout, slot, la)
         ts%slot_orb(slot) = key
         ts%orb_slot(key) = slot
         ts%ngathers = ts%ngathers + 1
      end if
      ts%slot_used(slot) = ts%tick
   end function ensure_slot

   subroutine build_slice(ints, amps, ts, o, layout, slot, la)
      !! Slot `slot` = A'_o (layout 1) or A2'_o (layout 2) from the fitted factors (module header, SLICE CACHE):
      !! the (o x|f y) columns are 2 nv^3 naux flops of GEMM, the t2 columns a gather
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: ts
      integer, intent(in) :: o, layout, slot
      type(trc_linalg_t), intent(in) :: la

      integer(int64) :: g0, other0
      integer :: nv, nv2, other

      nv = ts%nv
      nv2 = nv*nv
      g0 = 1_int64 + int(slot - 1, int64)*nv2*(nv + ts%no)
      ! the other layout of the same orbital, if cached: its integral columns are this one's with x and y exchanged
      other = ts%orb_slot(o + ts%no*(2 - layout))
      if (other /= 0) then
         ! a permuting copy moves 2 nv^3 doubles; the GEMM below costs 2 nv^3 naux flops
         other0 = 1_int64 + int(other - 1, int64)*nv2*(nv + ts%no)
         call kernel_copy_swapped(nv, size(ts%gs, kind=int64), other0, g0, ts%gs)
         ts%ncopy_builds = ts%ncopy_builds + 1
      else
         ts%ngemm_builds = ts%ngemm_builds + 1
         ! A'(x,(y,f)) = sum_P b_ov(o,x;P) b_vv((y,f),P) = (o x|y f) = (o x|f y): ONE GEMM, b_vv read in place
         call kernel_gather_bo(ts%no, nv, ints%naux, o, ints%b_ov, ts%bo)
         call la%gemm("N", "T", nv, nv2, ints%naux, 1.0_dp, ts%bo, nv, ints%b_vv, nv2, 0.0_dp, ts%gs(g0), nv)
         ! A2'(x,y,f) = A'(y,x,f): exchanged in place
         if (layout == 2) call kernel_swap_in_place(nv, ts%gs(g0))
      end if
      call kernel_gather_slice_t2(ts%no, nv, o, layout, amps%t2, ts%gs(g0))
   end subroutine build_slice

   function cache_slot(ts, key) result(slot)
      !! The slot for the slice `key`: the one holding it, else an empty one, else the least recently used one not
      !! stamped with the current tick (those hold the block triple's slices); 0 if there is none
      type(trc_rccsd_triples_ijk_scratch_t), intent(in) :: ts
      integer, intent(in) :: key
      integer :: slot
      integer :: k

      slot = ts%orb_slot(key)
      if (slot /= 0) return
      do k = 1, ts%nslot
         if (ts%slot_orb(k) == 0) then
            slot = k
            return
         end if
      end do
      do k = 1, ts%nslot
         if (ts%slot_used(k) == ts%tick) cycle
         if (slot == 0) then
            slot = k
         else if (ts%slot_used(k) < ts%slot_used(slot)) then
            slot = k
         end if
      end do
   end function cache_slot

   subroutine tile_energy(ints, amps, ts, la, nt)
      !! The energy kernels over the `nt` triples of the tile; the sums stay on the device
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_ijk_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: nt

      type(timer_type) :: clock

      call la%sync()
      call clock%start()
      call kernel_gather_vm(ts%no, ts%nv, nt, ts%ncap, ints%ovov, ts%trip, ts%vm)
      call kernel_gather_vt(ts%no, ts%nv, nt, ts%ncap, amps%t1, ts%trip, ts%vt)
      call kernel_energy_rows(ts%no, ts%nv, nt, ts%ncap, ts%q, ts%vm, ts%vt, ints%eps_o, ints%eps_v, ts%trip, &
                              ts%amap, ts%part)
      call kernel_energy_cols(ts%nv, nt, ts%ncap, ts%part, ts%col)
      call kernel_energy_sum(ts%nv, nt, ts%ncap, ts%col, ts%etrip)
      call la%sync()
      ts%seconds_energy = ts%seconds_energy + clock%get_elapsed_time()
   end subroutine tile_energy

   subroutine log_summary(ints, ts, e_t, seconds)
      !! The verbose table: tiles, triples, device time by stage, counted flops and the rates
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_triples_ijk_scratch_t), intent(in) :: ts
      real(dp), intent(in) :: e_t, seconds

      character(len=240) :: line

      write (line, "(a,i0,a,i0,a,f20.12,a,f9.2,a)") "  rccsd (T) device ijk: no ", ints%no, " nv ", ints%nv, &
         merge("   E(T) of this part ", "   E(T)              ", ts%parts > 1), e_t, "   ", seconds, " s"
      call global_logger%info(trim(line))
      if (ts%parts > 1) then
         write (line, "(a,i0,a,i0,a,i0,a,i0,a)") "    part ", ts%my_part + 1, " of ", ts%parts, ": ", ts%ncomputed, &
            " of the ", ts%ntotal, " triples"
         call global_logger%info(trim(line))
      end if
      write (line, "(a,i0,a,i0,a,i0,a,f9.1,a,i0,a,i0,a)") "    ", ts%ntiles, " tiles of up to ", ts%ncap, &
         " occupied triples (", ts%ntotal, " in all), device scratch ", real(ts%bytes, dp)/1.0e6_dp, " MB with ", &
         ts%nslot, " integral slices cached, blocks of ", ts%nblock, " occupied orbitals"
      call global_logger%info(trim(line))
      write (line, "(a,f9.3,a,f9.3,a,f9.3,a,f9.3,a,i0,a,i0,a,i0,a)") "    seconds: wide gemms ", ts%seconds_wide, &
         "  per-array gemms ", ts%seconds_per, "  energy kernels ", ts%seconds_energy, "  repacks and slice builds ", &
         ts%seconds_gather, " (", ts%ngathers, " slices built: ", ts%ngemm_builds, " by GEMM, ", ts%ncopy_builds, &
         " copies)"
      call global_logger%info(trim(line))
      write (line, "(a,es10.3,a,f9.1,a,es10.3,a,f9.1,a,f9.1,a)") "    flops: wide ", ts%flops_wide, " (", &
         1.0e-9_dp*ts%flops_wide/max(ts%seconds_wide, tiny(1.0_dp)), " GFLOP/s)  per-array ", ts%flops_per, " (", &
         1.0e-9_dp*ts%flops_per/max(ts%seconds_per, tiny(1.0_dp)), " GFLOP/s)  both ", &
         1.0e-9_dp*(ts%flops_wide + ts%flops_per)/max(ts%seconds_wide + ts%seconds_per, tiny(1.0_dp)), " GFLOP/s"
      call global_logger%info(trim(line))
   end subroutine log_summary

   subroutine kernel_pack_tvb(no, nv, t2, ooov, tvb)
      !! ---------------------------------------------------------------------
      !! the kernels. Each takes its arrays as explicit-shape dummies and runs one
      !! `do concurrent` over ONE flattened index, whose body calls an item routine
      !! of this module and assigns one array element. No `x%name` appears in any.
      !! ---------------------------------------------------------------------
      !! tvb(f,w,o2,o3) = t2(o3,o2,w,f) for f <= nv, tvb(nv+m,w,o2,o3) = -ooov(o2,m,o3,w): B' of every (o2,o3)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(int(no, int64)*no*nv*nv)
      real(dp), intent(in) :: ooov(int(no, int64)*no*no*nv)
      real(dp), intent(out) :: tvb(int(nv + no, int64)*nv*no*no)
      integer(int64) :: t, n

      n = int(nv + no, int64)*nv*no*no
      do concurrent(t=1:n)
         tvb(t) = tvb_item(TRC_DC_INDEX(t), no, nv, t2, ooov)
      end do
   end subroutine kernel_pack_tvb

   subroutine kernel_pack_tvb2(no, nv, tvb, tvb2)
      !! tvb2(:,:,a,b) = tvb(:,:,b,a)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tvb(int(nv + no, int64)*nv*no*no)
      real(dp), intent(out) :: tvb2(int(nv + no, int64)*nv*no*no)
      integer(int64) :: t, n

      n = int(nv + no, int64)*nv*no*no
      do concurrent(t=1:n)
         tvb2(t) = tvb2_item(TRC_DC_INDEX(t), no, nv, tvb)
      end do
   end subroutine kernel_pack_tvb2

   subroutine kernel_gather_bo(no, nv, naux, o, b_ov, bo)
      !! bo(x,P) = b_ov(o + no (x-1), P): the strided rows of one occupied orbital, contiguous for the slice GEMM
      integer, intent(in) :: no, nv, naux, o
      real(dp), intent(in) :: b_ov(int(no, int64)*nv*naux)
      real(dp), intent(out) :: bo(int(nv, int64)*naux)
      integer(int64) :: t, n

      n = int(nv, int64)*naux
      do concurrent(t=1:n)
         bo(t) = bo_item(TRC_DC_INDEX(t), no, nv, naux, o, b_ov)
      end do
   end subroutine kernel_gather_bo

   subroutine kernel_copy_swapped(nv, ng, src0, dst0, g)
      !! g(dst0 - 1 + (x,y,f)) = g(src0 - 1 + (y,x,f)) over the nv^3 integral columns of two slots of one array
      !! (one dummy for both: the slots are disjoint parts of the same actual)
      integer, intent(in) :: nv
      integer(int64), intent(in) :: ng, src0, dst0
      real(dp), intent(inout) :: g(ng)
      integer(int64) :: t, n

      n = int(nv, int64)*nv*nv
      do concurrent(t=1:n)
         g(dst0 - 1 + t) = g(src0 - 1 + swapped_pos(TRC_DC_INDEX(t), nv))
      end do
   end subroutine kernel_copy_swapped

   subroutine kernel_swap_in_place(nv, g)
      !! g(x,y,f) <-> g(y,x,f) over the nv^3 integral columns of one slot; one thread per x < y
      integer, intent(in) :: nv
      real(dp), intent(inout) :: g(int(nv, int64)*nv*nv)
      integer(int64) :: t, n

      n = int(nv, int64)*nv*nv
      do concurrent(t=1:n)
         call swap_item(TRC_DC_INDEX(t), nv, g)
      end do
   end subroutine kernel_swap_in_place

   subroutine kernel_gather_slice_t2(no, nv, o, layout, t2, g)
      !! The t2 columns of one slot: A'((x,y),nv+m) = t2(m,o,y,x) (layout 1), A2'((x,y),nv+m) = A'((y,x),nv+m) (layout 2)
      integer, intent(in) :: no, nv, o, layout
      real(dp), intent(in) :: t2(int(no, int64)*no*nv*nv)
      real(dp), intent(inout) :: g(int(nv, int64)*nv*(nv + no))
      integer(int64) :: t, n, off

      n = int(nv, int64)*nv*no
      off = int(nv, int64)*nv*nv
      do concurrent(t=1:n)
         g(off + t) = slice_t2_item(TRC_DC_INDEX(t), no, nv, o, layout, t2)
      end do
   end subroutine kernel_gather_slice_t2

   subroutine kernel_gather_vm(no, nv, nt, ncap, ovov, trip, vm)
      !! vm(:,:,1,s) = ovov(i,:,j,:), vm(:,:,2,s) = ovov(i,:,k,:), vm(:,:,3,s) = ovov(j,:,k,:) for (i,j,k) = trip(:,s)
      integer, intent(in) :: no, nv, nt
      integer(int64), intent(in) :: ncap
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      integer, intent(in) :: trip(3, ncap)
      real(dp), intent(out) :: vm(int(nv, int64)*nv*3*ncap)
      integer(int64) :: t, n

      n = int(nv, int64)*nv*3*nt
      do concurrent(t=1:n)
         vm(t) = vm_item(TRC_DC_INDEX(t), no, nv, ncap, ovov, trip)
      end do
   end subroutine kernel_gather_vm

   subroutine kernel_gather_vt(no, nv, nt, ncap, t1, trip, vt)
      !! vt(:,1,s) = t1(i,:), vt(:,2,s) = t1(j,:), vt(:,3,s) = t1(k,:) for (i,j,k) = trip(:,s)
      integer, intent(in) :: no, nv, nt
      integer(int64), intent(in) :: ncap
      real(dp), intent(in) :: t1(no, nv)
      integer, intent(in) :: trip(3, ncap)
      real(dp), intent(out) :: vt(int(nv, int64)*3*ncap)
      integer(int64) :: t, n

      n = int(nv, int64)*3*nt
      do concurrent(t=1:n)
         vt(t) = vt_item(TRC_DC_INDEX(t), no, nv, ncap, t1, trip)
      end do
   end subroutine kernel_gather_vt

   subroutine kernel_energy_rows(no, nv, nt, ncap, q, vm, vt, eps_o, eps_v, trip, amap, part)
      !! part(c,b,s) = sum_{a>=b} of the energy of (a,b,c) for occupied triple s (0 for c > b), s = 1..nt
      integer, intent(in) :: no, nv, nt
      integer(int64), intent(in) :: ncap
      real(dp), intent(in) :: q(3_int64*nv*nv*nv*ncap)
      real(dp), intent(in) :: vm(nv, nv, 3, ncap), vt(nv, 3, ncap), eps_o(no), eps_v(nv)
      integer, intent(in) :: trip(3, ncap), amap(4, ncap)
      real(dp), intent(out) :: part(int(nv, int64)*nv*ncap)
      integer(int64) :: t, n

      n = int(nv, int64)*nv*nt
      do concurrent(t=1:n)
         part(t) = row_energy_item(TRC_DC_INDEX(t), no, nv, ncap, q, vm, vt, eps_o, eps_v, trip, amap)
      end do
   end subroutine kernel_energy_rows

   subroutine kernel_energy_cols(nv, nt, ncap, part, col)
      !! col(b,s) = sum_c part(c,b,s), in index order
      integer, intent(in) :: nv, nt
      integer(int64), intent(in) :: ncap
      real(dp), intent(in) :: part(int(nv, int64)*nv*ncap)
      real(dp), intent(out) :: col(int(nv, int64)*ncap)
      integer(int64) :: t, n

      n = int(nv, int64)*nt
      do concurrent(t=1:n)
         col(t) = run_sum_item(TRC_DC_INDEX(t), nv, int(nv, int64)*nv*ncap, part)
      end do
   end subroutine kernel_energy_cols

   subroutine kernel_energy_sum(nv, nt, ncap, col, etrip)
      !! etrip(s) = sum_b col(b,s), in index order
      integer, intent(in) :: nv, nt
      integer(int64), intent(in) :: ncap
      real(dp), intent(in) :: col(int(nv, int64)*ncap)
      real(dp), intent(out) :: etrip(ncap)
      integer(int64) :: s

      do concurrent(s=1:int(nt, int64))
         etrip(s) = run_sum_item(TRC_DC_INDEX(s), nv, int(nv, int64)*ncap, col)
      end do
   end subroutine kernel_energy_sum

   pure subroutine split4(t, n1, n2, n3, i1, i2, i3, i4)
      !! ---------------------------------------------------------------------
      !! the item routines: pure, `routine seq`, explicit-shape dummies, scalars
      !! only. A flat 1-based index t is the position in the array that is written,
      !! first index fastest.
      !! A call in a `do concurrent` body passes the loop index through the TRC_DC_INDEX macro:
      !! see inc/trc_dc_index.inc (an ifx miscompilation of the plain variable).
      !! ---------------------------------------------------------------------
      !! t - 1 = (i1-1) + n1 ((i2-1) + n2 ((i3-1) + n3 (i4-1))), 64-bit t
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: n1, n2, n3
      integer, intent(out) :: i1, i2, i3, i4
      integer(int64) :: r

      r = t - 1
      i1 = int(mod(r, int(n1, int64))) + 1
      r = r/n1
      i2 = int(mod(r, int(n2, int64))) + 1
      r = r/n2
      i3 = int(mod(r, int(n3, int64))) + 1
      i4 = int(r/n3) + 1
   end subroutine split4

   pure function tvb_item(t, no, nv, t2, ooov) result(v)
      !! tvb(r,w,o2,o3): t2(o3,o2,w,r) for r <= nv, else -ooov(o2,r-nv,o3,w); position t of (nv + no, nv, no, no)
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(int(no, int64)*no*nv*nv)
      real(dp), intent(in) :: ooov(int(no, int64)*no*no*nv)
      real(dp) :: v
      integer :: r, w, o2, o3

      call split4(t, nv + no, nv, no, r, w, o2, o3)
      if (r <= nv) then
         v = t2(o3 + int(no, int64)*(o2 - 1) + int(no, int64)*no*(w - 1) + int(no, int64)*no*nv*(r - 1))
      else
         v = -ooov(o2 + int(no, int64)*(r - nv - 1) + int(no, int64)*no*(o3 - 1) + int(no, int64)*no*no*(w - 1))
      end if
   end function tvb_item

   pure function tvb2_item(t, no, nv, tvb) result(v)
      !! tvb2(r,w,a,b) = tvb(r,w,b,a), position t of (nv + no, nv, no, no)
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tvb(int(nv + no, int64)*nv*no*no)
      real(dp) :: v
      integer :: r, w, a, b

      call split4(t, nv + no, nv, no, r, w, a, b)
      v = tvb(r + int(nv + no, int64)*(w - 1) + int(nv + no, int64)*nv*((b - 1) + int(no, int64)*(a - 1)))
   end function tvb2_item

   pure function swapped_pos(t, nv) result(q)
      !! The position of (y,x,f) for the position t of (x,y,f) in (nv, nv, nv)
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: nv
      integer(int64) :: q
      integer(int64) :: x, y, f

      x = mod(t - 1, int(nv, int64))
      y = mod((t - 1)/nv, int(nv, int64))
      f = (t - 1)/(int(nv, int64)*nv)
      q = 1 + y + nv*(x + nv*f)
   end function swapped_pos

   pure subroutine swap_item(t, nv, g)
      !! g(x,y,f) <-> g(y,x,f) for the thread that owns the element of the pair with the smaller position: the other
      !! element of the pair has no thread that touches it, and a diagonal one (x = y) moves nowhere
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: nv
      real(dp), intent(inout) :: g(int(nv, int64)*nv*nv)
      integer(int64) :: q
      real(dp) :: held

      q = swapped_pos(t, nv)
      if (q > t) then
         held = g(t)
         g(t) = g(q)
         g(q) = held
      end if
   end subroutine swap_item

   pure function bo_item(t, no, nv, naux, o, b_ov) result(v)
      !! b_ov(o + no (x-1), P) at position t of (nv, naux)
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: no, nv, naux, o
      real(dp), intent(in) :: b_ov(int(no, int64)*nv*naux)
      real(dp) :: v
      integer(int64) :: x, p

      x = mod(t - 1, int(nv, int64)) + 1
      p = (t - 1)/nv + 1
      v = b_ov(o + int(no, int64)*(x - 1) + int(no, int64)*nv*(p - 1))
   end function bo_item

   pure function slice_t2_item(t, no, nv, o, layout, t2) result(v)
      !! A'((x,y),nv+m) = t2(m,o,y,x) of orbital o, with x and y exchanged for layout 2; position t of (nv, nv, no)
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: no, nv, o, layout
      real(dp), intent(in) :: t2(int(no, int64)*no*nv*nv)
      real(dp) :: v
      integer :: x, y, m, unused, u, w

      call split4(t, nv, nv, no, x, y, m, unused)
      if (layout == 1) then
         u = x
         w = y
      else
         u = y
         w = x
      end if
      v = t2(m + int(no, int64)*(o - 1) + int(no, int64)*no*(w - 1) + int(no, int64)*no*nv*(u - 1))
   end function slice_t2_item

   pure function run_sum_item(t, n, nx, x) result(e)
      !! x((t-1) n + 1) + ... + x(t n), added in index order
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer(int64), intent(in) :: nx
      integer, intent(in) :: n
      real(dp), intent(in) :: x(nx)
      real(dp) :: e
      integer :: r

      e = 0.0_dp
      do r = 1, n
         e = e + x((t - 1)*n + r)
      end do
   end function run_sum_item

   pure function w_at(p, q, r, nv, b1, b2, b3, qq, qn) result(w)
      !! W(p,q,r) = U1(p,q,r) + U2(p,r,q) + U3(q,r,p), (p,q,r) the virtuals of (i,j,k) and U1..U3 the arrays of the
      !! frames (i,j,k), (i,k,j), (j,k,i) starting `bn` elements into `qq` (module header)
      !$acc routine seq
      integer, intent(in) :: p, q, r, nv
      integer(int64), intent(in) :: b1, b2, b3, qn
      real(dp), intent(in) :: qq(qn)
      real(dp) :: w
      integer(int64) :: n1, n2

      n1 = int(nv, int64)
      n2 = n1*nv
      w = qq(b1 + p + n1*(q - 1) + n2*(r - 1)) &
          + qq(b2 + p + n1*(r - 1) + n2*(q - 1)) &
          + qq(b3 + q + n1*(r - 1) + n2*(p - 1))
   end function w_at

   pure function vm_item(t, no, nv, ncap, ovov, trip) result(v)
      !! vm(p,q,m,s), position t of (nv, nv, 3, nt): ovov(i,p,j,q), ovov(i,p,k,q), ovov(j,p,k,q) for m = 1, 2, 3
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer(int64), intent(in) :: ncap
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      integer, intent(in) :: trip(3, ncap)
      real(dp) :: v
      integer :: p, q, m, s

      call split4(t, nv, nv, 3, p, q, m, s)
      if (m == 1) then
         v = ovov(trip(1, s), p, trip(2, s), q)
      else if (m == 2) then
         v = ovov(trip(1, s), p, trip(3, s), q)
      else
         v = ovov(trip(2, s), p, trip(3, s), q)
      end if
   end function vm_item

   pure function vt_item(t, no, nv, ncap, t1, trip) result(v)
      !! vt(p,m,s), position t of (nv, 3, nt): t1(o,p) for o the m-th orbital of triple s
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer(int64), intent(in) :: ncap
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv)
      integer, intent(in) :: trip(3, ncap)
      real(dp) :: v
      integer :: p, m, s, unused

      call split4(t, nv, 3, 1, p, m, unused, s)
      v = t1(trip(m, s), p)
   end function vt_item

   pure function row_energy_item(t, no, nv, ncap, qq, vm, vt, eps_o, eps_v, trip, amap) result(e)
      !! sum_{a >= b} 1/(s_v D) sum_P X(P) R(P) for the virtual pair (b,c), b >= c, of one occupied triple: position t
      !! of (nv[c], nv[b], nt[s]); 0 when c > b
      !!
      !! V(p,q,r) = M1(p,q) tk(r) + M2(p,r) tj(q) + M3(q,r) ti(p) from the triple's compact copies (vm, vt), with
      !! everything that depends on (b,c) only read once outside the loop over a; the operations are those of the
      !! three-term V in the module header, in the same order
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer(int64), intent(in) :: ncap
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: qq(3_int64*nv*nv*nv*ncap)
      real(dp), intent(in) :: vm(nv, nv, 3, ncap), vt(nv, 3, ncap), eps_o(no), eps_v(nv)
      integer, intent(in) :: trip(3, ncap), amap(4, ncap)
      real(dp) :: e
      integer :: a, b, c, s, unused
      integer(int64) :: base, qn, n3, b1, b2, b3
      real(dp) :: w1, w2, w3, w4, w5, w6, s_even, s_odd, x, d, sv, eo, epsb, epsc
      real(dp) :: m1bc, m1cb, m2bc, m2cb, m3bc, m3cb, tib, tic, tjb, tjc, tkb, tkc
      real(dp) :: m1ab, m1ac, m1ba, m1ca, m2ab, m2ac, m2ba, m2ca, m3ab, m3ac, m3ba, m3ca, tia, tja, tka

      e = 0.0_dp
      call split4(t, nv, nv, 1, c, b, unused, s)
      if (c > b) return
      if (amap(4, s) == 0) return
      qn = 3_int64*nv*nv*nv*ncap
      n3 = int(nv, int64)*nv*nv
      ! the tile is q(nv^3, ncap, frame)
      base = int(s - 1, int64)*n3
      b1 = base + (amap(1, s) - 1)*ncap*n3
      b2 = base + (amap(2, s) - 1)*ncap*n3
      b3 = base + (amap(3, s) - 1)*ncap*n3
      eo = eps_o(trip(1, s)) + eps_o(trip(2, s)) + eps_o(trip(3, s))
      epsb = eps_v(b)
      epsc = eps_v(c)
      m1bc = vm(b, c, 1, s)
      m1cb = vm(c, b, 1, s)
      m2bc = vm(b, c, 2, s)
      m2cb = vm(c, b, 2, s)
      m3bc = vm(b, c, 3, s)
      m3cb = vm(c, b, 3, s)
      tib = vt(b, 1, s)
      tic = vt(c, 1, s)
      tjb = vt(b, 2, s)
      tjc = vt(c, 2, s)
      tkb = vt(b, 3, s)
      tkc = vt(c, 3, s)
      do a = b, nv
         ! the orderings: even (a,b,c), (b,c,a), (c,a,b); odd (a,c,b), (b,a,c), (c,b,a)
         w1 = w_at(a, b, c, nv, b1, b2, b3, qq, qn)
         w2 = w_at(b, c, a, nv, b1, b2, b3, qq, qn)
         w3 = w_at(c, a, b, nv, b1, b2, b3, qq, qn)
         w4 = w_at(a, c, b, nv, b1, b2, b3, qq, qn)
         w5 = w_at(b, a, c, nv, b1, b2, b3, qq, qn)
         w6 = w_at(c, b, a, nv, b1, b2, b3, qq, qn)
         m1ab = vm(a, b, 1, s)
         m1ac = vm(a, c, 1, s)
         m1ba = vm(b, a, 1, s)
         m1ca = vm(c, a, 1, s)
         m2ab = vm(a, b, 2, s)
         m2ac = vm(a, c, 2, s)
         m2ba = vm(b, a, 2, s)
         m2ca = vm(c, a, 2, s)
         m3ab = vm(a, b, 3, s)
         m3ac = vm(a, c, 3, s)
         m3ba = vm(b, a, 3, s)
         m3ca = vm(c, a, 3, s)
         tia = vt(a, 1, s)
         tja = vt(a, 2, s)
         tka = vt(a, 3, s)
         s_even = w1 + w2 + w3
         s_odd = w4 + w5 + w6
         ! V(p,q,r) = M1(p,q) tk(r) + M2(p,r) tj(q) + M3(q,r) ti(p) at each ordering
         x = (w1 + (m1ab*tkc + m2ac*tjb + m3bc*tia))*(3.0_dp*w1 + s_even - 2.0_dp*s_odd) &
             + (w2 + (m1bc*tka + m2ba*tjc + m3ca*tib))*(3.0_dp*w2 + s_even - 2.0_dp*s_odd) &
             + (w3 + (m1ca*tkb + m2cb*tja + m3ab*tic))*(3.0_dp*w3 + s_even - 2.0_dp*s_odd) &
             + (w4 + (m1ac*tkb + m2ab*tjc + m3cb*tia))*(3.0_dp*w4 + s_odd - 2.0_dp*s_even) &
             + (w5 + (m1ba*tkc + m2bc*tja + m3ac*tib))*(3.0_dp*w5 + s_odd - 2.0_dp*s_even) &
             + (w6 + (m1cb*tka + m2ca*tjb + m3ba*tic))*(3.0_dp*w6 + s_odd - 2.0_dp*s_even)
         sv = 1.0_dp
         if (a == c) then
            sv = 6.0_dp
         else if (a == b .or. b == c) then
            sv = 2.0_dp
         end if
         d = eo - eps_v(a) - epsb - epsc
         e = e + x/(sv*d)
      end do
   end function row_energy_item

end module trc_rccsd_triples_ijk
