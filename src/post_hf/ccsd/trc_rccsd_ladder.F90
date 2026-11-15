!
! The particle-particle ladder of the t2 equation of a closed-shell RI-CCSD(T)
! calculation, in its three forms, and the t1 dressing that goes with it.
!
! Ported from the same author's mqc_czt_rcc.f90 (ladder_accumulate,
! ri_ladder_pair, ri_ladder_direct, ladder_t1_dressing, the fitted branch of
! particle_ladder). Chemists' notation, spatial orbitals, REAL orbitals; the
! integral blocks, tau and the amplitudes are those of trc_rccsd_ints and
! trc_rccsd_inter.
!
! WHAT IS COMPUTED. With the fitted (ac|bd) = sum_P b_vv(a,c; P) b_vv(b,d; P)
! (b_vv(ab,P) has a fastest), mqc's ladder and its dressing add to the new t2,
! before any symmetrisation and before the denominators, exactly
!
!   ladder    L(i,j,a,b)  = sum_cd (ac|bd) tau(i,j,c,d)           einsum 'ijcd,acbd->ijab'
!   dressing  D(i,j,a,b)  = - sum_cd [ sum_k ovvv(k,d,a,c) t1(k,b)
!                                    + sum_k ovvv(k,c,b,d) t1(k,a) ] tau(i,j,c,d)
!   total     L + D       = sum_cd Wvvvv(a,b,c,d) tau(i,j,c,d)    einsum 'abcd,ijcd->ijab'
!
! with Wvvvv(a,b,c,d) = (ac|bd) - sum_k ovvv(k,d,a,c) t1(k,b) - sum_k ovvv(k,c,b,d) t1(k,a)
! (pyscf/cc/rintermediates.py cc_Wvvvv: einsum('kdac,kb->abcd', ovvv, -t1) and
! einsum('kcbd,ka->abcd', ovvv, t1) added to eris.vvvv.transpose(0,2,1,3)).
! So L + D is the `t2new += lib.einsum('abcd,ijcd->ijab', Wvvvv, tau)` line of
! pyscf/cc/rccsd.py update_amps: no factor, no symmetrisation. PySCF's
! `mycc._add_vvvv(t1, t2, eris)` (ccsd.py _add_vvvv_full, einsum
! 'ijcd,acbd->ijab' on tau = t2 + t1 t1) is only L, the UNDRESSED ladder: the
! ovvv terms are applied elsewhere in pyscf's update_amps. The test therefore
! compares L with the `ladder` probe, and L + D, D with two further probes
! (`ladder_w`, `dressing`) that generate.py writes from cc_Wvvvv.
!
! SYMMETRY. L, D and L + D satisfy R(i,j,a,b) = R(j,i,b,a) when tau does
! (tau(j,i,d,c) = tau(i,j,c,d), true for t2 = MP2 and every converged t2): for L
! swap the roles of c and d and use (bd|ac) = (ac|bd); the two terms of D are
! mapped into each other by the same swap.
!
! LAYOUT. Results are in mqc's (i,j,a,b) layout, first index fastest, in
! `trc_rccsd_ladder_t%h`. `rccsd_ladder_pp` OVERWRITES it with L,
! `rccsd_ladder_dressing` ADDS D to it (so on a freshly initialised, zeroed
! type it gives D alone), `rccsd_ladder_build` is the two in that order.
!
! THE THREE FORMS (o = active occupied, v = virtual, P = auxiliary functions; MACs
! are multiply-adds, flops are twice that). None ever holds (ab|cd) whole:
!
!   ASSEMBLED  for a batch of nb compound columns cd = c + v (d-1)
!                W(a,b,cd) = (ac|bd) = sum_P b_vv(a,c;P) b_vv(b,d;P)
!              is one strided batched GEMM per run of columns that share d
!              (m = n = v, k = naux, the c-th A item is the rows (c-1)v+1..c v of
!              b_vv, B is the rows of d for every item): the raw (ab|cd) = one GEMM
!              of b_vv against b_vv would give (ab|cd), not the (ac|bd) the ladder
!              needs, and re-sorting it would cost a pass over v^2 nb. Then
!                L(ij; ab) += tau(ij; cd-batch) . W(ab; cd-batch)^T
!              is ONE GEMM per batch (m = o^2, n = v^2, k = nb; tau and the
!              result are the matrices their memory already is).
!                MACs  naux v^4 + o^2 v^4     scratch  v^2 cd_batch
!              The default cd_batch is the largest the memory budget allows.
!   DIRECT     for each auxiliary function P: G_P(a,c) = b_vv(a,c;P),
!                L(ij; ab) += sum_P G_P tau_ij G_P^T,     tau_ij(c,d) = tau(i,j,c,d)
!              restructured for the GPU as three steps per batch of np auxiliary
!              functions instead of o^2 np small products:
!                1. taut(c, d ij) = tau(i,j,c,d) (once per call);
!                   Y(a, d ij; P) = sum_c G_P(a,c) taut(c, d ij): ONE strided
!                   batched GEMM over P with B shared (m = v, n = v o^2, k = v);
!                2. a permutation pack Y2(ij, a, d, P) = Y(a, d ij; P);
!                3. L(ij a; b) += sum_{d,P} Y2(ij a; d P) b_vv(b, d; P)
!                   -- ONE GEMM (m = o^2 v, n = v, k = v np), because b_vv's
!                   columns P..P+np-1 ARE the matrix (v; d P) with ld v.
!                MACs  2 naux o^2 v^3          scratch  taut + 2 o^2 v^2 np
!              (the pack moves o^2 v^2 np doubles twice per batch).
!   HALF       W(a,b,c,d) = (ac|bd) = W(b,a,d,c), so the c < d half of the sum is the
!              c > d half with (a,b) and the tau indices swapped. With w_cd = 1 for
!              c > d and 1/2 for c = d,
!                N(ij,ab)  = sum_{c>=d} w_cd W(a,b,c,d) tau(i,j,c,d)
!                N'(ij,ab) = sum_{c>=d} w_cd W(a,b,c,d) tau(i,j,d,c)
!                L(i,j,a,b) = N(i,j,a,b) + N'(i,j,b,a)              for any tau,
!              and when tau(i,j,d,c) = tau(j,i,c,d), N'(i,j,b,a) = N(j,i,b,a):
!                L(i,j,a,b) = N(i,j,a,b) + N(j,i,b,a)               one contraction.
!              Only the c >= d columns of W are built. Per pass of d = d0..d0+nd-1:
!                1. C((a,c'),(b,d')) = sum_P b_vv((a,c),P) b_vv((b,d),P), c = d0..v:
!                   ONE GEMM (m = v (v-d0+1), n = v nd, k = naux), both operands
!                   row ranges of b_vv read in place (lda = ldb = v^2). The rows
!                   d0 <= c < d are computed and dropped (6.6 % of the build at
!                   v 265, d_batch 8);
!                2. one kernel re-sorts the c >= d columns, weight folded in, to
!                   W2((a,b), col), and one gathers taus(ij, col) of tau;
!                3. N(ij; ab) += taus(ij; col) W2(ab; col)^T, ONE GEMM (m = o^2,
!                   n = v^2, k = ncol = the pass's c >= d columns, ~v nd).
!              After the last pass one kernel forms L from N (and N'). The
!              one-contraction path is taken only when the caller says tau has the
!              symmetry (`tau_symmetric`, default false); it gathers the symmetrised
!              (tau(i,j,c,d) + tau(j,i,d,c))/2, so for any tau it is the ladder of
!              that symmetrised tau. The general path gathers tau(i,j,d,c) as well
!              and runs a second contraction into N'. N accumulates in `taut`, which
!              the undressed ladder does not otherwise use (general path: N in the
!              result, N' in `taut`).
!                MACs  naux v^4 / 2 + o^2 v^4 / 2 (symmetric tau; general tau o^2 v^4)
!                scratch  v^3 nd (C) + v^2 ncol (W2) + 2 o^2 ncol (the tau gathers),
!                         ncol = nd v - nd (nd-1)/2 at d_batch; 299 MB per d at v 265
!              The first pass takes d_batch values of d; each later pass, having fewer
!              c rows, takes as many d as the same C and W2 hold (`half_pass_size`),
!              so the contraction's k stays near v d_batch instead of shrinking
!              through the tail (fewer, larger GEMMs; more dropped rows in the tail).
!              d_batch comes from its own default budget, 2.5 GiB, because the GEMMs
!              need n = v d_batch in the thousands (measurements at DEFAULT_HALF_BUDGET).
!   DRESSED BUILD (`rccsd_ladder_build` with the half form). With the dressed factor
!                bt(a,c;P) = b_vv(a,c;P) - sum_k t1(k,a) b_ov(k,c;P)
!              sum_P bt(a,c;P) bt(b,d;P) = Wvvvv(a,b,c,d) + sum_kl t1(k,a) t1(l,b) ovov(k,c,l,d)
!              (the fitted (kc|bd) = ovvv(k,c,b,d) and (kc|ld) = ovov(k,c,l,d) of the
!              same factors), and this Wt keeps Wt(a,b,c,d) = Wt(b,a,d,c), the one
!              symmetry the half form uses (bt itself is not symmetric in (a,c)). So
!                L + D = [the half form on bt] - sum_kl t1(k,a) t1(l,b) Y(k,l,i,j),
!                Y(k,l,i,j) = sum_cd ovov(k,c,l,d) tau(i,j,c,d)
!              and the dressing's 2 o^3 v^3 MACs are not spent. bt is a copy of b_vv and
!              ONE GEMM (m = v, n = v naux, k = o: b_ov((k,c),P) is the matrix (k; (c,P)));
!              Y is ONE GEMM (o^2, o^2, v^2) against ovov packed once by `init` as
!              ((k,l);(c,d)), then G(ij,k; b) = Y t1 (o^3, v, o) and X(ij; a,b) =
!              G_b t1 (o^2, v, o, batched over b), subtracted by one kernel (the
!              symmetric path subtracts the X of the symmetrised tau, matching its passes).
!                MACs  o v^2 naux + o^4 v^2 + o^4 v + o^3 v^2 on top of the passes
!                scratch  bt v^2 naux, ovov pack o^2 v^2, Y o^4, G o^3 v
!              `rccsd_ladder_pp` (L) and `rccsd_ladder_dressing` (D) keep their own
!              routes; the other forms build L + D as pp then dressing.
!   The choice is mqc's between direct and the rest: on counted MACs, (naux + o^2) v^4
!   against 2 naux o^2 v^3 (`rccsd_ladder_prefers_direct`, unchanged), so direct wins
!   when v > 2 naux o^2 / (naux + o^2), about 2 o^2 for naux >> o^2. Otherwise AUTO
!   takes the half form when one d of it (16 v^3 + 16 o^2 v bytes) fits its budget, else
!   the assembled one. The direct threshold was not re-measured against the half form.
!   `options%form` forces any of the three.
!
!   The dressing (`rccsd_ladder_dressing`, and the build of the assembled and direct forms) is O(o^3 v^3):
!     Z1(k,a; ij) = sum_cd ovvv(k,d,a,c) tau(i,j,c,d)     D1(ij,a,b) = - sum_k Z1(k,a; ij) t1(k,b)
!     Z2(k,b; ij) = sum_cd ovvv(k,c,b,d) tau(i,j,c,d)     D2(ij,a,b) = - sum_k Z2(k,b; ij) t1(k,a)
!   Per virtual a (b) the pack pk(k; cd) = ovvv(k,d,a,c) (ovvv(k,c,b,d)) is
!   o v^2, then Z = pk . taut (m = o, n = o^2, k = v^2) and cmat(ij; b) =
!   Z^T t1 (m = o^2, n = v, k = o); `virtual_batch` virtuals go in one strided
!   batched GEMM each (B shared), and a scatter kernel subtracts cmat. No (ia|bc)
!   block exists (trc_rccsd_ints): the packs of a batch are ONE GEMM of the factors,
!   pk(k,x,y,w) = sum_P b_ov((k,x),P) b_vv((y,w),P) = (kx|wy) (m = o v, n = v nvb,
!   k = naux; the rows (y,w) of b_vv for the batch's w are contiguous because w is
!   its slow index). Half 2 reads it as (kc|bd) at (k,c,d,b), half 1 as (kd|ac) at
!   (k,d,c,a): the pair comes in the order (d,c), so half 1 contracts it with tau
!   packed as taut(x,y,ij) = tau(i,j,y,x), exact for any tau. A pack for all
!   virtuals at once would be an o v^3 array, which is why it is chunked.
!                MACs  2 o^3 v^3 + 2 o^3 v^2 + 2 o v^3 naux   (both halves)
!                scratch  taut + (o v^2 + o^3 + o^2 v) virtual_batch
!   Flops are twice the MACs: `rccsd_ladder_flops` for the undressed form that ran (its counter, `rccsd_ladder_macs`,
!   matches the GEMMs the assembled and direct forms call for any batch plan; `ls%flops_pp` tallies the calls to show
!   it) and `rccsd_ladder_dressing_flops` for the dressing, whose MACs `rccsd_ladder_macs` does not include. The half
!   form's tally exceeds its `rccsd_ladder_flops` by the rows its passes compute and drop, which depend on the batch
!   plan (the step table reports the tally); its dressed build tallies bt and the Y correction as the dressing.
!
! NUMBERS for no = 25, nv = 130, naux = 630 (default budget 512 MiB):
!     assembled   717 GFLOP (build 360 + contraction 357),  wbatch 3970 columns = 537 MB
!     direct     3460 GFLOP (4.8 x the assembled; the selector picks assembled),  np = 3: 2 * 254 MB
!     dressing    138 GFLOP (19 % of the assembled), virtual_batch 129 = 537 MB
!     taut        84.5 MB, the result 84.5 MB.
! MEASURED on GH200 (H100), water pentamer cc-pVTZ / cc-pVTZ-RIFIT, no = 20, nv = 265,
! naux = 705, one undressed ladder per CCSD iteration (synchronised timers, GPU otherwise idle):
!     assembled   0.336 s (10.9 TFLOP, 32 TFLOP/s)
!     half        0.150 s at d_batch 16, symmetric tau: build 0.080 s (46 TFLOP/s), re-sort and
!                 gather 0.014 s, contraction 0.056 s (35 TFLOP/s); 0.161 s at d_batch 8, 0.240 s at 1
!     (fixed passes of d_batch d; the adaptive passes of `half_pass_size` took 0.152 s at d_batch 8,
!     contraction 0.056 s, and 0.153 s at 14, contraction 0.054 s = 36 TFLOP/s against 35 fixed)
!     the ladder stage of the 14-iteration CCSD, ladder + dressing:
!     5.30 s assembled; 2.71 s half (pp + dressing, 0.043 s per iteration of it); 2.12 to 2.13 s
!     half dressed build (bt 0.6 ms and the Y correction 0.9 ms per iteration); energies unchanged
!     to the printed 1e-12.
!
! DETERMINISM. The batches are visited in ascending order and each GEMM
! accumulates into the result with beta = 1, the scatter kernels write each
! element once: on the device the order of every sum is fixed. On the host a
! threaded BLAS regroups the sums inside a GEMM, which changes round-off only
! (assembled against direct: 1.7e-17 on one thread, 2.1e-17 on four). Changing a
! batch size regroups the sums over cd (assembled), over (d,P) (direct) or over the
! c >= d columns (half) inside the GEMM: round-off only. The half form sums the two
! halves (c > d and c < d) in a different order from the other forms, also round-off
! only (3e-17 against the assembled one in the test), and its symmetric path gives a
! result with R(i,j,a,b) = R(j,i,b,a) exactly in the ladder part.
!
! RESIDENCY. Every array of the types here exists on the host AND on the device
! (`enter data` in `init`, `exit data` in `release`, on the same `x%name`
! expression). The device copy is the one the kernels write; `pull` brings the
! result back, `push` sends the host copy up. Each component is mapped on its own
! and the parent type is never mapped, so a kernel gets an array only as an
! explicit-shape dummy of a `kernel_*` routine and never names `x%name` in a
! `do concurrent` body (an illegal address under mem:separate that the host build
! hides). The host copies of the scratch are never read.
!
! STALE SCRATCH. The scratch holds no integral-derived data, but its sizes and
! batch plan are those of one `ints%build`: it records `ints%stamp` and every
! routine refuses (message in `ls%status`, nothing written) after the integrals
! were rebuilt. Re-`init` it.
!
! SHARES (over the ranks of a run, trc_rccsd_driver RANKS). The half form's passes run over d = ls%d_lo..ls%d_hi
! (`init` sets 1..nv; an empty range, d_hi < d_lo, runs no pass), and the ladder is linear in that range: N, N' and the
! symmetrising or transposing kernel that finishes the passes (h = N + N^T) are linear in N, so the partial arrays of the
! ranks add to the whole one, and the t1 t1 correction, which uses only tau and t1, is subtracted once, from the sum.
! `rccsd_ladder_partial` is the bt and the passes of one rank's range (the symmetrising kernel included, before any sum, so
! a one-rank run, which has no sum, keeps its bits), `rccsd_ladder_correct` the correction;
! `rccsd_ladder_build` of the half form is the two in a row on the whole range. A pass never straddles two ranges: a
! range is cut into passes by `half_pass_size` from its own first d and each pass is capped at its last. Only the half
! form can be shared (`ready` refuses another form with a range other than 1..nv). Sums regrouped by a range boundary
! differ from the one-range sums by round-off, like a change of d_batch (TOL_BATCH in the check).
! Failures are reported in `ls%status` (a routine whose arguments differ in
! (no, nv, naux), are not built, or hold arrays of other sizes sets it and
! returns without writing anything); none stops. A block of more than 2^31-1
! elements is refused by `init`, not wrapped.
#include "inc/trc_dc_index.inc"

module trc_rccsd_ladder
   use pic_types, only: dp, int32, int64
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t
   use trc_rccsd_steps, only: RCCSD_ST_LADDER, RCCSD_ST_DRESSING, rccsd_steps_gemm_flops, rccsd_steps_open, rccsd_steps_close
   use pic_timer, only: timer_type
   use pic_strings, only: to_string
   use trc_step_table, only: step_table_t
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   character(len=*), parameter :: STALE_MSG = "trc_rccsd_ladder: the integrals were rebuilt after the scratch was "// &
                                              "initialised; re-init the scratch after ints%build"

   integer, parameter, public :: RCCSD_LADDER_AUTO = 0       !! `options%form`: pick by counted operations
   integer, parameter, public :: RCCSD_LADDER_ASSEMBLED = 1  !! `options%form`: build (ac|bd) in column batches
   integer, parameter, public :: RCCSD_LADDER_DIRECT = 2     !! `options%form`: contract G_P tau G_P^T, never form (ac|bd)
   integer, parameter, public :: RCCSD_LADDER_HALF = 3       !! `options%form`: build only the c >= d columns of (ac|bd)
   ! What the half form's tau gather reads for the column (c,d): tau(i,j,c,d), tau(i,j,d,c), or the symmetrised
   ! (tau(i,j,c,d) + tau(j,i,d,c))/2 of the one-contraction path
   integer, parameter :: HALF_TAU_CD = 0, HALF_TAU_DC = 1, HALF_TAU_SYM = 2
   integer(int64), parameter :: DEFAULT_BUDGET = 536870912_int64  !! bytes per working set when `budget_bytes` is 0
   ! The half form's C + W2 when `budget_bytes` is 0: its GEMMs are v^2 x v d_batch x naux and o^2 x v^2 x ~v d_batch.
   ! On GH200 (no 20, nv 265, naux 705, 299 MB per d) one dressed build with adaptive passes took 0.161 / 0.152 /
   ! 0.153 s at d_batch 4 / 8 / 14: 2.5 GiB gives d_batch 8 there, as fast as 14 in half the memory.
   integer(int64), parameter :: DEFAULT_HALF_BUDGET = 2684354560_int64

   public :: trc_rccsd_ladder_t, trc_rccsd_lscratch_t, trc_rccsd_ladder_options_t
   public :: rccsd_ladder_pp, rccsd_ladder_dressing, rccsd_ladder_build, rccsd_ladder_partial, rccsd_ladder_correct
   public :: rccsd_ladder_column_flops, rccsd_ladder_correction_flops
   public :: rccsd_ladder_prefers_direct, rccsd_ladder_macs, rccsd_ladder_flops
   public :: rccsd_ladder_dressing_macs, rccsd_ladder_dressing_flops

   type :: trc_rccsd_ladder_options_t
      !! What `init` of the scratch needs to know besides the integrals. Plain values only. A zero is "choose".
      integer :: form = RCCSD_LADDER_AUTO  !! RCCSD_LADDER_AUTO, _ASSEMBLED, _DIRECT or _HALF
      integer :: cd_batch = 0              !! compound (cd) columns per pass of the assembled form
      integer :: d_batch = 0               !! virtuals d per pass of the half form
      integer :: p_batch = 0               !! auxiliary functions per pass of the direct form
      integer :: virtual_batch = 0         !! virtuals per pass of the dressing
      integer(int64) :: budget_bytes = 0   !! bytes per working set (the form's, the dressing's); 0 means 512 MiB, 2.5 GiB
                                           !! for the half form's C + W2
   end type trc_rccsd_ladder_options_t

   type :: trc_rccsd_ladder_t
      !! The ladder contribution on host and device, valid after the matching routine
      type(error_t) :: status  !! set by `init`: a refusal leaves its message here and the type empty
      integer :: no = 0
      integer :: nv = 0
      real(dp), allocatable :: h(:, :, :, :)  !! h(i,j,a,b), zeroed by `init`
   contains
      procedure :: init => ladder_init
      procedure :: release => ladder_release
      procedure :: push => ladder_push
      procedure :: pull => ladder_pull
   end type trc_rccsd_ladder_t

   type :: trc_rccsd_lscratch_t
      !! The work arrays of the routines below, sizes in the module header
      type(error_t) :: status  !! set by `init` and by every routine that refuses its arguments
      integer :: no = 0
      integer :: nv = 0
      integer :: naux = 0
      integer :: stamp = 0          !! `ints%stamp` of the integrals `init` planned for
      integer :: form = 0           !! the form in use: RCCSD_LADDER_ASSEMBLED, _DIRECT or _HALF (never AUTO after `init`)
      integer :: cd_batch = 0       !! columns per pass, assembled form
      integer :: d_batch = 0        !! virtuals d per pass, half form
      integer :: p_batch = 0        !! auxiliary functions per pass, direct form
      integer :: virtual_batch = 0  !! virtuals per pass, dressing
      integer :: d_lo = 1           !! the half form runs the passes over d = d_lo..d_hi: this rank's SHARE of the ladder
      integer :: d_hi = 0           !! (module header, SHARES). `init` sets the whole range 1..nv; an empty range has d_hi < d_lo
      integer(int64) :: bytes = 0   !! bytes of the work arrays, host and device each
      integer(int64) :: flops_pp = 0        !! 2 m n k summed over the GEMMs the last undressed ladder called, batches counted
      integer(int64) :: flops_dressing = 0  !! the same for the last dressing: the tallies the step table reports
      real(dp), allocatable :: taut(:)    !! o^2 v^2, taut(c,d,ij) = tau(i,j,c,d), per call; the half form's N (or N')
      real(dp), allocatable :: ovpack(:)  !! o v^2 virtual_batch, (kx|wy) per virtual w of the batch (dressing)
      real(dp), allocatable :: zmat(:)    !! o^3 virtual_batch, Z(k; ij) per virtual
      real(dp), allocatable :: cmat(:)    !! o^2 v virtual_batch, Z^T t1 per virtual
      real(dp), allocatable :: wbatch(:)  !! v^2 cd_batch, (ac|bd) per column (assembled only)
      real(dp), allocatable :: gtau(:)    !! o^2 v^2 p_batch, Y(a, d ij; P) (direct only)
      real(dp), allocatable :: gpack(:)   !! o^2 v^2 p_batch, Y2(ij, a, d, P) (direct only)
      real(dp), allocatable :: cbuf(:)    !! v^3 d_batch, C(a,c',b,d') = (ac|bd) for c >= d0 (half only)
      real(dp), allocatable :: w2(:)      !! v^2 ncol, w_cd (ac|bd) on the c >= d columns of a pass (half only)
      real(dp), allocatable :: taus(:)    !! o^2 ncol, tau gathered on the same columns (half only)
      real(dp), allocatable :: taus2(:)   !! o^2 ncol, tau(i,j,d,c) on the same columns, general tau (half only)
      real(dp), allocatable :: bt(:)      !! v^2 naux, the t1-dressed factor bt(a,c;P) of the half form's build
      real(dp), allocatable :: ovovp(:)   !! o^2 v^2, ovov(k,c,l,d) as ((k,l);(c,d)), packed by `init` (half only)
      real(dp), allocatable :: ymat(:)    !! o^4, Y(ij; k,l) = sum_cd tau(i,j,c,d) ovov(k,c,l,d) (half only)
      real(dp), allocatable :: gmat(:)    !! o^3 v, G(ij,k; b) = sum_l Y(ij; k,l) t1(l,b) (half only)
   contains
      procedure :: init => lscratch_init
      procedure :: release => lscratch_release
   end type trc_rccsd_lscratch_t

contains

   subroutine ladder_init(this, no, nv)
      !! ---------------------------------------------------------------------
      !! the types
      !! ---------------------------------------------------------------------
      !! A zeroed result on host and device. A block above 2^31-1 elements is refused in `status`.
      class(trc_rccsd_ladder_t), intent(inout) :: this
      integer, intent(in) :: no, nv

      call this%release()
      call this%status%clear()
      if (no < 0 .or. nv < 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: negative orbital count")
         return
      end if
      if (elements(no, no, nv, nv) > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: a block has more than 2^31-1 elements")
         return
      end if
      this%no = no
      this%nv = nv
      allocate (this%h(no, no, nv, nv))
      this%h = 0.0_dp
      !$acc enter data copyin(this%h)
   end subroutine ladder_init

   subroutine ladder_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_ladder_t), intent(inout) :: this

      if (allocated(this%h)) then
         !$acc exit data delete(this%h)
         deallocate (this%h)
      end if
      this%no = 0
      this%nv = 0
   end subroutine ladder_release

   subroutine ladder_push(this)
      !! Host copy to device
      class(trc_rccsd_ladder_t), intent(inout) :: this

      if (.not. allocated(this%h)) return
      !$acc update device(this%h)
   end subroutine ladder_push

   subroutine ladder_pull(this)
      !! Device copy to host
      class(trc_rccsd_ladder_t), intent(inout) :: this

      if (.not. allocated(this%h)) return
      !$acc update self(this%h)
   end subroutine ladder_pull

   pure function plan_batch(requested, budget, item_elements, upper) result(n)
      !! Items per pass: `requested` if positive, else as many as `budget` bytes hold; never above `upper`
      !! or above what a default-integer index addresses; 0 when not even one item fits
      integer, intent(in) :: requested, upper
      integer(int64), intent(in) :: budget, item_elements
      integer :: n

      integer(int64) :: cap

      cap = int(huge(0_int32), int64)/max(item_elements, 1_int64)
      if (requested > 0) then
         cap = min(cap, int(requested, int64))
      else
         cap = min(cap, max(1_int64, budget/(8_int64*max(item_elements, 1_int64))))
      end if
      cap = min(cap, int(upper, int64))
      n = int(cap)
      if (cap < 1) n = 0
   end function plan_batch

   subroutine lscratch_init(this, ints, options)
      !! Plans the batches and allocates the work arrays for these integrals. Refusals go to `status`.
      class(trc_rccsd_lscratch_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_ladder_options_t), intent(in), optional :: options

      type(trc_rccsd_ladder_options_t) :: opts
      integer :: no, nv, naux
      integer(int64) :: budget, half_budget, half_item, no2, nv2

      call this%release()
      call this%status%clear()
      if (.not. allocated(ints%b_vv) .or. .not. allocated(ints%b_ov) .or. .not. allocated(ints%ovov)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the integral blocks are not built")
         return
      end if
      opts = trc_rccsd_ladder_options_t()
      if (present(options)) opts = options
      if (opts%form < RCCSD_LADDER_AUTO .or. opts%form > RCCSD_LADDER_HALF .or. opts%cd_batch < 0 .or. &
          opts%p_batch < 0 .or. opts%d_batch < 0 .or. opts%virtual_batch < 0 .or. opts%budget_bytes < 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: an unknown form or a negative batch size or budget")
         return
      end if
      no = ints%no
      nv = ints%nv
      naux = ints%naux
      no2 = int(no, int64)*no
      nv2 = int(nv, int64)*nv
      if (no < 1 .or. nv < 1) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: no active occupied or no virtual orbitals")
         return
      end if
      ! the kernels index o^2 v^2 arrays (tau, the result, the packs) with default integers
      if (no2*nv2 > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: an o^2 v^2 array has more than 2^31-1 elements "// &
                              "(default-integer kernel index)")
         return
      end if
      budget = merge(opts%budget_bytes, DEFAULT_BUDGET, opts%budget_bytes > 0)
      half_budget = merge(opts%budget_bytes, DEFAULT_HALF_BUDGET, opts%budget_bytes > 0)
      ! per d of the half form: C (v^3 at d0 = 1, the largest pass) + W2 (at most v^3) + the two tau gathers
      half_item = 2_int64*nv2*nv + 2_int64*no2*nv
      this%no = no
      this%nv = nv
      this%naux = naux
      this%form = opts%form
      if (this%form == RCCSD_LADDER_AUTO) then
         ! Where the counts do not favour the direct form, the half one (measured 2.2 x faster than the assembled one
         ! on GH200, module header) when one d of it fits the budget; the assembled one when it does not.
         if (rccsd_ladder_prefers_direct(no, nv, naux)) then
            this%form = RCCSD_LADDER_DIRECT
         else if (8_int64*half_item <= half_budget) then
            this%form = RCCSD_LADDER_HALF
         else
            this%form = RCCSD_LADDER_ASSEMBLED
         end if
      end if
      this%virtual_batch = plan_batch(opts%virtual_batch, budget, int(no, int64)*nv2 + no2*no + no2*nv, nv)
      if (this%form == RCCSD_LADDER_ASSEMBLED) then
         this%cd_batch = plan_batch(opts%cd_batch, budget, nv2, int(nv2))
      else if (this%form == RCCSD_LADDER_HALF) then
         this%d_batch = plan_batch(opts%d_batch, half_budget, half_item, nv)
      else
         this%p_batch = plan_batch(opts%p_batch, budget, 2_int64*no2*nv2, naux)
      end if
      if (this%virtual_batch < 1 .or. (this%form == RCCSD_LADDER_ASSEMBLED .and. this%cd_batch < 1) .or. &
          (this%form == RCCSD_LADDER_DIRECT .and. this%p_batch < 1) .or. &
          (this%form == RCCSD_LADDER_HALF .and. this%d_batch < 1)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: not even one batch item fits in a default integer index")
         call this%release()
         return
      end if
      this%stamp = ints%stamp
      this%d_lo = 1
      this%d_hi = nv
      call allocate_work(this)
      ! The dressed build's correction reads ovov as ((k,l);(c,d)): it depends on the integrals only, so once here.
      if (this%form == RCCSD_LADDER_HALF) call kernel_pack_ovov(no, nv, ints%ovov, this%ovovp)
   end subroutine lscratch_init

   subroutine allocate_work(this)
      !! The arrays of the plan in `this`, mapped to the device
      type(trc_rccsd_lscratch_t), intent(inout) :: this

      integer(int64) :: no, nv, no2, nv2, ncol

      no = this%no
      nv = this%nv
      no2 = no*no
      nv2 = nv*nv
      allocate (this%taut(no2*nv2))
      allocate (this%ovpack(no*nv2*this%virtual_batch), this%zmat(no2*no*this%virtual_batch), &
                this%cmat(no2*nv*this%virtual_batch))
      !$acc enter data create(this%taut, this%ovpack, this%zmat, this%cmat)
      this%bytes = 8_int64*(no2*nv2 + (no*nv2 + no2*no + no2*nv)*this%virtual_batch)
      if (this%form == RCCSD_LADDER_ASSEMBLED) then
         allocate (this%wbatch(nv2*this%cd_batch))
         !$acc enter data create(this%wbatch)
         this%bytes = this%bytes + 8_int64*nv2*this%cd_batch
      else if (this%form == RCCSD_LADDER_HALF) then
         ncol = half_ncol(this%d_batch, this%nv, 1)
         allocate (this%cbuf(nv2*nv*this%d_batch), this%w2(nv2*ncol), this%taus(no2*ncol), this%taus2(no2*ncol))
         allocate (this%bt(nv2*this%naux), this%ovovp(no2*nv2), this%ymat(no2*no2), this%gmat(no2*no*nv))
         !$acc enter data create(this%cbuf, this%w2, this%taus, this%taus2, this%bt, this%ovovp, this%ymat, this%gmat)
         this%bytes = this%bytes + 8_int64*(nv2*nv*this%d_batch + (nv2 + 2_int64*no2)*ncol + nv2*this%naux + &
                                            no2*nv2 + no2*no2 + no2*no*nv)
      else
         allocate (this%gtau(no2*nv2*this%p_batch), this%gpack(no2*nv2*this%p_batch))
         !$acc enter data create(this%gtau, this%gpack)
         this%bytes = this%bytes + 16_int64*no2*nv2*this%p_batch
      end if
   end subroutine allocate_work

   subroutine lscratch_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_lscratch_t), intent(inout) :: this

      if (allocated(this%taut)) then
         !$acc exit data delete(this%taut)
         deallocate (this%taut)
      end if
      if (allocated(this%ovpack)) then
         !$acc exit data delete(this%ovpack)
         deallocate (this%ovpack)
      end if
      if (allocated(this%zmat)) then
         !$acc exit data delete(this%zmat)
         deallocate (this%zmat)
      end if
      if (allocated(this%cmat)) then
         !$acc exit data delete(this%cmat)
         deallocate (this%cmat)
      end if
      if (allocated(this%wbatch)) then
         !$acc exit data delete(this%wbatch)
         deallocate (this%wbatch)
      end if
      if (allocated(this%gtau)) then
         !$acc exit data delete(this%gtau)
         deallocate (this%gtau)
      end if
      if (allocated(this%gpack)) then
         !$acc exit data delete(this%gpack)
         deallocate (this%gpack)
      end if
      if (allocated(this%cbuf)) then
         !$acc exit data delete(this%cbuf)
         deallocate (this%cbuf)
      end if
      if (allocated(this%w2)) then
         !$acc exit data delete(this%w2)
         deallocate (this%w2)
      end if
      if (allocated(this%taus)) then
         !$acc exit data delete(this%taus)
         deallocate (this%taus)
      end if
      if (allocated(this%taus2)) then
         !$acc exit data delete(this%taus2)
         deallocate (this%taus2)
      end if
      if (allocated(this%bt)) then
         !$acc exit data delete(this%bt)
         deallocate (this%bt)
      end if
      if (allocated(this%ovovp)) then
         !$acc exit data delete(this%ovovp)
         deallocate (this%ovovp)
      end if
      if (allocated(this%ymat)) then
         !$acc exit data delete(this%ymat)
         deallocate (this%ymat)
      end if
      if (allocated(this%gmat)) then
         !$acc exit data delete(this%gmat)
         deallocate (this%gmat)
      end if
      this%no = 0
      this%nv = 0
      this%naux = 0
      this%stamp = 0
      this%form = 0
      this%cd_batch = 0
      this%d_batch = 0
      this%p_batch = 0
      this%virtual_batch = 0
      this%d_lo = 1
      this%d_hi = 0
      this%bytes = 0
      this%flops_pp = 0_int64
      this%flops_dressing = 0_int64
   end subroutine lscratch_release

   pure function rccsd_ladder_macs(form, no, nv, naux) result(macs)
      !! ---------------------------------------------------------------------
      !! cost model and size checks
      !! ---------------------------------------------------------------------
      !! Multiply-adds of one undressed ladder: assembled (naux + o^2) v^4, direct 2 naux o^2 v^3, half
      !! (naux + o^2) v^4 / 2 (symmetric tau; the rows a pass computes and drops and the second contraction of a
      !! general tau are not counted, module header); 0 for another form
      integer, intent(in) :: form, no, nv, naux
      real(dp) :: macs

      macs = 0.0_dp
      if (form == RCCSD_LADDER_ASSEMBLED) macs = (real(naux, dp) + real(no, dp)**2)*real(nv, dp)**4
      if (form == RCCSD_LADDER_DIRECT) macs = 2.0_dp*real(naux, dp)*real(nv, dp)**3*real(no, dp)**2
      if (form == RCCSD_LADDER_HALF) macs = 0.5_dp*(real(naux, dp) + real(no, dp)**2)*real(nv, dp)**4
   end function rccsd_ladder_macs

   pure function rccsd_ladder_flops(form, no, nv, naux) result(flops)
      !! Flops of one undressed ladder, twice `rccsd_ladder_macs`: every GEMM of the assembled and direct forms is
      !! counted in it (the build of (ac|bd) and its contraction, or the two products around the pack), and
      !! `ls%flops_pp` tallies the calls actually made, which agree for any batch plan. For the half form it is the
      !! useful work only; the tally adds the dropped rows of its passes. 0 for another form.
      integer, intent(in) :: form, no, nv, naux
      integer(int64) :: flops

      flops = nint(2.0_dp*rccsd_ladder_macs(form, no, nv, naux), int64)
   end function rccsd_ladder_flops

   pure function rccsd_ladder_dressing_macs(no, nv, naux) result(macs)
      !! Multiply-adds of the dressing, both halves: the pack pk = b_ov b_vv^T is o v * v * naux per virtual,
      !! Z = pk . taut o * o^2 * v^2 per virtual, cmat = Z^T t1 o^2 * v * o per virtual, for v virtuals each:
      !! 2 o v^3 naux + 2 o^3 v^3 + 2 o^3 v^2 (the macs the header quotes)
      integer, intent(in) :: no, nv, naux
      real(dp) :: macs

      macs = 2.0_dp*real(no, dp)*real(nv, dp)**3*real(naux, dp) + 2.0_dp*real(no, dp)**3*real(nv, dp)**3 + &
             2.0_dp*real(no, dp)**3*real(nv, dp)**2
   end function rccsd_ladder_dressing_macs

   pure function rccsd_ladder_dressing_flops(no, nv, naux) result(flops)
      !! Flops of one dressing, twice `rccsd_ladder_dressing_macs`; `ls%flops_dressing` tallies the calls actually made
      integer, intent(in) :: no, nv, naux
      integer(int64) :: flops

      flops = nint(2.0_dp*rccsd_ladder_dressing_macs(no, nv, naux), int64)
   end function rccsd_ladder_dressing_flops

   pure function rccsd_ladder_column_flops(no, nv, naux) result(flops)
      !! Flops of the half form per (c,d) column with c >= d, symmetric tau: the build of C and its contraction,
      !! 2 v^2 naux + 2 o^2 v^2. A virtual d has nv - d + 1 of them, so the passes d_lo..d_hi cost
      !! this times the number of their columns; all of them, v (v + 1) / 2 columns, are
      !! `rccsd_ladder_flops` of the half form times (v + 1) / v
      integer, intent(in) :: no, nv, naux
      real(dp) :: flops

      flops = 2.0_dp*real(nv, dp)**2*(real(naux, dp) + real(no, dp)**2)
   end function rccsd_ladder_column_flops

   pure function rccsd_ladder_correction_flops(no, nv) result(flops)
      !! Flops of the t1 t1 ovov correction of the dressed build (Y, G and X of `dressed_correction`), which only the
      !! rank holding the summed ladder does; bt (2 o v^2 naux) is built by every rank and is not in it
      integer, intent(in) :: no, nv
      integer(int64) :: flops

      flops = rccsd_steps_gemm_flops(no*no, no*no, nv*nv) + rccsd_steps_gemm_flops(no*no*no, nv, no) + &
              rccsd_steps_gemm_flops(no*no, nv, no, nv)
   end function rccsd_ladder_correction_flops

   pure function half_pass_size(d_batch, nv, d0) result(nd)
      !! Virtuals d in the half-form pass starting at d0: as many as the buffers planned for d_batch passes from d0 = 1
      !! hold. C of a pass is v^2 (v - d0 + 1) nd and W2 v^2 half_ncol, so a later pass, with fewer c rows, takes more
      !! d and the contraction's k stays near v d_batch through the tail instead of shrinking with the c rows.
      integer, intent(in) :: d_batch, nv, d0
      integer :: nd
      integer :: nvc

      nvc = nv - d0 + 1
      nd = int(min(int(nvc, int64), (int(nv, int64)*d_batch)/nvc))
      do while (nd > 1 .and. half_ncol(nd, nv, d0) > half_ncol(d_batch, nv, 1))
         nd = nd - 1
      end do
      nd = max(nd, 1)
   end function half_pass_size

   pure function half_ncol(nd, nv, d0) result(n)
      !! Columns (c,d) with c >= d of a half-form pass over d = d0..d0+nd-1: the sum over those d of nv - d + 1
      integer, intent(in) :: nd, nv, d0
      integer(int64) :: n

      n = int(nd, int64)*(nv - d0 + 1) - (int(nd, int64)*(nd - 1))/2
   end function half_ncol

   pure function rccsd_ladder_prefers_direct(no, nv, naux) result(direct)
      !! Whether the direct form needs fewer multiply-adds than the assembled one (mqc's rule: counts only,
      !! because both forms are GEMMs of comparable shape, so the counts predict the order but not the constant)
      integer, intent(in) :: no, nv, naux
      logical :: direct

      direct = rccsd_ladder_macs(RCCSD_LADDER_DIRECT, no, nv, naux) < rccsd_ladder_macs(RCCSD_LADDER_ASSEMBLED, no, nv, naux)
   end function rccsd_ladder_prefers_direct

   pure function elements(n1, n2, n3, n4) result(n)
      !! n1 n2 n3 n4 as a 64-bit count
      integer, intent(in) :: n1, n2, n3, n4
      integer(int64) :: n

      n = int(n1, int64)*n2*n3*n4
   end function elements

   function ready(ints, lad, ls) result(ok)
      !! True when the integrals, the result and the scratch are built, agree on (no, nv, naux), hold arrays of
      !! those sizes and the scratch is not stale; else `ls%status` says which not
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_ladder_t), intent(in) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      logical :: ok

      integer :: no, nv
      integer(int64) :: no2, nv2
      character(len=200) :: text

      ok = .false.
      if (.not. allocated(ints%b_vv) .or. .not. allocated(ints%b_ov) .or. .not. allocated(ints%ovov)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the integral blocks are not built")
         return
      end if
      if (.not. allocated(lad%h)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the result is not initialised")
         return
      end if
      if (.not. allocated(ls%taut) .or. .not. allocated(ls%ovpack) .or. .not. allocated(ls%zmat) .or. &
          .not. allocated(ls%cmat)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the scratch is not initialised")
         return
      end if
      no = ints%no
      nv = ints%nv
      if (lad%no /= no .or. lad%nv /= nv .or. ls%no /= no .or. ls%nv /= nv .or. ls%naux /= ints%naux) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: (no, nv, naux) differ between integrals, result and scratch")
         return
      end if
      no2 = int(no, int64)*no
      nv2 = int(nv, int64)*nv
      if (size(ints%b_vv, kind=int64) /= nv2*ints%naux .or. &
          size(ints%b_ov, kind=int64) /= int(no, int64)*nv*ints%naux .or. &
          size(ints%ovov, kind=int64) /= elements(no, nv, no, nv)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: an integral block does not have the size (no, nv) says")
         return
      end if
      if (size(lad%h, kind=int64) /= no2*nv2) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the result does not have the size (no, nv) says")
         return
      end if
      if (.not. scratch_sized(ls, no2, nv2)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: a scratch array does not have the size its plan says")
         return
      end if
      if (ls%stamp /= ints%stamp) then
         call ls%status%set(ERROR_VALIDATION, STALE_MSG)
         return
      end if
      if (ls%d_lo < 1 .or. ls%d_hi > nv .or. ls%d_hi < ls%d_lo - 1) then
         text = "trc_rccsd_ladder: the share d = "//to_string(ls%d_lo)//".."//to_string(ls%d_hi)// &
                " is not a range inside 1.."//to_string(nv)
         call ls%status%set(ERROR_VALIDATION, trim(text))
         return
      end if
      if (ls%form /= RCCSD_LADDER_HALF .and. (ls%d_lo /= 1 .or. ls%d_hi /= nv)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: a share of the ladder (d range) needs the half form")
         return
      end if
      ok = .true.
   end function ready

   function scratch_sized(ls, no2, nv2) result(ok)
      !! True when every array of the plan is allocated with the size the plan gives
      type(trc_rccsd_lscratch_t), intent(in) :: ls
      integer(int64), intent(in) :: no2, nv2
      logical :: ok

      integer(int64) :: no, vb, ncol

      no = ls%no
      vb = ls%virtual_batch
      ok = size(ls%taut, kind=int64) == no2*nv2 .and. size(ls%ovpack, kind=int64) == no*nv2*vb .and. &
           size(ls%zmat, kind=int64) == no2*no*vb .and. size(ls%cmat, kind=int64) == no2*ls%nv*vb
      if (.not. ok) return
      if (ls%form == RCCSD_LADDER_ASSEMBLED) then
         ok = allocated(ls%wbatch)
         if (ok) ok = size(ls%wbatch, kind=int64) == nv2*ls%cd_batch
      else if (ls%form == RCCSD_LADDER_DIRECT) then
         ok = allocated(ls%gtau) .and. allocated(ls%gpack)
         if (ok) ok = size(ls%gtau, kind=int64) == no2*nv2*ls%p_batch .and. size(ls%gpack, kind=int64) == no2*nv2*ls%p_batch
      else if (ls%form == RCCSD_LADDER_HALF) then
         ncol = half_ncol(ls%d_batch, ls%nv, 1)
         ok = allocated(ls%cbuf) .and. allocated(ls%w2) .and. allocated(ls%taus) .and. allocated(ls%taus2) .and. &
              allocated(ls%bt) .and. allocated(ls%ovovp) .and. allocated(ls%ymat) .and. allocated(ls%gmat)
         if (ok) ok = size(ls%cbuf, kind=int64) == nv2*ls%nv*ls%d_batch .and. size(ls%w2, kind=int64) == nv2*ncol .and. &
                      size(ls%taus, kind=int64) == no2*ncol .and. size(ls%taus2, kind=int64) == no2*ncol .and. &
                      size(ls%bt, kind=int64) == nv2*ls%naux .and. size(ls%ovovp, kind=int64) == no2*nv2 .and. &
                      size(ls%ymat, kind=int64) == no2*no2 .and. size(ls%gmat, kind=int64) == no2*no*ls%nv
      else
         ok = .false.
      end if
   end function scratch_sized

   function ready_tau(ints, inter, ls) result(ok)
      !! True when `inter` holds a tau of the integrals' (no, nv); else `ls%status` says so
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      logical :: ok

      ok = .false.
      if (.not. allocated(inter%tau)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the intermediates (tau) are not initialised")
         return
      end if
      if (inter%no /= ints%no .or. inter%nv /= ints%nv .or. &
          size(inter%tau, kind=int64) /= elements(ints%no, ints%no, ints%nv, ints%nv)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: tau does not have the size (no, nv) says")
         return
      end if
      ok = .true.
   end function ready_tau

   function ready_amps(ints, amps, ls) result(ok)
      !! True when `amps` holds a t1 of the integrals' (no, nv); else `ls%status` says so
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      logical :: ok

      ok = .false.
      if (.not. allocated(amps%t1)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the amplitudes are not initialised")
         return
      end if
      if (amps%no /= ints%no .or. amps%nv /= ints%nv .or. &
          size(amps%t1, kind=int64) /= elements(ints%no, ints%nv, 1, 1)) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: t1 does not have the size (no, nv) says")
         return
      end if
      ok = .true.
   end function ready_amps

   subroutine rccsd_ladder_pp(ints, inter, lad, ls, la, tau_symmetric)
      !! ---------------------------------------------------------------------
      !! the ladder
      !! ---------------------------------------------------------------------
      !! L(i,j,a,b) = sum_cd (ac|bd) tau(i,j,c,d), overwriting `lad%h` (pyscf _add_vvvv, einsum 'ijcd,acbd->ijab')
      !!
      !! The form is `ls%form`. Needs tau at the current amplitudes (`rccsd_inter_build_tau`). `tau_symmetric`
      !! (default false) lets the half form take its one-contraction path, which computes the ladder of the
      !! symmetrised tau, (tau(i,j,c,d) + tau(j,i,d,c))/2: that is L only when tau(j,i,d,c) = tau(i,j,c,d), as
      !! for the tau = t2 + t1 t1 of a CCSD iteration. The other forms are exact for any tau and ignore it.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      logical, intent(in), optional :: tau_symmetric

      logical :: sym

      sym = .false.
      if (present(tau_symmetric)) sym = tau_symmetric
      call ls%status%clear()
      if (.not. ready(ints, lad, ls)) return
      if (.not. ready_tau(ints, inter, ls)) return
      call ladder_pp_run(ints, inter, lad, ls, la, sym)
   end subroutine rccsd_ladder_pp

   subroutine ladder_pp_run(ints, inter, lad, ls, la, tau_symmetric)
      !! `rccsd_ladder_pp` after its checks
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      logical, intent(in) :: tau_symmetric

      ls%flops_pp = 0_int64
      call kernel_zero(ints%no*ints%no*ints%nv*ints%nv, lad%h)
      if (ls%form == RCCSD_LADDER_ASSEMBLED) then
         call ladder_assembled(ints, inter, lad, ls, la)
      else if (ls%form == RCCSD_LADDER_HALF) then
         call ladder_half(ints, inter, lad, ls, la, tau_symmetric, ints%b_vv)
      else
         call ladder_direct(ints, inter, lad, ls, la)
      end if
   end subroutine ladder_pp_run

   subroutine ladder_assembled(ints, inter, lad, ls, la)
      !! The assembled form (module header): per batch of cd columns, W by batched GEMMs, then one GEMM with tau
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: shape
      integer :: no, nv, naux, nv2, cd0, cd1, cd, col, nb, c, d, nrun

      no = ints%no
      nv = ints%nv
      naux = ints%naux
      nv2 = nv*nv
      do cd0 = 1, nv2, ls%cd_batch
         cd1 = min(cd0 + ls%cd_batch - 1, nv2)
         nb = cd1 - cd0 + 1
         ! W(a,b,col) = (ac|bd) = sum_P b_vv(a,c;P) b_vv(b,d;P), one batched call per run of columns with the same d:
         ! item c reads the rows (c-1)v+1..c v of b_vv, B is the rows of d for all of them.
         cd = cd0
         col = 0
         do while (cd <= cd1)
            c = mod(cd - 1, nv) + 1
            d = (cd - 1)/nv + 1
            nrun = min(nv - c + 1, cd1 - cd + 1)
            shape = trc_gemm_batch(nv, nv, naux, nrun, lda=nv2, ldb=nv2, ldc=nv, &
                                   sa=int(nv, int64), sb=0_int64, sc=int(nv2, int64), transb="T")
            call la%gemm_batched(shape, ints%b_vv((c - 1)*nv + 1, 1), ints%b_vv((d - 1)*nv + 1, 1), &
                                 ls%wbatch(col*nv2 + 1))
            ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(shape%m, shape%n, shape%k, shape%count)
            col = col + nrun
            cd = cd + nrun
         end do
         ! L(ij; ab) += tau(ij; cd-batch) W(ab; cd-batch)^T. tau is (o^2, v^2) in memory and the batch is a
         ! contiguous column range of it, starting at the (c,d) of cd0.
         c = mod(cd0 - 1, nv) + 1
         d = (cd0 - 1)/nv + 1
         call la%gemm("N", "T", no*no, nv2, nb, 1.0_dp, inter%tau(1, 1, c, d), no*no, ls%wbatch, nv2, 1.0_dp, lad%h, no*no)
         ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(no*no, nv2, nb)
      end do
   end subroutine ladder_assembled

   subroutine ladder_half(ints, inter, lad, ls, la, tau_symmetric, fac)
      !! The half form (module header): per pass of d = d0..d1, C by one GEMM, its c >= d columns re-sorted with
      !! their weight into W2 and the tau columns gathered, then one GEMM with tau (two for a general tau).
      !! `fac` is the factor whose products sum_P fac(a,c;P) fac(b,d;P) are contracted: b_vv for L, the dressed bt
      !! for the dressed build; both are laid out as b_vv (rows (a,c), a fastest).
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      logical, intent(in) :: tau_symmetric
      real(dp), intent(in) :: fac(ints%nv*ints%nv, ints%naux)

      integer :: no, nv, naux, no2, nv2, d0, nd, nvc, ncol

      no = ints%no
      nv = ints%nv
      naux = ints%naux
      no2 = no*no
      nv2 = nv*nv
      ! N accumulates in taut, which nothing reads during the undressed ladder (the dressing repacks it). A
      ! general tau puts N in lad%h (zeroed by the caller) and N' in taut.
      call kernel_zero(no2*nv2, ls%taut)
      d0 = ls%d_lo
      do while (d0 <= ls%d_hi)
         nvc = nv - d0 + 1
         ! a pass never reaches past the end of this rank's share
         nd = min(half_pass_size(ls%d_batch, nv, d0), ls%d_hi - d0 + 1)
         ncol = int(half_ncol(nd, nv, d0))
         ! C((a,c'),(b,d')) = sum_P b_vv((a,c),P) b_vv((b,d),P), c = d0..nv, d = d0..d0+nd-1: both operands are
         ! contiguous row ranges of b_vv, read in place. The rows d0 <= c < d are computed and dropped below.
         call la%gemm("N", "T", nv*nvc, nv*nd, naux, 1.0_dp, fac((d0 - 1)*nv + 1, 1), nv2, &
                      fac((d0 - 1)*nv + 1, 1), nv2, 0.0_dp, ls%cbuf, nv*nvc)
         ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(nv*nvc, nv*nd, naux)
         call kernel_half_resort(nv, nvc, nd, d0, ncol, ls%cbuf, ls%w2)
         if (tau_symmetric) then
            call kernel_half_gather(no, nv, nvc, nd, d0, ncol, HALF_TAU_SYM, inter%tau, ls%taus)
            call la%gemm("N", "T", no2, nv2, ncol, 1.0_dp, ls%taus, no2, ls%w2, nv2, 1.0_dp, ls%taut, no2)
            ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(no2, nv2, ncol)
         else
            call kernel_half_gather(no, nv, nvc, nd, d0, ncol, HALF_TAU_CD, inter%tau, ls%taus)
            call la%gemm("N", "T", no2, nv2, ncol, 1.0_dp, ls%taus, no2, ls%w2, nv2, 1.0_dp, lad%h, no2)
            call kernel_half_gather(no, nv, nvc, nd, d0, ncol, HALF_TAU_DC, inter%tau, ls%taus2)
            call la%gemm("N", "T", no2, nv2, ncol, 1.0_dp, ls%taus2, no2, ls%w2, nv2, 1.0_dp, ls%taut, no2)
            ls%flops_pp = ls%flops_pp + 2_int64*rccsd_steps_gemm_flops(no2, nv2, ncol)
         end if
         d0 = d0 + nd
      end do
      if (tau_symmetric) then
         call kernel_half_sym(no, nv, ls%taut, lad%h)
      else
         call kernel_half_general(no, nv, ls%taut, lad%h)
      end if
   end subroutine ladder_half

   subroutine ladder_direct(ints, inter, lad, ls, la)
      !! The direct form (module header): per batch of auxiliary functions, GEMM against tau, a pack, GEMM against b_vv
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: shape
      integer :: no, nv, naux, p0, np

      no = ints%no
      nv = ints%nv
      naux = ints%naux
      call kernel_pack_taut(no, nv, inter%tau, ls%taut)
      do p0 = 1, naux, ls%p_batch
         np = min(ls%p_batch, naux - p0 + 1)
         ! Y(a, d ij; P) = sum_c G_P(a,c) taut(c, d ij): G_P is the (v, v) column P of b_vv, taut is shared
         shape = trc_gemm_batch(nv, nv*no*no, nv, np, lda=nv, ldb=nv, ldc=nv, &
                                sa=int(nv, int64)*nv, sb=0_int64, sc=int(nv, int64)*nv*no*no)
         call la%gemm_batched(shape, ints%b_vv(1, p0), ls%taut, ls%gtau)
         ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(shape%m, shape%n, shape%k, shape%count)
         call kernel_pack_y2(no, nv, np, ls%gtau, ls%gpack)
         ! L(ij a; b) += sum_{d,P} Y2(ij a; d P) b_vv(b, d; P): columns p0.. of b_vv are the matrix (v; d P), ld v
         call la%gemm("N", "T", no*no*nv, nv, nv*np, 1.0_dp, ls%gpack, no*no*nv, ints%b_vv(1, p0), nv, 1.0_dp, &
                      lad%h, no*no*nv)
         ls%flops_pp = ls%flops_pp + rccsd_steps_gemm_flops(no*no*nv, nv, nv*np)
      end do
   end subroutine ladder_direct

   subroutine rccsd_ladder_dressing(ints, amps, inter, lad, ls, la)
      !! Adds D(i,j,a,b) = - sum_cd [sum_k ovvv(k,d,a,c) t1(k,b) + sum_k ovvv(k,c,b,d) t1(k,a)] tau(i,j,c,d) to `lad%h`
      !!
      !! pyscf cc_Wvvvv's einsum('kdac,kb->abcd', ovvv, -t1) and einsum('kcbd,ka->abcd', ovvv, t1), contracted with
      !! tau. Needs tau at the current amplitudes and `amps` on the device (the caller pushes after a host edit). It
      !! ADDS: call it on a zeroed result for D alone, or after `rccsd_ladder_pp` for the total.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      call ls%status%clear()
      if (.not. ready(ints, lad, ls)) return
      if (.not. ready_tau(ints, inter, ls)) return
      if (.not. ready_amps(ints, amps, ls)) return
      call ladder_dressing_run(ints, amps, inter, lad, ls, la)
   end subroutine rccsd_ladder_dressing

   subroutine ladder_dressing_run(ints, amps, inter, lad, ls, la)
      !! `rccsd_ladder_dressing` after its checks
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      integer :: half

      ls%flops_dressing = 0_int64
      do half = 1, 2
         ! the pack of half 1 comes out with its (c,d) pair in the order (d,c): tau is packed to match
         if (half == 1) then
            call kernel_pack_taut_swapped(ints%no, ints%nv, inter%tau, ls%taut)
         else
            call kernel_pack_taut(ints%no, ints%nv, inter%tau, ls%taut)
         end if
         call dressing_half(half, ints, amps, lad, ls, la)
      end do
   end subroutine ladder_dressing_run

   subroutine dressing_half(half, ints, amps, lad, ls, la)
      !! half = 1: D1 (free virtual a is the first of the pair), 2: D2 (b is the second); virtual_batch at a time
      integer, intent(in) :: half
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: shape_z, shape_c
      integer :: no, nv, no2, nv2, nov, v0, nvb

      no = ints%no
      nv = ints%nv
      no2 = no*no
      nv2 = nv*nv
      nov = no*nv
      do v0 = 1, nv, ls%virtual_batch
         nvb = min(ls%virtual_batch, nv - v0 + 1)
         ! pk(k,x,y,w) = sum_P b_ov((k,x),P) b_vv((y,w),P) = (kx|wy) for the batch's w: the rows (y,w) of b_vv are
         ! contiguous because w is its SLOW index there. Half 2 reads it as (kc|bd) at (k,c,d,b); half 1 as
         ! (kd|ac) at (k,d,c,a), the pair in the order (d,c), against the tau packed to match
         call la%gemm("N", "T", nov, nv*nvb, ints%naux, 1.0_dp, ints%b_ov, nov, ints%b_vv(1 + nv*(v0 - 1), 1), nv2, &
                      0.0_dp, ls%ovpack, nov)
         ls%flops_dressing = ls%flops_dressing + rccsd_steps_gemm_flops(nov, nv*nvb, ints%naux)
         ! Z(k; ij) = sum_cd pk(k; cd) taut(cd; ij), one item per virtual, taut shared
         shape_z = trc_gemm_batch(no, no2, nv2, nvb, lda=no, ldb=nv2, ldc=no, &
                                  sa=int(no, int64)*nv2, sb=0_int64, sc=int(no, int64)*no2)
         call la%gemm_batched(shape_z, ls%ovpack, ls%taut, ls%zmat)
         ls%flops_dressing = ls%flops_dressing + rccsd_steps_gemm_flops(shape_z%m, shape_z%n, shape_z%k, shape_z%count)
         ! cmat(ij; x) = sum_k Z(k; ij) t1(k,x), t1 shared
         shape_c = trc_gemm_batch(no2, nv, no, nvb, transa="T", lda=no, ldb=no, ldc=no2, &
                                  sa=int(no, int64)*no2, sb=0_int64, sc=int(no2, int64)*nv)
         call la%gemm_batched(shape_c, ls%zmat, amps%t1, ls%cmat)
         ls%flops_dressing = ls%flops_dressing + rccsd_steps_gemm_flops(shape_c%m, shape_c%n, shape_c%k, shape_c%count)
         if (half == 1) then
            call kernel_sub_z1(no, nv, v0, nvb, ls%cmat, lad%h)
         else
            call kernel_sub_z2(no, nv, v0, nvb, ls%cmat, lad%h)
         end if
      end do
   end subroutine dressing_half

   subroutine rccsd_ladder_build(ints, amps, inter, lad, ls, la, stab, tau_symmetric)
      !! L + D = sum_cd Wvvvv(a,b,c,d) tau(i,j,c,d) in `lad%h`, overwriting it (pyscf rccsd.update_amps's
      !! einsum('abcd,ijcd->ijab', Wvvvv, tau)). Every check is made before the first write, so a refusal changes nothing.
      !! With `stab` the undressed ladder and the dressing are two steps of that table (RCCSD_ST_LADDER and
      !! RCCSD_ST_DRESSING), each with its counted flops; the half form's dressed build puts its passes on bt in the
      !! ladder step and bt with the t1 t1 ovov correction in the dressing step. `tau_symmetric` as for
      !! `rccsd_ladder_pp`; the dressing is exact for any tau.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      type(step_table_t), intent(inout), optional :: stab
      logical, intent(in), optional :: tau_symmetric

      type(timer_type) :: tc
      logical :: sym

      sym = .false.
      if (present(tau_symmetric)) sym = tau_symmetric
      call ls%status%clear()
      if (.not. ready(ints, lad, ls)) return
      if (.not. ready_tau(ints, inter, ls)) return
      if (.not. ready_amps(ints, amps, ls)) return
      if (ls%form == RCCSD_LADDER_HALF) then
         call ladder_build_dressed(ints, amps, inter, lad, ls, la, sym, stab)
      else
         call rccsd_steps_open(stab, la, RCCSD_ST_LADDER, tc)
         call ladder_pp_run(ints, inter, lad, ls, la, sym)
         ! The tallies of the GEMMs actually called, not the formulas: the table then stays right when the
         ! GEMMs change (the check compares the tallies with the formulas).
         call rccsd_steps_close(stab, la, RCCSD_ST_LADDER, tc, ls%flops_pp)
         call rccsd_steps_open(stab, la, RCCSD_ST_DRESSING, tc)
         call ladder_dressing_run(ints, amps, inter, lad, ls, la)
         call rccsd_steps_close(stab, la, RCCSD_ST_DRESSING, tc, ls%flops_dressing)
      end if
   end subroutine rccsd_ladder_build

   subroutine rccsd_ladder_partial(ints, amps, inter, lad, ls, la, stab, tau_symmetric)
      !! The first half of the dressed half-form build, for a RANK'S SHARE of the ladder (module header, SHARES): the
      !! dressed factor bt and the passes over d = ls%d_lo..ls%d_hi, finished by the symmetrising kernel, into
      !! `lad%h` (overwritten). That is sum over the share's d of Wt tau, linear in the d range, so the partial arrays of
      !! the ranks add up to the whole ladder of the half form on bt. The t1 t1 ovov correction is NOT in it:
      !! `rccsd_ladder_correct` subtracts it once, from the sum. Needs the half form. Steps and flops as `rccsd_ladder_build`.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      type(step_table_t), intent(inout), optional :: stab
      logical, intent(in) :: tau_symmetric

      call ls%status%clear()
      if (.not. ready(ints, lad, ls)) return
      if (.not. ready_tau(ints, inter, ls)) return
      if (.not. ready_amps(ints, amps, ls)) return
      if (ls%form /= RCCSD_LADDER_HALF) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: a share of the ladder needs the half form")
         return
      end if
      call dressed_partial(ints, amps, inter, lad, ls, la, tau_symmetric, stab)
   end subroutine rccsd_ladder_partial

   subroutine rccsd_ladder_correct(ints, amps, inter, lad, ls, la, stab, tau_symmetric)
      !! The second half of the dressed build: subtracts the t1 t1 ovov correction from the summed ladder in `lad%h`
      !! (it uses only tau and t1, so it is done once, on the rank that holds the sum). Needs
      !! `rccsd_ladder_partial` of the same amplitudes (it leaves bt's flops in `ls%flops_dressing`, closed here).
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      type(step_table_t), intent(inout), optional :: stab
      logical, intent(in) :: tau_symmetric

      call ls%status%clear()
      if (.not. ready(ints, lad, ls)) return
      if (.not. ready_tau(ints, inter, ls)) return
      if (.not. ready_amps(ints, amps, ls)) return
      if (ls%form /= RCCSD_LADDER_HALF) then
         call ls%status%set(ERROR_VALIDATION, "trc_rccsd_ladder: the dressed correction needs the half form")
         return
      end if
      call dressed_correction(ints, amps, inter, lad, ls, la, tau_symmetric, stab)
   end subroutine rccsd_ladder_correct

   subroutine ladder_build_dressed(ints, amps, inter, lad, ls, la, tau_symmetric, stab)
      !! L + D by the half form on the dressed factor (module header, DRESSED BUILD): the half form's passes on bt
      !! give sum_cd Wt tau, then the t1 t1 ovov part of Wt, which Wvvvv does not have, is subtracted
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      logical, intent(in) :: tau_symmetric
      type(step_table_t), intent(inout), optional :: stab

      call dressed_partial(ints, amps, inter, lad, ls, la, tau_symmetric, stab)
      call dressed_correction(ints, amps, inter, lad, ls, la, tau_symmetric, stab)
   end subroutine ladder_build_dressed

   subroutine dressed_partial(ints, amps, inter, lad, ls, la, tau_symmetric, stab)
      !! bt, then the passes of the share `ls%d_lo..ls%d_hi` into `lad%h` (the whole range in a one-rank run)
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      logical, intent(in) :: tau_symmetric
      type(step_table_t), intent(inout), optional :: stab

      type(timer_type) :: tc
      integer :: no, nv, naux, no2, nv2

      ls%flops_pp = 0_int64
      ls%flops_dressing = 0_int64
      no = ints%no
      nv = ints%nv
      naux = ints%naux
      no2 = no*no
      nv2 = nv*nv
      ! bt(a; c,P) = b_vv(a; c,P) - sum_k t1(k,a) b_ov(k; c,P): b_ov((k,c),P) is the matrix (k; (c,P)), ld o
      call rccsd_steps_open(stab, la, RCCSD_ST_DRESSING, tc)
      ! (a rank whose share is empty runs no pass and has no use for bt)
      if (ls%d_hi >= ls%d_lo) then
         call kernel_copy(nv2*naux, ints%b_vv, ls%bt)
         call la%gemm("T", "N", nv, nv*naux, no, -1.0_dp, amps%t1, no, ints%b_ov, no, 1.0_dp, ls%bt, nv)
         ls%flops_dressing = ls%flops_dressing + rccsd_steps_gemm_flops(nv, nv*naux, no)
      end if
      ! its flops go with the correction's at the close of `dressed_correction`
      call rccsd_steps_close(stab, la, RCCSD_ST_DRESSING, tc)
      call rccsd_steps_open(stab, la, RCCSD_ST_LADDER, tc)
      call kernel_zero(no2*nv2, lad%h)
      call ladder_half(ints, inter, lad, ls, la, tau_symmetric, ls%bt)
      call rccsd_steps_close(stab, la, RCCSD_ST_LADDER, tc, ls%flops_pp)
   end subroutine dressed_partial

   subroutine dressed_correction(ints, amps, inter, lad, ls, la, tau_symmetric, stab)
      !! Y, G and X of the module header, subtracted from `lad%h`
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      logical, intent(in) :: tau_symmetric
      type(step_table_t), intent(inout), optional :: stab

      type(trc_gemm_batch_t) :: shape
      type(timer_type) :: tc
      integer :: no, nv, no2, nv2

      no = ints%no
      nv = ints%nv
      no2 = no*no
      nv2 = nv*nv
      call rccsd_steps_open(stab, la, RCCSD_ST_DRESSING, tc)
      ! Y(ij; k,l) = sum_cd tau(ij; cd) ovov(k,c,l,d); G(ij,k; b) = sum_l Y(ij,k; l) t1(l,b);
      ! X(ij; a, b) = sum_k G(ij; k, b) t1(k, a), one item per b, into taut (free once the passes are done)
      call la%gemm("N", "T", no2, no2, nv2, 1.0_dp, inter%tau, no2, ls%ovovp, no2, 0.0_dp, ls%ymat, no2)
      ls%flops_dressing = ls%flops_dressing + rccsd_steps_gemm_flops(no2, no2, nv2)
      call la%gemm("N", "N", no2*no, nv, no, 1.0_dp, ls%ymat, no2*no, amps%t1, no, 0.0_dp, ls%gmat, no2*no)
      ls%flops_dressing = ls%flops_dressing + rccsd_steps_gemm_flops(no2*no, nv, no)
      shape = trc_gemm_batch(no2, nv, no, nv, lda=no2, ldb=no, ldc=no2, &
                             sa=int(no2, int64)*no, sb=0_int64, sc=int(no2, int64)*nv)
      call la%gemm_batched(shape, ls%gmat, amps%t1, ls%taut)
      ls%flops_dressing = ls%flops_dressing + rccsd_steps_gemm_flops(shape%m, shape%n, shape%k, shape%count)
      if (tau_symmetric) then
         call kernel_sub_xsym(no, nv, ls%taut, lad%h)
      else
         call kernel_sub_x(no2*nv2, ls%taut, lad%h)
      end if
      call rccsd_steps_close(stab, la, RCCSD_ST_DRESSING, tc, ls%flops_dressing)
   end subroutine dressed_correction

   subroutine kernel_zero(n, h)
      !! ---------------------------------------------------------------------
      !! the kernels. Each takes its arrays as explicit-shape dummies and runs one
      !! `do concurrent` over ONE flattened index, whose body calls an item routine
      !! of this module. No `x%name` appears in any of them.
      !! ---------------------------------------------------------------------
      !! h = 0
      integer, intent(in) :: n
      real(dp), intent(out) :: h(n)
      integer :: t

      do concurrent(t=1:n)
         h(t) = 0.0_dp
      end do
   end subroutine kernel_zero

   subroutine kernel_pack_taut(no, nv, tau, taut)
      !! taut(c,d,ij) = tau(i,j,c,d), ij = i + o (j-1)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(out) :: taut(nv*nv*no*no)
      integer :: t

      do concurrent(t=1:nv*nv*no*no)
         taut(t) = taut_item(TRC_DC_INDEX(t), no, nv, tau)
      end do
   end subroutine kernel_pack_taut

   subroutine kernel_half_resort(nv, nvc, nd, d0, ncol, cbuf, w2)
      !! w2(a,b,col(c,d)) = w_cd cbuf(a,c',b,d') for c >= d (w_cd = 1/2 on c = d, else 1); the rows c < d are skipped
      integer, intent(in) :: nv, nvc, nd, d0, ncol
      real(dp), intent(in) :: cbuf(nv*nvc*nv*nd)
      real(dp), intent(out) :: w2(nv*nv*ncol)
      integer :: t, dest

      do concurrent(t=1:nv*nvc*nv*nd) local(dest)
         dest = resort_dest(TRC_DC_INDEX(t), nv, nvc, d0)
         if (dest > 0) w2(dest) = cbuf(t)
         if (dest < 0) w2(-dest) = 0.5_dp*cbuf(t)
      end do
   end subroutine kernel_half_resort

   subroutine kernel_half_gather(no, nv, nvc, nd, d0, ncol, mode, tau, taus)
      !! taus(ij,col(c,d)) = tau(i,j,c,d), tau(i,j,d,c) or (tau(i,j,c,d) + tau(j,i,d,c))/2 by `mode`, for c >= d
      integer, intent(in) :: no, nv, nvc, nd, d0, ncol, mode
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(out) :: taus(no*no*ncol)
      integer :: t, dest

      do concurrent(t=1:no*no*nvc*nd) local(dest)
         dest = gather_dest(TRC_DC_INDEX(t), no, nv, nvc, d0)
         if (dest > 0) taus(dest) = gather_value(TRC_DC_INDEX(t), no, nv, nvc, d0, mode, tau)
      end do
   end subroutine kernel_half_gather

   subroutine kernel_half_sym(no, nv, nacc, h)
      !! h(i,j,a,b) = nacc(i,j,a,b) + nacc(j,i,b,a)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: nacc(no*no*nv*nv)
      real(dp), intent(out) :: h(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         h(t) = nacc(t) + nacc(swap_pair(TRC_DC_INDEX(t), no, nv))
      end do
   end subroutine kernel_half_sym

   subroutine kernel_half_general(no, nv, nprime, h)
      !! h(i,j,a,b) += nprime(i,j,b,a)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: nprime(no*no*nv*nv)
      real(dp), intent(inout) :: h(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         h(t) = h(t) + nprime(swap_ab(TRC_DC_INDEX(t), no, nv))
      end do
   end subroutine kernel_half_general

   subroutine kernel_copy(n, x, y)
      !! y = x
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: y(n)
      integer :: t

      do concurrent(t=1:n)
         y(t) = x(t)
      end do
   end subroutine kernel_copy

   subroutine kernel_pack_ovov(no, nv, ovov, ovovp)
      !! ovovp(k,l,c,d) = ovov(k,c,l,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp), intent(out) :: ovovp(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         ovovp(t) = ovovp_item(TRC_DC_INDEX(t), no, nv, ovov)
      end do
   end subroutine kernel_pack_ovov

   subroutine kernel_sub_x(n, x, h)
      !! h -= x
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(inout) :: h(n)
      integer :: t

      do concurrent(t=1:n)
         h(t) = h(t) - x(t)
      end do
   end subroutine kernel_sub_x

   subroutine kernel_sub_xsym(no, nv, x, h)
      !! h(i,j,a,b) -= (x(i,j,a,b) + x(j,i,b,a))/2: the correction of the symmetrised tau, which is what the symmetric
      !! path's passes contract
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: x(no*no*nv*nv)
      real(dp), intent(inout) :: h(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         h(t) = h(t) - 0.5_dp*(x(t) + x(swap_pair(TRC_DC_INDEX(t), no, nv)))
      end do
   end subroutine kernel_sub_xsym

   subroutine kernel_pack_y2(no, nv, np, y, y2)
      !! y2(ij,a,d,P) = y(a,d,ij,P)
      integer, intent(in) :: no, nv, np
      real(dp), intent(in) :: y(nv*nv*no*no*np)
      real(dp), intent(out) :: y2(no*no*nv*nv*np)
      integer :: t

      do concurrent(t=1:no*no*nv*nv*np)
         y2(t) = y2_item(TRC_DC_INDEX(t), no, nv, np, y)
      end do
   end subroutine kernel_pack_y2

   subroutine kernel_pack_taut_swapped(no, nv, tau, taut)
      !! taut(x,y,ij) = tau(i,j,y,x), ij = i + o (j-1): the pair in the order of the dressing's half-1 pack
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(out) :: taut(nv*nv*no*no)
      integer :: t

      do concurrent(t=1:nv*nv*no*no)
         taut(t) = taut_swapped_item(TRC_DC_INDEX(t), no, nv, tau)
      end do
   end subroutine kernel_pack_taut_swapped

   subroutine kernel_sub_z1(no, nv, v0, nvb, cmat, h)
      !! h(i,j,a,b) -= cmat(ij; b, al), a = v0 + al - 1
      integer, intent(in) :: no, nv, v0, nvb
      real(dp), intent(in) :: cmat(no*no*nv*nvb)
      real(dp), intent(inout) :: h(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nvb)
         h(sub_z1_dest(TRC_DC_INDEX(t), no, nv, v0)) = h(sub_z1_dest(TRC_DC_INDEX(t), no, nv, v0)) - cmat(t)
      end do
   end subroutine kernel_sub_z1

   subroutine kernel_sub_z2(no, nv, v0, nvb, cmat, h)
      !! h(i,j,a,b) -= cmat(ij; a, bl), b = v0 + bl - 1
      integer, intent(in) :: no, nv, v0, nvb
      real(dp), intent(in) :: cmat(no*no*nv*nvb)
      real(dp), intent(inout) :: h(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nvb)
         h(sub_z2_dest(TRC_DC_INDEX(t), no, nv, v0)) = h(sub_z2_dest(TRC_DC_INDEX(t), no, nv, v0)) - cmat(t)
      end do
   end subroutine kernel_sub_z2

   pure subroutine split4(t, n1, n2, n3, i1, i2, i3, i4)
      !! ---------------------------------------------------------------------
      !! the item routines: pure, `routine seq`, explicit-shape dummies, no
      !! arrays of their own. A flat 1-based index t is the position in the
      !! array that is written (or, for the scatters, read), first index fastest.
      !! A call in a `do concurrent` body passes the loop index through the TRC_DC_INDEX macro:
      !! see inc/trc_dc_index.inc (an ifx miscompilation of the plain variable).
      !! ---------------------------------------------------------------------
      !! t - 1 = (i1-1) + n1 ((i2-1) + n2 ((i3-1) + n3 (i4-1)))
      !$acc routine seq
      integer, intent(in) :: t, n1, n2, n3
      integer, intent(out) :: i1, i2, i3, i4
      integer :: q

      q = t - 1
      i1 = mod(q, n1) + 1
      q = q/n1
      i2 = mod(q, n2) + 1
      q = q/n2
      i3 = mod(q, n3) + 1
      i4 = q/n3 + 1
   end subroutine split4

   pure function taut_item(t, no, nv, tau) result(v)
      !! taut(c,d,ij) = tau(i,j,c,d), position t of (nv, nv, no, no)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp) :: v
      integer :: c, d, i, j

      call split4(t, nv, nv, no, c, d, i, j)
      v = tau(i, j, c, d)
   end function taut_item

   pure function half_col(c, d, d0, nv) result(col)
      !! Column of (c,d), c >= d, in a half-form pass from d0: the blocks of d0..d-1 (nv - d'' + 1 columns each), then c - d + 1
      !$acc routine seq
      integer, intent(in) :: c, d, d0, nv
      integer :: col
      integer :: k

      k = d - d0
      col = k*(nv - d0 + 1) - (k*(k - 1))/2 + (c - d) + 1
   end function half_col

   pure function resort_dest(t, nv, nvc, d0) result(dest)
      !! Position in w2(nv, nv, ncol) of cbuf position t of (nv[a], nvc[c - d0 + 1], nv[b], nd[d - d0 + 1]); 0 for
      !! c < d, negated for c = d (the columns of weight 1/2)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: nv, nvc, d0
      integer :: dest
      integer :: a, cp, b, dq

      call split4(t, nv, nvc, nv, a, cp, b, dq)
      dest = 0
      if (cp < dq) return
      dest = a + nv*(b - 1) + nv*nv*(half_col(d0 + cp - 1, d0 + dq - 1, d0, nv) - 1)
      if (cp == dq) dest = -dest
   end function resort_dest

   pure function gather_dest(t, no, nv, nvc, d0) result(dest)
      !! Position in taus(no*no, ncol) of position t of (no[i], no[j], nvc[c - d0 + 1], nd[d - d0 + 1]); 0 for c < d
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv, nvc, d0
      integer :: dest
      integer :: i, j, cp, dq

      call split4(t, no, no, nvc, i, j, cp, dq)
      dest = 0
      if (cp < dq) return
      dest = i + no*(j - 1) + no*no*(half_col(d0 + cp - 1, d0 + dq - 1, d0, nv) - 1)
   end function gather_dest

   pure function gather_value(t, no, nv, nvc, d0, mode, tau) result(v)
      !! The tau value `kernel_half_gather` stores for position t of (no[i], no[j], nvc[c - d0 + 1], nd[d - d0 + 1])
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv, nvc, d0, mode
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp) :: v
      integer :: i, j, cp, dq, c, d

      call split4(t, no, no, nvc, i, j, cp, dq)
      c = d0 + cp - 1
      d = d0 + dq - 1
      if (mode == HALF_TAU_CD) then
         v = tau(i, j, c, d)
      else if (mode == HALF_TAU_DC) then
         v = tau(i, j, d, c)
      else
         v = 0.5_dp*(tau(i, j, c, d) + tau(j, i, d, c))
      end if
   end function gather_value

   pure function ovovp_item(t, no, nv, ovov) result(v)
      !! ovov(k,c,l,d) at position t of (no[k], no[l], nv[c], nv[d])
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ovov(no, nv, no, nv)
      real(dp) :: v
      integer :: k, l, c, d

      call split4(t, no, no, nv, k, l, c, d)
      v = ovov(k, c, l, d)
   end function ovovp_item

   pure function swap_pair(t, no, nv) result(s)
      !! Position of (j,i,b,a) for position t of (i,j,a,b) in (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      integer :: s
      integer :: i, j, a, b

      call split4(t, no, no, nv, i, j, a, b)
      s = j + no*(i - 1) + no*no*(b - 1) + no*no*nv*(a - 1)
   end function swap_pair

   pure function swap_ab(t, no, nv) result(s)
      !! Position of (i,j,b,a) for position t of (i,j,a,b) in (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      integer :: s
      integer :: ij, a, b, unused

      call split4(t, no*no, nv, nv, ij, a, b, unused)
      s = ij + no*no*(b - 1) + no*no*nv*(a - 1)
   end function swap_ab

   pure function y2_item(t, no, nv, np, y) result(v)
      !! y2(ij,a,d,P) = y(a,d,ij,P), position t of (no*no, nv, nv, np)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv, np
      real(dp), intent(in) :: y(nv*nv*no*no*np)
      real(dp) :: v
      integer :: ij, a, d, p

      call split4(t, no*no, nv, nv, ij, a, d, p)
      v = y(a + nv*(d - 1) + nv*nv*(ij - 1) + nv*nv*no*no*(p - 1))
   end function y2_item

   pure function taut_swapped_item(t, no, nv, tau) result(v)
      !! taut(x,y,ij) = tau(i,j,y,x), position t of (nv, nv, no, no)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp) :: v
      integer :: x, y, i, j

      call split4(t, nv, nv, no, x, y, i, j)
      v = tau(i, j, y, x)
   end function taut_swapped_item

   pure function sub_z1_dest(t, no, nv, v0) result(dest)
      !! Position in h(i,j,a,b) of cmat position t of (no*no, nv[b], nvb[a - v0 + 1])
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv, v0
      integer :: dest
      integer :: ij, b, al, unused

      call split4(t, no*no, nv, 1, ij, b, unused, al)
      dest = ij + no*no*(v0 + al - 2) + no*no*nv*(b - 1)
   end function sub_z1_dest

   pure function sub_z2_dest(t, no, nv, v0) result(dest)
      !! Position in h(i,j,a,b) of cmat position t of (no*no, nv[a], nvb[b - v0 + 1])
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv, v0
      integer :: dest

      dest = t + no*no*nv*(v0 - 1)
   end function sub_z2_dest

end module trc_rccsd_ladder
