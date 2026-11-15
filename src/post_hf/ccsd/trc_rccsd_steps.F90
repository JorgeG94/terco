module trc_rccsd_steps
   !!
   !! The steps of an RI-CCSD run as the step table of the driver counts them, and
   !! the helpers every stage shares to feed it: the indices and names, the flop
   !! count of one (batched) GEMM, and `rccsd_steps_open` / `rccsd_steps_close`, which bracket a
   !! stage the way RI-MP2 brackets its own steps (a device sync, the clock and an
   !! NVTX range at both ends, see trc_step_table) and do nothing when the caller
   !! passed no table, so a stage can still be driven on its own.
   !!
   !! OVER RANKS the table has three more rows, RCCSD_ST_BCAST, _WAIT and _REDUCE (RCCSD_MULTI_STEP_NAME): they
   !! are the cost of the distributed ladder, 0 in a one-rank run, and a one-rank run does not print them.
   !!
   !! WHAT A ROW MEANS. The rows from `RCCSD_ST_METRIC` to `RCCSD_ST_BLOCKS` are
   !! inside the `RCCSD_ST_INTEGRALS` row (indented in the print): the integrals row
   !! is the whole `ints%build`, and the rows under it are the pieces that build
   !! times itself, so they do not add up to it (the orbital transforms to the
   !! Cartesian space, the allocations and the copy of the blocks back to the host
   !! are in the first only). Every other row is a stage of one iteration, summed
   !! over the iterations.
   !!
   !! WHERE THE FLOPS ARE COUNTED. A row's work is the sum of 2 m n k over the GEMMs
   !! (and the triangular solves and the Cholesky factorisation) of its stage, with
   !! the count of a batched call multiplied in. Kernel-only work (packs, scatters,
   !! the denominators, the per-element base kernels of W, Lvv, the small T1 terms,
   !! the energy, DIIS) has time but no count, so a row that mixes GEMMs and kernels
   !! (every stage from the intermediates to the T2 rings) shows a rate below the
   !! card's GEMM rate by the share of its time the kernels take; the ladder and
   !! the T2 rings are almost pure GEMM. The integral rows `3-index integrals` and
   !! `metric` are integral kernels plus the GEMMs of the Cartesian-to-spherical
   !! conversion, which are not counted (a few per cent of the build's GEMM work
   !! on water with cc-pVDZ-RIFIT). The counts are those of the stage
   !! modules' own functions (`rccsd_inter_flops`, `rccsd_w_flops`,
   !! `rccsd_ladder_flops`, `rccsd_ladder_dressing_flops`, `rccsd_update_*_flops`), and
   !! trc_rccsd_ints sums the calls it makes.
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_linalg, only: trc_linalg_t
   use trc_step_table, only: step_table_t
   implicit none
   private

   integer, parameter, public :: RCCSD_ST_INTEGRALS = 1  !! the whole `ints%build`
   integer, parameter, public :: RCCSD_ST_METRIC = 2     !! (P|Q) and its conversion to the spherical auxiliary basis
   integer, parameter, public :: RCCSD_ST_CHOLESKY = 3   !! the Cholesky factor of the metric
   integer, parameter, public :: RCCSD_ST_3C = 4         !! (mu nu|P), with its conversion to the spherical space
   integer, parameter, public :: RCCSD_ST_HALF = 5       !! the half transform and the three MO transforms per auxiliary block
   integer, parameter, public :: RCCSD_ST_TRSM = 6       !! B = X L^-T
   integer, parameter, public :: RCCSD_ST_BLOCKS = 7     !! the four blocks (oo|oo), (oo|ov), (oo|vv), (ov|ov) from B
   integer, parameter, public :: RCCSD_ST_INTER = 8      !! tau and the one-particle intermediates
   integer, parameter, public :: RCCSD_ST_W = 9          !! Woooo, Wvoov, Wvovo
   integer, parameter, public :: RCCSD_ST_LADDER = 10    !! the particle ladder; its name says the form that ran
   integer, parameter, public :: RCCSD_ST_DRESSING = 11  !! the t1 dressing of the ladder
   integer, parameter, public :: RCCSD_ST_T1 = 12        !! the T1 equation
   integer, parameter, public :: RCCSD_ST_T2 = 13        !! the T2 equation without its three ring terms
   integer, parameter, public :: RCCSD_ST_RING = 14      !! the three ring terms of the T2 equation
   integer, parameter, public :: RCCSD_ST_DIIS = 15      !! the step measure, the history and the extrapolation
   integer, parameter, public :: RCCSD_ST_ENERGY = 16    !! the energy of the new amplitudes
   integer, parameter, public :: RCCSD_NSTEP = 16
   integer, parameter, public :: RCCSD_ST_BCAST = 17     !! several ranks only: the amplitudes to the other ranks, once an iteration
   integer, parameter, public :: RCCSD_ST_WAIT = 18      !! several ranks only: the ranks arriving at the end of the ladder's passes
   integer, parameter, public :: RCCSD_ST_REDUCE = 19    !! several ranks only: the partial ladders to the first rank
   integer, parameter, public :: RCCSD_NSTEP_MULTI = 19  !! rows of the table of a run over several ranks
   integer, parameter, public :: RCCSD_NAME_LEN = 24     !! the step table's name width

   character(len=RCCSD_NAME_LEN), parameter, public :: RCCSD_STEP_NAME(RCCSD_NSTEP) = [character(len=RCCSD_NAME_LEN) :: &
      "integrals", "  metric", "  cholesky", "  3-index integrals", "  half transform", "  B = X L^-T", "  four blocks", &
      "tau, intermediates", "W intermediates", "ladder", "ladder dressing", "T1 equation", "T2 equation", &
      "T2 ring terms", "DIIS", "energy"]

   character(len=RCCSD_NAME_LEN), parameter, public :: RCCSD_MULTI_STEP_NAME(RCCSD_NSTEP_MULTI) = &
      [character(len=RCCSD_NAME_LEN) :: RCCSD_STEP_NAME, "amplitudes to ranks", "ladder wait (idle)", "ladder to first rank"]
   !! the rows of a run over several ranks: the single-rank ones, then the two transfers of an iteration (the amplitudes
   !! from the first rank to each other, the partial ladders from each other rank to the first; point to point, packed
   !! when pair-symmetric, bytes counted on the first rank) and the idle time at the ladder's barrier, which is where a
   !! share too large or too small for its rank shows

   public :: rccsd_steps_gemm_flops, rccsd_steps_trsm_flops, rccsd_steps_potrf_flops, rccsd_steps_open, rccsd_steps_close

contains

   pure function rccsd_steps_gemm_flops(m, n, k, count) result(flops)
      !! 2 m n k for each of `count` GEMMs of one shape; 0 when any extent is below one
      integer, intent(in) :: m, n, k
      integer, intent(in), optional :: count
      integer(int64) :: flops

      integer :: c

      c = 1
      if (present(count)) c = count
      flops = 0_int64
      if (m < 1 .or. n < 1 .or. k < 1 .or. c < 1) return
      flops = 2_int64*m*n*k*c
   end function rccsd_steps_gemm_flops

   pure function rccsd_steps_trsm_flops(m, n) result(flops)
      !! A triangular solve of an m by n block against an n by n factor: m n^2, as RI-MP2 counts it
      integer, intent(in) :: m, n
      integer(int64) :: flops

      flops = int(m, int64)*n*n
   end function rccsd_steps_trsm_flops

   pure function rccsd_steps_potrf_flops(n) result(flops)
      !! The Cholesky factorisation of an n by n matrix: n^3/3, as RI-MP2 counts it
      integer, intent(in) :: n
      integer(int64) :: flops

      flops = int(n, int64)**3/3
   end function rccsd_steps_potrf_flops

   subroutine rccsd_steps_open(stab, la, k, tc)
      !! Opens step `k` of `stab` (device sync, clock, NVTX range); nothing when there is no table
      type(step_table_t), intent(in), optional :: stab
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: k
      type(timer_type), intent(out) :: tc

      if (present(stab)) call stab%step_on(la, k, tc)
   end subroutine rccsd_steps_open

   subroutine rccsd_steps_close(stab, la, k, tc, flops)
      !! Closes step `k` after the same sync and adds the `flops` its stage executed; nothing when there is no table
      type(step_table_t), intent(inout), optional :: stab
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: k
      type(timer_type), intent(in) :: tc
      integer(int64), intent(in), optional :: flops

      if (.not. present(stab)) return
      call stab%step_off(la, k, tc)
      if (present(flops)) call stab%add_flop(k, flops)
   end subroutine rccsd_steps_close

end module trc_rccsd_steps
