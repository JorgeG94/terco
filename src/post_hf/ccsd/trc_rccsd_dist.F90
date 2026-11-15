module trc_rccsd_dist
   !!
   !! The iterations of an RI-CCSD run over the ranks of a communicator, one rank per device: what every rank does at the
   !! top of an iteration and what a rank other than the first does for the rest of it. The first rank runs the
   !! iteration itself (trc_rccsd_driver, trc_rccsd_update); module trc_rccsd_driver (RANKS) has the whole scheme.
   !!
   !! THE TOP OF AN ITERATION (`rccsd_dist_start`, collective, one call per iteration on every rank):
   !!   1. the control word, one integer from the first rank: 0 to stop, else 2 it + (1 if tau has the pair symmetry the
   !!      half ladder's one-contraction path needs, else 0). The ladder's path must be the same on every rank, so the
   !!      flag travels with the word instead of being decided twice;
   !!   2. when it is not 0, the amplitudes, device to device (trc_rccsd_xfer: point to point from the first rank to each
   !!      other, the unique half of t2 when tau has the pair symmetry, through the host or MPI device buffers by the
   !!      run-time switch of trc_mpi_device). The row RCCSD_ST_BCAST is this whole step, its bytes counted on the first
   !!      rank as the sum over the receivers.
   !! A follower that gets 0 leaves its loop. The first rank sends 0 when it leaves the loop for any reason of its
   !! own (converged, the limit, a NaN, a refusal in a stage after the ladder's sum); a failure the ranks agreed on inside
   !! the iteration (trc_rccsd_update, `ladder_over_ranks`) stopped them all already and no word is due.
   use pic_mpi_lib, only: comm_t
   use trc_linalg, only: trc_linalg_t
   use trc_step_table, only: step_table_t
   use pic_timer, only: timer_type
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_update, only: trc_rccsd_update_t
   use trc_rccsd_steps, only: RCCSD_ST_BCAST, rccsd_steps_open, rccsd_steps_close
   use trc_rccsd_mpi, only: rccsd_mpi_control
   use trc_rccsd_xfer, only: rccsd_xfer_amps_down, rccsd_xfer_amps_bytes
   use pic_types, only: int64
   implicit none
   private

   public :: rccsd_dist_start, rccsd_dist_follow, rccsd_dist_word

contains

   pure function rccsd_dist_word(iteration, tau_symmetric) result(word)
      !! The control word of iteration `iteration` (>= 1); 0 is the stop word
      integer, intent(in) :: iteration
      logical, intent(in) :: tau_symmetric
      integer :: word

      word = 2*iteration + merge(1, 0, tau_symmetric)
   end function rccsd_dist_word

   function rccsd_dist_start(comm, lead, word, upd, amps, la, stab) result(received)
      !! Collective. The top of an iteration on every rank; `word` is read on the first rank only. Returns the word and
      !! sets `upd%tau_symmetric` from it on a follower.
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: lead
      integer, intent(in) :: word
      type(trc_rccsd_update_t), intent(inout) :: upd
      type(trc_rccsd_amps_t), intent(inout) :: amps
      type(trc_linalg_t), intent(in) :: la
      type(step_table_t), intent(inout) :: stab
      integer :: received

      type(timer_type) :: tc
      logical :: active(comm%size())

      received = rccsd_mpi_control(comm, word)
      if (received == 0) return
      if (.not. lead) upd%tau_symmetric = mod(received, 2) == 1
      ! the ranks that get the amplitudes are decided from the one share table every rank holds: a follower whose range
      ! is empty takes the control word and the agreement and nothing else
      active = upd%range_hi >= upd%range_lo
      active(1) = .true.
      if (.not. active(comm%rank() + 1)) return
      call rccsd_steps_open(stab, la, RCCSD_ST_BCAST, tc)
      call rccsd_xfer_amps_down(comm, lead, upd%xf, amps, upd%tau_symmetric, active)
      call rccsd_steps_close(stab, la, RCCSD_ST_BCAST, tc)
      if (lead) call stab%add_bytes(RCCSD_ST_BCAST, int(count(active) - 1, int64)* &
                                    rccsd_xfer_amps_bytes(upd%no, upd%nv, upd%tau_symmetric))
   end function rccsd_dist_start

   subroutine rccsd_dist_follow(comm, upd, ints, amps, la, stab)
      !! A rank other than the first, from its first control word to the end of the iterations: receive the amplitudes,
      !! run the share, until the stop word or an agreed failure. `upd` was initialised as a follower and given its share.
      type(comm_t), intent(in) :: comm
      type(trc_rccsd_update_t), intent(inout) :: upd
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(inout) :: amps
      type(trc_linalg_t), intent(in) :: la
      type(step_table_t), intent(inout) :: stab

      integer :: word

      do
         word = rccsd_dist_start(comm, .false., 0, upd, amps, la, stab)
         if (word == 0) exit
         call upd%follow(ints, amps, la, stab)
         if (upd%agreed_stop) exit
      end do
   end subroutine rccsd_dist_follow

end module trc_rccsd_dist
