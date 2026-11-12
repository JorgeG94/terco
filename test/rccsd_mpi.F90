!
! RI-CCSD(T) over MPI ranks, one rank per device (src/post_hf/ccsd/trc_rccsd_driver.F90 and trc_rccsd_triples.F90,
! RANKS). Water and a C1 hydrogen peroxide, cc-pVDZ with cc-pVDZ-RIFIT (check_rccsd_triples' molecules and pyscf
! references, read through rccsd_ref_energies). The SCF runs over the communicator, as the CLI and the C API run it
! (test/scf_mpi.F90). The program is generic in the rank count; ctest runs it on 1 and 2, on a host MPI build too (the
! (T) routine is the occupied-driven one, which runs on the host).
!
! No argument: rank 0 runs the CCSD(T) alone first (no communicator, the single-rank code path), then every rank runs
! it over the communicator. The CCSD iterations are shared too (the ladder of every iteration is cut over the ranks), so
! the ranks accumulate the ladder in another order than one rank does: the numbers agree to round-off, NOT bitwise.
! Required, on every rank:
!   - converged, (T) done, no error;
!   - E_CCSD within TOL_ONE_RANK_E of the single-rank value, and the same number of iterations (rank 0 computes the
!     single-rank run alone, the others receive its numbers). Both stop on |dE| < 1e-12 and max |step| < 1e-10 at the
!     same iteration, so they differ by the round-off of the ladder sum carried through the contraction of the
!     iteration. Measured on water and h2o2, host and GPU builds, 2 to 4 ranks: 1.7e-16 at most. The bound is
!     TOL_ONE_RANK_E = 1e-14, two orders over that; a ladder lost, doubled or built from stale amplitudes moves the
!     energy by 1e-4 or more;
!   - EVERY element of t1 and t2, read back from the rank's device, within TOL_ONE_RANK_AMP = 1e-14 of rank 0's
!     single-rank amplitudes (the converged iterates; measured 4.2e-17 at most);
!   - E(T) within TOL_TRIPLES_ONE_RANK = 1e-14 of the single-rank value (measured 3.5e-18 at most; the amplitudes
!     differ by the above) -- and within TOL_TRIPLES_SAME = 1e-16, the original bound, where the amplitudes are
!     bitwise the single-rank ones (a communicator of one rank);
!   - the energy (and the largest step) after ONE and after TWO iterations (max_iterations 1 and 2, not converged)
!     within 1e-14 of the one-rank run's: a follower working from amplitudes one iteration old, or none, changes these
!     by far more, while the final energy of a converged run would absorb part of it (`prefix_case`);
!   - the energies against pyscf;
!   - E_CCSD, E(T), their sum and the iteration count equal over the ranks (MIN equals MAX), and the amplitudes bitwise
!     equal among the ranks (their sums are: the final broadcast copies the first rank's).
! Water is run once more with the sign of one occupied and one virtual orbital column flipped on every rank but the
! first (the GUARD case: a caller that hands ranks different orbitals). The same energies and amplitudes are required;
! without the orbital broadcast of trc_rccsd_solve the (T) of the other ranks is wrong and this case fails.
!
! Argument "fail": the failure paths, which must end on every rank without a rank left waiting (ctest puts a timeout
! on it), and with the TEXT of the refusing rank's reason on every rank: a CCSD that does not converge (every rank sees
! it; the followers stop with the first rank and (T) is refused on all of them); a refusal on ONE rank (the last, and
! rank 0) that makes every rank refuse; a refusal of one rank's LADDER share (a negative d_batch on that rank only);
! an explicit assembled ladder, which cannot be cut over ranks (every rank asks, then one rank); a rank whose
! ladder range is refused INSIDE an iteration (a range beyond nv, on the last rank and then on rank 0, through the
! public `set_share`/`step`/`follow`): every rank must stop with that rank's reason, and the same objects then run a
! healthy iteration on the converged amplitudes and reproduce them (the transfers of an iteration between the ranks,
! checked on their own), then the same with the last rank's range EMPTY (a hand-set table: the rule never leaves a
! follower empty for these molecules): that rank must be sent no amplitudes, build no tau and still take the control word
! and the agreement, with the first rank's send loop and the follower deciding from the same table; and a step table
! without the transfer rows, refused through the agreement. WHICH PATH THIS COVERS: the agreed stop inside an
! iteration is exercised through the update type's own API (`set_share`, `step`, `follow`), not through
! `trc_rccsd_solve`: the cut a solve uses is the pure rule (always valid), the integrals are checked and agreed before
! the loop, and nothing in a solve's inputs makes one rank's share fail later, so a solve cannot reach that branch
! without a test-only hook, which production code does not have. The solve-level stops that ARE tested are the
! agreed ones before the loop (below) and the control word's stop (not converged); a (T)
! routine that one rank asks for differently from the others; the host reference, which has no work share; one rank's
! state marked not converged; one rank's caller-set work share; each followed by a correct (T) on the same state.
!
! A collective-MPI program: no internal procedures (NVHPC 26.5 + HPC-X). The cases are module procedures.
module rccsd_mpi_cases
   use pic_types, only: dp
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_MAX, MPI_MIN
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_error, only: error_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rccsd_driver, only: trc_rccsd_options_t, trc_rccsd_state_t, trc_rccsd_solve
   use trc_rccsd_ladder, only: RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_HALF, trc_rccsd_ladder_options_t
   use trc_rccsd_update, only: trc_rccsd_update_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_share, only: rccsd_share_ranges
   use trc_step_table, only: step_table_t
   use trc_rccsd_steps, only: RCCSD_STEP_NAME
   use trc_rccsd_triples, only: trc_rccsd_triples_options_t, rccsd_triples_finish, RCCSD_TRIPLES_DEVICE_IJK, &
                                RCCSD_TRIPLES_DEVICE_ABC
   use rccsd_test_support, only: NMOL, MOL_LABEL, MOL_REF, MOL_NFROZEN, load_molecule
   use rccsd_ref_energies, only: rccsd_refs_t, read_rccsd_refs
   implicit none
   private
   public :: success_case, failure_cases, skew_start_case, prefix_case

   ! Mirrors of constants defined in the single-rank checks (test programs cannot share them without a module of their
   ! own): TOL_CORR, TOL_T, TOL_TOTAL are check_capi_rccsd's; TOL_TRIPLES is check_rccsd_triples' TOL_IJK_TILE.
   real(dp), parameter :: TOL_CORR = 1.0e-10_dp   !! E_CCSD against pyscf
   real(dp), parameter :: TOL_T = 3.0e-12_dp      !! E(T) against pyscf
   real(dp), parameter :: TOL_TOTAL = 3.0e-11_dp  !! the CCSD(T) total against pyscf
   real(dp), parameter :: TOL_ONE_RANK_E = 1.0e-14_dp    !! E_CCSD over the ranks against one rank; measured 1.7e-16 (header)
   real(dp), parameter :: TOL_ONE_RANK_AMP = 1.0e-14_dp  !! the largest amplitude difference; measured 4.2e-17 (header)
   real(dp), parameter :: TOL_TRIPLES_ONE_RANK = 1.0e-14_dp  !! E(T) over the ranks against one rank; measured 3.5e-18
   real(dp), parameter :: TOL_TRIPLES_SAME = 1.0e-16_dp  !! E(T) from bitwise the same amplitudes (one rank in the communicator)

contains

   function same_on_ranks(comm, x) result(same)
      !! True when `x` is the same number on every rank of `comm` (MIN equals MAX); collective
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: x
      logical :: same

      real(dp) :: lo, hi

      lo = x
      hi = x
      call allreduce(comm, lo, op=MPI_MIN)
      call allreduce(comm, hi, op=MPI_MAX)
      same = lo == hi
   end function same_on_ranks

   function all_ranks(comm, flag) result(all_true)
      !! True when `flag` is true on every rank of `comm`; collective
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: flag
      logical :: all_true

      real(dp) :: x

      x = merge(1.0_dp, 0.0_dp, flag)
      call allreduce(comm, x, op=MPI_MIN)
      all_true = x > 0.5_dp
   end function all_ranks

   subroutine scf_over_comm(comm, bas, nocc, res)
      !! RHF to energy 1e-13 and commutator 1e-11 (rccsd_test_support's run_scf), over the communicator as the CLI runs it
      type(comm_t), intent(in) :: comm
      type(trc_basis_t), intent(in) :: bas
      integer, intent(in) :: nocc
      type(trc_scf_result_t), intent(out) :: res

      type(trc_scf_options_t) :: opts

      opts = trc_scf_options_t()
      opts%conv_energy = 1.0e-13_dp
      opts%conv_density = 1.0e-11_dp
      opts%conv_diis = 1.0e-11_dp
      call trc_scf_run(bas, nocc, nocc, opts, res, comm=comm)
   end subroutine scf_over_comm

   subroutine converged_options(comm, nocc, nfz, opts)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: nocc, nfz
      type(trc_rccsd_options_t), intent(out) :: opts

      opts%nocc = nocc
      opts%nfrozen = nfz
      opts%conv_energy = 1.0e-12_dp
      opts%conv_amplitude = 1.0e-10_dp
      opts%max_iterations = 200
      opts%comm = comm
      opts%has_comm = .true.
   end subroutine converged_options

   subroutine amps_gap(comm, state, one_t1, one_t2, gap)
      !! The largest |difference| of any element of this rank's amplitudes (pulled back from its device) from the
      !! single-rank ones, which rank 0 broadcasts (`one_t1`, `one_t2` hold rank 0's on entry and everywhere on exit)
      type(comm_t), intent(in) :: comm
      type(trc_rccsd_state_t), intent(inout) :: state
      real(dp), intent(inout) :: one_t1(:, :), one_t2(:, :, :, :)
      real(dp), intent(out) :: gap

      real(dp), allocatable :: flat(:)

      flat = reshape(one_t1, [size(one_t1)])
      call bcast(comm, flat, size(flat), 0)
      one_t1 = reshape(flat, shape(one_t1))
      flat = reshape(one_t2, [size(one_t2)])
      call bcast(comm, flat, size(flat), 0)
      one_t2 = reshape(flat, shape(one_t2))
      call state%amps%pull()
      gap = max(maxval(abs(state%amps%t1 - one_t1)), maxval(abs(state%amps%t2 - one_t2)))
   end subroutine amps_gap

   subroutine run_over_ranks(comm, label, bas, aux, pl, cmo, eps, opts, topts, la, single, one_t1, one_t2, refs, good)
      !! One solve and one (T) over the communicator from the orbitals `cmo`, `eps` of THIS rank, against the single-rank
      !! numbers `single` = (E_CCSD, E(T), iterations) and amplitudes of rank 0 and against pyscf; `good` is the same on every rank
      type(comm_t), intent(in) :: comm
      character(len=*), intent(in) :: label
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: cmo(:, :), eps(:)
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_rccsd_triples_options_t), intent(in) :: topts
      type(trc_linalg_t), intent(inout) :: la
      real(dp), intent(in) :: single(3)
      real(dp), intent(inout) :: one_t1(:, :), one_t2(:, :, :, :)
      type(rccsd_refs_t), intent(in) :: refs
      logical, intent(out) :: good

      type(trc_rccsd_state_t) :: state
      type(error_t) :: status
      logical :: uniform
      real(dp) :: gap, e_total, e_gap
      integer :: rank

      rank = comm%rank()
      call trc_rccsd_solve(bas, aux, pl, cmo, eps, opts, la, state)
      call rccsd_triples_finish(state, status, topts)
      good = .not. state%result%refused .and. state%result%converged .and. .not. status%has_error() .and. &
             state%result%has_triples
      if (status%has_error()) print "(a,i0,a)", "rccsd_mpi: rank ", rank, ": (T) refused: "//status%get_message()
      call amps_gap(comm, state, one_t1, one_t2, gap)
      e_total = refs%e_scf + state%result%e_ccsd_t
      uniform = same_on_ranks(comm, state%result%e_ccsd) .and. same_on_ranks(comm, state%result%e_t) .and. &
                same_on_ranks(comm, state%result%e_ccsd_t) .and. same_on_ranks(comm, real(state%result%iterations, dp))
      ! the amplitudes bitwise among the ranks: their sums are the same number (the sum runs in the same order on each)
      uniform = uniform .and. same_on_ranks(comm, sum(state%amps%t2)) .and. same_on_ranks(comm, sum(abs(state%amps%t1)))
      e_gap = abs(state%result%e_ccsd - single(1))
      if (rank == 0) then
         print "(a,i0,a,f20.12,a,i0,a,es9.2,a,es9.2,a,es9.2,a,es9.2)", "rccsd_mpi: "//label//"ranks ", comm%size(), &
            "  E_CCSD ", state%result%e_ccsd, "  iterations ", state%result%iterations, "  |E_CCSD - one rank| ", e_gap, &
            "  max |amplitude - one rank| ", gap, "  |E(T) - one rank| ", abs(state%result%e_t - single(2)), &
            "  |E(T) - pyscf| ", abs(state%result%e_t - refs%e_t)
      end if
      good = good .and. e_gap <= TOL_ONE_RANK_E .and. gap <= TOL_ONE_RANK_AMP .and. &
             nint(single(3)) == state%result%iterations .and. &
             abs(state%result%e_t - single(2)) <= merge(TOL_TRIPLES_SAME, TOL_TRIPLES_ONE_RANK, comm%size() == 1) .and. &
             abs(state%result%e_ccsd - refs%e_ccsd_corr) < TOL_CORR .and. abs(state%result%e_t - refs%e_t) < TOL_T .and. &
             abs(e_total - refs%e_ccsd_t) < TOL_TOTAL .and. abs(state%result%e_ccsd_t - state%result%e_ccsd - &
                                                                 state%result%e_t) < 1.0e-15_dp
      good = all_ranks(comm, good .and. uniform)
      call state%release()
   end subroutine run_over_ranks

   subroutine success_case(comm, im, ok)
      !! Molecule `im`, no argument mode (program header); `ok` is the same on every rank
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: im
      logical, intent(inout) :: ok

      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: la
      type(trc_rccsd_state_t) :: solo
      type(trc_rccsd_options_t) :: opts
      type(trc_rccsd_triples_options_t) :: topts
      type(error_t) :: status
      type(rccsd_refs_t) :: refs
      integer, allocatable :: zint(:)
      real(dp), allocatable :: cmo(:, :), one_t1(:, :), one_t2(:, :, :, :)
      integer :: nocc, rank, no, nv
      logical :: supported, good
      real(dp) :: single(3)
      character(len=:), allocatable :: label

      rank = comm%rank()
      label = trim(MOL_LABEL(im))//": "
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) then
         if (rank == 0) print "(a)", "rccsd_mpi: SKIP, this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
         return
      end if
      call scf_over_comm(comm, bas, nocc, res)
      call la%init(2)
      call read_rccsd_refs(trim(MOL_REF(im)), refs)
      call converged_options(comm, nocc, MOL_NFROZEN(im), opts)
      topts%use_device = .true.
      topts%device_routine = RCCSD_TRIPLES_DEVICE_IJK

      ! the single-rank run, rank 0 only: no communicator, no MPI call
      single = 0.0_dp
      no = nocc - MOL_NFROZEN(im)
      nv = bas%nao_sph - nocc
      allocate (one_t1(no, nv), one_t2(no, no, nv, nv))
      one_t1 = 0.0_dp
      one_t2 = 0.0_dp
      if (rank == 0) then
         opts%has_comm = .false.
         call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, solo)
         call rccsd_triples_finish(solo, status, topts)
         if (status%has_error() .or. .not. solo%result%converged) error stop "rccsd_mpi: the single-rank run failed: "// &
            status%get_message()
         call solo%amps%pull()
         single = [solo%result%e_ccsd, solo%result%e_t, real(solo%result%iterations, dp)]
         one_t1 = solo%amps%t1
         one_t2 = solo%amps%t2
         opts%has_comm = .true.
         call solo%release()
      end if
      call bcast(comm, single, 3, 0)

      ! over the ranks, with the orbitals of the collective SCF
      topts%comm = comm
      topts%has_comm = .true.
      topts%verbose = .true.
      cmo = res%cmo(:, :, 1)
      call run_over_ranks(comm, label, bas, aux, pl, cmo, res%eps(:, 1), opts, topts, la, single, one_t1, one_t2, refs, good)
      if (rank == 0 .and. .not. good) print "(a)", "rccsd_mpi: "//label//"FAILED: converged/T/bitwise/pyscf/uniform check"
      ok = ok .and. good

      ! the guard: ranks after the first get orbitals with one occupied and one virtual column of the other sign
      if (im == 1) then
         topts%verbose = .false.
         if (rank > 0) then
            cmo(:, nocc) = -cmo(:, nocc)
            cmo(:, nocc + 1) = -cmo(:, nocc + 1)
         end if
         call run_over_ranks(comm, label//"sign-flipped orbitals on ranks > 0, ", bas, aux, pl, cmo, res%eps(:, 1), opts, &
                             topts, la, single, one_t1, one_t2, refs, good)
         if (rank == 0 .and. .not. good) print "(a)", "rccsd_mpi: "//label//"FAILED: the sign-flip guard case"
         ok = ok .and. good
         ! several passes inside every rank's share: the cut must not depend on d_batch beyond round-off
         cmo = res%cmo(:, :, 1)
         opts%ladder%d_batch = 2
         call run_over_ranks(comm, label//"d_batch 2, ", bas, aux, pl, cmo, res%eps(:, 1), opts, topts, la, single, &
                             one_t1, one_t2, refs, good)
         if (rank == 0 .and. .not. good) print "(a)", "rccsd_mpi: "//label//"FAILED: the d_batch 2 case"
         ok = ok .and. good
      end if
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine success_case

   subroutine prefix_case(comm, im, ok)
      !! The energy and the largest step after ONE and after TWO iterations (not converged), over the ranks against one
      !! rank, to 1e-14 (header). A follower that works from the previous iteration's amplitudes, or from none, is a
      !! different iteration 1 or 2.
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: im
      logical, intent(inout) :: ok

      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: la
      type(trc_rccsd_state_t) :: state
      type(trc_rccsd_options_t) :: opts
      integer, allocatable :: zint(:)
      integer :: nocc, rank, k
      logical :: supported, good
      real(dp) :: one(2), gap

      rank = comm%rank()
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) return
      call scf_over_comm(comm, bas, nocc, res)
      call la%init(2)
      do k = 1, 2
         call converged_options(comm, nocc, MOL_NFROZEN(im), opts)
         opts%max_iterations = k
         one = 0.0_dp
         if (rank == 0) then
            opts%has_comm = .false.
            call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
            one = [state%result%e_ccsd, state%result%step_max]
            call state%release()
            opts%has_comm = .true.
         end if
         call bcast(comm, one, 2, 0)
         call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
         gap = max(abs(state%result%e_ccsd - one(1)), abs(state%result%step_max - one(2)))
         good = all_ranks(comm, gap <= TOL_ONE_RANK_E .and. state%result%iterations == k)
         if (rank == 0) print "(a,i0,a,i0,a,es9.2,a,l1)", "rccsd_mpi: after ", k, " iteration(s) over ", comm%size(), &
            " ranks, max |E_CCSD, step - one rank| ", gap, ": ", good
         ok = ok .and. good
         call state%release()
      end do
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine prefix_case

   subroutine iteration_failure_case(comm, state, la, ok)
      !! Case 2d of the header, on the converged `state` of every rank (its integrals and amplitudes)
      type(comm_t), intent(in) :: comm
      type(trc_rccsd_state_t), intent(inout) :: state
      type(trc_linalg_t), intent(inout) :: la
      logical, intent(inout) :: ok

      type(trc_rccsd_update_t) :: upd
      type(trc_rccsd_amps_t) :: newamps
      type(trc_rccsd_ladder_options_t) :: lopts
      type(step_table_t) :: short_table
      integer, allocatable :: lo(:), hi(:)
      integer :: rank, n, no, nv, k, bad_rank
      logical :: good
      real(dp) :: gap

      rank = comm%rank()
      n = comm%size()
      if (n < 2) return
      no = state%ints%no
      nv = state%ints%nv
      lopts%form = RCCSD_LADDER_HALF
      call newamps%init(no, nv)
      allocate (lo(n), hi(n))
      call short_table%init("rccsd", RCCSD_STEP_NAME)
      ! misuse guard: set_share(lead) on an update initialised as a follower is refused through its status
      call upd%init(state%ints, lopts, follower=.true.)
      call rccsd_share_ranges(no, nv, state%ints%naux, n, lo, hi)
      call upd%set_share(comm, .true., lo, hi)
      good = upd%status%has_error() .and. .not. upd%distributed .and. &
             index(upd%status%get_message(), "initialised as a follower") > 0
      if (rank == 0) print "(a,l1)", "rccsd_mpi: fail 2d, set_share(lead) on a follower-initialised update refused: ", good
      ok = ok .and. good
      call upd%release()
      do k = 1, 5
         ! rounds 1 and 2: the last rank, then rank 0, is given a range beyond nv; round 3: every range correct
         bad_rank = merge(n - 1, 0, k == 1)
         call upd%init(state%ints, lopts, follower=rank /= 0)
         call rccsd_share_ranges(no, nv, state%ints%naux, n, lo, hi)
         if (k < 3) hi(bad_rank + 1) = nv + 7
         if (k == 4) then
            ! the last rank gets nothing: the others share 1..nv by the rule for n - 1 ranks
            call rccsd_share_ranges(no, nv, state%ints%naux, n - 1, lo(1:n - 1), hi(1:n - 1))
            lo(n) = nv + 1
            hi(n) = nv
         end if
         call upd%set_share(comm, rank == 0, lo, hi)
         if (k == 5) then
            if (rank == 0) then
               call upd%step(state%ints, state%amps, newamps, la, short_table)
            else
               call upd%follow(state%ints, state%amps, la, short_table)
            end if
         else if (rank == 0) then
            call upd%step(state%ints, state%amps, newamps, la)
         else
            call upd%follow(state%ints, state%amps, la)
         end if
         if (k == 5) then
            good = upd%agreed_stop .and. upd%status%has_error()
            good = all_ranks(comm, good) .and. has_text(comm, upd%status%get_message(), "no rows for the transfers")
            if (rank == 0) print "(a,l1)", "rccsd_mpi: fail 2d, a step table without the transfer rows refused on all ranks: ", good
         else if (k < 3) then
            good = upd%agreed_stop .and. upd%status%has_error()
            good = all_ranks(comm, good) .and. has_text(comm, upd%status%get_message(), "rank "//char(ichar("0") + bad_rank)// &
                                                        ": trc_rccsd_ladder: the share d")
            if (rank == 0) print "(a,i0,a,l1)", "rccsd_mpi: fail 2d, the share of rank ", bad_rank, &
               " refused inside an iteration, all ranks stopped with its reason: ", good
         else
            good = .not. upd%agreed_stop .and. .not. upd%status%has_error()
            gap = 0.0_dp
            if (rank == 0) then
               call newamps%pull()
               call state%amps%pull()
               gap = max(maxval(abs(newamps%t2 - state%amps%t2)), maxval(abs(newamps%t1 - state%amps%t1)))
               good = good .and. gap < 1.0e-9_dp
            end if
            good = all_ranks(comm, good)
            if (rank == 0) print "(a,i0,a,es9.2,a,l1)", "rccsd_mpi: fail 2d, a healthy iteration over the ranks (case ", k, &
               ": 3 the rule's cut, 4 the last rank's range empty) reproduces the converged amplitudes, max |change| ", gap, ": ", good
         end if
         ok = ok .and. good
         call upd%release()
      end do
      call newamps%release()
   end subroutine iteration_failure_case

   subroutine skew_start_case(comm, im, ok)
      !! A start whose t2 is NOT pair-symmetric: the whole run takes the general path (full arrays on the wire, two
      !! contractions per pass). Over the ranks it must reach the one-rank energy and iteration count of the same start.
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: im
      logical, intent(inout) :: ok

      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: la
      type(trc_rccsd_state_t) :: solo, multi
      type(trc_rccsd_options_t) :: opts
      integer, allocatable :: zint(:)
      integer :: nocc, rank, no, nv
      logical :: supported, good
      real(dp) :: e1(2), gap

      rank = comm%rank()
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) return
      call scf_over_comm(comm, bas, nocc, res)
      call la%init(2)
      call converged_options(comm, nocc, MOL_NFROZEN(im), opts)
      no = nocc - MOL_NFROZEN(im)
      nv = bas%nao_sph - nocc
      e1 = 0.0_dp
      if (rank == 0) then
         opts%has_comm = .false.
         call skew_amps(solo, no, nv)
         call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, solo)
         e1 = [solo%result%e_ccsd, real(solo%result%iterations, dp)]
         call solo%release()
         opts%has_comm = .true.
      end if
      call bcast(comm, e1, 2, 0)
      call skew_amps(multi, no, nv)
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, multi)
      gap = abs(multi%result%e_ccsd - e1(1))
      good = multi%result%converged .and. gap <= TOL_ONE_RANK_E .and. nint(e1(2)) == multi%result%iterations
      good = all_ranks(comm, good) .and. same_on_ranks(comm, multi%result%e_ccsd)
      if (rank == 0) print "(a,i0,a,es9.2,a,l1)", "rccsd_mpi: non-symmetric start, "//trim(MOL_LABEL(im))//", ranks ", comm%size(), &
         "  |E_CCSD - one rank| ", gap, ": ", good
      ok = ok .and. good
      call multi%release()
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine skew_start_case

   subroutine skew_amps(state, no, nv)
      !! A start with t1 = 0 and a t2 that has no pair symmetry (check_rccsd_run's)
      type(trc_rccsd_state_t), intent(inout) :: state
      integer, intent(in) :: no, nv

      integer :: i, j, a, b

      call state%amps%init(no, nv)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  state%amps%t2(i, j, a, b) = 1.0e-2_dp*sin(real(i + 2*j + 3*a + 5*b, dp))
               end do
            end do
         end do
      end do
      state%use_start = .true.
   end subroutine skew_amps

   function has_text(comm, message, text) result(everywhere)
      !! True when `message` contains `text` on every rank; collective
      type(comm_t), intent(in) :: comm
      character(len=*), intent(in) :: message, text
      logical :: everywhere

      everywhere = all_ranks(comm, index(message, text) > 0)
   end function has_text

   subroutine failure_cases(comm, im, ok)
      !! Argument "fail" mode (program header)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: im
      logical, intent(inout) :: ok

      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: la
      type(trc_rccsd_state_t) :: state
      type(trc_rccsd_options_t) :: opts
      type(trc_rccsd_triples_options_t) :: topts
      type(error_t) :: status
      type(rccsd_refs_t) :: refs
      integer, allocatable :: zint(:)
      integer :: nocc, rank, nranks, failing, k
      logical :: supported, good

      rank = comm%rank()
      nranks = comm%size()
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) return
      call scf_over_comm(comm, bas, nocc, res)
      call la%init(2)
      call read_rccsd_refs(trim(MOL_REF(im)), refs)
      topts%use_device = .true.
      topts%device_routine = RCCSD_TRIPLES_DEVICE_IJK
      topts%comm = comm
      topts%has_comm = .true.

      ! 1. does not converge: every rank sees it, and (T) is refused on every rank
      call converged_options(comm, nocc, MOL_NFROZEN(im), opts)
      opts%max_iterations = 2
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
      call rccsd_triples_finish(state, status, topts)
      good = .not. state%result%converged .and. .not. state%result%refused .and. status%has_error() .and. &
             .not. state%result%has_triples .and. state%result%iterations == 2 .and. allocated(state%amps%t2)
      good = all_ranks(comm, good) .and. same_on_ranks(comm, state%result%e_ccsd)
      good = good .and. has_text(comm, status%get_message(), "did not converge") .and. &
             has_text(comm, state%result%message, "not converged in 2 iterations")
      if (rank == 0) print "(a,l1)", "rccsd_mpi: fail 1, not converged on all ranks, (T) refused on all, with the reason: ", good
      ok = ok .and. good
      call state%release()

      ! 2. one rank is refused (the last, then rank 0): every rank returns refused with ITS reason, nothing left allocated
      do k = 1, 2
         failing = merge(nranks - 1, 0, k == 1)
         call converged_options(comm, nocc, MOL_NFROZEN(im), opts)
         if (rank == failing) opts%nocc = 0
         call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
         good = state%result%refused .and. .not. state%result%converged .and. state%ints%stamp == 0 .and. &
                .not. allocated(state%amps%t1)
         good = all_ranks(comm, good) .and. has_text(comm, state%result%message, "opts%nocc must be")
         if (rank == 0) print "(a,i0,a,l1)", "rccsd_mpi: fail 2, rank ", failing, &
            " refused, all ranks refused with its message: ", good
         ok = ok .and. good
         call state%release()
      end do

      ! 2b. the ladder share of one rank refuses (a negative d_batch on that rank only): every rank refuses with its reason
      do k = 1, 2
         failing = merge(nranks - 1, 0, k == 1)
         call converged_options(comm, nocc, MOL_NFROZEN(im), opts)
         if (rank == failing) opts%ladder%d_batch = -1
         call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
         good = state%result%refused .and. .not. state%result%converged .and. state%ints%stamp == 0 .and. &
                .not. allocated(state%amps%t1)
         good = all_ranks(comm, good) .and. has_text(comm, state%result%message, "negative batch size")
         if (rank == 0) print "(a,i0,a,l1)", "rccsd_mpi: fail 2b, the ladder share of rank ", failing, &
            " refused, all ranks refused with its message: ", good
         ok = ok .and. good
         call state%release()
      end do

      ! 2c. an explicit assembled ladder cannot be cut over ranks: refused on every rank, or (on one rank only) on all
      if (nranks > 1) then
         do k = 1, 2
            call converged_options(comm, nocc, MOL_NFROZEN(im), opts)
            if (k == 1 .or. rank == nranks - 1) opts%ladder%form = RCCSD_LADDER_ASSEMBLED
            call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
            good = state%result%refused .and. .not. state%result%converged .and. state%ints%stamp == 0 .and. &
                   .not. allocated(state%amps%t1)
            good = all_ranks(comm, good) .and. has_text(comm, state%result%message, "needs the half form")
            if (rank == 0) print "(a,i0,a,l1)", "rccsd_mpi: fail 2c, the assembled ladder refused on all ranks, case ", k, &
               " (1: every rank asks, 2: only the last), with the reason: ", good
            ok = ok .and. good
            call state%release()
         end do
      end if

      ! 3. converged state, then the (T) refusals
      call converged_options(comm, nocc, MOL_NFROZEN(im), opts)
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
      if (.not. state%result%converged) error stop "rccsd_mpi: the converged run did not converge"
      call iteration_failure_case(comm, state, la, ok)
      if (nranks > 1) then
         ! a (T) routine one rank asks for differently: all refuse
         topts%device_routine = merge(RCCSD_TRIPLES_DEVICE_ABC, RCCSD_TRIPLES_DEVICE_IJK, rank == nranks - 1)
         call rccsd_triples_finish(state, status, topts)
         good = all_ranks(comm, status%has_error() .and. .not. state%result%has_triples) .and. &
                has_text(comm, status%get_message(), "occupied-driven")
         if (rank == 0) print "(a,l1)", "rccsd_mpi: fail 3, differing (T) routine refused on all ranks, with the reason: ", good
         ok = ok .and. good
         ! the host reference has no work share: refused on every rank
         topts%device_routine = RCCSD_TRIPLES_DEVICE_IJK
         topts%use_device = .false.
         call rccsd_triples_finish(state, status, topts)
         good = all_ranks(comm, status%has_error() .and. .not. state%result%has_triples) .and. &
                has_text(comm, status%get_message(), "occupied-driven")
         if (rank == 0) print "(a,l1)", "rccsd_mpi: fail 4, host reference refused on all ranks, with the reason: ", good
         ok = ok .and. good
      end if
      topts%use_device = .true.
      topts%device_routine = RCCSD_TRIPLES_DEVICE_IJK
      ! one rank's state is marked not converged: all refuse with that reason
      if (rank == nranks - 1) state%result%converged = .false.
      call rccsd_triples_finish(state, status, topts)
      good = all_ranks(comm, status%has_error() .and. .not. state%result%has_triples) .and. &
             has_text(comm, status%get_message(), "did not converge")
      if (rank == 0) print "(a,l1)", "rccsd_mpi: fail 5, one rank's state not converged, all ranks refused, with the reason: ", good
      ok = ok .and. good
      state%result%converged = .true.
      ! one rank sets the work share itself
      if (rank == nranks - 1) topts%device%parts = 2
      call rccsd_triples_finish(state, status, topts)
      good = all_ranks(comm, status%has_error() .and. .not. state%result%has_triples) .and. &
             has_text(comm, status%get_message(), "device%parts")
      if (rank == 0) print "(a,l1)", "rccsd_mpi: fail 6, one rank sets parts, all ranks refused, with the reason: ", good
      ok = ok .and. good
      topts%device%parts = 1
      ! and the same state runs a correct (T) on every rank
      call rccsd_triples_finish(state, status, topts)
      good = all_ranks(comm, .not. status%has_error() .and. state%result%has_triples .and. &
                       abs(state%result%e_t - refs%e_t) < TOL_T)
      if (rank == 0) print "(a,l1)", "rccsd_mpi: after the refusals the same state runs a correct (T) on all ranks: ", good
      ok = ok .and. good
      call state%release()
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine failure_cases

end module rccsd_mpi_cases

program rccsd_mpi
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize
   use trc_api, only: trc_bind_device
   use rccsd_test_support, only: NMOL
   use rccsd_mpi_cases, only: success_case, failure_cases, skew_start_case, prefix_case
   implicit none

   type(comm_t) :: comm
   integer :: dev, im
   logical :: ok, fail_mode
   character(len=16) :: arg

   call pic_mpi_init()
   comm = comm_world()
   dev = trc_bind_device(comm%rank())
   fail_mode = .false.
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      fail_mode = trim(arg) == "fail"
   end if
   ok = .true.
   if (fail_mode) then
      call failure_cases(comm, 1, ok)
   else
      do im = 1, NMOL
         call success_case(comm, im, ok)
      end do
      call skew_start_case(comm, 1, ok)
      call prefix_case(comm, 1, ok)
   end if
   if (comm%rank() == 0) print "(a)", merge("rccsd_mpi: PASS", "rccsd_mpi: FAIL", ok)
   call pic_mpi_finalize()
   if (.not. ok) stop 1
end program rccsd_mpi
