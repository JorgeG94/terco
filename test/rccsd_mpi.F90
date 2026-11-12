!
! RI-CCSD(T) over MPI ranks, one rank per device (src/post_hf/ccsd/trc_rccsd_driver.F90 and trc_rccsd_triples.F90,
! RANKS). Water and a C1 hydrogen peroxide, cc-pVDZ with cc-pVDZ-RIFIT (check_rccsd_triples' molecules and pyscf
! references, read through rccsd_ref_energies). The SCF runs over the communicator, as the CLI and the C API run it
! (test/scf_mpi.F90). The program is generic in the rank count; ctest runs it on 1 and 2, on a host MPI build too (the
! (T) routine is the occupied-driven one, which runs on the host).
!
! No argument: rank 0 runs the CCSD(T) alone first (no communicator, the single-rank code path), then every rank runs
! it over the communicator. Required, on every rank:
!   - converged, (T) done, no error;
!   - E_CCSD equal to the single-rank value in every bit (rank 0 computes it alone, the others receive it);
!   - EVERY element of t1 and t2, read back from the rank's device, equal to rank 0's single-rank amplitudes in every
!     bit (the broadcast and the push to the device happened);
!   - E(T) within TOL_TRIPLES of the single-rank value: the per-triple energies of the ranks are those of one rank to
!     rounding, and the sum of the reduced vector differs from the one-rank sum only in rounding;
!   - the energies against pyscf;
!   - E_CCSD, E(T), their sum and the iteration count equal over the ranks (MIN equals MAX).
! Water is run once more with the sign of one occupied and one virtual orbital column flipped on every rank but the
! first (the GUARD case: a caller that hands ranks different orbitals). The same energies and amplitudes are required;
! without the orbital broadcast of trc_rccsd_solve the (T) of the other ranks is wrong and this case fails.
!
! Argument "fail": the failure paths, which must end on every rank without a rank left waiting (ctest puts a timeout
! on it), and with the TEXT of the refusing rank's reason on every rank: a CCSD that does not converge (every rank sees
! it; (T) is refused on all of them); a refusal on ONE rank (the last, and rank 0) that makes every rank refuse; a (T)
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
   use trc_rccsd_triples, only: trc_rccsd_triples_options_t, rccsd_triples_finish, RCCSD_TRIPLES_DEVICE_IJK, &
                                RCCSD_TRIPLES_DEVICE_ABC
   use rccsd_test_support, only: NMOL, MOL_LABEL, MOL_REF, MOL_NFROZEN, load_molecule
   use rccsd_ref_energies, only: rccsd_refs_t, read_rccsd_refs
   implicit none
   private
   public :: success_case, failure_cases

   ! Mirrors of constants defined in the single-rank checks (test programs cannot share them without a module of their
   ! own): TOL_CORR, TOL_T, TOL_TOTAL are check_capi_rccsd's; TOL_TRIPLES is check_rccsd_triples' TOL_IJK_TILE.
   real(dp), parameter :: TOL_CORR = 1.0e-10_dp   !! E_CCSD against pyscf
   real(dp), parameter :: TOL_T = 3.0e-12_dp      !! E(T) against pyscf
   real(dp), parameter :: TOL_TOTAL = 3.0e-11_dp  !! the CCSD(T) total against pyscf
   real(dp), parameter :: TOL_TRIPLES = 1.0e-16_dp  !! E(T) over the ranks against one rank, absolute

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
      !! numbers `single` = (E_CCSD, E(T)) and amplitudes of rank 0 and against pyscf; `good` is the same on every rank
      type(comm_t), intent(in) :: comm
      character(len=*), intent(in) :: label
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: cmo(:, :), eps(:)
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_rccsd_triples_options_t), intent(in) :: topts
      type(trc_linalg_t), intent(inout) :: la
      real(dp), intent(in) :: single(2)
      real(dp), intent(inout) :: one_t1(:, :), one_t2(:, :, :, :)
      type(rccsd_refs_t), intent(in) :: refs
      logical, intent(out) :: good

      type(trc_rccsd_state_t) :: state
      type(error_t) :: status
      logical :: uniform
      real(dp) :: gap, e_total
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
      if (rank == 0) then
         print "(a,i0,a,f20.12,a,f20.14,a,es9.2,a,es9.2,a,es9.2)", "rccsd_mpi: "//label//"ranks ", comm%size(), &
            "  E_CCSD ", state%result%e_ccsd, "  E(T) ", state%result%e_t, "  |E(T) - one rank| ", &
            abs(state%result%e_t - single(2)), "  |E(T) - pyscf| ", abs(state%result%e_t - refs%e_t), &
            "  max |amplitude - one rank| ", gap
      end if
      good = good .and. state%result%e_ccsd == single(1) .and. gap == 0.0_dp .and. &
             abs(state%result%e_t - single(2)) <= TOL_TRIPLES .and. &
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
      real(dp) :: single(2)
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
         single = [solo%result%e_ccsd, solo%result%e_t]
         one_t1 = solo%amps%t1
         one_t2 = solo%amps%t2
         opts%has_comm = .true.
         call solo%release()
      end if
      call bcast(comm, single, 2, 0)

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
      end if
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine success_case

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
      good = good .and. has_text(comm, status%get_message(), "did not converge")
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

      ! 3. converged state, then the (T) refusals
      call converged_options(comm, nocc, MOL_NFROZEN(im), opts)
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
      if (.not. state%result%converged) error stop "rccsd_mpi: the converged run did not converge"
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
   use rccsd_mpi_cases, only: success_case, failure_cases
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
   end if
   if (comm%rank() == 0) print "(a)", merge("rccsd_mpi: PASS", "rccsd_mpi: FAIL", ok)
   call pic_mpi_finalize()
   if (.not. ok) stop 1
end program rccsd_mpi
