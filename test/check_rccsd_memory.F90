module check_rccsd_memory_suite
   !!
   !! A RI-CCSD(T) run whose arrays do not fit on the device is REFUSED, with a message that names the allocation and the
   !! bytes it needs and the bytes free, instead of aborting inside the runtime (trc_device_memory). A host build has no
   !! device, so the free memory is capped by the test hook (`device_free_cap` of the CCSD and (T) options, or
   !! `device_memory_set_cap` for a type used on its own); every case below is a size check that the host build
   !! cannot otherwise reach.
   !!
   !! COVERAGE. The refusals of the integral build, of each type's `init` (called directly) and of the (T) are reached;
   !! the driver's forwarding of a LATER stage's refusal (`stages_setup`, `follower_amps`) is not, because the ints build is
   !! the largest stage and one static cap cannot pass it and refuse a later one.
   !!
   !! Water, cc-pVDZ with cc-pVDZ-RIFIT, frozen core 1 (o = 4, v = 19), the molecule and SCF of the other RI-CCSD checks.
   !!
   !! Cases:
   !!   1 `trc_rccsd_solve` with a cap of 1000 bytes: refused with the integral build's message, every type empty; the
   !!     same state then runs without the cap
   !!   2 the W scratch: refused one byte below the bytes it maps (5 o^2 v^2 + max(o,v) o naux doubles, formed here from
   !!     the sizes), accepted at that many, nothing allocated after the refusal. This is the allocation that aborted a
   !!     V100 at caffeine cc-pVTZ
   !!   3 the amplitudes, and the DIIS history: refused, empty; the history on the host fits where the one on the device
   !!     does not
   !!   4 the amplitude update (`init` of everything the iteration holds): refused at the ladder scratch, after the earlier parts are mapped, and every part released
   !!   5 (T): the converged CCSD, then `rccsd_triples_finish` with a cap below the scratch of the occupied-driven routine
   !!     and of the tiled one: the failure is the (T)'s, naming its scratch, the CCSD result is as it was and has_triples
   !!     is false; the tiled routine's dense (ia|bc) block is refused after a scratch that fit; with the cap lifted the
   !!     same state runs (T) and agrees with the host reference
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_error, only: error_t
   use trc_linalg, only: trc_linalg_t
   use trc_device_memory, only: device_memory_set_cap
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_scratch_t
   use trc_rccsd_w, only: trc_rccsd_wscratch_t
   use trc_rccsd_diis, only: trc_rccsd_diis_t
   use trc_rccsd_update, only: trc_rccsd_update_t
   use trc_rccsd_driver, only: trc_rccsd_options_t, trc_rccsd_state_t, trc_rccsd_solve
   use trc_rccsd_triples, only: trc_rccsd_triples_options_t, rccsd_triples_finish, RCCSD_TRIPLES_DEVICE_ABC, &
                                RCCSD_TRIPLES_DEVICE_IJK
   use rccsd_test_support, only: MOL_NFROZEN, load_molecule, run_scf
   implicit none
   private
   public :: collect_check_rccsd_memory

   integer, parameter :: WATER = 1
   logical :: computed = .false.
   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl
   type(trc_scf_result_t) :: res
   type(trc_linalg_t) :: la
   integer :: nocc = 0

contains

   subroutine collect_check_rccsd_memory(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("solve refused: the integral blocks do not fit", test_solve_refused), &
                  new_unittest("W scratch: refused one byte short, accepted at its size", test_w_scratch), &
                  new_unittest("amplitudes and DIIS history: refused, the host history fits", test_amps_diis), &
                  new_unittest("update: refused, every part released", test_update), &
                  new_unittest("(T): refused, the CCSD result kept, runs once the cap is lifted", test_triples) &
                  ]
   end subroutine collect_check_rccsd_memory

   subroutine ensure_computed()
      integer, allocatable :: zint(:)
      integer :: dev
      logical :: supported

      if (computed) return
      dev = trc_bind_device(0)
      call load_molecule(WATER, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) then
         print "(a)", "SKIP check_rccsd_memory: this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
         stop 0
      end if
      call run_scf(bas, nocc, res)
      call la%init(2)
      computed = .true.
   end subroutine ensure_computed

   function options() result(opts)
      type(trc_rccsd_options_t) :: opts

      opts%nocc = nocc
      opts%nfrozen = MOL_NFROZEN(WATER)
   end function options

   subroutine build_ints(ints)
      type(trc_rccsd_ints_t), intent(out) :: ints

      type(trc_rccsd_ints_options_t) :: bopts

      bopts%nalpha = nocc
      bopts%nbeta = nocc
      bopts%nfrozen = MOL_NFROZEN(WATER)
      call device_memory_set_cap(0_int64)
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
   end subroutine build_ints

   subroutine names(error, message, parts, label)
      !! Every text of `parts` is in `message`
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: message, parts(:), label

      integer :: k

      do k = 1, size(parts)
         call check(error, index(message, trim(parts(k))) > 0, label//": '"//trim(parts(k))//"' is not in '"//message//"'")
         if (allocated(error)) return
      end do
   end subroutine names

   subroutine test_solve_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_state_t) :: state
      type(trc_rccsd_options_t) :: opts

      call ensure_computed()
      opts = options()
      opts%device_free_cap = 1000_int64
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
      call check(error, state%result%refused, "a run with 1000 free bytes must be refused")
      if (allocated(error)) return
      call names(error, state%result%message, [character(len=60) :: "fitted integral blocks", "GB of device memory", &
                                                "GB are free", "smaller basis"], "solve")
      if (allocated(error)) return
      call check(error, .not. state%result%converged .and. state%ints%stamp == 0 .and. .not. allocated(state%ints%ovov) .and. &
                 .not. allocated(state%ints%b_vv) .and. .not. allocated(state%amps%t1), "a refusal leaves nothing allocated")
      if (allocated(error)) return
      call state%release()
      ! the hook is per call: the next run, without the cap, is not refused
      opts%device_free_cap = 0_int64
      opts%max_iterations = 1
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
      call check(error, .not. state%result%refused, "the cap of one call must not outlive it: "//trim(state%result%message))
      call state%release()
   end subroutine test_solve_refused

   subroutine test_w_scratch(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_wscratch_t) :: ws
      integer(int64) :: need

      call ensure_computed()
      call build_ints(ints)
      call check(error, .not. ints%status%has_error(), "the integrals build")
      if (allocated(error)) return
      need = 8_int64*(5_int64*ints%no*ints%no*ints%nv*ints%nv + int(max(ints%no, ints%nv), int64)*ints%no*ints%naux)
      call device_memory_set_cap(need - 1)
      call ws%init(ints)
      call device_memory_set_cap(0_int64)
      call check(error, ws%status%has_error(), "one byte short of the W scratch must be refused")
      if (allocated(error)) return
      call names(error, ws%status%get_message(), [character(len=60) :: "trc_rccsd_w: the W scratch", "GB of device memory"], "W")
      if (allocated(error)) return
      call check(error, .not. allocated(ws%akl) .and. .not. allocated(ws%f1) .and. ws%stamp == 0, "nothing was allocated")
      if (allocated(error)) return
      call device_memory_set_cap(need)
      call ws%init(ints)
      call device_memory_set_cap(0_int64)
      call check(error, .not. ws%status%has_error() .and. allocated(ws%akl), "the exact size is accepted")
      if (allocated(error)) return
      ! the formula is the bytes of the arrays that were mapped, counted from the arrays
      call check(error, need == ws_bytes(ws), "the bytes asked for are those of the arrays allocated: "//to_string(need)// &
                 " against "//to_string(ws_bytes(ws)))
      call ws%release()
      call ints%release()
   end subroutine test_w_scratch

   function ws_bytes(ws) result(bytes)
      type(trc_rccsd_wscratch_t), intent(in) :: ws
      integer(int64) :: bytes

      bytes = 8_int64*(size(ws%akl, kind=int64) + size(ws%a2, kind=int64) + size(ws%b1, kind=int64) + &
                       size(ws%b2, kind=int64) + size(ws%r, kind=int64) + size(ws%f1, kind=int64))
   end function ws_bytes

   function scr_bytes(scr) result(bytes)
      type(trc_rccsd_scratch_t), intent(in) :: scr
      integer(int64) :: bytes

      bytes = 8_int64*(size(scr%lov, kind=int64) + size(scr%pack1, kind=int64) + size(scr%lq, kind=int64) + &
                       size(scr%lam, kind=int64) + size(scr%part_s, kind=int64) + size(scr%part_d, kind=int64) + &
                       size(scr%n1, kind=int64) + size(scr%n1p, kind=int64) + size(scr%bovp, kind=int64) + &
                       size(scr%u, kind=int64))
   end function scr_bytes

   subroutine test_amps_diis(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_amps_t) :: amps
      type(trc_rccsd_diis_t) :: diis
      type(trc_rccsd_ints_t) :: ints
      integer(int64) :: need
      integer :: no, nv

      call ensure_computed()
      call build_ints(ints)
      no = ints%no
      nv = ints%nv
      call ints%release()
      ! the amplitudes: the bytes asked for are those of t1 and t2
      call amps%init(no, nv)
      need = 8_int64*(size(amps%t1, kind=int64) + size(amps%t2, kind=int64))
      call amps%release()
      call device_memory_set_cap(need - 1)
      call amps%init(no, nv)
      call device_memory_set_cap(0_int64)
      call check(error, amps%status%has_error() .and. .not. allocated(amps%t1) .and. .not. allocated(amps%t2), &
                 "amplitudes one byte short must be refused, nothing allocated")
      if (allocated(error)) return
      call names(error, amps%status%get_message(), &
                 [character(len=60) :: "trc_rccsd_inter: the amplitudes", "GB of device memory"], "amps")
      if (allocated(error)) return
      call device_memory_set_cap(need)
      call amps%init(no, nv)
      call device_memory_set_cap(0_int64)
      call check(error, .not. amps%status%has_error() .and. allocated(amps%t2), "the exact size of the amplitudes is accepted")
      call amps%release()
      if (allocated(error)) return
      ! the device history: the bytes asked for are those of err, vec, part and coef
      call diis%init(no, nv, 8, packed=.false., host=.false.)
      need = 8_int64*(size(diis%err, kind=int64) + size(diis%vec, kind=int64) + size(diis%part, kind=int64) + &
                      size(diis%coef, kind=int64))
      call diis%release()
      call device_memory_set_cap(need - 1)
      call diis%init(no, nv, 8, packed=.false., host=.false.)
      call device_memory_set_cap(0_int64)
      call check(error, diis%status%has_error() .and. .not. allocated(diis%err), "the device history one byte short is refused")
      if (allocated(error)) return
      call names(error, diis%status%get_message(), [character(len=60) :: "trc_rccsd_diis: the DIIS history", "--cc-diis"], &
                 "diis")
      if (allocated(error)) return
      call device_memory_set_cap(need)
      call diis%init(no, nv, 8, packed=.false., host=.false.)
      call device_memory_set_cap(0_int64)
      call check(error, .not. diis%status%has_error() .and. allocated(diis%err), "the exact size of the history is accepted")
      call diis%release()
      if (allocated(error)) return
      ! the history on the host asks for far less (two columns and the partials) and fits where the device one does not
      call device_memory_set_cap(need - 1)
      call diis%init(no, nv, 8, packed=.false., host=.true.)
      call device_memory_set_cap(0_int64)
      call check(error, .not. diis%status%has_error() .and. allocated(diis%err), "the history on the host fits where the "// &
                 "device history does not")
      call diis%release()
   end subroutine test_amps_diis

   subroutine test_update(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_update_t) :: upd
      integer(int64) :: need_scr, need_ws, need_w, need_lad, need_ls, cap

      call ensure_computed()
      call build_ints(ints)
      ! what each part asks for, counted from the arrays of an update built without a cap
      call upd%init(ints)
      call check(error, .not. upd%status%has_error(), "the update builds without a cap")
      if (allocated(error)) return
      need_scr = scr_bytes(upd%scr)
      need_ws = ws_bytes(upd%ws)
      need_w = 8_int64*(size(upd%w%woooo, kind=int64) + size(upd%w%wvoov, kind=int64) + size(upd%w%wvovo, kind=int64))
      need_lad = 8_int64*size(upd%lad%h, kind=int64)
      need_ls = upd%ls%bytes
      call upd%release()
      ! a cap that passes the scratch of the intermediates, W, the W scratch and the ladder result (each is mapped before
      ! the next is asked) and stops at the ladder's work arrays, the largest part
      cap = max(need_scr, need_ws, need_w, need_lad)
      call check(error, need_ls > cap, "the ladder's work arrays ("//to_string(need_ls)//" bytes) are the part that "// &
                 "does not fit under "//to_string(cap))
      if (allocated(error)) return
      call device_memory_set_cap(cap)
      call upd%init(ints)
      call device_memory_set_cap(0_int64)
      call check(error, upd%status%has_error(), "an update whose ladder scratch does not fit must be refused")
      if (allocated(error)) return
      call names(error, upd%status%get_message(), [character(len=60) :: "trc_rccsd_ladder: the ladder's work arrays", &
                                                    "GB of device memory", "GB are free"], "update")
      if (allocated(error)) return
      call check(error, .not. allocated(upd%scr%lov) .and. .not. allocated(upd%ws%akl) .and. .not. allocated(upd%w%woooo) .and. &
                 .not. allocated(upd%lad%h) .and. .not. allocated(upd%ls%taut) .and. .not. allocated(upd%inter%tau) .and. &
                 .not. allocated(upd%lamp) .and. upd%stamp == 0, "the parts mapped before the refusal are released")
      if (allocated(error)) return
      ! and a cap below the first part refuses before anything is mapped
      call device_memory_set_cap(1000_int64)
      call upd%init(ints)
      call device_memory_set_cap(0_int64)
      call check(error, upd%status%has_error() .and. .not. allocated(upd%scr%lov), "an update of 1000 free bytes is refused")
      if (allocated(error)) return
      call upd%init(ints)
      call check(error, .not. upd%status%has_error(), "without the cap the same update is built: "// &
                 merge(upd%status%get_message(), "", upd%status%has_error()))
      call upd%release()
      call ints%release()
   end subroutine test_update

   subroutine test_triples(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_state_t) :: state
      type(trc_rccsd_options_t) :: opts
      type(trc_rccsd_triples_options_t) :: topts
      type(error_t) :: status
      real(dp) :: e_ccsd, e_t_device, e_t_host

      call ensure_computed()
      opts = options()
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
      call check(error, state%result%converged, "the CCSD converges")
      if (allocated(error)) return
      e_ccsd = state%result%e_ccsd

      ! the occupied-driven routine: its scratch is a few MB here
      topts%use_device = .true.
      topts%device_routine = RCCSD_TRIPLES_DEVICE_IJK
      topts%device_free_cap = 1000_int64
      call rccsd_triples_finish(state, status, topts)
      call check(error, status%has_error() .and. .not. state%result%has_triples, "(T) must fail, not set an energy")
      if (allocated(error)) return
      call names(error, status%get_message(), &
                 [character(len=60) :: "trc_rccsd_triples_ijk: the (T) scratch", "GB of device memory"], "ijk")
      if (allocated(error)) return
      call check(error, state%result%converged .and. .not. state%result%refused .and. state%result%e_ccsd == e_ccsd .and. &
                 state%result%e_t == 0.0_dp .and. allocated(state%amps%t2) .and. state%ints%stamp /= 0, &
                 "the converged CCSD is untouched by the (T) failure")
      if (allocated(error)) return

      ! the tiled routine
      topts%device_routine = RCCSD_TRIPLES_DEVICE_ABC
      call rccsd_triples_finish(state, status, topts)
      call check(error, status%has_error() .and. .not. state%result%has_triples, "the tiled (T) must fail too")
      if (allocated(error)) return
      call names(error, status%get_message(), [character(len=60) :: "trc_rccsd_triples_dev: the (T) scratch"], "tiled")
      if (allocated(error)) return

      ! one triple a tile: the scratch (105 kB) fits 150 kB, the dense (ia|bc) block (219 kB) does not
      topts%device%tile_triples = 1
      topts%device_free_cap = 150000_int64
      call rccsd_triples_finish(state, status, topts)
      call check(error, status%has_error() .and. .not. state%result%has_triples, "the dense block must be refused")
      if (allocated(error)) return
      call names(error, status%get_message(), [character(len=60) :: "trc_rccsd_triples_dev: the dense (ia|bc) block"], "dense")
      if (allocated(error)) return

      ! lifted, the same state runs the (T) and agrees with the host reference
      topts%device%tile_triples = 0
      topts%device_free_cap = 0_int64
      topts%device_routine = RCCSD_TRIPLES_DEVICE_IJK
      call rccsd_triples_finish(state, status, topts)
      call check(error, .not. status%has_error() .and. state%result%has_triples, "(T) runs once the cap is lifted: "// &
                 merge(status%get_message(), "", status%has_error()))
      if (allocated(error)) return
      e_t_device = state%result%e_t
      topts%use_device = .false.
      call rccsd_triples_finish(state, status, topts)
      e_t_host = state%result%e_t
      call check(error, .not. status%has_error() .and. abs(e_t_device - e_t_host) < 1.0e-10_dp .and. e_t_host /= 0.0_dp, &
                 "E(T) of the device routine after a refusal equals the host reference")
      call state%release()
   end subroutine test_triples

end module check_rccsd_memory_suite

program check_rccsd_memory
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_memory_suite, only: collect_check_rccsd_memory
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_memory, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_memory
