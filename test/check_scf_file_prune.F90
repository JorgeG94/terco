!
! trc_scf_file (save/load restart) against a basis whose primitive pruning
! (trc_basis_json's prune_block, src/basis/trc_basis_json.F90) changed
! between the save and the load -- the scenario of an SCF saved by a build
! or a run before TRC_PRUNE_PRIMITIVES existed (or with it set to "0"),
! read back by a build or a run that prunes a general contraction (the
! default). The saved basis then has a different per-shell primitive count
! (and usually a different maxnp) from the one the loader builds, so it
! must refuse rather than hand back orbitals that do not belong to this AO
! space -- exactly what trc_scf_file's own header promises for "a different
! basis file of the same size".
!
! A test-drive suite.
!
module check_scf_file_prune_suite
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_file, only: trc_scf_save, trc_scf_load
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_scf_file_prune

   character(len=*), parameter :: SCRATCH = "check_scf_file_prune_tmp.bin"

   interface
      function c_setenv(name, value, overwrite) bind(c, name="setenv") result(r)
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: r
      end function c_setenv
   end interface

contains

   subroutine collect_scf_file_prune(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("loading an unpruned-basis SCF into a pruned basis is refused", &
                      test_mismatch_refused), &
         new_unittest("loading matches when both sides prune the same way", &
                      test_match_accepted) &
         ]
   end subroutine collect_scf_file_prune

   subroutine set_prune(v)
      character(*), intent(in) :: v
      integer(c_int) :: r
      r = c_setenv("TRC_PRUNE_PRIMITIVES"//c_null_char, trim(v)//c_null_char, 1_c_int)
      if (r /= 0) error stop 'check_scf_file_prune: setenv(TRC_PRUNE_PRIMITIVES) failed'
   end subroutine set_prune

   ! water/cc-pVDZ: oxygen's s and p blocks are general contractions
   ! (check_prune_primitives.F90's own header has the exponent counts), so
   ! TRC_PRUNE_PRIMITIVES actually changes this basis's shape, not just its
   ! numbers.
   subroutine build_and_run(prune_flag, bas, at_z, at_r, res)
      character(*), intent(in) :: prune_flag
      type(trc_basis_t), intent(out) :: bas
      real(dp), allocatable, intent(out) :: at_z(:), at_r(:, :)
      type(trc_scf_result_t), intent(out) :: res
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      integer :: natm
      integer, allocatable :: zint(:)

      call set_prune(prune_flag)
      call read_xyz('water.xyz', natm, zint, at_r)
      call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, bas, err, force_cart=.true.)
      if (err%has_error()) error stop 'check_scf_file_prune: basis build failed'
      at_z = real(zint, dp)
      call bas%to_device()
      opts%conv_energy = 1.0e-12_dp
      opts%conv_density = 1.0e-10_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      call trc_scf_run(bas, 5, 5, opts, res)   ! water, RHF, 10 electrons
      if (.not. res%converged) error stop 'check_scf_file_prune: RHF did not converge'
   end subroutine build_and_run

   !
   ! Save a run made with pruning off (the shape a pre-pruning build, or
   ! TRC_PRUNE_PRIMITIVES=0, would have written), then try to load it
   ! against a freshly built basis with pruning at its default (on). The
   ! two bases differ in oxygen's s and p shell primitive counts, so the
   ! load must come back with a non-empty refusal message and must not
   ! touch `res_loaded`.
   !
   subroutine test_mismatch_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas_saved, bas_now
      real(dp), allocatable :: z_saved(:), r_saved(:, :), z_now(:), r_now(:, :)
      type(trc_scf_result_t) :: res_saved, res_loaded
      character(len=256) :: msg

      call build_and_run("0", bas_saved, z_saved, r_saved, res_saved)
      call trc_scf_save(SCRATCH, bas_saved, z_saved, r_saved, 0, 1, 5, 5, "", "", 0.0_dp, res_saved, msg)
      call check(error, len_trim(msg) == 0, "trc_scf_save itself must not fail: "//trim(msg))
      if (allocated(error)) then
         call cleanup()
         return
      end if

      call build_and_run("1", bas_now, z_now, r_now, res_loaded)
      ! build_and_run already ran an SCF into res_loaded; trc_scf_load must
      ! refuse before touching it, so overwrite its energy with a sentinel
      ! and confirm the refusal leaves that sentinel alone.
      res_loaded%energy = -123.0_dp
      call trc_scf_load(SCRATCH, bas_now, z_now, r_now, 0, 1, res_saved%nspin, "", "", 0.0_dp, res_loaded, msg)
      call check(error, len_trim(msg) > 0, "a pruned-vs-unpruned basis mismatch must be refused, not accepted")
      if (.not. allocated(error)) then
         call check(error, res_loaded%energy, -123.0_dp, "a refused load must not touch the result")
      end if
      call cleanup()
   end subroutine test_mismatch_refused

   !
   ! The same two builds, both pruned (today's default both times) --
   ! confirms test_mismatch_refused's refusal is about the pruning
   ! difference and not some unrelated save/load bug.
   !
   subroutine test_match_accepted(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas_saved, bas_now
      real(dp), allocatable :: z_saved(:), r_saved(:, :), z_now(:), r_now(:, :)
      type(trc_scf_result_t) :: res_saved, res_loaded
      character(len=256) :: msg

      call build_and_run("1", bas_saved, z_saved, r_saved, res_saved)
      call trc_scf_save(SCRATCH, bas_saved, z_saved, r_saved, 0, 1, 5, 5, "", "", 0.0_dp, res_saved, msg)
      call check(error, len_trim(msg) == 0, "trc_scf_save itself must not fail: "//trim(msg))
      if (allocated(error)) then
         call cleanup()
         return
      end if

      call build_and_run("1", bas_now, z_now, r_now, res_loaded)
      call trc_scf_load(SCRATCH, bas_now, z_now, r_now, 0, 1, res_saved%nspin, "", "", 0.0_dp, res_loaded, msg)
      call check(error, len_trim(msg) == 0, "same pruning both sides must be accepted: "//trim(msg))
      if (.not. allocated(error)) then
         call check(error, abs(res_loaded%energy - res_saved%energy) < 1.0e-10_dp, &
                    "a matching load must hand back the saved energy")
      end if
      call cleanup()
   end subroutine test_match_accepted

   subroutine cleanup()
      integer :: u, ios
      logical :: there
      inquire (file=SCRATCH, exist=there)
      if (.not. there) return
      open (newunit=u, file=SCRATCH, status="old", iostat=ios)
      if (ios == 0) close (u, status="delete")
   end subroutine cleanup

end module check_scf_file_prune_suite

program check_scf_file_prune
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_scf_file_prune_suite, only: collect_scf_file_prune
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_file_prune, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_scf_file_prune
