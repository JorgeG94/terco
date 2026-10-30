!
! Pruning a general contraction's redundant primitives (trc_basis_json's
! prune_block, src/basis/trc_basis_json.F90) against pyscf's own
! "optimized general contraction" cc-pVDZ counts, and against the energy
! the unpruned (BSE-literal) basis gives.
!
! water/cc-pVDZ: oxygen's s block is 3 columns sharing 9 exponents (two
! real contractions, one single-primitive diffuse column); its p block is
! 2 columns sharing 4 exponents (one contraction, one single primitive).
! Both diffuse columns' exponent also shows up, with a small nonzero
! coefficient, inside the block's contracted column(s) -- the redundancy
! this module removes. Hydrogen's s block is the same shape one l lower.
! Counts below are `gto.basis.load('cc-pvdz', <symbol>)` (pyscf 2.14.0)
! against the BSE file metalquicha ships (test/basis_sets/cc-pvdz.json,
! identical to /home/jorge/dev/metalquicha/basis_sets/cc-pvdz.json).
!
! A test-drive suite.
!
module check_prune_primitives_suite
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_read_real, ref_close
   implicit none
   private
   public :: collect_prune_primitives

   interface
      function c_setenv(name, value, overwrite) bind(c, name="setenv") result(r)
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: r
      end function c_setenv
   end interface

contains

   subroutine collect_prune_primitives(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("pruned per-shell primitive counts match pyscf's cc-pVDZ", test_counts_match_pyscf), &
         new_unittest("the off switch reproduces the BSE-literal (padded) counts", test_off_switch_is_literal), &
         new_unittest("pruning leaves the RHF energy unchanged", test_energy_invariant), &
         new_unittest("the pruned energy matches the pyscf anchor", test_matches_anchor), &
         new_unittest("a column pruned to all zero is restored, not left NaN", &
                      test_allzero_column_restored) &
         ]
   end subroutine collect_prune_primitives

   subroutine set_prune(v)
      character(*), intent(in) :: v
      integer(c_int) :: r
      r = c_setenv("TRC_PRUNE_PRIMITIVES"//c_null_char, trim(v)//c_null_char, 1_c_int)
      if (r /= 0) error stop 'check_prune_primitives: setenv(TRC_PRUNE_PRIMITIVES) failed'
   end subroutine set_prune

   subroutine build_water_dz(bas)
      type(trc_basis_t), intent(out) :: bas
      type(error_t) :: err
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      call read_xyz('water.xyz', natm, zint, at_r)
      call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, bas, err, force_cart=.true.)
      if (err%has_error()) error stop 'check_prune_primitives: basis build failed'
   end subroutine build_water_dz

   !
   ! TRC_PRUNE_PRIMITIVES is read once per `trc_basis_from_json` call
   ! (trc_basis_json's prune_enabled), so toggling it between two builds in
   ! the same process is exactly the intended use -- test/check_rys_vs_hgp.F90
   ! does the same with TRC_RYS_CLASSES.
   !
   subroutine test_counts_match_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      ! O: s block (9 exponents, 3 columns: 2 contracted + 1 diffuse single
      ! primitive) -> 8, 8, 1; p block (4 exponents, 2 columns) -> 3, 1; the
      ! d is already a single primitive. H: s block (4 exponents, 2
      ! columns) -> 3, 1; the p is already a single primitive.
      integer, parameter :: L_EXP(12) = [0, 0, 0, 1, 1, 2, 0, 0, 1, 0, 0, 1]
      integer, parameter :: NP_EXP(12) = [8, 8, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1]
      integer :: i

      call set_prune("1")
      call build_water_dz(bas)
      call check(error, bas%nshell, 12, "water/cc-pVDZ has 12 shells before or after pruning")
      if (allocated(error)) return
      do i = 1, 12
         call check(error, bas%sh_l(i), L_EXP(i), "shell l unchanged by pruning")
         if (allocated(error)) return
         call check(error, bas%sh_np(i), NP_EXP(i), "pruned primitive count matches pyscf's cc-pVDZ")
         if (allocated(error)) return
      end do
   end subroutine test_counts_match_pyscf

   subroutine test_off_switch_is_literal(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      ! TRC_PRUNE_PRIMITIVES=0 only turns prune_block's elimination off;
      ! trc_api's basis_build drops a LITERALLY zero coefficient
      ! unconditionally (it always has -- "dropped here, once", that
      ! subroutine's own comment), so a column that was already a single
      ! primitive padded with zeros (BSE's own diffuse columns) still
      ! comes back as np=1 with pruning off. Only the two dense columns'
      ! own primitive count is what the switch controls: 9 and 4 here,
      ! not 8 and 3. This is today's behaviour before this module existed.
      integer, parameter :: NP_EXP(12) = [9, 9, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1]
      integer :: i

      call set_prune("0")
      call build_water_dz(bas)
      do i = 1, 12
         call check(error, bas%sh_np(i), NP_EXP(i), "TRC_PRUNE_PRIMITIVES=0 keeps the BSE-literal count")
         if (allocated(error)) return
      end do
      call set_prune("1")   ! restore the default for the tests that follow
   end subroutine test_off_switch_is_literal

   subroutine run_rhf(prune_flag, energy)
      character(*), intent(in) :: prune_flag
      real(dp), intent(out) :: energy
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer :: nocc

      call set_prune(prune_flag)
      call build_water_dz(bas)
      call bas%to_device()
      nocc = 5   ! water, 10 electrons
      opts%conv_energy = 1.0e-12_dp
      opts%conv_density = 1.0e-10_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      call trc_scf_run(bas, nocc, nocc, opts, res)
      if (.not. res%converged) error stop 'check_prune_primitives: RHF did not converge'
      energy = res%energy
   end subroutine run_rhf

   !
   ! The elimination only ever subtracts a zero -- see prune_block's own
   ! header -- so the pruned and unpruned bases span exactly the same
   ! functions and RHF must land on the same stationary point to round-off,
   ! not merely to engineering tolerance.
   !
   subroutine test_energy_invariant(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e_on, e_off

      call run_rhf("1", e_on)
      call run_rhf("0", e_off)
      call set_prune("1")
      call check(error, abs(e_on - e_off) < 1.0e-10_dp, &
                 "E_RHF(pruned) - E_RHF(unpruned) is at round-off, not " &
                 //"an approximation")
   end subroutine test_energy_invariant

   !
   ! Same construction (force_cart=.true.) as test/scf_cphf.F90, whose own
   ! cphf_water_cc-pvdz.ref carries this system's pyscf RHF energy
   ! (test/pyscf_refs/generate.py's gen_cphf). Reusing it here, rather than
   ! a second generated reference, needs no new "# files:" entry: this
   ! check only reads the committed anchor, same as check_ao.F90 and
   ! check_mult.F90 read more than one reference each.
   !
   subroutine test_matches_anchor(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e_on, e_ref
      integer :: u

      call run_rhf("1", e_on)
      call ref_require('cphf_water_cc-pvdz.ref', u)
      call ref_read_real(u, 'E_SCF', e_ref)
      call ref_close(u)
      call check(error, abs(e_on - e_ref) < 1.0e-9_dp, "E_RHF(pruned) matches the pyscf anchor")
   end subroutine test_matches_anchor

   !
   ! prune_block's NaN guard: a column that shares EVERY one of its
   ! primitives with the block's pivots comes out all zero after the
   ! elimination loop, not merely smaller. A hand-written one-block s
   ! shell (element 2 -- arbitrary, this basis is never used for anything
   ! but this test) with exponents 1, 2, 3 and three columns: two
   ! single-primitive pivots (at exponents 1 and 2) and a third column
   ! that is a contraction over EXACTLY those same two primitives and
   ! nothing else (coefficients 0.5, 0.5, 0 over exponents 1, 2, 3).
   ! prune_block zeroes the third column at both pivot positions, which
   ! is its every nonzero primitive -- without the guard, the column
   ! would reach `normalise` all zero, which divides by sqrt(overlap) and
   ! hands back NaN coefficients with no error status at all.
   !
   subroutine build_allzero_column_bas(bas, error)
      type(trc_basis_t), intent(out) :: bas
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: SCRATCH = "check_prune_primitives_allzero_tmp.json"
      type(error_t) :: err
      integer, parameter :: ELEMENTS(1) = [2]
      real(dp), parameter :: AT_R(3, 1) = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])
      integer :: u, ios

      open (newunit=u, file=SCRATCH, status="replace", action="write")
      write (u, "(a)") '{'
      write (u, "(a)") ' "elements": {'
      write (u, "(a)") '  "2": {'
      write (u, "(a)") '   "electron_shells": ['
      write (u, "(a)") '    {'
      write (u, "(a)") '     "angular_momentum": [0],'
      write (u, "(a)") '     "exponents": ["1.0", "2.0", "3.0"],'
      write (u, "(a)") '     "coefficients": ['
      write (u, "(a)") '      ["1.0", "0.0", "0.0"],'
      write (u, "(a)") '      ["0.0", "1.0", "0.0"],'
      write (u, "(a)") '      ["0.5", "0.5", "0.0"]'
      write (u, "(a)") '     ]'
      write (u, "(a)") '    }'
      write (u, "(a)") '   ]'
      write (u, "(a)") '  }'
      write (u, "(a)") ' }'
      write (u, "(a)") '}'
      close (u)

      call trc_basis_from_json(SCRATCH, 1, ELEMENTS, AT_R, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(prune=.true.))
      open (newunit=u, file=SCRATCH, status="old", iostat=ios)
      if (ios == 0) close (u, status="delete")
      call check(error, .not. err%has_error(), "trc_basis_from_json failed on the contrived block")
   end subroutine build_allzero_column_bas

   subroutine test_allzero_column_restored(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas

      call build_allzero_column_bas(bas, error)
      if (allocated(error)) return
      call check(error, bas%nshell, 3, "the contrived block keeps all 3 columns")
      if (allocated(error)) return
      call check(error, .not. any(ieee_is_nan(bas%sh_c)), &
                "a column pruned to all zero must not come back NaN")
      if (.not. allocated(error)) then
         ! Restored to its two real primitives (exponents 1 and 2), not
         ! compacted to the 0 the elimination would otherwise leave it at.
         call check(error, bas%sh_np(3), 2, &
                   "the restored column keeps its 2 real primitives, not 0")
      end if
   end subroutine test_allzero_column_restored

end module check_prune_primitives_suite

program check_prune_primitives
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_prune_primitives_suite, only: collect_prune_primitives
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_prune_primitives, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_prune_primitives
