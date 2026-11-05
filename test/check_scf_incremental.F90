!
! The incremental Fock build against full builds, run for run: the same
! molecule, basis and thresholds with `incremental` off and on must reach
! the same energy (to ETOL, below) in at most one more iteration, and the
! incremental run must actually have built from the density difference
! for all but a few iterations (MAX_FULL, below). Restricted and
! unrestricted, Hartree-Fock, a range-separated hybrid (the erf exchange
! increments too) and the continuum (rebuilt in full every iteration, on
! the full density). And a run that converges on its last allowed
! iteration is reported converged (test_last_iteration).
!
! Before trc_scf_incremental an open shell never incremented, and a closed
! shell waited for |FDS - SDF| < 1e-2 before its first increment.
!
module check_scf_incremental_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_scf_incremental_tests

   real(dp), parameter :: ETOL = 1.0e-10_dp
   !! Incremental against full builds. Measured, all five cases: at most
   !! 2.8e-13 Eh under gfortran 15 and 5.7e-13 under ifx 2026.2 (conv_energy
   !! 1e-10), so 1e-10 is the convergence test's own scale with two decades
   !! to spare.
   integer, parameter :: MAX_FULL = 4
   !! Full builds an incremental run may take, at least: the first, the one
   !! leaving the loose early phase (trc_scf_incremental's rule 0), the one
   !! predicted to converge (rule 2), and -- when that prediction narrowly
   !! misses, its energy change just above conv_energy -- the next build,
   !! full again, that converges. ifx 2025.0 took exactly that path on the
   !! UKS case: 4 of 15 builds full. A long run (an open shell's damped
   !! start) may rebuild on a stall more often, up to a quarter of its
   !! iterations.

contains

   subroutine collect_scf_incremental_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF water trimer/6-31G*", test_rhf), &
                  new_unittest("UHF OH/6-31G*", test_uhf), &
                  new_unittest("UKS wB97X OH/6-31G* (range-separated)", test_uks_rsh), &
                  new_unittest("RHF water/6-31G* C-PCM", test_rhf_pcm), &
                  new_unittest("converging on the last allowed iteration", test_last_iteration) &
                  ]
   end subroutine collect_scf_incremental_tests

   subroutine test_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_options_t) :: opts
      call compare(error, "w3.xyz", 15, 15, opts, "RHF w3/6-31G*")
   end subroutine test_rhf

   subroutine test_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_options_t) :: opts
      call compare(error, "oh.xyz", 5, 4, opts, "UHF OH/6-31G*")
   end subroutine test_uhf

   subroutine test_uks_rsh(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_options_t) :: opts
      opts%functional = "wb97x"
      call compare(error, "oh.xyz", 5, 4, opts, "UKS wB97X OH/6-31G*")
   end subroutine test_uks_rsh

   subroutine test_rhf_pcm(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_options_t) :: opts
      opts%pcm_model = "cpcm"
      opts%pcm_eps = 78.3553_dp
      call compare(error, "water.xyz", 5, 5, opts, "RHF water/6-31G* C-PCM")
   end subroutine test_rhf_pcm

   !
   ! A run that converges on its LAST allowed iteration must say so. With
   ! the prediction off (`incr_predict = 0`), the convergence test is first
   ! met on an increment, refused there, and retaken on the next build,
   ! which is full: the unbounded run converges at N with N - 1 an
   ! increment that met the test. Capped at N - 1 iterations, that last
   ! iteration has to be a full build (trc_scf_incremental's
   ! `next_build`), or the run ends unconverged although it met the test.
   !
   subroutine test_last_iteration(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: free, capped
      integer :: nat, dev
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :)
      type(trc_basis_t) :: bas
      type(error_t) :: err

      dev = trc_bind_device(0)
      call read_xyz("w3.xyz", nat, at_z, r0)
      call trc_basis_from_json("basis_sets/6-31gs.json", nat, at_z, r0, bas, err)
      call check(error, .not. err%has_error(), "last iteration: basis build failed: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-6_dp
      opts%max_iteration = 100
      opts%incremental = .true.
      opts%incr_predict = 0.0_dp
      call trc_scf_run(bas, 15, 15, opts, free)
      opts%max_iteration = free%iterations - 1
      call trc_scf_run(bas, 15, 15, opts, capped)
      call bas%release()
      print '(a,i4,a,l2,a,es10.2)', "last iteration: unbounded run converged in", free%iterations, &
         "; capped one below it, converged", capped%converged, ", dE ", capped%energy - free%energy

      call check(error, free%converged, "last iteration: the unbounded run did not converge")
      if (allocated(error)) return
      call check(error, capped%converged .and. capped%iterations == free%iterations - 1, &
                 "last iteration: capped at "//to_string(free%iterations - 1)//" iterations, the run "// &
                 "stopped after "//to_string(capped%iterations)//", converged "// &
                 trim(merge("yes", "no ", capped%converged)))
      if (allocated(error)) return
      call check(error, abs(capped%energy - free%energy) <= ETOL, "last iteration: capped E "// &
                 to_string(capped%energy, "f20.12")//" against "//to_string(free%energy, "f20.12"))
   end subroutine test_last_iteration

   subroutine compare(error, xyzfile, nalpha, nbeta, opts_in, label)
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: xyzfile, label
      integer, intent(in) :: nalpha, nbeta
      type(trc_scf_options_t), intent(in) :: opts_in
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: full, incr
      integer :: nat, dev, allowed
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :)
      type(trc_basis_t) :: bas
      type(error_t) :: err

      dev = trc_bind_device(0)
      call read_xyz(xyzfile, nat, at_z, r0)
      call trc_basis_from_json("basis_sets/6-31gs.json", nat, at_z, r0, bas, err)
      call check(error, .not. err%has_error(), label//": basis build failed: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()

      opts = opts_in
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-6_dp
      opts%max_iteration = 100
      opts%incremental = .false.
      call trc_scf_run(bas, nalpha, nbeta, opts, full)
      opts%incremental = .true.
      call trc_scf_run(bas, nalpha, nbeta, opts, incr)
      call bas%release()
      ! For the record (and for setting ETOL): the run's numbers, one line.
      print '(2a,i4,a,i4,a,i3,a,es10.2)', label, ": full ", full%iterations, " iterations, incremental ", &
         incr%iterations, " with ", incr%iterations - incr%n_increments, " full builds, dE ", incr%energy - full%energy

      call check(error, full%converged .and. incr%converged, label//": both SCFs must converge (full "// &
                 trim(merge("yes", "no ", full%converged))//", incremental "// &
                 trim(merge("yes", "no ", incr%converged))//")")
      if (allocated(error)) return
      call check(error, full%n_increments == 0, label//": incremental off took "// &
                 to_string(full%n_increments)//" increments")
      if (allocated(error)) return
      call check(error, abs(incr%energy - full%energy) <= ETOL, label//": incremental E "// &
                 to_string(incr%energy, "f20.12")//" against full "//to_string(full%energy, "f20.12")// &
                 ", diff "//to_string(incr%energy - full%energy, "es10.2"))
      if (allocated(error)) return
      call check(error, incr%iterations <= full%iterations + 1, label//": incremental took "// &
                 to_string(incr%iterations)//" iterations against "//to_string(full%iterations))
      if (allocated(error)) return
      allowed = max(MAX_FULL, incr%iterations/4)
      call check(error, incr%iterations - incr%n_increments <= allowed, label//": "// &
                 to_string(incr%iterations - incr%n_increments)//" of "//to_string(incr%iterations)// &
                 " builds were full, at most "//to_string(allowed)//" expected")
   end subroutine compare

end module check_scf_incremental_suite

program check_scf_incremental
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_scf_incremental_suite, only: collect_scf_incremental_tests
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_incremental_tests, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_scf_incremental
