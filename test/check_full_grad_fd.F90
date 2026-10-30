!
! The FULL analytic RHF/UHF nuclear gradient (trc_scf_grad%trc_scf_gradient:
! eri%gradient + trc_1e_gradient + nuclear repulsion) against a central
! finite difference of terco's own FULLY RECONVERGED SCF energy -- the
! Hellmann-Feynman/Pulay theorem's actual claim, extending
! check_rys_grad_fd.F90's technique (which only exercises the two-electron
! piece at a FIXED, non-SCF density) to the complete gradient at a
! converged one.
!
! RHF on water, UHF on the OH radical (a genuinely open-shell case, nalpha
! /= nbeta) -- both 6-31G*. Every displacement reconverges the SCF from
! scratch (GWH guess, tight conv_energy/conv_diis): the FD step is small
! enough (1e-3 bohr) that this is deterministic and costs nothing dramatic
! for these two-and-three-atom systems.
!
! Also checks translational invariance of the analytic gradient alone
! (sum over atoms ~ 0), the one property that holds exactly regardless of
! step size or SCF tolerance.
!
! NO INTERNAL PROCEDURES. `run_case`/`scf_energy_at` are genuine module
! procedures, not a `contains` block under the program, for the same
! reason app/main.F90's `trc_run` module gives for its own equivalent
! refactor: an internal procedure anywhere in the executable -- even one
! with no state in common with the GPU call at all -- corrupted a later
! collective Fock build on this exact NVHPC/OpenMPI toolchain,
! deterministically, regardless of optimisation level. Measured here too:
! as a `contains` block, this program's SECOND trc_scf_run (the UHF OH
! case, after water's reference + gradient + 18-point FD loop) silently
! failed to converge on GPU 0 -- not a crash, just wrong behaviour -- and
! the identical sequence, as plain module procedures, does not.
!
module check_full_grad_fd_cases
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none
   private
   public :: run_case

   real(dp), parameter :: H = 1.0e-3_dp
   real(dp), parameter, public :: RELTOL = 1.0e-5_dp
   !! Central difference at H = 1e-3 on a reconverged SCF energy: the
   !! truncation floor measured here is ~1e-7 absolute on components of
   !! order 1e-2 to 1e-1 -- a tight ABSOLUTE floor, not RELTOL's own
   !! bottleneck, which is why ABSTOL (below) is the gate that matters for
   !! the near-zero components translational invariance already covers.
   real(dp), parameter, public :: ABSTOL = 5.0e-7_dp
   real(dp), parameter, public :: INVTOL = 1.0e-10_dp

contains

   subroutine run_case(xyzfile, basis_json, nalpha, nbeta, label, any_fail)
      character(*), intent(in) :: xyzfile, basis_json, label
      integer, intent(in) :: nalpha, nbeta
      logical, intent(inout) :: any_fail

      integer :: nat, dev, a, d
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :), rtmp(:, :), grad(:, :)
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      real(dp) :: eplus, eminus, fd, worst_abs, worst_rel, denom
      real(dp) :: inv(3)

      dev = trc_bind_device(0)
      call read_xyz(xyzfile, nat, at_z, r0)
      call trc_basis_from_json(basis_json, nat, at_z, r0, bas, err)
      if (err%has_error()) then
         print '(a)', 'check_full_grad_fd: '//err%get_message(); stop 1
      end if
      call bas%to_device()

      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iteration = 200
      call trc_scf_run(bas, nalpha, nbeta, opts, res)
      if (.not. res%converged) then
         print '(a)', 'check_full_grad_fd: '//label//': reference SCF did not converge'
         any_fail = .true.; return
      end if

      allocate (grad(3, nat))
      call lp%init(bas%nao)
      call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, 1.0e-12_dp, grad)
      ! Released before `lp` goes out of scope: trc_linalg_t has no
      ! finalizer, and an unreleased one leaves its device workspace mapped
      ! after Fortran frees the host array, so a later allocation at that
      ! host address is "present" with the wrong size and device memory.
      call lp%release()
      call bas%release()

      inv = sum(grad, dim=2)
      if (maxval(abs(inv)) > INVTOL) any_fail = .true.

      allocate (rtmp(3, nat))
      worst_abs = 0.0_dp; worst_rel = 0.0_dp
      do a = 1, nat
         do d = 1, 3
            rtmp = r0; rtmp(d, a) = r0(d, a) + H
            eplus = scf_energy_at(basis_json, nat, at_z, rtmp, nalpha, nbeta, opts)
            rtmp = r0; rtmp(d, a) = r0(d, a) - H
            eminus = scf_energy_at(basis_json, nat, at_z, rtmp, nalpha, nbeta, opts)
            fd = (eplus - eminus)/(2.0_dp*H)
            denom = max(abs(grad(d, a)), abs(fd), 1.0e-10_dp)
            worst_abs = max(worst_abs, abs(grad(d, a) - fd))
            worst_rel = max(worst_rel, abs(grad(d, a) - fd)/denom)
         end do
      end do

      print '(a,a,a,es10.2,a,es10.2,a,3es10.2)', '  ', label, ': worst abs ', worst_abs, &
         '  worst rel ', worst_rel, '  translational-invariance residual: ', inv
      if (worst_abs > ABSTOL .and. worst_rel > RELTOL) any_fail = .true.

      deallocate (at_z, r0, rtmp, grad)
   end subroutine run_case

   real(dp) function scf_energy_at(basis_json, nat, at_z, r, nalpha, nbeta, opts) result(e)
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nalpha, nbeta
      real(dp), intent(in) :: r(3, nat)
      type(trc_scf_options_t), intent(in) :: opts
      type(trc_basis_t) :: b2
      type(trc_scf_result_t) :: r2
      type(error_t) :: err
      call trc_basis_from_json(basis_json, nat, at_z, r, b2, err)
      if (err%has_error()) then
         print '(a)', 'check_full_grad_fd: '//err%get_message(); stop 1
      end if
      call b2%to_device()
      call trc_scf_run(b2, nalpha, nbeta, opts, r2)
      e = r2%energy
      call b2%release()
   end function scf_energy_at

end module check_full_grad_fd_cases

program check_full_grad_fd
   use check_full_grad_fd_cases, only: run_case
   implicit none
   logical :: any_fail
   any_fail = .false.

   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'RHF water/6-31G*', any_fail)
   call run_case('oh.xyz', 'basis_sets/6-31gs.json', 5, 4, 'UHF OH/6-31G*', any_fail)
   ! SPHERICAL cases: cc-pVDZ is tagged gto_spherical in BSE, so
   ! trc_basis_from_json comes back with `bas%spherical` true here, and
   ! `run_case`'s own trc_scf_gradient call takes the spherical fold in
   ! trc_scf_grad.F90 (D and W folded nao_sph -> nao, T D T^T) -- these two
   ! are this program's own validation of that fold against a finite
   ! difference of terco's own reconverged energy, independent of pyscf.
   call run_case('water.xyz', 'basis_sets/cc-pvdz.json', 5, 5, 'RHF water/cc-pVDZ (spherical)', any_fail)
   call run_case('oh.xyz', 'basis_sets/cc-pvdz.json', 5, 4, 'UHF OH/cc-pVDZ (spherical)', any_fail)

   if (any_fail) then
      print '(a)', '  RESULT: FAIL'; stop 1
   end if
   print '(a)', '  RESULT: PASS'
end program check_full_grad_fd
