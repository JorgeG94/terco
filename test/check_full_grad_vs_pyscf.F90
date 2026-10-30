!
! The FULL analytic RHF/UHF nuclear gradient (trc_scf_grad%trc_scf_gradient:
! one-electron + eri%gradient's two-electron part + nuclear repulsion)
! against pyscf's own `Gradients` class AT A CONVERGED density, on
! SPHERICAL cc-pVDZ bases -- water/RHF, formamide/RHF, OH/UHF.
!
! WHY THIS PROGRAM, ON TOP OF check_full_grad_fd.F90 AND
! check_1e_grad_vs_pyscf.F90. Neither of those exercises trc_scf_gradient's
! own spherical fold (D and W folded nao_sph -> the Cartesian nao the
! kernels need, T D T^T -- src/scf/trc_scf_grad.F90's own header):
! check_full_grad_fd.F90 only runs 6-31G* (Cartesian, no fold at all), and
! check_1e_grad_vs_pyscf.F90 hands trc_1e_gradient a synthetic density
! directly at the Cartesian nao, bypassing the SCF (and the fold) entirely
! even on a basis BSE tags spherical. This program is the one path that
! builds a converged SPHERICAL SCF, then asks trc_scf_gradient for its
! gradient -- the fold runs for real here.
!
! FIXED INITIAL DENSITY. See test/pyscf_refs/generate.py's gen_full_grad
! for the full reasoning (its own docstring): water and formamide need
! none (a closed-shell ground state has no nearby competing stationary
! point); OH does, the same recipe rks.ref's own OH cases use, because its
! unpaired electron sits in a genuinely degenerate pair of pi orbitals.
! The committed DM0 (OH only) seeds `trc_scf_run`'s own `dguess`, so both
! sides converge toward the SAME orientation.
!
module check_full_grad_vs_pyscf_cases
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix, ref_read_matrix3, ref_report
   implicit none
   private
   public :: run_case

   real(dp), parameter, public :: ABSTOL = 1.0e-6_dp
   !! Analytic-vs-analytic at two independently converged (but, for water/
   !! formamide, not bit-identical) SCFs -- tighter than check_full_grad_fd's
   !! FD-truncation-limited 5e-7 is not expected; this is the same order of
   !! magnitude, on top of whatever residual the two SCFs' own, slightly
   !! different convergence paths leave in the density.

contains

   subroutine run_case(u, xyzfile, name, nalpha, nbeta, nspin, any_fail)
      integer, intent(in) :: u
      character(*), intent(in) :: xyzfile, name
      integer, intent(in) :: nalpha, nbeta, nspin
      logical, intent(inout) :: any_fail

      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :), grad(:, :), grad_ref(:, :), dm0(:, :, :)
      character(len=64) :: name_ref
      integer :: nat, nao, nao_ref, nat_ref, nspin_ref, dev, nfail
      real(dp) :: e_ref, gnorm_ref, worst

      dev = trc_bind_device(0)
      call read_xyz(xyzfile, nat, at_z, r0)
      call trc_basis_from_json('basis_sets/cc-pvdz.json', nat, at_z, r0, bas, err)
      if (err%has_error()) then
         print '(a)', 'check_full_grad_vs_pyscf: '//err%get_message(); stop 1
      end if
      if (.not. bas%spherical) then
         print '(a)', 'check_full_grad_vs_pyscf: '//trim(name)// &
            ': cc-pVDZ did not come back spherical -- nothing to fold here'
         any_fail = .true.; return
      end if
      call bas%to_device()
      nao = bas%nao_sph

      ! --- the committed reference, strictly in write order (see
      ! test/pyscf_refs/generate.py's gen_full_grad and pyscf_ref_io.F90's
      ! own header: tags are read sequentially, not by seeking a name) ---
      call ref_read_str(u, 'NAME_'//trim(name), name_ref)
      call ref_read_int(u, 'NSPIN_'//trim(name), nspin_ref)
      call ref_read_int(u, 'NAO_'//trim(name), nao_ref)
      call ref_read_int(u, 'NAT_'//trim(name), nat_ref)
      if (nspin_ref /= nspin .or. nao_ref /= nao .or. nat_ref /= nat) then
         print '(a)', 'check_full_grad_vs_pyscf: '//trim(name)//': reference shape does not match this build'
         any_fail = .true.; return
      end if
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iteration = 200
      if (nspin == 2) then
         allocate (dm0(nao, nao, 2))
         call ref_read_matrix3(u, 'DM0_'//trim(name), nao, nao, 2, dm0)
         call trc_scf_run(bas, nalpha, nbeta, opts, res, dguess=dm0)
      else
         call trc_scf_run(bas, nalpha, nbeta, opts, res)
      end if
      call ref_read_real(u, 'E_'//trim(name), e_ref)
      call ref_read_real(u, 'GRAD_NORM_'//trim(name), gnorm_ref)
      allocate (grad_ref(3, nat))
      call ref_read_matrix(u, 'GRAD_'//trim(name), 3, nat, grad_ref)

      if (.not. res%converged) then
         print '(a)', 'check_full_grad_vs_pyscf: '//trim(name)//': terco SCF did not converge'
         any_fail = .true.; return
      end if

      allocate (grad(3, nat))
      call lp%init(nao)
      call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, 1.0e-12_dp, grad)
      ! Released before `lp` goes out of scope -- trc_linalg_t has no
      ! finalizer (see check_full_grad_fd.F90's own comment on this).
      call lp%release()
      call bas%release()

      nfail = 0
      worst = maxval(abs(grad - grad_ref))
      call ref_report(nfail, worst < ABSTOL, trim(name)//': gradient vs pyscf, worst|diff|', worst)
      print '(a,a,a,es12.4,a,es12.4,a,es10.2)', '  ', trim(name), ': |E(terco)-E(pyscf)| ', &
         abs(res%energy - e_ref), '  pyscf |g| ', gnorm_ref, '  worst|diff| ', worst
      if (nfail > 0) any_fail = .true.

      deallocate (at_z, r0, grad, grad_ref)
      if (allocated(dm0)) deallocate (dm0)
   end subroutine run_case

end module check_full_grad_vs_pyscf_cases

program check_full_grad_vs_pyscf
   use check_full_grad_vs_pyscf_cases, only: run_case
   use pyscf_ref_io, only: ref_require, ref_close
   implicit none
   integer :: u
   logical :: any_fail
   any_fail = .false.

   call ref_require('full_grad_ccpvdz.ref', u)
   call run_case(u, 'water.xyz', 'water_ccpvdz', 5, 5, 1, any_fail)
   call run_case(u, 'formamide.xyz', 'formamide_ccpvdz', 12, 12, 1, any_fail)
   call run_case(u, 'oh.xyz', 'oh_ccpvdz', 5, 4, 2, any_fail)
   call ref_close(u)

   if (any_fail) then
      print '(a)', '  RESULT: FAIL'; stop 1
   end if
   print '(a)', '  RESULT: PASS'
end program check_full_grad_vs_pyscf
