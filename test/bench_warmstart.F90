!
! Iteration counts: the plain density warm start (a finite-difference
! Hessian's own, `dguess` at the PREVIOUS geometry's converged density)
! against the orbital-projection warm start (`cmo_project`, trc_scf_context's
! `project_density`), at the SAME displaced geometry -- water/def2-SVP,
! HF and B3LYP, at a finite-difference Hessian's own h=0.005 Bohr and a
! larger, geometry-optimiser-sized step. Informational only (prints a
! table, asserts nothing). Measured result, both warm starts needing the
! same iteration count at every step: HF 10/13/15, B3LYP 8/10/12 (steps
! 0.005, 0.2, 1.0 Bohr in that order).
!
program bench_warmstart
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none

   real(dp), parameter :: HSMALL = 0.005_dp, HLARGE = 0.2_dp, HHUGE = 1.0_dp
   integer, allocatable :: z(:)
   real(dp), allocatable :: r0(:, :)
   integer :: nat, dev, nalpha, nbeta

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', nat, z, r0)
   nalpha = sum(z)/2; nbeta = nalpha

   print '(a)', 'func      step(Bohr)  iters(density warm start)  iters(projected warm start)'
   call run_func('hf', HSMALL)
   call run_func('hf', HLARGE)
   call run_func('hf', HHUGE)
   call run_func('b3lyp', HSMALL)
   call run_func('b3lyp', HLARGE)
   call run_func('b3lyp', HHUGE)

contains

   subroutine run_func(func, h)
      character(*), intent(in) :: func
      real(dp), intent(in) :: h
      type(trc_basis_t) :: b0, b1
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: ref, res_density, res_projected
      type(error_t) :: err
      real(dp), allocatable :: r1(:, :), d0(:, :, :)

      call trc_basis_from_json('basis_sets/def2-svp.json', nat, z, r0, b0, err)
      if (err%has_error()) then
         print '(a)', 'bench_warmstart: '//err%get_message(); stop 1
      end if
      call b0%to_device()
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iter = 100
      if (func /= 'hf') opts%functional = func
      call trc_scf_run(b0, nalpha, nbeta, opts, ref)
      if (.not. ref%converged) then
         print '(a)', 'bench_warmstart: '//trim(func)//': the reference SCF did not converge'; stop 1
      end if
      allocate (d0, source=ref%dmat)
      allocate (r1, source=r0)
      r1(1, 1) = r1(1, 1) + h

      call trc_basis_from_json('basis_sets/def2-svp.json', nat, z, r1, b1, err)
      if (err%has_error()) then
         print '(a)', 'bench_warmstart: '//err%get_message(); stop 1
      end if
      call b1%to_device()

      call trc_scf_run(b1, nalpha, nbeta, opts, res_density, dguess=d0)
      call trc_scf_run(b1, nalpha, nbeta, opts, res_projected, cmo_project=ref%cmo)

      print '(a6,a,f10.3,a,i0,a,i0)', func, '    ', h, repeat(' ', 20), res_density%iterations, &
         repeat(' ', 20), res_projected%iterations
      if (.not. res_density%converged .or. .not. res_projected%converged) then
         print '(a,l1,a,l1)', '  (converged: density=', res_density%converged, ' projected=', res_projected%converged
      end if

      call b0%release(); call b1%release()
      deallocate (d0, r1)
   end subroutine run_func

end program bench_warmstart
