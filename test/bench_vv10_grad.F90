!
! Where the VV10 GRADIENT's time goes, next to the rest of the wB97X-V
! gradient, on a real molecule in a real basis -- bench_vv10.F90's own
! idea, applied to trc_scf_gradient instead of the energy path's
! trc_vv10_rks.
!
!   bench_vv10_grad <xyz> <basis.json> <nalpha> <nbeta>
!
! Converges wB97X-V once, then times trc_scf_gradient TWICE at the SAME
! converged orbitals: once with the functional exactly as converged
! (func%nlc_b/nlc_c nonzero -- the full gradient, semilocal GGA + erf
! long-range exchange + VV10), and once with a LOCAL COPY of that same
! functional with nlc_b/nlc_c zeroed (everything else identical: same
! kid/coef semilocal kernel list, same omega/exx_lr) -- "the rest" of the
! gradient at the exact same cost, minus VV10's own nonlocal term. The
! difference between the two wall times is VV10's own share.
!
! A timing is not a pass, so this is built and not registered as a test
! (test/CMakeLists.txt, same as bench_vv10).
!
program bench_vv10_grad
   use pic_types, only: dp
   use trc_wall_clock, only: trc_wall
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none

   character(len=256) :: xyzfile, bas_json, arg
   integer :: nat, nalpha, nbeta, dev
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: r(:, :), grad_full(:, :), grad_novv10(:, :)
   type(trc_basis_t) :: bas
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: res
   type(trc_xc_functional_t) :: func, func_novv10
   type(trc_linalg_t) :: lp
   type(error_t) :: err
   real(dp) :: t0, t1, t_full, t_novv10

   if (command_argument_count() < 4) then
      print '(a)', 'usage: bench_vv10_grad <xyz> <basis.json> <nalpha> <nbeta>'; stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, bas_json)
   call get_command_argument(3, arg); read (arg, *) nalpha
   call get_command_argument(4, arg); read (arg, *) nbeta

   if (.not. xc_functional_by_name('wb97x-v', func)) then
      print '(a)', 'bench_vv10_grad: wb97x-v not available'; stop 1
   end if

   dev = trc_bind_device(0)
   call read_xyz(trim(xyzfile), nat, at_z, r)
   call trc_basis_from_json(trim(bas_json), nat, at_z, r, bas, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', 'bench_vv10_grad: '//err%get_message(); stop 1
   end if
   call bas%to_device()

   opts%functional = 'wb97x-v'
   opts%grid_level = 3
   opts%conv_energy = 1.0e-10_dp; opts%conv_diis = 1.0e-8_dp; opts%max_iter = 100
   t0 = trc_wall()
   call trc_scf_run(bas, nalpha, nbeta, opts, res)
   t1 = trc_wall()
   if (.not. res%converged) then
      print '(a)', 'bench_vv10_grad: SCF did not converge'; stop 1
   end if
   print '(a,i0,a,i0,a,f8.3,a)', 'bench_vv10_grad: ', nat, ' atoms, ', bas%nao, ' AOs; SCF ', t1 - t0, ' s'

   ! "The rest": the same functional with VV10's two coefficients zeroed,
   ! so trc_xc_grad's semilocal GGA term and trc_scf_grad.F90's own
   ! long-range exchange pass run exactly as they did for the real
   ! functional, and only vv10_gradient_add's own call is skipped.
   func_novv10 = func
   func_novv10%nlc_b = 0.0_dp
   func_novv10%nlc_c = 0.0_dp

   allocate (grad_full(3, nat), grad_novv10(3, nat))
   call lp%init(bas%nao)

   t0 = trc_wall()
   call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, 1.0e-12_dp, grad_full, &
                         func=func, grid_level=opts%grid_level, grid_max_pts=opts%grid_max_pts, &
                         rho_tol=opts%rho_tol, nlc_grid_level=opts%nlc_grid_level)
   t1 = trc_wall()
   t_full = t1 - t0

   t0 = trc_wall()
   call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, 1.0e-12_dp, grad_novv10, &
                         func=func_novv10, grid_level=opts%grid_level, grid_max_pts=opts%grid_max_pts, &
                         rho_tol=opts%rho_tol, nlc_grid_level=opts%nlc_grid_level)
   t1 = trc_wall()
   t_novv10 = t1 - t0

   call lp%release()
   call bas%release()

   print '(a,f10.3,a)', '  full gradient (semilocal GGA + erf LR-x + VV10): ', t_full, ' s'
   print '(a,f10.3,a)', '  same, func%nlc_b/nlc_c zeroed (no VV10 call):    ', t_novv10, ' s'
   print '(a,f10.3,a,f6.2,a)', '  VV10 gradient piece, by difference:              ', t_full - t_novv10, ' s (', &
      100.0_dp*(t_full - t_novv10)/t_full, '% of the full gradient)'

   deallocate (at_z, r, grad_full, grad_novv10)
end program bench_vv10_grad
