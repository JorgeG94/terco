program time_tddft_grad
   !! End-to-end timing of TDA (CIS for Hartree-Fock) excited-state gradients:
   !! SCF, the TDA solve and trc_tddft_gradient_rks, each timed on its own,
   !! then the gradient norms of every root with ten significant digits so that two
   !! builds (host against device, or two roots_per_sweep) compare numerically.
   !!
   !! Usage: time_tddft_grad basis.json molecule.xyz functional nroots roots_per_sweep
   !! Run from the test/ source directory (the basis path is relative to it); the binary
   !! is built under the build tree's test/standalone. The functional is `hf` or any name
   !! xc_functional_by_name takes that the gradient accepts (LDA/GGA, global
   !! hybrids). The basis is always loaded Cartesian (force_cart), since the
   !! gradient refuses spherical shells. Same setup path as check_tddft_grad, with looser SCF thresholds.
   !!
   !! Environment: TDA_CONV_TOL overrides the Davidson residual tolerance (default
   !! 1e-10; 1e-6 is plenty for a timing run). The TDA solve runs verbose, printing one
   !! line per Davidson iteration, then the iteration and matvec counts; max_iteration
   !! is 200 so a hard system is not cut off mid-timing.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use pic_timer, only: timer_type
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_tddft, only: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_solve_rks, trc_tddft_reference_t
   use trc_tddft_grad, only: trc_tddft_grad_options_t, trc_tddft_gradient_rks
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_constants, only: HARTREE_TO_EV
   use trc_error, only: error_t
   implicit none

   real(dp), parameter :: ERI_THRESH = 1.0e-14_dp
   real(dp), parameter :: CONV_TOL = 1.0e-10_dp
   real(dp), parameter :: CPKS_TOL = 1.0e-11_dp

   character(len=256) :: bas_path, xyz_path, func_name
   character(len=32) :: arg
   integer :: nroots, rps, natm, nocc, dev, k
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   logical :: dft, ok
   type(error_t) :: err
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   type(trc_eri_t) :: eri
   type(trc_linalg_t) :: la
   type(trc_xc_grid_t) :: xg
   type(trc_xc_functional_t) :: func
   type(trc_scf_options_t) :: sopts
   type(trc_scf_result_t) :: scf
   type(trc_tddft_options_t) :: topts
   character(len=32) :: envbuf
   integer :: envstat
   type(trc_tddft_result_t) :: td
   type(trc_tddft_reference_t) :: ref
   type(trc_tddft_grad_options_t) :: gopts
   type(timer_type) :: tm
   real(dp) :: t_setup, t_scf, t_td, t_grad, g, rms, gmax

   if (command_argument_count() < 5) then
      call global_logger%error("usage: time_tddft_grad basis.json molecule.xyz functional nroots roots_per_sweep")
      error stop 2
   end if
   call get_command_argument(1, bas_path)
   call get_command_argument(2, xyz_path)
   call get_command_argument(3, func_name)
   call get_command_argument(4, arg)
   read (arg, *) nroots
   call get_command_argument(5, arg)
   read (arg, *) rps

   dev = trc_bind_device(0)
   dft = trim(func_name) /= 'hf'
   if (dft) then
      ok = xc_functional_by_name(trim(func_name), func)
      if (.not. ok) then
         call global_logger%error("unrecognised functional '"//trim(func_name)//"'")
         error stop 2
      end if
   end if

   call tm%start()
   call read_xyz(trim(xyz_path), natm, zint, at_r)
   call trc_basis_from_json(trim(bas_path), natm, zint, at_r, bas, err, force_cart=.true.)
   if (err%has_error()) then
      call global_logger%error(err%get_message())
      error stop 1
   end if
   call bas%to_device()
   call pl%build(bas, ERI_THRESH)
   call pl%to_device()
   call eri%build(bas, ERI_THRESH)
   call la%init(bas%nao)
   nocc = nint(sum(real(zint, dp)))/2
   call tm%stop()
   t_setup = tm%get_elapsed_time()
   call global_logger%info("time_tddft_grad: "//trim(xyz_path)//", "//trim(bas_path)//", "//trim(func_name)//", "// &
                           to_string(natm)//" atoms, "//to_string(bas%nao)//" basis functions, "//to_string(nocc)// &
                           " occupied, "//to_string(nroots)//" roots, roots_per_sweep "//to_string(rps))

   ! check_tddft_grad's 1e-12/1e-10/1e-11 sit at the arithmetic's floor and the
   ! SCF ran out of iterations on 6-31G* H2O2 (gfortran + MKL); a timing run
   ! does not need them, and a large molecule could never meet them.
   sopts%conv_energy = 1.0e-10_dp
   sopts%conv_density = 1.0e-8_dp
   sopts%conv_diis = 1.0e-9_dp
   sopts%max_iteration = 200
   sopts%grid_level = 3
   if (dft) sopts%functional = trim(func_name)

   call tm%start()
   call trc_scf_run(bas, nocc, nocc, sopts, scf)
   call tm%stop()
   t_scf = tm%get_elapsed_time()
   if (.not. scf%converged) then
      call global_logger%error("SCF did not converge")
      call release_all()
      error stop 1
   end if
   call global_logger%info("SCF energy "//to_string(scf%energy, "f20.10")//" Eh in "//to_string(scf%iterations)// &
                           " iterations")

   if (dft) then
      call xg%build(size(scf%grid_w), scf%grid_r, scf%grid_w, bas, max_pts=sopts%grid_max_pts)
      call xg%to_device()
      call ref%set(nocc, scf%cmo(:, :, 1), scf%eps(:, 1), scf%dmat(:, :, 1), func)
   else
      call ref%set(nocc, scf%cmo(:, :, 1), scf%eps(:, 1), scf%dmat(:, :, 1))
   end if

   topts%nroots = nroots
   topts%conv_tol = CONV_TOL
   topts%verbose = .true.
   topts%max_iteration = 200
   call get_environment_variable("TDA_CONV_TOL", envbuf, status=envstat)
   if (envstat == 0 .and. len_trim(envbuf) > 0) read (envbuf, *) topts%conv_tol
   call tm%start()
   call trc_tddft_solve_rks(bas, pl, eri, la, ref, topts, td, xg)
   call tm%stop()
   t_td = tm%get_elapsed_time()
   call global_logger%info("TDA: conv_tol "//to_string(topts%conv_tol, "es9.2")//", "//to_string(td%iterations)// &
                           " Davidson iterations, "//to_string(td%nmatvec)//" matvecs")
   if (.not. td%converged) then
      call global_logger%error("TDA solve failed: "//td%err%get_message())
      call release_all()
      error stop 1
   end if
   do k = 1, nroots
      call global_logger%info("root "//to_string(k)//": omega "//to_string(td%omega(k)*HARTREE_TO_EV, "f12.6")//" eV, "// &
                              to_string(td%omega(k), "f18.12")//" Eh")
   end do

   gopts%roots_per_sweep = rps
   gopts%cpks_tol = CPKS_TOL
   gopts%grid_response = .false.
   gopts%verbose = .true.
   call tm%start()
   call trc_tddft_gradient_rks(bas, pl, eri, la, ref, td, gopts, xg)
   call tm%stop()
   t_grad = tm%get_elapsed_time()
   if (gopts%err%has_error()) then
      call global_logger%error("gradient failed: "//gopts%err%get_message())
      call release_all()
      error stop 1
   end if

   do k = 1, nroots
      g = sum(gopts%grad(:, :, k)**2)
      rms = sqrt(g/real(size(gopts%grad(:, :, k)), dp))
      gmax = maxval(abs(gopts%grad(:, :, k)))
      call global_logger%info("root "//to_string(k)//" gradient: rms "//to_string(rms, "es18.10")//" max|g| "// &
                              to_string(gmax, "es18.10")//" Eh/bohr")
   end do
   call global_logger%info("sweeps "//to_string(gopts%nsweep)//", CPKS iterations max "// &
                           to_string(maxval(gopts%cpks_iterations)))
   call global_logger%info("wall times (s): setup "//to_string(t_setup, "f10.3")//"  SCF "//to_string(t_scf, "f10.3")// &
                           "  TDA "//to_string(t_td, "f10.3")//"  gradient "//to_string(t_grad, "f10.3"))
   call release_all()

contains

   subroutine release_all()
      !! Every object built above, on every path out.
      call eri%release()
      call la%release()
      if (dft) call xg%release()
      call pl%release()
      call bas%release()
   end subroutine release_all

end program time_tddft_grad
