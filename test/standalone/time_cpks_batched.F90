program time_cpks_batched
   !! Wall time of the coupled-perturbed solve for 1, 2, 4 and 8 right-hand
   !! sides, batched (opts%max_batch = 8: one Fock and f_xc pass per CG
   !! iteration for every column still iterating) against one column at a
   !! time (max_batch = 1, the pre-batching behaviour), interleaved, the
   !! minimum of nrep repeats. Not a ctest: it asserts nothing, and is built
   !! as `time_cpks_batched`.
   !!
   !!   time_cpks_batched <basis.json> <molecule.xyz> [functional=hf] [nrep=3]
   !!
   !! Run it from test/ (paths are relative). The right-hand sides are the
   !! occupied-virtual blocks of the first multipole matrices, a dipole-like
   !! set of unlike difficulty. "passes" counts the Fock builds: the sum of the
   !! columns' iterations one at a time, the largest of them batched (one
   !! group); both are read from res%iterations.
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
   use trc_integrals, only: trc_multipoles, TRC_NMULT
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_cphf_solve_rhf
   use trc_cpks, only: trc_cpks_solve_rks
   use trc_xc, only: trc_xc_fxc_cache_t, trc_xc_fxc_rks, trc_xc_fxc_rks_apply_many
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none

   integer, parameter :: NCOLS(4) = [1, 2, 4, 8]
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   type(trc_eri_t) :: eri
   type(trc_linalg_t) :: la
   type(trc_xc_grid_t) :: xg
   type(trc_xc_functional_t) :: func
   type(trc_scf_options_t) :: sopts
   type(trc_scf_result_t) :: sres
   type(trc_cphf_options_t) :: opts
   type(trc_cphf_result_t) :: res
   type(error_t) :: err
   type(timer_type) :: tm
   character(len=256) :: basis_json, xyz, funcname, arg
   character(len=512) :: line
   integer :: natm, nocc, nvirt, nao, dev, k, rep, nrep, mode, ncol, ios
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :), mmat(:, :, :), cmo(:, :), eps(:), dmat(:, :), tmp(:, :), bvec(:, :, :), u(:, :, :)
   real(dp) :: t(2), tbest(2), udiff
   integer :: passes(2)
   real(dp), allocatable :: ubatch(:, :, :)
   logical :: ks
   type(trc_xc_fxc_cache_t) :: cache
   real(dp) :: t_setup, t_apply1, t_apply8, t_once, t_each8
   real(dp), allocatable :: dens(:, :, :), vr(:, :, :)

   if (command_argument_count() < 2) then
      call global_logger%error("usage: time_cpks_batched <basis.json> <molecule.xyz> [functional=hf] [nrep=3]")
      error stop 1
   end if
   call get_command_argument(1, basis_json)
   call get_command_argument(2, xyz)
   funcname = "hf"
   if (command_argument_count() >= 3) call get_command_argument(3, funcname)
   nrep = 3
   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg)
      read (arg, *, iostat=ios) nrep
   end if
   ks = trim(funcname) /= "hf"

   dev = trc_bind_device(0)
   call read_xyz(trim(xyz), natm, zint, at_r)
   call trc_basis_from_json(trim(basis_json), natm, zint, at_r, bas, err, force_cart=.true.)
   if (err%has_error()) error stop 1
   call bas%to_device()
   call pl%build(bas, 1.0e-12_dp)
   call pl%to_device()
   nocc = nint(sum(real(zint, dp)))/2
   nao = bas%nao
   nvirt = nao - nocc
   sopts%conv_energy = 1.0e-10_dp
   sopts%conv_density = 1.0e-8_dp
   sopts%conv_diis = 1.0e-8_dp
   sopts%max_iteration = 100
   if (ks) then
      sopts%functional = trim(funcname)
      sopts%grid_level = 3
   end if
   call trc_scf_run(bas, nocc, nocc, sopts, sres)
   if (.not. sres%converged) error stop 2
   if (ks) then
      if (.not. xc_functional_by_name(trim(funcname), func)) error stop 3
      call xg%build(size(sres%grid_w), sres%grid_r, sres%grid_w, bas, max_pts=sopts%grid_max_pts)
      call xg%to_device()
   end if
   call eri%build(bas, 1.0e-12_dp)
   call la%init(nao)

   allocate (cmo(nao, nao), eps(nao), dmat(nao, nao), mmat(nao, nao, TRC_NMULT), tmp(nao, nocc))
   cmo = sres%cmo(:, :, 1)
   eps = sres%eps(:, 1)
   dmat = sres%dmat(:, :, 1)
   call trc_multipoles(bas, pl, [0.0_dp, 0.0_dp, 0.0_dp], mmat)
   allocate (bvec(nvirt, nocc, 8), u(nvirt, nocc, 8), ubatch(nvirt, nocc, 8))
   !$acc enter data copyin(cmo, mmat) create(tmp, bvec, u)
   do k = 1, 8
      call la%gemm_mat(mmat(:, :, k + 1), cmo(:, 1:nocc), tmp)
      call la%gemm_mat(cmo(:, nocc + 1:nao), tmp, bvec(:, :, k), transa="T")
   end do
   !$acc update self(bvec)

   if (ks) then
      ! The pieces of a KS matvec, alone: the one-off cache of the grid
      ! collocation the batched solve builds, and a response pass over 1 and 8
      ! densities against it (after a warm-up).
      allocate (dens(nao, nao, 8), vr(nao, nao, 8))
      dens = 0.0_dp
      do k = 1, 8
         dens(:, :, k) = dmat*(0.01_dp*real(k, dp))
      end do
      !$acc enter data copyin(dens) create(vr)
      call cache%setup(bas, xg, func, dmat, la)
      call cache%release()
      call tm%start()
      call cache%setup(bas, xg, func, dmat, la)
      t_setup = tm%get_elapsed_time()
      ! Each time is the minimum of five calls, after a warm-up: the card is
      ! often shared and a single call is noisy.
      !$acc enter data copyin(dmat)
      call trc_xc_fxc_rks_apply_many(bas, xg, cache, 8, dens, vr, la=la)
      call trc_xc_fxc_rks(bas, xg, func, dmat, dens(:, :, 1), vr(:, :, 1), la=la)
      t_apply1 = huge(1.0_dp)
      t_apply8 = huge(1.0_dp)
      t_once = huge(1.0_dp)
      t_each8 = huge(1.0_dp)
      do rep = 1, 5
         call tm%start()
         call trc_xc_fxc_rks_apply_many(bas, xg, cache, 1, dens, vr, la=la)
         t_apply1 = min(t_apply1, tm%get_elapsed_time())
         call tm%start()
         call trc_xc_fxc_rks_apply_many(bas, xg, cache, 8, dens, vr, la=la)
         t_apply8 = min(t_apply8, tm%get_elapsed_time())
         call tm%start()
         call trc_xc_fxc_rks(bas, xg, func, dmat, dens(:, :, 1), vr(:, :, 1), la=la)
         t_once = min(t_once, tm%get_elapsed_time())
         call tm%start()
         do k = 1, 8
            call trc_xc_fxc_rks(bas, xg, func, dmat, dens(:, :, k), vr(:, :, k), la=la)
         end do
         t_each8 = min(t_each8, tm%get_elapsed_time())
      end do
      !$acc exit data delete(dmat)
      call cache%release()
      !$acc exit data delete(dens, vr)
      write (line, "(a,f8.3,a,2f8.3,a,2f8.3)") "  f_xc: cache setup ", t_setup, " s; cached many-density 1 / 8 densities ", &
         t_apply1, t_apply8, " s; uncached one call / eight calls ", t_once, t_each8
      call global_logger%info(trim(line))
   end if

   opts%conv_residual = 1.0e-8_dp
   opts%max_iteration = 100
   write (line, "(a,a,a,a,a,i0,a,i0,a,i0)") "time_cpks_batched: ", trim(xyz), " ", trim(basis_json), " nao ", nao, &
      " nocc ", nocc, " nrep ", nrep
   call global_logger%info(trim(line))
   call global_logger%info("  ncol  seq_s  batch_s  speedup  passes_seq  passes_batch  max|u_seq-u_batch|")
   do k = 1, size(NCOLS)
      ncol = NCOLS(k)
      tbest = huge(1.0_dp)
      do rep = 0, nrep
         do mode = 1, 2
            ! mode 1: one at a time (max_batch 1); mode 2: batched.
            opts%max_batch = merge(1, 8, mode == 1)
            call tm%start()
            call run_solve(ncol)
            t(mode) = tm%get_elapsed_time()
            passes(mode) = merge(sum(res%iterations), maxval(res%iterations), mode == 1)
            if (mode == 1) then
               ubatch(:, :, 1:ncol) = u(:, :, 1:ncol)
            else
               udiff = maxval(abs(ubatch(:, :, 1:ncol) - u(:, :, 1:ncol)))
            end if
            ! rep 0 is the warm-up (device code, cache setup)
            if (rep > 0) tbest(mode) = min(tbest(mode), t(mode))
         end do
      end do
      write (line, "(i6,f8.3,f9.3,f9.2,i11,i14,es14.2)") ncol, tbest(1), tbest(2), tbest(1)/tbest(2), passes(1), &
         passes(2), udiff
      call global_logger%info(trim(line))
   end do

   !$acc exit data delete(cmo, mmat, tmp, bvec, u)
   call eri%release()
   call la%release()
   if (ks) call xg%release()
   call pl%release()
   call bas%release()

contains

   subroutine run_solve(nc)
      integer, intent(in) :: nc
      u = 0.0_dp
      !$acc update device(u)
      if (ks) then
         call trc_cpks_solve_rks(bas, eri, la, xg, func, nao, nocc, nvirt, cmo, eps, dmat, nc, bvec, u, opts, res)
      else
         call trc_cphf_solve_rhf(bas, eri, la, nao, nocc, nvirt, cmo, eps, nc, bvec, u, opts, res)
      end if
      !$acc update self(u)
   end subroutine run_solve

end program time_cpks_batched
