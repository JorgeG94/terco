program time_grad_pairs
   !! Wall time of eri%gradient_pairs against eri%gradient, per call, for one
   !! molecule and basis: the run-time measurement the pair-kernel form
   !! (PAIRS_RT_PLAN in scripts/gen_rys_grad.py) needs and no test makes.
   !! Not a ctest: it asserts nothing and is built as `time_grad_pairs`.
   !!
   !!   time_grad_pairs <basis.json> <molecule.xyz> [ncol=4] [nrep=3]
   !!
   !! Run it from test/ (the basis path is relative), e.g. on a GPU node:
   !!   ../build/test/standalone/time_grad_pairs basis_sets/cc-pvtz.json ../bench/md-engine/w5.xyz 1
   !!   ../build/test/standalone/time_grad_pairs basis_sets/cc-pvtz.json ../bench/md-engine/w5.xyz 4
   !! One RHF-like column is the pair table (1, 1, 1.0, 0.5), the same as
   !! eri%gradient's own density; further columns reuse it with a different
   !! density, so every column's work is of the same size. The first call of
   !! each kind is a warm-up and not counted; the table lists the others.
   !! Compare a build with a class in PAIRS_RT_PLAN against one without
   !! (delete the entry, reconfigure) at the same ncol; both builds print the
   !! eri%gradient time too, which the pair table does not change.
   !!
   !! Careful with that comparison on nvfortran: removing a class's
   !! PAIRS_RT_PLAN entry rebuilds its LITERAL pair kernels, each of which
   !! takes over 600 s to compile for (21|21), (31|30) and (33|10) (see the
   !! table in scripts/gen_rys_grad.py). `--rt-plan none` forces literal pair
   !! kernels for one class in a hand-run generator, but CMake cannot pass it.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use pic_timer, only: timer_type
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t, trc_grad_pair_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none

   type(trc_basis_t) :: bas
   type(trc_eri_t) :: eri
   type(error_t) :: err
   type(timer_type) :: tm
   character(len=256) :: basis_json, xyz, arg
   integer :: nat, ncol, nrep, c, k, ios
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: r(:, :), pool(:, :, :), g0(:, :), g1(:, :, :)
   type(trc_grad_pair_t), allocatable :: ps(:)
   real(dp) :: tg, tp, t0

   if (command_argument_count() < 2) then
      call global_logger%error("usage: time_grad_pairs <basis.json> <molecule.xyz> [ncol=4] [nrep=3]")
      error stop 1
   end if
   call get_command_argument(1, basis_json)
   call get_command_argument(2, xyz)
   ncol = 4
   nrep = 3
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg)
      read (arg, *, iostat=ios) ncol
      if (ios /= 0 .or. ncol < 1) then
         call global_logger%error("time_grad_pairs: ncol must be a positive integer")
         error stop 1
      end if
   end if
   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg)
      read (arg, *, iostat=ios) nrep
      if (ios /= 0 .or. nrep < 1) then
         call global_logger%error("time_grad_pairs: nrep must be a positive integer")
         error stop 1
      end if
   end if

   call read_xyz(trim(xyz), nat, at_z, r)
   call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err, &
                            overrides=trc_basis_json_overrides(prune=.false.))
   if (err%has_error()) then
      call global_logger%error("basis: "//err%get_message())
      error stop 1
   end if
   call bas%to_device()
   call eri%build(bas, 1.0e-12_dp)

   allocate (pool(bas%nao, bas%nao, ncol), g0(3, nat), g1(3, nat, ncol), ps(ncol))
   do c = 1, ncol
      call synth_density(bas%nao, c, pool(:, :, c))
      ps(c) = trc_grad_pair_t(left=c, right=c, column=c, coulomb=1.0_dp, exchange=0.5_dp)
   end do

   ! Warm-up, then nrep timed calls of each.
   call eri%gradient(bas, pool(:, :, 1:1), g0)
   call eri%gradient_pairs(bas, pool, ps, ncol, g1, err)
   if (err%has_error()) then
      call global_logger%error("gradient_pairs: "//err%get_message())
      error stop 1
   end if
   tg = 0.0_dp
   tp = 0.0_dp
   do k = 1, nrep
      call tm%start()
      call eri%gradient(bas, pool(:, :, 1:1), g0)
      call tm%stop()
      t0 = tm%get_elapsed_time()
      tg = tg + t0
      call tm%start()
      call eri%gradient_pairs(bas, pool, ps, ncol, g1, err)
      call tm%stop()
      t0 = tm%get_elapsed_time()
      tp = tp + t0
      call global_logger%info("  rep "//to_string(k)//": gradient "//to_string(tg/k, "f10.4")// &
                              " s (mean so far), gradient_pairs "//to_string(tp/k, "f10.4")//" s")
   end do
   call global_logger%info(trim(basis_json)//" "//trim(xyz)//": nao "//to_string(bas%nao)// &
                           ", ncol "//to_string(ncol)//", mean of "//to_string(nrep)// &
                           ": eri%gradient "//to_string(tg/nrep, "f10.4")//" s, eri%gradient_pairs "// &
                           to_string(tp/nrep, "f10.4")//" s")
   call eri%release()
   call bas%release()

contains

   subroutine synth_density(nao, phase, d)
      integer, intent(in) :: nao, phase
      real(dp), intent(out) :: d(nao, nao)
      integer :: i, j
      do j = 1, nao
         do i = 1, nao
            d(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j + 7*phase, dp))
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
   end subroutine synth_density

end program time_grad_pairs
