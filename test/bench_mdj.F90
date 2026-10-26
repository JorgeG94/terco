program bench_mdj
   !! Wall time of J by the MD engine (eri%coulomb) against the four-centre
   !! build's J-only branch (eri%fock at k_scale = 0), same density, same
   !! threshold, and their largest difference. Not a test: a measuring tool.
   !!
   !!   bench_mdj <xyz> <basis.json> [thresh] [repeats]
   !!
   !! The density decays with the distance between the AOs' centres, so the
   !! density-weighted screens of both engines have something to bite on.
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none

   character(len=256) :: xyz, basis, arg
   type(trc_basis_t) :: bas
   type(trc_eri_t) :: eri
   type(error_t) :: err
   type(timer_type) :: tm
   integer :: nat, nao, i, j, k, nrep, si, sj
   integer, allocatable :: at_z(:), aosh(:)
   real(dp), allocatable :: r(:, :), d(:, :), jref(:, :), jmd(:, :)
   real(dp) :: thresh, t_fock, t_md, t_build, r2, t_jk, t_k
   real(dp), allocatable :: g(:, :)

   call get_command_argument(1, xyz)
   call get_command_argument(2, basis)
   thresh = 1.0e-12_dp
   nrep = 3
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg)
      read (arg, *) thresh
   end if
   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg)
      read (arg, *) nrep
   end if

   call read_xyz(trim(xyz), nat, at_z, r)
   call trc_basis_from_json(trim(basis), nat, at_z, r, bas, err)
   if (err%has_error()) error stop "bench_mdj: basis"
   nao = bas%nao
   allocate (d(nao, nao), jref(nao, nao), jmd(nao, nao), aosh(nao))
   do k = 1, bas%nshell
      do i = bas%sh_ao(k), bas%sh_ao(k) + (bas%sh_l(k) + 1)*(bas%sh_l(k) + 2)/2 - 1
         aosh(i) = k
      end do
   end do
   do j = 1, nao
      sj = aosh(j)
      do i = 1, nao
         si = aosh(i)
         r2 = sum((bas%sh_r(:, si) - bas%sh_r(:, sj))**2)
         d(i, j) = (0.3_dp + 0.2_dp*cos(0.37_dp*real(i + j, dp)))*exp(-0.05_dp*r2)
      end do
   end do

   call bas%to_device()
   call tm%start()
   call eri%build(bas, thresh)
   t_build = tm%get_elapsed_time()
   print '(a,i0,a,i0,a,es8.1,a,f8.3,a)', "nao ", nao, "  shells ", bas%nshell, "  thresh ", thresh, &
      "  eri build ", t_build, " s"

   ! First call of each builds what it caches (the MD pair list); time it apart.
   call tm%start()
   call eri%coulomb(bas, d, jmd)
   print '(a,f9.4,a)', "md first call (incl. build) ", tm%get_elapsed_time(), " s"
   call tm%start()
   call eri%fock(bas, d, jref, j_scale=1.0_dp, k_scale=0.0_dp)
   print '(a,f9.4,a)', "fock J first call           ", tm%get_elapsed_time(), " s"

   t_fock = huge(1.0_dp)
   t_md = huge(1.0_dp)
   do k = 1, nrep
      call tm%start()
      call eri%fock(bas, d, jref, j_scale=1.0_dp, k_scale=0.0_dp)
      t_fock = min(t_fock, tm%get_elapsed_time())
      call tm%start()
      call eri%coulomb(bas, d, jmd)
      t_md = min(t_md, tm%get_elapsed_time())
   end do
   print '(a,f9.4,a,f9.4,a,f7.2)', "best of reps: fock J ", t_fock, " s   md J ", t_md, &
      " s   fock/md ", t_fock/t_md
   ! What the K engine has to beat: the whole HGP/Rys Fock build, and its
   ! exchange alone (j_scale = 0).
   allocate (g(nao, nao))
   t_jk = huge(1.0_dp)
   t_k = huge(1.0_dp)
   do k = 1, nrep
      call tm%start()
      call eri%fock(bas, d, g)
      t_jk = min(t_jk, tm%get_elapsed_time())
      call tm%start()
      call eri%fock(bas, d, g, j_scale=0.0_dp, k_scale=1.0_dp)
      t_k = min(t_k, tm%get_elapsed_time())
   end do
   print '(a,f9.4,a,f9.4,a)', "best of reps: fock J+K ", t_jk, " s   fock K only ", t_k, " s"
   print '(a,es10.3,a,es10.3)', "max|J_md - J_fock| ", maxval(abs(jmd - jref)), &
      "   max|J| ", maxval(abs(jref))
   print '(a,i0,a,i0)', "md primitive pairs ", eri%md%nhp, "   shell pairs ", eri%md%nsp
   call eri%release()
   call bas%release()
end program bench_mdj
