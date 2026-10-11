program bench_sph
   !! The spherical transform's cost against one Fock build, on any geometry
   !! and a spherical basis:  bench_sph <xyz> <basis.json> [reps]
   !!
   !! `sph_sph_to_cart`/`sph_cart_to_sph` run once per Fock build already (one
   !! of each, trc_scf.F90's per-iteration D_cart/G_sph steps); this times
   !! `reps` of BOTH back to back against one `eri%fock_resident`, on a
   !! basis large enough (1000+ AO) that the O(nao^4) build and the
   !! O(nshell^2) transform are not both too small to time honestly.
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_sph, only: sph_cart_to_sph, sph_sph_to_cart
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   character(len=256) :: xyzfile, basfile, arg
   integer :: natm, dev, reps, r
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   real(dp), allocatable :: d_c(:, :), g_c(:, :), d_s(:, :), g_s(:, :)
   real(dp) :: t0, t1, t_fock, t_c2s, t_s2c
   type(trc_basis_t) :: b
   type(trc_eri_t) :: eri
   type(error_t) :: err

   reps = 10
   if (command_argument_count() < 2) then
      print '(a)', 'usage: bench_sph <xyz> <basis.json> [reps]'
      stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basfile)
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg); read (arg, *) reps
   end if

   dev = trc_bind_device(0)
   call read_xyz(trim(xyzfile), natm, zint, at_r)
   call trc_basis_from_json(trim(basfile), natm, zint, at_r, b, err)
   if (err%has_error()) then
      print '(a)', 'bench_sph: '//err%get_message(); stop 1
   end if
   call b%to_device()
   print '(a,i0,a,i0,a,i0,a,l1)', '  nao_cart ', b%nao, '   nao_sph ', b%nao_sph, &
      '   shells ', b%nshell, '   spherical ', b%spherical

   allocate (d_c(b%nao, b%nao), g_c(b%nao, b%nao))
   allocate (d_s(b%nao_sph, b%nao_sph), g_s(b%nao_sph, b%nao_sph))
   call random_number(d_c); d_c = 0.5_dp*(d_c + transpose(d_c))
   call random_number(d_s); d_s = 0.5_dp*(d_s + transpose(d_s))
   !$acc enter data copyin(d_c, g_c, d_s, g_s)

   call eri%build(b, 1.0e-10_dp)

   t0 = wall()
   call eri%fock_resident(b, d_c, g_c, k_scale=1.0_dp)
   !$acc wait
   t1 = wall()
   t_fock = t1 - t0

   t0 = wall()
   do r = 1, reps
      call sph_sph_to_cart(b, d_s, d_c)
   end do
   !$acc wait
   t1 = wall()
   t_s2c = (t1 - t0)/real(reps, dp)

   t0 = wall()
   do r = 1, reps
      call sph_cart_to_sph(b, g_c, g_s)
   end do
   !$acc wait
   t1 = wall()
   t_c2s = (t1 - t0)/real(reps, dp)

   print '(a,f12.6,a)', '  one Fock build (fock_resident)   : ', t_fock, ' s'
   print '(a,f12.6,a)', '  one sph_sph_to_cart (D_sph->D_cart): ', t_s2c, ' s'
   print '(a,f12.6,a)', '  one sph_cart_to_sph (G_cart->G_sph): ', t_c2s, ' s'
   print '(a,es10.2)', '  (D transform + G transform) / Fock build = ', (t_s2c + t_c2s)/t_fock

   !$acc exit data delete(d_c, g_c, d_s, g_s)
   call eri%release(); call b%release()

contains

   function wall() result(t)
      real(dp) :: t
      integer(kind=8) :: cc, rate
      call system_clock(cc, rate)
      t = real(cc, dp)/real(rate, dp)
   end function wall

end program bench_sph
