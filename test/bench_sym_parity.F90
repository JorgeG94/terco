!
! Parity skipping, timed (docs/SYMMETRY_PARITY_SKIP.md): the skeleton J+K Fock
! build of one planar molecule with the dead-axis masks and the engine frame
! against the same build without them, and the full build without a point
! group. CPU-indicative only; not registered as a ctest, a timing is not a
! pass, same as this directory's other bench_*.
!
!   bench_sym_parity naphthalene_d2h.xyz basis_sets/cc-pvdz.json
!
program bench_sym_parity
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t
   use trc_parity, only: PARITY_OFF, PARITY_AUTO
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect
   use trc_test_basis, only: read_xyz
   implicit none

   integer, parameter :: NREP = 3
   character(len=256) :: xyzfile, basfile
   type(trc_basis_t) :: b
   type(trc_sym_geom_t) :: geom
   type(trc_eri_t) :: e_on, e_off, e_c1
   type(error_t) :: err
   type(timer_type) :: tm
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :), dc(:, :), g(:, :)
   integer :: natm, r, i, j
   real(dp) :: t_on, t_off, t_c1

   xyzfile = 'naphthalene_d2h.xyz'; basfile = 'basis_sets/cc-pvdz.json'
   if (command_argument_count() >= 1) call get_command_argument(1, xyzfile)
   if (command_argument_count() >= 2) call get_command_argument(2, basfile)

   call read_xyz(trim(xyzfile), natm, zint, at_r)
   geom = sym_detect(at_r, zint, err, basis_id=zint)
   if (err%has_error()) error stop 'bench_sym_parity: sym_detect'
   call trc_basis_from_json(trim(basfile), natm, zint, geom%r_std, b, err, force_cart=.true.)
   if (err%has_error()) error stop 'bench_sym_parity: basis'
   call b%to_device()
   allocate (dc(b%nao, b%nao), g(b%nao, b%nao))
   do j = 1, b%nao
      do i = 1, b%nao
         dc(i, j) = 0.1_dp*exp(-0.1_dp*abs(i - j))
      end do
   end do

   e_on%parity_mode = PARITY_AUTO
   call e_on%build(b, 1.0e-10_dp, sym=geom)
   e_off%parity_mode = PARITY_OFF
   call e_off%build(b, 1.0e-10_dp, sym=geom)
   call e_c1%build(b, 1.0e-10_dp)
   print '(a,a,a,i0,a,a,a,i0)', '  ', trim(xyzfile), '  bf ', b%nao, '  group ', trim(geom%group%name), &
      '  dead axis ', e_on%parity_axis

   call e_on%fock(b, dc, g, skeleton=.true.)
   call e_off%fock(b, dc, g, skeleton=.true.)
   call e_c1%fock(b, dc, g)
   call tm%start()
   do r = 1, NREP
      call e_c1%fock(b, dc, g)
   end do
   t_c1 = tm%get_elapsed_time()/real(NREP, dp)
   call tm%start()
   do r = 1, NREP
      call e_off%fock(b, dc, g, skeleton=.true.)
   end do
   t_off = tm%get_elapsed_time()/real(NREP, dp)
   call tm%start()
   do r = 1, NREP
      call e_on%fock(b, dc, g, skeleton=.true.)
   end do
   t_on = tm%get_elapsed_time()/real(NREP, dp)
   print '(a,f9.3,a)', '  no symmetry, full plan:          ', t_c1, ' s'
   print '(a,f9.3,a,f6.2,a)', '  skeleton, parity skipping off:   ', t_off, ' s   (', t_c1/t_off, 'x of no symmetry)'
   print '(a,f9.3,a,f6.2,a,f6.2,a)', '  skeleton, parity skipping on:    ', t_on, ' s   (', t_c1/t_on, &
      'x of no symmetry; ', t_off/t_on, 'x of parity off)'
   print '(a,i0,a,i0,a,f7.4)', '  variant launches ', e_on%nlaunch_variant, '  masked work ', e_on%nwork_masked, &
      '  modelled ops masked/unmasked ', e_on%parity_ops_masked/max(e_on%parity_ops_base, 1.0_dp)
end program bench_sym_parity
