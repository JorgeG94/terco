!
! Old vs new Schwarz-bound setup time, CPU-indicative only (see this
! change's own report for the GPU numbers): schwarz_bounds_hgp (one HGP
! thread per shell pair, whole (ab|ab) block) against schwarz_bounds_rys
! (one Rys thread per (pair, component pair), diagonal only), on the same
! hp_* primitive-pair data. Not registered as a ctest -- a timing is not a
! pass, same as this directory's other bench_*.
!
program bench_schwarz
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_tables, only: tables_init, LMAX
   use trc_cart, only: cart_init
   use trc_boys, only: boys_init
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_hgp, only: build_pairs_hgp
   use trc_screen, only: schwarz_bounds_hgp, schwarz_bounds_rys
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none

   integer, parameter :: NREP = 3
   character(len=256) :: xyzfile, basfile
   type(trc_basis_t) :: b
   type(error_t) :: err
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   real(dp), allocatable :: one(:)
   integer, allocatable :: hp_off(:), hp_n(:), hp_ki(:), hp_kj(:)
   real(dp), allocatable :: hp_p(:), hp_r(:, :), hp_ra(:, :), hp_rb(:, :), hp_c(:)
   real(dp), allocatable :: q(:)
   integer :: natm, nhpp, r
   type(timer_type) :: tm
   real(dp) :: t_hgp, t_rys

   xyzfile = 'water.xyz'; basfile = 'basis_sets/cc-pvdz.json'
   if (command_argument_count() >= 1) call get_command_argument(1, xyzfile)
   if (command_argument_count() >= 2) call get_command_argument(2, basfile)

   call tables_init(); call cart_init(); call boys_init()
   call read_xyz(trim(xyzfile), natm, zint, at_r)
   call trc_basis_from_json(trim(basfile), natm, zint, at_r, b, err)
   if (err%has_error()) error stop 'bench_schwarz: basis'

   allocate (one(b%nshell)); one = 1.0_dp
   call build_pairs_hgp(b%nshell, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                        one, b%sh_c, 1.0e-30_dp, &
                        hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, hp_ki, hp_kj, nhpp)
   deallocate (hp_ki, hp_kj)

   allocate (q(b%nshell*(b%nshell + 1)/2))
   print '(a,i0,a,i0,a,i0,a,i0)', '  shells ', b%nshell, '   bf ', b%nao, &
      '   primitive pairs ', nhpp, '   LMAX ', LMAX

   call tm%start()
   do r = 1, NREP
      call schwarz_bounds_hgp(b%nshell, nhpp, b%sh_l, hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, q)
   end do
   t_hgp = tm%get_elapsed_time()/real(NREP, dp)

   call tm%start()
   do r = 1, NREP
      call schwarz_bounds_rys(b%nshell, nhpp, b%sh_l, hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, q)
   end do
   t_rys = tm%get_elapsed_time()/real(NREP, dp)

   print '(a,es12.4,a)', '  schwarz_bounds_hgp : ', t_hgp, ' s/call'
   print '(a,es12.4,a)', '  schwarz_bounds_rys : ', t_rys, ' s/call'
   print '(a,f8.2,a)',   '  speedup            : ', t_hgp/t_rys, 'x'
end program bench_schwarz
