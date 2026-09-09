!! Segmented path against the primitive-shell view, timed, on any
!! geometry and basis:  bench_gc <xyz> <basis.json> [reps]
!!
!! One density, a decaying band as check_gc uses, so that neither path
!! gets to screen on a block-diagonal guess. Reports the build time, the
!! per-Fock time of each path, and how far the two Fock matrices differ.
program bench_gc
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_binkernel, only: trc_set_prim_margin
   implicit none
   character(len=256) :: xyzfile, basfile, arg
   character(len=512) :: dfile
   integer :: natm, n, i, j, dev, reps, r
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :), d(:, :), g1(:, :), g2(:, :), gref(:, :)
   real(dp) :: worst, t0, t1, tb1, tb2, tf1, tf2, thr, pmargin
   integer :: bres
   type(trc_basis_t) :: b
   type(trc_eri_t) :: e1, e2
   type(error_t) :: err

   reps = 3
   if (command_argument_count() < 2) then
      print '(a)', 'usage: bench_gc <xyz> <basis.json> [reps]'
      stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basfile)
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg); read (arg, *) reps
   end if
   thr = 1.0e-10_dp
   pmargin = 1.0e-3_dp
   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg); read (arg, *) thr
   end if

   if (command_argument_count() >= 5) then
      call get_command_argument(5, arg); read (arg, *) pmargin
   end if
   bres = 1
   if (command_argument_count() >= 6) then
      call get_command_argument(6, arg); read (arg, *) bres
   end if
   print '(a,i0)', '  batch res  ', bres
   call trc_set_prim_margin(pmargin)

   dev = trc_bind_device(0)
   call read_xyz(trim(xyzfile), natm, zint, at_r)
   call trc_basis_from_json(trim(basfile), natm, zint, at_r, b, err)
   if (err%has_error()) then
      print '(a)', 'bench_gc: '//err%get_message(); stop 1
   end if
   call b%to_device()
   n = b%nao
   print '(a,es9.2)', '  thresh     ', thr
   print '(a,i0,a,i0,a,i0,a,i0)', '  atoms ', natm, '  shells ', b%nshell, '  nao ', n, '  maxnp ', b%maxnp
   allocate (d(n, n), g1(n, n), g2(n, n), gref(n, n))
   ! A real density if one is offered, the model band otherwise. The band
   ! is small and short-ranged, so it screens far harder than any density an
   ! SCF actually meets and the timings come out about three times
   ! optimistic. TRC_DUMP_DENSITY on an SCF run writes one out; point
   ! TRC_BENCH_DENSITY at it here.
   call get_environment_variable("TRC_BENCH_DENSITY", dfile)
   if (len_trim(dfile) > 0) then
      call read_density(trim(dfile), n, d)
      print '(a,a)', '  density    : ', trim(dfile)
   else
      do j = 1, n
         do i = 1, n
            d(i, j) = 0.05_dp*exp(-0.3_dp*abs(i - j))
         end do
      end do
      print '(a)', '  density    : model band (optimistic; see TRC_BENCH_DENSITY)'
   end if
   !$acc enter data copyin(d) create(g1, g2, gref)

   call tick(t0)
   call e1%build(b, thr, general=.false., batch_res=bres)
   call tick(t1); tb1 = t1 - t0
   call e1%fock_resident(b, d, g1, k_scale=1.0_dp, count_survivors=.true.)   ! warm up, and count
   call tick(t0)
   do r = 1, reps
      call e1%fock_resident(b, d, g1, k_scale=1.0_dp)
   end do
   call tick(t1); tf1 = (t1 - t0)/reps

   call tick(t0)
   call e2%build(b, thr, general=.true., batch_res=bres)
   call tick(t1); tb2 = t1 - t0
   call e2%fock_resident(b, d, g2, k_scale=1.0_dp, count_survivors=.true.)
   call tick(t0)
   do r = 1, reps
      call e2%fock_resident(b, d, g2, k_scale=1.0_dp)
   end do
   call tick(t1); tf2 = (t1 - t0)/reps

   ! What the prescreen margin costs, against the same path run at a
   ! margin tight enough to be a reference. The margin is read per call,
   ! so this needs no rebuild of the object.
   call trc_set_prim_margin(1.0e-9_dp)
   call e1%fock_resident(b, d, gref, k_scale=1.0_dp)
   call trc_set_prim_margin(pmargin)
   !$acc update self(g1, g2, gref)
   print '(a,es10.2,a,es10.2)', '  margin error vs 1e-9   ', maxval(abs(g1 - gref)), &
      '   relative ', maxval(abs(g1 - gref))/max(1.0_dp, maxval(abs(gref)))
   worst = maxval(abs(g1 - g2))
   print '(a,f8.3,a,f9.3,a,i0,a,i0,a,i0)', '  segmented : build ', tb1, ' s   fock ', tf1, &
      ' s   launches ', e1%nlaunch, '   Mquartets enum ', int(e1%nwork/1000000_8), ' kept ', int(e1%nkept/1000000_8)
   print '(a,f8.3,a,f9.3,a,i0,a,i0,a,i0)', '  general   : build ', tb2, ' s   fock ', tf2, &
      ' s   launches ', e2%nlaunch, '   Mquartets enum ', int(e2%nwork/1000000_8), ' kept ', int(e2%nkept/1000000_8)
   print '(a,es10.2,a,es10.2)', '  worst |G_seg - G_gen| ', worst, '  scale ', maxval(abs(g1))
contains

   subroutine read_density(fn, nao, d)
      character(len=*), intent(in) :: fn
      integer, intent(in) :: nao
      real(dp), intent(out) :: d(nao, nao)
      integer :: u, ios
      integer(kind=8) :: nchk
      open (newunit=u, file=fn, access='stream', form='unformatted', &
            status='old', iostat=ios)
      if (ios /= 0) then
         print '(a,a)', '  cannot open ', trim(fn); stop 1
      end if
      read (u) nchk
      if (int(nchk) /= nao) then
         print '(a,i0,a,i0)', '  density nao ', nchk, ' /= ', nao; stop 1
      end if
      read (u) d
      close (u)
   end subroutine read_density
   subroutine tick(t)
      real(dp), intent(out) :: t
      integer(kind=8) :: c, rate
      call system_clock(c, rate)
      t = real(c, dp)/real(rate, dp)
   end subroutine tick
end program bench_gc
