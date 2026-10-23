!
! --opt-connect (DL-FIND task 1011): after a saddle search, displace along
! the imaginary mode both ways and minimise each to a downhill minimum.
! Works with --opt-saddle-method dimer; --opt-saddle-method prfo is
! refused at the CLI (app/main.F90) because it crashes INSIDE libdlfind
! itself -- confirmed with a debug (-fcheck=all) build of libdlfind: the
! fault is a bad array copy at dlf_formstep.f90:2362 (`array(:)=
! tsmode_r(:)` in dlf_formstep_get_ra, called from dlf_task.f90's
! task=1011 case right after the P-RFO TS sub-search converges). The size
! check one line above the fault passes, so `tsmode_r` itself is a
! corrupted allocation inside libdlfind's own TS-mode bookkeeping
! (dlf_formstep_set_tsmode's mode 10/11 path, which only P-RFO/NR ever
! exercise) by the time the downhill phase asks for it -- not a size
! mismatch or missing argument this bridge supplies or could catch. The
! dimer method's own TS-mode tracking never goes through that code path
! and does not reproduce the fault, which is what this check exercises
! end to end: --opt-target saddle --opt-saddle-method dimer --opt-connect
! from test/hcn_ts_guess.xyz finds the TS and both downhill minima.
!
! Runs the CLI with WORKING_DIRECTORY set to a scratch dir (not test/
! itself): DL-FIND writes path.xyz/path_active.xyz/path_tsmode.xyz/
! path_force.xyz unconditionally into its CWD once printf>=3 (dl-find.f90
! -- confirmed by reading it; the condition is `printf`, not `printl`),
! and a `connect` run needs printf=4 for the TS/minima callbacks to fire
! at all. This check also asserts those files do NOT land in test/
! itself.
!
! Validation: the two minima are HCN and HNC. Whichever is energetically
! lower (HCN; opt_hcn.ref, -91.67520899058249 Eh) is checked against that
! EXISTING anchor directly, since a connect run from this start should
! land on the SAME HCN minimum --opt alone already validates there. The
! higher one (HNC) is checked with --hessian for
! ZERO imaginary modes -- a genuine minimum, not a second saddle a bad
! connect run slid into.
!
program check_opt_connect
   use pic_types, only: dp
   use json_module, only: json_file
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   implicit none

   character(len=:), allocatable :: cli_exe, work_dir, cmd
   character(len=1024) :: arg
   integer :: rc, nfail
   logical :: converged, ts_found, min_minus_found, min_plus_found, found
   real(dp) :: e_ts, e_minus, e_plus
   type(json_file) :: json
   real(dp), parameter :: ANCHOR_ENERGY_TOL = 1.0e-7_dp   !! Eh -- same reasoning as check_opt_ts.F90's own

   nfail = 0
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      print '(a)', 'check_opt_connect: usage: check_opt_connect <terco CLI> <scratch dir>'
      stop 1
   end if
   ! Everything the subprocess needs is copied INTO work_dir, so `cd
   ! work_dir &&` below can use plain relative paths: DL-FIND's own
   ! path*.xyz (gated on `printf`, confirmed by reading dl-find.f90)
   ! then land there too, not in test/ itself (this program's own
   ! WORKING_DIRECTORY, ctest's convention) -- see this file's own
   ! header.
   call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f hcn_ts_guess.xyz '//work_dir//'/hcn_ts_guess.xyz', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                             exitstat=rc, wait=.true.)

   cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn_ts_guess.xyz -b basis_sets/sto-3g.json --opt '// &
         '--opt-target saddle --opt-saddle-method dimer '// &
         '--opt-connect --opt-maxiter 150 --opt-gtol 1e-5 --silent --json connect.json '// &
         '> connect.out 2>&1'
   call execute_command_line(cmd, exitstat=rc, wait=.true.)
   call expect(rc == 0, 'check_opt_connect: the CLI run itself failed')

   call json%initialize()
   call json%load(filename=work_dir//'/connect.json')
   call expect(.not. json%failed(), 'check_opt_connect: cannot read connect.json')
   call json%get('optimization.converged', converged, found)
   call expect(found, 'check_opt_connect: no optimization.converged')
   call report(nfail, converged, 'check_opt_connect: the TS sub-search converged')
   call json%get('optimization.ts_found', ts_found, found); call expect(found, 'check_opt_connect: no ts_found')
   call report(nfail, ts_found, 'check_opt_connect: a transition state was found')
   call json%get('optimization.ts_energy', e_ts, found)
   call report(nfail, found, 'check_opt_connect: ts_energy present')
   call json%get('optimization.minimum_minus_found', min_minus_found, found)
   call expect(found, 'check_opt_connect: no minimum_minus_found')
   call report(nfail, min_minus_found, 'check_opt_connect: the -distort downhill minimum was found')
   call json%get('optimization.minimum_minus_energy', e_minus, found)
   call report(nfail, found, 'check_opt_connect: minimum_minus_energy present')
   call json%get('optimization.minimum_plus_found', min_plus_found, found)
   call expect(found, 'check_opt_connect: no minimum_plus_found')
   call report(nfail, min_plus_found, 'check_opt_connect: the +distort downhill minimum was found')
   call json%get('optimization.minimum_plus_energy', e_plus, found)
   call report(nfail, found, 'check_opt_connect: minimum_plus_energy present')
   call json%destroy()

   print '(a,f20.12)', '  transition state  E = ', e_ts
   print '(a,f20.12)', '  minimum (-distort) E = ', e_minus
   print '(a,f20.12)', '  minimum (+distort) E = ', e_plus

   call report(nfail, e_ts > e_minus .and. e_ts > e_plus, &
              'check_opt_connect: the transition state sits above BOTH downhill minima')

   if (min_minus_found .and. min_plus_found) then
      if (e_minus < e_plus) then
         call check_hcn_side(e_minus, nfail)
         call check_hnc_side(cli_exe, work_dir, 'minus', nfail)
      else
         call check_hcn_side(e_plus, nfail)
         call check_hnc_side(cli_exe, work_dir, 'plus', nfail)
      end if
   end if

   ! DL-FIND's own path*.xyz/fort.32, and --opt's own
   ! <geometry>_opt.xyz/_traj.xyz/_ts.xyz/_minus.xyz/_plus.xyz
   ! (app/main.F90), must not have landed in test/ itself (this program's
   ! own WORKING_DIRECTORY, ctest's convention) -- the whole point of the
   ! `cd work_dir &&` above.
   block
      character(len=*), parameter :: stray(9) = &
         [character(len=22) :: 'path.xyz', 'path_active.xyz', 'path_tsmode.xyz', 'path_force.xyz', &
                                'hcn_ts_guess_opt.xyz', 'hcn_ts_guess_traj.xyz', 'hcn_ts_guess_ts.xyz', &
                                'hcn_ts_guess_minus.xyz', 'hcn_ts_guess_plus.xyz']
      integer :: k
      logical :: exists
      do k = 1, size(stray)
         inquire (file=trim(stray(k)), exist=exists)
         call report(nfail, .not. exists, 'check_opt_connect: '//trim(stray(k))// &
                    ' was not left in test/ itself')
      end do
   end block

   if (nfail > 0) then
      print '(a,i0,a)', 'check_opt_connect: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_opt_connect: RESULT: PASS'

contains

   !
   ! The lower-energy minimum must be HCN: opt_hcn.ref's own anchor
   ! (test/pyscf_refs/README.md's two-tier scheme -- the SAME anchor
   ! check_opt_vs_pyscf.F90 reads for a plain --opt run).
   !
   subroutine check_hcn_side(e_hcn, nf)
      real(dp), intent(in) :: e_hcn
      integer, intent(inout) :: nf
      integer :: uu, nat_ref
      real(dp), allocatable :: z_ref(:), coords_ref(:, :)
      real(dp) :: e_ref, gmax_ref
      call ref_require('opt_hcn.ref', uu)
      call ref_read_int(uu, 'NAT', nat_ref)
      allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
      call ref_read_vector(uu, 'Z', nat_ref, z_ref)
      call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
      call ref_read_real(uu, 'E_FINAL', e_ref)
      call ref_read_real(uu, 'GMAX_FINAL', gmax_ref)
      call ref_close(uu)
      call ref_report(nf, abs(e_hcn - e_ref) < ANCHOR_ENERGY_TOL, &
                      'check_opt_connect: the lower downhill minimum matches the opt_hcn anchor, Eh', &
                      abs(e_hcn - e_ref))
   end subroutine check_hcn_side

   !
   ! The higher-energy minimum (HNC) must be a genuine minimum: --hessian
   ! on its own geometry shows zero imaginary modes.
   !
   subroutine check_hnc_side(exe, wdir, side, nf)
      character(len=*), intent(in) :: exe, wdir, side
      integer, intent(inout) :: nf
      character(len=:), allocatable :: cmd2
      type(json_file) :: hj
      integer :: rc2, n_imag
      logical :: ok
      cmd2 = 'cd '//wdir//' && '//trim(exe)//' hcn_ts_guess_'//side//'.xyz -b basis_sets/sto-3g.json '// &
             '--hessian --silent --json hnc_hess.json > hnc_hess.out 2>&1'
      call execute_command_line(cmd2, exitstat=rc2, wait=.true.)
      call expect(rc2 == 0, 'check_opt_connect: the HNC --hessian re-check CLI run itself failed')
      call hj%initialize()
      call hj%load(filename=wdir//'/hnc_hess.json')
      call expect(.not. hj%failed(), 'check_opt_connect: cannot read hnc_hess.json')
      call hj%get('frequencies.n_imaginary', n_imag, ok)
      call expect(ok, 'check_opt_connect: no frequencies.n_imaginary')
      call hj%destroy()
      call report(nf, n_imag == 0, 'check_opt_connect: the higher downhill minimum (HNC) has zero '// &
                 'imaginary modes', real(n_imag, dp))
   end subroutine check_hnc_side

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', what; stop 1
      end if
   end subroutine expect

   subroutine report(nf, cond, what, val)
      integer, intent(inout) :: nf
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      real(dp), intent(in), optional :: val
      if (cond) then
         print '(a,a)', '  PASS  ', what
      else
         if (present(val)) then
            print '(a,a,a,es12.4)', '  FAIL  ', what, '   value=', val
         else
            print '(a,a)', '  FAIL  ', what
         end if
         nf = nf + 1
      end if
   end subroutine report

end program check_opt_connect
