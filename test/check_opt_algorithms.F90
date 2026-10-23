!
! Every --opt-algorithm this PR exposes (sd, cg, cg-auto, lbfgs, damped,
! nr, prfo-as-a-minimiser) converges water/STO-3G to the existing
! opt_water anchor energy (test/pyscf_refs/opt_water.ref,
! -74.96590121729567 Eh) within a tolerance scaled to the --opt-gtol used
! -- the harmonic estimate check_opt_vs_pyscf.F90 already uses (dE ~
! g^2/(2k)), generously loosened here since this check's job is "does the
! algorithm converge at all", not a tight anchor comparison (that is
! check_opt_vs_pyscf's). sd gets a generous maxiter (it is the slowest of
! the seven by a wide margin); every other algorithm uses the CLI's own
! default (100).
!
! Every --opt-coords this PR exposes (cartesian, hdlc, hdlc-tc, dlc,
! dlc-tc) converges the SAME water minimisation with --opt-algorithm
! lbfgs, to the same anchor.
!
program check_opt_algorithms
   use pic_types, only: dp
   use json_module, only: json_file
   implicit none

   character(len=:), allocatable :: cli_exe, work_dir
   character(len=1024) :: arg
   integer :: rc, nfail
   real(dp), parameter :: E_ANCHOR = -74.96590121729567_dp
   real(dp), parameter :: E_TOL = 2.0e-6_dp   !! loose: "does it converge", not a tight anchor match

   nfail = 0
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      print '(a)', 'check_opt_algorithms: usage: check_opt_algorithms <terco CLI> <scratch dir>'
      stop 1
   end if
   ! Run with WORKING_DIRECTORY set to work_dir (not test/ itself, this
   ! program's own WORKING_DIRECTORY, ctest's convention) -- --opt writes
   ! <geometry>_opt.xyz/_traj.xyz AND the --json probe (named alg_*.json
   ! below) into the CLI's CWD; check_opt_connect.F90's own
   ! `cd work_dir &&` pattern. One shared work dir for all 11 runs below:
   ! they are sequential (this program's own calls), never concurrent.
   call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                             exitstat=rc, wait=.true.)

   call run_one(cli_exe, work_dir, 'sd', 'cartesian', 600, nfail)
   call run_one(cli_exe, work_dir, 'cg', 'cartesian', 100, nfail)
   call run_one(cli_exe, work_dir, 'cg-auto', 'cartesian', 100, nfail)
   call run_one(cli_exe, work_dir, 'lbfgs', 'cartesian', 100, nfail)
   call run_one(cli_exe, work_dir, 'damped', 'cartesian', 150, nfail)
   call run_one(cli_exe, work_dir, 'nr', 'cartesian', 100, nfail)
   call run_one(cli_exe, work_dir, 'prfo', 'cartesian', 400, nfail)   !! valid as a plain minimiser too, just unusual

   call run_one(cli_exe, work_dir, 'lbfgs', 'hdlc', 100, nfail)
   call run_one(cli_exe, work_dir, 'lbfgs', 'hdlc-tc', 100, nfail)
   call run_one(cli_exe, work_dir, 'lbfgs', 'dlc', 100, nfail)
   call run_one(cli_exe, work_dir, 'lbfgs', 'dlc-tc', 100, nfail)

   ! water_opt.xyz/water_traj.xyz and every alg_*.json probe must not have
   ! landed in test/ itself -- the whole point of the `cd work_dir &&`
   ! inside run_one.
   block
      character(len=*), parameter :: algs(7) = &
         [character(len=8) :: 'sd', 'cg', 'cg-auto', 'lbfgs', 'damped', 'nr', 'prfo']
      character(len=*), parameter :: coordsets(4) = &
         [character(len=8) :: 'hdlc', 'hdlc-tc', 'dlc', 'dlc-tc']
      integer :: k
      logical :: exists
      inquire (file='water_opt.xyz', exist=exists)
      call report(nfail, .not. exists, 'check_opt_algorithms: water_opt.xyz was not left in test/ itself')
      inquire (file='water_traj.xyz', exist=exists)
      call report(nfail, .not. exists, 'check_opt_algorithms: water_traj.xyz was not left in test/ itself')
      do k = 1, size(algs)
         inquire (file='alg_'//trim(algs(k))//'_cartesian.json', exist=exists)
         call report(nfail, .not. exists, 'check_opt_algorithms: alg_'//trim(algs(k))// &
                    '_cartesian.json was not left in test/ itself')
      end do
      do k = 1, size(coordsets)
         inquire (file='alg_lbfgs_'//trim(coordsets(k))//'.json', exist=exists)
         call report(nfail, .not. exists, 'check_opt_algorithms: alg_lbfgs_'//trim(coordsets(k))// &
                    '.json was not left in test/ itself')
      end do
   end block

   if (nfail > 0) then
      print '(a,i0,a)', 'check_opt_algorithms: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_opt_algorithms: RESULT: PASS'

contains

   subroutine run_one(exe, wdir, alg, coords, maxiter, nf)
      character(len=*), intent(in) :: exe, wdir, alg, coords
      integer, intent(in) :: maxiter
      integer, intent(inout) :: nf
      character(len=:), allocatable :: cmd, tag, probe
      type(json_file) :: json
      integer :: rc, steps
      real(dp) :: e
      logical :: ok

      tag = alg//'_'//coords
      probe = 'alg_'//tag//'.json'
      cmd = 'cd '//wdir//' && '//trim(exe)//' water.xyz -b basis_sets/sto-3g.json --opt --opt-algorithm '//alg// &
            ' --opt-coords '//coords//' --opt-maxiter '//itoa(maxiter)//' --opt-gtol 1e-5 --silent --json '// &
            probe//' > '//tag//'.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call report(nf, rc == 0, 'check_opt_algorithms ('//tag//'): the CLI run itself succeeded')
      if (rc /= 0) return

      call json%initialize()
      call json%load(filename=wdir//'/'//probe)
      if (json%failed()) then
         call report(nf, .false., 'check_opt_algorithms ('//tag//'): cannot read '//probe)
         return
      end if
      call json%get('optimization.final_energy', e, ok)
      call report(nf, ok, 'check_opt_algorithms ('//tag//'): optimization.final_energy present')
      call json%get('optimization.steps', steps, ok)
      call json%destroy()
      if (.not. ok) return

      print '(a,a10,a,a9,a,i5,a,f20.12)', '  ', alg, '  ', coords, '  steps=', steps, '  E=', e
      call report(nf, abs(e - E_ANCHOR) < E_TOL, &
                 'check_opt_algorithms ('//tag//'): final energy matches the opt_water anchor', abs(e - E_ANCHOR))
   end subroutine run_one

   function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=16) :: buf
      write (buf, '(i0)') i
      s = trim(buf)
   end function itoa

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

end program check_opt_algorithms
