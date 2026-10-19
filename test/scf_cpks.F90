!
! RKS static dipole polarizability through trc_cpks, against pyscf's
! pyscf.prop.polarizability.rks (test/pyscf_refs/generate.py's gen_cpks).
!
! Converges the SCF independently (functional given on the command line),
! builds the grid trc_cpks's XC kernel response will use from the SAME
! `trc_scf_run` that converged the density (`res%grid_r`/`res%grid_w`),
! and feeds that grid to pyscf unchanged -- the same "match pyscf's grid
! to terco's" pattern VV10's rks_ref.py and scf_rks.F90's probe already
! use, so the only disagreement left to measure is CPKS's own residual,
! not two different quadratures. See test/scf_cphf.F90's own header for
! why converging independently (rather than reading back a committed
! pyscf density) is safe on these closed-shell cases.
!
! Builds its OWN pair list, ERI engine and linalg handle, exactly as
! scf_cphf.F90 does, for the same reason: trc_scf_run releases everything
! it built internally once it returns.
!
! Arguments: xyz [basis.json [molname [basname [funcname [probe]]]]].
! Defaults water/6-31G*/svwn.
!
program scf_cpks
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t
   use trc_cpks, only: trc_scf_polarizability_rks, cpks_refuse_reason, trc_cpks_solve_rks
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix, ref_report
   implicit none

   integer :: natm, nocc, unit, dev, x, y
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   type(trc_eri_t) :: eri
   type(trc_linalg_t) :: la
   type(trc_xc_grid_t) :: xg
   type(trc_xc_functional_t) :: func
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: res
   type(trc_cphf_options_t) :: cphf_opts
   type(trc_cphf_result_t) :: cphf_res
   type(error_t) :: err
   logical :: ok
   character(len=256) :: xyz, bas_json, probe
   character(len=16) :: molname, basname, funcname
   real(dp) :: alpha(3, 3)
   real(dp), parameter :: ERI_THRESH = 1.0e-12_dp

   dev = trc_bind_device(0)

   xyz = 'water.xyz'; bas_json = 'basis_sets/6-31gs.json'
   molname = 'water'; basname = '6-31gs'; funcname = 'svwn'; probe = 'cpks_probe.bin'
   if (command_argument_count() >= 1) call get_command_argument(1, xyz)
   if (command_argument_count() >= 2) call get_command_argument(2, bas_json)
   if (command_argument_count() >= 3) call get_command_argument(3, molname)
   if (command_argument_count() >= 4) call get_command_argument(4, basname)
   if (command_argument_count() >= 5) call get_command_argument(5, funcname)
   if (command_argument_count() >= 6) call get_command_argument(6, probe)

   call read_xyz(trim(xyz), natm, zint, at_r)
   call trc_basis_from_json(trim(bas_json), natm, zint, at_r, bas, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', "scf_cpks: "//err%get_message()
      stop 1
   end if
   call bas%to_device()
   call pl%build(bas, ERI_THRESH)
   call pl%to_device()

   nocc = nint(sum(real(zint, dp)))/2
   opts%conv_energy = 1.0e-12_dp
   opts%conv_density = 1.0e-10_dp
   opts%conv_diis = 1.0e-10_dp
   opts%max_iter = 200
   opts%functional = trim(funcname)
   opts%grid_level = 3
   block
      character(len=8) :: envv
      integer :: st
      call get_environment_variable("TERCO_VERBOSE", envv, status=st)
      opts%verbose = st == 0
      call get_environment_variable("TERCO_CPHF_VERBOSE", envv, status=st)
      cphf_opts%verbose = st == 0
   end block
   call trc_scf_run(bas, nocc, nocc, opts, res)
   ok = res%converged
   if (.not. ok) print '(a)', "scf_cpks: "//trim(res%message)

   if (.not. xc_functional_by_name(trim(funcname), func)) then
      print '(a)', "scf_cpks: unrecognised functional '"//trim(funcname)//"'"
      stop 1
   end if

   ! The SAME grid the SCF above converged its density on -- res%grid_r/
   ! grid_w, exactly as trc_scf_context.F90 built xg (max_pts =
   ! opts%grid_max_pts, default tol/w_tol). gen_cpks reads this grid back
   ! out of the probe file below and sets it as pyscf's own.
   call xg%build(size(res%grid_w), res%grid_r, res%grid_w, bas, max_pts=opts%grid_max_pts)
   call xg%to_device()

   call eri%build(bas, ERI_THRESH)
   call la%init(bas%nao)
   cphf_opts%conv_residual = 1.0e-8_dp
   cphf_opts%max_iter = 100

   block
      type(timer_type) :: tm
      real(dp) :: t_cpks
      integer :: niter_tot
      call tm%start()
      call trc_scf_polarizability_rks(bas, pl, eri, xg, func, nocc, res%cmo(:, :, 1), res%eps(:, 1), &
                                      res%dmat(:, :, 1), alpha, cphf_opts, cphf_res, la)
      t_cpks = tm%get_elapsed_time()
      niter_tot = sum(cphf_res%iterations)
      print '(a,f10.3,a,i0,a,f10.4,a)', "  CPKS wall time ", t_cpks, " s over ", niter_tot, &
         " total CG iterations: ", t_cpks/real(max(niter_tot, 1), dp), " s/iteration"
   end block
   if (.not. cphf_res%converged) then
      ok = .false.
      print '(a)', "scf_cpks: "//trim(cphf_res%message)
   end if

   print '(a,i0,a,i0)', "scf_cpks: nao ", bas%nao, "  nocc ", nocc
   print '(a,a,a,f20.12,a,i0)', "  E(", trim(funcname), ") = ", res%energy, "   SCF iter ", res%iterations
   print '(a,3(i0,1x))', "  CPKS iterations (x,y,z RHS): ", cphf_res%iterations
   print '(a)', "  polarizability tensor (a.u.):"
   do x = 1, 3
      print '(3f16.8)', (alpha(x, y), y=1, 3)
   end do

   open (newunit=unit, file=trim(probe), access='stream', form='unformatted', status='replace')
   write (unit) int(bas%nao, int64), int(natm, int64), int(nocc, int64)
   write (unit) molname, basname, funcname
   write (unit) res%energy
   write (unit) alpha
   write (unit) int(size(res%grid_w), int64)
   write (unit) res%grid_r, res%grid_w
   close (unit)

   block
      integer :: u, nfail
      real(dp) :: e_scf_ref, alpha_ref(3, 3)
      character(len=256) :: refname
      refname = 'cpks_'//trim(molname)//'_'//trim(basname)//'_'//trim(funcname)//'.ref'
      call ref_require(trim(refname), u)
      call ref_read_real(u, 'E_SCF', e_scf_ref)
      call ref_read_matrix(u, 'ALPHA', 3, 3, alpha_ref)
      call ref_close(u)
      nfail = 0
      call ref_report(nfail, abs(res%energy - e_scf_ref) < 1.0e-9_dp, 'E vs pyscf', &
                      abs(res%energy - e_scf_ref))
      do y = 1, 3
         do x = 1, 3
            call ref_report(nfail, abs(alpha(x, y) - alpha_ref(x, y)) < 1.0e-6_dp, &
                            'alpha vs pyscf', abs(alpha(x, y) - alpha_ref(x, y)))
         end do
      end do
      if (nfail > 0) ok = .false.
   end block

   ! --- cpks_refuse_reason: pure, so no basis/SCF/device state is needed
   ! -- four trc_xc_functional_t values built DIRECTLY (not through
   ! xc_functional_by_name: no named functional in today's table has
   ! nlc_b nonzero with omega zero, since every VV10 entry is also a
   ! range-separated hybrid, which would exercise the omega branch
   ! instead -- see this module's own header) so each of the four
   ! branches is exercised on its own.
   block
      type(trc_xc_functional_t) :: f_missing, f_rsh, f_vv10, f_clean
      character(len=256) :: msg

      f_missing%missing = "needs TPSS_X: meta-GGA kernels are not translated yet"
      msg = cpks_refuse_reason(f_missing)
      call check_refused(msg, "meta-GGA (missing)", trim(f_missing%missing))

      f_rsh%omega = 0.3_dp
      f_rsh%exx_lr = 0.833_dp
      msg = cpks_refuse_reason(f_rsh)
      call check_refused(msg, "range-separated hybrid (omega)", "range-separated")

      f_vv10%nlc_b = 6.0_dp
      f_vv10%nlc_c = 0.01_dp
      msg = cpks_refuse_reason(f_vv10)
      call check_refused(msg, "VV10 (nlc_b/nlc_c)", "VV10")

      ! f_clean is default-initialised: missing/omega/exx_lr/nlc_b/nlc_c
      ! all at their zero/empty defaults, same as a plain LDA/GGA.
      msg = cpks_refuse_reason(f_clean)
      call expect(len_trim(msg) == 0, "cpks_refuse_reason: a plain functional is not refused")
   end block

   ! --- a NaN right-hand side must report NOT converged, never a silent
   ! "converged" with zero iterations -- trc_cpks_solve_rks has its own
   ! separate CG loop (not a call into trc_cphf's), with the same
   ! `.not. (resid < conv_residual)` loop condition and pap breakdown
   ! guard trc_cphf_solve_rhf's own comment explains. Reuses the
   ! already-converged SCF's own cmo/eps/dmat and the functional already
   ! selected for this run; nothing here needs to be physically
   ! meaningful.
   block
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      real(dp) :: cmo_nan(bas%nao, bas%nao), dmat_nan(bas%nao, bas%nao)
      real(dp), allocatable :: bvec_nan(:, :, :), u_nan(:, :, :)
      type(trc_cphf_result_t) :: res_nan
      integer :: nvirt_t
      real(dp) :: nan_val

      nan_val = ieee_value(nan_val, ieee_quiet_nan)
      nvirt_t = bas%nao - nocc
      cmo_nan = res%cmo(:, :, 1)
      dmat_nan = res%dmat(:, :, 1)
      allocate (bvec_nan(nvirt_t, nocc, 1), u_nan(nvirt_t, nocc, 1))
      bvec_nan = nan_val
      !$acc enter data copyin(cmo_nan, dmat_nan, bvec_nan) create(u_nan)
      call trc_cpks_solve_rks(bas, eri, la, xg, func, bas%nao, nocc, nvirt_t, cmo_nan, res%eps(:, 1), &
                              dmat_nan, 1, bvec_nan, u_nan, cphf_opts, res_nan)
      !$acc exit data delete(cmo_nan, dmat_nan, bvec_nan, u_nan)
      if (res_nan%converged) then
         ok = .false.
         print '(a)', "scf_cpks: FAILED: a NaN right-hand side must not report converged"
      else
         print '(a)', "scf_cpks: NaN right-hand side correctly refused: "//trim(res_nan%message)
      end if
      deallocate (bvec_nan, u_nan)
   end block

   print '(a)', merge("scf_cpks: PASS", "scf_cpks: FAIL", ok)
   call eri%release(); call la%release(); call xg%release(); call pl%release(); call bas%release()
   if (.not. ok) stop 1

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "scf_cpks: FAILED: "//what
         ok = .false.
      end if
   end subroutine expect

   subroutine check_refused(msg, what, substring)
      character(len=*), intent(in) :: msg, what, substring
      call expect(len_trim(msg) > 0, "cpks_refuse_reason: "//what//" must be refused")
      call expect(index(msg, substring) > 0, &
                  "cpks_refuse_reason: "//what//" message must mention '"//trim(substring)//"'")
   end subroutine check_refused

end program scf_cpks
