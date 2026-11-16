module check_xc_kxc_suite
   !!
   !! trc_xc_kxc_rks_apply, the third functional derivative contracted with
   !! two AO densities (src/xc/trc_xc.F90), restricted LDA/GGA.
   !!
   !! Converges an RKS SCF (functional on the command line) and builds the
   !! grid it used, the same pattern as scf_cpks.F90. The probe densities
   !! A_k, B_k (k = 1, 2) and C are fixed symmetric matrices of a closed
   !! formula (`probe_matrix`, also coded in test/pyscf_refs/generate.py's
   !! gen_kxc), so the pyscf reference needs only the converged density and
   !! grid, which this program dumps as a probe file for generate.py.
   !!
   !! Checks: the traces Tr(C v[A_k,B_k]) against pyscf (libxc
   !! eval_xc(deriv=3) contracted by hand at terco's own density and grid);
   !! the symmetry E'''[A,B,C] under permutations; the density-direction
   !! finite difference of the existing trc_xc_fxc_rks; a batch of pairs
   !! against a loop of single pairs; and the refusals.
   !!
   !! Arguments: xyz [basis.json [molname [basname [funcname [probe]]]]].
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name, XC_FAMILY_MGGA
   use trc_xc, only: trc_xc_fxc_rks, trc_xc_kxc_rks_apply, trc_xc_kxc_refuse_reason
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_vector
   implicit none
   private
   public :: collect_check_xc_kxc, xyz_arg, bas_json_arg, molname_arg, basname_arg, funcname_arg, probe_arg

   integer, parameter :: NVEC = 2
   real(dp), parameter :: TOL_PYSCF = 1.0e-8_dp
   !! Relative to the largest |Tr| of the case. Measured at the committed
   !! step (XC_KXC_STEP): 2e-9 to 4e-9, the floor of the comparison.
   real(dp), parameter :: TOL_SYM = 1.0e-9_dp
   !! The permutation symmetry of the third derivative, relative to the
   !! largest |Tr|. Measured: 1.5e-11 to 1.6e-10 (A,B swap), 1e-12 to 3e-12
   !! (A,C swap).
   real(dp), parameter :: TOL_FD = 1.0e-8_dp
   !! Against the Richardson finite difference of trc_xc_fxc_rks (steps
   !! 1e-3 and 5e-4 of B). Measured 3e-10 to 1.7e-9.
   real(dp), parameter :: TOL_BATCH = 1.0e-13_dp

   character(len=256) :: xyz_arg = 'water.xyz'
   character(len=256) :: bas_json_arg = 'basis_sets/6-31gs.json'
   character(len=16) :: molname_arg = 'water'
   character(len=16) :: basname_arg = '6-31gs'
   character(len=16) :: funcname_arg = 'svwn'
   character(len=256) :: probe_arg = 'kxc_probe.bin'

   type :: kxc_case_t
      !! Everything the checks share, built once by `setup_case`.
      type(trc_basis_t) :: bas
      type(trc_xc_grid_t) :: xg
      type(trc_xc_functional_t) :: func
      type(trc_linalg_t) :: la
      real(dp), allocatable :: dmat(:, :)
      real(dp), allocatable :: a(:, :, :), bm(:, :, :), c(:, :)
   end type kxc_case_t

contains

   subroutine collect_check_xc_kxc(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("kxc traces vs pyscf, symmetry, fxc finite difference, batch == loop", test_kxc), &
                  new_unittest("kxc refuses meta-GGA, VV10, range separation, unavailable", test_refuse) &
                  ]
   end subroutine collect_check_xc_kxc

   pure function probe_matrix(n, p) result(m)
      !! The symmetric probe density number p: 0.05 cos(0.37 p (i+j) + 0.013 i j),
      !! 1-based i, j (test/pyscf_refs/generate.py's gen_kxc has the same).
      integer, intent(in) :: n, p
      real(dp) :: m(n, n)
      integer :: i, j
      do j = 1, n
         do i = 1, n
            m(i, j) = 0.05_dp*cos(0.37_dp*real(p, dp)*real(i + j, dp) + 0.013_dp*real(i*j, dp))
         end do
      end do
   end function probe_matrix

   function trace_prod(x, y) result(t)
      real(dp), intent(in) :: x(:, :), y(:, :)
      real(dp) :: t
      t = sum(x*y)
   end function trace_prod

   subroutine setup_case(kc, error)
      type(kxc_case_t), intent(inout) :: kc
      type(error_type), allocatable, intent(out) :: error
      integer :: natm, nocc, unit, dev, k
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_pairlist_t) :: pl
      type(error_t) :: err
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      logical :: func_ok

      dev = trc_bind_device(0)
      call read_xyz(trim(xyz_arg), natm, zint, at_r)
      call trc_basis_from_json(trim(bas_json_arg), natm, zint, at_r, kc%bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call kc%bas%to_device()

      nocc = nint(sum(real(zint, dp)))/2
      opts%conv_energy = 1.0e-12_dp
      opts%conv_density = 1.0e-10_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      opts%functional = trim(funcname_arg)
      opts%grid_level = 3
      call trc_scf_run(kc%bas, nocc, nocc, opts, res)
      call check(error, res%converged, "SCF did not converge")
      if (allocated(error)) return
      func_ok = xc_functional_by_name(trim(funcname_arg), kc%func)
      call check(error, func_ok, "unrecognised functional '"//trim(funcname_arg)//"'")
      if (allocated(error)) return

      call kc%xg%build(size(res%grid_w), res%grid_r, res%grid_w, kc%bas, max_pts=opts%grid_max_pts)
      call kc%xg%to_device()
      call kc%la%init(kc%bas%nao)

      kc%dmat = res%dmat(:, :, 1)
      allocate (kc%a(kc%bas%nao, kc%bas%nao, NVEC), kc%bm(kc%bas%nao, kc%bas%nao, NVEC))
      do k = 1, NVEC
         kc%a(:, :, k) = probe_matrix(kc%bas%nao, 2*k - 1)
         kc%bm(:, :, k) = probe_matrix(kc%bas%nao, 2*k)
      end do
      kc%c = probe_matrix(kc%bas%nao, 5)

      ! The overlap goes in the probe: generate.py needs terco's AO
      ! normalisation to carry the densities into pyscf's basis.
      allocate (smat(kc%bas%nao, kc%bas%nao), tmat(kc%bas%nao, kc%bas%nao), vmat(kc%bas%nao, kc%bas%nao))
      call pl%build(kc%bas, 1.0e-12_dp)
      call pl%to_device()
      !$acc enter data create(smat, tmat, vmat)
      call trc_1e(kc%bas, pl, smat, tmat, vmat)
      !$acc update self(smat)
      !$acc exit data delete(smat, tmat, vmat)
      call pl%release()

      open (newunit=unit, file=trim(probe_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(kc%bas%nao, int64), int(natm, int64), int(nocc, int64)
      write (unit) molname_arg, basname_arg, funcname_arg
      write (unit) res%energy
      write (unit) kc%dmat
      write (unit) smat
      write (unit) int(size(res%grid_w), int64)
      write (unit) res%grid_r, res%grid_w
      close (unit)
   end subroutine setup_case

   subroutine release_case(kc)
      type(kxc_case_t), intent(inout) :: kc
      call kc%la%release()
      call kc%xg%release()
      call kc%bas%release()
   end subroutine release_case

   subroutine kxc_apply(kc, nv, da, db, v, error)
      !! trc_xc_kxc_rks_apply with the failure turned into a test error.
      type(kxc_case_t), intent(inout) :: kc
      integer, intent(in) :: nv
      real(dp), intent(in) :: da(:, :, :), db(:, :, :)
      real(dp), intent(out) :: v(:, :, :)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      call trc_xc_kxc_rks_apply(kc%bas, kc%xg, kc%func, kc%dmat, da, db, nv, v, err, la=kc%la)
      call check(error, .not. err%has_error(), err%get_message())
   end subroutine kxc_apply

   subroutine test_kxc(error)
      type(error_type), allocatable, intent(out) :: error
      type(kxc_case_t) :: kc
      real(dp), allocatable :: v(:, :, :), v1(:, :, :), vs(:, :, :), tr(:), ref(:)
      real(dp) :: scale, dev, t_ab, t_ba, t_cb_a, h, fd1, fd2, fd_rich, t_alt
      integer :: u, n, k

      call setup_case(kc, error)
      if (allocated(error)) then
         call release_case(kc)
         return
      end if
      n = kc%bas%nao
      allocate (v(n, n, NVEC), v1(n, n, 1), vs(n, n, 1), tr(NVEC), ref(NVEC))

      ! --- 1. against pyscf ------------------------------------------------
      call kxc_apply(kc, NVEC, kc%a, kc%bm, v, error)
      if (allocated(error)) then
         call release_case(kc)
         return
      end if
      do k = 1, NVEC
         tr(k) = trace_prod(kc%c, v(:, :, k))
      end do
      block
         character(len=256) :: refname
         refname = 'kxc_'//trim(molname_arg)//'_'//trim(basname_arg)//'_'//trim(funcname_arg)//'.ref'
         call ref_require(trim(refname), u)
         call ref_read_vector(u, 'TR_CAB', NVEC, ref)
         call ref_close(u)
      end block
      scale = maxval(abs(ref))
      dev = maxval(abs(tr - ref))/scale
      call check(error, scale > 1.0e-6_dp, "reference traces are not vanishing (probe is not testing anything)")
      if (allocated(error)) then
         call release_case(kc)
         return
      end if
      call check(error, dev < TOL_PYSCF, "Tr(C v[A,B]) vs pyscf, worst relative deviation "// &
                 to_string(dev, "es10.2")//": terco "//to_string(tr(1), "es18.10")//" pyscf "//to_string(ref(1), "es18.10"))
      if (allocated(error)) then
         call release_case(kc)
         return
      end if

      ! --- 2. symmetry E'''[A,B,C] = E'''[B,A,C] = E'''[C,B,A] --------------
      call kxc_apply(kc, 1, kc%bm(:, :, 1:1), kc%a(:, :, 1:1), v1, error)
      if (allocated(error)) then
         call release_case(kc)
         return
      end if
      t_ab = tr(1)
      t_ba = trace_prod(kc%c, v1(:, :, 1))
      call kxc_apply(kc, 1, reshape(kc%c, [n, n, 1]), kc%bm(:, :, 1:1), vs, error)
      if (allocated(error)) then
         call release_case(kc)
         return
      end if
      t_cb_a = trace_prod(kc%a(:, :, 1), vs(:, :, 1))
      call check(error, abs(t_ab - t_ba) < TOL_SYM*scale, "Tr(C v[A,B]) == Tr(C v[B,A]): "// &
                 to_string(t_ab, "es18.10")//" vs "//to_string(t_ba, "es18.10"))
      if (allocated(error)) then
         call release_case(kc)
         return
      end if
      call check(error, abs(t_ab - t_cb_a) < TOL_SYM*scale, "Tr(C v[A,B]) == Tr(A v[C,B]): "// &
                 to_string(t_ab, "es18.10")//" vs "//to_string(t_cb_a, "es18.10"))
      if (allocated(error)) then
         call release_case(kc)
         return
      end if
      ! Negative control: the same comparison with a different B must NOT
      ! agree, or the tolerance above would not resolve a wrong operator.
      call kxc_apply(kc, 1, kc%a(:, :, 1:1), kc%a(:, :, 2:2), vs, error)
      if (allocated(error)) then
         call release_case(kc)
         return
      end if
      t_alt = trace_prod(kc%c, vs(:, :, 1))
      call check(error, abs(t_ab - t_alt) > 100.0_dp*TOL_SYM*scale, &
                 "control: a different B must change the trace by far more than the symmetry tolerance")
      if (allocated(error)) then
         call release_case(kc)
         return
      end if

      ! --- 3. finite difference of trc_xc_fxc_rks along B --------------------
      ! d/dt Tr(C f_xc[D + t B](A)) at 0, central, Richardson in the step.
      h = 1.0e-3_dp
      fd1 = fxc_trace_fd(kc, h)
      fd2 = fxc_trace_fd(kc, 0.5_dp*h)
      fd_rich = (4.0_dp*fd2 - fd1)/3.0_dp
      call check(error, abs(fd_rich - t_ab) < TOL_FD*scale, "kxc vs finite difference of trc_xc_fxc_rks: "// &
                 to_string(t_ab, "es18.10")//" vs "//to_string(fd_rich, "es18.10"))
      if (allocated(error)) then
         call release_case(kc)
         return
      end if

      ! --- 4. a batch of pairs equals a loop of single pairs ------------------
      do k = 1, NVEC
         call kxc_apply(kc, 1, kc%a(:, :, k:k), kc%bm(:, :, k:k), v1, error)
         if (allocated(error)) then
            call release_case(kc)
            return
         end if
         dev = maxval(abs(v(:, :, k) - v1(:, :, 1)))/max(maxval(abs(v(:, :, k))), tiny(1.0_dp))
         call check(error, dev < TOL_BATCH, "batch vs single pair "//to_string(k)//": relative "// &
                    to_string(dev, "es10.2"))
         if (allocated(error)) then
            call release_case(kc)
            return
         end if
      end do
      call release_case(kc)
   end subroutine test_kxc

   function fxc_trace_fd(kc, h) result(fd)
      !! (Tr(C f_xc[D + h B1](A1)) - Tr(C f_xc[D - h B1](A1)))/(2h).
      type(kxc_case_t), intent(inout) :: kc
      real(dp), intent(in) :: h
      real(dp) :: fd
      real(dp), allocatable :: dpl(:, :), dmi(:, :), rp(:, :), rm(:, :)
      integer :: n
      n = kc%bas%nao
      allocate (dpl(n, n), dmi(n, n), rp(n, n), rm(n, n))
      dpl = kc%dmat + h*kc%bm(:, :, 1)
      dmi = kc%dmat - h*kc%bm(:, :, 1)
      call trc_xc_fxc_rks(kc%bas, kc%xg, kc%func, dpl, kc%a(:, :, 1), rp, la=kc%la)
      call trc_xc_fxc_rks(kc%bas, kc%xg, kc%func, dmi, kc%a(:, :, 1), rm, la=kc%la)
      fd = (trace_prod(kc%c, rp) - trace_prod(kc%c, rm))/(2.0_dp*h)
   end function fxc_trace_fd

   subroutine test_refuse(error)
      !! No SCF needed: the refusal precedes any use of the basis or grid,
      !! so empty ones with zero-size arrays suffice.
      type(error_type), allocatable, intent(out) :: error
      type(trc_xc_functional_t) :: f_missing, f_mgga, f_rsh, f_vv10, f_clean
      type(trc_basis_t) :: bas
      type(trc_xc_grid_t) :: xg
      type(error_t) :: err
      real(dp) :: d(0, 0), da(0, 0, 1), db(0, 0, 1), v(0, 0, 1)

      f_missing%missing = "needs the density Laplacian"
      f_mgga%family = XC_FAMILY_MGGA
      f_rsh%omega = 0.3_dp
      f_rsh%exx_lr = 0.8_dp
      f_vv10%nlc_b = 6.0_dp
      f_vv10%nlc_c = 0.01_dp
      call check(error, index(trc_xc_kxc_refuse_reason(f_missing), "not available") > 0, "missing is refused")
      if (allocated(error)) return
      call check(error, index(trc_xc_kxc_refuse_reason(f_mgga), "meta-GGA") > 0, "meta-GGA is refused")
      if (allocated(error)) return
      call check(error, index(trc_xc_kxc_refuse_reason(f_rsh), "range-separated") > 0, "range separation is refused")
      if (allocated(error)) return
      call check(error, index(trc_xc_kxc_refuse_reason(f_vv10), "VV10") > 0, "VV10 is refused")
      if (allocated(error)) return
      call check(error, len_trim(trc_xc_kxc_refuse_reason(f_clean)) == 0, "a plain functional is not refused")
      if (allocated(error)) return
      v = 1.0_dp
      call trc_xc_kxc_rks_apply(bas, xg, f_mgga, d, da, db, 1, v, err)
      call check(error, err%has_error(), "trc_xc_kxc_rks_apply must report the refusal")
   end subroutine test_refuse

end module check_xc_kxc_suite

program check_xc_kxc
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_xc_kxc_suite, only: collect_check_xc_kxc, xyz_arg, bas_json_arg, molname_arg, basname_arg, &
                                  funcname_arg, probe_arg
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 2) call get_command_argument(2, bas_json_arg)
   if (command_argument_count() >= 3) call get_command_argument(3, molname_arg)
   if (command_argument_count() >= 4) call get_command_argument(4, basname_arg)
   if (command_argument_count() >= 5) call get_command_argument(5, funcname_arg)
   if (command_argument_count() >= 6) call get_command_argument(6, probe_arg)

   stat = 0
   call run_testsuite(collect_check_xc_kxc, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_xc_kxc
