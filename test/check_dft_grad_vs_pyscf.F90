!
! The FULL KS-DFT nuclear gradient (trc_scf_grad%trc_scf_gradient with
! `func` present) against pyscf's own grad.RKS/grad.UKS, at TERCO'S OWN
! converged orbitals and TERCO'S OWN quadrature -- the same "compare at
! terco's own converged density and grid" recipe check_xc_grad_vs_pyscf.F90
! uses, but WITHOUT that program's force_cart=.true.: the basis here loads
! in whatever space its own JSON tags it ("gto_spherical" or not), so
! cc-pVDZ comes back SPHERICAL and 6-31G* comes back CARTESIAN, and
! trc_scf_gradient's own spherical fold (D and W, nao_sph -> Cartesian nao,
! T D T^T -- see that module's header) runs for real on the spherical
! cases. This is the combination neither sibling branch tested alone:
! grad-xc's own check_xc_grad_vs_pyscf.F90 forced Cartesian throughout
! (its header says so explicitly, "the sibling branch's spherical density
! folding is what adds it"), and grad-cli's check_full_grad_vs_pyscf.F90
! exercises the fold only for the plain (functional-free) gradient.
!
! Also passes grid_level/grid_max_pts/rho_tol through from the SAME `opts`
! the SCF converged with -- trc_capi.F90's trc_run_gradient does the same,
! and this program is the one that checks the gradient's quadrature
! actually is the SCF's, not merely "a reasonable one".
!
! grid_response=.false. EXPLICITLY: this program and its committed
! dft_grad_*.ref were both written (and pyscf's own g.grid_response left
! at its False default -- test/pyscf_refs/generate.py's gen_dft_grad)
! before trc_xc_grad.F90's grid-response term existed at all, so every
! comparison here is Pulay-only, same as check_xc_grad_vs_pyscf.F90's own
! default (non-gridresp) cases. trc_scf_gradient's OWN default for this
! argument is now ON (the module it forwards to, trc_xc_grad.F90's own
! header) -- omitting it here would silently add the grid-response term
! into terco's side only and fail against this Pulay-only reference.
!
! CASES: water/formamide, RKS, SVWN/PBE/B3LYP (same closed-shell molecules
! as check_xc_grad_vs_pyscf.F90, now on cc-pVDZ *spherical* or 6-31G*
! *Cartesian* rather than cc-pVDZ forced Cartesian); and NH2, UKS PBE,
! spherical cc-pVDZ, FROM A FIXED INITIAL DENSITY -- `opts%guess = "core"`
! (trc_scf_settings.F90's TRC_GUESS_CORE equivalent on the Fortran side: a
! bare one-electron-Hamiltonian diagonalisation, fully determined by the
! molecule and basis alone, no atomic sub-SCF and no history-dependent
! heuristic) rather than this build's default GWH, so the SCF this
! program's own gradient is taken at is reproducible run to run without
! needing a committed initial density matrix at all.
!
! One ctest entry per case (test/CMakeLists.txt's own `foreach(case ...)`):
! each invocation of this executable gets a different xyz/basis/functional/
! grid_level/dump/case_name on its command line and runs exactly that one
! case -- never looped over internally.
!
module check_dft_grad_vs_pyscf_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use trc_sph, only: sph_cart_to_sph
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_report
   implicit none
   private
   public :: collect_dft_grad_vs_pyscf_tests, check_dft_grad_vs_pyscf_set_args

   real(dp), parameter :: THRESH = 1.0e-12_dp

   character(len=256) :: g_xyzfile, g_basis_json, g_functional, g_dump, g_case_name
   integer :: g_grid_level

contains

   subroutine check_dft_grad_vs_pyscf_set_args(xyzfile, basis_json, functional, grid_level, dump, case_name)
      character(*), intent(in) :: xyzfile, basis_json, functional, dump, case_name
      integer, intent(in) :: grid_level
      g_xyzfile = xyzfile; g_basis_json = basis_json; g_functional = functional
      g_grid_level = grid_level; g_dump = dump; g_case_name = case_name
   end subroutine check_dft_grad_vs_pyscf_set_args

   subroutine collect_dft_grad_vs_pyscf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("KS-DFT gradient vs pyscf grad.RKS/UKS", test_dft_grad_vs_pyscf)]
   end subroutine collect_dft_grad_vs_pyscf_tests

   subroutine test_dft_grad_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nat, nao, nao_sp, nspin, nalpha, nbeta, unit, dev
      logical :: is_sph
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), grad(:, :)
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_xc_functional_t) :: func
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :), s_sp(:, :)
      logical :: ok
      integer :: u, nat_ref, nao_ref, nfail
      real(dp), allocatable :: grad_ref(:, :)
      real(dp) :: worst
      logical :: shape_ok

      ok = xc_functional_by_name(trim(g_functional), func)
      if (.not. ok .or. len_trim(func%missing) > 0) then
         print '(a)', 'check_dft_grad_vs_pyscf: functional not available'; stop 1
      end if

      dev = trc_bind_device(0)
      call read_xyz(trim(g_xyzfile), nat, at_z, r)

      ! No force_cart here -- see this program's own header: the basis loads
      ! in whatever space basis_json tags it, which is the whole point.
      call trc_basis_from_json(trim(g_basis_json), nat, at_z, r, bas, err)
      if (err%has_error()) then
         print '(a)', 'check_dft_grad_vs_pyscf: '//err%get_message(); stop 1
      end if
      nao = bas%nao
      nao_sp = merge(bas%nao_sph, bas%nao, bas%spherical)
      ! Captured BEFORE bas%release() below, which resets bas%spherical
      ! (and bas%nao/bas%nao_sph) to its type default -- the dump's own
      ! "is_spherical" flag has to be read while `bas` is still alive.
      is_sph = bas%spherical
      call bas%to_device()

      ! The electron count and spin are implied by the molecule, the same
      ! case_name-sniffing convention check_xc_grad_vs_pyscf.F90 uses.
      ! `nh2` is UKS (doublet, 5 alpha/4 beta -- same count as the OH radical
      ! other tests already use); everything else here is closed-shell RKS.
      opts%guess = "gwh"
      if (index(trim(g_case_name), 'nh2') == 1) then
         nalpha = 5; nbeta = 4
         ! the FIXED initial density -- see this program's own header.
         opts%guess = "core"
      else if (index(trim(g_case_name), 'formamide') == 1) then
         nalpha = 12; nbeta = 12
      else
         nalpha = 5; nbeta = 5
      end if

      opts%functional = trim(g_functional)
      opts%grid_level = g_grid_level
      ! conv_diis tightened from 1e-9: CPU and GPU DIIS error floors differ
      ! enough at 1e-9 that the two hardware's converged densities themselves
      ! differ by ~1e-8/1e-9 (measured), which this program's comparison
      ! cannot tell apart from a real gradient bug -- see this program's
      ! own gate comment below.
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-11_dp; opts%max_iteration = 200
      call trc_scf_run(bas, nalpha, nbeta, opts, res)
      call check(error, res%converged, "check_dft_grad_vs_pyscf: reference SCF did not converge")
      if (allocated(error)) then
         call bas%release()
         return
      end if
      nspin = res%nspin

      allocate (grad(3, nat))
      call lp%init(nao_sp)
      ! `nlc_grid_level=grid_level`, NOT `opts%nlc_grid_level` (the SCF's own,
      ! deliberately coarse, default level 1): this comparison is against
      ! PYSCF's INDEPENDENT reference, at a grid PYSCF controls entirely on
      ! its own side (its own `nlcgrids`, not synced to terco's choice at
      ! all) -- there is no "must match the SCF's own grid" requirement the
      ! way check_xc_grad_fd.F90's FD-of-terco's-own-energy comparison has
      ! (that program's own comment). What matters HERE is evaluating the
      ! Pulay+response integral accurately; `opts%nlc_grid_level`'s coarse
      ! default measures 9.37e-7 against this exact pyscf reference (level 1
      ! quadrature error in the VV10 term, not a formula bug --
      ! check_vv10_grad_fixedd_fd.F90 independently confirms the formula to
      ! ~1e-11 at both level 1 AND level 2), while the semilocal grid's own
      ! level (3, same as the reference was generated with) reproduces the
      ! original ~1e-12/1e-13. Explicit here rather than left to
      ! vv10_grad_add's own default so a future change to that default
      ! cannot silently change what this test measures.
      call trc_scf_gradient(bas, lp, nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad, &
                            func=func, grid_level=g_grid_level, grid_max_pts=opts%grid_max_pts, &
                            rho_tol=opts%rho_tol, grid_response=.false., nlc_grid_level=g_grid_level)
      call lp%release()

      ! The overlap, in the SCF's OWN AO space (nao_sp), for the AO-
      ! permutation/normalisation check generate.py's gen_dft_grad runs --
      ! needed in BOTH the Cartesian and the spherical case: terco's
      ! Cartesian d/f off-diagonal components are not unit-normalised the
      ! way pyscf's cart=True ones are (every other *_vs_pyscf generator's
      ! own `ao_norm_from_overlap`), and `sph_cart_to_sph`'s transform
      ! inherits that same per-component scale into the spherical AOs it
      ! builds FROM those Cartesian ones (verified directly: terco's own
      ! S_sph diagonal for an l=2 shell here is ~0.398, not 1 -- a uniform,
      ! m-independent factor per shell, consistent with a Cartesian-input
      ! scale riding through the T^T S T contraction rather than a
      ! genuinely different spherical-harmonic convention). So this is NOT
      ! the "no renormalisation needed" case terco_to_pyscf_perm_sph's own
      ! docstring describes for a committed DM0 permutation -- that claim
      ! holds for permuting an ALREADY-pyscf-side density (gen_full_grad's
      ! own use), not for carrying a terco-side C into pyscf's normalisation,
      ! which is what this program's own comparison does.
      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      call pl%build(bas, THRESH)
      call pl%to_device()
      !$acc enter data create(smat, tmat, vmat)
      call trc_1e(bas, pl, smat, tmat, vmat)
      !$acc update self(smat)
      !$acc exit data delete(smat, tmat, vmat)
      call pl%release()
      allocate (s_sp(nao_sp, nao_sp))
      if (is_sph) then
         call sph_cart_to_sph(bas, smat, s_sp)
      else
         s_sp = smat
      end if
      call bas%release()

      open (newunit=unit, file=trim(g_dump), access='stream', form='unformatted', status='replace')
      write (unit) int(nat, int64), int(nao_sp, int64), int(nspin, int64), &
         int(nalpha, int64), int(nbeta, int64), int(merge(1, 0, is_sph), int64)
      write (unit) int(at_z, int64)
      write (unit) r
      write (unit) int(size(res%grid_r, 2), int64)
      write (unit) res%grid_r
      write (unit) res%grid_w
      write (unit) res%cmo
      write (unit) res%eps
      write (unit) grad
      write (unit) s_sp
      close (unit)
      print '(a)', 'check_dft_grad_vs_pyscf: wrote '//trim(g_dump)//' for the pyscf-reference comparison'

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      call ref_require('dft_grad_'//trim(g_case_name)//'.ref', u)
      call ref_read_int(u, 'NAT', nat_ref)
      call ref_read_int(u, 'NAO', nao_ref)
      shape_ok = nat_ref == nat .and. nao_ref == nao_sp
      call check(error, shape_ok, "check_dft_grad_vs_pyscf: pyscf reference shape does not match this build")
      if (allocated(error)) then
         call ref_close(u)
         deallocate (at_z, r, grad, smat, tmat, vmat, s_sp)
         return
      end if
      allocate (grad_ref(3, nat))
      call ref_read_matrix(u, 'GRAD', 3, nat, grad_ref)
      call ref_close(u)
      worst = maxval(abs(grad - grad_ref))
      ! 1e-9, not the old 3e-6: with conv_diis tightened above, CPU and
      ! GPU converge to the same density (opts%guess="core" for NH2's
      ! fixed initial density; "gwh" otherwise lands on the same minimum
      ! independently, same reasoning as every other *_vs_pyscf check in
      ! this tree) and measured worst|diff| across all 10 cases this
      ! program's own CMakeLists.txt registers, BOTH hardware (gfortran
      ! CPU, nvfortran/V100 GPU): 1.66e-12..4.14e-12 (every pure-GGA/
      ! global-hybrid case, RKS and UKS, cc-pVDZ and 6-31G*) and
      ! 1.99e-13..3.53e-12 (the two wB97X cases) -- CPU and GPU agree to
      ! within ~1e-12 of EACH OTHER on every case, not just within this
      ! gate. 1e-9 is ~250x the worst of those 20 numbers: tight enough
      ! that a genuine ~1e-7 gradient bug (this check's whole job) is
      ! not hidden in the margin the old 3e-6 gate left, but not so
      ! tight it chases hardware-dependent summation-order noise.
      nfail = 0
      call ref_report(nfail, worst < 1.0e-9_dp, 'KS-DFT gradient vs pyscf, worst|diff|', worst)
      call check(error, worst < 1.0e-9_dp, "KS-DFT gradient vs pyscf, worst|diff| = "// &
                 to_string(worst, "es10.2")//" exceeds gate 1.0e-9")
      deallocate (at_z, r, grad, smat, tmat, vmat, s_sp, grad_ref)
   end subroutine test_dft_grad_vs_pyscf

end module check_dft_grad_vs_pyscf_suite

program check_dft_grad_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_dft_grad_vs_pyscf_suite, only: collect_dft_grad_vs_pyscf_tests, check_dft_grad_vs_pyscf_set_args
   implicit none
   integer :: stat, grid_level
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=256) :: xyzfile, basis_json, functional, dump, case_name, glevel_str

   if (command_argument_count() < 6) then
      print '(a)', 'usage: check_dft_grad_vs_pyscf <xyz> <basis.json> <functional> '// &
         '<grid_level> <dump.bin> <case_name>'; stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basis_json)
   call get_command_argument(3, functional)
   call get_command_argument(4, glevel_str)
   call get_command_argument(5, dump)
   call get_command_argument(6, case_name)
   read (glevel_str, *) grid_level

   call check_dft_grad_vs_pyscf_set_args(trim(xyzfile), trim(basis_json), trim(functional), grid_level, &
                                         trim(dump), trim(case_name))

   stat = 0
   testsuites = [new_testsuite("check_dft_grad_vs_pyscf", collect_dft_grad_vs_pyscf_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop 1
   end if
end program check_dft_grad_vs_pyscf
