module check_tddft_suite
   !!
   !! TDA (CIS for Hartree-Fock) singlet excitation energies and oscillator
   !! strengths of a closed-shell reference through trc_tddft, against pyscf's
   !! tdscf TDA (test/pyscf_refs/generate.py's gen_tddft).
   !!
   !! One (molecule, basis, functional) per process, one SCF and one solve
   !! shared by the cases below (built by the first, released by the last):
   !!  1. omega and oscillator strengths against the committed pyscf
   !!     reference, every root; the measured deviations and their margins to
   !!     the tolerances are logged;
   !!  2. the residual ||A X - w X|| of every root, from a FRESH operator and
   !!     one fresh matvec (not the solver's own subspace bookkeeping), below
   !!     the solver's conv_tol;
   !!  3. the amplitudes are orthonormal across roots;
   !!  4. negative control (first root only): the same solve with the exact-exchange weight
   !!     corrupted (110% of the functional's, or 0.1 for a pure functional)
   !!     moves omega past the pyscf tolerance, so the comparison in 1 can
   !!     fail. Skipped for Hartree-Fock, which has no weight to corrupt;
   !!  5. every smaller nroots = 1..N-1 returns the first nroots reference roots
   !!     (a Davidson that skips a root returns a clean-looking wrong set). The
   !!     formamide Hartree-Fock sweep is THE GUARD for the skipped-root bug
   !!     (nroots = 2 once returned roots 1 and 3; the old solver fails exactly
   !!     there, and passes the methane cases). Methane (Td, run for HF and
   !!     B3LYP only, to save CI time) covers the degenerate-multiplet
   !!     handling instead: the cut falls inside degenerate triplets, whose
   !!     individual oscillator strengths are arbitrary, so there f is
   !!     compared only as the sum over the first (complete) triplet. For
   !!     formamide Kohn-Sham one functional (SVWN) sweeps nroots = 2, 3, 5;
   !!  6. refusals (triplet, RPA, UHF, meta-GGA, VV10, range-separated and
   !!     double hybrids) and a solve that is not given enough iterations
   !!     reports an error rather than returning.
   !!
   !! The grid is terco's own (the SCF's), written with terco's numbers to the
   !! probe file for pyscf.
   !!
   !! Tolerances: omega 1e-6 Eh and f 1e-5. The two codes agree to 1e-9 (omega)
   !! and 6e-8 (f) in every case measured, so the gates leave two to three
   !! orders of margin for the platform and grid-sum differences, while a wrong
   !! exchange or kernel weight moves omega by 1e-3 to 1e-1 (the negative
   !! control). The solver runs at conv_tol 1e-7 (omega error ~ conv_tol^2/gap);
   !! the reference SCFs at conv_tol 1e-12 with |g| < 1e-8 (measured 3e-10 to
   !! 8e-9), so SCF noise is far below the gates.
   !!
   !! Arguments: xyz basis.json molname basname funcname nroots probe.bin.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name, XC_FAMILY_MGGA
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_tddft, only: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_op_t, trc_tddft_solve_rks, &
                        trc_tddft_refuse_reason
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_vector
   implicit none
   private
   public :: collect_check_tddft, xyz_arg, bas_json_arg, molname_arg, basname_arg, funcname_arg, nroots_arg, probe_arg

   real(dp), parameter :: ERI_THRESH = 1.0e-12_dp
   real(dp), parameter :: TOL_OMEGA = 1.0e-6_dp
   real(dp), parameter :: TOL_OSC = 1.0e-5_dp
   real(dp), parameter :: CONV_TOL = 1.0e-7_dp

   character(len=256) :: xyz_arg = 'water.xyz'
   character(len=256) :: bas_json_arg = 'basis_sets/6-31gs.json'
   character(len=16) :: molname_arg = 'water'
   character(len=16) :: basname_arg = '6-31gs'
   character(len=16) :: funcname_arg = 'hf'
   integer :: nroots_arg = 4
   character(len=256) :: probe_arg = 'tddft_probe.bin'

   type :: fixture_t
      logical :: ready = .false.
      logical :: dft = .false.
      integer :: nocc = 0, natm = 0
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_linalg_t) :: la
      type(trc_xc_grid_t) :: xg
      type(trc_xc_functional_t) :: func
      type(trc_scf_result_t) :: scf
      type(trc_tddft_options_t) :: opts
      type(trc_tddft_result_t) :: td
      real(dp) :: t_solve = 0.0_dp
   end type fixture_t

   type(fixture_t), save :: fx

contains

   subroutine collect_check_tddft(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("TDA omega and f vs pyscf", test_vs_pyscf), &
                  new_unittest("residual of a fresh matvec below conv_tol", test_residual), &
                  new_unittest("amplitudes orthonormal across roots", test_orthonormal), &
                  new_unittest("negative control: corrupted exchange weight fails the comparison", test_control), &
                  new_unittest("every nroots = 1..N-1 gives the first roots of the reference", test_sweep), &
                  new_unittest("a tiny max_space (repeated collapse) gives the same roots", test_collapse), &
                  new_unittest("refusals and non-convergence are errors", test_refusals), &
                  new_unittest("release", test_release) &
                  ]
   end subroutine collect_check_tddft

   subroutine fixture_build(error)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(trc_scf_options_t) :: sopts
      integer :: dev, natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      logical :: ok
      type(timer_type) :: tm

      if (fx%ready) return
      dev = trc_bind_device(0)
      fx%dft = trim(funcname_arg) /= 'hf'
      if (fx%dft) then
         ok = xc_functional_by_name(trim(funcname_arg), fx%func)
         call check(error, ok, "unrecognised functional '"//trim(funcname_arg)//"'")
         if (allocated(error)) return
      end if
      call read_xyz(trim(xyz_arg), natm, zint, at_r)
      fx%natm = natm
      call trc_basis_from_json(trim(bas_json_arg), natm, zint, at_r, fx%bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call fx%bas%to_device()
      call fx%pl%build(fx%bas, ERI_THRESH)
      call fx%pl%to_device()
      call fx%eri%build(fx%bas, ERI_THRESH)
      call fx%la%init(fx%bas%nao)
      fx%ready = .true.

      fx%nocc = nint(sum(real(zint, dp)))/2
      sopts%conv_energy = 1.0e-12_dp
      sopts%conv_density = 1.0e-10_dp
      sopts%conv_diis = 1.0e-10_dp
      sopts%max_iteration = 200
      if (fx%dft) sopts%functional = trim(funcname_arg)
      sopts%grid_level = 3
      call trc_scf_run(fx%bas, fx%nocc, fx%nocc, sopts, fx%scf)
      call check(error, fx%scf%converged, "SCF did not converge")
      if (allocated(error)) return
      if (fx%dft) then
         call fx%xg%build(size(fx%scf%grid_w), fx%scf%grid_r, fx%scf%grid_w, fx%bas, max_pts=sopts%grid_max_pts)
         call fx%xg%to_device()
      end if

      fx%opts%nroots = nroots_arg
      fx%opts%conv_tol = CONV_TOL
      call tm%start()
      call solve(fx%opts, fx%td, err)
      fx%t_solve = tm%get_elapsed_time()
      call check(error, .not. err%has_error(), "TDA solve failed: "//err%get_message())
   end subroutine fixture_build

   subroutine solve(opts, td, err, func)
      !! The solver at the fixture's SCF; `func` overrides the fixture's
      !! functional (the negative control).
      type(trc_tddft_options_t), intent(in) :: opts
      type(trc_tddft_result_t), intent(out) :: td
      type(error_t), intent(out) :: err
      type(trc_xc_functional_t), intent(in), optional :: func
      type(trc_xc_functional_t) :: f
      f = fx%func
      if (present(func)) f = func
      if (fx%dft) then
         call trc_tddft_solve_rks(fx%bas, fx%pl, fx%eri, fx%la, fx%nocc, fx%scf%cmo(:, :, 1), fx%scf%eps(:, 1), &
                                  fx%scf%dmat(:, :, 1), opts, td, err, func=f, xg=fx%xg)
      else
         call trc_tddft_solve_rks(fx%bas, fx%pl, fx%eri, fx%la, fx%nocc, fx%scf%cmo(:, :, 1), fx%scf%eps(:, 1), &
                                  fx%scf%dmat(:, :, 1), opts, td, err)
      end if
   end subroutine solve

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: om_ref(:), osc_ref(:)
      real(dp) :: e_ref, dom, dosc
      integer :: u, nr, unit
      character(len=256) :: refname

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%td%converged, "TDA did not converge")
      if (allocated(error)) return
      nr = nroots_arg

      open (newunit=unit, file=trim(probe_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(fx%bas%nao, int64), int(fx%natm, int64), int(fx%nocc, int64), int(nr, int64)
      write (unit) molname_arg, basname_arg, funcname_arg
      write (unit) fx%scf%energy
      write (unit) fx%td%omega
      write (unit) fx%td%osc
      write (unit) int(size(fx%scf%grid_w), int64)
      write (unit) fx%scf%grid_r, fx%scf%grid_w
      close (unit)

      refname = 'tddft_'//trim(molname_arg)//'_'//trim(basname_arg)//'_'//trim(funcname_arg)//'.ref'
      call ref_require(trim(refname), u)
      allocate (om_ref(nr), osc_ref(nr))
      call ref_read_real(u, 'E_SCF', e_ref)
      call ref_read_vector(u, 'OMEGA', nr, om_ref)
      call ref_read_vector(u, 'OSC', nr, osc_ref)
      call ref_close(u)

      call check(error, abs(fx%scf%energy - e_ref) < 1.0e-9_dp, &
                 "E vs pyscf, |dE| = "//to_string(abs(fx%scf%energy - e_ref), "es10.2"))
      if (allocated(error)) return
      dom = maxval(abs(fx%td%omega - om_ref))
      dosc = maxval(abs(fx%td%osc - osc_ref))
      call global_logger%info("tddft "//trim(molname_arg)//"/"//trim(funcname_arg)//": "//to_string(nr)// &
                              " roots, max|d omega| "//to_string(dom, "es10.2")//" (tol "//to_string(TOL_OMEGA, "es8.1")// &
                              ", margin x"//to_string(TOL_OMEGA/max(dom, 1.0e-30_dp), "es8.1")//"), max|d f| "// &
                              to_string(dosc, "es10.2")//" (tol "//to_string(TOL_OSC, "es8.1")//", margin x"// &
                              to_string(TOL_OSC/max(dosc, 1.0e-30_dp), "es8.1")//"); "//to_string(fx%td%iterations)// &
                              " iterations, "//to_string(fx%td%nmatvec)//" matvecs, solve "// &
                              to_string(fx%t_solve, "f8.2")//" s; omega(1) "//to_string(fx%td%omega(1), "f12.8"))
      call check(error, dom < TOL_OMEGA, "omega vs pyscf, max|diff| = "//to_string(dom, "es10.2"))
      if (allocated(error)) return
      if (degenerate_molecule()) then
         ! A degenerate triplet's individual f depend on the arbitrary rotation
         ! inside it; the triplet's sum does not.
         dosc = abs(sum(fx%td%osc(1:3)) - sum(osc_ref(1:3)))
      end if
      call check(error, dosc < TOL_OSC, "oscillator strength vs pyscf, max|diff| = "//to_string(dosc, "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(fx%td%osc) > 1.0e-4_dp .or. nr < 1, "no root has a measurable oscillator strength")
   end subroutine test_vs_pyscf

   logical function degenerate_molecule()
      degenerate_molecule = trim(molname_arg) == 'ch4'
   end function degenerate_molecule

   logical function skip_extra()
      !! The sweep and the collapse solve repeat full solves: not for every
      !! functional of every molecule.
      skip_extra = .false.
      if (degenerate_molecule()) skip_extra = trim(funcname_arg) == 'svwn' .or. trim(funcname_arg) == 'pbe'
      if (fx%dft .and. fx%nocc > 6) skip_extra = trim(funcname_arg) /= 'svwn'
   end function skip_extra

   subroutine test_sweep(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_options_t) :: o
      type(trc_tddft_result_t) :: td
      type(error_t) :: err
      real(dp), allocatable :: om_ref(:), osc_ref(:)
      real(dp) :: e_ref, worst, d
      integer :: u, nr, k
      character(len=256) :: refname

      call fixture_build(error)
      if (allocated(error)) return
      nr = nroots_arg
      ! Cost: formamide Kohn-Sham solves take 10 to 20 s each, so only SVWN
      ! (the cheapest kernel) sweeps, over nroots = 2, 3, 5; methane runs the
      ! sweep for HF and B3LYP only.
      if (skip_extra()) then
         call global_logger%info("tddft: nroots sweep skipped for "//trim(molname_arg)//"/"//trim(funcname_arg)// &
                                 " (cost)")
         return
      end if
      refname = 'tddft_'//trim(molname_arg)//'_'//trim(basname_arg)//'_'//trim(funcname_arg)//'.ref'
      call ref_require(trim(refname), u)
      allocate (om_ref(nr), osc_ref(nr))
      call ref_read_real(u, 'E_SCF', e_ref)
      call ref_read_vector(u, 'OMEGA', nr, om_ref)
      call ref_read_vector(u, 'OSC', nr, osc_ref)
      call ref_close(u)
      worst = 0.0_dp
      do k = 1, nr - 1
         if (fx%dft .and. fx%nocc > 6 .and. all(k /= [2, 3, 5])) cycle
         o = fx%opts
         o%nroots = k
         call solve(o, td, err)
         call check(error, .not. err%has_error(), "nroots = "//to_string(k)//" failed: "//err%get_message())
         if (allocated(error)) return
         d = maxval(abs(td%omega - om_ref(1:k)))
         worst = max(worst, d)
         call check(error, d < TOL_OMEGA, "nroots = "//to_string(k)//" does not return the first "//to_string(k)// &
                    " reference roots, max|d omega| = "//to_string(d, "es10.2"))
         if (allocated(error)) return
      end do
      call global_logger%info("tddft "//trim(molname_arg)//"/"//trim(funcname_arg)//": nroots sweep 1.."// &
                              to_string(nr - 1)//", worst max|d omega| "//to_string(worst, "es10.2"))
   end subroutine test_sweep

   subroutine test_residual(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_op_t) :: op
      type(error_t) :: err
      real(dp), allocatable :: ax(:, :, :)
      real(dp) :: rk, worst
      integer :: k, nr

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%td%converged, "TDA did not converge")
      if (allocated(error)) return
      nr = nroots_arg
      if (fx%dft) then
         call op%setup(fx%bas, fx%la, fx%nocc, fx%scf%cmo(:, :, 1), fx%scf%eps(:, 1), fx%scf%dmat(:, :, 1), err, &
                       func=fx%func, xg=fx%xg)
      else
         call op%setup(fx%bas, fx%la, fx%nocc, fx%scf%cmo(:, :, 1), fx%scf%eps(:, 1), fx%scf%dmat(:, :, 1), err)
      end if
      call check(error, .not. err%has_error(), "operator setup: "//err%get_message())
      if (allocated(error)) return
      allocate (ax(fx%nocc, op%nvir, nr))
      call op%apply(fx%bas, fx%eri, fx%la, nr, fx%td%xamp, ax, fx%xg)
      worst = 0.0_dp
      do k = 1, nr
         rk = sqrt(sum((ax(:, :, k) - fx%td%omega(k)*fx%td%xamp(:, :, k))**2))
         worst = max(worst, rk)
      end do
      call op%release()
      call global_logger%info("tddft "//trim(molname_arg)//"/"//trim(funcname_arg)//": fresh-matvec max residual "// &
                              to_string(worst, "es10.2")//" (conv_tol "//to_string(CONV_TOL, "es8.1")//")")
      call check(error, worst < CONV_TOL, "fresh residual "//to_string(worst, "es10.2")//" >= conv_tol")
   end subroutine test_residual

   subroutine test_orthonormal(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: dev
      integer :: j, k
      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%td%converged, "TDA did not converge")
      if (allocated(error)) return
      dev = 0.0_dp
      do k = 1, nroots_arg
         do j = 1, nroots_arg
            dev = max(dev, abs(sum(fx%td%xamp(:, :, j)*fx%td%xamp(:, :, k)) - merge(1.0_dp, 0.0_dp, j == k)))
         end do
      end do
      call global_logger%info("tddft "//trim(molname_arg)//"/"//trim(funcname_arg)//": max|X^T X - 1| "// &
                              to_string(dev, "es10.2"))
      call check(error, dev < 1.0e-10_dp, "X not orthonormal, max dev "//to_string(dev, "es10.2"))
      if (allocated(error)) return
      do k = 2, nroots_arg
         call check(error, fx%td%omega(k) >= fx%td%omega(k - 1), "omega not ascending")
         if (allocated(error)) return
      end do
   end subroutine test_orthonormal

   subroutine test_control(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_xc_functional_t) :: bad
      type(trc_tddft_options_t) :: o1
      type(trc_tddft_result_t) :: td
      type(error_t) :: err
      real(dp), allocatable :: om_ref(:), osc_ref(:)
      real(dp) :: e_ref, dom
      integer :: u, nr
      character(len=256) :: refname

      call fixture_build(error)
      if (allocated(error)) return
      if (.not. fx%dft) then
         call global_logger%info("tddft: negative control skipped for Hartree-Fock (no exchange weight)")
         return
      end if
      nr = nroots_arg
      bad = fx%func
      if (fx%func%exx > 0.0_dp) then
         bad%exx = 1.1_dp*fx%func%exx
      else
         bad%exx = 0.1_dp
      end if
      o1 = fx%opts
      o1%nroots = 1
      call solve(o1, td, err, func=bad)
      call check(error, .not. err%has_error(), "corrupted solve failed: "//err%get_message())
      if (allocated(error)) return
      refname = 'tddft_'//trim(molname_arg)//'_'//trim(basname_arg)//'_'//trim(funcname_arg)//'.ref'
      call ref_require(trim(refname), u)
      allocate (om_ref(nr), osc_ref(nr))
      call ref_read_real(u, 'E_SCF', e_ref)
      call ref_read_vector(u, 'OMEGA', nr, om_ref)
      call ref_read_vector(u, 'OSC', nr, osc_ref)
      call ref_close(u)
      dom = abs(td%omega(1) - om_ref(1))
      call global_logger%info("tddft "//trim(molname_arg)//"/"//trim(funcname_arg)//": negative control (c_x "// &
                              to_string(fx%func%exx, "f6.3")//" -> "//to_string(bad%exx, "f6.3")//") moves omega by "// &
                              to_string(dom, "es10.2")//" (tolerance "//to_string(TOL_OMEGA, "es8.1")//")")
      call check(error, dom > 100.0_dp*TOL_OMEGA, "a corrupted exchange weight must fail the omega comparison, "// &
                 "moved only "//to_string(dom, "es10.2"))
   end subroutine test_control

   subroutine test_collapse(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_options_t) :: o
      type(trc_tddft_result_t) :: td
      type(error_t) :: err
      real(dp) :: dom
      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%td%converged, "TDA did not converge")
      if (allocated(error)) return
      ! A second full solve: not on formamide, and for methane HF and B3LYP only.
      if (nroots_arg > 5 .or. skip_extra()) return
      o = fx%opts
      o%max_space = 3*nroots_arg
      o%max_iter = 400
      call solve(o, td, err)
      call check(error, .not. err%has_error(), "collapsing solve failed: "//err%get_message())
      if (allocated(error)) return
      dom = maxval(abs(td%omega - fx%td%omega))
      call global_logger%info("tddft "//trim(molname_arg)//"/"//trim(funcname_arg)//": max_space "// &
                              to_string(3*nroots_arg)//" took "//to_string(td%iterations)//" iterations (vs "// &
                              to_string(fx%td%iterations)//"), max|d omega| "//to_string(dom, "es10.2"))
      call check(error, dom < 1.0e-8_dp, "collapsed solve differs, max|d omega| "//to_string(dom, "es10.2"))
   end subroutine test_collapse

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_options_t) :: o
      type(trc_xc_functional_t) :: f
      type(trc_tddft_result_t) :: td
      type(error_t) :: err

      ! Pure refusal logic.
      o%singlet = .false.
      call check_refused(error, trc_tddft_refuse_reason(o), "triplet", "singlet")
      if (allocated(error)) return
      o = trc_tddft_options_t()
      o%tda = .false.
      call check_refused(error, trc_tddft_refuse_reason(o), "RPA", "Tamm-Dancoff")
      if (allocated(error)) return
      o = trc_tddft_options_t()
      call check_refused(error, trc_tddft_refuse_reason(o, nspin=2), "UHF/UKS", "closed-shell")
      if (allocated(error)) return
      f%family = XC_FAMILY_MGGA
      call check_refused(error, trc_tddft_refuse_reason(o, f), "meta-GGA", "meta-GGA")
      if (allocated(error)) return
      f = trc_xc_functional_t()
      f%nlc_b = 6.0_dp
      call check_refused(error, trc_tddft_refuse_reason(o, f), "VV10", "VV10")
      if (allocated(error)) return
      f = trc_xc_functional_t()
      f%omega = 0.3_dp
      call check_refused(error, trc_tddft_refuse_reason(o, f), "range-separated hybrid", "range-separated")
      if (allocated(error)) return
      f = trc_xc_functional_t()
      f%pt2_cos = 0.25_dp
      call check_refused(error, trc_tddft_refuse_reason(o, f), "double hybrid", "double hybrid")
      if (allocated(error)) return
      f = trc_xc_functional_t()
      call check(error, len_trim(trc_tddft_refuse_reason(o, f)) == 0, "a plain LDA/GGA must not be refused")
      if (allocated(error)) return
      call check(error, len_trim(trc_tddft_refuse_reason(o)) == 0, "Hartree-Fock must not be refused")
      if (allocated(error)) return

      ! Through the solver: an error and a non-converged result, no work done.
      call fixture_build(error)
      if (allocated(error)) return
      o = fx%opts
      o%tda = .false.
      call solve(o, td, err)
      call check(error, err%has_error() .and. .not. td%converged, "tda=.false. must be an error")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "Tamm-Dancoff") > 0, "refusal message: "//err%get_message())
      if (allocated(error)) return

      ! More roots than occupied-virtual pairs.
      o = fx%opts
      o%nroots = fx%nocc*(size(fx%scf%eps, 1) - fx%nocc) + 1
      call solve(o, td, err)
      call check(error, err%has_error() .and. .not. td%converged, "nroots > nocc*nvir must be an error")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "exceeds") > 0, "nroots message: "//err%get_message())
      if (allocated(error)) return

      ! One iteration cannot converge: the failure must be reported.
      o = fx%opts
      o%nroots = 1
      o%max_iter = 1
      call solve(o, td, err)
      call check(error, err%has_error() .and. .not. td%converged, "running out of iterations must be an error")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "not converged") > 0, "non-convergence message: "//err%get_message())
   end subroutine test_refusals

   subroutine check_refused(error, msg, what, substring)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: msg, what, substring
      call check(error, len_trim(msg) > 0, what//" must be refused")
      if (allocated(error)) return
      call check(error, index(msg, substring) > 0, what//" message must mention '"//substring//"'")
   end subroutine check_refused

   subroutine test_release(error)
      type(error_type), allocatable, intent(out) :: error
      if (fx%ready) then
         call fx%eri%release()
         call fx%la%release()
         if (fx%dft) call fx%xg%release()
         call fx%pl%release()
         call fx%bas%release()
         fx%ready = .false.
      end if
      if (.false.) allocate (error)
   end subroutine test_release

end module check_tddft_suite

program check_tddft
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_tddft_suite, only: collect_check_tddft, xyz_arg, bas_json_arg, molname_arg, basname_arg, funcname_arg, &
                                 nroots_arg, probe_arg
   implicit none
   integer :: stat
   character(len=16) :: buf

   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 2) call get_command_argument(2, bas_json_arg)
   if (command_argument_count() >= 3) call get_command_argument(3, molname_arg)
   if (command_argument_count() >= 4) call get_command_argument(4, basname_arg)
   if (command_argument_count() >= 5) call get_command_argument(5, funcname_arg)
   if (command_argument_count() >= 6) then
      call get_command_argument(6, buf)
      read (buf, *) nroots_arg
   end if
   if (command_argument_count() >= 7) call get_command_argument(7, probe_arg)

   stat = 0
   call run_testsuite(collect_check_tddft, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_tddft
