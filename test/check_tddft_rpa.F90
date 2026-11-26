module check_tddft_rpa_suite
   !!
   !! Full RPA (TDHF for Hartree-Fock, Casida for Kohn-Sham) singlet excitation
   !! energies and oscillator strengths of a closed-shell reference through
   !! trc_tddft (opts%tda = .false.), against pyscf's tdscf TDHF / TDDFT
   !! (test/pyscf_refs/generate.py's gen_rpa). The TDA's own suite is
   !! check_tddft.F90; this is its RPA counterpart and shares no reference file
   !! with it.
   !!
   !! One (molecule, basis, functional) per process, one SCF and one RPA solve
   !! shared by the cases below (built by the first, released by the last):
   !!  1. omega, oscillator strengths and ||Y||^2 against the committed pyscf
   !!     reference, every root; the measured deviations and their margins to
   !!     the tolerances are logged. ||Y||^2 is compared with pyscf's (its
   !!     (X, Y) are normalised 2 (X.X - Y.Y) = 1, so ||Y||^2(terco) =
   !!     2 ||Y||^2(pyscf)); a degenerate multiplet is compared as a sum;
   !!  2. the residual of the (X, Y) equation [[A, B], [-B, -A]] of every
   !!     root, from a FRESH operator and a fresh (A+B) and (A-B) matvec on
   !!     X+Y and X-Y (not the solver's own subspace bookkeeping), below the
   !!     solver's conv_tol;
   !!  3. X and Y are normalised X.X - Y.Y = 1 for every pair of roots
   !!     (orthonormality in the RPA metric);
   !!  4. for a molecule of at most 200 pairs, A+B and A-B are built densely
   !!     from unit-vector matvecs: both are symmetric, A-B has no Coulomb or
   !!     kernel term (for a pure functional it is the diagonal e_a - e_i, with
   !!     the antisymmetric exchange half of the batch forced on), and the
   !!     lowest roots of (A-B)^{1/2} (A+B) (A-B)^{1/2} are the solver's w^2.
   !!     For every molecule, w(RPA) <= w(TDA) root by root, with the TDA roots
   !!     computed in the test by the solver itself (opts%tda = .true.);
   !!  5. negative control (first root only): the same solve with the exact-exchange weight
   !!     corrupted (110% of the functional's, or 0.1 for a pure functional)
   !!     moves omega past the pyscf tolerance, so the comparison in 1 can
   !!     fail. Skipped for Hartree-Fock, which has no weight to corrupt;
   !!  6. every smaller nroots = 1..N-1 returns the first nroots reference roots
   !!     (a Davidson that skips a root returns a clean-looking wrong set).
   !!     Methane (Td, run for HF and B3LYP only, to save CI time) covers the
   !!     degenerate-multiplet handling: the cut falls inside degenerate
   !!     triplets, whose individual oscillator strengths are arbitrary, so
   !!     there f is compared only as the sum over the first (complete) triplet.
   !!     For formamide Kohn-Sham one functional (SVWN) sweeps nroots = 2, 3, 5;
   !!  7. a tiny max_space gives the same roots, and repeated collapses
   !!     happened (the RPA subspace is rebuilt from p and m on a collapse);
   !!  8. refusals (triplet, UHF, meta-GGA, VV10, range-separated and double
   !!     hybrids; the RPA itself is not refused) and a solve that is not
   !!     given enough iterations report an error rather than returning;
   !!  9. the RPA error paths: a corrupted exact-exchange weight (c_x = 2.5)
   !!     makes A-B not positive definite in the Davidson subspace, and a
   !!     kernel scaled to a large negative value (pure functionals) makes
   !!     A+B not positive definite; each is an error naming the matrix and
   !!     never a silent return.
   !!
   !! Instability is detected only in the Davidson subspace, V^T (A-B) V and
   !! V^T (A+B) V, so a reference that is unstable along a direction the
   !! subspace never reaches is not flagged; the error paths above exercise the
   !! detection that exists, not a full-space stability analysis.
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
   use pic_types, only: dp, int64, default_int
   use pic_strings, only: to_string
   use pic_lapack_interfaces, only: pic_syev
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc, only: trc_xc_fxc_vec_chunk
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name, XC_FAMILY_MGGA
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_tddft, only: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_op_t, trc_tddft_solve_rks, &
                        trc_tddft_reference_t, &
                        trc_tddft_refuse_reason
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_vector
   implicit none
   private
   public :: collect_check_tddft_rpa, xyz_arg, bas_json_arg, molname_arg, basname_arg, funcname_arg, nroots_arg, probe_arg

   real(dp), parameter :: ERI_THRESH = 1.0e-12_dp
   real(dp), parameter :: TOL_OMEGA = 1.0e-6_dp
   real(dp), parameter :: TOL_OSC = 1.0e-5_dp
   real(dp), parameter :: CONV_TOL = 1.0e-7_dp
   real(dp), parameter :: TOL_YY = 1.0e-6_dp
   integer, parameter :: MAX_DENSE = 200
   !! Largest number of occupied-virtual pairs of the dense RPA check.

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
      type(trc_tddft_reference_t) :: ref
      real(dp) :: t_solve = 0.0_dp
   end type fixture_t

   type(fixture_t), save :: fx

contains

   subroutine collect_check_tddft_rpa(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("omega and f vs pyscf", test_vs_pyscf), &
                  new_unittest("residual of a fresh matvec below conv_tol", test_residual), &
                  new_unittest("amplitudes orthonormal across roots", test_orthonormal), &
                  new_unittest("w(RPA) <= w(TDA), the TDA computed by the solver", test_below_tda), &
                  new_unittest("dense (A+B) and (A-B) from unit-vector matvecs", test_rpa_dense), &
                  new_unittest("negative control: corrupted exchange weight fails the comparison", test_control), &
                  new_unittest("every nroots = 1..N-1 gives the first roots of the reference", test_sweep), &
                  new_unittest("a tiny max_space (repeated collapse) gives the same roots", test_collapse), &
                  new_unittest("refusals and non-convergence are errors", test_refusals), &
                  new_unittest("RPA: indefinite A-B and A+B are errors", test_rpa_errors), &
                  new_unittest("release", test_release) &
                  ]
   end subroutine collect_check_tddft_rpa

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

      if (fx%dft) then
         call fx%ref%set(fx%nocc, fx%scf%cmo(:, :, 1), fx%scf%eps(:, 1), fx%scf%dmat(:, :, 1), fx%func)
      else
         call fx%ref%set(fx%nocc, fx%scf%cmo(:, :, 1), fx%scf%eps(:, 1), fx%scf%dmat(:, :, 1))
      end if
      fx%opts%nroots = nroots_arg
      fx%opts%tda = .false.
      fx%opts%conv_tol = CONV_TOL
      call tm%start()
      call solve(fx%opts, fx%td, err)
      fx%t_solve = tm%get_elapsed_time()
      call check(error, .not. err%has_error(), "solve failed: "//err%get_message())
   end subroutine fixture_build

   subroutine solve(opts, td, err, func)
      !! The solver at the fixture's SCF; `func` overrides the fixture's
      !! functional (the negative control).
      type(trc_tddft_options_t), intent(in) :: opts
      type(trc_tddft_result_t), intent(out) :: td
      type(error_t), intent(out) :: err
      type(trc_xc_functional_t), intent(in), optional :: func
      type(trc_tddft_reference_t) :: r
      r = fx%ref
      if (present(func)) r%func = func
      call trc_tddft_solve_rks(fx%bas, fx%pl, fx%eri, fx%la, r, opts, td, fx%xg)
      err = td%err
   end subroutine solve

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: om_ref(:), osc_ref(:), yy_ref(:), yy(:)
      real(dp) :: e_ref, dom, dosc, dyy
      integer :: u, nr, unit, k
      character(len=256) :: refname

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%td%converged, "the solve did not converge")
      if (allocated(error)) return
      nr = nroots_arg

      open (newunit=unit, file=trim(probe_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(fx%bas%nao, int64), int(fx%natm, int64), int(fx%nocc, int64), int(nr, int64)
      write (unit) molname_arg, basname_arg, funcname_arg
      write (unit) fx%scf%energy
      write (unit) fx%td%omega
      write (unit) fx%td%osc
      allocate (yy(nr))
      do k = 1, nr
         yy(k) = sum(fx%td%yamp(:, :, k)**2)
      end do
      write (unit) yy
      ! Hartree-Fock leaves the grid unallocated: write a zero count and no arrays.
      if (allocated(fx%scf%grid_w)) then
         write (unit) int(size(fx%scf%grid_w), int64)
         write (unit) fx%scf%grid_r, fx%scf%grid_w
      else
         write (unit) 0_int64
      end if
      close (unit)

      refname = ref_name()
      call ref_require(trim(refname), u)
      allocate (om_ref(nr), osc_ref(nr))
      call ref_read_real(u, 'E_SCF', e_ref)
      call ref_read_vector(u, 'OMEGA', nr, om_ref)
      call ref_read_vector(u, 'OSC', nr, osc_ref)
      allocate (yy_ref(nr))
      call ref_read_vector(u, 'YY', nr, yy_ref)
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
      if (allocated(error)) return
      ! ||Y||^2 is invariant to sign; inside a degenerate multiplet only its sum is.
      if (degenerate_molecule()) then
         dyy = abs(sum(yy(1:3)) - sum(yy_ref(1:3)))
      else
         dyy = maxval(abs(yy - yy_ref))
      end if
      call global_logger%info("tddft rpa "//trim(molname_arg)//"/"//trim(funcname_arg)//": max|d ||Y||^2| "// &
                              to_string(dyy, "es10.2")//" (tol "//to_string(TOL_YY, "es8.1")//", margin x"// &
                              to_string(TOL_YY/max(dyy, 1.0e-30_dp), "es8.1")//"), largest ||Y||^2 "// &
                              to_string(maxval(yy), "es10.2"))
      call check(error, dyy < TOL_YY, "||Y||^2 vs pyscf, max|diff| = "//to_string(dyy, "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(yy) > 1.0e-6_dp, "Y is zero: not an RPA solution")
   end subroutine test_vs_pyscf

   function ref_name() result(name)
      !! The committed pyscf reference of this case.
      character(len=256) :: name
      name = 'rpa_'//trim(molname_arg)//'_'//trim(basname_arg)//'_'//trim(funcname_arg)//'.ref'
   end function ref_name

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
      refname = ref_name()
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
      real(dp), allocatable :: ax2(:, :, :, :), bx2(:, :, :, :), pp(:, :, :), mm(:, :, :)
      real(dp) :: rk, worst
      integer :: k, nr, nchunk

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%td%converged, "the solve did not converge")
      if (allocated(error)) return
      nr = nroots_arg
      call op%setup(fx%bas, fx%la, fx%ref, err, fx%xg)
      call check(error, .not. err%has_error(), "operator setup: "//err%get_message())
      if (allocated(error)) return
      ! (A+B)(X+Y) and (A-B)(X-Y), each from its own fresh batch of nr vectors.
      allocate (ax2(fx%nocc, op%nvir, nr, 2), bx2(fx%nocc, op%nvir, nr, 2), pp(fx%nocc, op%nvir, nr), &
                mm(fx%nocc, op%nvir, nr))
      pp = fx%td%xamp + fx%td%yamp
      mm = fx%td%xamp - fx%td%yamp
      call op%apply_rpa(fx%bas, fx%eri, fx%la, nr, pp, ax2, fx%xg)
      call op%apply_rpa(fx%bas, fx%eri, fx%la, nr, mm, bx2, fx%xg)
      if (fx%dft) then
         ! No budget is set, so the default must keep the batch in ONE fxc
         ! pass: a budget of one double once made it one vector per pass, with
         ! no visible change in the answers. The default is capped by memory
         ! (formamide holds 8 vectors, not 12, on the host), so the check is on
         ! a batch of the wanted roots, which is what every later iteration
         ! sends at most.
         nchunk = trc_xc_fxc_vec_chunk(fx%bas, fx%xg, op%cache, nr)
         call check(error, nchunk == nr, "the default fxc vec chunk is "//to_string(nchunk)//", not one pass for "// &
                    to_string(nr)//" vectors")
         if (allocated(error)) then
            call op%release()
            return
         end if
      end if
      worst = 0.0_dp
      do k = 1, nr
         ! [[A, B], [-B, -A]] (X, Y) - w (X, Y) from (A+B)p and (A-B)m:
         ! A X + B Y = ((A+B)p + (A-B)m)/2, B X + A Y = ((A+B)p - (A-B)m)/2.
         rk = sqrt(sum((0.5_dp*(ax2(:, :, k, 1) + bx2(:, :, k, 2)) - fx%td%omega(k)*fx%td%xamp(:, :, k))**2) &
                   + sum((0.5_dp*(ax2(:, :, k, 1) - bx2(:, :, k, 2)) + fx%td%omega(k)*fx%td%yamp(:, :, k))**2))
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
      call check(error, fx%td%converged, "the solve did not converge")
      if (allocated(error)) return
      dev = 0.0_dp
      do k = 1, nroots_arg
         do j = 1, nroots_arg
            ! The RPA's metric: X.X - Y.Y.
            dev = max(dev, abs(sum(fx%td%xamp(:, :, j)*fx%td%xamp(:, :, k)) - sum(fx%td%yamp(:, :, j)*fx%td%yamp(:, :, k)) &
                               - merge(1.0_dp, 0.0_dp, j == k)))
         end do
      end do
      call global_logger%info("tddft "//trim(molname_arg)//"/"//trim(funcname_arg)//": max|X^T X - Y^T Y - 1| "// &
                              to_string(dev, "es10.2"))
      call check(error, dev < 1.0e-10_dp, "X not orthonormal, max dev "//to_string(dev, "es10.2"))
      if (allocated(error)) return
      do k = 2, nroots_arg
         call check(error, fx%td%omega(k) >= fx%td%omega(k - 1), "omega not ascending")
         if (allocated(error)) return
      end do
   end subroutine test_orthonormal

   subroutine test_below_tda(error)
      !! w(RPA) <= w(TDA) root by root, with the TDA roots from the solver
      !! itself (opts%tda = .true.) at the same SCF: the RPA's Y couples the
      !! TDA's excitations to de-excitations and can only lower the roots.
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_options_t) :: o
      type(trc_tddft_result_t) :: td
      type(error_t) :: err
      real(dp) :: dtda

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%td%converged, "the solve did not converge")
      if (allocated(error)) return
      o = fx%opts
      o%tda = .true.
      call solve(o, td, err)
      call check(error, .not. err%has_error() .and. td%converged, "the TDA solve failed: "//err%get_message())
      if (allocated(error)) return
      dtda = maxval(fx%td%omega - td%omega)
      call global_logger%info("tddft rpa "//trim(molname_arg)//"/"//trim(funcname_arg)// &
                              ": max (w_RPA - w_TDA) "//to_string(dtda, "es10.2")//" (must be <= 0), smallest "// &
                              "TDA-RPA shift "//to_string(minval(td%omega - fx%td%omega), "es10.2"))
      call check(error, dtda <= 1.0e-9_dp, "an RPA root lies above its TDA root, by "//to_string(dtda, "es10.2"))
      if (allocated(error)) return
      call check(error, minval(td%omega - fx%td%omega) > 1.0e-7_dp, "no RPA root is below its TDA root: Y is not coupled in")
   end subroutine test_below_tda

   subroutine test_rpa_dense(error)
      !! At most MAX_DENSE pairs: A+B and A-B dense from unit-vector matvecs;
      !! their symmetry, A-B's structure for a pure functional and the lowest
      !! roots of (A-B)^{1/2} (A+B) (A-B)^{1/2} against the solver.
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_op_t) :: op
      type(error_t) :: err
      real(dp), allocatable :: unit_x(:, :, :), prod(:, :, :, :), apb(:, :), amb(:, :), amh(:, :), work(:, :)
      real(dp), allocatable :: sam(:), w2(:)
      real(dp) :: dsym, doff, dw
      integer :: nov, nocc, nvir, nr, i, j
      integer(default_int) :: linfo
      logical :: is_pure
      character(len=64) :: offmsg

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%td%converged, "the solve did not converge")
      if (allocated(error)) return
      nr = nroots_arg
      nocc = fx%nocc
      nvir = size(fx%scf%eps, 1) - nocc
      nov = nocc*nvir

      if (nov > MAX_DENSE) then
         call global_logger%info("tddft rpa: dense check skipped for "//trim(molname_arg)//" ("//to_string(nov)// &
                                 " pairs, cost)")
         return
      end if
      is_pure = fx%dft .and. fx%func%exx == 0.0_dp
      call op%setup(fx%bas, fx%la, fx%ref, err, fx%xg)
      call check(error, .not. err%has_error(), "operator setup: "//err%get_message())
      if (allocated(error)) return
      ! The antisymmetric half is forced on for a pure functional, so that code
      ! path of the matvec runs here too (its Fock results are zeros by construction).
      op%force_antisym = .true.
      allocate (unit_x(nocc, nvir, nov), prod(nocc, nvir, nov, 2), apb(nov, nov), amb(nov, nov))
      unit_x = 0.0_dp
      do j = 1, nov
         unit_x(mod(j - 1, nocc) + 1, (j - 1)/nocc + 1, j) = 1.0_dp
      end do
      call op%apply_rpa(fx%bas, fx%eri, fx%la, nov, unit_x, prod, fx%xg)
      apb = reshape(prod(:, :, :, 1), [nov, nov])
      amb = reshape(prod(:, :, :, 2), [nov, nov])
      call op%release()

      dsym = max(maxval(abs(apb - transpose(apb))), maxval(abs(amb - transpose(amb))))
      call check(error, dsym < 1.0e-8_dp, "A+B or A-B not symmetric, max asymmetry "//to_string(dsym, "es10.2"))
      if (allocated(error)) return
      doff = 0.0_dp
      if (is_pure) then
         ! A-B = diag(e_a - e_i): no J, no kernel, no exchange.
         do j = 1, nov
            do i = 1, nov
               if (i /= j) doff = max(doff, abs(amb(i, j)))
            end do
         end do
         call check(error, doff < 1.0e-10_dp, "A-B of a pure functional is not diagonal, max off-diagonal "// &
                    to_string(doff, "es10.2"))
         if (allocated(error)) return
      end if

      ! w^2 = eig of (A-B)^{1/2} (A+B) (A-B)^{1/2}.
      allocate (amh(nov, nov), work(nov, nov), sam(nov), w2(nov))
      amh = 0.5_dp*(amb + transpose(amb))
      call pic_syev(amh, sam, "V", "U", linfo)
      call check(error, linfo == 0 .and. sam(1) > 0.0_dp, "A-B is not positive definite")
      if (allocated(error)) return
      do j = 1, nov
         work(:, j) = amh(:, j)*sqrt(sam(j))
      end do
      amh = matmul(work, transpose(amh))
      work = matmul(amh, matmul(0.5_dp*(apb + transpose(apb)), amh))
      work = 0.5_dp*(work + transpose(work))
      call pic_syev(work, w2, "V", "U", linfo)
      call check(error, linfo == 0 .and. w2(1) > 0.0_dp, "(A-B)^{1/2} (A+B) (A-B)^{1/2} is not positive definite")
      if (allocated(error)) return
      dw = maxval(abs(sqrt(w2(1:nr)) - fx%td%omega))
      offmsg = ""
      if (is_pure) offmsg = ", A-B max off-diagonal "//to_string(doff, "es10.2")
      call global_logger%info("tddft rpa "//trim(molname_arg)//"/"//trim(funcname_arg)//": dense "// &
                              to_string(nov)//"x"//to_string(nov)//" problem, max|w(dense) - w(Davidson)| "// &
                              to_string(dw, "es10.2")//", max asymmetry "//to_string(dsym, "es10.2")//trim(offmsg))
      call check(error, dw < 1.0e-8_dp, "dense roots differ from the Davidson, max|dw| = "//to_string(dw, "es10.2"))
   end subroutine test_rpa_dense

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
      refname = ref_name()
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
      call check(error, fx%td%converged, "the solve did not converge")
      if (allocated(error)) return
      ! A second full solve: not on formamide, and for methane HF and B3LYP only.
      if (nroots_arg > 5 .or. skip_extra()) return
      o = fx%opts
      o%max_space = 3*nroots_arg
      o%max_iteration = 400
      call solve(o, td, err)
      call check(error, .not. err%has_error(), "collapsing solve failed: "//err%get_message())
      if (allocated(error)) return
      dom = maxval(abs(td%omega - fx%td%omega))
      call global_logger%info("tddft "//trim(molname_arg)//"/"//trim(funcname_arg)//": max_space "// &
                              to_string(3*nroots_arg)//" took "//to_string(td%iterations)//" iterations (vs "// &
                              to_string(fx%td%iterations)//"), "//to_string(td%ncollapse)//" collapses, max|d omega| "// &
                              to_string(dom, "es10.2"))
      call check(error, dom < 1.0e-8_dp, "collapsed solve differs, max|d omega| "//to_string(dom, "es10.2"))
      if (allocated(error)) return
      ! The RPA's subspace is rebuilt from p and m on a collapse, which the
      ! agreement above is only evidence for if a collapse happened.
      call check(error, td%ncollapse > 0, "the tiny max_space solve never collapsed")
   end subroutine test_collapse

   subroutine test_rpa_errors(error)
      !! The two RPA breakdown errors, provoked by corrupting the functional the
      !! solver sees (the SCF is untouched): a large exact-exchange weight makes
      !! A-B indefinite, a large negative kernel weight makes A+B indefinite. Each
      !! must be an error naming the matrix, with no result returned as converged.
      !! Both are detected in the Davidson subspace only.
      type(error_type), allocatable, intent(out) :: error
      type(trc_xc_functional_t) :: bad
      type(trc_tddft_options_t) :: o
      type(trc_tddft_result_t) :: td
      type(error_t) :: err

      call fixture_build(error)
      if (allocated(error)) return
      if (.not. fx%dft) then
         call global_logger%info("tddft rpa: error-path test skipped for Hartree-Fock (no functional to corrupt)")
         return
      end if
      o = fx%opts
      o%nroots = 1

      bad = fx%func
      bad%exx = 2.5_dp
      call solve(o, td, err, func=bad)
      call global_logger%info("tddft rpa "//trim(molname_arg)//"/"//trim(funcname_arg)//": c_x = 2.5 gives '"// &
                              trim(err%get_message())//"'")
      call check(error, err%has_error(), "c_x = 2.5 must be an error, not a silent return")
      if (allocated(error)) return
      call check(error, .not. td%converged, "c_x = 2.5 must not be reported as converged")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "A-B is not positive definite") > 0, &
                 "the error must name A-B, got: "//err%get_message())
      if (allocated(error)) return

      ! A pure functional's A-B is the diagonal e_a - e_i, so a kernel scaled to a
      ! large negative weight breaks A+B alone.
      if (fx%func%exx == 0.0_dp) then
         bad = fx%func
         bad%coef = 30.0_dp*bad%coef
         call solve(o, td, err, func=bad)
         call global_logger%info("tddft rpa "//trim(molname_arg)//"/"//trim(funcname_arg)// &
                                 ": a 30x kernel gives '"//trim(err%get_message())//"'")
         call check(error, err%has_error(), "a 30x kernel must be an error, not a silent return")
         if (allocated(error)) return
         call check(error, .not. td%converged, "a 30x kernel must not be reported as converged")
         if (allocated(error)) return
         call check(error, index(err%get_message(), "A+B is not positive definite") > 0, &
                    "the error must name A+B, got: "//err%get_message())
      end if
   end subroutine test_rpa_errors

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
      call check(error, len_trim(trc_tddft_refuse_reason(o)) == 0, "the RPA must not be refused")
      if (allocated(error)) return
      o%singlet = .false.
      call check_refused(error, trc_tddft_refuse_reason(o), "RPA triplet", "singlet")
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
      o%singlet = .false.
      call solve(o, td, err)
      call check(error, err%has_error() .and. .not. td%converged, "singlet=.false. must be an error")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "singlet") > 0, "refusal message: "//err%get_message())
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
      o%max_iteration = 1
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

end module check_tddft_rpa_suite

program check_tddft_rpa
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_tddft_rpa_suite, only: collect_check_tddft_rpa, xyz_arg, bas_json_arg, molname_arg, basname_arg, funcname_arg, &
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
   call run_testsuite(collect_check_tddft_rpa, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_tddft_rpa
