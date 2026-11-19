module check_tddft_basis_suite
   !!
   !! TDA and full RPA singlet excitation energies and oscillator strengths of a
   !! closed-shell reference through trc_tddft, across the shipped basis sets,
   !! against pyscf's tdscf TDA and TDHF / TDDFT (test/pyscf_refs/generate.py's
   !! gen_tddft_basis). check_tddft.F90 and check_tddft_rpa.F90 pin the solver
   !! on 6-31G* only; this suite is what says that the other basis families
   !! (minimal, split valence, diffuse, polarised, correlation consistent,
   !! def2, up to d functions) reach the solver intact: contraction
   !! and decontraction, general contractions, the Cartesian AO space, diffuse
   !! functions in the occupied-virtual space.
   !!
   !! One (molecule, basis, functional) per process, one SCF, one TDA and one
   !! RPA solve shared by the cases below (built by the first, released by the
   !! last):
   !!  1. TDA omega and f against the committed pyscf reference, every root;
   !!  2. RPA omega, f and ||Y||^2 against the reference (pyscf's (X, Y) are
   !!     normalised 2 (X.X - Y.Y) = 1, so ||Y||^2(terco) = 2 ||Y||^2(pyscf));
   !!  3. the residual of the solved equation of every root, TDA (A X - w X) and
   !!     RPA ([[A, B], [-B, -A]] (X, Y) - w (X, -Y)), from a FRESH operator
   !!     and fresh matvecs (not the solver's own subspace bookkeeping), below
   !!     the solver's conv_tol;
   !!  4. the normalisation: X.X = 1 (TDA), X.X - Y.Y = 1 (RPA), orthonormal in
   !!     that metric across roots, omega ascending;
   !!  5. w(RPA) <= w(TDA) root by root, the TDA roots from the solver itself;
   !!  6. negative controls, each of which must move omega past the tolerance:
   !!     the reference of a DIFFERENT basis of the same molecule and
   !!     functional (so the comparison in 1 is sensitive to the basis it
   !!     runs), and, for a Kohn-Sham functional, the same RPA solve with the
   !!     exact-exchange weight corrupted (110% of the functional's, or 0.1 for
   !!     a pure functional);
   !!  7. the lowest roots against the dense eigenvalues of A (TDA) and of
   !!     (A-B)^1/2 (A+B) (A-B)^1/2 (RPA), built from unit-vector matvecs when
   !!     there are at most MAX_DENSE pairs: no reference can say that a root is
   !!     MISSING from a converged, residual-clean answer, this can;
   !!     and again for nroots = 1 to 6 (Hartree-Fock);
   !!  8. a basis with a shell above MULT_LMAX is refused with an error that names
   !!     the limit, up front (any build), and after an SCF in an LMAX=3 build;
   !!  9. release.
   !!
   !! A case whose control argument is `none` has no reference and runs only 7.
   !!
   !! TDDFT is Cartesian only, so the basis is loaded with force_cart and the
   !! reference is pyscf's mol.cart = True. The grid is terco's own (the SCF's),
   !! written with terco's numbers to the probe file for pyscf.
   !!
   !! A basis whose highest angular momentum is above this build's LMAX cannot
   !! run an SCF: main prints "SKIP: <basis> needs l ..." and stops (ctest's
   !! anchored SKIP_REGULAR_EXPRESSION turns that line, and only that line, into
   !! "Skipped", never "Passed"; nothing else here prints "SKIP"). One above MULT_LMAX but within
   !! LMAX (an f basis in an LMAX=3 build: the dipole integrals of the
   !! oscillator strengths stop at d whatever TERCO_LMAX is) is an
   !! expected-refusal case: the SCF runs and both solves must return a clean
   !! error naming the limit, with no pyscf reference. A missing reference is a
   !! hard failure, never a skip.
   !!
   !! Tolerances are those of check_tddft.F90 and check_tddft_rpa.F90, not
   !! loosened for any basis: omega 1e-6 Eh and f 1e-5. The solver runs at
   !! conv_tol 1e-7; the reference SCFs at conv_tol 1e-12.
   !!
   !! Arguments: xyz basis.json molname basname funcname nroots probe.bin
   !! control_basname (the basis whose reference is the negative control).
   !!
   use pic_types, only: dp, int64, default_int
   use pic_lapack_interfaces, only: pic_syev
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_tddft, only: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_op_t, trc_tddft_solve_rks, &
                        trc_tddft_reference_t, trc_tddft_refuse_reason, tddft_window, tddft_space_t
   use trc_tables, only: LMAX
   use trc_mult_kernels, only: MULT_LMAX
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_vector
   implicit none
   private
   public :: collect_check_tddft_basis, collect_check_tddft_refusal, collect_check_tddft_unit
   public :: collect_check_tddft_dense, basis_beyond_lmax, refusal_mode
   public :: xyz_arg, bas_json_arg, molname_arg, basname_arg, funcname_arg, nroots_arg, probe_arg, control_arg

   real(dp), parameter :: ERI_THRESH = 1.0e-12_dp
   real(dp), parameter :: TOL_OMEGA = 1.0e-6_dp
   real(dp), parameter :: TOL_OSC = 1.0e-5_dp
   real(dp), parameter :: TOL_YY = 1.0e-6_dp
   real(dp), parameter :: CONV_TOL = 1.0e-7_dp
   integer, parameter :: MAX_DENSE = 400
   !! Largest number of occupied-virtual pairs for which the dense check runs.

   character(len=256) :: xyz_arg = 'water.xyz'
   character(len=256) :: bas_json_arg = '../basis_sets/6-31g.json'
   character(len=16) :: molname_arg = 'water'
   character(len=16) :: basname_arg = '6-31g'
   character(len=16) :: funcname_arg = 'hf'
   integer :: nroots_arg = 3
   character(len=256) :: probe_arg = 'tddft_basis_probe.bin'
   character(len=16) :: control_arg = 'sto-3g'
   logical :: refusal_mode = .false.
   !! Set by the program for a basis above MULT_LMAX: only the refusal is run.

   type :: ref_t
      real(dp) :: e_scf = 0.0_dp
      real(dp), allocatable :: om_tda(:), osc_tda(:), om_rpa(:), osc_rpa(:), yy(:)
   end type ref_t

   type :: fixture_t
      logical :: ready = .false.   ! set only once the SCF and both solves succeeded
      logical :: built = .false.   ! resources are allocated and test_release must free them
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
      type(trc_tddft_result_t) :: tda, rpa
      type(trc_tddft_reference_t) :: ref
      real(dp) :: t_scf = 0.0_dp, t_tda = 0.0_dp, t_rpa = 0.0_dp
      logical :: dense_ready = .false.
      real(dp), allocatable :: w_tda(:), w_rpa(:)
      !! The lowest dense eigenvalues of A and of the RPA problem (`dense_roots`).
   end type fixture_t

   type(fixture_t), save :: fx

contains

   subroutine collect_check_tddft_basis(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("lowest roots equal the dense eigenvalues of A and (A-B)^1/2 (A+B) (A-B)^1/2", test_dense_roots), &
                  new_unittest("the same for nroots = 1 to 6 (Hartree-Fock)", test_dense_nroots), &
                  new_unittest("TDA omega and f vs pyscf", test_tda_vs_pyscf), &
                  new_unittest("RPA omega, f and ||Y||^2 vs pyscf", test_rpa_vs_pyscf), &
                  new_unittest("residuals of fresh matvecs below conv_tol", test_residual), &
                  new_unittest("normalisation X.X = 1 (TDA), X.X - Y.Y = 1 (RPA)", test_normalisation), &
                  new_unittest("w(RPA) <= w(TDA), the TDA computed by the solver", test_below_tda), &
                  new_unittest("negative control: another basis's reference fails the comparison", test_control_basis), &
                  new_unittest("negative control: corrupted exchange weight fails the comparison", test_control_exx), &
                  new_unittest("release", test_release) &
                  ]
   end subroutine collect_check_tddft_basis

   subroutine collect_check_tddft_unit(testsuite)
      !! The checks that need no SCF and no reference, run once (the `unit` case),
      !! not in every case's process.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("a basis with f shells is refused up front, naming MULT_LMAX", test_f_refused_up_front), &
                  new_unittest("the window rule on a synthetic Ritz set: cap, count, open candidates", test_window_rule) &
                  ]
   end subroutine collect_check_tddft_unit

   subroutine collect_check_tddft_dense(testsuite)
      !! A case with no pyscf reference (control argument "none"): the solver
      !! against the dense eigenvalues only.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("lowest roots equal the dense eigenvalues of A and (A-B)^1/2 (A+B) (A-B)^1/2", test_dense_roots), &
                  new_unittest("the same for nroots = 1 to 6 (Hartree-Fock)", test_dense_nroots), &
                  new_unittest("release", test_release) &
                  ]
   end subroutine collect_check_tddft_dense

   subroutine collect_check_tddft_refusal(testsuite)
      !! The suite of a basis above MULT_LMAX (within LMAX): the refusal after a real SCF.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("TDA and RPA of an f basis return a clean error after the SCF", test_f_refused_after_scf), &
                  new_unittest("release", test_release) &
                  ]
   end subroutine collect_check_tddft_refusal

   function basis_beyond_lmax(highest) result(beyond)
      !! True when the basis has a shell above this build's LMAX; `highest`
      !! is its highest angular momentum. Read through the same loader the
      !! fixture uses, so the answer is the one the solver would meet.
      integer, intent(out) :: highest
      logical :: beyond
      type(trc_basis_t) :: bas
      type(error_t) :: err
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)

      call read_xyz(trim(xyz_arg), natm, zint, at_r)
      call trc_basis_from_json(trim(bas_json_arg), natm, zint, at_r, bas, err, force_cart=.true.)
      highest = -1
      if (.not. err%has_error()) highest = maxval(bas%sh_l)
      beyond = highest > LMAX
      call bas%release()
   end function basis_beyond_lmax

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
      call check(error, .not. fx%built, "fixture setup failed in an earlier test; nothing to compare")
      if (allocated(error)) return
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
      fx%built = .true.

      fx%nocc = nint(sum(real(zint, dp)))/2
      sopts%conv_energy = 1.0e-12_dp
      sopts%conv_density = 1.0e-10_dp
      sopts%conv_diis = 1.0e-10_dp
      sopts%max_iteration = 200
      if (fx%dft) sopts%functional = trim(funcname_arg)
      sopts%grid_level = 3
      call tm%start()
      call trc_scf_run(fx%bas, fx%nocc, fx%nocc, sopts, fx%scf)
      fx%t_scf = tm%get_elapsed_time()
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
      fx%opts%conv_tol = CONV_TOL
      if (refusal_mode) return

      fx%opts%tda = .true.
      call tm%start()
      call solve(fx%opts, fx%tda, err)
      fx%t_tda = tm%get_elapsed_time()
      call check(error, .not. err%has_error(), "TDA solve failed: "//err%get_message())
      if (allocated(error)) return
      call global_logger%info("tddft basis "//trim(tag())//": TDA "//to_string(fx%tda%iterations)//" iterations, "// &
                              to_string(fx%tda%nmatvec)//" matvecs, space "//to_string(fx%tda%space))

      fx%opts%tda = .false.
      call tm%start()
      call solve(fx%opts, fx%rpa, err)
      fx%t_rpa = tm%get_elapsed_time()
      call check(error, .not. err%has_error(), "RPA solve failed: "//err%get_message())
      if (allocated(error)) return
      call global_logger%info("tddft basis "//trim(tag())//": RPA "//to_string(fx%rpa%iterations)//" iterations, "// &
                              to_string(fx%rpa%nmatvec)//" matvecs, space "//to_string(fx%rpa%space))
      fx%ready = .true.
      call write_probe()
   end subroutine fixture_build

   subroutine write_probe()
      !! Terco's numbers and its grid, for generate.py's gen_tddft_basis.
      real(dp), allocatable :: yy(:)
      integer :: unit, k
      allocate (yy(nroots_arg))
      do k = 1, nroots_arg
         yy(k) = sum(fx%rpa%yamp(:, :, k)**2)
      end do
      open (newunit=unit, file=trim(probe_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(fx%bas%nao, int64), int(fx%natm, int64), int(fx%nocc, int64), int(nroots_arg, int64)
      write (unit) molname_arg, basname_arg, funcname_arg
      write (unit) fx%scf%energy
      write (unit) fx%tda%omega
      write (unit) fx%tda%osc
      write (unit) fx%rpa%omega
      write (unit) fx%rpa%osc
      write (unit) yy
      write (unit) int(size(fx%scf%grid_w), int64)
      write (unit) fx%scf%grid_r, fx%scf%grid_w
      close (unit)
   end subroutine write_probe

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

   function ref_name(basname) result(name)
      !! The committed pyscf reference of this molecule, functional and the
      !! given basis (the basis under test or the control's).
      character(len=*), intent(in) :: basname
      character(len=256) :: name
      name = 'tddft_basis_'//trim(molname_arg)//'_'//trim(basname)//'_'//trim(funcname_arg)//'.ref'
   end function ref_name

   subroutine read_ref(basname, r)
      character(len=*), intent(in) :: basname
      type(ref_t), intent(out) :: r
      integer :: u, nr
      nr = nroots_arg
      call ref_require(trim(ref_name(basname)), u)
      allocate (r%om_tda(nr), r%osc_tda(nr), r%om_rpa(nr), r%osc_rpa(nr), r%yy(nr))
      call ref_read_real(u, 'E_SCF', r%e_scf)
      call ref_read_vector(u, 'OMEGA_TDA', nr, r%om_tda)
      call ref_read_vector(u, 'OSC_TDA', nr, r%osc_tda)
      call ref_read_vector(u, 'OMEGA_RPA', nr, r%om_rpa)
      call ref_read_vector(u, 'OSC_RPA', nr, r%osc_rpa)
      call ref_read_vector(u, 'YY_RPA', nr, r%yy)
      call ref_close(u)
   end subroutine read_ref

   function tag() result(t)
      character(len=96) :: t
      t = trim(molname_arg)//"/"//trim(basname_arg)//"/"//trim(funcname_arg)
   end function tag

   subroutine check_refused_f(error, td, err, what)
      !! A refusal that names the limit and leaves no result behind.
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_result_t), intent(in) :: td
      type(error_t), intent(in) :: err
      character(len=*), intent(in) :: what
      call check(error, err%has_error() .and. .not. td%converged, what//": an f basis must be an error, not a result")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "MULT_LMAX") > 0 .and. index(err%get_message(), "l = 3") > 0, &
                 what//": the message must name the limit, got '"//err%get_message()//"'")
      if (allocated(error)) return
      call check(error, .not. allocated(td%omega) .and. .not. allocated(td%osc), what//": a refusal returned roots")
   end subroutine check_refused_f

   subroutine test_f_refused_up_front(error)
      !! Water in cc-pVTZ (f shells), no SCF: the refusal comes before any work, at any
      !! TERCO_LMAX, because the transition dipole's multipole integrals stop at
      !! MULT_LMAX = 2. Without the check the solve reached trc_multipoles and
      !! error-stopped the process, after the whole iteration. The operator, pair
      !! list and integral engine are never touched, so they stay unbuilt.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_linalg_t) :: la
      type(trc_tddft_reference_t) :: ref
      type(trc_tddft_options_t) :: o
      type(trc_tddft_result_t) :: td
      type(error_t) :: err
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)

      call read_xyz('water.xyz', natm, zint, at_r)
      call trc_basis_from_json('../basis_sets/cc-pvtz.json', natm, zint, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, maxval(bas%sh_l) == 3, "cc-pVTZ must have f shells")
      if (allocated(error)) return
      call check(error, len_trim(trc_tddft_refuse_reason(o, lmax=maxval(bas%sh_l))) > 0, "refuse_reason must refuse l = 3")
      if (allocated(error)) return
      call check(error, len_trim(trc_tddft_refuse_reason(o, lmax=MULT_LMAX)) == 0, "refuse_reason must accept l = MULT_LMAX")
      if (allocated(error)) return
      call trc_tddft_solve_rks(bas, pl, eri, la, ref, o, td)
      call check_refused_f(error, td, td%err, "TDA")
      if (allocated(error)) return
      o%tda = .false.
      call trc_tddft_solve_rks(bas, pl, eri, la, ref, o, td)
      call check_refused_f(error, td, td%err, "RPA")
      call bas%release()
   end subroutine test_f_refused_up_front

   subroutine test_window_rule(error)
      !! tddft_window on a synthetic Ritz set (nt = 3 always iterated, cap nkcap = 5,
      !! nr = 1, theta(1) = 0.30 Eh, theta(k) = 0.30 + 0.004 (k - 1)): a candidate
      !! is counted ONCE whatever its index (one at k = 25 takes one slot, so the
      !! collapse keeps nt + 1 vectors, not 25), the cap bounds the iterated ones,
      !! every candidate is counted in nopen, a pair with |r| >= theta/2 says too little
      !! and is not counted (|r| >= theta/2), a converged pair or one whose interval
      !! stays above theta(nr) is not a candidate.
      type(error_type), allocatable, intent(out) :: error
      type(tddft_space_t) :: sp
      real(dp) :: theta(30), rn(30)
      logical :: act(30)
      integer :: idx(30), nact, nopen, k

      sp%nt = 3
      sp%nkcap = 5
      do k = 1, 30
         theta(k) = 0.30_dp + 0.004_dp*real(k - 1, dp)
      end do
      rn = 1.0e-9_dp
      ! one candidate at a high index: |r| 0.15 < theta/2 and theta - r < 0.30
      rn(25) = 0.15_dp
      call tddft_window(sp, 1, theta, rn, 30, .true., 1.0e-3_dp, act, idx, nact, nopen)
      call check(error, nact == 4 .and. nopen == 1 .and. act(25) .and. idx(4) == 25 .and. all(idx(1:3) == [1, 2, 3]), &
                 "a single candidate at k = 25 must take one slot: nact "//to_string(nact)//", nopen "//to_string(nopen))
      if (allocated(error)) return
      ! eight candidates: the cap keeps the lowest two, all eight are counted
      rn(6) = 0.15_dp
      rn(8) = 0.15_dp
      rn(10) = 0.15_dp
      rn(12) = 0.15_dp
      rn(14) = 0.15_dp
      rn(16) = 0.15_dp
      rn(18) = 0.15_dp
      rn(20) = 0.15_dp
      call tddft_window(sp, 1, theta, rn, 30, .true., 1.0e-3_dp, act, idx, nact, nopen)
      call check(error, nact == 5 .and. nopen == 9 .and. all(idx(1:5) == [1, 2, 3, 6, 8]), &
                 "the cap must keep nkcap pairs, the lowest candidates first: nact "//to_string(nact)//", nopen "// &
                 to_string(nopen))
      if (allocated(error)) return
      ! no tracking (budget used up): nothing beyond nt is iterated, all are open
      call tddft_window(sp, 1, theta, rn, 30, .false., 1.0e-3_dp, act, idx, nact, nopen)
      call check(error, nact == 3 .and. nopen == 9, "without tracking only nt are kept and the candidates stay open: nact "// &
                 to_string(nact)//", nopen "//to_string(nopen))
      if (allocated(error)) return
      ! a pair with |r| >= theta/2 is not a candidate (iterated or counted); a converged pair
      ! and one whose interval stays above theta(nr) are not candidates either
      rn = 1.0e-9_dp
      rn(7) = 0.5_dp
      rn(9) = 0.2_dp
      rn(11) = 5.0e-4_dp
      rn(13) = 0.01_dp
      call tddft_window(sp, 1, theta, rn, 30, .true., 1.0e-3_dp, act, idx, nact, nopen)
      call check(error, nact == 3 .and. nopen == 0, &
                 "|r| >= theta/2, converged and clearly-above pairs must not be candidates: nact "//to_string(nact)// &
                 ", nopen "//to_string(nopen))
   end subroutine test_window_rule

   subroutine test_f_refused_after_scf(error)
      !! The SCF of an f basis (LMAX = 3 builds) converges, then both solves are
      !! refused: no error stop, a message that names the limit.
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_result_t) :: td
      type(error_t) :: err
      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%scf%converged .and. fx%scf%energy < 0.0_dp, "the f-basis SCF did not converge")
      if (allocated(error)) return
      fx%opts%tda = .true.
      call solve(fx%opts, td, err)
      call check_refused_f(error, td, err, "TDA")
      if (allocated(error)) return
      fx%opts%tda = .false.
      call solve(fx%opts, td, err)
      call check_refused_f(error, td, err, "RPA")
      if (allocated(error)) return
      call global_logger%info("tddft basis "//trim(tag())//": TDA and RPA refused after the SCF: "//trim(err%get_message()))
   end subroutine test_f_refused_after_scf

   subroutine test_tda_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(ref_t) :: r
      real(dp) :: dom, dosc

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%tda%converged, "TDA did not converge")
      if (allocated(error)) return
      call read_ref(basname_arg, r)
      call check(error, abs(fx%scf%energy - r%e_scf) < 1.0e-9_dp, &
                 "E vs pyscf, |dE| = "//to_string(abs(fx%scf%energy - r%e_scf), "es10.2"))
      if (allocated(error)) return
      dom = maxval(abs(fx%tda%omega - r%om_tda))
      dosc = maxval(abs(fx%tda%osc - r%osc_tda))
      call global_logger%info("tddft basis "//trim(tag())//": nao "//to_string(fx%bas%nao)//", SCF "// &
                              to_string(fx%t_scf, "f8.2")//" s (|dE| "// &
                              to_string(abs(fx%scf%energy - r%e_scf), "es8.1")//"); TDA max|d omega| "// &
                              to_string(dom, "es10.2")//" (tol "//to_string(TOL_OMEGA, "es8.1")//"), max|d f| "// &
                              to_string(dosc, "es10.2")//" (tol "//to_string(TOL_OSC, "es8.1")//"); "// &
                              to_string(fx%tda%iterations)//" iterations, solve "//to_string(fx%t_tda, "f8.2")// &
                              " s; omega(1) "//to_string(fx%tda%omega(1), "f12.8"))
      call check(error, dom < TOL_OMEGA, "TDA omega vs pyscf, max|diff| = "//to_string(dom, "es10.2"))
      if (allocated(error)) return
      call check(error, dosc < TOL_OSC, "TDA oscillator strength vs pyscf, max|diff| = "//to_string(dosc, "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(r%osc_tda) > 1.0e-4_dp, "no reference root has a measurable oscillator strength")
   end subroutine test_tda_vs_pyscf

   subroutine test_rpa_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(ref_t) :: r
      real(dp), allocatable :: yy(:)
      real(dp) :: dom, dosc, dyy
      integer :: k

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%rpa%converged, "RPA did not converge")
      if (allocated(error)) return
      call read_ref(basname_arg, r)
      allocate (yy(nroots_arg))
      do k = 1, nroots_arg
         yy(k) = sum(fx%rpa%yamp(:, :, k)**2)
      end do
      dom = maxval(abs(fx%rpa%omega - r%om_rpa))
      dosc = maxval(abs(fx%rpa%osc - r%osc_rpa))
      dyy = maxval(abs(yy - r%yy))
      call global_logger%info("tddft basis "//trim(tag())//": RPA max|d omega| "//to_string(dom, "es10.2")// &
                              " (tol "//to_string(TOL_OMEGA, "es8.1")//"), max|d f| "//to_string(dosc, "es10.2")// &
                              " (tol "//to_string(TOL_OSC, "es8.1")//"), max|d ||Y||^2| "//to_string(dyy, "es10.2")// &
                              " (tol "//to_string(TOL_YY, "es8.1")//"); "//to_string(fx%rpa%iterations)// &
                              " iterations, solve "//to_string(fx%t_rpa, "f8.2")//" s; omega(1) "// &
                              to_string(fx%rpa%omega(1), "f12.8"))
      call check(error, dom < TOL_OMEGA, "RPA omega vs pyscf, max|diff| = "//to_string(dom, "es10.2"))
      if (allocated(error)) return
      call check(error, dosc < TOL_OSC, "RPA oscillator strength vs pyscf, max|diff| = "//to_string(dosc, "es10.2"))
      if (allocated(error)) return
      call check(error, dyy < TOL_YY, "||Y||^2 vs pyscf, max|diff| = "//to_string(dyy, "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(yy) > 1.0e-6_dp, "Y is zero: not an RPA solution")
   end subroutine test_rpa_vs_pyscf

   subroutine test_residual(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_op_t) :: op
      type(error_t) :: err
      real(dp), allocatable :: ax(:, :, :), ax2(:, :, :, :), bx2(:, :, :, :), pp(:, :, :), mm(:, :, :)
      real(dp) :: rk, worst_tda, worst_rpa
      integer :: k, nr

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%tda%converged .and. fx%rpa%converged, "a solve did not converge")
      if (allocated(error)) return
      nr = nroots_arg
      call op%setup(fx%bas, fx%la, fx%ref, err, fx%xg)
      call check(error, .not. err%has_error(), "operator setup: "//err%get_message())
      if (allocated(error)) return
      allocate (ax(fx%nocc, op%nvir, nr), ax2(fx%nocc, op%nvir, nr, 2), bx2(fx%nocc, op%nvir, nr, 2), &
                pp(fx%nocc, op%nvir, nr), mm(fx%nocc, op%nvir, nr))
      call op%apply(fx%bas, fx%eri, fx%la, nr, fx%tda%xamp, ax, fx%xg)
      worst_tda = 0.0_dp
      do k = 1, nr
         rk = sqrt(sum((ax(:, :, k) - fx%tda%omega(k)*fx%tda%xamp(:, :, k))**2))
         worst_tda = max(worst_tda, rk)
      end do
      ! (A+B)(X+Y) and (A-B)(X-Y), each from its own fresh batch:
      ! A X + B Y = ((A+B)p + (A-B)m)/2, B X + A Y = ((A+B)p - (A-B)m)/2.
      pp = fx%rpa%xamp + fx%rpa%yamp
      mm = fx%rpa%xamp - fx%rpa%yamp
      call op%apply_rpa(fx%bas, fx%eri, fx%la, nr, pp, ax2, fx%xg)
      call op%apply_rpa(fx%bas, fx%eri, fx%la, nr, mm, bx2, fx%xg)
      call op%release()
      worst_rpa = 0.0_dp
      do k = 1, nr
         rk = sqrt(sum((0.5_dp*(ax2(:, :, k, 1) + bx2(:, :, k, 2)) - fx%rpa%omega(k)*fx%rpa%xamp(:, :, k))**2) &
                   + sum((0.5_dp*(ax2(:, :, k, 1) - bx2(:, :, k, 2)) + fx%rpa%omega(k)*fx%rpa%yamp(:, :, k))**2))
         worst_rpa = max(worst_rpa, rk)
      end do
      call global_logger%info("tddft basis "//trim(tag())//": fresh-matvec max residual TDA "// &
                              to_string(worst_tda, "es10.2")//", RPA "//to_string(worst_rpa, "es10.2")// &
                              " (conv_tol "//to_string(CONV_TOL, "es8.1")//")")
      call check(error, worst_tda < CONV_TOL, "TDA fresh residual "//to_string(worst_tda, "es10.2")//" >= conv_tol")
      if (allocated(error)) return
      call check(error, worst_rpa < CONV_TOL, "RPA fresh residual "//to_string(worst_rpa, "es10.2")//" >= conv_tol")
   end subroutine test_residual

   subroutine test_normalisation(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: dev_tda, dev_rpa
      integer :: j, k

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%tda%converged .and. fx%rpa%converged, "a solve did not converge")
      if (allocated(error)) return
      dev_tda = 0.0_dp
      dev_rpa = 0.0_dp
      do k = 1, nroots_arg
         do j = 1, nroots_arg
            dev_tda = max(dev_tda, abs(sum(fx%tda%xamp(:, :, j)*fx%tda%xamp(:, :, k)) - merge(1.0_dp, 0.0_dp, j == k)))
            ! The RPA's metric: X.X - Y.Y.
            dev_rpa = max(dev_rpa, abs(sum(fx%rpa%xamp(:, :, j)*fx%rpa%xamp(:, :, k)) &
                                       - sum(fx%rpa%yamp(:, :, j)*fx%rpa%yamp(:, :, k)) &
                                       - merge(1.0_dp, 0.0_dp, j == k)))
         end do
      end do
      call global_logger%info("tddft basis "//trim(tag())//": max|X^T X - 1| (TDA) "//to_string(dev_tda, "es10.2")// &
                              ", max|X^T X - Y^T Y - 1| (RPA) "//to_string(dev_rpa, "es10.2"))
      call check(error, dev_tda < 1.0e-10_dp, "TDA X not orthonormal, max dev "//to_string(dev_tda, "es10.2"))
      if (allocated(error)) return
      call check(error, dev_rpa < 1.0e-10_dp, "RPA (X, Y) not orthonormal in X.X - Y.Y, max dev "// &
                 to_string(dev_rpa, "es10.2"))
      if (allocated(error)) return
      do k = 2, nroots_arg
         call check(error, fx%tda%omega(k) >= fx%tda%omega(k - 1) .and. fx%rpa%omega(k) >= fx%rpa%omega(k - 1), &
                    "omega not ascending")
         if (allocated(error)) return
      end do
   end subroutine test_normalisation

   subroutine test_below_tda(error)
      !! The RPA's Y couples the TDA's excitations to de-excitations and can
      !! only lower the roots.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: dmax, dmin

      call fixture_build(error)
      if (allocated(error)) return
      call check(error, fx%tda%converged .and. fx%rpa%converged, "a solve did not converge")
      if (allocated(error)) return
      dmax = maxval(fx%rpa%omega - fx%tda%omega)
      dmin = minval(fx%tda%omega - fx%rpa%omega)
      call global_logger%info("tddft basis "//trim(tag())//": max (w_RPA - w_TDA) "//to_string(dmax, "es10.2")// &
                              " (must be <= 0), smallest TDA-RPA shift "//to_string(dmin, "es10.2"))
      call check(error, dmax <= 1.0e-9_dp, "an RPA root lies above its TDA root, by "//to_string(dmax, "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(fx%tda%omega - fx%rpa%omega) > 1.0e-7_dp, &
                 "no RPA root is below its TDA root: Y is not coupled in")
   end subroutine test_below_tda

   subroutine dense_roots(error)
      !! Fill fx%w_tda and fx%w_rpa (nroots_arg lowest, plus two) with the dense
      !! eigenvalues of A and of (A-B)^(1/2) (A+B) (A-B)^(1/2), built from nov
      !! unit-vector matvecs. Only for at most MAX_DENSE pairs; fx%dense_ready
      !! says whether it ran.
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_op_t) :: op
      type(error_t) :: err
      real(dp), allocatable :: unit_x(:, :, :), ax(:, :, :), prod(:, :, :, :), amat(:, :), apb(:, :), amb(:, :)
      real(dp), allocatable :: w(:), sam(:), amh(:, :), work(:, :), w2(:)
      integer :: nov, nocc, nvir, nl, j
      integer(default_int) :: linfo

      if (fx%dense_ready) return
      nocc = fx%nocc
      nvir = size(fx%scf%eps, 1) - nocc
      nov = nocc*nvir
      if (nov > MAX_DENSE) return
      nl = min(max(nroots_arg + 2, 6), nov)
      call op%setup(fx%bas, fx%la, fx%ref, err, fx%xg)
      call check(error, .not. err%has_error(), "operator setup: "//err%get_message())
      if (allocated(error)) return
      allocate (unit_x(nocc, nvir, nov), ax(nocc, nvir, nov), prod(nocc, nvir, nov, 2))
      unit_x = 0.0_dp
      do j = 1, nov
         unit_x(mod(j - 1, nocc) + 1, (j - 1)/nocc + 1, j) = 1.0_dp
      end do
      call op%apply(fx%bas, fx%eri, fx%la, nov, unit_x, ax, fx%xg)
      op%force_antisym = .true.
      call op%apply_rpa(fx%bas, fx%eri, fx%la, nov, unit_x, prod, fx%xg)
      call op%release()

      allocate (amat(nov, nov), w(nov))
      amat = reshape(ax, [nov, nov])
      amat = 0.5_dp*(amat + transpose(amat))
      call pic_syev(amat, w, "N", "U", linfo)
      call check(error, linfo == 0, "dense TDA eigenproblem failed")
      if (allocated(error)) return
      allocate (apb(nov, nov), amb(nov, nov), amh(nov, nov), work(nov, nov), sam(nov), w2(nov))
      apb = reshape(prod(:, :, :, 1), [nov, nov])
      amb = reshape(prod(:, :, :, 2), [nov, nov])
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
      call pic_syev(work, w2, "N", "U", linfo)
      call check(error, linfo == 0 .and. w2(1) > 0.0_dp, "(A-B)^1/2 (A+B) (A-B)^1/2 is not positive definite")
      if (allocated(error)) return
      fx%w_tda = w(1:nl)
      fx%w_rpa = sqrt(w2(1:nl))
      fx%dense_ready = .true.
   end subroutine dense_roots

   subroutine test_dense_roots(error)
      !! The solver's lowest roots against the dense problem: eigenvalues of A (TDA)
      !! and of (A-B)^(1/2) (A+B) (A-B)^(1/2) (RPA). This is the check no pyscf
      !! reference can give for a missed root: it fails when a root below the
      !! highest returned one is absent (the formaldehyde/6-31+G/HF case, where
      !! the 0.348172 Eh root was never returned). Cost: nov matvecs, for at most
      !! MAX_DENSE pairs.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: dtda, drpa
      character(len=512) :: roots

      call fixture_build(error)
      if (allocated(error)) return
      call dense_roots(error)
      if (allocated(error)) return
      if (.not. fx%dense_ready) then
         call global_logger%info("tddft basis "//trim(tag())//": dense root check skipped (more than "// &
                                 to_string(MAX_DENSE)//" pairs)")
         return
      end if
      call check(error, .not. (fx%tda%window_open .or. fx%rpa%window_open), "the missed-root check left candidates open ("// &
                 to_string(fx%tda%nwindow_open)//" TDA, "//to_string(fx%rpa%nwindow_open)//" RPA) on a case with no missing root")
      if (allocated(error)) return
      dtda = maxval(abs(fx%w_tda(1:nroots_arg) - fx%tda%omega))
      drpa = maxval(abs(fx%w_rpa(1:nroots_arg) - fx%rpa%omega))
      write (roots, "(*(f10.6))") fx%w_tda
      call global_logger%info("tddft basis "//trim(tag())//": dense A lowest roots "//trim(roots))
      write (roots, "(*(f10.6))") fx%tda%omega
      call global_logger%info("tddft basis "//trim(tag())//": TDA solver roots     "//trim(roots))
      call global_logger%info("tddft basis "//trim(tag())//": dense vs solver, max|dw| TDA "//to_string(dtda, "es10.2")// &
                              ", RPA "//to_string(drpa, "es10.2"))
      call check(error, dtda < 100.0_dp*CONV_TOL, "TDA: a root is missing or wrong, max|w(dense) - w(Davidson)| = "// &
                 to_string(dtda, "es10.2"))
      if (allocated(error)) return
      call check(error, drpa < 100.0_dp*CONV_TOL, "RPA: a root is missing or wrong, max|w(dense) - w(Davidson)| = "// &
                 to_string(drpa, "es10.2"))
   end subroutine test_dense_roots

   subroutine test_dense_nroots(error)
      !! The same comparison for the other numbers of roots of 1 to 6 (the case's own
      !! nroots_arg is the test above), Hartree-Fock only: a Kohn-Sham matvec costs
      !! far more. The lowest roots of formaldehyde/
      !! 6-31+G/HF were missed for nroots = 3, 4 and 5 alike, and what a solver
      !! returns must not depend on how many roots were asked for. A solve that
      !! returns k roots must equal the lowest k dense eigenvalues, for TDA and RPA.
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_options_t) :: o
      type(trc_tddft_result_t) :: td
      type(error_t) :: err
      real(dp) :: dw
      integer :: nr, pass, nl

      call fixture_build(error)
      if (allocated(error)) return
      if (fx%dft) return
      call dense_roots(error)
      if (allocated(error)) return
      if (.not. fx%dense_ready) return
      nl = size(fx%w_tda)
      do nr = 1, 6
         if (nr == nroots_arg .or. nr > nl) cycle
         do pass = 1, 2
            o = fx%opts
            o%nroots = nr
            o%tda = pass == 1
            call solve(o, td, err)
            call check(error, .not. err%has_error(), "solve failed: "//err%get_message())
            if (allocated(error)) return
            if (pass == 1) then
               dw = maxval(abs(fx%w_tda(1:nr) - td%omega))
            else
               dw = maxval(abs(fx%w_rpa(1:nr) - td%omega))
            end if
            call global_logger%info("tddft basis "//trim(tag())//": nroots "//to_string(nr)//" "// &
                                    merge("TDA", "RPA", pass == 1)//" "//to_string(td%nmatvec)//" matvecs, max|dw| "// &
                                    to_string(dw, "es10.2"))
            call check(error, dw < 100.0_dp*CONV_TOL, merge("TDA", "RPA", pass == 1)//" with nroots = "//to_string(nr)// &
                       ": a root is missing, max|w(dense) - w(Davidson)| = "//to_string(dw, "es10.2"))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_dense_nroots

   subroutine test_control_basis(error)
      !! The solved roots against another basis's reference must differ by far
      !! more than the tolerance, or the comparison could not tell the bases
      !! apart.
      type(error_type), allocatable, intent(out) :: error
      type(ref_t) :: r
      real(dp) :: dtda, drpa

      call fixture_build(error)
      if (allocated(error)) return
      call read_ref(control_arg, r)
      dtda = maxval(abs(fx%tda%omega - r%om_tda))
      drpa = maxval(abs(fx%rpa%omega - r%om_rpa))
      call global_logger%info("tddft basis "//trim(tag())//": negative control, the "//trim(control_arg)// &
                              " reference differs by TDA "//to_string(dtda, "es10.2")//", RPA "// &
                              to_string(drpa, "es10.2")//" (tolerance "//to_string(TOL_OMEGA, "es8.1")//")")
      call check(error, dtda > 100.0_dp*TOL_OMEGA, "another basis's TDA reference must fail the comparison, "// &
                 "differs only by "//to_string(dtda, "es10.2"))
      if (allocated(error)) return
      call check(error, drpa > 100.0_dp*TOL_OMEGA, "another basis's RPA reference must fail the comparison, "// &
                 "differs only by "//to_string(drpa, "es10.2"))
   end subroutine test_control_basis

   subroutine test_control_exx(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_xc_functional_t) :: bad
      type(trc_tddft_options_t) :: o1
      type(trc_tddft_result_t) :: td
      type(error_t) :: err
      type(ref_t) :: r
      real(dp) :: dom

      call fixture_build(error)
      if (allocated(error)) return
      if (.not. fx%dft) then
         call global_logger%info("tddft basis: exchange-weight control skipped for Hartree-Fock (no weight)")
         return
      end if
      bad = fx%func
      if (fx%func%exx > 0.0_dp) then
         bad%exx = 1.1_dp*fx%func%exx
      else
         bad%exx = 0.1_dp
      end if
      o1 = fx%opts
      o1%nroots = 1
      o1%tda = .false.
      call solve(o1, td, err, func=bad)
      call check(error, .not. err%has_error() .and. td%converged, "corrupted solve failed: "//err%get_message())
      if (allocated(error)) return
      call read_ref(basname_arg, r)
      dom = abs(td%omega(1) - r%om_rpa(1))
      call global_logger%info("tddft basis "//trim(tag())//": negative control (c_x "//to_string(fx%func%exx, "f6.3")// &
                              " -> "//to_string(bad%exx, "f6.3")//") moves the RPA omega by "// &
                              to_string(dom, "es10.2")//" (tolerance "//to_string(TOL_OMEGA, "es8.1")//")")
      call check(error, dom > 100.0_dp*TOL_OMEGA, "a corrupted exchange weight must fail the omega comparison, "// &
                 "moved only "//to_string(dom, "es10.2"))
   end subroutine test_control_exx

   subroutine test_release(error)
      type(error_type), allocatable, intent(out) :: error
      if (fx%built) then
         call fx%eri%release()
         call fx%la%release()
         if (fx%dft .and. fx%scf%converged) call fx%xg%release()
         call fx%pl%release()
         call fx%bas%release()
         fx%built = .false.
         fx%ready = .false.
      end if
   end subroutine test_release

end module check_tddft_basis_suite

program check_tddft_basis
   use, intrinsic :: iso_fortran_env, only: error_unit
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use testdrive, only: run_testsuite
   use trc_tables, only: LMAX
   use trc_mult_kernels, only: MULT_LMAX
   use check_tddft_basis_suite, only: collect_check_tddft_basis, collect_check_tddft_refusal, collect_check_tddft_unit, &
                                      collect_check_tddft_dense, basis_beyond_lmax, refusal_mode, xyz_arg, bas_json_arg, &
                                      molname_arg, basname_arg, funcname_arg, &
                                      nroots_arg, probe_arg, control_arg
   implicit none
   integer :: stat, highest
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
   if (command_argument_count() >= 8) call get_command_argument(8, control_arg)

   ! The `unit` case (no SCF, no reference): the checks that run once.
   stat = 0
   if (trim(basname_arg) == 'unit') then
      call run_testsuite(collect_check_tddft_unit, error_unit, stat, parallel=.false.)
      if (stat > 0) then
         write (error_unit, "(i0, a)") stat, " test(s) failed"
         error stop 1
      end if
      stop 0
   end if

   ! A basis above this build's LMAX has no SCF to start from: nothing to test
   ! yet, not an error. This is the only line the program prints with "SKIP", and
   ! test/CMakeLists.txt's SKIP_REGULAR_EXPRESSION matches exactly it.
   if (basis_beyond_lmax(highest)) then
      call global_logger%info("SKIP: "//trim(basname_arg)//" needs l = "//to_string(highest)// &
                              " and this build has LMAX = "//to_string(LMAX)//".")
      stop 0
   end if

   stat = 0
   ! Above MULT_LMAX the transition dipole cannot be built: the solver must refuse.
   refusal_mode = highest > MULT_LMAX
   if (refusal_mode) then
      call run_testsuite(collect_check_tddft_refusal, error_unit, stat, parallel=.false.)
   else if (trim(control_arg) == 'none') then
      call run_testsuite(collect_check_tddft_dense, error_unit, stat, parallel=.false.)
   else
      call run_testsuite(collect_check_tddft_basis, error_unit, stat, parallel=.false.)
   end if
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_tddft_basis
