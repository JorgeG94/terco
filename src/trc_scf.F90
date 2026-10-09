!
! The SCF: a basis and an electron count in, a converged energy and
! density out, with every matrix resident on the device.
!
! terco's contract used to stop at the matrices -- give me a density, get
! a Fock matrix -- and the iteration was the caller's. This module is the
! other half of that contract, so a host code can hand over a basis and a
! guess and receive a result: Hartree-Fock, or Kohn-Sham with any of the
! functionals in trc_xc_functional, restricted or unrestricted.
!
! WHAT CROSSES THE BUS
! --------------------
! Per iteration: two energies, a norm, and a DIIS matrix of at most ten by
! ten. The density, the Fock matrices, the orbitals, the DIIS history and
! the orthogonaliser are put on the device once and stay there; the Fock
! build is the resident driver, the XC integrator finds its density already
! present, and the products, the commutator and the eigenproblem are
! trc_linalg, which is cuBLAS and cuSOLVER under the OpenACC build. That
! is the design metalquicha's device SCF settled on after measuring host
! LAPACK plus copies at three times the cost of the GPU Fock build, and it
! is followed here rather than rediscovered.
!
! Under a host build the same source runs on host arrays with LAPACK, and
! the numbers agree to rounding; that is what makes the host suite a test
! of this module.
!
! THE ITERATION
! -------------
! A few damped steps from the guess, then Pulay's DIIS on the commutator
! FDS - SDF in the orthogonal basis, per spin with one set of coefficients
! for both. Converged when the energy change and the RMS density change
! are both below their thresholds. The DIIS system is at most eleven by
! eleven and is solved here by pivoted elimination, as metalquicha's
! mqc_diis does, with the error block scaled to order one against the
! constraint border first: near convergence the overlaps are 1e-16 and an
! unscaled elimination makes the weights noise, which shows not as a
! blow-up but as a density that drifts from the fixed point while the
! energy sits flat.
!
! UNRESTRICTED
! ------------
! F_s = H + J[D_a + D_b] - a_x K[D_s] + V_xc,s. The Fock driver returns
! j J[D] - k K[D]/2 for a density it takes to be a total, so J comes from
! one call on D_a + D_b with k = 0, and -a_x K[D_s] from one call on 2 D_s
! with j = 0 and k = a_x. Three builds per iteration where the restricted
! case needs one; a batched build over the three densities is the obvious
! improvement and is not done here.
!
module trc_scf_driver
   use, intrinsic :: iso_fortran_env, only: int64
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t
   use trc_dft_grid, only: dft_grid_t, build_dft_grid, build_dft_grid_block
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_xc, only: trc_xc_rks, trc_xc_uks, XC_RHO_TOL
   use trc_linalg, only: trc_linalg_t
   use pic_types, only: default_int
   use pic_mpi_lib, only: comm_t, bcast, allreduce, MPI_SUM
   implicit none
   private

   public :: trc_scf_options_t, trc_scf_result_t, trc_scf_run

   type :: trc_scf_options_t
      character(len=16) :: functional = ""   !! empty is Hartree-Fock
      !! Starting orbitals when no density is given: "gwh" (Wolfsberg-Helmholz
      !! over the overlap and the core Hamiltonian, the default) or "core".
      !! The core guess is hopeless past a few atoms -- cholesterol starts
      !! 660 hartree high and never recovers -- and GWH costs nothing.
      character(len=8) :: guess = "gwh"
      logical :: unrestricted = .false.      !! forced; an open shell is unrestricted regardless
      real(dp) :: conv_energy = 1.0e-10_dp
      real(dp) :: conv_density = 1.0e-7_dp   !! RMS change per element, printed; no longer gates
      !! Frobenius norm of the orthogonalised commutator X^T (FDS - SDF) X,
      !! the DIIS error, below which the SCF is converged once the energy
      !! change is too. The density change was the gate before, and on
      !! cholesterol/PBE it sat at the noise floor of the quadrature and the
      !! integral screen: energy converged to 1e-11 by iteration 40, then 37
      !! more iterations waiting for an RMS change under 1e-7 that wandered
      !! between 1e-7 and 1e-6. PySCF stops on the gradient norm at
      !! sqrt(conv_tol) for the same reason.
      real(dp) :: conv_diis = 1.0e-5_dp
      !! FRACTIONAL OCCUPATIONS, for the free-atom runs behind the SAD guess.
      !! Each iteration the orbitals of a spin channel are filled by aufbau on
      !! their own eigenvalues, and the degenerate set at the frontier shares
      !! what is left evenly, so the atom comes out spherical. Restricted
      !! (`unrestricted` off) fills `nelec_frac` electrons two per orbital;
      !! unrestricted fills `nalpha` and `nbeta` one per orbital in their own
      !! channels, a proper spin-polarised atom (carbon: two alpha thirds in
      !! each p, no beta p). `nalpha`/`nbeta` are ignored for the count only
      !! in the restricted case.
      logical :: frac_occ = .false.
      real(dp) :: nelec_frac = 0.0_dp
      integer :: max_iter = 100
      !> Shell-pair batch resolution: buckets per decade of Schwarz bound
      !> when the pairs are binned. 1 is decades, which is what this always
      !> did; higher values make the pre-launch screen finer, at more and
      !> smaller batches. Settable through trc_set_batching.
      integer :: batch_res = 1
      !> Rebuild the two-electron matrix from the DENSITY DIFFERENCE rather
      !> than the density, accumulating onto the last one. G is linear in D,
      !> so this is exact arithmetic; what it buys is screening, because the
      !> difference shrinks as the SCF converges while the density does not.
      !> On a 123-atom silica slice the later iterations are where the time
      !> goes, and they are exactly the ones the difference screens hardest.
      !> A full build is taken every `incr_reset` iterations as a backstop
      !> against long-run drift, though at the measured `incr_thresh_scale`
      !> and `incr_start` below a single increment's own error is already
      !> inside `conv_energy` and the backstop rarely has to do any work.
      !>
      !> WHY THIS WAS OFF, AND WHAT CHANGED. Handing the difference density
      !> to the SAME absolute screening threshold a full build uses was the
      !> bug, not the idea: the dropped quartets are bounded in absolute
      !> terms, which is right when the density hanging over them is O(1),
      !> and is not a bound that shrinks with a difference that is heading
      !> to zero. On w172/6-31G (450 atoms, 1950 AOs) that showed up as
      !> RMS(D) reaching 1e-10 while the energy kept bouncing by 1e-4 to
      !> 1e-6 Eh every iteration with every full rebuild agreeing with every
      !> other to 1e-9 -- consistent, so the dropped set IS the same one
      !> every full build, but each DIFFERENCE build drops a different set,
      !> because ddif's spatial pattern changes iteration to iteration the
      !> way D's does not, and those errors do not cancel. `incr_thresh_scale`
      !> and `incr_start` below are the fix, measured on the same case:
      !> 200 iterations, not converged, becomes 20, converged to the same
      !> energy as incremental-off (-11395.98363707 against -11395.98363707,
      !> both well inside 1e-8), in 23.3 s against the 24.4 s reference.
      !> Set TRC_INCREMENTAL=1 to turn it on without a rebuild.
      logical :: incremental = .false.
      integer :: incr_reset = 20
      !! Multiplies `eri_thresh` for an incremental Fock build only, so its
      !! screen is tighter than the full build's by this fixed factor
      !! rather than tracking how small the increment itself has become --
      !! that was tried first (scale by ||ddif||/||dref||) and is backwards:
      !! it loosens towards "keep everything" as ddif shrinks, which buys
      !! back the accuracy at the cost of the speed an increment build
      !! exists for (measured on w172: 25.6 s, slower than incremental off).
      !! A FIXED tighter floor keeps the per-quartet error small enough,
      !! every time, that an unlucky run of iterations whose dropped sets
      !! do not cancel still cannot add up past `conv_energy`. 1e-4,
      !! measured on w172/6-31G (450 atoms, 1950 AOs), reproduces the
      !! incremental-off energy to 2e-10 (comfortably inside `conv_energy`
      !! at 1e-8) while still costing less than a full build once the
      !! density difference is small enough for `incr_start` to have let
      !! incrementing start.
      real(dp) :: incr_thresh_scale = 1.0e-4_dp
      !! Do not increment until the PREVIOUS iteration's commutator norm
      !! is already below this. Early on, before CDIIS has settled, D
      !! swings by a large fraction of itself iteration to iteration --
      !! on w172 the energy jumps 1267 Eh between iterations 3 and 4 -- so
      !! ddif is no smaller than D itself, and `incr_thresh_scale`'s fixed
      !! floor then makes an increment survey MORE quartets than a full
      !! build would, not fewer (the floor is tight against a bound, dd,
      !! that is not small yet). A handful of full builds while the
      !! commutator is still large costs less than that: measured on
      !! w172, gating the fixed-scale increment on `errmax < incr_start`
      !! rather than starting it at iteration 2 turns a 34 s run (slower
      !! than incremental off) into a 23.3 s one (faster than the 24.4 s
      !! reference), at the same energy.
      real(dp) :: incr_start = 1.0e-2_dp
      integer :: ndiis = 10
      !! The damped start, -1 meaning "by spin": a closed shell goes to DIIS
      !! from the first iteration -- damping only delayed it until the
      !! commutator fell under `diis_start`, which took cholesterol from 29
      !! iterations to 51 and water from 10 to 19 -- while an open shell
      !! keeps four damped steps and a 0.5 start, because UHF on the OH
      !! radical from GWH oscillates for a hundred iterations without them.
      !! Any value >= 0 here, or a finite `diis_start`, is taken as given.
      integer :: ndamp = -1                  !! damped steps before DIIS, at least
      !! DIIS also waits until the largest commutator element is below this;
      !! extrapolating from a wild guess wanders, damping does not.
      real(dp) :: diis_start = -1.0_dp     !! commutator norm below which DIIS may start; -1 by spin
      real(dp) :: damp = 0.3_dp              !! weight of the new density while damping
      !! Level shift on the virtual space while DIIS is off, in hartree.
      !! From a Wolfsberg-Helmholz guess cholesterol sits 20 hartree above
      !! its answer and a damped Roothaan step alone throws it 300 hartree
      !! away; raising the virtuals keeps the occupied space from collapsing.
      !! ADIIS (Hu & Yang, JCP 132, 054109) while damping has ended but the
      !! commutator is still too large for CDIIS (see `diis_start`):
      !! minimise an energy functional over the convex hull of the
      !! Fock/density history instead of taking the level-shifted Roothaan
      !! step. Off by default -- measured against the existing suite,
      !! where `ndamp`/`level_shift` already reach the CDIIS window in the
      !! same iteration count, so turning it on there only pays for the
      !! history bookkeeping. Where it earns its keep is the hard cases
      !! (stretched bonds, small gaps, charged open shells) that plain
      !! damping oscillates on; see the robustness test set.
      logical :: adiis = .false.
      real(dp) :: level_shift = 1.0_dp
      !! Canonical orthogonalisation: an eigenvalue of S below this is
      !! LINEAR DEPENDENCE, and that direction is projected out of X rather
      !! than scaled by S^-1/4, which would divide by (nearly) zero. A
      !! compact double-zeta basis never gets near this; a diffuse one
      !! does -- def2-TZVPD and aug-cc-pVDZ put several S eigenvalues below
      !! 1e-8 on anything but a single atom, because two diffuse shells on
      !! neighbouring centres are nearly the same function. 0 restores the
      !! old behaviour: `error stop` the first time S is not numerically
      !! positive definite, which is what every existing test still gets
      !! since none of them is diffuse enough to reach this threshold.
      real(dp) :: lindep_thresh = 1.0e-8_dp
      real(dp) :: eri_thresh = 1.0e-12_dp
      integer :: grid_level = 3
      integer :: grid_max_pts = 512
      real(dp) :: rho_tol = XC_RHO_TOL
      logical :: verbose = .false.
   end type trc_scf_options_t

   type :: trc_scf_result_t
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: nspin = 1
      real(dp) :: energy = 0.0_dp      !! total
      real(dp) :: e_nuc = 0.0_dp
      real(dp) :: e_one = 0.0_dp       !! Tr(D H)
      real(dp) :: e_two = 0.0_dp       !! Coulomb and exact exchange
      real(dp) :: e_xc = 0.0_dp
      real(dp) :: nelec_grid = 0.0_dp  !! what the grid integrated the density to
      !! Functions dropped by canonical orthogonalisation: S eigenvalues
      !! below `opts%lindep_thresh`. 0 for every compact basis; see
      !! `lindep_thresh` for why a diffuse one is different.
      integer :: n_lindep = 0
      real(dp), allocatable :: dmat(:, :, :)   !! (nao, nao, nspin); restricted holds the total
      real(dp), allocatable :: cmo(:, :, :)    !! (nao, nao, nspin) orbitals
      real(dp), allocatable :: eps(:, :)       !! (nao, nspin) orbital energies
      real(dp), allocatable :: grid_r(:, :)    !! (3, npts) the grid used, for a same-grid comparison
      real(dp), allocatable :: grid_w(:)
      character(len=256) :: message = ""
   end type trc_scf_result_t

contains

   !> The density at a chosen iteration, written where a benchmark can read
   !> it back. Off unless TRC_DUMP_DENSITY names a file. Stream format:
   !> an 8-byte nao, then the matrix, which is what the test readers expect.
   subroutine dump_density(nao, d)
      integer, intent(in) :: nao
      real(dp), intent(in) :: d(nao, nao)
      character(len=512) :: path
      integer :: u, ios
      integer(int64) :: prev
      call get_environment_variable("TRC_DUMP_DENSITY", path)
      if (len_trim(path) == 0) return
      ! The atomic SCFs of a SAD guess come through here too, and would
      ! overwrite the molecular density with a 5x5 one. Keep the largest.
      open (newunit=u, file=trim(path), access='stream', form='unformatted', &
            status='old', iostat=ios)
      if (ios == 0) then
         read (u, iostat=ios) prev
         close (u)
         if (ios == 0) then
            if (prev >= int(nao, int64)) return
         end if
      end if
      open (newunit=u, file=trim(path), access='stream', form='unformatted', &
            status='replace', iostat=ios)
      if (ios /= 0) return
      write (u) int(nao, int64)
      write (u) d
      close (u)
      print '(a,a,a,i0,a)', "   [density dumped to ", trim(path), ", nao ", nao, "]"
   end subroutine dump_density

   !> Experiment hook for the incremental build, so its effect can be
   !> measured through a driver that does not expose the option yet.
   logical function incr_env()
      character(len=8) :: v
      call get_environment_variable("TRC_INCREMENTAL", v)
      incr_env = (trim(v) == "1")
   end function incr_env

   !> Full rebuild interval, and an integral threshold override. A
   !> difference build carries the same absolute screening error as a full
   !> one, so the error accumulates across increments and DIIS stalls
   !> before the energy converges; the cure is a tighter integral
   !> threshold on the difference builds, and resetting often enough that
   !> what does accumulate stays below it. Both are measured here rather
   !> than assumed.
   integer function incr_reset_env(dflt)
      integer, intent(in) :: dflt
      character(len=8) :: v
      integer :: ios, k
      call get_environment_variable("TRC_INCR_RESET", v)
      incr_reset_env = dflt
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) k
         if (ios == 0 .and. k > 0) incr_reset_env = k
      end if
   end function incr_reset_env

   !> Experiment hook, like `incr_reset_env`: TRC_INCR_THRESH_SCALE sweeps
   !> `incr_thresh_scale` without a rebuild.
   real(dp) function incr_thresh_scale_env(dflt)
      real(dp), intent(in) :: dflt
      character(len=16) :: v
      integer :: ios
      real(dp) :: t
      call get_environment_variable("TRC_INCR_THRESH_SCALE", v)
      incr_thresh_scale_env = dflt
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) t
         if (ios == 0 .and. t > 0.0_dp) incr_thresh_scale_env = t
      end if
   end function incr_thresh_scale_env

   !> Experiment hook, like `incr_reset_env`: TRC_INCR_START sweeps
   !> `incr_start` without a rebuild.
   real(dp) function incr_start_env(dflt)
      real(dp), intent(in) :: dflt
      character(len=16) :: v
      integer :: ios
      real(dp) :: t
      call get_environment_variable("TRC_INCR_START", v)
      incr_start_env = dflt
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) t
         if (ios == 0 .and. t > 0.0_dp) incr_start_env = t
      end if
   end function incr_start_env

   real(dp) function eri_thresh_env(dflt)
      real(dp), intent(in) :: dflt
      character(len=16) :: v
      integer :: ios
      real(dp) :: t
      call get_environment_variable("TRC_ERI_THRESH", v)
      eri_thresh_env = dflt
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) t
         if (ios == 0 .and. t > 0.0_dp) eri_thresh_env = t
      end if
   end function eri_thresh_env

   !
   ! Run the SCF. `nalpha` and `nbeta` are the electron counts; the case is
   ! restricted when they are equal and `unrestricted` is not forced. The
   ! basis must already be on the device. `dguess` is (nao, nao, nspin) in
   ! the result's convention; without it the guess is the core Hamiltonian.
   !
   subroutine trc_scf_run(b, nalpha, nbeta, opts, res, dguess, comm, hguess)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nalpha, nbeta
      type(trc_scf_options_t), intent(in) :: opts
      type(trc_scf_result_t), intent(out) :: res
      real(dp), intent(in), optional :: dguess(:, :, :)
      !! A guess HAMILTONIAN to occupy instead of hcore or the GWH matrix --
      !! how SAP gets in (see trc_sap), which cannot come through `dguess`
      !! because what it produces is a one-electron operator and not a
      !! density. Ignored if `dguess` is also given: a density the caller
      !! already has is more specific than a matrix to make one from.
      real(dp), intent(in), optional :: hguess(:, :)
      !! Ranks. Every rank runs the same iteration on the same matrices;
      !! the quartets are split across them and the Fock matrix summed
      !! back, so the density stays identical everywhere. Bind devices
      !! first (trc_bind_device). The XC batches are split across ranks too.
      type(comm_t), intent(in), optional :: comm

      integer :: nao, nspin, s, it, nocc(2), ndiis_used, i, j, k, n2, m
      logical :: ok
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :), hcore(:, :), x(:, :)
      real(dp), allocatable :: fock(:, :, :), gmat(:, :, :), vxc(:, :, :), dtot(:, :), jmat(:, :)
      real(dp), allocatable :: dref(:, :, :), gacc(:, :, :), ddif(:, :)
      logical :: incr_now
      real(dp), allocatable :: dold(:, :, :), errv(:, :, :), fstore(:, :, :, :), estore(:, :, :, :)
      real(dp), allocatable :: w1(:, :), w2(:, :), w3(:, :), bmat(:, :), rhs(:)
      !! ADIIS history (fock, density) and its own slot count, kept apart
      !! from the CDIIS `fstore`/`estore` above: ADIIS runs exactly in the
      !! window CDIIS has not opened yet, so the two never need the same
      !! entries, and keeping them separate means `opts%adiis = .false.`
      !! leaves the CDIIS path untouched, to the bit.
      real(dp), allocatable :: fstore_a(:, :, :, :), dstore_a(:, :, :, :), cad(:)
      integer :: nadiis_used
      logical :: adiis_active
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(trc_xc_functional_t) :: func
      type(trc_linalg_t) :: la
      type(error_t) :: err
      logical :: dft, talk, diis_on
      real(dp) :: exx, etot, eold, drms, e1, e2, exc, ngrid, occ, errmax
      real(dp), allocatable :: ofrac(:)
      integer :: nfrac, ndamp
      real(dp) :: diis_start
      real(dp) :: tw0, tw1, t_setup, t_grid, t_fock, t_xc, t_rest, tx_pts, tx_prs, tp1, tq1
      !> This iteration's share of the two totals, so the table says where
      !> the time went while it is going rather than only at the end. A
      !> Fock build is not the same cost twice: the SAD guess is atomic, so
      !> density screening throws nearly all of iteration 1 away, and the
      !> per-iteration cost then climbs to a plateau. Comparing two runs by
      !> total/iterations is therefore only valid at equal iteration counts.
      real(dp) :: tf_it, tx_it

      nao = b%nao
      n2 = nao*nao
      talk = opts%verbose
      if (present(comm)) talk = talk .and. comm%rank() == 0
      dft = len_trim(opts%functional) > 0
      exx = 1.0_dp
      if (dft) then
         if (.not. xc_functional_by_name(opts%functional, func)) then
            res%message = "trc_scf: unknown functional "//trim(opts%functional)
            return
         end if
         exx = func%exx
      end if
      nspin = 1
      if ((nalpha /= nbeta .or. opts%unrestricted) .and. .not. (opts%frac_occ .and. .not. opts%unrestricted)) nspin = 2
      nocc = [nalpha, nbeta]
      if (opts%frac_occ .and. nspin == 1) nocc = int(opts%nelec_frac/2.0_dp)   ! the guess fills whole orbitals; the loop refills
      occ = merge(2.0_dp, 1.0_dp, nspin == 1)
      res%nspin = nspin
      ndamp = opts%ndamp
      if (ndamp < 0) ndamp = merge(4, 0, nspin == 2)
      diis_start = opts%diis_start
      if (diis_start < 0.0_dp) diis_start = merge(0.5_dp, huge(1.0_dp), nspin == 2)

      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao), hcore(nao, nao), x(nao, nao))
      allocate (fock(nao, nao, nspin), gmat(nao, nao, nspin), vxc(nao, nao, nspin), dtot(nao, nao), jmat(nao, nao))
      allocate (dref(nao, nao, nspin), gacc(nao, nao, nspin), ddif(nao, nao))
      allocate (dold(nao, nao, nspin), errv(nao, nao, nspin))
      allocate (fstore(nao, nao, nspin, opts%ndiis), estore(nao, nao, nspin, opts%ndiis))
      allocate (w1(nao, nao), w2(nao, nao), w3(nao, nao))
      allocate (ofrac(nao)); ofrac = 0.0_dp
      allocate (bmat(opts%ndiis + 1, opts%ndiis + 1), rhs(opts%ndiis + 1))
      allocate (res%dmat(nao, nao, nspin), res%cmo(nao, nao, nspin), res%eps(nao, nspin))
      ! Only when asked for: a second ndiis-deep history the size of
      ! fstore/estore, which matters on a system the size of w172, where
      ! one nao^2*ndiis block is already hundreds of megabytes.
      if (opts%adiis) then
         allocate (fstore_a(nao, nao, nspin, opts%ndiis), dstore_a(nao, nao, nspin, opts%ndiis))
         allocate (cad(opts%ndiis))
         fstore_a = 0.0_dp; dstore_a = 0.0_dp
      end if

      t_setup = 0.0_dp; t_grid = 0.0_dp; t_fock = 0.0_dp; t_xc = 0.0_dp; t_rest = 0.0_dp
      tf_it = 0.0_dp; tx_it = 0.0_dp
      tx_pts = 0.0_dp; tx_prs = 0.0_dp
      tw0 = wall()
      ! --- once per geometry, on the host ------------------------------------
      call pl%build(b, eri_thresh_env(opts%eri_thresh))
      call pl%to_device()
      call trc_1e(b, pl, smat, tmat, vmat)
      hcore = tmat + vmat
      if (present(comm)) then
         call eri%build(b, eri_thresh_env(opts%eri_thresh), comm, batch_res=opts%batch_res)
      else
         call eri%build(b, eri_thresh_env(opts%eri_thresh), batch_res=opts%batch_res)
      end if
      res%e_nuc = nuclear_repulsion(b)
      call la%init(nao)
      call sym_orthog(la, nao, smat, x, opts%lindep_thresh, res%n_lindep)
      if (talk .and. res%n_lindep > 0) print '(a,i0,a,es9.2,a)', &
         "   canonical orthogonalisation: ", res%n_lindep, &
         " function(s) dropped, S eigenvalue below ", opts%lindep_thresh, ""
      t_setup = wall() - tw0

      if (dft) then
         tw0 = wall()
         ! The partition is split over the ranks; the summed cell weights are
         ! then rank 0's everywhere, as every reduced quantity is.
         if (present(comm)) then
            call build_dft_grid_block(b%at_r, nint(b%at_z), grid, err, comm%rank(), comm%size(), &
                                      level=opts%grid_level)
         else
            call build_dft_grid(b%at_r, nint(b%at_z), grid, err, level=opts%grid_level)
         end if
         if (err%has_error()) then
            res%message = "trc_scf: "//err%get_message()
            call eri%release(); call pl%release()
            return
         end if
         if (present(comm)) then
            if (comm%size() > 1) call sum_flat(comm, grid%weights, grid%n_points)
         end if
         call xg%build(grid%n_points, grid%coords, grid%weights, b, max_pts=opts%grid_max_pts)
         call xg%to_device()
         res%grid_r = grid%coords
         res%grid_w = grid%weights
         t_grid = wall() - tw0
      end if

      ! --- the guess, on the host, then everything goes up -------------------
      if (present(dguess)) then
         res%dmat = dguess
      else if (present(hguess)) then
         do s = 1, nspin
            call guess_from(hguess, nocc(s))
         end do
      else if (opts%guess == "core") then
         do s = 1, nspin
            call guess_from(hcore, nocc(s))
         end do
      else
         ! GWH: F_ij = 1.75/2 S_ij (H_ii + H_jj), diagonal H_ii.
         do j = 1, nao
            do i = 1, nao
               w1(i, j) = 0.875_dp*smat(i, j)*(hcore(i, i) + hcore(j, j))
            end do
            w1(j, j) = hcore(j, j)
         end do
         do s = 1, nspin
            call guess_from(w1, nocc(s))
         end do
      end if
      res%cmo = 0.0_dp; res%eps = 0.0_dp
      fock = 0.0_dp; gmat = 0.0_dp; vxc = 0.0_dp; dold = 0.0_dp; errv = 0.0_dp
      fstore = 0.0_dp; estore = 0.0_dp; dtot = 0.0_dp; jmat = 0.0_dp
      w1 = 0.0_dp; w2 = 0.0_dp; w3 = 0.0_dp

      !$acc enter data copyin(ofrac)
      !$acc enter data copyin(hcore, smat, x, res%dmat, res%cmo, res%eps) &
      !$acc            copyin(fock, gmat, vxc, dold, errv, fstore, estore, dtot, jmat, w1, w2, w3) &
      !$acc            create(dref, gacc, ddif)
      if (opts%adiis) then
         !$acc enter data copyin(fstore_a, dstore_a)
      end if

      eold = 0.0_dp
      ndiis_used = 0
      nadiis_used = 0
      diis_on = .false.
      exc = 0.0_dp; ngrid = 0.0_dp
      ! What the SCF is actually running on. Two callers can hand over the
      ! same molecule and basis name and different shells -- one had a
      ! reader that kept dead primitives -- and that only shows here.
      if (talk) print '(a,i0,a,i0,a,i0,a,i0)', "   basis: ", nao, " functions, ", b%nshell, &
         " shells, ", sum(b%sh_np), " primitives, max ", maxval(b%sh_np)
      if (talk) then
         if (dft) then
            print '(a)', "   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s)     xc(s)"
         else
            print '(a)', "   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s)"
         end if
      end if
      do it = 1, opts%max_iter
         tw0 = wall()
         ! --- two-electron part, resident ------------------------------------
         if (nspin == 1) then
            !
            ! INCREMENTAL. G(D) = G(D_ref) + G(D - D_ref), exactly, because
            ! G is linear in D. The point is the screening: the kernel
            ! screens on the density it is handed, and the difference is
            ! small and getting smaller while the density is neither.
            !
            incr_now = (opts%incremental .or. incr_env()) .and. it > 1 &
                       .and. errmax < incr_start_env(opts%incr_start) &
                       .and. mod(it - 1, max(1, incr_reset_env(opts%incr_reset))) /= 0
            if (incr_now) then
               do concurrent(i=1:nao, j=1:nao)
                  ddif(i, j) = res%dmat(i, j, 1) - dref(i, j, 1)
               end do
               !
               ! A FIXED, TIGHTER THRESHOLD ON THE INCREMENT -- not one
               ! that tracks how small ddif has become.
               !
               ! `eri`'s Schwarz/density screen drops a quartet whose
               ! bound, qq*dd, is below an ABSOLUTE threshold: right for a
               ! full build, where dd comes from a density of order one
               ! and the dropped piece really is negligible against the
               ! matrix being formed. Handed the same threshold on ddif
               ! instead, dd is the bound of a density that is getting
               ! small, so quartets the full build would never have kept
               ! survive here on sufferance, right at the edge of the
               ! test -- and WHICH quartets sit at that edge changes from
               ! one iteration's ddif to the next, because ddif's shape
               ! changes (it is a difference of two iterates, not a
               ! density that settles into one spatial pattern the way D
               ! does). So the quartets this build drops are not the same
               ! set the last one dropped, each drop is an independent
               ! truncation of a comparable size, and the errors do not
               ! average down across iterations the way the full build's
               ! -- which drops the SAME negligible quartets every time --
               ! does. Measured on w172/6-31G (450 atoms, 1950 AOs): RMS(D)
               ! reaches 1e-10 while the energy keeps bouncing by 1e-4 to
               ! 1e-6 Eh every iteration, with no downward trend, which is
               ! exactly what an error that does not shrink with ddif
               ! looks like.
               !
               ! The fix is not to chase ddif's shrinking scale -- that
               ! was tried, scaling `thresh` down by ||ddif||/||dref|| makes
               ! the screen loosen towards "keep everything" as ddif gets
               ! small, which is backwards: it buys back the accuracy at
               ! the cost of all of the speed an increment build exists
               ! for, measured on this case as 25.6 s against the 24.4 s
               ! baseline, i.e. no win at all. What works is a threshold
               ! that is simply TIGHTER, by a fixed factor, for every
               ! increment, so the per-quartet error floor itself is
               ! small enough that an unlucky run of iterations whose
               ! dropped sets do not cancel still cannot add up past
               ! `conv_energy`. `incr_thresh_scale` is that factor,
               ! default 1e-4: measured on the same case, tight enough
               ! that |dE| clears 1e-8 at every reset through iteration
               ! 20 (matching the incremental-off reference to the same
               ! energy), while the kept-quartet count stays low enough
               ! that each incremental build is still a few times
               ! cheaper than a full one, not equal to it.
               !
               block
                  real(dp) :: thresh0
                  thresh0 = eri%thresh
                  eri%thresh = thresh0*incr_thresh_scale_env(opts%incr_thresh_scale)
                  call eri%fock_resident(b, ddif, gmat(:, :, 1), k_scale=exx)
                  eri%thresh = thresh0
               end block
               do concurrent(i=1:nao, j=1:nao)
                  gacc(i, j, 1) = gacc(i, j, 1) + gmat(i, j, 1)
                  gmat(i, j, 1) = gacc(i, j, 1)
                  dref(i, j, 1) = res%dmat(i, j, 1)
               end do
            else
               call eri%fock_resident(b, res%dmat(:, :, 1), gmat(:, :, 1), k_scale=exx)
               do concurrent(i=1:nao, j=1:nao)
                  gacc(i, j, 1) = gmat(i, j, 1)
                  dref(i, j, 1) = res%dmat(i, j, 1)
               end do
            end if
         else
            do concurrent(i=1:nao, j=1:nao)
               dtot(i, j) = res%dmat(i, j, 1) + res%dmat(i, j, 2)
            end do
            call eri%fock_resident(b, dtot, jmat, k_scale=0.0_dp)
            do s = 1, 2
               if (exx /= 0.0_dp) then
                  do concurrent(i=1:nao, j=1:nao)
                     w1(i, j) = 2.0_dp*res%dmat(i, j, s)
                  end do
                  call eri%fock_resident(b, w1, gmat(:, :, s), j_scale=0.0_dp, k_scale=exx)
                  do concurrent(i=1:nao, j=1:nao)
                     gmat(i, j, s) = gmat(i, j, s) + jmat(i, j)
                  end do
               else
                  do concurrent(i=1:nao, j=1:nao)
                     gmat(i, j, s) = jmat(i, j)
                  end do
               end if
            end do
         end if
         tw1 = wall()
         t_fock = t_fock + (tw1 - tw0)
         tf_it = tw1 - tw0
         ! --- exchange-correlation: finds its density present ----------------
         if (dft) then
            if (nspin == 1) then
               call trc_xc_rks(b, xg, func, res%dmat(:, :, 1), vxc(:, :, 1), exc, ngrid, rho_tol=opts%rho_tol, &
                               comm=comm, t_points=tp1, t_pairs=tq1, la=la)
            else
               call trc_xc_uks(b, xg, func, res%dmat, vxc, exc, ngrid, rho_tol=opts%rho_tol, &
                               comm=comm, t_points=tp1, t_pairs=tq1, la=la)
            end if
            tx_pts = tx_pts + tp1; tx_prs = tx_prs + tq1
         end if
         tw0 = wall()
         t_xc = t_xc + (tw0 - tw1)
         tx_it = tw0 - tw1
         ! --- energy at this density -----------------------------------------
         e1 = 0.0_dp; e2 = 0.0_dp
         do s = 1, nspin
            e1 = e1 + la%dot(n2, res%dmat(:, :, s), hcore)
            e2 = e2 + 0.5_dp*la%dot(n2, res%dmat(:, :, s), gmat(:, :, s))
         end do
         etot = e1 + e2 + exc + res%e_nuc
         ! --- Fock and its commutator, per spin ------------------------------
         do s = 1, nspin
            if (dft) then
               do concurrent(i=1:nao, j=1:nao)
                  fock(i, j, s) = hcore(i, j) + gmat(i, j, s) + vxc(i, j, s)
               end do
            else
               do concurrent(i=1:nao, j=1:nao)
                  fock(i, j, s) = hcore(i, j) + gmat(i, j, s)
               end do
            end if
            ! w2 = F D S - S D F
            call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, fock(:, :, s), nao, res%dmat(:, :, s), nao, 0.0_dp, w1, nao)
            call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, w1, nao, smat, nao, 0.0_dp, w2, nao)
            call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, smat, nao, res%dmat(:, :, s), nao, 0.0_dp, w1, nao)
            call la%gemm('N', 'N', nao, nao, nao, -1.0_dp, w1, nao, fock(:, :, s), nao, 1.0_dp, w2, nao)
            ! errv = X^T w2 X
            call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, w2, nao, x, nao, 0.0_dp, w1, nao)
            call la%gemm('T', 'N', nao, nao, nao, 1.0_dp, x, nao, w1, nao, 0.0_dp, errv(:, :, s), nao)
         end do
         do concurrent(i=1:nao, j=1:nao, s=1:nspin)
            dold(i, j, s) = res%dmat(i, j, s)
         end do
         ! --- DIIS, once the commutator is small enough to extrapolate from --
         errmax = sqrt(la%dot(n2*nspin, errv, errv))
         if (.not. diis_on .and. it > ndamp) diis_on = errmax < diis_start
         !
         ! A/EDIIS WHILE DAMPING HAS ENDED BUT THE COMMUTATOR IS STILL TOO
         ! LARGE FOR CDIIS.
         !
         ! CDIIS minimises an error norm with no constraint on its
         ! coefficients, which is exactly why `diis_start` gates it: far
         ! from convergence the error vectors are not yet meaningful, and
         ! an unconstrained extrapolation from them can overshoot or head
         ! for a saddle rather than downhill. ADIIS (Hu & Yang, JCP 132,
         ! 054109) instead minimises an energy functional over the convex
         ! hull of the Fock/density history -- coefficients confined to
         ! the simplex {c >= 0, sum c = 1}, so the step can only move
         ! towards a density already visited, never away from all of
         ! them. That makes it reliable exactly where CDIIS is not, and
         ! slow exactly where CDIIS is fast, which is why the two SWITCH
         ! rather than blend: gmshpc measured a blended step turning a
         ! stretched-water case that converged in 21 CDIIS iterations into
         ! one that did not converge at all.
         !
         adiis_active = opts%adiis .and. .not. diis_on .and. it > ndamp
         if (adiis_active) then
            if (nadiis_used < opts%ndiis) then
               nadiis_used = nadiis_used + 1
            else
               do k = 1, opts%ndiis - 1
                  do concurrent(i=1:nao, j=1:nao, s=1:nspin)
                     fstore_a(i, j, s, k) = fstore_a(i, j, s, k + 1)
                     dstore_a(i, j, s, k) = dstore_a(i, j, s, k + 1)
                  end do
               end do
            end if
            k = nadiis_used
            do concurrent(i=1:nao, j=1:nao, s=1:nspin)
               fstore_a(i, j, s, k) = fock(i, j, s)
               dstore_a(i, j, s, k) = dold(i, j, s)
            end do
            ! One entry: the simplex is a point, c = 1 on it, a no-op.
            ! Skip the solve and leave `fock` as it is, which is also what
            ! `adiis_coeffs` would return from "start at the newest entry"
            ! after finding no improving step.
            if (nadiis_used > 1) then
               call adiis_coeffs(nadiis_used, cad(1:nadiis_used))
               do concurrent(i=1:nao, j=1:nao, s=1:nspin)
                  fock(i, j, s) = 0.0_dp
               end do
               do k = 1, nadiis_used
                  do concurrent(i=1:nao, j=1:nao, s=1:nspin)
                     fock(i, j, s) = fock(i, j, s) + cad(k)*fstore_a(i, j, s, k)
                  end do
               end do
            end if
         end if
         ! --- level shift while DIIS and A/EDIIS are off: F += mu (S - S D S / occ)
         if (.not. diis_on .and. .not. adiis_active .and. opts%level_shift > 0.0_dp) then
            do s = 1, nspin
               call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, smat, nao, res%dmat(:, :, s), nao, 0.0_dp, w1, nao)
               call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, w1, nao, smat, nao, 0.0_dp, w2, nao)
               do concurrent(i=1:nao, j=1:nao)
                  fock(i, j, s) = fock(i, j, s) + opts%level_shift*(smat(i, j) - w2(i, j)/occ)
               end do
            end do
         end if
         if (diis_on) then
            if (ndiis_used < opts%ndiis) then
               ndiis_used = ndiis_used + 1
            else
               ! Slot by slot, IN ORDER: slot k reads slot k+1 and slot k+1 is
               ! written next, so the slot loop carries a dependence. It was a
               ! `do concurrent` over k as well, and on the device that let a
               ! slot be overwritten before it was read: the history went
               ! wrong exactly when the subspace filled, and the commutator
               ! climbed for eight iterations on cholesterol before DIIS
               ! recovered. The host ran the loop in order and never saw it.
               do k = 1, opts%ndiis - 1
                  do concurrent(i=1:nao, j=1:nao, s=1:nspin)
                     fstore(i, j, s, k) = fstore(i, j, s, k + 1)
                     estore(i, j, s, k) = estore(i, j, s, k + 1)
                  end do
               end do
            end if
            k = ndiis_used
            do concurrent(i=1:nao, j=1:nao, s=1:nspin)
               fstore(i, j, s, k) = fock(i, j, s)
               estore(i, j, s, k) = errv(i, j, s)
            end do
            m = ndiis_used + 1
            bmat = -1.0_dp
            bmat(m, m) = 0.0_dp
            do j = 1, ndiis_used
               do i = 1, j
                  bmat(i, j) = la%dot(n2*nspin, estore(:, :, :, i), estore(:, :, :, j))
                  bmat(j, i) = bmat(i, j)
               end do
            end do
            call solve_diis(m, bmat(1:m, 1:m), rhs(1:m), ok)
            if (ok) then
               do concurrent(i=1:nao, j=1:nao, s=1:nspin)
                  fock(i, j, s) = 0.0_dp
               end do
               do k = 1, ndiis_used
                  do concurrent(i=1:nao, j=1:nao, s=1:nspin)
                     fock(i, j, s) = fock(i, j, s) + rhs(k)*fstore(i, j, s, k)
                  end do
               end do
            end if
         end if
         ! --- diagonalise, per spin, and the new density ---------------------
         do s = 1, nspin
            call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, fock(:, :, s), nao, x, nao, 0.0_dp, w1, nao)
            call la%gemm('T', 'N', nao, nao, nao, 1.0_dp, x, nao, w1, nao, 0.0_dp, w2, nao)
            call la%syev(nao, w2, nao, res%eps(:, s))
            call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, x, nao, w2, nao, 0.0_dp, res%cmo(:, :, s), nao)
            if (opts%frac_occ) then
               ! D = sum_i o_i c_i c_i^T through a column-scaled copy of C.
               !$acc update self(res%eps(:, s))
               if (nspin == 1) then
                  call frac_occupations(nao, res%eps(:, s), opts%nelec_frac, 2.0_dp, ofrac, nfrac)
               else
                  call frac_occupations(nao, res%eps(:, s), real(nocc(s), dp), 1.0_dp, ofrac, nfrac)
               end if
               !$acc update device(ofrac)
               do concurrent(i=1:nao, j=1:nao)
                  w1(i, j) = res%cmo(i, j, s)*sqrt(ofrac(j))
               end do
               call la%gemm('N', 'T', nao, nao, nfrac, 1.0_dp, w1, nao, w1, nao, 0.0_dp, res%dmat(:, :, s), nao)
            else if (nocc(s) > 0) then
               call la%gemm('N', 'T', nao, nao, nocc(s), occ, res%cmo(:, :, s), nao, res%cmo(:, :, s), nao, &
                            0.0_dp, res%dmat(:, :, s), nao)
            else
               do concurrent(i=1:nao, j=1:nao)
                  res%dmat(i, j, s) = 0.0_dp
               end do
            end if
         end do
         if (.not. diis_on .and. .not. adiis_active) then
            do concurrent(i=1:nao, j=1:nao, s=1:nspin)
               res%dmat(i, j, s) = opts%damp*res%dmat(i, j, s) + (1.0_dp - opts%damp)*dold(i, j, s)
            end do
         end if
         ! Every rank must hold the same density to the bit. The Fock sum is
         ! broadcast already; the XC potential is accumulated with atomics
         ! whose order differs per rank, 1e-13 in the energy on water, and a
         ! difference that small compounds through a diverging iteration.
         if (present(comm)) then
            if (comm%size() > 1) then
               !$acc update self(res%dmat)
               call bcast_flat(comm, res%dmat, n2*nspin)
               !$acc update device(res%dmat)
            end if
         end if
         do concurrent(i=1:nao, j=1:nao, s=1:nspin)
            errv(i, j, s) = res%dmat(i, j, s) - dold(i, j, s)
         end do
         drms = sqrt(la%dot(n2*nspin, errv, errv))/real(nao, dp)
         t_rest = t_rest + (wall() - tw0)
         if (talk) then
            if (dft) then
               print '(i5,f22.12,es14.4,es14.4,es14.4,f10.3,f10.3)', &
                  it, etot, etot - eold, drms, errmax, tf_it, tx_it
            else
               print '(i5,f22.12,es14.4,es14.4,es14.4,f10.3)', &
                  it, etot, etot - eold, drms, errmax, tf_it
            end if
         end if
         res%iterations = it
         ! A real density, on the way past, for the benchmarks. The model
         ! density bench_gc invents is small and banded, so it screens far
         ! harder than anything an SCF actually meets, and a Fock build
         ! timed on it flatters the code by about three times. Set
         ! TRC_DUMP_DENSITY to a path and the density leaving iteration 2 is
         ! written there, in the stream format read_density already reads.
         if (it == 2) then
            ! The density lives on the DEVICE -- every update to it since the
            ! guess was a gemm up there -- so the host copy is the SAD guess
            ! until it is pulled back. Dumping without this wrote a
            ! block-diagonal atomic density and called it iteration 2, which
            ! screens about three times harder than anything real and made
            ! every benchmark built on the file flattering and wrong.
            !$acc update self(res%dmat)
            call dump_density(nao, res%dmat(:, :, 1))
         end if
         if (it > 1 .and. abs(etot - eold) < opts%conv_energy .and. errmax < opts%conv_diis) then
            res%converged = .true.
            exit
         end if
         eold = etot
      end do

      ! The density the energy was evaluated at, and the orbitals that made it.
      !$acc update self(dold, res%cmo, res%eps)
      !$acc exit data delete(ofrac)
      !$acc exit data delete(hcore, smat, x, res%dmat, res%cmo, res%eps) &
      !$acc           delete(fock, gmat, vxc, dold, errv, fstore, estore, dtot, jmat, w1, w2, w3) &
      !$acc           delete(dref, gacc, ddif)
      if (opts%adiis) then
         !$acc exit data delete(fstore_a, dstore_a)
      end if
      res%dmat = dold
      res%energy = etot
      res%e_one = e1
      res%e_two = e2
      res%e_xc = exc
      res%nelec_grid = ngrid
      if (.not. res%converged) res%message = "trc_scf: not converged"
      if (talk) then
         print '(a)', "   time (s)   setup      grid      fock        xc   (points     pairs)   diag+diis"
         print '(a,7f10.3)', "          ", t_setup, t_grid, t_fock, t_xc, tx_pts, tx_prs, t_rest
      end if
      call la%release()
      call eri%release(); call pl%release()
      if (dft) then
         call xg%release(); call grid%destroy()
      end if

   contains

      function wall() result(t)
         real(dp) :: t
         integer(kind=8) :: cc, rate
         call system_clock(cc, rate)
         t = real(cc, dp)/real(rate, dp)
      end function wall

      ! Guess for one spin from a Fock-like h, into res%dmat(:,:,s): the
      ! eigenproblem and the products on the device, the result copied back
      ! for the resident copy made below.
      subroutine guess_from(h, nocc_s)
         real(dp), intent(in) :: h(nao, nao)
         integer, intent(in) :: nocc_s
         real(dp), allocatable :: fp(:, :), c(:, :), t(:, :), e(:), d(:, :)
         allocate (fp(nao, nao), c(nao, nao), t(nao, nao), e(nao), d(nao, nao))
         !$acc data copyin(h, x) create(fp, c, t, e) copyout(d)
         call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, h, nao, x, nao, 0.0_dp, t, nao)
         call la%gemm('T', 'N', nao, nao, nao, 1.0_dp, x, nao, t, nao, 0.0_dp, fp, nao)
         call la%syev(nao, fp, nao, e)
         call la%gemm('N', 'N', nao, nao, nao, 1.0_dp, x, nao, fp, nao, 0.0_dp, c, nao)
         if (nocc_s > 0) then
            call la%gemm('N', 'T', nao, nao, nocc_s, occ, c, nao, c, nao, 0.0_dp, d, nao)
         else
            call zero_dev(nao, d)
         end if
         !$acc end data
         res%dmat(:, :, s) = d
      end subroutine guess_from

      pure subroutine frac_occupations(n, eps, nelec, omax, o, nocc_out)
         !! Aufbau over the eigenvalues, `omax` per orbital (2 restricted, 1 per
         !! spin channel), the frontier's degenerate set (within FRAC_TOL)
         !! sharing the remainder evenly. `nocc_out` is the last orbital with
         !! any occupation, for the GEMM.
         integer, intent(in) :: n
         real(dp), intent(in) :: eps(n), nelec, omax
         real(dp), intent(out) :: o(n)
         integer, intent(out) :: nocc_out
         real(dp), parameter :: FRAC_TOL = 1.0e-3_dp
         integer :: nfull, i, g0, g1, ninside
         real(dp) :: rem, ef, share
         nfull = int(nelec/omax)
         rem = nelec - omax*real(nfull, dp)
         o = 0.0_dp
         if (nfull > 0) o(1:min(nfull, n)) = omax
         nocc_out = min(nfull, n)
         if (nfull >= n) return
         ! the frontier: the first not-full orbital when electrons are left,
         ! else the last full one (an open shell of even count, carbon 2p^2)
         if (rem > 0.0_dp .or. nfull == 0) then
            ef = eps(nfull + 1)
         else
            ef = eps(nfull)
         end if
         g0 = n + 1; g1 = 0
         do i = 1, n
            if (abs(eps(i) - ef) < FRAC_TOL) then
               g0 = min(g0, i); g1 = max(g1, i)
            end if
         end do
         if (g1 < g0) return
         ninside = count([(i <= nfull, i=g0, g1)])
         share = (omax*real(ninside, dp) + rem)/real(g1 - g0 + 1, dp)
         o(g0:g1) = share
         nocc_out = g1
      end subroutine frac_occupations

      subroutine zero_dev(n, a)
         integer, intent(in) :: n
         real(dp), intent(inout) :: a(n, n)
         integer :: i, j
         do concurrent(i=1:n, j=1:n)
            a(i, j) = 0.0_dp
         end do
      end subroutine zero_dev

      !
      ! ADIIS coefficients over the first `m` slots of fstore_a/dstore_a
      ! (host association: trc_scf_run's locals, same way `la`/`n2`/`nspin`
      ! are reached here). Slot `m` is always the newest -- the ADIIS
      ! history above is appended the same way CDIIS's is, at the top of
      ! the used range, not cycled -- so Hu & Yang's "entry n" is simply
      ! the last index.
      !
      ! Minimises, over the simplex {c >= 0, sum c = 1},
      !
      !     f(c) = 2 sum_i c_i <D_i - D_m, F_m>
      !            + sum_ij c_i c_j <D_i - D_m, F_j - F_m>
      !
      ! by projected gradient with backtracking: m <= ndiis, so this is a
      ! few hundred flops on the host once the dot products (device, via
      ! `la%dot`, the same call the CDIIS B matrix above uses) are in
      ! hand. Starting from all weight on the newest entry means a step
      ! that makes no progress degrades to the plain Roothaan step rather
      ! than to nonsense.
      !
      subroutine adiis_coeffs(m, c)
         integer, intent(in) :: m
         real(dp), intent(out) :: c(m)
         real(dp) :: dvec(m), amat(m, m), g(m), trial(m), f, ft, step
         integer :: ii, jj, itr, ls
         logical :: improved
         do ii = 1, m
            dvec(ii) = la%dot(n2*nspin, dstore_a(:, :, :, ii), fstore_a(:, :, :, m)) &
                       - la%dot(n2*nspin, dstore_a(:, :, :, m), fstore_a(:, :, :, m))
            do jj = 1, m
               amat(ii, jj) = la%dot(n2*nspin, dstore_a(:, :, :, ii), fstore_a(:, :, :, jj)) &
                              - la%dot(n2*nspin, dstore_a(:, :, :, ii), fstore_a(:, :, :, m)) &
                              - la%dot(n2*nspin, dstore_a(:, :, :, m), fstore_a(:, :, :, jj)) &
                              + la%dot(n2*nspin, dstore_a(:, :, :, m), fstore_a(:, :, :, m))
            end do
         end do
         c = 0.0_dp; c(m) = 1.0_dp
         f = adiis_value(m, dvec, amat, c)
         step = 1.0_dp
         do itr = 1, 200
            do ii = 1, m
               g(ii) = 2.0_dp*dvec(ii)
               do jj = 1, m
                  g(ii) = g(ii) + (amat(ii, jj) + amat(jj, ii))*c(jj)
               end do
            end do
            improved = .false.
            do ls = 1, 30
               trial = c - step*g
               call project_simplex(m, trial)
               ft = adiis_value(m, dvec, amat, trial)
               if (ft < f - 1.0e-14_dp) then
                  c = trial; f = ft; improved = .true.
                  exit
               end if
               step = step*0.5_dp
            end do
            if (.not. improved) exit
            step = step*1.5_dp
         end do
      end subroutine adiis_coeffs

      pure function adiis_value(m, dvec, amat, c) result(f)
         integer, intent(in) :: m
         real(dp), intent(in) :: dvec(m), amat(m, m), c(m)
         real(dp) :: f
         integer :: ii, jj
         f = 0.0_dp
         do ii = 1, m
            f = f + 2.0_dp*c(ii)*dvec(ii)
            do jj = 1, m
               f = f + c(ii)*c(jj)*amat(ii, jj)
            end do
         end do
      end function adiis_value

      ! Projection onto the simplex {c : sum c = 1, c >= 0}: the standard
      ! sort-and-threshold algorithm, which is what enforces ADIIS's
      ! non-negativity constraint, the thing that keeps it from wandering
      ! off far from convergence.
      pure subroutine project_simplex(m, c)
         integer, intent(in) :: m
         real(dp), intent(inout) :: c(m)
         real(dp) :: u(m), css, theta, t
         integer :: ii, rho
         u = c
         call sort_desc(m, u)
         css = 0.0_dp; theta = 0.0_dp; rho = 0
         do ii = 1, m
            css = css + u(ii)
            t = (css - 1.0_dp)/real(ii, dp)
            if (u(ii) - t > 0.0_dp) then
               rho = ii; theta = t
            end if
         end do
         if (rho == 0) then
            c = 1.0_dp/real(m, dp)
            return
         end if
         c = max(c - theta, 0.0_dp)
      end subroutine project_simplex

      ! Insertion sort, descending: m <= ndiis, so this is never worth more
      ! than a handful of compares.
      pure subroutine sort_desc(m, u)
         integer, intent(in) :: m
         real(dp), intent(inout) :: u(m)
         integer :: ii, jj
         real(dp) :: t
         do ii = 1, m - 1
            do jj = ii + 1, m
               if (u(jj) > u(ii)) then
                  t = u(ii); u(ii) = u(jj); u(jj) = t
               end if
            end do
         end do
      end subroutine sort_desc

   end subroutine trc_scf_run

   ! Sum a contiguous array over the ranks, then rank 0's copy everywhere.
   subroutine sum_flat(comm, g, n)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call allreduce(comm, g, op=MPI_SUM)
      call bcast(comm, g, n, 0)
   end subroutine sum_flat

   ! Broadcast rank 0's copy of a contiguous array of any rank.
   subroutine bcast_flat(comm, g, n)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call bcast(comm, g, n, 0)
   end subroutine bcast_flat

   real(dp) function nuclear_repulsion(b) result(e)
      type(trc_basis_t), intent(in) :: b
      integer :: a, c
      e = 0.0_dp
      do a = 1, b%natm
         do c = a + 1, b%natm
            e = e + b%at_z(a)*b%at_z(c)/norm2(b%at_r(:, a) - b%at_r(:, c))
         end do
      end do
   end function nuclear_repulsion

   ! X = S^(-1/2), symmetric orthogonalisation -- or, once `lindep_thresh`
   ! is positive and some S eigenvalue falls below it, CANONICAL
   ! orthogonalisation: those columns of U are set to zero rather than
   ! divided by sqrt(sqrt(w)), which is what symmetric orthogonalisation
   ! does to a near-null mode and is why it diverges first. A zeroed column
   ! of X contributes nothing to X^T F X either, so that mode decouples
   ! from the eigenproblem instead of blowing it up; the eigenvectors that
   ! come back from it are an exact, orthogonal, physically unused
   ! complement, not noise that has to be swept up afterwards. Still once
   ! per geometry, on the host side of the call.
   subroutine sym_orthog(la, n, s, x, lindep_thresh, ndropped)
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: n
      real(dp), intent(in) :: s(n, n)
      real(dp), intent(out) :: x(n, n)
      real(dp), intent(in), optional :: lindep_thresh
      integer, intent(out), optional :: ndropped
      real(dp), allocatable :: u(:, :), w(:)
      real(dp) :: thr
      integer :: i, k, ndrop
      allocate (u(n, n), w(n))
      u = s
      thr = 0.0_dp
      if (present(lindep_thresh)) thr = lindep_thresh
      !$acc data copy(u) create(w) copyout(x)
      call la%syev(n, u, n, w)
      !$acc update self(w)
      ! w is ascending (LAPACK/cuSOLVER convention), so the near-null modes
      ! -- the ones canonical orthogonalisation drops -- are the leading
      ! entries, the same ones the old unconditional check read as w(1).
      ndrop = 0
      if (thr > 0.0_dp) then
         do k = 1, n
            if (w(k) < thr) ndrop = ndrop + 1
         end do
      else if (w(1) <= 0.0_dp) then
         error stop "trc_scf: the overlap is not positive definite"
      end if
      if (present(ndropped)) ndropped = ndrop
      ! Columns scaled by w^-1/4 on the device, then X = U U^T; the first
      ! `ndrop` columns -- the near-null modes -- are zeroed instead, so
      ! they carry no weight into X^T F X and never has to divide by a
      ! near-zero w.
      !$acc update device(w)
      do concurrent(i=1:n, k=1:n)
         if (k <= ndrop) then
            u(i, k) = 0.0_dp
         else
            u(i, k) = u(i, k)/sqrt(sqrt(w(k)))
         end if
      end do
      call la%gemm('N', 'T', n, n, n, 1.0_dp, u, n, u, n, 0.0_dp, x, n)
      !$acc end data
   end subroutine sym_orthog

   !
   ! The DIIS system B c = (0, ..., 0, -1), pivoted elimination on the host,
   ! with the error block scaled to O(1) against the -1 border first. The
   ! weights are invariant under that scaling; only the Lagrange multiplier
   ! moves, and nothing reads it. `ok` is false on a pivot below 1e-14, in
   ! which case the caller keeps the unextrapolated Fock matrix.
   !
   subroutine solve_diis(m, b, c, ok)
      integer, intent(in) :: m
      real(dp), intent(in) :: b(m, m)
      real(dp), intent(out) :: c(m)
      logical, intent(out) :: ok
      real(dp), allocatable :: a(:, :)
      real(dp) :: scale, pivot, factor
      integer :: i, j, p
      allocate (a(m, m + 1))
      a(:, 1:m) = b
      a(:, m + 1) = 0.0_dp
      a(m, m + 1) = -1.0_dp
      scale = maxval(abs(b(1:m - 1, 1:m - 1)))
      if (scale > 0.0_dp) a(1:m - 1, 1:m - 1) = a(1:m - 1, 1:m - 1)/scale
      ok = .true.
      c = 0.0_dp
      do i = 1, m
         p = i
         do j = i + 1, m
            if (abs(a(j, i)) > abs(a(p, i))) p = j
         end do
         if (p /= i) a([i, p], :) = a([p, i], :)
         pivot = a(i, i)
         if (abs(pivot) < 1.0e-14_dp) then
            ok = .false.
            return
         end if
         do j = i + 1, m
            factor = a(j, i)/pivot
            a(j, :) = a(j, :) - factor*a(i, :)
         end do
      end do
      do i = m, 1, -1
         c(i) = (a(i, m + 1) - sum(a(i, i + 1:m)*c(i + 1:m)))/a(i, i)
      end do
   end subroutine solve_diis

end module trc_scf_driver
