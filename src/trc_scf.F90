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
      real(dp) :: level_shift = 1.0_dp
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
      real(dp), allocatable :: dmat(:, :, :)   !! (nao, nao, nspin); restricted holds the total
      real(dp), allocatable :: cmo(:, :, :)    !! (nao, nao, nspin) orbitals
      real(dp), allocatable :: eps(:, :)       !! (nao, nspin) orbital energies
      real(dp), allocatable :: grid_r(:, :)    !! (3, npts) the grid used, for a same-grid comparison
      real(dp), allocatable :: grid_w(:)
      character(len=256) :: message = ""
   end type trc_scf_result_t

contains

   !
   ! Run the SCF. `nalpha` and `nbeta` are the electron counts; the case is
   ! restricted when they are equal and `unrestricted` is not forced. The
   ! basis must already be on the device. `dguess` is (nao, nao, nspin) in
   ! the result's convention; without it the guess is the core Hamiltonian.
   !
   subroutine trc_scf_run(b, nalpha, nbeta, opts, res, dguess, comm)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nalpha, nbeta
      type(trc_scf_options_t), intent(in) :: opts
      type(trc_scf_result_t), intent(out) :: res
      real(dp), intent(in), optional :: dguess(:, :, :)
      !! Ranks. Every rank runs the same iteration on the same matrices;
      !! the quartets are split across them and the Fock matrix summed
      !! back, so the density stays identical everywhere. Bind devices
      !! first (trc_bind_device). The XC batches are split across ranks too.
      type(comm_t), intent(in), optional :: comm

      integer :: nao, nspin, s, it, nocc(2), ndiis_used, i, j, k, n2, m
      logical :: ok
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :), hcore(:, :), x(:, :)
      real(dp), allocatable :: fock(:, :, :), gmat(:, :, :), vxc(:, :, :), dtot(:, :), jmat(:, :)
      real(dp), allocatable :: dold(:, :, :), errv(:, :, :), fstore(:, :, :, :), estore(:, :, :, :)
      real(dp), allocatable :: w1(:, :), w2(:, :), w3(:, :), bmat(:, :), rhs(:)
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
      allocate (dold(nao, nao, nspin), errv(nao, nao, nspin))
      allocate (fstore(nao, nao, nspin, opts%ndiis), estore(nao, nao, nspin, opts%ndiis))
      allocate (w1(nao, nao), w2(nao, nao), w3(nao, nao))
      allocate (ofrac(nao)); ofrac = 0.0_dp
      allocate (bmat(opts%ndiis + 1, opts%ndiis + 1), rhs(opts%ndiis + 1))
      allocate (res%dmat(nao, nao, nspin), res%cmo(nao, nao, nspin), res%eps(nao, nspin))

      t_setup = 0.0_dp; t_grid = 0.0_dp; t_fock = 0.0_dp; t_xc = 0.0_dp; t_rest = 0.0_dp
      tx_pts = 0.0_dp; tx_prs = 0.0_dp
      tw0 = wall()
      ! --- once per geometry, on the host ------------------------------------
      call pl%build(b, opts%eri_thresh)
      call pl%to_device()
      call trc_1e(b, pl, smat, tmat, vmat)
      hcore = tmat + vmat
      if (present(comm)) then
         call eri%build(b, opts%eri_thresh, comm)
      else
         call eri%build(b, opts%eri_thresh)
      end if
      res%e_nuc = nuclear_repulsion(b)
      call la%init(nao)
      call sym_orthog(la, nao, smat, x)
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
      !$acc            copyin(fock, gmat, vxc, dold, errv, fstore, estore, dtot, jmat, w1, w2, w3)

      eold = 0.0_dp
      ndiis_used = 0
      diis_on = .false.
      exc = 0.0_dp; ngrid = 0.0_dp
      ! What the SCF is actually running on. Two callers can hand over the
      ! same molecule and basis name and different shells -- one had a
      ! reader that kept dead primitives -- and that only shows here.
      if (talk) print '(a,i0,a,i0,a,i0,a,i0)', "   basis: ", nao, " functions, ", b%nshell, &
         " shells, ", sum(b%sh_np), " primitives, max ", maxval(b%sh_np)
      if (talk) print '(a)', "   it        E(total)            dE          RMS(D)      |FDS-SDF|"
      do it = 1, opts%max_iter
         tw0 = wall()
         ! --- two-electron part, resident ------------------------------------
         if (nspin == 1) then
            call eri%fock_resident(b, res%dmat(:, :, 1), gmat(:, :, 1), k_scale=exx)
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
         ! --- level shift while DIIS is off: F += mu (S - S D S / occ) --------
         if (.not. diis_on .and. opts%level_shift > 0.0_dp) then
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
         if (.not. diis_on) then
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
         if (talk) print '(i5,f22.12,es14.4,es14.4,es14.4)', it, etot, etot - eold, drms, errmax
         res%iterations = it
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
      !$acc           delete(fock, gmat, vxc, dold, errv, fstore, estore, dtot, jmat, w1, w2, w3)
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

   ! X = S^(-1/2), symmetric orthogonalisation. Once per geometry, on the
   ! host: it is where near-null modes would be filtered, and it is not
   ! per iteration.
   subroutine sym_orthog(la, n, s, x)
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: n
      real(dp), intent(in) :: s(n, n)
      real(dp), intent(out) :: x(n, n)
      real(dp), allocatable :: u(:, :), w(:)
      integer :: i, k
      allocate (u(n, n), w(n))
      u = s
      !$acc data copy(u) create(w) copyout(x)
      call la%syev(n, u, n, w)
      !$acc update self(w)
      if (w(1) <= 0.0_dp) error stop "trc_scf: the overlap is not positive definite"
      ! Columns scaled by w^-1/4 on the device, then X = U U^T.
      !$acc update device(w)
      do concurrent(i=1:n, k=1:n)
         u(i, k) = u(i, k)/sqrt(sqrt(w(k)))
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
