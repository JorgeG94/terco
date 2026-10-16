!
!
! The SCF: a basis and an electron count in, a converged energy and
! density out, with every matrix resident on the device. Hartree-Fock or
! Kohn-Sham (any functional in trc_xc_functional), restricted or unrestricted.
!
! WHAT CROSSES THE BUS
! --------------------
! Per iteration: two energies, a norm, and a DIIS matrix of at most ten by
! ten. The density, Fock matrices, orbitals, DIIS history and orthogonaliser
! stay on the device; the Fock build is the resident driver, the XC integrator
! finds its density already present, and the products, commutator and
! eigenproblem are trc_linalg (cuBLAS/cuSOLVER under OpenACC). This follows
! metalquicha's device SCF, which measured host LAPACK plus copies at three
! times the cost of the GPU Fock build. A host build runs the same source with
! LAPACK, and the numbers agree to rounding, so the host suite tests this module.
!
! THE ITERATION
! -------------
! A few damped steps from the guess, then Pulay's DIIS on the commutator
! FDS - SDF in the orthogonal basis, per spin with one set of coefficients.
! Converged when the energy change and the RMS density change are both below
! threshold. The DIIS system (at most 11 by 11) is solved by pivoted
! elimination as in metalquicha's mqc_diis, with the error block scaled to
! order one against the constraint border first: near convergence the overlaps
! are 1e-16 and an unscaled elimination makes the weights noise, showing as a
! density that drifts from the fixed point while the energy sits flat.
!
! UNRESTRICTED
! ------------
! F_s = H + J[D_a + D_b] - a_x K[D_s] + V_xc,s. The Fock driver returns
! j J[D] - k K[D]/2 for a density it takes to be a total, so J is one call on
! D_a + D_b with k = 0 and -a_x K[D_s] one call on 2 D_s with j = 0, k = a_x:
! three builds per iteration against one restricted. A batched build over the
! three densities is the obvious improvement and is not done.
!
! SPHERICAL BASES
! ---------------
! `b%spherical` (trc_api.F90, from a JSON basis's function_type) changes the
! size of every matrix this module owns (S, H, D, F, C, eps, X, the DIIS
! history) from `b%nao` to `b%nao_sph`; everything operating on them is
! unchanged, just at a smaller `nao`. Everything touching a kernel stays
! Cartesian: `trc_1e`, `trc_eri_t%fock_resident`, and the XC integrator.
! Cartesian-sized scratch matrices (suffix `_c`: `s_c`, `h_c`, `d_c`,
! `gmat_c`, `vxc_c`, ...) sit between the kernels and the SCF's own matrices,
! built once per geometry (`s_c`/`h_c`, now `trc_scf_context`'s `smat`/`hcore`)
! or once per iteration (`d_c` from the current density; `gmat_c`/`vxc_c` from
! the Fock/XC calls, now `trc_scf_state`'s), and folded by `sph_cart_to_sph`.
! `sph` guards every call into trc_sph.F90, so the Cartesian path is the same
! code path as before spherical support existed.
!
! The XC integrator gets `d_c`, the Cartesian density rebuilt from `res%dmat`
! each iteration, not spherical AO values on the grid: the transform is shared
! with the Coulomb/exchange build, and `shell_collocate` evaluates Cartesian
! components straight off a flat table, so a second AO space there would
! duplicate trc_sph.F90's block structure.
!
! SAD and GWH reach this module already in the right space (SAD's atomic
! sub-SCFs run on a `subset` basis carrying `spherical`; GWH is built here from
! `hcore`). SAP stays Cartesian and is transformed once at the trc_capi.F90
! call site, as it has no basis handle to carry `spherical`.
!
! POLARISABLE CONTINUUM
! ---------------------
! With `pcm_model` set, trc_pcm builds the cavity, K and the nuclear potential
! once. Every iteration, at the density the energy is evaluated at, the surface
! potential of the total Cartesian density gives the charges, q.v/2 into the
! energy and the charges' full operator into each spin's Fock matrix (as
! pyscf.solvent.pcm's SCFWithSolvent, which also adds it before DIIS). It is
! rebuilt in full each iteration, independent of the incremental-ERI path.
!
! MODULES
! -------
! `trc_scf_settings` resolves `trc_scf_options_t` into what the loop reads
! (`scf_settings_t`, "cfg" below); `trc_scf_context` builds everything once
! per geometry (`scf_ctx_t`, "ctx": the pair list, S/H, the ERI engine, the
! linear-algebra handles, the orthogonaliser, and -- when asked for -- the
! PCM cavity and the DFT grid); `trc_scf_state` owns the per-iteration arrays
! (`scf_state_t`, "st") and the one enter/exit pair that puts them, `ctx`'s
! own matrices and `res` on the device together. `build_g`/`assemble_fock`
! and the DIIS push/combine helpers below stay in THIS file for now (they
! are not yet split into `trc_scf_fock`/`trc_scf_diis`): every one of them
! already takes its arrays as explicit-shape dummies rather than touching
! `ctx`/`st` inside a `do concurrent` body, which is what lets this compile
! correctly once those arrays are struct components instead of plain locals.
!
module trc_scf_driver
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t
   use trc_sph, only: sph_cart_to_sph, sph_sph_to_cart
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t
   use trc_vv10, only: trc_vv10_rks, trc_vv10_uks
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_xc, only: trc_xc_rks, trc_xc_uks, XC_RHO_TOL
   use trc_linalg, only: trc_linalg_t
   use trc_scf_settings, only: trc_scf_options_t, trc_scf_result_t, scf_settings_t
   use trc_scf_context, only: scf_ctx_t, bcast_flat
   use trc_scf_state, only: scf_state_t, scf_timers_t
   use pic_types, only: default_int
   use pic_logger, only: global_logger
   use pic_mpi_lib, only: comm_t
   implicit none
   private

   public :: trc_scf_options_t, trc_scf_result_t, trc_scf_run

contains

   subroutine dump_density(nao, d)
      !! The density at a chosen iteration, written where a benchmark can read
      !! it back. Off unless TRC_DUMP_DENSITY names a file. Stream format:
      !! an 8-byte nao, then the matrix, which is what the test readers expect.
      integer, intent(in) :: nao
      real(dp), intent(in) :: d(nao, nao)
      character(len=512) :: path
      integer :: u, ios
      integer(int64) :: prev
      call get_environment_variable("TRC_DUMP_DENSITY", path)
      if (len_trim(path) == 0) return
      ! The atomic SCFs of a SAD guess come through here too, and would
      ! overwrite the molecular density with a 5x5 one. Keep the largest.
      open (newunit=u, file=trim(path), access="stream", form="unformatted", &
            status="old", action="read", iostat=ios)
      if (ios == 0) then
         read (u, iostat=ios) prev
         close (u)
         if (ios == 0) then
            if (prev >= int(nao, int64)) return
         end if
      end if
      open (newunit=u, file=trim(path), access="stream", form="unformatted", &
            status="replace", action="write", iostat=ios)
      if (ios /= 0) return
      write (u) int(nao, int64)
      write (u) d
      close (u)
      print "(a,a,a,i0,a)", "   [density dumped to ", trim(path), ", nao ", nao, "]"
   end subroutine dump_density

   !
   ! Run the SCF. `nalpha` and `nbeta` are the electron counts; the case is
   ! restricted when they are equal and `unrestricted` is not forced. The
   ! basis must already be on the device. `dguess` is (nao, nao, nspin) in
   ! the result's convention; without it the guess is the core Hamiltonian.
   !
   subroutine trc_scf_run(b, nalpha, nbeta, opts, res, dguess, comm, hguess, la)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nalpha, nbeta
      type(trc_scf_options_t), intent(in) :: opts
      type(trc_scf_result_t), intent(out) :: res
      real(dp), intent(in), optional :: dguess(:, :, :)
      real(dp), intent(in), optional :: hguess(:, :)
      !! A guess HAMILTONIAN to occupy instead of hcore or the GWH matrix --
      !! how SAP gets in (see trc_sap), which cannot come through `dguess`
      !! because what it produces is a one-electron operator and not a
      !! density. Ignored if `dguess` is also given: a density the caller
      !! already has is more specific than a matrix to make one from.
      type(comm_t), intent(in), optional :: comm
      !! Ranks. Every rank runs the same iteration on the same matrices;
      !! the quartets are split across them and the Fock matrix summed
      !! back, so the density stays identical everywhere. Bind devices
      !! first (trc_bind_device). The XC batches are split across ranks too.
      type(trc_linalg_t), intent(inout), optional, target :: la
      !! cuBLAS/cuSOLVER handles made once by the caller and carried through
      !! the calculation (trc_linalg's init is create-once and grows the
      !! workspace to nao). Absent, the run makes its own and releases them.

      integer :: nao, nao_c, nspin, s, it, ndiis_used, i, j, k, n2, m
      integer :: nocc(2)
      integer :: nxc, nu_c, nlr_c, nin_c, npcm   !! full size, or 1 when unused (see trc_scf_state%alloc)
      logical :: ok, sph
      type(scf_settings_t) :: cfg
      type(scf_ctx_t), target :: ctx
      !! TARGET: ctx%setup points ctx%lp at ctx%la_own when the caller gives
      !! no linalg handles of its own; that pointer must stay valid for the
      !! rest of this routine, which needs ctx itself to carry TARGET.
      type(scf_state_t) :: st
      type(scf_timers_t) :: tmr
      character(len=256) :: msg
      logical :: incr_now
      real(dp), allocatable :: fstore(:, :, :, :), estore(:, :, :, :)
      real(dp), allocatable :: bmat(:, :), rhs(:)
      real(dp), allocatable :: fstore_a(:, :, :, :), dstore_a(:, :, :, :), cad(:)
      !! ADIIS history (fock, density) and its own slot count, kept apart
      !! from the CDIIS `fstore`/`estore` above: ADIIS runs exactly in the
      !! window CDIIS has not opened yet, so the two never need the same
      !! entries, and keeping them separate means `opts%adiis = .false.`
      !! leaves the CDIIS path untouched, to the bit.
      integer :: nadiis_used
      logical :: adiis_active
      logical :: nlc
      real(dp) :: e_nlc
      type(trc_xc_functional_t) :: func
      type(trc_linalg_t), pointer :: lp
      type(error_t) :: err
      logical :: pcm_on
      real(dp) :: epcm
      logical :: dft, talk, diis_on, root, incr_on
      real(dp) :: ctl(2)
      real(dp) :: exx_lr, omega
      !! A range-separated hybrid's extra long-range exchange fraction and
      !! its omega: K_total = exx*K + exx_lr*K_erf(omega). Zero otherwise.
      real(dp) :: exx, etot, eold, drms, e1, e2, exc, ngrid, occ, errmax
      integer :: nfrac, ndamp
      real(dp) :: diis_start
      real(dp) :: tw0, tw1
      ! The TRC_* experiment hooks (incr_*_env, eri_thresh_env), read once
      ! per run: they only reread the environment, and the loop asked for
      ! the same answer every iteration.
      real(dp) :: eri_thresh_run, incr_start_run, incr_scale_run
      integer :: incr_reset_run
      !! This iteration's share of the two totals, so the table says where
      !! the time went while it is going rather than only at the end. A
      !! Fock build is not the same cost twice: the SAD guess is atomic, so
      !! density screening throws nearly all of iteration 1 away, and the
      !! per-iteration cost then climbs to a plateau. Comparing two runs by
      !! total/iterations is therefore only valid at equal iteration counts.

      ! Everything `opts` resolves to before the loop starts: the functional
      ! lookup and its rejections, the spin rule, the by-spin damping/DIIS
      ! defaults, and the env-hook overrides (trc_scf_settings). `cfg` is
      ! immutable from here on; the local copies below are what the rest of
      ! this routine already read them as, unchanged by the move.
      ok = cfg%resolve(opts, nalpha, nbeta, msg=msg, comm=comm)
      if (.not. ok) then
         res%message = msg
         return
      end if
      talk = cfg%talk
      root = cfg%root
      dft = cfg%dft
      nlc = cfg%nlc
      exx = cfg%exx
      exx_lr = cfg%exx_lr
      omega = cfg%omega
      func = cfg%func
      nspin = cfg%nspin
      nocc = cfg%nocc
      occ = cfg%occ
      res%nspin = nspin
      ndamp = cfg%ndamp
      diis_start = cfg%diis_start
      pcm_on = cfg%pcm_on
      incr_on = cfg%incr_on
      eri_thresh_run = cfg%eri_thresh_run
      incr_start_run = cfg%incr_start_run
      incr_reset_run = cfg%incr_reset_run
      incr_scale_run = cfg%incr_scale_run
      ! The cavity, K and its LU are built and solved whole on every rank
      ! (replicated: cheap next to the integrals); pcm%operator splits the
      ! two integral passes over the ranks itself, given `comm`.
      epcm = 0.0_dp

      ! --- everything once per geometry: pairs, 1e, S/H fold, ERI, linalg,
      ! orthogonaliser, and (when asked for) the PCM cavity and DFT grid ---
      call ctx%setup(b, opts, cfg, err, comm=comm, la=la)
      if (err%has_error()) then
         res%message = "trc_scf: "//err%get_message()
         return
      end if
      nao = ctx%nao
      nao_c = ctx%nao_c
      sph = ctx%sph
      n2 = nao*nao
      lp => ctx%lp
      res%e_nuc = ctx%e_nuc
      res%n_lindep = ctx%n_lindep
      if (pcm_on) res%pcm_npoints = ctx%pcm%surface%n_points
      if (dft) then
         res%grid_r = ctx%grid%coords
         res%grid_w = ctx%grid%weights
         if (nlc) then
            res%grid_r_nlc = ctx%grid_nlc%coords
            res%grid_w_nlc = ctx%grid_nlc%weights
         end if
      end if

      ! ONLY WHAT THIS RUN USES IS FULL SIZE (trc_scf_state%alloc owns the
      ! array sizes; `nu_c`/`nlr_c`/`nin_c`/`nxc`/`npcm` below are the same
      ! sizes again, for the calls into `build_g`/`assemble_fock`, which take
      ! explicit-shape dummies and so must be told how big a 1 x 1 "unused"
      ! array really is). One that a run never touches is nao^2 of device
      ! memory for nothing: on ionic05 in aug-cc-pVTZ (4255 functions, 4950
      ! Cartesian) an HF run carried ~2.4 GB of XC, unrestricted,
      ! range-separated, incremental and PCM matrices it had no use for.
      nxc = merge(nao, 1, dft)
      nu_c = merge(nao_c, 1, nspin == 2)
      nlr_c = merge(nao_c, 1, exx_lr /= 0.0_dp)
      nin_c = merge(nao_c, 1, incr_on)
      npcm = merge(nao, 1, pcm_on)
      call st%alloc(ctx, cfg)
      allocate (fstore(nao, nao, nspin, opts%ndiis), estore(nao, nao, nspin, opts%ndiis))
      fstore = 0.0_dp
      estore = 0.0_dp
      allocate (bmat(opts%ndiis + 1, opts%ndiis + 1), rhs(opts%ndiis + 1))
      allocate (res%dmat(nao, nao, nspin), res%cmo(nao, nao, nspin), res%eps(nao, nspin))
      ! Only when asked for: a second ndiis-deep history the size of
      ! fstore/estore, which matters on a system the size of w172, where
      ! one nao^2*ndiis block is already hundreds of megabytes.
      if (opts%adiis) then
         allocate (fstore_a(nao, nao, nspin, opts%ndiis), dstore_a(nao, nao, nspin, opts%ndiis))
         allocate (cad(opts%ndiis))
         fstore_a = 0.0_dp
         dstore_a = 0.0_dp
      end if

      ! --- the guess, on the host, then everything goes up -------------------
      if (present(dguess)) then
         res%dmat = dguess
      else if (present(hguess)) then
         do s = 1, nspin
            call guess_from(hguess, nocc(s), res%dmat(:, :, s))
         end do
      else if (opts%guess == "core") then
         do s = 1, nspin
            call guess_from(ctx%hcore, nocc(s), res%dmat(:, :, s))
         end do
      else
         ! GWH: F_ij = 1.75/2 S_ij (H_ii + H_jj), diagonal H_ii.
         do j = 1, nao
            do i = 1, nao
               st%w1(i, j) = 0.875_dp*ctx%smat(i, j)*(ctx%hcore(i, i) + ctx%hcore(j, j))
            end do
            st%w1(j, j) = ctx%hcore(j, j)
         end do
         do s = 1, nspin
            call guess_from(st%w1, nocc(s), res%dmat(:, :, s))
         end do
      end if
      res%cmo = 0.0_dp
      res%eps = 0.0_dp

      call st%to_device(ctx, res)
      !$acc enter data copyin(fstore, estore)
      if (opts%adiis) then
         !$acc enter data copyin(fstore_a, dstore_a)
      end if

      eold = 0.0_dp
      ! Read after the loop (res%e_one/e_two) even when it never runs
      ! (max_iter=0): the loop body is the only other place that sets them.
      e1 = 0.0_dp
      e2 = 0.0_dp
      ! Read by the first iteration's `incr_now` before any commutator exists
      ! (Fortran's .and. does not short-circuit), so give it a value.
      errmax = huge(1.0_dp)
      ndiis_used = 0
      nadiis_used = 0
      diis_on = .false.
      exc = 0.0_dp
      ngrid = 0.0_dp
      ! What the SCF is actually running on. Two callers can hand over the
      ! same molecule and basis name and different shells -- one had a
      ! reader that kept dead primitives -- and that only shows here.
      if (talk) print "(a,i0,a,i0,a,i0,a,i0)", "   basis: ", nao, " functions, ", b%nshell, &
         " shells, ", sum(b%sh_np), " primitives, max ", maxval(b%sh_np)
      if (talk) then
         if (dft) then
            print "(a)", "   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s)     xc(s)"
         else
            print "(a)", "   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s)"
         end if
      end if
      do it = 1, opts%max_iter
         tw0 = wall()
         ! --- D_cart := T D_sph T^T, the one transform every Cartesian
         ! quantity below (the Fock build, the XC density) is built from.
         ! Not taken at all when `sph` is false -- see this module's own
         ! header on why that keeps the Cartesian path byte-for-byte
         ! unchanged rather than merely numerically equivalent.
         do s = 1, nspin
            call to_cart(res%dmat(:, :, s), st%d_c(:, :, s))
         end do
         ! --- two-electron part, resident, entirely in Cartesian space -------
         !
         ! INCREMENTAL. G(D) = G(D_ref) + G(D - D_ref), exactly, because G is
         ! linear in D. The point is the screening: the kernel screens on the
         ! density it is handed, and the difference is small and getting
         ! smaller while the density is neither. See `build_g`'s own header.
         !
         incr_now = incr_on .and. it > 1 &
                    .and. errmax < incr_start_run &
                    .and. mod(it - 1, max(1, incr_reset_run)) /= 0
         call build_g(b, ctx%eri, nao_c, nspin, nu_c, nlr_c, nin_c, incr_now, incr_on, &
                      exx, exx_lr, omega, incr_scale_run, &
                      st%d_c, st%gmat_c, st%dtot_c, st%jmat_c, st%w1_c, st%klr_c, &
                      st%gacc_c, st%dref_c, st%ddif_c)
         ! --- G_sph := T^T G_cart T, back up to the space the rest of the
         ! iteration runs in.
         do s = 1, nspin
            call to_sph(st%gmat_c(:, :, s), st%gmat(:, :, s))
         end do
         tw1 = wall()
         tmr%t_fock = tmr%t_fock + (tw1 - tw0)
         tmr%tf_it = tw1 - tw0
         ! --- exchange-correlation: finds its (Cartesian) density present ----
         if (dft) then
            if (nspin == 1) then
               call trc_xc_rks(b, ctx%xg, func, st%d_c(:, :, 1), st%vxc_c(:, :, 1), exc, ngrid, &
                               rho_tol=opts%rho_tol, comm=comm, t_points=tmr%tx_pts, t_pairs=tmr%tx_prs, la=lp)
            else
               call trc_xc_uks(b, ctx%xg, func, st%d_c, st%vxc_c, exc, ngrid, rho_tol=opts%rho_tol, &
                               comm=comm, t_points=tmr%tx_pts, t_pairs=tmr%tx_prs, la=lp)
            end if
            ! VV10 adds its energy to exc and its potential into vxc_c, on
            ! its own (coarser) grid, before the spherical transform below.
            if (nlc) then
               if (nspin == 1) then
                  call trc_vv10_rks(b, ctx%xg_nlc, st%d_c(:, :, 1), func%nlc_b, func%nlc_c, e_nlc, st%vxc_c(:, :, 1))
               else
                  call trc_vv10_uks(b, ctx%xg_nlc, st%d_c(:, :, 1), st%d_c(:, :, 2), func%nlc_b, func%nlc_c, e_nlc, &
                                    st%vxc_c(:, :, 1), st%vxc_c(:, :, 2))
               end if
               exc = exc + e_nlc
            end if
            do s = 1, nspin
               call to_sph(st%vxc_c(:, :, s), st%vxc(:, :, s))
            end do
         end if
         tw0 = wall()
         tmr%t_xc = tmr%t_xc + (tw0 - tw1)
         tmr%tx_it = tw0 - tw1
         ! --- the continuum, at the TOTAL density (pyscf's choice for UKS
         ! too: one set of charges, one operator added to both spins) -----
         if (pcm_on) then
            ! An absent `comm` passes through as absent: one rank, unsplit.
            if (nspin == 1) then
               call ctx%pcm%operator(b, ctx%pl, st%d_c(:, :, 1), st%vpcm_c, epcm, comm)
            else
               call ctx%pcm%operator(b, ctx%pl, st%dtot_c, st%vpcm_c, epcm, comm)
            end if
            call to_sph(st%vpcm_c, st%vpcm)
            tw1 = wall()
            tmr%t_pcm = tmr%t_pcm + (tw1 - tw0)
            tw0 = tw1
         end if
         ! --- energy at this density -----------------------------------------
         e1 = 0.0_dp
         e2 = 0.0_dp
         do s = 1, nspin
            e1 = e1 + lp%dot(n2, res%dmat(:, :, s), ctx%hcore)
            e2 = e2 + 0.5_dp*lp%dot(n2, res%dmat(:, :, s), st%gmat(:, :, s))
         end do
         etot = e1 + e2 + exc + res%e_nuc + epcm
         ! --- Fock and its commutator, per spin ------------------------------
         call assemble_fock(nao, nxc, npcm, nspin, dft, pcm_on, ctx%hcore, st%gmat, st%vxc, st%vpcm, st%fock)
         do s = 1, nspin
            ! w2 = F D S - S D F
            call lp%gemm("N", "N", nao, nao, nao, 1.0_dp, st%fock(:, :, s), nao, res%dmat(:, :, s), nao, &
                         0.0_dp, st%w1, nao)
            call lp%gemm("N", "N", nao, nao, nao, 1.0_dp, st%w1, nao, ctx%smat, nao, 0.0_dp, st%w2, nao)
            call lp%gemm("N", "N", nao, nao, nao, 1.0_dp, ctx%smat, nao, res%dmat(:, :, s), nao, 0.0_dp, st%w1, nao)
            call lp%gemm("N", "N", nao, nao, nao, -1.0_dp, st%w1, nao, st%fock(:, :, s), nao, 1.0_dp, st%w2, nao)
            ! errv = X^T w2 X
            call lp%gemm("N", "N", nao, nao, nao, 1.0_dp, st%w2, nao, ctx%x, nao, 0.0_dp, st%w1, nao)
            call lp%gemm("T", "N", nao, nao, nao, 1.0_dp, ctx%x, nao, st%w1, nao, 0.0_dp, st%errv(:, :, s), nao)
         end do
         ! `dold` must hold `res%dmat` from BEFORE DIIS/diagonalisation touch
         ! it -- see trc_scf_state's own header.
         call st%save_dold(nao, nspin, res%dmat)
         ! --- DIIS, once the commutator is small enough to extrapolate from --
         errmax = sqrt(lp%dot(n2*nspin, st%errv, st%errv))
         ! Over ranks the Fock matrices agree only to the XC atomics' 1e-13,
         ! so the numbers every branch below and the convergence test read
         ! are rank 0's, on every rank: one rank leaving the loop while
         ! another iterates on would hang the next collective.
         if (present(comm)) then
            if (comm%size() > 1) then
               ctl = [etot, errmax]
               call bcast_flat(comm, ctl, 2)
               etot = ctl(1)
               errmax = ctl(2)
            end if
         end if
         if (.not. diis_on .and. it > ndamp) diis_on = errmax < diis_start
         !
         ! A/EDIIS while damping has ended but the commutator is still too large for
         ! CDIIS. CDIIS has no constraint on its coefficients, hence the
         ! `diis_start` gate: far from convergence it can overshoot. ADIIS (Hu &
         ! Yang, JCP 132, 054109) minimises an energy functional over the convex
         ! hull of the history, so it only moves towards a visited density: reliable
         ! where CDIIS is not and slow where CDIIS is fast, hence the two switch
         ! rather than blend (gmshpc: a blended step turned a stretched-water case
         ! converging in 21 CDIIS iterations into one that never converged).
         !
         adiis_active = opts%adiis .and. .not. diis_on .and. it > ndamp
         ! The new density is rank 0's alone: extrapolation, eigensolve and density
         ! run on one device and the density is broadcast, since dsyevd need not
         ! return the same vectors on two devices. Orbitals and energies come over
         ! once, after the loop.
         if (root) then
            if (adiis_active) then
               if (nadiis_used < opts%ndiis) then
                  nadiis_used = nadiis_used + 1
               else
                  call shift_history(opts%ndiis, fstore_a, dstore_a)
               end if
               call diis_store(nao, nspin, nadiis_used, st%fock, st%dold, fstore_a, dstore_a)
               ! One entry: the simplex is a point, c = 1 on it, a no-op.
               ! Skip the solve and leave `fock` as it is, which is also what
               ! `adiis_coeffs` would return from "start at the newest entry"
               ! after finding no improving step.
               if (nadiis_used > 1) then
                  call adiis_coeffs(nadiis_used, cad(1:nadiis_used))
                  call diis_combine(nao, nspin, nadiis_used, cad(1:nadiis_used), fstore_a, st%fock)
               end if
            end if
            ! --- level shift while DIIS and A/EDIIS are off: F += mu (S - S D S / occ)
            if (.not. diis_on .and. .not. adiis_active .and. opts%level_shift > 0.0_dp) then
               call level_shift(lp, nao, nspin, opts%level_shift, occ, ctx%smat, res%dmat, st%w1, st%w2, st%fock)
            end if
            if (diis_on) then
               if (ndiis_used < opts%ndiis) then
                  ndiis_used = ndiis_used + 1
               else
                  call shift_history(opts%ndiis, fstore, estore)
               end if
               call diis_store(nao, nspin, ndiis_used, st%fock, st%errv, fstore, estore)
               !
               ! An ill-conditioned subspace loses its oldest entry, not the step.
               ! A failed solve followed by a plain Roothaan step is not neutral: on
               ! water/def2-SVP PBE0 the history holds errors of order 1 next to the
               ! newest at 1e-7, B's scaled pivots drop under solve_diis's 1e-14, and
               ! plain steps grow the commutator x1.28 per iteration along one
               ! direction, which keeps B singular: 65 iterations of runaway. Dropping
               ! the oldest entries until the system solves keeps the extrapolation;
               ! two distinct error vectors always solve, and that step cancels one
               ! growing mode.
               !
               do
                  m = ndiis_used + 1
                  bmat = -1.0_dp
                  bmat(m, m) = 0.0_dp
                  do j = 1, ndiis_used
                     do i = 1, j
                        bmat(i, j) = lp%dot(n2*nspin, estore(:, :, :, i), estore(:, :, :, j))
                        bmat(j, i) = bmat(i, j)
                     end do
                  end do
                  call solve_diis(m, bmat(1:m, 1:m), rhs(1:m), ok)
                  if (ok .or. ndiis_used <= 2) exit
                  call shift_history(ndiis_used, fstore, estore)
                  ndiis_used = ndiis_used - 1
               end do
               if (ok) call diis_combine(nao, nspin, ndiis_used, rhs(1:ndiis_used), fstore, st%fock)
            end if
            ! --- diagonalise, per spin, and the new density ---------------------
            do s = 1, nspin
               call orbitals_and_density(st%fock(:, :, s), nocc(s), opts%frac_occ, st%w1, st%w2, &
                                         res%eps(:, s), res%cmo(:, :, s), res%dmat(:, :, s), st%ofrac)
            end do
            if (.not. diis_on .and. .not. adiis_active) call damp(nao, nspin, opts%damp, st%dold, res%dmat)
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
         call st%update_errv(nao, nspin, res%dmat)
         drms = sqrt(lp%dot(n2*nspin, st%errv, st%errv))/real(nao, dp)
         tmr%t_rest = tmr%t_rest + (wall() - tw0)
         if (talk) then
            if (dft) then
               print "(i5,f22.12,es14.4,es14.4,es14.4,f10.3,f10.3)", &
                  it, etot, etot - eold, drms, errmax, tmr%tf_it, tmr%tx_it
            else
               print "(i5,f22.12,es14.4,es14.4,es14.4,f10.3)", &
                  it, etot, etot - eold, drms, errmax, tmr%tf_it
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
      !$acc update self(st%dold, res%cmo, res%eps)
      if (present(comm)) then
         if (comm%size() > 1) then
            call bcast_flat(comm, res%cmo, n2*nspin)
            call bcast_flat(comm, res%eps, nao*nspin)
         end if
      end if
      call st%release_device(ctx, res)
      !$acc exit data delete(fstore, estore)
      if (opts%adiis) then
         !$acc exit data delete(fstore_a, dstore_a)
      end if
      res%dmat = st%dold
      res%energy = etot
      res%e_one = e1
      res%e_two = e2
      res%e_xc = exc
      res%nelec_grid = ngrid
      res%e_pcm = epcm
      res%pcm_qtotal = ctx%pcm%q_total
      if (.not. res%converged) res%message = "trc_scf: not converged"
      if (talk) then
         call report_timing()
         if (pcm_on) print "(a,f20.12,a,f10.6,a,f10.3,a,f10.4,a)", "   pcm: E_diel = ", epcm, &
            "  q_total = ", ctx%pcm%q_total, "  setup ", ctx%t_pcm_setup, " s, ", &
            tmr%t_pcm/real(max(1, res%iterations), dp), " s per iteration"
      end if
      call ctx%release()

   contains

      subroutine report_timing()
         !! Where the run's time went, through pic's logger at PERFORMANCE,
         !! which its default INFO level shows: the phases, then setup and
         !! the ERI engine's build broken into their steps.
         character(len=160) :: line
         write (line, "(a,7f10.3)") "          ", ctx%t_setup, ctx%t_grid, tmr%t_fock, tmr%t_xc, &
            tmr%tx_pts, tmr%tx_prs, tmr%t_rest
         call global_logger%performance( &
            "   time (s)   setup      grid      fock        xc   (points     pairs)   diag+diis")
         call global_logger%performance(trim(line))
         call global_logger%performance( &
            "   setup (s)  pairs        1e   sph fold       eri    linalg    orthog")
         write (line, "(a,6f10.3)") "          ", ctx%ts
         call global_logger%performance(trim(line))
         call global_logger%performance( &
            "   eri (s) hgp pairs   schwarz      bins  quartets    upload")
         write (line, "(a,5f10.3)") "          ", ctx%eri%t_build
         call global_logger%performance(trim(line))
      end subroutine report_timing

      function wall() result(t)
         real(dp) :: t
         integer(int64) :: cc, rate
         call system_clock(cc, rate)
         t = real(cc, dp)/real(rate, dp)
      end function wall

      ! Guess for one spin from a Fock-like h, into `dout`: the eigenproblem
      ! and the products on the device, the result copied back for the
      ! resident copy made below. h and x are host arrays here, copied up for
      ! the region; the loop runs the same sequence on resident ones.
      subroutine guess_from(h, nocc_s, dout)
         real(dp), intent(in) :: h(nao, nao)
         integer, intent(in) :: nocc_s
         real(dp), intent(out) :: dout(nao, nao)
         real(dp), allocatable :: fp(:, :), c(:, :), t(:, :), e(:), d(:, :)
         allocate (fp(nao, nao), c(nao, nao), t(nao, nao), e(nao), d(nao, nao))
         !$acc data copyin(h, ctx%x) create(fp, c, t, e) copyout(d)
         call orbitals_and_density(h, nocc_s, .false., t, fp, e, c, d, st%ofrac)
         !$acc end data
         dout = d
      end subroutine guess_from

      subroutine orbitals_and_density(f, nocc_s, frac, t, fp, eps_s, cmo_s, d_s, ofrac)
         !! One spin's step from a Fock-like matrix to its density: F' = X^T F X,
         !! its eigenvectors, C = X C', then D = occ C(:,1:nocc) C(:,1:nocc)^T
         !! (or, with `frac`, the fractional occupations of the free-atom runs).
         !! Everything is on the device already (the loop's resident arrays, or
         !! the data region of `guess_from`); `t` and `fp` are the caller's
         !! n x n workspaces, which hold nothing on return that matters. `ofrac`
         !! is an explicit argument (not host-associated from `st%ofrac`
         !! directly) so the `frac` branch's `do concurrent` below never names
         !! a derived-type component.
         real(dp), intent(in) :: f(nao, nao)
         integer, intent(in) :: nocc_s
         logical, intent(in) :: frac
         real(dp), intent(inout) :: t(nao, nao), fp(nao, nao)
         real(dp), intent(inout) :: eps_s(nao), cmo_s(nao, nao), d_s(nao, nao)
         real(dp), intent(inout) :: ofrac(nao)
         integer :: i, j
         call lp%gemm("N", "N", nao, nao, nao, 1.0_dp, f, nao, ctx%x, nao, 0.0_dp, t, nao)
         call lp%gemm("T", "N", nao, nao, nao, 1.0_dp, ctx%x, nao, t, nao, 0.0_dp, fp, nao)
         call lp%syev(nao, fp, nao, eps_s)
         call lp%gemm("N", "N", nao, nao, nao, 1.0_dp, ctx%x, nao, fp, nao, 0.0_dp, cmo_s, nao)
         if (frac) then
            ! D = sum_i o_i c_i c_i^T through a column-scaled copy of C.
            !$acc update self(eps_s)
            if (nspin == 1) then
               call frac_occupations(nao, eps_s, opts%nelec_frac, 2.0_dp, ofrac, nfrac)
            else
               call frac_occupations(nao, eps_s, real(nocc_s, dp), 1.0_dp, ofrac, nfrac)
            end if
            !$acc update device(ofrac)
            do concurrent(i=1:nao, j=1:nao)
               t(i, j) = cmo_s(i, j)*sqrt(ofrac(j))
            end do
            call lp%gemm("N", "T", nao, nao, nfrac, 1.0_dp, t, nao, t, nao, 0.0_dp, d_s, nao)
         else if (nocc_s > 0) then
            call lp%gemm("N", "T", nao, nao, nocc_s, occ, cmo_s, nao, cmo_s, nao, 0.0_dp, d_s, nao)
         else
            call zero_dev(nao, d_s)
         end if
      end subroutine orbitals_and_density

      subroutine to_cart(a, a_c)
         !! A matrix from the SCF's space to the Cartesian one the kernels run
         !! in: T A T^T when the basis is spherical, a plain copy when it is not
         !! (no transform call at all on the Cartesian path, see this module's
         !! header). `a` and `a_c` are the same size then.
         real(dp), intent(in) :: a(nao, nao)
         real(dp), intent(out) :: a_c(nao_c, nao_c)
         integer :: i, j
         if (sph) then
            call sph_sph_to_cart(b, a, a_c)
         else
            do concurrent(i=1:nao_c, j=1:nao_c)
               a_c(i, j) = a(i, j)
            end do
         end if
      end subroutine to_cart

      subroutine to_sph(a_c, a)
         !! The reverse fold: T^T A_cart T into the SCF's space, or a plain copy.
         real(dp), intent(in) :: a_c(nao_c, nao_c)
         real(dp), intent(out) :: a(nao, nao)
         integer :: i, j
         if (sph) then
            call sph_cart_to_sph(b, a_c, a)
         else
            do concurrent(i=1:nao, j=1:nao)
               a(i, j) = a_c(i, j)
            end do
         end if
      end subroutine to_sph

      subroutine shift_history(n, a, bb)
         !! Slot k <- slot k+1 for k = 1..n-1 in two histories at once, which
         !! drops the oldest entry and frees slot n.
         integer, intent(in) :: n
         real(dp), intent(inout) :: a(nao, nao, nspin, *), bb(nao, nao, nspin, *)
         integer :: i, j, s, k
         ! Slot by slot, in order: slot k reads slot k+1, so the slot loop
         ! carries a dependence. As a `do concurrent` over k, a device slot
         ! was overwritten before being read and the commutator climbed for
         ! eight iterations on cholesterol; the host never saw it.
         do k = 1, n - 1
            do concurrent(i=1:nao, j=1:nao, s=1:nspin)
               a(i, j, s, k) = a(i, j, s, k + 1)
               bb(i, j, s, k) = bb(i, j, s, k + 1)
            end do
         end do
      end subroutine shift_history

      pure subroutine frac_occupations(n, eps, nelec, omax, o, nocc_out)
         !! Aufbau over the eigenvalues, `omax` per orbital (2 restricted, 1 per
         !! spin channel), the frontier's degenerate set (within FRAC_TOL)
         !! sharing the remainder evenly. `nocc_out` is the last orbital with
         !! any occupation, for the GEMM.
         integer, intent(in) :: n
         real(dp), intent(in) :: nelec, omax
         real(dp), intent(in) :: eps(n)
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
         g0 = n + 1
         g1 = 0
         do i = 1, n
            if (abs(eps(i) - ef) < FRAC_TOL) then
               g0 = min(g0, i)
               g1 = max(g1, i)
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
      ! ADIIS coefficients over the first `m` slots of fstore_a/dstore_a (host
      ! association). Slot `m` is the newest, so Hu & Yang's "entry n" is the last
      ! index. Minimises over the simplex {c >= 0, sum c = 1}
      !
      !     f(c) = 2 sum_i c_i <D_i - D_m, F_m>
      !            + sum_ij c_i c_j <D_i - D_m, F_j - F_m>
      !
      ! by projected gradient with backtracking (a few hundred host flops once the
      ! device dot products are in hand). Starting with all weight on the newest
      ! entry means no progress degrades to the plain Roothaan step.
      !
      subroutine adiis_coeffs(m, c)
         integer, intent(in) :: m
         real(dp), intent(out) :: c(m)
         real(dp) :: f, ft, step
         real(dp) :: dvec(m), amat(m, m), g(m), trial(m)
         integer :: ii, jj, itr, ls
         logical :: improved
         do ii = 1, m
            dvec(ii) = lp%dot(n2*nspin, dstore_a(:, :, :, ii), fstore_a(:, :, :, m)) &
                       - lp%dot(n2*nspin, dstore_a(:, :, :, m), fstore_a(:, :, :, m))
            do jj = 1, m
               amat(ii, jj) = lp%dot(n2*nspin, dstore_a(:, :, :, ii), fstore_a(:, :, :, jj)) &
                              - lp%dot(n2*nspin, dstore_a(:, :, :, ii), fstore_a(:, :, :, m)) &
                              - lp%dot(n2*nspin, dstore_a(:, :, :, m), fstore_a(:, :, :, jj)) &
                              + lp%dot(n2*nspin, dstore_a(:, :, :, m), fstore_a(:, :, :, m))
            end do
         end do
         c = 0.0_dp
         c(m) = 1.0_dp
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
                  c = trial
                  f = ft
                  improved = .true.
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
         real(dp) :: css, theta, t
         real(dp) :: u(m)
         integer :: ii, rho
         u = c
         call sort_desc(m, u)
         css = 0.0_dp
         theta = 0.0_dp
         rho = 0
         do ii = 1, m
            css = css + u(ii)
            t = (css - 1.0_dp)/real(ii, dp)
            if (u(ii) - t > 0.0_dp) then
               rho = ii
               theta = t
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
                  t = u(ii)
                  u(ii) = u(jj)
                  u(jj) = t
               end if
            end do
         end do
      end subroutine sort_desc

   end subroutine trc_scf_run

   !
   ! The two-electron build, resident and entirely in Cartesian space:
   ! g_rks (nspin = 1, a full build), g_rks_incremental (nspin = 1, from the
   ! density DIFFERENCE, see trc_scf_run's own header) and g_uks (nspin = 2),
   ! dispatched on `nspin`/`incr_now` rather than three separately-named
   ! entry points, since the caller already has exactly one of them to ask
   ! for. Every array here is an explicit-shape dummy -- `nu_c`/`nlr_c`/
   ! `nin_c` carry the TRUE size of `dtot_c`/`jmat_c`, `klr_c`, and
   ! `dref_c`/`gacc_c`/`ddif_c` respectively (1 when this run has no use for
   ! them) -- so this compiles correctly once the caller's arrays are
   ! `trc_scf_state` components instead of plain locals: a `do concurrent`
   ! body here never names `st%anything`, only its own dummies.
   !
   subroutine build_g(b, eri, nao_c, nspin, nu_c, nlr_c, nin_c, incr_now, incr_on, &
                      exx, exx_lr, omega, incr_scale_run, &
                      d_c, gmat_c, dtot_c, jmat_c, w1_c, klr_c, gacc_c, dref_c, ddif_c)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: nao_c, nspin, nu_c, nlr_c, nin_c
      logical, intent(in) :: incr_now, incr_on
      real(dp), intent(in) :: exx, exx_lr, omega, incr_scale_run
      real(dp), intent(in) :: d_c(nao_c, nao_c, nspin)
      real(dp), intent(inout) :: gmat_c(nao_c, nao_c, nspin)
      real(dp), intent(inout) :: dtot_c(nu_c, nu_c), jmat_c(nu_c, nu_c)
      real(dp), intent(inout) :: w1_c(nao_c, nao_c), klr_c(nlr_c, nlr_c)
      real(dp), intent(inout) :: gacc_c(nin_c, nin_c), dref_c(nin_c, nin_c), ddif_c(nin_c, nin_c)
      integer :: i, j, s

      if (nspin == 1) then
         if (incr_now) then
            do concurrent(i=1:nao_c, j=1:nao_c)
               ddif_c(i, j) = d_c(i, j, 1) - dref_c(i, j)
            end do
            !
            ! A fixed, tighter threshold on the increment, not one tracking ddif.
            !
            ! The Schwarz/density screen drops a quartet whose bound qq*dd is below
            ! an absolute threshold: fine for a full build, where dd is O(1). On
            ! ddif, quartets sit at the edge of the test, and which ones changes
            ! each iteration because ddif is a difference of iterates with a
            ! changing shape. The dropped sets differ and the errors do not average
            ! down as the full build's (always the same negligible quartets) do. On
            ! w172/6-31G (450 atoms, 1950 AOs) RMS(D) reaches 1e-10 while the
            ! energy bounces by 1e-4 to 1e-6 Eh per iteration.
            !
            ! Scaling `thresh` by ||ddif||/||dref|| is backwards: it loosens as
            ! ddif shrinks and loses the speed (25.6 s against 24.4 s baseline).
            ! A fixed tighter factor, `incr_thresh_scale` (default 1e-4), keeps
            ! the error floor low enough that |dE| clears 1e-8 at every reset
            ! through iteration 20, while each increment stays a few times cheaper
            ! than a full build.
            !
            block
               real(dp) :: thresh0
               thresh0 = eri%thresh
               eri%thresh = thresh0*incr_scale_run
               call eri%fock_resident(b, ddif_c, gmat_c(:, :, 1), k_scale=exx)
               if (exx_lr /= 0.0_dp) then
                  call eri%fock_resident(b, ddif_c, klr_c, j_scale=0.0_dp, k_scale=exx_lr, omega=omega)
                  do concurrent(i=1:nao_c, j=1:nao_c)
                     gmat_c(i, j, 1) = gmat_c(i, j, 1) + klr_c(i, j)
                  end do
               end if
               eri%thresh = thresh0
            end block
            do concurrent(i=1:nao_c, j=1:nao_c)
               gacc_c(i, j) = gacc_c(i, j) + gmat_c(i, j, 1)
               gmat_c(i, j, 1) = gacc_c(i, j)
               dref_c(i, j) = d_c(i, j, 1)
            end do
         else
            call eri%fock_resident(b, d_c(:, :, 1), gmat_c(:, :, 1), k_scale=exx)
            if (exx_lr /= 0.0_dp) then
               call eri%fock_resident(b, d_c(:, :, 1), klr_c, j_scale=0.0_dp, k_scale=exx_lr, omega=omega)
               do concurrent(i=1:nao_c, j=1:nao_c)
                  gmat_c(i, j, 1) = gmat_c(i, j, 1) + klr_c(i, j)
               end do
            end if
            ! The reference for the next increment, kept only when there
            ! will be one (gacc_c and dref_c are 1 x 1 otherwise).
            if (incr_on) then
               do concurrent(i=1:nao_c, j=1:nao_c)
                  gacc_c(i, j) = gmat_c(i, j, 1)
                  dref_c(i, j) = d_c(i, j, 1)
               end do
            end if
         end if
      else
         do concurrent(i=1:nao_c, j=1:nao_c)
            dtot_c(i, j) = d_c(i, j, 1) + d_c(i, j, 2)
         end do
         call eri%fock_resident(b, dtot_c, jmat_c, k_scale=0.0_dp)
         do s = 1, 2
            if (exx /= 0.0_dp .or. exx_lr /= 0.0_dp) then
               do concurrent(i=1:nao_c, j=1:nao_c)
                  w1_c(i, j) = 2.0_dp*d_c(i, j, s)
               end do
               call eri%fock_resident(b, w1_c, gmat_c(:, :, s), j_scale=0.0_dp, k_scale=exx)
               if (exx_lr /= 0.0_dp) then
                  call eri%fock_resident(b, w1_c, klr_c, j_scale=0.0_dp, k_scale=exx_lr, omega=omega)
                  do concurrent(i=1:nao_c, j=1:nao_c)
                     gmat_c(i, j, s) = gmat_c(i, j, s) + klr_c(i, j)
                  end do
               end if
               do concurrent(i=1:nao_c, j=1:nao_c)
                  gmat_c(i, j, s) = gmat_c(i, j, s) + jmat_c(i, j)
               end do
            else
               do concurrent(i=1:nao_c, j=1:nao_c)
                  gmat_c(i, j, s) = jmat_c(i, j)
               end do
            end if
         end do
      end if
   end subroutine build_g

   subroutine assemble_fock(nao, nxc, npcm, nspin, dft, pcm_on, hcore, gmat, vxc, vpcm, fock)
      !! F = H + G (+ V_xc, + the continuum's operator), per spin. `vxc` is
      !! 1 x 1 (unused) when `.not. dft`, `vpcm` 1 x 1 when `.not. pcm_on`;
      !! `nxc`/`npcm` are their true sizes, so the explicit-shape dummies
      !! below match whichever the caller actually allocated.
      integer, intent(in) :: nao, nxc, npcm, nspin
      logical, intent(in) :: dft, pcm_on
      real(dp), intent(in) :: hcore(nao, nao)
      real(dp), intent(in) :: gmat(nao, nao, nspin), vxc(nxc, nxc, nspin)
      real(dp), intent(in) :: vpcm(npcm, npcm)
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      integer :: i, j, s
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
         if (pcm_on) then
            do concurrent(i=1:nao, j=1:nao)
               fock(i, j, s) = fock(i, j, s) + vpcm(i, j)
            end do
         end if
      end do
   end subroutine assemble_fock

   subroutine diis_store(nao, nspin, k, fock, companion, fstore, cstore)
      !! One push, shared by CDIIS (companion = the commutator `errv`) and
      !! ADIIS (companion = the density `dold`): slot k <- the current
      !! iterate in both histories at once.
      integer, intent(in) :: nao, nspin, k
      real(dp), intent(in) :: fock(nao, nao, nspin), companion(nao, nao, nspin)
      real(dp), intent(inout) :: fstore(nao, nao, nspin, *), cstore(nao, nao, nspin, *)
      integer :: i, j, s
      do concurrent(i=1:nao, j=1:nao, s=1:nspin)
         fstore(i, j, s, k) = fock(i, j, s)
         cstore(i, j, s, k) = companion(i, j, s)
      end do
   end subroutine diis_store

   subroutine diis_combine(nao, nspin, m, coeff, fstore, fock)
      !! fock := sum_k coeff(k) * fstore(:,:,:,k), the extrapolated/combined
      !! step both CDIIS (coeff = the solved weights) and ADIIS (coeff = the
      !! simplex coefficients) take once they have more than one history slot.
      integer, intent(in) :: nao, nspin, m
      real(dp), intent(in) :: coeff(m)
      real(dp), intent(in) :: fstore(nao, nao, nspin, *)
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      integer :: i, j, s, k
      do concurrent(i=1:nao, j=1:nao, s=1:nspin)
         fock(i, j, s) = 0.0_dp
      end do
      do k = 1, m
         do concurrent(i=1:nao, j=1:nao, s=1:nspin)
            fock(i, j, s) = fock(i, j, s) + coeff(k)*fstore(i, j, s, k)
         end do
      end do
   end subroutine diis_combine

   subroutine level_shift(lp, nao, nspin, mu, occ, smat, dmat, w1, w2, fock)
      !! F += mu (S - S D S / occ) on the virtual space while DIIS is off,
      !! per spin: see trc_scf_options_t%level_shift's own comment.
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: mu, occ
      real(dp), intent(in) :: smat(nao, nao), dmat(nao, nao, nspin)
      real(dp), intent(inout) :: w1(nao, nao), w2(nao, nao)
      real(dp), intent(inout) :: fock(nao, nao, nspin)
      integer :: i, j, s
      do s = 1, nspin
         call lp%gemm("N", "N", nao, nao, nao, 1.0_dp, smat, nao, dmat(:, :, s), nao, 0.0_dp, w1, nao)
         call lp%gemm("N", "N", nao, nao, nao, 1.0_dp, w1, nao, smat, nao, 0.0_dp, w2, nao)
         do concurrent(i=1:nao, j=1:nao)
            fock(i, j, s) = fock(i, j, s) + mu*(smat(i, j) - w2(i, j)/occ)
         end do
      end do
   end subroutine level_shift

   subroutine damp(nao, nspin, alpha, dold, dmat)
      !! The damped step while DIIS/A-DIIS are off: dmat := alpha*dmat +
      !! (1 - alpha)*dold, the weight `alpha` given to the new density.
      integer, intent(in) :: nao, nspin
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: dold(nao, nao, nspin)
      real(dp), intent(inout) :: dmat(nao, nao, nspin)
      integer :: i, j, s
      do concurrent(i=1:nao, j=1:nao, s=1:nspin)
         dmat(i, j, s) = alpha*dmat(i, j, s) + (1.0_dp - alpha)*dold(i, j, s)
      end do
   end subroutine damp

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
