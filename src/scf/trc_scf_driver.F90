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
   use trc_error, only: error_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_settings, only: trc_scf_options_t, trc_scf_result_t, scf_settings_t
   use trc_scf_context, only: scf_ctx_t, bcast_flat
   use trc_scf_state, only: scf_state_t, scf_timers_t
   use trc_scf_fock, only: build_g, build_vxc, build_pcm, scf_energy, assemble_fock, &
                           commutator_error, to_cart, to_sph
   use trc_scf_df, only: df_jk_t
   use trc_scf_diis, only: cdiis_t, adiis_t, level_shift, damp
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
   subroutine trc_scf_run(b, nalpha, nbeta, opts, res, dguess, comm, hguess, la, dfjk)
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
      type(df_jk_t), intent(inout), optional, target :: dfjk
      !! The density-fitted J/K object for `opts%df` (`--df`), already fitted
      !! by the caller -- required when `opts%df` is true (trc_scf_context's
      !! `setup` refuses otherwise), ignored when it is false. Owned by the
      !! caller across more than one run (trc_capi.F90's context handle), the
      !! same reason `la` is a caller-owned, optional argument rather than
      !! something this routine builds for itself.

      integer :: nao, nao_c, s, it, i, j, n2
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
      type(cdiis_t) :: cdiis
      type(adiis_t) :: adiis
      !! ADIIS's history (fock, density) is kept apart from CDIIS's own (fock,
      !! commutator): ADIIS runs exactly in the window CDIIS has not opened
      !! yet, so the two never need the same entries, and keeping them
      !! separate means `opts%adiis = .false.` leaves the CDIIS path
      !! untouched, to the bit.
      logical :: adiis_active
      type(trc_linalg_t), pointer :: lp
      type(error_t) :: err
      real(dp) :: epcm
      logical :: diis_on
      real(dp) :: ctl(2)
      real(dp) :: etot, eold, drms, e1, e2, exc, ngrid, errmax
      integer :: nfrac
      !! This iteration's share of the two totals, so the table says where
      !! the time went while it is going rather than only at the end. A
      !! Fock build is not the same cost twice: the SAD guess is atomic, so
      !! density screening throws nearly all of iteration 1 away, and the
      !! per-iteration cost then climbs to a plateau. Comparing two runs by
      !! total/iterations is therefore only valid at equal iteration counts.

      ! Everything `opts` resolves to before the loop starts: the functional
      ! lookup and its rejections, the spin rule, the by-spin damping/DIIS
      ! defaults, and the env-hook overrides (trc_scf_settings). `cfg` is
      ! immutable from here on, and the rest of this routine reads it
      ! directly.
      ok = cfg%resolve(opts, nalpha, nbeta, msg=msg, comm=comm)
      if (.not. ok) then
         res%message = msg
         return
      end if
      res%nspin = cfg%nspin
      ! The cavity, K and its LU are built and solved whole on every rank
      ! (replicated: cheap next to the integrals); pcm%operator splits the
      ! two integral passes over the ranks itself, given `comm`.
      epcm = 0.0_dp

      ! --- everything once per geometry: pairs, 1e, S/H fold, ERI, linalg,
      ! orthogonaliser, and (when asked for) the PCM cavity and DFT grid ---
      call ctx%setup(b, opts, cfg, err, comm=comm, la=la, dfjk=dfjk)
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
      if (cfg%pcm_on) res%pcm_npoints = ctx%pcm%surface%n_points
      if (cfg%dft) then
         res%grid_r = ctx%grid%coords
         res%grid_w = ctx%grid%weights
         if (cfg%nlc) then
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
      nxc = merge(nao, 1, cfg%dft)
      nu_c = merge(nao_c, 1, cfg%nspin == 2)
      nlr_c = merge(nao_c, 1, cfg%exx_lr /= 0.0_dp)
      nin_c = merge(nao_c, 1, cfg%incr_on)
      npcm = merge(nao, 1, cfg%pcm_on)
      call st%alloc(ctx, cfg)
      call cdiis%alloc(nao, cfg%nspin, opts%ndiis)
      allocate (res%dmat(nao, nao, cfg%nspin), res%cmo(nao, nao, cfg%nspin), res%eps(nao, cfg%nspin))
      ! Only when asked for: a second ndiis-deep history the size of
      ! cdiis's own, which matters on a system the size of w172, where
      ! one nao^2*ndiis block is already hundreds of megabytes.
      if (opts%adiis) call adiis%alloc(nao, cfg%nspin, opts%ndiis)

      ! --- the guess, on the host, then everything goes up -------------------
      if (present(dguess)) then
         res%dmat = dguess
      else if (present(hguess)) then
         do s = 1, cfg%nspin
            call guess_from(hguess, cfg%nocc(s), res%dmat(:, :, s))
         end do
      else if (opts%guess == "core") then
         do s = 1, cfg%nspin
            call guess_from(ctx%hcore, cfg%nocc(s), res%dmat(:, :, s))
         end do
      else
         ! GWH: F_ij = 1.75/2 S_ij (H_ii + H_jj), diagonal H_ii.
         do j = 1, nao
            do i = 1, nao
               st%w1(i, j) = 0.875_dp*ctx%smat(i, j)*(ctx%hcore(i, i) + ctx%hcore(j, j))
            end do
            st%w1(j, j) = ctx%hcore(j, j)
         end do
         do s = 1, cfg%nspin
            call guess_from(st%w1, cfg%nocc(s), res%dmat(:, :, s))
         end do
      end if
      res%cmo = 0.0_dp
      res%eps = 0.0_dp

      call st%to_device(ctx, res)
      call cdiis%to_device()
      if (opts%adiis) call adiis%to_device()

      eold = 0.0_dp
      ! Read after the loop (res%e_one/e_two) even when it never runs
      ! (max_iter=0): the loop body is the only other place that sets them.
      e1 = 0.0_dp
      e2 = 0.0_dp
      ! Read by the first iteration's `incr_now` before any commutator exists
      ! (Fortran's .and. does not short-circuit), so give it a value.
      errmax = huge(1.0_dp)
      diis_on = .false.
      exc = 0.0_dp
      ngrid = 0.0_dp
      ! What the SCF is actually running on. Two callers can hand over the
      ! same molecule and basis name and different shells -- one had a
      ! reader that kept dead primitives -- and that only shows here.
      if (cfg%talk) print "(a,i0,a,i0,a,i0,a,i0)", "   basis: ", nao, " functions, ", b%nshell, &
         " shells, ", sum(b%sh_np), " primitives, max ", maxval(b%sh_np)
      ! Say so when increments are on: they are opt-in and change the timings.
      ! An `i` after an iteration's line marks a Fock matrix built from the
      ! density difference.
      if (cfg%talk .and. cfg%incr_on) print "(a,es8.1,a,i0,a)", &
         "   incremental Fock: from |FDS-SDF| < ", cfg%incr_start_run, &
         ", full build every ", max(1, cfg%incr_reset_run), " iterations (marked i)"
      if (cfg%talk) then
         if (cfg%dft) then
            print "(a)", "   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s)     xc(s)"
         else
            print "(a)", "   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s)"
         end if
      end if
      do it = 1, opts%max_iter
         call tmr%fock%start()
         incr_now = .false.
         ! --- D_cart := T D_sph T^T, the one transform every Cartesian
         ! quantity below (the Fock build, the XC density) is built from.
         ! Not taken at all when `sph` is false -- see this module's own
         ! header on why that keeps the Cartesian path byte-for-byte
         ! unchanged rather than merely numerically equivalent.
         do s = 1, cfg%nspin
            call to_cart(b, sph, nao, nao_c, res%dmat(:, :, s), st%d_c(:, :, s))
         end do
         ! --- two-electron part, resident, entirely in Cartesian space -------
         !
         ! INCREMENTAL. G(D) = G(D_ref) + G(D - D_ref), exactly, because G is
         ! linear in D. The point is the screening: the kernel screens on the
         ! density it is handed, and the difference is small and getting
         ! smaller while the density is neither. See `build_g`'s own header.
         !
         incr_now = cfg%incr_on .and. it > 1 &
                    .and. errmax < cfg%incr_start_run &
                    .and. mod(it - 1, max(1, cfg%incr_reset_run)) /= 0
         call build_g(b, ctx%eri, nao_c, cfg%nspin, nu_c, nlr_c, nin_c, incr_now, cfg%incr_on, &
                      cfg%exx, cfg%exx_lr, cfg%omega, cfg%incr_scale_run, &
                      st%d_c, st%gmat_c, st%dtot_c, st%jmat_c, st%w1_c, st%klr_c, &
                      st%gacc_c, st%dref_c, st%ddif_c, cfg%df, ctx%dfj, lp)
         ! --- G_sph := T^T G_cart T, back up to the space the rest of the
         ! iteration runs in.
         do s = 1, cfg%nspin
            call to_sph(b, sph, nao, nao_c, st%gmat_c(:, :, s), st%gmat(:, :, s))
         end do
         tmr%tf_it = tmr%fock%get_elapsed_time()
         tmr%t_fock = tmr%t_fock + tmr%tf_it
         call tmr%xc%start()
         ! --- exchange-correlation: finds its (Cartesian) density present ----
         if (cfg%dft) then
            call build_vxc(b, ctx%xg, ctx%xg_nlc, cfg%func, cfg%nlc, cfg%nspin, nao_c, nao, sph, opts%rho_tol, comm, lp, &
                           st%d_c, st%vxc_c, st%vxc, exc, ngrid, tmr%tx_pts, tmr%tx_prs)
         end if
         tmr%tx_it = tmr%xc%get_elapsed_time()
         tmr%t_xc = tmr%t_xc + tmr%tx_it
         call tmr%rest%start()
         ! --- the continuum, at the TOTAL density (pyscf's choice for UKS
         ! too: one set of charges, one operator added to both spins) -----
         if (cfg%pcm_on) then
            call tmr%pcm%start()
            call build_pcm(b, ctx%pcm, ctx%pl, cfg%nspin, nao_c, nu_c, nao, sph, &
                           st%d_c, st%dtot_c, st%vpcm_c, st%vpcm, epcm, comm)
            tmr%t_pcm = tmr%t_pcm + tmr%pcm%get_elapsed_time()
            ! the PCM step is its own column, not part of the rest
            call tmr%rest%start()
         end if
         ! --- energy at this density -----------------------------------------
         call scf_energy(lp, nao, cfg%nspin, res%dmat, ctx%hcore, st%gmat, res%e_nuc, exc, epcm, e1, e2, etot)
         ! --- Fock and its commutator, per spin ------------------------------
         call assemble_fock(nao, nxc, npcm, cfg%nspin, cfg%dft, cfg%pcm_on, ctx%hcore, st%gmat, st%vxc, st%vpcm, st%fock)
         call commutator_error(lp, nao, cfg%nspin, st%fock, res%dmat, ctx%smat, ctx%x, st%w1, st%w2, st%errv)
         ! `dold` must hold `res%dmat` from BEFORE DIIS/diagonalisation touch
         ! it -- see trc_scf_state's own header.
         call st%save_dold(nao, cfg%nspin, res%dmat)
         ! --- DIIS, once the commutator is small enough to extrapolate from --
         errmax = sqrt(lp%dot_3d(st%errv, st%errv))
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
         if (.not. diis_on .and. it > cfg%ndamp) diis_on = errmax < cfg%diis_start
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
         adiis_active = opts%adiis .and. .not. diis_on .and. it > cfg%ndamp
         ! The new density is rank 0's alone: extrapolation, eigensolve and density
         ! run on one device and the density is broadcast, since dsyevd need not
         ! return the same vectors on two devices. Orbitals and energies come over
         ! once, after the loop.
         if (cfg%root) then
            if (adiis_active) then
               call adiis%push(lp, nao, cfg%nspin, st%fock, st%dold)
               call adiis%extrapolate(nao, cfg%nspin, st%fock)
            end if
            ! --- level shift while DIIS and A/EDIIS are off: F += mu (S - S D S / occ)
            if (.not. diis_on .and. .not. adiis_active .and. opts%level_shift > 0.0_dp) then
               call level_shift(lp, nao, cfg%nspin, opts%level_shift, cfg%occ, ctx%smat, res%dmat, st%w1, st%w2, st%fock)
            end if
            if (diis_on) then
               call cdiis%push(lp, nao, cfg%nspin, st%fock, st%errv)
               call cdiis%extrapolate(nao, cfg%nspin, st%fock, ok)
            end if
            ! --- diagonalise, per spin, and the new density ---------------------
            do s = 1, cfg%nspin
               call orbitals_and_density(st%fock(:, :, s), cfg%nocc(s), opts%frac_occ, st%w1, st%w2, &
                                         res%eps(:, s), res%cmo(:, :, s), res%dmat(:, :, s), st%ofrac)
            end do
            if (.not. diis_on .and. .not. adiis_active) call damp(nao, cfg%nspin, opts%damp, st%dold, res%dmat)
         end if
         ! Every rank must hold the same density to the bit. The Fock sum is
         ! broadcast already; the XC potential is accumulated with atomics
         ! whose order differs per rank, 1e-13 in the energy on water, and a
         ! difference that small compounds through a diverging iteration.
         if (present(comm)) then
            if (comm%size() > 1) then
               !$acc update self(res%dmat)
               call bcast_flat(comm, res%dmat, n2*cfg%nspin)
               !$acc update device(res%dmat)
            end if
         end if
         call st%update_errv(nao, cfg%nspin, res%dmat)
         drms = sqrt(lp%dot_3d(st%errv, st%errv))/real(nao, dp)
         tmr%t_rest = tmr%t_rest + tmr%rest%get_elapsed_time()
         if (cfg%talk) then
            if (cfg%dft) then
               print "(i5,f22.12,es14.4,es14.4,es14.4,f10.3,f10.3,a)", &
                  it, etot, etot - eold, drms, errmax, tmr%tf_it, tmr%tx_it, merge(" i", "  ", incr_now)
            else
               print "(i5,f22.12,es14.4,es14.4,es14.4,f10.3,a)", &
                  it, etot, etot - eold, drms, errmax, tmr%tf_it, merge(" i", "  ", incr_now)
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
            call bcast_flat(comm, res%cmo, n2*cfg%nspin)
            call bcast_flat(comm, res%eps, nao*cfg%nspin)
         end if
      end if
      call st%release_device(ctx, res)
      call cdiis%release_device()
      if (opts%adiis) call adiis%release_device()
      res%dmat = st%dold
      res%energy = etot
      res%e_one = e1
      res%e_two = e2
      res%e_xc = exc
      res%nelec_grid = ngrid
      res%e_pcm = epcm
      res%pcm_qtotal = ctx%pcm%q_total
      if (.not. res%converged) res%message = "trc_scf: not converged"
      if (cfg%talk) then
         call report_timing()
         if (cfg%pcm_on) print "(a,f20.12,a,f10.6,a,f10.3,a,f10.4,a)", "   pcm: E_diel = ", epcm, &
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
         ! F' = X^T F X
         call lp%gemm_mat(f, ctx%x, t)
         call lp%gemm_mat(ctx%x, t, fp, transa="T")
         call lp%syev(fp, eps_s)
         call lp%gemm_mat(ctx%x, fp, cmo_s)
         if (frac) then
            ! D = sum_i o_i c_i c_i^T through a column-scaled copy of C.
            !$acc update self(eps_s)
            if (cfg%nspin == 1) then
               call frac_occupations(nao, eps_s, opts%nelec_frac, 2.0_dp, ofrac, nfrac)
            else
               call frac_occupations(nao, eps_s, real(nocc_s, dp), 1.0_dp, ofrac, nfrac)
            end if
            !$acc update device(ofrac)
            do concurrent(i=1:nao, j=1:nao)
               t(i, j) = cmo_s(i, j)*sqrt(ofrac(j))
            end do
            ! D = occ * C(:,1:k) C(:,1:k)^T with occ = 1, k = nfrac, on the
            ! column-scaled copy -- the same GEMM lp%density does on t with
            ! those arguments (A(1:nao,1:nfrac) is the same block either way).
            call lp%density(t, nfrac, 1.0_dp, d_s)
         else
            ! k <= 0 is lp%density's own zeroing branch (do concurrent),
            ! replacing zero_dev.
            call lp%density(cmo_s, nocc_s, cfg%occ, d_s)
         end if
      end subroutine orbitals_and_density

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

   end subroutine trc_scf_run

end module trc_scf_driver
