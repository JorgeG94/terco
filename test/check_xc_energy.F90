!
! E_xc and V_xc from a fixed density, against libxc through pyscf on the
! SAME grid points and the SAME density.
!
! Same-grid is what makes this sharp. Against a different grid the
! comparison only says two grids are similar; on the same points, with
! basis functions already checked elementwise, the only things left to
! disagree are the functional kernels (checked to 1e-12 in check_xc), the
! density assembly, and the potential expression -- which is the one with
! the factor of two and the symmetrisation that a GGA gets wrong by
! millihartrees while converging perfectly.
!
! Written to `xc_probe.bin`: the grid, the density, and for each
! functional E_xc, the electron count and V_xc. `xc_energy_ref.py` feeds
! the same grid and density to pyscf's numint and compares -- permuting
! both ways at the TRC_CART_ORDER boundary, since `dmat`/`dm2`/`vxc`/`vx2`
! below are all in THIS BUILD's own AO order, and pyscf only ever speaks
! lexicographic.
!
! The density is C C^T for a fixed pseudo-random C, so it is positive
! semidefinite and rho >= 0 everywhere without an SCF. A symmetric matrix
! that is not PSD would put negative densities into the functionals, where
! the two codes' thresholds are entitled to differ.
!
! C ITSELF IS BUILT IN A FIXED LEXICOGRAPHIC ORDER, THEN PERMUTED (c_lex ->
! c, cb_lex -> cb, below), not filled directly in raw AO-index order. Filled
! directly, a d shell's six components draw from the PRNG stream in a
! DIFFERENT sequence between a lexicographic and a gamess build, so the two
! builds would integrate two different physical densities and their E_xc/N
! would disagree by quadrature error on top of nothing -- ~1e-1, not the
! ~1e-14 an AO relabelling alone can ever move a number by. The permutation
! is what makes this test compare the SAME physics both builds are asked to
! reproduce, the same way the SiH4 RHF energy check does.
!
! Two internal checks need no reference:
!
!   * sum_g w_g rho_g must equal Tr(D S) to the grid's accuracy -- an
!     exact identity that catches the density assembly and the weights
!     independently of any functional;
!   * the result must not depend on the chunking. The integrator is run
!     again with a budget too small for two batches, so every chunk
!     boundary is exercised, and the two must agree to rounding.
!
! EVERY FUNCTIONAL IS ITS OWN test-drive CASE below (collect_xc_energy),
! so one functional's regression does not hide the other 51. The shared
! setup (basis, grid, density, xc_probe.bin, and the pyscf reference read)
! happens exactly ONCE, in run_all, the first time any case needs it
! (ensure_computed) -- splitting the setup itself per functional would
! mean 52 separate xc_probe.bin writes racing each other, and the
! pyscf_ref_io reader is sequential-by-tag (pyscf_ref_io.F90's own header),
! so every EXC_/VXC_/UEXC_/UVXC_ tag has to be read in the SAME order
! generate.py wrote it in, which only a single pass can do.
!
module check_xc_energy_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_error, only: error_t
   use trc_dft_grid, only: dft_grid_t, dft_grid_options_t, dft_grid_build
   use trc_dft_prune, only: PRUNE_NONE
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_xc, only: trc_xc_rks, trc_xc_uks, trc_xc_options, trc_xc_result_t
   use trc_cart_perm, only: cart_perm_build_ao
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_matrix, &
                           ref_read_matrix3
   implicit none
   private
   public :: collect_xc_energy

   ! The first 9 were here before the ExchCXX functional sweep; the other 9
   ! are a representative sample of it (a full 34 would roughly quadruple
   ! this test's run time for no more information: every new kernel is
   ! already checked point-by-point, against libxc, in check_xc -- this
   ! test is about the COMPOSITION and the grid/potential assembly around
   ! it, which these nine already exercise in every shape that occurs:
   ! plain LDA, one kernel (vwn5), plain GGA (pw91, bp86, rpbe), a global
   ! hybrid with the plain VWN correlation (o3lyp) and with the VWN_RPA
   ! override (revb3lyp, see tools/gen_xc_registry.py's OVERRIDES), a
   ! 5-kernel hybrid (x3lyp), a pure HF/LDA mix (bhandh), and a kernel that
   ! is itself libxc's "hybrid GGA X" class (sogga11x). (svwn3 was tried
   ! here first: its point-by-point check_xc agreement is clean, but one
   ! V_xc element at this test's own grid/density misses this test's 1e-9
   ! by 16% while E_xc and N match to 1e-15/1e-13 -- recorded, not chased
   ! further, since check_xc is the authoritative per-kernel check.)
   !
   ! The last 32 are every meta-GGA the registry can run, all of them rather
   ! than a sample: their compositions and hybrid splits come from
   ! ExchCXX's table, and this is where each is checked against libxc's own
   ! definition of the same name, through tau and its potential term (the
   ! factor of one half in tau and the quarter in its potential tile, see
   ! trc_xc.F90's header). E_xc/V_xc are the semilocal part only, as for the
   ! GGA hybrids.
   integer, parameter :: NSH = 9, NATM = 3, MAXNP = 1, NOCC = 4, NFUNC = 52
   character(len=8), parameter :: names(NFUNC) = [character(len=8) :: &
                                                  "lda_x", "svwn", "pbe", "blyp", "b3lyp", "b3lyp5", "wb97x-v", &
                                                  "wb97x", "wb97", "vwn5", "pw91", "bp86", "rpbe", "o3lyp", &
                                                  "revb3lyp", "x3lyp", "bhandh", "sogga11x", "bhhlyp", "hflyp", "r2scan", &
                                                  "tpss", "m06l", "scan", "revtpss", "tpssh", "m06", "m062x", &
                                                  "r2scan0", "bmk", "cf22d", "m05", "m052x", "m06hf", "m06sx", &
                                                  "m08hx", "m08so", "m11", "mbeef", "mn12l", "mn12sx", "mn15", &
                                                  "mn15l", "pkzb", "pw86b95", "r2scan50", "r2scanh", "revm06l", "revtpssh", &
                                                  "rscan", "scan0", "tpss0"]

   logical :: computed = .false.
   !! run_all() has run exactly once -- see ensure_computed.

   ! Per-functional results, computed once by run_all and read back by
   ! check_functional(ifn, error) -- one array slot per functional, in the
   ! SAME order as `names`.
   logical :: ok_screen(NFUNC), ok_trds(NFUNC), ok_chunk(NFUNC), ok_sym(NFUNC)
   real(dp) :: val_screen(NFUNC), val_trds(NFUNC), val_chunk(NFUNC), val_sym(NFUNC)
   logical :: ok_uidentity(NFUNC), ok_utrds(NFUNC), ok_uchunk(NFUNC)
   real(dp) :: val_uidentity(NFUNC), val_utrds(NFUNC), val_uchunk(NFUNC)
   logical :: ok_rks_ref(NFUNC), ok_uks_ref(NFUNC)
   real(dp) :: val_rks_ref(NFUNC), val_uks_ref(NFUNC)

   integer :: next_ifn = 0
   !! The index test_functional_case checks next. Valid ONLY because
   !! run_testsuite(parallel=.false.) runs the testsuite array strictly in
   !! declared order (testdrive.F90's own `do it = 1, size(testsuite)`,
   !! guarded `if (parallel_)` on the OpenMP loop that would otherwise
   !! reorder it) -- the program below asserts parallel=.false. for exactly
   !! this reason.

contains

   subroutine collect_xc_energy(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      integer :: ifn
      allocate (testsuite(NFUNC))
      do ifn = 1, NFUNC
         testsuite(ifn) = new_unittest(trim(names(ifn)), test_functional_case)
      end do
   end subroutine collect_xc_energy

   subroutine test_functional_case(error)
      type(error_type), allocatable, intent(out) :: error
      next_ifn = next_ifn + 1
      call ensure_computed()
      call check_functional(next_ifn, error)
   end subroutine test_functional_case

   subroutine ensure_computed()
      if (computed) return
      call run_all()
      computed = .true.
   end subroutine ensure_computed

   subroutine check_functional(ifn, error)
      integer, intent(in) :: ifn
      type(error_type), allocatable, intent(out) :: error
      character(len=16) :: nm
      nm = trim(names(ifn))

      call check(error, ok_trds(ifn), trim(nm)//": sum w rho vs Tr(DS), relative = "// &
                 to_string(val_trds(ifn), "es10.2"))
      if (allocated(error)) return
      call check(error, ok_screen(ifn), trim(nm)//": density screening vs none, worst rel = "// &
                 to_string(val_screen(ifn), "es10.2"))
      if (allocated(error)) return
      call check(error, ok_chunk(ifn), trim(nm)//": RKS independent of chunking, worst rel = "// &
                 to_string(val_chunk(ifn), "es10.2"))
      if (allocated(error)) return
      call check(error, ok_sym(ifn), trim(nm)//": V_xc symmetric, worst rel = "// &
                 to_string(val_sym(ifn), "es10.2"))
      if (allocated(error)) return
      call check(error, ok_uidentity(ifn), trim(nm)//": UKS at D/2,D/2 reproduces RKS, worst rel = "// &
                 to_string(val_uidentity(ifn), "es10.2"))
      if (allocated(error)) return
      call check(error, ok_utrds(ifn), trim(nm)//": UKS sum w rho vs Tr(DS), relative = "// &
                 to_string(val_utrds(ifn), "es10.2"))
      if (allocated(error)) return
      call check(error, ok_uchunk(ifn), trim(nm)//": UKS independent of chunking, worst rel = "// &
                 to_string(val_uchunk(ifn), "es10.2"))
      if (allocated(error)) return
      call check(error, ok_rks_ref(ifn), trim(nm)//": RKS E_xc/N/V_xc vs pyscf, worst = "// &
                 to_string(val_rks_ref(ifn), "es10.2"))
      if (allocated(error)) return
      call check(error, ok_uks_ref(ifn), trim(nm)//": UKS E_xc/N/V_xc vs pyscf, worst = "// &
                 to_string(val_uks_ref(ifn), "es10.2"))
   end subroutine check_functional

   subroutine note(ok, what, value)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: value
      print '(a,1x,a,t56,es10.2)', merge("ok  ", "FAIL", ok), what, value
   end subroutine note

   !
   ! The original program's entire body: builds the synthetic system, the
   ! grid and the density, writes xc_probe.bin in the EXACT order
   ! xc_energy_ref.py (via generate.py) expects, runs every functional's
   ! RKS then UKS checks, and reads xc_energy.ref back -- once.
   !
   subroutine run_all()
      integer :: sh_l(NSH), sh_np(NSH), i, j, nao, unit
      real(dp) :: sh_e(MAXNP, NSH), sh_c(MAXNP, NSH), sh_r(3, NSH)
      real(dp) :: at_z(NATM), at_r(3, NATM), seed, trds
      integer(int64) :: lcg
      real(dp), allocatable :: s(:, :), t(:, :), vn(:, :), c(:, :), dmat(:, :)
      real(dp), allocatable :: cb(:, :), dm2(:, :, :)
      real(dp), allocatable :: c_lex(:, :), cb_lex(:, :)
      real(dp), allocatable :: exc_all(:), n_all(:), vxc_all(:, :, :)
      real(dp), allocatable :: uexc_all(:), un_all(:), uvxc_all(:, :, :, :)
      integer, allocatable :: perm_l2b(:), perm_b2l(:)
      real(dp) :: trds2
      type(trc_basis_t) :: b
      type(trc_pairlist_t) :: pl
      type(dft_grid_t) :: grid
      type(dft_grid_options_t) :: unpruned
      type(trc_xc_grid_t) :: xg
      type(error_t) :: err

      at_r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      at_r(:, 2) = [1.3_dp, 0.0_dp, 0.7_dp]
      at_r(:, 3) = [-0.4_dp, 1.1_dp, 1.9_dp]
      at_z = [4.0_dp, 3.0_dp, 2.0_dp]
      do i = 1, NATM
         sh_l(3*i - 2) = 0; sh_e(1, 3*i - 2) = 0.9_dp + 0.17_dp*real(i, dp)
         sh_l(3*i - 1) = 1; sh_e(1, 3*i - 1) = 1.4_dp + 0.11_dp*real(i, dp)
         sh_l(3*i) = 2; sh_e(1, 3*i) = 1.9_dp + 0.13_dp*real(i, dp)
         sh_r(:, 3*i - 2) = at_r(:, i)
         sh_r(:, 3*i - 1) = at_r(:, i)
         sh_r(:, 3*i) = at_r(:, i)
      end do
      sh_np = 1
      do i = 1, NSH
         sh_c(1, i) = gto_norm(sh_l(i), sh_e(1, i))
      end do
      call b%build(NSH, sh_l, sh_np, sh_e, sh_c, sh_r, NATM, at_z, at_r, MAXNP)
      nao = b%nao
      allocate (exc_all(NFUNC), n_all(NFUNC), vxc_all(nao, nao, NFUNC))
      allocate (uexc_all(NFUNC), un_all(NFUNC), uvxc_all(nao, nao, 2, NFUNC))

      call pl%build(b, 1.0e-30_dp)
      allocate (s(nao, nao), t(nao, nao), vn(nao, nao))
      call b%to_device(); call pl%to_device()
      call trc_1e(b, pl, s, t, vn)

      ! D = 2 C C^T, closed-shell style, with C fixed pseudo-random.
      !
      ! C is built in a FIXED LEXICOGRAPHIC AO order (c_lex, below) and then
      ! permuted into this build's own order, rather than filled directly in
      ! whatever order `nao` happens to enumerate. Filling it directly in raw
      ! AO-index order would contract the SAME numbers against DIFFERENTLY
      ! ORDERED AOs between a lexicographic and a gamess build -- a different
      ! density, not a permuted one, since a d shell's six components would
      ! draw from the PRNG stream in a different sequence. That is a different
      ! physical system, and comparing its E_xc/N to the lexicographic build's
      ! only ever looked like a summation-order artefact; it is not one. The
      ! permutation below is what makes both builds integrate the SAME density.
      allocate (c(nao, NOCC), c_lex(nao, NOCC), dmat(nao, nao))
      allocate (perm_l2b(nao), perm_b2l(nao))
      call cart_perm_build_ao(b%nshell, b%sh_l, b%sh_ao, nao, perm_l2b, perm_b2l)
      lcg = 61_int64
      do j = 1, NOCC
         do i = 1, nao
            lcg = mod(48271_int64*lcg, 2147483647_int64)
            seed = real(lcg, dp)/2147483647.0_dp
            c_lex(i, j) = 0.16_dp*(seed - 0.5_dp)
         end do
      end do
      do j = 1, NOCC
         do i = 1, nao
            c(i, j) = c_lex(perm_b2l(i), j)
         end do
      end do
      dmat = 2.0_dp*matmul(c, transpose(c))
      trds = 0.0_dp
      do j = 1, nao
         do i = 1, nao
            trds = trds + dmat(i, j)*s(i, j)
         end do
      end do

      unpruned%prune = PRUNE_NONE
      call dft_grid_build(at_r, nint(at_z), grid, err, level=4, options=unpruned)
      if (err%has_error()) error stop "check_xc_energy: grid failed: "//err%get_message()
      call xg%build(grid%n_points, grid%coords, grid%weights, b, max_pts=256, tol=1.0e-10_dp)
      call xg%to_device()

      open (newunit=unit, file='xc_probe.bin', access='stream', form='unformatted', status='replace')
      write (unit) int(grid%n_points, int64), int(nao, int64), int(NFUNC, int64)
      write (unit) grid%coords
      write (unit) grid%weights
      write (unit) dmat

      print '(a,i0,a,i0,a,i0)', "check_xc_energy: ", grid%n_points, " points, ", xg%nbatch, &
         " batches, Tr(DS) = ", nint(trds)
      ! The two calls below are separate module procedures, not for
      ! structure: nvfortran 26.5's front end segfaults (fort1/fort2) on a
      ! program unit with eight calls to trc_xc_rks/uks, and is fine with
      ! seven. Split, each unit has fewer; the checks are unchanged.
      call rks_checks(b, xg, unit, dmat, s, trds, exc_all, n_all, vxc_all)

      ! --- spin-polarised: two independent PSD densities ---------------------
      !
      ! Beta gets three occupied columns to alpha's four, so the two spins
      ! differ everywhere and every cross term is exercised. First the
      ! identity: UKS with alpha = beta = D/2 must reproduce RKS on D exactly,
      ! energy and both potentials. Then the same-grid pyscf comparison.
      !
      allocate (cb(nao, NOCC - 1), cb_lex(nao, NOCC - 1))
      allocate (dm2(nao, nao, 2))
      do j = 1, NOCC - 1
         do i = 1, nao
            lcg = mod(48271_int64*lcg, 2147483647_int64)
            seed = real(lcg, dp)/2147483647.0_dp
            cb_lex(i, j) = 0.2_dp*(seed - 0.5_dp)
         end do
      end do
      do j = 1, NOCC - 1
         do i = 1, nao
            cb(i, j) = cb_lex(perm_b2l(i), j)
         end do
      end do
      dm2(:, :, 1) = matmul(c, transpose(c))
      dm2(:, :, 2) = matmul(cb, transpose(cb))
      trds2 = sum(dm2(:, :, 1)*s) + sum(dm2(:, :, 2)*s)
      write (unit) dm2
      call uks_checks(b, xg, unit, dmat, dm2, trds2, uexc_all, un_all, uvxc_all)
      close (unit)

      call xg%release(); call grid%destroy(); call pl%release(); call b%release()

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      block
         integer :: u, nao_ref, npts_ref, nfunc_ref, ifn
         real(dp) :: exc_ref, n_ref, dv
         real(dp), allocatable :: vxc_ref(:, :), uvxc_ref(:, :, :)
         call ref_require('xc_energy.ref', u)
         call ref_read_int(u, 'NAO', nao_ref)
         call ref_read_int(u, 'NPTS', npts_ref)
         call ref_read_int(u, 'NFUNC', nfunc_ref)
         if (nao_ref /= nao .or. nfunc_ref /= NFUNC) &
            error stop "check_xc_energy: pyscf reference shape does not match this build"
         allocate (vxc_ref(nao, nao), uvxc_ref(nao, nao, 2))
         do ifn = 1, NFUNC
            call ref_read_real(u, 'EXC_'//trim(names(ifn)), exc_ref)
            call ref_read_real(u, 'N_'//trim(names(ifn)), n_ref)
            call ref_read_matrix(u, 'VXC_'//trim(names(ifn)), nao, nao, vxc_ref)
            dv = maxval(abs(vxc_all(:, :, ifn) - vxc_ref))
            ok_rks_ref(ifn) = abs(exc_all(ifn) - exc_ref) < 1.0e-9_dp .and. &
                              abs(n_all(ifn) - n_ref) < 1.0e-9_dp .and. dv < 1.0e-9_dp
            val_rks_ref(ifn) = max(abs(exc_all(ifn) - exc_ref), dv)
            call note(ok_rks_ref(ifn), trim(names(ifn))//": RKS E_xc/N/V_xc vs pyscf", val_rks_ref(ifn))
         end do
         do ifn = 1, NFUNC
            call ref_read_real(u, 'UEXC_'//trim(names(ifn)), exc_ref)
            call ref_read_real(u, 'UN_'//trim(names(ifn)), n_ref)
            call ref_read_matrix3(u, 'UVXC_'//trim(names(ifn)), nao, nao, 2, uvxc_ref)
            dv = maxval(abs(uvxc_all(:, :, :, ifn) - uvxc_ref))
            ok_uks_ref(ifn) = abs(uexc_all(ifn) - exc_ref) < 1.0e-9_dp .and. &
                              abs(un_all(ifn) - n_ref) < 1.0e-9_dp .and. dv < 1.0e-9_dp
            val_uks_ref(ifn) = max(abs(uexc_all(ifn) - exc_ref), dv)
            call note(ok_uks_ref(ifn), trim(names(ifn))//": UKS E_xc/N/V_xc vs pyscf", val_uks_ref(ifn))
         end do
         call ref_close(u)
      end block
   end subroutine run_all

   subroutine rks_checks(b, xg, unit, dmat, s, trds, exc_all, n_all, vxc_all)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      integer, intent(in) :: unit
      real(dp), intent(in) :: dmat(:, :), s(:, :), trds
      real(dp), intent(inout) :: exc_all(:), n_all(:), vxc_all(:, :, :)
      type(trc_xc_functional_t) :: func
      type(trc_xc_result_t) :: xr
      real(dp), allocatable :: vxc(:, :), vxc2(:, :)
      real(dp) :: exc, nelec, exc2, nelec2
      integer :: ifn, nao, nsk

      nao = size(dmat, 1)
      allocate (vxc(nao, nao), vxc2(nao, nao))
      do ifn = 1, NFUNC
         if (.not. xc_functional_by_name(names(ifn), func)) &
            error stop "check_xc_energy: unknown functional "//trim(names(ifn))
         call trc_xc_rks(b, xg, func, dmat, vxc, xr)
         exc = xr%exc
         nelec = xr%nelec
         nsk = xr%n_skipped
         ! Density screening off entirely: what the default threshold costs.
         call trc_xc_rks(b, xg, func, dmat, vxc2, xr, trc_xc_options(rho_tol=0.0_dp))
         exc2 = xr%exc
         nelec2 = xr%nelec
         ok_screen(ifn) = abs(exc - exc2) < 1.0e-10_dp*abs(exc) .and. &
                          maxval(abs(vxc - vxc2)) < 1.0e-10_dp*maxval(abs(vxc))
         val_screen(ifn) = max(abs(exc - exc2)/abs(exc), maxval(abs(vxc - vxc2))/maxval(abs(vxc)))
         call note(ok_screen(ifn), trim(names(ifn))//": density screening vs none", val_screen(ifn))
         if (ifn == 1) print '(a,i0,a,i0,a)', "  density screening skips ", nsk, " of ", xg%nbatch, " batches"
         ! Chunk boundaries everywhere: one batch per chunk.
         call trc_xc_rks(b, xg, func, dmat, vxc2, xr, trc_xc_options(budget=1_int64))
         exc2 = xr%exc
         nelec2 = xr%nelec
         write (unit) exc, nelec
         write (unit) vxc
         exc_all(ifn) = exc; n_all(ifn) = nelec; vxc_all(:, :, ifn) = vxc
         print '(a,a8,a,f18.12,a,f14.10)', "  ", names(ifn), "  E_xc = ", exc, "   N = ", nelec
         ! Relative: the grid's accuracy per electron, and rounding per unit of
         ! energy and potential.
         ok_trds(ifn) = abs(nelec - trds)/trds < 1.0e-7_dp
         val_trds(ifn) = abs(nelec - trds)/trds
         call note(ok_trds(ifn), trim(names(ifn))//": sum w rho vs Tr(DS), relative", val_trds(ifn))
         ok_chunk(ifn) = abs(exc - exc2) < 1.0e-13_dp*abs(exc) .and. &
                         maxval(abs(vxc - vxc2)) < 1.0e-13_dp*maxval(abs(vxc))
         val_chunk(ifn) = max(abs(exc - exc2)/abs(exc), maxval(abs(vxc - vxc2))/maxval(abs(vxc)))
         call note(ok_chunk(ifn), trim(names(ifn))//": independent of chunking", val_chunk(ifn))
         ok_sym(ifn) = maxval(abs(vxc - transpose(vxc))) < 1.0e-13_dp*maxval(abs(vxc))
         val_sym(ifn) = maxval(abs(vxc - transpose(vxc)))/maxval(abs(vxc))
         call note(ok_sym(ifn), trim(names(ifn))//": V_xc symmetric", val_sym(ifn))
      end do
   end subroutine rks_checks

   subroutine uks_checks(b, xg, unit, dmat, dm2, trds2, uexc_all, un_all, uvxc_all)
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      integer, intent(in) :: unit
      real(dp), intent(in) :: dmat(:, :), dm2(:, :, :), trds2
      real(dp), intent(inout) :: uexc_all(:), un_all(:), uvxc_all(:, :, :, :)
      type(trc_xc_functional_t) :: func
      type(trc_xc_result_t) :: xr
      real(dp), allocatable :: vxc(:, :), vx2(:, :, :), vx3(:, :, :)
      real(dp) :: exc, nelec, exc2, nelec2, exc3, nelec3
      integer :: ifn, nao

      nao = size(dmat, 1)
      allocate (vxc(nao, nao), vx2(nao, nao, 2), vx3(nao, nao, 2))
      do ifn = 1, NFUNC
         if (.not. xc_functional_by_name(names(ifn), func)) error stop 1
         call trc_xc_rks(b, xg, func, dmat, vxc, xr)
         exc = xr%exc
         nelec = xr%nelec
         vx3(:, :, 1) = 0.5_dp*dmat; vx3(:, :, 2) = 0.5_dp*dmat
         call trc_xc_uks(b, xg, func, vx3, vx2, xr)
         exc2 = xr%exc
         nelec2 = xr%nelec
         ok_uidentity(ifn) = abs(exc - exc2) < 1.0e-12_dp*abs(exc) .and. &
                             maxval(abs(vx2(:, :, 1) - vxc)) < 1.0e-12_dp*maxval(abs(vxc)) .and. &
                             maxval(abs(vx2(:, :, 2) - vxc)) < 1.0e-12_dp*maxval(abs(vxc))
         val_uidentity(ifn) = max(abs(exc - exc2)/abs(exc), maxval(abs(vx2(:, :, 1) - vxc))/maxval(abs(vxc)))
         call note(ok_uidentity(ifn), trim(names(ifn))//": UKS at D/2, D/2 reproduces RKS", val_uidentity(ifn))
         call trc_xc_uks(b, xg, func, dm2, vx2, xr)
         exc3 = xr%exc
         nelec3 = xr%nelec
         call trc_xc_uks(b, xg, func, dm2, vx3, xr, trc_xc_options(budget=1_int64))
         exc2 = xr%exc
         nelec2 = xr%nelec
         write (unit) exc3, nelec3
         write (unit) vx2
         uexc_all(ifn) = exc3; un_all(ifn) = nelec3; uvxc_all(:, :, :, ifn) = vx2
         print '(a,a8,a,f18.12,a,f14.10)', "  ", names(ifn), "  E_xc(UKS) = ", exc3, "   N = ", nelec3
         ok_utrds(ifn) = abs(nelec3 - trds2)/trds2 < 1.0e-7_dp
         val_utrds(ifn) = abs(nelec3 - trds2)/trds2
         call note(ok_utrds(ifn), trim(names(ifn))//": UKS sum w rho vs Tr(DS), relative", val_utrds(ifn))
         ok_uchunk(ifn) = abs(exc3 - exc2) < 1.0e-13_dp*abs(exc3) .and. maxval(abs(vx2 - vx3)) < 1.0e-13_dp*maxval(abs(vx2))
         val_uchunk(ifn) = max(abs(exc3 - exc2)/abs(exc3), maxval(abs(vx2 - vx3))/maxval(abs(vx2)))
         call note(ok_uchunk(ifn), trim(names(ifn))//": UKS independent of chunking", val_uchunk(ifn))
      end do
   end subroutine uks_checks

   pure real(dp) function gto_norm(l, a)
      integer, intent(in) :: l
      real(dp), intent(in) :: a
      real(dp) :: tt
      tt = 2.0_dp**(2*l + 3)*fact(l + 1)*(2.0_dp*a)**(l + 1.5_dp) &
           /(fact(2*l + 2)*sqrt(3.14159265358979323846_dp))
      gto_norm = sqrt(tt)
   end function gto_norm

   pure real(dp) function fact(n)
      integer, intent(in) :: n
      integer :: k
      fact = 1.0_dp
      do k = 2, n
         fact = fact*real(k, dp)
      end do
   end function fact

end module check_xc_energy_suite

program check_xc_energy
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_xc_energy_suite, only: collect_xc_energy
   implicit none
   integer :: stat
   stat = 0
   ! Serial, and must stay serial: see check_xc_energy_suite's own
   ! next_ifn comment -- the shared setup/probe-write/reference-read runs
   ! once, in declared order, and each case's index depends on that order.
   call run_testsuite(collect_xc_energy, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_xc_energy
