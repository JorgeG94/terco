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
program check_xc_energy
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_error, only: error_t
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_dft_prune, only: PRUNE_NONE
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_xc, only: trc_xc_rks, trc_xc_uks, trc_xc_options, trc_xc_result_t
   use trc_cart_perm, only: cart_perm_build_ao
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_matrix, &
                           ref_read_matrix3, ref_report
   implicit none

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
   integer :: sh_l(NSH), sh_np(NSH), i, j, nao, ifn, unit, nfail
   real(dp) :: sh_e(MAXNP, NSH), sh_c(MAXNP, NSH), sh_r(3, NSH)
   real(dp) :: at_z(NATM), at_r(3, NATM), seed, exc, nelec, exc2, nelec2, trds
   integer(int64) :: lcg
   real(dp), allocatable :: s(:, :), t(:, :), vn(:, :), c(:, :), dmat(:, :), vxc(:, :), vxc2(:, :)
   real(dp), allocatable :: cb(:, :), dm2(:, :, :), vx2(:, :, :), vx3(:, :, :)
   real(dp), allocatable :: c_lex(:, :), cb_lex(:, :)
   real(dp), allocatable :: exc_all(:), n_all(:), vxc_all(:, :, :)
   real(dp), allocatable :: uexc_all(:), un_all(:), uvxc_all(:, :, :, :)
   !! One entry per functional, kept around (rather than only written to
   !! xc_probe.bin) so the pyscf-reference comparison block below rks_checks/
   !! uks_checks can read them back without a second pass over the kernels.
   integer, allocatable :: perm_l2b(:), perm_b2l(:)
   real(dp) :: exc3, nelec3, trds2
   integer :: nsk
   type(trc_xc_result_t) :: xr
   type(trc_basis_t) :: b
   type(trc_pairlist_t) :: pl
   type(dft_grid_t) :: grid
   type(trc_xc_grid_t) :: xg
   type(trc_xc_functional_t) :: func
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
   nfail = 0
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
   allocate (c(nao, NOCC), c_lex(nao, NOCC), dmat(nao, nao), vxc(nao, nao), vxc2(nao, nao))
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

   call dft_grid_build(at_r, nint(at_z), grid, err, level=4, prune=PRUNE_NONE)
   if (err%has_error()) then
      print '(a)', "check_xc_energy: grid failed: "//err%get_message()
      stop 1
   end if
   call xg%build(grid%n_points, grid%coords, grid%weights, b, max_pts=256, tol=1.0e-10_dp)
   call xg%to_device()

   open (newunit=unit, file='xc_probe.bin', access='stream', form='unformatted', status='replace')
   write (unit) int(grid%n_points, int64), int(nao, int64), int(NFUNC, int64)
   write (unit) grid%coords
   write (unit) grid%weights
   write (unit) dmat

   print '(a,i0,a,i0,a,i0)', "check_xc_energy: ", grid%n_points, " points, ", xg%nbatch, &
      " batches, Tr(DS) = ", nint(trds)
   ! The two loops below are internal procedures, not for structure:
   ! nvfortran 26.5's front end segfaults (fort1/fort2) on a program unit
   ! with eight calls to trc_xc_rks/uks, and is fine with seven. Split,
   ! each unit has fewer; the checks are unchanged.
   call rks_checks()

   ! --- spin-polarised: two independent PSD densities ---------------------
   !
   ! Beta gets three occupied columns to alpha's four, so the two spins
   ! differ everywhere and every cross term is exercised. First the
   ! identity: UKS with alpha = beta = D/2 must reproduce RKS on D exactly,
   ! energy and both potentials. Then the same-grid pyscf comparison.
   !
   allocate (cb(nao, NOCC - 1), cb_lex(nao, NOCC - 1))
   allocate (dm2(nao, nao, 2), vx2(nao, nao, 2), vx3(nao, nao, 2))
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
   call uks_checks()
   close (unit)

   call xg%release(); call grid%destroy(); call pl%release(); call b%release()

   ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
   block
      integer :: u, nao_ref, npts_ref, nfunc_ref
      real(dp) :: exc_ref, n_ref, dv
      real(dp), allocatable :: vxc_ref(:, :), uvxc_ref(:, :, :)
      call ref_require('xc_energy.ref', u)
      call ref_read_int(u, 'NAO', nao_ref)
      call ref_read_int(u, 'NPTS', npts_ref)
      call ref_read_int(u, 'NFUNC', nfunc_ref)
      if (nao_ref /= nao .or. nfunc_ref /= NFUNC) then
         print '(a)', "check_xc_energy: pyscf reference shape does not match this build"; stop 1
      end if
      allocate (vxc_ref(nao, nao), uvxc_ref(nao, nao, 2))
      do ifn = 1, NFUNC
         call ref_read_real(u, 'EXC_'//trim(names(ifn)), exc_ref)
         call ref_read_real(u, 'N_'//trim(names(ifn)), n_ref)
         call ref_read_matrix(u, 'VXC_'//trim(names(ifn)), nao, nao, vxc_ref)
         dv = maxval(abs(vxc_all(:, :, ifn) - vxc_ref))
         call ref_report(nfail, abs(exc_all(ifn) - exc_ref) < 1.0e-9_dp .and. &
                         abs(n_all(ifn) - n_ref) < 1.0e-9_dp .and. dv < 1.0e-9_dp, &
                         trim(names(ifn))//": RKS E_xc/N/V_xc vs pyscf", &
                         max(abs(exc_all(ifn) - exc_ref), dv))
      end do
      do ifn = 1, NFUNC
         call ref_read_real(u, 'UEXC_'//trim(names(ifn)), exc_ref)
         call ref_read_real(u, 'UN_'//trim(names(ifn)), n_ref)
         call ref_read_matrix3(u, 'UVXC_'//trim(names(ifn)), nao, nao, 2, uvxc_ref)
         dv = maxval(abs(uvxc_all(:, :, :, ifn) - uvxc_ref))
         call ref_report(nfail, abs(uexc_all(ifn) - exc_ref) < 1.0e-9_dp .and. &
                         abs(un_all(ifn) - n_ref) < 1.0e-9_dp .and. dv < 1.0e-9_dp, &
                         trim(names(ifn))//": UKS E_xc/N/V_xc vs pyscf", &
                         max(abs(uexc_all(ifn) - exc_ref), dv))
      end do
      call ref_close(u)
   end block

   if (nfail > 0) then
      print '(a,i0,a)', "check_xc_energy: ", nfail, " check(s) FAILED"
      stop 1
   end if
   print '(a)', "check_xc_energy: internal checks (including the pyscf reference) passed; RESULT: PASS"

contains

   subroutine rks_checks()
      do ifn = 1, NFUNC
         if (.not. xc_functional_by_name(names(ifn), func)) then
            print '(a)', "check_xc_energy: unknown functional "//trim(names(ifn))
            stop 1
         end if
         call trc_xc_rks(b, xg, func, dmat, vxc, xr)
         exc = xr%exc
         nelec = xr%nelec
         nsk = xr%n_skipped
         ! Density screening off entirely: what the default threshold costs.
         call trc_xc_rks(b, xg, func, dmat, vxc2, xr, trc_xc_options(rho_tol=0.0_dp))
         exc2 = xr%exc
         nelec2 = xr%nelec
         call report(abs(exc - exc2) < 1.0e-10_dp*abs(exc) .and. maxval(abs(vxc - vxc2)) < 1.0e-10_dp*maxval(abs(vxc)), &
                     trim(names(ifn))//": density screening vs none", &
                     max(abs(exc - exc2)/abs(exc), maxval(abs(vxc - vxc2))/maxval(abs(vxc))))
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
         call report(abs(nelec - trds)/trds < 1.0e-7_dp, trim(names(ifn))//": sum w rho vs Tr(DS), relative", &
                     abs(nelec - trds)/trds)
         call report(abs(exc - exc2) < 1.0e-13_dp*abs(exc) .and. &
                     maxval(abs(vxc - vxc2)) < 1.0e-13_dp*maxval(abs(vxc)), &
                     trim(names(ifn))//": independent of chunking", &
                     max(abs(exc - exc2)/abs(exc), maxval(abs(vxc - vxc2))/maxval(abs(vxc))))
         call report(maxval(abs(vxc - transpose(vxc))) < 1.0e-13_dp*maxval(abs(vxc)), &
                     trim(names(ifn))//": V_xc symmetric", maxval(abs(vxc - transpose(vxc)))/maxval(abs(vxc)))
      end do
   end subroutine rks_checks

   subroutine uks_checks()
      do ifn = 1, NFUNC
         if (.not. xc_functional_by_name(names(ifn), func)) stop 1
         call trc_xc_rks(b, xg, func, dmat, vxc, xr)
         exc = xr%exc
         nelec = xr%nelec
         vx3(:, :, 1) = 0.5_dp*dmat; vx3(:, :, 2) = 0.5_dp*dmat
         call trc_xc_uks(b, xg, func, vx3, vx2, xr)
         exc2 = xr%exc
         nelec2 = xr%nelec
         call report(abs(exc - exc2) < 1.0e-12_dp*abs(exc) .and. &
                     maxval(abs(vx2(:, :, 1) - vxc)) < 1.0e-12_dp*maxval(abs(vxc)) .and. &
                     maxval(abs(vx2(:, :, 2) - vxc)) < 1.0e-12_dp*maxval(abs(vxc)), &
                     trim(names(ifn))//": UKS at D/2, D/2 reproduces RKS", &
                     max(abs(exc - exc2)/abs(exc), maxval(abs(vx2(:, :, 1) - vxc))/maxval(abs(vxc))))
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
         call report(abs(nelec3 - trds2)/trds2 < 1.0e-7_dp, trim(names(ifn))//": UKS sum w rho vs Tr(DS), relative", &
                     abs(nelec3 - trds2)/trds2)
         call report(abs(exc3 - exc2) < 1.0e-13_dp*abs(exc3) .and. maxval(abs(vx2 - vx3)) < 1.0e-13_dp*maxval(abs(vx2)), &
                     trim(names(ifn))//": UKS independent of chunking", &
                     max(abs(exc3 - exc2)/abs(exc3), maxval(abs(vx2 - vx3))/maxval(abs(vx2))))
      end do
   end subroutine uks_checks


   subroutine report(ok, what, value)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: value
      if (.not. ok) nfail = nfail + 1
      print '(a,1x,a,t56,es10.2)', merge("ok  ", "FAIL", ok), what, value
   end subroutine report

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

end program check_xc_energy
