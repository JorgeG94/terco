!
! GAMESS's per-component AO normalisation, at the C API boundary.
!
! terco (like libcint) normalises a Cartesian shell per ANGULAR MOMENTUM: a
! d shell's six components share one constant, so <XY|XY> is 1/3, not 1.
! GAMESS normalises every component INDIVIDUALLY to 1. trc_set_ao_norm
! switches trc_capi's own AO-indexed entries between the two; see
! src/trc_ao_norm.F90's header for the closed form and src/trc_capi.F90's
! "AO NORMALISATION" section for exactly which entries convert.
!
! FORCED CARTESIAN. cc-pVDZ's d shells are tagged "gto_spherical" in BSE, so
! since the spherical AO space was added this basis defaults to running
! spherical (trc_api.F90's `b%spherical`). This check's entire subject is
! the per-COMPONENT GAMESS convention, which is defined on CARTESIAN
! components only and refuses to combine with a spherical basis at all
! (trc_set_ao_norm's own refusal, trc_capi.F90) -- so every basis load here
! forces Cartesian (`force_cart=.true.`/`trc_set_force_cartesian`), which is
! what keeps this check testing the same thing it always has rather than
! either failing outright on the refusal or silently changing subject.
!
! SIX CHECKS.
!
!   (e) THE FACTOR ITSELF, l = 0..4, against test/reference_ao_norm.dat --
!       sqrt(gamess_norm_sq(lx,ly,lz)) from scripts/cartorder.py, generated
!       once and recorded (the same reason check_binfock reads reference.dat
!       rather than calling libfint from Fortran: no Python at ctest time).
!       Indexed by the EXPONENT TRIPLE, not by within-shell offset, so the
!       one file checks whichever TRC_CART_ORDER this build was made with --
!       order and normalisation are independent on purpose (both headers
!       above say why), and this test would rather prove that than assume
!       it. Through trc_ao_norm's `gamess_norm` directly -- the pure ratio,
!       not `ao_norm_vector`'s full boundary conversion, which ALSO folds
!       in terco's own libcint-inherited `common_fac_sp` scale for s and p
!       (see trc_ao_norm.F90's header; that is a real but SEPARATE effect
!       this check is not the place to re-derive). Nothing here exercises
!       the actual kernels, because ONE_E_LMAX in this generated kernel set
!       is 2 and the formula does not need an integral to be checked
!       against itself -- see (a) below for where the kernels DO get
!       exercised, up to the d this build's one-electron path supports.
!
!   (a) THE OVERLAP DIAGONAL, water/cc-pVDZ (real d functions), through
!       trc_context_basis + trc_compute_1e: exactly 1 for every s, p and
!       "pure power" d component (XX/YY/ZZ), under EITHER mode -- terco's
!       own raw <s|s>/<p_x|p_x> already equals GAMESS's (and pyscf's)
!       unit convention, see trc_ao_norm.F90's header for why. Only the
!       three "mixed" d components (XY/XZ/YZ) differ by mode: 1/3 under
!       TRC_AO_NORM_SHELL (terco shares one normalisation constant across
!       a shell's components), 1 under TRC_AO_NORM_UNIT (GAMESS normalises
!       every component individually) -- found by looking up each AO's own
!       (px,py,pz) via cart_pow, the same way (e) does, so this check does
!       not depend on knowing TRC_CART_ORDER's slot assignment by hand
!       either.
!
!       AN EARLIER VERSION OF THIS CHECK asserted 1/(4*pi) and 3/(4*pi) for
!       SHELL-mode s and p. That number was the trc_basis_json bug this
!       basis (loaded through trc_basis_from_json, below) used to have --
!       folding `common_fac_sp` twice into an s or p shell's coefficients --
!       reverse-engineered into the test as if it were the convention,
!       without checking it against the OTHER basis constructor
!       (test/trc_test_basis.F90's `build_631g`) or against pyscf, both of
!       which give 1 directly. It stayed green because
!       trc_ao_norm.F90's own `ao_norm_vector` carried a SECOND
!       `common_fac_sp` divide that canceled it for (b)/(c)/(d)/(f), which
!       test self-consistency and cannot see an overall basis-wide scale
!       error; only this check's hard-coded numbers could, and did, the
!       moment the JSON reader's fold was fixed on its own.
!
!   (b) THE SCF ENERGY, SHELL and UNIT, must agree to ~1e-12 Eh: the whole
!       point of a change of AO NORMALISATION is that it is not a change of
!       physics, only of what a getter hands back.
!
!   (c) Tr(D S) = N_electron, in EACH mode's own D and S -- the physical
!       invariant a contravariant-covariant pair is built to preserve
!       (Tr(N^-1 D N^-1 . N S N) = Tr(D S) exactly, for any diagonal N), so
!       this is a self-consistency check on ao_scale_cov/ao_scale_contra
!       agreeing about which direction is which.
!
!   (d) ROUND TRIP: the density trc_density hands back in a mode is fed
!       straight into trc_set_guess in the SAME mode; the SCF must already
!       be converged (few iterations, same energy). A mismatched
!       conversion (covariant where it should be contravariant, or a sign
!       error in the exponent) would still pass (b) and (c) -- both are
!       symmetric in a way that a transposed mistake cannot break -- but
!       would feed trc_run_scf a density that is not actually the converged
!       one, and it would take real iterations to get back.
!
!   (f) THE FOCK BUILD (`check_fock`), the entry GAMESS's own driver calls
!       for J/K (run_Fockbuild's role in hpcchem_driver.src): ONE physical
!       density, D_shell from a converged SHELL-mode SCF, fed to trc_fock
!       in SHELL mode directly and, EXPLICITLY converted here (D_unit =
!       N^-1 D_shell N^-1, not from a second, independently-converged SCF
!       -- see the subroutine's own comment for why that distinction is
!       the whole point), to trc_fock in UNIT mode. Checks F_unit = N
!       F_shell N to ~1e-12 relative, and that Tr(D F) -- the two-electron
!       energy expression -- agrees between the two mode's own (D, F)
!       pairs, which holds by the SAME congruence identity (c) uses,
!       independent of any SCF.
!
program check_ao_norm
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char, c_loc
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_batch, only: cart_pow
   use trc_ao_norm, only: gamess_norm, ao_norm_vector, ao_scale_cov, ao_scale_contra
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_ao_norm, trc_get_ao_norm, trc_set_force_cartesian, &
                               trc_context_basis, trc_context_eri, &
                               trc_pairs_create, trc_pairs_destroy, trc_compute_1e, trc_fock, &
                               trc_set_method, trc_set_convergence, trc_set_guess, trc_run_scf, &
                               trc_energy, trc_density, trc_nao, trc_converged, trc_iterations, &
                               TRC_OK, TRC_GUESS_GIVEN, &
                               TRC_AO_NORM_SHELL, TRC_AO_NORM_UNIT
   implicit none

   type(trc_basis_t) :: b
   type(error_t) :: err
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   integer :: natm, dev
   real(c_double) :: e_shell, e_unit
   logical :: ok

   dev = trc_bind_device(0)
   ok = .true.

   call check_factors(ok)

   call read_xyz('water.xyz', natm, zint, at_r)
   call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, b, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', "check_ao_norm: "//err%get_message()
      stop 1
   end if

   call run_mode(TRC_AO_NORM_SHELL, "shell", e_shell, ok)
   call run_mode(TRC_AO_NORM_UNIT, "unit", e_unit, ok)

   print '(a,f20.12,a,f20.12,a,es9.2)', "check_ao_norm: E_scf shell ", e_shell, "  unit ", e_unit, &
      "  diff ", e_shell - e_unit
   call expect(ok, abs(e_shell - e_unit) < 1.0e-12_dp, "(b) SCF energy agrees between shell and unit")

   call check_fock(ok)

   call b%release()
   if (.not. ok) then
      print '(a)', "check_ao_norm: FAIL"
      stop 1
   end if
   print '(a)', "check_ao_norm: PASS"

contains

   subroutine expect(ok, cond, what)
      logical, intent(inout) :: ok
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "check_ao_norm: FAILED: "//what
         ok = .false.
      end if
   end subroutine expect

   subroutine check_factors(ok)
      !! (e) The closed form, l = 0..4, every component, against
      !! test/reference_ao_norm.dat -- indexed by (lx,ly,lz), not by
      !! within-shell offset, so this build's own TRC_CART_ORDER (via
      !! `cart_pow`) decides which offset is checked against which line, and
      !! the one file serves every build.
      logical, intent(inout) :: ok
      integer, parameter :: NREF = 35
      integer :: ref_p(3, NREF)
      real(dp) :: ref_n(NREF)
      integer :: u, i, l, ic, ncart, px, py, pz, j
      real(dp) :: got
      character(len=256) :: line
      logical :: found

      open (newunit=u, file='reference_ao_norm.dat', status='old', action='read')
      read (u, '(a)') line  ! header comment
      do i = 1, NREF
         read (u, *) ref_p(1, i), ref_p(2, i), ref_p(3, i), ref_n(i)
      end do
      close (u)

      do l = 0, 4
         ncart = (l + 1)*(l + 2)/2
         do ic = 0, ncart - 1
            call cart_pow(l, ic, px, py, pz)
            found = .false.
            do j = 1, NREF
               if (ref_p(1, j) == px .and. ref_p(2, j) == py .and. ref_p(3, j) == pz) then
                  found = .true.; exit
               end if
            end do
            call expect(ok, found, "(e) reference has l="//itoa(l)//" component")
            if (.not. found) cycle
            got = gamess_norm(px, py, pz)
            call expect(ok, abs(got - ref_n(j)) < 1.0e-12_dp, &
                        "(e) l="//itoa(l)//" ic="//itoa(ic)//" ("//itoa(px)//itoa(py)//itoa(pz)// &
                        ") factor matches scripts/cartorder.py")
         end do
      end do
   end subroutine check_factors

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=12) :: buf
      write (buf, '(i0)') i
      s = trim(buf)
   end function itoa

   subroutine run_mode(mode, label, e_out, ok)
      !! Runs (a), the SCF, (c) and (d) in one AO-normalisation mode; hands
      !! back the SCF energy for (b)'s comparison across modes.
      integer(c_int), intent(in) :: mode
      character(len=*), intent(in) :: label
      real(c_double), intent(out) :: e_out
      logical, intent(inout) :: ok
      type(c_ptr) :: h, basis, pairs
      integer(c_int) :: rc, nao, got_mode, flag
      real(c_double), allocatable, target :: smat(:), tmat(:), vmat(:), dmat(:)
      real(dp) :: trds
      integer :: i, ish, ic, ncart, base, px, py, pz
      integer(c_int) :: niter2

      rc = trc_create(h); call expect(ok, rc == TRC_OK, label//": create")
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(ok, rc == TRC_OK, label//": set_molecule")
      rc = trc_set_ao_norm(h, mode); call expect(ok, rc == TRC_OK, label//": set_ao_norm")
      rc = trc_get_ao_norm(h, got_mode)
      call expect(ok, rc == TRC_OK .and. got_mode == mode, label//": get_ao_norm reads back what was set")
      rc = trc_set_force_cartesian(h, 1_c_int)
      call expect(ok, rc == TRC_OK, label//": set_force_cartesian")
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      call expect(ok, rc == TRC_OK, label//": set_basis_json")
      rc = trc_nao(h, nao); call expect(ok, rc == TRC_OK .and. nao == b%nao, label//": nao")

      ! (a) THE OVERLAP DIAGONAL, through trc_context_basis + trc_compute_1e.
      rc = trc_context_basis(h, basis); call expect(ok, rc == TRC_OK, label//": context_basis")
      rc = trc_pairs_create(basis, 1.0e-14_c_double, pairs)
      call expect(ok, rc == TRC_OK, label//": pairs_create")
      allocate (smat(nao*nao), tmat(nao*nao), vmat(nao*nao))
      rc = trc_compute_1e(basis, pairs, smat, tmat, vmat)
      call expect(ok, rc == TRC_OK, label//": compute_1e")
      do ish = 1, b%nshell
         base = b%sh_ao(ish) - 1
         ncart = (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2
         do ic = 0, ncart - 1
            call cart_pow(b%sh_l(ish), ic, px, py, pz)
            i = base + ic + 1
            associate (diag => smat((i - 1)*nao + i))
               if (mode == TRC_AO_NORM_UNIT) then
                  call expect(ok, abs(diag - 1.0_c_double) < 1.0e-12_c_double, &
                              label//": (a) unit-mode overlap diagonal is 1, AO "//itoa(i))
               else if (b%sh_l(ish) == 0) then
                  ! s: shell mode's raw diagonal is 1, same as GAMESS and
                  ! pyscf -- see trc_ao_norm.F90's header for why.
                  call expect(ok, abs(diag - 1.0_c_double) < 1.0e-10_c_double, &
                              label//": (a) shell-mode s diagonal is 1, AO "//itoa(i))
               else if (b%sh_l(ish) == 1) then
                  ! p: same reason, same answer.
                  call expect(ok, abs(diag - 1.0_c_double) < 1.0e-10_c_double, &
                              label//": (a) shell-mode p diagonal is 1, AO "//itoa(i))
               else if (px == 2 .or. py == 2 .or. pz == 2) then
                  ! a "pure power" d/f/g component (XX-type): GAMESS's and
                  ! terco's conventions coincide there, diagonal is 1.
                  call expect(ok, abs(diag - 1.0_c_double) < 1.0e-10_c_double, &
                              label//": (a) shell-mode pure-power diagonal is 1, AO "//itoa(i))
               else if (b%sh_l(ish) == 2) then
                  ! a mixed d component (XY/XZ/YZ): shell normalisation
                  ! leaves <chi|chi> = 1/gamess_norm_sq = 1/3.
                  call expect(ok, abs(diag - 1.0_c_double/3.0_c_double) < 1.0e-10_c_double, &
                              label//": (a) shell-mode d XY/XZ/YZ diagonal is 1/3, AO "//itoa(i))
               end if
            end associate
         end do
      end do

      ! --- the SCF, SAD guess, default convergence ---------------------
      rc = trc_set_method(h, ''//c_null_char, 3_c_int); call expect(ok, rc == TRC_OK, label//": set_method")
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      call expect(ok, rc == TRC_OK, label//": set_convergence")
      rc = trc_run_scf(h); call expect(ok, rc == TRC_OK, label//": run_scf")
      rc = trc_converged(h, flag); call expect(ok, flag == 1, label//": converged")
      rc = trc_energy(h, e_out)

      ! (c) Tr(D S) = N_electron, THIS mode's own D and S.
      allocate (dmat(nao*nao))
      rc = trc_density(h, dmat); call expect(ok, rc == TRC_OK, label//": density")
      trds = 0.0_dp
      do i = 1, nao*nao
         trds = trds + real(dmat(i), dp)*real(smat(i), dp)
      end do
      ! S above is from BEFORE the SCF touched anything, so it is the same
      ! S under this same mode -- trc_compute_1e depends only on the basis
      ! and its ao_norm_kind, neither of which the SCF changes.
      call expect(ok, abs(trds - real(sum(zint), dp)) < 1.0e-8_dp, &
                  label//": (c) Tr(D S) = N_electron")

      ! (d) ROUND TRIP: this mode's own converged density, fed back as the
      ! guess in this SAME mode, must already be the answer.
      rc = trc_set_guess(h, TRC_GUESS_GIVEN, c_loc(dmat(1)), 1_c_int)
      call expect(ok, rc == TRC_OK, label//": (d) set_guess given")
      rc = trc_run_scf(h); call expect(ok, rc == TRC_OK, label//": (d) second run_scf")
      rc = trc_iterations(h, niter2)
      call expect(ok, niter2 <= 3, label//": (d) round trip converges in a few iterations, got "//itoa(int(niter2)))
      block
         real(c_double) :: e2
         rc = trc_energy(h, e2)
         call expect(ok, abs(e2 - e_out) < 1.0e-10_c_double, label//": (d) round-trip energy matches")
      end block

      rc = trc_pairs_destroy(pairs)
      rc = trc_destroy(h)
      deallocate (smat, tmat, vmat, dmat)
   end subroutine run_mode

   subroutine check_fock(ok)
      !! (f) THE FOCK BUILD, the entry GAMESS's own driver calls for J/K
      !! (run_Fockbuild's role in hpcchem_driver.src) -- see the module
      !! header's own paragraph on this check.
      !!
      !! ONE physical density: D_shell, from a converged SHELL-mode SCF, fed
      !! to trc_fock in SHELL mode directly, and D_unit = N^-1 D_shell N^-1,
      !! built HERE by explicit conversion -- NOT from a second, independently
      !! converged SCF in UNIT mode. Two independent SCF runs agree only to
      !! the SCF's own convergence threshold (~1e-8 to ~1e-10 here), which
      !! would cap this check's precision at that, not the ~1e-12 the
      !! transform itself is exact to; an explicit conversion of the SAME
      !! array sidesteps that entirely. h_unit therefore never runs an SCF at
      !! all -- it exists only to lend trc_context_basis/trc_context_eri
      !! handles for trc_fock's UNIT-mode path to run on.
      logical, intent(inout) :: ok
      type(c_ptr) :: h_shell, h_unit, basis_shell, basis_unit, eri_shell, eri_unit
      integer(c_int) :: rc, nao
      real(c_double), allocatable, target :: d_shell(:), f_shell(:), d_unit(:), f_unit(:)
      real(dp), allocatable :: nrm(:), tmp(:, :)
      real(dp) :: trdf_shell, trdf_unit, maxdiff, maxf

      rc = trc_create(h_shell); call expect(ok, rc == TRC_OK, "fock: shell create")
      rc = trc_set_molecule(h_shell, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                            0_c_int, 1_c_int)
      rc = trc_set_ao_norm(h_shell, TRC_AO_NORM_SHELL)
      rc = trc_set_force_cartesian(h_shell, 1_c_int)
      rc = trc_set_basis_json(h_shell, 'basis_sets/cc-pvdz.json'//c_null_char)
      call expect(ok, rc == TRC_OK, "fock: shell set_basis_json")
      rc = trc_nao(h_shell, nao)
      rc = trc_set_method(h_shell, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h_shell, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(h_shell); call expect(ok, rc == TRC_OK, "fock: shell run_scf")

      allocate (d_shell(nao*nao), f_shell(nao*nao))
      rc = trc_density(h_shell, d_shell); call expect(ok, rc == TRC_OK, "fock: shell density")
      rc = trc_context_basis(h_shell, basis_shell); call expect(ok, rc == TRC_OK, "fock: shell context_basis")
      rc = trc_context_eri(h_shell, eri_shell); call expect(ok, rc == TRC_OK, "fock: shell context_eri")
      rc = trc_fock(eri_shell, basis_shell, d_shell, f_shell, 1.0_c_double, 1.0_c_double, 1_c_int)
      call expect(ok, rc == TRC_OK, "fock: shell fock")

      ! D_unit = N^-1 D_shell N^-1: CONTRAVARIANT, the direction trc_density
      ! itself uses (see that entry's own comment) -- a density is built by
      ! CONTRACTING two AO indices, not by pairing two AO functions the way
      ! an integral does, so this is `ao_scale_contra`, not `ao_scale_cov`.
      allocate (nrm(nao), tmp(nao, nao))
      call ao_norm_vector(b%nshell, b%sh_l, b%sh_ao, nao, TRC_AO_NORM_UNIT, nrm)
      tmp = reshape(real(d_shell, dp), [nao, nao])
      call ao_scale_contra(nrm, tmp)
      allocate (d_unit(nao*nao), f_unit(nao*nao))

      rc = trc_create(h_unit); call expect(ok, rc == TRC_OK, "fock: unit create")
      rc = trc_set_molecule(h_unit, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                            0_c_int, 1_c_int)
      rc = trc_set_ao_norm(h_unit, TRC_AO_NORM_UNIT)
      rc = trc_set_force_cartesian(h_unit, 1_c_int)
      rc = trc_set_basis_json(h_unit, 'basis_sets/cc-pvdz.json'//c_null_char)
      call expect(ok, rc == TRC_OK, "fock: unit set_basis_json")
      rc = trc_context_basis(h_unit, basis_unit); call expect(ok, rc == TRC_OK, "fock: unit context_basis")
      rc = trc_context_eri(h_unit, eri_unit); call expect(ok, rc == TRC_OK, "fock: unit context_eri")
      d_unit = real(reshape(tmp, [nao*nao]), c_double)
      rc = trc_fock(eri_unit, basis_unit, d_unit, f_unit, 1.0_c_double, 1.0_c_double, 1_c_int)
      call expect(ok, rc == TRC_OK, "fock: unit fock")

      ! F_unit must equal N F_shell N.
      tmp = reshape(real(f_shell, dp), [nao, nao])
      call ao_scale_cov(nrm, tmp)
      maxf = maxval(abs(tmp))
      maxdiff = maxval(abs(tmp - reshape(real(f_unit, dp), [nao, nao])))
      print '(a,es9.2,a,es9.2)', "check_ao_norm: (f) Fock  max|F_unit - N F_shell N| ", maxdiff, &
         "  relative to max|F| ", maxdiff/maxf
      call expect(ok, maxdiff < 1.0e-12_dp*maxf, "(f) F_unit = N F_shell N to ~1e-12 relative")

      ! Tr(D F), the two-electron energy expression, agrees between the
      ! two modes' own (D, F) pairs by the same congruence identity (c)
      ! uses: Tr(N^-1 D N^-1 . N F N) = Tr(D F) exactly.
      trdf_shell = sum(real(d_shell, dp)*real(f_shell, dp))
      trdf_unit = sum(real(d_unit, dp)*real(f_unit, dp))
      print '(a,es9.2)', "check_ao_norm: (f) Tr(D F)  |shell - unit| ", abs(trdf_shell - trdf_unit)
      call expect(ok, abs(trdf_shell - trdf_unit) < 1.0e-10_dp*max(abs(trdf_shell), 1.0_dp), &
                  "(f) Tr(D F) agrees between shell and unit")

      rc = trc_destroy(h_shell); rc = trc_destroy(h_unit)
      deallocate (d_shell, f_shell, d_unit, f_unit, nrm, tmp)
   end subroutine check_fock

end program check_ao_norm
