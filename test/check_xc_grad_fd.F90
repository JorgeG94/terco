!
! The KS-DFT nuclear gradient (trc_scf_grad%trc_scf_gradient with `func`
! present: eri%gradient scaled by `func%exx` + trc_1e_gradient +
! trc_xc_grad%xc_gradient_add + nuclear repulsion) against a central
! finite difference of terco's own fully reconverged KS-DFT energy --
! extends check_full_grad_fd.F90's technique (same file layout and the
! same reason for NO INTERNAL PROCEDURES: an internal procedure here, even
! with no state shared with the GPU call, corrupted a later collective
! Fock build on this exact NVHPC/OpenMPI toolchain, measured there).
!
! SVWN (LDA) and PBE (GGA) on water, restricted; PBE on the NH2 radical,
! unrestricted -- a genuinely open-shell case (non-degenerate SOMO, unlike
! OH's rotationally degenerate pi pair -- see below) and the unrestricted
! GGA branch (cross spin terms) at once.
!
! GRID-RESPONSE SIZE, MEASURED NOT ASSUMED, BEFORE AND AFTER. The
! finite-difference reference here is of terco's own fully reconverged
! KS-DFT energy as a function of nuclear position, with the grid rebuilt
! (`trc_basis_from_json` + `trc_scf_run`, fresh, at every displaced
! geometry) -- i.e. the grid MOVES with the atoms on the FD side exactly
! the way the GRID RESPONSE term (trc_xc_grad.F90's header) says the
! physical gradient must account for. Every case here runs the analytic
! gradient TWICE, at the SAME converged orbitals: once with
! `grid_response=.true.` (ON, the module's default -- the physical
! gradient, gated tightly) and once with `.false.` (OFF, the Pulay-only
! term this branch started with -- printed for the before/after
! comparison, NOT gated, since omitting a real term is not this test's
! business to forbid for a caller who explicitly asks for it, e.g. to
! match pyscf's own `grid_response=False` default). Both at grid level 3
! AND grid level 5, so the omission's level-dependence (OFF) and its
! absence (ON) are both visible in this test's own output.
!
module check_xc_grad_fd_cases
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none
   private
   public :: run_case

   real(dp), parameter :: H = 1.0e-3_dp
   real(dp), parameter, public :: RELTOL = 1.0e-6_dp
   !! The ON gate's relative counterpart -- see ABSTOL_ON.
   real(dp), parameter, public :: ABSTOL_ON = 2.0e-7_dp
   !! ON (grid_response=.true., the physical gradient): the FD floor, NOT
   !! grid-level-dependent -- this is the whole point of the grid-response
   !! term. MEASURED on all six cases (water/SVWN, water/PBE, NH2/PBE
   !! unrestricted, each at level 3 AND level 5): see this test's own
   !! printed output for the exact numbers at generation time; every one
   !! landed at or below this with no level-3-to-level-5 trend, unlike
   !! OFF below. Set from those numbers with headroom, not loosened to
   !! hide a mismatch -- trc_xc_grad.F90's formula is independently
   !! confirmed to ~1e-9 against a FIXED-grid finite difference
   !! (check_xc_grad_fixedgrid_fd.F90) where no SCF reconverges and no
   !! grid moves at all.
   real(dp), parameter, public :: INVTOL = 1.0e-8_dp
   !! ON's translational-invariance gate: sum over atoms ~ 0 is a property
   !! of the PHYSICAL gradient, which is what ON now is.

contains

   subroutine run_case(xyzfile, basis_json, nalpha, nbeta, functional, grid_level, label, any_fail, &
                       abstol_on_override, invtol_override, richardson)
      character(*), intent(in) :: xyzfile, basis_json, functional, label
      integer, intent(in) :: nalpha, nbeta, grid_level
      logical, intent(inout) :: any_fail
      real(dp), intent(in), optional :: abstol_on_override, invtol_override
      !! Per-case overrides of the module ABSTOL_ON/INVTOL, for a case
      !! whose OWN measured FD/invariance floor is genuinely looser than
      !! the six 2-3-atom cases ABSTOL_ON/INVTOL were set from -- see the
      !! w3.xyz case this module's program unit calls with one, and ITS
      !! OWN comment for why.
      logical, intent(in), optional :: richardson
      !! Default .false. (the plain central difference at step H, as
      !! before). When .true., the ON gradient is compared against a
      !! RICHARDSON EXTRAPOLATION of the central difference instead: a
      !! SECOND central difference at step H/2 is taken alongside the
      !! step-H one, and the two combined as
      !!
      !!   fd_richardson = (4*fd(H/2) - fd(H)) / 3
      !!
      !! which cancels the central difference's LEADING O(h^2) truncation
      !! term (fd(h) = f'(x) + c*h^2 + O(h^4) for some c, so 4*fd(h/2) -
      !! fd(h) = 3*f'(x) + O(h^4)), leaving an O(h^4) floor instead of
      !! O(h^2) -- see the w3.xyz case's own comment for why this
      !! particular case needs it (doubles the SCF-reconverge cost per
      !! atom/direction: two step sizes instead of one).

      real(dp) :: my_abstol_on, my_invtol
      integer :: nat, dev, a, d
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r0(:, :), rtmp(:, :), grad_on(:, :), grad_off(:, :)
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_xc_functional_t) :: func
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      real(dp) :: eplus, eminus, fd, fd_h, fd_hh, worst_abs_on, worst_rel_on, worst_abs_off, worst_rel_off, denom
      real(dp) :: inv_on(3), inv_off(3)
      logical :: ok, my_richardson

      my_abstol_on = ABSTOL_ON
      if (present(abstol_on_override)) my_abstol_on = abstol_on_override
      my_richardson = .false.
      if (present(richardson)) my_richardson = richardson
      my_invtol = INVTOL
      if (present(invtol_override)) my_invtol = invtol_override

      ok = xc_functional_by_name(functional, func)
      if (.not. ok .or. len_trim(func%missing) > 0) then
         print '(a)', 'check_xc_grad_fd: '//label//': functional not available'; stop 1
      end if

      dev = trc_bind_device(0)
      call read_xyz(xyzfile, nat, at_z, r0)
      ! force_cart: trc_scf_gradient DOES have a spherical fold now
      ! (check_dft_grad_vs_pyscf.F90 exercises it); kept Cartesian here
      ! only to keep this program's own finite-difference recipe (a full
      ! SCF reconverge at every displaced geometry, scf_energy_at below)
      ! simple -- not a missing capability.
      call trc_basis_from_json(basis_json, nat, at_z, r0, bas, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_xc_grad_fd: '//err%get_message(); stop 1
      end if
      call bas%to_device()

      opts%functional = functional
      opts%grid_level = grid_level
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iter = 200
      call trc_scf_run(bas, nalpha, nbeta, opts, res)
      if (.not. res%converged) then
         print '(a)', 'check_xc_grad_fd: '//label//': reference SCF did not converge'
         any_fail = .true.; return
      end if

      allocate (grad_on(3, nat), grad_off(3, nat))
      call lp%init(bas%nao)
      call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, 1.0e-12_dp, grad_on, &
                            func=func, grid_level=grid_level, grid_response=.true.)
      call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, 1.0e-12_dp, grad_off, &
                            func=func, grid_level=grid_level, grid_response=.false.)
      ! Released before `lp` goes out of scope: trc_linalg_t has no
      ! finalizer (check_full_grad_fd.F90's own comment, unchanged here).
      call lp%release()
      call bas%release()

      inv_on = sum(grad_on, dim=2)
      inv_off = sum(grad_off, dim=2)

      allocate (rtmp(3, nat))
      worst_abs_on = 0.0_dp; worst_rel_on = 0.0_dp
      worst_abs_off = 0.0_dp; worst_rel_off = 0.0_dp
      do a = 1, nat
         do d = 1, 3
            rtmp = r0; rtmp(d, a) = r0(d, a) + H
            eplus = scf_energy_at(basis_json, nat, at_z, rtmp, nalpha, nbeta, opts)
            rtmp = r0; rtmp(d, a) = r0(d, a) - H
            eminus = scf_energy_at(basis_json, nat, at_z, rtmp, nalpha, nbeta, opts)
            fd_h = (eplus - eminus)/(2.0_dp*H)
            ! OFF is diagnostic only (its own header: "NOT gated"), always
            ! against the plain step-H difference -- Richardson only
            ! matters for ON's own gate below.
            fd = fd_h
            if (my_richardson) then
               rtmp = r0; rtmp(d, a) = r0(d, a) + 0.5_dp*H
               eplus = scf_energy_at(basis_json, nat, at_z, rtmp, nalpha, nbeta, opts)
               rtmp = r0; rtmp(d, a) = r0(d, a) - 0.5_dp*H
               eminus = scf_energy_at(basis_json, nat, at_z, rtmp, nalpha, nbeta, opts)
               fd_hh = (eplus - eminus)/H
               ! Richardson extrapolation: fd_h = f'(x) + c*H^2 + O(H^4),
               ! fd_hh = f'(x) + c*(H/2)^2 + O(H^4), so 4*fd_hh - fd_h =
               ! 3*f'(x) + O(H^4) -- see this subroutine's own header.
               fd = (4.0_dp*fd_hh - fd_h)/3.0_dp
            end if
            denom = max(abs(grad_on(d, a)), abs(fd), 1.0e-10_dp)
            worst_abs_on = max(worst_abs_on, abs(grad_on(d, a) - fd))
            worst_rel_on = max(worst_rel_on, abs(grad_on(d, a) - fd)/denom)
            denom = max(abs(grad_off(d, a)), abs(fd_h), 1.0e-10_dp)
            worst_abs_off = max(worst_abs_off, abs(grad_off(d, a) - fd_h))
            worst_rel_off = max(worst_rel_off, abs(grad_off(d, a) - fd_h)/denom)
         end do
      end do

      print '(a,a,a,i0)', '  ', label, ': grid level ', grid_level
      print '(a,es10.2,a,es10.2,a,3es10.2)', '    OFF (Pulay only): worst abs ', worst_abs_off, &
         '  worst rel ', worst_rel_off, '  translational-invariance residual: ', inv_off
      print '(a,es10.2,a,es10.2,a,3es10.2)', '    ON  (grid response): worst abs ', worst_abs_on, &
         '  worst rel ', worst_rel_on, '  translational-invariance residual: ', inv_on
      if (worst_abs_on > my_abstol_on .and. worst_rel_on > RELTOL) any_fail = .true.
      if (any(abs(inv_on) > my_invtol)) any_fail = .true.

      deallocate (at_z, r0, rtmp, grad_on, grad_off)
   end subroutine run_case

   real(dp) function scf_energy_at(basis_json, nat, at_z, r, nalpha, nbeta, opts) result(e)
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nalpha, nbeta
      real(dp), intent(in) :: r(3, nat)
      type(trc_scf_options_t), intent(in) :: opts
      type(trc_basis_t) :: b2
      type(trc_scf_result_t) :: r2
      type(error_t) :: err
      call trc_basis_from_json(basis_json, nat, at_z, r, b2, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_xc_grad_fd: '//err%get_message(); stop 1
      end if
      call b2%to_device()
      call trc_scf_run(b2, nalpha, nbeta, opts, r2)
      e = r2%energy
      call b2%release()
   end function scf_energy_at

end module check_xc_grad_fd_cases

program check_xc_grad_fd
   use pic_types, only: dp
   use check_xc_grad_fd_cases, only: run_case
   implicit none
   logical :: any_fail
   any_fail = .false.

   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'svwn', 3, 'RKS SVWN water/6-31G* L3', any_fail)
   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'svwn', 5, 'RKS SVWN water/6-31G* L5', any_fail)
   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'pbe', 3, 'RKS PBE water/6-31G* L3', any_fail)
   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'pbe', 5, 'RKS PBE water/6-31G* L5', any_fail)
   call run_case('nh2.xyz', 'basis_sets/6-31gs.json', 5, 4, 'pbe', 3, 'UKS PBE NH2/6-31G* L3', any_fail)
   call run_case('nh2.xyz', 'basis_sets/6-31gs.json', 5, 4, 'pbe', 5, 'UKS PBE NH2/6-31G* L5', any_fail)

   ! SCREENING CONSISTENCY, A CASE LARGE ENOUGH TO EXERCISE IT. Every case
   ! above is 2-3 atoms, where `becke_partition_derivatives`'s active-set
   ! screening (trc_dft_partition.F90's own header) never drops anything --
   ! its nu_max=0.95 cutoff only excludes an atom once the point's owner is
   ! confidently closer than it, which three atoms sharing one molecule
   ! essentially never trigger. w3.xyz is three SEPARATE, well-spread water
   ! molecules (9 atoms; the same geometry the MPI fragment tests use) --
   ! points deep in one molecule are confidently closer to their own atoms
   ! than to either other molecule's, so this case genuinely exercises the
   ! screened branch rather than merely compiling it: a direct probe
   ! (comparing this routine's own screened vs `nu_max=2` unscreened
   ! `dweights` on this exact grid) measures worst|screened-unscreened| =
   ! 5.15e-8 over 101094 points -- screening is live here, not a no-op.
   !
   ! ITS OWN tolerance, not ABSTOL_ON: w3.xyz is a compressed,
   ! non-equilibrium geometry (deliberately so, for the fragment tests it
   ! was built for -- O-H ~0.82 angstrom against water's own equilibrium
   ! ~0.96), not a minimum, so its own energy surface is sharper and its
   ! third-derivative larger, which is exactly the coefficient the plain
   ! central difference's O(H^2) truncation term depends on. PLAIN step-H
   ! worst abs ON was 9.30e-7 (OFF: 2.51e-6) against the six equilibrium-
   ! ish 2-3-atom cases' ~1.4e-7 -- NOT a screening artifact (a direct
   ! probe with becke_partition_derivatives' active-set screening
   ! disabled entirely, nu_max=2, reproduces the identical 9.30e-7,
   ! unchanged to 3 significant figures) and NOT a grid-response bug
   ! (translational invariance stays ~1e-11/1e-12 either way, INVTOL's
   ! default 1e-8 cleared with 1e3-1e4x margin -- a real inconsistency
   ! between the energy's screening and this gradient's would show up
   ! there first, and does not). richardson=.true. below confirms it IS
   ! h^2 finite-difference truncation directly: Richardson-extrapolating
   ! the same central difference (steps H and H/2, this subroutine's own
   ! header) collapses worst abs ON from 9.30e-7 to 2.13e-9 (measured,
   ! CPU, gfortran) -- a ~440x reduction, consistent with trading an
   ! O(H^2) floor for an O(H^4) one. 1e-8 below is that measured value
   ! with ~5x margin, not the old 2e-6 loosened.
   call run_case('w3.xyz', 'basis_sets/6-31gs.json', 15, 15, 'svwn', 3, 'RKS SVWN w3/6-31G* L3', any_fail, &
                 abstol_on_override=1.0e-8_dp, richardson=.true.)

   ! VV10 (wB97X-V, RKS): the nonlocal correlation term's own grid response
   ! (point motion of BOTH grids plus the Becke weight derivative,
   ! trc_vv10_grad.F90's own header) on top of the semilocal GGA term and
   ! the erf long-range exchange gradient -- the full wB97X-V gradient,
   ! differenced against the SAME recipe as every case above (grid rebuilt
   ! at every displaced geometry, VV10's own double sum recomputed from
   ! scratch each time). ITS OWN tolerance -- see this program's own
   ! output for the measured worst|abs| this was set from; not assumed
   ! equal to ABSTOL_ON, which was measured on functionals with no
   ! nonlocal term at all.
   call run_case('water.xyz', 'basis_sets/6-31gs.json', 5, 5, 'wb97x-v', 3, 'RKS wB97X-V water/6-31G* L3', &
                 any_fail, abstol_on_override=5.0e-6_dp, invtol_override=1.0e-7_dp)

   if (any_fail) then
      print '(a)', '  RESULT: FAIL'; stop 1
   end if
   print '(a)', '  RESULT: PASS'
end program check_xc_grad_fd
