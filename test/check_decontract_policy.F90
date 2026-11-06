!
! The split policy of trc_decontract: which groups it splits, and that the
! choice never shows in the answer.
!
! `trc_decontract_options_t` decides how the four-centre engine's working
! basis is cut up -- whether a general contraction is split at all (the shell
! price), and whether diffuse primitives leave a segmented shell (the
! diffuse split). The integrals over any of those working bases are exact,
! so the claims are:
!
! - each policy splits what it says: the shell counts below are worked out
!   by hand from the basis files;
! - the split basis still spans the caller's functions exactly: C^T S_p C = S
!   from trc_1e, which shares nothing with the four-centre path (the same
!   external answer check_decontract uses);
! - an RHF energy is the same under every policy to round-off, on a small
!   diffuse case where the diffuse split and the general-contraction split
!   both act.
!
! A test-drive suite.
!
module check_decontract_policy_suite
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_basis_json, only: trc_basis_from_json
   use trc_decontract, only: decon_t, decontract_basis, decontract_fold_host, decontract_release, &
                             trc_decontract_options_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_decontract_policy

   real(dp), parameter :: SPAN_TOL = 1.0e-12_dp
   real(dp), parameter :: ENERGY_TOL = 1.0e-10_dp
   real(dp), parameter :: WHOLE_TOL = 3.0e-8_dp
      !! w5/cc-pVDZ at an ERI threshold of 1e-8, gfortran 15: O's s kept
      !! whole lands 9.9e-09 Eh from it split, 8.8e-08 without the
      !! amplitude `screen_whole` gives it. Midway, a factor of 3 each side.
   real(dp), parameter :: NEVER = 1.0e6_dp
      !! A shell price no split can beat: the caller's basis, unsplit.

contains

   subroutine collect_decontract_policy(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("the shell price splits Si s and keeps C/N/O s whole", test_shell_price), &
                  new_unittest("the diffuse split moves diffuse primitives out", test_diffuse_split), &
                  new_unittest("the diffuse split spans the same functions", test_diffuse_span), &
                  new_unittest("RHF energy is the same under every policy", test_energy_invariant), &
                  new_unittest("kept whole screens as tightly as split", test_screen_invariant) &
                  ]
   end subroutine collect_decontract_policy

   subroutine load(xyz, json, b)
      character(*), intent(in) :: xyz, json
      type(trc_basis_t), intent(out) :: b
      type(error_t) :: err
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      call read_xyz(xyz, natm, zint, at_r)
      call trc_basis_from_json(json, natm, zint, at_r, b, err)
      if (err%has_error()) error stop 'check_decontract_policy: basis build failed'
   end subroutine load

   integer function nshell_split(b, opts, active) result(n)
      !! Shells in the working basis `opts` gives, and whether it split.
      type(trc_basis_t), intent(in) :: b
      type(trc_decontract_options_t), intent(in) :: opts
      logical, intent(out) :: active
      type(trc_basis_t) :: pb
      type(decon_t) :: m
      call decontract_basis(b, pb, m, opts)
      active = m%active
      n = b%nshell
      if (active) n = pb%nshell
      call decontract_release(m)
   end function nshell_split

   !
   ! After trc_basis_json's pruning, the only general contractions in
   ! cc-pVDZ are O's s (8 primitives under 2 columns, all compact) and Si's
   ! s (11 under 3) and p (7 under 2), both of which share a diffuse
   ! primitive between their columns. At the default price of 2 primitives
   ! a shell O's s is kept (16 + 2*2 = 20 against 8*3 = 24). Si's s and p
   ! split whatever the price, because their columns share a primitive
   ! below the 0.3 cutoff (3 -> 11 and 2 -> 7 shells); only a cutoff below
   ! it hands them back to the price. A price of zero splits everything, as
   ! the module did before the price existed.
   !
   subroutine test_shell_price(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(trc_decontract_options_t) :: opts
      logical :: active
      integer :: n

      call load('water.xyz', 'basis_sets/cc-pvdz.json', b)
      n = nshell_split(b, opts, active)
      call check(error, .not. active, "water/cc-pVDZ: the default price keeps O's s whole")
      if (allocated(error)) return
      opts%shell_cost = 0.0_dp
      n = nshell_split(b, opts, active)
      call check(error, active, "water/cc-pVDZ: price 0 splits O's s")
      if (allocated(error)) return
      call check(error, n, b%nshell - 2 + 8, "water/cc-pVDZ: price 0 gives O's s 8 shells")
      if (allocated(error)) return
      ! An MD-built Coulomb matrix splits whatever the price.
      opts = trc_decontract_options_t()
      opts%md_coulomb = .true.
      n = nshell_split(b, opts, active)
      call check(error, n, b%nshell - 2 + 8, "water/cc-pVDZ: md_coulomb splits O's s at the default price")
      if (allocated(error)) return

      call load('sih4.xyz', 'basis_sets/cc-pvdz.json', b)
      opts = trc_decontract_options_t()
      n = nshell_split(b, opts, active)
      call check(error, n, b%nshell - 3 + 11 - 2 + 7, "SiH4/cc-pVDZ: Si s and p, sharing a diffuse primitive, split")
      if (allocated(error)) return
      opts%shell_cost = NEVER
      n = nshell_split(b, opts, active)
      call check(error, n, b%nshell - 3 + 11 - 2 + 7, "SiH4/cc-pVDZ: a shared diffuse primitive splits at any price")
      if (allocated(error)) return
      ! Below every shared exponent (Si s 0.2517, p 0.2964), only the price
      ! decides: at 0.2 the s (33 + 3*2 = 39 against 33) splits, the p
      ! (14 + 2*2 = 18 against 21) stays whole.
      opts = trc_decontract_options_t()
      opts%cutoff = 0.2_dp
      n = nshell_split(b, opts, active)
      call check(error, n, b%nshell - 3 + 11, "SiH4/cc-pVDZ: at cutoff 0.2 Si s splits on price, p is kept")
      if (allocated(error)) return
      opts%shell_cost = NEVER
      n = nshell_split(b, opts, active)
      call check(error, .not. active, "SiH4/cc-pVDZ: at cutoff 0.2 a price nothing beats splits nothing")
   end subroutine test_shell_price

   !
   ! 6-31G* is segmented. Si's 3-primitive s and p shells (1.727, 0.5729,
   ! 0.2222) each hold one primitive below the 0.3 cutoff: the diffuse split
   ! makes each a 2-primitive segment plus a shell of one, two more shells.
   ! H's 3-primitive s goes down to 0.6401 and is left alone at 0.3. At a
   ! cutoff of 0.7 Si's s and p keep one compact primitive and lose two
   ! (four more shells), and every H s splits too (one more each).
   !
   subroutine test_diffuse_split(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(trc_decontract_options_t) :: opts
      logical :: active
      integer :: n

      call load('sih4.xyz', 'basis_sets/6-31gs.json', b)
      n = nshell_split(b, opts, active)
      call check(error, .not. active, "SiH4/6-31G*: no split with the diffuse split off")
      if (allocated(error)) return
      opts%diffuse = .true.
      n = nshell_split(b, opts, active)
      call check(error, active, "SiH4/6-31G*: the diffuse split acts")
      if (allocated(error)) return
      call check(error, n, b%nshell + 2, "SiH4/6-31G*: Si s and p each gain a diffuse shell")
      if (allocated(error)) return
      opts%cutoff = 0.7_dp
      n = nshell_split(b, opts, active)
      call check(error, n, b%nshell + 4 + 4, "SiH4/6-31G*: at 0.7 Si sp and every H s split")
   end subroutine test_diffuse_split

   subroutine test_diffuse_span(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b, pb
      type(trc_pairlist_t) :: p, pp
      type(decon_t) :: m
      type(trc_decontract_options_t) :: opts
      real(dp), allocatable :: s(:, :), t(:, :), v(:, :), sp(:, :), tp(:, :), vp(:, :), f(:, :)
      real(dp) :: worst

      call load('sih4.xyz', 'basis_sets/6-31gs.json', b)
      opts%diffuse = .true.
      opts%cutoff = 0.7_dp
      call decontract_basis(b, pb, m, opts)
      call check(error, m%active, "SiH4/6-31G*: the diffuse split acts")
      if (allocated(error)) return

      call b%to_device(); call pb%to_device()
      call p%build(b, 1.0e-14_dp); call pp%build(pb, 1.0e-14_dp)
      allocate (s(b%nao, b%nao), t(b%nao, b%nao), v(b%nao, b%nao), f(b%nao, b%nao))
      allocate (sp(pb%nao, pb%nao), tp(pb%nao, pb%nao), vp(pb%nao, pb%nao))
      call trc_1e(b, p, s, t, v)
      call trc_1e(pb, pp, sp, tp, vp)
      call decontract_fold_host(m, sp, f)
      worst = maxval(abs(f - s))
      call check(error, worst, 0.0_dp, thr=SPAN_TOL, message="worst |C^T S_p C - S| exceeds tolerance")
      if (allocated(error)) return
      call decontract_fold_host(m, vp, f)
      worst = maxval(abs(f - v))
      call check(error, worst, 0.0_dp, thr=SPAN_TOL*maxval(abs(v)), &
                 message="worst |C^T V_p C - V| exceeds tolerance")
      call decontract_release(m)
      call p%release(); call pp%release()
      call b%release(); call pb%release()
   end subroutine test_diffuse_span

   subroutine run_rhf(xyz, json, nocc, opts, energy, eri_thresh)
      character(*), intent(in) :: xyz, json
      integer, intent(in) :: nocc
      type(trc_decontract_options_t), intent(in) :: opts
      real(dp), intent(out) :: energy
      real(dp), intent(in), optional :: eri_thresh
      type(trc_basis_t) :: b
      type(trc_scf_options_t) :: scf
      type(trc_scf_result_t) :: res

      call load(xyz, json, b)
      call b%to_device()
      scf%conv_energy = 1.0e-12_dp
      scf%conv_density = 1.0e-10_dp
      scf%conv_diis = 1.0e-10_dp
      scf%max_iteration = 200
      scf%decontract = opts
      if (present(eri_thresh)) then
         ! A loose screen has a noise floor of its own: converge to it, not
         ! past it.
         scf%eri_thresh = eri_thresh
         scf%conv_diis = 1.0e-5_dp
         scf%conv_energy = 1.0e-9_dp
      end if
      call trc_scf_run(b, nocc, nocc, scf, res)
      call b%release()
      if (.not. res%converged) error stop 'check_decontract_policy: RHF did not converge'
      energy = res%energy
   end subroutine run_rhf

   !
   ! Water in d-aug-cc-pVDZ: two diffuse shells per l on top of cc-pVDZ.
   ! At a cutoff of 0.5 the diffuse split takes H's 0.4446 out of its s
   ! segment (every other diffuse primitive is already a shell of one), and
   ! a price of zero splits O's general s. The unsplit basis (a price
   ! nothing beats) is the reference.
   !
   subroutine test_energy_invariant(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_decontract_options_t) :: opts
      real(dp) :: e_ref, e

      opts%shell_cost = NEVER
      call run_rhf('water.xyz', 'basis_sets/d-aug-cc-pvdz.json', 5, opts, e_ref)

      opts = trc_decontract_options_t()
      opts%shell_cost = 0.0_dp
      call run_rhf('water.xyz', 'basis_sets/d-aug-cc-pvdz.json', 5, opts, e)
      call check(error, e, e_ref, thr=ENERGY_TOL, message="water/d-aug-cc-pVDZ: O s split moves E_RHF")
      if (allocated(error)) return

      opts = trc_decontract_options_t()
      opts%diffuse = .true.
      opts%cutoff = 0.5_dp
      call run_rhf('water.xyz', 'basis_sets/d-aug-cc-pvdz.json', 5, opts, e)
      call check(error, e, e_ref, thr=ENERGY_TOL, message="water/d-aug-cc-pVDZ: diffuse split moves E_RHF")
      if (allocated(error)) return

      opts%shell_cost = 0.0_dp
      call run_rhf('water.xyz', 'basis_sets/d-aug-cc-pvdz.json', 5, opts, e)
      call check(error, e, e_ref, thr=ENERGY_TOL, message="water/d-aug-cc-pVDZ: both splits move E_RHF")
      if (allocated(error)) return

      ! Si: the general s and p and the 6-31G* diffuse sp, all at once.
      opts = trc_decontract_options_t()
      opts%shell_cost = NEVER
      call run_rhf('sih4.xyz', 'basis_sets/6-31gs.json', 9, opts, e_ref)
      opts%shell_cost = 0.0_dp
      opts%diffuse = .true.
      call run_rhf('sih4.xyz', 'basis_sets/6-31gs.json', 9, opts, e)
      call check(error, e, e_ref, thr=ENERGY_TOL, message="SiH4/6-31G*: diffuse split moves E_RHF")
   end subroutine test_energy_invariant

   !
   ! The point of `screen_whole`: at a loose screening threshold the policy
   ! must not move the energy. A general contraction kept whole screens on
   ! its contracted bound, which drops far more than the split's
   ! amplitude-scaled primitive bounds; carrying the amplitude the split
   ! would have keeps the two together. Five waters in cc-pVDZ at 1e-8, so
   ! that screening, not round-off, sets the difference: the energies with
   ! O's s split and kept whole must agree far better than kept whole
   ! without the amplitude does.
   !
   subroutine test_screen_invariant(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_decontract_options_t) :: opts
      real(dp), parameter :: LOOSE = 1.0e-8_dp
      real(dp) :: e_split, e_whole, e_bare

      opts%shell_cost = 0.0_dp
      call run_rhf('w5.xyz', 'basis_sets/cc-pvdz.json', 25, opts, e_split, LOOSE)
      opts = trc_decontract_options_t()
      call run_rhf('w5.xyz', 'basis_sets/cc-pvdz.json', 25, opts, e_whole, LOOSE)
      opts%screen_whole = .false.
      call run_rhf('w5.xyz', 'basis_sets/cc-pvdz.json', 25, opts, e_bare, LOOSE)
      print '(a,es10.3,a,es10.3)', '  w5/cc-pVDZ at 1e-8: |E_whole - E_split| ', abs(e_whole - e_split), &
         '   without the amplitude ', abs(e_bare - e_split)
      call check(error, abs(e_whole - e_split) < WHOLE_TOL, &
                 "w5/cc-pVDZ: kept whole moves E_RHF at a loose threshold")
      if (allocated(error)) return
      ! The comparison above can fail: without the amplitude it does.
      call check(error, abs(e_bare - e_split) > WHOLE_TOL, &
                 "w5/cc-pVDZ: the loose threshold no longer exposes the bare whole-shell screen")
   end subroutine test_screen_invariant

end module check_decontract_policy_suite

program check_decontract_policy
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_bind_device
   use check_decontract_policy_suite, only: collect_decontract_policy
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, dev

   dev = trc_bind_device(0)
   stat = 0
   testsuites = [new_testsuite("check_decontract_policy", collect_decontract_policy)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_decontract_policy
