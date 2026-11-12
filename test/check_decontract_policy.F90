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
                             trc_decontract_options_t, decontract_whole_screen
   use trc_hgp, only: hgp_build_pairs, trc_pair_prune_t
   use trc_eri_data, only: trc_shell_pairs_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   implicit none
   private
   public :: collect_decontract_policy

   real(dp), parameter :: SPAN_TOL = 1.0e-12_dp
   real(dp), parameter :: ENERGY_TOL = 1.0e-10_dp
   real(dp), parameter :: WHOLE_SLACK = 1.1_dp, WHOLE_FLOOR = 1.0e-10_dp
      !! Kept whole may be at most 10% less accurate than split, plus the
      !! convergence floor of these runs (conv_energy 1e-11).
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
                  new_unittest("primitive pairs come in the kernels' |c|/p order", test_pair_order), &
                  new_unittest("a shell kept whole carries the split's primitive screen", test_whole_screen), &
                  new_unittest("kept whole gives the split's energy at a loose threshold", test_screen_invariant) &
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
   ! cc-pVDZ are O's s (8 primitives under 2 columns) and Si's s (11 under
   ! 3) and p (7 under 2). At the default price of 2 primitives a shell O's
   ! s is kept (16 + 2*2 = 20 against 8*3 = 24) and so is Si's p (14 + 4 =
   ! 18 against 21). Si's s splits whatever the price: its 1s column gives
   ! the outer primitives coefficients a thousand times smaller than the
   ! 3s column does (group_unevenness 998, above DECON_WHOLE_UNEVEN_MAX).
   ! A price of zero splits everything.
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
      call check(error, n, b%nshell - 3 + 11, "SiH4/cc-pVDZ: Si s split, p kept at the default price")
      if (allocated(error)) return
      opts%shell_cost = NEVER
      n = nshell_split(b, opts, active)
      call check(error, n, b%nshell - 3 + 11, "SiH4/cc-pVDZ: Si s, too uneven to keep, splits at any price")
      if (allocated(error)) return
      opts%shell_cost = 0.0_dp
      n = nshell_split(b, opts, active)
      call check(error, n, b%nshell - 3 + 11 - 2 + 7, "SiH4/cc-pVDZ: price 0 also splits Si's p")
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
         scf%conv_diis = 1.0e-6_dp
         scf%conv_energy = 1.0e-11_dp
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
   ! THE PRIMITIVE PAIRS OF A SHELL KEPT WHOLE ARE IN KERNEL ORDER.
   !
   ! The kernels walk a shell pair's primitive pairs and EXIT at the first
   ! one whose |c|/p bound fails, so hgp_build_pairs must hand them over in
   ! non-increasing |c|/p. It sorted by |c|^2/p, which no split basis could
   ! show (one primitive pair per shell pair) and a general contraction kept
   ! whole did: the silica slice was 9.2e-08 Eh off at 1e-11. Water in
   ! cc-pVDZ, O's s kept whole, every shell pair checked.
   !
   subroutine test_pair_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(trc_shell_pairs_t) :: pairs
      type(trc_pair_prune_t) :: prune
      real(dp), allocatable :: one(:)
      integer :: key, k, nbad

      call load('water.xyz', 'basis_sets/cc-pvdz.json', b)
      allocate (one(b%nshell))
      one = 1.0_dp
      ! Component by component: nvfortran 26.5 has an internal compiler
      ! error (search_arr) on this structure constructor with a basis-array argument.
      prune%camp = b%sh_c
      prune%cut = 1.0e-20_dp
      call hgp_build_pairs(b%sh_np, b%sh_e, b%sh_c, b%sh_r, one, prune, pairs)
      nbad = 0
      do key = 1, b%nshell*b%nshell
         do k = pairs%off(key) + 1, pairs%off(key) + pairs%n(key) - 1
            if (abs(pairs%c(k + 1))/pairs%p(k + 1) > abs(pairs%c(k))/pairs%p(k)*(1.0_dp + 1.0e-14_dp)) nbad = nbad + 1
         end do
      end do
      call check(error, nbad, 0, "hgp_build_pairs: primitive pairs out of |c|/p order, the kernels' exit "// &
                 "would drop live quartets")
   end subroutine test_pair_order

   !
   ! THE PRIMITIVE SCREEN OF A SHELL KEPT WHOLE IS THE SPLIT'S.
   !
   ! decontract_whole_screen gives every primitive of a general contraction
   ! kept whole the coefficient the split would (the largest any column
   ! gives it) as its prune amplitude, and reports how much smaller a
   ! column's own coefficient can be (`uneven`), which tightens the kernels'
   ! cut by uneven**4. Checked on the shapes the default policy keeps whole:
   ! O's cc-pVDZ s (water) and Si's p (SiH4), and on a segmented basis, which
   ! must come back untouched.
   !
   subroutine test_whole_screen(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b, pb
      type(decon_t) :: m
      type(trc_decontract_options_t) :: opts
      real(dp), allocatable :: camp(:, :)
      real(dp) :: uneven
      integer :: i, j, k, kk

      call load('water.xyz', 'basis_sets/cc-pvdz.json', b)
      allocate (camp(b%maxnp, b%nshell))
      call decontract_whole_screen(b, camp, uneven)
      call check(error, uneven > 65.0_dp .and. uneven < 75.0_dp, &
                 "water/cc-pVDZ: O's s unevenness is not the 71 its columns give: "//to_string(uneven))
      if (allocated(error)) return
      ! Every primitive of O's two s columns carries the larger of the two
      ! coefficients the columns give it, and nothing is ever below its own.
      do i = 1, b%nshell
         do j = 1, b%nshell
            if (i == j .or. b%sh_l(i) /= 0 .or. b%sh_l(j) /= 0 .or. b%sh_np(i) < 2 .or. b%sh_np(j) < 2) cycle
            if (b%sh_at(i) /= b%sh_at(j)) cycle
            do k = 1, b%sh_np(i)
               do kk = 1, b%sh_np(j)
                  if (b%sh_e(k, i) /= b%sh_e(kk, j)) cycle
                  call check(error, camp(k, i), max(abs(b%sh_c(k, i)), abs(b%sh_c(kk, j))), &
                             message="water/cc-pVDZ: a kept-whole primitive does not carry the split's coefficient")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
      deallocate (camp)

      ! SiH4: the default policy splits Si's s and keeps its p whole; the
      ! working basis's own groups are then exactly that p.
      call load('sih4.xyz', 'basis_sets/cc-pvdz.json', b)
      call decontract_basis(b, pb, m, opts)
      allocate (camp(pb%maxnp, pb%nshell))
      call decontract_whole_screen(pb, camp, uneven)
      call check(error, uneven > 48.0_dp .and. uneven < 58.0_dp, &
                 "SiH4/cc-pVDZ: the kept-whole Si p unevenness is not the 53 its columns give: "//to_string(uneven))
      if (allocated(error)) return
      call decontract_release(m)
      deallocate (camp)

      call load('water.xyz', 'basis_sets/6-31g.json', b)
      allocate (camp(b%maxnp, b%nshell))
      call decontract_whole_screen(b, camp, uneven)
      call check(error, uneven, 1.0_dp, message="water/6-31G: a segmented basis must not tighten the cut")
      if (allocated(error)) return
      call check(error, all(camp == b%sh_c), "water/6-31G: a segmented basis must keep its own amplitudes")
   end subroutine test_whole_screen

   !
   ! THE ENERGY OF A SHELL KEPT WHOLE. At a loose threshold, where screening
   ! and not round-off sets the error, keeping a general contraction whole
   ! must be no less accurate than splitting it, judged against the same
   ! system at 1e-14: O's s (w5), S's p ((H2S)4) and Si's p (disiloxane),
   ! all cc-pVDZ, at 1e-8. Measured on w5 (V100): split 6.9e-08 Eh off,
   ! kept whole 6.1e-08, and kept whole without `screen_whole` 1.4e-07 --
   ! which the last check asserts, so that this test can fail.
   !
   subroutine test_screen_invariant(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: bare, split
      call screen_case(error, 'w5.xyz', 25, 'w5/cc-pVDZ', bare, split)
      if (allocated(error)) return
      call check(error, bare > split*WHOLE_SLACK + WHOLE_FLOOR, "w5/cc-pVDZ: without screen_whole the "// &
                 "kept-whole s is no longer less accurate, so this test cannot fail: "//to_string(bare)// &
                 " against "//to_string(split))
      if (allocated(error)) return
      call screen_case(error, 'h2s4.xyz', 36, '(H2S)4/cc-pVDZ', bare, split)
      if (allocated(error)) return
      call screen_case(error, 'disiloxane.xyz', 21, 'disiloxane/cc-pVDZ', bare, split)
   end subroutine test_screen_invariant

   subroutine screen_case(error, xyz, nocc, label, bare, split)
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: xyz, label
      integer, intent(in) :: nocc
      real(dp), intent(out) :: bare, split
      !! |E - E(1e-14)| kept whole without `screen_whole`, and split.
      type(trc_decontract_options_t) :: opts
      real(dp), parameter :: LOOSE = 1.0e-8_dp, TIGHT = 1.0e-14_dp
      real(dp) :: e_ref, e_split, e_whole, e_bare, whole

      opts%shell_cost = 0.0_dp
      call run_rhf(xyz, 'basis_sets/cc-pvdz.json', nocc, opts, e_ref, TIGHT)
      call run_rhf(xyz, 'basis_sets/cc-pvdz.json', nocc, opts, e_split, LOOSE)
      opts = trc_decontract_options_t()
      call run_rhf(xyz, 'basis_sets/cc-pvdz.json', nocc, opts, e_whole, LOOSE)
      opts%screen_whole = .false.
      call run_rhf(xyz, 'basis_sets/cc-pvdz.json', nocc, opts, e_bare, LOOSE)
      split = abs(e_split - e_ref)
      whole = abs(e_whole - e_ref)
      bare = abs(e_bare - e_ref)
      call global_logger%info("  "//label//" at 1e-8, |E - E(1e-14)|: split "//to_string(split)// &
                              ", kept whole "//to_string(whole)//", without screen_whole "//to_string(bare))
      call check(error, whole <= split*WHOLE_SLACK + WHOLE_FLOOR, label//": kept whole is less accurate "// &
                 "than split at a loose threshold: "//to_string(whole)//" against "//to_string(split))
   end subroutine screen_case

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
