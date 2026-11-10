!
! The unique-pair build (docs/SYMMETRY_UNIQUE_PAIRS.md): with a point group
! `eri%build` forms the primitive pairs and the Schwarz bounds of the orbit
! representatives only and rebuilds the images from them (the product centre
! and the shell centres sign-flipped by the operation's mask). Checked here,
! per molecule, in cc-pVDZ (Cartesian caller space):
!
!   - the path ran: `eri%upair`, the number of shell pairs whose primitive
!     pairs were formed is below nshell**2, `pairs` holds fewer primitive pairs
!     than the pairs of the bins use, by at least the ratios measured when this
!     was written (measured 1.21, 2.21, 2.31 and 2.62 for water, B2H6, benzene and
!     naphthalene; the thresholds sit a little below) -- a silent fall back to the merged build
!     would give one;
!   - the tables are the full build's: the bin-interleaved primitive-pair
!     tables rebuilt from `pairs` of every shell pair (hgp_build_pairs without
!     the wanted mask) equal the unique build's, array by array, to 1e-13
!     (they are bit-identical when this was written), and the primitive counts
!     exactly;
!   - the paths that read `pairs` itself still work: the non-symmetric Fock
!     build on the unique engine (it completes the pairs first, and
!     `eri%upair` is false after) equals the one on an engine built without
!     symmetry, to the 1e-8 the merged check allows for the orbit-maximum Q.
!
! The G-equals-full-G oracle for the skeleton build on this path is
! check_sym_fock (its default engine is the unique-pair one); it also runs
! under TRC_UNIQUE_PAIRS=0 as check_sym_fock_nouniq.
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_unique_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_bins, only: bins_build_pair_quartets
   use trc_eri, only: trc_eri_t
   use trc_eri_data, only: trc_shell_pairs_t, trc_pair_quartets_t
   use trc_error, only: error_t
   use trc_hgp, only: hgp_build_pairs, trc_pair_prune_t
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_sym_unique

   real(dp), parameter :: THRESH = 1.0e-10_dp
   real(dp), parameter :: TOL_TABLE = 1.0e-13_dp
   real(dp), parameter :: TOL_REF = 1.0e-8_dp

contains

   subroutine collect_sym_unique(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water C2v: unique pairs, tables, non-symmetric build", test_water), &
                  new_unittest("B2H6 D2h (nonplanar)", test_b2h6), &
                  new_unittest("benzene D2h (via S12)", test_benzene), &
                  new_unittest("naphthalene D2h", test_naphthalene) &
                  ]
   end subroutine collect_sym_unique

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case("water", "C2v", 1.15_dp, error)
   end subroutine test_water

   subroutine test_b2h6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case("b2h6_d2h", "D2h", 2.1_dp, error)
   end subroutine test_b2h6

   subroutine test_benzene(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case("tops_benzene", "D2h", 2.2_dp, error)
   end subroutine test_benzene

   subroutine test_naphthalene(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case("naphthalene_d2h", "D2h", 2.5_dp, error)
   end subroutine test_naphthalene

   subroutine run_case(name, group, min_ratio, error)
      character(len=*), intent(in) :: name, group
      real(dp), intent(in) :: min_ratio
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_eri_t) :: eri, eri0
      type(error_t) :: err
      integer, allocatable :: elem(:)
      integer :: nat, dev, n
      real(dp), allocatable :: r(:, :), d(:, :), g1(:, :), g0(:, :)
      real(dp) :: ratio, e, gmax

      dev = trc_bind_device(0)
      call read_xyz(trim(name)//".xyz", nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, trim(group), "group of "//trim(name))
      if (allocated(error)) return
      call trc_basis_from_json("basis_sets/cc-pvdz.json", nat, elem, geom%r_std, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      bas%spherical = .false.
      call bas%to_device()

      call eri%build(bas, THRESH, sym=geom)
      call check(error, eri%upair, trim(name)//": the unique-pair build did not run (TRC_UNIQUE_PAIRS=0 in the environment?)")
      if (allocated(error)) return
      call check(error, eri%upair_nshell_pairs < eri%nbas*eri%nbas, trim(name)//": every shell pair got primitive pairs")
      if (allocated(error)) return
      ratio = real(eri%upair_npp_bins, dp)/real(eri%upair_npp_stored, dp)
      print "(a,a,a,i0,a,i0,a,i0,a,i0,a,f6.2)", "  ", trim(name), ": shell pairs formed ", eri%upair_nshell_pairs, " of ", &
         eri%nbas*eri%nbas, "; primitive pairs stored ", eri%upair_npp_stored, " for ", eri%upair_npp_bins, &
         " in the bins, ratio ", ratio
      call check(error, ratio >= min_ratio, trim(name)//": primitive pairs stored vs used, ratio "//to_string(ratio, "f6.2")// &
                 " below "//to_string(min_ratio, "f6.2"))
      if (allocated(error)) return

      call check_tables(eri, bas, trim(name), error)
      if (allocated(error)) return

      ! the non-symmetric build completes the pairs
      call eri0%build(bas, THRESH)
      n = bas%nao
      allocate (d(n, n), g1(n, n), g0(n, n))
      call density(n, d)
      call eri%fock_nosym(bas, d, g1)
      call eri0%fock_nosym(bas, d, g0)
      gmax = maxval(abs(g0))
      e = maxval(abs(g1 - g0))
      print "(a,a,a,es10.3,a,es10.3)", "  ", trim(name), " fock_nosym, unique vs plain: max|G| ", gmax, "  err ", e
      call check(error, e <= TOL_REF*max(1.0_dp, gmax), trim(name)//": fock_nosym on the unique engine differs by "// &
                 to_string(e, "es10.3"))
      if (allocated(error)) return
      call check(error, .not. eri%upair, trim(name)//": the pairs were not completed for the non-symmetric build")
      if (allocated(error)) return
      call check(error, eri%pairs%npp >= eri%upair_npp_bins, &
         trim(name)//": completed pairs hold fewer primitive pairs than the bins use")

      ! the gradient reads the pairs through the bins too (completed by now,
      ! so rebuild the unique engine to see the gradient complete them itself)
      call eri%build(bas, THRESH, sym=geom)
      call check(error, eri%upair, trim(name)//": the second unique-pair build did not run")
      if (allocated(error)) return
      block
         real(dp), allocatable :: dd(:, :, :), gr1(:, :), gr0(:, :)
         allocate (dd(n, n, 1), gr1(3, nat), gr0(3, nat))
         dd(:, :, 1) = d
         call eri%gradient(bas, dd, gr1)
         call eri0%gradient(bas, dd, gr0)
         e = maxval(abs(gr1 - gr0))
         print "(a,a,a,es10.3,a,es10.3)", "  ", trim(name), " gradient, unique vs plain: max|g| ", maxval(abs(gr0)), "  err ", e
         call check(error, e <= TOL_REF*max(1.0_dp, maxval(abs(gr0))), &
            trim(name)//": the gradient on the unique engine differs by "// &
                    to_string(e, "es10.3"))
         if (allocated(error)) return
         call check(error, .not. eri%upair, trim(name)//": the pairs were not completed for the gradient")
      end block

      call eri%release()
      call eri0%release()
      call bas%release()
   end subroutine run_case

   subroutine check_tables(eri, bas, label, error)
      !! The unique build's bin-interleaved tables against the ones made from
      !! the primitive pairs of every shell pair.
      type(trc_eri_t), intent(in) :: eri
      type(trc_basis_t), intent(in) :: bas
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error
      type(trc_shell_pairs_t) :: pf
      type(trc_pair_prune_t) :: prune
      type(trc_pair_quartets_t) :: qf
      type(trc_basis_t) :: wb
      real(dp), allocatable :: one(:)
      real(dp) :: dmax
      integer :: nb

      if (eri%dec%active) then
         wb = eri%pb
      else
         wb = bas
      end if
      nb = wb%nshell
      allocate (one(nb))
      one = 1.0_dp
      pf%nbas = nb
      pf%natm = wb%natm
      allocate (pf%sh_l(nb), pf%ao_off(nb), pf%sh_at(nb))
      pf%sh_l = wb%sh_l
      pf%ao_off = wb%sh_ao
      pf%sh_at = wb%sh_at
      ! Component by component: nvfortran 26.5 has an internal compiler
      ! error (search_arr) on this structure constructor with a basis-array argument.
      prune%camp = wb%sh_c
      prune%cut = THRESH*1.0e-3_dp
      call hgp_build_pairs(wb%sh_np, wb%sh_e, wb%sh_c, wb%sh_r, one, prune, pf)
      call bins_build_pair_quartets(eri%bins, pf, qf)
      call check(error, qf%nqp == eri%quartets%nqp, label//": table sizes differ "//to_string(qf%nqp)//" vs "// &
                 to_string(eri%quartets%nqp))
      if (allocated(error)) return
      call check(error, all(qf%n == eri%quartets%n), label//": primitive pair counts of the images differ from the direct ones")
      if (allocated(error)) return
      dmax = max(maxval(abs(qf%z - eri%quartets%z)), maxval(abs(qf%c - eri%quartets%c)), &
                 maxval(abs(qf%a - eri%quartets%a)), &
                 maxval(abs(qf%px - eri%quartets%px)), maxval(abs(qf%py - eri%quartets%py)), &
                 maxval(abs(qf%pz - eri%quartets%pz)), &
                 maxval(abs(qf%ax - eri%quartets%ax)), maxval(abs(qf%ay - eri%quartets%ay)), &
                 maxval(abs(qf%az - eri%quartets%az)), &
                 maxval(abs(qf%abx - eri%quartets%abx)), maxval(abs(qf%aby - eri%quartets%aby)), &
                 maxval(abs(qf%abz - eri%quartets%abz)))
      print "(a,a,a,es10.3)", "  ", label, " bin tables, unique vs full pairs: max deviation ", dmax
      call check(error, dmax <= TOL_TABLE, label//": the rebuilt bin tables deviate by "//to_string(dmax, "es10.3"))
   end subroutine check_tables

   subroutine density(n, d)
      !! A deterministic symmetric matrix; the non-symmetric build needs no
      !! invariance, and the density screen is off by symmetry-blind input.
      integer, intent(in) :: n
      real(dp), intent(out) :: d(n, n)
      integer :: i, j
      do j = 1, n
         do i = 1, n
            d(i, j) = 0.05_dp*(cos(real(2*i + 3*j, dp)) + cos(real(2*j + 3*i, dp)))
         end do
      end do
   end subroutine density

end module check_sym_unique_suite

program check_sym_unique
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_sym_unique_suite, only: collect_sym_unique
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_sym_unique, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_unique
