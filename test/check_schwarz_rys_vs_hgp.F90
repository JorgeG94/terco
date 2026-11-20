module check_schwarz_rys_vs_hgp_suite
   !!
   !! The Rys-quadrature diagonal Schwarz bound (trc_schwarz_rys.F90,
   !! trc_screen.F90's screen_schwarz_bounds_rys) against the HGP whole-block one it
   !! replaces (screen_schwarz_bounds_hgp): both read the SAME hp_* primitive-pair
   !! data and must agree to round-off, because they are two ways of computing
   !! the SAME number -- Q(i,j) = sqrt(max_{ab} |(ij|ij)_{ab,ab}|) -- not two
   !! different approximations to it (trc_screen.F90's own module header).
   !!
   !! Cases: water and formamide in cc-pVDZ (s/p/d, every build), and water in
   !! cc-pVQZ (g on O) when the build's LMAX reaches 4 -- the class
   !! screen_schwarz_bounds_hgp's header says it still has to bound even though the
   !! four-centre Fock path routes it to Rys instead of HGP.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_tables, only: tables_init
   use trc_cart, only: cart_init
   use trc_boys, only: boys_init
   use trc_rys_select, only: rys_select_tables_init
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_hgp, only: hgp_build_pairs, trc_pair_prune_t
   use trc_eri_data, only: trc_shell_pairs_t
   use trc_screen, only: screen_schwarz_bounds_hgp, screen_schwarz_bounds_rys
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_vector
   implicit none
   private
   public :: collect_schwarz_rys_vs_hgp

   ! Per pair, |Q_rys - Q_hgp| <= TOL_REL*Q_hgp + TOL_FLOOR. Measured worst
   ! relative difference: water/cc-pVDZ 1.1e-14, formamide/cc-pVDZ 1.0e-12,
   ! water/cc-pVQZ (g on O) 1.3e-11 -- HGP's recurrences against Rys
   ! quadrature on (gg|gg), round-off. Relative, not a fixed absolute gate:
   ! the tight QZ s/p pairs have Q of order 1-10, where a fixed 1e-12 is
   ! tighter than two methods can agree in double precision. A wrong formula
   ! still fails by many orders of magnitude; the floor only matters for a
   ! Q screened to ~0 by both.
   real(dp), parameter :: TOL_REL = 1.0e-10_dp
   real(dp), parameter :: TOL_FLOOR = 1.0e-15_dp

contains

   subroutine collect_schwarz_rys_vs_hgp(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water_dz", test_water_dz), &
                  new_unittest("formamide_dz", test_formamide_dz), &
                  new_unittest("water_qz", test_water_qz), &
                  new_unittest("water_5z", test_water_5z) &
                  ]
   end subroutine collect_schwarz_rys_vs_hgp

   subroutine test_water_dz(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_q(error, 'water.xyz', 'basis_sets/cc-pvdz.json')
   end subroutine test_water_dz

   subroutine test_formamide_dz(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_q(error, 'formamide.xyz', 'basis_sets/cc-pvdz.json')
   end subroutine test_formamide_dz

   subroutine test_water_qz(error)
      !! g on O: a build below TERCO_LMAX=4 cannot represent this basis at
      !! all (basis_build itself would refuse it), so this case is simply
      !! left out of a build whose own LMAX cannot reach it -- same
      !! reasoning check_binfock.F90 uses per reference system, hoisted to
      !! the whole case rather than tested per shell.
      use trc_tables, only: LMAX
      type(error_type), allocatable, intent(out) :: error
      if (LMAX < 4) then
         call skip_test(error, "SKIP: needs TERCO_LMAX >= 4 (g shells)")
         return
      end if
      call compare_q(error, 'water.xyz', 'basis_sets/cc-pvqz-pyscf17.json')
   end subroutine test_water_qz

   subroutine test_water_5z(error)
      !! h on O: the (hh| pair has nroots = la + lb + 1 = 11 Rys roots, which
      !! only the root tables of 10 and 11 roots (and the 11-way select case of
      !! the diagonal kernel) cover, and HGP's general recurrence reaches it
      !! through Boys order 20. Left out of a build below LMAX=5.
      use trc_tables, only: LMAX
      type(error_type), allocatable, intent(out) :: error
      if (LMAX < 5) then
         call skip_test(error, "SKIP: needs TERCO_LMAX >= 5 (h shells)")
         return
      end if
      call compare_q(error, 'water.xyz', 'basis_sets/cc-pv5z-pyscf17.json', ref="schwarz_water_cc5z.ref")
   end subroutine test_water_5z

   subroutine compare_q(error, xyzfile, basfile, ref)
      !!
      !! Build hp_* exactly as trc_eri.F90's build_structures does (cfac = 1,
      !! since the basis already folded common_fac_sp; camp = sh_c, the
      !! amplitude the kernel's own prune judges by), then compare the two
      !! Schwarz routines on that one primitive-pair set.
      !!
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyzfile, basfile
      character(len=*), intent(in), optional :: ref
      !! A pyscf reference of the exact diagonal bounds (generate.py's
      !! gen_h_schwarz), which both routines are then also held to.

      type(trc_basis_t) :: b
      type(error_t) :: err
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      real(dp), allocatable :: one(:)
      type(trc_shell_pairs_t) :: pairs
      type(trc_pair_prune_t) :: prune
      real(dp), allocatable :: q_hgp(:), q_rys(:)
      integer :: natm, npair, n
      real(dp) :: worst_abs, worst_rel, qmax, d
      logical :: within
      character(len=256) :: msg

      call tables_init(); call cart_init(); call boys_init(); call rys_select_tables_init()

      call read_xyz(xyzfile, natm, zint, at_r)
      ! prune=.false. for the pyscf case: its reference is built from the file's
      ! literal contraction columns, and pruning replaces them by combinations.
      if (present(ref)) then
         call trc_basis_from_json(basfile, natm, zint, at_r, b, err, overrides=trc_basis_json_overrides(prune=.false.))
      else
         call trc_basis_from_json(basfile, natm, zint, at_r, b, err)
      end if
      if (err%has_error()) then
         call check(error, .false., "basis load failed: "//xyzfile//" / "//basfile)
         return
      end if

      allocate (one(b%nshell)); one = 1.0_dp
      pairs%sh_l = b%sh_l
      ! Component by component: nvfortran 26.5 has an internal compiler
      ! error (search_arr) on this structure constructor with a basis-array argument.
      prune%camp = b%sh_c
      prune%cut = 1.0e-30_dp
      call hgp_build_pairs(b%sh_np, b%sh_e, b%sh_c, b%sh_r, one, prune, pairs)

      npair = b%nshell*(b%nshell + 1)/2
      allocate (q_hgp(npair), q_rys(npair))
      call screen_schwarz_bounds_hgp(pairs, q_hgp)
      call screen_schwarz_bounds_rys(pairs, q_rys)

      qmax = maxval(q_hgp)
      worst_abs = 0.0_dp; worst_rel = 0.0_dp
      within = .true.
      do n = 1, npair
         d = abs(q_rys(n) - q_hgp(n))
         worst_abs = max(worst_abs, d)
         if (d > TOL_REL*abs(q_hgp(n)) + TOL_FLOOR) within = .false.
         ! Reported relative difference: only above the noise floor of the
         ! largest Q, where it means something.
         if (q_hgp(n) > qmax*1.0e-10_dp) then
            worst_rel = max(worst_rel, d/q_hgp(n))
         end if
      end do

      write (msg, '(a,es10.2,a,es10.2,a,i0,a)') &
         '  worst |dQ| ', worst_abs, '   worst |dQ|/|Q| ', worst_rel, '   (', npair, ' pairs)'
      print '(a)', trim(msg)

      if (.not. present(ref)) then
         call check(error, within, &
                    message="screen_schwarz_bounds_rys disagrees with screen_schwarz_bounds_hgp past round-off "// &
                    "for "//xyzfile//" / "//basfile)
         return
      end if
      block
         real(dp), allocatable :: q_ref(:)
         real(dp) :: worst_rys, worst_hgp
         integer :: u
         logical :: ok
         call ref_open(ref, u, ok)
         call check(error, ok, "missing pyscf_refs/"//ref//" -- run test/pyscf_refs/generate.py schwarz_h")
         if (allocated(error)) return
         allocate (q_ref(npair))
         call ref_read_vector(u, "Q", npair, q_ref)
         call ref_close(u)
         worst_rys = maxval(abs(q_rys - q_ref)/(abs(q_ref) + 1.0e-12_dp*maxval(q_ref)))
         worst_hgp = maxval(abs(q_hgp - q_ref)/(abs(q_ref) + 1.0e-12_dp*maxval(q_ref)))
         write (msg, '(a,es10.2,a,es10.2)') '  against pyscf: worst relative |dQ|, Rys ', worst_rys, ', HGP ', worst_hgp
         print '(a)', trim(msg)
         call check(error, worst_rys <= TOL_REL, "the Rys Schwarz bound differs from pyscf's exact diagonal by "// &
                    to_string(worst_rys, "es10.3"))
         if (allocated(error)) return
         call check(error, worst_hgp <= TOL_REL, "the HGP Schwarz bound differs from pyscf's exact diagonal by "// &
                    to_string(worst_hgp, "es10.3"))
         if (allocated(error)) return
      end block
      call check(error, within, &
                 message="screen_schwarz_bounds_rys disagrees with screen_schwarz_bounds_hgp past round-off "// &
                 "for "//xyzfile//" / "//basfile)
   end subroutine compare_q

end module check_schwarz_rys_vs_hgp_suite

program check_schwarz_rys_vs_hgp
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, run_selected
   use check_schwarz_rys_vs_hgp_suite, only: collect_schwarz_rys_vs_hgp
   implicit none
   integer :: stat
   character(len=64) :: name
   stat = 0
   ! One case per ctest entry (the name is the first argument): a build that
   ! cannot reach g or h skips those entries without hiding the others.
   if (command_argument_count() >= 1) then
      call get_command_argument(1, name)
      call run_selected(collect_schwarz_rys_vs_hgp, trim(name), error_unit, stat)
   else
      call run_testsuite(collect_schwarz_rys_vs_hgp, error_unit, stat, parallel=.false.)
   end if
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_schwarz_rys_vs_hgp
