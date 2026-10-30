!
! The folded and enumerated digestions, on a SYMMETRIC density.
!
! They must agree for ANY symmetric density -- that is the whole premise of
! having two of them. `fock` folds the eight-fold permutational symmetry into
! six atomic updates and assumes D(mu,nu) = D(nu,mu); `fock_nosym` enumerates
! the permutations and assumes nothing. Where the assumption holds the two are
! the same operator, so a disagreement here is a bug in one of them and a
! disagreement on an ASYMMETRIC density means nothing until this passes.
!
! Phase 3 recorded this agreeing to 1.24e-14 and then nothing re-ran it. This
! is that control, as a standing test.
!
! IT IS ALSO THE ORACLE FOR THE HYBRID DISPATCH, which is why it takes a basis
! ceiling. `fock` routes each class to an unrolled kernel or, when none exists,
! to the general path over the same segments; `fock_nosym` sends every class
! through the non-symmetric kernel, one_bin_item_nosym. So running this with f functions tests both halves
! of the routing at once: a dropped f class shows up as a disagreement, and so
! does an s..d class digested twice.
!
! Without the ceiling argument it is 6-31G, which stops at d and therefore says
! nothing about any of that.
!
! Registered twice (ctest terco/check_nosym ARGS water.xyz; terco/check_nosym_f
! ARGS water.xyz 3, the uncontracted l=0..3 set that exercises the hybrid
! dispatch). Each invocation gets one argv pair, so this program unit reads
! argv and probes the basis's own angular-momentum ceiling BEFORE any
! test-drive call: when this build's LMAX cannot reach what argv asks for, it
! prints "SKIP: ..." and does `stop 0` right there (ctest matches "SKIP" via
! SKIP_REGULAR_EXPRESSION), never touching run_testsuite -- a clean early
! exit, not a failed check.
!
module check_nosym_suite
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_eri, only: trc_eri_t
   use trc_test_basis, only: read_xyz, build_aux
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_nosym_tests, build_probe_basis, g_xyzfile, g_auxl

   character(len=256) :: g_xyzfile = ''
   integer :: g_auxl = -1
      !! Set from argv by the program unit before run_testsuite is called
      !! (and before the self-skip probe, which reads them too).

contains

   subroutine collect_nosym_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("folded and enumerated Fock digestions agree on a symmetric density", &
                               test_nosym) &
                  ]
   end subroutine collect_nosym_tests

   subroutine build_probe_basis(xyzfile, auxl, bas, err)
      ! An uncontracted l = 0..auxl set, because 6-31G stops at d and the
      ! point of the second argv is to reach past it.
      character(len=*), intent(in) :: xyzfile
      integer, intent(in) :: auxl
      type(trc_basis_t), intent(out) :: bas
      type(error_t), intent(out) :: err
      integer :: nat, nsh, maxnp
      integer, allocatable :: at_z(:), sh_l(:), sh_np(:)
      real(dp), allocatable :: at_r(:, :), zat(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)

      call read_xyz(xyzfile, nat, at_z, at_r)
      if (auxl >= 0) then
         call build_aux(nat, at_z, at_r, auxl, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
         allocate (zat(nat)); zat = real(at_z, dp)
         call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, at_r, maxnp)
      else
         ! 6-31G, through the one basis-construction path.
         call trc_basis_from_json('basis_sets/6-31g.json', nat, at_z, at_r, bas, err)
      end if
   end subroutine build_probe_basis

   subroutine test_nosym(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(error_t) :: err
      real(dp), allocatable :: d(:, :), gs(:, :), gn(:, :)
      real(dp) :: worst, thresh
      integer :: nao

      call build_probe_basis(g_xyzfile, g_auxl, bas, err)
      call check(error,.not. err%has_error(), 'check_nosym: '//err%get_message())
      if (allocated(error)) return

      thresh = 1.0e-10_dp
      call pl%build(bas, thresh)
      nao = bas%nao
      call bas%to_device(); call pl%to_device()
      call eri%build(bas, thresh)

      ! Any symmetric density will do; a deterministic one that is not sparse
      ! exercises more of the digestion than a converged one would.
      allocate (d(nao, nao), gs(nao, nao), gn(nao, nao))
      call fill_density(nao, d)

      call eri%fock(bas, d, gs)
      call eri%fock_nosym(bas, d, gn)
      worst = maxval(abs(gs - gn))

      print '(a,i0,a,i0,a,i0)', '  ', nao, ' functions, shells ', bas%nshell, &
         ', max l ', maxval(bas%sh_l)
      print '(a,es12.4)', '  asymmetry of the density ', maxval(abs(d - transpose(d)))
      print '(a,es24.16)', '  sum folded     ', sum(gs)
      print '(a,es24.16)', '  sum enumerated ', sum(gn)
      print '(a,f18.10)', '  ratio          ', sum(gn)/sum(gs)
      print '(a,es12.4)', '  largest |G|                 ', maxval(abs(gn))
      print '(a,es12.4)', '  worst |folded - enumerated| ', worst
      print '(a,es12.4)', '  relative to largest |G|     ', worst/maxval(abs(gn))
      !
      ! RELATIVE, because the tolerance has to mean the same thing at both sizes
      ! it is run at. 6-31G on water reaches |G| ~ 20; the uncontracted l = 0..3
      ! set reaches ~1e3 over 140 functions, and holding that to the same ABSOLUTE
      ! 1e-10 demands 1e-13 relative -- below what summing that many terms in
      ! double precision returns. Scaling by the largest element asks both for the
      ! same number of correct digits.
      !
      call check(error, worst < 1.0e-11_dp*max(1.0_dp, maxval(abs(gn))), &
                 "folded vs enumerated digestion: worst |diff| relative to largest |G| exceeds tolerance")

      call eri%release(); call pl%release(); call bas%release()
   end subroutine test_nosym

   !
   ! In its own routine, explicit-shape, not inline in the main program. The
   ! one-liner 0.3*cos(real(2*i + 3*j)) + 0.1*real(min(i, j)) crashed
   ! nvfortran 26.5 and 26.9 at -fast -O3 ("bitcast with differing sizes",
   ! "unexpected type for intrinsic function"), and spelling it out inline
   ! only moved the crash ("Errors in Lowering") as the rest of the program
   ! changed. check_rys_fock_vs_pyscf.F90 builds its density the same way,
   ! in a module routine, and compiles. The values are the same to the bit.
   !
   subroutine fill_density(n, dd)
      integer, intent(in) :: n
      real(dp), intent(out) :: dd(n, n)
      integer :: ii, jj, kk
      real(dp) :: xx
      do jj = 1, n
         do ii = 1, n
            xx = real(2*ii + 3*jj, dp)
            kk = ii
            if (jj < ii) kk = jj
            dd(ii, jj) = 0.3_dp*cos(xx) + 0.1_dp*real(kk, dp)
         end do
      end do
      dd = 0.5_dp*(dd + transpose(dd))
   end subroutine fill_density

end module check_nosym_suite

program check_nosym
   use, intrinsic :: iso_fortran_env, only: error_unit
   use pic_types, only: dp
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_basis_t
   use trc_tables, only: LMAX
   use trc_error, only: error_t
   use check_nosym_suite, only: collect_nosym_tests, build_probe_basis, g_xyzfile, g_auxl
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   type(trc_basis_t) :: probe
   type(error_t) :: err
   character(len=256) :: arg
   integer :: stat

   if (command_argument_count() < 1) then
      print '(a)', 'usage: check_nosym <xyz> [aux-lmax]'; stop 1
   end if
   call get_command_argument(1, g_xyzfile)
   g_auxl = -1
   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg); read (arg, *) g_auxl
   end if

   ! The self-skip has to happen before any test-drive call, so this probes
   ! the basis's own ceiling exactly as check_nosym always did, releases it,
   ! and only then decides whether to run the suite at all.
   call build_probe_basis(g_xyzfile, g_auxl, probe, err)
   if (err%has_error()) then
      print '(a)', 'check_nosym: '//err%get_message(); stop 1
   end if
   if (maxval(probe%sh_l) > LMAX) then
      print '(a,i0,a,i0,a)', '  SKIP: basis needs l = ', maxval(probe%sh_l), &
         ' and this build has LMAX = ', LMAX, '.'
      call probe%release()
      stop 0
   end if
   call probe%release()

   stat = 0
   testsuites = [new_testsuite("check_nosym", collect_nosym_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_nosym
