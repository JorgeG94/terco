!
! The decisions of trc_sym_run, with no SCF: which group a request means for a
! detected geometry (subgroups, groups this build does not support yet), the
! electrons-per-irrep text (restricted and alpha/beta), and the occupation
! helpers the gradient guard and the orbital report use.
!
module check_sym_run_suite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_sym_group, only: trc_sym_group_t, sym_group_from_name
   use trc_sym_geom, only: trc_sym_geom_t
   use trc_sym_run, only: trc_sym_request_t, SYM_OFF, SYM_AUTO, SYM_EXPLICIT, sym_parse_request, sym_choose_group, &
                          sym_is_subgroup, sym_parse_irrep_nelec, sym_aufbau_matches, sym_fixed_occupation, &
                          sym_occupied_counts
   implicit none
   private
   public :: collect_sym_run

contains

   subroutine collect_sym_run(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("requests: off, auto, groups in any case, nonsense", test_request), &
                  new_unittest("subgroup table", test_subgroups), &
                  new_unittest("choosing the group for a detected geometry", test_choose), &
                  new_unittest("irrep_nelec text: restricted, unrestricted, errors", test_nelec), &
                  new_unittest("occupation helpers", test_occupations) &
                  ]
   end subroutine collect_sym_run

   subroutine test_request(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_request_t) :: req
      type(error_t) :: err
      call sym_parse_request('off', req, err)
      call check(error, req%mode, SYM_OFF, 'off')
      call sym_parse_request(' AUTO ', req, err)
      call check(error, req%mode, SYM_AUTO, 'auto, any case, padded')
      call sym_parse_request('c2V', req, err)
      call check(error, req%mode == SYM_EXPLICIT .and. req%group == 'C2v', 'c2V is C2v')
      call sym_parse_request('cs', req, err)
      call check(error, req%group, 'Cs', 'cs is Cs')
      call sym_parse_request('td', req, err)
      call check(error, err%has_error(), 'td is not an Abelian group of this code')
   end subroutine test_request

   subroutine test_subgroups(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, sym_is_subgroup('Cs', 'C2v') .and. sym_is_subgroup('Cs', 'D2h') .and. sym_is_subgroup('C2', 'D2'), &
                 'Cs < C2v, D2h; C2 < D2')
      if (allocated(error)) return
      call check(error, .not. sym_is_subgroup('Cs', 'D2') .and. .not. sym_is_subgroup('C2v', 'Cs') .and. &
                 .not. sym_is_subgroup('Ci', 'C2v') .and. .not. sym_is_subgroup('Cs', 'C1'), 'and not the reverse')
      if (allocated(error)) return
      call check(error, sym_is_subgroup('C1', 'Ci') .and. sym_is_subgroup('Cs', 'Cs'), 'C1 is in everything, a group in itself')
   end subroutine test_subgroups

   subroutine test_choose(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(trc_sym_request_t) :: req
      type(error_t) :: err, e2, e3, e4, e5
      character(len=3) :: g

      ! A water-like geometry: C2v found and reported (C2v is enabled since S8).
      geom%detected = 'C2v'
      geom%group = sym_group_from_name('C2v', err)
      req%mode = SYM_AUTO
      g = sym_choose_group(req, geom, err)
      call check(error, g == 'C2v' .and. .not. err%has_error(), 'auto uses the reported group')
      if (allocated(error)) return
      req%mode = SYM_EXPLICIT
      req%group = 'C2v'
      g = sym_choose_group(req, geom, e2)
      call check(error, g == 'C2v' .and. .not. e2%has_error(), 'C2v is the detected group and enabled')
      if (allocated(error)) return
      ! A subgroup below the reported group (here Cs of C2v) has no frame of its own
      ! yet: sym_detect reports one group, so it is refused rather than run in the wrong frame.
      req%group = 'Cs'
      g = sym_choose_group(req, geom, e3)
      call check(error, e3%has_error() .and. index(e3%get_message(), 'cannot be used for this molecule') > 0, &
                 'Cs of a C2v molecule is refused: '//e3%get_message())
      if (allocated(error)) return
      ! A detected group that is not enabled (D2h, S11) is refused as not yet supported.
      geom%detected = 'D2h'
      geom%group = sym_group_from_name('Cs', err)
      req%group = 'D2h'
      g = sym_choose_group(req, geom, e5)
      call check(error, e5%has_error() .and. index(e5%get_message(), 'not yet supported') > 0, &
                 'D2h is not supported yet: '//e5%get_message())
      if (allocated(error)) return
      geom%detected = 'C2v'
      geom%group = sym_group_from_name('C2v', err)
      req%group = 'D2'
      g = sym_choose_group(req, geom, e4)
      call check(error, e4%has_error() .and. index(e4%get_message(), 'subgroups') > 0, 'D2 is not in C2v: '//e4%get_message())
      if (allocated(error)) return
      req%group = 'C1'
      g = sym_choose_group(req, geom, err)
      call check(error, g == 'C1' .and. .not. err%has_error(), 'C1 is always allowed')
   end subroutine test_choose

   subroutine test_nelec(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: cs
      type(error_t) :: err, e2
      integer, allocatable :: n(:, :)
      cs = sym_group_from_name('Cs', err)
      call sym_parse_irrep_nelec("A'=20, A""=6", cs, 1, n, err)
      call check(error, .not. err%has_error() .and. all(n(:, 1) == [10, 3]), 'restricted: electrons halve')
      if (allocated(error)) return
      call sym_parse_irrep_nelec("a'=4/3", cs, 2, n, err)
      call check(error, .not. err%has_error() .and. all(n(:, 1) == [4, 0]) .and. all(n(:, 2) == [3, 0]), &
                 'unrestricted: alpha/beta, unlisted irreps empty, case-insensitive')
      if (allocated(error)) return
      call sym_parse_irrep_nelec("A'=3", cs, 1, n, err)
      call check(error, err%has_error(), 'an odd count in a restricted run')
      call sym_parse_irrep_nelec("A'=4", cs, 2, n, err)
      call check(error, err%has_error(), 'one count in an unrestricted run')
      call sym_parse_irrep_nelec("A'=2/2", cs, 1, n, err)
      call check(error, err%has_error(), 'alpha/beta in a restricted run')
      call sym_parse_irrep_nelec("A'=2,A'=4", cs, 1, n, err)
      call check(error, err%has_error(), 'an irrep twice')
      call sym_parse_irrep_nelec("A1=2", cs, 1, n, e2)
      call check(error, e2%has_error() .and. index(e2%get_message(), "A', A""") > 0, 'a wrong name lists the right ones')
   end subroutine test_nelec

   subroutine test_occupations(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: lab(6, 1)
      logical :: occ(6)
      lab(:, 1) = [1, 2, 1, 1, 2, 1]
      call check(error, sym_aufbau_matches(lab, [3], reshape([2, 1], [2, 1])), 'lowest 3 are 2 + 1')
      call check(error, .not. sym_aufbau_matches(lab, [3], reshape([3, 0], [2, 1])), 'but not 3 + 0')
      if (allocated(error)) return
      occ = sym_fixed_occupation(lab(:, 1), [3, 0])
      call check(error, all(occ .eqv. [.true., .false., .true., .true., .false., .false.]), 'fixed 3 + 0 skips the irrep-2 orbitals')
      if (allocated(error)) return
      call check(error, all(sym_occupied_counts(lab(:, 1), occ, 2) == [3, 0]), 'counts of the mask')
   end subroutine test_occupations

end module check_sym_run_suite

program check_sym_run
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_sym_run_suite, only: collect_sym_run
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_sym_run, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_run
