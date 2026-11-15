module check_sym_run_suite
   !!
   !! The decisions of trc_sym_run, with no SCF: which group a request means for a
   !! detected geometry (subgroups, groups this build does not support yet), the
   !! electrons-per-irrep text (restricted and alpha/beta), and the occupation
   !! helpers the gradient guard and the orbital report use.
   !!
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_sym_group, only: trc_sym_group_t, sym_group_from_name
   use trc_sym_geom, only: trc_sym_geom_t
   use trc_sym_run, only: sym_run_skeleton_request, SYM_SKEL_AUTO, SYM_SKEL_ON, SYM_SKEL_OFF
   use trc_sym_run, only: trc_sym_request_t, SYM_OFF, SYM_AUTO, SYM_EXPLICIT, sym_run_parse_request, sym_run_choose_group, &
                          sym_run_is_subgroup, sym_run_request_text, sym_run_parse_irrep_nelec, &
                          sym_run_aufbau_matches, sym_run_fixed_occupation, &
                          sym_run_occupied_counts
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
                  new_unittest("occupation helpers", test_occupations), &
                  new_unittest("trc_set_symmetry_skeleton's values: -1 the rule, 0 off, 1 on", test_skeleton_request) &
                  ]
   end subroutine collect_sym_run

   subroutine test_request(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_request_t) :: req
      type(error_t) :: err
      call sym_run_parse_request('off', req, err)
      call check(error, req%mode, SYM_OFF, 'off')
      call sym_run_parse_request(' AUTO ', req, err)
      call check(error, req%mode, SYM_AUTO, 'auto, any case, padded')
      call sym_run_parse_request('c2V', req, err)
      call check(error, req%mode == SYM_EXPLICIT .and. req%group == 'C2v', 'c2V is C2v')
      call sym_run_parse_request('cs', req, err)
      call check(error, req%group, 'Cs', 'cs is Cs')
      call sym_run_parse_request('td', req, err)
      call check(error, err%has_error(), 'td is not an Abelian group of this code')
      if (allocated(error)) return
      ! Symmetry is opt-in: `default` is off; `on`/`auto` is lenient (input
      ! frame, C1 fallback); a group name is strict.
      call sym_run_parse_request('default', req, err)
      call check(error, req%mode == SYM_OFF .and. .not. req%lenient .and. sym_run_request_text(req) == 'off', &
                 'default is off')
      if (allocated(error)) return
      call sym_run_parse_request('on', req, err)
      call check(error, req%mode == SYM_AUTO .and. req%lenient .and. sym_run_request_text(req) == 'auto', 'on is auto')
      if (allocated(error)) return
      call sym_run_parse_request('auto', req, err)
      call check(error, req%mode == SYM_AUTO .and. req%lenient .and. sym_run_request_text(req) == 'auto', &
                 'auto is lenient')
      if (allocated(error)) return
      call sym_run_parse_request('on', req, err)
      call sym_run_parse_request('cs', req, err)
      call check(error, req%mode == SYM_EXPLICIT .and. .not. req%lenient, 'a group after on is strict')
      call sym_run_parse_request('on', req, err)
      call sym_run_parse_request('off', req, err)
      call check(error, req%mode == SYM_OFF .and. .not. req%lenient, 'off after on')
      ! the type itself starts at off: the Fortran entry has no default symmetry
      req = trc_sym_request_t()
      call check(error, req%mode == SYM_OFF .and. .not. req%lenient, 'a fresh request is off')
   end subroutine test_request

   subroutine test_subgroups(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, sym_run_is_subgroup('Cs', 'C2v') .and. sym_run_is_subgroup('Cs', &
         'D2h') .and. sym_run_is_subgroup('C2', 'D2'), &
                 'Cs < C2v, D2h; C2 < D2')
      if (allocated(error)) return
      call check(error, .not. sym_run_is_subgroup('Cs', 'D2') .and. .not. sym_run_is_subgroup('C2v', 'Cs') .and. &
                 .not. sym_run_is_subgroup('Ci', 'C2v') .and. .not. sym_run_is_subgroup('Cs', 'C1'), 'and not the reverse')
      if (allocated(error)) return
      call check(error, sym_run_is_subgroup('C1', 'Ci') .and. sym_run_is_subgroup('Cs', 'Cs'), &
         'C1 is in everything, a group in itself')
   end subroutine test_subgroups

   subroutine test_choose(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(trc_sym_request_t) :: req
      type(error_t) :: err, e2, e3, e4, e5, e6, e7
      character(len=3) :: g

      ! A water-like geometry: C2v found and reported.
      geom%detected = 'C2v'
      geom%group = sym_group_from_name('C2v', err)
      req%mode = SYM_AUTO
      g = sym_run_choose_group(req, geom, err)
      call check(error, g == 'C2v' .and. .not. err%has_error(), 'auto uses the reported group')
      if (allocated(error)) return
      req%mode = SYM_EXPLICIT
      req%group = 'C2v'
      g = sym_run_choose_group(req, geom, e2)
      call check(error, g == 'C2v' .and. .not. e2%has_error(), 'C2v is the detected group and enabled')
      if (allocated(error)) return
      ! A subgroup of the detected group: sym_geom_detect(requested=) reports it, and it is accepted.
      geom%group = sym_group_from_name('Cs', err)
      req%group = 'Cs'
      g = sym_run_choose_group(req, geom, e3)
      call check(error, g == 'Cs' .and. .not. e3%has_error(), 'Cs of a C2v molecule, reported as Cs, is accepted')
      if (allocated(error)) return
      ! A subgroup the detection did not report (no frame rule for the pair, or the
      ! request not passed to sym_geom_detect) is refused rather than run in the wrong frame.
      geom%group = sym_group_from_name('C2v', err)
      g = sym_run_choose_group(req, geom, e6)
      call check(error, e6%has_error() .and. index(e6%get_message(), 'cannot be used for this molecule') > 0, &
                 'Cs of a C2v molecule not reported as Cs is refused: '//e6%get_message())
      if (allocated(error)) return
      req%group = 'Ci'
      g = sym_run_choose_group(req, geom, e7)
      call check(error, e7%has_error() .and. index(e7%get_message(), 'subgroups') > 0, 'Ci is not in C2v: '//e7%get_message())
      if (allocated(error)) return
      ! All eight groups are enabled, so no request is "not yet supported"; a D2h
      ! molecule is accepted as D2h, and D2 (pyscf has no D2 subgroup of D2h) is refused for want of a frame.
      geom%detected = 'D2h'
      geom%group = sym_group_from_name('D2h', err)
      req%group = 'D2h'
      g = sym_run_choose_group(req, geom, e5)
      call check(error, g == 'D2h' .and. .not. e5%has_error(), 'D2h is accepted: '//e5%get_message())
      if (allocated(error)) return
      geom%detected = 'D2h'
      geom%group = sym_group_from_name('D2h', err)
      req%group = 'D2'
      g = sym_run_choose_group(req, geom, e5)
      call check(error, e5%has_error() .and. index(e5%get_message(), 'no frame is defined') > 0, &
                 'D2 of a D2h molecule is refused: '//e5%get_message())
      if (allocated(error)) return
      geom%detected = 'C2v'
      geom%group = sym_group_from_name('C2v', err)
      req%group = 'D2'
      g = sym_run_choose_group(req, geom, e4)
      call check(error, e4%has_error() .and. index(e4%get_message(), 'subgroups') > 0, 'D2 is not in C2v: '//e4%get_message())
      if (allocated(error)) return
      req%group = 'C1'
      g = sym_run_choose_group(req, geom, err)
      call check(error, g == 'C1' .and. .not. err%has_error(), 'C1 is always allowed')
   end subroutine test_choose

   subroutine test_nelec(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_group_t) :: cs
      type(error_t) :: err, e2
      integer, allocatable :: n(:, :)
      cs = sym_group_from_name('Cs', err)
      call sym_run_parse_irrep_nelec("A'=20, A""=6", cs, 1, n, err)
      call check(error, .not. err%has_error() .and. all(n(:, 1) == [10, 3]), 'restricted: electrons halve')
      if (allocated(error)) return
      call sym_run_parse_irrep_nelec("a'=4/3", cs, 2, n, err)
      call check(error, .not. err%has_error() .and. all(n(:, 1) == [4, 0]) .and. all(n(:, 2) == [3, 0]), &
                 'unrestricted: alpha/beta, unlisted irreps empty, case-insensitive')
      if (allocated(error)) return
      call sym_run_parse_irrep_nelec("A'=3", cs, 1, n, err)
      call check(error, err%has_error(), 'an odd count in a restricted run')
      call sym_run_parse_irrep_nelec("A'=4", cs, 2, n, err)
      call check(error, err%has_error(), 'one count in an unrestricted run')
      call sym_run_parse_irrep_nelec("A'=2/2", cs, 1, n, err)
      call check(error, err%has_error(), 'alpha/beta in a restricted run')
      call sym_run_parse_irrep_nelec("A'=2,A'=4", cs, 1, n, err)
      call check(error, err%has_error(), 'an irrep twice')
      call sym_run_parse_irrep_nelec("A1=2", cs, 1, n, e2)
      call check(error, e2%has_error() .and. index(e2%get_message(), "A', A""") > 0, 'a wrong name lists the right ones')
   end subroutine test_nelec

   subroutine test_occupations(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: lab(6, 1)
      logical :: occ(6)
      lab(:, 1) = [1, 2, 1, 1, 2, 1]
      call check(error, sym_run_aufbau_matches(lab, [3], reshape([2, 1], [2, 1])), 'lowest 3 are 2 + 1')
      call check(error, .not. sym_run_aufbau_matches(lab, [3], reshape([3, 0], [2, 1])), 'but not 3 + 0')
      if (allocated(error)) return
      occ = sym_run_fixed_occupation(lab(:, 1), [3, 0])
      call check(error, all(occ .eqv. [.true., .false., .true., .true., .false., .false.]), &
         'fixed 3 + 0 skips the irrep-2 orbitals')
      if (allocated(error)) return
      call check(error, all(sym_run_occupied_counts(lab(:, 1), occ, 2) == [3, 0]), 'counts of the mask')
   end subroutine test_occupations

   subroutine test_skeleton_request(error)
      !! -1 must hand the choice back to the size rule: mapped to "on" (any
      !! nonzero), it would force the skeleton build on the molecules the rule
      !! declines, with no value left to undo it.
      type(error_type), allocatable, intent(out) :: error
      call check(error, sym_run_skeleton_request(-1) == SYM_SKEL_AUTO, "-1 is the size rule")
      if (allocated(error)) return
      call check(error, sym_run_skeleton_request(0) == SYM_SKEL_OFF, "0 is off")
      if (allocated(error)) return
      call check(error, sym_run_skeleton_request(1) == SYM_SKEL_ON, "1 is on")
      if (allocated(error)) return
      call check(error, sym_run_skeleton_request(2) < 0 .and. sym_run_skeleton_request(-2) < 0, "other values are invalid")
   end subroutine test_skeleton_request

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
