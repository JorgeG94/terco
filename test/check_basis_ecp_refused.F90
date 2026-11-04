module check_basis_ecp_refused_suite
   !! The basis reader refuses an element whose basis is written for an
   !! effective core potential (`ecp_potentials` in the BSE JSON): terco has
   !! no ECPs, and reading the valence-only shells all-electron gave HI/def2-SVP
   !! -1996.90 Eh against -297.23 (pyscf, mqc). Elements of the same file
   !! without an ECP still load.
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_basis_ecp_refused

   character(len=*), parameter :: DEF2_SVP = 'basis_sets/def2-svp.json'

contains

   subroutine collect_basis_ecp_refused(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("HI/def2-SVP is refused, naming the ECP", test_hi_refused), &
                  new_unittest("water/def2-SVP still loads", test_water_loads) &
                  ]
   end subroutine collect_basis_ecp_refused

   subroutine test_hi_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(error_t) :: err
      real(dp) :: r(3, 2)
      r = 0.0_dp
      r(3, 2) = 3.04_dp   ! H-I, Bohr
      call trc_basis_from_json(DEF2_SVP, 2, [1, 53], r, bas, err)
      call check(error, err%has_error(), "HI/def2-SVP loaded although iodine's def2-SVP needs an ECP")
      if (allocated(error)) then
         call bas%release()
         return
      end if
      call check(error, index(err%get_message(), "effective core potential") > 0 .and. &
                 index(err%get_message(), "28 core electrons") > 0, &
                 "refusal does not name the ECP: "//err%get_message())
   end subroutine test_hi_refused

   subroutine test_water_loads(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(error_t) :: err
      real(dp) :: r(3, 3)
      r = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.43_dp, 1.11_dp, 0.0_dp, -1.43_dp, 1.11_dp], [3, 3])
      call trc_basis_from_json(DEF2_SVP, 3, [8, 1, 1], r, bas, err)
      call check(error, .not. err%has_error(), "water/def2-SVP: "//err%get_message())
      if (allocated(error)) return
      call check(error, bas%nao_sph == 24, "water/def2-SVP should have 24 spherical functions")
      call bas%release()
   end subroutine test_water_loads

end module check_basis_ecp_refused_suite

program check_basis_ecp_refused
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_basis_ecp_refused_suite, only: collect_basis_ecp_refused
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_basis_ecp_refused, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_basis_ecp_refused
