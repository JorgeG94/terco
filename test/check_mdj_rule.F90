module check_mdj_rule_suite
   !! The J-engine choice of a resident Fock build (trc_eri's eri_mdj_size_rule,
   !! read through mdj_requested): an extended cc-pVDZ chain takes MD J and
   !! a compact water cluster keeps the bins' shared J+K pass, and the two
   !! engines give the same G on the same density, so the choice moves
   !! the time and not the numbers.
   !!
   !! The rule's verdict is what `build` stores in eri%mdj_extended; the
   !! engines are forced through trc_fock_options_t's j_engine, which a
   !! TRC_MDJ in the environment would override (and the equality then
   !! compares an engine with itself, which still passes).
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t, trc_fock_options_t, TRC_J_ENGINE_BINS, TRC_J_ENGINE_MDJ
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_mdj_rule

   real(dp), parameter :: THRESH = 1.0e-12_dp
   real(dp), parameter :: TOL = 1.0e-11_dp
   !! max |G_mdj - G_bins| relative to max |G|; measured 9.3e-13 (gfortran 15, host)
   real(dp), parameter :: DECAY = 0.3_dp
   !! of the synthetic density with the distance between shell centres, bohr^-1

contains

   subroutine collect_mdj_rule(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("gly10 cc-pVDZ (chain): MD J", test_chain), &
                  new_unittest("w5 cc-pVDZ (compact): bins", test_cluster), &
                  new_unittest("gly3 cc-pVDZ: MD J and bins give the same G", test_same_g) &
                  ]
   end subroutine collect_mdj_rule

   subroutine verdict(error, xyz, basis, expect)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis
      logical, intent(in) :: expect
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      logical :: got

      call read_xyz(xyz, nat, at_z, r)
      call trc_basis_from_json(basis, nat, at_z, r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call eri%build(bas, 1.0e-11_dp)
      got = eri%mdj_extended
      call eri%release()
      call bas%release()
      call check(error, got .eqv. expect, xyz//" "//basis//": eri_mdj_size_rule chose the wrong J engine")
   end subroutine verdict

   subroutine test_chain(error)
      type(error_type), allocatable, intent(out) :: error
      call verdict(error, "gly10.xyz", "basis_sets/cc-pvdz.json", .true.)
   end subroutine test_chain

   subroutine test_cluster(error)
      type(error_type), allocatable, intent(out) :: error
      call verdict(error, "w5.xyz", "basis_sets/cc-pvdz.json", .false.)
   end subroutine test_cluster

   subroutine test_same_g(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      integer :: nat, nao, i, j, si, sj, ni, nj
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d(:, :), gb(:, :), gm(:, :), work(:)
      real(dp) :: worst, scale, dist
      character(len=160) :: msg

      call read_xyz("gly3.xyz", nat, at_z, r)
      call trc_basis_from_json("basis_sets/cc-pvdz.json", nat, at_z, r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao
      ! A symmetric density that decays with distance, so both engines'
      ! density screens have something to drop.
      allocate (d(nao, nao), gb(nao, nao), gm(nao, nao))
      do sj = 1, bas%nshell
         nj = (bas%sh_l(sj) + 1)*(bas%sh_l(sj) + 2)/2
         do si = 1, bas%nshell
            ni = (bas%sh_l(si) + 1)*(bas%sh_l(si) + 2)/2
            dist = norm2(bas%sh_r(:, si) - bas%sh_r(:, sj))
            do j = bas%sh_ao(sj), bas%sh_ao(sj) + nj - 1
               do i = bas%sh_ao(si), bas%sh_ao(si) + ni - 1
                  d(i, j) = 0.4_dp*exp(-DECAY*dist)*(1.0_dp + 0.1_dp*cos(real(i + j, dp)))
               end do
            end do
         end do
      end do
      d = 0.5_dp*(d + transpose(d))

      call bas%to_device()
      call eri%build(bas, THRESH)
      allocate (work(eri%fock_work_size(1)))
      !$acc enter data copyin(d) create(gb, gm, work)
      call eri%fock_resident_lent(bas, d, gb, trc_fock_options_t(j_engine=TRC_J_ENGINE_BINS), work)
      call eri%fock_resident_lent(bas, d, gm, trc_fock_options_t(j_engine=TRC_J_ENGINE_MDJ), work)
      !$acc update self(gb, gm)
      !$acc exit data delete(d, gb, gm, work)
      call eri%release()
      call bas%release()

      scale = maxval(abs(gb))
      worst = maxval(abs(gm - gb))/scale
      write (msg, '(a,es9.2,a,es9.2)') "max |G_mdj - G_bins|/max|G| = ", worst, ", max|G| = ", scale
      call check(error, worst <= TOL, trim(msg))
   end subroutine test_same_g

end module check_mdj_rule_suite

program check_mdj_rule
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_mdj_rule_suite, only: collect_mdj_rule
   implicit none
   integer :: stat
   stat = 0
   ! Serial: every case builds device-resident engines.
   call run_testsuite(collect_mdj_rule, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, '(i0,a)') stat, " test(s) failed"
      error stop 1
   end if
end program check_mdj_rule
