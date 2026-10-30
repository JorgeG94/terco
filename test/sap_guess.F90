module sap_guess_suite
   !! terco's SAP guess: the screening potential it builds, and whether starting
   !! from it is better than starting from GWH.
   !!
   !! TWO CHECKS THAT CAN FAIL INDEPENDENTLY OF THE SCF, each its own
   !! test-drive case below.
   !!
   !! The radial construction is checked against two things it must satisfy and
   !! which nothing else in the guess would reveal:
   !!
   !!   4 pi int r^2 rho dr  =  N_A        the atom's own electron count
   !!   r Vscr(r) -> N_A     as r -> inf   the screened nucleus seen from afar
   !!
   !! Both come out of the same radial quadrature that builds the potential, so
   !! if the collocation, the mesh or either cumulative integral is wrong they
   !! move -- whereas the SCF converges to the same energy from any guess and
   !! would say nothing. The second is the sharper of the two: it is the only
   !! statement about the OUTWARD integral, which is the one accumulated
   !! backwards and the one a plausible-looking implementation gets wrong.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build, trc_atom_ranges
   use trc_sap, only: trc_sap_build, trc_sap_screen_table
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_sap_guess

contains

   subroutine collect_sap_guess(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("radial screening table, per atom", test_radial_table), &
                  new_unittest("SAP as a guess, against GWH and SAD", test_sap_vs_gwh_sad) &
                  ]
   end subroutine collect_sap_guess

   subroutine build_water_631g(bas, pl, smat, tmat, vmat, zint, natm)
      type(trc_basis_t), intent(out) :: bas
      type(trc_pairlist_t), intent(out) :: pl
      real(dp), allocatable, intent(out) :: smat(:, :), tmat(:, :), vmat(:, :)
      integer, allocatable, intent(out) :: zint(:)
      integer, intent(out) :: natm
      real(dp), allocatable :: at_r(:, :)
      type(error_t) :: err
      integer :: dev
      dev = trc_bind_device(0)
      call read_xyz('water.xyz', natm, zint, at_r)
      call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, bas, err)
      if (err%has_error()) error stop "sap_guess: "//err%get_message()
      call bas%to_device()
      call pl%build(bas, 1.0e-12_dp); call pl%to_device()
      allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
      call trc_1e(bas, pl, smat, tmat, vmat)
   end subroutine build_water_631g

   subroutine test_radial_table(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      integer, allocatable :: zint(:), atom_of(:), first(:), last(:)
      integer :: natm, ia, a0, a1, nr
      real(dp), allocatable :: dsad(:, :), rtab(:), vtab(:)
      type(error_t) :: err
      real(dp) :: nel, tail, worst_n, worst_t

      call build_water_631g(bas, pl, smat, tmat, vmat, zint, natm)

      call trc_sad_build(bas, dsad, err)
      if (.not. err%has_error()) call trc_atom_ranges(bas, atom_of, first, last, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) then
         call pl%release(); call bas%release()
         return
      end if
      worst_n = 0.0_dp; worst_t = 0.0_dp
      do ia = 1, natm
         a0 = bas%sh_ao(first(ia))
         a1 = bas%sh_ao(last(ia)) + (bas%sh_l(last(ia)) + 1)*(bas%sh_l(last(ia)) + 2)/2 - 1
         call trc_sap_screen_table(bas, first(ia), last(ia), bas%at_r(:, ia), zint(ia), &
                                   dsad(a0:a1, a0:a1), rtab, vtab, nel, err)
         if (err%has_error()) exit
         nr = size(rtab)
         tail = rtab(nr)*vtab(nr)
         print '(a,i0,a,i0,a,f14.9,a,i0,a,f14.9)', "sap_guess: Z = ", zint(ia), &
            "  want ", zint(ia), " electrons, radial density gives ", nel, &
            "   r*Vscr at r = ", nint(rtab(nr)), " is ", tail
         worst_n = max(worst_n, abs(nel - real(zint(ia), dp)))
         worst_t = max(worst_t, abs(tail - real(zint(ia), dp)))
         deallocate (rtab, vtab)
      end do
      call pl%release(); call bas%release()

      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      ! The radial mesh is a quadrature, so these are quadrature-accurate and not
      ! exact. A part in 1e-6 of an electron is far tighter than a guess needs
      ! and loose enough that the mesh size is not being asserted.
      print '(a,es12.4,a,es12.4)', "sap_guess: worst |N - Z| ", worst_n, "   worst |r*Vscr - Z| ", worst_t
      call check(error, worst_n <= 1.0e-6_dp, "worst |N - Z| = "//to_string(worst_n, "es12.4"))
      if (allocated(error)) return
      call check(error, worst_t <= 1.0e-6_dp, "worst |r*Vscr - Z| = "//to_string(worst_t, "es12.4"))
   end subroutine test_radial_table

   subroutine test_sap_vs_gwh_sad(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      integer, allocatable :: zint(:)
      integer :: natm, nz
      real(dp), allocatable :: dsad(:, :), hsap(:, :), dg3(:, :, :)
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: r_gwh, r_sad, r_sap
      type(error_t) :: err
      real(dp) :: sym_err

      call build_water_631g(bas, pl, smat, tmat, vmat, zint, natm)
      nz = sum(zint)
      opts%conv_energy = 1.0e-10_dp

      call trc_sad_build(bas, dsad, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) then
         call pl%release(); call bas%release()
         return
      end if

      allocate (hsap(bas%nao, bas%nao))
      call trc_sap_build(bas, pl, tmat, vmat, hsap, err, verbose=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) then
         call pl%release(); call bas%release()
         return
      end if
      sym_err = maxval(abs(hsap - transpose(hsap)))

      opts%guess = "gwh"
      call trc_scf_run(bas, nz/2, nz/2, opts, r_gwh)
      allocate (dg3(bas%nao, bas%nao, 1)); dg3(:, :, 1) = dsad
      call trc_scf_run(bas, nz/2, nz/2, opts, r_sad, dguess=dg3)
      call trc_scf_run(bas, nz/2, nz/2, opts, r_sap, hguess=hsap)
      print '(a,f18.12,a,i0,a)', "sap_guess: gwh E = ", r_gwh%energy, " (", r_gwh%iterations, ")"
      print '(a,f18.12,a,i0,a)', "sap_guess: sad E = ", r_sad%energy, " (", r_sad%iterations, ")"
      print '(a,f18.12,a,i0,a)', "sap_guess: sap E = ", r_sap%energy, " (", r_sap%iterations, ")"

      call pl%release(); call bas%release()

      call check(error, sym_err <= 1.0e-10_dp, "H_sap is not symmetric")
      if (allocated(error)) return
      call check(error, r_sap%converged, "SAP-guess SCF did not converge")
      if (allocated(error)) return
      ! One SCF solution from three guesses. Only the energy has to agree; the
      ! iteration counts are reported rather than asserted, because water in
      ! 6-31G converges from almost anything and a count that happened to be
      ! equal here would not be evidence about SAP.
      call check(error, abs(r_sap%energy - r_gwh%energy) <= 1.0e-8_dp, "SAP energy differs from GWH's")
      if (allocated(error)) return
      call check(error, abs(r_sap%energy - r_sad%energy) <= 1.0e-8_dp, "SAP energy differs from SAD's")
   end subroutine test_sap_vs_gwh_sad

end module sap_guess_suite

program sap_guess
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use sap_guess_suite, only: collect_sap_guess
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_sap_guess, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program sap_guess
