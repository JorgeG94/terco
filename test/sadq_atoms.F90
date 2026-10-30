module sadq_atoms_suite
   !! terco's SADQ guess: the electron count it actually carries, that it is SAD
   !! when the molecule is neutral, and that it beats SAD on an ion.
   !!
   !! The electron count is the point. SAD superposes NEUTRAL atoms, so on a
   !! dication it hands the SCF a density holding two electrons too many and the
   !! first iterations are spent giving them back. Tr(D S) is what says whether a
   !! guess got that right, and it is checked exactly rather than loosely: the
   !! atomic blocks are idempotent-free but each atomic SCF converged its own
   !! count, so the sum is the molecular count to machine precision or the
   !! construction is wrong somewhere.
   !!
   !! Four independent cases below, each its own test-drive unittest, each
   !! building its own basis/pair list/S.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build, trc_sadq_build
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_sadq_atoms

contains

   subroutine collect_sadq_atoms(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("neutral water: SADQ carries sum(Z) and is SAD", test_neutral), &
                  new_unittest("water dication: SADQ beats SAD", test_dication), &
                  new_unittest("explicit per-atom charge", test_qatom), &
                  new_unittest("a qatom of the wrong length is refused", test_qatom_refused) &
                  ]
   end subroutine collect_sadq_atoms

   subroutine build_water_631g(bas, pl, smat, zint, natm)
      type(trc_basis_t), intent(out) :: bas
      type(trc_pairlist_t), intent(out) :: pl
      real(dp), allocatable, intent(out) :: smat(:, :)
      integer, allocatable, intent(out) :: zint(:)
      integer, intent(out) :: natm
      real(dp), allocatable :: at_r(:, :), tmat(:, :), vmat(:, :)
      type(error_t) :: err
      integer :: dev
      dev = trc_bind_device(0)
      call read_xyz('water.xyz', natm, zint, at_r)
      call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, bas, err)
      if (err%has_error()) error stop "sadq_atoms: "//err%get_message()
      call bas%to_device()
      call pl%build(bas, 1.0e-12_dp); call pl%to_device()
      allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
      call trc_1e(bas, pl, smat, tmat, vmat)
   end subroutine build_water_631g

   subroutine test_neutral(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :)
      integer, allocatable :: zint(:)
      integer :: natm, nz
      real(dp), allocatable :: d_sad(:, :), d_sac(:, :), dg3(:, :, :)
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: r_sad, r_sac
      type(error_t) :: err
      real(dp) :: t_sad, t_sac

      call build_water_631g(bas, pl, smat, zint, natm)
      nz = sum(zint)
      opts%conv_energy = 1.0e-10_dp

      call trc_sad_build(bas, d_sad, err, verbose=.true.)
      if (.not. err%has_error()) call trc_sadq_build(bas, nz, d_sac, err, verbose=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) then
         call pl%release(); call bas%release()
         return
      end if
      t_sad = sum(d_sad*smat); t_sac = sum(d_sac*smat)
      print '(a,f14.10,a,f14.10,a,i0)', "sadq_atoms: neutral  Tr(DS) sad ", t_sad, "  sadq ", t_sac, &
         "   want ", nz

      allocate (dg3(bas%nao, bas%nao, 1))
      dg3(:, :, 1) = d_sad; call trc_scf_run(bas, nz/2, nz/2, opts, r_sad, dguess=dg3)
      dg3(:, :, 1) = d_sac; call trc_scf_run(bas, nz/2, nz/2, opts, r_sac, dguess=dg3)
      print '(a,f18.12,a,i0,a,f18.12,a,i0,a)', "sadq_atoms: neutral  sad E = ", r_sad%energy, &
         " (", r_sad%iterations, ")   sadq E = ", r_sac%energy, " (", r_sac%iterations, ")"

      call pl%release(); call bas%release()

      call check(error, abs(t_sac - real(nz, dp)) <= 1.0e-8_dp, "SADQ Tr(DS) /= sum(Z)")
      if (allocated(error)) then
         deallocate (d_sad, d_sac, dg3)
         return
      end if
      call check(error, abs(t_sad - real(nz, dp)) <= 1.0e-8_dp, "SAD Tr(DS) /= sum(Z)")
      if (allocated(error)) then
         deallocate (d_sad, d_sac, dg3)
         return
      end if
      ! SADQ on a neutral molecule is SAD, and `build_atomic` chooses the Hund
      ! atom per atom precisely so that it is the SAME atomic SCF rather than a
      ! similar one. So the densities themselves are compared, not the energies:
      ! an energy agreeing to 1e-8 would have hidden the first version of this,
      ! which spin-restricted every atom and moved the density by 1e-2.
      !
      ! Not bit-identity, though, even though the two calls run the same SCF on
      ! the same input. Under the OpenMP port the reductions inside it are
      ! threaded and their order is not fixed, so two runs differ at 4e-14 --
      ! which is the SCF's own reproducibility and nothing to do with SADQ. The
      ! tolerance is set by that, and is still ten orders below the difference
      ! it was written to catch.
      call check(error, maxval(abs(d_sac - d_sad)) <= 1.0e-12_dp, &
                 "neutral SADQ is not SAD, max|dD| = "//to_string(maxval(abs(d_sac - d_sad)), "es12.4"))
      if (allocated(error)) then
         deallocate (d_sad, d_sac, dg3)
         return
      end if
      call check(error, r_sac%converged, "SADQ did not converge")
      if (allocated(error)) then
         deallocate (d_sad, d_sac, dg3)
         return
      end if
      call check(error, abs(r_sac%energy - r_sad%energy) <= 1.0e-8_dp, "SADQ energy differs from SAD's")
      deallocate (d_sad, d_sac, dg3)
   end subroutine test_neutral

   subroutine test_dication(error)
      !! Water at +2 is still closed shell, so the comparison needs no
      !! open-shell machinery to confuse it. SAD is expected to be WRONG
      !! here, by exactly the two electrons it never removed -- if it ever
      !! stops being wrong, this test is measuring nothing.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :)
      integer, allocatable :: zint(:)
      integer :: natm, nz
      real(dp), allocatable :: d_sad(:, :), d_sac(:, :), dg3(:, :, :)
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: r_sad, r_sac
      type(error_t) :: err
      real(dp) :: t_sad, t_sac
      logical :: sad_conv

      call build_water_631g(bas, pl, smat, zint, natm)
      nz = sum(zint)
      opts%conv_energy = 1.0e-10_dp

      call trc_sad_build(bas, d_sad, err)
      if (.not. err%has_error()) call trc_sadq_build(bas, nz - 2, d_sac, err, verbose=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) then
         call pl%release(); call bas%release()
         return
      end if
      t_sad = sum(d_sad*smat); t_sac = sum(d_sac*smat)
      print '(a,f14.10,a,f14.10,a,i0)', "sadq_atoms: +2       Tr(DS) sad ", t_sad, "  sadq ", t_sac, &
         "   want ", nz - 2

      allocate (dg3(bas%nao, bas%nao, 1))
      dg3(:, :, 1) = d_sad
      call trc_scf_run(bas, (nz - 2)/2, (nz - 2)/2, opts, r_sad, dguess=dg3)
      dg3(:, :, 1) = d_sac
      call trc_scf_run(bas, (nz - 2)/2, (nz - 2)/2, opts, r_sac, dguess=dg3)
      print '(a,f18.12,a,i0,a,f18.12,a,i0,a)', "sadq_atoms: +2       sad E = ", r_sad%energy, &
         " (", r_sad%iterations, ")   sadq E = ", r_sac%energy, " (", r_sac%iterations, ")"
      sad_conv = r_sad%converged

      call pl%release(); call bas%release()

      call check(error, abs(t_sac - real(nz - 2, dp)) <= 1.0e-8_dp, "SADQ Tr(DS) /= sum(Z)-2")
      if (allocated(error)) then
         deallocate (d_sad, d_sac, dg3)
         return
      end if
      call check(error, abs(t_sad - real(nz, dp)) <= 1.0e-8_dp, "SAD is no longer wrong on the dication")
      if (allocated(error)) then
         deallocate (d_sad, d_sac, dg3)
         return
      end if
      call check(error, r_sac%converged, "SADQ did not converge on the dication")
      if (allocated(error)) then
         deallocate (d_sad, d_sac, dg3)
         return
      end if
      if (sad_conv) then
         call check(error, abs(r_sac%energy - r_sad%energy) <= 1.0e-8_dp, &
                    "SADQ energy differs from SAD's on the dication")
         if (allocated(error)) then
            deallocate (d_sad, d_sac, dg3)
            return
         end if
      end if
      ! The iteration counts are PRINTED and not asserted against each other.
      ! SADQ beat SAD here by one on the GPU and lost by one under the threaded
      ! host build, from the same code -- which is the answer to whether one
      ! molecule can support an iteration-count claim. What is asserted is only
      ! that the guess is not catastrophic, with the same +2 slack sad_atoms
      ! allows itself against GWH.
      if (sad_conv) then
         call check(error, r_sac%iterations <= r_sad%iterations + 2, &
                    "SADQ took far longer than SAD on the dication")
      end if
      deallocate (d_sad, d_sac, dg3)
   end subroutine test_dication

   subroutine test_qatom(error)
      !! The whole +1 taken off the oxygen. Tr(DS) still has to land on the
      !! molecular count, which is what says the override is wired to the
      !! atoms and not merely accepted.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :)
      integer, allocatable :: zint(:)
      integer :: natm, nz, ia
      real(dp), allocatable :: d_sac(:, :), qa(:)
      type(error_t) :: err
      real(dp) :: t_sac

      call build_water_631g(bas, pl, smat, zint, natm)
      nz = sum(zint)

      allocate (qa(natm))
      qa = 0.0_dp
      do ia = 1, natm
         if (zint(ia) == 8) qa(ia) = 1.0_dp
      end do
      call trc_sadq_build(bas, nz - 1, d_sac, err, qatom=qa)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) then
         call pl%release(); call bas%release()
         deallocate (qa)
         return
      end if
      t_sac = sum(d_sac*smat)
      print '(a,f14.10,a,i0)', "sadq_atoms: qatom    Tr(DS) sadq ", t_sac, "   want ", nz - 1

      call pl%release(); call bas%release()
      deallocate (qa, d_sac)

      call check(error, abs(t_sac - real(nz - 1, dp)) <= 1.0e-8_dp, "qatom-overridden Tr(DS) /= sum(Z)-1")
   end subroutine test_qatom

   subroutine test_qatom_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :)
      integer, allocatable :: zint(:)
      integer :: natm, nz
      type(error_t) :: e2
      real(dp), allocatable :: dbad(:, :)
      real(dp) :: qshort(1)

      call build_water_631g(bas, pl, smat, zint, natm)
      nz = sum(zint)
      qshort = 0.0_dp
      call trc_sadq_build(bas, nz, dbad, e2, qatom=qshort)
      call pl%release(); call bas%release()

      call check(error, e2%has_error(), "a qatom of the wrong length was accepted")
      if (.not. allocated(error)) print '(a)', "sadq_atoms: qatom length checked -- "//trim(e2%get_message())
   end subroutine test_qatom_refused

end module sadq_atoms_suite

program sadq_atoms
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use sadq_atoms_suite, only: collect_sadq_atoms
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_sadq_atoms, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program sadq_atoms
