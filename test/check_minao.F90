!
! The minao initial guess (TRC_GUESS_MINAO / trc_minao's trc_minao_build)
! against PySCF's init_guess_by_minao, plus the algorithm's own closed-form
! checks (the occupation table, and self-projection, where the answer is
! exactly the diagonal occupation matrix whatever the projection code does).
!
module check_minao_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_basis_json, only: trc_basis_from_json
   use trc_sph, only: sph_cart_to_sph
   use trc_minao, only: trc_minao_build, minao_occupation
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_str, ref_read_real, &
                           ref_read_matrix
   implicit none
   private
   public :: collect_minao_tests

   real(dp), parameter :: PYSCF_TOL = 1.0e-9_dp
      !! Same order as mqc's own 1e-9 on the identical algorithm (measured
      !! here: worst |dD| 1.1e-10 on water/cc-pVDZ). terco's cross overlap
      !! goes through one extra Cartesian->spherical fold (sph_cart_to_sph
      !! on the merged basis) that mqc's libcint-native spherical
      !! integrals never pay, and the orthogonaliser is `sym_orthog`'s
      !! eigendecomposition rather than mqc's own, so the two pick up
      !! round-off differently. Far tighter than a guess needs to be
      !! useful.

contains

   subroutine collect_minao_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("occupations follow frac_occ", test_occupations), &
                  new_unittest("self projection is the occupation", test_self_projection), &
                  new_unittest("heavy element is refused", test_heavy_refused), &
                  new_unittest("cartesian occupied-d is refused", test_cartesian_refused), &
                  new_unittest("water and silane match pyscf", test_vs_pyscf) &
                  ]
   end subroutine collect_minao_tests

   !
   ! Tr(D S), the electron count a density carries. `d` is in the SCF's own
   ! space (nao_sph if b%spherical, else nao, same as trc_minao_build's
   ! own output) -- basis_spherical_overlap's fold applies the identity
   ! when the basis is Cartesian throughout (nao_sph == nao then).
   !
   function electrons(b, d) result(n)
      type(trc_basis_t), intent(inout) :: b
      real(dp), intent(in) :: d(:, :)
      real(dp) :: n
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: s(:, :), t(:, :), v(:, :), s_sph(:, :)
      allocate (s(b%nao, b%nao), t(b%nao, b%nao), v(b%nao, b%nao))
      call pl%build(b, 1.0e-12_dp)
      call b%to_device(); call pl%to_device()
      call trc_1e(b, pl, s, t, v)
      call pl%release()
      if (b%spherical) then
         allocate (s_sph(b%nao_sph, b%nao_sph))
         call sph_cart_to_sph(b, s, s_sph)
         n = sum(d*s_sph)
      else
         n = sum(d*s)
      end if
   end function electrons

   subroutine test_occupations(error)
      !! A few rows of PySCF's frac_occ, including the table's oddities --
      !! same cases test_mqc_czt_minao pins, since the table is a verbatim
      !! port of the same NRSRHF_CONFIGURATION data.
      type(error_type), allocatable, intent(out) :: error
      integer :: nd
      real(dp) :: f

      call minao_occupation(8, 1, nd, f)      ! O 2p4: the open shell, nothing doubly occupied
      call check(error, nd == 0 .and. abs(f - 4.0_dp/3.0_dp) < 1.0e-14_dp, "O p")
      if (allocated(error)) return
      call minao_occupation(10, 1, nd, f)     ! Ne: 2p closed, no open column
      call check(error, nd == 1 .and. f == 0.0_dp, "Ne p")
      if (allocated(error)) return
      call minao_occupation(1, 0, nd, f)      ! H 1s1
      call check(error, nd == 0 .and. abs(f - 1.0_dp) < 1.0e-14_dp, "H s")
      if (allocated(error)) return
      call minao_occupation(21, 1, nd, f)     ! Sc: PySCF puts 13 electrons in p
      call check(error, nd == 2 .and. abs(f - 1.0_dp/3.0_dp) < 1.0e-14_dp, "Sc p")
      if (allocated(error)) return
      call minao_occupation(26, 2, nd, f)     ! Fe 3d8 in this table
      call check(error, nd == 0 .and. abs(f - 1.6_dp) < 1.0e-14_dp, "Fe d")
      if (allocated(error)) return
      call minao_occupation(37, 0, nd, f)     ! Rb: outside the data
      call check(error, nd == 0 .and. f == 0.0_dp, "Rb is not covered")
   end subroutine test_occupations

   subroutine test_self_projection(error)
      !! In its own basis the projection is the identity: P = S^+ S is one
      !! (up to the scf_context_orthogonalize floor, which the minimal basis's own
      !! overlap never approaches), so the guess must come back as the
      !! diagonal occupation matrix, and its trace against S as exactly
      !! water's ten neutral electrons.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(error_t) :: err
      real(dp), allocatable :: d(:, :)
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: off
      integer :: i, j, natm_

      call read_xyz('water.xyz', natm_, z, r)
      call trc_basis_from_json('basis_sets/minao.json', natm_, z, r, b, err)
      call check(error,.not. err%has_error(), "water must build in the minimal basis: "//err%get_message())
      if (allocated(error)) return

      call trc_minao_build(b, d, err)
      call check(error,.not. err%has_error(), "minao must build: "//err%get_message())
      if (allocated(error)) return

      ! O: 1s 2s 2p(x3) at 2, 2, 4/3 each; H: 1s at 1.
      call check(error, size(d, 1) == 7, "water's minimal basis is seven functions")
      if (allocated(error)) return
      off = 0.0_dp
      do j = 1, 7
         do i = 1, 7
            if (i /= j) off = max(off, abs(d(i, j)))
         end do
      end do
      call check(error, off < 1.0e-8_dp, "the self-projected density must be diagonal")
      if (allocated(error)) return
      call check(error, abs(d(1, 1) - 2.0_dp) < 1.0e-8_dp, "O 1s")
      if (allocated(error)) return
      call check(error, abs(d(2, 2) - 2.0_dp) < 1.0e-8_dp, "O 2s")
      if (allocated(error)) return
      call check(error, abs(d(3, 3) - 4.0_dp/3.0_dp) < 1.0e-8_dp, "O 2p")
      if (allocated(error)) return
      call check(error, abs(d(6, 6) - 1.0_dp) < 1.0e-8_dp, "H 1s")
      if (allocated(error)) return
      call check(error, abs(electrons(b, d) - 10.0_dp) < 1.0e-6_dp, &
                 "the self-projected density must carry ten electrons")
   end subroutine test_self_projection

   subroutine test_heavy_refused(error)
      !! An element without data is an error, not a silent zero block.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(error_t) :: err
      real(dp), allocatable :: d(:, :)
      real(dp) :: r(3, 1)

      r = 0.0_dp
      call trc_basis_from_json('basis_sets/def2-svp.json', 1, [40], r, b, err)   ! Zr, Z=40
      call check(error,.not. err%has_error(), "Zr must build: "//err%get_message())
      if (allocated(error)) return
      call trc_minao_build(b, d, err)
      call check(error, err%has_error(), "minao must refuse zirconium")
   end subroutine test_heavy_refused

   subroutine test_cartesian_refused(error)
      !! A Cartesian target basis cannot host the five spherical components
      !! minao's occupation table describes once an atom has occupied d.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(error_t) :: err
      real(dp), allocatable :: d(:, :)
      real(dp) :: r(3, 1)

      r = 0.0_dp
      call trc_basis_from_json('basis_sets/def2-svp.json', 1, [22], r, b, err, force_cart=.true.)  ! Ti, Z=22
      call check(error,.not. err%has_error(), "Ti must build: "//err%get_message())
      if (allocated(error)) return
      call check(error, .not. b%spherical, "the probe basis must actually be Cartesian")
      if (allocated(error)) return
      call trc_minao_build(b, d, err)
      call check(error, err%has_error(), "minao must refuse a Cartesian basis with occupied d")
   end subroutine test_cartesian_refused

   subroutine test_vs_pyscf(error)
      !! test/pyscf_refs/minao.ref: water/cc-pVDZ and silane/def2-SVP
      !! against PySCF's init_guess_by_minao, generated by
      !! test/pyscf_refs/generate.py's gen_minao (see its own docstring
      !! for the AO-order/normalisation argument).
      !!
      !! Also dumps "minao_probe.bin" (terco's own spherical S for both
      !! cases), UNCONDITIONALLY, before the reference is even opened --
      !! scf_sph_corr.F90's own convention (never gated on the reference
      !! existing, since generate.py's gen_minao needs this dump to apply
      !! the same per-AO normalisation rescale gen_sph_corr's group B
      !! does for an l>=2 spherical shell: terco's own spherical S
      !! diagonal is not 1 there, see terco_to_pyscf_perm_sph's own
      !! docstring).
      type(error_type), allocatable, intent(out) :: error
      character(len=16), parameter :: MOLNAME(2) = ["water           ", "sih4            "]
      character(len=16), parameter :: BASNAME(2) = ["cc-pvdz.json    ", "def2-svp.json   "]
      type(trc_basis_t) :: b
      type(error_t) :: err
      real(dp), allocatable :: d(:, :), dm_ref(:, :), s_sph(:, :)
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      character(len=64) :: name
      integer :: u, pu, ncase, ic, natm_, nao_ref, i, j
      real(dp) :: nelec_ref, worst

      open (newunit=pu, file="minao_probe.bin", access="stream", form="unformatted", &
            status="replace", action="write")
      do ic = 1, 2
         call read_xyz(trim(MOLNAME(ic))//'.xyz', natm_, z, r)
         call trc_basis_from_json('basis_sets/'//trim(BASNAME(ic)), natm_, z, r, b, err)
         call check(error,.not. err%has_error(), trim(MOLNAME(ic))//" must build: "//err%get_message())
         if (allocated(error)) then
            close (pu)
            return
         end if
         call basis_spherical_overlap(b, s_sph)
         write (pu) int(b%nao_sph, int64)
         write (pu) s_sph
         deallocate (s_sph)
      end do
      close (pu)

      call ref_require("minao.ref", u)
      call ref_read_int(u, "NCASE", ncase)
      call check(error, ncase == 2, "minao.ref must carry two cases")
      if (allocated(error)) return

      do ic = 0, ncase - 1
         call ref_read_str(u, "NAME_"//to_string(ic), name)
         call ref_read_int(u, "NAO_"//to_string(ic), nao_ref)
         call ref_read_real(u, "NELEC_"//to_string(ic), nelec_ref)
         allocate (dm_ref(nao_ref, nao_ref))
         call ref_read_matrix(u, "DM_"//to_string(ic), nao_ref, nao_ref, dm_ref)

         call read_xyz(trim(MOLNAME(ic + 1))//'.xyz', natm_, z, r)
         call trc_basis_from_json('basis_sets/'//trim(BASNAME(ic + 1)), natm_, z, r, b, err)
         call check(error,.not. err%has_error(), trim(name)//" must build: "//err%get_message())
         if (allocated(error)) then
            deallocate (dm_ref)
            call ref_close(u)
            return
         end if
         call check(error, b%nao_sph == nao_ref, trim(name)//": nao mismatch against pyscf")
         if (allocated(error)) then
            deallocate (dm_ref)
            call ref_close(u)
            return
         end if

         call trc_minao_build(b, d, err)
         call check(error,.not. err%has_error(), trim(name)//": minao must build: "//err%get_message())
         if (allocated(error)) then
            deallocate (dm_ref)
            call ref_close(u)
            return
         end if

         worst = 0.0_dp
         do j = 1, nao_ref
            do i = 1, nao_ref
               worst = max(worst, abs(d(i, j) - dm_ref(i, j)))
            end do
         end do
         call check(error, worst < PYSCF_TOL, trim(name)//": density differs from PySCF's, worst |dD| = "// &
                    to_string(worst, "es12.4"))
         if (allocated(error)) then
            deallocate (dm_ref)
            call ref_close(u)
            return
         end if

         call check(error, abs(electrons(b, d) - nelec_ref) < 1.0e-6_dp, &
                    trim(name)//": Tr(D S) differs from PySCF's")
         if (allocated(error)) then
            deallocate (dm_ref)
            call ref_close(u)
            return
         end if
         deallocate (dm_ref)
      end do
      call ref_close(u)
   end subroutine test_vs_pyscf

   !
   ! terco's own overlap, in its spherical AO space (identity fold if the
   ! basis is Cartesian throughout) -- the same S trc_scf_driver itself
   ! builds every geometry.
   !
   subroutine basis_spherical_overlap(b, s_sph)
      type(trc_basis_t), intent(inout) :: b
      real(dp), allocatable, intent(out) :: s_sph(:, :)
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: s(:, :), t(:, :), v(:, :)
      allocate (s(b%nao, b%nao), t(b%nao, b%nao), v(b%nao, b%nao))
      call pl%build(b, 1.0e-12_dp)
      call b%to_device(); call pl%to_device()
      call trc_1e(b, pl, s, t, v)
      call pl%release()
      allocate (s_sph(b%nao_sph, b%nao_sph))
      call sph_cart_to_sph(b, s, s_sph)
   end subroutine basis_spherical_overlap

end module check_minao_suite

program check_minao
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_minao_suite, only: collect_minao_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("check_minao", collect_minao_tests)]
   write (error_unit, fmt) "Testing:", testsuites(1)%name
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop 1
   end if
end program check_minao
