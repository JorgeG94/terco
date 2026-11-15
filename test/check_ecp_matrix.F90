module check_ecp_matrix_suite
   !!
   !! The scalar ECP one-electron matrix (trc_ecp_matrix) against pyscf's
   !! intor('ECPscalar'), fed the SAME BSE JSON terco reads (def2-svp.json) --
   !! not pyscf's own bundled ECP tables, which can differ (ECP_DESIGN.md's own
   !! note, and mqc's validation approach for the same reason).
   !!
   !! Five molecules, one per invocation (test/CMakeLists.txt's foreach), the
   !! case name selecting a hardcoded geometry below:
   !!
   !!   hi_svp   HI, bonded              -- the bug this feature exists to fix
   !!            (test/basis_sets/def2-svp.json's own ecp_potentials; an
   !!            all-electron read of I gave -1996.90 Eh against -297.23)
   !!   ih_far   I...H at 25 bohr        -- large separation (ECP_DESIGN.md)
   !!   xe       Xe, a single ECP atom with no orbital partner nearby
   !!   csh      CsH, bonded
   !!   snh4     SnH4, tetrahedral, one ECP atom with four light partners
   !!   i2       I2 at 5.04 bohr         -- two ECP atoms: a shell on one against
   !!            the other's channels. The cases above all pair the ECP atom
   !!            with hydrogen, so they never reach the semi-local blocks that
   !!            overflowed the type-2 scratch or that never show two equal
   !!            quadrature levels (#194 review); this one does both.
   !!   ch3i     CH3I                    -- the same through a carbon partner
   !!
   !! Each case also checks the SPHERICAL output: the same molecule read
   !! without force_cart must give exactly T^T U_cart T, folded here on the
   !! host from the Cartesian matrix above. That runs trc_ecp_matrix's own
   !! spherical transform, which on the GPU build needs its own data mapping.
   !!
   !! The reference is lexicographic (pyscf's order); it is read through
   !! cart_perm_b2l, as check_1e does, so a TERCO_CART_ORDER=gamess build
   !! compares like with like.
   !!
   !! Cartesian throughout (force_cart=.true. here, cart=True in pyscf): this
   !! build's raw Cartesian AO convention already matches pyscf's int1e_ovlp_cart
   !! with no rescale for an UNCONTRACTED probe (check_mult.F90's own finding),
   !! but def2-svp's shells are contracted, and gen_rys_fock's own history shows
   !! a contracted shell's normalisation can still differ measurably -- so this
   !! reads terco's own overlap (trc_1e's S) alongside U and lets
   !! test/pyscf_refs/generate.py's gen_ecp_ints rescale pyscf's U by the SAME
   !! ao_norm_from_overlap factor every other AO-level pyscf comparison here
   !! uses, rather than assume Cartesian needs none.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_ecp_types, only: trc_ecp_t
   use trc_ecp_json, only: trc_ecp_from_json
   use trc_ecp, only: trc_ecp_matrix
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_report
   use trc_cart_perm, only: cart_perm_b2l
   use trc_1e_kernels, only: ONE_E_LMAX
   use trc_sph_tables, only: SPH_T
   implicit none
   private
   public :: collect_ecp_matrix_tests

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: TOL = 1.0e-9_dp
   character(len=*), parameter :: DEF2_SVP = 'basis_sets/def2-svp.json'

contains

   subroutine collect_ecp_matrix_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("ECP matrix vs pyscf", test_case)]
   end subroutine collect_ecp_matrix_tests

   subroutine case_geometry(case_name, nat, elements, r)
      !!
      !! The five hardcoded geometries, Bohr. `case_name` picks one; `generate.py`'s
      !! gen_ecp_ints carries the SAME literal numbers (this module's own header
      !! on why: no xyz file, so the two have to agree by inspection rather than
      !! by a shared reader).
      !!
      character(len=*), intent(in) :: case_name
      integer, intent(out) :: nat
      integer, allocatable, intent(out) :: elements(:)
      real(dp), allocatable, intent(out) :: r(:, :)

      select case (trim(case_name))
      case ("hi_svp")
         nat = 2
         allocate (elements(nat), r(3, nat))
         elements = [1, 53]
         r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
         r(:, 2) = [0.0_dp, 0.0_dp, 3.04_dp]
      case ("ih_far")
         nat = 2
         allocate (elements(nat), r(3, nat))
         elements = [53, 1]
         r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
         r(:, 2) = [0.0_dp, 0.0_dp, 25.0_dp]
      case ("xe")
         nat = 1
         allocate (elements(nat), r(3, nat))
         elements = [54]
         r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      case ("csh")
         nat = 2
         allocate (elements(nat), r(3, nat))
         elements = [55, 1]
         r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
         r(:, 2) = [0.0_dp, 0.0_dp, 4.71_dp]
      case ("snh4")
         nat = 5
         allocate (elements(nat), r(3, nat))
         elements = [50, 1, 1, 1, 1]
         r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
         r(:, 2) = [1.8652_dp, 1.8652_dp, 1.8652_dp]
         r(:, 3) = [1.8652_dp, -1.8652_dp, -1.8652_dp]
         r(:, 4) = [-1.8652_dp, 1.8652_dp, -1.8652_dp]
         r(:, 5) = [-1.8652_dp, -1.8652_dp, 1.8652_dp]
      case ("i2")
         nat = 2
         allocate (elements(nat), r(3, nat))
         elements = [53, 53]
         r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
         r(:, 2) = [0.0_dp, 0.0_dp, 5.04_dp]
      case ("ch3i")
         ! C-I 4.04 bohr, C-H 2.06 bohr, H-C-I 107.5 degrees.
         nat = 5
         allocate (elements(nat), r(3, nat))
         elements = [6, 53, 1, 1, 1]
         r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
         r(:, 2) = [0.0_dp, 0.0_dp, 4.04_dp]
         r(:, 3) = [1.9647_dp, 0.0_dp, -0.6195_dp]
         r(:, 4) = [-0.98235_dp, 1.70148_dp, -0.6195_dp]
         r(:, 5) = [-0.98235_dp, -1.70148_dp, -0.6195_dp]
      case default
         nat = 0
      end select
   end subroutine case_geometry

   subroutine test_case(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: case_name, dump
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_ecp_t) :: ecp
      type(error_t) :: err
      integer :: nat, nao, unit
      integer, allocatable :: elements(:)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :), umat(:, :)
      integer, allocatable :: sh_l(:), sh_ao(:)

      call check(error, command_argument_count() >= 2, &
                 "usage: check_ecp_matrix <case_name> <dump.bin>")
      if (allocated(error)) return
      call get_command_argument(1, case_name)
      call get_command_argument(2, dump)

      call case_geometry(trim(case_name), nat, elements, r)
      call check(error, nat > 0, "unknown case_name: "//trim(case_name))
      if (allocated(error)) return

      ! allow_ecp=.true.: this test reads an element's valence-only shells
      ! itself knowing full well it also carries ecp_potentials (it builds
      ! trc_ecp_t from the SAME file right below) -- see
      ! trc_basis_overrides_t%allow_ecp's own comment. Every other caller
      ! of trc_basis_from_json still gets fix/refuse-ecp-basis's refusal.
      call trc_basis_from_json(DEF2_SVP, nat, elements, r, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(allow_ecp=.true.))
      call check(error,.not. err%has_error(), trim(case_name)//" basis must build: "//err%get_message())
      if (allocated(error)) return

      call trc_ecp_from_json(DEF2_SVP, nat, elements, ecp, err)
      call check(error,.not. err%has_error(), trim(case_name)//" ECP must read: "//err%get_message())
      if (allocated(error)) then
         call bas%release()
         return
      end if

      nao = bas%nao
      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      ! trc_1e only visits SURVIVING pairs from the screened pairlist; a
      ! pair this build's threshold drops is simply never written, so
      ! these start at zero rather than whatever the allocator handed back
      ! (seen directly: an uninitialised diagonal entry read back as
      ! huge(1.0_dp)).
      smat = 0.0_dp; tmat = 0.0_dp; vmat = 0.0_dp
      ! S only feeds generate.py's normalisation check, and trc_1e has no
      ! kernels above ONE_E_LMAX: Xe's f shell in an LMAX=2 build wrote past
      ! its block (seen under -fcheck=bounds). Such a build dumps a zero S;
      ! generate the reference from an LMAX=3 build (TERCO_PYSCF_REF_BUILD).
      if (maxval(bas%sh_l) <= ONE_E_LMAX) then
         call pl%build(bas, THRESH)
         call trc_1e(bas, pl, smat, tmat, vmat)
         call pl%release()
      end if

      call trc_ecp_matrix(bas, ecp, umat, err)
      call check(error,.not. err%has_error(), trim(case_name)//" ECP matrix must build: "//err%get_message())
      if (allocated(error)) then
         call bas%release()
         call ecp%release()
         return
      end if
      sh_l = bas%sh_l
      sh_ao = bas%sh_ao
      call bas%release()
      call ecp%release()

      call check_spherical(error, trim(case_name), nat, elements, r, umat)
      if (allocated(error)) return

      open (newunit=unit, file=trim(dump), access='stream', form='unformatted', status='replace')
      write (unit) int(nat, int64), int(nao, int64)
      write (unit) int(elements, int64)
      write (unit) r
      write (unit) smat
      write (unit) umat
      close (unit)

      block
         integer :: u, nat_ref, nao_ref, nfail
         real(dp), allocatable :: u_ref(:, :)
         real(dp) :: worst

         call ref_require('ecp_ints_'//trim(case_name)//'.ref', u)
         call ref_read_int(u, 'NAT', nat_ref)
         call ref_read_int(u, 'NAO', nao_ref)
         call check(error, nat_ref == nat .and. nao_ref == nao, &
                    trim(case_name)//": pyscf reference shape does not match this build")
         if (allocated(error)) then
            call ref_close(u)
            return
         end if
         allocate (u_ref(nao, nao))
         call ref_read_matrix(u, 'U', nao, nao, u_ref)
         call ref_close(u)

         ! The reference is lexicographic; this build's slot (i, j) holds
         ! the reference's (b2l(i), b2l(j)) -- the identity at the default order.
         worst = 0.0_dp
         block
            integer :: i, j, li_, ish, ic
            integer, allocatable :: b2l(:)
            allocate (b2l(nao))
            do ish = 1, size(sh_l)
               li_ = sh_l(ish)
               do ic = 0, (li_ + 1)*(li_ + 2)/2 - 1
                  b2l(sh_ao(ish) + ic) = sh_ao(ish) + cart_perm_b2l(li_, ic)
               end do
            end do
            do j = 1, nao
               do i = 1, nao
                  worst = max(worst, abs(umat(i, j) - u_ref(b2l(i), b2l(j))))
               end do
            end do
         end block
         nfail = 0
         call ref_report(nfail, worst < TOL, 'ECP matrix vs pyscf, max|diff|', worst)
         call check(error, worst < TOL, trim(case_name)// &
                    ": ECP matrix vs pyscf, max|diff| = "//to_string(worst, "es12.4"))
         deallocate (u_ref)
      end block

      deallocate (elements, r, smat, tmat, vmat, umat)
   end subroutine test_case

   subroutine check_spherical(error, case_name, nat, elements, r, ucart)
      !!
      !! The spherical output of trc_ecp_matrix against T^T U_cart T, folded
      !! here on the host, shell block by shell block, from the Cartesian
      !! matrix the caller already checked against pyscf.
      !!
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: case_name
      integer, intent(in) :: nat
      integer, intent(in) :: elements(:)
      real(dp), intent(in) :: r(:, :), ucart(:, :)
      type(trc_basis_t) :: bas
      type(trc_ecp_t) :: ecp
      type(error_t) :: err
      real(dp), allocatable :: usph(:, :), ufold(:, :)
      integer :: ish, jsh, li, lj, ai, aj, bi, bj, p, q, rr, ss, nfail
      real(dp) :: acc, worst

      call trc_basis_from_json(DEF2_SVP, nat, elements, r, bas, err, &
                               overrides=trc_basis_json_overrides(allow_ecp=.true.))
      call check(error,.not. err%has_error(), case_name//" spherical basis must build: "//err%get_message())
      if (allocated(error)) return
      call check(error, bas%spherical, case_name//": def2-SVP must read as spherical")
      if (allocated(error)) then
         call bas%release()
         return
      end if
      call trc_ecp_from_json(DEF2_SVP, nat, elements, ecp, err)
      if (.not. err%has_error()) call trc_ecp_matrix(bas, ecp, usph, err)
      call check(error,.not. err%has_error(), case_name//" spherical ECP matrix must build: "//err%get_message())
      if (allocated(error)) then
         call bas%release()
         call ecp%release()
         return
      end if

      allocate (ufold(bas%nao_sph, bas%nao_sph))
      ufold = 0.0_dp
      do jsh = 1, bas%nshell
         lj = bas%sh_l(jsh); aj = bas%sh_ao(jsh); bj = bas%sh_ao_sph(jsh)
         do ish = 1, bas%nshell
            li = bas%sh_l(ish); ai = bas%sh_ao(ish); bi = bas%sh_ao_sph(ish)
            do ss = 1, 2*lj + 1
               do rr = 1, 2*li + 1
                  acc = 0.0_dp
                  do q = 1, (lj + 1)*(lj + 2)/2
                     do p = 1, (li + 1)*(li + 2)/2
                        acc = acc + SPH_T(p, rr, li)*ucart(ai + p - 1, aj + q - 1)*SPH_T(q, ss, lj)
                     end do
                  end do
                  ufold(bi + rr - 1, bj + ss - 1) = acc
               end do
            end do
         end do
      end do
      worst = maxval(abs(usph - ufold))
      nfail = 0
      call ref_report(nfail, worst < 1.0e-12_dp, 'ECP spherical vs T^T U_cart T, max|diff|', worst)
      call check(error, worst < 1.0e-12_dp, case_name// &
                 ": spherical ECP matrix vs T^T U_cart T, max|diff| = "//to_string(worst, "es12.4"))
      call bas%release()
      call ecp%release()
   end subroutine check_spherical

end module check_ecp_matrix_suite

program check_ecp_matrix
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_ecp_matrix_suite, only: collect_ecp_matrix_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("check_ecp_matrix", collect_ecp_matrix_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_ecp_matrix
