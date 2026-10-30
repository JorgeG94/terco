!
! The SCF driver end to end, against pyscf on the same grid.
!
! trc_scf is the half of terco's contract that used to be the caller's:
! basis and electron count in, converged energy and density out, every
! matrix resident on the device. This runs it on two molecules in 6-31G
! UNCONTRACTED to its primitives -- every primitive its own shell with
! coefficient one, because pyscf renormalises a contracted shell to unit
! self-overlap and terco does not, and a total energy compared at 1e-8
! cannot carry the 1e-6 that leaves in the basis -- and writes what pyscf
! needs to reproduce each case exactly: the shells as built, the grid, the
! energy and the final density.
!
!   water, closed shell:  RHF, then SVWN, PBE, BLYP, B3LYP, wB97X-V, wB97X,
!                         wB97, wB97M-V restricted
!   OH radical, doublet:  UHF, then PBE, wB97X-V, wB97X, wB97M-V unrestricted
!
! wB97M-V is the meta-GGA case: same RSH/VV10 machinery as wB97X-V plus
! tau, so it exercises res%grid_w_nlc (VV10's own, coarser grid) and the
! meta-GGA dispatch together, the one combination nothing else here covers.
!
! So the restricted and unrestricted paths, Hartree-Fock and Kohn-Sham,
! LDA, GGA, global and range-separated hybrid, with and without VV10's
! nonlocal correlation, all go through the one driver. `rks_ref.py` rebuilds
! the same basis in the same order and runs pyscf's RHF/UHF/RKS/UKS on the
! same points from terco's converged density -- VV10's own (coarser) grid
! too, for a functional that has it, set as pyscf's `nlcgrids` exactly as
! the semilocal grid is set as `grids`.
!
! Fourteen independent (molecule, functional) cases below, each its own
! test-drive unittest; the shared setup (per-molecule uncontracted basis,
! the probe write, the golden dump and the sequential pyscf-reference
! read) runs once -- see check_xc_energy_suite's own `computed` comment
! for why, and why this suite stays `parallel=.false.`.
!
module scf_rks_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz, uncontract_basis
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix3
   use golden_io, only: golden_begin, golden_put, golden_end
   implicit none
   private
   public :: collect_scf_rks

   integer, parameter :: NCASE = 14
   character(len=8), parameter :: mol(NCASE) = [character(len=8) :: &
                                                "water", "water", "water", "water", "water", "water", "water", &
                                                "water", "water", "oh", "oh", "oh", "oh", "oh"]
   character(len=8), parameter :: fn(NCASE) = [character(len=8) :: &
                                               "", "svwn", "pbe", "blyp", "b3lyp", "wb97x-v", "wb97x", "wb97", &
                                               "wb97m-v", "", "pbe", "wb97x-v", "wb97x", "wb97m-v"]
   integer, parameter :: nalpha(NCASE) = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
   integer, parameter :: nbeta(NCASE) = [5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4]

   logical :: computed = .false.
   logical :: ok_case(NCASE)
   character(len=256) :: msg_case(NCASE)

contains

   subroutine collect_scf_rks(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest(case_name(1), case_01), &
                  new_unittest(case_name(2), case_02), &
                  new_unittest(case_name(3), case_03), &
                  new_unittest(case_name(4), case_04), &
                  new_unittest(case_name(5), case_05), &
                  new_unittest(case_name(6), case_06), &
                  new_unittest(case_name(7), case_07), &
                  new_unittest(case_name(8), case_08), &
                  new_unittest(case_name(9), case_09), &
                  new_unittest(case_name(10), case_10), &
                  new_unittest(case_name(11), case_11), &
                  new_unittest(case_name(12), case_12), &
                  new_unittest(case_name(13), case_13), &
                  new_unittest(case_name(14), case_14) &
                ]
   end subroutine collect_scf_rks

   function case_name(ic) result(name)
      integer, intent(in) :: ic
      character(len=:), allocatable :: name
      if (len_trim(fn(ic)) > 0) then
         name = trim(mol(ic))//"/"//trim(fn(ic))
      else
         name = trim(mol(ic))//"/hf"
      end if
   end function case_name

   ! One wrapper per case, each passing its own index: a unittest must not
   ! depend on which unittests ran before it (a name-selected run, or a
   ! reordered suite, would otherwise check the wrong case).
   subroutine case_01(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 1)
   end subroutine case_01

   subroutine case_02(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 2)
   end subroutine case_02

   subroutine case_03(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 3)
   end subroutine case_03

   subroutine case_04(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 4)
   end subroutine case_04

   subroutine case_05(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 5)
   end subroutine case_05

   subroutine case_06(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 6)
   end subroutine case_06

   subroutine case_07(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 7)
   end subroutine case_07

   subroutine case_08(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 8)
   end subroutine case_08

   subroutine case_09(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 9)
   end subroutine case_09

   subroutine case_10(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 10)
   end subroutine case_10

   subroutine case_11(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 11)
   end subroutine case_11

   subroutine case_12(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 12)
   end subroutine case_12

   subroutine case_13(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 13)
   end subroutine case_13

   subroutine case_14(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 14)
   end subroutine case_14

   subroutine test_case(error, ic)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ic
      call ensure_computed()
      call check(error, ok_case(ic), trim(msg_case(ic)))
   end subroutine test_case

   subroutine ensure_computed()
      if (computed) return
      call run_all()
      computed = .true.
   end subroutine ensure_computed

   subroutine run_all()
      integer :: ic, unit, natm, nsh, maxnp, ish
      integer, allocatable :: sh_l(:), sh_np(:), zint(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), at_z(:), at_r(:, :)
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts, opts2
      type(trc_scf_result_t) :: res, res2
      character(len=8) :: current
      integer :: refu, ncase_ref
      character(len=16) :: tag
      character(len=32) :: gkey
      logical :: has_ref

      open (newunit=unit, file="rks_probe.bin", access="stream", form="unformatted", status="replace")
      write (unit) int(NCASE, int64)
      current = ""

      ! Opened (not required) before any case runs: a missing reference must
      ! not stop this suite before rks_probe.bin is written, or
      ! test/pyscf_refs/generate.py has nothing to bootstrap from. The FAIL
      ! happens below, after the dump is safely written.
      call ref_open("rks.ref", refu, has_ref)
      call golden_begin("scf_rks")
      if (has_ref) then
         call ref_read_int(refu, "NCASE", ncase_ref)
         if (ncase_ref /= NCASE) error stop "scf_rks: pyscf reference case count does not match this build"
      end if
      ok_case = .false.
      msg_case = "pyscf reference missing"

      do ic = 1, NCASE
         if (mol(ic) /= current) then
            if (len_trim(current) > 0) call bas%release()
            current = mol(ic)
            call read_xyz(trim(current)//".xyz", natm, zint, at_r)
            if (allocated(at_z)) deallocate (at_z)
            allocate (at_z(natm))
            at_z = real(zint, dp)
            block
               type(trc_basis_t) :: jbas
               type(error_t) :: err
               call trc_basis_from_json("basis_sets/6-31g.json", natm, zint, at_r, jbas, err)
               if (err%has_error()) error stop "scf_rks: "//err%get_message()
               call uncontract_basis(jbas, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
            end block
            call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, at_z, at_r, maxnp)
            call bas%to_device()
            print "(a)", ""
            print "(a,a,a,i0,a,i0,a)", "scf_rks: ", trim(current), ", 6-31G uncontracted, ", nsh, " shells, ", &
               bas%nao, " AOs"
         end if
         opts = trc_scf_options_t()
         opts%functional = fn(ic)
         opts%grid_level = 3
         opts%conv_energy = 1.0e-11_dp
         opts%conv_density = 1.0e-8_dp
         ! conv_density no longer gates (trc_scf_driver.F90); the commutator norm does.
         ! At its default 1e-5 the slow UKS wB97X-V case could stop up to ~2e-7
         ! above its minimum (seen under ifx), past the 1e-8 gate below, while
         ! terco's energy AT pyscf's density agreed to 2.5e-11. Converge to the
         ! gate's own precision instead.
         opts%conv_diis = 1.0e-8_dp
         ! OH's unpaired electron sits in one of a rotationally degenerate
         ! pair of pi orbitals; starting terco's OWN SCF from the same
         ! symmetry-broken density test/pyscf_refs/generate.py's gen_rks
         ! commits (DM0_i, read here before the NAME_i/DM_i the at-density
         ! control reads later -- same per-case order the file was written
         ! in) means both terco and pyscf converge toward the SAME
         ! orientation rather than two that merely share an energy. Only for
         ! the open-shell cases (nalpha /= nbeta) that record carries; a
         ! missing reference (the first, bootstrap run) falls back to the
         ! guess trc_scf_run would pick on its own, same as always.
         if (has_ref .and. nalpha(ic) /= nbeta(ic)) then
            block
               character(len=16) :: name_chk
               integer :: nspin_chk, nao_chk
               real(dp), allocatable :: dm0(:, :, :)
               write (tag, "(i0)") ic - 1
               call ref_read_str(refu, "NAME_"//trim(tag), name_chk)
               call ref_read_int(refu, "NSPIN_"//trim(tag), nspin_chk)
               call ref_read_int(refu, "NAO_"//trim(tag), nao_chk)
               if (nspin_chk /= 2 .or. nao_chk /= bas%nao) &
                  error stop "scf_rks: pyscf reference shape does not match this build"
               allocate (dm0(bas%nao, bas%nao, 2))
               call ref_read_matrix3(refu, "DM0_"//trim(tag), bas%nao, bas%nao, 2, dm0)
               call trc_scf_run(bas, nalpha(ic), nbeta(ic), opts, res, dguess=dm0)
            end block
         else
            call trc_scf_run(bas, nalpha(ic), nbeta(ic), opts, res)
         end if
         if (.not. res%converged) print "(a,a)", "  NOT CONVERGED: ", trim(res%message)
         if (len_trim(fn(ic)) > 0) then
            print "(a,a8,a,i1,a,i3,a,f22.12,a,f14.10)", "  ", fn(ic), " nspin ", res%nspin, " iter ", res%iterations, &
               "  E = ", res%energy, "  N(grid) = ", res%nelec_grid
         else
            print "(a,a8,a,i1,a,i3,a,f22.12)", "  ", "hf", " nspin ", res%nspin, " iter ", res%iterations, &
               "  E = ", res%energy
         end if
         gkey = trim(mol(ic))//"_"//merge(fn(ic), "hf      ", len_trim(fn(ic)) > 0)
         call golden_put("E_"//trim(gkey), res%energy)

         ! Everything pyscf needs to reproduce this case.
         write (unit) int(bas%nao, int64), int(natm, int64), int(nsh, int64), int(res%nspin, int64), &
            int(nalpha(ic), int64), int(nbeta(ic), int64)
         write (unit) at_z, at_r
         do ish = 1, nsh
            write (unit) int(centre_of(sh_r, at_r, natm, ish), int64), int(sh_l(ish), int64), sh_e(1, ish)
         end do
         write (unit) fn(ic)
         if (len_trim(fn(ic)) > 0) then
            write (unit) int(size(res%grid_w), int64)
            write (unit) res%grid_r, res%grid_w
         else
            write (unit) 0_int64
         end if
         if (allocated(res%grid_w_nlc)) then
            write (unit) int(size(res%grid_w_nlc), int64)
            write (unit) res%grid_r_nlc, res%grid_w_nlc
         else
            write (unit) 0_int64
         end if
         write (unit) res%energy
         write (unit) res%dmat

         ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
         if (has_ref) then
            block
               character(len=16) :: name_ref
               integer :: nspin_ref, nao_ref
               real(dp) :: e_ref, de, dat, de_tol, pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm
               real(dp), allocatable :: dm_ref(:, :, :)
               write (tag, "(i0)") ic - 1
               ! NAME_/NSPIN_/NAO_ for this case were already read above, before
               ! the live SCF, for every open-shell case (nalpha /= nbeta) --
               ! DM0_'s own position in the file, right after them -- so only a
               ! closed-shell case still has them here next in sequence.
               if (nalpha(ic) == nbeta(ic)) then
                  call ref_read_str(refu, "NAME_"//trim(tag), name_ref)
                  call ref_read_int(refu, "NSPIN_"//trim(tag), nspin_ref)
                  call ref_read_int(refu, "NAO_"//trim(tag), nao_ref)
               else
                  name_ref = fn(ic)
                  if (len_trim(name_ref) == 0) name_ref = "hf"
                  nspin_ref = res%nspin
                  nao_ref = bas%nao
               end if
               call ref_read_real(refu, "E_"//trim(tag), e_ref)
               ! Provenance only (not asserted on here): the pyscf convergence
               ! this energy was committed at -- see generate.py's converge_tight,
               ! which refuses a case whose gradient norm does not clear
               ! conv_tol_grad, tighter than pyscf's own sqrt(conv_tol) default.
               call ref_read_real(refu, "CONV_TOL_"//trim(tag), pyscf_conv_tol)
               call ref_read_real(refu, "CONV_TOL_GRAD_"//trim(tag), pyscf_conv_tol_grad)
               call ref_read_real(refu, "GRAD_NORM_"//trim(tag), pyscf_grad_norm)
               if (nspin_ref /= res%nspin .or. nao_ref /= bas%nao) &
                  error stop "scf_rks: pyscf reference shape does not match this build"
               allocate (dm_ref(bas%nao, bas%nao, res%nspin))
               call ref_read_matrix3(refu, "DM_"//trim(tag), bas%nao, bas%nao, res%nspin, dm_ref)
               de = abs(res%energy - e_ref)
               ! The at-density control: terco's OWN energy at pyscf's committed,
               ! converged density -- max_iteration=1 with dguess runs exactly one
               ! Roothaan iteration, whose etot (src/scf/trc_scf_driver.F90's own "energy at
               ! this density" comment) is computed from dguess BEFORE the
               ! diagonalisation updates it, i.e. exactly energy_tot(dm0)'s
               ! equivalent -- never an iterated answer.
               opts2 = opts
               opts2%max_iteration = 1
               call trc_scf_run(bas, nalpha(ic), nbeta(ic), opts2, res2, dguess=dm_ref)
               dat = abs(res2%energy - e_ref)
               call golden_put("E_AT_PYSCF_D_"//trim(gkey), res2%energy)
               print "(a,a,a,es10.2,a,es10.2)", "    ", trim(name_ref), &
                  ": de(converged) ", de, "  de(at pyscf D) ", dat
               ! The CONVERGED energy of an open-shell DFT case gets a looser bound
               ! than everything else. OH's 2-Pi hole can sit anywhere in the
               ! degenerate pi plane, and with a spherically averaged SAD guess
               ! nothing picks the orientation, so round-off does -- which differs
               ! with the order of an OpenMP reduction. The grid is not invariant
               ! under an arbitrary rotation about the bond, so orientations that
               ! are equivalent exactly differ by grid error, measured up to
               ! 6.1e-7 (pbe) on the OpenMP build. The at-density check is what
               ! pins correctness, and it stays at 1e-8 for every case, as does the
               ! converged energy for HF (no grid) and every closed-shell case.
               if (res%nspin == 2 .and. len_trim(fn(ic)) > 0) then
                  de_tol = 2.0e-6_dp
               else
                  de_tol = 1.0e-8_dp
               end if
               ok_case(ic) = res%converged .and. de < de_tol .and. dat < 1.0e-8_dp
               write (msg_case(ic), "(a,es10.2)") trim(name_ref)// &
                  ": E(converged) and E(terco at pyscf D) vs pyscf, worst = ", max(de, dat)
               if (.not. res%converged) msg_case(ic) = "SCF did not converge ("//trim(res%message)//"); "// &
                                                       trim(msg_case(ic))
               deallocate (dm_ref)
            end block
         end if
      end do
      call golden_end()
      close (unit)
      call bas%release()
      if (has_ref) call ref_close(refu)
      if (.not. has_ref) error stop "FAIL: missing reference pyscf_refs/rks.ref -- run " &
         //"test/pyscf_refs/generate.py (see test/pyscf_refs/README.md)"
   end subroutine run_all

   integer function centre_of(sh_r, at_r, natm, ish)
      real(dp), intent(in) :: sh_r(:, :), at_r(:, :)
      integer, intent(in) :: natm, ish
      integer :: a
      centre_of = 0
      do a = 1, natm
         if (all(abs(sh_r(:, ish) - at_r(:, a)) < 1.0e-12_dp)) centre_of = a
      end do
      if (centre_of == 0) error stop "scf_rks: a shell is on no atom"
   end function centre_of

end module scf_rks_suite

program scf_rks
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_rks_suite, only: collect_scf_rks
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_rks, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_rks
