!
! The FULL KS-DFT nuclear gradient (trc_scf_grad%trc_scf_gradient with
! `func` present) against pyscf's own grad.RKS/grad.UKS, evaluated at
! TERCO'S OWN converged orbitals rather than a re-converged pyscf SCF
! (test/pyscf_refs/generate.py's gen_xc_grad: the same "compare at terco's
! own converged density" rule rks.ref already uses for VV10) and on
! TERCO'S OWN quadrature (`res%grid_r`/`res%grid_w`, dumped below and
! handed to pyscf's `mf.grids` verbatim -- test/pyscf_refs/generate.py's
! gen_rks does the same), so neither an SCF difference nor a quadrature
! difference can hide inside the comparison: what differs is this
! module's formula (and the grid-response term it omits by default --
! see trc_xc_grad.F90's header and check_xc_grad_fd.F90 for that size).
!
! pyscf's `grid_response` defaults to False -- the basis-function term
! only, no travelling grid points or responding partition weights -- the
! same omission trc_xc_grad.F90 makes. An optional 7th argument flips
! BOTH sides to the grid-response gradient: terco's own
! (trc_scf_gradient's `grid_response=.true.`, this module's default) and
! pyscf's (`g%grid_response = .true.`, which needs pyscf's OWN grid, NOT
! the dumped one -- see the block below this header's own comment on
! why). Absent (or 0), this is the original Pulay-only comparison,
! unchanged, against its own already-committed reference.
!
module check_xc_grad_vs_pyscf_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix
   implicit none
   private
   public :: collect_xc_grad_vs_pyscf, xyzfile_arg, basis_json_arg, functional_arg, &
              dump_arg, case_name_arg, grid_level_arg, gresp_arg

   real(dp), parameter :: THRESH = 1.0e-12_dp

   ! All set by the program from argv before run_testsuite.
   character(len=256) :: xyzfile_arg = ''
   !! The molecule, .xyz.
   character(len=256) :: basis_json_arg = ''
   !! The orbital basis, MolSSI BSE JSON.
   character(len=256) :: functional_arg = ''
   !! The functional, by registry name.
   character(len=256) :: dump_arg = ''
   !! The probe file this run writes for generate.py.
   character(len=256) :: case_name_arg = ''
   !! Which case of the reference this run is.
   integer :: grid_level_arg = 3
   !! DFT grid level.
   logical :: gresp_arg = .false.
   !! Whether the grid-response term is included.

contains

   subroutine collect_xc_grad_vs_pyscf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("KS-DFT gradient vs pyscf", test_vs_pyscf) &
                  ]
   end subroutine collect_xc_grad_vs_pyscf

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nat, nao, nspin, nalpha, nbeta, unit, dev
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), grad(:, :)
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_xc_functional_t) :: func
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      logical :: ok
      character(len=256) :: ref_name

      ok = xc_functional_by_name(trim(functional_arg), func)
      call check(error, ok .and. len_trim(func%missing) == 0, "functional not available")
      if (allocated(error)) return

      dev = trc_bind_device(0)
      call read_xyz(trim(xyzfile_arg), nat, at_z, r)
      ! trc_scf_gradient DOES fold a spherical density back to Cartesian now
      ! (check_dft_grad_vs_pyscf.F90 exercises it, and trc_scf_grad.F90's
      ! own header has the formula); force_cart=.true. here is this
      ! suite's own choice, kept so every case matches pyscf's own
      ! cart=True mol below exactly -- not a missing capability.
      call trc_basis_from_json(trim(basis_json_arg), nat, at_z, r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis build: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao
      call bas%to_device()

      ! The electron count is implied by the molecule, not a command-line
      ! argument: water is the closed-shell 5/5, formamide 12/12, the OH
      ! radical the open-shell 5/4 check_full_grad_fd.F90 also uses.
      ! `case_name_arg` already encodes which molecule this is.
      if (index(trim(case_name_arg), 'oh_') == 1 .or. trim(case_name_arg) == 'oh') then
         nalpha = 5; nbeta = 4
      else if (index(trim(case_name_arg), 'nh2_') == 1 .or. trim(case_name_arg) == 'nh2') then
         nalpha = 5; nbeta = 4
      else if (index(trim(case_name_arg), 'formamide') == 1) then
         nalpha = 12; nbeta = 12
      else
         nalpha = 5; nbeta = 5
      end if

      opts%functional = trim(functional_arg)
      opts%grid_level = grid_level_arg
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iteration = 200
      call trc_scf_run(bas, nalpha, nbeta, opts, res)
      call check(error, res%converged, "reference SCF did not converge")
      if (allocated(error)) then
         call bas%release()
         return
      end if
      nspin = res%nspin

      allocate (grad(3, nat))
      call lp%init(nao)
      call trc_scf_gradient(bas, lp, nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad, &
                            func=func, grid_level=grid_level_arg, grid_response=gresp_arg)
      call lp%release()

      ! The overlap, for the AO-permutation/normalisation check generate.py's
      ! gen_1e_grad/gen_rys_grad already run (terco's Cartesian d/f are not
      ! all unit-normalised the way pyscf's cart=True ones are).
      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      call pl%build(bas, THRESH)
      call pl%to_device()
      !$acc enter data create(smat, tmat, vmat)
      call trc_1e(bas, pl, smat, tmat, vmat)
      !$acc update self(smat)
      !$acc exit data delete(smat, tmat, vmat)
      call pl%release()
      call bas%release()

      open (newunit=unit, file=trim(dump_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(nat, int64), int(nao, int64), int(nspin, int64), int(nalpha, int64), int(nbeta, int64)
      write (unit) int(at_z, int64)
      write (unit) r
      write (unit) int(size(res%grid_r, 2), int64)
      write (unit) res%grid_r
      write (unit) res%grid_w
      write (unit) res%cmo
      write (unit) res%eps
      write (unit) grad
      write (unit) smat
      close (unit)
      print '(a)', 'check_xc_grad_vs_pyscf: wrote '//trim(dump_arg)//' for the pyscf-reference comparison'

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      block
         integer :: u, nat_ref, nao_ref
         real(dp), allocatable :: grad_ref(:, :)
         real(dp) :: worst
         if (gresp_arg) then
            ref_name = 'xc_grad_gridresp_'//trim(case_name_arg)//'.ref'
         else
            ref_name = 'xc_grad_'//trim(case_name_arg)//'.ref'
         end if
         call ref_require(trim(ref_name), u)
         call ref_read_int(u, 'NAT', nat_ref)
         call ref_read_int(u, 'NAO', nao_ref)
         call check(error, nat_ref == nat .and. nao_ref == nao, &
                    "pyscf reference shape does not match this build")
         if (allocated(error)) then
            call ref_close(u)
            deallocate (at_z, r, grad, smat, tmat, vmat)
            return
         end if
         allocate (grad_ref(3, nat))
         call ref_read_matrix(u, 'GRAD', 3, nat, grad_ref)
         call ref_close(u)
         worst = maxval(abs(grad - grad_ref))
         if (gresp_arg) then
            ! The grid_response=True cases are two INDEPENDENT
            ! implementations of the same point-motion/partition-weight
            ! physics (trc_xc_grad.F90's header, GRID RESPONSE), not a
            ! Pulay-term-only comparison muddied by a quadrature truncation
            ! -- measured worst|terco-pyscf| at generation time was 1.501e-10
            ! (water/svwn), 1.615e-10 (water/pbe), 1.563e-10 (nh2/pbe, UKS,
            ! the cross-spin Q2/rmatx terms) -- with dft_partition_becke_derivatives'
            ! active-set screening ON (its own default, matching the energy's);
            ! an EARLIER, unscreened measurement here was 5.151e-13/6.470e-13/
            ! 2.481e-11 on the same three cases, so screening costs about two
            ! orders of magnitude of the agreement without threatening it. 1e-9
            ! is ~6x the worst of the three current numbers -- not the
            ! Pulay-only tolerance below loosened: a 30% error in either new
            ! term (B/C), or a screening inconsistency between this gradient and
            ! the energy, would be orders of magnitude above this and caught.
            call check(error, worst < 1.0e-9_dp, &
                       "KS-DFT gradient (grid response) vs pyscf, worst|diff| = "//to_string(worst, "es10.2"))
         else
            call check(error, worst < 3.0e-6_dp, &
                       "KS-DFT gradient vs pyscf, worst|diff| = "//to_string(worst, "es10.2"))
         end if
      end block

      deallocate (at_z, r, grad, smat, tmat, vmat)
   end subroutine test_vs_pyscf

end module check_xc_grad_vs_pyscf_suite

program check_xc_grad_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_xc_grad_vs_pyscf_suite, only: collect_xc_grad_vs_pyscf, xyzfile_arg, basis_json_arg, &
                                            functional_arg, dump_arg, case_name_arg, grid_level_arg, gresp_arg
   implicit none
   integer :: stat, gresp_i
   character(len=256) :: glevel_str, gresp_str

   if (command_argument_count() < 6) then
      write (error_unit, '(a)') 'usage: check_xc_grad_vs_pyscf <xyz> <basis.json> <functional> '// &
         '<grid_level> <dump.bin> <case_name> [grid_response 0|1]'
      error stop 1
   end if
   call get_command_argument(1, xyzfile_arg)
   call get_command_argument(2, basis_json_arg)
   call get_command_argument(3, functional_arg)
   call get_command_argument(4, glevel_str)
   call get_command_argument(5, dump_arg)
   call get_command_argument(6, case_name_arg)
   read (glevel_str, *) grid_level_arg
   gresp_i = 0
   if (command_argument_count() >= 7) then
      call get_command_argument(7, gresp_str)
      read (gresp_str, *) gresp_i
   end if
   gresp_arg = gresp_i /= 0

   stat = 0
   call run_testsuite(collect_xc_grad_vs_pyscf, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_xc_grad_vs_pyscf
