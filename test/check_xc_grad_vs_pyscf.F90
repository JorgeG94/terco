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
program check_xc_grad_vs_pyscf
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_report
   implicit none

   real(dp), parameter :: THRESH = 1.0e-12_dp
   character(len=256) :: xyzfile, basis_json, functional, dump, case_name, glevel_str, gresp_str
   integer :: nat, nao, nspin, nalpha, nbeta, grid_level, unit, dev, gresp_i
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
   logical :: ok, gresp
   character(len=256) :: ref_name

   if (command_argument_count() < 6) then
      print '(a)', 'usage: check_xc_grad_vs_pyscf <xyz> <basis.json> <functional> '// &
         '<grid_level> <dump.bin> <case_name> [grid_response 0|1]'; stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basis_json)
   call get_command_argument(3, functional)
   call get_command_argument(4, glevel_str)
   call get_command_argument(5, dump)
   call get_command_argument(6, case_name)
   read (glevel_str, *) grid_level
   gresp_i = 0
   if (command_argument_count() >= 7) then
      call get_command_argument(7, gresp_str)
      read (gresp_str, *) gresp_i
   end if
   gresp = gresp_i /= 0

   ok = xc_functional_by_name(trim(functional), func)
   if (.not. ok .or. len_trim(func%missing) > 0) then
      print '(a)', 'check_xc_grad_vs_pyscf: functional not available'; stop 1
   end if

   dev = trc_bind_device(0)
   call read_xyz(trim(xyzfile), nat, at_z, r)
   ! trc_scf_gradient (this branch) takes cmo/eps straight from
   ! trc_scf_run's result at `bas%nao`: a spherical basis (cc-pVDZ's JSON
   ! tags "gto_spherical") would make that SMALLER than bas%nao, which the
   ! gradient has no fold-back for yet -- the sibling branch's "spherical
   ! density folding" (this program's own header, and trc_scf_grad.F90's)
   ! is what adds it. force_cart=.true. keeps every basis here Cartesian,
   ! matching pyscf's own cart=True mol below exactly.
   call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', 'check_xc_grad_vs_pyscf: '//err%get_message(); stop 1
   end if
   nao = bas%nao
   call bas%to_device()

   ! The electron count is implied by the molecule, not a command-line
   ! argument: water is the closed-shell 5/5, formamide 12/12, the OH
   ! radical the open-shell 5/4 check_full_grad_fd.F90 also uses.
   ! `case_name` already encodes which molecule this is.
   if (index(trim(case_name), 'oh_') == 1 .or. trim(case_name) == 'oh') then
      nalpha = 5; nbeta = 4
   else if (index(trim(case_name), 'nh2_') == 1 .or. trim(case_name) == 'nh2') then
      nalpha = 5; nbeta = 4
   else if (index(trim(case_name), 'formamide') == 1) then
      nalpha = 12; nbeta = 12
   else
      nalpha = 5; nbeta = 5
   end if

   opts%functional = trim(functional)
   opts%grid_level = grid_level
   opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-9_dp; opts%max_iter = 200
   call trc_scf_run(bas, nalpha, nbeta, opts, res)
   if (.not. res%converged) then
      print '(a)', 'check_xc_grad_vs_pyscf: reference SCF did not converge'; stop 1
   end if
   nspin = res%nspin

   allocate (grad(3, nat))
   call lp%init(nao)
   call trc_scf_gradient(bas, lp, nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad, &
                         func=func, grid_level=grid_level, grid_response=gresp)
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

   open (newunit=unit, file=trim(dump), access='stream', form='unformatted', status='replace')
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
   print '(a)', 'check_xc_grad_vs_pyscf: wrote '//trim(dump)//' for the pyscf-reference comparison'

   ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
   block
      integer :: u, nat_ref, nao_ref, nfail
      real(dp), allocatable :: grad_ref(:, :)
      real(dp) :: worst
      if (gresp) then
         ref_name = 'xc_grad_gridresp_'//trim(case_name)//'.ref'
      else
         ref_name = 'xc_grad_'//trim(case_name)//'.ref'
      end if
      call ref_require(trim(ref_name), u)
      call ref_read_int(u, 'NAT', nat_ref)
      call ref_read_int(u, 'NAO', nao_ref)
      if (nat_ref /= nat .or. nao_ref /= nao) then
         print '(a)', 'check_xc_grad_vs_pyscf: pyscf reference shape does not match this build'; stop 1
      end if
      allocate (grad_ref(3, nat))
      call ref_read_matrix(u, 'GRAD', 3, nat, grad_ref)
      call ref_close(u)
      nfail = 0
      worst = maxval(abs(grad - grad_ref))
      call ref_report(nfail, worst < 3.0e-6_dp, 'KS-DFT gradient vs pyscf, worst|diff|', worst)
      if (nfail > 0) then
         print '(a)', 'check_xc_grad_vs_pyscf: FAILED against the pyscf reference'; stop 1
      end if
      print '(a)', 'check_xc_grad_vs_pyscf: RESULT: PASS'
   end block

   deallocate (at_z, r, grad, smat, tmat, vmat)

end program check_xc_grad_vs_pyscf
