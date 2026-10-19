!
! The FULL density-fitted (`--df`) RHF/RKS nuclear gradient
! (trc_scf_grad%trc_scf_gradient with `dfjk`/`aux`: fitted J and K, 1e, XC,
! nuclear repulsion) against pyscf's own density-fitted gradient
! (pyscf.df.grad.rhf / pyscf.df.grad.rks, auxiliary basis
! def2-universal-jkfit, mol.cart=True) AT TERCO'S OWN converged orbitals and,
! for a functional, TERCO'S OWN grid -- gen_dft_grad's recipe
! (test/pyscf_refs/generate.py's gen_df_grad), so the comparison is of the
! gradient formula at one fixed state, not of two SCFs.
!
! pyscf's DF gradient includes the auxiliary-basis response (both the
! (mu nu|P) aux-centre derivative and the metric's) by default; terco's
! always does. Grid response is OFF on both sides for this comparison
! (pyscf's default); check_df_grad_fd.F90 checks the grid-response-ON
! gradient against a finite difference instead.
!
! Writes the dump gen_df_grad reads: geometry, grid, orbitals, terco's
! gradient and the Cartesian overlap (for the per-AO normalisation).
!
program check_df_grad_vs_pyscf
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_scf_df, only: df_jk_t
   use trc_df_kernels, only: DF_LMAX_AUX
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_report
   implicit none

   real(dp), parameter :: THRESH = 1.0e-12_dp
   real(dp) :: tol
   !! Per case, from test/CMakeLists.txt (the 8th argument): 1e-10 on
   !! water, 5e-9 on formamide. The formamide floor is the FIT's
   !! conditioning, not the gradient: def2-universal-jkfit's Cartesian
   !! metric on formamide has cond(M) = 1.6e14, and a 1e-14 relative
   !! perturbation of pyscf's OWN (mu nu|P) and (P|Q) moves pyscf's own
   !! fitted gradient by 4e-10 (measured), which is the level terco's
   !! integrals and libcint's differ at. terco's derivative integrals agree
   !! with pyscf's to 2e-13 on the same molecule (check_df_grad_ints).
   character(len=256) :: xyzfile, basis_json, aux_json, functional, dump, case_name, glevel_str, tol_str
   integer :: nat, nao, nalpha, nbeta, grid_level, unit, dev, npts
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: r(:, :), grad(:, :)
   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: res
   type(trc_xc_functional_t) :: func
   type(trc_linalg_t) :: lp
   type(df_jk_t) :: dfjk
   type(error_t) :: err
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
   logical :: is_dft

   if (command_argument_count() < 8) then
      print '(a)', 'usage: check_df_grad_vs_pyscf <xyz> <basis.json> <aux.json> <functional|hf> '// &
         '<grid_level> <dump.bin> <case_name> <tol>'; stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basis_json)
   call get_command_argument(3, aux_json)
   call get_command_argument(4, functional)
   call get_command_argument(5, glevel_str)
   call get_command_argument(6, dump)
   call get_command_argument(7, case_name)
   call get_command_argument(8, tol_str)
   read (glevel_str, *) grid_level
   read (tol_str, *) tol
   is_dft = trim(functional) /= 'hf'
   if (is_dft) then
      if (.not. xc_functional_by_name(trim(functional), func)) then
         print '(a)', 'check_df_grad_vs_pyscf: functional not available'; stop 1
      end if
   end if

   dev = trc_bind_device(0)
   call read_xyz(trim(xyzfile), nat, at_z, r)
   nalpha = sum(at_z)/2
   nbeta = nalpha
   call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', 'check_df_grad_vs_pyscf: '//err%get_message(); stop 1
   end if
   call trc_basis_from_json(trim(aux_json), nat, at_z, r, aux, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', 'check_df_grad_vs_pyscf: '//err%get_message(); stop 1
   end if
   if (maxval(aux%sh_l) > DF_LMAX_AUX) then
      print '(a)', '  SKIP: the auxiliary basis needs more angular momentum than this build''s DF_LMAX_AUX.'
      stop 0
   end if
   nao = bas%nao
   call bas%to_device()
   call aux%to_device()
   call pl%build(bas, THRESH)
   call pl%to_device()
   call lp%init(nao)
   call dfjk%setup(bas, aux, pl, lp, err, nocc_hint=nalpha)
   if (err%has_error()) then
      print '(a)', 'check_df_grad_vs_pyscf: '//err%get_message(); stop 1
   end if

   opts%df = .true.
   opts%guess = "gwh"
   if (is_dft) opts%functional = trim(functional)
   opts%grid_level = grid_level
   opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-10_dp; opts%max_iter = 200
   call trc_scf_run(bas, nalpha, nbeta, opts, res, la=lp, dfjk=dfjk)
   if (.not. res%converged) then
      print '(a)', 'check_df_grad_vs_pyscf: reference SCF did not converge'; stop 1
   end if

   allocate (grad(3, nat))
   if (is_dft) then
      call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad, &
                            func=func, grid_level=grid_level, grid_max_pts=opts%grid_max_pts, &
                            rho_tol=opts%rho_tol, grid_response=.false., dfjk=dfjk, aux=aux)
   else
      call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad, &
                            dfjk=dfjk, aux=aux)
   end if

   allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
   !$acc enter data create(smat, tmat, vmat)
   call trc_1e(bas, pl, smat, tmat, vmat)
   !$acc update self(smat)
   !$acc exit data delete(smat, tmat, vmat)
   call dfjk%release()
   call pl%release()
   call lp%release()
   call aux%release()
   call bas%release()

   npts = 0
   if (is_dft) npts = size(res%grid_r, 2)
   open (newunit=unit, file=trim(dump), access='stream', form='unformatted', status='replace')
   write (unit) int(nat, int64), int(nao, int64), int(nalpha, int64)
   write (unit) int(at_z, int64)
   write (unit) r
   write (unit) int(npts, int64)
   if (is_dft) then
      write (unit) res%grid_r
      write (unit) res%grid_w
   end if
   write (unit) res%cmo
   write (unit) res%eps
   write (unit) grad
   write (unit) smat
   close (unit)
   print '(a)', 'check_df_grad_vs_pyscf: wrote '//trim(dump)

   block
      integer :: u, nat_ref, nao_ref, nfail
      real(dp), allocatable :: grad_ref(:, :)
      real(dp) :: worst
      call ref_require('dfg_'//trim(case_name)//'.ref', u)
      call ref_read_int(u, 'NAT', nat_ref)
      call ref_read_int(u, 'NAO', nao_ref)
      if (nat_ref /= nat .or. nao_ref /= nao) then
         print '(a)', 'check_df_grad_vs_pyscf: pyscf reference shape does not match this build'; stop 1
      end if
      allocate (grad_ref(3, nat))
      call ref_read_matrix(u, 'GRAD', 3, nat, grad_ref)
      call ref_close(u)
      nfail = 0
      worst = maxval(abs(grad - grad_ref))
      call ref_report(nfail, worst < tol, '--df gradient vs pyscf df.grad, worst|diff|', worst)
      print '(a,es10.2)', '  translational-invariance residual: ', maxval(abs(sum(grad, dim=2)))
      if (nfail > 0) then
         print '(a)', 'check_df_grad_vs_pyscf: FAILED against the pyscf reference'; stop 1
      end if
      print '(a)', 'check_df_grad_vs_pyscf: RESULT: PASS'
   end block

   deallocate (at_z, r, grad, smat, tmat, vmat)

end program check_df_grad_vs_pyscf
