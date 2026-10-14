!
! trc_df_2c and trc_df_3c against pyscf's int2c2e/int3c2e, same geometry,
! orbital and auxiliary bases read from the same JSON files pyscf reads --
! the density-fitting analogue of check_rys_fock_vs_pyscf.F90, and the one
! that can exercise a g SHELL (l=4) on the auxiliary side, not just an f
! orbital shell.
!
! SKIPS CLEANLY (stop 0, after printing SKIP) when either basis needs more
! angular momentum than this build's DF_LMAX/DF_LMAX_AUX reach, same
! convention as check_rys_fock_vs_pyscf.F90's TRC_LMAX skip -- a build at
! the shipped ceiling (DF_LMAX=2, DF_LMAX_AUX=3) has nothing with a g
! shell to run this against, and that is not a failure.
!
! PER-AO NORMALISATION WITHOUT trc_1e. check_rys_fock_vs_pyscf.F90 derives
! its normalisation factors from the overlap matrix (trc_1e), but trc_1e's
! own ceiling is ONE_E_LMAX -- independent of DF_LMAX, and not raised by
! TERCO_DF_LMAX alone (this build can have DF_LMAX=3 with ONE_E_LMAX still
! at 2, and one_e_dispatch has no case above its own ceiling: calling it
! here would read an uninitialised block, not a wrong number caught at
! source). Used instead: trc_df_2c on the ORBITAL basis treated as its own
! auxiliary set, i.e. the (i|j) two-centre Coulomb metric over the SAME
! shells trc_df_3c's bra uses. Its diagonal (i|i) is, like any two-centre
! integral's diagonal, strictly positive and scales as n_i^2 under a
! per-AO rescaling exactly the way the overlap diagonal does -- n is a
! property of the represented basis function, not of which operator the
! integral happens to be -- so it is as valid a normalisation probe as S,
! and it only needs DF_LMAX, which this program has already checked.
!
! Writes the dump df_vs_pyscf_ref.py reads: geometry, the orbital (i|j)
! metric, trc_df_2c's (P|Q), and trc_df_3c's (mu nu|P).
!
program check_df_vs_pyscf
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_2c, trc_df_3c
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_basis_json, only: trc_basis_from_json
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none

   real(dp), parameter :: THRESH = 1.0e-14_dp
   character(len=256) :: xyzfile, basis_json, aux_json, dump
   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl
   type(error_t) :: err
   integer :: nat, unit, nao, naux
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: r(:, :)
   real(dp), allocatable :: omat(:, :), jmat(:, :), tens(:, :, :)

   if (command_argument_count() < 4) then
      print '(a)', 'usage: check_df_vs_pyscf <xyz> <basis.json> <aux.json> <dump.bin>'
      stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basis_json)
   call get_command_argument(3, aux_json)
   call get_command_argument(4, dump)

   call read_xyz(trim(xyzfile), nat, at_z, r)
   call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', 'check_df_vs_pyscf: '//err%get_message(); stop 1
   end if
   call trc_basis_from_json(trim(aux_json), nat, at_z, r, aux, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', 'check_df_vs_pyscf: '//err%get_message(); stop 1
   end if

   ! See this program's own header: a basis above this build's ceiling is
   ! nothing to test yet, not an error.
   if (maxval(bas%sh_l) > DF_LMAX) then
      print '(a,i0,a,i0,a)', '  SKIP: orbital basis needs l = ', maxval(bas%sh_l), &
         ' and this build has DF_LMAX = ', DF_LMAX, '.'
      stop 0
   end if
   if (maxval(aux%sh_l) > DF_LMAX_AUX) then
      print '(a,i0,a,i0,a)', '  SKIP: auxiliary basis needs l = ', maxval(aux%sh_l), &
         ' and this build has DF_LMAX_AUX = ', DF_LMAX_AUX, '.'
      stop 0
   end if

   call bas%to_device()
   call aux%to_device()
   call pl%build(bas, THRESH)
   call pl%to_device()

   nao = bas%nao; naux = aux%nao
   allocate (omat(nao, nao), jmat(naux, naux), tens(nao, nao, naux))
   !$acc enter data create(omat, jmat, tens)
   call trc_df_2c(bas, omat)
   call trc_df_2c(aux, jmat)
   call trc_df_3c(bas, pl, aux, tens)
   !$acc update self(omat, jmat, tens)
   !$acc exit data delete(omat, jmat, tens)

   call pl%release()
   call bas%release()
   call aux%release()

   open (newunit=unit, file=trim(dump), access='stream', &
         form='unformatted', status='replace')
   write (unit) int(nat, kind=8), int(nao, kind=8), int(naux, kind=8)
   write (unit) int(at_z, kind=8)
   write (unit) r
   write (unit) omat
   write (unit) jmat
   write (unit) tens
   close (unit)
   print '(a,i0,a,i0,a,i0,a)', 'check_df_vs_pyscf: wrote '//trim(dump)//' (', &
      nat, ' atoms, ', nao, ' orbital AOs, ', naux, ' auxiliary AOs)'

   deallocate (at_z, r, omat, jmat, tens)

end program check_df_vs_pyscf
