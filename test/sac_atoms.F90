!! terco's SAC guess: the electron count it actually carries, that it is SAD
!! when the molecule is neutral, and that it beats SAD on an ion.
!!
!! The electron count is the point. SAD superposes NEUTRAL atoms, so on a
!! dication it hands the SCF a density holding two electrons too many and the
!! first iterations are spent giving them back. Tr(D S) is what says whether a
!! guess got that right, and it is checked exactly rather than loosely: the
!! atomic blocks are idempotent-free but each atomic SCF converged its own
!! count, so the sum is the molecular count to machine precision or the
!! construction is wrong somewhere.
program sac_atoms
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build, trc_sac_build
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   integer :: natm, nz, dev, ia
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :), smat(:, :), tmat(:, :), vmat(:, :)
   real(dp), allocatable :: d_sad(:, :), d_sac(:, :), dg3(:, :, :), qa(:)
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: r_sad, r_sac
   type(error_t) :: err
   real(dp) :: t_sad, t_sac
   logical :: ok

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   nz = sum(zint)
   opts%conv_energy = 1.0e-10_dp
   ok = .true.

   call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, bas, err)
   if (err%has_error()) then
      print '(a)', "sac_atoms: "//err%get_message(); stop 1
   end if
   call bas%to_device()
   call pl%build(bas, 1.0e-12_dp); call pl%to_device()
   allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
   call trc_1e(bas, pl, smat, tmat, vmat)

   ! --- 1. neutral: SAC carries sum(Z), and reaches SAD's answer ------------
   call trc_sad_build(bas, d_sad, err, verbose=.true.)
   if (.not. err%has_error()) call trc_sac_build(bas, nz, d_sac, err, verbose=.true.)
   if (err%has_error()) then
      print '(a)', "sac_atoms: "//err%get_message(); stop 1
   end if
   t_sad = sum(d_sad*smat); t_sac = sum(d_sac*smat)
   print '(a,f14.10,a,f14.10,a,i0)', "sac_atoms: neutral  Tr(DS) sad ", t_sad, "  sac ", t_sac, &
      "   want ", nz
   if (abs(t_sac - real(nz, dp)) > 1.0e-8_dp) ok = .false.
   if (abs(t_sad - real(nz, dp)) > 1.0e-8_dp) ok = .false.

   allocate (dg3(bas%nao, bas%nao, 1))
   dg3(:, :, 1) = d_sad; call trc_scf_run(bas, nz/2, nz/2, opts, r_sad, dguess=dg3)
   dg3(:, :, 1) = d_sac; call trc_scf_run(bas, nz/2, nz/2, opts, r_sac, dguess=dg3)
   print '(a,f18.12,a,i0,a,f18.12,a,i0,a)', "sac_atoms: neutral  sad E = ", r_sad%energy, &
      " (", r_sad%iterations, ")   sac E = ", r_sac%energy, " (", r_sac%iterations, ")"
   ! SAC on a neutral molecule is SAD, and `build_atomic` chooses the Hund
   ! atom per atom precisely so that it is EXACTLY SAD rather than something
   ! close to it. So the densities themselves are compared, not the energies:
   ! an energy agreeing to 1e-8 would have hidden the earlier version, which
   ! spin-restricted every atom and cost an iteration.
   if (maxval(abs(d_sac - d_sad)) > 1.0e-14_dp) then
      print '(a,es12.4)', "sac_atoms: neutral SAC is not SAD, max|dD| = ", &
         maxval(abs(d_sac - d_sad))
      ok = .false.
   end if
   if (.not. r_sac%converged) ok = .false.
   if (abs(r_sac%energy - r_sad%energy) > 1.0e-8_dp) ok = .false.
   if (r_sac%iterations /= r_sad%iterations) ok = .false.
   deallocate (d_sad, d_sac)

   ! --- 2. the dication, where the two guesses differ ----------------------
   ! Water at +2 is still closed shell, so the comparison needs no open-shell
   ! machinery to confuse it.
   call trc_sad_build(bas, d_sad, err)
   if (.not. err%has_error()) call trc_sac_build(bas, nz - 2, d_sac, err, verbose=.true.)
   if (err%has_error()) then
      print '(a)', "sac_atoms: "//err%get_message(); stop 1
   end if
   t_sad = sum(d_sad*smat); t_sac = sum(d_sac*smat)
   print '(a,f14.10,a,f14.10,a,i0)', "sac_atoms: +2       Tr(DS) sad ", t_sad, "  sac ", t_sac, &
      "   want ", nz - 2
   if (abs(t_sac - real(nz - 2, dp)) > 1.0e-8_dp) ok = .false.
   ! SAD is expected to be WRONG here, by exactly the two electrons it never
   ! removed. If it ever stops being wrong, this test is measuring nothing.
   if (abs(t_sad - real(nz, dp)) > 1.0e-8_dp) ok = .false.

   dg3(:, :, 1) = d_sad
   call trc_scf_run(bas, (nz - 2)/2, (nz - 2)/2, opts, r_sad, dguess=dg3)
   dg3(:, :, 1) = d_sac
   call trc_scf_run(bas, (nz - 2)/2, (nz - 2)/2, opts, r_sac, dguess=dg3)
   print '(a,f18.12,a,i0,a,f18.12,a,i0,a)', "sac_atoms: +2       sad E = ", r_sad%energy, &
      " (", r_sad%iterations, ")   sac E = ", r_sac%energy, " (", r_sac%iterations, ")"
   if (.not. r_sac%converged) ok = .false.
   if (r_sad%converged .and. abs(r_sac%energy - r_sad%energy) > 1.0e-8_dp) ok = .false.
   ! The claim SAC is here to make. Not "fewer iterations on average" -- one
   ! molecule cannot show that -- but that starting from the right electron
   ! count does not cost iterations, which is the weakest thing worth
   ! asserting and the only one a single case supports.
   if (r_sad%converged .and. r_sac%iterations > r_sad%iterations) then
      print '(a)', "sac_atoms: SAC took longer than SAD on the dication"
      ok = .false.
   end if

   ! --- 3. explicit per-atom charges --------------------------------------
   ! The whole +1 taken off the oxygen. Tr(DS) still has to land on the
   ! molecular count, which is what says the override is wired to the atoms
   ! and not merely accepted.
   allocate (qa(natm))
   qa = 0.0_dp
   do ia = 1, natm
      if (zint(ia) == 8) qa(ia) = 1.0_dp
   end do
   deallocate (d_sac)
   call trc_sac_build(bas, nz - 1, d_sac, err, qatom=qa)
   if (err%has_error()) then
      print '(a)', "sac_atoms: "//err%get_message(); stop 1
   end if
   t_sac = sum(d_sac*smat)
   print '(a,f14.10,a,i0)', "sac_atoms: qatom    Tr(DS) sac ", t_sac, "   want ", nz - 1
   if (abs(t_sac - real(nz - 1, dp)) > 1.0e-8_dp) ok = .false.

   ! --- 4. a bad override is refused, not absorbed -------------------------
   block
      type(error_t) :: e2
      real(dp), allocatable :: dbad(:, :)
      real(dp) :: qshort(1)
      qshort = 0.0_dp
      call trc_sac_build(bas, nz, dbad, e2, qatom=qshort)
      if (.not. e2%has_error()) then
         print '(a)', "sac_atoms: a qatom of the wrong length was accepted"
         ok = .false.
      else
         print '(a)', "sac_atoms: qatom length checked -- "//trim(e2%get_message())
      end if
   end block

   call pl%release(); call bas%release()
   if (.not. ok) then
      print '(a)', "sac_atoms: FAIL"
      stop 1
   end if
   print '(a)', "sac_atoms: PASS"
end program sac_atoms
