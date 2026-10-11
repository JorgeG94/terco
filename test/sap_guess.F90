!! terco's SAP guess: the screening potential it builds, and whether starting
!! from it is better than starting from GWH.
!!
!! TWO CHECKS THAT CAN FAIL INDEPENDENTLY OF THE SCF.
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
program sap_guess
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build, trc_atom_ranges
   use trc_sap, only: trc_sap_build, trc_sap_screen_table
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   integer :: natm, nz, dev, ia, a0, a1, nr
   integer, allocatable :: zint(:), atom_of(:), first(:), last(:)
   real(dp), allocatable :: at_r(:, :), smat(:, :), tmat(:, :), vmat(:, :)
   real(dp), allocatable :: dsad(:, :), hsap(:, :), rtab(:), vtab(:), dg3(:, :, :)
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: r_gwh, r_sad, r_sap
   type(error_t) :: err
   real(dp) :: nel, tail, worst_n, worst_t
   logical :: ok

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   nz = sum(zint)
   opts%conv_energy = 1.0e-10_dp
   ok = .true.

   call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, bas, err)
   if (err%has_error()) then
      print '(a)', "sap_guess: "//err%get_message(); stop 1
   end if
   call bas%to_device()
   call pl%build(bas, 1.0e-12_dp); call pl%to_device()
   allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
   call trc_1e(bas, pl, smat, tmat, vmat)

   ! --- 1. the radial screening table, per atom ---------------------------
   call trc_sad_build(bas, dsad, err)
   if (.not. err%has_error()) call trc_atom_ranges(bas, atom_of, first, last, err)
   if (err%has_error()) then
      print '(a)', "sap_guess: "//err%get_message(); stop 1
   end if
   worst_n = 0.0_dp; worst_t = 0.0_dp
   do ia = 1, natm
      a0 = bas%sh_ao(first(ia))
      a1 = bas%sh_ao(last(ia)) + (bas%sh_l(last(ia)) + 1)*(bas%sh_l(last(ia)) + 2)/2 - 1
      call trc_sap_screen_table(bas, first(ia), last(ia), bas%at_r(:, ia), zint(ia), &
                                dsad(a0:a1, a0:a1), rtab, vtab, nel, err)
      if (err%has_error()) then
         print '(a)', "sap_guess: "//err%get_message(); stop 1
      end if
      nr = size(rtab)
      tail = rtab(nr)*vtab(nr)
      print '(a,i0,a,i0,a,f14.9,a,i0,a,f14.9)', "sap_guess: Z = ", zint(ia), &
         "  want ", zint(ia), " electrons, radial density gives ", nel, &
         "   r*Vscr at r = ", nint(rtab(nr)), " is ", tail
      worst_n = max(worst_n, abs(nel - real(zint(ia), dp)))
      worst_t = max(worst_t, abs(tail - real(zint(ia), dp)))
      deallocate (rtab, vtab)
   end do
   ! The radial mesh is a quadrature, so these are quadrature-accurate and not
   ! exact. A part in 1e-6 of an electron is far tighter than a guess needs
   ! and loose enough that the mesh size is not being asserted.
   print '(a,es12.4,a,es12.4)', "sap_guess: worst |N - Z| ", worst_n, "   worst |r*Vscr - Z| ", worst_t
   if (worst_n > 1.0e-6_dp) ok = .false.
   if (worst_t > 1.0e-6_dp) ok = .false.

   ! --- 2. SAP as a guess, against GWH and SAD ----------------------------
   allocate (hsap(bas%nao, bas%nao))
   call trc_sap_build(bas, pl, tmat, vmat, hsap, err, verbose=.true.)
   if (err%has_error()) then
      print '(a)', "sap_guess: "//err%get_message(); stop 1
   end if
   if (maxval(abs(hsap - transpose(hsap))) > 1.0e-10_dp) then
      print '(a)', "sap_guess: H_sap is not symmetric"
      ok = .false.
   end if

   opts%guess = "gwh"
   call trc_scf_run(bas, nz/2, nz/2, opts, r_gwh)
   allocate (dg3(bas%nao, bas%nao, 1)); dg3(:, :, 1) = dsad
   call trc_scf_run(bas, nz/2, nz/2, opts, r_sad, dguess=dg3)
   call trc_scf_run(bas, nz/2, nz/2, opts, r_sap, hguess=hsap)
   print '(a,f18.12,a,i0,a)', "sap_guess: gwh E = ", r_gwh%energy, " (", r_gwh%iterations, ")"
   print '(a,f18.12,a,i0,a)', "sap_guess: sad E = ", r_sad%energy, " (", r_sad%iterations, ")"
   print '(a,f18.12,a,i0,a)', "sap_guess: sap E = ", r_sap%energy, " (", r_sap%iterations, ")"
   if (.not. r_sap%converged) ok = .false.
   ! One SCF solution from three guesses. Only the energy has to agree; the
   ! iteration counts are reported rather than asserted, because water in
   ! 6-31G converges from almost anything and a count that happened to be
   ! equal here would not be evidence about SAP.
   if (abs(r_sap%energy - r_gwh%energy) > 1.0e-8_dp) ok = .false.
   if (abs(r_sap%energy - r_sad%energy) > 1.0e-8_dp) ok = .false.

   call pl%release(); call bas%release()
   if (.not. ok) then
      print '(a)', "sap_guess: FAIL"
      stop 1
   end if
   print '(a)', "sap_guess: PASS"
end program sap_guess
