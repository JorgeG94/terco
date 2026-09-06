!! The context surface against the Fortran drivers: water/cc-pVDZ RHF and
!! RI-MP2 through one handle, set up in steps, must reproduce trc_scf_run
!! and trc_rimp2_run to the last digit, and every misuse must be a status.
program context_scf
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char, c_char
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_method, trc_set_convergence, trc_set_rimp2, trc_run_scf, trc_run_rimp2, &
                               trc_nao, trc_nspin, trc_energy, trc_converged, trc_iterations, trc_density, &
                               trc_mo_energies, trc_rimp2_energy, trc_message, trc_set_guess, &
                               TRC_OK, TRC_ERR_STATE, TRC_ERR_NULL, TRC_GUESS_GWH
   implicit none
   integer :: natm, dev, nocc
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :), dsad(:, :), dg3(:, :, :), dlib(:, :, :), elib(:, :)
   type(trc_basis_t) :: b, aux
   type(trc_pairlist_t) :: pl
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: ref
   type(trc_rimp2_result_t) :: mp2
   type(error_t) :: err
   type(c_ptr) :: h
   integer(c_int) :: rc, n, ns, it, flag
   real(c_double) :: e, e_os, e_ss
   character(kind=c_char) :: msg(256)
   logical :: ok

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   nocc = sum(zint)/2
   ok = .true.

   ! --- the reference, through the Fortran drivers, from terco's own SAD ----
   call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, b, err)
   if (.not. err%has_error()) call trc_basis_from_json('basis_sets/cc-pvdz-rifit.json', natm, zint, at_r, aux, err)
   if (.not. err%has_error()) call trc_sad_build(b, dsad, err)
   if (err%has_error()) then
      print '(a)', "context_scf: "//err%get_message()
      stop 1
   end if
   call b%to_device(); call aux%to_device()
   opts%conv_energy = 1.0e-10_dp; opts%conv_diis = 1.0e-6_dp
   allocate (dg3(b%nao, b%nao, 1)); dg3(:, :, 1) = dsad
   call trc_scf_run(b, nocc, nocc, opts, ref, dguess=dg3)
   call pl%build(b, 1.0e-12_dp); call pl%to_device()
   call trc_rimp2_run(b, aux, pl, nocc, ref%cmo(:, :, 1), ref%eps(:, 1), mp2, nfrozen=1, aux_block=20)

   ! --- misuse first: every one a status, none a crash ---------------------
   rc = trc_destroy(c_null_ptr);            call expect(rc == TRC_ERR_NULL, "destroy(null) is NULL")
   rc = trc_create(h);                      call expect(rc == TRC_OK, "create")
   rc = trc_run_scf(h);                     call expect(rc == TRC_ERR_STATE, "run before setup is STATE")
   rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
   call expect(rc == TRC_ERR_STATE, "basis before molecule is STATE")
   rc = trc_run_rimp2(h);                   call expect(rc == TRC_ERR_STATE, "rimp2 before scf is STATE")
   rc = trc_message(h, msg, 256_c_int);     call expect(rc == TRC_OK .and. msg(1) /= c_null_char, "a refusal leaves a message")
   rc = trc_set_method(h, 'pbe99'//c_null_char, 3_c_int)
   call expect(rc /= TRC_OK, "an unknown functional is refused at the setter")
   rc = trc_message(h, msg, 256_c_int);     call expect(index(transfer(msg, repeat(' ', 256)), 'pbe99') > 0, "and named")
   rc = trc_set_method(h, ''//c_null_char, -1_c_int)
   call expect(rc /= TRC_OK, "a negative grid level is refused")

   ! --- the same calculation through the handle ----------------------------
   rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
   call expect(rc == TRC_OK, "set_molecule")
   rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char);       call expect(rc == TRC_OK, "set_basis_json")
   rc = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char);   call expect(rc == TRC_OK, "set_aux_json")
   rc = trc_set_method(h, ''//c_null_char, 3_c_int);                          call expect(rc == TRC_OK, "set_method")
   rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
   call expect(rc == TRC_OK, "set_convergence")
   rc = trc_set_rimp2(h, 1_c_int, 20_c_int);                                  call expect(rc == TRC_OK, "set_rimp2")
   rc = trc_nao(h, n);                      call expect(rc == TRC_OK .and. n == b%nao, "nao")
   rc = trc_run_scf(h);                     call expect(rc == TRC_OK, "run_scf")
   rc = trc_nspin(h, ns);                   call expect(rc == TRC_OK .and. ns == 1, "nspin")
   rc = trc_converged(h, flag);             call expect(flag == 1, "converged")
   rc = trc_iterations(h, it);              call expect(it == ref%iterations, "same iteration count as the driver")
   rc = trc_energy(h, e)
   print '(a,f20.12,a,f20.12,a,es9.2)', "context_scf: E_scf handle ", e, "  driver ", ref%energy, "  diff ", e - ref%energy
   call expect(abs(e - ref%energy) < 1.0e-10_dp, "SCF energy matches the driver")
   allocate (dlib(n, n, 1), elib(n, 1))
   rc = trc_density(h, dlib); rc = trc_mo_energies(h, elib)
   call expect(maxval(abs(dlib(:, :, 1) - ref%dmat(:, :, 1))) < 1.0e-8_dp, "density matches")
   call expect(maxval(abs(elib(:, 1) - ref%eps(:, 1))) < 1.0e-8_dp, "orbital energies match")
   ! the auxiliary basis set AFTER the SCF must leave the SCF standing: that
   ! is the order a driver naturally uses
   rc = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char);   call expect(rc == TRC_OK, "aux after scf")
   rc = trc_converged(h, flag);             call expect(rc == TRC_OK .and. flag == 1, "SCF survives a new aux basis")
   rc = trc_run_rimp2(h);                   call expect(rc == TRC_OK, "run_rimp2")
   rc = trc_rimp2_energy(h, e_os, e_ss)
   print '(a,f20.12,a,f20.12,a,es9.2)', "context_scf: E_corr handle ", e_os + e_ss, "  driver ", mp2%e_os + mp2%e_ss, &
      "  diff ", e_os + e_ss - mp2%e_os - mp2%e_ss
   call expect(abs(e_os - mp2%e_os) < 1.0e-9_dp .and. abs(e_ss - mp2%e_ss) < 1.0e-9_dp, "RI-MP2 matches the driver")
   rc = trc_energy(h, e)
   call expect(abs(e - (ref%energy + mp2%e_os + mp2%e_ss)) < 1.0e-9_dp, "total includes the correlation")

   ! --- a setter after the run invalidates, and a second run from GWH agrees -
   rc = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int); call expect(rc == TRC_OK, "set_guess gwh")
   rc = trc_rimp2_energy(h, e_os, e_ss);    call expect(rc == TRC_ERR_STATE, "rimp2 result gone after a setter")
   rc = trc_run_scf(h);                     call expect(rc == TRC_OK, "second run")
   rc = trc_energy(h, e);                   call expect(abs(e - ref%energy) < 1.0e-9_dp, "GWH start, same energy")
   rc = trc_message(h, msg, 256_c_int);     call expect(rc == TRC_OK, "message")
   rc = trc_destroy(h);                     call expect(rc == TRC_OK, "destroy")

   call pl%release(); call b%release(); call aux%release()
   if (.not. ok) then
      print '(a)', "context_scf: FAIL"
      stop 1
   end if
   print '(a)', "context_scf: PASS"

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "context_scf: FAILED: "//what
         ok = .false.
      end if
   end subroutine expect

end program context_scf
