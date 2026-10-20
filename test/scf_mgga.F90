!
! Meta-GGA SCF end to end, against pyscf on the same grid: scf_rks.F90's
! harness, its 6-31G-uncontracted water and OH radical, its same-grid
! comparison and its at-density control, for the functionals that need the
! kinetic-energy density tau.
!
!   water, closed shell:  r2SCAN, TPSS, M06-L, SCAN, and the hybrids TPSSh,
!                         r2SCAN0, M06, M06-2X, restricted
!   OH radical, doublet:  r2SCAN, TPSS, M06-L, TPSSh unrestricted
!
! tau = 1/2 sum_i |grad psi_i|^2 of the total density in RKS and of each
! spin in UKS, pyscf's and libxc's convention (src/xc/trc_xc.F90's header).
! test/pyscf_refs/generate.py's gen_scf_mgga rebuilds the same basis in the
! same order and runs pyscf's RKS/UKS on terco's own grid points and
! weights, converging each case tightly and committing its energy and
! density; this program then compares terco's converged energy and,
! sharper, terco's energy AT pyscf's density (one iteration from it).
!
program scf_mgga
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz, uncontract_basis
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix3, ref_report
   use golden_io, only: golden_begin, golden_put, golden_end
   implicit none

   integer, parameter :: NCASE = 12
   character(len=8), parameter :: mol(NCASE) = [character(len=8) :: &
                                                "water", "water", "water", "water", "water", "water", "water", &
                                                "water", "oh", "oh", "oh", "oh"]
   character(len=8), parameter :: fn(NCASE) = [character(len=8) :: &
                                               "r2scan", "tpss", "m06l", "scan", "tpssh", "r2scan0", "m06", &
                                               "m062x", "r2scan", "tpss", "m06l", "tpssh"]
   integer, parameter :: nalpha(NCASE) = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
   integer, parameter :: nbeta(NCASE) = [5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4]
   integer :: ic, unit, natm, nsh, maxnp, ish
   integer, allocatable :: sh_l(:), sh_np(:), zint(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), at_z(:), at_r(:, :)
   type(trc_basis_t) :: bas
   type(trc_scf_options_t) :: opts, opts2
   type(trc_scf_result_t) :: res, res2
   character(len=8) :: current
   logical :: all_ok
   integer :: refu, ncase_ref, nfail
   character(len=16) :: tag
   character(len=32) :: gkey
   logical :: has_ref

   open (newunit=unit, file="mgga_probe.bin", access="stream", form="unformatted", status="replace")
   write (unit) int(NCASE, int64)
   all_ok = .true.
   current = ""

   ! Opened (not required) before any case runs: a missing reference must
   ! not stop this program before rks_probe.bin is written, or
   ! test/pyscf_refs/generate.py has nothing to bootstrap from. The FAIL
   ! happens below, after the dump is safely written.
   call ref_open("mgga.ref", refu, has_ref)
   call golden_begin("scf_mgga")
   nfail = 0
   if (has_ref) then
      call ref_read_int(refu, "NCASE", ncase_ref)
      if (ncase_ref /= NCASE) then
         print "(a)", "scf_mgga: pyscf reference case count does not match this build"
         stop 1
      end if
   end if
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
            if (err%has_error()) then
               print "(a)", "scf_mgga: "//err%get_message()
               stop 1
            end if
            call uncontract_basis(jbas, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
         end block
         call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, at_z, at_r, maxnp)
         call bas%to_device()
         print "(a)", ""
         print "(a,a,a,i0,a,i0,a)", "scf_mgga: ", trim(current), ", 6-31G uncontracted, ", nsh, " shells, ", bas%nao, " AOs"
      end if
      opts = trc_scf_options_t()
      opts%functional = fn(ic)
      opts%grid_level = 3
      opts%conv_energy = 1.0e-11_dp
      opts%conv_density = 1.0e-8_dp
      ! conv_density no longer gates (trc_scf_driver.F90); the commutator norm does.
      ! At its default 1e-5 scf_rks.F90's slow UKS wB97X-V case could stop up to ~2e-7
      ! above its minimum (seen under ifx), past the 1e-8 gate below, while
      ! terco's energy AT pyscf's density agreed to 2.5e-11. Converge to the
      ! gate's own precision instead.
      opts%conv_diis = 1.0e-8_dp
      ! OH's unpaired electron sits in one of a rotationally degenerate
      ! pair of pi orbitals; starting terco's OWN SCF from the same
      ! symmetry-broken density test/pyscf_refs/generate.py's gen_scf_mgga
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
            if (nspin_chk /= 2 .or. nao_chk /= bas%nao) then
               print "(a)", "scf_mgga: pyscf reference shape does not match this build"
               stop 1
            end if
            allocate (dm0(bas%nao, bas%nao, 2))
            call ref_read_matrix3(refu, "DM0_"//trim(tag), bas%nao, bas%nao, 2, dm0)
            call trc_scf_run(bas, nalpha(ic), nbeta(ic), opts, res, dguess=dm0)
         end block
      else
         call trc_scf_run(bas, nalpha(ic), nbeta(ic), opts, res)
      end if
      if (.not. res%converged) then
         print "(a,a)", "  NOT CONVERGED: ", trim(res%message)
         all_ok = .false.
      end if
      print "(a,a8,a,i1,a,i3,a,f22.12,a,f14.10)", "  ", fn(ic), " nspin ", res%nspin, " iter ", res%iterations, &
         "  E = ", res%energy, "  N(grid) = ", res%nelec_grid
      gkey = trim(mol(ic))//"_"//trim(fn(ic))
      call golden_put("E_"//trim(gkey), res%energy)

      ! Everything pyscf needs to reproduce this case.
      write (unit) int(bas%nao, int64), int(natm, int64), int(nsh, int64), int(res%nspin, int64), &
         int(nalpha(ic), int64), int(nbeta(ic), int64)
      write (unit) at_z, at_r
      do ish = 1, nsh
         write (unit) int(centre_of(ish), int64), int(sh_l(ish), int64), sh_e(1, ish)
      end do
      write (unit) fn(ic)
      write (unit) int(size(res%grid_w), int64)
      write (unit) res%grid_r, res%grid_w
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
         if (nspin_ref /= res%nspin .or. nao_ref /= bas%nao) then
            print "(a)", "scf_mgga: pyscf reference shape does not match this build"
            stop 1
         end if
         allocate (dm_ref(bas%nao, bas%nao, res%nspin))
         call ref_read_matrix3(refu, "DM_"//trim(tag), bas%nao, bas%nao, res%nspin, dm_ref)
         de = abs(res%energy - e_ref)
         ! The at-density control: terco's OWN energy at pyscf's committed,
         ! converged density -- max_iter=1 with dguess runs exactly one
         ! Roothaan iteration, whose etot (src/scf/trc_scf_driver.F90's own "energy at
         ! this density" comment) is computed from dguess BEFORE the
         ! diagonalisation updates it, i.e. exactly energy_tot(dm0)'s
         ! equivalent -- never an iterated answer.
         opts2 = opts
         opts2%max_iter = 1
         call trc_scf_run(bas, nalpha(ic), nbeta(ic), opts2, res2, dguess=dm_ref)
         dat = abs(res2%energy - e_ref)
         call golden_put("E_AT_PYSCF_D_"//trim(gkey), res2%energy)
         print "(a,a,a,es10.2,a,es10.2)", "    ", trim(name_ref), &
            ": de(converged) ", de, "  de(at pyscf D) ", dat
         ! The CONVERGED energy of an open-shell case gets a looser bound than
         ! everything else, for scf_rks.F90's reason: OH's 2-Pi hole can sit
         ! anywhere in the degenerate pi plane, and orientations that are
         ! equivalent exactly differ by grid error. The at-density check is
         ! what pins correctness, at 1e-9 for every case, as is the converged
         ! energy of every closed-shell case: measured 6.0e-11 at worst
         ! (M06-L, both), against scf_rks.F90's 1e-8.
         if (res%nspin == 2) then
            de_tol = 2.0e-6_dp
         else
            de_tol = 1.0e-9_dp
         end if
         call ref_report(nfail, de < de_tol .and. dat < 1.0e-9_dp, &
                         trim(name_ref)//": E(converged) and E(terco at pyscf D) vs pyscf", max(de, dat))
         deallocate (dm_ref)
      end block
      end if
   end do
   call golden_end()
   close (unit)
   call bas%release()
   if (has_ref) call ref_close(refu)
   if (.not. has_ref) then
      print "(a)", "FAIL: missing reference pyscf_refs/mgga.ref -- run test/pyscf_refs/generate.py " &
         //"(see test/pyscf_refs/README.md)"
      stop 1
   end if
   if (.not. all_ok .or. nfail > 0) stop 1
   print "(a)", ""
   print "(a)", "scf_mgga: all cases converged and agree with the pyscf reference; RESULT: PASS"

contains

   integer function centre_of(ish)
      integer, intent(in) :: ish
      integer :: a
      centre_of = 0
      do a = 1, natm
         if (all(abs(sh_r(:, ish) - at_r(:, a)) < 1.0e-12_dp)) centre_of = a
      end do
      if (centre_of == 0) then
         print "(a)", "scf_mgga: a shell is on no atom"
         stop 1
      end if
   end function centre_of

end program scf_mgga
