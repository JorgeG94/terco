!
! Density-fitted RHF/RKS (`--df`'s library entry, trc_scf_df.F90) against
! pyscf's own `density_fit`, same auxiliary basis (def2-universal-jkfit) on
! both sides -- water and formamide, cc-pVDZ and 6-31G*, Hartree-Fock, PBE
! and B3LYP. Converged energy only (no at-density control: see
! test/pyscf_refs/generate.py's gen_scf_df for why that is the one thing
! this leaves out of rks.ref's own pattern).
!
program scf_df
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_df, only: df_jk_t
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none

   integer, parameter :: NGROUP = 4, NMETHOD = 3
   character(len=16), parameter :: GMOL(NGROUP) = [character(len=16) :: &
                                                   "water", "water", "formamide", "formamide"]
   character(len=16), parameter :: GBASIS(NGROUP) = [character(len=16) :: &
                                                      "cc-pvdz", "6-31gs", "cc-pvdz", "6-31gs"]
   character(len=8), parameter :: METHOD(NMETHOD) = [character(len=8) :: "", "pbe", "b3lyp"]

   type(trc_basis_t) :: bas, auxb
   type(trc_pairlist_t) :: pl
   type(trc_linalg_t) :: la
   type(df_jk_t) :: dfjk
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: res
   type(error_t) :: err
   integer :: ig, im, unit, natm, nsh, nsha, maxnp, nalpha, nbeta, nelec
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   logical :: all_ok

   open (newunit=unit, file="scf_df_probe.bin", access="stream", form="unformatted", status="replace", &
         action="write")
   write (unit) int(NGROUP*NMETHOD, int64)
   all_ok = .true.

   do ig = 1, NGROUP
      call read_xyz(trim(GMOL(ig))//".xyz", natm, zint, at_r)
      nelec = sum(zint)
      nalpha = (nelec + mod(nelec, 2))/2
      nbeta = nelec - nalpha

      call trc_basis_from_json("basis_sets/"//trim(GBASIS(ig))//".json", natm, zint, at_r, bas, err, &
                               force_cart=.true.)
      if (err%has_error()) then
         print "(a)", "scf_df: "//err%get_message()
         stop 1
      end if
      call bas%to_device()
      call trc_basis_from_json("basis_sets/def2-universal-jkfit.json", natm, zint, at_r, auxb, err, &
                               force_cart=.true.)
      if (err%has_error()) then
         print "(a)", "scf_df: "//err%get_message()
         stop 1
      end if
      call auxb%to_device()

      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()
      call la%init(bas%nao)
      call dfjk%setup(bas, auxb, pl, la, err, nocc_hint=max(nalpha, 1))
      if (err%has_error()) then
         print "(a)", "scf_df: "//err%get_message()
         stop 1
      end if
      call pl%release()

      print "(a)", ""
      print "(a,a,a,a,a,i0,a)", "scf_df: ", trim(GMOL(ig)), "/", trim(GBASIS(ig)), ", ", bas%nao, " AOs"

      do im = 1, NMETHOD
         opts = trc_scf_options_t()
         opts%df = .true.
         opts%functional = METHOD(im)
         opts%grid_level = 3
         opts%conv_energy = 1.0e-11_dp
         opts%conv_diis = 1.0e-8_dp
         call trc_scf_run(bas, nalpha, nbeta, opts, res, la=la, dfjk=dfjk)
         if (.not. res%converged) then
            print "(a,a)", "  NOT CONVERGED: ", trim(res%message)
            all_ok = .false.
         end if
         if (len_trim(METHOD(im)) > 0) then
            print "(a,a8,a,i3,a,f22.12)", "  ", METHOD(im), " iter ", res%iterations, "  E = ", res%energy
         else
            print "(a,a8,a,i3,a,f22.12)", "  ", "hf", " iter ", res%iterations, "  E = ", res%energy
         end if

         write (unit) int(bas%nao, int64), int(natm, int64), int(bas%nshell, int64), int(auxb%nshell, int64), &
            int(nalpha, int64), int(nbeta, int64)
         write (unit) real(zint, dp), at_r
         call dump_shells(unit, bas)
         call dump_shells(unit, auxb)
         write (unit) METHOD(im)
         write (unit) res%energy
      end do

      call dfjk%release()
      call auxb%release()
      call bas%release()
      call la%release()
      deallocate (zint, at_r)
   end do
   close (unit)
   if (.not. all_ok) stop 1
   print "(a)", ""
   print "(a)", "scf_df: every case converged; RESULT: PASS (compare against pyscf with generate.py)"

contains

   subroutine dump_shells(u, b)
      !! One shell per record: its atom (1-based), angular momentum,
      !! primitive count, then its exponents and (already normalised, as
      !! terco stores them) contraction coefficients -- exactly the shape
      !! test/pyscf_refs/generate.py's rimp2 reader (`_gen_rimp2_one`'s
      !! `shells()`) already knows how to turn into a pyscf basis dict.
      integer, intent(in) :: u
      type(trc_basis_t), intent(in) :: b
      integer :: ish
      do ish = 1, b%nshell
         if (b%sh_at(ish) <= 0) then
            print "(a)", "scf_df: a shell is on no atom"
            stop 1
         end if
         write (u) int(b%sh_at(ish), int64), int(b%sh_l(ish), int64), int(b%sh_np(ish), int64)
         write (u) b%sh_e(1:b%sh_np(ish), ish), b%sh_c(1:b%sh_np(ish), ish)
      end do
   end subroutine dump_shells

end program scf_df
