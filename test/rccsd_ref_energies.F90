module rccsd_ref_energies
   !!
   !! The converged-energy records of test/pyscf_refs/rccsd_ints_*.ref, for the checks that drive RI-CCSD(T) through
   !! the C ABI and the CLI (check_capi_rccsd, check_cli_rccsd). The reader is sequential, so the probe records of the
   !! integral-level checks that sit between E_FIXED and E_CCSD_CORR are read and dropped; their names and lengths are
   !! those of gen_rccsd_ints in generate.py. Keeping this reader apart from rccsd_test_support leaves the hashed
   !! helper of the integral-level checks untouched.
   !!
   use pic_types, only: dp
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector
   implicit none
   private
   public :: rccsd_refs_t, read_rccsd_refs

   type :: rccsd_refs_t
      !! pyscf's DF-CCSD(T) of one molecule, cc-pVDZ / cc-pVDZ-RIFIT, converged to 1e-12 / 1e-10
      integer :: nfrozen = 0
      real(dp) :: e_scf = 0.0_dp         !! RHF total energy
      real(dp) :: e_mp2 = 0.0_dp         !! DF-MP2 correlation energy
      real(dp) :: e_ccsd_corr = 0.0_dp   !! DF-CCSD correlation energy
      real(dp) :: e_t = 0.0_dp           !! the (T) correction
      real(dp) :: e_ccsd_t = 0.0_dp      !! the CCSD(T) TOTAL energy
      real(dp) :: t1_diag = 0.0_dp
   end type rccsd_refs_t

contains

   subroutine read_rccsd_refs(name, refs)
      !! Reads `pyscf_refs/<name>`; a missing file stops the program (run generate.py)
      character(len=*), intent(in) :: name
      type(rccsd_refs_t), intent(out) :: refs

      character(len=10), parameter :: NAMES(22) = [character(len=10) :: "oooo", "ooov", "oovv", "ovov", "ovvv", "t1", "t2", &
                                                   "Foo", "Fvv", "Fov", "Loo", "Lvv", "Woooo", "Wvoov", "Wvovo", "ladder", &
                                                   "ladder_w", "dressing", "t1_new", "t2_new", "t1_conv", "t2_conv"]
      integer, parameter :: LENGTHS(22) = [4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4]
      real(dp) :: probe(4), skipped
      integer :: u, i, nskip

      call ref_require(name, u)
      call ref_read_int(u, "NO", nskip)
      call ref_read_int(u, "NV", nskip)
      call ref_read_int(u, "NAUX", nskip)
      call ref_read_int(u, "NFROZEN", refs%nfrozen)
      call ref_read_real(u, "E_SCF", refs%e_scf)
      call ref_read_real(u, "E_MP2", refs%e_mp2)
      call ref_read_real(u, "E_FIXED", skipped)
      do i = 1, size(LENGTHS)
         call ref_read_vector(u, "P_"//trim(NAMES(i)), LENGTHS(i), probe)
      end do
      call ref_read_real(u, "E_CCSD_CORR", refs%e_ccsd_corr)
      call ref_read_real(u, "E_T", refs%e_t)
      call ref_read_real(u, "E_CCSD_T", refs%e_ccsd_t)
      call ref_read_real(u, "E_UPDATE1", skipped)
      call ref_read_real(u, "T1_DIAG", refs%t1_diag)
      call ref_close(u)
   end subroutine read_rccsd_refs

end module rccsd_ref_energies
