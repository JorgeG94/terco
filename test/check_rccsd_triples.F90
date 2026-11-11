!
! The perturbative triples (T) of RI-CCSD (src/post_hf/ccsd/trc_rccsd_triples.F90, the host reference) against the
! DF-RCCSD(T) of pyscf, on a water and on a C1 hydrogen peroxide, cc-pVDZ with cc-pVDZ-RIFIT, spherical, frozen core
! 1 and 2 (the molecules, SCF and dipole blocks of check_rccsd_ints, shared through rccsd_test_support). Water is C2v;
! the C1 peroxide is the molecule in which a transposed W pairing table cannot hide behind symmetry.
!
! References (test/pyscf_refs/rccsd_ints_*.ref, gen_rccsd_ints of generate.py):
!   E_T, E_CCSD_T       pyscf ccsd_t() on the converged DF-CCSD (conv_tol 1e-12, normt 1e-10), and the total energy
!   E_T_FIXED           ccsd_t.kernel at the FIXED amplitudes of the other checks (t2 = MP2, t1_ia from the dipoles):
!                       the triples alone, independent of the CCSD solver, with a non-trivial t1
!   E_T_FIXED_NO_T1     the same with t1 = 0
!
! Cases, per molecule:
!   1 E(T) at the fixed amplitudes against E_T_FIXED
!   2 a converged `trc_rccsd_solve` (conv_energy 1e-12, conv_amplitude 1e-10), then `rccsd_triples_finish`: E(T)
!     against E_T, and E_SCF + E_CCSD + E(T) against E_CCSD_T
!   3 the fixed amplitudes against an independent host evaluation written from pyscf ccsd_t_slow einsums over ALL
!     (a,b,c), with no a >= b >= c restriction and no W pairing table: the oracle for a wrong table
!   4 t1 = 0: against the same evaluation and against E_T_FIXED_NO_T1; E(T) <= 0
!   5 the result is the same for c batches of 1, 7 and 32, and bitwise the same for 1 to 6 threads (OpenMP builds)
!   6 refusals: an empty type, mismatched amplitudes, a NaN amplitude, a non-converged state, a batch of 0
!   7 the DEVICE version (trc_rccsd_triples_dev) against the host reference to 1e-16, against pyscf with the
!     tolerances of cases 1, 2 and 4, and against the independent evaluation, at the fixed amplitudes (with and
!     without t1) and at the converged ones
!   8 tile independence: forced tiles of 1, 7 and 37 triples (37 cuts the (b,c) slabs unevenly), the default
!     budget, and a budget that holds exactly one triple, all within 1e-16
!   9 the energy of every triple of the device against the host routine's per-(a,b)-pair energies summed over c
!  10 `rccsd_triples_finish` with the device and with the reference: the same result fields
!  11 device refusals: a stale scratch after a second `ints%build`, amplitudes of another size, empty integrals, a
!     budget below one triple, a NaN amplitude, unallocated amplitudes
!  12 the occupied-driven DEVICE version (trc_rccsd_triples_ijk) against the host reference, pyscf and the
!     independent evaluation, at the fixed amplitudes (with and without t1) and at the converged ones; its i = j = k
!     triples, skipped by default, contribute round-off zeros when computed
!  13 its tiles of 1, 2 and 7 occupied triples, all of them, the default budget, a one-triple budget and a slice
!     cache of 4 (evicting) agree; a budget
!     below one triple is refused; `rccsd_triples_finish` with it agrees with the reference
! and, once, without a molecule:
!   - two electrons (one occupied orbital): E(T) = 0 for any integrals and amplitudes
!   - random integrals with the symmetries of the real ones and more virtuals than one c batch (nv = 40): the
!     routine against the independent evaluation, with batches of 32 and 5, and the device version with tiles of 1,
!     7 and the default
!   - the 36-term sum of the device item routine, unrolled, against the table-driven host sum on random W
!   - no occupied or no virtual orbitals is a refusal
!
module rccsd_host_triples
   !! E(T) as plain loops from the ccsd_t_slow of pyscf: the six W of a triple and its z = r3(W + V/2)/D, summed over
   !! ALL ordered (a,b,c) with weight 1. Of the 36 terms the sorted loop of pyscf writes per unordered triple (with
   !! the 1/2 or 1/6 of coincident indices in its denominator), only the Z_abc row is needed here: the other five
   !! rows are this row for the other five orderings of (a,b,c), which the ordered sum visits itself, and a
   !! coincident index makes the ordered visits repeat exactly as the weights say.
   use pic_types, only: dp
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   implicit none
   private
   public :: host_triples

contains

   subroutine host_triples(ints, t1, t2, e_t)
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), intent(out) :: e_t

      real(dp), allocatable :: w(:, :, :, :), v(:, :, :), z(:, :, :)
      real(dp) :: d3, s
      integer :: no, nv, a, b, c, i, j, k

      no = ints%no
      nv = ints%nv
      allocate (w(no, no, no, 6), v(no, no, no), z(no, no, no))
      e_t = 0.0_dp
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               call get_w(ints, t2, a, b, c, w(:, :, :, 1))
               call get_w(ints, t2, a, c, b, w(:, :, :, 2))
               call get_w(ints, t2, b, a, c, w(:, :, :, 3))
               call get_w(ints, t2, b, c, a, w(:, :, :, 4))
               call get_w(ints, t2, c, a, b, w(:, :, :, 5))
               call get_w(ints, t2, c, b, a, w(:, :, :, 6))
               ! v(i,j,k) = (ia|jb) t1(k,c); x = w + v/2; z = r3(x)/D with D = e_i + e_j + e_k - e_a - e_b - e_c
               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        v(i, j, k) = w(i, j, k, 1) + 0.5_dp*ints%ovov(i, a, j, b)*t1(k, c)
                     end do
                  end do
               end do
               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        d3 = ints%eps_o(i) + ints%eps_o(j) + ints%eps_o(k) - ints%eps_v(a) - ints%eps_v(b) - ints%eps_v(c)
                        z(i, j, k) = (4.0_dp*v(i, j, k) + v(k, i, j) + v(j, k, i) - 2.0_dp*v(k, j, i) &
                                      - 2.0_dp*v(i, k, j) - 2.0_dp*v(j, i, k))/d3
                     end do
                  end do
               end do
               ! et += einsum('ijk,ijk', wabc, z) + 'ikj,ijk' wacb + 'jik,ijk' wbac + 'jki,ijk' wbca + 'kij,ijk' wcab
               !       + 'kji,ijk' wcba
               s = 0.0_dp
               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        s = s + z(i, j, k)*(w(i, j, k, 1) + w(i, k, j, 2) + w(j, i, k, 3) + w(j, k, i, 4) &
                                            + w(k, i, j, 5) + w(k, j, i, 6))
                     end do
                  end do
               end do
               e_t = e_t + s
            end do
         end do
      end do
      e_t = 2.0_dp*e_t
   end subroutine host_triples

   subroutine get_w(ints, t2, a, b, c, w)
      !! pyscf get_w: einsum('if,fkj->ijk', vvov[a,b], t2T[c]) - einsum('ijm,mk->ijk', vooo[a], t2T[b,c]) with
      !! vvov[a,b][i,f] = ovvv(i,a,f,b), t2T[c,f,k,j] = t2(k,j,c,f), vooo[a][i,j,m] = ooov(j,m,i,a)
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t2(:, :, :, :)
      integer, intent(in) :: a, b, c
      real(dp), intent(out) :: w(:, :, :)

      integer :: i, j, k, f, m
      real(dp) :: s

      do k = 1, ints%no
         do j = 1, ints%no
            do i = 1, ints%no
               s = 0.0_dp
               do f = 1, ints%nv
                  s = s + ints%ovvv(i, a, f, b)*t2(k, j, c, f)
               end do
               do m = 1, ints%no
                  s = s - ints%ooov(j, m, i, a)*t2(m, k, b, c)
               end do
               w(i, j, k) = s
            end do
         end do
      end do
   end subroutine get_w

end module rccsd_host_triples

module check_rccsd_triples_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_error, only: error_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_driver, only: trc_rccsd_options_t, trc_rccsd_state_t, trc_rccsd_solve
   use trc_rccsd_triples, only: trc_rccsd_triples_options_t, rccsd_triples_reference, rccsd_triples_finish, &
                                RCCSD_TRIPLES_WMAP, rccsd_triples_energy_block, RCCSD_TRIPLES_DEVICE_ABC, &
                                RCCSD_TRIPLES_DEVICE_IJK
   use trc_rccsd_triples_ijk, only: trc_rccsd_triples_ijk_scratch_t, rccsd_triples_ijk_run, rccsd_triples_ijk_run_scratch
   use trc_rccsd_triples_dev, only: trc_rccsd_triples_dev_options_t, trc_rccsd_triples_scratch_t, rccsd_triples_dev_run, &
                                    rccsd_triples_dev_run_scratch, rccsd_triples_dev_point_energy
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: NMOL, MOL_LABEL, MOL_REF, MOL_NFROZEN, load_molecule, run_scf, mo_dipoles
   use rccsd_host_triples, only: host_triples
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector
   !$ use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private
   public :: collect_check_rccsd_triples

   integer, parameter :: NCASE = 13
   integer, parameter :: REF_E_SCF = 1, REF_E_CCSD_CORR = 2, REF_E_T = 3, REF_E_CCSD_T = 4, REF_E_T_FIXED = 5, &
                         REF_E_T_NO_T1 = 6, NREF = 6  !! indices into the scalars read from a .ref
   ! Tolerances sit about one order above what is measured (host gfortran build and a V100), so a regression at the
   ! 1e-10 level cannot pass. Measured values are in brackets.
   real(dp), parameter :: TOL_PYSCF_FIXED = 1.0e-13_dp
   !! E(T) at fixed amplitudes against pyscf [6e-15 host, 1.3e-14 GPU]: the same arithmetic on the same inputs, so the
   !! limit is summation order over about 1e6 triples
   real(dp), parameter :: TOL_PYSCF_CONV = 3.0e-12_dp
   !! E(T) at converged amplitudes against pyscf [3.0e-13]: limited by the CCSD amplitudes' own 1e-10-level
   !! convergence (conv_amplitude) in both codes, not by the (T)
   real(dp), parameter :: TOL_PYSCF_TOTAL = 3.0e-11_dp
   !! E_SCF + E_CCSD + E(T) against pyscf [4.7e-12]: the same amplitude convergence plus the SCF and the energy
   !! functional, which carry |E| of about 1e2
   real(dp), parameter :: TOL_HOST = 1.0e-14_dp
   !! E(T) against the independent all-(a,b,c) evaluation, absolute, for the molecules [8.5e-17]: only summation order
   real(dp), parameter :: TOL_HOST_REL = 1.0e-13_dp
   !! the same, relative, for the random symmetric case where |E(T)| is about 471
   real(dp), parameter :: TOL_DEV_REF = 1.0e-16_dp
   !! device E(T) against the host reference on the same inputs, absolute, molecules: the same arithmetic, only the
   !! gemm shapes and the order of the sum over (i,j,k) differ [3.5e-18 host, 3.5e-18 V100]
   real(dp), parameter :: TOL_DEV_REL = 1.0e-13_dp
   !! the same, relative, for the random symmetric case where |E(T)| is about 471
   real(dp), parameter :: TOL_DEV_TILE = 1.0e-16_dp
   !! device E(T) for another tile size, absolute: the tile only cuts the gemms [0 host, 1.7e-18 V100: cuBLAS may pick
   !! other kernels as the shapes change]
   real(dp), parameter :: TOL_DEV_PAIR = 1.0e-17_dp
   !! the energy of a pair (a,b) from the device's per-triple energies against the host routine's, largest difference
   !! [2.2e-19 host, 7.1e-20 V100]
   real(dp), parameter :: TOL_IJK_REF = 1.0e-16_dp
   !! the occupied-driven device E(T) against the host reference, absolute, molecules: the same quantity summed in
   !! another order (over occupied triples, the 36 terms regrouped) [1.7e-18 GH200]
   real(dp), parameter :: TOL_IJK_TILE = 1.0e-16_dp
   !! the occupied-driven E(T) for another tile size, absolute: the tile only cuts the work [0 GH200]
   real(dp), parameter :: TOL_BATCH = 1.0e-15_dp
   !! the same E(T) for another c batch, absolute [0]: the batch only cuts the gemms, the sums are the same

   logical :: computed = .false.
   logical :: ok_case(NCASE, NMOL)
   character(len=1200) :: msg_case(NCASE, NMOL)

contains

   subroutine collect_check_rccsd_triples(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water: (T) at fixed amplitudes vs pyscf", w1), &
                  new_unittest("water: (T) of the converged run vs pyscf", w2), &
                  new_unittest("water: (T) vs the independent all-(a,b,c) evaluation", w3), &
                  new_unittest("water: t1 = 0", w4), &
                  new_unittest("water: c batch 1, 7, 32 and thread counts", w5), &
                  new_unittest("water: refusals", w6), &
                  new_unittest("h2o2: (T) at fixed amplitudes vs pyscf", h1), &
                  new_unittest("h2o2: (T) of the converged run vs pyscf", h2), &
                  new_unittest("h2o2: (T) vs the independent all-(a,b,c) evaluation", h3), &
                  new_unittest("h2o2: t1 = 0", h4), &
                  new_unittest("h2o2: c batch 1, 7, 32 and thread counts", h5), &
                  new_unittest("h2o2: refusals", h6), &
                  new_unittest("water: device (T) vs reference, pyscf and the independent evaluation", w7), &
                  new_unittest("water: tile sizes and the one-triple budget", w8), &
                  new_unittest("water: per-triple energies vs the host pair energies", w9), &
                  new_unittest("water: finish with device and with reference", w10), &
                  new_unittest("water: device refusals", w11), &
                  new_unittest("h2o2: device (T) vs reference, pyscf and the independent evaluation", h7), &
                  new_unittest("h2o2: tile sizes and the one-triple budget", h8), &
                  new_unittest("h2o2: per-triple energies vs the host pair energies", h9), &
                  new_unittest("h2o2: finish with device and with reference", h10), &
                  new_unittest("h2o2: device refusals", h11), &
                  new_unittest("water: occupied-driven device (T) vs reference, pyscf and the independent evaluation", w12), &
                  new_unittest("water: occupied-driven tiles, budgets and finish", w13), &
                  new_unittest("h2o2: occupied-driven device (T) vs reference, pyscf and the independent evaluation", h12), &
                  new_unittest("h2o2: occupied-driven tiles, budgets and finish", h13), &
                  new_unittest("two electrons: (T) is zero", two_electrons), &
                  new_unittest("random symmetric integrals, nv above the c batch", random_symmetric), &
                  new_unittest("an empty space is a refusal", empty_space), &
                  new_unittest("device item routine: unrolled 36-term sum equals the WMAP table", unrolled_vs_table) &
                  ]
   end subroutine collect_check_rccsd_triples

   subroutine w1(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1, 1)
   end subroutine w1
   subroutine w2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 2, 1)
   end subroutine w2
   subroutine w3(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 3, 1)
   end subroutine w3
   subroutine w4(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 4, 1)
   end subroutine w4
   subroutine w5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5, 1)
   end subroutine w5
   subroutine w6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 6, 1)
   end subroutine w6
   subroutine h1(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1, 2)
   end subroutine h1
   subroutine h2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 2, 2)
   end subroutine h2
   subroutine h3(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 3, 2)
   end subroutine h3
   subroutine h4(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 4, 2)
   end subroutine h4
   subroutine h5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5, 2)
   end subroutine h5
   subroutine h6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 6, 2)
   end subroutine h6

   subroutine w7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 7, 1)
   end subroutine w7
   subroutine w8(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 8, 1)
   end subroutine w8
   subroutine w9(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 9, 1)
   end subroutine w9
   subroutine w10(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 10, 1)
   end subroutine w10
   subroutine w11(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 11, 1)
   end subroutine w11
   subroutine h7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 7, 2)
   end subroutine h7
   subroutine h8(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 8, 2)
   end subroutine h8
   subroutine h9(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 9, 2)
   end subroutine h9
   subroutine h10(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 10, 2)
   end subroutine h10
   subroutine h11(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 11, 2)
   end subroutine h11

   subroutine w12(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 12, 1)
   end subroutine w12
   subroutine w13(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 13, 1)
   end subroutine w13
   subroutine h12(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 12, 2)
   end subroutine h12
   subroutine h13(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 13, 2)
   end subroutine h13

   subroutine run_case(error, ic, im)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ic, im
      call ensure_computed()
      call check(error, ok_case(ic, im), trim(msg_case(ic, im)))
   end subroutine run_case

   subroutine ensure_computed()
      integer :: im, dev
      if (computed) return
      dev = trc_bind_device(0)
      ok_case = .false.
      msg_case = "not computed"
      do im = 1, NMOL
         call compute_molecule(im)
      end do
      computed = .true.
   end subroutine ensure_computed

   subroutine record(im, ic, passed, text)
      !! One (ok, message) per case, also printed so a run shows the measured numbers
      integer, intent(in) :: im, ic
      logical, intent(in) :: passed
      character(len=*), intent(in) :: text
      ok_case(ic, im) = passed
      msg_case(ic, im) = text
      print "(a,l1,a,a)", "  check_rccsd_triples ", passed, " ", text
   end subroutine record

   function ieee_nan() result(x)
      !! A quiet NaN without a compile-time constant folding it away
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      real(dp) :: x
      x = ieee_value(1.0_dp, ieee_quiet_nan)
   end function ieee_nan

   subroutine skip_to_energies(ref_unit)
      !! The reader is sequential: the probes of the other checks, in the order of the file, up to E_CCSD_CORR
      integer, intent(in) :: ref_unit

      character(len=10), parameter :: NAMES(22) = [character(len=10) :: "oooo", "ooov", "oovv", "ovov", "ovvv", "t1", "t2", &
                                                   "Foo", "Fvv", "Fov", "Loo", "Lvv", "Woooo", "Wvoov", "Wvovo", "ladder", &
                                                   "ladder_w", "dressing", "t1_new", "t2_new", "t1_conv", "t2_conv"]
      !! the length of each probe: 3 for a 2-index tensor, 4 for a 4-index one
      integer, parameter :: LENGTHS(22) = [4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4]
      real(dp) :: probe(4)
      integer :: i

      do i = 1, size(LENGTHS)
         call ref_read_vector(ref_unit, "P_"//trim(NAMES(i)), LENGTHS(i), probe)
      end do
   end subroutine skip_to_energies

   subroutine read_refs(ref_unit, refs)
      integer, intent(in) :: ref_unit
      real(dp), intent(out) :: refs(NREF)

      real(dp) :: skipped
      integer :: nskip

      call ref_read_int(ref_unit, "NO", nskip)
      call ref_read_int(ref_unit, "NV", nskip)
      call ref_read_int(ref_unit, "NAUX", nskip)
      call ref_read_int(ref_unit, "NFROZEN", nskip)
      call ref_read_real(ref_unit, "E_SCF", refs(REF_E_SCF))
      call ref_read_real(ref_unit, "E_MP2", skipped)
      call ref_read_real(ref_unit, "E_FIXED", skipped)
      call skip_to_energies(ref_unit)
      call ref_read_real(ref_unit, "E_CCSD_CORR", refs(REF_E_CCSD_CORR))
      call ref_read_real(ref_unit, "E_T", refs(REF_E_T))
      call ref_read_real(ref_unit, "E_CCSD_T", refs(REF_E_CCSD_T))
      call ref_read_real(ref_unit, "E_UPDATE1", skipped)
      call ref_read_real(ref_unit, "T1_DIAG", skipped)
      call ref_read_real(ref_unit, "E_T_FIXED", refs(REF_E_T_FIXED))
      call ref_read_real(ref_unit, "E_T_FIXED_NO_T1", refs(REF_E_T_NO_T1))
   end subroutine read_refs

   subroutine compute_molecule(im)
      !! SCF, the converged run, the fixed amplitudes, and every number the cases compare
      integer, intent(in) :: im

      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: la
      type(trc_rccsd_state_t) :: state
      type(trc_rccsd_options_t) :: opts
      type(trc_rccsd_amps_t) :: fixed
      integer, allocatable :: zint(:)
      real(dp), allocatable :: dip(:, :, :)
      real(dp) :: refs(NREF)
      integer :: nocc, nfz, nmo, no, nv, ref_unit
      logical :: supported
      character(len=:), allocatable :: label

      label = trim(MOL_LABEL(im))//": "
      nfz = MOL_NFROZEN(im)
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) then
         print "(a)", "SKIP check_rccsd_triples: this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
         stop 0
      end if
      nmo = bas%nao_sph
      call run_scf(bas, nocc, res)
      call la%init(2)

      opts%nocc = nocc
      opts%nfrozen = nfz
      opts%conv_energy = 1.0e-12_dp
      opts%conv_amplitude = 1.0e-10_dp
      opts%max_iterations = 200
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
      if (state%result%refused) error stop "check_rccsd_triples: the converged run was refused: "//trim(state%result%message)
      no = state%result%no
      nv = state%result%nv
      allocate (dip(no + nv, no + nv, 4))
      call mo_dipoles(bas, pl, res%cmo(:, nfz + 1:nmo, 1), dip)
      call ref_require(trim(MOL_REF(im)), ref_unit)
      call read_refs(ref_unit, refs)
      call ref_close(ref_unit)

      ! the fixed amplitudes of the other checks: t2 = MP2, t1_ia = (0.03 x + 0.02 y + 0.01 z)_ia / (e_i - e_a)
      call fixed%init(no, nv)
      call fixed_amplitudes(state%ints, dip, fixed)
      call fixed%push()  ! the device version reads the device copy

      call fixed_cases(im, label, state%ints, fixed, refs)
      call converged_case(im, label, state, refs)
      call batch_case(im, label, state%ints, fixed)
      call refusal_cases(im, label, state, fixed)
      call device_cases(im, label, state, fixed, refs)
      call device_tile_cases(im, label, state%ints, fixed)
      call device_finish_case(im, label, state)
      call ijk_cases(im, label, state, fixed, refs)
      call ijk_tile_cases(im, label, state, fixed)
      call device_refusal_cases(im, label, state, fixed, bas, aux, pl, res, nocc, nfz)
      call fixed%release()
      call state%release()
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine compute_molecule

   ! ---------------------------------------------------------------------
   ! The device version (cases 7 to 11)
   ! ---------------------------------------------------------------------

   subroutine device_cases(im, label, state, fixed, refs)
      !! Cases 7 and 9: the device E(T) against the reference, pyscf and the independent evaluation, at the fixed
      !! amplitudes with and without t1 and at the converged ones; and its per-triple energies against the pair energies
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_state_t), intent(in) :: state
      type(trc_rccsd_amps_t), intent(in) :: fixed
      real(dp), intent(in) :: refs(:)

      type(trc_rccsd_amps_t) :: no_t1
      type(error_t) :: status
      real(dp), allocatable :: e_trip(:), e_pair(:)
      real(dp) :: e_dev, e_ref, e_host, e_dev0, e_ref0, e_host0, e_devc, e_refc, secs
      real(dp) :: d_ref, d_py, d_host, d_ref0, d_py0, d_host0, d_refc, d_pyc, worst_pair, pair
      integer(int64) :: n
      integer :: a, b, c

      call rccsd_triples_dev_run(state%ints, fixed, e_dev, status, seconds=secs, triple_energies=e_trip)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      call rccsd_triples_reference(state%ints, fixed, e_ref, status, pair_energies=e_pair)
      call host_triples(state%ints, fixed%t1, fixed%t2, e_host)
      call no_t1%init(state%ints%no, state%ints%nv)
      no_t1%t2 = fixed%t2
      call no_t1%push()
      call rccsd_triples_dev_run(state%ints, no_t1, e_dev0, status)
      call rccsd_triples_reference(state%ints, no_t1, e_ref0, status)
      call host_triples(state%ints, no_t1%t1, no_t1%t2, e_host0)
      call no_t1%release()
      call rccsd_triples_dev_run(state%ints, state%amps, e_devc, status)
      call rccsd_triples_reference(state%ints, state%amps, e_refc, status)
      d_ref = abs(e_dev - e_ref)
      d_py = abs(e_dev - refs(REF_E_T_FIXED))
      d_host = abs(e_dev - e_host)
      d_ref0 = abs(e_dev0 - e_ref0)
      d_py0 = abs(e_dev0 - refs(REF_E_T_NO_T1))
      d_host0 = abs(e_dev0 - e_host0)
      d_refc = abs(e_devc - e_refc)
      d_pyc = abs(e_devc - refs(REF_E_T))
      call record(im, 7, d_ref < TOL_DEV_REF .and. d_py < TOL_PYSCF_FIXED .and. d_host < TOL_HOST .and. &
                  d_ref0 < TOL_DEV_REF .and. d_py0 < TOL_PYSCF_FIXED .and. d_host0 < TOL_HOST .and. &
                  d_refc < TOL_DEV_REF .and. d_pyc < TOL_PYSCF_CONV .and. e_dev < 0.0_dp .and. e_dev0 < 0.0_dp .and. &
                  abs(e_dev - e_dev0) > 1.0e-6_dp, &
                  label//"device E(T) fixed "//to_string(e_dev, "f20.14")//": |reference| "//to_string(d_ref, "es9.2")// &
                  ", |pyscf| "//to_string(d_py, "es9.2")//", |independent| "//to_string(d_host, "es9.2")//"; t1 = 0 "// &
                  to_string(e_dev0, "f20.14")//": "//to_string(d_ref0, "es9.2")//", "//to_string(d_py0, "es9.2")//", "// &
                  to_string(d_host0, "es9.2")//"; converged "//to_string(e_devc, "f20.14")//": |reference| "// &
                  to_string(d_refc, "es9.2")//", |pyscf| "//to_string(d_pyc, "es9.2")//"; wall "//to_string(secs, "f9.3")//" s")

      ! the energy of pair (a,b) is the sum over c of the triples' (the order of the host routine's e_pair)
      worst_pair = 0.0_dp
      n = 0
      do a = 1, state%ints%nv
         do b = 1, a
            pair = 0.0_dp
            do c = 1, b
               n = n + 1
               pair = pair + e_trip(n)
            end do
            worst_pair = max(worst_pair, abs(pair - e_pair((a*(a - 1))/2 + b)))
         end do
      end do
      call record(im, 9, worst_pair < TOL_DEV_PAIR .and. n == size(e_trip, kind=int64) .and. &
                  n == int(state%ints%nv, int64)*(state%ints%nv + 1)*(state%ints%nv + 2)/6 .and. &
                  abs(sum(e_pair)*2.0_dp - e_ref) < TOL_DEV_REF, &
                  label//to_string(n)//" triples in "//to_string(size(e_pair))//" pairs: the largest difference of a pair "// &
                  "energy, device against reference, is "//to_string(worst_pair, "es9.2"))
   end subroutine device_cases

   subroutine device_tile_cases(im, label, ints, fixed)
      !! Case 8: the tile changes how the products are cut, not the energy
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: fixed

      integer(int64), parameter :: TILES(5) = [1_int64, 7_int64, 37_int64, 100000_int64, 0_int64]
      type(trc_rccsd_triples_dev_options_t) :: dopts
      type(trc_rccsd_triples_scratch_t) :: ts
      type(trc_linalg_t) :: la
      type(error_t) :: status
      real(dp) :: e(size(TILES)), e_default, e_one, e_again, spread
      integer(int64) :: ncap(size(TILES)), ntiles(size(TILES)), ntotal
      integer :: k

      call la%init(1)
      ntotal = int(ints%nv, int64)*(ints%nv + 1)*(ints%nv + 2)/6
      do k = 1, size(TILES)
         dopts%tile_triples = TILES(k)
         call ts%init(ints, dopts)
         if (ts%status%has_error()) error stop "check_rccsd_triples: "//ts%status%get_message()
         call rccsd_triples_dev_run_scratch(ints, fixed, ts, la, e(k), status)
         if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
         ncap(k) = ts%ncap
         ntiles(k) = ts%ntiles
         call ts%release()
      end do
      e_default = e(size(TILES))
      ! a budget that holds exactly one triple (6 no^3 doubles), and a second call with the same plan
      dopts%tile_triples = 0
      dopts%budget_bytes = 48_int64*int(ints%no, int64)**3
      call ts%init(ints, dopts)
      if (ts%status%has_error()) error stop "check_rccsd_triples: "//ts%status%get_message()
      call rccsd_triples_dev_run_scratch(ints, fixed, ts, la, e_one, status)
      call rccsd_triples_dev_run_scratch(ints, fixed, ts, la, e_again, status)
      k = int(ts%ncap)
      call ts%release()
      call la%release()
      spread = max(maxval(abs(e - e_default)), abs(e_one - e_default))
      call record(im, 8, spread < TOL_DEV_TILE .and. ncap(1) == 1 .and. ntiles(1) == ntotal .and. ncap(2) == 7 .and. &
                  ncap(3) == 37 .and. ntiles(3) == (ntotal + 36)/37 .and. ncap(4) == ntotal .and. ntiles(4) == 1 .and. &
                  k == 1 .and. e_again == e_one, &
                  label//"tiles of 1, 7, 37, all and the default budget ("//to_string(ncap(5))//" triples, "// &
                  to_string(ntiles(5))//" tile) and a one-triple budget: E(T) "//to_string(e_default, "f20.14")// &
                  ", largest difference "//to_string(spread, "es9.2")//"; a second call with the same scratch differs by "// &
                  to_string(abs(e_again - e_one), "es9.2"))
   end subroutine device_tile_cases

   subroutine device_finish_case(im, label, state)
      !! Case 10: `rccsd_triples_finish` with the device (set explicitly: the default follows the build) and with the reference
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_state_t), intent(inout) :: state

      type(trc_rccsd_triples_options_t) :: topts
      type(error_t) :: status
      real(dp) :: e_dev, e_total_dev, e_ref, e_total_ref
      logical :: flag_dev, flag_ref

      state%result%has_triples = .false.
      topts%use_device = .true.
      topts%device_routine = RCCSD_TRIPLES_DEVICE_ABC
      call rccsd_triples_finish(state, status, topts)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      e_dev = state%result%e_t
      e_total_dev = state%result%e_ccsd_t
      flag_dev = state%result%has_triples .and. state%result%seconds_triples > 0.0_dp
      state%result%has_triples = .false.
      topts%use_device = .false.
      call rccsd_triples_finish(state, status, topts)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      e_ref = state%result%e_t
      e_total_ref = state%result%e_ccsd_t
      flag_ref = state%result%has_triples .and. state%result%seconds_triples > 0.0_dp
      ! The sums e_ccsd + e_t (about -0.2 to -0.4) are held to a few units in the last place, not to
      ! TOL_DEV_REF: one rounding of the sum is already 5.5e-17.
      call record(im, 10, flag_dev .and. flag_ref .and. abs(e_dev - e_ref) < TOL_DEV_REF .and. &
                  abs(e_total_dev - e_total_ref) <= 4.0_dp*spacing(e_total_ref) .and. e_dev /= 0.0_dp .and. &
                  abs(e_total_dev - state%result%e_ccsd - e_dev) < 1.0e-15_dp, &
                  label//"finish: E(T) device "//to_string(e_dev, "f20.14")//", reference "//to_string(e_ref, "f20.14")// &
                  ", |difference| "//to_string(abs(e_dev - e_ref), "es9.2")//", E_CCSD(T) |difference| "// &
                  to_string(abs(e_total_dev - e_total_ref), "es9.2"))
   end subroutine device_finish_case

   ! ---------------------------------------------------------------------
   ! The occupied-driven device version (cases 12 and 13)
   ! ---------------------------------------------------------------------

   subroutine ijk_cases(im, label, state, fixed, refs)
      !! Case 12: the occupied-driven E(T) against the reference, pyscf and the independent evaluation, at the fixed
      !! amplitudes with and without t1 and at the converged ones; the i = j = k triples give round-off zeros
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_state_t), intent(in) :: state
      type(trc_rccsd_amps_t), intent(in) :: fixed
      real(dp), intent(in) :: refs(:)

      type(trc_rccsd_amps_t) :: no_t1
      type(trc_rccsd_triples_dev_options_t) :: dopts
      type(error_t) :: status
      real(dp), allocatable :: e_trip(:), e_trip_diag(:)
      real(dp) :: e_ijk, e_ijk_diag, e_ref, e_host, e_ijk0, e_ref0, e_host0, e_ijkc, e_refc, total, diag
      real(dp) :: d_ref, d_py, d_host, d_ref0, d_py0, d_host0, d_refc, d_pyc
      integer(int64) :: n
      integer :: i, j, k

      call rccsd_triples_ijk_run(state%ints, fixed, e_ijk, status, triple_energies=e_trip)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      ! the i = j = k triples are skipped by default (their energy is zero); computed on request they are round-off zeros
      dopts%diagonal_triples = .true.
      call rccsd_triples_ijk_run(state%ints, fixed, e_ijk_diag, status, dopts, triple_energies=e_trip_diag)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      call rccsd_triples_reference(state%ints, fixed, e_ref, status)
      call host_triples(state%ints, fixed%t1, fixed%t2, e_host)
      call no_t1%init(state%ints%no, state%ints%nv)
      no_t1%t2 = fixed%t2
      call no_t1%push()
      call rccsd_triples_ijk_run(state%ints, no_t1, e_ijk0, status)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      call rccsd_triples_reference(state%ints, no_t1, e_ref0, status)
      call host_triples(state%ints, no_t1%t1, no_t1%t2, e_host0)
      call no_t1%release()
      call rccsd_triples_ijk_run(state%ints, state%amps, e_ijkc, status)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      call rccsd_triples_reference(state%ints, state%amps, e_refc, status)
      ! the per-triple energies add up to E(T) in their order; the i = j = k ones are zero up to round-off
      total = 0.0_dp
      diag = 0.0_dp
      n = 0
      do i = 1, state%ints%no
         do j = 1, i
            do k = 1, j
               n = n + 1
               total = total + e_trip(n)
               if (i == k) diag = max(diag, abs(e_trip_diag(n)) + abs(e_trip(n)))
            end do
         end do
      end do
      d_ref = abs(e_ijk - e_ref)
      d_py = abs(e_ijk - refs(REF_E_T_FIXED))
      d_host = abs(e_ijk - e_host)
      d_ref0 = abs(e_ijk0 - e_ref0)
      d_py0 = abs(e_ijk0 - refs(REF_E_T_NO_T1))
      d_host0 = abs(e_ijk0 - e_host0)
      d_refc = abs(e_ijkc - e_refc)
      d_pyc = abs(e_ijkc - refs(REF_E_T))
      call record(im, 12, d_ref < TOL_IJK_REF .and. d_py < TOL_PYSCF_FIXED .and. d_host < TOL_HOST .and. &
                  d_ref0 < TOL_IJK_REF .and. d_py0 < TOL_PYSCF_FIXED .and. d_host0 < TOL_HOST .and. &
                  d_refc < TOL_IJK_REF .and. d_pyc < TOL_PYSCF_CONV .and. e_ijk < 0.0_dp .and. e_ijk0 < 0.0_dp .and. &
                  abs(e_ijk - e_ijk0) > 1.0e-6_dp .and. n == size(e_trip, kind=int64) .and. abs(total - e_ijk) < 1.0e-16_dp .and. &
                  diag < 1.0e-17_dp .and. abs(e_ijk_diag - e_ijk) < TOL_IJK_TILE, &
                  label//"ijk E(T) fixed "//to_string(e_ijk, "f20.14")//": |reference| "//to_string(d_ref, "es9.2")// &
                  ", |pyscf| "//to_string(d_py, "es9.2")//", |independent| "//to_string(d_host, "es9.2")//"; t1 = 0 "// &
                  to_string(e_ijk0, "f20.14")//": "//to_string(d_ref0, "es9.2")//", "//to_string(d_py0, "es9.2")//", "// &
                  to_string(d_host0, "es9.2")//"; converged "//to_string(e_ijkc, "f20.14")//": |reference| "// &
                  to_string(d_refc, "es9.2")//", |pyscf| "//to_string(d_pyc, "es9.2")//"; largest i = j = k triple (computed "// &
                  "on request) "//to_string(diag, "es9.2")//", E(T) with them "//to_string(abs(e_ijk_diag - e_ijk), "es9.2"))
   end subroutine ijk_cases

   subroutine ijk_tile_cases(im, label, state, fixed)
      !! Case 13: tiles, budgets and the finish with the occupied-driven routine
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_state_t), intent(inout) :: state
      type(trc_rccsd_amps_t), intent(in) :: fixed

      integer(int64), parameter :: TILES(5) = [1_int64, 2_int64, 7_int64, 100000_int64, 0_int64]
      type(trc_rccsd_triples_dev_options_t) :: dopts
      type(trc_rccsd_triples_ijk_scratch_t) :: ts
      type(trc_rccsd_triples_options_t) :: topts
      type(trc_linalg_t) :: la
      type(error_t) :: status, st_small
      real(dp) :: e(size(TILES)), e_one, e_again, e_small, spread, e_fin, e_ref, e_lru
      integer :: nslot_lru, nslot_default
      integer(int64) :: ncap(size(TILES)), ntiles(size(TILES)), ntotal, per_triple
      integer :: k, ncap_one
      logical :: ok

      call la%init(1)
      ntotal = int(state%ints%no, int64)*(state%ints%no + 1)*(state%ints%no + 2)/6
      per_triple = 24_int64*int(state%ints%nv, int64)**3
      do k = 1, size(TILES)
         dopts%tile_triples = TILES(k)
         call ts%init(state%ints, dopts)
         if (ts%status%has_error()) error stop "check_rccsd_triples: "//ts%status%get_message()
         call rccsd_triples_ijk_run_scratch(state%ints, fixed, ts, la, e(k), status)
         if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
         ncap(k) = ts%ncap
         ntiles(k) = ts%ntiles
         call ts%release()
      end do
      ! a budget that holds exactly one occupied triple (3 nv^3 doubles), and a second call with the same plan
      dopts%tile_triples = 0
      dopts%budget_bytes = per_triple
      call ts%init(state%ints, dopts)
      if (ts%status%has_error()) error stop "check_rccsd_triples: "//ts%status%get_message()
      call rccsd_triples_ijk_run_scratch(state%ints, fixed, ts, la, e_one, status)
      call rccsd_triples_ijk_run_scratch(state%ints, fixed, ts, la, e_again, status)
      ncap_one = int(ts%ncap)
      call ts%release()
      ! a slice cache of four (the least it takes), so slices are evicted and gathered again, against the default
      dopts%budget_bytes = 0
      dopts%slice_bytes = 1
      call la%init(1)
      call ts%init(state%ints, dopts)
      if (ts%status%has_error()) error stop "check_rccsd_triples: "//ts%status%get_message()
      call rccsd_triples_ijk_run_scratch(state%ints, fixed, ts, la, e_lru, status)
      nslot_lru = ts%nslot
      call ts%release()
      dopts%slice_bytes = 0
      call ts%init(state%ints, dopts)
      nslot_default = ts%nslot
      call ts%release()
      call la%release()
      ! a budget one byte short of a triple
      dopts%budget_bytes = per_triple - 1_int64
      call rccsd_triples_ijk_run(state%ints, fixed, e_small, st_small, dopts)
      ! the finish with the occupied-driven routine against the reference
      state%result%has_triples = .false.
      topts%use_device = .true.
      topts%device_routine = RCCSD_TRIPLES_DEVICE_IJK
      call rccsd_triples_finish(state, status, topts)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      e_fin = state%result%e_t
      ok = state%result%has_triples
      call rccsd_triples_reference(state%ints, state%amps, e_ref, status)
      spread = max(maxval(abs(e - e(size(TILES)))), abs(e_one - e(size(TILES))), abs(e_lru - e(size(TILES))))
      ! a tile is a run of k at one (i,j), cut at capacity: the run of (i,j) has j triples
      ok = ok .and. spread < TOL_IJK_TILE .and. ncap(1) == 1 .and. ntiles(1) == ntotal .and. ncap(2) == min(2_int64, ntotal) &
           .and. ncap(3) == min(7_int64, ntotal) .and. ntiles(3) == run_tiles(state%ints%no, ncap(3)) .and. &
           ncap(4) == ntotal .and. ntiles(4) == run_tiles(state%ints%no, ncap(4)) .and. &
           ntiles(2) == run_tiles(state%ints%no, ncap(2)) .and. ncap_one == 1 .and. e_again == e_one .and. &
           st_small%has_error() .and. e_small == 0.0_dp .and. index(st_small%get_message(), "one occupied triple") > 0 .and. &
           abs(e_fin - e_ref) < TOL_IJK_REF .and. nslot_lru == 4 .and. nslot_default == 2*state%ints%no
      call record(im, 13, ok, label//"ijk tiles of 1, 2, 7, all ("//to_string(ntotal)//" triples) and the default ("// &
                  to_string(ncap(5))//" a tile), a one-triple budget and a cache of 4 slices (default "// &
                  to_string(nslot_default)//"): E(T) "//to_string(e(size(TILES)), "f20.14")// &
                  ", largest difference "//to_string(spread, "es9.2")//"; a budget one byte short: ["// &
                  st_small%get_message()//"]; finish against the reference "//to_string(abs(e_fin - e_ref), "es9.2"))
   end subroutine ijk_tile_cases

   pure function run_tiles(no, cap) result(n)
      !! Tiles of the occupied-driven routine for tiles of up to `cap` triples: ceil(j / cap) for each pair i >= j
      integer, intent(in) :: no
      integer(int64), intent(in) :: cap
      integer(int64) :: n
      integer :: i, j

      n = 0
      do i = 1, no
         do j = 1, i
            n = n + (j + cap - 1)/cap
         end do
      end do
   end function run_tiles

   subroutine device_refusal_cases(im, label, state, fixed, bas, aux, pl, res, nocc, nfz)
      !! Case 11: every device refusal returns E(T) = 0 with a message; a scratch planned for an earlier build of the
      !! integrals is refused after a second `ints%build`
      integer, intent(in) :: im, nocc, nfz
      character(len=*), intent(in) :: label
      type(trc_rccsd_state_t), intent(inout) :: state
      type(trc_rccsd_amps_t), intent(in) :: fixed
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res

      type(trc_rccsd_ints_t) :: empty
      type(trc_rccsd_amps_t) :: bad
      type(trc_rccsd_triples_dev_options_t) :: dopts
      type(trc_rccsd_triples_scratch_t) :: ts
      type(trc_rccsd_ints_options_t) :: bopts
      type(trc_linalg_t) :: la
      type(error_t) :: st(8)
      real(dp) :: e(8), e_before, e_after
      logical :: ok
      integer :: k
      character(len=1200) :: text

      call rccsd_triples_dev_run(state%ints, fixed, e_before, st(8))
      ! 1 empty integrals
      call rccsd_triples_dev_run(empty, fixed, e(1), st(1))
      ! 2 amplitudes of another size
      call bad%init(state%ints%no, state%ints%nv + 1)
      call rccsd_triples_dev_run(state%ints, bad, e(2), st(2))
      call bad%release()
      ! 3 a NaN amplitude on the device: reported, not returned as a number
      call bad%init(state%ints%no, state%ints%nv)
      bad%t1 = fixed%t1
      bad%t2 = fixed%t2
      bad%t2(1, 1, 1, 1) = ieee_nan()
      call bad%push()
      call rccsd_triples_dev_run(state%ints, bad, e(3), st(3))
      call bad%release()
      ! 4 a budget that does not hold one triple
      dopts%budget_bytes = 48_int64*int(state%ints%no, int64)**3 - 1_int64
      call rccsd_triples_dev_run(state%ints, fixed, e(4), st(4), dopts)
      ! 5 amplitudes that were never allocated
      call rccsd_triples_dev_run(state%ints, bad, e(5), st(5))
      ! 6 a negative tile size
      dopts%budget_bytes = 0
      dopts%tile_triples = -1
      call rccsd_triples_dev_run(state%ints, fixed, e(6), st(6), dopts)
      ! 7 a scratch planned before a second build of the integrals
      dopts%tile_triples = 0
      call la%init(1)
      call ts%init(state%ints, dopts)
      bopts%nalpha = nocc
      bopts%nbeta = nocc
      bopts%nfrozen = nfz
      call state%ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      call rccsd_triples_dev_run_scratch(state%ints, fixed, ts, la, e(7), st(7))
      call ts%release()
      call la%release()
      call rccsd_triples_dev_run(state%ints, fixed, e_after, st(8))
      ok = .true.
      text = ""
      do k = 1, 7
         ok = ok .and. st(k)%has_error() .and. e(k) == 0.0_dp
         if (st(k)%has_error()) text = trim(text)//" ["//st(k)%get_message()//"]"
      end do
      ok = ok .and. index(st(3)%get_message(), "NaN") > 0 .and. index(st(4)%get_message(), "one triple") > 0 .and. &
           index(st(7)%get_message(), "stale") > 0 .and. .not. st(8)%has_error() .and. abs(e_after - e_before) < TOL_DEV_REF
      call record(im, 11, ok, label//"device refusals (empty integrals, amplitudes of another size, a NaN amplitude, a "// &
                  "budget below one triple, unallocated amplitudes, a negative tile, a stale scratch after a second "// &
                  "build; the rebuilt integrals give E(T) back to "//to_string(abs(e_after - e_before), "es9.2")//"): "// &
                  trim(text))
   end subroutine device_refusal_cases

   subroutine fixed_amplitudes(ints, dip, amps)
      !! t2(i,j,a,b) = (ia|jb) / (e_i + e_j - e_a - e_b); t1 from the dipole blocks
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: dip(:, :, :)
      type(trc_rccsd_amps_t), intent(inout) :: amps

      integer :: i, j, a, b, no

      no = ints%no
      do b = 1, ints%nv
         do a = 1, ints%nv
            do j = 1, no
               do i = 1, no
                  amps%t2(i, j, a, b) = ints%ovov(i, a, j, b)/(ints%eps_o(i) + ints%eps_o(j) - ints%eps_v(a) - ints%eps_v(b))
               end do
            end do
         end do
      end do
      do a = 1, ints%nv
         do i = 1, no
            amps%t1(i, a) = (0.03_dp*dip(i, no + a, 1) + 0.02_dp*dip(i, no + a, 2) + 0.01_dp*dip(i, no + a, 3))/ &
                            (ints%eps_o(i) - ints%eps_v(a))
         end do
      end do
   end subroutine fixed_amplitudes

   subroutine fixed_cases(im, label, ints, fixed, refs)
      !! Cases 1, 3 and 4: the fixed amplitudes against pyscf, against the independent evaluation, and with t1 = 0
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: fixed
      real(dp), intent(in) :: refs(:)

      type(trc_rccsd_amps_t) :: no_t1
      type(error_t) :: status
      real(dp) :: e_t, e_host, e_t0, e_host0, secs, d_py, d_host, d_py0, d_host0

      call rccsd_triples_reference(ints, fixed, e_t, status, seconds=secs)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      d_py = abs(e_t - refs(REF_E_T_FIXED))
      call record(im, 1, d_py < TOL_PYSCF_FIXED .and. e_t < 0.0_dp, label//"E(T) at the fixed amplitudes "// &
                  to_string(e_t, "f20.14")//", pyscf "//to_string(refs(REF_E_T_FIXED), "f20.14")//", |difference| "// &
                  to_string(d_py, "es9.2")//", wall "//to_string(secs, "f9.3")//" s")

      call host_triples(ints, fixed%t1, fixed%t2, e_host)
      d_host = abs(e_t - e_host)
      call record(im, 3, d_host < TOL_HOST .and. abs(e_host) > 1.0e-5_dp, label//"E(T) "//to_string(e_t, "f20.14")// &
                  ", independent all-(a,b,c) evaluation "//to_string(e_host, "f20.14")//", |difference| "// &
                  to_string(d_host, "es9.2"))

      call no_t1%init(ints%no, ints%nv)
      no_t1%t2 = fixed%t2
      call rccsd_triples_reference(ints, no_t1, e_t0, status)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      call host_triples(ints, no_t1%t1, no_t1%t2, e_host0)
      d_host0 = abs(e_t0 - e_host0)
      d_py0 = abs(e_t0 - refs(REF_E_T_NO_T1))
      call record(im, 4, d_host0 < TOL_HOST .and. d_py0 < TOL_PYSCF_FIXED .and. e_t0 <= 0.0_dp .and. abs(e_t0 - e_t) > 1.0e-6_dp, &
                  label//"t1 = 0: E(T) "//to_string(e_t0, "f20.14")//", independent "//to_string(e_host0, "f20.14")// &
                  " (|difference| "//to_string(d_host0, "es9.2")//"), pyscf "//to_string(refs(REF_E_T_NO_T1), "f20.14")// &
                  " (|difference| "//to_string(d_py0, "es9.2")//"); the disconnected term moves E(T) by "// &
                  to_string(e_t - e_t0, "es9.2"))
      call no_t1%release()
   end subroutine fixed_cases

   subroutine converged_case(im, label, state, refs)
      !! Case 2: `rccsd_triples_finish` on the converged state, against pyscf E_T and E_CCSD_T
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_state_t), intent(inout) :: state
      real(dp), intent(in) :: refs(:)

      type(error_t) :: status
      type(trc_rccsd_triples_options_t) :: topts
      real(dp) :: d_t, d_total, e_total
      logical :: has_before

      topts%verbose = .true.
      has_before = state%result%has_triples
      call rccsd_triples_finish(state, status, topts)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      d_t = abs(state%result%e_t - refs(REF_E_T))
      e_total = refs(REF_E_SCF) + state%result%e_ccsd_t
      d_total = abs(e_total - refs(REF_E_CCSD_T))
      call record(im, 2, .not. has_before .and. state%result%has_triples .and. state%result%converged .and. &
                  d_t < TOL_PYSCF_CONV .and. d_total < TOL_PYSCF_TOTAL .and. &
                  abs(state%result%e_ccsd_t - state%result%e_ccsd - state%result%e_t) < 1.0e-15_dp, &
                  label//"converged in "//to_string(state%result%iterations)//" iterations, E(T) "// &
                  to_string(state%result%e_t, "f20.14")//", pyscf "//to_string(refs(REF_E_T), "f20.14")//", |difference| "// &
                  to_string(d_t, "es9.2")//"; E_SCF + E_CCSD + E(T) "//to_string(e_total, "f20.12")//" vs pyscf "// &
                  to_string(refs(REF_E_CCSD_T), "f20.12")//", |difference| "//to_string(d_total, "es9.2")//"; (T) wall "// &
                  to_string(state%result%seconds_triples, "f9.3")//" s")
   end subroutine converged_case

   subroutine batch_case(im, label, ints, fixed)
      !! Case 5: the batch over c changes how the products are cut, not the energy
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: fixed

      integer, parameter :: BATCHES(4) = [1, 7, 32, 100000]  !! the last is far above nv: it must not over-allocate or fail
      type(trc_rccsd_triples_options_t) :: topts
      type(error_t) :: status
      real(dp) :: e(4), e_again, spread, e_threads(6)
      integer :: k, nt, nt_saved

      do k = 1, size(BATCHES)
         topts%c_batch = BATCHES(k)
         call rccsd_triples_reference(ints, fixed, e(k), status, topts)
         if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      end do
      call rccsd_triples_reference(ints, fixed, e_again, status, topts)
      spread = maxval(e) - minval(e)
      ! The pairs are summed in a fixed order, so the energy is bitwise the same for any thread count. That relies on
      ! the BLAS returning thread-invariant results for these gemms (true with MKL gnu_thread, as measured), and it is
      ! only meaningful in an OpenMP build: without OpenMP these are six identical serial calls.
      nt_saved = 1
      !$ nt_saved = omp_get_max_threads()
      topts%c_batch = 32
      do nt = 1, size(e_threads)
         !$ call omp_set_num_threads(nt)
         call rccsd_triples_reference(ints, fixed, e_threads(nt), status, topts)
      end do
      !$ call omp_set_num_threads(nt_saved)
      call record(im, 5, spread < TOL_BATCH .and. e_again == e(3) .and. all(e_threads == e(3)), &
                  label//"c batch 1, 7, 32, 100000 and 1 to 6 threads (bitwise, worst "// &
                  to_string(maxval(abs(e_threads - e(3))), "es9.2")//"): E(T) "// &
                  to_string(e(1), "f20.14")//" "//to_string(e(2), "f20.14")//" "//to_string(e(3), "f20.14")//" "// &
                  to_string(e(4), "f20.14")// &
                  ", spread "//to_string(spread, "es9.2")//"; a second call with batch 32 differs by "// &
                  to_string(abs(e_again - e(3)), "es9.2"))
   end subroutine batch_case

   subroutine refusal_cases(im, label, state, fixed)
      !! Case 6: every refusal returns E(T) = 0 with a message; a non-converged state is not finished
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_state_t), intent(inout) :: state
      type(trc_rccsd_amps_t), intent(in) :: fixed

      type(trc_rccsd_ints_t) :: empty
      type(trc_rccsd_amps_t) :: bad
      type(trc_rccsd_triples_options_t) :: topts
      type(trc_rccsd_state_t) :: stopped
      type(error_t) :: st(6)
      real(dp) :: e(6)
      logical :: ok, has_ok
      integer :: k
      character(len=1200) :: text

      ! 1 empty integrals
      call rccsd_triples_reference(empty, fixed, e(1), st(1))
      ! 2 amplitudes of another size
      call bad%init(state%ints%no, state%ints%nv + 1)
      call rccsd_triples_reference(state%ints, bad, e(2), st(2))
      call bad%release()
      ! 3 a NaN amplitude: reported, not returned as a number
      call bad%init(state%ints%no, state%ints%nv)
      bad%t1 = fixed%t1
      bad%t2 = fixed%t2
      bad%t2(1, 1, 1, 1) = ieee_nan()
      call rccsd_triples_reference(state%ints, bad, e(3), st(3))
      call bad%release()
      ! 4 a batch of zero
      topts%c_batch = 0
      call rccsd_triples_reference(state%ints, fixed, e(4), st(4), topts)
      ! 5 amplitudes that were never allocated
      call rccsd_triples_reference(state%ints, bad, e(5), st(5))
      ! 6 (T) of a state that did not converge
      stopped%result%converged = .false.
      call rccsd_triples_finish(stopped, st(6))
      e(6) = stopped%result%e_t
      has_ok = .not. stopped%result%has_triples
      ok = .true.
      text = ""
      do k = 1, 6
         ok = ok .and. st(k)%has_error() .and. e(k) == 0.0_dp
         if (st(k)%has_error()) text = trim(text)//" ["//st(k)%get_message()//"]"
      end do
      ok = ok .and. index(st(3)%get_message(), "NaN") > 0 .and. has_ok .and. &
           index(st(6)%get_message(), "did not converge") > 0 .and. index(st(6)%get_message(), "refused") == 0 .and. &
           index(st(6)%get_message(), "non-finite") == 0
      call record(im, 6, ok .and. state%result%e_t /= 0.0_dp, label//"refusals (empty integrals, amplitudes of another size, "// &
                  "a NaN amplitude, a batch of 0, unallocated amplitudes, a state that did not converge): "//trim(text))
   end subroutine refusal_cases

   ! ---------------------------------------------------------------------
   ! Without a molecule
   ! ---------------------------------------------------------------------

   function lcg(seed) result(x)
      !! A deterministic number in (-1, 1)
      integer(int64), intent(inout) :: seed
      real(dp) :: x
      seed = mod(seed*1103515245_int64 + 12345_int64, 2147483648_int64)
      x = 2.0_dp*real(seed, dp)/2147483648.0_dp - 1.0_dp
   end function lcg

   subroutine random_problem(no, nv, seed, symmetric, ints, amps)
      !! Random blocks, orbital energies in canonical order (occupied below -0.5, virtual above 0.3). With
      !! `symmetric`, the blocks have the permutation symmetry of real integrals: (ij|ka) = (ji|ka),
      !! (ia|jb) = (jb|ia), (ia|bc) = (ia|cb).
      integer, intent(in) :: no, nv
      integer(int64), intent(inout) :: seed
      logical, intent(in) :: symmetric
      type(trc_rccsd_ints_t), intent(out) :: ints
      type(trc_rccsd_amps_t), intent(out) :: amps

      integer :: i, j, k, a, b, c

      ints%no = no
      ints%nv = nv
      ints%stamp = 1
      allocate (ints%ooov(no, no, no, nv), ints%ovov(no, nv, no, nv), ints%ovvv(no, nv, nv, nv))
      allocate (ints%eps_o(no), ints%eps_v(nv))
      do i = 1, no
         ints%eps_o(i) = -2.0_dp + 1.4_dp*real(i - 1, dp)/real(max(no, 2) - 1, dp) + 0.05_dp*lcg(seed)
      end do
      do a = 1, nv
         ints%eps_v(a) = 0.3_dp + 2.0_dp*real(a - 1, dp)/real(max(nv, 2) - 1, dp) + 0.05_dp*lcg(seed)
      end do
      do a = 1, nv
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  ints%ooov(i, j, k, a) = 0.1_dp*lcg(seed)
               end do
            end do
         end do
      end do
      do b = 1, nv
         do j = 1, no
            do a = 1, nv
               do i = 1, no
                  ints%ovov(i, a, j, b) = 0.1_dp*lcg(seed)
               end do
            end do
         end do
      end do
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               do i = 1, no
                  ints%ovvv(i, a, b, c) = 0.1_dp*lcg(seed)
               end do
            end do
         end do
      end do
      if (symmetric) then
         do a = 1, nv
            do k = 1, no
               do j = 1, no
                  do i = 1, j
                     ints%ooov(j, i, k, a) = ints%ooov(i, j, k, a)
                  end do
               end do
            end do
         end do
         do b = 1, nv
            do j = 1, no
               do a = 1, nv
                  do i = 1, no
                     if (j < i .or. (j == i .and. b < a)) ints%ovov(i, a, j, b) = ints%ovov(j, b, i, a)
                  end do
               end do
            end do
         end do
         do c = 1, nv
            do b = 1, c
               do a = 1, nv
                  do i = 1, no
                     ints%ovvv(i, a, c, b) = ints%ovvv(i, a, b, c)
                  end do
               end do
            end do
         end do
      end if
      call amps%init(no, nv)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  amps%t2(i, j, a, b) = 0.1_dp*lcg(seed)
               end do
            end do
         end do
      end do
      do a = 1, nv
         do i = 1, no
            amps%t1(i, a) = 0.1_dp*lcg(seed)
         end do
      end do
      if (symmetric) then
         ! t2(i,j,a,b) = t2(j,i,b,a)
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     if (j < i .or. (j == i .and. b < a)) amps%t2(i, j, a, b) = amps%t2(j, i, b, a)
                  end do
               end do
            end do
         end do
      end if
   end subroutine random_problem

   subroutine two_electrons(error)
      !! One occupied orbital: the antisymmetriser of the (T) vanishes, whatever the integrals and amplitudes
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_amps_t) :: amps
      type(error_t) :: status
      integer(int64) :: seed
      real(dp) :: e_t, worst
      integer :: trial

      worst = 0.0_dp
      do trial = 1, 3
         seed = 12345_int64 + int(trial, int64)
         call random_problem(1, 6 + trial, seed, .false., ints, amps)
         call rccsd_triples_reference(ints, amps, e_t, status)
         if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
         worst = max(worst, abs(e_t))
         call amps%release()
      end do
      call check(error, worst < 1.0e-12_dp, "two electrons: |E(T)| = "//to_string(worst, "es9.2")//" (random integrals)")
   end subroutine two_electrons

   subroutine random_symmetric(error)
      !! nv = 40 is above one c batch of 32 and a batch of 5 cuts every b: the routine against the independent
      !! evaluation, on integrals nothing about which is molecular; and the device version on the same inputs with
      !! tiles of 1 and 7 triples and the default budget
      type(error_type), allocatable, intent(out) :: error

      integer(int64), parameter :: RANDOM_TILES(3) = [1_int64, 7_int64, 0_int64]
      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_amps_t) :: amps
      type(trc_rccsd_triples_options_t) :: topts
      type(trc_rccsd_triples_dev_options_t) :: dopts
      type(error_t) :: status
      integer(int64) :: seed
      real(dp) :: e_32, e_5, e_host, e_dev(3), worst_dev, e_ijk(3), worst_ijk
      integer :: k

      seed = 987654321_int64
      call random_problem(3, 40, seed, .true., ints, amps)
      call map_ints(ints)
      call amps%push()
      call rccsd_triples_reference(ints, amps, e_32, status)
      if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      topts%c_batch = 5
      call rccsd_triples_reference(ints, amps, e_5, status, topts)
      call host_triples(ints, amps%t1, amps%t2, e_host)
      do k = 1, 3
         dopts%tile_triples = RANDOM_TILES(k)
         call rccsd_triples_dev_run(ints, amps, e_dev(k), status, dopts)
         if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      end do
      worst_dev = max(maxval(abs(e_dev - e_host)), maxval(abs(e_dev - e_32)))
      ! the occupied-driven routine: no = 3 has distinct, two-equal and three-equal occupied triples
      do k = 1, 3
         dopts%tile_triples = RANDOM_TILES(k)
         call rccsd_triples_ijk_run(ints, amps, e_ijk(k), status, dopts)
         if (status%has_error()) error stop "check_rccsd_triples: "//status%get_message()
      end do
      worst_ijk = max(maxval(abs(e_ijk - e_host)), maxval(abs(e_ijk - e_32)))
      call amps%release()
      call ints%release()
      call check(error, abs(e_32 - e_host) < TOL_HOST_REL*abs(e_host) .and. abs(e_5 - e_host) < TOL_HOST_REL*abs(e_host) .and. &
                 abs(e_host) > 1.0e-6_dp .and. worst_dev < TOL_DEV_REL*abs(e_host) .and. &
                 maxval(abs(e_dev - e_dev(3))) < TOL_DEV_REL*abs(e_host) .and. worst_ijk < TOL_DEV_REL*abs(e_host), &
                 "random symmetric no = 3, nv = 40: batch 32 "//to_string(e_32, "f20.14")//", batch 5 "// &
                 to_string(e_5, "f20.14")//", independent "//to_string(e_host, "f20.14")//", |difference| "// &
                 to_string(max(abs(e_32 - e_host), abs(e_5 - e_host)), "es9.2")//"; device tiles 1, 7, default "// &
                 to_string(e_dev(1), "f20.14")//" "//to_string(e_dev(2), "f20.14")//" "//to_string(e_dev(3), "f20.14")// &
                 ", largest difference from the independent evaluation and the reference "//to_string(worst_dev, "es9.2")// &
                 "; occupied-driven tiles 1, 7, default "//to_string(e_ijk(1), "f20.14")//" "//to_string(e_ijk(2), "f20.14")// &
                 " "//to_string(e_ijk(3), "f20.14")//", largest difference "//to_string(worst_ijk, "es9.2"))
   end subroutine random_symmetric

   subroutine map_ints(ints)
      !! Puts the blocks the (T) reads on the device, as `ints%build` leaves them; `ints%release` takes them off
      type(trc_rccsd_ints_t), intent(inout) :: ints

      !$acc enter data copyin(ints%ooov, ints%ovov, ints%ovvv, ints%eps_o, ints%eps_v)
   end subroutine map_ints

   subroutine unrolled_vs_table(error)
      !! The 36 terms of the device item routine, written out row by row, against the table-driven sum of the host
      !! routine (RCCSD_TRIPLES_WMAP) on random W, V inputs: distinct, two-equal and all-equal (a,b,c)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NO = 4, NV = 5
      integer, parameter :: ABC(3, 5) = reshape([4, 2, 1, 3, 3, 1, 5, 2, 2, 2, 2, 2, 1, 1, 1], [3, 5])
      real(dp), allocatable :: w(:, :, :, :), x(:, :, :), z(:, :, :, :)
      real(dp) :: wflat(6*NO**3)
      real(dp) :: ovov(NO, NV, NO, NV), t1(NO, NV), eps_o(NO), eps_v(NV)
      real(dp) :: e_table, e_unrolled, worst, d3, scale
      integer(int64) :: seed
      integer :: trial, p, i, j, k, a, b, c, perm(3, 6)

      seed = 20240607_int64
      allocate (w(NO, NO, NO, 6), x(NO, NO, NO), z(NO, NO, NO, 6))
      worst = 0.0_dp
      do trial = 1, size(ABC, 2)
         do p = 1, 6
            do k = 1, NO
               do j = 1, NO
                  do i = 1, NO
                     w(i, j, k, p) = lcg(seed)
                  end do
               end do
            end do
         end do
         do k = 1, NV
            do i = 1, NO
               t1(i, k) = lcg(seed)
            end do
         end do
         do b = 1, NV
            do j = 1, NO
               do a = 1, NV
                  do i = 1, NO
                     ovov(i, a, j, b) = lcg(seed)
                  end do
               end do
            end do
         end do
         do i = 1, NO
            eps_o(i) = -2.0_dp + 0.3_dp*real(i, dp) + 0.01_dp*lcg(seed)
         end do
         do a = 1, NV
            eps_v(a) = 1.0_dp + 0.4_dp*real(a, dp) + 0.01_dp*lcg(seed)
         end do
         a = ABC(1, trial)
         b = ABC(2, trial)
         c = ABC(3, trial)
         perm(:, 1) = [a, b, c]
         perm(:, 2) = [a, c, b]
         perm(:, 3) = [b, a, c]
         perm(:, 4) = [b, c, a]
         perm(:, 5) = [c, a, b]
         perm(:, 6) = [c, b, a]
         scale = 1.0_dp
         if (a == c) then
            scale = 6.0_dp
         else if (a == b .or. b == c) then
            scale = 2.0_dp
         end if
         ! the host routine's steps: X = W + V/2, Z = r3(X)/D, then the table-driven sum
         do p = 1, 6
            do k = 1, NO
               do j = 1, NO
                  do i = 1, NO
                     x(i, j, k) = w(i, j, k, p) + 0.5_dp*(ovov(i, perm(1, p), j, perm(2, p))*t1(k, perm(3, p)))
                  end do
               end do
            end do
            do k = 1, NO
               do j = 1, NO
                  do i = 1, NO
                     d3 = scale*(eps_o(i) + eps_o(j) + eps_o(k) - eps_v(a) - eps_v(b) - eps_v(c))
                     z(i, j, k, p) = (4.0_dp*x(i, j, k) + x(k, i, j) + x(j, k, i) - 2.0_dp*x(k, j, i) &
                                      - 2.0_dp*x(i, k, j) - 2.0_dp*x(j, i, k))/d3
                  end do
               end do
            end do
         end do
         e_table = rccsd_triples_energy_block(NO, RCCSD_TRIPLES_WMAP, w, z)
         wflat = reshape(w, [6*NO**3])
         e_unrolled = 0.0_dp
         do k = 1, NO
            do j = 1, NO
               do i = 1, NO
                  e_unrolled = e_unrolled + rccsd_triples_dev_point_energy(i, j, k, NO, NV, 0_int64, int(NO**3, int64), a, b, &
                                                                           c, wflat, 6_int64*NO**3, ovov, t1, &
                                                                           eps_o, eps_v)
               end do
            end do
         end do
         worst = max(worst, abs(e_table - e_unrolled)/max(abs(e_table), 1.0e-30_dp))
         if (abs(e_table) < 1.0e-8_dp) worst = max(worst, 1.0_dp)
      end do
      call check(error, worst < 1.0e-13_dp, "unrolled 36-term sum against the table-driven sum over 5 random triples "// &
                 "(distinct, two equal, three equal): worst relative difference "//to_string(worst, "es9.2"))
   end subroutine unrolled_vs_table

   subroutine empty_space(error)
      !! No occupied or no virtual orbitals is not a zero energy: it is refused, with the E(T) left at 0
      type(error_type), allocatable, intent(out) :: error

      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_amps_t) :: amps
      type(error_t) :: st_o, st_v
      integer(int64) :: seed
      real(dp) :: e_o, e_v
      logical :: ok

      seed = 4242_int64
      call random_problem(2, 4, seed, .false., ints, amps)
      ints%no = 0
      call rccsd_triples_reference(ints, amps, e_o, st_o)
      ints%no = 2
      ints%nv = 0
      call rccsd_triples_reference(ints, amps, e_v, st_v)
      ok = st_o%has_error() .and. st_v%has_error() .and. e_o == 0.0_dp .and. e_v == 0.0_dp
      ints%nv = 4
      call amps%release()
      call check(error, ok, "an empty occupied or virtual space is refused with E(T) = 0")
   end subroutine empty_space

end module check_rccsd_triples_suite

program check_rccsd_triples
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_triples_suite, only: collect_check_rccsd_triples
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_triples, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_triples
