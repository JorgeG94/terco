module check_h_suite
   !! Shells of angular momentum h (l = 5), against pyscf.
   !!
   !! Needs TERCO_LMAX >= 5: below it every case skips. h classes are the
   !! McMurchie-Davidson engine's (trc_mdg for exchange, trc_mdj for J), so
   !! these cases are the end-to-end statement that the generated MD kernels,
   !! the Boys function up to order 20, the Schwarz bounds, the one-electron
   !! kernels, the spherical transform and the XC grid evaluation agree with
   !! libcint's at h.
   !!
   !! Two kinds of case, each reading a reference test/pyscf_refs/generate.py
   !! wrote from pyscf (never run here):
   !!
   !!   ints_*   S, T, V and the Coulomb and exchange matrices of a fixed,
   !!            non-SCF density, in the Cartesian AO space, as row sums
   !!            against a fixed integer-hash weight matrix (weights below).
   !!            A matrix element that disagrees moves its row's sum, so the
   !!            row index says which shell is wrong without committing the
   !!            matrices themselves.
   !!   scf_*    a full SCF energy: RHF, UHF, and GGA, meta-GGA, hybrid and
   !!            range-separated hybrid DFT, spherical or Cartesian. DFT cases run on terco's own
   !!            grid, which they write to a probe for generate.py (the grid
   !!            is the one input the two codes would otherwise build
   !!            differently).
   !!
   !! The bases are pyscf's own cc-pV5Z data exported at full precision
   !! (scripts/pyscf_basis_to_json.py): the real set for water, and two
   !! uncontracted cut-downs that keep one h shell on C, N and O for the
   !! small systems (`min`: one primitive of each of s..h; `prims`: several per
   !! l). The cases on `min` take seconds to a minute on the host; the `prims`
   !! and cc-pV5Z ones take minutes and are labelled LONG. Test-drive cases run
   !! one at a time through `run_selected` (the first argument), which is how
   !! ctest and generate.py call them.
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, PCM_LMAX
   use trc_eri, only: trc_eri_t, trc_fock_options_t, eri_lmax_refusal, eri_grad_lmax_refusal
   use trc_tables, only: LMAX
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_pcm, only: trc_pcm_t, PCM_MODEL_CPCM
   use trc_chelpg, only: trc_chelpg_potential
   use trc_linalg, only: trc_linalg_t
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_scf_polarizability_rhf
   use trc_cpks, only: trc_scf_polarizability_rks
   use trc_tddft, only: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_reference_t, trc_tddft_solve_rks
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect
   use trc_sym_ao, only: sym_ao_build
   use trc_sym_scf, only: trc_sym_scf_t
   use trc_sph, only: sph_cart_to_sph
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_real, ref_read_vector, ref_read_int
   implicit none
   private
   public :: collect_h

   real(dp), parameter :: TOL_E = 1.0e-9_dp
   !! Eh, terco against pyscf for a full SCF (the g checks use the same)
   real(dp), parameter :: TOL_ROW = 1.0e-13_dp
   !! relative to the row's own sum of |weighted elements|, for the matrix rows.
   !! Measured worst row 6.1e-15 on the host (gfortran 15.1) and 5.7e-15 on a
   !! V100 (nvfortran 26.5), over S, T, V, J and K of ints_hco and the pure-
   !! harmonic overlap; the gate sits at about 16 times the larger.
   real(dp), parameter :: TOL_SYM = 1.0e-10_dp
   !! Eh, a symmetric run against the same molecule as C1 (the HF gate the
   !! point-group checks use, check_sym_c2v.F90's TOL_E_HF)
   real(dp), parameter :: THRESH = 1.0e-14_dp
   !! two-electron screening threshold of the matrix cases
   character(len=*), parameter :: B_MIN = "basis_sets/h-min-pyscf17.json"
   character(len=*), parameter :: B_PRIMS = "basis_sets/h-prims-pyscf17.json"
   character(len=*), parameter :: B_5Z = "basis_sets/cc-pv5z-pyscf17.json"

contains

   subroutine collect_h(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("refusals", test_refusals), &
                  new_unittest("ints_hco", test_ints_hco), &
                  new_unittest("sph_overlap_hco", test_sph_overlap_hco), &
                  new_unittest("rhf_water_min", test_rhf_water_min), &
                  new_unittest("rhf_water_prims", test_rhf_water_prims), &
                  new_unittest("rhf_n2_prims", test_rhf_n2_prims), &
                  new_unittest("uhf_hco_min", test_uhf_hco_min), &
                  new_unittest("rhf_water_min_cart", test_rhf_water_min_cart), &
                  new_unittest("rhf_water_cc5z", test_rhf_water_cc5z), &
                  new_unittest("lda_water_min", test_lda_water_min), &
                  new_unittest("pbe_n2_min", test_pbe_n2_min), &
                  new_unittest("pbe_water_prims", test_pbe_water_prims), &
                  new_unittest("b3lyp_hco_min", test_b3lyp_hco_min), &
                  new_unittest("r2scan_water_min", test_r2scan_water_min), &
                  new_unittest("wb97x_water_min", test_wb97x_water_min), &
                  new_unittest("rhf_water_min_sym", test_rhf_water_min_sym) &
                  ]
   end subroutine collect_h

   subroutine one_shell_basis(b, l)
      !! An s shell on each of two hydrogens and one shell of angular
      !! momentum `l` on the first (check_lmax_refusal.F90's basis).
      type(trc_basis_t), intent(out) :: b
      integer, intent(in) :: l
      real(dp) :: sh_e(1, 3), sh_c(1, 3), sh_r(3, 3), at_r(3, 2)
      at_r = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.4_dp], [3, 2])
      sh_e = reshape([1.0_dp, 1.0_dp, 0.8_dp], [1, 3])
      sh_c = 1.0_dp
      sh_r(:, 1) = at_r(:, 1)
      sh_r(:, 2) = at_r(:, 2)
      sh_r(:, 3) = at_r(:, 1)
      call b%build(3, [0, 0, l], [1, 1, 1], sh_e, sh_c, sh_r, 2, [1.0_dp, 1.0_dp], at_r, 1)
   end subroutine one_shell_basis

   subroutine test_refusals(error)
      !! What an h basis is told it cannot do: the four-centre integrals take
      !! it (at TERCO_LMAX=5), the analytic gradient names its own ceiling, an
      !! i shell is refused with the rebuild limit, and each module whose
      !! kernels stop below h answers with its own status and message: the PCM
      !! cavity and the CHELPG potential (PCM_LMAX = 4 < LMAX here), and the
      !! polarizability (CPHF, CPKS) and TDDFT, which need the dipole
      !! integrals (MULT_LMAX = 2). (The density-fitting and one-electron
      !! ceilings are the check_df_*_ceiling and check_one_e*_ceiling
      !! programs'.) Each refusal is asserted by its message, so a run that
      !! stops for another reason does not pass.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bh, bi
      if (LMAX < 5) then
         call skip_test(error, "SKIP: needs TERCO_LMAX >= 5 (h shells)")
         return
      end if
      call one_shell_basis(bh, 5)
      call check(error, eri_lmax_refusal(bh), "", "an h basis must pass the four-centre ceiling at TERCO_LMAX=5")
      if (.not. allocated(error)) then
         call check(error, index(eri_grad_lmax_refusal(bh), "analytic gradient stops at l = ") > 0, &
                    "the analytic gradient must refuse an h basis by name: "//eri_grad_lmax_refusal(bh))
      end if
      if (.not. allocated(error)) call refusal_pcm(error, bh)
      if (.not. allocated(error)) call refusal_chelpg(error, bh)
      if (.not. allocated(error)) call refusal_response(error, bh)
      call bh%release()
      if (allocated(error)) return
      call one_shell_basis(bi, 6)
      call check(error, index(eri_lmax_refusal(bi), "go up to h") > 0, &
                 "an i basis must be refused with the h ceiling: "//eri_lmax_refusal(bi))
      call bi%release()
   end subroutine test_refusals

   subroutine refusal_pcm(error, bh)
      !! The PCM build stops at PCM_LMAX with a validation error.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bh
      type(trc_pcm_t) :: pcm
      type(trc_linalg_t) :: la
      type(error_t) :: err
      call check(error, PCM_LMAX < 5, "this check needs PCM_LMAX below h (it is "//to_string(PCM_LMAX)//")")
      if (allocated(error)) return
      call pcm%init(PCM_MODEL_CPCM, 78.3553_dp, err)
      call check(error, .not. err%has_error(), "pcm init: "//err%get_message())
      if (.not. allocated(error)) then
         call pcm%build(la, bh, 302, 1.2_dp, err)
         call check(error, err%has_error(), "PCM must refuse an h basis")
      end if
      if (.not. allocated(error)) then
         call check(error, index(err%get_message(), "PCM: the basis has angular momentum above") > 0, &
                    "the PCM refusal names the surface-integral kernels: '"//err%get_message()//"'")
      end if
      call pcm%release()
      call la%release()
   end subroutine refusal_pcm

   subroutine refusal_chelpg(error, bh)
      !! The CHELPG potential stops at PCM_LMAX (its kernel is PCM's).
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bh
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), potential(:)
      allocate (density(bh%nao, bh%nao))
      density = 0.0_dp
      call trc_chelpg_potential(bh, density, reshape([0.0_dp, 0.0_dp, 3.0_dp], [3, 1]), potential, err)
      call check(error, err%has_error(), "CHELPG must refuse an h basis")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "above this build's potential-integral ceiling PCM_LMAX = "// &
                              to_string(PCM_LMAX)) > 0, &
                 "the CHELPG refusal names PCM_LMAX: '"//err%get_message()//"'")
   end subroutine refusal_chelpg

   subroutine refusal_response(error, bh)
      !! CPHF, CPKS and TDDFT come back unconverged with the multipole
      !! ceiling in their message, before any integral is built.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bh
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_linalg_t) :: la
      type(trc_cphf_options_t) :: copts
      type(trc_cphf_result_t) :: cres
      type(trc_tddft_options_t) :: topts
      type(trc_tddft_result_t) :: tres
      type(trc_tddft_reference_t) :: ref
      type(trc_xc_grid_t) :: xg
      type(trc_xc_functional_t) :: func
      real(dp) :: alpha(3, 3)
      real(dp), allocatable :: cmo(:, :), eps(:), dmat(:, :)
      allocate (cmo(bh%nao, bh%nao), eps(bh%nao), dmat(bh%nao, bh%nao))
      cmo = 0.0_dp
      eps = 0.0_dp
      dmat = 0.0_dp

      call trc_scf_polarizability_rhf(bh, pl, eri, 1, cmo, eps, alpha, copts, cres)
      call check(error, .not. cres%converged, "CPHF must refuse an h basis")
      if (.not. allocated(error)) then
         call check(error, index(cres%message, "multipole") > 0 .and. index(cres%message, "(l = 5)") > 0, &
                    "the CPHF refusal names the multipole ceiling and the shell: '"//trim(cres%message)//"'")
      end if
      if (allocated(error)) return

      call trc_scf_polarizability_rks(bh, pl, eri, xg, func, 1, cmo, eps, dmat, alpha, copts, cres)
      call check(error, .not. cres%converged, "CPKS must refuse an h basis")
      if (.not. allocated(error)) then
         call check(error, index(cres%message, "multipole") > 0 .and. index(cres%message, "(l = 5)") > 0, &
                    "the CPKS refusal names the multipole ceiling and the shell: '"//trim(cres%message)//"'")
      end if
      if (allocated(error)) return

      call trc_tddft_solve_rks(bh, pl, eri, la, ref, topts, tres)
      call check(error, .not. tres%converged .and. tres%err%has_error(), "TDDFT must refuse an h basis")
      if (.not. allocated(error)) then
         call check(error, index(tres%err%get_message(), "multipole") > 0 .and. &
                    index(tres%err%get_message(), "(l = 5)") > 0, &
                    "the TDDFT refusal names the multipole ceiling and the shell: '"//tres%err%get_message()//"'")
      end if
      call la%release()
   end subroutine refusal_response

   subroutine test_ints_hco(error)
      type(error_type), allocatable, intent(out) :: error
      call run_ints(error, "ints_hco", "hco.xyz", B_MIN)
   end subroutine test_ints_hco

   subroutine test_sph_overlap_hco(error)
      !! The spherical transform at h (and every l below it): the pure-harmonic
      !! overlap T^T S T of the h-min basis, unit-normalised, against pyscf's
      !! cart=False overlap. Unit normalisation takes terco's per-shell
      !! constant (its spherical S diagonal is not 1 above p) out of the
      !! comparison, leaving the transform's structure: orthogonality of the
      !! harmonics, their order within a shell (m = -l..l) and their signs.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(error_t) :: err
      integer :: natm, nao, nsph, refu, nref, i, j
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), scart(:, :), tmat(:, :), vmat(:, :), ssph(:, :), ref(:), dg(:)
      logical :: ok
      if (LMAX < 5) then
         call skip_test(error, "SKIP: needs TERCO_LMAX >= 5 (h shells)")
         return
      end if
      call read_xyz("hco.xyz", natm, zint, at_r)
      call trc_basis_from_json(B_MIN, natm, zint, at_r, bas, err, overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call check(error, bas%spherical, "the basis should come back spherical")
      if (allocated(error)) return
      nao = bas%nao
      nsph = bas%nao_sph
      allocate (scart(nao, nao), tmat(nao, nao), vmat(nao, nao), ssph(nsph, nsph), dg(nsph))
      call bas%to_device()
      call pl%build(bas, THRESH)
      call pl%to_device()
      !$acc enter data create(scart, tmat, vmat, ssph)
      call trc_1e(bas, pl, scart, tmat, vmat)
      call sph_cart_to_sph(bas, scart, ssph)
      !$acc update self(ssph)
      !$acc exit data delete(scart, tmat, vmat, ssph)
      call pl%release()
      call bas%release()
      do i = 1, nsph
         dg(i) = sqrt(ssph(i, i))
      end do
      do j = 1, nsph
         do i = 1, nsph
            ssph(i, j) = ssph(i, j)/(dg(i)*dg(j))
         end do
      end do
      call ref_open("h_sph_overlap_hco.ref", refu, ok)
      call check(error, ok, "missing pyscf_refs/h_sph_overlap_hco.ref -- run test/pyscf_refs/generate.py h_sph_overlap_hco")
      if (allocated(error)) return
      call ref_read_int(refu, "NAO", nref)
      call check(error, nref == nsph, "spherical AO count "//to_string(nsph)//" against the reference's "//to_string(nref))
      if (allocated(error)) then
         call ref_close(refu)
         return
      end if
      allocate (ref(nsph))
      call ref_read_vector(refu, "S_ROWS", nsph, ref)
      call ref_close(refu)
      call compare_rows(error, "S (spherical, unit-normalised)", ssph, ref)
   end subroutine test_sph_overlap_hco

   subroutine test_rhf_water_min(error)
      !! One primitive of each of s..h per heavy atom: every h class, small.
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "rhf_water_min", "water.xyz", B_MIN, "", 0, .false.)
   end subroutine test_rhf_water_min

   subroutine test_rhf_water_prims(error)
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "rhf_water_prims", "water.xyz", B_PRIMS, "", 0, .false.)
   end subroutine test_rhf_water_prims

   subroutine test_rhf_n2_prims(error)
      !! Two h centres: (hh|hh) classes between different atoms.
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "rhf_n2_prims", "tops_n2.xyz", B_PRIMS, "", 0, .false.)
   end subroutine test_rhf_n2_prims

   subroutine test_uhf_hco_min(error)
      !! HCO (2A'): open shell, unequal centres, no symmetry beyond the plane.
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "uhf_hco_min", "hco.xyz", B_MIN, "", 1, .false.)
   end subroutine test_uhf_hco_min

   subroutine test_rhf_water_min_cart(error)
      !! The Cartesian AO space (pyscf cart=True): 21 components per h shell.
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "rhf_water_min_cart", "water.xyz", B_MIN, "", 0, .true.)
   end subroutine test_rhf_water_min_cart

   subroutine test_rhf_water_cc5z(error)
      !! The real cc-pV5Z (contracted, 201 spherical AOs): the one case with
      !! contracted h-bearing basis data.
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "rhf_water_cc5z", "water.xyz", B_5Z, "", 0, .false.)
   end subroutine test_rhf_water_cc5z

   subroutine test_lda_water_min(error)
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "lda_water_min", "water.xyz", B_MIN, "svwn", 0, .false.)
   end subroutine test_lda_water_min

   subroutine test_pbe_n2_min(error)
      !! GGA on two h centres.
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "pbe_n2_min", "tops_n2.xyz", B_MIN, "pbe", 0, .false.)
   end subroutine test_pbe_n2_min

   subroutine test_pbe_water_prims(error)
      !! The same functional with several primitives per l (the slower check).
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "pbe_water_prims", "water.xyz", B_PRIMS, "pbe", 0, .false.)
   end subroutine test_pbe_water_prims

   subroutine test_b3lyp_hco_min(error)
      !! UKS hybrid on the open-shell radical.
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "b3lyp_hco_min", "hco.xyz", B_MIN, "b3lyp", 1, .false.)
   end subroutine test_b3lyp_hco_min

   subroutine test_r2scan_water_min(error)
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "r2scan_water_min", "water.xyz", B_MIN, "r2scan", 0, .false.)
   end subroutine test_r2scan_water_min

   subroutine test_wb97x_water_min(error)
      !! A range-separated hybrid: the erfc-attenuated exchange (omega > 0)
      !! through the MD exchange engine at h, with its own XC functional.
      type(error_type), allocatable, intent(out) :: error
      call run_scf(error, "wb97x_water_min", "water.xyz", B_MIN, "wb97x", 0, .false.)
   end subroutine test_wb97x_water_min

   subroutine test_rhf_water_min_sym(error)
      !! Point-group symmetry with h shells (the shell transform blocks and the
      !! AO permutations are per l, so h is one more l): water on the C2v
      !! standard frame, run with the irrep blocks and as C1, which must agree
      !! to TOL_SYM, and against pyscf's rhf_water_min energy (a rotation of
      !! the molecule does not change it).
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_scf_t) :: sc
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: rs, r1
      type(error_t) :: err
      integer :: natm, refu, nocc
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      real(dp) :: e_ref
      logical :: ok
      if (LMAX < 5) then
         call skip_test(error, "SKIP: needs TERCO_LMAX >= 5 (h shells)")
         return
      end if
      call read_xyz("water.xyz", natm, zint, at_r)
      geom = sym_geom_detect(at_r, zint, err, basis_id=zint)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, "C2v", "the group of water")
      if (allocated(error)) return
      call trc_basis_from_json(B_MIN, natm, zint, geom%r_std, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call check(error, maxval(bas%sh_l) == 5 .and. bas%spherical, "the basis should be spherical and reach h")
      if (allocated(error)) return
      call bas%to_device()
      nocc = sum(zint)/2
      opts%conv_energy = 1.0e-11_dp
      opts%conv_diis = 1.0e-8_dp
      opts%max_iteration = 200
      sc%ao = sym_ao_build(bas, geom, err)
      call check(error, .not. err%has_error(), "sym_ao_build: "//err%get_message())
      if (.not. allocated(error)) then
         call trc_scf_run(bas, nocc, nocc, opts, rs, sym=sc)
         call check(error, rs%converged, "the symmetric SCF did not converge: "//trim(rs%message))
      end if
      if (.not. allocated(error)) then
         call sc%release()
         call trc_scf_run(bas, nocc, nocc, opts, r1)
         call check(error, r1%converged, "the C1 SCF did not converge: "//trim(r1%message))
      end if
      call bas%release()
      if (allocated(error)) return
      print "(4x,a,es9.2)", "rhf_water_min_sym: |E(C2v) - E(C1)| ", abs(rs%energy - r1%energy)
      call check(error, abs(rs%energy - r1%energy) <= TOL_SYM, "|E(C2v) - E(C1)| = "// &
                 to_string(abs(rs%energy - r1%energy), "es10.2")//" exceeds "//to_string(TOL_SYM, "es10.2"))
      if (allocated(error)) return
      call ref_open("h_rhf_water_min.ref", refu, ok)
      call check(error, ok, "missing pyscf_refs/h_rhf_water_min.ref")
      if (allocated(error)) return
      call ref_read_real(refu, "E", e_ref)
      call ref_close(refu)
      call check(error, abs(rs%energy - e_ref) <= TOL_E, "C2v energy "//to_string(rs%energy, "es20.12")// &
                 " against pyscf's "//to_string(e_ref, "es20.12"))
   end subroutine test_rhf_water_min_sym

   pure function weight(i, j) result(w)
      !! The row-sum weights, exact in both languages: integer arithmetic,
      !! one division. generate.py's `row_weights` is the same formula.
      integer, intent(in) :: i, j
      real(dp) :: w
      w = real(modulo(i*7919 + j*104729 + i*j*31, 1009), dp)/1009.0_dp - 0.5_dp
   end function weight

   subroutine row_sums(m, rs, ra)
      !! rs(i) = sum_j w(i,j) m(i,j); ra(i) = sum_j |w(i,j) m(i,j)|.
      real(dp), intent(in) :: m(:, :)
      real(dp), intent(out) :: rs(size(m, 1)), ra(size(m, 1))
      integer :: i, j
      real(dp) :: x
      rs = 0.0_dp
      ra = 0.0_dp
      do j = 1, size(m, 2)
         do i = 1, size(m, 1)
            x = weight(i, j)*m(i, j)
            rs(i) = rs(i) + x
            ra(i) = ra(i) + abs(x)
         end do
      end do
   end subroutine row_sums

   subroutine synth_density(n, dd)
      !! The same deterministic, dense, symmetric density as
      !! check_rys_fock_vs_pyscf.F90 and check_erf_fock_vs_pyscf.F90.
      integer, intent(in) :: n
      real(dp), intent(out) :: dd(n, n)
      integer :: ii, jj
      do jj = 1, n
         do ii = 1, n
            dd(ii, jj) = 0.4_dp + 0.1_dp*cos(real(3*ii + 5*jj, dp))
         end do
      end do
      dd = 0.5_dp*(dd + transpose(dd))
   end subroutine synth_density

   subroutine compare_rows(error, label, m, refrs)
      !! The worst row of `m` against its reference row sum, relative to the
      !! row's own absolute sum.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: m(:, :), refrs(:)
      real(dp), allocatable :: rs(:), ra(:)
      real(dp) :: worst, rel
      integer :: i, iw
      allocate (rs(size(m, 1)), ra(size(m, 1)))
      call row_sums(m, rs, ra)
      worst = 0.0_dp
      iw = 0
      do i = 1, size(m, 1)
         rel = abs(rs(i) - refrs(i))/(1.0_dp + ra(i))
         if (rel > worst) then
            worst = rel
            iw = i
         end if
      end do
      print "(4x,a,a,es10.2,a,i0)", label, ": worst relative row-sum difference ", worst, " at row ", iw
      call check(error, worst <= TOL_ROW, label//": row "//to_string(iw)//" differs from pyscf by "// &
                 to_string(worst, "es10.3")//" (relative), above "//to_string(TOL_ROW, "es10.3"))
   end subroutine compare_rows

   subroutine run_ints(error, name, xyz, basis)
      !! S, T, V, J and K on the fixed density, Cartesian, against pyscf.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name, xyz, basis
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(trc_pairlist_t) :: pl
      type(error_t) :: err
      integer :: natm, nao, unit, refu, nao_ref
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), d(:, :), gj(:, :), gk(:, :), smat(:, :), tmat(:, :), vmat(:, :)
      real(dp), allocatable :: ref(:)
      logical :: ok

      if (LMAX < 5) then
         call skip_test(error, "SKIP: needs TERCO_LMAX >= 5 (h shells)")
         return
      end if
      call read_xyz(xyz, natm, zint, at_r)
      ! force_cart: the Cartesian AO space pyscf's cart=True has. prune=.false.:
      ! the reference is built from the file's literal columns.
      call trc_basis_from_json(basis, natm, zint, at_r, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call check(error, maxval(bas%sh_l) == 5, "the basis should reach h")
      if (allocated(error)) return
      nao = bas%nao
      allocate (d(nao, nao), gj(nao, nao), gk(nao, nao), smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      call synth_density(nao, d)

      call bas%to_device()
      call eri%build(bas, THRESH)
      ! K = -2 g from a second call with j_scale = 0, k_scale = 1 (the kernels are
      ! linear in the scales; check_erf_fock_vs_pyscf.F90 does the same).
      call eri%fock(bas, d, gj, opts=trc_fock_options_t(density_screen=.false., j_scale=1.0_dp, k_scale=0.0_dp))
      call eri%fock(bas, d, gk, opts=trc_fock_options_t(density_screen=.false., j_scale=0.0_dp, k_scale=1.0_dp))
      call eri%release()
      gk = -2.0_dp*gk
      call pl%build(bas, THRESH)
      call pl%to_device()
      !$acc enter data create(smat, tmat, vmat)
      call trc_1e(bas, pl, smat, tmat, vmat)
      !$acc update self(smat, tmat, vmat)
      !$acc exit data delete(smat, tmat, vmat)
      call pl%release()
      call bas%release()

      ! The probe generate.py reads: geometry, density, and terco's matrices
      ! (S for the per-AO normalisation rescale; the rest for its printout).
      open (newunit=unit, file="h_probe_"//name//".bin", access="stream", form="unformatted", status="replace")
      write (unit) int(natm, int64), int(nao, int64)
      write (unit) int(zint, int64), at_r
      write (unit) d
      write (unit) smat
      write (unit) tmat
      write (unit) vmat
      write (unit) gj
      write (unit) gk
      close (unit)

      call ref_open("h_"//name//".ref", refu, ok)
      call check(error, ok, "missing pyscf_refs/h_"//name//".ref -- run test/pyscf_refs/generate.py h_"//name)
      if (allocated(error)) return
      call ref_read_int(refu, "NAO", nao_ref)
      call check(error, nao_ref == nao, "AO count "//to_string(nao)//" against the reference's "//to_string(nao_ref))
      if (allocated(error)) then
         call ref_close(refu)
         return
      end if
      allocate (ref(nao))
      call ref_read_vector(refu, "S_ROWS", nao, ref)
      call compare_rows(error, "S", smat, ref)
      if (.not. allocated(error)) then
         call ref_read_vector(refu, "T_ROWS", nao, ref)
         call compare_rows(error, "T", tmat, ref)
      end if
      if (.not. allocated(error)) then
         call ref_read_vector(refu, "V_ROWS", nao, ref)
         call compare_rows(error, "V", vmat, ref)
      end if
      if (.not. allocated(error)) then
         call ref_read_vector(refu, "J_ROWS", nao, ref)
         call compare_rows(error, "J", gj, ref)
      end if
      if (.not. allocated(error)) then
         call ref_read_vector(refu, "K_ROWS", nao, ref)
         call compare_rows(error, "K", gk, ref)
      end if
      call ref_close(refu)
   end subroutine run_ints

   subroutine run_scf(error, name, xyz, basis, functional, nunpaired, cartesian)
      !! A converged SCF energy against pyscf's, on the same exponents,
      !! contractions, geometry (Bohr, terco's constant) and DFT grid.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name, xyz, basis, functional
      integer, intent(in) :: nunpaired
      logical, intent(in) :: cartesian
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer :: natm, nelec, na, nb, unit, refu
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      real(dp) :: e_ref, de
      logical :: ok
      character(len=200) :: msg
      character(len=8) :: fn8

      if (LMAX < 5) then
         call skip_test(error, "SKIP: needs TERCO_LMAX >= 5 (h shells)")
         return
      end if
      call read_xyz(xyz, natm, zint, at_r)
      nelec = sum(zint)
      nb = (nelec - nunpaired)/2
      na = nelec - nb
      call trc_basis_from_json(basis, natm, zint, at_r, bas, err, force_cart=cartesian, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call check(error, maxval(bas%sh_l) == 5, "the basis should reach h")
      if (allocated(error)) return
      call check(error, bas%spherical .neqv. cartesian, "AO space is not the requested one")
      if (allocated(error)) return
      call bas%to_device()

      opts = trc_scf_options_t()
      opts%functional = functional
      opts%grid_level = 3
      opts%conv_energy = 1.0e-11_dp
      opts%conv_diis = 1.0e-8_dp
      opts%max_iteration = 200
      opts%force_cartesian = cartesian
      call trc_scf_run(bas, na, nb, opts, res)
      call bas%release()
      call check(error, res%converged, name//": terco SCF did not converge: "//trim(res%message))
      if (allocated(error)) return

      ! The probe generate.py reads: molecule, spin, AO space, functional,
      ! the grid (DFT) and terco's energy (for its own printout only).
      fn8 = functional
      open (newunit=unit, file="h_probe_"//name//".bin", access="stream", form="unformatted", status="replace")
      write (unit) int(natm, int64), int(res%nspin, int64), int(na, int64), int(nb, int64)
      write (unit) int(zint, int64), at_r
      write (unit) fn8
      if (allocated(res%grid_w)) then
         write (unit) int(size(res%grid_w), int64)
         write (unit) res%grid_r, res%grid_w
      else
         write (unit) 0_int64
      end if
      write (unit) res%energy
      close (unit)

      call ref_open("h_"//name//".ref", refu, ok)
      call check(error, ok, "missing pyscf_refs/h_"//name//".ref -- run test/pyscf_refs/generate.py h_"//name)
      if (allocated(error)) return
      call ref_read_real(refu, "E", e_ref)
      call ref_close(refu)
      de = abs(res%energy - e_ref)
      write (msg, "(a,a,f20.12,a,f20.12,a,es9.2)") name, ": E(terco) ", res%energy, "  E(pyscf) ", e_ref, &
         "  |dE| ", de
      print "(4x,a)", trim(msg)
      call check(error, de <= TOL_E, trim(msg))
   end subroutine run_scf

end module check_h_suite

program check_h_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, run_selected
   use check_h_suite, only: collect_h
   implicit none
   integer :: stat
   character(len=64) :: name
   stat = 0
   ! One case per ctest entry (and per generate.py case): the name is the
   ! first argument; none runs them all.
   if (command_argument_count() >= 1) then
      call get_command_argument(1, name)
      call run_selected(collect_h, trim(name), error_unit, stat)
   else
      call run_testsuite(collect_h, error_unit, stat, parallel=.false.)
   end if
   if (stat > 0) then
      write (error_unit, '(a)') "check_h_vs_pyscf: failed"
      error stop 1
   end if
end program check_h_vs_pyscf
