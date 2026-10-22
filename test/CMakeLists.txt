#
# New checks are test-drive suites (a module of named unit tests, `check`
# with a message, run_testsuite); see AGENTS.md. The older checks are plain
# programs and move over when they are next touched.
#
find_package("test-drive" REQUIRED)

#
# The checks.
#
# Every one compares against something outside terco -- libfint through a
# recorded oracle, pyscf, or another of terco's own paths that must agree by
# construction. None of them assert a number this library produced and then
# froze, because that tests only that nothing changed.
#
# They run on the GPU, so they cannot run in CI on a machine without one.
# `CUDA_VISIBLE_DEVICES=0` pins them to one device: terco is single-GPU today,
# and a four-GPU machine otherwise reports timings that are not what they look
# like.
#

#
# The test programs supply their own basis. terco ships no basis set -- a
# consumer passes exponents and coefficients -- so the tests have to build one,
# and they may as well share the code.
#
set(TERCO_TEST_DEFS TRC_LMAX=${TERCO_LMAX} TRC_FOCK6)

add_library(terco_test_basis STATIC trc_test_basis.F90)
target_compile_definitions(terco_test_basis PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(terco_test_basis PUBLIC ${main_lib})
set_target_properties(terco_test_basis PROPERTIES Fortran_MODULE_DIRECTORY
                                                  ${CMAKE_BINARY_DIR}/modules)

#
# The reader for the committed pyscf references under test/pyscf_refs/ (see
# that directory's README.md and pyscf_ref_io.F90's own header). Every check
# that used to hand a dump to a `*_ref.py` companion links this instead, and
# compares in Fortran against a file test/pyscf_refs/generate.py wrote
# offline -- pyscf itself never runs in CI.
#
add_library(terco_pyscf_ref_io STATIC pyscf_ref_io.F90 golden_io.F90)
target_compile_definitions(terco_pyscf_ref_io PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(terco_pyscf_ref_io PUBLIC ${main_lib})
set_target_properties(terco_pyscf_ref_io PROPERTIES Fortran_MODULE_DIRECTORY
                                                     ${CMAKE_BINARY_DIR}/modules)

#
# The library's definitions are PRIVATE -- TRC_FOCK6 is an
# implementation choice and no consumer should need them. A test can, though:
# `check_binfock` compares against a raw oracle and has to know whether the
# digestion folded the eight-fold symmetry into six updates or enumerated it,
# because the two return different things to combine.
#
# Passing them here rather than making them PUBLIC keeps that knowledge where
# it belongs. Getting this wrong is not a build failure -- the test compiles,
# runs, takes the unfolded branch and reports a worst error of 5e+01 as though
# the kernels were broken.
#
function(terco_add_check name source)
  cmake_parse_arguments(ARG "NEEDS_BASIS;NEEDS_PYSCF_REF" "" "ARGS" ${ARGN})
  add_executable(${name} ${source})
  target_compile_definitions(${name} PRIVATE ${TERCO_TEST_DEFS})
  target_link_libraries(${name} PRIVATE ${main_lib})
  if(ARG_NEEDS_BASIS)
    target_link_libraries(${name} PRIVATE terco_test_basis)
  endif()
  if(ARG_NEEDS_PYSCF_REF)
    target_link_libraries(${name} PRIVATE terco_pyscf_ref_io)
  endif()
  target_include_directories(${name} PRIVATE ${CMAKE_BINARY_DIR}/modules)
  add_test(NAME terco/${name} COMMAND ${name} ${ARG_ARGS})
  # The oracles (reference.dat/reference1e.dat/reference_df.dat/
  # reference_xc.dat/pyscf_refs/*.ref) are read by relative path, so the
  # tests run from test/. TERCO_CART_ORDER is not needed by the Fortran
  # side (it already knows its own TRC_CART_ORDER at compile time), but is
  # set here anyway for uniformity with test/pyscf_refs/generate.py, which
  # DOES need it when regenerating a reference by hand outside this build.
  set_tests_properties(
    terco/${name} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                             ENVIRONMENT "CUDA_VISIBLE_DEVICES=0;TERCO_CART_ORDER=${TERCO_CART_ORDER}")
endfunction()

#
# Correctness with NOTHING external: closed forms and identities.
#
# Runs first because it is the one that can run anywhere. The oracle-backed
# checks below are stronger evidence -- an independent implementation catches
# a shared misunderstanding that an identity cannot -- but they need a 6 MB
# recorded file, and pyscf, and someone with neither should still be able to
# tell whether this library works.
#
terco_add_check(check_selftest check_selftest.F90)

#
# trc_linalg's shape-inferred forms (gemm_mat, dot_2d/3d, syev_mat, sandwich,
# density) against the raw forms they must stay bitwise identical to -- see
# check_linalg_mat.F90's own header. No basis, no pyscf: it is dense linear
# algebra, not integrals.
#
terco_add_check(check_linalg_mat check_linalg_mat.F90)

# The step tables' SI-prefixed counts (2.45 PFLOP, not 2453056.2 GFLOP).
terco_add_check(check_si check_si.F90)
target_link_libraries(check_si PRIVATE test-drive::test-drive)

# A basis above this build's four-centre ceiling is refused with a message,
# not run into NaN integrals and an eigensolver failure.
terco_add_check(check_lmax_refusal check_lmax_refusal.F90)
target_link_libraries(check_lmax_refusal PRIVATE test-drive::test-drive terco_c_interfaces)

# trc_constants' derived relations (ANGSTROM_TO_BOHR vs BOHR_TO_ANGSTROM,
# KB_HARTREE, HARTREE_TO_KCALMOL, HARTREE_TO_CM1, all against the OTHER
# CODATA 2018 values on the same file) and the atomic-mass tables, H/C/O/P
# and the out-of-range refusal. No basis, no pyscf: these are closed-form
# identities and a lookup table.
terco_add_check(check_constants check_constants.F90)
target_link_libraries(check_constants PRIVATE test-drive::test-drive)

#
# The staleness guard over every committed pyscf reference under
# test/pyscf_refs/ (see that directory's README.md and check_staleness.py's
# own header). Standard-library Python only -- no pyscf, no numpy -- so it
# runs in CI like every other check here, and fails loudly (not silently)
# when a reference's own recorded inputs (geometry, basis JSON, the
# Fortran source that builds a fixed density or grid, or generate.py
# itself) have changed since it was written.
#
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  add_test(NAME terco/check_pyscf_refs_fresh
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pyscf_refs/check_staleness.py)
  set_tests_properties(
    terco/check_pyscf_refs_fresh
    PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  # refhash.py's own tests: an edit to one generator restales only its
  # references, a shared helper exactly its users', a comment none.
  add_test(NAME terco/check_refhash
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pyscf_refs/test_refhash.py)
endif()

#
# The translated XC kernels against libxc, through the values
# test/pyscf_refs/generate.py's gen_xc_kernels() records in reference_xc.dat.
# Every derivative to second order, both spin cases, eleven decades of
# density. Needs no GPU and no pyscf to run; pyscf only to regenerate it.
#
terco_add_check(check_xc check_xc.F90)

# The quadrature against closed forms: Lebedev moments, the partition of
# unity, an exact radial Gaussian, and sums of Gaussians and Slaters over
# water integrating to the atom count.
terco_add_check(check_grid check_grid.F90 NEEDS_BASIS)

# becke_partition_derivatives' O(n_atoms^2) production routine against a
# hand-kept copy of the O(n_atoms^3) implementation it replaced, over
# every point of a real water and formamide grid, both partition/adjust
# combinations trc_dft_partition.F90 offers.
terco_add_check(check_dft_partition_deriv check_dft_partition_deriv.F90 NEEDS_BASIS)

#
# Basis functions on grid points: the numerical overlap against trc_1e's
# analytic one, and the batched, screened evaluation against the unscreened
# one. (check_ao itself -- the pyscf comparison on top of it -- is
# registered further down, alongside check_mult and check_xc_energy's.)
#

# Four-centre integrals and the Fock digestion, against libfint's own numbers
# recorded in reference.dat.
terco_add_check(check_binfock check_binfock.F90)
# The same check with every bin split into chunks of at most 3 pairs
# (TRC_BIN_CAP; trc_bins.F90's split_large_bins). A real system only splits
# past 32768 pairs in a bin, around 14,000 basis functions; this is that
# path at a size CI can run, against the same oracle and tolerance.
add_test(NAME terco/check_binfock_split COMMAND check_binfock)
set_tests_properties(terco/check_binfock_split PROPERTIES
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  ENVIRONMENT "CUDA_VISIBLE_DEVICES=0;TERCO_CART_ORDER=${TERCO_CART_ORDER};TRC_BIN_CAP=3")
# The decontraction transform against an external answer: C^T S_p C = S and
# C^T T_p C = T, from trc_1e, which shares nothing with the four-centre path.
terco_add_check(check_decontract check_decontract.F90 NEEDS_BASIS ARGS water.xyz basis_sets/cc-pvdz.json)
# Second row: twelve silicon s primitives under four contraction columns, the
# case the whole exercise exists for.
add_test(NAME terco/check_decontract_si COMMAND check_decontract sih4.xyz basis_sets/cc-pvdz.json)
set_tests_properties(terco/check_decontract_si PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                   ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")
# build_631g against trc_basis_from_json, on water through a five-water
# cluster -- the one check allowed to call build_631g for its own basis
# (see that subroutine's comment in trc_test_basis.F90). Water through
# five waters, not one geometry: a per-shell coefficient scale error can
# still agree with the reference to several digits at a single small
# geometry, which is what let the common_fac_sp regression this check
# guards against go unnoticed before.
terco_add_check(check_basis_631g_ref check_basis_631g_ref.F90 NEEDS_BASIS)
# The spherical (pure) AO space against pyscf's cart=False: water in
# def2-SVP and cc-pVDZ, a three-water cluster in def2-SVP, SiH4 in
# cc-pVDZ, and the OH radical open-shell in cc-pVDZ. See check_spherical.F90's
# own header for the reference numbers' provenance and def2-svp.json's own
# (copied from gmshpc's data/basis, itself the Basis Set Exchange's JSON
# bundle -- same source and format as every other file in basis_sets/).
terco_add_check(check_spherical check_spherical.F90 NEEDS_BASIS)
terco_add_check(sad_atoms sad_atoms.F90 NEEDS_BASIS)
terco_add_check(sadq_atoms sadq_atoms.F90 NEEDS_BASIS)
# The spherical averaging itself (trc_sad.F90's `average_atom_density`):
# determinism, block proportionality and trace preservation on a free
# O/cc-pVTZ atom's s/p/d/f shells, plus Tr(D S) = Z through the full
# trc_sad_build pipeline when this build's one-/two-electron ceiling
# reaches f (see sad_sphavg.F90's own header).
terco_add_check(sad_sphavg sad_sphavg.F90 NEEDS_BASIS)
terco_add_check(sap_guess sap_guess.F90 NEEDS_BASIS)
terco_add_check(sapcmp sapcmp.F90 NEEDS_BASIS)

# The public API end to end: overlap, kinetic, nuclear attraction, the two- and
# three-centre Coulomb integrals, and the multipoles, all against libfint.
terco_add_check(check_api check_api.F90)

# The one-electron kernels on their own, block by block, against libfint.
#
# There is deliberately no separate density-fitting check here: `check_df.F90`
# called the two- and three-centre dispatchers directly and went stale when
# they started taking pre-computed primitive-pair data instead of raw shells.
# `check_api` compares the same integrals against the same oracle through the
# public API, which is the only way a consumer can reach them anyway.
terco_add_check(check_1e check_1e.F90)

#
# The folded and enumerated digestions, on a symmetric density.
#
# They must agree for ANY symmetric density -- that is the premise of having
# two of them -- so this is a control rather than a comparison. It exists
# because they did NOT agree at LMAX=2 for a month, and nothing was looking:
# only the non-symmetric build reaches one_bin_item_nosym, so every other
# check passes without executing a line of it.
#
terco_add_check(check_nosym check_nosym.F90 NEEDS_BASIS ARGS water.xyz)
# 6-31G needs only l <= 1, so this only self-skips (SKIP + stop 0, same
# check_nosym.F90 branch check_nosym_f's comment below describes) at
# TERCO_LMAX=0 -- a build that drops even p functions. Covered for the
# same reason check_nosym_f is, below.
set_tests_properties(terco/check_nosym PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")

#
# The same comparison with f functions, which is the only test that exercises
# the hybrid dispatch: at LMAX=3 with the unrolled set generated at 2, `fock`
# splits the work between the two paths and `fock_nosym` does not, so they
# agree only if the split is exact.
#
# Registered at every ceiling and skipped from inside when the build cannot
# reach l = 3, so it never silently stops testing anything: a build that CAN
# run it always does.
#
terco_add_check(check_nosym_f check_nosym.F90 NEEDS_BASIS ARGS water.xyz 3)
# At TERCO_LMAX < 3 this prints "SKIP: ..." and `stop 0` instead of running
# anything (see check_nosym.F90's own header) -- a pass by exit code, with
# nothing checked. SKIP_REGULAR_EXPRESSION makes ctest report that as
# "Skipped" rather than "Passed" so a CI matrix that never raises TERCO_LMAX
# cannot mistake silence for coverage.
set_tests_properties(terco/check_nosym_f PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")

#
# Rys against HGP, same density and screening, on water/6-31G* and
# water/cc-pVDZ. GPU 2, not the fleet default of
# 0 `terco_add_check` otherwise sets, to stay off the GPUs other work on
# this host uses.
#
terco_add_check(check_rys_vs_hgp_631gs check_rys_vs_hgp.F90 NEEDS_BASIS
                ARGS basis_sets/6-31gs.json)
terco_add_check(check_rys_vs_hgp_ccpvdz check_rys_vs_hgp.F90 NEEDS_BASIS
                ARGS basis_sets/cc-pvdz.json)
set_tests_properties(
  terco/check_rys_vs_hgp_631gs terco/check_rys_vs_hgp_ccpvdz
  PROPERTIES ENVIRONMENT "CUDA_VISIBLE_DEVICES=0;TERCO_CART_ORDER=${TERCO_CART_ORDER}")

#
# eri%gradient (Rys, the only scheme it has) against a central difference
# of eri%fock's energy through HGP, RHF and UHF, water/6-31G* and
# water/cc-pVDZ. TRC_RYS_CLASSES=none, deliberately -- the default routes
# some d classes to the Rys energy kernels, and the finite-difference
# reference must go through HGP, so this cross-validates the gradient's
# algebra against an UNRELATED energy path, not Rys against itself.
#
terco_add_check(check_rys_grad_fd_631gs check_rys_grad_fd.F90 NEEDS_BASIS
                ARGS basis_sets/6-31gs.json)
terco_add_check(check_rys_grad_fd_ccpvdz check_rys_grad_fd.F90 NEEDS_BASIS
                ARGS basis_sets/cc-pvdz.json)
set_tests_properties(
  terco/check_rys_grad_fd_631gs terco/check_rys_grad_fd_ccpvdz
  PROPERTIES ENVIRONMENT "CUDA_VISIBLE_DEVICES=0;TERCO_CART_ORDER=${TERCO_CART_ORDER};TRC_RYS_CLASSES=none")

#
# eri%gradient's long-range (erf(omega*r12)/r12) exchange term alone
# against a central difference of E_K^lr(D) -- the omega > 0 analogue of
# check_rys_grad_fd, RHF only (see check_lr_exchange_grad_fd.F90's own
# header for why). water/6-31G* and water/cc-pVDZ, omega = 0.3 and 0.4
# (wB97X-V/wB97X and wB97's own values, src/xc/trc_xc_functional.F90) in
# one binary, plus a bit-identity check at omega = 0.
#
terco_add_check(check_lr_exchange_grad_fd_631gs check_lr_exchange_grad_fd.F90 NEEDS_BASIS
                ARGS basis_sets/6-31gs.json)
terco_add_check(check_lr_exchange_grad_fd_ccpvdz check_lr_exchange_grad_fd.F90 NEEDS_BASIS
                ARGS basis_sets/cc-pvdz.json)
set_tests_properties(
  terco/check_lr_exchange_grad_fd_631gs terco/check_lr_exchange_grad_fd_ccpvdz
  PROPERTIES ENVIRONMENT "CUDA_VISIBLE_DEVICES=0;TERCO_CART_ORDER=${TERCO_CART_ORDER}")

# The FULL analytic gradient (eri%gradient + trc_1e_gradient + nuclear
# repulsion) against finite differences of terco's own fully reconverged
# SCF energy -- RHF water, UHF OH -- see check_full_grad_fd.F90's header.
terco_add_check(check_full_grad_fd check_full_grad_fd.F90 NEEDS_BASIS)

# The same full analytic gradient, on SPHERICAL cc-pVDZ bases (water/RHF,
# formamide/RHF, OH/UHF), against pyscf's own Gradients class at a
# converged density -- the one check that exercises trc_scf_gradient's own
# spherical fold (see check_full_grad_vs_pyscf.F90's own header).
terco_add_check(check_full_grad_vs_pyscf check_full_grad_vs_pyscf.F90 NEEDS_BASIS NEEDS_PYSCF_REF)

# The FULL KS-DFT gradient, SAME recipe as check_xc_grad_vs_pyscf.F90
# (terco's own converged orbitals and grid, no pyscf SCF) but with NO
# force_cart: cc-pVDZ comes back SPHERICAL here, exercising
# trc_scf_gradient's fold TOGETHER with the XC term (see
# check_dft_grad_vs_pyscf.F90's own header) -- water/formamide RKS
# PBE/B3LYP on spherical cc-pVDZ and Cartesian 6-31G*, and UKS PBE/B3LYP
# NH2 from a fixed (core-guess) initial density -- B3LYP NH2 is the
# GLOBAL-hybrid UKS case (nonzero func%exx, omega/exx_lr still zero).
foreach(case water_pbe_ccpvdz_sph water_b3lyp_ccpvdz_sph
             formamide_pbe_ccpvdz_sph formamide_b3lyp_ccpvdz_sph
             formamide_pbe_631gs formamide_b3lyp_631gs
             nh2_pbe_ccpvdz_sph nh2_b3lyp_ccpvdz_sph
             water_wb97x_ccpvdz_sph nh2_wb97x_ccpvdz_sph
             water_wb97xv_ccpvdz_sph nh2_wb97xv_ccpvdz_sph)
  if(case STREQUAL "water_pbe_ccpvdz_sph")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func pbe)
  elseif(case STREQUAL "water_b3lyp_ccpvdz_sph")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func b3lyp)
  elseif(case STREQUAL "formamide_pbe_ccpvdz_sph")
    set(xyz formamide.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func pbe)
  elseif(case STREQUAL "formamide_b3lyp_ccpvdz_sph")
    set(xyz formamide.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func b3lyp)
  elseif(case STREQUAL "formamide_pbe_631gs")
    set(xyz formamide.xyz)
    set(json basis_sets/6-31gs.json)
    set(func pbe)
  elseif(case STREQUAL "formamide_b3lyp_631gs")
    set(xyz formamide.xyz)
    set(json basis_sets/6-31gs.json)
    set(func b3lyp)
  elseif(case STREQUAL "nh2_b3lyp_ccpvdz_sph")
    set(xyz nh2.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func b3lyp)
  elseif(case STREQUAL "water_wb97x_ccpvdz_sph")
    # Range-separated hybrid, spherical, RKS -- the full assembled
    # wB97X gradient (semilocal GGA term + trc_scf_grad.F90's own
    # long-range exchange pass) through trc_scf_gradient, SPHERICAL
    # fold included, guarded against pyscf's own analytic gradient.
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func wb97x)
  elseif(case STREQUAL "nh2_wb97x_ccpvdz_sph")
    # Same, UKS (ndens=2): the erf long-range exchange gradient's
    # unrestricted (per-spin-channel) pass, guarded for real rather
    # than only the RHF path check_lr_exchange_grad_vs_pyscf.F90 covers.
    set(xyz nh2.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func wb97x)
  elseif(case STREQUAL "water_wb97xv_ccpvdz_sph")
    # VV10 nonlocal correlation (wB97X-V), RKS, spherical -- the semilocal
    # GGA term, the erf long-range exchange gradient, AND
    # trc_vv10_grad%vv10_gradient_add's nonlocal term, all three assembled
    # by trc_scf_gradient. grid_response left at trc_scf_gradient's own
    # default (OFF, this program's own `.false.` below) -- pyscf's own
    # nlc gradient has no grid-response term to match against, see
    # trc_vv10_grad.F90's own header.
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func wb97x-v)
  elseif(case STREQUAL "nh2_wb97xv_ccpvdz_sph")
    # Same, UKS (ndens=2, fixed core-guess density) -- VV10 reads the
    # TOTAL density only (trc_vv10_grad.F90's own header), so this is the
    # one case that checks the UKS reduction (dtot = dens_a + dens_b)
    # for real.
    set(xyz nh2.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func wb97x-v)
  else()
    set(xyz nh2.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func pbe)
  endif()
  set(glev 3)
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/dft_grad_probe_${case}.bin)
  terco_add_check(check_dft_grad_vs_pyscf_${case} check_dft_grad_vs_pyscf.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${xyz} ${json} ${func} ${glev} ${dump} ${case})
endforeach()

# The CLI's --grad against a direct trc_run_gradient/trc_gradient call on
# the same molecule/basis/guess -- see check_cli_grad.F90's own header.
# Needs the terco_cli executable itself (app/CMakeLists.txt, added AFTER
# this directory by the root CMakeLists.txt, but still within the same
# configure run -- $<TARGET_FILE:...> is a generator expression, resolved
# at generate time, so the forward reference is fine); TERCO_BUILD_APP is
# a separate switch from TERCO_ENABLE_TESTING, so this is guarded on it
# rather than assumed.
if(TERCO_BUILD_APP)
  set(_cli_grad_probe ${CMAKE_CURRENT_BINARY_DIR}/cli_grad_probe.json)
  terco_add_check(check_cli_grad check_cli_grad.F90 NEEDS_BASIS
                  ARGS $<TARGET_FILE:terco_cli> ${_cli_grad_probe})
  # `add_dependencies(check_cli_grad terco_cli)` belongs in
  # app/CMakeLists.txt, not here: that target does not exist yet in THIS
  # directory's own configure pass (app/ is added after test/ by the root
  # CMakeLists.txt), only by generate time -- which is exactly why the
  # ARGS above is a generator expression and the explicit build-order
  # dependency is not.
  unset(_cli_grad_probe)

  # The same comparison, for a KS-DFT gradient (--func pbe): the CLI's
  # --grad is NOT refused for an LDA/GGA functional any more (this
  # branch's own lift of trc_run_gradient's DFT refusal) -- see
  # check_cli_grad_dft.F90's own header.
  set(_cli_grad_dft_probe ${CMAKE_CURRENT_BINARY_DIR}/cli_grad_dft_probe.json)
  terco_add_check(check_cli_grad_dft check_cli_grad_dft.F90 NEEDS_BASIS
                  ARGS $<TARGET_FILE:terco_cli> ${_cli_grad_dft_probe})
  unset(_cli_grad_dft_probe)

  # The same comparison for --df --grad (B3LYP: fitted J, fitted K, XC) --
  # see check_cli_grad_df.F90's own header.
  set(_cli_grad_df_probe ${CMAKE_CURRENT_BINARY_DIR}/cli_grad_df_probe.json)
  terco_add_check(check_cli_grad_df check_cli_grad_df.F90 NEEDS_BASIS
                  ARGS $<TARGET_FILE:terco_cli> ${_cli_grad_df_probe})
  unset(_cli_grad_df_probe)
endif()

# The FULL KS-DFT gradient (the above, plus trc_xc_grad's LDA/GGA
# exchange-correlation term) against finite differences of terco's own
# fully reconverged KS-DFT energy, at two grid levels -- RKS SVWN/PBE on
# water, UKS PBE on the OH radical -- see check_xc_grad_fd.F90's header.
terco_add_check(check_xc_grad_fd check_xc_grad_fd.F90 NEEDS_BASIS)
# LONG: over half of the whole suite's time (1091 of 1973 s on a gcc-15 CI
# runner, 1920 s on a slow one), so pull requests skip it (`ctest -LE LONG`)
# and the nightly scheduled CI run on main runs it with everything else.
set_tests_properties(terco/check_xc_grad_fd PROPERTIES LABELS LONG)

# The VV10 gradient's own formula, isolated from the SCF entirely: a
# fixed-density finite difference of E_nlc[D] alone -- see
# check_vv10_grad_fixedd_fd.F90's own header for why this is a sharper
# gate than check_xc_grad_fd's reconverged-SCF wB97X-V case.
terco_add_check(check_vv10_grad_fixedd_fd check_vv10_grad_fixedd_fd.F90 NEEDS_BASIS)

# trc_xc_grad's basis-function (Pulay) term ALONE, against a finite
# difference that holds the grid itself fixed too -- no SCF, no grid
# response on either side of the comparison -- see
# check_xc_grad_fixedgrid_fd.F90's own header for why this is the test
# that pins the formula down to ~1e-9, independent of check_xc_grad_fd's
# SCF-based (and therefore noisier) comparison above.
terco_add_check(check_xc_grad_fixedgrid_fd check_xc_grad_fixedgrid_fd.F90 NEEDS_BASIS)

#
# Cartesian multipoles through the octupole, against a reference
# test/pyscf_refs/generate.py recorded from pyscf. No pyscf in CI, ever --
# see test/pyscf_refs/README.md: the Fortran program below reads
# test/pyscf_refs/mult.ref back and does the comparison itself.
#
terco_add_check(check_mult check_mult.F90 NEEDS_PYSCF_REF)

#
# Basis function values and gradients, elementwise against a committed
# pyscf eval_ao reference (test/pyscf_refs/ao.ref).
#
terco_add_check(check_ao check_ao.F90 NEEDS_PYSCF_REF)

# E_xc and V_xc on the same grid and density, twenty functionals, RKS and
# UKS, against test/pyscf_refs/xc_energy.ref.
terco_add_check(check_xc_energy check_xc_energy.F90 NEEDS_PYSCF_REF)

# VV10's E_nlc and potential on the same grid and density, against
# test/pyscf_refs/vv10.ref: water/6-31G* and water/cc-pVDZ (RKS), the OH
# radical/cc-pVDZ (UKS).
terco_add_check(check_vv10_vs_pyscf check_vv10_vs_pyscf.F90 NEEDS_BASIS NEEDS_PYSCF_REF)

# eri%gradient against a committed pyscf int2e_ip1 contraction (pyscf's
# own Gradients class assumes a converged SCF density; this is a fixed,
# non-SCF density, same reasoning as check_rys_grad_fd). RHF and UHF,
# water/6-31G* and water/cc-pVDZ.
foreach(b 631gs ccpvdz)
  if(b STREQUAL "631gs")
    set(json basis_sets/6-31gs.json)
  else()
    set(json basis_sets/cc-pvdz.json)
  endif()
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/rys_grad_probe_${b}.bin)
  terco_add_check(check_rys_grad_vs_pyscf_${b} check_rys_grad_vs_pyscf.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${json} ${dump} ${b})
endforeach()

# The FULL KS-DFT gradient against pyscf's own grad.RKS/grad.UKS, at
# terco's own converged orbitals AND terco's own grid (see
# check_xc_grad_vs_pyscf.F90's header) -- SVWN/PBE/B3LYP on water and
# formamide, UKS PBE on the OH radical.
foreach(case water_svwn_631gs water_pbe_631gs water_b3lyp_631gs
             water_svwn_ccpvdz water_pbe_ccpvdz water_b3lyp_ccpvdz
             formamide_pbe_ccpvdz formamide_b3lyp_ccpvdz
             oh_pbe_631gs)
  if(case STREQUAL "water_svwn_631gs")
    set(xyz water.xyz)
    set(json basis_sets/6-31gs.json)
    set(func svwn)
    set(glev 3)
  elseif(case STREQUAL "water_pbe_631gs")
    set(xyz water.xyz)
    set(json basis_sets/6-31gs.json)
    set(func pbe)
    set(glev 3)
  elseif(case STREQUAL "water_b3lyp_631gs")
    set(xyz water.xyz)
    set(json basis_sets/6-31gs.json)
    set(func b3lyp)
    set(glev 3)
  elseif(case STREQUAL "water_svwn_ccpvdz")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func svwn)
    set(glev 3)
  elseif(case STREQUAL "water_pbe_ccpvdz")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func pbe)
    set(glev 3)
  elseif(case STREQUAL "water_b3lyp_ccpvdz")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func b3lyp)
    set(glev 3)
  elseif(case STREQUAL "formamide_pbe_ccpvdz")
    set(xyz formamide.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func pbe)
    set(glev 3)
  elseif(case STREQUAL "formamide_b3lyp_ccpvdz")
    set(xyz formamide.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func b3lyp)
    set(glev 3)
  else()
    set(xyz oh.xyz)
    set(json basis_sets/6-31gs.json)
    set(func pbe)
    set(glev 3)
  endif()
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/xc_grad_probe_${case}.bin)
  terco_add_check(check_xc_grad_vs_pyscf_${case} check_xc_grad_vs_pyscf.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${xyz} ${json} ${func} ${glev} ${dump} ${case})
endforeach()

# The SAME full KS-DFT gradient, grid_response=.true. on BOTH sides (see
# check_xc_grad_vs_pyscf.F90's header, and generate.py's gen_xc_grad
# docstring on why pyscf's own grid is built fresh here rather than
# overridden) -- RKS SVWN/PBE on water, UKS PBE on NH2 (a genuinely
# open-shell, non-degenerate-SOMO case -- check_xc_grad_fd.F90's own
# header explains why NH2 and not OH).
foreach(case water_svwn_631gs water_pbe_631gs nh2_pbe_631gs
             water_wb97x_631gs water_wb97x_ccpvdz)
  if(case STREQUAL "water_svwn_631gs")
    set(xyz water.xyz)
    set(json basis_sets/6-31gs.json)
    set(func svwn)
    set(glev 3)
  elseif(case STREQUAL "water_pbe_631gs")
    set(xyz water.xyz)
    set(json basis_sets/6-31gs.json)
    set(func pbe)
    set(glev 3)
  elseif(case STREQUAL "water_wb97x_631gs")
    # Range-separated hybrid (wB97X: nonzero omega/exx_lr -- see
    # src/xc/trc_xc_grad.F90's own header on why lifting the refusal
    # there is safe): the semilocal GGA part is this grid-response term,
    # the erf(omega r)/r long-range exchange is trc_scf_grad.F90's own
    # (grad-rsh's) eri%gradient second pass, grid-independent.
    set(xyz water.xyz)
    set(json basis_sets/6-31gs.json)
    set(func wb97x)
    set(glev 3)
  elseif(case STREQUAL "water_wb97x_ccpvdz")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
    set(func wb97x)
    set(glev 3)
  else()
    set(xyz nh2.xyz)
    set(json basis_sets/6-31gs.json)
    set(func pbe)
    set(glev 3)
  endif()
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/xc_grad_probe_gridresp_${case}.bin)
  terco_add_check(check_xc_grad_vs_pyscf_gridresp_${case} check_xc_grad_vs_pyscf.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${xyz} ${json} ${func} ${glev} ${dump} ${case} 1)
endforeach()

# The long-range exchange gradient term alone (j_scale=0, k_scale=1, omega)
# against a committed pyscf int2e_ip1 contraction under
# mol.with_range_coulomb(omega) -- the omega > 0 analogue of
# check_rys_grad_vs_pyscf, RHF only. Two omega values (0.3, 0.4) in one
# binary/reference -- see check_lr_exchange_grad_vs_pyscf.F90's own header.
foreach(b 631gs ccpvdz)
  if(b STREQUAL "631gs")
    set(json basis_sets/6-31gs.json)
  else()
    set(json basis_sets/cc-pvdz.json)
  endif()
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/lr_exchange_grad_probe_${b}.bin)
  terco_add_check(check_lr_exchange_grad_vs_pyscf_${b} check_lr_exchange_grad_vs_pyscf.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${json} ${dump} ${b})
endforeach()

# trc_1e_gradient against committed pyscf int1e_ipovlp/int1e_ipkin/
# int1e_ipnuc(+Hellmann-Feynman) contractions, each term separated (see
# check_1e_grad_vs_pyscf.F90's own header), at a fixed non-SCF density --
# water/6-31G*, water/cc-pVDZ, formamide/cc-pVDZ.
foreach(case water_631gs water_ccpvdz formamide_ccpvdz)
  if(case STREQUAL "water_631gs")
    set(xyz water.xyz)
    set(json basis_sets/6-31gs.json)
  elseif(case STREQUAL "water_ccpvdz")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
  else()
    set(xyz formamide.xyz)
    set(json basis_sets/cc-pvdz.json)
  endif()
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/onee_grad_probe_${case}.bin)
  terco_add_check(check_1e_grad_vs_pyscf_${case} check_1e_grad_vs_pyscf.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${xyz} ${json} ${dump} ${case})
endforeach()

#
# eri%fock against a committed pyscf raw int2e reference, on a fixed
# density -- the first comparison that can exercise an f shell (l=3).
# check_rys_fock_vs_pyscf.F90 skips itself (stop 0, after printing SKIP)
# when the basis needs more l than this build's TRC_LMAX reaches, exactly
# the way check_nosym_f does, so it is registered unconditionally and at
# every ceiling: a build that CAN reach f always runs it -- this exercises the ordinary eri%fock dispatch,
# whichever scheme ends up routing an f class once TERCO_LMAX reaches 3.
#
# water/cc-pVDZ has no f on O and passes TODAY, at LMAX=2: the proof that
# the AO permutation and the per-AO normalisation this shares with the
# rys_grad case (pyscf's cart=True does not keep a JSON's literal shell
# order, and its Cartesian d/f components are not all unit-normalised the
# way this build's are) are right before any f kernel exists to test.
# water/cc-pVTZ (f on O) and formamide/def2-TZVP (f on C, N and O --
# formamide.xyz, an idealised planar geometry) are the f cases: at
# TERCO_LMAX=3 they run for real against test/pyscf_refs/rys_fock_*.ref,
# which the LMAX=3 CI job requires to exist (see check_rys_fock_vs_pyscf.F90's
# own FAIL-not-SKIP branch for a missing reference at that ceiling).
#
foreach(case water_ccpvdz water_ccpvtz formamide_def2tzvp)
  # The gate is 1e-9 by default, tightened per case where measurement
  # leaves ample margin: water/cc-pVTZ measures 8e-13, four orders of
  # magnitude under even 1e-10, so it is passed explicitly here.
  # Formamide/def2-TZVP measures 2.6e-11 -- a real size effect (more atoms,
  # more AOs, not noise) -- and keeps the 1e-9 default rather than risk a
  # false red on a larger system tested later.
  set(tol 1e-9)
  if(case STREQUAL "water_ccpvdz")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
  elseif(case STREQUAL "water_ccpvtz")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvtz.json)
    set(tol 1e-10)
  else()
    set(xyz formamide.xyz)
    set(json basis_sets/def2-tzvp.json)
  endif()
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/rys_fock_probe_${case}.bin)
  terco_add_check(check_rys_fock_vs_pyscf_${case} check_rys_fock_vs_pyscf.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${xyz} ${json} ${dump} ${tol})
  # The build-side LMAX self-skip still prints "SKIP" and exits 0 at
  # TERCO_LMAX < 3 -- ctest calls that "Skipped", not "Passed". At
  # TERCO_LMAX>=3 a MISSING reference is a hard FAIL instead (see
  # ref_require in pyscf_ref_io.F90), never silently skipped.
  set_tests_properties(terco/check_rys_fock_vs_pyscf_${case}
    PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")
endforeach()

#
# omega > 0 (erf(0.3*r12)/r12, wB97X-V's long-range exchange operator)
# against a committed pyscf range-separated J/K
# (test/pyscf_refs/erf_fock_*.ref, from mol.with_range_coulomb/scf.hf.get_jk
# -- see check_erf_fock_vs_pyscf.F90's own header for the J/K-isolation
# trick). water/6-31G* and water/cc-pVDZ exercise s/p only; formamide/cc-pVDZ
# puts d functions on every heavy atom. Each case runs three ways:
# eri%fock through HGP (TRC_RYS_CLASSES=none), eri%fock through Rys
# (TRC_RYS_CLASSES=all), and eri%fock_nosym (the enumerated/general path) --
# the same substitution has to be right in all three, and all three compare
# against the SAME reference file (the fixed density and pyscf's J/K do not
# depend on which terco kernel path produced the number under test).
#
foreach(case water_631gs water_ccpvdz formamide_ccpvdz)
  if(case STREQUAL "water_631gs")
    set(xyz water.xyz)
    set(json basis_sets/6-31gs.json)
  elseif(case STREQUAL "water_ccpvdz")
    set(xyz water.xyz)
    set(json basis_sets/cc-pvdz.json)
  else()
    set(xyz formamide.xyz)
    set(json basis_sets/cc-pvdz.json)
  endif()
  foreach(path hgp nosym)
    set(dump ${CMAKE_CURRENT_BINARY_DIR}/erf_fock_probe_${case}_${path}.bin)
    if(path STREQUAL "nosym")
      set(nosym_flag 1)
    else()
      set(nosym_flag 0)
    endif()
    terco_add_check(check_erf_fock_vs_pyscf_${case}_${path}
                    check_erf_fock_vs_pyscf.F90
                    NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${xyz} ${json} ${dump} ${nosym_flag} ${case})
  endforeach()
  # "hgp" is HGP for every class: the default sends some d classes to Rys.
  set_tests_properties(terco/check_erf_fock_vs_pyscf_${case}_hgp
    PROPERTIES ENVIRONMENT "CUDA_VISIBLE_DEVICES=0;TERCO_CART_ORDER=${TERCO_CART_ORDER};TRC_RYS_CLASSES=none")
  # Same as the "hgp" case above, but routed through the Rys kernels
  # instead of HGP's (TRC_RYS_CLASSES=all).
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/erf_fock_probe_${case}_rys.bin)
  terco_add_check(check_erf_fock_vs_pyscf_${case}_rys
                  check_erf_fock_vs_pyscf.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${xyz} ${json} ${dump} 0 ${case})
  set_tests_properties(terco/check_erf_fock_vs_pyscf_${case}_rys
    PROPERTIES ENVIRONMENT "CUDA_VISIBLE_DEVICES=0;TERCO_CART_ORDER=${TERCO_CART_ORDER};TRC_RYS_CLASSES=all")
endforeach()

#
# trc_df_2c ((P|Q)) and trc_df_3c ((mu nu|P)) against a committed pyscf
# int2c2e/aux_e2(int3c2e) reference (test/pyscf_refs/df_vs_pyscf.ref) -- the
# density-fitting analogue of the eri%fock comparison above, and the first
# one that can exercise a g SHELL (l=4), on the auxiliary basis.
# check_df_vs_pyscf.F90 skips itself (stop 0, after printing SKIP) when
# either basis needs more l than this build's DF_LMAX/DF_LMAX_AUX reach
# (gen_df.py's two ceilings, src/CMakeLists.txt's TERCO_DF_LMAX), same
# convention as check_rys_fock_vs_pyscf.F90's TRC_LMAX skip -- registered
# unconditionally, nothing to run at the shipped ceiling (DF_LMAX=2,
# DF_LMAX_AUX=3), and running water/cc-pVTZ (f on O) against cc-pVTZ-RI
# (g on O, scripts/pyscf_basis_to_json.py's cc-pvtz-rifit.json) for real
# once DF_LMAX reaches 3 -- a missing reference there is a hard FAIL.
#
set(dump ${CMAKE_CURRENT_BINARY_DIR}/df_vs_pyscf_probe.bin)
terco_add_check(check_df_vs_pyscf check_df_vs_pyscf.F90
                NEEDS_BASIS NEEDS_PYSCF_REF ARGS water.xyz basis_sets/cc-pvtz.json
                basis_sets/cc-pvtz-rifit.json ${dump})
# Same self-skip convention as check_rys_fock_vs_pyscf above, one level up.
set_tests_properties(terco/check_df_vs_pyscf
  PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")

#
# The density-fitting DERIVATIVE integrals (scripts/gen_df_grad.py) against
# pyscf's int3c2e_ip1/ip2 and int2c2e_ip1: raw one-hot tensors from the host
# class dispatch, and the device drivers' Y/Z contraction against both those
# and test/pyscf_refs/dfg_ints_<case>.ref -- see check_df_grad_ints.F90's
# header. def2-universal-jkfit carries g auxiliary shells, so this covers the
# shipped DFG_LMAX_AUX = 4 ceiling.
#
foreach(case IN ITEMS water_631gs water_ccpvdz formamide_ccpvdz)
  string(REGEX MATCH "^[a-z]+" _mol "${case}")
  if(case MATCHES "631gs$")
    set(_bas basis_sets/6-31gs.json)
  else()
    set(_bas basis_sets/cc-pvdz.json)
  endif()
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/dfg_ints_${case}_probe.bin)
  terco_add_check(check_df_grad_ints_${case} check_df_grad_ints.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${_mol}.xyz ${_bas}
                  basis_sets/def2-universal-jkfit.json ${dump} ${case})
  set_tests_properties(terco/check_df_grad_ints_${case}
    PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")
  # Its own module directory: both cases compile the same source, whose
  # helper module would otherwise race into one shared directory.
  set_target_properties(check_df_grad_ints_${case} PROPERTIES Fortran_MODULE_DIRECTORY
                        ${CMAKE_CURRENT_BINARY_DIR}/mod_check_df_grad_ints_${case})
endforeach()

# The --df analytic gradient (fitted J and K, trc_scf_grad_df.F90) against a
# Richardson central finite difference of terco's own reconverged --df energy
# -- RHF, RKS PBE, RKS B3LYP on water/6-31G*, see check_df_grad_fd.F90.
terco_add_check(check_df_grad_fd check_df_grad_fd.F90 NEEDS_BASIS)
set_tests_properties(terco/check_df_grad_fd PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")
# The same, RHF on formamide/cc-pVDZ (cond(M) = 8.6e13), independent of
# pyscf: ~25 s on 4 CPU cores.
add_test(NAME terco/check_df_grad_fd_formamide COMMAND check_df_grad_fd formamide)
set_tests_properties(
  terco/check_df_grad_fd_formamide
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} SKIP_REGULAR_EXPRESSION "SKIP"
             ENVIRONMENT "CUDA_VISIBLE_DEVICES=0;TERCO_CART_ORDER=${TERCO_CART_ORDER}")

# The full --df gradient against pyscf's own df.grad.rhf/rks at terco's
# converged orbitals (and grid), def2-universal-jkfit, Cartesian --
# test/pyscf_refs/dfg_<case>.ref, see check_df_grad_vs_pyscf.F90.
foreach(case IN ITEMS water_hf_631gs formamide_hf_ccpvdz water_pbe_631gs water_b3lyp_631gs
                      formamide_b3lyp_ccpvdz)
  string(REGEX MATCH "^[a-z]+" _mol "${case}")
  string(REGEX MATCH "_(hf|pbe|b3lyp)_" _fn "${case}")
  string(REPLACE "_" "" _fn "${_fn}")
  if(case MATCHES "631gs$")
    set(_bas basis_sets/6-31gs.json)
  else()
    set(_bas basis_sets/cc-pvdz.json)
  endif()
  # The formamide gate is the fit's conditioning floor (measured: the gap is
  # the (P|Q) values through cond(M) = 8.6e13), see the program's `tol`
  # comment; water's is the gradient's own.
  if(_mol STREQUAL "formamide")
    set(_tol 5e-9)
  else()
    set(_tol 1e-10)
  endif()
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/dfg_${case}_probe.bin)
  terco_add_check(check_df_grad_vs_pyscf_${case} check_df_grad_vs_pyscf.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS ${_mol}.xyz ${_bas}
                  basis_sets/def2-universal-jkfit.json ${_fn} 3 ${dump} ${case} ${_tol})
  set_tests_properties(terco/check_df_grad_vs_pyscf_${case}
    PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")
endforeach()
# The --grad --df refusals (UHF, spherical, range-separated, no fit) at the C
# API, by status and message, never the library's error stop -- see
# check_grad_df_refusals.F90.
terco_add_check(check_grad_df_refusals check_grad_df_refusals.F90 NEEDS_BASIS
                ARGS ${CMAKE_CURRENT_BINARY_DIR}/grad_df_refusals.scf)

# The occupied-chunked path (trc_scf_grad_df.F90's MEMORY note), forced on
# with TERCO_DFG_OCC_CHUNK=5: formamide's 12 occupied orbitals in chunks of
# 5, 5 and 2, against the same reference.
add_test(NAME terco/check_df_grad_vs_pyscf_formamide_b3lyp_ccpvdz_chunked
         COMMAND check_df_grad_vs_pyscf_formamide_b3lyp_ccpvdz formamide.xyz basis_sets/cc-pvdz.json
                 basis_sets/def2-universal-jkfit.json b3lyp 3
                 ${CMAKE_CURRENT_BINARY_DIR}/dfg_chunked_probe.bin formamide_b3lyp_ccpvdz 5e-9)
set_tests_properties(
  terco/check_df_grad_vs_pyscf_formamide_b3lyp_ccpvdz_chunked
  PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} SKIP_REGULAR_EXPRESSION "SKIP"
             ENVIRONMENT "CUDA_VISIBLE_DEVICES=0;TERCO_CART_ORDER=${TERCO_CART_ORDER};TERCO_DFG_OCC_CHUNK=5")

#
# The polarisable continuum's cavity, S and D, and charge solve, with no
# integrals: Born-ion and sphere identities, the refusals, and pyscf's own
# surface arrays from test/pyscf_refs/pcm_cavity_water.ref.
#
terco_add_check(check_pcm_cavity check_pcm_cavity.F90 NEEDS_BASIS NEEDS_PYSCF_REF)
#
# Its solute-surface integrals (gen_1e.py's --pcm kernels) on fixed inputs
# against pyscf's int3c2e/fakemol path: water/6-31G* (d) at every ceiling,
# water/cc-pVTZ (f on O) for real at TERCO_LMAX=3 and skipped below it, the
# same convention as check_rys_fock_vs_pyscf.
#
foreach(case water_631gs water_ccpvtz)
  if(case STREQUAL "water_631gs")
    set(json basis_sets/6-31gs.json)
  else()
    set(json basis_sets/cc-pvtz.json)
  endif()
  set(dump ${CMAKE_CURRENT_BINARY_DIR}/pcm_ints_probe_${case}.bin)
  terco_add_check(check_pcm_ints_${case} check_pcm_ints.F90
                  NEEDS_BASIS NEEDS_PYSCF_REF ARGS water.xyz ${json} ${dump} ${case})
  set_tests_properties(terco/check_pcm_ints_${case} PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")
endforeach()

#
# Negative tests: trc_df_2c and trc_df_3c must refuse a basis above the
# ceiling the DF dispatchers were generated for, by `error stop`, rather
# than read an unwritten block. Nothing above exercises this -- every
# check builds within the ceiling on purpose -- so without these nothing
# notices if the guard in trc_api.F90 is ever lost. Each hands the entry
# point a shell one l above THIS build's DF_LMAX/DF_LMAX_AUX
# (trc_df_kernels), so the test is valid whatever the ceiling is, not only
# at the shipped one. PASS_REGULAR_EXPRESSION matches the error-stop
# message itself rather than WILL_FAIL, so a different crash -- a segfault
# before the guard runs -- does not count as this test passing.
#
terco_add_check(check_df_2c_ceiling check_df_2c_ceiling.F90 NEEDS_BASIS)
set_tests_properties(terco/check_df_2c_ceiling PROPERTIES
  PASS_REGULAR_EXPRESSION "trc_df_2c: auxiliary basis above DF_LMAX_AUX")
terco_add_check(check_df_3c_ceiling check_df_3c_ceiling.F90 NEEDS_BASIS)
set_tests_properties(terco/check_df_3c_ceiling PROPERTIES
  PASS_REGULAR_EXPRESSION "trc_df_3c: orbital basis above DF_LMAX")

#
# A complete RHF through the public API, and the makeEFP pathway on top of it.
#
# These need LAPACK and terco does not -- the eigenproblem is the CALLER's
# half of the contract. terco's own contract is: given a density, return a
# Fock matrix; given a basis, return one-electron matrices. Everything above
# that, including the diagonalisation and DIIS, stays with whoever is driving.
#
find_package(LAPACK REQUIRED)

add_executable(scf_rhf scf_rhf.F90)
target_compile_definitions(scf_rhf PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_rhf PRIVATE ${main_lib} terco_test_basis LAPACK::LAPACK)
target_include_directories(scf_rhf PRIVATE ${CMAKE_BINARY_DIR}/modules)

#
# The library's own SCF driver, resident on the device: water RHF and four
# functionals restricted, the OH radical UHF and PBE unrestricted, in 6-31G
# uncontracted to primitives. Writes rks_probe.bin for the same-grid pyscf
# comparison registered above.
#
add_executable(scf_rks scf_rks.F90)
target_compile_definitions(scf_rks PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_rks PRIVATE ${main_lib} terco_test_basis terco_pyscf_ref_io LAPACK::LAPACK)
target_include_directories(scf_rks PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_rks COMMAND scf_rks)
set_tests_properties(
  terco/scf_rks PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                           ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

#
# The same harness for the meta-GGAs (tau): water r2SCAN, TPSS, M06-L, SCAN,
# TPSSh, r2SCAN0, M06, M06-2X restricted, OH r2SCAN, TPSS, M06-L, TPSSh
# unrestricted, against pyscf on terco's own grid (pyscf_refs/mgga.ref).
#
add_executable(scf_mgga scf_mgga.F90)
target_compile_definitions(scf_mgga PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_mgga PRIVATE ${main_lib} terco_test_basis terco_pyscf_ref_io LAPACK::LAPACK)
target_include_directories(scf_mgga PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_mgga COMMAND scf_mgga)
set_tests_properties(
  terco/scf_mgga PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                            ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

#
# Density-fitted RHF/RKS (trc_scf_df.F90, `--df`'s library entry): water and
# formamide, cc-pVDZ and 6-31G*, Hartree-Fock/PBE/B3LYP, against pyscf's own
# density_fit on the SAME auxiliary basis (def2-universal-jkfit), at rks.ref's
# closed-shell tolerance (test/pyscf_refs/scf_df.ref, generate.py's
# gen_scf_df). def2-universal-jkfit needs g auxiliary functions on O/C/N,
# within the shipped DF_LMAX_AUX=4 ceiling since that became the default
# (see the top-level CMakeLists.txt's TERCO_DF_LMAX_AUX), so this runs for
# real in a plain build; scf_df.F90 still self-skips (stop 0, after
# printing SKIP), same convention as check_df_vs_pyscf's, for anyone who
# lowers TERCO_DF_LMAX_AUX back down.
#
add_executable(scf_df scf_df.F90)
target_compile_definitions(scf_df PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_df PRIVATE ${main_lib} terco_test_basis terco_pyscf_ref_io LAPACK::LAPACK)
target_include_directories(scf_df PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_df COMMAND scf_df)
set_tests_properties(
  terco/scf_df PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                          ENVIRONMENT "CUDA_VISIBLE_DEVICES=0"
                          SKIP_REGULAR_EXPRESSION "SKIP")

#
# The SCF with a polarisable continuum against pyscf's .PCM(): converged
# energies and E_diel, and both again at pyscf's own converged density
# (test/pyscf_refs/pcm_scf_main.ref, pcm_scf_f.ref). `main` is water and
# formamide (HF, B3LYP, wB97X-V; C-PCM and IEF-PCM), the OH radical
# unrestricted, and one spherical case; `f` is water/cc-pVTZ, run for real at
# TERCO_LMAX=3 and skipped below it.
#
add_executable(scf_pcm scf_pcm.F90)
target_compile_definitions(scf_pcm PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_pcm PRIVATE ${main_lib} terco_test_basis terco_pyscf_ref_io)
target_include_directories(scf_pcm PRIVATE ${CMAKE_BINARY_DIR}/modules)
foreach(set main f)
  add_test(NAME terco/scf_pcm_${set} COMMAND scf_pcm ${set})
  set_tests_properties(
    terco/scf_pcm_${set} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                    ENVIRONMENT "CUDA_VISIBLE_DEVICES=0"
                                    SKIP_REGULAR_EXPRESSION "SKIP")
endforeach()

#
# The hard-case set for SCF robustness (M2, Part B): A/EDIIS ahead of CDIIS
# and canonical orthogonalisation at their measured defaults, on four small
# systems chosen because plain damped DIIS is weak on each for a different
# reason -- see scf_robust.F90's header. Compares against pyscf's own
# converged energy, committed in test/pyscf_refs/scf_robust.ref.
#
add_executable(scf_robust scf_robust.F90)
target_compile_definitions(scf_robust PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_robust PRIVATE ${main_lib} terco_test_basis terco_pyscf_ref_io LAPACK::LAPACK)
target_include_directories(scf_robust PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_robust COMMAND scf_robust)
set_tests_properties(
  terco/scf_robust PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                              ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

# The SAD guess read from mqc's JSON blocks, against terco's own start.
add_executable(scf_sad scf_sad.F90)
target_compile_definitions(scf_sad PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_sad PRIVATE ${main_lib} terco_test_basis LAPACK::LAPACK)
target_include_directories(scf_sad PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_sad COMMAND scf_sad)
set_tests_properties(
  terco/scf_sad PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                           ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

# RI-MP2 on the RHF reference, against pyscf's DF-MP2 on the same sets.
# The shipped C-interface declarations, compiled once: two targets compiling
# the same source both write trc_c_interfaces.mod, and Ninja refuses that.
add_library(terco_c_interfaces STATIC ${PROJECT_SOURCE_DIR}/include/trc_c_interfaces.f90)
set_target_properties(terco_c_interfaces PROPERTIES Fortran_MODULE_DIRECTORY
                                                    ${CMAKE_CURRENT_BINARY_DIR}/c_interfaces)
target_include_directories(terco_c_interfaces INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/c_interfaces)

add_executable(scf_rimp2 scf_rimp2.F90)
target_compile_definitions(scf_rimp2 PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_rimp2 PRIVATE ${main_lib} terco_test_basis terco_c_interfaces terco_pyscf_ref_io LAPACK::LAPACK)

# The context surface against the Fortran drivers, through the compiled
# interface declarations like scf_rimp2.
add_executable(context_scf context_scf.F90)
target_compile_definitions(context_scf PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(context_scf PRIVATE ${main_lib} terco_test_basis terco_c_interfaces LAPACK::LAPACK)
add_test(NAME terco/context_scf COMMAND context_scf)
set_tests_properties(terco/context_scf PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                  ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

# trc_set_basis_libcint/trc_set_aux_libcint/trc_set_jkfit_libcint with
# cartesian = 0 (spherical), against the equivalent trc_set_*_json path.
add_executable(check_capi_libcint_sph check_capi_libcint_sph.F90)
target_compile_definitions(check_capi_libcint_sph PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(check_capi_libcint_sph PRIVATE ${main_lib} terco_test_basis terco_c_interfaces LAPACK::LAPACK)
add_test(NAME terco/check_capi_libcint_sph COMMAND check_capi_libcint_sph)
set_tests_properties(terco/check_capi_libcint_sph PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                  ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

#
# trc_run_polarizability/trc_polarizability (CPHF/CPKS) and trc_nmo/
# trc_mo_occupations, through the context API, against the SAME pyscf
# references test/scf_cphf.F90 and test/scf_cpks.F90 check their own
# Fortran drivers with.
add_executable(check_capi_polarizability check_capi_polarizability.F90)
target_compile_definitions(check_capi_polarizability PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(check_capi_polarizability PRIVATE ${main_lib} terco_test_basis terco_c_interfaces
                                                        terco_pyscf_ref_io)
add_test(NAME terco/check_capi_polarizability COMMAND check_capi_polarizability)
set_tests_properties(terco/check_capi_polarizability PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                  ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

#
# trc_set_unrestricted, trc_set_level_shift and trc_set_charge, through
# the context API: forced UHF on a closed shell matches RHF, a level
# shift only changes the path (not the converged energy), and an
# ionisation series on one --df handle matches fresh handles while
# reusing the fit.
add_executable(check_capi_charge_unrestricted check_capi_charge_unrestricted.F90)
target_compile_definitions(check_capi_charge_unrestricted PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(check_capi_charge_unrestricted PRIVATE ${main_lib} terco_test_basis terco_c_interfaces)
add_test(NAME terco/check_capi_charge_unrestricted COMMAND check_capi_charge_unrestricted)
set_tests_properties(terco/check_capi_charge_unrestricted PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                  ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

#
# trc_set_pcm_radii, through the context API: an override equal to
# terco's own table is bit-identical to the table path; an override at
# mqc's own radii convention (plain Bondi x 1.2) matches an independent
# pyscf reference built with the same convention.
add_executable(check_capi_pcm_radii check_capi_pcm_radii.F90)
target_compile_definitions(check_capi_pcm_radii PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(check_capi_pcm_radii PRIVATE ${main_lib} terco_test_basis terco_c_interfaces
                                                    terco_pyscf_ref_io)
add_test(NAME terco/check_capi_pcm_radii COMMAND check_capi_pcm_radii)
set_tests_properties(terco/check_capi_pcm_radii PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                  ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

#
# terco driven from C, through include/terco.h only -- the C-side proof
# trc_c_interfaces.f90 is for a Fortran caller (test/link_gfortran.F90).
# Needs the C compiler, which nothing else in this build does; enabled
# here, scoped to this one target, rather than at the top level.
#
enable_language(C)
add_executable(link_c link_c.c)
# Linked against libterco's own FILE, not the CMake target name: going
# through the target name pulls in terco's PUBLIC usage requirements,
# including -cudalib=cublas,cusolver under TERCO_ENABLE_GPU (top-level
# CMakeLists.txt) -- an nvfortran flag the C compiler linking this
# executable does not understand. Linking nvfortran itself as the
# driver instead does understand that flag, but then supplies its OWN
# `main` (expecting a Fortran program's `MAIN_` entry), which collides
# with link_c.c's own `main` -- worse, not better. The file path sidesteps
# both: the shared library already carries cuBLAS/cuSOLVER as its own
# linked dependencies, so nothing here needs to ask for them again.
add_dependencies(link_c ${main_lib})
target_link_libraries(link_c PRIVATE $<TARGET_FILE:${main_lib}>)
target_include_directories(link_c PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set_target_properties(link_c PROPERTIES BUILD_RPATH "$<TARGET_FILE_DIR:${main_lib}>")
add_test(NAME terco/link_c COMMAND link_c)
set_tests_properties(terco/link_c PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                             ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

# GAMESS's per-component AO normalisation (trc_set_ao_norm): the factor
# itself against scripts/cartorder.py's recorded reference, l = 0..4, then
# the overlap diagonal, the SCF energy, Tr(D S) and a guess round trip, all
# in both TRC_AO_NORM_SHELL and TRC_AO_NORM_UNIT, on water/cc-pVDZ (real d
# functions). See check_ao_norm.F90's own header for why each check is
# there and reference_ao_norm.dat's header for how it was generated.
add_executable(check_ao_norm check_ao_norm.F90)
target_compile_definitions(check_ao_norm PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(check_ao_norm PRIVATE ${main_lib} terco_test_basis terco_c_interfaces)
target_include_directories(check_ao_norm PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/check_ao_norm COMMAND check_ao_norm)
set_tests_properties(terco/check_ao_norm PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                    ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

target_include_directories(scf_rimp2 PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_rimp2 COMMAND scf_rimp2)
set_tests_properties(
  terco/scf_rimp2 PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                             ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

# The same with both sets from MolSSI BSE JSON, the files mqc ships: cc-pVDZ
# in its Cartesian form with cc-pVDZ-RI, general contractions split per
# column as the mqc bridge splits them.
add_test(NAME terco/scf_rimp2_ccpvdz
         COMMAND scf_rimp2 water.xyz basis_sets/cc-pvdz.json basis_sets/cc-pvdz-rifit.json
                 rimp2_ccpvdz_probe.bin)
set_tests_properties(
  terco/scf_rimp2_ccpvdz PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                    ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

# UNRESTRICTED RI-MP2 (trc_urimp2_run): OH, triplet O2, NH2 in 6-31G* and
# cc-pVDZ with cc-pVDZ-RIFIT against pyscf's own UHF + DF-UMP2, plus a
# closed-shell-as-UHF self-check against the restricted path.
add_executable(scf_rimp2_uhf scf_rimp2_uhf.F90)
target_compile_definitions(scf_rimp2_uhf PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_rimp2_uhf PRIVATE ${main_lib} terco_test_basis terco_pyscf_ref_io LAPACK::LAPACK)
target_include_directories(scf_rimp2_uhf PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_rimp2_uhf COMMAND scf_rimp2_uhf)
set_tests_properties(
  terco/scf_rimp2_uhf PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                 ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

# Double hybrids (b2plyp, b2gp-plyp -- trc_xc_functional.F90): the SCF
# (four-centre AND --df, cc-pVDZ-RIFIT doing double duty as --df's jkfit
# too) then RI-MP2 (plain and chained off --df's fit) on water and
# formamide, cc-pVDZ and 6-31G*, plus OH/b2plyp UKS from a fixed density,
# against pyscf's own custom-hybrid KS ('0.53*HF + 0.47*B88, 0.73*LYP')
# plus DF-MP2/DF-UMP2 (test/pyscf_refs/dh.ref, generate.py's gen_dh).
add_executable(scf_dh scf_dh.F90)
target_compile_definitions(scf_dh PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_dh PRIVATE ${main_lib} terco_test_basis terco_pyscf_ref_io LAPACK::LAPACK)
target_include_directories(scf_dh PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_dh COMMAND scf_dh)
set_tests_properties(
  terco/scf_dh PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                          ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

# SPHERICAL-basis RI-MP2 and `--df`: restricted RI-MP2 (water, formamide,
# cc-pVDZ + cc-pVDZ-RIFIT), unrestricted RI-MP2 (OH, from a fixed initial
# density), `--df` HF/B3LYP (water, formamide, def2-universal-jkfit) and
# B2PLYP (water, four-centre AND --df, chained) -- all spherical, against
# pyscf's own `cart=False` on the identical sets (test/pyscf_refs/
# sph_corr.ref, generate.py's gen_sph_corr). The Cartesian counterpart of
# each of these is scf_rimp2(_uhf)/scf_df/scf_dh, unchanged by this file.
add_executable(scf_sph_corr scf_sph_corr.F90)
target_compile_definitions(scf_sph_corr PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_sph_corr PRIVATE ${main_lib} terco_test_basis terco_pyscf_ref_io LAPACK::LAPACK)
target_include_directories(scf_sph_corr PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_sph_corr COMMAND scf_sph_corr)
set_tests_properties(
  terco/scf_sph_corr PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

#
# Coupled-perturbed HF (src/scf/trc_cphf.F90): the static dipole
# polarizability, RHF on water/formamide and UHF on the OH radical, each
# against test/pyscf_refs/generate.py's gen_cphf/gen_cphf_uhf
# (pyscf.prop.polarizability). Four RHF cases share scf_cphf.F90, one
# executable per (molecule, basis) pair -- a different name/ARGS per
# terco_add_check call, same source -- and the UHF case is its own
# program (trc_cphf_solve_uhf's coupled spin channels).
#
terco_add_check(scf_cphf scf_cphf.F90 NEEDS_BASIS NEEDS_PYSCF_REF)
terco_add_check(scf_cphf_ccpvdz scf_cphf.F90 NEEDS_BASIS NEEDS_PYSCF_REF
                 ARGS water.xyz basis_sets/cc-pvdz.json water cc-pvdz cphf_ccpvdz_probe.bin)
terco_add_check(scf_cphf_formamide scf_cphf.F90 NEEDS_BASIS NEEDS_PYSCF_REF
                 ARGS formamide.xyz basis_sets/6-31gs.json formamide 6-31gs cphf_formamide_probe.bin)
terco_add_check(scf_cphf_formamide_ccpvdz scf_cphf.F90 NEEDS_BASIS NEEDS_PYSCF_REF
                 ARGS formamide.xyz basis_sets/cc-pvdz.json formamide cc-pvdz cphf_formamide_ccpvdz_probe.bin)
terco_add_check(scf_cphf_uhf scf_cphf_uhf.F90 NEEDS_BASIS NEEDS_PYSCF_REF)

#
# Coupled-perturbed KS (src/scf/trc_cpks.F90): the static dipole
# polarizability, RKS on water and formamide for SVWN (LDA), PBE and
# B3LYP (GGA/hybrid-GGA), 6-31G*, against test/pyscf_refs/generate.py's
# gen_cpks (pyscf.prop.polarizability.rks, grid matched to terco's own --
# see scf_cpks.F90's own header). One executable, one ARGS set per case.
#
terco_add_check(scf_cpks scf_cpks.F90 NEEDS_BASIS NEEDS_PYSCF_REF
                 ARGS water.xyz basis_sets/6-31gs.json water 6-31gs svwn cpks_svwn_probe.bin)
terco_add_check(scf_cpks_pbe scf_cpks.F90 NEEDS_BASIS NEEDS_PYSCF_REF
                 ARGS water.xyz basis_sets/6-31gs.json water 6-31gs pbe cpks_pbe_probe.bin)
terco_add_check(scf_cpks_b3lyp scf_cpks.F90 NEEDS_BASIS NEEDS_PYSCF_REF
                 ARGS water.xyz basis_sets/6-31gs.json water 6-31gs b3lyp cpks_b3lyp_probe.bin)
terco_add_check(scf_cpks_formamide_pbe scf_cpks.F90 NEEDS_BASIS NEEDS_PYSCF_REF
                 ARGS formamide.xyz basis_sets/6-31gs.json formamide 6-31gs pbe cpks_formamide_pbe_probe.bin)

#
# The SCF over ranks. Built in every configuration -- pic-mpi's serial
# backend makes it rank 0 of 1 -- and registered on one, two and four ranks
# when MPI is on, each run asserting the single-rank energies. No
# CUDA_VISIBLE_DEVICES here: the point is one card per rank.
#
# step_table_t%reduce: a step table's time, work and bytes over the ranks
# (RI-MP2's table is the job's, not rank 0's share). A test-drive suite;
# registered beside scf_mpi, on the same launcher.
add_executable(check_step_table_mpi check_step_table_mpi.F90)
target_link_libraries(check_step_table_mpi PRIVATE ${main_lib} test-drive::test-drive)
target_include_directories(check_step_table_mpi PRIVATE ${CMAKE_BINARY_DIR}/modules)

add_executable(scf_mpi scf_mpi.F90)
target_compile_definitions(scf_mpi PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_mpi PRIVATE ${main_lib} terco_test_basis terco_c_interfaces)
target_include_directories(scf_mpi PRIVATE ${CMAKE_BINARY_DIR}/modules)
if(TERCO_ENABLE_MPI)
  # The launcher beside the compiler wrapper, not the first mpiexec on PATH:
  # an executable linked against NVHPC's OpenMPI launched by another
  # OpenMPI's mpiexec fails at load with an undefined MPI symbol, because
  # the launcher points it at the wrong runtime.
  get_filename_component(_mpi_bindir "${CMAKE_Fortran_COMPILER}" DIRECTORY)
  if(EXISTS "${_mpi_bindir}/mpiexec")
    set(_mpiexec "${_mpi_bindir}/mpiexec")
    set(_npflag "-np")
  else()
    find_package(MPI REQUIRED COMPONENTS Fortran)
    set(_mpiexec "${MPIEXEC_EXECUTABLE}")
    set(_npflag "${MPIEXEC_NUMPROC_FLAG}")
  endif()
  foreach(np 1 2 4)
    add_test(NAME terco/scf_mpi_np${np}
             COMMAND ${_mpiexec} ${_npflag} ${np} $<TARGET_FILE:scf_mpi>)
    set_tests_properties(terco/scf_mpi_np${np} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endforeach()
  foreach(np 1 2 4)
    add_test(NAME terco/check_step_table_mpi_np${np}
             COMMAND ${_mpiexec} ${_npflag} ${np} $<TARGET_FILE:check_step_table_mpi>)
  endforeach()
else()
  add_test(NAME terco/scf_mpi COMMAND scf_mpi)
  add_test(NAME terco/check_step_table_mpi COMMAND check_step_table_mpi)
  set_tests_properties(terco/scf_mpi PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")
endif()

#
# GAMESS's DDI and a fragment method, imitated: see scf_mpi_frag.F90's own
# header for what each of the four scenarios checks. Needs exactly 4 ranks
# -- two DDI compute ranks and two data servers, or two 2-rank fragments --
# so it is only meaningful, and only registered, under MPI; it still builds
# in a non-MPI configuration (uniformly with every other test here) and
# just says so and exits. A TIMEOUT, not only the barrier inside the test
# itself: a regression that brings back the MPI_COMM_WORLD hang this test
# exists to catch should fail ctest, not wedge it.
#
add_executable(scf_mpi_frag scf_mpi_frag.F90)
target_compile_definitions(scf_mpi_frag PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_mpi_frag PRIVATE ${main_lib} terco_test_basis terco_c_interfaces)
target_include_directories(scf_mpi_frag PRIVATE ${CMAKE_BINARY_DIR}/modules)
if(TERCO_ENABLE_MPI)
  add_test(NAME terco/scf_mpi_frag
           COMMAND ${_mpiexec} ${_npflag} 4 $<TARGET_FILE:scf_mpi_frag>)
  set_tests_properties(terco/scf_mpi_frag PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                     TIMEOUT 300)
endif()

add_executable(demo_makefp demo_makefp.F90)
target_compile_definitions(demo_makefp PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(demo_makefp PRIVATE ${main_lib} terco_test_basis LAPACK::LAPACK)
target_include_directories(demo_makefp PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/demo_makefp COMMAND demo_makefp water.xyz)
set_tests_properties(
  terco/demo_makefp PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                               ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")

# Benchmarks are built but not registered as tests: a timing is not a pass.
add_executable(bench_binwater bench_binwater.F90)
target_compile_definitions(bench_binwater PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_binwater PRIVATE ${main_lib})
add_executable(bench_many bench_many.F90)
target_compile_definitions(bench_many PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_many PRIVATE ${main_lib} terco_test_basis)
add_executable(bench_xc bench_xc.F90)
target_compile_definitions(bench_xc PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_xc PRIVATE ${main_lib} terco_test_basis)
add_executable(bench_vv10 bench_vv10.F90)
target_compile_definitions(bench_vv10 PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_vv10 PRIVATE ${main_lib} terco_test_basis)
# The VV10 GRADIENT's own time next to the rest of the wB97X-V gradient --
# see bench_vv10_grad.F90's own header. Not a pass/fail, like its siblings.
add_executable(bench_vv10_grad bench_vv10_grad.F90)
target_compile_definitions(bench_vv10_grad PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_vv10_grad PRIVATE ${main_lib} terco_test_basis)
add_executable(bench_gc bench_gc.F90)
target_compile_definitions(bench_gc PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_gc PRIVATE ${main_lib} terco_test_basis)
# The spherical transform's cost against one Fock build -- see
# bench_sph.F90's own header. Not a pass/fail: a timing, like its siblings.
add_executable(bench_sph bench_sph.F90)
target_compile_definitions(bench_sph PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_sph PRIVATE ${main_lib} terco_test_basis)
# The one-electron gradient's own time next to eri%gradient's two-electron
# part -- see bench_grad_1e.F90's own header. Not a pass/fail, like its
# siblings.
add_executable(bench_grad_1e bench_grad_1e.F90)
target_compile_definitions(bench_grad_1e PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_grad_1e PRIVATE ${main_lib} terco_test_basis)
foreach(b bench_binwater bench_many bench_xc bench_gc bench_sph bench_grad_1e)
  target_include_directories(${b} PRIVATE ${CMAKE_BINARY_DIR}/modules)
endforeach()

#
# The golden tier ("Two tiers" in test/pyscf_refs/README.md): terco checked
# against its own stored numbers, at a tolerance far tighter than the pyscf
# anchors'. Each check below records its quantities (golden_io.F90) into
# golden_out/ while it runs its pyscf comparison, and terco/golden/<case>
# then diffs them against test/goldens/<case>.gold. The fixture orders the two; the check
# itself never runs twice.
#
if(Python3_FOUND)
  set(TERCO_GOLDEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/golden_out)
  file(MAKE_DIRECTORY ${TERCO_GOLDEN_OUT})
  function(terco_add_golden producer case)
    set_property(TEST terco/${producer} APPEND PROPERTY ENVIRONMENT "TERCO_GOLDEN_OUT=${TERCO_GOLDEN_OUT}")
    set_property(TEST terco/${producer} APPEND PROPERTY FIXTURES_SETUP golden_${case})
    add_test(NAME terco/golden/${case}
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/goldens/compare_goldens.py
                     ${case} ${TERCO_GOLDEN_OUT})
    set_tests_properties(terco/golden/${case} PROPERTIES FIXTURES_REQUIRED golden_${case}
                                                         SKIP_REGULAR_EXPRESSION "^SKIP")
  endfunction()
  terco_add_golden(scf_rks scf_rks)
  terco_add_golden(scf_mgga scf_mgga)
  terco_add_golden(scf_df scf_df)
  terco_add_golden(scf_robust scf_robust)
  terco_add_golden(scf_rimp2 scf_rimp2)
  terco_add_golden(scf_rimp2_ccpvdz scf_rimp2_ccpvdz)
  foreach(set main f)
    terco_add_golden(scf_pcm_${set} scf_pcm_${set})
  endforeach()
  foreach(case water_ccpvdz water_ccpvtz formamide_def2tzvp)
    terco_add_golden(check_rys_fock_vs_pyscf_${case} rys_fock_${case})
  endforeach()
  foreach(case water_631gs water_ccpvdz formamide_ccpvdz)
    foreach(path hgp nosym rys)
      terco_add_golden(check_erf_fock_vs_pyscf_${case}_${path} erf_fock_${case}_${path})
    endforeach()
  endforeach()
endif()
