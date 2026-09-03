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
set(TERCO_TEST_DEFS TRC_LMAX=${TERCO_LMAX} TRC_FOCK6 TRC_PERCLASS)

add_library(terco_test_basis STATIC trc_test_basis.F90)
target_compile_definitions(terco_test_basis PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(terco_test_basis PUBLIC ${main_lib})
set_target_properties(terco_test_basis PROPERTIES Fortran_MODULE_DIRECTORY
                                                  ${CMAKE_BINARY_DIR}/modules)

#
# The library's definitions are PRIVATE -- TRC_FOCK6 and TRC_PERCLASS are
# implementation choices and no consumer should need them. A test can, though:
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
  cmake_parse_arguments(ARG "NEEDS_BASIS" "" "ARGS" ${ARGN})
  add_executable(${name} ${source})
  target_compile_definitions(${name} PRIVATE ${TERCO_TEST_DEFS})
  target_link_libraries(${name} PRIVATE ${main_lib})
  if(ARG_NEEDS_BASIS)
    target_link_libraries(${name} PRIVATE terco_test_basis)
  endif()
  target_include_directories(${name} PRIVATE ${CMAKE_BINARY_DIR}/modules)
  add_test(NAME terco/${name} COMMAND ${name} ${ARG_ARGS})
  # The oracles are read by relative path, so the tests run from test/.
  set_tests_properties(
    terco/${name} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                             ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")
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
# The translated XC kernels against libxc, through the values `xc_ref.py`
# recorded in reference_xc.dat. Every derivative to second order, both spin
# cases, eleven decades of density. Needs no GPU and no pyscf to run; pyscf
# only to regenerate the reference.
#
terco_add_check(check_xc check_xc.F90)

# The quadrature against closed forms: Lebedev moments, the partition of
# unity, an exact radial Gaussian, and sums of Gaussians and Slaters over
# water integrating to the atom count.
terco_add_check(check_grid check_grid.F90 NEEDS_BASIS)

#
# Basis functions on grid points: the numerical overlap against trc_1e's
# analytic one, and the batched, screened evaluation against the unscreened
# one. It also writes ao_probe.bin for the pyscf comparison registered
# further down.
#
terco_add_check(check_ao check_ao.F90)

#
# The integrator from a fixed density: the electron count against Tr(DS),
# independence of the chunking, symmetry of V_xc -- and xc_probe.bin for
# the same-grid pyscf comparison of E_xc and V_xc registered further down.
#
terco_add_check(check_xc_energy check_xc_energy.F90)

# Four-centre integrals and the Fock digestion, against libfint's own numbers
# recorded in reference.dat.
terco_add_check(check_binfock check_binfock.F90)

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
# under TRC_PERCLASS only the non-symmetric build reaches the shared kernel,
# so every other check passes without executing a line of it.
#
terco_add_check(check_nosym check_nosym.F90 NEEDS_BASIS ARGS water.xyz)

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

#
# Cartesian multipoles through the octupole. The Fortran half writes a binary
# dump; `mult_ref.py` compares all 39 components against pyscf and refuses to
# judge any of them until the overlap agrees first.
#
terco_add_check(check_mult check_mult.F90)

#
# The pyscf comparison is registered only if pyscf is actually importable.
#
# terco does not depend on pyscf and neither should its test suite: a missing
# optional oracle should not turn the suite red, because a red suite that is
# always red stops being read. `check_mult` above already validates the
# multipole kernels against libfint through `check_api`; this adds an
# independent second opinion when one is available.
#
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy, pyscf"
    RESULT_VARIABLE _pyscf_rc
    OUTPUT_QUIET ERROR_QUIET)
  if(_pyscf_rc EQUAL 0)
    add_test(NAME terco/check_mult_vs_pyscf
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/mult_ref.py)
    set_tests_properties(
      terco/check_mult_vs_pyscf
      PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS
                 terco/check_mult)
    # Basis function values and gradients, elementwise against eval_ao.
    add_test(NAME terco/check_ao_vs_pyscf
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ao_ref.py)
    set_tests_properties(
      terco/check_ao_vs_pyscf
      PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS
                 terco/check_ao)
    # E_xc and V_xc on the same grid and density, six functionals.
    add_test(NAME terco/check_xc_energy_vs_pyscf
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/xc_energy_ref.py)
    set_tests_properties(
      terco/check_xc_energy_vs_pyscf
      PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS
                 terco/check_xc_energy)
  else()
    message(STATUS "terco: pyscf not importable, skipping the multipole cross-check")
  endif()
endif()

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
target_link_libraries(scf_rks PRIVATE ${main_lib} terco_test_basis LAPACK::LAPACK)
target_include_directories(scf_rks PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/scf_rks COMMAND scf_rks)
set_tests_properties(
  terco/scf_rks PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                           ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")
if(Python3_FOUND AND _pyscf_rc EQUAL 0)
  add_test(NAME terco/scf_rks_vs_pyscf
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rks_ref.py)
  set_tests_properties(
    terco/scf_rks_vs_pyscf PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                      DEPENDS terco/scf_rks)
endif()

#
# The SCF over ranks. Built in every configuration -- pic-mpi's serial
# backend makes it rank 0 of 1 -- and registered on one, two and four ranks
# when MPI is on, each run asserting the single-rank energies. No
# CUDA_VISIBLE_DEVICES here: the point is one card per rank.
#
add_executable(scf_mpi scf_mpi.F90)
target_compile_definitions(scf_mpi PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_mpi PRIVATE ${main_lib} terco_test_basis)
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
else()
  add_test(NAME terco/scf_mpi COMMAND scf_mpi)
  set_tests_properties(terco/scf_mpi PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                                ENVIRONMENT "CUDA_VISIBLE_DEVICES=0")
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
foreach(b bench_binwater bench_many bench_xc)
  target_include_directories(${b} PRIVATE ${CMAKE_BINARY_DIR}/modules)
endforeach()
