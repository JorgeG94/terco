#
# The programs that are not test-drive suites: benchmarks, reference
# generators, demos, the MPI programs (a program with collective MPI calls
# keeps no `contains`, so these stay plain) and the SCF driver. Moved out of
# test/ so the top level holds only suites. The programs read their data files
# (xyz, basis_sets, references) relative to test/, so their tests keep test/
# as WORKING_DIRECTORY: TERCO_TEST_DIR, set by test/CMakeLists.txt, stands
# where CMAKE_CURRENT_SOURCE_DIR did. _mpiexec and _npflag come from the launcher
# detection there too. The LONG-label disable loop in test/CMakeLists.txt reads
# the TESTS property of test/ only, so a LONG-labelled test added here would
# need its own handling.
#

# A complete RHF through the public API (scf_rhf.F90's own header). Built, not
# registered: it prints an energy and exits on convergence only.
add_executable(scf_rhf scf_rhf.F90)
target_compile_definitions(scf_rhf PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_rhf PRIVATE ${main_lib} terco_test_basis LAPACK::LAPACK)
target_include_directories(scf_rhf PRIVATE ${CMAKE_BINARY_DIR}/modules)

# check_mdj's timing harness: not a test, a tool (bench_mdj <xyz> <basis> [thresh] [repeats]).
add_executable(bench_mdj bench_mdj.F90)
target_compile_definitions(bench_mdj PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_mdj PRIVATE ${main_lib} terco_test_basis)
target_include_directories(bench_mdj PRIVATE ${CMAKE_BINARY_DIR}/modules)

#
# The SCF over ranks. Built in every configuration -- pic-mpi's serial
# backend makes it rank 0 of 1 -- and registered on one, two and four ranks
# when MPI is on, each run asserting the single-rank energies. No
# CUDA_VISIBLE_DEVICES here: the point is one card per rank.
#
add_executable(scf_mpi scf_mpi.F90)
target_compile_definitions(scf_mpi PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(scf_mpi PRIVATE ${main_lib} terco_test_basis terco_c_interfaces)
target_include_directories(scf_mpi PRIVATE ${CMAKE_BINARY_DIR}/modules)
if(TERCO_ENABLE_MPI)
  foreach(np 1 2 4)
    add_test(NAME terco/scf_mpi_np${np}
             COMMAND ${_mpiexec} ${_npflag} ${np} $<TARGET_FILE:scf_mpi>)
    set_tests_properties(terco/scf_mpi_np${np} PROPERTIES WORKING_DIRECTORY ${TERCO_TEST_DIR})
  endforeach()
else()
  add_test(NAME terco/scf_mpi COMMAND scf_mpi)
  set_tests_properties(terco/scf_mpi PROPERTIES WORKING_DIRECTORY ${TERCO_TEST_DIR}
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
  set_tests_properties(terco/scf_mpi_frag PROPERTIES WORKING_DIRECTORY ${TERCO_TEST_DIR}
                                                     TIMEOUT 300)
endif()

# RI-CCSD(T) over ranks, one device each (rccsd_mpi.F90's header): the energies, E_CCSD and amplitudes (to round-off:
# the ladder is cut over the ranks) and the (T) sum against the single-rank run and pyscf, the same result on every rank, and with the argument "fail" the
# failure paths (a rank refused, a CCSD that does not converge, a (T) routine that differs, ...), which must end on all
# ranks. The TIMEOUT is the test of the second: a rank left in a collective call fails ctest, it does not wedge it.
#
add_executable(rccsd_mpi rccsd_mpi.F90)
target_compile_definitions(rccsd_mpi PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(rccsd_mpi PRIVATE ${main_lib} terco_rccsd_test_support terco_rccsd_ref_energies)
target_include_directories(rccsd_mpi PRIVATE ${CMAKE_BINARY_DIR}/modules)
if(TERCO_ENABLE_MPI)
  # One rank per device is the rule, and a CI box may have fewer GPUs than ranks: a GPU build registers np 1 and 2
  # only; a build without a GPU (ranks are host processes) also runs 3 and 4, which cut the ladder three and four ways.
  if(TERCO_ENABLE_GPU)
    set(_rccsd_mpi_nps 1 2)
  else()
    set(_rccsd_mpi_nps 1 2 3 4)
  endif()
  foreach(np ${_rccsd_mpi_nps})
    add_test(NAME terco/rccsd_mpi_np${np}
             COMMAND ${_mpiexec} ${_npflag} ${np} $<TARGET_FILE:rccsd_mpi>)
    set_tests_properties(terco/rccsd_mpi_np${np} PROPERTIES WORKING_DIRECTORY ${TERCO_TEST_DIR}
                                                           TIMEOUT 900)
  endforeach()
  if(TERCO_ENABLE_GPU)
    # the transfers through the host (trc_mpi_device: TERCO_MPI_DEVICE=0); the default run uses MPI device buffers
    # where the build says the MPI takes them
    add_test(NAME terco/rccsd_mpi_np2_hostbuf
             COMMAND ${_mpiexec} ${_npflag} 2 $<TARGET_FILE:rccsd_mpi>)
    set_tests_properties(terco/rccsd_mpi_np2_hostbuf PROPERTIES WORKING_DIRECTORY ${TERCO_TEST_DIR}
                                                               TIMEOUT 900 ENVIRONMENT "TERCO_MPI_DEVICE=0")
  endif()
  add_test(NAME terco/rccsd_mpi_fail_np2
           COMMAND ${_mpiexec} ${_npflag} 2 $<TARGET_FILE:rccsd_mpi> fail)
  set_tests_properties(terco/rccsd_mpi_fail_np2 PROPERTIES WORKING_DIRECTORY ${TERCO_TEST_DIR}
                                                           TIMEOUT 600)
endif()

# The makeEFP pathway on top of a complete RHF (demo_makefp.F90's own header).
add_executable(demo_makefp demo_makefp.F90)
target_compile_definitions(demo_makefp PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(demo_makefp PRIVATE ${main_lib} terco_test_basis LAPACK::LAPACK)
target_include_directories(demo_makefp PRIVATE ${CMAKE_BINARY_DIR}/modules)
add_test(NAME terco/demo_makefp COMMAND demo_makefp water.xyz)
set_tests_properties(
  terco/demo_makefp PROPERTIES WORKING_DIRECTORY ${TERCO_TEST_DIR}
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
# Warm-start iteration counts, density vs orbital-projection -- see
# bench_warmstart.F90's own header. Not a pass/fail, like its siblings.
add_executable(bench_warmstart bench_warmstart.F90)
target_compile_definitions(bench_warmstart PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_warmstart PRIVATE ${main_lib} terco_test_basis)
add_executable(bench_schwarz bench_schwarz.F90)
target_compile_definitions(bench_schwarz PRIVATE ${TERCO_TEST_DEFS})
target_link_libraries(bench_schwarz PRIVATE ${main_lib} terco_test_basis)
foreach(b bench_binwater bench_many bench_xc bench_gc bench_sph bench_grad_1e bench_warmstart bench_schwarz)
  target_include_directories(${b} PRIVATE ${CMAKE_BINARY_DIR}/modules)
endforeach()

# --opt's own warm-start measurement (bench_opt_warmstart.F90's own
# header): SCF iterations per step, mqc's own baseline vs the density and
# projected warm starts. Both benches link terco::terco_opt (src/opt/
# CMakeLists.txt) for module `trc_dlfind_bridge` now, rather than
# recompiling its source here -- that library is exactly what terco_cli
# itself links. TERCO_ENABLE_DLFIND only: the stub build has no optimiser
# to measure (and the library target does not exist).
if(TERCO_ENABLE_DLFIND)
  add_executable(bench_opt_warmstart bench_opt_warmstart.F90)
  target_compile_definitions(bench_opt_warmstart PRIVATE ${TERCO_TEST_DEFS})
  target_link_libraries(bench_opt_warmstart PRIVATE ${main_lib} terco_test_basis terco::terco_opt)
  target_include_directories(bench_opt_warmstart PRIVATE ${CMAKE_BINARY_DIR}/modules)
  set_target_properties(bench_opt_warmstart PROPERTIES Fortran_MODULE_DIRECTORY
                                                        ${CMAKE_CURRENT_BINARY_DIR}/modules_opt_bench)

  # The per-image-vs-shared-density NEB warm start measurement
  # (bench_opt_neb_warmstart.F90's own header) -- the --opt-neb analogue
  # of bench_opt_warmstart above.
  add_executable(bench_opt_neb_warmstart bench_opt_neb_warmstart.F90)
  target_compile_definitions(bench_opt_neb_warmstart PRIVATE ${TERCO_TEST_DEFS})
  target_link_libraries(bench_opt_neb_warmstart PRIVATE ${main_lib} terco_test_basis terco::terco_opt)
  target_include_directories(bench_opt_neb_warmstart PRIVATE ${CMAKE_BINARY_DIR}/modules)
  set_target_properties(bench_opt_neb_warmstart PROPERTIES Fortran_MODULE_DIRECTORY
                                                            ${CMAKE_CURRENT_BINARY_DIR}/modules_opt_neb_bench)
endif()

