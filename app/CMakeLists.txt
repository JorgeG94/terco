#
# The `terco` command (main.F90): an SCF on an xyz file, from the shell.
#
# It drives the library through its context API, so it adds no physics of its
# own; what it adds is argument handling (pic_cli), an xyz reader and the
# lookup that turns `--basis def2-svp` into a file. The target is terco_cli
# because `terco` is already the library's; the executable is still `terco`.
#
include(GNUInstallDirs)

set(TERCO_BASIS_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/terco/basis_sets")

#
# --opt (DL-FIND), behind TERCO_ENABLE_DLFIND -- see that option's own
# comment in the root CMakeLists.txt for the licensing reason it defaults
# OFF, and src/opt/CMakeLists.txt for libterco_opt itself, which is where
# libdlfind is now fetched and linked (never into ${main_lib}, which is
# what keeps the LGPL relink obligation and DL-FIND's `error stop` out of
# the library).
#
# With the option ON, terco_cli gets module `trc_dlfind_bridge` from
# libterco_opt (linked below); with it OFF, from trc_dlfind_stub.F90,
# compiled straight into terco_cli, same public interface
# (trc_dlfind_optimize) either way, so main.F90 never branches on which
# build this is: the stub refuses --opt with a message instead of failing
# to link.
#
if(TERCO_ENABLE_DLFIND)
  set(_dlfind_bridge_src)
else()
  set(_dlfind_bridge_src trc_dlfind_stub.F90)
endif()

add_executable(
  terco_cli
  main.F90
  trc_app.F90
  trc_run.F90
  trc_run_report.F90
  trc_run_opt.F90
  trc_run_json.F90
  trc_run_input.F90
  ${_dlfind_bridge_src})
set_target_properties(
  terco_cli PROPERTIES OUTPUT_NAME terco Fortran_MODULE_DIRECTORY
                                         ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_compile_definitions(
  terco_cli
  PRIVATE TERCO_VERSION="${PROJECT_VERSION}"
          TERCO_BASIS_DIR="${CMAKE_INSTALL_PREFIX}/${TERCO_BASIS_INSTALL_DIR}"
          TERCO_SOURCE_BASIS_DIR="${PROJECT_SOURCE_DIR}/test/basis_sets")
# pic-mpi: the CLI initialises and finalises MPI, and asks for the world
# communicator's Fortran handle, itself (app/main.F90) -- the library only
# ever sees the handle trc_set_comm hands it.
# toml-f reads the TOML input files (app/trc_run_input.F90); only the CLI
# links it.
find_package("toml-f" REQUIRED)
target_link_libraries(terco_cli PRIVATE ${main_lib} pic::pic pic-mpi::pic-mpi toml-f::toml-f)
if(TERCO_ENABLE_DLFIND)
  target_compile_definitions(terco_cli PRIVATE TERCO_ENABLE_DLFIND)
  # terco::terco_opt (src/opt/CMakeLists.txt, added by the root
  # CMakeLists.txt before this directory) carries trc_dlfind_bridge's own
  # .mod, libdlfind's link requirement and libdlfind's own install() rule
  # -- terco_cli needs none of that spelled out again here.
  target_link_libraries(terco_cli PRIVATE terco::terco_opt)
endif()
unset(_dlfind_bridge_src)
# Find libterco next to the install's bin/, wherever the prefix is moved to.
# The library installs to lib/ (the root CMakeLists.txt), not LIBDIR.
if(APPLE)
  set_target_properties(terco_cli PROPERTIES INSTALL_RPATH "@loader_path/../lib")
else()
  set_target_properties(terco_cli PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
endif()

install(TARGETS terco_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# test/check_cli_grad (test/CMakeLists.txt, added before this directory by
# the root CMakeLists.txt) runs this executable as a subprocess, given its
# path through a generator expression -- but ctest itself does not know
# that check_cli_grad's own build needs terco_cli built first unless told,
# so that ordering is declared here, now that both targets exist.
if(TERCO_ENABLE_TESTING AND TARGET check_cli_grad)
  add_dependencies(check_cli_grad terco_cli)
endif()
if(TERCO_ENABLE_TESTING AND TARGET check_cli_grad_df)
  add_dependencies(check_cli_grad_df terco_cli)
endif()
if(TERCO_ENABLE_TESTING AND TARGET check_cli_lowercase)
  add_dependencies(check_cli_lowercase terco_cli)
endif()
if(TERCO_ENABLE_TESTING AND TARGET check_cli_hessian)
  add_dependencies(check_cli_hessian terco_cli)
endif()
if(TERCO_ENABLE_TESTING AND TARGET check_cli_opt)
  add_dependencies(check_cli_opt terco_cli)
endif()
if(TERCO_ENABLE_TESTING AND TARGET check_cli_opt_logging)
  add_dependencies(check_cli_opt_logging terco_cli)
endif()
if(TERCO_ENABLE_TESTING AND TARGET check_opt_vs_mqc)
  add_dependencies(check_opt_vs_mqc terco_cli)
endif()
foreach(_opt_pyscf_case water hcn)
  if(TERCO_ENABLE_TESTING AND TARGET check_opt_vs_pyscf_${_opt_pyscf_case})
    add_dependencies(check_opt_vs_pyscf_${_opt_pyscf_case} terco_cli)
  endif()
endforeach()
unset(_opt_pyscf_case)

# The basis sets the tests use, so `--basis def2-svp` works out of the box
# after an install. sap_helfem_large.json is not a basis but the fitted atomic
# potentials `--guess sap` reads; it is found the same way.
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/test/basis_sets/
  DESTINATION ${TERCO_BASIS_INSTALL_DIR}
  FILES_MATCHING
  PATTERN "*.json")

if(TERCO_ENABLE_TESTING)
  # From the source tree, not the install: TERCO_BASIS_PATH points at the
  # test basis sets, and each case checks an energy terco already pins
  # elsewhere.
  set(_env "TERCO_BASIS_PATH=${PROJECT_SOURCE_DIR}/test/basis_sets"
           "CUDA_VISIBLE_DEVICES=0")
  set(_xyz ${PROJECT_SOURCE_DIR}/test)

  # Water RHF/6-31G: scf_mpi's E_RHF (-75.983974465718), here from the SAD
  # guess rather than GWH, through a basis given by name.
  add_test(NAME terco/cli_rhf COMMAND terco_cli ${_xyz}/water.xyz --basis 6-31g)
  set_tests_properties(
    terco/cli_rhf PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION
                             "TOTAL ENERGY +-75\\.98397446")

  # A basis named with "*": 6-31G* resolves to the shipped 6-31gs.json (BSE's
  # own 6-31g_st_.json first, when a basis directory has it; BSE keeps "+"). Water RHF/6-31G*
  # is test/pyscf_refs/cphf_water_6-31gs.ref's E_SCF, -76.0105049954.
  add_test(NAME terco/cli_basis_star COMMAND terco_cli ${_xyz}/water.xyz --basis 6-31G*)
  set_tests_properties(
    terco/cli_basis_star PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION
                                    "TOTAL ENERGY +-76\\.010504")

  # The same through every other guess, so each one stays reachable from the
  # command line and converges to the same state.
  foreach(g sap gwh core minao)
    add_test(NAME terco/cli_guess_${g} COMMAND terco_cli ${_xyz}/water.xyz -b
                                                6-31g --guess ${g})
    set_tests_properties(
      terco/cli_guess_${g} PROPERTIES ENVIRONMENT "${_env}"
                                      PASS_REGULAR_EXPRESSION
                                      "TOTAL ENERGY +-75\\.98397446")
  endforeach()

  # --lindep reaches the orthogonalisation: at 0.3 water/6-31G loses two
  # functions and its energy rises (terco's own number for that basis).
  add_test(NAME terco/cli_lindep COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g
                                         --lindep 0.3)
  set_tests_properties(
    terco/cli_lindep
    PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION
               "2 function\\(s\\) dropped.*TOTAL ENERGY +-75\\.96606683")

  # A NaN threshold is refused, not passed through as "no dependence" (the
  # CLI's real parser already rejects it; trc_set_lindep refuses it too,
  # check_fukui).
  add_test(NAME terco/cli_lindep_nan COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g
                                             --lindep nan)
  set_tests_properties(
    terco/cli_lindep_nan
    PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION
               "lindep.*(not a real|must be in)")

  # Hirshfeld charges from outside the source tree: the free-atom table has
  # to be found the way a basis is, not relative to the working directory.
  # The O charge is terco's own (check_hirshfeld validates the method).
  add_test(NAME terco/cli_charges_hirshfeld
           COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g --charges hirshfeld)
  set_tests_properties(
    terco/cli_charges_hirshfeld
    PROPERTIES ENVIRONMENT "${_env}" WORKING_DIRECTORY
               ${CMAKE_CURRENT_BINARY_DIR} PASS_REGULAR_EXPRESSION
               "1  O +-0\\.342803")

  # CHELPG charges through the CLI, alone and beside Hirshfeld's, and an
  # unknown scheme refused before the SCF. The O charge is terco's own
  # (check_chelpg validates the method against pyscf and mqc).
  add_test(NAME terco/cli_charges_chelpg
           COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g --charges chelpg)
  set_tests_properties(
    terco/cli_charges_chelpg
    PROPERTIES ENVIRONMENT "${_env}" WORKING_DIRECTORY
               ${CMAKE_CURRENT_BINARY_DIR} PASS_REGULAR_EXPRESSION
               "chelpg charges.*1  O +-0\\.944615")
  add_test(NAME terco/cli_charges_hirshfeld_chelpg
           COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g --charges hirshfeld,chelpg)
  set_tests_properties(
    terco/cli_charges_hirshfeld_chelpg
    PROPERTIES ENVIRONMENT "${_env}" WORKING_DIRECTORY
               ${CMAKE_CURRENT_BINARY_DIR} PASS_REGULAR_EXPRESSION
               "hirshfeld charges.*1  O +-0\\.342803.*chelpg charges.*1  O +-0\\.944615")
  add_test(NAME terco/cli_charges_unknown
           COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g --charges mulliken)
  set_tests_properties(
    terco/cli_charges_unknown
    PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION
               "unknown --charges 'mulliken'")

  # --fukui without --guess runs from minao (the default it sets, which the
  # ions inherit; the NA-040 anion needs it), and f+, f- and f0 each sum to
  # one electron.
  add_test(NAME terco/cli_fukui_hirshfeld
           COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g --fukui hirshfeld)
  set_tests_properties(
    terco/cli_fukui_hirshfeld
    PROPERTIES ENVIRONMENT "${_env}" WORKING_DIRECTORY
               ${CMAKE_CURRENT_BINARY_DIR} PASS_REGULAR_EXPRESSION
               "guess +: minao.*sum +1\\.000000 +1\\.000000 +1\\.000000")

  # The CHELPG-condensed variant: each of N, N+1 and N-1 is fitted on its own
  # density (check_fukui validates the three charge sets and the indices
  # against pyscf). The O values are terco's own at water/6-31G; O's f+ is
  # negative, which Hirshfeld's (0.295041, above) never is.
  add_test(NAME terco/cli_fukui_chelpg
           COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g --fukui chelpg)
  set_tests_properties(
    terco/cli_fukui_chelpg
    PROPERTIES ENVIRONMENT "${_env}" WORKING_DIRECTORY
               ${CMAKE_CURRENT_BINARY_DIR} PASS_REGULAR_EXPRESSION
               "guess +: minao.*CHELPG-condensed Fukui.*1  O +-0\\.159473 +0\\.776391 +0\\.308459.*sum +1\\.000000 +1\\.000000 +1\\.000000")
  add_test(NAME terco/cli_fukui_unknown
           COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g --fukui mulliken)
  set_tests_properties(
    terco/cli_fukui_unknown
    PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION
               "unknown --fukui 'mulliken': use hirshfeld or chelpg")

  # Every convergence aid reaches the SCF and converges to the same state --
  # but NOT in the same number of iterations as the defaults, which is what
  # this test actually has to show: an energy match alone would still pass
  # if main.F90 silently dropped every one of these options on the floor,
  # since DIIS converges this system to the same state either way. Measured
  # with this binary: water/6-31G/SAD takes 9 iterations with no flags
  # (terco/cli_rhf) and 11 with --diis 20 --adiis --damp 0.5 --damp-steps 2
  # --diis-start 0.1 --level-shift 0.5, so pin the "SCF : N iterations" line
  # to 11, not just the energy.
  add_test(NAME terco/cli_scf_accel
           COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g --diis 20 --adiis --damp 0.5
                   --damp-steps 2 --diis-start 0.1 --level-shift 0.5)
  set_tests_properties(
    terco/cli_scf_accel
    PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION
               "SCF +: +11 iterations, converged.*TOTAL ENERGY +-75\\.98397446")
  # And a value outside its range is refused, not passed on.
  add_test(NAME terco/cli_scf_accel_refused
           COMMAND terco_cli ${_xyz}/water.xyz -b 6-31g --damp 0)
  set_tests_properties(
    terco/cli_scf_accel_refused PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION
                                           "--damp must be in")

  # A missing geometry is an error with a message, not a crash.
  # A double hybrid whose SCF stops unconverged never runs RI-MP2, so its
  # TOTAL ENERGY is the SCF energy, never the unset double-hybrid total (0).
  add_test(NAME terco/cli_dh_unconverged COMMAND terco_cli ${_xyz}/water.xyz -b
                                                 6-31g -f b2plyp --iterations 1)
  set_tests_properties(
    terco/cli_dh_unconverged
    PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION
               "TOTAL ENERGY +-7[0-9]\\." FAIL_REGULAR_EXPRESSION
               "TOTAL ENERGY +-?0\\.0+ ")

  # A meta-GGA's gradient runs through the CLI (its numbers are checked by
  # test/check_xc_grad_vs_pyscf and check_xc_grad_fd; this is the wiring).
  add_test(NAME terco/cli_grad_mgga COMMAND terco_cli ${_xyz}/water.xyz -b
                                            6-31g -f r2scan --grad)
  set_tests_properties(
    terco/cli_grad_mgga
    PROPERTIES ENVIRONMENT "${_env}" PASS_REGULAR_EXPRESSION "max \\|gradient\\| +: +[0-9]"
               FAIL_REGULAR_EXPRESSION "not implemented")

  add_test(NAME terco/cli_missing_file COMMAND terco_cli no_such_file.xyz)
  set_tests_properties(terco/cli_missing_file PROPERTIES WILL_FAIL TRUE)

  # --save-scf, then RI-MP2 alone from the file with --restart: the same
  # energies as the run that saved it, and a refusal for a run it does not
  # fit (another basis, method or charge).
  add_test(
    NAME terco/cli_restart
    COMMAND
      ${CMAKE_COMMAND} -DTERCO=$<TARGET_FILE:terco_cli> -DXYZ=${_xyz}/water.xyz
      -DWORK=${CMAKE_CURRENT_BINARY_DIR}/cli_restart -P
      ${CMAKE_CURRENT_SOURCE_DIR}/cli_restart_check.cmake)
  set_tests_properties(terco/cli_restart PROPERTIES ENVIRONMENT "${_env}"
                                                    PASS_REGULAR_EXPRESSION "restart: RI-MP2 identical")

  #
  # The CLI itself over ranks -- mirroring test/scf_mpi.F90's own launcher
  # logic (same reason for the compiler-wrapper's own mpiexec, not the
  # first one on PATH), since this is a different driver over the same
  # trc_set_comm/trc_run_scf and needs its own check: np 1, np 2 and np 3
  # must agree on water RHF, B3LYP and wB97X-V (RKS, VV10), and on NH2
  # radical wB97X-V (UKS, VV10) -- see cli_mpi_check.py's own header.
  # Registered only when this build can run it.
  #
  if(TERCO_ENABLE_MPI)
    find_package(Python3 COMPONENTS Interpreter)
    get_filename_component(_mpi_bindir "${CMAKE_Fortran_COMPILER}" DIRECTORY)
    if(EXISTS "${_mpi_bindir}/mpiexec")
      set(_cli_mpiexec "${_mpi_bindir}/mpiexec")
      set(_cli_npflag "-np")
    else()
      find_package(MPI REQUIRED COMPONENTS Fortran)
      set(_cli_mpiexec "${MPIEXEC_EXECUTABLE}")
      set(_cli_npflag "${MPIEXEC_NUMPROC_FLAG}")
    endif()
    if(Python3_FOUND)
      add_test(
        NAME terco/cli_mpi
        COMMAND
          ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cli_mpi_check.py
          ${_cli_mpiexec} ${_cli_npflag} $<TARGET_FILE:terco_cli> ${_xyz}/water.xyz
          ${PROJECT_SOURCE_DIR}/test/basis_sets ${_xyz}/nh2.xyz)
      set_tests_properties(terco/cli_mpi PROPERTIES WORKING_DIRECTORY
                                                     ${CMAKE_CURRENT_SOURCE_DIR}
                                                     TIMEOUT 300)
    endif()
    unset(_cli_mpiexec)
    unset(_cli_npflag)
  endif()
endif()
