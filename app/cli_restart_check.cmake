# --save-scf then --restart: RI-MP2 from the saved SCF must be the RI-MP2 of
# the run that saved it, to the last printed digit, and a restart into a
# run it does not fit must be refused rather than used.
#
#   cmake -DTERCO=<terco> -DXYZ=<water.xyz> -DWORK=<scratch dir> -P cli_restart_check.cmake
#
# TERCO_BASIS_PATH comes from the test's environment.

file(MAKE_DIRECTORY ${WORK})
set(_scf ${WORK}/water.scf)
file(REMOVE ${_scf} ${WORK}/full.json ${WORK}/restart.json)

function(terco_run out_var rc_var)
  execute_process(COMMAND ${TERCO} ${XYZ} ${ARGN} WORKING_DIRECTORY ${WORK}
                  OUTPUT_VARIABLE _out ERROR_VARIABLE _err RESULT_VARIABLE _rc)
  set(${out_var} "${_out}${_err}" PARENT_SCOPE)
  set(${rc_var} ${_rc} PARENT_SCOPE)
endfunction()

function(json_number path key var)
  file(READ ${path} _j)
  string(JSON _v GET "${_j}" "${key}")
  set(${var} ${_v} PARENT_SCOPE)
endfunction()

# The full run, saving its SCF.
terco_run(out rc -b 6-31g --rimp2 --save-scf ${_scf} --json ${WORK}/full.json)
if(NOT rc EQUAL 0 OR NOT EXISTS ${_scf})
  message(FATAL_ERROR "the full run failed (${rc}):\n${out}")
endif()

# RI-MP2 alone, from the file.
terco_run(out rc -b 6-31g --rimp2 --restart ${_scf} --json ${WORK}/restart.json)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "the restarted run failed (${rc}):\n${out}")
endif()
if(NOT out MATCHES "read from")
  message(FATAL_ERROR "the restarted run does not say it read the SCF:\n${out}")
endif()
if(out MATCHES "fock\\(s\\)")
  message(FATAL_ERROR "the restarted run ran SCF iterations:\n${out}")
endif()
foreach(key "scf energy" "mp2 opposite-spin energy" "mp2 same-spin energy"
            "two-electron energy" "iterations")
  json_number(${WORK}/full.json "${key}" a)
  json_number(${WORK}/restart.json "${key}" b)
  if(NOT a STREQUAL b)
    message(FATAL_ERROR "${key}: ${a} from the full run, ${b} from the restart")
  endif()
endforeach()

# Runs the file does not fit, each differing from the saved one in one thing
# only. Each must fail, for that reason.
foreach(case "different basis;-b;cc-pvdz" "made with;-b;6-31g;-f;b3lyp" "charge 0;-b;6-31g;-c;2")
  list(POP_FRONT case why)
  terco_run(out rc ${case} --rimp2 --restart ${_scf})
  if(rc EQUAL 0 OR NOT out MATCHES "does not fit this run: [^\n]*${why}")
    message(FATAL_ERROR "a restart with ${case} was not refused for '${why}' (${rc}):\n${out}")
  endif()
endforeach()

message(STATUS "restart: RI-MP2 identical, three mismatches refused")
