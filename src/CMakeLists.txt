#
# The library: a hand-written core, and kernels that are generated.
#
# Roughly nine tenths of the compiled lines here come out of `scripts/`. That
# is the design and not an accident -- see docs/GENERATORS.md -- so this file
# is mostly about deciding WHICH generated set to use.
#

# --- hand written -----------------------------------------------------------
target_sources(
  ${main_lib}
  PRIVATE trc_boys.F90       # Boys function, tabulated and asymptotic
          trc_tables.F90     # Cartesian index tables
          trc_cart.F90
          trc_batch.F90      # primitive-pair batching
          trc_hgp.F90        # Head-Gordon-Pople driver
          trc_screen.F90     # Schwarz bounds
          trc_bins.F90       # the binned shell-pair container
          trc_binkernel.F90  # the shared quartet kernel and the digestion
          trc_api.F90        # basis, pair list, 1e, multipoles, DF
          trc_fock.F90       # the four-centre Fock builds
          trc_capi.F90       # the bind(c) skin
          trc_acc_shim.F90)  # see the file: this one is not optional

# `trc_binkernel.F90` textually includes the quartet body rather than calling
# it. NVHPC will not inline an `!$acc routine seq` across a module boundary --
# `v` has to be materialised in memory to cross the call and 35 kB of spill
# appears -- so the recurrence is pasted in and the compiler sees all of it.
target_include_directories(${main_lib} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# --- generated --------------------------------------------------------------
#
# The repository ships ONE generated configuration: TERCO_LMAX=2, in
# src/generated. That is what `fpm build` uses, since fpm cannot run a
# generator, and it is what a consumer gets without needing Python.
#
# Any other ceiling -- or an explicit TERCO_REGENERATE -- runs the generators.
# At configure time rather than build time, because the file list is not known
# until they have run and CMake needs it to construct the target.
#
set(_gen_dir "${CMAKE_CURRENT_SOURCE_DIR}/generated")
set(_regen OFF)
if(NOT TERCO_LMAX STREQUAL "2")
  set(_gen_dir "${CMAKE_BINARY_DIR}/generated")
  set(_regen ON)
elseif(TERCO_REGENERATE)
  set(_regen ON)
endif()

if(_regen)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  file(MAKE_DIRECTORY "${_gen_dir}")
  set(_scripts "${CMAKE_SOURCE_DIR}/scripts")

  message(STATUS "terco: generating kernels for LMAX=${TERCO_LMAX} in ${_gen_dir}")

  # Four-centre, one module per (la,lb,lc,ld) class. `--split` matters: at
  # LMAX=2 a single module is a gigabyte of NVVM that does not finish.
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_perclass.py --lmax
            ${TERCO_LMAX} --split ${_gen_dir}
    RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_perclass.py failed (${_rc})")
  endif()

  # One-electron, multipoles and density fitting carry angular momentum on two
  # or three centres rather than four, so they grow as L^2 or L^3 and can
  # afford a higher ceiling than the four-centre kernels. Fixed at 2 here
  # because that is what the shipped set was built with; raising it is a
  # one-line change and a longer compile.
  execute_process(COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py --lmax 2
                          -o ${_gen_dir}/trc_1e_kernels.F90 RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_1e.py failed (${_rc})")
  endif()
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py --lmax 2 --multipole
            ${_gen_dir}/trc_mult_kernels.F90 RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_1e.py --multipole failed (${_rc})")
  endif()
  execute_process(COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_df.py --lmax 2
                          -o ${_gen_dir}/trc_df_kernels.F90 RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_df.py failed (${_rc})")
  endif()

  # The unrolled recurrences are textual includes, and only the shared kernel
  # uses them -- at LMAX<=1, for the reason the top-level CMakeLists gives.
  # They carry a `#if TRC_LMAX != n / #error` guard, so a stale one is a
  # compile error rather than a wrong number. That guard exists because a
  # stale one WAS a wrong number, silently, for a month.
  if(TERCO_LMAX LESS_EQUAL 1)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_vrr.py --lmax ${TERCO_LMAX}
              --mode include -o ${_gen_dir}/trc_vrr_cases.inc OUTPUT_QUIET)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_hrr.py --lmax ${TERCO_LMAX}
              -o ${_gen_dir}/trc_hrr_cases.inc OUTPUT_QUIET)
    target_include_directories(${main_lib} PRIVATE ${_gen_dir})
  endif()
endif()

file(GLOB _generated CONFIGURE_DEPENDS "${_gen_dir}/*.F90")
if(NOT _generated)
  message(
    FATAL_ERROR
      "No generated kernels in ${_gen_dir}. Configure with -DTERCO_REGENERATE=ON "
      "(needs Python 3), or check out the committed set.")
endif()
list(LENGTH _generated _ngen)
message(STATUS "terco: ${_ngen} generated kernel modules from ${_gen_dir}")
target_sources(${main_lib} PRIVATE ${_generated})
