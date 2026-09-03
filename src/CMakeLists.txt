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

# --- the molecular quadrature -------------------------------------------------
#
# Lebedev spheres, the Treutler-Ahlrichs radial mesh, NWChem-style pruning,
# Becke and Stratmann partitions, and the assembled grid. Lifted from
# metalquicha's src/methods with the module prefix renamed and nothing else,
# plus the minimal error type they report through.
#
target_sources(
  ${main_lib}
  PRIVATE grid/trc_error.F90
          grid/trc_lebedev_data.F90
          grid/trc_lebedev.F90
          grid/trc_dft_radial_data.F90
          grid/trc_dft_radial.F90
          grid/trc_dft_partition.F90
          grid/trc_dft_prune.F90
          grid/trc_dft_grid.F90)

# --- exchange-correlation ---------------------------------------------------
#
# The XC functionals. One hand-written helper module and, per functional
# kernel, one module TRANSLATED from ExchCXX's Maple-generated C++ by
# tools/exchcxx_to_fortran.py. Translated rather than generated: the source
# is another project's checkout, so the output is committed and the tool is
# rerun by hand when a kernel is added. See src/xc/trc_xc_util.F90.
#
target_sources(
  ${main_lib}
  PRIVATE xc/bit_repro.f90          # bit-reproducible exp/log/pow/erfc/cbrt, vendored
          xc/bit_repro_helpers.f90  # its submodule
          xc/trc_xc_util.F90
          xc/trc_xc_slater_exchange.F90
          xc/trc_xc_vwn.F90
          xc/trc_xc_vwn_rpa.F90
          xc/trc_xc_pbe_x.F90
          xc/trc_xc_pbe_c.F90
          xc/trc_xc_b88.F90
          xc/trc_xc_lyp.F90
          xc/trc_xc_functional.F90  # compositions, and the per-point dispatch
          xc/trc_collocation.F90    # basis functions and gradients at a point
          xc/trc_xc_batch.F90       # the grid, batched and screened
          xc/trc_xc.F90)            # the integrator: E_xc and V_xc

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
  # afford the same ceiling as the four-centre kernels for far less compile.
  #
  # They follow TERCO_LMAX rather than a literal, because a one-electron path
  # that stops below the four-centre one is not a limitation -- it is a wrong
  # Fock matrix. H would be built from whatever one_e_dispatch returns for a
  # key it has no case for, and an SCF converges on that without complaint.
  execute_process(COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py --lmax ${TERCO_LMAX}
                          -o ${_gen_dir}/trc_1e_kernels.F90 RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_1e.py failed (${_rc})")
  endif()
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py --lmax ${TERCO_LMAX} --multipole
            ${_gen_dir}/trc_mult_kernels.F90 RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_1e.py --multipole failed (${_rc})")
  endif()
  execute_process(COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_df.py --lmax ${TERCO_LMAX}
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
