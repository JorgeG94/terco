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
  PRIVATE core/trc_boys.F90       # Boys function, tabulated and asymptotic
          core/trc_tables.F90     # Cartesian index tables
          core/trc_cart.F90
          eri/trc_batch.F90      # primitive-pair batching
          eri/trc_hgp.F90        # Head-Gordon-Pople driver
          eri/trc_screen.F90     # Schwarz bounds
          eri/trc_bins.F90       # the binned shell-pair container
          eri/trc_binkernel.F90  # the shared quartet kernel and the digestion
          eri/trc_grad.F90       # two-electron gradient launch plan, Rys only
          api/trc_api.F90        # basis, pair list, 1e, multipoles, DF
          basis/trc_decontract.F90 # splitting a general contraction into primitives
          eri/trc_mdj.F90       # the McMurchie-Davidson J engine (no quartets)
          eri/trc_mdg.F90       # MD exchange through batched GEMMs, for the classes TRC_MDG_CLASSES names
          eri/trc_eri.F90       # the four-centre Fock builds
          core/trc_linalg.F90     # dense linear algebra on device-resident matrices
          core/trc_si.F90         # a quantity printed with its SI prefix (2.45 PFLOP)
          core/trc_step_table.F90 # shared timed-step table (RI-MP2, DF-JK): time/FLOP/B/s, NVTX
          core/trc_mpi_device.F90 # whether MPI gets device buffers directly (CUDA-aware)
          core/trc_constants.F90  # the single source of truth for math/CODATA constants and atomic masses
          solvation/trc_pcm_radii.F90  # cavity radii for the polarisable continuum
          solvation/trc_pcm.F90        # the polarisable continuum: cavity, charges, operator
          scf/trc_scf_settings.F90      # trc_scf_options_t/result_t and the options a run resolves them to
          scf/trc_scf_context.F90       # everything an SCF run builds once per geometry
          scf/trc_scf_state.F90         # the SCF's per-iteration arrays, and the one enter/exit pair for them
          scf/trc_scf_fock.F90          # the Fock build: two-electron, XC, PCM, energy, commutator, folds
          scf/trc_scf_df.F90            # density-fitted J/K (--df): metric, Cholesky, (mn|P), B, per-iteration J/K
          scf/trc_scf_diis.F90          # CDIIS and ADIIS, the level shift and the damped step
          scf/trc_scf_driver.F90        # the SCF: HF and Kohn-Sham, R and U, resident
          scf/trc_scf_file.F90          # a converged SCF to a file and back, for a restart
          scf/trc_scf_grad_df.F90       # the density-fitted (--df) two-electron gradient: fitted J and K
          scf/trc_cphf.F90               # coupled-perturbed HF: the orbital-response solve and the static polarizability
          scf/trc_cpks.F90               # coupled-perturbed KS: trc_cphf plus the XC kernel response, RKS LDA/GGA
          scf/trc_scf_grad.F90          # the full analytic gradient: 1e + eri%gradient + nuclear repulsion
          scf/trc_vib.F90                # harmonic vibrational analysis of a Cartesian Hessian, host only
          scf/trc_thermo.F90             # ideal-gas RRHO thermochemistry from trc_vib's frequencies
          post_hf/rimp2/trc_rimp2_driver.F90      # RI-MP2 on the RHF reference, resident
          basis/trc_basis_json.F90 # a basis from a MolSSI BSE JSON file
          scf/trc_guess_json.F90 # a SAD guess from mqc's per-element JSON blocks
          scf/trc_sad.F90 # the SAD guess built here, one atomic SCF per element
          scf/trc_sap.F90 # SAP: free-atom potentials, superposed on a quadrature
          api/trc_capi.F90       # the bind(c) skin
          core/trc_cart_perm.F90  # AO-index permutation, lexicographic <-> this build's order
          core/trc_ao_norm.F90    # AO-index normalisation, this build's <-> GAMESS's per-component convention
          core/trc_sph.F90        # the spherical (pure) AO space: block-sparse Cartesian <-> spherical transforms
          core/trc_acc_shim.F90)  # see the file: this one is not optional

# --- the molecular quadrature -------------------------------------------------
#
# Lebedev spheres, the Treutler-Ahlrichs radial mesh, NWChem-style pruning,
# Becke and Stratmann partitions, and the assembled grid. Lifted from
# metalquicha's src/methods with the module prefix renamed and nothing else,
# plus the minimal error type they report through.
#
target_sources(
  ${main_lib}
  PRIVATE core/trc_error.F90
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
# trc_xc_functional.F90 (the kernel ids, the per-point dispatch, and the
# name -> composition registry) is GENERATED from the same ExchCXX checkout
# by tools/gen_xc_registry.py, likewise rerun by hand; see
# docs/functionals.md for what it produced.
#
target_sources(
  ${main_lib}
  PRIVATE xc/vendor/bit_repro.f90          # bit-reproducible exp/log/pow/erfc/cbrt, vendored
          xc/vendor/bit_repro_helpers.f90  # its submodule
          xc/trc_xc_util.F90
          xc/functionals/trc_xc_slater_exchange.F90
          xc/functionals/trc_xc_vwn.F90
          xc/functionals/trc_xc_vwn_rpa.F90
          xc/functionals/trc_xc_pbe_x.F90
          xc/functionals/trc_xc_pbe_c.F90
          xc/functionals/trc_xc_b88.F90
          xc/functionals/trc_xc_lyp.F90
          xc/functionals/trc_xc_wb97x_v_xc.F90
          xc/functionals/trc_xc_wb97x_xc.F90
          xc/functionals/trc_xc_wb97_xc.F90
          xc/functionals/trc_xc_vwn3.F90
          xc/functionals/trc_xc_pz81.F90
          xc/functionals/trc_xc_pz81_mod.F90
          xc/functionals/trc_xc_pw91_lda.F90
          xc/functionals/trc_xc_pw91_lda_mod.F90
          xc/functionals/trc_xc_pw91_lda_rpa.F90
          xc/functionals/trc_xc_rev_pbe_x.F90
          xc/functionals/trc_xc_rpbe_x.F90
          xc/functionals/trc_xc_optx_x.F90
          xc/functionals/trc_xc_mpw91_x.F90
          xc/functionals/trc_xc_pw86_x.F90
          xc/functionals/trc_xc_pw91_x.F90
          xc/functionals/trc_xc_p86_c.F90
          xc/functionals/trc_xc_p86vwn_ft_c.F90
          xc/functionals/trc_xc_pw91_c.F90
          xc/functionals/trc_xc_pbe_sol_c.F90
          xc/functionals/trc_xc_bmk_c.F90
          xc/functionals/trc_xc_n12_c.F90
          xc/functionals/trc_xc_n12_sx_c.F90
          xc/functionals/trc_xc_sogga11_x_c.F90
          xc/functionals/trc_xc_sogga11_x_x.F90
          xc/functionals/trc_xc_b97_d.F90
          xc/functionals/trc_xc_ityh_x.F90
          xc/functionals/trc_xc_ityh_x_015.F90
          xc/functionals/trc_xc_ityh_x_033.F90
          xc/functionals/trc_xc_hjs_pbe_x.F90
          xc/functionals/trc_xc_lrcwpbe_hjs_pbe_x.F90
          xc/functionals/trc_xc_lrcwpbeh_hjs_pbe_x.F90
          xc/functionals/trc_xc_wb97x_d_xc.F90
          xc/functionals/trc_xc_wb97x_d3_xc.F90
          # meta-GGA kernels: exc and exc_vxc only (see tools/exchcxx_to_fortran.py)
          xc/functionals/trc_xc_bc95_c.F90
          xc/functionals/trc_xc_bmk_x.F90
          xc/functionals/trc_xc_cf22d_c.F90
          xc/functionals/trc_xc_cf22d_x.F90
          xc/functionals/trc_xc_m05_2x_c.F90
          xc/functionals/trc_xc_m05_2x_x.F90
          xc/functionals/trc_xc_m05_c.F90
          xc/functionals/trc_xc_m05_x.F90
          xc/functionals/trc_xc_m06_2x_c.F90
          xc/functionals/trc_xc_m06_2x_x.F90
          xc/functionals/trc_xc_m06_c.F90
          xc/functionals/trc_xc_m06_hf_c.F90
          xc/functionals/trc_xc_m06_hf_x.F90
          xc/functionals/trc_xc_m06_l_c.F90
          xc/functionals/trc_xc_m06_l_x.F90
          xc/functionals/trc_xc_m06_sx_c.F90
          xc/functionals/trc_xc_m06_sx_x.F90
          xc/functionals/trc_xc_m06_x.F90
          xc/functionals/trc_xc_m08_hx_c.F90
          xc/functionals/trc_xc_m08_hx_x.F90
          xc/functionals/trc_xc_m08_so_c.F90
          xc/functionals/trc_xc_m08_so_x.F90
          xc/functionals/trc_xc_m11_c.F90
          xc/functionals/trc_xc_m11_x.F90
          xc/functionals/trc_xc_mbeef_x.F90
          xc/functionals/trc_xc_mn12_l_c.F90
          xc/functionals/trc_xc_mn12_l_x.F90
          xc/functionals/trc_xc_mn12_sx_c.F90
          xc/functionals/trc_xc_mn12_sx_x.F90
          xc/functionals/trc_xc_mn15_c.F90
          xc/functionals/trc_xc_mn15_l_c.F90
          xc/functionals/trc_xc_mn15_l_x.F90
          xc/functionals/trc_xc_mn15_x.F90
          xc/functionals/trc_xc_pkzb_c.F90
          xc/functionals/trc_xc_pkzb_x.F90
          xc/functionals/trc_xc_r2scan_c.F90
          xc/functionals/trc_xc_r2scan_x.F90
          xc/functionals/trc_xc_rev_m06_l_c.F90
          xc/functionals/trc_xc_rev_m06_l_x.F90
          xc/functionals/trc_xc_revtpss_c.F90
          xc/functionals/trc_xc_revtpss_x.F90
          xc/functionals/trc_xc_rscan_c.F90
          xc/functionals/trc_xc_rscan_x.F90
          xc/functionals/trc_xc_scan_c.F90
          xc/functionals/trc_xc_scan_x.F90
          xc/functionals/trc_xc_tpss_c.F90
          xc/functionals/trc_xc_tpss_x.F90
          xc/functionals/trc_xc_wb97m_v_xc.F90  # hand-written, not translated: see its own header
          xc/trc_xc_functional.F90  # compositions, and the per-point dispatch
          xc/trc_collocation.F90    # basis functions and gradients at a point
          xc/trc_xc_batch.F90       # the grid, batched and screened
          xc/trc_xc.F90             # the integrator: E_xc and V_xc
          xc/trc_xc_grad.F90        # the exchange-correlation nuclear gradient, LDA/GGA
          xc/trc_vv10.F90           # VV10 nonlocal correlation: E_nlc and its potential
          xc/trc_vv10_grad.F90)     # VV10's nuclear gradient, with grid response

# ifx's default -fp-model=fast reassociates and contracts the translated
# kernels' Maple power series, which were written assuming IEEE-strict
# evaluation order. Measured on the host (ifx 2026.0, check_xc): with the
# default flags, several kernels' own four evaluation routines
# (eval_exc/eval_exc_vxc/eval_fxc/eval_vxc_fxc) disagree with EACH OTHER
# past check_xc's 1e-12 "self" tolerance (p86vwn_ft_c, bmk_c, n12_c,
# n12_sx_c, sogga11_x_c, pw91_c) even though every one of them is close to
# libxc -- a sign of reassociated, path-dependent rounding, not a wrong
# formula. wb97_xc and wb97x_d_xc go further: at small s (the reduced
# density gradient, where their B97 power series has heavy cancellation),
# -ffast-math-style contraction loses enough precision that vxc/fxc are off
# by 1e-1 to 1 relative, not merely past a tight tolerance. Compiling these
# translated kernels with -fp-model=precise (strict IEEE evaluation order,
# no reassociation or contraction) fixed every one of them at check_xc's
# UNCHANGED tolerances; gfortran and nvfortran need nothing extra because
# neither reassociates by default.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set_source_files_properties(
    xc/functionals/trc_xc_slater_exchange.F90
    xc/functionals/trc_xc_vwn.F90
    xc/functionals/trc_xc_vwn_rpa.F90
    xc/functionals/trc_xc_pbe_x.F90
    xc/functionals/trc_xc_pbe_c.F90
    xc/functionals/trc_xc_b88.F90
    xc/functionals/trc_xc_lyp.F90
    xc/functionals/trc_xc_wb97x_v_xc.F90
    xc/functionals/trc_xc_wb97x_xc.F90
    xc/functionals/trc_xc_wb97_xc.F90
    xc/functionals/trc_xc_vwn3.F90
    xc/functionals/trc_xc_pz81.F90
    xc/functionals/trc_xc_pz81_mod.F90
    xc/functionals/trc_xc_pw91_lda.F90
    xc/functionals/trc_xc_pw91_lda_mod.F90
    xc/functionals/trc_xc_pw91_lda_rpa.F90
    xc/functionals/trc_xc_rev_pbe_x.F90
    xc/functionals/trc_xc_rpbe_x.F90
    xc/functionals/trc_xc_optx_x.F90
    xc/functionals/trc_xc_mpw91_x.F90
    xc/functionals/trc_xc_pw86_x.F90
    xc/functionals/trc_xc_pw91_x.F90
    xc/functionals/trc_xc_p86_c.F90
    xc/functionals/trc_xc_p86vwn_ft_c.F90
    xc/functionals/trc_xc_pw91_c.F90
    xc/functionals/trc_xc_pbe_sol_c.F90
    xc/functionals/trc_xc_bmk_c.F90
    xc/functionals/trc_xc_n12_c.F90
    xc/functionals/trc_xc_n12_sx_c.F90
    xc/functionals/trc_xc_sogga11_x_c.F90
    xc/functionals/trc_xc_sogga11_x_x.F90
    xc/functionals/trc_xc_b97_d.F90
    xc/functionals/trc_xc_ityh_x.F90
    xc/functionals/trc_xc_ityh_x_015.F90
    xc/functionals/trc_xc_ityh_x_033.F90
    xc/functionals/trc_xc_hjs_pbe_x.F90
    xc/functionals/trc_xc_lrcwpbe_hjs_pbe_x.F90
    xc/functionals/trc_xc_lrcwpbeh_hjs_pbe_x.F90
    xc/functionals/trc_xc_wb97x_d_xc.F90
    xc/functionals/trc_xc_wb97x_d3_xc.F90
    DIRECTORY ${PROJECT_SOURCE_DIR}
    PROPERTIES COMPILE_OPTIONS "-fp-model=precise")
endif()

# The translated kernels' Maple literals go in instructions, not constant bank
# 2. Every functional is reachable from the trc_xc_*_points kernels, so their
# literals and jump tables add up in those kernels' c[0x2] bank: 35808 of its
# 64 KiB on main, and over 64 KiB with the 47 meta-GGA kernels, where
# CUDA 12.9's nvlink silently mis-links (see TERCO_CHECK_DEVICE_LINK in the
# top-level CMakeLists.txt). With the literals inlined the largest bank is
# 1080 bytes (an ERI kernel), the energies are bit-identical, and the point
# kernel's time is unchanged within 1% (bench_xc gly10/cc-pVDZ level 3, V100:
# PBE here; PBE, B3LYP, wB97X-V, r2SCAN, TPSS and M06 on the meta-GGA branch).
# Device code grows by 0.5 MB of 33 MB (1.2 MB of 37 MB with the meta-GGAs).
if(TERCO_ENABLE_GPU)
  get_target_property(_trc_srcs ${main_lib} SOURCES)
  list(FILTER _trc_srcs INCLUDE REGEX "xc/functionals/")
  set_source_files_properties(
    ${_trc_srcs} DIRECTORY ${PROJECT_SOURCE_DIR}
    PROPERTIES COMPILE_OPTIONS "-Wptxas,--disable-optimizer-constants")
endif()

# `trc_binkernel.F90` textually includes the quartet body rather than calling
# it. NVHPC will not inline an `!$acc routine seq` across a module boundary --
# `v` has to be materialised in memory to cross the call and 35 kB of spill
# appears -- so the recurrence is pasted in and the compiler sees all of it.
target_include_directories(${main_lib} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# --- generated --------------------------------------------------------------
#
# The repository ships ONE generated configuration: TERCO_LMAX=2,
# TERCO_CART_ORDER=lexicographic, in src/generated. That is what `fpm build`
# uses, since fpm cannot run a generator, and it is what a consumer gets
# without needing Python.
#
# Any other ceiling, any other Cartesian order, or an explicit
# TERCO_REGENERATE runs the generators. At configure time rather than build
# time, because the file list is not known until they have run and CMake
# needs it to construct the target.
#
# "libcint" is accepted as an alias for TERCO_CART_ORDER, same as
# cartorder.py's own ORDERS dict -- it costs one line and nothing else needs
# to know the alias exists.
if(TERCO_CART_ORDER STREQUAL "libcint")
  set(TERCO_CART_ORDER "lexicographic")
endif()

# Generators write into a staging directory, and only files whose contents
# changed are copied to where they are compiled from. Written straight into
# place, every configure gave every generated file a new timestamp, and
# make recompiled all of them -- hours of nvfortran at TERCO_LMAX=3 for a
# reconfigure that changed one file. `prune` also removes files the
# generators no longer write; it is never set for src/generated, whose
# committed files are not all regenerated by every configuration.
function(terco_sync_generated stage dest prune)
  file(MAKE_DIRECTORY "${dest}")
  file(GLOB _new LIST_DIRECTORIES false RELATIVE "${stage}" "${stage}/*")
  foreach(_f ${_new})
    file(COPY_FILE "${stage}/${_f}" "${dest}/${_f}" ONLY_IF_DIFFERENT)
  endforeach()
  if(prune)
    file(GLOB _old LIST_DIRECTORIES false RELATIVE "${dest}" "${dest}/*")
    foreach(_f ${_old})
      if(NOT _f IN_LIST _new)
        file(REMOVE "${dest}/${_f}")
      endif()
    endforeach()
  endif()
endfunction()

# GENERATION STAMPS. The sync above stops an unchanged output from being
# RECOMPILED; these stop it from being regenerated at all. Every configure
# (`cmake ..`, or CMake's own re-run after a CMakeLists edit) used to run
# every generator below -- the Rys template, ~150 gen_rys.py/gen_rys_grad.py
# runs and both dispatchers, and with _regen the HGP, 1e, DF, Cartesian and
# spherical sets too -- 40-90 s of Python for a reconfigure that changed
# nothing.
#
# A block's generators are written once, as terco_gen calls inside a macro,
# and the macro is run twice: first with TERCO_GEN_MODE "plan", which only
# records each command line, then with "run", which executes them -- but
# only if the stamp says it must. The stamp is a SHA256 over
#   - every generator command line, exactly as it would run: the scripts'
#     paths, the interpreter (and its version), every option that reaches a
#     generator (TERCO_LMAX, TERCO_DF_LMAX, TERCO_DF_LMAX_AUX, TERCO_CART_ORDER,
#     TERCO_RYS_UNROLL_MAX, TERCO_DF_FORCE_RYS) and the class lists
#     computed below, so no option can be left out of it by forgetting to
#     list it; and
#   - the contents of every script those commands run or load (passed as
#     SCRIPTS: the scripts each block invokes plus what they import or load
#     with importlib -- none of them reads any other committed file).
# The stamp file also records the SHA256 of every file the block left in its
# destination. A configure skips the block only if the stamp matches AND
# that destination holds exactly those files with exactly those contents;
# a deleted, edited or stray output regenerates.
#
# Editing a generator re-runs the configure on the next build
# (CMAKE_CONFIGURE_DEPENDS), and the changed hash regenerates. Nothing here
# applies to TERCO_REGENERATE, which always runs its generators: its
# in-place refresh of src/generated is how the committed set is checked to
# reproduce. And TERCO_FORCE_REGENERATE ignores every stamp for one
# configure.
function(terco_gen what)
  get_property(_mode GLOBAL PROPERTY TERCO_GEN_MODE)
  if(_mode STREQUAL "plan")
    string(JOIN " " _line ${ARGN})
    set_property(GLOBAL APPEND PROPERTY TERCO_GEN_PLAN "${_line}")
  else()
    execute_process(COMMAND ${ARGN} RESULT_VARIABLE _rc OUTPUT_QUIET)
    if(NOT _rc EQUAL 0)
      message(FATAL_ERROR "${what} failed (${_rc})")
    endif()
  endif()
endfunction()

# The stamp's manifest: line 1 the stamp, then one "<sha256> <file>" per
# file in `dest`.
function(terco_gen_write_stamp stamp_file stamp dest)
  file(GLOB _files LIST_DIRECTORIES false RELATIVE "${dest}" "${dest}/*")
  list(SORT _files)
  set(_text "${stamp}\n")
  foreach(_f ${_files})
    file(SHA256 "${dest}/${_f}" _h)
    string(APPEND _text "${_h} ${_f}\n")
  endforeach()
  file(WRITE "${stamp_file}" "${_text}")
endfunction()

function(terco_gen_uptodate out stamp_file stamp dest)
  set(${out} FALSE PARENT_SCOPE)
  if(NOT EXISTS "${stamp_file}")
    return()
  endif()
  file(STRINGS "${stamp_file}" _lines)
  list(POP_FRONT _lines _first)
  if(NOT _first STREQUAL "${stamp}")
    return()
  endif()
  set(_want "")
  foreach(_l ${_lines})
    if(NOT _l MATCHES "^([0-9a-f]+) (.+)$")
      return()
    endif()
    set(_h "${CMAKE_MATCH_1}")
    set(_f "${CMAKE_MATCH_2}")
    if(NOT EXISTS "${dest}/${_f}")
      return()
    endif()
    file(SHA256 "${dest}/${_f}" _have)
    if(NOT _have STREQUAL _h)
      return()
    endif()
    list(APPEND _want "${_f}")
  endforeach()
  file(GLOB _files LIST_DIRECTORIES false RELATIVE "${dest}" "${dest}/*")
  list(SORT _files)
  list(SORT _want)
  if(NOT _files STREQUAL _want)
    return()
  endif()
  set(${out} TRUE PARENT_SCOPE)
endfunction()

# Plan, stamp, then run `gen_macro` (a macro of terco_gen calls) and sync
# `stage` into `dest` -- or skip both. `stampable` OFF always runs (the
# in-place src/generated refresh). The stamp is kept beside the stage, as
# `${stage}.stamp`; `name` is only for the status lines.
function(terco_generate name gen_macro stage dest prune stampable)
  set(_tg_gen_deps ${ARGN})
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_tg_gen_deps})

  set_property(GLOBAL PROPERTY TERCO_GEN_MODE "plan")
  set_property(GLOBAL PROPERTY TERCO_GEN_PLAN "")
  cmake_language(CALL ${gen_macro})
  get_property(_tg_plan GLOBAL PROPERTY TERCO_GEN_PLAN)

  set(_tg_sig "dest ${dest}\nprune ${prune}\npython ${Python3_VERSION}\n")
  foreach(_tg_s ${_tg_gen_deps})
    file(SHA256 "${_tg_s}" _tg_h)
    string(APPEND _tg_sig "script ${_tg_s} ${_tg_h}\n")
  endforeach()
  foreach(_tg_c ${_tg_plan})
    string(APPEND _tg_sig "run ${_tg_c}\n")
  endforeach()
  string(SHA256 _tg_stamp "${_tg_sig}")
  string(SUBSTRING "${_tg_stamp}" 0 12 _tg_short)
  set(_tg_stamp_file "${stage}.stamp")

  if(stampable AND NOT TERCO_FORCE_REGENERATE)
    terco_gen_uptodate(_tg_ok "${_tg_stamp_file}" "${_tg_stamp}" "${dest}")
    if(_tg_ok)
      message(STATUS "terco: ${name} up to date (stamp ${_tg_short})")
      return()
    endif()
  endif()

  message(STATUS "terco: generating ${name} into ${dest} (stamp ${_tg_short})")
  file(REMOVE "${_tg_stamp_file}")
  file(REMOVE_RECURSE "${stage}")
  file(MAKE_DIRECTORY "${stage}")
  set_property(GLOBAL PROPERTY TERCO_GEN_MODE "run")
  cmake_language(CALL ${gen_macro})
  terco_sync_generated("${stage}" "${dest}" ${prune})
  if(stampable)
    terco_gen_write_stamp("${_tg_stamp_file}" "${_tg_stamp}" "${dest}")
  endif()
endfunction()

set(_gen_dir "${CMAKE_CURRENT_SOURCE_DIR}/generated")
set(_regen OFF)
if(NOT TERCO_LMAX STREQUAL "2")
  set(_gen_dir "${CMAKE_BINARY_DIR}/generated")
  set(_regen ON)
elseif(NOT TERCO_CART_ORDER STREQUAL "lexicographic")
  set(_gen_dir "${CMAKE_BINARY_DIR}/generated")
  set(_regen ON)
elseif((NOT "${TERCO_DF_LMAX}" STREQUAL "" AND NOT TERCO_DF_LMAX STREQUAL "2")
       OR NOT TERCO_DF_LMAX_AUX STREQUAL "4")
  # A DF orbital or auxiliary ceiling above the shipped ones, with the
  # four-centre ceiling left at its default: into the build tree, same as
  # the TERCO_LMAX branch above -- never into src/generated, which
  # TERCO_REGENERATE below writes into DELIBERATELY (regenerating the
  # shipped --lmax 2 --lmax-aux 4 set in place to verify it reproduces the
  # committed files). Without this branch, TERCO_DF_LMAX/TERCO_DF_LMAX_AUX
  # alone did nothing (_regen stayed OFF) unless paired with
  # TERCO_REGENERATE, which then overwrote the committed DF kernels with a
  # higher-ceiling set -- the mistake this branch exists to make impossible.
  set(_gen_dir "${CMAKE_BINARY_DIR}/generated")
  set(_regen ON)
elseif(NOT "${TERCO_DF_FORCE_RYS}" STREQUAL "")
  # TESTING ONLY (see the top-level CMakeLists.txt's own comment on this
  # cache var): forcing the Rys body onto classes that would otherwise keep
  # the VRR+HRR still has to go through the generator, so it needs the same
  # build-tree routing as every other non-default DF configuration above --
  # never src/generated, for the same reason.
  set(_gen_dir "${CMAKE_BINARY_DIR}/generated")
  set(_regen ON)
elseif(TERCO_REGENERATE)
  set(_regen ON)
endif()

# `--order` is passed to every generator call below ONLY when it is not the
# default, exactly like `--lmax` is always passed (TERCO_LMAX has no
# "leave it out" state). Suppressing the default keeps a lexicographic
# regeneration's invocation line -- and so its GENERATED header -- identical
# to what every generator emitted before this option existed, which is what
# "regenerates byte for byte" in docs/GENERATORS.md means for the shipped
# set.
set(_order_flag "")
if(NOT TERCO_CART_ORDER STREQUAL "lexicographic")
  set(_order_flag "--order" "${TERCO_CART_ORDER}")
endif()

if(_regen)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  set(_gen_out "${CMAKE_BINARY_DIR}/stage/generated")
  set(_scripts "${CMAKE_SOURCE_DIR}/scripts")

  set(_hgp_lmax ${TERCO_LMAX})
  if(_hgp_lmax GREATER 2)
    set(_hgp_lmax 2)
  endif()
  set(_df_lmax ${TERCO_LMAX})
  if(NOT "${TERCO_DF_LMAX}" STREQUAL "")
    set(_df_lmax ${TERCO_DF_LMAX})
  endif()
  set(_df_lmax_aux ${TERCO_DF_LMAX_AUX})
  set(_df_force_rys_flag "")
  if(NOT "${TERCO_DF_FORCE_RYS}" STREQUAL "")
    set(_df_force_rys_flag "--force-rys" "${TERCO_DF_FORCE_RYS}")
  endif()

  # The generators, in a macro so terco_generate can run it once to plan
  # (record each command line for the stamp) and once to execute.
  macro(_terco_gen_main)
    # Four-centre, one module per (la,lb,lc,ld) class. `--split` matters: at
    # LMAX=2 a single module is a gigabyte of NVVM that does not finish. HGP
    # stops at d whatever TERCO_LMAX is: classes with an f shell are Rys only
    # (see the Rys section below), and HGP's d kernels do not depend on the
    # ceiling -- generated at --lmax 3 they differ from --lmax 2 only in the
    # GENERATED line.
    terco_gen("gen_perclass.py" ${Python3_EXECUTABLE} ${_scripts}/gen_perclass.py
              --lmax ${_hgp_lmax} ${_order_flag} --split ${_gen_out})

    # One-electron, multipoles and density fitting carry angular momentum on two
    # or three centres rather than four, so they grow as L^2 or L^3 and can
    # afford the same ceiling as the four-centre kernels for far less compile.
    #
    # They follow TERCO_LMAX rather than a literal, because a one-electron path
    # that stops below the four-centre one is not a limitation -- it is a wrong
    # Fock matrix. H would be built from whatever one_e_dispatch returns for a
    # key it has no case for, and an SCF converges on that without complaint.
    # DF is the exception that proved too expensive at f: it follows
    # TERCO_DF_LMAX when that is set (see its comment in the top-level
    # CMakeLists.txt), and its entry points stop on a basis above it.
    terco_gen("gen_1e.py" ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py --lmax ${TERCO_LMAX}
              ${_order_flag} -o ${_gen_out}/trc_1e_kernels.F90)
    # One module per class (gen_1e.py's split emitters): as one file the
    # multipoles took 1165 s to compile at LMAX=3 and the gradient kernels
    # over 90 minutes at LMAX=4, single-threaded.
    terco_gen("gen_1e.py --multipole-dir" ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py
              --lmax ${TERCO_LMAX} --multipole-dir ${_gen_out} ${_order_flag})
    # The PCM smeared-charge surface integrals follow TERCO_LMAX for the
    # same reason the 1e kernels do: a continuum that stops below the
    # four-centre ceiling would have no operator for the basis it is given.
    terco_gen("gen_1e.py --pcm" ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py
              --lmax ${TERCO_LMAX} --pcm ${_gen_out}/trc_pcm_kernels.F90 ${_order_flag})
    # The one-electron GRADIENT kernels follow TERCO_LMAX for the same
    # reason trc_1e's own do: a ceiling below the orbital basis is not a
    # limitation, it is a wrong force.
    terco_gen("gen_1e.py --deriv-dir" ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py
              --lmax ${TERCO_LMAX} --deriv-dir ${_gen_out} ${_order_flag})

    # One module per (la,lb|lc)/(la|lc) class, same reason as gen_perclass.py's
    # --split above: at --lmax 3 --lmax-aux 4 (f orbitals, g auxiliary) a single
    # module did not finish compiling on the GPU back end. Classes above
    # gen_df.py's own SHIPPED_L/SHIPPED_LA get a Rys quadrature body instead
    # of the written-out VRR+HRR -- same reason, see that script's module
    # docstring -- so it takes the SAME unroll threshold gen_rys.py's
    # four-centre Rys classes do (TERCO_RYS_UNROLL_MAX), not a new one; the
    # header this embeds (MOD_HEADER) never includes either flag, so this
    # does not disturb the shipped set's byte-for-byte regeneration.
    # TERCO_DF_FORCE_RYS is TESTING ONLY and empty in every normal build (see
    # its own comment, top-level CMakeLists.txt).
    terco_gen("gen_df.py" ${Python3_EXECUTABLE} ${_scripts}/gen_df.py --lmax ${_df_lmax}
              --lmax-aux ${_df_lmax_aux} ${_order_flag} --unroll-max-comp ${TERCO_RYS_UNROLL_MAX}
              ${_df_force_rys_flag} --split ${_gen_out})
    # The DF integrals' nuclear derivatives, contracted (the --df gradient,
    # trc_scf_grad_df.F90): Rys for every class, since a derivative raises an
    # index past the written-out VRR+HRR's compile wall (gen_df_grad.py's
    # docstring). Same ceilings as gen_df.py above; their own unroll
    # threshold, TERCO_DFG_UNROLL_MAX (top-level CMakeLists.txt), with the
    # root loop never unrolled.
    terco_gen("gen_df_grad.py" ${Python3_EXECUTABLE} ${_scripts}/gen_df_grad.py --lmax ${_df_lmax}
              --lmax-aux ${_df_lmax_aux} ${_order_flag} --unroll-max-comp ${TERCO_DFG_UNROLL_MAX}
              --unroll-roots-max-comp 0 --split ${_gen_out})

    # The (l,ic)->(px,py,pz) table shared by trc_hgp.F90's/trc_binkernel.F90's
    # `comp` and trc_batch.F90's `cart_pow`, the flat table
    # trc_collocation.F90's `shell_collocate` reads, and the lexicographic<->
    # this-build permutation trc_cart_perm.F90 wraps for every reader of
    # externally-authored AO data. Needed at every TERCO_LMAX, unlike the
    # unrolled VRR/HRR includes below, because these call sites are not
    # conditional on TRC_UNROLL_VRR/HRR.
    terco_gen("gen_cartcomp.py" ${Python3_EXECUTABLE} ${_scripts}/gen_cartcomp.py ${_order_flag}
              --comp-out ${_gen_out}/trc_cart_comp.inc
              --colloc-out ${_gen_out}/trc_collocation_table.inc
              --perm-l2b-out ${_gen_out}/trc_cart_perm_l2b.inc
              --perm-b2l-out ${_gen_out}/trc_cart_perm_b2l.inc)

    # The Cartesian-to-pure-spherical-harmonic transform blocks, l = 0..4
    # (src/core/trc_sph.F90's only consumer). Fixed at l<=4 regardless of
    # TERCO_LMAX or TRC_MAXL, for the same reason gen_cartcomp.py's CEIL is 4
    # and not TERCO_LMAX: f support (milestone M3) needs this table ready
    # before any four-centre kernel reaches f, and generating a small table
    # four levels higher than today's d ceiling costs nothing.
    terco_gen("gen_sph.py" ${Python3_EXECUTABLE} ${_scripts}/gen_sph.py ${_order_flag}
              -o ${_gen_out}/trc_sph_tables.F90)

    # The unrolled recurrences are textual includes, and only the non-symmetric
    # kernel (one_bin_item_nosym) uses them -- at LMAX<=1, for the reason the top-level CMakeLists gives.
    # They carry a `#if TRC_LMAX != n / #error` guard, so a stale one is a
    # compile error rather than a wrong number. That guard exists because a
    # stale one WAS a wrong number, silently, for a month. gen_hrr.py's output
    # additionally depends on the Cartesian order (its digestion writes AO
    # slots directly), so it also carries a `#if TRC_CART_ORDER != n` guard and
    # takes `--order`; gen_vrr.py's VRR intermediates are addressed only
    # internally and never reach an AO index, so it has no `--order` at all.
    if(TERCO_LMAX LESS_EQUAL 1)
      terco_gen("gen_vrr.py" ${Python3_EXECUTABLE} ${_scripts}/gen_vrr.py --lmax ${TERCO_LMAX}
                --mode include -o ${_gen_out}/trc_vrr_cases.inc)
      terco_gen("gen_hrr.py" ${Python3_EXECUTABLE} ${_scripts}/gen_hrr.py --lmax ${TERCO_LMAX}
                ${_order_flag} -o ${_gen_out}/trc_hrr_cases.inc)
    endif()
  endmacro()

  # src/generated (TERCO_REGENERATE) is refreshed in place and never
  # pruned. TERCO_REGENERATE is never stamped either, in place or (with a
  # non-default ceiling or order) in the build tree: it always runs its
  # generators, see GENERATION STAMPS above.
  set(_prune ON)
  if(_gen_dir STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/generated")
    set(_prune OFF)
  endif()
  set(_stampable ON)
  if(TERCO_REGENERATE)
    set(_stampable OFF)
  endif()
  terco_generate(
    "kernels for LMAX=${TERCO_LMAX}, ${TERCO_CART_ORDER} order" _terco_gen_main
    "${_gen_out}" "${_gen_dir}" ${_prune} ${_stampable}
    ${_scripts}/gen_perclass.py ${_scripts}/gen_vrr.py ${_scripts}/gen_hrr.py
    ${_scripts}/gen_1e.py ${_scripts}/gen_df.py ${_scripts}/gen_df_grad.py ${_scripts}/rys_inline.py
    ${_scripts}/gen_cartcomp.py ${_scripts}/gen_sph.py ${_scripts}/cartorder.py)
endif()

# Always on the include path, not just when regenerating for LMAX<=1: the
# shipped src/generated also carries trc_cart_comp.inc and
# trc_collocation_table.inc now, which trc_hgp.F90, trc_binkernel.F90,
# trc_batch.F90 and trc_collocation.F90 #include unconditionally.
target_include_directories(${main_lib} PRIVATE ${_gen_dir})

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

# --- Rys ----------------------------------------------------------------------
#
# Always built; see the comment above the TERCO_RYS_UNROLL_MAX setting in the
# top-level CMakeLists.txt.
#
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# The root/weight tables: committed, like the Boys table, because they
# need mpmath to regenerate and this project's generators otherwise run
# on a bare Python. In their own subdirectory so the GLOB above -- which
# is not recursive -- never picks them up, ON or OFF.
target_sources(
  ${main_lib}
  PRIVATE eri/trc_rys_select.F90 # TRC_RYS_CLASSES, read once per Fock build
          generated/rys/trc_rys_table1.F90
          generated/rys/trc_rys_table2.F90
          generated/rys/trc_rys_table3.F90
          generated/rys/trc_rys_table4.F90
          generated/rys/trc_rys_table5.F90
          generated/rys/trc_rys_table6.F90
          generated/rys/trc_rys_table7.F90
          generated/rys/trc_rys_table8.F90
          generated/rys/trc_rys_table9.F90)

# Every class with at least one d shell, plus (pp|pp) and (pp|ps) as the
# light end of the HGP/Rys comparison. Seven heavy d classes take their
# chunked Rys kernel by default (gen_rys.py's FOCK_RYS_DEFAULT, measured);
# HGP stays faster everywhere else it has a kernel, including (ss|ss),
# (ps|ss), (ps|ps) and (pp|ss), which is why those have no Rys energy
# kernel at all. The gradient family covers all 81. Regenerated
# at configure time into the build tree, never committed: unlike the root
# tables above, a per-class kernel needs no mpmath, and it is only ever as
# current as the HGP scalar kernel it is spliced from, which this project
# already keeps in-tree-regenerable for exactly this kind of drift.
set(_rys_dir "${CMAKE_BINARY_DIR}/generated_rys")
set(_rys_stage "${CMAKE_BINARY_DIR}/stage/generated_rys")

# The template the Rys generators splice from: a private SKELETON copy of
# the HGP per-class kernels (gen_perclass.py --skeleton: every marker and
# size the splicers read, without HGP's VRR/HRR body they discard),
# written here rather than read out of ${_gen_dir},
# because gen_rys.py and gen_rys_grad.py slice those files on literal
# anchors (`pure subroutine pcsi...`, `!$acc atomic update`) and anything
# that edits src/ in place -- CI's OpenMP port strips `pure` and rewrites
# `!$acc` before configuring -- removes them. Same flags as the ${_gen_dir}
# set, so in an unedited tree the Rys kernels spliced from it are
# byte-identical to ones spliced from the full HGP text (checked for all
# 81 classes up to d, energy and gradient). The skeleton is what makes f
# affordable: the full HGP text at --lmax 3 takes over six minutes and
# 136 MB to write, the skeleton well under a second.
#
# The Rys kernels themselves are deliberately NOT converted: under
# -fopenmp their `do concurrent` is a serial loop and their `!$acc` lines
# comments, and that is correct because the Rys dispatchers are only ever
# called from the sequential class loops in trc_binkernel.F90 and
# trc_grad.F90, never from inside a converted parallel region.
# Rys covers every class up to f. A class with a g shell is the
# McMurchie-Davidson engine's (trc_mdj for J, trc_mdg for K) unless its total
# angular momentum is at most TERCO_RYS_G_LTOT: there Rys has few roots and
# MD's Hermite expansion buys nothing (w20 cc-pVQZ, V100: the g classes up to
# L = 9 ran 2.8-3.8x behind GPU4PySCF's Rys in MD, those from L = 13 ahead).
set(TERCO_RYS_G_LTOT 9 CACHE STRING
    "g classes with la+lb+lc+ld up to this take a Rys exchange kernel (TERCO_LMAX=4)")
set(_rys_lmax ${TERCO_LMAX})
if(_rys_lmax GREATER 3)
  set(_rys_lmax 3)
endif()
set(_all_classes "")
foreach(la RANGE 0 ${TERCO_LMAX})
  foreach(lb RANGE 0 ${TERCO_LMAX})
    foreach(lc RANGE 0 ${TERCO_LMAX})
      foreach(ld RANGE 0 ${TERCO_LMAX})
        math(EXPR _ltot "${la} + ${lb} + ${lc} + ${ld}")
        if(la GREATER _rys_lmax OR lb GREATER _rys_lmax OR lc GREATER _rys_lmax
           OR ld GREATER _rys_lmax)
          # g: only the angular-momentum-canonical orderings the bins form
          # (trc_bins: la >= lb, lc >= ld, bra >= ket), as for trc_mdg below.
          if(_ltot GREATER TERCO_RYS_G_LTOT OR lb GREATER la OR ld GREATER lc
             OR lc GREATER la OR (lc EQUAL la AND ld GREATER lb))
            continue()
          endif()
        endif()
        list(APPEND _all_classes "${la}${lb}${lc}${ld}")
      endforeach()
    endforeach()
  endforeach()
endforeach()
list(JOIN _all_classes "," _all_classes_csv)

set(_hgp_template "${CMAKE_BINARY_DIR}/hgp_template")
file(MAKE_DIRECTORY "${_hgp_template}")

# Every class with a shell of d or above, plus (pp|pp) and (pp|ps). Above
# d these are the only kernels a class has (trc_rys_select.F90 routes
# them unconditionally); up to d, HGP is the default except for
# FOCK_RYS_DEFAULT's classes, and TRC_RYS_CLASSES picks either per class.
set(_rys_classes "")
foreach(_tag ${_all_classes})
  string(REGEX MATCH "[2-9]" _has_d "${_tag}")
  if(_has_d
     OR _tag STREQUAL "1111"
     OR _tag STREQUAL "1110")
    list(APPEND _rys_classes "${_tag}")
  endif()
endforeach()

set(_rys_sources "")
foreach(_tag ${_rys_classes})
  list(APPEND _rys_sources "${_rys_dir}/trc_pc_rys${_tag}.F90")
endforeach()
list(JOIN _rys_classes "," _rys_classes_csv)

list(LENGTH _rys_sources _nrys)
message(STATUS "terco: ${_nrys} Rys class kernels in ${_rys_dir}")
target_sources(${main_lib} PRIVATE ${_rys_sources} ${_rys_dir}/trc_rys_kernels.F90)

# --- Rys GRADIENT kernels: all 81 classes up to (dd|dd), no HGP fallback --
#
# Unlike the energy kernel above (generated only where Rys measures faster
# than HGP), the gradient has no HGP counterpart to fall back to -- an HGP
# gradient is out of scope -- so every class a basis up to LMAX=2 can form
# needs one. Not above d yet: gen_rys_grad.py stops at five roots, so an f
# gradient stops in grad_bins (rysg_has_class) rather than running.
#
set(_grad_lmax ${TERCO_LMAX})
if(_grad_lmax GREATER 2)
  set(_grad_lmax 2)
endif()
set(_rys_grad_classes "")
foreach(la RANGE 0 ${_grad_lmax})
  foreach(lb RANGE 0 ${_grad_lmax})
    foreach(lc RANGE 0 ${_grad_lmax})
      foreach(ld RANGE 0 ${_grad_lmax})
        list(APPEND _rys_grad_classes "${la}${lb}${lc}${ld}")
      endforeach()
    endforeach()
  endforeach()
endforeach()

set(_rys_grad_sources "")
foreach(_tag ${_rys_grad_classes})
  list(APPEND _rys_grad_sources "${_rys_dir}/trc_pc_rysg${_tag}.F90")
endforeach()
list(JOIN _rys_grad_classes "," _rys_grad_classes_csv)

list(LENGTH _rys_grad_sources _nrysg)
message(STATUS "terco: ${_nrysg} Rys gradient class kernels in ${_rys_dir}")
target_sources(${main_lib} PRIVATE ${_rys_grad_sources}
                                   ${_rys_dir}/trc_rys_grad_kernels.F90)

# The Rys generators, planned and run (or skipped) by terco_generate: the
# template first, since every splicer reads it, then one energy and one
# gradient kernel per class and the two dispatchers. The template is not
# part of the stamp's manifest -- it is an intermediate, rewritten
# whenever the Rys kernels are -- and is not needed when they are skipped.
macro(_terco_gen_rys)
  terco_gen("gen_perclass.py (Rys template)" ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/gen_perclass.py --lmax ${TERCO_LMAX} ${_order_flag}
            --classes ${_all_classes_csv} --skeleton --split ${_hgp_template})
  foreach(_tag ${_rys_classes})
    terco_gen("gen_rys.py ${_tag}" ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/gen_rys.py
              ${_tag} ${_order_flag} --src-dir ${_hgp_template}
              --unroll-max-comp ${TERCO_RYS_UNROLL_MAX} -o ${_rys_stage}/trc_pc_rys${_tag}.F90)
  endforeach()
  terco_gen("gen_rys_dispatch.py" ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/gen_rys_dispatch.py --classes ${_rys_classes_csv}
            -o ${_rys_stage}/trc_rys_kernels.F90)
  foreach(_tag ${_rys_grad_classes})
    terco_gen("gen_rys_grad.py ${_tag}" ${Python3_EXECUTABLE}
              ${CMAKE_SOURCE_DIR}/scripts/gen_rys_grad.py ${_tag} ${_order_flag}
              --src-dir ${_hgp_template} --unroll-max-comp ${TERCO_RYS_UNROLL_MAX}
              -o ${_rys_stage}/trc_pc_rysg${_tag}.F90)
  endforeach()
  terco_gen("gen_rys_grad_dispatch.py" ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/gen_rys_grad_dispatch.py
            --classes ${_rys_grad_classes_csv} -o ${_rys_stage}/trc_rys_grad_kernels.F90)
endmacro()
set(_s "${CMAKE_SOURCE_DIR}/scripts")
terco_generate(
  "Rys kernels" _terco_gen_rys "${_rys_stage}" "${_rys_dir}" ON ON
  ${_s}/gen_perclass.py ${_s}/gen_vrr.py ${_s}/gen_hrr.py ${_s}/gen_rys.py
  ${_s}/gen_rys_grad.py ${_s}/rys_inline.py ${_s}/gen_rys_dispatch.py
  ${_s}/gen_rys_grad_dispatch.py ${_s}/cartorder.py)

# --- McMurchie-Davidson J engine ----------------------------------------------
#
# Generated into the build tree at every TERCO_LMAX, like the Rys kernels:
# the ket and bra stages per (la, lb), the Hermite stage per
# (la + lb, lc + ld) with its R batches, and the host dispatcher
# (scripts/gen_md_j.py). Nothing is spliced from another generator, and the
# whole set at LMAX=3 is about 110k lines, a small fraction of the HGP and
# Rys sets, because the engine never forms a Cartesian quartet.
# TERCO_MDJ_BATCHES overrides the R batch count per total L ("7:6,8:9").
set(TERCO_MDJ_BATCHES "" CACHE STRING "R batches per total L for the MD J engine, L:n,...")
set(_mdj_dir "${CMAKE_BINARY_DIR}/generated_mdj")
set(_mdj_stage "${CMAKE_BINARY_DIR}/stage/generated_mdj")
set(_mdj_batch_flag "")
if(NOT "${TERCO_MDJ_BATCHES}" STREQUAL "")
  set(_mdj_batch_flag "--batches" "${TERCO_MDJ_BATCHES}")
endif()
set(TERCO_MDG_CLASSES "2121,1010" CACHE STRING "classes with an MD GEMM Fock prototype (TRC_MDG_CLASSES routes them)")
set(TERCO_MDK_CLASSES "2121,1010" CACHE STRING "classes with an MD Fock prototype kernel (TRC_MDK_CLASSES routes them)")
# Every class with a shell above f goes to the MD GEMM pass, whatever
# TERCO_MDG_CLASSES says: nothing else in terco forms it.
set(_mdg_classes "${TERCO_MDG_CLASSES}")
if(TERCO_LMAX GREATER 3)
  foreach(la RANGE 0 ${TERCO_LMAX})
    foreach(lb RANGE 0 ${la})
      foreach(lc RANGE 0 ${la})
        foreach(ld RANGE 0 ${lc})
          if(lc EQUAL la AND ld GREATER lb)
            continue()
          endif()
          if(la GREATER 3 OR lc GREATER 3)
            string(APPEND _mdg_classes ",${la}${lb}${lc}${ld}")
          endif()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endif()
macro(_terco_gen_mdj)
  terco_gen("gen_md_j.py" ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/gen_md_j.py
            --lmax ${TERCO_LMAX} ${_order_flag} ${_mdj_batch_flag} -o ${_mdj_stage})
  terco_gen("gen_md_k.py" ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/gen_md_k.py
            --lmax ${TERCO_LMAX} --classes ${TERCO_MDK_CLASSES} ${_order_flag} -o ${_mdj_stage})
  terco_gen("gen_md_kgemm.py" ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/gen_md_kgemm.py
            --lmax ${TERCO_LMAX} --classes ${_mdg_classes} ${_order_flag} -o ${_mdj_stage})
endmacro()
terco_generate(
  "MD J kernels" _terco_gen_mdj "${_mdj_stage}" "${_mdj_dir}" ON ON
  ${_s}/gen_md_j.py ${_s}/gen_md_k.py ${_s}/gen_md_kgemm.py ${_s}/md/__init__.py ${_s}/md/hermite.py ${_s}/md/rtree.py ${_s}/cartorder.py)
file(GLOB _mdj_sources "${_mdj_dir}/*.F90")
list(LENGTH _mdj_sources _nmdj)
message(STATUS "terco: ${_nmdj} MD J engine modules in ${_mdj_dir}")
target_sources(${main_lib} PRIVATE ${_mdj_sources})

# TERCO_FORCE_REGENERATE is one-shot: having ignored the stamps for this
# configure, it turns itself off so the next one trusts them again.
if(TERCO_FORCE_REGENERATE)
  set(TERCO_FORCE_REGENERATE OFF CACHE BOOL
      "Regenerate the build-tree kernels once, ignoring the generation stamps" FORCE)
endif()
