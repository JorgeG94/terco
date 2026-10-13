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
          trc_grad.F90       # two-electron gradient launch plan, Rys only
          trc_api.F90        # basis, pair list, 1e, multipoles, DF
          trc_decontract.F90 # splitting a general contraction into primitives
          trc_fock.F90       # the four-centre Fock builds
          trc_linalg.F90     # dense linear algebra on device-resident matrices
          trc_scf.F90        # the SCF: HF and Kohn-Sham, R and U, resident
          trc_rimp2.F90      # RI-MP2 on the RHF reference, resident
          trc_basis_json.F90 # a basis from a MolSSI BSE JSON file
          trc_guess_json.F90 # a SAD guess from mqc's per-element JSON blocks
          trc_sad.F90 # the SAD guess built here, one atomic SCF per element
          trc_sap.F90 # SAP: free-atom potentials, superposed on a quadrature
          trc_capi.F90       # the bind(c) skin
          trc_cart_perm.F90  # AO-index permutation, lexicographic <-> this build's order
          trc_ao_norm.F90    # AO-index normalisation, this build's <-> GAMESS's per-component convention
          trc_sph.F90        # the spherical (pure) AO space: block-sparse Cartesian <-> spherical transforms
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

set(_gen_dir "${CMAKE_CURRENT_SOURCE_DIR}/generated")
set(_regen OFF)
if(NOT TERCO_LMAX STREQUAL "2")
  set(_gen_dir "${CMAKE_BINARY_DIR}/generated")
  set(_regen ON)
elseif(NOT TERCO_CART_ORDER STREQUAL "lexicographic")
  set(_gen_dir "${CMAKE_BINARY_DIR}/generated")
  set(_regen ON)
elseif(NOT "${TERCO_DF_LMAX}" STREQUAL "" AND NOT TERCO_DF_LMAX STREQUAL "2")
  # A DF ceiling above the shipped one, with the four-centre ceiling left at
  # its default: into the build tree, same as the TERCO_LMAX branch above --
  # never into src/generated, which TERCO_REGENERATE below writes into
  # DELIBERATELY (regenerating the shipped --lmax 2 --lmax-aux 3 set in
  # place to verify it reproduces the committed files). Without this branch,
  # TERCO_DF_LMAX alone did nothing (_regen stayed OFF) unless paired with
  # TERCO_REGENERATE, which then overwrote the committed DF kernels with a
  # higher-ceiling set -- the mistake this branch exists to make impossible.
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
  file(REMOVE_RECURSE "${_gen_out}")
  file(MAKE_DIRECTORY "${_gen_out}")
  set(_scripts "${CMAKE_SOURCE_DIR}/scripts")

  message(STATUS "terco: generating kernels for LMAX=${TERCO_LMAX} "
                 "CART_ORDER=${TERCO_CART_ORDER} in ${_gen_dir}")

  # Four-centre, one module per (la,lb,lc,ld) class. `--split` matters: at
  # LMAX=2 a single module is a gigabyte of NVVM that does not finish. HGP
  # stops at d whatever TERCO_LMAX is: classes with an f shell are Rys only
  # (see the Rys section below), and HGP's d kernels do not depend on the
  # ceiling -- generated at --lmax 3 they differ from --lmax 2 only in the
  # GENERATED line.
  set(_hgp_lmax ${TERCO_LMAX})
  if(_hgp_lmax GREATER 2)
    set(_hgp_lmax 2)
  endif()
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_perclass.py --lmax
            ${_hgp_lmax} ${_order_flag} --split ${_gen_out}
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
  # DF is the exception that proved too expensive at f: it follows
  # TERCO_DF_LMAX when that is set (see its comment in the top-level
  # CMakeLists.txt), and its entry points stop on a basis above it.
  execute_process(COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py --lmax ${TERCO_LMAX}
                          ${_order_flag} -o ${_gen_out}/trc_1e_kernels.F90
                  RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_1e.py failed (${_rc})")
  endif()
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py --lmax ${TERCO_LMAX} --multipole
            ${_gen_out}/trc_mult_kernels.F90 ${_order_flag} RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_1e.py --multipole failed (${_rc})")
  endif()
  set(_df_lmax ${TERCO_LMAX})
  if(NOT "${TERCO_DF_LMAX}" STREQUAL "")
    set(_df_lmax ${TERCO_DF_LMAX})
  endif()
  # One module per (la,lb|lc)/(la|lc) class, same reason as gen_perclass.py's
  # --split above: at --lmax 3 --lmax-aux 4 (f orbitals, g auxiliary) a single
  # module did not finish compiling on the GPU back end.
  execute_process(COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_df.py --lmax ${_df_lmax}
                          ${_order_flag} --split ${_gen_out}
                  RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_df.py failed (${_rc})")
  endif()

  # The (l,ic)->(px,py,pz) table shared by trc_hgp.F90's/trc_binkernel.F90's
  # `comp` and trc_batch.F90's `cart_pow`, the flat table
  # trc_collocation.F90's `shell_collocate` reads, and the lexicographic<->
  # this-build permutation trc_cart_perm.F90 wraps for every reader of
  # externally-authored AO data. Needed at every TERCO_LMAX, unlike the
  # unrolled VRR/HRR includes below, because these call sites are not
  # conditional on TRC_UNROLL_VRR/HRR.
  execute_process(
    COMMAND
      ${Python3_EXECUTABLE} ${_scripts}/gen_cartcomp.py ${_order_flag}
      --comp-out ${_gen_out}/trc_cart_comp.inc --colloc-out
      ${_gen_out}/trc_collocation_table.inc --perm-l2b-out
      ${_gen_out}/trc_cart_perm_l2b.inc --perm-b2l-out
      ${_gen_out}/trc_cart_perm_b2l.inc
    RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_cartcomp.py failed (${_rc})")
  endif()

  # The Cartesian-to-pure-spherical-harmonic transform blocks, l = 0..4
  # (src/trc_sph.F90's only consumer). Fixed at l<=4 regardless of
  # TERCO_LMAX or TRC_MAXL, for the same reason gen_cartcomp.py's CEIL is 4
  # and not TERCO_LMAX: f support (milestone M3) needs this table ready
  # before any four-centre kernel reaches f, and generating a small table
  # four levels higher than today's d ceiling costs nothing.
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_sph.py ${_order_flag}
            -o ${_gen_out}/trc_sph_tables.F90
    RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_sph.py failed (${_rc})")
  endif()

  # The unrolled recurrences are textual includes, and only the shared kernel
  # uses them -- at LMAX<=1, for the reason the top-level CMakeLists gives.
  # They carry a `#if TRC_LMAX != n / #error` guard, so a stale one is a
  # compile error rather than a wrong number. That guard exists because a
  # stale one WAS a wrong number, silently, for a month. gen_hrr.py's output
  # additionally depends on the Cartesian order (its digestion writes AO
  # slots directly), so it also carries a `#if TRC_CART_ORDER != n` guard and
  # takes `--order`; gen_vrr.py's VRR intermediates are addressed only
  # internally and never reach an AO index, so it has no `--order` at all.
  if(TERCO_LMAX LESS_EQUAL 1)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_vrr.py --lmax ${TERCO_LMAX}
              --mode include -o ${_gen_out}/trc_vrr_cases.inc OUTPUT_QUIET)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} ${_scripts}/gen_hrr.py --lmax ${TERCO_LMAX}
              ${_order_flag} -o ${_gen_out}/trc_hrr_cases.inc OUTPUT_QUIET)
  endif()

  set(_prune ON)
  if(_gen_dir STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/generated")
    set(_prune OFF)
  endif()
  terco_sync_generated("${_gen_out}" "${_gen_dir}" ${_prune})
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

# --- Rys (on by default) -----------------------------------------------------
#
# See the TERCO_RYS option's own comment (top-level CMakeLists.txt). OFF
# adds nothing: no root tables, no per-class kernel, no dispatcher, and
# nothing above this line is affected either way.
#
if(TERCO_RYS)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)

  # The root/weight tables: committed, like the Boys table, because they
  # need mpmath to regenerate and this project's generators otherwise run
  # on a bare Python. In their own subdirectory so the GLOB above -- which
  # is not recursive -- never picks them up, ON or OFF.
  target_sources(
    ${main_lib}
    PRIVATE trc_rys_select.F90 # TRC_RYS_CLASSES, read once per Fock build
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
  # light end of the HGP/Rys comparison: measured on a 450-atom water cluster
  # in 6-31G*, Rys only beats HGP at (dd|dd) (1.11x) and is 1.1-7.4x slower
  # everywhere below, (pp|pp) the slowest, which is why HGP stays the energy
  # default for every class. The gradient family covers all 81. Regenerated
  # at configure time into the build tree, never committed: unlike the root
  # tables above, a per-class kernel needs no mpmath, and it is only ever as
  # current as the HGP scalar kernel it is spliced from, which this project
  # already keeps in-tree-regenerable for exactly this kind of drift.
  set(_rys_dir "${CMAKE_BINARY_DIR}/generated_rys")
  set(_rys_stage "${CMAKE_BINARY_DIR}/stage/generated_rys")
  file(REMOVE_RECURSE "${_rys_stage}")
  file(MAKE_DIRECTORY "${_rys_stage}")

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
  set(_all_classes "")
  foreach(la RANGE 0 ${TERCO_LMAX})
    foreach(lb RANGE 0 ${TERCO_LMAX})
      foreach(lc RANGE 0 ${TERCO_LMAX})
        foreach(ld RANGE 0 ${TERCO_LMAX})
          list(APPEND _all_classes "${la}${lb}${lc}${ld}")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
  list(JOIN _all_classes "," _all_classes_csv)

  set(_hgp_template "${CMAKE_BINARY_DIR}/hgp_template")
  file(MAKE_DIRECTORY "${_hgp_template}")
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/gen_perclass.py
            --lmax ${TERCO_LMAX} ${_order_flag} --classes ${_all_classes_csv}
            --skeleton --split ${_hgp_template}
    RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_perclass.py (Rys template) failed (${_rc})")
  endif()

  # Every class with a shell of d or above, plus (pp|pp) and (pp|ps). Above
  # d these are the only kernels a class has (trc_rys_select.F90 routes
  # them unconditionally); up to d, HGP stays the default and these exist
  # for TRC_RYS_CLASSES.
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
    set(_rys_out "${_rys_stage}/trc_pc_rys${_tag}.F90")
    execute_process(
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/gen_rys.py
              ${_tag} ${_order_flag} --src-dir ${_hgp_template}
              --unroll-max-comp ${TERCO_RYS_UNROLL_MAX} -o ${_rys_out}
      RESULT_VARIABLE _rc OUTPUT_QUIET)
    if(NOT _rc EQUAL 0)
      message(FATAL_ERROR "gen_rys.py ${_tag} failed (${_rc})")
    endif()
    list(APPEND _rys_sources "${_rys_dir}/trc_pc_rys${_tag}.F90")
  endforeach()

  list(JOIN _rys_classes "," _rys_classes_csv)
  set(_rys_disp "${_rys_stage}/trc_rys_kernels.F90")
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/gen_rys_dispatch.py
            --classes ${_rys_classes_csv} -o ${_rys_disp}
    RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_rys_dispatch.py failed (${_rc})")
  endif()

  list(LENGTH _rys_sources _nrys)
  message(STATUS "terco: ${_nrys} Rys class kernels generated into ${_rys_dir}")
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
    set(_rys_grad_out "${_rys_stage}/trc_pc_rysg${_tag}.F90")
    execute_process(
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/gen_rys_grad.py
              ${_tag} ${_order_flag} --src-dir ${_hgp_template}
              --unroll-max-comp ${TERCO_RYS_UNROLL_MAX} -o ${_rys_grad_out}
      RESULT_VARIABLE _rc OUTPUT_QUIET)
    if(NOT _rc EQUAL 0)
      message(FATAL_ERROR "gen_rys_grad.py ${_tag} failed (${_rc})")
    endif()
    list(APPEND _rys_grad_sources "${_rys_dir}/trc_pc_rysg${_tag}.F90")
  endforeach()

  list(JOIN _rys_grad_classes "," _rys_grad_classes_csv)
  set(_rys_grad_disp "${_rys_stage}/trc_rys_grad_kernels.F90")
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/gen_rys_grad_dispatch.py
            --classes ${_rys_grad_classes_csv} -o ${_rys_grad_disp}
    RESULT_VARIABLE _rc OUTPUT_QUIET)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "gen_rys_grad_dispatch.py failed (${_rc})")
  endif()

  list(LENGTH _rys_grad_sources _nrysg)
  message(STATUS "terco: ${_nrysg} Rys gradient class kernels generated into ${_rys_dir}")
  target_sources(${main_lib} PRIVATE ${_rys_grad_sources}
                                     ${_rys_dir}/trc_rys_grad_kernels.F90)

  terco_sync_generated("${_rys_stage}" "${_rys_dir}" ON)
endif()
