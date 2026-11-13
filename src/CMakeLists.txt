#
# The library: a hand-written core, and kernels that are generated.
#
# Roughly nine tenths of the compiled lines here come out of `scripts/`. That
# is the design and not an accident -- see docs/GENERATORS.md -- so this file
# is mostly about deciding WHICH generated set to use.
#

add_subdirectory(core)

# Source registrations live beside their Fortran sources.
add_subdirectory(eri)
add_subdirectory(api)
add_subdirectory(basis)
add_subdirectory(solvation)
add_subdirectory(scf)
add_subdirectory(post_hf)
add_subdirectory(symmetry)
add_subdirectory(grid)
add_subdirectory(xc)
include("${PROJECT_SOURCE_DIR}/cmake/terco_generation.cmake")

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
  set(_scripts "${PROJECT_SOURCE_DIR}/scripts")

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
    # The multipoles stop at d whatever TERCO_LMAX is: only the response
    # code reads them (dipole perturbations), the f and g classes cost
    # 20 and 41 minutes of nvfortran, and trc_multipoles stops on a basis
    # above MULT_LMAX rather than returning a wrong matrix.
    set(_mult_lmax ${TERCO_LMAX})
    if(_mult_lmax GREATER 2)
      set(_mult_lmax 2)
    endif()
    terco_gen("gen_1e.py --multipole-dir" ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py
              --lmax ${_mult_lmax} --multipole-dir ${_gen_out} ${_order_flag})
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
    # The PCM surface-integral GRADIENT kernels, same ceiling as --pcm
    # above for the same reason: a continuum force that stops below the
    # orbital basis is a wrong force. One module PER CLASS, same reason
    # --deriv-dir above is split rather than single-module: the monolithic
    # dispatcher sized every thread's registers for the heaviest compiled
    # class regardless of which one it ran (ncu, see the PR that split
    # this; measured 1.08 s/call against the no-derivative PCM Fock term's
    # 0.08 s on the same integral family, local-memory spills).
    terco_gen("gen_1e.py --pcm-grad-dir" ${Python3_EXECUTABLE} ${_scripts}/gen_1e.py
              --lmax ${TERCO_LMAX} --pcm-grad-dir ${_gen_out} ${_order_flag})

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
    # `comp` and trc_batch.F90's `batch_cart_pow`, the flat table
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

    # The ECP angular path's own Cartesian-to-pure-spherical-harmonic
    # coefficients, l = 0..8 (src/eri/trc_ecp_angular.F90's only consumer).
    # Fixed, like gen_sph.py's table above, regardless of TERCO_LMAX: an
    # ECP channel's own angular momentum (up to ECP_LMAX = 5,
    # trc_ecp_quadrature.F90) can push the total past TERCO_LMAX's own
    # ceiling even at TERCO_LMAX=0.
    terco_gen("gen_ecp_tables.py" ${Python3_EXECUTABLE} ${_scripts}/gen_ecp_tables.py
              -o ${_gen_out}/trc_ecp_tables.F90)

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
    ${_scripts}/gen_cartcomp.py ${_scripts}/gen_sph.py ${_scripts}/gen_ecp_tables.py ${_scripts}/cartorder.py)
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
          generated/schwarz/trc_schwarz_rys.F90 # Rys diagonal Schwarz kernel (scripts/gen_schwarz.py); no LMAX dependence, own subdir like the root tables below so the GLOB never picks it up twice
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
#
# Only in the PAIR-CANONICAL orderings la >= lb and lc >= ld. Every bin holds
# its pairs with the higher-l shell first (`lorder` in trc_bins.F90) and the
# plans copy a bin's (bin_la, bin_lb) straight into (sLA, sLB) and (sLC, sLD)
# (bins_plan_emit), so no launch site can present lb > la or ld > lc and a
# kernel for one of those orderings can never run. That is 38 of 67 classes at
# TERCO_LMAX=2 and 149 of 242 at TERCO_LMAX=3, each worth up to 60 kB of
# per-thread stack and 255 registers of compile time.
#
# NOT the bra >= ket half of the gradient filter below, which would be wrong
# here: the skeleton (symmetry) Fock plan deliberately launches a bin pair with
# its bra and ket SWAPPED when that side has fewer representatives
# (bins_skel_cut_first, up to SKEL_SWAP_LMAX = 3 precisely because "the kernels
# of every ordering exist up to f"), which is a class with (lc, ld) > (la, lb).
# Removing those would silently send such a launch to HGP up to d and
# `error stop` above it ("a class above d has no Rys kernel in this build",
# trc_rys_select.F90).
#
# The gradient filter below keeps that third condition, and what makes it
# safe is that EVERY gradient plan is emitted ANGULAR-MOMENTUM CANONICALLY:
# bins_plan_emit's `canon` reorders the two sides of each segment so that
# (sLA, sLB) >= (sLC, sLD), whichever side the walk cut (trc_bins.F90).
# grad_bins asks for it from both plans it builds, the full one
# (bins_plan_build's `canonical`) and the skeleton one; the energy plans keep
# the walk order, which the bra<ket energy classes above serve.
#
# Ascending bin index alone is not enough. plan%build iterates b%live in
# ascending index, and bin_key packs (la, lb) as its most significant fields,
# so ascending index USUALLY means (la, lb) >= (lc, ld). But split_large_bins
# and arrange_orbits append a split bin's extra chunks PAST every original
# index while keeping their parent's (la, lb), so a chunk of a low-l bin can
# be the bra against a high-l ket -- the same mechanism that makes the 42
# bra<ket energy classes above reachable -- and grad_classes would `error
# stop` on a class it has no kernel for. At TERCO_LMAX=2 all 81 gradient
# classes are built, so nothing can fire; above d it takes an f shell and a
# bin of more than 32768 pairs, which test/check_sym_grad.F90's capped case
# builds at TERCO_LMAX=3.
# The 30 gradient classes that filter drops (NOT ~42: that is the energy
# side's bra<ket count above, over a different class set) are the
# f-containing ones with (lc, ld) > (la, lb):
#   0030 0031 0032 0033 1030 1031 1032 1033 1130 1131 1132 1133
#   2030 2031 2032 2033 2130 2131 2132 2133 2230 2231 2232 2233
#   3031 3032 3033 3132 3133 3233
set(_rys_classes "")
foreach(_tag ${_all_classes})
  string(SUBSTRING "${_tag}" 0 1 _tag_la)
  string(SUBSTRING "${_tag}" 1 1 _tag_lb)
  string(SUBSTRING "${_tag}" 2 1 _tag_lc)
  string(SUBSTRING "${_tag}" 3 1 _tag_ld)
  if(_tag_lb GREATER _tag_la OR _tag_ld GREATER _tag_lc)
    continue()
  endif()
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

# Parity skipping: the z-masked variants of the
# Rys-routed classes that get one (gen_perclass.py's PARITY_RYS_TOP, whose
# comment has the coverage measurement that chose them; the HGP ones are in
# ${_gen_dir}), where this build has them. Keep this list and PARITY_RYS_TOP
# in step: a tag missing here only means the launch falls back to the unmasked
# kernel, which is the right answer and the slower one.
set(_rys_masked_classes "")
foreach(_tag 2111 2120 2121 2210 2211 2221 2222)
  if(_tag IN_LIST _rys_classes)
    list(APPEND _rys_masked_classes "${_tag}")
    list(APPEND _rys_sources "${_rys_dir}/trc_pc_rys${_tag}z.F90")
  endif()
endforeach()
list(JOIN _rys_masked_classes "," _rys_masked_classes_csv)

list(LENGTH _rys_sources _nrys)
message(STATUS "terco: ${_nrys} Rys class kernels in ${_rys_dir}")
target_sources(${main_lib} PRIVATE ${_rys_sources} ${_rys_dir}/trc_rys_kernels.F90)

# --- Rys GRADIENT kernels: every class the pair bins launch -------------
#
# Unlike the energy kernel above (generated only where Rys measures faster
# than HGP), the gradient has no HGP counterpart to fall back to -- an HGP
# gradient is out of scope -- so every class a basis can form needs one.
# Up to d that is all 81 tags, as before. Above d only the
# ANGULAR-MOMENTUM CANONICAL ones (la >= lb, lc >= ld, (la, lb) >= (lc,
# ld); trc_bins.F90) are built: grad_bins never launches the others, and
# at LMAX=3 they would be ~175 dead kernels of compile time. Not above f:
# no g gradient has been validated, so a TERCO_LMAX=4 build stops at f and
# trc_capi refuses a g basis up front (RYSG_LMAX).
#
set(_grad_lmax ${TERCO_LMAX})
if(_grad_lmax GREATER 3)
  set(_grad_lmax 3)
endif()
set(_rys_grad_classes "")
foreach(la RANGE 0 ${_grad_lmax})
  foreach(lb RANGE 0 ${_grad_lmax})
    foreach(lc RANGE 0 ${_grad_lmax})
      foreach(ld RANGE 0 ${_grad_lmax})
        set(_keep TRUE)
        if(la GREATER 2 OR lb GREATER 2 OR lc GREATER 2 OR ld GREATER 2)
          math(EXPR _ab "${la}*16 + ${lb}")
          math(EXPR _cd "${lc}*16 + ${ld}")
          if(lb GREATER la OR ld GREATER lc OR _cd GREATER _ab)
            set(_keep FALSE)
          endif()
        endif()
        if(_keep)
          list(APPEND _rys_grad_classes "${la}${lb}${lc}${ld}")
        endif()
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
            ${PROJECT_SOURCE_DIR}/scripts/gen_perclass.py --lmax ${TERCO_LMAX} ${_order_flag}
            --classes ${_all_classes_csv} --skeleton --split ${_hgp_template})
  foreach(_tag ${_rys_classes})
    terco_gen("gen_rys.py ${_tag}" ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/gen_rys.py
              ${_tag} ${_order_flag} --src-dir ${_hgp_template}
              --unroll-max-comp ${TERCO_RYS_UNROLL_MAX} -o ${_rys_stage}/trc_pc_rys${_tag}.F90)
  endforeach()
  foreach(_tag ${_rys_masked_classes})
    terco_gen("gen_rys.py ${_tag} --dead z" ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/gen_rys.py
              ${_tag} ${_order_flag} --src-dir ${_hgp_template} --dead z
              --unroll-max-comp ${TERCO_RYS_UNROLL_MAX} -o ${_rys_stage}/trc_pc_rys${_tag}z.F90)
  endforeach()
  terco_gen("gen_rys_dispatch.py" ${Python3_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/scripts/gen_rys_dispatch.py --classes ${_rys_classes_csv}
            --masked-classes "${_rys_masked_classes_csv}"
            -o ${_rys_stage}/trc_rys_kernels.F90)
  foreach(_tag ${_rys_grad_classes})
    terco_gen("gen_rys_grad.py ${_tag}" ${Python3_EXECUTABLE}
              ${PROJECT_SOURCE_DIR}/scripts/gen_rys_grad.py ${_tag} ${_order_flag}
              --src-dir ${_hgp_template} --unroll-max-comp ${TERCO_RYS_UNROLL_MAX}
              -o ${_rys_stage}/trc_pc_rysg${_tag}.F90)
  endforeach()
  terco_gen("gen_rys_grad_dispatch.py" ${Python3_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/scripts/gen_rys_grad_dispatch.py
            --classes ${_rys_grad_classes_csv} -o ${_rys_stage}/trc_rys_grad_kernels.F90)
endmacro()
set(_s "${PROJECT_SOURCE_DIR}/scripts")
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
# The P <-> Q symmetric Hermite kernels (see src/eri/trc_mdj.F90) are the
# default for every class (99 is above any lp + lq); -1 builds the gather
# kernels alone. TERCO_MDJ_SYM_LMAX >= 0 uses them for every class with lp + lq
# up to it;
# TERCO_MDJ_GATHER_ALL also keeps every gather kernel, so TRC_MDJ_SYM_LMAX can
# race the two forms at run time (that build compiles both sets).
set(TERCO_MDJ_SYM_LMAX "99" CACHE STRING "MD J: largest lp+lq on the symmetric Hermite kernels (99: all, -1: none)")
option(TERCO_MDJ_GATHER_ALL "MD J: keep every gather Hermite kernel beside the symmetric ones" OFF)
# Always passed: the generator's own default is 99, so leaving the flag out
# for -1 built the symmetric kernels anyway.
list(APPEND _mdj_batch_flag "--sym-lmax" "${TERCO_MDJ_SYM_LMAX}")
if(TERCO_MDJ_GATHER_ALL)
  list(APPEND _mdj_batch_flag "--gather-all")
endif()
# MD exchange (trc_mdg, scripts/gen_md_k.py) takes the g classes Rys does not:
# total L above TERCO_RYS_G_LTOT (V100, w20 cc-pVQZ: MD won 19 of the 25 g
# classes at L 10-12 against the best Rys chunk shape, and every class above).
# Only the angular-momentum-canonical orderings the bins form are generated.
# Routing d and f classes here was measured and dropped: it won only four
# (ff|xx) classes on one system and none on a contracted one.
set(_mdg_classes "")
if(TERCO_LMAX GREATER 3)
  foreach(la RANGE 0 ${TERCO_LMAX})
    foreach(lb RANGE 0 ${la})
      foreach(lc RANGE 0 ${la})
        foreach(ld RANGE 0 ${lc})
          if(lc EQUAL la AND ld GREATER lb)
            continue()
          endif()
          math(EXPR _ltot "${la} + ${lb} + ${lc} + ${ld}")
          if((la GREATER 3 OR lc GREATER 3) AND _ltot GREATER TERCO_RYS_G_LTOT)
            list(APPEND _mdg_classes "${la}${lb}${lc}${ld}")
          endif()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endif()
list(JOIN _mdg_classes "," _mdg_classes)
macro(_terco_gen_mdj)
  terco_gen("gen_md_j.py" ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/gen_md_j.py
            --lmax ${TERCO_LMAX} ${_order_flag} ${_mdj_batch_flag} -o ${_mdj_stage})
  terco_gen("gen_md_k.py" ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/gen_md_k.py
            --lmax ${TERCO_LMAX} --classes=${_mdg_classes} ${_order_flag} -o ${_mdj_stage})
endmacro()
terco_generate(
  "MD J kernels" _terco_gen_mdj "${_mdj_stage}" "${_mdj_dir}" ON ON
  ${_s}/gen_md_j.py ${_s}/gen_md_k.py ${_s}/md/__init__.py ${_s}/md/hermite.py ${_s}/md/rtree.py ${_s}/cartorder.py)
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
