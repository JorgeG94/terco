target_sources(
  ${main_lib}
  PRIVATE
          trc_sym_group.F90    # the eight Abelian point groups: operations, characters, irrep names
          trc_sym_tops.F90     # pyscf's axis search for symmetric/spherical tops, linear molecules, atoms
          trc_sym_geom.F90     # point-group detection, orientation, symmetrization of a geometry
          trc_sym_ao.F90       # AO signed permutations and the SO basis
          trc_sym_pairs.F90    # shell-pair orbits under the group: images, representatives, weights, orbit-invariant Q and dsh
          trc_sym_scf.F90      # SO-blocked SCF step: per-irrep diagonalisation, labels, symmetrisation
          trc_sym_frame.F90    # exact AO rotation between the standard frame and the user's input frame
          trc_sym_run.F90     # what a run asks for: off/auto/default/group, subgroup rules, per-irrep occupations
)
