target_sources(
  ${main_lib}
  PRIVATE
          trc_guess_json.F90  # a SAD guess from mqc's per-element JSON blocks
          trc_sad.F90  # the SAD guess built here, one atomic SCF per element
          trc_sap.F90  # SAP: free-atom potentials, superposed on a quadrature
          trc_minao.F90  # the minao guess: a minimal-ANO density projected into the target basis
)
