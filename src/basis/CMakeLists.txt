target_sources(
  ${main_lib}
  PRIVATE
          trc_decontract.F90  # splitting a general contraction into primitives
          trc_ecp_types.F90  # trc_ecp_t: the ECP data model, flat rows, separate from the orbital shells
          trc_ecp_json.F90  # an effective core potential from a MolSSI BSE JSON file
          trc_basis_json.F90  # a basis from a MolSSI BSE JSON file
)
