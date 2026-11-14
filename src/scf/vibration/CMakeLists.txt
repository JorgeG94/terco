target_sources(
  ${main_lib}
  PRIVATE
          trc_vib.F90                 # harmonic vibrational analysis of a Cartesian Hessian, host only
          trc_thermo.F90              # ideal-gas RRHO thermochemistry from trc_vib's frequencies
)
