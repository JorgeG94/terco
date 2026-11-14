target_sources(
  ${main_lib}
  PRIVATE
          trc_rimp2_driver.F90       # RI-MP2 on the RHF reference, resident
          trc_rimp2_relax.F90        # the RI-MP2 gradient's fitted side: amplitudes, Lagrangian, 3c/2c terms
          trc_rimp2_grad.F90                 # the analytic RI-MP2 gradient: Z-vector, relaxed densities, assembly
)
