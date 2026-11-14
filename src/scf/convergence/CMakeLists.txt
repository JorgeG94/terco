target_sources(
  ${main_lib}
  PRIVATE
          trc_scf_diis.F90       # CDIIS and ADIIS, the level shift and the damped step
          trc_scf_gdm.F90        # geometric direct minimisation: L-BFGS history, preconditioner, step control
          trc_scf_gdm_step.F90   # GDM as a step mode of the SCF loop: accept-on-descent orbital rotations
          trc_scf_soscf.F90      # second-order SCF (--soscf): trust-region Newton on the orbital rotations, RHF
)
