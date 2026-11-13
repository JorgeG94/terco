target_sources(
  ${main_lib}
  PRIVATE
          trc_scf_grad_df.F90        # the density-fitted (--df) two-electron gradient: fitted J and K
          trc_scf_grad.F90           # the full analytic gradient: 1e + eri%gradient + nuclear repulsion
)
