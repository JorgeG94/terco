target_sources(
  ${main_lib}
  PRIVATE
          trc_pcm_radii.F90   # cavity radii for the polarisable continuum
          trc_pcm.F90         # the polarisable continuum: cavity, charges, operator
          trc_pcm_grad.F90    # the PCM analytic gradient: grad_qv + grad_nuc + grad_solver
)
