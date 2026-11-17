target_sources(
  ${main_lib}
  PRIVATE
          trc_cphf.F90                # coupled-perturbed HF: the orbital-response solve and the static polarizability
          trc_cpks.F90                # coupled-perturbed KS: trc_cphf plus the XC kernel response, RKS LDA/GGA
          trc_tddft.F90               # TDA/CIS and full-RPA singlet excitation energies, block Davidson over fock_nosym_many + fxc_apply_many
)
