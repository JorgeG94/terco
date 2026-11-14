target_sources(
  ${main_lib}
  PRIVATE
          trc_api.F90         # the Fortran API: re-exports the three above
          trc_capi.F90        # the bind(c) skin
)
