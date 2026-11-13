target_sources(
  ${main_lib}
  PRIVATE
          bit_repro.f90           # bit-reproducible exp/log/pow/erfc/cbrt, vendored
          bit_repro_helpers.f90   # its submodule
)
