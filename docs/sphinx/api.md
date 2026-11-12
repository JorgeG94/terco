# API reference

Generated from the public API sources themselves
(`src/api/trc_api.F90`, `src/api/trc_capi.F90`) by
[sphinx-fortran-domain](https://github.com/fortran-lang/sphinx-fortran-domain)'s
FORD-backed lexer, so this page and the source cannot drift apart silently.
`src/generated/` (the per-class integral kernels) is never parsed here --
see `AGENTS.md`: it is generator output, not something anyone reads as API.

## `trc_api` -- the Fortran API

`trc_api` is a thin facade: it `use`s and re-exports the public names of
three lower-layer modules (`trc_basis_types`, `trc_device_binding`,
`trc_integrals`) so that a Fortran caller needs only `use trc_api`. Because
every type and procedure it exposes is *defined* in one of those other
modules, not in `trc_api.F90` itself, the autodoc view below is sparse: it
shows the module's own header comment and nothing else, since the domain
only has FORD-level doc comments to show for entities declared in the file
it parses. The module's own header comment (`src/api/trc_api.F90`, lines
1-45) explains the facade and its design in prose; the public names it
re-exports are listed there too.

```{f:module} trc_api
```

## `trc_capi` -- the C ABI (Fortran declaration)

`trc_capi.F90` is where the `bind(c)` entry points actually live, each one
documented with a FORD `!!` block placed, per terco's convention, as the
first lines of the procedure body (so it documents the procedure, not its
first argument). This is the module sphinx-fortran-domain's FORD lexer
exercises most: hundreds of `bind(c)` functions, most with a doc comment.

```{f:module} trc_capi
```

The RI-CCSD(T) entry points (`trc_set_rccsd`, `trc_set_rccsd_convergence`,
`trc_set_rccsd_algorithm`, `trc_run_rccsd`, `trc_rccsd_energy`,
`trc_rccsd_dimensions`; ABI 1.28) are among them; {doc}`library` has an example
and the list of what each call returns.
