# Examples

Each directory is one run: a TOML input (the format is
[docs/INPUT.md](../docs/INPUT.md)) and the geometries it names. Paths in a
file are relative to the file, so a run can start from anywhere:

```sh
terco examples/water/water.toml
```

Basis names (`def2-svp`, `sto-3g`, ...) are looked up as on the command
line: in `$TERCO_BASIS_PATH`, then in the basis sets an install puts in
`share/terco/basis_sets`. From a build tree, point at the repository's own:

```sh
export TERCO_BASIS_PATH=$PWD/test/basis_sets
```

| Example | What it shows | Needs |
|---|---|---|
| [`water/`](water/water.toml) | the minimal input: RHF/def2-SVP single point | |
| [`formamide_water/`](formamide_water/formamide.toml) | B3LYP, density fitting, C-PCM water, analytic gradient, Hirshfeld charges | |
| [`water_opt_freq/`](water_opt_freq/water_opt_freq.toml) | geometry optimisation, then frequencies and thermochemistry at the minimum | DL-FIND build |
| [`hcn_neb/`](hcn_neb/hcn_neb.toml) | HCN to HNC climbing-image NEB with a seeded initial path | DL-FIND build |

"DL-FIND build" means CMake with `-DTERCO_ENABLE_DLFIND=ON` (or the fpm opt-in
in the top-level README); without it `--opt` is refused with a message.

Every run writes a JSON summary whose `"input"` object lists every setting
it used, defaults included, in the same sections as the file.
