# Basis sets shipped with terco

Every `.json` file here is the Basis Set Exchange's own JSON export, unmodified,
from `basis_set_exchange` **0.12**, written by `tools/fetch_basis.py` (which
holds the list) and checked by CI (`tools/fetch_basis.py --check`). File names
are BSE's (`transform_basis_name`): lower case, `*` as `_st_`, `/` as `_sl_`, so
`terco water.xyz -b 6-31+G*` finds `6-31+g_st_.json`.

`cmake --install` copies them to `share/terco/basis_sets`, which `-b NAME`
searches after `$TERCO_BASIS_PATH`. For an fpm build, point
`$TERCO_BASIS_PATH` at this directory.

Which build a family needs: d functions (6-31G*, def2-SVP, cc-pVDZ) run in the
default `TERCO_LMAX=2`; f (def2-TZVP, cc-pVTZ, 6-311++G(3df,3pd)) needs
`TERCO_LMAX=3`; g (def2-QZVP, cc-pVQZ) needs `TERCO_LMAX=4`. A build refuses a
basis above its ceiling and says which to rebuild with. The def2 files carry
the def2 ECPs for the elements that use them.

To add a basis, add its BSE name to `BASES` in `tools/fetch_basis.py` and rerun
it. `test/basis_sets/` is separate: the tests' frozen inputs, some converted
from pyscf with a different contraction layout.

## Citing

Please cite the Basis Set Exchange:

> B. P. Pritchard, D. Altarawy, B. Didier, T. D. Gibson, T. L. Windus,
> "New Basis Set Exchange: An Open, Up-to-Date Resource for the Molecular
> Sciences Community", *J. Chem. Inf. Model.* **59**, 4814-4820 (2019),
> doi:10.1021/acs.jcim.9b00725

and the original paper of each basis set you use: every file lists them under
each element's `references` (or run `bse get-refs <name> txt`).
