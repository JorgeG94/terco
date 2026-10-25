#!/usr/bin/env python3
"""Tests for refhash.py: which references an edit to generate.py restales.

Standard library only, like check_staleness.py; registered as the ctest
terco/check_refhash. Each test edits a copy of the real generate.py text,
re-parses it, and compares every CASES entry's code digest against the
unedited one.

    python3 test/pyscf_refs/test_refhash.py
"""
from __future__ import annotations

import ast
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refhash  # noqa: E402

SOURCE = (HERE / "generate.py").read_text()
BASE = refhash.Generator(SOURCE)
BASE_DIGESTS = {c: BASE.code_digest(c) for c in BASE.cases}


def restaled(source: str) -> set[str]:
    gen = refhash.Generator(source)
    assert set(gen.cases) == set(BASE.cases)
    return {c for c in gen.cases if gen.code_digest(c) != BASE_DIGESTS[c]}


def edit_once(source: str, old: str, new: str) -> str:
    assert source.count(old) == 1, f"expected exactly one {old!r} in generate.py"
    return source.replace(old, new)


def body_line(source: str, func: str) -> str:
    """The first line of `func`'s body after its docstring, as written."""
    for node in ast.parse(source).body:
        if isinstance(node, ast.FunctionDef) and node.name == func:
            stmt = node.body[1] if ast.get_docstring(node) else node.body[0]
            return source.splitlines(keepends=True)[stmt.lineno - 1]
    raise AssertionError(f"no def {func} in generate.py")


class RestaleTests(unittest.TestCase):

    def test_unedited_restales_nothing(self):
        self.assertEqual(restaled(SOURCE), set())

    def test_generator_body_restales_only_its_case(self):
        line = body_line(SOURCE, "gen_scf_df")
        src = edit_once(SOURCE, line, line + "    _refhash_probe = 1\n")
        self.assertEqual(restaled(src), {"scf_df"})

    def test_shared_helper_restales_exactly_its_users(self):
        # converge_tight is called by gen_rks, gen_scf_mgga, gen_scf_robust, gen_full_grad,
        # gen_rimp2_uhf, gen_dh, gen_sph_corr, gen_vib_thermo (so all three
        # vib_thermo_* cases, which share that one generator), gen_opt_ts
        # (opt_ts_hcn), and through _gen_rimp2_one / _gen_scf_df_one /
        # _gen_cphf_one / _gen_cpks_one by gen_rimp2 / gen_scf_df / gen_cphf /
        # gen_cphf_uhf / gen_cpks. Update this set when a new generator
        # starts calling it.
        src = edit_once(SOURCE, "conv_tol_grad=1.0e-8, max_cycle=200,",
                        "conv_tol_grad=1.0e-8, max_cycle=201,")
        self.assertEqual(restaled(src),
                         {"rks", "scf_mgga", "scf_robust", "rimp2", "scf_df", "full_grad",
                          "rimp2_uhf", "dh", "sph_corr", "cphf", "cphf_uhf", "cpks",
                          "vib_thermo_water", "vib_thermo_hcn", "vib_thermo_nh3_planar",
                          "opt_ts_hcn"})

    def test_case_arguments_restale_only_that_case(self):
        src = edit_once(SOURCE, '"basis_sets/cc-pvdz.json", 1e-9),',
                        '"basis_sets/cc-pvdz.json", 1e-8),')
        self.assertEqual(restaled(src), {"rys_fock_water_ccpvdz"})

    def test_comments_docstrings_and_layout_restale_nothing(self):
        src = SOURCE
        # A comment on its own line and at the end of a code line.
        src = edit_once(src, "def converge_tight(", "# a new comment\ndef converge_tight(")
        line = body_line(src, "gen_scf_df")
        src = edit_once(src, line, line.rstrip("\n") + "  # trailing comment\n")
        # A docstring rewritten, and blank lines added between definitions.
        src = edit_once(src, '    """Refuse an mf whose own `converged` flag is unset.',
                        '    """Refuse an mf whose `converged` flag is unset (reworded).')
        src = edit_once(src, "\n\ndef gen_rks():", "\n\n\n\ndef gen_rks():")
        self.assertEqual(restaled(src), set())

    def test_header_names_the_generator_function(self):
        self.assertEqual(BASE.function("scf_df"), "gen_scf_df")
        self.assertEqual(BASE.function("rys_fock_water_ccpvdz"), "gen_rys_fock")

    def test_imports_hash_the_same_on_every_python(self):
        # 3.15 gave Import/ImportFrom an `is_lazy` field (PEP 810); a plain
        # import must serialise as it did before, or every generator with a
        # function-local import restales when CI's Python moves.
        tree = ast.parse("import numpy\nfrom pyscf import gto")
        self.assertEqual(refhash._canon(tree.body[0]),
                         "Import(names=[alias(name='numpy')])")
        self.assertEqual(refhash._canon(tree.body[1]),
                         "ImportFrom(module='pyscf',names=[alias(name='gto')],level=0)")


if __name__ == "__main__":
    unittest.main(verbosity=2)
