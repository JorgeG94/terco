"""Unit tests for the do concurrent -> OpenMP converter.

Run with: pytest tools/test_dc_to_omp.py
"""
import importlib.util
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
_SPEC = importlib.util.spec_from_file_location(
    "dc_to_omp", pathlib.Path(__file__).parent / "dc_to_omp.py")
dc_to_omp = importlib.util.module_from_spec(_SPEC)
# dc_to_omp defines a @dataclass, which looks its module up in sys.modules.
sys.modules["dc_to_omp"] = dc_to_omp
_SPEC.loader.exec_module(dc_to_omp)


def test_masked_header_becomes_guarded_body():
    # A skipped loop is not harmless: amdflang runs an unconverted do
    # concurrent on the host, over host copies of device-resident arrays
    # (four loops in trc_scf_df.F90 did). The mask guards the innermost body
    # so the loops stay perfectly nested for collapse.
    src = [
        "      do concurrent(mu=1:n, nu=1:n, nu <= mu) local(t)\n",
        "         a(mu, nu) = 1\n",
        "      end do\n",
    ]
    out, n = dc_to_omp.transform_lines(src, "gpu")
    text = "".join(out)
    assert n == 1
    assert "do concurrent" not in text
    assert "!$omp target teams distribute parallel do collapse(2) private(t)" in text
    lines = [l.strip() for l in out]
    i = lines.index("do nu = 1, n")
    assert lines[i + 1] == "if (nu <= mu) then"
    assert lines[i + 3] == "end if"


def test_unmasked_header_unchanged_shape():
    src = ["      do concurrent(i=1:n)\n", "         a(i) = 1\n", "      end do\n"]
    out, n = dc_to_omp.transform_lines(src, "gpu")
    assert n == 1
    assert "if (" not in "".join(out)
