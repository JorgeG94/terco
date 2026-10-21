#!/usr/bin/env python3
"""Generate src/xc/trc_xc_citations.F90: every runnable functional's papers.

Inputs, all read from source rather than retyped:

  - libxc's src/references.c: every reference's text and DOI, keyed by
    libxc's own name for it (`xc_ref_Perdew1996_3865` -> `Perdew1996_3865`).
  - libxc's functional info tables (`xc_func_info_type xc_func_info_<name>`
    in src/*.c): which references each libxc functional cites.
  - ExchCXX's src/libxc.cxx (Kernel:: enum -> libxc id) together with
    tools/gen_xc_registry.py's own enum -> terco-stem map, so each kernel a
    terco functional is built from resolves to the libxc functional it was
    translated from.
  - tools/gen_xc_registry.py's registry itself: the functionals and their
    kernels.

A functional cites, in order: the libxc HYBRID or XC entry of the same name
(B3LYP's Stephens 1994, which none of its B88/LYP/VWN kernels carries), then
each kernel's own references, deduplicated. Only `hyb_*` and `*_xc_*` libxc
entries are matched by name: a plain `x`/`c` name match is a kernel already
reached through the kernel map, or a different functional that happens to
share the name (libxc's MGGA_C_B88 is not terco's b88).

Two small local tables, kept here because nothing upstream states them:
NAME_TO_LIBXC (terco names libxc spells differently: pbe0 is libxc's PBEH)
and LOCAL_REFS (the double hybrids, which libxc does not have; their text
is the registry's own comment for each, with no DOI).

Usage:
    python tools/gen_xc_citations.py --libxc ~/dev/libxc --exchcxx ~/dev/ExchCXX \\
        -o src/xc/trc_xc_citations.F90
"""

from __future__ import annotations

import argparse
import re
import sys
import unicodedata
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import gen_xc_registry  # noqa: E402

LINE_WIDTH = 96
KEY_LEN = 32

# terco names libxc spells differently, for the name match above.
NAME_TO_LIBXC = {"pbe0": "XC_HYB_GGA_XC_PBEH"}

# Not in libxc. Text as the registry's own comments give it; no DOI.
LOCAL_REFS = {
    "Grimme2006_034108": ("S. Grimme, J. Chem. Phys. 124, 034108 (2006)", ""),
    "Karton2008_12868": ("A. Karton, A. Tarnopolsky, J.-F. Lamere, G. C. Schatz, and J. M. L. Martin, "
                         "J. Phys. Chem. A 112, 12868 (2008)", ""),
    "Schwabe2006_4398": ("T. Schwabe and S. Grimme, Phys. Chem. Chem. Phys. 8, 4398 (2006)", ""),
}
LOCAL_FUNCTIONAL_REFS = {
    "b2plyp": ["Grimme2006_034108"],
    "b2gp-plyp": ["Karton2008_12868"],
    "mpw2plyp": ["Schwabe2006_4398"],
}

# Always emitted, for trc_citations' method table: libxc itself, and VV10,
# whose nonlocal kernel terco evaluates on its own grid.
EXTRA_KEYS = ["Lehtola2018_1", "Vydrov2010_244103"]


def parse_references(libxc: Path) -> dict[str, tuple[str, str]]:
    text = (libxc / "src/references.c").read_text()
    pat = re.compile(r'func_reference_type xc_ref_(\w+) = \{\s*"((?:[^"\\]|\\.)*)",\s*"((?:[^"\\]|\\.)*)"')
    return {m.group(1): (clean(m.group(2)), m.group(3)) for m in pat.finditer(text)}


def clean(text: str) -> str:
    # libxc's text ends its author list "A. B. and C. D.,  Journal ..."
    text = text.replace(".,  ", ", ").replace(",  ", ", ")
    text = re.sub(r"\s+", " ", text).strip()
    # Fortran source stays ASCII: en/em dashes become "-", accents drop.
    text = text.replace("\u2013", "-").replace("\u2014", "-")
    text = text.translate(str.maketrans({"\u00f8": "o", "\u00d8": "O", "\u00df": "ss",
                                         "\u00e6": "ae", "\u0142": "l", "\u0111": "d"}))
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    return text


def parse_func_info(libxc: Path) -> dict[str, tuple[str, list[str]]]:
    """libxc id (XC_GGA_X_PBE) -> (info name gga_x_pbe, [reference keys])."""
    pat = re.compile(r'xc_func_info_type xc_func_info_(\w+) = \{\s*(XC_\w+),\s*\w+,\s*"(?:[^"\\]|\\.)*",'
                     r'\s*\w+,\s*\{([^}]*)\}')
    out = {}
    for path in sorted((libxc / "src").glob("*.c")):
        for m in pat.finditer(path.read_text(errors="replace")):
            keys = [r.strip().lstrip("&").removeprefix("xc_ref_") for r in m.group(3).split(",")]
            out[m.group(2)] = (m.group(1), [k for k in keys if k and k != "NULL"])
    return out


def functional_keys(exchcxx: Path, libxc: Path, refs: dict) -> dict[str, list[str]]:
    info = parse_func_info(libxc)
    lx = gen_xc_registry.strip_comments((exchcxx / "src/libxc.cxx").read_text())
    enum_to_xc = dict(re.findall(r"\{\s*Kernel::(\w+)\s*,\s*(XC_\w+)\s*\}", lx))
    stem_to_xc: dict[str, str] = {}
    for enum, stem in gen_xc_registry.parse_enum_to_stem(exchcxx).items():
        if enum in enum_to_xc:
            stem_to_xc.setdefault(stem, enum_to_xc[enum])

    by_name: dict[str, list[str]] = {}
    for xc_id, (name, _) in info.items():
        m = re.match(r"(hyb_)?(?:lda|gga|mgga)_(x|c|xc|k)_(\w+)$", name)
        if m and (m.group(1) or m.group(2) == "xc"):
            by_name.setdefault(m.group(3).replace("_", ""), []).append(xc_id)

    entries, _ = gen_xc_registry.build_registry(exchcxx)
    out = {}
    for e in entries:
        if e.missing:
            continue
        ids = [NAME_TO_LIBXC[e.name]] if e.name in NAME_TO_LIBXC else \
            sorted(by_name.get(e.name.replace("-", "").replace("_", ""), []))
        ids += [stem_to_xc[s] for _, s in e.kernels if s in stem_to_xc]
        keys = list(LOCAL_FUNCTIONAL_REFS.get(e.name, []))
        for xc_id in ids:
            for k in info.get(xc_id, ("", []))[1]:
                if k not in keys and k in refs:
                    keys.append(k)
        if not keys:
            raise SystemExit(f"{e.name}: no reference found")
        out[e.name] = keys
    return out


def fstr(s: str) -> str:
    return '"' + s.replace('"', '""') + '"'


def emit_string(lhs: str, s: str, indent: str) -> list[str]:
    """`lhs = "..."`, split over continuation lines at LINE_WIDTH."""
    first = f"{indent}{lhs} = "
    room = LINE_WIDTH - len(first) - 6
    pieces, cur = [], ""
    for word in s.split(" "):
        add = (" " if cur else "") + word
        if cur and len(cur) + len(add) > room:
            pieces.append(cur + " ")
            cur = word
        else:
            cur += add
    pieces.append(cur)
    lines = [first + fstr(pieces[0]) + (" // &" if len(pieces) > 1 else "")]
    for i, p in enumerate(pieces[1:], 1):
        lines.append(f"{indent}   {fstr(p)}" + (" // &" if i < len(pieces) - 1 else ""))
    return lines


def emit(fkeys: dict[str, list[str]], refs: dict[str, tuple[str, str]], version: str) -> str:
    used = []
    for keys in fkeys.values():
        used += [k for k in keys if k not in used]
    used += [k for k in EXTRA_KEYS if k not in used]
    nmax = max(len(k) for k in fkeys.values())
    L = []
    L += ["!",
          f"! GENERATED by tools/gen_xc_citations.py from libxc {version}'s src/references.c and its",
          "! functionals' info tables, through ExchCXX's kernel map -- do not edit.",
          "!",
          "! Every runnable functional's papers: libxc's own hybrid/XC entry of the same name",
          "! first, then each translated kernel's, deduplicated. The text and DOI are libxc's.",
          "!",
          "module trc_xc_citations",
          "   implicit none",
          "   private",
          "   public :: xc_functional_ref_keys, xc_ref_text",
          f"   integer, parameter, public :: XC_REF_KEY_LEN = {KEY_LEN}",
          f"   integer, parameter, public :: XC_MAX_REFS = {nmax}",
          "   !! The most references any one functional cites.",
          "",
          "contains",
          "",
          "   pure subroutine xc_functional_ref_keys(name, keys, nkeys)",
          "      !! The reference keys a functional cites, by its registry name",
          "      !! (`trc_xc_functional_t%name`); `nkeys` is 0 for an unknown name.",
          "      character(len=*), intent(in) :: name",
          "      character(len=XC_REF_KEY_LEN), intent(out) :: keys(XC_MAX_REFS)",
          "      integer, intent(out) :: nkeys",
          '      keys = ""',
          "      select case (name)"]
    for name, keys in fkeys.items():
        L.append(f"      case ({fstr(name)})")
        L.append(f"         nkeys = {len(keys)}")
        for i, k in enumerate(keys, 1):
            L.append(f"         keys({i}) = {fstr(k)}")
    L += ["      case default",
          "         nkeys = 0",
          "      end select",
          "   end subroutine xc_functional_ref_keys",
          "",
          "   pure subroutine xc_ref_text(key, text, doi, found)",
          "      !! A reference's text and DOI (\"\" when it has none), by key.",
          "      character(len=*), intent(in) :: key",
          "      character(len=:), allocatable, intent(out) :: text, doi",
          "      logical, intent(out) :: found",
          "      found = .true.",
          "      select case (key)"]
    allrefs = dict(refs)
    allrefs.update(LOCAL_REFS)
    for k in used:
        text, doi = allrefs[k]
        L.append(f"      case ({fstr(k)})")
        L += emit_string("text", text, "         ")
        L.append(f"         doi = {fstr(doi)}")
    L += ["      case default",
          "         found = .false.",
          '         text = ""',
          '         doi = ""',
          "      end select",
          "   end subroutine xc_ref_text",
          "",
          "end module trc_xc_citations",
          ""]
    for k in used:
        if len(k) > KEY_LEN:
            raise SystemExit(f"reference key {k} longer than {KEY_LEN}")
    return "\n".join(L)


def libxc_version(libxc: Path) -> str:
    for f in ("VERSION", "configure.ac", "CMakeLists.txt"):
        p = libxc / f
        if p.exists():
            m = re.search(r"(\d+\.\d+\.\d+)", p.read_text())
            if m:
                return m.group(1)
    return "unknown"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--libxc", required=True, type=Path)
    ap.add_argument("--exchcxx", required=True, type=Path)
    ap.add_argument("-o", "--output", required=True, type=Path)
    args = ap.parse_args(argv)
    refs = parse_references(args.libxc)
    fkeys = functional_keys(args.exchcxx, args.libxc, refs)
    args.output.write_text(emit(fkeys, refs, libxc_version(args.libxc)))
    nrefs = len({k for v in fkeys.values() for k in v})
    print(f"{args.output}: {len(fkeys)} functionals, {nrefs} references")
    return 0


if __name__ == "__main__":
    sys.exit(main())
