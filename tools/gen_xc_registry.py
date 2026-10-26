#!/usr/bin/env python3
"""Generate src/xc/trc_xc_functional.F90 and docs/functionals.md.

Inputs, all read from source rather than retyped:

  - tools/xc_kernel_discovery.py: which ExchCXX kernel headers
    tools/exchcxx_to_fortran.py can translate (and why it refuses the
    rest), and which translated kernels disagree with libxc at check_xc's
    tolerances despite translating cleanly.
  - ExchCXX's own functional table, ~/dev/ExchCXX/src/xc_functional.cxx:
    every named functional's kernel mix and HybCoeffs {alpha, beta, omega}.
  - ExchCXX's own kernel enum -> libxc-kernel-and-builtin-class mapping,
    ~/dev/ExchCXX/src/builtin_interface.cxx and the kernel headers, to turn
    a `Kernel::FOO` reference in the functional table into the terco stem
    (`foo`) its translated module uses.

Plus three small local tables, kept here rather than derived because
nothing upstream states them: LOCAL (terco-only entries -- svwn_rpa,
b3lyp/b3lyp5's libxc-not-ExchCXX split, wb97x-v's VV10 coefficients --
carried over from the hand-written registry this replaces, unchanged),
OVERRIDES (one entry: ExchCXX's own revB3LYP composition names
Kernel::VWN5, which its own builtin_interface.cxx says is the VWN_RPA
formula, but matches libxc's actual revB3LYP/revB3LYP5 only through the
plain VWN formula -- checked numerically against pyscf; see the commit
that added this generator), and LOCAL_KERNELS (hand-written kernel
modules ExchCXX does not define at all, e.g. wb97m_v_xc -- a LOCAL entry
may list one exactly like a translated stem; see the module's own header
for why it exists and where its coefficients come from).

HybCoeffs -> terco's (exx, exx_lr, omega): ExchCXX's alpha is the fraction
of LONG-range (erf(omega r)/r) exact exchange, beta the amount by which the
SHORT-range fraction differs from it, so alpha+beta is the FULL-range exact
exchange fraction (what libxc calls hybrid_coeff) and -beta is the extra
long-range fraction terco calls exx_lr. omega=0 hybrids (beta=0) reduce to
exx=alpha, exx_lr=0 exactly as the hand-written registry had them.

A functional any of whose kernels failed to translate, disagrees with
libxc, or needs a dispersion correction (B97-D and the wB97X-D/-D3 pair:
ExchCXX's kernel is the semilocal part only, same as wb97x_v_xc excludes
VV10's nonlocal term, but the empirical correction itself is not something
terco computes) is still registered, with `missing` set to why.

Usage:
    python tools/gen_xc_registry.py --exchcxx ~/dev/ExchCXX \
        -o src/xc/trc_xc_functional.F90 --docs docs/functionals.md
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import xc_kernel_discovery  # noqa: E402

INDENT = "   "
LINE_WIDTH = 96
USE_WIDTH = 112  # a `use` line that fits stays one line

# --------------------------------------------------------------------------
# terco-only entries, unchanged from the hand-written registry. `kernels`
# is [(coef, stem), ...]; `exx`/`exx_lr`/`omega`/`nlc_b`/`nlc_c` default to
# 0. `comment` is emitted verbatim above the case, FORD style.
# --------------------------------------------------------------------------


@dataclass
class Entry:
    name: str
    aliases: tuple[str, ...]
    family: str  # "lda" | "gga" | "mgga"
    kernels: list[tuple[float, str]]
    exx: float = 0.0
    exx_lr: float = 0.0
    omega: float = 0.0
    nlc_b: float = 0.0
    nlc_c: float = 0.0
    pt2_cos: float = 0.0
    pt2_css: float = 0.0
    comment: list[str] = field(default_factory=list)
    missing: str | None = None
    source: str = "local"  # local | exchcxx
    after: str | None = None  # emit right after this (ExchCXX) entry, not in LOCAL's place


LOCAL: list[Entry] = [
    Entry("lda_x", ("slater",), "lda", [(1.0, "slater_exchange")]),
    Entry("svwn", ("lda", "svwn5"), "lda", [(1.0, "slater_exchange"), (1.0, "vwn")]),
    Entry("svwn_rpa", (), "lda", [(1.0, "slater_exchange"), (1.0, "vwn_rpa")]),
    Entry("pbe", (), "gga", [(1.0, "pbe_x"), (1.0, "pbe_c")]),
    Entry("blyp", (), "gga", [(1.0, "b88"), (1.0, "lyp")]),
    Entry("pbe0", (), "gga", [(0.75, "pbe_x"), (1.0, "pbe_c")], exx=0.25),
    Entry(
        "b3lyp", (), "gga",
        [(0.08, "slater_exchange"), (0.72, "b88"), (0.19, "vwn_rpa"), (0.81, "lyp")],
        exx=0.20,
        comment=["libxc HYB_GGA_XC_B3LYP: VWN_RPA in the correlation"],
    ),
    Entry(
        "b3lyp5", (), "gga",
        [(0.08, "slater_exchange"), (0.72, "b88"), (0.19, "vwn"), (0.81, "lyp")],
        exx=0.20,
        comment=["libxc HYB_GGA_XC_B3LYP5 and ExchCXX's B3LYP: VWN5"],
    ),
    Entry(
        "bhandhlyp", ("bhhlyp",), "gga", [(0.5, "b88"), (1.0, "lyp")], exx=0.5,
        comment=["libxc HYB_GGA_XC_BHANDHLYP (Gaussian's BHandHLYP): 0.5 B88 + LYP,",
                 "half exact exchange. Not in ExchCXX's functional table."],
    ),
    # Double hybrids: a_x exact exchange plus (1 - a_x) of a GGA exchange
    # kernel, (1 - c) of a GGA correlation kernel plus c of RI-MP2
    # correlation on the converged KS orbitals (c_os = c_ss = c for all
    # three -- no spin-component scaling, unlike a DSD-type double
    # hybrid). Not in libxc (metalquicha's mqc_xc_spec.f90 carries the
    # same citations and the same a_x/c split; terco's own PT2 part is
    # driven by trc_capi.F90's trc_run_rimp2/trc_energy, never libxc).
    Entry(
        "b2plyp", (), "gga", [(0.47, "b88"), (0.73, "lyp")], exx=0.53,
        pt2_cos=0.27, pt2_css=0.27,
        comment=["Grimme, J. Chem. Phys. 124, 034108 (2006): a_x = 0.53, c = 0.27.",
                 "Not in libxc (no double hybrid is); the hybrid part is 0.53 HF +",
                 "0.47 B88 exchange, 0.73 LYP correlation, and c is the RI-MP2",
                 "fraction trc_run_rimp2/trc_energy add on top (pt2_cos = pt2_css,",
                 "no spin-component scaling)."],
    ),
    Entry(
        "b2gp-plyp", ("b2gpplyp",), "gga", [(0.35, "b88"), (0.64, "lyp")], exx=0.65,
        pt2_cos=0.36, pt2_css=0.36,
        comment=["Karton, Tarnopolsky, Lamere, Schatz and Martin, J. Phys. Chem. A",
                 "112, 12868 (2008): a_x = 0.65, c = 0.36."],
    ),
    Entry(
        "mpw2plyp", ("mpw2-plyp",), "gga", [(0.45, "mpw91_x"), (0.75, "lyp")], exx=0.55,
        pt2_cos=0.25, pt2_css=0.25,
        comment=["Schwabe and Grimme, Phys. Chem. Chem. Phys. 8, 4398 (2006):",
                 "a_x = 0.55, c = 0.25; the GGA exchange kernel is mPW91's, not",
                 "B88's -- the other difference from b2plyp/b2gp-plyp above."],
        after="mpw1k",
    ),
    Entry(
        "hflyp", ("hf-lyp", "hf_lyp"), "gga", [(1.0, "lyp")], exx=1.0,
        comment=["Full exact exchange plus LYP correlation; pyscf's 'HF,LYP'",
                 "(libxc GGA_C_LYP with hybrid_coeff 1). Not in ExchCXX's table."],
    ),
    Entry(
        "wb97x-v", ("wb97xv", "wb97x_v"), "gga", [(1.0, "wb97x_v_xc")],
        exx=0.167, exx_lr=0.833, omega=0.3, nlc_b=6.0, nlc_c=0.01,
        comment=[
            "libxc HYB_GGA_XC_WB97X_V. rsh_coeff (omega, alpha, beta) =",
            "(0.3, 1.0, -0.833): alpha is the LONG-RANGE exact-exchange",
            "fraction, hyb = alpha + beta the FULL-range one that multiplies",
            "SR_HFX = K_full - K_erf. pyscf/libxc build",
            "  vk = hyb*K_full + (alpha - hyb)*K_erf(omega)",
            "which is exx*K_full + exx_lr*K_erf with exx = hyb, exx_lr = alpha",
            "- hyb. VV10 (nlc_b, nlc_c) from libxc's nlc_coeff, ((6.0, 0.01), 1).",
        ],
    ),
    Entry(
        "wb97x", (), "gga", [(1.0, "wb97x_xc")], exx=0.157706, exx_lr=0.842294, omega=0.3,
        comment=[
            "libxc HYB_GGA_XC_WB97X. rsh_coeff (omega, alpha, beta) =",
            "(0.3, 1.0, -0.842294); hyb = alpha + beta = 0.157706, exx_lr =",
            "alpha - hyb = 0.842294. No VV10 (nlc_coeff returns empty). The",
            "kernel's own translated `beta` parameter is -0.842294 exactly,",
            "the same direct cross-check wb97x-v's gave.",
        ],
    ),
    Entry(
        "wb97", (), "gga", [(1.0, "wb97_xc")], exx=0.0, exx_lr=1.0, omega=0.4,
        comment=[
            "libxc HYB_GGA_XC_WB97. rsh_coeff (omega, alpha, beta) =",
            "(0.4, 1.0, -1.0); hyb = alpha + beta = 0, exx_lr = alpha - hyb =",
            "1 -- a pure range-separated functional, no global HF mixing. No",
            "VV10. The kernel's own translated `beta` is -1.0 exactly.",
        ],
    ),
    Entry(
        "wb97m-v", ("wb97mv", "wb97m_v"), "mgga", [(1.0, "wb97m_v_xc")],
        exx=0.15, exx_lr=0.85, omega=0.3, nlc_b=6.0, nlc_c=0.01,
        comment=[
            "libxc HYB_MGGA_XC_WB97M_V (id 531). rsh_coeff (omega, alpha,",
            "beta) = (0.3, 1.0, -0.85): hyb = alpha + beta = 0.15 = exx,",
            "exx_lr = alpha - hyb = 0.85, same reading as wb97x-v above.",
            "VV10 (nlc_b, nlc_c) = (6.0, 0.01), wb97x-v's own values, left",
            "untouched by this functional's fit (Mardirossian and",
            "Head-Gordon, J. Chem. Phys. 144, 214110 (2016), Section 5).",
            "The kernel itself is HAND-WRITTEN (trc_xc_wb97m_v_xc.F90's own",
            "header): ExchCXX has no wB97M-V, so there is nothing to",
            "translate here, unlike every other entry in this table.",
        ],
    ),
]

# Kernels nothing upstream defines: hand-written, not translated (see each
# module's own header for why and its source). `available_kernels` treats
# these exactly like an ExchCXX-discovered stem -- same family cache, same
# kernel-id numbering, same `use trc_xc_<stem>` line -- so a LOCAL Entry
# above can list one in its `kernels` just like a translated stem.
LOCAL_KERNELS: list[tuple[str, str]] = [
    ("wb97m_v_xc", "mgga"),
]

# ExchCXX names its own revB3LYP kernel Kernel::VWN5 -- which its own
# src/libxc.cxx maps to XC_LDA_C_VWN_RPA -- but the coefficients (0.13,
# 0.67, 0.16, 0.84) are libxc's revB3LYP5's, and libxc's revB3LYP5
# evaluated against pyscf matches LDA_C_VWN (terco's plain `vwn`), not
# LDA_C_VWN_RPA, to machine precision; the VWN_RPA substitution is off by
# ~1e-2 relative at the same points. ExchCXX's X3LYP, which also uses
# Kernel::VWN5, DOES match LDA_C_VWN_RPA exactly, so this is a property of
# revB3LYP specifically, not a wrong enum->kernel mapping in general.
OVERRIDES: dict[tuple[str, str], str] = {("revB3LYP", "VWN5"): "vwn"}

# Needs an empirical dispersion correction (Grimme D2/D3) terco does not
# have: ExchCXX's kernel is the bare semilocal XC term these functionals
# were parametrised assuming a dispersion correction sits alongside, the
# same way wb97x_v_xc excludes VV10's nonlocal piece.
DISPERSION_ONLY = {
    "B97D": "needs an empirical dispersion correction (Grimme D2), not in terco",
    "B97D3ZERO": "needs an empirical dispersion correction (Grimme D3(0)), not in terco",
    "wB97XD": "needs an empirical dispersion correction (Grimme D2), not in terco",
    "wB97XD3": "needs an empirical dispersion correction (Grimme D3), not in terco",
}


# --------------------------------------------------------------------------
# Parsing ExchCXX's sources.
# --------------------------------------------------------------------------


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def parse_enum_to_stem(exchcxx: Path) -> dict[str, str]:
    """ExchCXX's Kernel::FOO -> terco's stem (the hpp file's own name)."""
    bi = strip_comments((exchcxx / "src/builtin_interface.cxx").read_text())
    pat = re.compile(r"Kernel::(\w+)\s*\)\s*\n?\s*return std::make_unique<Builtin(\w+)>")
    enum_to_builtin = dict(pat.findall(bi))

    kdir = exchcxx / xc_kernel_discovery.x2f.KERNEL_DIR
    builtin_to_stem: dict[str, str] = {}
    for path in kdir.glob("*.hpp"):
        text = path.read_text()
        for m in re.finditer(r"struct Builtin(\w+)\s*:", text):
            builtin_to_stem[m.group(1)] = path.stem

    return {enum: builtin_to_stem[b] for enum, b in enum_to_builtin.items() if b in builtin_to_stem}


@dataclass
class ExchCXXFunctional:
    name: str
    kernels: list[tuple[float, str]]  # (coef, Kernel:: enum name)
    hyb: tuple[float, float, float] | None  # (alpha, beta, omega)


def parse_functional_table(exchcxx: Path) -> list[ExchCXXFunctional]:
    text = strip_comments((exchcxx / "src/xc_functional.cxx").read_text())
    start = text.index("functional_factory")
    body = text[start:]
    pat = re.compile(r"(?:else\s+)?if\s*\(\s*func\s*==\s*Functional::(\w+)\s*\)\s*\{?")
    matches = list(pat.finditer(body))
    out = []
    for i, m in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(body)
        block = body[m.end():end]
        kerns = re.findall(
            r"\{\s*([\d.eE+-]+)\s*,\s*XCKernel\(\s*backend\s*,\s*Kernel::(\w+)\s*,\s*polar\s*\)\s*\}", block)
        hyb_m = re.search(r"hyb_coefs\s*=\s*\{([^}]+)\}", block)
        hyb = tuple(float(x) for x in hyb_m.group(1).split(",")) if hyb_m else None
        out.append(ExchCXXFunctional(m.group(1), [(float(c), k) for c, k in kerns], hyb))
    return out


# --------------------------------------------------------------------------
# Building the registry.
# --------------------------------------------------------------------------


def build_registry(exchcxx: Path) -> tuple[list[Entry], list[str]]:
    """Return (entries, notes-for-docs)."""
    translated, refused = xc_kernel_discovery.discover(exchcxx)
    translated_stems = {s for s, _ in translated}
    disagrees = xc_kernel_discovery.DISAGREES_WITH_LIBXC
    enum_to_stem = parse_enum_to_stem(exchcxx)
    functionals = parse_functional_table(exchcxx)

    entries = list(LOCAL)
    local_names = {n.lower() for e in LOCAL for n in (e.name, *e.aliases)}
    notes = []

    for fn in functionals:
        lname = fn.name.lower()
        if lname in local_names:
            notes.append(f"- `{fn.name}`: name already used by terco's own entry above "
                         f"(content may differ -- libxc, not ExchCXX, is terco's reference).")
            continue

        if fn.name in DISPERSION_ONLY:
            # Still resolve the kernel so the semilocal part is on record,
            # but a dispersion correction is the one and only reason given.
            entries.append(Entry(lname, (), "gga", [], missing=DISPERSION_ONLY[fn.name], source="exchcxx"))
            continue

        kernel_stems = []
        problems: dict[str, str] = {}  # stem_or_enum -> reason
        for coef, enum in fn.kernels:
            stem = OVERRIDES.get((fn.name, enum)) or enum_to_stem.get(enum)
            if stem is None:
                problems[enum] = "no built-in kernel implementation"
            elif stem in refused:
                problems[enum] = refused[stem]
            elif stem in disagrees:
                problems[enum] = disagrees[stem]
            elif stem not in translated_stems:
                problems[enum] = "not translated"
            else:
                kernel_stems.append((coef, stem))

        family = "gga"
        if kernel_stems:
            fams = {f for s, f in translated if s in {st for _, st in kernel_stems}}
            if fams == {"lda"}:
                family = "lda"
            elif "mgga" in fams:
                # Any meta-GGA kernel makes the whole functional one (PW86B95
                # is a GGA exchange under BC95's meta-GGA correlation): the
                # integrator then evaluates tau and calls xc_eval_point_mgga,
                # which dispatches the LDA/GGA kernels too.
                family = "mgga"

        exx = exx_lr = omega = 0.0
        if fn.hyb:
            alpha, beta, omega = fn.hyb
            exx = alpha + beta
            exx_lr = -beta

        if problems:
            # Group enums by identical reason so "needs R2SCANL_C, R2SCANL_X:
            # needs the density Laplacian, ..." prints once.
            by_reason: dict[str, list[str]] = {}
            for enum, reason in problems.items():
                by_reason.setdefault(reason, []).append(enum)
            missing = "; ".join(f"needs {', '.join(sorted(enums))}: {reason}"
                                 for reason, enums in by_reason.items())
            entries.append(Entry(lname, (), family, kernel_stems, exx, exx_lr, omega,
                                  missing=missing, source="exchcxx"))
        else:
            entries.append(Entry(lname, (), family, kernel_stems, exx, exx_lr, omega, source="exchcxx"))

    for e in [e for e in entries if e.after]:
        entries.remove(e)
        entries.insert(next(i for i, x in enumerate(entries) if x.name == e.after) + 1, e)
    return entries, notes


# --------------------------------------------------------------------------
# Emitting Fortran.
# --------------------------------------------------------------------------


def wrap(line: str, indent: str) -> list[str]:
    if len(line) <= LINE_WIDTH:
        return [line]
    out, cur, cont = [], line, indent + INDENT
    while len(cur) > LINE_WIDTH:
        k = cur.rfind(",", 0, LINE_WIDTH - 2)
        if k <= len(indent):
            k = cur.rfind(" ", 0, LINE_WIDTH - 2)
        if k <= len(indent):
            break
        out.append(cur[:k + 1].rstrip() + " &")
        cur = cont + cur[k + 1:].lstrip()
    out.append(cur)
    return out


def fortran_real(x: float) -> str:
    if x == int(x) and abs(x) < 1e15:
        return f"{x:.1f}_dp"
    return f"{x!r}_dp"


def emit_fortran(entries: list[Entry], available_kernels: list[str]) -> str:
    kid_name = {stem: f"XC_K_{stem.upper()}" for stem in available_kernels}
    semilocal = [st for st in available_kernels if _family_of(st) != "mgga"]
    # What xc_eval_point_mgga dispatches: every meta-GGA kernel, plus the
    # LDA/GGA kernels some available meta-GGA composition puts beside them.
    mgga_used = {st for e in entries if e.family == "mgga" and e.missing is None for _, st in e.kernels}
    mgga_dispatch = [st for st in semilocal if st in mgga_used] + \
        [st for st in available_kernels if _family_of(st) == "mgga"]
    lines = [
        "!",
        "! A density functional, as the integrator sees it: a short list of",
        "! (coefficient, kernel) pairs, a family, and a fraction of exact exchange.",
        "!",
        "! GENERATED by tools/gen_xc_registry.py from ExchCXX's functional table",
        "! (~/dev/ExchCXX/src/xc_functional.cxx) plus the generator's own small",
        "! local table of terco-only entries -- do not edit.",
        "!",
        "! The kernels are the translated ExchCXX routines in the sibling modules.",
        "! This module knows their ids and how to call each one at a point; the",
        "! compositions below are libxc's -- libxc is the reference every energy is",
        "! checked against, so its definition of B3LYP is the one that counts, VWN",
        "! variant included. ExchCXX's own B3LYP uses VWN5 and is `b3lyp5` here.",
        "!",
        "! A functional whose kernels did not all translate, disagree with libxc,",
        "! or need a capability terco does not have (an empirical dispersion",
        "! correction, say) is still registered here: `missing` names exactly what",
        "! is missing, and xc_functional_by_name still returns it with `ok = .true.`",
        "! -- the caller (trc_capi's trc_set_method) is the one that refuses it, by",
        "! checking `missing`, rather than reporting an unrecognised name. See",
        "! docs/functionals.md for the full list.",
        "!",
        "! `xc_eval_point` is the only thing that runs inside a kernel. It takes the",
        "! composition as plain arrays rather than the type, because a derived type",
        "! with allocatable components is exactly what a `do concurrent` body should",
        "! not touch.",
        "!",
        "! A meta-GGA has its own pair, `xc_eval_point_mgga` and its polarised twin:",
        "! tau in, vtau out, dispatching every kernel a meta-GGA composition uses",
        "! (its meta-GGA kernels and any LDA/GGA ones beside them). The LDA/GGA",
        "! evaluators dispatch only LDA and GGA kernels, so the meta-GGA family",
        "! leaves what they compile to unchanged. Kernel ids are LDA/GGA first and",
        "! meta-GGA after, for the same reason.",
        "!",
        "! A double hybrid (b2plyp, b2gp-plyp, mpw2plyp) is still just a (coef,",
        "! kernel) composition plus `exx` for its SCF -- `pt2_cos`/`pt2_css` are",
        "! the RI-MP2 fraction on top, which this module only STORES; running it",
        "! (always RI, on the converged KS orbitals) is trc_capi.F90's",
        "! trc_run_rimp2/trc_energy, never anything in here.",
        "!",
        "module trc_xc_functional",
        f"{INDENT}use pic_types, only: dp",
    ]
    for stem in available_kernels:
        if _family_of(stem) == "mgga":
            # No second derivatives translated for a meta-GGA.
            lines.append(f"{INDENT}use trc_xc_{stem}, only: {stem}_exc_vxc_unpolar, {stem}_exc_vxc_polar")
            continue
        one = (f"{INDENT}use trc_xc_{stem}, only: {stem}_exc_vxc_unpolar, {stem}_exc_vxc_polar, "
               f"{stem}_vxc_fxc_unpolar")
        if len(one) <= USE_WIDTH:
            lines.append(one)
        else:
            lines.append(f"{INDENT}use trc_xc_{stem}, only: {stem}_exc_vxc_unpolar, {stem}_exc_vxc_polar, &")
            lines.append(f"{INDENT}{INDENT}{stem}_vxc_fxc_unpolar")
    lines += [
        f"{INDENT}implicit none",
        f"{INDENT}private",
        "",
        f"{INDENT}public :: trc_xc_functional_t, xc_eval_point, xc_eval_point_polar, xc_functional_by_name",
        f"{INDENT}public :: xc_eval_point_fxc",
        f"{INDENT}public :: xc_eval_point_mgga, xc_eval_point_mgga_polar",
        f"{INDENT}public :: xc_is_double_hybrid",
        f"{INDENT}public :: XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_MAX_KERNELS",
        f"{INDENT}public :: XC_FUNCTIONAL_NAMES, XC_AVAILABLE_FUNCTIONALS",
        "",
        f"{INDENT}integer, parameter :: XC_FAMILY_LDA = 1, XC_FAMILY_GGA = 2, XC_FAMILY_MGGA = 3",
        f"{INDENT}integer, parameter :: XC_MAX_KERNELS = 8",
        "",
    ]
    for i, stem in enumerate(available_kernels, start=1):
        lines.append(f"{INDENT}integer, parameter :: {kid_name[stem]} = {i}")
    lines.append("")
    names_sorted = sorted({e.name for e in entries} | {a for e in entries for a in e.aliases})
    namelen = max(len(n) for n in names_sorted)
    lines.append(f"{INDENT}!! Every selectable name, available or `missing`: the CLI's --func")
    lines.append(f"{INDENT}!! help and trc_capi's unknown-functional message are built from this")
    lines.append(f"{INDENT}!! rather than keeping their own copy of the list.")
    lines.append(f"{INDENT}integer, parameter :: XC_FUNCTIONAL_COUNT = {len(names_sorted)}")
    lines.append(f"{INDENT}character(len={namelen}), parameter :: XC_FUNCTIONAL_NAMES(XC_FUNCTIONAL_COUNT) = [ &")
    for i in range(0, len(names_sorted), 6):
        chunk = names_sorted[i:i + 6]
        items = ", ".join(f'"{n:<{namelen}}"' for n in chunk)
        suffix = "]" if i + 6 >= len(names_sorted) else ", &"
        lines.append(f"{INDENT}{INDENT}{items}{suffix}")
    available_names = sorted({e.name for e in entries if e.missing is None} |
                              {a for e in entries if e.missing is None for a in e.aliases})
    joined = ", ".join(available_names)
    lines.append(f'{INDENT}character(len=*), parameter :: XC_AVAILABLE_FUNCTIONALS = "{joined}"')
    lines.append("")

    lines.append(f"{INDENT}type :: trc_xc_functional_t")
    lines.append(f'{INDENT}{INDENT}character(len=24) :: name = ""')
    lines.append(f"{INDENT}{INDENT}integer :: family = XC_FAMILY_LDA")
    lines.append(f"{INDENT}{INDENT}integer :: nk = 0")
    lines.append(f"{INDENT}{INDENT}integer :: kid(XC_MAX_KERNELS) = 0")
    lines.append(f"{INDENT}{INDENT}real(dp) :: coef(XC_MAX_KERNELS) = 0.0_dp")
    lines.append(f"{INDENT}{INDENT}real(dp) :: exx = 0.0_dp   !! fraction of exact exchange the caller must add")
    lines.append(f"{INDENT}{INDENT}!")
    lines.append(f"{INDENT}{INDENT}! Range-separated hybrids. `exx` keeps its meaning above unchanged:")
    lines.append(f"{INDENT}{INDENT}! the fraction of FULL-range exact exchange, 1/r12. These four are")
    lines.append(f"{INDENT}{INDENT}! zero for every functional that is not one.")
    lines.append(f"{INDENT}{INDENT}!")
    lines.append(f"{INDENT}{INDENT}real(dp) :: exx_lr = 0.0_dp  !! extra fraction of LONG-RANGE exact exchange, erf(omega r)/r")
    lines.append(f"{INDENT}{INDENT}real(dp) :: omega  = 0.0_dp  !! range-separation parameter, bohr^-1; 0 = none")
    lines.append(f"{INDENT}{INDENT}real(dp) :: nlc_b  = 0.0_dp  !! VV10 b; 0 = no nonlocal correlation")
    lines.append(f"{INDENT}{INDENT}real(dp) :: nlc_c  = 0.0_dp  !! VV10 C")
    lines.append(f"{INDENT}{INDENT}!")
    lines.append(f"{INDENT}{INDENT}! Double hybrids: the RI-MP2 fraction on top of the semilocal/exact-")
    lines.append(f"{INDENT}{INDENT}! exchange composition above, split opposite-/same-spin so a")
    lines.append(f"{INDENT}{INDENT}! spin-component-scaled (DSD-type) double hybrid can set them")
    lines.append(f"{INDENT}{INDENT}! unequally; both zero for every functional that is not one, and")
    lines.append(f"{INDENT}{INDENT}! equal (c_os = c_ss = c) for every double hybrid terco has today.")
    lines.append(f"{INDENT}{INDENT}! trc_capi.F90's trc_run_rimp2/trc_energy read these, never a kernel.")
    lines.append(f"{INDENT}{INDENT}!")
    lines.append(f"{INDENT}{INDENT}! WHOEVER ADDS THE FIRST DSD-TYPE ENTRY (pt2_cos /= pt2_css) MUST")
    lines.append(f"{INDENT}{INDENT}! ALSO ADD A TEST THAT EXERCISES UNEQUAL WEIGHTS: every double")
    lines.append(f"{INDENT}{INDENT}! hybrid and every test of its weighted total today has pt2_cos ==")
    lines.append(f"{INDENT}{INDENT}! pt2_css, so a cos/css coefficient swap anywhere in the chain")
    lines.append(f"{INDENT}{INDENT}! (here, trc_capi's trc_energy, or a test's own oracle) is")
    lines.append(f"{INDENT}{INDENT}! undetectable -- see context_scf.F90's check_double_hybrid for")
    lines.append(f"{INDENT}{INDENT}! where that gap is called out and why literal-coefficient")
    lines.append(f"{INDENT}{INDENT}! assertions alone do not close it.")
    lines.append(f"{INDENT}{INDENT}!")
    lines.append(f"{INDENT}{INDENT}real(dp) :: pt2_cos = 0.0_dp  !! opposite-spin RI-MP2 fraction")
    lines.append(f"{INDENT}{INDENT}real(dp) :: pt2_css = 0.0_dp  !! same-spin RI-MP2 fraction")
    lines.append(f"{INDENT}{INDENT}character(len=160) :: missing = \"\"")
    lines.append(f"{INDENT}{INDENT}!! Empty when the functional runs; otherwise exactly what it needs")
    lines.append(f"{INDENT}{INDENT}!! that terco does not have, e.g. \"needs R2SCANL_C: needs the")
    lines.append(f'{INDENT}{INDENT}!! density Laplacian, ...\". Selecting it by name still succeeds;')
    lines.append(f"{INDENT}{INDENT}!! running it is what the caller must refuse.")
    lines.append(f"{INDENT}end type trc_xc_functional_t")
    lines.append("")
    lines.append("contains")
    lines.append("")

    # xc_eval_point / xc_eval_point_polar
    lines.append(f"{INDENT}!")
    lines.append(f"{INDENT}! Energy per particle and first derivatives of the composed functional at")
    lines.append(f"{INDENT}! one point of a closed-shell density. `sigma` is |grad rho|^2 and is")
    lines.append(f"{INDENT}! ignored by an LDA kernel.")
    lines.append(f"{INDENT}!")
    lines.append(f"{INDENT}pure subroutine xc_eval_point(nk, kid, coef, rho, sigma, eps, vrho, vsigma)")
    lines.append(f"{INDENT}{INDENT}!$acc routine seq")
    lines.extend(_kid_decls(["rho", "sigma"]))
    lines.append(f"{INDENT}{INDENT}real(dp), intent(out) :: eps, vrho, vsigma")
    lines.append(f"{INDENT}{INDENT}real(dp) :: e, vr, vs")
    lines.append(f"{INDENT}{INDENT}integer :: k")
    lines.append("")
    lines.extend(_zero(2, ["eps", "vrho", "vsigma"]))
    lines.append(f"{INDENT}{INDENT}do k = 1, nk")
    lines.extend(_zero(3, ["vs"]))
    lines.append(f"{INDENT}{INDENT}{INDENT}select case (kid(k))")
    for stem in semilocal:
        call = f"call {stem}_exc_vxc_unpolar(rho, sigma, e, vr, vs)" if \
            _family_of(stem) == "gga" else f"call {stem}_exc_vxc_unpolar(rho, e, vr)"
        lines.extend(_case(kid_name[stem], call))
    lines.extend(_case_default(["e", "vr"]))
    lines.append(f"{INDENT}{INDENT}{INDENT}end select")
    lines.append(f"{INDENT}{INDENT}{INDENT}eps = eps + coef(k)*e")
    lines.append(f"{INDENT}{INDENT}{INDENT}vrho = vrho + coef(k)*vr")
    lines.append(f"{INDENT}{INDENT}{INDENT}vsigma = vsigma + coef(k)*vs")
    lines.append(f"{INDENT}{INDENT}end do")
    lines.append(f"{INDENT}end subroutine xc_eval_point")
    lines.append("")

    # xc_eval_point_fxc: the restricted-spin second derivatives CPKS needs
    # (v2rho2, v2rhosigma, v2sigma2), plus vsigma -- the GGA response term
    # needs both the kernel's second derivative AND its first, since the
    # grad-rho factor the energy path's own vsigma multiplies is itself
    # unperturbed; only the coefficient in front of it (vsigma) and the
    # trial density's own gradient respond. Same `_vxc_fxc_unpolar` routine
    # check_xc.F90 already validates against libxc (tol_f = 1e-10) -- no
    # new kernel arithmetic, only a new way of calling what is there.
    lines.append(f"{INDENT}!")
    lines.append(f"{INDENT}! The second derivatives of the composed functional at one point of a")
    lines.append(f"{INDENT}! closed-shell density, plus vsigma -- the GGA response term needs both")
    lines.append(f"{INDENT}! the kernel's second derivative AND its first derivative, since the")
    lines.append(f"{INDENT}! grad-rho factor the energy path's vsigma multiplies is itself")
    lines.append(f"{INDENT}! unperturbed; only vsigma's own value and the trial density's gradient")
    lines.append(f"{INDENT}! respond. Same `_vxc_fxc_unpolar` routine check_xc.F90 already validates")
    lines.append(f"{INDENT}! against libxc (tol_f = 1e-10 relative, see that program's header) --")
    lines.append(f"{INDENT}! this calls what is there, with no new kernel arithmetic. An LDA kernel")
    lines.append(f"{INDENT}! leaves vsigma, v2rhosigma and v2sigma2 at zero.")
    lines.append(f"{INDENT}!")
    lines.append(f"{INDENT}pure subroutine xc_eval_point_fxc(nk, kid, coef, rho, sigma, vsigma, v2rho2, v2rhosigma, v2sigma2)")
    lines.append(f"{INDENT}{INDENT}!$acc routine seq")
    lines.extend(_kid_decls(["rho", "sigma"]))
    lines.append(f"{INDENT}{INDENT}real(dp), intent(out) :: vsigma, v2rho2, v2rhosigma, v2sigma2")
    lines.append(f"{INDENT}{INDENT}real(dp) :: vr, vs, f2rr, f2rs, f2ss")
    lines.append(f"{INDENT}{INDENT}integer :: k")
    lines.append("")
    lines.extend(_zero(2, ["vsigma", "v2rho2", "v2rhosigma", "v2sigma2"]))
    lines.append(f"{INDENT}{INDENT}do k = 1, nk")
    lines.extend(_zero(3, ["vs", "f2rs", "f2ss"]))
    lines.append(f"{INDENT}{INDENT}{INDENT}select case (kid(k))")
    for stem in semilocal:
        if _family_of(stem) == "gga":
            call = f"call {stem}_vxc_fxc_unpolar(rho, sigma, vr, vs, f2rr, f2rs, f2ss)"
        else:
            call = f"call {stem}_vxc_fxc_unpolar(rho, vr, f2rr)"
        lines.extend(_case(kid_name[stem], call))
    lines.extend(_case_default(["f2rr"]))
    lines.append(f"{INDENT}{INDENT}{INDENT}end select")
    lines.append(f"{INDENT}{INDENT}{INDENT}vsigma = vsigma + coef(k)*vs")
    lines.append(f"{INDENT}{INDENT}{INDENT}v2rho2 = v2rho2 + coef(k)*f2rr")
    lines.append(f"{INDENT}{INDENT}{INDENT}v2rhosigma = v2rhosigma + coef(k)*f2rs")
    lines.append(f"{INDENT}{INDENT}{INDENT}v2sigma2 = v2sigma2 + coef(k)*f2ss")
    lines.append(f"{INDENT}{INDENT}end do")
    lines.append(f"{INDENT}end subroutine xc_eval_point_fxc")
    lines.append("")

    lines.append(f"{INDENT}!")
    lines.append(f"{INDENT}! The same for a spin-polarised density: two densities, three gradient")
    lines.append(f"{INDENT}! invariants (aa, ab, bb), and the derivatives with respect to each.")
    lines.append(f"{INDENT}! `eps` is still per particle of the TOTAL density.")
    lines.append(f"{INDENT}!")
    sig = ("pure subroutine xc_eval_point_polar(nk, kid, coef, rho_a, rho_b, s_aa, s_ab, s_bb, &")
    lines.append(f"{INDENT}{sig}")
    lines.append(f"{INDENT}{INDENT}{INDENT}{INDENT}eps, vrho_a, vrho_b, vs_aa, vs_ab, vs_bb)")
    lines.append(f"{INDENT}{INDENT}!$acc routine seq")
    lines.extend(_kid_decls(["rho_a", "rho_b", "s_aa", "s_ab", "s_bb"]))
    lines.append(f"{INDENT}{INDENT}real(dp), intent(out) :: eps, vrho_a, vrho_b, vs_aa, vs_ab, vs_bb")
    lines.append(f"{INDENT}{INDENT}real(dp) :: e, va, vb, saa, sab, sbb")
    lines.append(f"{INDENT}{INDENT}integer :: k")
    lines.append("")
    lines.extend(_zero(2, ["eps", "vrho_a", "vrho_b", "vs_aa", "vs_ab", "vs_bb"]))
    lines.append(f"{INDENT}{INDENT}do k = 1, nk")
    lines.extend(_zero(3, ["saa", "sab", "sbb"]))
    lines.append(f"{INDENT}{INDENT}{INDENT}select case (kid(k))")
    for stem in semilocal:
        if _family_of(stem) == "gga":
            call = f"call {stem}_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)"
        else:
            call = f"call {stem}_exc_vxc_polar(rho_a, rho_b, e, va, vb)"
        lines.extend(_case(kid_name[stem], call, wrapped=True))
    lines.extend(_case_default(["e", "va", "vb"]))
    lines.append(f"{INDENT}{INDENT}{INDENT}end select")
    lines.append(f"{INDENT}{INDENT}{INDENT}eps = eps + coef(k)*e")
    lines.append(f"{INDENT}{INDENT}{INDENT}vrho_a = vrho_a + coef(k)*va")
    lines.append(f"{INDENT}{INDENT}{INDENT}vrho_b = vrho_b + coef(k)*vb")
    lines.append(f"{INDENT}{INDENT}{INDENT}vs_aa = vs_aa + coef(k)*saa")
    lines.append(f"{INDENT}{INDENT}{INDENT}vs_ab = vs_ab + coef(k)*sab")
    lines.append(f"{INDENT}{INDENT}{INDENT}vs_bb = vs_bb + coef(k)*sbb")
    lines.append(f"{INDENT}{INDENT}end do")
    lines.append(f"{INDENT}end subroutine xc_eval_point_polar")
    lines.append("")

    lines.extend(emit_mgga_eval(mgga_dispatch, kid_name))

    # xc_functional_by_name
    lines.append(f"{INDENT}!")
    lines.append(f"{INDENT}! libxc's compositions where terco has one; ExchCXX's own table for")
    lines.append(f"{INDENT}! everything else. Names are case-insensitive. A functional this build")
    lines.append(f"{INDENT}! cannot run yet is still returned, ok = .true., with `missing` set.")
    lines.append(f"{INDENT}!")
    lines.append(f"{INDENT}function xc_functional_by_name(name, f) result(ok)")
    lines.append(f"{INDENT}{INDENT}character(len=*), intent(in) :: name")
    lines.append(f"{INDENT}{INDENT}type(trc_xc_functional_t), intent(out) :: f")
    lines.append(f"{INDENT}{INDENT}logical :: ok")
    lines.append(f"{INDENT}{INDENT}character(len=24) :: key")
    lines.append(f"{INDENT}{INDENT}integer :: i")
    lines.append("")
    lines.append(f"{INDENT}{INDENT}key = name")
    lines.append(f"{INDENT}{INDENT}do i = 1, len_trim(key)")
    lines.append(f'{INDENT}{INDENT}{INDENT}if (key(i:i) >= "A" .and. key(i:i) <= "Z") key(i:i) = achar(iachar(key(i:i)) + 32)')
    lines.append(f"{INDENT}{INDENT}end do")
    lines.append(f"{INDENT}{INDENT}f%name = key")
    lines.append(f"{INDENT}{INDENT}ok = .true.")
    lines.append(f"{INDENT}{INDENT}select case (trim(key))")
    for e in entries:
        case_names = ", ".join(f'"{n}"' for n in (e.name, *e.aliases))
        for ln in wrap(f"{INDENT}{INDENT}case ({case_names})", INDENT*2):
            lines.append(ln)
        for c in e.comment:
            lines.append(f"{INDENT}{INDENT}{INDENT}! {c}")
        if e.missing:
            # No `set` call: a functional short of a kernel gets nk = 0 (the
            # type's own default), not a partial composition that would
            # quietly run the wrong physics if a caller ever skipped the
            # `missing` check. Available vs not is a sharp line here.
            esc = e.missing.replace('"', "'")
            lines.append(f'{INDENT}{INDENT}{INDENT}f%missing = "{esc}"')
        else:
            family_param = {"lda": "XC_FAMILY_LDA", "gga": "XC_FAMILY_GGA", "mgga": "XC_FAMILY_MGGA"}[e.family]
            kids = ", ".join(kid_name[s] for _, s in e.kernels)
            coefs = ", ".join(fortran_real(c) for c, _ in e.kernels)
            for ln in wrap(f"{INDENT}{INDENT}{INDENT}call set(f, {family_param}, [{kids}], [{coefs}])", INDENT*3):
                lines.append(ln)
            if e.exx:
                lines.append(f"{INDENT}{INDENT}{INDENT}f%exx = {fortran_real(e.exx)}")
            if e.exx_lr:
                lines.append(f"{INDENT}{INDENT}{INDENT}f%exx_lr = {fortran_real(e.exx_lr)}")
            if e.omega:
                lines.append(f"{INDENT}{INDENT}{INDENT}f%omega = {fortran_real(e.omega)}")
            if e.nlc_b:
                lines.append(f"{INDENT}{INDENT}{INDENT}f%nlc_b = {fortran_real(e.nlc_b)}")
            if e.nlc_c:
                lines.append(f"{INDENT}{INDENT}{INDENT}f%nlc_c = {fortran_real(e.nlc_c)}")
            if e.pt2_cos:
                lines.append(f"{INDENT}{INDENT}{INDENT}f%pt2_cos = {fortran_real(e.pt2_cos)}")
            if e.pt2_css:
                lines.append(f"{INDENT}{INDENT}{INDENT}f%pt2_css = {fortran_real(e.pt2_css)}")
    lines.append(f"{INDENT}{INDENT}case default")
    lines.append(f"{INDENT}{INDENT}{INDENT}ok = .false.")
    lines.append(f"{INDENT}{INDENT}end select")
    lines.append(f"{INDENT}contains")
    lines.append(f"{INDENT}{INDENT}subroutine set(f, family, kid, coef)")
    lines.append(f"{INDENT}{INDENT}{INDENT}type(trc_xc_functional_t), intent(inout) :: f")
    lines.append(f"{INDENT}{INDENT}{INDENT}integer, intent(in) :: family")
    lines.append(f"{INDENT}{INDENT}{INDENT}integer, intent(in) :: kid(:)")
    lines.append(f"{INDENT}{INDENT}{INDENT}real(dp), intent(in) :: coef(:)")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%family = family")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%nk = size(kid)")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%kid = 0")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%coef = 0.0_dp")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%kid(1:f%nk) = kid")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%coef(1:f%nk) = coef")
    lines.append(f"{INDENT}{INDENT}end subroutine set")
    lines.append(f"{INDENT}end function xc_functional_by_name")
    lines.append("")
    lines.append(f"{INDENT}!")
    lines.append(f"{INDENT}! Whether `f` needs an RI-MP2 pass on top of its SCF: either PT2")
    lines.append(f"{INDENT}! fraction nonzero. The one place this is asked, rather than reading")
    lines.append(f"{INDENT}! pt2_cos/pt2_css directly, so a future third PT2 field cannot leave a")
    lines.append(f"{INDENT}! caller checking only the first two of three.")
    lines.append(f"{INDENT}!")
    lines.append(f"{INDENT}pure function xc_is_double_hybrid(f) result(is_dh)")
    lines.append(f"{INDENT}{INDENT}type(trc_xc_functional_t), intent(in) :: f")
    lines.append(f"{INDENT}{INDENT}logical :: is_dh")
    lines.append(f"{INDENT}{INDENT}is_dh = f%pt2_cos /= 0.0_dp .or. f%pt2_css /= 0.0_dp")
    lines.append(f"{INDENT}end function xc_is_double_hybrid")
    lines.append("")
    lines.append("end module trc_xc_functional")
    lines.append("")
    return "\n".join(lines)


_FAMILY_CACHE: dict[str, str] = {}


def emit_mgga_eval(dispatch: list[str], kid_name: dict[str, str]) -> list[str]:
    """xc_eval_point_mgga and xc_eval_point_mgga_polar, the meta-GGA point
    evaluators. An LDA or GGA kernel among them ignores tau and contributes
    no vtau, the way libxc evaluates a mixed composition kernel by kernel."""
    I = INDENT
    out = [
        f"{I}!",
        f"{I}! Energy per particle and first derivatives of a composed meta-GGA at one",
        f"{I}! point of a closed-shell density: `sigma` is |grad rho|^2 and `tau` the",
        f"{I}! kinetic-energy density 1/2 sum_i |grad psi_i|^2 of the total density,",
        f"{I}! libxc's and pyscf's convention. `vtau` is d(rho eps)/d tau. An LDA or",
        f"{I}! GGA kernel in the composition ignores tau.",
        f"{I}!",
        f"{I}pure subroutine xc_eval_point_mgga(nk, kid, coef, rho, sigma, tau, eps, vrho, vsigma, vtau)",
        f"{I}{I}!$acc routine seq",
    ]
    out += _kid_decls(["rho", "sigma", "tau"])
    out += [
        f"{I}{I}real(dp), intent(out) :: eps, vrho, vsigma, vtau",
        f"{I}{I}real(dp) :: e, vr, vs, vt",
        f"{I}{I}integer :: k",
        "",
    ]
    out += _zero(2, ["eps", "vrho", "vsigma", "vtau"])
    out.append(f"{I}{I}do k = 1, nk")
    out += _zero(3, ["vs", "vt"])
    out.append(f"{I}{I}{I}select case (kid(k))")
    for stem in dispatch:
        fam = _family_of(stem)
        if fam == "mgga":
            call = f"call {stem}_exc_vxc_unpolar(rho, sigma, tau, e, vr, vs, vt)"
        elif fam == "gga":
            call = f"call {stem}_exc_vxc_unpolar(rho, sigma, e, vr, vs)"
        else:
            call = f"call {stem}_exc_vxc_unpolar(rho, e, vr)"
        out += _case(kid_name[stem], call, wrapped=True)
    out += _case_default(["e", "vr"])
    out += [
        f"{I}{I}{I}end select",
        f"{I}{I}{I}eps = eps + coef(k)*e",
        f"{I}{I}{I}vrho = vrho + coef(k)*vr",
        f"{I}{I}{I}vsigma = vsigma + coef(k)*vs",
        f"{I}{I}{I}vtau = vtau + coef(k)*vt",
        f"{I}{I}end do",
        f"{I}end subroutine xc_eval_point_mgga",
        "",
        f"{I}!",
        f"{I}! The same for a spin-polarised density: per-spin densities, the three",
        f"{I}! gradient invariants and per-spin tau_s = 1/2 sum_i |grad psi_is|^2.",
        f"{I}!",
        f"{I}pure subroutine xc_eval_point_mgga_polar(nk, kid, coef, rho_a, rho_b, s_aa, s_ab, s_bb, &",
        f"{I}{I}{I}{I}tau_a, tau_b, eps, vrho_a, vrho_b, vs_aa, vs_ab, vs_bb, vtau_a, vtau_b)",
        f"{I}{I}!$acc routine seq",
    ]
    out += _kid_decls(["rho_a", "rho_b", "s_aa", "s_ab", "s_bb", "tau_a", "tau_b"])
    out += [
        f"{I}{I}real(dp), intent(out) :: eps, vrho_a, vrho_b, vs_aa, vs_ab, vs_bb, vtau_a, vtau_b",
        f"{I}{I}real(dp) :: e, va, vb, saa, sab, sbb, ta, tb",
        f"{I}{I}integer :: k",
        "",
    ]
    out += _zero(2, ["eps", "vrho_a", "vrho_b", "vs_aa", "vs_ab", "vs_bb", "vtau_a", "vtau_b"])
    out.append(f"{I}{I}do k = 1, nk")
    out += _zero(3, ["saa", "sab", "sbb", "ta", "tb"])
    out.append(f"{I}{I}{I}select case (kid(k))")
    for stem in dispatch:
        fam = _family_of(stem)
        if fam == "mgga":
            call = (f"call {stem}_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, tau_a, tau_b, "
                    f"e, va, vb, saa, sab, sbb, ta, tb)")
        elif fam == "gga":
            call = f"call {stem}_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)"
        else:
            call = f"call {stem}_exc_vxc_polar(rho_a, rho_b, e, va, vb)"
        out += _case(kid_name[stem], call, wrapped=True)
    out += _case_default(["e", "va", "vb"])
    out += [
        f"{I}{I}{I}end select",
        f"{I}{I}{I}eps = eps + coef(k)*e",
        f"{I}{I}{I}vrho_a = vrho_a + coef(k)*va",
        f"{I}{I}{I}vrho_b = vrho_b + coef(k)*vb",
        f"{I}{I}{I}vs_aa = vs_aa + coef(k)*saa",
        f"{I}{I}{I}vs_ab = vs_ab + coef(k)*sab",
        f"{I}{I}{I}vs_bb = vs_bb + coef(k)*sbb",
        f"{I}{I}{I}vtau_a = vtau_a + coef(k)*ta",
        f"{I}{I}{I}vtau_b = vtau_b + coef(k)*tb",
        f"{I}{I}end do",
        f"{I}end subroutine xc_eval_point_mgga_polar",
        "",
    ]
    return out


# One statement per line and one kind of declaration per line: the style
# fortitude's S-rules hold src/xc/ to (see fpm.toml).
def _kid_decls(scalars_in: list[str]) -> list[str]:
    return [f"{INDENT}{INDENT}integer, intent(in) :: nk",
            f"{INDENT}{INDENT}integer, intent(in) :: kid(XC_MAX_KERNELS)",
            f"{INDENT}{INDENT}real(dp), intent(in) :: {', '.join(scalars_in)}",
            f"{INDENT}{INDENT}real(dp), intent(in) :: coef(XC_MAX_KERNELS)"]


def _zero(depth: int, names: list[str]) -> list[str]:
    return [f"{INDENT * depth}{n} = 0.0_dp" for n in names]


def _case(kid: str, call: str, wrapped: bool = False) -> list[str]:
    """`case (kid)` and its call, at the select's own indent. A wrapped call
    breaks as if it still followed `case (kid); ` on the same line."""
    ind = INDENT * 3
    out = wrap(f"{ind}case ({kid}); {call}", ind) if wrapped else [f"{ind}case ({kid}); {call}"]
    head, first = out[0].split("; ", 1)
    return [head, ind + first] + out[1:]


def _case_default(names: list[str]) -> list[str]:
    return [f"{INDENT * 3}case default"] + _zero(3, names)


def _family_of(stem: str) -> str:
    return _FAMILY_CACHE[stem]


# --------------------------------------------------------------------------
# Emitting docs/functionals.md.
# --------------------------------------------------------------------------


NUMBER_WORDS = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five"}


def emit_docs(entries: list[Entry], notes: list[str], nkernel_translated: int,
              nkernel_refused: int, nkernel_disagree: int) -> str:
    avail = [e for e in entries if e.missing is None]
    miss = [e for e in entries if e.missing is not None]
    lines = [
        "# terco's XC functionals",
        "",
        "GENERATED by `tools/gen_xc_registry.py` -- do not edit.",
        "",
        f"Every functional [ExchCXX](https://github.com/wavefunction91/ExchCXX) "
        f"defines, plus terco's own local entries, is selectable by name through "
        f"`xc_functional_by_name`. {len(avail)} run; {len(miss)} are registered but "
        f"refuse to run, each for a stated reason.",
        "",
        f"Kernels: {nkernel_translated} translated by `tools/exchcxx_to_fortran.py`, "
        f"of which {nkernel_translated - nkernel_disagree} agree with libxc at "
        f"check_xc's tolerances and {nkernel_disagree} do not (listed under their "
        f"functionals below); {nkernel_refused} refused outright.",
        "",
    ]
    dh = [e for e in avail if e.pt2_cos or e.pt2_css]
    if dh:
        names = ", ".join(e.name for e in dh)
        lines += [
            f"{NUMBER_WORDS.get(len(dh), str(len(dh)))} of the {len(avail)} are double hybrids ({names}): the table",
            "below lists their SCF composition only (a_x exact exchange plus a GGA",
            "exchange/correlation mix), not the RI-MP2 fraction on top, which is never a",
            "kernel -- see `trc_xc_functional_t`'s own `pt2_cos`/`pt2_css` and",
            "`trc_capi.F90`'s `trc_run_rimp2`/`trc_energy`.",
            "",
        ]
    lines += [
        "## Available",
        "",
        "| name | aliases | family |",
        "| --- | --- | --- |",
    ]
    # The double hybrids first, in their own order, then everything else.
    for e in dh + sorted((e for e in avail if e not in dh), key=lambda e: e.name):
        lines.append(f"| `{e.name}` | {', '.join(e.aliases) or '--'} | {e.family} |")
    lines += ["", "## Missing", "", "| name | reason |", "| --- | --- |"]
    for e in sorted(miss, key=lambda e: e.name):
        lines.append(f"| `{e.name}` | {e.missing} |")
    if notes:
        lines += ["", "## Notes", ""]
        lines.extend(notes)
    lines.append("")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exchcxx", required=True, type=Path)
    ap.add_argument("-o", "--output", required=True, type=Path)
    ap.add_argument("--docs", required=True, type=Path)
    args = ap.parse_args(argv)

    translated, refused = xc_kernel_discovery.discover(args.exchcxx)
    # LOCAL_KERNELS joins the same family cache and id numbering as a
    # translated stem (see its own comment above `LOCAL`), but is kept out
    # of `translated` itself: emit_docs' "N translated by
    # exchcxx_to_fortran.py" count must stay exactly what came from ExchCXX.
    for stem, fam in translated + LOCAL_KERNELS:
        _FAMILY_CACHE[stem] = fam
    available_kernels = [s for s, _ in translated if s not in xc_kernel_discovery.DISAGREES_WITH_LIBXC]
    available_kernels += [s for s, _ in LOCAL_KERNELS]

    # LDA/GGA ids first, meta-GGA after: adding the meta-GGA family must not
    # renumber the LDA/GGA dispatch.
    ordered = sorted(st for st in available_kernels if _FAMILY_CACHE[st] != "mgga") + \
        sorted(st for st in available_kernels if _FAMILY_CACHE[st] == "mgga")

    entries, notes = build_registry(args.exchcxx)
    args.output.write_text(emit_fortran(entries, ordered))
    args.docs.write_text(emit_docs(entries, notes, len(translated), len(refused),
                                    len(xc_kernel_discovery.DISAGREES_WITH_LIBXC)))
    navail = sum(1 for e in entries if e.missing is None)
    print(f"{args.output}: {len(entries)} functionals ({navail} available), "
          f"{len(available_kernels)} kernels")
    print(f"{args.docs}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
