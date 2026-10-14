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

Plus two small local tables, kept here rather than derived because nothing
upstream states them: LOCAL (terco-only entries -- svwn_rpa, b3lyp/b3lyp5's
libxc-not-ExchCXX split, wb97x-v's VV10 coefficients -- carried over from
the hand-written registry this replaces, unchanged) and OVERRIDES (one
entry: ExchCXX's own revB3LYP composition names Kernel::VWN5, which its own
builtin_interface.cxx says is the VWN_RPA formula, but matches libxc's
actual revB3LYP/revB3LYP5 only through the plain VWN formula -- checked
numerically against pyscf; see the commit that added this generator).

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

# --------------------------------------------------------------------------
# terco-only entries, unchanged from the hand-written registry. `kernels`
# is [(coef, stem), ...]; `exx`/`exx_lr`/`omega`/`nlc_b`/`nlc_c` default to
# 0. `comment` is emitted verbatim above the case, FORD style.
# --------------------------------------------------------------------------


@dataclass
class Entry:
    name: str
    aliases: tuple[str, ...]
    family: str  # "lda" | "gga"
    kernels: list[tuple[float, str]]
    exx: float = 0.0
    exx_lr: float = 0.0
    omega: float = 0.0
    nlc_b: float = 0.0
    nlc_c: float = 0.0
    comment: list[str] = field(default_factory=list)
    missing: str | None = None
    source: str = "local"  # local | exchcxx


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

        exx = exx_lr = omega = 0.0
        if fn.hyb:
            alpha, beta, omega = fn.hyb
            exx = alpha + beta
            exx_lr = -beta

        if problems:
            # Group enums by identical reason so "needs TPSS_X, TPSS_C:
            # meta-GGA kernels are not translated yet" prints once.
            by_reason: dict[str, list[str]] = {}
            for enum, reason in problems.items():
                by_reason.setdefault(reason, []).append(enum)
            missing = "; ".join(f"needs {', '.join(sorted(enums))}: {reason}"
                                 for reason, enums in by_reason.items())
            entries.append(Entry(lname, (), family, kernel_stems, exx, exx_lr, omega,
                                  missing=missing, source="exchcxx"))
        else:
            entries.append(Entry(lname, (), family, kernel_stems, exx, exx_lr, omega, source="exchcxx"))

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
        "module trc_xc_functional",
        f"{INDENT}use trc_boys, only: dp",
    ]
    for stem in available_kernels:
        lines.append(f"{INDENT}use trc_xc_{stem}, only: {stem}_exc_vxc_unpolar, {stem}_exc_vxc_polar")
    lines += [
        f"{INDENT}implicit none",
        f"{INDENT}private",
        "",
        f"{INDENT}public :: trc_xc_functional_t, xc_eval_point, xc_eval_point_polar, xc_functional_by_name",
        f"{INDENT}public :: XC_FAMILY_LDA, XC_FAMILY_GGA, XC_MAX_KERNELS",
        f"{INDENT}public :: XC_FUNCTIONAL_NAMES, XC_AVAILABLE_FUNCTIONALS",
        "",
        f"{INDENT}integer, parameter :: XC_FAMILY_LDA = 1, XC_FAMILY_GGA = 2",
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
        suffix = " ]" if i + 6 >= len(names_sorted) else ", &"
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
    lines.append(f"{INDENT}{INDENT}real(dp) :: exx_lr = 0.0_dp !! extra fraction of LONG-RANGE exact exchange, erf(omega r)/r")
    lines.append(f"{INDENT}{INDENT}real(dp) :: omega  = 0.0_dp !! range-separation parameter, bohr^-1; 0 = none")
    lines.append(f"{INDENT}{INDENT}real(dp) :: nlc_b  = 0.0_dp !! VV10 b; 0 = no nonlocal correlation")
    lines.append(f"{INDENT}{INDENT}real(dp) :: nlc_c  = 0.0_dp !! VV10 C")
    lines.append(f"{INDENT}{INDENT}!! Empty when the functional runs; otherwise exactly what it needs")
    lines.append(f"{INDENT}{INDENT}!! that terco does not have, e.g. \"needs TPSS_X: meta-GGA kernels")
    lines.append(f'{INDENT}{INDENT}!! are not translated yet\". Selecting it by name still succeeds;')
    lines.append(f"{INDENT}{INDENT}!! running it is what the caller must refuse.")
    lines.append(f"{INDENT}{INDENT}character(len=160) :: missing = \"\"")
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
    lines.append(f"{INDENT}{INDENT}integer, intent(in) :: nk, kid(XC_MAX_KERNELS)")
    lines.append(f"{INDENT}{INDENT}real(dp), intent(in) :: coef(XC_MAX_KERNELS), rho, sigma")
    lines.append(f"{INDENT}{INDENT}real(dp), intent(out) :: eps, vrho, vsigma")
    lines.append(f"{INDENT}{INDENT}real(dp) :: e, vr, vs")
    lines.append(f"{INDENT}{INDENT}integer :: k")
    lines.append("")
    lines.append(f"{INDENT}{INDENT}eps = 0.0_dp; vrho = 0.0_dp; vsigma = 0.0_dp")
    lines.append(f"{INDENT}{INDENT}do k = 1, nk")
    lines.append(f"{INDENT}{INDENT}{INDENT}vs = 0.0_dp")
    lines.append(f"{INDENT}{INDENT}{INDENT}select case (kid(k))")
    for stem in available_kernels:
        call = f"call {stem}_exc_vxc_unpolar(rho, sigma, e, vr, vs)" if \
            _family_of(stem) == "gga" else f"call {stem}_exc_vxc_unpolar(rho, e, vr)"
        lines.append(f"{INDENT}{INDENT}{INDENT}case ({kid_name[stem]}); {call}")
    lines.append(f"{INDENT}{INDENT}{INDENT}case default; e = 0.0_dp; vr = 0.0_dp")
    lines.append(f"{INDENT}{INDENT}{INDENT}end select")
    lines.append(f"{INDENT}{INDENT}{INDENT}eps = eps + coef(k)*e")
    lines.append(f"{INDENT}{INDENT}{INDENT}vrho = vrho + coef(k)*vr")
    lines.append(f"{INDENT}{INDENT}{INDENT}vsigma = vsigma + coef(k)*vs")
    lines.append(f"{INDENT}{INDENT}end do")
    lines.append(f"{INDENT}end subroutine xc_eval_point")
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
    lines.append(f"{INDENT}{INDENT}integer, intent(in) :: nk, kid(XC_MAX_KERNELS)")
    lines.append(f"{INDENT}{INDENT}real(dp), intent(in) :: coef(XC_MAX_KERNELS), rho_a, rho_b, s_aa, s_ab, s_bb")
    lines.append(f"{INDENT}{INDENT}real(dp), intent(out) :: eps, vrho_a, vrho_b, vs_aa, vs_ab, vs_bb")
    lines.append(f"{INDENT}{INDENT}real(dp) :: e, va, vb, saa, sab, sbb")
    lines.append(f"{INDENT}{INDENT}integer :: k")
    lines.append("")
    lines.append(f"{INDENT}{INDENT}eps = 0.0_dp; vrho_a = 0.0_dp; vrho_b = 0.0_dp")
    lines.append(f"{INDENT}{INDENT}vs_aa = 0.0_dp; vs_ab = 0.0_dp; vs_bb = 0.0_dp")
    lines.append(f"{INDENT}{INDENT}do k = 1, nk")
    lines.append(f"{INDENT}{INDENT}{INDENT}saa = 0.0_dp; sab = 0.0_dp; sbb = 0.0_dp")
    lines.append(f"{INDENT}{INDENT}{INDENT}select case (kid(k))")
    for stem in available_kernels:
        if _family_of(stem) == "gga":
            call = f"call {stem}_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)"
        else:
            call = f"call {stem}_exc_vxc_polar(rho_a, rho_b, e, va, vb)"
        for ln in wrap(f"{INDENT}{INDENT}{INDENT}case ({kid_name[stem]}); {call}", INDENT*3):
            lines.append(ln)
    lines.append(f"{INDENT}{INDENT}{INDENT}case default; e = 0.0_dp; va = 0.0_dp; vb = 0.0_dp")
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
    lines.append(f"{INDENT}{INDENT}{INDENT}if (key(i:i) >= 'A' .and. key(i:i) <= 'Z') key(i:i) = achar(iachar(key(i:i)) + 32)")
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
            family_param = "XC_FAMILY_GGA" if e.family == "gga" else "XC_FAMILY_LDA"
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
    lines.append(f"{INDENT}{INDENT}case default")
    lines.append(f"{INDENT}{INDENT}{INDENT}ok = .false.")
    lines.append(f"{INDENT}{INDENT}end select")
    lines.append(f"{INDENT}contains")
    lines.append(f"{INDENT}{INDENT}subroutine set(f, family, kid, coef)")
    lines.append(f"{INDENT}{INDENT}{INDENT}type(trc_xc_functional_t), intent(inout) :: f")
    lines.append(f"{INDENT}{INDENT}{INDENT}integer, intent(in) :: family, kid(:)")
    lines.append(f"{INDENT}{INDENT}{INDENT}real(dp), intent(in) :: coef(:)")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%family = family")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%nk = size(kid)")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%kid = 0; f%coef = 0.0_dp")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%kid(1:f%nk) = kid")
    lines.append(f"{INDENT}{INDENT}{INDENT}f%coef(1:f%nk) = coef")
    lines.append(f"{INDENT}{INDENT}end subroutine set")
    lines.append(f"{INDENT}end function xc_functional_by_name")
    lines.append("")
    lines.append("end module trc_xc_functional")
    lines.append("")
    return "\n".join(lines)


_FAMILY_CACHE: dict[str, str] = {}


def _family_of(stem: str) -> str:
    return _FAMILY_CACHE[stem]


# --------------------------------------------------------------------------
# Emitting docs/functionals.md.
# --------------------------------------------------------------------------


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
        "## Available",
        "",
        "| name | aliases | family |",
        "| --- | --- | --- |",
    ]
    for e in sorted(avail, key=lambda e: e.name):
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
    for stem, fam in translated:
        _FAMILY_CACHE[stem] = fam
    available_kernels = [s for s, _ in translated if s not in xc_kernel_discovery.DISAGREES_WITH_LIBXC]

    entries, notes = build_registry(args.exchcxx)
    args.output.write_text(emit_fortran(entries, sorted(available_kernels)))
    args.docs.write_text(emit_docs(entries, notes, len(translated), len(refused),
                                    len(xc_kernel_discovery.DISAGREES_WITH_LIBXC)))
    navail = sum(1 for e in entries if e.missing is None)
    print(f"{args.output}: {len(entries)} functionals ({navail} available), "
          f"{len(available_kernels)} kernels")
    print(f"{args.docs}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
