#!/usr/bin/env python3
"""Translate ExchCXX built-in XC kernels into pure Fortran.

ExchCXX (https://github.com/wavefunction91/ExchCXX, MIT) carries, under
`include/exchcxx/impl/builtin/kernels/`, one header per exchange or
correlation kernel. Each is Maple output: a `kernel_traits` struct holding a
few constants and eight member functions

    eval_{exc,exc_vxc,fxc,vxc_fxc}_{unpolar,polar}_impl

whose bodies are straight-line scalar arithmetic -- `const double tN = ...;`
sequences ending in assignments to the reference outputs. No loops, no
branches, no pointers. The only C++ in them is

    safe_math::cbrt(x) and friends      -> a Fortran intrinsic or a helper
    piecewise_functor_3(c, a, b)         -> merge(a, b, c)
    piecewise_functor_5(b, x, c, y, z)   -> merge(x, merge(y, z, c), b)
    constants::m_cbrt_3                  -> a parameter from trc_xc_util
    0.1e1                                -> 0.1e1_dp

`piecewise_functor_*` are ordinary functions taking VALUES, so C++ evaluates
both arms before the call exactly as `merge` does. The translation is
semantically identical, including for any NaN an unselected arm produces.

The Maple output also stores comparisons in `double`s and passes them where
a `bool` is wanted (`const double t2 = rho / 0.2e1 <= dens_tol;`). Those
become `logical` here; the translator infers the type from the expression.

The screening wrappers in `screening_interface.hpp` -- zero the outputs, and
only above `dens_tol` clamp the inputs and call the body -- are emitted from
a template rather than translated, since they are the same for every kernel
of a family.

Meta-GGA kernels are refused. Their wrapper has one more rule (the Fermi
hole curvature bound) and nothing needs them yet.

Usage:
    python tools/exchcxx_to_fortran.py --exchcxx ~/dev/ExchCXX \
        -o src/xc slater_exchange vwn vwn_rpa pbe_x pbe_c b88 lyp
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

KERNEL_DIR = Path("include/exchcxx/impl/builtin/kernels")

# C++ callee -> (Fortran callee, argument reorder)
# A reorder of None keeps the arguments as they are.
FUNC_MAP: dict[str, str] = {
    "safe_math::cbrt": "xc_cbrt",
    "safe_math::sqrt": "sqrt",
    "safe_math::log": "log",
    "safe_math::exp": "exp",
    "safe_math::atan": "atan",
    "safe_math::erf": "erf",
    "safe_math::xc_erfcx": "xc_erfcx",
    "fabs": "abs",
    "square": "xc_square",
    "pow_3_2": "pow_3_2",
    "pow_1_4": "pow_1_4",
    "safe_max": "max",
    "safe_min": "min",
    "enforce_polar_sigma_constraints": "enforce_polar_sigma_constraints",
    "enforce_fermi_hole_curvature": "enforce_fermi_hole_curvature",
}
UTIL_NAMES = {
    "xc_cbrt", "xc_square", "pow_3_2", "pow_1_4", "xc_erfcx",
    "enforce_polar_sigma_constraints", "enforce_fermi_hole_curvature",
}
CONSTANT_NAMES = {
    "m_cbrt_2", "m_cbrt_3", "m_cbrt_4", "m_cbrt_6", "m_cbrt_pi", "m_pi",
    "m_one_ov_pi", "m_pi_sq", "m_cbrt_one_ov_pi", "m_cbrt_pi_sq",
    "x2s", "x_factor_c",
}

LINE_WIDTH = 96
INDENT = "   "


class Unsupported(Exception):
    pass


# --------------------------------------------------------------------------
# Expression translation
# --------------------------------------------------------------------------

def _split_top_level(s: str) -> list[str]:
    """Split on commas at parenthesis depth zero."""
    out, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return [a.strip() for a in out]


_IDENT_CALL = re.compile(r"([A-Za-z_][A-Za-z_0-9:]*)\s*\(")


def _translate_calls(expr: str) -> str:
    """Rewrite every function call, innermost first, by matching parentheses."""
    i = 0
    while True:
        m = _IDENT_CALL.search(expr, i)
        if not m:
            return expr
        name = m.group(1)
        start = m.end()  # just past '('
        depth, j = 1, start
        while depth:
            if j >= len(expr):
                raise Unsupported(f"unbalanced parentheses in: {expr}")
            if expr[j] == "(":
                depth += 1
            elif expr[j] == ")":
                depth -= 1
            j += 1
        inner = _translate_calls(expr[start:j - 1])
        args = _split_top_level(inner)
        rep = _render_call(name, args)
        expr = expr[:m.start()] + rep + expr[j:]
        i = m.start() + len(rep)


def _render_call(name: str, args: list[str]) -> str:
    if name == "piecewise_functor_3":
        c, a, b = args
        return f"merge({a}, {b}, {c})"
    if name == "piecewise_functor_5":
        b, x, c, y, z = args
        return f"merge({x}, merge({y}, {z}, {c}), {b})"
    if name == "safe_math::pow":
        x, e = args
        return f"({x})**({e})"
    if name in FUNC_MAP:
        return f"{FUNC_MAP[name]}({', '.join(args)})"
    raise Unsupported(f"function {name}")


# A floating literal: 0.1e1, 1.0, .5, 1e-24, 2.1544e-43 -- with an optional
# long-double L that C++ constants carry. Not preceded by a word character
# (so t12 stays an identifier) nor followed by one (so 1e1_dp is not doubled).
_FLOAT = re.compile(
    r"(?<![\w.])((?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?|\d+[eE][+-]?\d+)L?(?![\w.])"
)
_INT_DIV = re.compile(r"(?<![\w.])\d+\s*/\s*\d+(?![\w.])")
_BARE_INT = re.compile(r"(?<![\w.])\d+(?![\w.eE])")


def translate_expr(expr: str) -> str:
    expr = expr.strip()
    if _INT_DIV.search(expr):
        raise Unsupported(f"integer division in: {expr}")
    expr = _translate_calls(expr)
    expr = expr.replace("constants::", "")
    expr = re.sub(r"\bX2S\b", "x2s", expr)
    expr = re.sub(r"\bX_FACTOR_C\b", "x_factor_c", expr)
    expr = _FLOAT.sub(r"\1_dp", expr)
    # Bare integers are promoted rather than left to Fortran's integer
    # arithmetic. None appear in the kernels today; this is for the day one does.
    expr = _BARE_INT.sub(lambda m: m.group(0) + ".0_dp", expr)
    expr = expr.replace("!=", "/=").replace("&&", ".and.").replace("||", ".or.")
    if re.search(r"(?<![/<>=])!(?!=)", expr):
        expr = re.sub(r"(?<![/<>=])!(?!=)", ".not.", expr)
    expr = re.sub(r"\s+", " ", expr)
    expr = expr.replace("( ", "(").replace(" )", ")")
    return expr


_CMP = re.compile(r"<=|>=|==|/=|<|>")


def is_logical(expr: str) -> bool:
    """True when the translated expression is a comparison at depth zero."""
    depth = 0
    i = 0
    while i < len(expr):
        ch = expr[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif depth == 0:
            m = _CMP.match(expr, i)
            if m:
                return True
        i += 1
    return False


# --------------------------------------------------------------------------
# Parsing the C++ header
# --------------------------------------------------------------------------

@dataclass
class Function:
    name: str          # eval_exc_unpolar_impl
    inputs: list[str]  # by-value doubles
    outputs: list[str]  # reference doubles
    stmts: list[tuple[str, str, bool]] = field(default_factory=list)  # (lhs, rhs, is_new)


@dataclass
class Kernel:
    name: str
    family: str  # lda | gga
    params: list[tuple[str, str]]  # (name, fortran literal)
    functions: list[Function]


_STATIC = re.compile(r"static\s+constexpr\s+(double|bool)\s+(\w+)\s*=\s*([^;]+);")
_FUNC = re.compile(r"BUILTIN_KERNEL_EVAL_RETURN\s+(eval_\w+_impl)\s*\(([^)]*)\)\s*\{")
_STMT = re.compile(r"^(?:(constexpr|const)\s+(?:double|auto)\s+)?(\w+)\s*=\s*(.+)$", re.S)


def parse_kernel(path: Path) -> Kernel:
    text = path.read_text()
    # Drop block and line comments.
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//[^\n]*", "", text)

    statics: dict[str, str] = {}
    for typ, name, val in _STATIC.findall(text):
        statics[name] = val.strip()

    is_lda = statics.get("is_lda") == "true"
    is_gga = statics.get("is_gga") == "true"
    is_mgga = statics.get("is_mgga") == "true"
    if is_mgga:
        raise Unsupported("meta-GGA kernels are not translated yet")
    if statics.get("is_kedf") == "true" or statics.get("is_epc") == "true":
        raise Unsupported("kinetic-energy and electron-proton kernels are not translated")
    family = "lda" if is_lda else "gga" if is_gga else None
    if family is None:
        raise Unsupported("kernel declares neither LDA nor GGA")

    params: list[tuple[str, str]] = []
    for name, val in statics.items():
        if name.startswith(("is_", "needs_")):
            continue
        if name == "tau_tol":
            val = "1e-20"  # `is_kedf ? 0.0 : 1e-20`, and is_kedf is false here
        params.append((name, translate_expr(val)))

    functions: list[Function] = []
    for m in _FUNC.finditer(text):
        fname = m.group(1)
        sig = m.group(2)
        inputs, outputs = [], []
        for arg in sig.split(","):
            arg = arg.strip()
            if not arg:
                continue
            am = re.match(r"double\s*(&?)\s*(\w+)", arg)
            if not am:
                raise Unsupported(f"argument {arg!r} in {fname}")
            (outputs if am.group(1) else inputs).append(am.group(2))
        # Body: up to the matching brace. The bodies have no nested braces.
        start = m.end()
        depth, j = 1, start
        while depth:
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
            j += 1
        body = text[start:j - 1]
        fn = Function(fname, inputs, outputs)
        for raw in body.split(";"):
            raw = raw.strip()
            if not raw or raw.startswith("(void)"):
                continue  # `(void)(sigma_ab)` silences an unused-argument warning
            sm = _STMT.match(raw)
            if not sm:
                raise Unsupported(f"statement {raw!r} in {fname}")
            fn.stmts.append((sm.group(2), sm.group(3).strip(), sm.group(1) is not None))
        functions.append(fn)

    if len(functions) != 8:
        raise Unsupported(f"expected 8 eval functions, found {len(functions)}")
    return Kernel(path.stem, family, params, functions)


# --------------------------------------------------------------------------
# Emitting Fortran
# --------------------------------------------------------------------------

def wrap(line: str, indent: str) -> list[str]:
    """Break a long statement at spaces, with `&` continuations."""
    if len(line) <= LINE_WIDTH:
        return [line]
    out = []
    cur = line
    cont = indent + INDENT
    while len(cur) > LINE_WIDTH:
        k = cur.rfind(" ", 0, LINE_WIDTH - 2)
        if k <= len(indent):
            k = cur.find(" ", LINE_WIDTH)
            if k < 0:
                break
        out.append(cur[:k] + " &")
        cur = cont + cur[k + 1:]
    out.append(cur)
    return out


def emit_impl(k: Kernel, fn: Function, used: set[str]) -> list[str]:
    sub = f"{k.name}_{fn.name[5:]}"  # drop "eval_"
    ind = INDENT
    lines = wrap(f"{ind}pure subroutine {sub}({', '.join(fn.inputs + fn.outputs)})", ind)
    lines.append(f"{ind}{INDENT}!$acc routine seq")
    lines.extend(wrap(f"{ind}{INDENT}real(dp), intent(in) :: {', '.join(fn.inputs)}", ind + INDENT))
    lines.extend(wrap(f"{ind}{INDENT}real(dp), intent(out) :: {', '.join(fn.outputs)}", ind + INDENT))
    reals, logicals, body = [], [], []
    seen = set()
    for lhs, rhs, is_new in fn.stmts:
        f_rhs = translate_expr(rhs)
        for nm in re.findall(r"[A-Za-z_]\w*", f_rhs):
            if nm in UTIL_NAMES or nm in CONSTANT_NAMES:
                used.add(nm)
        if is_new:
            if lhs in seen:
                raise Unsupported(f"{lhs} declared twice in {fn.name}")
            seen.add(lhs)
            (logicals if is_logical(f_rhs) else reals).append(lhs)
        elif lhs not in fn.outputs:
            raise Unsupported(f"assignment to undeclared {lhs} in {fn.name}")
        body.extend(wrap(f"{ind}{INDENT}{lhs} = {f_rhs}", ind + INDENT))
    for group, typ in ((logicals, "logical"), (reals, "real(dp)")):
        for i in range(0, len(group), 8):
            lines.append(f"{ind}{INDENT}{typ} :: {', '.join(group[i:i + 8])}")
    lines.extend(body)
    lines.append(f"{ind}end subroutine {sub}")
    lines.append("")
    return lines


def emit_wrapper(k: Kernel, fn: Function) -> list[str]:
    """The screening_interface.hpp rule, one template for LDA and GGA."""
    pub = f"{k.name}_{fn.name[5:-5]}"  # drop "eval_" and "_impl"
    impl = f"{k.name}_{fn.name[5:]}"
    polar = fn.name.endswith("polar_impl") and "unpolar" not in fn.name
    ind = INDENT
    lines = wrap(f"{ind}pure subroutine {pub}({', '.join(fn.inputs + fn.outputs)})", ind)
    lines.append(f"{ind}{INDENT}!$acc routine seq")
    lines.extend(wrap(f"{ind}{INDENT}real(dp), intent(in) :: {', '.join(fn.inputs)}", ind + INDENT))
    lines.extend(wrap(f"{ind}{INDENT}real(dp), intent(out) :: {', '.join(fn.outputs)}", ind + INDENT))
    clamped = [x + "_c" for x in fn.inputs]
    lines.append(f"{ind}{INDENT}real(dp) :: {', '.join(clamped)}")
    for o in fn.outputs:
        lines.append(f"{ind}{INDENT}{o} = 0.0_dp")
    total = "rho_a + rho_b" if polar else "rho"
    lines.append(f"{ind}{INDENT}if ({total} > dens_tol) then")
    b = ind + INDENT + INDENT
    if polar:
        lines.append(f"{b}rho_a_c = max(rho_a, dens_tol)")
        lines.append(f"{b}rho_b_c = max(rho_b, dens_tol)")
        if k.family == "gga":
            lines.append(f"{b}sigma_aa_c = max(sigma_aa, sigma_tol*sigma_tol)")
            lines.append(f"{b}sigma_bb_c = max(sigma_bb, sigma_tol*sigma_tol)")
            lines.append(f"{b}sigma_ab_c = enforce_polar_sigma_constraints(sigma_aa_c, sigma_ab, sigma_bb_c)")
    else:
        lines.append(f"{b}rho_c = rho")
        if k.family == "gga":
            lines.append(f"{b}sigma_c = max(sigma, sigma_tol*sigma_tol)")
    call = f"{b}call {impl}({', '.join(clamped + fn.outputs)})"
    lines.extend(wrap(call, b))
    lines.append(f"{ind}{INDENT}end if")
    lines.append(f"{ind}end subroutine {pub}")
    lines.append("")
    return lines


def emit_module(k: Kernel, src: Path) -> str:
    used: set[str] = set()
    impls, wrappers, publics = [], [], []
    for fn in k.functions:
        impls.extend(emit_impl(k, fn, used))
        wrappers.extend(emit_wrapper(k, fn))
        publics.append(f"{k.name}_{fn.name[5:-5]}")
    if k.family == "gga":
        used.add("enforce_polar_sigma_constraints")
    mod = f"trc_xc_{k.name}"
    head = [
        "!",
        f"! {k.name.upper()}: {k.family.upper()} exchange-correlation kernel.",
        "!",
        f"! GENERATED by tools/exchcxx_to_fortran.py from ExchCXX's",
        f"! {src.relative_to(src.parents[5])} -- do not edit.",
        "!",
        "! The arithmetic is ExchCXX's (MIT; Copyright (c) 2020-2024, The Regents",
        "! of the University of California, through Lawrence Berkeley National",
        "! Laboratory; portions Copyright (c) Microsoft Corporation), itself",
        "! generated by Maple from the functional's definition. See",
        "! https://github.com/wavefunction91/ExchCXX.",
        "!",
        "! Every routine is scalar, pure, and callable from `do concurrent`.",
        "! The public `<name>_<what>_<polar>` routines apply the density and",
        "! gradient thresholds the way ExchCXX's screening_interface.hpp does;",
        "! the `_impl` bodies below them are the bare Maple expressions.",
        "!",
        "! Unpolarised routines take the TOTAL density and its gradient",
        "! invariant; polarised ones take the two spin densities and the three",
        "! gradient invariants. `eps` is the energy per particle, so the",
        "! integrand is rho * eps.",
        "!",
        f"module {mod}",
        f"{INDENT}use trc_boys, only: dp",
    ]
    imports = sorted(used)
    if imports:
        head.extend(wrap(f"{INDENT}use trc_xc_util, only: {', '.join(imports)}", INDENT))
    head.append(f"{INDENT}implicit none")
    head.append(f"{INDENT}private")
    head.append("")
    for p in publics:
        head.append(f"{INDENT}public :: {p}")
    head.append("")
    for name, val in k.params:
        head.append(f"{INDENT}real(dp), parameter :: {name} = {val}")
    head.append("")
    head.append("contains")
    head.append("")
    body = wrappers + impls
    tail = [f"end module {mod}", ""]
    return "\n".join(head + body + tail)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exchcxx", required=True, type=Path, help="ExchCXX checkout")
    ap.add_argument("-o", "--outdir", required=True, type=Path)
    ap.add_argument("kernels", nargs="+", help="kernel header stems, e.g. pbe_x")
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    rc = 0
    for stem in args.kernels:
        src = args.exchcxx / KERNEL_DIR / f"{stem}.hpp"
        if not src.exists():
            print(f"{stem}: no such kernel at {src}", file=sys.stderr)
            rc = 1
            continue
        try:
            k = parse_kernel(src)
            out = args.outdir / f"trc_xc_{stem}.F90"
            out.write_text(emit_module(k, src.resolve()))
            print(f"{stem}: {k.family}, {len(k.params)} parameters -> {out}")
        except Unsupported as e:
            print(f"{stem}: unsupported: {e}", file=sys.stderr)
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
