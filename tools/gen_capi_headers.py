#!/usr/bin/env python3
"""
Parse the `bind(c, name="...")` procedures and the public `integer(c_int),
parameter` declarations out of src/api/trc_capi.F90 and emit:

  include/terco.h              C prototypes, `extern "C"` guards, the
                                status/guess/AO-norm/ABI-version constants,
                                an opaque handle typedef.
  include/trc_c_interfaces.f90  the Fortran interface module host codes
                                compile instead of restating the ABI by hand.

WHY GENERATE RATHER THAN HAND-MAINTAIN
---------------------------------------
Both files are a second copy of the same contract trc_capi.F90 already
states once, in iso_c_binding. Two of them, kept in step by hand, has
nothing to catch drift: add an argument in trc_capi.F90, forget one of the
two copies, and the forgotten one still compiles AND still links -- a C
symbol carries no signature -- and fails only at run time, by reading
whatever was in the register. Deriving both from trc_capi.F90 removes the
hand-copying step the drift needs.

WHAT IS PARSED
---------------
Only two shapes, both at module scope in trc_capi.F90:

  * a FUNCTION whose header carries `bind(c, name="...")`. Its leading
    `!!` doc-comment lines (the contiguous run immediately after the
    header, before the first declaration) are carried into both outputs
    as that procedure's documentation. Every other procedure in the file
    -- anything without a `bind(c, name=...)` clause -- is implementation
    and is not part of the ABI this tool restates.

  * a module-level `integer(c_int), parameter, public :: NAME = EXPR[,
    ...]` declaration whose NAME starts with `TRC_`. `EXPR` is either an
    integer literal or a bare identifier that is itself a renamed `use`
    of another module's own named parameter (TRC_AO_NORM_SHELL/UNIT's
    `ao_norm_shell_ => TRC_AO_NORM_SHELL` -- see `resolve_param`); it is
    followed there and resolved to a literal for both outputs, which have
    no `use` of their own to inherit it from.

SUPPORTED DUMMY-ARGUMENT SHAPES
---------------------------------
  integer(c_int)          value, or intent(in)/(out)/(inout) -- scalar, an
  real(c_double)           assumed-size array NAME(*), or an explicit-shape
                            array bound by another dummy, NAME(n) -- both
                            compile to a flat pointer at the C boundary, with
                            no array descriptor, which is the property that
                            matters here; the Fortran interface file copies
                            whichever shape trc_capi.F90 itself declared.
  type(c_ptr)              value (a handle in) or intent(out) (a handle out)
  character(kind=c_char)   intent(in)/(out), an assumed-size array, NAME(*)

Anything else -- `optional`, an assumed-shape array (NAME(:)), `allocatable`,
`pointer`, a derived type crossing the C boundary, a kind this list does not
name -- is refused with a hard error rather than guessed at: a header
generator that invents a plausible C type for a shape it has not been told
about would ship a silently wrong ABI. (Most of these the Fortran standard
already refuses inside a `bind(c)` procedure, so this tool's own check is
the backstop for the rest: an unfamiliar KIND, or a shape it cannot read.)

USAGE
    tools/gen_capi_headers.py
        regenerate both files in place, from src/api/trc_capi.F90.
    tools/gen_capi_headers.py --check
        regenerate into a scratch location and diff against what is
        committed, without touching the working tree (the CI step).
"""

import argparse
import difflib
import os
import re
import sys
from collections import OrderedDict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAPI_SOURCE = os.path.join(REPO_ROOT, "src", "api", "trc_capi.F90")
HEADER_OUT = os.path.join(REPO_ROOT, "include", "terco.h")
INTERFACE_OUT = os.path.join(REPO_ROOT, "include", "trc_c_interfaces.f90")


class GeneratorError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Free-form source -> logical units (comments and continuation-joined
# statements), so a header spanning two or three physical lines reads as one.
# ---------------------------------------------------------------------------

def strip_inline_comment(line):
    """Cut a code line at the first `!` that is not inside a quoted string."""
    quote = None
    for idx, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            return line[:idx]
    return line


def read_logical_units(path):
    with open(path, "r") as f:
        lines = f.readlines()
    units = []
    i = 0
    n = len(lines)
    while i < n:
        raw = lines[i].rstrip("\n")
        stripped = raw.strip()
        if stripped == "":
            i += 1
            continue
        if stripped.startswith("!"):
            units.append(("comment", i + 1, stripped))
            i += 1
            continue
        parts = []
        first = True
        start_lineno = i + 1
        while True:
            text = strip_inline_comment(lines[i].rstrip("\n")).strip()
            if not first and text.startswith("&"):
                text = text[1:].lstrip()
            cont = text.endswith("&")
            if cont:
                text = text[:-1].rstrip()
            parts.append(text)
            i += 1
            first = False
            if not cont:
                break
        units.append(("code", start_lineno, " ".join(p for p in parts if p)))
    return units


def split_toplevel_commas(text):
    out = []
    depth = 0
    buf = []
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    out.append("".join(buf))
    return [s.strip() for s in out if s.strip()]


DECL_RE = re.compile(
    r"^(integer|real|character|logical|type|class)\s*(\([^)]*\))?\s*"
    r"(?:,\s*(.+?))?\s*::\s*(.+)$",
    re.IGNORECASE,
)
FUNC_RE = re.compile(
    r"^(?:pure\s+|elemental\s+|recursive\s+)*function\s+(\w+)\s*\((.*?)\)\s*"
    r"(?:result\s*\(\s*\w+\s*\))?\s*bind\s*\(\s*c\s*,\s*name\s*=\s*\"([^\"]+)\"\s*\)\s*$",
    re.IGNORECASE,
)
END_FUNC_RE = re.compile(r"^end\s+function\b", re.IGNORECASE)
PARAM_RE = re.compile(
    r"^integer\s*\(\s*c_int\s*\)\s*,\s*parameter\s*,\s*public\s*::\s*(.+)$",
    re.IGNORECASE,
)
USE_RENAME_RE = re.compile(
    r"^use\s+(\w+)\s*,\s*only\s*:\s*(.+)$", re.IGNORECASE
)


class Arg:
    def __init__(self, name, typename, kindspec, attrs, entity_text):
        self.name = name
        self.typename = typename.lower()
        self.kindspec = kindspec
        self.attrs = attrs.lower() if attrs else ""
        self.entity_text = entity_text

    @property
    def kind(self):
        if not self.kindspec:
            return None
        inner = self.kindspec.strip()[1:-1]
        if "=" in inner:
            inner = inner.split("=", 1)[1]
        return inner.strip().lower()

    @property
    def is_value(self):
        return "value" in self.attrs

    @property
    def is_array(self):
        return "(" in self.entity_text


class Proc:
    def __init__(self, fname, cname, args, doc):
        self.fname = fname
        self.cname = cname
        self.args = args
        self.doc = doc


def parse_capi(path):
    units = read_logical_units(path)
    procs = []
    raw_params = OrderedDict()  # name -> raw RHS expression text
    rename_map = {}  # local name -> (module, remote name)

    i = 0
    n = len(units)
    while i < n:
        kind, lineno, text = units[i]
        if kind != "code":
            i += 1
            continue

        m_use = USE_RENAME_RE.match(text)
        if m_use:
            module_name = m_use.group(1)
            for item in split_toplevel_commas(m_use.group(2)):
                if "=>" in item:
                    local, remote = [s.strip() for s in item.split("=>", 1)]
                    rename_map[local] = (module_name, remote)
            i += 1
            continue

        m_param = PARAM_RE.match(text)
        if m_param:
            for ent in split_toplevel_commas(m_param.group(1)):
                if "=" not in ent:
                    raise GeneratorError(
                        f"{path}:{lineno}: public TRC_* parameter without a value: {ent!r}"
                    )
                name, expr = [s.strip() for s in ent.split("=", 1)]
                if name.startswith("TRC_"):
                    raw_params[name] = expr
            i += 1
            continue

        m_func = FUNC_RE.match(text)
        if m_func:
            fname, argtext, cname = m_func.groups()
            argnames = [a.strip() for a in argtext.split(",") if a.strip()]

            j = i + 1
            doc = []
            while j < n and units[j][0] == "comment" and units[j][2].lstrip().startswith("!!"):
                line = units[j][2].lstrip()
                line = line[2:]
                if line.startswith(" "):
                    line = line[1:]
                doc.append(line)
                j += 1

            found = {}
            k = j
            while k < n:
                ukind, ulineno, utext = units[k]
                if ukind == "comment":
                    k += 1
                    continue
                m_decl = DECL_RE.match(utext)
                if not m_decl:
                    break
                typename, kindspec, attrs, entities = m_decl.groups()
                attrs = attrs or ""
                if "value" in attrs.lower() or "intent" in attrs.lower():
                    for ent in split_toplevel_commas(entities):
                        base = re.match(r"^(\w+)", ent)
                        base = base.group(1) if base else ent
                        if base in argnames and base not in found:
                            found[base] = Arg(base, typename, kindspec, attrs, ent)
                k += 1

            missing = [a for a in argnames if a not in found]
            if missing:
                raise GeneratorError(
                    f"{path}:{lineno}: {fname} ({cname}): no dummy declaration "
                    f"found for argument(s) {missing} -- supported shapes only "
                    f"(see this tool's own module docstring)"
                )
            args = [found[a] for a in argnames]
            procs.append(Proc(fname, cname, args, doc))

            m_end = k
            while m_end < n:
                ukind, ulineno, utext = units[m_end]
                if ukind == "code" and END_FUNC_RE.match(utext):
                    break
                m_end += 1
            i = m_end + 1
            continue

        i += 1

    params = OrderedDict()
    for name, expr in raw_params.items():
        params[name] = resolve_param(name, expr, raw_params, rename_map)
    return procs, params


def resolve_param(name, expr, known, rename_map, _stack=None):
    expr = expr.strip()
    if re.fullmatch(r"-?\d+", expr):
        return int(expr)
    _stack = _stack or set()
    if expr in _stack:
        raise GeneratorError(f"circular constant reference resolving {name} = {expr}")
    if expr in known:
        return resolve_param(expr, str(known[expr]), known, rename_map, _stack | {expr})
    if expr in rename_map:
        module_name, remote_name = rename_map[expr]
        value = find_module_constant(module_name, remote_name)
        if value is not None:
            return value
    raise GeneratorError(
        f"cannot resolve public parameter {name} = {expr}: not an integer "
        f"literal, not another TRC_* parameter in this file, and not a "
        f"renamed `use` of a module this tool knows how to search"
    )


def find_module_constant(module_name, remote_name):
    """Search src/ for `module <module_name>` and, within it, a
    `parameter :: ... remote_name = <int> ...` declaration."""
    src_root = os.path.join(REPO_ROOT, "src")
    mod_re = re.compile(r"^\s*module\s+" + re.escape(module_name) + r"\b", re.IGNORECASE)
    for dirpath, _dirs, files in os.walk(src_root):
        for fn in files:
            if not fn.endswith(".F90") and not fn.endswith(".f90"):
                continue
            path = os.path.join(dirpath, fn)
            with open(path, "r") as f:
                text = f.read()
            if not any(mod_re.match(ln) for ln in text.splitlines()):
                continue
            for _kind, _lineno, line in read_logical_units(path):
                mm = re.match(
                    r"^(integer|real)\s*(\([^)]*\))?\s*,\s*parameter[^:]*::\s*(.+)$",
                    line, re.IGNORECASE,
                )
                if not mm:
                    continue
                for ent in split_toplevel_commas(mm.group(3)):
                    if "=" not in ent:
                        continue
                    ename, eexpr = [s.strip() for s in ent.split("=", 1)]
                    if ename == remote_name:
                        eexpr = eexpr.strip()
                        if re.fullmatch(r"-?\d+", eexpr):
                            return int(eexpr)
    return None


# ---------------------------------------------------------------------------
# C header
# ---------------------------------------------------------------------------

_INT_KINDS = {"c_int"}
_REAL_KINDS = {"c_double"}


def check_shape(arg, proc_name):
    if ":" in arg.entity_text:
        raise GeneratorError(
            f"{proc_name}: dummy {arg.name} is an assumed-shape array "
            f"({arg.entity_text!r}) -- not interoperable with C, and not a "
            f"shape this tool supports"
        )
    if "allocatable" in arg.attrs or "pointer" in arg.attrs:
        raise GeneratorError(f"{proc_name}: dummy {arg.name} is allocatable/pointer, not supported")
    if "optional" in arg.attrs:
        raise GeneratorError(f"{proc_name}: dummy {arg.name} is optional, not supported at the C boundary")


HANDLE_NAMES = {"handle", "basis", "pairs", "aux", "eri"}
"""by-value type(c_ptr) dummies that are terco handles; every other by-value
c_ptr is a nullable data pointer (const void* in terco.h)."""


def c_type_for(arg, proc_name):
    check_shape(arg, proc_name)
    kind = arg.kind
    if arg.typename == "type":
        if kind != "c_ptr":
            raise GeneratorError(f"{proc_name}: unsupported derived type dummy {arg.entity_text!r}")
        if arg.is_value:
            # A by-value c_ptr is either an opaque handle or a nullable data
            # pointer (an optional array: dguess, radii). Only the names
            # terco uses for its handles map to trc_handle_t; anything
            # else is data, so a C caller passes its double* uncast.
            if arg.name in HANDLE_NAMES:
                return "trc_handle_t"
            return "const void*"
        if "intent(out)" in arg.attrs:
            return "trc_handle_t*"
        raise GeneratorError(f"{proc_name}: type(c_ptr) dummy {arg.name} is neither value nor intent(out)")
    if arg.typename == "character":
        if kind != "c_char":
            raise GeneratorError(f"{proc_name}: unsupported character kind for {arg.name}")
        if "intent(in)" in arg.attrs:
            return "const char*"
        if "intent(out)" in arg.attrs or "intent(inout)" in arg.attrs:
            return "char*"
        raise GeneratorError(f"{proc_name}: character dummy {arg.name} has no intent")
    if arg.typename == "integer":
        if kind not in _INT_KINDS:
            raise GeneratorError(f"{proc_name}: unsupported integer kind for {arg.name}: {kind}")
        base = "int"
    elif arg.typename == "real":
        if kind not in _REAL_KINDS:
            raise GeneratorError(f"{proc_name}: unsupported real kind for {arg.name}: {kind}")
        base = "double"
    else:
        raise GeneratorError(f"{proc_name}: unsupported dummy type for {arg.name}: {arg.typename}")
    if arg.is_value:
        return base
    if "intent(in)" in arg.attrs:
        return f"const {base}*"
    if "intent(out)" in arg.attrs or "intent(inout)" in arg.attrs:
        return f"{base}*"
    raise GeneratorError(f"{proc_name}: dummy {arg.name} has neither `value` nor an `intent`")


def wrap_doc_c(doc):
    if not doc:
        return ""
    lines = ["/**"] + [" * " + d if d else " *" for d in doc] + [" */"]
    return "\n".join(lines) + "\n"


def render_header(procs, params):
    out = []
    out.append("/*")
    out.append(" * terco's C ABI -- GENERATED, do not edit by hand.")
    out.append(" *")
    out.append(" * tools/gen_capi_headers.py parses src/api/trc_capi.F90's own")
    out.append(" * `bind(c, name=...)` procedures and public TRC_* parameters and")
    out.append(" * writes this file and include/trc_c_interfaces.f90 from the same")
    out.append(" * pass, so the two cannot drift from each other or from the")
    out.append(" * library. Regenerate with that script; do not edit the output.")
    out.append(" */")
    out.append("#ifndef TERCO_H")
    out.append("#define TERCO_H")
    out.append("")
    out.append("#ifdef __cplusplus")
    out.append('extern "C" {')
    out.append("#endif")
    out.append("")
    out.append("/* An opaque handle to any object this API hands back -- a basis, a")
    out.append(" * shell-pair list, an ERI context, or a whole calculation context.")
    out.append(" * The Fortran side owns the memory; a caller only ever holds this")
    out.append(" * and passes it back. */")
    out.append("typedef void* trc_handle_t;")
    out.append("")
    out.append("/* Status codes. Zero is success; everything else is a reason. */")
    for name, value in params.items():
        out.append(f"#define {name} {value}")
    out.append("")
    for proc in procs:
        doc = wrap_doc_c(proc.doc)
        if doc:
            out.append(doc.rstrip())
        args_c = ", ".join(f"{c_type_for(a, proc.cname)} {a.name}" for a in proc.args) or "void"
        out.append(f"int {proc.cname}({args_c});")
        out.append("")
    out.append("#ifdef __cplusplus")
    out.append("}")
    out.append("#endif")
    out.append("")
    out.append("#endif /* TERCO_H */")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Fortran interface module
# ---------------------------------------------------------------------------

def f_decl_for(arg, proc_name):
    kind = arg.kind
    if arg.typename == "type":
        if kind != "c_ptr":
            raise GeneratorError(f"{proc_name}: unsupported derived type dummy {arg.entity_text!r}")
        attr = "value" if arg.is_value else "intent(out)"
        return f"type(c_ptr), {attr} :: {arg.name}", {"c_ptr"}
    if arg.typename == "character":
        if kind != "c_char":
            raise GeneratorError(f"{proc_name}: unsupported character kind for {arg.name}")
        if "intent(in)" in arg.attrs:
            attr = "intent(in)"
        elif "intent(out)" in arg.attrs:
            attr = "intent(out)"
        else:
            attr = "intent(inout)"
        return f"character(kind=c_char), {attr} :: {arg.entity_text}", {"c_char"}
    if arg.typename == "integer":
        if kind not in _INT_KINDS:
            raise GeneratorError(f"{proc_name}: unsupported integer kind for {arg.name}: {kind}")
        fkind = "c_int"
        ftype = "integer"
    elif arg.typename == "real":
        if kind not in _REAL_KINDS:
            raise GeneratorError(f"{proc_name}: unsupported real kind for {arg.name}: {kind}")
        fkind = "c_double"
        ftype = "real"
    else:
        raise GeneratorError(f"{proc_name}: unsupported dummy type for {arg.name}: {arg.typename}")
    if arg.is_value:
        attr = "value"
    elif "intent(in)" in arg.attrs:
        attr = "intent(in)"
    elif "intent(out)" in arg.attrs:
        attr = "intent(out)"
    elif "intent(inout)" in arg.attrs:
        attr = "intent(inout)"
    else:
        raise GeneratorError(f"{proc_name}: dummy {arg.name} has neither `value` nor an `intent`")
    decl = f"{ftype}({fkind}), {attr} :: {arg.entity_text}"
    return decl, {fkind}


def wrap_fortran_line(line, maxlen):
    """Split an over-long `... NAME(args...) ...` line across Fortran `&`
    continuations, breaking only between arguments (never inside a nested
    paren), re-indented to the opening paren's column."""
    if len(line) <= maxlen:
        return [line]
    # The argument list's own opening paren -- NOT the first `(` in the
    # line, which is `integer(c_int)`'s.
    fn = line.index("function")
    start = line.index("(", fn)
    depth = 0
    end = None
    for idx in range(start, len(line)):
        if line[idx] == "(":
            depth += 1
        elif line[idx] == ")":
            depth -= 1
            if depth == 0:
                end = idx
                break
    prefix = line[: start + 1]
    inner = line[start + 1 : end]
    suffix = line[end:]
    items = split_toplevel_commas(inner)
    indent = " " * len(prefix)
    lines = []
    cur = prefix
    for idx, item in enumerate(items):
        piece = item + (", " if idx < len(items) - 1 else "")
        if cur != prefix and len(cur) + len(piece) + 2 > maxlen:
            lines.append(cur.rstrip() + " &")
            cur = indent
        cur += piece
    if len(cur) + len(suffix) <= maxlen:
        lines.append(cur + suffix)
    else:
        lines.append(cur.rstrip() + " &")
        lines.append(indent + suffix.lstrip())
    return lines


def render_interface(procs, params):
    out = []
    out.append("!")
    out.append("! The authoritative Fortran declaration of terco's C ABI -- GENERATED,")
    out.append("! do not edit by hand.")
    out.append("!")
    out.append("! terco is nvfortran and its callers usually are not, so they reach it")
    out.append("! through `bind(c)`. tools/gen_capi_headers.py parses src/api/trc_capi.F90's")
    out.append("! own `bind(c, name=...)` procedures and public TRC_* parameters and")
    out.append("! writes this file and include/terco.h from the same pass, so a Fortran")
    out.append("! caller, a C caller and the library itself cannot drift from each other.")
    out.append("! Regenerate with that script; do not edit the output.")
    out.append("!")
    out.append("! `test/link_gfortran.F90` compiles this file with gfortran and asserts an")
    out.append("! SCF energy, so a drift anywhere in the chain is a wrong number in terco's")
    out.append("! own test suite before any host code sees it.")
    out.append("!")
    out.append("module trc_c_interfaces")
    out.append("   use, intrinsic :: iso_c_binding")
    out.append("   implicit none")
    out.append("   public")
    out.append("")
    for name, value in params.items():
        out.append(f"   integer(c_int), parameter :: {name} = {value}")
    out.append("")
    out.append("   interface")
    out.append("")
    for proc in procs:
        resname = "res_" if "res" in [a.name for a in proc.args] else "res"
        argnames = ", ".join(a.name for a in proc.args)
        header = f"      integer(c_int) function {proc.cname}({argnames}) bind(c) result({resname})"
        out.extend(wrap_fortran_line(header, 132))
        for d in proc.doc:
            out.append(f"         !! {d}" if d else "         !!")
        decls = []
        kinds = {"c_int"}
        for arg in proc.args:
            decl, used_kinds = f_decl_for(arg, proc.cname)
            decls.append(decl)
            kinds |= used_kinds
        import_list = ", ".join(sorted(kinds))
        out.append(f"         import :: {import_list}")
        out.append("         implicit none")
        for d in decls:
            out.append(f"         {d}")
        out.append(f"      end function {proc.cname}")
        out.append("")
    out.append("   end interface")
    out.append("")
    out.append("end module trc_c_interfaces")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                     help="diff against the committed files instead of writing them")
    args = ap.parse_args()

    procs, params = parse_capi(CAPI_SOURCE)
    header_text = render_header(procs, params)
    interface_text = render_interface(procs, params)

    if args.check:
        ok = True
        for path, text in ((HEADER_OUT, header_text), (INTERFACE_OUT, interface_text)):
            with open(path, "r") as f:
                current = f.read()
            if current != text:
                ok = False
                sys.stderr.write(f"{path} is out of date with {CAPI_SOURCE}:\n")
                diff = difflib.unified_diff(
                    current.splitlines(keepends=True), text.splitlines(keepends=True),
                    fromfile=path, tofile=path + " (regenerated)",
                )
                sys.stderr.writelines(diff)
        if not ok:
            sys.exit(1)
        return

    with open(HEADER_OUT, "w") as f:
        f.write(header_text)
    with open(INTERFACE_OUT, "w") as f:
        f.write(interface_text)
    print(f"wrote {HEADER_OUT}")
    print(f"wrote {INTERFACE_OUT}")


if __name__ == "__main__":
    main()
