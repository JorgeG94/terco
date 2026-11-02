#!/usr/bin/env python3
"""Map terco's modules, their contents and who calls what.

    python3 tools/modmap.py graph  [--dir scf] [--no-pic] [--no-external] > g.mmd
    python3 tools/modmap.py layers
    python3 tools/modmap.py md     > MODULES.md
    python3 tools/modmap.py calls  <procedure|module> [--depth 2] [--reverse]
    python3 tools/modmap.py lint
    python3 tools/modmap.py json

Scans app/ and src/ (src/generated/ only with --generated). Regex based: it
relies on terco's regular style (one statement per line, `use m, only: ...`),
not on a full Fortran parser. The `calls` view matches `call name(` and
`name(` against known procedure names, so it is approximate; generic
interfaces and type-bound calls are matched by their bound/specific name.
Output goes to stdout; nothing is written to the repository.
"""
import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXTS = {".f90", ".F90", ".inc"}
INTRINSIC = {"iso_c_binding", "iso_fortran_env", "ieee_arithmetic",
             "ieee_exceptions", "ieee_features", "openacc", "omp_lib", "mpi",
             "mpi_f08", "omp_lib_kinds"}

RE_MOD = re.compile(r"^\s*(module|program)\s+(?!procedure\b)(\w+)\s*$", re.I)
RE_SUBMOD = re.compile(r"^\s*submodule\s*\(\s*(\w+)\s*\)\s*(\w+)", re.I)
RE_END_MOD = re.compile(r"^\s*end\s+(module|program|submodule)\b", re.I)
RE_USE = re.compile(
    r"^\s*use\s*(?:,\s*(?:non_)?intrinsic\s*)?(?:::)?\s*(\w+)\s*(?:,\s*only\s*:\s*(.*))?$",
    re.I)
RE_PROC = re.compile(
    r"^\s*((?:(?:pure|elemental|recursive|impure|module)\s+)*)"
    r"(?:(?:integer|real|logical|complex|character|type)\b[^!]*?\s+)?"
    r"(subroutine|function)\s+(\w+)", re.I)
RE_END_PROC = re.compile(r"^\s*end\s+(subroutine|function)\b", re.I)
RE_TYPE = re.compile(r"^\s*type\s*(?:,[^:]*)?::\s*(\w+)", re.I)
RE_TYPE_ALT = re.compile(r"^\s*type\s+(?!is\b|\()(\w+)\s*$", re.I)
RE_END_TYPE = re.compile(r"^\s*end\s+type\b", re.I)
RE_IFACE = re.compile(r"^\s*(abstract\s+)?interface\s+(\w+)\s*$", re.I)
RE_BOUND = re.compile(
    r"^\s*(?:procedure|generic)\s*(?:,[^:]*)?::\s*(\w+)\s*(?:=>\s*(\w+))?", re.I)
RE_PUBLIC = re.compile(r"^\s*public\s*(?:::)?\s*(.*)$", re.I)
RE_PRIVATE_ALL = re.compile(r"^\s*private\s*$", re.I)
RE_CALL = re.compile(r"\bcall\s+(?:\w+(?:\([^)]*\))?%)*(\w+)", re.I)
RE_WORD_CALL = re.compile(r"\b(\w+)\s*\(")
RE_ACC_SEQ = re.compile(r"^\s*!\$acc\s+routine\s+seq", re.I)
RE_FORD = re.compile(r"^\s*!!\s?(.*)$")


def strip(line):
    """Drop a trailing comment (outside strings) but keep `!$` directives."""
    if line.lstrip().startswith("!"):
        return line
    out, q = [], None
    for ch in line:
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "!":
            break
        out.append(ch)
    return "".join(out)


def logical_lines(text):
    """Yield (lineno, statement) joining `&` continuations."""
    buf, start = "", 0
    for i, raw in enumerate(text.splitlines(), 1):
        s = strip(raw).rstrip()
        if not buf:
            start = i
        cont = s.endswith("&")
        s = s[:-1] if cont else s
        if buf and s.lstrip().startswith("&"):
            s = s.lstrip()[1:]
        buf += s
        if not cont:
            for part in (buf.split(";") if not buf.lstrip().startswith("!") else [buf]):
                yield start, part
            buf = ""


class Proc:
    def __init__(self, name, kind, attrs, line):
        self.name, self.kind, self.line = name, kind, line
        self.attrs = sorted(a for a in attrs.lower().split() if a)
        self.device = False
        self.calls = set()
        self.owner = None
        self.public = True
        self.doc = ""


class Unit:
    def __init__(self, name, kind, path, line):
        self.name, self.kind, self.path, self.line = name, kind, path, line
        self.uses = {}            # module -> set(only names) or None
        self.procs = {}
        self.types = {}           # name -> [bound names]
        self.ifaces = set()
        self.public_names = None  # explicit public list
        self.default_private = False
        self.doc = ""
        self.words = set()

    @property
    def layer(self):
        p = self.path.parts
        if p[0] == "src":
            return p[1] if len(p) > 2 else "src"
        return p[0]


def parse_file(path):
    rel = path.relative_to(ROOT)
    text = path.read_text(errors="replace")
    raw_lines = text.splitlines()
    units, cur, cur_proc, cur_type = [], None, None, None
    depth_iface = 0
    pending_acc = False
    for ln, st in logical_lines(text):
        if st.lstrip().startswith("!"):
            if RE_ACC_SEQ.match(st):
                pending_acc = True
            continue
        if not st.strip():
            continue
        m = RE_SUBMOD.match(st)
        if m:
            cur = Unit(m.group(2), "submodule", rel, ln)
            cur.uses[m.group(1)] = None
            units.append(cur)
            continue
        m = RE_MOD.match(st)
        if m and cur is None:
            cur = Unit(m.group(2), m.group(1).lower(), rel, ln)
            units.append(cur)
            # FORD summary: first `!!` after the header
            for r in raw_lines[ln:ln + 6]:
                d = RE_FORD.match(r)
                if d:
                    cur.doc = d.group(1).strip()
                    break
                if r.strip() and not r.lstrip().startswith("!"):
                    break
            continue
        if cur is None:
            continue
        if RE_END_MOD.match(st):
            cur = None
            continue
        words = re.findall(r"\w+", st.lower())
        cur.words.update(words)
        m = RE_USE.match(st)
        if m:
            mod = m.group(1).lower()
            only = m.group(2)
            if only is None:
                cur.uses.setdefault(mod, None)
            else:
                names = set()
                for n in only.split(","):
                    n = n.strip().lower()
                    if not n:
                        continue
                    if "=>" in n:  # rename: keep the original
                        n = n.split("=>")[1].strip()
                    names.add(n)
                if cur.uses.get(mod) is not None or mod not in cur.uses:
                    cur.uses.setdefault(mod, set()).update(names)
            continue
        if RE_PRIVATE_ALL.match(st) and cur_proc is None and cur_type is None:
            cur.default_private = True
            continue
        m = RE_PUBLIC.match(st)
        if m and cur_proc is None and cur_type is None:
            names = [n.strip().lower() for n in m.group(1).split(",") if n.strip()]
            cur.public_names = (cur.public_names or set()) | set(names)
            continue
        m = RE_IFACE.match(st)
        if m:
            cur.ifaces.add(m.group(2).lower())
            depth_iface += 1
            continue
        if re.match(r"^\s*end\s+interface\b", st, re.I):
            depth_iface = max(0, depth_iface - 1)
            continue
        if depth_iface:
            continue  # interface bodies are not definitions
        m = RE_TYPE.match(st) or RE_TYPE_ALT.match(st)
        if m and not re.match(r"^\s*type\s*\(", st, re.I):
            cur_type = m.group(1).lower()
            cur.types[cur_type] = []
            continue
        if cur_type and RE_END_TYPE.match(st):
            cur_type = None
            continue
        if cur_type:
            b = RE_BOUND.match(st)
            if b:
                cur.types[cur_type].append((b.group(2) or b.group(1)).lower())
            continue
        m = RE_PROC.match(st)
        if m and not re.match(r"^\s*end\b", st, re.I) and cur_proc is None:
            if re.match(r"^\s*module\s+procedure\b", st, re.I):
                continue
            p = Proc(m.group(3).lower(), m.group(2).lower(), m.group(1), ln)
            p.owner = cur.name
            p.device = pending_acc
            pending_acc = False
            for r in raw_lines[ln:ln + 4]:
                d = RE_FORD.match(r)
                if d:
                    p.doc = d.group(1).strip()
                    break
                if r.strip() and not r.lstrip().startswith("!"):
                    break
            cur.procs[p.name] = p
            cur_proc = p
            continue
        if cur_proc is not None:
            if RE_END_PROC.match(st):
                cur_proc = None
                continue
            if RE_ACC_SEQ.match(st):
                cur_proc.device = True
            for c in RE_CALL.findall(st):
                cur_proc.calls.add(c.lower())
            for c in RE_WORD_CALL.findall(st):
                cur_proc.calls.add(c.lower())
    # visibility
    for u in units:
        for p in u.procs.values():
            if u.public_names is not None:
                p.public = p.name in u.public_names or u.public_names == set() and not u.default_private
            if u.default_private:
                p.public = p.name in (u.public_names or set())
    return units


def collect(generated=False):
    files = []
    for base in ("app", "src"):
        for p in sorted((ROOT / base).rglob("*")):
            if p.suffix not in EXTS or not p.is_file():
                continue
            parts = p.relative_to(ROOT).parts
            if "vendor" in parts or (not generated and "generated" in parts):
                continue
            files.append(p)
    units = []
    for f in files:
        units.extend(parse_file(f))
    return {u.name.lower(): u for u in units}


def classify(mod, units):
    if mod in units:
        return "terco"
    if mod.startswith("pic_"):
        return "pic"
    if mod in INTRINSIC or mod.startswith(("iso_", "ieee_")):
        return "intrinsic"
    return "external"


def sanitize(n):
    return re.sub(r"\W", "_", n)


# ---------------------------------------------------------------- views
def cmd_graph(units, a):
    keep = {n: u for n, u in units.items() if not a.dir or u.layer in a.dir}
    out = ["flowchart LR"]
    by_layer = defaultdict(list)
    for n, u in keep.items():
        by_layer[u.layer].append(n)
    for layer, names in sorted(by_layer.items()):
        out.append(f"  subgraph {sanitize(layer)}[{layer}]")
        for n in sorted(names):
            out.append(f"    {sanitize(n)}[{n}]")
        out.append("  end")
    extras, edges = set(), set()
    for n, u in keep.items():
        for m in u.uses:
            k = classify(m, units)
            if k == "terco" and m in keep:
                edges.add((n, m))
            elif k == "terco" and a.dir:
                extras.add(m)
                edges.add((n, m))
            elif k == "pic" and not a.no_pic:
                extras.add(m)
                edges.add((n, m))
            elif k == "external" and not a.no_external:
                extras.add(m)
                edges.add((n, m))
    for m in sorted(extras):
        out.append(f"  {sanitize(m)}([{m}])")
    for s, d in sorted(edges):
        out.append(f"  {sanitize(s)} --> {sanitize(d)}")
    print("\n".join(out))


def layer_edges(units):
    e = defaultdict(set)
    for n, u in units.items():
        for m in u.uses:
            if m in units and units[m].layer != u.layer:
                e[(u.layer, units[m].layer)].add((n, m))
    return e


# dependency rank: lower may not use higher
LAYER_ORDER = ["core", "basis", "grid", "eri", "xc", "solvation", "scf",
               "post_hf", "api", "src", "app"]


def cmd_layers(units, a):
    e = layer_edges(units)
    layers = sorted({u.layer for u in units.values()},
                    key=lambda x: LAYER_ORDER.index(x) if x in LAYER_ORDER else 99)
    w = max(len(x) for x in layers) + 1
    print("rows USE columns (count of module pairs)\n")
    print(" " * w + " ".join(f"{c[:7]:>7}" for c in layers))
    for r in layers:
        print(f"{r:<{w}}" + " ".join(
            f"{len(e.get((r, c), ())) or '.':>7}" for c in layers))
    print("\nupward dependencies (lower layer uses a higher one; order: "
          + " < ".join(LAYER_ORDER) + "):")
    found = False
    for (s, d), pairs in sorted(e.items()):
        if s in LAYER_ORDER and d in LAYER_ORDER and \
                LAYER_ORDER.index(s) < LAYER_ORDER.index(d):
            found = True
            for x, y in sorted(pairs):
                print(f"  {s}/{x} -> {d}/{y}")
    if not found:
        print("  none")


def vis(u, p):
    return "public" if p.public and not u.default_private or \
        (u.default_private and p.name in (u.public_names or set())) else "private"


def cmd_md(units, a):
    print("# terco module reference (generated by tools/modmap.py)\n")
    cur = None
    for n, u in sorted(units.items(), key=lambda kv: (kv[1].layer, kv[0])):
        if u.layer != cur:
            cur = u.layer
            print(f"\n## {cur}\n")
        print(f"### `{u.name}` ({u.kind})\n")
        print(f"`{u.path}:{u.line}`" + (f" - {u.doc}" if u.doc else "") + "\n")
        t = sorted(m for m in u.uses if classify(m, units) == "terco")
        o = sorted(m for m in u.uses if classify(m, units) != "terco")
        if t:
            print("- uses (terco): " + ", ".join(f"`{x}`" for x in t))
        if o:
            print("- uses (other): " + ", ".join(f"`{x}`" for x in o))
        if u.types:
            print("- types: " + ", ".join(
                f"`{k}`" + (f" ({', '.join(v)})" if v else "")
                for k, v in sorted(u.types.items())))
        if u.ifaces:
            print("- generics: " + ", ".join(f"`{x}`" for x in sorted(u.ifaces)))
        pub = [p for p in u.procs.values() if vis(u, p) == "public"]
        if pub:
            print("\n| procedure | kind | attrs | line | summary |\n|---|---|---|---|---|")
            for p in sorted(pub, key=lambda x: x.line):
                at = " ".join(p.attrs + (["acc-seq"] if p.device else []))
                print(f"| `{p.name}` | {p.kind} | {at} | {p.line} | {p.doc} |")
        npriv = len(u.procs) - len(pub)
        if npriv:
            print(f"\n_{npriv} private procedure(s)_")
        print()


def build_proc_index(units):
    idx = defaultdict(list)
    for u in units.values():
        for p in u.procs.values():
            idx[p.name].append(p)
        for t, bound in u.types.items():
            for b in bound:
                if b in u.procs:
                    idx[b].append(u.procs[b])
    return idx


def resolve(p, units, idx):
    """Known procedures p may call, preferring ones visible through its uses."""
    u = units[p.owner]
    vis_mods = {u.name}
    vis_mods.update(m for m in u.uses if m in units)
    out = set()
    for c in p.calls:
        if c == p.name:
            continue
        for q in idx.get(c, ()):
            if q.owner in vis_mods:
                out.add(q)
    return out


def cmd_calls(units, a):
    idx = build_proc_index(units)
    target = a.target.lower()
    roots = list(idx.get(target, []))
    if not roots and target in units:
        roots = [p for p in units[target].procs.values() if vis(units[target], p) == "public"]
    if not roots:
        sys.exit(f"no procedure or module named {a.target!r}")
    rev = defaultdict(set)
    fwd = {}
    for u in units.values():
        for p in u.procs.values():
            fwd[p] = resolve(p, units, idx)
            for q in fwd[p]:
                rev[q].add(p)
    graph = rev if a.reverse else fwd
    seen = set()
    lines = ["flowchart LR"] if a.mermaid else []

    def label(p):
        return f"{p.owner}::{p.name}"

    def walk(p, depth, indent):
        if a.mermaid:
            for q in sorted(graph.get(p, ()), key=label):
                lines.append(f'  {sanitize(label(p))}["{label(p)}"] --> '
                             f'{sanitize(label(q))}["{label(q)}"]'
                             if not a.reverse else
                             f'  {sanitize(label(q))}["{label(q)}"] --> '
                             f'{sanitize(label(p))}["{label(p)}"]')
                if q not in seen and depth < a.depth:
                    seen.add(q)
                    walk(q, depth + 1, indent)
            return
        for q in sorted(graph.get(p, ()), key=label):
            again = q in seen
            lines.append("  " * indent + label(q) + (" (...)" if again and graph.get(q) else ""))
            if not again and depth < a.depth:
                seen.add(q)
                walk(q, depth + 1, indent + 1)

    for r in roots:
        seen.add(r)
        if not a.mermaid:
            lines.append(label(r) + f"   [{units[r.owner].path}:{r.line}]")
        walk(r, 1, 1)
    print("\n".join(lines))


# ---------------------------------------------------------------- lint
def find_cycles(units):
    g = {n: [m for m in u.uses if m in units and m != n] for n, u in units.items()}
    idx, low, onst, st, sccs, c = {}, {}, set(), [], [], [0]

    def sc(v):
        idx[v] = low[v] = c[0]; c[0] += 1
        st.append(v); onst.add(v)
        for w in g[v]:
            if w not in idx:
                sc(w); low[v] = min(low[v], low[w])
            elif w in onst:
                low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp = []
            while True:
                w = st.pop(); onst.discard(w); comp.append(w)
                if w == v:
                    break
            if len(comp) > 1:
                sccs.append(sorted(comp))
    sys.setrecursionlimit(10000)
    for v in g:
        if v not in idx:
            sc(v)
    return sccs


def tokens(name):
    return [t for t in name.lower().split("_") if t]


def cmd_lint(units, a):
    issues = defaultdict(list)
    idx = build_proc_index(units)
    for n, u in sorted(units.items()):
        if u.kind == "program":
            continue
        if not n.startswith("trc_"):
            issues["module name lacks trc_ prefix"].append(f"{n}  ({u.path})")
        stem = u.path.stem.lower()
        if stem != n and not stem.startswith(n):
            issues["module name differs from file name"].append(f"{n}  ({u.path})")
        pubs = [p for p in u.procs.values() if vis(u, p) == "public"]
        if len(pubs) > a.max_public:
            issues[f"public surface > {a.max_public} procedures"].append(
                f"{n}: {len(pubs)}")
        short = n[4:] if n.startswith("trc_") else n
        for p in pubs:
            if not p.name.startswith("trc_") and short not in p.name and \
                    p.name not in {b for bl in u.types.values() for b in bl}:
                issues["public procedure carries neither trc_ nor module name"].append(
                    f"{n}::{p.name}  ({u.path}:{p.line})")
        # unused imports (explicit only-lists)
        for m, names in u.uses.items():
            if names:
                for nm in sorted(names):
                    if sum(1 for w in [nm] if w in u.words) and \
                            not any(nm in p_calls for p_calls in
                                    [pp.calls for pp in u.procs.values()]) \
                            and not _used(u, nm):
                        issues["imported name never referenced"].append(
                            f"{n}: {m}: only {nm}")
    for cyc in find_cycles(units):
        issues["module dependency cycle"].append(" <-> ".join(cyc))
    # duplicate procedure names across modules
    for name, ps in sorted(idx.items()):
        owners = sorted({p.owner for p in ps if vis(units[p.owner], p) == "public"})
        if len(owners) > 1 and not name.startswith(("init", "release")):
            issues["same public procedure name in several modules"].append(
                f"{name}: {', '.join(owners)}")
    # near-duplicate names: same token set, different order or synonyms
    syn = {"calc": "compute", "calculate": "compute", "comp": "compute",
           "build": "make", "form": "make", "get": "fetch", "free": "release",
           "destroy": "release", "finalize": "release", "delete": "release",
           "alloc": "init", "allocate": "init", "setup": "init",
           "initialize": "init", "create": "init", "new": "init"}
    groups = defaultdict(set)
    for name, ps in idx.items():
        if not any(vis(units[p.owner], p) == "public" for p in ps):
            continue
        key = tuple(sorted(syn.get(t, t) for t in tokens(name) if t != "trc"))
        groups[key].add(name)
    for key, names in sorted(groups.items()):
        if len(names) > 1 and key:
            issues["near-duplicate names (reordered/synonym verbs)"].append(
                ", ".join(sorted(names)))
    verbs = defaultdict(set)
    for name in idx:
        t = tokens(name)
        if len(t) > 1 and t[0] != "trc":
            verbs[t[0]].add(name)
    for v in ("calc", "calculate", "compute", "build", "make", "form", "free",
              "destroy", "release", "finalize", "setup", "init", "initialize"):
        if len(verbs.get(v, ())) >= 3:
            issues["verb usage tally (informational)"].append(
                f"{v}_*: {len(verbs[v])}")
    # dead public API: never referenced in any other unit
    all_words = {n: u.words for n, u in units.items()}
    for n, u in sorted(units.items()):
        for p in u.procs.values():
            if vis(u, p) != "public" or u.kind == "program":
                continue
            if not any(p.name in w for m, w in all_words.items() if m != n):
                issues["public procedure not referenced outside its module (src/app only)"].append(
                    f"{n}::{p.name}")
    if not issues:
        print("no findings")
        return
    for k, v in issues.items():
        print(f"\n## {k} ({len(v)})")
        for x in v[:a.limit]:
            print("  " + x)
        if len(v) > a.limit:
            print(f"  ... {len(v) - a.limit} more (use --limit)")


def _used(u, name):
    """Name appears on a line other than its own use statement."""
    return u._count.get(name, 0) > 1 if hasattr(u, "_count") else True


def cmd_json(units, a):
    d = {}
    for n, u in units.items():
        d[n] = {"path": str(u.path), "kind": u.kind, "layer": u.layer,
                "uses": {m: (sorted(v) if v else None) for m, v in u.uses.items()},
                "types": u.types, "generics": sorted(u.ifaces),
                "procedures": {p.name: {"kind": p.kind, "line": p.line,
                                        "attrs": p.attrs, "acc_seq": p.device,
                                        "public": vis(u, p) == "public",
                                        "calls": sorted(p.calls & set(build_proc_index(units)))}
                               for p in u.procs.values()}}
    print(json.dumps(d, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--generated", action="store_true",
                    help="also scan src/generated/")
    sp = ap.add_subparsers(dest="cmd", required=True)
    g = sp.add_parser("graph")
    g.add_argument("--dir", action="append", help="layer to focus on (repeatable)")
    g.add_argument("--no-pic", action="store_true")
    g.add_argument("--no-external", action="store_true")
    sp.add_parser("layers")
    sp.add_parser("md")
    c = sp.add_parser("calls")
    c.add_argument("target")
    c.add_argument("--depth", type=int, default=2)
    c.add_argument("--reverse", action="store_true", help="callers instead of callees")
    c.add_argument("--mermaid", action="store_true")
    l = sp.add_parser("lint")
    l.add_argument("--max-public", type=int, default=40)
    l.add_argument("--limit", type=int, default=40)
    sp.add_parser("json")
    a = ap.parse_args()
    units = collect(a.generated)
    # per-unit word counts, used by the unused-import check
    for u in units.values():
        u._count = defaultdict(int)
    for n, u in units.items():
        text = (ROOT / u.path).read_text(errors="replace").lower()
        for nm in {x for v in u.uses.values() if v for x in v}:
            u._count[nm] = len(re.findall(r"\b" + re.escape(nm) + r"\b", text))
    {"graph": cmd_graph, "layers": cmd_layers, "md": cmd_md, "calls": cmd_calls,
     "lint": cmd_lint, "json": cmd_json}[a.cmd](units, a)


if __name__ == "__main__":
    main()
