"""Per-case staleness hash for the pyscf references. Standard library only.

generate.py (which runs pyscf, offline) and check_staleness.py (which runs
in CI, where pyscf is never installed) both import this module, so the hash
a reference is written with and the hash CI recomputes come from the same
code and cannot drift apart.

A reference's hash covers exactly what produced it:

  * the CASES entry in generate.py that wrote it (the dict value, so a
    lambda's arguments count), plus the transitive closure of every
    module-level function, class and constant of generate.py that entry
    refers to -- found by walking the parsed source with `ast`, not by
    reading names off the page;
  * the bytes of the input files on its "files:" line (geometry, basis
    JSON, the Fortran source that builds a fixed density or grid).

Each definition is hashed as a canonical dump of its parsed node (see
_canon) with docstrings removed, so comments, docstrings, blank lines and
line positions never restale a reference; a change to code that the case
actually runs always does. Editing one generator restales only the references it writes, and
editing a shared helper restales exactly the references whose generators
reach it.

Names generate.py imports (numpy, pyscf, cart_perm) are not part of the
closure: the pyscf/libxc versions are recorded on their own header lines.
"""
from __future__ import annotations

import ast
import hashlib
from pathlib import Path

#: Recorded on every reference written with this scheme; check_staleness.py
#: refuses a reference that does not carry it.
SCHEME = "case-ast-v1"

#: Never part of any case's closure: the registry itself and the driver.
_EXCLUDED = {"CASES", "main"}


def _strip_docstrings(node: ast.AST) -> ast.AST:
    """Drop the leading string-constant statement of every module, class
    and function body under `node` (in place) and return `node`."""
    for n in ast.walk(node):
        if isinstance(n, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = n.body
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                n.body = body[1:] or [ast.Pass()]
    return node


def _defined_names(stmt: ast.stmt) -> list[str]:
    """Module-level names a top-level statement binds."""
    if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
        return [stmt.name]
    targets = []
    if isinstance(stmt, ast.Assign):
        targets = stmt.targets
    elif isinstance(stmt, (ast.AnnAssign, ast.AugAssign)):
        targets = [stmt.target]
    names = []
    for t in targets:
        for n in ast.walk(t):
            if isinstance(n, ast.Name):
                names.append(n.id)
    return names


def _canon(node) -> str:
    """A serialisation of an AST that is the same on every Python 3.9+.

    `ast.dump` is not: 3.13 stopped printing empty optional fields, and
    3.12 added fields (type_params) older versions lack. Fields that are
    None or empty lists are skipped here, so all of those agree; position
    attributes (lineno, col_offset) are never fields, so moving code does
    not change it."""
    if isinstance(node, ast.AST):
        parts = []
        for name, value in ast.iter_fields(node):
            if value is None or (isinstance(value, list) and not value):
                continue
            parts.append(f"{name}={_canon(value)}")
        return f"{type(node).__name__}({','.join(parts)})"
    if isinstance(node, list):
        return "[" + ",".join(_canon(v) for v in node) + "]"
    return repr(node)


def _referenced(node: ast.AST) -> set[str]:
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


class Generator:
    """generate.py parsed once: its module-level definitions and its CASES."""

    def __init__(self, source: str):
        tree = _strip_docstrings(ast.parse(source))
        self.defs: dict[str, list[ast.stmt]] = {}
        self.cases: dict[str, ast.expr] = {}
        for stmt in tree.body:
            for name in _defined_names(stmt):
                self.defs.setdefault(name, []).append(stmt)
            if (isinstance(stmt, ast.Assign) and len(stmt.targets) == 1
                    and isinstance(stmt.targets[0], ast.Name)
                    and stmt.targets[0].id == "CASES"
                    and isinstance(stmt.value, ast.Dict)):
                for k, v in zip(stmt.value.keys, stmt.value.values):
                    if isinstance(k, ast.Constant) and isinstance(k.value, str):
                        self.cases[k.value] = v

    @classmethod
    def from_path(cls, path: Path) -> "Generator":
        return cls(Path(path).read_text())

    def closure(self, case: str) -> list[str]:
        """Sorted names of every module-level definition `case` reaches."""
        if case not in self.cases:
            raise KeyError(f"no CASES entry {case!r} in generate.py")
        seen: set[str] = set()
        todo = sorted(_referenced(self.cases[case]))
        while todo:
            name = todo.pop()
            if name in seen or name in _EXCLUDED or name not in self.defs:
                continue
            seen.add(name)
            for stmt in self.defs[name]:
                todo.extend(_referenced(stmt) - seen)
        return sorted(seen)

    def function(self, case: str) -> str:
        """The gen_* function a CASES entry calls, for the header."""
        v = self.cases[case]
        if isinstance(v, ast.Lambda):
            v = v.body
        if isinstance(v, ast.Call):
            v = v.func
        return v.id if isinstance(v, ast.Name) else ast.unparse(v)

    def code_digest(self, case: str) -> str:
        """sha256 of the case's own CASES entry and its closure, normalised."""
        h = hashlib.sha256()
        h.update(b"case\0" + case.encode() + b"\0")
        h.update(_canon(self.cases[case]).encode() + b"\0")
        for name in self.closure(case):
            h.update(b"def\0" + name.encode() + b"\0")
            for stmt in self.defs[name]:
                h.update(_canon(stmt).encode() + b"\0")
        return h.hexdigest()


def case_hash(gen: Generator, case: str, files: list[Path]) -> str:
    """The value of a reference's "hash:" line: its code digest followed by
    the bytes of each input file, in "files:" order."""
    h = hashlib.sha256()
    h.update(SCHEME.encode() + b"\0" + gen.code_digest(case).encode() + b"\0")
    for p in files:
        data = Path(p).read_bytes()
        h.update(str(len(data)).encode() + b"\0" + data)
    return h.hexdigest()
