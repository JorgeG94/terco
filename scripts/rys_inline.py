"""Inline Fortran text for the Rys roots lookup, the 2D recurrence and the
HRR, shared by gen_rys.py (energy) and gen_rys_grad.py (gradient) so every
generated Rys kernel is call-free.

WHY
---
Both generators splice the HGP scalar kernel and then replace the
scheme-specific middle with a Rys quadrature. Before this module existed,
that middle issued four `call`s per primitive quartet: one into
`trc_rys_roots{n}_tab` (module `trc_rys`), three into a per-axis
`rys2d_{tag}` helper, and three into `rys_hrr_{tag}` (energy) or
`rys_hrr_g_{tag}` (gradient). Even with the callee in the SAME module,
nvfortran did not inline them: Nsight Compute showed 17 non-inlined CALL
instructions and 255 registers on the (ps|ss) gradient kernel, with the
spilled registers showing up as local-memory load/store traffic. The fix is
textual, the same one the HGP kernels already use for the Boys evaluation
(`trc_boys`'s `boys_table`): put the callee's body directly in the caller,
reading only the table DATA from its module.

THE ROOTS LOOKUP
-----------------
Below is the same table lookup as `trc_rys.F90`'s `trc_rys_roots{n}_tab` --
same interval, same Chebyshev evaluation, same asymptotic branch -- with its
early `return` (illegal once inlined: it would return from the whole item
routine, not just the lookup) restructured into an if/else. Every operation
and its order is unchanged; only the control-flow shape is. It reads the
root/weight tables (`RYS_C{n}` etc.) from `trc_rys_table{n}` as module data,
exactly as the HGP kernels read `boys_table` from `trc_boys` -- the tables
stay where they are; only the evaluation code moves into the kernel body.

NAME COLLISIONS
----------------
Every pcsi{tag}/rysi{tag}/rysgi{tag} item routine gen_perclass.py's --split
template emits declares the SAME set of local names regardless of class (the
array bounds change, the names do not -- checked against (10|00) and
(22|22), the smallest and largest classes). Of trc_rys_roots{n}_tab's own
locals (`t`, `t2`, `d0`, `d1`, `tmp`, `fr`, `fw`, `xc`, `iv`, `s`, `j`,
`base`) and rys2d/rys_hrr's (`i`, `k`, `a`, `b`, `c`, `ab2`, `cd2`, and the
HRR's own intermediate array `t`), only the roots lookup's `t` collides (the
decoded primitive-pair index is already called `t`, as an integer) and the
subroutine dummies `u`/`w` collide (an int64 scalar and the digestion's
per-component weight are already called that -- renamed to `ru`/`rw` below,
which is also what both generators already call the arrays they pass the
lookup's result into). The HRR intermediate array is renamed `ht` to dodge
the same `t` collision. Everything else keeps its original name.

UNROLLED MODE
--------------
Inlining (above) removed the CALL boundary; it did not remove the runtime
array indexing inside the class kernel's component loops -- `hx(cnx0(aa),
cnx1(bb), ...)` with `aa`/`bb`/... as Fortran `do`-loop induction variables,
`gam(idx)`, the 2D-recurrence table `gx/gy/gz`, the HRR table `hx/hy/hz` and,
for the gradient, the per-primitive density blocks. A runtime index into a
small per-thread array forces nvfortran to keep that array in LOCAL memory
rather than registers, however small it is.

For classes with few enough Cartesian components (`--unroll-max-comp`, a
generator option; see gen_rys.py/gen_rys_grad.py), every (root, component)
combination is instead emitted as its own straight-line block, Python-side,
with literal indices: `rys2d_scalar_assignments`/`hrr_scalar_assignments`
emit the SAME algebra as `rys2d_inline_body`/`hrr_inline_body` but write a
named scalar (`nm`, e.g. `hx_1_0_2_0`) in place of each array element, so the
element a given `do rr`/`do aa` iteration used to visit at runtime is instead
baked into a distinct variable NAME at generation time -- there is no array
left for any index to be runtime. `scalar_decl_lines`/`unrolled_table_decls`
declare them. Every table element this project's own tables actually need is
emitted (checked against the array form's bounds, which were already exact,
not padded -- see each generator's own module docstring for the class-
specific accumulation this feeds). Classes above the threshold keep the
array/loop form (`rys2d_inline_body`/`hrr_inline_body`), unchanged.
"""

#: Binomial row b, for the horizontal transfer written out rather than
#: looped -- LMAX stays 2, so a shift of 0, 1 or 2 is everything this ever
#: needs. Index [b][t] is C(b, t).
BINOM = [[1], [1, 1], [1, 2, 1]]


def roots_table_use(nroots):
    """`use trc_rys_table{n}, only: ...` for the table data the inlined
    roots lookup below reads -- replaces `use trc_rys, only:
    trc_rys_roots{n}_tab`."""
    n = nroots
    return (f"   use trc_rys_table{n}, only: RYS_C{n}, RYS_ASY_R{n}, RYS_ASY_W{n}, "
            f"RYS_NROOTS{n}, &\n"
            f"                             RYS_DEGREE{n}, RYS_NINT{n}, RYS_WIDTH{n}, "
            f"RYS_XMAX{n}\n")


def roots_decls():
    """Declarations for the inlined roots lookup's scratch locals. `rlt` is
    `trc_rys_roots{n}_tab`'s `t`, renamed to dodge the decoded primitive-pair
    index already called `t`; `s`, `j`, `iv`, `base`, `xc`, `d0`, `d1`,
    `tmp`, `fr`, `fw`, `t2` are free in every class and kept as-is."""
    return ("      real(dp) :: rlt, t2, d0, d1, tmp, fr, fw, xc\n"
            "      integer  :: iv, s, j, base\n")


def roots_inline_body(nroots, x_expr, indent):
    """The body of `trc_rys_roots{n}_tab`, inlined: same table, same
    interval, same Chebyshev evaluation, same asymptotic branch, with the
    early `return` turned into an if/else (see module docstring). Writes
    `ru(1:nroots)`, `rw(1:nroots)` from `x_expr`."""
    n = nroots
    pad = " "*indent
    return f"""\
{pad}! Rys root/weight lookup, inlined from trc_rys_roots{n}_tab -- same
{pad}! table, interval and Chebyshev evaluation; the early return there
{pad}! becomes this if/else (a return here would exit the whole kernel).
{pad}if ({x_expr} >= RYS_XMAX{n}) then
{pad}   fr = 1.0_dp/{x_expr}
{pad}   fw = sqrt(fr)
{pad}   do s = 1, RYS_NROOTS{n}
{pad}      ru(s) = RYS_ASY_R{n}(s)*fr
{pad}      rw(s) = RYS_ASY_W{n}(s)*fw
{pad}   end do
{pad}else
{pad}   xc = max({x_expr}, 0.0_dp)
{pad}   iv = min(int(xc/RYS_WIDTH{n}) + 1, RYS_NINT{n})
{pad}   rlt = (2.0_dp*xc - (2*iv - 1)*RYS_WIDTH{n})/RYS_WIDTH{n}
{pad}   t2 = 2.0_dp*rlt
{pad}   do s = 1, 2*RYS_NROOTS{n}
{pad}      base = ((iv - 1)*2*RYS_NROOTS{n} + (s - 1))*(RYS_DEGREE{n} + 1)
{pad}      d0 = 0.0_dp
{pad}      d1 = 0.0_dp
{pad}      do j = RYS_DEGREE{n}, 1, -1
{pad}         tmp = d0
{pad}         d0 = t2*d0 - d1 + RYS_C{n}(base + j + 1)
{pad}         d1 = tmp
{pad}      end do
{pad}      tmp = rlt*d0 - d1 + RYS_C{n}(base + 1)
{pad}      if (s <= RYS_NROOTS{n}) then
{pad}         ru(s) = tmp
{pad}      else
{pad}         rw(s - RYS_NROOTS{n}) = tmp
{pad}      end if
{pad}   end do
{pad}end if
"""


def rys2d_inline_body(ni, nk, c00, cp00, b00, b10, bp01, g00, g, indent):
    """One axis of the 2D Rys recurrence, inlined -- same algebra as
    rys2d_{tag}, writing directly into `g` (dims 0:ni, 0:nk) instead of a
    dummy argument. `i`/`k` are the loop variables, declared once in the
    caller and reused for all three axes (sequential, not nested)."""
    pad = " "*indent
    lines = [f"{pad}{g}(0, 0) = {g00}"]
    if ni >= 1:
        lines.append(f"{pad}{g}(1, 0) = {c00}*{g00}")
    if ni >= 2:
        lines.append(f"{pad}do i = 1, {ni - 1}")
        lines.append(f"{pad}   {g}(i + 1, 0) = real(i, dp)*{b10}*{g}(i - 1, 0) + {c00}*{g}(i, 0)")
        lines.append(f"{pad}end do")
    if nk >= 1:
        lines.append(f"{pad}{g}(0, 1) = {cp00}*{g00}")
    if nk >= 1 and ni >= 1:
        lines.append(f"{pad}do i = 1, {ni}")
        lines.append(f"{pad}   {g}(i, 1) = real(i, dp)*{b00}*{g}(i - 1, 0) + {cp00}*{g}(i, 0)")
        lines.append(f"{pad}end do")
    if nk >= 2:
        lines.append(f"{pad}do k = 1, {nk - 1}")
        lines.append(f"{pad}   {g}(0, k + 1) = real(k, dp)*{bp01}*{g}(0, k - 1) + {cp00}*{g}(0, k)")
        if ni >= 1:
            lines.append(f"{pad}   do i = 1, {ni}")
            lines.append(f"{pad}      {g}(i, k + 1) = real(i, dp)*{b00}*{g}(i - 1, k) &")
            lines.append(f"{pad}                  + real(k, dp)*{bp01}*{g}(i, k - 1) + {cp00}*{g}(i, k)")
            lines.append(f"{pad}   end do")
        lines.append(f"{pad}end do")
    return "\n".join(lines) + "\n"


def hrr_inline_body(ni, nk, la, lb, lc, ld, g, ab, cd, h, indent):
    """The horizontal transfer on both indices, inlined -- same algebra as
    rys_hrr_{tag}/rys_hrr_g_{tag}: (a,b+1| = (a+1,b| + (A-B)(a,b|, binomials
    for a shift of 0, 1 or 2 written out. Writes directly into `h` from `g`
    (dims 0:ni, 0:nk), through the shared scratch array `ht` (dims 0:la,
    0:lb, 0:nk) -- renamed from rys_hrr_{tag}'s local `t` to dodge the
    decoded primitive-pair index already called `t` in the caller. `ht`,
    `ab2`, `cd2` and the loop variables `a`/`b`/`c` are declared once in the
    caller and reused for all three axes (sequential, not nested)."""
    pad = " "*indent
    lines = [f"{pad}ab2 = {ab}*{ab}", f"{pad}cd2 = {cd}*{cd}"]
    lines.append(f"{pad}do a = 0, {la}")
    lines.append(f"{pad}   do c = 0, {nk}")
    for b in range(lb + 1):
        terms = []
        for t in range(b + 1):
            coef, p = BINOM[b][t], b - t
            pre = "" if coef == 1 else f"{float(coef)}_dp*"
            mul = "" if p == 0 else (f"{ab}*" if p == 1 else "ab2*")
            terms.append(f"{pre}{mul}{g}(a + {t}, c)")
        lines.append(f"{pad}      ht(a, {b}, c) = " + " + ".join(terms))
    lines.append(f"{pad}   end do")
    lines.append(f"{pad}end do")
    lines.append(f"{pad}do a = 0, {la}")
    lines.append(f"{pad}   do b = 0, {lb}")
    lines.append(f"{pad}      do c = 0, {lc}")
    for d in range(ld + 1):
        terms = []
        for t in range(d + 1):
            coef, p = BINOM[d][t], d - t
            pre = "" if coef == 1 else f"{float(coef)}_dp*"
            mul = "" if p == 0 else (f"{cd}*" if p == 1 else "cd2*")
            terms.append(f"{pre}{mul}ht(a, b, c + {t})")
        lines.append(f"{pad}         {h}(a, b, c, {d}) = " + " + ".join(terms))
    lines.append(f"{pad}      end do")
    lines.append(f"{pad}   end do")
    lines.append(f"{pad}end do")
    return "\n".join(lines) + "\n"


def nm(prefix, *idx):
    """A named-scalar Fortran identifier for a table element -- e.g.
    `nm("hx", 1, 0, 2, 0)` is `hx_1_0_2_0`. Used by the unrolled emission
    mode (see the module docstring's UNROLLED MODE section below) in place
    of an array element `hx(1, 0, 2, 0)`: the index is baked into the NAME
    at generation time, so there is no runtime index for nvfortran to force
    into local memory."""
    return prefix + "".join(f"_{i}" for i in idx)


def rys2d_scalar_names(ni, nk, prefix):
    """Every named scalar `rys2d_scalar_assignments` below assigns, in no
    particular order -- for declarations."""
    return [nm(prefix, i, k) for k in range(nk + 1) for i in range(ni + 1)]


def rys2d_scalar_assignments(ni, nk, c00, cp00, b00, b10, bp01, g00, prefix, indent):
    """The unrolled-mode twin of `rys2d_inline_body`: same algebra, same
    traversal order (so the dependency chain a flat sequence of assignments
    needs is respected automatically -- each g(i,k) is assigned only after
    every g it reads), but every `g(i, k)` becomes a NAMED SCALAR
    `{prefix}_{i}_{k}` instead of an array element addressed by the runtime
    loop variables `i`/`k`. `real(i, dp)`/`real(k, dp)` fold to their literal
    value here, since i and k are Python ints, not Fortran loop variables."""
    pad = " "*indent
    g = lambda i, k: nm(prefix, i, k)
    lines = [f"{pad}{g(0, 0)} = {g00}"]
    if ni >= 1:
        lines.append(f"{pad}{g(1, 0)} = {c00}*{g(0, 0)}")
    if ni >= 2:
        for i in range(1, ni):
            lines.append(f"{pad}{g(i + 1, 0)} = {float(i)}_dp*{b10}*{g(i - 1, 0)} + {c00}*{g(i, 0)}")
    if nk >= 1:
        lines.append(f"{pad}{g(0, 1)} = {cp00}*{g(0, 0)}")
    if nk >= 1 and ni >= 1:
        for i in range(1, ni + 1):
            lines.append(f"{pad}{g(i, 1)} = {float(i)}_dp*{b00}*{g(i - 1, 0)} + {cp00}*{g(i, 0)}")
    if nk >= 2:
        for k in range(1, nk):
            lines.append(f"{pad}{g(0, k + 1)} = {float(k)}_dp*{bp01}*{g(0, k - 1)} + {cp00}*{g(0, k)}")
            if ni >= 1:
                for i in range(1, ni + 1):
                    lines.append(f"{pad}{g(i, k + 1)} = {float(i)}_dp*{b00}*{g(i - 1, k)} "
                                 f"+ {float(k)}_dp*{bp01}*{g(i, k - 1)} + {cp00}*{g(i, k)}")
    return "\n".join(lines) + "\n"


def hrr_scalar_names(la, lb, lc, ld, prefix):
    """Every named `h` scalar the HRR unroll below assigns -- for
    declarations."""
    return [nm(prefix, a, b, c, d)
            for a in range(la + 1) for b in range(lb + 1)
            for c in range(lc + 1) for d in range(ld + 1)]


def ht_scalar_names(la, lb, nk):
    """Every named `ht` intermediate the HRR unroll below assigns -- shared
    across the x/y/z axis calls (same bounds every time), exactly as the
    array form shares one `ht` array across all three."""
    return [nm("ht", a, b, c)
            for a in range(la + 1) for b in range(lb + 1) for c in range(nk + 1)]


def hrr_scalar_assignments(ni, nk, la, lb, lc, ld, g, ab, cd, h, indent):
    """The unrolled-mode twin of `hrr_inline_body`: same two-stage binomial
    transfer, same traversal order, but `g(i, k)`/`ht(a, b, c)`/`h(a, b, c,
    d)` all become named scalars (`nm`) instead of array elements. `ab2`/
    `cd2` stay true scalars, recomputed per axis call exactly as the array
    form does -- declared once by the caller, assigned here, reused x/y/z
    sequentially."""
    pad = " "*indent
    gf = lambda i, k: nm(g, i, k)
    ht = lambda a, b, c: nm("ht", a, b, c)
    hf = lambda a, b, c, d: nm(h, a, b, c, d)
    lines = [f"{pad}ab2 = {ab}*{ab}", f"{pad}cd2 = {cd}*{cd}"]
    for a in range(la + 1):
        for c in range(nk + 1):
            for b in range(lb + 1):
                terms = []
                for t in range(b + 1):
                    coef, p = BINOM[b][t], b - t
                    pre = "" if coef == 1 else f"{float(coef)}_dp*"
                    mul = "" if p == 0 else (f"{ab}*" if p == 1 else "ab2*")
                    terms.append(f"{pre}{mul}{gf(a + t, c)}")
                lines.append(f"{pad}{ht(a, b, c)} = " + " + ".join(terms))
    for a in range(la + 1):
        for b in range(lb + 1):
            for c in range(lc + 1):
                for d in range(ld + 1):
                    terms = []
                    for t in range(d + 1):
                        coef, p = BINOM[d][t], d - t
                        pre = "" if coef == 1 else f"{float(coef)}_dp*"
                        mul = "" if p == 0 else (f"{cd}*" if p == 1 else "cd2*")
                        terms.append(f"{pre}{mul}{ht(a, b, c + t)}")
                    lines.append(f"{pad}{hf(a, b, c, d)} = " + " + ".join(terms))
    return "\n".join(lines) + "\n"


def scalar_decl_lines(names, indent):
    """Fortran declarations for a flat list of scalar names, batched so no
    single `real(dp) ::` statement grows unboundedly: at most 6 names per
    continuation line, at most 10 lines (60 names) per statement. Returns ''
    for an empty list (no `real(dp) ::` with nothing after it)."""
    if not names:
        return ""
    pad = " "*indent
    per_line, stmt_lines = 6, 10
    per_stmt = per_line*stmt_lines
    out = []
    for s in range(0, len(names), per_stmt):
        chunk = names[s:s + per_stmt]
        rows = [chunk[i:i + per_line] for i in range(0, len(chunk), per_line)]
        stmt = []
        for i, row in enumerate(rows):
            row_txt = ", ".join(row)
            cont = ", &" if i < len(rows) - 1 else ""
            head = f"{pad}real(dp) :: " if i == 0 else f"{pad}             "
            stmt.append(f"{head}{row_txt}{cont}")
        out.append("\n".join(stmt))
    return "\n".join(out) + "\n"


def unrolled_table_decls(ni, nk, la, lb, lc, ld, indent):
    """Declarations for the unrolled mode's 2D-recurrence and HRR scratch:
    `ab2`/`cd2` (true scalars, as in the array form) plus every named gx/gy/
    gz, hx/hy/hz and the shared ht scalar, for the table sizes (ni, nk, la,
    lb, lc, ld) the caller is using -- already raised for the gradient's
    extended table, same as `rys2d_hrr_decls`."""
    pad = " "*indent
    g_names, h_names = [], []
    for ax in ("x", "y", "z"):
        g_names += rys2d_scalar_names(ni, nk, f"g{ax}")
        h_names += hrr_scalar_names(la, lb, lc, ld, f"h{ax}")
    ht_names = ht_scalar_names(la, lb, nk)
    return (f"{pad}real(dp) :: ab2, cd2\n"
            + scalar_decl_lines(g_names, indent)
            + scalar_decl_lines(h_names, indent)
            + scalar_decl_lines(ht_names, indent))


def rys2d_hrr_decls(la, lb, nk_h):
    """Declarations shared by the inlined 2D recurrence and HRR: the loop
    variables and the HRR's scratch array `ht`, sized for one axis -- the
    same array is reused for x, y, z in turn. `nk_h` is the HRR's own
    second dimension bound (equal to `nk` for the energy kernel, `nk + 1`
    for the gradient kernel's extended table)."""
    return (f"      integer  :: i, k, a, b, c\n"
            f"      real(dp) :: ab2, cd2\n"
            f"      real(dp) :: ht(0:{la}, 0:{lb}, 0:{nk_h})\n")
