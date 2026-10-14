#!/usr/bin/env python3
"""
Emit two- and three-centre Coulomb kernels for density fitting.

    (a|c)   = int int  chi_a(r1) (1/r12) chi_c(r2)
    (ab|c)  = int int  phi_a(r1) phi_b(r1) (1/r12) chi_c(r2)

WHY THIS IS THE CHEAP ONE
-------------------------
Nothing new is needed. The four-centre VRR already builds [a0|c0]^(m), which
IS the two-centre quantity; the four-centre form only differs by having a
Gaussian PRODUCT on each side instead of a single function. A single shell is
formally a product with a unit function of zero exponent, so:

    2-centre:  zeta = alpha,      P = A,  K_ab = 1
               eta  = gamma,      Q = C,  K_cd = 1
    3-centre:  zeta = alpha+beta, P as usual, K_ab = exp(-xi |AB|^2)
               eta  = gamma,      Q = C,  K_cd = 1

and the base case is the same one the ERI uses,

    [00|00]^(m) = 2 pi^(5/2) / (zeta eta sqrt(zeta+eta)) K_ab K_cd F_m(rho R_PQ^2)

with rho = zeta*eta/(zeta+eta). The recurrence coefficients are untouched. The
only structural difference is that the ket carries no angular momentum to
transfer, so the ket HRR disappears and only the bra HRR runs (3-centre), or
neither does (2-centre).

So this reuses gen_vrr's sieved emitter verbatim -- including Gill's sieve --
and supplies its own target set: the ket index is pinned to |c| = lc rather
than ranging over [lc, lc+ld].

These produce INTEGRALS, not a Fock matrix. Density fitting wants the (mn|P)
tensor itself, so there is no digestion here.

ONE MODULE PER CLASS
---------------------
At --lmax 2 --lmax-aux 3 (the shipped set) every class in one module was
0.9 MB; at --lmax 3 --lmax-aux 4 (f orbitals, g auxiliary) the same scheme is
7.1 MB and nvfortran's GPU back end does not finish it in any reasonable time.
gen_perclass.py hit the identical wall for the four-centre kernels and the
fix there is the one used here: one module per (la,lb|lc) or (la|lc) class,
each owning its OWN `do concurrent` launch with the class body inlined in the
same translation unit (stdpar does not inline `!$acc routine seq` across a
module boundary -- see gen_perclass.py's docstring and trc_api.F90's -- so the
item routine that calls the class kernel has to share its file). The host
side groups shell pairs by (la,lb) and auxiliary shells by lc into index
LISTS (scripts/gen_perclass.py's --split groups by PHYSICALLY PERMUTING its
segment list instead; a list is simpler here because the primitive-pair pool
is already addressed by offset, so nothing downstream needs reordering), then
a small dispatcher -- host code, free to cross modules -- calls each
non-empty class's own driver once.

ABOVE THE CEILING: RYS
-----------------------
One module per class moves the wall; it does not remove it. At --lmax 3
--lmax-aux 4, single-file compiles with real GPU flags measured (g|g) 2c at
1m43s, (fd|g) 3c at 22m, and (ff|f)/(ff|g) 3c killed unfinished past 30-50m
-- the written-out VRR+HRR's own code size, not the one-module-per-class
split, is what nvfortran chokes on past a shell or two of f/g. Classes past
SHIPPED_L/SHIPPED_LA (below this module's own `main`) use a RYS QUADRATURE
body instead (`emit_class_rys`): the same algebra gen_rys.py already uses
for the four-centre kernels above d, with the ket pinned to a single
primitive (eta = gamma, Q = C, ld = 0 -- a three- or two-centre integral's
ket is a lone shell, same observation the VRR/HRR body above already made).
Reusing scripts/rys_inline.py rather than writing a second quadrature keeps
the two Rys bodies (four-centre and this one) sharing one set of inlined
roots/2D-recurrence/HRR emitters -- a bug fixed in one is fixed in both.
Classes within the ceiling are untouched: `default_use_rys` returns False
for every one of them regardless of --lmax/--lmax-aux, so the shipped
--lmax 2 --lmax-aux 3 set keeps its written-out VRR+HRR, byte for byte.
`--force-rys` (TESTING ONLY, never passed by src/CMakeLists.txt's normal
configuration) overrides the policy per class, for comparing the two bodies
against each other on classes where both exist.

USAGE
    gen_df.py --lmax 3 --lmax-aux 4 --split src/generated
"""

import argparse
import importlib.util
import os

import cartorder


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HERE = os.path.dirname(os.path.abspath(__file__))
gv = _load("gv", os.path.join(HERE, "gen_vrr.py"))
gh = _load("gh", os.path.join(HERE, "gen_hrr.py"))
ri = _load("ri", os.path.join(HERE, "rys_inline.py"))


def cart(l, order="lexicographic"):
    """(nx,ny,nz) for one shell, in `order`.

    Thin wrapper: the order itself is the single definition in cartorder.py.
    Called with NO explicit order (the default) wherever the result only
    feeds `I`, the VRR sieve's internal index map (built from `gv.cart_list`,
    itself always lexicographic) -- a dict lookup by triple does not care what
    order produced the triple, so building `targets` this way is safe
    regardless. Only `emit_bra_hrr`'s and the two-centre body's OUTPUT
    enumerations (the AO/tensor slot order) take `order` explicitly.
    """
    return cartorder.components(l, order)


def binom(n, k):
    r = 1
    for i in range(1, k + 1):
        r = r*(n - k + i)//i
    return r


def powexpr(var, n):
    return "*".join([var]*n) if n else ""


def emit_bra_hrr(la, lb, lc, I, nca, order="lexicographic"):
    """(ab|c) from [e0|c0], e in [la, la+lb]. Ket untouched: ld = 0.

    OUTPUT: `slot` below is the (ab|c) tensor's own AO axis, so the three
    outer loops follow `order`. `I[...]` addresses `g`, the VRR's internal
    intermediate, and stays lexicographic regardless (see `cart`'s
    docstring).
    """
    o, w = [], None
    w = o.append
    na, nb = len(cart(la)), len(cart(lb))
    for ic, cp in enumerate(cart(lc, order)):
        jc = I[cp]
        for ib, (bx, by, bz) in enumerate(cart(lb, order)):
            for ia, (ax, ay, az) in enumerate(cart(la, order)):
                terms = []
                for px in range(bx + 1):
                    for py in range(by + 1):
                        for pz in range(bz + 1):
                            cbin = binom(bx, px)*binom(by, py)*binom(bz, pz)
                            ja = I[(ax + px, ay + py, az + pz)]
                            fac = [f for f in (powexpr("abx", bx - px),
                                               powexpr("aby", by - py),
                                               powexpr("abz", bz - pz)) if f]
                            head = "" if cbin == 1 else f"{float(cbin)}_dp*"
                            body = "*".join(fac)
                            if body:
                                body += "*"
                            terms.append(f"{head}{body}g({ja + (jc - 1)*nca})")
                expr = terms[0]
                for t in terms[1:]:
                    expr += f" &\n                  + {t}"
                slot = 1 + ia + na*(ib + nb*ic)
                w(f"               out({slot}) = out({slot}) + {expr}")
    return "\n".join(o)


def emit_class(la, lb, lc, three, order="lexicographic"):
    """One class. `three` selects (ab|c); otherwise (a|c)."""
    lab = la + lb
    lt = lab + lc
    nca, ncc = gv.ncum(lab), gv.ncum(lc)
    na, nb, nc = len(cart(la)), len(cart(lb)), len(cart(lc))
    nout = na*nb*nc if three else na*nc
    nv = nca*ncc

    P = [None] + gv.cart_list(max(lab, lc))
    I = {p: i for i, p in enumerate(P) if p is not None}

    # Targets: the ket index is PINNED to |c| = lc, since there is no ket
    # transfer. That is the whole difference from the four-centre target set,
    # and it lets the sieve prune harder.
    kets = [I[p] for p in cart(lc)]
    if three:
        bras = [I[p] for l in range(la, lab + 1) for p in cart(l)]
    else:
        bras = [I[p] for p in cart(la)]
    targets = {(b, c, 0) for b in bras for c in kets}

    vrr, fin, nlive = gv.emit_class_sieved(lab, lc, targets, P, I, nca)
    vrr = "\n".join("         " + ln if ln.strip() else ln
                    for ln in vrr.split("\n"))

    if three:
        tag = f"3c{la}{lb}{lc}"
        body = emit_bra_hrr(la, lb, lc, I, nca, order)
        #
        # The bra arrives as PRECOMPUTED primitive-pair data, not as two sets
        # of exponents. This kernel is called once per (shell pair, auxiliary
        # shell); rebuilding zeta, P and K_ab from the exponents here would
        # redo the bra product once per auxiliary SHELL, which for a real
        # fitting basis is hundreds of times over. The pair list computes it
        # once per shell pair and this reads it.
        #
        sig = ("npp, ppz, ppk, ppp, ppa, abx, aby, abz, "
               "nprim_c, ec, cc_, rc, out")
        decls = """      integer,  intent(in)  :: npp, nprim_c
      real(dp), intent(in)  :: ppz(npp), ppk(npp), ppp(3, npp), ppa(3, npp)
      real(dp), intent(in)  :: abx, aby, abz
      real(dp), intent(in)  :: ec(nprim_c), cc_(nprim_c), rc(3)"""
        bra_loop_open = """      do ip = 1, npp
         zeta = ppz(ip)
         kab  = ppk(ip)
         px = ppp(1, ip); py = ppp(2, ip); pz = ppp(3, ip)
         pax = ppa(1, ip); pay = ppa(2, ip); paz = ppa(3, ip)"""
        bra_loop_close = ""
        abset = ""
    else:
        tag = f"2c{la}{lc}"
        # No bra transfer either: the surviving [a0|c0] IS the answer.
        # OUTPUT: `slot` is the (a|c) tensor's own AO axis, so this follows
        # `order`; `I[...]` addresses `g` and stays lexicographic.
        o = []
        for ic, cp in enumerate(cart(lc, order)):
            jc = I[cp]
            for ia, ap in enumerate(cart(la, order)):
                ja = I[ap]
                slot = 1 + ia + na*ic
                o.append(f"               out({slot}) = out({slot})"
                         f" + g({ja + (jc - 1)*nca})")
        body = "\n".join(o)
        sig_extra = ""
        decl_extra = ""
        # Two-centre keeps the plain exponent interface: there is no pair to
        # precompute, and it is called once per (shell, shell) rather than
        # once per (pair, shell), so nothing is recomputed.
        sig = "nprim_a, ea, ca, ra, nprim_c, ec, cc_, rc, out"
        decls = """      integer,  intent(in)  :: nprim_a, nprim_c
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a), ra(3)
      real(dp), intent(in)  :: ec(nprim_c), cc_(nprim_c), rc(3)"""
        bra_loop_open = """      do ip = 1, nprim_a
         al = ea(ip)
         zeta = al
         kab = ca(ip)
         px = ra(1); py = ra(2); pz = ra(3)
         pax = 0.0_dp; pay = 0.0_dp; paz = 0.0_dp"""
        bra_loop_close = ""
        abset = ""   # no HRR in the two-centre case, so no AB vector

    return f"""
   pure subroutine df_{tag}({sig})
      !! {'(ab|c)' if three else '(a|c)'} class ({la}{lb}|{lc}) -- {nlive} VRR statements after the sieve.
      !$acc routine seq
{decls}
      real(dp), intent(out) :: out({nout})

      integer  :: ip, kp, k, cur
      real(dp) :: al, gam, zeta, eta, zpe, rho, kab, pref
      real(dp) :: px, py, pz, pqx, pqy, pqz, tval, wc
      real(dp) :: pax, pay, paz, qcx, qcy, qcz
      real(dp) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re
      real(dp) :: v({nv},0:1), g({nv})
      real(dp) :: f(0:BOYS_MMAX)
{abset}

      do k = 1, {nout}
         out(k) = 0.0_dp
      end do

{bra_loop_open}

         do k = 1, {nv}
            g(k) = 0.0_dp
         end do

         do kp = 1, nprim_c
            gam = ec(kp)
            eta = gam
            zpe = zeta + eta
            rho = zeta*eta/zpe
            ! The ket is a single shell, so Q = C and QC = 0. Kept as
            ! variables rather than folded away because the VRR is emitted
            ! by the four-centre generator and expects them.
            qcx = 0.0_dp; qcy = 0.0_dp; qcz = 0.0_dp
            pqx = px - rc(1); pqy = py - rc(2); pqz = pz - rc(3)
            wc = (zeta*px + eta*rc(1))/zpe; wpx = wc - px; wqx = wc - rc(1)
            wc = (zeta*py + eta*rc(2))/zpe; wpy = wc - py; wqy = wc - rc(2)
            wc = (zeta*pz + eta*rc(3))/zpe; wpz = wc - pz; wqz = wc - rc(3)
            oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = 0.5_dp/zpe
            rz = rho/zeta; re = rho/eta
            tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)
            ! K_cd = 1: a single Gaussian has no product exponential.
            pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*kab*cc_(kp)
            call boys_eval({lt}, tval, f)

{vrr}

            do k = 1, {nv}
               g(k) = g(k) + v(k,{fin})
            end do
         end do

{body}
      end do
{bra_loop_close}
   end subroutine df_{tag}
"""


#
# The shipped ceiling: the (--lmax, --lmax-aux) the committed src/generated
# DF set was built at. NOT args.lmax/args.lmax_aux below, which is whatever
# a non-default invocation asks for -- these two stay fixed so the classes
# the "byte for byte" promise covers (docs/GENERATORS.md) never depend on a
# LATER ceiling. A 3-centre class needs Rys once la, lb or lc goes past this;
# a 2-centre one (both sides drawn from the AUXILIARY ceiling -- main()'s
# loop runs la, lc over 0..LA for 2-centre) once la or lc does. See this
# module's docstring for the nvfortran compile-time wall that makes this
# necessary: (fd|g) 3c took 22 minutes to compile at LMAX=3/LMAX_AUX=4;
# (ff|f) and (ff|g) did not finish in 30-50 minutes.
#
SHIPPED_L = 2
SHIPPED_LA = 3


def default_use_rys(la, lb, lc, three):
    """The un-forced policy -- True once a class needs a shell this
    project's GPU back end cannot compile a written-out VRR+HRR for in any
    reasonable time. `--force-rys` (main()) can override this per class,
    for the Rys-vs-HGP correctness comparison and other testing; the
    shipped configuration (src/CMakeLists.txt) never passes it, so this
    function is the only thing that decides for a normal build."""
    if three:
        return la > SHIPPED_L or lb > SHIPPED_L or lc > SHIPPED_LA
    return la > SHIPPED_LA or lc > SHIPPED_LA


def emit_class_rys(la, lb, lc, three, order, unroll_max):
    """One class, RYS QUADRATURE body -- the above-ceiling twin of
    `emit_class`. Same public signature (df_{tag}, same dummy arguments)
    and the same digestion convention (`out` slot = 1 + ia + na*(ib + nb*ic)
    for (ab|c), 1 + ia + na*ic for (a|c), both cartorder.py's own ordering),
    so a caller cannot tell which body a given df_{tag} runs; only the
    integral itself differs.

    ALGEBRA: gen_rys.py's four-centre Rys middle (scripts/rys_inline.py's
    inlined roots lookup, 2D recurrence and horizontal transfer), with the
    ket pinned to a single primitive -- eta = gamma, Q = C, ld = 0
    throughout (see gen_df.py's module docstring), so there is no ket HRR
    and the ket's own "shell" is always the one-component s set (cart(0) =
    [(0,0,0)]), addressed by the literal index 0 rather than a cnx3/cny3/
    cnz3 table -- unlike the four-centre generator, which keeps ld generic.
    nroots = (la+lb+lc)//2 + 1, the general nroots = (la+lb+lc+ld)//2 + 1
    with ld = 0.

    `out` IS the four-centre splice's `vbuf`, not `emit_class`'s `g`: the
    horizontal transfer is linear in each axis's 2D table, but the final
    digestion is the PRODUCT hx*hy*hz across the three axes, so the
    per-root, per-primitive-quartet contribution has to land directly in
    the output block. There is no g-then-HRR two-stage deferral the way
    `emit_class`'s fully-3D VRR intermediate allows (its `g`, summed over
    the auxiliary primitive kp, then transferred once per bra primitive ip
    by `emit_bra_hrr`) -- a per-axis 1D quadrature only becomes additive
    AFTER the x/y/z product, which is why this zeroes `out` once (like
    `emit_class`) but accumulates into it inside BOTH primitive loops and
    the root loop, like gen_rys.py's four-centre middle does into `vbuf`.

    Heavy classes -- this function's entire reason to exist -- stay in the
    array/loop form (rys_inline.py's UNROLLED MODE); `unroll_max` is the
    SAME na*nb*nc <= threshold gen_rys.py's `--unroll-max-comp` already
    uses for na*nb*nc*nd (nd = 1 here, ld = 0), not a new one -- see this
    module's main() for where it comes from (CMake's TERCO_RYS_UNROLL_MAX,
    threaded through unchanged)."""
    ni, nk = la + lb, lc
    nroots = (la + lb + lc)//2 + 1
    na, nb, nc = len(cart(la)), len(cart(lb)), len(cart(lc))
    nout = na*nb*nc if three else na*nc
    unroll = na*nb*nc <= unroll_max

    if three:
        tag = f"3c{la}{lb}{lc}"
        sig = ("npp, ppz, ppk, ppp, ppa, abx, aby, abz, "
               "nprim_c, ec, cc_, rc, out")
        decls = """      integer,  intent(in)  :: npp, nprim_c
      real(dp), intent(in)  :: ppz(npp), ppk(npp), ppp(3, npp), ppa(3, npp)
      real(dp), intent(in)  :: abx, aby, abz
      real(dp), intent(in)  :: ec(nprim_c), cc_(nprim_c), rc(3)"""
        bra_loop_open = """      do ip = 1, npp
         zeta = ppz(ip)
         kab  = ppk(ip)
         px = ppp(1, ip); py = ppp(2, ip); pz = ppp(3, ip)
         pax = ppa(1, ip); pay = ppa(2, ip); paz = ppa(3, ip)"""
        ab_axes = ("abx", "aby", "abz")
    else:
        tag = f"2c{la}{lc}"
        sig = "nprim_a, ea, ca, ra, nprim_c, ec, cc_, rc, out"
        decls = """      integer,  intent(in)  :: nprim_a, nprim_c
      real(dp), intent(in)  :: ea(nprim_a), ca(nprim_a), ra(3)
      real(dp), intent(in)  :: ec(nprim_c), cc_(nprim_c), rc(3)"""
        bra_loop_open = """      do ip = 1, nprim_a
         zeta = ea(ip)
         kab = ca(ip)
         px = ra(1); py = ra(2); pz = ra(3)
         pax = 0.0_dp; pay = 0.0_dp; paz = 0.0_dp"""
        ab_axes = ("0.0_dp", "0.0_dp", "0.0_dp")  # lb = 0: never multiplied

    decl_new = (
        # `oz`, not `k`: rys2d_hrr_decls (below, loop form) declares its own
        # `k` for the 2D-recurrence axis loop, and a second `integer :: k`
        # here would be a duplicate declaration in the same scope -- caught
        # by compiling a forced-Rys loop-form class, which is also the ONLY
        # form that reaches rys2d_hrr_decls (the unrolled branch never
        # declares `k`, so it would not have caught this).
        "      integer  :: ip, kp, oz\n"
        "      real(dp) :: zeta, eta, zpe, rho, kab, pref\n"
        "      real(dp) :: px, py, pz, pqx, pqy, pqz, tval\n"
        "      real(dp) :: pax, pay, paz, qcx, qcy, qcz\n"
        f"      real(dp) :: ru({nroots}), rw({nroots})\n"
        "      real(dp) :: u2, dinv, b00r, b10r, bp01r, f00\n"
        "      real(dp) :: c00x, c00y, c00z, cp0x, cp0y, cp0z\n"
        + ri.roots_decls())
    tabs = ""
    if unroll:
        decl_new += ri.unrolled_table_decls(ni, nk, la, lb, lc, 0, 6)
    else:
        decl_new += (
            f"      real(dp) :: gx(0:{ni}, 0:{nk}), gy(0:{ni}, 0:{nk}), gz(0:{ni}, 0:{nk})\n"
            f"      real(dp) :: hx(0:{la}, 0:{lb}, 0:{lc}, 0:0)\n"
            f"      real(dp) :: hy(0:{la}, 0:{lb}, 0:{lc}, 0:0)\n"
            f"      real(dp) :: hz(0:{la}, 0:{lb}, 0:{lc}, 0:0)\n"
            "      integer  :: rr, aa, bb, cc, idx\n"
            + ri.rys2d_hrr_decls(la, lb, 0, nk))
        for cnm, ax in (("cnx", 0), ("cny", 1), ("cnz", 2)):
            for which, l in enumerate((la, lb, lc)):
                comps = cart(l, order)
                vals = ", ".join(str(c[ax]) for c in comps)
                tabs += f"   integer, parameter :: {cnm}{which}_{tag}({len(comps)}) = [{vals}]\n"

    kp_head = f"""         do kp = 1, nprim_c
            eta = ec(kp)
            zpe = zeta + eta
            rho = zeta*eta/zpe
            ! The ket is a single shell, so Q = C and QC = 0 -- same as
            ! emit_class's HGP body; ld = 0 throughout, so there is no ket
            ! HRR either.
            qcx = 0.0_dp; qcy = 0.0_dp; qcz = 0.0_dp
            pqx = px - rc(1); pqy = py - rc(2); pqz = pz - rc(3)
            tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)
            ! K_cd = 1: a single Gaussian has no product exponential.
            pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*kab*cc_(kp)
"""

    if unroll:
        comps = [cart(l, order) for l in (la, lb, lc)]
        pad = " "*12
        mid = [ri.roots_inline_body(nroots, "tval", 12)]
        for rr in range(1, nroots + 1):
            mid.append(f"""\
            u2 = ru({rr})/(1.0_dp - ru({rr}))*rho
            dinv = 1.0_dp/(zeta*eta + u2*(zeta + eta))
            b00r  = 0.5_dp*u2*dinv
            b10r  = 0.5_dp*(eta + u2)*dinv
            bp01r = 0.5_dp*(zeta + u2)*dinv
            f00 = pref*rw({rr})
            c00x = pax - u2*eta*dinv*pqx
            c00y = pay - u2*eta*dinv*pqy
            c00z = paz - u2*eta*dinv*pqz
            cp0x = qcx + u2*zeta*dinv*pqx
            cp0y = qcy + u2*zeta*dinv*pqy
            cp0z = qcz + u2*zeta*dinv*pqz
""")
            mid.append(ri.rys2d_scalar_assignments(ni, nk, "c00x", "cp0x", "b00r", "b10r", "bp01r", "1.0_dp", "gx", 12))
            mid.append(ri.rys2d_scalar_assignments(ni, nk, "c00y", "cp0y", "b00r", "b10r", "bp01r", "1.0_dp", "gy", 12))
            mid.append(ri.rys2d_scalar_assignments(ni, nk, "c00z", "cp0z", "b00r", "b10r", "bp01r", "f00",    "gz", 12))
            mid.append(ri.hrr_scalar_assignments(ni, nk, la, lb, lc, 0, "gx", ab_axes[0], "0.0_dp", "hx", 12))
            mid.append(ri.hrr_scalar_assignments(ni, nk, la, lb, lc, 0, "gy", ab_axes[1], "0.0_dp", "hy", 12))
            mid.append(ri.hrr_scalar_assignments(ni, nk, la, lb, lc, 0, "gz", ab_axes[2], "0.0_dp", "hz", 12))
            idx = 0
            acc = []
            for icc in range(1, nc + 1):
                for ibb in range(1, nb + 1):
                    for iaa in range(1, na + 1):
                        idx += 1
                        eix, eiy, eiz = comps[0][iaa - 1]
                        ejx, ejy, ejz = comps[1][ibb - 1]
                        ekx, eky, ekz = comps[2][icc - 1]
                        hxn = ri.nm("hx", eix, ejx, ekx, 0)
                        hyn = ri.nm("hy", eiy, ejy, eky, 0)
                        hzn = ri.nm("hz", eiz, ejz, ekz, 0)
                        acc.append(f"{pad}out({idx}) = out({idx}) + {hxn}*{hyn}*{hzn}")
            mid.append("\n".join(acc) + "\n")
        kp_body = "".join(mid)
    else:
        kp_body = f"""\
{ri.roots_inline_body(nroots, "tval", 12)}
            do rr = 1, {nroots}
               u2 = ru(rr)/(1.0_dp - ru(rr))*rho
               dinv = 1.0_dp/(zeta*eta + u2*(zeta + eta))
               b00r  = 0.5_dp*u2*dinv
               b10r  = 0.5_dp*(eta + u2)*dinv
               bp01r = 0.5_dp*(zeta + u2)*dinv
               f00 = pref*rw(rr)
               c00x = pax - u2*eta*dinv*pqx
               c00y = pay - u2*eta*dinv*pqy
               c00z = paz - u2*eta*dinv*pqz
               cp0x = qcx + u2*zeta*dinv*pqx
               cp0y = qcy + u2*zeta*dinv*pqy
               cp0z = qcz + u2*zeta*dinv*pqz

{ri.rys2d_inline_body(ni, nk, "c00x", "cp0x", "b00r", "b10r", "bp01r", "1.0_dp", "gx", 15)}
{ri.rys2d_inline_body(ni, nk, "c00y", "cp0y", "b00r", "b10r", "bp01r", "1.0_dp", "gy", 15)}
{ri.rys2d_inline_body(ni, nk, "c00z", "cp0z", "b00r", "b10r", "bp01r", "f00",    "gz", 15)}
{ri.hrr_inline_body(ni, nk, la, lb, lc, 0, "gx", ab_axes[0], "0.0_dp", "hx", 15)}
{ri.hrr_inline_body(ni, nk, la, lb, lc, 0, "gy", ab_axes[1], "0.0_dp", "hy", 15)}
{ri.hrr_inline_body(ni, nk, la, lb, lc, 0, "gz", ab_axes[2], "0.0_dp", "hz", 15)}

               idx = 0
               do cc = 1, {nc}
                  do bb = 1, {nb}
                     do aa = 1, {na}
                        idx = idx + 1
                        out(idx) = out(idx) &
                           + hx(cnx0_{tag}(aa), cnx1_{tag}(bb), cnx2_{tag}(cc), 0) &
                           * hy(cny0_{tag}(aa), cny1_{tag}(bb), cny2_{tag}(cc), 0) &
                           * hz(cnz0_{tag}(aa), cnz1_{tag}(bb), cnz2_{tag}(cc), 0)
                     end do
                  end do
               end do
            end do
"""

    body = f"""
   pure subroutine df_{tag}({sig})
      !! {'(ab|c)' if three else '(a|c)'} class ({la}{lb}|{lc}) -- Rys quadrature,
      !! {nroots} root{'s' if nroots != 1 else ''}, {'loop' if not unroll else 'unrolled'} form.
      !$acc routine seq
{decls}
      real(dp), intent(out) :: out({nout})

{decl_new}
      do oz = 1, {nout}
         out(oz) = 0.0_dp
      end do

{bra_loop_open}

{kp_head}
{kp_body}
         end do
      end do
   end subroutine df_{tag}
"""
    return body, tabs


MOD_HEADER = """!
! {title}
!
! GENERATED by scripts/gen_df.py --lmax {L} --lmax-aux {LA}{order_arg} --split -- do not edit.
!
! One module per class so no single compilation unit is enormous -- see
! gen_df.py's module docstring. The item routine and the launch stay in the
! same file as the class kernel, which is what the inlining rule requires;
! the dispatcher in trc_df_kernels.F90 is host code and may cross modules.
!
! TRC_CART_ORDER GUARD -- see gen_perclass.py's non-split header for the full
! explanation; same guard, one per class file.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
"""

BOYS_USE = ("   use trc_boys, only: dp, BOYS_MMAX, boys_table, BOYS_NCHEB, "
            "BOYS_NGRID, BOYS_TMAX, BOYS_DTINV, BOYS_DT\n")


def emit_2c_module(la, lc, L, LA, order, order_val, order_arg, rys, unroll_max):
    """trc_df_2c_k{la}{lc}: the (la|lc) class, its item routine and launch.

    `rys` picks the integral body: False (shipped ceiling, every class
    below SHIPPED_LA on both sides) keeps today's written-out VRR+HRR,
    byte-for-byte; True (above the ceiling, or `--force-rys`) uses
    `emit_class_rys` instead -- see that function and this module's
    docstring. Everything else about the module -- the launch, the item
    routine, the digestion -- is identical either way."""
    na, nc = len(cart(la)), len(cart(lc))
    mod = f"trc_df_2c_k{la}{lc}"
    drv = f"df_2c_k{la}{lc}"
    itm = f"df_2c_i{la}{lc}"
    kern = f"df_2c{la}{lc}"
    if rys:
        kernel_body, tabs = emit_class_rys(la, 0, lc, three=False, order=order,
                                           unroll_max=unroll_max)
        nroots = (la + lc)//2 + 1
        use_block = "   use trc_boys, only: dp\n" + ri.roots_table_use(nroots)
        boys_incl = ""
    else:
        kernel_body = emit_class(la, 0, lc, three=False, order=order)
        tabs = ""
        use_block = BOYS_USE
        boys_incl = '#include "inc/trc_boys_eval.inc"\n\n'

    text = MOD_HEADER.format(
        title=f"Class ({la}|{lc}) two-centre Coulomb kernel.",
        L=L, LA=LA, order_arg=order_arg, order_val=order_val)
    text += f"""module {mod}
{use_block}   implicit none
   private
   public :: {drv}

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp
{tabs}
contains

{boys_incl}   subroutine {drv}(nshell, maxnp, sh_np, sh_e, sh_c, sh_r, &
                    na, ord_a, nc, ord_c, nao, sh_ao, jmat)
      !! ({la}|{lc}) launch: `do concurrent` over the na x nc cross product of
      !! two shell-index LISTS the caller built by grouping on angular
      !! momentum (see trc_df_kernels.F90's df_2c_driver) -- not a runtime
      !! (la, lc) pair, so no select case runs inside the loop.
      integer,  intent(in)    :: nshell, maxnp, na, nc, nao
      integer,  intent(in)    :: sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      integer,  intent(in)    :: ord_a(na), ord_c(nc)
      real(dp), intent(inout) :: jmat(nao, nao)

      integer(kind=8) :: t

      do concurrent(t = 0:int(na, 8)*nc - 1)
         call {itm}(t, nshell, maxnp, sh_np, sh_e, sh_c, sh_r, &
                   na, ord_a, nc, ord_c, nao, sh_ao, jmat)
      end do
   end subroutine {drv}

   pure subroutine {itm}(t, nshell, maxnp, sh_np, sh_e, sh_c, sh_r, &
                         na, ord_a, nc, ord_c, nao, sh_ao, jmat)
      !$acc routine seq
      integer(kind=8), intent(in) :: t
      integer,  intent(in)    :: nshell, maxnp, na, nc, nao
      integer,  intent(in)    :: sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      integer,  intent(in)    :: ord_a(na), ord_c(nc)
      real(dp), intent(inout) :: jmat(nao, nao)

      integer  :: ii, jj, a, c, mu, nu
      real(dp) :: blk({na*nc})

      ii = ord_a(int(mod(t, int(na, 8))) + 1)
      jj = ord_c(int(t/na) + 1)
      call {kern}(sh_np(ii), sh_e(:, ii), sh_c(:, ii), sh_r(:, ii), &
                 sh_np(jj), sh_e(:, jj), sh_c(:, jj), sh_r(:, jj), blk)
      ! `blk` is addressed by the SAME affine formula the kernel's own
      ! digestion used to fill it (slot = 1 + ia + na*ic, emit_class's
      ! else-branch above) rather than a running counter incremented once
      ! per iteration: nvfortran 26.5 fails "Missing increment value" on
      ! `!$acc routine seq` translation of a nested loop that increments a
      ! counter and uses it as another array's subscript while ALSO writing
      ! a two-index array through loop-derived indices (mu, nu) -- confirmed
      ! with a four-line reproducer outside this project, bisected down to
      ! exactly that combination. An affine index from the loop variables
      ! sidesteps it and costs nothing: na is a compile-time literal here.
      do c = 0, {nc - 1}
         nu = sh_ao(jj) + c
         do a = 0, {na - 1}
            mu = sh_ao(ii) + a
            jmat(mu, nu) = blk(1 + a + {na}*c)
         end do
      end do
   end subroutine {itm}
{kernel_body}
end module {mod}
"""
    return mod, drv, text


def emit_3c_module(la, lb, lc, L, LA, order, order_val, order_arg, rys, unroll_max):
    """trc_df_3c_k{la}{lb}{lc}: the (la,lb|lc) class, item routine, launch.

    `rys` picks the integral body -- see `emit_2c_module`'s identical note.
    """
    na, nb, nc = len(cart(la)), len(cart(lb)), len(cart(lc))
    mod = f"trc_df_3c_k{la}{lb}{lc}"
    drv = f"df_3c_k{la}{lb}{lc}"
    itm = f"df_3c_i{la}{lb}{lc}"
    kern = f"df_3c{la}{lb}{lc}"
    if rys:
        kernel_body, tabs = emit_class_rys(la, lb, lc, three=True, order=order,
                                           unroll_max=unroll_max)
        nroots = (la + lb + lc)//2 + 1
        use_block = "   use trc_boys, only: dp\n" + ri.roots_table_use(nroots)
        boys_incl = ""
    else:
        kernel_body = emit_class(la, lb, lc, three=True, order=order)
        tabs = ""
        use_block = BOYS_USE
        boys_incl = '#include "inc/trc_boys_eval.inc"\n\n'

    text = MOD_HEADER.format(
        title=f"Class ({la}{lb}|{lc}) three-centre Coulomb kernel.",
        L=L, LA=LA, order_arg=order_arg, order_val=order_val)
    text += f"""module {mod}
{use_block}   implicit none
   private
   public :: {drv}

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp
{tabs}
contains

{boys_incl}   subroutine {drv}(nsh, nao, sh_ao, sh_r, &
                    npair, npp, p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                    npc, ord_p, &
                    nsx, nax, mnx, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                    nac, ord_a, &
                    tens)
      !! ({la}{lb}|{lc}) launch: `do concurrent` over the npc x nac cross
      !! product of a pair-index LIST and an auxiliary-shell-index LIST the
      !! caller built by grouping on angular momentum (see
      !! trc_df_kernels.F90's df_3c_driver) -- not a runtime (la, lb, lc)
      !! triple, so no select case runs inside the loop.
      integer,  intent(in)    :: nsh, nao, npair, npp, npc, nsx, nax, mnx, nac
      integer,  intent(in)    :: sh_ao(nsh)
      real(dp), intent(in)    :: sh_r(3, nsh)
      integer,  intent(in)    :: p_i(npair), p_j(npair)
      integer,  intent(in)    :: pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      integer,  intent(in)    :: ord_p(npc)
      integer,  intent(in)    :: ax_np(nsx), ax_ao(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: ord_a(nac)
      real(dp), intent(inout) :: tens(nao, nao, nax)

      integer(kind=8) :: t

      if (npc <= 0 .or. nac <= 0) return
      do concurrent(t = 0:int(npc, 8)*nac - 1)
         call {itm}(t, nsh, nao, sh_ao, sh_r, &
                   npair, npp, p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                   npc, ord_p, &
                   nsx, nax, mnx, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                   nac, ord_a, &
                   tens)
      end do
   end subroutine {drv}

   pure subroutine {itm}(t, nsh, nao, sh_ao, sh_r, &
                         npair, npp, p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                         npc, ord_p, &
                         nsx, nax, mnx, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                         nac, ord_a, &
                         tens)
      !$acc routine seq
      integer(kind=8), intent(in) :: t
      integer,  intent(in)    :: nsh, nao, npair, npp, npc, nsx, nax, mnx, nac
      integer,  intent(in)    :: sh_ao(nsh)
      real(dp), intent(in)    :: sh_r(3, nsh)
      integer,  intent(in)    :: p_i(npair), p_j(npair)
      integer,  intent(in)    :: pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      integer,  intent(in)    :: ord_p(npc)
      integer,  intent(in)    :: ax_np(nsx), ax_ao(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: ord_a(nac)
      real(dp), intent(inout) :: tens(nao, nao, nax)

      integer  :: ip, kk, ii, jj, a, c, e, mu, nu, pq, o
      real(dp) :: blk({na*nb*nc}), abx, aby, abz

      ip = ord_p(int(mod(t, int(npc, 8))) + 1)
      kk = ord_a(int(t/npc) + 1)
      ii = p_i(ip); jj = p_j(ip)
      abx = sh_r(1, ii) - sh_r(1, jj)
      aby = sh_r(2, ii) - sh_r(2, jj)
      abz = sh_r(3, ii) - sh_r(3, jj)
      o = pp_off(ip)

      call {kern}(pp_n(ip), pp_z(o + 1:o + pp_n(ip)), pp_k(o + 1:o + pp_n(ip)), &
                 pp_p(:, o + 1:o + pp_n(ip)), pp_a(:, o + 1:o + pp_n(ip)), &
                 abx, aby, abz, &
                 ax_np(kk), ax_e(:, kk), ax_c(:, kk), ax_r(:, kk), blk)

      ! (mu nu|P) is symmetric in mu <-> nu, so each thread fills both. On a
      ! diagonal shell pair the second write repeats the first with the same
      ! value -- same convention as the old single-dispatch df_3c_item.
      !
      ! `blk` is addressed by the SAME affine formula the kernel's own
      ! digestion used to fill it (slot = 1 + ia + na*(ib + nb*ic),
      ! emit_bra_hrr above) rather than a running counter -- see
      ! emit_2c_module's identical comment for the nvfortran bug this
      ! avoids. na and nb are compile-time literals here.
      do e = 0, {nc - 1}
         pq = ax_ao(kk) + e
         do c = 0, {nb - 1}
            nu = sh_ao(jj) + c
            do a = 0, {na - 1}
               mu = sh_ao(ii) + a
               tens(mu, nu, pq) = blk(1 + a + {na}*(c + {nb}*e))
               tens(nu, mu, pq) = blk(1 + a + {na}*(c + {nb}*e))
            end do
         end do
      end do
   end subroutine {itm}
{kernel_body}
end module {mod}
"""
    return mod, drv, text


def _parse_unroll_max(s):
    """Same parsing as gen_rys.py's --unroll-max-comp: an int, or 'inf'/
    'all' for no ceiling."""
    if s.strip().lower() in ("inf", "all"):
        return float("inf")
    return int(s)


def _parse_force_rys(s):
    """--force-rys: '' (default, no forcing), 'all' (every class), or a
    comma-separated set of raw class-digit tags -- see main()'s --force-rys
    help for the two/three-digit convention."""
    s = s.strip()
    if not s:
        return set()
    if s.lower() == "all":
        return "all"
    return {t.strip() for t in s.split(",") if t.strip()}


def emit_handoff(two, three, L, LA, order_arg, order_val):
    """trc_df_kernels.F90: DF_LMAX/DF_LMAX_AUX, df_2c_driver, df_3c_driver.

    Host code only -- grouping shells/pairs into index lists by angular
    momentum, and the two small dispatchers that pick a class module's
    driver. Free to cross module boundaries: nothing here runs inside a
    `do concurrent`, that is entirely within the per-class modules this
    imports.
    """
    use2 = "\n".join(f"   use trc_df_2c_k{la}{lc}, only: df_2c_k{la}{lc}"
                     for la, lc in two)
    use3 = "\n".join(f"   use trc_df_3c_k{la}{lb}{lc}, only: df_3c_k{la}{lb}{lc}"
                     for la, lb, lc in three)

    disp2 = ["""
   subroutine df_2c_dispatch(la, lc, nshell, maxnp, sh_np, sh_e, sh_c, sh_r, &
                             na, ord_a, nc, ord_c, nao, sh_ao, jmat)
      integer,  intent(in)    :: la, lc, nshell, maxnp, na, nc, nao
      integer,  intent(in)    :: sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      integer,  intent(in)    :: ord_a(na), ord_c(nc)
      real(dp), intent(inout) :: jmat(nao, nao)
      select case (la*(DF_LMAX_AUX + 1) + lc)"""]
    for la, lc in two:
        disp2.append(f"""      case ({la*(LA + 1) + lc}); call df_2c_k{la}{lc}(nshell, maxnp, sh_np, sh_e, sh_c, sh_r, &
                          na, ord_a, nc, ord_c, nao, sh_ao, jmat)""")
    disp2.append("""      end select
   end subroutine df_2c_dispatch
""")

    disp3 = ["""
   subroutine df_3c_dispatch(la, lb, lc, nsh, nao, sh_ao, sh_r, &
                             npair, npp, p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                             npc, ord_p, &
                             nsx, nax, mnx, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                             nac, ord_a, &
                             tens)
      integer,  intent(in)    :: la, lb, lc, nsh, nao, npair, npp, npc, nsx, nax, mnx, nac
      integer,  intent(in)    :: sh_ao(nsh)
      real(dp), intent(in)    :: sh_r(3, nsh)
      integer,  intent(in)    :: p_i(npair), p_j(npair)
      integer,  intent(in)    :: pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      integer,  intent(in)    :: ord_p(npc)
      integer,  intent(in)    :: ax_np(nsx), ax_ao(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: ord_a(nac)
      real(dp), intent(inout) :: tens(nao, nao, nax)
      select case ((la*(DF_LMAX + 1) + lb)*(DF_LMAX_AUX + 1) + lc)"""]
    for la, lb, lc in three:
        key = (la*(L + 1) + lb)*(LA + 1) + lc
        disp3.append(
            f"""      case ({key}); call df_3c_k{la}{lb}{lc}(nsh, nao, sh_ao, sh_r, &
                          npair, npp, p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                          npc, ord_p, nsx, nax, mnx, ax_np, ax_ao, ax_e, ax_c, ax_r, nac, ord_a, &
                          tens)""")
    disp3.append("""      end select
   end subroutine df_3c_dispatch
""")

    text = f"""!
! Density-fitting hand-off: DF_LMAX, DF_LMAX_AUX, and the two entry points
! trc_api.F90 calls.
!
! GENERATED by scripts/gen_df.py --lmax {L} --lmax-aux {LA}{order_arg} --split -- do not edit.
!
! The actual kernels live one class per module (trc_df_2c_k*.F90,
! trc_df_3c_k*.F90) -- see gen_df.py's module docstring for why. This file is
! HOST code only: `df_2c_driver`/`df_3c_driver` group shells (or shell pairs)
! into index lists by angular momentum, then a small `select case` dispatcher
! calls the one non-empty class module whose driver covers that group. None
! of that runs inside a `do concurrent`, so crossing a module boundary here
! costs nothing -- same reasoning as `pc_dispatch` in
! scripts/gen_perclass.py's --split output.
!
! No TRC_CART_ORDER guard here: this file only groups work and dispatches,
! and neither depends on Cartesian order -- the guard lives in each class
! module, where the order-dependent digestion actually is.
!
module trc_df_kernels
   use trc_boys, only: dp
{use2}
{use3}
   implicit none
   private
   public :: df_2c_driver, df_3c_driver, DF_LMAX, DF_LMAX_AUX

   integer, parameter :: DF_LMAX = {L}, DF_LMAX_AUX = {LA}

contains

   subroutine df_2c_driver(nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                           sh_e, sh_c, sh_r, jmat)
      integer,  intent(in)  :: nshell, nao, maxnp
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      real(dp), intent(out) :: jmat(nao, nao)

      integer :: i, l, la, lc
      integer :: cnt(0:DF_LMAX_AUX), off(0:DF_LMAX_AUX + 1), cur(0:DF_LMAX_AUX)
      integer, allocatable :: ord(:)

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(l=1:nao, i=1:nao)
         jmat(i, l) = 0.0_dp
      end do

      ! Counting sort of the shells by angular momentum, on the host: every
      ! class kernel below takes a fixed (la, lc) and a LIST, not a runtime
      ! l, so the select case in df_2c_dispatch runs once per non-empty
      ! class pair rather than once per shell pair.
      allocate (ord(nshell))
      cnt = 0
      do i = 1, nshell
         cnt(sh_l(i)) = cnt(sh_l(i)) + 1
      end do
      off(0) = 0
      do l = 0, DF_LMAX_AUX
         off(l + 1) = off(l) + cnt(l)
      end do
      cur = off(0:DF_LMAX_AUX)
      do i = 1, nshell
         l = sh_l(i)
         cur(l) = cur(l) + 1
         ord(cur(l)) = i
      end do

      do la = 0, DF_LMAX_AUX
         if (off(la + 1) <= off(la)) cycle
         do lc = 0, DF_LMAX_AUX
            if (off(lc + 1) <= off(lc)) cycle
            call df_2c_dispatch(la, lc, nshell, maxnp, sh_np, sh_e, sh_c, sh_r, &
                                off(la + 1) - off(la), ord(off(la) + 1:off(la + 1)), &
                                off(lc + 1) - off(lc), ord(off(lc) + 1:off(lc + 1)), &
                                nao, sh_ao, jmat)
         end do
      end do
      deallocate (ord)
   end subroutine df_2c_driver

   subroutine df_3c_driver(nsh, nao, mnp, sh_l, sh_ao, sh_r, &
                           npair, npp, p_i, p_j, pp_off, pp_n, &
                           pp_z, pp_k, pp_p, pp_a, &
                           nsx, nax, mnx, ax_l, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                           tens)
      integer,  intent(in)  :: nsh, nao, mnp, npair, npp, nsx, nax, mnx
      integer,  intent(in)  :: sh_l(nsh), sh_ao(nsh)
      real(dp), intent(in)  :: sh_r(3, nsh)
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      integer,  intent(in)  :: pp_off(npair), pp_n(npair)
      real(dp), intent(in)  :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      integer,  intent(in)  :: ax_l(nsx), ax_np(nsx), ax_ao(nsx)
      real(dp), intent(in)  :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      real(dp), intent(out) :: tens(nao, nao, nax)

      integer :: i, j, k, l, la, lb, lc
      integer :: pcnt(0:(DF_LMAX + 1)*(DF_LMAX + 1) - 1)
      integer :: poff(0:(DF_LMAX + 1)*(DF_LMAX + 1))
      integer :: pcur(0:(DF_LMAX + 1)*(DF_LMAX + 1) - 1)
      integer, allocatable :: ordp(:)
      integer :: acnt(0:DF_LMAX_AUX), aoff(0:DF_LMAX_AUX + 1), acur(0:DF_LMAX_AUX)
      integer, allocatable :: orda(:)

      ! `do concurrent` rather than an OpenACC parallel loop -- same note as
      ! df_2c_driver above.
      do concurrent(k=1:nax, j=1:nao, i=1:nao)
         tens(i, j, k) = 0.0_dp
      end do

      ! Counting sort of the pairs by (la, lb), and of the auxiliary shells
      ! by lc, on the host -- see df_2c_driver's comment for why.
      allocate (ordp(npair))
      pcnt = 0
      do i = 1, npair
         l = sh_l(p_i(i))*(DF_LMAX + 1) + sh_l(p_j(i))
         pcnt(l) = pcnt(l) + 1
      end do
      poff(0) = 0
      do l = 0, (DF_LMAX + 1)*(DF_LMAX + 1) - 1
         poff(l + 1) = poff(l) + pcnt(l)
      end do
      pcur = poff(0:(DF_LMAX + 1)*(DF_LMAX + 1) - 1)
      do i = 1, npair
         l = sh_l(p_i(i))*(DF_LMAX + 1) + sh_l(p_j(i))
         pcur(l) = pcur(l) + 1
         ordp(pcur(l)) = i
      end do

      allocate (orda(nsx))
      acnt = 0
      do i = 1, nsx
         acnt(ax_l(i)) = acnt(ax_l(i)) + 1
      end do
      aoff(0) = 0
      do l = 0, DF_LMAX_AUX
         aoff(l + 1) = aoff(l) + acnt(l)
      end do
      acur = aoff(0:DF_LMAX_AUX)
      do i = 1, nsx
         l = ax_l(i)
         acur(l) = acur(l) + 1
         orda(acur(l)) = i
      end do

      do la = 0, DF_LMAX
         do lb = 0, DF_LMAX
            k = la*(DF_LMAX + 1) + lb
            if (poff(k + 1) <= poff(k)) cycle
            do lc = 0, DF_LMAX_AUX
               if (aoff(lc + 1) <= aoff(lc)) cycle
               call df_3c_dispatch(la, lb, lc, nsh, nao, sh_ao, sh_r, &
                                   npair, npp, p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                                   poff(k + 1) - poff(k), ordp(poff(k) + 1:poff(k + 1)), &
                                   nsx, nax, mnx, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                                   aoff(lc + 1) - aoff(lc), orda(aoff(lc) + 1:aoff(lc + 1)), &
                                   tens)
            end do
         end do
      end do
      deallocate (ordp, orda)
   end subroutine df_3c_driver
{"\n".join(disp2)}
{"\n".join(disp3)}
end module trc_df_kernels
"""
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=2,
                    help="max l on the orbital (bra) shells")
    # Default one above the orbital ceiling: a fitting basis normally
    # outranks the orbital basis by one, and the shipped set (--lmax 2) was
    # built with f on the auxiliary shell. CI and CMake regenerate with
    # --lmax alone, so the default IS the shipped configuration; a default
    # of 2 here once made CI's regeneration disagree with the committed
    # kernels by nine thousand lines while nothing was actually wrong.
    ap.add_argument("--lmax-aux", type=int, default=None,
                    help="max l on the auxiliary shell (default: lmax + 1)")
    ap.add_argument("--order", choices=("lexicographic", "gamess"),
                    default="lexicographic",
                    help="Cartesian order of the AO axes of `out`, the "
                         "(a|c)/(ab|c) tensor. INTERNAL VRR indices "
                         "(gv.cart_list, `I`) stay lexicographic regardless.")
    ap.add_argument("--split", metavar="DIR", required=True,
                    help="emit one module per class into DIR, plus the "
                         "trc_df_kernels.F90 hand-off module. The only mode: "
                         "see gen_df.py's module docstring for why a single "
                         "file does not scale past (dd|f).")
    ap.add_argument("--unroll-max-comp", default="0",
                    help="unroll a Rys class's root/component loops (see "
                         "rys_inline.py's UNROLLED MODE) when its na*nb*nc "
                         "is at most this; 'inf'/'all' unrolls every Rys "
                         "class. Default 0: every Rys class keeps the "
                         "array/loop form. The SAME threshold gen_rys.py's "
                         "--unroll-max-comp takes (CMake's "
                         "TERCO_RYS_UNROLL_MAX) -- pass that value here too, "
                         "do not pick a different one.")
    ap.add_argument("--force-rys", default="",
                    help="TESTING ONLY, never passed by the shipped "
                         "configuration: force the Rys body for classes "
                         "that would otherwise keep the written-out "
                         "VRR+HRR. 'all' forces every class; otherwise a "
                         "comma-separated list of raw class-digit tags -- "
                         "two digits 'la lc' for a two-centre class (e.g. "
                         "33), three digits 'la lb lc' for a three-centre "
                         "one (e.g. 223) -- unambiguous since every "
                         "angular momentum here is one decimal digit.")
    args = ap.parse_args()
    L = args.lmax
    LA = args.lmax_aux if args.lmax_aux is not None else L + 1
    order = args.order
    order_val = 0 if order == "lexicographic" else 1
    order_arg = "" if order == "lexicographic" else f" --order {order}"

    os.makedirs(args.split, exist_ok=True)

    unroll_max = _parse_unroll_max(args.unroll_max_comp)
    force = _parse_force_rys(args.force_rys)

    def use_rys(digits, la, lb, lc, three):
        if force == "all" or digits in force:
            return True
        return default_use_rys(la, lb, lc, three)

    two, three = [], []
    n_rys2 = n_rys3 = 0
    # Both sides of (a|c) are auxiliary shells: the metric (P|Q) needs the
    # bra up to LA as well, and a missing (f|f) case left the block unset and
    # the metric indefinite the first time an RI set with f functions was
    # tried.
    for la in range(LA + 1):
        for lc in range(LA + 1):
            rys = use_rys(f"{la}{lc}", la, 0, lc, False)
            n_rys2 += rys
            mod, drv, text = emit_2c_module(la, lc, L, LA, order, order_val,
                                            order_arg, rys, unroll_max)
            with open(os.path.join(args.split, mod + ".F90"), "w") as fh:
                fh.write(text)
            two.append((la, lc))
    for la in range(L + 1):
        for lb in range(L + 1):
            for lc in range(LA + 1):
                rys = use_rys(f"{la}{lb}{lc}", la, lb, lc, True)
                n_rys3 += rys
                mod, drv, text = emit_3c_module(la, lb, lc, L, LA, order,
                                                order_val, order_arg, rys,
                                                unroll_max)
                with open(os.path.join(args.split, mod + ".F90"), "w") as fh:
                    fh.write(text)
                three.append((la, lb, lc))

    handoff = emit_handoff(two, three, L, LA, order_arg, order_val)
    with open(os.path.join(args.split, "trc_df_kernels.F90"), "w") as fh:
        fh.write(handoff)

    total_lines = handoff.count("\n")
    print(f"wrote {len(two)} two-centre ({n_rys2} Rys) + {len(three)} "
          f"three-centre ({n_rys3} Rys) class modules + trc_df_kernels.F90 "
          f"into {args.split}")


if __name__ == "__main__":
    main()
