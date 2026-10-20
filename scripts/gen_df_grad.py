#!/usr/bin/env python3
"""
Emit the nuclear-DERIVATIVE kernels of the density-fitting integrals:

    d(ab|c)/dR   three-centre, a and b orbital shells, c an auxiliary shell
    d(a|c)/dR    two-centre, both auxiliary shells

already CONTRACTED against a coefficient block -- Y^c_ab for the three-centre
integral, Z_ac for the two-centre one -- so each kernel call returns a few
gradient numbers rather than a derivative-integral block. These are what the
fitted (`--df`) SCF gradient needs (src/scf/trc_scf_grad_df.F90 builds Y and Z;
this file only differentiates the integrals).

WHY RYS, FOR EVERY CLASS
-------------------------
A derivative raises one index by one, so differentiating a (dd|g) integral
needs an (fd|g)-sized table. That is past the written-out VRR+HRR's compile
wall at every class with a g auxiliary shell (see gen_df.py's own docstring:
(fd|g) took 22 minutes). The four-centre gradient is Rys-only for the same
reason (gen_rys_grad.py). So every class here is a Rys quadrature, reusing
scripts/rys_inline.py's inlined roots lookup, 2D recurrence and horizontal
transfer -- the same text gen_df.py's above-ceiling classes and the
four-centre Rys kernels use.

THE ALGEBRA
-----------
The derivative of a Cartesian primitive x_A^i exp(-alpha x_A^2) with respect
to A_x is 2 alpha x_A^(i+1) e - i x_A^(i-1) e. Per axis, with h(i, j, k) the
horizontal-transfer table of the 1D Rys integrals:

    fi = 2 alpha h(i+1, j, k) - i h(i-1, j, k)                      d/dA
    fj = 2 beta (h(i+1, j, k) + (A-B) h(i, j, k)) - j h(i, j-1, k)  d/dB
                                     (h(i, j+1) = h(i+1, j) + (A-B) h(i, j))

and the three-dimensional derivative integral is fi_x * h_y * h_z (and its y,
z twins). The table is built to (la+1, lb, lc) -- only the bra-combined index
is raised, never lb and never lc -- so the 2D recurrence runs to
nig = la+lb+1 on the bra and nkg = lc on the ket. The auxiliary centre's
derivative is never computed: the integral is translation-invariant, so

    three-centre:  dC = -(dA + dB)        two-centre:  dC = -dA

which is also why the ket needs no raise at all. No Hellmann-Feynman term
exists (no nucleus appears in a Coulomb integral between basis functions).

nroots = (la + lb + lc + 1)//2 + 1: the polynomial degree actually summed is
la+lb+lc+1 (one index raised by one), which the existing root tables cover
up to (dd|g) (5 roots) and (g|g) (5 roots).

THE CONTRACTION. Each kernel takes the coefficient block y(na*nb*nc) (or
z(na*nc)), in the SAME slot order gen_df.py's energy kernels write their
integrals in (slot = 1 + ia + na*(ib + nb*ic), cartorder.py's order), and
folds it into six running sums (A x,y,z and B x,y,z) inside the root loop
-- the four-centre gradient's own "fold, do not materialise" pattern
(gen_rys_grad.py), with no derivative block ever held.

THE PAIR PLAN, NOT THE PAIR LIST. The fold needs alpha and beta
separately; trc_pairlist_t only keeps their sum, and its layout stays as it
is (trc_scf_df.F90 records an nvfortran miscompile from changing it). So
the hand-off's df_grad_3c_plan builds the kernels' own view of the pair
list, once per gradient, on the host: every pair turned so its first shell
has the higher l (CANONICAL CLASSES, below), grouped by class, sorted by
primitive-pair count inside a class -- one warp, one primitive trip count;
the pairs came interleaved 1x1 / 9x9 before, and the warp ran the longest
-- and each pair's primitive pairs (zeta, K_ab, P, alpha, beta) computed
once into pd(8, npp) instead of once per (pair, auxiliary shell) item. The
auxiliary shells are sorted by primitive count the same way. Measured on a
V100, gly10/cc-pVDZ, with the J part below: 0.49 -> 0.25 s.

CANONICAL CLASSES. Only (la >= lb | lc) and (la >= lc) are generated; a pair
the other way round is handed in turned (y transposed to match, which Y's
and Z's symmetry make exact). 30 three-centre + 15 two-centre modules
instead of 45 + 25, and the horizontal transfer always runs towards the
lower l.

THE J PART IN THE KERNEL. Y = j c_P D_mn + Y_K: the three-centre item forms
the J part from c and D itself (`jc`, `cq`, `dens`), so a J-only gradient
(a pure functional) never stores an (n, n, naux) tensor at all and runs
one pass over the whole auxiliary basis, and with K the dense input is
Y_K alone (`ky` = 1; the GEMM writes it with beta = 0, no fill pass) -- or,
when trc_scf_grad_df.F90 has to chunk the occupied orbitals and so cannot
rebuild X's rows j > i, as M with Y_K = M + M^T (`ky` = 2), the kernel
reading both.

HOST-CALLABLE CLASS KERNELS. Every class kernel is public, and the hand-off
module (trc_df_grad_kernels.F90) has a host `select case` over them
(df_grad_3c_block / df_grad_2c_block), so a test can call one class kernel
with a one-hot coefficient block and read back the raw derivative integrals
(test/check_df_grad_ints.F90) -- the exact code the drivers run.

UNROLLING. --unroll-max-comp unrolls a class's components, and
--unroll-roots-max-comp (default: the same number) its root loop too. The
shipped set (CMake's TERCO_DFG_UNROLL_MAX) unrolls the components up to 135
and never the roots: a run-time root loop around unrolled components ran
as fast as the fully unrolled form on a V100 (within 3% per class) in a
third of the text, and 2-3x faster than the loop form (whose run-time
indexed tables live in local memory) for the classes above 81.

USAGE
    gen_df_grad.py --lmax 2 --lmax-aux 4 --unroll-max-comp 135 --unroll-roots-max-comp 0 --split src/generated
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

PD_W = 8
# Doubles per primitive pair in the three-centre plan's pd(PD_W, npp):
# zeta, K_ab, P(3), alpha, beta, and one pad -- 64 bytes, two whole sectors.
ri = _load("ri", os.path.join(HERE, "rys_inline.py"))


def cart(l, order="lexicographic"):
    return cartorder.components(l, order)


def ncart(l):
    return (l + 1)*(l + 2)//2


def _parse_unroll_max(s):
    if s.strip().lower() in ("inf", "all"):
        return float("inf")
    return int(s)


def _fold_unrolled(la, lb, lc, three, order, pad):
    """The (aa, bb, cc) contraction for one root, unrolled Python-side:
    every table read is a literal-indexed named scalar, and a `- e*h(e-1)`
    term whose literal exponent is 0 is simply not emitted."""
    comps = [cart(l, order) for l in (la, lb, lc)]
    na, nb, nc = len(comps[0]), len(comps[1]), len(comps[2])
    out = []
    for ic in range(nc):
        for ib in range(nb):
            for ia in range(na):
                slot = 1 + ia + na*(ib + nb*ic)
                e = (comps[0][ia], comps[1][ib], comps[2][ic])
                lines = [f"{pad}yv = y({slot})"]
                v, fi, fj = {}, {}, {}
                for ax, xn in enumerate("xyz"):
                    ei, ej, ek = e[0][ax], e[1][ax], e[2][ax]
                    h = lambda i, j, k, xn=xn: ri.nm(f"h{xn}", i, j, k, 0)
                    v[xn] = h(ei, ej, ek)
                    s = f"2.0_dp*al*{h(ei + 1, ej, ek)}"
                    if ei > 0:
                        s += f" - {float(ei)}_dp*{h(ei - 1, ej, ek)}"
                    fi[xn] = s
                    if three:
                        s = f"2.0_dp*be*({h(ei + 1, ej, ek)} + ab{xn}*{h(ei, ej, ek)})"
                        if ej > 0:
                            s += f" - {float(ej)}_dp*{h(ei, ej - 1, ek)}"
                        fj[xn] = s
                lines.append(f"{pad}pyz = yv*{v['y']}*{v['z']}; pxz = yv*{v['x']}*{v['z']}; "
                             f"pxy = yv*{v['x']}*{v['y']}")
                lines.append(f"{pad}gax = gax + ({fi['x']})*pyz")
                lines.append(f"{pad}gay = gay + ({fi['y']})*pxz")
                lines.append(f"{pad}gaz = gaz + ({fi['z']})*pxy")
                if three:
                    lines.append(f"{pad}gbx = gbx + ({fj['x']})*pyz")
                    lines.append(f"{pad}gby = gby + ({fj['y']})*pxz")
                    lines.append(f"{pad}gbz = gbz + ({fj['z']})*pxy")
                out.append("\n".join(lines))
    return "\n".join(out) + "\n"


def _fold_loop(la, lb, lc, three, tag, pad):
    """The loop-form contraction, for classes above --unroll-max-comp: the
    component exponents come from the cn{x,y,z}{0,1,2}_{tag} tables and the
    `- e*h(e-1)` terms are guarded at run time. `yv` is addressed by the
    same affine slot formula the energy kernels' digestion uses -- never a
    running counter (see gen_df.py's emit_2c_module for the nvfortran bug
    a counter-plus-two-index-write trips)."""
    na, nb, nc = ncart(la), ncart(lb), ncart(lc)
    p = pad
    s = [f"{p}do ce = 1, {nc}",
         f"{p}   do bb = 1, {nb}",
         f"{p}      do aa = 1, {na}",
         f"{p}         yv = y(aa + {na}*(bb - 1 + {nb}*(ce - 1)))"]
    q = p + "         "
    for xn in "xyz":
        s.append(f"{q}i{xn} = cn{xn}0_{tag}(aa); j{xn} = cn{xn}1_{tag}(bb); k{xn} = cn{xn}2_{tag}(ce)")
    for xn in "xyz":
        s.append(f"{q}v{xn} = h{xn}(i{xn}, j{xn}, k{xn}, 0)")
        s.append(f"{q}fi{xn} = 2.0_dp*al*h{xn}(i{xn} + 1, j{xn}, k{xn}, 0)")
        s.append(f"{q}if (i{xn} > 0) fi{xn} = fi{xn} - real(i{xn}, dp)*h{xn}(i{xn} - 1, j{xn}, k{xn}, 0)")
        if three:
            s.append(f"{q}fj{xn} = 2.0_dp*be*(h{xn}(i{xn} + 1, j{xn}, k{xn}, 0) + ab{xn}*v{xn})")
            s.append(f"{q}if (j{xn} > 0) fj{xn} = fj{xn} - real(j{xn}, dp)*h{xn}(i{xn}, j{xn} - 1, k{xn}, 0)")
    s.append(f"{q}pyz = yv*vy*vz; pxz = yv*vx*vz; pxy = yv*vx*vy")
    s.append(f"{q}gax = gax + fix*pyz; gay = gay + fiy*pxz; gaz = gaz + fiz*pxy")
    if three:
        s.append(f"{q}gbx = gbx + fjx*pyz; gby = gby + fjy*pxz; gbz = gbz + fjz*pxy")
    s += [f"{p}      end do", f"{p}   end do", f"{p}end do"]
    return "\n".join(s) + "\n"


def emit_kernel(la, lb, lc, three, order, unroll_max, roots_max):
    """One class's contracted derivative kernel. Returns (text, tables).

    Three forms, by the class's component count nout = na*nb*nc: nout <=
    roots_max (and <= unroll_max) unrolls the roots and the components;
    nout <= unroll_max unrolls the components inside a run-time root loop
    (one copy of the 2D/HRR/fold text instead of nroots); above that, the
    loop form (run-time-indexed tables, which live in local memory)."""
    nig, nkg = la + lb + 1, lc
    nroots = (la + lb + lc + 1)//2 + 1
    na, nb, nc = ncart(la), ncart(lb), ncart(lc)
    nout = na*nb*nc
    unroll = nout <= unroll_max
    rloop = unroll and nout > roots_max
    tag = f"3c{la}{lb}{lc}" if three else f"2c{la}{lc}"
    form = "root-loop unrolled" if rloop else "unrolled" if unroll else "loop"

    if three:
        sig = "npp, poff, pd, ra, rb, npc, ec, cc_, rc, y, ga, gb"
        decls = (
            "      integer,  intent(in)  :: npp, poff, npc\n"
            "      real(dp), intent(in)  :: pd(PD_W, *), ra(3), rb(3)\n"
            "      real(dp), intent(in)  :: ec(npc), cc_(npc), rc(3)\n"
            f"      real(dp), intent(in)  :: y({nout})\n"
            "      real(dp), intent(out) :: ga(3), gb(3)\n")
        doc = (f"      !! d(ab|c)/dA and d(ab|c)/dB of class ({la}{lb}|{lc}), contracted with\n"
               f"      !! y(slot): ga = sum y d/dA, gb = sum y d/dB (dC = -(ga+gb)). Rys,\n"
               f"      !! {nroots} root{'s' if nroots != 1 else ''}, {form} form.\n"
               f"      !! The bra's npp primitive pairs are columns poff+1..poff+npp of pd\n"
               f"      !! (zeta, K_ab, P, alpha, beta: trc_df_grad_kernels' df_grad_3c_plan).\n")
        bra_open = """\
      abx = ra(1) - rb(1); aby = ra(2) - rb(2); abz = ra(3) - rb(3)
      do ip = poff + 1, poff + npp
         zeta = pd(1, ip)
         kab = pd(2, ip)
         px = pd(3, ip); py = pd(4, ip); pz = pd(5, ip)
         al = pd(6, ip)
         be = pd(7, ip)
         pax = px - ra(1); pay = py - ra(2); paz = pz - ra(3)
"""
        bra_close = "      end do\n"
        ab_axes = ("abx", "aby", "abz")
        ind = 12
    else:
        sig = "npa, ea, ca, ra, npc, ec, cc_, rc, y, ga"
        decls = (
            "      integer,  intent(in)  :: npa, npc\n"
            "      real(dp), intent(in)  :: ea(npa), ca(npa), ra(3)\n"
            "      real(dp), intent(in)  :: ec(npc), cc_(npc), rc(3)\n"
            f"      real(dp), intent(in)  :: y({nout})\n"
            "      real(dp), intent(out) :: ga(3)\n")
        doc = (f"      !! d(a|c)/dA of class ({la}|{lc}), contracted with y(slot) (dC = -ga).\n"
               f"      !! Rys, {nroots} root{'s' if nroots != 1 else ''}, {form} form.\n")
        bra_open = """\
      do ip = 1, npa
         al = ea(ip)
         zeta = al
         kab = ca(ip)
         px = ra(1); py = ra(2); pz = ra(3)
         pax = 0.0_dp; pay = 0.0_dp; paz = 0.0_dp
"""
        bra_close = "      end do\n"
        ab_axes = ("0.0_dp", "0.0_dp", "0.0_dp")
        ind = 12

    pad = " "*ind
    loc = ("      integer  :: ip, kp\n"
           "      real(dp) :: al, be, zeta, eta, zpe, rho, kab, pref\n"
           "      real(dp) :: px, py, pz, pqx, pqy, pqz, tval\n"
           "      real(dp) :: pax, pay, paz, qcx, qcy, qcz, abx, aby, abz\n"
           f"      real(dp) :: ru({nroots}), rw({nroots})\n"
           "      real(dp) :: u2, dinv, b00r, b10r, bp01r, f00\n"
           "      real(dp) :: c00x, c00y, c00z, cp0x, cp0y, cp0z\n"
           "      real(dp) :: gax, gay, gaz, gbx, gby, gbz, yv, pyz, pxz, pxy\n"
           + ri.roots_decls())
    tabs = ""
    if unroll:
        loc += ri.unrolled_table_decls(nig, nkg, la + 1, lb, lc, 0, 6)
        if rloop:
            loc += "      integer  :: rr\n"
    else:
        loc += (f"      real(dp) :: gx(0:{nig}, 0:{nkg}), gy(0:{nig}, 0:{nkg}), gz(0:{nig}, 0:{nkg})\n"
                f"      real(dp) :: hx(0:{la + 1}, 0:{lb}, 0:{lc}, 0:0)\n"
                f"      real(dp) :: hy(0:{la + 1}, 0:{lb}, 0:{lc}, 0:0)\n"
                f"      real(dp) :: hz(0:{la + 1}, 0:{lb}, 0:{lc}, 0:0)\n"
                "      integer  :: rr, aa, bb, ce, ix, iy, iz, jx, jy, jz, kx, ky, kz\n"
                "      real(dp) :: vx, vy, vz, fix, fiy, fiz, fjx, fjy, fjz\n"
                + ri.rys2d_hrr_decls(la + 1, lb, 0, nkg))
        for cnm, ax in (("cnx", 0), ("cny", 1), ("cnz", 2)):
            for which, l in enumerate((la, lb, lc)):
                comps = cart(l, order)
                vals = ", ".join(str(c[ax]) for c in comps)
                tabs += f"   integer, parameter :: {cnm}{which}_{tag}({len(comps)}) = [{vals}]\n"

    kp_head = f"""\
{pad[:-3]}do kp = 1, npc
{pad}eta = ec(kp)
{pad}zpe = zeta + eta
{pad}rho = zeta*eta/zpe
{pad}! The ket is one auxiliary shell: Q = C, QC = 0, K_cd = 1.
{pad}qcx = 0.0_dp; qcy = 0.0_dp; qcz = 0.0_dp
{pad}pqx = px - rc(1); pqy = py - rc(2); pqz = pz - rc(3)
{pad}tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)
{pad}pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*kab*cc_(kp)
"""

    def root_setup(r, p):
        return f"""\
{p}u2 = ru({r})/(1.0_dp - ru({r}))*rho
{p}dinv = 1.0_dp/(zeta*eta + u2*(zeta + eta))
{p}b00r  = 0.5_dp*u2*dinv
{p}b10r  = 0.5_dp*(eta + u2)*dinv
{p}bp01r = 0.5_dp*(zeta + u2)*dinv
{p}f00 = pref*rw({r})
{p}c00x = pax - u2*eta*dinv*pqx
{p}c00y = pay - u2*eta*dinv*pqy
{p}c00z = paz - u2*eta*dinv*pqz
{p}cp0x = qcx + u2*zeta*dinv*pqx
{p}cp0y = qcy + u2*zeta*dinv*pqy
{p}cp0z = qcz + u2*zeta*dinv*pqz
"""

    if rloop:
        p2 = pad + "   "
        i2 = ind + 3
        kp_body = "".join(
            [ri.roots_inline_body(nroots, "tval", ind),
             f"{pad}do rr = 1, {nroots}\n",
             root_setup("rr", p2),
             ri.rys2d_scalar_assignments(nig, nkg, "c00x", "cp0x", "b00r", "b10r", "bp01r", "1.0_dp", "gx", i2),
             ri.rys2d_scalar_assignments(nig, nkg, "c00y", "cp0y", "b00r", "b10r", "bp01r", "1.0_dp", "gy", i2),
             ri.rys2d_scalar_assignments(nig, nkg, "c00z", "cp0z", "b00r", "b10r", "bp01r", "f00", "gz", i2)]
            + [ri.hrr_scalar_assignments(nig, nkg, la + 1, lb, lc, 0, f"g{xn}", abn, "0.0_dp", f"h{xn}", i2)
               for xn, abn in zip("xyz", ab_axes)]
            + [_fold_unrolled(la, lb, lc, three, order, p2), f"{pad}end do\n"])
    elif unroll:
        mid = [ri.roots_inline_body(nroots, "tval", ind)]
        for r in range(1, nroots + 1):
            mid.append(root_setup(r, pad))
            mid.append(ri.rys2d_scalar_assignments(nig, nkg, "c00x", "cp0x", "b00r", "b10r", "bp01r",
                                                   "1.0_dp", "gx", ind))
            mid.append(ri.rys2d_scalar_assignments(nig, nkg, "c00y", "cp0y", "b00r", "b10r", "bp01r",
                                                   "1.0_dp", "gy", ind))
            mid.append(ri.rys2d_scalar_assignments(nig, nkg, "c00z", "cp0z", "b00r", "b10r", "bp01r",
                                                   "f00", "gz", ind))
            for xn, abn in zip("xyz", ab_axes):
                mid.append(ri.hrr_scalar_assignments(nig, nkg, la + 1, lb, lc, 0, f"g{xn}", abn, "0.0_dp",
                                                     f"h{xn}", ind))
            mid.append(_fold_unrolled(la, lb, lc, three, order, pad))
        kp_body = "".join(mid)
    else:
        p2 = pad + "   "
        kp_body = (ri.roots_inline_body(nroots, "tval", ind)
                   + f"{pad}do rr = 1, {nroots}\n"
                   + root_setup("rr", p2)
                   + ri.rys2d_inline_body(nig, nkg, "c00x", "cp0x", "b00r", "b10r", "bp01r", "1.0_dp", "gx",
                                          ind + 3)
                   + ri.rys2d_inline_body(nig, nkg, "c00y", "cp0y", "b00r", "b10r", "bp01r", "1.0_dp", "gy",
                                          ind + 3)
                   + ri.rys2d_inline_body(nig, nkg, "c00z", "cp0z", "b00r", "b10r", "bp01r", "f00", "gz",
                                          ind + 3)
                   + "".join(ri.hrr_inline_body(nig, nkg, la + 1, lb, lc, 0, f"g{xn}", abn, "0.0_dp",
                                                f"h{xn}", ind + 3)
                             for xn, abn in zip("xyz", ab_axes))
                   + _fold_loop(la, lb, lc, three, tag, p2)
                   + f"{pad}end do\n")

    tail = "      ga(1) = gax; ga(2) = gay; ga(3) = gaz\n"
    if three:
        tail += "      gb(1) = gbx; gb(2) = gby; gb(3) = gbz\n"

    text = f"""
   pure subroutine dfg_{tag}({sig})
{doc}      !$acc routine seq
{decls}
{loc}
      gax = 0.0_dp; gay = 0.0_dp; gaz = 0.0_dp
      gbx = 0.0_dp; gby = 0.0_dp; gbz = 0.0_dp
      abx = 0.0_dp; aby = 0.0_dp; abz = 0.0_dp
      be = 0.0_dp

{bra_open}
{kp_head}
{kp_body}{pad[:-3]}end do
{bra_close}{tail}   end subroutine dfg_{tag}
"""
    return text, tabs, nroots


MOD_HEADER = """!
! {title}
!
! GENERATED by scripts/gen_df_grad.py --lmax {L} --lmax-aux {LA}{order_arg} --split -- do not edit.
!
! One module per class, for the same compile-time reason as gen_df.py's own
! split; the item routine and the launch stay in the class kernel's file
! (stdpar inlines only within a translation unit). The hand-off module
! trc_df_grad_kernels.F90 is host code and may cross modules.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
"""


def emit_3c_module(la, lb, lc, L, LA, order, order_val, order_arg, unroll_max, roots_max):
    na, nb, nc = ncart(la), ncart(lb), ncart(lc)
    tag = f"{la}{lb}{lc}"
    mod, drv, itm, kern = f"trc_dfg_3c_k{tag}", f"dfg_3c_k{tag}", f"dfg_3c_i{tag}", f"dfg_3c{tag}"
    body, tabs, nroots = emit_kernel(la, lb, lc, True, order, unroll_max, roots_max)
    text = MOD_HEADER.format(title=f"Class ({la}{lb}|{lc}) three-centre Coulomb derivative, contracted.",
                             L=L, LA=LA, order_arg=order_arg, order_val=order_val)
    text += f"""module {mod}
   use pic_types, only: dp, int64
{ri.roots_table_use(nroots)}   implicit none
   private
   public :: {drv}, {kern}

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp
   integer, parameter :: PD_W = {PD_W}
{tabs}
contains

   subroutine {drv}(nsh, nao, sh_ao, sh_at, sh_r, npq, q_a, q_b, q_o, q_n, npp, pd, &
                    nsx, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, nac, ord_a, &
                    natm, nlane, ky, ny1, ny3, ytens, jc, nqc, qoff, cq, dens, gp)
      !! ({la}{lb}|{lc}) launch: `do concurrent` over the npq x nac cross product of
      !! this class's shell pairs (q_a/q_b, la >= lb, sorted by primitive
      !! count: trc_df_grad_kernels.F90's df_grad_3c_plan) and an auxiliary-
      !! shell index LIST (df_grad_3c_driver).
      integer,  intent(in)    :: nsh, nao, npq, npp, nsx, mnx, nac, natm, nlane, ny1, ny3, nqc, qoff
      integer,  intent(in)    :: sh_ao(nsh), sh_at(nsh)
      real(dp), intent(in)    :: sh_r(3, nsh)
      integer,  intent(in)    :: q_a(npq), q_b(npq), q_o(npq), q_n(npq)
      real(dp), intent(in)    :: pd(PD_W, npp)
      integer,  intent(in)    :: ax_np(nsx), ax_ao(nsx), ax_at(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: ord_a(nac)
      integer,  intent(in)    :: ky
      real(dp), intent(in)    :: ytens(ny1, ny1, ny3), jc, cq(nqc), dens(nao, nao)
      real(dp), intent(inout) :: gp(3, natm, nlane)

      integer(int64) :: t

      if (npq <= 0 .or. nac <= 0) return
      do concurrent(t = 0:int(npq, int64)*nac - 1)
         call {itm}(t, nsh, nao, sh_ao, sh_at, sh_r, npq, q_a, q_b, q_o, q_n, npp, pd, &
                   nsx, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, nac, ord_a, &
                   natm, nlane, ky, ny1, ny3, ytens, jc, nqc, qoff, cq, dens, gp)
      end do
   end subroutine {drv}

   pure subroutine {itm}(t, nsh, nao, sh_ao, sh_at, sh_r, npq, q_a, q_b, q_o, q_n, npp, pd, &
                         nsx, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, nac, ord_a, &
                         natm, nlane, ky, ny1, ny3, ytens, jc, nqc, qoff, cq, dens, gp)
      !! One (shell pair, auxiliary shell) triple: its Y block, the J part
      !! jc c_P D_mn formed here from c and D (never stored as a tensor) plus
      !! the K part read from `ytens` (ky = 1: Y_K itself; ky = 2: its one-
      !! sided half M, Y_K = M + M^T; ky = 0: none) -- weighted 2 for an
      !! off-diagonal pair, which stands for both (mu,nu) and (nu,mu) -- the
      !! contracted derivative, and nine atomic adds into lane mod(t, nlane)+1
      !! of `gp`. A triple on one atom is skipped: its three derivatives sum
      !! to zero on that atom.
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer,  intent(in)    :: nsh, nao, npq, npp, nsx, mnx, nac, natm, nlane, ny1, ny3, nqc, qoff
      integer,  intent(in)    :: sh_ao(nsh), sh_at(nsh)
      real(dp), intent(in)    :: sh_r(3, nsh)
      integer,  intent(in)    :: q_a(npq), q_b(npq), q_o(npq), q_n(npq)
      real(dp), intent(in)    :: pd(PD_W, npp)
      integer,  intent(in)    :: ax_np(nsx), ax_ao(nsx), ax_at(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: ord_a(nac)
      integer,  intent(in)    :: ky
      real(dp), intent(in)    :: ytens(ny1, ny1, ny3), jc, cq(nqc), dens(nao, nao)
      real(dp), intent(inout) :: gp(3, natm, nlane)

      integer  :: ip, kk, ii, jj, a, c, e, mu, nu, pq, lane, ia, ib, ic
      real(dp) :: yb({na*nb*nc}), ga(3), gb(3), wt, yv

      ip = int(mod(t, int(npq, int64))) + 1
      kk = ord_a(int(t/npq) + 1)
      ii = q_a(ip); jj = q_b(ip)
      ia = sh_at(ii); ib = sh_at(jj); ic = ax_at(kk)
      if (ia == ib .and. ia == ic) return
      wt = merge(1.0_dp, 2.0_dp, ii == jj)
      do e = 0, {nc - 1}
         pq = ax_ao(kk) + e
         do c = 0, {nb - 1}
            nu = sh_ao(jj) + c
            do a = 0, {na - 1}
               mu = sh_ao(ii) + a
               yv = jc*cq(qoff + pq)*dens(mu, nu)
               if (ky == 1) then
                  yv = yv + ytens(mu, nu, pq)
               else if (ky == 2) then
                  yv = yv + ytens(mu, nu, pq) + ytens(nu, mu, pq)
               end if
               yb(1 + a + {na}*(c + {nb}*e)) = wt*yv
            end do
         end do
      end do

      call {kern}(q_n(ip), q_o(ip), pd, sh_r(:, ii), sh_r(:, jj), &
                 ax_np(kk), ax_e(:, kk), ax_c(:, kk), ax_r(:, kk), yb, ga, gb)

      lane = int(mod(t, int(nlane, int64))) + 1
      !$acc atomic update
      gp(1, ia, lane) = gp(1, ia, lane) + ga(1)
      !$acc atomic update
      gp(2, ia, lane) = gp(2, ia, lane) + ga(2)
      !$acc atomic update
      gp(3, ia, lane) = gp(3, ia, lane) + ga(3)
      !$acc atomic update
      gp(1, ib, lane) = gp(1, ib, lane) + gb(1)
      !$acc atomic update
      gp(2, ib, lane) = gp(2, ib, lane) + gb(2)
      !$acc atomic update
      gp(3, ib, lane) = gp(3, ib, lane) + gb(3)
      !$acc atomic update
      gp(1, ic, lane) = gp(1, ic, lane) - (ga(1) + gb(1))
      !$acc atomic update
      gp(2, ic, lane) = gp(2, ic, lane) - (ga(2) + gb(2))
      !$acc atomic update
      gp(3, ic, lane) = gp(3, ic, lane) - (ga(3) + gb(3))
   end subroutine {itm}
{body}
end module {mod}
"""
    return mod, text


def emit_2c_module(la, lc, L, LA, order, order_val, order_arg, unroll_max, roots_max):
    na, nc = ncart(la), ncart(lc)
    tag = f"{la}{lc}"
    mod, drv, itm, kern = f"trc_dfg_2c_k{tag}", f"dfg_2c_k{tag}", f"dfg_2c_i{tag}", f"dfg_2c{tag}"
    body, tabs, nroots = emit_kernel(la, 0, lc, False, order, unroll_max, roots_max)
    text = MOD_HEADER.format(title=f"Class ({la}|{lc}) two-centre Coulomb derivative, contracted.",
                             L=L, LA=LA, order_arg=order_arg, order_val=order_val)
    text += f"""module {mod}
   use pic_types, only: dp, int64
{ri.roots_table_use(nroots)}   implicit none
   private
   public :: {drv}, {kern}

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp
{tabs}
contains

   subroutine {drv}(nsx, nax, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, &
                    nq, q_i, q_j, natm, nlane, zmat, gp)
      !! ({la}|{lc}) launch: `do concurrent` over a LIST of auxiliary shell
      !! pairs (i >= j) of this class, built on the host
      !! (trc_df_grad_kernels.F90's df_grad_2c_driver).
      integer,  intent(in)    :: nsx, nax, mnx, nq, natm, nlane
      integer,  intent(in)    :: ax_np(nsx), ax_ao(nsx), ax_at(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: q_i(nq), q_j(nq)
      real(dp), intent(in)    :: zmat(nax, nax)
      real(dp), intent(inout) :: gp(3, natm, nlane)

      integer :: t

      if (nq <= 0) return
      do concurrent(t = 1:nq)
         call {itm}(t, nsx, nax, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, &
                   nq, q_i, q_j, natm, nlane, zmat, gp)
      end do
   end subroutine {drv}

   pure subroutine {itm}(t, nsx, nax, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, &
                         nq, q_i, q_j, natm, nlane, zmat, gp)
      !! One auxiliary shell pair: its Z block (weighted 2 off the diagonal),
      !! the contracted derivative on the bra centre, and its negative on the
      !! ket centre. A pair on one atom is skipped.
      !$acc routine seq
      integer,  intent(in)    :: t, nsx, nax, mnx, nq, natm, nlane
      integer,  intent(in)    :: ax_np(nsx), ax_ao(nsx), ax_at(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: q_i(nq), q_j(nq)
      real(dp), intent(in)    :: zmat(nax, nax)
      real(dp), intent(inout) :: gp(3, natm, nlane)

      integer  :: ii, jj, a, c, mu, nu, lane, ia, ic
      real(dp) :: zb({na*nc}), ga(3), wt

      ii = q_i(t); jj = q_j(t)
      ia = ax_at(ii); ic = ax_at(jj)
      if (ia == ic) return
      wt = merge(1.0_dp, 2.0_dp, ii == jj)
      do c = 0, {nc - 1}
         nu = ax_ao(jj) + c
         do a = 0, {na - 1}
            mu = ax_ao(ii) + a
            zb(1 + a + {na}*c) = wt*zmat(mu, nu)
         end do
      end do

      call {kern}(ax_np(ii), ax_e(:, ii), ax_c(:, ii), ax_r(:, ii), &
                 ax_np(jj), ax_e(:, jj), ax_c(:, jj), ax_r(:, jj), zb, ga)

      lane = mod(t - 1, nlane) + 1
      !$acc atomic update
      gp(1, ia, lane) = gp(1, ia, lane) + ga(1)
      !$acc atomic update
      gp(2, ia, lane) = gp(2, ia, lane) + ga(2)
      !$acc atomic update
      gp(3, ia, lane) = gp(3, ia, lane) + ga(3)
      !$acc atomic update
      gp(1, ic, lane) = gp(1, ic, lane) - ga(1)
      !$acc atomic update
      gp(2, ic, lane) = gp(2, ic, lane) - ga(2)
      !$acc atomic update
      gp(3, ic, lane) = gp(3, ic, lane) - ga(3)
   end subroutine {itm}
{body}
end module {mod}
"""
    return mod, text


def emit_handoff(two, three, L, LA, order_arg):
    use2 = "\n".join(f"   use trc_dfg_2c_k{a}{c}, only: dfg_2c_k{a}{c}, dfg_2c{a}{c}" for a, c in two)
    use3 = "\n".join(f"   use trc_dfg_3c_k{a}{b}{c}, only: dfg_3c_k{a}{b}{c}, dfg_3c{a}{b}{c}"
                     for a, b, c in three)

    # Canonical classes only: la >= lb for the orbital pair, la >= lc for the
    # auxiliary pair; tri(a, b) = a(a+1)/2 + b numbers them.
    def tri(a, b):
        return a*(a + 1)//2 + b

    def k3(a, b, c):
        return tri(a, b)*(LA + 1) + c

    case2 = "\n".join(
        f"      case ({tri(a, c)}); call dfg_2c_k{a}{c}(nsx, nax, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, &\n"
        f"                          nq, q_i, q_j, natm, nlane, zmat, gp)"
        for a, c in two)
    case3 = "\n".join(
        f"      case ({k3(a, b, c)}); call dfg_3c_k{a}{b}{c}(nsh, nao, sh_ao, sh_at, sh_r, npq, q_a, q_b, q_o, q_n, &\n"
        f"                          npp, pd, nsx, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, nac, ord_a, &\n"
        f"                          natm, nlane, ky, ny1, ny3, ytens, jc, nqc, qoff, cq, dens, gp)"
        for a, b, c in three)
    blk2 = "\n".join(
        f"      case ({tri(a, c)}); call dfg_2c{a}{c}(npa, ea, ca, ra, npc, ec, cc_, rc, y, ga)"
        for a, c in two)
    blk3 = "\n".join(
        f"      case ({k3(a, b, c)}); call dfg_3c{a}{b}{c}(npp, 0, pd, ra, rb, npc, ec, cc_, rc, y, ga, gb)"
        for a, b, c in three)
    npcls = tri(L, L) + 1
    nacls = tri(LA, LA) + 1

    return f"""!
! Density-fitting DERIVATIVE hand-off: the pair plan and the two gradient
! drivers trc_scf_grad_df.F90 calls, and two host-side class dispatchers for
! tests.
!
! GENERATED by scripts/gen_df_grad.py --lmax {L} --lmax-aux {LA}{order_arg} --split -- do not edit.
!
! The kernels live one class per module (trc_dfg_2c_k*.F90, trc_dfg_3c_k*.F90),
! see gen_df_grad.py's docstring. This file is HOST code only: it groups work
! by angular momentum and dispatches to the one class module covering each
! group, exactly like trc_df_kernels.F90 does for the energy.
!
! CANONICAL CLASSES. Only la >= lb (three-centre) and la >= lc (two-centre)
! are generated: a pair whose first shell has the lower l is handed to the
! kernel the other way round (the higher-l shell first), which the
! integrals' symmetry in the two shells -- and Y's and Z's -- makes exact.
!
module trc_df_grad_kernels
   use pic_types, only: dp
{use2}
{use3}
   implicit none
   private
   public :: df_grad_3c_plan, df_grad_2c_driver, df_grad_3c_driver, df_grad_2c_block, df_grad_3c_block
   public :: DFG_LMAX, DFG_LMAX_AUX, DFG_PD_W, DFG_NPCLS

   integer, parameter :: DFG_LMAX = {L}, DFG_LMAX_AUX = {LA}
   integer, parameter :: DFG_PD_W = {PD_W}
   !! Doubles per primitive pair in the plan's pd: zeta, K_ab, P(3), alpha,
   !! beta, one pad (64 bytes).
   integer, parameter :: DFG_NPCLS = {npcls}
   !! Canonical orbital pair classes (la >= lb), numbered la(la+1)/2 + lb.
   integer, parameter :: NACLS = {nacls}
   !! Canonical auxiliary pair classes (la >= lc), numbered the same way.
   integer, parameter :: NLANE = 16
   !! Lane-private gradient copies (gp(3, natm, NLANE)): a work item adds
   !! into lane mod(t, NLANE)+1, so neighbouring threads rarely contend for
   !! the same atomic; the lanes are summed once at the end of each driver.

contains

   subroutine df_grad_3c_plan(nsh, mnp, sh_l, sh_np, sh_e, sh_c, sh_r, npair, p_i, p_j, &
                              coff, q_a, q_b, q_o, q_n, pd)
      !! HOST, once per gradient: the three-centre kernels' view of the
      !! energy's screened pair list (p_i >= p_j). Each pair is turned so its
      !! first shell q_a has the higher l (its canonical class), the pairs
      !! are grouped by class -- class k is q_*(coff(k)+1:coff(k+1)) -- and,
      !! within a class, sorted by primitive-pair count, so the threads of
      !! one warp run the same primitive trip count. pd(:, q_o(i)+1 ..
      !! q_o(i)+q_n(i)) holds pair i's primitive pairs, computed once here
      !! rather than once per (pair, auxiliary shell) item. The caller copies
      !! q_* and pd to the device.
      integer,  intent(in)  :: nsh, mnp, npair
      integer,  intent(in)  :: sh_l(nsh), sh_np(nsh), p_i(npair), p_j(npair)
      real(dp), intent(in)  :: sh_e(mnp, nsh), sh_c(mnp, nsh), sh_r(3, nsh)
      integer,  intent(out) :: coff(0:DFG_NPCLS)
      integer,  allocatable, intent(out) :: q_a(:), q_b(:), q_o(:), q_n(:)
      real(dp), allocatable, intent(out) :: pd(:, :)

      integer :: i, a, b, la, lb, key, nkey, nmax, k, npp
      integer, allocatable :: kcnt(:), kidx(:)

      nmax = mnp*mnp
      nkey = DFG_NPCLS*(nmax + 1)
      allocate (kcnt(0:nkey), kidx(npair), q_a(npair), q_b(npair), q_o(npair), q_n(npair))
      kcnt = 0
      do i = 1, npair
         a = p_i(i); b = p_j(i)
         if (sh_l(a) < sh_l(b)) then
            a = p_j(i); b = p_i(i)
         end if
         la = sh_l(a); lb = sh_l(b)
         kidx(i) = (la*(la + 1)/2 + lb)*(nmax + 1) + sh_np(a)*sh_np(b)
         kcnt(kidx(i) + 1) = kcnt(kidx(i) + 1) + 1
      end do
      do key = 1, nkey
         kcnt(key) = kcnt(key) + kcnt(key - 1)
      end do
      do k = 0, DFG_NPCLS
         coff(k) = kcnt(k*(nmax + 1))
      end do
      do i = 1, npair
         a = p_i(i); b = p_j(i)
         if (sh_l(a) < sh_l(b)) then
            a = p_j(i); b = p_i(i)
         end if
         kcnt(kidx(i)) = kcnt(kidx(i)) + 1
         k = kcnt(kidx(i))
         q_a(k) = a; q_b(k) = b
         q_n(k) = sh_np(a)*sh_np(b)
      end do
      npp = 0
      do k = 1, npair
         q_o(k) = npp
         npp = npp + q_n(k)
      end do
      allocate (pd(DFG_PD_W, max(npp, 1)))
      do k = 1, npair
         call pair_prims(sh_np(q_a(k)), sh_e(:, q_a(k)), sh_c(:, q_a(k)), sh_r(:, q_a(k)), &
                         sh_np(q_b(k)), sh_e(:, q_b(k)), sh_c(:, q_b(k)), sh_r(:, q_b(k)), &
                         pd(:, q_o(k) + 1:q_o(k) + q_n(k)))
      end do
      deallocate (kcnt, kidx)
   end subroutine df_grad_3c_plan

   subroutine pair_prims(npa, ea, ca, ra, npb, eb, cb, rb, pd)
      !! HOST: one shell pair's primitive pairs, a outer, b inner -- the
      !! plan's and df_grad_3c_block's single source of these numbers.
      integer,  intent(in)  :: npa, npb
      real(dp), intent(in)  :: ea(npa), ca(npa), ra(3), eb(npb), cb(npb), rb(3)
      real(dp), intent(out) :: pd(DFG_PD_W, npa*npb)
      integer  :: ip, jp, k
      real(dp) :: al, be, zeta, rab2
      rab2 = (ra(1) - rb(1))**2 + (ra(2) - rb(2))**2 + (ra(3) - rb(3))**2
      k = 0
      do ip = 1, npa
         al = ea(ip)
         do jp = 1, npb
            be = eb(jp)
            zeta = al + be
            k = k + 1
            pd(1, k) = zeta
            pd(2, k) = ca(ip)*cb(jp)*exp(-al*be/zeta*rab2)
            pd(3, k) = (al*ra(1) + be*rb(1))/zeta
            pd(4, k) = (al*ra(2) + be*rb(2))/zeta
            pd(5, k) = (al*ra(3) + be*rb(3))/zeta
            pd(6, k) = al
            pd(7, k) = be
            pd(8, k) = 0.0_dp
         end do
      end do
   end subroutine pair_prims

   subroutine df_grad_3c_driver(nsh, nao, sh_ao, sh_at, sh_r, coff, npair, q_a, q_b, q_o, q_n, npp, pd, &
                                nsx, mnx, ax_l, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, natm, &
                                ky, ny1, ny3, ytens, jc, nqc, qoff, cq, dens, grad)
      !! grad += sum_(mu nu P) Y(mu,nu,P) d(mu nu|P)/dR over the planned
      !! shell pairs (df_grad_3c_plan) and every shell of the auxiliary
      !! (sub)basis, with
      !!
      !!     Y(mu,nu,P) = jc cq(qoff+P) dens(mu,nu)
      !!                  + ytens(mu,nu,P)                    (ky = 1)
      !!                  + ytens(mu,nu,P) + ytens(nu,mu,P)   (ky = 2)
      !!
      !! -- the J part formed inside the kernels from c and D, never stored,
      !! so a J-only gradient (ky = 0; ytens then a dummy, ny1 = ny3 = 1)
      !! needs no (n, n, naux) tensor at all. P is the auxiliary (sub)basis's
      !! OWN AO numbering (a subset basis from trc_basis_t%subset numbers
      !! from 1); qoff maps it into cq. With ky = 1 `ytens` is the K part
      !! itself (symmetric in mu, nu); with ky = 2 its one-sided half M,
      !! which trc_scf_grad_df.F90 builds from X's packed rows j <= i when
      !! the occupied orbitals come in more than one chunk. Every array but
      !! `coff` is device-resident; `grad` is added into.
      integer,  intent(in)    :: nsh, nao, npair, npp, nsx, mnx, natm, ny1, ny3, nqc, qoff
      integer,  intent(in)    :: sh_ao(nsh), sh_at(nsh)
      real(dp), intent(in)    :: sh_r(3, nsh)
      integer,  intent(in)    :: coff(0:DFG_NPCLS)
      integer,  intent(in)    :: q_a(npair), q_b(npair), q_o(npair), q_n(npair)
      real(dp), intent(in)    :: pd(DFG_PD_W, npp)
      integer,  intent(in)    :: ax_l(nsx), ax_np(nsx), ax_ao(nsx), ax_at(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: ky
      real(dp), intent(in)    :: ytens(ny1, ny1, ny3), jc, cq(nqc), dens(nao, nao)
      real(dp), intent(inout) :: grad(3, natm)

      integer :: i, l, la, lb, lc, k, ga_, gc_, gs, key, nkey
      integer :: aoff(0:DFG_LMAX_AUX + 1)
      integer, allocatable :: orda(:), kcnt(:)
      real(dp) :: gsum
      real(dp), allocatable :: gp(:, :, :)

      ! The auxiliary shells grouped by l (one launch per class each) and,
      ! within an l, sorted by primitive count -- the plan's warp-uniformity
      ! reason again.
      nkey = (DFG_LMAX_AUX + 1)*(mnx + 1)
      allocate (orda(nsx), kcnt(0:nkey))
      kcnt = 0
      do i = 1, nsx
         key = ax_l(i)*(mnx + 1) + ax_np(i)
         kcnt(key + 1) = kcnt(key + 1) + 1
      end do
      do key = 1, nkey
         kcnt(key) = kcnt(key) + kcnt(key - 1)
      end do
      do l = 0, DFG_LMAX_AUX + 1
         aoff(l) = kcnt(l*(mnx + 1))
      end do
      do i = 1, nsx
         key = ax_l(i)*(mnx + 1) + ax_np(i)
         kcnt(key) = kcnt(key) + 1
         orda(kcnt(key)) = i
      end do
      deallocate (kcnt)
      !$acc enter data copyin(orda)

      allocate (gp(3, natm, NLANE))
      !$acc enter data create(gp)
      do concurrent(gs=1:NLANE, ga_=1:natm, gc_=1:3)
         gp(gc_, ga_, gs) = 0.0_dp
      end do

      do la = 0, DFG_LMAX
         do lb = 0, la
            k = la*(la + 1)/2 + lb
            if (coff(k + 1) <= coff(k)) cycle
            do lc = 0, DFG_LMAX_AUX
               if (aoff(lc + 1) <= aoff(lc)) cycle
               call df_grad_3c_dispatch(k, lc, nsh, nao, sh_ao, sh_at, sh_r, coff(k + 1) - coff(k), &
                                        q_a(coff(k) + 1:coff(k + 1)), q_b(coff(k) + 1:coff(k + 1)), &
                                        q_o(coff(k) + 1:coff(k + 1)), q_n(coff(k) + 1:coff(k + 1)), npp, pd, &
                                        nsx, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, &
                                        aoff(lc + 1) - aoff(lc), orda(aoff(lc) + 1:aoff(lc + 1)), &
                                        natm, NLANE, ky, ny1, ny3, ytens, jc, nqc, qoff, cq, dens, gp)
            end do
         end do
      end do

      do concurrent(ga_=1:natm, gc_=1:3) local(gsum, gs)
         gsum = 0.0_dp
         do gs = 1, NLANE
            gsum = gsum + gp(gc_, ga_, gs)
         end do
         grad(gc_, ga_) = grad(gc_, ga_) + gsum
      end do
      !$acc exit data delete(gp, orda)
      deallocate (gp, orda)
   end subroutine df_grad_3c_driver

   subroutine df_grad_2c_driver(nsx, nax, mnx, ax_l, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, &
                                natm, zmat, grad)
      !! grad += sum_(P Q) Z(P,Q) d(P|Q)/dR over every auxiliary shell pair
      !! (each once, weighted inside the item, the higher-l shell first: its
      !! canonical class). `zmat` symmetric, device-resident.
      integer,  intent(in)    :: nsx, nax, mnx, natm
      integer,  intent(in)    :: ax_l(nsx), ax_np(nsx), ax_ao(nsx), ax_at(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      real(dp), intent(in)    :: zmat(nax, nax)
      real(dp), intent(inout) :: grad(3, natm)

      integer :: i, j, a, c, l, la, lc, nq, ga_, gc_, gs
      integer :: cnt(0:NACLS - 1), off(0:NACLS), cur(0:NACLS - 1)
      integer, allocatable :: qi(:), qj(:)
      real(dp) :: gsum
      real(dp), allocatable :: gp(:, :, :)

      cnt = 0
      do i = 1, nsx
         do j = 1, i
            if (ax_at(i) == ax_at(j)) cycle
            a = max(ax_l(i), ax_l(j)); c = min(ax_l(i), ax_l(j))
            l = a*(a + 1)/2 + c
            cnt(l) = cnt(l) + 1
         end do
      end do
      off(0) = 0
      do l = 0, NACLS - 1
         off(l + 1) = off(l) + cnt(l)
      end do
      nq = off(NACLS)
      allocate (qi(max(nq, 1)), qj(max(nq, 1)))
      cur = off(0:NACLS - 1)
      do i = 1, nsx
         do j = 1, i
            if (ax_at(i) == ax_at(j)) cycle
            a = i; c = j
            if (ax_l(i) < ax_l(j)) then
               a = j; c = i
            end if
            l = ax_l(a)*(ax_l(a) + 1)/2 + ax_l(c)
            cur(l) = cur(l) + 1
            qi(cur(l)) = a
            qj(cur(l)) = c
         end do
      end do
      !$acc enter data copyin(qi, qj)

      allocate (gp(3, natm, NLANE))
      !$acc enter data create(gp)
      do concurrent(gs=1:NLANE, ga_=1:natm, gc_=1:3)
         gp(gc_, ga_, gs) = 0.0_dp
      end do

      do la = 0, DFG_LMAX_AUX
         do lc = 0, la
            l = la*(la + 1)/2 + lc
            if (off(l + 1) <= off(l)) cycle
            call df_grad_2c_dispatch(l, nsx, nax, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, &
                                     off(l + 1) - off(l), qi(off(l) + 1:off(l + 1)), qj(off(l) + 1:off(l + 1)), &
                                     natm, NLANE, zmat, gp)
         end do
      end do

      do concurrent(ga_=1:natm, gc_=1:3) local(gsum, gs)
         gsum = 0.0_dp
         do gs = 1, NLANE
            gsum = gsum + gp(gc_, ga_, gs)
         end do
         grad(gc_, ga_) = grad(gc_, ga_) + gsum
      end do
      !$acc exit data delete(gp, qi, qj)
      deallocate (gp, qi, qj)
   end subroutine df_grad_2c_driver

   subroutine df_grad_3c_dispatch(k, lc, nsh, nao, sh_ao, sh_at, sh_r, npq, q_a, q_b, q_o, q_n, npp, pd, &
                                  nsx, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, nac, ord_a, &
                                  natm, nlane, ky, ny1, ny3, ytens, jc, nqc, qoff, cq, dens, gp)
      integer,  intent(in)    :: k, lc, nsh, nao, npq, npp, nsx, mnx, nac, natm, nlane, ny1, ny3, nqc, qoff
      integer,  intent(in)    :: sh_ao(nsh), sh_at(nsh)
      real(dp), intent(in)    :: sh_r(3, nsh)
      integer,  intent(in)    :: q_a(npq), q_b(npq), q_o(npq), q_n(npq)
      real(dp), intent(in)    :: pd(DFG_PD_W, npp)
      integer,  intent(in)    :: ax_np(nsx), ax_ao(nsx), ax_at(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: ord_a(nac)
      integer,  intent(in)    :: ky
      real(dp), intent(in)    :: ytens(ny1, ny1, ny3), jc, cq(nqc), dens(nao, nao)
      real(dp), intent(inout) :: gp(3, natm, nlane)
      select case (k*(DFG_LMAX_AUX + 1) + lc)
{case3}
      case default
         error stop "trc_df_grad_kernels: no three-centre derivative class for this angular momentum"
      end select
   end subroutine df_grad_3c_dispatch

   subroutine df_grad_2c_dispatch(l, nsx, nax, mnx, ax_np, ax_ao, ax_at, ax_e, ax_c, ax_r, &
                                  nq, q_i, q_j, natm, nlane, zmat, gp)
      integer,  intent(in)    :: l, nsx, nax, mnx, nq, natm, nlane
      integer,  intent(in)    :: ax_np(nsx), ax_ao(nsx), ax_at(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      integer,  intent(in)    :: q_i(nq), q_j(nq)
      real(dp), intent(in)    :: zmat(nax, nax)
      real(dp), intent(inout) :: gp(3, natm, nlane)
      select case (l)
{case2}
      case default
         error stop "trc_df_grad_kernels: no two-centre derivative class for this angular momentum"
      end select
   end subroutine df_grad_2c_dispatch

   subroutine df_grad_3c_block(la, lb, lc, npa, ea, ca, ra, npb, eb, cb, rb, npc, ec, cc_, rc, y, ga, gb)
      !! HOST, for tests: one three-centre class kernel, by angular momentum.
      !! With a one-hot `y` this returns one raw derivative integral
      !! (d/dA in ga, d/dB in gb) -- the same code the driver runs, including
      !! its canonical turn for la < lb (B handed in first, y transposed to
      !! match, ga and gb swapped back).
      integer,  intent(in)  :: la, lb, lc, npa, npb, npc
      real(dp), intent(in)  :: ea(npa), ca(npa), ra(3), eb(npb), cb(npb), rb(3)
      real(dp), intent(in)  :: ec(npc), cc_(npc), rc(3), y(*)
      real(dp), intent(out) :: ga(3), gb(3)
      integer :: na, nb, nc, ia, ib, ic, npp
      real(dp), allocatable :: pd(:, :), yt(:)
      real(dp) :: g1(3), g2(3)
      npp = npa*npb
      allocate (pd(DFG_PD_W, npp))
      if (la >= lb) then
         call pair_prims(npa, ea, ca, ra, npb, eb, cb, rb, pd)
         call block3(la, lb, lc, npp, pd, ra, rb, npc, ec, cc_, rc, y, ga, gb)
      else
         na = (la + 1)*(la + 2)/2; nb = (lb + 1)*(lb + 2)/2; nc = (lc + 1)*(lc + 2)/2
         allocate (yt(na*nb*nc))
         do ic = 0, nc - 1
            do ib = 0, nb - 1
               do ia = 0, na - 1
                  yt(1 + ib + nb*(ia + na*ic)) = y(1 + ia + na*(ib + nb*ic))
               end do
            end do
         end do
         call pair_prims(npb, eb, cb, rb, npa, ea, ca, ra, pd)
         call block3(lb, la, lc, npp, pd, rb, ra, npc, ec, cc_, rc, yt, g1, g2)
         ga = g2
         gb = g1
      end if
   end subroutine df_grad_3c_block

   subroutine block3(la, lb, lc, npp, pd, ra, rb, npc, ec, cc_, rc, y, ga, gb)
      !! HOST: the canonical (la >= lb) three-centre class kernel by l.
      integer,  intent(in)  :: la, lb, lc, npp, npc
      real(dp), intent(in)  :: pd(DFG_PD_W, npp), ra(3), rb(3), ec(npc), cc_(npc), rc(3), y(*)
      real(dp), intent(out) :: ga(3), gb(3)
      select case ((la*(la + 1)/2 + lb)*(DFG_LMAX_AUX + 1) + lc)
{blk3}
      case default
         error stop "trc_df_grad_kernels: no three-centre derivative class for this angular momentum"
      end select
   end subroutine block3

   subroutine df_grad_2c_block(la, lc, npa, ea, ca, ra, npc, ec, cc_, rc, y, ga)
      !! HOST, for tests: one two-centre class kernel (d/dA, contracted with
      !! y), with the driver's canonical turn for la < lc: C handed in first,
      !! y transposed, and d/dA = -d/dC by translational invariance.
      integer,  intent(in)  :: la, lc, npa, npc
      real(dp), intent(in)  :: ea(npa), ca(npa), ra(3), ec(npc), cc_(npc), rc(3), y(*)
      real(dp), intent(out) :: ga(3)
      integer :: na, nc, ia, ic
      real(dp), allocatable :: yt(:)
      if (la >= lc) then
         call block2(la, lc, npa, ea, ca, ra, npc, ec, cc_, rc, y, ga)
      else
         na = (la + 1)*(la + 2)/2; nc = (lc + 1)*(lc + 2)/2
         allocate (yt(na*nc))
         do ic = 0, nc - 1
            do ia = 0, na - 1
               yt(1 + ic + nc*ia) = y(1 + ia + na*ic)
            end do
         end do
         call block2(lc, la, npc, ec, cc_, rc, npa, ea, ca, ra, yt, ga)
         ga = -ga
      end if
   end subroutine df_grad_2c_block

   subroutine block2(la, lc, npa, ea, ca, ra, npc, ec, cc_, rc, y, ga)
      !! HOST: the canonical (la >= lc) two-centre class kernel by l.
      integer,  intent(in)  :: la, lc, npa, npc
      real(dp), intent(in)  :: ea(npa), ca(npa), ra(3), ec(npc), cc_(npc), rc(3), y(*)
      real(dp), intent(out) :: ga(3)
      select case (la*(la + 1)/2 + lc)
{blk2}
      case default
         error stop "trc_df_grad_kernels: no two-centre derivative class for this angular momentum"
      end select
   end subroutine block2

end module trc_df_grad_kernels
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=2, help="max l on the orbital shells")
    ap.add_argument("--lmax-aux", type=int, default=None, help="max l on the auxiliary shells (default lmax + 2)")
    ap.add_argument("--order", choices=("lexicographic", "gamess"), default="lexicographic")
    ap.add_argument("--split", metavar="DIR", required=True,
                    help="emit one module per class into DIR, plus trc_df_grad_kernels.F90")
    ap.add_argument("--unroll-max-comp", default="0",
                    help="unroll a class's root/component loops when na*nb*nc is at most this "
                         "(the same TERCO_RYS_UNROLL_MAX gen_df.py takes)")
    ap.add_argument("--unroll-roots-max-comp", default=None,
                    help="also unroll the root loop only when na*nb*nc is at most this (default: "
                         "--unroll-max-comp); classes between the two keep a run-time root loop around "
                         "unrolled components")
    args = ap.parse_args()
    L = args.lmax
    LA = args.lmax_aux if args.lmax_aux is not None else L + 2
    order = args.order
    order_val = 0 if order == "lexicographic" else 1
    order_arg = "" if order == "lexicographic" else f" --order {order}"
    unroll_max = _parse_unroll_max(args.unroll_max_comp)
    roots_max = unroll_max if args.unroll_roots_max_comp is None else _parse_unroll_max(args.unroll_roots_max_comp)
    os.makedirs(args.split, exist_ok=True)

    two, three = [], []
    for la in range(LA + 1):
        for lc in range(la + 1):
            mod, text = emit_2c_module(la, lc, L, LA, order, order_val, order_arg, unroll_max, roots_max)
            with open(os.path.join(args.split, mod + ".F90"), "w") as fh:
                fh.write(text)
            two.append((la, lc))
    for la in range(L + 1):
        for lb in range(la + 1):
            for lc in range(LA + 1):
                mod, text = emit_3c_module(la, lb, lc, L, LA, order, order_val, order_arg, unroll_max,
                                           roots_max)
                with open(os.path.join(args.split, mod + ".F90"), "w") as fh:
                    fh.write(text)
                three.append((la, lb, lc))
    with open(os.path.join(args.split, "trc_df_grad_kernels.F90"), "w") as fh:
        fh.write(emit_handoff(two, three, L, LA, order_arg))
    print(f"wrote {len(two)} two-centre + {len(three)} three-centre derivative class modules "
          f"+ trc_df_grad_kernels.F90 into {args.split}")


if __name__ == "__main__":
    main()
