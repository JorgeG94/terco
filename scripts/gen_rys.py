#!/usr/bin/env python3
"""Emit a Rys kernel for one class, by splicing terco's own scalar kernel.

WHY SPLICE RATHER THAN GENERATE WHOLE
--------------------------------------
The comparison this exists to serve is HGP against Rys with everything else
held fixed. Reusing the scalar HGP kernel's prologue (segment decode,
screening, the bin-interleaved primitive-pair loads, the primitive-quartet
prescreen) and its folded six-atomic digestion VERBATIM is what makes that
true: only the middle -- the integral itself -- differs, so a timing
difference is the scheme and nothing else.

This reads `src/generated/trc_pc_k{tag}.F90`, the file `gen_perclass.py
--split` already produced for the scalar (decontracted, single-column)
kernel `pcs{tag}`/`pcsi{tag}`. General contraction is not this file's
concern at all: `trc_decontract` splits it away before any kernel runs, so
every class file -- HGP or Rys -- only ever needs the scalar form.

THE CUT POINTS
---------------
`pcsi{tag}` falls into three parts, found by text search rather than by
counting lines, so a reflow of the surrounding generator does not silently
misalign the splice:

  1. HEAD, up to and including `tval = rho*(pqx*pqx + ...)`. Segment decode,
     screens, the pair-table offsets, `abx..cdz`, and the bra/ket primitive
     loop through everything two shell-pair centres and one primitive
     quartet need. Identical for HGP and Rys, kept verbatim.
  2. The INTEGRAL itself -- HGP's Boys evaluation, VRR and the accumulation
     into its pre-HRR block. Replaced by the Rys root lookup, the 2D
     recurrence and the horizontal transfer, run per primitive quartet
     (HGP contracts over primitives and only then transforms once; Rys
     builds the shifted class per primitive quartet instead -- see
     `docs/RYS.md`'s "the cause"). This writes the final Cartesian block
     straight into `vbuf`, so there is no separate HRR pass afterward.
  3. TAIL: the degeneracy/AO-offset decode (`! Single column on every
     side:`), the digestion (`#ifdef TRC_NO_DIGEST ... #endif`), and
     `end subroutine`. Scheme-independent, kept verbatim. HGP's own
     "--- HRR ---" text, which sits between the decode and the digestion in
     the source file, is dropped: Rys's vbuf is already the transferred
     block by the time the decode runs.

THE ALGEBRA
-----------
From gamess-libERI's int*_rysgen, reduced: writing P-A and Q-C through PQ
makes the c00/cp00 coefficients collapse to one term each, so the kernel
needs only pax/qcx/pqx -- all of which the prologue already has, since HGP
needs them too.  `tval` IS the Rys argument and `pref` already carries
2*pi^(5/2)/(zeta*eta*sqrt(zpe))*(coefficients). This algebra was validated
against libcint through (dd|dd) on explore/rys-4c (not carried over here);
what checks it in this tree is `test/check_rys_vs_hgp.F90`, comparing the
Fock matrix this produces against HGP's on the same density.
"""
import argparse
import importlib.util
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


co = _load("co", os.path.join(HERE, "cartorder.py"))

#: Binomial row b, for the horizontal transfer written out rather than
#: looped -- LMAX stays 2, so a shift of 0, 1 or 2 is everything this ever
#: needs. Index [b][t] is C(b, t).
_BINOM = [[1], [1, 1], [1, 2, 1]]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("class_tag", help="four digits, e.g. 1111 or 2222")
    ap.add_argument("--order", choices=("lexicographic", "gamess"),
                    default="lexicographic",
                    help="must match the --order the scalar HGP file in "
                         "--src-dir was generated with; the digestion it is "
                         "spliced from already scatters into this order, "
                         "and the component tables below must agree with it")
    ap.add_argument("--src-dir", default="src/generated",
                    help="where trc_pc_k{tag}.F90 (the scalar HGP kernel) "
                         "already lives")
    ap.add_argument("-o", "--output", required=True,
                    help="required, and never a path under src/generated: "
                         "unlike the Rys root tables, these per-class "
                         "kernels are not committed -- CMake regenerates "
                         "them at configure time into the build directory, "
                         "same as a non-default TERCO_LMAX's HGP kernels")
    args = ap.parse_args()
    tag = args.class_tag
    la, lb, lc, ld = (int(c) for c in tag)
    ni, nk = la + lb, lc + ld
    nroots = (la + lb + lc + ld)//2 + 1
    dims = [len(co.components(l, args.order)) for l in (la, lb, lc, ld)]
    nv = dims[0]*dims[1]*dims[2]*dims[3]

    src_path = os.path.join(args.src_dir, f"trc_pc_k{tag}.F90")
    src = open(src_path).read()

    # Scope every search to pcsi{tag} itself (not pcs{tag}, its driver, and
    # not any pcr{tag}/pcri{tag} warp-role kernel the same file may also
    # carry for HGP -- ROLE_SPLITS in gen_perclass.py -- which would also
    # contain a `tval = rho*(pqx...` line and silently splice the wrong one).
    i_item = src.index(f"   pure subroutine pcsi{tag}(")
    i_end = src.index(f"end subroutine pcsi{tag}\n") + len(f"end subroutine pcsi{tag}\n")
    whole = src[:i_end]
    assert i_item < i_end

    tval_line = "tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)\n"
    i_tval = whole.index(tval_line) + len(tval_line)
    head = whole[:i_tval]

    decode_marker = "         ! Single column on every side:"
    hrr_marker = "         ! --- HRR ---\n"
    dig_marker = "#ifdef TRC_NO_DIGEST"
    i_decode = whole.index(decode_marker)
    i_hrr = whole.index(hrr_marker)
    i_dig = whole.index(dig_marker)
    decode = whole[i_decode:i_hrr]
    tail = whole[i_dig:]

    # The pre-HRR VRR accumulator's size, read back from its own
    # declaration rather than recomputed -- it is whatever gen_perclass.py
    # actually emitted, and asserting it here catches a mismatch as a clear
    # error instead of a wrong array bound three edits later.
    decl_re = re.compile(
        r"      real\(dp\) :: v\((\d+), 0:1\), g1\(\1\), vbuf\((\d+)\)\n")
    m = decl_re.search(head)
    assert m, "could not find the HGP v/g1/vbuf declaration to replace"
    nv_vrr, nv_in_file = int(m.group(1)), int(m.group(2))
    assert nv_in_file == nv, (
        f"vbuf({nv_in_file}) in {src_path} disagrees with "
        f"na*nb*nc*nd = {nv} computed here for --order {args.order}; "
        "the --order argument likely does not match the source file")

    decl_new = (
        f"      real(dp) :: vbuf({nv})\n"
        f"      real(dp) :: ru({nroots}), rw({nroots})\n"
        f"      real(dp) :: gx(0:{ni}, 0:{nk}), gy(0:{ni}, 0:{nk}), gz(0:{ni}, 0:{nk})\n"
        f"      real(dp) :: hx(0:{la}, 0:{lb}, 0:{lc}, 0:{ld})\n"
        f"      real(dp) :: hy(0:{la}, 0:{lb}, 0:{lc}, 0:{ld})\n"
        f"      real(dp) :: hz(0:{la}, 0:{lb}, 0:{lc}, 0:{ld})\n"
        f"      real(dp) :: u2, dinv, b00r, b10r, bp01r, f00\n"
        f"      real(dp) :: c00x, c00y, c00z, cp0x, cp0y, cp0z\n"
        f"      integer  :: rr, aa, bb, cc, dd2\n")
    head = head[:m.start()] + decl_new + head[m.end():]

    # Boys scratch, now unused.
    for pat in (r"      real\(dp\) :: f\(0:BOYS_MMAX\)\n",
                r"      integer  :: bi, bj, bbase\n",
                r"      real\(dp\) :: bx, bx2, b0, b1, b2, btt, bet\n"):
        head, n1 = re.subn(pat, "", head, count=1)
        assert n1 == 1, pat

    # HGP zeroes its pre-HRR accumulator once before the primitive loop;
    # Rys zeroes `vbuf` there instead, since nothing defers the transfer.
    zero_re = re.compile(r"(            do x = 1, )" + str(nv_vrr) +
                         r"(\n               )g1\(x\)( = 0\.0_dp\n            end do\n)")
    head, n1 = zero_re.subn(lambda g: f"{g.group(1)}{nv}{g.group(2)}vbuf(x){g.group(3)}", head, 1)
    assert n1 == 1, "could not find the g1 zero loop to retarget at vbuf"

    # No Boys function; the root lookup replaces it.
    head, n1 = re.subn(r"   use trc_boys, only:.*?BOYS_DT, BOYS_DTINV\n",
                       f"   use trc_boys, only: dp\n"
                       f"   use trc_rys, only: trc_rys_roots{nroots}_tab\n",
                       head, count=1, flags=re.S)
    assert n1 == 1

    def ren(t):
        t = t.replace(f"trc_pc_k{tag}", f"trc_pc_rys{tag}")
        t = t.replace(f"pcsi{tag}", f"rysi{tag}")
        t = t.replace(f"pcs{tag}", f"rys{tag}")
        return t

    head = re.sub(r"public :: pcs" + tag + r".*\n", f"public :: rys{tag}\n", head, count=1)
    head, decode, tail = ren(head), ren(decode), ren(tail)

    # The file-top comment still describes the HGP file it was read from;
    # replace it with one that says what this file actually is.
    old_header_re = re.compile(r"\A!\n.*?#endif\n", re.S)
    assert old_header_re.search(head)
    order_val = 0 if args.order == "lexicographic" else 1
    new_header = f"""!
! Class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}) kernel, Rys quadrature.
!
! GENERATED by scripts/gen_rys.py {tag} -- do not edit.
!
! Spliced from the scalar HGP kernel in src/generated/trc_pc_k{tag}.F90: the
! prologue (segment decode, screens, the primitive-pair loads, the
! primitive-quartet prescreen) and the folded six-atomic digestion are that
! file's text, verbatim; only the integral itself differs. See this
! generator's module docstring for the cut points and docs/RYS_PLAN.md for
! why the splice is done this way.
!
! TRC_CART_ORDER GUARD -- see gen_perclass.py's non-split header for the
! full explanation; same guard, carried over from the file this was spliced
! from, since the digestion below scatters into the same AO slot order.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
"""
    head = old_header_re.sub(new_header, head, count=1)

    # The Rys integral itself: roots, the 2D recurrence per axis, the
    # horizontal transfer per axis, and the outer product straight into
    # vbuf -- no deferred HRR, see the module docstring.
    middle = f"""
                  call trc_rys_roots{nroots}_tab(tval, ru, rw)

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

                     call rys2d_{tag}(c00x, cp0x, b00r, b10r, bp01r, 1.0_dp, gx)
                     call rys2d_{tag}(c00y, cp0y, b00r, b10r, bp01r, 1.0_dp, gy)
                     call rys2d_{tag}(c00z, cp0z, b00r, b10r, bp01r, f00,    gz)
                     call rys_hrr_{tag}(gx, abx, cdx, hx)
                     call rys_hrr_{tag}(gy, aby, cdy, hy)
                     call rys_hrr_{tag}(gz, abz, cdz, hz)

                     idx = 0
                     do dd2 = 1, {dims[3]}
                        do cc = 1, {dims[2]}
                           do bb = 1, {dims[1]}
                              do aa = 1, {dims[0]}
                                 idx = idx + 1
                                 vbuf(idx) = vbuf(idx) &
                                    + hx(cnx0_{tag}(aa), cnx1_{tag}(bb), cnx2_{tag}(cc), cnx3_{tag}(dd2)) &
                                    * hy(cny0_{tag}(aa), cny1_{tag}(bb), cny2_{tag}(cc), cny3_{tag}(dd2)) &
                                    * hz(cnz0_{tag}(aa), cnz1_{tag}(bb), cnz2_{tag}(cc), cnz3_{tag}(dd2))
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
"""

    # Component tables: aa/bb/cc/dd2 above index into these to pick up each
    # shell's (lx,ly,lz) triple in the BUILD'S OWN Cartesian order -- the
    # same order gen_perclass.py's digestion scattered into when it wrote
    # the tail this file reused verbatim.
    tabs = ""
    for nm, ax in (("cnx", 0), ("cny", 1), ("cnz", 2)):
        for which, l in enumerate((la, lb, lc, ld)):
            comps = co.components(l, args.order)
            vals = ", ".join(str(c[ax]) for c in comps)
            tabs += f"   integer, parameter :: {nm}{which}_{tag}({len(comps)}) = [{vals}]\n"

    helpers = f"""
{tabs}
contains

   pure subroutine rys2d_{tag}(c00, cp00, b00, b10, bp01, g00, g)
      !! One axis of 2D Rys integrals, all angular momentum on A and C.
      !$acc routine seq
      real(dp), intent(in)  :: c00, cp00, b00, b10, bp01, g00
      real(dp), intent(out) :: g(0:{ni}, 0:{nk})
      integer :: i, k

      g(0, 0) = g00
{"      g(1, 0) = c00*g00" if ni >= 1 else ""}
{f"      do i = 1, {ni - 1}" if ni >= 2 else ""}
{"         g(i + 1, 0) = real(i, dp)*b10*g(i - 1, 0) + c00*g(i, 0)" if ni >= 2 else ""}
{"      end do" if ni >= 2 else ""}
{"      g(0, 1) = cp00*g00" if nk >= 1 else ""}
{f"      do i = 1, {ni}" if nk >= 1 and ni >= 1 else ""}
{"         g(i, 1) = real(i, dp)*b00*g(i - 1, 0) + cp00*g(i, 0)" if nk >= 1 and ni >= 1 else ""}
{"      end do" if nk >= 1 and ni >= 1 else ""}
{f"      do k = 1, {nk - 1}" if nk >= 2 else ""}
{"         g(0, k + 1) = real(k, dp)*bp01*g(0, k - 1) + cp00*g(0, k)" if nk >= 2 else ""}
{f"         do i = 1, {ni}" if nk >= 2 and ni >= 1 else ""}
{"            g(i, k + 1) = real(i, dp)*b00*g(i - 1, k) &" if nk >= 2 and ni >= 1 else ""}
{"                        + real(k, dp)*bp01*g(i, k - 1) + cp00*g(i, k)" if nk >= 2 and ni >= 1 else ""}
{"         end do" if nk >= 2 and ni >= 1 else ""}
{"      end do" if nk >= 2 else ""}
   end subroutine rys2d_{tag}

   pure subroutine rys_hrr_{tag}(g, ab, cd, h)
      !! Horizontal transfer on both indices: (a,b+1| = (a+1,b| + (A-B)(a,b|.
      !! Binomials for a shift of 0, 1 or 2, written out -- LMAX stays 2.
      !$acc routine seq
      real(dp), intent(in)  :: g(0:{ni}, 0:{nk}), ab, cd
      real(dp), intent(out) :: h(0:{la}, 0:{lb}, 0:{lc}, 0:{ld})
      real(dp) :: t(0:{la}, 0:{lb}, 0:{nk})
      integer  :: a, b, c
      real(dp) :: ab2, cd2

      ab2 = ab*ab
      cd2 = cd*cd
      do a = 0, {la}
         do c = 0, {nk}
{chr(10).join(f"            t(a, {b}, c) = " + " + ".join(
    (("" if co_ == 1 else f"{float(co_)}_dp*") +
     ("" if p == 0 else (f"ab*" if p == 1 else f"ab2*")) + f"g(a + {t}, c)")
    for t, (co_, p) in enumerate((_BINOM[b][t], b - t) for t in range(b + 1)))
    for b in range(lb + 1))}
         end do
      end do
      do a = 0, {la}
         do b = 0, {lb}
            do c = 0, {lc}
{chr(10).join(f"               h(a, b, c, {d}) = " + " + ".join(
    (("" if co_ == 1 else f"{float(co_)}_dp*") +
     ("" if p == 0 else (f"cd*" if p == 1 else f"cd2*")) + f"t(a, b, c + {t})")
    for t, (co_, p) in enumerate((_BINOM[d][t], d - t) for t in range(d + 1)))
    for d in range(ld + 1))}
            end do
         end do
      end do
   end subroutine rys_hrr_{tag}
"""

    # `ren` renamed the module's own "contains" line nowhere -- it is bare
    # text -- so insert the component tables and the two helpers at the
    # module's original `contains` (the one right after its TWO_PI_2_5
    # parameter, in the specification part, before any subroutine; the
    # tables are `parameter` arrays and may not follow one).
    assert head.count("\ncontains\n") == 1
    head = head.replace("\ncontains\n", helpers, 1)

    out = head + middle + decode + tail + f"\nend module trc_pc_rys{tag}\n"
    dest = args.output
    with open(dest, "w") as fh:
        fh.write(out)
    print(f"  {dest}: ({la}{lb}|{lc}{ld}), {nroots} roots, vbuf({nv}), {len(out.splitlines())} lines")


if __name__ == "__main__":
    main()
