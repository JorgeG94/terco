#!/usr/bin/env python3
"""Emit a Rys GRADIENT kernel for one class, by splicing terco's own scalar
HGP kernel -- the same technique gen_rys.py uses for the energy kernel, with
the same reasons (see that file's docstring for why splicing rather than
generating whole).

WHAT IS DIFFERENT FROM THE ENERGY SPLICE
-----------------------------------------
The energy kernel defers the density fold to AFTER the primitive/root loop:
it sums the bare integral into `vbuf` and only then multiplies by the density
blocks, once per quartet. The gradient cannot do that -- the derivative
weight `fix`/`fjx`/`fkx`/`flx` depends on the PRIMITIVE's own exponent, so it
has to be formed, and folded against the density, INSIDE the primitive/root
loop. That means the density blocks have to be loaded BEFORE that loop
starts rather than after it ends, so the HEAD/DECODE/TAIL split gen_rys.py
uses does not carry over unchanged:

  - HEAD is still the segment decode, screens, bin-interleaved primitive
    loads and primitive-quartet prescreen, verbatim -- but with the density
    blocks, the degeneracy weight `wq`, and a precomputed `gam(idx)` per
    output Cartesian component folded in right where the HGP file already
    has every piece they need (si/sj/sk/sl, same_ab/same_cd/same_pair, the
    AO offsets) and right before the primitive loop starts. This is new
    text, not copied from HGP's tail -- the energy digestion's "Single
    column..." decode and `#ifdef TRC_NO_DIGEST` block are not reused at
    all.
  - MIDDLE replaces the Rys integral build exactly like the energy splice,
    except it builds ONE extended horizontal-transfer table per axis,
    `h(0:la+1, 0:lb, 0:lc+1, 0:ld)` -- i and k (the bra- and ket-combined
    totals) raised by one each, j and l never raised -- and folds the
    derivative directly into 12 running sums (A, B, C, D x, y, z) instead of
    an output integral block.
  - TAIL is new too: four atomic updates into `grad(3, natm)`, through the
    shell -> atom map `sh_at`, in place of HGP's six-block Fock digestion.

THE ALGEBRA (GPU4PySCF's, not re-derived)
------------------------------------------
Checked against gpu4pyscf's actual production kernel
(gvhf-rys/rys_contract_jk_ip1.cu, the BUILD_4C_GXYZ/fix/fjx/fkx/flx block),
not only its description:

  - Raise ONLY the bra- and ket-combined 2D-integral indices (i, k), never
    j or l. The j derivative comes from the horizontal-transfer identity
    instead of a third raised table:
        fix = 2*alpha_i*h(i+1,j,k,l) [- i*h(i-1,j,k,l) if i>0]
        fkx = 2*alpha_k*h(i,j,k+1,l) [- k*h(i,j,k-1,l) if k>0]
        fjx = 2*alpha_j*(h(i+1,j,k,l) + (Ax-Bx)*h(i,j,k,l)) [- j*h(i,j-1,k,l) if j>0]
        flx = -fix - fjx - fkx                      (translational invariance)
  - nroots = (L + 1)//2 + 1 with L = la+lb+lc+ld (integer division), NOT
    (L+2)//2+1: every derivative TERM raises only one of i or k, never both
    in the same product, so the polynomial degree actually summed is L+1,
    not L+2, even though the table's nominal range covers both raised by
    one. Confirmed against gpu4pyscf's own bounds (gvhf-rys/...ip1.cu
    around its `order = li+lj+lk+ll; nroots = (order+1)/2+1`). This means
    the existing n = 1..5 root tables cover every class up to (dd|dd); no
    new table is needed.
  - Gamma, from the standard RHF/UHF two-electron gradient identity (not
    copied from the Fock digestion's internal 4x/divide-by-8 bookkeeping,
    which is a Fock-MATRIX convention that does not transfer to a gradient
    SCALAR):
        RHF (ndens=1, dmat(:,:,1) = total D):
            Gamma = jfac*D_ab*D_cd - kfac/4*(D_ac*D_bd + D_ad*D_bc)
        UHF (ndens=2, dmat(:,:,1)=Da, dmat(:,:,2)=Db):
            Gamma = jfac*(Da+Db)_ab*(Da+Db)_cd
                    - kfac*(Da_ac*Da_bd+Da_ad*Da_bc + Db_ac*Db_bd+Db_ad*Db_bc)
    jfac/kfac default to 1 (pure HF); a hybrid functional's exchange
    fraction multiplies kfac, same meaning as eri_fock's k_scale.
  - A quartet whose four shells all sit on the same atom contributes
    nothing (A=B=C=D already satisfies translational invariance on its
    own) and is skipped before the primitive loop.

Never hand-edit the output; change this generator and regenerate.
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
#: needs. Index [b][t] is C(b, t). Same table as gen_rys.py's.
_BINOM = [[1], [1, 1], [1, 2, 1]]


def hrr_sub(name, ni, nk, la, lb, lc, ld, ab, cd):
    """Text for `pure subroutine {name}(g, {ab}, {cd}, h)`: the horizontal
    transfer from a combined-index 2D table g(0:ni, 0:nk) to
    h(0:la, 0:lb, 0:lc, 0:ld). Same algebra as gen_rys.py's rys_hrr_{tag},
    factored out so the gradient kernel can call it with la+1/lc+1/ni+1/nk+1
    for the extended table without duplicating the derivation."""
    t_rows = []
    for b in range(lb + 1):
        terms = []
        for t in range(b + 1):
            coef, p = _BINOM[b][t], b - t
            pre = "" if coef == 1 else f"{float(coef)}_dp*"
            mul = "" if p == 0 else (f"{ab}*" if p == 1 else f"{ab}2*")
            terms.append(f"{pre}{mul}g(a + {t}, c)")
        t_rows.append(f"            t(a, {b}, c) = " + " + ".join(terms))
    h_rows = []
    for d in range(ld + 1):
        terms = []
        for t in range(d + 1):
            coef, p = _BINOM[d][t], d - t
            pre = "" if coef == 1 else f"{float(coef)}_dp*"
            mul = "" if p == 0 else (f"{cd}*" if p == 1 else f"{cd}2*")
            terms.append(f"{pre}{mul}t(a, b, c + {t})")
        h_rows.append(f"               h(a, b, c, {d}) = " + " + ".join(terms))
    return f"""
   pure subroutine {name}(g, {ab}, {cd}, h)
      !! Horizontal transfer on both indices, to the gradient-extended
      !! table: (a,b+1| = (a+1,b| + (A-B)(a,b|. Binomials for a shift of
      !! 0, 1 or 2, written out -- LMAX stays 2.
      !$acc routine seq
      real(dp), intent(in)  :: g(0:{ni}, 0:{nk}), {ab}, {cd}
      real(dp), intent(out) :: h(0:{la}, 0:{lb}, 0:{lc}, 0:{ld})
      real(dp) :: t(0:{la}, 0:{lb}, 0:{nk})
      integer  :: a, b, c
      real(dp) :: {ab}2, {cd}2

      {ab}2 = {ab}*{ab}
      {cd}2 = {cd}*{cd}
      do a = 0, {la}
         do c = 0, {nk}
{chr(10).join(t_rows)}
         end do
      end do
      do a = 0, {la}
         do b = 0, {lb}
            do c = 0, {lc}
{chr(10).join(h_rows)}
            end do
         end do
      end do
   end subroutine {name}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("class_tag", help="four digits, e.g. 1111 or 2222")
    ap.add_argument("--order", choices=("lexicographic", "gamess"),
                    default="lexicographic",
                    help="must match the --order the scalar HGP file in "
                         "--src-dir was generated with")
    ap.add_argument("--src-dir", default="src/generated",
                    help="where trc_pc_k{tag}.F90 (the scalar HGP kernel) "
                         "already lives")
    ap.add_argument("-o", "--output", required=True,
                    help="never a path under src/generated: regenerated at "
                         "configure time into the build directory")
    args = ap.parse_args()
    tag = args.class_tag
    la, lb, lc, ld = (int(c) for c in tag)
    ni, nk = la + lb, lc + ld
    nig, nkg = ni + 1, nk + 1       # the gradient-extended table's dims
    L = la + lb + lc + ld
    nroots = (L + 1)//2 + 1
    if nroots > 5:
        raise SystemExit(f"class {tag}: nroots={nroots} exceeds the 5 "
                          "committed root tables")
    dims = [len(co.components(l, args.order)) for l in (la, lb, lc, ld)]
    nv = dims[0]*dims[1]*dims[2]*dims[3]

    src_path = os.path.join(args.src_dir, f"trc_pc_k{tag}.F90")
    src = open(src_path).read()

    i_item = src.index(f"   pure subroutine pcsi{tag}(")
    i_end = src.index(f"end subroutine pcsi{tag}\n") + len(f"end subroutine pcsi{tag}\n")
    whole = src[:i_end]
    assert i_item < i_end

    tval_line = "tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)\n"
    i_tval = whole.index(tval_line) + len(tval_line)
    head = whole[:i_tval]

    # --- 1. thread qp_a (the bra/ket primitive's own exponent) through the
    #     signatures and the one call between them, exactly like qp_z. ---
    head, n1 = re.subn(r"qp_z\(nqp\), qp_c\(nqp\)",
                       "qp_z(nqp), qp_c(nqp), qp_a(nqp)", head)
    assert n1 == 2, n1
    head, n1 = re.subn(r"nqp, qp_z, qp_c, qp_px, qp_py, qp_pz",
                       "nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz", head)
    assert n1 == 3, n1

    # --- 2. thread sh_at/natm/grad in place of jmat, next to sh_l/ao_off. ---
    head, n1 = re.subn(r"sh_l, ao_off,", "sh_l, ao_off, sh_at, natm,", head)
    assert n1 == 3, n1
    head, n1 = re.subn(r"sh_l\(nbas\), ao_off\(nbas\)",
                       "sh_l(nbas), ao_off(nbas)\n"
                       "      integer,  intent(in)    :: sh_at(nbas), natm",
                       head)
    assert n1 == 2, n1
    head, n1 = re.subn(r"dmat, jmat", "dmat, grad", head)
    assert n1 == 3, n1
    head, n1 = re.subn(r"real\(dp\), intent\(inout\) :: jmat\(ndens, nao, nao\)",
                       "real(dp), intent(inout) :: grad(3, natm)", head)
    assert n1 == 2, n1

    # --- 3. alpha_i/alpha_k: the individual primitive exponent, out of the
    #     SAME bin-interleaved slot qp_z/qp_px read -- not a new primitive
    #     loop, just one more field of the pair already being visited. ---
    head, n1 = re.subn(
        r"(rpx = qp_px\(kp\); rpy = qp_py\(kp\); rpz = qp_pz\(kp\)\n)",
        r"\1               alpha_i = qp_a(kp)\n"
        r"               alpha_j = zeta - alpha_i\n", head)
    assert n1 == 1, n1
    head, n1 = re.subn(
        r"(eta = qp_z\(kq\)\n)",
        r"\1                  alpha_k = qp_a(kq)\n", head)
    assert n1 == 1, n1

    # --- 4. drop the pre-HRR VRR declaration and the Boys scratch, exactly
    #     like gen_rys.py -- the gradient needs neither. ---
    decl_re = re.compile(
        r"      real\(dp\) :: v\((\d+), 0:1\), g1\(\1\), vbuf\((\d+)\)\n")
    m = decl_re.search(head)
    assert m, "could not find the HGP v/g1/vbuf declaration to replace"
    nv_vrr = int(m.group(1))
    nv_in_file = int(m.group(2))
    assert nv_in_file == nv, (
        f"vbuf({nv_in_file}) disagrees with na*nb*nc*nd = {nv}; "
        "the --order argument likely does not match the source file")

    decl_new = (
        f"      real(dp) :: ru({nroots}), rw({nroots})\n"
        f"      real(dp) :: gx(0:{nig}, 0:{nkg}), gy(0:{nig}, 0:{nkg}), gz(0:{nig}, 0:{nkg})\n"
        f"      real(dp) :: hx(0:{la + 1}, 0:{lb}, 0:{lc + 1}, 0:{ld})\n"
        f"      real(dp) :: hy(0:{la + 1}, 0:{lb}, 0:{lc + 1}, 0:{ld})\n"
        f"      real(dp) :: hz(0:{la + 1}, 0:{lb}, 0:{lc + 1}, 0:{ld})\n"
        f"      real(dp) :: u2, dinv, b00r, b10r, bp01r, f00\n"
        f"      real(dp) :: c00x, c00y, c00z, cp0x, cp0y, cp0z\n"
        f"      real(dp) :: alpha_i, alpha_j, alpha_k, kfac_eff\n"
        f"      integer  :: rr, aa, bb, cc, dd2\n"
        f"      integer  :: eix, eiy, eiz, ejx, ejy, ejz\n"
        f"      integer  :: ekx, eky, ekz, elx, ely, elz\n"
        f"      real(dp) :: vix, viy, viz, pyz, pxz, pxy\n"
        f"      real(dp) :: fix, fiy, fiz, fjx, fjy, fjz\n"
        f"      real(dp) :: fkx, fky, fkz, flx, fly, flz\n"
        f"      real(dp) :: vax, vay, vaz, vbx, vby, vbz\n"
        f"      real(dp) :: vcx, vcy, vcz, vdx, vdy, vdz\n"
        f"      real(dp) :: dab({dims[0]*dims[1]}), dcd({dims[2]*dims[3]})\n"
        f"      real(dp) :: dac(2, {dims[0]*dims[2]}), dad(2, {dims[0]*dims[3]})\n"
        f"      real(dp) :: dbc(2, {dims[1]*dims[2]}), dbd(2, {dims[1]*dims[3]})\n"
        f"      real(dp) :: gam({nv})\n")
    head = head[:m.start()] + decl_new + head[m.end():]

    for pat in (r"      real\(dp\) :: f\(0:BOYS_MMAX\)\n",
                r"      integer  :: bi, bj, bbase\n",
                r"      real\(dp\) :: bx, bx2, b0, b1, b2, btt, bet\n"):
        head, n1 = re.subn(pat, "", head, count=1)
        assert n1 == 1, pat

    # The six Fock-digestion block declarations (jab/jcd/kac/kad/kbc/kbd and
    # dab/dcd/dac/dad/dbc/dbd) are dead for the gradient -- `dab`/`dcd` are
    # redeclared above with the shape this file needs, and the rest are
    # never read at all once `tail` is replaced. Drop all four lines.
    for pat in (r"      real\(dp\) :: jab\(\d+\), jcd\(\d+\), kac\(\d+\)\n",
                r"      real\(dp\) :: kad\(\d+\), kbc\(\d+\), kbd\(\d+\)\n",
                r"      real\(dp\) :: dab\(\d+\), dcd\(\d+\), dac\(\d+\)\n",
                r"      real\(dp\) :: dad\(\d+\), dbc\(\d+\), dbd\(\d+\)\n"):
        head, n1 = re.subn(pat, "", head, count=1)
        assert n1 == 1, pat

    # No Boys function; the root lookup replaces it, same as the energy splice.
    head, n1 = re.subn(r"   use trc_boys, only:.*?BOYS_DT, BOYS_DTINV\n",
                       f"   use trc_boys, only: dp\n"
                       f"   use trc_rys, only: trc_rys_roots{nroots}_tab\n",
                       head, count=1, flags=re.S)
    assert n1 == 1

    # --- 5. the zero-the-VRR-accumulator loop becomes: same-atom skip, the
    #     degeneracy weight, the AO offsets, the density blocks and the
    #     per-component Gamma, all computed ONCE per quartet before the
    #     primitive loop -- see the module docstring for why this has to
    #     move here instead of staying in a deferred digestion. ---
    zero_re = re.compile(r"            do x = 1, " + str(nv_vrr) +
                         r"\n               g1\(x\) = 0\.0_dp\n            end do\n")
    m2 = zero_re.search(head)
    assert m2, "could not find the g1 zero loop to replace"

    d0, d1, d2, d3 = dims
    inject = f"""
         ! A quartet whose four shells sit on one atom is translation-
         ! invariant on its own: dGrad/dR = 0 for every R trivially, so skip
         ! it before the primitive loop rather than compute and cancel it.
         if (sh_at(si) == sh_at(sj) .and. sh_at(si) == sh_at(sk) &
             .and. sh_at(si) == sh_at(sl)) return

         mui = ao_off(si); nuj = ao_off(sj); lamk = ao_off(sk); sigl = ao_off(sl)
         wq = 1.0_dp
         if (same_ab) wq = wq*0.5_dp
         if (same_cd) wq = wq*0.5_dp
         if (same_pair) wq = wq*0.5_dp
         ! The 8 below is NOT part of wq: terco visits one representative
         ! quartet per (ab|cd)/(cd|ab)/(ba|cd)/... 8-fold permutation class
         ! (wq already halves it where two of those labels coincide), but
         ! dE2/dR = sum Gamma (ab|cd)/dR runs over the FULL unrestricted
         ! index sum, where a generic (non-degenerate) quartet's integral
         ! appears in all 8 equivalent terms. Found empirically against the
         ! finite-difference oracle (check_rys_grad_fd.F90): without it the
         ! analytic gradient was a uniform 1/8 of the correct one on every
         ! component, for every system and direction tried.
         !
         ! kfac_eff: RHF's Gamma = D_ab D_cd - kfac/4 (D_ac D_bd + D_ad D_bc)
         ! using the TOTAL density: the 1/4 is standard (Szabo & Ostlund),
         ! checked against the FD oracle directly. UHF's K term is NOT the
         ! same coefficient reused per spin -- it is 1/2, not 1/4 or 1 --
         ! found the same way the 8 above was: Da = Db is the RHF solution
         ! with D = 2 Da, so eri%gradient(Da, Da) must equal eri%gradient
         ! (2 Da) exactly; only kfac_eff = 1/2 per spin reproduces that
         ! identity (1/4 and 1 both failed it before the FD oracle confirmed
         ! the full density was actually wrong, not just this constant).
         kfac_eff = merge(0.25_dp, 0.5_dp, ndens == 1)*kfac

         do ib = 0, {d1 - 1}
            do ia = 0, {d0 - 1}
               dab(1 + ia + {d0}*ib) = 0.0_dp
               do idens = 1, ndens
                  dab(1 + ia + {d0}*ib) = dab(1 + ia + {d0}*ib) + dmat(idens, mui + ia, nuj + ib)
               end do
            end do
         end do
         do id = 0, {d3 - 1}
            do ic = 0, {d2 - 1}
               dcd(1 + ic + {d2}*id) = 0.0_dp
               do idens = 1, ndens
                  dcd(1 + ic + {d2}*id) = dcd(1 + ic + {d2}*id) + dmat(idens, lamk + ic, sigl + id)
               end do
            end do
         end do
         do idens = 1, ndens
            do ic = 0, {d2 - 1}
               do ia = 0, {d0 - 1}
                  dac(idens, 1 + ia + {d0}*ic) = dmat(idens, mui + ia, lamk + ic)
               end do
            end do
            do id = 0, {d3 - 1}
               do ia = 0, {d0 - 1}
                  dad(idens, 1 + ia + {d0}*id) = dmat(idens, mui + ia, sigl + id)
               end do
            end do
            do ic = 0, {d2 - 1}
               do ib = 0, {d1 - 1}
                  dbc(idens, 1 + ib + {d1}*ic) = dmat(idens, nuj + ib, lamk + ic)
               end do
            end do
            do id = 0, {d3 - 1}
               do ib = 0, {d1 - 1}
                  dbd(idens, 1 + ib + {d1}*id) = dmat(idens, nuj + ib, sigl + id)
               end do
            end do
         end do

         idx = 0
         do dd2 = 1, {d3}
            do cc = 1, {d2}
               do bb = 1, {d1}
                  do aa = 1, {d0}
                     idx = idx + 1
                     gam(idx) = 8.0_dp*wq*jfac*dab(aa + {d0}*(bb - 1))*dcd(cc + {d2}*(dd2 - 1))
                     do idens = 1, ndens
                        gam(idx) = gam(idx) - 8.0_dp*wq*kfac_eff*( &
                           dac(idens, aa + {d0}*(cc - 1))*dbd(idens, bb + {d1}*(dd2 - 1)) &
                           + dad(idens, aa + {d0}*(dd2 - 1))*dbc(idens, bb + {d1}*(cc - 1)))
                     end do
                  end do
               end do
            end do
         end do

         vax = 0.0_dp; vay = 0.0_dp; vaz = 0.0_dp
         vbx = 0.0_dp; vby = 0.0_dp; vbz = 0.0_dp
         vcx = 0.0_dp; vcy = 0.0_dp; vcz = 0.0_dp
         vdx = 0.0_dp; vdy = 0.0_dp; vdz = 0.0_dp
"""
    head = head[:m2.start()] + inject + head[m2.end():]

    def ren(t):
        t = t.replace(f"trc_pc_k{tag}", f"trc_pc_rysg{tag}")
        t = t.replace(f"pcsi{tag}", f"rysgi{tag}")
        t = t.replace(f"pcs{tag}", f"rysg{tag}")
        return t

    head = re.sub(r"public :: pcs" + tag + r".*\n", f"public :: rysg{tag}\n", head, count=1)
    head = ren(head)

    old_header_re = re.compile(r"\A!\n.*?#endif\n", re.S)
    assert old_header_re.search(head)
    order_val = 0 if args.order == "lexicographic" else 1
    new_header = f"""!
! Class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}) GRADIENT kernel, Rys quadrature.
!
! GENERATED by scripts/gen_rys_grad.py {tag} -- do not edit.
!
! Spliced from the scalar HGP kernel in src/generated/trc_pc_k{tag}.F90: the
! segment decode, screens and primitive-pair loads are that file's text,
! verbatim. The density fold and the digestion are NOT HGP's -- see this
! generator's module docstring for why they could not be.
!
! TRC_CART_ORDER GUARD -- see gen_perclass.py's non-split header for the
! full explanation; same guard, carried over from the file this was spliced
! from.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
"""
    head = old_header_re.sub(new_header, head, count=1)

    # --- the Rys gradient integral: roots, the 2D recurrence per axis, the
    #     extended horizontal transfer, then the derivative fold straight
    #     into the 12 running sums -- no output integral block is ever
    #     materialised. ---
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
                     call rys_hrr_g_{tag}(gx, abx, cdx, hx)
                     call rys_hrr_g_{tag}(gy, aby, cdy, hy)
                     call rys_hrr_g_{tag}(gz, abz, cdz, hz)

                     idx = 0
                     do dd2 = 1, {d3}
                        do cc = 1, {d2}
                           do bb = 1, {d1}
                              do aa = 1, {d0}
                                 idx = idx + 1
                                 eix = cnx0_{tag}(aa); eiy = cny0_{tag}(aa); eiz = cnz0_{tag}(aa)
                                 ejx = cnx1_{tag}(bb); ejy = cny1_{tag}(bb); ejz = cnz1_{tag}(bb)
                                 ekx = cnx2_{tag}(cc); eky = cny2_{tag}(cc); ekz = cnz2_{tag}(cc)
                                 elx = cnx3_{tag}(dd2); ely = cny3_{tag}(dd2); elz = cnz3_{tag}(dd2)

                                 vix = hx(eix, ejx, ekx, elx); viy = hy(eiy, ejy, eky, ely)
                                 viz = hz(eiz, ejz, ekz, elz)
                                 pyz = viy*viz; pxz = vix*viz; pxy = vix*viy

                                 fix = 2.0_dp*alpha_i*hx(eix + 1, ejx, ekx, elx)
                                 if (eix > 0) fix = fix - real(eix, dp)*hx(eix - 1, ejx, ekx, elx)
                                 fkx = 2.0_dp*alpha_k*hx(eix, ejx, ekx + 1, elx)
                                 if (ekx > 0) fkx = fkx - real(ekx, dp)*hx(eix, ejx, ekx - 1, elx)
                                 fjx = 2.0_dp*alpha_j*(hx(eix + 1, ejx, ekx, elx) + abx*hx(eix, ejx, ekx, elx))
                                 if (ejx > 0) fjx = fjx - real(ejx, dp)*hx(eix, ejx - 1, ekx, elx)
                                 flx = -fix - fjx - fkx

                                 fiy = 2.0_dp*alpha_i*hy(eiy + 1, ejy, eky, ely)
                                 if (eiy > 0) fiy = fiy - real(eiy, dp)*hy(eiy - 1, ejy, eky, ely)
                                 fky = 2.0_dp*alpha_k*hy(eiy, ejy, eky + 1, ely)
                                 if (eky > 0) fky = fky - real(eky, dp)*hy(eiy, ejy, eky - 1, ely)
                                 fjy = 2.0_dp*alpha_j*(hy(eiy + 1, ejy, eky, ely) + aby*hy(eiy, ejy, eky, ely))
                                 if (ejy > 0) fjy = fjy - real(ejy, dp)*hy(eiy, ejy - 1, eky, ely)
                                 fly = -fiy - fjy - fky

                                 fiz = 2.0_dp*alpha_i*hz(eiz + 1, ejz, ekz, elz)
                                 if (eiz > 0) fiz = fiz - real(eiz, dp)*hz(eiz - 1, ejz, ekz, elz)
                                 fkz = 2.0_dp*alpha_k*hz(eiz, ejz, ekz + 1, elz)
                                 if (ekz > 0) fkz = fkz - real(ekz, dp)*hz(eiz, ejz, ekz - 1, elz)
                                 fjz = 2.0_dp*alpha_j*(hz(eiz + 1, ejz, ekz, elz) + abz*hz(eiz, ejz, ekz, elz))
                                 if (ejz > 0) fjz = fjz - real(ejz, dp)*hz(eiz, ejz - 1, ekz, elz)
                                 flz = -fiz - fjz - fkz

                                 vax = vax + fix*pyz*gam(idx); vay = vay + fiy*pxz*gam(idx); vaz = vaz + fiz*pxy*gam(idx)
                                 vbx = vbx + fjx*pyz*gam(idx); vby = vby + fjy*pxz*gam(idx); vbz = vbz + fjz*pxy*gam(idx)
                                 vcx = vcx + fkx*pyz*gam(idx); vcy = vcy + fky*pxz*gam(idx); vcz = vcz + fkz*pxy*gam(idx)
                                 vdx = vdx + flx*pyz*gam(idx); vdy = vdy + fly*pxz*gam(idx); vdz = vdz + flz*pxy*gam(idx)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do

         !$acc atomic update
         grad(1, sh_at(si)) = grad(1, sh_at(si)) + vax
         !$acc atomic update
         grad(2, sh_at(si)) = grad(2, sh_at(si)) + vay
         !$acc atomic update
         grad(3, sh_at(si)) = grad(3, sh_at(si)) + vaz
         !$acc atomic update
         grad(1, sh_at(sj)) = grad(1, sh_at(sj)) + vbx
         !$acc atomic update
         grad(2, sh_at(sj)) = grad(2, sh_at(sj)) + vby
         !$acc atomic update
         grad(3, sh_at(sj)) = grad(3, sh_at(sj)) + vbz
         !$acc atomic update
         grad(1, sh_at(sk)) = grad(1, sh_at(sk)) + vcx
         !$acc atomic update
         grad(2, sh_at(sk)) = grad(2, sh_at(sk)) + vcy
         !$acc atomic update
         grad(3, sh_at(sk)) = grad(3, sh_at(sk)) + vcz
         !$acc atomic update
         grad(1, sh_at(sl)) = grad(1, sh_at(sl)) + vdx
         !$acc atomic update
         grad(2, sh_at(sl)) = grad(2, sh_at(sl)) + vdy
         !$acc atomic update
         grad(3, sh_at(sl)) = grad(3, sh_at(sl)) + vdz
      end subroutine rysgi{tag}
"""

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
      !! One axis of 2D Rys integrals, gradient-extended by one unit of
      !! angular momentum on each of the bra- and ket-combined indices.
      !$acc routine seq
      real(dp), intent(in)  :: c00, cp00, b00, b10, bp01, g00
      real(dp), intent(out) :: g(0:{nig}, 0:{nkg})
      integer :: i, k

      g(0, 0) = g00
{"      g(1, 0) = c00*g00" if nig >= 1 else ""}
{f"      do i = 1, {nig - 1}" if nig >= 2 else ""}
{"         g(i + 1, 0) = real(i, dp)*b10*g(i - 1, 0) + c00*g(i, 0)" if nig >= 2 else ""}
{"      end do" if nig >= 2 else ""}
{"      g(0, 1) = cp00*g00" if nkg >= 1 else ""}
{f"      do i = 1, {nig}" if nkg >= 1 and nig >= 1 else ""}
{"         g(i, 1) = real(i, dp)*b00*g(i - 1, 0) + cp00*g(i, 0)" if nkg >= 1 and nig >= 1 else ""}
{"      end do" if nkg >= 1 and nig >= 1 else ""}
{f"      do k = 1, {nkg - 1}" if nkg >= 2 else ""}
{"         g(0, k + 1) = real(k, dp)*bp01*g(0, k - 1) + cp00*g(0, k)" if nkg >= 2 else ""}
{f"         do i = 1, {nig}" if nkg >= 2 and nig >= 1 else ""}
{"            g(i, k + 1) = real(i, dp)*b00*g(i - 1, k) &" if nkg >= 2 and nig >= 1 else ""}
{"                        + real(k, dp)*bp01*g(i, k - 1) + cp00*g(i, k)" if nkg >= 2 and nig >= 1 else ""}
{"         end do" if nkg >= 2 and nig >= 1 else ""}
{"      end do" if nkg >= 2 else ""}
   end subroutine rys2d_{tag}
{hrr_sub(f"rys_hrr_g_{tag}", nig, nkg, la + 1, lb, lc + 1, ld, "ab", "cd")}
"""

    assert head.count("\ncontains\n") == 1
    head = head.replace("\ncontains\n", helpers, 1)

    out = head + middle + f"\nend module trc_pc_rysg{tag}\n"
    with open(args.output, "w") as fh:
        fh.write(out)
    print(f"  {args.output}: ({la}{lb}|{lc}{ld}), {nroots} roots, "
          f"{len(out.splitlines())} lines")


if __name__ == "__main__":
    main()
