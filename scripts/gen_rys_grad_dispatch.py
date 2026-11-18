#!/usr/bin/env python3
"""Emit the dispatcher over the Rys per-class GRADIENT kernels
`gen_rys_grad.py` wrote. Same role as gen_rys_dispatch.py, for the gradient
driver (`trc_grad.F90`'s `grad_bins`) instead of `trc_rys_select.F90`: there
is no HGP fallback for the gradient, so `grad_bins` runs the screen and
compact dispatchers unconditionally for every class run and relies on every
class up to RYSG_LMAX that the pair bins launch having a kernel --
`rysg_has_class` exists anyway, as a cheap defensive check against a class
CMake's loop did not actually cover. RYSG_LMAX (the highest l in the class
list) is what the C API's gradient refusal reads, so a basis above it is
refused before the SCF rather than stopped inside `grad_bins`.
"""
import argparse
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gp = _load("gp", os.path.join(HERE, "gen_perclass.py"))
gg = _load("gg", os.path.join(HERE, "gen_rys_grad.py"))
PAIRS_NCB = gg.PAIRS_NCB

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--classes", required=True,
                    help="comma-separated four-digit class tags")
    # The options that decide a class's pair widths, from the generator's own
    # parser: pass this the same ones the generator and its --pair-widths
    # query got, or the widths this names are not the ones that were written.
    gg.add_plan_options(ap)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    tags = args.classes.split(",")
    widths = {tag: gg.class_pair_widths(tag, args) for tag in tags}

    keys = []
    use_lines = []
    lmax = max(int(c) for tag in tags for c in tag)
    for tag in tags:
        la, lb, lc, ld = (int(c) for c in tag)
        key = ((la*gp.CLASS_RADIX + lb)*gp.CLASS_RADIX + lc)*gp.CLASS_RADIX + ld
        keys.append((key, tag))
        use_lines.append(f"   use trc_pc_rysg{tag}, only: rysg{tag}_screen, rysg{tag}_compact")
        use_lines.extend(f"   use {gg.pair_module_name(tag, nc)}, only: {gg.pair_driver_name(tag, nc)}"
                         for nc in widths[tag])
    keys.sort()

    # Wrapped, not one line: at TERCO_LMAX=3 the list runs to 242 keys
    # (~1400 characters), which nvfortran rejects as a source line too long.
    case_list = ", &\n               ".join(
        ", ".join(str(k) for k, _ in keys[i:i + 16])
        for i in range(0, len(keys), 16))
    screen_cases = "\n".join(
        f"""      case ({key}); call rysg{tag}_screen(lo, hi, w0, w1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                          npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, sh_at, &
                          out_seg, out_iab, out_icd, nmax, ctr)"""
        for key, tag in keys)

    compact_cases = "\n".join(
        f"""      case ({key}); call rysg{tag}_compact(nout, out_seg, out_iab, out_icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &
                          npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, omega, nbas, nao, ao_off, sh_at, natm, &
                          nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, ndens, dmat, npriv, gp)"""
        for key, tag in keys)

    def pairs_case(key, tag):
        """One class's case: walk the columns in blocks of the widest width
        the class has, and run each block through the narrowest width that
        holds it (the columns past a block of a wider width would only be
        skipped). The block is the table's own sub-table: pstart(c0) and
        gp(.., c0) start the dummies there."""
        ws = widths[tag]
        if not ws:
            return (f"      case ({key})\n"
                    f'         error stop "terco: gradient_pairs: the class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}) '
                    'is never launched"')
        vmax = ws[-1]
        txt = (f"      case ({key})\n"
               f"         do c0 = 1, ncb, {vmax}\n"
               f"            nblk = min({vmax}, ncb - c0 + 1)\n")
        for n, nc in enumerate(ws):
            kw = ("if (nblk <= %d) then" % nc if n == 0 else
                  "else if (nblk <= %d) then" % nc if nc != vmax else "else")
            txt += (f"            {kw}\n"
                    f"               call {gg.pair_driver_name(tag, nc)}(nout, out_seg, out_iab, out_icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &\n"
                    "                          npair, sp_i, sp_j, sp_q, pcut, omega, nbas, nao, ao_off, sh_at, natm, &\n"
                    "                          nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &\n"
                    "                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &\n"
                    "                          npool, dmat, nblk, pstart(c0), pl, pr, pj, pk, npriv, gp(1, 1, 1, c0))\n")
        return txt + "            end if\n         end do"

    pairs_cases = "\n".join(pairs_case(key, tag) for key, tag in keys)

    out = f"""!
! Dispatcher over the Rys per-class GRADIENT kernels gen_rys_grad.py wrote.
!
! GENERATED by scripts/gen_rys_grad_dispatch.py --classes {args.classes} -- do not edit.
!
! The ORDERED SURVIVOR COMPACTION pair `rysg_screen_dispatch` /
! `rysg_compact_dispatch` that `grad_bins` runs per class -- see
! scripts/gen_rys_grad.py's screen_and_compact_text for the technique.
!
module trc_rys_grad_kernels
   use pic_types, only: dp, int64
   use trc_tables, only: LMAX
{chr(10).join(use_lines)}
   implicit none
   private
   public :: rysg_has_class
   public :: rysg_screen_dispatch, rysg_compact_dispatch, rysg_pairs_dispatch

   integer, parameter, public :: RYSG_PAIRS_NCB = {PAIRS_NCB}
      !! Columns one rysg_pairs_dispatch call carries (gen_rys_grad.py's PAIRS_NCB).

   integer, parameter, public :: RYSG_LMAX = {lmax}
      !! Highest l with gradient kernels: above it the analytic gradient is refused.

contains

   pure logical function rysg_has_class(key)
      !! Is there a generated Rys gradient kernel for this class key? There
      !! is no HGP gradient fallback, so a class `grad_bins` meets without
      !! one is a configuration error, not a silently skipped contribution --
      !! callers that care check this before dispatching.
      integer, intent(in) :: key
      select case (key)
      case ({case_list})
         rysg_has_class = .true.
      case default
         rysg_has_class = .false.
      end select
   end function rysg_has_class

   subroutine rysg_screen_dispatch(key, lo, hi, w0, w1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                          npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, sh_at, &
                          out_seg, out_iab, out_icd, nmax, ctr)
      !! ORDERED SURVIVOR COMPACTION, screen pass: decode + the three
      !! quartet-level screens only, over work items [w0, w1] of segments
      !! [lo, hi], appending survivors to a dense list. See
      !! scripts/gen_rys_grad.py's screen_and_compact_text.
      integer, intent(in) :: key
      integer,  intent(in)    :: lo, hi, nseg, npair, nbas, nmax
      integer(int64), intent(in) :: w0, w1
      integer(int64), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_at(nbas)
      integer,  intent(inout) :: out_seg(nmax), out_iab(nmax), out_icd(nmax)
      integer(int64), intent(inout) :: ctr(1)
      select case (key)
{screen_cases}
      end select
   end subroutine rysg_screen_dispatch

   subroutine rysg_compact_dispatch(key, nout, out_seg, out_iab, out_icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &
                          npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, omega, nbas, nao, ao_off, sh_at, natm, &
                          nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                          ndens, dmat, npriv, gp)
      !! ORDERED SURVIVOR COMPACTION, compact pass: the integral body, run
      !! over the dense survivor
      !! list `rysg_screen_dispatch` produced, adding into npriv
      !! lane-private gradient copies gp (see gen_rys_grad.py's compact
      !! driver) that the caller sums into grad. `omega` switches the
      !! operator to erf(omega*r12)/r12, same meaning as eri_fock's own --
      !! the screen pass above takes no `omega`, since the Schwarz bound it
      !! tests against stays the plain Coulomb one (see eri_fock's own
      !! comment: pointwise larger, so it bounds the attenuated integral
      !! too).
      integer, intent(in) :: key
      integer(int64), intent(in) :: nout
      integer,  intent(in)    :: out_seg(*), out_iab(*), out_icd(*)
      integer,  intent(in)    :: nseg, npair, nbas, nao
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      integer,  intent(in)    :: sSA(nseg), sSB(nseg)
      !! Pairs in the whole bra / ket BIN: the stride of the bin-interleaved
      !! primitive tables, apart from the segment's own extent sA/sNB. They
      !! differ exactly when a segment covers part of a bin, which every
      !! segment of a skeleton gradient plan does.
      logical,  intent(in)    :: sD(nseg)
      real(dp), intent(in)    :: sW(nseg)
      !! The segment's orbit weight (one unless the plan is a skeleton one),
      !! which the kernels fold into their degeneracy factor.
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), pcut, jfac, kfac, omega
      integer,  intent(in)    :: ao_off(nbas), sh_at(nbas), natm
      integer,  intent(in)    :: nqp
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp), qp_a(nqp)
      real(dp), intent(in)    :: qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(npair), sq_ay(npair), sq_az(npair)
      real(dp), intent(in)    :: sq_abx(npair), sq_aby(npair), sq_abz(npair)
      integer,  intent(in)    :: sq_n(npair)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      integer,  intent(in)    :: npriv
      real(dp), intent(inout) :: gp(3, natm, npriv)
      select case (key)
{compact_cases}
      end select
   end subroutine rysg_compact_dispatch

   subroutine rysg_pairs_dispatch(key, nout, out_seg, out_iab, out_icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &
                          npair, sp_i, sp_j, sp_q, pcut, omega, nbas, nao, ao_off, sh_at, natm, &
                          nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                          npool, dmat, ncb, pstart, pl, pr, pj, pk, npriv, gp)
      !! The compact pass over a COLUMN PAIR TABLE: every survivor's derivative
      !! integrals are contracted against all `ncb` columns (at most
      !! RYSG_PAIRS_NCB) at once, into gp(:, :, :, column). Column c owns pairs
      !! pstart(c) .. pstart(c+1)-1 of (pl, pr, pj, pk); see
      !! gen_rys_grad.py's pairs_modules for what a pair contributes.
      integer, intent(in) :: key
      integer(int64), intent(in) :: nout
      integer,  intent(in)    :: out_seg(*), out_iab(*), out_icd(*)
      integer,  intent(in)    :: nseg, npair, nbas, nao
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      integer,  intent(in)    :: sSA(nseg), sSB(nseg)
      logical,  intent(in)    :: sD(nseg)
      real(dp), intent(in)    :: sW(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), pcut, omega
      integer,  intent(in)    :: ao_off(nbas), sh_at(nbas), natm
      integer,  intent(in)    :: nqp
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp), qp_a(nqp)
      real(dp), intent(in)    :: qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(npair), sq_ay(npair), sq_az(npair)
      real(dp), intent(in)    :: sq_abx(npair), sq_aby(npair), sq_abz(npair)
      integer,  intent(in)    :: sq_n(npair)
      integer,  intent(in)    :: npool, ncb
      real(dp), intent(in)    :: dmat(npool, nao, nao)
      integer,  intent(in)    :: pstart(ncb + 1), pl(*), pr(*)
      real(dp), intent(in)    :: pj(*), pk(*)
      integer,  intent(in)    :: npriv
      real(dp), intent(inout) :: gp(3, natm, npriv, ncb)
      integer :: c0, nblk
      select case (key)
{pairs_cases}
      end select
   end subroutine rysg_pairs_dispatch

end module trc_rys_grad_kernels
"""
    with open(args.output, "w") as fh:
        fh.write(out)
    print(f"  {args.output}: dispatcher for {len(keys)} Rys gradient classes")


if __name__ == "__main__":
    main()
