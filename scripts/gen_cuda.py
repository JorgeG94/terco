"""CUDA Fortran cooperative kernels for generally contracted shells.

Loaded by gen_perclass.py (--cuda) the way gen_vrr/gen_hrr are, and fed the
same per-class pieces: the sieved VRR body, the HRR body, and the text of the
`do concurrent` block kernel, from which the item decode/screen head and the
HRR-plus-digestion tail are sliced verbatim so the three kernels of a class
cannot drift apart.

WHY A THIRD KERNEL
------------------
A generally contracted shell (cc-pVDZ silicon: twelve s primitives under three
contractions) puts the same primitive quartet under nine column combinations
of an (ss|ss) and 81 of a (pp|pp). One thread per quartet with registers only
-- which is what `do concurrent` is -- cannot hold the contracted blocks for
all of them, so the blocked kernel repeats the primitive VRR per block and on
every class with p on both sides it is the segmented cost again. Reusing the
VRR needs threads that cooperate through shared memory. This is that kernel,
and it is the only one in terco that is not plain `do concurrent`.

SHAPE
-----
One thread block per primitive-shell quartet item of a general segment run,
PG_T threads. The block walks the bra primitive pairs; for each, the ket pairs
go in chunks of PG_T with one lane per ket pair computing the primitive VRR
block v(nv) in registers, exactly the scalar kernel's code. The ket-column
contraction tc(x, cd) = sum_kq c_c c_d v(x) is a reduction across the lanes,
done through shared memory in x-chunks of PG_XC; its nv*ncdc entries are
owned round-robin by the threads (tcown). After the last ket chunk each owner
folds the bra columns, gown(x, ab, cd) += c_a c_b tc(x, cd), so the factorised
cost of the blocked kernel is kept: nv*nccd FMAs per primitive quartet,
nv*ncab*nccd per bra pair, and the VRR once. The column grid is chunked at
run time to a fixed register budget (PG_NTC tc entries, PG_NG g entries per
thread); a class whose grid does not fit repeats its primitive loops per
chunk, never more often than the blocked kernel did.

Digestion is one thread per column combination, in batches of PG_GB: the
owners scatter their g into shared memory, the thread copies its column to a
local g1 and runs the generated HRR and six-block folded digestion unchanged,
with `!$acc atomic update` become `atomicadd`.

Every `syncthreads` sits outside every lane-divergent branch: the loops over
ket chunks, x-chunks and digestion batches have block-uniform trip counts, and
the early returns of the decode/screen head are taken by all threads alike.
"""
import re

PG_T = 64      # threads per block: one ket primitive pair per lane
PG_XC = 32     # VRR components staged through shared memory at a time
PG_NTC = 8     # ket-contracted entries a thread owns, at most
PG_NG = 32     # bra-folded entries a thread owns, at most
PG_NCC = 16    # ket column pairs per chunk, at most (PS_NCOL_MAX**2)
PG_NCAB = 16   # bra column pairs per chunk, at most
SHM_G = 24576  # bytes of shared memory for the digestion gather

CUDA_PROLOGUE = """      integer, value :: nranks, lo, hi, nseg, npair, nbas, npp, nao, ncoltot, ncoef, ndens
      integer(kind=8), value :: g0, gend
      real(dp), value :: thresh, jfac, kfac
      integer(kind=8), device :: sOff(nseg + 1)
      integer,  device :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      logical,  device :: sD(nseg)
      integer,  device :: sp_i(npair), sp_j(npair)
      real(dp), device :: sp_q(npair)
      real(dp), device :: dsh(nbas, nbas)
      integer,  device :: sh_l(nbas), ao_off(nbas)
      integer,  device :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), device :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3)
      real(dp), device :: pp_rb(npp, 3), pp_c(npp)
      real(dp), device :: pp_cs(npp)
      integer,  device :: pp_ki(npp), pp_kj(npp)
      integer,  device :: ps_np(nbas), ps_ncol(nbas), ps_soff(nbas), ps_coff(nbas)
      integer,  device :: col_ao(ncoltot)
      real(dp), device :: ps_coef(ncoef)
      real(dp), device :: boys_d(*)
      real(dp), device :: dmat(ndens, nao, nao)
      real(dp), device :: jmat(ndens, nao, nao)
"""

ARGS = ("g0, gend, nranks, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sD, "
        "npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, "
        "pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, pp_cs, pp_ki, pp_kj, "
        "ncoltot, ncoef, ps_np, ps_ncol, ps_soff, ps_coff, col_ao, ps_coef, boys_d, "
        "ndens, dmat, jmat")


def _atomics_to_cuda(text):
    """`!$acc atomic update` + `a = a + x` (two-line form) -> `aold = atomicadd(a, x)`."""
    pat = re.compile(
        r"^(\s*)!\$acc atomic update\n"
        r"\s*jmat\((.*?)\) = jmat\((.*?)\) &\n"
        r"\s*\+ (.*?)\n", re.M)
    def rep(m):
        assert m.group(2) == m.group(3), m.group(0)
        return f"{m.group(1)}aold = atomicadd(jmat({m.group(2)}), {m.group(4)})\n"
    out, n = pat.subn(rep, text)
    assert "!$acc atomic" not in out, "an atomic update escaped the rewrite"
    return out


def emit_class(la, lb, lc, ld, vrr_body, hrr_body, block_txt, ncum, L):
    lab, lcd = la + lb, lc + ld
    lt = lab + lcd
    nca, ncc = ncum(lab), ncum(lcd)
    na, nb = (la + 1)*(la + 2)//2, (lb + 1)*(lb + 2)//2
    nc, nd = (lc + 1)*(lc + 2)//2, (ld + 1)*(ld + 2)//2
    nv = nca*ncc
    n4 = na*nb*nc*nd
    tag = f"{la}{lb}{lc}{ld}"
    gb = max(1, min(PG_T, SHM_G//(8*nv)))
    cur = lt % 2
    # A thread must be able to own one full ket column of the class even at
    # a single column per chunk: nv/PG_T entries, 20 at (dd|dd).
    ntc_max = max(PG_NTC, -(-nv // PG_T))
    ng_max = max(PG_NG, ntc_max)

    # --- verbatim slices of the do concurrent block kernel ---------------
    ihead = block_txt.index("      ! locate the segment")
    igc = block_txt.index("         ! GENERAL CONTRACTION.")
    head = block_txt[ihead:igc]
    head = head.replace("      ! locate the segment", "         ! locate the segment (every thread alike)")
    ihrr = block_txt.index("         ! --- HRR ---\n")
    iend = block_txt.index("         end do   ! qab\n")
    tail = _atomics_to_cuda(block_txt[ihrr:iend])
    assert "cycle" not in tail and "return" not in tail
    hrr_cuda = re.sub(r"\bg\(", "g1(", hrr_body)
    assert hrr_cuda in tail or hrr_body in tail

    boys = block_txt[block_txt.index("                  if (tval >= BOYS_TMAX) then"):
                     block_txt.index("                  oo2z = 0.5_dp/zeta")]
    boys = boys.replace("boys_table(", "boys_d(")

    mod = f"trc_pg_k{tag}"
    return f"""!
! Class ({la}{lb}|{lc}{ld}) cooperative kernel, CUDA Fortran.
!
! GENERATED by scripts/gen_perclass.py --lmax {L} --split --cuda -- do not edit.
! See scripts/gen_cuda.py for the shape and the reasons.
!
! Compiled to an empty module unless TRC_CUDAF is defined, so every other
! configuration -- host CI, fpm, the OpenMP port -- never sees CUDA Fortran.
!
module {mod}
#ifdef TRC_CUDAF
   use cudafor
   use trc_boys, only: dp, BOYS_MMAX, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV
   implicit none
   private
   public :: pg{tag}

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp
   integer, parameter :: PG_T = {PG_T}, PG_XC = {PG_XC}, PG_NTC = {ntc_max}, PG_NG = {ng_max}
   integer, parameter :: PG_NCC = {PG_NCC}, PG_NCAB = {PG_NCAB}, PG_GB = {gb}

contains

   attributes(global) subroutine pg{tag}({ARGS})
{CUDA_PROLOGUE}
      real(dp), shared :: Vs(PG_XC, PG_T), Ws(PG_T, PG_NCC), Gs({nv}, PG_GB)
      integer(kind=8) :: gt
      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl
      integer(kind=8) :: nsa, u, kx
      real(dp) :: qcut, pcut
      integer :: keyab, keycd, offab, offcd, nab, ncd
      integer :: kp, kq, kq0, d, x, x0, xc, cur, ia, ib, ic, id, idx, idens
      integer :: mu, nu, lam, sig, mui, nuj, lamk, sigl
      logical :: dij, dkl, dpq, have, ok
      real(dp) :: zeta, eta, zpe, rho, tval, pref, wc, wmax, acc, aold
      real(dp) :: pqx, pqy, pqz, pax, pay, paz, qcx, qcy, qcz
      real(dp) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re, sc, vv
      real(dp) :: abx, aby, abz, cdx, cdy, cdz
      real(dp) :: f(0:BOYS_MMAX)
      integer  :: bi, bj, bbase
      real(dp) :: bx, bx2, b0, b1, b2, btt, bet
      real(dp) :: v({nv}, 0:1), g1({nv}), vbuf({n4})
      real(dp) :: wq
      integer  :: nca, ncb, nccl, ncdl, npi, npj, npk, npl, ncab, nccd
      integer  :: ki, kj, kk, kl, ia2, ib2, ic2, id2, iabc, icdc
      logical  :: same_ab, same_cd, same_pair
      integer  :: tau, cdc, ncdc, ntc, ntcl, abc, nabc, ab0, cd0, a, c, j, e, l, ncomb, cb0, nbt, combo
      real(dp) :: cab(PG_NCAB), tcown(PG_NTC), gown(PG_NG)
      integer  :: offa(PG_NCAB), offb(PG_NCAB), offc(PG_NCC), offd(PG_NCC)
      real(dp) :: jab({na*nb}), jcd({nc*nd}), kac({na*nc})
      real(dp) :: kad({na*nd}), kbc({nb*nc}), kbd({nb*nd})
      real(dp) :: dab({na*nb}), dcd({nc*nd}), dac({na*nc})
      real(dp) :: dad({na*nd}), dbc({nb*nc}), dbd({nb*nd})

      tau = threadIdx%x
      gt = g0 + int(blockIdx%x - 1, 8)*int(nranks, 8)
      if (gt > gend) return
{head}
         ncab = nca*ncb; nccd = nccl*ncdl
         !
         ! COLUMN GRID CHUNKS, to the register budget: a thread owns at most
         ! PG_NTC ket-contracted entries and PG_NG bra-folded ones. Most
         ! classes fit in one pass; the rest repeat the primitive loops per
         ! chunk, which is still far fewer VRR evaluations than one per
         ! column combination.
         !
         cdc = max(1, min(nccd, (PG_NTC*PG_T)/{nv}))
         ntcl = ({nv}*cdc + PG_T - 1)/PG_T
         abc = max(1, min(ncab, PG_NG/ntcl))
         cur = {cur}
         do cd0 = 1, nccd, cdc
            ncdc = min(cdc, nccd - cd0 + 1)
            ntc = {nv}*ncdc
            ntcl = (ntc + PG_T - 1)/PG_T
            do c = 1, ncdc
               icdc = cd0 + c - 1
               ic2 = mod(icdc - 1, nccl) + 1; id2 = (icdc - 1)/nccl + 1
               offc(c) = ps_coff(sk) + (ic2 - 1)*npk
               offd(c) = ps_coff(sl) + (id2 - 1)*npl
            end do
            do ab0 = 1, ncab, abc
               nabc = min(abc, ncab - ab0 + 1)
               do a = 1, nabc
                  iabc = ab0 + a - 1
                  ia2 = mod(iabc - 1, nca) + 1; ib2 = (iabc - 1)/nca + 1
                  offa(a) = ps_coff(si) + (ia2 - 1)*npi
                  offb(a) = ps_coff(sj) + (ib2 - 1)*npj
               end do
               do j = 1, PG_NG
                  gown(j) = 0.0_dp
               end do

               do kp = offab + 1, offab + nab
                  zeta = pp_p(kp)
                  ki = pp_ki(kp); kj = pp_kj(kp)
                  wmax = 0.0_dp
                  do a = 1, nabc
                     cab(a) = ps_coef(offa(a) + ki)*ps_coef(offb(a) + kj)
                     wmax = max(wmax, abs(cab(a)))
                  end do
                  do j = 1, PG_NTC
                     tcown(j) = 0.0_dp
                  end do

                  do kq0 = offcd + 1, offcd + ncd, PG_T
                     ! The previous chunk's reduction reads Ws; nobody may
                     ! overwrite a row until every thread is past it.
                     call syncthreads()
                     kq = kq0 + tau - 1
                     have = kq <= offcd + ncd
                     if (have) then
                        eta = pp_p(kq)
                        kk = pp_ki(kq); kl = pp_kj(kq)
                        zpe = zeta + eta
                        pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)
                     else
                        eta = 1.0_dp; kk = 1; kl = 1; zpe = 1.0_dp; pref = 0.0_dp
                     end if
                     ! ket column weights of this lane's pair, and the
                     ! prescreen on the largest weight either side can give.
                     acc = 0.0_dp
                     do c = 1, ncdc
                        wc = ps_coef(offc(c) + kk)*ps_coef(offd(c) + kl)
                        if (.not. have) wc = 0.0_dp
                        Ws(tau, c) = wc
                        acc = max(acc, abs(wc))
                     end do
                     if (abs(pref)*wmax*acc <= pcut) have = .false.
                     if (have) then
                        rho = zeta*eta/zpe
                        pqx = pp_r(kp, 1) - pp_r(kq, 1)
                        pqy = pp_r(kp, 2) - pp_r(kq, 2)
                        pqz = pp_r(kp, 3) - pp_r(kq, 3)
                        pax = pp_r(kp, 1) - pp_ra(kp, 1)
                        pay = pp_r(kp, 2) - pp_ra(kp, 2)
                        paz = pp_r(kp, 3) - pp_ra(kp, 3)
                        qcx = pp_r(kq, 1) - pp_ra(kq, 1)
                        qcy = pp_r(kq, 2) - pp_ra(kq, 2)
                        qcz = pp_r(kq, 3) - pp_ra(kq, 3)
                        wc = (zeta*pp_r(kp, 1) + eta*pp_r(kq, 1))/zpe
                        wpx = wc - pp_r(kp, 1); wqx = wc - pp_r(kq, 1)
                        wc = (zeta*pp_r(kp, 2) + eta*pp_r(kq, 2))/zpe
                        wpy = wc - pp_r(kp, 2); wqy = wc - pp_r(kq, 2)
                        wc = (zeta*pp_r(kp, 3) + eta*pp_r(kq, 3))/zpe
                        wpz = wc - pp_r(kp, 3); wqz = wc - pp_r(kq, 3)
                        tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)
{boys}
                        oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = 0.5_dp/zpe
                        rz = rho/zeta; re = rho/eta
{vrr_body}
                     end if

                     ! KET CONTRACTION ACROSS THE LANES: tc(x, c) += sum over
                     ! the chunk's pairs of W(pair, c) v(x), staged through
                     ! shared memory PG_XC components at a time.
                     do x0 = 1, {nv}, PG_XC
                        xc = min(PG_XC, {nv} - x0 + 1)
                        call syncthreads()
                        do x = 1, xc
                           if (have) then
                              Vs(x, tau) = v(x0 + x - 1, cur)
                           else
                              Vs(x, tau) = 0.0_dp
                           end if
                        end do
                        call syncthreads()
                        do j = 1, ntcl
                           e = tau + (j - 1)*PG_T
                           if (e > ntc) exit
                           x = mod(e - 1, {nv}) + 1
                           c = (e - 1)/{nv} + 1
                           if (x >= x0 .and. x < x0 + xc) then
                              acc = 0.0_dp
                              do l = 1, PG_T
                                 acc = acc + Vs(x - x0 + 1, l)*Ws(l, c)
                              end do
                              tcown(j) = tcown(j) + acc
                           end if
                        end do
                     end do
                  end do   ! kq0

                  ! BRA FOLD, by the owners: g(x, a, c) += c_a c_b tc(x, c).
                  do j = 1, ntcl
                     e = tau + (j - 1)*PG_T
                     if (e > ntc) exit
                     do a = 1, nabc
                        gown(j + (a - 1)*ntcl) = gown(j + (a - 1)*ntcl) + cab(a)*tcown(j)
                     end do
                  end do
               end do   ! kp

               !
               ! DIGESTION, one thread per column combination of the chunk,
               ! in batches of PG_GB: the owners scatter into shared memory,
               ! the thread takes its column and runs the HRR and the folded
               ! digestion exactly as the do concurrent kernels do.
               !
               ncomb = nabc*ncdc
               do cb0 = 1, ncomb, PG_GB
                  nbt = min(PG_GB, ncomb - cb0 + 1)
                  call syncthreads()
                  do j = 1, ntcl
                     e = tau + (j - 1)*PG_T
                     if (e > ntc) exit
                     x = mod(e - 1, {nv}) + 1
                     c = (e - 1)/{nv} + 1
                     do a = 1, nabc
                        combo = a + (c - 1)*nabc
                        if (combo >= cb0 .and. combo < cb0 + nbt) then
                           Gs(x, combo - cb0 + 1) = gown(j + (a - 1)*ntcl)
                        end if
                     end do
                  end do
                  call syncthreads()
                  if (tau <= nbt) then
                     combo = cb0 + tau - 1
                     a = mod(combo - 1, nabc) + 1
                     c = (combo - 1)/nabc + 1
                     iabc = ab0 + a - 1
                     icdc = cd0 + c - 1
                     ia2 = mod(iabc - 1, nca) + 1; ib2 = (iabc - 1)/nca + 1
                     ic2 = mod(icdc - 1, nccl) + 1; id2 = (icdc - 1)/nccl + 1
                     ! Only the canonical contracted quartets, as the
                     ! contracted path enumerated them.
                     ok = .not. (same_ab .and. ia2 < ib2)
                     if (same_cd .and. ic2 < id2) ok = .false.
                     if (same_pair .and. iabc < icdc) ok = .false.
                     if (ok) then
                        dij = .not. (same_ab .and. ia2 == ib2)
                        dkl = .not. (same_cd .and. ic2 == id2)
                        dpq = .not. (same_pair .and. iabc == icdc)
                        mui = col_ao(ps_soff(si) + ia2); nuj = col_ao(ps_soff(sj) + ib2)
                        lamk = col_ao(ps_soff(sk) + ic2); sigl = col_ao(ps_soff(sl) + id2)
                        do x = 1, {nv}
                           g1(x) = Gs(x, tau)
                        end do
{tail}
                     end if
                  end if
               end do   ! cb0
            end do   ! ab0
         end do   ! cd0
   end subroutine pg{tag}
#endif
end module {mod}
"""


def emit_dispatch(names, L):
    """names: list of (key, tag). Returns the text of trc_pg_kernels.F90."""
    uses = "\n".join(f"   use trc_pg_k{tag}, only: pg{tag}" for _k, tag in names)
    cases = "\n".join(
        f"      case ({key}); call pg{tag}<<<grid, PG_T, 0, strm>>>(" + ARGS + ")"
        for key, tag in names)
    return f"""!
! Launcher for the CUDA Fortran cooperative kernels.
!
! GENERATED by scripts/gen_perclass.py --lmax {L} --split --cuda -- do not edit.
!
! The dummies are C pointers so that trc_binkernel, which is not compiled as
! CUDA Fortran, can call this from inside an `!$acc host_data use_device`
! region with `c_loc` of every OpenACC-resident array. Here they become
! device pointer arrays and the class kernel is launched on OpenACC's own
! stream, so it is ordered against the `do concurrent` kernels and the data
! transfers around it without a synchronisation of its own.
!
module trc_pg_kernels
#ifdef TRC_CUDAF
   use, intrinsic :: iso_c_binding, only: c_ptr
   use cudafor
   use openacc, only: acc_get_cuda_stream, acc_async_sync
   use trc_boys, only: dp, boys_table
{uses}
   implicit none
   private
   public :: pg_init, pg_dispatch

   integer, parameter :: PG_T = {PG_T}
   !> The Boys table, held by CUDA Fortran rather than reached through
   !> OpenACC: `declare create` module data has no `host_data` handle.
   real(dp), device, allocatable :: boys_d(:)

contains

   subroutine pg_init()
      if (.not. allocated(boys_d)) then
         allocate (boys_d(size(boys_table)))
         boys_d = boys_table
      end if
   end subroutine pg_init

   subroutine pg_dispatch(key, lo, hi, nseg, rank, nranks, npair, nbas, npp, nao, ncoltot, ncoef, ndens, &
                          thresh, jfac, kfac, p_sOff, p_sA, p_sNB, p_sOA, p_sOB, p_sD, p_sp_i, p_sp_j, p_sp_q, &
                          p_dsh, p_sh_l, p_ao_off, p_pp_off, p_pp_n, p_pp_p, p_pp_r, p_pp_ra, p_pp_rb, p_pp_c, &
                          p_pp_cs, p_pp_ki, p_pp_kj, p_ps_np, p_ps_ncol, p_ps_soff, p_ps_coff, p_col_ao, &
                          p_ps_coef, p_dmat, p_jmat, h_sOff_lo, h_sOff_hi1)
      integer, intent(in) :: key, lo, hi, nseg, rank, nranks, npair, nbas, npp, nao, ncoltot, ncoef, ndens
      real(dp), intent(in) :: thresh, jfac, kfac
      type(c_ptr), intent(in) :: p_sOff, p_sA, p_sNB, p_sOA, p_sOB, p_sD, p_sp_i, p_sp_j, p_sp_q
      type(c_ptr), intent(in) :: p_dsh, p_sh_l, p_ao_off, p_pp_off, p_pp_n, p_pp_p, p_pp_r, p_pp_ra, p_pp_rb
      type(c_ptr), intent(in) :: p_pp_c, p_pp_cs, p_pp_ki, p_pp_kj, p_ps_np, p_ps_ncol, p_ps_soff, p_ps_coff
      type(c_ptr), intent(in) :: p_col_ao, p_ps_coef, p_dmat, p_jmat
      !> sOff(lo) and sOff(hi+1) on the host: the launch geometry.
      integer(kind=8), intent(in) :: h_sOff_lo, h_sOff_hi1

      integer(kind=8), device, pointer :: sOff(:)
      integer,  device, pointer :: sA(:), sNB(:), sOA(:), sOB(:), sp_i(:), sp_j(:), sh_l(:), ao_off(:)
      logical,  device, pointer :: sD(:)
      real(dp), device, pointer :: sp_q(:), dsh(:, :), pp_p(:), pp_r(:, :), pp_ra(:, :), pp_rb(:, :)
      real(dp), device, pointer :: pp_c(:), pp_cs(:), ps_coef(:), dmat(:, :, :), jmat(:, :, :)
      integer,  device, pointer :: pp_off(:), pp_n(:), pp_ki(:), pp_kj(:)
      integer,  device, pointer :: ps_np(:), ps_ncol(:), ps_soff(:), ps_coff(:), col_ao(:)
      integer(kind=8) :: g0, gend, nr
      integer :: grid
      integer(kind=cuda_stream_kind) :: strm
      type(c_devptr) :: dp0

      g0 = h_sOff_lo + 1 + rank
      gend = h_sOff_hi1
      nr = 0
      if (gend >= g0) nr = (gend - g0)/nranks + 1
      if (nr <= 0) return
      grid = int(nr)

      call c_f_pointer(transfer(p_sOff, dp0), sOff, [nseg + 1])
      call c_f_pointer(transfer(p_sA, dp0), sA, [nseg])
      call c_f_pointer(transfer(p_sNB, dp0), sNB, [nseg])
      call c_f_pointer(transfer(p_sOA, dp0), sOA, [nseg])
      call c_f_pointer(transfer(p_sOB, dp0), sOB, [nseg])
      call c_f_pointer(transfer(p_sD, dp0), sD, [nseg])
      call c_f_pointer(transfer(p_sp_i, dp0), sp_i, [npair])
      call c_f_pointer(transfer(p_sp_j, dp0), sp_j, [npair])
      call c_f_pointer(transfer(p_sp_q, dp0), sp_q, [npair])
      call c_f_pointer(transfer(p_dsh, dp0), dsh, [nbas, nbas])
      call c_f_pointer(transfer(p_sh_l, dp0), sh_l, [nbas])
      call c_f_pointer(transfer(p_ao_off, dp0), ao_off, [nbas])
      call c_f_pointer(transfer(p_pp_off, dp0), pp_off, [nbas*nbas])
      call c_f_pointer(transfer(p_pp_n, dp0), pp_n, [nbas*nbas])
      call c_f_pointer(transfer(p_pp_p, dp0), pp_p, [npp])
      call c_f_pointer(transfer(p_pp_r, dp0), pp_r, [npp, 3])
      call c_f_pointer(transfer(p_pp_ra, dp0), pp_ra, [npp, 3])
      call c_f_pointer(transfer(p_pp_rb, dp0), pp_rb, [npp, 3])
      call c_f_pointer(transfer(p_pp_c, dp0), pp_c, [npp])
      call c_f_pointer(transfer(p_pp_cs, dp0), pp_cs, [npp])
      call c_f_pointer(transfer(p_pp_ki, dp0), pp_ki, [npp])
      call c_f_pointer(transfer(p_pp_kj, dp0), pp_kj, [npp])
      call c_f_pointer(transfer(p_ps_np, dp0), ps_np, [nbas])
      call c_f_pointer(transfer(p_ps_ncol, dp0), ps_ncol, [nbas])
      call c_f_pointer(transfer(p_ps_soff, dp0), ps_soff, [nbas])
      call c_f_pointer(transfer(p_ps_coff, dp0), ps_coff, [nbas])
      call c_f_pointer(transfer(p_col_ao, dp0), col_ao, [ncoltot])
      call c_f_pointer(transfer(p_ps_coef, dp0), ps_coef, [ncoef])
      call c_f_pointer(transfer(p_dmat, dp0), dmat, [ndens, nao, nao])
      call c_f_pointer(transfer(p_jmat, dp0), jmat, [ndens, nao, nao])

      strm = acc_get_cuda_stream(acc_async_sync)
      select case (key)
{cases}
      case default
         error stop "trc_pg_kernels: no cooperative kernel for this class"
      end select
   end subroutine pg_dispatch
#else
   implicit none
   private
#endif
end module trc_pg_kernels
"""
