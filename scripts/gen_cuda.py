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
One thread block of one warp per primitive-shell quartet item of a general
segment run. The warp walks the bra primitive pairs; for each, the ket pairs
go in chunks of 32 with one lane per ket pair computing the primitive
prefactor, the prescreen, the Boys function and the sieved VRR block v(nv)
in registers -- the scalar kernel's code -- and accumulating its own
quartets' ket-column contraction tcl(x, c) += c_c c_d v(x) in registers,
all of it literal-indexed so nothing leaves the register file. Once per
bra pair the nv*CDC entries are summed across the lanes by shuffle
butterflies and the owner lane folds the bra columns in, gown(x, a, c) +=
c_a c_b tc(x, c). So the factorised cost of the blocked kernel is kept --
nv*CDC FMAs per primitive quartet, the reduction and nv*ABC*CDC per bra
pair -- with the VRR computed once. The column grid is chunked at
generation time to a register budget (PG_ACC accumulators beside the VRR
block); a class whose grid does not fit repeats its primitive loops per
chunk, never more often than the blocked kernel did.

Digestion is one lane per column combination, in batches of PG_GB: the
owners scatter their g into shared memory, the lane copies its column to a
local g1 and runs the generated HRR and six-block folded digestion
unchanged, with `!$acc atomic update` become `atomicadd`. Those two
barriers are the only ones in the kernel.
"""
import re

PG_T = 32      # one warp per block: a ket primitive pair per lane
PG_ACC = 64    # registers of per-lane ket-contracted accumulator, at most
PG_CDC = 4     # ket column pairs per chunk, at most: shell pairs have 1-9
PG_NCAB = 16   # bra column pairs per chunk (PS_NCOL_MAX**2), one per lane
SHM_G = 24576  # bytes of shared memory for the digestion gather

CUDA_PROLOGUE = """      integer, value :: nranks, lo, hi, nseg, npair, nbas, npp, nao, ncoltot, ncoef, ndens
      integer(int64), value :: g0, gend
      real(dp), value :: thresh, jfac, kfac
      integer(int64), device :: sOff(nseg + 1)
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
    # Ket columns accumulated per lane at once: nv*cdc registers of
    # accumulator beside the nv of the VRR block, literal-indexed so they
    # stay registers. Bra columns folded per chunk to a fixed local budget.
    cdc = max(1, min(PG_CDC, PG_ACC//nv))
    ntc = nv*cdc                       # entries reduced per bra pair
    abc = PG_NCAB                      # every bra column at once: lane a owns column a

    # --- verbatim slices of the do concurrent block kernel ---------------
    ihead = block_txt.index("      ! locate the segment")
    igc = block_txt.index("         ! GENERAL CONTRACTION.")
    head = block_txt[ihead:igc]
    head = head.replace("      ! locate the segment", "         ! locate the segment (every lane alike)")
    ihrr = block_txt.index("         ! --- HRR ---\n")
    iend = block_txt.index("         end do   ! qab\n")
    tail = _atomics_to_cuda(block_txt[ihrr:iend])
    assert "cycle" not in tail and "return" not in tail
    boys = block_txt[block_txt.index("                  if (tval >= BOYS_TMAX) then"):
                     block_txt.index("                  oo2z = 0.5_dp/zeta")]
    boys = boys.replace("boys_table(", "boys_d(")

    # Per-lane accumulation, fully unrolled with literal indices.
    accum = "".join(
        f"                        tcl({x + (c - 1)*nv}) = tcl({x + (c - 1)*nv}) + wcol({c})*v({x}, {cur})\n"
        for c in range(1, cdc + 1) for x in range(1, nv + 1))
    tcl_zero = "".join(f"                     tcl({e}) = 0.0_dp\n" for e in range(1, ntc + 1))
    # Warp reduction of every entry, once per bra pair; the owner folds
    # the bra columns in. Butterfly, so every lane ends with the total.
    reduce = "".join(
        (f"                  if ({(e - 1)//nv + 1} <= ncdc) then\n" if (e - 1) % nv == 0 else "") +
        f"""                  val = tcl({e})
                  val = val + __shfl_xor(val, 16)
                  val = val + __shfl_xor(val, 8)
                  val = val + __shfl_xor(val, 4)
                  val = val + __shfl_xor(val, 2)
                  val = val + __shfl_xor(val, 1)
                  gown({e}) = gown({e}) + cabl*val
""" + ("                  end if\n" if e % nv == 0 else "")
        for e in range(1, ntc + 1))

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
   use, intrinsic :: iso_fortran_env, only: int64
   use cudafor
   use trc_boys, only: dp, BOYS_MMAX, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV
   implicit none
   private
   public :: pg{tag}

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp
   integer, parameter :: PG_T = {PG_T}, PG_GB = {gb}
   !> ket columns per chunk, bra columns per chunk, owned entries per lane
   integer, parameter :: CDC = {cdc}, ABC = {abc}, NTC = {ntc}

contains

   attributes(global) subroutine pg{tag}({ARGS})
{CUDA_PROLOGUE}
      real(dp), shared :: Gs({nv}, PG_GB)
      integer(int64) :: gt
      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl
      integer(int64) :: nsa, u, kx
      real(dp) :: qcut, pcut
      integer :: keyab, keycd, offab, offcd, nab, ncd
      integer :: kp, kq, kq0, d, x, cur, ia, ib, ic, id, idx, idens
      integer :: mu, nu, lam, sig, mui, nuj, lamk, sigl
      logical :: dij, dkl, dpq, have, ok
      real(dp) :: zeta, eta, zpe, rho, tval, pref, wc, wmax, acc, aold, val
      real(dp) :: pqx, pqy, pqz, pax, pay, paz, qcx, qcy, qcz
      real(dp) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re, sc, vv
      real(dp) :: abx, aby, abz, cdx, cdy, cdz
      real(dp) :: f(0:BOYS_MMAX)
      integer  :: bi, bj, bbase
      real(dp) :: bx, bx2, b0, b1, b2, btt, bet
      real(dp) :: v({nv}, 0:1), g1({nv}), vbuf({n4})
      real(dp) :: tcl(NTC), wcol(CDC), gown(NTC), cabl
      real(dp) :: wq
      integer  :: nca, ncb, nccl, ncdl, npi, npj, npk, npl, ncab, nccd
      integer  :: ki, kj, kk, kl, ia2, ib2, ic2, id2, iabc, icdc
      logical  :: same_ab, same_cd, same_pair
      integer  :: tau, ncdc, nabc, ab0, cd0, a, c, j, ncomb, cb0, nbt, combo
      integer  :: offal, offbl, offc(CDC), offd(CDC)
      real(dp) :: jab({na*nb}), jcd({nc*nd}), kac({na*nc})
      real(dp) :: kad({na*nd}), kbc({nb*nc}), kbd({nb*nd})
      real(dp) :: dab({na*nb}), dcd({nc*nd}), dac({na*nc})
      real(dp) :: dad({na*nd}), dbc({nb*nc}), dbd({nb*nd})

      tau = threadIdx%x
      gt = g0 + int(blockIdx%x - 1, int64)*int(nranks, int64)
      if (gt > gend) return
      pcut = thresh*1.0e-3_dp
{head}
         ncab = nca*ncb; nccd = nccl*ncdl
         cur = {cur}
         !
         ! COLUMN GRID in chunks of ABC x CDC. A lane accumulates its own
         ! quartets' ket-contracted block in registers over the ket loop;
         ! the sum across the warp's lanes is taken once per bra pair, by
         ! shuffles, and the owner folds the bra columns in. Neither needs
         ! shared memory or a barrier. A class whose grid exceeds a chunk
         ! repeats its primitive loops per chunk, still far fewer VRR
         ! evaluations than one per column combination.
         !
         do cd0 = 1, nccd, CDC
            ncdc = min(CDC, nccd - cd0 + 1)
            do c = 1, CDC
               icdc = min(cd0 + c - 1, nccd)
               ic2 = mod(icdc - 1, nccl) + 1; id2 = (icdc - 1)/nccl + 1
               offc(c) = ps_coff(sk) + (ic2 - 1)*npk
               offd(c) = ps_coff(sl) + (id2 - 1)*npl
            end do
            do ab0 = 1, ncab, ABC
               nabc = min(ABC, ncab - ab0 + 1)
               ! Lane a owns bra column a of the chunk: its coefficient
               ! offsets, its weight per bra pair, and its slice of g.
               iabc = min(ab0 + tau - 1, ncab)
               ia2 = mod(iabc - 1, nca) + 1; ib2 = (iabc - 1)/nca + 1
               offal = ps_coff(si) + (ia2 - 1)*npi
               offbl = ps_coff(sj) + (ib2 - 1)*npj
               do j = 1, NTC
                  gown(j) = 0.0_dp
               end do

               do kp = offab + 1, offab + nab
                  zeta = pp_p(kp)
                  ki = pp_ki(kp); kj = pp_kj(kp)
                  cabl = 0.0_dp
                  if (tau <= nabc) cabl = ps_coef(offal + ki)*ps_coef(offbl + kj)
                  wmax = abs(cabl)
                  wmax = max(wmax, __shfl_xor(wmax, 16))
                  wmax = max(wmax, __shfl_xor(wmax, 8))
                  wmax = max(wmax, __shfl_xor(wmax, 4))
                  wmax = max(wmax, __shfl_xor(wmax, 2))
                  wmax = max(wmax, __shfl_xor(wmax, 1))
{tcl_zero}
                  do kq0 = offcd + 1, offcd + ncd, PG_T
                     kq = kq0 + tau - 1
                     have = kq <= offcd + ncd
                     if (have) then
                        eta = pp_p(kq)
                        kk = pp_ki(kq); kl = pp_kj(kq)
                        zpe = zeta + eta
                        pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)
                        acc = 0.0_dp
                        do c = 1, CDC
                           wc = ps_coef(offc(c) + kk)*ps_coef(offd(c) + kl)
                           if (c > ncdc) wc = 0.0_dp
                           wcol(c) = wc
                           acc = max(acc, abs(wc))
                        end do
                        if (abs(pref)*wmax*acc <= pcut) have = .false.
                     end if
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
{accum}                     end if
                  end do   ! kq0

                  ! WARP REDUCTION, once per bra pair, and the bra fold.
{reduce}               end do   ! kp

               !
               ! DIGESTION, one lane per column combination of the chunk,
               ! in batches of PG_GB: the owners scatter into shared memory,
               ! the lane takes its column and runs the HRR and the folded
               ! digestion exactly as the do concurrent kernels do.
               !
               ncomb = nabc*ncdc
               do cb0 = 1, ncomb, PG_GB
                  nbt = min(PG_GB, ncomb - cb0 + 1)
                  call syncthreads()
                  if (tau <= nabc) then
                     do c = 1, ncdc
                        combo = tau + (c - 1)*nabc
                        if (combo >= cb0 .and. combo < cb0 + nbt) then
                           do x = 1, {nv}
                              Gs(x, combo - cb0 + 1) = gown(x + (c - 1)*{nv})
                           end do
                        end if
                     end do
                  end if
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
      integer(int64), intent(in) :: h_sOff_lo, h_sOff_hi1

      integer(int64), device, pointer :: sOff(:)
      integer,  device, pointer :: sA(:), sNB(:), sOA(:), sOB(:), sp_i(:), sp_j(:), sh_l(:), ao_off(:)
      logical,  device, pointer :: sD(:)
      real(dp), device, pointer :: sp_q(:), dsh(:, :), pp_p(:), pp_r(:, :), pp_ra(:, :), pp_rb(:, :)
      real(dp), device, pointer :: pp_c(:), pp_cs(:), ps_coef(:), dmat(:, :, :), jmat(:, :, :)
      integer,  device, pointer :: pp_off(:), pp_n(:), pp_ki(:), pp_kj(:)
      integer,  device, pointer :: ps_np(:), ps_ncol(:), ps_soff(:), ps_coff(:), col_ao(:)
      integer(int64) :: g0, gend, nr
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
