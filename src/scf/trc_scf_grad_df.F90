!
! The two-electron part of the density-fitted (`--df`) RHF/RKS nuclear
! gradient: fitted J and fitted K, both through the same auxiliary basis and
! the same metric Cholesky factor the DF-SCF (trc_scf_df.F90) converged with.
!
! THE ENERGY. With M = (P|Q) = L L^T, t_P = sum_mn D_mn (mn|P) and, for the
! occupied orbitals C (D = occ_s C C^T, occ_s = 2), h^{ij}_P = sum_mn C_mi
! (mn|P) C_nj, the fitted two-electron energy terco converges is
!
!     E_2 = (j/2) t^T M^-1 t  -  (k/4) occ_s^2 sum_ij h^{ij} M^-1 h^{ij}
!
! (j = 1, k the exact-exchange fraction: 1 for HF, func%exx for a hybrid,
! the same G = j J - (k/2) K split trc_scf_df's build uses).
!
! THE GRADIENT. Stationary in C (the SCF converged), so only the integrals
! move. With c = M^-1 t and X^{ij} = M^-1 h^{ij} -- one triangular solve
! each against the RETAINED L (df_jk_t%cholesky), since df_jk_t already
! holds B = T L^-T, whose contractions are L^-1 t and L^-1 h:
!
!     dE_2/dR = sum_(mn P) Y^P_mn d(mn|P)/dR + sum_(P Q) Z_PQ d(P|Q)/dR
!
!     Y^P_mn = j c_P D_mn - (k/2) occ_s^2 [C X^P C^T]_mn
!     Z_PQ   = -(j/2) c_P c_Q + (k/4) occ_s^2 sum_ij X^{ij}_P X^{ij}_Q
!
! Both coefficient tensors multiply derivatives taken with respect to the
! nuclear centres (scripts/gen_df_grad.py computes d/dA and d/dB and sets dC
! from translational invariance), so no sign convention of any integral
! library enters: Y and Z go in exactly as written. Verified numerically
! against finite differences before this was written (the K 2c term's sign
! in particular, which a library's nabla-on-the-electron convention flips).
!
! THE ORDER OF WORK
!   1. c: one GEMV on packed B (t' = B^T w D, the energy's own J step) and
!      one triangular solve, L^T c = t'.
!   2. K only (k /= 0; a pure functional skips the K work outright, not
!      by scaling by zero): H(ij, P) = C^T B_P C one auxiliary block at a
!      time (unpack, then one GEMM and one strided batch), kept only for
!      j <= i (X^{ij} = X^{ji}; the off-diagonal rows scaled by sqrt 2),
!      then X = H L^-1, one triangular solve over every kept row at once.
!   3. Z: an outer product plus, for K, X^T X over the kept rows, its lower
!      triangle only, one column panel per GEMM. Then the two-centre
!      derivative kernels (trc_df_grad_kernels' df_grad_2c_driver).
!   4. The three-centre derivative kernels, on the pair plan built once
!      (df_grad_3c_plan) and Y = j c_P D_mn + Y_K, whose J part the kernels
!      form themselves from c and D (df_grad_3c_driver). J only: one pass
!      over the whole auxiliary basis, no Y tensor. With K: Y_K one
!      auxiliary block at a time from the kept rows of X (unpacked to the
!      full X^P with one chunk; as M with Y_K = M + M^T with several), one
!      strided batch and one GEMM, and the kernels on that block's own
!      auxiliary subset.
!
! MEMORY. X is nocc^2 naux doubles. When the device cannot hold it next to
! everything else (choose_grad_blocks), steps 2-4 run once per CHUNK of
! occupied orbitals i -- K is a sum over i, so each chunk's rows of X give
! its own share of Z and of Y exactly -- and the three-centre kernels run
! once per chunk. A system too large for even one orbital per chunk stops
! with the sizes, not at a failed allocation.
!
! DEVICE RESIDENCY. Every argument array (cmo, dens, grad) and both bases
! and the pair list must already be device-resident, as must `dfjk`'s own
! B and L (they are, from its setup). Everything else is this routine's own.
!
module trc_scf_grad_df
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_scf_df, only: df_jk_t
   use trc_linalg, only: trc_linalg_t
   use trc_step_table, only: step_table_t
   use pic_timer, only: timer_type
   use trc_df_grad_kernels, only: df_grad_2c_driver, df_grad_3c_driver, df_grad_3c_plan, DFG_LMAX, &
                                  DFG_LMAX_AUX, DFG_NPCLS
#ifdef _OPENACC
   use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, &
                      acc_property_free_memory
#endif
   implicit none
   private
   public :: df_two_electron_gradient

   ! The timed steps, reported (global_logger, PERFORMANCE) after every
   ! gradient in the same table shape DF-JK and RI-MP2 use: device-synced
   ! seconds, and GFLOP/s for the dense algebra or GB/s for the steps that
   ! only move memory. The derivative kernels report time only -- their
   ! operation count is a property of each generated class, measured with
   ! a profiler, not modelled here.
   integer, parameter :: ST_JGAM = 1, ST_CSOL = 2, ST_UNPK = 3, ST_HGEMM = 4, ST_XPACK = 5, ST_XTRSM = 6, &
                         ST_ZGEMM = 7, ST_PLAN = 8, ST_YK = 9, ST_D3C = 10, ST_D2C = 11, NSTEP = 11
   real(dp), parameter :: SQRT2 = 1.4142135623730951_dp
   integer, parameter :: ZPANEL = 1024
   !! Columns per panel of Z's lower-triangle GEMMs.
   character(len=22), parameter :: NAMES(NSTEP) = [character(len=22) :: &
                                   "J gamma (gemv)", "c = L^-T gamma (trsm)", "unpack B", &
                                   "H = C^T B C (gemm)", "pack/unpack X (i>=j)", "X = H L^-1 (trsm)", &
                                   "Z += X^T X (gemm)", &
                                   "3c pair plan (host)", "Y = K part (gemm)", "3c derivative kernels", &
                                   "2c derivative kernels"]

contains

   subroutine df_two_electron_gradient(b, aux, pl, dfjk, lp, nocc, cmo, dens, j_scale, k_scale, grad)
      !! grad += dE_2/dR for the density-fitted J (scaled by j_scale) and K
      !! (scaled by k_scale), restricted, Cartesian -- see this module's
      !! header. `cmo` the converged orbitals (their first `nocc` columns are
      !! the occupied ones), `dens` the total density 2 C_occ C_occ^T,
      !! `pl` the pair list the fit's (mu nu|P) was built over.
      type(trc_basis_t),    intent(in)    :: b
      type(trc_basis_t),    intent(in)    :: aux
      type(trc_pairlist_t), intent(in)    :: pl
      type(df_jk_t),        intent(in), target :: dfjk
      type(trc_linalg_t),   intent(inout) :: lp
      integer,              intent(in)    :: nocc
      real(dp),             intent(in)    :: cmo(b%nao, b%nao)
      real(dp),             intent(in)    :: dens(b%nao, b%nao)
      real(dp),             intent(in)    :: j_scale, k_scale
      real(dp),             intent(inout) :: grad(3, b%natm)

      real(dp), pointer :: bpk(:, :), chol(:, :)
      real(dp), allocatable :: dpk(:), cq(:), zmat(:, :), cocc(:, :), xp(:), hb(:, :, :)
      real(dp), allocatable :: tens(:, :, :), half(:, :, :), ydum(:, :, :)
      type(trc_basis_t) :: blk
      integer :: n, naux, npack, ablk, nci, nchunk, ic, i0, i1, nc, s0, s1, p0, np, mu, nu, q, i, p
      integer(int64) :: n8, o8, c8, t
      integer :: q0, w, nrp, rbase, il, j, kmode
      integer :: coff(0:DFG_NPCLS)
      integer, allocatable :: q_a(:), q_b(:), q_o(:), q_n(:)
      real(dp), allocatable :: pd(:, :)
      logical :: do_k
      real(dp) :: occ_s, jc
      type(timer_type) :: tc
      type(step_table_t) :: tbl
      type(timer_type) :: tall
      logical :: dense(NSTEP)

      if (b%spherical) error stop "trc_scf_grad_df: the orbital basis is spherical; the DF gradient is Cartesian only"
      if (.not. dfjk%is_ready()) error stop "trc_scf_grad_df: the density fit has not been set up"
      if (maxval(b%sh_l) > DFG_LMAX .or. maxval(aux%sh_l) > DFG_LMAX_AUX) then
         error stop "trc_scf_grad_df: a basis is above the DF derivative kernels' ceiling (gen_df_grad.py)"
      end if
      call dfjk%dims(n, naux, npack)
      if (n /= b%nao .or. naux /= aux%nao) error stop "trc_scf_grad_df: the fit does not match these bases"
      if (aux%natm /= b%natm) error stop "trc_scf_grad_df: the auxiliary basis is on a different molecule"
      call tbl%init("dfgrad", NAMES, bytes_steps=[ST_JGAM, ST_UNPK, ST_XPACK])
      call lp%sync()
      call tall%start()
      bpk => dfjk%packed_b()
      chol => dfjk%cholesky()
      occ_s = 2.0_dp
      do_k = k_scale /= 0.0_dp .and. nocc > 0
      n8 = int(n, int64)
      o8 = int(max(nocc, 1), int64)
      ! --- 1. c = M^-1 t: the energy's own J gamma, then L^T c = gamma -----
      allocate (dpk(npack), cq(naux), zmat(naux, naux))
      !$acc enter data create(dpk, cq, zmat)
      ! One flat index (a masked two-index form left one loop sequential
      ! on the device): the lower triangle of D, off-diagonal doubled.
      do concurrent(t=0:n8*n8 - 1) local(mu, nu)
         mu = int(mod(t, n8)) + 1
         nu = int(t/n8) + 1
         if (nu <= mu) dpk(pack_idx(mu, nu)) = merge(1.0_dp, 2.0_dp, mu == nu)*dens(mu, nu)
      end do
      call tbl%step_on(lp, ST_JGAM, tc)
      call lp%gemv("T", npack, naux, 1.0_dp, bpk, npack, dpk, 0.0_dp, cq)
      call tbl%step_off(lp, ST_JGAM, tc)
      call tbl%add_bytes(ST_JGAM, 8_int64*int(npack, int64)*int(naux, int64))
      call tbl%step_on(lp, ST_CSOL, tc)
      call lp%trsm("L", "L", "T", naux, 1, 1.0_dp, chol, naux, cq, naux)
      call tbl%step_off(lp, ST_CSOL, tc)
      call tbl%add_flop(ST_CSOL, int(naux, int64)**2)

      do concurrent(q=1:naux, p=1:naux)
         zmat(p, q) = -0.5_dp*j_scale*cq(p)*cq(q)
      end do
      if (do_k) then
         allocate (cocc(n, nocc))
         !$acc enter data create(cocc)
         do concurrent(mu=1:n, i=1:nocc)
            cocc(mu, i) = cmo(mu, i)
         end do
      end if

      ! The three-centre kernels' pair plan: every screened pair turned to
      ! its canonical class, sorted by primitive count, its primitive pairs
      ! precomputed (trc_df_grad_kernels' df_grad_3c_plan), once.
      call tbl%step_on(lp, ST_PLAN, tc)
      call df_grad_3c_plan(b%nshell, b%maxnp, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, pl%npair, pl%p_i, &
                           pl%p_j, coff, q_a, q_b, q_o, q_n, pd)
      !$acc enter data copyin(q_a, q_b, q_o, q_n, pd)
      call tbl%step_off(lp, ST_PLAN, tc)

      ! The blocking, chosen only now that everything resident above is
      ! already allocated (the plan included): auxiliary functions per
      ! block, and occupied orbitals i per chunk of X (packed, X is
      ! nocc(nocc+1)/2 naux doubles, GBs at a few hundred occupied and
      ! several thousand auxiliary functions). One chunk is the usual case.
      call choose_grad_blocks(n, naux, merge(nocc, 0, do_k), ablk, nci)
      nchunk = 1
      if (do_k) nchunk = (nocc + nci - 1)/nci

      ! --- J only: Y = j c_P D_mn is formed inside the three-centre kernels
      ! from c and D -- no (n, n, naux) tensor, no auxiliary blocks, one pass.
      if (.not. do_k) then
         allocate (ydum(1, 1, 1))
         ydum = 0.0_dp
         !$acc enter data copyin(ydum)
         call tbl%step_on(lp, ST_D3C, tc)
         call df_grad_3c_driver(b%nshell, n, b%sh_ao, b%sh_at, b%sh_r, coff, pl%npair, q_a, q_b, q_o, q_n, &
                                size(pd, 2), pd, aux%nshell, aux%maxnp, aux%sh_l, aux%sh_np, aux%sh_ao, &
                                aux%sh_at, aux%sh_e, aux%sh_c, aux%sh_r, b%natm, 0, 1, 1, ydum, j_scale, &
                                naux, 0, cq, dens, grad)
         call tbl%step_off(lp, ST_D3C, tc)
         !$acc exit data delete(ydum)
         deallocate (ydum)
         nchunk = 0
      end if

      ! --- 2-4, with K, per occupied chunk [i0, i1] ------------------------
      ! K is a sum over occupied i, so X's rows i in [i0, i1] give that
      ! chunk's share of Z (X_c^T X_c) and of Y (C_c X_c^P C^T) on their
      ! own; the three-centre derivative is linear in Y, so it runs once
      ! per chunk on that chunk's Y. The J part (formed in the kernels)
      ! rides with the first chunk.
      do ic = 1, nchunk
         i0 = 1 + (ic - 1)*nci
         i1 = min(nocc, ic*nci)
         nc = max(i1 - i0 + 1, 0)
         c8 = int(max(nc, 1), int64)
         jc = merge(j_scale, 0.0_dp, ic == 1)

         ! --- 2. K: H(ij, P) = C_c^T B_P C for i in the chunk, X = H L^-1 --
         ! X^{ij} = X^{ji}, so only the rows j <= i are kept: xp(r, P), r the
         ! chunk-local packed index of (i, j), scaled by w_ij = sqrt(2) off
         ! the diagonal so that xp^T xp is Z's full sum over (i, j). Half the
         ! rows for the triangular solve and for Z (w60/6-31G: 5.7 and 11.5
         ! TFLOP before). Each auxiliary block's H lands in a block-sized
         ! buffer and is packed straight into xp: no full nc x nocc x naux H.
         nrp = (i1*(i1 + 1) - (i0 - 1)*i0)/2
         rbase = (i0 - 1)*i0/2
         allocate (xp(int(nrp, int64)*int(naux, int64)))
         !$acc enter data create(xp)
         s0 = 1
         do while (s0 <= aux%nshell)
            call next_block(aux, s0, ablk, s1, p0, np)
            allocate (tens(n, n, np), half(nc, n, np), hb(nc, nocc, np))
            !$acc enter data create(tens, half, hb)
            call tbl%step_on(lp, ST_UNPK, tc)
            ! One flat index over the whole (n, n, np) block: the three-
            ! index form ran its q loop sequentially inside each thread
            ! (nvfortran -Minfo), 62 GB/s on a V100.
            do concurrent(t=0:n8*n8*int(np, int64) - 1) local(mu, nu, q)
               mu = int(mod(t, n8)) + 1
               nu = int(mod(t/n8, n8)) + 1
               q = int(t/(n8*n8)) + 1
               tens(mu, nu, q) = bpk(pack_idx(max(mu, nu), min(mu, nu)), p0 - 1 + q)
            end do
            call tbl%step_off(lp, ST_UNPK, tc)
            call tbl%add_bytes(ST_UNPK, 8_int64*n8*n8*int(np, int64))
            call tbl%step_on(lp, ST_HGEMM, tc)
            call lp%gemm("T", "N", nc, n*np, n, 1.0_dp, cocc(1, i0), n, tens, n, 0.0_dp, half, nc)
            call lp%gemm_strided("N", "N", nc, nocc, n, 1.0_dp, half, nc, c8*n8, cocc, n, 0_int64, &
                                 0.0_dp, hb, nc, c8*o8, np)
            call tbl%step_off(lp, ST_HGEMM, tc)
            call tbl%add_flop(ST_HGEMM, 2_int64*c8*n8*int(np, int64)*(n8 + o8))
            call tbl%step_on(lp, ST_XPACK, tc)
            do concurrent(t=0:c8*o8*int(np, int64) - 1) local(il, j, q, i)
               il = int(mod(t, c8)) + 1
               j = int(mod(t/c8, o8)) + 1
               q = int(t/(c8*o8)) + 1
               i = i0 - 1 + il
               if (j <= i) xp(int(i*(i - 1)/2 + j - rbase, int64) + int(nrp, int64)*int(p0 - 2 + q, int64)) = &
                  merge(1.0_dp, SQRT2, j == i)*hb(il, j, q)
            end do
            call tbl%step_off(lp, ST_XPACK, tc)
            call tbl%add_bytes(ST_XPACK, 8_int64*c8*o8*int(np, int64))
            !$acc exit data delete(tens, half, hb)
            deallocate (tens, half, hb)
            s0 = s1 + 1
         end do
         call tbl%step_on(lp, ST_XTRSM, tc)
         call lp%trsm("R", "L", "N", nrp, naux, 1.0_dp, chol, naux, xp, nrp)
         call tbl%step_off(lp, ST_XTRSM, tc)
         call tbl%add_flop(ST_XTRSM, int(nrp, int64)*int(naux, int64)**2)
         ! --- 3a. Z's K share, this chunk: Z is symmetric, so only its lower
         ! triangle, one column panel at a time (each panel's own diagonal
         ! block in full), mirrored once after the last chunk -- half the
         ! flops of the one square GEMM, which also ran at only ~51% of peak
         ! (cuBLAS's 64x64 TN tile at this shape, gly10/cc-pVDZ).
         call tbl%step_on(lp, ST_ZGEMM, tc)
         do q0 = 1, naux, ZPANEL
            w = min(ZPANEL, naux - q0 + 1)
            call lp%gemm("T", "N", naux - q0 + 1, w, nrp, 0.25_dp*k_scale*occ_s*occ_s, &
                         xp(1 + int(q0 - 1, int64)*int(nrp, int64)), nrp, &
                         xp(1 + int(q0 - 1, int64)*int(nrp, int64)), nrp, 1.0_dp, zmat(q0, q0), naux)
            call tbl%add_flop(ST_ZGEMM, 2_int64*int(nrp, int64)*int(naux - q0 + 1, int64)*int(w, int64))
         end do
         call tbl%step_off(lp, ST_ZGEMM, tc)

         ! --- 4. Y, one auxiliary block at a time, and the three-centre one -
         ! One chunk (every i): V^P = X^P in full, rows j > i read from their
         ! (j, i) partner, and Y_K = C V C^T is symmetric (kmode 1). More
         ! than one chunk: the partner of a row j > i may sit in another
         ! chunk, so V^P_ij = X^P_ij below the diagonal, X^P_ii/2 on it and 0
         ! above, M = C_c V C^T and Y_K = M + M^T (X symmetric) -- the
         ! kernels read both M(mu,nu) and M(nu,mu) (kmode 2), a second,
         ! uncoalesced read (3% of the kernels' time at gly10/cc-pVDZ).
         kmode = merge(1, 2, nchunk == 1)
         s0 = 1
         do while (s0 <= aux%nshell)
            call next_block(aux, s0, ablk, s1, p0, np)
            call aux%subset(s0, s1, blk)
            call blk%to_device()
            allocate (tens(n, n, np), half(nc, n, np), hb(nc, nocc, np))
            !$acc enter data create(tens, half, hb)
            call tbl%step_on(lp, ST_XPACK, tc)
            do concurrent(t=0:c8*o8*int(np, int64) - 1) local(il, j, q, i)
               il = int(mod(t, c8)) + 1
               j = int(mod(t/c8, o8)) + 1
               q = int(t/(c8*o8)) + 1
               i = i0 - 1 + il
               if (j < i) then
                  hb(il, j, q) = xp(int(i*(i - 1)/2 + j - rbase, int64) + int(nrp, int64)*int(p0 - 2 + q, int64)) &
                                 /SQRT2
               else if (j == i) then
                  hb(il, j, q) = merge(1.0_dp, 0.5_dp, kmode == 1)* &
                                 xp(int(i*(i - 1)/2 + j - rbase, int64) + int(nrp, int64)*int(p0 - 2 + q, int64))
               else if (kmode == 1) then
                  hb(il, j, q) = xp(int(j*(j - 1)/2 + i - rbase, int64) + int(nrp, int64)*int(p0 - 2 + q, int64)) &
                                 /SQRT2
               else
                  hb(il, j, q) = 0.0_dp
               end if
            end do
            call tbl%step_off(lp, ST_XPACK, tc)
            call tbl%add_bytes(ST_XPACK, 8_int64*c8*o8*int(np, int64))
            ! half_q = V_q C^T (nc x n), then M = -(k/2) occ_s^2 C_c half,
            ! every q of the block in one GEMM; the J part is the kernels'.
            call tbl%step_on(lp, ST_YK, tc)
            call lp%gemm_strided("N", "T", nc, n, nocc, 1.0_dp, hb, nc, c8*o8, cocc, n, 0_int64, 0.0_dp, &
                                 half, nc, c8*n8, np)
            call lp%gemm("N", "N", n, n*np, nc, -0.5_dp*k_scale*occ_s*occ_s, cocc(1, i0), n, half, nc, &
                         0.0_dp, tens, n)
            call tbl%step_off(lp, ST_YK, tc)
            call tbl%add_flop(ST_YK, 2_int64*c8*n8*int(np, int64)*(n8 + o8))
            !$acc exit data delete(half, hb)
            deallocate (half, hb)
            call tbl%step_on(lp, ST_D3C, tc)
            call df_grad_3c_driver(b%nshell, n, b%sh_ao, b%sh_at, b%sh_r, coff, pl%npair, q_a, q_b, q_o, q_n, &
                                   size(pd, 2), pd, blk%nshell, blk%maxnp, blk%sh_l, blk%sh_np, blk%sh_ao, &
                                   blk%sh_at, blk%sh_e, blk%sh_c, blk%sh_r, b%natm, kmode, n, np, tens, jc, &
                                   naux, p0 - 1, cq, dens, grad)
            call tbl%step_off(lp, ST_D3C, tc)
            !$acc exit data delete(tens)
            deallocate (tens)
            call blk%release()
            s0 = s1 + 1
         end do
         !$acc exit data delete(xp)
         deallocate (xp)
      end do
      !$acc exit data delete(q_a, q_b, q_o, q_n, pd)
      deallocate (q_a, q_b, q_o, q_n, pd)
      if (do_k) then
         ! Z's upper triangle from its lower one (the panels above wrote the
         ! lower triangle and each diagonal block): reads only row > column,
         ! writes only row < column.
         do concurrent(t=0:int(naux, int64)**2 - 1) local(p, q)
            p = int(mod(t, int(naux, int64))) + 1
            q = int(t/int(naux, int64)) + 1
            if (p < q) zmat(p, q) = zmat(q, p)
         end do
      end if

      ! --- 3b. the two-centre derivative, once Z is complete ---------------
      call tbl%step_on(lp, ST_D2C, tc)
      call df_grad_2c_driver(aux%nshell, naux, aux%maxnp, aux%sh_l, aux%sh_np, aux%sh_ao, aux%sh_at, &
                             aux%sh_e, aux%sh_c, aux%sh_r, b%natm, zmat, grad)
      call tbl%step_off(lp, ST_D2C, tc)

      if (do_k) then
         !$acc exit data delete(cocc)
         deallocate (cocc)
      end if
      !$acc exit data delete(dpk, cq, zmat)
      deallocate (dpk, cq, zmat)
      dense = .false.
      dense([ST_CSOL, ST_HGEMM, ST_XTRSM, ST_ZGEMM, ST_YK]) = .true.
      call lp%sync()
      call tbl%report(tall%get_elapsed_time(), 0, 1, subtotal_mask=dense, subtotal_label="all dense algebra")
   end subroutine df_two_electron_gradient

   subroutine next_block(aux, s0, ablk, s1, p0, np)
      !! The run of whole auxiliary shells starting at s0 whose functions fit
      !! in `ablk` (at least one shell): last shell s1, first function p0,
      !! function count np. Same rule as trc_scf_df's setup blocks.
      type(trc_basis_t), intent(in) :: aux
      integer, intent(in) :: s0, ablk
      integer, intent(out) :: s1, p0, np
      integer :: ish
      s1 = s0
      do ish = s0 + 1, aux%nshell
         if (aux%sh_ao(ish) + ncart(aux%sh_l(ish)) - aux%sh_ao(s0) > ablk) exit
         s1 = ish
      end do
      p0 = aux%sh_ao(s0)
      np = aux%sh_ao(s1) + ncart(aux%sh_l(s1)) - p0
   end subroutine next_block

   subroutine choose_grad_blocks(n, naux, nocc, ablk, nci)
      !! The gradient's two block sizes, from the device's free memory AFTER
      !! everything resident is allocated (B, L, Z, C_occ, the caller's own
      !! arrays), half a gigabyte kept back for the libraries:
      !!
      !!   ablk  auxiliary functions per block, up to naux/4: each costs a
      !!         dense n x n slab, an (nocc x n) half-transform and an
      !!         (nocc x nocc) block of H or V;
      !!   nci   occupied orbitals per chunk of X (the packed rows j <= i of
      !!         the chunk's i, at most nci nocc naux doubles).
      !!
      !! When all of packed X (nocc(nocc+1)/2 naux doubles) fits next to a
      !! block of at least MIN_BLOCK functions -- the usual case -- X is
      !! reserved first and the blocks get everything else: every block is
      !! a full pass of the three-centre kernels over all the pairs, and a
      !! small one starves the GPU (w40/cc-pVDZ on a V100: the old even
      !! split of the free memory gave 25 blocks and 0.45 s of kernels, this
      !! rule ~15 blocks and 0.35 s, one pass over everything 0.18 s).
      !! Otherwise, the old split: blocks within half of what is free, X
      !! chunked in the rest.
      !!
      !! `nocc` = 0 (no exchange: no X at all) budgets the slab alone. A
      !! system for which not even one auxiliary function and one occupied
      !! orbital fit stops here, with the sizes, rather than at an
      !! allocation. TERCO_DFG_OCC_CHUNK=<k> caps nci at k (TESTING: forces
      !! the chunked path on a small molecule). The host build has no
      !! device to budget: naux/4 and every orbital, unless capped.
      integer, intent(in) :: n, naux, nocc
      integer, intent(out) :: ablk, nci
      integer(int64) :: per_q, per_i
      character(len=16) :: v
      integer :: st, cap
#ifdef _OPENACC
      integer, parameter :: MIN_BLOCK = 256
      integer(int64) :: free, avail, xall
      character(len=160) :: msg
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
      avail = free - 512_int64*1024*1024
#endif
      per_q = 8_int64*(int(n, int64)*(int(n, int64) + int(max(nocc, 1), int64)) + int(max(nocc, 1), int64)**2)
      per_i = 8_int64*int(max(nocc, 1), int64)*int(naux, int64)
      ablk = max(naux/4, 1)
      nci = max(nocc, 1)
#ifdef _OPENACC
      xall = 8_int64*int(naux, int64)*(int(nocc, int64)*int(nocc + 1, int64)/2)
      if (xall + int(min(MIN_BLOCK, naux), int64)*per_q <= avail) then
         ablk = int(min(int(ablk, int64), (avail - xall)/per_q))
      else
         ablk = int(min(int(ablk, int64), max(0_int64, (avail/2)/per_q)))
         if (nocc > 0) nci = int(min(int(nci, int64), max(0_int64, (avail - ablk*per_q)/per_i)))
      end if
      if (ablk < 1 .or. nci < 1) then
         write (msg, "(a,i0,a,i0,a,i0,a)") "trc_scf_grad_df: not enough device memory for the --df gradient (", &
            avail/1024/1024, " MiB free; one auxiliary function needs ", per_q/1024/1024, &
            " MiB, one occupied orbital's share of X ", per_i/1024/1024, " MiB)"
         error stop msg
      end if
#endif
      call get_environment_variable("TERCO_DFG_OCC_CHUNK", v, status=st)
      if (st == 0) then
         read (v, *, iostat=st) cap
         if (st == 0 .and. cap >= 1) nci = min(nci, cap)
      end if
      ablk = max(1, min(ablk, naux))
   end subroutine choose_grad_blocks

   pure integer function ncart(l) result(res)
      integer, intent(in) :: l
      res = (l + 1)*(l + 2)/2
   end function ncart

   pure integer function pack_idx(mu, nu) result(p)
      !! trc_scf_df's packed index, mu >= nu: mu(mu-1)/2 + nu.
      !$acc routine seq
      integer, intent(in) :: mu, nu
      p = mu*(mu - 1)/2 + nu
   end function pack_idx

end module trc_scf_grad_df
