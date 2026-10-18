!
! RI-MP2 for a closed-shell reference, on the device.
!
! With the auxiliary basis {P} and its metric J_PQ = (P|Q) = L L^T,
!
!     B_ia^Q = sum_P (ia|P) (L^-T)_PQ        so that  (ia|jb) = sum_Q B_ia^Q B_jb^Q
!
! and the correlation energy splits into its opposite- and same-spin parts,
!
!     E_os = sum_ijab (ia|jb)^2 / D_ijab
!     E_ss = sum_ijab [(ia|jb)^2 - (ia|jb)(ib|ja)] / D_ijab,   D = e_i + e_j - e_a - e_b
!
! E_MP2 = E_os + E_ss; SCS and SOS are weights on the two and cost nothing.
!
! THE STEPS
! ---------
!   1. (P|Q) by trc_df_2c, Cholesky on the device (trc_linalg%potrf).
!   2. (mu nu|P) by trc_df_3c, a block of auxiliary shells at a time so the
!      tensor never exceeds `aux_block` functions deep; each block is an
!      auxiliary basis of its own, built from a shell range, and the
!      kernels see nothing unusual.
!   3. (ia|P) = C_occ^T (mu nu|P) C_vir for this rank's occupied block
!      only, into X(ia, P) with a fastest: one GEMM and one strided batch
!      per auxiliary block.
!   4. B = X L^-T: one triangular solve from the right, on this rank's rows.
!   5. K(ia, jb) = sum_Q B(ia, Q) B(jb, Q) for the pairs this rank owns, a
!      group of occupied i per GEMM; then a kernel over (j, a) summing b
!      with the denominators. Both (ia|jb) and (ib|ja) are in K. The pair
!      energy is symmetric -- swapping i with j and a with b leaves
!      (ia|jb)^2/D and (ia|jb)(ib|ja)/D as they were -- so each unordered
!      pair is computed once and an i /= j pair counted twice, which takes
!      the GEMMs, nearly all of RI-MP2, to about half the work.
!
! RANKS: B IS DISTRIBUTED
! -----------------------
! The correlated occupied orbitals are split into one contiguous block per
! rank, O_r, and rank r holds only B_r, the rows ia with i in O_r: nocc nvir
! naux / nranks doubles instead of all of them, which is what bounds the
! system on one device. The integrals of step 2 are computed on every rank
! (they are the cheap part, and the transform needs every P of them); the
! transform and the solve are each rank's own rows.
!
! The pairs then come from a ring. At step s = 0 .. nranks/2 rank r holds
! the block of rank r + s (its own at s = 0, then each step passes the held
! block one rank to the left) and computes the pairs (i in O_r, j in
! O_{r+s}):
!   s = 0             j >= i inside its own block;
!   0 < s < n/2       every pair of the two blocks;
!   s = n/2, n even   the blocks r and r + n/2 meet twice, once on each
!                     rank, so the lower rank takes the first half of its
!                     own rows and the higher rank the second half of the
!                     visiting block's: each pair once.
! Every unordered pair (i, j) is then computed by exactly one rank, and
! E_os and E_ss are summed over the ranks at the end. On one rank only
! s = 0 runs, which is the single-device algorithm.
!
! A CUDA-aware MPI is handed the device buffers directly, which between
! GPUs joined by NVLink moves a block in a fraction of the time; any other
! MPI gets them through the host (trc_mpi_device decides, at run time). This
! is the distribution of Stocks, Palethorpe and Barca, JCTC 2024, 20, 7503.
!
! Steps 3-5 are cuBLAS from OpenACC through trc_linalg, pic-blas and LAPACK
! on the host; the two kernels are `do concurrent`.
!
! Module trc_rimp2_driver, not trc_rimp2: the C entry is bound to the
! symbol trc_rimp2, and a binding label may not equal a module name
! (F2018 19.2); gfortran silently resolves the call to the wrong thing.
module trc_rimp2_driver
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_2c, trc_df_3c
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t
   use pic_mpi_lib, only: comm_t, allreduce, bcast, isend, recv, request_t, waitall, MPI_SUM, &
                          comm_isend_real_dp_array_n, comm_recv_real_dp_array_n
   use trc_mpi_device, only: trc_mpi_device_direct
   use pic_flop_recorder, only: flop_recorder_type
   use trc_step_table, only: step_table_t
   use trc_scf_df, only: df_jk_t
   implicit none
   private

   public :: trc_rimp2_run, trc_rimp2_result_t

   !
   ! THE STEPS, AS TIMED, through the shared trc_step_table (src/core/
   ! trc_step_table.F90): each is device time to completion -- the clock is
   ! read after a device sync, since cuBLAS returns before its kernels end --
   ! and the floating-point operations it executes, so verbose output can say
   ! how close each step runs to the device's peak. The integral steps are
   ! time only: their operation counts are not GEMM-shaped and depend on the
   ! screening. Each step is also an NVTX range of the same name (a no-op
   ! unless built with TERCO_NVTX), opened and closed at the same device
   ! syncs as the timer, so on an Nsight Systems timeline a range covers its
   ! kernels and not only their launches.
   !
   integer, parameter :: NSTEP = 8
   integer, parameter :: ST_METRIC = 1, ST_CHOL = 2, ST_3C = 3, ST_TRANS = 4, ST_TRSM = 5, &
                         ST_KGEMM = 6, ST_PAIR = 7, ST_RING = 8
   character(len=18), parameter :: STEP_NAME(NSTEP) = [character(len=18) :: &
                                                       "metric (P|Q)", "Cholesky", "(mu nu|P)", &
                                                       "transform (ia|P)", "B = X L^-T", "K_i GEMMs", &
                                                       "pair energies", "B exchange"]

   !
   ! The largest message the ring sends in one call: MPI counts are 32-bit,
   ! and a block of B passes 2^31 doubles near w100 cc-pVDZ.
   !
   integer(int64), parameter :: MSG_CHUNK = 2_int64**30

   type :: trc_rimp2_result_t
      real(dp) :: e_os = 0.0_dp     !! opposite-spin pair energy
      real(dp) :: e_ss = 0.0_dp     !! same-spin pair energy
      real(dp) :: e_corr = 0.0_dp   !! E_os + E_ss
      integer :: naux = 0, nocc = 0, nvir = 0
      character(len=200) :: message = ""
   end type trc_rimp2_result_t

contains

   !
   ! `cmo` (nao, nao) and `eps` (nao) are the RHF orbitals, occupied first;
   ! `nocc` doubly occupied, `nfrozen` of them left out of the correlation.
   ! `pl` is the orbital basis's pair list. `la` optional, made here if
   ! absent; `comm` optional, one rank without it.
   !
   ! `dfjk`, optional: a density-fitted J/K object (trc_scf_df.F90) already
   ! fitted on an auxiliary basis the CALLER has checked is the same one
   ! as `aux` (same JSON path/contents -- trc_capi.F90's trc_run_rimp2 is
   ! the one caller that does; this routine only checks the dimensions
   ! agree, not that the basis itself is really the same one). When
   ! present, steps 1 (metric, Cholesky), 2 (the (mu nu|P) build) and 4
   ! (B = X L^-T) are skipped entirely: `dfjk%packed_b()` is already B
   ! (packed mu >= nu, L^-T applied), and the occupied/virtual half-
   ! transform of step 3 run on it directly IS RI-MP2's own B -- L^-T acts
   ! only on the auxiliary index, the half-transform only on (mu, nu), and
   ! the two commute. `aux` is still used for its shell structure (the
   ! block loop below only needs auxiliary function COUNTS and offsets,
   ! identical either way).
   !
   subroutine trc_rimp2_run(b, aux, pl, nocc, cmo, eps, res, nfrozen, aux_block, la, comm, verbose, dfjk)
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(b%nao, b%nao), eps(b%nao)
      type(trc_rimp2_result_t), intent(out) :: res
      integer, intent(in), optional :: nfrozen, aux_block
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(comm_t), intent(in), optional :: comm
      logical, intent(in), optional :: verbose   !! one line per step with its wall time, rank 0
      type(df_jk_t), intent(in), optional, target :: dfjk

      real(dp) :: tw0, tw1
      logical :: talk, chained
      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      type(trc_basis_t) :: blk
      real(dp), pointer :: bchain(:, :)
      integer :: nao, naux, nfz, no, nv, nov, nblk, s0, s1, p0, np, ish, i, rank, nranks, ng, ldk, i0, i1
      integer :: lo, nb, nmax, novr, step, nstep_ring, qv, cur, ncolmax, ldv
      integer :: ilo, ihi, jlo, jhi, jbase, jstart, ncol, roff, coff, jfirst, ncol_i
      integer :: n_chk, naux_chk, npack_chk, mu, nu, q
      integer, allocatable :: blo(:), bn(:)
      logical :: diag, direct
      real(dp), allocatable :: jm(:, :), cocc(:, :), cvir(:, :), tens(:, :, :), half(:, :, :)
      real(dp), allocatable :: xr(:), vis(:, :)
      real(dp), allocatable :: ki(:, :), kt(:, :), part(:), eo(:), ev(:)
      integer(int64) :: novr8
      real(dp) :: e_os, e_ss
      real(dp) :: es(2)
      real(dp) :: gf_all(1)
      real(dp) :: tc
      type(step_table_t) :: stab
      logical :: dense_mask(NSTEP)
      integer(int64) :: nv8, nao8, naux8

      nao = b%nao
      naux = aux%nao
      ! Density fitting stays Cartesian-only for this milestone (see
      ! trc_capi.F90's AO classification header, and its own refusal on
      ! the same basis at the C boundary): `cmo` above is declared
      ! (b%nao, b%nao), and a caller whose SCF ran spherical hands this a
      ! `cmo` actually shaped (b%nao_sph, b%nao_sph) -- smaller, so this
      ! would read past the end of the caller's array under sequence
      ! association rather than merely giving a wrong number. Refused
      ! before that read happens, not after.
      if (b%spherical) then
         res%message = "trc_rimp2: the orbital basis is spherical; RI-MP2 only "// &
                       "supports a Cartesian basis today"
         return
      end if
      nfz = 0
      if (present(nfrozen)) nfz = max(nfrozen, 0)
      no = nocc - nfz
      nv = nao - nocc
      nov = no*nv
      res%naux = naux
      res%nocc = no
      res%nvir = nv
      if (no < 1 .or. nv < 1) then
         res%message = "trc_rimp2: nothing to correlate"
         return
      end if
      if (maxval(b%sh_l) > DF_LMAX) then
         res%message = "trc_rimp2: the orbital basis carries angular momentum above what the "// &
                       "density-fitting kernels were generated for (see gen_df.py --lmax)"
         return
      end if
      if (maxval(aux%sh_l) > DF_LMAX_AUX) then
         res%message = "trc_rimp2: the auxiliary basis carries angular momentum above what the "// &
                       "density-fitting kernels were generated for (see gen_df.py --lmax-aux)"
         return
      end if
      chained = .false.
      bchain => null()
      if (present(dfjk)) then
         if (dfjk%is_ready()) then
            call dfjk%dims(n_chk, naux_chk, npack_chk)
            if (n_chk /= nao .or. naux_chk /= naux) then
               res%message = "trc_rimp2: the pre-fitted B does not match this orbital/auxiliary basis"
               return
            end if
            bchain => dfjk%packed_b()
            chained = .true.
         end if
      end if
      nblk = naux
      if (present(aux_block)) nblk = max(aux_block, 1)
      rank = 0
      nranks = 1
      if (present(comm)) then
         rank = comm%rank()
         nranks = comm%size()
      end if
      if (present(la)) then
         lp => la
      else
         lp => la_own
      end if
      ! This rank's occupied block, O_rank = lo .. lo + nb - 1, and every
      ! rank's, which the ring needs to size what it receives.
      allocate (blo(0:nranks - 1), bn(0:nranks - 1))
      do i = 0, nranks - 1
         bn(i) = no/nranks + merge(1, 0, i < mod(no, nranks))
         blo(i) = i*(no/nranks) + min(i, mod(no, nranks)) + 1
      end do
      lo = blo(rank)
      nb = bn(rank)
      nmax = maxval(bn)
      novr = nb*nv
      novr8 = novr
      ! Create-once: a caller's handles are kept, the workspace only grows.
      call lp%init(naux)

      talk = .false.
      if (present(verbose)) talk = verbose .and. rank == 0
      nv8 = nv
      nao8 = nao
      naux8 = naux
      call stab%init("rimp2", STEP_NAME)
      tw0 = tick(lp)
      ! --- 1. the metric and its Cholesky factor (skipped when chained) ----
      if (.not. chained) then
         allocate (jm(naux, naux))
         !$acc enter data create(jm)
         call stab%step_on(lp, ST_METRIC, tc)
         call trc_df_2c(aux, jm)
         call stab%step_off(lp, ST_METRIC, tc)
         call stab%step_on(lp, ST_CHOL, tc)
         call lp%potrf(naux, jm, naux)
         call stab%step_off(lp, ST_CHOL, tc)
         call stab%add_flop(ST_CHOL, naux8**3/3)
         if (talk) print "(a,i0,a,f8.2,a)", "  rimp2: metric (", naux, ") and Cholesky ", tick(lp) - tw0, " s"
      else if (talk) then
         print "(a)", "  rimp2: metric, Cholesky and (mu nu|P) skipped -- chained off --df's fit"
      end if

      ! --- 3. (ia|P) for i in this rank's block, an auxiliary block at a time
      allocate (cocc(nao, max(nb, 1)), cvir(nao, nv), xr(max(novr8*naux8, 1_int64)))
      if (nb > 0) cocc(:, 1:nb) = cmo(:, nfz + lo:nfz + lo + nb - 1)
      cvir = cmo(:, nocc + 1:nao)
      !$acc enter data copyin(cocc, cvir) create(xr)
      s0 = 1
      do while (s0 <= aux%nshell .and. nb > 0)
         ! Shells s0..s1 hold at most nblk functions (and at least one shell).
         s1 = s0
         do ish = s0 + 1, aux%nshell
            if (aux%sh_ao(ish) + ncart(aux%sh_l(ish)) - aux%sh_ao(s0) > nblk) exit
            s1 = ish
         end do
         p0 = aux%sh_ao(s0)
         np = aux%sh_ao(s1) + ncart(aux%sh_l(s1)) - p0
         call aux%subset(s0, s1, blk)
         call blk%to_device()
         allocate (tens(nao, nao, np), half(nb, nao, np))
         !$acc enter data create(tens, half)
         call stab%step_on(lp, ST_3C, tc)
         if (chained) then
            ! B, already L^-T applied, packed mu >= nu: unpack this block's
            ! columns straight into the dense tensor the half-transform
            ! below reads either way -- no (mu nu|P) build at all.
            do concurrent(q=1:np, mu=1:nao, nu=1:nao)
               tens(mu, nu, q) = bchain(rimp2_pack_idx(max(mu, nu), min(mu, nu)), p0 - 1 + q)
            end do
         else
            call trc_df_3c(b, pl, blk, tens)
         end if
         call stab%step_off(lp, ST_3C, tc)
         call stab%step_on(lp, ST_TRANS, tc)
         ! Two library calls per block, not two per auxiliary function: the
         ! occupied side of the whole block is one GEMM,
         !   half(i, nu q) = sum_mu C_occ(mu, i) T(mu, nu q)   (nb x nao np),
         ! and the virtual side one strided batch over q with C_vir shared,
         !   X(:, P) = C_vir^T half_q^T                        (nv x nb, a fastest).
         ! 2 x 3840 calls at w40 cc-pVDZ ran at 41% of the device's peak.
         call lp%gemm("T", "N", nb, nao*np, nao, 1.0_dp, cocc, nao, tens, nao, 0.0_dp, half, nb)
         call lp%gemm_strided("T", "T", nv, nb, nao, 1.0_dp, cvir, nao, 0_int64, half, nb, &
                              int(nb, int64)*nao8, 0.0_dp, xr(1 + int(p0 - 1, int64)*novr8), nv, novr8, np)
         call stab%step_off(lp, ST_TRANS, tc)
         call stab%add_flop(ST_TRANS, int(np, int64)*(2*int(nb, int64)*nao8*nao8 + 2*nv8*int(nb, int64)*nao8))
         !$acc exit data delete(tens, half)
         deallocate (tens, half)
         call blk%release()
         s0 = s1 + 1
      end do
      if (talk) print "(a,f8.2,a,f8.2,a)", "  rimp2: (mu nu|P) ", stab%t(ST_3C), " s, transform to (ia|P) ", &
         stab%t(ST_TRANS), " s"

      ! --- 4. B = X L^-T: solve B L^T = X from the right, this rank's rows ---
      ! (skipped when chained: `xr` is already built from B, which has
      ! L^-T applied already -- see this routine's own `dfjk` comment)
      if (.not. chained) then
         call stab%step_on(lp, ST_TRSM, tc)
         if (nb > 0) call lp%trsm("R", "L", "T", novr, naux, 1.0_dp, jm, naux, xr, novr)
         call stab%step_off(lp, ST_TRSM, tc)
         call stab%add_flop(ST_TRSM, novr8*naux8*naux8)
         if (talk) print "(a,f8.2,a)", "  rimp2: B = X L^-T ", stab%t(ST_TRSM), " s"
      end if
      tw1 = tick(lp)

      ! --- 5. the pair energies, over the ring ------------------------------
      ! The visiting blocks need two buffers, the one in use and the one
      ! arriving; one rank needs neither.
      nstep_ring = nranks/2
      direct = trc_mpi_device_direct()
      if (talk .and. nranks > 1) print "(a,a)", "  rimp2: B exchange ", &
         trim(merge("device to device (CUDA-aware MPI)", "through the host                 ", direct))
      if (nranks > 1) then
         allocate (vis(int(nmax, int64)*nv8*naux8, 2))
         !$acc enter data create(vis)
      end if
      ncolmax = nmax*nv
      ng = k_group(nmax, nv, ncolmax)
      ldk = ng*nv
      if (talk) print "(a,i0,a,i0,a)", "  rimp2: K GEMMs over ", ng, " occupied orbitals at a time, ", &
         nranks, " rank(s)"
      allocate (ki(ldk, ncolmax), kt(ncolmax, nv), part(2*ncolmax), eo(no), ev(nv))
      eo = eps(nfz + 1:nocc)
      ev = eps(nocc + 1:nao)
      !$acc enter data create(ki, kt, part) copyin(eo, ev)
      e_os = 0.0_dp
      e_ss = 0.0_dp
      cur = 0
      do step = 0, nstep_ring
         qv = mod(rank + step, nranks)
         if (step > 0) then
            call stab%step_on(lp, ST_RING, tc)
            call ring_shift(comm, rank, nranks, xr, vis, cur, &
                            int(bn(mod(rank + step - 1, nranks)), int64)*nv8*naux8, &
                            int(bn(qv), int64)*nv8*naux8, step == 1, direct)
            call stab%step_off(lp, ST_RING, tc)
         end if
         ! Rows i in ilo..ihi of this rank's block against columns j in
         ! jlo..jhi of block qv; `diag` keeps j >= i.
         diag = step == 0
         ilo = lo
         ihi = lo + nb - 1
         jbase = blo(qv)
         jlo = jbase
         jhi = jbase + bn(qv) - 1
         if (2*step == nranks) then
            if (rank < qv) then
               ihi = lo + nb/2 - 1
            else
               jlo = jbase + bn(qv)/2
            end if
         end if
         if (step == 0) then
            ldv = novr
         else
            ldv = bn(qv)*nv
         end if
         ! Groups of `ng` occupied orbitals: one GEMM per group, K(i a, jb)
         ! for i in the group and j from the group's first column on, so
         ! cuBLAS gets (ng nv) rows instead of nv. On the diagonal block it
         ! also computes the j < i corner inside the group, which the pair
         ! kernel skips.
         do i0 = ilo, ihi, ng
            i1 = min(ihi, i0 + ng - 1)
            jstart = jlo
            if (diag) jstart = i0
            ncol = (jhi - jstart + 1)*nv
            if (ncol <= 0) cycle
            call stab%step_on(lp, ST_KGEMM, tc)
            if (step == 0) then
               call lp%gemm("N", "T", (i1 - i0 + 1)*nv, ncol, naux, 1.0_dp, xr(1 + (i0 - lo)*nv), novr, &
                            xr(1 + (jstart - jbase)*nv), ldv, 0.0_dp, ki, ldk)
            else
               call lp%gemm("N", "T", (i1 - i0 + 1)*nv, ncol, naux, 1.0_dp, xr(1 + (i0 - lo)*nv), novr, &
                            vis(1 + (jstart - jbase)*nv, cur), ldv, 0.0_dp, ki, ldk)
            end if
            call stab%step_off(lp, ST_KGEMM, tc)
            call stab%add_flop(ST_KGEMM, 2*int(i1 - i0 + 1, int64)*nv8*int(ncol, int64)*naux8)
            do i = i0, i1
               ! K_i(a, jb) is the block of rows roff + a and, on the
               ! diagonal, columns from j = i; off it, every column.
               roff = (i - i0)*nv
               if (diag) then
                  coff = (i - i0)*nv
                  jfirst = i
               else
                  coff = 0
                  jfirst = jstart
               end if
               ncol_i = (jhi - jfirst + 1)*nv
               call stab%step_on(lp, ST_PAIR, tc)
               ! K_i^T, so the kernel reads (ia|jb) and (ib|ja) both coalesced
               call lp%transpose(nv, ncol_i, ki(roff + 1, coff + 1), ldk, kt, ncolmax)
               call pair_energy(i, jfirst, no, nv, ldk, roff, coff, ncolmax, ncol_i, eo, ev, ki, kt, part)
               !$acc update self(part(1:2*ncol_i))
               call stab%step_off(lp, ST_PAIR, tc)
               call stab%add_flop(ST_PAIR, 10*nv8*int(ncol_i, int64))
               e_os = e_os + sum(part(1:2*ncol_i:2))
               e_ss = e_ss + sum(part(2:2*ncol_i:2))
            end do
         end do
      end do
      !$acc exit data delete(ki, kt, part, eo, ev, xr, cocc, cvir)
      if (.not. chained) then
         !$acc exit data delete(jm)
      end if
      if (nranks > 1) then
         !$acc exit data delete(vis)
      end if
      if (nranks > 1) then
         es = [e_os, e_ss]
         call allreduce(comm, es, op=MPI_SUM)
         call bcast(comm, es, 2, 0)
         e_os = es(1)
         e_ss = es(2)
      end if
      ! Every rank's operations, for the rate of the whole job: the table
      ! below is rank 0's own share.
      gf_all(1) = 0.0_dp
      do i = 1, NSTEP
         gf_all(1) = gf_all(1) + real(stab%fl(i)%get(), dp)*1.0e-9_dp
      end do
      if (nranks > 1) call allreduce(comm, gf_all, op=MPI_SUM)
      if (talk) then
         print "(a,f8.2,a,f8.2,a)", "  rimp2: pair energies ", tick(lp) - tw1, " s;  total ", tick(lp) - tw0, " s"
         dense_mask = .false.
         dense_mask(ST_CHOL) = .true.
         dense_mask(ST_TRANS) = .true.
         dense_mask(ST_TRSM) = .true.
         dense_mask(ST_KGEMM) = .true.
         call stab%report(tick(lp) - tw0, rank, nranks, subtotal_mask=dense_mask, subtotal_label="all dense algebra")
         if (nranks > 1) print "(a,a24,f12.3,f13.1,f12.1)", "  rimp2:   ", "all ranks", tick(lp) - tw0, &
            gf_all(1), gf_all(1)/(tick(lp) - tw0)
      end if
      res%e_os = e_os
      res%e_ss = e_ss
      res%e_corr = e_os + e_ss
      if (.not. present(la)) call la_own%release()
   end subroutine trc_rimp2_run

   function wall() result(t)
      real(dp) :: t
      integer(int64) :: cc, rate
      call system_clock(cc, rate)
      t = real(cc, dp)/real(rate, dp)
   end function wall

   !
   ! The clock after everything queued on the device has finished.
   !
   function tick(lp) result(t)
      type(trc_linalg_t), intent(in) :: lp
      real(dp) :: t
      call lp%sync()
      t = wall()
   end function tick

   !
   ! Occupied orbitals per K GEMM. One i gives cuBLAS nvir rows, which it
   ! tiles small; four give it a matrix tall enough for its larger tiles.
   ! Measured at w40 cc-pVDZ on a V100, K GEMM seconds for 1, 2, 4, 8, 16
   ! per group: 9.37, 9.09, 8.78, 9.24, 9.38 -- four is the best, and
   ! larger groups lose to the j < i corner each one computes and throws
   ! away. K is that many times larger, so on the device the group shrinks
   ! to what the free memory holds next to K^T and the energies, keeping
   ! 512 MB back for the libraries' own workspace.
   !
   !
   ! One step of the ring: pass the block this rank holds to the rank on its
   ! left and take the one on its right, through the host. At the first
   ! step the held block is this rank's own (`xr`), after that the visiting
   ! buffer `cur`, and what arrives goes into the other buffer, which
   ! becomes `cur`. Blocking calls in an order that cannot deadlock: even
   ! ranks send first and odd ranks receive first, so every send meets a
   ! receive that is already posted or about to be. In pieces of at most
   ! MSG_CHUNK doubles, since MPI counts are 32-bit.
   !
   subroutine ring_shift(comm, rank, nranks, xr, vis, cur, nsend, nrecv, first, direct)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: rank, nranks
      real(dp), intent(inout) :: xr(:), vis(:, :)
      integer, intent(inout) :: cur
      integer(int64), intent(in) :: nsend, nrecv
      logical, intent(in) :: first, direct
      integer :: left, right, nxt, nreq
      type(request_t), allocatable :: reqs(:)
      left = mod(rank - 1 + nranks, nranks)
      right = mod(rank + 1, nranks)
      nxt = merge(1, 3 - cur, first)
      nreq = int((nsend + MSG_CHUNK - 1)/MSG_CHUNK)
      allocate (reqs(max(nreq, 1)))
      if (direct) then
         ! The device addresses straight to MPI: inside host_data, through
         ! pic-mpi's explicit-count procedures, which pass the buffer on
         ! without a Fortran temporary (a temporary would be a host copy of
         ! a device address).
         !$acc host_data use_device(xr, vis)
         call send_block_dev(comm, xr, vis, cur, nsend, left, first, reqs)
         call recv_block_dev(comm, vis, nxt, nrecv, right)
         !$acc end host_data
         if (nreq > 0) call waitall(reqs(1:nreq))
         cur = nxt
         return
      end if
      if (first) then
         !$acc update self(xr(1:nsend))
      else
         !$acc update self(vis(1:nsend, cur))
      end if
      ! The send is posted first and non-blocking, so both directions of
      ! every link move at once; the receive then blocks, and the send is
      ! waited for before its buffer can be reused.
      call send_block(comm, xr, vis, cur, nsend, left, first, reqs)
      call recv_block(comm, vis, nxt, nrecv, right)
      if (nreq > 0) call waitall(reqs(1:nreq))
      !$acc update device(vis(1:nrecv, nxt))
      cur = nxt
   end subroutine ring_shift

   !
   ! The held block to `dest`, from the host copy, in MSG_CHUNK pieces.
   ! Module procedures rather than internal ones: see app/main.F90's
   ! terco_run for what a sibling internal procedure did under NVHPC's MPI.
   !
   subroutine send_block(comm, xr, vis, cur, n, dest, first, reqs)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in), target :: xr(:), vis(:, :)
      integer, intent(in) :: cur, dest
      integer(int64), intent(in) :: n
      logical, intent(in) :: first
      type(request_t), intent(inout) :: reqs(:)
      integer(int64) :: k, m
      integer :: r
      k = 1
      r = 0
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         r = r + 1
         if (first) then
            call isend(comm, xr(k:k + m - 1), dest, 7, reqs(r))
         else
            call isend(comm, vis(k:k + m - 1, cur), dest, 7, reqs(r))
         end if
         k = k + m
      end do
   end subroutine send_block

   !
   ! The next block from `source` into buffer `c` of `vis`, in MSG_CHUNK pieces.
   !
   subroutine recv_block(comm, vis, c, n, source)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: vis(:, :)
      integer, intent(in) :: c, source
      integer(int64), intent(in) :: n
      integer(int64) :: k, m
      k = 1
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         call recv(comm, vis(k:k + m - 1, c), int(m), source, 7)
         k = k + m
      end do
   end subroutine recv_block

   !
   ! send_block and recv_block for device addresses, called inside
   ! `host_data use_device`: pic-mpi's explicit-count procedures by their
   ! specific names, since its generic interfaces misbehave there under
   ! nvfortran.
   !
   subroutine send_block_dev(comm, xr, vis, cur, n, dest, first, reqs)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in), target :: xr(:), vis(:, :)
      integer, intent(in) :: cur, dest
      integer(int64), intent(in) :: n
      logical, intent(in) :: first
      type(request_t), intent(inout) :: reqs(:)
      integer(int64) :: k, m
      integer :: r
      k = 1
      r = 0
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         r = r + 1
         if (first) then
            call comm_isend_real_dp_array_n(comm, xr(k:k + m - 1), int(m), dest, 7, reqs(r))
         else
            call comm_isend_real_dp_array_n(comm, vis(k:k + m - 1, cur), int(m), dest, 7, reqs(r))
         end if
         k = k + m
      end do
   end subroutine send_block_dev

   subroutine recv_block_dev(comm, vis, c, n, source)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout), target :: vis(:, :)
      integer, intent(in) :: c, source
      integer(int64), intent(in) :: n
      integer(int64) :: k, m
      k = 1
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         call comm_recv_real_dp_array_n(comm, vis(k:k + m - 1, c), int(m), source, 7)
         k = k + m
      end do
   end subroutine recv_block_dev

   integer function k_group(no, nv, ncol) result(ng)
#ifdef _OPENACC
      use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, &
                         acc_property_free_memory
#endif
      integer, intent(in) :: no, nv, ncol
      integer, parameter :: NG_BEST = 4
#ifdef _OPENACC
      integer(int64) :: free, per_i, fixed
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
      per_i = 8_int64*int(nv, int64)*int(ncol, int64)
      fixed = per_i + 16_int64*int(ncol, int64) + 512_int64*1024*1024
      ng = int(min(int(NG_BEST, int64), max(1_int64, (free - fixed)/per_i)))
#else
      ng = NG_BEST
#endif
      ng = max(1, min(ng, no))
   end function k_group

   pure integer function ncart(l) result(res)
      integer, intent(in) :: l
      res = (l + 1)*(l + 2)/2
   end function ncart

   !! idx(mu, nu) = mu(mu-1)/2 + nu, mu >= nu, 1-based: trc_scf_df.F90's own
   !! `pack_idx`, duplicated here rather than used across the module
   !! boundary -- a `do concurrent` body may not call a procedure from
   !! another module (the same reason `ncart` above is its own copy and
   !! not trc_scf_df.F90's).
   pure integer function rimp2_pack_idx(mu, nu) result(p)
      !$acc routine seq
      integer, intent(in) :: mu, nu
      p = mu*(mu - 1)/2 + nu
   end function rimp2_pack_idx

   !
   ! For occupied i and the pairs j = jfirst.., one thread per (j, a), with a
   ! fastest and b summed inside: part(2 ja - 1) = w sum_b (ia|jb)^2 / D and
   ! part(2 ja) = w sum_b [(ia|jb)^2 - (ia|jb)(ib|ja)] / D, w = 1 for j = i
   ! and 2 otherwise -- the (j, i) pair is the same numbers again and is
   ! never computed. No reduction crosses threads; the host adds 2 nji
   ! numbers. Diagonal terms are in the sum: (ia|ia)^2 - (ia|ia)^2 = 0 for
   ! the same-spin part, as it must.
   !
   ! (ia|jb) is K(roff + a, coff + (j - jfirst) nv + b) and (ib|ja) is
   ! K(roff + b, coff + (j - jfirst) nv + a) = K_i^T(ja, b): with a fastest
   ! across threads, both are unit-stride reads. From K alone one of the
   ! two is always strided by nv across a warp, which held this kernel to a
   ! fifth of the device's bandwidth.
   !
   subroutine pair_energy(i, jfirst, no, nv, ldk, roff, coff, ldt, nji, eo, ev, ki, kt, part)
      integer, intent(in) :: i, jfirst, no, nv, ldk, roff, coff, ldt, nji
      real(dp), intent(in) :: eo(no), ev(nv), ki(ldk, *), kt(ldt, nv)
      real(dp), intent(out) :: part(*)
      integer :: ja
      do concurrent(ja=1:nji)
         call pair_energy_body(i, jfirst, ja, no, nv, ldk, roff, coff, ldt, eo, ev, ki, kt, &
                               part(2*ja - 1), part(2*ja))
      end do
   end subroutine pair_energy

   pure subroutine pair_energy_body(i, jfirst, ja, no, nv, ldk, roff, coff, ldt, eo, ev, ki, kt, pos, pss)
      !$acc routine seq
      integer, intent(in) :: i, jfirst, ja, no, nv, ldk, roff, coff, ldt
      real(dp), intent(in) :: eo(no), ev(nv), ki(ldk, *), kt(ldt, nv)
      real(dp), intent(out) :: pos, pss
      integer :: jj, j, a, b
      real(dp) :: t, u, r, eija, w
      jj = (ja - 1)/nv               ! j - jfirst
      j = jfirst + jj
      a = ja - jj*nv
      eija = eo(i) + eo(j) - ev(a)
      w = merge(1.0_dp, 2.0_dp, j == i)
      pos = 0.0_dp
      pss = 0.0_dp
      do b = 1, nv
         t = ki(roff + a, coff + jj*nv + b)  ! (ia|jb)
         u = kt(ja, b)                       ! (ib|ja)
         r = 1.0_dp/(eija - ev(b))
         pos = pos + t*t*r
         pss = pss + (t*t - t*u)*r
      end do
      pos = w*pos
      pss = w*pss
   end subroutine pair_energy_body

end module trc_rimp2_driver
