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
!   3. (ia|P) = C_occ^T (mu nu|P) C_vir: two GEMMs per P, into X(ia, P)
!      with a fastest.
!   4. B = X L^-T: one triangular solve from the right.
!   5. Per occupied i, K_i(a, jb) = sum_Q B(ia, Q) B(jb, Q) for j >= i only:
!      one GEMM (nvir x naux) (naux x (nocc - i + 1) nvir); then a kernel
!      over (j, b) summing a with the denominators, an i /= j pair counted
!      twice. Both (ia|jb) and (ib|ja) are in K_i. The pair energy is
!      symmetric -- swapping i with j and a with b leaves (ia|jb)^2/D and
!      (ia|jb)(ib|ja)/D as they were -- so the j < i half is the j > i
!      half again, and skipping it takes the GEMMs, which are nearly all
!      of RI-MP2, to (nocc + 1)/(2 nocc) of the work.
!   6. Ranks: i dealt round-robin, E_os and E_ss allreduced. Everything
!      before step 5 is replicated -- it is the cheap part -- so the ranks
!      never communicate a tensor.
!
! Steps 3-5 are cuBLAS from OpenACC through trc_linalg, pic-blas and LAPACK
! on the host; the two kernels are `do concurrent`. This is the structure of
! Stocks, Palethorpe and Barca, JCTC 2024, 20, 7503, without their
! distribution of the transformed tensor across GPUs: the ranks here each
! hold the whole of B, which bounds the system size at what one device
! holds, nocc x nvir x naux doubles.
!
! Module trc_rimp2_driver, not trc_rimp2: the C entry is bound to the
! symbol trc_rimp2, and a binding label may not equal a module name
! (F2018 19.2); gfortran silently resolves the call to the wrong thing.
module trc_rimp2_driver
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_2c, trc_df_3c
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   use pic_flop_recorder, only: flop_recorder_type
   use pic_profiler, only: profiler_start, profiler_stop
   implicit none
   private

   public :: trc_rimp2_run, trc_rimp2_result_t

   !
   ! THE STEPS, AS TIMED. Each is device time to completion -- the clock is
   ! read after a device sync, since cuBLAS returns before its kernels end --
   ! and the floating-point operations it executes, so verbose output can say
   ! how close each step runs to the device's peak. The integral steps are
   ! time only: their operation counts are not GEMM-shaped and depend on the
   ! screening.
   !
   integer, parameter :: NSTEP = 7
   integer, parameter :: ST_METRIC = 1, ST_CHOL = 2, ST_3C = 3, ST_TRANS = 4, ST_TRSM = 5, &
                         ST_KGEMM = 6, ST_PAIR = 7
   !
   ! Each step is also an NVTX range of the same name (pic_profiler; a no-op
   ! unless built with TERCO_NVTX), opened and closed at the same device
   ! syncs as the timer, so on an Nsight Systems timeline a range covers its
   ! kernels and not only their launches.
   !
   character(len=18), parameter :: STEP_NAME(NSTEP) = [character(len=18) :: &
                                                       "metric (P|Q)", "Cholesky", "(mu nu|P)", &
                                                       "transform (ia|P)", "B = X L^-T", "K_i GEMMs", &
                                                       "pair energies"]

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
   subroutine trc_rimp2_run(b, aux, pl, nocc, cmo, eps, res, nfrozen, aux_block, la, comm, verbose)
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(b%nao, b%nao), eps(b%nao)
      type(trc_rimp2_result_t), intent(out) :: res
      integer, intent(in), optional :: nfrozen, aux_block
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(comm_t), intent(in), optional :: comm
      logical, intent(in), optional :: verbose   !! one line per step with its wall time, rank 0

      real(dp) :: tw0, tw1
      logical :: talk
      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      type(trc_basis_t) :: blk
      integer :: nao, naux, nfz, no, nv, nov, nblk, s0, s1, p0, np, ish, q, i, rank, nranks, nji, ng, ldk, ig, i0, i1, off
      real(dp), allocatable :: jm(:, :), cocc(:, :), cvir(:, :), x(:, :), tens(:, :, :), half(:, :, :)
      real(dp), allocatable :: ki(:, :), kt(:, :), part(:), eo(:), ev(:)
      real(dp) :: e_os, e_ss, es(2)
      real(dp) :: t_step(NSTEP), tc
      type(flop_recorder_type) :: fl(NSTEP)
      integer(int64) :: no8, nv8, nao8, naux8, nov8

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
      res%naux = naux; res%nocc = no; res%nvir = nv
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
      nblk = naux
      if (present(aux_block)) nblk = max(aux_block, 1)
      rank = 0; nranks = 1
      if (present(comm)) then
         rank = comm%rank(); nranks = comm%size()
      end if
      if (present(la)) then
         lp => la
      else
         lp => la_own
      end if
      ! Create-once: a caller's handles are kept, the workspace only grows.
      call lp%init(naux)

      talk = .false.
      if (present(verbose)) talk = verbose .and. rank == 0
      no8 = no; nv8 = nv; nao8 = nao; naux8 = naux; nov8 = nov
      t_step = 0.0_dp
      tw0 = tick(lp)
      ! --- 1. the metric and its Cholesky factor ---------------------------
      allocate (jm(naux, naux))
      !$acc enter data create(jm)
      call step_on(lp, ST_METRIC, tc)
      call trc_df_2c(aux, jm)
      call step_off(lp, ST_METRIC, tc, t_step)
      call step_on(lp, ST_CHOL, tc)
      call lp%potrf(naux, jm, naux)
      call step_off(lp, ST_CHOL, tc, t_step)
      call fl(ST_CHOL)%add(naux8**3/3)
      if (talk) print '(a,i0,a,f8.2,a)', "  rimp2: metric (", naux, ") and Cholesky ", tick(lp) - tw0, " s"

      ! --- 3. (ia|P) for every P, a block of auxiliary shells at a time ----
      allocate (cocc(nao, no), cvir(nao, nv), x(nov, naux))
      cocc = cmo(:, nfz + 1:nocc)
      cvir = cmo(:, nocc + 1:nao)
      !$acc enter data copyin(cocc, cvir) create(x)
      s0 = 1
      do while (s0 <= aux%nshell)
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
         allocate (tens(nao, nao, np), half(no, nao, np))
         !$acc enter data create(tens, half)
         call step_on(lp, ST_3C, tc)
         call trc_df_3c(b, pl, blk, tens)
         call step_off(lp, ST_3C, tc, t_step)
         call step_on(lp, ST_TRANS, tc)
         ! Two library calls per block, not two per auxiliary function: the
         ! occupied side of the whole block is one GEMM,
         !   half(i, nu q) = sum_mu C_occ(mu, i) T(mu, nu q)   (no x nao np),
         ! and the virtual side one strided batch over q with C_vir shared,
         !   X(:, P) = C_vir^T half_q^T                        (nv x no, a fastest).
         ! 2 x 3840 calls at w40 cc-pVDZ ran at 41% of the device's peak.
         call lp%gemm('T', 'N', no, nao*np, nao, 1.0_dp, cocc, nao, tens, nao, 0.0_dp, half, no)
         call lp%gemm_strided('T', 'T', nv, no, nao, 1.0_dp, cvir, nao, 0_int64, half, no, no8*nao8, &
                              0.0_dp, x(1, p0), nv, nov8, np)
         call step_off(lp, ST_TRANS, tc, t_step)
         call fl(ST_TRANS)%add(int(np, int64)*(2*no8*nao8*nao8 + 2*nv8*no8*nao8))
         !$acc exit data delete(tens, half)
         deallocate (tens, half)
         call blk%release()
         s0 = s1 + 1
      end do
      if (talk) print '(a,f8.2,a,f8.2,a)', "  rimp2: (mu nu|P) ", t_step(ST_3C), " s, transform to (ia|P) ", &
         t_step(ST_TRANS), " s"

      ! --- 4. B = X L^-T: solve B L^T = X from the right ---------------------
      call step_on(lp, ST_TRSM, tc)
      call lp%trsm('R', 'L', 'T', nov, naux, 1.0_dp, jm, naux, x, nov)
      call step_off(lp, ST_TRSM, tc, t_step)
      call fl(ST_TRSM)%add(nov8*naux8*naux8)
      if (talk) print '(a,f8.2,a)', "  rimp2: B = X L^-T ", t_step(ST_TRSM), " s"
      tw1 = tick(lp)

      ! --- 5. the pair energies, one occupied i at a time ------------------
      ng = k_group(no, nv, nov)
      ldk = ng*nv
      if (talk) print '(a,i0,a)', "  rimp2: K GEMMs over ", ng, " occupied orbitals at a time"
      allocate (ki(ldk, nov), kt(nov, nv), part(2*nov), eo(no), ev(nv))
      eo = eps(nfz + 1:nocc)
      ev = eps(nocc + 1:nao)
      !$acc enter data create(ki, kt, part) copyin(eo, ev)
      e_os = 0.0_dp; e_ss = 0.0_dp
      ! Groups of `ng` occupied orbitals, dealt round-robin over the ranks.
      ! One GEMM per group, K(i a, jb) for i in the group and j >= i0, the
      ! group's first: (ng nv) rows instead of nv, so cuBLAS gets a taller
      ! matrix. It also computes the j < i corner inside the group, which
      ! the pair kernel then skips -- (ng - 1)/(2 (no - i0 + 1)) of the group's
      ! work on average.
      do ig = 1 + rank, (no + ng - 1)/ng, nranks
         i0 = (ig - 1)*ng + 1
         i1 = min(no, i0 + ng - 1)
         nji = (no - i0 + 1)*nv
         call step_on(lp, ST_KGEMM, tc)
         call lp%gemm('N', 'T', (i1 - i0 + 1)*nv, nji, naux, 1.0_dp, x((i0 - 1)*nv + 1, 1), nov, &
                      x((i0 - 1)*nv + 1, 1), nov, 0.0_dp, ki, ldk)
         call step_off(lp, ST_KGEMM, tc, t_step)
         call fl(ST_KGEMM)%add(2*int(i1 - i0 + 1, int64)*nv8*int(nji, int64)*naux8)
         do i = i0, i1
            ! K_i(a, jb), j >= i, is the block of rows (i - i0) nv + a and
            ! columns from (i - i0) nv + 1.
            off = (i - i0)*nv
            nji = (no - i + 1)*nv
            call step_on(lp, ST_PAIR, tc)
            ! K_i^T, so the kernel reads (ia|jb) and (ib|ja) both coalesced
            call lp%transpose(nv, nji, ki(off + 1, off + 1), ldk, kt, nov)
            call pair_energy(i, no, nv, ldk, off, nov, nji, eo, ev, ki, kt, part)
            !$acc update self(part(1:2*nji))
            call step_off(lp, ST_PAIR, tc, t_step)
            call fl(ST_PAIR)%add(10*nv8*int(nji, int64))
            e_os = e_os + sum(part(1:2*nji:2))
            e_ss = e_ss + sum(part(2:2*nji:2))
         end do
      end do
      !$acc exit data delete(ki, kt, part, eo, ev, x, cocc, cvir, jm)
      if (nranks > 1) then
         es = [e_os, e_ss]
         call allreduce(comm, es, op=MPI_SUM)
         call bcast(comm, es, 2, 0)
         e_os = es(1); e_ss = es(2)
      end if
      if (talk) then
         print '(a,f8.2,a,f8.2,a)', "  rimp2: pair energies ", tick(lp) - tw1, " s;  total ", tick(lp) - tw0, " s"
         call report_steps(t_step, fl, tick(lp) - tw0, rank, nranks)
      end if
      res%e_os = e_os
      res%e_ss = e_ss
      res%e_corr = e_os + e_ss
      if (.not. present(la)) call la_own%release()
   end subroutine trc_rimp2_run

   function wall() result(t)
      real(dp) :: t
      integer(kind=8) :: cc, rate
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
   ! Open step `k`: its start time in `tc`, and its NVTX range.
   !
   subroutine step_on(lp, k, tc)
      type(trc_linalg_t), intent(in) :: lp
      integer, intent(in) :: k
      real(dp), intent(out) :: tc
      tc = tick(lp)
      call profiler_start("rimp2: "//trim(STEP_NAME(k)), nvtx_only=.true.)
   end subroutine step_on

   !
   ! Close step `k`: the time since `tc` added to its total, its range closed.
   !
   subroutine step_off(lp, k, tc, t_step)
      type(trc_linalg_t), intent(in) :: lp
      integer, intent(in) :: k
      real(dp), intent(in) :: tc
      real(dp), intent(inout) :: t_step(NSTEP)
      t_step(k) = t_step(k) + (tick(lp) - tc)
      call profiler_stop("rimp2: "//trim(STEP_NAME(k)))
   end subroutine step_off

   !
   ! One line per step: device time, operations executed and the rate, then
   ! the GEMM-shaped steps together. Rank 0's own share when there are
   ! several ranks: the pair steps are split over them, the rest replicated.
   !
   subroutine report_steps(t_step, fl, t_total, rank, nranks)
      real(dp), intent(in) :: t_step(NSTEP), t_total
      type(flop_recorder_type), intent(in) :: fl(NSTEP)
      integer, intent(in) :: rank, nranks
      integer :: k
      real(dp) :: gf, tg, fg
      print '(a)', "  rimp2: step                   time (s)        GFLOP     GFLOP/s"
      do k = 1, NSTEP
         gf = real(fl(k)%get(), dp)*1.0e-9_dp
         if (fl(k)%get() > 0 .and. t_step(k) > 0.0_dp) then
            print '(a,a18,f14.3,f13.1,f12.1)', "  rimp2:   ", STEP_NAME(k), t_step(k), gf, gf/t_step(k)
         else
            print '(a,a18,f14.3,a13,a12)', "  rimp2:   ", STEP_NAME(k), t_step(k), "-", "-"
         end if
      end do
      tg = t_step(ST_CHOL) + t_step(ST_TRANS) + t_step(ST_TRSM) + t_step(ST_KGEMM)
      fg = real(fl(ST_CHOL)%get() + fl(ST_TRANS)%get() + fl(ST_TRSM)%get() + fl(ST_KGEMM)%get(), dp)*1.0e-9_dp
      if (tg > 0.0_dp) print '(a,a18,f14.3,f13.1,f12.1)', "  rimp2:   ", label("all dense algebra"), tg, fg, fg/tg
      print '(a,a18,f14.3,a,i0,a,i0)', "  rimp2:   ", label("total"), t_total, "   rank ", rank, " of ", nranks
   end subroutine report_steps

   pure function label(t) result(l)
      !! `t` left-justified in the table's column: an `a18` edit descriptor
      !! right-justifies anything shorter than 18.
      character(len=*), intent(in) :: t
      character(len=18) :: l
      l = t
   end function label

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
   integer function k_group(no, nv, nov) result(ng)
#ifdef _OPENACC
      use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, &
                         acc_property_free_memory
#endif
      integer, intent(in) :: no, nv, nov
      integer, parameter :: NG_BEST = 4
#ifdef _OPENACC
      integer(int64) :: free, per_i, fixed
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
      per_i = 8_int64*int(nv, int64)*int(nov, int64)
      fixed = per_i + 16_int64*int(nov, int64) + 512_int64*1024*1024
      ng = int(min(int(NG_BEST, int64), max(1_int64, (free - fixed)/per_i)))
#else
      ng = NG_BEST
#endif
      ng = max(1, min(ng, no))
   end function k_group

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   !
   ! For occupied i and the pairs j >= i, one thread per (j, a), with a
   ! fastest and b summed inside: part(2 ja - 1) = w sum_b (ia|jb)^2 / D and
   ! part(2 ja) = w sum_b [(ia|jb)^2 - (ia|jb)(ib|ja)] / D, ja counting
   ! (j, a) from j = i, w = 2 for j > i and 1 for j = i -- the j < i pairs
   ! are the same numbers again. No reduction crosses threads; the host adds
   ! 2 nji numbers. Diagonal terms are in the sum: (ia|ia)^2 - (ia|ia)^2 = 0
   ! for the same-spin part, as it must.
   !
   ! (ia|jb) is K_i(a, jb) and (ib|ja) is K_i(b, ja) = K_i^T(ja, b): with a
   ! fastest across threads, both are unit-stride reads. From K_i alone one
   ! of the two is always strided by nv across a warp, which held this
   ! kernel to a fifth of the device's bandwidth.
   !
   subroutine pair_energy(i, no, nv, ldk, off, ldt, nji, eo, ev, ki, kt, part)
      integer, intent(in) :: i, no, nv, ldk, off, ldt, nji
      real(dp), intent(in) :: eo(no), ev(nv), ki(ldk, *), kt(ldt, nv)
      real(dp), intent(out) :: part(*)
      integer :: ja
      do concurrent(ja=1:nji)
         call pair_energy_body(i, ja, no, nv, ldk, off, ldt, eo, ev, ki, kt, part(2*ja - 1), part(2*ja))
      end do
   end subroutine pair_energy

   pure subroutine pair_energy_body(i, ja, no, nv, ldk, off, ldt, eo, ev, ki, kt, pos, pss)
      !$acc routine seq
      integer, intent(in) :: i, ja, no, nv, ldk, off, ldt
      real(dp), intent(in) :: eo(no), ev(nv), ki(ldk, *), kt(ldt, nv)
      real(dp), intent(out) :: pos, pss
      integer :: jj, j, a, b
      real(dp) :: t, u, r, eija, w
      jj = (ja - 1)/nv               ! j - i
      j = i + jj
      a = ja - jj*nv
      eija = eo(i) + eo(j) - ev(a)
      w = merge(1.0_dp, 2.0_dp, jj == 0)
      pos = 0.0_dp
      pss = 0.0_dp
      do b = 1, nv
         t = ki(off + a, off + jj*nv + b)    ! (ia|jb)
         u = kt(ja, b)                       ! (ib|ja)
         r = 1.0_dp/(eija - ev(b))
         pos = pos + t*t*r
         pss = pss + (t*t - t*u)*r
      end do
      pos = w*pos
      pss = w*pss
   end subroutine pair_energy_body

end module trc_rimp2_driver
