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
!   5. Per occupied i, K_i(a, jb) = sum_Q B(ia, Q) B(jb, Q): one GEMM
!      (nvir x naux) (naux x nocc nvir); then a kernel over (j, b) summing
!      a with the denominators. Both (ia|jb) and (ib|ja) are in K_i.
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
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_2c, trc_df_3c
   use trc_linalg, only: trc_linalg_t
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   implicit none
   private

   public :: trc_rimp2_run, trc_rimp2_result_t

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
   subroutine trc_rimp2_run(b, aux, pl, nocc, cmo, eps, res, nfrozen, aux_block, la, comm)
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(b%nao, b%nao), eps(b%nao)
      type(trc_rimp2_result_t), intent(out) :: res
      integer, intent(in), optional :: nfrozen, aux_block
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(comm_t), intent(in), optional :: comm

      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      type(trc_basis_t) :: blk
      integer :: nao, naux, nfz, no, nv, nov, nblk, s0, s1, p0, np, ish, q, i, rank, nranks
      real(dp), allocatable :: jm(:, :), cocc(:, :), cvir(:, :), x(:, :), tens(:, :, :), tmp(:, :)
      real(dp), allocatable :: ki(:, :), part(:), eo(:), ev(:)
      real(dp) :: e_os, e_ss, es(2)

      nao = b%nao
      naux = aux%nao
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
      nblk = naux
      if (present(aux_block)) nblk = max(aux_block, 1)
      rank = 0; nranks = 1
      if (present(comm)) then
         rank = comm%rank(); nranks = comm%size()
      end if
      if (present(la)) then
         lp => la
      else
         call la_own%init(naux)
         lp => la_own
      end if

      ! --- 1. the metric and its Cholesky factor ---------------------------
      allocate (jm(naux, naux))
      !$acc enter data create(jm)
      call trc_df_2c(aux, jm)
      call lp%potrf(naux, jm, naux)

      ! --- 3. (ia|P) for every P, a block of auxiliary shells at a time ----
      allocate (cocc(nao, no), cvir(nao, nv), x(nov, naux), tmp(no, nao))
      cocc = cmo(:, nfz + 1:nocc)
      cvir = cmo(:, nocc + 1:nao)
      !$acc enter data copyin(cocc, cvir) create(x, tmp)
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
         allocate (tens(nao, nao, np))
         !$acc enter data create(tens)
         call trc_df_3c(b, pl, blk, tens)
         do q = 1, np
            ! tmp = C_occ^T T_q  (no x nao), then X(:, P) = (C_vir^T tmp^T) = (nv x no), a fastest
            call lp%gemm('T', 'N', no, nao, nao, 1.0_dp, cocc, nao, tens(1, 1, q), nao, 0.0_dp, tmp, no)
            call lp%gemm('T', 'T', nv, no, nao, 1.0_dp, cvir, nao, tmp, no, 0.0_dp, x(1, p0 + q - 1), nv)
         end do
         !$acc exit data delete(tens)
         deallocate (tens)
         call blk%release()
         s0 = s1 + 1
      end do

      ! --- 4. B = X L^-T: solve B L^T = X from the right ---------------------
      call lp%trsm('R', 'L', 'T', nov, naux, 1.0_dp, jm, naux, x, nov)

      ! --- 5. the pair energies, one occupied i at a time ------------------
      allocate (ki(nv, nov), part(2*nov), eo(no), ev(nv))
      eo = eps(nfz + 1:nocc)
      ev = eps(nocc + 1:nao)
      !$acc enter data create(ki, part) copyin(eo, ev)
      e_os = 0.0_dp; e_ss = 0.0_dp
      do i = 1 + rank, no, nranks
         ! K_i(a, jb) = sum_Q B(ia, Q) B(jb, Q); rows ia of X start at (i-1) nv + 1
         call lp%gemm('N', 'T', nv, nov, naux, 1.0_dp, x((i - 1)*nv + 1, 1), nov, x, nov, 0.0_dp, ki, nv)
         call pair_energy(i, no, nv, nov, eo, ev, ki, part)
         !$acc update self(part)
         e_os = e_os + sum(part(1:2*nov:2))
         e_ss = e_ss + sum(part(2:2*nov:2))
      end do
      !$acc exit data delete(ki, part, eo, ev, x, tmp, cocc, cvir, jm)
      if (nranks > 1) then
         es = [e_os, e_ss]
         call allreduce(comm, es, op=MPI_SUM)
         call bcast(comm, es, 2, 0)
         e_os = es(1); e_ss = es(2)
      end if
      res%e_os = e_os
      res%e_ss = e_ss
      res%e_corr = e_os + e_ss
      if (.not. present(la)) call la_own%release()
   end subroutine trc_rimp2_run

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   !
   ! For occupied i: part(2 jb - 1) = sum_a (ia|jb)^2 / D, part(2 jb) =
   ! sum_a [(ia|jb)^2 - (ia|jb)(ib|ja)] / D -- one thread per (j, b), a summed
   ! inside, so no reduction crosses threads and the host adds nov numbers.
   ! Diagonal terms are in the sum: (ia|ib)^2 - (ia|ib)(ib|ia) = 0 for the
   ! same-spin part, as it must, and i = j counts once as in the closed formula.
   !
   subroutine pair_energy(i, no, nv, nov, eo, ev, ki, part)
      integer, intent(in) :: i, no, nv, nov
      real(dp), intent(in) :: eo(no), ev(nv), ki(nv, nov)
      real(dp), intent(out) :: part(2*nov)
      integer :: jb
      do concurrent(jb=1:nov)
         call pair_energy_body(i, jb, no, nv, nov, eo, ev, ki, part(2*jb - 1), part(2*jb))
      end do
   end subroutine pair_energy

   pure subroutine pair_energy_body(i, jb, no, nv, nov, eo, ev, ki, pos, pss)
      !$acc routine seq
      integer, intent(in) :: i, jb, no, nv, nov
      real(dp), intent(in) :: eo(no), ev(nv), ki(nv, nov)
      real(dp), intent(out) :: pos, pss
      integer :: j, bb, a
      real(dp) :: t, u, d, eij
      j = (jb - 1)/nv + 1
      bb = jb - (j - 1)*nv
      eij = eo(i) + eo(j)
      pos = 0.0_dp; pss = 0.0_dp
      do a = 1, nv
         t = ki(a, jb)                       ! (ia|jb)
         u = ki(bb, (j - 1)*nv + a)          ! (ib|ja)
         d = eij - ev(a) - ev(bb)
         pos = pos + t*t/d
         pss = pss + (t*t - t*u)/d
      end do
   end subroutine pair_energy_body

end module trc_rimp2_driver
