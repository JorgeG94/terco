!
! The density-fitting DERIVATIVE integrals (scripts/gen_df_grad.py's kernels)
! against pyscf's int3c2e_ip1, int3c2e_ip2 and int2c2e_ip1, before any of the
! fitted gradient's own contraction is involved -- "is d(mu nu|P)/dR right",
! the same kernel-check role check_df_vs_pyscf.F90 plays for the integrals
! themselves.
!
! TWO HALVES.
!
!   RAW. Every class kernel is called on the HOST (trc_df_grad_kernels'
!   df_grad_3c_block / df_grad_2c_block) once per Cartesian component with
!   a one-hot coefficient block, which returns that one derivative integral
!   -- the exact code the device drivers run. The full tensors go into the
!   dump; test/pyscf_refs/generate.py's gen_df_grad_ints compares them
!   element by element against pyscf at generation time and prints the
!   worst difference (the tensors are too big to commit). Checked here
!   too: the B-centre derivative of (mu nu|P) equals the A-centre one of
!   (nu mu|P), and translational invariance of every derivative block.
!
!   CONTRACTED. The device drivers (df_grad_3c_driver over the energy's
!   own triangular pair list, df_grad_2c_driver over auxiliary shell
!   pairs) contract a fixed synthetic symmetric Y(mu,nu,P) and Z(P,Q) into
!   a per-atom gradient. That is compared against (a) the same contraction
!   of the raw host tensors above -- the drivers' weights, pair list,
!   lanes and atomics -- and (b) pyscf's own integrals contracted with the
!   same Y/Z, committed as dfg_ints_<case>.ref.
!
! PER-AO NORMALISATION as check_df_vs_pyscf.F90: from the (i|j) metric of
! the orbital basis and the (P|Q) metric of the auxiliary one, both
! trc_df_2c, dumped for generate.py.
!
! SKIPS CLEANLY when a basis needs more angular momentum than the
! derivative kernels were generated for (DFG_LMAX/DFG_LMAX_AUX).
!
module check_df_grad_ints_m
   use pic_types, only: dp
   use trc_api, only: trc_basis_t
   use trc_df_grad_kernels, only: df_grad_3c_block, df_grad_2c_block
   implicit none
   private
   public :: raw_3c, raw_2c, synth_y, synth_z, ao_atoms

contains

   subroutine raw_3c(b, aux, ip1, gbsw)
      !! ip1(x, mu, nu, P) = d(mu nu|P)/dA_x, A the centre of mu -- the
      !! one-hot A-derivative of every ordered shell pair. `gbsw` is the
      !! B-derivative of the same calls, stored at (x, nu, mu, P) so it
      !! should equal ip1 exactly.
      type(trc_basis_t), intent(in) :: b, aux
      real(dp), intent(out) :: ip1(3, b%nao, b%nao, aux%nao), gbsw(3, b%nao, b%nao, aux%nao)
      integer :: ii, jj, kk, la, lb, lc, na, nb, nc, a, c, e, slot
      real(dp), allocatable :: y(:)
      real(dp) :: ga(3), gb(3)
      ! Sized for the largest class of THESE bases, not a fixed ceiling: a
      ! DF_LMAX=3 build reaches (ff|g), 10*10*15 components.
      la = maxval(b%sh_l); lc = maxval(aux%sh_l)
      allocate (y(((la + 1)*(la + 2)/2)**2*(lc + 1)*(lc + 2)/2))
      do kk = 1, aux%nshell
         lc = aux%sh_l(kk); nc = (lc + 1)*(lc + 2)/2
         do jj = 1, b%nshell
            lb = b%sh_l(jj); nb = (lb + 1)*(lb + 2)/2
            do ii = 1, b%nshell
               la = b%sh_l(ii); na = (la + 1)*(la + 2)/2
               do e = 0, nc - 1
                  do c = 0, nb - 1
                     do a = 0, na - 1
                        slot = 1 + a + na*(c + nb*e)
                        y = 0.0_dp
                        y(slot) = 1.0_dp
                        call df_grad_3c_block(la, lb, lc, b%sh_np(ii), b%sh_e(:, ii), b%sh_c(:, ii), b%sh_r(:, ii), &
                                              b%sh_np(jj), b%sh_e(:, jj), b%sh_c(:, jj), b%sh_r(:, jj), &
                                              aux%sh_np(kk), aux%sh_e(:, kk), aux%sh_c(:, kk), aux%sh_r(:, kk), &
                                              y, ga, gb)
                        ip1(:, b%sh_ao(ii) + a, b%sh_ao(jj) + c, aux%sh_ao(kk) + e) = ga
                        gbsw(:, b%sh_ao(jj) + c, b%sh_ao(ii) + a, aux%sh_ao(kk) + e) = gb
                     end do
                  end do
               end do
            end do
         end do
      end do
      deallocate (y)
   end subroutine raw_3c

   subroutine raw_2c(aux, ip1)
      !! ip1(x, P, Q) = d(P|Q)/dA_x, A the centre of P.
      type(trc_basis_t), intent(in) :: aux
      real(dp), intent(out) :: ip1(3, aux%nao, aux%nao)
      integer :: ii, jj, la, lc, na, nc, a, c
      real(dp), allocatable :: y(:)
      real(dp) :: ga(3)
      la = maxval(aux%sh_l)
      allocate (y(((la + 1)*(la + 2)/2)**2))
      do jj = 1, aux%nshell
         lc = aux%sh_l(jj); nc = (lc + 1)*(lc + 2)/2
         do ii = 1, aux%nshell
            la = aux%sh_l(ii); na = (la + 1)*(la + 2)/2
            do c = 0, nc - 1
               do a = 0, na - 1
                  y = 0.0_dp
                  y(1 + a + na*c) = 1.0_dp
                  call df_grad_2c_block(la, lc, aux%sh_np(ii), aux%sh_e(:, ii), aux%sh_c(:, ii), aux%sh_r(:, ii), &
                                        aux%sh_np(jj), aux%sh_e(:, jj), aux%sh_c(:, jj), aux%sh_r(:, jj), y, ga)
                  ip1(:, aux%sh_ao(ii) + a, aux%sh_ao(jj) + c) = ga
               end do
            end do
         end do
      end do
      deallocate (y)
   end subroutine raw_2c

   subroutine synth_y(n, naux, y)
      !! A fixed, dense coefficient tensor, symmetric in its first two
      !! indices (the drivers rely on that), O(1) entries of both signs.
      integer, intent(in) :: n, naux
      real(dp), intent(out) :: y(n, n, naux)
      integer :: i, j, p
      do p = 1, naux
         do j = 1, n
            do i = 1, n
               y(i, j, p) = 0.3_dp*cos(real(3*i + 5*j + 7*p, dp)) + 0.3_dp*cos(real(5*i + 3*j + 7*p, dp)) &
                            + 0.05_dp*sin(real(i + j + 2*p, dp))
            end do
         end do
      end do
   end subroutine synth_y

   subroutine synth_z(naux, z)
      integer, intent(in) :: naux
      real(dp), intent(out) :: z(naux, naux)
      integer :: i, j
      do j = 1, naux
         do i = 1, naux
            z(i, j) = 0.2_dp*cos(real(2*i + 3*j, dp)) + 0.2_dp*cos(real(3*i + 2*j, dp)) + 0.1_dp
         end do
      end do
   end subroutine synth_z

   subroutine ao_atoms(b, at)
      !! at(mu): the atom AO mu sits on, through its shell.
      type(trc_basis_t), intent(in) :: b
      integer, intent(out) :: at(b%nao)
      integer :: ish, l
      do ish = 1, b%nshell
         l = b%sh_l(ish)
         at(b%sh_ao(ish):b%sh_ao(ish) + (l + 1)*(l + 2)/2 - 1) = b%sh_at(ish)
      end do
   end subroutine ao_atoms

end module check_df_grad_ints_m

program check_df_grad_ints
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_2c
   use trc_df_grad_kernels, only: df_grad_3c_plan, df_grad_3c_driver, df_grad_2c_driver, DFG_LMAX, DFG_LMAX_AUX, &
                                  DFG_NPCLS
   use trc_basis_json, only: trc_basis_from_json
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_report
   use check_df_grad_ints_m, only: raw_3c, raw_2c, synth_y, synth_z, ao_atoms
   implicit none

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: REL_TOL = 1.0e-10_dp
   !! Relative to the largest gradient component (O(10-100) for these
   !! synthetic O(1) coefficients summed over every integral).
   character(len=256) :: xyzfile, basis_json, aux_json, dump, case_name
   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl
   type(error_t) :: err
   integer :: nat, nao, naux, unit, mu, nu, p, q, nfail
   integer, allocatable :: at_z(:), aoat(:), axat(:)
   real(dp), allocatable :: r(:, :), omat(:, :), jmat(:, :)
   real(dp), allocatable :: ip1(:, :, :, :), gbsw(:, :, :, :), ip2c(:, :, :)
   real(dp), allocatable :: y(:, :, :), z(:, :), g3(:, :), g2(:, :), g3h(:, :), g2h(:, :)
   real(dp) :: dsw, inv3, inv2, d3, d2, scale3, scale2, dmax(3)
   integer :: coff(0:DFG_NPCLS)
   integer, allocatable :: q_a(:), q_b(:), q_o(:), q_n(:)
   real(dp), allocatable :: pd(:, :), cq0(:), d0(:, :)

   if (command_argument_count() < 5) then
      print '(a)', 'usage: check_df_grad_ints <xyz> <basis.json> <aux.json> <dump.bin> <case>'; stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basis_json)
   call get_command_argument(3, aux_json)
   call get_command_argument(4, dump)
   call get_command_argument(5, case_name)

   call read_xyz(trim(xyzfile), nat, at_z, r)
   call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', 'check_df_grad_ints: '//err%get_message(); stop 1
   end if
   call trc_basis_from_json(trim(aux_json), nat, at_z, r, aux, err, force_cart=.true.)
   if (err%has_error()) then
      print '(a)', 'check_df_grad_ints: '//err%get_message(); stop 1
   end if
   if (maxval(bas%sh_l) > DFG_LMAX .or. maxval(aux%sh_l) > DFG_LMAX_AUX) then
      print '(a)', '  SKIP: a basis needs more angular momentum than this build''s DF derivative kernels.'
      stop 0
   end if
   nao = bas%nao; naux = aux%nao
   call bas%to_device()
   call aux%to_device()
   call pl%build(bas, THRESH)
   call pl%to_device()

   ! --- normalisation probes, as check_df_vs_pyscf.F90 ---
   allocate (omat(nao, nao), jmat(naux, naux))
   !$acc enter data create(omat, jmat)
   call trc_df_2c(bas, omat)
   call trc_df_2c(aux, jmat)
   !$acc update self(omat, jmat)
   !$acc exit data delete(omat, jmat)

   ! --- raw derivative tensors, host, one-hot ---
   allocate (ip1(3, nao, nao, naux), gbsw(3, nao, nao, naux), ip2c(3, naux, naux))
   call raw_3c(bas, aux, ip1, gbsw)
   call raw_2c(aux, ip2c)
   nfail = 0
   dsw = maxval(abs(ip1 - gbsw))
   call ref_report(nfail, dsw < 1.0e-13_dp*max(1.0_dp, maxval(abs(ip1))), &
                   'd/dB (mu nu|P) == d/dA (nu mu|P), max|diff|', dsw)
   ! Translational invariance per element: dA + dB + dC = 0 makes dC itself
   ! a derived quantity, so check the two-centre one, (P|Q)'s dA = -dC.
   inv2 = 0.0_dp
   do q = 1, naux
      do p = 1, naux
         inv2 = max(inv2, maxval(abs(ip2c(:, p, q) + ip2c(:, q, p))))
      end do
   end do
   call ref_report(nfail, inv2 < 1.0e-13_dp*max(1.0_dp, maxval(abs(ip2c))), &
                   'd/dA (P|Q) == -d/dA (Q|P), max|diff|', inv2)

   ! --- contracted, on the device drivers ---
   allocate (y(nao, nao, naux), z(naux, naux), g3(3, nat), g2(3, nat), g3h(3, nat), g2h(3, nat))
   call synth_y(nao, naux, y)
   call synth_z(naux, z)
   g3 = 0.0_dp; g2 = 0.0_dp
   ! Y alone (the drivers' K-style dense input): the J part's c and D zero.
   allocate (cq0(naux), d0(nao, nao))
   cq0 = 0.0_dp; d0 = 0.0_dp
   call df_grad_3c_plan(bas%nshell, bas%maxnp, bas%sh_l, bas%sh_np, bas%sh_e, bas%sh_c, bas%sh_r, &
                        pl%npair, pl%p_i, pl%p_j, coff, q_a, q_b, q_o, q_n, pd)
   !$acc enter data copyin(y, z, g3, g2, cq0, d0, q_a, q_b, q_o, q_n, pd)
   call df_grad_3c_driver(bas%nshell, nao, bas%sh_ao, bas%sh_at, bas%sh_r, coff, pl%npair, q_a, q_b, q_o, q_n, &
                          size(pd, 2), pd, aux%nshell, aux%maxnp, aux%sh_l, aux%sh_np, aux%sh_ao, aux%sh_at, &
                          aux%sh_e, aux%sh_c, aux%sh_r, nat, .true., nao, naux, y, 0.0_dp, naux, 0, cq0, d0, g3)
   call df_grad_2c_driver(aux%nshell, naux, aux%maxnp, aux%sh_l, aux%sh_np, aux%sh_ao, aux%sh_at, &
                          aux%sh_e, aux%sh_c, aux%sh_r, nat, z, g2)
   !$acc update self(g3, g2)
   !$acc exit data delete(y, z, g3, g2, cq0, d0, q_a, q_b, q_o, q_n, pd)

   ! --- the same contraction from the raw host tensors ---
   allocate (aoat(nao), axat(naux))
   call ao_atoms(bas, aoat)
   call ao_atoms(aux, axat)
   g3h = 0.0_dp; g2h = 0.0_dp
   do p = 1, naux
      do nu = 1, nao
         do mu = 1, nao
            ! A on mu's atom, B on nu's (via the symmetric partner), C on P's.
            g3h(:, aoat(mu)) = g3h(:, aoat(mu)) + y(mu, nu, p)*ip1(:, mu, nu, p)
            g3h(:, aoat(nu)) = g3h(:, aoat(nu)) + y(mu, nu, p)*ip1(:, nu, mu, p)
            g3h(:, axat(p)) = g3h(:, axat(p)) - y(mu, nu, p)*(ip1(:, mu, nu, p) + ip1(:, nu, mu, p))
         end do
      end do
   end do
   do q = 1, naux
      do p = 1, naux
         g2h(:, axat(p)) = g2h(:, axat(p)) + z(p, q)*ip2c(:, p, q)
         g2h(:, axat(q)) = g2h(:, axat(q)) - z(p, q)*ip2c(:, p, q)
      end do
   end do
   scale3 = max(1.0_dp, maxval(abs(g3h)))
   scale2 = max(1.0_dp, maxval(abs(g2h)))
   d3 = maxval(abs(g3 - g3h))
   d2 = maxval(abs(g2 - g2h))
   call ref_report(nfail, d3 < 1.0e-12_dp*scale3, 'device 3c driver == host raw contraction', d3)
   call ref_report(nfail, d2 < 1.0e-12_dp*scale2, 'device 2c driver == host raw contraction', d2)
   inv3 = maxval(abs(sum(g3, dim=2)))
   inv2 = maxval(abs(sum(g2, dim=2)))
   call ref_report(nfail, inv3 < 1.0e-12_dp*scale3, '3c translational invariance', inv3)
   call ref_report(nfail, inv2 < 1.0e-12_dp*scale2, '2c translational invariance', inv2)

   open (newunit=unit, file=trim(dump), access='stream', form='unformatted', status='replace')
   write (unit) int(nat, int64), int(nao, int64), int(naux, int64)
   write (unit) int(at_z, int64)
   write (unit) r
   write (unit) omat
   write (unit) jmat
   write (unit) ip1
   write (unit) ip2c
   write (unit) y
   write (unit) z
   write (unit) g3
   write (unit) g2
   close (unit)
   print '(a)', 'check_df_grad_ints: wrote '//trim(dump)

   ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
   block
      integer :: u, nat_ref, nao_ref, naux_ref
      real(dp), allocatable :: g3_ref(:, :), g2_ref(:, :)
      call ref_require('dfg_ints_'//trim(case_name)//'.ref', u)
      call ref_read_int(u, 'NAT', nat_ref)
      call ref_read_int(u, 'NAO', nao_ref)
      call ref_read_int(u, 'NAUX', naux_ref)
      if (nat_ref /= nat .or. nao_ref /= nao .or. naux_ref /= naux) then
         print '(a)', 'check_df_grad_ints: pyscf reference shape does not match this build'; stop 1
      end if
      allocate (g3_ref(3, nat), g2_ref(3, nat))
      call ref_read_matrix(u, 'G3', 3, nat, g3_ref)
      call ref_read_matrix(u, 'G2', 3, nat, g2_ref)
      call ref_close(u)
      dmax(1) = maxval(abs(g3 - g3_ref))
      dmax(2) = maxval(abs(g2 - g2_ref))
      call ref_report(nfail, dmax(1) < REL_TOL*scale3, '3c Y:d(mu nu|P)/dR vs pyscf ip1/ip2, worst', dmax(1))
      call ref_report(nfail, dmax(2) < REL_TOL*scale2, '2c Z:d(P|Q)/dR vs pyscf ip1, worst', dmax(2))
      print '(a,2es10.2)', '  (gradient scales, 3c and 2c: ', scale3, scale2
   end block

   call pl%release()
   call bas%release()
   call aux%release()
   if (nfail > 0) then
      print '(a)', 'check_df_grad_ints: RESULT: FAIL'; stop 1
   end if
   print '(a)', 'check_df_grad_ints: RESULT: PASS'

end program check_df_grad_ints
