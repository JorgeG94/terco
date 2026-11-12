!
! The amplitude update, DIIS and the iteration driver of RI-CCSD
! (src/post_hf/ccsd/trc_rccsd_update.F90, _diis.F90, _driver.F90) against the DF-RCCSD of pyscf, on a
! water and on a C1 hydrogen peroxide, cc-pVDZ with cc-pVDZ-RIFIT, spherical, frozen core 1 and 2 (the
! molecules, SCF and probes of check_rccsd_ints, shared through rccsd_test_support).
!
! References (test/pyscf_refs/rccsd_ints_*.ref, generate.py's gen_rccsd_ints):
!   P_t1_new, P_t2_new  pyscf's update_amps applied ONCE to the FIXED amplitudes of the other checks
!                       (t2 = MP2, t1_ia = (0.03 x + 0.02 y + 0.01 z)_ia / (e_i - e_a))
!   E_UPDATE1           the correlation energy of those updated amplitudes
!   E_MP2, E_CCSD_CORR  the MP2 and the converged DF-CCSD correlation energies (conv_tol 1e-12, normt 1e-10)
!   P_t1_conv, P_t2_conv, T1_DIAG   the converged amplitudes (as the sign-invariant probes) and get_t1_diagnostic
!
! Cases, per molecule:
!   1 one update from the fixed amplitudes against the pyscf probes and E_UPDATE1
!   2 the same update against an independent host implementation of update_amps written term by term from the
!     pyscf einsums (explicit (ac|bd) from b_vv; nothing of the modules under test), elementwise
!   3 a converged run (conv_energy 1e-12, conv_amplitude 1e-10): E_CCSD_CORR, E_MP2, the probes of the converged
!     amplitudes, the T1 diagnostic, the converged flag and the iteration count
!   4 the converged amplitudes are a fixed point: one more update changes them by less than 1e-9 (max abs)
!   5 DIIS on and off (diis_size 0) reach the same energy; DIIS takes no more iterations
!   6 ladder forced assembled and forced direct reach the same energy
!   7 non-convergence is reported (max_iterations 2: not converged, a message, a finite energy), and a NaN in
!     the starting amplitudes stops the run as not converged
!   8 refusals (an ECP, an open shell, two ranks, no virtuals, a frozen core as large as the occupied space,
!     an unsized orbital array): refused, a message, nothing left allocated, and a normal run afterwards
!   9 a loose energy tolerance does not stop the run while the amplitudes still move (the criterion is the AND)
! and, once, without a molecule:
!   - DIIS on a linear fixed-point problem x = A x + b reaches the solution in at most n + 1 steps; a restart on
!     parallel error vectors; a NaN step is measured, flagged and never appended; the step measure
!     (max and 2-norm, across more than one chunk) against a known vector
!
module rccsd_host_update
   !! update_amps of pyscf/cc/rccsd.py (non-CC2 branch, canonical RHF orbitals, level_shift 0) as plain loops,
   !! each einsum written out with the integral blocks and t1/t2 only: (ac|bd) is summed from b_vv here, the Fock
   !! matrix is diag(e) and is added to Foo, Fvv, Loo, Lvv and subtracted again exactly as pyscf does. Nothing of
   !! the modules under test is used except the arrays of the integrals type.
   use pic_types, only: dp
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   implicit none
   private
   public :: host_update

contains

   subroutine host_update(ints, t1, t2, t1n, t2n)
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), intent(out) :: t1n(:, :), t2n(:, :, :, :)

      real(dp), allocatable :: foo(:, :), fvv(:, :), fov(:, :), loo(:, :), lvv(:, :), tau(:, :, :, :)
      real(dp), allocatable :: woooo(:, :, :, :), wvoov(:, :, :, :), wvovo(:, :, :, :), wvvvv(:, :, :, :)
      integer :: no, nv

      no = ints%no
      nv = ints%nv
      allocate (foo(no, no), fvv(nv, nv), fov(no, nv), loo(no, no), lvv(nv, nv), tau(no, no, nv, nv))
      allocate (woooo(no, no, no, no), wvoov(nv, no, no, nv), wvovo(nv, no, nv, no), wvvvv(nv, nv, nv, nv))
      tau = t2
      call add_t1t1(no, nv, t1, tau)
      call host_f_and_l(ints, t1, t2, foo, fvv, fov, loo, lvv)
      call host_w(ints, t1, t2, woooo, wvoov, wvovo, wvvvv)
      call host_t1_equation(ints, t1, t2, foo, fvv, fov, t1n)
      call host_t2_equation(ints, t1, t2, tau, loo, lvv, woooo, wvoov, wvovo, wvvvv, t2n)
   end subroutine host_update

   subroutine add_t1t1(no, nv, t1, tau)
      !! tau += einsum('ia,jb->ijab', t1, t1)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv)
      real(dp), intent(inout) :: tau(no, no, nv, nv)
      integer :: i, j, a, b

      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  tau(i, j, a, b) = tau(i, j, a, b) + t1(i, a)*t1(j, b)
               end do
            end do
         end do
      end do
   end subroutine add_t1t1

   subroutine host_f_and_l(ints, t1, t2, foo, fvv, fov, loo, lvv)
      !! imd.cc_Foo, cc_Fvv, cc_Fov, Loo, Lvv, then `Foo[diag] -= mo_e_o`, `Fvv[diag] -= mo_e_v` and the same for L
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), intent(out) :: foo(:, :), fvv(:, :), fov(:, :), loo(:, :), lvv(:, :)
      real(dp), allocatable :: ovvv(:, :, :, :)

      integer :: no, nv, i, k, l, a, c, d
      real(dp) :: s

      no = ints%no
      nv = ints%nv
      allocate (ovvv(no, nv, nv, nv))
      ! (ia|bc) from the factors: trc_rccsd_ints no longer holds the dense block
      call ints%fill_ovvv(ovvv)
      ! cc_Foo: 2 'kcld,ilcd->ki' - 'kdlc,ilcd->ki' + 2 'kcld,ic,ld->ki' - 'kdlc,ic,ld->ki' + foo
      do i = 1, no
         do k = 1, no
            s = 0.0_dp
            do d = 1, nv
               do c = 1, nv
                  do l = 1, no
                     s = s + 2.0_dp*ints%ovov(k, c, l, d)*t2(i, l, c, d) - ints%ovov(k, d, l, c)*t2(i, l, c, d) &
                         + 2.0_dp*ints%ovov(k, c, l, d)*t1(i, c)*t1(l, d) - ints%ovov(k, d, l, c)*t1(i, c)*t1(l, d)
                  end do
               end do
            end do
            foo(k, i) = s
            if (k == i) foo(k, i) = foo(k, i) + ints%eps_o(i)
         end do
      end do
      ! cc_Fvv: -2 'kcld,klad->ac' + 'kdlc,klad->ac' - 2 'kcld,ka,ld->ac' + 'kdlc,ka,ld->ac' + fvv
      do c = 1, nv
         do a = 1, nv
            s = 0.0_dp
            do d = 1, nv
               do l = 1, no
                  do k = 1, no
                     s = s - 2.0_dp*ints%ovov(k, c, l, d)*t2(k, l, a, d) + ints%ovov(k, d, l, c)*t2(k, l, a, d) &
                         - 2.0_dp*ints%ovov(k, c, l, d)*t1(k, a)*t1(l, d) + ints%ovov(k, d, l, c)*t1(k, a)*t1(l, d)
                  end do
               end do
            end do
            fvv(a, c) = s
            if (a == c) fvv(a, c) = fvv(a, c) + ints%eps_v(a)
         end do
      end do
      ! cc_Fov: 2 'kcld,ld->kc' - 'kdlc,ld->kc' + fov (zero)
      do c = 1, nv
         do k = 1, no
            s = 0.0_dp
            do d = 1, nv
               do l = 1, no
                  s = s + 2.0_dp*ints%ovov(k, c, l, d)*t1(l, d) - ints%ovov(k, d, l, c)*t1(l, d)
               end do
            end do
            fov(k, c) = s
         end do
      end do
      ! Loo = cc_Foo + 2 'lcki,lc->ki' - 'kcli,lc->ki' with ovoo(l,c,k,i) = (lc|ki) = ooov(k,i,l,c)
      do i = 1, no
         do k = 1, no
            s = foo(k, i)
            do c = 1, nv
               do l = 1, no
                  s = s + 2.0_dp*ints%ooov(k, i, l, c)*t1(l, c) - ints%ooov(l, i, k, c)*t1(l, c)
               end do
            end do
            loo(k, i) = s
         end do
      end do
      ! Lvv = cc_Fvv + 2 'kdac,kd->ac' - 'kcad,kd->ac'
      do c = 1, nv
         do a = 1, nv
            s = fvv(a, c)
            do d = 1, nv
               do k = 1, no
                  s = s + 2.0_dp*ovvv(k, d, a, c)*t1(k, d) - ovvv(k, c, a, d)*t1(k, d)
               end do
            end do
            lvv(a, c) = s
         end do
      end do
      do i = 1, no
         foo(i, i) = foo(i, i) - ints%eps_o(i)
         loo(i, i) = loo(i, i) - ints%eps_o(i)
      end do
      do a = 1, nv
         fvv(a, a) = fvv(a, a) - ints%eps_v(a)
         lvv(a, a) = lvv(a, a) - ints%eps_v(a)
      end do
   end subroutine host_f_and_l

   subroutine host_w(ints, t1, t2, woooo, wvoov, wvovo, wvvvv)
      !! imd.cc_Woooo, cc_Wvoov, cc_Wvovo, cc_Wvvvv with the explicit (ac|bd) = sum_P b_vv(ac,P) b_vv(bd,P)
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), intent(out) :: woooo(:, :, :, :), wvoov(:, :, :, :), wvovo(:, :, :, :), wvvvv(:, :, :, :)
      real(dp), allocatable :: ovvv(:, :, :, :)

      integer :: no, nv, i, j, k, l, a, b, c, d, p
      real(dp) :: s

      no = ints%no
      nv = ints%nv
      allocate (ovvv(no, nv, nv, nv))
      ! (ia|bc) from the factors: trc_rccsd_ints no longer holds the dense block
      call ints%fill_ovvv(ovvv)
      ! cc_Woooo: 'lcki,jc->klij' + 'kclj,ic->klij' + 'kcld,ijcd->klij' + 'kcld,ic,jd->klij' + oooo.transpose(0,2,1,3)
      do j = 1, no
         do i = 1, no
            do l = 1, no
               do k = 1, no
                  s = ints%oooo(k, i, l, j)
                  do c = 1, nv
                     s = s + ints%ooov(k, i, l, c)*t1(j, c) + ints%ooov(l, j, k, c)*t1(i, c)
                     do d = 1, nv
                        s = s + ints%ovov(k, c, l, d)*(t2(i, j, c, d) + t1(i, c)*t1(j, d))
                     end do
                  end do
                  woooo(k, l, i, j) = s
               end do
            end do
         end do
      end do
      ! cc_Wvoov
      do c = 1, nv
         do i = 1, no
            do k = 1, no
               do a = 1, nv
                  ! ovvo.transpose(2,0,3,1): ovvo(k,c,a,i) = (kc|ai) = ovov(k,c,i,a)
                  s = ints%ovov(k, c, i, a)
                  do d = 1, nv
                     s = s + ovvv(k, c, a, d)*t1(i, d)            ! 'kcad,id->akic'
                  end do
                  do l = 1, no
                     s = s - ints%ooov(l, i, k, c)*t1(l, a)            ! 'kcli,la->akic', ovoo(k,c,l,i) = ooov(l,i,k,c)
                     do d = 1, nv
                        s = s - 0.5_dp*ints%ovov(l, d, k, c)*t2(i, l, d, a) - 0.5_dp*ints%ovov(l, c, k, d)*t2(i, l, a, d) &
                            - ints%ovov(l, d, k, c)*t1(i, d)*t1(l, a) + ints%ovov(l, d, k, c)*t2(i, l, a, d)
                     end do
                  end do
                  wvoov(a, k, i, c) = s
               end do
            end do
         end do
      end do
      ! cc_Wvovo
      do i = 1, no
         do c = 1, nv
            do k = 1, no
               do a = 1, nv
                  s = ints%oovv(k, i, a, c)                              ! oovv.transpose(2,0,3,1)
                  do d = 1, nv
                     s = s + ovvv(k, d, a, c)*t1(i, d)            ! 'kdac,id->akci'
                  end do
                  do l = 1, no
                     s = s - ints%ooov(k, i, l, c)*t1(l, a)            ! 'lcki,la->akci', ovoo(l,c,k,i) = ooov(k,i,l,c)
                     do d = 1, nv
                        s = s - 0.5_dp*ints%ovov(l, c, k, d)*t2(i, l, d, a) - ints%ovov(l, c, k, d)*t1(i, d)*t1(l, a)
                     end do
                  end do
                  wvovo(a, k, c, i) = s
               end do
            end do
         end do
      end do
      ! cc_Wvvvv: 'kdac,kb->abcd' (-t1) + 'kcbd,ka->abcd' (-t1) + vvvv.transpose(0,2,1,3), vvvv(a,c,b,d) = (ac|bd)
      do d = 1, nv
         do c = 1, nv
            do b = 1, nv
               do a = 1, nv
                  s = 0.0_dp
                  do p = 1, ints%naux
                     s = s + ints%b_vv(a + nv*(c - 1), p)*ints%b_vv(b + nv*(d - 1), p)
                  end do
                  do k = 1, no
                     s = s - ovvv(k, d, a, c)*t1(k, b) - ovvv(k, c, b, d)*t1(k, a)
                  end do
                  wvvvv(a, b, c, d) = s
               end do
            end do
         end do
      end do
   end subroutine host_w

   subroutine host_t1_equation(ints, t1, t2, foo, fvv, fov, t1n)
      !! The t1new lines of update_amps, in pyscf's order (fov = 0, so the '-2 kc,ka,ic' line is absent)
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :), foo(:, :), fvv(:, :), fov(:, :)
      real(dp), intent(out) :: t1n(:, :)
      real(dp), allocatable :: ovvv(:, :, :, :)

      integer :: no, nv, i, k, l, a, c, d
      real(dp) :: s

      no = ints%no
      nv = ints%nv
      allocate (ovvv(no, nv, nv, nv))
      ! (ia|bc) from the factors: trc_rccsd_ints no longer holds the dense block
      call ints%fill_ovvv(ovvv)
      do a = 1, nv
         do i = 1, no
            s = 0.0_dp
            do c = 1, nv
               s = s + fvv(a, c)*t1(i, c)                                         ! 'ac,ic->ia'
            end do
            do k = 1, no
               s = s - foo(k, i)*t1(k, a)                                         ! '-ki,ka->ia'
            end do
            do c = 1, nv
               do k = 1, no
                  s = s + 2.0_dp*fov(k, c)*t2(k, i, c, a) - fov(k, c)*t2(i, k, c, a) + fov(k, c)*t1(i, c)*t1(k, a)
                  ! '2 kc,kica->ia', '-kc,ikca->ia', 'kc,ic,ka->ia'
                  s = s + 2.0_dp*ints%ovov(k, c, i, a)*t1(k, c) - ints%oovv(k, i, a, c)*t1(k, c)
                  ! '2 kcai,kc->ia' with ovvo(k,c,a,i) = ovov(k,c,i,a); '-kiac,kc->ia'
               end do
            end do
            do d = 1, nv
               do c = 1, nv
                  do k = 1, no
                     s = s + 2.0_dp*ovvv(k, d, a, c)*t2(i, k, c, d) - ovvv(k, c, a, d)*t2(i, k, c, d) &
                         + 2.0_dp*ovvv(k, d, a, c)*t1(k, d)*t1(i, c) - ovvv(k, c, a, d)*t1(k, d)*t1(i, c)
                     ! '2 kdac,ikcd->ia', '-kcad,ikcd->ia', '2 kdac,kd,ic->ia', '-kcad,kd,ic->ia'
                  end do
               end do
            end do
            do c = 1, nv
               do l = 1, no
                  do k = 1, no
                     s = s - 2.0_dp*ints%ooov(k, i, l, c)*t2(k, l, a, c) + ints%ooov(l, i, k, c)*t2(k, l, a, c) &
                         - 2.0_dp*ints%ooov(k, i, l, c)*t1(l, c)*t1(k, a) + ints%ooov(l, i, k, c)*t1(l, c)*t1(k, a)
                     ! '-2 lcki,klac->ia', 'kcli,klac->ia', '-2 lcki,lc,ka->ia', 'kcli,lc,ka->ia'
                  end do
               end do
            end do
            t1n(i, a) = s/(ints%eps_o(i) - ints%eps_v(a))
         end do
      end do
   end subroutine host_t1_equation

   subroutine host_t2_equation(ints, t1, t2, tau, loo, lvv, woooo, wvoov, wvovo, wvvvv, t2n)
      !! The t2new lines of update_amps, in pyscf's order
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :), tau(:, :, :, :), loo(:, :), lvv(:, :)
      real(dp), intent(in) :: woooo(:, :, :, :), wvoov(:, :, :, :), wvovo(:, :, :, :), wvvvv(:, :, :, :)
      real(dp), intent(out) :: t2n(:, :, :, :)
      real(dp), allocatable :: ovvv(:, :, :, :)

      real(dp), allocatable :: tmp(:, :, :, :), tmp2(:, :, :, :)
      integer :: no, nv, i, j, k, l, a, b, c, d
      real(dp) :: s

      no = ints%no
      nv = ints%nv
      allocate (ovvv(no, nv, nv, nv))
      ! (ia|bc) from the factors: trc_rccsd_ints no longer holds the dense block
      call ints%fill_ovvv(ovvv)
      allocate (tmp(no, no, nv, nv))
      ! tmp2 = 'kibc,ka->abic' (-t1) + ovvv.transpose(1,3,0,2), i.e. tmp2(a,b,i,c) = ovvv(i,a,c,b) - sum_k oovv(k,i,b,c) t1(k,a)
      ! tmp = 'abic,jc->ijab'; t2new = tmp + tmp.transpose(1,0,3,2)
      allocate (tmp2(nv, nv, no, nv))
      do c = 1, nv
         do i = 1, no
            do b = 1, nv
               do a = 1, nv
                  s = ovvv(i, a, c, b)
                  do k = 1, no
                     s = s - ints%oovv(k, i, b, c)*t1(k, a)
                  end do
                  tmp2(a, b, i, c) = s
               end do
            end do
         end do
      end do
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  s = 0.0_dp
                  do c = 1, nv
                     s = s + tmp2(a, b, i, c)*t1(j, c)
                  end do
                  tmp(i, j, a, b) = s
               end do
            end do
         end do
      end do
      call sym_add(no, nv, 1.0_dp, tmp, t2n, .true.)
      deallocate (tmp2)
      ! tmp2 = 'kcai,jc->akij' + ovoo.transpose(1,3,0,2): tmp2(a,k,i,j) = sum_c ovvo(k,c,a,i) t1(j,c) + (ia|jk)
      ! tmp = 'akij,kb->ijab'; t2new -= tmp + tmp.transpose(1,0,3,2)
      allocate (tmp2(nv, no, no, no))
      do j = 1, no
         do i = 1, no
            do k = 1, no
               do a = 1, nv
                  s = ints%ooov(j, k, i, a)
                  do c = 1, nv
                     s = s + ints%ovov(k, c, i, a)*t1(j, c)
                  end do
                  tmp2(a, k, i, j) = s
               end do
            end do
         end do
      end do
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  s = 0.0_dp
                  do k = 1, no
                     s = s + tmp2(a, k, i, j)*t1(k, b)
                  end do
                  tmp(i, j, a, b) = s
               end do
            end do
         end do
      end do
      call sym_add(no, nv, -1.0_dp, tmp, t2n, .false.)
      ! t2new += eris.ovov.transpose(0,2,1,3)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  t2n(i, j, a, b) = t2n(i, j, a, b) + ints%ovov(i, a, j, b)
               end do
            end do
         end do
      end do
      ! t2new += 'klij,klab->ijab' (Woooo, tau) and 'abcd,ijcd->ijab' (Wvvvv, tau)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  s = 0.0_dp
                  do l = 1, no
                     do k = 1, no
                        s = s + woooo(k, l, i, j)*tau(k, l, a, b)
                     end do
                  end do
                  do d = 1, nv
                     do c = 1, nv
                        s = s + wvvvv(a, b, c, d)*tau(i, j, c, d)
                     end do
                  end do
                  t2n(i, j, a, b) = t2n(i, j, a, b) + s
               end do
            end do
         end do
      end do
      ! tmp = 'ac,ijcb->ijab' (Lvv, t2); t2new += tmp + tmp.T
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  s = 0.0_dp
                  do c = 1, nv
                     s = s + lvv(a, c)*t2(i, j, c, b)
                  end do
                  tmp(i, j, a, b) = s
               end do
            end do
         end do
      end do
      call sym_add(no, nv, 1.0_dp, tmp, t2n, .false.)
      ! tmp = 'ki,kjab->ijab' (Loo, t2); t2new -= tmp + tmp.T
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  s = 0.0_dp
                  do k = 1, no
                     s = s + loo(k, i)*t2(k, j, a, b)
                  end do
                  tmp(i, j, a, b) = s
               end do
            end do
         end do
      end do
      call sym_add(no, nv, -1.0_dp, tmp, t2n, .false.)
      ! tmp = 2 'akic,kjcb->ijab' - 'akci,kjcb->ijab'; t2new += tmp + tmp.T
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  s = 0.0_dp
                  do c = 1, nv
                     do k = 1, no
                        s = s + 2.0_dp*wvoov(a, k, i, c)*t2(k, j, c, b) - wvovo(a, k, c, i)*t2(k, j, c, b)
                     end do
                  end do
                  tmp(i, j, a, b) = s
               end do
            end do
         end do
      end do
      call sym_add(no, nv, 1.0_dp, tmp, t2n, .false.)
      ! tmp = 'akic,kjbc->ijab'; t2new -= tmp + tmp.T
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  s = 0.0_dp
                  do c = 1, nv
                     do k = 1, no
                        s = s + wvoov(a, k, i, c)*t2(k, j, b, c)
                     end do
                  end do
                  tmp(i, j, a, b) = s
               end do
            end do
         end do
      end do
      call sym_add(no, nv, -1.0_dp, tmp, t2n, .false.)
      ! tmp = 'bkci,kjac->ijab'; t2new -= tmp + tmp.T
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  s = 0.0_dp
                  do c = 1, nv
                     do k = 1, no
                        s = s + wvovo(b, k, c, i)*t2(k, j, a, c)
                     end do
                  end do
                  tmp(i, j, a, b) = s
               end do
            end do
         end do
      end do
      call sym_add(no, nv, -1.0_dp, tmp, t2n, .false.)
      ! t2new /= eijab
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  t2n(i, j, a, b) = t2n(i, j, a, b)/(ints%eps_o(i) + ints%eps_o(j) - ints%eps_v(a) - ints%eps_v(b))
               end do
            end do
         end do
      end do
   end subroutine host_t2_equation

   subroutine sym_add(no, nv, factor, tmp, t2n, first)
      !! t2n (+)= factor (tmp + tmp.transpose(1,0,3,2)); `first` starts t2n instead of adding
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: factor, tmp(no, no, nv, nv)
      real(dp), intent(inout) :: t2n(no, no, nv, nv)
      logical, intent(in) :: first
      integer :: i, j, a, b

      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  if (first) t2n(i, j, a, b) = 0.0_dp
                  t2n(i, j, a, b) = t2n(i, j, a, b) + factor*(tmp(i, j, a, b) + tmp(j, i, b, a))
               end do
            end do
         end do
      end do
   end subroutine sym_add

end module rccsd_host_update

module check_rccsd_run_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, rccsd_inter_energy
   use trc_rccsd_ladder, only: RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_DIRECT, RCCSD_LADDER_HALF, trc_rccsd_ladder_options_t
   use trc_rccsd_update, only: trc_rccsd_update_t
   use trc_rccsd_steps, only: RCCSD_ST_INTEGRALS, RCCSD_ST_METRIC, RCCSD_ST_CHOLESKY, RCCSD_ST_3C, RCCSD_ST_HALF, &
                              RCCSD_ST_TRSM, RCCSD_ST_BLOCKS, RCCSD_ST_INTER, RCCSD_ST_W, RCCSD_ST_LADDER, &
                              RCCSD_ST_DRESSING, RCCSD_ST_T1, RCCSD_ST_T2, RCCSD_ST_RING, RCCSD_ST_DIIS, RCCSD_ST_ENERGY, &
                              RCCSD_NSTEP, RCCSD_STEP_NAME
   use trc_step_table, only: step_table_t
   use trc_rccsd_diis, only: trc_rccsd_diis_t
   use trc_rccsd_driver, only: trc_rccsd_options_t, trc_rccsd_result_t, trc_rccsd_state_t, trc_rccsd_run, trc_rccsd_solve, &
                               RCCSD_TIME_INTEGRALS, RCCSD_TIME_INTER, RCCSD_TIME_W, RCCSD_TIME_LADDER, RCCSD_TIME_AMPS, &
                               RCCSD_TIME_DIIS, RCCSD_TIME_TOTAL
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: NMOL, MOL_LABEL, MOL_REF, MOL_NFROZEN, TOL_PROBE, PROBE_FLOOR, load_molecule, &
                                 run_scf, mo_dipoles, probe2, probe4
   use rccsd_host_update, only: host_update
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector
   implicit none
   private
   public :: collect_check_rccsd_run

   integer, parameter :: NCASE = 11
   integer, parameter :: REF_E_CCSD = 1, REF_E_UPDATE1 = 2, REF_T1_DIAG = 3, REF_E_MP2 = 4  !! indices into the scalars read from a .ref
   character(len=4), parameter :: SKIP_TAG(5) = [character(len=4) :: "oooo", "ooov", "oovv", "ovov", "ovvv"]
   character(len=4), parameter :: SKIP_ONE(7) = [character(len=4) :: "t1", "t2", "Foo", "Fvv", "Fov", "Loo", "Lvv"]
   character(len=6), parameter :: SKIP_W(3) = [character(len=6) :: "Woooo", "Wvoov", "Wvovo"]
   character(len=10), parameter :: SKIP_LADDER(3) = [character(len=10) :: "ladder", "ladder_w", "dressing"]
   real(dp), parameter :: TOL_ENERGY = 1.0e-10_dp     !! converged E_CCSD_CORR against pyscf, absolute
   real(dp), parameter :: TOL_CONV_PROBE = 1.0e-9_dp  !! probes of the converged amplitudes against pyscf's, relative
   real(dp), parameter :: TOL_HOST = 1.0e-11_dp       !! one update against the host loops, max abs
   real(dp), parameter :: TOL_FIXED = 1.0e-9_dp       !! one more update from the converged amplitudes, max abs
   real(dp), parameter :: TOL_DIIS_ENERGY = 1.0e-10_dp
   real(dp), parameter :: TOL_LADDER_ENERGY = 1.0e-11_dp

   logical :: computed = .false.
   logical :: ok_case(NCASE, NMOL)
   character(len=900) :: msg_case(NCASE, NMOL)

contains

   subroutine collect_check_rccsd_run(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water: one update vs pyscf", w1), &
                  new_unittest("water: one update vs host einsum loops", w2), &
                  new_unittest("water: converged run vs pyscf", w3), &
                  new_unittest("water: converged amplitudes are a fixed point", w4), &
                  new_unittest("water: DIIS on and off", w5), &
                  new_unittest("water: ladder forms", w6), &
                  new_unittest("water: non-convergence and NaN", w7), &
                  new_unittest("water: refusals", w8), &
                  new_unittest("water: loose energy, tight amplitudes", w9), &
                  new_unittest("water: step table counts the flops of every GEMM", w10), &
                  new_unittest("water: step table times add up", w11), &
                  new_unittest("h2o2: one update vs pyscf", h1), &
                  new_unittest("h2o2: one update vs host einsum loops", h2), &
                  new_unittest("h2o2: converged run vs pyscf", h3), &
                  new_unittest("h2o2: converged amplitudes are a fixed point", h4), &
                  new_unittest("h2o2: DIIS on and off", h5), &
                  new_unittest("h2o2: ladder forms", h6), &
                  new_unittest("h2o2: non-convergence and NaN", h7), &
                  new_unittest("h2o2: refusals", h8), &
                  new_unittest("h2o2: loose energy, tight amplitudes", h9), &
                  new_unittest("h2o2: step table counts the flops of every GEMM", h10), &
                  new_unittest("h2o2: step table times add up", h11), &
                  new_unittest("diis: linear fixed-point problem in at most n + 1 steps", diis_linear), &
                  new_unittest("diis: packed and full histories agree on a symmetric problem", diis_packed), &
                  new_unittest("diis: host and device histories agree, with a ring that wraps", diis_host), &
                  new_unittest("diis: parallel errors restart the history", diis_singular), &
                  new_unittest("diis: a NaN step is flagged and never appended", diis_nan), &
                  new_unittest("diis: step measure against a known vector", diis_step) &
                  ]
   end subroutine collect_check_rccsd_run

   subroutine w1(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1, 1)
   end subroutine w1

   subroutine w2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 2, 1)
   end subroutine w2

   subroutine w3(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 3, 1)
   end subroutine w3

   subroutine w4(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 4, 1)
   end subroutine w4

   subroutine w5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5, 1)
   end subroutine w5

   subroutine w6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 6, 1)
   end subroutine w6

   subroutine w7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 7, 1)
   end subroutine w7

   subroutine w8(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 8, 1)
   end subroutine w8

   subroutine w9(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 9, 1)
   end subroutine w9

   subroutine w10(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 10, 1)
   end subroutine w10

   subroutine w11(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 11, 1)
   end subroutine w11

   subroutine h10(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 10, 2)
   end subroutine h10

   subroutine h11(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 11, 2)
   end subroutine h11

   subroutine h1(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1, 2)
   end subroutine h1

   subroutine h2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 2, 2)
   end subroutine h2

   subroutine h3(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 3, 2)
   end subroutine h3

   subroutine h4(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 4, 2)
   end subroutine h4

   subroutine h5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5, 2)
   end subroutine h5

   subroutine h6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 6, 2)
   end subroutine h6

   subroutine h7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 7, 2)
   end subroutine h7

   subroutine h8(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 8, 2)
   end subroutine h8

   subroutine h9(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 9, 2)
   end subroutine h9

   subroutine run_case(error, ic, im)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ic, im
      call ensure_computed()
      call check(error, ok_case(ic, im), trim(msg_case(ic, im)))
   end subroutine run_case

   subroutine ensure_computed()
      integer :: im, dev
      if (computed) return
      dev = trc_bind_device(0)
      ok_case = .false.
      msg_case = "not computed"
      do im = 1, NMOL
         call compute_molecule(im)
      end do
      computed = .true.
   end subroutine ensure_computed

   subroutine record(im, ic, passed, text)
      !! One (ok, message) per case, also printed so a run shows the measured numbers
      integer, intent(in) :: im, ic
      logical, intent(in) :: passed
      character(len=*), intent(in) :: text
      ok_case(ic, im) = passed
      msg_case(ic, im) = text
      print "(a,l1,a,a)", "  check_rccsd_run ", passed, " ", text
   end subroutine record

   function rel_probe(p, ref, floor) result(w)
      !! Worst |p - ref| / max(|ref|, floor)
      real(dp), intent(in) :: p(:), ref(:), floor
      real(dp) :: w
      w = maxval(abs(p - ref)/max(abs(ref), floor))
   end function rel_probe

   function options_for(im, nocc) result(opts)
      !! The tight options of the converged runs
      integer, intent(in) :: im, nocc
      type(trc_rccsd_options_t) :: opts
      opts%nocc = nocc
      opts%nfrozen = MOL_NFROZEN(im)
      opts%conv_energy = 1.0e-12_dp
      opts%conv_amplitude = 1.0e-10_dp
      opts%max_iterations = 200
   end function options_for

   subroutine skip_probes(ref_unit)
      !! The reader is sequential: the probes of the other checks, in the order of the file, up to P_t1_new
      integer, intent(in) :: ref_unit

      integer :: i
      real(dp) :: ref4(4), ref2(3)

      do i = 1, size(SKIP_TAG)
         call ref_read_vector(ref_unit, "P_"//SKIP_TAG(i), 4, ref4)
      end do
      do i = 1, size(SKIP_ONE)
         if (SKIP_ONE(i) == "t2") then
            call ref_read_vector(ref_unit, "P_"//trim(SKIP_ONE(i)), 4, ref4)
         else
            call ref_read_vector(ref_unit, "P_"//trim(SKIP_ONE(i)), 3, ref2)
         end if
      end do
      do i = 1, size(SKIP_W)
         call ref_read_vector(ref_unit, "P_"//trim(SKIP_W(i)), 4, ref4)
      end do
      do i = 1, size(SKIP_LADDER)
         call ref_read_vector(ref_unit, "P_"//trim(SKIP_LADDER(i)), 4, ref4)
      end do
   end subroutine skip_probes

   subroutine read_refs(ref_unit, p_t1_new, p_t2_new, p_t1_conv, p_t2_conv, refs)
      !! The numbers of the cases, in the order of the file; `refs` is indexed by the REF_* parameters
      integer, intent(in) :: ref_unit
      real(dp), intent(out) :: p_t1_new(3), p_t2_new(4), p_t1_conv(3), p_t2_conv(4), refs(:)

      real(dp) :: skipped
      integer :: nskip

      call ref_read_int(ref_unit, "NO", nskip)
      call ref_read_int(ref_unit, "NV", nskip)
      call ref_read_int(ref_unit, "NAUX", nskip)
      call ref_read_int(ref_unit, "NFROZEN", nskip)
      call ref_read_real(ref_unit, "E_SCF", skipped)
      call ref_read_real(ref_unit, "E_MP2", refs(REF_E_MP2))
      call ref_read_real(ref_unit, "E_FIXED", skipped)
      call skip_probes(ref_unit)
      call ref_read_vector(ref_unit, "P_t1_new", 3, p_t1_new)
      call ref_read_vector(ref_unit, "P_t2_new", 4, p_t2_new)
      call ref_read_vector(ref_unit, "P_t1_conv", 3, p_t1_conv)
      call ref_read_vector(ref_unit, "P_t2_conv", 4, p_t2_conv)
      call ref_read_real(ref_unit, "E_CCSD_CORR", refs(REF_E_CCSD))
      call ref_read_real(ref_unit, "E_T", skipped)
      call ref_read_real(ref_unit, "E_CCSD_T", skipped)
      call ref_read_real(ref_unit, "E_UPDATE1", refs(REF_E_UPDATE1))
      call ref_read_real(ref_unit, "T1_DIAG", refs(REF_T1_DIAG))
   end subroutine read_refs

   subroutine compute_molecule(im)
      !! SCF, the integrals, and every number the cases compare
      integer, intent(in) :: im

      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_linalg_t) :: la
      integer, allocatable :: zint(:)
      real(dp), allocatable :: dip(:, :, :)
      integer :: nocc, nfz, nmo, no, nv, ref_unit
      logical :: supported
      character(len=:), allocatable :: label
      real(dp) :: refs(REF_E_MP2), p_t1_new(3), p_t2_new(4), p_t1_conv(3), p_t2_conv(4)
      type(trc_rccsd_state_t) :: state

      label = trim(MOL_LABEL(im))//": "
      nfz = MOL_NFROZEN(im)
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) then
         print "(a)", "SKIP check_rccsd_run: this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
         stop 0
      end if
      nmo = bas%nao_sph
      call run_scf(bas, nocc, res)
      call la%init(2)

      ! ---- the converged run first: it also fixes (no, nv) and the dipole blocks ------------
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), options_for(im, nocc), la, state)
      if (state%result%refused) error stop "check_rccsd_run: the converged run was refused: "//trim(state%result%message)
      no = state%result%no
      nv = state%result%nv
      allocate (dip(no + nv, no + nv, 4))
      call mo_dipoles(bas, pl, res%cmo(:, nfz + 1:nmo, 1), dip)

      call ref_require(trim(MOL_REF(im)), ref_unit)
      call read_refs(ref_unit, p_t1_new, p_t2_new, p_t1_conv, p_t2_conv, refs)
      call ref_close(ref_unit)

      call update_cases(im, label, bas, aux, pl, res, la, nocc, dip, p_t1_new, p_t2_new, refs(REF_E_UPDATE1))
      call converged_cases(im, label, state, dip, p_t1_conv, p_t2_conv, refs(REF_E_CCSD), refs(REF_E_MP2), refs(REF_T1_DIAG))
      call fixed_point_case(im, label, la, state)
      call variant_cases(im, label, bas, aux, pl, res, la, nocc, state%result)
      call failure_cases(im, label, bas, aux, pl, res, la, nocc, no, nv)
      call step_cases(im, label, bas, aux, pl, res, la, nocc)
      call state%release()
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine compute_molecule

   subroutine update_cases(im, label, bas, aux, pl, res, la, nocc, dip, p_t1_new, p_t2_new, e_update1)
      !! Cases 1 and 2: one update from the fixed amplitudes, against pyscf and against the host loops
      integer, intent(in) :: im, nocc
      character(len=*), intent(in) :: label
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res
      type(trc_linalg_t), intent(inout) :: la
      real(dp), intent(in) :: dip(:, :, :), p_t1_new(:), p_t2_new(:), e_update1

      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_ints_options_t) :: bopts
      type(trc_rccsd_amps_t) :: amps, newamps
      type(trc_rccsd_update_t) :: upd
      real(dp), allocatable :: t1fix(:, :), t2fix(:, :, :, :), t1h(:, :), t2h(:, :, :, :)
      real(dp) :: d1, d2, de, es, ed, host_t1, host_t2, scale
      integer :: no, nv, a, i
      logical :: ok

      bopts%nalpha = nocc
      bopts%nbeta = nocc
      bopts%nfrozen = MOL_NFROZEN(im)
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      if (ints%status%has_error()) error stop "check_rccsd_run: build: "//ints%status%get_message()
      no = ints%no
      nv = ints%nv
      allocate (t1fix(no, nv), t2fix(no, no, nv, nv))
      do a = 1, nv
         do i = 1, no
            t1fix(i, a) = (0.03_dp*dip(i, no + a, 1) + 0.02_dp*dip(i, no + a, 2) + 0.01_dp*dip(i, no + a, 3))/ &
                          (ints%eps_o(i) - ints%eps_v(a))
         end do
      end do
      call mp2_t2(ints, t2fix)
      call amps%init(no, nv)
      call newamps%init(no, nv)
      amps%t1 = t1fix
      amps%t2 = t2fix
      call amps%push()
      call upd%init(ints)
      if (upd%status%has_error()) error stop "check_rccsd_run: update init: "//upd%status%get_message()
      call upd%step(ints, amps, newamps, la)
      if (upd%status%has_error()) error stop "check_rccsd_run: update: "//upd%status%get_message()
      call newamps%pull()
      call rccsd_inter_energy(newamps, upd%scr, es, ed)

      d1 = rel_probe(probe2(newamps%t1, "ov", dip, no), p_t1_new, PROBE_FLOOR)
      d2 = rel_probe(probe4(newamps%t2, "oovv", dip, no), p_t2_new, PROBE_FLOOR)
      de = abs(es + ed - e_update1)
      call record(im, 1, max(d1, d2) < TOL_PROBE .and. de < TOL_ENERGY .and. maxval(abs(newamps%t2)) > 1.0e-3_dp, &
                  label//"one update vs pyscf, worst relative (floor 0.1): t1 "//to_string(d1, "es9.2")//", t2 "// &
                  to_string(d2, "es9.2")//", energy after the update, absolute "//to_string(de, "es9.2"))

      allocate (t1h(no, nv), t2h(no, no, nv, nv))
      call host_update(ints, t1fix, t2fix, t1h, t2h)
      host_t1 = maxval(abs(newamps%t1 - t1h))
      host_t2 = maxval(abs(newamps%t2 - t2h))
      scale = max(1.0_dp, maxval(abs(t2h)))
      ok = host_t1 < TOL_HOST*scale .and. host_t2 < TOL_HOST*scale .and. maxval(abs(t1h)) > 1.0e-3_dp
      call record(im, 2, ok, label//"one update vs the host einsum loops, max abs: t1 "//to_string(host_t1, "es9.2")// &
                  ", t2 "//to_string(host_t2, "es9.2")//" (|t1|max "//to_string(maxval(abs(t1h)), "es9.2")//", |t2|max "// &
                  to_string(maxval(abs(t2h)), "es9.2")//")")
      call upd%release()
      call amps%release()
      call newamps%release()
      call ints%release()
   end subroutine update_cases

   subroutine mp2_t2(ints, t2)
      !! t2(i,j,a,b) = (ia|jb) / (e_i + e_j - e_a - e_b)
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(out) :: t2(:, :, :, :)
      integer :: i, j, a, b

      do b = 1, ints%nv
         do a = 1, ints%nv
            do j = 1, ints%no
               do i = 1, ints%no
                  t2(i, j, a, b) = ints%ovov(i, a, j, b)/(ints%eps_o(i) + ints%eps_o(j) - ints%eps_v(a) - ints%eps_v(b))
               end do
            end do
         end do
      end do
   end subroutine mp2_t2

   subroutine converged_cases(im, label, state, dip, p_t1_conv, p_t2_conv, e_ref, e_mp2_ref, t1_diag_ref)
      !! Case 3: the converged run against pyscf
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_state_t), intent(in) :: state
      real(dp), intent(in) :: dip(:, :, :), p_t1_conv(:), p_t2_conv(:), e_ref, e_mp2_ref, t1_diag_ref

      real(dp) :: d_e, d_mp2, d1, d2, d_diag, d_split
      integer :: no

      no = state%result%no
      d_e = abs(state%result%e_ccsd - e_ref)
      d_mp2 = abs(state%result%e_mp2 - e_mp2_ref)
      d1 = rel_probe(probe2(state%amps%t1, "ov", dip, no), p_t1_conv, PROBE_FLOOR)
      d2 = rel_probe(probe4(state%amps%t2, "oovv", dip, no), p_t2_conv, PROBE_FLOOR)
      d_diag = abs(state%result%t1_diagnostic - t1_diag_ref)
      d_split = abs(state%result%e_singles + state%result%e_doubles - state%result%e_ccsd)
      call record(im, 3, state%result%converged .and. d_e < TOL_ENERGY .and. d_mp2 < TOL_ENERGY .and. &
                  max(d1, d2) < TOL_CONV_PROBE .and. d_diag < 1.0e-9_dp .and. d_split < 1.0e-14_dp, &
                  label//"converged in "//to_string(state%result%iterations)//" iterations, converged = "// &
                  merge("T", "F", state%result%converged)//", E_CCSD(corr) "//to_string(state%result%e_ccsd, "f20.14")// &
                  " minus pyscf "//to_string(d_e, "es9.2")//", E_MP2 minus pyscf "//to_string(d_mp2, "es9.2")// &
                  ", worst probe (relative, floor 0.1) t1 "//to_string(d1, "es9.2")//" t2 "//to_string(d2, "es9.2")// &
                  ", T1 diagnostic "//to_string(state%result%t1_diagnostic, "f9.6")//" minus pyscf "// &
                  to_string(d_diag, "es9.2")//", last step max "//to_string(state%result%step_max, "es9.2")//" 2-norm "// &
                  to_string(state%result%step_norm2, "es9.2"))
   end subroutine converged_cases

   subroutine fixed_point_case(im, label, la, state)
      !! Case 4: one more update from the converged amplitudes changes them by less than TOL_FIXED
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_state_t), intent(inout) :: state

      type(trc_rccsd_amps_t) :: newamps
      type(trc_rccsd_update_t) :: upd
      real(dp) :: d1, d2, es, ed, d_ret

      call upd%init(state%ints)
      if (upd%status%has_error()) error stop "check_rccsd_run: fixed point: init: "//upd%status%get_message()
      call rccsd_inter_energy(state%amps, upd%scr, es, ed)
      d_ret = abs(es + ed - state%result%e_ccsd)
      call newamps%init(state%ints%no, state%ints%nv)
      call upd%step(state%ints, state%amps, newamps, la)
      if (upd%status%has_error()) error stop "check_rccsd_run: fixed point: step: "//upd%status%get_message()
      call newamps%pull()
      d1 = maxval(abs(newamps%t1 - state%amps%t1))
      d2 = maxval(abs(newamps%t2 - state%amps%t2))
      ! the (T) step reads state%amps: the energy reported must be the functional of exactly those amplitudes
      call record(im, 4, max(d1, d2) < TOL_FIXED .and. maxval(abs(state%amps%t2)) > 1.0e-3_dp .and. d_ret < 1.0e-13_dp, &
                  label//"one more update from the converged amplitudes, max abs change: t1 "//to_string(d1, "es9.2")// &
                  ", t2 "//to_string(d2, "es9.2")//" (|t1|max "//to_string(maxval(abs(state%amps%t1)), "es9.2")//"); energy functional of the returned amplitudes "// &
                  "minus result%e_ccsd "//to_string(d_ret, "es9.2"))
      call upd%release()
      call newamps%release()
   end subroutine fixed_point_case

   subroutine variant_cases(im, label, bas, aux, pl, res, la, nocc, base)
      !! Cases 5 and 6: DIIS off, the ladder forms
      integer, intent(in) :: im, nocc
      character(len=*), intent(in) :: label
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_result_t), intent(in) :: base

      type(trc_rccsd_options_t) :: opts
      type(trc_rccsd_result_t) :: plain, assembled, direct

      opts = options_for(im, nocc)
      opts%diis_size = 0
      call trc_rccsd_run(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), plain, opts, la)
      call record(im, 5, plain%converged .and. base%converged .and. abs(plain%e_ccsd - base%e_ccsd) < TOL_DIIS_ENERGY .and. &
                  base%iterations <= plain%iterations, &
                  label//"DIIS 8: "//to_string(base%iterations)//" iterations, no DIIS: "//to_string(plain%iterations)// &
                  " iterations, energies differ by "//to_string(abs(plain%e_ccsd - base%e_ccsd), "es9.2")// &
                  ", restarts "//to_string(base%diis_restarts))

      opts = options_for(im, nocc)
      opts%ladder%form = RCCSD_LADDER_ASSEMBLED
      call trc_rccsd_run(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), assembled, opts, la)
      opts%ladder%form = RCCSD_LADDER_DIRECT
      call trc_rccsd_run(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), direct, opts, la)
      call record(im, 6, assembled%converged .and. direct%converged .and. &
                  abs(assembled%e_ccsd - direct%e_ccsd) < TOL_LADDER_ENERGY .and. &
                  abs(assembled%e_ccsd - base%e_ccsd) < TOL_LADDER_ENERGY, &
                  label//"ladder assembled vs direct, converged energies differ by "// &
                  to_string(abs(assembled%e_ccsd - direct%e_ccsd), "es9.2")//" ("//to_string(assembled%iterations)// &
                  " and "//to_string(direct%iterations)//" iterations)")
   end subroutine variant_cases

   subroutine failure_cases(im, label, bas, aux, pl, res, la, nocc, no, nv)
      !! Cases 7, 8 and 9: non-convergence and NaN, refusals, the AND of the criteria
      integer, intent(in) :: im, nocc, no, nv
      character(len=*), intent(in) :: label
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res
      type(trc_linalg_t), intent(inout) :: la

      type(trc_rccsd_options_t) :: opts
      type(trc_rccsd_result_t) :: short, loose, refused
      type(trc_rccsd_state_t) :: state, nan_state, early, skew
      real(dp) :: d_early
      real(dp) :: e_skew(2), gap_skew(2)
      logical :: ok_skew
      integer :: i, j, a, b
      logical :: ok_early
      type(trc_rccsd_result_t) :: normal
      logical :: ok_short, ok_nan, ok_refuse
      integer :: nmo, k
      character(len=900) :: refuse_text
      real(dp) :: bad_eps(1)

      opts = options_for(im, nocc)
      opts%max_iterations = 2
      call trc_rccsd_run(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), short, opts, la)
      ! stopped early with DIIS on, where the extrapolated and the plain amplitudes still differ visibly
      opts%max_iterations = 3
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, early)
      d_early = returned_gap(early, la)
      ok_early = d_early < 1.0e-13_dp .and. early%result%iterations == 3 .and. .not. early%result%converged
      call early%release()
      ok_short = .not. short%converged .and. short%iterations == 2 .and. len_trim(short%message) > 0 .and. &
                 abs(short%e_ccsd) < 10.0_dp .and. short%e_ccsd == short%e_ccsd .and. short%finite .and. .not. short%refused

      ! a NaN in the starting amplitudes: not converged, not finite, a message
      opts = options_for(im, nocc)
      call nan_state%amps%init(no, nv)
      nan_state%amps%t2 = 0.01_dp
      nan_state%amps%t2(1, 1, 1, 1) = ieee_nan()
      nan_state%use_start = .true.
      call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, nan_state)
      ok_nan = .not. nan_state%result%converged .and. .not. nan_state%result%finite .and. &
               index(nan_state%result%message, "NaN") > 0 .and. nan_state%result%iterations == 0 .and. &
               nan_state%result%e_ccsd == 0.0_dp .and. index(nan_state%result%message, "no finite energy") > 0
      ! a start WITHOUT the pair symmetry t2(j,i,b,a) = t2(i,j,a,b): two iterations from it must not depend on the
      ! ladder form. The half form's fast path assumes the symmetry (it would return the ladder of the symmetrised
      ! tau); the driver has to notice and take the general path. Two iterations, not one: the energy weights are
      ! pair-symmetric, so the first energy cannot see a symmetrised ladder; the second update can.
      do k = 1, 2
         opts = options_for(im, nocc)
         opts%max_iterations = 2
         opts%ladder%form = merge(RCCSD_LADDER_HALF, RCCSD_LADDER_ASSEMBLED, k == 1)
         call skew%amps%init(no, nv)
         skew%amps%t1 = 0.0_dp
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     skew%amps%t2(i, j, a, b) = 1.0e-2_dp*sin(real(i + 2*j + 3*a + 5*b, dp))
                  end do
               end do
            end do
         end do
         skew%use_start = .true.
         call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, skew)
         e_skew(k) = skew%result%e_ccsd
         ! The second iterate is a DIIS extrapolation. A history that kept half of t2 (valid only for symmetric
         ! amplitudes) would hand back an exactly pair-symmetric t2; the full history keeps the start's asymmetry.
         gap_skew(k) = pair_gap(skew%amps%t2)
         call skew%release()
      end do
      ok_skew = abs(e_skew(1) - e_skew(2)) < 1.0e-12_dp .and. e_skew(1) == e_skew(1) .and. minval(gap_skew) > 1.0e-6_dp
      call record(im, 7, ok_short .and. ok_nan .and. ok_early .and. ok_skew, label//"non-symmetric start, two "// &
                  "iterations, half minus assembled ladder "//to_string(e_skew(1) - e_skew(2), "es9.2")// &
                  ", max |t2(i,j,a,b) - t2(j,i,b,a)| after them "//to_string(minval(gap_skew), "es9.2")// &
                  " (E "//to_string(e_skew(2), "f16.10")//"); max_iterations 2: converged = "// &
                  merge("T", "F", short%converged)//", iterations "//to_string(short%iterations)//", E "// &
                  to_string(short%e_ccsd, "f16.10")//", message '"//trim(short%message)//"'; NaN start: converged = "// &
                  merge("T", "F", nan_state%result%converged)//", finite = "//merge("T", "F", nan_state%result%finite)// &
                  ", message '"//trim(nan_state%result%message)//"'; after 3 iterations the energy functional of the "// &
                  "returned amplitudes minus result%e_ccsd "//to_string(d_early, "es9.2"))
      call nan_state%release()

      ! refusals
      nmo = bas%nao_sph
      ok_refuse = .true.
      refuse_text = ""
      do k = 1, 7
         opts = options_for(im, nocc)
         state%use_start = .false.
         select case (k)
         case (1); opts%has_ecp = .true.
         case (2); opts%nbeta = nocc - 1
         case (3); opts%max_iterations = 0
         case (4); opts%nfrozen = nocc
         case (5); opts%nocc = nmo
         case (6); opts%nocc = 0
         case (7); opts%ladder%form = 7
         end select
         call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
         ok_refuse = ok_refuse .and. state%result%refused .and. .not. state%result%converged .and. &
                     len_trim(state%result%message) > 0 .and. .not. allocated(state%amps%t1) .and. &
                     .not. allocated(state%ints%ovov) .and. state%ints%stamp == 0
         refuse_text = trim(refuse_text)//" ["//trim(state%result%message)//"]"
      end do
      ! an orbital array of the wrong size is refused by the integrals
      opts = options_for(im, nocc)
      bad_eps = 0.0_dp
      call trc_rccsd_run(bas, aux, pl, res%cmo(:, :, 1), bad_eps, refused, opts, la)
      ok_refuse = ok_refuse .and. refused%refused .and. len_trim(refused%message) > 0
      refuse_text = trim(refuse_text)//" ["//trim(refused%message)//"]"
      ! and the same state, the same device handle, runs normally afterwards
      opts = options_for(im, nocc)
      opts%conv_energy = 1.0e-8_dp
      opts%conv_amplitude = 1.0e-6_dp
      call trc_rccsd_run(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), normal, opts, la)
      ok_refuse = ok_refuse .and. normal%converged .and. .not. normal%refused
      call record(im, 8, ok_refuse, label//"refusals (ECP, open shell, max_iterations 0, frozen = occupied, no virtuals, "// &
                  "nocc 0, invalid ladder form, wrong eps size), then a default run converged = "//merge("T", "F", normal%converged)//": "// &
                  trim(refuse_text))
      call state%release()

      ! the criteria are an AND: a loose energy tolerance must not stop a run whose amplitudes still move
      opts = options_for(im, nocc)
      opts%conv_energy = 1.0_dp
      call trc_rccsd_run(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), loose, opts, la)
      call record(im, 9, loose%converged .and. loose%step_max < opts%conv_amplitude .and. loose%iterations > 3, &
                  label//"conv_energy 1, conv_amplitude 1e-10: converged = "//merge("T", "F", loose%converged)// &
                  " after "//to_string(loose%iterations)//" iterations, last step max "// &
                  to_string(loose%step_max, "es9.2")//", E "//to_string(loose%e_ccsd, "f20.14"))
   end subroutine failure_cases

   pure function pass_flops(k, form, o, v, naux, nao, d_batch) result(f)
      !! Flops of ONE pass of step `k`, written here from the contractions themselves and not from the library's
      !! counters: 2 m n k of each GEMM, as polynomials in o, v, the auxiliary count and the Cartesian AO count.
      !! A step that has no GEMM (the kernels, the integral steps, DIIS, the energy) is 0. The half form's ladder
      !! also depends on its batch plan, `d_batch` (`half_ladder_flops`).
      integer, intent(in) :: k, form, o, v, naux, nao, d_batch
      integer(int64) :: f

      integer(int64) :: oo, vv, px, ao

      oo = o
      vv = v
      px = naux
      ao = nao
      select case (k)
      case (RCCSD_ST_CHOLESKY)
         f = px**3/3
      case (RCCSD_ST_HALF)
         ! C^T (mu nu|P) over all active orbitals, then the oo, ov and vv MO transforms of each (p nu|P)
         f = 2*(oo + vv)*ao*ao*px + 2*ao*px*(oo*oo + oo*vv + vv*vv)
      case (RCCSD_ST_TRSM)
         f = (oo*oo + oo*vv + vv*vv)*px*px
      case (RCCSD_ST_BLOCKS)
         ! oooo, ooov, oovv, ovov: (ia|bc) is not built
         f = 2*px*(oo**4 + oo**3*vv + 2*oo**2*vv**2)
      case (RCCSD_ST_INTER)
         ! Foo, Fvv, Fov, Loo; then Lvv from the factors: N1 = t1 b_vv (o, v naux, v), u = b_ov^T t1, 2 b_vv u,
         ! and the exchange N1 b_ov over (k,P) (v, v, o naux)
         f = 2*oo**3*vv**2 + 2*oo**2*vv**3 + 2*oo**2*vv**2 + 2*oo**3*vv + 4*oo*vv**2*px + 2*oo*vv*px + 2*vv**2*px
      case (RCCSD_ST_W)
         ! Woooo, three rings, two ooov t1 batches; the (ia|bc) t1 terms from the factors: N1a (v, o, v) x naux and
         ! b_ov N1a^T (o v, v o, naux) for Wvoov, O1 (o, o, v) x naux and O1 b_vv^T (o^2, v^2, naux) for Wvovo
         f = 2*oo**4*vv**2 + 6*oo**3*vv**3 + 4*oo**3*vv**2 + 2*oo*vv**2*px + 4*oo**2*vv**2*px + 2*oo**2*vv*px
      case (RCCSD_ST_LADDER)
         if (form == RCCSD_LADDER_ASSEMBLED) then
            f = 2*(px + oo**2)*vv**4
         else if (form == RCCSD_LADDER_HALF) then
            f = half_ladder_flops(o, v, naux, d_batch)
         else
            f = 4*px*oo**2*vv**3
         end if
      case (RCCSD_ST_DRESSING)
         if (form == RCCSD_LADDER_HALF) then
            ! the dressed build: bt = b_vv - t1^T b_ov (v, v naux, o), Y = tau ovov^T (o^2, o^2, v^2),
            ! G = Y t1 (o^3, v, o), X_b = G_b t1 (o^2, v, o) for each of v values of b
            f = 2*vv*vv*px*oo + 2*oo**4*vv**2 + 2*oo**4*vv + 2*oo**3*vv**2
         else
            ! Z and cmat per virtual, both halves, and the pack of each virtual batch from the factors,
            ! b_ov b_vv^T (o v, v nvb, naux), both halves: 4 o v^3 naux over the batches
            f = 4*oo**3*vv**3 + 4*oo**3*vv**2 + 4*oo*vv**3*px
         end if
      case (RCCSD_ST_T1)
         ! Y = ycomb b_ov (o v, naux, o v), Y b_vv^T (o, v, v naux), the lam term (o, v, o^2 v)
         f = 2*oo**2*vv**2*px + 2*oo*vv**2*px + 2*oo**3*vv**2
      case (RCCSD_ST_T2)
         ! term A's (ia|bc) t1 is b_ov N1^T (o v, o v, naux); the rest as before
         f = 2*oo**2*vv**3 + 2*oo**2*vv**2*px + 10*oo**3*vv**2 + 2*oo**4*vv**2
      case (RCCSD_ST_RING)
         f = 6*oo**3*vv**3
      case default
         f = 0
      end select
   end function pass_flops

   pure function half_ladder_flops(o, v, naux, d_batch) result(f)
      !! The half form's passes on a symmetric tau, written from its contractions: the first pass covers d_batch
      !! values of d; a pass from d0 has nc = v - d0 + 1 rows c and takes the most d, nd <= nc, with nc nd <= v d_batch
      !! (the C buffer) and at most the first pass's number of (c >= d) columns (the W2 buffer). Each pass is
      !! C = (v nc, v nd, naux) and one contraction (o^2, v^2, ncol), ncol = nd nc - nd (nd - 1)/2.
      integer, intent(in) :: o, v, naux, d_batch
      integer(int64) :: f

      integer(int64) :: d0, nc, nd, cap

      f = 0
      cap = int(d_batch, int64)*v - (int(d_batch, int64)*(d_batch - 1))/2
      d0 = 1
      do while (d0 <= v)
         nc = v - d0 + 1
         nd = min(nc, (int(v, int64)*d_batch)/nc)
         do while (nd > 1 .and. nd*nc - (nd*(nd - 1))/2 > cap)
            nd = nd - 1
         end do
         nd = max(nd, 1_int64)
         f = f + 2*(v*nc)*(v*nd)*naux + 2*(int(o, int64)*o)*(int(v, int64)*v)*(nd*nc - (nd*(nd - 1))/2)
         d0 = d0 + nd
      end do
   end function half_ladder_flops

   subroutine step_cases(im, label, bas, aux, pl, res, la, nocc)
      !! Cases 10 and 11: the step table of a run counts the flops of its GEMMs, for both ladder forms and a split
      !! auxiliary basis, and its times add up
      integer, intent(in) :: im, nocc
      character(len=*), intent(in) :: label
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res
      type(trc_linalg_t), intent(inout) :: la

      type(trc_rccsd_state_t) :: state
      type(trc_rccsd_options_t) :: opts
      integer(int64) :: want, got
      integer :: form, k, n, no, nv, naux, nbad
      logical :: ok_times
      character(len=900) :: bad
      real(dp) :: rows, pieces, tot

      nbad = 0
      bad = ""
      ok_times = .true.
      tot = 0.0_dp
      do form = RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_HALF
         opts = options_for(im, nocc)
         opts%conv_energy = 1.0e-8_dp
         opts%conv_amplitude = 1.0e-6_dp
         opts%ladder%form = form
         opts%ladder%d_batch = 3
         ! a split auxiliary basis in one of the runs: the integral counts must not depend on the blocks
         opts%aux_block = merge(0, 30, form == RCCSD_LADDER_ASSEMBLED)
         call trc_rccsd_solve(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), opts, la, state)
         if (state%result%refused) error stop "check_rccsd_run: step table: refused: "//trim(state%result%message)
         n = state%result%iterations
         no = state%result%no
         nv = state%result%nv
         naux = state%result%naux
         do k = 1, RCCSD_NSTEP
            want = pass_flops(k, form, no, nv, naux, bas%nao, 3)
            if (k >= RCCSD_ST_INTER) want = want*n
            got = state%steps%fl(k)%get()
            if (got /= want) then
               nbad = nbad + 1
               if (nbad <= 4) bad = trim(bad)//" ["//trim(adjustl(RCCSD_STEP_NAME(k)))//", form "//to_string(form)//": counted "// &
                                    to_string(got)//" expected "//to_string(want)//"]"
            end if
         end do
         if (form == RCCSD_LADDER_ASSEMBLED) then
            ! the times of the first run
            ok_times = all(state%steps%t >= 0.0_dp)
            rows = state%steps%t(RCCSD_ST_INTEGRALS) + sum(state%steps%t(RCCSD_ST_INTER:RCCSD_ST_ENERGY))
            pieces = sum(state%steps%t(RCCSD_ST_METRIC:RCCSD_ST_BLOCKS))
            tot = state%result%seconds(RCCSD_TIME_TOTAL)
            ok_times = ok_times .and. rows <= tot*(1.0_dp + 1.0e-9_dp) .and. pieces <= state%steps%t(RCCSD_ST_INTEGRALS)* &
                       (1.0_dp + 1.0e-9_dp) .and. state%steps%t(RCCSD_ST_LADDER) > 0.0_dp .and. &
                       state%steps%t(RCCSD_ST_RING) > 0.0_dp
            ok_times = ok_times .and. state%result%seconds(RCCSD_TIME_INTEGRALS) == state%steps%t(RCCSD_ST_INTEGRALS) .and. &
                       state%result%seconds(RCCSD_TIME_INTER) == state%steps%t(RCCSD_ST_INTER) .and. &
                       state%result%seconds(RCCSD_TIME_W) == state%steps%t(RCCSD_ST_W) .and. &
                       state%result%seconds(RCCSD_TIME_LADDER) == state%steps%t(RCCSD_ST_LADDER) + &
                       state%steps%t(RCCSD_ST_DRESSING) .and. &
                       state%result%seconds(RCCSD_TIME_AMPS) == state%steps%t(RCCSD_ST_T1) + state%steps%t(RCCSD_ST_T2) + &
                       state%steps%t(RCCSD_ST_RING) .and. state%result%seconds(RCCSD_TIME_DIIS) == state%steps%t(RCCSD_ST_DIIS)
         end if
         call state%release()
      end do
      nbad = nbad + step_pass_cases(im, bas, aux, pl, res, la, nocc, bad)
      call record(im, 10, nbad == 0, label//"step table flops, one run per ladder form (the direct one with a split auxiliary "// &
                  "basis) and one update per batch plan: "//to_string(nbad)//" counts differ from the independent "// &
                  "formulas"//trim(bad))
      call record(im, 11, ok_times, label//"step table times: rows non-negative, stage rows within the total "// &
                  to_string(tot, "f9.4")//" s, build pieces within the integrals row, result%seconds equal the table's totals")
   end subroutine step_cases

   function step_pass_cases(im, bas, aux, pl, res, la, nocc, bad) result(nbad)
      !! One update per ladder form and batch plan, with the tiny batches that make the batched GEMM calls many: what
      !! the ladder tallied from its own calls (`flops_pp`, `flops_dressing`) and what each step of the table holds
      !! are the formulas of `pass_flops`
      integer, intent(in) :: im, nocc
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res
      type(trc_linalg_t), intent(inout) :: la
      character(len=*), intent(inout) :: bad
      integer :: nbad

      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_ints_options_t) :: bopts
      type(trc_rccsd_amps_t) :: amps, newamps
      type(trc_rccsd_update_t) :: upd
      type(trc_rccsd_ladder_options_t) :: lopts
      type(step_table_t) :: steps
      integer(int64) :: want, got
      integer :: form, plan, k, no, nv, naux
      character(len=16) :: tag

      nbad = 0
      bopts%nalpha = nocc
      bopts%nbeta = nocc
      bopts%nfrozen = MOL_NFROZEN(im)
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      if (ints%status%has_error()) error stop "check_rccsd_run: step pass: build: "//ints%status%get_message()
      no = ints%no
      nv = ints%nv
      naux = ints%naux
      call amps%init(no, nv)
      call newamps%init(no, nv)
      call mp2_t2(ints, amps%t2)
      amps%t1 = 0.0_dp
      amps%t1(1, 1) = 0.01_dp
      call amps%push()
      do form = RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_HALF
         do plan = 1, 2
            lopts = trc_rccsd_ladder_options_t()
            lopts%form = form
            if (plan == 2) then
               lopts%cd_batch = 7
               lopts%d_batch = 2
               lopts%p_batch = 5
               lopts%virtual_batch = 3
            end if
            call upd%init(ints, lopts)
            if (upd%status%has_error()) error stop "check_rccsd_run: step pass: init: "//upd%status%get_message()
            call steps%init("rccsd", RCCSD_STEP_NAME)
            call upd%step(ints, amps, newamps, la, steps)
            if (upd%status%has_error()) error stop "check_rccsd_run: step pass: step: "//upd%status%get_message()
            tag = " ("//to_string(form)//","//to_string(plan)//")"
            ! the half form's plan 1 takes d_batch from the budget: the plan is an input of the formula
            want = pass_flops(RCCSD_ST_LADDER, form, no, nv, naux, bas%nao, upd%ls%d_batch)
            if (upd%ls%flops_pp /= want) call flag("ladder GEMMs made"//trim(tag), upd%ls%flops_pp, want)
            want = pass_flops(RCCSD_ST_DRESSING, form, no, nv, naux, bas%nao, upd%ls%d_batch)
            if (upd%ls%flops_dressing /= want) call flag("dressing GEMMs made"//trim(tag), upd%ls%flops_dressing, want)
            do k = RCCSD_ST_INTER, RCCSD_ST_RING
               want = pass_flops(k, form, no, nv, naux, bas%nao, upd%ls%d_batch)
               got = steps%fl(k)%get()
               if (got /= want) call flag(trim(adjustl(RCCSD_STEP_NAME(k)))//trim(tag), got, want)
            end do
            call upd%release()
         end do
      end do
      call amps%release()
      call newamps%release()
      call ints%release()

   contains

      subroutine flag(what, got, want)
         character(len=*), intent(in) :: what
         integer(int64), intent(in) :: got, want
         nbad = nbad + 1
         if (nbad <= 4) bad = trim(bad)//" ["//what//": counted "//to_string(got)//" expected "//to_string(want)//"]"
      end subroutine flag

   end function step_pass_cases

   function returned_gap(state, la) result(gap)
      !! |energy functional of state%amps - state%result%e_ccsd|: the reported energy belongs to the returned amplitudes
      type(trc_rccsd_state_t), intent(in) :: state
      type(trc_linalg_t), intent(inout) :: la
      real(dp) :: gap

      type(trc_rccsd_update_t) :: upd
      real(dp) :: es, ed

      call upd%init(state%ints)
      if (upd%status%has_error()) error stop "check_rccsd_run: returned gap: init: "//upd%status%get_message()
      call rccsd_inter_energy(state%amps, upd%scr, es, ed)
      gap = abs(es + ed - state%result%e_ccsd)
      call upd%release()
   end function returned_gap

   function pair_gap(t2) result(gap)
      !! max |t2(i,j,a,b) - t2(j,i,b,a)|
      real(dp), intent(in) :: t2(:, :, :, :)
      real(dp) :: gap
      integer :: i, j, a, b

      gap = 0.0_dp
      do b = 1, size(t2, 4)
         do a = 1, size(t2, 3)
            do j = 1, size(t2, 2)
               do i = 1, size(t2, 1)
                  gap = max(gap, abs(t2(i, j, a, b) - t2(j, i, b, a)))
               end do
            end do
         end do
      end do
   end function pair_gap

   function ieee_nan() result(x)
      !! A quiet NaN without a compile-time constant folding it away
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      real(dp) :: x
      x = ieee_value(1.0_dp, ieee_quiet_nan)
   end function ieee_nan

   ! ---------------------------------------------------------------------
   ! DIIS on its own, no molecule
   ! ---------------------------------------------------------------------

   function lcg(seed) result(x)
      !! A deterministic number in (-1, 1)
      integer(int64), intent(inout) :: seed
      real(dp) :: x
      seed = mod(seed*1103515245_int64 + 12345_int64, 2147483648_int64)
      x = 2.0_dp*real(seed, dp)/2147483648.0_dp - 1.0_dp
   end function lcg

   subroutine unpack_vector(amps, x)
      !! (t1, t2) of the host copy as one vector
      type(trc_rccsd_amps_t), intent(in) :: amps
      real(dp), intent(out) :: x(:)
      integer :: n1
      n1 = size(amps%t1)
      x(1:n1) = reshape(amps%t1, [n1])
      x(n1 + 1:) = reshape(amps%t2, [size(amps%t2)])
   end subroutine unpack_vector

   subroutine pack_vector(x, amps)
      !! One vector into the host and device copies of (t1, t2)
      real(dp), intent(in) :: x(:)
      type(trc_rccsd_amps_t), intent(inout) :: amps
      integer :: n1
      n1 = size(amps%t1)
      amps%t1 = reshape(x(1:n1), shape(amps%t1))
      amps%t2 = reshape(x(n1 + 1:), shape(amps%t2))
      call amps%push()
   end subroutine pack_vector

   subroutine diis_linear(error)
      !! x = A x + b with a random contraction A: the solution from plain iteration (400 steps), then DIIS from x = 0
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NO1 = 1, NV1 = 2, N = NO1*NV1 + NO1*NO1*NV1*NV1
      type(trc_linalg_t) :: la
      type(trc_rccsd_diis_t) :: diis
      type(trc_rccsd_amps_t) :: cur, nxt
      real(dp) :: a(N, N), b(N), x(N), xs(N), xn(N), worst
      integer :: i, j, k, reached
      integer(int64) :: seed
      logical :: ex

      seed = 12345_int64
      do j = 1, N
         do i = 1, N
            a(i, j) = 0.7_dp*lcg(seed)/real(N, dp)
         end do
         b(j) = lcg(seed)
      end do
      xs = 0.0_dp
      do k = 1, 400
         xs = matmul(a, xs) + b
      end do
      call la%init(2)
      call diis%init(NO1, NV1, N + 2)
      call cur%init(NO1, NV1)
      call nxt%init(NO1, NV1)
      x = 0.0_dp
      reached = -1
      do k = 1, N + 4
         xn = matmul(a, x) + b
         call pack_vector(x, cur)
         call pack_vector(xn, nxt)
         call diis%push(la, cur, nxt, .true.)
         call diis%next(la, nxt, cur, ex)
         call cur%pull()
         call unpack_vector(cur, x)
         worst = maxval(abs(x - xs))
         if (reached < 0 .and. worst < 1.0e-9_dp) reached = k
      end do
      call check(error, reached > 0 .and. reached <= N + 1, "diis linear: the solution to 1e-9 after "// &
                 to_string(reached)//" steps (n = "//to_string(N)//", limit n + 1), error after "//to_string(N + 4)// &
                 " steps "//to_string(worst, "es9.2")//", restarts "//to_string(diis%restarts))
      call cur%release()
      call nxt%release()
      call diis%release()
      call la%release()
   end subroutine diis_linear

   subroutine symmetrize_vector(no, nv, x)
      !! The t2 block of the concatenated (t1, t2) made pair-symmetric: t2(i,j,a,b) <- (t2(i,j,a,b) + t2(j,i,b,a))/2
      integer, intent(in) :: no, nv
      real(dp), intent(inout) :: x(:)
      real(dp), allocatable :: t2(:, :, :, :)
      integer :: n1

      n1 = no*nv
      t2 = reshape(x(n1 + 1:), [no, no, nv, nv])
      t2 = 0.5_dp*(t2 + reshape(t2, [no, no, nv, nv], order=[2, 1, 4, 3]))
      x(n1 + 1:) = reshape(t2, [size(t2)])
   end subroutine symmetrize_vector

   function two_point(w1, w2) result(x)
      !! The DIIS combination of two entries whose old amplitudes were zero (each error is its vector)
      real(dp), intent(in) :: w1(:), w2(:)
      real(dp) :: x(size(w1))
      real(dp) :: c1

      c1 = (dot_product(w2, w2) - dot_product(w1, w2))/(dot_product(w1, w1) - 2.0_dp*dot_product(w1, w2) + &
                                                           dot_product(w2, w2))
      x = c1*w1 + (1.0_dp - c1)*w2
   end function two_point

   subroutine diis_packed(error)
      !! A pair-symmetric linear fixed-point problem driven through a full and a packed history in lockstep: the
      !! same weights and the same extrapolated vector to round-off, the packed one exactly pair-symmetric, and
      !! the solution within (packed length) + 1 steps
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NO1 = 2, NV1 = 2, NOV = NO1*NV1, N = NOV + NOV*NOV, NP = NOV + NOV*(NOV + 1)/2
      type(trc_linalg_t) :: la
      type(trc_rccsd_diis_t) :: full, half
      type(trc_rccsd_amps_t) :: cur, nxt, out_full, out_half
      real(dp) :: a(N, N), b(N), x(N), xs(N), xn(N), xf(N)
      real(dp) :: dvec, dw, dnorm, dmax, gap, worst
      integer :: i, j, k, reached
      integer(int64) :: seed
      logical :: ex_full, ex_half, ok

      seed = 4242_int64
      do j = 1, N
         do i = 1, N
            a(i, j) = 0.7_dp*lcg(seed)/real(N, dp)
         end do
         b(j) = lcg(seed)
      end do
      xs = 0.0_dp
      do k = 1, 400
         xs = matmul(a, xs) + b
         call symmetrize_vector(NO1, NV1, xs)
      end do
      call la%init(2)
      call full%init(NO1, NV1, N + 2)
      call half%init(NO1, NV1, N + 2, packed=.true.)
      call cur%init(NO1, NV1)
      call nxt%init(NO1, NV1)
      call out_full%init(NO1, NV1)
      call out_half%init(NO1, NV1)
      ok = half%n == NP .and. full%n == N .and. half%packed .and. .not. full%packed
      x = 0.0_dp
      reached = -1
      dvec = 0.0_dp
      dw = 0.0_dp
      dnorm = 0.0_dp
      dmax = 0.0_dp
      gap = 0.0_dp
      worst = huge(1.0_dp)
      do k = 1, NP + 2
         xn = matmul(a, x) + b
         call symmetrize_vector(NO1, NV1, xn)
         call pack_vector(x, cur)
         call pack_vector(xn, nxt)
         call full%push(la, cur, nxt, .true.)
         call half%push(la, cur, nxt, .true.)
         dnorm = max(dnorm, abs(full%step_norm2 - half%step_norm2)/full%step_norm2)
         dmax = max(dmax, abs(full%step_max - half%step_max))
         call full%next(la, nxt, out_full, ex_full)
         call half%next(la, nxt, out_half, ex_half)
         ! at the solution B turns singular and both restart, at the same step
         ok = ok .and. (ex_full .eqv. ex_half) .and. full%restarts == half%restarts
         call out_full%pull()
         call out_half%pull()
         call unpack_vector(out_full, xf)
         call unpack_vector(out_half, x)
         dvec = max(dvec, maxval(abs(x - xf)))
         gap = max(gap, pair_gap(out_half%t2))
         ! compare the weights while the errors are well above round-off: at the solution B is singular
         if (ex_full .and. worst > 1.0e-6_dp) then
            dw = max(dw, maxval(abs(full%coef - half%coef))/maxval(abs(full%coef)))
         end if
         worst = maxval(abs(x - xs))
         if (reached < 0 .and. worst < 1.0e-9_dp) reached = k
      end do
      ok = ok .and. dvec < 1.0e-12_dp .and. dw < 1.0e-10_dp .and. dnorm < 1.0e-13_dp .and. dmax < 1.0e-15_dp .and. &
           gap == 0.0_dp .and. reached > 0 .and. reached <= NP + 1
      call check(error, ok, "diis packed: packed minus full, extrapolated vector "//to_string(dvec, "es9.2")// &
                 ", weights (relative) "//to_string(dw, "es9.2")//", ||step||_2 (relative) "//to_string(dnorm, "es9.2")// &
                 ", max|step| "//to_string(dmax, "es9.2")//"; packed output pair asymmetry "//to_string(gap, "es9.2")// &
                 "; solution to 1e-9 after "//to_string(reached)//" steps (packed length "//to_string(NP)//"), restarts "// &
                 to_string(full%restarts)//" and "//to_string(half%restarts))
      call cur%release()
      call nxt%release()
      call out_full%release()
      call out_half%release()
      call full%release()
      call half%release()
      call la%release()
   end subroutine diis_packed

   subroutine diis_host(error)
      !! The pair-symmetric linear fixed-point problem of `diis_packed` driven through a device and a host history in
      !! lockstep, packed with room for every step and full with a ring of 3 that wraps: the step measures bitwise
      !! equal (both computed on the device), the weights and the extrapolated vector the same to round-off (dots on
      !! the same vectors; the host sums the extrapolation), the same restarts, and the packed pair converges
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NO1 = 2, NV1 = 2, NOV = NO1*NV1, N = NOV + NOV*NOV, NP = NOV + NOV*(NOV + 1)/2
      type(trc_linalg_t) :: la
      type(trc_rccsd_diis_t) :: dev, host
      type(trc_rccsd_amps_t) :: cur, nxt, out_dev, out_host
      real(dp) :: a(N, N), b(N), x(N), xs(N), xn(N), xd(N)
      real(dp) :: dvec, dw, worst
      integer :: i, j, k, cfg, nsize, reached
      integer(int64) :: seed
      logical :: ex_dev, ex_host, ok, same_step, packed
      character(len=:), allocatable :: text

      seed = 4242_int64
      do j = 1, N
         do i = 1, N
            a(i, j) = 0.7_dp*lcg(seed)/real(N, dp)
         end do
         b(j) = lcg(seed)
      end do
      xs = 0.0_dp
      do k = 1, 400
         xs = matmul(a, xs) + b
         call symmetrize_vector(NO1, NV1, xs)
      end do
      call la%init(2)
      call cur%init(NO1, NV1)
      call nxt%init(NO1, NV1)
      call out_dev%init(NO1, NV1)
      call out_host%init(NO1, NV1)
      ok = .true.
      text = "diis host:"
      do cfg = 1, 2
         packed = cfg == 1
         nsize = merge(N + 2, 3, packed)
         call dev%init(NO1, NV1, nsize, packed=packed)
         call host%init(NO1, NV1, nsize, packed=packed, host=.true.)
         ok = ok .and. host%host .and. .not. dev%host
         x = 0.0_dp
         reached = -1
         dvec = 0.0_dp
         dw = 0.0_dp
         same_step = .true.
         worst = huge(1.0_dp)
         do k = 1, NP + 2
            xn = matmul(a, x) + b
            call symmetrize_vector(NO1, NV1, xn)
            call pack_vector(x, cur)
            call pack_vector(xn, nxt)
            call dev%push(la, cur, nxt, .true.)
            call host%push(la, cur, nxt, .true.)
            same_step = same_step .and. dev%step_norm2 == host%step_norm2 .and. dev%step_max == host%step_max
            call dev%next(la, nxt, out_dev, ex_dev)
            call host%next(la, nxt, out_host, ex_host)
            ok = ok .and. (ex_dev .eqv. ex_host) .and. dev%restarts == host%restarts .and. dev%nstored == host%nstored
            call out_dev%pull()
            call out_host%pull()
            call unpack_vector(out_dev, xd)
            call unpack_vector(out_host, x)
            dvec = max(dvec, maxval(abs(x - xd)))
            ! the weights while the errors are well above round-off: at the solution B is singular
            if (ex_dev .and. worst > 1.0e-6_dp) dw = max(dw, maxval(abs(dev%coef - host%coef))/maxval(abs(dev%coef)))
            worst = maxval(abs(x - xs))
            if (reached < 0 .and. worst < 1.0e-9_dp) reached = k
         end do
         ok = ok .and. same_step .and. dvec < 1.0e-12_dp .and. dw < 1.0e-10_dp
         if (packed) ok = ok .and. reached > 0 .and. reached <= NP + 1
         text = text//merge(" packed, ring of ", " full, ring of   ", packed)//to_string(nsize)// &
                ": host minus device, extrapolated vector "//to_string(dvec, "es9.2")//", weights (relative) "// &
                to_string(dw, "es9.2")//", step measures bitwise "//merge("yes", "no ", same_step)// &
                ", solution to 1e-9 after "//to_string(reached)//" steps, restarts "//to_string(host%restarts)//";"
         call dev%release()
         call host%release()
      end do
      call check(error, ok, text)
      call cur%release()
      call nxt%release()
      call out_dev%release()
      call out_host%release()
      call la%release()
   end subroutine diis_host

   subroutine diis_singular(error)
      !! Two identical steps make B singular: next() must restart, report, and return the plain new amplitudes. The
      !! restarted history then extrapolates again from two new entries. Full and packed histories.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NO1 = 1, NV1 = 2, N = NO1*NV1 + NO1*NO1*NV1*NV1
      type(trc_linalg_t) :: la
      type(trc_rccsd_diis_t) :: diis
      type(trc_rccsd_amps_t) :: cur, nxt
      real(dp) :: v(N), w(N), x(N)
      integer :: k, nstored_before, mode
      logical :: ex, ok, ok_after
      real(dp) :: d_after
      character(len=:), allocatable :: text

      ok = .true.
      text = ""
      call la%init(2)
      do mode = 1, 2
         do k = 1, N
            v(k) = 0.1_dp*real(k, dp)
            w(k) = 0.1_dp*real(k, dp) + 0.03_dp*real(k*k, dp)
         end do
         if (mode == 2) then
            call symmetrize_vector(NO1, NV1, v)
            call symmetrize_vector(NO1, NV1, w)
         end if
         call diis%init(NO1, NV1, 4, packed=mode == 2)
         call cur%init(NO1, NV1)
         call nxt%init(NO1, NV1)
         call pack_vector(v, nxt)
         call diis%push(la, cur, nxt, .true.)
         call diis%next(la, nxt, cur, ex)
         ok = ok .and. .not. ex .and. diis%nstored == 1
         call cur%init(NO1, NV1)                ! back to the same old amplitudes
         call diis%push(la, cur, nxt, .true.)   ! the same step again: parallel to the first
         nstored_before = diis%nstored
         call diis%next(la, nxt, cur, ex)
         call cur%pull()
         call unpack_vector(cur, x)
         ok = ok .and. nstored_before == 2 .and. .not. ex .and. diis%restarts == 1 .and. diis%nstored == 0 .and. &
              diis%status%has_error() .and. maxval(abs(x - v)) < 1.0e-14_dp
         text = text//" ["//merge("packed", "full  ", mode == 2)//": restarts = "//to_string(diis%restarts)// &
                ", history length after "//to_string(diis%nstored)//", output differs from the plain new amplitudes by "// &
                to_string(maxval(abs(x - v)), "es9.2")
         ! after the restart: two entries that are not parallel extrapolate again, to the two-point combination
         call cur%init(NO1, NV1)
         call diis%push(la, cur, nxt, .true.)
         call diis%next(la, nxt, cur, ex)
         ok_after = .not. ex .and. diis%nstored == 1
         call cur%init(NO1, NV1)
         call pack_vector(w, nxt)
         call diis%push(la, cur, nxt, .true.)
         call diis%next(la, nxt, cur, ex)
         call cur%pull()
         call unpack_vector(cur, x)
         d_after = maxval(abs(x - two_point(v, w)))
         ok_after = ok_after .and. ex .and. diis%nstored == 2 .and. diis%restarts == 1 .and. d_after < 1.0e-13_dp
         ok = ok .and. ok_after
         text = text//"; then two new entries extrapolate = "//merge("T", "F", ex)//", off the two-point "// &
                "combination by "//to_string(d_after, "es9.2")//"]"
         call cur%release()
         call nxt%release()
         call diis%release()
      end do
      call check(error, ok, "diis singular:"//text)
      call la%release()
   end subroutine diis_singular

   subroutine diis_nan(error)
      !! A NaN in the new amplitudes: measured as not finite, never appended; the history it left still extrapolates.
      !! Full and packed histories.
      type(error_type), allocatable, intent(out) :: error

      type(trc_linalg_t) :: la
      type(trc_rccsd_diis_t) :: diis
      type(trc_rccsd_amps_t) :: cur, nxt
      real(dp) :: v(2*3 + 36), w(2*3 + 36), x(2*3 + 36)
      real(dp) :: d_after
      integer :: mode
      logical :: ok, ex
      character(len=:), allocatable :: text

      ok = .true.
      text = ""
      call la%init(2)
      do mode = 1, 2
         call diis%init(2, 3, 4, packed=mode == 2)
         call cur%init(2, 3)
         call nxt%init(2, 3)
         nxt%t2 = 0.5_dp
         call nxt%push()
         call unpack_vector(nxt, v)
         call diis%push(la, cur, nxt, .true.)
         ok = ok .and. diis%finite .and. diis%nstored == 1
         nxt%t2(2, 1, 3, 2) = ieee_nan()
         call nxt%push()
         call diis%push(la, cur, nxt, .true.)
         ok = ok .and. .not. diis%finite .and. diis%status%has_error() .and. diis%nstored == 1 .and. &
              diis%step_max > 1.0e300_dp
         nxt%t2(2, 1, 3, 2) = huge(1.0_dp)   ! an overflow of the squared norm is just as unusable
         call nxt%push()
         call diis%push(la, cur, nxt, .true.)
         ok = ok .and. .not. diis%finite .and. diis%nstored == 1
         ! a finite step after them: the entry kept from before the bad steps is intact
         nxt%t2 = 0.5_dp
         nxt%t1 = 0.1_dp
         call nxt%push()
         call unpack_vector(nxt, w)
         call diis%push(la, cur, nxt, .true.)
         call diis%next(la, nxt, cur, ex)
         call cur%pull()
         call unpack_vector(cur, x)
         d_after = maxval(abs(x - two_point(v, w)))
         ok = ok .and. diis%finite .and. ex .and. diis%nstored == 2 .and. d_after < 1.0e-13_dp
         text = text//" ["//merge("packed", "full  ", mode == 2)//": "//to_string(diis%nstored)// &
                " entries after a finite step, extrapolation off the two-point combination by "// &
                to_string(d_after, "es9.2")//"]"
         call cur%release()
         call nxt%release()
         call diis%release()
      end do
      call check(error, ok, "diis NaN: a NaN or an overflowing step is flagged not finite and never appended:"//text)
      call la%release()
   end subroutine diis_nan

   subroutine diis_step(error)
      !! The step measure of a known vector, spread over several chunks; full and packed histories (the packed one on
      !! a pair-symmetric vector whose largest element is off the diagonal, where the stored step is scaled by
      !! sqrt(2)). A packed vector above 2^31-1 elements is refused.
      type(error_type), allocatable, intent(out) :: error

      type(trc_linalg_t) :: la
      type(trc_rccsd_diis_t) :: diis
      type(trc_rccsd_amps_t) :: cur, nxt
      real(dp), allocatable :: v(:)
      integer :: n, k, mode
      integer(int64) :: seed
      real(dp) :: vmax, vnorm
      logical :: ok, ok_mode
      character(len=:), allocatable :: text

      ok = .true.
      text = ""
      call la%init(2)
      n = 4*19 + 4*4*19*19
      allocate (v(n))
      do mode = 1, 2
         call diis%init(4, 19, 0, packed=mode == 2)
         call cur%init(4, 19)
         call nxt%init(4, 19)
         seed = 777_int64
         do k = 1, n
            v(k) = 0.001_dp*lcg(seed)
         end do
         if (mode == 1) then
            v(n) = -0.75_dp   ! the maximum sits in the last chunk
         else
            call symmetrize_vector(4, 19, v)
            ! t2(3,4,19,19) and its partner t2(4,3,19,19), off the diagonal, in the last chunk
            v(4*19 + 3 + 4*3 + 16*18 + 16*19*18) = -0.75_dp
            v(4*19 + 4 + 4*2 + 16*18 + 16*19*18) = -0.75_dp
         end if
         call pack_vector(v, nxt)
         call diis%push(la, cur, nxt, .true.)
         vmax = maxval(abs(v))
         vnorm = sqrt(sum(v*v))
         ok_mode = diis%finite .and. abs(diis%step_max - vmax) < 1.0e-15_dp .and. &
                   abs(diis%step_norm2 - vnorm) < 1.0e-13_dp .and. diis%nstored == 0 .and. diis%npart > 1
         ok = ok .and. ok_mode
         text = text//" ["//merge("packed", "full  ", mode == 2)//": max "//to_string(diis%step_max, "es12.5")// &
                " (expected "//to_string(vmax, "es12.5")//"), 2-norm "//to_string(diis%step_norm2, "es12.5")// &
                " (expected "//to_string(vnorm, "es12.5")//"), "//to_string(diis%npart)//" chunks]"
         call cur%release()
         call nxt%release()
         call diis%release()
      end do
      ! o v = 70000: o v (o v + 1)/2 alone is 2.45e9
      call diis%init(1, 70000, 0, packed=.true.)
      ok_mode = diis%status%has_error() .and. .not. allocated(diis%err) .and. &
                index(diis%status%get_message(), "2^31-1") > 0
      ok = ok .and. ok_mode
      text = text//" [packed o 1 v 70000 refused: "//merge("T", "F", ok_mode)//"]"
      call diis%release()
      call check(error, ok, "diis step, no history with size 0:"//text)
      call la%release()
   end subroutine diis_step

end module check_rccsd_run_suite

program check_rccsd_run
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_run_suite, only: collect_check_rccsd_run
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_run, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_run
