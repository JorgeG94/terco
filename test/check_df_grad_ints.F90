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
! derivative kernels were generated for (DFG_LMAX/DFG_LMAX_AUX): prints
! "SKIP: ..." and `stop 0` before run_testsuite is even entered, same
! convention as check_nosym.F90 -- ctest's SKIP_REGULAR_EXPRESSION
! (test/CMakeLists.txt) turns that into "Skipped", not "Passed".
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

module check_df_grad_ints_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_2c
   use trc_df_grad_kernels, only: df_grad_3c_plan, df_grad_3c_driver, df_grad_2c_driver, DFG_LMAX, DFG_LMAX_AUX, &
                                  DFG_NPCLS
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix
   use check_df_grad_ints_m, only: raw_3c, raw_2c, synth_y, synth_z, ao_atoms
   implicit none
   private
   public :: collect_df_grad_ints_tests, check_df_grad_ints_prepare, check_df_grad_ints_cleanup

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: REL_TOL = 1.0e-10_dp
   !! Relative to the largest gradient component (O(10-100) for these
   !! synthetic O(1) coefficients summed over every integral).

   character(len=256) :: g_xyzfile, g_basis_json, g_aux_json, g_dump, g_case_name
   type(trc_basis_t) :: g_bas, g_aux
   type(trc_pairlist_t) :: g_pl
   logical :: g_prepared = .false.
   logical :: g_computed = .false.
   integer :: g_nat, g_nao, g_naux
   integer, allocatable :: g_at_z(:)
   real(dp), allocatable :: g_r(:, :)
   real(dp), allocatable :: g_ip1(:, :, :, :), g_gbsw(:, :, :, :), g_ip2c(:, :, :)
   real(dp), allocatable :: g_y(:, :, :), g_z(:, :), g_g3(:, :), g_g2(:, :), g_g3h(:, :), g_g2h(:, :)
   real(dp) :: g_dsw, g_inv2anti, g_d3, g_d2, g_inv3, g_inv2trans, g_scale3, g_scale2
   real(dp) :: g_dmax1, g_dmax2
   logical :: g_shape_ok

contains

   subroutine check_df_grad_ints_prepare(xyzfile, basis_json, aux_json, dump, case_name, skip)
      !! Reads the command-line-selected molecule/bases and checks the
      !! DFG_LMAX/DFG_LMAX_AUX ceiling, before any test-drive suite is
      !! built -- the program calls this and self-skips (print SKIP, stop
      !! 0) on skip=.true. without ever entering run_testsuite.
      character(*), intent(in) :: xyzfile, basis_json, aux_json, dump, case_name
      logical, intent(out) :: skip
      type(error_t) :: err

      g_xyzfile = xyzfile; g_basis_json = basis_json; g_aux_json = aux_json
      g_dump = dump; g_case_name = case_name

      call read_xyz(trim(g_xyzfile), g_nat, g_at_z, g_r)
      ! prune=.false.: the pyscf reference is built from this file's literal
      ! (unpruned) columns, so the AO-level comparison stays on them.
      call trc_basis_from_json(trim(g_basis_json), g_nat, g_at_z, g_r, g_bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(prune=.false.))
      if (err%has_error()) then
         print '(a)', 'check_df_grad_ints: '//err%get_message(); stop 1
      end if
      ! prune=.false.: the pyscf reference is built from this file's literal
      ! (unpruned) columns, so the AO-level comparison stays on them.
      call trc_basis_from_json(trim(g_aux_json), g_nat, g_at_z, g_r, g_aux, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(prune=.false.))
      if (err%has_error()) then
         print '(a)', 'check_df_grad_ints: '//err%get_message(); stop 1
      end if
      skip = maxval(g_bas%sh_l) > DFG_LMAX .or. maxval(g_aux%sh_l) > DFG_LMAX_AUX
      g_prepared = .not. skip
   end subroutine check_df_grad_ints_prepare

   subroutine ensure_computed()
      !! The expensive, shared computation (raw tensors, device-driver
      !! contraction, the pyscf dump and comparison), memoised: the first
      !! test to run triggers it, every later test in this process reuses
      !! the cached results -- test-drive's `check()` deallocates its
      !! intent(out) `error` on entry, so several independent assertions
      !! over the SAME computed state must be separate new_unittest
      !! entries, not sequential checks in one, or an earlier failure
      !! would be silently overwritten by the next check() call.
      integer :: mu, nu, p, q, unit
      integer, allocatable :: aoat(:), axat(:)
      real(dp), allocatable :: omat(:, :), jmat(:, :)
      real(dp), allocatable :: cq0(:), d0(:, :)
      integer :: coff(0:DFG_NPCLS)
      integer, allocatable :: q_a(:), q_b(:), q_o(:), q_n(:)
      real(dp), allocatable :: pd(:, :)

      if (g_computed) return

      g_nao = g_bas%nao; g_naux = g_aux%nao
      call g_bas%to_device()
      call g_aux%to_device()
      call g_pl%build(g_bas, THRESH)
      call g_pl%to_device()

      ! --- normalisation probes, as check_df_vs_pyscf.F90 ---
      allocate (omat(g_nao, g_nao), jmat(g_naux, g_naux))
      !$acc enter data create(omat, jmat)
      call trc_df_2c(g_bas, omat)
      call trc_df_2c(g_aux, jmat)
      !$acc update self(omat, jmat)
      !$acc exit data delete(omat, jmat)

      ! --- raw derivative tensors, host, one-hot ---
      allocate (g_ip1(3, g_nao, g_nao, g_naux), g_gbsw(3, g_nao, g_nao, g_naux), g_ip2c(3, g_naux, g_naux))
      call raw_3c(g_bas, g_aux, g_ip1, g_gbsw)
      call raw_2c(g_aux, g_ip2c)
      g_dsw = maxval(abs(g_ip1 - g_gbsw))
      g_inv2anti = 0.0_dp
      do q = 1, g_naux
         do p = 1, g_naux
            g_inv2anti = max(g_inv2anti, maxval(abs(g_ip2c(:, p, q) + g_ip2c(:, q, p))))
         end do
      end do

      ! --- contracted, on the device drivers ---
      allocate (g_y(g_nao, g_nao, g_naux), g_z(g_naux, g_naux))
      allocate (g_g3(3, g_nat), g_g2(3, g_nat), g_g3h(3, g_nat), g_g2h(3, g_nat))
      call synth_y(g_nao, g_naux, g_y)
      call synth_z(g_naux, g_z)
      g_g3 = 0.0_dp; g_g2 = 0.0_dp
      ! Y alone (the drivers' K-style dense input, read as M + M^T, so y/2
      ! for a symmetric y): the J part's c and D zero.
      allocate (cq0(g_naux), d0(g_nao, g_nao))
      cq0 = 0.0_dp; d0 = 0.0_dp
      call df_grad_3c_plan(g_bas%nshell, g_bas%maxnp, g_bas%sh_l, g_bas%sh_np, g_bas%sh_e, g_bas%sh_c, g_bas%sh_r, &
                           g_pl%npair, g_pl%p_i, g_pl%p_j, coff, q_a, q_b, q_o, q_n, pd)
      g_y = 0.5_dp*g_y
      !$acc enter data copyin(g_y, g_z, g_g3, g_g2, cq0, d0, q_a, q_b, q_o, q_n, pd)
      call df_grad_3c_driver(g_bas%nshell, g_nao, g_bas%sh_ao, g_bas%sh_at, g_bas%sh_r, coff, g_pl%npair, &
                             q_a, q_b, q_o, q_n, size(pd, 2), pd, g_aux%nshell, g_aux%maxnp, g_aux%sh_l, &
                             g_aux%sh_np, g_aux%sh_ao, g_aux%sh_at, g_aux%sh_e, g_aux%sh_c, g_aux%sh_r, &
                             g_nat, 2, g_nao, g_naux, g_y, 0.0_dp, g_naux, 0, cq0, d0, g_g3)
      call df_grad_2c_driver(g_aux%nshell, g_naux, g_aux%maxnp, g_aux%sh_l, g_aux%sh_np, g_aux%sh_ao, g_aux%sh_at, &
                             g_aux%sh_e, g_aux%sh_c, g_aux%sh_r, g_nat, g_z, g_g2)
      !$acc update self(g_g3, g_g2)
      !$acc exit data delete(g_y, g_z, g_g3, g_g2, cq0, d0, q_a, q_b, q_o, q_n, pd)
      g_y = 2.0_dp*g_y

      ! --- the same contraction from the raw host tensors ---
      allocate (aoat(g_nao), axat(g_naux))
      call ao_atoms(g_bas, aoat)
      call ao_atoms(g_aux, axat)
      g_g3h = 0.0_dp; g_g2h = 0.0_dp
      do p = 1, g_naux
         do nu = 1, g_nao
            do mu = 1, g_nao
               ! A on mu's atom, B on nu's (via the symmetric partner), C on P's.
               g_g3h(:, aoat(mu)) = g_g3h(:, aoat(mu)) + g_y(mu, nu, p)*g_ip1(:, mu, nu, p)
               g_g3h(:, aoat(nu)) = g_g3h(:, aoat(nu)) + g_y(mu, nu, p)*g_ip1(:, nu, mu, p)
               g_g3h(:, axat(p)) = g_g3h(:, axat(p)) - g_y(mu, nu, p)*(g_ip1(:, mu, nu, p) + g_ip1(:, nu, mu, p))
            end do
         end do
      end do
      do q = 1, g_naux
         do p = 1, g_naux
            g_g2h(:, axat(p)) = g_g2h(:, axat(p)) + g_z(p, q)*g_ip2c(:, p, q)
            g_g2h(:, axat(q)) = g_g2h(:, axat(q)) - g_z(p, q)*g_ip2c(:, p, q)
         end do
      end do
      g_scale3 = max(1.0_dp, maxval(abs(g_g3h)))
      g_scale2 = max(1.0_dp, maxval(abs(g_g2h)))
      g_d3 = maxval(abs(g_g3 - g_g3h))
      g_d2 = maxval(abs(g_g2 - g_g2h))
      g_inv3 = maxval(abs(sum(g_g3, dim=2)))
      g_inv2trans = maxval(abs(sum(g_g2, dim=2)))

      open (newunit=unit, file=trim(g_dump), access='stream', form='unformatted', status='replace')
      write (unit) int(g_nat, int64), int(g_nao, int64), int(g_naux, int64)
      write (unit) int(g_at_z, int64)
      write (unit) g_r
      write (unit) omat
      write (unit) jmat
      write (unit) g_ip1
      write (unit) g_ip2c
      write (unit) g_y
      write (unit) g_z
      write (unit) g_g3
      write (unit) g_g2
      close (unit)
      print '(a)', 'check_df_grad_ints: wrote '//trim(g_dump)

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      block
         integer :: u, nat_ref, nao_ref, naux_ref
         real(dp), allocatable :: g3_ref(:, :), g2_ref(:, :)
         call ref_require('dfg_ints_'//trim(g_case_name)//'.ref', u)
         call ref_read_int(u, 'NAT', nat_ref)
         call ref_read_int(u, 'NAO', nao_ref)
         call ref_read_int(u, 'NAUX', naux_ref)
         g_shape_ok = nat_ref == g_nat .and. nao_ref == g_nao .and. naux_ref == g_naux
         if (.not. g_shape_ok) then
            call ref_close(u)
            g_computed = .true.
            return
         end if
         allocate (g3_ref(3, g_nat), g2_ref(3, g_nat))
         call ref_read_matrix(u, 'G3', 3, g_nat, g3_ref)
         call ref_read_matrix(u, 'G2', 3, g_nat, g2_ref)
         call ref_close(u)
         g_dmax1 = maxval(abs(g_g3 - g3_ref))
         g_dmax2 = maxval(abs(g_g2 - g2_ref))
         print '(a,2es10.2)', '  (gradient scales, 3c and 2c: ', g_scale3, g_scale2
      end block
      g_computed = .true.
   end subroutine ensure_computed

   subroutine collect_df_grad_ints_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("raw derivative identities (B==A swap, 2c antisymmetry)", &
                               test_raw_identities), &
                  new_unittest("device 3c/2c drivers match the host raw contraction", &
                               test_device_vs_host), &
                  new_unittest("3c/2c contracted gradients are translationally invariant", &
                               test_translational_invariance), &
                  new_unittest("contracted gradients match pyscf's int3c2e_ip/int2c2e_ip1", &
                               test_vs_pyscf) &
                  ]
   end subroutine collect_df_grad_ints_tests

   subroutine test_raw_identities(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, g_dsw < 1.0e-13_dp*max(1.0_dp, maxval(abs(g_ip1))), &
                 "d/dB (mu nu|P) == d/dA (nu mu|P), max|diff| = "//to_string(g_dsw, "es10.2"))
      if (allocated(error)) return
      call check(error, g_inv2anti < 1.0e-13_dp*max(1.0_dp, maxval(abs(g_ip2c))), &
                 "d/dA (P|Q) == -d/dA (Q|P), max|diff| = "//to_string(g_inv2anti, "es10.2"))
   end subroutine test_raw_identities

   subroutine test_device_vs_host(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, g_d3 < 1.0e-12_dp*g_scale3, &
                 "device 3c driver != host raw contraction, max|diff| = "//to_string(g_d3, "es10.2"))
      if (allocated(error)) return
      call check(error, g_d2 < 1.0e-12_dp*g_scale2, &
                 "device 2c driver != host raw contraction, max|diff| = "//to_string(g_d2, "es10.2"))
   end subroutine test_device_vs_host

   subroutine test_translational_invariance(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, g_inv3 < 1.0e-12_dp*g_scale3, &
                 "3c translational invariance residual = "//to_string(g_inv3, "es10.2"))
      if (allocated(error)) return
      call check(error, g_inv2trans < 1.0e-12_dp*g_scale2, &
                 "2c translational invariance residual = "//to_string(g_inv2trans, "es10.2"))
   end subroutine test_translational_invariance

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, g_shape_ok, "check_df_grad_ints: pyscf reference shape does not match this build")
      if (allocated(error)) return
      call check(error, g_dmax1 < REL_TOL*g_scale3, &
                 "3c Y:d(mu nu|P)/dR vs pyscf ip1/ip2, worst = "//to_string(g_dmax1, "es10.2"))
      if (allocated(error)) return
      call check(error, g_dmax2 < REL_TOL*g_scale2, &
                 "2c Z:d(P|Q)/dR vs pyscf ip1, worst = "//to_string(g_dmax2, "es10.2"))
   end subroutine test_vs_pyscf

   subroutine check_df_grad_ints_cleanup()
      if (.not. g_prepared) return
      call g_pl%release()
      call g_bas%release()
      call g_aux%release()
   end subroutine check_df_grad_ints_cleanup

end module check_df_grad_ints_suite

program check_df_grad_ints
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_df_grad_ints_suite, only: collect_df_grad_ints_tests, check_df_grad_ints_prepare, &
                                       check_df_grad_ints_cleanup
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=256) :: xyzfile, basis_json, aux_json, dump, case_name
   logical :: skip

   if (command_argument_count() < 5) then
      print '(a)', 'usage: check_df_grad_ints <xyz> <basis.json> <aux.json> <dump.bin> <case>'; stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basis_json)
   call get_command_argument(3, aux_json)
   call get_command_argument(4, dump)
   call get_command_argument(5, case_name)

   call check_df_grad_ints_prepare(trim(xyzfile), trim(basis_json), trim(aux_json), trim(dump), &
                                   trim(case_name), skip)
   if (skip) then
      print '(a)', '  SKIP: a basis needs more angular momentum than this build''s DF derivative kernels.'
      stop 0
   end if

   stat = 0
   testsuites = [new_testsuite("check_df_grad_ints", collect_df_grad_ints_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   call check_df_grad_ints_cleanup()
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop 1
   end if
end program check_df_grad_ints
