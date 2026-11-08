!
! trc_sym_frame: the exact AO rotation between the standard frame of a
! symmetric run and the user's input frame, checked against brute force.
!
! The rotated basis is built, not derived: the same shells (every l from 0 to
! this build's TERCO_LMAX, uncontracted, several exponents) on the atoms of a
! nonplanar molecule, once at the standard-frame positions r and once at
! s + rot^T r, the positions the input frame sees. For a random rotation
! `rot` and shift `s`, Cartesian and spherical:
!   - the overlap satisfies S_std = D^T S_in D to 1e-12 of its largest
!     element, and the same for the core Hamiltonian T + V of the molecule's
!     own nuclei (so the identity holds for operators that know the atoms,
!     not only for the overlap);
!   - a density transforms as D D_std D^T and keeps Tr(D S) and Tr(D H);
!   - orbitals: C_in = D C_std keeps C^T S C (the MO overlap) unchanged;
! and the cases where the answer is known:
!   - the identity rotation gives the identity matrix;
!   - a quarter turn about z gives a signed permutation of the components
!     (every entry 0 or +-1), exactly;
!   - D (rot) D (rot^T) is the identity (the blocks of the inverse rotation).
!
! A test-drive suite, serial (`parallel=.false.`). The l range follows the
! build, so an LMAX=4 build also covers g shells.
!
module check_sym_frame_suite
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_error, only: error_t
   use trc_sph, only: sph_cart_to_sph
   use trc_tables, only: LMAX
   use trc_sym_frame, only: sym_ao_rot_blocks, sym_ao_rotate_rows, sym_ao_rotate_dmat, &
                            sym_frame_to_input, sym_frame_to_standard
   use trc_test_basis, only: read_xyz, build_aux
   implicit none
   private
   public :: collect_sym_frame

   real(dp), parameter :: TOL = 1.0e-12_dp  !! relative to the largest element; measured ~1e-14
   character(len=4), parameter :: KIND_NAME(2) = [character(len=4) :: "cart", "sph"]

contains

   subroutine collect_sym_frame(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("overlap and core Hamiltonian of a rotated basis: S = D^T S_in D", test_rotated_basis), &
                  new_unittest("density and orbitals keep their traces", test_invariants), &
                  new_unittest("identity, quarter turn and inverse rotations", test_known_rotations), &
                  new_unittest("frame maps of points are inverse", test_points) &
                  ]
   end subroutine collect_sym_frame

   function random_rotation(seed) result(rot)
      !! A proper rotation from a unit quaternion built out of `seed`.
      real(dp), intent(in) :: seed
      real(dp) :: rot(3, 3)
      real(dp) :: q(4), a, b, c, d
      q = [1.0_dp + 0.3_dp*sin(seed), 0.7_dp*cos(2.0_dp*seed), -0.4_dp + 0.2_dp*sin(3.0_dp*seed), 0.55_dp]
      q = q/sqrt(sum(q**2))
      a = q(1); b = q(2); c = q(3); d = q(4)
      rot(1, :) = [a*a + b*b - c*c - d*d, 2.0_dp*(b*c - a*d), 2.0_dp*(b*d + a*c)]
      rot(2, :) = [2.0_dp*(b*c + a*d), a*a - b*b + c*c - d*d, 2.0_dp*(c*d - a*b)]
      rot(3, :) = [2.0_dp*(b*d - a*c), 2.0_dp*(c*d + a*b), a*a - b*b - c*c + d*d]
   end function random_rotation

   subroutine two_bases(rot, shift, spherical, bas_std, bas_in, error)
      !! The same shells on the molecule's standard-frame atoms and on the
      !! atoms the input frame sees.
      real(dp), intent(in) :: rot(3, 3), shift(3)
      logical, intent(in) :: spherical
      type(trc_basis_t), intent(out) :: bas_std, bas_in
      type(error_type), allocatable, intent(out) :: error
      integer :: nat, nsh, maxnp
      integer, allocatable :: elem(:), sh_l(:), sh_np(:)
      real(dp), allocatable :: r(:, :), r_in(:, :), sh_e(:, :), sh_c(:, :), sh_r(:, :)

      call read_xyz('chfcl2.xyz', nat, elem, r)
      r_in = sym_frame_to_input(rot, shift, r)
      call build_aux(nat, elem, r, LMAX, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
      call bas_std%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, real(elem, dp), r, maxnp)
      bas_std%spherical = spherical
      call build_aux(nat, elem, r_in, LMAX, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
      call bas_in%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, real(elem, dp), r_in, maxnp)
      bas_in%spherical = spherical
      call check(error, bas_std%nshell == bas_in%nshell, "the two bases have different shells")
   end subroutine two_bases

   subroutine one_electron(bas, s, h)
      !! S and H = T + V in the AO space of `bas`.
      type(trc_basis_t), intent(inout) :: bas
      real(dp), allocatable, intent(out) :: s(:, :), h(:, :)
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: sc(:, :), tc(:, :), vc(:, :), hc(:, :)
      integer :: n, ns
      n = bas%nao
      ns = merge(bas%nao_sph, bas%nao, bas%spherical)
      call pl%build(bas, 1.0e-14_dp)
      call bas%to_device(); call pl%to_device()
      allocate (sc(n, n), tc(n, n), vc(n, n), hc(n, n), s(ns, ns), h(ns, ns))
      call trc_1e(bas, pl, sc, tc, vc)
      hc = tc + vc
      if (bas%spherical) then
         call sph_cart_to_sph(bas, sc, s)
         call sph_cart_to_sph(bas, hc, h)
      else
         s = sc
         h = hc
      end if
      call pl%release()
   end subroutine one_electron

   function dense_d(bas, dl) result(d)
      !! The dense AO rotation D (columns: the image of each std AO).
      type(trc_basis_t), intent(in) :: bas
      real(dp), intent(in) :: dl(:, :, 0:)
      real(dp), allocatable :: d(:, :)
      integer :: n, i
      n = merge(bas%nao_sph, bas%nao, bas%spherical)
      allocate (d(n, n))
      d = 0.0_dp
      do i = 1, n
         d(i, i) = 1.0_dp
      end do
      call sym_ao_rotate_rows(bas, dl, d)
   end function dense_d

   subroutine test_rotated_basis(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bs, bi
      real(dp), allocatable :: s_std(:, :), h_std(:, :), s_in(:, :), h_in(:, :), dl(:, :, :), d(:, :)
      real(dp) :: rot(3, 3), shift(3), es, eh
      type(error_t) :: err
      integer :: ik, itry
      do itry = 1, 2
         rot = random_rotation(0.9_dp*itry)
         shift = [0.31_dp, -0.77_dp, 1.2_dp*itry]
         do ik = 1, 2
            call two_bases(rot, shift, ik == 2, bs, bi, error)
            if (allocated(error)) return
            call one_electron(bs, s_std, h_std)
            call one_electron(bi, s_in, h_in)
            call sym_ao_rot_blocks(rot, ik == 2, LMAX, dl, err)
            call check(error, .not. err%has_error(), "rot blocks: "//err%get_message())
            if (allocated(error)) return
            d = dense_d(bs, dl)
            es = maxval(abs(matmul(transpose(d), matmul(s_in, d)) - s_std))/maxval(abs(s_std))
            eh = maxval(abs(matmul(transpose(d), matmul(h_in, d)) - h_std))/maxval(abs(h_std))
            print '(a,a,i2,a,es10.2,a,es10.2)', '  ', KIND_NAME(ik), itry, ': |D^T S_in D - S|/|S| ', es, &
               '   |D^T H_in D - H|/|H| ', eh
            call check(error, es < TOL, KIND_NAME(ik)//": overlap is not rotated exactly: "//fmt(es))
            if (allocated(error)) return
            call check(error, eh < TOL, KIND_NAME(ik)//": core Hamiltonian is not rotated exactly: "//fmt(eh))
            if (allocated(error)) return
            call bs%release(); call bi%release()
         end do
      end do
   end subroutine test_rotated_basis

   subroutine test_invariants(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bs, bi
      real(dp), allocatable :: s_std(:, :), h_std(:, :), s_in(:, :), h_in(:, :), dl(:, :, :)
      real(dp), allocatable :: c(:, :), c_in(:, :), p(:, :), p_in(:, :), q(:, :)
      real(dp) :: rot(3, 3), shift(3), t1, t2, t3, t4
      type(error_t) :: err
      integer :: ik, n, i, j
      rot = random_rotation(2.1_dp)
      shift = [-0.5_dp, 0.25_dp, 0.8_dp]
      do ik = 1, 2
         call two_bases(rot, shift, ik == 2, bs, bi, error)
         if (allocated(error)) return
         call one_electron(bs, s_std, h_std)
         call one_electron(bi, s_in, h_in)
         call sym_ao_rot_blocks(rot, ik == 2, LMAX, dl, err)
         n = size(s_std, 1)
         ! a deterministic non-symmetric orbital matrix with a few columns
         allocate (c(n, 7), p(n, n))
         do j = 1, 7
            do i = 1, n
               c(i, j) = sin(0.37_dp*i*j + 0.11_dp*i) + 0.2_dp*cos(1.3_dp*i - 0.7_dp*j)
            end do
         end do
         p = matmul(c, transpose(c))
         c_in = c
         call sym_ao_rotate_rows(bs, dl, c_in)
         p_in = p
         call sym_ao_rotate_dmat(bs, dl, p_in)
         t1 = sum(p*transpose(s_std)); t2 = sum(p_in*transpose(s_in))
         t3 = sum(p*transpose(h_std)); t4 = sum(p_in*transpose(h_in))
         call check(error, abs(t1 - t2) < 1.0e-9_dp*abs(t1), KIND_NAME(ik)//": Tr(D S) changed: "//fmt(abs(t1 - t2)))
         if (allocated(error)) return
         call check(error, abs(t3 - t4) < 1.0e-9_dp*abs(t3), KIND_NAME(ik)//": Tr(D H) changed: "//fmt(abs(t3 - t4)))
         if (allocated(error)) return
         ! the MO overlap C^T S C is the same in both frames
         q = matmul(transpose(c), matmul(s_std, c)) - matmul(transpose(c_in), matmul(s_in, c_in))
         call check(error, maxval(abs(q)) < 1.0e-9_dp*maxval(abs(matmul(transpose(c), matmul(s_std, c)))), &
                    KIND_NAME(ik)//": C^T S C changed: "//fmt(maxval(abs(q))))
         if (allocated(error)) return
         ! the density equals C_in C_in^T
         q = p_in - matmul(c_in, transpose(c_in))
         call check(error, maxval(abs(q)) < 1.0e-9_dp*maxval(abs(p_in)), KIND_NAME(ik)//": D C C^T D^T /= D (C C^T) D^T")
         if (allocated(error)) return
         print '(a,a,a,es10.2,a,es10.2)', '  ', KIND_NAME(ik), ': |Tr(PS) - Tr(P_in S_in)| ', abs(t1 - t2), &
            '   |Tr(PH) - Tr(P_in H_in)| ', abs(t3 - t4)
         deallocate (c, p)
         call bs%release(); call bi%release()
      end do
   end subroutine test_invariants

   subroutine test_known_rotations(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bs, bi
      real(dp), allocatable :: dl(:, :, :), dl_inv(:, :, :), d(:, :), e(:, :)
      real(dp) :: rot(3, 3), id(3, 3), qt(3, 3), tol_exact
      type(error_t) :: err
      integer :: ik, n, i
      id = 0.0_dp
      do i = 1, 3
         id(i, i) = 1.0_dp
      end do
      qt = reshape([0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])  ! 90 degrees about z
      tol_exact = 1.0e-13_dp
      do ik = 1, 2
         call two_bases(id, [0.0_dp, 0.0_dp, 0.0_dp], ik == 2, bs, bi, error)
         if (allocated(error)) return
         n = merge(bs%nao_sph, bs%nao, ik == 2)
         ! identity
         call sym_ao_rot_blocks(id, ik == 2, LMAX, dl, err)
         d = dense_d(bs, dl)
         e = d
         do i = 1, n
            e(i, i) = e(i, i) - 1.0_dp
         end do
         call check(error, maxval(abs(e)) < tol_exact, KIND_NAME(ik)//": the identity rotation is not the identity matrix")
         if (allocated(error)) return
         ! a quarter turn: a signed permutation in the Cartesian space; in the
         ! spherical space too (real harmonics map to +-harmonics, or mix m, -m)
         call sym_ao_rot_blocks(qt, ik == 2, LMAX, dl, err)
         d = dense_d(bs, dl)
         if (ik == 1) then
            call check(error, all(abs(abs(d) - 1.0_dp) < tol_exact .or. abs(d) < tol_exact), &
                       "cart: a quarter turn about z is not a signed permutation")
            if (allocated(error)) return
         end if
         ! the inverse rotation is the transpose of the rotation matrix; its blocks undo D
         call sym_ao_rot_blocks(transpose(qt), ik == 2, LMAX, dl_inv, err)
         e = matmul(dense_d(bs, dl_inv), d)
         do i = 1, n
            e(i, i) = e(i, i) - 1.0_dp
         end do
         call check(error, maxval(abs(e)) < tol_exact, KIND_NAME(ik)//": D(rot^T) D(rot) is not the identity")
         if (allocated(error)) return
         ! a generic rotation likewise
         rot = random_rotation(0.4_dp)
         call sym_ao_rot_blocks(rot, ik == 2, LMAX, dl, err)
         call sym_ao_rot_blocks(transpose(rot), ik == 2, LMAX, dl_inv, err)
         e = matmul(dense_d(bs, dl_inv), dense_d(bs, dl))
         do i = 1, n
            e(i, i) = e(i, i) - 1.0_dp
         end do
         call check(error, maxval(abs(e)) < 1.0e-12_dp, KIND_NAME(ik)//": D(rot^T) D(rot) is not the identity: "// &
                    fmt(maxval(abs(e))))
         if (allocated(error)) return
         call bs%release(); call bi%release()
      end do
   end subroutine test_known_rotations

   subroutine test_points(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: rot(3, 3), shift(3), x(3, 4), y(3, 4)
      rot = random_rotation(1.7_dp)
      shift = [1.0_dp, 2.0_dp, -3.0_dp]
      x = reshape([0.1_dp, 0.2_dp, 0.3_dp, -1.0_dp, 0.5_dp, 2.0_dp, 3.0_dp, -2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [3, 4])
      y = sym_frame_to_standard(rot, shift, sym_frame_to_input(rot, shift, x))
      call check(error, maxval(abs(y - x)) < 1.0e-14_dp, "to_standard(to_input(x)) /= x")
      if (allocated(error)) return
      y = sym_frame_to_input(rot, shift, x)
      call check(error, abs(sum((y(:, 1) - y(:, 2))**2) - sum((x(:, 1) - x(:, 2))**2)) < 1.0e-13_dp, "distances not preserved")
   end subroutine test_points

   function fmt(x) result(t)
      real(dp), intent(in) :: x
      character(len=16) :: t
      write (t, '(es10.2)') x
   end function fmt

end module check_sym_frame_suite

program check_sym_frame
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_sym_frame_suite, only: collect_sym_frame
   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(collect_sym_frame, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_frame
