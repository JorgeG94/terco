!
! Cartesian multipole moments against pyscf.
!
! WHY ITS OWN SYSTEM, UNCONTRACTED
! --------------------------------
! The first attempt reused the libfint reference system 2 and compared against
! a pyscf rebuild of it. Every component disagreed -- and so did the OVERLAP,
! which terco reproduces from libfint to 2.2e-15. The fault was the rebuild:
! that system has a contracted s shell, and pyscf renormalises contracted
! shells so their self-overlap is one while `gen_reference` does not. Two
! different basis sets, compared.
!
! So this builds its own system with UNCONTRACTED shells only. For a single
! primitive there is nothing to renormalise: pyscf's coefficient and
! `gto_norm` are the same number, and the two codes are then looking at the
! same functions. The overlap is dumped alongside as the control -- if that
! disagrees, the systems differ and nothing else in the file means anything.
!
! Three centres, off-axis, so no accidental symmetry can hide a component
! swap: with everything on the z axis the x and y moments coincide and half
! the possible ordering errors are invisible.
!
program check_mult
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_multipoles, &
                      TRC_NMULT
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_read_matrix3, ref_report
   implicit none

   integer, parameter :: NSH = 9, NATM = 3, MAXNP = 1
   integer  :: sh_l(NSH), sh_np(NSH), i
   real(dp) :: sh_e(MAXNP, NSH), sh_c(MAXNP, NSH), sh_r(3, NSH)
   real(dp) :: at_z(NATM), at_r(3, NATM)
   real(dp), allocatable :: mm(:, :, :), s(:, :), t(:, :), v(:, :)
   type(trc_basis_t) :: b
   type(trc_pairlist_t) :: pl
   integer :: u, nao

   ! deliberately off-axis and unequal
   at_r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
   at_r(:, 2) = [1.3_dp, 0.0_dp, 0.7_dp]
   at_r(:, 3) = [-0.4_dp, 1.1_dp, 1.9_dp]
   at_z = [4.0_dp, 3.0_dp, 2.0_dp]

   ! s, p and d on each centre, one primitive each
   do i = 1, NATM
      sh_l(3*i - 2) = 0; sh_e(1, 3*i - 2) = 0.9_dp + 0.17_dp*real(i, dp)
      sh_l(3*i - 1) = 1; sh_e(1, 3*i - 1) = 1.4_dp + 0.11_dp*real(i, dp)
      sh_l(3*i)     = 2; sh_e(1, 3*i)     = 1.9_dp + 0.13_dp*real(i, dp)
      sh_r(:, 3*i - 2) = at_r(:, i)
      sh_r(:, 3*i - 1) = at_r(:, i)
      sh_r(:, 3*i)     = at_r(:, i)
   end do
   sh_np = 1
   do i = 1, NSH
      sh_c(1, i) = gto_norm(sh_l(i), sh_e(1, i))
   end do

   call b%build(NSH, sh_l, sh_np, sh_e, sh_c, sh_r, NATM, at_z, at_r, MAXNP)
   call pl%build(b, 1.0e-30_dp)
   nao = b%nao

   allocate (mm(nao, nao, TRC_NMULT), s(nao, nao), t(nao, nao), v(nao, nao))
   call b%to_device(); call pl%to_device()
   !$acc enter data create(mm, s, t, v)
   call trc_1e(b, pl, s, t, v)
   call trc_multipoles(b, pl, [0.0_dp, 0.0_dp, 0.0_dp], mm)
   !$acc update self(mm, s)
   !$acc exit data delete(mm, s, t, v)

   open (newunit=u, file='mult_probe.bin', access='stream', &
         form='unformatted', status='replace')
   write (u) int(nao, int64), int(TRC_NMULT, int64)
   write (u) s
   write (u) mm
   close (u)

   print '(a,i0,a,i0)', '  wrote mult_probe.bin: nao ', nao, &
      '   components ', TRC_NMULT
   print '(a,es20.12)', '  sum(S)  ', sum(s)
   call pl%release(); call b%release()

   ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
   block
      integer :: u, nao_ref, nmult_ref, nfail, k
      real(dp), allocatable :: s_ref(:, :), mm_ref(:, :, :)
      real(dp) :: ds, worst
      call ref_require('mult.ref', u)
      call ref_read_int(u, 'NAO', nao_ref)
      call ref_read_int(u, 'NMULT', nmult_ref)
      if (nao_ref /= nao .or. nmult_ref /= TRC_NMULT) then
         print '(a)', 'check_mult: pyscf reference shape does not match this build'; stop 1
      end if
      allocate (s_ref(nao, nao), mm_ref(nao, nao, TRC_NMULT))
      call ref_read_matrix(u, 'S', nao, nao, s_ref)
      call ref_read_matrix3(u, 'MM', nao, nao, TRC_NMULT, mm_ref)
      call ref_close(u)
      nfail = 0
      ds = maxval(abs(s - s_ref))
      call ref_report(nfail, ds < 1.0e-12_dp, 'CONTROL overlap vs pyscf, max|diff|', ds)
      if (ds >= 1.0e-12_dp) then
         print '(a)', '  -> systems differ; nothing below means anything'; stop 1
      end if
      worst = 0.0_dp
      do k = 1, TRC_NMULT
         worst = max(worst, maxval(abs(mm(:, :, k) - mm_ref(:, :, k))))
      end do
      call ref_report(nfail, worst < 1.0e-10_dp, 'multipoles vs pyscf, max|diff| over all 39 components', worst)
      if (nfail > 0) then
         print '(a)', 'check_mult: FAILED against the pyscf reference'; stop 1
      end if
      print '(a)', 'check_mult: multipoles vs pyscf RESULT: PASS'
   end block

contains

   pure real(dp) function gto_norm(l, a)
      !! libcint's primitive normalisation. For ONE primitive this is also what
      !! pyscf uses, which is the whole reason this system is uncontracted.
      integer,  intent(in) :: l
      real(dp), intent(in) :: a
      real(dp) :: nn, gi
      nn = real(2*l + 2, dp)
      gi = 0.5_dp*(2.0_dp*a)**(-(nn + 1.0_dp)/2.0_dp)*gamma((nn + 1.0_dp)/2.0_dp)
      gto_norm = 1.0_dp/sqrt(gi)
   end function gto_norm

end program check_mult
