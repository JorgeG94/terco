!
! Harmonic vibrational analysis of a Cartesian Hessian: mass-weight it,
! remove translation and rotation EXACTLY (not by diagonalising the full
! 3N x 3N matrix and guessing which near-zero eigenvalues are TR -- see
! below), and diagonalise what is left.
!
! THE TR SUBSPACE. Six mass-weighted generators (3 translations, 3
! rotations about the centre of mass -- trc_vib_tr_generators), built in
! the input Cartesian frame: translations are trivially independent, and
! a rotation about any one lab-frame axis is a linear combination of
! rotations about the molecule's own principal axes, so the lab-frame
! generators span the SAME subspace a principal-axis construction would
! (pyscf's hessian/thermo.py `_get_TR` uses principal axes; this uses lab
! axes instead, for the same subspace with less bookkeeping). For a
! linear molecule, rotation about the molecular axis moves nothing, so
! one combination of the three lab-frame rotation generators is
! identically zero: the generator MATRIX's rank drops from 6 to 5, which
! an SVD finds directly (a relative threshold on the singular values --
! trc_vib_tr_generators) rather than asking the caller whether the
! molecule is linear.
!
! Once the TR space's own orthonormal basis U is in hand (the SVD's left
! singular vectors above the threshold), its orthogonal complement is an
! eigendecomposition of the projector P = I - U U^T: P's eigenvalues are
! exactly 0 (multiplicity = TR's own dimension) or exactly 1 (multiplicity
! 3N - that), so the eigenvectors past the TR count are a complement basis
! by construction, not by a second threshold (pyscf's own
! `harmonic_analysis` does the equivalent projector-eigh step, with QR in
! place of the SVD above). Diagonalising the mass-weighted Hessian
! RESTRICTED to that complement (its dimension is exactly 3N-5 or 3N-6)
! gives the vibrations alone: no translation/rotation residual ever
! reaches the eigenvalue list, so there is nothing to miscount as
! vibrational or imaginary (mqc's defects D4/D6/D7 -- see
! HESSIAN_VIB_THERMO_SPEC.md -- do not arise because there is only one
! mode set here, not a 3N-long one a caller has to filter).
!
! UNITS. mass in amu, coords_bohr in Bohr, hess in Hartree/Bohr^2 (the
! finite-difference convention row/column = 3*(atom-1) + {x,y,z}, 1
! indexed, same as a later CLI's --hessian). Frequencies in cm^-1, signed:
! negative is imaginary (an unconverged geometry's "frequency", never
! literally imaginary). Reduced masses in amu; force constants in
! Hartree/Bohr^2 (the PHYSICAL constant lambda*mu, not the mass-weighted
! eigenvalue lambda itself -- pyscf's own `force_const_au` field is
! lambda, a different quantity despite the name; this module's
! `force_const_hartree_bohr2` is mu*lambda, matching what a test comparing
! against pyscf must multiply pyscf's `force_const_au` by `reduced_mass`
! to get). Displacements are pyscf's own `norm_mode` convention: L/sqrt(m)
! per atom, NOT further renormalised, so a mode's reduced mass completes
! its normalisation (||displacement(:,:,k)||^2 = 1/reduced_mass_amu(k)).
!
! EIGENSOLVER. pic_syev/pic_gesvd (pic-blas's LAPACK wrappers), the same
! ones trc_linalg.F90 wraps for the device build -- but this module calls
! them directly rather than through trc_linalg_t: every array here is a
! small, purely HOST object (a finite-difference Hessian assembled on the
! host), and trc_linalg_t's own contract requires every argument resident
! on the device under -acc. Routing a host array through it would be a
! silent illegal address under mem:separate, not a wrong number, so this
! module stays off that path entirely -- there is nothing for a GPU build
! to do here.
!
module trc_vib
   use pic_types, only: dp, default_int
   use pic_lapack_interfaces, only: pic_syev, pic_gesvd
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use trc_constants, only: HARTREE_PER_AMU_BOHR2_TO_CM1
   use trc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   implicit none
   private

   public :: trc_vibrational_analysis
   public :: trc_vib_log_table

   real(dp), parameter :: TR_RANK_REL_TOL = 1.0e-8_dp
      !! Relative threshold on the generator SVD's singular values (against
      !! the largest one) that separates a genuine translation/rotation
      !! direction from the zero a linear molecule's third rotation
      !! generator collapses to. Not a physical constant: a numerical
      !! rank-detection tolerance, chosen with wide margin against both the
      !! TR singular values (O(1) after mass-weighting) and double
      !! precision's own noise floor (~1e-16).
   real(dp), parameter :: PROJECTOR_EIGVAL_SPLIT = 0.5_dp
      !! P = I - U U^T is an exact projector: its eigenvalues are 0 (the TR
      !! directions) or 1 (everything else), with no eigenvalue genuinely
      !! between the two. 0.5 just has to land strictly between roundoff
      !! noise around 0 and around 1; it is not tuned to anything.

contains

   !
   ! natm = size(mass). hess is (3*natm, 3*natm), upper triangle read
   ! (pic_syev's convention); the caller need not symmetrise it first.
   ! Every output array is allocated here to exactly n_modes = 3*natm -
   ! n_tr columns (0 for a single atom), overwriting whatever the actual
   ! argument held.
   !
   subroutine trc_vibrational_analysis(mass, coords_bohr, hess, freq_cm1, reduced_mass_amu, &
                                        force_const_hartree_bohr2, displacement, is_linear, n_tr, error)
      real(dp), intent(in) :: mass(:)
      real(dp), intent(in) :: coords_bohr(:, :)
      real(dp), intent(in) :: hess(:, :)
      real(dp), allocatable, intent(out) :: freq_cm1(:)
      real(dp), allocatable, intent(out) :: reduced_mass_amu(:)
      real(dp), allocatable, intent(out) :: force_const_hartree_bohr2(:)
      real(dp), allocatable, intent(out) :: displacement(:, :, :)
      logical, intent(out) :: is_linear
      integer, intent(out) :: n_tr
      type(error_t), intent(inout) :: error

      integer :: natm, ndof, n_modes, k, a, c, info
      integer(default_int) :: linfo
      real(dp), allocatable :: hmw(:, :), u_tr(:, :), bvec(:, :), h_reduced(:, :)
      real(dp), allocatable :: mode_full(:, :), lambda(:)
      real(dp) :: msum

      is_linear = .false.
      n_tr = 0
      natm = size(mass)
      ndof = 3*natm

      if (natm < 1) then
         call error%set(ERROR_VALIDATION, "trc_vibrational_analysis: no atoms")
         return
      end if
      if (size(coords_bohr, 1) /= 3 .or. size(coords_bohr, 2) /= natm) then
         call error%set(ERROR_VALIDATION, "trc_vibrational_analysis: coords_bohr must be (3, "// &
                        to_string(natm)//")")
         return
      end if
      if (size(hess, 1) /= ndof .or. size(hess, 2) /= ndof) then
         call error%set(ERROR_VALIDATION, "trc_vibrational_analysis: hess must be ("// &
                        to_string(ndof)//", "//to_string(ndof)//")")
         return
      end if
      if (any(mass <= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "trc_vibrational_analysis: every mass must be positive")
         return
      end if

      allocate (hmw(ndof, ndof))
      call trc_vib_mass_weight(mass, hess, hmw)

      call trc_vib_tr_generators(mass, coords_bohr, u_tr, n_tr, info)
      if (info /= 0) then
         call error%set(ERROR_GENERIC, "trc_vibrational_analysis: the TR-generator SVD did not converge")
         return
      end if
      is_linear = (natm >= 2 .and. n_tr == 5)

      n_modes = ndof - n_tr
      allocate (freq_cm1(n_modes), reduced_mass_amu(n_modes), force_const_hartree_bohr2(n_modes))
      allocate (displacement(3, natm, n_modes))
      if (n_modes == 0) return   ! a single atom: no vibrations at all

      call trc_vib_complement_basis(u_tr, bvec, info)
      if (info /= 0) then
         call error%set(ERROR_GENERIC, "trc_vibrational_analysis: the TR-projector eigh did not converge")
         return
      end if
      if (size(bvec, 2) /= n_modes) then
         call error%set(ERROR_GENERIC, "trc_vibrational_analysis: TR complement has the wrong dimension "// &
                        "("//to_string(size(bvec, 2))//", wanted "//to_string(n_modes)// &
                        ") -- the generator SVD and the projector eigh disagree")
         return
      end if

      ! h_reduced = bvec^T Hmw bvec, the mass-weighted Hessian restricted to
      ! the vibrational complement; matmul is intrinsic and these matrices
      ! are at most a few hundred wide (a post-processing step on a
      ! finite-difference Hessian, not a kernel), so there is no reduction
      ! to a device routine to reuse here.
      allocate (h_reduced(n_modes, n_modes))
      h_reduced = matmul(transpose(bvec), matmul(hmw, bvec))

      allocate (lambda(n_modes))
      call pic_syev(h_reduced, lambda, jobz="V", uplo="U", info=linfo)
      info = int(linfo)
      if (info /= 0) then
         call error%set(ERROR_GENERIC, "trc_vibrational_analysis: pic_syev did not converge "// &
                        "on the projected Hessian")
         return
      end if

      allocate (mode_full(ndof, n_modes))
      mode_full = matmul(bvec, h_reduced)   ! h_reduced now holds the reduced-space eigenvectors

      do k = 1, n_modes
         freq_cm1(k) = sign(1.0_dp, lambda(k))*sqrt(abs(lambda(k)))*HARTREE_PER_AMU_BOHR2_TO_CM1
         msum = 0.0_dp
         do a = 1, natm
            do c = 1, 3
               msum = msum + mode_full(3*(a - 1) + c, k)**2/mass(a)
            end do
         end do
         reduced_mass_amu(k) = 1.0_dp/msum
         force_const_hartree_bohr2(k) = lambda(k)*reduced_mass_amu(k)
         do a = 1, natm
            do c = 1, 3
               displacement(c, a, k) = mode_full(3*(a - 1) + c, k)/sqrt(mass(a))
            end do
         end do
      end do
   end subroutine trc_vibrational_analysis

   ! Hmw(i,j) = hess(i,j) / sqrt(mass(atom(i)) * mass(atom(j))), atom(i) =
   ! the atom that Cartesian row/column i belongs to.
   subroutine trc_vib_mass_weight(mass, hess, hmw)
      real(dp), intent(in) :: mass(:)
      real(dp), intent(in) :: hess(:, :)
      real(dp), intent(out) :: hmw(:, :)
      integer :: ndof, i, j
      real(dp), allocatable :: inv_sqrt_m(:)
      ndof = size(hess, 1)
      allocate (inv_sqrt_m(ndof))
      do i = 1, ndof
         inv_sqrt_m(i) = 1.0_dp/sqrt(mass((i - 1)/3 + 1))
      end do
      do j = 1, ndof
         do i = 1, ndof
            hmw(i, j) = hess(i, j)*inv_sqrt_m(i)*inv_sqrt_m(j)
         end do
      end do
   end subroutine trc_vib_mass_weight

   !
   ! The six mass-weighted translation/rotation generators (pyscf's `_get_TR`
   ! convention, lab-frame rotations -- see this module's own header), as
   ! columns of D(3*natm, 6):
   !
   !   translation c:  D(3*(a-1)+c, c)      = sqrt(m_a)
   !   rotation about lab axis x (col 4):   (0, sqrt(m_a) r_z, -sqrt(m_a) r_y)
   !   rotation about lab axis y (col 5):   (-sqrt(m_a) r_z, 0, sqrt(m_a) r_x)
   !   rotation about lab axis z (col 6):   (sqrt(m_a) r_y, -sqrt(m_a) r_x, 0)
   !
   ! r_a = coords_bohr(:,a) - centre of mass. Returns u_tr(3*natm, n_tr), the
   ! SVD's left singular vectors above TR_RANK_REL_TOL * the largest
   ! singular value -- n_tr is 6 for a generic (nonlinear, >1-atom)
   ! molecule, 5 for linear (the third rotation column collapses to zero),
   ! 3 for a single atom (all three rotation columns are identically zero,
   ! r_a = 0).
   !
   subroutine trc_vib_tr_generators(mass, coords_bohr, u_tr, n_tr, info)
      real(dp), intent(in) :: mass(:)
      real(dp), intent(in) :: coords_bohr(:, :)
      real(dp), allocatable, intent(out) :: u_tr(:, :)
      integer, intent(out) :: n_tr
      integer, intent(out) :: info
      integer :: natm, ndof, a, nsing
      real(dp) :: rx, ry, rz, sm
      real(dp) :: com(3)
      real(dp), allocatable :: d_mat(:, :), s_vals(:), u_full(:, :)
      integer(default_int) :: linfo

      natm = size(mass)
      ndof = 3*natm
      nsing = min(ndof, 6)

      com = 0.0_dp
      do a = 1, natm
         com = com + mass(a)*coords_bohr(:, a)
      end do
      com = com/sum(mass)

      allocate (d_mat(ndof, 6))
      d_mat = 0.0_dp
      do a = 1, natm
         sm = sqrt(mass(a))
         rx = coords_bohr(1, a) - com(1)
         ry = coords_bohr(2, a) - com(2)
         rz = coords_bohr(3, a) - com(3)
         d_mat(3*(a - 1) + 1, 1) = sm
         d_mat(3*(a - 1) + 2, 2) = sm
         d_mat(3*(a - 1) + 3, 3) = sm
         d_mat(3*(a - 1) + 2, 4) = sm*rz
         d_mat(3*(a - 1) + 3, 4) = -sm*ry
         d_mat(3*(a - 1) + 1, 5) = -sm*rz
         d_mat(3*(a - 1) + 3, 5) = sm*rx
         d_mat(3*(a - 1) + 1, 6) = sm*ry
         d_mat(3*(a - 1) + 2, 6) = -sm*rx
      end do

      allocate (s_vals(nsing), u_full(ndof, nsing))
      call pic_gesvd(d_mat, s_vals, u_full, jobu="S", jobvt="N", info=linfo)
      info = int(linfo)
      if (info /= 0) return

      n_tr = count(s_vals > TR_RANK_REL_TOL*s_vals(1))
      allocate (u_tr(ndof, n_tr))
      u_tr = u_full(:, 1:n_tr)
   end subroutine trc_vib_tr_generators

   !
   ! bvec(3*natm, 3*natm - size(u_tr,2)): an orthonormal basis of u_tr's
   ! orthogonal complement, from the eigendecomposition of the projector
   ! P = I - u_tr u_tr^T. P's eigenvalues are exactly 0 or 1 (to roundoff);
   ! pic_syev returns them ascending, so the complement is simply the
   ! eigenvectors past the TR count -- no second threshold, see this
   ! module's own header.
   !
   subroutine trc_vib_complement_basis(u_tr, bvec, info)
      real(dp), intent(in) :: u_tr(:, :)
      real(dp), allocatable, intent(out) :: bvec(:, :)
      integer, intent(out) :: info
      integer :: ndof, k_tr, n_modes
      real(dp), allocatable :: p_mat(:, :), eigval(:)
      integer :: i, j
      integer(default_int) :: linfo

      ndof = size(u_tr, 1)
      k_tr = size(u_tr, 2)
      n_modes = ndof - k_tr

      allocate (p_mat(ndof, ndof))
      p_mat = -matmul(u_tr, transpose(u_tr))
      do i = 1, ndof
         p_mat(i, i) = p_mat(i, i) + 1.0_dp
      end do

      allocate (eigval(ndof))
      call pic_syev(p_mat, eigval, jobz="V", uplo="U", info=linfo)
      info = int(linfo)
      if (info /= 0) return

      allocate (bvec(ndof, n_modes))
      j = 0
      do i = 1, ndof
         if (eigval(i) > PROJECTOR_EIGVAL_SPLIT) then
            j = j + 1
            if (j <= n_modes) bvec(:, j) = p_mat(:, i)
         end if
      end do
   end subroutine trc_vib_complement_basis

   !
   ! A frequency/reduced-mass/force-constant table through global_logger,
   ! for the CLI's later --hessian report. Imaginary modes (freq_cm1 < 0)
   ! print with an "i" suffix on the magnitude, Gaussian's own convention.
   !
   subroutine trc_vib_log_table(freq_cm1, reduced_mass_amu, force_const_hartree_bohr2)
      real(dp), intent(in) :: freq_cm1(:)
      real(dp), intent(in) :: reduced_mass_amu(:)
      real(dp), intent(in) :: force_const_hartree_bohr2(:)
      character(len=256) :: line
      integer :: k, n_imag

      n_imag = count(freq_cm1 < 0.0_dp)
      write (line, "(a,i0,a,i0,a)") "trc_vib: ", size(freq_cm1), " vibrational mode(s), ", &
         n_imag, " imaginary"
      call global_logger%info(trim(line))
      write (line, "(a)") "  mode   freq/cm-1   reduced mass/amu   force const/(Eh/a0^2)"
      call global_logger%info(trim(line))
      do k = 1, size(freq_cm1)
         if (freq_cm1(k) < 0.0_dp) then
            write (line, "(2x,i4,4x,f10.2,a,2x,f12.6,10x,es12.4)") &
               k, abs(freq_cm1(k)), "i", reduced_mass_amu(k), force_const_hartree_bohr2(k)
         else
            write (line, "(2x,i4,4x,f10.2,1x,2x,f12.6,10x,es12.4)") &
               k, freq_cm1(k), reduced_mass_amu(k), force_const_hartree_bohr2(k)
         end if
         call global_logger%info(trim(line))
      end do
   end subroutine trc_vib_log_table

end module trc_vib
