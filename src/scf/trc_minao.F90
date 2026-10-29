!
! The minimal-ANO starting density, PySCF's `minao` guess.
!
! Ported from metalquicha's mqc_czt_minao (and its own project_density,
! mqc_czt_projection), with mqc_ -> trc_ renames.
!
module trc_minao
   !! A starting density built in a minimal atomic natural orbital basis and
   !! projected into the calculation's own.
   !!
   !! Each atom contributes the leading ANO-RCC contractions for its occupied
   !! shells, occupied with the spherically averaged ground-state configuration:
   !! for each l, floor(n_l / 2(2l+1)) doubly occupied columns and then one open
   !! column holding the remainder spread evenly over its 2l+1 components. That
   !! diagonal density is carried into the target basis by `project_density`:
   !!
   !!     D_B = S_BB^+ S_BA D_A S_AB S_BB^+
   !!
   !! with S_BB^+ the pseudo-inverse `scf_context_orthogonalize` builds from the target
   !! overlap's own eigendecomposition, near-null directions dropped rather
   !! than divided by -- the same object the SCF's own orthogonaliser is. This
   !! is PySCF's `project_dm_nr2nr` with nothing dropped.
   !!
   !! There is no free-atom SCF, so nothing to converge or cache, and the atoms
   !! are neutral: the guess carries sum(Z) electrons whatever the molecule's
   !! charge, and the first diagonalisation places the rest. PySCF does not
   !! renormalise the electron count either.
   !!
   !! The basis is `basis_sets/minao.json`, H-Kr, copied verbatim from
   !! metalquicha's `basis_sets/minao/minao.json` (itself written by
   !! `tools/minao/gen_minao_basis.py` from PySCF's `ano` data). Shipped as an
   !! ordinary Basis Set Exchange JSON file -- `trc_basis_json` already reads
   !! this schema and every other basis set here ships the same way -- rather
   !! than as a generated Fortran table: it is basis data, not a kernel.
   !!
   !! terco's own integral engine is Cartesian throughout; a spherical AO
   !! space is a transform applied once at the boundary (`trc_sph`). This
   !! module follows the same convention `trc_sap` does for its own one-
   !! electron Hamiltonian: everything through the cross overlap is built
   !! Cartesian, over a MERGED basis (the target's shells, then the minimal
   !! one's, over the same atoms), and the one spherical transform -- if the
   !! target basis is spherical -- is applied to that merged overlap in one
   !! shot before slicing out the target-target and target-minimal blocks.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_basis_json, only: trc_basis_from_json, trc_basis_overrides_t
   use trc_sph, only: sph_cart_to_sph
   use trc_linalg, only: trc_linalg_t
   use trc_scf_context, only: scf_context_orthogonalize
   implicit none
   private

   public :: trc_minao_build    !! Target-basis density from the minimal ANO atoms
   public :: minao_occupation     !! Doubly occupied columns and open-shell fraction

   character(len=*), parameter :: MINAO_BASIS_DEFAULT = "basis_sets/minao.json"
      !! Overridable with TRC_MINAO_BASIS (app/main.F90's own --guess minao
      !! route), same convention as trc_sap's SAP_BASIS_DEFAULT/TRC_SAP_BASIS.

   integer, parameter :: MINAO_Z_MAX = 36

   real(dp), parameter :: MINAO_LINDEP_DEFAULT = 1.0e-8_dp
      !! Matches trc_scf_settings's own `lindep_thresh` default: a caller
      !! that passes nothing gets the same floor the SCF's own orthogonaliser
      !! would apply to this same target overlap.

   ! Electrons per angular momentum (s, p, d, f) in the spherically averaged
   ! ground state, from `pyscf.data.elements.NRSRHF_CONFIGURATION` (PySCF
   ! 2.14), copied rather than re-derived so the guess is PySCF's. Its oddities
   ! are kept on purpose: Sc puts its third valence electron in p, Cr is
   ! 4s2 3d4, and Mn-Ni empty 4s.
   integer, parameter :: CONFIGURATION(0:3, MINAO_Z_MAX) = reshape([ &
                                                                   1, 0, 0, 0, &   ! H
                                                                   2, 0, 0, 0, &   ! He
                                                                   3, 0, 0, 0, &   ! Li
                                                                   4, 0, 0, 0, &   ! Be
                                                                   4, 1, 0, 0, &   ! B
                                                                   4, 2, 0, 0, &   ! C
                                                                   4, 3, 0, 0, &   ! N
                                                                   4, 4, 0, 0, &   ! O
                                                                   4, 5, 0, 0, &   ! F
                                                                   4, 6, 0, 0, &   ! Ne
                                                                   5, 6, 0, 0, &   ! Na
                                                                   6, 6, 0, 0, &   ! Mg
                                                                   6, 7, 0, 0, &   ! Al
                                                                   6, 8, 0, 0, &   ! Si
                                                                   6, 9, 0, 0, &   ! P
                                                                   6, 10, 0, 0, &  ! S
                                                                   6, 11, 0, 0, &  ! Cl
                                                                   6, 12, 0, 0, &  ! Ar
                                                                   7, 12, 0, 0, &  ! K
                                                                   8, 12, 0, 0, &  ! Ca
                                                                   8, 13, 0, 0, &  ! Sc
                                                                   8, 12, 2, 0, &  ! Ti
                                                                   8, 12, 3, 0, &  ! V
                                                                   8, 12, 4, 0, &  ! Cr
                                                                   6, 12, 7, 0, &  ! Mn
                                                                   6, 12, 8, 0, &  ! Fe
                                                                   6, 12, 9, 0, &  ! Co
                                                                   6, 12, 10, 0, &  ! Ni
                                                                   7, 12, 10, 0, &  ! Cu
                                                                   8, 12, 10, 0, &  ! Zn
                                                                   8, 13, 10, 0, &  ! Ga
                                                                   8, 14, 10, 0, &  ! Ge
                                                                   8, 15, 10, 0, &  ! As
                                                                   8, 16, 10, 0, &  ! Se
                                                                   8, 17, 10, 0, &  ! Br
                                                                   8, 18, 10, 0 &  ! Kr
                                                                   ], [4, MINAO_Z_MAX])

contains

   pure subroutine minao_occupation(atomic_number, l, n_double, open_fraction)
      !! The occupation pattern of one angular momentum, as PySCF's `frac_occ`
      !!
      !! `n_double` columns hold two electrons per function; the next holds
      !! `open_fraction` per function, which may be zero. Zero and zero for an
      !! element outside H-Kr or an l above f.
      integer, intent(in) :: atomic_number
      integer, intent(in) :: l
      integer, intent(out) :: n_double
      real(dp), intent(out) :: open_fraction   !! Electrons per function, in [0, 2)

      integer :: n_elec, per_shell

      n_double = 0
      open_fraction = 0.0_dp
      if (atomic_number < 1 .or. atomic_number > MINAO_Z_MAX) return
      if (l < 0 .or. l > 3) return
      n_elec = CONFIGURATION(l, atomic_number)
      per_shell = 2*(2*l + 1)
      n_double = n_elec/per_shell
      open_fraction = 2.0_dp*real(n_elec - n_double*per_shell, dp)/real(per_shell, dp)
   end subroutine minao_occupation

   subroutine trc_minao_build(b, density, error, la, threshold, eri_thresh, path)
      !! The minao starting density in `b`'s basis, both spins together
      !!
      !! Refused for an element above Kr, and for a Cartesian target basis
      !! (b%spherical false) on an atom with occupied d, whose six Cartesian
      !! functions are not the five the occupations describe. Ghost centres
      !! (b%at_z <= 0) are given no electrons.
      type(trc_basis_t), intent(in) :: b
      real(dp), allocatable, intent(out) :: density(:, :)   !! (n, n), n = scf_space_nao(b)
      type(error_t), intent(inout) :: error
      type(trc_linalg_t), intent(inout), optional, target :: la
         !! The calculation's cuBLAS/cuSOLVER handles. Absent, a local one is
         !! created and released here.
      real(dp), intent(in), optional :: threshold
         !! The SCF's linear-dependence threshold -- see `scf_context_orthogonalize`.
         !! Absent or non-positive takes MINAO_LINDEP_DEFAULT.
      real(dp), intent(in), optional :: eri_thresh
         !! Shell-pair screening for the merged basis's own overlap pass.
         !! Absent takes 1e-12 (trc_scf_settings's own eri_thresh default).
      character(len=*), intent(in), optional :: path
         !! Overrides MINAO_BASIS_DEFAULT/TRC_MINAO_BASIS, for a caller (a
         !! test) that wants a specific file regardless of environment.

      type(trc_basis_t) :: minimal, merged
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: s_cart(:, :), t_cart(:, :), v_cart(:, :), s_sph(:, :)
      real(dp), allocatable :: occ(:), d_small(:, :)
      real(dp), allocatable :: s_target(:, :), s_cross(:, :)
      real(dp), allocatable :: x(:, :), y(:, :), p(:, :), pd(:, :)
      character(len=:), allocatable :: bpath
      integer :: ia, nb, na, n_double
      logical :: owns_la, use_sph
      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      real(dp) :: lindep_threshold, pairlist_threshold, open_fraction

      if (error%has_error()) return

      do ia = 1, b%natm
         if (b%at_elem(ia) < 1 .or. b%at_elem(ia) > MINAO_Z_MAX) then
            call error%set(ERROR_VALIDATION, "minao guess: the minimal ANO data covers "// &
                           "H-Kr and atom "//to_string(ia)//" has Z="//to_string(b%at_elem(ia)))
            return
         end if
         if (.not. b%spherical) then
            call minao_occupation(b%at_elem(ia), 2, n_double, open_fraction)
            if (n_double > 0 .or. open_fraction > 0.0_dp) then
               call error%set(ERROR_VALIDATION, "minao guess: a Cartesian target basis on an "// &
                              "element with occupied d (Z="//to_string(b%at_elem(ia))//") is "// &
                              "not supported")
               return
            end if
         end if
      end do

      bpath = MINAO_BASIS_DEFAULT
      block
         character(len=512) :: pe
         pe = " "
         call get_environment_variable("TRC_MINAO_BASIS", pe)
         if (len_trim(pe) > 0) bpath = trim(pe)
      end block
      if (present(path)) bpath = path

      ! Same cartesian-ness as the target, so the minimal basis's own d/f
      ! shells (if any atom has them occupied) carry 2l+1 spherical
      ! components, matching what `minao_occupation` describes -- the
      ! refusal above already rules out the one case that would mismatch
      ! (a Cartesian target with occupied d).
      call trc_basis_from_json(bpath, b%natm, b%at_elem, b%at_r, minimal, error, &
                               force_cart=.not. b%spherical, overrides=trc_basis_overrides_t(charges=b%at_z))
      if (error%has_error()) then
         call error%set(error%get_code(), "minao guess: "//error%get_message())
         return
      end if

      use_sph = b%spherical
      call minimal_occupation_vector(minimal, b%at_z, use_sph, occ, error)
      if (error%has_error()) then
         call minimal%release()
         return
      end if
      na = size(occ)
      allocate (d_small(na, na))
      d_small = 0.0_dp
      do ia = 1, na
         d_small(ia, ia) = occ(ia)
      end do

      call merge_bases(b, minimal, merged)

      pairlist_threshold = 1.0e-12_dp
      if (present(eri_thresh)) pairlist_threshold = eri_thresh
      call pl%build(merged, pairlist_threshold)
      call pl%to_device()
      call merged%to_device()
      allocate (s_cart(merged%nao, merged%nao), t_cart(merged%nao, merged%nao), v_cart(merged%nao, merged%nao))
      !$acc enter data create(s_cart, t_cart, v_cart)
      call trc_1e(merged, pl, s_cart, t_cart, v_cart)
      !$acc update self(s_cart)
      !$acc exit data delete(s_cart, t_cart, v_cart)
      call pl%release()

      nb = merge(b%nao_sph, b%nao, use_sph)
      if (use_sph) then
         allocate (s_sph(merged%nao_sph, merged%nao_sph))
         call sph_cart_to_sph(merged, s_cart, s_sph)
         s_target = s_sph(1:nb, 1:nb)
         s_cross = s_sph(1:nb, nb + 1:nb + na)
      else
         s_target = s_cart(1:nb, 1:nb)
         s_cross = s_cart(1:nb, nb + 1:nb + na)
      end if
      call merged%release()
      call minimal%release()

      owns_la = .not. present(la)
      if (present(la)) then
         lp => la
      else
         lp => la_own
      end if
      call lp%init(nb)

      lindep_threshold = MINAO_LINDEP_DEFAULT
      if (present(threshold)) then
         if (threshold > 0.0_dp) lindep_threshold = threshold
      end if

      allocate (x(nb, nb))
      call scf_context_orthogonalize(lp, nb, s_target, x, error, lindep_thresh=lindep_threshold)
      if (error%has_error()) then
         call error%set(error%get_code(), "minao guess: "//error%get_message())
         if (owns_la) call lp%release()
         return
      end if

      ! P = S_BB^+ S_BA = X X^T S_BA; D_B = P D_A P^T.
      allocate (y(nb, na), p(nb, na), pd(nb, na), density(nb, nb))
      y = matmul(transpose(x), s_cross)
      p = matmul(x, y)
      pd = matmul(p, d_small)
      density = matmul(pd, transpose(p))

      if (owns_la) call lp%release()
   end subroutine trc_minao_build

   subroutine merge_bases(target_b, minimal, merged)
      !! One basis carrying both shell lists, target first -- Cartesian
      !! throughout (both `sh_e`/`sh_c` are always Cartesian GTO data
      !! regardless of either basis's own `spherical` flag; `nao_sph`/
      !! `sh_ao_sph` are likewise always computed, see trc_api's
      !! `basis_build`), field by field so the already-normalised
      !! coefficients are not folded through `common_fac_sp` a second time
      !! (the same reason `basis_subset` copies rather than rebuilds).
      type(trc_basis_t), intent(in) :: target_b, minimal
      type(trc_basis_t), intent(out) :: merged

      integer :: nt, ns, maxnp

      nt = target_b%nshell
      ns = minimal%nshell
      maxnp = max(target_b%maxnp, minimal%maxnp)

      merged%nshell = nt + ns
      merged%natm = target_b%natm
      merged%nao = target_b%nao + minimal%nao
      merged%nao_sph = target_b%nao_sph + minimal%nao_sph
      merged%maxnp = maxnp
      merged%spherical = target_b%spherical
      merged%on_device = .false.

      allocate (merged%sh_l(nt + ns), merged%sh_np(nt + ns), merged%sh_ao(nt + ns), merged%sh_ao_sph(nt + ns))
      allocate (merged%sh_e(maxnp, nt + ns), merged%sh_c(maxnp, nt + ns), merged%sh_r(3, nt + ns))
      allocate (merged%sh_at(nt + ns))
      merged%sh_e = 0.0_dp
      merged%sh_c = 0.0_dp

      merged%sh_l(1:nt) = target_b%sh_l
      merged%sh_np(1:nt) = target_b%sh_np
      merged%sh_ao(1:nt) = target_b%sh_ao
      merged%sh_ao_sph(1:nt) = target_b%sh_ao_sph
      merged%sh_e(1:target_b%maxnp, 1:nt) = target_b%sh_e
      merged%sh_c(1:target_b%maxnp, 1:nt) = target_b%sh_c
      merged%sh_r(:, 1:nt) = target_b%sh_r
      merged%sh_at(1:nt) = target_b%sh_at

      merged%sh_l(nt + 1:nt + ns) = minimal%sh_l
      merged%sh_np(nt + 1:nt + ns) = minimal%sh_np
      merged%sh_ao(nt + 1:nt + ns) = minimal%sh_ao + target_b%nao
      merged%sh_ao_sph(nt + 1:nt + ns) = minimal%sh_ao_sph + target_b%nao_sph
      merged%sh_e(1:minimal%maxnp, nt + 1:nt + ns) = minimal%sh_e
      merged%sh_c(1:minimal%maxnp, nt + 1:nt + ns) = minimal%sh_c
      merged%sh_r(:, nt + 1:nt + ns) = minimal%sh_r
      merged%sh_at(nt + 1:nt + ns) = minimal%sh_at

      allocate (merged%at_z(merged%natm), merged%at_r(3, merged%natm), merged%at_elem(merged%natm))
      merged%at_z = target_b%at_z
      merged%at_r = target_b%at_r
      merged%at_elem = target_b%at_elem
   end subroutine merge_bases

   subroutine minimal_occupation_vector(minimal, target_at_z, use_sph, occupation, error)
      !! Per-function occupation of the minimal basis, atom by atom
      !!
      !! Checks that the minao.json file holds exactly the columns the
      !! configuration occupies for each atom, per l -- the generator and
      !! this module's own CONFIGURATION table are two copies of one rule,
      !! and this is what keeps them from drifting apart silently.
      type(trc_basis_t), intent(in) :: minimal
      real(dp), intent(in) :: target_at_z(:)
      logical, intent(in) :: use_sph
      real(dp), allocatable, intent(out) :: occupation(:)   !! (nao or nao_sph)
      type(error_t), intent(inout) :: error

      integer :: ish, iatom, l, z, n_double, first, ncomp, expected
      integer :: seen(0:3)
      real(dp) :: open_fraction, occ

      if (use_sph) then
         allocate (occupation(minimal%nao_sph))
      else
         allocate (occupation(minimal%nao))
      end if
      occupation = 0.0_dp

      do iatom = 1, minimal%natm
         z = minimal%at_elem(iatom)
         seen = 0
         do ish = 1, minimal%nshell
            if (minimal%sh_at(ish) /= iatom) cycle
            l = minimal%sh_l(ish)
            if (l > 3) then
               call error%set(ERROR_VALIDATION, "minao guess: the minimal basis has an l > 3 shell")
               return
            end if
            seen(l) = seen(l) + 1
            call minao_occupation(z, l, n_double, open_fraction)
            if (seen(l) <= n_double) then
               occ = 2.0_dp
            else
               occ = open_fraction
            end if
            ncomp = 2*l + 1
            if (use_sph) then
               first = minimal%sh_ao_sph(ish)
            else
               first = minimal%sh_ao(ish)
            end if
            if (target_at_z(iatom) > 0.0_dp) occupation(first:first + ncomp - 1) = occ
         end do

         do l = 0, 3
            call minao_occupation(z, l, n_double, open_fraction)
            expected = n_double
            if (open_fraction > 0.0_dp) expected = expected + 1
            if (seen(l) /= expected) then
               call error%set(ERROR_VALIDATION, "minao guess: the minimal basis has "// &
                              to_string(seen(l))//" columns of l="//to_string(l)// &
                              " for Z="//to_string(z)//" where the configuration occupies "// &
                              to_string(expected))
               return
            end if
         end do
      end do
   end subroutine minimal_occupation_vector

end module trc_minao
