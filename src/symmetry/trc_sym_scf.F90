module trc_sym_scf
   !! The SCF in the symmetry-adapted basis: per-irrep orthogonalisation and
   !! diagonalisation, orbital labels, symmetrisation of D and F, and of the
   !! gradient. The Fock BUILD is unchanged (full C1, in the AO basis); this
   !! module only changes how a Fock matrix becomes orbitals.
   !!
   !! One block per irrep i, with the SO columns U_i of `trc_sym_ao`:
   !!   S_i = U_i^T S U_i,  D_i = diag(S_i)^(-1/2),  Sn_i = D_i S_i D_i = V w V^T
   !!   X_i = D_i V_k w_k^(-1/2)        (the kept eigenpairs, w_k >= lindep_thresh)
   !!   T_i = U_i X_i                   (nao x nk_i), so T_i^T S T_i = 1
   !! which is canonical orthogonalisation with the same normalise-first rule
   !! and threshold as `scf_context_orthogonalize`, applied block by block. An
   !! iteration is then, per block, Fp_i = T_i^T F T_i, Fp_i = C' e C'^T and
   !! C_i = T_i C'. The matrix `T` is built once per geometry (`prepare`).
   !!
   !! Nothing is dropped in a compact basis, so the eigenvalues are those of
   !! the C1 problem exactly. With functions dropped, the number dropped per
   !! block can differ from the C1 count (the normalisation is by the SO
   !! diagonal, not the AO one); the MOs past `nmo` get a large energy and a
   !! zero vector, as the null shift of the C1 path places them last.
   !!
   !! The merged MO list is sorted by energy (a k-way merge of the ascending
   !! block spectra, ties to the lower irrep), so aufbau is the first `nocc`
   !! columns exactly as in C1. `mo_irrep(m, spin)` is the irrep of MO m.
   !! Fixed per-irrep occupations (`irrep_nocc`, the `irrep_nelec` of pyscf)
   !! occupy the lowest orbitals of each irrep instead.
   !!
   !! One `trc_linalg_t` per block, so each eigensolve has the workspace of
   !! its own size; `release` releases every one, and `prepare` calls it
   !! first and on each of its error paths.
   !!
   !! DEVICE DATA (GPU build). `prepare` maps `t`, `cblk`, `wk`, `fk`,
   !! `epsb`, `imap` and `isgn` with `enter data`, and `release` unmaps them;
   !! `diagonalize` and `symmetrize` expect F, D and the scratch they are
   !! given to be present. The merge index `ord` is mapped for the length of
   !! one call (`data copyin`). Setup (S_i, X_i, T_i) is host code, except
   !! the eigensolve of Sn_i, which runs in an `acc data` region.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use trc_linalg, only: trc_linalg_t
   use trc_sym_group, only: trc_sym_group_t
   use trc_sym_ao, only: trc_sym_ao_t
   use trc_sym_geom, only: trc_sym_geom_t
   use trc_sym_run, only: SYM_SKEL_AUTO
   use trc_basis_types, only: trc_basis_t
   implicit none
   private

   public :: trc_sym_scf_t, sym_scf_grad_symmetrize

   real(dp), parameter :: PAD_ENERGY = 1.0e3_dp
      !! eigenvalue floor given to the MOs that canonical orthogonalisation
      !! dropped, above any real level and never occupied

   type :: trc_sym_scf_t
      type(trc_sym_ao_t) :: ao  !! AO maps, signs and SO basis: set by the caller
      logical :: decline_lindep = .false.
      !! Decline (`trc_scf_result_t%sym_declined`) instead of running when canonical
      !! orthogonalisation drops functions: the per-irrep drop can differ from the
      !! C1 one (see above), and an on/auto request, which promises C1 wherever
      !! symmetry cannot run, must then give the C1 numbers. Read by `trc_scf_run`.
      integer :: skeleton = SYM_SKEL_AUTO
      !! Whether the Fock build is the SKELETON build: the integral engine is
      !! built with the point
      !! group, launches only the orbit representatives, and the SCF
      !! symmetrizes G after the fold. SYM_SKEL_OFF is the full Fock build,
      !! then symmetrization of F (the S4 path), with an engine built without
      !! the group. SYM_SKEL_AUTO, the default, leaves the choice to
      !! `sym_run_skeleton_size_rule` (trc_sym_run), which declines the skeleton
      !! build on the molecules it was measured to lose on. Only read by
      !! `trc_scf_run` with a group of more than one operation, one rank and
      !! four-centre integrals (density fitting has its own path).
      type(trc_sym_geom_t) :: geom
      !! Optional: the detected geometry (group and atom permutations) the
      !! skeleton build needs. Unset (no `geom%perm`), `geometry` derives the
      !! permutations of `ao%group` from the basis centres, so a caller that
      !! only sets `ao` still gets the skeleton build.
      integer, allocatable :: irrep_nocc(:, :)
      !! Optional fixed occupations, (nirrep, nspin): orbitals per irrep, per
      !! spin channel (restricted: doubly occupied orbitals). Unallocated is
      !! aufbau across irreps. Set by the caller before the run. The gradient
      !! routine rebuilds the density by aufbau from the orbitals, so it is
      !! right only when these occupations are the aufbau ones.
      integer, allocatable :: mo_irrep(:, :)
      !! (nao, nspin) irrep of each MO of the last diagonalisation, 0 for a
      !! dropped one. Kept after `release`, for the result.
      integer :: nmo = 0  !! MOs kept (nao - dropped)
      integer :: nlast = 0
      !! Output of `diagonalize` with fixed occupations: the last occupied MO
      !! (its `occ` argument holds the occupations). Kept on the object so the
      !! call stays within the argument limit.
      integer :: nspin = 0
      logical :: ready = .false.
      logical :: fixed = .false.  !! `irrep_nocc` is in force
      integer, allocatable :: nkeep(:)  !! kept orbitals per irrep
      integer, allocatable :: koff(:)  !! first column of irrep i in `t` is koff(i)+1
      real(dp), allocatable :: t(:, :)  !! (nao, nmo) T = U X per block, columns grouped by irrep
      real(dp), allocatable :: cblk(:, :)  !! (nao, nmo) block eigenvectors, same columns
      real(dp), allocatable :: wk(:), fk(:)  !! flat scratch: F T_i and T_i^T F T_i
      real(dp), allocatable :: epsb(:)  !! (nmo) block eigenvalues
      integer, allocatable :: imap(:, :), isgn(:, :)  !! device copies of ao_map, ao_sign
      type(trc_linalg_t), allocatable :: la(:)  !! one handle per irrep block
   contains
      procedure :: is_active => sym_scf_is_active
      procedure :: prepare => sym_scf_prepare
      procedure :: release => sym_scf_release
      procedure :: check_occupations => sym_scf_check_occupations
      procedure :: diagonalize => sym_scf_diagonalize
      procedure :: symmetrize => sym_scf_symmetrize
      procedure :: geometry => sym_scf_geometry
   end type trc_sym_scf_t

contains

   pure logical function sym_scf_is_active(this) result(on)
      !! True when the group has more than the identity: C1 takes the old path.
      class(trc_sym_scf_t), intent(in) :: this
      on = this%ao%nop > 1
   end function sym_scf_is_active

   subroutine sym_scf_geometry(this, bas, geom, err)
      !! The point-group geometry of the skeleton Fock build: `this%geom` when
      !! the caller set it, else the group of `this%ao` with the atom
      !! permutations found by applying each operation's sign mask to the
      !! basis centres (`bas` sits on the standard frame, as for `sym_ao_build`).
      class(trc_sym_scf_t), intent(in) :: this
      type(trc_basis_t), intent(in) :: bas
      type(trc_sym_geom_t), intent(out) :: geom
      type(error_t), intent(out) :: err
      real(dp), parameter :: TOL = 1.0e-6_dp
      real(dp) :: x(3), s(3)
      integer :: k, ia, ja, ic
      logical :: found

      if (allocated(this%geom%perm)) then
         geom = this%geom
         return
      end if
      geom%group = this%ao%group
      allocate (geom%perm(bas%natm, this%ao%nop))
      do k = 1, this%ao%nop
         do ic = 1, 3
            s(ic) = merge(-1.0_dp, 1.0_dp, btest(this%ao%group%ops(k), ic - 1))
         end do
         do ia = 1, bas%natm
            x = s*bas%at_r(:, ia)
            found = .false.
            do ja = 1, bas%natm
               if (maxval(abs(x - bas%at_r(:, ja))) <= TOL) then
                  geom%perm(ia, k) = ja
                  found = .true.
                  exit
               end if
            end do
            if (.not. found) then
               call err%set(ERROR_VALIDATION, "trc_sym_scf: the basis centres are not invariant under the group "// &
                            "operations (atom "//to_string(ia)//", operation "//to_string(k)//")")
               return
            end if
         end do
      end do
   end subroutine sym_scf_geometry

   subroutine sym_scf_prepare(this, smat, nspin, lindep_thresh, err)
      !! Per-block orthogonalisers and device data for a run of `nspin` spin
      !! channels in the overlap `smat` (host, AO space of `this%ao`).
      class(trc_sym_scf_t), intent(inout) :: this
      real(dp), intent(in) :: smat(:, :)
      integer, intent(in) :: nspin
      real(dp), intent(in) :: lindep_thresh
      type(error_t), intent(out) :: err

      integer :: nao, ir, nd, nk, i, j, k, c0, info, nkmax, off
      real(dp), allocatable :: ui(:, :), si(:, :), w(:), d(:), xr(:, :)

      call this%release()
      nao = this%ao%nao
      if (.not. allocated(this%ao%u) .or. size(smat, 1) /= nao .or. size(smat, 2) /= nao) then
         call err%set(ERROR_VALIDATION, "trc_sym_scf: the SO basis does not match the overlap matrix")
         return
      end if
      this%nspin = nspin
      allocate (this%nkeep(this%ao%nirrep), this%koff(this%ao%nirrep), this%la(this%ao%nirrep))
      allocate (this%t(nao, nao))
      this%nkeep = 0
      this%koff = 0
      this%t = 0.0_dp
      c0 = 0
      nkmax = 0
      do ir = 1, this%ao%nirrep
         nd = this%ao%irrep_dim(ir)
         off = this%ao%irrep_off(ir)
         this%koff(ir) = c0
         if (nd == 0) cycle
         allocate (ui(nao, nd), si(nd, nd), w(nd), d(nd), xr(nd, nd))
         ui = this%ao%u(:, off + 1:off + nd)
         si = matmul(transpose(ui), matmul(smat, ui))
         si = 0.5_dp*(si + transpose(si))
         do i = 1, nd
            d(i) = 1.0_dp/sqrt(si(i, i))
         end do
         do j = 1, nd
            do i = 1, nd
               si(i, j) = d(i)*si(i, j)*d(j)
            end do
         end do
         call this%la(ir)%init(nd)
         !$acc data copy(si) create(w)
         call this%la(ir)%syev(si, w, info)
         !$acc update self(w)
         !$acc end data
         if (info /= 0) then
            call err%set(ERROR_GENERIC, "trc_sym_scf: the overlap of irrep "//to_string(ir)// &
                         " did not diagonalise (info="//to_string(info)//")")
            call this%release()
            return
         end if
         nk = 0
         if (lindep_thresh > 0.0_dp) then
            nk = count(w >= lindep_thresh)
         else if (w(1) > 0.0_dp) then
            nk = nd
         end if
         if (nk == 0 .and. lindep_thresh <= 0.0_dp) then
            call err%set(ERROR_GENERIC, "trc_sym_scf: the overlap of irrep "//to_string(ir)// &
                         " is not positive definite")
            call this%release()
            return
         end if
         ! the kept eigenpairs are the last nk (w is ascending)
         k = 0
         do j = nd - nk + 1, nd
            k = k + 1
            do i = 1, nd
               xr(i, k) = d(i)*si(i, j)/sqrt(w(j))
            end do
         end do
         this%nkeep(ir) = nk
         if (nk > 0) this%t(:, c0 + 1:c0 + nk) = matmul(ui, xr(:, 1:nk))
         c0 = c0 + nk
         nkmax = max(nkmax, nk)
         deallocate (ui, si, w, d, xr)
      end do
      this%nmo = c0

      allocate (this%cblk(nao, max(this%nmo, 1)), this%wk(nao*max(nkmax, 1)), this%fk(max(nkmax, 1)**2), &
                this%epsb(max(this%nmo, 1)))
      this%cblk = 0.0_dp
      this%wk = 0.0_dp
      this%fk = 0.0_dp
      this%epsb = 0.0_dp
      this%imap = this%ao%ao_map
      this%isgn = this%ao%ao_sign
      if (allocated(this%mo_irrep)) deallocate (this%mo_irrep)
      allocate (this%mo_irrep(nao, nspin))
      this%mo_irrep = 0
      !$acc enter data copyin(this%t, this%imap, this%isgn) create(this%cblk, this%wk, this%fk, this%epsb)
      this%ready = .true.
   end subroutine sym_scf_prepare

   subroutine sym_scf_release(this)
      !! Every block handle, and the device data. Idempotent. Keeps `ao`,
      !! `irrep_nocc` and `mo_irrep`.
      class(trc_sym_scf_t), intent(inout) :: this
      integer :: i
      if (allocated(this%la)) then
         do i = 1, size(this%la)
            call this%la(i)%release()
         end do
         deallocate (this%la)
      end if
      if (this%ready) then
         !$acc exit data delete(this%t, this%imap, this%isgn, this%cblk, this%wk, this%fk, this%epsb)
      end if
      this%ready = .false.
      if (allocated(this%t)) deallocate (this%t)
      if (allocated(this%imap)) deallocate (this%imap)
      if (allocated(this%isgn)) deallocate (this%isgn)
      if (allocated(this%cblk)) deallocate (this%cblk)
      if (allocated(this%wk)) deallocate (this%wk)
      if (allocated(this%fk)) deallocate (this%fk)
      if (allocated(this%epsb)) deallocate (this%epsb)
      if (allocated(this%nkeep)) deallocate (this%nkeep)
      if (allocated(this%koff)) deallocate (this%koff)
   end subroutine sym_scf_release

   subroutine sym_scf_check_occupations(this, nocc, err)
      !! Fixed per-irrep occupations must add up to the electron count of
      !! each spin and fit in their irrep. Sets `fixed`.
      class(trc_sym_scf_t), intent(inout) :: this
      integer, intent(in) :: nocc(:)  !! occupied orbitals per spin channel
      type(error_t), intent(out) :: err
      integer :: s
      this%fixed = allocated(this%irrep_nocc)
      if (.not. this%fixed) return
      if (size(this%irrep_nocc, 1) /= this%ao%nirrep .or. size(this%irrep_nocc, 2) < this%nspin) then
         call err%set(ERROR_VALIDATION, "trc_sym_scf: irrep_nocc must be (nirrep, nspin)")
         return
      end if
      do s = 1, this%nspin
         if (any(this%irrep_nocc(:, s) < 0) .or. any(this%irrep_nocc(:, s) > this%nkeep)) then
            call err%set(ERROR_VALIDATION, "trc_sym_scf: irrep_nocc is negative or exceeds the irrep size")
            return
         end if
         if (sum(this%irrep_nocc(:, s)) /= nocc(s)) then
            call err%set(ERROR_VALIDATION, "trc_sym_scf: irrep_nocc of spin "//to_string(s)//" sums to "// &
                         to_string(sum(this%irrep_nocc(:, s)))//", not the "//to_string(nocc(s))//" occupied orbitals")
            return
         end if
      end do
   end subroutine sym_scf_check_occupations

   subroutine sym_scf_diagonalize(this, f, spin, omax, eps_s, cmo_s, info, occ)
      !! The orbitals of one spin from the Fock matrix `f` (device resident in the
      !! GPU build): block eigensolves, the energy-ordered merge into `eps_s`
      !! and `cmo_s` (both present), and the labels. With fixed occupations
      !! `occ(m)` is `omax` for an occupied MO and 0 otherwise (host array,
      !! the caller updates the device copy) and `this%nlast` the last occupied.
      class(trc_sym_scf_t), intent(inout) :: this
      real(dp), intent(in) :: f(:, :)
      integer, intent(in) :: spin
      real(dp), intent(in) :: omax
      real(dp), intent(inout) :: eps_s(:), cmo_s(:, :)
      integer, intent(out) :: info
      real(dp), intent(inout) :: occ(:)

      integer :: nao, ir, nk, c0, m, ib, nmo, ni
      integer, allocatable :: ord(:), ptr(:), cnt(:)
      real(dp) :: best

      best = 0.0_dp
      nao = this%ao%nao
      nmo = this%nmo
      info = 0
      this%nlast = 0
      do ir = 1, this%ao%nirrep
         nk = this%nkeep(ir)
         if (nk == 0) cycle
         c0 = this%koff(ir)
         call block_step(this%la(ir), f, this%t(1, c0 + 1), this%wk(1), this%fk(1), &
                         this%cblk(1, c0 + 1), this%epsb(c0 + 1:c0 + nk), ni)
         if (ni /= 0 .and. info == 0) info = ni
      end do
      if (info /= 0) return

      !$acc update self(this%epsb)
      allocate (ord(nmo), ptr(this%ao%nirrep), cnt(this%ao%nirrep))
      ptr = 0
      cnt = 0
      this%mo_irrep(:, spin) = 0
      do m = 1, nmo
         ib = 0
         best = huge(1.0_dp)
         do ir = 1, this%ao%nirrep
            if (ptr(ir) >= this%nkeep(ir)) cycle
            if (ib == 0 .or. this%epsb(this%koff(ir) + ptr(ir) + 1) < best) then
               ib = ir
               best = this%epsb(this%koff(ir) + ptr(ir) + 1)
            end if
         end do
         ptr(ib) = ptr(ib) + 1
         ord(m) = this%koff(ib) + ptr(ib)
         eps_s(m) = best
         this%mo_irrep(m, spin) = ib
         if (this%fixed) then
            cnt(ib) = cnt(ib) + 1
            if (cnt(ib) <= this%irrep_nocc(ib, spin)) then
               occ(m) = omax
               this%nlast = m
            else
               occ(m) = 0.0_dp
            end if
         end if
      end do
      if (nmo < nao) then
         eps_s(nmo + 1:nao) = PAD_ENERGY + abs(best)
         if (this%fixed) occ(nmo + 1:nao) = 0.0_dp
      end if
      !$acc update device(eps_s)
      call merge_columns(nao, nmo, this%cblk, ord, cmo_s)
   end subroutine sym_scf_diagonalize

   subroutine block_step(la, f, tb, wk, fk, cb, eb, info)
      !! One irrep block: Fp = T^T F T, its eigenpairs, C = T C'. The block
      !! size is that of `eb` (nk) and `f` (nao), so they need not travel as
      !! arguments; the matrices have explicit shapes, so the callers can hand
      !! it the start of a contiguous run of columns of the packed arrays.
      type(trc_linalg_t), intent(inout) :: la
      real(dp), intent(in), contiguous :: f(:, :)
      real(dp), intent(inout) :: eb(:)
      real(dp), intent(in) :: tb(size(f, 1), size(eb))
      real(dp), intent(inout) :: wk(size(f, 1), size(eb)), fk(size(eb), size(eb)), cb(size(f, 1), size(eb))
      integer, intent(out) :: info
      call la%gemm_mat(f, tb, wk)
      call la%gemm_mat(tb, wk, fk, transa="T")
      call la%syev(fk, eb, info)
      if (info /= 0) return
      call la%gemm_mat(tb, fk, cb)
   end subroutine block_step

   subroutine merge_columns(nao, nmo, cblk, ord, cmo)
      !! cmo(:, m) = cblk(:, ord(m)) for m <= nmo, zero after. `ord` lives on
      !! the device for this call only.
      integer, intent(in) :: nao, nmo
      real(dp), intent(in) :: cblk(nao, *)
      integer, intent(in) :: ord(nmo)
      real(dp), intent(inout) :: cmo(nao, nao)
      integer :: i, j
      !$acc data copyin(ord)
      do concurrent(j=1:nao, i=1:nao)
         if (j <= nmo) then
            cmo(i, j) = cblk(i, ord(j))
         else
            cmo(i, j) = 0.0_dp
         end if
      end do
      !$acc end data
   end subroutine merge_columns

   subroutine sym_scf_symmetrize(this, n, a, w)
      !! a <- (1/|G|) sum_g R(g) a R(g)^T for a symmetric `a` (a density or a
      !! Fock matrix); `w` is scratch of the same size. Both present.
      class(trc_sym_scf_t), intent(in) :: this
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(n, n), w(n, n)
      call symmetrize_dev(n, this%ao%nop, this%imap, this%isgn, a, w)
   end subroutine sym_scf_symmetrize

   subroutine symmetrize_dev(n, nop, imap, isgn, a, w)
      integer, intent(in) :: n, nop
      integer, intent(in) :: imap(n, nop), isgn(n, nop)
      real(dp), intent(inout) :: a(n, n), w(n, n)
      integer :: i, j, g
      ! R a R^T moves element (i, j) to (map i, map j) with the signs of
      ! i and j; each operation squares to the identity, so the average is
      ! (1/|G|) sum_g sgn(i,g) sgn(j,g) a(map(i,g), map(j,g)). One kernel per
      ! operation, accumulating into `w`.
      do concurrent(i=1:n, j=1:n)
         w(i, j) = a(i, j)
      end do
      do g = 2, nop
         do concurrent(i=1:n, j=1:n)
            w(i, j) = w(i, j) + real(isgn(i, g)*isgn(j, g), dp)*a(imap(i, g), imap(j, g))
         end do
      end do
      do concurrent(i=1:n, j=1:n)
         a(i, j) = w(i, j)/real(nop, dp)
      end do
   end subroutine symmetrize_dev

   subroutine sym_scf_grad_symmetrize(group, perm, grad)
      !! The gradient of a symmetric molecule is equivariant: the force on
      !! the image atom is the sign-flipped force. Averaging over the group,
      !! grad(:, A) <- (1/|G|) sum_g s_g * grad(:, perm(A, g)), removes the
      !! components symmetry forces to zero and the numerical noise in the
      !! rest. `grad` is in the standard frame of the geometry. Host code.
      type(trc_sym_group_t), intent(in) :: group
      integer, intent(in) :: perm(:, :)
      real(dp), intent(inout) :: grad(:, :)
      real(dp), allocatable :: g0(:, :)
      real(dp) :: s(3)
      integer :: a, k, ic
      g0 = grad
      grad = 0.0_dp
      do k = 1, group%order
         do ic = 1, 3
            s(ic) = merge(-1.0_dp, 1.0_dp, btest(group%ops(k), ic - 1))
         end do
         do a = 1, size(grad, 2)
            grad(:, a) = grad(:, a) + s*g0(:, perm(a, k))
         end do
      end do
      grad = grad/real(group%order, dp)
   end subroutine sym_scf_grad_symmetrize

end module trc_sym_scf
