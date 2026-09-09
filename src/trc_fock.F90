!
! The four-centre Fock build, behind one object.
!
! WHY THIS EXISTS
! ---------------
! The four-centre path and the one-electron/DF API grew separately and speak
! different languages: `fock_bins` wants a `pair_bins_t`, the `hp_*` primitive
! arrays and a Schwarz vector, all built by three different routines with
! different argument conventions, while the API wants an `trc_basis_t`. A
! caller wanting both had to build two containers over one geometry and know
! which routine fills which array.
!
! This closes that. `trc_eri_t` is built from a basis and a threshold, owns
! every four-centre structure, stays device-resident, and exposes one call:
!
!     call eri%build(basis, thresh)
!     do while (.not. converged)
!        call eri%fock(basis, dmat, gmat)      ! G = 2J - K
!        ! diagonalise, DIIS, new density -- the caller's
!     end do
!     call eri%release()
!
! WHAT `fock` RETURNS
! -------------------
! The two-electron part of the RHF Fock matrix,
!
!     vhf = J - K/2        for the TOTAL density D (trace = nelec)
!
! so the caller writes `F = Hcore + vhf` with no factor to remember. NOT the
! full Fock matrix: Hcore does not change between iterations, and folding it
! in here would mean either recomputing it every time or storing it in an
! object with no reason to hold it.
!
! The convention matters and is easy to get wrong. `2J - K` is equally common
! and is what the internal kernel and the benchmarks produce -- it is twice
! this. The first SCF written against this API used it and converged smoothly
! to -98.50 for water instead of -75.99, which is the failure mode worth
! designing out: a factor of two in the two-electron energy does not diverge,
! it converges to the wrong answer.
!
! Internally the kernel accumulates a folded half-matrix at twice scale, so
! vhf = (raw + raw^T)/8. Both factors are applied here; neither is the
! caller's business.
!
module trc_eri
   use trc_boys, only: dp, boys_init
   use trc_tables, only: tables_init
   use trc_cart, only: cart_init
   use trc_batch, only: build_pairs
   use trc_hgp, only: build_pairs_hgp
   use trc_screen, only: schwarz_bounds
   use trc_bins, only: pair_bins_t, build_binned_pairs
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   use trc_binkernel, only: fock_bins
   use trc_api, only: trc_basis_t
   use trc_decontract, only: decon_t, decontract_basis, decon_expand, decon_fold, decon_release, &
                             DECON_CUTOFF_DEFAULT, &
                             decon_expand_host, decon_fold_host
   implicit none
   private
   public :: trc_eri_t

   type :: trc_eri_t
      real(dp) :: thresh = 1.0e-10_dp
      integer :: nbas = 0, nao = 0, nhpp = 0
      type(pair_bins_t) :: bins
      !! HGP primitive-pair data: what the kernels actually read.
      integer,  allocatable :: hp_off(:), hp_n(:)
      real(dp), allocatable :: hp_p(:), hp_r(:, :), hp_ra(:, :), hp_rb(:, :)
      real(dp), allocatable :: hp_c(:)
      !! Shell-block density bound, refreshed per iteration from D.
      real(dp), allocatable :: dsh(:, :)
      integer,  allocatable :: sh_l(:), ao_off(:)
      !! The basis the integrals are actually evaluated over, and the
      !! transform back to the one the caller speaks. On a generally
      !! contracted basis `dec%active` is true, `pb` is the split basis and
      !! `nao_w` is its (larger) function count; the public `nao` stays the
      !! contracted one so no caller sees the difference. On a segmented
      !! basis nothing here is allocated and `nao_w == nao`.
      type(trc_basis_t) :: pb
      type(decon_t) :: dec
      integer :: nao_w = 0
      logical :: on_device = .false.
      integer :: nlaunch = 0
      !> Quartets that survive both screens in the last build, counted only
      !> when fock_resident is asked for it. Diagnostic: it costs a full
      !> extra enumeration pass.
      integer(kind=8) :: nkept = 0
      integer(kind=8) :: nwork = 0
      !! Ranks. Every rank builds the same bins and takes every nranks-th
      !! item of the sorted work list from its rank; the Fock matrix is then
      !! summed across ranks. One rank is rank 0 of 1 and no communicator.
      type(comm_t) :: comm
      integer :: rank = 0, nranks = 1
      logical :: distributed = .false.
   contains
      procedure :: build   => eri_build
      procedure :: fock      => eri_fock
      procedure :: fock_many  => eri_fock_many
      procedure :: fock_nosym => eri_fock_nosym
      procedure :: fock_resident => eri_fock_resident
      procedure :: release => eri_release
   end type trc_eri_t

contains

   !
   ! Everything that depends on the geometry and not on the density.
   !
   ! Three structures, all of which the four-centre path needs and none of
   ! which a caller should have to assemble: the MMD primitive pairs (only so
   ! Schwarz bounds can be formed from them), the Schwarz vector itself, and
   ! the HGP primitive pairs the kernels read. The MMD set is discarded once
   ! the bounds exist.
   !
   subroutine eri_build(this, b, thresh, comm, batch_res)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      !! Ranks to split the quartets over. Bind each rank to its device
      !! (trc_bind_device) BEFORE any object goes to the device.
      type(comm_t), intent(in), optional :: comm
      !! Buckets per decade of Schwarz bound in the pair bins; 1 is decades.
      integer, intent(in), optional :: batch_res

      integer :: bres
      character(len=8) :: bt_env

      this%rank = 0; this%nranks = 1; this%distributed = .false.
      if (present(comm)) then
         this%comm = comm
         this%rank = comm%rank()
         this%nranks = comm%size()
         this%distributed = this%nranks > 1
      end if
      bres = 1
      if (present(batch_res)) bres = batch_res

      call this%release()
      this%thresh = thresh
      this%nao = b%nao

      !
      ! SPLIT A GENERAL CONTRACTION BEFORE ANYTHING IS BUILT.
      !
      ! Everything below -- pairs, Schwarz, the HGP primitive pairs, the bins,
      ! the density bound -- is built over `wb`, which is the split basis when
      ! there is one and the caller's basis otherwise. Only the transform of D
      ! on the way in and G on the way out knows the difference, so every
      ! entry point and every caller gets this for free. See trc_decontract
      ! for why it is worth 2.4x on cc-pVDZ.
      !
      ! TRC_NO_DECONTRACT=1 forces the old path, for bisecting a discrepancy.
      !
      bt_env = ' '
      call get_environment_variable('TRC_NO_DECONTRACT', bt_env)
      if (len_trim(bt_env) == 0) then
         ! TRC_DECON_CUTOFF overrides the diffuse threshold the split rule
         ! uses, so the rule can be swept without a rebuild.
         call decontract_basis(b, this%pb, this%dec, cutoff=decon_cutoff())
      end if
      if (this%dec%active) then
         ! The transform runs inside `fock_resident`, where everything is
         ! device-resident; the maps have to be up there with it.
         ! The containing type first, then its arrays. Mapping a component
         ! without the type that holds it is what nvfortran reports as
         ! "partially present", and the same rule is why this%bins is in the
         ! clause below alongside this%bins%sp_i.
         !$acc enter data copyin(this%dec)
         !$acc enter data copyin(this%dec%f_off, this%dec%f_c, this%dec%f_a, &
         !$acc                   this%dec%t_off, this%dec%t_p, this%dec%t_a)
         this%dec%on_device = .true.
      end if
      if (this%dec%active) then
         call build_structures(this, this%pb, thresh, bres)
      else
         call build_structures(this, b, thresh, bres)
      end if
   end subroutine eri_build

   !
   ! Everything eri_build does once the basis to build over has been chosen.
   !
   !> Diffuse threshold for the decontraction rule; DECON_CUTOFF_DEFAULT
   !> unless TRC_DECON_CUTOFF says otherwise.
   real(dp) function decon_cutoff()
      character(len=32) :: e
      integer :: ios
      e = ' '
      call get_environment_variable('TRC_DECON_CUTOFF', e)
      decon_cutoff = DECON_CUTOFF_DEFAULT
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) decon_cutoff
         if (ios /= 0) decon_cutoff = DECON_CUTOFF_DEFAULT
      end if
   end function decon_cutoff

   subroutine build_structures(this, b, thresh, bres)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      integer, intent(in) :: bres

      integer :: npp, i
      integer,  allocatable :: pp_off(:), pp_n(:)
      real(dp), allocatable :: pp_p(:), pp_r(:, :), pp_c(:), pp_e(:, :)
      real(dp), allocatable :: qs(:), one(:)
      integer,  allocatable :: hp_ki(:), hp_kj(:)
      ! Stage timing, printed only when TRC_BUILD_TIMING is set. `build` is
      ! four routines with very different costs and no way to tell from the
      ! outside which one is the slow one.
      logical :: btime
      integer :: bt_c0, bt_c1, bt_rate
      character(len=8) :: bt_env

      bt_env = ' '
      call get_environment_variable('TRC_BUILD_TIMING', bt_env)
      btime = len_trim(bt_env) > 0
      call system_clock(bt_c0, bt_rate)

      ! All three, not just Boys. The Schwarz bounds go through the MMD path,
      ! which reads the Hermite and Cartesian index tables; without them the
      ! bounds come back identically zero, every pair is pre-screened away and
      ! the Fock matrix is silently empty. That is what happened the first
      ! time this ran.
      call boys_init()
      call tables_init()
      call cart_init()
      this%nbas = b%nshell
      this%nao_w = b%nao

      allocate (this%sh_l(b%nshell), this%ao_off(b%nshell))
      this%sh_l   = b%sh_l
      this%ao_off = b%sh_ao

      ! The basis already folded common_fac_sp into its coefficients, so the
      ! separate per-shell factor these builders expect is unity. Passing the
      ! factor twice is the classic way to be wrong by a smooth number here.
      allocate (one(b%nshell))
      one = 1.0_dp

      call build_pairs(b%nshell, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, thresh, &
                       pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, npp)
      if (btime) then
         call system_clock(bt_c1)
         print '(a,a,f9.3,a)', '  [build] ', 'build_pairs       ', real(bt_c1 - bt_c0, dp)/real(bt_rate, dp), ' s'
         bt_c0 = bt_c1
      end if


      allocate (qs(b%nshell*(b%nshell + 1)/2))
      call schwarz_bounds(b%nshell, npp, b%sh_l, pp_off, pp_n, pp_p, pp_r, &
                          pp_c, pp_e, qs)
      if (btime) then
         call system_clock(bt_c1)
         print '(a,a,f9.3,a)', '  [build] ', 'schwarz_bounds    ', real(bt_c1 - bt_c0, dp)/real(bt_rate, dp), ' s'
         bt_c0 = bt_c1
      end if

      deallocate (pp_off, pp_n, pp_p, pp_r, pp_c, pp_e)

      ! The kernel multiplies each pair by sh_c, so sh_c is the amplitude
      ! the prune judges by. The cutoff sits three orders below the
      ! screening threshold, as the pair prefilter in build_pairs does.
      call build_pairs_hgp(b%nshell, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                           one, b%sh_c, thresh*1.0e-3_dp, &
                           this%hp_off, this%hp_n, this%hp_p, this%hp_r, &
                           this%hp_ra, this%hp_rb, this%hp_c, hp_ki, hp_kj, this%nhpp)
      if (btime) then
         call system_clock(bt_c1)
         print '(a,a,f9.3,a)', '  [build] ', 'build_pairs_hgp   ', real(bt_c1 - bt_c0, dp)/real(bt_rate, dp), ' s'
         bt_c0 = bt_c1
      end if

      ! The contracted kernels fold the coefficients into hp_c and never
      ! ask which primitive a pair came from.
      deallocate (hp_ki, hp_kj)

      call build_binned_pairs(b%nshell, b%sh_l, b%sh_np, b%sh_r, qs, thresh, &
                              this%bins, pp_n=this%hp_n, res=bres)
      if (btime) then
         call system_clock(bt_c1)
         print '(a,a,f9.3,a)', '  [build] ', 'build_binned_pairs', real(bt_c1 - bt_c0, dp)/real(bt_rate, dp), ' s'
         bt_c0 = bt_c1
      end if


      allocate (this%dsh(b%nshell, b%nshell))
      this%dsh = huge(1.0_dp)*1.0e-30_dp

      deallocate (qs, one)

      !$acc enter data copyin(this%sh_l, this%ao_off, this%dsh, &
      !$acc                   this%bins, this%bins%sp_i, this%bins%sp_j, &
      !$acc                   this%bins%sp_q, &
      !$acc                   this%hp_off, this%hp_n, this%hp_p, this%hp_r, &
      !$acc                   this%hp_ra, this%hp_rb, this%hp_c)
      this%on_device = .true.
   end subroutine build_structures

   !
   ! G = 2J - K for a given density.
   !
   ! `dmat` and `gmat` are host arrays; the density goes up and the result
   ! comes back. A device-resident variant belongs here for a caller whose
   ! density already lives on the GPU, but it is a separate entry point rather
   ! than a flag -- confusing the two is a segfault, not a wrong number.
   !
   subroutine eri_fock(this, b, dmat, gmat, k_scale, j_scale, hcore, &
                       density_screen)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)  :: dmat(this%nao, this%nao)
      !> j*J - k*K/2 for the total density, plus Hcore if given.
      real(dp), intent(out) :: gmat(this%nao, this%nao)
      !> Fraction of exact exchange. One is Hartree-Fock and the default; zero
      !> is pure DFT exchange; a hybrid is between. Present so a Kohn-Sham
      !> build needs no second routine, and because the coupled-perturbed
      !> equations pass it explicitly.
      real(dp), intent(in), optional :: k_scale
      !> Fraction of the Coulomb term, default one. Zero is what the second
      !> pass of a range-separated hybrid wants: long-range exchange and
      !> nothing else.
      real(dp), intent(in), optional :: j_scale
      !> Core Hamiltonian. If present the result is the full Fock matrix
      !> rather than its two-electron part -- the caller almost always wants
      !> that, and making them add it themselves is friction for nothing.
      real(dp), intent(in), optional :: hcore(this%nao, this%nao)
      !> Weight the Schwarz bound by the density it multiplies before
      !> screening. Default true, which is what an SCF wants: a quartet whose
      !> density elements are all negligible contributes nothing however large
      !> its integral.
      !>
      !> A COUPLED-PERTURBED SOLVE MUST PASS FALSE. Its density is a trial
      !> vector the solver drives towards zero, so a screen keyed on the
      !> density's magnitude tightens as the solve proceeds -- the operator
      !> stops being the same linear map from one matvec to the next and the
      !> iteration cannot converge. False recovers the plain Schwarz screen,
      !> which depends on the basis alone and leaves the operator fixed.
      logical, intent(in), optional :: density_screen

      real(dp), allocatable :: jmat(:, :, :), kmat(:, :, :), dwork(:, :, :)
      real(dp), allocatable :: dp_w(:, :), gp_w(:, :)
      integer :: n, nw, i, j
      real(dp) :: jfac, kfac

      jfac = 1.0_dp; kfac = 1.0_dp
      if (present(j_scale)) jfac = j_scale
      if (present(k_scale)) kfac = k_scale

      n = this%nao
      nw = this%nao_w
      ! On a split basis the kernels work in the primitive basis: the density
      ! goes up through C D C^T and the result comes back through C^T G_p C,
      ! a few tens of megaflops each. See trc_decontract.
      allocate (jmat(1, nw, nw), kmat(1, nw, nw), dwork(1, nw, nw))
      allocate (dp_w(nw, nw), gp_w(nw, nw))
      if (this%dec%active) then
         call decon_expand_host(this%dec, dmat, dp_w)
      else
         dp_w = dmat
      end if
      dwork(1, :, :) = dp_w
      jmat = 0.0_dp
      kmat = 0.0_dp

      if (screen_on(density_screen)) then
         call density_blocks(this, dwork, 1)
      else
         this%dsh = 1.0_dp
      end if

      !$acc enter data copyin(dwork, jmat, kmat)
      !$acc update device(this%dsh)
      call fock_bins(this%bins, this%nbas, this%nhpp, nw, this%sh_l, &
                     this%ao_off, this%thresh, .false., &
                     jfac, kfac, .false., this%dsh, &
                     this%hp_off, this%hp_n, this%hp_p, this%hp_r, &
                     this%hp_ra, this%hp_rb, this%hp_c, &
                     1, dwork, jmat, kmat, this%rank, this%nranks, this%nlaunch, this%nwork)
      !$acc wait
      !$acc update self(jmat)
      !$acc exit data delete(dwork, jmat, kmat)

      ! Symmetrise in the basis the integrals were formed in, THEN fold. The
      ! other order is a different matrix -- C is not square.
      do j = 1, nw
         do i = 1, nw
            gp_w(i, j) = 0.5_dp*(jmat(1, i, j) + jmat(1, j, i))
         end do
      end do
      if (this%dec%active) then
         call decon_fold_host(this%dec, gp_w, gmat)
      else
         do j = 1, n
            do i = 1, n
               gmat(i, j) = gp_w(i, j)
            end do
         end do
      end if

      ! vhf = (raw + raw^T)/2.
      !
      ! Determined by measurement, not derivation, and the two disagreed: the
      ! benchmark path reports 2J - K as (raw + raw^T)/4, which would make
      ! J - K/2 equal (raw + raw^T)/8. Against gpu4pyscf on the same density
      ! that is uniformly a factor of four too small -- elementwise ratio
      ! 0.2500 at every element, min = median = max, so a pure scale and not
      ! a structural error. The benchmark and this path must differ by two
      ! somewhere in how the folded halves are counted; the constant here is
      ! the one checked against an external code on more than one system.
      if (present(hcore)) then
         do j = 1, n
            do i = 1, n
               gmat(i, j) = gmat(i, j) + hcore(i, j)
            end do
         end do
      end if
      ! hcore, if it was added, is on every rank: sum the two-electron part only.
      if (this%distributed) then
         if (present(hcore)) gmat = gmat - hcore
         call sum_ranks(this, gmat, n*n)
         if (present(hcore)) gmat = gmat + hcore
      end if

      deallocate (jmat, kmat, dwork, dp_w, gp_w)
   end subroutine eri_fock

   !> Sum a matrix over the ranks, in place, every rank receiving the total.
   !>
   !> Then rank 0's copy is broadcast. An allreduce is not obliged to leave
   !> the same bits on every rank -- 1e-13 apart in practice -- and an SCF
   !> iterated from bit-different matrices on different ranks drifts apart,
   !> slowly when it converges and catastrophically when it does not. The
   !> broadcast costs one more nao^2 message and makes the ranks identical.
   subroutine sum_ranks(this, g, n)
      class(trc_eri_t), intent(in) :: this
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call allreduce(this%comm, g, op=MPI_SUM)
      call bcast(this%comm, g, n, 0)
   end subroutine sum_ranks

   !
   ! Largest |D| in each shell-pair block, for the in-kernel density screen.
   !
   ! Over a BATCH the bound is the maximum across every density, so all of them
   ! see exactly the same quartets. That is deliberate and matches mqc: a batch
   ! then gives bit-comparable results to the same densities passed one at a
   ! time, which is what makes the single-density wrapper safe. Screening each
   ! set on its own bound would be tighter and would make the batch disagree
   ! with the loop.
   !
   !
   ! J - K/2 for a density that is NOT symmetric.
   !
   ! Routes through the enumerated kernel, which reads D_mn and D_nm as
   ! distinct and writes J and K separately, so no fold is applied and the
   ! antisymmetric part of K survives. This is what the coupled-perturbed
   ! Hessian needs; mqc calls it `build_fock_direct_nosym`.
   !
   subroutine eri_fock_nosym(this, b, dmat, gmat, k_scale, j_scale, &
                             density_screen)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)  :: dmat(this%nao, this%nao)
      real(dp), intent(out) :: gmat(this%nao, this%nao)
      real(dp), intent(in), optional :: k_scale, j_scale
      !> Weight the Schwarz bound by the density it multiplies before
      !> screening. Default true, which is what an SCF wants: a quartet whose
      !> density elements are all negligible contributes nothing however large
      !> its integral.
      !>
      !> A COUPLED-PERTURBED SOLVE MUST PASS FALSE. Its density is a trial
      !> vector the solver drives towards zero, so a screen keyed on the
      !> density's magnitude tightens as the solve proceeds -- the operator
      !> stops being the same linear map from one matvec to the next and the
      !> iteration cannot converge. False recovers the plain Schwarz screen,
      !> which depends on the basis alone and leaves the operator fixed.
      logical, intent(in), optional :: density_screen

      real(dp), allocatable :: jmat(:, :, :), kmat(:, :, :), dwork(:, :, :)
      real(dp), allocatable :: dp_w(:, :), gp_w(:, :)
      integer :: n, nw, i, j
      real(dp) :: jfac, kfac

      jfac = 1.0_dp; kfac = 1.0_dp
      if (present(j_scale)) jfac = j_scale
      if (present(k_scale)) kfac = k_scale

      n = this%nao
      nw = this%nao_w
      allocate (jmat(1, nw, nw), kmat(1, nw, nw), dwork(1, nw, nw))
      allocate (dp_w(nw, nw), gp_w(nw, nw))
      if (this%dec%active) then
         call decon_expand_host(this%dec, dmat, dp_w)
      else
         dp_w = dmat
      end if
      dwork(1, :, :) = dp_w
      jmat = 0.0_dp; kmat = 0.0_dp
      if (screen_on(density_screen)) then
         call density_blocks(this, dwork, 1)
      else
         this%dsh = 1.0_dp
      end if

      !$acc enter data copyin(dwork, jmat, kmat)
      !$acc update device(this%dsh)
      call fock_bins(this%bins, this%nbas, this%nhpp, nw, this%sh_l, &
                     this%ao_off, this%thresh, .false., &
                     jfac, kfac, .true., this%dsh, &
                     this%hp_off, this%hp_n, this%hp_p, this%hp_r, &
                     this%hp_ra, this%hp_rb, this%hp_c, &
                     1, dwork, jmat, kmat, this%rank, this%nranks, this%nlaunch, this%nwork)
      !$acc wait
      !$acc update self(jmat, kmat)
      !$acc exit data delete(dwork, jmat, kmat)

      ! J and K come out in full and unfolded, so this is the definition.
      do j = 1, nw
         do i = 1, nw
            gp_w(i, j) = jfac*jmat(1, i, j) - 0.5_dp*kfac*kmat(1, i, j)
         end do
      end do
      if (this%dec%active) then
         call decon_fold_host(this%dec, gp_w, gmat)
      else
         do j = 1, n
            do i = 1, n
               gmat(i, j) = gp_w(i, j)
            end do
         end do
      end if
      if (this%distributed) call sum_ranks(this, gmat, n*n)

      deallocate (jmat, kmat, dwork, dp_w, gp_w)
   end subroutine eri_fock_nosym

   !
   ! Fock build with the density and result ALREADY ON THE DEVICE.
   !
   ! `fock` copies the density up and the result back every call. At 832
   ! functions that is 11 MB each way, and the coupled-perturbed equations run
   ! a hundred of them for matrices that never need to be on the host at all.
   !
   ! Here the caller owns the residency:
   !
   !     !$acc enter data copyin(dmat) create(gmat)
   !     do ...
   !        call eri%fock_resident(b, dmat, gmat)   ! nothing crosses
   !     end do
   !     !$acc update self(gmat)
   !     !$acc exit data delete(dmat, gmat)
   !
   ! A SEPARATE ENTRY POINT rather than a flag on `fock`. Passing a host array
   ! to this routine is a segfault, not a wrong number, and a flag makes that
   ! one typo away; two names make it a compile-time choice.
   !
   ! The scratch it still needs -- the folded accumulator and the density-block
   ! bound -- is allocated per call. Those are nao^2, the same order as the
   ! matrices, so a caller doing many builds should expect this to be improved
   ! by hoisting them into the object; it is the next thing to do here.
   !
   subroutine eri_fock_resident(this, b, dmat, gmat, k_scale, j_scale, &
                                density_screen, count_survivors)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      !! Both must already be present on the device.
      real(dp), intent(in)  :: dmat(this%nao, this%nao)
      real(dp), intent(out) :: gmat(this%nao, this%nao)
      real(dp), intent(in), optional :: k_scale, j_scale
      !> Weight the Schwarz bound by the density before screening (default
      !> true). A COUPLED-PERTURBED SOLVE MUST PASS FALSE -- see `eri_fock`.
      logical, intent(in), optional :: density_screen
      !! Also count the quartets that survive screening, into this%nkept.
      !! A full extra pass over the work list, so it is opt-in.
      logical, intent(in), optional :: count_survivors

      real(dp), allocatable :: jmat(:, :, :), kmat(:, :, :), dwork(:, :, :)
      real(dp), allocatable :: dp_w(:, :), gp_w(:, :)
      integer :: n, nw, i, j
      real(dp) :: jfac, kfac

      jfac = 1.0_dp; kfac = 1.0_dp
      if (present(j_scale)) jfac = j_scale
      if (present(k_scale)) kfac = k_scale
      n = this%nao
      nw = this%nao_w

      allocate (jmat(1, nw, nw), kmat(1, nw, nw), dwork(1, nw, nw))
      allocate (dp_w(nw, nw), gp_w(nw, nw))
      !$acc enter data create(jmat, kmat, dwork, dp_w, gp_w)

      ! The density into the basis the integrals are formed in, on the device.
      ! Nothing crosses the bus here -- that is the whole point of `resident`.
      if (this%dec%active) then
         call decon_expand(this%dec, dmat, dp_w)
      else
         do concurrent(j=1:n, i=1:n)
            dp_w(i, j) = dmat(i, j)
         end do
      end if

      ! Reshape on the device: the kernel wants (ndens, nao, nao) and the
      ! caller has (nao, nao). This is a device-to-device copy, not a transfer.
      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(j=1:nw, i=1:nw)
         dwork(1, i, j) = dp_w(i, j)
         jmat(1, i, j) = 0.0_dp
         kmat(1, i, j) = 0.0_dp
      end do

      ! The screening bound is computed ON THE DEVICE. Doing it host-side, as
      ! `fock` does, would mean pulling the whole density down and would undo
      ! the entire point of this routine.
      if (screen_on(density_screen)) then
         call density_blocks_device(this%nbas, nw, this%sh_l, this%ao_off, dp_w, this%dsh)
      else
         do concurrent(j=1:this%nbas, i=1:this%nbas)
            this%dsh(i, j) = 1.0_dp
         end do
      end if

      if (present(count_survivors)) then
         if (count_survivors) then
            call fock_bins(this%bins, this%nbas, this%nhpp, nw, this%sh_l, &
                           this%ao_off, this%thresh, .false., &
                           jfac, kfac, .false., this%dsh, &
                           this%hp_off, this%hp_n, this%hp_p, this%hp_r, &
                           this%hp_ra, this%hp_rb, this%hp_c, &
                           1, dwork, jmat, kmat, this%rank, this%nranks, this%nlaunch, this%nwork, &
                           nkept=this%nkept)
         end if
      end if
      call fock_bins(this%bins, this%nbas, this%nhpp, nw, this%sh_l, &
                     this%ao_off, this%thresh, .false., &
                     jfac, kfac, .false., this%dsh, &
                     this%hp_off, this%hp_n, this%hp_p, this%hp_r, &
                     this%hp_ra, this%hp_rb, this%hp_c, &
                     1, dwork, jmat, kmat, this%rank, this%nranks, this%nlaunch, this%nwork)
      !$acc wait

      ! Fold on the device, straight into the caller's resident result.
      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      !
      ! Symmetrise in the basis the integrals were formed in, THEN fold back.
      ! The other order is a different matrix, because C is not square.
      do concurrent(j=1:nw, i=1:nw)
         gp_w(i, j) = 0.5_dp*(jmat(1, i, j) + jmat(1, j, i))
      end do
      if (this%dec%active) then
         call decon_fold(this%dec, gp_w, gmat)
      else
         do concurrent(j=1:n, i=1:n)
            gmat(i, j) = gp_w(i, j)
         end do
      end if

      !$acc exit data delete(jmat, kmat, dwork, dp_w, gp_w)
      ! Resident on every rank; the sum goes through the host for now.
      if (this%distributed) then
         !$acc update self(gmat)
         call sum_ranks(this, gmat, n*n)
         !$acc update device(gmat)
      end if
      deallocate (jmat, kmat, dwork, dp_w, gp_w)
   end subroutine eri_fock_resident

   !
   ! Shell-block density bound, computed on the device.
   !
   ! The host version reads a density that has to be there anyway. This one
   ! exists so `fock_resident` never pulls the density down: one thread per
   ! shell pair, each scanning its own block, writing `dsh` in place while it
   ! is device-resident.
   !
   !> `density_screen` with its default applied, in one place.
   !>
   !> Written as a function rather than repeated three times because the
   !> default is the DANGEROUS one for a response solve: a caller that forgets
   !> the argument gets a screen that changes with the trial vector, and a
   !> CPHF that will not converge for a reason nothing in the output points
   !> at. One place to read the default is one place to get it wrong.
   pure logical function screen_on(flag)
      logical, intent(in), optional :: flag
      screen_on = .true.
      if (present(flag)) screen_on = flag
   end function screen_on

   subroutine density_blocks_device(nshell, nao, sh_l, sh_ao, dmat, dsh)
      integer,  intent(in)    :: nshell, nao
      integer,  intent(in)    :: sh_l(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: dmat(nao, nao)
      real(dp), intent(inout) :: dsh(nshell, nshell)

      integer :: ij, i, j, a, c, ni, nj
      real(dp) :: m

      ! Every per-iteration scalar named. This one writes the screening bound
      ! that decides which quartets are computed at all, so stating the
      ! locality rather than leaving it to be inferred is worth the line.
      do concurrent(j=1:nshell, i=1:nshell) local(ni, nj, m, a, c)
         ni = (sh_l(i) + 1)*(sh_l(i) + 2)/2
         nj = (sh_l(j) + 1)*(sh_l(j) + 2)/2
         m = 0.0_dp
         do c = 0, nj - 1
            do a = 0, ni - 1
               m = max(m, abs(dmat(sh_ao(i) + a, sh_ao(j) + c)))
            end do
         end do
         dsh(i, j) = m
      end do
   end subroutine density_blocks_device

   !
   ! Shell-block density bound, over the basis the integrals are built in.
   !
   ! That is `this%sh_l`/`this%ao_off`, not the caller's basis: on a split
   ! general contraction the two differ, and bounding the wrong one indexes
   ! past the end of the density it was handed.
   !
   subroutine density_blocks(this, dmat, ndens)
      class(trc_eri_t), intent(inout) :: this
      integer, intent(in) :: ndens
      real(dp), intent(in) :: dmat(ndens, this%nao_w, this%nao_w)
      integer :: i, j, ni, nj, d
      this%dsh = 0.0_dp
      do d = 1, ndens
         do j = 1, this%nbas
            nj = (this%sh_l(j) + 1)*(this%sh_l(j) + 2)/2
            do i = 1, this%nbas
               ni = (this%sh_l(i) + 1)*(this%sh_l(i) + 2)/2
               this%dsh(i, j) = max(this%dsh(i, j), maxval(abs( &
                  dmat(d, this%ao_off(i):this%ao_off(i) + ni - 1, &
                       this%ao_off(j):this%ao_off(j) + nj - 1))))
            end do
         end do
      end do
   end subroutine density_blocks

   !
   ! Many densities, one pass over the integrals.
   !
   ! The point of the whole exercise: in a direct scheme the integral
   ! evaluation dominates and the contractions against it are nearly free, so
   ! contracting one quartet against every density in hand is the difference
   ! between one integral pass and N of them. The coupled-perturbed equations
   ! for dynamic polarizabilities need about a hundred right-hand sides.
   !
   ! The ceiling is one over the digestion fraction -- measured at 24 to 37%
   ! here, so between 2.7x and 4.2x, approached slowly.
   !
   subroutine eri_fock_many(this, b, ndens, dmats, gmats, k_scale, j_scale, &
                            density_screen)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: ndens
      !! (ndens, nao, nao): the density index is FASTEST on purpose. With it
      !! slowest the N atomic updates for one (mu,nu) are nao^2 apart and the
      !! batch degrades past N=4; adjacent, they coalesce.
      real(dp), intent(in)  :: dmats(ndens, this%nao, this%nao)
      real(dp), intent(out) :: gmats(ndens, this%nao, this%nao)
      real(dp), intent(in), optional :: k_scale, j_scale
      !> Weight the Schwarz bound by the density before screening (default
      !> true). A COUPLED-PERTURBED SOLVE MUST PASS FALSE -- see `eri_fock`.
      !>
      !> Over a batch the bound is the maximum across every density, so all of
      !> them see the same quartets whichever way this is set; that is what
      !> keeps a batch bit-comparable with the same densities passed one at a
      !> time.
      logical, intent(in), optional :: density_screen

      real(dp), allocatable :: jmat(:, :, :), kmat(:, :, :), dwork(:, :, :)
      real(dp), allocatable :: dp_w(:, :), gp_w(:, :)
      integer :: n, nw, i, j, d
      real(dp) :: jfac, kfac

      jfac = 1.0_dp; kfac = 1.0_dp
      if (present(j_scale)) jfac = j_scale
      if (present(k_scale)) kfac = k_scale

      n = this%nao
      nw = this%nao_w
      allocate (jmat(ndens, nw, nw), kmat(ndens, nw, nw), dwork(ndens, nw, nw))
      allocate (dp_w(nw, nw), gp_w(nw, nw))
      ! Every density transformed on the way in, every result folded on the
      ! way out. Both are per-density, and both are negligible against the
      ! one pass over the integrals they share.
      do d = 1, ndens
         if (this%dec%active) then
            call decon_expand_host(this%dec, dmats(d, :, :), dp_w)
         else
            dp_w = dmats(d, :, :)
         end if
         dwork(d, :, :) = dp_w
      end do
      jmat = 0.0_dp
      kmat = 0.0_dp

      if (screen_on(density_screen)) then
         call density_blocks(this, dwork, ndens)
      else
         this%dsh = 1.0_dp
      end if

      !$acc enter data copyin(dwork, jmat, kmat)
      !$acc update device(this%dsh)
      call fock_bins(this%bins, this%nbas, this%nhpp, nw, this%sh_l, &
                     this%ao_off, this%thresh, .false., &
                     jfac, kfac, .false., this%dsh, &
                     this%hp_off, this%hp_n, this%hp_p, this%hp_r, &
                     this%hp_ra, this%hp_rb, this%hp_c, &
                     ndens, dwork, jmat, kmat, this%rank, this%nranks, this%nlaunch, this%nwork)
      !$acc wait
      !$acc update self(jmat)
      !$acc exit data delete(dwork, jmat, kmat)

      do d = 1, ndens
         do j = 1, nw
            do i = 1, nw
               gp_w(i, j) = 0.5_dp*(jmat(d, i, j) + jmat(d, j, i))
            end do
         end do
         if (this%dec%active) then
            call decon_fold_host(this%dec, gp_w, gmats(d, :, :))
         else
            do j = 1, n
               do i = 1, n
                  gmats(d, i, j) = gp_w(i, j)
               end do
            end do
         end if
      end do

      if (this%distributed) call sum_ranks(this, gmats, ndens*n*n)
      deallocate (jmat, kmat, dwork, dp_w, gp_w)
   end subroutine eri_fock_many

   subroutine eri_release(this)
      class(trc_eri_t), intent(inout) :: this
      if (this%on_device) then
         !$acc exit data delete(this%sh_l, this%ao_off, this%dsh, &
         !$acc                  this%bins%sp_i, this%bins%sp_j, &
         !$acc                  this%bins%sp_q, this%bins, &
         !$acc                  this%hp_off, this%hp_n, this%hp_p, this%hp_r, &
         !$acc                  this%hp_ra, this%hp_rb, this%hp_c)
         this%on_device = .false.
      end if
      if (allocated(this%hp_off)) deallocate (this%hp_off)
      if (allocated(this%hp_n))   deallocate (this%hp_n)
      if (allocated(this%hp_p))   deallocate (this%hp_p)
      if (allocated(this%hp_r))   deallocate (this%hp_r)
      if (allocated(this%hp_ra))  deallocate (this%hp_ra)
      if (allocated(this%hp_rb))  deallocate (this%hp_rb)
      if (allocated(this%hp_c))   deallocate (this%hp_c)
      if (allocated(this%dsh))    deallocate (this%dsh)
      if (allocated(this%sh_l))   deallocate (this%sh_l)
      if (allocated(this%ao_off)) deallocate (this%ao_off)
      call decon_release(this%dec)
      call this%pb%release()
      this%nbas = 0; this%nao = 0; this%nao_w = 0; this%nhpp = 0
   end subroutine eri_release

end module trc_eri
