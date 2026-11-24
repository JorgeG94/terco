module check_mdj_grad_suite
   !!
   !! The McMurchie-Davidson Coulomb gradient (trc_mdj `gradient`, reached through
   !! eri%gradient with eri%grad_j_engine = TRC_J_ENGINE_MDJ) against terco's own Rys
   !! gradient, which forms exactly the same quantity, d/dX Tr(D J[D]) of the
   !! total density, from quartet derivative integrals. The Rys kernels are the
   !! reference (a pyscf anchor for them exists already, check_full_grad_vs_pyscf);
   !! the numpy oracle of the MD algebra, scripts/md/reference_grad.py, agrees
   !! with pyscf's int2e_ip1 to 1e-14.
   !!
   !! Every call names its engine (eri%grad_j_engine), so the checks stay what they say when the
   !! whole suite is run with TRC_MDJ_GRAD=0 (which only changes callers that
   !! leave the engine to the default, now the MD route).
   !!
   !! The basis comes from the command line: argument 1 is a basis JSON, argument 2
   !! (optional) is "split", which forces the general contraction to be split
   !! into primitives (md_coulomb), the working basis the MD engines prefer.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t
   use trc_tables, only: LMAX
   use trc_eri, only: trc_eri_t, trc_fock_options_t, TRC_J_ENGINE_MDJ, TRC_J_ENGINE_BINS, eri_mdj_grad_force_of, &
                      eri_mdj_grad_pays, MDJ_GRAD_MIN_NAO
   use trc_decontract, only: trc_decontract_options_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_mdj_grad_tests

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: THRESH_LOOSE = 1.0e-7_dp
   real(dp), parameter :: AGREE = 3.0e-13_dp
   !! MD against Rys, absolute, on gradients of order 1e1-1e3, at THRESH 1e-14,
   !! i.e. 1e-14 relative to the largest gradient component's scale (the engines
   !! screen against that threshold, not per integral). Observed over
   !! all three bases: MD - Rys <= 1.6e-11 (cc-pVDZ; 8.2e-12 for 6-31G*, UHF
   !! CH2Cl2, max |g| 1.5e3); the tolerance is 6x the worst. Measured with the
   !! engines' own error: at THRESH 1e-12 most of the gap is the Rys gradient's
   !! (CHFCl2 6-31G*, Coulomb only, against each engine at 1e-15: MD 4.3e-11,
   !! Rys 4.0e-10; at 1e-10 MD 1.2e-10, Rys 1.7e-7), so the engines are compared
   !! at a threshold where an MD mistake below 1e-9 would show.
   real(dp), parameter :: SCREENTOL = 1.0e-7_dp
   !! The (P, Q) screen at THRESH_LOOSE against the same build without it
   !! (test_screen): observed <= 5.7e-8 over all molecules and bases (the cut is
   !! thresh*1e-3 per pair, summed over many pairs); the mutations "energy bound
   !! mp_ip instead of mp_gip" and "I'_P scaled by 0.1" gave 1.3e-7 and 2.0e-7
   !! and fail it.
   real(dp), parameter :: INVTOL = 1.0e-9_dp
   real(dp), parameter :: ROTTOL = 1.0e-9_dp
   real(dp), parameter :: FDTOL = 1.0e-7_dp
   real(dp), parameter :: H = 1.0e-4_dp
   integer, parameter :: NMOL = 3
   !! (cc-pVTZ's JSON has no Cl: nmol_for takes the first two molecules there.)
   character(len=*), parameter :: MOLS(NMOL) = [character(len=14) :: "water.xyz", "h2o2.xyz", "ch2cl2.xyz"]

contains

   subroutine collect_mdj_grad_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("MD J gradient equals the Rys J gradient, RHF dense", test_rhf), &
                  new_unittest("MD J gradient equals the Rys J gradient, RHF banded", test_banded), &
                  new_unittest("MD J gradient equals the Rys J gradient, UHF", test_uhf), &
                  new_unittest("derivative screening bound: error scales with the threshold", test_screen), &
                  new_unittest("switch on: full gradient equals the fused one", test_full), &
                  new_unittest("unsupported calls fall back to the fused kernels", test_fallback), &
                  new_unittest("MD J gradient is translation invariant and rotates", test_symmetry), &
                  new_unittest("MD J gradient matches a finite difference of the J energy", test_fd), &
                  new_unittest("TRC_MDJ_GRAD parser accepts only 0 and 1", test_env_parse), &
                  new_unittest("default route: MD without exchange, by size with it", test_size_rule), &
                 new_unittest("rank shares of the MD route sum to the whole gradient", test_rank_shares) &
                  ]
   end subroutine collect_mdj_grad_tests

   subroutine test_size_rule(error)
      !! The default's size rule (eri_mdj_grad_pays): a pure functional takes the
      !! MD route at any size, exchange only from MDJ_GRAD_MIN_NAO functions.
      type(error_type), allocatable, intent(out) :: error
      call check(error, eri_mdj_grad_pays(10, 0.0_dp), "no exchange: MD at any size")
      if (allocated(error)) return
      call check(error, .not. eri_mdj_grad_pays(MDJ_GRAD_MIN_NAO - 1, 1.0_dp), "exchange below the size: fused")
      if (allocated(error)) return
      call check(error, eri_mdj_grad_pays(MDJ_GRAD_MIN_NAO, 0.25_dp), "exchange at the size: MD")
   end subroutine test_size_rule

   subroutine test_env_parse(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, eri_mdj_grad_force_of("0") == 0, 'TRC_MDJ_GRAD="0" must force the fused route')
      if (allocated(error)) return
      call check(error, eri_mdj_grad_force_of("1") == 1, 'TRC_MDJ_GRAD="1" must force the MD route')
      if (allocated(error)) return
      call check(error, eri_mdj_grad_force_of(" 0") == 0, 'TRC_MDJ_GRAD=" 0" is "0" after adjustl')
      if (allocated(error)) return
      call check(error, eri_mdj_grad_force_of("off") == -1, 'TRC_MDJ_GRAD="off" must be ignored')
      if (allocated(error)) return
      call check(error, eri_mdj_grad_force_of("01") == -1, 'TRC_MDJ_GRAD="01" must be ignored')
      if (allocated(error)) return
      call check(error, eri_mdj_grad_force_of("") == -1, "unset TRC_MDJ_GRAD is the default route")
   end subroutine test_env_parse

   integer function nmol_for(basis_json) result(n)
      character(len=*), intent(in) :: basis_json
      n = NMOL
      if (index(basis_json, 'cc-pvtz') > 0) n = 2
   end function nmol_for

   subroutine basis_args(basis_json, split)
      character(len=*), intent(out) :: basis_json
      logical, intent(out) :: split
      character(len=16) :: a2
      basis_json = ''
      a2 = ''
      split = .false.
      if (command_argument_count() >= 1) call get_command_argument(1, basis_json)
      if (command_argument_count() >= 2) call get_command_argument(2, a2)
      split = trim(a2) == "split"
   end subroutine basis_args

   subroutine synth_density(nao, phase, banded, d)
      !! A deterministic symmetric density: dense (every element 0.3-0.5, so the
      !! density screens drop nothing) or banded (decaying with the AO
      !! distance, the shape a real density has, so they do).
      integer, intent(in) :: nao, phase
      logical, intent(in) :: banded
      real(dp), intent(out) :: d(nao, nao)
      integer :: i, j
      do j = 1, nao
         do i = 1, nao
            if (banded) then
               d(i, j) = (0.6_dp + 0.2_dp*cos(real(3*i + 5*j + 7*phase, dp)))*exp(-0.12_dp*real(abs(i - j), dp))
            else
               d(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j + 7*phase, dp))
            end if
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
   end subroutine synth_density

   subroutine setup(xyz, basis_json, split, bas, eri, nat, at_z, r, error, thresh, rot)
      !! Geometry, basis and engine; SKIPs a basis the build has no kernels for.
      character(*), intent(in) :: xyz, basis_json
      logical, intent(in) :: split
      type(trc_basis_t), intent(out) :: bas
      type(trc_eri_t), intent(out) :: eri
      integer, intent(out) :: nat
      integer, allocatable, intent(out) :: at_z(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: thresh
      real(dp), intent(in), optional :: rot(3, 3)
      type(error_t) :: err
      type(trc_decontract_options_t) :: policy
      call read_xyz(xyz, nat, at_z, r)
      if (present(rot)) r = matmul(rot, r)
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "check_mdj_grad: "//err%get_message())
      if (allocated(error)) return
      if (maxval(bas%sh_l) > LMAX) then
         call skip_test(error, "SKIP: basis needs l = "//to_string(maxval(bas%sh_l))// &
                        ", this build has LMAX = "//to_string(LMAX))
         call bas%release()
         return
      end if
      call bas%to_device()
      policy%md_coulomb = split
      call eri%build(bas, thresh, decontract=policy)
   end subroutine setup

   subroutine gradient_pair(eri, bas, dmat, k_scale, j_scale, omega, gmd, grys, nwork_md, nwork_rys)
      !! The same call twice, once per engine.
      type(trc_eri_t), intent(inout) :: eri
      type(trc_basis_t), intent(in) :: bas
      real(dp), intent(in) :: dmat(:, :, :), k_scale, j_scale, omega
      real(dp), intent(out) :: gmd(:, :), grys(:, :)
      integer(int64), intent(out) :: nwork_md, nwork_rys
      eri%grad_j_engine = TRC_J_ENGINE_MDJ
      call eri%gradient(bas, dmat, gmd, k_scale=k_scale, j_scale=j_scale, omega=omega)
      nwork_md = eri%nwork
      eri%grad_j_engine = TRC_J_ENGINE_BINS
      call eri%gradient(bas, dmat, grys, k_scale=k_scale, j_scale=j_scale, omega=omega)
      nwork_rys = eri%nwork
   end subroutine gradient_pair

   subroutine j_only(tag, ndens, banded, error)
      !! MD J gradient alone against Rys with kfac = 0, jfac = 1, every molecule.
      character(*), intent(in) :: tag
      integer, intent(in) :: ndens
      logical, intent(in) :: banded
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      logical :: split
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      integer :: nat, im, id
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), dmat(:, :, :), gmd(:, :), grys(:, :)
      integer(int64) :: nmd, nrys
      real(dp) :: dmax

      call basis_args(basis_json, split)
      call check(error, len_trim(basis_json) > 0, "usage: check_mdj_grad <basis.json> [split]")
      if (allocated(error)) return
      do im = 1, nmol_for(basis_json)
         call setup(trim(MOLS(im)), basis_json, split, bas, eri, nat, at_z, r, error, THRESH)
         if (allocated(error)) return
         allocate (dmat(bas%nao, bas%nao, ndens), gmd(3, nat), grys(3, nat))
         do id = 1, ndens
            call synth_density(bas%nao, id, banded, dmat(:, :, id))
         end do
         call gradient_pair(eri, bas, dmat, 0.0_dp, 1.0_dp, 0.0_dp, gmd, grys, nmd, nrys)
         dmax = maxval(abs(gmd - grys))
         print '(a,a,a,a,a,es10.2,a,es10.2)', '  ', tag, ' ', trim(MOLS(im)), ': max |MD - Rys| ', dmax, &
            '  max |g| ', maxval(abs(grys))
         call check(error, nmd == 0 .and. nrys > 0, trim(MOLS(im))//": the MD route did not engage")
         ! The split variant must really run over the split working basis.
         if (.not. allocated(error) .and. split) &
            call check(error, eri%dec%active, trim(MOLS(im))//": the split basis was not split")
         if (.not. allocated(error)) &
            call check(error, dmax <= AGREE*max(1.0_dp, maxval(abs(grys))), tag//" "//trim(MOLS(im))//": max |MD - Rys| "// &
                       to_string(dmax, "es10.2"))
         deallocate (at_z, r, dmat, gmd, grys)
         call eri%release()
         call bas%release()
         if (allocated(error)) return
      end do
   end subroutine j_only

   subroutine test_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call j_only("RHF dense", 1, .false., error)
   end subroutine test_rhf

   subroutine test_banded(error)
      type(error_type), allocatable, intent(out) :: error
      call j_only("RHF banded", 1, .true., error)
   end subroutine test_banded

   subroutine test_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      call j_only("UHF", 2, .false., error)
   end subroutine test_uhf

   subroutine test_screen(error)
      !! The derivative pair distribution's screening bound, isolated: the MD
      !! gradient at a loose threshold with the (P, Q) screen on, against the
      !! same build with the screen off (md%cut = 0; same primitive pairs, so
      !! the primitive prune cancels). What is left is exactly what the screen
      !! dropped, which a correct bound keeps near the cut (thresh*1e-3 per
      !! pair); a bound too small (the energy's I_P in its place, or the
      !! derivative bound scaled down) drops real pairs and fails it. Dense
      !! synthetic density, so no density screen helps.
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      logical :: split
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      integer :: nat, im
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), dmat(:, :, :), g(:, :), g0(:, :)
      real(dp) :: dmax

      call basis_args(basis_json, split)
      do im = 1, nmol_for(basis_json)
         call setup(trim(MOLS(im)), basis_json, split, bas, eri, nat, at_z, r, error, THRESH_LOOSE)
         if (allocated(error)) return
         allocate (dmat(bas%nao, bas%nao, 1), g(3, nat), g0(3, nat))
         call synth_density(bas%nao, 1, .false., dmat(:, :, 1))
         eri%grad_j_engine = TRC_J_ENGINE_MDJ
         call eri%gradient(bas, dmat, g, k_scale=0.0_dp)
         call check(error, eri%md%cut > 0.0_dp, "the MD engine was not built by the gradient call")
         if (.not. allocated(error)) then
            eri%md%cut = 0.0_dp
            call eri%gradient(bas, dmat, g0, k_scale=0.0_dp)
            dmax = maxval(abs(g - g0))
            print '(a,a,a,es10.2,a,es10.2)', '  ', trim(MOLS(im)), ': max |g_MD(screen on) - g_MD(screen off)| ', &
               dmax, '  max |g| ', maxval(abs(g0))
            call check(error, dmax <= SCREENTOL, trim(MOLS(im))//": the (P,Q) screen dropped "// &
                       to_string(dmax, "es10.2")//" of the gradient")
         end if
         deallocate (at_z, r, dmat, g, g0)
         call eri%release()
         call bas%release()
         if (allocated(error)) return
      end do
   end subroutine test_screen

   subroutine test_full(error)
      !! The sum of the MD Coulomb part and the exchange-only Rys part against the
      !! fused kernels: HF, a pure functional (no exchange), a global hybrid, and
      !! a range-separated hybrid (the long-range exchange is its own call with
      !! jfac = 0, as the SCF makes it).
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      logical :: split
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      integer :: nat, im, id, ic, ndens
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), dmat(:, :, :), gmd(:, :), grys(:, :), gm2(:, :), gr2(:, :)
      integer(int64) :: nmd, nrys
      real(dp) :: dmax
      character(len=24), parameter :: NAMES(4) = [character(len=24) :: "HF", "pure", "global hybrid", "RSH (SR + LR calls)"]

      call basis_args(basis_json, split)
      call check(error, len_trim(basis_json) > 0, "usage: check_mdj_grad <basis.json> [split]")
      if (allocated(error)) return
      do im = 1, nmol_for(basis_json)
         call setup(trim(MOLS(im)), basis_json, split, bas, eri, nat, at_z, r, error, THRESH)
         if (allocated(error)) return
         allocate (gmd(3, nat), grys(3, nat), gm2(3, nat), gr2(3, nat))
         do ndens = 1, 2
            allocate (dmat(bas%nao, bas%nao, ndens))
            do id = 1, ndens
               call synth_density(bas%nao, id, .true., dmat(:, :, id))
            end do
            do ic = 1, 4
               select case (ic)
               case (1)
                  call gradient_pair(eri, bas, dmat, 1.0_dp, 1.0_dp, 0.0_dp, gmd, grys, nmd, nrys)
               case (2)
                  call gradient_pair(eri, bas, dmat, 0.0_dp, 1.0_dp, 0.0_dp, gmd, grys, nmd, nrys)
               case (3)
                  call gradient_pair(eri, bas, dmat, 0.25_dp, 1.0_dp, 0.0_dp, gmd, grys, nmd, nrys)
               case (4)
                  ! The short-range part with the Coulomb J and a fraction of K, then
                  ! the long-range exchange alone (jfac = 0, omega > 0).
                  call gradient_pair(eri, bas, dmat, 0.2_dp, 1.0_dp, 0.0_dp, gmd, grys, nmd, nrys)
                  call gradient_pair(eri, bas, dmat, 0.6_dp, 0.0_dp, 0.33_dp, gm2, gr2, nmd, nrys)
                  gmd = gmd + gm2
                  grys = grys + gr2
               end select
               dmax = maxval(abs(gmd - grys))
               print '(a,a,a,i0,a,a,a,es10.2,a,es10.2)', '  ', trim(MOLS(im)), ' ndens ', ndens, ' ', &
                  trim(NAMES(ic)), ': max |MD - fused| ', dmax, '  max |g| ', maxval(abs(grys))
               call check(error, dmax <= AGREE*max(1.0_dp, maxval(abs(grys))), trim(MOLS(im))//" "//trim(NAMES(ic))//" ndens "// &
                          to_string(ndens)//": max |MD - fused| "//to_string(dmax, "es10.2"))
               if (allocated(error)) exit
            end do
            deallocate (dmat)
            if (allocated(error)) exit
         end do
         deallocate (at_z, r, gmd, grys, gm2, gr2)
         call eri%release()
         call bas%release()
         if (allocated(error)) return
      end do
   end subroutine test_full

   subroutine test_rank_shares(error)
      !! What a run over several ranks adds up: each rank's `gradient` is its
      !! share (every nranks-th ket chunk of the MD Coulomb part, its segments
      !! of the exchange plan) and the shares are summed. Here one process
      !! plays every rank in turn, with no communicator, and the sum of the
      !! shares must be the gradient of rank 0 of 1; a chunk or segment taken
      !! twice or by no rank shows at the size of the gradient. HF (both
      !! parts) and a pure functional (the MD part alone), two and three ranks.
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      logical :: split
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      integer :: nat, nr, ir, ic
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), dmat(:, :, :), gall(:, :), gsum(:, :), gshare(:, :)
      real(dp) :: kfac, dmax

      call basis_args(basis_json, split)
      call check(error, len_trim(basis_json) > 0, "usage: check_mdj_grad <basis.json> [split]")
      if (allocated(error)) return
      call setup(trim(MOLS(NMOL)), basis_json, split, bas, eri, nat, at_z, r, error, THRESH)
      if (allocated(error)) return
      allocate (dmat(bas%nao, bas%nao, 1), gall(3, nat), gsum(3, nat), gshare(3, nat))
      call synth_density(bas%nao, 1, .false., dmat(:, :, 1))
      eri%grad_j_engine = TRC_J_ENGINE_MDJ
      do ic = 1, 2
         kfac = merge(1.0_dp, 0.0_dp, ic == 1)
         call eri%gradient(bas, dmat, gall, k_scale=kfac)
         do nr = 2, 3
            gsum = 0.0_dp
            eri%nranks = nr
            do ir = 0, nr - 1
               eri%rank = ir
               call eri%gradient(bas, dmat, gshare, k_scale=kfac)
               gsum = gsum + gshare
            end do
            eri%rank = 0
            eri%nranks = 1
            dmax = maxval(abs(gsum - gall))
            call check(error, dmax <= AGREE*max(1.0_dp, maxval(abs(gall))) .and. maxval(abs(gshare)) > 0.0_dp, &
                       "k_scale "//to_string(kfac, "f3.1")//", "//to_string(nr)//" ranks: max |sum of shares - whole| "// &
                       to_string(dmax, "es10.2"))
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
      end do
      deallocate (at_z, r, dmat, gall, gsum, gshare)
      call eri%release()
      call bas%release()
   end subroutine test_rank_shares

   subroutine test_fallback(error)
      !! The calls the MD route does not take (a range-separated J, a skeleton
      !! gradient is covered by check_sym_grad) run the fused kernels whatever
      !! the engine argument says: same numbers, same work.
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      logical :: split
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), dmat(:, :, :), gmd(:, :), grys(:, :)
      integer(int64) :: nmd, nrys

      call basis_args(basis_json, split)
      call setup("water.xyz", basis_json, split, bas, eri, nat, at_z, r, error, THRESH)
      if (allocated(error)) return
      allocate (dmat(bas%nao, bas%nao, 1), gmd(3, nat), grys(3, nat))
      call synth_density(bas%nao, 1, .true., dmat(:, :, 1))
      call gradient_pair(eri, bas, dmat, 0.5_dp, 1.0_dp, 0.3_dp, gmd, grys, nmd, nrys)
      call check(error, maxval(abs(gmd - grys)) <= 1.0e-12_dp .and. nmd == nrys .and. nrys > 0, &
                 "omega > 0 with a Coulomb J must take the fused kernels")
      deallocate (at_z, r, dmat, gmd, grys)
      call eri%release()
      call bas%release()
   end subroutine test_fallback

   subroutine test_symmetry(error)
      !! Translational invariance on every molecule, and rotation on a basis whose
      !! shells are all s and p: with a density that is a weighted identity on each
      !! shell (it rotates with the Cartesian s and p functions, which d does not
      !! under an identity density), the gradient of the rotated molecule is the
      !! rotated gradient.
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      logical :: split
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      integer :: nat, im, is, i, nl
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), dmat(:, :, :), g(:, :), g0(:, :)
      real(dp) :: inv, rot(3, 3), ca, sa, cb, sb, resid

      call basis_args(basis_json, split)
      do im = 1, nmol_for(basis_json)
         call setup(trim(MOLS(im)), basis_json, split, bas, eri, nat, at_z, r, error, THRESH)
         if (allocated(error)) return
         allocate (dmat(bas%nao, bas%nao, 1), g(3, nat))
         call synth_density(bas%nao, 1, .false., dmat(:, :, 1))
         eri%grad_j_engine = TRC_J_ENGINE_MDJ
      call eri%gradient(bas, dmat, g, k_scale=0.0_dp)
         inv = maxval(abs(sum(g, dim=2)))
         print '(a,a,a,es10.2)', '  ', trim(MOLS(im)), ': translational-invariance residual ', inv
         call check(error, inv <= INVTOL, trim(MOLS(im))//": translational-invariance residual "// &
                    to_string(inv, "es10.2"))
         deallocate (at_z, r, dmat, g)
         call eri%release()
         call bas%release()
         if (allocated(error)) return
      end do

      ! Rotation. The density is diagonal and constant on each shell.
      ca = cos(0.7_dp); sa = sin(0.7_dp); cb = cos(1.1_dp); sb = sin(1.1_dp)
      rot = matmul(reshape([ca, sa, 0.0_dp, -sa, ca, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3]), &
                   reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, cb, sb, 0.0_dp, -sb, cb], [3, 3]))
      do nl = 1, 2
         if (nl == 1) then
            call setup("water.xyz", "basis_sets/sto-3g.json", .false., bas, eri, nat, at_z, r, error, THRESH)
         else
            call setup("water.xyz", "basis_sets/sto-3g.json", .false., bas, eri, nat, at_z, r, error, THRESH, rot)
         end if
         if (allocated(error)) return
         allocate (dmat(bas%nao, bas%nao, 1), g(3, nat))
         dmat = 0.0_dp
         do is = 1, bas%nshell
            do i = 1, (bas%sh_l(is) + 1)*(bas%sh_l(is) + 2)/2
               dmat(bas%sh_ao(is) + i - 1, bas%sh_ao(is) + i - 1, 1) = 1.0_dp + 0.1_dp*real(is, dp)
            end do
         end do
         eri%grad_j_engine = TRC_J_ENGINE_MDJ
      call eri%gradient(bas, dmat, g, k_scale=0.0_dp)
         if (nl == 1) then
            allocate (g0(3, nat))
            g0 = g
         else
            resid = maxval(abs(g - matmul(rot, g0)))
            print '(a,es10.2,a,es10.2)', '  rotation: max |g(R x) - R g(x)| ', resid, '  max |g| ', maxval(abs(g0))
            call check(error, resid <= ROTTOL, "rotated gradient differs from the rotated gradient by "// &
                       to_string(resid, "es10.2"))
         end if
         deallocate (at_z, r, dmat, g)
         call eri%release()
         call bas%release()
         if (allocated(error)) return
      end do
      deallocate (g0)
   end subroutine test_symmetry

   real(dp) function j_energy(basis_json, split, nat, at_z, r, nao, dmat, error)
      !! Tr(D J[D]) at the geometry r, the density held fixed.
      character(*), intent(in) :: basis_json
      logical, intent(in) :: split
      integer, intent(in) :: nat, at_z(nat), nao
      real(dp), intent(in) :: r(3, nat), dmat(nao, nao)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      type(trc_decontract_options_t) :: policy
      real(dp) :: vhf(nao, nao)
      j_energy = 0.0_dp
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "check_mdj_grad: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      policy%md_coulomb = split
      call eri%build(bas, THRESH, decontract=policy)
      call eri%fock(bas, dmat, vhf, opts=trc_fock_options_t(k_scale=0.0_dp, density_screen=.false.))
      j_energy = sum(dmat*vhf)
      call eri%release()
      call bas%release()
   end function j_energy

   subroutine test_fd(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      logical :: split
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      integer :: nat, ia, d, nbad
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), rp(:, :), rm(:, :), dmat(:, :, :), g(:, :)
      real(dp) :: ep, em, fd, worst, rel

      call basis_args(basis_json, split)
      call setup("water.xyz", basis_json, split, bas, eri, nat, at_z, r, error, THRESH)
      if (allocated(error)) return
      allocate (dmat(bas%nao, bas%nao, 1), g(3, nat), rp(3, nat), rm(3, nat))
      call synth_density(bas%nao, 1, .false., dmat(:, :, 1))
      eri%grad_j_engine = TRC_J_ENGINE_MDJ
      call eri%gradient(bas, dmat, g, k_scale=0.0_dp)
      worst = 0.0_dp
      nbad = 0
      do ia = 1, nat
         do d = 1, 3
            rp = r
            rm = r
            rp(d, ia) = rp(d, ia) + H
            rm(d, ia) = rm(d, ia) - H
            ep = j_energy(trim(basis_json), split, nat, at_z, rp, bas%nao, dmat(:, :, 1), error)
            if (allocated(error)) exit
            em = j_energy(trim(basis_json), split, nat, at_z, rm, bas%nao, dmat(:, :, 1), error)
            if (allocated(error)) exit
            fd = (ep - em)/(2.0_dp*H)
            rel = abs(fd - g(d, ia))/max(abs(g(d, ia)), abs(fd), 1.0e-10_dp)
            worst = max(worst, rel)
            if (.not. (rel <= FDTOL)) nbad = nbad + 1
         end do
         if (allocated(error)) exit
      end do
      if (.not. allocated(error)) then
         print '(a,es10.2,a,i0)', '  finite difference of Tr(D J[D]): worst rel ', worst, '  outside tol ', nbad
         call check(error, nbad == 0, "MD J gradient vs finite difference: worst rel "// &
                    to_string(worst, "es10.2"))
      end if
      deallocate (at_z, r, dmat, g, rp, rm)
      call eri%release()
      call bas%release()
   end subroutine test_fd

end module check_mdj_grad_suite

program check_mdj_grad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_mdj_grad_suite, only: collect_mdj_grad_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("check_mdj_grad", collect_mdj_grad_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_mdj_grad
