!
! The perturbative triples correction (T) of closed-shell RI-CCSD(T): the HOST reference implementation.
!
! JCP 94, 442 (1991) (Rendell, Lee, Komornicki), as PySCF's cc/ccsd_t_slow.py writes it, ported verbatim from the
! metalquicha `triples_correction` and its helpers. It is the oracle the device version is checked against, so it
! keeps the source's loop structure on purpose: do not restructure it for the GPU here.
!
! WHAT IT COMPUTES. With the converged (or any fixed) amplitudes t1(i,a), t2(i,j,a,b), the integral blocks of
! `trc_rccsd_ints` and canonical orbital energies, for each virtual triple a >= b >= c (weight 1/6 when all three
! coincide, 1/2 when two do, 1 otherwise, folded into the denominator) and each of the six permutations p of (a,b,c):
!
!     W_p(i,j,k) = sum_f (i a_p|f b_p) t2(k,j,c_p,f)  -  sum_m (i a_p|j m) t2(m,k,b_p,c_p)    ((a_p,b_p,c_p) = perm_p(a,b,c))
!     V_p(i,j,k) = (i a_p|j b_p) t1(k,c_p)                 (the disconnected term; the t2 f_ov term vanishes: canonical)
!     Z_p = [4 X(i,j,k) + X(k,i,j) + X(j,k,i) - 2 X(k,j,i) - 2 X(i,k,j) - 2 X(j,i,k)] / D3,   X = W_p + V_p / 2,
!     D3 = e_i + e_j + e_k - e_a - e_b - e_c                                                          (`triples_r3`)
!     E(T) = 2 sum_{abc} sum_p sum_ijk Z_p(i,j,k) [ sum_q W_WMAP(p,q)(perm_q(i,j,k)) ]
!
! the two-electron integrals being the fitted ones of `ints`. PySCF's `ccsd_t_slow.kernel` is the same sum with
! the loops written out. The paper's Eq. (1) restriction is wrong; the loops here restrict a >= b >= c over all
! (i,j,k), which is the form that sums every distinct virtual triple once.
!
! DEVIATIONS FROM THE SOURCE. The two `!$omp parallel do` repacks (the t2 repack into tt and the ooov repack into
! oop) run serially here, so the first-touch placement of those arrays across sockets, which the source relied on, is
! lost; the values are identical.
!
! THREE INDEX IDENTITIES let the blocks be read straight out of the ones held, so no transposed copy of ovvv
! (the largest block) is made:
!
!     vvov(a,b,i,f) = ovvv(i,a,f,b)     vooo(a,i,j,m) = ooov(j,m,i,a)     vvoo(a,b,i,j) = ovov(i,a,j,b)
!
! COST AND MEMORY. The particle half of W is O(no^3 nv^4) flops (about 2 no^3 nv^4 in all), the hole half O(no^4 nv^3).
! No t3 is ever held. Shared arrays: tt (nv, no^2 nv) the repacked t2 and oop (no^3, nv) the repacked ooov. Per
! thread, with cb the c batch (capped at nv for the allocation): W 6 no^3 cb, Z 6 no^3, X and D3 2 no^3, the gemm
! operands 2 no^3 cb + no nv + nv cb no, all doubles. For no = 25, nv = 130 and cb = 32 that is 85 MB + 16 MB shared
! and about 34 MB per thread.
!
! HOST ONLY (this routine). `rccsd_triples_finish` runs a device version in a device build (`opts%use_device` defaults
! to true there, false in a host build) and this reference otherwise: `opts%device_routine` picks the one driven by
! occupied triples (trc_rccsd_triples_ijk, the default) or the one driven by virtual triples (trc_rccsd_triples_dev).
! The reference reads `ints%ovvv`, `ints%ooov`, `ints%ovov`, `ints%eps_o`,
! `ints%eps_v` and `amps%t1`, `amps%t2` as they are in host memory, and never touches device data: the host copies must
! be current. They are after `trc_rccsd_solve` (`amps%pull`), and after `ints%build`; a caller that rewrote one of them
! on the device must `update self` first. Nothing here is a `do concurrent`, so under the GPU build it runs on the host
! unchanged. The device version needs the opposite: the device copies current.
!
! THREADS. The loop over the pairs (a,b) is an OpenMP loop, each thread with its own W and scratch. The energy of
! each pair is stored in `e_pair` and summed in pair order afterwards, so the result does not depend on the thread
! count or the schedule. The BLAS is the host BLAS through pic_gemm (the host path of trc_linalg_t, which has no
! entry for host operands in the GPU build), called from the threads.
module trc_rccsd_triples
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_blas_interfaces, only: pic_gemm
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_driver, only: trc_rccsd_state_t
   use trc_rccsd_triples_dev, only: trc_rccsd_triples_dev_options_t, rccsd_triples_dev_run
   use trc_rccsd_triples_ijk, only: rccsd_triples_ijk_run
   implicit none
   private

   public :: trc_rccsd_triples_options_t, RCCSD_TRIPLES_C_BATCH
   public :: RCCSD_TRIPLES_DEVICE_ABC, RCCSD_TRIPLES_DEVICE_IJK
   public :: rccsd_triples_reference, rccsd_triples_finish
   public :: RCCSD_TRIPLES_WMAP, rccsd_triples_energy_block

   integer, parameter :: N_TRIPLE_PERMS = 6
   integer, parameter :: RCCSD_TRIPLES_C_BATCH = 32
   !! The default number of c per batched gemm; sized so a thread's W (6 no^3 doubles per c) stays in cache

   integer, parameter :: RCCSD_TRIPLES_WMAP(N_TRIPLE_PERMS, N_TRIPLE_PERMS) = reshape([ &
                                                                        1, 2, 3, 4, 5, 6, &
                                                                        2, 1, 5, 6, 3, 4, &
                                                                        3, 4, 1, 2, 6, 5, &
                                                                        4, 3, 6, 5, 1, 2, &
                                                                        5, 6, 2, 1, 4, 3, &
                                                                        6, 5, 4, 3, 2, 1], &
                                                                        [N_TRIPLE_PERMS, N_TRIPLE_PERMS], order=[2, 1])
      !! Which W pairs with which occupied permutation, for each Z. Row is the
      !! Z's virtual permutation, column the occupied one, entry the W's. The
      !! S3 multiplication table, written out rather than derived: getting the
      !! positions-versus-labels convention backwards transposes half of it
      !! into a still-symmetric, still-plausible, wrong answer.

   integer, parameter :: RCCSD_TRIPLES_DEVICE_ABC = 1
   !! the device routine driven by virtual triples a >= b >= c (trc_rccsd_triples_dev)
   integer, parameter :: RCCSD_TRIPLES_DEVICE_IJK = 2
   !! the device routine driven by occupied triples i >= j >= k (trc_rccsd_triples_ijk)
   integer, parameter :: DEFAULT_DEVICE_ROUTINE = RCCSD_TRIPLES_DEVICE_IJK
   !! the occupied-driven one: GEMMs sized by the virtual space (H10O5 cc-pVTZ, no 20, nv 265, GH200: see that module)

#ifdef TRC_DEVICE
   logical, parameter :: DEFAULT_USE_DEVICE = .true.   !! the GPU build runs the device version by default
#else
   logical, parameter :: DEFAULT_USE_DEVICE = .false.  !! a host build runs the reference by default: it is faster there
#endif

   type :: trc_rccsd_triples_options_t
      !! Plain values only
      integer :: c_batch = RCCSD_TRIPLES_C_BATCH  !! c per batched gemm; any value >= 1 gives the same energy
      logical :: verbose = .false.                !! one line with E(T) and the wall time
      logical :: use_device = DEFAULT_USE_DEVICE  !! `rccsd_triples_finish` runs the device version; false: this host reference
      !! The default follows the build: true under `TRC_DEVICE` (the GPU build), false in a host build, where the device
      !! version has no OpenMP over pairs and its batched gemm is a serial loop (h2o2 0.46 s against 0.17 s)
      integer :: device_routine = DEFAULT_DEVICE_ROUTINE
      !! which device routine `use_device` runs: RCCSD_TRIPLES_DEVICE_IJK (occupied triples) or _ABC (virtual triples)
      type(trc_rccsd_triples_dev_options_t) :: device  !! budget, forced tile size and verbosity of the device version
   end type trc_rccsd_triples_options_t

contains

   subroutine rccsd_triples_reference(ints, amps, e_t, status, opts, seconds, pair_energies)
      !! E(T) in Hartree from the host copies of `ints` and `amps`; `status` says why not when it is refused
      !!
      !! Refused, with `e_t` = 0: empty integrals, amplitudes missing or of another size than the integrals, no
      !! occupied or no virtual orbitals (an empty space has no triples; that is not a zero), a batch below 1, and a
      !! result that is NaN or infinite (a zero denominator, a NaN amplitude).
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      real(dp), intent(out) :: e_t
      type(error_t), intent(out) :: status
      type(trc_rccsd_triples_options_t), intent(in), optional :: opts
      real(dp), intent(out), optional :: seconds  !! wall time of the triples
      real(dp), intent(out), allocatable, optional :: pair_energies(:)
      !! the energy of each pair (a,b), a >= b, in the order a, b (b fastest), before the factor 2; for the device check

      type(trc_rccsd_triples_options_t) :: o
      type(timer_type) :: clock
      character(len=200) :: line
      character(len=:), allocatable :: why

      e_t = 0.0_dp
      if (present(seconds)) seconds = 0.0_dp
      if (present(opts)) o = opts
      why = refusal(ints, amps, o)
      if (len(why) > 0) then
         call status%set(ERROR_VALIDATION, "trc_rccsd_triples: "//why)
         return
      end if
      call clock%start()
      call triples_correction(ints, amps%t1, amps%t2, ints%no, ints%nv, o%c_batch, e_t, pair_energies)
      if (present(seconds)) seconds = clock%get_elapsed_time()
      if (.not. abs(e_t) <= huge(1.0_dp)) then
         e_t = 0.0_dp
         call status%set(ERROR_VALIDATION, "trc_rccsd_triples: E(T) is NaN or infinite (a zero denominator or a "// &
                         "non-finite amplitude or integral)")
         return
      end if
      if (o%verbose) then
         write (line, "(a,i0,a,i0,a,f20.12,a,f9.2,a)") "  rccsd (T): no ", ints%no, " nv ", ints%nv, "   E(T) ", e_t, &
            "   ", clock%get_elapsed_time(), " s"
         call global_logger%info(trim(line))
      end if
   end subroutine rccsd_triples_reference

   subroutine rccsd_triples_finish(state, status, opts)
      !! (T) for a state `trc_rccsd_solve` left: fills `state%result%e_t`, `e_ccsd_t` and `seconds_triples`
      !!
      !! Refused, with the result untouched, when the run was refused, did not converge or is not finite: (T) of
      !! amplitudes that are not a CCSD solution is not the (T) correction.
      type(trc_rccsd_state_t), intent(inout) :: state
      type(error_t), intent(out) :: status
      type(trc_rccsd_triples_options_t), intent(in), optional :: opts

      real(dp) :: e_t, seconds

      if (state%result%refused) then
         call status%set(ERROR_VALIDATION, "trc_rccsd_triples: the CCSD run was refused; there are no amplitudes for (T)")
         return
      end if
      if (.not. state%result%finite) then
         call status%set(ERROR_VALIDATION, "trc_rccsd_triples: the CCSD run hit a non-finite value; (T) needs its amplitudes")
         return
      end if
      if (.not. state%result%converged) then
         call status%set(ERROR_VALIDATION, "trc_rccsd_triples: the CCSD run did not converge; (T) needs its amplitudes")
         return
      end if
      if (use_device(opts) .and. device_routine(opts) == RCCSD_TRIPLES_DEVICE_IJK) then
         call rccsd_triples_ijk_run(state%ints, state%amps, e_t, status, device_options(opts), seconds)
      else if (use_device(opts)) then
         call rccsd_triples_dev_run(state%ints, state%amps, e_t, status, device_options(opts), seconds)
      else
         call rccsd_triples_reference(state%ints, state%amps, e_t, status, opts, seconds)
      end if
      if (status%has_error()) return
      state%result%e_t = e_t
      state%result%e_ccsd_t = state%result%e_ccsd + e_t
      state%result%seconds_triples = seconds
      state%result%has_triples = .true.
   end subroutine rccsd_triples_finish

   pure function use_device(opts) result(device)
      !! The build's default (device version under `TRC_DEVICE`, else the reference) unless the caller chose
      type(trc_rccsd_triples_options_t), intent(in), optional :: opts
      logical :: device

      device = DEFAULT_USE_DEVICE
      if (present(opts)) device = opts%use_device
   end function use_device

   pure function device_routine(opts) result(routine)
      !! The default device routine unless the caller chose
      type(trc_rccsd_triples_options_t), intent(in), optional :: opts
      integer :: routine

      routine = DEFAULT_DEVICE_ROUTINE
      if (present(opts)) routine = opts%device_routine
   end function device_routine

   function device_options(opts) result(d)
      !! The device options inside `opts`, with its `verbose` joined to the host one; the defaults when absent
      type(trc_rccsd_triples_options_t), intent(in), optional :: opts
      type(trc_rccsd_triples_dev_options_t) :: d

      if (present(opts)) then
         d = opts%device
         d%verbose = opts%device%verbose .or. opts%verbose
      end if
   end function device_options

   function refusal(ints, amps, o) result(why)
      !! Empty when the inputs can be used, else the reason
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_options_t), intent(in) :: o
      character(len=:), allocatable :: why

      why = ""
      if (ints%stamp == 0 .or. .not. allocated(ints%ovvv) .or. .not. allocated(ints%ooov) .or. &
          .not. allocated(ints%ovov) .or. .not. allocated(ints%eps_o) .or. .not. allocated(ints%eps_v)) then
         why = "the integral blocks are empty (no successful build)"
      else if (ints%no < 1 .or. ints%nv < 1) then
         why = "no occupied or no virtual orbitals: there are no triples"
      else if (.not. allocated(amps%t1) .or. .not. allocated(amps%t2)) then
         why = "the amplitudes are not allocated"
      else if (amps%no /= ints%no .or. amps%nv /= ints%nv) then
         why = "the amplitudes are of another size than the integrals"
      else if (any(shape(amps%t1) /= [ints%no, ints%nv]) .or. any(shape(amps%t2) /= [ints%no, ints%no, ints%nv, ints%nv])) then
         why = "the amplitude arrays do not have the shape (no, nv) and (no, no, nv, nv)"
      else if (any(shape(ints%ovvv) /= [ints%no, ints%nv, ints%nv, ints%nv]) .or. &
               any(shape(ints%ooov) /= [ints%no, ints%no, ints%no, ints%nv]) .or. &
               any(shape(ints%ovov) /= [ints%no, ints%nv, ints%no, ints%nv]) .or. &
               size(ints%eps_o) /= ints%no .or. size(ints%eps_v) /= ints%nv) then
         why = "the integral blocks do not have the shapes of (no, nv)"
      else if (o%c_batch < 1) then
         why = "c_batch must be >= 1"
      end if
      if (len(why) > 0) why = trim(why)//" ("//to_string(ints%no)//" occupied, "//to_string(ints%nv)//" virtual)"
   end function refusal

   subroutine triples_term2(oop_u, t2_vw, no, w)
      !! w(i,j,k) -= sum_m (iu|jm) t2(m,k,v,w3), the hole half of W
      !!
      !! One gemm, (n_occ^2 x n_occ)(n_occ x n_occ): `oop_u` is `ooov` for this
      !! u repacked so (i,j) leads, `t2_vw` the (m,k) block of t2 at (v,w3), and
      !! `w` takes the product as the (i,j)-by-k matrix it already is.
      integer, intent(in) :: no
      real(dp), intent(in) :: oop_u(no*no, no)
      real(dp), intent(in) :: t2_vw(no, no)
      real(dp), intent(inout) :: w(no*no, no)

      call pic_gemm(oop_u, t2_vw, w, alpha=-1.0_dp, beta=1.0_dp)
   end subroutine triples_term2

   subroutine triples_v(ints, t1, no, a, b, c, v)
      !! v(i,j,k) = (ia|jb) t1(k,c)
      !!
      !! The reference's second term, `t2(i,j,a,b) f_vo(c,k)`, is absent: the
      !! occupied-virtual Fock block vanishes for canonical orbitals, which this
      !! module assumes throughout.
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :)
      integer, intent(in) :: no, a, b, c
      real(dp), intent(out) :: v(:, :, :)

      integer :: i, j, k

      do k = 1, no
         do j = 1, no
            do i = 1, no
               v(i, j, k) = ints%ovov(i, a, j, b)*t1(k, c)
            end do
         end do
      end do
   end subroutine triples_v

   subroutine triples_r3(x, d3, no, z)
      !! z = [4x(i,j,k) + x(k,i,j) + x(j,k,i) - 2x(k,j,i) - 2x(i,k,j) - 2x(j,i,k)] / d3
      real(dp), intent(in) :: x(:, :, :), d3(:, :, :)
      integer, intent(in) :: no
      real(dp), intent(out) :: z(:, :, :)

      integer :: i, j, k

      do k = 1, no
         do j = 1, no
            do i = 1, no
               z(i, j, k) = (4.0_dp*x(i, j, k) + x(k, i, j) + x(j, k, i) &
                             - 2.0_dp*x(k, j, i) - 2.0_dp*x(i, k, j) &
                             - 2.0_dp*x(j, i, k))/d3(i, j, k)
            end do
         end do
      end do
   end subroutine triples_r3

   subroutine triples_pack_t2(t2, no, nv, tt)
      !! tt(f, (c,j,k)) = t2(k,j,c,f), built once for the whole (T)
      !!
      !! The particle half of W contracts over a virtual `f` against t2 and wants
      !! it with `f` leading. Built whole rather than per triple, so the batched
      !! products below take contiguous slices instead of gathering.
      real(dp), intent(in) :: t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: tt(:, :)

      integer :: c, f, j, k

      do c = 1, nv
         do f = 1, nv
            do j = 1, no
               do k = 1, no
                  tt(f, (c - 1)*no*no + (j - 1)*no + k) = t2(k, j, c, f)
               end do
            end do
         end do
      end do
   end subroutine triples_pack_t2

   subroutine triples_pack_ovvv(ints, no, nv, c0, nc, fixed, fixed_first, mmt)
      !! mmt(f, (c,i)) = (iu|fv) for the four permutations whose pair carries c
      !!
      !! `fixed_first` says which side of the pair the loop's fixed virtual sits
      !! on: true gives (i,fixed|f,c), false gives (i,c|f,fixed). Transposed --
      !! `f` leading rather than `i` -- so that the column range for this batch
      !! is contiguous and the product can take it without a copy.
      type(trc_rccsd_ints_t), intent(in) :: ints
      integer, intent(in) :: no, nv, c0, nc, fixed
      logical, intent(in) :: fixed_first
      real(dp), intent(out) :: mmt(:, :)

      integer :: cl, c, i, f

      do cl = 1, nc
         c = c0 + cl - 1
         do f = 1, nv
            do i = 1, no
               if (fixed_first) then
                  mmt(f, (cl - 1)*no + i) = ints%ovvv(i, fixed, f, c)
               else
                  mmt(f, (cl - 1)*no + i) = ints%ovvv(i, c, f, fixed)
               end if
            end do
         end do
      end do
   end subroutine triples_pack_ovvv

   subroutine triples_correction(ints, t1, t2, no, nv, c_batch, e_triples, pair_energies)
      !! The (T) correction, batched over the innermost virtual
      !!
      !! The outer loop runs a >= b >= c and weights the denominator to recover
      !! the terms it skipped -- six-fold when all three coincide, two-fold when
      !! exactly two do.
      !!
      !! The c loop is not innermost. The particle half of W dominates the
      !! triples at O(n_occ^3 n_vir^4), but per (a,b,c) it is only n_occ n_vir by
      !! n_occ^2 -- too small a product to be worth a call, and called n_vir^3
      !! times. Batching the whole c range of a given (a,b) into one product
      !! fixes the shape without changing the arithmetic. The six permutations
      !! split two ways by where c appears:
      !!
      !!   p = 1,3      c is the third index, so the integral half is fixed and
      !!                the amplitude half stacks   -> (n_occ x n_vir) (n_vir x n_occ^2 nc)
      !!   p = 2,4,5,6  c is in the pair, so the amplitude half is fixed and the
      !!                integral half stacks        -> (nc n_occ x n_vir) (n_vir x n_occ^2)
      !!
      !! Everything downstream -- the denominators, r3, and the thirty-six-term
      !! sum -- is unchanged and still runs one (a,b,c) at a time.
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :)
      integer, intent(in) :: no, nv, c_batch
      ! Explicit shape, so a (m,k) block of it can be handed to the hole-half
      ! gemm by its first element.
      real(dp), intent(in) :: t2(no, no, nv, nv)
      real(dp), intent(out) :: e_triples
      real(dp), intent(out), allocatable, optional :: pair_energies(:)

      real(dp), allocatable :: w(:, :, :, :, :), z(:, :, :, :), x(:, :, :), d3(:, :, :)
      real(dp), allocatable :: tt(:, :), mfix(:, :), mmt(:, :), r1(:, :), r2(:, :)
      integer, allocatable :: alist(:), blist(:)
      real(dp), allocatable :: e_pair(:)
      real(dp), allocatable :: oop(:, :, :, :)
         !! (n_occ, n_occ, n_occ, n_vir) ooov(j,m,i,u) stored as (i,j,m,u), so the
         !! hole half of W is a gemm over m with (i,j) as one row index
      integer :: a, b, c, cl, c0, c1, i, j, k, m, p, nc, no2, cb
      integer :: ipair, npair
      integer :: trip(3), perm(3, N_TRIPLE_PERMS)
      real(dp) :: scale

      no2 = no*no
      ! A batch never exceeds nv, so a larger c_batch would only over-allocate (W is 6 no^3 cb per thread).
      cb = min(c_batch, nv)

      ! Only the shared arrays are allocated here. W and the scratch that feeds
      ! it are allocated per thread inside the parallel region below, which is
      ! also where they are sized -- an allocatable named in a `private` clause
      ! arrives unallocated, and allocating it out here would be both wrong and
      ! a per-thread copy of the wrong size.
      allocate (tt(nv, no2*nv))

      call triples_pack_t2(t2, no, nv, tt)

      allocate (oop(no, no, no, nv))
      do a = 1, nv
         do m = 1, no
            do j = 1, no
               do i = 1, no
                  oop(i, j, m, a) = ints%ooov(j, m, i, a)
               end do
            end do
         end do
      end do

      ! The (a,b) loop is triangular, so it is flattened before being handed to
      ! OpenMP: a triangular nest cannot be collapsed, and the work per outer
      ! index runs from nothing at a = 1 to the whole inner range at a = n_vir,
      ! so splitting on `a` alone would leave the first threads idle. One flat
      ! index over the pairs, taken dynamically, balances itself.
      npair = nv*(nv + 1)/2
      allocate (alist(npair), blist(npair), e_pair(npair))
      ipair = 0
      do a = 1, nv
         do b = 1, a
            ipair = ipair + 1
            alist(ipair) = a
            blist(ipair) = b
         end do
      end do

      e_pair = 0.0_dp

      ! Everything above this line is shared and read-only for the duration:
      ! the packed amplitudes, the permutation table and the pair list. Each
      ! thread keeps its own W and its own scratch, which is what c_batch is
      ! sized against.
      !$omp parallel default(none) &
      !$omp    shared(ints, t1, t2, no, nv, no2, c_batch, cb, tt, alist, blist, npair, oop) &
      !$omp    shared(e_pair) &
      !$omp    private(w, z, x, d3, mfix, mmt, r1, r2, &
      !$omp            a, b, c, cl, c0, c1, nc, i, j, k, p, trip, perm, scale)
      allocate (w(no, no, no, N_TRIPLE_PERMS, cb))
      allocate (z(no, no, no, N_TRIPLE_PERMS))
      allocate (x(no, no, no), d3(no, no, no))
      allocate (mfix(no, nv), mmt(nv, cb*no))
      allocate (r1(no, no2*cb), r2(cb*no, no2))

      !$omp do schedule(dynamic)
      do ipair = 1, npair
         a = alist(ipair)
         b = blist(ipair)

         do c0 = 1, b, c_batch
            c1 = min(c0 + c_batch - 1, b)
            nc = c1 - c0 + 1

            ! ---- p = 1 and 3: the pair is fixed, the amplitudes stack -------
            do k = 1, nv
               do i = 1, no
                  mfix(i, k) = ints%ovvv(i, a, k, b)
               end do
            end do
            call pic_gemm(mfix, tt(:, (c0 - 1)*no2 + 1:c1*no2), r1(:, 1:no2*nc))
            call scatter_w_third(r1, no, nc, 1, w)

            do k = 1, nv
               do i = 1, no
                  mfix(i, k) = ints%ovvv(i, b, k, a)
               end do
            end do
            call pic_gemm(mfix, tt(:, (c0 - 1)*no2 + 1:c1*no2), r1(:, 1:no2*nc))
            call scatter_w_third(r1, no, nc, 3, w)

            ! ---- p = 2, 4, 5, 6: the amplitudes are fixed, the pair stacks --
            ! The third index is b for p = 2 and 5, a for p = 4 and 6, so each
            ! takes the matching contiguous block of `tt`.
            call triples_pack_ovvv(ints, no, nv, c0, nc, a, .true., mmt)
            call pic_gemm(mmt(:, 1:nc*no), tt(:, (b - 1)*no2 + 1:b*no2), &
                          r2(1:nc*no, :), transa="T")
            call scatter_w_pair(r2, no, nc, 2, w)

            call triples_pack_ovvv(ints, no, nv, c0, nc, b, .true., mmt)
            call pic_gemm(mmt(:, 1:nc*no), tt(:, (a - 1)*no2 + 1:a*no2), &
                          r2(1:nc*no, :), transa="T")
            call scatter_w_pair(r2, no, nc, 4, w)

            call triples_pack_ovvv(ints, no, nv, c0, nc, a, .false., mmt)
            call pic_gemm(mmt(:, 1:nc*no), tt(:, (b - 1)*no2 + 1:b*no2), &
                          r2(1:nc*no, :), transa="T")
            call scatter_w_pair(r2, no, nc, 5, w)

            call triples_pack_ovvv(ints, no, nv, c0, nc, b, .false., mmt)
            call pic_gemm(mmt(:, 1:nc*no), tt(:, (a - 1)*no2 + 1:a*no2), &
                          r2(1:nc*no, :), transa="T")
            call scatter_w_pair(r2, no, nc, 6, w)

            do cl = 1, nc
               c = c0 + cl - 1
               perm(:, 1) = [a, b, c]
               perm(:, 2) = [a, c, b]
               perm(:, 3) = [b, a, c]
               perm(:, 4) = [b, c, a]
               perm(:, 5) = [c, a, b]
               perm(:, 6) = [c, b, a]

               ! The hole half, one small gemm per permutation.
               do p = 1, N_TRIPLE_PERMS
                  trip = perm(:, p)
                  call triples_term2(oop(1, 1, 1, trip(1)), t2(1, 1, trip(2), trip(3)), &
                                     no, w(1, 1, 1, p, cl))
               end do

               scale = 1.0_dp
               if (a == c) then
                  scale = 6.0_dp          ! a == b == c
               else if (a == b .or. b == c) then
                  scale = 2.0_dp
               end if

               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        d3(i, j, k) = scale*(ints%eps_o(i) + ints%eps_o(j) + ints%eps_o(k) &
                                             - ints%eps_v(a) - ints%eps_v(b) - ints%eps_v(c))
                     end do
                  end do
               end do

               do p = 1, N_TRIPLE_PERMS
                  trip = perm(:, p)
                  call triples_v(ints, t1, no, trip(1), trip(2), trip(3), x)
                  ! explicit loops, not array syntax: nvfortran -stdpar may offload an array statement
                  do k = 1, no
                     do j = 1, no
                        do i = 1, no
                           x(i, j, k) = w(i, j, k, p, cl) + 0.5_dp*x(i, j, k)
                        end do
                     end do
                  end do
                  call triples_r3(x, d3, no, z(:, :, :, p))
               end do

               ! The thirty-six terms: every Z against every occupied
               ! permutation of the W the table pairs it with.
               e_pair(ipair) = e_pair(ipair) + rccsd_triples_energy_block(no, RCCSD_TRIPLES_WMAP, w(1, 1, 1, 1, cl), z)
            end do
         end do
      end do
      !$omp end do

      deallocate (w, z, x, d3, mfix, mmt, r1, r2)
      !$omp end parallel

      ! Summed here, in pair order, rather than by an OpenMP reduction. A
      ! reduction combines the per-thread partials in whatever order the threads
      ! finish, and dynamic scheduling changes that order between runs; a fixed
      ! array summed in a fixed order makes the correction independent of the
      ! thread count. Worth a n_vir^2/2 array, because this energy is differenced
      ! against others in a many-body expansion, where a wandering last digit is
      ! a term and not noise.
      e_triples = 0.0_dp
      do ipair = 1, npair
         e_triples = e_triples + e_pair(ipair)
      end do
      e_triples = 2.0_dp*e_triples

      if (present(pair_energies)) then
         if (allocated(pair_energies)) deallocate (pair_energies)
         allocate (pair_energies(npair))
         pair_energies = e_pair
      end if
      deallocate (tt, alist, blist, e_pair, oop)
   end subroutine triples_correction

   pure function rccsd_triples_energy_block(no, wmap, w, z) result(acc)
      !! sum_pq sum_ijk W_wmap(p,q)(perm_q(i,j,k)) Z_p(i,j,k), one triple's share
      !!
      !! The thirty-six-term sum, taken one (i,j,k) at a time: the six Z there
      !! are read once and each meets the six occupied permutations of its
      !! paired W, addressed directly rather than through a table of permuted
      !! offsets. The occupied permutations run ijk, ikj, jik, jki, kij, kji --
      !! the order the virtual ones take in `perm` -- which is the column order
      !! of `wmap`.
      integer, intent(in) :: no
      integer, intent(in) :: wmap(N_TRIPLE_PERMS, N_TRIPLE_PERMS)
      real(dp), intent(in) :: w(no, no, no, N_TRIPLE_PERMS)
      real(dp), intent(in) :: z(no, no, no, N_TRIPLE_PERMS)
      real(dp) :: acc

      integer :: i, j, k, p
      real(dp) :: s

      acc = 0.0_dp
      do k = 1, no
         do j = 1, no
            do i = 1, no
               do p = 1, N_TRIPLE_PERMS
                  s = w(i, j, k, wmap(p, 1)) + w(i, k, j, wmap(p, 2)) &
                      + w(j, i, k, wmap(p, 3)) + w(j, k, i, wmap(p, 4)) &
                      + w(k, i, j, wmap(p, 5)) + w(k, j, i, wmap(p, 6))
                  acc = acc + z(i, j, k, p)*s
               end do
            end do
         end do
      end do
   end function rccsd_triples_energy_block

   subroutine scatter_w_third(r, no, nc, p, w)
      !! R(i, (c,j,k)) -> w(i,j,k,p,c), for the permutations where c is third
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: no, nc, p
      real(dp), intent(inout) :: w(:, :, :, :, :)

      integer :: c, i, j, k

      do c = 1, nc
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  w(i, j, k, p, c) = r(i, (c - 1)*no*no + (j - 1)*no + k)
               end do
            end do
         end do
      end do
   end subroutine scatter_w_third

   subroutine scatter_w_pair(r, no, nc, p, w)
      !! R((c,i), (j,k)) -> w(i,j,k,p,c), for the permutations where c is paired
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: no, nc, p
      real(dp), intent(inout) :: w(:, :, :, :, :)

      integer :: c, i, j, k

      do c = 1, nc
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  w(i, j, k, p, c) = r((c - 1)*no + i, (j - 1)*no + k)
               end do
            end do
         end do
      end do
   end subroutine scatter_w_pair

end module trc_rccsd_triples
