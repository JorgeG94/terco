!
! The analytic gradient over MPI ranks: the frame every rank enters
! (trc_run_gradient, trc_grad_par, trc_comm_agree).
!
! Four things are asserted, none of which needs a GPU count or a timing:
!
!   FRAME     grad_par_finish with a partial gradient per rank: the sum lands
!             on every rank with the same bits; a failure injected on one rank,
!             on rank 0, on all but rank 0, or on every rank returns on EVERY
!             rank with the status and message of the lowest-numbered failing
!             rank, leaves the gradient alone, and does not hang.
!   CHEMISTRY water/6-31G*, split SCF then trc_run_gradient on all ranks:
!             RHF, B3LYP and C-PCM RHF each match a one-rank context's
!             gradient to 1e-8 and are bitwise equal on every rank. (PBE is
!             scf_mpi's.) Every term is rank 0's today, so these also pin that
!             the other ranks add nothing: a rank that added its share of a
!             term would show up as a doubled term.
!   REFUSAL   a gradient asked for before any SCF over the ranks returns the
!             same non-OK status and message on every rank. That refusal, like
!             every refusal decided from the replicated settings, is taken
!             before the first collective call; no refusal reachable through
!             the public API is decided on one rank only (a failing PCM build
!             would be, but the SCF fails first on the same settings), so the
!             failure agreement itself is tested under FRAME.
!
! What guards what, honestly: the FRAME reduction check is the only thing that
! guards the SUM (every term is rank 0's today, so the chemistry cases would
! still pass with the allreduce removed); the CHEMISTRY spread-over-ranks check
! guards that the finish (agreement + broadcast) is reached on every rank and
! leaves identical bits, and the one-rank comparison guards that no rank adds
! a term twice; the REFUSAL case only pins "no hang on a pre-communicator
! return". The failure agreement is tested by FRAME injection alone.
! Cases beyond RHF/B3LYP/C-PCM cover the branches the frame restructured:
! --df (the DF term's enter/exit data on ranks that skip it), UHF (two spin
! densities) and wB97X-V (the VV10 block).
!
! A plain program, not test-drive: it makes collective MPI calls (scf_mpi's
! reason). Its helpers live in a module for the same reason: a procedure that
! makes a collective call must have no internal procedures (NVHPC 26.5 +
! HPC-X).
!
module grad_mpi_support
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr, c_char, c_null_char
   use pic_types, only: dp, int32
   use pic_mpi_lib, only: comm_t, allreduce, MPI_MAX, MPI_MIN
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_method, &
                       trc_set_df, trc_set_jkfit_json, trc_set_force_cartesian, trc_set_guess, trc_set_convergence, trc_set_verbose, trc_set_comm, trc_set_pcm, &
                       trc_run_scf, trc_run_gradient, trc_gradient, trc_message, TRC_OK, TRC_GUESS_GWH
   use trc_test_basis, only: read_xyz
   use trc_grad_par, only: trc_grad_par_t, grad_par_finish
   implicit none
   private
   public :: frame_checks, gradient_case, gradient_spread, refusal_case, case_t

   type :: case_t
      character(len=12) :: tag = ""
      character(len=12) :: func = ""
      character(len=8) :: pcm = ""
      integer :: charge = 0
      integer :: mult = 1
      logical :: df = .false.
   end type case_t

contains

   function cstr(s) result(c)
      character(*), intent(in) :: s
      character(kind=c_char) :: c(len(s) + 1)
      integer :: i
      do i = 1, len(s)
         c(i) = s(i:i)
      end do
      c(len(s) + 1) = c_null_char
   end function cstr

   function context_message(h) result(m)
      type(c_ptr), intent(in) :: h
      character(len=:), allocatable :: m
      character(kind=c_char) :: buf(512)
      integer(c_int) :: rc
      integer :: i, n
      buf = c_null_char
      rc = trc_message(h, buf, 512_c_int)
      n = 0
      do i = 1, 512
         if (buf(i) == c_null_char) exit
         n = i
      end do
      m = repeat(" ", n)
      do i = 1, n
         m(i:i) = buf(i)
      end do
   end function context_message

   subroutine new_context(h, fcomm, cs)
      !! water/6-31G* for the case `cs` (functional, continuum, charge and
      !! multiplicity, --df), over `fcomm` (-1 = one rank).
      type(c_ptr), intent(out) :: h
      integer(c_int), intent(in) :: fcomm
      type(case_t), intent(in) :: cs
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      integer(c_int) :: rc
      call read_xyz('water.xyz', nat, at_z, r)
      rc = trc_create(h)
      rc = trc_set_verbose(h, 0_c_int)
      rc = trc_set_molecule(h, int(nat, c_int), real(at_z, c_double), r, int(cs%charge, c_int), int(cs%mult, c_int))
      rc = trc_set_basis_json(h, cstr('basis_sets/6-31gs.json'))
      if (cs%df) then
         rc = trc_set_force_cartesian(h, 1_c_int)
         rc = trc_set_jkfit_json(h, cstr('basis_sets/def2-universal-jkfit.json'))
         rc = trc_set_df(h, 1_c_int)
      end if
      rc = trc_set_method(h, cstr(trim(cs%func)), 3_c_int)
      rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
      rc = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
      if (len_trim(cs%pcm) > 0) rc = trc_set_pcm(h, cstr(trim(cs%pcm)), 78.3553_c_double, 0_c_int, 0.0_c_double)
      if (fcomm /= -1_c_int) rc = trc_set_comm(h, fcomm)
   end subroutine new_context

   subroutine gradient_case(fcomm, cs, g, ok)
      !! SCF then gradient of one context; `g` is the (3, 3) gradient, `ok`
      !! false when any step returned a non-OK status.
      integer(c_int), intent(in) :: fcomm
      type(case_t), intent(in) :: cs
      real(dp), intent(out) :: g(9)
      logical, intent(inout) :: ok
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: gc(9)
      g = 0.0_dp
      call new_context(h, fcomm, cs)
      rc = trc_run_scf(h)
      if (rc == TRC_OK) rc = trc_run_gradient(h)
      if (rc == TRC_OK) rc = trc_gradient(h, gc)
      if (rc /= TRC_OK) then
         ok = .false.
      else
         g = gc
      end if
      rc = trc_destroy(h)
   end subroutine gradient_case

   function gradient_spread(comm, g) result(spread)
      !! Largest max-min over the ranks of any component of `g`.
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: g(:)
      real(dp) :: spread
      real(dp) :: gmax, gmin
      integer :: k
      spread = 0.0_dp
      do k = 1, size(g)
         gmax = g(k); gmin = g(k)
         if (comm%size() > 1) then
            call allreduce(comm, gmax, op=MPI_MAX)
            call allreduce(comm, gmin, op=MPI_MIN)
         end if
         spread = max(spread, gmax - gmin)
      end do
   end function gradient_spread

   subroutine refusal_case(comm, fcomm, ok)
      !! trc_run_gradient before any SCF, over the ranks: every rank must get
      !! the same non-OK status and the same text.
      type(comm_t), intent(in) :: comm
      integer(c_int), intent(in) :: fcomm
      logical, intent(inout) :: ok
      type(c_ptr) :: h
      integer(c_int) :: rc
      character(len=:), allocatable :: msg
      integer(int32) :: nbad
      call new_context(h, fcomm, case_t())
      rc = trc_run_gradient(h)
      msg = context_message(h)
      nbad = 0
      if (rc == TRC_OK .or. len(msg) == 0) nbad = 1
      ! The text is compared by its hash-like checksum and length across ranks.
      call same_everywhere(comm, int(rc, int32), nbad)
      call same_everywhere(comm, int(len(msg), int32), nbad)
      call same_everywhere(comm, int(checksum(msg), int32), nbad)
      if (comm%rank() == 0) print '(a,i0,a,a)', "grad_mpi: refusal before the SCF, status ", rc, "  -- ", msg
      if (nbad /= 0) ok = .false.
      rc = trc_destroy(h)
   end subroutine refusal_case

   integer function checksum(s)
      character(*), intent(in) :: s
      integer :: i
      checksum = 0
      do i = 1, len(s)
         checksum = mod(checksum*31 + iachar(s(i:i)), 1000003)
      end do
   end function checksum

   subroutine same_everywhere(comm, v, nbad)
      !! Raise `nbad` when `v` is not the same on every rank.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: v
      integer(int32), intent(inout) :: nbad
      integer(int32) :: vmax, vmin
      vmax = v; vmin = v
      if (comm%size() > 1) then
         call allreduce(comm, vmax, op=MPI_MAX)
         call allreduce(comm, vmin, op=MPI_MIN)
      end if
      if (vmax /= vmin) nbad = nbad + 1
   end subroutine same_everywhere

   subroutine frame_checks(comm, ok)
      !! grad_par_finish: the reduction and the failure agreement.
      type(comm_t), intent(in) :: comm
      logical, intent(inout) :: ok
      type(trc_grad_par_t) :: par
      real(dp) :: grad(3, 2), before(3, 2), expect(3, 2)
      integer(int32) :: st, want_st
      character(len=:), allocatable :: msg, want_msg
      integer :: nr, r, k, pat, lowest
      logical :: fails, good
      nr = comm%size()
      r = comm%rank()
      par%comm = comm
      par%has_comm = .true.
      par%rank = r
      par%nranks = nr

      ! --- the reduction: rank-dependent partials, the sum on every rank
      do k = 1, 6
         before(mod(k - 1, 3) + 1, (k - 1)/3 + 1) = real(r + 1, dp)*0.1_dp + 0.25_dp*real(k, dp)
         expect(mod(k - 1, 3) + 1, (k - 1)/3 + 1) = real(nr*(nr + 1)/2, dp)*0.1_dp + real(nr, dp)*0.25_dp*real(k, dp)
      end do
      grad = before
      st = 0
      msg = ""
      call grad_par_finish(par, grad, st, msg)
      good = st == 0 .and. maxval(abs(grad - expect)) < 1.0e-13_dp .and. gradient_spread(comm, reshape(grad, [6])) == 0.0_dp
      if (comm%rank() == 0) print '(a,l1,a,es9.2)', "grad_mpi: frame reduction ok ", good, "  max |sum - expected| ", &
         maxval(abs(grad - expect))
      if (.not. good) ok = .false.

      ! --- the agreement: pattern 1 the last rank, 2 rank 0, 3 every rank but 0, 4 every rank
      do pat = 1, 4
         if (pat == 3 .and. nr < 2) cycle   ! no rank but 0 to fail
         select case (pat)
         case (1); fails = r == nr - 1; lowest = nr - 1
         case (2); fails = r == 0; lowest = 0
         case (3); fails = r >= 1; lowest = min(1, nr - 1)
         case default; fails = .true.; lowest = 0
         end select
         grad = before
         st = 0
         msg = ""
         if (fails) then
            st = int(10 + r, int32)
            msg = "injected failure on rank "//achar(iachar('0') + r)
         end if
         call grad_par_finish(par, grad, st, msg)
         want_st = int(10 + lowest, int32)
         want_msg = "injected failure on rank "//achar(iachar('0') + lowest)
         good = st == want_st .and. msg == want_msg .and. all(grad == before)
         if (.not. good) print '(a,i0,a,i0,a,i0,a,a)', "grad_mpi: FRAME FAIL rank ", r, " pattern ", pat, &
            " status ", st, " message ", msg
         if (r == 0) print '(a,i0,a,l1,a,i0,a,a)', "grad_mpi: frame failure pattern ", pat, " ok ", good, &
            "  agreed status ", st, "  ", msg
         if (.not. good) ok = .false.
      end do
   end subroutine frame_checks

end module grad_mpi_support

program grad_mpi
   use, intrinsic :: iso_c_binding, only: c_int
   use pic_types, only: dp
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize
   use trc_api, only: trc_bind_device
   use grad_mpi_support, only: frame_checks, gradient_case, gradient_spread, refusal_case, case_t
   implicit none

   type(comm_t) :: comm
   integer :: dev, ic
   integer(c_int) :: fcomm
   logical :: ok
   real(dp) :: g_split(9), g_one(9), gspread, gdiff
   integer, parameter :: NCASE = 6
   type(case_t) :: cases(NCASE)
   call pic_mpi_init()
   comm = comm_world()
   dev = trc_bind_device(comm%rank())
   fcomm = int(comm%fortran_handle(), c_int)
   ok = .true.

   call frame_checks(comm, ok)

   cases(1) = case_t("rhf")
   cases(2) = case_t("b3lyp", "b3lyp")
   cases(3) = case_t("rhf + cpcm", "", "cpcm")
   cases(4) = case_t("rhf --df", "", "", 0, 1, .true.)
   cases(5) = case_t("uhf h2o+", "", "", 1, 2)
   cases(6) = case_t("wb97x-v", "wb97x-v")
   do ic = 1, NCASE
      call gradient_case(fcomm, cases(ic), g_split, ok)
      call gradient_case(-1_c_int, cases(ic), g_one, ok)
      gspread = gradient_spread(comm, g_split)
      gdiff = maxval(abs(g_split - g_one))
      if (comm%rank() == 0) print '(a,i0,a,a,a,es9.2,a,es9.2)', "grad_mpi: ", comm%size(), " rank(s)  ", cases(ic)%tag, &
         "  gradient spread over ranks ", gspread, "  vs one rank ", gdiff
      if (gspread /= 0.0_dp .or. gdiff > 1.0e-8_dp) ok = .false.
   end do

   call refusal_case(comm, fcomm, ok)

   if (comm%rank() == 0) print '(a)', merge("grad_mpi: PASS", "grad_mpi: FAIL", ok)
   call pic_mpi_finalize()
   if (.not. ok) stop 1
end program grad_mpi
