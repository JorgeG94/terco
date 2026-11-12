!
! How the analytic gradient is shared over the ranks of a communicator.
!
! Every rank enters the gradient together and accumulates a PARTIAL
! 3 x natom gradient. Each term is either
!
!   SPLIT      its work items are dealt so that each runs on exactly one rank,
!   ROOT-ONLY  rank 0 adds the whole term and every other rank adds nothing.
!
! GRAD_TERM_SPLIT below is the one table that says which; a term is moved from
! root-only to split by setting its entry and teaching its own code to deal
! the items out (par%computes(term) is the guard each term already sits
! behind). Today every entry is false, so rank 0's numbers are exactly the
! single-rank ones.
!
! The partial gradients are summed by grad_par_finish: ONE allreduce, then a
! broadcast of rank 0's copy so that every rank holds the same bits (an
! allreduce need not leave identical bits; test/scf_mpi.F90 asserts they are
! identical, and the eri_t%sum_ranks note says why this matters). Before the
! sum the ranks agree on failure (trc_comm_agree), so a refusal on one rank
! returns on all of them with the same status and message.
!
module trc_grad_par
   use pic_types, only: dp, int32
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   use trc_comm_agree, only: comm_agree_failure
   implicit none
   private
   public :: trc_grad_par_t, grad_par_finish
   public :: GRAD_TERM_FOUR_CENTRE, GRAD_TERM_DF, GRAD_TERM_ONE_ELECTRON, GRAD_TERM_NUC_REP, &
             GRAD_TERM_XC, GRAD_TERM_VV10, GRAD_TERM_PCM, GRAD_NTERMS, GRAD_TERM_SPLIT

   integer, parameter :: GRAD_TERM_FOUR_CENTRE = 1   !! eri%gradient, the full-range and long-range passes
   integer, parameter :: GRAD_TERM_DF = 2            !! the fitted (--df) two-electron term
   integer, parameter :: GRAD_TERM_ONE_ELECTRON = 3  !! overlap Pulay, kinetic, nuclear attraction
   integer, parameter :: GRAD_TERM_NUC_REP = 4       !! nuclear repulsion
   integer, parameter :: GRAD_TERM_XC = 5            !! semilocal exchange-correlation, with weight derivatives
   integer, parameter :: GRAD_TERM_VV10 = 6          !! VV10 nonlocal correlation
   integer, parameter :: GRAD_TERM_PCM = 7           !! the polarisable continuum
   integer, parameter :: GRAD_NTERMS = 7

   logical, parameter :: GRAD_TERM_SPLIT(GRAD_NTERMS) = [ &
                         .false., &   ! four-centre
                         .false., &   ! df
                         .false., &   ! one-electron
                         .false., &   ! nuclear repulsion
                         .false., &   ! xc
                         .false., &   ! vv10
                         .false.]     ! pcm
   !! .true. where a term is dealt over the ranks, .false. where rank 0 alone
   !! computes it.

   type :: trc_grad_par_t
      !! The communicator a gradient runs over and the rank layout, carried as
      !! one value so that no argument list grows. The default is one rank
      !! with no communicator, which computes every term.
      type(comm_t) :: comm                 !! valid only when `has_comm`
      logical :: has_comm = .false.
      integer :: rank = 0
      integer :: nranks = 1
   contains
      procedure :: computes => grad_par_computes
   end type trc_grad_par_t

contains

   logical function grad_par_computes(this, term) result(yes)
      !! Whether this rank computes `term`: always on rank 0, elsewhere only for
      !! a split term.
      class(trc_grad_par_t), intent(in) :: this
      integer, intent(in) :: term
      yes = this%rank == 0 .or. GRAD_TERM_SPLIT(term)
   end function grad_par_computes

   subroutine grad_par_finish(par, grad, status, message)
      !! Collective over `par%comm` when it has one, nothing otherwise. Agrees
      !! on failure (`status` 0 is success, `message` the text of a failure;
      !! both come back as the lowest-numbered failing rank's), and when no rank
      !! failed sums the partial `grad` over the ranks and leaves rank 0's copy
      !! on every rank. After a failure `grad` is not touched.
      type(trc_grad_par_t), intent(in) :: par
      real(dp), intent(inout) :: grad(:, :)
      integer(int32), intent(inout) :: status
      character(len=:), allocatable, intent(inout) :: message
      real(dp), allocatable :: gbuf(:)
      integer :: n
      if (.not. par%has_comm) return
      call comm_agree_failure(par%comm, status, message)
      if (status /= 0_int32) return
      n = size(grad)
      gbuf = reshape(grad, [n])
      call allreduce(par%comm, gbuf, op=MPI_SUM)
      call bcast(par%comm, gbuf, n, 0)
      grad = reshape(gbuf, shape(grad))
   end subroutine grad_par_finish

end module trc_grad_par
