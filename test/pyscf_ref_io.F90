!
! A tiny reader for the committed pyscf references under test/pyscf_refs/.
!
! Jorge's rule, from the scf_rks_vs_pyscf/UKS red build: pyscf never runs in
! CI. test/pyscf_refs/generate.py runs it once, offline, on a developer
! machine with pyscf installed, and writes one plain-text file per case --
! the same "recorded oracle" shape as reference.dat/reference1e.dat/
! reference_df.dat/reference_xc.dat, just with pyscf as the oracle instead of
! libfint or libxc. Every Fortran check below reads its own file back and
! compares; none of them import anything, because none of them are Python.
!
! FORMAT. Comment lines start with '#' and may appear anywhere; everything
! else is TAG value(s), one tag per logical record:
!   INT   <tag> <i>
!   REAL  <tag> <r>                (%.17e from the generator; exact round trip)
!   STR   <tag> <string, no spaces>
!   VECTOR <tag> <n>
!   <n reals, list-directed, any layout>
!   MATRIX <tag> <nr> <nc>
!   <nr*nc reals, column-major (Fortran order), list-directed>
!   MATRIX3 <tag> <n1> <n2> <n3>
!   <n1*n2*n3 reals, column-major, list-directed>
! A reader asks for a tag by name; a mismatch is a hard stop (the file and
! the check have drifted, which is exactly what should not pass silently).
!
! Every matrix/vector here is already in THIS BUILD's own AO order and
! normalisation convention -- generate.py does the permutation and the
! per-AO rescaling once, offline, the same transforms rys_grad_ref.py and
! its siblings used to do at test time. A Fortran check compares its own
! array to what it reads here directly, no permutation needed.
!
module pyscf_ref_io
   use pic_types, only: dp
   implicit none
   private

   public :: ref_open, ref_require, ref_close
   public :: ref_read_int, ref_read_real, ref_read_str
   public :: ref_read_vector, ref_read_matrix, ref_read_matrix3
   public :: ref_report, REF_DIR

   character(len=*), parameter :: REF_DIR = 'pyscf_refs/'

contains

   !
   ! Open a reference file; `ok` is false (not a crash) if it is missing, so
   ! a caller still inside its own LMAX self-skip can decide what that means.
   !
   subroutine ref_open(name, u, ok)
      character(len=*), intent(in) :: name
      integer, intent(out) :: u
      logical, intent(out) :: ok
      integer :: ios
      open (newunit=u, file=REF_DIR//trim(name), status='old', action='read', iostat=ios)
      ok = (ios == 0)
   end subroutine ref_open

   !
   ! Open or FAIL LOUDLY. For a case already known to be in scope (past any
   ! LMAX/DF_LMAX guard) -- a missing file here means "run generate.py", not
   ! "nothing to test yet", so this stops rather than skips.
   !
   subroutine ref_require(name, u)
      character(len=*), intent(in) :: name
      integer, intent(out) :: u
      logical :: ok
      call ref_open(name, u, ok)
      if (.not. ok) then
         print '(a)', 'FAIL: missing reference '//REF_DIR//trim(name)// &
            ' -- run test/pyscf_refs/generate.py (see test/pyscf_refs/README.md)'
         stop 1
      end if
   end subroutine ref_require

   subroutine ref_close(u)
      integer, intent(in) :: u
      close (u)
   end subroutine ref_close

   subroutine next_line(u, line)
      integer, intent(in) :: u
      character(len=*), intent(out) :: line
      do
         read (u, '(a)') line
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         exit
      end do
   end subroutine next_line

   subroutine expect_tag(line, tag, rest)
      character(len=*), intent(in) :: line, tag
      character(len=*), intent(out) :: rest
      integer :: n
      n = len_trim(tag)
      if (line(1:n) /= tag(1:n)) then
         print '(a)', 'pyscf_ref_io: expected tag "'//trim(tag)//'", got: '//trim(line)
         stop 1
      end if
      rest = adjustl(line(n + 1:))
   end subroutine expect_tag

   subroutine ref_read_int(u, tag, ival)
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      integer, intent(out) :: ival
      character(len=512) :: line, rest
      call next_line(u, line)
      call expect_tag(line, tag, rest)
      read (rest, *) ival
   end subroutine ref_read_int

   subroutine ref_read_real(u, tag, rval)
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      real(dp), intent(out) :: rval
      character(len=512) :: line, rest
      call next_line(u, line)
      call expect_tag(line, tag, rest)
      read (rest, *) rval
   end subroutine ref_read_real

   subroutine ref_read_str(u, tag, sval)
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      character(len=*), intent(out) :: sval
      character(len=512) :: line, rest
      call next_line(u, line)
      call expect_tag(line, tag, rest)
      sval = trim(rest)
   end subroutine ref_read_str

   subroutine ref_read_vector(u, tag, n, v)
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      integer, intent(in) :: n
      real(dp), intent(out) :: v(n)
      character(len=512) :: line, rest
      integer :: m
      call next_line(u, line)
      call expect_tag(line, tag, rest)
      read (rest, *) m
      if (m /= n) then
         print '(a,a,a,i0,a,i0)', 'pyscf_ref_io: ', trim(tag), ' length mismatch: file has ', m, ', wanted ', n
         stop 1
      end if
      read (u, *) v
   end subroutine ref_read_vector

   subroutine ref_read_matrix(u, tag, nr, nc, a)
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      integer, intent(in) :: nr, nc
      real(dp), intent(out) :: a(nr, nc)
      character(len=512) :: line, rest
      integer :: mr, mc
      call next_line(u, line)
      call expect_tag(line, tag, rest)
      read (rest, *) mr, mc
      if (mr /= nr .or. mc /= nc) then
         print '(a,a,a,i0,1x,i0,a,i0,1x,i0)', 'pyscf_ref_io: ', trim(tag), ' shape mismatch: file has ', mr, mc, &
            ', wanted ', nr, nc
         stop 1
      end if
      read (u, *) a
   end subroutine ref_read_matrix

   subroutine ref_read_matrix3(u, tag, n1, n2, n3, a)
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      integer, intent(in) :: n1, n2, n3
      real(dp), intent(out) :: a(n1, n2, n3)
      character(len=512) :: line, rest
      integer :: m1, m2, m3
      call next_line(u, line)
      call expect_tag(line, tag, rest)
      read (rest, *) m1, m2, m3
      if (m1 /= n1 .or. m2 /= n2 .or. m3 /= n3) then
         print '(a)', 'pyscf_ref_io: '//trim(tag)//' shape mismatch'
         stop 1
      end if
      read (u, *) a
   end subroutine ref_read_matrix3

   !
   ! One line of PASS/FAIL bookkeeping, shared across every check that reads
   ! a reference here -- the same shape check_ao.F90/check_xc_energy.F90's
   ! own `report` already used, lifted out so the pyscf-comparison checks
   ! share one copy of it too.
   !
   subroutine ref_report(nfail, ok, what, value)
      integer, intent(inout) :: nfail
      logical, intent(in) :: ok
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: value
      if (.not. ok) nfail = nfail + 1
      print '(a,1x,a,t64,es10.2)', merge('ok  ', 'FAIL', ok), what, value
   end subroutine ref_report

end module pyscf_ref_io
