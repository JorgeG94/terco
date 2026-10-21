! The papers behind what a run actually did, cited as it does it.
!
! A feature calls `trc_cite(key)` when it runs (the Rys quadrature when a
! Fock build takes a Rys class, the Becke partition when a grid is built,
! ADIIS when it extrapolates, ...), and a functional's own papers come in
! through `trc_cite_functional`. Each key is recorded once, in the order it
! was first cited. `trc_citations_report` logs the list at the end of a
! run, and `trc_citations_text` hands the same lines to the C API, so a
! paper is printed only when the run used the method it describes.
!
! Method references live in `method_ref` below; functional references are
! libxc's own, generated into trc_xc_citations by tools/gen_xc_citations.py.
! A key in neither is a programming error and stops the run, so a typo in a
! `trc_cite` call cannot silently drop a reference.
!
! The list is per process and survives across runs in one process (the C
! API's contexts share it): `trc_citations_reset` empties it. Host only:
! nothing here is called from device code or inside an OpenMP region.
module trc_citations
   use pic_logger, only: global_logger
   use trc_xc_citations, only: xc_functional_ref_keys, xc_ref_text, XC_MAX_REFS, XC_REF_KEY_LEN
   implicit none
   private

   public :: trc_cite, trc_cite_functional, trc_cited, trc_citations_reset
   public :: trc_citations_report, trc_citations_text

   integer, parameter :: MAX_CITED = 128
   character(len=XC_REF_KEY_LEN), save :: cited(MAX_CITED) = ""
   integer, save :: ncited = 0

contains

   subroutine trc_cite(key)
      !! Record one reference by key, once. Unknown keys stop the run.
      character(len=*), intent(in) :: key
      character(len=:), allocatable :: text, doi
      logical :: found
      if (trc_cited(key)) return
      call reference(key, text, doi, found)
      if (.not. found) error stop "trc_cite: no reference with key '"//key//"'"
      if (ncited == MAX_CITED) error stop "trc_cite: more than MAX_CITED references"
      ncited = ncited + 1
      cited(ncited) = key
   end subroutine trc_cite

   subroutine trc_cite_functional(name)
      !! Record every paper a functional cites, by its registry name
      !! (`trc_xc_functional_t%name`), plus libxc's own, which its kernels'
      !! definitions come from.
      character(len=*), intent(in) :: name
      character(len=XC_REF_KEY_LEN) :: keys(XC_MAX_REFS)
      integer :: nkeys, i
      call xc_functional_ref_keys(trim(name), keys, nkeys)
      do i = 1, nkeys
         call trc_cite(trim(keys(i)))
      end do
      call trc_cite("Lehtola2018_1")
   end subroutine trc_cite_functional

   function trc_cited(key) result(is_cited)
      !! Whether `key` has been cited since the last reset.
      character(len=*), intent(in) :: key
      logical :: is_cited
      is_cited = any(cited(1:ncited) == key)
   end function trc_cited

   subroutine trc_citations_reset()
      cited = ""
      ncited = 0
   end subroutine trc_citations_reset

   function trc_citations_text() result(lines)
      !! The cited references, one per line, numbered in the order first
      !! cited, each line ending in a newline; "" when nothing was cited.
      character(len=:), allocatable :: lines
      character(len=:), allocatable :: text, doi
      character(len=8) :: num
      logical :: found
      integer :: i
      lines = ""
      do i = 1, ncited
         call reference(trim(cited(i)), text, doi, found)
         write (num, "(a,i2,a)") "[", i, "] "
         lines = lines//num//text
         if (len(doi) > 0) lines = lines//", doi:"//doi
         lines = lines//new_line("a")
      end do
   end function trc_citations_text

   subroutine trc_citations_report()
      !! Log the cited references through `global_logger%info`.
      character(len=:), allocatable :: lines
      integer :: i0, i1
      if (ncited == 0) return
      call global_logger%info("")
      call global_logger%info("   references for the methods this run used:")
      lines = trc_citations_text()
      i0 = 1
      do while (i0 <= len(lines))
         i1 = index(lines(i0:), new_line("a")) + i0 - 1
         call global_logger%info("   "//lines(i0:i1 - 1))
         i0 = i1 + 1
      end do
   end subroutine trc_citations_report

   subroutine reference(key, text, doi, found)
      character(len=*), intent(in) :: key
      character(len=:), allocatable, intent(out) :: text, doi
      logical, intent(out) :: found
      call method_ref(key, text, doi, found)
      if (.not. found) call xc_ref_text(key, text, doi, found)
   end subroutine reference

   subroutine method_ref(key, text, doi, found)
      !! The papers behind terco's own methods, as the code that implements
      !! each one cites it. A DOI is given only where the README's own
      !! reference list gives one.
      character(len=*), intent(in) :: key
      character(len=:), allocatable, intent(out) :: text, doi
      logical, intent(out) :: found
      found = .true.
      doi = ""
      select case (key)
      case ("GalvezVallejo2022_e2112987")  ! binning, Boys interpolation, VRR, G = 2J - K
         text = "J. L. Galvez Vallejo, G. M. J. Barca, and M. S. Gordon, Mol. Phys. 121, " // &
                "e2112987 (2022)"
         doi = "10.1080/00268976.2022.2112987"
      case ("Schlegel1995_83")  ! Cartesian -> spherical coefficients (scripts/gen_sph.py)
         text = "H. B. Schlegel and M. J. Frisch, Int. J. Quantum Chem. 54, 83 (1995)"
      case ("Becke1988_2547")  ! atomic partition of the molecular grid
         text = "A. D. Becke, J. Chem. Phys. 88, 2547 (1988)"
      case ("Stratmann1996_213")  ! the screenable partition, PARTITION_STRATMANN
         text = "R. E. Stratmann, G. E. Scuseria, and M. J. Frisch, Chem. Phys. Lett. 257, 213 (1996)"
      case ("Treutler1995_346")  ! M4 radial grid and the size adjustment
         text = "O. Treutler and R. Ahlrichs, J. Chem. Phys. 102, 346 (1995)"
      case ("Lebedev1999_477")  ! angular grids, DFT and PCM
         text = "V. I. Lebedev and D. N. Laikov, Dokl. Math. 59, 477 (1999)"
      case ("Stocks2025_10263")  ! batched XC integration (Batch Dense D)
         text = "R. Stocks and G. M. J. Barca, J. Chem. Theory Comput. 21, 10263 (2025)"
      case ("Hu2010_054109")  ! ADIIS
         text = "X. Hu and W. Yang, J. Chem. Phys. 132, 054109 (2010)"
      case ("Lehtola2019_1593")  ! SAP guess
         text = "S. Lehtola, J. Chem. Theory Comput. 15, 1593 (2019)"
      case ("Lange2010_244111")  ! SWIG discretisation of the PCM cavity
         text = "A. W. Lange and J. M. Herbert, J. Chem. Phys. 133, 244111 (2010)"
      case ("Bondi1964_441")  ! PCM radii
         text = "A. Bondi, J. Phys. Chem. 68, 441 (1964)"
      case ("Mantina2009_5806")  ! PCM radii Bondi does not cover
         text = "M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer, and D. G. Truhlar, " // &
                "J. Phys. Chem. A 113, 5806 (2009)"
      case ("Pople1979_225")  ! CPHF
         text = "J. A. Pople, R. Krishnan, H. B. Schlegel, and J. S. Binkley, " // &
                "Int. J. Quantum Chem. 16, 225 (1979)"
      case ("Stocks2024_7503")  ! RI-MP2 multi-GPU pair distribution
         text = "R. Stocks, E. Palethorpe, and G. M. J. Barca, J. Chem. Theory Comput. 20, 7503 (2024)"
      case default
         found = .false.
         text = ""
      end select
   end subroutine method_ref

end module trc_citations
