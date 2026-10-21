program check_citations
   !! The reference list a run prints: trc_citations, through the C API
   !! entries trc_citation_list and trc_citation_reset.
   !!
   !! Static: every functional terco runs resolves to at least one paper,
   !! and every key it cites has text. A functional added to the registry
   !! without a regenerated trc_xc_citations fails here.
   !!
   !! Per run, water/6-31G*: an RHF run cites the ERI paper and no grid or
   !! functional paper; a PBE run adds the grid (Becke, Treutler-Ahlrichs,
   !! Lebedev), the batched XC integration, PBE's own paper and libxc's;
   !! a reset empties the list. Each claim is one a wrong cite (or a cite
   !! in the wrong place, run whatever the method) breaks.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char
   use pic_types, only: dp
   use trc_api, only: trc_bind_device
   use trc_test_basis, only: read_xyz
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name, XC_AVAILABLE_FUNCTIONALS
   use trc_xc_citations, only: xc_functional_ref_keys, xc_ref_text, XC_MAX_REFS, XC_REF_KEY_LEN
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_force_cartesian, trc_set_method, trc_set_convergence, &
                               trc_run_scf, trc_citation_list, trc_citation_reset, TRC_OK, TRC_ERR_BADARG
   implicit none

   integer :: natm, dev
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   logical :: ok
   character(len=:), allocatable :: refs

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   ok = .true.

   call check_every_functional()

   call expect(trc_citation_reset() == TRC_OK, "reset")
   call run_water("")
   refs = citation_list()
   call expect(index(refs, "Galvez Vallejo") > 0, "RHF cites the ERI paper")
   call expect(index(refs, "Becke, J. Chem. Phys. 88") == 0, "RHF cites no grid paper")
   call expect(index(refs, "SoftwareX") == 0, "RHF cites no libxc paper")
   call check_buffer_protocol()

   call expect(trc_citation_reset() == TRC_OK, "reset")
   call run_water("pbe")
   refs = citation_list()
   call expect(index(refs, "Galvez Vallejo") > 0, "PBE cites the ERI paper")
   call expect(index(refs, "Becke, J. Chem. Phys. 88") > 0, "PBE cites the Becke partition")
   call expect(index(refs, "Treutler and R. Ahlrichs") > 0, "PBE cites the radial grid")
   call expect(index(refs, "Lebedev and D. N. Laikov") > 0, "PBE cites the angular grid")
   call expect(index(refs, "Stocks and G. M. J. Barca") > 0, "PBE cites the batched XC integration")
   call expect(index(refs, "Phys. Rev. Lett. 77, 3865") > 0, "PBE cites PBE")
   call expect(index(refs, "SoftwareX") > 0, "PBE cites libxc")
   call expect(index(refs, "Vydrov") == 0, "PBE cites no VV10")
   call expect(count_lines(refs) == count_unique_lines(refs), "no reference printed twice")

   call expect(trc_citation_reset() == TRC_OK, "reset")
   call expect(len(citation_list()) == 0, "reset empties the list")

   if (.not. ok) then
      print '(a)', "check_citations: FAIL"
      stop 1
   end if
   print '(a)', "check_citations: PASS"

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "check_citations: FAILED: "//what
         ok = .false.
      end if
   end subroutine expect

   subroutine check_every_functional()
      type(trc_xc_functional_t) :: f
      character(len=XC_REF_KEY_LEN) :: keys(XC_MAX_REFS)
      character(len=:), allocatable :: text, doi, name
      integer :: i0, i1, nkeys, k, nfunc
      logical :: found
      nfunc = 0
      i0 = 1
      do while (i0 <= len_trim(XC_AVAILABLE_FUNCTIONALS))
         i1 = index(XC_AVAILABLE_FUNCTIONALS(i0:), ",") + i0 - 1
         if (i1 < i0) i1 = len_trim(XC_AVAILABLE_FUNCTIONALS) + 1
         name = trim(adjustl(XC_AVAILABLE_FUNCTIONALS(i0:i1 - 1)))
         i0 = i1 + 1
         if (.not. xc_functional_by_name(name, f)) then
            call expect(.false., "listed functional "//name//" does not resolve")
            cycle
         end if
         if (len_trim(f%missing) > 0) cycle
         nfunc = nfunc + 1
         call xc_functional_ref_keys(trim(f%name), keys, nkeys)
         call expect(nkeys >= 1, "functional "//name//" cites no paper")
         do k = 1, nkeys
            call xc_ref_text(trim(keys(k)), text, doi, found)
            call expect(found .and. len(text) > 0, "functional "//name//": reference "//trim(keys(k))//" has no text")
         end do
      end do
      print '(a,i0,a)', "check_citations: ", nfunc, " functional names checked"
      call expect(nfunc > 60, "fewer runnable functionals than expected")
   end subroutine check_every_functional

   subroutine run_water(func)
      character(len=*), intent(in) :: func
      type(c_ptr) :: h
      integer(c_int) :: rc
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_force_cartesian(h, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_set_method(h, func//c_null_char, 3_c_int)
      rc = trc_set_convergence(h, 1.0e-8_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(h)
      call expect(rc == TRC_OK, "run_scf "//func)
      rc = trc_destroy(h)
   end subroutine run_water

   function citation_list() result(s)
      character(len=:), allocatable :: s
      character(kind=c_char), allocatable :: buf(:)
      character(kind=c_char) :: unused(1)
      integer(c_int) :: reqlen, rc
      integer :: i
      rc = trc_citation_list(unused, 0_c_int, reqlen)
      call expect(rc == TRC_OK .and. reqlen >= 1, "citation_list: length query")
      allocate (buf(reqlen))
      rc = trc_citation_list(buf, reqlen, reqlen)
      call expect(rc == TRC_OK, "citation_list: fill")
      call expect(buf(reqlen) == c_null_char, "citation_list: NUL-terminated")
      allocate (character(len=reqlen - 1) :: s)
      do i = 1, reqlen - 1
         s(i:i) = buf(i)
      end do
   end function citation_list

   subroutine check_buffer_protocol()
      character(kind=c_char) :: small(4)
      integer(c_int) :: reqlen, rc
      rc = trc_citation_list(small, 4_c_int, reqlen)
      call expect(rc == TRC_ERR_BADARG .and. reqlen > 4, "a buffer too small is refused, with the length needed")
   end subroutine check_buffer_protocol

   integer function count_lines(s)
      character(len=*), intent(in) :: s
      integer :: i
      count_lines = 0
      do i = 1, len(s)
         if (s(i:i) == new_line("a")) count_lines = count_lines + 1
      end do
   end function count_lines

   integer function count_unique_lines(s)
      !! Lines with the "[ n] " number stripped, counted once each.
      character(len=*), intent(in) :: s
      character(len=1024) :: seen(256)
      integer :: i0, i1, n, j
      logical :: dup
      n = 0
      i0 = 1
      do while (i0 <= len(s))
         i1 = index(s(i0:), new_line("a")) + i0 - 1
         dup = .false.
         do j = 1, n
            if (seen(j) == s(i0 + 5:i1 - 1)) dup = .true.
         end do
         if (.not. dup) then
            n = n + 1
            seen(n) = s(i0 + 5:i1 - 1)
         end if
         i0 = i1 + 1
      end do
      count_unique_lines = n
   end function count_unique_lines

end program check_citations
