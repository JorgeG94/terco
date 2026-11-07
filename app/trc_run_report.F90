submodule (trc_run) trc_run_report
   !! Console output of the `terco` command: the banner and the energy,
   !! gradient, charge, Fukui and Hessian summaries (bodies of `trc_run`).
   use pic_strings, only: to_string
   use trc_app, only: app_element_symbol, BOHR_PER_ANGSTROM
   use trc_vib, only: trc_vib_log_table
   use trc_thermo, only: trc_thermo_log_summary
   use trc_sym_run, only: sym_occupied_counts
   implicit none

contains

   module function method(r) result(m)
      type(run_t), intent(in) :: r
      character(len=:), allocatable :: m
      if (len(r%func) == 0) then
         m = merge("RHF", "UHF", r%mult == 1)
      else
         m = merge("RKS", "UKS", r%mult == 1)//"/"//r%func
      end if
   end function method

   module subroutine run_banner(r)
      type(run_t), intent(in) :: r
      call global_logger%info("")
      call global_logger%info("      -----------------------------------------------------------")
      call global_logger%info("     |                  t e r c o   "//VERSION)
      call global_logger%info("     |         GPU accelerated Quantum Chemistry ")
      call global_logger%info("     |     (c) 2026, Jorge Luis Galvez Vallejo, NCI  ")
      call global_logger%info("     |               Powered by Fortran        ")
      call global_logger%info("      -----------------------------------------------------------")
      call global_logger%info("")
      call global_logger%info("   geometry         : "//r%geom)
      if (len_trim(r%title) > 0) call global_logger%info("   title            : "//trim(r%title))
      call global_logger%info("   atoms            : "//to_string(r%natm)//"  ("//formula(r%z, r%natm)//")")
      block
         character(len=1024) :: logline
         write (logline, "(a,i0,a,i0,a,i0)") "   charge           : ", r%chrg, "   unpaired: ", r%uhf, &
            "   electrons: ", r%nelec
         call global_logger%info(trim(logline))
      end block
      if (r%nranks > 1) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a)") "   ranks            : ", r%nranks, " (one GPU each)"
            call global_logger%info(trim(logline))
         end block
      end if
      call global_logger%info("   method           : "//method(r))
      if (r%symmetry /= "off") call global_logger%info("   symmetry         : "//r%symmetry)
      call global_logger%info("   basis            : "//r%basis_path)
      if (r%rimp2) call global_logger%info("   auxiliary basis  : "//r%aux_path)
      ! The AO space itself (spherical or Cartesian) is not known until the
      ! basis is loaded, below -- `summary`'s own "AO basis" line reports
      ! it, for every run, once it is.
      if (len(r%restart) > 0) then
         call global_logger%info("   SCF              : read from "//r%restart//", not run")
      else
         call global_logger%info("   guess            : "//r%guess)
      end if
      block
         character(len=1024) :: logline
         write (logline, "(a,es8.1,a,es8.1,a,es8.1)") "   convergence      : dE ", r%econv, &
            "   |FDS-SDF| ", r%gconv, "   screening ", r%thresh
         call global_logger%info(trim(logline))
      end block
      if (len(r%func) > 0) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0)") "   grid level       : ", r%grid
            if (r%vv10) write (logline, "(a,i0,a,i0)") "   grid level       : ", r%grid, &
               "   VV10: ", r%nlc_grid
            call global_logger%info(trim(logline))
         end block
      end if
      if (len(r%pcm) > 0) then
         block
            character(len=1024) :: logline
            write (logline, "(a,a,a,f0.4)") "   solvation        : ", r%pcm, ", eps ", r%eps
            call global_logger%info(trim(logline))
         end block
      end if
   end subroutine run_banner

   function formula(z, natm) result(f)
      !! Hill order: with carbon, C then H then the rest alphabetically;
      !! without, everything alphabetically.
      !!
      !! Counted once into a per-element histogram with plain loops. The
      !! previous `any(z == ia)` test came back true for every element on
      !! nvfortran/aarch64 (a GH200, 1155-atom rna9), which printed all 118
      !! symbols while `count` still gave the right numbers.
      integer, intent(in) :: natm
      integer, intent(in) :: z(:)
      integer, parameter :: NELEM = 118
      integer :: ia, k, n
      integer :: cnt(NELEM), order(NELEM)
      character(len=:), allocatable :: f
      cnt = 0
      do ia = 1, natm
         if (z(ia) >= 1 .and. z(ia) <= NELEM) cnt(z(ia)) = cnt(z(ia)) + 1
      end do
      n = 0
      if (cnt(6) > 0) then
         n = n + 1
         order(n) = 6
         if (cnt(1) > 0) then
            n = n + 1
            order(n) = 1
         end if
      end if
      k = n
      do ia = 1, NELEM
         if (cnt(ia) == 0) cycle
         if (k > 0 .and. (ia == 6 .or. (ia == 1 .and. k == 2))) cycle
         n = n + 1
         order(n) = ia
      end do
      call sort_by_symbol(order(k + 1:n))
      f = ""
      do ia = 1, n
         f = f//trim(app_element_symbol(order(ia)))
         if (cnt(order(ia)) > 1) f = f//to_string(cnt(order(ia)))
      end do
   end function formula

   subroutine sort_by_symbol(a)
      integer, intent(inout) :: a(:)
      integer :: i, j, t
      do i = 2, size(a)
         t = a(i)
         j = i - 1
         do while (j >= 1)
            if (app_element_symbol(a(j)) <= app_element_symbol(t)) exit
            a(j + 1) = a(j)
            j = j - 1
         end do
         a(j + 1) = t
      end do
   end subroutine sort_by_symbol

   module subroutine run_summary(r)
      type(run_t), intent(in) :: r
      character(len=*), parameter :: bar = "          ---------------------------------------------------"
      call global_logger%info("")
      block
         character(len=1024) :: logline
         write (logline, "(a,i0,a,a)") "   AO basis         : ", r%nao, " functions, ", &
            merge("spherical", "Cartesian", r%sph /= 0)
         call global_logger%info(trim(logline))
      end block
      if (len(r%restart) > 0) then
         call global_logger%info("   SCF              : "//to_string(r%niter) &
                                 //trim(merge(" iterations, converged    ", " iterations, NOT converged", r%conv /= 0)) &
                                 //" (read from "//r%restart//")")
      else
         call global_logger%info("   SCF              : "//to_string(r%niter) &
                                 //trim(merge(" iterations, converged    ", " iterations, NOT converged", r%conv /= 0)))
      end if
      call global_logger%info("")
      call global_logger%info(bar)
      ! "TOTAL ENERGY"/"total energy" means the SAME thing in text and
      ! JSON: the double hybrid's weighted total when one ran, the SCF
      ! total otherwise -- JSON already worked this way; this merge
      ! makes the text banner/summary agree with it, rather than
      ! renaming either key (every external reader of terco's output,
      ! metalquicha included, looks for "total energy"/"TOTAL ENERGY"
      ! by that name). "SCF ENERGY"/"scf energy" is new: always
      ! present, always the SCF-level value alone, so a caller wanting
      ! that specifically no longer has to infer it from context.
      call line("TOTAL ENERGY", merge(r%etot_dh, r%etot, r%dh .and. r%conv /= 0))
      call line("SCF ENERGY", r%etot)
      call frontier_lines(r)
      if (r%rimp2 .and. r%conv /= 0) then
         if (r%dh) then
            call line("PT2 CORRELATION", r%eos + r%ess)
            call line("DOUBLE HYBRID TOTAL", r%etot_dh)
         else
            call line("MP2 CORRELATION", r%eos + r%ess)
            call line("TOTAL MP2 ENERGY", r%etot + r%eos + r%ess)
         end if
      end if
      call global_logger%info(bar)
      call global_logger%info("")
      call global_logger%info("   ::::::::::::::::::::::::::::::::::::::::::::::::::::::")
      call global_logger%info("   ::"//repeat(" ", 21)//"SUMMARY"//repeat(" ", 22)//"::")
      call global_logger%info("   ::::::::::::::::::::::::::::::::::::::::::::::::::::::")
      call row("total energy", merge(r%etot_dh, r%etot, r%dh .and. r%conv /= 0))
      call row("scf energy", r%etot)
      call row("nuclear repulsion", r%enuc)
      call row("one-electron energy", r%eone)
      call row("two-electron energy", r%etwo)
      if (len(r%func) > 0) call row("exchange-correlation", r%exc)
      if (len(r%pcm) > 0) call row("solvation (E_diel)", r%epcm)
      if (r%rimp2 .and. r%conv /= 0) then
         call row(merge("PT2 opposite spin", "MP2 opposite spin", r%dh), r%eos)
         call row(merge("PT2 same spin    ", "MP2 same spin    ", r%dh), r%ess)
         if (r%dh) call row("double hybrid total", r%etot_dh)
      end if
      call global_logger%info("   ::::::::::::::::::::::::::::::::::::::::::::::::::::::")
   end subroutine run_summary

   module subroutine run_symmetry_summary(r)
      !! What the symmetry request did: the group detected and used, how far the
      !! atoms moved, the geometry the SCF actually ran on (the standard frame,
      !! not the input's), the electrons per irrep, and a label on each orbital.
      !! Printed only when a symmetry request ran, so a run with symmetry off is
      !! unchanged.
      type(run_t), intent(in) :: r
      character(len=1024) :: logline
      integer :: ia, ir, s, m, nmo, nocc, m0, m1, ns
      integer, allocatable :: cnt(:, :)
      character(len=64) :: tag

      call global_logger%info("")
      call global_logger%info("   symmetry         : requested "//r%symmetry//", detected "//trim(r%sym_detected)// &
                              ", used "//trim(r%sym_used))
      if (.not. r%sym_active) then
         call global_logger%info("                      no symmetry in the SCF (the group used is C1)")
         return
      end if
      write (logline, "(a,es9.2,a)") "   symmetrisation   : largest atom displacement ", r%sym_residual, " Bohr"
      call global_logger%info(trim(logline))
      call global_logger%info("   NOTE: the coordinates below, the gradient and the orbitals are in the STANDARD FRAME of "// &
                              trim(r%sym_used)//":")
      call global_logger%info("         the molecule is symmetrised, moved to its centre of nuclear charge and rotated "// &
                              "to the group's axes")
      call global_logger%info("         (--json symmetry.frame rotation and symmetry.frame origin map the input to it).")
      call global_logger%info("   geometry, standard frame (Angstrom)")
      do ia = 1, r%natm
         write (logline, "(a,i4,2x,a2,3f16.8)") "   ", ia, app_element_symbol(r%z(ia)), r%geom_std(:, ia)/BOHR_PER_ANGSTROM
         call global_logger%info(trim(logline))
      end do
      if (.not. allocated(r%sym_mo_irrep)) return
      ns = r%nspin
      nmo = size(r%sym_occ, 1)
      allocate (cnt(r%sym_nirrep, ns))
      do s = 1, ns
         cnt(:, s) = sym_occupied_counts(r%sym_mo_irrep(1:nmo, s), r%sym_occ(:, s) > 0.0_dp, r%sym_nirrep)
      end do
      if (ns == 1) then
         call global_logger%info("   irrep  SO functions  occupied orbitals")
      else
         call global_logger%info("   irrep  SO functions  occupied orbitals (alpha / beta)")
      end if
      do ir = 1, r%sym_nirrep
         if (ns == 1) then
            write (logline, "(a,a6,i12,i12)") "   ", trim(r%sym_irrep_names(ir)), r%sym_irrep_dim(ir), cnt(ir, 1)
         else
            write (logline, "(a,a6,i12,i8,a,i0)") "   ", trim(r%sym_irrep_names(ir)), r%sym_irrep_dim(ir), cnt(ir, 1), &
               " / ", cnt(ir, 2)
         end if
         call global_logger%info(trim(logline))
      end do
      do s = 1, ns
         nocc = count(r%sym_occ(:, s) > 0.0_dp)
         ! Every orbital of a small molecule; a window around the frontier otherwise.
         m0 = 1
         m1 = nmo
         if (nmo > 50) then
            m0 = max(1, nocc - 14)
            m1 = min(nmo, nocc + 10)
         end if
         call global_logger%info("")
         tag = "   orbitals"
         if (ns == 2) tag = merge("   orbitals (alpha)", "   orbitals (beta) ", s == 1)
         if (m0 > 1 .or. m1 < nmo) tag = trim(tag)//", window around the HOMO"
         call global_logger%info(trim(tag)//": index  irrep  energy (Eh)  occupation")
         do m = m0, m1
            write (logline, "(a,i8,2x,a6,f18.10,f8.1)") "   ", m, trim(r%sym_irrep_names(max(1, r%sym_mo_irrep(m, s)))), &
               r%sym_eps(m, s), r%sym_occ(m, s)
            call global_logger%info(trim(logline))
         end do
      end do
   end subroutine run_symmetry_summary

   module subroutine run_gradient_summary(r)
      !! Per-atom analytic nuclear gradient (Eh/Bohr), plus the RMS and the
      !! largest per-atom norm -- the two scalars an optimiser's own
      !! convergence check usually gates on. `r%gradient` is already
      !! filled by the time this is called (trc_main's own job, right
      !! after `trc_gradient`).
      type(run_t), intent(in) :: r
      real(dp) :: rms, gmax
      integer :: ia
      call global_logger%info("")
      call global_logger%info("   gradient (Eh/Bohr)")
      call global_logger%info("   ----------------------------------------------------------")
      do ia = 1, r%natm
         block
            character(len=1024) :: logline
            write (logline, "(a,i4,2x,a2,3f16.8)") "   ", ia, app_element_symbol(r%z(ia)), r%gradient(:, ia)
            call global_logger%info(trim(logline))
         end block
      end do
      rms = sqrt(sum(r%gradient**2)/real(3*r%natm, dp))
      gmax = sqrt(maxval(sum(r%gradient**2, dim=1)))
      call global_logger%info("   ----------------------------------------------------------")
      block
         character(len=1024) :: logline
         write (logline, "(a,es12.4)") "   RMS gradient     : ", rms
         call global_logger%info(trim(logline))
      end block
      block
         character(len=1024) :: logline
         write (logline, "(a,es12.4)") "   max |gradient|   : ", gmax
         call global_logger%info(trim(logline))
      end block
   end subroutine run_gradient_summary

   module subroutine run_charges_summary(r)
      !! Per-atom charges (e) of each scheme --charges named, and their
      !! sum, which matches the molecular charge to round-off: Hirshfeld's
      !! are renormalised to Tr(D S) (trc_hirshfeld.F90), CHELPG's are fitted
      !! under that constraint (trc_chelpg.F90). Both arrays are filled by
      !! the time this is called (trc_main's own job, right after the
      !! trc_*_charges calls).
      type(run_t), intent(in) :: r
      if (allocated(r%hirshfeld_charges)) call charges_table(r, "hirshfeld", r%hirshfeld_charges)
      if (allocated(r%chelpg_charges)) call charges_table(r, "chelpg", r%chelpg_charges)
   end subroutine run_charges_summary

   subroutine charges_table(r, scheme, charges)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: scheme
      real(dp), intent(in) :: charges(:)
      character(len=1024) :: logline
      integer :: ia
      call global_logger%info("")
      call global_logger%info("   "//scheme//" charges (e)")
      call global_logger%info("   ----------------------------------------------------------")
      do ia = 1, r%natm
         write (logline, "(a,i4,2x,a2,f16.8)") "   ", ia, app_element_symbol(r%z(ia)), charges(ia)
         call global_logger%info(trim(logline))
      end do
      call global_logger%info("   ----------------------------------------------------------")
      write (logline, "(a,es12.4)") "   sum              : ", sum(charges)
      call global_logger%info(trim(logline))
   end subroutine charges_table

   module subroutine run_fukui_summary(r)
      !! Per-atom Hirshfeld-condensed Fukui indices (f+, f-, f0) and their
      !! sums, each of which should be one electron to round-off --
      !! `r%fukui_f_plus`/`fukui_f_minus`/`fukui_f_zero` are already
      !! filled by the time this is called (trc_main's own job, right
      !! after `trc_fukui_hirshfeld`).
      type(run_t), intent(in) :: r
      integer :: ia
      call global_logger%info("")
      call global_logger%info("   Hirshfeld-condensed Fukui indices")
      call global_logger%info("   ----------------------------------------------------------")
      call global_logger%info("     atom            f+          f-          f0")
      do ia = 1, r%natm
         block
            character(len=1024) :: logline
            write (logline, "(a,i4,2x,a2,3f12.6)") "   ", ia, app_element_symbol(r%z(ia)), &
               r%fukui_f_plus(ia), r%fukui_f_minus(ia), r%fukui_f_zero(ia)
            call global_logger%info(trim(logline))
         end block
      end do
      call global_logger%info("   ----------------------------------------------------------")
      block
         character(len=1024) :: logline
         write (logline, "(a,3f12.6)") "   sum             ", sum(r%fukui_f_plus), &
            sum(r%fukui_f_minus), sum(r%fukui_f_zero)
         call global_logger%info(trim(logline))
      end block
      call global_logger%info("")
      call global_logger%info("   f+ is where the molecule accepts charge, f- where it gives "// &
                              "it up; f0 is their average.")
   end subroutine run_fukui_summary

   !
   ! The frequency table and thermochemistry summary for --hessian,
   ! through trc_vib's/trc_thermo's own logger helpers -- `r%freq_cm1`
   ! etc. and `r%thermo` are already filled by trc_main itself, right
   ! after trc_vibrational_analysis/trc_thermochemistry succeed.
   !
   module subroutine run_hessian_summary(r)
      type(run_t), intent(in) :: r
      call global_logger%info("")
      call global_logger%info("   vibrational analysis and thermochemistry (--hessian)")
      call global_logger%info("   ----------------------------------------------------------")
      block
         character(len=1024) :: logline
         write (logline, "(a,f0.4,a,a,a,i0,a,i0)") "   step ", r%hess_step, " Bohr   masses: ", r%masses_kind, &
            "   symmetry number: ", r%sym_number, " (not determined automatically)   mult: ", r%mult
         call global_logger%info(trim(logline))
      end block
      call trc_vib_log_table(r%freq_cm1, r%redmass_amu, r%force_const)
      call trc_thermo_log_summary(r%thermo)
   end subroutine run_hessian_summary

   subroutine frontier_lines(r)
      !! HOMO/LUMO rows of the energy box: one pair for a restricted run,
      !! an alpha and a beta pair for an unrestricted one.
      type(run_t), intent(in) :: r
      character(len=8), parameter :: tag(2) = ["(alpha) ", "(beta)  "]
      integer :: s
      do s = 1, r%nspin
         if (r%nspin == 1) then
            if (r%has_homo(s)) call line("HOMO", real(r%homo(s), c_double))
            if (r%has_lumo(s)) call line("LUMO", real(r%lumo(s), c_double))
         else
            if (r%has_homo(s)) call line("HOMO "//trim(tag(s)), real(r%homo(s), c_double))
            if (r%has_lumo(s)) call line("LUMO "//trim(tag(s)), real(r%lumo(s), c_double))
         end if
      end do
   end subroutine frontier_lines

   subroutine line(label, e)
      character(len=*), intent(in) :: label
      real(c_double), intent(in) :: e
      block
         character(len=1024) :: logline
         write (logline, "(a,a20,f22.12,a)") "          | ", label, e, " Eh    |"
         call global_logger%info(trim(logline))
      end block
   end subroutine line

   subroutine row(label, e)
      character(len=*), intent(in) :: label
      real(c_double), intent(in) :: e
      block
         character(len=1024) :: logline
         write (logline, "(a,a22,f22.12,a)") "   :: ", label, e, " Eh  ::"
         call global_logger%info(trim(logline))
      end block
   end subroutine row

end submodule trc_run_report
