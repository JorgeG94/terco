!
! RI-CCSD and RI-CCSD(T) through the real `terco` binary: the flags --riccsd, --riccsd-t, --cc-econv, --cc-tconv,
! --cc-iterations, --cc-diis, --cc-ladder, --cc-triples, --cc-triples-budget, the [correlation] method and [riccsd] keys of a TOML
! file, the summary and the --json output, and every combination the command line refuses.
!
! Reference. Water, cc-pVDZ with cc-pVDZ-RIFIT, frozen core 1 (the geometry of test/hirshfeld_water.xyz): pyscf's
! DF-CCSD(T) in test/pyscf_refs/rccsd_ints_w1.ref, as check_capi_rccsd uses it, with the same tolerances:
! E_CCSD correlation 1e-10, E(T) 3e-12, the CCSD(T) total energy 3e-11. The runs converge the SCF with
! --econv 1e-13 --gconv 1e-11 and the coupled-cluster iteration with --cc-econv 1e-12 --cc-tconv 1e-10.
!
! The summary prints the energies with 12 decimals, so its values are compared with the tolerance above plus the
! 5e-13 of that rounding. The JSON carries every digit.
!
! Cases:
!   - --riccsd-t: JSON energies against pyscf, and the summary's lines against pyscf and the JSON;
!   - --riccsd alone has no (T) keys; --riccsd-t alone is the same CCSD as --riccsd (it implies it);
!   - --frozen -1 (the default: the core counted by element) is --frozen 1 for water;
!   - a TOML file ([correlation] method = "ri-ccsd(t)", [riccsd] keys) gives the numbers of its flags, "ri-ccsd" gives
!     no (T), and the plain "ccsd" and "ccsd(t)" are refused with the accepted names,
!     and the JSON "input" echo names the method;
!   - the convergence flags reach the solver: --cc-iterations 2 stops after two iterations, a loose --cc-econv /
!     --cc-tconv stops earlier and lands elsewhere, --cc-diis 1 takes more iterations than the default;
!   - --cc-ladder and --cc-triples-budget do not change the energy;
!   - a CCSD that does not converge: exit status 2, "NOT CONVERGED" in the summary, "ri-ccsd converged": false and no
!     total energy in the JSON, and the message on the error stream;
!   - a (T) that refuses after a converged RI-CCSD (the tiled routine, a budget of 1e-6 GiB = 1073 bytes, below the
!     3072 bytes one water triple needs), by flags and by a TOML file: exit status 1, the RI-CCSD energies in the box,
!     the summary and the JSON, a "(T) FAILED" line, the reason on the error stream and in the JSON key
!     "ri-ccsd triples error", and no (T) correction or RI-CCSD(T) total anywhere; --cc-triples reference with the
!     same budget succeeds (the reference routine has no budget); the routine is always named, so the case means the
!     same in a host build and in a device build;
!   - --silent prints the summary but not the iteration table;
!   - each refusal: --rimp2, a functional, --uhf, --grad, --hessian, --opt, an ECP, the --cc-* options without
!     --riccsd, a bad --cc-ladder, --cc-triples-budget without (T), bad values: a non-zero exit status and the
!     reason on the error stream.
!
! Arguments: the terco executable and a scratch directory. A build whose density-fitting ceilings are below
! the cc-pVDZ / cc-pVDZ-RIFIT shells prints SKIP. A test-drive suite.
!
module check_cli_rccsd_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use rccsd_ref_energies, only: rccsd_refs_t, read_rccsd_refs
   implicit none
   private
   public :: collect_cli_rccsd, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns; the runs happen in its `case/` subdirectory.
   real(dp), parameter :: TOL_CORR = 1.0e-10_dp
   real(dp), parameter :: TOL_T = 3.0e-12_dp
   real(dp), parameter :: TOL_TOTAL = 3.0e-11_dp
   real(dp), parameter :: PRINT_ROUNDING = 5.0e-13_dp
   !! half a unit of the 12th decimal the summary prints
   real(dp), parameter :: TOL_SAME = 1.0e-12_dp
   !! two runs of the same numbers reached by different settings (a ladder form, a tile size, a TOML file)
   character(len=*), parameter :: REF = "rccsd_ints_w1.ref"
   character(len=*), parameter :: COMMON = ' water.xyz -b basis_sets/cc-pvdz.json --aux basis_sets/cc-pvdz-rifit.json'// &
                                          ' --econv 1e-13 --gconv 1e-11'
   character(len=*), parameter :: TIGHT = ' --cc-econv 1e-12 --cc-tconv 1e-10'

contains

   subroutine collect_cli_rccsd(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("--riccsd-t: JSON and summary match pyscf", test_triples), &
                  new_unittest("--riccsd alone has no (T); --riccsd-t implies --riccsd", test_implies), &
                  new_unittest("--frozen -1 counts the water core", test_frozen_default), &
                  new_unittest("a TOML file gives the numbers of its flags", test_toml), &
                  new_unittest("the convergence flags reach the solver", test_convergence_flags), &
                  new_unittest("--cc-ladder and --cc-triples-budget do not change the energy", test_algorithm_flags), &
                  new_unittest("a CCSD that does not converge is unmistakable and exits 2", test_noconv), &
                  new_unittest("a (T) that refuses after a converged RI-CCSD keeps the RI-CCSD and exits 1", &
                               test_triples_failed), &
                  new_unittest("--silent prints the summary and no iteration table", test_silent), &
                  new_unittest("combinations the method does not cover are refused", test_refusals), &
                  new_unittest("bad --cc-* values are refused", test_bad_values) &
                  ]
   end subroutine collect_cli_rccsd

   function supported() result(ok)
      logical :: ok
      ok = DF_LMAX >= 2 .and. DF_LMAX_AUX >= 3
      if (.not. ok) print "(a)", "SKIP: this build's density-fitting ceilings are below cc-pVDZ-RIFIT"
   end function supported

   subroutine setup()
      !! case/ with the geometry and the two basis files
      integer :: rc
      call run('mkdir -p '//work_dir//'/case/basis_sets', rc)
      call run('cp -f hirshfeld_water.xyz '//work_dir//'/case/water.xyz', rc)
      call run('cp -f hi.xyz '//work_dir//'/case/hi.xyz', rc)
      call run('cp -f basis_sets/cc-pvdz.json '//work_dir//'/case/basis_sets/cc-pvdz.json', rc)
      call run('cp -f basis_sets/cc-pvdz-rifit.json '//work_dir//'/case/basis_sets/cc-pvdz-rifit.json', rc)
   end subroutine setup

   subroutine run(cmd, rc)
      character(len=*), intent(in) :: cmd
      integer, intent(out) :: rc
      rc = 0
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
   end subroutine run

   subroutine terco(args, name, rc)
      !! `terco args` in case/, stdout and stderr to case/<name>.out, the JSON to case/<name>.json
      character(len=*), intent(in) :: args, name
      integer, intent(out) :: rc
      call run('cd '//work_dir//'/case && '//cli_exe//args//' --json '//name//'.json > '//name//'.out 2>&1', rc)
   end subroutine terco

   function slurp(path) result(text)
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: text
      integer :: u, ios, n
      text = ""
      open (newunit=u, file=path, status="old", action="read", access="stream", form="unformatted", iostat=ios)
      if (ios /= 0) return
      inquire (unit=u, size=n)
      deallocate (text)
      allocate (character(len=n) :: text)
      read (u, iostat=ios) text
      close (u)
   end function slurp

   function out_text(name) result(text)
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: text
      text = slurp(work_dir//'/case/'//name//'.out')
   end function out_text

   subroutine json_real(name, key, value, found)
      character(len=*), intent(in) :: name, key
      real(dp), intent(out) :: value
      logical, intent(out) :: found
      type(json_file) :: json
      value = 0.0_dp
      call json%initialize()
      call json%load(filename=work_dir//'/case/'//name//'.json')
      call json%get(key, value, found)
      call json%destroy()
   end subroutine json_real

   subroutine json_int(name, key, value, found)
      character(len=*), intent(in) :: name, key
      integer, intent(out) :: value
      logical, intent(out) :: found
      type(json_file) :: json
      value = 0
      call json%initialize()
      call json%load(filename=work_dir//'/case/'//name//'.json')
      call json%get(key, value, found)
      call json%destroy()
   end subroutine json_int

   subroutine json_logical(name, key, value, found)
      character(len=*), intent(in) :: name, key
      logical, intent(out) :: value, found
      type(json_file) :: json
      value = .false.
      call json%initialize()
      call json%load(filename=work_dir//'/case/'//name//'.json')
      call json%get(key, value, found)
      call json%destroy()
   end subroutine json_logical

   subroutine json_string(name, key, value, found)
      character(len=*), intent(in) :: name, key
      character(len=:), allocatable, intent(out) :: value
      logical, intent(out) :: found
      type(json_file) :: json
      call json%initialize()
      call json%load(filename=work_dir//'/case/'//name//'.json')
      call json%get(key, value, found)
      call json%destroy()
   end subroutine json_string

   subroutine summary_value(name, label, value, found)
      !! The number after `label` on a line of the energy box of case/<name>.out
      character(len=*), intent(in) :: name, label
      real(dp), intent(out) :: value
      logical, intent(out) :: found
      character(len=:), allocatable :: text
      integer :: i, j, ios
      value = 0.0_dp
      found = .false.
      text = out_text(name)
      i = index(text, "| "//repeat(" ", max(0, 24 - len(label)))//label)
      if (i == 0) return
      j = i + 2 + max(len(label), 24)
      read (text(j:), *, iostat=ios) value
      found = ios == 0
   end subroutine summary_value

   function has_text(name, needle) result(found)
      character(len=*), intent(in) :: name, needle
      logical :: found
      found = index(out_text(name), needle) > 0
   end function has_text

   function summary_rows_aligned(name) result(aligned)
      !! Every row of the SUMMARY block (the lines that start "   ::") is as wide as the first one: a row with a
      !! longer label field sticks out of the box.
      character(len=*), intent(in) :: name
      logical :: aligned
      character(len=:), allocatable :: text
      integer :: pos, eol, width, first_width, nrows

      text = out_text(name)
      aligned = .true.
      first_width = -1
      nrows = 0
      pos = 1
      do while (pos <= len(text))
         eol = index(text(pos:), new_line("a"))
         if (eol == 0) then
            eol = len(text) - pos + 2
         end if
         if (eol - 1 >= 5) then
            if (text(pos:pos + 4) == "   ::") then
               width = len_trim(text(pos:pos + eol - 2))
               nrows = nrows + 1
               if (first_width < 0) first_width = width
               if (width /= first_width) aligned = .false.
            end if
         end if
         pos = pos + eol
      end do
      if (nrows < 10) aligned = .false.
   end function summary_rows_aligned

   subroutine require_run(error, args, name, ok)
      !! The CLI run succeeded (exit 0)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: args, name
      logical, intent(out) :: ok
      integer :: rc
      call terco(args, name, rc)
      ok = rc == 0
      if (.not. ok) call test_failed(error, "terco"//args//" exited "//to_string(rc)//" (see case/"//name//".out)")
   end subroutine require_run

   subroutine compare_to_reference(error, name, triples, label)
      !! The JSON of case/<name>.json against pyscf: correlation, (T) when `triples`, totals
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: name, label
      logical, intent(in) :: triples
      type(rccsd_refs_t) :: refs
      real(dp) :: e_ccsd, e_t, e_total, e_total_t
      logical :: f1, f2, f3, f4

      call read_rccsd_refs(REF, refs)
      call json_real(name, "ri-ccsd correlation energy", e_ccsd, f1)
      call json_real(name, "ri-ccsd total energy", e_total, f2)
      call check(error, f1 .and. f2, label//": no ccsd correlation/total energy in case/"//name//".json")
      if (allocated(error)) return
      call check(error, abs(e_ccsd - refs%e_ccsd_corr) < TOL_CORR, label//": E_CCSD "//to_string(e_ccsd, "f20.12")// &
                 " vs pyscf "//to_string(refs%e_ccsd_corr, "f20.12")//", |difference| "// &
                 to_string(abs(e_ccsd - refs%e_ccsd_corr), "es9.2"))
      if (allocated(error)) return
      call check(error, abs(e_total - (refs%e_scf + refs%e_ccsd_corr)) < TOL_TOTAL, label//": CCSD total "// &
                 to_string(e_total, "f20.12"))
      if (allocated(error) .or. .not. triples) return
      call json_real(name, "ri-ccsd triples correction", e_t, f3)
      call json_real(name, "ri-ccsd triples total energy", e_total_t, f4)
      call check(error, f3 .and. f4, label//": no (T) keys in case/"//name//".json")
      if (allocated(error)) return
      call check(error, abs(e_t - refs%e_t) < TOL_T, label//": E(T) "//to_string(e_t, "f20.14")//" vs pyscf "// &
                 to_string(refs%e_t, "f20.14")//", |difference| "//to_string(abs(e_t - refs%e_t), "es9.2"))
      if (allocated(error)) return
      call check(error, abs(e_total_t - refs%e_ccsd_t) < TOL_TOTAL, label//": CCSD(T) total "// &
                 to_string(e_total_t, "f20.12")//" vs pyscf "//to_string(refs%e_ccsd_t, "f20.12")//", |difference| "// &
                 to_string(abs(e_total_t - refs%e_ccsd_t), "es9.2"))
   end subroutine compare_to_reference

   subroutine test_triples(error)
      type(error_type), allocatable, intent(out) :: error
      type(rccsd_refs_t) :: refs
      real(dp) :: v, e_json
      logical :: ok, found
      integer :: iterations

      if (.not. supported()) return
      call setup()
      call require_run(error, COMMON//' --frozen 1 --riccsd-t'//TIGHT, 'triples', ok)
      if (.not. ok) return
      call compare_to_reference(error, 'triples', .true., "--riccsd-t")
      if (allocated(error)) return
      call read_rccsd_refs(REF, refs)
      call summary_value('triples', "RI-CCSD CORRELATION", v, found)
      call json_real('triples', "ri-ccsd correlation energy", e_json, ok)
      call check(error, found .and. abs(v - e_json) < PRINT_ROUNDING + 1.0e-14_dp, "the summary's RI-CCSD CORRELATION "// &
                 to_string(v, "f20.12")//" is not the JSON's "//to_string(e_json, "f20.12"))
      if (allocated(error)) return
      call summary_value('triples', "RI-CCSD (T) CORRECTION", v, found)
      call check(error, found .and. abs(v - refs%e_t) < TOL_T + PRINT_ROUNDING, "the summary's (T) CORRECTION "// &
                 to_string(v, "f20.12")//" vs pyscf "//to_string(refs%e_t, "f20.12"))
      if (allocated(error)) return
      call summary_value('triples', "TOTAL RI-CCSD(T) ENERGY", v, found)
      call check(error, found .and. abs(v - refs%e_ccsd_t) < TOL_TOTAL + PRINT_ROUNDING, &
                 "the summary's TOTAL RI-CCSD(T) ENERGY "//to_string(v, "f20.12")//" vs pyscf "//to_string(refs%e_ccsd_t, "f20.12"))
      if (allocated(error)) return
      call summary_value('triples', "RI-MP2 CORRELATION", v, found)
      call check(error, found .and. abs(v - refs%e_mp2) < 1.0e-10_dp + PRINT_ROUNDING, "the summary's RI-MP2 CORRELATION "// &
                 to_string(v, "f20.12")//" vs pyscf "//to_string(refs%e_mp2, "f20.12"))
      if (allocated(error)) return
      ! The MP2 of an RI-CCSD run is the density-fitted one and the JSON names it so.
      call json_real('triples', "ri-mp2 correlation energy", v, found)
      call check(error, found .and. abs(v - refs%e_mp2) < 1.0e-10_dp, 'the JSON "ri-mp2 correlation energy" '// &
                 to_string(v, "f20.12")//" vs pyscf "//to_string(refs%e_mp2, "f20.12"))
      if (allocated(error)) return
      call json_int('triples', "ri-ccsd iterations", iterations, found)
      call json_logical('triples', "ri-ccsd converged", ok, found)
      call check(error, found .and. ok .and. iterations > 3, "ri-ccsd converged/iterations in the JSON: "//to_string(iterations))
      if (allocated(error)) return
      call json_real('triples', "ri-ccsd t1 diagnostic", v, found)
      call check(error, found .and. abs(v - refs%t1_diag) < 1.0e-9_dp, "ri-ccsd t1 diagnostic "//to_string(v, "f12.8")// &
                 " vs pyscf "//to_string(refs%t1_diag, "f12.8"))
      if (allocated(error)) return
      ! the spaces: water cc-pVDZ is 24 orbitals, 5 occupied with the oxygen 1s frozen; cc-pVDZ-RIFIT is 84 functions
      call json_int('triples', "ri-ccsd frozen core orbitals", iterations, found)
      call check(error, found .and. iterations == 1, "ri-ccsd frozen core orbitals is "//to_string(iterations))
      if (allocated(error)) return
      call json_int('triples', "ri-ccsd active occupied", iterations, found)
      call check(error, found .and. iterations == 4, "ri-ccsd active occupied is "//to_string(iterations))
      if (allocated(error)) return
      call json_int('triples', "ri-ccsd virtual orbitals", iterations, found)
      call check(error, found .and. iterations == 19, "ri-ccsd virtual orbitals is "//to_string(iterations))
      if (allocated(error)) return
      call json_int('triples', "ri-ccsd auxiliary functions", iterations, found)
      call check(error, found .and. iterations == 84, "ri-ccsd auxiliary functions is "//to_string(iterations))
      if (allocated(error)) return
      ! the summary rows carry the JSON's names for the same quantities
      call check(error, has_text('triples', "ri-ccsd(t) total") .and. has_text('triples', "ri-ccsd active occ.") &
                 .and. has_text('triples', "ri-ccsd auxiliary") .and. summary_rows_aligned('triples'), &
                 "the SUMMARY lacks the RI-CCSD rows, or they are not as wide as the other rows of the block")
   end subroutine test_triples

   subroutine test_implies(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e_ccsd_only, e_ccsd_t, v
      logical :: ok, found

      if (.not. supported()) return
      call setup()
      call require_run(error, COMMON//' --frozen 1 --riccsd'//TIGHT, 'ccsd_only', ok)
      if (.not. ok) return
      call compare_to_reference(error, 'ccsd_only', .false., "--riccsd")
      if (allocated(error)) return
      call json_real('ccsd_only', "ri-ccsd triples correction", v, found)
      call check(error, .not. found, "--riccsd wrote a (T) correction")
      if (allocated(error)) return
      call check(error, .not. has_text('ccsd_only', "(T) CORRECTION"), "--riccsd printed a (T) line")
      if (allocated(error)) return
      call require_run(error, COMMON//' --frozen 1 --riccsd-t'//TIGHT, 'ccsd_t', ok)
      if (.not. ok) return
      call json_real('ccsd_only', "ri-ccsd correlation energy", e_ccsd_only, found)
      call json_real('ccsd_t', "ri-ccsd correlation energy", e_ccsd_t, ok)
      call check(error, found .and. ok .and. abs(e_ccsd_only - e_ccsd_t) < TOL_SAME, &
                 "--riccsd-t is not the same CCSD as --riccsd: "//to_string(e_ccsd_t - e_ccsd_only, "es9.2"))
   end subroutine test_implies

   subroutine test_frozen_default(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: ok

      if (.not. supported()) return
      call setup()
      call require_run(error, COMMON//' --riccsd-t'//TIGHT, 'byelement', ok)
      if (.not. ok) return
      call compare_to_reference(error, 'byelement', .true., "--riccsd-t, the default --frozen -1")
      if (allocated(error)) return
      call require_run(error, COMMON//' --frozen 0 --riccsd'//TIGHT, 'nofrozen', ok)
      if (.not. ok) return
      call check_not_reference(error, 'nofrozen')
   end subroutine test_frozen_default

   subroutine check_not_reference(error, name)
      !! --frozen 0 must move the energy: an ignored --frozen would still match the frozen-core reference
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: name
      type(rccsd_refs_t) :: refs
      real(dp) :: e
      logical :: found
      call read_rccsd_refs(REF, refs)
      call json_real(name, "ri-ccsd correlation energy", e, found)
      call check(error, found .and. abs(e - refs%e_ccsd_corr) > 1.0e-4_dp, "--frozen 0 gives the frozen-core energy")
   end subroutine check_not_reference

   subroutine write_toml(path, lines)
      character(len=*), intent(in) :: path, lines(:)
      integer :: u, i
      open (newunit=u, file=path, status='replace', action='write')
      do i = 1, size(lines)
         write (u, '(a)') trim(lines(i))
      end do
      close (u)
   end subroutine write_toml

   subroutine test_toml(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=64), allocatable :: lines(:)
      character(len=:), allocatable :: method
      real(dp) :: e_flags, e_file, t_flags, t_file, v
      logical :: ok, found, found2
      integer :: iterations

      if (.not. supported()) return
      call setup()
      lines = [character(len=64) :: '[molecule]', 'geometry = "water.xyz"', '[basis]', &
               'main = "basis_sets/cc-pvdz.json"', 'auxiliary = "basis_sets/cc-pvdz-rifit.json"', &
               '[scf]', 'energy_tol = 1e-13', 'commutator_tol = 1e-11', &
               '[correlation]', 'method = "ri-ccsd(t)"', 'frozen_core = -1', 'aux_block = 40', &
               '[riccsd]', 'energy_tol = 1e-12', 'amplitude_tol = 1e-10', 'max_iterations = 50', 'diis_size = 8', &
               'ladder = "direct"', 'triples_budget = 0.5']
      call write_toml(work_dir//'/case/run.toml', lines)
      call terco_file('run.toml', 'toml', ok)
      if (.not. ok) then
         call test_failed(error, "terco run.toml failed (see case/toml.out)")
         return
      end if
      call compare_to_reference(error, 'toml', .true., "TOML method = ri-ccsd(t)")
      if (allocated(error)) return
      call require_run(error, COMMON//' --frozen -1 --riccsd-t'//TIGHT, 'flags', ok)
      if (.not. ok) return
      call json_real('flags', "ri-ccsd correlation energy", e_flags, found)
      call json_real('toml', "ri-ccsd correlation energy", e_file, found2)
      call json_real('flags', "ri-ccsd triples correction", t_flags, ok)
      call json_real('toml', "ri-ccsd triples correction", t_file, ok)
      call check(error, found .and. found2 .and. abs(e_flags - e_file) < TOL_SAME .and. abs(t_flags - t_file) < TOL_SAME, &
                 "the TOML file and its flags differ: E_CCSD "//to_string(e_file - e_flags, "es9.2")//", E(T) "// &
                 to_string(t_file - t_flags, "es9.2"))
      if (allocated(error)) return
      call json_string('toml', "input.correlation.method", method, found)
      call check(error, found, "no input.correlation.method in the JSON echo")
      if (allocated(error)) return
      call check(error, method == "ri-ccsd(t)", "the echo says method '"//method//"'")
      if (allocated(error)) return
      call json_int('toml', "input.riccsd.max_iterations", iterations, found)
      call check(error, found .and. iterations == 50, "the echo lacks the [riccsd] keys")
      if (allocated(error)) return
      ! method = "ri-ccsd": no (T)
      lines(10) = 'method = "ri-ccsd"'
      call write_toml(work_dir//'/case/run_ccsd.toml', lines(1:18))
      call terco_file('run_ccsd.toml', 'toml_ccsd', ok)
      if (.not. ok) then
         call test_failed(error, "terco run_ccsd.toml failed (see case/toml_ccsd.out)")
         return
      end if
      call json_real('toml_ccsd', "ri-ccsd correlation energy", v, found)
      call json_real('toml_ccsd', "ri-ccsd triples correction", t_file, found2)
      call check(error, found .and. .not. found2, "method = ri-ccsd: the (T) keys are "//merge("present", "absent ", found2))
      if (allocated(error)) return
      ! the plain names are refused, naming the accepted ones: nobody should mistake this for a conventional CCSD
      lines(10) = 'method = "ccsd"'
      call write_toml(work_dir//'/case/run_plain.toml', lines(1:18))
      call terco_file('run_plain.toml', 'toml_plain', ok)
      call check(error, .not. ok .and. has_text('toml_plain', "'ri-ccsd', 'ri-ccsd(t)'"), &
                 "[correlation] method = ccsd was accepted or the refusal does not name the RI names (case/toml_plain.out)")
      if (allocated(error)) return
      lines(10) = 'method = "ccsd(t)"'
      call write_toml(work_dir//'/case/run_plain_t.toml', lines(1:18))
      call terco_file('run_plain_t.toml', 'toml_plain_t', ok)
      call check(error, .not. ok .and. has_text('toml_plain_t', "'ri-ccsd', 'ri-ccsd(t)'"), &
                 "[correlation] method = ccsd(t) was accepted or the refusal does not name the RI names")
   end subroutine test_toml

   subroutine terco_file(file, name, ok)
      character(len=*), intent(in) :: file, name
      logical, intent(out) :: ok
      integer :: rc
      call run('cd '//work_dir//'/case && '//cli_exe//' '//file//' --json '//name//'.json > '//name//'.out 2>&1', rc)
      ok = rc == 0
   end subroutine terco_file

   subroutine test_convergence_flags(error)
      type(error_type), allocatable, intent(out) :: error
      type(rccsd_refs_t) :: refs
      integer :: n_default, n_loose, n_nodiis, n_two
      real(dp) :: e_loose
      logical :: ok, found

      if (.not. supported()) return
      call setup()
      call read_rccsd_refs(REF, refs)
      call require_run(error, COMMON//' --frozen 1 --riccsd'//TIGHT, 'conv_default', ok)
      if (.not. ok) return
      call json_int('conv_default', "ri-ccsd iterations", n_default, found)
      call require_run(error, COMMON//' --frozen 1 --riccsd --cc-econv 1e-3 --cc-tconv 1e-2', 'conv_loose', ok)
      if (.not. ok) return
      call json_int('conv_loose', "ri-ccsd iterations", n_loose, found)
      call json_real('conv_loose', "ri-ccsd correlation energy", e_loose, found)
      call check(error, n_loose < n_default .and. abs(e_loose - refs%e_ccsd_corr) > 1.0e-6_dp, &
                 "loose --cc-econv/--cc-tconv: "//to_string(n_loose)//" iterations against "//to_string(n_default)// &
                 " and E_CCSD off the converged one by "//to_string(abs(e_loose - refs%e_ccsd_corr), "es9.2"))
      if (allocated(error)) return
      ! the stopping rule is energy AND amplitude: tightening either one alone costs iterations
      call require_run(error, COMMON//' --frozen 1 --riccsd --cc-econv 1e-12 --cc-tconv 1e-2', 'conv_econv', ok)
      if (.not. ok) return
      call json_int('conv_econv', "ri-ccsd iterations", n_two, found)
      call check(error, n_two > n_loose, "--cc-econv 1e-12 alone took "//to_string(n_two)//" iterations, the loose run "// &
                 to_string(n_loose))
      if (allocated(error)) return
      call require_run(error, COMMON//' --frozen 1 --riccsd --cc-econv 1e-3 --cc-tconv 1e-10', 'conv_tconv', ok)
      if (.not. ok) return
      call json_int('conv_tconv', "ri-ccsd iterations", n_two, found)
      call check(error, n_two > n_loose, "--cc-tconv 1e-10 alone took "//to_string(n_two)//" iterations, the loose run "// &
                 to_string(n_loose))
      if (allocated(error)) return
      call require_run(error, COMMON//' --frozen 1 --riccsd'//TIGHT//' --cc-diis 1', 'conv_nodiis', ok)
      if (.not. ok) return
      call json_int('conv_nodiis', "ri-ccsd iterations", n_nodiis, found)
      call check(error, n_nodiis > n_default, "--cc-diis 1 took "//to_string(n_nodiis)//" iterations, the default "// &
                 to_string(n_default))
      if (allocated(error)) return
      ! --cc-iterations is exercised by the non-convergence case; a value that suffices changes nothing
      call require_run(error, COMMON//' --frozen 1 --riccsd'//TIGHT//' --cc-iterations 40', 'conv_forty', ok)
      if (.not. ok) return
      call json_int('conv_forty', "ri-ccsd iterations", n_two, found)
      call check(error, n_two == n_default, "--cc-iterations 40 changed the iteration count: "//to_string(n_two))
   end subroutine test_convergence_flags

   subroutine test_algorithm_flags(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: e_auto, t_auto, e, t
      character(len=16), parameter :: FORMS(2) = [character(len=16) :: "assembled", "direct"]
      logical :: ok, found
      integer :: k

      if (.not. supported()) return
      call setup()
      call require_run(error, COMMON//' --frozen 1 --riccsd-t'//TIGHT, 'alg_auto', ok)
      if (.not. ok) return
      call json_real('alg_auto', "ri-ccsd correlation energy", e_auto, found)
      call json_real('alg_auto', "ri-ccsd triples correction", t_auto, found)
      do k = 1, 2
         call require_run(error, COMMON//' --frozen 1 --riccsd-t'//TIGHT//' --cc-ladder '//trim(FORMS(k)), 'alg_'// &
                          trim(FORMS(k)), ok)
         if (.not. ok) return
         call json_real('alg_'//trim(FORMS(k)), "ri-ccsd correlation energy", e, found)
         call check(error, found .and. abs(e - e_auto) < 1.0e-10_dp, "--cc-ladder "//trim(FORMS(k))//" moved E_CCSD by "// &
                    to_string(e - e_auto, "es9.2"))
         if (allocated(error)) return
      end do
      call require_run(error, COMMON//' --frozen 1 --riccsd-t'//TIGHT//' --cc-triples-budget 0.001', 'alg_budget', ok)
      if (.not. ok) return
      call json_real('alg_budget', "ri-ccsd triples correction", t, found)
      call check(error, found .and. abs(t - t_auto) < TOL_SAME, "a 1 MiB --cc-triples-budget moved E(T) by "// &
                 to_string(t - t_auto, "es9.2"))
   end subroutine test_algorithm_flags

   subroutine test_noconv(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: text
      logical :: converged, found, f2
      real(dp) :: e, v
      integer :: rc, iterations

      if (.not. supported()) return
      call setup()
      call terco(COMMON//' --frozen 1 --riccsd-t --cc-iterations 2', 'noconv', rc)
      call check(error, rc == 2, "a CCSD that does not converge exited "//to_string(rc)//", not 2 (case/noconv.out)")
      if (allocated(error)) return
      text = out_text('noconv')
      call check(error, index(text, "RI-CCSD DID NOT CONVERGE") > 0, "the energy box does not say the CCSD did not converge")
      if (allocated(error)) return
      call check(error, index(text, "NOT CONVERGED") > 0, "the SUMMARY block does not say NOT CONVERGED")
      if (allocated(error)) return
      call check(error, index(text, "not converged in 2 iterations") > 0, "no error message with the count (the library's)")
      if (allocated(error)) return
      call check(error, index(text, "TOTAL RI-MP2 ENERGY") > 0, "the RI-MP2 total is missing although the MP2 part is finished")
      if (allocated(error)) return
      call check(error, index(text, "TOTAL RI-CCSD ENERGY") == 0 .and. index(text, "(T) CORRECTION") == 0, &
                 "a total or a (T) correction is printed for an unconverged CCSD")
      if (allocated(error)) return
      call json_logical('noconv', "ri-ccsd converged", converged, found)
      call check(error, found .and. .not. converged, "the JSON does not say ccsd converged: false")
      if (allocated(error)) return
      call json_int('noconv', "ri-ccsd iterations", iterations, found)
      call check(error, found .and. iterations == 2, "the JSON ccsd iterations is "//to_string(iterations))
      if (allocated(error)) return
      call json_real('noconv', "ri-ccsd total energy", v, found)
      call json_real('noconv', "ri-ccsd triples correction", e, f2)
      call check(error, .not. found .and. .not. f2, "the JSON carries a total energy or a (T) correction of a failed CCSD")
      if (allocated(error)) return
      call json_real('noconv', "ri-ccsd correlation energy", e, found)
      call check(error, found .and. e < 0.0_dp, "the last iterate is not in the JSON")
   end subroutine test_noconv

   subroutine check_t_failed(error, name, label)
      !! case/<name>.out and .json are those of a converged RI-CCSD whose (T) refused: the RI-CCSD stands, no (T) anywhere
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: name, label
      character(len=:), allocatable :: reason
      real(dp) :: v, e_json
      logical :: found, ok, f2

      if (allocated(error)) return
      call compare_to_reference(error, name, .false., label)
      if (allocated(error)) return
      call json_logical(name, "ri-ccsd converged", ok, found)
      call check(error, found .and. ok, label//': the JSON does not say "ri-ccsd converged": true')
      if (allocated(error)) return
      call json_real(name, "ri-ccsd triples correction", v, found)
      call json_real(name, "ri-ccsd triples total energy", v, f2)
      call check(error, .not. found .and. .not. f2, label//": the JSON carries a (T) correction or total of a failed (T)")
      if (allocated(error)) return
      call json_string(name, "ri-ccsd triples error", reason, found)
      call check(error, found, label//': no "ri-ccsd triples error" in the JSON')
      if (allocated(error)) return
      call check(error, index(reason, "budget") > 0 .and. index(reason, "(T)") > 0, label//": the JSON's reason is '"// &
                 reason//"'")
      if (allocated(error)) return
      ! the energy box and the summary keep the RI-CCSD, and the box says what happened
      call json_real(name, "ri-ccsd correlation energy", e_json, found)
      call summary_value(name, "RI-CCSD CORRELATION", v, ok)
      call check(error, ok .and. found .and. abs(v - e_json) < PRINT_ROUNDING + 1.0e-14_dp, label// &
                 ": the box's RI-CCSD CORRELATION "//to_string(v, "f20.12")//" is not the JSON's "//to_string(e_json, "f20.12"))
      if (allocated(error)) return
      call json_real(name, "ri-ccsd total energy", e_json, found)
      call summary_value(name, "TOTAL RI-CCSD ENERGY", v, ok)
      call check(error, ok .and. found .and. abs(v - e_json) < PRINT_ROUNDING + 1.0e-14_dp, label// &
                 ": the box's TOTAL RI-CCSD ENERGY "//to_string(v, "f20.12")//" is not the JSON's "//to_string(e_json, "f20.12"))
      if (allocated(error)) return
      call check(error, has_text(name, "(T) FAILED"), label//': no "(T) FAILED" line in case/'//name//".out")
      if (allocated(error)) return
      call check(error, has_text(name, "ri-ccsd correlation") .and. has_text(name, "ri-ccsd total energy") .and. &
                 has_text(name, "ri-ccsd (t) error"), label//": the SUMMARY lacks the RI-CCSD rows or the (T) error row")
      if (allocated(error)) return
      call check(error, .not. has_text(name, "(T) CORRECTION") .and. .not. has_text(name, "TOTAL RI-CCSD(T) ENERGY") .and. &
                 .not. has_text(name, "ri-ccsd (t) correction") .and. .not. has_text(name, "ri-ccsd(t) total"), &
                 label//": a (T) correction or an RI-CCSD(T) total is printed for a failed (T)")
      if (allocated(error)) return
      call check(error, has_text(name, "does not hold one triple"), label//": the reason is not on the error stream")
   end subroutine check_t_failed

   subroutine test_triples_failed(error)
      !! 1e-6 GiB is 1073 bytes; one water triple needs 48*4^3 = 3072 bytes of W
      type(error_type), allocatable, intent(out) :: error
      character(len=64), allocatable :: lines(:)
      character(len=:), allocatable :: text
      integer :: rc
      logical :: ok

      if (.not. supported()) return
      call setup()
      call terco(COMMON//' --frozen 1 --riccsd-t --cc-triples tiled --cc-triples-budget 1e-6'//TIGHT, 't_flags', rc)
      call check(error, rc == 1, "a refused (T) exited "//to_string(rc)//", not 1 (case/t_flags.out)")
      call check_t_failed(error, 't_flags', "flags")
      if (allocated(error)) return
      lines = [character(len=64) :: '[molecule]', 'geometry = "water.xyz"', '[basis]', &
               'main = "basis_sets/cc-pvdz.json"', 'auxiliary = "basis_sets/cc-pvdz-rifit.json"', &
               '[scf]', 'energy_tol = 1e-13', 'commutator_tol = 1e-11', &
               '[correlation]', 'method = "ri-ccsd(t)"', 'frozen_core = 1', &
               '[riccsd]', 'energy_tol = 1e-12', 'amplitude_tol = 1e-10', 'triples = "tiled"', 'triples_budget = 1e-6']
      call write_toml(work_dir//'/case/t_run.toml', lines)
      call run('cd '//work_dir//'/case && '//cli_exe//' t_run.toml --json t_toml.json > t_toml.out 2>&1', rc)
      call check(error, rc == 1, "a refused (T) from a TOML file exited "//to_string(rc)//", not 1 (case/t_toml.out)")
      call check_t_failed(error, 't_toml', "TOML")
      if (allocated(error)) return
      call json_string('t_toml', "input.riccsd.triples", text, ok)
      call check(error, ok, 'the JSON echo lacks input.riccsd.triples')
      if (allocated(error)) return
      call check(error, text == "tiled", "the echo says triples '"//text//"'")
      if (allocated(error)) return
      ! the reference routine has no budget: the same tiny budget is harmless, and (T) is pyscf's
      call require_run(error, COMMON//' --frozen 1 --riccsd-t --cc-triples reference --cc-triples-budget 1e-6'//TIGHT, &
                       't_reference', ok)
      if (.not. ok) return
      call compare_to_reference(error, 't_reference', .true., "--cc-triples reference")
      if (allocated(error)) return
      call require_run(error, COMMON//' --frozen 1 --riccsd-t --cc-triples tiled'//TIGHT, 't_tiled', ok)
      if (.not. ok) return
      call compare_to_reference(error, 't_tiled', .true., "--cc-triples tiled")
   end subroutine test_triples_failed

   subroutine test_silent(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: ok

      if (.not. supported()) return
      call setup()
      call require_run(error, COMMON//' --frozen 1 --riccsd-t'//TIGHT//' --silent', 'silent', ok)
      if (.not. ok) return
      call check(error, has_text('silent', "TOTAL RI-CCSD(T) ENERGY"), "--silent lost the summary")
      if (allocated(error)) return
      call check(error, .not. has_text('silent', "E_CCSD(corr)") .and. .not. has_text('silent', "iter  ") .and. &
                 .not. has_text('silent', "rccsd: step"), "--silent printed the iteration table or the stage timings")
      if (allocated(error)) return
      call require_run(error, COMMON//' --frozen 1 --riccsd-t'//TIGHT, 'verbose', ok)
      if (.not. ok) return
      call check(error, has_text('verbose', "E_CCSD(corr)") .and. has_text('verbose', "rccsd: step"), &
                 "the default run lacks the iteration table or the stage timings")
   end subroutine test_silent

   subroutine refused(error, args, name, reason)
      !! `terco args` exits non-zero and says `reason`
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: args, name, reason
      integer :: rc
      if (allocated(error)) return
      call terco(args, name, rc)
      call check(error, rc /= 0, "terco"//args//" was not refused (exit 0, case/"//name//".out)")
      if (allocated(error)) return
      call check(error, has_text(name, reason), "terco"//args//": the refusal does not contain '"//reason//"' (case/"// &
                 name//".out)")
      if (allocated(error)) return
      ! A message handed over with a wrong length prints a null and stray text after it (nvfortran 26.9 did this
      ! to one refusal, and crashed on another machine).
      call check(error, .not. has_text(name, achar(0)), "terco"//args//": the output holds a null character (case/"// &
                 name//".out): a message was passed with a wrong length")
   end subroutine refused

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error

      if (.not. supported()) return
      call setup()
      call refused(error, COMMON//' --riccsd --rimp2', 'r_rimp2', "cannot be combined")
      call refused(error, COMMON//' --riccsd -f pbe', 'r_func', "Hartree-Fock reference")
      call refused(error, COMMON//' --riccsd-t --uhf 2', 'r_uhf', "--uhf 2 asks for unpaired electrons")
      call refused(error, COMMON//' --riccsd-t -c 1', 'r_charge', "check --charge")
      call refused(error, COMMON//' --riccsd --grad', 'r_grad', "has no gradient")
      call refused(error, COMMON//' --riccsd --hessian', 'r_hessian', "has no gradient")
      call refused(error, COMMON//' --riccsd --opt', 'r_opt', "has no gradient")
      call refused(error, ' hi.xyz -b def2-svp --riccsd', 'r_ecp', "effective core potential")
      call refused(error, COMMON//' --cc-econv 1e-8', 'r_nocc', "need --riccsd")
      call refused(error, COMMON//' --cc-ladder direct', 'r_nocc2', "need --riccsd")
      call refused(error, COMMON//' --riccsd --cc-triples-budget 1', 'r_budget', "needs --riccsd-t")
      call refused(error, COMMON//' --riccsd --cc-triples tiled', 'r_triples', "needs --riccsd-t")
      call refused(error, COMMON//' --cc-triples tiled', 'r_triples2', "need --riccsd")
   end subroutine test_refusals

   subroutine test_bad_values(error)
      type(error_type), allocatable, intent(out) :: error

      if (.not. supported()) return
      call setup()
      call refused(error, COMMON//' --riccsd --cc-econv 0', 'b_econv', "--cc-econv ([riccsd] energy_tol) must be positive")
      call refused(error, COMMON//' --riccsd --cc-econv -1e-8', 'b_econv2', "[riccsd] energy_tol")
      call refused(error, COMMON//' --riccsd --cc-tconv 0', 'b_tconv', "--cc-tconv ([riccsd] amplitude_tol) must be positive")
      call refused(error, COMMON//' --riccsd --cc-iterations 0', 'b_iter', &
                   "--cc-iterations ([riccsd] max_iterations) must be at least 1")
      call refused(error, COMMON//' --riccsd --cc-diis 0', 'b_diis', "--cc-diis ([riccsd] diis_size) must be between 1")
      call refused(error, COMMON//' --riccsd --cc-diis 51', 'b_diis2', "and 50")
      call refused(error, COMMON//' --riccsd --cc-ladder fast', 'b_ladder', &
                   "--cc-ladder ([riccsd] ladder) must be auto, assembled or direct")
      call refused(error, COMMON//' --riccsd-t --cc-triples fast', 'b_triples', &
                   "--cc-triples ([riccsd] triples) must be auto, reference or tiled")
      call refused(error, COMMON//' --riccsd-t --cc-triples-budget 0', 'b_budget', "must be positive and at most 1e6 GiB")
      call refused(error, COMMON//' --riccsd-t --cc-triples-budget -1', 'b_budget2', "must be positive and at most 1e6 GiB")
      call refused(error, COMMON//' --riccsd-t --cc-triples-budget 2e6', 'b_budget3', "must be positive and at most 1e6 GiB")
      call refused(error, COMMON//' --riccsd --frozen 5', 'b_frozen', "leaves no correlated orbital")
   end subroutine test_bad_values

end module check_cli_rccsd_suite

program check_cli_rccsd
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_rccsd_suite, only: collect_cli_rccsd, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_cli_rccsd: usage: check_cli_rccsd <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_cli_rccsd, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_rccsd
