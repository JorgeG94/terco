module check_cli_tddft_suite
   !!
   !! The CLI's `--tddft`: the whole path a user runs (app/main.F90 -> trc_set_tddft, trc_run_scf, trc_run_tddft,
   !! trc_run_gradient), read back from its --json summary:
   !!  1. TDA excitation energies, oscillator strengths and the total energy of the target root of water PBE/6-31G*
   !!     against pyscf's tdscf (tddft_water_6-31gs_pbe.ref), and the same run in eV;
   !!  2. RPA energies, oscillator strengths and ||Y||^2 against pyscf's TDDFT (rpa_water_6-31gs_pbe.ref);
   !!  3. the same settings given in a TOML file ([tddft] table, Hartree-Fock: CIS against tddft_water_6-31gs_hf.ref);
   !!  4. --grad: the gradient of root 2 of water PBE/6-31G with --no-tddft-grid-response against pyscf's
   !!     grad.tdrks (tddft_grad_water_6-31g_pbe.ref, Pulay only), and the default run (grid response) differing
   !!     from it by the grid's own derivative, about 1e-5, with the JSON saying which energy it differentiates;
   !!  5. what is refused before the SCF, with the reason: a triplet, a meta-GGA, a range-separated hybrid, a
   !!     double hybrid, an open shell, --hessian, a root out of range, a --tddft-* option without --tddft,
   !!     contradicting grid response flags and a spherical basis.
   !! Arguments: the terco executable, a scratch path for its --json output and optionally the mode (see `mode_arg`).
   !! Run from test/.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_vector, ref_read_matrix3
   implicit none
   private
   public :: collect_cli_tddft_tests, set_cli_tddft_args, mode_arg

   character(len=16) :: mode_arg = ""
   !! "grad": the Pulay-only gradient against pyscf; "grid": the default (grid response) gradient; "refusals": the
   !! refusals; empty: the energies and the TOML input. Four runs of one program keep each around ten seconds.

   real(dp), parameter :: TOL_OMEGA = 1.0e-6_dp
   real(dp), parameter :: TOL_OSC = 1.0e-5_dp
   real(dp), parameter :: TOL_YY = 1.0e-6_dp
   !! The library checks' gates against pyscf (check_tddft, check_tddft_rpa).
   real(dp), parameter :: TOL_GRAD = 2.0e-8_dp
   !! The CLI's own SCF stops at --gconv 1e-9, so the gradient is held to the CLI checks' gate (check_cli_grad_rimp2).
   real(dp), parameter :: HARTREE_TO_EV = 27.211386245988_dp

   character(len=:), allocatable :: g_cli_exe, g_probe

contains

   subroutine set_cli_tddft_args(cli_exe, probe)
      character(*), intent(in) :: cli_exe, probe
      g_cli_exe = cli_exe
      g_probe = probe
   end subroutine set_cli_tddft_args

   subroutine collect_cli_tddft_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (trim(mode_arg) == "grad") then
         testsuite = [new_unittest("CLI --tddft --grad --no-tddft-grid-response: the gradient matches pyscf", test_gradient)]
      else if (trim(mode_arg) == "grid") then
         testsuite = [new_unittest("CLI --tddft --grad: the grid response is the default and moves the gradient", &
                                   test_grid_response)]
      else if (trim(mode_arg) == "refusals") then
         testsuite = [new_unittest("CLI --tddft refusals, before the SCF", test_refusals)]
      else
         testsuite = [ &
                     new_unittest("CLI --tddft TDA energies match pyscf, in Eh and eV", test_tda), &
                     new_unittest("CLI --tddft --tddft-approx rpa: energies and ||Y||^2 match pyscf", test_rpa), &
                     new_unittest("CLI [tddft] table of a TOML file: CIS matches pyscf", test_toml), &
                     new_unittest("CLI --tddft-missed-root-iterations and its TOML key reach the solver", test_missed_root), &
                     new_unittest("CLI --no-tddft-density-screen and its TOML key reach the solver", test_density_screen) &
                     ]
      end if
   end subroutine collect_cli_tddft_tests

   subroutine run_cli(error, args, probe, json)
      !! The CLI with `args` plus --json `probe` must succeed; its summary is read into `json`.
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: args, probe
      type(json_file), intent(out) :: json
      integer :: rc

      call execute_command_line(g_cli_exe//' '//args//' --econv 1e-12 --gconv 1e-9 --json '//probe//' > '//probe//'.out 2>&1', &
                                exitstat=rc, wait=.true.)
      call check(error, rc == 0, "the CLI run failed, exit code "//to_string(rc)//" (see "//probe//".out)")
      if (allocated(error)) return
      call json%initialize()
      call json%load(filename=probe)
      call check(error, .not. json%failed(), "cannot read "//probe)
   end subroutine run_cli

   subroutine read_roots(error, json, nr, omega, osc, ynorm)
      type(error_type), allocatable, intent(out) :: error
      type(json_file), intent(inout) :: json
      integer, intent(in) :: nr
      real(dp), allocatable, intent(out) :: omega(:), osc(:), ynorm(:)
      logical :: found_omega, found_osc, found_y

      call json%get("tddft.excitation energies", omega, found_omega)
      call json%get("tddft.oscillator strengths", osc, found_osc)
      call json%get("tddft.y norms", ynorm, found_y)
      call check(error, found_omega .and. found_osc .and. found_y, "the JSON has no tddft block")
      if (allocated(error)) return
      call check(error, size(omega) == nr .and. size(osc) == nr .and. size(ynorm) == nr, &
                 "the JSON holds "//to_string(size(omega))//" roots, "//to_string(nr)//" were asked for")
   end subroutine read_roots

   subroutine compare_with_ref(error, ref_name, omega, osc, ynorm, with_y)
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: ref_name
      real(dp), intent(in) :: omega(:), osc(:), ynorm(:)
      logical, intent(in) :: with_y
      real(dp), allocatable :: om_ref(:), osc_ref(:), yy_ref(:)
      real(dp) :: e_ref
      integer :: u, nr

      nr = size(omega)
      allocate (om_ref(nr), osc_ref(nr), yy_ref(nr))
      call ref_require(ref_name, u)
      call ref_read_real(u, "E_SCF", e_ref)
      call ref_read_vector(u, "OMEGA", nr, om_ref)
      call ref_read_vector(u, "OSC", nr, osc_ref)
      if (with_y) call ref_read_vector(u, "YY", nr, yy_ref)
      call ref_close(u)
      call check(error, maxval(abs(omega - om_ref)) < TOL_OMEGA, "omega vs pyscf, "// &
                 to_string(maxval(abs(omega - om_ref)), "es10.2"))
      if (allocated(error)) return
      call check(error, maxval(abs(osc - osc_ref)) < TOL_OSC, "oscillator strengths vs pyscf, "// &
                 to_string(maxval(abs(osc - osc_ref)), "es10.2"))
      if (allocated(error)) return
      if (with_y) then
         call check(error, maxval(abs(ynorm - yy_ref)) < TOL_YY, "||Y||^2 vs pyscf, "// &
                    to_string(maxval(abs(ynorm - yy_ref)), "es10.2"))
      else
         call check(error, maxval(abs(ynorm)) == 0.0_dp, "the TDA has Y = 0")
      end if
   end subroutine compare_with_ref

   subroutine test_tda(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp), allocatable :: omega(:), osc(:), ynorm(:), ev(:)
      real(dp) :: etot, escf
      logical :: found
      integer :: root

      call run_cli(error, 'water.xyz -b basis_sets/6-31gs.json -f pbe --tddft --tddft-roots 4 --tddft-root 3', g_probe, json)
      if (allocated(error)) return
      call read_roots(error, json, 4, omega, osc, ynorm)
      if (.not. allocated(error)) call compare_with_ref(error, "tddft_water_6-31gs_pbe.ref", omega, osc, ynorm, .false.)
      if (.not. allocated(error)) then
         call json%get("tddft.excitation energies ev", ev, found)
         call check(error, found, "no excitation energies in eV")
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(ev - omega*HARTREE_TO_EV)) < 1.0e-9_dp, "eV is Eh times 27.2113862")
      end if
      if (.not. allocated(error)) then
         call json%get("tddft.target root", root, found)
         call json%get("tddft.target root total energy", etot, found)
         call json%get("scf energy", escf, found)
         call check(error, root == 3, "the target root is the one asked for")
         if (.not. allocated(error)) call check(error, abs(etot - (escf + omega(3))) < 1.0e-12_dp, &
                                                "the target root's energy is E_SCF + omega")
      end if
      call json%destroy()
   end subroutine test_tda

   subroutine test_rpa(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp), allocatable :: omega(:), osc(:), ynorm(:)

      call run_cli(error, 'water.xyz -b basis_sets/6-31gs.json -f pbe --tddft --tddft-roots 4 --tddft-approx rpa', &
                   g_probe//'.rpa.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 4, omega, osc, ynorm)
      if (.not. allocated(error)) call compare_with_ref(error, "rpa_water_6-31gs_pbe.ref", omega, osc, ynorm, .true.)
      call json%destroy()
   end subroutine test_rpa

   subroutine test_toml(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp), allocatable :: omega(:), osc(:), ynorm(:)

      ! test/tddft_cis_water.toml: its paths are relative to the file.
      call run_cli(error, 'tddft_cis_water.toml', g_probe//'.toml.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 4, omega, osc, ynorm)
      if (.not. allocated(error)) call compare_with_ref(error, "tddft_water_6-31gs_hf.ref", omega, osc, ynorm, .false.)
      call json%destroy()
   end subroutine test_toml

   subroutine test_missed_root(error)
      !! Formaldehyde 6-31+G Hartree-Fock, three roots: the third is reached only by iterating the candidate pairs
      !! after the probe. The default is a budget of 0, which loses it (a different third root, 0.353367 Eh); 24, given
      !! as a flag or as the TOML key, recovers it (0.348172 Eh).
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: BASE = 'formaldehyde.xyz -b ../basis_sets/6-31+g.json --tddft --tddft-roots 3'
      type(json_file) :: json
      real(dp), allocatable :: omega(:), osc(:), ynorm(:)
      real(dp) :: w_default, w_0, w_flag, w_toml, w_grad, w_grad0
      integer :: n_default, n_0, n_flag, n_toml, n_grad

      call run_cli(error, BASE, g_probe//'.default.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 3, omega, osc, ynorm)
      if (allocated(error)) return
      w_default = omega(3)
      n_default = open_count(json)
      call json%destroy()
      call run_cli(error, BASE//' --tddft-missed-root-iterations 0', g_probe//'.m0.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 3, omega, osc, ynorm)
      if (allocated(error)) return
      w_0 = omega(3)
      n_0 = open_count(json)
      call json%destroy()
      call run_cli(error, BASE//' --tddft-missed-root-iterations 24', g_probe//'.m24.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 3, omega, osc, ynorm)
      if (allocated(error)) return
      w_flag = omega(3)
      n_flag = open_count(json)
      call json%destroy()
      call run_cli(error, 'tddft_missed_root_formaldehyde.toml', g_probe//'.m24toml.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 3, omega, osc, ynorm)
      if (allocated(error)) return
      w_toml = omega(3)
      n_toml = open_count(json)
      call json%destroy()
      ! Separate runs: device atomics give run-to-run noise near 1e-12, far below the 5e-3 a lost root moves.
      call check(error, abs(w_0 - w_default) < 1.0e-9_dp, "--tddft-missed-root-iterations 0 is not the default")
      if (allocated(error)) return
      call check(error, abs(w_flag - w_default) > 1.0e-3_dp, "--tddft-missed-root-iterations 24 left the third root alone")
      if (allocated(error)) return
      call check(error, abs(w_flag - 0.348172_dp) < 2.0e-6_dp, "--tddft-missed-root-iterations 24 did not recover the root")
      if (allocated(error)) return
      call check(error, abs(w_toml - w_flag) < 1.0e-9_dp, "the [tddft] missed_root_iterations key and the flag disagree")
      if (allocated(error)) return
      ! The JSON's "candidates open": the missed-root check's open Ritz pairs, with the roots that lost one.
      call check(error, n_default > 0 .and. n_0 == n_default, "the default (budget 0) left no candidate open in the JSON")
      if (allocated(error)) return
      call check(error, n_flag == 0 .and. n_toml == 0, "a budget of 24 left candidates open: flag "// &
                 to_string(n_flag)//", TOML "//to_string(n_toml))
      if (allocated(error)) return
      ! --grad follows a root by index, so it keeps the guard (24) unless a budget is given: the third root's gradient
      ! must not be that of the wrong state. An explicit 0 is honoured.
      call run_cli(error, BASE//' --tddft-root 3 --grad', g_probe//'.grad.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 3, omega, osc, ynorm)
      if (allocated(error)) return
      w_grad = omega(3)
      n_grad = open_count(json)
      call json%destroy()
      call run_cli(error, BASE//' --tddft-root 3 --grad --tddft-missed-root-iterations 0', g_probe//'.grad0.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 3, omega, osc, ynorm)
      if (allocated(error)) return
      w_grad0 = omega(3)
      call json%destroy()
      call check(error, abs(w_grad - w_flag) < 1.0e-9_dp .and. n_grad == 0, &
                 "--grad did not keep the missed-root guard: third root "//to_string(w_grad, "f10.6")// &
                 ", candidates open "//to_string(n_grad))
      if (allocated(error)) return
      call check(error, abs(w_grad0 - w_default) < 1.0e-9_dp, "--grad with an explicit budget of 0 did not honour it")
   end subroutine test_missed_root

   integer function open_count(json) result(n)
      !! The `candidates open` count of the run's `tddft` block (-1 when the key is absent).
      type(json_file), intent(inout) :: json
      logical :: found
      call json%get("tddft.candidates open", n, found)
      if (.not. found) n = -1
   end function open_count

   function warned(file, words) result(n)
      !! How many lines of `file` contain `words`.
      character(*), intent(in) :: file, words
      integer :: n
      character(len=1024) :: line
      integer :: u, ios
      n = 0
      open (newunit=u, file=file, status="old", action="read", iostat=ios)
      do while (ios == 0)
         read (u, "(a)", iostat=ios) line
         if (ios == 0 .and. index(line, words) > 0) n = n + 1
      end do
      close (u)
   end function warned

   subroutine test_density_screen(error)
      !! Water STO-3G Hartree-Fock, two roots: the default run screens the matvec and says nothing; --no-tddft-density-screen
      !! and the [tddft] density_screen = false key each switch it off, which the run logs ONCE (that is what shows the
      !! setting reached the solver: the roots do not change measurably). The roots agree to twice the Davidson tolerance
      !! (--tddft-tol 1e-8 by default): the two operators differ by quartets below the integral threshold. Given twice in
      !! opposite senses, or without --tddft, the flag is refused.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: BASE = 'water.xyz -b basis_sets/sto-3g.json --tddft --tddft-roots 2'
      character(len=*), parameter :: WORDS = 'density screen of the excited-state matvec is OFF'
      type(json_file) :: json
      real(dp), allocatable :: omega(:), osc(:), ynorm(:)
      real(dp) :: w_default(2), w_flag(2), w_toml(2)

      call run_cli(error, BASE, g_probe//'.screen_default.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 2, omega, osc, ynorm)
      w_default = omega
      call json%destroy()
      call run_cli(error, BASE//' --no-tddft-density-screen', g_probe//'.screen_flag.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 2, omega, osc, ynorm)
      w_flag = omega
      call json%destroy()
      call run_cli(error, 'tddft_density_screen_off.toml', g_probe//'.screen_toml.json', json)
      if (allocated(error)) return
      call read_roots(error, json, 2, omega, osc, ynorm)
      w_toml = omega
      call json%destroy()
      call check(error, warned(g_probe//'.screen_default.json.out', WORDS) == 0, "the default run warned about the screen")
      if (allocated(error)) return
      call check(error, warned(g_probe//'.screen_flag.json.out', WORDS) == 1, &
                 "--no-tddft-density-screen did not log its warning once")
      if (allocated(error)) return
      call check(error, warned(g_probe//'.screen_toml.json.out', WORDS) == 1, &
                 "[tddft] density_screen = false did not log its warning once")
      if (allocated(error)) return
      call check(error, maxval(abs(w_flag - w_default)) < 2.0e-8_dp, "the unscreened roots differ from the default's by "// &
                 to_string(maxval(abs(w_flag - w_default)), "es10.3"))
      if (allocated(error)) return
      call check(error, maxval(abs(w_toml - w_flag)) < 2.0e-8_dp, "the TOML key and the flag disagree")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' --tddft-density-screen --no-tddft-density-screen', '.screen_both.out', "contradict")
      if (allocated(error)) return
      call expect_refusal(error, 'water.xyz -b basis_sets/sto-3g.json --no-tddft-density-screen', '.screen_alone.out', &
                          "need --tddft")
   end subroutine test_density_screen

   subroutine read_gradient(error, json, grad, what)
      type(error_type), allocatable, intent(out) :: error
      type(json_file), intent(inout) :: json
      real(dp), intent(out) :: grad(3, 3)
      character(len=:), allocatable, intent(out) :: what
      real(dp), allocatable :: row(:)
      logical :: found
      integer :: i

      grad = 0.0_dp
      call json%get("gradient of", what, found)
      if (.not. found) what = ""
      do i = 1, 3
         call json%get("gradient("//to_string(i)//")", row, found)
         call check(error, found, "no gradient row "//to_string(i)//" in the JSON")
         if (allocated(error)) return
         grad(:, i) = row
      end do
   end subroutine read_gradient

   subroutine test_gradient(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: pulay(3, 3), gref(3, 3, 3), e_ref
      character(len=:), allocatable :: what
      logical :: found, response
      integer :: u

      call run_cli(error, 'water.xyz -b basis_sets/6-31g.json -f pbe --tddft --tddft-roots 3 --tddft-root 2 '// &
                   '--no-tddft-grid-response --tddft-tol 1e-10 --tddft-cpks-tol 1e-11 --grad', g_probe//'.pulay.json', json)
      if (allocated(error)) return
      call read_gradient(error, json, pulay, what)
      if (.not. allocated(error)) then
         call json%get("tddft.grid response", response, found)
         call check(error, found .and. .not. response, "the JSON says the grid response was off")
      end if
      if (.not. allocated(error)) then
         call check(error, index(what, "tddft root 2") > 0, "the JSON does not say the gradient is of root 2: '"//what//"'")
      end if
      call json%destroy()
      if (allocated(error)) return
      call ref_require("tddft_grad_water_6-31g_pbe.ref", u)
      call ref_read_real(u, "E_SCF", e_ref)
      call ref_read_matrix3(u, "GRAD", 3, 3, 3, gref)
      call ref_close(u)
      call check(error, maxval(abs(pulay - gref(:, :, 2))) < TOL_GRAD, "Pulay-only gradient of root 2 vs pyscf, "// &
                 to_string(maxval(abs(pulay - gref(:, :, 2))), "es10.2"))
   end subroutine test_gradient

   subroutine test_grid_response(error)
      !! The default gradient has the grid response and differs from pyscf's Pulay-only one by the grid's derivative
      !! (about 1e-5): the reference is that Pulay-only gradient, which test_gradient shows the CLI reproduces with
      !! --no-tddft-grid-response to 2e-8.
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: exact(3, 3), gref(3, 3, 3), e_ref, diff
      character(len=:), allocatable :: what
      logical :: found, response
      integer :: u

      call run_cli(error, 'water.xyz -b basis_sets/6-31g.json -f pbe --tddft --tddft-roots 3 --tddft-root 2 '// &
                   '--tddft-tol 1e-10 --tddft-cpks-tol 1e-11 --grad', g_probe//'.exact.json', json)
      if (allocated(error)) return
      call read_gradient(error, json, exact, what)
      call json%get("tddft.grid response", response, found)
      call json%destroy()
      call check(error, found .and. response, "the grid response is on by default")
      if (allocated(error)) return
      call ref_require("tddft_grad_water_6-31g_pbe.ref", u)
      call ref_read_real(u, "E_SCF", e_ref)
      call ref_read_matrix3(u, "GRAD", 3, 3, 3, gref)
      call ref_close(u)
      diff = maxval(abs(exact - gref(:, :, 2)))
      call check(error, diff > 1.0e-6_dp .and. diff < 1.0e-4_dp, "the default gradient should differ from the "// &
                 "Pulay-only one by the grid's derivative, about 1e-5; got "//to_string(diff, "es10.2"))
   end subroutine test_grid_response

   subroutine expect_refusal(error, args, suffix, reason)
      !! The CLI with `args` must exit non-zero and print `reason`.
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: args, suffix, reason
      character(len=:), allocatable :: out
      character(len=512) :: line
      integer :: rc, u, ios
      logical :: seen

      out = g_probe//suffix
      call execute_command_line(g_cli_exe//' '//args//' > '//out//' 2>&1', exitstat=rc, wait=.true.)
      call check(error, rc /= 0, "not refused (exit code 0, see "//out//")")
      if (allocated(error)) return
      seen = .false.
      open (newunit=u, file=out, status="old", action="read", iostat=ios)
      do while (ios == 0)
         read (u, "(a)", iostat=ios) line
         if (ios == 0 .and. index(line, reason) > 0) seen = .true.
      end do
      close (u)
      call check(error, seen, "the refusal does not say '"//reason//"' (see "//out//")")
   end subroutine expect_refusal

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: BASE = 'water.xyz -b basis_sets/sto-3g.json --tddft'

      call expect_refusal(error, BASE//' --tddft-spin triplet', '.triplet.out', "triplet excitations are not implemented")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' -f tpss', '.mgga.out', "meta-GGA")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' -f wb97x', '.rsh.out', "range-separated")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' -f wb97x-v', '.vv10.out', "range-separated")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' -f b2plyp', '.dh.out', "double-hybrid")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' --uhf 2', '.uhf.out', "closed-shell only")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' --hessian', '.hessian.out', "no Hessian")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' --df', '.df.out', "--df")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' --tddft-roots 2 --tddft-root 3', '.root.out', "must not exceed")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' --tddft-approx cis', '.approx.out', "must be tda")
      if (allocated(error)) return
      call expect_refusal(error, 'water.xyz -b basis_sets/sto-3g.json --tddft-roots 2', '.alone.out', "need --tddft")
      if (allocated(error)) return
      call expect_refusal(error, BASE//' -f pbe --tddft-grid-response --no-tddft-grid-response', '.flags.out', "contradict")
      if (allocated(error)) return
      call expect_refusal(error, 'water.xyz -b basis_sets/cc-pvdz.json --tddft', '.spherical.out', "spherical basis")
      if (allocated(error)) return
      ! The gap window may be zero (the warning off) but not negative; the C ABI (trc_opt_set_tddft) agrees.
      call expect_refusal(error, BASE//' --tddft-gap-warn -0.001', '.gap.out', "zero or more")
      if (allocated(error)) return
      call expect_refusal(error, 'water.xyz -b basis_sets/sto-3g.json --tddft-strict-root', '.strict.out', "need --tddft")
      if (allocated(error)) return
      ! The missed-root budget may be zero but not negative; given without --tddft it is refused like its siblings.
      call expect_refusal(error, BASE//' --tddft-missed-root-iterations -1', '.missed.out', "zero or more")
      if (allocated(error)) return
      call expect_refusal(error, 'water.xyz -b basis_sets/sto-3g.json --tddft-missed-root-iterations 0', '.missedalone.out', &
                          "need --tddft")
   end subroutine test_refusals

end module check_cli_tddft_suite

program check_cli_tddft
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_tddft_suite, only: collect_cli_tddft_tests, set_cli_tddft_args, mode_arg
   implicit none
   character(len=1024) :: cli_exe, probe
   integer :: stat

   call get_command_argument(1, cli_exe)
   call get_command_argument(2, probe)
   call get_command_argument(3, mode_arg)
   call set_cli_tddft_args(trim(cli_exe), trim(probe))
   stat = 0
   call run_testsuite(collect_cli_tddft_tests, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_tddft
