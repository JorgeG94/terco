module check_tddft_df_suite
   !!
   !! TDA and full-RPA singlet excitation energies of a density-fitted (--df) RHF/RKS reference through
   !! `trc_tddft_run_solve`, the excitation operator applied with the same fit the SCF converged with
   !! (trc_tddft's THE FITTED ROUTE, trc_tddft_df.F90):
   !!  1. TDA energies and oscillator strengths of the lowest three singlets against pyscf's TDA on
   !!     `mf.density_fit(auxbasis=def2-universal-jkfit)` (test/pyscf_refs/generate.py, `tddft_df_*` cases);
   !!  2. RPA energies, oscillator strengths and ||Y||^2 (pyscf's (X, Y) are normalised 2 (X.X - Y.Y) = 1,
   !!     so ||Y||^2(terco) = 2 ||Y||^2(pyscf));
   !!  3. the fitted excitation energies against a FOUR-CENTRE SCF and TDDFT of the same molecule: they differ by
   !!     the fitting error and must stay within FIT_BOUND, a sanity check that the operator is the fitted
   !!     one's and not an unrelated one (it is not a substitute for 1 and 2);
   !!  4. the blocked paths against the unblocked one: chunks of eight auxiliary functions (a partial last chunk) with every block
   !!     resident, with none resident (every chunk rebuilt from the packed B inside each application) and with
   !!     some resident, to round-off (`tddft_df_set_blocks`);
   !!  5. the refusals: the excited-state gradient over a fit (a TDA solve, whose message must survive the approximation
   !!     check -- it was overwritten by an operator-precedence slip), a fit-less solve, an unset fit, a fit of another
   !!     basis, and the refusal plumbing for a device cap below the reserve (`device_memory_set_cap`), after which
   !!     the fit is intact and solves again. A second ctest entry of the water Cartesian HF case runs the anchors with
   !!     TERCO_TDDFT_DF_BLOCK=8 and TERCO_TDDFT_DF_RESIDENT=0 in the environment (the blocked, rebuilt path against
   !!     the pyscf anchors).
   !! The spaces are as check_df_grad_vs_pyscf.F90's: "cart" (the orbital and auxiliary bases Cartesian) and
   !! "sph" (both as their JSON tags them, spherical for cc-pVDZ and def2-universal-jkfit), the latter running
   !! over the Cartesian orbitals `trc_tddft_run` folds from the spherical ones. One (molecule, space,
   !! functional) per process; the SCF is run to pyscf's tight convergence and terco's numbers and its grid go
   !! to the probe file for generate.py.
   !!
   !! Arguments: xyz basis.json aux.json molname space(cart|sph) funcname probe.bin [mode]. Run from test/.
   !! Mode "anchors" (the default) runs 1-3 and the release; mode "extras" runs 4, 5 and the release. Items 4
   !! and 5 do not depend on the molecule, the functional or the space, so ctest runs them once for HF (exchange on)
   !! and once for PBE (exchange off) on water Cartesian instead of in every program (each costs two SCFs and a dozen solves).
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_df, only: df_jk_t
   use trc_tddft_df, only: tddft_df_set_blocks
   use trc_device_memory, only: device_memory_set_cap
   use trc_df_kernels, only: DF_LMAX_AUX
   use trc_tddft_run, only: trc_tddft_run_options_t, trc_tddft_run_result_t, trc_tddft_run_solve, &
                            trc_tddft_run_solve_fitted, trc_tddft_run_gradient
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_vector
   implicit none
   private
   public :: collect_check_tddft_df, tddft_df_needs_more_aux_lmax
   public :: xyz_arg, basis_arg, aux_arg, molname_arg, space_arg, funcname_arg, probe_arg, mode_arg

   real(dp), parameter :: TOL_OMEGA = 1.0e-6_dp
   real(dp), parameter :: TOL_OSC = 1.0e-5_dp
   real(dp), parameter :: TOL_YY = 1.0e-6_dp
   ! Above TOL_OMEGA, TOL_OSC, TOL_YY: the Cartesian TDDFT checks' gates against pyscf (check_tddft_basis,
   ! check_tddft_sph).
   ! Below, TOL_BLOCK_OMEGA and TOL_BLOCK_OSC: blocked against unblocked, the same operator summed in another
   ! order and solved to conv_tol 1e-10.
   real(dp), parameter :: TOL_BLOCK_OMEGA = 1.0e-11_dp
   real(dp), parameter :: TOL_BLOCK_OSC = 1.0e-7_dp
   real(dp), parameter :: FIT_BOUND = 1.0e-4_dp
   !! Fitted against four-centre excitation energies (hartree): the fitting error of def2-universal-jkfit in the
   !! lowest roots is 1.0e-5 to 1.6e-5 (measured, water and formamide, HF/PBE/B3LYP); a wrong operator (a missing
   !! exchange term, a wrong weight) moves them by 1e-2 or more.
   integer, parameter :: NROOTS = 3
   real(dp), parameter :: THRESH = 1.0e-12_dp

   character(len=256) :: xyz_arg = "water.xyz"
   character(len=256) :: basis_arg = "basis_sets/6-31gs.json"
   character(len=256) :: aux_arg = "basis_sets/def2-universal-jkfit.json"
   character(len=16) :: molname_arg = "water"
   character(len=16) :: space_arg = "cart"
   character(len=16) :: funcname_arg = "hf"
   character(len=256) :: probe_arg = "tddft_df_probe.bin"
   character(len=16) :: mode_arg = "anchors"

   type :: fixture_t
      logical :: ready = .false.   ! the bases and the fit exist: test_release frees them
      logical :: solved = .false.  ! the SCF and both solves succeeded
      integer :: nocc = 0, natm = 0
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: la
      type(df_jk_t) :: dfjk
      type(trc_scf_options_t) :: sopts
      type(trc_scf_result_t) :: scf
      type(trc_tddft_run_options_t) :: opts
      type(trc_tddft_run_result_t) :: tda, rpa
   end type fixture_t

   type(fixture_t), save :: fx

contains

   function tddft_df_needs_more_aux_lmax() result(skip)
      !! The auxiliary basis needs more angular momentum than this build's DF_LMAX_AUX: the program skips cleanly,
      !! as check_df_grad_vs_pyscf.F90 does.
      logical :: skip
      type(trc_basis_t) :: aux
      type(error_t) :: err
      integer :: nat, dev
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      dev = trc_bind_device(0)
      call read_xyz(trim(xyz_arg), nat, at_z, r)
      call trc_basis_from_json(trim(aux_arg), nat, at_z, r, aux, err, force_cart=trim(space_arg) == "cart")
      if (err%has_error()) then
         call global_logger%error("check_tddft_df: "//err%get_message())
         error stop 1
      end if
      skip = maxval(aux%sh_l) > DF_LMAX_AUX
      call aux%release()
   end function tddft_df_needs_more_aux_lmax

   subroutine collect_check_tddft_df(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (trim(mode_arg) /= "extras" .and. trim(mode_arg) /= "anchors") then
         call global_logger%error("check_tddft_df: unknown mode '"//trim(mode_arg)//"' (anchors or extras)")
         error stop 1
      end if
      if (trim(mode_arg) == "extras") then
         testsuite = [ &
                     new_unittest("blocked and rebuilt chunks equal the unblocked solve", test_blocks), &
                     new_unittest("refusals", test_refusals), &
                     new_unittest("release", test_release) &
                     ]
      else
         testsuite = [ &
                     new_unittest("fitted TDA energies vs pyscf", test_tda), &
                     new_unittest("fitted RPA energies and ||Y||^2 vs pyscf", test_rpa), &
                     new_unittest("fitted roots within the fitting error of the four-centre ones", test_vs_exact), &
                     new_unittest("release", test_release) &
                     ]
      end if
   end subroutine collect_check_tddft_df

   function ref_name(prefix) result(name)
      character(len=*), intent(in) :: prefix
      character(len=256) :: name
      name = prefix//trim(molname_arg)//"_"//trim(bas_tag())//"_"//trim(space_arg)//"_"//trim(funcname_arg)//".ref"
   end function ref_name

   function bas_tag() result(tag)
      character(len=16) :: tag
      tag = merge("631gs  ", "ccpvdz ", trim(space_arg) == "cart")
   end function bas_tag

   subroutine scf_settings(sopts, df)
      !! pyscf's tight convergence; Kohn-Sham on terco's level-3 grid.
      type(trc_scf_options_t), intent(inout) :: sopts
      logical, intent(in) :: df
      sopts%conv_energy = 1.0e-12_dp
      sopts%conv_density = 1.0e-10_dp
      sopts%conv_diis = 1.0e-11_dp
      sopts%max_iteration = 200
      sopts%grid_level = 3
      sopts%df = df
      if (trim(funcname_arg) /= "hf") sopts%functional = trim(funcname_arg)
   end subroutine scf_settings

   subroutine fixture_build(error)
      !! The bases, the fit, one tightly converged fitted SCF and the TDA and RPA solves over it.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp), allocatable :: at_r(:, :)
      integer, allocatable :: zint(:)
      logical :: cart
      integer :: dev

      if (fx%solved) return
      call check(error, .not. fx%ready, "fixture setup failed in an earlier test; nothing to compare")
      if (allocated(error)) return
      cart = trim(space_arg) == "cart"
      dev = trc_bind_device(0)
      call read_xyz(trim(xyz_arg), fx%natm, zint, at_r)
      call trc_basis_from_json(trim(basis_arg), fx%natm, zint, at_r, fx%bas, err, force_cart=cart)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call trc_basis_from_json(trim(aux_arg), fx%natm, zint, at_r, fx%aux, err, force_cart=cart)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, fx%bas%spherical .eqv. (.not. cart), "the orbital basis is not in the requested space")
      if (allocated(error)) return
      call fx%bas%to_device()
      call fx%aux%to_device()
      call fx%pl%build(fx%bas, THRESH)
      call fx%pl%to_device()
      call fx%la%init(fx%bas%nao)
      fx%nocc = nint(sum(real(zint, dp)))/2
      call fx%dfjk%setup(fx%bas, fx%aux, fx%pl, fx%la, err, nocc_hint=fx%nocc)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      fx%ready = .true.
      call scf_settings(fx%sopts, .true.)
      call trc_scf_run(fx%bas, fx%nocc, fx%nocc, fx%sopts, fx%scf, la=fx%la, scf_dfjk=fx%dfjk)
      call check(error, fx%scf%converged, "the fitted SCF did not converge")
      if (allocated(error)) return

      fx%opts%nroots = NROOTS
      fx%opts%conv_tol = 1.0e-10_dp
      fx%opts%tda = .true.
      call trc_tddft_run_solve_fitted(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, fx%opts, fx%tda, fx%dfjk)
      call check(error, fx%tda%converged, "TDA solve failed: "//fx%tda%err%get_message())
      if (allocated(error)) return
      fx%opts%tda = .false.
      call trc_tddft_run_solve_fitted(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, fx%opts, fx%rpa, fx%dfjk)
      call check(error, fx%rpa%converged, "RPA solve failed: "//fx%rpa%err%get_message())
      if (allocated(error)) return
      fx%solved = .true.
      call write_probe()
   end subroutine fixture_build

   subroutine write_probe()
      !! Terco's numbers and its grid, for generate.py's tddft_df cases.
      integer :: unit
      open (newunit=unit, file=trim(probe_arg), access="stream", form="unformatted", status="replace")
      write (unit) int(merge(fx%bas%nao_sph, fx%bas%nao, fx%bas%spherical), int64), int(fx%natm, int64), &
         int(fx%nocc, int64), int(NROOTS, int64)
      write (unit) molname_arg, bas_tag(), funcname_arg, space_arg
      write (unit) fx%scf%energy
      write (unit) fx%tda%omega
      write (unit) fx%tda%osc
      write (unit) fx%rpa%omega
      write (unit) fx%rpa%osc
      write (unit) fx%rpa%ynorm
      ! Hartree-Fock has no grid; the generator reads zero points.
      if (allocated(fx%scf%grid_w)) then
         write (unit) int(size(fx%scf%grid_w), int64)
         write (unit) fx%scf%grid_r, fx%scf%grid_w
      else
         write (unit) 0_int64
      end if
      close (unit)
   end subroutine write_probe

   subroutine compare_roots(error, res, prefix, with_y, tag)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_run_result_t), intent(in) :: res
      character(len=*), intent(in) :: prefix, tag
      logical, intent(in) :: with_y
      real(dp), allocatable :: om_ref(:), osc_ref(:), yy_ref(:)
      real(dp) :: dom, dosc, dyy, e_ref
      integer :: u

      call fixture_build(error)
      if (allocated(error)) return
      call ref_require(ref_name(prefix), u)
      allocate (om_ref(NROOTS), osc_ref(NROOTS), yy_ref(NROOTS))
      call ref_read_real(u, "E_SCF", e_ref)
      call ref_read_vector(u, "OMEGA", NROOTS, om_ref)
      call ref_read_vector(u, "OSC", NROOTS, osc_ref)
      if (with_y) call ref_read_vector(u, "YY", NROOTS, yy_ref)
      call ref_close(u)
      dom = maxval(abs(res%omega - om_ref))
      dosc = maxval(abs(res%osc - osc_ref))
      call global_logger%info("tddft df "//tag//" "//trim(molname_arg)//"/"//trim(space_arg)//"/"//trim(funcname_arg)// &
                              ": max|d omega| "//to_string(dom, "es10.2")//", max|d f| "//to_string(dosc, "es10.2"))
      call check(error, abs(fx%scf%energy - e_ref) < 1.0e-8_dp, "E_SCF vs pyscf, "// &
                 to_string(abs(fx%scf%energy - e_ref), "es10.2"))
      if (allocated(error)) return
      call check(error, dom < TOL_OMEGA, "omega vs pyscf, "//to_string(dom, "es10.2"))
      if (allocated(error)) return
      call check(error, dosc < TOL_OSC, "oscillator strength vs pyscf, "//to_string(dosc, "es10.2"))
      if (allocated(error)) return
      if (with_y) then
         dyy = maxval(abs(res%ynorm - yy_ref))
         call global_logger%info("tddft df "//tag//": max|d ||Y||^2| "//to_string(dyy, "es10.2"))
         call check(error, dyy < TOL_YY, "||Y||^2 vs pyscf, "//to_string(dyy, "es10.2"))
      else
         call check(error, maxval(res%ynorm) == 0.0_dp, "the TDA has Y = 0")
      end if
   end subroutine compare_roots

   subroutine test_tda(error)
      type(error_type), allocatable, intent(out) :: error
      call fixture_build(error)
      if (allocated(error)) return
      call compare_roots(error, fx%tda, "tddft_df_", .false., "TDA")
   end subroutine test_tda

   subroutine test_rpa(error)
      type(error_type), allocatable, intent(out) :: error
      call fixture_build(error)
      if (allocated(error)) return
      call compare_roots(error, fx%rpa, "tddft_df_rpa_", .true., "RPA")
   end subroutine test_rpa

   subroutine test_vs_exact(error)
      !! The same molecule through a four-centre SCF and TDDFT: the fitted roots must be close, not equal.
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_options_t) :: sopts4
      type(trc_scf_result_t) :: scf4
      type(trc_tddft_run_result_t) :: ex
      type(trc_tddft_run_options_t) :: o
      real(dp) :: d_tda, d_rpa, d_min

      call fixture_build(error)
      if (allocated(error)) return
      call scf_settings(sopts4, .false.)
      call trc_scf_run(fx%bas, fx%nocc, fx%nocc, sopts4, scf4, la=fx%la)
      call check(error, scf4%converged, "the four-centre SCF did not converge")
      if (allocated(error)) return
      o = fx%opts
      o%tda = .true.
      call trc_tddft_run_solve(fx%bas, fx%la, scf4, fx%nocc, sopts4, o, ex)
      call check(error, ex%converged, "four-centre TDA failed: "//ex%err%get_message())
      if (allocated(error)) return
      d_tda = maxval(abs(fx%tda%omega - ex%omega))
      o%tda = .false.
      call trc_tddft_run_solve(fx%bas, fx%la, scf4, fx%nocc, sopts4, o, ex)
      call check(error, ex%converged, "four-centre RPA failed: "//ex%err%get_message())
      if (allocated(error)) return
      d_rpa = maxval(abs(fx%rpa%omega - ex%omega))
      call global_logger%info("tddft df vs four-centre "//trim(molname_arg)//"/"//trim(space_arg)//"/"// &
                              trim(funcname_arg)//": max|d omega| TDA "//to_string(d_tda, "es10.2")//", RPA "// &
                              to_string(d_rpa, "es10.2")//" (bound "//to_string(FIT_BOUND, "es8.1")//")")
      call check(error, d_tda < FIT_BOUND, "TDA fitted vs four-centre, "//to_string(d_tda, "es10.2"))
      if (allocated(error)) return
      call check(error, d_rpa < FIT_BOUND, "RPA fitted vs four-centre, "//to_string(d_rpa, "es10.2"))
      if (allocated(error)) return
      ! The fit is not exact: identical roots would mean the four-centre operator ran.
      d_min = min(d_tda, d_rpa)
      call check(error, d_min > 1.0e-9_dp, "fitted and four-centre roots are identical to 1e-9: no fit in the operator")
   end subroutine test_vs_exact

   subroutine test_blocks(error)
      !! Chunks of eight auxiliary functions: every block resident, none resident (each chunk rebuilt in every
      !! application), and the first three chunks resident.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: QBLK = 8
      integer, parameter :: RESIDENT(3) = [-1, 0, 3*QBLK]
      type(trc_tddft_run_result_t) :: tda, rpa
      type(trc_tddft_run_options_t) :: o
      real(dp) :: dom, dosc, dyy
      integer :: k, naux, nchunk

      call fixture_build(error)
      if (allocated(error)) return
      ! The last chunk must be partial (water Cartesian def2-universal-jkfit has 133 = 16*8 + 5 functions), and
      ! the third plan keeps three chunks resident and rebuilds the rest: a mix of both.
      naux = merge(fx%aux%nao_sph, fx%aux%nao, fx%aux%spherical)
      nchunk = (naux + QBLK - 1)/QBLK
      call global_logger%info("tddft df blocks: "//to_string(naux)//" auxiliary functions, "//to_string(nchunk)// &
                              " chunks of "//to_string(QBLK)//", last chunk "//to_string(naux - (nchunk - 1)*QBLK)// &
                              "; resident plan "//to_string(RESIDENT(3)/QBLK)//" resident, "// &
                              to_string(nchunk - RESIDENT(3)/QBLK)//" rebuilt")
      call check(error, mod(naux, QBLK) /= 0, "the chunk size leaves no partial last chunk, so the tail is untested")
      if (allocated(error)) return
      call check(error, RESIDENT(3) > 0 .and. RESIDENT(3) < naux, "the resident count does not mix resident and rebuilt")
      if (allocated(error)) return
      o = fx%opts
      do k = 1, size(RESIDENT)
         call tddft_df_set_blocks(QBLK, RESIDENT(k))
         o%tda = .true.
         call trc_tddft_run_solve_fitted(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, o, tda, fx%dfjk)
         o%tda = .false.
         call trc_tddft_run_solve_fitted(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, o, rpa, fx%dfjk)
         call tddft_df_set_blocks(0, -1)
         call check(error, tda%converged .and. rpa%converged, "a blocked solve failed: "//tda%err%get_message()// &
                    rpa%err%get_message())
         if (allocated(error)) return
         dom = max(maxval(abs(tda%omega - fx%tda%omega)), maxval(abs(rpa%omega - fx%rpa%omega)))
         dosc = max(maxval(abs(tda%osc - fx%tda%osc)), maxval(abs(rpa%osc - fx%rpa%osc)))
         dyy = maxval(abs(rpa%ynorm - fx%rpa%ynorm))
         call global_logger%info("tddft df blocks (7, resident "//to_string(RESIDENT(k))//"): max|d omega| "// &
                                 to_string(dom, "es10.2")//", max|d f| "//to_string(dosc, "es10.2")// &
                                 ", max|d ||Y||^2| "//to_string(dyy, "es10.2"))
         call check(error, dom < TOL_BLOCK_OMEGA, "blocked omega vs unblocked, "//to_string(dom, "es10.2"))
         if (allocated(error)) return
         call check(error, dosc < TOL_BLOCK_OSC, "blocked f vs unblocked, "//to_string(dosc, "es10.2"))
         if (allocated(error)) return
         call check(error, dyy < TOL_BLOCK_OSC, "blocked ||Y||^2 vs unblocked, "//to_string(dyy, "es10.2"))
         if (allocated(error)) return
      end do
   end subroutine test_blocks

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_run_result_t) :: r
      type(trc_tddft_run_options_t) :: o
      type(df_jk_t) :: unready

      call fixture_build(error)
      if (allocated(error)) return
      ! The gradient over a fit is not implemented: refused, before any work, with its reason.
      r = fx%tda
      o = fx%opts
      o%tda = .true.
      o%target_root = 1
      call trc_tddft_run_gradient(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, o, r)
      call check(error, r%err%has_error(), "the gradient over a --df reference was not refused")
      if (allocated(error)) return
      call check(error, index(r%err%get_message(), "not implemented yet") > 0 .and. &
                 index(r%err%get_message(), "--df") > 0, "the gradient refusal: "//r%err%get_message())
      if (allocated(error)) return
      call check(error, .not. allocated(r%grad), "a refused gradient left a gradient")
      if (allocated(error)) return
      call check(error, index(r%err%get_message(), "approximation") == 0, &
                 "the refusal was replaced by the approximation message: "//r%err%get_message())
      if (allocated(error)) return
      ! A fitted reference without its fit, and with a fit that is not set up.
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, o, r)
      call check(error, r%err%has_error() .and. index(r%err%get_message(), "trc_tddft_run_solve_fitted") > 0, &
                 "a --df solve without its fit was not refused: "//r%err%get_message())
      if (allocated(error)) return
      call trc_tddft_run_solve_fitted(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, o, r, unready)
      call check(error, r%err%has_error() .and. index(r%err%get_message(), "gone") > 0, &
                 "a --df solve with an unset fit was not refused: "//r%err%get_message())
      if (allocated(error)) return
      call refuse_other_basis(error)
      if (allocated(error)) return
      ! The refusal plumbing for a device without room: the cap is below the fixed reserve (RESERVE, 512 MiB, plus
      ! the 256 MiB cuBLAS workspace, both in tddft_df_reserve), so any molecule is refused, with the sizes in the message and the fit left intact (it solves again after). That the
      ! refusal threshold is sized right is not checked here.
      call device_memory_set_cap(600000000_int64)
      call trc_tddft_run_solve_fitted(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, o, r, fx%dfjk)
      call device_memory_set_cap(0_int64)
      call check(error, r%err%has_error() .and. index(r%err%get_message(), "device memory") > 0, &
                 "a device cap below the reserve was not refused with the sizes: "//r%err%get_message())
      if (allocated(error)) return
      call check(error, fx%dfjk%is_ready(), "the fit was left lent after a refused setup")
      if (allocated(error)) return
      call trc_tddft_run_solve_fitted(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, o, r, fx%dfjk)
      call check(error, r%converged, "the fit does not solve again after the refusal: "//r%err%get_message())
   end subroutine test_refusals

   subroutine refuse_other_basis(error)
      !! A fit built for another (smaller) basis is refused, not read past its end.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas2, aux2
      type(trc_pairlist_t) :: pl2
      type(df_jk_t) :: fit2
      type(trc_tddft_run_result_t) :: r
      type(error_t) :: err
      real(dp), allocatable :: at_r(:, :)
      integer, allocatable :: zint(:)
      integer :: natm

      call read_xyz(trim(xyz_arg), natm, zint, at_r)
      call trc_basis_from_json("basis_sets/sto-3g.json", natm, zint, at_r, bas2, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call trc_basis_from_json(trim(aux_arg), natm, zint, at_r, aux2, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call bas2%to_device()
      call aux2%to_device()
      call pl2%build(bas2, THRESH)
      call pl2%to_device()
      call fit2%setup(bas2, aux2, pl2, fx%la, err, nocc_hint=fx%nocc)
      call check(error, .not. err%has_error(), err%get_message())
      if (.not. allocated(error)) then
         call trc_tddft_run_solve_fitted(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, fx%opts, r, fit2)
         call check(error, r%err%has_error() .and. index(r%err%get_message(), "another basis") > 0, &
                    "a fit of another basis was not refused: "//r%err%get_message())
      end if
      call fit2%release()
      call pl2%release()
      call aux2%release()
      call bas2%release()
   end subroutine refuse_other_basis

   subroutine test_release(error)
      type(error_type), allocatable, intent(out) :: error
      if (fx%ready) then
         call fx%dfjk%release()
         call fx%pl%release()
         call fx%la%release()
         call fx%aux%release()
         call fx%bas%release()
         fx%ready = .false.
         fx%solved = .false.
      end if
      call check(error, .not. fx%ready, "released")
   end subroutine test_release

end module check_tddft_df_suite

program check_tddft_df
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use pic_logger, only: global_logger
   use check_tddft_df_suite, only: collect_check_tddft_df, tddft_df_needs_more_aux_lmax, xyz_arg, basis_arg, aux_arg, &
                                   molname_arg, space_arg, funcname_arg, probe_arg, mode_arg
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 2) call get_command_argument(2, basis_arg)
   if (command_argument_count() >= 3) call get_command_argument(3, aux_arg)
   if (command_argument_count() >= 4) call get_command_argument(4, molname_arg)
   if (command_argument_count() >= 5) call get_command_argument(5, space_arg)
   if (command_argument_count() >= 6) call get_command_argument(6, funcname_arg)
   if (command_argument_count() >= 7) call get_command_argument(7, probe_arg)
   if (command_argument_count() >= 8) call get_command_argument(8, mode_arg)
   if (trim(mode_arg) /= "anchors" .and. trim(mode_arg) /= "extras") then
      call global_logger%error("check_tddft_df: unknown mode '"//trim(mode_arg)//"' (anchors or extras)")
      error stop 1
   end if
   stat = 0
   if (tddft_df_needs_more_aux_lmax()) then
      call global_logger%info("SKIP check_tddft_df: the auxiliary basis needs more angular momentum than DF_LMAX_AUX")
      stop 0
   end if
   call run_testsuite(collect_check_tddft_df, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_tddft_df
