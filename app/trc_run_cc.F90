submodule (trc_run) trc_run_cc
   !! The CLI side of closed-shell RI-CCSD and RI-CCSD(T) (`--riccsd`,
   !! `--riccsd-t`, `--cc-*`): the options, what the method refuses before any
   !! work, and the call into the library (bodies of `trc_run`).
   use pic_strings, only: to_string
   use trc_capi, only: trc_set_rccsd, trc_set_rccsd_convergence, trc_set_rccsd_algorithm, trc_run_rccsd, &
                       trc_rccsd_energy, trc_rccsd_dimensions, trc_energy, TRC_ERR_NOCONV
   implicit none

   ! An option as the user met it: the flag, and the TOML key that gives the same setting. Literals, so that no
   ! character function result is involved (see ladder_code).
   character(len=*), parameter :: N_ECONV = "--cc-econv ([riccsd] energy_tol)"
   character(len=*), parameter :: N_TCONV = "--cc-tconv ([riccsd] amplitude_tol)"
   character(len=*), parameter :: N_ITER = "--cc-iterations ([riccsd] max_iterations)"
   character(len=*), parameter :: N_DIIS = "--cc-diis ([riccsd] diis_size)"
   character(len=*), parameter :: N_LADDER = "--cc-ladder ([riccsd] ladder)"
   character(len=*), parameter :: N_TRIPLES = "--cc-triples ([riccsd] triples)"
   character(len=*), parameter :: N_BUDGET = "--cc-triples-budget ([riccsd] triples_budget)"
   integer, parameter :: MAX_DIIS = 50
   !! DIIS keeps two copies of the amplitudes per vector, o v + o^2 v^2 doubles each: 50 is far past what helps

contains

   pure integer function ladder_code(word) result(code)
      !! The `--cc-ladder` word as trc_set_rccsd_algorithm's code (0 auto, 1 assembled, 2 direct), -1 for anything
      !! else, in any case. It lower-cases into a fixed-size local and returns an integer. The same job done with
      !! the allocatable-character `run_to_lower_alloc` gave a result 100 characters long with a garbage tail in
      !! this unit on the GPU build (nvfortran 26.5), at a call that is correct and that works in `trc_main`; the
      !! failure came and went with unrelated statements around the call, so no form of it is safe to rely on.
      character(len=*), intent(in) :: word
      character(len=16) :: low
      integer :: i, c

      code = -1
      if (len_trim(word) > len(low)) return
      low = word
      do i = 1, len(low)
         c = iachar(low(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) low(i:i) = achar(c + 32)
      end do
      if (low == "auto") then
         code = 0
      else if (low == "assembled") then
         code = 1
      else if (low == "direct") then
         code = 2
      end if
   end function ladder_code

   pure integer function triples_code(word) result(code)
      !! The `--cc-triples` word as trc_set_rccsd_algorithm's `triples_routine` (0 auto, 1 reference, 2 tiled, 3 ijk), -1 for
      !! anything else, in any case. Fixed-size local and an integer result, for the reason given at `ladder_code`.
      character(len=*), intent(in) :: word
      character(len=16) :: low
      integer :: i, c

      code = -1
      if (len_trim(word) > len(low)) return
      low = word
      do i = 1, len(low)
         c = iachar(low(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) low(i:i) = achar(c + 32)
      end do
      if (low == "auto") then
         code = 0
      else if (low == "reference") then
         code = 1
      else if (low == "tiled") then
         code = 2
      else if (low == "ijk") then
         code = 3
      end if
   end function triples_code

   module subroutine run_cc_parse(r, nfrozen)
      !! Reads `--riccsd`, `--riccsd-t` (which implies `--riccsd`) and the `--cc-*` options. Every combination the
      !! method does not cover is refused here, before the SCF: a correlated step that fails after an hour of SCF is
      !! worse than one that fails on the command line. The library refuses the same cases again for a caller that
      !! does not come through here.
      type(run_t), intent(inout) :: r
      integer, intent(in) :: nfrozen   !! --frozen
      character(len=:), allocatable :: sval
      real(dp) :: gib
      logical :: any_option

      r%cc_triples = r%cli%is_set("riccsd-t")
      r%cc = r%cli%is_set("riccsd") .or. r%cc_triples
      any_option = r%cli%is_set("cc-econv") .or. r%cli%is_set("cc-tconv") .or. r%cli%is_set("cc-iterations") .or. &
                   r%cli%is_set("cc-diis") .or. r%cli%is_set("cc-ladder") .or. r%cli%is_set("cc-triples") .or. &
                   r%cli%is_set("cc-triples-budget")
      if (.not. r%cc) then
         if (any_option) then
            call die_options_alone(r)
         end if
         return
      end if
      if (r%cli%is_set("rimp2")) call run_die(r, "--riccsd and --rimp2 cannot be combined: the RI-MP2 correlation "// &
                                          "energy of the same reference is part of the --riccsd report")
      if (len(r%func) > 0) call run_die(r, "--riccsd needs a Hartree-Fock reference; --func "//r%func// &
                                    " is a Kohn-Sham one (a double hybrid has its own MP2 part)")
      if (r%uhf /= 0) then
         if (r%cli%is_set("uhf")) then
            call run_die(r, "--riccsd is closed-shell only; --uhf "//to_string(r%uhf)//" asks for unpaired electrons")
         else
            call run_die(r, "--riccsd is closed-shell only; the "//to_string(r%nelec)//" electrons of this molecule "// &
                         "(charge "//to_string(r%chrg)//") cannot all be paired: check --charge")
         end if
      end if
      if (r%grad .or. r%hessian .or. r%opt) call run_die(r, "--riccsd has no gradient: it cannot be combined with "// &
                                                      "--grad, --hessian or --opt")
      if (r%ecp_core > 0) call run_die(r, "--riccsd does not support an effective core potential")
      if (r%nranks > 1) call run_die(r, "--riccsd runs on one rank only ("//to_string(r%nranks)//" ranks given)")
      if (nfrozen >= r%nelec/2) call run_die(r, "--frozen "//to_string(nfrozen)//" leaves no correlated orbital: this "// &
                                         "molecule has "//to_string(r%nelec/2)//" occupied orbitals")
      if (len(r%pcm) > 0 .and. r%rank0) then
         call global_logger%warning("--riccsd over a --pcm reference: the correlation is computed on the orbitals of "// &
                                    "the solvated SCF and has no solvent response of its own (as --rimp2)")
      end if
      if (r%cli%is_set("cc-econv")) then
         call r%cli%get("cc-econv", r%cc_econv, r%err)
         call run_check(r, N_ECONV)
         if (.not. (r%cc_econv > 0.0_dp .and. r%cc_econv <= huge(1.0_dp))) call run_die(r, N_ECONV//" must be positive")
      end if
      if (r%cli%is_set("cc-tconv")) then
         call r%cli%get("cc-tconv", r%cc_tconv, r%err)
         call run_check(r, N_TCONV)
         if (.not. (r%cc_tconv > 0.0_dp .and. r%cc_tconv <= huge(1.0_dp))) call run_die(r, N_TCONV//" must be positive")
      end if
      if (r%cli%is_set("cc-iterations")) then
         call r%cli%get("cc-iterations", r%cc_iterations, r%err)
         call run_check(r, N_ITER)
         if (r%cc_iterations < 1) call run_die(r, N_ITER//" must be at least 1")
      end if
      if (r%cli%is_set("cc-diis")) then
         call r%cli%get("cc-diis", r%cc_diis, r%err)
         call run_check(r, N_DIIS)
         if (r%cc_diis < 1 .or. r%cc_diis > MAX_DIIS) call run_die(r, N_DIIS//" must be between 1 (DIIS off) and "// &
                                                                to_string(MAX_DIIS)//": each vector costs two copies "// &
                                                                "of the amplitudes")
      end if
      if (r%cli%is_set("cc-ladder")) then
         call r%cli%get("cc-ladder", sval, r%err)
         call run_check(r, N_LADDER)
         r%cc_ladder = ladder_code(sval)
         if (r%cc_ladder < 0) call run_die(r, N_LADDER//" must be auto, assembled or direct, got '"//trim(sval)//"'")
      end if
      if (r%cli%is_set("cc-triples")) then
         if (.not. r%cc_triples) call run_die(r, N_TRIPLES//" needs --riccsd-t ([correlation] method = 'ri-ccsd(t)')")
         call r%cli%get("cc-triples", sval, r%err)
         call run_check(r, N_TRIPLES)
         r%cc_triples_routine = triples_code(sval)
         if (r%cc_triples_routine < 0) call run_die(r, N_TRIPLES//" must be auto, reference, tiled or ijk, got '"// &
                                                    trim(sval)//"'")
      end if
      if (r%cli%is_set("cc-triples-budget")) then
         if (.not. r%cc_triples) call run_die(r, N_BUDGET//" needs --riccsd-t ([correlation] method = 'ri-ccsd(t)')")
         call r%cli%get("cc-triples-budget", gib, r%err)
         call run_check(r, N_BUDGET)
         ! Any positive size is a budget: one below a single triple (48 o^3 bytes) is how the (T) refusal is reached.
         if (.not. (gib > 0.0_dp .and. gib <= 1.0e6_dp)) call run_die(r, N_BUDGET//" must be positive and at "// &
                                                                      "most 1e6 GiB")
         r%cc_budget_mib = gib*1024.0_dp
      end if
   end subroutine run_cc_parse

   module subroutine run_cc_configure(r, nfrozen, aux_block)
      !! The library's RI-CCSD settings from the parsed options. The solver settings are passed only when one was
      !! given, so the library's defaults (1e-8, 1e-6, 100 iterations, DIIS 8) stand otherwise.
      type(run_t), intent(in) :: r
      integer, intent(in) :: nfrozen, aux_block

      call run_api(r, trc_set_rccsd(r%h, merge(1_c_int, 0_c_int, r%cc_triples), int(nfrozen, c_int), &
                                    int(aux_block, c_int)))
      if (r%cc_econv > 0.0_dp .or. r%cc_tconv > 0.0_dp .or. r%cc_iterations > 0 .or. r%cc_diis > 0) then
         call run_api(r, trc_set_rccsd_convergence(r%h, real(r%cc_econv, c_double), real(r%cc_tconv, c_double), &
                                                   int(r%cc_iterations, c_int), int(r%cc_diis, c_int)))
      end if
      if (r%cc_ladder /= 0 .or. r%cc_triples_routine /= 0 .or. r%cc_budget_mib > 0.0_dp) then
         call run_api(r, trc_set_rccsd_algorithm(r%h, int(r%cc_ladder, c_int), int(r%cc_triples_routine, c_int), &
                                                 real(r%cc_budget_mib, c_double)))
      end if
   end subroutine run_cc_configure

   subroutine die_options_alone(r)
      !! The refusal of `--cc-*` options given without `--riccsd`, in a routine of its own.
      !!
      !! Written in `run_cc_parse`, where it is the first `run_die` call, nvfortran 26.9 passed `run_die` a wrong
      !! length for the message whatever the message was (a concatenation, one named constant, a local variable of
      !! the exact length, a thirteen-character literal): about twenty characters too long on x86 (the text, a null,
      !! neighbouring constants), an uninitialised value on aarch64 (757101167, then FTN-F-STR_MALLOC or a crash in
      !! the string copy). The later `run_die` calls of that routine are not affected, and 26.5 did neither. Called
      !! from here the length is right on both. The CLI check looks for the text and for a null in every refusal.
      type(run_t), intent(in) :: r

      call run_die(r, "the --cc-* options (the [riccsd] table) need --riccsd or --riccsd-t ([correlation] method = "// &
                   "'ri-ccsd' or 'ri-ccsd(t)')")
   end subroutine die_options_alone

   subroutine keep_message(r)
      !! The library's last message into `r%cc_message`, copied by hand (no character function result)
      type(run_t), intent(inout) :: r
      integer, parameter :: MSG_LEN = 512
      character(kind=c_char) :: buffer(MSG_LEN)
      integer(c_int) :: rc
      integer :: i

      r%cc_message = ""
      rc = trc_message(r%h, buffer, int(size(buffer), c_int))
      if (rc /= TRC_OK) return
      do i = 1, size(buffer)
         if (buffer(i) == achar(0)) exit
         r%cc_message(i:i) = buffer(i)
      end do
   end subroutine keep_message

   module subroutine run_cc_execute(r)
      !! Runs RI-CCSD(T) on the SCF the handle holds. A refusal stops the run with the library's message. A CCSD
      !! that did not converge is not a refusal: its last iterate is read and kept, `r%cc_converged` is false, and
      !! the summary and the JSON say so before the run ends with a failure status. A (T) that failed after a
      !! converged CCSD is kept the same way (`r%cc_t_failed`): the CCSD energies are still a result.
      type(run_t), intent(inout) :: r
      integer(c_int) :: status, rc, flag_converged, flag_triples

      status = trc_run_rccsd(r%h)
      call keep_message(r)
      flag_converged = 0
      flag_triples = 0
      rc = trc_rccsd_energy(r%h, r%cc_e_mp2, r%cc_e_ccsd, r%cc_e_t, r%cc_t1_diagnostic, r%cc_niter, &
                            flag_converged, flag_triples)
      if (status == TRC_OK .or. status == TRC_ERR_NOCONV) then
         call run_api(r, rc)
      else if (rc /= TRC_OK .or. flag_converged == 0) then
         call run_api(r, status)
      else
         r%cc_t_failed = .true.
      end if
      rc = trc_rccsd_dimensions(r%h, r%cc_nfrozen, r%cc_nactive, r%cc_nvirtual, r%cc_naux)
      call run_api(r, rc)
      r%cc_ran = .true.
      r%cc_converged = flag_converged /= 0
      r%cc_has_triples = flag_triples /= 0
      if (r%cc_converged) call run_api(r, trc_energy(r%h, r%etot_total))
   end subroutine run_cc_execute

end submodule trc_run_cc
