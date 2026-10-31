submodule (trc_run) trc_run_json
   !! The `--json` summary writer of the `terco` command (bodies of `trc_run`).
   use json_module, only: json_core, json_value
   use trc_app, only: BOHR_PER_ANGSTROM
   implicit none

contains

   module subroutine run_write_json(r, path)
      !! The run's JSON summary, built with json-fortran, which escapes the
      !! strings and writes the reals at full precision. The key names and
      !! their order are an interface: scripts and the check_cli_grad tests
      !! read them.
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: path
      type(json_core) :: json
      type(json_value), pointer :: root
      call json%initialize()
      call json%create_object(root, "")
      call json%add(root, "program", "terco")
      call json%add(root, "version", VERSION)
      call json%add(root, "geometry", r%geom)
      call json%add(root, "method", method(r))
      call json%add(root, "basis", r%basis_path)
      call json%add(root, "guess", r%guess)
      call json%add(root, "charge", r%chrg)
      call json%add(root, "unpaired electrons", r%uhf)
      call json%add(root, "number of basis functions", int(r%nao))
      call json%add(root, "spherical", r%sph /= 0)
      call json%add(root, "iterations", int(r%niter))
      call json%add(root, "converged", r%conv /= 0)
      call json%add(root, "nuclear repulsion energy", real(r%enuc, dp))
      call json%add(root, "one-electron energy", real(r%eone, dp))
      call json%add(root, "two-electron energy", real(r%etwo, dp))
      call json%add(root, "exchange-correlation energy", real(r%exc, dp))
      if (len(r%pcm) > 0) then
         call json%add(root, "solvation model", r%pcm)
         call json%add(root, "dielectric constant", r%eps)
         call json%add(root, "solvation energy", real(r%epcm, dp))
         call json%add(root, "solvation surface points", int(r%npcm))
      end if
      if (r%rimp2 .and. r%conv /= 0) then
         call json%add(root, "mp2 opposite-spin energy", real(r%eos, dp))
         call json%add(root, "mp2 same-spin energy", real(r%ess, dp))
      end if
      if (r%rimp2 .and. r%conv /= 0 .and. r%dh) then
         call json%add(root, "pt2 opposite-spin fraction", real(r%pt2_cos, dp))
         call json%add(root, "pt2 same-spin fraction", real(r%pt2_css, dp))
         call json%add(root, "scf energy", real(r%etot, dp))
         call json%add(root, "total energy", real(r%etot_dh, dp))
      else
         call json%add(root, "scf energy", real(r%etot, dp))
         call json%add(root, "total energy", real(r%etot, dp))
         if (r%grad .and. allocated(r%gradient)) call add_gradient(json, root, r)
      end if
      if (r%hessian .and. r%hess_ok) call add_hessian(json, root, r)
      if (r%opt) call add_optimization(json, root, r)
      if (r%charges == "hirshfeld" .and. allocated(r%hirshfeld_charges)) call add_charges(json, root, r)
      if (r%fukui == "hirshfeld" .and. allocated(r%fukui_f_plus)) call add_fukui(json, root, r)
      call json%print(root, path)
      if (json%failed()) then
         call json%destroy(root)
         call run_die(r, "cannot write "//path)
      end if
      call json%destroy(root)
   end subroutine run_write_json

   subroutine add_gradient(json, root, r)
      !! The analytic gradient as an array of [x, y, z] triples, one per atom
      !! in `r%gradient`'s own (and `trc_set_molecule`'s) atom order, plus the
      !! RMS and max per-atom norm `gradient_summary` also prints.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: root
      type(run_t), intent(in) :: r
      type(json_value), pointer :: rows
      integer :: ia
      call json%add(root, "gradient status", "validated: single rank (rank 0 alone on a multi-rank run)")
      call json%create_array(rows, "gradient")
      do ia = 1, r%natm
         call json%add(rows, "", r%gradient(:, ia))
      end do
      call json%add(root, rows)
      call json%add(root, "gradient rms", sqrt(sum(r%gradient**2)/real(3*r%natm, dp)))
      call json%add(root, "gradient max", sqrt(maxval(sum(r%gradient**2, dim=1))))
   end subroutine add_gradient

   subroutine add_charges(json, root, r)
      !! Hirshfeld charges as an array, one per atom in `r%hirshfeld_charges`'s
      !! own (and `trc_set_molecule`'s) atom order, plus their sum.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: root
      type(run_t), intent(in) :: r
      call json%add(root, "charges scheme", r%charges)
      call json%add(root, "hirshfeld charges", r%hirshfeld_charges)
      call json%add(root, "hirshfeld charges sum", sum(r%hirshfeld_charges))
   end subroutine add_charges

   subroutine add_fukui(json, root, r)
      !! Hirshfeld-condensed Fukui indices, one array per index in
      !! `r%fukui_f_plus`'s own (and `trc_set_molecule`'s) atom order, plus
      !! each array's sum and the three states' own Hirshfeld charges
      !! (N, N+1, N-1) they were differenced from.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: root
      type(run_t), intent(in) :: r
      type(json_value), pointer :: section
      call json%create_object(section, "fukui")
      call json%add(root, section)
      call json%add(section, "population scheme", r%fukui)
      call json%add(section, "f+", r%fukui_f_plus)
      call json%add(section, "f-", r%fukui_f_minus)
      call json%add(section, "f0", r%fukui_f_zero)
      call json%add(section, "f+ sum", sum(r%fukui_f_plus))
      call json%add(section, "f- sum", sum(r%fukui_f_minus))
      call json%add(section, "charges neutral (N)", r%fukui_q_n)
      call json%add(section, "charges anion (N+1)", r%fukui_q_np1)
      call json%add(section, "charges cation (N-1)", r%fukui_q_nm1)
   end subroutine add_fukui

   !
   ! --hessian's JSON: three blocks, as asked for -- "hessian" (the
   ! matrix itself and the step it was built at), "frequencies"
   ! (trc_vibrational_analysis's own per-mode output) and
   ! "thermochemistry" (trc_thermo_result_t, trc_thermo.F90's own
   ! fields, nothing recomputed here). Row/column k = 3*(atom-1)+
   ! {x,y,z}, the same order `trc_hessian` and `trc_set_molecule` use.
   !
   subroutine add_hessian(json, root, r)
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: root
      type(run_t), intent(in) :: r
      type(json_value), pointer :: hobj, hrows, fobj, tobj, tsub
      integer :: ia, ndof

      ndof = 3*r%natm

      call json%create_object(hobj, "hessian")
      call json%add(hobj, "step_bohr", r%hess_step)
      call json%create_array(hrows, "matrix")
      do ia = 1, ndof
         call json%add(hrows, "", r%hess(ia, :))
      end do
      call json%add(hobj, hrows)
      call json%add(root, hobj)

      call json%create_object(fobj, "frequencies")
      call json%add(fobj, "n_modes", size(r%freq_cm1))
      call json%add(fobj, "n_imaginary", r%thermo%n_imaginary)
      call json%add(fobj, "is_linear", r%is_linear)
      call json%add(fobj, "frequencies_cm1", r%freq_cm1)
      call json%add(fobj, "reduced_masses_amu", r%redmass_amu)
      call json%add(fobj, "force_constants_hartree_bohr2", r%force_const)
      call json%add(root, fobj)

      call json%create_object(tobj, "thermochemistry")
      call json%add(tobj, "temperature_K", r%thermo%temperature_k)
      call json%add(tobj, "pressure_atm", r%pressure_atm)
      call json%add(tobj, "symmetry_number", r%sym_number)
      call json%add(tobj, "spin_multiplicity", r%mult)
      call json%add(tobj, "masses", r%masses_kind)
      call json%add(tobj, "zero_point_energy_hartree", r%thermo%zpe_hartree)

      call json%create_object(tsub, "thermal_corrections_hartree")
      call json%add(tsub, "to_energy", r%thermo%delta_e_hartree)
      call json%add(tsub, "to_enthalpy", r%thermo%delta_h_hartree)
      call json%add(tsub, "to_gibbs", r%thermo%delta_g_hartree)
      call json%add(tobj, tsub)

      call json%create_object(tsub, "total_energies_hartree")
      call json%add(tsub, "electronic", r%thermo%e_elec_hartree)
      call json%add(tsub, "electronic_plus_zpe", r%thermo%e_zpe_hartree)
      call json%add(tsub, "electronic_plus_thermal_E", r%thermo%e_tot_hartree)
      call json%add(tsub, "electronic_plus_thermal_H", r%thermo%h_tot_hartree)
      call json%add(tsub, "electronic_plus_thermal_G", r%thermo%g_tot_hartree)
      call json%add(tobj, tsub)

      call json%create_object(tsub, "entropy_cal_mol_K")
      call json%add(tsub, "translational", r%thermo%s_trans_calmolk)
      call json%add(tsub, "rotational", r%thermo%s_rot_calmolk)
      call json%add(tsub, "vibrational", r%thermo%s_vib_calmolk)
      call json%add(tsub, "electronic", r%thermo%s_elec_calmolk)
      call json%add(tsub, "total", r%thermo%s_tot_calmolk)
      call json%add(tobj, tsub)

      call json%create_object(tsub, "heat_capacity_cal_mol_K")
      call json%add(tsub, "translational", r%thermo%cv_trans_calmolk)
      call json%add(tsub, "rotational", r%thermo%cv_rot_calmolk)
      call json%add(tsub, "vibrational", r%thermo%cv_vib_calmolk)
      call json%add(tsub, "electronic", r%thermo%cv_elec_calmolk)
      call json%add(tsub, "total", r%thermo%cv_tot_calmolk)
      call json%add(tobj, tsub)

      call json%add(root, tobj)
   end subroutine add_hessian

   !
   ! --opt's JSON block: converged, step/evaluation counts, the final
   ! energy and gradient, the coordinate system, the trajectory's
   ! per-step energies and the final geometry. The final geometry is in
   ! Angstrom, named accordingly so a reader never has to guess the unit
   ! the way it would from a bare "final_geometry".
   !
   subroutine add_optimization(json, root, r)
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: root
      type(run_t), intent(in) :: r
      type(json_value), pointer :: oobj, rows, csub
      integer :: ia

      call json%create_object(oobj, "optimization")
      call json%add(oobj, "converged", r%opt_converged)
      call json%add(oobj, "steps", r%opt_steps)
      call json%add(oobj, "energy_evaluations", r%opt_evals)
      call json%add(oobj, "final_energy", real(r%etot, dp))
      call json%add(oobj, "final_gradient_max", r%opt_gmax)
      call json%add(oobj, "final_gradient_rms", r%opt_grms)
      call json%add(oobj, "gradient_tolerance", r%opt_gtol)
      call json%add(oobj, "coordinates", opt_coords_name(r%opt_coords_kind))
      call json%add(oobj, "guess", trim(r%opt_guess_arg))
      call json%add(oobj, "target", opt_target_name(r%opt_target_kind))
      if (r%opt_target_kind == 1) then
         call json%add(oobj, "saddle_method", opt_saddle_method_name(r%opt_saddle_method_kind))
      else
         call json%add(oobj, "algorithm", opt_algorithm_name(r%opt_algorithm_kind))
      end if
      call json%add(oobj, "connect", r%opt_connect)
      call json%add(oobj, "hessian_evaluations", r%opt_hess_evals)
      if (allocated(r%opt_traj_energy)) then
         call json%add(oobj, "trajectory_energies", r%opt_traj_energy)
      end if
      if (allocated(r%opt_final_xyz)) then
         call json%create_array(rows, "final_geometry_angstrom")
         do ia = 1, r%natm
            call json%add(rows, "", r%opt_final_xyz(:, ia)/BOHR_PER_ANGSTROM)
         end do
         call json%add(oobj, rows)
      end if
      if (r%opt_connect) then
         ! Each sub-search's OWN convergence verdict and gradient, never
         ! another sub-search's -- the top-level converged/final_gradient_max
         ! above are the TS search's own (trc_dlfind_bridge's own override),
         ! but a reader who wants the two downhill minima's separately
         ! finds them here instead of having to infer them.
         call json%add(oobj, "ts_found", r%opt_have_ts)
         if (r%opt_have_ts) then
            call json%create_object(csub, "ts")
            call json%add(csub, "energy", r%opt_ts_energy)
            call json%add(csub, "converged", r%opt_ts_converged)
            call json%add(csub, "final_gradient_max", r%opt_ts_gmax)
            call json%add(oobj, csub)
         end if
         call json%add(oobj, "minimum_minus_found", r%opt_have_min_minus)
         if (r%opt_have_min_minus) then
            call json%create_object(csub, "minimum_minus")
            call json%add(csub, "energy", r%opt_min_minus_energy)
            call json%add(csub, "converged", r%opt_min_minus_converged)
            call json%add(csub, "final_gradient_max", r%opt_min_minus_gmax)
            call json%add(oobj, csub)
         end if
         call json%add(oobj, "minimum_plus_found", r%opt_have_min_plus)
         if (r%opt_have_min_plus) then
            call json%create_object(csub, "minimum_plus")
            call json%add(csub, "energy", r%opt_min_plus_energy)
            call json%add(csub, "converged", r%opt_min_plus_converged)
            call json%add(csub, "final_gradient_max", r%opt_min_plus_gmax)
            call json%add(oobj, csub)
         end if
      end if
      ! --opt-neb's own block: the image energies, the max-energy
      ! (climbing) image and whether the run converged, plus the
      ! per-image geometry so a JSON reader never has to go back to the
      ! multi-frame xyz for it.
      if (r%opt_neb_have) call add_neb(json, oobj, r)
      call json%add(root, oobj)
   end subroutine add_optimization

   !
   ! --opt-neb's own JSON block, nested under "optimization": image
   ! energies, the max-energy ("climbing") image's index and geometry, and
   ! whether this run converged by the SAME criterion the rest of --opt's
   ! JSON uses (the max Cartesian gradient component of the optimised
   ! images against --opt-gtol -- trc_dlfind_bridge's own result assembly
   ! excludes frozen endpoints from that test, see its own comment).
   !
   subroutine add_neb(json, oobj, r)
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: oobj
      type(run_t), intent(in) :: r
      type(json_value), pointer :: nobj, rows, geomrows
      integer :: k, ia

      call json%create_object(nobj, "neb")
      call json%add(nobj, "n_images", r%opt_neb_n_images)
      call json%add(nobj, "endpoints", trim(r%opt_neb_endpoints_arg))
      call json%add(nobj, "climbing_image_requested", r%opt_neb_climb)
      call json%add(nobj, "climbing_image_index", r%opt_neb_climb_image)
      call json%add(nobj, "climbing_image_energy", r%opt_neb_image_energy(r%opt_neb_climb_image))
      call json%add(nobj, "image_energies", r%opt_neb_image_energy)
      if (allocated(r%opt_neb_scf_iters)) call json%add(nobj, "scf_iterations_per_image", r%opt_neb_scf_iters)
      call json%create_array(rows, "climbing_image_geometry_angstrom")
      do ia = 1, r%natm
         call json%add(rows, "", r%opt_neb_image_xyz(:, ia, r%opt_neb_climb_image)/BOHR_PER_ANGSTROM)
      end do
      call json%add(nobj, rows)
      call json%create_array(geomrows, "path_geometries_angstrom")
      do k = 1, r%opt_neb_n_images
         block
            type(json_value), pointer :: framerows
            call json%create_array(framerows, "")
            do ia = 1, r%natm
               call json%add(framerows, "", r%opt_neb_image_xyz(:, ia, k)/BOHR_PER_ANGSTROM)
            end do
            call json%add(geomrows, framerows)
         end block
      end do
      call json%add(nobj, geomrows)
      call json%add(oobj, nobj)
   end subroutine add_neb

end submodule trc_run_json
