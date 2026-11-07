/*
 * terco driven from C, through terco.h only.
 *
 * WHY THIS TEST EXISTS
 * ---------------------
 * trc_c_interfaces.f90 proves the ABI to a Fortran caller that cannot read
 * terco's own `.mod` files (test/link_gfortran.F90); this is the same
 * proof for a C or C++ caller, against the OTHER file
 * tools/gen_capi_headers.py derives from trc_capi.F90: terco.h. It compiles
 * with the project's C compiler (not nvfortran, not gfortran) and links
 * libterco, so a mismatch between the generated header and the library's
 * real ABI -- an argument dropped, a type that does not match, a constant
 * that drifted -- is a build or run failure here, in terco's own test
 * suite, before any host code ever sees it.
 *
 * It runs a minimal context: water/6-31G RHF through trc_create ...
 * trc_run_scf ... trc_energy, and checks the energy against the SAME
 * value test/link_gfortran.F90 and test/scf_rhf.F90 get (trc_basis_from_
 * json's own reference, 2026-10-08). It also calls trc_abi_version (and
 * checks it against the header's own TRC_ABI_MAJOR/MINOR) and
 * trc_functional_info("b3lyp"), since those are the two handle-free
 * entries this task added alongside the context ones.
 */
#include "terco.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int must(int rc, const char *what) {
   if (rc != TRC_OK) {
      fprintf(stderr, "link_c: %s failed with status %d\n", what, rc);
      return 1;
   }
   return 0;
}

int main(void) {
   int rc, major, minor, flags, dev;
   double exx, c_os, c_ss;
   trc_handle_t h = NULL;
   double z[3] = {8.0, 1.0, 1.0};
   /* test/water.xyz, Angstrom -> Bohr (CODATA 2018, AGENTS.md). */
   const double ang = 1.8897261254578281;
   double xyz_ang[9] = {
       0.0000000,  0.0000000,  0.1173000,
       0.0000000,  0.7572000, -0.4692000,
       0.0000000, -0.7572000, -0.4692000};
   double xyz[9];
   double energy;
   int i;

   /* --- ABI version: the header's own constants must be what the ----------
    * library just reported, or this header does not describe this .so. */
   rc = trc_abi_version(&major, &minor);
   if (must(rc, "trc_abi_version")) return 1;
   printf("link_c: terco ABI %d.%d\n", major, minor);
   if (major != TRC_ABI_MAJOR || minor != TRC_ABI_MINOR) {
      fprintf(stderr, "link_c: header says ABI %d.%d, library says %d.%d\n",
              TRC_ABI_MAJOR, TRC_ABI_MINOR, major, minor);
      return 1;
   }

   /* Rank 0 must bind device 0 when there is a GPU, and report -1 on a
      host build: mod(0, ndev) either way. */
   rc = trc_select_device(0, &dev);
   if (must(rc, "trc_select_device")) return 1;
   printf("link_c: trc_select_device(0) chose device %d\n", dev);
   if (dev != 0 && dev != -1) {
      fprintf(stderr, "link_c: rank 0 should bind device 0 (or -1 without a GPU), got %d\n", dev);
      return 1;
   }

   /* --- a functional lookup needing no handle at all ----------------------*/
   rc = trc_functional_info("b3lyp", &flags, &exx, &c_os, &c_ss);
   if (must(rc, "trc_functional_info(b3lyp)")) return 1;
   if (!(flags & 1)) {
      fprintf(stderr, "link_c: b3lyp should be known/runnable (flags=%d)\n", flags);
      return 1;
   }
   if (!(flags & 2)) {
      fprintf(stderr, "link_c: b3lyp should be a hybrid (flags=%d)\n", flags);
      return 1;
   }
   if (fabs(exx - 0.2) > 1.0e-12) {
      fprintf(stderr, "link_c: b3lyp exx should be 0.2, got %.12f\n", exx);
      return 1;
   }
   /* an unknown name must come back non-OK with flags = 0 */
   rc = trc_functional_info("not-a-functional", &flags, &exx, &c_os, &c_ss);
   if (rc == TRC_OK || flags != 0) {
      fprintf(stderr, "link_c: an unknown functional name should refuse with flags = 0\n");
      return 1;
   }

   /* --- the minimal context: water/6-31G RHF -------------------------------*/
   for (i = 0; i < 9; i++) xyz[i] = xyz_ang[i] * ang;

   rc = trc_create(&h);
   if (must(rc, "trc_create")) return 1;
   rc = trc_set_molecule(h, 3, z, xyz, 0, 1);
   if (must(rc, "trc_set_molecule")) return 1;
   rc = trc_set_basis_json(h, "basis_sets/6-31g.json");
   if (must(rc, "trc_set_basis_json")) return 1;

   /* --- shell-data accessors, on the basis handle this context holds ----*/
   {
      trc_handle_t basis;
      int nshell, maxnp, spherical;
      /* basis_sets/6-31g.json's own shell list, water atom order O/H/H
       * (test/water.xyz): O is s(6), s(3)/p(3), s(1)/p(1); each H is
       * s(3), s(1). */
      int want_atom[9] = {1, 1, 1, 1, 1, 2, 2, 3, 3};
      int want_l[9] = {0, 0, 1, 0, 1, 0, 0, 0, 0};
      int want_np[9] = {6, 3, 3, 1, 1, 3, 1, 3, 1};
      int atom[9], l[9], nprim[9];
      double exps[54], coefs[54];
      int i;

      rc = trc_context_basis(h, &basis);
      if (must(rc, "trc_context_basis")) return 1;

      rc = trc_basis_shell_count(basis, &nshell);
      if (must(rc, "trc_basis_shell_count")) return 1;
      if (nshell != 9) {
         fprintf(stderr, "link_c: water/6-31g should have 9 shells, got %d\n", nshell);
         return 1;
      }
      rc = trc_basis_shell_max_primitives(basis, &maxnp);
      if (must(rc, "trc_basis_shell_max_primitives")) return 1;
      if (maxnp != 6) {
         fprintf(stderr, "link_c: water/6-31g's widest shell has 6 primitives, got %d\n", maxnp);
         return 1;
      }
      rc = trc_basis_shell_info(basis, atom, l, nprim);
      if (must(rc, "trc_basis_shell_info")) return 1;
      for (i = 0; i < 9; i++) {
         if (atom[i] != want_atom[i] || l[i] != want_l[i] || nprim[i] != want_np[i]) {
            fprintf(stderr, "link_c: shell %d: atom/l/nprim = %d/%d/%d, want %d/%d/%d\n",
                    i, atom[i], l[i], nprim[i], want_atom[i], want_l[i], want_np[i]);
            return 1;
         }
      }
      rc = trc_basis_shell_primitives(basis, exps, coefs);
      if (must(rc, "trc_basis_shell_primitives")) return 1;
      /* exponents are untouched by trc_basis_json's own renormalisation
       * (trc_basis_shell_primitives' own doc), so these ARE the JSON
       * file's literal numbers for oxygen's tightest and loosest 1s
       * primitives. */
      if (fabs(exps[0] - 5484.671660) > 1.0e-4) {
         fprintf(stderr, "link_c: shell 0 primitive 0 exponent %.6f, want 5484.671660\n", exps[0]);
         return 1;
      }
      if (fabs(exps[5] - 5.799635340) > 1.0e-6) {
         fprintf(stderr, "link_c: shell 0 primitive 5 exponent %.6f, want 5.799635340\n", exps[5]);
         return 1;
      }
      rc = trc_basis_is_spherical(basis, &spherical);
      if (must(rc, "trc_basis_is_spherical")) return 1;
      if (spherical != 0) {
         fprintf(stderr, "link_c: 6-31g is Cartesian, trc_basis_is_spherical should be 0, got %d\n", spherical);
         return 1;
      }
   }

   /* a basis built DIRECTLY (trc_basis_create), not through JSON: here
    * trc_basis_shell_primitives must hand back exactly what went in --
    * oxygen's own 1s shell, basis_sets/6-31g.json's literal exponents and
    * coefficients, BSE's published values (trc_basis_shell_primitives'
    * own doc: no renormalisation on this path). */
   {
      trc_handle_t basis2;
      int l2[1] = {0};
      int nprim2[1] = {6};
      double exps2[6] = {5484.671660, 825.2349460, 188.0469580, 52.96450000, 16.89757040, 5.799635340};
      double coefs2[6] = {0.0018310744, 0.0139501722, 0.0684450781, 0.2327143360, 0.4701928980, 0.3585208530};
      double centres2[3] = {0.0, 0.0, 0.0};
      double zatm2[1] = {8.0};
      double ratm2[3] = {0.0, 0.0, 0.0};
      double exps_out[6], coefs_out[6];
      int i;

      rc = trc_basis_create(1, 6, l2, nprim2, exps2, coefs2, centres2, 1, zatm2, ratm2, &basis2);
      if (must(rc, "trc_basis_create")) return 1;
      rc = trc_basis_shell_primitives(basis2, exps_out, coefs_out);
      if (must(rc, "trc_basis_shell_primitives (direct)")) return 1;
      for (i = 0; i < 6; i++) {
         if (fabs(exps_out[i] - exps2[i]) > 1.0e-10) {
            fprintf(stderr, "link_c: direct basis primitive %d exponent %.10f, want %.10f\n",
                    i, exps_out[i], exps2[i]);
            return 1;
         }
         if (fabs(coefs_out[i] - coefs2[i]) > 1.0e-8) {
            fprintf(stderr, "link_c: direct basis primitive %d coefficient %.10f, want %.10f\n",
                    i, coefs_out[i], coefs2[i]);
            return 1;
         }
      }
      rc = trc_basis_destroy(basis2);
      if (must(rc, "trc_basis_destroy")) return 1;
   }

   /* --- physical constants, against trc_constants' CODATA 2018 own ------*/
   {
      double value;
      rc = trc_physical_constant(TRC_CONST_BOHR_TO_ANGSTROM, &value);
      if (must(rc, "trc_physical_constant(BOHR_TO_ANGSTROM)")) return 1;
      if (fabs(value - 0.529177210903) > 1.0e-12) {
         fprintf(stderr, "link_c: BOHR_TO_ANGSTROM %.12f, want 0.529177210903\n", value);
         return 1;
      }
      rc = trc_physical_constant(TRC_CONST_HARTREE_TO_EV, &value);
      if (must(rc, "trc_physical_constant(HARTREE_TO_EV)")) return 1;
      if (fabs(value - 27.211386245988) > 1.0e-9) {
         fprintf(stderr, "link_c: HARTREE_TO_EV %.12f, want 27.211386245988\n", value);
         return 1;
      }
      rc = trc_physical_constant(TRC_CONST_HARTREE_TO_KCALMOL, &value);
      if (must(rc, "trc_physical_constant(HARTREE_TO_KCALMOL)")) return 1;
      if (fabs(value - 627.5094740631) > 1.0e-6) {
         fprintf(stderr, "link_c: HARTREE_TO_KCALMOL %.10f, want 627.5094740631\n", value);
         return 1;
      }
      rc = trc_physical_constant(TRC_CONST_AMU_TO_ME, &value);
      if (must(rc, "trc_physical_constant(AMU_TO_ME)")) return 1;
      if (fabs(value - 1822.888486209) > 1.0e-6) {
         fprintf(stderr, "link_c: AMU_TO_ME %.10f, want 1822.888486209\n", value);
         return 1;
      }
      rc = trc_physical_constant(TRC_CONST_KB_SI, &value);
      if (must(rc, "trc_physical_constant(KB_SI)")) return 1;
      if (fabs(value - 1.380649e-23) > 1.0e-30) {
         fprintf(stderr, "link_c: KB_SI %.6e, want 1.380649e-23\n", value);
         return 1;
      }
      rc = trc_physical_constant(TRC_CONST_NA_SI, &value);
      if (must(rc, "trc_physical_constant(NA_SI)")) return 1;
      if (fabs(value - 6.02214076e23) > 1.0e10) {
         fprintf(stderr, "link_c: NA_SI %.6e, want 6.02214076e23\n", value);
         return 1;
      }
      /* an id this header does not know must refuse with value left at 0 */
      value = -1.0;
      rc = trc_physical_constant(999, &value);
      if (rc == TRC_OK || value != 0.0) {
         fprintf(stderr, "link_c: an unknown constant id should refuse with value = 0\n");
         return 1;
      }
   }

   /* ABI 1.24: a fresh context has symmetry on. This file's first run is the
    * plain C1 one its symmetry block below compares against, so it asks for
    * "off"; the default has its own block at the end. */
   rc = trc_set_symmetry(h, "off", 0.0);
   if (must(rc, "trc_set_symmetry(off) before the C1 run")) return 1;
   rc = trc_set_method(h, "", 3); /* empty functional name: Hartree-Fock */
   if (must(rc, "trc_set_method")) return 1;
   rc = trc_set_convergence(h, 1.0e-10, 1.0e-7, 100, 10);
   if (must(rc, "trc_set_convergence")) return 1;
   rc = trc_run_scf(h);
   if (must(rc, "trc_run_scf")) return 1;
   rc = trc_energy(h, &energy);
   if (must(rc, "trc_energy")) return 1;

   printf("link_c: RHF/6-31G water energy %.10f\n", energy);
   /* trc_basis_from_json's own 6-31G water RHF value -- the same one
    * test/link_gfortran.F90 and test/scf_rhf.F90 assert. */
   if (fabs(energy + 75.983974465717) > 1.0e-8) {
      fprintf(stderr, "link_c: energy %.12f does not match -75.983974465717\n", energy);
      return 1;
   }

   /* --- symmetry (ABI 1.22): water is C2v (A1 A2 B1 B2). The C1 gradient first, then the same context with
    * symmetry = auto: the energy and the norm of the gradient (which the move
    * to the standard frame rotates but cannot change) must agree. ----------*/
   {
      double g1[9], g2[9], n1 = 0.0, n2 = 0.0, e_sym, resid, rot[9], shift[3], geom[9];
      char name[16], det[16];
      int nirrep, dims[4], irreps[13], nao, k, nocc[4] = {0, 0, 0, 0};

      rc = trc_run_gradient(h);
      if (must(rc, "trc_run_gradient (C1)")) return 1;
      rc = trc_gradient(h, g1);
      if (must(rc, "trc_gradient (C1)")) return 1;
      /* before any symmetric run the getters have nothing to say */
      if (trc_symmetry_group(h, name, 16) != TRC_ERR_STATE) {
         fprintf(stderr, "link_c: trc_symmetry_group should be TRC_ERR_STATE before a symmetric run\n");
         return 1;
      }
      if (trc_set_symmetry(h, "bogus", 0.0) != TRC_ERR_BADARG) {
         fprintf(stderr, "link_c: trc_set_symmetry(bogus) should be refused\n");
         return 1;
      }
      rc = trc_set_symmetry(h, "auto", 0.0);
      if (must(rc, "trc_set_symmetry(auto)")) return 1;
      rc = trc_run_scf(h);
      if (must(rc, "trc_run_scf (symmetric)")) return 1;
      rc = trc_energy(h, &e_sym);
      if (must(rc, "trc_energy (symmetric)")) return 1;
      if (fabs(e_sym - energy) > 1.0e-8) {
         fprintf(stderr, "link_c: symmetric energy %.12f differs from C1 %.12f\n", e_sym, energy);
         return 1;
      }
      rc = trc_symmetry_group(h, name, 16);
      if (must(rc, "trc_symmetry_group")) return 1;
      rc = trc_symmetry_detected(h, det, 16);
      if (must(rc, "trc_symmetry_detected")) return 1;
      printf("link_c: water symmetry: detected %s, used %s\n", det, name);
      if (strcmp(name, "C2v") != 0 || strcmp(det, "C2v") != 0) {
         fprintf(stderr, "link_c: water should be detected C2v and run in C2v, got %s / %s\n", det, name);
         return 1;
      }
      rc = trc_symmetry_residual(h, &resid);
      if (must(rc, "trc_symmetry_residual")) return 1;
      if (resid > 1.0e-6) {
         fprintf(stderr, "link_c: the water geometry is symmetric, residual %.3e\n", resid);
         return 1;
      }
      rc = trc_symmetry_frame(h, rot, shift);
      if (must(rc, "trc_symmetry_frame")) return 1;
      /* the rotation is orthogonal */
      for (k = 0; k < 3; k++) {
         double d = rot[3*k]*rot[3*k] + rot[3*k+1]*rot[3*k+1] + rot[3*k+2]*rot[3*k+2];
         if (fabs(d - 1.0) > 1.0e-12) {
            fprintf(stderr, "link_c: frame row %d has norm^2 %.14f\n", k, d);
            return 1;
         }
      }
      rc = trc_symmetry_nirrep(h, &nirrep);
      if (must(rc, "trc_symmetry_nirrep")) return 1;
      rc = trc_symmetry_irrep_name(h, 1, name, 16);
      if (must(rc, "trc_symmetry_irrep_name(1)")) return 1;
      if (nirrep != 4 || strcmp(name, "A1") != 0) {
         fprintf(stderr, "link_c: C2v has four irreps starting with A1, got %d / %s\n", nirrep, name);
         return 1;
      }
      rc = trc_symmetry_irrep_name(h, 4, name, 16);
      if (must(rc, "trc_symmetry_irrep_name(4)")) return 1;
      if (strcmp(name, "B2") != 0) {
         fprintf(stderr, "link_c: irrep 4 of C2v should be B2, got %s\n", name);
         return 1;
      }
      if (trc_symmetry_irrep_name(h, 5, name, 16) != TRC_ERR_BADARG) {
         fprintf(stderr, "link_c: irrep 5 of C2v should be refused\n");
         return 1;
      }
      rc = trc_symmetry_irrep_dim(h, dims);
      if (must(rc, "trc_symmetry_irrep_dim")) return 1;
      rc = trc_nao(h, &nao);
      if (must(rc, "trc_nao")) return 1;
      if (dims[0] + dims[1] + dims[2] + dims[3] != nao) {
         fprintf(stderr, "link_c: SO dimensions %d + %d + %d + %d do not sum to nao %d\n", dims[0], dims[1], dims[2], dims[3], nao);
         return 1;
      }
      rc = trc_mo_irreps(h, irreps);
      if (must(rc, "trc_mo_irreps")) return 1;
      for (k = 0; k < nao; k++) {
         if (irreps[k] < 1 || irreps[k] > 4) {
            fprintf(stderr, "link_c: orbital %d has irrep %d\n", k + 1, irreps[k]);
            return 1;
         }
         if (k < 5) nocc[irreps[k] - 1]++;   /* water: 5 doubly occupied orbitals */
      }
      printf("link_c: water occupied orbitals per irrep %d %d %d %d (dims %d %d %d %d)\n", nocc[0], nocc[1], nocc[2], nocc[3],
             dims[0], dims[1], dims[2], dims[3]);
      if (nocc[0] != 3 || nocc[1] != 0 || nocc[2] != 1 || nocc[3] != 1) {
         fprintf(stderr, "link_c: water in C2v has 3 a1, 1 b1 (out of plane) and 1 b2 (in plane) occupied orbitals\n");
         return 1;
      }
      rc = trc_geometry(h, geom);
      if (must(rc, "trc_geometry")) return 1;
      rc = trc_run_gradient(h);
      if (must(rc, "trc_run_gradient (symmetric)")) return 1;
      rc = trc_gradient(h, g2);
      if (must(rc, "trc_gradient (symmetric)")) return 1;
      for (k = 0; k < 9; k++) { n1 += g1[k]*g1[k]; n2 += g2[k]*g2[k]; }
      if (fabs(sqrt(n1) - sqrt(n2)) > 1.0e-6) {
         fprintf(stderr, "link_c: gradient norm %.9f (C1) vs %.9f (symmetric)\n", sqrt(n1), sqrt(n2));
         return 1;
      }
      /* what is not supported is refused, not run in C1 */
      if (trc_run_rimp2(h) != TRC_ERR_UNSUPPORTED || trc_run_hessian(h) != TRC_ERR_UNSUPPORTED) {
         fprintf(stderr, "link_c: RI-MP2 and the Hessian should be refused with symmetry on\n");
         return 1;
      }
      /* a subgroup of the detected group: Cs on water runs, in Cs, with the C2v energy */
      rc = trc_set_symmetry(h, "cs", 0.0);
      if (must(rc, "trc_set_symmetry(cs)")) return 1;
      rc = trc_run_scf(h);
      if (must(rc, "trc_run_scf (Cs on water)")) return 1;
      rc = trc_energy(h, &e_sym);
      if (must(rc, "trc_energy (Cs on water)")) return 1;
      if (fabs(e_sym - energy) > 1.0e-8) {
         fprintf(stderr, "link_c: Cs energy %.12f differs from C1 %.12f\n", e_sym, energy);
         return 1;
      }
      rc = trc_symmetry_group(h, name, 16);
      if (must(rc, "trc_symmetry_group (Cs)")) return 1;
      rc = trc_symmetry_detected(h, det, 16);
      if (must(rc, "trc_symmetry_detected (Cs)")) return 1;
      if (strcmp(name, "Cs") != 0 || strcmp(det, "C2v") != 0) {
         fprintf(stderr, "link_c: Cs requested on water should be used as Cs of a detected C2v, got %s / %s\n", name, det);
         return 1;
      }
      /* back to off: a later run is the plain C1 run */
      rc = trc_set_symmetry(h, "off", 0.0);
      if (must(rc, "trc_set_symmetry(off)")) return 1;
      rc = trc_run_scf(h);
      if (must(rc, "trc_run_scf (off again)")) return 1;
      if (trc_symmetry_group(h, name, 16) != TRC_ERR_STATE) {
         fprintf(stderr, "link_c: no symmetry request, no symmetry result\n");
         return 1;
      }

      /* the default request (a fresh context's own, ABI 1.24): symmetric, and
       * reported in the input frame -- the geometry is the input, not the
       * standard frame, and the energy and the gradient norm are the C1 ones */
      rc = trc_set_symmetry(h, "default", 0.0);
      if (must(rc, "trc_set_symmetry(default)")) return 1;
      rc = trc_run_scf(h);
      if (must(rc, "trc_run_scf (default)")) return 1;
      rc = trc_energy(h, &e_sym);
      if (must(rc, "trc_energy (default)")) return 1;
      if (fabs(e_sym - energy) > 1.0e-8) {
         fprintf(stderr, "link_c: default-symmetry energy %.12f differs from C1 %.12f\n", e_sym, energy);
         return 1;
      }
      rc = trc_symmetry_group(h, name, 16);
      if (must(rc, "trc_symmetry_group (default)")) return 1;
      if (strcmp(name, "C2v") != 0) {
         fprintf(stderr, "link_c: the default request should run water in C2v, got %s\n", name);
         return 1;
      }
      rc = trc_symmetry_report_frame(h, &k);
      if (must(rc, "trc_symmetry_report_frame")) return 1;
      if (k != 1) {
         fprintf(stderr, "link_c: the default request should report the input frame\n");
         return 1;
      }
      rc = trc_geometry(h, geom);
      if (must(rc, "trc_geometry (default)")) return 1;
      for (k = 0; k < 9; k++) {
         if (fabs(geom[k] - xyz[k]) > 1.0e-5) {
            fprintf(stderr, "link_c: the default request must report the input geometry, coordinate %d is %.8f not %.8f\n",
                    k, geom[k], xyz[k]);
            return 1;
         }
      }
      rc = trc_set_symmetry_report_frame(h, 0);
      if (must(rc, "trc_set_symmetry_report_frame(0)")) return 1;
      rc = trc_symmetry_report_frame(h, &k);
      if (must(rc, "trc_symmetry_report_frame (standard)")) return 1;
      if (k != 0) {
         fprintf(stderr, "link_c: report frame 0 should read back as the standard frame\n");
         return 1;
      }
   }

   rc = trc_destroy(h);
   if (must(rc, "trc_destroy")) return 1;

   printf("link_c: PASS\n");
   return 0;
}
