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

   rc = trc_destroy(h);
   if (must(rc, "trc_destroy")) return 1;

   printf("link_c: PASS\n");
   return 0;
}
