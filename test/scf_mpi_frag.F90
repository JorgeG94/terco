!
! GAMESS/DDI and fragment methods, imitated, on 4 ranks.
!
! trc_capi's `fcomm` handles used to accept only MPI_COMM_WORLD (see
! trc_build_comm in trc_capi.F90). That is not what either caller this test
! stands in for ever hands terco:
!
!   * GAMESS's DDI puts DATA SERVER ranks in MPI_COMM_WORLD that never call
!     into any library at all -- they sit in DDI's own server loop. GAMESS
!     hands libraries DDI's *compute* communicator, a strict sub-communicator
!     of MPI_COMM_WORLD. The old code either rejected that handle outright,
!     or -- had it been passed MPI_COMM_WORLD instead, to get past the
!     check -- would have called MPI_Comm_dup(MPI_COMM_WORLD) from the
!     compute ranks only, which never returns: a dup is collective over
!     every rank of the communicator it names, and the data servers never
!     call it.
!   * A fragment method hands terco a different sub-communicator per
!     fragment, several of which are alive on the same MPI_COMM_WORLD at
!     once, none of them the world communicator either.
!
! Four scenarios, all against a build with TERCO_ENABLE_MPI=ON and run under
! `mpiexec -np 4`:
!
!   a. DATA SERVERS. World splits into "compute" (ranks 0-1) and "servers"
!      (ranks 2-3). The servers never call trc_scf_mpi; they only reach a
!      barrier on MPI_COMM_WORLD, which the compute ranks reach in turn only
!      after their own run (through trc_scf_mpi, on the COMPUTE
!      communicator's own handle) is done. That barrier is what a hang would
!      show up as -- the ctest TIMEOUT below is the backstop if it does not.
!   b. FRAGMENTS. World splits into two 2-rank groups running DIFFERENT
!      molecules AT THE SAME TIME, each on its own sub-communicator's
!      handle: water on ranks 0-1, the OH radical (open-shell, so this also
!      covers UHF over a communicator) on ranks 2-3.
!   c. OWNERSHIP. Every sub-communicator this test hands to terco in (a) and
!      (b) is still usable for an ordinary MPI collective afterwards, and
!      THIS PROGRAM frees it, not terco -- proof that trc_build_comm
!      duplicated it rather than adopting it.
!   d. GPU BINDING. Each rank prints the device trc_scf_mpi bound it to;
!      rank 0 checks that the two ranks of a communicator that shares a node
!      (every rank here does -- one node, four V100s) land on distinct
!      devices, now that binding goes by NODE-LOCAL rank rather than by
!      comm%rank() in whatever communicator the caller happened to hand in.
!
! Every collective energy is checked against this same rank's OWN
! single-rank run (fcomm = -1) on the same molecule, to ~1e-10 Eh: that is
! the property trc_build_comm exists to preserve -- an arbitrary
! communicator must reproduce trc_scf on one rank -- not a frozen literal.
!
program scf_mpi_frag
   use trc_boys, only: dp
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize, allreduce, &
                          allgather, MPI_MAX, MPI_MIN
   use trc_test_basis, only: read_xyz, unpack_basis
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char
   use trc_c_interfaces, only: trc_basis_create, trc_basis_destroy, trc_basis_nao, trc_scf_mpi, TRC_OK
   ! This test plays the part of the CALLER (GAMESS's DDI, or a fragment
   ! driver), and a caller reaches terco only through the plain Fortran
   ! handle of a communicator it built itself. It builds those communicators
   ! with pic-mpi rather than `use mpi_f08`, and hands terco
   ! comm%fortran_handle(): the same integer MPI_Comm_c2f or the legacy
   ! `mpi` module would give, and the only way to build on every backend --
   ! some MPIs (Cray MPICH under nvfortran) offer no mpi_f08 module at all.
#ifdef _OPENACC
   use openacc, only: acc_get_device_num, acc_device_nvidia
#endif
   implicit none

   integer :: natm_w, natm_o, nsh_w, nsh_o, maxnp_w, maxnp_o
   integer, allocatable :: z_w(:), z_o(:), l_w(:), np_w(:), l_o(:), np_o(:)
   real(dp), allocatable :: r_w(:, :), r_o(:, :), e_w(:, :), c_w(:, :), sr_w(:, :)
   real(dp), allocatable :: e_o(:, :), c_o(:, :), sr_o(:, :)
   ! Pre-converted to the C API's own kinds and held in named variables,
   ! not passed as expressions: `l`/`nprim` are assumed-size C dummies
   ! (trc_c_interfaces.f90), and an assumed-size actual argument must be an
   ! array variable, not the result of an inline int()/real() conversion.
   integer(c_int), allocatable :: lc_w(:), npc_w(:), lc_o(:), npc_o(:)
   real(c_double), allocatable :: zat_w(:), zat_o(:)
   type(c_ptr) :: hbas_w, hbas_o
   integer(c_int) :: nao_w, nao_o, rc
   real(c_double), allocatable :: dmat_w(:), eps_w(:), dmat_o(:), eps_o(:)
   type(trc_basis_t) :: jbas
   type(error_t) :: err
   logical :: ok

   ok = .true.
   call pic_mpi_init()

   ! The two molecules, 6-31G: water (closed shell, the reference RHF every
   ! other MPI test in this suite already uses) and the OH radical
   ! (open-shell, nalpha /= nbeta -- UHF), so scenario (b) is not running
   ! the same physics on both fragments.
   call read_xyz('water.xyz', natm_w, z_w, r_w)
   call trc_basis_from_json('basis_sets/6-31g.json', natm_w, z_w, r_w, jbas, err)
   if (err%has_error()) then
      print '(a)', 'scf_mpi_frag: '//err%get_message(); call pic_mpi_finalize(); stop 1
   end if
   call unpack_basis(jbas, nsh_w, l_w, np_w, e_w, c_w, sr_w, maxnp_w)
   allocate (lc_w(nsh_w), npc_w(nsh_w), zat_w(natm_w))
   lc_w = int(l_w, c_int); npc_w = int(np_w, c_int); zat_w = real(z_w, c_double)
   rc = trc_basis_create(int(nsh_w, c_int), int(maxnp_w, c_int), lc_w, npc_w, &
                         e_w, c_w, sr_w, int(natm_w, c_int), zat_w, r_w, hbas_w)
   if (rc /= TRC_OK) ok = .false.
   rc = trc_basis_nao(hbas_w, nao_w)
   if (rc /= TRC_OK) ok = .false.
   allocate (dmat_w(nao_w*nao_w), eps_w(nao_w))

   call read_xyz('oh.xyz', natm_o, z_o, r_o)
   call trc_basis_from_json('basis_sets/6-31g.json', natm_o, z_o, r_o, jbas, err)
   if (err%has_error()) then
      print '(a)', 'scf_mpi_frag: '//err%get_message(); call pic_mpi_finalize(); stop 1
   end if
   call unpack_basis(jbas, nsh_o, l_o, np_o, e_o, c_o, sr_o, maxnp_o)
   allocate (lc_o(nsh_o), npc_o(nsh_o), zat_o(natm_o))
   lc_o = int(l_o, c_int); npc_o = int(np_o, c_int); zat_o = real(z_o, c_double)
   rc = trc_basis_create(int(nsh_o, c_int), int(maxnp_o, c_int), lc_o, npc_o, &
                         e_o, c_o, sr_o, int(natm_o, c_int), zat_o, r_o, hbas_o)
   if (rc /= TRC_OK) ok = .false.
   rc = trc_basis_nao(hbas_o, nao_o)
   if (rc /= TRC_OK) ok = .false.
   ! The OH radical is run unrestricted (nalpha=5, nbeta=4, nspin=2 -- see
   ! trc_scf_mpi's own comment: "the caller sizes dmat and eps as (nao, nao,
   ! nspin) and (nao, nspin)"), so both hold TWO spin blocks, not one. Sizing
   ! these as nao_o*nao_o/nao_o alone -- easy to do by analogy with water's
   ! closed-shell nspin=1 arrays above -- overruns the buffer by exactly one
   ! spin block and corrupts the heap; that surfaced far away and confusingly,
   ! as an MPI_Comm_dup/MPI_Comm_free crash deep in OpenMPI on a later call.
   allocate (dmat_o(nao_o*nao_o*2), eps_o(nao_o*2))

   ! With no MPI (pic-mpi's serial backend) there is one rank, so the size
   ! check below skips: the DDI and fragment topologies need real ranks to
   ! mean anything. Not registered as a ctest in that configuration either
   ! (see test/CMakeLists.txt) -- it only keeps the source building
   ! everywhere, in the house style.
   block
      type(comm_t) :: world, sub_comm
      integer :: wrank, wsize, color, sum_ranks
      integer(c_int) :: fcomm
      real(c_double) :: e_ref, e_xc, e_a, e_b, spread_max, spread_min
      integer :: niter
      integer :: dev_a, dev_b
      integer :: all_dev_a(4), all_dev_b(4)

      world = comm_world()
      wrank = world%rank()
      wsize = world%size()
      if (wsize /= 4) then
         if (wrank == 0) print '(a,i0,a)', "scf_mpi_frag: needs 4 ranks, got ", wsize, " -- SKIP"
      else
         dev_a = -1; dev_b = -1

         !
         ! (a) DATA SERVERS. Ranks 0-1 are "compute", 2-3 are "servers" that
         ! never touch terco.
         !
         color = merge(0, 1, wrank < 2)
         sub_comm = world%split_by(color)
         if (sub_comm%is_null()) ok = .false.

         if (wrank < 2) then
            fcomm = int(sub_comm%fortran_handle(), c_int)
            ! This rank's own single-rank reference, for the same molecule.
            rc = trc_scf_mpi(-1_c_int, hbas_w, 5_c_int, 5_c_int, c_null_char, 3_c_int, &
                             1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, 0_c_int, &
                             e_ref, e_xc, dmat_w, eps_w, niter)
            if (rc /= TRC_OK) ok = .false.
            ! The collective run, on the COMPUTE communicator's own handle --
            ! not MPI_COMM_WORLD. If trc_build_comm still special-cased
            ! MPI_COMM_WORLD, or still duplicated it regardless of what
            ! `fcomm` names, this would either be refused or hang here,
            ! because rank 2 and 3 are on their way to a barrier that has
            ! nothing to do with this call and never make an MPI_Comm_dup
            ! of anything.
            rc = trc_scf_mpi(fcomm, hbas_w, 5_c_int, 5_c_int, c_null_char, 3_c_int, &
                             1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, 0_c_int, &
                             e_a, e_xc, dmat_w, eps_w, niter)
            if (rc /= TRC_OK) ok = .false.
            if (abs(e_a - e_ref) > 1.0e-10_dp) ok = .false.
#ifdef _OPENACC
            dev_a = acc_get_device_num(acc_device_nvidia)
#endif
            ! (c) OWNERSHIP: the compute communicator terco just duplicated
            ! and finalized its own copy of is still this program's, and an
            ! ordinary collective on it still works.
            sum_ranks = wrank
            call allreduce(sub_comm, sum_ranks)
            if (sum_ranks /= 1) ok = .false.  ! ranks 0 + 1
            spread_max = e_a
            call allreduce(sub_comm, spread_max, MPI_MAX)
            spread_min = e_a
            call allreduce(sub_comm, spread_min, MPI_MIN)
            if (spread_max - spread_min > 1.0e-12_dp) ok = .false.
            ! This program frees the communicator it made, itself; terco
            ! freed only the duplicate it made from the handle.
            call sub_comm%finalize()
         else
            ! The "servers": a real DDI server loop would be off doing its
            ! own thing right now. This one has nothing to do, so it frees
            ! its half of the split immediately and heads straight for the
            ! barrier below, where it waits for compute to catch up.
            call sub_comm%finalize()
         end if
         ! Every rank reaches this together only once compute's run above
         ! has returned -- the servers were already here. A hang anywhere
         ! in (a) shows up as this barrier never completing; the ctest
         ! TIMEOUT is the backstop.
         call world%barrier()

         !
         ! (b) FRAGMENTS. Two 2-rank groups, two different molecules, at the
         ! same time: water on 0-1 again, the OH radical (UHF) on 2-3.
         !
         color = merge(0, 1, wrank < 2)
         sub_comm = world%split_by(color)
         if (sub_comm%is_null()) ok = .false.
         fcomm = int(sub_comm%fortran_handle(), c_int)

         if (wrank < 2) then
            rc = trc_scf_mpi(-1_c_int, hbas_w, 5_c_int, 5_c_int, c_null_char, 3_c_int, &
                             1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, 0_c_int, &
                             e_ref, e_xc, dmat_w, eps_w, niter)
            if (rc /= TRC_OK) ok = .false.
            rc = trc_scf_mpi(fcomm, hbas_w, 5_c_int, 5_c_int, c_null_char, 3_c_int, &
                             1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, 0_c_int, &
                             e_b, e_xc, dmat_w, eps_w, niter)
         else
            rc = trc_scf_mpi(-1_c_int, hbas_o, 5_c_int, 4_c_int, c_null_char, 3_c_int, &
                             1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, 0_c_int, &
                             e_ref, e_xc, dmat_o, eps_o, niter)
            if (rc /= TRC_OK) ok = .false.
            rc = trc_scf_mpi(fcomm, hbas_o, 5_c_int, 4_c_int, c_null_char, 3_c_int, &
                             1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, 0_c_int, &
                             e_b, e_xc, dmat_o, eps_o, niter)
         end if
         if (rc /= TRC_OK) ok = .false.
         if (abs(e_b - e_ref) > 1.0e-10_dp) ok = .false.
#ifdef _OPENACC
         dev_b = acc_get_device_num(acc_device_nvidia)
#endif

         ! (c) OWNERSHIP again, for the fragment communicator.
         spread_max = e_b
         call allreduce(sub_comm, spread_max, MPI_MAX)
         spread_min = e_b
         call allreduce(sub_comm, spread_min, MPI_MIN)
         if (spread_max - spread_min > 1.0e-12_dp) ok = .false.
         call sub_comm%finalize()

         !
         ! (d) GPU BINDING. Gather every rank's device for both scenarios
         ! and let rank 0 print and check them: dev_a is only meaningful for
         ! ranks 0-1 (the compute pair of (a)), dev_b for all four (every
         ! rank belongs to a fragment in (b)). All four ranks are on the one
         ! node this machine has, so a communicator that is node-local
         ! (both sub_comm's here are) must bind its members to DISTINCT
         ! devices -- the whole point of binding by node-local rank instead
         ! of by comm%rank() in whatever communicator the caller handed in.
         !
         call allgather(world, dev_a, all_dev_a)
         call allgather(world, dev_b, all_dev_b)
         if (wrank == 0) then
            block
               integer :: i
               do i = 1, 4
                  print '(a,i0,a,i0,a,i0)', "scf_mpi_frag: world rank ", i - 1, "  device(a) ", &
                     all_dev_a(i), "  device(b) ", all_dev_b(i)
               end do
            end block
#ifdef _OPENACC
            if (all_dev_a(1) == all_dev_a(2)) ok = .false.  ! (a)'s compute pair, one node
            if (all_dev_b(1) == all_dev_b(2)) ok = .false.  ! (b)'s water fragment, one node
            if (all_dev_b(3) == all_dev_b(4)) ok = .false.  ! (b)'s OH fragment, one node
#endif
         end if

         if (wrank == 0) then
            print '(a,es9.2,a,es9.2)', "scf_mpi_frag: (a) data servers  |e-e_ref| ", &
               abs(e_a - e_ref), "  (b) fragments  spread ", spread_max - spread_min
         end if
      end if
      call world%finalize()
   end block

   rc = trc_basis_destroy(hbas_w)
   if (rc /= TRC_OK) ok = .false.
   rc = trc_basis_destroy(hbas_o)
   if (rc /= TRC_OK) ok = .false.

   print '(a)', merge("scf_mpi_frag: PASS", "scf_mpi_frag: FAIL", ok)
   call pic_mpi_finalize()
   if (.not. ok) stop 1
end program scf_mpi_frag
