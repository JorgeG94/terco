!
! The consumer-facing API: a device-resident basis object and batched drivers.
!
! WHAT THIS IS FOR
! ----------------
! Everything below the kernels is per-class `!$acc routine seq` code that
! computes ONE shell block. That is the right shape for a generator and the
! wrong shape for a caller. This module is the boundary: a caller hands over a
! basis once, the data goes to the device once and stays there, and each
! driver is a single `do concurrent` over shell pairs.
!
! SHAPE OF THE INTERFACE
! ----------------------
! The basis is accepted in libcint's `atm`/`bas`/`env` layout, because that is
! what mqc already has -- it carries a libcint backend, so no caller-side
! conversion is needed to try this. It is converted ONCE into a
! structure-of-arrays container that the kernels can read coalesced, rather
! than being walked in that packed integer form on the device.
!
!     type(trc_basis_t) :: b
!     call b%from_libcint(atm, natm, bas, nshell, env)
!     call b%to_device()
!     call trc_1e(b, smat, tmat, vmat)      ! nao x nao each
!     call b%release()
!
! WHY ONE CALL FOR S, T AND V
! ---------------------------
! They share every shell-pair quantity -- zeta, P, the Gaussian product
! prefactor -- and a caller forming a core Hamiltonian wants all three. Asking
! for them separately would read the primitive data three times to save
! nothing.
!
! NO ATOMICS HERE
! ---------------
! Distinct shell pairs write disjoint blocks of the output matrices, so each
! thread owns its output outright. That is the one respect in which the
! one-electron drivers are simpler than the Fock build, which has to scatter
! into six overlapping blocks per quartet.
!
! The API is a facade over three lower-layer modules, re-exported here so a
! caller needs only `use trc_api`: trc_basis_types (core: the basis and the
! pair list), trc_device_binding (core: one GPU per rank) and trc_integrals
! (eri: the 1e, multipole, PCM and DF drivers and their LMAX ceilings).
! Library code below the API layer uses those modules directly.
!
module trc_api
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_device_binding, only: trc_bind_device, trc_bind_device_env, trc_local_rank_env
   use trc_integrals, only: trc_1e, trc_1e_gradient, trc_multipoles, trc_df_2c, trc_df_3c, &
                            trc_pcm_potential, trc_pcm_fock, PCM_LMAX, TRC_NMULT, TRC_MAXL
   implicit none
   private
   public :: trc_basis_t, trc_pairlist_t
   public :: trc_1e, trc_1e_gradient, trc_multipoles, trc_df_2c, trc_df_3c
   public :: trc_pcm_potential, trc_pcm_fock, PCM_LMAX
   public :: trc_bind_device, trc_bind_device_env, trc_local_rank_env
   public :: TRC_NMULT, TRC_MAXL
end module trc_api
