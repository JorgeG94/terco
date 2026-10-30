!
! The batched one-electron, multipole, PCM and density-fitting integral
! drivers over a `trc_basis_t` and its `trc_pairlist_t`: the host-side entry
! points that check the build's angular-momentum ceilings and hand explicit
! arrays to the generated kernel drivers, each a single `do concurrent` over
! shell pairs. Every argument must already be device-resident.
!
! An eri-layer module, next to the generated kernels it launches, so that the
! SCF, PCM, RI-MP2 and the guesses can call it without depending on the API.
! `trc_api` re-exports every public name here.
!
! RYS ROOT TABLES. The density-fitting kernels read the Rys root tables
! (trc_rys_table1..9). Their `declare create` storage does not carry the host
! initialiser to the device, so trc_df_2c and trc_df_3c push them
! (rys_select_tables_init) before every launch: a few hundred kB of `update
! device`, nothing beside a 3c build, and correct whichever device is
! current at the time. A host build compiles the push away. The other
! drivers here read the Boys table only, which trc_basis_t's `build` pushes.
!
module trc_integrals
   use pic_types, only: dp
   use trc_tables, only: LMAX
   use trc_rys_select, only: rys_select_tables_init
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_1e_kernels, only: one_e_driver, ONE_E_LMAX
   use trc_1e_grad_kernels, only: one_e_grad_driver, ONE_E_GRAD_LMAX
   use trc_mult_kernels, only: mult_driver, TRC_NMULT, MULT_LMAX
   use trc_df_kernels, only: df_2c_driver, df_3c_driver, DF_LMAX, DF_LMAX_AUX
   use trc_pcm_kernels, only: pcm_potential_driver, pcm_fock_driver, PCM_LMAX
   implicit none
   private
   public :: trc_1e, trc_1e_gradient, trc_multipoles, trc_df_2c, trc_df_3c
   public :: trc_pcm_potential, trc_pcm_fock, PCM_LMAX
   public :: TRC_NMULT

   !
   ! Angular-momentum limits, taken FROM THE GENERATORS rather than restated.
   !
   ! The four-centre ERIs and everything else have different ceilings on
   ! purpose. A four-centre class carries momentum on all four centres, so its
   ! kernels grow as L^4 and (dd|dd) is already 20k lines; the one-electron and
   ! density-fitting classes carry it on two or three, so f and g are
   ! affordable there while they would not be for the ERIs. Splitting the two
   ! is what lets the fitting basis outrank the orbital basis, which is how
   ! fitting bases are normally built.
   !
   ! TRC_MAXL is the limit for THIS module's paths (1e and DF). The
   ! four-centre limit is TRC_LMAX, a preprocessor symbol consumed by
   ! trc_binkernel and its generated kernels; the two do not have to agree.
   !
   integer, parameter, public :: TRC_MAXL = &
      max(ONE_E_LMAX, max(DF_LMAX, DF_LMAX_AUX))

   !
   ! COMPILE-TIME ASSERTION: the one-electron kernels must reach as high as the
   ! four-centre ones.
   !
   ! "The two do not have to agree" above is true only in the direction where
   ! 1e is HIGHER. The other direction is not a limitation, it is a wrong
   ! answer: a build whose four-centre path outranks ONE_E_LMAX accepts a basis
   ! it cannot do one-electron integrals for, and builds H from whatever
   ! one_e_dispatch returns for a key it has no case for. That is a converged
   ! SCF on a wrong Hamiltonian, with no diagnostic anywhere.
   !
   ! It holds today because both ceilings are 2. It is asserted so that raising
   ! one without the other stops the build instead of the science.
   !
   ! Against trc_tables' LMAX and not the TRC_LMAX preprocessor symbol: CMake
   ! defines TRC_LMAX and fpm does not, so referring to it here compiles under
   ! one build system and not the other. trc_tables already carries the
   ! `#ifndef TRC_LMAX / #define TRC_LMAX 2` default and exports the result as
   ! a parameter, which is the value both builds actually use.
   !
   ! Division rather than a zero-size array because a negative extent is legal
   ! Fortran and silently gives an empty array, whereas 1/0 in a constant
   ! expression is a hard error in gfortran, ifx and nvfortran alike.
   !
   integer, parameter, private :: ASSERT_1E_COVERS_ERI = &
      1/merge(1, 0, ONE_E_LMAX >= LMAX)
   ! The 1e, multipole and DF drivers and their per-item routines live in
   ! the four *_driver.inc files under src/inc, included into the kernel modules they dispatch
   ! into: called from here they crossed a module boundary, and stdpar only
   ! inlines within a translation unit.

contains

   !
   ! Cartesian multipole moments about `orig`, through the octupole.
   !
   ! (nao, nao, 39): 3 dipole, 9 quadrupole, 27 octupole, full tensors in
   ! libcint's order with the last index fastest -- what makeEFP's distributed
   ! multipoles consume.
   !
   ! A separate call from `trc_1e` even though the two share 1-D overlap
   ! tables. Fusing them was tried and put 39 components per class into one
   ! device routine, which NVVM chewed on for over thirty minutes at LMAX=2
   ! against two and a half apart. Both are computed once per geometry, so the
   ! shared tables were never worth it.
   !
   subroutine trc_multipoles(b, p, orig, mmat)
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: p
      real(dp), intent(in)  :: orig(3)
      real(dp), intent(out) :: mmat(b%nao, b%nao, TRC_NMULT)
      if (maxval(b%sh_l) > MULT_LMAX) error stop &
         "terco: multipole integrals stop at d (MULT_LMAX); this basis has a higher shell"
      call mult_driver(b%nshell, b%nao, b%maxnp, b%sh_l, b%sh_np, b%sh_ao, &
                       b%sh_e, b%sh_c, b%sh_r, p%npair, p%p_i, p%p_j, orig, mmat)
   end subroutine trc_multipoles


   subroutine trc_1e(b, p, smat, tmat, vmat)
      !! Overlap, kinetic and nuclear attraction, all three at once.
      type(trc_basis_t), intent(in)  :: b
      type(trc_pairlist_t), intent(in) :: p
      real(dp), intent(out) :: smat(b%nao, b%nao)
      real(dp), intent(out) :: tmat(b%nao, b%nao)
      real(dp), intent(out) :: vmat(b%nao, b%nao)
      call one_e_driver(b%nshell, b%nao, b%natm, b%maxnp, b%sh_l, b%sh_np, &
                        b%sh_ao, b%sh_e, b%sh_c, b%sh_r, b%at_z, b%at_r, &
                        p%npair, p%p_i, p%p_j, smat, tmat, vmat)
   end subroutine trc_1e

   !
   ! The one-electron GRADIENT: -W:dS/dR + D:dT/dR + D:dV/dR (basis
   ! derivative on both centres, Hellmann-Feynman on every nucleus),
   ! ADDED into `grad` -- a caller assembling the full analytic gradient
   ! calls this after eri%gradient (which starts `grad` from zero; see its
   ! own header) and before the trivial nuclear-repulsion term, same
   ! convention as trc_scf_grad%trc_scf_gradient. `wmat`/`dmat`/`grad` must
   ! already be device-resident, same as trc_1e's own smat/tmat/vmat.
   !
   subroutine trc_1e_gradient(b, p, wmat, dmat, grad)
      type(trc_basis_t), intent(in)  :: b
      type(trc_pairlist_t), intent(in) :: p
      real(dp), intent(in)    :: wmat(b%nao, b%nao)
      !! The energy-weighted density, W_mu,nu = sum_i n_i eps_i c_mu,i c_nu,i
      !! (trc_scf_grad.F90's scf_grad_energy_weighted_density builds it).
      real(dp), intent(in)    :: dmat(b%nao, b%nao)
      !! The ordinary density -- the TOTAL one (Da + Db for UHF): T and V
      !! are both linear in the density, so the spin split eri%gradient
      !! needs for its exchange term does not apply here.
      real(dp), intent(inout) :: grad(3, b%natm)
      call one_e_grad_driver(b%nshell, b%nao, b%natm, b%maxnp, b%sh_l, b%sh_np, &
                             b%sh_ao, b%sh_at, b%sh_e, b%sh_c, b%sh_r, b%at_z, b%at_r, &
                             p%npair, p%p_i, p%p_j, wmat, dmat, grad)
   end subroutine trc_1e_gradient

   !
   ! The polarisable continuum's solute-surface integrals. A surface point
   ! is a normalised s Gaussian of exponent zeta_k**2 centred at pts(:,k);
   ! its potential is erf(zeta_k r)/r, so these are the nuclear-attraction
   ! integrals with that operator in place of 1/r (gen_1e.py's --pcm
   ! kernels). Both Cartesian, both device-resident like trc_1e: b, p, pts,
   ! zeta and the matrix/vector arguments must already be present.
   !
   ! The pair list's primitive-pair data is what the kernels read, so `p`
   ! must have been built from `b`.
   !
   subroutine trc_pcm_potential(b, p, npts, pts, zeta, dmat, v)
      !! `v(k) = sum_{mu nu} D_mu,nu <mu| erf(zeta_k|r - pts_k|)/|r - pts_k| |nu>`
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: p
      integer, intent(in) :: npts
      real(dp), intent(in)  :: pts(3, npts), zeta(npts)
      real(dp), intent(in)  :: dmat(b%nao, b%nao)
      real(dp), intent(out) :: v(npts)
      if (maxval(b%sh_l) > PCM_LMAX) error stop &
         "terco: trc_pcm_potential: basis above PCM_LMAX (gen_1e.py --pcm --lmax)"
      call pcm_potential_driver(b%nshell, b%nao, b%sh_l, b%sh_ao, b%sh_r, p%npair, p%npp, &
                                p%p_i, p%p_j, p%pp_off, p%pp_n, p%pp_z, p%pp_k, p%pp_p, p%pp_a, &
                                npts, pts, zeta, dmat, v)
   end subroutine trc_pcm_potential

   subroutine trc_pcm_fock(b, p, npts, pts, zeta, q, vmat)
      !! `vmat(mu,nu) += sum_k q_k <mu| erf(zeta_k|r - pts_k|)/|r - pts_k| |nu>` -- ADDS
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: p
      integer, intent(in) :: npts
      real(dp), intent(in)    :: pts(3, npts), zeta(npts), q(npts)
      real(dp), intent(inout) :: vmat(b%nao, b%nao)
      if (maxval(b%sh_l) > PCM_LMAX) error stop &
         "terco: trc_pcm_fock: basis above PCM_LMAX (gen_1e.py --pcm --lmax)"
      call pcm_fock_driver(b%nshell, b%nao, b%sh_l, b%sh_ao, b%sh_r, p%npair, p%npp, &
                           p%p_i, p%p_j, p%pp_off, p%pp_n, p%pp_z, p%pp_k, p%pp_p, p%pp_a, &
                           npts, pts, zeta, q, vmat)
   end subroutine trc_pcm_fock

   !
   ! Explicit-shape dummies with the extents passed in, and the extents
   ! declared first. An assumed-shape dummy makes nvfortran walk a descriptor
   ! per launch; reaching through b%component inside the loop does the same.
   !

   subroutine trc_df_2c(aux, jmat)
      !! (P|Q) over an auxiliary basis.
      type(trc_basis_t), intent(in) :: aux
      real(dp), intent(out) :: jmat(aux%nao, aux%nao)
      ! The DF dispatchers have no case for a class above the ceilings they
      ! were generated for and would leave its block unwritten, silently.
      if (maxval(aux%sh_l) > DF_LMAX_AUX) error stop &
         "terco: trc_df_2c: auxiliary basis above DF_LMAX_AUX (gen_df.py --lmax-aux)"
      call rys_select_tables_init()
      call df_2c_driver(aux%nshell, aux%nao, aux%maxnp, aux%sh_l, aux%sh_np, &
                        aux%sh_ao, aux%sh_e, aux%sh_c, aux%sh_r, jmat)
   end subroutine trc_df_2c


   subroutine trc_df_3c(b, p, aux, tens)
      !! (mu nu|P), the three-index tensor density fitting is built on.
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: p
      real(dp), intent(out) :: tens(b%nao, b%nao, aux%nao)
      ! Same reason as trc_df_2c: no case means an unwritten block, not a stop.
      if (maxval(b%sh_l) > DF_LMAX) error stop &
         "terco: trc_df_3c: orbital basis above DF_LMAX (gen_df.py --lmax)"
      if (maxval(aux%sh_l) > DF_LMAX_AUX) error stop &
         "terco: trc_df_3c: auxiliary basis above DF_LMAX_AUX (gen_df.py --lmax-aux)"
      call rys_select_tables_init()
      call df_3c_driver(b%nshell, b%nao, b%maxnp, b%sh_l, b%sh_ao, b%sh_r, &
                        p%npair, p%npp, p%p_i, p%p_j, p%pp_off, p%pp_n, &
                        p%pp_z, p%pp_k, p%pp_p, p%pp_a, &
                        aux%nshell, aux%nao, aux%maxnp, aux%sh_l, aux%sh_np, &
                        aux%sh_ao, aux%sh_e, aux%sh_c, aux%sh_r, tens)
   end subroutine trc_df_3c

end module trc_integrals
