!
! The sizes of the ECP block kernels of trc_ecp_block.inc, for the module's
! specification part. The including module defines first:
!
!    ECP_BLOCK_LBRA   highest l of the bra shell the instance handles
!    ECP_BLOCK_LKET   highest l of the ket shell
!
! and uses dp, int64, ECP_LMAX and ECP_NRS. Two instances exist: trc_ecp (f
! bra, f ket: the energy matrix and every gradient task whose derivative bra
! is f or below) and trc_ecp_gbra (g bra, f ket: the derivative of an f
! shell). One instance sized for g would make every task carry g's scratch:
! the K1b columns grow by a third and the task kernel's per-thread stack, and
! with it the local memory the GPU reserves, by more than half.
!
   integer, parameter :: NPMAX = 16       !! primitives per orbital shell or ECP row
   real(dp), parameter :: ECP_NEGLIGIBLE = 1.0e-20_dp
   !! A contribution below this, to a matrix whose real entries run from
   !! ~1e-3 to a few hundred (measured: Xe/def2-SVP's U diagonal reaches
   !! 206), can never matter past the 1e-9 tolerance anything here is
   !! validated to. Used to stop the adaptive radial loop from chasing a
   !! RELATIVE 1e-12 match on a term that is already ABSOLUTELY negligible
   !! -- which a large separation (ecp_ints_ih_far's 25 bohr) or a tight
   !! core primitive pair reaches often, and which the relative test alone
   !! can fail to confirm within this port's ECP_LEVEL_MAX even though the
   !! term's actual contribution to U is already settled at zero.
   integer, parameter :: NFBRA = (ECP_BLOCK_LBRA + 1)*(ECP_BLOCK_LBRA + 2)/2   !! Cartesian bra functions
   integer, parameter :: NFKET = (ECP_BLOCK_LKET + 1)*(ECP_BLOCK_LKET + 2)/2   !! Cartesian ket functions
   integer, parameter :: D1MAX = ECP_BLOCK_LBRA + ECP_BLOCK_LKET + 1
   integer, parameter :: D2MAX = D1MAX*D1MAX
   integer, parameter :: D3MAX = D2MAX*D1MAX
   integer, parameter :: LC1BRA = ECP_BLOCK_LBRA + ECP_LMAX + 1
   integer, parameter :: LC1KET = ECP_BLOCK_LKET + ECP_LMAX + 1
   integer, parameter :: D2MAX2 = LC1BRA*LC1KET
   integer, parameter :: D3MAX2 = D1MAX*D2MAX2
   integer, parameter :: BUFMAX = NFBRA*(2*ECP_LMAX + 1)*LC1KET
   ! K1b's per-(shell, channel) global caches, sized to the instance's own
   ! ceilings and read at the caller's actual, smaller extent -- exactly
   ! `ecp_radial_type1_static_facs`'s, `ecp_radial_type2_facs_rad`'s and
   ! `ecp_type2_facs_ang`'s output extents at l = the ceiling.
   integer, parameter :: FAC1SIZE_BRA = NFBRA*(ECP_BLOCK_LBRA + 1)**3
   integer, parameter :: FAC1SIZE_KET = NFKET*(ECP_BLOCK_LKET + 1)**3
   integer, parameter :: FAC2SIZE_BRA = ECP_NRS*LC1BRA
   integer, parameter :: FAC2SIZE_KET = ECP_NRS*LC1KET
   integer, parameter :: ANG2SIZE_BRA = (ECP_BLOCK_LBRA + 1)*NFBRA*(2*ECP_LMAX + 1)*LC1BRA
   integer, parameter :: ANG2SIZE_KET = (ECP_BLOCK_LKET + 1)*NFKET*(2*ECP_LMAX + 1)*LC1KET
   integer(int64), parameter :: BRA_COL_BYTES = 8_int64*(FAC1SIZE_BRA + FAC2SIZE_BRA + ANG2SIZE_BRA)
   !! One (shell, channel group) column of the three K1b caches, bra sizes.
