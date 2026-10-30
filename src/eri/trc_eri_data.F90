!
! The data the four-centre launchers share.
!
! Three containers, no procedures: the shells and their primitive pairs
! (hgp_build_pairs fills the pairs), the bin-interleaved copy of those pairs
! (bins_build_pair_quartets), and one launch over the bins (what to build and
! what came back). They exist so the launchers -- binkernel_fock_bins,
! grad_bins, the Schwarz bounds -- take a handful of objects instead of forty
! arrays.
!
! NONE OF THESE IS EVER MAPPED WHOLE. Their arrays go to the device one by
! one (`enter data copyin(pairs%p, ...)` in trc_eri), and a launcher hands
! them to its kernels as explicit-shape actual arguments, so a `do concurrent`
! body never reads a derived-type component. Growing a type that IS mapped
! whole has crashed nvfortran before (see trc_bins' bins_build_pair_quartets);
! that is why these are separate types rather than more fields of
! `pair_bins_t`.
!
module trc_eri_data
   use pic_types, only: dp, int64
   implicit none
   private

   type, public :: trc_shell_pairs_t
      !! The shells the integrals run over, and every kept primitive pair of
      !! every shell pair: what the Fock, gradient and Schwarz kernels read
      !! about the basis. The caller fills the shell part; hgp_build_pairs
      !! fills the pairs.
      integer :: nbas = 0
      !! Shells.
      integer :: natm = 0
      !! Atoms: the range of `sh_at`, and the gradient's second extent.
      integer, allocatable :: sh_l(:)
      !! (nbas) angular momentum.
      integer, allocatable :: ao_off(:)
      !! (nbas) first Cartesian AO of the shell, 1-based.
      integer, allocatable :: sh_at(:)
      !! (nbas) atom of the shell; only the gradient reads it.
      integer :: npp = 0
      !! Kept primitive pairs.
      integer, allocatable :: off(:), n(:)
      !! (nbas*nbas) first primitive pair - 1 and how many, rectangular key
      !! (i-1)*nbas + j.
      real(dp), allocatable :: p(:)
      !! (npp) zeta = alpha + beta.
      real(dp), allocatable :: r(:, :), ra(:, :), rb(:, :)
      !! (npp, 3) the product centre P and the two shell centres A and B.
      real(dp), allocatable :: c(:)
      !! (npp) the combined coefficient, K_ab = exp(-mu |AB|^2) folded in.
      real(dp), allocatable :: a(:)
      !! (npp) the i-side primitive's own exponent; beta is p - a. Only the
      !! Rys gradient kernels read it.
   end type trc_shell_pairs_t

   type, public :: trc_pair_quartets_t
      !! The primitive pairs again, laid out bin by bin and interleaved by
      !! primitive index, so a warp's threads read neighbouring slots (see
      !! bins_build_pair_quartets). Slot ix(b, k, s) = bin_pbase(b) +
      !! (k-1)*n_b + s for bin b, primitive k, pair s.
      integer :: nqp = 0
      !! Slots.
      integer, allocatable :: bin_pbase(:), bin_k(:)
      !! (nbin) first slot - 1 of a bin, and its largest kept count.
      real(dp), allocatable :: z(:), c(:), a(:)
      !! (nqp) zeta, coefficient and the i-side exponent, as in the pairs.
      real(dp), allocatable :: px(:), py(:), pz(:)
      !! (nqp) the product centre.
      real(dp), allocatable :: ax(:), ay(:), az(:)
      !! (npair) per bin slot: the A centre.
      real(dp), allocatable :: abx(:), aby(:), abz(:)
      !! (npair) per bin slot: A - B.
      integer, allocatable :: n(:)
      !! (npair) per bin slot: the pair's real kept count.
   end type trc_pair_quartets_t

   type, public :: trc_bin_launch_t
      !! One launch over the bins: the operator, the operands' shape, this
      !! rank's share, and what the launch reports back.
      real(dp) :: thresh = 1.0e-10_dp
      !! Integral (Schwarz) threshold.
      real(dp) :: jfac = 1.0_dp, kfac = 1.0_dp
      !! Coulomb and exchange scalings; 1 and 1 is Hartree-Fock.
      real(dp) :: omega = 0.0_dp
      !! Range-separation parameter, bohr^-1; zero is the Coulomb operator.
      integer :: ndens = 1
      !! Densities in the batch: the first extent of dmat, jmat, kmat.
      integer :: nao = 0
      !! Functions: the other two extents.
      integer :: rank = 0, nranks = 1
      !! This rank takes every nranks-th work item from `rank`; the caller
      !! sums the ranks.
      logical :: use_density = .false.
      !! Bin-level density screen (not armed; see density_reject in trc_bins).
      logical :: nosym = .false.
      !! The density is not symmetric: route through the enumerated kernel.
      logical :: count_kept = .false.
      !! Also count the quartets that survive both screens, into `nkept`.
      !! A full extra pass, so it is opt-in.
      integer :: nlaunch = 0
      !! Out: kernel launches.
      integer(int64) :: nwork = 0
      !! Out: work items (quartets before the in-kernel screens).
      integer(int64) :: nkept = 0
      !! Out: quartets that survived, when `count_kept` asked for it.
   end type trc_bin_launch_t

end module trc_eri_data
