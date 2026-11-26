module trc_ecp_grad_types
   !!
   !! The data of the ECP gradient (trc_ecp's THE GRADIENT), for the two
   !! modules whose kernels read it: trc_ecp and trc_ecp_gbra. Everything a
   !! kernel reads is mapped array by array by `trc_ecp`'s gradient driver.
   !!
   use pic_types, only: dp
   use trc_ecp_quadrature, only: ECP_NRS
   implicit none
   private
   public :: trc_ecp_channels_t, trc_ecp_deriv_basis_t, trc_ecp_grad_tasks_t
   public :: trc_ecp_k1b_cache_t, trc_ecp_chunk_t

   type :: trc_ecp_channels_t
      !! The ECP as the kernels see it: its (atom, l) channel groups, the atoms
      !! that carry them, the radial grid and each group's U(r) on it.
      integer :: ngroup = 0, natm_ecp = 0
      integer, allocatable :: grp_atom(:), grp_l(:), grp_row0(:), grp_row1(:)
      integer, allocatable :: atm_ecp(:), atm_grp0(:), atm_grp1(:)
      real(dp) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1)
      real(dp), allocatable :: ur(:, :)
   end type trc_ecp_channels_t

   type :: trc_ecp_deriv_basis_t
      !! The derivative basis: the real shells, then per real shell the
      !! raised one and (l > 0) the lowered one, as flat arrays (this file's
      !! header), with the fold arrays the task item contracts the density
      !! through.
      integer :: nreal = 0                !! real shells, 1 .. nreal
      integer :: nshell = 0               !! real and derivative shells
      integer :: natm = 0, nao = 0, maxnp = 0   !! the real basis's own
      real(dp), allocatable :: at_r(:, :)       !! its atom coordinates
      integer, allocatable :: l(:), np(:)
      real(dp), allocatable :: e(:, :), c(:, :), r(:, :)
      integer, allocatable :: atom(:)     !! the atom of derivative shell nreal + f
      integer :: nfold = 0
      !! Functions `parent` and `coef` hold per derivative shell: NFMAX (f),
      !! or 15 when a derivative shell is g.
      integer, allocatable :: parent(:, :, :)
      !! (component, lexicographic function of derivative shell f, f): the AO
      !! of the real basis that function is a derivative term of.
      real(dp), allocatable :: coef(:, :, :)
      !! Its coefficient: 1 (raised), -l_c (lowered), 0 where there is none.
      integer, allocatable :: ket_ao(:, :)
      !! (lexicographic function, real shell): its AO, in this build's order.
   end type trc_ecp_deriv_basis_t

   type :: trc_ecp_grad_tasks_t
      !! The gradient's kept tasks and the K1b cache columns they read. First
      !! the tasks whose bra is f or below, ECP atom by ECP atom (`first`),
      !! then those with a g derivative bra, again by atom (`gfirst`): the two
      !! task kernels.
      integer :: ntask = 0
      integer :: ng = 0                   !! of them, with a g bra (the last ng)
      integer, allocatable :: bra(:)      !! derivative bra shell, in the derivative basis
      integer, allocatable :: ket(:)      !! real ket shell
      integer, allocatable :: atm(:)      !! ECP atom, an index into the atm_* lists
      integer, allocatable :: first(:)    !! atom a's f-bra tasks are first(a) + 1 .. first(a + 1)
      integer, allocatable :: gfirst(:)   !! atom a's g-bra tasks are gfirst(a) + 1 .. gfirst(a + 1)
      integer, allocatable :: col(:, :)
      !! (shell, atom): the 0-based column, within the atom's own block of
      !! columns, of the shell's first channel group; -1 when no task of this
      !! atom reads the shell. A g shell's column is in the g-sized caches,
      !! every other shell's in the f-sized ones.
      integer, allocatable :: ncol(:)     !! f-sized cache columns atom a needs
      integer, allocatable :: ncol_g(:)   !! g-sized cache columns atom a needs
   end type trc_ecp_grad_tasks_t

   type :: trc_ecp_k1b_cache_t
      !! K1b's three caches, one (shell, channel group) per column, the column
      !! sizes those of the module that fills it.
      integer :: ncol = 0
      real(dp), allocatable :: fac1(:, :), fac2(:, :), ang2(:, :)
   end type trc_ecp_k1b_cache_t

   type :: trc_ecp_chunk_t
      !! One chunk of ECP atoms: its task slots, first the f-bra tasks 1..nt
      !! (tasks t0 + 1 on of the list) and then the g-bra ones nt + 1..nt + ntg
      !! (tasks tg0 + 1 on), and the g shells' cache columns.
      integer :: nt = 0, ntg = 0, t0 = 0, tg0 = 0
      integer, allocatable :: ci(:), cj(:)            !! each slot's bra and ket cache column
      integer, allocatable :: gcol_sh(:), gcol_g(:)   !! each g column's shell and channel group
      real(dp), allocatable :: w(:, :)                !! (3, nt + ntg): each task's three numbers
      integer, allocatable :: conv_fail(:)
   end type trc_ecp_chunk_t

end module trc_ecp_grad_types
