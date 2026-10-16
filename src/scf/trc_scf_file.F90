!
! A converged SCF, written to a file and read back, so a correlated step
! can start from it without repeating the SCF.
!
! WHAT IS IN THE FILE
! -------------------
! Binary stream, native byte order, 32-bit integers and 64-bit reals:
!
!   "TERCOSCF", version                         the format, checked first
!   natm, nshell, maxnp, nao, n, spherical,     the shape of everything below;
!   cart_order, nspin, nalpha, nbeta, charge,   `n` is the SCF's own space
!   multiplicity, converged, iterations,        (nao_sph when spherical)
!   n_lindep, pcm_npoints
!   functional (16 chars), pcm_model (8)        the method that made it
!   pcm_eps, energy, e_nuc, e_one, e_two,       the scalars of the result
!   e_xc, e_pcm, pcm_qtotal, nelec_grid
!   at_z(natm), at_r(3, natm)                   the molecule
!   sh_l, sh_np, sh_e, sh_c, sh_r               the basis, exactly as built
!   eps(n, nspin), cmo(n, n, nspin),            the orbitals and the density
!   dmat(n, n, nspin)
!
! WHAT A READ REFUSES
! -------------------
! Orbitals are only meaningful against the molecule, the basis and the AO
! convention they were made in, and nothing about a wrong file would show
! up as an error later: RI-MP2 would contract it and return a number. So
! everything that defines the AO space is stored and compared on the way
! in -- the atoms and positions, every shell's l, primitive count,
! exponents, coefficients and centre (that also catches --cart against a
! spherical run, and a different basis file of the same size), the
! Cartesian component order this build was made with, the charge and spin,
! the functional and the solvent. Nothing is converted: a mismatch is a
! refusal, with the first difference named.
!
! Run settings that change only how closely the result was converged (the
! screening threshold, the grid, the convergence criteria) are not
! compared; the file carries the converged flag and the iteration count,
! and the reader decides what to accept.
!
module trc_scf_file
   use, intrinsic :: iso_fortran_env, only: int32, real64
   use pic_types, only: dp
   use trc_api, only: trc_basis_t
   use trc_scf_driver, only: trc_scf_result_t
   use trc_cart_perm, only: trc_cart_order_id
   implicit none
   private

   public :: trc_scf_save, trc_scf_load

   character(len=8), parameter :: MAGIC = "TERCOSCF"
   integer(int32), parameter :: VERSION = 1
   integer, parameter :: N_SCALARS = 8   !! energy, e_nuc, e_one, e_two, e_xc, e_pcm, pcm_qtotal, nelec_grid

   ! The molecule is read from the same xyz the original run read, through
   ! the same conversion, so its positions come back bit-identical; the
   ! tolerances only allow for a different compiler's last bit.
   real(dp), parameter :: TOL_POS = 1.0e-8_dp    !! bohr
   real(dp), parameter :: TOL_REL = 1.0e-12_dp   !! exponents, coefficients, eps

contains

   !
   ! Write `res` (an SCF on `b`, the molecule `at_z`/`at_r` with `charge` and
   ! `mult`, `functional` and the solvent `pcm_model`/`pcm_eps`) to `path`.
   ! `msg` is empty on success.
   !
   subroutine trc_scf_save(path, b, at_z, at_r, charge, mult, nalpha, nbeta, functional, &
                           pcm_model, pcm_eps, res, msg)
      character(len=*), intent(in) :: path
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: at_z(:), at_r(:, :)
      integer, intent(in) :: charge, mult, nalpha, nbeta
      character(len=*), intent(in) :: functional, pcm_model
      real(dp), intent(in) :: pcm_eps
      type(trc_scf_result_t), intent(in) :: res
      character(len=*), intent(out) :: msg
      integer :: u, ios, n
      character(len=16) :: fname
      character(len=8) :: pname

      msg = ""
      if (.not. (allocated(res%cmo) .and. allocated(res%eps) .and. allocated(res%dmat))) then
         msg = "trc_scf_save: there is no SCF result to save"
         return
      end if
      n = size(res%cmo, 1)
      fname = functional
      pname = pcm_model
      open (newunit=u, file=path, access="stream", form="unformatted", status="replace", &
            action="write", iostat=ios)
      if (ios /= 0) then
         msg = "trc_scf_save: cannot write "//trim(path)
         return
      end if
      write (u, iostat=ios) MAGIC, VERSION, &
         i32(size(at_z)), i32(b%nshell), i32(b%maxnp), i32(b%nao), i32(n), &
         i32(merge(1, 0, b%spherical)), i32(trc_cart_order_id()), i32(res%nspin), &
         i32(nalpha), i32(nbeta), i32(charge), i32(mult), &
         i32(merge(1, 0, res%converged)), i32(res%iterations), i32(res%n_lindep), i32(res%pcm_npoints), &
         fname, pname, &
         real(pcm_eps, real64), real(res%energy, real64), real(res%e_nuc, real64), &
         real(res%e_one, real64), real(res%e_two, real64), real(res%e_xc, real64), &
         real(res%e_pcm, real64), real(res%pcm_qtotal, real64), real(res%nelec_grid, real64)
      if (ios == 0) write (u, iostat=ios) real(at_z, real64), real(at_r, real64)
      if (ios == 0) write (u, iostat=ios) int(b%sh_l, int32), int(b%sh_np, int32), &
         real(b%sh_e, real64), real(b%sh_c, real64), real(b%sh_r, real64)
      if (ios == 0) write (u, iostat=ios) real(res%eps, real64)
      if (ios == 0) write (u, iostat=ios) real(res%cmo, real64)
      if (ios == 0) write (u, iostat=ios) real(res%dmat, real64)
      close (u)
      if (ios /= 0) msg = "trc_scf_save: writing "//trim(path)//" failed"
   end subroutine trc_scf_save

   !
   ! Read `path` into `res`, after checking it was made for exactly this
   ! molecule (`at_z`, `at_r`, `charge`, `mult`), basis `b`, method and
   ! solvent. `nspin` is the spin count this run would use. On a refusal
   ! `msg` says why and `res` is untouched.
   !
   subroutine trc_scf_load(path, b, at_z, at_r, charge, mult, nspin, functional, &
                           pcm_model, pcm_eps, res, msg)
      character(len=*), intent(in) :: path
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: at_z(:), at_r(:, :)
      integer, intent(in) :: charge, mult, nspin
      character(len=*), intent(in) :: functional, pcm_model
      real(dp), intent(in) :: pcm_eps
      type(trc_scf_result_t), intent(inout) :: res
      character(len=*), intent(out) :: msg

      character(len=8) :: magic_in, f_pname
      character(len=16) :: f_fname
      integer(int32) :: ver, f_natm, f_nshell, f_maxnp, f_nao, f_n, f_sph, f_order, f_nspin, &
                        f_nalpha, f_nbeta, f_charge, f_mult, f_conv, f_iter, f_lindep, f_npcm
      real(real64) :: f_pcm_eps
      real(real64) :: sc(N_SCALARS)
      real(real64), allocatable :: f_z(:), f_r(:, :), f_e(:, :), f_c(:, :), f_sr(:, :)
      integer(int32), allocatable :: f_l(:), f_np(:)
      real(dp), allocatable :: eps(:, :), cmo(:, :, :), dmat(:, :, :)
      integer :: u, ios, n, natm
      logical :: there

      msg = ""
      inquire (file=path, exist=there)
      if (.not. there) then
         msg = "trc_scf_load: no file "//trim(path)
         return
      end if
      open (newunit=u, file=path, access="stream", form="unformatted", status="old", &
            action="read", iostat=ios)
      if (ios /= 0) then
         msg = "trc_scf_load: cannot open "//trim(path)
         return
      end if
      read (u, iostat=ios) magic_in, ver
      if (ios /= 0 .or. magic_in /= MAGIC) then
         call refuse(u, path, msg, trim(path)//" is not a terco SCF file")
         return
      end if
      if (ver /= VERSION) then
         call refuse(u, path, msg, trim(path)//" is format version "//itoa(int(ver))//", this build reads "// &
                     itoa(int(VERSION)))
         return
      end if
      read (u, iostat=ios) f_natm, f_nshell, f_maxnp, f_nao, f_n, f_sph, f_order, f_nspin, &
         f_nalpha, f_nbeta, f_charge, f_mult, f_conv, f_iter, f_lindep, f_npcm, &
         f_fname, f_pname, f_pcm_eps, sc
      if (ios /= 0) then
         call refuse(u, path, msg, trim(path)//" is truncated")
         return
      end if

      natm = size(at_z)
      n = merge(b%nao_sph, b%nao, b%spherical)
      if (f_natm /= natm) then
         call refuse(u, path, msg, "it holds "//itoa(int(f_natm))//" atoms, this molecule has "//itoa(natm))
         return
      end if
      if (f_charge /= charge .or. f_mult /= mult) then
         call refuse(u, path, msg, "it was made with charge "//itoa(int(f_charge))//" and multiplicity "// &
                     itoa(int(f_mult))//", this run has "//itoa(charge)//" and "//itoa(mult))
         return
      end if
      if (f_order /= trc_cart_order_id()) then
         call refuse(u, path, msg, "it was written by a build with a different Cartesian component order "// &
                     "(TERCO_CART_ORDER)")
         return
      end if
      if ((f_sph /= 0) .neqv. b%spherical) then
         call refuse(u, path, msg, "it was made in "//merge("spherical", "Cartesian", f_sph /= 0)// &
                     " AOs and this run uses "//merge("spherical", "Cartesian", b%spherical)// &
                     " ones (--cart, or --rimp2, which forces Cartesian)")
         return
      end if
      if (f_nshell /= b%nshell .or. f_maxnp /= b%maxnp .or. f_nao /= b%nao .or. f_n /= n) then
         call refuse(u, path, msg, "it was made in a different basis: "//itoa(int(f_nshell))//" shells and "// &
                     itoa(int(f_n))//" functions, this run has "//itoa(b%nshell)//" and "//itoa(n))
         return
      end if
      if (f_nspin /= nspin) then
         call refuse(u, path, msg, "it holds "//trim(merge("restricted  ", "unrestricted", f_nspin == 1))// &
                     " orbitals and this run is "//trim(merge("restricted  ", "unrestricted", nspin == 1)))
         return
      end if
      if (trim(f_fname) /= trim(functional)) then
         call refuse(u, path, msg, "it was made with "//method_name(f_fname)//" and this run asks for "// &
                     method_name(functional))
         return
      end if
      if (trim(f_pname) /= trim(pcm_model) .or. .not. close_rel(real(f_pcm_eps, dp), pcm_eps)) then
         call refuse(u, path, msg, "it was made with a different solvent (--pcm/--eps)")
         return
      end if

      allocate (f_z(natm), f_r(3, natm))
      read (u, iostat=ios) f_z, f_r
      if (ios /= 0) then
         call refuse(u, path, msg, trim(path)//" is truncated")
         return
      end if
      if (any(abs(f_z - at_z) > 0.0_dp)) then
         call refuse(u, path, msg, "its atoms are different elements, or in a different order")
         return
      end if
      if (any(abs(f_r - at_r) > TOL_POS)) then
         call refuse(u, path, msg, "its atoms are at different positions (largest difference "// &
                     rtoa(maxval(abs(f_r - at_r)))//" bohr)")
         return
      end if

      allocate (f_l(b%nshell), f_np(b%nshell))
      allocate (f_e(b%maxnp, b%nshell), f_c(b%maxnp, b%nshell), f_sr(3, b%nshell))
      read (u, iostat=ios) f_l, f_np, f_e, f_c, f_sr
      if (ios /= 0) then
         call refuse(u, path, msg, trim(path)//" is truncated")
         return
      end if
      if (any(f_l /= b%sh_l) .or. any(f_np /= b%sh_np) .or. &
          .not. all_close(f_e, b%sh_e) .or. .not. all_close(f_c, b%sh_c) .or. &
          any(abs(f_sr - b%sh_r) > TOL_POS)) then
         call refuse(u, path, msg, "it was made in a different basis with the same number of shells and functions")
         return
      end if

      allocate (eps(n, nspin), cmo(n, n, nspin), dmat(n, n, nspin))
      read (u, iostat=ios) eps
      if (ios == 0) read (u, iostat=ios) cmo
      if (ios == 0) read (u, iostat=ios) dmat
      if (ios /= 0) then
         call refuse(u, path, msg, trim(path)//" is truncated")
         return
      end if
      close (u)

      res%converged = f_conv /= 0
      res%iterations = int(f_iter)
      res%nspin = int(f_nspin)
      res%energy = sc(1)
      res%e_nuc = sc(2)
      res%e_one = sc(3)
      res%e_two = sc(4)
      res%e_xc = sc(5)
      res%e_pcm = sc(6)
      res%pcm_qtotal = sc(7)
      res%nelec_grid = sc(8)
      res%pcm_npoints = int(f_npcm)
      res%n_lindep = int(f_lindep)
      res%message = ""
      call move_alloc(eps, res%eps)
      call move_alloc(cmo, res%cmo)
      call move_alloc(dmat, res%dmat)
   end subroutine trc_scf_load

   !
   ! A module procedure rather than an internal one: see app/main.F90's
   ! terco_run for what a sibling internal procedure did to NVHPC's
   ! collective SCF.
   !
   subroutine refuse(u, path, msg, why)
      integer, intent(in) :: u
      character(len=*), intent(in) :: path, why
      character(len=*), intent(out) :: msg
      msg = "trc_scf_load: "//trim(path)//" does not fit this run: "//why
      close (u)
   end subroutine refuse

   pure function i32(i) result(res)
      integer, intent(in) :: i
      integer(int32) :: res
      res = int(i, int32)
   end function i32

   pure function close_rel(a, b) result(res)
      real(dp), intent(in) :: a, b
      logical :: res
      res = abs(a - b) <= TOL_REL*max(abs(a), abs(b))
   end function close_rel

   pure function all_close(a, b) result(res)
      real(real64), intent(in) :: a(:, :)
      real(dp), intent(in) :: b(:, :)
      logical :: res
      res = all(abs(a - b) <= TOL_REL*max(abs(a), abs(b)))
   end function all_close

   pure function method_name(f) result(s)
      character(len=*), intent(in) :: f
      character(len=:), allocatable :: s
      if (len_trim(f) == 0) then
         s = "Hartree-Fock"
      else
         s = trim(f)
      end if
   end function method_name

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=16) :: t
      write (t, "(i0)") i
      s = trim(t)
   end function itoa

   function rtoa(x) result(s)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=16) :: t
      write (t, "(es10.3)") x
      s = trim(adjustl(t))
   end function rtoa

end module trc_scf_file
