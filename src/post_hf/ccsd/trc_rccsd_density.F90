!
! The unrelaxed correlation densities of closed-shell, density-fitted CCSD, from the converged amplitudes and
! Lambda multipliers, in the active MO space: the one-particle density and the three-index two-particle density.
! What pyscf/cc/ccsd_rdm.py `_gamma1_intermediates` and `_gamma2_outcore` define (REAL orbitals, spatial, no
! reference part), without ever forming an o v^3 or v^4 array. Orbital relaxation, the gradient and (T) are not
! here.
!
! WHAT IS COMPUTED. With i, j active occupied, a, b virtual, p, q either, and gamma and Gamma the correlation
! one- and two-particle densities symmetrised so that E_corr = sum_pq gamma(pq) f(pq) + 1/2 sum_pqrs Gamma(pqrs)
! (pq|rs), Gamma carrying the eightfold symmetry of the integrals:
!
!     gamma_oo(i,j), gamma_ov(i,a), gamma_vv(a,b)          gamma_ov is also the vo block transposed
!     g3(pq,P) = sum_rs Gamma(pqrs) b(rs,P)
!
! with b the three-index factors of trc_rccsd_ints, so (pq|rs) = sum_P b(pq,P) b(rs,P). g3 comes in the layout of
! those factors, first orbital fastest: g3_oo(i + (j-1) no, P), g3_ov(i + (a-1) no, P), g3_vv(a + (b-1) nv, P);
! g3_oo and g3_vv are symmetric in their pair, and g3_vo is g3_ov transposed. Hence
!
!     E_corr = sum_pp gamma(pp) eps(p) + 1/2 [ g3_oo . b_oo + 2 g3_ov . b_ov + g3_vv . b_vv ]
!
! at converged t and l. The auxiliary index P is that of `ints%b_*` (the Cholesky factor of the metric), so a
! consumer that wants densities over the raw auxiliary functions transforms with that factor.
!
! NOTHING ELSE OF THE TWO-PARTICLE DENSITY IS NEEDED LATER. The orbital Lagrangian of a gradient is
! X(t,p) = ... + sum_qrs (tq|rs) Gamma(pqrs) = ... + sum_qP b(tq,P) g3(pq,P) for any t (frozen orbitals
! included; p, q active), and the three- and two-centre derivative densities are g3 transformed to the AO and
! auxiliary bases and contracted with the fit. Gamma itself is never required, and is not kept.
!
! BLOCKS. pyscf's names, with Gamma assembled as pyscf's `_make_rdm2` does and then symmetrised: dovov, doooo,
! doovv, dovvo and dooov are formed (they are o^2 v^2 at the largest) and contracted with b; dovvv and dvvvv
! are not. Their contractions with b are re-ordered:
!   - the six terms of dovvv that are a product through one occupied index end in GEMMs of an o^2 v^2 array with
!     a three-index one;
!   - dvvvv, and the one term of dovvv that contains it, go through
!         C'(s,q,r,d) = sum_ij [tau(ij,sq) - 1/2 tau(ij,qs)] l2(ij,rd),
!     built for a slab of `slab_batch` virtuals d at a time (v^3 slab_batch doubles, never v^4) and contracted
!     with b_vv and with Bt(a,d;P) = sum_i t1(i,d) b_ov(i,a,P): with E(sr,qd) = C'(s,q,r,d),
!         y(d,q;P) = sum_sr E(sr,qd) b_vv(sr,P),      g3_vv(s r) -= 4 sum_ad E(sr,ad) Bt(a,d;P),
!         g3_vv(q,d) += 4 y(d,q),                      g3_ov(i,a) -= 2 sum_d t1(i,d) y(d,a).
!
! INPUTS. t2 and l2 must be pair-symmetric, t(j,i,b,a) = t(i,j,a,b), as the CCSD and Lambda solvers leave them
! (not checked here; trc_rccsd_lambda checks t2). The integrals enter through b_oo, b_ov and b_vv only.
!
! MEMORY (o active occupied, v virtual, x auxiliary, doubles, host and device each). The result is
! (o^2 + o v + v^2) x doubles of g3. During the build: 7 o^2 v^2, 2 o^3 v, 3 o^4, o v x, o^2 x, v^2 x, the
! contraction plans' scratch (trc_contract; `work_doubles`) and the slab, 2 v^3 slab_batch + v slab_batch x.
! `slab_batch` is `options%slab_batch`, else as many virtuals as SLAB_BUDGET holds, never more than a quarter of
! the virtual space; the environment variable TERCO_CCDENS_D_CHUNK=<k> caps it (TESTING: forces several slabs on
! a small molecule). Everything is sized and checked against the free device memory before anything is
! allocated: a refusal states the bytes and leaves the result empty.
!
! FLOPS. The slabs are 2 o^2 v^4 + 4 v^4 x, the rest o^3 v^3 and o^2 v^2 x terms.
!
! RESIDENCY. The result arrays exist on host and device; the build writes the device copies, `pull` brings them
! back. The build reads the device copies of the amplitudes, multipliers and factors.
!
! Failures are reported in `status`; none stops. No MPI call: one rank.
#include "inc/trc_dc_index.inc"

module trc_rccsd_density
   use pic_types, only: dp, int32, int64
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use pic_timer, only: timer_type
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_driver, only: trc_rccsd_state_t
   use trc_rccsd_lambda, only: trc_rccsd_lambda_state_t
   use trc_rccsd_steps, only: rccsd_steps_gemm_flops, rccsd_steps_open, rccsd_steps_close
   use trc_contract, only: trc_contract_plan_t, contract_permute
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_step_table, only: step_table_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_device_memory, only: device_memory_set_cap, device_memory_free_bytes, device_memory_fits
   implicit none
   private

   public :: trc_rccsd_density_t, trc_rccsd_density_options_t
   public :: trc_rccsd_density_solve, rccsd_density_build, rccsd_density_bytes

   integer, parameter, public :: RCCSD_DENSITY_ST_ONE = 1    !! gamma, tau and the two o^2 v^2 products of l2 and t2
   integer, parameter, public :: RCCSD_DENSITY_ST_OOOO = 2   !! doooo and dooov
   integer, parameter, public :: RCCSD_DENSITY_ST_OVVO = 3   !! doovv and dovvo
   integer, parameter, public :: RCCSD_DENSITY_ST_OVOV = 4   !! dovov
   integer, parameter, public :: RCCSD_DENSITY_ST_OVVV = 5   !! the factorised terms of dovvv
   integer, parameter, public :: RCCSD_DENSITY_ST_VVVV = 6   !! the slabs: dvvvv and the dovvv term that holds it
   integer, parameter, public :: RCCSD_DENSITY_NSTEP = 6
   character(len=24), parameter, public :: RCCSD_DENSITY_STEP_NAME(RCCSD_DENSITY_NSTEP) = [character(len=24) :: &
      "gamma, tau, l2 t2", "oooo, ooov blocks", "oovv, ovvo blocks", "ovov block", "ovvv terms", "vvvv slabs"]
   !! rows of `trc_rccsd_density_t%steps`

   integer(int64), parameter :: SLAB_BUDGET = 1073741824_int64  !! bytes of a slab when `options%slab_batch` is 0
   integer, parameter :: SLAB_FRACTION = 4  !! a slab holds at most nv/SLAB_FRACTION virtuals (at least one)
   integer, parameter :: NSQUARE = 7        !! o^2 v^2 work arrays

   ! the contraction plans, by what they form
   integer, parameter :: P_DOO = 1, P_DVV = 2, P_GOO = 3, P_GVV = 4, P_XTV = 5, P_GOV_A = 6, P_GOV_B = 7, &
                         P_GOV_C = 8, P_GOV_D = 9, P_PA = 10, P_PB_A = 11, P_PB_B = 12, P_MIA = 13, P_MAB = 14, &
                         P_MIJ = 15, P_GOOOO = 16, P_G3_A = 17, P_G3_B = 18, P_G3_C = 19, P_G3_D = 20, P_G3_E = 21, &
                         P_K_OO = 22, P_OOOV_OO = 23, P_OOOV_OV = 24, P_U = 25, P_UP = 26, P_U_T = 27, &
                         P_OOVV_OO = 28, P_OOVV_VV = 29, P_OOVV_A = 30, P_OOVV_B = 31, P_OOVV_C = 32, &
                         P_OOVV_D = 33, P_RING_A = 34, P_RING_B = 35, P_RING_C = 36, P_M_B = 37, P_H_T = 38, &
                         P_H_L = 39, P_K_H = 40, P_K_B = 41, P_Z_T = 42, P_Z_L = 43, P_R = 44, P_R_T = 45, &
                         P_MOO = 46, P_P2_MOO = 47, P_N = 48, P_N_T = 49, P_BT = 50, P_M_BT = 51, NPLAN = 51

   type :: trc_rccsd_density_options_t
      !! Plain values only
      integer :: slab_batch = 0  !! virtuals per vvvv slab; 0 chooses (module header, MEMORY)
      integer(int64) :: extra_work_limit = -1
      !! the contraction plans' `extra_work_limit` (doubles); negative keeps their default
      integer(int64) :: device_free_cap = 0_int64
      !! TEST HOOK, as trc_rccsd_options_t's: above 0, the free device memory every size check of the build sees
      logical :: verbose = .false.  !! the step table
   end type trc_rccsd_density_options_t

   type :: trc_rccsd_density_t
      !! The densities on host and device (module header), valid when `built`
      type(error_t) :: status  !! set by the build: a refusal leaves its message here and the type empty
      integer :: no = 0
      integer :: nv = 0
      integer :: naux = 0
      logical :: built = .false.
      integer :: slab_batch = 0           !! virtuals per slab the last build took
      integer(int64) :: work_doubles = 0  !! doubles of contraction scratch of the last build
      integer(int64) :: temp_bytes = 0    !! bytes of every temporary of the last build
      integer(int64) :: flops(RCCSD_DENSITY_NSTEP) = 0  !! 2 m n k of the GEMMs of the last build, per step
      type(step_table_t) :: steps         !! rows RCCSD_DENSITY_STEP_NAME of the last build
      real(dp), allocatable :: gamma_oo(:, :)  !! gamma(i,j)
      real(dp), allocatable :: gamma_ov(:, :)  !! gamma(i,a) = gamma(a,i)
      real(dp), allocatable :: gamma_vv(:, :)  !! gamma(a,b)
      real(dp), allocatable :: g3_oo(:, :)     !! g3(ij,P), no^2 by naux
      real(dp), allocatable :: g3_ov(:, :)     !! g3(ia,P), no nv by naux
      real(dp), allocatable :: g3_vv(:, :)     !! g3(ab,P), nv^2 by naux
   contains
      procedure :: release => density_release
      procedure :: pull => density_pull
   end type trc_rccsd_density_t

   type :: work_t
      !! The plans and temporaries of one build; names of the arrays are of their first use
      type(trc_contract_plan_t) :: plan(NPLAN)
      type(error_t) :: failure        !! the first contraction that refused
      integer(int64) :: flops = 0     !! of the open step
      real(dp), allocatable :: work(:)    !! the plans' scratch
      real(dp), allocatable :: tau(:)     !! o^2 v^2, tau(i,j,a,b); the slabs turn it into tau - tau swapped/2
      real(dp), allocatable :: pa(:)      !! o^2 v^2, pa(i,a,j,b) = pyscf pvOOv(a,i,j,b)
      real(dp), allocatable :: pb(:)      !! o^2 v^2, pb(i,a,j,b) = pyscf pvoOV(a,i,j,b)
      real(dp), allocatable :: m(:)       !! o^2 v^2, the matrix that multiplies b_ov in g3_ov
      real(dp), allocatable :: w1(:), w2(:), w3(:)  !! o^2 v^2
      real(dp), allocatable :: s1(:), s2(:)         !! o^3 v
      real(dp), allocatable :: goooo(:), d4(:), k4(:)  !! o^4
      real(dp), allocatable :: hov(:)     !! o v naux
      real(dp), allocatable :: hoo(:)     !! o^2 naux
      real(dp), allocatable :: hvv(:)     !! v^2 naux: N, then Bt(a,d;P)
      real(dp), allocatable :: u(:)       !! naux
      real(dp), allocatable :: cbuf(:)    !! v^3 slab_batch, C'
      real(dp), allocatable :: ebuf(:)    !! v^3 slab_batch, E
      real(dp), allocatable :: ybuf(:)    !! v slab_batch naux, y
      real(dp), allocatable :: doo(:, :), dvv(:, :)  !! sum l2 theta over (k,a,b) and over (i,j,c)
      real(dp), allocatable :: xtv(:, :), moo(:, :), mij(:, :)     !! o^2
      real(dp), allocatable :: mvv(:, :), mab(:, :), m2(:, :)      !! v^2
      real(dp), allocatable :: mia(:, :)                           !! o v
      real(dp), allocatable :: d1oo(:, :), d1vv(:, :)  !! pyscf doo, dvv
   end type work_t

contains

   subroutine density_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_density_t), intent(inout) :: this

      if (allocated(this%gamma_oo)) then
         !$acc exit data delete(this%gamma_oo, this%gamma_ov, this%gamma_vv, this%g3_oo, this%g3_ov, this%g3_vv)
         deallocate (this%gamma_oo, this%gamma_ov, this%gamma_vv, this%g3_oo, this%g3_ov, this%g3_vv)
      end if
      this%no = 0
      this%nv = 0
      this%naux = 0
      this%built = .false.
      this%slab_batch = 0
      this%work_doubles = 0
      this%temp_bytes = 0
      this%flops = 0_int64
   end subroutine density_release

   subroutine density_pull(this)
      !! Device copy to host
      class(trc_rccsd_density_t), intent(inout) :: this

      if (.not. allocated(this%gamma_oo)) return
      !$acc update self(this%gamma_oo, this%gamma_ov, this%gamma_vv, this%g3_oo, this%g3_ov, this%g3_vv)
   end subroutine density_pull

   subroutine trc_rccsd_density_solve(cc, lam, options, la, den)
      !! The densities of the converged RI-CCSD `cc` and its converged Lambda solution `lam` (module header).
      !! Refuses, with the message in `den%status`, a CCSD or Lambda that was refused, did not converge, was
      !! released or ran over several ranks. `cc` and `lam` are read, their device copies current.
      type(trc_rccsd_state_t), intent(in) :: cc
      type(trc_rccsd_lambda_state_t), intent(in) :: lam
      type(trc_rccsd_density_options_t), intent(in) :: options
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_density_t), intent(inout) :: den

      call den%release()
      call den%status%clear()
      if (cc%result%refused .or. .not. cc%result%converged) then
         call refuse(den, "trc_rccsd_density: the CCSD was refused or did not converge")
      else if (allocated(cc%rank_seconds)) then
         call refuse(den, "trc_rccsd_density: the CCSD ran over several ranks: the densities are one rank only")
      else if (.not. allocated(cc%ints%b_ov) .or. .not. allocated(cc%amps%t2)) then
         call refuse(den, "trc_rccsd_density: the CCSD state holds no integrals or amplitudes (released?)")
      else if (lam%result%refused .or. .not. lam%result%converged) then
         call refuse(den, "trc_rccsd_density: the Lambda equations were refused or did not converge")
      else if (.not. allocated(lam%l%t2)) then
         call refuse(den, "trc_rccsd_density: the Lambda state holds no multipliers (released?)")
      else
         call rccsd_density_build(cc%ints, cc%amps, lam%l, options, la, den)
      end if
   end subroutine trc_rccsd_density_solve

   subroutine refuse(den, text)
      !! Leaves `den` empty with `text` in its status
      type(trc_rccsd_density_t), intent(inout) :: den
      character(len=*), intent(in) :: text

      call den%release()
      call den%status%set(ERROR_VALIDATION, text)
      call global_logger%warning(text)
   end subroutine refuse

   subroutine rccsd_density_build(ints, amps, l, options, la, den)
      !! The densities (module header) of the amplitudes `amps` and the multipliers `l` (l1 in `%t1`, l2 in
      !! `%t2`), whatever they are: nothing checks convergence. t2 and l2 pair-symmetric; all device copies
      !! current. Of `ints` only (no, nv, naux) and b_oo, b_ov, b_vv are read. A refusal is in `den%status` and
      !! leaves `den` empty.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(trc_rccsd_density_options_t), intent(in) :: options
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_density_t), intent(inout) :: den

      type(work_t) :: wk

      call den%release()
      call den%status%clear()
      if (.not. inputs_ok(ints, amps, l, options, den)) return
      call device_memory_set_cap(options%device_free_cap)
      call plan_all(ints%no, ints%nv, ints%naux, options%extra_work_limit, wk%plan, den%status)
      if (.not. den%status%has_error()) call allocate_all(options, wk, den, ints%no, ints%nv, ints%naux)
      call device_memory_set_cap(0_int64)
      if (den%status%has_error()) then
         call global_logger%warning(den%status%get_message())
         return
      end if
      call den%steps%init("rccsd density", RCCSD_DENSITY_STEP_NAME)
      call run_steps(ints, amps, l, wk, den, la)
      call free_work(wk)
      if (wk%failure%has_error()) then
         call refuse(den, "trc_rccsd_density: "//wk%failure%get_message())
         return
      end if
      den%built = .true.
      if (options%verbose) call den%steps%report(sum(den%steps%t), 0, 1)
   end subroutine rccsd_density_build

   function inputs_ok(ints, amps, l, options, den) result(ok)
      !! True when the factors, amplitudes and multipliers are there with the sizes of `ints`
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(trc_rccsd_density_options_t), intent(in) :: options
      type(trc_rccsd_density_t), intent(inout) :: den
      logical :: ok

      ok = .false.
      if (.not. allocated(ints%b_oo) .or. .not. allocated(ints%b_ov) .or. .not. allocated(ints%b_vv)) then
         call refuse(den, "trc_rccsd_density: the three-index factors are not built")
         return
      end if
      if (ints%no < 1 .or. ints%nv < 1 .or. ints%naux < 1) then
         call refuse(den, "trc_rccsd_density: no active occupied, virtual or auxiliary functions")
         return
      end if
      if (.not. allocated(amps%t1) .or. .not. allocated(amps%t2) .or. .not. allocated(l%t1) .or. &
          .not. allocated(l%t2)) then
         call refuse(den, "trc_rccsd_density: the amplitudes or the multipliers are not initialised")
         return
      end if
      if (amps%no /= ints%no .or. amps%nv /= ints%nv .or. l%no /= ints%no .or. l%nv /= ints%nv) then
         call refuse(den, "trc_rccsd_density: (no, nv) differ between integrals, amplitudes and multipliers")
         return
      end if
      if (options%slab_batch < 0) then
         call refuse(den, "trc_rccsd_density: a negative slab_batch")
         return
      end if
      ok = .true.
   end function inputs_ok

   subroutine rccsd_density_bytes(no, nv, naux, slab_batch, result_bytes, temporary_bytes, status)
      !! Device bytes of a build at these sizes with `slab_batch` virtuals per slab: of the result, and of every
      !! temporary (the host mirrors are as much again). Nothing is allocated. `status` says when a size is
      !! beyond what the build handles.
      integer, intent(in) :: no, nv, naux, slab_batch
      integer(int64), intent(out) :: result_bytes, temporary_bytes
      type(error_t), intent(inout) :: status

      type(trc_contract_plan_t) :: plan(NPLAN)
      integer(int64) :: work

      result_bytes = 0
      temporary_bytes = 0
      call plan_all(no, nv, naux, -1_int64, plan, status)
      if (status%has_error()) return
      work = largest_work(plan)
      result_bytes = result_size(no, nv, naux)
      temporary_bytes = fixed_size(no, nv, naux, work) + int(max(slab_batch, 1), int64)*slab_size(nv, naux)
   end subroutine rccsd_density_bytes

   pure function result_size(no, nv, naux) result(bytes)
      !! Bytes of gamma and g3
      integer, intent(in) :: no, nv, naux
      integer(int64) :: bytes
      integer(int64) :: o, v, x

      o = no
      v = nv
      x = naux
      bytes = 8_int64*(o*o + o*v + v*v)*(x + 1_int64)
   end function result_size

   pure function fixed_size(no, nv, naux, work) result(bytes)
      !! Bytes of the temporaries that do not depend on the slab
      integer, intent(in) :: no, nv, naux
      integer(int64), intent(in) :: work
      integer(int64) :: bytes
      integer(int64) :: o, v, x

      o = no
      v = nv
      x = naux
      bytes = 8_int64*(work + NSQUARE*o*o*v*v + 2_int64*o*o*o*v + 3_int64*o*o*o*o + (o*v + o*o + v*v)*x + x &
                       + 5_int64*(o*o + v*v) + o*v)
   end function fixed_size

   pure function slab_size(nv, naux) result(bytes)
      !! Bytes of one virtual of a slab
      integer, intent(in) :: nv, naux
      integer(int64) :: bytes
      integer(int64) :: v

      v = nv
      bytes = 8_int64*(2_int64*v*v*v + v*naux)
   end function slab_size

   function largest_work(plan) result(work)
      !! The scratch the plans need, at least one double
      type(trc_contract_plan_t), intent(in) :: plan(NPLAN)
      integer(int64) :: work
      integer :: k

      work = 1
      do k = 1, NPLAN
         work = max(work, plan(k)%work_size)
      end do
   end function largest_work

   subroutine plan_one(plan, spec, extent_a, extent_b, error)
      !! Builds `plan`, keeping the first error
      type(trc_contract_plan_t), intent(inout) :: plan
      character(len=*), intent(in) :: spec
      integer, intent(in) :: extent_a(:), extent_b(:)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      call plan%build(spec, extent_a, extent_b, error)
   end subroutine plan_one

   subroutine plan_all(no, nv, naux, extra_work_limit, plan, error)
      !! Every contraction of the build outside the slabs, before anything is allocated
      integer, intent(in) :: no, nv, naux
      integer(int64), intent(in) :: extra_work_limit
      type(trc_contract_plan_t), intent(inout) :: plan(NPLAN)
      type(error_t), intent(inout) :: error

      integer(int64) :: o, v, x, big
      integer :: k

      o = no
      v = nv
      x = naux
      ! every array below is indexed with default integers
      big = max(o*o*v*v, o*o*o*max(o, v), v*v*x, o*v*x, o*o*x, v*v*v)
      if (big > int(huge(0_int32), int64)) then
         call error%set(ERROR_VALIDATION, "trc_rccsd_density: a block has more than 2^31-1 elements (o "// &
                        to_string(no)//", v "//to_string(nv)//", naux "//to_string(naux)//")")
         return
      end if
      if (extra_work_limit >= 0) then
         do k = 1, NPLAN
            plan(k)%extra_work_limit = extra_work_limit
         end do
      end if
      call plan_one_particle(no, nv, plan, error)
      call plan_blocks(no, nv, naux, plan, error)
      call plan_three_index(no, nv, naux, plan, error)
   end subroutine plan_all

   subroutine plan_one_particle(no, nv, plan, error)
      !! The plans of gamma, of the two products of l2 and t2 and of the small intermediates
      integer, intent(in) :: no, nv
      type(trc_contract_plan_t), intent(inout) :: plan(NPLAN)
      type(error_t), intent(inout) :: error

      integer :: t1(2), t2(4), oo(2), vv(2)

      t1 = [no, nv]
      t2 = [no, no, nv, nv]
      oo = [no, no]
      vv = [nv, nv]
      call plan_one(plan(P_DOO), "ikab,jkab->ij", t2, t2, error)
      call plan_one(plan(P_DVV), "jica,jicb->ab", t2, t2, error)
      call plan_one(plan(P_GOO), "ia,ja->ij", t1, t1, error)
      call plan_one(plan(P_GVV), "ia,ib->ab", t1, t1, error)
      call plan_one(plan(P_XTV), "ie,me->im", t1, t1, error)
      call plan_one(plan(P_GOV_A), "im,ma->ia", oo, t1, error)
      call plan_one(plan(P_GOV_B), "imae,me->ia", t2, t1, error)
      call plan_one(plan(P_GOV_C), "mi,ma->ia", oo, t1, error)
      call plan_one(plan(P_GOV_D), "ie,ae->ia", t1, vv, error)
      call plan_one(plan(P_PA), "ikca,jkcb->iajb", t2, t2, error)
      call plan_one(plan(P_PB_A), "ikca,jkbc->iajb", t2, t2, error)
      call plan_one(plan(P_PB_B), "ikac,jkbc->iajb", t2, t2, error)
      call plan_one(plan(P_MIA), "ikac,kc->ia", t2, t1, error)
      call plan_one(plan(P_MAB), "kc,kb->cb", t1, t1, error)
      call plan_one(plan(P_MIJ), "jc,kc->jk", t1, t1, error)
   end subroutine plan_one_particle

   subroutine plan_blocks(no, nv, naux, plan, error)
      !! The plans of the o^2 v^2 and smaller blocks of the two-particle density
      integer, intent(in) :: no, nv, naux
      type(trc_contract_plan_t), intent(inout) :: plan(NPLAN)
      type(error_t), intent(inout) :: error

      integer :: t1(2), t2(4), ovov(4), o4(4), ooov(4), oovo(4), boo(3), bov(3), bvv(3), oo(2), vv(2)

      t1 = [no, nv]
      t2 = [no, no, nv, nv]
      ovov = [no, nv, no, nv]
      o4 = [no, no, no, no]
      ooov = [no, no, no, nv]
      oovo = [no, no, nv, no]
      oo = [no, no]
      vv = [nv, nv]
      boo = [no, no, naux]
      bov = [no, nv, naux]
      bvv = [nv, nv, naux]
      call plan_one(plan(P_GOOOO), "ijab,klab->ijkl", t2, t2, error)
      call plan_one(plan(P_G3_A), "kc,icja->jkia", t1, ovov, error)
      call plan_one(plan(P_G3_B), "jc,icka->jkia", t1, ovov, error)
      call plan_one(plan(P_G3_C), "jkil,la->jkia", o4, t1, error)
      call plan_one(plan(P_G3_D), "ib,jkba->jkia", t1, t2, error)
      call plan_one(plan(P_G3_E), "jkba,ib->jkia", t2, t1, error)
      call plan_one(plan(P_K_OO), "pqrs,rsP->pqP", o4, boo, error)
      call plan_one(plan(P_OOOV_OO), "pqra,raP->pqP", ooov, bov, error)
      call plan_one(plan(P_OOOV_OV), "jkia,jkP->iaP", ooov, boo, error)
      call plan_one(plan(P_U), "ikac,jc->ikaj", t2, t1, error)
      call plan_one(plan(P_UP), "ikca,jc->ikaj", t2, t1, error)
      call plan_one(plan(P_U_T), "ikaj,kb->iajb", oovo, t1, error)
      call plan_one(plan(P_OOVV_OO), "ijab,abP->ijP", t2, bvv, error)
      call plan_one(plan(P_OOVV_VV), "ijab,ijP->abP", t2, boo, error)
      call plan_one(plan(P_OOVV_A), "ijkl,klab->ijab", o4, t2, error)
      call plan_one(plan(P_OOVV_B), "ikab,jk->ijab", t2, oo, error)
      call plan_one(plan(P_OOVV_C), "ijac,cb->ijab", t2, vv, error)
      call plan_one(plan(P_OOVV_D), "ijad,bd->ijab", t2, vv, error)
      call plan_one(plan(P_RING_A), "ldib,jlda->ijab", ovov, t2, error)
      call plan_one(plan(P_RING_B), "ldia,jldb->ijab", ovov, t2, error)
      call plan_one(plan(P_RING_C), "ldia,jlbd->ijab", ovov, t2, error)
      call plan_one(plan(P_M_B), "iajb,jbP->iaP", ovov, bov, error)
   end subroutine plan_blocks

   subroutine plan_three_index(no, nv, naux, plan, error)
      !! The plans of the factorised terms of dovvv
      integer, intent(in) :: no, nv, naux
      type(trc_contract_plan_t), intent(inout) :: plan(NPLAN)
      type(error_t), intent(inout) :: error

      integer :: t1(2), ovov(4), boo(3), bov(3), bvv(3), vv(2)

      t1 = [no, nv]
      ovov = [no, nv, no, nv]
      vv = [nv, nv]
      boo = [no, no, naux]
      bov = [no, nv, naux]
      bvv = [nv, nv, naux]
      call plan_one(plan(P_H_T), "kc,cbP->kbP", t1, bvv, error)
      call plan_one(plan(P_H_L), "jb,bcP->jcP", t1, bvv, error)
      call plan_one(plan(P_K_H), "iakb,kbP->iaP", ovov, bov, error)
      call plan_one(plan(P_K_B), "iakb,iaP->kbP", ovov, bov, error)
      call plan_one(plan(P_Z_T), "kb,kcP->bcP", t1, bov, error)
      call plan_one(plan(P_Z_L), "jb,jcP->bcP", t1, bov, error)
      call plan_one(plan(P_R), "kbic,bcP->kiP", ovov, bvv, error)
      call plan_one(plan(P_R_T), "kiP,ka->iaP", boo, t1, error)
      call plan_one(plan(P_MOO), "ka,iaP->kiP", t1, bov, error)
      call plan_one(plan(P_P2_MOO), "kbic,kiP->bcP", ovov, boo, error)
      call plan_one(plan(P_N), "ba,bcP->acP", vv, bvv, error)
      call plan_one(plan(P_N_T), "ic,acP->iaP", t1, bvv, error)
      call plan_one(plan(P_BT), "iaP,ic->acP", bov, t1, error)
      call plan_one(plan(P_M_BT), "ba,acP->bcP", vv, bvv, error)
   end subroutine plan_three_index

   subroutine allocate_all(options, wk, den, no, nv, naux)
      !! Chooses the slab, checks that result and temporaries fit, then allocates and maps them (module header,
      !! MEMORY). Not fitting sets `den%status` and allocates nothing.
      type(trc_rccsd_density_options_t), intent(in) :: options
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      integer, intent(in) :: no, nv, naux

      integer(int64) :: work, held, per_slab, free, budget
      integer :: nb, cap

      work = largest_work(wk%plan)
      held = result_size(no, nv, naux) + fixed_size(no, nv, naux, work)
      per_slab = slab_size(nv, naux)
      cap = max(1, nv/SLAB_FRACTION)
      if (options%slab_batch > 0) then
         nb = min(options%slab_batch, cap)
      else
         budget = SLAB_BUDGET
         free = device_memory_free_bytes()
         if (free >= 0) budget = min(budget, max(0_int64, free - held))
         nb = int(min(int(cap, int64), max(1_int64, budget/per_slab)))
      end if
      nb = forced_slab(nb)
      ! v^3 nb and v nb naux are indexed with default integers
      nb = int(min(int(nb, int64), int(huge(0_int32), int64)/max(int(nv, int64)**3, int(nv, int64)*naux)))
      nb = max(nb, 1)
      if (.not. device_memory_fits(den%status, "trc_rccsd_density: the densities and their work arrays (o "// &
                                   to_string(no)//", v "//to_string(nv)//", naux "//to_string(naux)//", slab "// &
                                   to_string(nb)//")", held + per_slab*nb)) return
      den%no = no
      den%nv = nv
      den%naux = naux
      den%slab_batch = nb
      den%work_doubles = work
      den%temp_bytes = fixed_size(no, nv, naux, work) + per_slab*nb
      allocate (den%gamma_oo(no, no), den%gamma_ov(no, nv), den%gamma_vv(nv, nv))
      allocate (den%g3_oo(no*no, naux), den%g3_ov(no*nv, naux), den%g3_vv(nv*nv, naux))
      !$acc enter data create(den%gamma_oo, den%gamma_ov, den%gamma_vv, den%g3_oo, den%g3_ov, den%g3_vv)
      call allocate_work(wk, work, no, nv, naux, nb)
   end subroutine allocate_all

   subroutine allocate_work(wk, work, no, nv, naux, nb)
      !! The temporaries of `work_t`, mapped
      type(work_t), intent(inout) :: wk
      integer(int64), intent(in) :: work
      integer, intent(in) :: no, nv, naux, nb

      integer :: n2

      n2 = no*no*nv*nv
      allocate (wk%work(work), wk%tau(n2), wk%pa(n2), wk%pb(n2), wk%m(n2), wk%w1(n2), wk%w2(n2), wk%w3(n2))
      allocate (wk%s1(no*no*no*nv), wk%s2(no*no*no*nv), wk%goooo(no**4), wk%d4(no**4), wk%k4(no**4))
      allocate (wk%hov(no*nv*naux), wk%hoo(no*no*naux), wk%hvv(nv*nv*naux), wk%u(naux))
      allocate (wk%cbuf(nv*nv*nv*nb), wk%ebuf(nv*nv*nv*nb), wk%ybuf(nv*nb*naux))
      allocate (wk%doo(no, no), wk%dvv(nv, nv), wk%xtv(no, no), wk%moo(no, no), wk%mij(no, no))
      allocate (wk%mvv(nv, nv), wk%mab(nv, nv), wk%m2(nv, nv), wk%mia(no, nv), wk%d1oo(no, no), wk%d1vv(nv, nv))
      !$acc enter data create(wk%work, wk%tau, wk%pa, wk%pb, wk%m, wk%w1, wk%w2, wk%w3)
      !$acc enter data create(wk%s1, wk%s2, wk%goooo, wk%d4, wk%k4, wk%hov, wk%hoo, wk%hvv, wk%u)
      !$acc enter data create(wk%cbuf, wk%ebuf, wk%ybuf)
      !$acc enter data create(wk%doo, wk%dvv, wk%xtv, wk%moo, wk%mij, wk%mvv, wk%mab, wk%m2, wk%mia)
      !$acc enter data create(wk%d1oo, wk%d1vv)
   end subroutine allocate_work

   subroutine free_work(wk)
      !! Unmaps and frees what `allocate_work` made
      type(work_t), intent(inout) :: wk

      if (.not. allocated(wk%work)) return
      !$acc exit data delete(wk%work, wk%tau, wk%pa, wk%pb, wk%m, wk%w1, wk%w2, wk%w3)
      !$acc exit data delete(wk%s1, wk%s2, wk%goooo, wk%d4, wk%k4, wk%hov, wk%hoo, wk%hvv, wk%u)
      !$acc exit data delete(wk%cbuf, wk%ebuf, wk%ybuf)
      !$acc exit data delete(wk%doo, wk%dvv, wk%xtv, wk%moo, wk%mij, wk%mvv, wk%mab, wk%m2, wk%mia)
      !$acc exit data delete(wk%d1oo, wk%d1vv)
      deallocate (wk%work, wk%tau, wk%pa, wk%pb, wk%m, wk%w1, wk%w2, wk%w3)
      deallocate (wk%s1, wk%s2, wk%goooo, wk%d4, wk%k4, wk%hov, wk%hoo, wk%hvv, wk%u)
      deallocate (wk%cbuf, wk%ebuf, wk%ybuf)
      deallocate (wk%doo, wk%dvv, wk%xtv, wk%moo, wk%mij, wk%mvv, wk%mab, wk%m2, wk%mia, wk%d1oo, wk%d1vv)
   end subroutine free_work

   function forced_slab(nb) result(capped)
      !! `nb` capped by TERCO_CCDENS_D_CHUNK when that is a positive integer
      integer, intent(in) :: nb
      integer :: capped

      character(len=16) :: text
      integer :: stat, cap

      capped = nb
      call get_environment_variable("TERCO_CCDENS_D_CHUNK", text, status=stat)
      if (stat /= 0) return
      read (text, *, iostat=stat) cap
      if (stat == 0 .and. cap >= 1) capped = min(nb, cap)
   end function forced_slab

   subroutine run_steps(ints, amps, l, wk, den, la)
      !! The six steps, each a row of `den%steps`; the order matters (each leaves arrays for the next)
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      type(timer_type) :: tc
      integer :: k

      call kernel_zero(den%no*den%no*den%naux, den%g3_oo)
      call kernel_zero(den%no*den%nv*den%naux, den%g3_ov)
      call kernel_zero(den%nv*den%nv*den%naux, den%g3_vv)
      do k = 1, RCCSD_DENSITY_NSTEP
         wk%flops = 0
         call rccsd_steps_open(den%steps, la, k, tc)
         select case (k)
         case (RCCSD_DENSITY_ST_ONE)
            call step_one_particle(amps, l, wk, den, la)
         case (RCCSD_DENSITY_ST_OOOO)
            call step_oooo(ints, amps, l, wk, den, la)
         case (RCCSD_DENSITY_ST_OVVO)
            call step_ovvo(ints, amps, l, wk, den, la)
         case (RCCSD_DENSITY_ST_OVOV)
            call step_ovov(ints, amps, l, wk, den, la)
         case (RCCSD_DENSITY_ST_OVVV)
            call step_ovvv(ints, amps, l, wk, den, la)
         case default
            call step_vvvv(ints, amps, l, wk, den, la)
         end select
         den%flops(k) = wk%flops
         call rccsd_steps_close(den%steps, la, k, tc, wk%flops)
      end do
   end subroutine run_steps

   subroutine run(wk, k, la, a, b, c, alpha, beta)
      !! C = alpha contraction_k(A, B) + beta C; a refusal is kept in `wk%failure`
      type(work_t), intent(inout) :: wk
      integer, intent(in) :: k
      type(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: a(*), b(*)
      real(dp), intent(inout) :: c(*)
      real(dp), intent(in) :: alpha, beta

      call wk%plan(k)%apply(la, a, b, c, wk%work, alpha, beta)
      if (wk%plan(k)%status%has_error() .and. .not. wk%failure%has_error()) then
         call wk%failure%set(ERROR_VALIDATION, wk%plan(k)%status%get_message())
      end if
      wk%flops = wk%flops + rccsd_steps_gemm_flops(wk%plan(k)%gemm%m, wk%plan(k)%gemm%n, wk%plan(k)%gemm%k, &
                                                   wk%plan(k)%gemm%count)
   end subroutine run

   subroutine permute(wk, spec, extent, source, destination, alpha, beta)
      !! `contract_permute`, a refusal kept in `wk%failure`
      type(work_t), intent(inout) :: wk
      character(len=*), intent(in) :: spec
      integer, intent(in) :: extent(:)
      real(dp), intent(in) :: source(*)
      real(dp), intent(inout) :: destination(*)
      real(dp), intent(in) :: alpha, beta

      type(error_t) :: err

      call contract_permute(spec, extent, source, destination, err, alpha, beta)
      if (err%has_error() .and. .not. wk%failure%has_error()) call wk%failure%set(ERROR_VALIDATION, err%get_message())
   end subroutine permute

   subroutine step_one_particle(amps, l, wk, den, la)
      !! gamma (pyscf _gamma1_intermediates), then tau, pa and pb and the small intermediates mia, mab, mij, moo,
      !! mvv of the two-particle blocks
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv

      no = den%no
      nv = den%nv
      ! w1 = theta(t2) = 2 t2 - t2 swapped
      call permute(wk, "ijab->ijab", [no, no, nv, nv], amps%t2, wk%w1, 2.0_dp, 0.0_dp)
      call permute(wk, "ijab->ijba", [no, no, nv, nv], amps%t2, wk%w1, -1.0_dp, 1.0_dp)
      ! doo(ij) = sum_kab l2(ikab) theta(jkab),  dvv(ab) = sum_ijc theta(jica) l2(jicb)
      call run(wk, P_DOO, la, l%t2, wk%w1, wk%doo, 1.0_dp, 0.0_dp)
      call run(wk, P_DVV, la, wk%w1, l%t2, wk%dvv, 1.0_dp, 0.0_dp)
      ! pyscf's doo and dvv, then gamma = d + d^T
      call kernel_axpby(no*no, -1.0_dp, wk%doo, 0.0_dp, wk%d1oo)
      call run(wk, P_GOO, la, l%t1, amps%t1, wk%d1oo, -1.0_dp, 1.0_dp)
      call kernel_add_transpose(no, wk%d1oo, den%gamma_oo)
      call kernel_axpby(nv*nv, 1.0_dp, wk%dvv, 0.0_dp, wk%d1vv)
      call run(wk, P_GVV, la, amps%t1, l%t1, wk%d1vv, 1.0_dp, 1.0_dp)
      call kernel_add_transpose(nv, wk%d1vv, den%gamma_vv)
      ! gamma(ov) = l1 + dvo^T
      call kernel_axpby(no*nv, 1.0_dp, l%t1, 0.0_dp, den%gamma_ov)
      call kernel_axpby(no*nv, 1.0_dp, amps%t1, 1.0_dp, den%gamma_ov)
      call run(wk, P_XTV, la, amps%t1, l%t1, wk%xtv, 1.0_dp, 0.0_dp)
      call run(wk, P_GOV_A, la, wk%xtv, amps%t1, den%gamma_ov, -1.0_dp, 1.0_dp)
      call run(wk, P_GOV_B, la, wk%w1, l%t1, den%gamma_ov, 1.0_dp, 1.0_dp)
      call run(wk, P_GOV_C, la, wk%doo, amps%t1, den%gamma_ov, -1.0_dp, 1.0_dp)
      call run(wk, P_GOV_D, la, amps%t1, wk%dvv, den%gamma_ov, -1.0_dp, 1.0_dp)
      ! mia(ia) = sum_kc theta(ikac) l1(kc)
      call run(wk, P_MIA, la, wk%w1, l%t1, wk%mia, 1.0_dp, 0.0_dp)

      call kernel_tau(no, nv, amps%t1, amps%t2, wk%tau)
      ! pa(iajb) = sum l2(ikca) t2(jkcb); pb(iajb) = -sum l2(ikca) t2(jkbc) + sum l2(ikac) theta(jkbc)
      call run(wk, P_PA, la, l%t2, amps%t2, wk%pa, 1.0_dp, 0.0_dp)
      call run(wk, P_PB_A, la, l%t2, amps%t2, wk%pb, -1.0_dp, 0.0_dp)
      call run(wk, P_PB_B, la, l%t2, wk%w1, wk%pb, 1.0_dp, 1.0_dp)
      ! moo(jl) = sum_d 2 pa(ldjd) + pb(ldjd);  mvv(db) = sum_l 2 pa(lbld) + pb(lbld)
      call kernel_trace_oo(no, nv, wk%pa, wk%pb, wk%moo)
      call kernel_trace_vv(no, nv, wk%pa, wk%pb, wk%mvv)
      call run(wk, P_MAB, la, l%t1, amps%t1, wk%mab, 1.0_dp, 0.0_dp)
      call kernel_axpby(no*no, 0.5_dp, wk%moo, 0.0_dp, wk%mij)
      call run(wk, P_MIJ, la, amps%t1, l%t1, wk%mij, 1.0_dp, 1.0_dp)
   end subroutine step_one_particle

   subroutine step_oooo(ints, amps, l, wk, den, la)
      !! goooo, the dooov and doooo blocks and their contractions with b
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv
      integer :: o4(4), ooov(4)

      no = den%no
      nv = den%nv
      o4 = [no, no, no, no]
      ooov = [no, no, no, nv]
      call run(wk, P_GOOOO, la, wk%tau, l%t2, wk%goooo, 0.5_dp, 0.0_dp)
      ! gooov(jkia) in s1
      call run(wk, P_G3_A, la, amps%t1, wk%pa, wk%s1, 1.0_dp, 0.0_dp)
      call run(wk, P_G3_B, la, amps%t1, wk%pb, wk%s1, -1.0_dp, 1.0_dp)
      call kernel_ooov_outer(no, nv, wk%moo, amps%t1, wk%s1)
      call run(wk, P_G3_C, la, wk%goooo, amps%t1, wk%s1, 2.0_dp, 1.0_dp)
      call run(wk, P_G3_D, la, l%t1, wk%tau, wk%s1, -1.0_dp, 1.0_dp)
      call run(wk, P_G3_E, la, l%t2, amps%t1, wk%s1, -1.0_dp, 1.0_dp)
      ! dooov(pqra) = 2 gooov(prqa) - gooov(rpqa), in s2
      call permute(wk, "jkia->jika", ooov, wk%s1, wk%s2, 2.0_dp, 0.0_dp)
      call permute(wk, "jkia->kija", ooov, wk%s1, wk%s2, -1.0_dp, 1.0_dp)
      call run(wk, P_OOOV_OO, la, wk%s2, ints%b_ov, den%g3_oo, 2.0_dp, 1.0_dp)
      call run(wk, P_OOOV_OV, la, wk%s2, ints%b_oo, den%g3_ov, 1.0_dp, 1.0_dp)
      ! doooo(pqrs) = 2 goooo(prqs) - goooo(prsq), then K = d(pqrs) + d(qpsr) + d(rspq) + d(srqp)
      call permute(wk, "ijkl->ikjl", o4, wk%goooo, wk%d4, 2.0_dp, 0.0_dp)
      call permute(wk, "ijkl->iljk", o4, wk%goooo, wk%d4, -1.0_dp, 1.0_dp)
      call permute(wk, "pqrs->pqrs", o4, wk%d4, wk%k4, 1.0_dp, 0.0_dp)
      call permute(wk, "pqrs->qpsr", o4, wk%d4, wk%k4, 1.0_dp, 1.0_dp)
      call permute(wk, "pqrs->rspq", o4, wk%d4, wk%k4, 1.0_dp, 1.0_dp)
      call permute(wk, "pqrs->srqp", o4, wk%d4, wk%k4, 1.0_dp, 1.0_dp)
      call run(wk, P_K_OO, la, wk%k4, ints%b_oo, den%g3_oo, 1.0_dp, 1.0_dp)
   end subroutine step_oooo

   subroutine step_ovvo(ints, amps, l, wk, den, la)
      !! pyscf's gOvvO and govVO; doovv and its contractions with b; dovvo's share of `wk%m`
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, nov
      integer :: ovov(4)

      no = den%no
      nv = den%nv
      nov = no*nv
      ovov = [no, nv, no, nv]
      ! w1(iajb) = gOvvO(iabj) = pa + sum_k [sum_c l2(ikca) t1(jc)] t1(kb)
      call kernel_axpby(nov*nov, 1.0_dp, wk%pa, 0.0_dp, wk%w1)
      call run(wk, P_UP, la, l%t2, amps%t1, wk%s1, 1.0_dp, 0.0_dp)
      call run(wk, P_U_T, la, wk%s1, amps%t1, wk%w1, 1.0_dp, 1.0_dp)
      ! w2(iajb) = govVO(iabj) = pb + l1(ia) t1(jb) - sum_k [sum_c l2(ikac) t1(jc)] t1(kb)
      call kernel_axpby(nov*nov, 1.0_dp, wk%pb, 0.0_dp, wk%w2)
      call la%gemm("N", "N", nov, nov, 1, 1.0_dp, l%t1, nov, amps%t1, 1, 1.0_dp, wk%w2, nov)
      wk%flops = wk%flops + rccsd_steps_gemm_flops(nov, nov, 1)
      call run(wk, P_U, la, l%t2, amps%t1, wk%s2, 1.0_dp, 0.0_dp)
      call run(wk, P_U_T, la, wk%s2, amps%t1, wk%w2, -1.0_dp, 1.0_dp)
      ! m(iajb) = dovvo(iabj) so far = 2 govVO + gOvvO;  w3 = Y = 2 gOvvO + govVO
      call kernel_axpby(nov*nov, 1.0_dp, wk%w1, 0.0_dp, wk%m)
      call kernel_axpby(nov*nov, 2.0_dp, wk%w2, 1.0_dp, wk%m)
      call kernel_axpby(nov*nov, 2.0_dp, wk%w1, 0.0_dp, wk%w3)
      call kernel_axpby(nov*nov, 1.0_dp, wk%w2, 1.0_dp, wk%w3)
      ! w1(ijab) = doovv(ijab) + doovv(jiba) = - Y(jaib) - Y(ibja)
      call permute(wk, "jaib->ijab", ovov, wk%w3, wk%w1, -1.0_dp, 0.0_dp)
      call permute(wk, "ibja->ijab", ovov, wk%w3, wk%w1, -1.0_dp, 1.0_dp)
      call run(wk, P_OOVV_OO, la, wk%w1, ints%b_vv, den%g3_oo, 1.0_dp, 1.0_dp)
      call run(wk, P_OOVV_VV, la, wk%w1, ints%b_oo, den%g3_vv, 1.0_dp, 1.0_dp)
   end subroutine step_ovvo

   subroutine step_ovov(ints, amps, l, wk, den, la)
      !! pyscf's goovv in w1; adds dovov to `wk%m`, symmetrises it and contracts it with b_ov into g3_ov
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, nov
      integer :: oovv(4)

      no = den%no
      nv = den%nv
      nov = no*nv
      oovv = [no, no, nv, nv]
      call kernel_oovv_start(no, nv, wk%mia, amps%t1, wk%tau, l%t2, wk%w1)
      call run(wk, P_OOVV_A, la, wk%goooo, wk%tau, wk%w1, 1.0_dp, 1.0_dp)
      call run(wk, P_OOVV_B, la, wk%tau, wk%mij, wk%w1, -1.0_dp, 1.0_dp)
      call run(wk, P_OOVV_C, la, amps%t2, wk%mab, wk%w1, -1.0_dp, 1.0_dp)
      call run(wk, P_OOVV_D, la, wk%tau, wk%mvv, wk%w1, -0.5_dp, 1.0_dp)
      ! w2 = tau - t2/2
      call kernel_axpby(nov*nov, 1.0_dp, wk%tau, 0.0_dp, wk%w2)
      call kernel_axpby(nov*nov, -0.5_dp, amps%t2, 1.0_dp, wk%w2)
      call run(wk, P_RING_A, la, wk%pa, wk%w2, wk%w1, 1.0_dp, 1.0_dp)
      call run(wk, P_RING_B, la, wk%pb, wk%w2, wk%w1, -1.0_dp, 1.0_dp)
      ! w3 = pb + pa/2
      call kernel_axpby(nov*nov, 1.0_dp, wk%pb, 0.0_dp, wk%w3)
      call kernel_axpby(nov*nov, 0.5_dp, wk%pa, 1.0_dp, wk%w3)
      call run(wk, P_RING_C, la, wk%w3, amps%t2, wk%w1, 1.0_dp, 1.0_dp)
      ! m(iajb) += dovov(iajb) = 2 goovv(ijab) - goovv(jiab), then m + m^T
      call permute(wk, "ijab->iajb", oovv, wk%w1, wk%m, 2.0_dp, 1.0_dp)
      call permute(wk, "jiab->iajb", oovv, wk%w1, wk%m, -1.0_dp, 1.0_dp)
      call kernel_symmetrise_sum(nov, wk%m)
      call run(wk, P_M_B, la, wk%m, ints%b_ov, den%g3_ov, 1.0_dp, 1.0_dp)
   end subroutine step_ovov

   subroutine step_ovvv(ints, amps, l, wk, den, la)
      !! The terms of dovvv that factor through one occupied index, contracted with b_vv (into g3_ov) and with
      !! b_ov (into g3_vv). Leaves wk%hvv free.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      call ovvv_pair_terms(ints, amps, l, wk, den, la)
      call ovvv_exchange_terms(ints, amps, wk, den, la)
      call ovvv_single_terms(ints, amps, l, wk, den, la)
   end subroutine step_ovvv

   subroutine ovvv_pair_terms(ints, amps, l, wk, den, la)
      !! The terms l1 t2, t1 l2 and the half of the pa/pb terms that contracts t1 with b_vv:
      !!     K4(ia,jc) = 2 t2(ijac) - t2(ijca)
      !!     K7(ia,kb) = 2 l2(ikab) - l2(ikba) + 2 pb(kbia) + pa(kbia)
      !!     g3_ov += K7 H7 + K4 H4,   H7(kb;P) = sum_c t1(kc) b_vv(cb,P),  H4 with l1
      !!     g3_vv += 2 [t1^T (K7^T b_ov) + l1^T (K4^T b_ov)]
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv
      integer :: oovv(4), ovov(4)

      no = den%no
      nv = den%nv
      oovv = [no, no, nv, nv]
      ovov = [no, nv, no, nv]
      call permute(wk, "ijac->iajc", oovv, l%t2, wk%w1, 2.0_dp, 0.0_dp)
      call permute(wk, "ijca->iajc", oovv, l%t2, wk%w1, -1.0_dp, 1.0_dp)
      call permute(wk, "kbia->iakb", ovov, wk%pb, wk%w1, 2.0_dp, 1.0_dp)
      call permute(wk, "kbia->iakb", ovov, wk%pa, wk%w1, 1.0_dp, 1.0_dp)
      call run(wk, P_H_T, la, amps%t1, ints%b_vv, wk%hov, 1.0_dp, 0.0_dp)
      call run(wk, P_K_H, la, wk%w1, wk%hov, den%g3_ov, 1.0_dp, 1.0_dp)
      call run(wk, P_K_B, la, wk%w1, ints%b_ov, wk%hov, 1.0_dp, 0.0_dp)
      call run(wk, P_Z_T, la, amps%t1, wk%hov, den%g3_vv, 2.0_dp, 1.0_dp)

      call permute(wk, "ijac->iajc", oovv, amps%t2, wk%w1, 2.0_dp, 0.0_dp)
      call permute(wk, "ijca->iajc", oovv, amps%t2, wk%w1, -1.0_dp, 1.0_dp)
      call run(wk, P_H_L, la, l%t1, ints%b_vv, wk%hov, 1.0_dp, 0.0_dp)
      call run(wk, P_K_H, la, wk%w1, wk%hov, den%g3_ov, 1.0_dp, 1.0_dp)
      call run(wk, P_K_B, la, wk%w1, ints%b_ov, wk%hov, 1.0_dp, 0.0_dp)
      call run(wk, P_Z_L, la, l%t1, wk%hov, den%g3_vv, 2.0_dp, 1.0_dp)
   end subroutine ovvv_pair_terms

   subroutine ovvv_exchange_terms(ints, amps, wk, den, la)
      !! The other half of the pa/pb terms: with P2(kbic) = pb(kbic) + 2 pa(kbic),
      !!     g3_ov(ia) -= sum_k R(ki) t1(ka),     R(ki;P) = sum_bc P2(kbic) b_vv(bc,P)
      !!     g3_vv(bc) -= 2 sum_ki P2(kbic) M(ki),  M(ki;P) = sum_a t1(ka) b_ov(ia,P)
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      integer :: n

      n = den%no*den%nv*den%no*den%nv
      call kernel_axpby(n, 1.0_dp, wk%pb, 0.0_dp, wk%w1)
      call kernel_axpby(n, 2.0_dp, wk%pa, 1.0_dp, wk%w1)
      call run(wk, P_R, la, wk%w1, ints%b_vv, wk%hoo, 1.0_dp, 0.0_dp)
      call run(wk, P_R_T, la, wk%hoo, amps%t1, den%g3_ov, -1.0_dp, 1.0_dp)
      call run(wk, P_MOO, la, amps%t1, ints%b_ov, wk%hoo, 1.0_dp, 0.0_dp)
      call run(wk, P_P2_MOO, la, wk%w1, wk%hoo, den%g3_vv, -2.0_dp, 1.0_dp)
   end subroutine ovvv_exchange_terms

   subroutine ovvv_single_terms(ints, amps, l, wk, den, la)
      !! The terms l1 t1 t1 and mvv t1: with M(xy) = sum_j l1(jx) t1(jy) + mvv(yx)/2,
      !!     g3_ov(ia) += 2 t1(ia) u - sum_c t1(ic) N(ac),   u = sum_bc M(bc) b_vv(bc,P), N(ac) = sum_b M(ba) b_vv(bc,P)
      !!     g3_vv(bc) += 4 M(bc) w - 2 sum_a M(ba) Bt(ac),   w = sum_ia t1(ia) b_ov(ia,P)
      !! Leaves Bt(a,d;P) = sum_i b_ov(ia,P) t1(id) in wk%hvv for the slabs.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      integer :: no, nv, nov, nv2, naux

      no = den%no
      nv = den%nv
      naux = den%naux
      nov = no*nv
      nv2 = nv*nv
      call permute(wk, "yx->xy", [nv, nv], wk%mvv, wk%m2, 0.5_dp, 0.0_dp)
      call la%gemm("T", "N", nv, nv, no, 1.0_dp, l%t1, no, amps%t1, no, 1.0_dp, wk%m2, nv)
      call la%gemv("T", nv2, naux, 1.0_dp, ints%b_vv, nv2, wk%m2, 0.0_dp, wk%u)
      call la%gemm("N", "N", nov, naux, 1, 2.0_dp, amps%t1, nov, wk%u, 1, 1.0_dp, den%g3_ov, nov)
      call run(wk, P_N, la, wk%m2, ints%b_vv, wk%hvv, 1.0_dp, 0.0_dp)
      call run(wk, P_N_T, la, amps%t1, wk%hvv, den%g3_ov, -1.0_dp, 1.0_dp)
      call la%gemv("T", nov, naux, 1.0_dp, ints%b_ov, nov, amps%t1, 0.0_dp, wk%u)
      call la%gemm("N", "N", nv2, naux, 1, 4.0_dp, wk%m2, nv2, wk%u, 1, 1.0_dp, den%g3_vv, nv2)
      call run(wk, P_BT, la, ints%b_ov, amps%t1, wk%hvv, 1.0_dp, 0.0_dp)
      call run(wk, P_M_BT, la, wk%m2, wk%hvv, den%g3_vv, -2.0_dp, 1.0_dp)
      wk%flops = wk%flops + rccsd_steps_gemm_flops(nv, nv, no) + rccsd_steps_gemm_flops(nov, naux, 1) &
                 + rccsd_steps_gemm_flops(nv2, naux, 1)
   end subroutine ovvv_single_terms

   subroutine step_vvvv(ints, amps, l, wk, den, la)
      !! The slabs over the virtual d (module header, BLOCKS), then the pair symmetrisation of g3_oo and g3_vv.
      !! Reads Bt in wk%hvv.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps, l
      type(work_t), intent(inout) :: wk
      type(trc_rccsd_density_t), intent(inout) :: den
      type(trc_linalg_t), intent(in) :: la

      type(trc_gemm_batch_t) :: per_p
      integer :: no, nv, naux, no2, nv2, d0, nd, nvd

      no = den%no
      nv = den%nv
      naux = den%naux
      no2 = no*no
      nv2 = nv*nv
      ! w1(ijsq) = tau(ijsq) - tau(ijqs)/2
      call permute(wk, "ijab->ijab", [no, no, nv, nv], wk%tau, wk%w1, 1.0_dp, 0.0_dp)
      call permute(wk, "ijab->ijba", [no, no, nv, nv], wk%tau, wk%w1, -0.5_dp, 1.0_dp)
      do d0 = 1, nv, den%slab_batch
         nd = min(den%slab_batch, nv - d0 + 1)
         nvd = nv*nd
         ! C'(sq,rd) = sum_ij w1(ij,sq) l2(ij,rd) for the slab's d
         call la%gemm("T", "N", nv2, nvd, no2, 1.0_dp, wk%w1, no2, l%t2(1, 1, 1, d0), no2, 0.0_dp, wk%cbuf, nv2)
         ! E(s,r,q,d) = C'(s,q,r,d)
         call permute(wk, "sqrd->srqd", [nv, nv, nv, nd], wk%cbuf, wk%ebuf, 1.0_dp, 0.0_dp)
         ! y(qd,P) = sum_sr E(sr,qd) b_vv(sr,P), which is y(d,q) of the module header
         call la%gemm("T", "N", nvd, naux, nv2, 1.0_dp, wk%ebuf, nv2, ints%b_vv, nv2, 0.0_dp, wk%ybuf, nvd)
         call kernel_add_rows(nv2, nvd, naux, nv*(d0 - 1), 4.0_dp, wk%ybuf, den%g3_vv)
         ! g3_vv(sr,P) -= 4 sum_ad E(sr,ad) Bt(ad,P): the rows of Bt for the slab's d
         call la%gemm("N", "N", nv2, naux, nvd, -4.0_dp, wk%ebuf, nv2, wk%hvv(nv*(d0 - 1) + 1), nv2, 1.0_dp, &
                      den%g3_vv, nv2)
         ! g3_ov(i,a,P) -= 2 sum_d t1(i,d) y(ad,P), one item per P, t1 shared
         per_p = trc_gemm_batch(no, nv, nd, naux, lda=no, ldb=nv, ldc=no, sa=0_int64, sb=int(nvd, int64), &
                                sc=int(no, int64)*nv, transb="T")
         call la%gemm_batched(per_p, amps%t1(1, d0), wk%ybuf, den%g3_ov, -2.0_dp, 1.0_dp)
         wk%flops = wk%flops + rccsd_steps_gemm_flops(nv2, nvd, no2) + 2_int64*rccsd_steps_gemm_flops(nvd, naux, nv2) &
                    + rccsd_steps_gemm_flops(no, nv, nd, naux)
      end do
      call kernel_symmetrise_pairs(no, naux, den%g3_oo)
      call kernel_symmetrise_pairs(nv, naux, den%g3_vv)
   end subroutine step_vvvv

   subroutine kernel_zero(n, x)
      !! ---------------------------------------------------------------------
      !! the kernels. Each takes its arrays as explicit-shape dummies and runs one
      !! `do concurrent` over ONE flattened index. The pure reorderings are
      !! trc_contract's `contract_permute`.
      !! ---------------------------------------------------------------------
      !! x = 0
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n)
      integer :: t

      do concurrent(t=1:n)
         x(t) = 0.0_dp
      end do
   end subroutine kernel_zero

   subroutine kernel_axpby(n, a, x, b, y)
      !! y = a x + b y; y is not read when b = 0
      integer, intent(in) :: n
      real(dp), intent(in) :: a, b
      real(dp), intent(in) :: x(n)
      real(dp), intent(inout) :: y(n)
      integer :: t

      if (abs(b) <= 0.0_dp) then
         do concurrent(t=1:n)
            y(t) = a*x(t)
         end do
      else
         do concurrent(t=1:n)
            y(t) = a*x(t) + b*y(t)
         end do
      end if
   end subroutine kernel_axpby

   subroutine kernel_add_transpose(n, d, g)
      !! g = d + d^T
      integer, intent(in) :: n
      real(dp), intent(in) :: d(n, n)
      real(dp), intent(out) :: g(n*n)
      integer :: t

      do concurrent(t=1:n*n)
         g(t) = d(mod(t - 1, n) + 1, (t - 1)/n + 1) + d((t - 1)/n + 1, mod(t - 1, n) + 1)
      end do
   end subroutine kernel_add_transpose

   subroutine kernel_symmetrise_sum(n, a)
      !! a = a + a^T in place: the iteration of the upper triangle writes both elements of its pair
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(n, n)
      integer :: t, i, j
      real(dp) :: s

      do concurrent(t=0:n*n - 1) local(i, j, s)
         i = mod(t, n) + 1
         j = t/n + 1
         if (i <= j) then
            s = a(i, j) + a(j, i)
            a(i, j) = s
            a(j, i) = s
         end if
      end do
   end subroutine kernel_symmetrise_sum

   subroutine kernel_symmetrise_pairs(n, m, a)
      !! a(ij,P) = (a(ij,P) + a(ji,P))/2 in place, for each of the m columns P
      integer, intent(in) :: n, m
      real(dp), intent(inout) :: a(n, n, m)
      integer :: t, i, j, p
      real(dp) :: s

      do concurrent(t=0:n*n*m - 1) local(i, j, p, s)
         i = mod(t, n) + 1
         j = mod(t/n, n) + 1
         p = t/(n*n) + 1
         if (i < j) then
            s = 0.5_dp*(a(i, j, p) + a(j, i, p))
            a(i, j, p) = s
            a(j, i, p) = s
         end if
      end do
   end subroutine kernel_symmetrise_pairs

   subroutine kernel_add_rows(nrows, n, m, offset, alpha, x, y)
      !! y(offset + r, P) += alpha x(r, P) for r = 1..n and the m columns P
      integer, intent(in) :: nrows, n, m, offset
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(n, m)
      real(dp), intent(inout) :: y(nrows, m)
      integer :: t, r, p

      do concurrent(t=0:n*m - 1) local(r, p)
         r = mod(t, n) + 1
         p = t/n + 1
         y(offset + r, p) = y(offset + r, p) + alpha*x(r, p)
      end do
   end subroutine kernel_add_rows

   subroutine kernel_tau(no, nv, t1, t2, tau)
      !! tau(ijab) = t2(ijab) + t1(ia) t1(jb)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv)
      real(dp), intent(out) :: tau(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         tau(t) = tau_item(TRC_DC_INDEX(t), no, nv, t1, t2)
      end do
   end subroutine kernel_tau

   subroutine kernel_trace_oo(no, nv, pa, pb, moo)
      !! moo(j,l) = sum_d 2 pa(l,d,j,d) + pb(l,d,j,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: pa(no, nv, no, nv), pb(no, nv, no, nv)
      real(dp), intent(out) :: moo(no*no)
      integer :: t

      do concurrent(t=1:no*no)
         moo(t) = trace_oo_item(TRC_DC_INDEX(t), no, nv, pa, pb)
      end do
   end subroutine kernel_trace_oo

   subroutine kernel_trace_vv(no, nv, pa, pb, mvv)
      !! mvv(d,b) = sum_l 2 pa(l,b,l,d) + pb(l,b,l,d)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: pa(no, nv, no, nv), pb(no, nv, no, nv)
      real(dp), intent(out) :: mvv(nv*nv)
      integer :: t

      do concurrent(t=1:nv*nv)
         mvv(t) = trace_vv_item(TRC_DC_INDEX(t), no, nv, pa, pb)
      end do
   end subroutine kernel_trace_vv

   subroutine kernel_ooov_outer(no, nv, moo, t1, g)
      !! g(j,k,i,a) -= moo(j,i) t1(k,a)/2
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: moo(no, no), t1(no, nv)
      real(dp), intent(inout) :: g(no*no*no*nv)
      integer :: t

      do concurrent(t=1:no*no*no*nv)
         g(t) = g(t) + ooov_outer_item(TRC_DC_INDEX(t), no, nv, moo, t1)
      end do
   end subroutine kernel_ooov_outer

   subroutine kernel_oovv_start(no, nv, mia, t1, tau, l2, g)
      !! g(ijab) = mia(ia) t1(jb) + (tau(ijab) + l2(ijab))/2
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: mia(no, nv), t1(no, nv), tau(no*no*nv*nv), l2(no*no*nv*nv)
      real(dp), intent(out) :: g(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         g(t) = outer_item(TRC_DC_INDEX(t), no, nv, mia, t1) + 0.5_dp*(tau(t) + l2(t))
      end do
   end subroutine kernel_oovv_start

   pure function tau_item(t, no, nv, t1, t2) result(w)
      !! tau at position t of (no, no, nv, nv)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv)
      real(dp) :: w
      integer :: q, i, j, a, b

      q = t - 1
      i = mod(q, no) + 1
      q = q/no
      j = mod(q, no) + 1
      q = q/no
      a = mod(q, nv) + 1
      b = q/nv + 1
      w = t2(i, j, a, b) + t1(i, a)*t1(j, b)
   end function tau_item

   pure function trace_oo_item(t, no, nv, pa, pb) result(s)
      !! The sum of `kernel_trace_oo` for t the position of (j,l)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: pa(no, nv, no, nv), pb(no, nv, no, nv)
      real(dp) :: s
      integer :: j, l, d

      j = mod(t - 1, no) + 1
      l = (t - 1)/no + 1
      s = 0.0_dp
      do d = 1, nv
         s = s + 2.0_dp*pa(l, d, j, d) + pb(l, d, j, d)
      end do
   end function trace_oo_item

   pure function trace_vv_item(t, no, nv, pa, pb) result(s)
      !! The sum of `kernel_trace_vv` for t the position of (d,b)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: pa(no, nv, no, nv), pb(no, nv, no, nv)
      real(dp) :: s
      integer :: d, b, l

      d = mod(t - 1, nv) + 1
      b = (t - 1)/nv + 1
      s = 0.0_dp
      do l = 1, no
         s = s + 2.0_dp*pa(l, b, l, d) + pb(l, b, l, d)
      end do
   end function trace_vv_item

   pure function ooov_outer_item(t, no, nv, moo, t1) result(w)
      !! -moo(j,i) t1(k,a)/2 at position t of (j,k,i,a)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: moo(no, no), t1(no, nv)
      real(dp) :: w
      integer :: q, j, k, i, a

      q = t - 1
      j = mod(q, no) + 1
      q = q/no
      k = mod(q, no) + 1
      q = q/no
      i = mod(q, no) + 1
      a = q/no + 1
      w = -0.5_dp*moo(j, i)*t1(k, a)
   end function ooov_outer_item

   pure function outer_item(t, no, nv, x, y) result(w)
      !! x(i,a) y(j,b) at position t of (i,j,a,b)
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: x(no, nv), y(no, nv)
      real(dp) :: w
      integer :: q, i, j, a, b

      q = t - 1
      i = mod(q, no) + 1
      q = q/no
      j = mod(q, no) + 1
      q = q/no
      a = mod(q, nv) + 1
      b = q/nv + 1
      w = x(i, a)*y(j, b)
   end function outer_item

end module trc_rccsd_density
