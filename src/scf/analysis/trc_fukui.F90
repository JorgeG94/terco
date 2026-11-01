!
! Hirshfeld-condensed Fukui indices, by finite difference in the electron
! count. Ported from metalquicha's mqc_czt_fukui, mqc_ -> trc_ renamed: the
! arithmetic below is finish_indices/check_sum_rule verbatim. mqc's own
! three-SCF orchestration (fukui_indices) has no counterpart here -- that
! belongs with the context box that runs the SCFs, trc_capi.F90's
! trc_fukui_hirshfeld, which calls trc_hirshfeld_populations at each of the
! three states and then this module to turn the three charge sets into
! indices.
!
module trc_fukui
   !! The Fukui function is the derivative of the density with respect to
   !! the number of electrons at fixed nuclei. Condensed onto atoms:
   !!
   !!     f_A+ = q_A(N) - q_A(N+1)      f_A- = q_A(N-1) - q_A(N)
   !!
   !! `f_A+` is where the molecule accepts charge (a nucleophile attacks
   !! there), `f_A-` where it gives charge up (an electrophile attacks
   !! there), and `f0 = (f+ + f-)/2` is the radical average. All three come
   !! from Hirshfeld charges of the SAME partition (`trc_hirshfeld_populations`,
   !! fixed by geometry and element alone) at three electron counts, so an
   !! index is the stockholder share of a density DIFFERENCE rather than a
   !! difference of two independent partitions.
   use pic_types, only: dp
   use trc_error, only: error_t, ERROR_VALIDATION
   use pic_strings, only: to_string
   implicit none
   private

   public :: trc_fukui_indices

   real(dp), parameter :: TOL_SUM_RULE = 1.0e-6_dp
      !! mqc_czt_fukui.f90's own TOL_SUM_RULE, verbatim: how far `sum(f+)`
      !! and `sum(f-)` may drift from one before the three states are
      !! refused rather than quoted. Looser than round-off because the
      !! charges are renormalised to Tr(D S) at each state independently
      !! (trc_hirshfeld's own comment) and a difference of two such
      !! renormalisations is not exact.

contains

   subroutine trc_fukui_indices(q_neutral, q_anion, q_cation, f_plus, f_minus, f_zero, error)
      !! The condensed indices from the three states' Hirshfeld charges,
      !! one atom order shared by all three (the caller's -- same geometry,
      !! same basis, same partition).
      real(dp), intent(in) :: q_neutral(:), q_anion(:), q_cation(:)
         !! Hirshfeld charges of the N, N+1 and N-1 electron states, same
         !! atom order and length.
      real(dp), allocatable, intent(out) :: f_plus(:), f_minus(:), f_zero(:)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      if (size(q_anion) /= size(q_neutral) .or. size(q_cation) /= size(q_neutral)) then
         call error%set(ERROR_VALIDATION, "trc_fukui_indices: the three states carry different "// &
                        "atom counts ("//trim(to_string(size(q_neutral)))//", "// &
                        trim(to_string(size(q_anion)))//", "//trim(to_string(size(q_cation)))// &
                        "); they must share one partition")
         return
      end if

      allocate (f_plus(size(q_neutral)), f_minus(size(q_neutral)), f_zero(size(q_neutral)))
      f_plus = q_neutral - q_anion
      f_minus = q_cation - q_neutral
      f_zero = 0.5_dp*(f_plus + f_minus)

      call check_sum_rule(f_plus, f_minus, error)
   end subroutine trc_fukui_indices

   subroutine check_sum_rule(f_plus, f_minus, error)
      !! Both indices must sum to one over the molecule: the three states
      !! differ by exactly one electron and the charges account for all of
      !! it. What this catches is a density that was not the total, an ion
      !! run at the wrong charge, or a partition that lost part of the
      !! molecule -- mqc_czt_fukui.f90's own check_sum_rule, verbatim.
      real(dp), intent(in) :: f_plus(:), f_minus(:)
      type(error_t), intent(inout) :: error

      real(dp) :: plus, minus

      if (error%has_error()) return
      plus = sum(f_plus)
      minus = sum(f_minus)
      if (abs(plus - 1.0_dp) > TOL_SUM_RULE .or. abs(minus - 1.0_dp) > TOL_SUM_RULE) then
         call error%set(ERROR_VALIDATION, "trc_fukui_indices: the condensed Fukui indices do "// &
                        "not sum to one: f+ gives "//trim(to_string(plus))//" and f- gives "// &
                        trim(to_string(minus))//". One electron was added and one removed, so "// &
                        "the difference has to account for exactly that much charge.")
      end if
   end subroutine check_sum_rule

end module trc_fukui
