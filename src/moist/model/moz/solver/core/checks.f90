!> Range predicates shared by the solver parameter checks
!>
!> - Each predicate is false for NaN and infinite values
module moist_model_moz_solver_core_checks
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: finite_positive, finite_nonnegative, unit_fraction, step_fraction

contains

   !> Whether a value is finite and positive
   !>
   !> @param[in] value  value to check
   elemental function finite_positive(value) result(valid)
      !> Value to check
      real(wp), intent(in) :: value
      !> True for a finite value above zero
      logical :: valid

      valid = ieee_is_finite(value) .and. value > 0.0_wp
   end function finite_positive

   !> Whether a value is finite and nonnegative
   !>
   !> @param[in] value  value to check
   elemental function finite_nonnegative(value) result(valid)
      !> Value to check
      real(wp), intent(in) :: value
      !> True for a finite value at or above zero
      logical :: valid

      valid = ieee_is_finite(value) .and. value >= 0.0_wp
   end function finite_nonnegative

   !> Whether a value is a fraction strictly between zero and one
   !>
   !> @param[in] value  fraction to check
   elemental function unit_fraction(value) result(valid)
      !> Fraction to check
      real(wp), intent(in) :: value
      !> True for a value in (0, 1)
      logical :: valid

      valid = ieee_is_finite(value) .and. value > 0.0_wp .and. value < 1.0_wp
   end function unit_fraction

   !> Whether a value is a step length in (0, 1]
   !>
   !> @param[in] value  step length to check
   elemental function step_fraction(value) result(valid)
      !> Step length to check
      real(wp), intent(in) :: value
      !> True for a value in (0, 1]
      logical :: valid

      valid = ieee_is_finite(value) .and. value > 0.0_wp .and. value <= 1.0_wp
   end function step_fraction

end module moist_model_moz_solver_core_checks
