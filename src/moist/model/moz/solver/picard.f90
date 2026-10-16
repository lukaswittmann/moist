!> Picard iteration with damped fixed-point updates
module moist_model_moz_solver_picard
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_moz_solver_base, only: &
      & moist_model_moz_solver_type, validate_solver_arrays
   use moist_model_moz_solver_evaluator, only: moist_model_moz_solver_evaluator_type
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_picard_type, new_picard_solver

   !> Damped fixed-point mixing
   type, extends(moist_model_moz_solver_type) :: moist_model_moz_solver_picard_type
      !> Fraction of the fixed-point update to apply
      real(wp) :: damping = 0.1_wp
   contains
      procedure :: outer_step => picard_outer_step
      procedure :: name => picard_name
   end type moist_model_moz_solver_picard_type

contains

   !> Display name of the solver
   !>
   !> @param[in] self Solver to name
   pure function picard_name(self) result(text)
      !> Solver to name
      class(moist_model_moz_solver_picard_type), intent(in) :: self
      !> Display name
      character(len=:), allocatable :: text
      text = "Picard"
   end function picard_name

   !> Construct a Picard solver with optional damping and verbosity
   !>
   !> @param[in] damping Fraction of the fixed-point update
   !> @param[in] verbosity Diagnostic verbosity
   function new_picard_solver(damping, verbosity) result(self)
      !> Fraction of the fixed-point update
      real(wp), intent(in), optional :: damping
      !> Diagnostic verbosity
      integer, intent(in), optional :: verbosity
      !> Constructed solver
      type(moist_model_moz_solver_picard_type) :: self

      if (present(damping)) self%damping = damping
      if (present(verbosity)) self%verbosity = verbosity
   end function new_picard_solver

   !> One Picard outer iteration: evaluate residual, apply damped linear mixing
   !>
   !> @param[inout] self Solver and step diagnostics
   !> @param[in] residual_fn Model residual evaluator
   !> @param[in] gamma Current field (npts, npair)
   !> @param[out] gamma_new Damped update (npts, npair)
   !> @param[out] rmse RMS residual of the input field
   !> @param[out] feasible Whether the residual is usable
   !> @param[out] error Evaluation failure
   subroutine picard_outer_step(self, residual_fn, gamma, gamma_new, rmse, feasible, error)
      !> Solver and step diagnostics
      class(moist_model_moz_solver_picard_type), intent(inout) :: self
      !> Model residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_fn
      !> Current field (npts, npair)
      real(wp), intent(in) :: gamma(:, :)
      !> Damped update (npts, npair)
      real(wp), intent(out) :: gamma_new(:, :)
      !> RMS residual of the input field
      real(wp), intent(out) :: rmse
      !> Whether the residual is usable
      logical, intent(out) :: feasible
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: F(:, :)

      feasible = .false.
      call validate_solver_arrays(gamma, gamma_new, error)
      if (allocated(error)) return
      allocate(F, mold=gamma)
      call residual_fn%residual(gamma, F, feasible, error)
      if (allocated(error)) feasible = .false.
      if (.not. feasible) return

      rmse = sqrt(sum(F**2) / real(size(F), wp))
      if (.not. ieee_is_finite(rmse)) then
         call fatal_error(error, "Solver residual norm must be finite")
         feasible = .false.
         return
      end if
      gamma_new = gamma - self%damping * F
      self%prev_step_maxabs = maxval(abs(gamma_new - gamma))
      self%prev_alpha = self%damping
   end subroutine picard_outer_step

end module moist_model_moz_solver_picard
