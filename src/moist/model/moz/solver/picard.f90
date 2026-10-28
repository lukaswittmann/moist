!> Picard iteration with damped fixed-point updates
module moist_model_moz_solver_picard
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_model_moz_solver_core_base, only: &
      & moist_model_moz_solver_type, validate_solver_arrays
   use moist_model_moz_solver_core_checks, only: finite_positive
   use moist_model_moz_solver_core_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_parameters, only: moist_model_parameters_type
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_picard_type, moist_model_moz_solver_picard_parameters_type
   public :: new_picard_solver

   !> Picard controls
   !>
   !> JSON/TOML keys: damping
   type, extends(moist_model_parameters_type) :: moist_model_moz_solver_picard_parameters_type
      !> Fraction of the fixed-point update, also used while DIIS fills its history
      real(wp) :: damping = 0.1_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => picard_init_defaults
      !> Declare fields for JSON/TOML input, output and printing
      procedure :: register_entries => picard_register_entries
      !> Check the damping range
      procedure :: validate => picard_validate
   end type moist_model_moz_solver_picard_parameters_type

   !> Damped fixed-point mixing
   type, extends(moist_model_moz_solver_type) :: moist_model_moz_solver_picard_type
      !> Damping controls
      type(moist_model_moz_solver_picard_parameters_type) :: params
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

   !> Restore compiled defaults
   !>
   !> @param[inout] self Parameter values
   subroutine picard_init_defaults(self)
      !> Parameter values
      class(moist_model_moz_solver_picard_parameters_type), intent(inout) :: self
      type(moist_model_moz_solver_picard_parameters_type) :: defaults

      self%damping = defaults%damping
   end subroutine picard_init_defaults

   !> Declare parameter fields for JSON/TOML input, output and printing
   !>
   !> @param[inout] self Parameter values
   subroutine picard_register_entries(self)
      !> Parameter values
      class(moist_model_moz_solver_picard_parameters_type), intent(inout), target :: self

      call self%register_real_scalar("damping", self%damping, level=2)
   end subroutine picard_register_entries

   !> Check the damping range
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error Invalid damping
   subroutine picard_validate(self, error)
      !> Parameter values
      class(moist_model_moz_solver_picard_parameters_type), intent(inout) :: self
      !> Invalid damping
      type(error_type), allocatable, intent(out) :: error

      if (.not. finite_positive(self%damping)) then
         call fatal_error(error, "Picard damping must be finite and positive")
      end if
   end subroutine picard_validate

   !> Construct a Picard solver from its controls
   !>
   !> @param[in] ctx Run context, borrowed for the lifetime of the solver
   !> @param[in] parameters Picard controls; defaults when absent
   function new_picard_solver(ctx, parameters) result(self)
      !> Run context, borrowed for the lifetime of the solver
      type(moist_context_type), intent(in), target :: ctx
      !> Picard controls
      type(moist_model_moz_solver_picard_parameters_type), intent(in), optional :: parameters
      !> Constructed solver
      type(moist_model_moz_solver_picard_type) :: self

      self%ctx => ctx
      if (present(parameters)) self%params = parameters
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
      gamma_new = gamma - self%params%damping * F
      self%prev_step_maxabs = maxval(abs(gamma_new - gamma))
      self%prev_alpha = self%params%damping
   end subroutine picard_outer_step

end module moist_model_moz_solver_picard
