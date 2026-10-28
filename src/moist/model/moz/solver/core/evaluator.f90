!> Residual evaluator for MOZ fixed-point and Newton-Krylov solvers
!>
!> - Concrete evaluators define the field representation
!> - Concrete evaluators own or borrow the problem and scratch
module moist_model_moz_solver_core_evaluator
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_moz_solver_jfnk_deflation, only: moist_model_moz_solver_deflation_type
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_evaluator_type

   !> Residual and optional Jacobian-vector operations supplied by the model
   type, abstract :: moist_model_moz_solver_evaluator_type
   contains
      !> MOZ residual R(gamma) = gamma - gamma_map(gamma)
      procedure(moz_evaluator_residual_iface), deferred :: residual
      !> Whether `linearize`/`jvp` provide exact Jacobian-vector products
      !> (default: no; a Newton-Krylov step with `use_exact_jvp` then fails)
      procedure :: has_jvp => moz_evaluator_has_jvp_default
      !> Prepare linearization data at gamma for later `jvp` calls (default: not supported)
      procedure :: linearize => moz_evaluator_linearize_default
      !> Exact product J v at the last `linearize` point (default: not supported)
      procedure :: jvp => moz_evaluator_jvp_default
      !> Whether a transpose product is available at the frozen linearization
      procedure :: has_adjoint => moz_evaluator_has_adjoint_default
      !> Euclidean transpose product J**T v (default: not supported)
      procedure :: adjoint => moz_evaluator_adjoint_default
      !> Fill and factorize the coarse matrix of a deflation object at the last
      !> `linearize` point (default: not supported; the deflation stays inactive)
      procedure :: deflation_setup => moz_evaluator_deflation_setup_default
   end type moist_model_moz_solver_evaluator_type

   !> Abstract interface of the deferred residual binding
   abstract interface
      !> Evaluate the MOZ residual at the current iterate
      !>
      !> @param[in]  self      Evaluator (read-only operator; carries the problem)
      !> @param[in]  gamma     Current field iterate (npts, npair)
      !> @param[out] residual  Residual R(gamma) = gamma - gamma_map (npts, npair)
      !> @param[out] feasible  Whether the evaluation produced a usable state
      !> @param[out] error     Error handle on failure
      subroutine moz_evaluator_residual_iface(self, gamma, residual, feasible, error)
         import :: moist_model_moz_solver_evaluator_type, wp, error_type
         implicit none(type, external)
         !> Evaluator (carries the problem context)
         class(moist_model_moz_solver_evaluator_type), intent(in) :: self
         !> Current iterate (npts, npair)
         real(wp), intent(in) :: gamma(:, :)
         !> Residual R(gamma) (npts, npair)
         real(wp), intent(out) :: residual(:, :)
         !> Whether the residual is usable
         logical, intent(out) :: feasible
         !> Error handle
         type(error_type), allocatable, intent(out) :: error
      end subroutine moz_evaluator_residual_iface
   end interface

contains

   !> Default: no exact Jacobian-vector products; result always false
   !>
   !> @param[in] self  evaluator
   logical function moz_evaluator_has_jvp_default(self) result(available)
      !> Evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: self

      associate (s => self); end associate
      available = .false.
   end function moz_evaluator_has_jvp_default

   !> Default linearization: not supported by this evaluator
   !>
   !> @param[in]  self      Evaluator
   !> @param[in]  gamma     Linearization point (npts, npair)
   !> @param[out] feasible  Always .false. here
   !> @param[out] error     Fatal error: exact Jacobian not implemented
   subroutine moz_evaluator_linearize_default(self, gamma, feasible, error)
      !> Evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: self
      !> Linearization point (npts, npair)
      real(wp), intent(in) :: gamma(:, :)
      !> Feasibility flag
      logical, intent(out) :: feasible
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      associate (s => self, g => gamma); end associate
      feasible = .false.
      call fatal_error(error, "This MOZ evaluator does not implement an exact Jacobian")
   end subroutine moz_evaluator_linearize_default

   !> Default Jacobian-vector product: not supported by this evaluator
   !>
   !> @param[in]  self       Evaluator
   !> @param[in]  direction  Direction v, flat (npts*npair)
   !> @param[out] jv         Product J v (set to 0 on the error path)
   !> @param[out] feasible   Always .false. here
   !> @param[out] error      Fatal error: exact Jacobian not implemented
   subroutine moz_evaluator_jvp_default(self, direction, jv, feasible, error)
      !> Evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: self
      !> Direction v, flat (npts*npair)
      real(wp), intent(in) :: direction(:)
      !> Product J v, flat (npts*npair)
      real(wp), intent(out) :: jv(:)
      !> Feasibility flag
      logical, intent(out) :: feasible
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      associate (s => self, d => direction); end associate
      jv = 0.0_wp
      feasible = .false.
      call fatal_error(error, "This MOZ evaluator does not implement an exact Jacobian")
   end subroutine moz_evaluator_jvp_default

   !> Default coarse-space setup: not supported, the deflation stays inactive
   !>
   !> - A singular coarse matrix leaves the deflation inactive without an error
   !>
   !> @param[in]    self       Evaluator
   !> @param[inout] deflation  Deflation object; left inactive
   !> @param[out]   error      Coarse-matrix failure; never set here
   subroutine moz_evaluator_deflation_setup_default(self, deflation, error)
      !> Evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: self
      !> Deflation object
      type(moist_model_moz_solver_deflation_type), intent(inout) :: deflation
      !> Coarse-matrix failure
      type(error_type), allocatable, intent(out) :: error

      associate (s => self); end associate
      call deflation%deactivate()
   end subroutine moz_evaluator_deflation_setup_default

   !> Default: the evaluator does not supply transpose products
   !>
   !> @param[in] self  evaluator
   logical function moz_evaluator_has_adjoint_default(self) result(available)
      !> Evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: self
      available = .false.
   end function moz_evaluator_has_adjoint_default

   !> Report a missing transpose operation
   !>
   !> @param[in]  self       evaluator
   !> @param[in]  direction  flat Euclidean direction
   !> @param[out] jtv        transpose product; defined only on success
   !> @param[out] feasible   always false here
   !> @param[out] error      fatal error: transpose product not provided
   subroutine moz_evaluator_adjoint_default(self, direction, jtv, feasible, error)
      !> Evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: self
      !> Flat Euclidean direction
      real(wp), intent(in) :: direction(:)
      !> Transpose product; defined only on success
      real(wp), intent(out) :: jtv(:)
      !> Product usability
      logical, intent(out) :: feasible
      !> Unsupported operation
      type(error_type), allocatable, intent(out) :: error
      feasible = .false.
      call fatal_error(error, "This evaluator does not provide a transpose Jacobian product")
   end subroutine moz_evaluator_adjoint_default

end module moist_model_moz_solver_core_evaluator
