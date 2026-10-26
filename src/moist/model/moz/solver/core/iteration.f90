!> Shared outer-iteration control, independent of grid and field representation
module moist_model_moz_solver_core_iteration
   use mctc_env, only: wp, error_type
   use moist_context, only: moist_context_type
   use moist_utils_timer, only: cat_solve
   implicit none(type, external)
   private
   public :: moist_model_moz_solver_iteration_type, run_iterations

   !> Solve-local hooks; concrete contexts own scratch and borrow their problem
   type, abstract :: moist_model_moz_solver_iteration_type
   contains
      procedure(advance_interface), deferred :: advance
      procedure(observe_interface), deferred :: observe
      procedure(checkpoint_interface), deferred :: checkpoint
   end type moist_model_moz_solver_iteration_type

   !> Abstract interfaces of the deferred iteration hooks
   abstract interface
      !> Compute and accept one update, including driver-specific postprocessing
      !>
      !> @param[inout] self Solve-local context
      !> @param[in] iter Stage-local attempted iteration number
      !> @param[out] rmse Solver's unchanged convergence measure
      !> @param[out] feasible Whether a usable update was accepted
      !> @param[out] error Failure information
      subroutine advance_interface(self, iter, rmse, feasible, error)
         import :: moist_model_moz_solver_iteration_type, wp, error_type
         implicit none(type, external)
         !> Solve-local context
         class(moist_model_moz_solver_iteration_type), intent(inout) :: self
         !> Stage-local iteration number
         integer, intent(in) :: iter
         !> Convergence measure
         real(wp), intent(out) :: rmse
         !> Whether the step succeeded
         logical, intent(out) :: feasible
         !> Failure information
         type(error_type), allocatable, intent(out) :: error
      end subroutine advance_interface

      !> Report an accepted update before testing convergence
      !>
      !> @param[inout] self Solve-local diagnostic state
      !> @param[in] iter Stage-local iteration number
      !> @param[in] rmse Convergence measure from advance
      !> @param[in] elapsed Advance time, s; excludes diagnostics and I/O
      subroutine observe_interface(self, iter, rmse, elapsed)
         import :: moist_model_moz_solver_iteration_type, wp
         implicit none(type, external)
         !> Solve-local context
         class(moist_model_moz_solver_iteration_type), intent(inout) :: self
         !> Stage-local iteration number
         integer, intent(in) :: iter
         !> Convergence measure
         real(wp), intent(in) :: rmse
         !> Advance time, s; excludes diagnostics and I/O
         real(wp), intent(in) :: elapsed
      end subroutine observe_interface

      !> Write a checkpoint and restore the solver representation even on error
      !>
      !> @param[inout] self Solve-local checkpoint context
      !> @param[in] iter Stage-local iteration number
      !> @param[in] rmse Convergence measure from advance
      !> @param[out] error Failure information
      subroutine checkpoint_interface(self, iter, rmse, error)
         import :: moist_model_moz_solver_iteration_type, wp, error_type
         implicit none(type, external)
         !> Solve-local context
         class(moist_model_moz_solver_iteration_type), intent(inout) :: self
         !> Stage-local iteration number
         integer, intent(in) :: iter
         !> Convergence measure
         real(wp), intent(in) :: rmse
         !> Failure information
         type(error_type), allocatable, intent(out) :: error
      end subroutine checkpoint_interface
   end interface
contains
   !> Run one stage with driver-selected checkpoint/convergence ordering
   !>
   !> - Failed attempts count as iterations but are not observed or checkpointed
   !> - Only advance is timed, on the context timer node `Iterations` under the
   !>   open node; observe receives the last interval
   !> - The node is open on the timer stack during advance, so the solver's
   !>   residual and Jacobian-product timers nest under it
   !>
   !> @param[inout] context Driver-specific operations and borrowed problem
   !> @param[inout] ctx Run context whose timer records the advance time
   !> @param[in] max_iter Maximum attempted steps in this stage
   !> @param[in] tolerance Strict upper bound for convergence
   !> @param[in] checkpoint_interval Positive frequency; nonpositive disables I/O
   !> @param[inout] rmse Last residual; unchanged when no steps are attempted
   !> @param[out] niter Number of attempted steps in this stage
   !> @param[out] converged Whether an accepted step met the tolerance
   !> @param[out] error Step or checkpoint failure information
   !> @param[in] checkpoint_before_convergence Also checkpoint converged steps; checkpoint failure takes precedence
   recursive subroutine run_iterations(context, ctx, max_iter, tolerance, checkpoint_interval, &
                                       rmse, niter, converged, error, checkpoint_before_convergence)
      !> Solve-local context; no field arrays are owned by the controller
      class(moist_model_moz_solver_iteration_type), intent(inout) :: context
      !> Run context whose timer records the advance time
      type(moist_context_type), intent(inout) :: ctx
      !> Stage limit and checkpoint frequency
      integer, intent(in) :: max_iter, checkpoint_interval
      !> Strict convergence threshold
      real(wp), intent(in) :: tolerance
      !> Last computed residual, retained on an empty stage
      real(wp), intent(inout) :: rmse
      !> Attempted iteration count
      integer, intent(out) :: niter
      !> Successful convergence flag
      logical, intent(out) :: converged
      !> Failure information
      type(error_type), allocatable, intent(out) :: error
      !> Preserve drivers that checkpoint before testing convergence
      logical, optional, intent(in) :: checkpoint_before_convergence
      !> Stage-local iteration
      integer :: iter
      !> Advance timer node while it is open; 0 when the timer is inactive
      integer :: id
      !> Step feasibility, convergence and checkpoint ordering flags
      logical :: feasible, reached_tolerance, checkpoint_first

      niter = 0
      converged = .false.
      checkpoint_first = .false.
      if (present(checkpoint_before_convergence)) checkpoint_first = checkpoint_before_convergence
      do iter = 1, max_iter
         niter = iter
         call ctx%timer%start("Iterations", category=cat_solve)
         id = ctx%timer%current()
         call context%advance(iter, rmse, feasible, error)
         ! Stop before the failure exit so the node is never left running
         call ctx%timer%stop("Iterations")
         if (allocated(error) .or. .not. feasible) exit
         call context%observe(iter, rmse, ctx%timer%last(id))
         reached_tolerance = rmse < tolerance
         ! Fortran does not guarantee short-circuit evaluation of logical AND
         if (checkpoint_interval > 0 .and. (.not. reached_tolerance .or. checkpoint_first)) then
            if (mod(iter, checkpoint_interval) == 0) then
               call context%checkpoint(iter, rmse, error)
               if (allocated(error)) exit
            end if
         end if
         if (reached_tolerance) then
            converged = .true.
            exit
         end if
      end do
   end subroutine run_iterations
end module moist_model_moz_solver_core_iteration
