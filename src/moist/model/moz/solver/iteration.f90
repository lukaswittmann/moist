!> Shared outer-iteration control, independent of grid and field representation
module moist_model_moz_solver_iteration
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, error_type
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
      !> @param[in] elapsed Advance time in seconds, excluding diagnostics/I/O
      subroutine observe_interface(self, iter, rmse, elapsed)
         import :: moist_model_moz_solver_iteration_type, wp
         implicit none(type, external)
         !> Solve-local context
         class(moist_model_moz_solver_iteration_type), intent(inout) :: self
         !> Stage-local iteration number
         integer, intent(in) :: iter
         !> Convergence measure and advance wall time
         real(wp), intent(in) :: rmse, elapsed
      end subroutine observe_interface

      !> Write a checkpoint and restore the solver representation even on error
      !>
      !> @param[inout] self Solve-local checkpoint context
      !> @param[out] error Failure information
      subroutine checkpoint_interface(self, error)
         import :: moist_model_moz_solver_iteration_type, error_type
         implicit none(type, external)
         !> Solve-local context
         class(moist_model_moz_solver_iteration_type), intent(inout) :: self
         !> Failure information
         type(error_type), allocatable, intent(out) :: error
      end subroutine checkpoint_interface
   end interface
contains
   !> Run one stage with driver-selected checkpoint/convergence ordering
   !> Failed attempts count as iterations but are not observed or checkpointed
   !>
   !> @param[inout] context Driver-specific operations and borrowed problem
   !> @param[in] max_iter Maximum attempted steps in this stage
   !> @param[in] tolerance Strict upper bound for convergence
   !> @param[in] checkpoint_interval Positive frequency; nonpositive disables I/O
   !> @param[inout] rmse Last residual; unchanged when no steps are attempted
   !> @param[out] niter Number of attempted steps in this stage
   !> @param[out] converged Whether an accepted step met the tolerance
   !> @param[out] error Step or checkpoint failure information
   !> @param[in] checkpoint_before_convergence Also checkpoint converged steps;
   !> checkpoint failure takes precedence over convergence
   recursive subroutine run_iterations(context, max_iter, tolerance, checkpoint_interval, &
                                       rmse, niter, converged, error, checkpoint_before_convergence)
      !> Solve-local context; no field arrays are owned by the controller
      class(moist_model_moz_solver_iteration_type), intent(inout) :: context
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
      !> Stage-local iteration and monotonic wall-clock counters
      integer :: iter
      integer(int64) :: clock_start, clock_end, clock_rate
      !> Time spent advancing and step feasibility
      real(wp) :: elapsed
      logical :: feasible, reached_tolerance, checkpoint_first

      niter = 0
      converged = .false.
      checkpoint_first = .false.
      if (present(checkpoint_before_convergence)) checkpoint_first = checkpoint_before_convergence
      call system_clock(count_rate=clock_rate)
      do iter = 1, max_iter
         call system_clock(count=clock_start)
         niter = iter
         call context%advance(iter, rmse, feasible, error)
         if (allocated(error) .or. .not. feasible) exit
         call system_clock(count=clock_end)
         elapsed = 0.0_wp
         if (clock_rate > 0_int64) elapsed = real(clock_end - clock_start, wp)/real(clock_rate, wp)
         call context%observe(iter, rmse, elapsed)
         reached_tolerance = rmse < tolerance
         ! Fortran does not guarantee short-circuit evaluation of logical AND
         if (checkpoint_interval > 0 .and. (.not. reached_tolerance .or. checkpoint_first)) then
            if (mod(iter, checkpoint_interval) == 0) then
               call context%checkpoint(error)
               if (allocated(error)) exit
            end if
         end if
         if (reached_tolerance) then
            converged = .true.
            exit
         end if
      end do
   end subroutine run_iterations
end module moist_model_moz_solver_iteration
