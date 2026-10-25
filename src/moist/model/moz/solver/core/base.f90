!> Shared interface for MOZ fixed-point and Newton-Krylov solvers
module moist_model_moz_solver_core_base
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_model_moz_solver_core_adjoint, only: moist_model_moz_solver_jfnk_adjoint_parameters_type
   use moist_model_moz_solver_core_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_moz_solver_core_printer, only: moist_model_moz_solver_printer_type
   use moist_model_moz_solver_core_preconditioner, only: moist_model_moz_solver_preconditioner_type
   use moist_model_moz_solver_jfnk_workspace, only: moist_model_moz_solver_workspace_type, exchange_solver_workspace
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_type, moist_model_moz_solver_statistics_type
   public :: validate_solver_arrays

   !> Solver counters; negative values mark unavailable metrics
   ! TODO: this IMO can be added to the abstract solver type (?)
   ! TODO: the krylov/gmres specific ones should then be moved to the jfnk solver type (?)
   type :: moist_model_moz_solver_statistics_type
      !> Stored DIIS vectors and GMRES work of the last update (including recycled refreshes)
      integer :: ndiis = 0, ngmres = 0
      !> Cumulative GMRES work, residual evaluations and exact Jacobian-product counts
      integer :: nkrylov = -1, nresidual = -1, njvp = -1
   end type moist_model_moz_solver_statistics_type

   !> One numerical update of an iterate supplied by a model driver
   type, abstract :: moist_model_moz_solver_type
      !> Borrowed run context, valid while the owning model lives
      type(moist_context_type), pointer :: ctx => null()
      !> Maximum absolute change and effective step fraction of the last update
      real(wp) :: prev_step_maxabs = 0.0_wp, prev_alpha = 0.0_wp
      !> Owned scratch; also caches transpose storage for non-Krylov primal solvers
      type(moist_model_moz_solver_workspace_type), allocatable :: workspace
      !> Controls of transpose solves that borrow this solver's scratch
      !>
      !> - The only copy: new_solver sets it from jfnk.adjoint for every solver kind
      !> - Defaults for solvers built by their own constructors
      type(moist_model_moz_solver_jfnk_adjoint_parameters_type) :: adjoint
   contains
      procedure(outer_step_interface), deferred :: outer_step
      procedure(name_interface), deferred :: name
      procedure :: allows_postprocessing => solver_allows_postprocessing
      procedure :: begin_stage => solver_begin_stage
      procedure :: set_iteration_printer => solver_set_iteration_printer
      procedure :: uses_preconditioner => solver_uses_preconditioner
      procedure :: set_preconditioner => solver_set_preconditioner
      procedure :: set_fallback_filter => solver_set_fallback_filter
      procedure :: get_statistics => solver_get_statistics
      !> One-line note on the fallbacks of the current stage; empty when there were none
      procedure :: stage_summary => solver_stage_summary
      !> Print one outer row with the solver's step, DIIS and GMRES columns
      procedure, non_overridable :: print_row => solver_print_row
      !> Print the work counts since a statistics snapshot below the table
      procedure, non_overridable :: print_statistics => solver_print_statistics
      !> OpenMP team size for this solver's own parallel regions (>= 1)
      procedure :: team_size => solver_team_size
      !> Transfer scratch ownership without copying arrays or solver settings
      procedure :: exchange_workspace => solver_exchange_workspace
   end type moist_model_moz_solver_type

   !> Abstract interfaces of the deferred solver bindings
   abstract interface
      !> Evaluate the residual and propose an update; the driver accepts the result
      !>
      !> @param[inout] self Solver state and workspace
      !> @param[in] residual_fn Model residual evaluator
      !> @param[in] gamma Current field (npts, npair)
      !> @param[out] gamma_new Proposed field of the same shape
      !> @param[out] rmse RMS residual of the input field
      !> @param[out] feasible Whether the proposed update is usable
      !> @param[out] error Failure information
      subroutine outer_step_interface(self, residual_fn, gamma, gamma_new, rmse, feasible, error)
         import :: moist_model_moz_solver_type, moist_model_moz_solver_evaluator_type, wp, error_type
         implicit none(type, external)
         !> Solver state and workspace
         class(moist_model_moz_solver_type), intent(inout) :: self
         !> Model residual evaluator
         class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_fn
         !> Current field (npts, npair)
         real(wp), intent(in) :: gamma(:, :)
         !> Proposed field (npts, npair)
         real(wp), intent(out) :: gamma_new(:, :)
         !> RMS residual of the input field
         real(wp), intent(out) :: rmse
         !> Whether the proposed update is usable
         logical, intent(out) :: feasible
         !> Failure information
         type(error_type), allocatable, intent(out) :: error
      end subroutine outer_step_interface

      !> Display name of the solver
      !>
      !> @param[in] self Solver to name
      pure function name_interface(self) result(text)
         import :: moist_model_moz_solver_type
         implicit none(type, external)
         !> Solver to name
         class(moist_model_moz_solver_type), intent(in) :: self
         !> Display name
         character(len=:), allocatable :: text
      end function name_interface
   end interface

contains

   !> Exchange scratch between completed solve stages; retain shared capacity
   !>
   !> @param[inout] self Solver with its settings and counters unchanged
   !> @param[inout] workspace Scratch owned by a distinct solver or solve context
   subroutine solver_exchange_workspace(self, workspace)
      !> Solver with its settings and counters unchanged
      class(moist_model_moz_solver_type), intent(inout) :: self
      !> Scratch owned by a distinct solver or solve context
      type(moist_model_moz_solver_workspace_type), allocatable, intent(inout) :: workspace
      call exchange_solver_workspace(self%workspace, workspace)
   end subroutine solver_exchange_workspace

   !> Check field shapes before evaluating a residual or allocating workspace
   !>
   !> @param[in] gamma Nonempty input field (npts, npair)
   !> @param[out] gamma_new Output storage; only its shape is inspected
   !> @param[out] error Empty input or mismatched output shape
   subroutine validate_solver_arrays(gamma, gamma_new, error)
      !> Input field
      real(wp), intent(in) :: gamma(:, :)
      !> Output field of the same shape
      real(wp), intent(out) :: gamma_new(:, :)
      !> Invalid field shape
      type(error_type), allocatable, intent(out) :: error

      if (size(gamma) == 0) then
         call fatal_error(error, "Solver fields must be nonempty")
      else if (any(shape(gamma) /= shape(gamma_new))) then
         call fatal_error(error, "Solver fields must have matching shapes")
      end if
   end subroutine validate_solver_arrays

   !> Whether the driver may filter or damp the returned update
   !>
   !> @param[in] self Solver to query
   pure logical function solver_allows_postprocessing(self) result(allowed)
      !> Solver to query
      class(moist_model_moz_solver_type), intent(in) :: self
      allowed = .true.
   end function solver_allows_postprocessing

   !> Begin a new equation or cold solve, clearing stage state and reusing workspace
   !>
   !> @param[inout] self Solver; the default has no stage-local state
   subroutine solver_begin_stage(self)
      !> Solver with no stage-local state by default
      class(moist_model_moz_solver_type), intent(inout) :: self
   end subroutine solver_begin_stage

   !> Borrow a printer whose target outlives the attachment; omit to detach
   !>
   !> - The default keeps no printer and shows neither the DIIS column nor the GMRES block
   !>
   !> @param[inout] self Solver to configure
   !> @param[inout] printer Borrowed printer
   subroutine solver_set_iteration_printer(self, printer)
      !> Solver to configure
      class(moist_model_moz_solver_type), intent(inout) :: self
      !> Borrowed printer
      type(moist_model_moz_solver_printer_type), target, optional, intent(inout) :: printer

      associate (s => self); end associate
      if (present(printer)) call printer%set_columns(diis=.false., gmres=.false.)
   end subroutine solver_set_iteration_printer

   !> Whether this solver supports a borrowed preconditioner
   !>
   !> @param[in] self Solver to query
   logical function solver_uses_preconditioner(self) result(uses)
      !> Solver to query
      class(moist_model_moz_solver_type), intent(in) :: self
      uses = .false.
   end function solver_uses_preconditioner

   !> Borrow a left preconditioner whose target outlives the attachment; omit to detach
   !>
   !> @param[inout] self Solver to configure; the default ignores preconditioning
   !> @param[inout] precond Borrowed operator and its scratch
   subroutine solver_set_preconditioner(self, precond)
      !> Solver to configure
      class(moist_model_moz_solver_type), intent(inout) :: self
      !> Borrowed operator and its scratch
      class(moist_model_moz_solver_preconditioner_type), target, optional, intent(inout) :: precond
   end subroutine solver_set_preconditioner

   !> Borrow an operator for the Picard fallback; omit to detach
   !>
   !> @param[inout] self Solver to configure; the default ignores the filter
   !> @param[inout] filter Borrowed operator and its scratch
   subroutine solver_set_fallback_filter(self, filter)
      !> Solver to configure
      class(moist_model_moz_solver_type), intent(inout) :: self
      !> Borrowed operator and its scratch
      class(moist_model_moz_solver_preconditioner_type), target, optional, intent(inout) :: filter
   end subroutine solver_set_fallback_filter

   !> Read solver counters without allocating diagnostic storage
   !>
   !> @param[in] self Solver to inspect
   function solver_get_statistics(self) result(stats)
      !> Solver to inspect
      class(moist_model_moz_solver_type), intent(in) :: self
      !> Statistics snapshot
      type(moist_model_moz_solver_statistics_type) :: stats
      stats = moist_model_moz_solver_statistics_type()
   end function solver_get_statistics

   !> Note on the fallbacks of the current stage, for the iteration table
   !>
   !> - Drivers print a nonempty note once after each stage
   !> - The default has no stage counters and returns an empty string
   !>
   !> @param[in] self Solver after a stage
   pure function solver_stage_summary(self) result(text)
      !> Solver after a stage
      class(moist_model_moz_solver_type), intent(in) :: self
      !> Note text; empty without fallbacks
      character(len=:), allocatable :: text

      associate (s => self); end associate
      text = ""
   end function solver_stage_summary

   !> Print one outer row: the model's diagnostics with the solver's columns
   !>
   !> @param[in]    self     Solver after an accepted step
   !> @param[inout] printer  Iteration table of the stage
   !> @param[in]    iter     Stage-local iteration
   !> @param[in]    elapsed  Advance time, s
   !> @param[in]    rmse     RMS residual
   !> @param[in]    hmax     max(h) column; absent leaves it blank
   !> @param[in]    energy   energy column; absent leaves it and dF blank
   !> @param[in]    lambda   lambda column of staged solves; absent leaves it blank
   subroutine solver_print_row(self, printer, iter, elapsed, rmse, hmax, energy, lambda)
      !> Solver after an accepted step
      class(moist_model_moz_solver_type), intent(in) :: self
      !> Iteration table of the stage
      type(moist_model_moz_solver_printer_type), intent(inout) :: printer
      !> Stage-local iteration
      integer, intent(in) :: iter
      !> Advance time and RMS residual
      real(wp), intent(in) :: elapsed, rmse
      !> Optional model diagnostics
      real(wp), intent(in), optional :: hmax, energy, lambda
      type(moist_model_moz_solver_statistics_type) :: stats

      stats = self%get_statistics()
      call printer%print_outer_row(iter, elapsed, rmse, self%prev_step_maxabs, hmax=hmax, lambda=lambda, &
         & energy=energy, alpha=self%prev_alpha, ndiis=stats%ndiis, gmres_total=stats%ngmres)
   end subroutine solver_print_row

   !> Print the Krylov, residual and Jacobian-product counts since a snapshot
   !>
   !> - Nothing for a solver without cumulative counters
   !>
   !> @param[in]    self     Solver after the solve
   !> @param[inout] printer  Iteration table of the solve
   !> @param[in]    since    Statistics taken before the solve
   subroutine solver_print_statistics(self, printer, since)
      !> Solver after the solve
      class(moist_model_moz_solver_type), intent(in) :: self
      !> Iteration table of the solve
      type(moist_model_moz_solver_printer_type), intent(inout) :: printer
      !> Statistics taken before the solve
      type(moist_model_moz_solver_statistics_type), intent(in) :: since
      type(moist_model_moz_solver_statistics_type) :: stats

      stats = self%get_statistics()
      if (stats%nresidual < 0) return
      call printer%stat_line("Krylov iterations", stats%nkrylov - since%nkrylov)
      call printer%stat_line("Residual evaluations", stats%nresidual - since%nresidual)
      call printer%stat_line("Exact Jacobian products", stats%njvp - since%njvp)
      call printer%blank()
   end subroutine solver_print_statistics

   !> OpenMP team size from the borrowed run context
   !>
   !> @param[in] self Solver to query
   function solver_team_size(self) result(nt)
      !> Solver to query
      class(moist_model_moz_solver_type), intent(in) :: self
      !> Team size (>= 1)
      integer :: nt
      nt = self%ctx%get_num_threads()
   end function solver_team_size

end module moist_model_moz_solver_core_base
