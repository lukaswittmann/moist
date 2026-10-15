!> MDIIS startup followed by Newton-Krylov refinement
!>
!> The RMS residual of the input iterate selects a one-way switch to JFNK.
!> Each solve or coupling stage starts in MDIIS. Both phases use the shared
!> evaluator interface; field representation and closure evaluation belong to the model.
module moist_model_moz_solver_hybrid
   use mctc_env, only: wp, error_type
   use moist_model_moz_solver_base, only: &
      & moist_model_moz_solver_type, moist_model_moz_solver_statistics_type, validate_solver_arrays
   use moist_model_moz_solver_diis, only: moist_model_moz_solver_diis_type
   use moist_model_moz_solver_jfnk, only: moist_model_moz_solver_jfnk_type
   use moist_model_moz_solver_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_moz_solver_printer, only: moist_model_moz_solver_printer_type
   use moist_model_moz_solver_preconditioner, only: moist_model_moz_solver_preconditioner_type
   use moist_model_moz_solver_workspace, only: moist_model_moz_solver_workspace_type
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_hybrid_type, new_hybrid_solver

   !> Compose the existing MDIIS and JFNK implementations with a residual switch
   type, extends(moist_model_moz_solver_type) :: moist_model_moz_solver_hybrid_type
      private
      !> Solver state and workspaces owned by this solve context
      type(moist_model_moz_solver_diis_type) :: mdiis
      type(moist_model_moz_solver_jfnk_type) :: jfnk
      !> Switch when the input RMS residual is at or below this value
      real(wp) :: switch_rmse = 1.0e-3_wp
      !> Remain in JFNK until begin_stage is called
      logical :: jfnk_active = .false.
      !> Residual evaluations made by MDIIS, including the handoff evaluation
      integer :: nmdiis_residual = 0
   contains
      procedure :: outer_step => hybrid_outer_step
      procedure :: name => hybrid_name
      procedure :: begin_stage => hybrid_begin_stage
      procedure :: allows_postprocessing => hybrid_allows_postprocessing
      procedure :: get_statistics => hybrid_get_statistics
      procedure :: uses_preconditioner => hybrid_uses_preconditioner
      procedure :: set_preconditioner => hybrid_set_preconditioner
      procedure :: set_fallback_filter => hybrid_set_fallback_filter
      procedure :: set_iteration_printer => hybrid_set_iteration_printer
      procedure :: exchange_workspace => hybrid_exchange_workspace
   end type moist_model_moz_solver_hybrid_type

contains

   !> Lend the embedded Newton solver's scratch between completed solve stages
   !>
   !> @param[inout] self Composite with mixing state and numerical controls retained
   !> @param[inout] workspace Owned scratch of a distinct solve context
   subroutine hybrid_exchange_workspace(self, workspace)
      class(moist_model_moz_solver_hybrid_type), intent(inout) :: self
      type(moist_model_moz_solver_workspace_type), allocatable, intent(inout) :: workspace
      call self%jfnk%exchange_workspace(workspace)
   end subroutine hybrid_exchange_workspace

   !> Construct a hybrid from configured MDIIS and JFNK solvers
   !>
   !> @param[in] mdiis Startup solver with its mixing controls
   !> @param[in] jfnk Refinement solver with its Krylov and line-search controls
   !> @param[in] switch_rmse Positive finite RMS threshold; defaults to 1e-3
   !> @return Hybrid solver starting in MDIIS
   function new_hybrid_solver(mdiis, jfnk, switch_rmse) result(self)
      !> Configured startup solver
      type(moist_model_moz_solver_diis_type), intent(in) :: mdiis
      !> Configured refinement solver
      type(moist_model_moz_solver_jfnk_type), intent(in) :: jfnk
      !> Positive finite switch threshold
      real(wp), optional, intent(in) :: switch_rmse
      !> Composed solver with fresh stage state
      type(moist_model_moz_solver_hybrid_type) :: self

      self%mdiis = mdiis
      self%jfnk = jfnk
      self%verbosity = mdiis%verbosity
      if (present(switch_rmse)) self%switch_rmse = switch_rmse
      call self%begin_stage()
   end function new_hybrid_solver

   !> Display name shared by both phases
   !>
   !> @param[in] self Solver to name
   pure function hybrid_name(self) result(text)
      !> Solver to name
      class(moist_model_moz_solver_hybrid_type), intent(in) :: self
      !> Display name
      character(len=:), allocatable :: text

      text = "HYBRID"
   end function hybrid_name

   !> Propose an MDIIS update or hand the current iterate to JFNK
   !>
   !> The switch uses the residual of gamma. On the handoff call, JFNK
   !> evaluates the same gamma and its accepted update replaces the MDIIS proposal.
   !> This costs one additional residual evaluation only at the handoff.
   !>
   !> @param[inout] self Solver phase, histories and workspaces
   !> @param[in] residual_fn Model residual evaluator
   !> @param[in] gamma Current field (npts, npair)
   !> @param[out] gamma_new Proposed field with the same shape
   !> @param[out] rmse RMS residual of the input field
   !> @param[out] feasible Whether the proposed update is usable
   !> @param[out] error Evaluation failure
   subroutine hybrid_outer_step(self, residual_fn, gamma, gamma_new, rmse, feasible, error)
      !> Solver phase, histories and workspaces
      class(moist_model_moz_solver_hybrid_type), intent(inout) :: self
      !> Model residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_fn
      !> Input field (npts, npair)
      real(wp), intent(in) :: gamma(:, :)
      !> Proposed field with the same shape
      real(wp), intent(out) :: gamma_new(:, :)
      !> RMS residual of gamma
      real(wp), intent(out) :: rmse
      !> Whether the update is usable
      logical, intent(out) :: feasible
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      feasible = .false.
      call validate_solver_arrays(gamma, gamma_new, error)
      if (allocated(error)) return
      ! The private components follow the composite's thread count, also when
      ! a model sets solver%nthreads after construction
      self%mdiis%nthreads = self%nthreads
      self%jfnk%nthreads = self%nthreads
      if (.not. self%jfnk_active) then
         call self%mdiis%outer_step(residual_fn, gamma, gamma_new, rmse, feasible, error)
         self%nmdiis_residual = self%nmdiis_residual + 1
         if (.not. feasible .or. allocated(error)) return
         if (rmse <= self%switch_rmse) then
            self%jfnk_active = .true.
            if (self%verbosity >= 2 .and. associated(self%mdiis%iter_printer)) then
               call self%mdiis%iter_printer%note("Switching from MDIIS to JFNK")
            end if
            ! Release the MDIIS history before JFNK allocates its Krylov workspace.
            call self%mdiis%reset()
            call self%jfnk%begin_stage()
         else
            self%prev_step_maxabs = self%mdiis%prev_step_maxabs
            self%prev_alpha = self%mdiis%prev_alpha
            return
         end if
      end if

      call self%jfnk%outer_step(residual_fn, gamma, gamma_new, rmse, feasible, error)
      if (.not. feasible .or. allocated(error)) return
      self%prev_step_maxabs = self%jfnk%prev_step_maxabs
      self%prev_alpha = self%jfnk%prev_alpha
   end subroutine hybrid_outer_step

   !> Start again in MDIIS with an empty history and reset JFNK stage state
   !>
   !> JFNK workspace and borrowed operators are retained across stages.
   !>
   !> @param[inout] self Solver starting a new stage
   subroutine hybrid_begin_stage(self)
      !> Solver to reset
      class(moist_model_moz_solver_hybrid_type), intent(inout) :: self

      call self%mdiis%reset()
      call self%jfnk%begin_stage()
      self%jfnk_active = .false.
      self%prev_step_maxabs = 0.0_wp
      self%prev_alpha = 0.0_wp
   end subroutine hybrid_begin_stage

   !> Allow driver damping/filtering only while the returned update uses MDIIS
   !>
   !> @param[in] self Solver after proposing an update
   pure logical function hybrid_allows_postprocessing(self) result(allowed)
      !> Solver phase to query
      class(moist_model_moz_solver_hybrid_type), intent(in) :: self

      allowed = .not. self%jfnk_active
   end function hybrid_allows_postprocessing

   !> Report the current phase and cumulative residual work of both solvers
   !>
   !> @param[in] self Solver to inspect
   !> @return Combined statistics
   function hybrid_get_statistics(self) result(stats)
      !> Solver to inspect
      class(moist_model_moz_solver_hybrid_type), intent(in) :: self
      !> Current history and cumulative counters
      type(moist_model_moz_solver_statistics_type) :: stats

      stats = self%jfnk%get_statistics()
      stats%nresidual = stats%nresidual + self%nmdiis_residual
      stats%ndiis = self%mdiis%nstored
      if (.not. self%jfnk_active) stats%ngmres = 0
   end function hybrid_get_statistics

   !> Expose JFNK's capability before the switch so the driver can attach operators
   !>
   !> @param[in] self Solver to query
   logical function hybrid_uses_preconditioner(self) result(uses)
      !> Composed solver
      class(moist_model_moz_solver_hybrid_type), intent(in) :: self

      uses = self%jfnk%uses_preconditioner()
   end function hybrid_uses_preconditioner

   !> Borrow a left preconditioner for JFNK; omit to detach
   !>
   !> The target must outlive its attachment to this solver.
   !>
   !> @param[inout] self Solver to configure
   !> @param[inout] precond Borrowed operator and its scratch
   subroutine hybrid_set_preconditioner(self, precond)
      !> Solver to configure
      class(moist_model_moz_solver_hybrid_type), intent(inout) :: self
      !> Borrowed operator
      class(moist_model_moz_solver_preconditioner_type), target, optional, intent(inout) :: precond

      call self%jfnk%set_preconditioner(precond)
   end subroutine hybrid_set_preconditioner

   !> Borrow a filter for JFNK's Picard fallback; omit to detach
   !>
   !> The target must outlive its attachment to this solver.
   !>
   !> @param[inout] self Solver to configure
   !> @param[inout] filter Borrowed filter and its scratch
   subroutine hybrid_set_fallback_filter(self, filter)
      !> Solver to configure
      class(moist_model_moz_solver_hybrid_type), intent(inout) :: self
      !> Borrowed fallback filter
      class(moist_model_moz_solver_preconditioner_type), target, optional, intent(inout) :: filter

      call self%jfnk%set_fallback_filter(filter)
   end subroutine hybrid_set_fallback_filter

   !> Borrow the same printer for both phases; omit to detach both
   !>
   !> The target must outlive its attachment to this solver.
   !>
   !> @param[inout] self Solver to configure
   !> @param[in] printer Borrowed printer
   subroutine hybrid_set_iteration_printer(self, printer)
      !> Solver to configure
      class(moist_model_moz_solver_hybrid_type), intent(inout) :: self
      !> Borrowed printer
      type(moist_model_moz_solver_printer_type), target, optional, intent(in) :: printer

      call self%mdiis%set_iteration_printer(printer)
      call self%jfnk%set_iteration_printer(printer)
   end subroutine hybrid_set_iteration_printer

end module moist_model_moz_solver_hybrid
