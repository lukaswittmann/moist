!> Jacobian-free Newton-Krylov solver with restarted GMRES
!> Uses finite differences or evaluator-supplied products without forming a Jacobian
module moist_model_moz_solver_jfnk
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_moz_solver_base, only: moist_model_moz_solver_type, validate_solver_arrays
   use moist_model_moz_solver_base, only: moist_model_moz_solver_statistics_type
   use moist_model_moz_solver_printer, only: moist_model_moz_solver_printer_type
   use moist_model_moz_solver_preconditioner, only: moist_model_moz_solver_preconditioner_type
   use moist_math_orthogonalization, only: orthogonalize_mgs => moist_orthogonalize_mgs
   use moist_math_blas, only: gemv, gemv_rows, dot, gemm_rows_tn, gemm_rows_nn
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_lapack_potrf, only: wrap_potrf
   use moist_math_lapack_ggev, only: dggev
   use moist_model_moz_solver_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_moz_solver_deflation, only: moist_model_moz_solver_deflation_type
   use moist_model_moz_solver_workspace, only: moist_model_moz_solver_workspace_type, exchange_solver_workspace
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_jfnk_type, new_jfnk_solver

   !> Avoid OpenMP team overhead for small vectors
   integer, parameter :: parallel_vector_size = 65536
   !> Fixed reduction partitions, independent of the OpenMP team size
   integer, parameter :: reduction_block_size = 8192

   !> Jacobian-free Newton-Krylov solver with restarted GMRES
   !> Operates on fields stored as (npts, npair) arrays
   type, extends(moist_model_moz_solver_type) :: moist_model_moz_solver_jfnk_type
      !> Maximum Krylov dimension of one cycle; convergence or max_krylov may end it early
      !> Shared workspace may retain capacity for the larger of primal/transpose restarts
      integer :: restart = 50
      !> Inner iterations per Newton step; the default is one restart cycle
      !> The driver sets the outer iteration limit separately
      integer :: max_krylov = 50
      !> Minimum relative residual reduction per restart cycle to continue
      !>
      !> - below this fraction, the inner solve bails out early
      !> - default 0.0 disables the early-out: GMRES runs to the
      !>   inexact-Newton tolerance eta*||F|| or the iteration cap, avoiding a
      !>   premature bail on non-symmetric MOZ Jacobian stagnation
      !> - set > 0 (e.g. 0.1) to re-enable the early-out
      real(wp) :: restart_stagnation_tol = 0.0_wp
      real(wp) :: eta0 = 0.5_wp
      real(wp) :: eta_min = 1.0e-3_wp
      !> Upper bound of the forcing term from the third update of a stage on
      real(wp) :: eta_max = 0.1_wp
      real(wp) :: forcing_beta = 0.9_wp
      real(wp) :: forcing_power = 2.0_wp
      real(wp) :: line_search_c = 1.0e-4_wp
      real(wp) :: line_search_tau = 0.5_wp
      real(wp) :: line_search_min_alpha = 1.0e-4_wp
      !> Hard cap on the number of trial steps (residual evaluations) taken
      !> by the backtracking line search
      !>
      !> Rejected trials shrink alpha to [0.1, line_search_tau]*alpha.
      !> This cap may stop the search before line_search_min_alpha is reached;
      !> step extension and the Picard fallback have their own trial loops
      integer :: max_line_search_iter = 10
      real(wp) :: line_search_short_alpha = 0.25_wp
      !> Longest step the line search may try beyond the full Newton step
      real(wp) :: alpha_max = 2.0_wp
      !> Extend only when the full step keeps more than this fraction of the
      !> merit phi = 0.5*||F||^2; a stronger decrease is the fast Newton regime,
      !> where a longer step overshoots and its trial would be wasted
      real(wp) :: line_search_expand_ratio = 0.25_wp
      real(wp) :: fd_eps = 0.0_wp
      !> Use the evaluator's exact Jacobian-vector products where it provides
      !> them; .false. forces finite differences
      logical :: exact_jvp = .true.
      !> Whether the evaluator holds a linearization at the current Newton point
      logical :: jvp_ready = .false.
      !> Whether the last update was a zero step at eta_min; a second one is a stall
      logical :: stalled_at_min = .false.
      !> Requested recycled subspace dimension k of GCRO-DR (0 = off); the solver
      !> also caps it by restart and max_krylov to leave at least one Arnoldi step
      integer :: recycle = 0
      !> Number of recycled vectors currently held in ws_u (0 = none yet)
      integer :: nrecycled = 0
      !> Whether the recycled subspace was built for the preconditioned operator
      logical :: recycle_precond = .false.
      !> Use CGS2 (classical Gram-Schmidt + reorthogonalization) for the
      !> Arnoldi step instead of modified Gram-Schmidt
      !>
      !> - fewer vector sweeps, with speed depending on field size and threading
      !> - rounding can change the nonlinear convergence path, so it defaults off
      !> - selected with the cgs2 constructor argument
      logical :: use_cgs2 = .false.
      !> Newton state shared by outer_step calls within one solve stage
      !> A negative forcing term selects eta0 on the next update
      real(wp) :: eta_current = -1.0_wp
      logical :: precond_valid_state = .false.
      real(wp) :: rmse_prev_state = huge(1.0_wp)
      !> Last inner work count: Arnoldi steps plus recycled-image refresh products
      integer :: gmres_iters_last = 0
      !> Cumulative inner work in this stage, including recycled refreshes and failed solves
      integer :: gmres_iters_total = 0
      !> Cumulative inner work, residual evaluations (Newton residuals,
      !> line-search trials and finite-difference probes) and exact
      !> Jacobian-vector products since construction
      integer :: nkrylov_total = 0, nresidual_total = 0, njvp_total = 0
      !> Attempted outer_step calls in the current stage
      integer :: newton_iter_count = 0
      !> Pointer to the outer iteration printer (set by driver, nullified after use)
      type(moist_model_moz_solver_printer_type), pointer :: iter_printer => null()
      !> Left preconditioner lent by the driver for the current solve, if any
      class(moist_model_moz_solver_preconditioner_type), pointer :: precond => null()
      !> Filter lent by the driver to smooth the Picard fallback step, if any
      class(moist_model_moz_solver_preconditioner_type), pointer :: smoother => null()
      !> Coarse-space deflation (DEF1), filled by evaluators that support it
      type(moist_model_moz_solver_deflation_type) :: deflation
      !> Last reported preconditioner status and deflation state, for the notes
      character(len=:), allocatable :: precond_status
      logical :: deflation_status = .false.
      !> Cumulative wall seconds allocating and initializing workspace
      real(wp) :: workspace_wall = 0.0_wp
      !> Cumulative linear-solve algebra wall seconds, excluding JVP/preconditioner calls
      real(wp) :: vector_wall = 0.0_wp
   contains
      procedure :: get_statistics => jfnk_get_statistics
      procedure :: allows_postprocessing => jfnk_allows_postprocessing
      procedure :: begin_stage => jfnk_begin_stage
      procedure :: set_iteration_printer => jfnk_set_iteration_printer
      procedure :: uses_preconditioner => jfnk_uses_preconditioner
      procedure :: set_preconditioner => jfnk_set_preconditioner
      procedure :: set_fallback_filter => jfnk_set_fallback_filter
      procedure :: outer_step => jfnk_outer_step
      procedure :: name => jfnk_name
      procedure :: ensure_workspace => jfnk_ensure_workspace
      procedure :: release_workspace => jfnk_release_workspace
      procedure :: exchange_workspace => jfnk_exchange_workspace
      procedure :: linear_solve => gmres_linear_solve
      procedure :: apply_precond => jfnk_apply_precond
      procedure :: apply_jacobian => jfnk_apply_jacobian
      procedure :: recycle_dim => gmres_recycle_dim
      procedure :: recycle_refresh => gmres_recycle_refresh
      procedure :: recycle_harvest => gmres_recycle_harvest
      procedure :: deflate => gmres_deflate
      procedure :: deflation_finish => gmres_deflation_finish
      procedure, private :: note_status => jfnk_note_status
      procedure, private :: ensure_fd_workspace => jfnk_ensure_fd_workspace
      procedure, private :: ensure_recycle_workspace => jfnk_ensure_recycle_workspace
   end type moist_model_moz_solver_jfnk_type

contains

   !> Move scratch between phases and invalidate operator-dependent Krylov data
   !>
   !> @param[inout] self Solver; numerical controls and counters are retained
   !> @param[inout] workspace Owned scratch of a distinct solve context
   subroutine jfnk_exchange_workspace(self, workspace)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      type(moist_model_moz_solver_workspace_type), allocatable, intent(inout) :: workspace

      call exchange_solver_workspace(self%workspace, workspace)
      self%nrecycled = 0
      self%recycle_precond = .false.
      self%jvp_ready = .false.
      call self%deflation%deactivate()
   end subroutine jfnk_exchange_workspace

   !> Lazily allocate FD fields, including failed-exact-linearization fallback
   !>
   !> @param[inout] self Solver with an allocated workspace
   !> @param[in] gamma Field shape; products overwrite every scratch value
   subroutine jfnk_ensure_fd_workspace(self, gamma)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      real(wp), intent(in) :: gamma(:, :)
      integer(int64) :: tick

      if (allocated(self%workspace%ws_gamma_pert)) return
      call system_clock(tick)
      allocate (self%workspace%ws_gamma_pert, self%workspace%ws_f_pert, mold=gamma)
      call accumulate_wall(tick, self%workspace_wall)
   end subroutine jfnk_ensure_fd_workspace

   !> Grow optional recycling storage without replacing a shared Arnoldi basis
   !>
   !> @param[inout] self Solver with shape-compatible shared scratch
   !> @param[in] kmax Requested recycled dimension; zero does not consume storage
   subroutine jfnk_ensure_recycle_workspace(self, kmax)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      integer, intent(in) :: kmax
      integer :: i, j, nflat, m
      integer(int64) :: tick

      ! Drop an oversized active basis while retaining shared storage.
      if (self%nrecycled > kmax) self%nrecycled = 0
      if (kmax == 0) return
      if (allocated(self%workspace%ws_u)) then
         if (size(self%workspace%ws_u, 2) >= kmax) return
      end if
      call system_clock(tick)
      nflat = self%workspace%ws_nflat
      m = size(self%workspace%ws_v, 2) - 1
      if (allocated(self%workspace%ws_u)) deallocate (self%workspace%ws_u, self%workspace%ws_c)
      deallocate (self%workspace%ws_bproj)
      allocate (self%workspace%ws_u(nflat, kmax), self%workspace%ws_c(nflat, kmax))
      allocate (self%workspace%ws_bproj(kmax, m), source=0.0_wp)
      self%nrecycled = 0
      !$omp parallel do collapse(2) schedule(static) num_threads(merge(self%team_size(), 1, &
      !$omp & nflat >= parallel_vector_size))
      do j = 1, kmax
         do i = 1, nflat
            self%workspace%ws_u(i, j) = 0.0_wp
            self%workspace%ws_c(i, j) = 0.0_wp
         end do
      end do
      !$omp end parallel do
      call accumulate_wall(tick, self%workspace_wall)
   end subroutine jfnk_ensure_recycle_workspace

   !> Display name of the solver
   !>
   !> @param[in] self Solver to name
   pure function jfnk_name(self) result(text)
      !> Solver to name
      class(moist_model_moz_solver_jfnk_type), intent(in) :: self
      !> Display name
      character(len=:), allocatable :: text
      text = "JFNK"
   end function jfnk_name

   !> Preserve the update accepted by the Newton line search
   !>
   !> @param[in] self Solver to query
   pure logical function jfnk_allows_postprocessing(self) result(allowed)
      !> Solver to query
      class(moist_model_moz_solver_jfnk_type), intent(in) :: self
      allowed = .false.
   end function jfnk_allows_postprocessing

   !> Print a note when the preconditioner status or the deflation state changes
   !>
   !> - The first step of a table reports the preconditioner in any case
   !> - Needs trace verbosity >= 5 and an attached printer; the status is tracked anyway
   !>
   !> @param[inout] self Solver after deciding this step's operator and coarse space
   subroutine jfnk_note_status(self)
      !> Solver with its last reported status
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Status of this step
      character(len=:), allocatable :: status
      !> Formatted coarse-space size
      character(len=16) :: buf
      !> Whether notes are printed
      logical :: show

      show = self%verbosity >= 5 .and. associated(self%iter_printer)
      if (associated(self%precond)) then
         status = self%precond%describe(self%rmse_prev_state)
      else
         status = "none"
      end if
      if (allocated(self%precond_status)) then
         if (status /= self%precond_status) then
            if (show) call self%iter_printer%note("preconditioner "//status)
         end if
      else
         if (show) call self%iter_printer%note("preconditioner "//status)
      end if
      self%precond_status = status

      if (self%deflation%active .eqv. self%deflation_status) return
      self%deflation_status = self%deflation%active
      if (.not. show) return
      if (self%deflation%active) then
         write (buf, "(i0)") self%deflation%kz
         call self%iter_printer%note("deflation on, "//trim(buf)//" coarse vectors")
      else
         call self%iter_printer%note("deflation off")
      end if
   end subroutine jfnk_note_status

   !> Invalidate stage forcing, recycling and stall state; retain workspace and borrowed operators
   !> The next outer_step reinitializes the remaining Newton state; cumulative counters persist
   !>
   !> @param[inout] self Solver starting a solve or coupling stage
   subroutine jfnk_begin_stage(self)
      !> Solver retaining its workspaces
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      self%eta_current = -1.0_wp
      self%gmres_iters_total = 0
      self%nrecycled = 0
      self%stalled_at_min = .false.
   end subroutine jfnk_begin_stage

   !> Attach or detach a borrowed inner-iteration printer
   !>
   !> - A new printer starts a new table, so the next step reports the status again
   !>
   !> @param[inout] self Solver to configure
   !> @param[in] printer Target must outlive the attachment; omit to detach
   subroutine jfnk_set_iteration_printer(self, printer)
      !> Solver to configure
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Optional borrowed printer
      type(moist_model_moz_solver_printer_type), target, optional, intent(in) :: printer
      nullify(self%iter_printer)
      if (present(printer)) self%iter_printer => printer
      if (allocated(self%precond_status)) deallocate (self%precond_status)
      self%deflation_status = .false.
   end subroutine jfnk_set_iteration_printer

   !> Whether this solver supports a lent preconditioner; does not mean one is active
   !>
   !> @param[in] self Solver to query
   !> @return Always true; attachment and readiness are checked for each linear solve
   logical function jfnk_uses_preconditioner(self) result(uses)
      !> Solver to query
      class(moist_model_moz_solver_jfnk_type), intent(in) :: self
      uses = .true.
   end function jfnk_uses_preconditioner

   !> Lend or detach the left preconditioner of the linear solves
   !>
   !> @param[inout] self Solver to configure
   !> @param[inout] precond Target must outlive the attachment; omit to detach
   subroutine jfnk_set_preconditioner(self, precond)
      !> Solver to configure
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Optional borrowed operator
      class(moist_model_moz_solver_preconditioner_type), target, optional, intent(inout) :: precond
      nullify(self%precond)
      if (present(precond)) self%precond => precond
   end subroutine jfnk_set_preconditioner

   !> Lend or detach the filter smoothing the Picard fallback step
   !>
   !> @param[inout] self Solver to configure
   !> @param[inout] filter Target must outlive the attachment; omit to detach
   subroutine jfnk_set_fallback_filter(self, filter)
      !> Solver to configure
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Optional borrowed operator
      class(moist_model_moz_solver_preconditioner_type), target, optional, intent(inout) :: filter
      nullify(self%smoother)
      if (present(filter)) self%smoother => filter
   end subroutine jfnk_set_fallback_filter

   !> Return an allocation-free snapshot of solver statistics
   !>
   !> @param[in] self Solver to inspect
   !> @return stats Counters and cumulative timings, or unavailable sentinels
   function jfnk_get_statistics(self) result(stats)
      !> Solver to inspect
      class(moist_model_moz_solver_jfnk_type), intent(in) :: self
      !> Statistics snapshot
      type(moist_model_moz_solver_statistics_type) :: stats
      stats = moist_model_moz_solver_statistics_type()
      stats%ngmres = self%gmres_iters_last
      stats%nkrylov = self%nkrylov_total
      stats%nresidual = self%nresidual_total
      stats%njvp = self%njvp_total
      stats%workspace_wall = self%workspace_wall
      stats%vector_wall = self%vector_wall
   end function jfnk_get_statistics

   !> Construct a JFNK solver with optional inner-solve and line-search controls
   !>
   !> @param[in] restart Krylov dimension per restart cycle
   !> @param[in] max_krylov Maximum inner iterations per Newton update
   !> @param[in] eta0 Initial inexact-Newton forcing term
   !> @param[in] eta_min Lower bound of subsequent forcing terms
   !> @param[in] eta_max Upper bound of subsequent forcing terms
   !> @param[in] fd_eps Finite-difference displacement; zero selects automatic scaling
   !> @param[in] verbosity Diagnostic verbosity
   !> @param[in] restart_stagnation_tol Minimum restart reduction; zero disables early exit
   !> @param[in] cgs2 Use reorthogonalized classical Gram-Schmidt
   !> @param[in] exact_jvp Prefer evaluator-supplied matrix-free products
   !> @param[in] recycle Maximum recycled Krylov dimension; zero disables recycling
   !> @param[in] deflation Enable evaluator-supplied coarse-space deflation
   !> @param[in] deflation_rmse Residual threshold for deflation
   !> @param[in] alpha_max Maximum line-search step fraction
   function new_jfnk_solver(restart, max_krylov, &
                             eta0, eta_min, eta_max, fd_eps, &
                             verbosity, restart_stagnation_tol, &
                             cgs2, exact_jvp, recycle, &
                             deflation, deflation_rmse, alpha_max) &
      result(self)
      !> Krylov dimension per restart cycle
      integer, intent(in), optional :: restart
      !> Maximum inner iterations per Newton update
      integer, intent(in), optional :: max_krylov
      !> Initial inexact-Newton forcing term
      real(wp), intent(in), optional :: eta0
      !> Lower bound of subsequent forcing terms
      real(wp), intent(in), optional :: eta_min
      !> Upper bound of subsequent forcing terms
      real(wp), intent(in), optional :: eta_max
      !> Finite-difference displacement; zero selects automatic scaling
      real(wp), intent(in), optional :: fd_eps
      !> Diagnostic verbosity
      integer, intent(in), optional :: verbosity
      !> Minimum restart reduction; zero disables early exit
      real(wp), intent(in), optional :: restart_stagnation_tol
      !> Use reorthogonalized classical Gram-Schmidt
      logical, intent(in), optional :: cgs2
      !> Use exact Jacobian-vector products where the evaluator provides them
      logical, intent(in), optional :: exact_jvp
      !> Recycled Krylov subspace dimension (GCRO-DR); 0 disables recycling
      integer, intent(in), optional :: recycle
      !> Enable the coarse-space deflation where the evaluator supports it
      logical, intent(in), optional :: deflation
      !> Residual RMS below which the deflation is used
      real(wp), intent(in), optional :: deflation_rmse
      !> Longest line-search step beyond the full Newton step, 1 = off
      real(wp), intent(in), optional :: alpha_max
      !> Constructed solver
      type(moist_model_moz_solver_jfnk_type) :: self

      if (present(restart)) self%restart = restart
      if (present(max_krylov)) self%max_krylov = max_krylov
      if (present(eta0)) self%eta0 = eta0
      if (present(eta_min)) self%eta_min = eta_min
      if (present(eta_max)) self%eta_max = eta_max
      if (present(fd_eps)) self%fd_eps = fd_eps
      if (present(verbosity)) self%verbosity = verbosity
      if (present(restart_stagnation_tol)) self%restart_stagnation_tol = restart_stagnation_tol
      if (present(cgs2)) self%use_cgs2 = cgs2
      if (present(exact_jvp)) self%exact_jvp = exact_jvp
      if (present(recycle)) self%recycle = max(0, recycle)
      if (present(deflation)) self%deflation%enabled = deflation
      if (present(deflation_rmse)) self%deflation%rmse_activate = deflation_rmse
      if (present(alpha_max)) self%alpha_max = max(1.0_wp, alpha_max)
   end function new_jfnk_solver

   !> One Newton-Krylov outer step: evaluate residual, GMRES linear solve, Armijo line search
   !> Called once per Newton iteration, in a loop driven by the caller
   !>
   !> @param[inout] self Solver state and workspace
   !> @param[in] residual_fn Model residual evaluator
   !> @param[in] gamma Current field (npts, npair)
   !> @param[out] gamma_new Line-search update (npts, npair)
   !> @param[out] rmse RMS residual of the input field
   !> @param[out] feasible Whether the update is usable
   !> @param[out] error Evaluation failure
   subroutine jfnk_outer_step(self, residual_fn, gamma, gamma_new, rmse, feasible, error)
      !> Solver state and workspace
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Model residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_fn
      !> Current field (npts, npair)
      real(wp), intent(in) :: gamma(:, :)
      !> Line-search update (npts, npair)
      real(wp), intent(out) :: gamma_new(:, :)
      !> RMS residual of the input field
      real(wp), intent(out) :: rmse
      !> Whether the update is usable
      logical, intent(out) :: feasible
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      integer :: npts, npair, nflat, i, j
      real(wp) :: alpha, armijo_lhs, armijo_rhs, slope
      real(wp) :: phi0, phi_a, abs_slope, denom, alpha_new
      !> Number of trial steps taken by the backtracking line search
      integer :: ls_iter
      !> .true. once a trial step satisfies the Armijo sufficient-decrease test
      logical :: accepted, any_feasible
      !> Merit of the best extended step and feasibility of an extension trial
      real(wp) :: phi_best
      logical :: expand_feasible
      logical :: gmres_feasible
      type(error_type), allocatable :: gmres_error

      feasible = .false.
      call validate_solver_arrays(gamma, gamma_new, error)
      if (allocated(error)) return
      self%gmres_iters_last = 0
      self%jvp_ready = .false.
      call self%ensure_workspace(gamma)
      npts = size(gamma, 1)
      npair = size(gamma, 2)
      nflat = npts*npair

      ! Initialize Newton state on first call (sentinel: eta_current < 0)
      if (self%eta_current < 0.0_wp) then
         self%eta_current = self%eta0
         self%precond_valid_state = .false.
         self%rmse_prev_state = huge(1.0_wp)
         self%newton_iter_count = 0
         self%gmres_iters_total = 0
         self%nrecycled = 0
      end if
      self%newton_iter_count = self%newton_iter_count + 1

      ! Evaluate residual F(gamma)
      call residual_fn%residual(gamma, self%workspace%ws_f_curr, feasible, error)
      self%nresidual_total = self%nresidual_total + 1
      if (allocated(error)) feasible = .false.
      if (.not. feasible) return

      rmse = sqrt(field_square_sum(self%workspace%ws_f_curr, self%team_size())/real(nflat, wp))
      if (.not. ieee_is_finite(rmse)) then
         feasible = .false.
         call fatal_error(error, "Solver residual norm must be finite")
         return
      end if
      if (rmse == 0.0_wp) then
         gamma_new = gamma
         self%prev_step_maxabs = 0.0_wp
         self%prev_alpha = 0.0_wp
         return
      end if

      ! Linearize at gamma when the evaluator provides exact Jacobian-vector
      ! products; every product of this Newton step (Arnoldi, restart
      ! residual, Armijo slope) is taken at gamma
      ! A failed linearization falls back to finite differences for this
      ! step instead of aborting it
      self%jvp_ready = .false.
      if (self%exact_jvp) then
         if (residual_fn%has_jvp()) then
            call residual_fn%linearize(gamma, self%jvp_ready, gmres_error)
            if (allocated(gmres_error)) then
               self%jvp_ready = .false.
               deallocate (gmres_error)
            end if
         end if
      end if

      ! Coarse-space deflation of this step's linear solve, once the
      ! residual is below the activation threshold
      ! Evaluator fills the coarse matrix at this linearization point;
      ! evaluators without support (and a singular coarse matrix) leave it
      ! inactive
      call self%deflation%deactivate()
      if (self%jvp_ready .and. self%deflation%wanted(rmse)) then
         call residual_fn%deflation_setup(self%deflation)
      end if

      ! The operator decides from the previous step's residual whether it is used
      self%precond_valid_state = .false.
      if (associated(self%precond)) self%precond_valid_state = self%precond%ready(self%rmse_prev_state)
      call self%note_status()

      ! Set up linear solve
      !$omp parallel do collapse(2) schedule(static) num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
      do j = 1, npair
         do i = 1, npts
            self%workspace%ws_rhs_flat(i + (j - 1)*npts) = self%workspace%ws_f_curr(i, j)
         end do
      end do
      !$omp end parallel do

      call self%linear_solve(residual_fn, gamma, self%workspace%ws_f_curr, &
                             self%workspace%ws_rhs_flat, self%workspace%ws_delta_flat, &
                             self%eta_current, self%precond_valid_state, &
                             gmres_feasible, gmres_error, &
                             gmres_iters=self%gmres_iters_last)
      self%gmres_iters_total = self%gmres_iters_total + self%gmres_iters_last
      self%nkrylov_total = self%nkrylov_total + self%gmres_iters_last
      if (allocated(gmres_error)) gmres_feasible = .false.
      if (.not. gmres_feasible) then
         feasible = .false.
         if (allocated(gmres_error)) call move_alloc(gmres_error, error)
         return
      end if

      call unpack_vector(self%workspace%ws_delta_flat, self%workspace%ws_delta, self%team_size())

      ! Compute true directional derivative for Armijo condition via one extra JVP
      ! As in pysofft, an unusable product falls back to slope = -1, and a
      ! direction that is not a descent direction is replaced by the residual,
      ! searched with slope -||F||^2 before the Picard fallback below
      call self%apply_jacobian(residual_fn, gamma, self%workspace%ws_f_curr, self%workspace%ws_delta_flat, &
                               self%workspace%ws_w, gmres_feasible, gmres_error)
      if (allocated(gmres_error)) then
         gmres_feasible = .false.
         deallocate (gmres_error)
      end if
      slope = -1.0_wp
      if (gmres_feasible) slope = -blocked_dot(self%workspace%ws_rhs_flat, self%workspace%ws_w, self%team_size())
      if (.not. ieee_is_finite(slope)) slope = -1.0_wp
      if (slope >= 0.0_wp) then
         call copy_vector(self%workspace%ws_rhs_flat, self%workspace%ws_delta_flat, self%team_size())
         call unpack_vector(self%workspace%ws_delta_flat, self%workspace%ws_delta, self%team_size())
         slope = -blocked_dot(self%workspace%ws_rhs_flat, self%workspace%ws_rhs_flat, self%team_size())
      end if

      ! Backtracking line search with safeguarded quadratic interpolation of
      ! the merit function phi(a) = 0.5*||F(gamma - a*delta)||^2
      !
      ! - a rejected step's next trial minimises the quadratic model through
      !   phi(0), phi'(0) = slope and the latest phi(a), clamped to
      !   [0.1, tau]*a (Dennis-Schnabel safeguards)
      ! - the interpolation can shrink faster than fixed tau-halving
      ! - max_line_search_iter caps rejected backtracking trials, independently
      !   of the alpha floor; accepted full steps may be extended below
      phi0 = 0.5_wp*field_square_sum(self%workspace%ws_f_curr, self%team_size())
      abs_slope = abs(slope)
      armijo_rhs = self%line_search_c*abs_slope
      alpha = 1.0_wp
      accepted = .false.
      ls_iter = 0
      ! slope < 0 here unless the residual underflows to zero
      do while (slope < 0.0_wp)
         ls_iter = ls_iter + 1
         call trial_field(gamma, self%workspace%ws_delta, alpha, self%workspace%ws_gamma_trial, self%team_size())
         call residual_fn%residual(self%workspace%ws_gamma_trial, self%workspace%ws_f_trial, feasible, gmres_error)
         self%nresidual_total = self%nresidual_total + 1
         if (allocated(gmres_error)) feasible = .false.
         if (feasible) then
            phi_a = 0.5_wp*field_square_sum(self%workspace%ws_f_trial, self%team_size())
            feasible = ieee_is_finite(phi_a)
         end if
         if (feasible) then
            armijo_lhs = phi0 - phi_a
            if (armijo_lhs >= alpha*armijo_rhs) then
               ! Sufficient decrease achieved: keep this Newton step regardless of
               ! how small alpha is, so a valid short step is not thrown away in
               ! favour of the weaker Picard fallback below
               accepted = .true.
               exit
            end if
            ! Minimiser of the quadratic model through phi(0), phi'(0), phi(a)
            ! Armijo failure implies denom > 0, but guard against roundoff
            denom = 2.0_wp*(phi_a - phi0 + alpha*abs_slope)
            if (denom > 0.0_wp) then
               alpha_new = abs_slope*alpha*alpha/denom
               alpha_new = min(max(alpha_new, 0.1_wp*alpha), self%line_search_tau*alpha)
            else
               alpha_new = self%line_search_tau*alpha
            end if
            alpha = alpha_new
         else
            ! Infeasible trial: no merit value to interpolate; shrink by tau
            alpha = self%line_search_tau*alpha
         end if
         if (alpha < self%line_search_min_alpha) exit
         if (ls_iter >= self%max_line_search_iter) exit
      end do

      ! Extend a full step that still leaves much of the merit: double alpha up
      ! to alpha_max while each longer step is feasible, meets the Armijo test,
      ! and lowers the merit below the best step so far; the first worse trial
      ! ends the extension and the best step is kept
      if (accepted .and. ls_iter == 1 .and. self%alpha_max > 1.0_wp) then
         if (phi_a > self%line_search_expand_ratio*phi0) then
            phi_best = phi_a
            do while (alpha < self%alpha_max)
               alpha_new = min(2.0_wp*alpha, self%alpha_max)
               call trial_field(gamma, self%workspace%ws_delta, alpha_new, self%workspace%ws_gamma_trial, self%team_size())
               call residual_fn%residual(self%workspace%ws_gamma_trial, self%workspace%ws_f_trial, expand_feasible, &
                  & gmres_error)
               self%nresidual_total = self%nresidual_total + 1
               if (allocated(gmres_error)) then
                  expand_feasible = .false.
                  deallocate (gmres_error)
               end if
               if (.not. expand_feasible) exit
               phi_a = 0.5_wp*field_square_sum(self%workspace%ws_f_trial, self%team_size())
               if (.not. ieee_is_finite(phi_a)) exit
               if (phi_a >= phi_best .or. phi0 - phi_a < alpha_new*armijo_rhs) exit
               alpha = alpha_new
               phi_best = phi_a
            end do
         end if
      end if

      ! If the GMRES line search never met the Armijo test (floor reached, trial
      ! budget exhausted, or every trial infeasible), fall back to a
      ! spectrally-filtered Picard step
      if (.not. accepted) then
         call copy_vector(self%workspace%ws_rhs_flat, self%workspace%ws_delta_flat, self%team_size())
         ! The lent fallback filter, not the active preconditioner
         if (associated(self%smoother)) call self%smoother%apply(self%workspace%ws_delta_flat)
         call unpack_vector(self%workspace%ws_delta_flat, self%workspace%ws_delta, self%team_size())
         ! Backtrack until a trial is both feasible and decreasing; take no
         ! step if none is. Only an evaluated, decreasing trial is accepted
         ! `feasible` is restored afterwards because the trial evaluation
         ! overwrites this routine's intent(out) flag
         alpha = self%line_search_short_alpha
         any_feasible = .false.
         do
            call trial_field(gamma, self%workspace%ws_delta, alpha, self%workspace%ws_gamma_trial, self%team_size())
            call residual_fn%residual(self%workspace%ws_gamma_trial, self%workspace%ws_f_trial, feasible, gmres_error)
            self%nresidual_total = self%nresidual_total + 1
            if (allocated(gmres_error)) feasible = .false.
            if (feasible) then
               phi_a = 0.5_wp*field_square_sum(self%workspace%ws_f_trial, self%team_size())
               feasible = ieee_is_finite(phi_a)
            end if
            if (feasible) then
               any_feasible = .true.
               if (phi_a < phi0) then
                  accepted = .true.
                  exit
               end if
            end if
            alpha = alpha*self%line_search_tau
            if (alpha < self%line_search_min_alpha) exit
         end do
         if (.not. accepted) alpha = 0.0_wp
         ! Every trial non-finite: report infeasible (with the last error, if
         ! any) instead of pretending the unchanged iterate is fine
         feasible = any_feasible
         if (.not. any_feasible) then
            if (allocated(gmres_error)) call move_alloc(gmres_error, error)
            if (.not. allocated(error)) then
               call fatal_error(error, "Newton and Picard trial evaluations were infeasible")
            end if
            return
         end if
         ! A non-decreasing search returns a zero step and the next update
         ! retries with a tightened forcing term. At eta_min the first zero
         ! step is still returned, so a driver can accept an input that meets
         ! its tolerance; a second one in a row would repeat the same solve and
         ! is reported as infeasible without an error. The accepted-step
         ! diagnostics are kept as on any failure
         if (.not. accepted .and. self%eta_current <= self%eta_min) then
            if (self%stalled_at_min) then
               feasible = .false.
               return
            end if
            self%stalled_at_min = .true.
         end if
      end if
      if (accepted) self%stalled_at_min = .false.

      self%prev_step_maxabs = alpha*maxval(abs(self%workspace%ws_delta_flat))
      self%prev_alpha = alpha
      call trial_field(gamma, self%workspace%ws_delta, alpha, gamma_new, self%team_size())

      ! Forcing term for the next update: Eisenstat-Walker Choice 2 on this and
      ! the previous input residual, so it lags one update
      ! After a zero step the retry starts from the same iterate; tighten the
      ! tolerance so it gets a more accurate Newton direction
      if (alpha <= 0.0_wp) then
         self%eta_current = max(self%eta_min, 0.1_wp*self%eta_current)
      else if (self%rmse_prev_state < huge(1.0_wp)*0.5_wp) then
         block
            real(wp) :: eta_trial
            eta_trial = self%forcing_beta*(rmse/self%rmse_prev_state)**self%forcing_power
            self%eta_current = max(self%eta_min, min(self%eta_max, eta_trial))
         end block
      end if
      self%rmse_prev_state = rmse
   end subroutine jfnk_outer_step

   !> Inexact Newton linear solve J delta = F by left-preconditioned GMRES
   !>
   !> Restarted GMRES(m) with modified Gram-Schmidt (or CGS2) Arnoldi
   !>
   !> - `recycle` = k > 0 makes it GCRO-DR(m, k) (Parks, de Sturler,
   !>   Mackey, Johnson and Maiti, SIAM J Sci Comput 28, 1651 (2006))
   !> - a recycled pair (U, C = A U, C^T C = I) carries across restarts
   !>   and Newton steps
   !> - each solve rebuilds C for the current Jacobian (k products)
   !> - takes the minimal-residual correction in span U
   !> - runs each cycle's m - k Arnoldi steps on (I - C C^T) A, so
   !>   directions kept in U need not be rediscovered every cycle
   !> - every cycle replaces U by the harmonic Ritz vectors of its
   !>   smallest eigenvalues
   !> - with k = 0 this is exactly plain restarted GMRES
   !>
   !> With an active coarse space (`deflation%active`) the solve is DEF1
   !> deflated (see moist_model_moz_solver_deflation)
   !>
   !> - GMRES runs on the operator
   !>   P (I - J Q) J x = P (I - J Q) b
   !> - every return adds the coarse part, delta = x + Q (b - J x)
   !> - each Arnoldi step then costs two Jacobian products,
   !>   J (v - Q J v), and recycling is skipped for the solve
   !>
   !> @param[inout] self          Solver (workspace, recycled subspace)
   !> @param[in]    residual_cb   Residual evaluator
   !> @param[in]    gamma_curr    Current Newton point
   !> @param[in]    f_curr        Residual at the current Newton point
   !> @param[in]    rhs_flat      Right-hand side F (flat)
   !> @param[out]   delta_flat    Newton direction (flat)
   !> @param[in]    eta           Forcing term (relative tolerance)
   !> @param[inout] precond_valid Whether the left preconditioner is applied
   !> @param[out]   feasible      Whether every product was usable
   !> @param[out]   error         Error handle
   !> @param[out]   gmres_iters   Arnoldi iterations plus recycled-image products
   subroutine gmres_linear_solve(self, residual_cb, gamma_curr, f_curr, rhs_flat, delta_flat, &
                                 eta, precond_valid, feasible, error, gmres_iters)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      real(wp), intent(in) :: gamma_curr(:, :)
      real(wp), intent(in) :: f_curr(:, :)
      real(wp), intent(in) :: rhs_flat(:)
      real(wp), intent(out) :: delta_flat(:)
      real(wp), intent(in) :: eta
      logical, intent(inout) :: precond_valid
      logical, intent(out) :: feasible
      type(error_type), allocatable, intent(out) :: error
      integer, intent(out), optional :: gmres_iters

      integer :: nflat, m, k, j, i, gmres_iter, nrestart, last_iter
      !> Recycled vectors in use, cycle length (Arnoldi steps), maximum recycle dimension,
      !> Jacobian products of the recycle refresh
      integer :: kr, mk, kmax, njv
      !> Cycle outcome
      logical :: converged, last_cycle, happy_breakdown
      !> Start of the current algebra-only timing interval
      integer(int64) :: tick
      real(wp) :: beta, rho, h_jj, tol, rho_prev_restart
      !> B y, the recycled-space part of the cycle correction
      real(wp), allocatable :: by(:)
      !> Whether this solve is deflated (DEF1)
      logical :: defl
      !> |P b| of the undeflated right-hand side
      real(wp) :: beta_rhs
      !> Coarse coordinates E^{-1} Z^T Omega w, (kz, 1)
      real(wp), allocatable :: ycoarse(:, :)

      feasible = .true.
      call system_clock(tick)
      nflat = size(rhs_flat)
      !> Krylov dimension of one cycle; matches `jfnk_ensure_workspace`
      m = max(1, self%restart)
      kmax = self%recycle_dim()
      last_iter = 0
      gmres_iter = 0
      kr = 0
      if (present(gmres_iters)) gmres_iters = 0

      ! A deflated solve skips recycling and drops the recycled subspace
      defl = self%deflation%active
      if (defl) then
         kmax = 0
         self%nrecycled = 0
         allocate (ycoarse(self%deflation%kz, 1))
         if (allocated(self%workspace%ws_z)) then
            if (size(self%workspace%ws_z) /= nflat) deallocate (self%workspace%ws_z)
         end if
         if (.not. allocated(self%workspace%ws_z)) allocate (self%workspace%ws_z(nflat))
      end if

      !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
      do i = 1, nflat
         delta_flat(i) = 0.0_wp
         self%workspace%ws_r(i) = rhs_flat(i)
      end do
      !$omp end parallel do
      if (defl) then
         ! Newton residual of the final direction equals the deflated GMRES
         ! residual, P (b - J delta) = P (I - J Q)(b - J x),
         ! so the inexact-Newton tolerance stays eta*|P b| of the
         ! undeflated right-hand side
         ! The oblique projector can inflate the right-hand side, so keep
         ! the tolerance relative to the undeflated preconditioned residual
         call copy_vector(rhs_flat, self%workspace%ws_w, self%team_size())
         call accumulate_wall(tick, self%vector_wall)
         if (precond_valid) call self%apply_precond(self%workspace%ws_w)
         call system_clock(tick)
         beta_rhs = sqrt(blocked_dot(self%workspace%ws_w, self%workspace%ws_w, self%team_size()))
         if (.not. ieee_is_finite(beta_rhs)) then
            feasible = .false.
            call fatal_error(error, "Preconditioned residual norm must be finite")
            return
         end if
         if (beta_rhs == 0.0_wp .and. any(rhs_flat /= 0.0_wp)) then
            feasible = .false.
            call fatal_error(error, "Preconditioner annihilated a nonzero residual")
            return
         end if
         ! r := (I - J Q) b
         call accumulate_wall(tick, self%vector_wall)
         call self%deflate(residual_cb, gamma_curr, f_curr, self%workspace%ws_r, ycoarse, feasible, error)
         call system_clock(tick)
         if (.not. feasible) then
            call accumulate_wall(tick, self%vector_wall)
            return
         end if
      end if
      if (precond_valid) then
         call accumulate_wall(tick, self%vector_wall)
         call self%apply_precond(self%workspace%ws_r)
         call system_clock(tick)
      end if
      beta = sqrt(blocked_dot(self%workspace%ws_r, self%workspace%ws_r, self%team_size()))
      if (.not. ieee_is_finite(beta)) then
         feasible = .false.
         call fatal_error(error, "Preconditioned residual norm must be finite")
         return
      end if
      if (beta == 0.0_wp) then
         ! An exact coarse-space correction can legitimately remove the RHS.
         ! An invertible preconditioner alone cannot remove a nonzero RHS.
         if (precond_valid .and. .not. defl .and. any(rhs_flat /= 0.0_wp)) then
            feasible = .false.
            call fatal_error(error, "Preconditioner annihilated a nonzero residual")
            return
         end if
         call accumulate_wall(tick, self%vector_wall)
         if (defl) call self%deflation_finish(residual_cb, gamma_curr, f_curr, rhs_flat, &
                                              delta_flat, .true., ycoarse, feasible, error)
         return
      end if

      tol = eta*beta
      if (defl) tol = eta*beta_rhs
      ! Preserve relative accuracy when the left operator is rescaled.
      tol = max(tiny(1.0_wp), tol)

      ! Recycled subspace: rebuild C = A U for this Jacobian, then take
      ! the minimal-residual correction in span U,
      ! delta = U C^T r, r := r - C C^T r
      ! A subspace built for the other preconditioning state is dropped
      if (kmax > 0 .and. self%nrecycled > 0) then
         if (precond_valid .eqv. self%recycle_precond) then
            call accumulate_wall(tick, self%vector_wall)
            call self%recycle_refresh(residual_cb, gamma_curr, f_curr, precond_valid, njv, &
                                      feasible, error)
            call system_clock(tick)
            gmres_iter = njv
            last_iter = njv
            if (.not. feasible) then
               if (present(gmres_iters)) gmres_iters = last_iter
               return
            end if
            kr = self%nrecycled
            if (kr > 0) then
               call project_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_r, &
                  & self%workspace%ws_proj(1:kr), self%team_size())
               call update_basis(nflat, kr, self%workspace%ws_u(:, 1:kr), self%workspace%ws_proj(1:kr), 1.0_wp, &
                                 delta_flat, self%team_size())
               call update_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_proj(1:kr), -1.0_wp, &
                                 self%workspace%ws_r, self%team_size())
               beta = sqrt(blocked_dot(self%workspace%ws_r, self%workspace%ws_r, self%team_size()))
            end if
         else
            self%nrecycled = 0
         end if
      end if
      if (kmax > 0) allocate (by(kmax))

      if (.not. ieee_is_finite(beta)) then
         feasible = .false.
         call fatal_error(error, "Recycled residual norm must be finite")
         return
      end if
      if (beta < tol) then
         if (present(gmres_iters)) gmres_iters = last_iter
         call accumulate_wall(tick, self%vector_wall)
         if (defl) call self%deflation_finish(residual_cb, gamma_curr, f_curr, rhs_flat, &
                                              delta_flat, .true., ycoarse, feasible, error)
         return
      end if

      call divide_vector(self%workspace%ws_r, beta, self%workspace%ws_v(:, 1), self%team_size())
      self%workspace%ws_g = 0.0_wp
      self%workspace%ws_g(1) = beta
      self%workspace%ws_h = 0.0_wp
      if (kmax > 0) then
         self%workspace%ws_hbar = 0.0_wp
         self%workspace%ws_bproj = 0.0_wp
      end if
      nrestart = 0
      rho_prev_restart = beta
      !> Arnoldi steps per cycle: the recycled vectors take part of the m budget
      mk = max(1, m - kr)
      k = 0

      do while (gmres_iter < self%max_krylov)
         gmres_iter = gmres_iter + 1
         last_iter = gmres_iter
         k = k + 1

         ! Arnoldi step: compute w = M^{-1} J v_k using preallocated workspace
         call accumulate_wall(tick, self%vector_wall)
         call self%apply_jacobian(residual_cb, gamma_curr, f_curr, self%workspace%ws_v(:, k), self%workspace%ws_w, &
                                  feasible, error)
         if (.not. feasible) then
            if (present(gmres_iters)) gmres_iters = last_iter
            return
         end if

         ! DEF1: w := (I - J Q) J v,
         ! computed as J (v - Q J v), one more product
         if (defl) then
            call self%deflation%coordinates(self%workspace%ws_w, ycoarse, nthreads=self%team_size())
            call copy_vector(self%workspace%ws_v(:, k), self%workspace%ws_z, self%team_size())
            call self%deflation%expand(ycoarse, -1.0_wp, self%workspace%ws_z, nthreads=self%team_size())
            call self%apply_jacobian(residual_cb, gamma_curr, f_curr, self%workspace%ws_z, self%workspace%ws_w, &
                                     feasible, error)
            if (.not. feasible) then
               if (present(gmres_iters)) gmres_iters = last_iter
               return
            end if
         end if

         ! Apply left preconditioner
         if (precond_valid) call self%apply_precond(self%workspace%ws_w)
         call system_clock(tick)

         ! GCRO deflation: w := (I - C C^T) w,
         ! keeping B(:, k) = C^T A v_k
         if (kr > 0) then
            call project_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_w, &
               & self%workspace%ws_bproj(1:kr, k), self%team_size())
            call update_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_bproj(1:kr, k), -1.0_wp, &
                              self%workspace%ws_w, self%team_size())
         end if

         ! Arnoldi orthogonalization of w against v_1..v_k
         if (self%use_cgs2) then
            ! Classical Gram-Schmidt with one reorthogonalization pass (CGS2)
            !
            ! - each pass streams the large w vector only twice: a BLAS-2
            !   projection s = V^T w and a back-update w := w - V s, instead
            !   of the 2*k separate sweeps modified Gram-Schmidt makes
         ! - the second pass reduces loss of orthogonality
         ! - fewer sweeps do not guarantee a faster solve; rounding can also
         !   change the Newton path. CGS2 is opt-in and MGS is the default
            call project_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_w, &
               & self%workspace%ws_proj(1:k), self%team_size())
            call update_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_proj(1:k), -1.0_wp, &
               & self%workspace%ws_w, self%team_size())
            call project_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_w, &
               & self%workspace%ws_proj2(1:k), self%team_size())
            call update_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_proj2(1:k), -1.0_wp, &
               & self%workspace%ws_w, self%team_size())
            self%workspace%ws_h(1:k, k) = self%workspace%ws_proj(1:k) + self%workspace%ws_proj2(1:k)
         else
            call orthogonalize_mgs(self%workspace%ws_v(:, 1:k), self%workspace%ws_w, self%workspace%ws_h(1:k, k), &
               & nthreads=self%team_size())
         end if

         h_jj = sqrt(blocked_dot(self%workspace%ws_w, self%workspace%ws_w, self%team_size()))
         if (.not. ieee_is_finite(h_jj) .or. .not. all(ieee_is_finite(self%workspace%ws_h(1:k, k)))) then
            feasible = .false.
            if (present(gmres_iters)) gmres_iters = last_iter
            call fatal_error(error, "GMRES Jacobian product and preconditioner must be finite")
            return
         end if
         self%workspace%ws_h(k + 1, k) = h_jj

         ! Compare against this operator column, rather than an absolute
         ! cutoff that changes meaning when the preconditioner is rescaled.
         happy_breakdown = h_jj <= 10.0_wp*epsilon(1.0_wp)*maxval(abs(self%workspace%ws_h(1:k, k)))
         if (.not. happy_breakdown) then
            call divide_vector(self%workspace%ws_w, h_jj, self%workspace%ws_v(:, k + 1), self%team_size())
         else
            ! Happy breakdown: the Krylov space is invariant
            ! Continuing would otherwise reuse whatever stale content
            ! ws_v(:,k+1) holds as an orthonormal basis vector
            !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
            do i = 1, nflat
               self%workspace%ws_v(i, k + 1) = 0.0_wp
            end do
            !$omp end parallel do
         end if
         ! Unrotated Hessenberg column for the harmonic Ritz harvest
         if (kmax > 0) self%workspace%ws_hbar(1:k + 1, k) = self%workspace%ws_h(1:k + 1, k)

         do j = 1, k - 1
            call apply_givens(self%workspace%ws_cs(j), self%workspace%ws_sn(j), self%workspace%ws_h(j, k), &
               & self%workspace%ws_h(j + 1, k))
         end do
         call generate_givens(self%workspace%ws_h(k, k), self%workspace%ws_h(k + 1, k), &
            & self%workspace%ws_cs(k), self%workspace%ws_sn(k))
         call apply_givens(self%workspace%ws_cs(k), self%workspace%ws_sn(k), self%workspace%ws_h(k, k), &
            & self%workspace%ws_h(k + 1, k))
         call apply_givens(self%workspace%ws_cs(k), self%workspace%ws_sn(k), self%workspace%ws_g(k), self%workspace%ws_g(k + 1))

         rho = abs(self%workspace%ws_g(k + 1))

         if (self%verbosity >= 5 .and. associated(self%iter_printer)) then
            if (gmres_iter <= 4 .or. mod(gmres_iter, 5) == 0 .or. rho < tol) then
               call self%iter_printer%print_inner_gmres_row(gmres_iter, k, rho, tol)
            end if
         end if

         ! Least-squares solution is only formed here; reaching the iteration
         ! cap also ends the cycle, so a capped solve (max_krylov < restart)
         ! still returns its Newton direction
         if (rho < tol .or. happy_breakdown .or. k == mk .or. gmres_iter == self%max_krylov) then
            call upper_backsolve(self%workspace%ws_h(1:k, 1:k), self%workspace%ws_g(1:k), self%workspace%ws_y(1:k), feasible)
            if (.not. feasible) then
               if (present(gmres_iters)) gmres_iters = last_iter
               call fatal_error(error, "GMRES encountered a singular or non-finite Krylov system")
               return
            end if
            call update_basis(nflat, k, self%workspace%ws_v(:, 1:k), self%workspace%ws_y(1:k), 1.0_wp, delta_flat, self%team_size())
            ! GCRO: the recycled coordinates of the correction are -B y
            if (kr > 0) then
               by(1:kr) = matmul(self%workspace%ws_bproj(1:kr, 1:k), self%workspace%ws_y(1:k))
               call update_basis(nflat, kr, self%workspace%ws_u(:, 1:kr), by(1:kr), -1.0_wp, delta_flat, self%team_size())
            end if
            converged = rho < tol
            last_cycle = converged .or. happy_breakdown .or. gmres_iter == self%max_krylov

            ! Harvest the recycled subspace from this cycle: always once one
            ! exists, otherwise from the first cycle that filled the restart
            ! length or the iteration budget, so a solve capped below the
            ! restart length (max_krylov < restart) still starts recycling
            ! C is only needed for another cycle of this solve; the next
            ! Newton step rebuilds it for its own Jacobian
            if (kmax > 0 .and. (kr > 0 .or. k == mk .or. gmres_iter == self%max_krylov) &
               & .and. kr + k > kmax) then
               call self%recycle_harvest(kr, k, .not. last_cycle)
               self%recycle_precond = precond_valid
            end if

            if (last_cycle) exit

            ! Check for restart stagnation: bail out if this cycle barely reduced residual
            if (rho_prev_restart > 0.0_wp .and. self%restart_stagnation_tol > 0.0_wp) then
               if ((rho_prev_restart - rho)/rho_prev_restart < self%restart_stagnation_tol) then
                  exit
               end if
            end if
            rho_prev_restart = rho

            ! Restart: recompute preconditioned residual r = M^{-1}(b - J*delta)
            nrestart = nrestart + 1
            call copy_vector(rhs_flat, self%workspace%ws_r, self%team_size())
            call accumulate_wall(tick, self%vector_wall)
            call self%apply_jacobian(residual_cb, gamma_curr, f_curr, delta_flat, self%workspace%ws_w, &
                                     feasible, error)
            if (.not. feasible) then
               if (present(gmres_iters)) gmres_iters = last_iter
               return
            end if
            call system_clock(tick)
            !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
            do i = 1, nflat
               self%workspace%ws_r(i) = self%workspace%ws_r(i) - self%workspace%ws_w(i)
            end do
            !$omp end parallel do
            if (defl) then
               call accumulate_wall(tick, self%vector_wall)
               call self%deflate(residual_cb, gamma_curr, f_curr, self%workspace%ws_r, ycoarse, feasible, &
                                 error)
               call system_clock(tick)
               if (.not. feasible) then
                  if (present(gmres_iters)) gmres_iters = last_iter
                  return
               end if
            end if
            ! Apply left preconditioner to restarted residual
            if (precond_valid) then
               call accumulate_wall(tick, self%vector_wall)
               call self%apply_precond(self%workspace%ws_r)
               call system_clock(tick)
            end if
            ! Deflate the restarted residual against the harvested C
            kr = self%nrecycled
            if (kr > 0) then
               call project_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_r, &
                  & self%workspace%ws_proj(1:kr), self%team_size())
               call update_basis(nflat, kr, self%workspace%ws_u(:, 1:kr), self%workspace%ws_proj(1:kr), 1.0_wp, &
                                 delta_flat, self%team_size())
               call update_basis(nflat, kr, self%workspace%ws_c(:, 1:kr), self%workspace%ws_proj(1:kr), -1.0_wp, &
                                 self%workspace%ws_r, self%team_size())
            end if
            beta = sqrt(blocked_dot(self%workspace%ws_r, self%workspace%ws_r, self%team_size()))
            if (.not. ieee_is_finite(beta)) then
               feasible = .false.
               if (present(gmres_iters)) gmres_iters = last_iter
               call fatal_error(error, "Restarted residual norm must be finite")
               return
            end if
            if (beta < tol) exit

            call divide_vector(self%workspace%ws_r, beta, self%workspace%ws_v(:, 1), self%team_size())
            self%workspace%ws_g = 0.0_wp
            self%workspace%ws_g(1) = beta
            self%workspace%ws_h = 0.0_wp
            if (kmax > 0) then
               self%workspace%ws_hbar = 0.0_wp
               self%workspace%ws_bproj = 0.0_wp
            end if
            mk = max(1, m - kr)
            k = 0
         end if
      end do

      if (present(gmres_iters)) gmres_iters = last_iter
      call accumulate_wall(tick, self%vector_wall)
      if (defl) call self%deflation_finish(residual_cb, gamma_curr, f_curr, rhs_flat, delta_flat, &
                                           .false., ycoarse, feasible, error)

   end subroutine gmres_linear_solve

   !> Deflate a residual in place, vec := (I - J Q) vec,
   !> Q = Z E^{-1} Z^T Omega
   !> Costs one Jacobian product; uses ws_z and ws_w, so vec must be neither
   !>
   !> @param[inout] self        Solver with an active deflation
   !> @param[in]    residual_cb Residual evaluator
   !> @param[in]    gamma_curr  Current Newton point
   !> @param[in]    f_curr      Residual at the current Newton point
   !> @param[inout] vec         Flat vector to deflate
   !> @param[out]   y           Coarse coordinates E^{-1} Z^T Omega vec, (kz, 1)
   !> @param[out]   feasible    Whether the product was usable
   !> @param[out]   error       Error handle
   subroutine gmres_deflate(self, residual_cb, gamma_curr, f_curr, vec, y, feasible, error)
      !> Solver
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      !> Current Newton point
      real(wp), intent(in) :: gamma_curr(:, :)
      !> Residual at the current Newton point
      real(wp), intent(in) :: f_curr(:, :)
      !> Flat vector, deflated in place
      real(wp), intent(inout) :: vec(:)
      !> Coarse coordinates
      real(wp), intent(out) :: y(:, :)
      !> Whether the product was usable
      logical, intent(out) :: feasible
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Vector index
      integer :: i

      call self%deflation%coordinates(vec, y, nthreads=self%team_size())
      !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, size(vec) >= parallel_vector_size))
      do i = 1, size(vec)
         self%workspace%ws_z(i) = 0.0_wp
      end do
      !$omp end parallel do
      call self%deflation%expand(y, 1.0_wp, self%workspace%ws_z, nthreads=self%team_size())
      call self%apply_jacobian(residual_cb, gamma_curr, f_curr, self%workspace%ws_z, self%workspace%ws_w, feasible, &
                               error)
      if (.not. feasible) return
      !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, size(vec) >= parallel_vector_size))
      do i = 1, size(vec)
         vec(i) = vec(i) - self%workspace%ws_w(i)
      end do
      !$omp end parallel do
   end subroutine gmres_deflate

   !> DEF1 coarse correction of a deflated solve:
   !> delta := x + Q (b - J x),
   !> which is Q b + (I - Q J) x
   !> Costs one Jacobian product unless x = 0
   !>
   !> @param[inout] self        Solver with an active deflation
   !> @param[in]    residual_cb Residual evaluator
   !> @param[in]    gamma_curr  Current Newton point
   !> @param[in]    f_curr      Residual at the current Newton point
   !> @param[in]    rhs_flat    Right-hand side b (flat)
   !> @param[inout] delta_flat  GMRES iterate x in, Newton direction out
   !> @param[in]    x_zero      Whether x = 0 (no Arnoldi step taken)
   !> @param[out]   y           Coarse coordinate scratch, (kz, 1)
   !> @param[out]   feasible    Whether the product was usable
   !> @param[out]   error       Error handle
   subroutine gmres_deflation_finish(self, residual_cb, gamma_curr, f_curr, rhs_flat, delta_flat, &
                                     x_zero, y, feasible, error)
      !> Solver
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      !> Current Newton point
      real(wp), intent(in) :: gamma_curr(:, :)
      !> Residual at the current Newton point
      real(wp), intent(in) :: f_curr(:, :)
      !> Right-hand side b
      real(wp), intent(in) :: rhs_flat(:)
      !> GMRES iterate x in, Newton direction out
      real(wp), intent(inout) :: delta_flat(:)
      !> Whether x = 0
      logical, intent(in) :: x_zero
      !> Coarse coordinate scratch
      real(wp), intent(out) :: y(:, :)
      !> Whether the product was usable
      logical, intent(out) :: feasible
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Vector index
      integer :: i

      feasible = .true.
      if (x_zero) then
         call copy_vector(rhs_flat, self%workspace%ws_w, self%team_size())
      else
         call self%apply_jacobian(residual_cb, gamma_curr, f_curr, delta_flat, self%workspace%ws_w, &
                                  feasible, error)
         if (.not. feasible) return
         !$omp parallel do schedule(static) num_threads(merge(self%team_size(), 1, size(rhs_flat) >= parallel_vector_size))
         do i = 1, size(rhs_flat)
            self%workspace%ws_w(i) = rhs_flat(i) - self%workspace%ws_w(i)
         end do
         !$omp end parallel do
      end if
      call self%deflation%coordinates(self%workspace%ws_w, y, nthreads=self%team_size())
      call self%deflation%expand(y, 1.0_wp, delta_flat, nthreads=self%team_size())
   end subroutine gmres_deflation_finish

   !> Recycled subspace dimension in use: `recycle`, capped so that every
   !> cycle and every capped solve keep at least one Arnoldi step
   !>
   !> @param[in] self Solver
   !> @return kmax Recycled vectors (0 = plain restarted GMRES)
   pure function gmres_recycle_dim(self) result(kmax)
      !> Solver
      class(moist_model_moz_solver_jfnk_type), intent(in) :: self
      !> Recycled vectors in use
      integer :: kmax
      kmax = max(0, min(self%recycle, self%restart - 1, self%max_krylov - 1))
   end function gmres_recycle_dim

   !> Rebuild the recycled image C = A U for the current Jacobian (GCRO-DR)
   !>
   !> - each recycled vector costs one Jacobian-vector product (and one
   !>   preconditioner application)
   !> - C is then orthonormalized by two passes of Cholesky QR,
   !>   C := C R^{-1}, with the same R^{-1} applied to U
   !>   so that A U = C still holds
   !> - a C too ill-conditioned for Cholesky QR (breakdown, or a diagonal
   !>   ratio of R below sqrt(eps)) drops the recycled subspace for this
   !>   Newton step; the solve then proceeds as plain GMRES
   !>
   !> @param[inout] self          Solver holding U (in) and C (out)
   !> @param[in]    residual_cb   Residual evaluator
   !> @param[in]    gamma_curr    Current Newton point
   !> @param[in]    f_curr        Residual at the current Newton point
   !> @param[in]    precond_valid Whether the left preconditioner is applied
   !> @param[out]   njv           Jacobian-vector products spent
   !> @param[out]   feasible      Whether every product was usable
   !> @param[out]   error         Error handle
   subroutine gmres_recycle_refresh(self, residual_cb, gamma_curr, f_curr, precond_valid, njv, &
                                    feasible, error)
      !> Solver
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      !> Current Newton point and its residual
      real(wp), intent(in) :: gamma_curr(:, :), f_curr(:, :)
      !> Whether the left preconditioner is applied
      logical, intent(in) :: precond_valid
      !> Jacobian-vector products spent
      integer, intent(out) :: njv
      !> Whether every product was usable
      logical, intent(out) :: feasible
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Recycled count, column and pass indices, Cholesky status
      integer :: kr, i, pass, info
      !> Gram matrix / Cholesky factor and its inverse
      real(wp), allocatable :: g(:, :), rinv(:, :)
      !> Smallest and largest diagonal entry of the Cholesky factor
      real(wp) :: dmin, dmax

      feasible = .true.
      njv = 0
      kr = self%nrecycled
      if (kr == 0) return
      do i = 1, kr
         call self%apply_jacobian(residual_cb, gamma_curr, f_curr, self%workspace%ws_u(:, i), &
                                  self%workspace%ws_c(:, i), feasible, error)
         njv = njv + 1
         if (.not. feasible) then
            self%nrecycled = 0
            return
         end if
         if (precond_valid) then
            call self%apply_precond(self%workspace%ws_c(:, i))
         end if
      end do

      allocate (g(kr, kr), rinv(kr, kr))
      do pass = 1, 2
         call gram_basis(self%workspace%ws_c, kr, self%workspace%ws_c, kr, g, self%team_size())
         call wrap_potrf(g, info, "U")
         if (info == 0) then
            do i = 1, kr
               g(i + 1:kr, i) = 0.0_wp
            end do
            dmin = minval([(g(i, i), i=1, kr)])
            dmax = maxval([(g(i, i), i=1, kr)])
            if (dmin < sqrt(epsilon(1.0_wp))*dmax) info = 1
         end if
         if (info /= 0) then
            self%nrecycled = 0
            return
         end if
         call upper_inverse(g, rinv)
         call combine_basis(self%workspace%ws_c, kr, rinv, kr, self%team_size())
         call combine_basis(self%workspace%ws_u, kr, rinv, kr, self%team_size())
      end do
   end subroutine gmres_recycle_refresh

   !> Harvest the recycled subspace from a finished cycle (GCRO-DR)
   !>
   !> With U~ = U D (D scales the recycled columns to unit norm) the
   !> cycle satisfies A [U~ V_j] = [C V_{j+1}] G with
   !>
   !>    G = [D B; 0 Hbar], B = C^T A V_j, Hbar the unrotated Hessenberg
   !>
   !> Harmonic Ritz pairs of A on span[U~ V_j] solve
   !> G^T G z = theta G^T W^T Vh z with W = [C V_{j+1}],
   !> Vh = [U~ V_j]
   !>
   !> - k of smallest |theta| approximate the eigenvalues nearest zero
   !>   (on the many-site VV, the deep-core modes that stall restarted
   !>   GMRES)
   !> - with G P = Q R the new pair is
   !>   U := Vh P R^{-1}, C := W Q, which keeps A U = C
   !>   and C^T C = I without a Jacobian-vector product
   !> - a complex conjugate pair contributes its real and imaginary
   !>   parts together, or not at all
   !> - on failure the previous pair is kept; it still satisfies A U = C
   !>
   !> @param[inout] self   Solver holding U, C, V, Hbar and B of the cycle
   !> @param[in]    kr     Recycled vectors used in the cycle
   !> @param[in]    j      Arnoldi steps of the cycle
   !> @param[in]    need_c Whether C is needed (another cycle follows)
   subroutine gmres_recycle_harvest(self, kr, j, need_c)
      !> Solver
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Recycled vectors used in the cycle and its Arnoldi steps
      integer, intent(in) :: kr, j
      !> Whether C must be updated too
      logical, intent(in) :: need_c

      !> Search-space and image dimensions, targets, counters, workspace size
      integer :: p, q, kmax, nsel, nkeep, i, ii, js, lwork
      !> LAPACK status
      integer(lapack_ik) :: lapack_info
      !> Eigenvalue ordering, kept QR columns
      integer, allocatable :: order(:), keep(:)
      logical, allocatable :: taken(:)
      !> Column scaling D, harmonic Ritz magnitudes, generalized eigenvalues
      real(wp), allocatable :: dsc(:), theta(:), alphar(:), alphai(:), betav(:), work(:)
      !> Small matrices of the harvest
      real(wp), allocatable :: gbar(:, :), wtv(:, :), aeig(:, :), beig(:, :), vr(:, :), vl(:, :)
      real(wp), allocatable :: psel(:, :), gp(:, :), qmat(:, :), rmat(:, :), rinv(:, :), pr(:, :)
      real(wp), allocatable :: t1(:, :), t2(:, :), pu(:, :)
      real(wp) :: tmp

      kmax = self%recycle_dim()
      p = kr + j
      q = p + 1
      if (p <= kmax .or. kmax == 0) return

      allocate (dsc(kr), gbar(q, p), wtv(q, p), source=0.0_wp)
      do i = 1, kr
         dsc(i) = 1.0_wp/sqrt(blocked_dot(self%workspace%ws_u(:, i), self%workspace%ws_u(:, i), self%team_size()))
         gbar(i, i) = dsc(i)
      end do
      if (kr > 0) gbar(1:kr, kr + 1:p) = self%workspace%ws_bproj(1:kr, 1:j)
      gbar(kr + 1:q, kr + 1:p) = self%workspace%ws_hbar(1:j + 1, 1:j)
      if (kr > 0) then
         allocate (t1(kr, kr), t2(j + 1, kr))
         call gram_basis(self%workspace%ws_c, kr, self%workspace%ws_u, kr, t1, self%team_size())
         call gram_basis(self%workspace%ws_v, j + 1, self%workspace%ws_u, kr, t2, self%team_size())
         do i = 1, kr
            wtv(1:kr, i) = t1(:, i)*dsc(i)
            wtv(kr + 1:q, i) = t2(:, i)*dsc(i)
         end do
      end if
      do i = 1, j
         wtv(kr + i, kr + i) = 1.0_wp
      end do

      ! Harmonic Ritz pairs: G^T G z = theta G^T (W^T Vh) z
      aeig = matmul(transpose(gbar), gbar)
      beig = matmul(transpose(gbar), wtv)
      lwork = 16*p + 16
      allocate (alphar(p), alphai(p), betav(p), vr(p, p), vl(1, 1), work(lwork), source=0.0_wp)
      call dggev("N", "V", int(p, lapack_ik), aeig, int(p, lapack_ik), beig, int(p, lapack_ik), &
                 alphar, alphai, betav, vl, 1_lapack_ik, vr, int(p, lapack_ik), work, &
                 int(lwork, lapack_ik), lapack_info)
      if (lapack_info /= 0) return

      ! Order by |theta|; infinite eigenvalues (beta = 0) last
      allocate (theta(p), order(p), taken(p))
      do i = 1, p
         tmp = hypot(alphar(i), alphai(i))
         if (abs(betav(i)) > tiny(1.0_wp)*max(1.0_wp, tmp)) then
            theta(i) = tmp/abs(betav(i))
         else
            theta(i) = huge(1.0_wp)
         end if
         order(i) = i
      end do
      do i = 2, p
         ii = order(i)
         js = i - 1
         do while (js >= 1)
            if (theta(order(js)) <= theta(ii)) exit
            order(js + 1) = order(js)
            js = js - 1
         end do
         order(js + 1) = ii
      end do

      ! Select up to k vectors; a conjugate pair enters as its real and imaginary parts
      allocate (psel(p, kmax), source=0.0_wp)
      taken = .false.
      nsel = 0
      do ii = 1, p
         i = order(ii)
         if (taken(i)) cycle
         if (theta(i) >= huge(1.0_wp)) exit
         if (alphai(i) == 0.0_wp) then
            if (nsel + 1 > kmax) exit
            nsel = nsel + 1
            psel(:, nsel) = vr(:, i)
            taken(i) = .true.
         else
            js = i
            if (alphai(i) < 0.0_wp) js = i - 1
            if (nsel + 2 > kmax) exit
            psel(:, nsel + 1) = vr(:, js)
            psel(:, nsel + 2) = vr(:, js + 1)
            nsel = nsel + 2
            taken(js) = .true.
            taken(js + 1) = .true.
         end if
      end do
      if (nsel == 0) return

      ! G P = Q R, dropping directions that G P does not resolve
      gp = matmul(gbar, psel(:, 1:nsel))
      allocate (qmat(q, nsel), rmat(nsel, nsel), keep(nsel))
      call small_qr(gp, qmat, rmat, keep, nkeep)
      if (nkeep == 0) return
      allocate (rinv(nkeep, nkeep))
      call upper_inverse(rmat(1:nkeep, 1:nkeep), rinv)
      pr = matmul(psel(:, keep(1:nkeep)), rinv)

      ! U := U~ P_u R^{-1} + V_j P_v R^{-1},
      ! C := C Q_c + V_{j+1} Q_v (in place)
      allocate (pu(kr, nkeep))
      do i = 1, kr
         pu(i, :) = dsc(i)*pr(i, :)
      end do
      call combine_basis(self%workspace%ws_u, kr, pu, nkeep, self%team_size(), self%workspace%ws_v, j, pr(kr + 1:p, :))
      if (need_c) then
         call combine_basis(self%workspace%ws_c, kr, qmat(1:kr, 1:nkeep), nkeep, self%team_size(), self%workspace%ws_v, j + 1, &
                            qmat(kr + 1:q, 1:nkeep))
      end if
      self%nrecycled = nkeep
   end subroutine gmres_recycle_harvest

   !> Gram block G = X(:, 1:kx)^T Y(:, 1:ky) of two tall bases, reduced over
   !> fixed row blocks in a fixed order (independent of the thread count)
   !>
   !> @param[in]  x,y   Tall bases with the same row count
   !> @param[in]  kx,ky Columns used
   !> @param[out] g     Gram block (kx, ky)
   !> @param[in] nt OpenMP team size
   subroutine gram_basis(x, kx, y, ky, g, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Tall bases
      real(wp), intent(in), contiguous :: x(:, :), y(:, :)
      !> Columns used
      integer, intent(in) :: kx, ky
      !> Gram block
      real(wp), intent(out) :: g(:, :)
      !> Per-block partial products
      real(wp), allocatable :: partial(:, :, :)
      !> Row count, block count, block index and bounds
      integer :: n, nb, b, first, count

      g = 0.0_wp
      n = size(x, 1)
      if (kx == 0 .or. ky == 0 .or. n == 0) return
      nb = (n - 1)/reduction_block_size + 1
      allocate (partial(kx, ky, nb))
      !$omp parallel do schedule(static) private(first, count) num_threads(merge(nt, 1, n >= parallel_vector_size))
      do b = 1, nb
         first = (b - 1)*reduction_block_size + 1
         count = min(reduction_block_size, n - first + 1)
         call gemm_rows_tn(x(:, 1:kx), y(:, 1:ky), first, count, partial(:, :, b), &
            & 1.0_wp, 0.0_wp)
      end do
      !$omp end parallel do
      do b = 1, nb
         g(1:kx, 1:ky) = g(1:kx, 1:ky) + partial(:, :, b)
      end do
   end subroutine gram_basis

   !> In-place basis combination X(:, 1:kout) := X(:, 1:kx) mx [+ Y(:, 1:ky) my],
   !> row block by row block through a small per-thread buffer
   !>
   !> @param[inout] x     Tall basis, overwritten in its first kout columns
   !> @param[in]    kx    Columns of x used as input (may be 0)
   !> @param[in]    mx    Coefficients (kx, kout)
   !> @param[in]    kout  Output columns
   !> @param[in]    y     Optional second tall basis
   !> @param[in]    ky    Columns of y used
   !> @param[in]    my    Coefficients (ky, kout)
   !> @param[in] nt OpenMP team size
   subroutine combine_basis(x, kx, mx, kout, nt, y, ky, my)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Tall basis
      real(wp), intent(inout), contiguous :: x(:, :)
      !> Input and output column counts
      integer, intent(in) :: kx, kout
      !> Coefficients of x
      real(wp), intent(in) :: mx(:, :)
      !> Second basis
      real(wp), intent(in), contiguous, optional :: y(:, :)
      !> Columns of y used
      integer, intent(in), optional :: ky
      !> Coefficients of y
      real(wp), intent(in), optional :: my(:, :)
      !> Rows per block of the buffer
      integer, parameter :: blk = 2048
      !> Per-thread result block
      real(wp), allocatable :: buf(:, :)
      !> Contiguous copies of the coefficients
      real(wp), allocatable :: cx(:, :), cy(:, :)
      !> Row count, columns of y, block index and bounds
      integer :: n, nyc, b, first, count
      !> BLAS accumulation factor: retain the x term if present, otherwise overwrite
      real(wp) :: beta_y

      n = size(x, 1)
      if (kout == 0 .or. n == 0) return
      nyc = 0
      if (present(y) .and. present(ky) .and. present(my)) nyc = ky
      cx = mx(1:kx, 1:kout)
      if (nyc > 0) cy = my(1:nyc, 1:kout)
      beta_y = merge(1.0_wp, 0.0_wp, kx > 0)
      !$omp parallel private(buf, b, first, count) num_threads(merge(nt, 1, n >= parallel_vector_size))
      allocate (buf(blk, kout))
      !$omp do schedule(static)
      do b = 1, (n - 1)/blk + 1
         first = (b - 1)*blk + 1
         count = min(blk, n - first + 1)
         if (kx > 0) call gemm_rows_nn(x(:, 1:kx), first, count, cx, buf, 1.0_wp, 0.0_wp)
         if (nyc > 0) call gemm_rows_nn(y(:, 1:nyc), first, count, cy, buf, 1.0_wp, beta_y)
         x(first:first + count - 1, 1:kout) = buf(1:count, 1:kout)
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel
   end subroutine combine_basis

   !> Thin QR of a small matrix by modified Gram-Schmidt with one
   !> reorthogonalization pass
   !> A column whose remaining norm is below 1e-10 of the largest column
   !> norm is dropped, so a(:, keep(1:nkeep)) = Q R
   !>
   !> @param[in]  a     Matrix (nrow, ncol)
   !> @param[out] q     Orthonormal columns (nrow, ncol), first nkeep valid
   !> @param[out] r     Upper triangular factor (ncol, ncol), leading nkeep block valid
   !> @param[out] keep  Indices of the kept columns
   !> @param[out] nkeep Number of kept columns
   pure subroutine small_qr(a, q, r, keep, nkeep)
      !> Matrix to factor
      real(wp), intent(in) :: a(:, :)
      !> Orthonormal factor
      real(wp), intent(out) :: q(:, :)
      !> Triangular factor
      real(wp), intent(out) :: r(:, :)
      !> Kept column indices
      integer, intent(out) :: keep(:)
      !> Number of kept columns
      integer, intent(out) :: nkeep
      !> Working column
      real(wp) :: w(size(a, 1))
      !> Largest column norm, projection, remaining norm
      real(wp) :: amax, s, nrm
      !> Column, pass and basis indices
      integer :: c, pass, l

      q = 0.0_wp
      r = 0.0_wp
      keep = 0
      nkeep = 0
      amax = 0.0_wp
      do c = 1, size(a, 2)
         amax = max(amax, norm2(a(:, c)))
      end do
      if (amax <= 0.0_wp) return
      do c = 1, size(a, 2)
         w = a(:, c)
         do pass = 1, 2
            do l = 1, nkeep
               s = dot_product(q(:, l), w)
               r(l, nkeep + 1) = r(l, nkeep + 1) + s
               w = w - s*q(:, l)
            end do
         end do
         nrm = norm2(w)
         if (nrm <= 1.0e-10_wp*amax) then
            r(:, nkeep + 1) = 0.0_wp
            cycle
         end if
         nkeep = nkeep + 1
         r(nkeep, nkeep) = nrm
         q(:, nkeep) = w/nrm
         keep(nkeep) = c
      end do
   end subroutine small_qr

   !> Inverse of a nonsingular upper triangular matrix by back substitution
   !>
   !> @param[in]  r    Upper triangular matrix (n, n)
   !> @param[out] rinv Its inverse (n, n), upper triangular
   pure subroutine upper_inverse(r, rinv)
      !> Upper triangular matrix
      real(wp), intent(in) :: r(:, :)
      !> Inverse
      real(wp), intent(out) :: rinv(:, :)
      !> Row, column and summation indices
      integer :: i, c, l
      !> Accumulated right-hand side
      real(wp) :: s

      rinv = 0.0_wp
      do c = 1, size(r, 2)
         do i = c, 1, -1
            s = merge(1.0_wp, 0.0_wp, i == c)
            do l = i + 1, c
               s = s - r(i, l)*rinv(l, c)
            end do
            rinv(i, c) = s/r(i, i)
         end do
      end do
   end subroutine upper_inverse

   !> Jacobian-vector product at the current Newton point
   !>
   !> Exact when the evaluator was linearized there (`jvp_ready`, set by
   !> `jfnk_outer_step`), a finite-difference probe otherwise
   !> Finite-difference accuracy depends on the probe size and residual noise.
   !>
   !> @param[inout] self        Solver (counters, finite-difference scratch)
   !> @param[in]    residual_cb Residual evaluator
   !> @param[in]    gamma_curr  Current state x
   !> @param[in]    f_curr      Residual F(x) at the current state
   !> @param[in]    direction   Direction d (flat)
   !> @param[out]   jv          J(x)*d (flat)
   !> @param[out]   feasible    Whether the product is usable
   !> @param[out]   error       Error handle
   subroutine jfnk_apply_jacobian(self, residual_cb, gamma_curr, f_curr, direction, jv, &
                                   feasible, error)
      !> Solver
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      !> Current state x
      real(wp), intent(in) :: gamma_curr(:, :)
      !> Residual F(x)
      real(wp), intent(in) :: f_curr(:, :)
      !> Direction d
      real(wp), intent(in) :: direction(:)
      !> Product J(x)*d
      real(wp), intent(out) :: jv(:)
      !> Whether the product is usable
      logical, intent(out) :: feasible
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      if (self%jvp_ready) then
         call residual_cb%jvp(direction, jv, feasible, error)
         self%njvp_total = self%njvp_total + 1
      else
         call self%ensure_fd_workspace(gamma_curr)
         call jfnk_fd_jvp(residual_cb, gamma_curr, f_curr, direction, jv, feasible, error, &
                           self%fd_eps, self%workspace%ws_gamma_pert, self%workspace%ws_f_pert, self%nresidual_total, &
                           self%team_size())
      end if
      if (allocated(error)) feasible = .false.
   end subroutine jfnk_apply_jacobian

   !> Finite-difference Jacobian-vector product: jv ~= (F(x + sigma*d) - F(x))/sigma
   !>
   !> Perturbation magnitude follows the JFNK (NITSOL / Brown-Saad) scaling
   !> sigma = sqrt(eps_mach)*(1 + ||x||)/||d||, but with ||x|| taken as the
   !> root-mean-square norm sqrt(sum(x**2)/n) rather than the raw L2 norm
   !>
   !> RMS scaling keeps the L2 displacement ||sigma*d|| independent of the field size
   !> for a fixed RMS state. Per-entry increments and accuracy still depend on the
   !> direction and residual noise; exact products avoid this probe-size choice
   !>
   !> @param[in]  residual_cb Residual evaluator
   !> @param[in]  gamma_curr  current state x
   !> @param[in]  f_curr      residual F(x) at the current state
   !> @param[in]  direction   Krylov direction d (flat, un-normalized)
   !> @param[out] jv          approximation to J(x)*d
   !> @param[out] feasible    Whether the probe succeeded
   !> @param[out] error       Probe failure information
   !> @param[in]  fd_eps      manual step override; if > 0 the probe is a fixed
   !>                         step fd_eps along the unit direction
   !> @param[inout] gamma_pert Perturbed field scratch (npts, npair)
   !> @param[inout] f_pert     Perturbed residual scratch (npts, npair)
   !> @param[inout] nresidual residual-evaluation counter, incremented per probe
   !> @param[in] nt OpenMP team size
   subroutine jfnk_fd_jvp(residual_cb, gamma_curr, f_curr, direction, jv, &
                           feasible, error, fd_eps, gamma_pert, f_pert, nresidual, nt)
      !> Model residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      !> Current point (npts, npair)
      real(wp), intent(in) :: gamma_curr(:, :)
      !> Residual at the current point (npts, npair)
      real(wp), intent(in) :: f_curr(:, :)
      !> Flat direction
      real(wp), intent(in) :: direction(:)
      !> Flat approximate product
      real(wp), intent(out) :: jv(:)
      !> Whether the probe succeeded
      logical, intent(out) :: feasible
      !> Probe failure information
      type(error_type), allocatable, intent(out) :: error
      !> Positive manual step or zero for automatic scaling
      real(wp), intent(in) :: fd_eps
      !> Perturbed field scratch
      real(wp), intent(inout) :: gamma_pert(:, :)
      !> Perturbed residual scratch
      real(wp), intent(inout) :: f_pert(:, :)
      !> Residual-evaluation counter
      integer, intent(inout) :: nresidual
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt

      integer :: npts, npair, i, j
      !> Perturbation scale applied to `direction`, L2 norm of `direction`
      real(wp) :: sigma, norm_d
      !> Root-mean-square norm of the current state vector
      real(wp) :: rms_gamma

      feasible = .true.
      npts = size(gamma_curr, 1)
      npair = size(gamma_curr, 2)

      norm_d = sqrt(blocked_dot(direction, direction, nt))
      if (.not. ieee_is_finite(norm_d)) then
         feasible = .false.
         call fatal_error(error, "Jacobian direction norm must be finite")
         return
      end if
      if (norm_d == 0.0_wp) then
         !$omp parallel do schedule(static) num_threads(merge(nt, 1, size(jv) >= parallel_vector_size))
         do i = 1, size(jv)
            jv(i) = 0.0_wp
         end do
         !$omp end parallel do
         return
      end if

      !> sigma scales the (un-normalized) direction, so
      !> jv = (F(x+sigma*d)-F(x))/sigma recovers J*direction
      !> A positive fd_eps overrides the auto scaling with a fixed step
      !> along the unit direction (sigma = fd_eps/||d||), so that
      !> jv = ||d||*(F(x + fd_eps*d/||d||) - F(x))/fd_eps as before
      if (fd_eps > 0.0_wp) then
         sigma = fd_eps/norm_d
      else
         !> Grid-size-independent RMS norm of the current state
         rms_gamma = sqrt(field_square_sum(gamma_curr, nt)/real(size(gamma_curr), wp))
         sigma = sqrt(epsilon(1.0_wp))*(1.0_wp + rms_gamma)/norm_d
      end if

      !$omp parallel do collapse(2) schedule(static) &
      !$omp& num_threads(merge(nt, 1, size(direction) >= parallel_vector_size))
      do j = 1, npair
         do i = 1, npts
            gamma_pert(i, j) = gamma_curr(i, j) + sigma*direction(i + (j - 1)*npts)
         end do
      end do
      !$omp end parallel do
      call residual_cb%residual(gamma_pert, f_pert, feasible, error)
      nresidual = nresidual + 1
      if (allocated(error)) feasible = .false.
      if (.not. feasible) return

      !$omp parallel do collapse(2) schedule(static) num_threads(merge(nt, 1, size(jv) >= parallel_vector_size))
      do j = 1, npair
         do i = 1, npts
            jv(i + (j - 1)*npts) = (f_pert(i, j) - f_curr(i, j))/sigma
         end do
      end do
      !$omp end parallel do
   end subroutine jfnk_fd_jvp

   !> Reuse compatible scratch, retaining capacity shared by completed phases
   !>
   !> @param[inout] self Solver workspace
   !> @param[in] gamma_mold Field defining the (npts, npair) shape
   subroutine jfnk_ensure_workspace(self, gamma_mold)
      !> Solver workspace
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Field shape
      real(wp), intent(in) :: gamma_mold(:, :)

      integer :: npts, npair, nflat, m, i, j
      logical :: shared
      !> Recycled subspace dimension
      integer :: kmax
      !> Workspace wall-clock start
      integer(int64) :: tick

      npts = size(gamma_mold, 1)
      npair = size(gamma_mold, 2)
      nflat = npts*npair
      !> A standalone workspace follows restart exactly. Shared scratch retains
      !> the larger capacity requested by either phase, avoiding repeated resize.
      !> The logical cycle bound remains self%restart; convergence, recycling
      !> and the remaining max_krylov budget can reduce its Arnoldi steps.
      m = max(1, self%restart)
      kmax = self%recycle_dim()

      shared = .false.
      if (allocated(self%workspace)) then
         shared = self%workspace%shared
         if (self%workspace%ws_nflat == nflat .and. self%workspace%ws_npts == npts .and. &
            & self%workspace%ws_npair == npair) then
            if (allocated(self%workspace%ws_v)) then
               if (shared .and. size(self%workspace%ws_v, 2) >= m + 1) then
                  call self%ensure_recycle_workspace(kmax)
                  return
               end if
               if (size(self%workspace%ws_v, 2) == m + 1) then
                  if (kmax == 0) then
                     if (.not. allocated(self%workspace%ws_u)) return
                  else if (allocated(self%workspace%ws_u)) then
                     if (size(self%workspace%ws_u, 2) == kmax) return
                  end if
               end if
            end if
         end if
      end if

      call system_clock(tick)
      call self%release_workspace()
      allocate (self%workspace)
      self%workspace%shared = shared

      self%workspace%ws_nflat = nflat
      self%workspace%ws_npts = npts
      self%workspace%ws_npair = npair

      allocate (self%workspace%ws_gamma_trial(npts, npair))
      allocate (self%workspace%ws_delta(npts, npair))
      allocate (self%workspace%ws_f_curr(npts, npair))
      allocate (self%workspace%ws_f_trial(npts, npair))

      allocate (self%workspace%ws_rhs_flat(nflat))
      allocate (self%workspace%ws_delta_flat(nflat))

      allocate (self%workspace%ws_r(nflat))
      allocate (self%workspace%ws_w(nflat))

      allocate (self%workspace%ws_v(nflat, m + 1))
      allocate (self%workspace%ws_h(m + 1, m))
      allocate (self%workspace%ws_cs(m))
      allocate (self%workspace%ws_sn(m))
      allocate (self%workspace%ws_g(m + 1))
      allocate (self%workspace%ws_y(m))
      allocate (self%workspace%ws_proj(m))
      allocate (self%workspace%ws_proj2(m))
      allocate (self%workspace%ws_hbar(m + 1, m), source=0.0_wp)
      allocate (self%workspace%ws_bproj(max(1, kmax), m), source=0.0_wp)
      if (kmax > 0) allocate (self%workspace%ws_u(nflat, kmax), self%workspace%ws_c(nflat, kmax))
      self%nrecycled = 0

      !$omp parallel num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size)) private(i, j)
      !$omp do collapse(2) schedule(static)
      do j = 1, npair
         do i = 1, npts
            self%workspace%ws_gamma_trial(i, j) = 0.0_wp
            self%workspace%ws_delta(i, j) = 0.0_wp
            self%workspace%ws_f_curr(i, j) = 0.0_wp
            self%workspace%ws_f_trial(i, j) = 0.0_wp
         end do
      end do
      !$omp end do
      !$omp do schedule(static)
      do i = 1, nflat
         self%workspace%ws_rhs_flat(i) = 0.0_wp
         self%workspace%ws_delta_flat(i) = 0.0_wp
         self%workspace%ws_r(i) = 0.0_wp
         self%workspace%ws_w(i) = 0.0_wp
      end do
      !$omp end do
      ! Match first touch of every Krylov column to its later vector sweeps
      do j = 1, m + 1
         !$omp do schedule(static)
         do i = 1, nflat
            self%workspace%ws_v(i, j) = 0.0_wp
         end do
         !$omp end do
      end do
      do j = 1, kmax
         !$omp do schedule(static)
         do i = 1, nflat
            self%workspace%ws_u(i, j) = 0.0_wp
            self%workspace%ws_c(i, j) = 0.0_wp
         end do
         !$omp end do
      end do
      !$omp end parallel
      self%workspace%ws_h = 0.0_wp
      self%workspace%ws_g = 0.0_wp
      self%workspace%ws_y = 0.0_wp
      call accumulate_wall(tick, self%workspace_wall)
   end subroutine jfnk_ensure_workspace

   !> Add elapsed monotonic wall seconds to a cumulative counter
   !>
   !> @param[in] tick Start clock count
   !> @param[inout] total Accumulated seconds
   subroutine accumulate_wall(tick, total)
      !> Starting system clock count
      integer(int64), intent(in) :: tick
      !> Cumulative elapsed wall seconds
      real(wp), intent(inout) :: total
      !> Current clock count and counts per second
      integer(int64) :: now, rate
      call system_clock(now, rate)
      if (rate > 0_int64) total = total + real(now - tick, wp)/real(rate, wp)
   end subroutine accumulate_wall

   !> Dot product with fixed block and combination order across thread counts
   !>
   !> @param[in] x First vector
   !> @param[in] y Second vector, of the same length
   !> @return Scalar dot product; small vectors retain the intrinsic path
   !> @param[in] nt OpenMP team size
   function blocked_dot(x, y, nt) result(value)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Input vectors
      real(wp), intent(in) :: x(:), y(:)
      !> Dot product
      real(wp) :: value
      !> One partial sum per fixed block
      real(wp), allocatable :: partial(:)
      !> Block index, block count, and first/last element
      integer :: b, nb, first, last

      if (size(x) < parallel_vector_size) then
         value = dot_product(x, y)
         return
      end if
      nb = (size(x) - 1)/reduction_block_size + 1
      allocate (partial(nb))
      !$omp parallel do schedule(static) private(first, last) num_threads(nt)
      do b = 1, nb
         first = (b - 1)*reduction_block_size + 1
         last = min(size(x), first + reduction_block_size - 1)
         partial(b) = dot(x(first:last), y(first:last))
      end do
      !$omp end parallel do
      value = sum(partial)
   end function blocked_dot

   !> Parallel BLAS projection without packing strided matrix tiles
   !>
   !> @param[in] n,m Basis dimensions
   !> @param[in] v,w Basis and vector to project
   !> @param[out] h Projection coefficients, combined in fixed block order
   !> @param[in] nt OpenMP team size
   subroutine project_basis(n, m, v, w, h, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Explicit dimensions preserve the original leading dimension for BLAS
      integer, intent(in) :: n, m
      !> Contiguous basis and vector
      real(wp), intent(in) :: v(n, m), w(n)
      !> Projection coefficients
      real(wp), intent(out) :: h(m)
      !> Small block products, not field-sized packing buffers
      real(wp), allocatable :: partial(:, :)
      !> Block bounds and coefficient index
      integer :: nb, b, first, count, j
      if (n < parallel_vector_size) then
         call gemv(v, w, h, trans="t")
         return
      end if
      nb = (n - 1)/reduction_block_size + 1
      allocate (partial(m, nb))
      !$omp parallel do schedule(static) private(first, count) num_threads(nt)
      do b = 1, nb
         first = (b - 1)*reduction_block_size + 1
         count = min(reduction_block_size, n - first + 1)
         call gemv_rows(v, first, count, w(first:first + count - 1), &
            & partial(:, b), 1.0_wp, 0.0_wp, "T")
      end do
      !$omp end parallel do
      do j = 1, m
         h(j) = sum(partial(j, :))
      end do
   end subroutine project_basis

   !> Cache-blocked BLAS basis update, sharing the destination across columns in cache
   !>
   !> @param[in] n,m Basis dimensions
   !> @param[in] v,h Basis and coefficients
   !> @param[in] alpha Scale factor
   !> @param[inout] w Updated vector w + alpha*V*h
   !> @param[in] nt OpenMP team size
   subroutine update_basis(n, m, v, h, alpha, w, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Explicit dimensions preserve BLAS leading dimensions without packing
      integer, intent(in) :: n, m
      !> Basis, coefficients, and scale factor
      real(wp), intent(in) :: v(n, m), h(m), alpha
      !> Output vector
      real(wp), intent(inout) :: w(n)
      !> Block index and bounds
      integer :: b, first, count, j
      if (n < parallel_vector_size) then
         do j = 1, m
            w = w + alpha*h(j)*v(:, j)
         end do
         return
      end if
      !$omp parallel do schedule(static) private(first, count) num_threads(nt)
      do b = 1, (n - 1)/reduction_block_size + 1
         first = (b - 1)*reduction_block_size + 1
         count = min(reduction_block_size, n - first + 1)
         call gemv_rows(v, first, count, h, w(first:first + count - 1), alpha, 1.0_wp, "N")
      end do
      !$omp end parallel do
   end subroutine update_basis

   !> Squared field norm without flattening or allocating a field-sized temporary
   !>
   !> @param[in] x Input field in native Fortran storage order
   !> @return Sum of squares with fixed block and combination order
   !> @param[in] nt OpenMP team size
   function field_square_sum(x, nt) result(value)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Input field
      real(wp), intent(in) :: x(:, :)
      !> Sum of squares
      real(wp) :: value
      !> Partial sums in spatial-block then column order
      real(wp), allocatable :: partial(:, :)
      !> Block/column indices, block count, and element bounds
      integer :: b, j, nb, first, last

      if (size(x) < parallel_vector_size) then
         value = sum(x**2)
         return
      end if
      nb = (size(x, 1) - 1)/reduction_block_size + 1
      allocate (partial(nb, size(x, 2)))
      !$omp parallel do collapse(2) schedule(static) private(first, last) num_threads(nt)
      do j = 1, size(x, 2)
         do b = 1, nb
            first = (b - 1)*reduction_block_size + 1
            last = min(size(x, 1), first + reduction_block_size - 1)
            partial(b, j) = sum(x(first:last, j)**2)
         end do
      end do
      !$omp end parallel do
      value = sum(partial)
   end function field_square_sum

   !> Copy a flat vector with the same static ownership as Krylov updates
   !>
   !> @param[in] source Input vector
   !> @param[out] dest Copied vector
   !> @param[in] nt OpenMP team size
   subroutine copy_vector(source, dest, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Input vector
      real(wp), intent(in) :: source(:)
      !> Output vector of the same length
      real(wp), intent(out) :: dest(:)
      !> Element index
      integer :: i
      !$omp parallel do schedule(static) num_threads(merge(nt, 1, size(source) >= parallel_vector_size))
      do i = 1, size(source)
         dest(i) = source(i)
      end do
      !$omp end parallel do
   end subroutine copy_vector

   !> Normalize a vector without changing division to reciprocal multiplication
   !>
   !> @param[in] source Input vector
   !> @param[in] divisor Normalization divisor
   !> @param[out] dest Normalized vector
   !> @param[in] nt OpenMP team size
   subroutine divide_vector(source, divisor, dest, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Input vector and divisor
      real(wp), intent(in) :: source(:), divisor
      !> Output vector of the same length
      real(wp), intent(out) :: dest(:)
      !> Element index
      integer :: i
      !$omp parallel do schedule(static) num_threads(merge(nt, 1, size(source) >= parallel_vector_size))
      do i = 1, size(source)
         dest(i) = source(i)/divisor
      end do
      !$omp end parallel do
   end subroutine divide_vector

   !> Unpack a vector in Fortran storage order without a reshape temporary
   !>
   !> @param[in] source Packed vector
   !> @param[out] dest Unpacked field
   !> @param[in] nt OpenMP team size
   subroutine unpack_vector(source, dest, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Packed input vector
      real(wp), intent(in) :: source(:)
      !> Output field with size equal to the vector length
      real(wp), intent(out) :: dest(:, :)
      !> Spatial and pair indices
      integer :: i, j
      !$omp parallel do collapse(2) schedule(static) num_threads(merge(nt, 1, size(source) >= parallel_vector_size))
      do j = 1, size(dest, 2)
         do i = 1, size(dest, 1)
            dest(i, j) = source(i + (j - 1)*size(dest, 1))
         end do
      end do
      !$omp end parallel do
   end subroutine unpack_vector

   !> Form a trial iterate with independent elementwise updates
   !>
   !> @param[in] gamma Current field
   !> @param[in] delta Newton direction
   !> @param[in] alpha Step length
   !> @param[out] trial Trial field
   !> @param[in] nt OpenMP team size
   subroutine trial_field(gamma, delta, alpha, trial, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Current field, direction, and step length
      real(wp), intent(in) :: gamma(:, :), delta(:, :), alpha
      !> Trial field of matching shape
      real(wp), intent(out) :: trial(:, :)
      !> Spatial and pair indices
      integer :: i, j
      !$omp parallel do collapse(2) schedule(static) num_threads(merge(nt, 1, size(gamma) >= parallel_vector_size))
      do j = 1, size(gamma, 2)
         do i = 1, size(gamma, 1)
            trial(i, j) = gamma(i, j) - alpha*delta(i, j)
         end do
      end do
      !$omp end parallel do
   end subroutine trial_field

   !> Release solver workspace while retaining borrowed operators
   !>
   !> @param[inout] self Solver workspace to release
   subroutine jfnk_release_workspace(self)
      !> Solver whose owned scratch is released
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      if (allocated(self%workspace)) deallocate (self%workspace)
      self%nrecycled = 0
   end subroutine jfnk_release_workspace

   !> Apply the lent preconditioner to vec_flat in place
   !>
   !> Only called for a step whose `ready` check passed, so the operator is attached
   !>
   !> @param[inout] self        Solver instance
   !> @param[inout] vec_flat    Flat (npts*npair) vector, preconditioned in place
   subroutine jfnk_apply_precond(self, vec_flat)
      !> Solver instance
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Flat vector, modified in place
      real(wp), contiguous, target, intent(inout) :: vec_flat(:)

      call self%precond%apply(vec_flat)
   end subroutine jfnk_apply_precond

   !> Form a Givens rotation that removes b without an absolute scale cutoff
   !>
   !> @param[in] a First entry
   !> @param[in] b Entry to eliminate
   !> @param[out] c Rotation cosine
   !> @param[out] s Rotation sine
   subroutine generate_givens(a, b, c, s)
      !> Entries to rotate
      real(wp), intent(in) :: a, b
      !> Rotation cosine and sine
      real(wp), intent(out) :: c, s
      real(wp) :: r

      if (b == 0.0_wp) then
         c = 1.0_wp
         s = 0.0_wp
      else if (abs(b) > abs(a)) then
         r = -a/b
         s = 1.0_wp/sqrt(1.0_wp + r**2)
         c = s*r
      else
         r = -b/a
         c = 1.0_wp/sqrt(1.0_wp + r**2)
         s = c*r
      end if
   end subroutine generate_givens

   subroutine apply_givens(c, s, a, b)
      real(wp), intent(in) :: c, s
      real(wp), intent(inout) :: a, b
      real(wp) :: tmp

      tmp = c*a - s*b
      b = s*a + c*b
      a = tmp
   end subroutine apply_givens

   !> Solve the small triangular Krylov system without dividing by a zero pivot
   !>
   !> @param[in] r Upper triangular matrix
   !> @param[in] rhs Rotated least-squares right-hand side
   !> @param[out] x Solution; only usable when feasible
   !> @param[out] feasible Nonsingular finite system and finite solution
   subroutine upper_backsolve(r, rhs, x, feasible)
      !> Upper triangular matrix
      real(wp), intent(in) :: r(:, :)
      !> Rotated right-hand side
      real(wp), intent(in) :: rhs(:)
      !> Solution coefficients
      real(wp), intent(out) :: x(:)
      !> Whether the solution is usable
      logical, intent(out) :: feasible

      integer :: n, i, j
      real(wp) :: sum_val

      n = size(rhs)
      feasible = .false.
      if (.not. all(ieee_is_finite(r)) .or. .not. all(ieee_is_finite(rhs))) return
      do i = n, 1, -1
         if (r(i, i) == 0.0_wp) return
         sum_val = rhs(i)
         do j = i + 1, n
            sum_val = sum_val - r(i, j)*x(j)
         end do
         x(i) = sum_val/r(i, i)
         if (.not. ieee_is_finite(x(i))) return
      end do
      feasible = .true.
   end subroutine upper_backsolve

end module moist_model_moz_solver_jfnk
