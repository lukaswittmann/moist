!> Jacobian-free Newton-Krylov solver with restarted GMRES
!>
!> - Uses finite differences or evaluator-supplied products without forming a Jacobian
module moist_model_moz_solver_jfnk
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_output_citations, only: citation
   use moist_model_moz_solver_core_adjoint, only: moist_model_moz_solver_jfnk_adjoint_parameters_type
   use moist_model_moz_solver_core_base, only: moist_model_moz_solver_type, moist_model_moz_solver_statistics_type
   use moist_model_moz_solver_core_checks, only: finite_positive, finite_nonnegative, unit_fraction, step_fraction
   use moist_model_moz_solver_core_printer, only: moist_model_moz_solver_printer_type
   use moist_model_moz_solver_core_preconditioner, only: moist_model_moz_solver_preconditioner_type
   use moist_model_moz_solver_core_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_moz_solver_jfnk_deflation, only: moist_model_moz_solver_deflation_type
   use moist_model_moz_solver_jfnk_workspace, only: moist_model_moz_solver_workspace_type, exchange_solver_workspace
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_model_moz_solver_jfnk_kernels, only: parallel_vector_size
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_jfnk_type, moist_model_moz_solver_jfnk_parameters_type
   public :: moist_model_moz_solver_jfnk_adjoint_parameters_type, new_jfnk_solver

   !> Inner-solve, forcing and line-search controls
   !>
   !> JSON/TOML keys: restart, max_krylov, cgs2, recycle, stagnation_tol, eta0, eta_min, eta_max, forcing_beta,
   !> forcing_power, jvp, fd_eps, deflation, deflation_rmse, line_search_c, line_search_tau, line_search_min_alpha,
   !> line_search_max_iter, line_search_expand_ratio, alpha_max, fallback_alpha
   type, extends(moist_model_parameters_type) :: moist_model_moz_solver_jfnk_parameters_type

      !* ------------------------------ GMRES inner solve ----------------------------- *!

      !> GMRES restart length: maximum Krylov dimension of one cycle
      !>
      !> - Convergence or max_krylov may end a cycle early
      !> - Shared workspace may retain capacity for the larger of primal/transpose restarts
      integer :: restart = 50
      !> Inner GMRES iterations per Newton update; the driver limits outer iterations
      integer :: max_krylov = 50
      !> Classical Gram-Schmidt with one reorthogonalization instead of MGS
      !>
      !> - Fewer vector sweeps, with speed depending on field size and threading
      !> - Rounding can change the nonlinear convergence path, so it defaults off
      logical :: cgs2 = .false.
      !> Recycled Krylov vectors (GCRO-DR); zero disables recycling
      !>
      !> - The solver caps it by restart and max_krylov to leave at least one Arnoldi step
      integer :: recycle = 0
      !> Relative decrease below which a restart cycle ends the inner solve; zero disables it
      !>
      !> - Zero lets GMRES run to eta*||F|| or the iteration cap, avoiding a premature
      !>   bail on non-symmetric MOZ Jacobian stagnation
      real(wp) :: stagnation_tol = 0.0_wp

      !* --------------------------- Inexact Newton forcing --------------------------- *!

      !> Forcing term of the first two Newton updates of a stage
      real(wp) :: eta0 = 0.5_wp
      !> Lower bound of later forcing terms
      real(wp) :: eta_min = 1.0e-3_wp
      !> Upper bound of later forcing terms
      real(wp) :: eta_max = 0.1_wp
      !> Scale beta of the Eisenstat-Walker forcing term beta*(rmse/rmse_prev)**forcing_power
      real(wp) :: forcing_beta = 0.9_wp
      !> Exponent of the Eisenstat-Walker forcing term
      real(wp) :: forcing_power = 2.0_wp

      !* -------------------------- Jacobian-vector products -------------------------- *!

      !> Exact evaluator Jacobian-vector products instead of finite differences
      !>
      !> - true: exact evaluator products; a step without them fails
      !> - false: finite-difference products
      logical :: use_exact_jvp = .true.
      !> Finite-difference displacement; zero selects an automatic value
      real(wp) :: fd_eps = 0.0_wp

      !* --------------------------- Coarse-space deflation --------------------------- *!

      !> Coarse-space (DEF1) deflation of the inner solves
      logical :: deflation = .false.
      !> Deflation is used below this RMS residual
      real(wp) :: deflation_rmse = 1.0e-3_wp

      !* --------------------------------- Line search -------------------------------- *!

      !> Armijo sufficient-decrease constant of the line search
      real(wp) :: line_search_c = 1.0e-4_wp
      !> Backtracking factor applied to a rejected step length
      real(wp) :: line_search_tau = 0.5_wp
      !> Step length below which the line search stops backtracking
      real(wp) :: line_search_min_alpha = 1.0e-4_wp
      !> Hard cap on trial steps (residual evaluations) of the backtracking line search
      !>
      !> - Rejected trials shrink alpha to [0.1, line_search_tau]*alpha
      !> - Cap may stop the search before line_search_min_alpha is reached
      !> - Step extension and the Picard fallback have their own trial loops
      integer :: line_search_max_iter = 10
      !> Extend only when the full step keeps more than this merit fraction
      !>
      !> - Merit is phi = 0.5*||F||^2
      !> - A stronger decrease is the fast Newton regime
      !> - There a longer step overshoots and its trial would be wasted
      real(wp) :: line_search_expand_ratio = 0.25_wp
      !> Longest line-search fraction beyond the full Newton step; one disables step extension
      real(wp) :: alpha_max = 2.0_wp
      !> Initial step length of the Picard fallback after a failed line search
      real(wp) :: fallback_alpha = 0.25_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => jfnk_init_defaults
      !> Declare fields for JSON/TOML input, output and printing
      procedure :: register_entries => jfnk_register_entries
      !> Check the ranges of all controls
      procedure :: validate => jfnk_validate
   end type moist_model_moz_solver_jfnk_parameters_type

   !> Jacobian-free Newton-Krylov solver with restarted GMRES
   !>
   !> - Operates on fields stored as (npts, npair) arrays
   type, extends(moist_model_moz_solver_type) :: moist_model_moz_solver_jfnk_type
      !> Inner-solve, forcing and line-search controls
      type(moist_model_moz_solver_jfnk_parameters_type) :: params
      !> Whether the evaluator holds a linearization at the current Newton point
      logical :: jvp_ready = .false.
      !> Last update was a zero step at eta_min; a second one is a stall
      logical :: stalled_at_min = .false.
      !> Number of recycled vectors currently held in ws_u (0 = none yet)
      integer :: nrecycled = 0
      !> Whether the recycled subspace was built for the preconditioned operator
      logical :: recycle_precond = .false.
      !> Forcing term for the next update, shared by outer_step calls within one solve stage
      !>
      !> - Negative value selects eta0 on the next update
      real(wp) :: eta_current = -1.0_wp
      !> Whether the preconditioner reported ready for the current update
      logical :: precond_valid_state = .false.
      !> Input residual RMS of the previous update; huge at the start of a stage
      real(wp) :: rmse_prev_state = huge(1.0_wp)
      !> Last inner work count: Arnoldi steps plus recycled-image refresh products
      integer :: gmres_iters_last = 0
      !> Cumulative inner work in this stage, including recycled refreshes and failed solves
      integer :: gmres_iters_total = 0
      !> Cumulative inner work, residual evaluations and exact Jacobian-vector products since construction
      !>
      !> - Residual evaluations cover Newton residuals, line-search trials and finite-difference probes
      integer :: nkrylov_total = 0, nresidual_total = 0, njvp_total = 0
      !> Attempted outer_step calls in the current stage
      integer :: newton_iter_count = 0
      !> Stage counts for the summary note
      !>
      !> - Newton steps that took the Picard fallback, and those that returned a zero step
      !> - Inner solves capped at max_krylov short of their tolerance, and stagnation exits
      !> - Steps whose requested deflation stayed inactive, and deflated solves that skipped recycling
      integer :: nfallback_stage = 0, nzero_stage = 0, ncapped_stage = 0, nstagnated_stage = 0
      integer :: nundeflated_stage = 0, nunrecycled_stage = 0
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
      !> Last reported deflation state
      logical :: deflation_status = .false.
   contains
      procedure :: get_statistics => jfnk_get_statistics
      procedure :: stage_summary => jfnk_stage_summary
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
   end type moist_model_moz_solver_jfnk_type

   interface

      !* ============================================================================== *!
      !*                             Newton step (step.f90)                             *!
      !* ============================================================================== *!

      !> One Newton-Krylov outer step: evaluate residual, GMRES linear solve, Armijo line search
      !>
      !> - Called once per Newton iteration, in a loop driven by the caller
      !>
      !> @param[inout] self Solver state and workspace
      !> @param[in] residual_fn Model residual evaluator
      !> @param[in] gamma Current field (npts, npair)
      !> @param[out] gamma_new Line-search update (npts, npair)
      !> @param[out] rmse RMS residual of the input field
      !> @param[out] feasible Whether the update is usable
      !> @param[out] error Evaluation failure
      module subroutine jfnk_outer_step(self, residual_fn, gamma, gamma_new, rmse, feasible, error)
         implicit none(type, external)
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
      end subroutine jfnk_outer_step

      !* ============================================================================== *!
      !*                        Jacobian products (jacobian.f90)                        *!
      !* ============================================================================== *!

      !> Jacobian-vector product at the current Newton point
      !>
      !> Exact when the evaluator was linearized there (`jvp_ready`, set by
      !> `jfnk_outer_step`), a finite-difference probe otherwise
      !> Finite-difference accuracy depends on the probe size and residual noise
      !>
      !> @param[inout] self        Solver (counters, finite-difference scratch)
      !> @param[in]    residual_cb Residual evaluator
      !> @param[in]    gamma_curr  Current state x
      !> @param[in]    f_curr      Residual F(x) at the current state
      !> @param[in]    direction   Direction d (flat)
      !> @param[out]   jv          J(x)*d (flat)
      !> @param[out]   feasible    Whether the product is usable
      !> @param[out]   error       Error handle
      module subroutine jfnk_apply_jacobian(self, residual_cb, gamma_curr, f_curr, direction, jv, &
                                            feasible, error)
         implicit none(type, external)
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
      end subroutine jfnk_apply_jacobian

      !* ============================================================================== *!
      !*                        GMRES and deflation (gmres.f90)                         *!
      !* ============================================================================== *!

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
      !> With an active coarse space (`deflation%active`) the solve is DEF1 deflated
      !> (see moist_model_moz_solver_jfnk_deflation)
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
      module subroutine gmres_linear_solve(self, residual_cb, gamma_curr, f_curr, rhs_flat, delta_flat, &
                                           eta, precond_valid, feasible, error, gmres_iters)
         implicit none(type, external)
         !> Solver (workspace, recycled subspace)
         class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
         !> Residual evaluator
         class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
         !> Current Newton point (npts, npair)
         real(wp), intent(in) :: gamma_curr(:, :)
         !> Residual at the current Newton point (npts, npair)
         real(wp), intent(in) :: f_curr(:, :)
         !> Right-hand side F, flat (npts*npair)
         real(wp), intent(in) :: rhs_flat(:)
         !> Newton direction, flat (npts*npair)
         real(wp), intent(out) :: delta_flat(:)
         !> Forcing term (relative tolerance)
         real(wp), intent(in) :: eta
         !> Whether the left preconditioner is applied
         logical, intent(inout) :: precond_valid
         !> Whether every product was usable
         logical, intent(out) :: feasible
         !> Error handle
         type(error_type), allocatable, intent(out) :: error
         !> Arnoldi iterations plus recycled-image products
         integer, intent(out), optional :: gmres_iters
      end subroutine gmres_linear_solve

      !> Deflate a residual in place, vec := (I - J Q) vec
      !>
      !> - Q = Z E^{-1} Z^T Omega
      !> - Costs one Jacobian product
      !> - Uses ws_z and ws_w, so vec must be neither
      !>
      !> @param[inout] self        Solver with an active deflation
      !> @param[in]    residual_cb Residual evaluator
      !> @param[in]    gamma_curr  Current Newton point
      !> @param[in]    f_curr      Residual at the current Newton point
      !> @param[inout] vec         Flat vector to deflate
      !> @param[out]   y           Coarse coordinates E^{-1} Z^T Omega vec, (kz, 1)
      !> @param[out]   feasible    Whether the product was usable
      !> @param[out]   error       Error handle
      module subroutine gmres_deflate(self, residual_cb, gamma_curr, f_curr, vec, y, feasible, error)
         implicit none(type, external)
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
      end subroutine gmres_deflate

      !> DEF1 coarse correction of a deflated solve
      !>
      !> - delta := x + Q (b - J x)
      !> - Equals Q b + (I - Q J) x
      !> - Costs one Jacobian product unless x = 0
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
      module subroutine gmres_deflation_finish(self, residual_cb, gamma_curr, f_curr, rhs_flat, delta_flat, &
                                               x_zero, y, feasible, error)
         implicit none(type, external)
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
      end subroutine gmres_deflation_finish

      !* ============================================================================== *!
      !*                         Krylov recycling (recycle.f90)                         *!
      !* ============================================================================== *!

      !> Recycled subspace dimension in use
      !>
      !> - `recycle`, capped to keep one Arnoldi step per cycle and capped solve
      !> - Zero means plain restarted GMRES
      !>
      !> @param[in] self  solver
      pure module function gmres_recycle_dim(self) result(kmax)
         implicit none(type, external)
         !> Solver
         class(moist_model_moz_solver_jfnk_type), intent(in) :: self
         !> Recycled vectors in use
         integer :: kmax
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
      module subroutine gmres_recycle_refresh(self, residual_cb, gamma_curr, f_curr, precond_valid, njv, &
                                              feasible, error)
         implicit none(type, external)
         !> Solver
         class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
         !> Residual evaluator
         class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
         !> Current Newton point
         real(wp), intent(in) :: gamma_curr(:, :)
         !> Residual at the current Newton point
         real(wp), intent(in) :: f_curr(:, :)
         !> Whether the left preconditioner is applied
         logical, intent(in) :: precond_valid
         !> Jacobian-vector products spent
         integer, intent(out) :: njv
         !> Whether every product was usable
         logical, intent(out) :: feasible
         !> Error handle
         type(error_type), allocatable, intent(out) :: error
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
      module subroutine gmres_recycle_harvest(self, kr, j, need_c)
         implicit none(type, external)
         !> Solver
         class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
         !> Recycled vectors used in the cycle
         integer, intent(in) :: kr
         !> Arnoldi steps of the cycle
         integer, intent(in) :: j
         !> Whether C must be updated too
         logical, intent(in) :: need_c
      end subroutine gmres_recycle_harvest

   end interface

contains

   !> Move scratch between phases and invalidate operator-dependent Krylov data
   !>
   !> @param[inout] self       solver; numerical controls and counters are retained
   !> @param[inout] workspace  owned scratch of a distinct solve context
   subroutine jfnk_exchange_workspace(self, workspace)
      !> Solver; numerical controls and counters are retained
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Owned scratch of a distinct solve context
      type(moist_model_moz_solver_workspace_type), allocatable, intent(inout) :: workspace

      call exchange_solver_workspace(self%workspace, workspace)
      self%nrecycled = 0
      self%recycle_precond = .false.
      self%jvp_ready = .false.
      call self%deflation%deactivate()
   end subroutine jfnk_exchange_workspace

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
   !> - Needs an attached printer with inner rows; the status is tracked anyway
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

      show = .false.
      if (associated(self%iter_printer)) show = self%iter_printer%inner_rows
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

   !> Invalidate stage forcing, recycling and stall state
   !>
   !> - Workspace and borrowed operators are retained
   !> - Next outer_step reinitializes the remaining Newton state
   !> - Stage counts and once-per-stage warnings are reset
   !> - Cumulative counters persist
   !>
   !> @param[inout] self Solver starting a solve or coupling stage
   subroutine jfnk_begin_stage(self)
      !> Solver retaining its workspaces
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      self%eta_current = -1.0_wp
      self%gmres_iters_total = 0
      self%nrecycled = 0
      self%stalled_at_min = .false.
      self%newton_iter_count = 0
      self%nfallback_stage = 0
      self%nzero_stage = 0
      self%ncapped_stage = 0
      self%nstagnated_stage = 0
      self%nundeflated_stage = 0
      self%nunrecycled_stage = 0
   end subroutine jfnk_begin_stage

   !> Note on the fallbacks of the current stage, for the iteration table
   !>
   !> - Picard-fallback steps out of the attempted Newton steps, with the zero steps among them
   !> - Inner solves capped at max_krylov short of their tolerance, and stagnation exits
   !> - Steps left undeflated by a failed or singular coarse matrix, and deflated
   !>   solves that skipped recycling
   !> - Empty when no step or inner solve fell back
   !>
   !> @param[in] self Solver after a stage
   pure function jfnk_stage_summary(self) result(text)
      !> Solver after a stage
      class(moist_model_moz_solver_jfnk_type), intent(in) :: self
      !> Note text, e.g. "JFNK: 7 of 40 steps used the Picard fallback"
      character(len=:), allocatable :: text
      !> One formatted count
      character(len=64) :: part

      text = ""
      if (self%nfallback_stage > 0) then
         write (part, "(i0,' of ',i0,' steps used the Picard fallback')") self%nfallback_stage, &
            & self%newton_iter_count
         text = trim(part)
         if (self%nzero_stage > 0) then
            write (part, "(', ',i0,' with a zero step')") self%nzero_stage
            text = text//trim(part)
         end if
      end if
      if (self%ncapped_stage > 0) then
         write (part, "(i0,1x,a,' hit max_krylov')") self%ncapped_stage, inner_solves(self%ncapped_stage)
         call append_part(text, trim(part))
      end if
      if (self%nstagnated_stage > 0) then
         write (part, "(i0,1x,a,' stagnated')") self%nstagnated_stage, inner_solves(self%nstagnated_stage)
         call append_part(text, trim(part))
      end if
      if (self%nundeflated_stage > 0) then
         write (part, "('deflation inactive in ',i0,' of ',i0,' steps')") self%nundeflated_stage, &
            & self%newton_iter_count
         call append_part(text, trim(part))
      end if
      if (self%nunrecycled_stage > 0) then
         write (part, "(i0,1x,a,' skipped recycling under deflation')") self%nunrecycled_stage, &
            & inner_solves(self%nunrecycled_stage)
         call append_part(text, trim(part))
      end if
      if (len(text) > 0) text = self%name()//": "//text

   contains

      !> Singular or plural noun of an inner-solve count
      !>
      !> @param[in] n  count
      pure function inner_solves(n) result(noun)
         !> Count
         integer, intent(in) :: n
         !> Noun phrase
         character(len=:), allocatable :: noun

         if (n == 1) then
            noun = "inner solve"
         else
            noun = "inner solves"
         end if
      end function inner_solves

      !> Join a part to the note, separated by a semicolon
      !>
      !> @param[inout] note  note so far
      !> @param[in]    item  part to append
      pure subroutine append_part(note, item)
         !> Note so far
         character(len=:), allocatable, intent(inout) :: note
         !> Part to append
         character(len=*), intent(in) :: item

         if (len(note) > 0) then
            note = note//"; "//item
         else
            note = item
         end if
      end subroutine append_part

   end function jfnk_stage_summary

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
      type(moist_model_moz_solver_printer_type), target, optional, intent(inout) :: printer
      nullify(self%iter_printer)
      if (present(printer)) then
         self%iter_printer => printer
         call printer%set_columns(diis=.false., gmres=.true.)
      end if
      if (allocated(self%precond_status)) deallocate (self%precond_status)
      self%deflation_status = .false.
   end subroutine jfnk_set_iteration_printer

   !> Whether this solver supports a lent preconditioner; always true, does not mean one is active
   !>
   !> - Attachment and readiness are checked for each linear solve
   !>
   !> @param[in] self Solver to query
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
   function jfnk_get_statistics(self) result(stats)
      !> Solver to inspect
      class(moist_model_moz_solver_jfnk_type), intent(in) :: self
      !> Counters, or unavailable sentinels
      type(moist_model_moz_solver_statistics_type) :: stats
      stats = moist_model_moz_solver_statistics_type()
      stats%ngmres = self%gmres_iters_last
      stats%nkrylov = self%nkrylov_total
      stats%nresidual = self%nresidual_total
      stats%njvp = self%njvp_total
   end function jfnk_get_statistics

   !> Restore compiled defaults
   !>
   !> @param[inout] self Parameter values
   subroutine jfnk_init_defaults(self)
      !> Parameter values
      class(moist_model_moz_solver_jfnk_parameters_type), intent(inout) :: self
      type(moist_model_moz_solver_jfnk_parameters_type) :: defaults

      self%restart = defaults%restart
      self%max_krylov = defaults%max_krylov
      self%cgs2 = defaults%cgs2
      self%recycle = defaults%recycle
      self%stagnation_tol = defaults%stagnation_tol
      self%eta0 = defaults%eta0
      self%eta_min = defaults%eta_min
      self%eta_max = defaults%eta_max
      self%forcing_beta = defaults%forcing_beta
      self%forcing_power = defaults%forcing_power
      self%use_exact_jvp = defaults%use_exact_jvp
      self%fd_eps = defaults%fd_eps
      self%deflation = defaults%deflation
      self%deflation_rmse = defaults%deflation_rmse
      self%line_search_c = defaults%line_search_c
      self%line_search_tau = defaults%line_search_tau
      self%line_search_min_alpha = defaults%line_search_min_alpha
      self%line_search_max_iter = defaults%line_search_max_iter
      self%line_search_expand_ratio = defaults%line_search_expand_ratio
      self%alpha_max = defaults%alpha_max
      self%fallback_alpha = defaults%fallback_alpha
   end subroutine jfnk_init_defaults

   !> Declare parameter fields for JSON/TOML input, output and printing
   !>
   !> @param[inout] self Parameter values
   subroutine jfnk_register_entries(self)
      !> Parameter values
      class(moist_model_moz_solver_jfnk_parameters_type), intent(inout), target :: self

      call self%register_int_scalar("restart", self%restart, level=2)
      call self%register_int_scalar("max_krylov", self%max_krylov, level=2)
      call self%register_logical("cgs2", self%cgs2)
      call self%register_int_scalar("recycle", self%recycle)
      call self%register_real_scalar("stagnation_tol", self%stagnation_tol)
      call self%register_real_scalar("eta0", self%eta0)
      call self%register_real_scalar("eta_min", self%eta_min)
      call self%register_real_scalar("eta_max", self%eta_max)
      call self%register_real_scalar("forcing_beta", self%forcing_beta)
      call self%register_real_scalar("forcing_power", self%forcing_power)
      call self%register_logical("use_exact_jvp", self%use_exact_jvp)
      call self%register_real_scalar("fd_eps", self%fd_eps)
      call self%register_logical("deflation", self%deflation, level=2)
      call self%register_real_scalar("deflation_rmse", self%deflation_rmse)
      call self%register_real_scalar("line_search_c", self%line_search_c)
      call self%register_real_scalar("line_search_tau", self%line_search_tau)
      call self%register_real_scalar("line_search_min_alpha", self%line_search_min_alpha)
      call self%register_int_scalar("line_search_max_iter", self%line_search_max_iter)
      call self%register_real_scalar("line_search_expand_ratio", self%line_search_expand_ratio)
      call self%register_real_scalar("alpha_max", self%alpha_max)
      call self%register_real_scalar("fallback_alpha", self%fallback_alpha)
   end subroutine jfnk_register_entries

   !> Check the ranges of all controls
   !>
   !> - Deflation controls are checked only when enabled
   !> - Deflation needs exact Jacobian-vector products
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error First invalid control
   subroutine jfnk_validate(self, error)
      !> Parameter values
      class(moist_model_moz_solver_jfnk_parameters_type), intent(inout) :: self
      !> First invalid control
      type(error_type), allocatable, intent(out) :: error

      if (self%restart < 1 .or. self%max_krylov < 1 .or. self%recycle < 0) then
         call fatal_error(error, "JFNK needs positive Krylov limits and a nonnegative recycling dimension")
      else if (.not. unit_fraction(self%eta0) .or. .not. unit_fraction(self%eta_min) .or. &
         & .not. unit_fraction(self%eta_max) .or. self%eta_min > self%eta_max) then
         call fatal_error(error, "JFNK forcing terms must lie between zero and one, with eta_min <= eta_max")
      else if (.not. step_fraction(self%forcing_beta) .or. .not. ieee_is_finite(self%forcing_power) .or. &
         & self%forcing_power <= 1.0_wp .or. self%forcing_power > 2.0_wp) then
         call fatal_error(error, "JFNK forcing_beta must lie in (0, 1] and forcing_power in (1, 2]")
      else if (.not. finite_nonnegative(self%fd_eps) .or. .not. finite_nonnegative(self%stagnation_tol) .or. &
         & self%stagnation_tol >= 1.0_wp) then
         call fatal_error(error, "Invalid JFNK finite-difference or stagnation controls")
      else if (.not. unit_fraction(self%line_search_c) .or. .not. unit_fraction(self%line_search_tau) .or. &
         & .not. step_fraction(self%line_search_min_alpha) .or. self%line_search_max_iter < 1 .or. &
         & .not. finite_nonnegative(self%line_search_expand_ratio) .or. self%line_search_expand_ratio > 1.0_wp .or. &
         & .not. ieee_is_finite(self%alpha_max) .or. self%alpha_max < 1.0_wp .or. &
         & .not. step_fraction(self%fallback_alpha)) then
         call fatal_error(error, "Invalid JFNK line-search controls")
      else if (self%deflation .and. .not. finite_positive(self%deflation_rmse)) then
         call fatal_error(error, "JFNK deflation threshold must be finite and positive")
      else if (self%deflation .and. .not. self%use_exact_jvp) then
         call fatal_error(error, "jfnk.deflation needs jfnk.use_exact_jvp")
      end if
   end subroutine jfnk_validate

   !> Construct a JFNK solver from its controls
   !>
   !> - Controls are expected to pass `validate`; no clamping here
   !> - Warns once when restart or max_krylov caps `recycle`
   !>
   !> @param[in] ctx Run context, borrowed for the lifetime of the solver
   !> @param[in] parameters Inner-solve, forcing and line-search controls; defaults when absent
   function new_jfnk_solver(ctx, parameters) result(self)
      !> Run context, borrowed for the lifetime of the solver
      type(moist_context_type), intent(in), target :: ctx
      !> Inner-solve, forcing and line-search controls
      type(moist_model_moz_solver_jfnk_parameters_type), intent(in), optional :: parameters
      !> Constructed solver
      type(moist_model_moz_solver_jfnk_type) :: self
      !> Recycled vectors in use
      integer :: kmax
      !> Formatted requested and used dimensions
      character(len=16) :: requested, used

      self%ctx => ctx
      call self%ctx%cite(citation%jfnk)
      call self%ctx%cite(citation%gmres)
      if (present(parameters)) self%params = parameters

      ! Recycling leaves at least one Arnoldi step per cycle and per capped solve
      kmax = self%recycle_dim()
      if (kmax < self%params%recycle) then
         write (requested, "(i0)") self%params%recycle
         write (used, "(i0)") kmax
         if (kmax == 0) then
            call ctx%warning(self%name()//": restart and max_krylov leave no room for recycling, "// &
               & "ignoring jfnk.recycle = "//trim(requested))
         else
            call ctx%warning(self%name()//": jfnk.recycle = "//trim(requested)//" capped at "//trim(used)// &
               & " by restart and max_krylov")
         end if
      end if
   end function new_jfnk_solver

   !> Reuse compatible scratch, retaining capacity shared by completed phases
   !>
   !> @param[inout] self Solver workspace
   !> @param[in] gamma_mold Field defining the (npts, npair) shape
   subroutine jfnk_ensure_workspace(self, gamma_mold)
      !> Solver workspace
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      !> Field shape
      real(wp), intent(in) :: gamma_mold(:, :)

      !> Field shape, flat length and restart dimension
      integer :: npts, npair, nflat, m
      !> Recycled subspace dimension
      integer :: kmax
      !> OpenMP team size for the first touch
      integer :: nthreads
      !> Whether the workspace replaced storage
      logical :: reallocated

      npts = size(gamma_mold, 1)
      npair = size(gamma_mold, 2)
      nflat = npts*npair
      m = max(1, self%params%restart)
      kmax = self%recycle_dim()

      nthreads = merge(self%team_size(), 1, nflat >= parallel_vector_size)

      if (.not. allocated(self%workspace)) allocate (self%workspace)
      reallocated = self%workspace%ensure(npts, npair, m, kmax, nthreads)
      ! Drop an invalidated or oversized active basis
      if (reallocated .or. self%nrecycled > kmax) self%nrecycled = 0
   end subroutine jfnk_ensure_workspace

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

end module moist_model_moz_solver_jfnk
