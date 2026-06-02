
!> Project anchor points onto an implicit level-set-function (LSF) surface.
!>
!> Given an anchor point, this module finds the nearest point(s) on the
!> surface LSF(x) = 0 by minimizing the quadratic anchor objective
!> phi(x) = phi_alpha/2 * ||x - anchor||^2 subject to LSF(x) = 0. Multiple
!> surface branches (e.g. near concave triple-junctions of the LSF) are
!> enumerated, Newton-refined, classified, and deduplicated.
!>
!> The procedures below are grouped by the section banners they sit under:
!>   1. Augmented-Lagrangian solver callbacks (Newton / KKT system)
!>   2. Initialization and SSD screening
!>   3. Solver callbacks (diagnostics, displacement control, SLSQP)
!>   4. Seed generation: multi-start solver drivers
!>   5. Projection pipeline (seed refinement and dispatch)
!>   6. Candidate filtering and deduplication
!>   7. Newton refinement of a single candidate
!>   8. Check constrained optimality
!>   9. Riemannian Newton escape from saddle points
!>  10. Surface retraction
!>  11. Cleanup and finalization
module moist_cavity_drop_projector
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use iso_fortran_env, only: output_unit
   use mctc_io, only: structure_type
   use mctc_io_convert, only: aatoau
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_objective_phi, only: moist_cavity_drop_objective_phi_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_type, only: solver_base_type
   use moist_math_solver_newton, only: new_newton_solver
   use moist_math_solver_slsqp, only: new_slsqp_solver
   use moist_math_solver_slsqp_multi_tangent, only: new_slsqp_multi_tangent_solver, &
                                                    moist_math_solver_slsqp_multi_tangent_type
   use moist_math_solver_slsqp_multistart, only: new_slsqp_multistart_solver, &
                                                 moist_math_solver_slsqp_multistart_type
   use moist_math_solver_slsqp_deflation, only: new_slsqp_deflation_solver, &
                                                moist_math_solver_slsqp_deflation_type
   use moist_math_solver_newton_deflation, only: new_newton_deflation_solver, &
                                                 moist_math_solver_newton_deflation_type
   use moist_cavity_drop_types, only: projection_workspace_type
   use moist_math_cell_grid, only: moist_cell_grid_type
   use moist_math_linalg, only: setup_tangent_frame, eig_2x2_symmetric
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none
   private

   public :: drop_projector_type
   !> Projector type for efficient grid point projection onto LSF surface

   !> Context type for thread-safe projection callbacks
   !> Contains all data needed by solver callbacks without module-level variables
   type :: projection_context_type
      type(drop_projector_type), pointer :: projector => null()  !> Pointer to parent projector
      class(moist_cavity_drop_lsf_type), pointer :: lsf => null() !> Polymorphic LSF
      real(wp) :: anchor(3)                                     !> Anchor point for projection
      integer :: owner                                          !> Owner index for projection
      real(wp) :: threshold                                     !> Maximum displacement threshold
   end type projection_context_type

   type :: drop_projector_type

      !*------------------------ Solver configuration------------------------ *!

      integer :: proj_level = 3

      !> Tolerances
      real(wp) :: multistart_tol = 1.0e-9_wp
      real(wp) :: slsqp_tol = 1.0e-10_wp
      real(wp) :: newton_tol = 1.0e-10_wp

      !> Maximum iterations
      integer :: multistart_max_iter = 100
      integer :: slsqp_max_iter = 100
      integer :: newton_max_iter = 100

      !> Verbosity level (0=silent, 1=default, 2=user events, 3=diagnostics)
      integer :: verbosity = 1

      !> Enable full debug output (solver iterations, stability tables, etc.)
      logical :: debug = .false.

      !> Maximum allowed displacement from anchor during projection (Bohr)
      real(wp) :: max_displacement_threshold = 15.0_wp

      !*-------------- Multi-start SLSQP sampling configuration-------------- *!

      !> Sampling radii for the multi-tangent seed rings (fraction of displacement bound)
      real(wp) :: multistart_radius(3) = [0.25_wp, 0.75_wp, 4.0_wp]
      !> Number of tangent-ring points per radius (historical name kept for compatibility)
      integer :: multistart_leb_num(3) = [6, 6, 14]

      !> Minimum point-to-point separation for deduplicating projected candidates (Bohr)
      real(wp) :: branch_sep_cut = 1.0E-10_wp
      !> Maximum allowed rho difference from the closest candidate to keep branch alternatives (Bohr)
      real(wp) :: branch_rho_cut = 2.0_wp

      !*------------------- Deflation solver configuration------------------- *!

      !> Maximum number of roots enumerated by either deflation solver per anchor
      integer :: deflation_max_roots = 12
      !> Deflation exponent p in ||x - x*||^{-p}
      integer :: deflation_p_power = 2
      !> Additive shift alpha in M(x) = prod (||x - x*||^{-p} + alpha)
      real(wp) :: deflation_alpha = 1.0_wp
      !> L2 tolerance for considering two enumerated roots identical (Bohr)
      real(wp) :: deflation_root_tol = 1.0e-6_wp
      !> Multiplier on inner-solver tolerances (tol/toldx/toldf for SLSQP,
      !> tol/tolx for Newton). Deflation only needs to identify the basin
      !> of each root; tight convergence is left to the downstream Newton
      !> refinement step. 100 = two orders of magnitude looser than the
      !> regular multistart tolerance.
      real(wp) :: deflation_tol_relax = 100.0_wp

      !*--------------------- Riemannian escape settings--------------------- *!

      !> Maximum iterations
      integer :: riemann_max_iter = 150
      !> Maximum iterations if line-search is used
      integer :: riemann_max_ls_iter = 20
      !> Initial line-search step factor for non-saddle Riemannian updates
      real(wp) :: riemann_alpha_init = 2.0_wp
      !> Geometric reduction factor applied when line-search trial is rejected
      real(wp) :: riemann_alpha_reduce = 0.75_wp
      !> Armijo sufficient-decrease coefficient for accepting line-search steps
      real(wp) :: riemann_armijo_c = 1.0e-4_wp
      !> Fixed step factor used to move along negative-curvature direction when escaping saddles
      real(wp) :: riemann_escape_step = 0.5_wp
      !> Eigenvalue threshold used to classify tangent reduced-Hessian curvature as saddle-like
      real(wp) :: riemann_saddle_eig_tol = 1.0e-8_wp
      !> Retraction settings
      integer :: retract_max_iter = 50

      !*----------------------- Optimality verification---------------------- *!

      !> Enable saddle-point check
      logical :: verify_optimality = .true.
      !> Abort/error if there is no strict minimum found
      logical :: strict_minimum_required = .true.

      !*------------------------- Working variables------------------------- *!

      !> Primitive functions (owned by projector)
      type(moist_cavity_drop_objective_phi_type) :: phi
      class(moist_cavity_drop_lsf_type), allocatable :: lsf

      !> Per-cell candidate atom lists for SSD screening
      type(moist_cell_grid_type) :: mol_cell_grid

      !> SSD point cache: skip recomputation when evaluation point unchanged.
      !> Eliminates redundant compute_ssd calls within SLSQP (objective+constraint
      !> share x) and Newton (residual+Jacobian share z) iterations.
      real(wp) :: cached_ssd_point(3) = [huge(0.0_wp), huge(0.0_wp), huge(0.0_wp)]
      logical :: ssd_cache_valid = .false.

      !> Cache for SLSQP callbacks (reuse phi computation across objective/gradient/constraint)
      real(wp) :: cached_phi0 = 0.0_wp
      real(wp) :: cached_phi1_r(3)

      !> Cache for LSF computations
      real(wp) :: cached_lsf0 = 0.0_wp
      real(wp) :: cached_lsf1_r(3)
      real(wp) :: cached_lsf2_rr(3, 3)

      !> Workspace for tangent-space saddle check
      real(wp) :: work_tangent_t1(3)
      real(wp) :: work_tangent_t2(3)

   contains
      procedure :: init => projector_init
      procedure :: init_primitives => projector_init_primitives
      procedure :: compute_ssd => projector_compute_ssd
      procedure :: project_point => projector_project_point
      procedure, private :: run_multistart_solver => projector_run_multistart_solver
      procedure, private :: run_single_solver => projector_run_single_solver
      procedure, private :: refine_seeds => projector_refine_seeds
      procedure :: filter_candidates
      procedure :: refine_point => projector_refine_point
      procedure :: check_constrained_optimality
      procedure :: riemannian_newton_escape
      procedure :: retract_to_surface
      procedure :: destroy => projector_destroy
      !> Finalizer
      final :: finalize_projector
   end type drop_projector_type

contains

   !* ================================================================================= *!
   !*            Augmented-Lagrangian solver callbacks (Newton / KKT system)            *!
   !* ================================================================================= *!

   !> Residual function for augmented Lagrangian system (context-aware)
   !> F(x, lambda) = [ phi_alpha*(x - anchor) - lambda*grad_LSF(x) ]
   !>                [            -LSF(x)                           ]
   !> Phi is the pure-quadratic anchor objective: phi_alpha/2 * ||x - anchor||^2.
   !> Also eagerly computes and caches lsf2_rr (Hessian) for the subsequent
   !> Jacobian call at the same point: the z012 accumulation always computes
   !> hessZ internally, so requesting lsf2_rr adds only the trivial quotient
   !> rule (~20 FLOPs), while saving a full z012 recomputation in the Jacobian.
   subroutine projection_residual(z, f, context)
      real(wp), dimension(:), intent(in) :: z   ! [x1, x2, x3, lambda]
      real(wp), dimension(:), intent(out) :: f  ! [f1, f2, f3, f4]
      class(*), intent(in) :: context           ! projection context

      real(wp) :: x(3), lambda, phi_alpha
      type(projection_context_type) :: ctx

      ! Extract context
      select type (context)
      type is (projection_context_type)
         ctx = context
      class default
         return
      end select

      if (.not. associated(ctx%projector)) return

      ! Extract variables
      x(:) = z(1:3)
      lambda = z(4)
      phi_alpha = ctx%projector%phi%param%phi_alpha

      ! Compute SSD (cached: skips if point unchanged)
      call ctx%projector%compute_ssd(x)

      ! Compute screened LSF value, gradient, AND Hessian.
      ! The Hessian is cached for the subsequent Jacobian call at the same point.
      call ctx%lsf%f012_r_screened( &
         lsf0=ctx%projector%cached_lsf0, &
         lsf1_r=ctx%projector%cached_lsf1_r, &
         lsf2_rr=ctx%projector%cached_lsf2_rr)

      ! Lagrangian gradient: grad_phi - lambda * grad_LSF
      ! grad_phi = phi_alpha * (x - anchor) for pure-quadratic anchor term
      f(1:3) = phi_alpha*(x(:) - ctx%anchor(:)) - lambda*ctx%projector%cached_lsf1_r(:)

      ! Constraint: -LSF(x) = 0
      f(4) = -ctx%projector%cached_lsf0

   end subroutine projection_residual

   !> Jacobian function for augmented Lagrangian system (context-aware)
   !> J = [ phi_alpha*I - lambda*Hess_LSF   -grad_LSF ]
   !>     [         -grad_LSF^T                  0     ]
   !> All LSF derivatives are already cached from the preceding residual call
   !> at the same point (nlesolver always calls func then grad at identical x).
   !> Hess_phi = phi_alpha * I (constant, pure-quadratic anchor).
   subroutine projection_jacobian(z, jac, context)
      real(wp), dimension(:), intent(in) :: z       ! [x1, x2, x3, lambda]
      real(wp), dimension(:, :), intent(out) :: jac  ! (4,4)
      class(*), intent(in) :: context               ! projection context

      real(wp) :: lambda, phi_alpha
      integer :: i
      type(projection_context_type) :: ctx

      ! Extract context
      select type (context)
      type is (projection_context_type)
         ctx = context
      class default
         return
      end select

      if (.not. associated(ctx%projector)) return

      lambda = z(4)
      phi_alpha = ctx%projector%phi%param%phi_alpha

      ! All LSF derivatives (lsf0, lsf1_r, lsf2_rr) were cached by the
      ! preceding projection_residual call at the same point. No SSD
      ! computation or LSF accumulation needed - pure matrix assembly.

      ! Upper-left 3x3: Lagrangian Hessian = phi_alpha*I - lambda*Hess_LSF
      jac(1:3, 1:3) = -lambda*ctx%projector%cached_lsf2_rr(:, :)
      do i = 1, 3
         jac(i, i) = jac(i, i) + phi_alpha
      end do

      ! Upper-right 3x1: -grad_LSF
      jac(1:3, 4) = -ctx%projector%cached_lsf1_r(:)

      ! Lower-left 1x3: -grad_LSF^T
      jac(4, 1:3) = -ctx%projector%cached_lsf1_r(:)

      ! Lower-right 1x1: 0
      jac(4, 4) = 0.0_wp

   end subroutine projection_jacobian

   !> Build a warm-started seed z = (x_init, lambda_init) for the
   !> Newton-deflation 4-D KKT system by taking one tangent-plane Newton
   !> step toward the level-set surface from the anchor.
   !>
   !> At z = (anchor, 0) the residual reduces to F = (0, 0, 0, -L(anchor)),
   !> so the inner solver's first Newton step has magnitude
   !> ||dx|| = |L(anchor)| / ||grad L(anchor)||, which becomes huge near
   !> branched anchors that sit close to triple-junctions of L (where
   !> ||grad L|| -> 0). Pre-computing the same step here lets us *cap*
   !> its magnitude and pass a moderate, well-conditioned seed to the
   !> inner Newton, so deflation iterations start in a regime where each
   !> step is small.
   !>
   !> Formula (assuming x_in = anchor):
   !>   step = -L * grad_L / ||grad_L||^2
   !>   x_init = anchor + step                     (tangent-plane projection)
   !>   lambda_init = -phi_alpha * L / ||grad_L||^2 (KKT condition at x_init)
   !>
   !> The step is capped to 0.5 * branch_rho_cut so the warm-started seed
   !> stays well inside the post-first-root ball cap. If ||grad_L||^2 is
   !> below an absolute floor we fall back to (x_in, 0) and let the inner
   !> Newton try its own first step.
   !>
   !> @param[in]  x_in           Caller's spatial seed (typically the anchor)
   !> @param[in]  anchor         Anchor point (defines the phi objective)
   !> @param[in]  phi_alpha      Coefficient in phi(x) = phi_alpha/2 * ||x - anchor||^2
   !> @param[in]  branch_rho_cut Maximum rho excursion allowed for branches (Bohr)
   !> @param[in]  context        Projection context (carries projector pointer)
   !> @param[out] z_seed         Warm-started 4-D seed (x_init, lambda_init)
   subroutine newton_warm_start_seed(x_in, anchor, phi_alpha, branch_rho_cut, &
                                     context, z_seed)
      !> Spatial seed (3 components)
      real(wp), intent(in) :: x_in(3)
      !> Anchor point (3 components)
      real(wp), intent(in) :: anchor(3)
      !> phi objective coefficient
      real(wp), intent(in) :: phi_alpha
      !> Maximum rho excursion (Bohr); used to cap the warm-start step
      real(wp), intent(in) :: branch_rho_cut
      !> Projection context (forwarded to projection_residual / _jacobian)
      class(*), intent(in) :: context
      !> Warm-started seed for the 4-D KKT solve
      real(wp), intent(out) :: z_seed(4)

      real(wp) :: z0(4), f0(4), jac0(4, 4)
      real(wp) :: lsf_val, grad_lsf(3), grad_norm_sq, step(3), step_norm, scale
      real(wp), parameter :: grad_floor_sq = 1.0e-12_wp
      real(wp), parameter :: step_cap_frac = 0.5_wp

      ! Default fallback: zero-lambda seed at the caller's spatial point.
      z_seed(1:3) = x_in
      z_seed(4) = 0.0_wp

      ! Probe (anchor, lambda=0): residual gives -L(anchor), Jacobian gives
      ! -grad L(anchor) in column 4 (or row 4 - they're symmetric).
      z0(1:3) = anchor
      z0(4) = 0.0_wp
      call projection_residual(z0, f0, context)
      call projection_jacobian(z0, jac0, context)

      lsf_val = -f0(4)
      grad_lsf = -jac0(1:3, 4)
      grad_norm_sq = dot_product(grad_lsf, grad_lsf)

      ! If ||grad L|| is essentially zero we cannot form the tangent-plane
      ! projection. Hand the un-warmed seed to Newton and let line search
      ! do what it can.
      if (grad_norm_sq < grad_floor_sq) return

      step = -lsf_val*grad_lsf/grad_norm_sq
      step_norm = norm2(step)

      ! Cap the step magnitude so the warm-started seed stays well inside
      ! the post-first-root ball cap (|x - anchor| <= phi_min + rho_cut).
      scale = 1.0_wp
      if (step_norm > step_cap_frac*branch_rho_cut .and. step_norm > 0.0_wp) then
         scale = step_cap_frac*branch_rho_cut/step_norm
         step = scale*step
      end if

      z_seed(1:3) = anchor + step
      z_seed(4) = scale*(-phi_alpha*lsf_val/grad_norm_sq)
   end subroutine newton_warm_start_seed

   !* ================================================================================= *!
   !*                          Initialization and SSD screening                         *!
   !* ================================================================================= *!

   !> Initialize the projector with molecular structure and solver parameters
   subroutine projector_init(self, param, lsf_model, branch_sep_cut, branch_rho_cut, &
                             tol, maxiter, verbosity, debug)
      class(drop_projector_type), intent(inout) :: self
      type(moist_cavity_drop_parameters_type), intent(in) :: param
      !> Polymorphic LSF template; the projector source-allocates its own clone.
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_model
      !> Minimum point-to-point separation for deduplicating projected candidates (Bohr)
      real(wp), intent(in) :: branch_sep_cut
      !> Maximum allowed rho difference from the closest candidate (Bohr)
      real(wp), intent(in) :: branch_rho_cut
      real(wp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      integer, intent(in), optional :: verbosity
      logical, intent(in), optional :: debug

      ! Store branching cutoffs
      self%branch_sep_cut = branch_sep_cut
      self%branch_rho_cut = branch_rho_cut

      ! Store solver configuration
      if (present(tol)) then
         self%slsqp_tol = tol
         self%newton_tol = tol
         self%multistart_tol = 10.0_wp*tol
      end if
      if (present(maxiter)) then
         self%slsqp_max_iter = maxiter
         self%newton_max_iter = maxiter
         self%multistart_max_iter = maxiter
      end if
      if (present(verbosity)) self%verbosity = verbosity
      if (present(debug)) self%debug = debug

      ! Set parameters for primitives. The LSF is cloned from the caller's
      ! template; its constructor was already invoked at cavity setup.
      call self%phi%set_parameters(param)
      if (allocated(self%lsf)) deallocate (self%lsf)
      allocate (self%lsf, source=lsf_model)

   end subroutine projector_init
   !> Initialize the primitives

   subroutine projector_init_primitives(self, mol, radii, mol_cell_grid)
      class(drop_projector_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      type(moist_cell_grid_type), intent(in), optional :: mol_cell_grid

      ! Initialize primitives
      call self%phi%set_input(mol, radii)
      call self%lsf%update(mol, radii)

      ! Invalidate SSD point cache (geometry changed)
      self%ssd_cache_valid = .false.

      if (present(mol_cell_grid)) then
         self%mol_cell_grid = mol_cell_grid
      end if

   end subroutine projector_init_primitives

   !> Compute SSD data for an evaluation point using per-cell screening.
   !>
   !> Queries the molecular cell grid for the candidate atom list of the cell
   !> containing `point` (strict clamp for points outside the atom bounding
   !> box) and passes it to the SSD subset routine. Zero allocation on the
   !> hot path - `cell_nlat(start+1:start+n)` is a contiguous slice.
   !>
   !> @param[inout] self  Projector instance (cell grid must be built)
   !> @param[in]    point Evaluation point (3)
   subroutine projector_compute_ssd(self, point)
      class(drop_projector_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)

      integer :: start, n

      !> Skip if we already computed SSD for this exact point
      if (self%ssd_cache_valid) then
         if (point(1) == self%cached_ssd_point(1) .and. &
             point(2) == self%cached_ssd_point(2) .and. &
             point(3) == self%cached_ssd_point(3)) return
      end if

      call self%mol_cell_grid%query(point, start, n)
      call self%lsf%prepare_subset(point, &
                                   self%mol_cell_grid%cell_nlat(start + 1:start + n))

      self%cached_ssd_point(:) = point(:)
      self%ssd_cache_valid = .true.
   end subroutine projector_compute_ssd

   !* ================================================================================= *!
   !*            Solver callbacks (diagnostics, displacement control, SLSQP)            *!
   !* ================================================================================= *!

   !> Debug callback - shows physical quantities for projection (context-aware)
   !> Variables: x(1:3) = point coordinates, x(4) = lambda
   !> Residuals: f(1:3) = dL (Lagrangian gradient), f(4) = -G (LSF constraint)
   subroutine newton_debug_callback(iter, x, f, context, jac, jac_sparse)
      integer, intent(in) :: iter
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(in) :: f
      class(*), intent(in) :: context
      real(wp), dimension(:, :), intent(in), optional :: jac
      real(wp), dimension(:), intent(in), optional :: jac_sparse

      real(wp) :: point(3), lambda_val
      real(wp) :: S_val, phi_val, grad_S(3)
      real(wp) :: norm_grad_S, norm_grad_L, constraint_viol
      real(wp) :: delta_phi, rho
      type(projection_context_type) :: ctx

      ! Extract context
      select type (context)
      type is (projection_context_type)
         ctx = context
      class default
         return
      end select

      ! Extract solution variables: [x, y, z, lambda ]
      point(:) = x(1:3)
      lambda_val = x(4)

      ! Evaluate physical quantities
      if (associated(ctx%projector)) then
         ! Phi objective and its gradient
         call ctx%projector%compute_ssd(point)
         phi_val = ctx%projector%phi%f0(point, ctx%anchor, ctx%owner)

         ! LSF constraint G(x) and its gradient (cached from residual/Jacobian)
         S_val = ctx%projector%cached_lsf0
         grad_S(:) = ctx%projector%cached_lsf1_r(:)
         norm_grad_S = sqrt(sum(grad_S**2))

         ! Change in phi from anchor point
         call ctx%projector%compute_ssd(ctx%anchor)
         delta_phi = phi_val - ctx%projector%phi%f0(ctx%anchor, ctx%anchor, ctx%owner)

         ! Distance from anchor (rho)
         rho = sqrt(sum((point - ctx%anchor)**2))
      else
         S_val = 0.0_wp
         norm_grad_S = 0.0_wp
         delta_phi = 0.0_wp
         rho = 0.0_wp
      end if

      ! Lagrangian gradient norm || dL|| (optimality condition)
      norm_grad_L = sqrt(sum(f(1:3)**2))

      ! Constraint violation |G(x)|
      constraint_viol = abs(f(4))

      ! Print header once
      if (iter == 1) then
         write (output_unit, '(a6,1x,a14,1x,a14,1x,a14,1x,a14,1x,a14,1x,a14,1x,a14)') &
            'Iter', 'LSF', '||dLSF||', 'dPhi', 'rho', 'lambda', '||dL||', '|G|'
         write (output_unit, '(A6,1X,7(A14,1X))') "-----", "--------------", &
            "--------------", "--------------", "--------------", "--------------", &
            "--------------", "--------------"
      end if

      ! Print iteration data with physically meaningful quantities
      write (output_unit, '(i6,1x,7(e14.6,1x),3f12.4)') &
         iter, S_val, norm_grad_S, delta_phi, rho, lambda_val, norm_grad_L, constraint_viol

   end subroutine newton_debug_callback

   !> Check if displacement exceeds threshold (user_input_check callback, context-aware)
   !> Returns .true. to stop the solver if rho > threshold
   function displacement_check(x, context) result(stop_solver)
      real(wp), dimension(:), intent(in) :: x
      class(*), intent(in) :: context
      logical :: stop_solver
      real(wp) :: point(3), rho
      type(projection_context_type) :: ctx

      stop_solver = .false.

      ! Extract context
      select type (context)
      type is (projection_context_type)
         ctx = context
      class default
         return
      end select

      ! Only check if threshold is set (> 0)
      if (ctx%threshold <= 0.0_wp) return

      ! Extract point coordinates
      point(:) = x(1:3)

      ! Compute displacement from anchor
      rho = sqrt(sum((point - ctx%anchor)**2))

      ! Stop if threshold exceeded
      if (rho > ctx%threshold) then
         stop_solver = .true.
         if (associated(ctx%projector)) then
            if (ctx%projector%debug) then
               write (output_unit, '(a,f12.6,a,f12.6,a)') &
                  '[displacement_check] Stopping: rho = ', rho, ' Bohr > threshold = ', ctx%threshold, ' Bohr'
            end if
         end if
      end if

   end function displacement_check
   !> SLSQP iteration callback for debugging (context-aware)

   subroutine slsqp_debug_callback(iter, x, f, c, context)
      integer, intent(in) :: iter
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(in) :: f
      real(wp), dimension(:), intent(in) :: c
      class(*), intent(in) :: context

      real(wp) :: S_val, grad_S(3), norm_grad_S
      real(wp) :: grad_phi(3), delta_phi, rho
      type(projection_context_type) :: ctx

      ! Extract context
      select type (context)
      type is (projection_context_type)
         ctx = context
      class default
         return
      end select

      if (.not. associated(ctx%projector)) return

      ! LSF constraint G(x) and its gradient (using screened)
      ! Note: slsqp_constraint/grad called before this, so values should be in cache
      S_val = ctx%projector%cached_lsf0
      grad_S(:) = ctx%projector%cached_lsf1_r(:)
      norm_grad_S = sqrt(sum(grad_S**2))

      ! Phi objective gradient
      call ctx%projector%compute_ssd(x)
      grad_phi(:) = ctx%projector%phi%f1_r(x, ctx%anchor, ctx%owner)

      ! Change in phi from anchor
      delta_phi = f  ! SLSQP gives us objective value directly

      ! Distance from anchor
      rho = sqrt(sum((x - ctx%anchor)**2))

      ! Print header on first iteration
      if (iter == 0) then
         write (output_unit, '(a6,1x,a14,1x,a15,1x,a16,1x,a16,1x,a16,1x,a14)') &
            'Iter', 'LSF', '||∇LSF||', 'ΔΦ (obj)', 'ρ', '||∇Φ||', '|G| (con)'
         write (output_unit, '(A6,1X,6(A14,1X))') "-----", "--------------", &
            "--------------", "--------------", "--------------", "--------------", &
            "--------------"
      end if

      ! Print iteration data
      write (output_unit, '(i6,1x,6(e14.6,1x))') &
         iter, S_val, norm_grad_S, delta_phi, rho, sqrt(sum(grad_phi**2)), abs(c(1))

   end subroutine slsqp_debug_callback

   !*---------------- SLSQP objective and constraint callbacks--------------- *!

   !> SLSQP objective function: Phi (x) = full phi objective (context-aware)
   !> Computes phi value and gradient, caches them for subsequent gradient/constraint calls
   subroutine slsqp_objective(x, f, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      class(*), intent(in) :: context
      type(projection_context_type) :: ctx

      ! Extract context
      select type (context)
      type is (projection_context_type)
         ctx = context
      class default
         return
      end select

      if (.not. associated(ctx%projector)) return

      ! Phi is pure-quadratic (anchor term only) and does not use SSD data,
      ! so compute_ssd is deferred to the constraint callback where it is needed.
      call ctx%projector%phi%f012_r(x, ctx%anchor, ctx%owner, &
         & ctx%projector%cached_phi0, ctx%projector%cached_phi1_r)

      f = ctx%projector%cached_phi0

   end subroutine slsqp_objective

   !> SLSQP objective gradient: dPhi (x) = full phi gradient (context-aware)
   !> Uses cached gradient computed by slsqp_objective
   subroutine slsqp_objective_grad(x, grad_f, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: grad_f
      class(*), intent(in) :: context
      type(projection_context_type) :: ctx

      ! Extract context
      select type (context)
      type is (projection_context_type)
         ctx = context
      class default
         return
      end select

      if (.not. associated(ctx%projector)) return

      ! Return cached gradient from slsqp_objective
      grad_f = ctx%projector%cached_phi1_r

   end subroutine slsqp_objective_grad

   !> SLSQP constraint: g(x) = LSF(x) = 0 (context-aware)
   !> Computes and caches screened LSF value
   subroutine slsqp_constraint(x, g, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: g
      class(*), intent(in) :: context
      type(projection_context_type) :: ctx

      ! Extract context
      select type (context)
      type is (projection_context_type)
         ctx = context
      class default
         return
      end select

      if (.not. associated(ctx%projector)) return

      ! Compute screened LSF value (and grad) and cache them in projector
      ! Using ctx%lsf which comes from the projector context
      call ctx%projector%compute_ssd(x)

      call ctx%lsf%f012_r_screened( &
         lsf0=ctx%projector%cached_lsf0, &
         lsf1_r=ctx%projector%cached_lsf1_r)

      g(1) = ctx%projector%cached_lsf0

   end subroutine slsqp_constraint
   !> SLSQP constraint gradient: dg(x) = dLSF(x) (context-aware)
   !> Uses cached LSF gradient computed by slsqp_constraint

   subroutine slsqp_constraint_grad(x, grad_g, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:, :), intent(out) :: grad_g
      class(*), intent(in) :: context
      type(projection_context_type) :: ctx

      ! Extract context
      select type (context)
      type is (projection_context_type)
         ctx = context
      class default
         return
      end select

      if (.not. associated(ctx%projector)) return

      ! Return cached LSF gradient from projector's cache.
      ! Assumes slsqp_constraint was called for the same x before this
      ! (standard SLSQP order: constraints evaluated, then gradients).

      grad_g(1, :) = ctx%projector%cached_lsf1_r(:)

   end subroutine slsqp_constraint_grad

   !* ================================================================================= *!
   !*                    Seed generation: multi-start solver drivers                    *!
   !* ================================================================================= *!

   !> Extract candidate seeds from a multi-start or multi-tangent solver.
   !> Unifies the type-dispatch + filter + packing logic shared by both solver types.
   !> @param[inout] self       Projector (needed for filter_candidates)
   !> @param[inout] solver     Solved multi-start or multi-tangent solver
   !> @param[in]    anchor     Anchor point for filtering
   !> @param[out]   x_seeds    Filtered seed points (3, n_seeds)
   !> @param[out]   n_seeds    Number of retained seeds
   !> @param[out]   error      Error descriptor
   subroutine extract_seeds_from_solver(self, solver, anchor, x_seeds, n_seeds, error)
      class(drop_projector_type), intent(inout) :: self
      class(solver_base_type), intent(inout) :: solver
      real(wp), intent(in) :: anchor(3)
      real(wp), allocatable, intent(out) :: x_seeds(:, :)
      integer, intent(out) :: n_seeds
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: candidates(:, :)
      logical, allocatable :: keep_mask(:)
      integer :: n_candidates, n_unique, i

      n_seeds = 0

      ! Dispatch to the type-specific get_raw_candidates call. Newton-deflation
      ! returns 4-row (x, lambda) candidates; we strip lambda right here so the
      ! downstream filter only sees xyz.
      select type (ms_solver => solver)
      type is (moist_math_solver_slsqp_multistart_type)
         call ms_solver%get_raw_candidates(candidates, n_candidates)
      type is (moist_math_solver_slsqp_multi_tangent_type)
         call ms_solver%get_raw_candidates(candidates, n_candidates)
      type is (moist_math_solver_slsqp_deflation_type)
         call ms_solver%get_raw_candidates(candidates, n_candidates)
      type is (moist_math_solver_newton_deflation_type)
         block
            real(wp), allocatable :: raw4(:, :)
            integer :: ic
            call ms_solver%get_raw_candidates(raw4, n_candidates)
            if (n_candidates > 0) then
               allocate (candidates(3, n_candidates))
               do ic = 1, n_candidates
                  candidates(:, ic) = raw4(1:3, ic)
               end do
            else
               allocate (candidates(3, 0))
            end if
            if (allocated(raw4)) deallocate (raw4)
         end block
      class default
         call fatal_error(error, "Unexpected solver type in multi-start projection")
         return
      end select

      ! Validate candidate count
      if (n_candidates <= 0) then
         call fatal_error(error, "Multi-start solver returned no candidates")
         if (allocated(candidates)) deallocate (candidates)
         return
      end if

      ! Filter duplicates and outliers
      call self%filter_candidates(candidates, n_candidates, keep_mask, error, anchor)
      if (allocated(error)) then
         if (allocated(candidates)) deallocate (candidates)
         if (allocated(keep_mask)) deallocate (keep_mask)
         return
      end if

      n_seeds = count(keep_mask)
      if (n_seeds <= 0) then
         call fatal_error(error, "No unique candidates retained after filtering")
         if (allocated(candidates)) deallocate (candidates)
         if (allocated(keep_mask)) deallocate (keep_mask)
         return
      end if

      ! Pack retained candidates into x_seeds
      allocate (x_seeds(3, n_seeds), source=0.0_wp)
      n_unique = 0
      do i = 1, n_candidates
         if (keep_mask(i)) then
            n_unique = n_unique + 1
            x_seeds(:, n_unique) = candidates(:, i)
         end if
      end do

      if (allocated(candidates)) deallocate (candidates)
      if (allocated(keep_mask)) deallocate (keep_mask)

   end subroutine extract_seeds_from_solver

   !> Run a multi-start or multi-tangent SLSQP solver and return filtered seed points.
   !> @param[inout] self         Projector instance
   !> @param[in]    level        Projection level (8 fine reference, 7 multistart, 6 Newton-defl, 4-5 SLSQP-defl, <4 multi-tangent)
   !> @param[in]    anchor       Anchor point
   !> @param[in]    xl           Lower bounds
   !> @param[in]    xu           Upper bounds
   !> @param[in]    proj_context Callback context
   !> @param[in]    owner        Atom owner index
   !> @param[in]    index        Grid point index (for diagnostics)
   !> @param[inout] x_slsqp     Initial guess (overwritten by solver)
   !> @param[out]   x_seeds      Filtered seed points (3, n_seeds)
   !> @param[out]   n_seeds      Number of retained seeds
   !> @param[in]    retry_radius Adaptive perturbation radius for SLSQP-deflation retries
   !> @param[out]   error        Error descriptor
   subroutine projector_run_multistart_solver(self, level, anchor, xl, xu, proj_context, &
                                              owner, index, x_slsqp, x_seeds, n_seeds, &
                                              retry_radius, error)
      class(drop_projector_type), intent(inout), target :: self
      integer, intent(in) :: level
      real(wp), intent(in) :: anchor(3), xl(3), xu(3)
      type(projection_context_type), intent(in) :: proj_context
      integer, intent(in) :: owner, index
      real(wp), intent(inout) :: x_slsqp(3)
      real(wp), allocatable, intent(out) :: x_seeds(:, :)
      integer, intent(out) :: n_seeds
      real(wp), intent(in) :: retry_radius
      type(error_type), allocatable, intent(out) :: error

      class(solver_base_type), allocatable :: solver
      real(wp) :: z_seed(4), lxl(4), lxu(4)
      real(wp), allocatable :: fine_radii(:)
      integer, allocatable :: fine_n_points(:)
      integer :: req_max_deriv

      n_seeds = 0

      ! Choose solver and assemble it. The SLSQP variants work in 3-D xyz;
      ! Newton-deflation works on the full 4-D KKT system z = (x, lambda),
      ! so its branch needs separate seed handling.
      !
      !   level == 8 : fine SLSQP multistart reference                    -- 3D
      !   level == 7 : regular SLSQP multistart baseline                  -- 3D
      !   level == 6 : Newton-deflation on the 4D KKT system               -- 4D
      !   level 4-5  : SLSQP-deflation (iterated Farrell deflation)        -- 3D
      !   level <  4 : multi-tangent (cheap local rings)                  -- 3D
      if (level >= 8) then
         call build_fine_multistart_profile(fine_radii, fine_n_points)
         call new_slsqp_multistart_solver( &
            anchor=anchor, &
            solver=solver, &
            n=3, &
            m=1, &
            meq=1, &
            obj_ctx=slsqp_objective, &
            obj_grad_ctx=slsqp_objective_grad, &
            con_ctx=slsqp_constraint, &
            con_grad_ctx=slsqp_constraint_grad, &
            context=proj_context, &
            xl=xl, &
            xu=xu, &
            max_iter=self%multistart_max_iter, &
            tol=self%multistart_tol, &
            toldx=self%multistart_tol, &
            toldf=self%multistart_tol, &
            radii=fine_radii, &
            n_points=fine_n_points, &
            normal=anchor - self%lsf%mol%xyz(:, owner), &
            error=error &
            )
      else if (level >= 7) then
         call new_slsqp_multistart_solver( &
            anchor=anchor, &
            solver=solver, &
            n=3, &
            m=1, &
            meq=1, &
            obj_ctx=slsqp_objective, &
            obj_grad_ctx=slsqp_objective_grad, &
            con_ctx=slsqp_constraint, &
            con_grad_ctx=slsqp_constraint_grad, &
            context=proj_context, &
            xl=xl, &
            xu=xu, &
            max_iter=self%multistart_max_iter, &
            tol=self%multistart_tol, &
            toldx=self%multistart_tol, &
            toldf=self%multistart_tol, &
            radii=self%multistart_radius, &
            n_points=self%multistart_leb_num, &
            normal=anchor - self%lsf%mol%xyz(:, owner), &
            error=error &
            )
      else if (level == 6) then
         ! Newton solves the augmented Lagrangian system, so we need
         ! second-order LSF derivatives for the Hessian block.
         lxl(1:3) = xl
         lxu(1:3) = xu
         lxl(4) = -1.0e6_wp
         lxu(4) = 1.0e6_wp
         call new_newton_deflation_solver( &
            solver=solver, &
            n=4, m=4, &
            func_ctx=projection_residual, &
            grad_ctx=projection_jacobian, &
            context=proj_context, &
            max_iter=self%newton_max_iter, &
            tol=self%newton_tol, &
            tolx=self%newton_tol*0.1_wp, &
            tol_relax_factor=self%deflation_tol_relax, &
            alpha=0.01_wp, &
            use_broyden=.false., &
            max_roots=self%deflation_max_roots, &
            p_power=self%deflation_p_power, &
            alpha_shift=self%deflation_alpha, &
            dedup_tol=self%deflation_root_tol, &
            bounds_mode=3, &
            xlow=lxl, xupp=lxu, &
            anchor=anchor, &
            branch_rho_cut=self%branch_rho_cut, &
            error=error &
            )
      else if (level >= 4) then
         call new_slsqp_deflation_solver( &
            solver=solver, &
            n=3, &
            m=1, &
            meq=1, &
            obj_ctx=slsqp_objective, &
            obj_grad_ctx=slsqp_objective_grad, &
            con_ctx=slsqp_constraint, &
            con_grad_ctx=slsqp_constraint_grad, &
            context=proj_context, &
            xl=xl, &
            xu=xu, &
            max_iter=self%multistart_max_iter, &
            tol=self%multistart_tol, &
            toldx=self%multistart_tol, &
            toldf=self%multistart_tol, &
            tol_relax_factor=self%deflation_tol_relax, &
            max_roots=self%deflation_max_roots, &
            p_power=self%deflation_p_power, &
            alpha_shift=self%deflation_alpha, &
            dedup_tol=self%deflation_root_tol, &
            retry_radius=retry_radius, &
            anchor=anchor, &
            branch_rho_cut=self%branch_rho_cut, &
            error=error &
            )
      else
         call new_slsqp_multi_tangent_solver( &
            anchor=anchor, &
            solver=solver, &
            n=3, &
            m=1, &
            meq=1, &
            obj_ctx=slsqp_objective, &
            obj_grad_ctx=slsqp_objective_grad, &
            con_ctx=slsqp_constraint, &
            con_grad_ctx=slsqp_constraint_grad, &
            context=proj_context, &
            xl=xl, &
            xu=xu, &
            owner_pos=self%lsf%mol%xyz(:, owner), &
            max_iter=self%multistart_max_iter, &
            tol=self%multistart_tol, &
            toldx=self%multistart_tol, &
            toldf=self%multistart_tol, &
            radii=self%multistart_radius, &
            n_points=self%multistart_leb_num, &
            error=error &
            )
      end if

      if (allocated(error)) then
         call fatal_error(error, "Failed to initialize multi-start solver for projection")
         return
      end if

      ! SLSQP-style solvers only need value+gradient (max_deriv=1, skips
      ! per-atom Hessians); Newton-deflation needs the LSF Hessian for the
      ! Lagrangian block, so it requires max_deriv=2.
      if (level == 6) then
         req_max_deriv = 2
      else
         req_max_deriv = 1
      end if
      call self%lsf%set_max_deriv(req_max_deriv)
      self%ssd_cache_valid = .false.

      ! Solve and extract candidates. Newton-deflation operates on z=(x,lambda),
      ! so it gets its own 4-D seed; SLSQP variants share x_slsqp.
      if (level == 6) then
         ! Warm-start the Newton seed with one tangent-plane step toward
         ! LSF=0 (and the corresponding lambda). Seeding at (anchor, 0)
         ! makes the first deflated Newton step's amplitude scale like
         ! L(anchor)/||grad L(anchor)||, which blows up at branched
         ! anchors near triple-junctions of LSF where the gradient is
         ! nearly stationary. A bounded warm-start sidesteps this.
         call newton_warm_start_seed(x_slsqp, anchor, &
                                     self%phi%param%phi_alpha, &
                                     self%branch_rho_cut, &
                                     proj_context, z_seed)
         call solver%solve(z_seed, error)
         x_slsqp = z_seed(1:3)
      else
         call solver%solve(x_slsqp, error)
      end if

      ! Restore second-order for subsequent Newton refinement
      call self%lsf%set_max_deriv(2)
      self%ssd_cache_valid = .false.

      if (allocated(error)) then
         if (self%verbosity >= 3) then
            write (output_unit, '(x,i12,x,a,a)') index, &
               '[multi-SLSQP] Failed: ', trim(error%message)
         end if
         call solver%destroy()
         deallocate (solver)
         return
      end if

      call extract_seeds_from_solver(self, solver, anchor, x_seeds, n_seeds, error)

      call solver%destroy()
      deallocate (solver)

   end subroutine projector_run_multistart_solver

   !> Build the deterministic high-density SLSQP multistart profile used by
   !> projection level 8 reference tests.
   subroutine build_fine_multistart_profile(radii, n_points)
      real(wp), allocatable, intent(out) :: radii(:)
      integer, allocatable, intent(out) :: n_points(:)

      integer, parameter :: n_shells = 20
      integer :: i
      real(wp) :: radius

      allocate (radii(n_shells), n_points(n_shells))
      do i = 1, n_shells
         radius = 0.10_wp + 0.35_wp*real(i - 1, wp)
         radii(i) = radius
         if (radius <= 0.70_wp) then
            n_points(i) = 50
         else if (radius <= 1.50_wp) then
            n_points(i) = 86
         else if (radius <= 2.50_wp) then
            n_points(i) = 110
         else if (radius <= 3.50_wp) then
            n_points(i) = 146
         else
            n_points(i) = 170
         end if
      end do
   end subroutine build_fine_multistart_profile

   !> Run a single local SLSQP solver and return its result as a one-element seed array.
   !> @param[inout] self         Projector instance
   !> @param[in]    anchor       Anchor point
   !> @param[in]    xl           Lower bounds
   !> @param[in]    xu           Upper bounds
   !> @param[in]    proj_context Callback context
   !> @param[in]    index        Grid point index (for diagnostics)
   !> @param[inout] x_slsqp     Initial guess (overwritten by solver)
   !> @param[out]   x_seeds      Seed point (3, 1)
   !> @param[out]   n_seeds      Always 1 on success
   !> @param[out]   error        Error descriptor
   subroutine projector_run_single_solver(self, anchor, xl, xu, proj_context, &
                                          index, x_slsqp, x_seeds, n_seeds, error)
      class(drop_projector_type), intent(inout), target :: self
      real(wp), intent(in) :: anchor(3), xl(3), xu(3)
      type(projection_context_type), intent(in) :: proj_context
      integer, intent(in) :: index
      real(wp), intent(inout) :: x_slsqp(3)
      real(wp), allocatable, intent(out) :: x_seeds(:, :)
      integer, intent(out) :: n_seeds
      type(error_type), allocatable, intent(out) :: error

      class(solver_base_type), allocatable :: solver

      call new_slsqp_solver( &
         solver=solver, &
         n=3, &
         m=1, &
         meq=1, &
         error=error, &
         obj_ctx=slsqp_objective, &
         obj_grad_ctx=slsqp_objective_grad, &
         con_ctx=slsqp_constraint, &
         con_grad_ctx=slsqp_constraint_grad, &
         context=proj_context, &
         max_iter=self%slsqp_max_iter, &
         tol=self%slsqp_tol, &
         toldx=self%slsqp_tol, &
         toldf=self%slsqp_tol, &
         xl=xl, &
         xu=xu &
         )

      if (allocated(error)) then
         call fatal_error(error, "Failed to initialize SLSQP solver for projection")
         return
      end if

      ! SLSQP only needs value + gradient (max_deriv=1 skips per-atom Hessians)
      call self%lsf%set_max_deriv(1)
      self%ssd_cache_valid = .false.

      call solver%solve(x_slsqp, error)

      ! Restore second-order for subsequent Newton refinement
      call self%lsf%set_max_deriv(2)
      self%ssd_cache_valid = .false.

      if (allocated(error)) then
         if (self%verbosity >= 3) then
            write (output_unit, '(x,i12,x,a,a)') index, &
               '[single-SLSQP] Failed: ', trim(error%message)
            write (output_unit, '(x,i12,x,a)') index, &
               '[single-SLSQP] Resorting to Newton ...'
         end if
         deallocate (error)
      end if

      call solver%destroy()
      deallocate (solver)

      n_seeds = 1
      allocate (x_seeds(3, 1), source=0.0_wp)
      x_seeds(:, 1) = x_slsqp(:)

   end subroutine projector_run_single_solver

   !* ================================================================================= *!
   !*                 Projection pipeline (seed refinement and dispatch)                *!
   !* ================================================================================= *!

   !> Refine seed points onto the LSF surface and populate the projection workspace.
   !> At level 1, uses fast SLSQP-only diagnostics; at level >= 2, applies Newton refinement.
   !> @param[inout] self      Projector instance
   !> @param[in]    level     Projection level
   !> @param[in]    anchor    Anchor point
   !> @param[in]    owner     Atom owner index
   !> @param[in]    index     Grid point index
   !> @param[in]    xl        Lower bounds
   !> @param[in]    xu        Upper bounds
   !> @param[in]    rho_max   Maximum displacement threshold
   !> @param[in]    x_seeds   Seed points (3, n_seeds)
   !> @param[in]    n_seeds   Number of seeds
   !> @param[inout] work      Projection workspace
   !> @param[inout] n_points  Running count of accepted points
   !> @param[out]   error     Error descriptor
   subroutine projector_refine_seeds(self, level, anchor, owner, index, xl, xu, &
                                     rho_max, x_seeds, n_seeds, work, n_points, error)
      class(drop_projector_type), intent(inout), target :: self
      integer, intent(in) :: level, owner, index, n_seeds
      real(wp), intent(in) :: anchor(3), xl(3), xu(3), rho_max
      real(wp), intent(in) :: x_seeds(:, :)
      type(projection_workspace_type), intent(inout) :: work
      integer, intent(inout) :: n_points
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: x_refined(3), lambda_refined, rho_refined, phi_refined, normal_refined(3)
      real(wp) :: grad_phi_seed(3), norm_grad2
      integer :: i_seed, n_failed
      character(len=512) :: seed_error_msgs(n_seeds)
      integer :: failed_indices(n_seeds)

      !> When n_seeds == 0 the caller (project_point, solver-failure path)
      !> intentionally passes an empty seed list so the n_points == 0 onion
      !> fallback can run. reserve(0) is a legal no-op in that case and the
      !> workspace may still be unallocated on its first use, so skip the
      !> allocation sanity check.
      if (n_seeds > 0) then
         call work%reserve(n_seeds)
         if (.not. allocated(work%points) .or. .not. allocated(work%normals) .or. &
             .not. allocated(work%rho) .or. .not. allocated(work%lambda) .or. &
             .not. allocated(work%phi) .or. .not. allocated(work%converged) .or. &
             .not. allocated(work%branch_weights)) then
            call fatal_error(error, "Projection workspace allocation failure after reserve")
            return
         end if
         if (size(work%points, dim=1) /= 3 .or. size(work%normals, dim=1) /= 3) then
            call fatal_error(error, "Projection workspace leading dimension corruption detected")
            return
         end if
      end if

      if (level == 1) then
         ! Fast path: keep pure SLSQP result, no Newton/KKT refinement
         ! Still populate diagnostics (rho, phi, normal, lambda estimate)
         do i_seed = 1, n_seeds
            n_points = n_points + 1
            if (n_points > work%capacity) call work%reserve(n_points)
            work%points(:, n_points) = x_seeds(:, i_seed)
            work%rho(n_points) = sqrt(dot_product( &
                                      x_seeds(:, i_seed) - anchor, x_seeds(:, i_seed) - anchor))

            call self%compute_ssd(x_seeds(:, i_seed))
            work%phi(n_points) = self%phi%f0( &
                                 x_seeds(:, i_seed), anchor, owner)
            grad_phi_seed = self%phi%f1_r( &
                            x_seeds(:, i_seed), anchor, owner)

            call self%lsf%f012_r_screened( &
               lsf0=self%cached_lsf0, &
               lsf1_r=self%cached_lsf1_r)
            norm_grad2 = dot_product(self%cached_lsf1_r, self%cached_lsf1_r)
            if (norm_grad2 > 1.0e-14_wp) then
               work%normals(:, n_points) = self%cached_lsf1_r/sqrt(norm_grad2)
               ! KKT condition: grad_phi = lambda * grad_S
               work%lambda(n_points) = dot_product(grad_phi_seed, self%cached_lsf1_r)/norm_grad2
            else
               work%normals(:, n_points) = 0.0_wp
               work%lambda(n_points) = 0.0_wp
            end if
            work%converged(n_points) = .true.
         end do
      else
         ! Accurate path: refine each seed with coupled Newton solve
         n_failed = 0
         do i_seed = 1, n_seeds
            call self%refine_point( &
               anchor=anchor, &
               x_slsqp=x_seeds(:, i_seed), &
               xl=xl, &
               xu=xu, &
               displacement_threshold=rho_max, &
               owner=owner, &
               index=index, &
               error=error, &
               gridpoint=x_refined, &
               lambda_out=lambda_refined, &
               rho_out=rho_refined, &
               phi_out=phi_refined, &
               normal_out=normal_refined &
               )

            if (allocated(error)) then
               n_failed = n_failed + 1
               failed_indices(n_failed) = i_seed
               seed_error_msgs(n_failed) = error%message
               if (self%verbosity >= 3) then
                  write (output_unit, '(x,i12,x,a,i0,a,a)') index, &
                     '[refinement] Dropping branch ', i_seed, &
                     ': ', trim(error%message)
               end if
               deallocate (error)
               cycle
            end if
            if (.not. allocated(work%points) .or. .not. allocated(work%normals)) then
               call fatal_error(error, "Projection workspace deallocated during refinement")
               return
            end if
            if (size(work%points, dim=1) /= 3 .or. size(work%normals, dim=1) /= 3) then
               call fatal_error(error, "Projection workspace leading dimension corrupted during refinement")
               return
            end if

            n_points = n_points + 1
            if (n_points > work%capacity) call work%reserve(n_points)
            work%points(:, n_points) = x_refined(:)
            work%lambda(n_points) = lambda_refined
            work%rho(n_points) = rho_refined
            work%phi(n_points) = phi_refined
            work%normals(:, n_points) = normal_refined(:)
            work%converged(n_points) = .true.
         end do

         ! Print collected errors when all seeds failed (no error propagation
         ! here - let the caller handle fallback via n_points == 0)
         if (n_points == 0 .and. n_failed > 0) then
            if (self%verbosity >= 3) then
               write (output_unit, '(x,i12,x,a,i0,a,i0,a)') &
                  index, '[refinement] All ', n_failed, ' of ', n_seeds, &
                  ' seed(s) failed:'
               do i_seed = 1, n_failed
                  write (output_unit, '(x,i12,x,a,i0,a,a)') &
                     index, '[refinement]   seed ', failed_indices(i_seed), ': ', &
                     trim(seed_error_msgs(i_seed))
               end do
            end if
         end if
      end if

   end subroutine projector_refine_seeds

   !*-------------------- Top-level projection entry point------------------- *!

   !> Project an anchor point onto the LSF surface and return one or more branches.
   !>
   !> Dispatches to a single-SLSQP or multi-start solver depending on
   !> projection level and local gradient quality, refines the resulting seeds,
   !> deduplicates, and populates the workspace and optional output arrays.
   !> If the single-solver path yields no surviving points, a multistart
   !> fallback is attempted before giving up.
   subroutine projector_project_point(self, anchor, gridpoints, n_points, owner, index, &
                                      proj_level, &
                                      initial_guess, error, &
                                      lambda_out, rho_out, normal_out, phi_out, &
                                      max_displacement, work)
      class(drop_projector_type), intent(inout), target :: self
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in), optional :: initial_guess(3)
      integer, intent(in), optional :: proj_level
      real(wp), allocatable, intent(out), optional :: gridpoints(:, :)
      integer, intent(out) :: n_points
      integer, intent(in) :: owner, index
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable, intent(out), optional :: lambda_out(:)
      real(wp), allocatable, intent(out), optional :: normal_out(:, :)
      real(wp), allocatable, intent(out), optional :: rho_out(:)
      real(wp), allocatable, intent(out), optional :: phi_out(:)
      real(wp), intent(in), optional :: max_displacement
      type(projection_workspace_type), intent(inout) :: work

      real(wp) :: S_grad(3), tmp_val(3), tmp_grad(3, 3)
      real(wp) :: x_slsqp(3), xl(3), xu(3), rho_max, adaptive_retry_rad
      real(wp), allocatable :: x_seeds(:, :)
      logical, allocatable :: keep_mask(:)
      integer :: n_seeds, i_keep, n_unique, level
      logical :: use_multistart
      type(projection_context_type) :: proj_context

      !> Setup

      n_points = 0
      call work%clear()
      self%ssd_cache_valid = .false.

      level = self%proj_level
      if (present(proj_level)) level = proj_level

      ! Initial guess
      if (present(initial_guess)) then
         x_slsqp(:) = initial_guess(:)
      else
         x_slsqp(:) = anchor(:)
      end if

      ! Callback context
      proj_context%projector => self
      proj_context%lsf => self%lsf
      proj_context%anchor(:) = anchor(:)
      proj_context%owner = owner

      ! Displacement bounds
      if (present(max_displacement)) then
         rho_max = max_displacement
      else
         rho_max = self%max_displacement_threshold
      end if
      proj_context%threshold = rho_max
      xl(:) = anchor(:) - rho_max
      xu(:) = anchor(:) + rho_max

      !> Evaluate LSF gradient at anchor

      call slsqp_constraint(anchor, tmp_val, proj_context)
      call slsqp_constraint_grad(anchor, tmp_grad, proj_context)
      S_grad(:) = tmp_grad(1, :)

      !> Solver dispatch

      ! Adaptive retry radius: mean of |S| and the regularized Newton step S/sqrt(|grad S|^2 + 0.3^2)
      adaptive_retry_rad = 0.5_wp*( &
                           abs(tmp_val(1)) + &
                           abs(tmp_val(1))/sqrt(dot_product(S_grad, S_grad) + 0.09_wp) &
                           )
      adaptive_retry_rad = max(adaptive_retry_rad, 0.05_wp)

      use_multistart = ( &
                       ((level >= 3) .and. ((norm2(S_grad) < 0.25_wp))) .or. &
                       ((level >= 4) .and. ((norm2(S_grad) < 0.60_wp))) .or. &
                       ((level >= 5)) &
                       )

      if (use_multistart) then
         call self%run_multistart_solver(level, anchor, xl, xu, proj_context, &
                                         owner, index, x_slsqp, x_seeds, n_seeds, &
                                         adaptive_retry_rad, error)
      else
         call self%run_single_solver(anchor, xl, xu, proj_context, &
                                     index, x_slsqp, x_seeds, n_seeds, error)
      end if

      ! A solver error is treated like "0 candidates"
      if (allocated(error)) then
         if (level >= 7) return
         if (self%verbosity >= 3) then
            write (output_unit, '(x,i12,x,a,a)') index, &
               '[solver-fail] Falling back to multistart: ', trim(error%message)
         end if
         deallocate (error)
         n_seeds = 0
         if (allocated(x_seeds)) deallocate (x_seeds)
         allocate (x_seeds(3, 0))
      end if

      !> Refine seeds

      call self%refine_seeds(level, anchor, owner, index, xl, xu, &
                             rho_max, x_seeds, n_seeds, work, n_points, error)
      if (allocated(error)) return

      !> Fallback: retry with regular multistart (level 7) if the chosen solver produced 0 points

      if (n_points == 0 .and. level < 7) then
         if (self%verbosity >= 3) then
            write (output_unit, '(x,i12,x,a)') &
               index, '[fallback] Previous run produced 0 points, retrying with multistart'
         end if

         if (allocated(x_seeds)) deallocate (x_seeds)
         call self%run_multistart_solver(7, anchor, xl, xu, proj_context, &
                                         owner, index, x_slsqp, x_seeds, n_seeds, &
                                         adaptive_retry_rad, error)
         if (allocated(error)) then
            deallocate (error)
            call fatal_error(error, &
                             "Projection failed for all candidate branches (including fallback)")
            return
         end if

         call work%clear()
         n_points = 0
         call self%refine_seeds(level, anchor, owner, index, xl, xu, &
                                rho_max, x_seeds, n_seeds, work, n_points, error)
         if (allocated(error)) return
         if (self%verbosity >= 3) then
            write (output_unit, '(1x,i12,x,a,i0,a)') &
               index, "[fallback] Multistart produced ", n_points, " candidate points"
         end if
      end if

      !> Final deduplication

      if (n_points > 1) then
         call self%filter_candidates(work%points(:, 1:n_points), n_points, &
                                     keep_mask, error, anchor)
         if (allocated(error)) return
         n_unique = 0
         do i_keep = 1, n_points
            if (keep_mask(i_keep)) then
               n_unique = n_unique + 1
               if (n_unique /= i_keep) then
                  work%points(:, n_unique) = work%points(:, i_keep)
                  work%lambda(n_unique) = work%lambda(i_keep)
                  work%rho(n_unique) = work%rho(i_keep)
                  work%phi(n_unique) = work%phi(i_keep)
                  work%normals(:, n_unique) = work%normals(:, i_keep)
                  work%converged(n_unique) = work%converged(i_keep)
               end if
            end if
         end do
         n_points = n_unique
         if (allocated(keep_mask)) deallocate (keep_mask)
      end if

      if (n_points <= 0) then
         call fatal_error(error, "Projection failed for all candidate branches")
         return
      end if

      !> Populate outputs

      work%n_points = n_points
      work%branch_weights(1:n_points) = 1.0_wp

      if (present(gridpoints)) then
         allocate (gridpoints(3, n_points), source=0.0_wp)
         gridpoints(:, :) = work%points(:, 1:n_points)
      end if
      if (present(lambda_out)) then
         allocate (lambda_out(n_points), source=0.0_wp)
         lambda_out(:) = work%lambda(1:n_points)
      end if
      if (present(rho_out)) then
         allocate (rho_out(n_points), source=0.0_wp)
         rho_out(:) = work%rho(1:n_points)
      end if
      if (present(phi_out)) then
         allocate (phi_out(n_points), source=0.0_wp)
         phi_out(:) = work%phi(1:n_points)
      end if
      if (present(normal_out)) then
         allocate (normal_out(3, n_points), source=0.0_wp)
         normal_out(:, :) = work%normals(:, 1:n_points)
      end if

   end subroutine projector_project_point

   !* ================================================================================= *!
   !*                       Candidate filtering and deduplication                       *!
   !* ================================================================================= *!

   !> Filter points and return keep flags for unique representatives
   !>
   !> @param[in]  candidates   Candidate points (3,n)
   !> @param[in]  n_candidates Number of candidates
   !> @param[out] keep_mask    Keep flags (.true.=keep, .false.=duplicate)
   !> @param[out] error        Error if no valid candidate exists
   !> @param[in]  anchor       Optional anchor for debug diagnostics
   !> @param[in]  pp_sep_cut_inp Optional minimum separation for point removal
   !> @param[in]  pm_sep_cut_inp Optional rho separation for point removal
   subroutine filter_candidates(self, candidates, n_candidates, keep_mask, error, anchor, &
                                pp_sep_cut_inp, pm_sep_cut_inp)
      class(drop_projector_type), intent(inout) :: self
      real(wp), intent(in) :: candidates(:, :)
      integer, intent(in) :: n_candidates
      logical, allocatable, intent(out) :: keep_mask(:)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in), optional :: pp_sep_cut_inp
      real(wp), intent(in), optional :: pm_sep_cut_inp

      real(wp) :: dist2_i, pp_sep_cut, pm_sep_cut, rho_min
      real(wp), allocatable :: unique(:, :)
      integer :: i, j, n_selected
      type(prettylistprinter) :: minima_tbl, dist_tbl
      integer, allocatable :: dist_widths(:)
      character(:), allocatable :: dist_headers(:)
      character(32) :: idx_label

      if (n_candidates <= 0) then
         call fatal_error(error, "filter_candidates: no candidates")
         return
      end if
      if (size(candidates, dim=1) /= 3) then
         call fatal_error(error, "filter_candidates: candidates must have leading dimension 3")
         return
      end if
      if (size(candidates, dim=2) /= n_candidates) then
         call fatal_error(error, "filter_candidates: inconsistent candidate count")
         return
      end if

      if (present(pp_sep_cut_inp)) then
         pp_sep_cut = pp_sep_cut_inp
      else
         pp_sep_cut = self%branch_sep_cut
      end if

      if (present(pm_sep_cut_inp)) then
         pm_sep_cut = pm_sep_cut_inp
      else
         pm_sep_cut = self%branch_rho_cut
      end if

      n_selected = 0
      allocate (keep_mask(n_candidates), source=.false.)
      allocate (unique(3, n_candidates), source=0.0_wp)

      !> Check distance between points (deduplication)
      do i = 1, n_candidates
         if (is_new_minimum(candidates(:, i), unique, n_selected, pp_sep_cut)) then
            n_selected = n_selected + 1
            unique(:, n_selected) = candidates(:, i)
            keep_mask(i) = .true.
         end if
      end do

      !> Discard points whose anchor-point distance exceeds the minimum by more than pm_sep_cut
      rho_min = huge(1.0_wp)
      do i = 1, n_candidates
         if (keep_mask(i)) then
            dist2_i = norm2(candidates(:, i) - anchor)
            if (dist2_i < rho_min) rho_min = dist2_i
         end if
      end do
      n_selected = 0
      do i = 1, n_candidates
         if (keep_mask(i)) then
            if (norm2(candidates(:, i) - anchor) - rho_min > pm_sep_cut) then
               keep_mask(i) = .false.
            else
               n_selected = n_selected + 1
               unique(:, n_selected) = candidates(:, i)
            end if
         end if
      end do

      if (n_selected <= 0) then
         call fatal_error(error, "filter_candidates: no unique candidates after filtering")
         deallocate (unique)
         return
      end if

      if (self%debug) then
         if (n_selected > 1) then
            write (output_unit, '(x,a,i0)') 'Number of unique minima found: ', n_selected
            write (output_unit, '(x,a)')
            write (output_unit, '(x,a,3f22.14)') 'Anchor point:', anchor(1), anchor(2), anchor(3)
            write (output_unit, '(x,a)')

            minima_tbl = new_prettylistprinter( &
                         widths=[10, 22, 22, 22, 22], &
                         headers=[character(len=22) :: 'Minimum', 'Rho', 'x', 'y', 'z'], &
                         unit=output_unit, offset=1, fmt_len=22, column_gap=0)
            call minima_tbl%print_header()
            call minima_tbl%separator()
            do j = 1, n_selected
               dist2_i = norm2(unique(:, j) - anchor)
               write (idx_label, '(I0)') j
               call minima_tbl%begin_row()
               call minima_tbl%add(trim(idx_label))
               call minima_tbl%add(dist2_i, 'F22.14')
               call minima_tbl%add(unique(1, j), 'F22.14')
               call minima_tbl%add(unique(2, j), 'F22.14')
               call minima_tbl%add(unique(3, j), 'F22.14')
               call minima_tbl%end_row()
            end do

            write (output_unit, '(x,a)')
            write (output_unit, '(x,a)') 'Distance matrix between unique minima:'
            write (output_unit, '(x,a)')

            allocate (dist_widths(n_selected + 1), source=10)
            allocate (character(len=10) :: dist_headers(n_selected + 1))
            dist_headers(1) = 'i\j'
            do j = 1, n_selected
               write (idx_label, '(I0)') j
               dist_headers(j + 1) = trim(idx_label)
            end do

            dist_tbl = new_prettylistprinter( &
                       widths=dist_widths, headers=dist_headers, unit=output_unit, offset=1, &
                       fmt_len=10, fmt_real='F10.2', fmt_exp='ES10.2', column_gap=0)
            call dist_tbl%print_header()
            do i = 1, n_selected
               write (idx_label, '(I0)') i
               call dist_tbl%begin_row()
               call dist_tbl%add(trim(idx_label))
               do j = 1, n_selected
                  if (i == j) then
                     call dist_tbl%skip()
                  else
                     dist2_i = sqrt(sum((unique(:, i) - unique(:, j))**2))
                     call dist_tbl%add(dist2_i, 'ES10.2')
                  end if
               end do
               call dist_tbl%end_row()
            end do
            if (allocated(dist_widths)) deallocate (dist_widths)
            if (allocated(dist_headers)) deallocate (dist_headers)
         else
            write (output_unit, '(x,a)') 'Only one unique minimum found.'
         end if
         write (output_unit, '(x,a)')
      end if

      deallocate (unique)
   end subroutine filter_candidates

   !> Check if a minimum is distinct from previous ones
   pure logical function is_new_minimum(x_min, minima, n_unique, min_sep)
      real(wp), intent(in) :: x_min(3)
      real(wp), intent(in) :: minima(:, :)
      integer, intent(in) :: n_unique
      real(wp), intent(in) :: min_sep
      integer :: j
      real(wp) :: d2

      do j = 1, n_unique
         d2 = sum((x_min - minima(:, j))**2)
         if (d2 < min_sep**2) then
            is_new_minimum = .false.
            return
         end if
      end do
      is_new_minimum = .true.
   end function is_new_minimum

   !* ================================================================================= *!
   !*                      Newton refinement of a single candidate                      *!
   !* ================================================================================= *!

   !> Refine a single SLSQP candidate with Newton and optimality checks
   !>
   !> @param[in]     anchor      Anchor point for projection
   !> @param[in]     x_slsqp     Candidate point from SLSQP stage
   !> @param[in]     xl          Lower bounds for xyz
   !> @param[in]     xu          Upper bounds for xyz
   !> @param[in]     displacement_threshold  Maximum allowed displacement
   !> @param[in]     owner       Owner atom index
   !> @param[in]     index       Grid-point index (for diagnostics)
   !> @param[out]    error       Error if refinement fails
   !> @param[out]    gridpoint   Refined projected point
   !> @param[out]    lambda_out  Optional lagrange multiplier
   !> @param[out]    rho_out     Optional displacement norm
   !> @param[out]    phi_out     Optional final objective value at converged point
   !> @param[out]    normal_out  Optional unit surface normal
   subroutine projector_refine_point(self, anchor, x_slsqp, xl, xu, displacement_threshold, &
                                     owner, index, error, gridpoint, lambda_out, rho_out, phi_out, normal_out)
      class(drop_projector_type), intent(inout), target :: self
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: x_slsqp(3)
      real(wp), intent(in) :: xl(3), xu(3)
      real(wp), intent(in) :: displacement_threshold
      integer, intent(in) :: owner, index
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(out) :: gridpoint(3)
      real(wp), intent(out), optional :: lambda_out
      real(wp), intent(out), optional :: rho_out
      real(wp), intent(out), optional :: phi_out
      real(wp), intent(out), optional :: normal_out(3)

      logical :: is_saddle, is_kkt
      real(wp) :: lambda_curr
      real(wp) :: phi_curr
      real(wp) :: z_solution(4)
      real(wp) :: rho_actual
      class(solver_base_type), allocatable :: solver
      real(wp) :: lxl(4), lxu(4)
      type(projection_context_type) :: proj_context

      proj_context%projector => self
      proj_context%lsf => self%lsf
      proj_context%anchor(:) = anchor(:)
      proj_context%owner = owner
      proj_context%threshold = displacement_threshold

      z_solution(1:3) = x_slsqp(:)
      z_solution(4) = 0.0_wp

      ! Keep xyz inside displacement box and lambda in a wide finite range
      lxl(1:3) = xl(:)
      lxl(4) = -1.0e6_wp
      lxu(1:3) = xu(:)
      lxu(4) = 1.0e6_wp

      if (self%debug) then
         call new_newton_solver( &
            solver=solver, &
            n=4, m=4, &
            error=error, &
            func_ctx=projection_residual, &
            grad_ctx=projection_jacobian, &
            context=proj_context, &
            max_iter=self%newton_max_iter, &
            tol=self%newton_tol, &
            tolx=self%newton_tol*0.1_wp, &
            use_broyden=.false., &
            alpha=0.01_wp, &
            verbose=.true., &
            bounds_mode=3, &
            xlow=lxl, &
            xupp=lxu, &
            debug_callback_ctx=newton_debug_callback, &
            user_input_check_ctx=displacement_check &
            )
      else
         call new_newton_solver( &
            solver=solver, &
            n=4, m=4, &
            error=error, &
            func_ctx=projection_residual, &
            grad_ctx=projection_jacobian, &
            context=proj_context, &
            max_iter=self%newton_max_iter, &
            tol=self%newton_tol, &
            tolx=self%newton_tol*0.1_wp, &
            alpha=0.01_wp, &
            use_broyden=.false., &
            bounds_mode=3, &
            xlow=lxl, &
            xupp=lxu &
            )
      end if
      if (allocated(error)) then
         call fatal_error(error, "Failed to initialize Newton solver for projection refinement")
         return
      end if

      ! Solve the coupled KKT nonlinear system for [x, lambda]
      call solver%solve(z_solution, error)
      if (allocated(error)) then
         if (self%verbosity >= 3) then
            write (output_unit, '(x,i12,x,a,a)') index, '[Newton] ', trim(error%message)
         end if
         return
      end if

      if (self%verify_optimality) then

         ! Classify converged point via reduced Hessian on tangent space
         call self%check_constrained_optimality( &
            r_star=z_solution(1:3), &
            anchor=anchor, &
            owner=owner, &
            index=index, &
            lambda=z_solution(4), &
            error=error, &
            is_saddle=is_saddle, &
            is_kkt=is_kkt &
            )

         if (allocated(error)) then
            if (self%verbosity >= 3) then
               write (output_unit, '(x,i12,x,a,a)') index, '[optimality] ', trim(error%message)
            end if
            deallocate (error)
         end if

         ! Drop non-KKT points so the caller can trigger multistart fallback
         if (.not. is_kkt) then
            if (allocated(solver)) then
               call solver%destroy()
               deallocate (solver)
            end if
            call fatal_error(error, "Converged point is not a KKT point; dropping candidate")
            return
         end if

         if (is_saddle) then
            if (self%verbosity >= 3) then
               write (output_unit, '(1x,i12,x,a)') &
                  index, "[optimality] Saddle point detected, attempting escape ..."
            end if

            ! Escape negative curvature region on the surface and re-polish
            call self%riemannian_newton_escape( &
               r_init=z_solution(1:3), &
               anchor=anchor, &
               owner=owner, &
               lambda_init=z_solution(4), &
               r_out=z_solution(1:3), &
               lambda_out=z_solution(4), &
               error=error &
               )

            if (.not. allocated(solver)) then
               lxl(1:3) = xl(:)
               lxl(4) = -1.0e6_wp
               lxu(1:3) = xu(:)
               lxu(4) = 1.0e6_wp
               call new_newton_solver( &
                  solver=solver, &
                  n=4, m=4, &
                  error=error, &
                  func_ctx=projection_residual, &
                  grad_ctx=projection_jacobian, &
                  context=proj_context, &
                  max_iter=self%newton_max_iter, &
                  tol=self%newton_tol, &
                  tolx=self%newton_tol*0.1_wp, &
                  alpha=0.01_wp, &
                  use_broyden=.false., &
                  bounds_mode=3, &
                  xlow=lxl, &
                  xupp=lxu &
                  )
               if (allocated(error)) return
            end if

            call solver%solve(z_solution, error)

            ! Recheck whether escape produced a proper local minimum
            call self%check_constrained_optimality( &
               r_star=z_solution(1:3), &
               anchor=anchor, &
               owner=owner, &
               index=index, &
               lambda=z_solution(4), &
               error=error, &
               is_saddle=is_saddle, &
               is_kkt=is_kkt &
               )

            if (self%verbosity >= 3) then
               if (is_saddle .and. is_kkt) then
                  write (output_unit, '(1x,i12,x,a)') &
                     index, '[optimality] Failed to escape saddle point!'
               else
                  write (output_unit, '(1x,i12,x,a)') &
                     index, '[optimality] Escape successful!'
               end if
            end if

            if (allocated(error)) then
               if (self%verbosity >= 3) then
                  write (output_unit, '(x,i12,x,a,a)') index, &
                     '[optimality] Newton after escape failed: ', trim(error%message)
               end if
               return
            end if

            ! Drop non-KKT points after escape (mirrors first-call guard)
            if (.not. is_kkt) then
               if (allocated(solver)) then
                  call solver%destroy()
                  deallocate (solver)
               end if
               call fatal_error(error, "Escaped point is not a KKT point; dropping candidate")
               return
            end if
         end if
      end if

      if (allocated(solver)) then
         call solver%destroy()
         deallocate (solver)
      end if

      gridpoint(:) = z_solution(1:3)
      lambda_curr = z_solution(4)
      rho_actual = sqrt(dot_product(gridpoint - anchor, gridpoint - anchor))

      if (self%debug) then
         write (output_unit, '(a)') '[projector] Converged! Projection results:'
         write (output_unit, '(a,i12)') '  anchor ID   = ', index
         write (output_unit, '(a,3f12.6)') '  anchor      = ', anchor
         write (output_unit, '(a,3f12.6)') '  gridpoint   = ', gridpoint
         write (output_unit, '(a,es12.4)') '  lambda      = ', lambda_curr
         write (output_unit, '(a,f12.6,a)') '  rho         = ', rho_actual, ' Bohr'
      end if

      if (present(normal_out)) then
         normal_out = self%cached_lsf1_r
         normal_out = normal_out/norm2(normal_out)
      end if

      if (present(phi_out)) then
         call self%compute_ssd(gridpoint)
         phi_curr = self%phi%f0(gridpoint, anchor, owner)
         phi_out = phi_curr
      end if

      if (present(lambda_out)) lambda_out = lambda_curr
      if (present(rho_out)) rho_out = rho_actual

   end subroutine projector_refine_point

   !* ================================================================================= *!
   !*                            Check constrained optimality                           *!
   !* ================================================================================= *!

   !> Check whether a converged KKT solution is a constrained local minimum or saddle
   !> Uses reduced Hessian eigenvalues in the tangent space (not bordered KKT)
   !>
   !> @param[in]     r_star     Candidate solution point
   !> @param[in]     anchor     Anchor point for projection (r0)
   !> @param[in]     owner      Owner atom index
   !> @param[in]     index      Grid point index (for diagnostics)
   !> @param[in]     lambda     Lagrange multiplier at solution
   !> @param[out]    error      Error if strict validation fails
   !> @param[out]    is_saddle  Optional flag: .true. if saddle detected
   subroutine check_constrained_optimality(self, r_star, anchor, owner, index, lambda, error, &
                                           is_saddle, is_kkt, is_converged)
      class(drop_projector_type), intent(inout), target :: self
      real(wp), intent(in) :: r_star(3)
      real(wp), intent(in) :: anchor(3)
      integer, intent(in) :: owner, index
      real(wp), intent(in) :: lambda
      type(error_type), allocatable, intent(out) :: error
      logical, intent(out), optional :: is_saddle
      logical, intent(out), optional :: is_kkt
      logical, intent(out), optional :: is_converged

      ! Local variables
      real(wp) :: S_val, grad_S(3), hess_S(3, 3)
      real(wp) :: grad_phi(3), hess_phi(3, 3)
      real(wp) :: phi_val
      real(wp) :: normal_vec(3), norm_grad_S
      real(wp) :: res_feas, res_stat, grad_L(3)
      real(wp) :: H_lagrangian(3, 3)
      real(wp) :: R_reduced(2, 2)
      real(wp) :: T_mat(3, 2)  ! Tangent basis [t1 | t2]
      real(wp) :: tmp(3, 2)    ! Temporary for matrix products
      real(wp) :: mu_min, mu_max
      real(wp) :: v_min(2), v_max(2)  ! Eigenvectors (unused, but needed for interface)
      real(wp) :: eig_tol, HL_norm
      real(wp), allocatable :: lsf3_S(:, :, :)  ! Third derivative of LSF (3,3,3)
      real(wp) :: d_degen(3)    ! Degenerate direction lifted to 3D
      real(wp) :: D3_ddd        ! Third-order directional derivative D^3 L[d,d,d]
      integer :: i, j, k, i_retract

      ! Tolerances for KKT verification
      real(wp), parameter :: tol_D3 = 1.0e-8_wp  ! Tolerance for third-order test
      real(wp), parameter :: tol_feas = 1.0e-6_wp
      real(wp), parameter :: tol_stat = 1.0e-6_wp
      real(wp), parameter :: tol_grad = 1.0e-10_wp
      real(wp), parameter :: tol_eig_base = 1.0e-10_wp

      ! Initialize output
      if (present(is_saddle)) is_saddle = .false.
      if (present(is_kkt)) is_kkt = .true.
      if (present(is_converged)) is_converged = .false.

      if (self%debug) then
         write (output_unit, '(x,a,2x,a,2x,a)') &
            '===================================', &
            'Stability analysis', &
            '==================================='
      end if

      ! Step 0: Verify KKT conditions
      ! Compute LSF constraint and derivatives at solution
      call self%compute_ssd(r_star)
      call self%lsf%f012_r_screened( &
         lsf0=S_val, lsf1_r=grad_S, lsf2_rr=hess_S)

      ! Check feasibility: |S(r*)| < tol_feas
      res_feas = abs(S_val)

      ! Check constraint qualification: ||grad_S|| > tol_grad
      norm_grad_S = sqrt(dot_product(grad_S, grad_S))
      if (norm_grad_S < tol_grad) then
         if (self%verbosity >= 2) then
            write (output_unit, '(x,i12,x,a,es12.4)') index, &
               '[KKT check] Degenerate LSF gradient ||grad_S|| = ', norm_grad_S
            write (output_unit, '(x,i12,x,a)') index, &
               '[KKT check] Constraint qualification violated - classification unreliable'
         end if
         return ! Cannot reliably classify
      end if

      ! Compute stationarity residual: ||grad_phi - lambda*grad_S||
      call self%phi%f012_r(r_star, anchor, owner, phi_val, grad_phi, hess_phi)
      grad_L = grad_phi - lambda*grad_S
      res_stat = sqrt(dot_product(grad_L, grad_L))

      if (self%debug) then
         write (output_unit, '(x,a)') 'KKT verification:'
         write (output_unit, '(3x,a,es12.4,a,es12.4,a)') &
            'Feasibility   |S(r*)| = ', res_feas, '  (tol: ', tol_feas, ')'
         write (output_unit, '(3x,a,es12.4,a,es12.4,a)') &
            'Stationarity  ||dL||  = ', res_stat, '  (tol: ', tol_stat, ')'
         write (output_unit, '(3x,a,es12.4)') &
            'Gradient norm ||dS||  = ', norm_grad_S
      end if

      if (res_feas > tol_feas .or. res_stat > tol_stat) then
         if (self%verbosity >= 3 .or. self%debug) then
            write (output_unit, '(x,i12,x,a)') index, &
               '[KKT check] Not at KKT point - skipping classification'
            write (output_unit, '(x,i12,x,a,es12.4,a,es12.4,a)') index, &
               '[KKT check]   Feasibility |S(r*)| = ', res_feas, ' (tol: ', tol_feas, ')'
            write (output_unit, '(x,i12,x,a,es12.4,a,es12.4,a)') index, &
               '[KKT check]   Stationarity ||dL|| = ', res_stat, ' (tol: ', tol_stat, ')'
         end if
         if (present(is_kkt)) is_kkt = .false.
         return ! Not converged to KKT point
      end if

      ! Step 1: Build orthonormal tangent basis
      normal_vec = grad_S/norm_grad_S
      call setup_tangent_frame(normal_vec, self%work_tangent_t1, self%work_tangent_t2)

      ! Construct tangent basis matrix T = [t1 | t2] (3x2)
      T_mat(:, 1) = self%work_tangent_t1(:)
      T_mat(:, 2) = self%work_tangent_t2(:)

      ! Step 2: Form Lagrangian Hessian H_L = I - lambda * hess_S
      ! (For this specific problem: Phi = 0.5*||r - r0||^2, so hess_phi = I)
      H_lagrangian = hess_phi - lambda*hess_S

      ! Step 3: Compute reduced Hessian R = T^T * H_L * T (2x2) using intrinsic matmul
      tmp = matmul(H_lagrangian, T_mat)     ! tmp = H_L * T (3x2)
      R_reduced = matmul(transpose(T_mat), tmp)  ! R = T^T * tmp (2x2)

      ! Step 4: Compute eigenvalues of R (2x2) analytically
      call eig_2x2_symmetric(R_reduced(1, 1), R_reduced(1, 2), R_reduced(2, 2), &
                             mu_min, mu_max, v_min, v_max)

      ! Adaptive tolerance based on Lagrangian Hessian norm
      HL_norm = 0.0_wp
      do i = 1, 3
         do j = 1, 3
            HL_norm = HL_norm + H_lagrangian(i, j)**2
         end do
      end do
      HL_norm = sqrt(HL_norm)  ! Frobenius norm
      eig_tol = max(tol_eig_base, 1.0e-8_wp*HL_norm)

      if (self%debug) then
         write (output_unit, '(x,a)') 'Constrained optimality classification:'
         write (output_unit, '(3x,a,3es14.6)') 'Reduced Hessian eigenvalues: ', mu_min, mu_max
      end if

      ! Classify based on smallest eigenvalue
      if (mu_min < -eig_tol) then
         ! SADDLE or MAXIMUM on the surface
         if (self%debug) then
            write (output_unit, '(3x,a)') '=> Negative curvature in tangent direction'
            write (output_unit, '(5x,a)') '*** YIKES! SADDLE/MAXIMUM detected ***'
         end if

         ! Set saddle flag for caller to handle
         if (present(is_saddle)) is_saddle = .true.

         ! Only error if strict mode and caller doesn't handle saddle
         if (self%strict_minimum_required .and. .not. present(is_saddle)) then
            call fatal_error(error, 'Saddle point detected but strict minimum required')
            return
         end if

      else if (mu_min > eig_tol) then
         ! LOCAL MINIMUM on the surface
         if (self%debug) then
            write (output_unit, '(3x,a)') '=> Positive curvature in all tangent directions'
            write (output_unit, '(5x,a)') '*** HURRAY! CONSTRAINED LOCAL MINIMUM ***'
         end if
         if (present(is_converged)) is_converged = .true.

      else
         ! DEGENERATE / FLAT: second-order test inconclusive, try third-order
         ! The degenerate direction in tangent coords is v_min; lift to 3D:
         !   d = T * v_min
         ! Third-order directional derivative of Lagrangian:
         !   D^3 L[d,d,d] = -lambda * sum_ijk (d3S/dr_i dr_j dr_k) d_i d_j d_k
         ! (phi is quadratic, so d3(phi) = 0)
         d_degen = matmul(T_mat, v_min)

         ! Temporarily upgrade SSD to third-order derivatives for degenerate test.
         ! The SSD is normally max_deriv=2; we allocate f3_rrr_arr on the fly,
         ! recompute for the current point, then restore.
         call self%lsf%set_max_deriv(3)
         self%ssd_cache_valid = .false.
         call self%compute_ssd(r_star)

         call self%lsf%f3_rrr_screened(lsf3_rrr=lsf3_S)

         ! Restore to second-order
         call self%lsf%set_max_deriv(2)
         self%ssd_cache_valid = .false.

         D3_ddd = 0.0_wp
         do i = 1, 3
            do j = 1, 3
               do k = 1, 3
                  D3_ddd = D3_ddd + lsf3_S(i, j, k)*d_degen(i)*d_degen(j)*d_degen(k)
               end do
            end do
         end do
         D3_ddd = -lambda*D3_ddd
         deallocate (lsf3_S)

         if (self%debug) then
            write (output_unit, '(3x,a)') '=> Second-order test degenerate; checking third-order condition'
            write (output_unit, '(3x,a,3es16.8)') 'Degenerate point r* = ', r_star
            write (output_unit, '(3x,a,3es16.8)') 'Anchor         r0 = ', anchor
            write (output_unit, '(3x,a,es16.8)') 'Lambda            = ', lambda
            write (output_unit, '(3x,a,2es16.8)') 'Reduced eigs      = ', mu_min, mu_max
            write (output_unit, '(3x,a,es16.8)') 'Eig tolerance     = ', eig_tol
            write (output_unit, '(3x,a,3es16.8)') 'Degen direction d = ', d_degen
            write (output_unit, '(3x,a,es16.8)') 'D3 L[d,d,d]       = ', D3_ddd
         end if

         if (abs(D3_ddd) > tol_D3) then
            ! Nonzero third derivative => inflection point on the surface
            ! This is NOT a local minimum; treat as saddle (escapable)
            if (self%debug) then
               write (output_unit, '(5x,a)') '*** INFLECTION: nonzero D3 => treating as saddle ***'
            end if
            if (present(is_saddle)) is_saddle = .true.

            if (self%strict_minimum_required .and. .not. present(is_saddle)) then
               call fatal_error(error, 'Inflection point detected but strict minimum required')
               return
            end if
         else
            if (self%debug) then
               write (output_unit, '(5x,a)') '*** YIKES! DEGENERATE to third order (flat D2 and D3) ***'
            end if
            call fatal_error(error, 'Degenerate to third order; classification unreliable')
            return
         end if
      end if

      if (self%debug) then
         write (output_unit, '(a)') ''
      end if

   end subroutine check_constrained_optimality

   !* ================================================================================= *!
   !*                  Riemannian (?) Newton Escape from Saddle Points                  *!
   !*      (Surface walker along the constraint in the direction of lowest penalty)     *!
   !* ================================================================================= *!

   !> Riemannian Newton method to escape saddle points on constrained surface
   !> Uses negative curvature direction to escape, then descends to minimum
   !>
   !> Strategy:
   !>   1. At saddle: identify negative curvature direction in tangent space
   !>   2. Take step along negative eigenvector to escape saddle
   !>   3. Once escaped, use standard Riemannian gradient descent to reach minimum
   !>
   !> @param[in]     r_init      Initial point
   !> @param[in]     anchor      Anchor point for projection (r0)
   !> @param[in]     owner       Owner atom index
   !> @param[in]     lambda_init Initial Lagrange multiplier
   !> @param[out]    r_out       Final point (escaped from saddle)
   !> @param[out]    lambda_out  Final Lagrange multiplier
   !> @param[out]    error       Error if escape fails
   subroutine riemannian_newton_escape(self, r_init, anchor, owner, lambda_init, r_out, lambda_out, error)
      class(drop_projector_type), intent(inout), target :: self
      !> Initial point (can be at saddle)
      real(wp), intent(in) :: r_init(3)
      !> Anchor point for projection
      real(wp), intent(in) :: anchor(3)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Initial Lagrange multiplier at r_init
      real(wp), intent(in) :: lambda_init
      !> Final point after escape (hopefully a local minimum on the surface)
      real(wp), intent(out) :: r_out(3)
      !> Final Lagrange multiplier at r_out
      real(wp), intent(out) :: lambda_out
      type(error_type), allocatable, intent(out) :: error

      ! Algorithm parameters
      real(wp), parameter :: tol_grad_tan = 1.0e-8_wp

      ! Local variables
      real(wp) :: r_curr(3), r_trial(3), r_retracted(3)
      real(wp) :: lambda_curr, lambda_trial
      real(wp) :: S_val, grad_S(3), hess_S(3, 3)
      real(wp) :: phi_val, grad_phi(3), hess_phi(3, 3)
      real(wp) :: phi_curr, phi_trial
      real(wp) :: normal_vec(3), norm_grad_S
      real(wp) :: T_mat(3, 2) ! Tangent basis [t1 | t2]
      real(wp) :: g_tan(2)    ! Tangent gradient T^T * grad_phi
      real(wp) :: H_lagrangian(3, 3), R_reduced(2, 2)
      real(wp) :: tmp(3, 2)   ! Temporary for matrix products
      real(wp) :: delta_r(3)  ! Full-space step = T * s_tan
      real(wp) :: alpha, directional_deriv
      real(wp) :: norm_g_tan
      real(wp) :: rho_curr, rho_trial ! Distance from anchor
      real(wp) :: mu_min, mu_max      ! Eigenvalues of reduced Hessian
      real(wp) :: v_min(2), v_max(2)  ! Eigenvectors (in tangent coords)
      real(wp) :: s_tan(2)            ! Tangent step
      integer :: iter, ls_iter
      logical :: converged, ls_dropcess, at_saddle

      ! Initialize
      r_curr = r_init
      lambda_curr = lambda_init
      converged = .false.

      ! Get initial phi value and distance
      call self%compute_ssd(r_curr)
      phi_curr = self%phi%f0(r_curr, anchor, owner)
      rho_curr = sqrt(dot_product(r_curr - anchor, r_curr - anchor))

      if (self%debug) then
         write (output_unit, '(a)') ''
         write (output_unit, '(a)') ' ====================== Rienmannian Newton Escape ======================'
         write (output_unit, '(a)') ''
         ! TODO: re-enable per-iteration r_init/anchor/rho/Phi debug dump when needed
         write (output_unit, '(a4,1x,a12,1x,a12,1x,a12,1x,a10,1x,a12,1x,a12)') &
            'Iter', 'Phi', 'rho', '||g_tan||', 'alpha', 'mu_min', 'Status'
         write (output_unit, '(a)') repeat('-', 85)
      end if

      ! Main Riemannian optimization loop
      do iter = 0, self%riemann_max_iter

         ! Step 1: Compute derivatives at current point
         call self%compute_ssd(r_curr)
         call self%lsf%f012_r_screened( &
            lsf0=S_val, lsf1_r=grad_S, lsf2_rr=hess_S)
         call self%phi%f012_r(r_curr, anchor, owner, phi_val, grad_phi, hess_phi)

         ! Build normal and tangent basis
         norm_grad_S = sqrt(dot_product(grad_S, grad_S))
         if (norm_grad_S < 1.0e-10_wp) then
            call fatal_error(error, 'Riemannian escape: degenerate LSF gradient')
            return
         end if
         normal_vec = grad_S/norm_grad_S
         call setup_tangent_frame(normal_vec, self%work_tangent_t1, self%work_tangent_t2)
         T_mat(:, 1) = self%work_tangent_t1
         T_mat(:, 2) = self%work_tangent_t2

         ! Step 2: Compute tangent gradient g_tan = T^T * grad_phi
         g_tan(1) = dot_product(T_mat(:, 1), grad_phi)
         g_tan(2) = dot_product(T_mat(:, 2), grad_phi)
         norm_g_tan = sqrt(g_tan(1)**2 + g_tan(2)**2)

         ! Step 3: Form reduced Hessian R = T^T * H_L * T
         H_lagrangian = hess_phi - lambda_curr*hess_S

         ! R = T^T * H_L * T (2x2) using intrinsic matmul
         tmp = matmul(H_lagrangian, T_mat)     ! tmp = H_L * T (3x2)
         R_reduced = matmul(transpose(T_mat), tmp)  ! R = T^T * tmp (2x2)

         ! Step 4: Compute eigenvalues/eigenvectors of R (2x2) analytically
         ! For symmetric 2x2: R = [a b; b c]
         call eig_2x2_symmetric(R_reduced(1, 1), R_reduced(1, 2), R_reduced(2, 2), &
                                mu_min, mu_max, v_min, v_max)

         ! Note: mu_min <= mu_max by construction from eig_2x2_symmetric

         ! Determine if we're at a saddle (negative curvature)
         at_saddle = (mu_min < -self%riemann_saddle_eig_tol)

         ! Check convergence: small gradient AND positive curvature (local minimum)
         if (norm_g_tan < tol_grad_tan .and. .not. at_saddle) then
            converged = .true.
            if (self%debug) then
               write (output_unit, '(a)') ''
               write (output_unit, '(a)') '[riemannian] Converged to local minimum'
            end if
            exit
         end if

         ! Step 5: Compute search direction
         if (at_saddle) then
            ! At saddle: move along negative curvature direction
            ! Choose sign to decrease phi (descent direction)
            s_tan = v_min  ! Eigenvector for negative eigenvalue

            ! Check which direction decreases phi
            delta_r = T_mat(:, 1)*s_tan(1) + T_mat(:, 2)*s_tan(2)
            directional_deriv = dot_product(grad_phi, delta_r)

            ! If gradient component is positive, flip direction
            if (directional_deriv > 0.0_wp) then
               s_tan = -s_tan
               delta_r = -delta_r
               directional_deriv = -directional_deriv
            end if

            ! Use fixed step for saddle escape (don't use Newton)
            alpha = self%riemann_escape_step
         else
            ! Not at saddle: use Newton direction (or gradient if Hessian not useful)
            if (mu_min > self%riemann_saddle_eig_tol) then
               ! Positive definite: use Newton direction s = -R^{-1} g
               ! For 2x2: transform to eigenbasis, solve, transform back
               ! s = -V * diag(1/mu) * V^T * g
               s_tan(1) = -(v_min(1)*g_tan(1) + v_min(2)*g_tan(2))/mu_min*v_min(1) &
                          - (v_max(1)*g_tan(1) + v_max(2)*g_tan(2))/mu_max*v_max(1)
               s_tan(2) = -(v_min(1)*g_tan(1) + v_min(2)*g_tan(2))/mu_min*v_min(2) &
                          - (v_max(1)*g_tan(1) + v_max(2)*g_tan(2))/mu_max*v_max(2)
            else
               ! Near-singular: use steepest descent
               s_tan = -g_tan/max(norm_g_tan, 1.0e-10_wp)
            end if

            delta_r = T_mat(:, 1)*s_tan(1) + T_mat(:, 2)*s_tan(2)
            directional_deriv = dot_product(grad_phi, delta_r)
            alpha = self%riemann_alpha_init
         end if

         ! Step 6: Line search with retraction
         ls_dropcess = .false.

         do ls_iter = 1, self%riemann_max_ls_iter
            ! Trial point
            r_trial = r_curr + alpha*delta_r

            ! Retraction maps trial step back to S(r)=0 before acceptance test
            call self%retract_to_surface(r_trial, anchor, owner, r_retracted, lambda_trial, error)
            if (allocated(error)) then
               deallocate (error)
               alpha = alpha*self%riemann_alpha_reduce
               cycle
            end if

            ! Evaluate phi at retracted point
            call self%compute_ssd(r_retracted)
            phi_trial = self%phi%f0(r_retracted, anchor, owner)
            rho_trial = sqrt(dot_product(r_retracted - anchor, r_retracted - anchor))

            ! Accept if phi decreased (simple sufficient decrease)
            if (phi_trial < phi_curr - 1.0e-10_wp) then
               ls_dropcess = .true.
               exit
            end if

            ! Also accept via Armijo if directional derivative is negative
            if (directional_deriv < 0.0_wp) then
               if (phi_trial <= phi_curr + self%riemann_armijo_c*alpha*directional_deriv) then
                  ls_dropcess = .true.
                  exit
               end if
            end if

            alpha = alpha*self%riemann_alpha_reduce
         end do

         if (.not. ls_dropcess) then
            if (self%debug) then
               write (output_unit, '(a)') ''
               write (output_unit, '(a)') '[riemannian] Line search failed'
            end if
            exit
         end if

         ! Accept step
         r_curr = r_retracted
         lambda_curr = lambda_trial
         phi_curr = phi_trial
         rho_curr = rho_trial

         ! Print iteration info
         if (self%debug) then
            if (at_saddle) then
               write (output_unit, '(i4,1x,es12.4,1x,f12.6,1x,es12.4,1x,es10.2,1x,es12.4,1x,a12)') &
                  iter, phi_curr, rho_curr, norm_g_tan, alpha, mu_min, 'SADDLE'
            else
               write (output_unit, '(i4,1x,es12.4,1x,f12.6,1x,es12.4,1x,es10.2,1x,es12.4,1x,a12)') &
                  iter, phi_curr, rho_curr, norm_g_tan, alpha, mu_min, 'MIN-TYPE'
            end if
         end if
      end do

      ! Output final result
      r_out = r_curr
      lambda_out = lambda_curr

      if (self%debug) then
         write (output_unit, '(a)') ''
         write (output_unit, '(a)') '─────────────────────────────────────────────────────────────────'
         write (output_unit, '(a,i0)') '  Iterations:    ', iter
         write (output_unit, '(a,es14.6)') '  Final Phi:     ', phi_curr
         write (output_unit, '(a,f12.6,a)') '  Final rho:     ', rho_curr, ' Bohr'
         write (output_unit, '(a,l1)') '  Converged:     ', converged
         write (output_unit, '(a)') ''
      end if

      if (.not. converged .and. iter >= self%riemann_max_iter) then
         call fatal_error(error, 'Riemannian Newton did not converge within max iterations')
      end if

   end subroutine riemannian_newton_escape

   !* ================================================================================= *!
   !*                                 Surface retraction                                *!
   !* ================================================================================= *!

   !> Retract a point to the constraint surface S(r) = 0 using 1D Newton along normal
   !>
   !> @param[in]     r_init     Initial point
   !> @param[in]     anchor     Anchor point for projection (r0)
   !> @param[in]     owner      Owner atom index
   !> @param[out]    r_out      Retracted point on surface
   !> @param[out]    lambda_out Lagrange multiplier at retracted point
   !> @param[out]    error      Error if retraction fails
   subroutine retract_to_surface(self, r_init, anchor, owner, r_out, lambda_out, error)
      class(drop_projector_type), intent(inout), target :: self
      !> Initial point (can be off surface)
      real(wp), intent(in) :: r_init(3)
      !> Anchor point for projection
      real(wp), intent(in) :: anchor(3)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Retracted point on surface |S(r)|<tol_S
      real(wp), intent(out) :: r_out(3)
      !> Lagrange multiplier at retracted point (from KKT stationarity)
      real(wp), intent(out) :: lambda_out
      type(error_type), allocatable, intent(out) :: error

      ! Parameters TODO: should be actual projection tolerance
      real(wp), parameter :: tol_S = 1.0e-12_wp

      ! Local variables
      real(wp) :: r_curr(3), S_val, grad_S(3), grad_phi(3), norm_grad_S
      real(wp) :: step_size
      integer :: iter

      r_curr = r_init
      lambda_out = 0.0_wp

      do iter = 1, self%retract_max_iter
         ! Evaluate constraint
         call self%compute_ssd(r_curr)
         call self%lsf%f012_r_screened(lsf0=S_val, lsf1_r=grad_S)

         ! Check convergence
         if (abs(S_val) < tol_S) then
            r_out = r_curr

            ! Compute lambda from KKT stationarity: grad_phi = lambda * grad_S
            ! Use least-squares projection to handle small numerical inconsistency.
            norm_grad_S = dot_product(grad_S, grad_S)
            if (norm_grad_S < 1.0e-14_wp) then
               call fatal_error(error, 'Retraction failed: degenerate LSF gradient at converged point')
               return
            end if
            grad_phi = self%phi%f1_r(r_curr, anchor, owner)
            lambda_out = dot_product(grad_phi, grad_S)/norm_grad_S

            self%cached_lsf0 = S_val
            self%cached_lsf1_r = grad_S
            return
         end if

         ! Newton step along normal: r_new = r - S / ||grad_S||^2 * grad_S
         norm_grad_S = dot_product(grad_S, grad_S)
         if (norm_grad_S < 1.0e-14_wp) then
            call fatal_error(error, 'Retraction failed: degenerate LSF gradient')
            return
         end if

         step_size = S_val/norm_grad_S
         r_curr = r_curr - step_size*grad_S
      end do

      ! Failed to converge
      call fatal_error(error, 'Retraction to surface did not converge')

   end subroutine retract_to_surface

   !* ================================================================================= *!
   !*                              Cleanup and finalization                             *!
   !* ================================================================================= *!

   !> Clean up projector resources
   subroutine projector_destroy(self)
      class(drop_projector_type), intent(inout) :: self
      call finalize_projector(self)
   end subroutine projector_destroy

   !> Finalizer for projector type to properly deallocate all allocatable components
   !> The nested phi type has its own finalizer that will be called automatically
   subroutine finalize_projector(self)
      type(drop_projector_type), intent(inout) :: self
      call self%mol_cell_grid%destroy()
      if (allocated(self%lsf)) deallocate (self%lsf)

   end subroutine finalize_projector

end module moist_cavity_drop_projector
