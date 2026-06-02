!> Curvature-guided multi-start SLSQP solver
!>
!> Analyzes the constraint surface curvature at the anchor point to
!> generate targeted seed directions for multi-start SLSQP optimization.
!>
!> Instead of scattering seeds uniformly on concentric Lebedev shells
!> (O(50) SLSQP solves), this solver:
!>
!>   1. Computes the tangent-plane Hessian of the constraint at the anchor
!>   2. Identifies principal curvature directions via 2x2 eigenanalysis
!>   3. Seeds along directions with small curvature ("flat" directions)
!>      at multiple radii, where projection branch ambiguity can arise
!>   4. Falls back to the atom-anchor direction when the gradient is near-zero
!>
!> Mathematical background:
!>   For a level-set constraint S(x) = 0, the principal curvatures are
!>   eigenvalues of the shape operator W = -P nabla^2 S P / ||nabla S||,
!>   where P = I - n_hat n_hat^T is the tangent-plane projector and
!>   n_hat = nabla S / ||nabla S||.  When |kappa_i| is small, the surface
!>   is locally flat in the i-th principal direction, creating ambiguity
!>   in the closest-point projection (multiple branches may exist).
!>
!> Typical seed count: 4-12 (vs 46 for uniform Lebedev), targeting only
!> the directions where projection ambiguity actually arises.
module moist_math_solver_slsqp_curvature
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use iso_fortran_env, only: output_unit
   use moist_type, only: solver_base_type
   use moist_math_solver_slsqp, only: new_slsqp_solver
   implicit none
   private

   public :: moist_math_solver_slsqp_curvature_type
   public :: new_slsqp_curvature_solver

   !> Maximum number of perturbation radii
   integer, parameter :: max_n_radii = 10

   ! Context-aware user function interfaces (mirrors moist_math_solver_slsqp)
   abstract interface
      subroutine objective_context_interface(x, f, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x
         real(wp), intent(out) :: f
         class(*), intent(in) :: context
      end subroutine objective_context_interface

      subroutine objective_grad_context_interface(x, df, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x
         real(wp), dimension(:), intent(out) :: df
         class(*), intent(in) :: context
      end subroutine objective_grad_context_interface

      subroutine constraints_context_interface(x, c, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x
         real(wp), dimension(:), intent(out) :: c
         class(*), intent(in) :: context
      end subroutine constraints_context_interface

      subroutine constraints_grad_context_interface(x, dc, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x
         real(wp), dimension(:, :), intent(out) :: dc
         class(*), intent(in) :: context
      end subroutine constraints_grad_context_interface

      subroutine iteration_callback_context_interface(iter, x, f, c, context)
         import :: wp
         integer, intent(in) :: iter
         real(wp), dimension(:), intent(in) :: x
         real(wp), intent(in) :: f
         real(wp), dimension(:), intent(in) :: c
         class(*), intent(in) :: context
      end subroutine iteration_callback_context_interface
   end interface

   !> Curvature-guided multi-start SLSQP solver
   type, extends(solver_base_type) :: moist_math_solver_slsqp_curvature_type
      private
      !> Anchor point for closest-point projection
      real(wp) :: anchor(3)
      !> Seed starting points generated from curvature analysis
      real(wp), allocatable :: seeds(:, :)
      !> Raw converged candidates from all seed runs
      real(wp), allocatable :: raw_candidates(:, :)
      !> Number of seed points
      integer :: n_seeds = 0
      !> Number of raw converged candidates
      integer :: n_raw_candidates = 0
      !> Underlying SLSQP solver (reused for each seed)
      class(solver_base_type), allocatable :: slsqp_solver
      !> User-provided objective function
      procedure(objective_context_interface), pointer, nopass :: user_obj_ctx => null()
      !> User-provided objective gradient
      procedure(objective_grad_context_interface), pointer, nopass :: user_obj_grad_ctx => null()
      !> User-provided constraint function
      procedure(constraints_context_interface), pointer, nopass :: user_con_ctx => null()
      !> User-provided constraint gradient
      procedure(constraints_grad_context_interface), pointer, nopass :: user_con_grad_ctx => null()
      !> User-provided iteration callback
      procedure(iteration_callback_context_interface), pointer, nopass :: &
         user_iter_callback_ctx => null()
      !> User context data for callbacks
      class(*), allocatable :: user_context
      !> Number of optimization variables
      integer :: n = 0
      !> Number of constraints
      integer :: m = 0
      !> Number of equality constraints
      integer :: meq = 0
      !> Lower bounds on variables
      real(wp), allocatable :: xl(:)
      !> Upper bounds on variables
      real(wp), allocatable :: xu(:)
      !> SLSQP convergence tolerance
      real(wp) :: tol = 1.0e-2_wp
      !> SLSQP step-size stagnation tolerance
      real(wp) :: toldx = 1.0e-2_wp
      !> SLSQP objective stagnation tolerance
      real(wp) :: toldf = 1.0e-2_wp
      !> Maximum SLSQP iterations per seed
      integer :: max_iter = 50
      !> Enable debug output
      logical :: debug = .false.
      !> Curvature threshold below which a tangent direction is "flat"
      real(wp) :: kappa_thr = 0.1_wp
      !> Gradient norm threshold below which to use atom-anchor fallback normal
      real(wp) :: grad_thr = 0.01_wp
      !> Perturbation radii for seed placement
      real(wp) :: radii(max_n_radii) = 0.0_wp
      !> Number of active perturbation radii
      integer :: n_radii = 3
   contains
      !> Run multi-start SLSQP from curvature-guided seeds
      procedure :: solve => slsqp_curvature_solve
      !> Retrieve raw converged candidates
      procedure :: get_raw_candidates => slsqp_curvature_get_raw_candidates
      !> Free resources
      procedure :: destroy => slsqp_curvature_destroy
   end type moist_math_solver_slsqp_curvature_type

contains

   !> Build an orthonormal tangent-plane basis for a given unit normal.
   !>
   !> Uses Gram-Schmidt on the coordinate axis least parallel to the normal
   !> to ensure numerical stability.
   !> @param[in]  normal  Unit normal vector (3)
   !> @param[out] t1      First tangent vector (3)
   !> @param[out] t2      Second tangent vector (3), equals normal x t1
   pure subroutine build_tangent_basis(normal, t1, t2)
      real(wp), intent(in) :: normal(3)
      real(wp), intent(out) :: t1(3), t2(3)

      real(wp) :: ax(3), dot_n, inv_norm
      integer :: k

      !> Pick coordinate axis least aligned with the normal
      k = 1
      if (abs(normal(2)) < abs(normal(k))) k = 2
      if (abs(normal(3)) < abs(normal(k))) k = 3

      ax = 0.0_wp
      ax(k) = 1.0_wp

      !> Gram-Schmidt: t1 = normalize(ax - (ax . normal) * normal)
      dot_n = ax(1)*normal(1) + ax(2)*normal(2) + ax(3)*normal(3)
      t1 = ax - dot_n*normal
      inv_norm = 1.0_wp/sqrt(t1(1)**2 + t1(2)**2 + t1(3)**2)
      t1 = t1*inv_norm

      !> t2 = normal x t1
      t2(1) = normal(2)*t1(3) - normal(3)*t1(2)
      t2(2) = normal(3)*t1(1) - normal(1)*t1(3)
      t2(3) = normal(1)*t1(2) - normal(2)*t1(1)
   end subroutine build_tangent_basis

   !> Eigendecompose a 2x2 symmetric matrix analytically.
   !>
   !> For [[a11, a12], [a12, a22]], computes eigenvalues and unit
   !> eigenvectors.  Returns eigenvalues sorted by magnitude
   !> (|lam1| >= |lam2|) with corresponding eigenvectors.
   !> @param[in]  a11, a12, a22  Matrix elements
   !> @param[out] lam1           Eigenvalue with larger magnitude
   !> @param[out] lam2           Eigenvalue with smaller magnitude
   !> @param[out] v1             Unit eigenvector for lam1 (2-component)
   !> @param[out] v2             Unit eigenvector for lam2 (2-component)
   pure subroutine eig_2x2_sym(a11, a12, a22, lam1, lam2, v1, v2)
      real(wp), intent(in) :: a11, a12, a22
      real(wp), intent(out) :: lam1, lam2
      real(wp), intent(out) :: v1(2), v2(2)

      real(wp) :: avg, diff, disc, e1, e2, inv_norm

      avg = 0.5_wp*(a11 + a22)
      diff = 0.5_wp*(a11 - a22)
      disc = sqrt(diff**2 + a12**2)

      e1 = avg + disc
      e2 = avg - disc

      !> Sort by magnitude: |lam1| >= |lam2|
      if (abs(e1) >= abs(e2)) then
         lam1 = e1
         lam2 = e2
      else
         lam1 = e2
         lam2 = e1
      end if

      !> Compute eigenvectors using (A - lam*I)v = 0
      !> From first row: (a11 - lam)*v1 + a12*v2 = 0  =>  v = [a12, lam - a11]
      if (abs(a12) > 1.0e-14_wp) then
         v1 = [a12, lam1 - a11]
         inv_norm = 1.0_wp/sqrt(v1(1)**2 + v1(2)**2)
         v1 = v1*inv_norm

         v2 = [a12, lam2 - a11]
         inv_norm = 1.0_wp/sqrt(v2(1)**2 + v2(2)**2)
         v2 = v2*inv_norm
      else
         !> Already diagonal
         if (abs(a11) >= abs(a22)) then
            v1 = [1.0_wp, 0.0_wp]
            v2 = [0.0_wp, 1.0_wp]
         else
            v1 = [0.0_wp, 1.0_wp]
            v2 = [1.0_wp, 0.0_wp]
         end if
      end if
   end subroutine eig_2x2_sym

   !> Generate curvature-guided seed points around the anchor.
   !>
   !> Analyzes the constraint Hessian projected onto the tangent plane
   !> of the constraint surface at the anchor.  Seeds are placed along
   !> principal curvature directions where the curvature is small (the
   !> surface is locally flat and projection branches may exist).
   !> Directions with large curvature receive only a minimal perturbation.
   !> Falls back to a small Lebedev shell when the gradient norm is too
   !> small to define a reliable surface normal.
   !>
   !> @param[in]    anchor     Anchor point (3)
   !> @param[in]    grad_s     Constraint gradient at anchor (3)
   !> @param[in]    hess_s     Constraint Hessian at anchor (3,3)
   !> @param[in]    radii      Perturbation radii for seed placement
   !> @param[in]    kappa_thr  Curvature threshold for "flat" classification
   !> @param[in]    grad_thr   Gradient norm threshold for Lebedev fallback
   !> @param[out]   seeds      Generated seed points (3, n_seeds)
   !> @param[out]   n_seeds    Number of seeds generated
   !> @param[in]    debug      Enable debug output
   !> @param[out]   error      Error status
   subroutine generate_curvature_seeds(anchor, grad_s, hess_s, radii, &
                                       kappa_thr, grad_thr, owner_pos, seeds, n_seeds, debug, error)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: grad_s(3)
      real(wp), intent(in) :: hess_s(3, 3)
      real(wp), intent(in) :: radii(:)
      real(wp), intent(in) :: kappa_thr
      real(wp), intent(in) :: grad_thr
      real(wp), intent(in) :: owner_pos(3)
      real(wp), allocatable, intent(out) :: seeds(:, :)
      integer, intent(out) :: n_seeds
      logical, intent(in) :: debug
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: grad_norm, normal(3), t1(3), t2(3)
      real(wp) :: R11, R12, R22, lam1, lam2, v1(2), v2(2)
      real(wp) :: kappa1, kappa2, dir1(3), dir2(3)
      real(wp) :: Ht1(3), Ht2(3)
      real(wp) :: fallback_dir(3), fallback_norm
      integer :: max_seeds, n_radii, ir, offset
      logical :: flat1, flat2

      n_radii = size(radii)
      grad_norm = sqrt(grad_s(1)**2 + grad_s(2)**2 + grad_s(3)**2)

      ! ----------------------------------------------------------------
      ! Case 1: Gradient too small - use atom-anchor direction as
      !         fallback normal, seed along both tangent directions
      ! ----------------------------------------------------------------
      if (grad_norm < grad_thr) then
         fallback_dir = anchor - owner_pos
         fallback_norm = sqrt(sum(fallback_dir**2))
         normal = fallback_dir/fallback_norm
         call build_tangent_basis(normal, t1, t2)

         if (debug) then
            write (output_unit, '(x,a,es10.3,a,es10.3)') &
               '[curvature] Gradient norm ', grad_norm, ' < threshold ', grad_thr
            write (output_unit, '(x,a)') &
               '[curvature] Using atom-anchor fallback normal'
            write (output_unit, '(x,a,3f10.5)') &
               '[curvature] Fallback normal: ', normal
         end if

         !> Both tangent directions are ambiguous - seed along both at all radii
         n_seeds = 4*n_radii
         allocate (seeds(3, n_seeds))
         offset = 0
         do ir = 1, n_radii
            offset = offset + 1
            seeds(:, offset) = anchor + radii(ir)*t1
            offset = offset + 1
            seeds(:, offset) = anchor - radii(ir)*t1
            offset = offset + 1
            seeds(:, offset) = anchor + radii(ir)*t2
            offset = offset + 1
            seeds(:, offset) = anchor - radii(ir)*t2
         end do
         return
      end if

      ! ----------------------------------------------------------------
      ! Case 2: Gradient OK - curvature-guided seeding
      ! ----------------------------------------------------------------
      normal = grad_s/grad_norm
      call build_tangent_basis(normal, t1, t2)

      !> Project Hessian onto tangent plane: R_ij = t_i . H . t_j
      Ht1 = matmul(hess_s, t1)
      Ht2 = matmul(hess_s, t2)
      R11 = dot_product(t1, Ht1)
      R12 = dot_product(t1, Ht2)
      R22 = dot_product(t2, Ht2)

      !> Eigendecompose the 2x2 projected Hessian
      call eig_2x2_sym(R11, R12, R22, lam1, lam2, v1, v2)

      !> Convert to principal curvatures: kappa = |lambda| / ||grad S||
      kappa1 = abs(lam1)/grad_norm
      kappa2 = abs(lam2)/grad_norm

      !> Convert eigenvectors to 3D principal directions
      dir1 = v1(1)*t1 + v1(2)*t2
      dir2 = v2(1)*t1 + v2(2)*t2

      flat1 = (kappa1 < kappa_thr)
      flat2 = (kappa2 < kappa_thr)

      if (debug) then
         write (output_unit, '(x,a)') &
            '========== Curvature-guided seeding ========='
         write (output_unit, '(x,a,es12.4)') &
            'Gradient norm:  ', grad_norm
         write (output_unit, '(x,a,3f10.5)') &
            'Normal:         ', normal
         write (output_unit, '(x,a,es12.4,a,l3)') &
            'kappa_1:        ', kappa1, '  flat:', flat1
         write (output_unit, '(x,a,3f10.5)') &
            'Direction 1:    ', dir1
         write (output_unit, '(x,a,es12.4,a,l3)') &
            'kappa_2:        ', kappa2, '  flat:', flat2
         write (output_unit, '(x,a,3f10.5)') &
            'Direction 2:    ', dir2
      end if

      !> Compute seed count:
      !>   - Flat directions: +/- at each radius  =>  2 * n_radii seeds
      !>   - Curved directions: +/- at smallest radius only  =>  2 seeds
      max_seeds = 0
      if (flat1) then
         max_seeds = max_seeds + 2*n_radii
      else
         max_seeds = max_seeds + 2
      end if
      if (flat2) then
         max_seeds = max_seeds + 2*n_radii
      else
         max_seeds = max_seeds + 2
      end if

      allocate (seeds(3, max_seeds))
      offset = 0

      !> Direction 1: flat -> all radii, curved -> smallest radius only
      if (flat1) then
         do ir = 1, n_radii
            offset = offset + 1
            seeds(:, offset) = anchor + radii(ir)*dir1
            offset = offset + 1
            seeds(:, offset) = anchor - radii(ir)*dir1
         end do
      else
         offset = offset + 1
         seeds(:, offset) = anchor + radii(1)*dir1
         offset = offset + 1
         seeds(:, offset) = anchor - radii(1)*dir1
      end if

      !> Direction 2: flat -> all radii, curved -> smallest radius only
      if (flat2) then
         do ir = 1, n_radii
            offset = offset + 1
            seeds(:, offset) = anchor + radii(ir)*dir2
            offset = offset + 1
            seeds(:, offset) = anchor - radii(ir)*dir2
         end do
      else
         offset = offset + 1
         seeds(:, offset) = anchor + radii(1)*dir2
         offset = offset + 1
         seeds(:, offset) = anchor - radii(1)*dir2
      end if

      n_seeds = offset

      if (debug) then
         write (output_unit, '(x,a,i0,a)') &
            '[curvature] Generated ', n_seeds, ' seeds'
         write (output_unit, '(x,a)') &
            '============================================='
      end if
   end subroutine generate_curvature_seeds

   !> Factory function to create a curvature-guided multi-start SLSQP solver.
   !>
   !> Accepts the constraint gradient and Hessian at the anchor to perform
   !> curvature analysis and generate targeted seeds.  Interface mirrors
   !> new_slsqp_multistart_solver but replaces Lebedev shell parameters
   !> with grad_s, hess_s, and optional curvature thresholds.
   !>
   !> @param[in]    anchor     Anchor point for projection (3)
   !> @param[out]   solver     Allocated solver instance
   !> @param[in]    n          Number of optimization variables (must be 3)
   !> @param[in]    m          Number of constraints
   !> @param[in]    meq        Number of equality constraints
   !> @param[in]    obj_ctx    Objective function callback (context-aware)
   !> @param[in]    obj_grad_ctx  Objective gradient callback (context-aware)
   !> @param[in]    con_ctx    Constraint function callback (context-aware)
   !> @param[in]    con_grad_ctx  Constraint gradient callback (context-aware)
   !> @param[in]    context    User context for callbacks
   !> @param[in]    xl         Lower bounds on variables (n)
   !> @param[in]    xu         Upper bounds on variables (n)
   !> @param[in]    grad_s     Constraint gradient at anchor (3)
   !> @param[in]    hess_s     Constraint Hessian at anchor (3,3)
   !> @param[in]    owner_pos  Position of owner atom (3), used as fallback normal
   !> @param[in]    max_iter   Maximum SLSQP iterations per seed (optional)
   !> @param[in]    tol        SLSQP convergence tolerance (optional)
   !> @param[in]    toldx      SLSQP step tolerance (optional)
   !> @param[in]    toldf      SLSQP function tolerance (optional)
   !> @param[in]    iter_callback_ctx  Iteration callback (optional)
   !> @param[in]    radii      Perturbation radii (optional, overrides type default)
   !> @param[in]    kappa_thr  Curvature threshold (optional, overrides type default)
   !> @param[in]    grad_thr   Gradient threshold (optional, overrides type default)
   !> @param[in]    debug      Enable debug output (optional)
   !> @param[out]   error      Error status (optional)
   subroutine new_slsqp_curvature_solver(anchor, solver, &
                                         n, m, meq, obj_ctx, obj_grad_ctx, con_ctx, con_grad_ctx, context, &
                                         xl, xu, grad_s, hess_s, owner_pos, &
                                         max_iter, tol, toldx, toldf, iter_callback_ctx, &
                                         radii, kappa_thr, grad_thr, debug, error)
      real(wp), intent(in) :: anchor(3)
      class(solver_base_type), allocatable, intent(out) :: solver
      integer, intent(in) :: n
      integer, intent(in) :: m
      integer, intent(in) :: meq
      procedure(objective_context_interface) :: obj_ctx
      procedure(objective_grad_context_interface) :: obj_grad_ctx
      procedure(constraints_context_interface) :: con_ctx
      procedure(constraints_grad_context_interface) :: con_grad_ctx
      class(*), intent(in) :: context
      real(wp), dimension(n), intent(in) :: xl
      real(wp), dimension(n), intent(in) :: xu
      real(wp), intent(in) :: grad_s(3)
      real(wp), intent(in) :: hess_s(3, 3)
      real(wp), intent(in) :: owner_pos(3)
      integer, intent(in), optional :: max_iter
      real(wp), intent(in), optional :: tol
      real(wp), intent(in), optional :: toldx
      real(wp), intent(in), optional :: toldf
      procedure(iteration_callback_context_interface), optional :: iter_callback_ctx
      real(wp), intent(in), optional :: radii(:)
      real(wp), intent(in), optional :: kappa_thr
      real(wp), intent(in), optional :: grad_thr
      logical, intent(in), optional :: debug
      type(error_type), allocatable, intent(out), optional :: error

      type(moist_math_solver_slsqp_curvature_type), allocatable :: tmp
      type(error_type), allocatable :: seed_error, solver_error
      logical :: debug_use

      allocate (tmp)
      tmp%anchor = anchor
      tmp%n = n
      tmp%m = m
      tmp%meq = meq
      tmp%user_obj_ctx => obj_ctx
      tmp%user_obj_grad_ctx => obj_grad_ctx
      tmp%user_con_ctx => con_ctx
      tmp%user_con_grad_ctx => con_grad_ctx
      if (present(iter_callback_ctx)) tmp%user_iter_callback_ctx => iter_callback_ctx
      allocate (tmp%user_context, source=context)
      allocate (tmp%xl(n), tmp%xu(n))
      tmp%xl = xl
      tmp%xu = xu

      if (present(tol)) tmp%tol = tol
      if (present(toldx)) tmp%toldx = toldx
      if (present(toldf)) tmp%toldf = toldf
      if (present(max_iter)) tmp%max_iter = max_iter

      debug_use = .false.
      if (present(debug)) debug_use = debug
      tmp%debug = debug_use

      !> Override curvature defaults from optional arguments
      if (present(kappa_thr)) tmp%kappa_thr = kappa_thr
      if (present(grad_thr)) tmp%grad_thr = grad_thr
      if (present(radii)) then
         tmp%n_radii = min(size(radii), max_n_radii)
         tmp%radii = 0.0_wp
         tmp%radii(:tmp%n_radii) = radii(:tmp%n_radii)
      else
         tmp%radii(1) = 0.2_wp
         tmp%radii(2) = 0.5_wp
         tmp%radii(3) = 0.8_wp
      end if

      if (debug_use) then
         write (output_unit, '(x,a)') &
            '========== Curvature SLSQP startup =========='
      end if

      !> Generate curvature-guided seeds
      call generate_curvature_seeds(anchor, grad_s, hess_s, &
                                    tmp%radii(:tmp%n_radii), tmp%kappa_thr, tmp%grad_thr, owner_pos, &
                                    tmp%seeds, tmp%n_seeds, debug_use, seed_error)

      if (allocated(seed_error)) then
         if (present(error)) then
            call fatal_error(error, seed_error%message)
         else
            error stop seed_error%message
         end if
         return
      end if

      !> Create underlying SLSQP solver (reused for all seed runs)
      call new_slsqp_solver( &
         solver=tmp%slsqp_solver, &
         n=tmp%n, m=tmp%m, meq=tmp%meq, &
         error=solver_error, &
         obj_ctx=tmp%user_obj_ctx, &
         obj_grad_ctx=tmp%user_obj_grad_ctx, &
         con_ctx=tmp%user_con_ctx, &
         con_grad_ctx=tmp%user_con_grad_ctx, &
         context=tmp%user_context, &
         max_iter=tmp%max_iter, &
         tol=tmp%tol, &
         toldx=tmp%toldx, &
         toldf=tmp%toldf, &
         xl=tmp%xl, &
         xu=tmp%xu &
         )

      if (allocated(solver_error)) then
         if (present(error)) then
            call fatal_error(error, solver_error%message)
         else
            error stop solver_error%message
         end if
         return
      end if

      call move_alloc(tmp, solver)
   end subroutine new_slsqp_curvature_solver

   !> Solve the constrained projection using curvature-guided multi-start SLSQP.
   !>
   !> Iterates over all seed points, runs SLSQP from each, and returns the
   !> feasible solution closest to the anchor.  All converged results are
   !> stored as raw candidates for downstream filtering.
   !> @param[inout] x      Initial guess in, best solution out (3)
   !> @param[out]   error  Error status
   subroutine slsqp_curvature_solve(self, x, error)
      class(moist_math_solver_slsqp_curvature_type), intent(inout), target :: self
      real(wp), dimension(:), intent(inout) :: x
      type(error_type), allocatable, intent(out) :: error

      type(error_type), allocatable :: solver_error
      real(wp) :: x_trial(3), best_x(3), best_dist2, dist2
      integer :: i, n_converged
      real(wp), allocatable :: converged(:, :)

      if (size(x) /= 3) then
         call fatal_error(error, "Curvature SLSQP solver: x must have size 3")
         return
      end if
      if (.not. allocated(self%seeds)) then
         call fatal_error(error, "Curvature SLSQP solver: seeds not initialized")
         return
      end if
      if (self%n_seeds <= 0 .or. size(self%seeds, dim=2) /= self%n_seeds) then
         call fatal_error(error, "Curvature SLSQP solver: invalid seed count")
         return
      end if
      if (.not. allocated(self%slsqp_solver)) then
         call fatal_error(error, "Curvature SLSQP solver: SLSQP not initialized")
         return
      end if

      best_dist2 = huge(1.0_wp)
      best_x = self%anchor

      allocate (converged(3, self%n_seeds))
      n_converged = 0

      do i = 1, self%n_seeds
         x_trial = self%seeds(:, i)

         call self%slsqp_solver%solve(x_trial, solver_error)
         if (allocated(solver_error)) then
            if (self%debug) then
               write (output_unit, '(x,a,i0,a,a)') &
                  'Seed ', i, ' failed: ', trim(solver_error%message)
            end if
            deallocate (solver_error)
            cycle
         end if

         dist2 = sum((x_trial - self%anchor)**2)
         if (dist2 < best_dist2) then
            best_dist2 = dist2
            best_x = x_trial
         end if

         n_converged = n_converged + 1
         converged(:, n_converged) = x_trial
      end do

      if (best_dist2 >= huge(1.0_wp)*0.5_wp) then
         call fatal_error(error, "Curvature SLSQP solver: no successful starts")
         return
      end if

      if (allocated(self%raw_candidates)) deallocate (self%raw_candidates)
      self%n_raw_candidates = n_converged
      allocate (self%raw_candidates(3, n_converged))
      self%raw_candidates(:, :) = converged(:, 1:n_converged)

      x = best_x

      if (self%debug) then
         write (output_unit, '(x,a,i0,a,i0,a)') &
            '[curvature] ', n_converged, '/', self%n_seeds, ' seeds converged'
         write (output_unit, '(x,a,es12.4)') &
            '[curvature] Best distance: ', sqrt(best_dist2)
      end if

      deallocate (converged)
   end subroutine slsqp_curvature_solve

   !> Get raw SLSQP candidates converged from curvature-guided seeds.
   !> @param[out] candidates   Raw candidate points (3, n_candidates)
   !> @param[out] n_candidates Number of available candidates
   subroutine slsqp_curvature_get_raw_candidates(self, candidates, n_candidates)
      class(moist_math_solver_slsqp_curvature_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: candidates(:, :)
      integer, intent(out) :: n_candidates

      n_candidates = self%n_raw_candidates
      if (n_candidates <= 0 .or. .not. allocated(self%raw_candidates)) then
         allocate (candidates(3, 0))
         return
      end if
      allocate (candidates(3, n_candidates))
      candidates(:, :) = self%raw_candidates(:, :)
   end subroutine slsqp_curvature_get_raw_candidates

   !> Clean up resources
   subroutine slsqp_curvature_destroy(self)
      class(moist_math_solver_slsqp_curvature_type), intent(inout), target :: self

      if (allocated(self%slsqp_solver)) then
         call self%slsqp_solver%destroy()
         deallocate (self%slsqp_solver)
      end if
      if (allocated(self%seeds)) deallocate (self%seeds)
      if (allocated(self%raw_candidates)) deallocate (self%raw_candidates)
      self%n_raw_candidates = 0
      if (allocated(self%user_context)) deallocate (self%user_context)
      if (allocated(self%xl)) deallocate (self%xl)
      if (allocated(self%xu)) deallocate (self%xu)
      self%user_obj_ctx => null()
      self%user_obj_grad_ctx => null()
      self%user_con_ctx => null()
      self%user_con_grad_ctx => null()
      self%user_iter_callback_ctx => null()
   end subroutine slsqp_curvature_destroy

end module moist_math_solver_slsqp_curvature
