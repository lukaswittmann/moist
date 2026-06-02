!> SLSQP solver wrapped with Farrell deflation for enumerating multiple
!> local minima of a constrained optimisation problem from a single seed.
!>
!> Given the caller's objective f(x) and its gradient g(x), plus a constraint
!> h(x) = 0, this solver runs SLSQP repeatedly on wrapped callbacks:
!>
!>    f_def(x) = M(x) * f(x)
!>    grad f_def(x) = M(x) * grad f(x) + f(x) * grad M(x)
!>
!> where M(x) = prod_i ( ||x - x*_i||^{-p} + alpha ) accumulates the roots
!> already found. The constraint (and its gradient) pass through unchanged, so
!> SLSQP's active-set machinery and Lagrange multipliers are unaffected and no
!> additional second-derivative information is required.
!>
!> Seeding policy is "single seed + iterated deflation": each outer iteration
!> starts SLSQP from the same caller-supplied seed; the deflation multiplier
!> guarantees the solver cannot return to an already-discovered root.
!> Termination happens when the inner SLSQP either fails or returns a point
!> within `dedup_tol` of an existing root, or when `max_roots` has been reached.
!>
!> Thread safety: the wrapped context carries pointers to the owning solver's
!> deflation operator and caller-context. Each solver instance owns its own
!> operator, so parallel callers (e.g. OMP projection) see no cross-thread
!> contention as long as each thread constructs its own solver.
module moist_math_solver_slsqp_deflation
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use iso_fortran_env, only: output_unit
   use moist_type, only: solver_base_type
   use moist_math_solver_slsqp, only: new_slsqp_solver
   use moist_math_solver_deflation, only: moist_deflation_operator_type
   implicit none
   private

   public :: moist_math_solver_slsqp_deflation_type
   public :: new_slsqp_deflation_solver

   !> Context-aware user function interfaces (mirror moist_math_solver_slsqp).
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
   end interface

   !> Internal context passed into the inner SLSQP solver. The inner solver
   !> stores a copy (via allocate/source=), which preserves the pointer
   !> components - they continue to target the owning solver's fields, so
   !> mutations to the deflation operator between outer iterations are seen
   !> by the wrapped callbacks on the next inner solve.
   type :: deflation_slsqp_context_type
      !> Caller's original objective
      procedure(objective_context_interface), pointer, nopass :: user_obj_ctx => null()
      !> Caller's original objective gradient
      procedure(objective_grad_context_interface), pointer, nopass :: user_obj_grad_ctx => null()
      !> Caller's original constraint (pass-through)
      procedure(constraints_context_interface), pointer, nopass :: user_con_ctx => null()
      !> Caller's original constraint gradient (pass-through)
      procedure(constraints_grad_context_interface), pointer, nopass :: user_con_grad_ctx => null()
      !> Pointer to the owning solver's deflation operator
      type(moist_deflation_operator_type), pointer :: deflation => null()
      !> Pointer to the caller's original user-context (owned by the solver)
      class(*), pointer :: user_context => null()
      !> Number of user (passed-through) constraints
      integer :: m_user = 0
      !> True if a ball constraint is appended (one extra inequality)
      logical :: has_ball = .false.
      !> Anchor point for the ball constraint (only used when has_ball)
      real(wp), pointer :: anchor(:) => null()
      !> Squared cap radius (phi_min + branch_rho_cut)^2; mutated between
      !> outer iterations. Set to huge() until the first root is accepted,
      !> so the inequality is trivially satisfied during the un-deflated solve.
      real(wp), pointer :: phi_max_sq => null()
   end type deflation_slsqp_context_type

   !> SLSQP-deflation solver.
   type, extends(solver_base_type) :: moist_math_solver_slsqp_deflation_type
      private
      !> Problem size
      integer :: n = 0
      integer :: m = 0
      integer :: meq = 0
      !> Variable bounds
      real(wp), allocatable :: xl(:)
      real(wp), allocatable :: xu(:)

      !> Accumulated known-root operator
      type(moist_deflation_operator_type) :: deflation

      !> Caller's original callbacks (not wrapped)
      procedure(objective_context_interface), pointer, nopass :: user_obj_ctx => null()
      procedure(objective_grad_context_interface), pointer, nopass :: user_obj_grad_ctx => null()
      procedure(constraints_context_interface), pointer, nopass :: user_con_ctx => null()
      procedure(constraints_grad_context_interface), pointer, nopass :: user_con_grad_ctx => null()
      !> Caller's original context (copied in by allocate/source=).
      class(*), allocatable :: user_context

      !> Inner SLSQP solver (polymorphic so the standard factory can be reused)
      class(solver_base_type), allocatable :: slsqp_solver

      !> Inner SLSQP tolerances / iteration limits
      real(wp) :: tol = 1.0e-8_wp
      real(wp) :: toldx = 1.0e-8_wp
      real(wp) :: toldf = 1.0e-8_wp
      integer  :: max_iter = 100

      !> Upper limit on outer deflation iterations.
      integer :: max_roots = 8

      !> Number of perturbed-seed retries if plain SLSQP fails on the first
      !> (un-deflated) attempt. Each retry offsets the seed along a fixed
      !> small set of directions scaled by retry_radius.
      integer :: max_retries = 6
      !> Offset magnitude for perturbed retries (in the units of x).
      real(wp) :: retry_radius = 0.25_wp

      !> Multiplier applied to (tol, toldx, toldf) when constructing the
      !> inner SLSQP. Deflation only needs to identify the *basin* of each
      !> root; tight convergence is left to a downstream refinement step.
      real(wp) :: tol_relax_factor = 1.0_wp

      !> Optional ball cap on subsequent roots: ||x - anchor||^2 <=
      !> (phi_min + branch_rho_cut)^2 where phi_min is the displacement
      !> norm of the first accepted root. Active only when both anchor
      !> and a positive branch_rho_cut are provided.
      logical :: has_ball = .false.
      real(wp), allocatable :: anchor(:)
      real(wp) :: branch_rho_cut = 0.0_wp
      !> Mutated between outer iterations. Reached through the wrapped
      !> context's pointer (subobject of a TARGET parent is itself a valid
      !> pointer target per F2008 16.4.1.4).
      real(wp) :: phi_max_sq = huge(1.0_wp)

      !> Converged candidate points (n_dim, n_raw_candidates)
      real(wp), allocatable :: raw_candidates(:, :)
      !> Number of successfully enumerated roots
      integer :: n_raw_candidates = 0

      !> Debug printing
      logical :: debug = .false.
   contains
      procedure :: solve => slsqp_deflation_solve
      procedure :: get_raw_candidates => slsqp_deflation_get_raw_candidates
      procedure :: destroy => slsqp_deflation_destroy
   end type moist_math_solver_slsqp_deflation_type

contains

   !> Factory: construct a new SLSQP-deflation solver.
   !>
   !> @param[out] solver         Allocated polymorphic solver handle
   !> @param[in]  n              Number of variables
   !> @param[in]  m              Total number of constraints
   !> @param[in]  meq            Number of equality constraints (meq <= m)
   !> @param[in]  obj_ctx        Caller's objective f(x; ctx)
   !> @param[in]  obj_grad_ctx   Caller's objective gradient grad f(x; ctx)
   !> @param[in]  con_ctx        Caller's constraints h(x; ctx)
   !> @param[in]  con_grad_ctx   Caller's constraint Jacobian dh/dx(x; ctx)
   !> @param[in]  context        Caller's user-context (passed back to callbacks)
   !> @param[in]  xl,xu          Variable bounds
   !> @param[in]  max_iter       Inner SLSQP iteration cap (optional)
   !> @param[in]  tol            Inner SLSQP convergence tolerance (optional)
   !> @param[in]  toldx,toldf    Inner SLSQP stagnation tolerances (optional)
   !> @param[in]  max_roots      Outer deflation iteration cap (optional, default 8)
   !> @param[in]  p_power        Deflation exponent (optional, default 2)
   !> @param[in]  alpha_shift    Deflation additive shift (optional, default 1.0)
   !> @param[in]  dedup_tol      Root-identity tolerance (optional, default 1e-6)
   !> @param[in]  max_retries    If the first (un-deflated) SLSQP from the
   !>                            caller seed fails, try this many perturbed
   !>                            seeds before giving up (optional, default 6).
   !> @param[in]  retry_radius   Magnitude of the perturbation for each retry
   !>                            (in units of x) (optional, default 0.25).
   !> @param[in]  tol_relax_factor   Multiplier on (tol, toldx, toldf) for the
   !>                            inner SLSQP. Default 1.0 (no relaxation).
   !>                            Larger values trade root accuracy for fewer
   !>                            inner iterations; deflation only needs to
   !>                            land in the correct basin, so 100x is typical.
   !> @param[in]  anchor         Optional anchor point used to enforce a ball
   !>                            cap once the first root is found. Must have
   !>                            length n if supplied.
   !> @param[in]  branch_rho_cut Maximum allowed displacement *beyond* the
   !>                            first-root distance: subsequent SLSQP solves
   !>                            see an extra inequality
   !>                            (phi_min + branch_rho_cut)^2 - ||x - anchor||^2 >= 0.
   !>                            Inactive while branch_rho_cut <= 0 or anchor
   !>                            is absent.
   !> @param[in]  debug          If true, print per-iteration diagnostics
   !> @param[out] error          Error descriptor
   subroutine new_slsqp_deflation_solver(solver, n, m, meq, &
                                         obj_ctx, obj_grad_ctx, con_ctx, con_grad_ctx, context, &
                                         xl, xu, max_iter, tol, toldx, toldf, &
                                         max_roots, p_power, alpha_shift, dedup_tol, &
                                         max_retries, retry_radius, tol_relax_factor, &
                                         anchor, branch_rho_cut, &
                                         debug, error)
      class(solver_base_type), allocatable, intent(out) :: solver
      integer, intent(in) :: n, m, meq
      procedure(objective_context_interface)         :: obj_ctx
      procedure(objective_grad_context_interface)    :: obj_grad_ctx
      procedure(constraints_context_interface)       :: con_ctx
      procedure(constraints_grad_context_interface)  :: con_grad_ctx
      class(*), intent(in) :: context
      real(wp), dimension(n), intent(in) :: xl, xu
      integer, intent(in), optional :: max_iter
      real(wp), intent(in), optional :: tol, toldx, toldf
      integer, intent(in), optional :: max_roots
      integer, intent(in), optional :: p_power
      real(wp), intent(in), optional :: alpha_shift
      real(wp), intent(in), optional :: dedup_tol
      integer, intent(in), optional :: max_retries
      real(wp), intent(in), optional :: retry_radius
      real(wp), intent(in), optional :: tol_relax_factor
      real(wp), dimension(:), intent(in), optional :: anchor
      real(wp), intent(in), optional :: branch_rho_cut
      logical, intent(in), optional :: debug
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_solver_slsqp_deflation_type), allocatable, target :: tmp
      type(error_type), allocatable :: solver_error
      type(deflation_slsqp_context_type) :: wrapped_ctx
      integer :: m_eff
      real(wp) :: tol_use, toldx_use, toldf_use

      allocate (tmp)
      tmp%n = n
      tmp%m = m
      tmp%meq = meq
      allocate (tmp%xl(n), tmp%xu(n))
      tmp%xl = xl
      tmp%xu = xu

      tmp%user_obj_ctx => obj_ctx
      tmp%user_obj_grad_ctx => obj_grad_ctx
      tmp%user_con_ctx => con_ctx
      tmp%user_con_grad_ctx => con_grad_ctx
      allocate (tmp%user_context, source=context)

      if (present(tol)) tmp%tol = tol
      if (present(toldx)) tmp%toldx = toldx
      if (present(toldf)) tmp%toldf = toldf
      if (present(max_iter)) tmp%max_iter = max_iter
      if (present(max_roots)) tmp%max_roots = max_roots
      if (present(max_retries)) tmp%max_retries = max_retries
      if (present(retry_radius)) tmp%retry_radius = retry_radius
      if (present(tol_relax_factor)) tmp%tol_relax_factor = tol_relax_factor
      if (present(debug)) tmp%debug = debug

      if (tmp%tol_relax_factor <= 0.0_wp) then
         call fatal_error(error, "SLSQP-deflation: tol_relax_factor must be positive")
         return
      end if

      ! Activate the ball cap only when the caller provides BOTH a finite
      ! positive radius and an anchor of the right size.
      tmp%has_ball = .false.
      if (present(branch_rho_cut)) tmp%branch_rho_cut = branch_rho_cut
      if (present(anchor) .and. tmp%branch_rho_cut > 0.0_wp) then
         if (size(anchor) /= n) then
            call fatal_error(error, "SLSQP-deflation: anchor size /= n")
            return
         end if
         allocate (tmp%anchor(n))
         tmp%anchor = anchor
         tmp%has_ball = .true.
      end if
      tmp%phi_max_sq = huge(1.0_wp)   ! cap inactive until first root lands

      if (tmp%max_roots <= 0) then
         call fatal_error(error, "SLSQP-deflation: max_roots must be positive")
         return
      end if

      call tmp%deflation%init(n_dim=n, max_roots=tmp%max_roots, error=solver_error, &
                              p_power=p_power, alpha_shift=alpha_shift, dedup_tol=dedup_tol)
      if (allocated(solver_error)) then
         call fatal_error(error, solver_error%message)
         return
      end if

      ! Wire the wrapped context's pointers to the solver's own state. The
      ! inner SLSQP will source-copy this context, but pointer components in
      ! the copy still target our fields.
      wrapped_ctx%user_obj_ctx => tmp%user_obj_ctx
      wrapped_ctx%user_obj_grad_ctx => tmp%user_obj_grad_ctx
      wrapped_ctx%user_con_ctx => tmp%user_con_ctx
      wrapped_ctx%user_con_grad_ctx => tmp%user_con_grad_ctx
      wrapped_ctx%deflation => tmp%deflation
      wrapped_ctx%user_context => tmp%user_context
      wrapped_ctx%m_user = tmp%m
      wrapped_ctx%has_ball = tmp%has_ball
      if (tmp%has_ball) then
         wrapped_ctx%anchor => tmp%anchor
         wrapped_ctx%phi_max_sq => tmp%phi_max_sq
      end if

      m_eff = tmp%m
      if (tmp%has_ball) m_eff = m_eff + 1
      tol_use = tmp%tol*tmp%tol_relax_factor
      toldx_use = tmp%toldx*tmp%tol_relax_factor
      toldf_use = tmp%toldf*tmp%tol_relax_factor

      call new_slsqp_solver( &
         solver=tmp%slsqp_solver, &
         n=tmp%n, m=m_eff, meq=tmp%meq, error=solver_error, &
         obj_ctx=deflated_objective, &
         obj_grad_ctx=deflated_objective_grad, &
         con_ctx=passthrough_constraint, &
         con_grad_ctx=passthrough_constraint_grad, &
         context=wrapped_ctx, &
         max_iter=tmp%max_iter, &
         tol=tol_use, toldx=toldx_use, toldf=toldf_use, &
         xl=tmp%xl, xu=tmp%xu, &
         verbose=.false.)
      if (allocated(solver_error)) then
         call fatal_error(error, solver_error%message)
         return
      end if

      call move_alloc(tmp, solver)
   end subroutine new_slsqp_deflation_solver

   !> Solve the constrained problem, enumerating up to `max_roots` distinct
   !> minima via iterated deflation. On exit `x` holds the first root
   !> discovered (the one reached from the caller's seed before any deflation
   !> is applied). All enumerated roots are stored in raw_candidates.
   !>
   !> Seeding strategy: iter=1 first tries the un-perturbed caller seed; if
   !> that fails, perturbed seeds are tried up to `max_retries` times. From
   !> iter=2 onward the un-perturbed seed is *skipped* (re-running it almost
   !> always rediscovers the previous root); each iter immediately tries up
   !> to `max_retries` perturbed seeds, accepting the first one that converges
   !> to a *new* (non-duplicate) root. This breaks the "all restarts converge
   !> to the same root" failure mode for high-symmetry geometries where
   !> deflation alone cannot redirect the inner solver.
   subroutine slsqp_deflation_solve(self, x, error)
      !> Solver instance
      class(moist_math_solver_slsqp_deflation_type), intent(inout), target :: self
      !> Initial guess on entry; first root on exit
      real(wp), dimension(:), intent(inout) :: x
      !> Error descriptor
      type(error_type), allocatable, intent(out) :: error

      type(error_type), allocatable :: inner_error
      real(wp), dimension(size(x)) :: x_seed, x_trial, x_first, perturb
      real(wp), allocatable :: converged(:, :)
      integer :: iter, attempt, k, start_attempt
      logical :: accepted, first_root_found, found_new_root

      if (size(x) /= self%n) then
         call fatal_error(error, "SLSQP-deflation: x has wrong dimension")
         return
      end if
      if (.not. allocated(self%slsqp_solver)) then
         call fatal_error(error, "SLSQP-deflation: inner SLSQP not initialized")
         return
      end if

      call self%deflation%reset()
      ! Disarm the ball cap before the first un-deflated solve. It is
      ! re-armed below as soon as the first root is accepted.
      !
      ! Use the squared box diagonal as the "inert" cap (not huge(1.0_wp)):
      ! the inner SLSQP keeps x inside [xl, xu] via its own bound
      ! constraints, so the inequality phi_max_sq - ||x - anchor||^2 >= 0
      ! is trivially satisfied for any feasible point. But SLSQP's QP
      ! subproblem normalizes constraint values with their gradients
      ! (~ 1/sqrt(c^2 + g^Tg)); a literal huge(1.0_wp) overflows c^2 and
      ! taints the line search with NaN, which on real molecules causes
      ! "no roots enumerated" failures even when the un-deflated problem
      ! is well-posed.
      if (self%has_ball) then
         self%phi_max_sq = sum((self%xu - self%xl)**2)
         if (self%phi_max_sq <= 0.0_wp) then
            self%phi_max_sq = (100.0_wp*max(self%branch_rho_cut, 1.0_wp))**2
         end if
      end if

      allocate (converged(self%n, self%max_roots), source=0.0_wp)
      self%n_raw_candidates = 0

      x_seed = x
      x_first = x
      first_root_found = .false.

      outer: do iter = 1, self%max_roots
         found_new_root = .false.
         ! iter 1: try un-perturbed seed first (attempt=0), then perturbations.
         ! iter >= 2: skip attempt=0 entirely - the un-perturbed seed already
         ! produced the first root, so re-running it just yields a duplicate.
         start_attempt = 0
         if (iter > 1) start_attempt = 1

         attempts: do attempt = start_attempt, self%max_retries
            if (attempt == 0) then
               x_trial = x_seed
            else
               ! Globally unique offset index: iter=1 uses k=1..max_retries,
               ! iter=2 uses k=max_retries+1..2*max_retries, etc.  This keeps
               ! later iterations probing fresh octants instead of recycling
               ! the same axis-aligned offsets the first iteration already
               ! tried.
               k = (iter - 1)*self%max_retries + attempt
               perturb = retry_offset(size(x), k, self%retry_radius)
               x_trial = x_seed + perturb
               call clip_to_box(x_trial, self%xl, self%xu)
            end if

            call self%slsqp_solver%solve(x_trial, inner_error)
            if (allocated(inner_error)) then
               if (self%debug) then
                  write (output_unit, '(x,a,i0,a,i0,a,a)') &
                     '[deflation] iter ', iter, ' attempt ', attempt, &
                     ' inner SLSQP failed: ', trim(inner_error%message)
               end if
               deallocate (inner_error)
               cycle attempts
            end if

            call self%deflation%append_root(x_trial, accepted)
            if (accepted) then
               found_new_root = .true.
               exit attempts
            end if
            if (self%debug) then
               write (output_unit, '(x,a,i0,a,i0,a)') &
                  '[deflation] iter ', iter, ' attempt ', attempt, &
                  ' converged to a known root (try next perturbation)'
            end if
         end do attempts

         if (.not. found_new_root) then
            if (self%debug) then
               write (output_unit, '(x,a,i0,a)') &
                  '[deflation] iter ', iter, ' exhausted all attempts (stop)'
            end if
            exit outer
         end if

         self%n_raw_candidates = self%n_raw_candidates + 1
         converged(:, self%n_raw_candidates) = x_trial
         if (.not. first_root_found) then
            x_first = x_trial
            first_root_found = .true.
            ! Activate the ball cap from iteration 2 onwards: subsequent
            ! SLSQP solves must stay inside ||x - anchor|| <= phi_min + rho_cut.
            ! The wrapped constraint sees this through context%phi_max_sq,
            ! which the inner solver's source-copied context still points
            ! to via this solver's `phi_max_sq` field.
            if (self%has_ball) then
               self%phi_max_sq = (norm2(x_trial - self%anchor) + self%branch_rho_cut)**2
               if (self%debug) then
                  write (output_unit, '(x,a,es12.4)') &
                     '[deflation] ball cap (phi_max) = ', sqrt(self%phi_max_sq)
               end if
            end if
         end if

         if (self%debug) then
            write (output_unit, '(x,a,i0,a,3(es12.4,x))') &
               '[deflation] iter ', iter, ' accepted root: ', x_trial
         end if
      end do outer

      if (allocated(self%raw_candidates)) deallocate (self%raw_candidates)
      if (self%n_raw_candidates > 0) then
         allocate (self%raw_candidates(self%n, self%n_raw_candidates))
         self%raw_candidates(:, :) = converged(:, 1:self%n_raw_candidates)
      else
         allocate (self%raw_candidates(self%n, 0))
      end if

      if (self%n_raw_candidates == 0) then
         call fatal_error(error, "SLSQP-deflation: no roots enumerated")
         deallocate (converged)
         return
      end if

      x = x_first
      deallocate (converged)
   end subroutine slsqp_deflation_solve

   !> Retry-offset pattern for `n_dim == 3`: cycle through +/- x, +/- y, +/- z
   !> scaled by `r`, so successive retries sample distinct octants. For other
   !> dimensions a deterministic pseudo-random offset is used.
   pure function retry_offset(n, k, r) result(off)
      !> Space dimension
      integer, intent(in) :: n
      !> Retry index (1-based)
      integer, intent(in) :: k
      !> Radius of the offset
      real(wp), intent(in) :: r
      real(wp), dimension(n) :: off

      integer :: axis, sign_k, i, kk
      real(wp) :: u

      off = 0.0_wp
      if (n == 3) then
         ! 6-point pattern (axis aligned) for the common 3D projection case.
         axis = mod(k - 1, 3) + 1
         sign_k = 1
         if (mod((k - 1)/3, 2) == 1) sign_k = -1
         off(axis) = real(sign_k, wp)*r
      else
         ! Deterministic spiral for other dimensions: offset_i = r*sin(k*i)/sqrt(n)
         kk = k
         do i = 1, n
            u = sin(real(kk*i, wp))
            off(i) = r*u/sqrt(real(n, wp))
         end do
      end if
   end function retry_offset

   !> Clip `x` into [xl, xu] element-wise (no error on out-of-bounds input).
   pure subroutine clip_to_box(x, xl, xu)
      !> Point to clip (modified in place)
      real(wp), dimension(:), intent(inout) :: x
      !> Lower bounds
      real(wp), dimension(:), intent(in) :: xl
      !> Upper bounds
      real(wp), dimension(:), intent(in) :: xu

      integer :: i
      do i = 1, size(x)
         if (x(i) < xl(i)) x(i) = xl(i)
         if (x(i) > xu(i)) x(i) = xu(i)
      end do
   end subroutine clip_to_box

   !> Return the full list of converged roots.
   subroutine slsqp_deflation_get_raw_candidates(self, candidates, n_candidates)
      class(moist_math_solver_slsqp_deflation_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: candidates(:, :)
      integer, intent(out) :: n_candidates

      n_candidates = self%n_raw_candidates
      if (n_candidates <= 0 .or. .not. allocated(self%raw_candidates)) then
         allocate (candidates(self%n, 0))
         return
      end if
      allocate (candidates(self%n, n_candidates))
      candidates(:, :) = self%raw_candidates(:, :)
   end subroutine slsqp_deflation_get_raw_candidates

   !> Release resources.
   subroutine slsqp_deflation_destroy(self)
      class(moist_math_solver_slsqp_deflation_type), intent(inout), target :: self

      if (allocated(self%slsqp_solver)) then
         call self%slsqp_solver%destroy()
         deallocate (self%slsqp_solver)
      end if
      call self%deflation%destroy()
      if (allocated(self%raw_candidates)) deallocate (self%raw_candidates)
      if (allocated(self%user_context)) deallocate (self%user_context)
      if (allocated(self%xl)) deallocate (self%xl)
      if (allocated(self%xu)) deallocate (self%xu)
      if (allocated(self%anchor)) deallocate (self%anchor)
      self%has_ball = .false.
      self%user_obj_ctx => null()
      self%user_obj_grad_ctx => null()
      self%user_con_ctx => null()
      self%user_con_grad_ctx => null()
      self%n = 0; self%m = 0; self%meq = 0
      self%n_raw_candidates = 0
   end subroutine slsqp_deflation_destroy

   !>==================================================================
   !> Callback wrappers installed on the inner SLSQP solver.
   !>==================================================================

   !> Deflated objective: f_def(x) = M(x) * f(x).
   subroutine deflated_objective(x, f, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      class(*), intent(in) :: context

      real(wp) :: f_raw, m_val

      f = 0.0_wp
      select type (ctx => context)
      type is (deflation_slsqp_context_type)
         if (.not. associated(ctx%user_obj_ctx)) return
         if (.not. associated(ctx%deflation)) return
         if (.not. associated(ctx%user_context)) return
         call ctx%user_obj_ctx(x, f_raw, ctx%user_context)
         m_val = ctx%deflation%multiplier(x)
         f = m_val*f_raw
      end select
   end subroutine deflated_objective

   !> Deflated objective gradient: grad f_def(x) = M*grad f + f*grad M.
   subroutine deflated_objective_grad(x, df, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: df
      class(*), intent(in) :: context

      real(wp) :: f_raw, m_val
      real(wp), dimension(size(x)) :: df_raw, grad_m

      df = 0.0_wp
      select type (ctx => context)
      type is (deflation_slsqp_context_type)
         if (.not. associated(ctx%user_obj_ctx)) return
         if (.not. associated(ctx%user_obj_grad_ctx)) return
         if (.not. associated(ctx%deflation)) return
         if (.not. associated(ctx%user_context)) return
         call ctx%user_obj_ctx(x, f_raw, ctx%user_context)
         call ctx%user_obj_grad_ctx(x, df_raw, ctx%user_context)
         m_val = ctx%deflation%multiplier(x)
         call ctx%deflation%gradient(x, grad_m)
         df = m_val*df_raw + f_raw*grad_m
      end select
   end subroutine deflated_objective_grad

   !> Constraint callback: passes the user constraints through unchanged in
   !> rows 1..m_user, and appends a single inequality
   !>     phi_max^2 - ||x - anchor||^2 >= 0
   !> in row m_user+1 when the ball cap is active. Until the first root
   !> lands, phi_max_sq == huge so the inequality is trivially satisfied.
   subroutine passthrough_constraint(x, c, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: c
      class(*), intent(in) :: context

      real(wp), allocatable :: c_user(:)
      real(wp) :: dx2

      c = 0.0_wp
      select type (ctx => context)
      type is (deflation_slsqp_context_type)
         if (.not. associated(ctx%user_con_ctx)) return
         if (.not. associated(ctx%user_context)) return
         if (ctx%m_user > 0) then
            allocate (c_user(ctx%m_user))
            call ctx%user_con_ctx(x, c_user, ctx%user_context)
            c(1:ctx%m_user) = c_user
            deallocate (c_user)
         end if
         if (ctx%has_ball .and. associated(ctx%anchor) .and. &
             associated(ctx%phi_max_sq) .and. size(c) > ctx%m_user) then
            dx2 = sum((x - ctx%anchor)**2)
            c(ctx%m_user + 1) = ctx%phi_max_sq - dx2
         end if
      end select
   end subroutine passthrough_constraint

   !> Constraint Jacobian callback: rows 1..m_user are the user Jacobian;
   !> row m_user+1 (when active) is d/dx ( phi_max^2 - ||x - anchor||^2 )
   !> = -2 (x - anchor).
   subroutine passthrough_constraint_grad(x, dc, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:, :), intent(out) :: dc
      class(*), intent(in) :: context

      real(wp), allocatable :: dc_user(:, :)

      dc = 0.0_wp
      select type (ctx => context)
      type is (deflation_slsqp_context_type)
         if (.not. associated(ctx%user_con_grad_ctx)) return
         if (.not. associated(ctx%user_context)) return
         if (ctx%m_user > 0) then
            allocate (dc_user(ctx%m_user, size(x)))
            call ctx%user_con_grad_ctx(x, dc_user, ctx%user_context)
            dc(1:ctx%m_user, :) = dc_user
            deallocate (dc_user)
         end if
         if (ctx%has_ball .and. associated(ctx%anchor) .and. &
             size(dc, 1) > ctx%m_user) then
            dc(ctx%m_user + 1, :) = -2.0_wp*(x - ctx%anchor)
         end if
      end select
   end subroutine passthrough_constraint_grad

end module moist_math_solver_slsqp_deflation
