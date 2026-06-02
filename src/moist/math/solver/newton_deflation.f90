!> Newton-Raphson solver wrapped with Farrell deflation for enumerating
!> multiple roots of a nonlinear system F(x) = 0 from a single seed.
!>
!> Given the caller's residual F(x) and its Jacobian J(x) = dF/dx, this solver
!> runs Newton's method repeatedly on wrapped callbacks:
!>
!>    F_def(x) = M(x) * F(x)
!>    J_def(x) = M(x) * J(x) + F(x) * grad M(x)^T
!>
!> where M(x) = prod_i ( ||x - x*_i||^{-p} + alpha ) accumulates the roots
!> already found, and the Jacobian correction is a rank-1 outer-product
!> update - each row of the correction is F_i * grad M.
!>
!> Seeding policy is "single seed + iterated deflation". Termination when the
!> inner Newton solver fails or returns a point within `dedup_tol` of a known
!> root, or when `max_roots` has been reached.
!>
!> Thread safety: context carries pointers to owning solver state; each solver
!> instance owns its own deflation operator, so multiple instances running in
!> parallel do not collide.
module moist_math_solver_newton_deflation
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use iso_fortran_env, only: output_unit
   use moist_type, only: solver_base_type
   use moist_math_solver_newton, only: new_newton_solver
   use moist_math_solver_deflation, only: moist_deflation_operator_type
   implicit none
   private

   public :: moist_math_solver_newton_deflation_type
   public :: new_newton_deflation_solver

   !> Context-aware user function interfaces (mirror moist_math_solver_newton).
   abstract interface
      subroutine func_context_interface(x, f, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x
         real(wp), dimension(:), intent(out) :: f
         class(*), intent(in) :: context
      end subroutine func_context_interface

      subroutine grad_context_interface(x, jac, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x
         real(wp), dimension(:, :), intent(out) :: jac
         class(*), intent(in) :: context
      end subroutine grad_context_interface
   end interface

   !> Internal context passed into the inner Newton solver. The inner solver
   !> stores a copy via allocate/source=; pointer components in the copy still
   !> target the owning solver's fields, so mutations to the deflation
   !> operator between outer iterations are visible to the wrapped callbacks.
   type :: deflation_newton_context_type
      procedure(func_context_interface), pointer, nopass :: user_func_ctx => null()
      procedure(grad_context_interface), pointer, nopass :: user_grad_ctx => null()
      type(moist_deflation_operator_type), pointer :: deflation => null()
      class(*), pointer :: user_context => null()
   end type deflation_newton_context_type

   !> Newton-deflation solver.
   type, extends(solver_base_type) :: moist_math_solver_newton_deflation_type
      private
      !> Problem size (n variables, m residuals). For square systems n == m.
      integer :: n = 0
      integer :: m = 0

      !> Accumulated known-root operator
      type(moist_deflation_operator_type) :: deflation

      !> Caller's original callbacks
      procedure(func_context_interface), pointer, nopass :: user_func_ctx => null()
      procedure(grad_context_interface), pointer, nopass :: user_grad_ctx => null()
      !> Caller's original context
      class(*), allocatable :: user_context

      !> Inner Newton solver
      class(solver_base_type), allocatable :: newton_solver

      !> Inner-solver tolerances
      real(wp) :: tol = 1.0e-8_wp
      real(wp) :: tolx = 1.0e-10_wp
      integer  :: max_iter = 100

      !> Multiplier on (tol, tolx) for the inner Newton; deflation only
      !> needs to land in the correct basin. Default 1.0 (no relaxation).
      real(wp) :: tol_relax_factor = 1.0_wp

      !> Outer deflation iteration cap
      integer :: max_roots = 8

      !> Inner Newton damping. Mirrors `alpha` of new_newton_solver: smaller
      !> values trade speed for staying inside the basin on stiff problems
      !> (the projection KKT system needs alpha ~ 0.01).
      real(wp) :: alpha = 1.0_wp
      logical  :: use_broyden = .false.

      !> Number of perturbed-seed retries on the first (un-deflated) Newton
      !> call. Mirrors the SLSQP-deflation retry pattern; needed because a
      !> single-seed solver loses to Lebedev multistart on pathological anchors.
      integer :: max_retries = 6
      real(wp) :: retry_radius = 0.25_wp

      !> Converged roots (n, n_raw_candidates)
      real(wp), allocatable :: raw_candidates(:, :)
      !> Number of successfully enumerated roots
      integer :: n_raw_candidates = 0

      !> Bounds plumbing for Newton (ignored when bounds_mode == 0).
      integer :: bounds_mode = 0
      real(wp), allocatable :: xl_init(:)
      real(wp), allocatable :: xu_init(:)

      !> Optional ball cap on subsequent roots: once the first root is
      !> accepted, the inner Newton is rebuilt with axis-aligned bounds
      !> that intersect [xl_init, xu_init] with [anchor - phi_max, anchor + phi_max]
      !> on the first n_anchor components, where
      !>     phi_max = phi_min + branch_rho_cut.
      !> Useful when the residual is the projection KKT system on
      !> z = (x, lambda): n_anchor == 3 and lambda keeps its wide bounds.
      logical :: has_ball = .false.
      real(wp), allocatable :: anchor(:)
      integer  :: n_anchor = 0
      real(wp) :: branch_rho_cut = 0.0_wp

      logical :: debug = .false.
   contains
      procedure :: solve => newton_deflation_solve
      procedure :: get_raw_candidates => newton_deflation_get_raw_candidates
      procedure :: destroy => newton_deflation_destroy
      procedure, private :: build_inner_newton => newton_deflation_build_inner
   end type moist_math_solver_newton_deflation_type

contains

   !> Factory: construct a new Newton-deflation solver.
   !>
   !> @param[out] solver         Allocated polymorphic solver handle
   !> @param[in]  n              Number of variables
   !> @param[in]  m              Number of residual components (often == n)
   !> @param[in]  func_ctx       Caller's residual F(x; ctx)
   !> @param[in]  grad_ctx       Caller's Jacobian dF/dx(x; ctx)
   !> @param[in]  context        Caller's user-context
   !> @param[in]  max_iter       Inner Newton iteration cap (optional)
   !> @param[in]  tol            Inner Newton residual tolerance (optional)
   !> @param[in]  tolx           Inner Newton step tolerance (optional)
   !> @param[in]  max_roots         Outer deflation iteration cap (optional, default 8)
   !> @param[in]  p_power           Deflation exponent (optional, default 2)
   !> @param[in]  alpha_shift       Deflation additive shift (optional, default 1.0)
   !> @param[in]  dedup_tol         Root-identity tolerance (optional, default 1e-6)
   !> @param[in]  tol_relax_factor  Multiplier on (tol, tolx). Default 1.0.
   !> @param[in]  bounds_mode       Newton bounds mode (forwarded to inner Newton)
   !> @param[in]  xlow              Initial lower bounds for inner Newton (length n)
   !> @param[in]  xupp              Initial upper bounds for inner Newton (length n)
   !> @param[in]  anchor            Anchor used for the post-first-root box cap.
   !>                               Must have length n_anchor (= size(anchor)).
   !>                               Only the first n_anchor entries of x are
   !>                               capped; remaining entries (e.g. lambda in
   !>                               a 4-D KKT system) keep [xlow, xupp].
   !> @param[in]  branch_rho_cut    Maximum allowed displacement beyond the
   !>                               first root: phi_max = phi_min + branch_rho_cut.
   !>                               Inactive when <= 0 or anchor is absent.
   !> @param[in]  debug             If true, print per-iteration diagnostics
   !> @param[out] error             Error descriptor
   subroutine new_newton_deflation_solver(solver, n, m, &
                                          func_ctx, grad_ctx, context, &
                                          max_iter, tol, tolx, tol_relax_factor, &
                                          alpha, use_broyden, max_retries, retry_radius, &
                                          max_roots, p_power, alpha_shift, dedup_tol, &
                                          bounds_mode, xlow, xupp, &
                                          anchor, branch_rho_cut, &
                                          debug, error)
      class(solver_base_type), allocatable, intent(out) :: solver
      integer, intent(in) :: n, m
      procedure(func_context_interface) :: func_ctx
      procedure(grad_context_interface) :: grad_ctx
      class(*), intent(in) :: context
      integer, intent(in), optional :: max_iter
      real(wp), intent(in), optional :: tol, tolx
      real(wp), intent(in), optional :: tol_relax_factor
      real(wp), intent(in), optional :: alpha
      logical, intent(in), optional :: use_broyden
      integer, intent(in), optional :: max_retries
      real(wp), intent(in), optional :: retry_radius
      integer, intent(in), optional :: max_roots
      integer, intent(in), optional :: p_power
      real(wp), intent(in), optional :: alpha_shift
      real(wp), intent(in), optional :: dedup_tol
      integer, intent(in), optional :: bounds_mode
      real(wp), dimension(:), intent(in), optional :: xlow, xupp
      real(wp), dimension(:), intent(in), optional :: anchor
      real(wp), intent(in), optional :: branch_rho_cut
      logical, intent(in), optional :: debug
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_solver_newton_deflation_type), allocatable, target :: tmp
      type(error_type), allocatable :: solver_error

      allocate (tmp)
      tmp%n = n
      tmp%m = m
      tmp%user_func_ctx => func_ctx
      tmp%user_grad_ctx => grad_ctx
      allocate (tmp%user_context, source=context)

      if (present(tol)) tmp%tol = tol
      if (present(tolx)) tmp%tolx = tolx
      if (present(tol_relax_factor)) tmp%tol_relax_factor = tol_relax_factor
      if (present(alpha)) tmp%alpha = alpha
      if (present(use_broyden)) tmp%use_broyden = use_broyden
      if (present(max_retries)) tmp%max_retries = max_retries
      if (present(retry_radius)) tmp%retry_radius = retry_radius
      if (present(max_iter)) tmp%max_iter = max_iter
      if (present(max_roots)) tmp%max_roots = max_roots
      if (present(debug)) tmp%debug = debug

      if (tmp%tol_relax_factor <= 0.0_wp) then
         call fatal_error(error, "Newton-deflation: tol_relax_factor must be positive")
         return
      end if

      if (tmp%max_roots <= 0) then
         call fatal_error(error, "Newton-deflation: max_roots must be positive")
         return
      end if

      ! Bounds: default to "ignore" mode if none provided.
      if (present(bounds_mode)) tmp%bounds_mode = bounds_mode
      if (present(xlow)) then
         if (size(xlow) /= n) then
            call fatal_error(error, "Newton-deflation: xlow size /= n")
            return
         end if
         allocate (tmp%xl_init(n)); tmp%xl_init = xlow
      end if
      if (present(xupp)) then
         if (size(xupp) /= n) then
            call fatal_error(error, "Newton-deflation: xupp size /= n")
            return
         end if
         allocate (tmp%xu_init(n)); tmp%xu_init = xupp
      end if

      ! Ball cap setup: anchor + positive branch_rho_cut, plus bounds present.
      tmp%has_ball = .false.
      if (present(branch_rho_cut)) tmp%branch_rho_cut = branch_rho_cut
      if (present(anchor) .and. tmp%branch_rho_cut > 0.0_wp) then
         if (size(anchor) > n) then
            call fatal_error(error, "Newton-deflation: anchor size > n")
            return
         end if
         tmp%n_anchor = size(anchor)
         allocate (tmp%anchor(tmp%n_anchor))
         tmp%anchor = anchor
         tmp%has_ball = .true.
      end if

      call tmp%deflation%init(n_dim=n, max_roots=tmp%max_roots, error=solver_error, &
                              p_power=p_power, alpha_shift=alpha_shift, dedup_tol=dedup_tol)
      if (allocated(solver_error)) then
         call fatal_error(error, solver_error%message)
         return
      end if

      ! Build the initial inner Newton with the caller's bounds.
      call tmp%build_inner_newton(tmp%xl_init, tmp%xu_init, error)
      if (allocated(error)) return

      call move_alloc(tmp, solver)
   end subroutine new_newton_deflation_solver

   !> Construct (or reconstruct) the inner Newton solver with the given
   !> bounds. Tears down a previous instance if present. Used at factory
   !> time and again whenever the ball cap tightens after a root is found.
   subroutine newton_deflation_build_inner(self, xl, xu, error)
      class(moist_math_solver_newton_deflation_type), intent(inout), target :: self
      real(wp), dimension(:), intent(in), optional :: xl, xu
      type(error_type), allocatable, intent(out) :: error

      type(deflation_newton_context_type) :: wrapped_ctx
      type(error_type), allocatable :: solver_error
      real(wp) :: tol_use, tolx_use

      if (allocated(self%newton_solver)) then
         call self%newton_solver%destroy()
         deallocate (self%newton_solver)
      end if

      wrapped_ctx%user_func_ctx => self%user_func_ctx
      wrapped_ctx%user_grad_ctx => self%user_grad_ctx
      wrapped_ctx%deflation => self%deflation
      wrapped_ctx%user_context => self%user_context

      tol_use = self%tol*self%tol_relax_factor
      tolx_use = self%tolx*self%tol_relax_factor

      if (present(xl) .and. present(xu)) then
         call new_newton_solver( &
            solver=self%newton_solver, &
            n=self%n, m=self%m, error=solver_error, &
            func_ctx=deflated_residual, &
            grad_ctx=deflated_jacobian, &
            context=wrapped_ctx, &
            max_iter=self%max_iter, &
            tol=tol_use, tolx=tolx_use, &
            alpha=self%alpha, use_broyden=self%use_broyden, &
            bounds_mode=self%bounds_mode, &
            xlow=xl, xupp=xu, &
            verbose=.false.)
      else
         call new_newton_solver( &
            solver=self%newton_solver, &
            n=self%n, m=self%m, error=solver_error, &
            func_ctx=deflated_residual, &
            grad_ctx=deflated_jacobian, &
            context=wrapped_ctx, &
            max_iter=self%max_iter, &
            tol=tol_use, tolx=tolx_use, &
            alpha=self%alpha, use_broyden=self%use_broyden, &
            verbose=.false.)
      end if
      if (allocated(solver_error)) then
         call fatal_error(error, solver_error%message)
         return
      end if
   end subroutine newton_deflation_build_inner

   !> Solve the nonlinear system, enumerating up to `max_roots` distinct
   !> roots via iterated deflation. On exit `x` holds the first root
   !> discovered. The full list is available via get_raw_candidates.
   !>
   !> Seeding strategy (combined continuation + anchor fallback):
   !>   iter=1     -- try the un-perturbed caller seed (warm-started by the
   !>                 caller for the projection problem). If it fails, try
   !>                 up to `max_retries` perturbations of that seed.
   !>   iter >= 2  -- two-phase probe:
   !>                 phase A (continuation): perturbations of the most
   !>                   recently accepted root `x_base`. Cheap and works
   !>                   for asymmetric branch arrangements where the next
   !>                   sibling root is geometrically near the previous one.
   !>                 phase B (anchor fallback): only entered if phase A
   !>                   exhausts without finding a new root. Perturbations
   !>                   of the original (warm-started) caller seed `x_seed`.
   !>                   Required for radially symmetric arrangements where
   !>                   sibling roots sit equidistant *around* the anchor
   !>                   but far *from each other* (e.g. octahedral 4-fold
   !>                   anchors), so a small perturbation from any one root
   !>                   cannot reach the next sibling.
   !>
   !> Perturbations only touch the leading n_anchor components when a ball
   !> cap is active so lambda flows through (the previous root's lambda is
   !> a much better guess than 0 for the neighboring root's lambda).
   subroutine newton_deflation_solve(self, x, error)
      class(moist_math_solver_newton_deflation_type), intent(inout), target :: self
      real(wp), dimension(:), intent(inout) :: x
      type(error_type), allocatable, intent(out) :: error

      type(error_type), allocatable :: inner_error, build_error
      real(wp), dimension(size(x)) :: x_seed, x_trial, x_first, perturb, x_base
      real(wp), allocatable :: converged(:, :), xl_tight(:), xu_tight(:)
      integer :: iter, attempt, k, n_perturb, start_attempt
      logical :: accepted, ball_armed, first_root_found, found_new_root
      real(wp) :: phi_min, phi_max

      if (size(x) /= self%n) then
         call fatal_error(error, "Newton-deflation: x has wrong dimension")
         return
      end if
      if (.not. allocated(self%newton_solver)) then
         call fatal_error(error, "Newton-deflation: inner Newton not initialized")
         return
      end if

      call self%deflation%reset()

      allocate (converged(self%n, self%max_roots), source=0.0_wp)
      self%n_raw_candidates = 0
      ball_armed = .false.
      first_root_found = .false.
      ! Perturb only the leading n_anchor components (xyz) when has_ball;
      ! otherwise the whole vector. lambda is left at the seed value
      ! (typically 0 or a warm-started tangent-plane estimate from the
      ! caller) because the right multiplier is a function of xyz, not a
      ! free knob.
      n_perturb = self%n
      if (self%has_ball) n_perturb = self%n_anchor

      x_seed = x
      x_first = x
      ! x_base is the seed for the *current* outer iteration. It starts at
      ! the caller's (warm-started) seed and is updated to the most recent
      ! accepted root after each successful iter, so subsequent iterations
      ! probe the neighborhood of the previous basin rather than restarting
      ! from the anchor every time.
      x_base = x_seed

      outer: do iter = 1, self%max_roots
         found_new_root = .false.
         ! Phase A (continuation): perturbations of x_base = previous root.
         ! iter 1: try un-perturbed seed first (attempt=0), then perturbations.
         ! iter >= 2: skip attempt=0 - x_base IS the previous accepted root,
         ! so the un-perturbed solve would just rediscover it.
         start_attempt = 0
         if (iter > 1) start_attempt = 1

         attempts_cont: do attempt = start_attempt, self%max_retries
            if (attempt == 0) then
               x_trial = x_base
            else
               ! Globally unique offset index across (iter, attempt) so each
               ! iteration probes a fresh octant rather than recycling.
               k = (iter - 1)*self%max_retries + attempt
               perturb = 0.0_wp
               perturb(1:n_perturb) = retry_offset(n_perturb, k, self%retry_radius)
               x_trial = x_base + perturb
            end if

            call self%newton_solver%solve(x_trial, inner_error)
            if (allocated(inner_error)) then
               if (self%debug) then
                  write (output_unit, '(x,a,i0,a,i0,a,a)') &
                     '[newton-deflation] iter ', iter, ' contA attempt ', attempt, &
                     ' inner Newton failed: ', trim(inner_error%message)
               end if
               deallocate (inner_error)
               cycle attempts_cont
            end if

            call self%deflation%append_root(x_trial, accepted)
            if (accepted) then
               found_new_root = .true.
               exit attempts_cont
            end if
            if (self%debug) then
               write (output_unit, '(x,a,i0,a,i0,a)') &
                  '[newton-deflation] iter ', iter, ' contA attempt ', attempt, &
                  ' converged to a known root (try next perturbation)'
            end if
         end do attempts_cont

         ! Phase B (anchor fallback): only entered when phase A exhausted
         ! and the original anchor seed is genuinely different from x_base
         ! (i.e. iter >= 2). For iter=1 phase A already perturbed x_seed.
         if (.not. found_new_root .and. iter > 1) then
            attempts_anchor: do attempt = 1, self%max_retries
               ! Offset range disjoint from phase A so phase B probes
               ! different perturbations even after the cyclic axis pattern
               ! would have repeated.
               k = (iter - 1)*self%max_retries + attempt &
                   + self%max_roots*self%max_retries
               perturb = 0.0_wp
               perturb(1:n_perturb) = retry_offset(n_perturb, k, self%retry_radius)
               x_trial = x_seed + perturb

               call self%newton_solver%solve(x_trial, inner_error)
               if (allocated(inner_error)) then
                  if (self%debug) then
                     write (output_unit, '(x,a,i0,a,i0,a,a)') &
                        '[newton-deflation] iter ', iter, ' anchorB attempt ', &
                        attempt, ' inner Newton failed: ', trim(inner_error%message)
                  end if
                  deallocate (inner_error)
                  cycle attempts_anchor
               end if

               call self%deflation%append_root(x_trial, accepted)
               if (accepted) then
                  found_new_root = .true.
                  exit attempts_anchor
               end if
               if (self%debug) then
                  write (output_unit, '(x,a,i0,a,i0,a)') &
                     '[newton-deflation] iter ', iter, ' anchorB attempt ', &
                     attempt, ' converged to a known root (try next perturbation)'
               end if
            end do attempts_anchor
         end if

         if (.not. found_new_root) then
            if (self%debug) then
               write (output_unit, '(x,a,i0,a)') &
                  '[newton-deflation] iter ', iter, ' exhausted all attempts (stop)'
            end if
            exit outer
         end if

         self%n_raw_candidates = self%n_raw_candidates + 1
         converged(:, self%n_raw_candidates) = x_trial
         ! Continuation: next iter's perturbations are around the just-found
         ! root, not the original anchor seed.
         x_base = x_trial
         if (self%n_raw_candidates == 1) then
            x_first = x_trial
            first_root_found = .true.

            ! After the first accepted root, tighten Newton's bounds to
            ! anchor +/- (phi_min + branch_rho_cut) on the n_anchor leading
            ! components. lambda (and any other tail entries) keep their
            ! initial bounds.
            if (self%has_ball .and. .not. ball_armed) then
               phi_min = norm2(x_trial(1:self%n_anchor) - self%anchor)
               phi_max = phi_min + self%branch_rho_cut
               allocate (xl_tight(self%n))
               allocate (xu_tight(self%n))
               if (allocated(self%xl_init)) then
                  xl_tight = self%xl_init
               else
                  xl_tight = -huge(1.0_wp)
               end if
               if (allocated(self%xu_init)) then
                  xu_tight = self%xu_init
               else
                  xu_tight = huge(1.0_wp)
               end if
               xl_tight(1:self%n_anchor) = max(xl_tight(1:self%n_anchor), &
                                               self%anchor - phi_max)
               xu_tight(1:self%n_anchor) = min(xu_tight(1:self%n_anchor), &
                                               self%anchor + phi_max)
               call self%build_inner_newton(xl_tight, xu_tight, build_error)
               deallocate (xl_tight, xu_tight)
               if (allocated(build_error)) then
                  if (self%debug) then
                     write (output_unit, '(x,a,a)') &
                        '[newton-deflation] bounds rebuild failed: ', &
                        trim(build_error%message)
                  end if
                  ! Non-fatal: stop the deflation search rather than aborting
                  ! the whole solve.
                  deallocate (build_error)
                  exit outer
               end if
               ball_armed = .true.
               if (self%debug) then
                  write (output_unit, '(x,a,es12.4)') &
                     '[newton-deflation] ball cap (phi_max) = ', phi_max
               end if
            end if
         end if

         if (self%debug) then
            write (output_unit, '(x,a,i0,a)') &
               '[newton-deflation] iter ', iter, ' accepted root'
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
         call fatal_error(error, "Newton-deflation: no roots enumerated")
         deallocate (converged)
         return
      end if

      x = x_first
      deallocate (converged)
   end subroutine newton_deflation_solve

   !> Retry-offset pattern: cycle through axis-aligned +/- directions for
   !> n=3 (the common xyz case), deterministic spiral otherwise. Mirrors
   !> the SLSQP-deflation helper so both deflation paths sample the same
   !> octants on the projection problem.
   pure function retry_offset(n, k, r) result(off)
      integer, intent(in) :: n
      integer, intent(in) :: k
      real(wp), intent(in) :: r
      real(wp), dimension(n) :: off

      integer :: axis, sign_k, i, kk
      real(wp) :: u

      off = 0.0_wp
      if (n == 3) then
         axis = mod(k - 1, 3) + 1
         sign_k = 1
         if (mod((k - 1)/3, 2) == 1) sign_k = -1
         off(axis) = real(sign_k, wp)*r
      else
         kk = k
         do i = 1, n
            u = sin(real(kk*i, wp))
            off(i) = r*u/sqrt(real(n, wp))
         end do
      end if
   end function retry_offset

   !> Return the full list of converged roots.
   subroutine newton_deflation_get_raw_candidates(self, candidates, n_candidates)
      class(moist_math_solver_newton_deflation_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: candidates(:, :)
      integer, intent(out) :: n_candidates

      n_candidates = self%n_raw_candidates
      if (n_candidates <= 0 .or. .not. allocated(self%raw_candidates)) then
         allocate (candidates(self%n, 0))
         return
      end if
      allocate (candidates(self%n, n_candidates))
      candidates(:, :) = self%raw_candidates(:, :)
   end subroutine newton_deflation_get_raw_candidates

   !> Release resources.
   subroutine newton_deflation_destroy(self)
      class(moist_math_solver_newton_deflation_type), intent(inout), target :: self

      if (allocated(self%newton_solver)) then
         call self%newton_solver%destroy()
         deallocate (self%newton_solver)
      end if
      call self%deflation%destroy()
      if (allocated(self%raw_candidates)) deallocate (self%raw_candidates)
      if (allocated(self%user_context)) deallocate (self%user_context)
      if (allocated(self%xl_init)) deallocate (self%xl_init)
      if (allocated(self%xu_init)) deallocate (self%xu_init)
      if (allocated(self%anchor)) deallocate (self%anchor)
      self%has_ball = .false.
      self%n_anchor = 0
      self%user_func_ctx => null()
      self%user_grad_ctx => null()
      self%n = 0; self%m = 0
      self%n_raw_candidates = 0
   end subroutine newton_deflation_destroy

   !>==================================================================
   !> Callback wrappers installed on the inner Newton solver.
   !>==================================================================

   !> Deflated residual: F_def(x) = M(x) * F(x).
   subroutine deflated_residual(x, f, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: f
      class(*), intent(in) :: context

      real(wp), dimension(size(f)) :: f_raw
      real(wp) :: m_val

      f = 0.0_wp
      select type (ctx => context)
      type is (deflation_newton_context_type)
         if (.not. associated(ctx%user_func_ctx)) return
         if (.not. associated(ctx%deflation)) return
         if (.not. associated(ctx%user_context)) return
         call ctx%user_func_ctx(x, f_raw, ctx%user_context)
         m_val = ctx%deflation%multiplier(x)
         f = m_val*f_raw
      end select
   end subroutine deflated_residual

   !> Deflated Jacobian: J_def(x) = M*J + F * grad_M^T (rank-1 outer product).
   subroutine deflated_jacobian(x, jac, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:, :), intent(out) :: jac
      class(*), intent(in) :: context

      real(wp), dimension(size(jac, 1)) :: f_raw
      real(wp), dimension(size(jac, 1), size(jac, 2)) :: jac_raw
      real(wp), dimension(size(x)) :: grad_m
      real(wp) :: m_val
      integer :: i, j

      jac = 0.0_wp
      select type (ctx => context)
      type is (deflation_newton_context_type)
         if (.not. associated(ctx%user_func_ctx)) return
         if (.not. associated(ctx%user_grad_ctx)) return
         if (.not. associated(ctx%deflation)) return
         if (.not. associated(ctx%user_context)) return
         call ctx%user_func_ctx(x, f_raw, ctx%user_context)
         call ctx%user_grad_ctx(x, jac_raw, ctx%user_context)
         m_val = ctx%deflation%multiplier(x)
         call ctx%deflation%gradient(x, grad_m)
         do j = 1, size(jac, 2)
            do i = 1, size(jac, 1)
               jac(i, j) = m_val*jac_raw(i, j) + f_raw(i)*grad_m(j)
            end do
         end do
      end select
   end subroutine deflated_jacobian

end module moist_math_solver_newton_deflation
