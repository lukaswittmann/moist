!> Multi-tangent SLSQP solver using tangent-plane ring seeds
!>
!> Builds a tangent basis at the anchor from the owner-anchor direction,
!> generates ring seeds in that tangent plane for configurable radii and
!> points-per-ring, runs SLSQP from each seed, and returns the solution
!> closest to the anchor. Uses the same objective/constraint callbacks as
!> the regular SLSQP setup.
module moist_math_solver_slsqp_multi_tangent
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use iso_fortran_env, only: output_unit
   use moist_type, only: solver_base_type
   use moist_math_solver_slsqp, only: new_slsqp_solver
   implicit none
   private

   public :: moist_math_solver_slsqp_multi_tangent_type
   public :: new_slsqp_multi_tangent_solver

   integer, parameter :: default_n_layers = 3
   real(wp), parameter :: zero_tol = 1.0e-12_wp
   real(wp), parameter :: default_radii(default_n_layers) = [0.2_wp, 0.5_wp, 0.8_wp]
   integer, parameter :: default_n_points(default_n_layers) = [6, 14, 26]

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

   !> Multi-tangent SLSQP solver
   type, extends(solver_base_type) :: moist_math_solver_slsqp_multi_tangent_type
      private
      real(wp) :: anchor(3)
      real(wp) :: owner_pos(3)
      real(wp), allocatable :: seeds(:, :)
      real(wp), allocatable :: raw_candidates(:, :)
      integer :: n_seeds = 0
      integer :: n_raw_candidates = 0
      class(solver_base_type), allocatable :: slsqp_solver
      procedure(objective_context_interface), pointer, nopass :: user_obj_ctx => null()
      procedure(objective_grad_context_interface), pointer, nopass :: user_obj_grad_ctx => null()
      procedure(constraints_context_interface), pointer, nopass :: user_con_ctx => null()
      procedure(constraints_grad_context_interface), pointer, nopass :: user_con_grad_ctx => null()
      procedure(iteration_callback_context_interface), pointer, nopass :: user_iter_callback_ctx => null()
      class(*), allocatable :: user_context
      integer :: n = 0
      integer :: m = 0
      integer :: meq = 0
      real(wp), allocatable :: xl(:)
      real(wp), allocatable :: xu(:)
      real(wp), allocatable :: radii(:)
      integer, allocatable :: n_points(:)
      integer :: n_layers = 0
      real(wp) :: tol = 1.0e-8_wp
      real(wp) :: toldx = 1.0e-8_wp
      real(wp) :: toldf = 1.0e-8_wp
      integer :: max_iter = 50
      logical :: debug = .false.
   contains
      procedure :: solve => slsqp_multi_tangent_solve
      procedure :: get_raw_candidates => slsqp_multi_tangent_get_raw_candidates
      procedure :: destroy => slsqp_multi_tangent_destroy
   end type moist_math_solver_slsqp_multi_tangent_type

contains

   !> Build an orthonormal tangent-plane basis for a given unit normal.
   !> @param[in]  normal  Unit normal vector (3)
   !> @param[out] t1      First tangent vector (3)
   !> @param[out] t2      Second tangent vector (3), equals normal x t1
   pure subroutine build_tangent_basis(normal, t1, t2)
      real(wp), intent(in) :: normal(3)
      real(wp), intent(out) :: t1(3), t2(3)

      real(wp) :: ax(3), dot_n, inv_norm
      integer :: k

      k = 1
      if (abs(normal(2)) < abs(normal(k))) k = 2
      if (abs(normal(3)) < abs(normal(k))) k = 3

      ax = 0.0_wp
      ax(k) = 1.0_wp

      dot_n = ax(1)*normal(1) + ax(2)*normal(2) + ax(3)*normal(3)
      t1 = ax - dot_n*normal
      inv_norm = 1.0_wp/sqrt(t1(1)**2 + t1(2)**2 + t1(3)**2)
      t1 = t1*inv_norm

      t2(1) = normal(2)*t1(3) - normal(3)*t1(2)
      t2(2) = normal(3)*t1(1) - normal(1)*t1(3)
      t2(3) = normal(1)*t1(2) - normal(2)*t1(1)
   end subroutine build_tangent_basis

   !> Generate seed points on a single tangent-plane ring around the anchor.
   !> @param[in]    anchor    Centre of the ring
   !> @param[in]    t1        First tangent basis vector
   !> @param[in]    t2        Second tangent basis vector
   !> @param[in]    radius    Ring radius
   !> @param[in]    n_points  Number of points on this ring
   !> @param[inout] seeds     Seed array to populate (3, total_seeds)
   !> @param[inout] offset    Current write offset; advanced by n_points on exit
   !> @param[out]   error    Error status
   subroutine generate_layer_seeds(anchor, t1, t2, radius, n_points, seeds, offset, error)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: t1(3), t2(3)
      real(wp), intent(in) :: radius
      integer, intent(in) :: n_points
      real(wp), intent(inout) :: seeds(:, :)
      integer, intent(inout) :: offset
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: two_pi = 2.0_wp*acos(-1.0_wp)
      real(wp) :: theta
      integer :: i

      if (n_points <= 0) then
         call fatal_error(error, "Multi-tangent SLSQP solver: points-per-ring must be > 0")
         return
      end if

      do i = 1, n_points
         theta = two_pi*real(i - 1, wp)/real(n_points, wp)
         seeds(:, offset + i) = anchor + radius*(cos(theta)*t1 + sin(theta)*t2)
      end do
      offset = offset + n_points
   end subroutine generate_layer_seeds

   !> Factory function to create a multi-tangent SLSQP solver.
   subroutine new_slsqp_multi_tangent_solver(anchor, owner_pos, solver, &
                                             n, m, meq, obj_ctx, obj_grad_ctx, con_ctx, con_grad_ctx, context, &
                                             xl, xu, max_iter, tol, toldx, toldf, verbose, iter_callback_ctx, &
                                             radii, n_points, debug, error)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: owner_pos(3)
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
      integer, intent(in), optional :: max_iter
      real(wp), intent(in), optional :: tol
      real(wp), intent(in), optional :: toldx
      real(wp), intent(in), optional :: toldf
      logical, intent(in), optional :: verbose
      procedure(iteration_callback_context_interface), optional :: iter_callback_ctx
      real(wp), intent(in), optional :: radii(:)
      integer, intent(in), optional :: n_points(:)
      logical, intent(in), optional :: debug
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_solver_slsqp_multi_tangent_type), allocatable :: tmp
      type(error_type), allocatable :: solver_error
      integer :: offset, total_seeds, ilayer
      real(wp) :: normal(3), norm_normal, t1(3), t2(3)
      logical :: verbose_use

      allocate (tmp)
      tmp%anchor = anchor
      tmp%owner_pos = owner_pos
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
      if (present(debug)) tmp%debug = debug

      if (present(radii) .neqv. present(n_points)) then
         call fatal_error(error, &
                          "Multi-tangent SLSQP solver: radii and n_points must be provided together")
         return
      end if
      if (present(radii)) then
         if (size(radii) <= 0) then
            call fatal_error(error, "Multi-tangent SLSQP solver: radii array must be non-empty")
            return
         end if
         if (size(n_points) <= 0) then
            call fatal_error(error, "Multi-tangent SLSQP solver: n_points array must be non-empty")
            return
         end if
         if (size(radii) /= size(n_points)) then
            call fatal_error(error, &
                             "Multi-tangent SLSQP solver: radii and n_points must have same size")
            return
         end if
         tmp%n_layers = size(radii)
         allocate (tmp%radii(tmp%n_layers), tmp%n_points(tmp%n_layers))
         tmp%radii = radii
         tmp%n_points = n_points
      else
         tmp%n_layers = default_n_layers
         allocate (tmp%radii(tmp%n_layers), tmp%n_points(tmp%n_layers))
         tmp%radii = default_radii
         tmp%n_points = default_n_points
      end if

      if (tmp%n /= 3) then
         call fatal_error(error, "Multi-tangent SLSQP solver: n must be 3 for tangent seeding")
         return
      end if

      normal = tmp%anchor - tmp%owner_pos
      norm_normal = sqrt(sum(normal**2))
      if (norm_normal <= zero_tol) then
         call fatal_error(error, &
                          "Multi-tangent SLSQP solver: owner position coincides with anchor")
         return
      end if
      normal = normal/norm_normal
      call build_tangent_basis(normal, t1, t2)

      ! Count active layers (skip layers with radius <= 0)
      total_seeds = 0
      do ilayer = 1, tmp%n_layers
         if (tmp%radii(ilayer) > 0.0_wp) then
            if (tmp%n_points(ilayer) <= 0) then
               call fatal_error(error, &
                                "Multi-tangent SLSQP solver: n_points must be > 0 for active radii")
               return
            end if
            total_seeds = total_seeds + tmp%n_points(ilayer)
         end if
      end do
      if (total_seeds <= 0) then
         call fatal_error(error, "Multi-tangent SLSQP solver: no active seed layers")
         return
      end if
      tmp%n_seeds = total_seeds
      allocate (tmp%seeds(3, total_seeds))

      offset = 0

      ! Generate seed points on each active tangent-plane ring
      do ilayer = 1, tmp%n_layers
         if (tmp%radii(ilayer) > 0.0_wp) then
            call generate_layer_seeds(tmp%anchor, t1, t2, tmp%radii(ilayer), tmp%n_points(ilayer), &
                                      tmp%seeds, offset, solver_error)
            if (allocated(solver_error)) then
               call fatal_error(error, solver_error%message)
               return
            end if
         end if
      end do

      verbose_use = .false.
      if (present(verbose)) verbose_use = verbose

      if (present(iter_callback_ctx)) then
         call new_slsqp_solver( &
            solver=tmp%slsqp_solver, &
            n=tmp%n, m=tmp%m, meq=tmp%meq, &
            error=solver_error, &
            obj_ctx=tmp%user_obj_ctx, &
            obj_grad_ctx=tmp%user_obj_grad_ctx, &
            con_ctx=tmp%user_con_ctx, &
            con_grad_ctx=tmp%user_con_grad_ctx, &
            iter_callback_ctx=tmp%user_iter_callback_ctx, &
            context=tmp%user_context, &
            max_iter=tmp%max_iter, &
            tol=tmp%tol, &
            toldx=tmp%toldx, &
            toldf=tmp%toldf, &
            xl=tmp%xl, &
            xu=tmp%xu, &
            verbose=verbose_use, &
            iunit=output_unit &
            )
      else
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
            xu=tmp%xu, &
            verbose=verbose_use, &
            iunit=output_unit &
            )
      end if

      if (allocated(solver_error)) then
         call fatal_error(error, solver_error%message)
         return
      end if

      call move_alloc(tmp, solver)
   end subroutine new_slsqp_multi_tangent_solver

   !> Solve the constrained projection using multi-tangent SLSQP.
   subroutine slsqp_multi_tangent_solve(self, x, error)
      class(moist_math_solver_slsqp_multi_tangent_type), intent(inout), target :: self
      real(wp), dimension(:), intent(inout) :: x
      type(error_type), allocatable, intent(out) :: error

      type(error_type), allocatable :: solver_error
      real(wp) :: x_trial(3), best_x(3), best_dist2, dist2
      integer :: i, n_converged
      real(wp), allocatable :: converged(:, :)

      if (size(x) /= 3) then
         call fatal_error(error, "Multi-tangent SLSQP solver: x must have size 3")
         return
      end if
      if (.not. allocated(self%seeds)) then
         call fatal_error(error, "Multi-tangent SLSQP solver: seeds not initialized")
         return
      end if
      if (self%n_seeds <= 0 .or. size(self%seeds, dim=2) /= self%n_seeds) then
         call fatal_error(error, "Multi-tangent SLSQP solver: invalid seed count")
         return
      end if
      if (.not. allocated(self%slsqp_solver)) then
         call fatal_error(error, "Multi-tangent SLSQP solver: SLSQP not initialized")
         return
      end if

      best_dist2 = huge(1.0_wp)
      best_x = self%anchor

      if (self%debug) then
         write (output_unit, '(x,a)') &
            '========== Multi-tangent SLSQP ========='
      end if

      allocate (converged(3, self%n_seeds))
      n_converged = 0

      do i = 1, self%n_seeds
         x_trial = self%seeds(:, i)

         call self%slsqp_solver%solve(x_trial, solver_error)
         if (allocated(solver_error)) then
            if (self%debug) then
               write (output_unit, '(x,a,i0,a,a)') 'Seed ', i, ' failed: ', trim(solver_error%message)
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
         call fatal_error(error, "Multi-tangent SLSQP solver: no successful starts")
         return
      end if

      if (allocated(self%raw_candidates)) deallocate (self%raw_candidates)
      self%n_raw_candidates = n_converged
      allocate (self%raw_candidates(3, n_converged))
      self%raw_candidates(:, :) = converged(:, 1:n_converged)

      x = best_x

      deallocate (converged)
   end subroutine slsqp_multi_tangent_solve

   !> Get raw SLSQP candidates converged from multi-start seeds.
   !> @param[out] candidates  Raw candidate points (3,n)
   !> @param[out] n_candidates Number of available candidates
   subroutine slsqp_multi_tangent_get_raw_candidates(self, candidates, n_candidates)
      class(moist_math_solver_slsqp_multi_tangent_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: candidates(:, :)
      integer, intent(out) :: n_candidates

      n_candidates = self%n_raw_candidates
      if (n_candidates <= 0 .or. .not. allocated(self%raw_candidates)) then
         allocate (candidates(3, 0))
         return
      end if
      allocate (candidates(3, n_candidates))
      candidates(:, :) = self%raw_candidates(:, :)
   end subroutine slsqp_multi_tangent_get_raw_candidates

   !> Clean up resources
   subroutine slsqp_multi_tangent_destroy(self)
      class(moist_math_solver_slsqp_multi_tangent_type), intent(inout), target :: self
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
      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%n_points)) deallocate (self%n_points)
      self%user_obj_ctx => null()
      self%user_obj_grad_ctx => null()
      self%user_con_ctx => null()
      self%user_con_grad_ctx => null()
      self%user_iter_callback_ctx => null()
   end subroutine slsqp_multi_tangent_destroy

end module moist_math_solver_slsqp_multi_tangent
