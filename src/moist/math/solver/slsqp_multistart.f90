!> Multi-start SLSQP solver using a small Lebedev seed cloud
!>
!> Generates a set of initial guesses around the anchor using configurable
!> Lebedev shells, runs SLSQP from each seed, and returns the
!> solution closest to the anchor. Uses the same objective/constraint
!> callbacks as the regular SLSQP setup.
module moist_math_solver_slsqp_multistart
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use iso_fortran_env, only: output_unit
   use moist_type, only: solver_base_type
   use moist_math_solver_slsqp, only: new_slsqp_solver
   use moist_math_grid_lebedev, only: lebedev_order_from_num, get_angular_grid
   use moist_math_trigonometry, only: rotation_z_to_n
   implicit none
   private

   public :: moist_math_solver_slsqp_multistart_type
   public :: new_slsqp_multistart_solver

   integer, parameter :: default_n_layers = 3
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

   !> Multi-start SLSQP solver
   type, extends(solver_base_type) :: moist_math_solver_slsqp_multistart_type
      private
      real(wp) :: anchor(3)
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
      !> Maximum number of retry attempts with expanded radius when no seeds converge
      integer :: max_retries = 3
      !> Radius increment (bohr) for each retry attempt
      real(wp) :: radius_increment = 1.0_wp
      !> Rotation matrix aligning Lebedev grid to surface normal (identity if unset)
      real(wp) :: rot(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
      logical :: debug = .false.
   contains
      procedure :: solve => slsqp_multistart_solve
      procedure :: get_raw_candidates => slsqp_multistart_get_raw_candidates
      procedure :: destroy => slsqp_multistart_destroy
   end type moist_math_solver_slsqp_multistart_type

contains

   !> Generate seed points on a single Lebedev shell around the anchor.
   !> @param[in]    anchor   Centre of the shell
   !> @param[in]    radius   Shell radius
   !> @param[in]    num_leb  Number of Lebedev points on this shell
   !> @param[inout] seeds    Seed array to populate (3, total_seeds)
   !> @param[inout] offset   Current write offset; advanced by num_leb on exit
   !> @param[in]    rot      Optional 3x3 rotation matrix applied to grid points
   !> @param[out]   error    Error status
   subroutine generate_layer_seeds(anchor, radius, num_leb, seeds, offset, error, rot)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: radius
      integer, intent(in) :: num_leb
      real(wp), intent(inout) :: seeds(:, :)
      integer, intent(inout) :: offset
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: rot(3, 3)

      real(wp), allocatable :: grid(:, :), weights(:)
      real(wp) :: pt(3)
      integer :: order, i

      call lebedev_order_from_num(num_leb, order, error)
      if (allocated(error)) return

      allocate (grid(3, num_leb), weights(num_leb))
      call get_angular_grid(order, grid, weights, error)
      if (allocated(error)) return

      do i = 1, num_leb
         if (present(rot)) then
            pt = matmul(rot, grid(:, i))
         else
            pt = grid(:, i)
         end if
         seeds(:, offset + i) = anchor + radius*pt
      end do
      offset = offset + num_leb
   end subroutine generate_layer_seeds

   !> Factory function to create a multi-start SLSQP solver.
   subroutine new_slsqp_multistart_solver(anchor, solver, &
                                          n, m, meq, obj_ctx, obj_grad_ctx, con_ctx, con_grad_ctx, context, &
                                          xl, xu, max_iter, tol, toldx, toldf, verbose, iter_callback_ctx, &
                                          radii, n_points, debug, max_retries, radius_increment, normal, error)
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
      integer, intent(in), optional :: max_iter
      real(wp), intent(in), optional :: tol
      real(wp), intent(in), optional :: toldx
      real(wp), intent(in), optional :: toldf
      logical, intent(in), optional :: verbose
      procedure(iteration_callback_context_interface), optional :: iter_callback_ctx
      real(wp), intent(in), optional :: radii(:)
      integer, intent(in), optional :: n_points(:)
      logical, intent(in), optional :: debug
      !> Maximum number of retry attempts with expanded radius (default: 3)
      integer, intent(in), optional :: max_retries
      !> Radius increment per retry attempt in bohr (default: 1.0)
      real(wp), intent(in), optional :: radius_increment
      !> Surface normal vector (atom -> anchor); orients the Lebedev grid
      real(wp), intent(in), optional :: normal(3)
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_solver_slsqp_multistart_type), allocatable :: tmp
      type(error_type), allocatable :: solver_error
      integer :: offset, total_seeds, ilayer
      logical :: verbose_use

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
      if (present(debug)) tmp%debug = debug
      if (present(max_retries)) tmp%max_retries = max_retries
      if (present(radius_increment)) tmp%radius_increment = radius_increment

      ! Build rotation matrix from surface normal
      if (present(normal)) then
         if (norm2(normal) > 0.0_wp) then
            call rotation_z_to_n(normal/norm2(normal), tmp%rot)
         end if
      end if

      if (present(radii) .neqv. present(n_points)) then
         call fatal_error(error, "Multi-start SLSQP solver: radii and n_points must be provided together")
         return
      end if
      if (present(radii)) then
         if (size(radii) <= 0) then
            call fatal_error(error, "Multi-start SLSQP solver: radii array must be non-empty")
            return
         end if
         if (size(n_points) <= 0) then
            call fatal_error(error, "Multi-start SLSQP solver: n_points array must be non-empty")
            return
         end if
         if (size(radii) /= size(n_points)) then
            call fatal_error(error, "Multi-start SLSQP solver: radii and n_points must have same size")
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
         call fatal_error(error, "Multi-start SLSQP solver: n must be 3")
         return
      end if

      ! Count active layers (skip layers with radius <= 0)
      total_seeds = 0
      do ilayer = 1, tmp%n_layers
         if (tmp%radii(ilayer) > 0.0_wp) then
            if (tmp%n_points(ilayer) <= 0) then
               call fatal_error(error, &
                                "Multi-start SLSQP solver: n_points must be > 0 for active radii")
               return
            end if
            total_seeds = total_seeds + tmp%n_points(ilayer)
         end if
      end do
      if (total_seeds <= 0) then
         call fatal_error(error, "Multi-start SLSQP solver: no active seed layers")
         return
      end if
      tmp%n_seeds = total_seeds
      allocate (tmp%seeds(3, total_seeds))

      offset = 0

      ! Generate seed points on each active Lebedev shell
      do ilayer = 1, tmp%n_layers
         if (tmp%radii(ilayer) > 0.0_wp) then
            call generate_layer_seeds(tmp%anchor, tmp%radii(ilayer), tmp%n_points(ilayer), &
                                      tmp%seeds, offset, solver_error, rot=tmp%rot)
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
   end subroutine new_slsqp_multistart_solver

   !> Solve the constrained projection using multi-start SLSQP.
   subroutine slsqp_multistart_solve(self, x, error)
      class(moist_math_solver_slsqp_multistart_type), intent(inout), target :: self
      real(wp), dimension(:), intent(inout) :: x
      type(error_type), allocatable, intent(out) :: error

      type(error_type), allocatable :: solver_error, layer_error
      real(wp) :: x_trial(3), best_x(3), best_dist2, dist2
      integer :: i, n_converged, iretry, retry_npts, retry_offset
      real(wp) :: retry_radius
      real(wp), allocatable :: converged(:, :), retry_seeds(:, :)

      if (size(x) /= 3) then
         call fatal_error(error, "Multi-start SLSQP solver: x must have size 3")
         return
      end if
      if (.not. allocated(self%seeds)) then
         call fatal_error(error, "Multi-start SLSQP solver: seeds not initialized")
         return
      end if
      if (self%n_seeds <= 0 .or. size(self%seeds, dim=2) /= self%n_seeds) then
         call fatal_error(error, "Multi-start SLSQP solver: invalid seed count")
         return
      end if
      if (.not. allocated(self%slsqp_solver)) then
         call fatal_error(error, "Multi-start SLSQP solver: SLSQP not initialized")
         return
      end if

      best_dist2 = huge(1.0_wp)
      best_x = self%anchor

      if (self%debug) then
         write (output_unit, '(x,a)') &
            '========== Multi-start SLSQP ========='
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

      ! Retry with progressively larger radii if no seeds converged
      if (n_converged == 0 .and. self%max_retries > 0) then
         retry_npts = maxval(self%n_points)

         do iretry = 1, self%max_retries
            retry_radius = maxval(self%radii) + iretry*self%radius_increment

            if (self%debug) then
               write (output_unit, '(x,a,i0,a,f8.3)') &
                  'Retry ', iretry, ': expanding radius to ', retry_radius
            end if

            allocate (retry_seeds(3, retry_npts))
            retry_offset = 0
            call generate_layer_seeds(self%anchor, retry_radius, retry_npts, &
                                      retry_seeds, retry_offset, layer_error, rot=self%rot)
            if (allocated(layer_error)) then
               deallocate (retry_seeds)
               deallocate (layer_error)
               cycle
            end if

            deallocate (converged)
            allocate (converged(3, retry_npts))

            do i = 1, retry_npts
               x_trial = retry_seeds(:, i)

               call self%slsqp_solver%solve(x_trial, solver_error)
               if (allocated(solver_error)) then
                  if (self%debug) then
                     write (output_unit, '(x,a,i0,a,i0,a,a)') &
                        'Retry ', iretry, ' seed ', i, ' failed: ', &
                        trim(solver_error%message)
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

            deallocate (retry_seeds)
            if (n_converged > 0) exit
         end do
      end if

      if (best_dist2 >= huge(1.0_wp)*0.5_wp) then
         call fatal_error(error, "Multi-start SLSQP solver: no successful starts")
         return
      end if

      if (allocated(self%raw_candidates)) deallocate (self%raw_candidates)
      self%n_raw_candidates = n_converged
      allocate (self%raw_candidates(3, n_converged))
      self%raw_candidates(:, :) = converged(:, 1:n_converged)

      x = best_x

      deallocate (converged)
   end subroutine slsqp_multistart_solve

   !> Get raw SLSQP candidates converged from multi-start seeds.
   !> @param[out] candidates  Raw candidate points (3,n)
   !> @param[out] n_candidates Number of available candidates
   subroutine slsqp_multistart_get_raw_candidates(self, candidates, n_candidates)
      class(moist_math_solver_slsqp_multistart_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: candidates(:, :)
      integer, intent(out) :: n_candidates

      n_candidates = self%n_raw_candidates
      if (n_candidates <= 0 .or. .not. allocated(self%raw_candidates)) then
         allocate (candidates(3, 0))
         return
      end if
      allocate (candidates(3, n_candidates))
      candidates(:, :) = self%raw_candidates(:, :)
   end subroutine slsqp_multistart_get_raw_candidates

   !> Clean up resources
   subroutine slsqp_multistart_destroy(self)
      class(moist_math_solver_slsqp_multistart_type), intent(inout), target :: self
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
   end subroutine slsqp_multistart_destroy

end module moist_math_solver_slsqp_multistart
