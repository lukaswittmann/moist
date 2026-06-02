!> Hierarchical Lebedev onion-shell global search solver for LSF surface
!>
!> Implements a hierarchical bisection algorithm that couples radial and angular
!> refinement to find points where LSF ~= 0, starting from an anchor point.
!> The Lebedev angular grid resolution is interpolated based on the current
!> search radius-coarse grids for large radii, fine grids near the surface.
!>
!> Algorithm:
!>   1. Start from small radius (near anchor) with coarse step, scan outward
!>   2. Lebedev grid size interpolated from radius (194 near anchor -> 26 far away)
!>   3. Stop immediately when any positive LSF detected at radius R_upper
!>   4. R_lower = last radius with all negative LSF, R_upper = first with any positive
!>   5. Bisect: halve step, refine grid, re-scan [R_lower, R_upper]
!>   6. Repeat until R_upper - R_lower < min_step or convergence
!>   7. Return best point found within converged bracket
!>
!> Uses screened LSF evaluation for O(N) complexity with many atoms.
module moist_math_solver_onion
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_type, only: solver_base_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_math_grid_lebedev, only: get_angular_grid, lebedev_order_from_num
   use iso_fortran_env, only: output_unit
   implicit none
   private

   public :: moist_math_solver_onion_type
   public :: new_onion_solver

   !> Onion-shell global search solver
   type, extends(solver_base_type) :: moist_math_solver_onion_type
      private
      !> LSF primitive for function evaluation
      class(moist_cavity_drop_lsf_type), pointer :: lsf => null()
      !> Anchor point (starting position for search)
      real(wp) :: anchor(3)
      !> Initial Lebedev order (number of angular points)
      integer :: initial_leb_num
      !> Maximum Lebedev order for refinement
      integer :: max_leb_num
      !> Initial radial step size (bohr)
      real(wp) :: initial_step
      !> Minimum radial step size (bohr)
      real(wp) :: min_step
      !> Maximum number of radial steps per direction
      integer :: max_steps
      !> LSF convergence tolerance
      real(wp) :: lsf_tol
      !> Step tolerance  (bohr)
      real(wp) :: step_tol
      !> Maximum search radius (bohr)
      real(wp) :: max_radius
      !> Debug flag
      logical :: debug
   contains
      !> Solve the global search problem
      procedure :: solve => onion_solve
      !> Clean up resources
      procedure :: destroy => onion_destroy
   end type moist_math_solver_onion_type

contains

   !> Factory function to create and initialize an onion solver
   !>
   !> Creates a hierarchical Lebedev global search solver for finding points
   !> on the LSF=0 surface starting from an anchor point.
   !>
   !> @param[in]  lsf           LSF primitive for function evaluation
   !> @param[in]  anchor         Starting point for radial search [3]
   !> @param[out] solver         Allocated solver instance
   !> @param[in]  initial_leb_num Initial Lebedev grid size (default: 26)
   !> @param[in]  max_leb_num    Maximum Lebedev grid size (default: 194)
   !> @param[in]  initial_step   Initial radial step size in bohr (default: 0.5)
   !> @param[in]  min_step       Minimum radial step size in bohr (default: 0.01)
   !> @param[in]  max_steps      Maximum radial steps per direction (default: 100)
   !> @param[in]  lsf_tol       LSF convergence tolerance (default: 1e-6)
   !> @param[in]  max_radius     Maximum search radius in bohr (default: 20.0)
   !> @param[in]  verbose        Verbosity level (default: 0)
   subroutine new_onion_solver(lsf, anchor, solver, &
                               initial_leb_num, max_leb_num, initial_step, min_step, &
                               max_steps, lsf_tol, step_tol, max_radius, debug)
      !> LSF primitive
      class(moist_cavity_drop_lsf_type), target, intent(inout) :: lsf
      !> Anchor point
      real(wp), intent(in) :: anchor(3)
      !> Output solver
      class(solver_base_type), allocatable, intent(out) :: solver
      !> Initial Lebedev grid size
      integer, intent(in), optional :: initial_leb_num
      !> Maximum Lebedev grid size
      integer, intent(in), optional :: max_leb_num
      !> Initial radial step size
      real(wp), intent(in), optional :: initial_step
      !> Minimum radial step size
      real(wp), intent(in), optional :: min_step
      !> Maximum radial steps
      integer, intent(in), optional :: max_steps
      !> LSF tolerance
      real(wp), intent(in), optional :: lsf_tol
      !> Step size tolerance
      real(wp), intent(in), optional :: step_tol
      !> Maximum search radius
      real(wp), intent(in), optional :: max_radius
      !> Debug flag
      logical, intent(in), optional :: debug

      type(moist_math_solver_onion_type), allocatable :: tmp

      allocate (tmp)

      ! Store pointers to context
      tmp%lsf => lsf
      tmp%anchor = anchor

      ! Set defaults
      tmp%initial_leb_num = 110
      tmp%max_leb_num = 302
      tmp%initial_step = 0.5_wp
      tmp%min_step = 0.001_wp
      tmp%max_steps = 100
      tmp%lsf_tol = 1.0e-6_wp
      tmp%step_tol = 1.0e-6_wp
      tmp%max_radius = 20.0_wp
      tmp%debug = .false.

      ! Override with user values
      if (present(initial_leb_num)) tmp%initial_leb_num = initial_leb_num
      if (present(max_leb_num)) tmp%max_leb_num = max_leb_num
      if (present(initial_step)) tmp%initial_step = initial_step
      if (present(min_step)) tmp%min_step = min_step
      if (present(max_steps)) tmp%max_steps = max_steps
      if (present(lsf_tol)) tmp%lsf_tol = lsf_tol
      if (present(max_radius)) tmp%max_radius = max_radius
      if (present(debug)) tmp%debug = debug

      ! Move to polymorphic output
      call move_alloc(tmp, solver)

   end subroutine new_onion_solver

   !> Solve the global search problem
   !>
   !> Performs hierarchical bisection with radius-dependent Lebedev refinement
   !> by scanning forward from anchor and stopping at first positive LSF.
   !>
   !> @param[inout] self   Solver instance
   !> @param[inout] x      On input: ignored. On output: best point found [3]
   !> @param[out]   error  Error handling
   subroutine onion_solve(self, x, error)
      class(moist_math_solver_onion_type), intent(inout), target :: self
      real(wp), dimension(:), intent(inout) :: x
      type(error_type), allocatable, intent(out) :: error

      integer :: oleb, num_leb, i_ang, i_shell, ierr
      real(wp), allocatable :: ang_grid(:, :), ang_weight(:)
      real(wp) :: r_curr, r_lower, r_upper, step_size
      real(wp) :: best_radius, best_lsf, shell_max_lsf, shell_max_point(3)
      integer :: best_ang, total_evals, level, iter, step
      integer, parameter :: max_refinements = 10
      logical :: found_bracket, found_positive

      if (.not. associated(self%lsf)) then
         call fatal_error(error, "Onion solver: LSF primitive not associated")
         return
      end if
      if (size(x) /= 3) then
         call fatal_error(error, "Onion solver: x must have size 3")
         return
      end if

      if (self%debug) then
         write (output_unit, '(x,a)') &
            '==================================  Starting Onion Solver  =================================='
      end if

      ! Initialize
      best_radius = huge(1.0_wp)
      best_lsf = huge(1.0_wp)
      best_ang = 1
      x = self%anchor
      total_evals = 0
      r_lower = 0.0_wp
      r_upper = self%max_radius
      step_size = self%initial_step
      found_positive = .false.

      num_leb = self%initial_leb_num
      call lebedev_order_from_num(num_leb, oleb, error)
      if (allocated(error)) return

      allocate (ang_grid(3, num_leb), ang_weight(num_leb))
      call get_angular_grid(oleb, ang_grid, ang_weight, error)
      if (allocated(error)) return

      !> First, do a scan until we find a positive LSF to establish initial bracket
      do step = 0, self%max_steps
         r_curr = step*step_size
         r_upper = r_curr
         call angular_scan_at_radius(self, self%anchor, r_curr, ang_grid, &
                                     shell_max_lsf, shell_max_point, error)
         if (shell_max_lsf > 0.0_wp) exit
         r_lower = r_curr

         if (self%debug) call onion_debug_callback(step + 1, i_ang, num_leb, r_curr, step_size, &
                                                   shell_max_lsf, r_lower, r_upper, found_positive)

      end do
      do iter = 1, self%max_steps

         r_curr = (r_upper + r_lower)/2.0_wp

         call angular_scan_at_radius(self, self%anchor, r_curr, ang_grid, &
                                     shell_max_lsf, shell_max_point, error)

         if (self%debug) call onion_debug_callback(iter, i_ang, num_leb, r_curr, step_size, &
                                                   shell_max_lsf, r_lower, r_upper, found_positive)

         ! Check convergence
         if (((step_size < self%step_tol) .or. (abs(shell_max_lsf) < self%lsf_tol)) &
             .and. (shell_max_lsf > 0.0_wp)) then
            x = shell_max_point
            exit
         end if

         if (shell_max_lsf < 0.0_wp) then
            r_lower = r_curr
            found_positive = .false.
         else if (shell_max_lsf > 0.0_wp) then
            r_upper = r_curr
            found_positive = .true.
         end if

         step_size = r_upper - r_lower

      end do

      if (self%debug) then
         write (output_unit, '(x,a)') 'Onion solver finished!'
      end if

   end subroutine onion_solve

   !> Perform single angular scan at given radius from anchor
   !>
   !> Evaluates LSF at all angular directions on a spherical shell and returns
   !> the maximum LSF value and the point where it occurs.
   !>
   !> @param[in]     self         Solver instance
   !> @param[in]     anchor       Center point for radial scan [3]
   !> @param[in]     radius       Radial distance from anchor
   !> @param[in]     ang_grid     Angular directions (unit vectors) [3, n_ang]
   !> @param[out]    max_lsf     Maximum LSF value found
   !> @param[out]    max_point    Point where maximum LSF occurs [3]
   !> @param[out]    max_ang_idx  Angular grid index of maximum
   !> @param[out]    n_evals      Number of LSF evaluations performed
   !> @param[out]    error        Error handling
   subroutine angular_scan_at_radius(self, anchor, radius, ang_grid, max_lsf, &
                                     max_point, error)
      type(moist_math_solver_onion_type), intent(inout) :: self
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: radius
      real(wp), intent(in) :: ang_grid(:, :)
      real(wp), intent(out) :: max_lsf
      real(wp), intent(out) :: max_point(3)
      integer :: max_ang_idx
      integer :: n_evals
      type(error_type), allocatable, intent(out) :: error

      integer :: i_ang, n_ang, ierr
      real(wp) :: direction(3), current_pos(3), lsf_val

      n_ang = size(ang_grid, 2)
      max_lsf = -huge(1.0_wp)
      max_ang_idx = 1
      max_point = anchor
      n_evals = 0

      ! Loop over all angular directions
      do i_ang = 1, n_ang
         direction = ang_grid(:, i_ang)
         current_pos = anchor + radius*direction

         call evaluate_lsf_at_point(self, current_pos, lsf_val, ierr)
         n_evals = n_evals + 1

         ! Skip failed evaluations
         if (ierr /= 0) cycle

         ! Track maximum LSF
         if (lsf_val > max_lsf) then
            max_lsf = lsf_val
            max_point = current_pos
            max_ang_idx = i_ang
         end if
      end do

   end subroutine angular_scan_at_radius

   !> Evaluate LSF at a given point using screened routine
   !>
   !> @param[in]     self      Solver instance
   !> @param[in]     point     Point to evaluate [3]
   !> @param[out]    lsf_val  LSF function value
   !> @param[out]    ierr      Error flag (0=success, nonzero=failure)
   subroutine evaluate_lsf_at_point(self, point, lsf_val, ierr)
      type(moist_math_solver_onion_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      real(wp), intent(out) :: lsf_val
      integer, intent(out) :: ierr

      real(wp) :: lsf0

      ierr = 0
      lsf_val = 0.0_wp

      ! Refresh per-point screening, then evaluate value-only.
      call self%lsf%prepare(point)
      call self%lsf%f012_r_screened(lsf0=lsf0)

      lsf_val = lsf0

      ! Check for NaN or Inf
      if (.not. (lsf_val == lsf_val)) then  ! NaN check
         ierr = 1
         return
      end if
      if (abs(lsf_val) > huge(1.0_wp)*0.1_wp) then  ! Inf check
         ierr = 2
         return
      end if

   end subroutine evaluate_lsf_at_point

   !> Debug callback for printing iteration progress
   !>
   !> @param[in] iter           Current step number
   !> @param[in] i_ang          Current angular direction index
   !> @param[in] num_leb        Total number of Lebedev directions
   !> @param[in] radius         Current radius from anchor
   !> @param[in] step_size      Current radial step size
   !> @param[in] shell_max_lsf Maximum LSF on current shell
   !> @param[in] r_lower        Lower bracket radius (all negative)
   !> @param[in] r_upper        Upper bracket radius (contains positive)
   !> @param[in] found_positive Whether positive LSF was found
   subroutine onion_debug_callback(iter, i_ang, num_leb, radius, step_size, &
                                   shell_max_lsf, r_lower, r_upper, found_positive)
      integer, intent(in) :: iter
      integer, intent(in) :: i_ang
      integer, intent(in) :: num_leb
      real(wp), intent(in) :: radius
      real(wp), intent(in) :: step_size
      real(wp), intent(in) :: shell_max_lsf
      real(wp), intent(in) :: r_lower
      real(wp), intent(in) :: r_upper
      logical, intent(in) :: found_positive

      character(len=1) :: status_flag

      ! Flag for positive found
      if (shell_max_lsf > 0.0_wp) then
         status_flag = 'U'
      else if (shell_max_lsf < 0.0_wp) then
         status_flag = 'L'
      end if

      ! Print header on first iteration
      if (iter == 1) then
         write (output_unit, '(x,a6,1x,a7,1x,a14,1x,a14,1x,a14,1x,a14,1x,a14,1x,a3)') &
            'Step', 'N_leb', 'R', 'maxLSF', 'R_lower', 'R_upper', 'Bracket', 'St'
         write (output_unit, '(x,a6,1x,a7,1x,a14,1x,a14,1x,a14,1x,a14,1x,a14,1x,a3)') &
            '------', '-------', '--------------', '--------------', '--------------', &
            '--------------', '--------------', '---'
      end if

      ! Print iteration data
      write (output_unit, '(x,i6,1x,i7,1x,e14.4,1x,e14.4,1x,f14.8,1x,f14.8,1x,'// &
             'es14.4,1x,a3)') &
         iter, num_leb, radius, shell_max_lsf, r_lower, r_upper, r_upper - r_lower, &
         status_flag

   end subroutine onion_debug_callback

   !> Clean up solver resources
   !>
   !> @param[inout] self  Solver instance
   subroutine onion_destroy(self)
      class(moist_math_solver_onion_type), intent(inout), target :: self

      ! Nullify pointer (do not deallocate - we don't own this)
      nullify (self%lsf)

   end subroutine onion_destroy

end module moist_math_solver_onion
