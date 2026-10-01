!> Reference quadrature rules on [-1, 1], the first layer of a radial grid
!>
!> A rule generates nodes and ordinary `dx` weights for a requested count;
!> optional finite bounds map the canonical interval affinely onto [lower, upper]
!>
!> - Chebyshev second kind: nodes x_i = cos(i*pi/(n+1)), i = 1..n, in
!>   descending order; weights pi/(n+1)*sqrt((1-x_i)*(1+x_i)) represent plain
!>   dx and are evaluated at the stored node, so every later factor sees the
!>   same rounded x_i; exact for sqrt(1-x^2)*p(x) with deg(p) <= 2n-1
!> - Midpoint: nodes x_i = -1 + (2i - 1)/n, i = 1..n, in ascending order;
!>   equal weights 2/n; exact through degree 1
!> - Gauss-Legendre: nodes are the roots of P_n in ascending order, found by
!>   Newton iteration on the three-term recurrence; weights
!>   2/((1 - x^2) P_n'(x)^2); exact through degree 2n-1
module moist_math_grid_radial_rule
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   implicit none(type, external)
   private

   public :: moist_math_grid_radial_rule_type
   public :: check_rule_request
   public :: scale_rule_interval
   public :: moist_math_grid_radial_rule_chebyshev2_type
   public :: new_chebyshev2_rule
   public :: moist_math_grid_radial_rule_midpoint_type
   public :: new_midpoint_rule
   public :: moist_math_grid_radial_rule_gauss_legendre_type
   public :: new_gauss_legendre_rule

   !> Abstract reference quadrature rule; carries only rule-specific parameters
   type, abstract :: moist_math_grid_radial_rule_type
   contains
      !> Generate n nodes and `dx` weights on [-1, 1] or on [lower, upper]
      procedure(rule_generate_i), deferred :: generate
   end type moist_math_grid_radial_rule_type

   !> Deferred generator shared by every concrete rule
   abstract interface
      !> Generate n nodes and `dx` weights
      !>
      !> @param[in]  self   Rule instance
      !> @param[in]  n      Number of nodes (n >= 1)
      !> @param[out] x      Nodes, shape (n)
      !> @param[out] w      Weights representing dx, shape (n)
      !> @param[out] error  Set on an invalid count or interval
      !> @param[in]  lower  Lower interval bound, finite (default -1)
      !> @param[in]  upper  Upper interval bound, finite, > lower (default 1)
      subroutine rule_generate_i(self, n, x, w, error, lower, upper)
         import :: moist_math_grid_radial_rule_type, wp, error_type
         implicit none(type, external)
         !> Rule instance
         class(moist_math_grid_radial_rule_type), intent(in) :: self
         !> Number of nodes
         integer, intent(in) :: n
         !> Nodes
         real(wp), allocatable, intent(out) :: x(:)
         !> Weights representing dx
         real(wp), allocatable, intent(out) :: w(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Lower interval bound
         real(wp), intent(in), optional :: lower
         !> Upper interval bound
         real(wp), intent(in), optional :: upper
      end subroutine rule_generate_i
   end interface

   !> Chebyshev second-kind rule; no parameters
   type, extends(moist_math_grid_radial_rule_type) :: moist_math_grid_radial_rule_chebyshev2_type
   contains
      !> Generate nodes and dx weights
      procedure :: generate => chebyshev2_generate
   end type moist_math_grid_radial_rule_chebyshev2_type

   !> Composite midpoint rule; no parameters
   type, extends(moist_math_grid_radial_rule_type) :: moist_math_grid_radial_rule_midpoint_type
   contains
      !> Generate nodes and dx weights
      procedure :: generate => midpoint_generate
   end type moist_math_grid_radial_rule_midpoint_type

   !> Maximum Newton iterations per root
   integer, parameter :: max_newton = 100

   !> Gauss-Legendre rule; no parameters
   type, extends(moist_math_grid_radial_rule_type) :: moist_math_grid_radial_rule_gauss_legendre_type
   contains
      !> Generate nodes and dx weights
      procedure :: generate => gauss_legendre_generate
   end type moist_math_grid_radial_rule_gauss_legendre_type

contains

   !> Validate a rule request: node count and optional interval bounds
   !>
   !> @param[in]  n      Number of nodes (n >= 1)
   !> @param[in]  label  Error-message label
   !> @param[out] error  Set on an invalid count or interval
   !> @param[in]  lower  Lower interval bound, finite
   !> @param[in]  upper  Upper interval bound, finite, > lower
   subroutine check_rule_request(n, label, error, lower, upper)
      !> Number of nodes
      integer, intent(in) :: n
      !> Error-message label
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Lower interval bound
      real(wp), intent(in), optional :: lower
      !> Upper interval bound
      real(wp), intent(in), optional :: upper

      real(wp) :: lo, hi

      if (n < 1) then
         call fatal_error(error, label//": number of nodes must be >= 1")
         return
      end if
      lo = -1.0_wp
      hi = 1.0_wp
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      if (.not. (ieee_is_finite(lo) .and. ieee_is_finite(hi))) then
         call fatal_error(error, label//": interval bounds must be finite")
         return
      end if
      if (.not. lo < hi) then
         call fatal_error(error, label//": lower bound must be below upper bound")
         return
      end if
   end subroutine check_rule_request

   !> Map canonical nodes and weights affinely onto [lower, upper]
   !>
   !> x -> (lower + upper)/2 + (upper - lower)/2 * x, w -> (upper - lower)/2 * w;
   !> no-op when both bounds are absent, so the canonical bits are kept
   !>
   !> @param[in,out] x      Nodes on [-1, 1] on entry, on [lower, upper] on exit
   !> @param[in,out] w      Weights for [-1, 1] on entry, for [lower, upper] on exit
   !> @param[in]     lower  Lower interval bound (default -1)
   !> @param[in]     upper  Upper interval bound (default 1)
   pure subroutine scale_rule_interval(x, w, lower, upper)
      !> Nodes
      real(wp), intent(inout) :: x(:)
      !> Weights
      real(wp), intent(inout) :: w(:)
      !> Lower interval bound
      real(wp), intent(in), optional :: lower
      !> Upper interval bound
      real(wp), intent(in), optional :: upper

      real(wp) :: lo, hi, mid, half

      if (.not. (present(lower) .or. present(upper))) return
      lo = -1.0_wp
      hi = 1.0_wp
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      mid = 0.5_wp*(lo + hi)
      half = 0.5_wp*(hi - lo)
      x(:) = mid + half*x(:)
      w(:) = half*w(:)
   end subroutine scale_rule_interval

   !> Create a Chebyshev second-kind rule
   !>
   !> @param[out] rule  New rule
   pure subroutine new_chebyshev2_rule(rule)
      !> New rule
      type(moist_math_grid_radial_rule_chebyshev2_type), intent(out) :: rule
   end subroutine new_chebyshev2_rule

   !> Generate n Chebyshev second-kind nodes (descending) and dx weights
   !>
   !> @param[in]  self   Rule instance
   !> @param[in]  n      Number of nodes (n >= 1)
   !> @param[out] x      Nodes, shape (n), descending
   !> @param[out] w      Weights representing dx, shape (n)
   !> @param[out] error  Set on an invalid count or interval
   !> @param[in]  lower  Lower interval bound, finite (default -1)
   !> @param[in]  upper  Upper interval bound, finite, > lower (default 1)
   subroutine chebyshev2_generate(self, n, x, w, error, lower, upper)
      !> Rule instance
      class(moist_math_grid_radial_rule_chebyshev2_type), intent(in) :: self
      !> Number of nodes
      integer, intent(in) :: n
      !> Nodes
      real(wp), allocatable, intent(out) :: x(:)
      !> Weights representing dx
      real(wp), allocatable, intent(out) :: w(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Lower interval bound
      real(wp), intent(in), optional :: lower
      !> Upper interval bound
      real(wp), intent(in), optional :: upper

      integer :: i, stat

      call check_rule_request(n, "Chebyshev-II rule", error, lower, upper)
      if (allocated(error)) return

      allocate (x(n), stat=stat)
      if (stat == 0) allocate (w(n), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Chebyshev-II rule: failed to allocate nodes")
         return
      end if
      do i = 1, n
         x(i) = cos(real(i, wp)*pi/real(n + 1, wp))
         w(i) = pi/real(n + 1, wp)*sqrt((1.0_wp - x(i))*(1.0_wp + x(i)))
      end do

      call scale_rule_interval(x, w, lower, upper)
   end subroutine chebyshev2_generate

   !> Create a composite midpoint rule
   !>
   !> @param[out] rule  New rule
   pure subroutine new_midpoint_rule(rule)
      !> New rule
      type(moist_math_grid_radial_rule_midpoint_type), intent(out) :: rule
   end subroutine new_midpoint_rule

   !> Generate n midpoint nodes (ascending) and dx weights
   !>
   !> @param[in]  self   Rule instance
   !> @param[in]  n      Number of nodes (n >= 1)
   !> @param[out] x      Nodes, shape (n), ascending
   !> @param[out] w      Weights representing dx, shape (n)
   !> @param[out] error  Set on an invalid count or interval
   !> @param[in]  lower  Lower interval bound, finite (default -1)
   !> @param[in]  upper  Upper interval bound, finite, > lower (default 1)
   subroutine midpoint_generate(self, n, x, w, error, lower, upper)
      !> Rule instance
      class(moist_math_grid_radial_rule_midpoint_type), intent(in) :: self
      !> Number of nodes
      integer, intent(in) :: n
      !> Nodes
      real(wp), allocatable, intent(out) :: x(:)
      !> Weights representing dx
      real(wp), allocatable, intent(out) :: w(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Lower interval bound
      real(wp), intent(in), optional :: lower
      !> Upper interval bound
      real(wp), intent(in), optional :: upper

      integer :: i, stat
      real(wp) :: wx

      call check_rule_request(n, "Midpoint rule", error, lower, upper)
      if (allocated(error)) return

      allocate (x(n), stat=stat)
      if (stat == 0) allocate (w(n), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Midpoint rule: failed to allocate nodes")
         return
      end if
      wx = 2.0_wp/real(n, wp)
      do i = 1, n
         x(i) = -1.0_wp + (2.0_wp*real(i, wp) - 1.0_wp)/real(n, wp)
         w(i) = wx
      end do

      call scale_rule_interval(x, w, lower, upper)
   end subroutine midpoint_generate

   !> Create a Gauss-Legendre rule
   !>
   !> @param[out] rule  New rule
   pure subroutine new_gauss_legendre_rule(rule)
      !> New rule
      type(moist_math_grid_radial_rule_gauss_legendre_type), intent(out) :: rule
   end subroutine new_gauss_legendre_rule

   !> Generate n Gauss-Legendre nodes (ascending) and dx weights
   !>
   !> Roots are found for x >= 0 and mirrored, so nodes are exactly symmetric;
   !> the center node of an odd rule is exactly 0
   !>
   !> @param[in]  self   Rule instance
   !> @param[in]  n      Number of nodes (n >= 1)
   !> @param[out] x      Nodes, shape (n), ascending
   !> @param[out] w      Weights representing dx, shape (n)
   !> @param[out] error  Set on an invalid count or interval, or a Newton failure
   !> @param[in]  lower  Lower interval bound, finite (default -1)
   !> @param[in]  upper  Upper interval bound, finite, > lower (default 1)
   subroutine gauss_legendre_generate(self, n, x, w, error, lower, upper)
      !> Rule instance
      class(moist_math_grid_radial_rule_gauss_legendre_type), intent(in) :: self
      !> Number of nodes
      integer, intent(in) :: n
      !> Nodes
      real(wp), allocatable, intent(out) :: x(:)
      !> Weights representing dx
      real(wp), allocatable, intent(out) :: w(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Lower interval bound
      real(wp), intent(in), optional :: lower
      !> Upper interval bound
      real(wp), intent(in), optional :: upper

      integer :: i, iter, stat
      logical :: converged
      real(wp) :: z, dz, p, dp, tol

      call check_rule_request(n, "Gauss-Legendre rule", error, lower, upper)
      if (allocated(error)) return

      allocate (x(n), stat=stat)
      if (stat == 0) allocate (w(n), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Gauss-Legendre rule: failed to allocate nodes")
         return
      end if

      ! Newton converges quadratically, so one step below this bound leaves
      ! an error far under round-off
      tol = 100.0_wp*epsilon(1.0_wp)
      do i = 1, n/2
         ! Tricomi initial guess for the i-th largest root
         z = cos(pi*(real(i, wp) - 0.25_wp)/(real(n, wp) + 0.5_wp))
         converged = .false.
         do iter = 1, max_newton
            call legendre_eval(n, z, p, dp)
            dz = p/dp
            z = z - dz
            if (abs(dz) <= tol) then
               converged = .true.
               exit
            end if
         end do
         if (.not. converged) then
            call fatal_error(error, "Gauss-Legendre rule: Newton iteration did not converge")
            return
         end if
         call legendre_eval(n, z, p, dp)
         x(n + 1 - i) = z
         x(i) = -z
         w(n + 1 - i) = 2.0_wp/((1.0_wp - z)*(1.0_wp + z)*dp*dp)
         w(i) = w(n + 1 - i)
      end do
      if (mod(n, 2) == 1) then
         z = 0.0_wp
         call legendre_eval(n, z, p, dp)
         x(n/2 + 1) = z
         w(n/2 + 1) = 2.0_wp/(dp*dp)
      end if

      call scale_rule_interval(x, w, lower, upper)
   end subroutine gauss_legendre_generate

   !> Evaluate P_n and P_n' at z in (-1, 1) by the three-term recurrence
   !>
   !> @param[in]  n   Polynomial degree (n >= 1)
   !> @param[in]  z   Evaluation point, |z| < 1
   !> @param[out] p   P_n(z)
   !> @param[out] dp  P_n'(z)
   pure subroutine legendre_eval(n, z, p, dp)
      !> Polynomial degree
      integer, intent(in) :: n
      !> Evaluation point
      real(wp), intent(in) :: z
      !> P_n(z)
      real(wp), intent(out) :: p
      !> P_n'(z)
      real(wp), intent(out) :: dp

      integer :: j
      real(wp) :: pm1, pm2

      pm1 = 1.0_wp
      p = z
      do j = 2, n
         pm2 = pm1
         pm1 = p
         p = (real(2*j - 1, wp)*z*pm1 - real(j - 1, wp)*pm2)/real(j, wp)
      end do
      ! P_n' = n*(P_{n-1} - z*P_n)/(1 - z^2)
      dp = real(n, wp)*(pm1 - z*p)/((1.0_wp - z)*(1.0_wp + z))
   end subroutine legendre_eval

end module moist_math_grid_radial_rule
