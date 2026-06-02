!> Generic deflation operator for multi-root enumeration (Farrell, Birkisson
!> & Funke 2015, arXiv:1410.5620).
!>
!> The operator encodes a list of already-known roots x*_i and produces the
!> scalar multiplier
!>
!>    M(x) = prod_i ( ||x - x*_i||^{-p} + alpha )
!>
!> which tends to infinity as x -> x*_i and smoothly to alpha far from any
!> known root. Multiplying a target residual (or objective) by M redirects a
!> local solver away from the already-found roots so iterating the same local
!> solver from the same seed can enumerate new roots.
!>
!> This module is deliberately solver-agnostic: it is consumed by both the
!> SLSQP-deflation and Newton-deflation solver wrappers, which differ only in
!> *which* quantity they multiply by M (objective scalar vs residual vector).
module moist_math_solver_deflation
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   implicit none
   private

   public :: moist_deflation_operator_type

   !> Operator carrying the accumulated known roots and the scalar M(x), grad_M(x).
   type :: moist_deflation_operator_type
      !> Dimension of the search space (length of x).
      integer :: n_dim = 0
      !> Maximum number of roots that can be stored.
      integer :: max_roots = 0
      !> Current number of known roots.
      integer :: n_known = 0
      !> Exponent p in ||x - x*_i||^{-p}. Farrell default: 2.
      integer :: p_power = 2
      !> Additive shift alpha; prevents M from vanishing far from known roots.
      real(wp) :: alpha_shift = 1.0_wp
      !> L2 tolerance for considering two roots to be the same point.
      real(wp) :: dedup_tol = 1.0e-6_wp
      !> Known-root list (n_dim, max_roots); first n_known columns are used.
      real(wp), allocatable :: known_roots(:, :)
   contains
      procedure :: init => deflation_init
      procedure :: reset => deflation_reset
      procedure :: append_root => deflation_append_root
      procedure :: multiplier => deflation_multiplier
      procedure :: gradient => deflation_gradient
      procedure :: is_near_known => deflation_is_near_known
      procedure :: destroy => deflation_destroy
   end type moist_deflation_operator_type

contains

   !> Initialise the operator with given dimensions and parameters.
   !>
   !> @param[inout] self        Operator instance
   !> @param[in]    n_dim       Length of the search-space vector x
   !> @param[in]    max_roots   Maximum number of roots to store
   !> @param[in]    p_power     Exponent p in ||x - x*||^{-p} (optional, default 2)
   !> @param[in]    alpha_shift Additive shift alpha (optional, default 1.0)
   !> @param[in]    dedup_tol   L2 tolerance for root-identity test (optional, default 1e-6)
   !> @param[out]   error       Error status
   subroutine deflation_init(self, n_dim, max_roots, error, p_power, alpha_shift, dedup_tol)
      !> Operator instance
      class(moist_deflation_operator_type), intent(inout) :: self
      !> Length of the search-space vector
      integer, intent(in) :: n_dim
      !> Maximum number of roots that can be stored
      integer, intent(in) :: max_roots
      !> Error status
      type(error_type), allocatable, intent(out) :: error
      !> Deflation exponent p
      integer, intent(in), optional :: p_power
      !> Additive shift alpha
      real(wp), intent(in), optional :: alpha_shift
      !> Root-identity L2 tolerance
      real(wp), intent(in), optional :: dedup_tol

      if (n_dim <= 0) then
         call fatal_error(error, "deflation: n_dim must be positive")
         return
      end if
      if (max_roots <= 0) then
         call fatal_error(error, "deflation: max_roots must be positive")
         return
      end if

      self%n_dim = n_dim
      self%max_roots = max_roots
      self%n_known = 0
      if (present(p_power)) self%p_power = p_power
      if (present(alpha_shift)) self%alpha_shift = alpha_shift
      if (present(dedup_tol)) self%dedup_tol = dedup_tol

      if (self%p_power <= 0) then
         call fatal_error(error, "deflation: p_power must be positive")
         return
      end if
      if (self%alpha_shift < 0.0_wp) then
         call fatal_error(error, "deflation: alpha_shift must be non-negative")
         return
      end if

      if (allocated(self%known_roots)) deallocate (self%known_roots)
      allocate (self%known_roots(n_dim, max_roots), source=0.0_wp)
   end subroutine deflation_init

   !> Forget all known roots (operator behaves as identity on M until next append).
   !>
   !> @param[inout] self Operator instance
   subroutine deflation_reset(self)
      !> Operator instance
      class(moist_deflation_operator_type), intent(inout) :: self
      self%n_known = 0
   end subroutine deflation_reset

   !> Append a newly-found root if it is not a duplicate of any known root.
   !>
   !> @param[inout] self     Operator instance
   !> @param[in]    x        New root candidate (length n_dim)
   !> @param[out]   accepted .true. if the point was added, .false. on duplicate or full buffer
   subroutine deflation_append_root(self, x, accepted)
      !> Operator instance
      class(moist_deflation_operator_type), intent(inout) :: self
      !> Candidate root
      real(wp), dimension(:), intent(in) :: x
      !> True on successful append, false if duplicate or buffer full
      logical, intent(out) :: accepted

      accepted = .false.
      if (size(x) /= self%n_dim) return
      if (self%n_known >= self%max_roots) return
      if (self%is_near_known(x)) return

      self%n_known = self%n_known + 1
      self%known_roots(:, self%n_known) = x
      accepted = .true.
   end subroutine deflation_append_root

   !> Evaluate M(x) = prod_i ( ||x - x*_i||^{-p} + alpha ).
   !>
   !> When n_known == 0 this returns 1.0 (no deflation yet), so the inner
   !> solver sees the original problem on the first invocation.
   !>
   !> @param[in] self Operator instance
   !> @param[in] x    Evaluation point
   !> @return    Scalar multiplier M(x)
   function deflation_multiplier(self, x) result(m)
      !> Operator instance
      class(moist_deflation_operator_type), intent(in) :: self
      !> Evaluation point
      real(wp), dimension(:), intent(in) :: x
      !> Result: M(x)
      real(wp) :: m

      integer :: i
      real(wp) :: r
      real(wp) :: r_floor

      m = 1.0_wp
      if (self%n_known == 0) return

      r_floor = max(self%dedup_tol, epsilon(1.0_wp))
      do i = 1, self%n_known
         r = max(norm2(x - self%known_roots(1:self%n_dim, i)), r_floor)
         m = m*(r**(-self%p_power) + self%alpha_shift)
      end do
   end function deflation_multiplier

   !> Evaluate grad_M(x) via the product rule.
   !>
   !> grad_M = sum_i [ -p * (x - x*_i) / r_i^{p+2} ] * prod_{j /= i} (r_j^{-p} + alpha)
   !>
   !> Implemented as a single pass that builds the full product M, then for
   !> each i scales the per-i gradient contribution by M / factor_i. This keeps
   !> the cost at O(n_known * n_dim) without the explicit double loop.
   !>
   !> @param[in]  self   Operator instance
   !> @param[in]  x      Evaluation point
   !> @param[out] grad_m Gradient of M at x (length n_dim)
   subroutine deflation_gradient(self, x, grad_m)
      !> Operator instance
      class(moist_deflation_operator_type), intent(in) :: self
      !> Evaluation point
      real(wp), dimension(:), intent(in) :: x
      !> Output gradient
      real(wp), dimension(:), intent(out) :: grad_m

      integer :: i
      real(wp) :: r, factor_i, m_total
      real(wp) :: r_floor
      real(wp), allocatable :: factor(:)

      grad_m = 0.0_wp
      if (self%n_known == 0) return

      r_floor = max(self%dedup_tol, epsilon(1.0_wp))
      allocate (factor(self%n_known))
      m_total = 1.0_wp
      do i = 1, self%n_known
         r = max(norm2(x - self%known_roots(1:self%n_dim, i)), r_floor)
         factor_i = r**(-self%p_power) + self%alpha_shift
         factor(i) = factor_i
         m_total = m_total*factor_i
      end do

      do i = 1, self%n_known
         r = max(norm2(x - self%known_roots(1:self%n_dim, i)), r_floor)
         ! Per-term local gradient:  -p * (x - x*_i) / r^{p+2}
         ! Multiply by the product of the *other* factors = m_total / factor(i).
         grad_m = grad_m + (m_total/factor(i))* &
                  (-real(self%p_power, wp))* &
                  (x - self%known_roots(1:self%n_dim, i))/r**(self%p_power + 2)
      end do

      deallocate (factor)
   end subroutine deflation_gradient

   !> Test whether x is within dedup_tol of any stored root.
   !>
   !> @param[in] self Operator instance
   !> @param[in] x    Query point
   !> @return    .true. if x is near (L2) an existing root, .false. otherwise
   function deflation_is_near_known(self, x) result(hit)
      !> Operator instance
      class(moist_deflation_operator_type), intent(in) :: self
      !> Query point
      real(wp), dimension(:), intent(in) :: x
      !> True if close to any stored root
      logical :: hit

      integer :: i

      hit = .false.
      do i = 1, self%n_known
         if (norm2(x - self%known_roots(1:self%n_dim, i)) < self%dedup_tol) then
            hit = .true.
            return
         end if
      end do
   end function deflation_is_near_known

   !> Release storage and reset state.
   !>
   !> @param[inout] self Operator instance
   subroutine deflation_destroy(self)
      !> Operator instance
      class(moist_deflation_operator_type), intent(inout) :: self
      if (allocated(self%known_roots)) deallocate (self%known_roots)
      self%n_dim = 0
      self%max_roots = 0
      self%n_known = 0
   end subroutine deflation_destroy

end module moist_math_solver_deflation
