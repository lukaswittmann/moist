!> Bounded flexible right GMRES for approximate inverses inside a preconditioner
!>
!> - Solves A x = b approximately, starting from x0 = b (the identity action), with at
!>   most `m` Arnoldi steps and no restart
!> - Flexible: each direction z_j = M^{-1} v_j is kept, so M may change from step to step
!>   (it may itself run an inner iteration); x = b + Z y
!> - Classical Gram-Schmidt with a second pass, Givens rotations
!> - Stops when the residual estimate reaches `rtol` times |b - A b|, or at breakdown;
!>   a singular projected system drops its last directions
!> - Accepted only when x is finite and the recomputed residual |b - A x| is not above
!>   |b - A b| (relative slack `slack`); otherwise x = b
!> - Reductions in fixed blocks, so results do not depend on the thread count
module moist_model_moz_3d_preconditioner_krylov
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use moist_model_moz_solver_jfnk_kernels, only: blocked_dot, project_basis, update_basis, generate_givens, &
      & apply_givens, upper_backsolve
   implicit none(type, external)
   private

   public :: bounded_operator_type, bounded_gmres_type, bounded_gmres_result_type

   !> Relative slack of the acceptance test against the start residual
   real(wp), parameter :: slack = 1.0e-10_wp

   !> Linear operator with an optional (possibly varying) right preconditioner
   type, abstract :: bounded_operator_type
   contains
      !> y = A x
      procedure(matvec_interface), deferred :: matvec
      !> z = M^{-1} q; the identity unless overridden
      procedure :: precondition => operator_identity
   end type bounded_operator_type

   abstract interface
      !> y = A x
      !>
      !> @param[inout] self  operator and its scratch
      !> @param[in]    x     vector
      !> @param[out]   y     product
      !> @param[out]   ok    false when the product could not be formed
      subroutine matvec_interface(self, x, y, ok)
         import :: bounded_operator_type, wp
         implicit none(type, external)
         !> Operator
         class(bounded_operator_type), intent(inout) :: self
         !> Vector
         real(wp), contiguous, intent(in) :: x(:)
         !> Product
         real(wp), contiguous, intent(out) :: y(:)
         !> Whether the product was formed
         logical, intent(out) :: ok
      end subroutine matvec_interface
   end interface

   !> Outcome of one bounded solve
   type :: bounded_gmres_result_type
      !> Whether the solution was kept; false returns x = b
      logical :: accepted = .false.
      !> Arnoldi steps taken, also when the solve is not accepted
      integer :: iterations = 0
      !> Products with A, including the start and the acceptance check
      integer :: nmatvec = 0
   end type bounded_gmres_result_type

   !> Work space of the bounded solve, allocated once by `setup`
   type :: bounded_gmres_type
      private
      !> Most Arnoldi steps
      integer :: m = 0
      !> Vector length
      integer :: n = 0
      !> OpenMP team size
      integer :: nthreads = 1
      !> Orthonormal basis (n, m+1) and preconditioned directions (n, m)
      real(wp), allocatable :: v(:, :), z(:, :)
      !> Product and residual scratch (n)
      real(wp), allocatable :: w(:), r(:)
      !> Rotated Hessenberg matrix (m+1, m), rotated right-hand side (m+1)
      real(wp), allocatable :: h(:, :), g(:)
      !> Rotations, solution coefficients and projections (m)
      real(wp), allocatable :: cs(:), sn(:), y(:), proj(:), proj2(:)
   contains
      !> Allocate the work space
      procedure :: setup => gmres_setup
      !> Approximate A^{-1} b
      procedure :: solve => gmres_solve
   end type bounded_gmres_type

contains

   !> Identity preconditioner: z = q
   !>
   !> @param[inout] self  operator
   !> @param[in]    q     vector
   !> @param[out]   z     preconditioned vector
   !> @param[out]   ok    always true
   subroutine operator_identity(self, q, z, ok)
      !> Operator
      class(bounded_operator_type), intent(inout) :: self
      !> Vector
      real(wp), contiguous, intent(in) :: q(:)
      !> Preconditioned vector
      real(wp), contiguous, intent(out) :: z(:)
      !> Whether z was formed
      logical, intent(out) :: ok
      z = q
      ok = .true.
   end subroutine operator_identity

   !> Allocate the work space for vectors of length n and at most m steps
   !>
   !> @param[inout] self      work space
   !> @param[in]    n         vector length
   !> @param[in]    m         most Arnoldi steps, at least 1
   !> @param[in]    nthreads  OpenMP team size
   subroutine gmres_setup(self, n, m, nthreads)
      !> Work space
      class(bounded_gmres_type), intent(inout) :: self
      !> Vector length
      integer, intent(in) :: n
      !> Most Arnoldi steps
      integer, intent(in) :: m
      !> OpenMP team size
      integer, intent(in) :: nthreads

      self%n = n
      self%m = max(1, m)
      self%nthreads = max(1, nthreads)
      if (allocated(self%v)) deallocate (self%v)
      if (allocated(self%z)) deallocate (self%z)
      if (allocated(self%w)) deallocate (self%w)
      if (allocated(self%r)) deallocate (self%r)
      if (allocated(self%h)) deallocate (self%h)
      if (allocated(self%g)) deallocate (self%g)
      if (allocated(self%cs)) deallocate (self%cs)
      if (allocated(self%sn)) deallocate (self%sn)
      if (allocated(self%y)) deallocate (self%y)
      if (allocated(self%proj)) deallocate (self%proj)
      if (allocated(self%proj2)) deallocate (self%proj2)
      allocate (self%v(n, self%m + 1), self%z(n, self%m), self%w(n), self%r(n))
      allocate (self%h(self%m + 1, self%m), self%g(self%m + 1))
      allocate (self%cs(self%m), self%sn(self%m), self%y(self%m), self%proj(self%m), self%proj2(self%m))
   end subroutine gmres_setup

   !> Approximate A^{-1} b from the start b; x = b when not accepted
   !>
   !> @param[inout] self    work space after `setup`
   !> @param[inout] op      operator and preconditioner
   !> @param[in]    b       right-hand side (n)
   !> @param[out]   x       approximate solution (n)
   !> @param[in]    rtol    residual reduction that ends the iteration
   !> @param[out]   result  acceptance, steps and products
   !> @param[in]    nsteps  most steps of this solve, at most `m`; default `m`
   subroutine gmres_solve(self, op, b, x, rtol, result, nsteps)
      !> Work space
      class(bounded_gmres_type), intent(inout) :: self
      !> Operator
      class(bounded_operator_type), intent(inout) :: op
      !> Right-hand side
      real(wp), contiguous, intent(in) :: b(:)
      !> Approximate solution
      real(wp), contiguous, intent(out) :: x(:)
      !> Residual reduction that ends the iteration
      real(wp), intent(in) :: rtol
      !> Outcome
      type(bounded_gmres_result_type), intent(out) :: result
      !> Most steps of this solve
      integer, intent(in), optional :: nsteps
      real(wp) :: beta, hnorm, rnew
      integer :: j, k, i, mk, nt
      logical :: ok, feasible

      nt = self%nthreads
      mk = self%m
      if (present(nsteps)) mk = max(1, min(self%m, nsteps))
      x = b
      call op%matvec(x, self%w, ok)
      result%nmatvec = 1
      if (.not. ok) return
      self%r = b - self%w
      beta = sqrt(blocked_dot(self%r, self%r, nt))
      if (.not. ieee_is_finite(beta)) return
      if (beta == 0.0_wp) then
         result%accepted = .true.
         return
      end if
      self%v(:, 1) = self%r/beta
      self%g = 0.0_wp
      self%g(1) = beta
      ! Clear the entries of the previous solve
      self%h = 0.0_wp
      k = 0
      do j = 1, mk
         result%iterations = j
         call op%precondition(self%v(:, j), self%z(:, j), ok)
         if (.not. ok) self%z(:, j) = self%v(:, j)
         if (.not. all(ieee_is_finite(self%z(:, j)))) self%z(:, j) = self%v(:, j)
         call op%matvec(self%z(:, j), self%w, ok)
         result%nmatvec = result%nmatvec + 1
         if (.not. ok) exit
         ! Classical Gram-Schmidt, second pass against loss of orthogonality
         call project_basis(self%n, j, self%v(:, 1:j), self%w, self%proj(1:j), nt)
         call update_basis(self%n, j, self%v(:, 1:j), self%proj(1:j), -1.0_wp, self%w, nt)
         call project_basis(self%n, j, self%v(:, 1:j), self%w, self%proj2(1:j), nt)
         call update_basis(self%n, j, self%v(:, 1:j), self%proj2(1:j), -1.0_wp, self%w, nt)
         self%h(1:j, j) = self%proj(1:j) + self%proj2(1:j)
         hnorm = sqrt(blocked_dot(self%w, self%w, nt))
         if (.not. ieee_is_finite(hnorm) .or. .not. all(ieee_is_finite(self%h(1:j, j)))) exit
         self%h(j + 1, j) = hnorm
         do i = 1, j - 1
            call apply_givens(self%cs(i), self%sn(i), self%h(i, j), self%h(i + 1, j))
         end do
         call generate_givens(self%h(j, j), self%h(j + 1, j), self%cs(j), self%sn(j))
         call apply_givens(self%cs(j), self%sn(j), self%h(j, j), self%h(j + 1, j))
         call apply_givens(self%cs(j), self%sn(j), self%g(j), self%g(j + 1))
         k = j
         if (abs(self%g(j + 1)) <= rtol*beta) exit
         ! Breakdown: the next basis vector is undefined
         if (hnorm <= 10.0_wp*epsilon(1.0_wp)*maxval(abs(self%h(1:j, j)))) exit
         if (j < mk) self%v(:, j + 1) = self%w/hnorm
      end do
      ! A projected system singular up to rounding drops its last directions
      do while (k > 0)
         call upper_backsolve(self%h(1:k, 1:k), self%g(1:k), self%y(1:k), feasible, &
            & rank_tol=10.0_wp*epsilon(1.0_wp))
         if (feasible) exit
         k = k - 1
      end do
      if (k == 0) return
      call update_basis(self%n, k, self%z(:, 1:k), self%y(1:k), 1.0_wp, x, nt)
      if (.not. all(ieee_is_finite(x))) then
         x = b
         return
      end if
      call op%matvec(x, self%w, ok)
      result%nmatvec = result%nmatvec + 1
      if (ok) then
         self%r = b - self%w
         rnew = sqrt(blocked_dot(self%r, self%r, nt))
         ok = ieee_is_finite(rnew) .and. rnew <= (1.0_wp + slack)*beta
      end if
      if (.not. ok) then
         x = b
         return
      end if
      result%accepted = .true.
   end subroutine gmres_solve

end module moist_model_moz_3d_preconditioner_krylov
