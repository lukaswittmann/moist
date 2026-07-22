!> PCM-specific solver wrappers
!> This module provides solver routines for the PCM linear system A*q = rhs.
!> It wraps general moist linear algebra routines (LAPACK) and can be extended
!> with iterative solvers in the future.
module moist_model_component_pcm_solvers
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use moist_math_lapack, only: getrf, getrs, getri, potrf, potrs
   use moist_math_blas, only: dot, gemv
   implicit none
   private

   public :: solve_pcm_lu
   public :: solve_pcm_cholesky
   public :: solve_pcm_iterative
   public :: solve_pcm_inversion

contains

   !> Debug routine to check matrix properties
   !> Validates symmetry, positive definiteness, and detects NaN/Inf values.
   !> Not called during normal operation - use for debugging only.
   subroutine check_amat(amat, error)
      !> System matrix to validate
      real(wp), intent(in) :: amat(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n, i, j
      real(wp) :: sym_error, diag_min, diag_max, off_diag_max

      n = size(amat, 1)

      ! Check for square matrix
      if (size(amat, 2) /= n) then
         call fatal_error(error, "[check_amat] Matrix is not square")
         return
      end if

      ! Check for NaN or Inf values
      do i = 1, n
         do j = 1, n
            if (isnan(amat(i, j))) then
               call fatal_error(error, "[check_amat] Matrix contains NaN values")
               return
            end if
            if (.not. (abs(amat(i, j)) < huge(1.0_wp))) then
               call fatal_error(error, "[check_amat] Matrix contains Inf values")
               return
            end if
         end do
      end do

      ! Check symmetry
      sym_error = 0.0_wp
      do i = 1, n
         do j = i + 1, n
            sym_error = max(sym_error, abs(amat(i, j) - amat(j, i)))
         end do
      end do

      if (sym_error > 1.0e-10_wp) then
         write (*, '(A,ES12.4)') &
            "[check_amat] WARNING: Matrix not symmetric, max error = ", sym_error
      end if

      ! Check diagonal properties
      diag_min = minval([(amat(i, i), i=1, n)])
      diag_max = maxval([(amat(i, i), i=1, n)])
      off_diag_max = 0.0_wp
      do i = 1, n
         do j = 1, n
            if (i /= j) off_diag_max = max(off_diag_max, abs(amat(i, j)))
         end do
      end do

      if (diag_min <= 0.0_wp) then
         write (*, '(A,ES12.4)') "[check_amat] WARNING: Non-positive diagonal, min = ", diag_min
      end if

      ! Check for potential diagonal dominance issues
      do i = 1, n
         if (abs(amat(i, i)) < sum(abs(amat(i, :))) - abs(amat(i, i))) then
            write (*, '(A,I0)') "[check_amat] WARNING: Row ", i, " not diagonally dominant"
            exit
         end if
      end do

      ! Print summary
      write (*, '(A)') "[check_amat] Matrix validation summary:"
      write (*, '(A,I0)') "  Size: ", n
      write (*, '(A,ES12.4)') "  Symmetry error: ", sym_error
      write (*, '(A,ES12.4)') "  Min diagonal: ", diag_min
      write (*, '(A,ES12.4)') "  Max diagonal: ", diag_max
      write (*, '(A,ES12.4)') "  Max off-diagonal: ", off_diag_max

   end subroutine check_amat

   !> Solve PCM system using LU factorization
   !> Solves A*q = rhs via LAPACK's LU decomposition (DGETRF + DGETRS).
   subroutine solve_pcm_lu(amat, rhs, q, error)
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side vector (ngrid)
      real(wp), intent(in) :: rhs(:)
      !> Solution vector - surface charges (ngrid)
      real(wp), intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n, info
      integer, allocatable :: ipiv(:)
      real(wp), allocatable :: amat_copy(:, :), q_mat(:, :)

      n = size(amat, 1)
      ! Copy matrix (LAPACK overwrites input)
      allocate (amat_copy(n, n))
      amat_copy = amat

      ! Reshape RHS into 2D matrix for getrs (n, 1)
      allocate (q_mat(n, 1))
      q_mat(:, 1) = rhs

      ! Allocate pivot indices
      allocate (ipiv(n))

      ! LU factorization
      call getrf(amat_copy, ipiv, info)
      if (info /= 0) then
         write (*, '(A,I0)') "[solve_pcm_lu] LAPACK getrf failed with info = ", info
         write (*, '(A,I0)') "[solve_pcm_lu] Matrix size n = ", n
         call fatal_error(error, "[solve_pcm_lu] LAPACK getrf failed")
         return
      end if

      ! Solve using factorization (getrs expects 2D matrix)
      call getrs(amat_copy, q_mat, ipiv, info)
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_lu] LAPACK getrs failed")
         return
      end if

      ! Extract solution from 2D matrix
      q = q_mat(:, 1)

   end subroutine solve_pcm_lu

   !> Solve PCM system using Cholesky factorization
   !> Solves A*q = rhs via LAPACK's Cholesky decomposition (DPOTRF + DPOTRS).
   !> Assumes A is symmetric positive definite - faster than LU for such matrices.
   subroutine solve_pcm_cholesky(amat, rhs, q, error)
      !> System matrix (ngrid, ngrid) - must be symmetric positive definite
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side vector (ngrid)
      real(wp), intent(in) :: rhs(:)
      !> Solution vector - surface charges (ngrid)
      real(wp), intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n, info
      real(wp), allocatable :: amat_copy(:, :), q_mat(:, :)

      n = size(amat, 1)

      ! Copy matrix (LAPACK overwrites input)
      allocate (amat_copy(n, n))
      amat_copy = amat

      ! Reshape RHS into 2D matrix for potrs (n, 1)
      allocate (q_mat(n, 1))
      q_mat(:, 1) = rhs

      ! Cholesky factorization: A = L*L^T (lower triangular)
      call potrf(amat_copy, info, uplo='l')
      if (info /= 0) then
         if (info > 0) then
            call fatal_error(error, "[solve_pcm_cholesky] Matrix not positive definite")
         else
            call fatal_error(error, "[solve_pcm_cholesky] LAPACK potrf failed")
         end if
         return
      end if

      ! Solve using factorization (potrs expects 2D matrix)
      call potrs(amat_copy, q_mat, info, uplo='l')
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_cholesky] LAPACK potrs failed")
         return
      end if

      ! Extract solution from 2D matrix
      q = q_mat(:, 1)

   end subroutine solve_pcm_cholesky

   !> Solve PCM system using matrix inversion
   !> Computes A^(-1) and then q = A^(-1)*rhs. Efficient if matrix is reused many times.
   subroutine solve_pcm_inversion(amat, rhs, q, error)
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side vector (ngrid)
      real(wp), intent(in) :: rhs(:)
      !> Solution vector - surface charges (ngrid)
      real(wp), intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n, info
      integer, allocatable :: ipiv(:)
      real(wp), allocatable :: amat_inv(:, :)

      n = size(amat, 1)

      ! Copy matrix for inversion
      allocate (amat_inv(n, n))
      amat_inv = amat

      ! Allocate pivot indices
      allocate (ipiv(n))

      ! LU factorization
      call getrf(amat_inv, ipiv, info)
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_inversion] LAPACK getrf failed")
         return
      end if

      ! Compute inverse
      call getri(amat_inv, ipiv, info=info)
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_inversion] LAPACK getri failed")
         return
      end if

      ! Multiply: q = A^(-1) * rhs
      call gemv(amat_inv, rhs, q)

   end subroutine solve_pcm_inversion

   !> Solve PCM system via preconditioned Conjugate Gradient
   subroutine solve_pcm_iterative(amat, rhs, q, tol, maxiter, error)
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side vector (ngrid)
      real(wp), intent(in) :: rhs(:)
      !> Solution vector - surface charges (ngrid)
      real(wp), intent(out) :: q(:)
      !> Convergence tolerance
      real(wp), intent(in) :: tol
      !> Maximum iterations
      integer, intent(in) :: maxiter
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n, iter, restart, i
      real(wp) :: alpha, beta, rho_old, rho_new, res_norm, pAp, rhs_norm
      real(wp) :: tol_eff
      logical :: converged
      real(wp), allocatable :: r(:), p(:), Ap(:), z(:), diag_inv(:)
      real(wp), parameter :: diag_tol = 100.0_wp*epsilon(1.0_wp)

      n = size(amat, 1)

      ! Allocate work arrays
      allocate (r(n), p(n), Ap(n), z(n), diag_inv(n))

      ! Jacobi preconditioner (diagonal)
      do i = 1, n
         if (abs(amat(i, i)) <= diag_tol) then
            call fatal_error(error, "[CG] Jacobi preconditioner failed: near-zero diagonal entry in A")
            return
         end if
         diag_inv(i) = 1.0_wp/amat(i, i)
      end do

      rhs_norm = sqrt(dot(rhs, rhs))
      tol_eff = tol*max(1.0_wp, rhs_norm)

      ! Initial guess q = 0
      q = 0.0_wp

      do restart = 0, 1

         ! True residual: r = b - A*q (first pass: r = b since q = 0)
         if (restart == 0) then
            r = rhs
         else
            call gemv(amat, q, Ap)
            r = rhs - Ap
         end if
         res_norm = sqrt(dot(r, r))
         if (res_norm < tol_eff) return

         ! Preconditioned direction: z = M^(-1) r, p = z, rho = (r, z)
         z = r*diag_inv
         p = z
         rho_old = dot(r, z)

         converged = .false.
         do iter = 1, maxiter

            call gemv(amat, p, Ap)

            ! Step size: alpha = (r, z) / (p, A*p).
            ! For SPD A and p /= 0, pAp > 0; a non-positive value means the
            ! matrix is not positive definite (or rounding destroyed it).
            pAp = dot(p, Ap)
            if (pAp <= 0.0_wp) then
               call fatal_error(error, "[CG] Matrix is not positive definite")
               return
            end if
            alpha = rho_old/pAp

            q = q + alpha*p
            r = r - alpha*Ap

            z = r*diag_inv
            rho_new = dot(r, z)
            res_norm = sqrt(dot(r, r))

            if (isnan(res_norm)) then
               call fatal_error(error, "[CG] NaN detected in residual")
               return
            end if

            if (res_norm < tol_eff) then
               converged = .true.
               exit
            end if

            beta = rho_new/rho_old
            p = z + beta*p
            rho_old = rho_new

         end do

         if (.not. converged) then
            call fatal_error(error, "[CG] Failed to converge within maximum iterations")
            return
         end if

         ! Verify convergence against the residual b - A*q
         call gemv(amat, q, Ap)
         r = rhs - Ap
         res_norm = sqrt(dot(r, r))
         if (res_norm < tol_eff) return

      end do

      call fatal_error(error, &
         & "[CG] Recurrence converged but true residual remains above tolerance")
      return

   end subroutine solve_pcm_iterative

end module moist_model_component_pcm_solvers
