!> PCM-specific solver wrappers
!> This module provides solver routines for the PCM linear system A*q = rhs.
!> It wraps general moist linear algebra routines (LAPACK) and can be extended
!> with iterative solvers in the future.
module moist_model_component_pcm_solvers
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use moist_math_lapack, only: getrf, getrs, getri, potrf, potrs
   use moist_math_blas, only: dot, gemv, gemm
   implicit none (type, external)
   private

   public :: solve_pcm_lu
   public :: solve_pcm_cholesky
   public :: solve_pcm_iterative
   public :: solve_pcm_inversion
   public :: factorize_pcm_lu, solve_pcm_lu_factored
   public :: factorize_pcm_cholesky, solve_pcm_cholesky_factored
   public :: invert_pcm_matrix, apply_pcm_inverse

contains

   !> Solve PCM system using LU factorization
   !> Solves A*q = rhs via LAPACK's LU decomposition (DGETRF + DGETRS).
   subroutine solve_pcm_lu(amat, rhs, q, error, unit)
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side vector (ngrid)
      real(wp), intent(in) :: rhs(:)
      !> Solution vector - surface charges (ngrid)
      real(wp), intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Output unit for the diagnostics; defaults to standard output
      integer, intent(in), optional :: unit

      integer :: n, info
      integer, allocatable :: ipiv(:)
      real(wp), allocatable :: amat_copy(:, :), q_mat(:, :)
      integer :: iunit

      iunit = output_unit
      if (present(unit)) iunit = unit

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
         write (iunit, "(A,I0)") "[solve_pcm_lu] LAPACK getrf failed with info = ", info
         write (iunit, "(A,I0)") "[solve_pcm_lu] Matrix size n = ", n
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
      call potrf(amat_copy, info, uplo="l")
      if (info /= 0) then
         if (info > 0) then
            call fatal_error(error, "[solve_pcm_cholesky] Matrix not positive definite")
         else
            call fatal_error(error, "[solve_pcm_cholesky] LAPACK potrf failed")
         end if
         return
      end if

      ! Solve using factorization (potrs expects 2D matrix)
      call potrs(amat_copy, q_mat, info, uplo="l")
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

   !> LU-factorize the PCM matrix once, for repeated multi-column solves
   !>
   !> The factor and the pivots are what [[solve_pcm_lu_factored]] consumes.
   !> This is the same DGETRF the single-column [[solve_pcm_lu]] runs; it is
   !> kept apart so that the charge solve of the energy path is never routed
   !> through a cache.
   !>
   !> @param[in]  amat   System matrix (ngrid, ngrid)
   !> @param[out] factor LU factor (ngrid, ngrid)
   !> @param[out] ipiv   Pivot indices (ngrid)
   !> @param[out] error  Error handling
   !> @param[in]  unit   Output unit for the diagnostics; defaults to standard output
   subroutine factorize_pcm_lu(amat, factor, ipiv, error, unit)
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> LU factor
      real(wp), allocatable, intent(out) :: factor(:, :)
      !> Pivot indices
      integer, allocatable, intent(out) :: ipiv(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Output unit for the diagnostics; defaults to standard output
      integer, intent(in), optional :: unit

      integer :: n, info, iunit

      iunit = output_unit
      if (present(unit)) iunit = unit

      n = size(amat, 1)
      allocate (factor(n, n), ipiv(n))
      factor = amat
      call getrf(factor, ipiv, info)
      if (info /= 0) then
         write (iunit, "(A,I0)") "[factorize_pcm_lu] LAPACK getrf failed with info = ", info
         write (iunit, "(A,I0)") "[factorize_pcm_lu] Matrix size n = ", n
         call fatal_error(error, "[factorize_pcm_lu] LAPACK getrf failed")
         return
      end if

   end subroutine factorize_pcm_lu

   !> Solve several right-hand sides against an LU factor
   !>
   !> @param[in]  factor LU factor from [[factorize_pcm_lu]]
   !> @param[in]  ipiv   Pivot indices from [[factorize_pcm_lu]]
   !> @param[in]  rhs    Right-hand sides (ngrid, nrhs)
   !> @param[out] sol    Solutions (ngrid, nrhs)
   !> @param[out] error  Error handling
   subroutine solve_pcm_lu_factored(factor, ipiv, rhs, sol, error)
      !> LU factor
      real(wp), intent(in) :: factor(:, :)
      !> Pivot indices
      integer, intent(in) :: ipiv(:)
      !> Right-hand sides
      real(wp), intent(in) :: rhs(:, :)
      !> Solutions
      real(wp), intent(out) :: sol(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: info

      sol = rhs
      call getrs(factor, sol, ipiv, info)
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_lu_factored] LAPACK getrs failed")
         return
      end if

   end subroutine solve_pcm_lu_factored

   !> Cholesky-factorize the PCM matrix once, for repeated multi-column solves
   !>
   !> @param[in]  amat   System matrix (ngrid, ngrid), symmetric positive definite
   !> @param[out] factor Lower Cholesky factor (ngrid, ngrid)
   !> @param[out] error  Error handling
   subroutine factorize_pcm_cholesky(amat, factor, error)
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> Lower Cholesky factor
      real(wp), allocatable, intent(out) :: factor(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n, info

      n = size(amat, 1)
      allocate (factor(n, n))
      factor = amat
      call potrf(factor, info, uplo="l")
      if (info /= 0) then
         if (info > 0) then
            call fatal_error(error, "[factorize_pcm_cholesky] Matrix not positive definite")
         else
            call fatal_error(error, "[factorize_pcm_cholesky] LAPACK potrf failed")
         end if
         return
      end if

   end subroutine factorize_pcm_cholesky

   !> Solve several right-hand sides against a Cholesky factor
   !>
   !> @param[in]  factor Lower Cholesky factor from [[factorize_pcm_cholesky]]
   !> @param[in]  rhs    Right-hand sides (ngrid, nrhs)
   !> @param[out] sol    Solutions (ngrid, nrhs)
   !> @param[out] error  Error handling
   subroutine solve_pcm_cholesky_factored(factor, rhs, sol, error)
      !> Lower Cholesky factor
      real(wp), intent(in) :: factor(:, :)
      !> Right-hand sides
      real(wp), intent(in) :: rhs(:, :)
      !> Solutions
      real(wp), intent(out) :: sol(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: info

      sol = rhs
      call potrs(factor, sol, info, uplo="l")
      if (info /= 0) then
         call fatal_error(error, "[solve_pcm_cholesky_factored] LAPACK potrs failed")
         return
      end if

   end subroutine solve_pcm_cholesky_factored

   !> Invert the PCM matrix once, for repeated multi-column products
   !>
   !> @param[in]  amat   System matrix (ngrid, ngrid)
   !> @param[out] ainv   Inverse matrix (ngrid, ngrid)
   !> @param[out] error  Error handling
   subroutine invert_pcm_matrix(amat, ainv, error)
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> Inverse matrix
      real(wp), allocatable, intent(out) :: ainv(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: n, info
      integer, allocatable :: ipiv(:)

      n = size(amat, 1)
      allocate (ainv(n, n), ipiv(n))
      ainv = amat
      call getrf(ainv, ipiv, info)
      if (info /= 0) then
         call fatal_error(error, "[invert_pcm_matrix] LAPACK getrf failed")
         return
      end if
      call getri(ainv, ipiv, info=info)
      if (info /= 0) then
         call fatal_error(error, "[invert_pcm_matrix] LAPACK getri failed")
         return
      end if

   end subroutine invert_pcm_matrix

   !> Apply the inverse matrix to several right-hand sides
   !>
   !> @param[in]  ainv Inverse matrix from [[invert_pcm_matrix]]
   !> @param[in]  rhs  Right-hand sides (ngrid, nrhs)
   !> @param[out] sol  Solutions (ngrid, nrhs)
   subroutine apply_pcm_inverse(ainv, rhs, sol)
      !> Inverse matrix
      real(wp), intent(in) :: ainv(:, :)
      !> Right-hand sides
      real(wp), intent(in) :: rhs(:, :)
      !> Solutions
      real(wp), intent(out) :: sol(:, :)

      sol = 0.0_wp
      call gemm(ainv, rhs, sol)

   end subroutine apply_pcm_inverse

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
