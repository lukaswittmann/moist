!> Small dense matrix inversion and analytic eigen-decomposition.
!>
!> Analytical routines for small fixed-size matrices that are not efficiently
!> handled by BLAS/LAPACK due to call overhead.
module moist_math_linalg_decomp
   use mctc_env_accuracy, only: wp
   implicit none
   private

   public :: mat3x3_inv
   public :: eig_2x2_symmetric

contains

   !> Compute analytical inverse of a 3x3 matrix
   !>
   !> Uses the adjugate matrix formula: A^{-1} = adj(A) / det(A)
   !> Returns .false. if matrix is singular (within tolerance).
   !> For 3x3 matrices, this analytical approach is faster than
   !> LAPACK LU decomposition due to lower overhead.
   !>
   !> @param[in]  A       Input 3x3 matrix to invert
   !> @param[out] Ainv    Output 3x3 inverse matrix
   !> @param[out] det     Determinant of A (optional)
   !> @param[in]  tol     Singularity tolerance (optional, default: 1.0e-30)
   !> @returns    success .true. if inversion succeeded, .false. if singular
   function mat3x3_inv(A, Ainv, det, tol) result(success)
      !> Input matrix [3,3]
      real(wp), intent(in) :: A(3, 3)
      !> Output inverse matrix [3,3]
      real(wp), intent(out) :: Ainv(3, 3)
      !> Determinant of A (optional output)
      real(wp), intent(out), optional :: det
      !> Singularity tolerance (optional, default 1.0e-30)
      real(wp), intent(in), optional :: tol
      !> Success flag
      logical :: success

      real(wp) :: a11, a12, a13, a21, a22, a23, a31, a32, a33
      real(wp) :: det_val, inv_det, tol_val
      real(wp) :: c11, c12, c13, c21, c22, c23, c31, c32, c33

      ! Set tolerance
      if (present(tol)) then
         tol_val = tol
      else
         tol_val = 1.0e-30_wp
      end if

      ! Extract matrix elements for clarity
      a11 = A(1, 1); a12 = A(1, 2); a13 = A(1, 3)
      a21 = A(2, 1); a22 = A(2, 2); a23 = A(2, 3)
      a31 = A(3, 1); a32 = A(3, 2); a33 = A(3, 3)

      ! Compute cofactors (each used for both determinant and inverse)
      c11 = a22*a33 - a23*a32
      c12 = a23*a31 - a21*a33
      c13 = a21*a32 - a22*a31
      c21 = a13*a32 - a12*a33
      c22 = a11*a33 - a13*a31
      c23 = a12*a31 - a11*a32
      c31 = a12*a23 - a13*a22
      c32 = a13*a21 - a11*a23
      c33 = a11*a22 - a12*a21

      ! Compute determinant via cofactor expansion along first row
      det_val = a11*c11 + a12*c12 + a13*c13

      ! Check for singularity
      if (abs(det_val) < tol_val) then
         success = .false.
         Ainv = 0.0_wp
         if (present(det)) det = det_val
         return
      end if

      success = .true.
      if (present(det)) det = det_val
      inv_det = 1.0_wp/det_val

      ! Compute adjugate matrix (transposed cofactors) divided by determinant
      ! Row 1
      Ainv(1, 1) = c11*inv_det
      Ainv(1, 2) = c21*inv_det
      Ainv(1, 3) = c31*inv_det
      ! Row 2
      Ainv(2, 1) = c12*inv_det
      Ainv(2, 2) = c22*inv_det
      Ainv(2, 3) = c32*inv_det
      ! Row 3
      Ainv(3, 1) = c13*inv_det
      Ainv(3, 2) = c23*inv_det
      Ainv(3, 3) = c33*inv_det

   end function mat3x3_inv

   !> Analytic eigenvalue decomposition for 2x2 symmetric matrix
   !>
   !> For a symmetric 2x2 matrix R = [a b; b c], computes eigenvalues and
   !> eigenvectors analytically using the characteristic polynomial. This is
   !> significantly faster than calling LAPACK for such small matrices due to
   !> elimination of function call overhead and optimized branch prediction.
   !>
   !> The eigenvalues are computed via the quadratic formula applied to
   !> the characteristic equation: lambda^2 - trace*lambda + det = 0
   !>
   !> Eigenvectors are computed using the standard formula (R - lambda*I)v = 0,
   !> with special handling for diagonal and degenerate cases.
   !>
   !> @param[in]  a          R(1,1) element
   !> @param[in]  b          R(1,2) = R(2,1) off-diagonal element
   !> @param[in]  c          R(2,2) element
   !> @param[out] lambda_min Smallest eigenvalue
   !> @param[out] lambda_max Largest eigenvalue
   !> @param[out] v_min      Eigenvector for smallest eigenvalue (normalized)
   !> @param[out] v_max      Eigenvector for largest eigenvalue (normalized)
   pure subroutine eig_2x2_symmetric(a, b, c, lambda_min, lambda_max, v_min, v_max)
      !> R(1,1) matrix element
      real(wp), intent(in)  :: a
      !> R(1,2) = R(2,1) off-diagonal element
      real(wp), intent(in)  :: b
      !> R(2,2) matrix element
      real(wp), intent(in)  :: c
      !> Smallest eigenvalue
      real(wp), intent(out) :: lambda_min
      !> Largest eigenvalue
      real(wp), intent(out) :: lambda_max
      !> Eigenvector for smallest eigenvalue
      real(wp), intent(out) :: v_min(2)
      !> Eigenvector for largest eigenvalue
      real(wp), intent(out) :: v_max(2)

      real(wp) :: trace, det, disc, sqrt_disc, norm

      ! Characteristic polynomial: lambda^2 - trace*lambda + det = 0
      trace = a + c
      det = a*c - b*b
      disc = trace*trace - 4.0_wp*det

      ! Eigenvalues (guaranteed real for symmetric matrix)
      if (disc < 0.0_wp) then
         ! Numerically negative discriminant; clamp to zero
         sqrt_disc = 0.0_wp
      else
         sqrt_disc = sqrt(disc)
      end if

      lambda_min = 0.5_wp*(trace - sqrt_disc)
      lambda_max = 0.5_wp*(trace + sqrt_disc)

      ! Eigenvectors: For eigenvalue lambda, solve (R - lambda*I)*v = 0
      ! First row: (a - lambda)*v1 + b*v2 = 0  =>  v2 = -(a - lambda)*v1 / b
      ! Then normalize. Handle special cases when b ~= 0 (diagonal matrix)

      if (abs(b) > 1.0e-14_wp) then
         ! Non-diagonal case: use standard formula
         ! For lambda_min:
         v_min(1) = b
         v_min(2) = lambda_min - a
         norm = sqrt(v_min(1)**2 + v_min(2)**2)
         if (norm > 1.0e-14_wp) then
            v_min = v_min/norm
         else
            ! Degenerate: set to canonical basis
            v_min = [1.0_wp, 0.0_wp]
         end if

         ! For lambda_max:
         v_max(1) = b
         v_max(2) = lambda_max - a
         norm = sqrt(v_max(1)**2 + v_max(2)**2)
         if (norm > 1.0e-14_wp) then
            v_max = v_max/norm
         else
            v_max = [0.0_wp, 1.0_wp]
         end if
      else
         ! Diagonal (or nearly diagonal) case
         if (abs(a - c) < 1.0e-14_wp) then
            ! Scalar multiple of identity: eigenspaces are degenerate
            v_min = [1.0_wp, 0.0_wp]
            v_max = [0.0_wp, 1.0_wp]
         else if (a < c) then
            ! a is smaller eigenvalue
            v_min = [1.0_wp, 0.0_wp]
            v_max = [0.0_wp, 1.0_wp]
         else
            ! c is smaller eigenvalue
            v_min = [0.0_wp, 1.0_wp]
            v_max = [1.0_wp, 0.0_wp]
         end if
      end if

   end subroutine eig_2x2_symmetric

end module moist_math_linalg_decomp
