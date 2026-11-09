!> Provides wrappers for LU factorization routines

!> Wrapper routines for LU factorization
module moist_math_lapack_getrf
   use mctc_env, only: dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none(type, external)
   private

   public :: wrap_getrf

   !> LU factorization with partial pivoting of a general M-by-N matrix A
   !>
   !> The factorization has the form
   !>    A = P * L * U
   !> where P is a permutation matrix, L is lower triangular with unit
   !> diagonal elements (lower trapezoidal if m > n), and U is upper
   !> triangular (upper trapezoidal if m < n). Real and complex matrices.
   interface wrap_getrf
      module procedure :: wrap_dgetrf
      module procedure :: wrap_zgetrf
   end interface wrap_getrf

   !> LAPACK LU factorization with partial pivoting of a general matrix
   interface lapack_getrf
      !> Double precision LU factorization
      !>
      !> @param[in]     m     number of rows of a
      !> @param[in]     n     number of columns of a
      !> @param[in,out] a     matrix on input, factors on output
      !> @param[in]     lda   leading dimension of a
      !> @param[out]    ipiv  pivot indices
      !> @param[out]    info  zero on success, positive if u is exactly singular
      pure subroutine dgetrf(m, n, a, lda, ipiv, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         !> Number of rows of a
         integer(lapack_ik), intent(in) :: m
         !> Number of columns of a
         integer(lapack_ik), intent(in) :: n
         !> Leading dimension of a
         integer(lapack_ik), intent(in) :: lda
         !> Matrix and factors
         real(dp), intent(inout) :: a(lda, *)
         !> Pivot indices
         integer(lapack_ik), intent(out) :: ipiv(*)
         !> Status
         integer(lapack_ik), intent(out) :: info
      end subroutine dgetrf
      !> Double complex LU factorization
      !>
      !> @param[in]     m     number of rows of a
      !> @param[in]     n     number of columns of a
      !> @param[in,out] a     matrix on input, factors on output
      !> @param[in]     lda   leading dimension of a
      !> @param[out]    ipiv  pivot indices
      !> @param[out]    info  zero on success, positive if u is exactly singular
      pure subroutine zgetrf(m, n, a, lda, ipiv, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         !> Number of rows of a
         integer(lapack_ik), intent(in) :: m
         !> Number of columns of a
         integer(lapack_ik), intent(in) :: n
         !> Leading dimension of a
         integer(lapack_ik), intent(in) :: lda
         !> Matrix and factors
         complex(dp), intent(inout) :: a(lda, *)
         !> Pivot indices
         integer(lapack_ik), intent(out) :: ipiv(*)
         !> Status
         integer(lapack_ik), intent(out) :: info
      end subroutine zgetrf
   end interface lapack_getrf

contains

!> Double precision LU factorization
!>
!> @param[in,out] amat  matrix on input, factors on output, (m, n)
!> @param[out]    ipiv  pivot indices, (min(m, n))
!> @param[out]    info  zero on success, positive if u is exactly singular
   subroutine wrap_dgetrf(amat, ipiv, info)
      !> Matrix and factors
      real(dp), intent(inout) :: amat(:, :)
      !> Pivot indices
      integer, intent(out) :: ipiv(:)
      !> Status
      integer, intent(out) :: info
      integer(lapack_ik) :: m, n, lda, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      lda = int(max(1, size(amat, 1)), lapack_ik)
      m = int(size(amat, 1), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      call lapack_getrf(m, n, amat, lda, ipiv_lapack, info_lapack)
      ipiv = int(ipiv_lapack)
      info = int(info_lapack)
   end subroutine wrap_dgetrf

!> Double complex LU factorization
!>
!> @param[in,out] amat  matrix on input, factors on output, (m, n)
!> @param[out]    ipiv  pivot indices, (min(m, n))
!> @param[out]    info  zero on success, positive if u is exactly singular
   subroutine wrap_zgetrf(amat, ipiv, info)
      !> Matrix and factors
      complex(dp), intent(inout) :: amat(:, :)
      !> Pivot indices
      integer, intent(out) :: ipiv(:)
      !> Status
      integer, intent(out) :: info
      integer(lapack_ik) :: m, n, lda, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      lda = int(max(1, size(amat, 1)), lapack_ik)
      m = int(size(amat, 1), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      call lapack_getrf(m, n, amat, lda, ipiv_lapack, info_lapack)
      ipiv = int(ipiv_lapack)
      info = int(info_lapack)
   end subroutine wrap_zgetrf

end module moist_math_lapack_getrf
