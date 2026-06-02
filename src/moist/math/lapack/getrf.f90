!> Provides wrappers for LU factorization routines

!> Wrapper rountines for LU factorization
module moist_math_lapack_getrf
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: wrap_getrf

   !> Computes an LU factorization of a general M-by-N matrix A
   !> using partial pivoting with row interchanges.
   !>
   !> The factorization has the form
   !>    A = P * L * U
   !> where P is a permutation matrix, L is lower triangular with unit
   !> diagonal elements (lower trapezoidal if m > n), and U is upper
   !> triangular (upper trapezoidal if m < n).
   interface wrap_getrf
      module procedure :: wrap_sgetrf
      module procedure :: wrap_dgetrf
   end interface wrap_getrf

   !> Computes an LU factorization of a general M-by-N matrix A
   !> using partial pivoting with row interchanges.
   !>
   !> The factorization has the form
   !>    A = P * L * U
   !> where P is a permutation matrix, L is lower triangular with unit
   !> diagonal elements (lower trapezoidal if m > n), and U is upper
   !> triangular (upper trapezoidal if m < n).
   interface lapack_getrf
      pure subroutine sgetrf(m, n, a, lda, ipiv, info)
         import :: sp, lapack_ik
         integer(lapack_ik), intent(in) :: m
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(out) :: ipiv(*)
         integer(lapack_ik), intent(out) :: info
      end subroutine sgetrf
      pure subroutine dgetrf(m, n, a, lda, ipiv, info)
         import :: dp, lapack_ik
         integer(lapack_ik), intent(in) :: m
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(out) :: ipiv(*)
         integer(lapack_ik), intent(out) :: info
      end subroutine dgetrf
   end interface lapack_getrf

contains

   subroutine wrap_sgetrf(amat, ipiv, info)
      real(sp), intent(inout) :: amat(:, :)
      integer, intent(out) :: ipiv(:)
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
   end subroutine wrap_sgetrf

   subroutine wrap_dgetrf(amat, ipiv, info)
      real(dp), intent(inout) :: amat(:, :)
      integer, intent(out) :: ipiv(:)
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

end module moist_math_lapack_getrf
