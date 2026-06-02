!> Provides interface to general linear system solver

!> Interface to LAPACK general linear system solver
module moist_math_lapack_gesv
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: lapack_gesv
   public :: dgesv, sgesv

   !> Computes the solution to a real system of linear equations
   !>    A * X = B,
   !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
   !> The LU decomposition with partial pivoting and row interchanges is
   !> used to factor A as
   !>    A = P * L * U,
   !> where P is a permutation matrix, L is unit lower triangular, and U is
   !> upper triangular.  The factored form of A is then used to solve the
   !> system of equations A * X = B.
   interface lapack_gesv
      pure subroutine sgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp, lapack_ik
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: nrhs
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldb
         real(sp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(out) :: ipiv(*)
         real(sp), intent(inout) :: b(ldb, *)
         integer(lapack_ik), intent(out) :: info
      end subroutine sgesv
      pure subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, lapack_ik
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: nrhs
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldb
         real(dp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(out) :: ipiv(*)
         real(dp), intent(inout) :: b(ldb, *)
         integer(lapack_ik), intent(out) :: info
      end subroutine dgesv
   end interface lapack_gesv

end module moist_math_lapack_gesv
