!> Computes the Cholesky factorization of a real symmetric positive definite matrix A.
module moist_math_lapack_potrf
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: wrap_potrf

   !> Computes the Cholesky factorization of a real symmetric
   !> positive definite matrix A.
   !>
   !> The factorization has the form
   !>    A = U**T * U,  if UPLO = 'U', or
   !>    A = L  * L**T,  if UPLO = 'L',
   !> where U is an upper triangular matrix and L is lower triangular.
   !>
   !> This is the block version of the algorithm, calling Level 3 BLAS.
   interface wrap_potrf
      module procedure :: wrap_spotrf
      module procedure :: wrap_dpotrf
   end interface wrap_potrf

   !> Computes the Cholesky factorization of a real symmetric
   !> positive definite matrix A.
   !>
   !> The factorization has the form
   !>    A = U**T * U,  if UPLO = 'U', or
   !>    A = L  * L**T,  if UPLO = 'L',
   !> where U is an upper triangular matrix and L is lower triangular.
   !>
   !> This is the block version of the algorithm, calling Level 3 BLAS.
   interface lapack_potrf
      pure subroutine spotrf(uplo, n, a, lda, info)
         import :: sp, lapack_ik
         character(len=1), intent(in) :: uplo
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(out) :: info
      end subroutine spotrf
      pure subroutine dpotrf(uplo, n, a, lda, info)
         import :: dp, lapack_ik
         character(len=1), intent(in) :: uplo
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(out) :: info
      end subroutine dpotrf
   end interface lapack_potrf

contains

   subroutine wrap_spotrf(amat, info, uplo)
      real(sp), intent(inout) :: amat(:, :)
      integer, intent(out) :: info
      character(len=1), intent(in), optional :: uplo
      integer(lapack_ik) :: n, lda, info_lapack
      character(len=1) :: ula

      ula = 'u'
      if (present(uplo)) ula = uplo
      lda = int(max(1, size(amat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      call lapack_potrf(ula, n, amat, lda, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_spotrf

   subroutine wrap_dpotrf(amat, info, uplo)
      real(dp), intent(inout) :: amat(:, :)
      integer, intent(out) :: info
      character(len=1), intent(in), optional :: uplo
      integer(lapack_ik) :: n, lda, info_lapack
      character(len=1) :: ula

      ula = 'u'
      if (present(uplo)) ula = uplo
      lda = int(max(1, size(amat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      call lapack_potrf(ula, n, amat, lda, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_dpotrf

end module moist_math_lapack_potrf
