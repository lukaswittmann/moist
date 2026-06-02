!> Provides wrappers to solve a linear equation system using Cholesky factorization

!> Wrappers to solve a system of linear equations using Cholesky factorization
module moist_math_lapack_potrs
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: wrap_potrs

   !> Solves a system of linear equations
   !>    A * X = B
   !> where A is a real symmetric positive definite matrix using the
   !> Cholesky factorization computed by ?POTRF.
   interface wrap_potrs
      module procedure :: wrap_spotrs
      module procedure :: wrap_dpotrs
   end interface wrap_potrs

   !> Solves a system of linear equations
   !>    A * X = B
   !> where A is a real symmetric positive definite matrix using the
   !> Cholesky factorization A = U**T * U or A = L * L**T computed by ?POTRF.
   interface lapack_potrs
      pure subroutine spotrs(uplo, n, nrhs, a, lda, b, ldb, info)
         import :: sp, lapack_ik
         character(len=1), intent(in) :: uplo
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: nrhs
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldb
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         integer(lapack_ik), intent(out) :: info
      end subroutine spotrs
      pure subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
         import :: dp, lapack_ik
         character(len=1), intent(in) :: uplo
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: nrhs
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldb
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         integer(lapack_ik), intent(out) :: info
      end subroutine dpotrs
   end interface lapack_potrs

contains

   subroutine wrap_spotrs(amat, bmat, info, uplo)
      real(sp), intent(in) :: amat(:, :)
      real(sp), intent(inout) :: bmat(:, :)
      integer, intent(out) :: info
      character(len=1), intent(in), optional :: uplo
      character(len=1) :: ula
      integer(lapack_ik) :: n, nrhs, lda, ldb, info_lapack
      if (present(uplo)) then
         ula = uplo
      else
         ula = 'u'
      end if
      lda = int(max(1, size(amat, 1)), lapack_ik)
      ldb = int(max(1, size(bmat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      nrhs = int(size(bmat, 2), lapack_ik)
      call lapack_potrs(ula, n, nrhs, amat, lda, bmat, ldb, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_spotrs

   subroutine wrap_dpotrs(amat, bmat, info, uplo)
      real(dp), intent(in) :: amat(:, :)
      real(dp), intent(inout) :: bmat(:, :)
      integer, intent(out) :: info
      character(len=1), intent(in), optional :: uplo
      character(len=1) :: ula
      integer(lapack_ik) :: n, nrhs, lda, ldb, info_lapack
      if (present(uplo)) then
         ula = uplo
      else
         ula = 'u'
      end if
      lda = int(max(1, size(amat, 1)), lapack_ik)
      ldb = int(max(1, size(bmat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      nrhs = int(size(bmat, 2), lapack_ik)
      call lapack_potrs(ula, n, nrhs, amat, lda, bmat, ldb, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_dpotrs

end module moist_math_lapack_potrs
