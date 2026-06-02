!> Provides wrappers to solve a linear equation system

!> Wrappers to solve a system of linear equations
module moist_math_lapack_getrs
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: wrap_getrs

   !> Solves a system of linear equations
   !>    A * X = B  or  A**T * X = B
   !> with a general N-by-N matrix A using the LU factorization computed
   !> by ?GETRF.
   interface wrap_getrs
      module procedure :: wrap_sgetrs
      module procedure :: wrap_dgetrs
   end interface wrap_getrs

   !> Solves a system of linear equations
   !>    A * X = B  or  A**T * X = B
   !> with a general N-by-N matrix A using the LU factorization computed
   !> by ?GETRF.
   interface lapack_getrs
      pure subroutine sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp, lapack_ik
         character(len=1), intent(in) :: trans
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: nrhs
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldb
         real(sp), intent(in) :: a(lda, *)
         integer(lapack_ik), intent(in) :: ipiv(*)
         real(sp), intent(inout) :: b(ldb, *)
         integer(lapack_ik), intent(out) :: info
      end subroutine sgetrs
      pure subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, lapack_ik
         character(len=1), intent(in) :: trans
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: nrhs
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldb
         real(dp), intent(in) :: a(lda, *)
         integer(lapack_ik), intent(in) :: ipiv(*)
         real(dp), intent(inout) :: b(ldb, *)
         integer(lapack_ik), intent(out) :: info
      end subroutine dgetrs
   end interface lapack_getrs

contains

   subroutine wrap_sgetrs(amat, bmat, ipiv, info, trans)
      real(sp), intent(in) :: amat(:, :)
      real(sp), intent(inout) :: bmat(:, :)
      integer, intent(in) :: ipiv(:)
      integer, intent(out) :: info
      character(len=1), intent(in), optional :: trans
      character(len=1) :: tra
      integer(lapack_ik) :: n, nrhs, lda, ldb, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      if (present(trans)) then
         tra = trans
      else
         tra = 'n'
      end if
      lda = int(max(1, size(amat, 1)), lapack_ik)
      ldb = int(max(1, size(bmat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      nrhs = int(size(bmat, 2), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      ipiv_lapack = int(ipiv, lapack_ik)
      call lapack_getrs(tra, n, nrhs, amat, lda, ipiv_lapack, bmat, ldb, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_sgetrs

   subroutine wrap_dgetrs(amat, bmat, ipiv, info, trans)
      real(dp), intent(in) :: amat(:, :)
      real(dp), intent(inout) :: bmat(:, :)
      integer, intent(in) :: ipiv(:)
      integer, intent(out) :: info
      character(len=1), intent(in), optional :: trans
      character(len=1) :: tra
      integer(lapack_ik) :: n, nrhs, lda, ldb, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      if (present(trans)) then
         tra = trans
      else
         tra = 'n'
      end if
      lda = int(max(1, size(amat, 1)), lapack_ik)
      ldb = int(max(1, size(bmat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      nrhs = int(size(bmat, 2), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      ipiv_lapack = int(ipiv, lapack_ik)
      call lapack_getrs(tra, n, nrhs, amat, lda, ipiv_lapack, bmat, ldb, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_dgetrs

end module moist_math_lapack_getrs
