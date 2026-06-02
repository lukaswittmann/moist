!> Provides wrappers for computing a matrix inverse

!> Wrappers to obtain the inverse of a matrix
module moist_math_lapack_getri
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: wrap_getri

   !> Computes the inverse of a matrix using the LU factorization
   !> computed by ?GETRF.
   !>
   !> This method inverts U and then computes inv(A) by solving the system
   !> inv(A)*L = inv(U) for inv(A).
   interface wrap_getri
      module procedure :: wrap_sgetri
      module procedure :: wrap_dgetri
   end interface wrap_getri

   !> Computes the inverse of a matrix using the LU factorization
   !> computed by ?GETRF.
   !>
   !> This method inverts U and then computes inv(A) by solving the system
   !> inv(A)*L = inv(U) for inv(A).
   interface lapack_getri
      pure subroutine sgetri(n, a, lda, ipiv, work, lwork, info)
         import :: sp, lapack_ik
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(in) :: ipiv(*)
         real(sp), intent(inout) :: work(*)
         integer(lapack_ik), intent(in) :: lwork
         integer(lapack_ik), intent(out) :: info
      end subroutine sgetri
      pure subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
         import :: dp, lapack_ik
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(in) :: ipiv(*)
         real(dp), intent(inout) :: work(*)
         integer(lapack_ik), intent(in) :: lwork
         integer(lapack_ik), intent(out) :: info
      end subroutine dgetri
   end interface lapack_getri

contains

   subroutine wrap_sgetri(amat, ipiv, info)
      real(sp), intent(inout) :: amat(:, :)
      integer, intent(in) :: ipiv(:)
      integer, intent(out) :: info
      integer(lapack_ik) :: n, lda, lwork, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      real(sp), allocatable :: work(:)
      real(sp) :: test(1)
      lda = int(max(1, size(amat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      lwork = int(-1, lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      ipiv_lapack = int(ipiv, lapack_ik)
      call lapack_getri(n, amat, lda, ipiv_lapack, test, lwork, info_lapack)
      info = int(info_lapack)
      if (info == 0) then
         lwork = nint(test(1), kind=lapack_ik)
         allocate (work(int(lwork)))
         call lapack_getri(n, amat, lda, ipiv_lapack, work, lwork, info_lapack)
         info = int(info_lapack)
      end if
   end subroutine wrap_sgetri

   subroutine wrap_dgetri(amat, ipiv, info)
      real(dp), intent(inout) :: amat(:, :)
      integer, intent(in) :: ipiv(:)
      integer, intent(out) :: info
      integer(lapack_ik) :: n, lda, lwork, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      real(dp), allocatable :: work(:)
      real(dp) :: test(1)
      lda = int(max(1, size(amat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      lwork = int(-1, lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      ipiv_lapack = int(ipiv, lapack_ik)
      call lapack_getri(n, amat, lda, ipiv_lapack, test, lwork, info_lapack)
      info = int(info_lapack)
      if (info == 0) then
         lwork = nint(test(1), kind=lapack_ik)
         allocate (work(int(lwork)))
         call lapack_getri(n, amat, lda, ipiv_lapack, work, lwork, info_lapack)
         info = int(info_lapack)
      end if
   end subroutine wrap_dgetri

end module moist_math_lapack_getri
