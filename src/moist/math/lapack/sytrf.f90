!> Provides wrappers for symmetric indefinite factorization routines

!> Wrapper routines for the Bunch-Kaufman factorization of a symmetric matrix
module moist_math_lapack_sytrf
   use mctc_env, only: dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none(type, external)
   private

   public :: wrap_sytrf

   !> Bunch-Kaufman factorization of a real or complex symmetric matrix A
   !>
   !> Complex matrices are symmetric (A = A**T), not Hermitian.
   !> The factorization has the form
   !>    A = U**T * D * U  or  A = L * D * L**T
   !> where U (or L) is a product of permutation and unit upper (lower)
   !> triangular matrices, and D is symmetric and block diagonal with
   !> 1-by-1 and 2-by-2 diagonal blocks.
   !>
   !> The workspace is queried and allocated internally.
   interface wrap_sytrf
      module procedure :: wrap_dsytrf
      module procedure :: wrap_zsytrf
   end interface wrap_sytrf

   !> LAPACK Bunch-Kaufman factorization of a real or complex symmetric matrix
   interface lapack_sytrf
      !> Double precision symmetric indefinite factorization
      !>
      !> @param[in]     uplo   "U" or "L"
      !> @param[in]     n      order of a
      !> @param[in,out] a      matrix on input, factors on output
      !> @param[in]     lda    leading dimension of a
      !> @param[out]    ipiv   pivot indices
      !> @param[out]    work   workspace; work(1) is the optimal lwork on return
      !> @param[in]     lwork  size of work, -1 for a workspace query
      !> @param[out]    info   zero on success, positive if d is exactly singular
      pure subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         !> Triangle of a that is stored
         character(len=1), intent(in) :: uplo
         !> Order of a
         integer(lapack_ik), intent(in) :: n
         !> Leading dimension of a
         integer(lapack_ik), intent(in) :: lda
         !> Matrix and factors
         real(dp), intent(inout) :: a(lda, *)
         !> Pivot indices
         integer(lapack_ik), intent(out) :: ipiv(*)
         !> Workspace
         real(dp), intent(out) :: work(*)
         !> Size of work
         integer(lapack_ik), intent(in) :: lwork
         !> Status
         integer(lapack_ik), intent(out) :: info
      end subroutine dsytrf
      !> Double complex symmetric indefinite factorization
      !>
      !> @param[in]     uplo   "U" or "L"
      !> @param[in]     n      order of a
      !> @param[in,out] a      matrix on input, factors on output
      !> @param[in]     lda    leading dimension of a
      !> @param[out]    ipiv   pivot indices
      !> @param[out]    work   workspace; work(1) is the optimal lwork on return
      !> @param[in]     lwork  size of work, -1 for a workspace query
      !> @param[out]    info   zero on success, positive if d is exactly singular
      pure subroutine zsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         !> Triangle of a that is stored
         character(len=1), intent(in) :: uplo
         !> Order of a
         integer(lapack_ik), intent(in) :: n
         !> Leading dimension of a
         integer(lapack_ik), intent(in) :: lda
         !> Matrix and factors
         complex(dp), intent(inout) :: a(lda, *)
         !> Pivot indices
         integer(lapack_ik), intent(out) :: ipiv(*)
         !> Workspace
         complex(dp), intent(out) :: work(*)
         !> Size of work
         integer(lapack_ik), intent(in) :: lwork
         !> Status
         integer(lapack_ik), intent(out) :: info
      end subroutine zsytrf
   end interface lapack_sytrf

contains

!> Double precision Bunch-Kaufman factorization with internal workspace query
!>
!> @param[in,out] amat  symmetric matrix on input, factors on output, (n, n)
!> @param[out]    ipiv  pivot indices, (n)
!> @param[out]    info  zero on success, positive if d is exactly singular
!> @param[in]     uplo  triangle of amat that is stored, "u" (default) or "l"
   subroutine wrap_dsytrf(amat, ipiv, info, uplo)
      !> Symmetric matrix and factors
      real(dp), intent(inout) :: amat(:, :)
      !> Pivot indices
      integer, intent(out) :: ipiv(:)
      !> Status
      integer, intent(out) :: info
      !> Triangle of amat that is stored
      character(len=1), intent(in), optional :: uplo
      character(len=1) :: ula
      integer(lapack_ik) :: n, lda, lwork, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      real(dp) :: query(1)
      real(dp), allocatable :: work(:)

      ula = "u"
      if (present(uplo)) ula = uplo
      lda = int(max(1, size(amat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      call lapack_sytrf(ula, n, amat, lda, ipiv_lapack, query, -1_lapack_ik, info_lapack)
      if (info_lapack /= 0) then
         info = int(info_lapack)
         return
      end if
      lwork = int(max(1, int(query(1))), lapack_ik)
      allocate (work(lwork))
      call lapack_sytrf(ula, n, amat, lda, ipiv_lapack, work, lwork, info_lapack)
      ipiv = int(ipiv_lapack)
      info = int(info_lapack)
   end subroutine wrap_dsytrf

!> Double complex Bunch-Kaufman factorization with internal workspace query
!>
!> @param[in,out] amat  symmetric matrix on input, factors on output, (n, n)
!> @param[out]    ipiv  pivot indices, (n)
!> @param[out]    info  zero on success, positive if d is exactly singular
!> @param[in]     uplo  triangle of amat that is stored, "u" (default) or "l"
   subroutine wrap_zsytrf(amat, ipiv, info, uplo)
      !> Symmetric matrix and factors
      complex(dp), intent(inout) :: amat(:, :)
      !> Pivot indices
      integer, intent(out) :: ipiv(:)
      !> Status
      integer, intent(out) :: info
      !> Triangle of amat that is stored
      character(len=1), intent(in), optional :: uplo
      character(len=1) :: ula
      integer(lapack_ik) :: n, lda, lwork, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      complex(dp) :: query(1)
      complex(dp), allocatable :: work(:)

      ula = "u"
      if (present(uplo)) ula = uplo
      lda = int(max(1, size(amat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      call lapack_sytrf(ula, n, amat, lda, ipiv_lapack, query, -1_lapack_ik, info_lapack)
      if (info_lapack /= 0) then
         info = int(info_lapack)
         return
      end if
      lwork = int(max(1, int(real(query(1)))), lapack_ik)
      allocate (work(lwork))
      call lapack_sytrf(ula, n, amat, lda, ipiv_lapack, work, lwork, info_lapack)
      ipiv = int(ipiv_lapack)
      info = int(info_lapack)
   end subroutine wrap_zsytrf

end module moist_math_lapack_sytrf
