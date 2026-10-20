!> Provides wrappers to solve a linear equation system

!> Wrappers to solve a system of linear equations
module moist_math_lapack_getrs
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none(type, external)
   private

   public :: wrap_getrs

   !> Solves A * X = B or A**T * X = B with LU factors from ?GETRF
   !>
   !> A is a general N-by-N real or complex matrix. trans = "t" is the plain
   !> transpose also for complex A; "c" gives the conjugate transpose.
   interface wrap_getrs
      module procedure :: wrap_sgetrs
      module procedure :: wrap_dgetrs
      module procedure :: wrap_zgetrs
   end interface wrap_getrs

   !> LAPACK solve with LU factors of a general matrix
   interface lapack_getrs
      !> Single precision solve with LU factors
      !>
      !> @param[in]     trans  "N", "T" or "C"
      !> @param[in]     n      order of a
      !> @param[in]     nrhs   number of right-hand sides
      !> @param[in]     a      factors from sgetrf
      !> @param[in]     lda    leading dimension of a
      !> @param[in]     ipiv   pivot indices from sgetrf
      !> @param[in,out] b      right-hand sides on input, solutions on output
      !> @param[in]     ldb    leading dimension of b
      !> @param[out]    info   zero on success
      pure subroutine sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp, lapack_ik
         implicit none(type, external)
         !> Orientation of the system
         character(len=1), intent(in) :: trans
         !> Order of a
         integer(lapack_ik), intent(in) :: n
         !> Number of right-hand sides
         integer(lapack_ik), intent(in) :: nrhs
         !> Leading dimension of a
         integer(lapack_ik), intent(in) :: lda
         !> Leading dimension of b
         integer(lapack_ik), intent(in) :: ldb
         !> Factors
         real(sp), intent(in) :: a(lda, *)
         !> Pivot indices
         integer(lapack_ik), intent(in) :: ipiv(*)
         !> Right-hand sides and solutions
         real(sp), intent(inout) :: b(ldb, *)
         !> Status
         integer(lapack_ik), intent(out) :: info
      end subroutine sgetrs
      !> Double precision solve with LU factors
      !>
      !> @param[in]     trans  "N", "T" or "C"
      !> @param[in]     n      order of a
      !> @param[in]     nrhs   number of right-hand sides
      !> @param[in]     a      factors from dgetrf
      !> @param[in]     lda    leading dimension of a
      !> @param[in]     ipiv   pivot indices from dgetrf
      !> @param[in,out] b      right-hand sides on input, solutions on output
      !> @param[in]     ldb    leading dimension of b
      !> @param[out]    info   zero on success
      pure subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         !> Orientation of the system
         character(len=1), intent(in) :: trans
         !> Order of a
         integer(lapack_ik), intent(in) :: n
         !> Number of right-hand sides
         integer(lapack_ik), intent(in) :: nrhs
         !> Leading dimension of a
         integer(lapack_ik), intent(in) :: lda
         !> Leading dimension of b
         integer(lapack_ik), intent(in) :: ldb
         !> Factors
         real(dp), intent(in) :: a(lda, *)
         !> Pivot indices
         integer(lapack_ik), intent(in) :: ipiv(*)
         !> Right-hand sides and solutions
         real(dp), intent(inout) :: b(ldb, *)
         !> Status
         integer(lapack_ik), intent(out) :: info
      end subroutine dgetrs
      !> Double complex solve with LU factors
      !>
      !> @param[in]     trans  "N", "T" (plain transpose) or "C" (conjugate transpose)
      !> @param[in]     n      order of a
      !> @param[in]     nrhs   number of right-hand sides
      !> @param[in]     a      factors from zgetrf
      !> @param[in]     lda    leading dimension of a
      !> @param[in]     ipiv   pivot indices from zgetrf
      !> @param[in,out] b      right-hand sides on input, solutions on output
      !> @param[in]     ldb    leading dimension of b
      !> @param[out]    info   zero on success
      pure subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         !> Orientation of the system
         character(len=1), intent(in) :: trans
         !> Order of a
         integer(lapack_ik), intent(in) :: n
         !> Number of right-hand sides
         integer(lapack_ik), intent(in) :: nrhs
         !> Leading dimension of a
         integer(lapack_ik), intent(in) :: lda
         !> Leading dimension of b
         integer(lapack_ik), intent(in) :: ldb
         !> Factors
         complex(dp), intent(in) :: a(lda, *)
         !> Pivot indices
         integer(lapack_ik), intent(in) :: ipiv(*)
         !> Right-hand sides and solutions
         complex(dp), intent(inout) :: b(ldb, *)
         !> Status
         integer(lapack_ik), intent(out) :: info
      end subroutine zgetrs
   end interface lapack_getrs

contains

!> Single precision solve with factors from wrap_getrf
!>
!> @param[in]     amat   factors from wrap_getrf, (n, n)
!> @param[in,out] bmat   right-hand sides on input, solutions on output, (n, nrhs)
!> @param[in]     ipiv   pivot indices from wrap_getrf, (n)
!> @param[out]    info   zero on success
!> @param[in]     trans  "n" (default), "t" or "c"
   subroutine wrap_sgetrs(amat, bmat, ipiv, info, trans)
      !> Factors
      real(sp), intent(in) :: amat(:, :)
      !> Right-hand sides and solutions
      real(sp), intent(inout) :: bmat(:, :)
      !> Pivot indices
      integer, intent(in) :: ipiv(:)
      !> Status
      integer, intent(out) :: info
      !> Orientation of the system
      character(len=1), intent(in), optional :: trans
      character(len=1) :: tra
      integer(lapack_ik) :: n, nrhs, lda, ldb, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      if (present(trans)) then
         tra = trans
      else
         tra = "n"
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

!> Double precision solve with factors from wrap_getrf
!>
!> @param[in]     amat   factors from wrap_getrf, (n, n)
!> @param[in,out] bmat   right-hand sides on input, solutions on output, (n, nrhs)
!> @param[in]     ipiv   pivot indices from wrap_getrf, (n)
!> @param[out]    info   zero on success
!> @param[in]     trans  "n" (default), "t" or "c"
   subroutine wrap_dgetrs(amat, bmat, ipiv, info, trans)
      !> Factors
      real(dp), intent(in) :: amat(:, :)
      !> Right-hand sides and solutions
      real(dp), intent(inout) :: bmat(:, :)
      !> Pivot indices
      integer, intent(in) :: ipiv(:)
      !> Status
      integer, intent(out) :: info
      !> Orientation of the system
      character(len=1), intent(in), optional :: trans
      character(len=1) :: tra
      integer(lapack_ik) :: n, nrhs, lda, ldb, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      if (present(trans)) then
         tra = trans
      else
         tra = "n"
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

!> Double complex solve with factors from wrap_getrf
!>
!> @param[in]     amat   factors from wrap_getrf, (n, n)
!> @param[in,out] bmat   right-hand sides on input, solutions on output, (n, nrhs)
!> @param[in]     ipiv   pivot indices from wrap_getrf, (n)
!> @param[out]    info   zero on success
!> @param[in]     trans  "n" (default), "t" (plain transpose) or "c" (conjugate)
   subroutine wrap_zgetrs(amat, bmat, ipiv, info, trans)
      !> Factors
      complex(dp), intent(in) :: amat(:, :)
      !> Right-hand sides and solutions
      complex(dp), intent(inout) :: bmat(:, :)
      !> Pivot indices
      integer, intent(in) :: ipiv(:)
      !> Status
      integer, intent(out) :: info
      !> Orientation of the system
      character(len=1), intent(in), optional :: trans
      character(len=1) :: tra
      integer(lapack_ik) :: n, nrhs, lda, ldb, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)
      if (present(trans)) then
         tra = trans
      else
         tra = "n"
      end if
      lda = int(max(1, size(amat, 1)), lapack_ik)
      ldb = int(max(1, size(bmat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      nrhs = int(size(bmat, 2), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      ipiv_lapack = int(ipiv, lapack_ik)
      call lapack_getrs(tra, n, nrhs, amat, lda, ipiv_lapack, bmat, ldb, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_zgetrs

end module moist_math_lapack_getrs
