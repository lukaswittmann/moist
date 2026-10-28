!> Provides wrappers to solve a symmetric indefinite linear equation system

!> Wrappers to solve a system of linear equations with Bunch-Kaufman factors
module moist_math_lapack_sytrs
   use mctc_env, only: dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none(type, external)
   private

   public :: wrap_sytrs

   !> Solves A * X = B for a real or complex symmetric matrix A with factors from ?SYTRF
   !>
   !> The factorization is A = U**T * D * U or A = L * D * L**T; complex
   !> matrices are symmetric (A = A**T), not Hermitian.
   interface wrap_sytrs
      module procedure :: wrap_dsytrs
      module procedure :: wrap_zsytrs
   end interface wrap_sytrs

   !> LAPACK solve with Bunch-Kaufman factors of a real or complex symmetric matrix
   interface lapack_sytrs
      !> Double precision solve with symmetric indefinite factors
      !>
      !> @param[in]     uplo  "U" or "L", as used in dsytrf
      !> @param[in]     n     order of a
      !> @param[in]     nrhs  number of right-hand sides
      !> @param[in]     a     factors from dsytrf
      !> @param[in]     lda   leading dimension of a
      !> @param[in]     ipiv  pivot indices from dsytrf
      !> @param[in,out] b     right-hand sides on input, solutions on output
      !> @param[in]     ldb   leading dimension of b
      !> @param[out]    info  zero on success
      pure subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         !> Triangle of the factorization
         character(len=1), intent(in) :: uplo
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
      end subroutine dsytrs
      !> Double complex solve with symmetric indefinite factors
      !>
      !> @param[in]     uplo  "U" or "L", as used in dsytrf
      !> @param[in]     n     order of a
      !> @param[in]     nrhs  number of right-hand sides
      !> @param[in]     a     factors from dsytrf
      !> @param[in]     lda   leading dimension of a
      !> @param[in]     ipiv  pivot indices from dsytrf
      !> @param[in,out] b     right-hand sides on input, solutions on output
      !> @param[in]     ldb   leading dimension of b
      !> @param[out]    info  zero on success
      pure subroutine zsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         !> Triangle of the factorization
         character(len=1), intent(in) :: uplo
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
      end subroutine zsytrs
   end interface lapack_sytrs

contains

!> Double precision solve with factors from wrap_sytrf
!>
!> @param[in]     amat  factors from wrap_sytrf, (n, n)
!> @param[in,out] bmat  right-hand sides on input, solutions on output, (n, nrhs)
!> @param[in]     ipiv  pivot indices from wrap_sytrf, (n)
!> @param[out]    info  zero on success
!> @param[in]     uplo  triangle used in wrap_sytrf, "u" (default) or "l"
   subroutine wrap_dsytrs(amat, bmat, ipiv, info, uplo)
      !> Factors
      real(dp), intent(in) :: amat(:, :)
      !> Right-hand sides and solutions
      real(dp), intent(inout) :: bmat(:, :)
      !> Pivot indices
      integer, intent(in) :: ipiv(:)
      !> Status
      integer, intent(out) :: info
      !> Triangle of the factorization
      character(len=1), intent(in), optional :: uplo
      character(len=1) :: ula
      integer(lapack_ik) :: n, nrhs, lda, ldb, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)

      ula = "u"
      if (present(uplo)) ula = uplo
      lda = int(max(1, size(amat, 1)), lapack_ik)
      ldb = int(max(1, size(bmat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      nrhs = int(size(bmat, 2), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      ipiv_lapack = int(ipiv, lapack_ik)
      call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv_lapack, bmat, ldb, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_dsytrs

!> Double complex solve with factors from wrap_sytrf
!>
!> @param[in]     amat  factors from wrap_sytrf, (n, n)
!> @param[in,out] bmat  right-hand sides on input, solutions on output, (n, nrhs)
!> @param[in]     ipiv  pivot indices from wrap_sytrf, (n)
!> @param[out]    info  zero on success
!> @param[in]     uplo  triangle used in wrap_sytrf, "u" (default) or "l"
   subroutine wrap_zsytrs(amat, bmat, ipiv, info, uplo)
      !> Factors
      complex(dp), intent(in) :: amat(:, :)
      !> Right-hand sides and solutions
      complex(dp), intent(inout) :: bmat(:, :)
      !> Pivot indices
      integer, intent(in) :: ipiv(:)
      !> Status
      integer, intent(out) :: info
      !> Triangle of the factorization
      character(len=1), intent(in), optional :: uplo
      character(len=1) :: ula
      integer(lapack_ik) :: n, nrhs, lda, ldb, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)

      ula = "u"
      if (present(uplo)) ula = uplo
      lda = int(max(1, size(amat, 1)), lapack_ik)
      ldb = int(max(1, size(bmat, 1)), lapack_ik)
      n = int(size(amat, 2), lapack_ik)
      nrhs = int(size(bmat, 2), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      ipiv_lapack = int(ipiv, lapack_ik)
      call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv_lapack, bmat, ldb, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_zsytrs

end module moist_math_lapack_sytrs
