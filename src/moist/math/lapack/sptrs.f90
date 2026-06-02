!> Provides wrappers for solving symmetric packed linear systems.
module moist_math_lapack_sptrs
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: wrap_dsptrf, wrap_dsptrs

   interface lapack_sptrf
      pure subroutine dsptrf(uplo, n, ap, ipiv, info)
         import :: dp, lapack_ik
         character(len=1), intent(in) :: uplo
         integer(lapack_ik), intent(in) :: n
         real(dp), intent(inout) :: ap(*)
         integer(lapack_ik), intent(out) :: ipiv(*)
         integer(lapack_ik), intent(out) :: info
      end subroutine dsptrf
   end interface lapack_sptrf

   interface lapack_sptrs
      pure subroutine dsptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
         import :: dp, lapack_ik
         character(len=1), intent(in) :: uplo
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: nrhs
         real(dp), intent(in) :: ap(*)
         integer(lapack_ik), intent(in) :: ipiv(*)
         integer(lapack_ik), intent(in) :: ldb
         real(dp), intent(inout) :: b(ldb, *)
         integer(lapack_ik), intent(out) :: info
      end subroutine dsptrs
   end interface lapack_sptrs

contains

   subroutine wrap_dsptrf(ap, ipiv, info, uplo)
      real(dp), intent(inout) :: ap(:)
      integer, intent(out) :: ipiv(:)
      integer, intent(out) :: info
      character(len=1), intent(in), optional :: uplo
      character(len=1) :: ula
      integer(lapack_ik) :: n, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)

      if (present(uplo)) then
         ula = uplo
      else
         ula = 'u'
      end if

      n = int((sqrt(real(8*size(ap) + 1, dp)) - 1.0_dp)/2.0_dp, lapack_ik)
      allocate (ipiv_lapack(max(1, size(ipiv))))
      call lapack_sptrf(ula, n, ap, ipiv_lapack, info_lapack)
      ipiv = int(ipiv_lapack, kind=kind(ipiv))
      info = int(info_lapack)
   end subroutine wrap_dsptrf

   subroutine wrap_dsptrs(ap, bmat, ipiv, info, uplo)
      real(dp), intent(in) :: ap(:)
      real(dp), intent(inout) :: bmat(:, :)
      integer, intent(in) :: ipiv(:)
      integer, intent(out) :: info
      character(len=1), intent(in), optional :: uplo
      character(len=1) :: ula
      integer(lapack_ik) :: n, nrhs, ldb, info_lapack
      integer(lapack_ik), allocatable :: ipiv_lapack(:)

      if (present(uplo)) then
         ula = uplo
      else
         ula = 'u'
      end if

      n = int((sqrt(real(8*size(ap) + 1, dp)) - 1.0_dp)/2.0_dp, lapack_ik)
      nrhs = int(size(bmat, 2), lapack_ik)
      ldb = int(max(1, size(bmat, 1)), lapack_ik)
      allocate (ipiv_lapack(size(ipiv)))
      ipiv_lapack = int(ipiv, lapack_ik)
      call lapack_sptrs(ula, n, nrhs, ap, ipiv_lapack, bmat, ldb, info_lapack)
      info = int(info_lapack)
   end subroutine wrap_dsptrs

end module moist_math_lapack_sptrs
