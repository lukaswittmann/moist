!> Provides interface to symmetric eigenvalue solver

!> Interface to LAPACK symmetric eigenvalue solver
module moist_math_lapack_syev
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: lapack_syev
   public :: dsyev, ssyev

   !> Computes all eigenvalues and, optionally, eigenvectors of a
   !> real symmetric matrix A.
   interface lapack_syev
      pure subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import :: sp, lapack_ik
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         real(sp), intent(out) :: work(*)
         integer(lapack_ik), intent(in) :: lwork
         integer(lapack_ik), intent(out) :: info
      end subroutine ssyev
      pure subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import :: dp, lapack_ik
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         real(dp), intent(out) :: work(*)
         integer(lapack_ik), intent(in) :: lwork
         integer(lapack_ik), intent(out) :: info
      end subroutine dsyev
   end interface lapack_syev

end module moist_math_lapack_syev
