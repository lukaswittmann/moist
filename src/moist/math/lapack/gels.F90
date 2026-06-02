!> Provides interface to LAPACK least-squares solver (?gels)

module moist_math_lapack_gels
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: lapack_gels
   public :: sgels, dgels

   interface lapack_gels
      pure subroutine sgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
         import :: sp, lapack_ik
         character(len=1), intent(in) :: trans
         integer(lapack_ik), intent(in) :: m
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: nrhs
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldb
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(inout) :: work(*)
         integer(lapack_ik), intent(in) :: lwork
         integer(lapack_ik), intent(out) :: info
      end subroutine sgels
      pure subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
         import :: dp, lapack_ik
         character(len=1), intent(in) :: trans
         integer(lapack_ik), intent(in) :: m
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: nrhs
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldb
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(inout) :: work(*)
         integer(lapack_ik), intent(in) :: lwork
         integer(lapack_ik), intent(out) :: info
      end subroutine dgels
   end interface lapack_gels

end module moist_math_lapack_gels
