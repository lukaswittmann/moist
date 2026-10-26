!> Provides interface to LAPACK least-squares solver (?gels)

module moist_math_lapack_gels
   use mctc_env, only: dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none(type, external)
   private

   public :: lapack_gels
   public :: dgels

   interface lapack_gels
      pure subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
         import :: dp, lapack_ik
         implicit none(type, external)
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
