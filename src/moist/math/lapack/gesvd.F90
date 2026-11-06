!> Provides interface to LAPACK singular value decomposition solver (?gesvd)

module moist_math_lapack_gesvd
   use mctc_env, only: dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none(type, external)
   private

   public :: lapack_gesvd
   public :: dgesvd

   interface lapack_gesvd
      pure subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         character(len=1), intent(in) :: jobu
         character(len=1), intent(in) :: jobvt
         integer(lapack_ik), intent(in) :: m
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldu
         integer(lapack_ik), intent(in) :: ldvt
         integer(lapack_ik), intent(in) :: lwork
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: s(*)
         real(dp), intent(out) :: u(ldu, *)
         real(dp), intent(out) :: vt(ldvt, *)
         real(dp), intent(inout) :: work(*)
         integer(lapack_ik), intent(out) :: info
      end subroutine dgesvd
   end interface lapack_gesvd

end module moist_math_lapack_gesvd
