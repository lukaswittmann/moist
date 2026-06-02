!> Provides interface to LAPACK singular value decomposition solver (?gesvd)

module moist_math_lapack_gesvd
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none
   private

   public :: lapack_gesvd
   public :: sgesvd, dgesvd

   interface lapack_gesvd
      pure subroutine sgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
         import :: sp, lapack_ik
         character(len=1), intent(in) :: jobu
         character(len=1), intent(in) :: jobvt
         integer(lapack_ik), intent(in) :: m
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         integer(lapack_ik), intent(in) :: ldu
         integer(lapack_ik), intent(in) :: ldvt
         integer(lapack_ik), intent(in) :: lwork
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: s(*)
         real(sp), intent(out) :: u(ldu, *)
         real(sp), intent(out) :: vt(ldvt, *)
         real(sp), intent(inout) :: work(*)
         integer(lapack_ik), intent(out) :: info
      end subroutine sgesvd
      pure subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
         import :: dp, lapack_ik
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
