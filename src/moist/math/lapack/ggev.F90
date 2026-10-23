!> Provides interface to the generalized nonsymmetric eigenvalue driver

!> Interface to LAPACK generalized nonsymmetric eigenvalue driver
module moist_math_lapack_ggev
   use mctc_env, only: sp, dp
   use moist_math_lapack_kinds, only: lapack_ik
   implicit none(type, external)
   private

   public :: lapack_ggev
   public :: sggev, dggev

   !> Generalized eigenvalues and optional eigenvectors of a real nonsymmetric pair (A, B)
   !>
   !> - A and B are N-by-N
   !> - Eigenvalue lambda = alpha/beta, such that A - lambda*B is singular
   !> - Returned as the pair (alpha, beta), which also covers beta = 0 and alpha = beta = 0
   !> - Right eigenvector v(j) satisfies A*v(j) = lambda(j)*B*v(j)
   !> - Complex conjugate pairs occupy consecutive columns of VR:
   !>   v(j) = VR(:, j) + i*VR(:, j+1), v(j+1) = VR(:, j) - i*VR(:, j+1)
   interface lapack_ggev

      pure subroutine sggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, &
                            vl, ldvl, vr, ldvr, work, lwork, info)
         import :: sp, lapack_ik
         implicit none(type, external)
         character(len=1), intent(in) :: jobvl
         character(len=1), intent(in) :: jobvr
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(in) :: ldb
         real(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: alphar(*)
         real(sp), intent(out) :: alphai(*)
         real(sp), intent(out) :: beta(*)
         integer(lapack_ik), intent(in) :: ldvl
         real(sp), intent(out) :: vl(ldvl, *)
         integer(lapack_ik), intent(in) :: ldvr
         real(sp), intent(out) :: vr(ldvr, *)
         real(sp), intent(inout) :: work(*)
         integer(lapack_ik), intent(in) :: lwork
         integer(lapack_ik), intent(out) :: info
      end subroutine sggev

      pure subroutine dggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, &
                            vl, ldvl, vr, ldvr, work, lwork, info)
         import :: dp, lapack_ik
         implicit none(type, external)
         character(len=1), intent(in) :: jobvl
         character(len=1), intent(in) :: jobvr
         integer(lapack_ik), intent(in) :: n
         integer(lapack_ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         integer(lapack_ik), intent(in) :: ldb
         real(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(out) :: alphar(*)
         real(dp), intent(out) :: alphai(*)
         real(dp), intent(out) :: beta(*)
         integer(lapack_ik), intent(in) :: ldvl
         real(dp), intent(out) :: vl(ldvl, *)
         integer(lapack_ik), intent(in) :: ldvr
         real(dp), intent(out) :: vr(ldvr, *)
         real(dp), intent(inout) :: work(*)
         integer(lapack_ik), intent(in) :: lwork
         integer(lapack_ik), intent(out) :: info
      end subroutine dggev

   end interface lapack_ggev

end module moist_math_lapack_ggev
