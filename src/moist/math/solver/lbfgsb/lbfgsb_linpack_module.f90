
!*******************************************************************************
!> license: BSD
!
!  LINPACK support routines for LBFGSB.
!  These have been refactored into modern Fortran.

module lbfgsb_linpack_module

   use moist_math_blas_legacy, only: daxpy, ddot
   use mctc_env_accuracy, only: wp

   implicit none

   private

   public :: dpofa, dtrsl

contains
!*******************************************************************************

!*******************************************************************************
!>
!  dpofa factors a real symmetric positive definite matrix.
!
!### History
!  * linpack.  this version dated 08/14/78 .
!    cleve moler, university of new mexico, argonne national lab.

   subroutine dpofa(a, lda, n, info)

      integer, intent(in) :: lda !! the leading dimension of the array `a`.
      integer, intent(in) :: n !! the order of the matrix a.
      integer, intent(out) :: info !!  * `info = 0` for normal return.
                               !!  * `info = k` signals an error condition.  the leading minor
                               !!     of order `k` is not positive definite.
      real(wp), intent(inout) :: a(lda, *) !! Dimension `(lda, n)`:
                                      !!
                                      !!  * On entry: the symmetric matrix to be factored.  only the
                                      !!    diagonal and upper triangle are used.
                                      !!  * On return: an upper triangular matrix `r` so that `a = trans(r)*r`
                                      !!    where `trans(r)` is the transpose.
                                      !!    the strict lower triangle is unaltered.
                                      !!    if `info /= 0`, the factorization is not complete.

      real(wp) :: t
      real(wp) :: s
      integer :: j, jm1, k

      do j = 1, n
         info = j
         s = 0.0_wp
         jm1 = j - 1
         if (jm1 >= 1) then
            do k = 1, jm1
               t = a(k, j) - ddot(k - 1, a(1, k), 1, a(1, j), 1)
               t = t/a(k, k)
               a(k, j) = t
               s = s + t*t
            end do
         end if
         s = a(j, j) - s
         if (s <= 0.0_wp) return
         a(j, j) = sqrt(s)
      end do
      info = 0
   end subroutine dpofa
!*******************************************************************************

!*******************************************************************************
!>
!  dtrsl solves systems of the form
!
!   `t * x = b`
!
!  or
!
!   `trans(t) * x = b`
!
!  where t is a triangular matrix of order n. here trans(t)
!  denotes the transpose of the matrix t.
!
!### History
!  * linpack. this version dated 08/14/78 .
!    g. w. stewart, university of maryland, argonne national lab.

   subroutine dtrsl(t, ldt, n, b, job, info)

      integer, intent(in) :: ldt !! the leading dimension of the array t.
      integer, intent(in) :: n !! the order of the system.
      integer, intent(in) :: job !! job specifies what kind of system is to be solved.
                          !! if job is:
                          !!
                          !!  * `00`   solve `t*x=b`, t lower triangular,
                          !!  * `01`   solve `t*x=b`, t upper triangular,
                          !!  * `10`   solve `trans(t)*x=b`, t lower triangular,
                          !!  * `11`   solve `trans(t)*x=b`, t upper triangular.
      integer, intent(out) :: info !! info contains zero if the system is nonsingular.
                            !! otherwise info contains the index of
                            !! the first zero diagonal element of t.
      real(wp), intent(in) :: t(ldt, *) !! t contains the matrix of the system. the zero
                                !! elements of the matrix are not referenced, and
                                !! the corresponding elements of the array can be
                                !! used to store other information.
      real(wp), intent(inout) :: b(*) !! On entry: the right hand side of the system.
                               !! On return: the solution, if info == 0. otherwise b is unaltered.

      real(wp) :: temp
      integer :: case, j, jj

      ! check for zero diagonal elements.
      do info = 1, n
         if (t(info, info) == 0.0_wp) return
      end do
      info = 0

      ! determine the task and go to it.
      case = 1
      if (mod(job, 10) /= 0) case = 2
      if (mod(job, 100)/10 /= 0) case = case + 2

      select case (case)
      case (1) ! solve t*x=b for t lower triangular
         b(1) = b(1)/t(1, 1)
         if (n >= 2) then
            do j = 2, n
               temp = -b(j - 1)
               call daxpy(n - j + 1, temp, t(j, j - 1), 1, b(j), 1)
               b(j) = b(j)/t(j, j)
            end do
         end if

      case (2) ! solve t*x=b for t upper triangular.
         b(n) = b(n)/t(n, n)
         if (n >= 2) then
            do jj = 2, n
               j = n - jj + 1
               temp = -b(j + 1)
               call daxpy(j, temp, t(1, j + 1), 1, b(1), 1)
               b(j) = b(j)/t(j, j)
            end do
         end if

      case (3) ! solve trans(t)*x=b for t lower triangular.
         b(n) = b(n)/t(n, n)
         if (n >= 2) then
            do jj = 2, n
               j = n - jj + 1
               b(j) = b(j) - ddot(jj - 1, t(j + 1, j), 1, b(j + 1), 1)
               b(j) = b(j)/t(j, j)
            end do
         end if

      case (4) ! solve trans(t)*x=b for t upper triangular.
         b(1) = b(1)/t(1, 1)
         if (n >= 2) then
            do j = 2, n
               b(j) = b(j) - ddot(j - 1, t(1, j), 1, b(1), 1)
               b(j) = b(j)/t(j, j)
            end do
         end if

      end select

   end subroutine dtrsl

end module lbfgsb_linpack_module
