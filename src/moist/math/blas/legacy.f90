!> Classic-signature level-1 BLAS routines bound to the linked BLAS library.
!>
!> moist always links an external BLAS/LAPACK backend. Several vendored solvers
!> (slsqp, lbfgsb) call BLAS with the classic Fortran-77 convention (explicit
!> length and increments, array offsets, non-unit/zero strides) that the
!> whole-array, unit-stride wrappers in moist_math_blas cannot express. This
!> module exposes the double-precision routines they need with that classic
!> signature so the solvers use the optimized library instead of bundling their
!> own copies.
module moist_math_blas_legacy
   use mctc_env_accuracy, only: wp
   implicit none
   private

   public :: daxpy, dcopy, ddot, dnrm2, dscal

   interface
      !> Constant times a vector plus a vector: dy := dy + da*dx.
      pure subroutine daxpy(n, da, dx, incx, dy, incy)
         import :: wp
         implicit none
         !> Number of elements to process
         integer, intent(in) :: n
         !> Scalar multiplier
         real(wp), intent(in) :: da
         !> Source vector
         real(wp), dimension(*), intent(in) :: dx
         !> Stride of dx
         integer, intent(in) :: incx
         !> Destination vector, updated in place
         real(wp), dimension(*), intent(inout) :: dy
         !> Stride of dy
         integer, intent(in) :: incy
      end subroutine daxpy

      !> Copy a vector: dy := dx.
      pure subroutine dcopy(n, dx, incx, dy, incy)
         import :: wp
         implicit none
         !> Number of elements to process
         integer, intent(in) :: n
         !> Source vector
         real(wp), dimension(*), intent(in) :: dx
         !> Stride of dx
         integer, intent(in) :: incx
         !> Destination vector
         real(wp), dimension(*), intent(out) :: dy
         !> Stride of dy
         integer, intent(in) :: incy
      end subroutine dcopy

      !> Dot product of two vectors.
      pure real(wp) function ddot(n, dx, incx, dy, incy)
         import :: wp
         implicit none
         !> Number of elements to process
         integer, intent(in) :: n
         !> First vector
         real(wp), dimension(*), intent(in) :: dx
         !> Stride of dx
         integer, intent(in) :: incx
         !> Second vector
         real(wp), dimension(*), intent(in) :: dy
         !> Stride of dy
         integer, intent(in) :: incy
      end function ddot

      !> Euclidean (2-)norm of a vector.
      pure function dnrm2(n, x, incx) result(norm)
         import :: wp
         implicit none
         !> Number of elements to process
         integer, intent(in) :: n
         !> Input vector
         real(wp), dimension(*), intent(in) :: x
         !> Stride of x
         integer, intent(in) :: incx
         !> Resulting norm
         real(wp) :: norm
      end function dnrm2

      !> Scale a vector by a constant: dx := da*dx.
      pure subroutine dscal(n, da, dx, incx)
         import :: wp
         implicit none
         !> Number of elements to process
         integer, intent(in) :: n
         !> Scalar multiplier
         real(wp), intent(in) :: da
         !> Vector scaled in place
         real(wp), dimension(*), intent(inout) :: dx
         !> Stride of dx
         integer, intent(in) :: incx
      end subroutine dscal
   end interface

end module moist_math_blas_legacy
