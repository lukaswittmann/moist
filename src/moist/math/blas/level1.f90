!> @file moist/blas/level1.f90
!> Provides interfactes to level 1 BLAS routines

!> High-level interface to level 1 basic linear algebra subprogram operations
module moist_math_blas_level1
   use mctc_env, only: sp, dp
   implicit none
   private

   public :: wrap_dot
   public :: wrap_nrm2
   public :: wrap_scal
   public :: wrap_copy
   public :: wrap_axpy

   !> Forms the dot product of two vectors.
   interface wrap_dot
      module procedure :: wrap_sdot
      module procedure :: wrap_ddot
      module procedure :: wrap_sdot12
      module procedure :: wrap_sdot21
      module procedure :: wrap_sdot22
      module procedure :: wrap_ddot12
      module procedure :: wrap_ddot21
      module procedure :: wrap_ddot22
   end interface wrap_dot

   !> Euclidean (2-)norm of a vector.
   interface wrap_nrm2
      module procedure :: wrap_snrm2
      module procedure :: wrap_dnrm2
   end interface wrap_nrm2

   !> Scale a vector in place by a scalar.
   interface wrap_scal
      module procedure :: wrap_sscal
      module procedure :: wrap_dscal
   end interface wrap_scal

   !> Copy a vector into another.
   interface wrap_copy
      module procedure :: wrap_scopy
      module procedure :: wrap_dcopy
   end interface wrap_copy

   !> Constant times a vector plus a vector (y := y + a*x).
   interface wrap_axpy
      module procedure :: wrap_saxpy
      module procedure :: wrap_daxpy
   end interface wrap_axpy

   !> Forms the dot product of two vectors.
   !> Uses unrolled loops for increments equal to one.
   interface blas_dot
      pure function sdot(n, x, incx, y, incy)
         import :: sp
         real(sp) :: sdot
         real(sp), intent(in) :: x(*)
         real(sp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function sdot
      pure function ddot(n, x, incx, y, incy)
         import :: dp
         real(dp) :: ddot
         real(dp), intent(in) :: x(*)
         real(dp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function ddot
   end interface blas_dot

   !> Euclidean norm of a vector (level-1 BLAS).
   interface blas_nrm2
      pure function snrm2(n, x, incx)
         import :: sp
         real(sp) :: snrm2
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function snrm2
      pure function dnrm2(n, x, incx)
         import :: dp
         real(dp) :: dnrm2
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function dnrm2
   end interface blas_nrm2

   !> Scales a vector by a constant (level-1 BLAS).
   interface blas_scal
      pure subroutine sscal(n, a, x, incx)
         import :: sp
         real(sp), intent(in) :: a
         real(sp), intent(inout) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine sscal
      pure subroutine dscal(n, a, x, incx)
         import :: dp
         real(dp), intent(in) :: a
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine dscal
   end interface blas_scal

   !> Copies a vector into another (level-1 BLAS).
   interface blas_copy
      pure subroutine scopy(n, x, incx, y, incy)
         import :: sp
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine scopy
      pure subroutine dcopy(n, x, incx, y, incy)
         import :: dp
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine dcopy
   end interface blas_copy

   !> Constant times a vector plus a vector (level-1 BLAS).
   interface blas_axpy
      pure subroutine saxpy(n, a, x, incx, y, incy)
         import :: sp
         real(sp), intent(in) :: a
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine saxpy
      pure subroutine daxpy(n, a, x, incx, y, incy)
         import :: dp
         real(dp), intent(in) :: a
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine daxpy
   end interface blas_axpy

contains

   function wrap_sdot(xvec, yvec) result(dot)
      real(sp) :: dot
      real(sp), intent(in) :: xvec(:)
      real(sp), intent(in) :: yvec(:)
      integer :: incx, incy, n
      incx = 1
      incy = 1
      n = size(xvec)
      dot = blas_dot(n, xvec, incx, yvec, incy)
   end function wrap_sdot

   function wrap_ddot(xvec, yvec) result(dot)
      real(dp) :: dot
      real(dp), intent(in) :: xvec(:)
      real(dp), intent(in) :: yvec(:)
      integer :: incx, incy, n
      incx = 1
      incy = 1
      n = size(xvec)
      dot = blas_dot(n, xvec, incx, yvec, incy)
   end function wrap_ddot

   function wrap_sdot12(xvec, yvec) result(dot)
      real(sp) :: dot
      real(sp), intent(in) :: xvec(:)
      real(sp), intent(in), contiguous, target :: yvec(:, :)
      real(sp), pointer :: yptr(:)
      yptr(1:size(yvec)) => yvec
      dot = wrap_dot(xvec, yptr)
   end function wrap_sdot12

   function wrap_sdot21(xvec, yvec) result(dot)
      real(sp) :: dot
      real(sp), intent(in), contiguous, target :: xvec(:, :)
      real(sp), intent(in) :: yvec(:)
      real(sp), pointer :: xptr(:)
      xptr(1:size(xvec)) => xvec
      dot = wrap_dot(xptr, yvec)
   end function wrap_sdot21

   function wrap_sdot22(xvec, yvec) result(dot)
      real(sp) :: dot
      real(sp), intent(in), contiguous, target :: xvec(:, :)
      real(sp), intent(in), contiguous, target :: yvec(:, :)
      real(sp), pointer :: xptr(:), yptr(:)
      xptr(1:size(xvec)) => xvec
      yptr(1:size(yvec)) => yvec
      dot = wrap_dot(xptr, yptr)
   end function wrap_sdot22

   function wrap_ddot12(xvec, yvec) result(dot)
      real(dp) :: dot
      real(dp), intent(in) :: xvec(:)
      real(dp), intent(in), contiguous, target :: yvec(:, :)
      real(dp), pointer :: yptr(:)
      yptr(1:size(yvec)) => yvec
      dot = wrap_dot(xvec, yptr)
   end function wrap_ddot12

   function wrap_ddot21(xvec, yvec) result(dot)
      real(dp) :: dot
      real(dp), intent(in), contiguous, target :: xvec(:, :)
      real(dp), intent(in) :: yvec(:)
      real(dp), pointer :: xptr(:)
      xptr(1:size(xvec)) => xvec
      dot = wrap_dot(xptr, yvec)
   end function wrap_ddot21

   function wrap_ddot22(xvec, yvec) result(dot)
      real(dp) :: dot
      real(dp), intent(in), contiguous, target :: xvec(:, :)
      real(dp), intent(in), contiguous, target :: yvec(:, :)
      real(dp), pointer :: xptr(:), yptr(:)
      xptr(1:size(xvec)) => xvec
      yptr(1:size(yvec)) => yvec
      dot = wrap_dot(xptr, yptr)
   end function wrap_ddot22

   function wrap_snrm2(xvec) result(nrm)
      real(sp) :: nrm
      real(sp), intent(in) :: xvec(:)
      nrm = blas_nrm2(size(xvec), xvec, 1)
   end function wrap_snrm2

   function wrap_dnrm2(xvec) result(nrm)
      real(dp) :: nrm
      real(dp), intent(in) :: xvec(:)
      nrm = blas_nrm2(size(xvec), xvec, 1)
   end function wrap_dnrm2

   subroutine wrap_sscal(xvec, a)
      real(sp), intent(inout) :: xvec(:)
      real(sp), intent(in) :: a
      call blas_scal(size(xvec), a, xvec, 1)
   end subroutine wrap_sscal

   subroutine wrap_dscal(xvec, a)
      real(dp), intent(inout) :: xvec(:)
      real(dp), intent(in) :: a
      call blas_scal(size(xvec), a, xvec, 1)
   end subroutine wrap_dscal

   subroutine wrap_scopy(xvec, yvec)
      real(sp), intent(in) :: xvec(:)
      real(sp), intent(inout) :: yvec(:)
      call blas_copy(size(xvec), xvec, 1, yvec, 1)
   end subroutine wrap_scopy

   subroutine wrap_dcopy(xvec, yvec)
      real(dp), intent(in) :: xvec(:)
      real(dp), intent(inout) :: yvec(:)
      call blas_copy(size(xvec), xvec, 1, yvec, 1)
   end subroutine wrap_dcopy

   subroutine wrap_saxpy(a, xvec, yvec)
      real(sp), intent(in) :: a
      real(sp), intent(in) :: xvec(:)
      real(sp), intent(inout) :: yvec(:)
      call blas_axpy(size(xvec), a, xvec, 1, yvec, 1)
   end subroutine wrap_saxpy

   subroutine wrap_daxpy(a, xvec, yvec)
      real(dp), intent(in) :: a
      real(dp), intent(in) :: xvec(:)
      real(dp), intent(inout) :: yvec(:)
      call blas_axpy(size(xvec), a, xvec, 1, yvec, 1)
   end subroutine wrap_daxpy

end module moist_math_blas_level1
