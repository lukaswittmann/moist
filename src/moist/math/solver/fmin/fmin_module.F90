!*****************************************************************************************
!>
!  Module for [[fmin]] 1D derative-free function minimizer.
!
!### License
!  * [BSD-3](https://github.com/jacobwilliams/fmin/blob/master/LICENSE)
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:
#ifdef REAL32
!      `real(kind=real32)` [4 bytes]
#elif REAL64
!      `real(kind=real64)` [8 bytes]
#elif REAL128
!      `real(kind=real128)` [16 bytes]
#else
!      `real(kind=real64)` [8 bytes]
#endif

module moist_math_solver_fmin

   use mctc_env_accuracy, only: wp

   implicit none

   private

   abstract interface
      function func(x) result(f)
        !! interface for user function
         import :: wp
         implicit none
         real(wp), intent(in) :: x  !! indep. variable
         real(wp)            :: f  !! function value `f(x)`
      end function func
   end interface

   public :: fmin

contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  An approximation x to the point where `f` attains a minimum on
!  the interval `(ax,bx)` is determined.
!
!  the method used is a combination of golden section search and
!  successive parabolic interpolation. convergence is never much slower
!  than that for a fibonacci search. if `f` has a continuous second
!  derivative which is positive at the minimum (which is not at `ax` or
!  `bx`), then convergence is superlinear, and usually of the order of
!  about 1.324.
!
!  the function `f` is never evaluated at two points closer together
!  than `eps*abs(fmin) + (tol/3)`, where `eps` is approximately the square
!  root of the relative machine precision. if `f` is a unimodal
!  function and the computed values of `f` are always unimodal when
!  separated by at least `eps*abs(x) + (tol/3)`, then fmin approximates
!  the abcissa of the global minimum of `f` on the interval `ax,bx` with
!  an error less than `3*eps*abs(fmin) + tol`. if `f` is not unimodal,
!  then `fmin` may approximate a local, but perhaps non-global, minimum to
!  the same accuracy.
!
!### Reference
!  * Richard brent, "algorithms for minimization without derivatives",
!    prentice - hall, inc. (1973).
!
!### See also
!  * [fmin from Netlib](http://www.netlib.org/fmm/fmin.f)

   function fmin(f, ax, bx, tol) result(xmin)

      implicit none

      procedure(func)     :: f    !! the function to minimize
      real(wp), intent(in) :: ax   !! left endpoint of initial interval
      real(wp), intent(in) :: bx   !! right endpoint of initial interval
      real(wp), intent(in) :: tol  !! desired length of the interval of
                                !! uncertainty of the final result (>=0)
      real(wp)            :: xmin !! abcissa approximating the point where
                                !! f attains a minimum

      real(wp) :: a, b, d, e, xm, p, q, r, tol1, tol2, u, v, w
      real(wp) :: fu, fv, fw, fx, x

      real(wp), parameter :: c = (3.0_wp - sqrt(5.0_wp))/2.0_wp  !! squared inverse of golden ratio
      real(wp), parameter :: half = 0.5_wp
      real(wp), parameter :: sqrteps = sqrt(epsilon(1.0_wp))

      ! initialization

      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0_wp
      fx = f(x)
      fv = fx
      fw = fx

      do    !  main loop starts here

         xm = half*(a + b)
         tol1 = sqrteps*abs(x) + tol/3.0_wp
         tol2 = 2.0_wp*tol1

         !  check stopping criterion

         if (abs(x - xm) <= (tol2 - half*(b - a))) then
            ! write(*,*) 'x             = ', x
            ! write(*,*) 'xm            = ', xm
            ! write(*,*) 'abs(x - xm)   = ', abs(x - xm)
            ! write(*,*) 'tol2          = ', tol2
            ! write(*,*) 'half*(b - a)  = ', half*(b - a)
            exit
         end if

         ! is golden-section necessary

         if (abs(e) <= tol1) then

            !  a golden-section step

            if (x >= xm) then
               e = a - x
            else
               e = b - x
            end if
            d = c*e

         else

            !  fit parabola

            r = (x - w)*(fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v)*q - (x - w)*r
            q = 2.0_wp*(q - r)
            if (q > 0.0_wp) p = -p
            q = abs(q)
            r = e
            e = d

            !  is parabola acceptable

            if ((abs(p) >= abs(half*q*r)) .or. (p <= q*(a - x)) .or. (p >= q*(b - x))) then

               !  a golden-section step

               if (x >= xm) then
                  e = a - x
               else
                  e = b - x
               end if
               d = c*e

            else

               !  a parabolic interpolation step

               d = p/q
               u = x + d

               !  f must not be evaluated too close to ax or bx

               if (((u - a) < tol2) .or. ((b - u) < tol2)) d = sign(tol1, xm - x)

            end if

         end if

         !  f must not be evaluated too close to x

         if (abs(d) >= tol1) then
            u = x + d
         else
            u = x + sign(tol1, d)
         end if
         fu = f(u)

         !  update a, b, v, w, and x

         if (fu <= fx) then
            if (u >= x) a = x
            if (u < x) b = x
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
         else
            if (u < x) a = u
            if (u >= x) b = u
            if (fu <= fw .or. w == x) then
               v = w
               fv = fw
               w = u
               fw = fu
            else if (fu <= fv .or. v == x .or. v == w) then
               v = u
               fv = fu
            end if
         end if

      end do    !  end of main loop

      xmin = x

   end function fmin
!*****************************************************************************************

!*****************************************************************************************
end module moist_math_solver_fmin
!*****************************************************************************************