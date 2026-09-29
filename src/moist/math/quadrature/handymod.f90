!> HandyMod finite radial quadratures for molecular integration grids
module moist_math_quadrature_handymod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   implicit none(type, external)
   private

   public :: handymod_radii
   public :: handymod_midpoint_radii

contains

   !> Generate Chebyshev-II HandyMod radial nodes and weights on a finite interval
   !>
   !> @param[in]  nr      Number of radial points (nr >= 1)
   !> @param[in]  rmin    Minimum radial distance in bohr
   !> @param[in]  rmax    Maximum radial distance in bohr (must exceed rmin)
   !> @param[in]  m       HandyMod map exponent (m >= 1 for this legacy path)
   !> @param[out] radii   Radial nodes in bohr, shape at least (nr)
   !> @param[out] weights Radial weights in bohr^3, including r^2 dr/dx
   !> @param[out] error   Propagated validation error
   subroutine handymod_radii(nr, rmin, rmax, m, radii, weights, error)
      !> Number of radial points
      integer, intent(in)  :: nr
      !> HandyMod map exponent
      integer, intent(in)  :: m
      !> Minimum radial distance in bohr
      real(wp), intent(in)  :: rmin
      !> Maximum radial distance in bohr
      real(wp), intent(in)  :: rmax
      !> Radial nodes in bohr
      real(wp), intent(out) :: radii(:)
      !> Radial weights in bohr^3 (carrying the r^2 dr Jacobian)
      real(wp), intent(out) :: weights(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: ir
      real(wp) :: sx, theta, x

      if (.not. validate_handymod_inputs(nr, rmin, rmax, real(m, wp), radii, weights, &
            & "HandyMod quadrature", error)) return
      if (m < 1) then
         call fatal_error(error, "HandyMod quadrature: m must be >= 1")
         return
      end if

      do ir = 1, nr
         theta = real(ir, wp)*pi/real(nr + 1, wp)
         x = cos(theta)
         sx = sin(theta)
         call transform_handymod_node(x, sx*pi/real(nr + 1, wp), rmin, rmax, real(m, wp), &
            & radii(ir), weights(ir), error)
         if (allocated(error)) return
      end do
   end subroutine handymod_radii

   !> Generate midpoint HandyMod radial nodes and weights
   !>
   !> @param[in]  nr      Number of midpoint nodes on [-1, 1] (nr >= 1)
   !> @param[in]  rmin    Minimum radial distance in bohr
   !> @param[in]  rmax    Maximum radial distance in bohr (must exceed rmin)
   !> @param[in]  m       Real-valued HandyMod transform parameter (m > 0)
   !> @param[out] radii   Radial nodes in bohr, shape at least (nr)
   !> @param[out] weights Radial weights in bohr^3, including r^2 dr/dx
   !> @param[out] error   Propagated validation error
   subroutine handymod_midpoint_radii(nr, rmin, rmax, m, radii, weights, error)
      !> Number of radial midpoint nodes
      integer, intent(in)  :: nr
      !> Minimum radial distance in bohr
      real(wp), intent(in)  :: rmin
      !> Maximum radial distance in bohr
      real(wp), intent(in)  :: rmax
      !> HandyMod transform parameter
      real(wp), intent(in)  :: m
      !> Radial nodes in bohr
      real(wp), intent(out) :: radii(:)
      !> Radial weights in bohr^3 (carrying the r^2 dr Jacobian)
      real(wp), intent(out) :: weights(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: ir
      real(wp) :: wx, x

      if (.not. validate_handymod_inputs(nr, rmin, rmax, m, radii, weights, &
            & "HandyMod midpoint quadrature", error)) return

      wx = 2.0_wp/real(nr, wp)
      do ir = 1, nr
         x = -1.0_wp + (2.0_wp*real(ir, wp) - 1.0_wp)/real(nr, wp)
         call transform_handymod_node(x, wx, rmin, rmax, m, radii(ir), weights(ir), error)
         if (allocated(error)) return
      end do
   end subroutine handymod_midpoint_radii

   !> Validate common HandyMod transform inputs
   !>
   !> @param[in]  nr      Number of radial nodes (nr >= 1)
   !> @param[in]  rmin    Minimum radial distance in bohr
   !> @param[in]  rmax    Maximum radial distance in bohr (must exceed rmin)
   !> @param[in]  m       HandyMod transform parameter (m > 0)
   !> @param[out] radii   Radial node output array, zeroed
   !> @param[out] weights Radial weight output array, zeroed
   !> @param[in]  label   Error-message label
   !> @param[out] error   Propagated validation error
   function validate_handymod_inputs(nr, rmin, rmax, m, radii, weights, label, error) result(valid)
      !> Number of radial nodes
      integer, intent(in) :: nr
      !> Minimum radial distance in bohr
      real(wp), intent(in) :: rmin
      !> Maximum radial distance in bohr
      real(wp), intent(in) :: rmax
      !> HandyMod transform parameter
      real(wp), intent(in) :: m
      !> Radial node output array
      real(wp), intent(out) :: radii(:)
      !> Radial weight output array
      real(wp), intent(out) :: weights(:)
      !> Error-message label
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Whether inputs are valid
      logical :: valid

      real(wp) :: c, d

      valid = .false.
      radii = 0.0_wp
      weights = 0.0_wp

      if (nr < 1) then
         call fatal_error(error, label//": nr must be >= 1")
         return
      end if
      if (.not. ieee_is_finite(rmin) .or. .not. ieee_is_finite(rmax)) then
         call fatal_error(error, label//": radial bounds must be finite")
         return
      end if
      if (.not. ieee_is_finite(m)) then
         call fatal_error(error, label//": m must be finite")
         return
      end if
      if (m <= 0.0_wp) then
         call fatal_error(error, label//": m must be > 0")
         return
      end if
      if (rmax <= rmin) then
         call fatal_error(error, label//": rmax must be greater than rmin")
         return
      end if
      if (size(radii) < nr) then
         call fatal_error(error, label//": radii array is too small")
         return
      end if
      if (size(weights) < nr) then
         call fatal_error(error, label//": weights array is too small")
         return
      end if

      d = rmax - rmin
      c = 2.0_wp**m
      if (d <= c - 1.0_wp) then
         call fatal_error(error, &
            & label//": rmax - rmin must be greater than 2^m - 1")
         return
      end if

      valid = .true.
   end function validate_handymod_inputs

   !> Apply the HandyMod radial transform to one node and quadrature weight
   !>
   !> @param[in]  x       Base-rule node on [-1, 1]
   !> @param[in]  wx      Base-rule weight on [-1, 1]
   !> @param[in]  rmin    Minimum radial distance in bohr
   !> @param[in]  rmax    Maximum radial distance in bohr (must exceed rmin)
   !> @param[in]  m       HandyMod transform parameter (m > 0)
   !> @param[out] radius  Transformed radial node in bohr
   !> @param[out] weight  Transformed radial weight in bohr^3
   !> @param[out] error   Propagated validation error
   subroutine transform_handymod_node(x, wx, rmin, rmax, m, radius, weight, error)
      !> Base-rule node on [-1, 1]
      real(wp), intent(in) :: x
      !> Base-rule weight on [-1, 1]
      real(wp), intent(in) :: wx
      !> Minimum radial distance in bohr
      real(wp), intent(in) :: rmin
      !> Maximum radial distance in bohr
      real(wp), intent(in) :: rmax
      !> HandyMod transform parameter
      real(wp), intent(in) :: m
      !> Transformed radial node in bohr
      real(wp), intent(out) :: radius
      !> Transformed radial weight in bohr^3
      real(wp), intent(out) :: weight
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: a, c, d, den, drdx, tol, u

      radius = 0.0_wp
      weight = 0.0_wp
      d = rmax - rmin
      c = 2.0_wp**m
      a = c*(1.0_wp - c + d)
      if (a <= 0.0_wp .or. .not. ieee_is_finite(a)) then
         call fatal_error(error, &
            & "HandyMod quadrature: invalid denominator scale")
         return
      end if

      u = (1.0_wp + x)**m
      den = a - u*(d - c)
      if (den <= 0.0_wp .or. .not. ieee_is_finite(den)) then
         call fatal_error(error, "HandyMod quadrature: invalid denominator")
         return
      end if

      radius = rmin + d*u/den
      drdx = d*m*(1.0_wp + x)**(m - 1.0_wp)*a/den**2
      weight = wx*radius**2*drdx

      tol = 100.0_wp*epsilon(1.0_wp)*max(1.0_wp, abs(rmin), abs(rmax))
      if (.not. ieee_is_finite(radius) .or. .not. ieee_is_finite(weight)) then
         call fatal_error(error, &
            & "HandyMod quadrature: generated non-finite radius or weight")
         return
      end if
      if (weight <= 0.0_wp) then
         call fatal_error(error, &
            & "HandyMod quadrature: generated non-positive weight")
         return
      end if
      if (radius < rmin - tol .or. radius > rmax + tol) then
         call fatal_error(error, &
            & "HandyMod quadrature: generated radius outside requested interval")
         return
      end if
   end subroutine transform_handymod_node

end module moist_math_quadrature_handymod
