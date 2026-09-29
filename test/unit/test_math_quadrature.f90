!> Test suite for the moist math/quadrature submodule
!>
!>   - Lebedev angular grid normalisation and unit-sphere placement
!>   - Chebyshev-2 radial quadrature of a Gaussian
!>   - HandyMod radial quadrature of finite ball/shell volumes
!>   - Midpoint-HandyMod radial quadrature
!>   - Becke partition-of-unity
!>
!> These exercise the bare quadrature rules in `moist_math_quadrature_*`,
!> independent of any grid abstraction (which is covered by `test_math_grid`)
module test_math_quadrature
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: center_at_origin
   use moist_math_quadrature_lebedev, only: get_angular_grid
   use moist_math_quadrature_chebyshev, only: chebyshev2_radii
   use moist_math_quadrature_handymod, only: handymod_radii, handymod_midpoint_radii
   use moist_math_quadrature_becke, only: becke_weights
   implicit none(type, external)
   private

   public :: collect_math_quadrature

contains

   !> Collect all math_quadrature tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_quadrature(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("angular_weights_sum_to_one", test_angular_weights_sum), &
         new_unittest("angular_points_on_unit_sphere", test_angular_unit_sphere), &
         new_unittest("radial_chebyshev_gauss", test_radial_chebyshev_gauss), &
         new_unittest("radial_handymod_ball_volume", test_radial_handymod_ball_volume), &
         new_unittest("radial_handymod_shell_volume", test_radial_handymod_shell_volume), &
         new_unittest("radial_handymod_midpoint_ball_volume", test_radial_handymod_midpoint_ball), &
         new_unittest("radial_handymod_midpoint_shell_volume", test_radial_handymod_midpoint_shell), &
         new_unittest("radial_handymod_midpoint_rejects_m", test_radial_handymod_midpoint_rejects_m), &
         new_unittest("becke_partition_of_unity", test_becke_partition_of_unity) &
         ]
   end subroutine collect_math_quadrature

   !> Lebedev weights should sum to 1 on the unit sphere
   !>
   !> @param[out] error  Test failure
   subroutine test_angular_weights_sum(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: npts_list(3) = [74, 230, 434]
      integer, parameter :: orders(3)    = [ 6,  12,  16]
      real(wp), allocatable :: xyz(:, :), w(:)
      type(mctc_error), allocatable :: merr
      integer :: k

      do k = 1, 3
         allocate(xyz(3, npts_list(k)), w(npts_list(k)))
         call get_angular_grid(orders(k), xyz, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, abs(sum(w) - 1.0_wp) < 1.0e-12_wp, &
            & "Lebedev weights do not sum to 1")
         deallocate(xyz, w)
         if (allocated(error)) return
      end do
   end subroutine test_angular_weights_sum

   !> Lebedev points should lie exactly on the unit sphere
   !>
   !> @param[out] error  Test failure
   subroutine test_angular_unit_sphere(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: npts_list(3) = [74, 230, 434]
      integer, parameter :: orders(3)    = [ 6,  12,  16]
      real(wp), allocatable :: xyz(:, :), w(:)
      type(mctc_error), allocatable :: merr
      integer :: k, i
      real(wp) :: r2

      do k = 1, 3
         allocate(xyz(3, npts_list(k)), w(npts_list(k)))
         call get_angular_grid(orders(k), xyz, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         do i = 1, npts_list(k)
            r2 = xyz(1, i)**2 + xyz(2, i)**2 + xyz(3, i)**2
            call check(error, abs(r2 - 1.0_wp) < 1.0e-12_wp, &
               & "Lebedev point off unit sphere")
            if (allocated(error)) exit
         end do
         deallocate(xyz, w)
         if (allocated(error)) return
      end do
   end subroutine test_angular_unit_sphere

   !> Chebyshev-2 radial quadrature (with r^2 dr Jacobian folded into w)
   !>
   !> integral_0^inf exp(-r^2) r^2 dr = sqrt(pi)/4; the chebyshev2_radii
   !> weights already include r^2 dr, so the sum reduces to
   !> sum_i w_i exp(-r_i^2)
   !>
   !> Measured absolute error 4.2e-13 at nr = 80
   !>
   !> @param[out] error  Test failure
   subroutine test_radial_chebyshev_gauss(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nr = 80
      real(wp) :: radii(nr), weights(nr)
      real(wp) :: result, expected
      integer :: ir

      call chebyshev2_radii(nr, 1.0_wp, radii, weights)
      result = 0.0_wp
      do ir = 1, nr
         result = result + weights(ir) * exp(-radii(ir)**2)
      end do
      expected = sqrt(pi) / 4.0_wp
      call check(error, abs(result - expected) < 3.0e-12_wp, &
         & "Chebyshev-2 quadrature of exp(-r^2) r^2 deviates from sqrt(pi)/4")
   end subroutine test_radial_chebyshev_gauss

   !> HandyMod radial weights should integrate the volume of a ball
   !>
   !> The bare radial weights include r^2 dr; multiplying their sum by 4*pi
   !> gives the 3D volume of a sphere when normalized Lebedev weights sum to 1
   !>
   !> @param[out] error  Test failure
   subroutine test_radial_handymod_ball_volume(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nr = 400
      integer, parameter :: m = 2
      real(wp), parameter :: rmin = 0.0_wp
      real(wp), parameter :: rmax = 6.0_wp
      real(wp), parameter :: rel_tol = 1.0e-4_wp
      real(wp) :: radii(nr), weights(nr)
      real(wp) :: radial_result, radial_expected
      real(wp) :: volume_result, volume_expected
      type(mctc_error), allocatable :: merr

      call handymod_radii(nr, rmin, rmax, m, radii, weights, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call check(error, all(radii >= rmin .and. radii <= rmax), &
         & "HandyMod radii lie outside the requested interval")
      if (allocated(error)) return
      call check(error, all(weights > 0.0_wp), &
         & "HandyMod radial weights should be positive")
      if (allocated(error)) return

      radial_result = sum(weights)
      radial_expected = (rmax**3 - rmin**3) / 3.0_wp
      call check(error, abs(radial_result - radial_expected) < rel_tol * radial_expected, &
         & "HandyMod radial weights do not reproduce the ball radial volume")
      if (allocated(error)) return

      volume_result = 4.0_wp * pi * radial_result
      volume_expected = 4.0_wp * pi * rmax**3 / 3.0_wp
      call check(error, abs(volume_result - volume_expected) < rel_tol * volume_expected, &
         & "HandyMod radial weights do not reproduce the ball volume")
   end subroutine test_radial_handymod_ball_volume

   !> HandyMod radial weights should integrate the volume of a finite shell
   !>
   !> A nonzero lower bound exercises the rmin offset in the radial transform
   !>
   !> @param[out] error  Test failure
   subroutine test_radial_handymod_shell_volume(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nr = 200
      integer, parameter :: m = 2
      real(wp), parameter :: rmin = 0.5_wp
      real(wp), parameter :: rmax = 4.0_wp
      real(wp), parameter :: rel_tol = 1.0e-4_wp
      real(wp) :: radii(nr), weights(nr)
      real(wp) :: radial_result, radial_expected
      real(wp) :: volume_result, volume_expected
      type(mctc_error), allocatable :: merr

      call handymod_radii(nr, rmin, rmax, m, radii, weights, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call check(error, all(radii >= rmin .and. radii <= rmax), &
         & "HandyMod shell radii lie outside the requested interval")
      if (allocated(error)) return
      call check(error, all(weights > 0.0_wp), &
         & "HandyMod shell radial weights should be positive")
      if (allocated(error)) return

      radial_result = sum(weights)
      radial_expected = (rmax**3 - rmin**3) / 3.0_wp
      call check(error, abs(radial_result - radial_expected) < rel_tol * radial_expected, &
         & "HandyMod radial weights do not reproduce the shell radial volume")
      if (allocated(error)) return

      volume_result = 4.0_wp * pi * radial_result
      volume_expected = 4.0_wp * pi * (rmax**3 - rmin**3) / 3.0_wp
      call check(error, abs(volume_result - volume_expected) < rel_tol * volume_expected, &
         & "HandyMod radial weights do not reproduce the shell volume")
   end subroutine test_radial_handymod_shell_volume

   !> Midpoint HandyMod weights should integrate a ball volume
   !>
   !> The fractional m = 0.1 transform has an endpoint-singular derivative;
   !> checks the expected midpoint convergence at the requested nr = 2000
   !>
   !> @param[out] error  Test failure
   subroutine test_radial_handymod_midpoint_ball(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nr = 2000
      real(wp), parameter :: rmin = 0.0_wp
      real(wp), parameter :: rmax = 6.0_wp
      real(wp), parameter :: m = 0.1_wp
      real(wp), parameter :: rel_tol = 1.0e-3_wp
      real(wp) :: radii(nr), weights(nr)
      real(wp) :: result, expected
      type(mctc_error), allocatable :: merr

      call handymod_midpoint_radii(nr, rmin, rmax, m, radii, weights, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call check(error, all(radii >= rmin .and. radii <= rmax), &
         & "Midpoint HandyMod radii lie outside the requested ball interval")
      if (allocated(error)) return
      call check(error, all(weights > 0.0_wp), &
         & "Midpoint HandyMod ball weights should be positive")
      if (allocated(error)) return

      result = sum(weights)
      expected = (rmax**3 - rmin**3) / 3.0_wp
      call check(error, abs(result - expected) < rel_tol * expected, &
         & "Midpoint HandyMod radial weights do not reproduce the ball volume")
   end subroutine test_radial_handymod_midpoint_ball

   !> Midpoint HandyMod weights should integrate a shell volume
   !>
   !> With rmin > 0 and m = 0.1, the endpoint-singular derivative is more
   !> pronounced, so the tolerance reflects the requested nr = 2000 midpoint rule
   !>
   !> @param[out] error  Test failure
   subroutine test_radial_handymod_midpoint_shell(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nr = 2000
      real(wp), parameter :: rmin = 0.5_wp
      real(wp), parameter :: rmax = 4.0_wp
      real(wp), parameter :: m = 0.1_wp
      real(wp), parameter :: rel_tol = 2.0e-2_wp
      real(wp) :: radii(nr), weights(nr)
      real(wp) :: result, expected
      type(mctc_error), allocatable :: merr

      call handymod_midpoint_radii(nr, rmin, rmax, m, radii, weights, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call check(error, all(radii >= rmin .and. radii <= rmax), &
         & "Midpoint HandyMod radii lie outside the requested shell interval")
      if (allocated(error)) return
      call check(error, all(weights > 0.0_wp), &
         & "Midpoint HandyMod shell weights should be positive")
      if (allocated(error)) return

      result = sum(weights)
      expected = (rmax**3 - rmin**3) / 3.0_wp
      call check(error, abs(result - expected) < rel_tol * expected, &
         & "Midpoint HandyMod radial weights do not reproduce the shell volume")
   end subroutine test_radial_handymod_midpoint_shell

   !> Midpoint HandyMod must reject nonpositive m via error_type
   !>
   !> @param[out] error  Test failure
   subroutine test_radial_handymod_midpoint_rejects_m(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nr = 16
      real(wp) :: radii(nr), weights(nr)
      type(mctc_error), allocatable :: merr

      call handymod_midpoint_radii(nr, 0.0_wp, 6.0_wp, 0.0_wp, radii, weights, merr)
      call check(error, allocated(merr), &
         & "Midpoint HandyMod should reject m = 0 through error_type")
      if (allocated(merr)) deallocate(merr)
      if (allocated(error)) return

      call handymod_midpoint_radii(nr, 0.0_wp, 6.0_wp, -0.1_wp, radii, weights, merr)
      call check(error, allocated(merr), &
         & "Midpoint HandyMod should reject negative m through error_type")
   end subroutine test_radial_handymod_midpoint_rejects_m

   !> Becke partition weights must sum to 1 at every sample point
   !>
   !> The molecule is MB16-43/H2 (any 2-atom system would do; the
   !> partition-of-unity property is geometry-independent)
   !>
   !> @param[out] error  Test failure
   subroutine test_becke_partition_of_unity(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :)
      integer,  allocatable :: numbers(:)
      real(wp) :: samples(3, 5)
      real(wp) :: weights(2)
      integer  :: k

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      xyz = mol%xyz
      allocate(numbers(mol%nat))
      do k = 1, mol%nat
         numbers(k) = mol%num(mol%id(k))
      end do

      samples(:, 1) = [ 0.1_wp,  0.2_wp, -0.3_wp]
      samples(:, 2) = [ 1.0_wp,  0.0_wp,  0.0_wp]
      samples(:, 3) = [-0.5_wp,  0.4_wp,  0.7_wp]
      samples(:, 4) = [ 0.0_wp,  0.0_wp,  3.0_wp]
      samples(:, 5) = [ 2.0_wp, -1.0_wp,  0.5_wp]

      do k = 1, size(samples, 2)
         call becke_weights(samples(:, k), mol%nat, xyz, numbers, weights)
         call check(error, abs(sum(weights) - 1.0_wp) < 1.0e-12_wp, &
            & "Becke weights do not sum to 1")
         if (allocated(error)) return
         call check(error, weights(1) >= -1.0e-14_wp .and. weights(2) >= -1.0e-14_wp, &
            & "Becke weights unexpectedly negative")
         if (allocated(error)) return
      end do
   end subroutine test_becke_partition_of_unity

end module test_math_quadrature
