!> Rigid-motion invariance of scalar integrals on all four 3D grid variants
!>
!> - Every integrand and grid setting from the integration suite
!> - Translate and rotate both carrier nuclei and the scalar field
!> - Update the grid from the moved molecule, retaining lab-frame grid directions
!> - Compare both integration APIs with the original integral and its reference
!> - Local callback context; no mutable module state under parallel test execution
module test_math_grid_3d_invariance
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use test_helpers, only: check_moist_error
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use test_math_grid_3d_integration, only: integration_case_type, make_cases, make_grid, &
      & grid_names, integral_tol
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_invariance

   !> Identity rotation
   real(wp), parameter :: identity(3, 3) = reshape([ &
      & 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
   !> Rotation pivot distinct from the origin, nuclei and integrand centers (bohr)
   real(wp), parameter :: pivot(3) = [1.13_wp, -0.71_wp, 0.37_wp]

   !> Active rigid motion r' = pivot + rotation*(r-pivot) + shift
   type :: rigid_motion_type
      !> Motion label in assertion failures
      character(len=32) :: name
      !> Proper orthogonal rotation, shape (3, 3)
      real(wp) :: rotation(3, 3)
      !> Translation after rotation (bohr)
      real(wp) :: shift(3)
   end type rigid_motion_type

contains

   !> Collect geometries; each runs every motion, integrand and grid variant
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_math_grid_3d_invariance(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("single_atom", test_single_atom), &
         & new_unittest("asymmetric_trimer", test_asymmetric_trimer)]
   end subroutine collect_math_grid_3d_invariance

   !> Off-center fields on a translated single-atom carrier
   subroutine test_single_atom(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [8], reshape([0.19_wp, -0.13_wp, 0.23_wp], [3, 1]))
      call check_invariance(mol, error)
   end subroutine test_single_atom

   !> Partitioned quadrature on an asymmetric heteronuclear trimer
   subroutine test_asymmetric_trimer(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [8, 1, 6], reshape([-0.6_wp, -0.2_wp, 0.1_wp, &
         & 0.7_wp, 0.3_wp, -0.25_wp, 0.1_wp, 1.2_wp, 0.4_wp], [3, 3]))
      call check_invariance(mol, error)
   end subroutine test_asymmetric_trimer

   !> Deterministic motions away from Cartesian and Lebedev rotation symmetries
   !>
   !> - Fractional shift has no component commensurate with the Cartesian spacing
   !> - Large shift moves the fields outside the original Cartesian box
   !> - Two oblique rotations and their noncommuting composition about an offset pivot
   !> - Final identity update checks recovery after the combined motion
   !>
   !> @param[out] motions Ordered translations, rotations and restoration
   subroutine make_motions(motions)
      !> Rigid motions applied independently to the original geometry
      type(rigid_motion_type), allocatable, intent(out) :: motions(:)
      real(wp) :: first(3, 3), second(3, 3)
      !> Fractional-grid translation (bohr)
      real(wp), parameter :: small_shift(3) = [0.317_wp, -0.431_wp, 0.283_wp]
      !> Translation larger than the original box half-width (bohr)
      real(wp), parameter :: large_shift(3) = [43.17_wp, -37.41_wp, 29.63_wp]
      !> Zero translation
      real(wp), parameter :: zero(3) = 0.0_wp

      first = axis_rotation([2.0_wp, -1.0_wp, 3.0_wp], 0.63_wp)
      second = axis_rotation([-1.0_wp, 4.0_wp, 2.0_wp], -0.91_wp)
      motions = [ &
         & rigid_motion_type("fractional translation", identity, small_shift), &
         & rigid_motion_type("large translation", identity, large_shift), &
         & rigid_motion_type("first oblique rotation", first, zero), &
         & rigid_motion_type("second oblique rotation", second, zero), &
         & rigid_motion_type("rotation and translation", matmul(second, first), large_shift), &
         & rigid_motion_type("restored geometry", identity, zero) &
         & ]
   end subroutine make_motions

   !> Rodrigues rotation about a nonzero axis
   !>
   !> @param[in] axis Rotation axis, shape (3)
   !> @param[in] angle Rotation angle (radians)
   pure function axis_rotation(axis, angle) result(rotation)
      !> Nonzero direction of the rotation axis
      real(wp), intent(in) :: axis(3)
      !> Rotation angle
      real(wp), intent(in) :: angle
      !> Proper orthogonal rotation matrix
      real(wp) :: rotation(3, 3)
      real(wp) :: u(3), cross(3, 3)
      integer :: j

      u = axis/sqrt(sum(axis**2))
      cross = reshape([0.0_wp, u(3), -u(2), -u(3), 0.0_wp, u(1), u(2), -u(1), 0.0_wp], [3, 3])
      rotation = cos(angle)*identity + sin(angle)*cross
      do j = 1, 3
         rotation(:, j) = rotation(:, j) + (1.0_wp - cos(angle))*u*u(j)
      end do
   end function axis_rotation

   !> Update each grid through the rigid motions and compare scalar integrals
   !>
   !> Cartesian and Lebedev grids keep fixed lab-frame directions, so arbitrary
   !> rotations preserve the continuum integral only to quadrature accuracy
   !>
   !> @param[in] mol Original carrier geometry
   !> @param[out] error Test failure
   subroutine check_invariance(mol, error)
      !> Original carrier geometry
      type(structure_type), intent(in) :: mol
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      class(moist_math_grid_3d_type), allocatable :: grid
      type(integration_case_type), allocatable :: cases(:)
      type(rigid_motion_type), allocatable :: motions(:)
      type(rigid_motion_type) :: original
      type(structure_type) :: moved
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: baseline(:, :)
      real(wp) :: integrals(2), gram(3, 3)
      character(len=:), allocatable :: label
      integer :: kind, icase, imotion, i, j

      call make_cases(cases)
      call make_motions(motions)
      original = rigid_motion_type("original geometry", identity, [0.0_wp, 0.0_wp, 0.0_wp])
      allocate (baseline(2, size(cases)))
      do imotion = 1, size(motions)
         gram = matmul(transpose(motions(imotion)%rotation), motions(imotion)%rotation)
         do j = 1, 3
            do i = 1, 3
               call check(error, gram(i, j), identity(i, j), thr=1.0e-14_wp, &
                  & more=trim(motions(imotion)%name)//": rotation orthogonality")
               if (allocated(error)) return
            end do
         end do
      end do

      do kind = 1, 4
         call make_grid(kind, mol, grid, error)
         if (allocated(error)) return
         do icase = 1, size(cases)
            label = trim(grid_names(kind))//": original: "//trim(cases(icase)%name)
            call integrate_moved_field(grid, cases(icase), original, baseline(:, icase), label, error)
            if (allocated(error)) return
         end do
         do imotion = 1, size(motions)
            moved = mol
            do i = 1, mol%nat
               moved%xyz(:, i) = pivot + matmul(motions(imotion)%rotation, mol%xyz(:, i) - pivot) &
                  & + motions(imotion)%shift
            end do
            call grid%update(moved, merr)
            label = trim(grid_names(kind))//": "//trim(motions(imotion)%name)
            call check_moist_error(error, merr, label//": update")
            if (allocated(error)) return
            do icase = 1, size(cases)
               label = trim(grid_names(kind))//": "//trim(motions(imotion)%name)//": "//trim(cases(icase)%name)
               call integrate_moved_field(grid, cases(icase), motions(imotion), integrals, label, error)
               if (allocated(error)) return
               call check(error, integrals(1), baseline(1, icase), thr_abs=integral_tol, thr_rel=integral_tol, &
                  & more=label//": integrate invariance")
               if (allocated(error)) return
               call check(error, integrals(2), baseline(2, icase), thr_abs=integral_tol, thr_rel=integral_tol, &
                  & more=label//": integrate_field invariance")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine check_invariance

   !> Integrate the transported field f'(r') = f(pivot + transpose(R)*(r'-pivot-t))
   !>
   !> Proper rotations have unit Jacobian; the whole-space reference is unchanged
   !>
   !> @param[in] grid Updated grid for the moved molecule
   !> @param[in] field Original scalar field and its independent integral
   !> @param[in] motion Active rigid motion
   !> @param[out] integrals Callback and tabulated integrals, shape (2)
   !> @param[in] label Assertion context
   !> @param[out] error Test failure
   subroutine integrate_moved_field(grid, field, motion, integrals, label, error)
      !> Grid regenerated from the moved molecule
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Original function and whole-space reference
      type(integration_case_type), intent(in) :: field
      !> Motion of the nuclei and scalar field
      type(rigid_motion_type), intent(in) :: motion
      !> Results from integrate and integrate_field
      real(wp), intent(out) :: integrals(2)
      !> Grid, motion and function labels
      character(len=*), intent(in) :: label
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: values(:)
      type(mctc_error), allocatable :: merr
      integer :: i

      call check(error, ieee_is_finite(field%reference), label//": non-finite reference")
      if (allocated(error)) return
      allocate (values(grid%ngrid))
      do i = 1, grid%ngrid
         values(i) = moved_integrand(grid%point(i))
      end do
      call check(error, all(ieee_is_finite(values)), label//": non-finite samples")
      if (allocated(error)) return
      call grid%integrate(moved_integrand, integrals(1))
      call grid%integrate_field(values, integrals(2), merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, integrals(1), field%reference, thr_abs=integral_tol, thr_rel=integral_tol, &
         & more=label//": integrate reference")
      if (allocated(error)) return
      call check(error, integrals(2), field%reference, thr_abs=integral_tol, thr_rel=integral_tol, &
         & more=label//": integrate_field reference")

   contains

      !> Pull the lab-frame evaluation point back into the original field frame
      !>
      !> @param[in] r Position on the updated grid (bohr), shape (3)
      pure function moved_integrand(r) result(value)
         !> Position in the transformed frame
         real(wp), intent(in) :: r(3)
         !> Scalar field value
         real(wp) :: value

         value = field%f(pivot + matmul(transpose(motion%rotation), r - pivot - motion%shift))
      end function moved_integrand

   end subroutine integrate_moved_field

end module test_math_grid_3d_invariance
