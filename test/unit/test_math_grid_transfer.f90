!> Transfer of fields between discretizations: math_grid_transfer
!>
!> - Radial spline on finer, coarser and non-uniform nodes, below the first stored node,
!>   and coverage beyond the last; negative nodes evaluated at |r|
!> - Rigid fit: a rotated and translated molecule, a pure translation (exact identity),
!>   one atom, a diatomic and a collinear triatomic (no rotation about the axis),
!>   nearly and exactly opposite axes, a bent structure fitted onto a straight one
!> - Coverage of a Cartesian box
!> - Malformed point layouts (two mode counts, a zero mode spacing, a NaN phase
!>   reference, too few weights, a NaN coordinate or a negative owner) found by
!>   `validate`, without a Fourier basis and covering nothing; each refused by the
!>   spectral transfer and by `prepare` in a test of its own, as are targets of two
!>   coordinates and a NaN radial node
!> - Spectral transfer: exact on a Cartesian lattice, off-lattice points of a periodic
!>   field, a larger box sharing the lattice
!> - Grid transfer: copies onto a larger box with the caller's values kept outside, an
!>   incomplete transfer without a coverage mask refused, several columns as each alone,
!>   the source layout no longer needed after `prepare`, radial copy and spline, exact
!>   copies only; each refused input is a test of its own
!> - Cartesian sources: trilinear onto point targets and Fourier onto lattices by default,
!>   either on request; trilinear reproduces affine fields and covers what Fourier covers
!> - Exact copies: k-d tree neighbours equal brute force; identical, shuffled and subset
!>   point sets copied bit for bit; coinciding source points, and lattice points closer
!>   than the tolerance, not copied and counted;
!>   lattice points copied and the others interpolated; the copy tolerance far below the
!>   smallest point spacing of atom-centered grids
!> - Point layouts (moving least squares): Gaussians on a molecular grid onto a Cartesian
!>   lattice, closer than the direct Fourier transform; affine fields reproduced; the
!>   error falls under refinement; flat and far neighbourhoods and strongly amplifying
!>   extrapolation not covered, mild extrapolation covered; coinciding source points;
!>   the same values for any thread count and with stencils fitted again when applied;
!>   the check of coverage that stops at the first stencil
module test_math_grid_transfer
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_grid_layout, only: moist_math_grid_layout_type
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   use moist_math_grid_transfer, only: rigid_fit_type, fit_rigid, transfer_radial, layout_covers, &
      & transfer_fourier, has_fourier_basis, grid_transfer_type, grid_transfer_options_type, &
      & grid_transfer_report_type, copy_tolerance
   use moist_math_kdtree, only: moist_kdtree_type
   use test_helpers, only: get_qc_handymod_recipe, check_moist_error
   implicit none(type, external)
   private

   public :: collect_math_grid_transfer

contains

   !> Collect the transfer tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_math_grid_transfer(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("radial_spline", test_radial_spline), &
         & new_unittest("radial_nonuniform", test_radial_nonuniform), &
         & new_unittest("fit_rotation", test_fit_rotation), &
         & new_unittest("fit_translation", test_fit_translation), &
         & new_unittest("fit_collinear", test_fit_collinear), &
         & new_unittest("fit_opposite", test_fit_opposite), &
         & new_unittest("fit_axes", test_fit_axes), &
         & new_unittest("fit_bent_straight", test_fit_bent_straight), &
         & new_unittest("covers_box", test_covers_box), &
         & new_unittest("layout_validate", test_layout_validate), &
         & new_unittest("fourier_malformed_modes", test_fourier_malformed_modes, should_fail=.true.), &
         & new_unittest("fourier_malformed_kspacing", test_fourier_malformed_kspacing, should_fail=.true.), &
         & new_unittest("fourier_malformed_kref", test_fourier_malformed_kref, should_fail=.true.), &
         & new_unittest("fourier_malformed_weights", test_fourier_malformed_weights, should_fail=.true.), &
         & new_unittest("fourier_malformed_xyz", test_fourier_malformed_xyz, should_fail=.true.), &
         & new_unittest("fourier_malformed_owner", test_fourier_malformed_owner, should_fail=.true.), &
         & new_unittest("prepare_malformed_modes", test_prepare_malformed_modes, should_fail=.true.), &
         & new_unittest("prepare_malformed_kspacing", test_prepare_malformed_kspacing, should_fail=.true.), &
         & new_unittest("prepare_malformed_kref", test_prepare_malformed_kref, should_fail=.true.), &
         & new_unittest("prepare_malformed_weights", test_prepare_malformed_weights, should_fail=.true.), &
         & new_unittest("prepare_malformed_xyz", test_prepare_malformed_xyz, should_fail=.true.), &
         & new_unittest("prepare_malformed_owner", test_prepare_malformed_owner, should_fail=.true.), &
         & new_unittest("fourier_two_coordinates", test_fourier_two_coordinates, should_fail=.true.), &
         & new_unittest("radial_nan_node", test_radial_nan_node, should_fail=.true.), &
         & new_unittest("fourier_cartesian", test_fourier_cartesian), &
         & new_unittest("transfer_coverage", test_transfer_coverage), &
         & new_unittest("transfer_columns", test_transfer_columns), &
         & new_unittest("transfer_kept_geometry", test_transfer_kept_geometry), &
         & new_unittest("transfer_radial", test_transfer_radial), &
         & new_unittest("radial_negative", test_radial_negative), &
         & new_unittest("lattice_ambiguous", test_lattice_ambiguous), &
         & new_unittest("transfer_copy_only", test_transfer_copy_only), &
         & new_unittest("apply_other_columns", test_apply_other_columns, should_fail=.true.), &
         & new_unittest("apply_other_rows", test_apply_other_rows, should_fail=.true.), &
         & new_unittest("prepare_mixed_kinds", test_prepare_mixed_kinds, should_fail=.true.), &
         & new_unittest("prepare_unknown_method", test_prepare_unknown_method, should_fail=.true.), &
         & new_unittest("prepare_box_without_corner", test_prepare_box_without_corner, should_fail=.true.), &
         & new_unittest("apply_after_failed_prepare", test_apply_after_failed_prepare, should_fail=.true.), &
         & new_unittest("trilinear", test_trilinear), &
         & new_unittest("kdtree_brute", test_kdtree_brute), &
         & new_unittest("copy_points", test_copy_points), &
         & new_unittest("copy_duplicates", test_copy_duplicates), &
         & new_unittest("copy_mixed", test_copy_mixed), &
         & new_unittest("copy_tolerance", test_copy_tolerance), &
         & new_unittest("points_molecular", test_points_molecular), &
         & new_unittest("points_affine", test_points_affine), &
         & new_unittest("points_refinement", test_points_refinement), &
         & new_unittest("points_uncovered", test_points_uncovered), &
         & new_unittest("points_threads", test_points_threads), &
         & new_unittest("points_check_only", test_points_check_only), &
         & new_unittest("points_duplicates", test_points_duplicates), &
         & new_unittest("points_extrapolation", test_points_extrapolation) &
         & ]
   end subroutine collect_math_grid_transfer

   !> Malformed layouts are found by `validate`, give no Fourier basis and cover nothing
   subroutine test_layout_validate(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_layout_type) :: good, bad
      logical, allocatable :: covered(:)
      character(len=:), allocatable :: reason
      integer :: n, variant

      call valid_points(good)
      call good%validate(n, reason)
      call check(error, n == 4 .and. len(reason) == 0 .and. has_fourier_basis(good), "a complete point layout refused")
      if (allocated(error)) return
      do variant = 1, 6
         call malformed_points(good, variant, bad)
         call bad%validate(n, reason)
         call check(error, len(reason) > 0 .and. n == -1 .and. .not. has_fourier_basis(bad), &
            & more="malformed point layout accepted, variant "//char(48 + variant))
         if (allocated(error)) return
      end do
      ! Cartesian boxes: a short shape, a NaN corner and a zero spacing cover nothing
      do variant = 1, 3
         call malformed_box(variant, bad)
         call layout_covers(bad, good%xyz, covered)
         call check(error, size(covered) == 4 .and. .not. any(covered), &
            & more="a malformed box covers points, variant "//char(48 + variant))
         if (allocated(error)) return
      end do
   end subroutine test_layout_validate

   !> Four points with weights, owners and a Fourier basis
   !>
   !> @param[out] layout  point layout
   subroutine valid_points(layout)
      !> Point layout
      type(moist_math_grid_layout_type), intent(out) :: layout
      layout%kind = "points"
      layout%xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 1.0_wp], [3, 4])
      layout%w = [0.25_wp, 0.25_wp, 0.25_wp, 0.25_wp]
      layout%owner = [1, 1, 1, 1]
      layout%kmodes = [8, 8, 8]
      layout%kspacing = [0.5_wp, 0.5_wp, 0.5_wp]
      layout%kref = [0.0_wp, 0.0_wp, 0.0_wp]
   end subroutine valid_points

   !> A point layout with one defect: two mode counts, a zero mode spacing, a NaN phase
   !> reference, three weights, a NaN coordinate or a negative owner
   !>
   !> @param[in]  good     valid point layout
   !> @param[in]  variant  defect, 1 to 6
   !> @param[out] bad      malformed layout
   subroutine malformed_points(good, variant, bad)
      !> Valid point layout
      type(moist_math_grid_layout_type), intent(in) :: good
      !> Defect
      integer, intent(in) :: variant
      !> Malformed layout
      type(moist_math_grid_layout_type), intent(out) :: bad
      bad = good
      select case (variant)
      case (1)
         bad%kmodes = [8, 8]
      case (2)
         bad%kspacing(2) = 0.0_wp
      case (3)
         bad%kref(1) = ieee_value(1.0_wp, ieee_quiet_nan)
      case (4)
         bad%w = bad%w(1:3)
      case (5)
         bad%xyz(2, 3) = ieee_value(1.0_wp, ieee_quiet_nan)
      case default
         bad%owner(4) = -1
      end select
   end subroutine malformed_points

   !> A box of 4x4x4 points with one defect: a short shape, a NaN corner or a zero spacing
   !>
   !> @param[in]  variant  defect, 1 to 3
   !> @param[out] bad      malformed box
   subroutine malformed_box(variant, bad)
      !> Defect
      integer, intent(in) :: variant
      !> Malformed box
      type(moist_math_grid_layout_type), intent(out) :: bad
      bad%kind = "cartesian"
      bad%shape = [4, 4, 4]
      bad%spacing = 0.5_wp
      bad%origin = [-1.0_wp, -1.0_wp, -1.0_wp]
      select case (variant)
      case (1)
         bad%shape = [4, 4]
      case (2)
         bad%origin(3) = ieee_value(1.0_wp, ieee_quiet_nan)
      case default
         bad%spacing = 0.0_wp
      end select
   end subroutine malformed_box

   !> Spectral transfer from a malformed point layout
   !>
   !> @param[out] error    the refusal
   !> @param[in]  variant  defect of `malformed_points`
   subroutine fourier_malformed(error, variant)
      !> The refusal
      type(error_type), allocatable, intent(out) :: error
      !> Defect
      integer, intent(in) :: variant
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: good, bad
      real(wp), allocatable :: f(:, :), g(:, :)

      call valid_points(good)
      call malformed_points(good, variant, bad)
      allocate (f(4, 1), source=1.0_wp)
      call transfer_fourier(bad, f, good%xyz, g, err)
      call check_moist_error(error, err)
   end subroutine fourier_malformed

   !> Grid transfer prepared from a malformed point layout
   !>
   !> @param[out] error    the refusal
   !> @param[in]  variant  defect of `malformed_points`
   subroutine prepare_malformed(error, variant)
      !> The refusal
      type(error_type), allocatable, intent(out) :: error
      !> Defect
      integer, intent(in) :: variant
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: good, bad
      type(grid_transfer_type) :: transfer

      call valid_points(good)
      call malformed_points(good, variant, bad)
      call transfer%prepare(bad, good, err)
      call check_moist_error(error, err)
   end subroutine prepare_malformed

   !> Refused: spectral transfer from a layout with two mode counts
   subroutine test_fourier_malformed_modes(error)
      type(error_type), allocatable, intent(out) :: error
      call fourier_malformed(error, 1)
   end subroutine test_fourier_malformed_modes

   !> Refused: spectral transfer from a layout with a zero mode spacing
   subroutine test_fourier_malformed_kspacing(error)
      type(error_type), allocatable, intent(out) :: error
      call fourier_malformed(error, 2)
   end subroutine test_fourier_malformed_kspacing

   !> Refused: spectral transfer from a layout with a NaN phase reference
   subroutine test_fourier_malformed_kref(error)
      type(error_type), allocatable, intent(out) :: error
      call fourier_malformed(error, 3)
   end subroutine test_fourier_malformed_kref

   !> Refused: spectral transfer from a layout with three weights for four points
   subroutine test_fourier_malformed_weights(error)
      type(error_type), allocatable, intent(out) :: error
      call fourier_malformed(error, 4)
   end subroutine test_fourier_malformed_weights

   !> Refused: spectral transfer from a layout with a NaN coordinate
   subroutine test_fourier_malformed_xyz(error)
      type(error_type), allocatable, intent(out) :: error
      call fourier_malformed(error, 5)
   end subroutine test_fourier_malformed_xyz

   !> Refused: spectral transfer from a layout with a negative owner
   subroutine test_fourier_malformed_owner(error)
      type(error_type), allocatable, intent(out) :: error
      call fourier_malformed(error, 6)
   end subroutine test_fourier_malformed_owner

   !> Refused: transfer prepared from a layout with two mode counts
   subroutine test_prepare_malformed_modes(error)
      type(error_type), allocatable, intent(out) :: error
      call prepare_malformed(error, 1)
   end subroutine test_prepare_malformed_modes

   !> Refused: transfer prepared from a layout with a zero mode spacing
   subroutine test_prepare_malformed_kspacing(error)
      type(error_type), allocatable, intent(out) :: error
      call prepare_malformed(error, 2)
   end subroutine test_prepare_malformed_kspacing

   !> Refused: transfer prepared from a layout with a NaN phase reference
   subroutine test_prepare_malformed_kref(error)
      type(error_type), allocatable, intent(out) :: error
      call prepare_malformed(error, 3)
   end subroutine test_prepare_malformed_kref

   !> Refused: transfer prepared from a layout with three weights for four points
   subroutine test_prepare_malformed_weights(error)
      type(error_type), allocatable, intent(out) :: error
      call prepare_malformed(error, 4)
   end subroutine test_prepare_malformed_weights

   !> Refused: transfer prepared from a layout with a NaN coordinate
   subroutine test_prepare_malformed_xyz(error)
      type(error_type), allocatable, intent(out) :: error
      call prepare_malformed(error, 5)
   end subroutine test_prepare_malformed_xyz

   !> Refused: transfer prepared from a layout with a negative owner
   subroutine test_prepare_malformed_owner(error)
      type(error_type), allocatable, intent(out) :: error
      call prepare_malformed(error, 6)
   end subroutine test_prepare_malformed_owner

   !> Refused: spectral transfer to points of two coordinates
   subroutine test_fourier_two_coordinates(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: good
      real(wp), allocatable :: f(:, :), g(:, :)

      call valid_points(good)
      allocate (f(4, 1), source=1.0_wp)
      call transfer_fourier(good, f, good%xyz(1:2, :), g, err)
      call check_moist_error(error, err)
   end subroutine test_fourier_two_coordinates

   !> Refused: radial transfer from a NaN node
   subroutine test_radial_nan_node(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp), allocatable :: fr(:, :)
      logical, allocatable :: covered(:)

      call transfer_radial([0.1_wp, ieee_value(1.0_wp, ieee_quiet_nan), 0.3_wp], &
         & reshape([1.0_wp, 2.0_wp, 3.0_wp], [3, 1]), [0.2_wp], fr, covered, err)
      call check_moist_error(error, err)
   end subroutine test_radial_nan_node

   !> Smooth radial test function, even in r
   !>
   !> @param[in] r  radius
   elemental real(wp) function smooth(r)
      !> Radius
      real(wp), intent(in) :: r
      smooth = exp(-0.5_wp*r**2)*cos(r)
   end function smooth

   !> Midpoint nodes refined three times and coarsened twice reproduce a smooth function
   subroutine test_radial_spline(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp), allocatable :: r_stored(:), r(:), f(:, :)
      logical, allocatable :: covered(:)
      integer :: i

      r_stored = [((real(i, wp) - 0.5_wp)*0.1_wp, i=1, 120)]
      r = [((real(i, wp) - 0.5_wp)*0.1_wp/3.0_wp, i=1, 400)]
      call transfer_radial(r_stored, reshape([smooth(r_stored), 2.0_wp*smooth(r_stored)], [120, 2]), r, f, covered, err)
      call check(error, .not. allocated(err), "radial transfer failed")
      if (allocated(error)) return
      call check(error, all(covered(:359)) .and. .not. any(covered(361:)), "coverage ends at the last stored node")
      if (allocated(error)) return
      call check(error, maxval(abs(f(:359, 1) - smooth(r(:359)))) < 1.0e-5_wp &
         & .and. maxval(abs(f(:359, 2) - 2.0_wp*smooth(r(:359)))) < 2.0e-5_wp, "refined nodes off the function")
      if (allocated(error)) return
      call check(error, maxval(abs(f(:3, 1) - smooth(r(:3)))) < 1.0e-5_wp, "nodes below the first stored node off")
      if (allocated(error)) return
      call check(error, all(f(361:, :) == 0.0_wp), "uncovered nodes not zero")
      if (allocated(error)) return
      r = [((real(i, wp) - 0.5_wp)*0.2_wp, i=1, 50)]
      call transfer_radial(r_stored, reshape(smooth(r_stored), [120, 1]), r, f, covered, err)
      call check(error, .not. allocated(err) .and. all(covered), "coarser nodes not covered")
      if (allocated(error)) return
      call check(error, maxval(abs(f(:, 1) - smooth(r))) < 1.0e-5_wp, "coarser nodes off the function")
   end subroutine test_radial_spline

   !> Non-uniform stored nodes, denser near the origin
   subroutine test_radial_nonuniform(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp), allocatable :: r_stored(:), r(:), f(:, :)
      logical, allocatable :: covered(:)
      integer :: i

      r_stored = [(10.0_wp*(real(i, wp)/200.0_wp)**2, i=1, 200)]
      r = [((real(i, wp) - 0.5_wp)*0.05_wp, i=1, 200)]
      call transfer_radial(r_stored, reshape(smooth(r_stored), [200, 1]), r, f, covered, err)
      call check(error, .not. allocated(err) .and. all(covered), "non-uniform transfer failed")
      if (allocated(error)) return
      call check(error, maxval(abs(f(:, 1) - smooth(r))) < 1.0e-4_wp, "non-uniform nodes off the function")
      if (allocated(error)) return
      call transfer_radial([1.0_wp, 0.5_wp], reshape([1.0_wp, 2.0_wp], [2, 1]), r, f, covered, err)
      call check(error, allocated(err), "decreasing stored nodes accepted")
   end subroutine test_radial_nonuniform

   !> Rotation matrix of a unit quaternion
   !>
   !> @param[in] q  unit quaternion
   pure function rotation(q) result(rot)
      !> Unit quaternion
      real(wp), intent(in) :: q(4)
      !> Rotation
      real(wp) :: rot(3, 3)
      rot(1, :) = [q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2, 2.0_wp*(q(2)*q(3) - q(1)*q(4)), &
         & 2.0_wp*(q(2)*q(4) + q(1)*q(3))]
      rot(2, :) = [2.0_wp*(q(2)*q(3) + q(1)*q(4)), q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2, &
         & 2.0_wp*(q(3)*q(4) - q(1)*q(2))]
      rot(3, :) = [2.0_wp*(q(2)*q(4) - q(1)*q(3)), 2.0_wp*(q(3)*q(4) + q(1)*q(2)), &
         & q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2]
   end function rotation

   !> A rotated and translated molecule maps back onto the stored one
   subroutine test_fit_rotation(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rigid_fit_type) :: fit
      real(wp) :: stored(3, 4), target(3, 4), q(4)
      real(wp), allocatable :: mapped(:, :)

      stored = reshape([0.0_wp, 0.0_wp, 0.2_wp, 1.4_wp, 0.1_wp, -0.9_wp, -1.4_wp, 0.3_wp, -0.8_wp, &
         & 0.2_wp, 2.1_wp, 0.5_wp], [3, 4])
      q = [0.8_wp, 0.3_wp, -0.4_wp, 0.33_wp]
      q = q/norm2(q)
      target = matmul(rotation(q), stored) + spread([1.5_wp, -0.7_wp, 3.2_wp], 2, 4)
      call fit_rigid(stored, target, fit, err)
      call check(error, .not. allocated(err) .and. .not. fit%is_identity(), "rotation not fitted")
      if (allocated(error)) return
      call fit%map(target, mapped)
      call check(error, maxval(abs(mapped - stored)) < 1.0e-10_wp, "fitted map misses the stored atoms")
      if (allocated(error)) return

      ! A map built by hand turns by its rotation: 90 degrees about z takes x to y
      fit = rigid_fit_type()
      fit%rot = reshape([0.0_wp, 1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      fit%c_target = [1.0_wp, 0.0_wp, 0.0_wp]
      fit%c_stored = [0.0_wp, 2.0_wp, 0.0_wp]
      call fit%map(reshape([2.0_wp, 0.0_wp, 0.0_wp], [3, 1]), mapped)
      call check(error, .not. fit%is_identity() .and. maxval(abs(mapped(:, 1) - [0.0_wp, 3.0_wp, 0.0_wp])) < 1.0e-15_wp, &
         & "a rotation set by hand was not applied")
   end subroutine test_fit_rotation

   !> A translation keeps the identity exactly; one atom is translated only
   subroutine test_fit_translation(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rigid_fit_type) :: fit
      real(wp) :: stored(3, 3), shift(3)
      real(wp), allocatable :: mapped(:, :)

      stored = reshape([0.0_wp, 0.0_wp, 0.2_wp, 1.4_wp, 0.1_wp, -0.9_wp, -1.4_wp, 0.3_wp, -0.8_wp], [3, 3])
      shift = [0.25_wp, -1.0_wp, 0.5_wp]
      call fit_rigid(stored, stored + spread(shift, 2, 3), fit, err)
      call check(error, .not. allocated(err) .and. fit%is_identity(), "translation not the identity")
      if (allocated(error)) return
      call fit%map(stored + spread(shift, 2, 3), mapped)
      call check(error, all(abs(mapped - stored) <= 1.0e-15_wp*4.0_wp), "translated atoms not mapped back")
      if (allocated(error)) return
      call fit_rigid(stored(:, 1:1), stored(:, 2:2), fit, err)
      call check(error, .not. allocated(err) .and. fit%is_identity(), "one atom rotated")
      if (allocated(error)) return
      call fit%map(stored(:, 2:2), mapped)
      call check(error, maxval(abs(mapped(:, 1) - stored(:, 1))) < 1.0e-14_wp, "one atom not translated")
   end subroutine test_fit_translation

   !> Collinear atoms: the axes are aligned by the smallest rotation, none about the axis
   subroutine test_fit_collinear(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rigid_fit_type) :: fit
      real(wp) :: stored(3, 3), target(3, 3), axis_s(3), axis_t(3), normal(3)
      real(wp), allocatable :: mapped(:, :)
      integer :: n

      axis_s = [0.0_wp, 0.0_wp, 1.0_wp]
      axis_t = [0.6_wp, 0.0_wp, 0.8_wp]
      ! The rotation about y turns axis_t into axis_s and leaves y alone
      normal = [0.0_wp, 1.0_wp, 0.0_wp]
      do n = 2, 3
         stored(:, 1) = -1.1_wp*axis_s
         stored(:, 2) = 1.1_wp*axis_s
         stored(:, 3) = 2.5_wp*axis_s
         target(:, 1) = -1.1_wp*axis_t + [1.0_wp, 2.0_wp, 3.0_wp]
         target(:, 2) = 1.1_wp*axis_t + [1.0_wp, 2.0_wp, 3.0_wp]
         target(:, 3) = 2.5_wp*axis_t + [1.0_wp, 2.0_wp, 3.0_wp]
         call fit_rigid(stored(:, :n), target(:, :n), fit, err)
         call check(error, .not. allocated(err), "collinear fit failed")
         if (allocated(error)) return
         call fit%map(target(:, :n), mapped)
         call check(error, maxval(abs(mapped - stored(:, :n))) < 1.0e-10_wp, "collinear atoms not mapped")
         if (allocated(error)) return
         call check(error, maxval(abs(matmul(fit%rot, normal) - normal)) < 1.0e-12_wp, "rotation about the axis")
         if (allocated(error)) return
      end do
   end subroutine test_fit_collinear

   !> Opposite and nearly opposite axes of a diatomic give a finite proper rotation
   subroutine test_fit_opposite(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rigid_fit_type) :: fit
      real(wp) :: stored(3, 2), target(3, 2), axis(3)
      real(wp), parameter :: tilts(3) = [0.0_wp, 1.0e-9_wp, 1.0e-4_wp]
      real(wp), allocatable :: mapped(:, :)
      integer :: k

      do k = 1, 3
         axis = [tilts(k), 0.0_wp, -1.0_wp]
         axis = axis/norm2(axis)
         stored(:, 1) = [0.0_wp, 0.0_wp, -0.7_wp]
         stored(:, 2) = [0.0_wp, 0.0_wp, 0.7_wp]
         target(:, 1) = -0.7_wp*axis + [0.3_wp, 0.0_wp, 0.0_wp]
         target(:, 2) = 0.7_wp*axis + [0.3_wp, 0.0_wp, 0.0_wp]
         ! The first atom of the target sits where the second stored one points
         call fit_rigid(stored, target, fit, err)
         call check(error, .not. allocated(err) .and. all(abs(fit%rot) <= 1.0_wp + 1.0e-12_wp), &
            & "opposite axes not fitted")
         if (allocated(error)) return
         call check(error, proper(fit%rot), "opposite axes give no proper rotation")
         if (allocated(error)) return
         call fit%map(target, mapped)
         call check(error, maxval(abs(mapped - stored)) < 1.0e-8_wp, "opposite diatomic not mapped")
         if (allocated(error)) return
      end do
   end subroutine test_fit_opposite

   !> Collinear structures at many angles: the turn about the common normal
   !>
   !> - Polar angles from 30 to 179 degrees at several azimuths
   !> - The target atoms land on the stored ones, the rotation is proper, the common
   !>   normal of the two axes stays in place
   subroutine test_fit_axes(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rigid_fit_type) :: fit
      real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
      real(wp), parameter :: polar(5) = [30.0_wp, 90.0_wp, 120.0_wp, 150.0_wp, 179.0_wp]*pi/180.0_wp
      real(wp), parameter :: azimuth(4) = [0.0_wp, 1.1_wp, 2.5_wp, 4.0_wp]
      real(wp) :: stored(3, 3), target(3, 3), u(3), normal(3)
      real(wp), allocatable :: mapped(:, :)
      integer :: i, j

      stored = reshape([0.0_wp, 0.0_wp, -1.2_wp, 0.0_wp, 0.0_wp, 0.4_wp, 0.0_wp, 0.0_wp, 2.0_wp], [3, 3])
      do i = 1, size(polar)
         do j = 1, size(azimuth)
            u = [sin(polar(i))*cos(azimuth(j)), sin(polar(i))*sin(azimuth(j)), cos(polar(i))]
            target(:, 1) = -1.2_wp*u
            target(:, 2) = 0.4_wp*u
            target(:, 3) = 2.0_wp*u
            call fit_rigid(stored, target, fit, err)
            call check(error, .not. allocated(err) .and. proper(fit%rot), "collinear fit gives no proper rotation")
            if (allocated(error)) return
            call fit%map(target, mapped)
            call check(error, maxval(abs(mapped - stored)) < 1.0e-10_wp, "collinear atoms not mapped")
            if (allocated(error)) return
            normal = [u(2), -u(1), 0.0_wp]
            normal = normal/norm2(normal)
            call check(error, maxval(abs(matmul(fit%rot, normal) - normal)) < 1.0e-10_wp, &
               & "the common normal of the axes was turned")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_fit_axes

   !> A bent structure fitted onto a straight one takes the full least-squares fit
   subroutine test_fit_bent_straight(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rigid_fit_type) :: fit, turned
      real(wp) :: stored(3, 3), target(3, 3), best, tilt(3, 3), angle
      real(wp), allocatable :: mapped(:, :)
      integer :: axis, sgn

      stored = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.4_wp, 1.1_wp, 0.0_wp, -1.4_wp, 1.1_wp, 0.0_wp], [3, 3])
      ! A straight chain with the first atom at one end
      target = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.08_wp, 0.0_wp, 1.44_wp, 2.16_wp, 0.0_wp, 2.88_wp], [3, 3])
      call fit_rigid(stored, target, fit, err)
      call check(error, .not. allocated(err) .and. proper(fit%rot), "bent onto straight gives no proper rotation")
      if (allocated(error)) return
      call fit%map(target, mapped)
      best = sum((mapped - stored)**2)
      ! No small extra rotation about the stored centroid improves the fit
      angle = 1.0e-3_wp
      do axis = 1, 3
         do sgn = -1, 1, 2
            tilt = small_rotation(axis, sgn*angle)
            turned = fit
            turned%rot = matmul(tilt, fit%rot)
            call turned%map(target, mapped)
            call check(error, sum((mapped - stored)**2) >= best - 1.0e-12_wp, "the fit is not a least-squares optimum")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_fit_bent_straight

   !> Whether a matrix is a proper rotation
   !>
   !> @param[in] rot  matrix
   pure logical function proper(rot)
      !> Matrix
      real(wp), intent(in) :: rot(3, 3)
      real(wp) :: det
      det = rot(1, 1)*(rot(2, 2)*rot(3, 3) - rot(2, 3)*rot(3, 2)) &
         & - rot(1, 2)*(rot(2, 1)*rot(3, 3) - rot(2, 3)*rot(3, 1)) &
         & + rot(1, 3)*(rot(2, 1)*rot(3, 2) - rot(2, 2)*rot(3, 1))
      proper = maxval(abs(matmul(transpose(rot), rot) - reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
         & 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3]))) < 1.0e-12_wp .and. abs(det - 1.0_wp) < 1.0e-12_wp
   end function proper

   !> Rotation by `angle` about the Cartesian axis `axis`
   !>
   !> @param[in] axis   Cartesian axis, 1 to 3
   !> @param[in] angle  angle, radian
   pure function small_rotation(axis, angle) result(rot)
      !> Cartesian axis, 1 to 3
      integer, intent(in) :: axis
      !> Angle, radian
      real(wp), intent(in) :: angle
      real(wp) :: rot(3, 3)
      integer :: i, j
      rot = reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      i = modulo(axis, 3) + 1
      j = modulo(axis + 1, 3) + 1
      rot(i, i) = cos(angle)
      rot(j, j) = cos(angle)
      rot(i, j) = -sin(angle)
      rot(j, i) = sin(angle)
   end function small_rotation

   !> A Cartesian layout covers its box
   subroutine test_covers_box(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_layout_type) :: layout
      logical, allocatable :: covered(:)

      layout%kind = "cartesian"
      layout%shape = [4, 4, 2]
      layout%spacing = 0.5_wp
      layout%origin = [-1.0_wp, -1.0_wp, 0.0_wp]
      call layout_covers(layout, reshape([0.0_wp, 0.0_wp, 0.5_wp, 0.99_wp, -0.99_wp, 0.01_wp, &
         & 1.2_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.0_wp, 1.1_wp], [3, 4]), covered)
      call check(error, all(covered .eqv. [.true., .true., .false., .false.]), "box coverage")
   end subroutine test_covers_box

   !> Periodic test field on the box [origin, origin + n dr], resolved by its lattice
   !>
   !> @param[in] x       point
   !> @param[in] origin  box corner
   !> @param[in] box     box edge lengths
   pure real(wp) function periodic(x, origin, box) result(f)
      !> Point
      real(wp), intent(in) :: x(3)
      !> Box corner and edge lengths
      real(wp), intent(in) :: origin(3), box(3)
      real(wp), parameter :: two_pi = 8.0_wp*atan(1.0_wp)
      real(wp) :: u(3)
      u = two_pi*(x - origin)/box
      f = 0.3_wp + cos(u(1) + 0.4_wp) + 0.5_wp*sin(2.0_wp*u(2))*cos(u(3)) - 0.25_wp*cos(3.0_wp*u(1) - u(3))
   end function periodic

   !> A Cartesian lattice transfers exactly: lattice points, off-lattice points, a larger box
   subroutine test_fourier_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: layout, larger
      real(wp), allocatable :: xs(:, :), f(:, :), y(:, :), g(:, :), exact(:)
      logical, allocatable :: covered(:)
      real(wp) :: box(3)
      integer :: i

      layout%kind = "cartesian"
      layout%shape = [16, 12, 10]
      layout%spacing = 0.4_wp
      layout%origin = [-3.2_wp, -2.4_wp, -2.0_wp]
      box = real(layout%shape, wp)*layout%spacing
      call check(error, has_fourier_basis(layout), "Cartesian layout without Fourier basis")
      if (allocated(error)) return
      call layout%points(xs)
      allocate (f(size(xs, 2), 2))
      do i = 1, size(xs, 2)
         f(i, 1) = periodic(xs(:, i), layout%origin, box)
      end do
      f(:, 2) = -2.0_wp*f(:, 1)
      ! Off-lattice points inside the box
      allocate (y(3, 50))
      do i = 1, 50
         y(:, i) = layout%origin + box*[modulo(0.618_wp*i, 1.0_wp), modulo(0.414_wp*i, 1.0_wp), modulo(0.732_wp*i, 1.0_wp)]
      end do
      call transfer_fourier(layout, f, y, g, err, tol=1.0e-12_wp)
      call check(error, .not. allocated(err), "Cartesian transfer failed")
      if (allocated(error)) return
      allocate (exact(50))
      do i = 1, 50
         exact(i) = periodic(y(:, i), layout%origin, box)
      end do
      call check(error, maxval(abs(g(:, 1) - exact)) < 1.0e-9_wp .and. maxval(abs(g(:, 2) + 2.0_wp*exact)) < 2.0e-9_wp, &
         & "off-lattice values of a resolved field")
      if (allocated(error)) return
      ! A larger box two cells wider on each side shares the lattice inside the old one
      larger%kind = "cartesian"
      larger%shape = layout%shape + 4
      larger%spacing = layout%spacing
      larger%origin = layout%origin - 2.0_wp*layout%spacing
      call larger%points(y)
      call layout_covers(layout, y, covered)
      call check(error, count(covered) == product(layout%shape), "the old box covers its own lattice")
      if (allocated(error)) return
      y = y(:, pack([(i, i=1, size(y, 2))], covered))
      call transfer_fourier(layout, f, y, g, err, tol=1.0e-12_wp)
      call check(error, .not. allocated(err), "transfer to the larger box failed")
      if (allocated(error)) return
      call check(error, maxval(abs(g(:, 1) - f(:, 1))) < 1.0e-9_wp, "shared lattice points changed")
   end subroutine test_fourier_cartesian

   !> Cartesian test box of 16 x 12 x 10 points, spacing 0.4 bohr
   !>
   !> @param[out] layout  layout of the box
   pure subroutine test_box(layout)
      !> Layout
      type(moist_math_grid_layout_type), intent(out) :: layout
      layout%kind = "cartesian"
      layout%shape = [16, 12, 10]
      layout%spacing = 0.4_wp
      layout%origin = [-3.2_wp, -2.4_wp, -2.0_wp]
   end subroutine test_box

   !> Columns f, -2 f and f**2 of the periodic test field at points of a box
   !>
   !> @param[in] layout  box
   !> @param[in] x       points (3, n)
   pure function box_field(layout, x) result(f)
      !> Box
      type(moist_math_grid_layout_type), intent(in) :: layout
      !> Points
      real(wp), intent(in) :: x(:, :)
      !> Columns (n, 3)
      real(wp), allocatable :: f(:, :)
      integer :: i
      allocate (f(size(x, 2), 3))
      do i = 1, size(x, 2)
         f(i, 1) = periodic(x(:, i), layout%origin, real(layout%shape, wp)*layout%spacing)
      end do
      f(:, 2) = -2.0_wp*f(:, 1)
      f(:, 3) = f(:, 1)**2
   end function box_field

   !> Off-lattice points inside the test box
   !>
   !> @param[in] layout  box
   !> @param[in] n       point count
   pure function inner_points(layout, n) result(y)
      !> Box
      type(moist_math_grid_layout_type), intent(in) :: layout
      !> Point count
      integer, intent(in) :: n
      !> Points (3, n)
      real(wp), allocatable :: y(:, :)
      integer :: i
      allocate (y(3, n))
      do i = 1, n
         y(:, i) = layout%origin + real(layout%shape, wp)*layout%spacing &
            & *[modulo(0.618_wp*i, 1.0_wp), modulo(0.414_wp*i, 1.0_wp), modulo(0.732_wp*i, 1.0_wp)]
      end do
   end function inner_points

   !> A larger box sharing the lattice: the old box is copied, the caller's values stay
   !> outside, and without a coverage mask the incomplete transfer is refused
   subroutine test_transfer_coverage(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: xs(:, :), f(:, :), g(:, :)
      logical, allocatable :: covered(:)
      integer :: nold

      call test_box(old)
      call old%points(xs)
      f = box_field(old, xs)
      nold = size(xs, 2)
      new = old
      new%shape = old%shape + 2
      new%origin = old%origin - old%spacing
      call transfer%prepare(old, new, err)
      call check(error, .not. allocated(err), "prepare failed")
      if (allocated(error)) return
      allocate (g(product(new%shape), 3), source=-7.0_wp)
      call transfer%apply(f, g, err)
      call check(error, allocated(err), "an incomplete transfer without a mask was accepted")
      if (allocated(error)) return
      call check(error, all(g == -7.0_wp), "a refused transfer changed the field")
      if (allocated(error)) return
      call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err), "apply failed")
      if (allocated(error)) return
      call check(error, report%method == "copy" .and. report%ncopied == nold .and. &
         & report%nuncovered == size(g, 1) - nold .and. report%reason == "outside the source grid", &
         & "report "//report%method//", "//report%reason)
      if (allocated(error)) return
      call check(error, count(covered) == nold, "covered points")
      if (allocated(error)) return
      ! The shared lattice keeps its values exactly, in the order of the old box
      call check(error, all(reshape(pack(g, spread(covered, 2, 3)), [nold, 3]) == f), "copied values changed")
      if (allocated(error)) return
      call check(error, all(pack(g, spread(.not. covered, 2, 3)) == -7.0_wp), "the caller's values outside changed")
   end subroutine test_transfer_coverage

   !> Several columns move as each column alone, and off-lattice values match the field
   subroutine test_transfer_columns(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_options_type) :: options
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: xs(:, :), f(:, :), g(:, :), g1(:, :), exact(:, :)
      integer :: col

      call test_box(old)
      call old%points(xs)
      f = box_field(old, xs)
      new%kind = "points"
      new%xyz = inner_points(old, 40)
      exact = box_field(old, new%xyz)
      options%nthreads = 1
      options%cartesian = "fourier"
      call transfer%prepare(old, new, err, options=options)
      call check(error, .not. allocated(err), "prepare failed")
      if (allocated(error)) return
      allocate (g(40, 3), g1(40, 1), source=0.0_wp)
      call transfer%apply(f, g, err, report=report)
      call check(error, .not. allocated(err), "apply failed")
      if (allocated(error)) return
      call check(error, report%method == "fourier" .and. report%ninterpolated == 40 .and. report%nuncovered == 0, &
         & "report "//report%method)
      if (allocated(error)) return
      call check(error, maxval(abs(g(:, 1:2) - exact(:, 1:2))) < 2.0e-8_wp, "off-lattice values of the field")
      if (allocated(error)) return
      do col = 1, 3
         call transfer%apply(f(:, col:col), g1, err)
         call check(error, .not. allocated(err), "apply of one column failed")
         if (allocated(error)) return
         call check(error, maxval(abs(g1(:, 1) - g(:, col))) <= 1.0e-13_wp*maxval(abs(g(:, col))), &
            & "a column moved alone differs")
         if (allocated(error)) return
      end do
   end subroutine test_transfer_columns

   !> After `prepare` the source layout may change or go away
   subroutine test_transfer_kept_geometry(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer, fresh
      real(wp), allocatable :: xs(:, :), f(:, :), g(:, :), expected(:, :)

      call test_box(old)
      call old%points(xs)
      f = box_field(old, xs)
      new%kind = "points"
      new%xyz = inner_points(old, 25)
      call fresh%prepare(old, new, err)
      if (.not. allocated(err)) call transfer%prepare(old, new, err)
      call check(error, .not. allocated(err), "prepare failed")
      if (allocated(error)) return
      allocate (expected(25, 3), g(25, 3), source=0.0_wp)
      call fresh%apply(f, expected, err)
      call fresh%destroy()
      deallocate (old%origin, old%shape)
      old%kind = "radial"
      deallocate (new%xyz)
      if (.not. allocated(err)) call transfer%apply(f, g, err)
      call check(error, .not. allocated(err), "apply after the layouts changed failed")
      if (allocated(error)) return
      call check(error, all(g == expected), "the transfer depends on the layouts after prepare")
      if (allocated(error)) return
      call transfer%destroy()
      call transfer%apply(f, g, err)
      call check(error, allocated(err), "a destroyed transfer was applied")
   end subroutine test_transfer_kept_geometry

   !> Radial nodes: the same nodes are copied, other nodes splined up to the last stored node
   subroutine test_transfer_radial(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g(:, :), expected(:, :)
      logical, allocatable :: covered(:), inside(:)
      integer :: i

      old%kind = "radial"
      old%r = [((real(i, wp) - 0.5_wp)*0.1_wp, i=1, 80)]
      f = reshape([smooth(old%r), 3.0_wp*smooth(old%r)], [80, 2])
      call transfer%prepare(old, old, err)
      call check(error, .not. allocated(err), "prepare on the same nodes failed")
      if (allocated(error)) return
      allocate (g(80, 2), source=0.0_wp)
      call transfer%apply(f, g, err, report=report)
      call check(error, .not. allocated(err) .and. report%method == "copy" .and. all(g == f), "same nodes not copied")
      if (allocated(error)) return
      new%kind = "radial"
      new%r = [((real(i, wp) - 0.5_wp)*0.07_wp, i=1, 130)]
      call transfer%prepare(old, new, err)
      call check(error, .not. allocated(err), "prepare on other nodes failed")
      if (allocated(error)) return
      deallocate (g)
      allocate (g(130, 2), source=5.0_wp)
      call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err) .and. report%method == "spline", "other nodes not splined")
      if (allocated(error)) return
      call transfer_radial(old%r, f, new%r, expected, inside, err)
      call check(error, all(covered .eqv. inside) .and. .not. all(covered), "radial coverage")
      if (allocated(error)) return
      call check(error, all(pack(g, spread(covered, 2, 2)) == pack(expected, spread(covered, 2, 2))) &
         & .and. all(pack(g, spread(.not. covered, 2, 2)) == 5.0_wp), "splined values")
   end subroutine test_transfer_radial

   !> Negative radial targets are evaluated at |r| and covered only within the stored range
   subroutine test_radial_negative(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp), allocatable :: f(:, :), mirrored(:, :)
      logical, allocatable :: covered(:), inside(:)

      call transfer_radial([1.0_wp, 2.0_wp], reshape([1.0_wp, 3.0_wp], [2, 1]), [-3.0_wp, -1.5_wp], f, covered, err)
      call check(error, .not. allocated(err), "radial transfer to negative nodes failed")
      if (allocated(error)) return
      call check(error, .not. covered(1) .and. f(1, 1) == 0.0_wp, "a node beyond -r_max was covered")
      if (allocated(error)) return
      call transfer_radial([1.0_wp, 2.0_wp], reshape([1.0_wp, 3.0_wp], [2, 1]), [1.5_wp], mirrored, inside, err)
      call check(error, .not. allocated(err) .and. covered(2) .and. inside(1), "a node within the range not covered")
      if (allocated(error)) return
      call check(error, f(2, 1) == mirrored(1, 1), "-r and r differ")
   end subroutine test_radial_negative

   !> A target within the copy tolerance of two lattice points is interpolated, not copied
   !>
   !> - Spacing 0.5 bohr: the target on a lattice point is copied
   !> - Spacing 1e-10 bohr: its axis neighbours lie within the tolerance too
   subroutine test_lattice_ambiguous(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g(:, :)
      real(wp) :: spacing(2)
      integer :: i

      spacing = [0.5_wp, 1.0e-10_wp]
      allocate (f(64, 1), source=1.0_wp)
      do i = 1, 2
         old = moist_math_grid_layout_type()
         old%kind = "cartesian"
         old%shape = [4, 4, 4]
         old%spacing = spacing(i)
         old%origin = [0.0_wp, 0.0_wp, 0.0_wp]
         new = moist_math_grid_layout_type()
         new%kind = "points"
         new%xyz = reshape([1.5_wp, 1.5_wp, 1.5_wp]*spacing(i), [3, 1])
         call transfer%prepare(old, new, err)
         call check(error, .not. allocated(err), "prepare onto a lattice point failed")
         if (allocated(error)) return
         if (allocated(g)) deallocate (g)
         allocate (g(1, 1), source=0.0_wp)
         call transfer%apply(f, g, err, report=report)
         call check(error, .not. allocated(err), "transfer onto a lattice point failed")
         if (allocated(error)) return
         if (i == 1) then
            call check(error, report%ncopied == 1 .and. report%nduplicate == 0, "a lattice point was not copied")
         else
            call check(error, report%ncopied == 0 .and. report%nduplicate == 1 .and. report%ninterpolated == 1, &
               & "an ambiguous lattice point was copied")
         end if
         if (allocated(error)) return
         call check(error, abs(g(1, 1) - 1.0_wp) < 1.0e-12_wp, "constant field not reproduced")
         if (allocated(error)) return
      end do
   end subroutine test_lattice_ambiguous

   !> Exact copies only leave off-lattice targets uncovered and the caller's values kept
   subroutine test_transfer_copy_only(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g(:, :)
      logical, allocatable :: covered(:)

      call copy_only_fixture(transfer, f, g, err)
      call check(error, .not. allocated(err), "prepare with copies only failed")
      if (allocated(error)) return
      call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err) .and. .not. any(covered) .and. report%nexcluded == 10 .and. &
         & report%reason == "only exact copies requested" .and. all(g == 1.0_wp), "copies only")
   end subroutine test_transfer_copy_only

   !> Copies only from a 3D box onto ten points between its lattice points
   !>
   !> - Refusal tests return no test error for a failing fixture, which then shows as an
   !>   unexpected pass
   !>
   !> @param[out] transfer  prepared transfer
   !> @param[out] f         field on the box (npoint, 3)
   !> @param[out] g         target field of ones (10, 3)
   !> @param[out] err       failed prepare
   subroutine copy_only_fixture(transfer, f, g, err)
      !> Prepared transfer
      type(grid_transfer_type), intent(out) :: transfer
      !> Field on the box
      real(wp), allocatable, intent(out) :: f(:, :)
      !> Target field
      real(wp), allocatable, intent(out) :: g(:, :)
      !> Failed prepare
      type(moist_error), allocatable, intent(out) :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_options_type) :: options
      real(wp), allocatable :: xs(:, :)

      call test_box(old)
      call old%points(xs)
      f = box_field(old, xs)
      new%kind = "points"
      new%xyz = inner_points(old, 10)
      options%copy_only = .true.
      allocate (g(10, 3), source=1.0_wp)
      call transfer%prepare(old, new, err, options=options)
   end subroutine copy_only_fixture

   !> Refused: a target field with other columns than the source field
   subroutine test_apply_other_columns(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(grid_transfer_type) :: transfer
      real(wp), allocatable :: f(:, :), g(:, :)
      logical, allocatable :: covered(:)

      call copy_only_fixture(transfer, f, g, err)
      if (allocated(err)) return
      call transfer%apply(f(:, 1:2), g, err, covered=covered)
      call check_moist_error(error, err)
   end subroutine test_apply_other_columns

   !> Refused: a source field of other rows than the source layout
   subroutine test_apply_other_rows(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(grid_transfer_type) :: transfer
      real(wp), allocatable :: f(:, :), g(:, :)
      logical, allocatable :: covered(:)

      call copy_only_fixture(transfer, f, g, err)
      if (allocated(err)) return
      call transfer%apply(f(2:, :), g, err, covered=covered)
      call check_moist_error(error, err)
   end subroutine test_apply_other_rows

   !> Refused: a transfer from a 3D box onto radial nodes
   subroutine test_prepare_mixed_kinds(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, radial
      type(grid_transfer_type) :: transfer

      call test_box(old)
      radial%kind = "radial"
      radial%r = [0.5_wp, 1.0_wp]
      call transfer%prepare(old, radial, err)
      call check_moist_error(error, err)
   end subroutine test_prepare_mixed_kinds

   !> Refused: an unknown method for Cartesian sources
   subroutine test_prepare_unknown_method(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_options_type) :: options

      call test_box(old)
      new%kind = "points"
      new%xyz = inner_points(old, 10)
      options%cartesian = "cubic"
      call transfer%prepare(old, new, err, options=options)
      call check_moist_error(error, err)
   end subroutine test_prepare_unknown_method

   !> Refused: a box without its corner
   subroutine test_prepare_box_without_corner(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer

      call test_box(old)
      new%kind = "points"
      new%xyz = inner_points(old, 10)
      deallocate (old%origin)
      call transfer%prepare(old, new, err)
      call check_moist_error(error, err)
   end subroutine test_prepare_box_without_corner

   !> Refused: applying a transfer whose prepare failed
   subroutine test_apply_after_failed_prepare(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      real(wp), allocatable :: f(:, :), g(:, :)
      logical, allocatable :: covered(:)

      call copy_only_fixture(transfer, f, g, err)
      if (allocated(err)) return
      call test_box(old)
      new%kind = "points"
      new%xyz = inner_points(old, 10)
      deallocate (old%origin)
      call transfer%prepare(old, new, err)
      ! A prepare that succeeds leaves nothing to refuse: an unexpected pass
      if (.not. allocated(err)) return
      deallocate (err)
      call transfer%apply(f, g, err, covered=covered)
      call check_moist_error(error, err)
   end subroutine test_apply_after_failed_prepare

   !> Cartesian sources: trilinear by default onto point targets, Fourier onto lattices,
   !> either on request; trilinear reproduces affine fields between the lattice points,
   !> copies the lattice points and covers what Fourier covers
   subroutine test_trilinear(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new, lattice
      type(grid_transfer_type) :: transfer
      type(grid_transfer_options_type) :: options
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: xs(:, :), f(:, :), g(:, :), exact(:, :), y(:, :)
      logical, allocatable :: covered(:), covered_fourier(:)
      real(wp) :: hull_lo(3), hull_hi(3), edge(3)
      integer :: i, n

      call test_box(old)
      call old%points(xs)
      f = affine(xs)
      ! Off-lattice points between the outermost lattice points, lattice points, and points
      ! in the half cell at the faces and outside the box
      hull_lo = old%origin + 0.5_wp*old%spacing
      hull_hi = old%origin + (real(old%shape, wp) - 0.5_wp)*old%spacing
      edge = real(old%shape, wp)*old%spacing
      allocate (y(3, 60))
      do i = 1, 40
         y(:, i) = hull_lo + (hull_hi - hull_lo)*[modulo(0.618_wp*i, 1.0_wp), modulo(0.414_wp*i, 1.0_wp), &
            & modulo(0.732_wp*i, 1.0_wp)]
      end do
      y(:, 41:50) = xs(:, 1:300:30)
      do i = 51, 60
         y(:, i) = old%origin + edge*[modulo(0.3_wp*i, 1.0_wp), 0.5_wp, 1.0_wp] + [0.0_wp, 0.0_wp, 0.3_wp*real(i - 55, wp)]
      end do
      new%kind = "points"
      new%xyz = y
      exact = affine(y)
      allocate (g(60, 2), source=0.0_wp)
      call transfer%prepare(old, new, err)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err), "trilinear transfer failed")
      if (allocated(error)) return
      call check(error, report%method == "copy+trilinear" .and. report%ncopied == 10, "methods "//report%method)
      if (allocated(error)) return
      call check(error, all(covered(1:50)) .and. maxval(abs(g(1:40, :) - exact(1:40, :))) <= 1.0e-12_wp*maxval(abs(exact)), &
         & "affine field not reproduced between the lattice points")
      if (allocated(error)) return
      call check(error, all(g(41:50, :) == f(1:300:30, :)), "lattice points not copied")
      if (allocated(error)) return
      ! Coverage equals the coverage of the Fourier transfer
      options%cartesian = "fourier"
      call transfer%prepare(old, new, err, options=options)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered_fourier, report=report)
      call check(error, .not. allocated(err) .and. report%method == "copy+fourier", "Fourier on request")
      if (allocated(error)) return
      call check(error, all(covered .eqv. covered_fourier) .and. .not. all(covered), "coverage differs from Fourier")
      if (allocated(error)) return
      ! A lattice target takes Fourier by default and trilinear on request
      lattice = old
      lattice%origin = old%origin + 0.3_wp
      n = product(lattice%shape)
      deallocate (g)
      allocate (g(n, 2), source=0.0_wp)
      call transfer%prepare(old, lattice, err)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err) .and. report%method == "fourier", "lattice target: "//report%method)
      if (allocated(error)) return
      options%cartesian = "trilinear"
      call transfer%prepare(old, lattice, err, options=options)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err) .and. report%method == "trilinear", "trilinear on request: "//report%method)
   contains
      !> Two affine columns
      !>
      !> @param[in] x  points (3, n)
      pure function affine(x) result(v)
         !> Points (3, n)
         real(wp), intent(in) :: x(:, :)
         real(wp), allocatable :: v(:, :)
         allocate (v(size(x, 2), 2))
         v(:, 1) = 0.5_wp - 0.7_wp*x(1, :) + 0.2_wp*x(2, :) + 1.3_wp*x(3, :)
         v(:, 2) = x(3, :)
      end function affine
   end subroutine test_trilinear

   !> Pseudo-random points in a box of edge `edge` around the origin, the same every run
   !>
   !> @param[in]  n     point count
   !> @param[in]  edge  box edge, bohr
   !> @param[out] x     points (3, n)
   subroutine scattered_points(n, edge, x)
      !> Point count
      integer, intent(in) :: n
      !> Box edge, bohr
      real(wp), intent(in) :: edge
      !> Points (3, n)
      real(wp), allocatable, intent(out) :: x(:, :)
      integer :: i
      allocate (x(3, n))
      do i = 1, n
         x(:, i) = edge*([modulo(0.7548776662_wp*i, 1.0_wp), modulo(0.5698402910_wp*i, 1.0_wp), &
            & modulo(0.3247179572_wp*real(i, wp)**1.5_wp, 1.0_wp)] - 0.5_wp)
      end do
   end subroutine scattered_points

   !> The k nearest points from the tree equal a brute-force search, ties by index,
   !> also with repeated points and with fewer points than asked for
   subroutine test_kdtree_brute(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_kdtree_type) :: tree
      real(wp), allocatable :: x(:, :), y(:, :), d2(:), dist2(:)
      integer, allocatable :: order(:), nearest(:)
      integer :: kk(3), ik, k, iq, i, j, t

      call scattered_points(700, 6.0_wp, x)
      ! Repeated points: equal distances, ordered by index
      x(:, 650:700) = x(:, 1:51)
      call scattered_points(60, 7.0_wp, y)
      y(:, 1:5) = x(:, 10:14)
      call tree%build(x)
      kk = [1, 7, 40]
      allocate (d2(size(x, 2)), order(size(x, 2)))
      do ik = 1, 3
         k = kk(ik)
         allocate (nearest(k), dist2(k))
         do iq = 1, size(y, 2)
            call tree%knn(y(:, iq), k, nearest, dist2)
            do i = 1, size(x, 2)
               d2(i) = sum((x(:, i) - y(:, iq))**2)
               order(i) = i
            end do
            ! Insertion sort by distance, then index
            do i = 2, size(x, 2)
               t = order(i)
               j = i - 1
               do while (j >= 1)
                  if (d2(order(j)) < d2(t) .or. (d2(order(j)) == d2(t) .and. order(j) < t)) exit
                  order(j + 1) = order(j)
                  j = j - 1
               end do
               order(j + 1) = t
            end do
            call check(error, all(nearest == order(1:k)) .and. all(dist2 == d2(order(1:k))), &
               & "k-d tree neighbours differ from brute force")
            if (allocated(error)) return
         end do
         deallocate (nearest, dist2)
      end do
      call tree%build(x(:, 1:3))
      allocate (nearest(5), dist2(5))
      call tree%knn(y(:, 7), 5, nearest, dist2)
      call check(error, all(nearest(4:) == 0) .and. all(dist2(4:) == huge(1.0_wp)) .and. all(nearest(:3) > 0), &
         & "fewer points than asked for")
   end subroutine test_kdtree_brute

   !> Identical, shuffled and subset point sets are copied bit for bit; a point moved by
   !> more than the tolerance is not, by less it is
   subroutine test_copy_points(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g(:, :)
      logical, allocatable :: covered(:)
      integer, allocatable :: order(:)
      integer :: i, n

      n = 300
      old%kind = "points"
      call scattered_points(n, 8.0_wp, old%xyz)
      allocate (f(n, 2))
      f(:, 1) = [(sin(real(i, wp)), i=1, n)]
      f(:, 2) = [(real(i, wp), i=1, n)]
      ! Identical
      new = old
      call transfer%prepare(old, new, err)
      allocate (g(n, 2), source=0.0_wp)
      if (.not. allocated(err)) call transfer%apply(f, g, err, report=report)
      call check(error, .not. allocated(err), "copy of identical points failed")
      if (allocated(error)) return
      call check(error, report%method == "copy" .and. report%ncopied == n .and. all(g == f), "identical points")
      if (allocated(error)) return
      ! Shuffled
      order = [(modulo(37*i, n) + 1, i=0, n - 1)]
      new%xyz = old%xyz(:, order)
      call transfer%prepare(old, new, err)
      if (.not. allocated(err)) call transfer%apply(f, g, err)
      call check(error, .not. allocated(err) .and. all(g == f(order, :)), "shuffled points")
      if (allocated(error)) return
      ! Subset, with one point moved beyond the tolerance and one within it; copies only
      new%xyz = old%xyz(:, 1:n:3)
      new%xyz(1, 2) = new%xyz(1, 2) + 1.0e-7_wp
      new%xyz(2, 3) = new%xyz(2, 3) + 0.2_wp*copy_tolerance
      deallocate (g)
      allocate (g(size(new%xyz, 2), 2), source=-1.0_wp)
      call transfer%prepare(old, new, err, options=grid_transfer_options_type(copy_only=.true.))
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err), "copy of a subset failed")
      if (allocated(error)) return
      call check(error, count(.not. covered) == 1 .and. .not. covered(2) .and. all(g(2, :) == -1.0_wp), &
         & "a point beyond the tolerance was copied")
      if (allocated(error)) return
      covered(2) = .true.
      g(2, :) = f(4, :)
      call check(error, all(g == f(1:n:3, :)), "subset values")
      if (allocated(error)) return
      ! Without copies only, the moved point is interpolated
      call transfer%prepare(old, new, err)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err) .and. all(covered) .and. report%method == "copy+points" &
         & .and. report%ninterpolated == 1, "the moved point was not interpolated: "//report%method)
   end subroutine test_copy_points

   !> Targets at two coinciding source points are not copied and are counted, on point
   !> layouts and on radial nodes; on point layouts the fit takes both points
   subroutine test_copy_duplicates(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g(:, :)
      logical, allocatable :: covered(:)
      integer :: i

      old%kind = "points"
      call scattered_points(50, 5.0_wp, old%xyz)
      old%xyz(:, 50) = old%xyz(:, 7) + [0.3_wp, 0.0_wp, 0.0_wp]*copy_tolerance
      allocate (f(50, 1))
      f(:, 1) = [(real(i, wp), i=1, 50)]
      new%kind = "points"
      new%xyz = old%xyz(:, 1:20)
      allocate (g(20, 1), source=0.0_wp)
      call transfer%prepare(old, new, err)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err), "transfer with coinciding source points failed")
      if (allocated(error)) return
      call check(error, report%nduplicate == 1 .and. report%ncopied == 19 .and. report%ninterpolated == 1 &
         & .and. all(covered), "coinciding source points copied: "//report%method)
      if (allocated(error)) return
      ! Radial nodes closer than the tolerance
      old%kind = "radial"
      old%r = [0.5_wp, 1.0_wp, 1.0_wp + 0.5_wp*copy_tolerance, 2.0_wp]
      new%kind = "radial"
      new%r = [0.5_wp, 1.0_wp, 2.0_wp]
      deallocate (g, f)
      allocate (f(4, 1), g(3, 1), source=0.0_wp)
      f(:, 1) = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      options_radial: block
         type(grid_transfer_options_type) :: options
         options%copy_only = .true.
         call transfer%prepare(old, new, err, options=options)
      end block options_radial
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err), "radial transfer with coinciding nodes failed")
      if (allocated(error)) return
      call check(error, report%nduplicate == 1 .and. all(covered .eqv. [.true., .false., .true.]) &
         & .and. g(1, 1) == 1.0_wp .and. g(3, 1) == 4.0_wp, "coinciding radial nodes copied")
   end subroutine test_copy_duplicates

   !> Targets on the source lattice are copied, the others interpolated
   subroutine test_copy_mixed(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      type(grid_transfer_options_type) :: options
      real(wp), allocatable :: xs(:, :), f(:, :), g(:, :), exact(:, :)

      call test_box(old)
      call old%points(xs)
      f = box_field(old, xs)
      new%kind = "points"
      new%xyz = reshape([xs(:, 1:400:7), inner_points(old, 30)], [3, 88])
      exact = box_field(old, new%xyz)
      allocate (g(88, 3), source=0.0_wp)
      options%cartesian = "fourier"
      call transfer%prepare(old, new, err, options=options)
      if (.not. allocated(err)) call transfer%apply(f, g, err, report=report)
      call check(error, .not. allocated(err), "mixed transfer failed")
      if (allocated(error)) return
      call check(error, report%method == "copy+fourier" .and. report%ncopied == 58 .and. report%ninterpolated == 30, &
         & "methods "//report%method)
      if (allocated(error)) return
      call check(error, all(g(1:58, :) == f(1:400:7, :)), "lattice points not copied exactly")
      if (allocated(error)) return
      call check(error, maxval(abs(g(59:, 1:2) - exact(59:, 1:2))) < 2.0e-8_wp, "interpolated points off the field")
   end subroutine test_copy_mixed

   !> The copy tolerance lies below 1e-3 of the smallest point spacing of atom-centered
   !> grids: the default molecular settings, a fine grid, and the largest Lebedev table
   subroutine test_copy_tolerance(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_kdtree_type) :: tree
      integer, parameter :: nrad(3) = [50, 100, 100], degree(3) = [35, 77, 131]
      real(wp), parameter :: rmax(3) = [10.0_wp, 32.0_wp, 32.0_wp]
      real(wp) :: spacing, dist2(2)
      integer :: nearest(2), ig, i
      character(len=80) :: text

      call new(mol, [10], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      do ig = 1, 3
         call get_qc_handymod_recipe(recipe, err, nrad=nrad(ig), degree=degree(ig), rmax=rmax(ig))
         if (.not. allocated(err)) call new_molecular_point_grid(grid, err, recipe=recipe, pruning_threshold=0.0_wp, &
            & dr=2.0_wp)
         if (.not. allocated(err)) call grid%update(mol, err)
         call check(error, .not. allocated(err), "atom-centered grid not built")
         if (allocated(error)) return
         call tree%build(grid%xyz)
         spacing = huge(1.0_wp)
         do i = 1, grid%ngrid
            call tree%knn(grid%xyz(:, i), 2, nearest, dist2)
            spacing = min(spacing, sqrt(dist2(2)))
         end do
         write (text, "(a,i0,a,i0,a,es9.2,a)") "nrad ", nrad(ig), ", degree ", degree(ig), ": smallest spacing ", &
            & spacing, " bohr"
         call check(error, copy_tolerance <= 1.0e-3_wp*spacing, trim(text))
         if (allocated(error)) return
         call grid%destroy()
      end do
   end subroutine test_copy_tolerance

   !> Water molecule, bohr
   !>
   !> @param[out] mol  molecule
   subroutine water(mol)
      !> Molecule
      type(structure_type), intent(out) :: mol
      call new(mol, [8, 1, 1], reshape([0.0_wp, 0.0_wp, 0.12_wp, 1.43_wp, 0.0_wp, -0.98_wp, &
         & -1.43_wp, 0.0_wp, -0.98_wp], [3, 3]))
   end subroutine water

   !> Molecular grid of a water molecule and its layout
   !>
   !> @param[in]  mol     molecule
   !> @param[in]  nrad    radial shells
   !> @param[in]  degree  Lebedev degree
   !> @param[out] layout  layout of the grid
   !> @param[out] err     error handling
   subroutine water_grid(mol, nrad, degree, layout, err)
      !> Molecule
      type(structure_type), intent(in) :: mol
      !> Radial shells and Lebedev degree
      integer, intent(in) :: nrad, degree
      !> Layout of the grid
      type(moist_math_grid_layout_type), intent(out) :: layout
      !> Error handling
      type(moist_error), allocatable, intent(out) :: err
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid

      call get_qc_handymod_recipe(recipe, err, nrad=nrad, degree=degree, rmax=8.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(grid, err, recipe=recipe, pruning_threshold=0.1_wp, &
         & dr=0.8_wp)
      if (.not. allocated(err)) call grid%update(mol, err)
      if (.not. allocated(err)) call grid%get_layout(layout)
      call grid%destroy()
   end subroutine water_grid

   !> Cartesian lattice of n^3 points around a molecule
   !>
   !> @param[in]  mol      molecule
   !> @param[in]  n        points per axis
   !> @param[in]  spacing  spacing, bohr
   !> @param[out] layout   layout of the box
   !> @param[out] err      error handling
   subroutine water_box(mol, n, spacing, layout, err)
      !> Molecule
      type(structure_type), intent(in) :: mol
      !> Points per axis
      integer, intent(in) :: n
      !> Spacing, bohr
      real(wp), intent(in) :: spacing
      !> Layout of the box
      type(moist_math_grid_layout_type), intent(out) :: layout
      !> Error handling
      type(moist_error), allocatable, intent(out) :: err
      type(moist_math_grid_3d_cartesian_type) :: cart

      call new_cartesian_point_grid(cart, err, n, n, n, spacing, margin=0.0_wp)
      if (.not. allocated(err)) call cart%update(mol, err)
      if (.not. allocated(err)) call cart%get_layout(layout)
   end subroutine water_box

   !> Gaussians at the atoms of a molecule
   !>
   !> @param[in] mol  molecule
   !> @param[in] x    points (3, n)
   pure function gaussians(mol, x) result(v)
      !> Molecule
      type(structure_type), intent(in) :: mol
      !> Points (3, n)
      real(wp), intent(in) :: x(:, :)
      !> Values (n)
      real(wp), allocatable :: v(:)
      integer :: i
      allocate (v(size(x, 2)))
      do i = 1, size(x, 2)
         v(i) = exp(-sum((x(:, i) - mol%xyz(:, 1))**2)/(2.0_wp*1.2_wp**2)) &
            & + 0.4_wp*exp(-sum((x(:, i) - mol%xyz(:, 2))**2)/(2.0_wp*0.8_wp**2)) &
            & + 0.4_wp*exp(-sum((x(:, i) - mol%xyz(:, 3))**2)/(2.0_wp*0.8_wp**2))
      end do
   end function gaussians

   !> Relative RMS error of values against a reference
   !>
   !> @param[in] v    values
   !> @param[in] ref  reference
   pure real(wp) function relrms(v, ref)
      !> Values and reference
      real(wp), intent(in) :: v(:), ref(:)
      relrms = sqrt(sum((v - ref)**2)/max(sum(ref**2), tiny(1.0_wp)))
   end function relrms

   !> Gaussians on a water molecular grid onto a Cartesian lattice: moving least squares is
   !> closer than the direct Fourier transform of the grid (4.0e-2 against 0.26 relative
   !> RMS error); the grid onto itself is copied
   subroutine test_points_molecular(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_layout_type) :: grid, box
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g(:, :), direct(:, :), y(:, :), exact(:)
      logical, allocatable :: covered(:)
      real(wp) :: err_points, err_direct
      character(len=96) :: text
      integer :: i

      call water(mol)
      call water_grid(mol, 24, 11, grid, err)
      if (.not. allocated(err)) call water_box(mol, 16, 0.75_wp, box, err)
      call check(error, .not. allocated(err), "grids not built")
      if (allocated(error)) return
      f = reshape(gaussians(mol, grid%xyz), [size(grid%xyz, 2), 1])
      call transfer%prepare(grid, box, err)
      allocate (g(product(box%shape), 1), source=0.0_wp)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err), "transfer from the molecular grid failed")
      if (allocated(error)) return
      call check(error, report%method == "points" .and. count(covered) > 3*size(covered)/4, &
         & "method "//report%method)
      if (allocated(error)) return
      call box%points(y)
      y = y(:, pack([(i, i=1, size(covered))], covered))
      exact = gaussians(mol, y)
      err_points = relrms(pack(g(:, 1), covered), exact)
      call transfer_fourier(grid, f, y, direct, err, nthreads=2)
      call check(error, .not. allocated(err), "direct transform failed")
      if (allocated(error)) return
      err_direct = relrms(direct(:, 1), exact)
      write (text, "(a,es9.2,a,es9.2)") "points ", err_points, ", direct Fourier ", err_direct
      call check(error, err_points < 0.06_wp .and. err_points < 0.25_wp*err_direct, trim(text))
      if (allocated(error)) return
      ! The grid onto itself
      call transfer%prepare(grid, grid, err)
      deallocate (g)
      allocate (g(size(f, 1), 1), source=0.0_wp)
      if (.not. allocated(err)) call transfer%apply(f, g, err, report=report)
      call check(error, .not. allocated(err) .and. report%method == "copy" .and. all(g == f), "grid onto itself")
   end subroutine test_points_molecular

   !> Affine fields from a molecular grid are reproduced at every interpolated target, on
   !> a Cartesian lattice and on a molecular grid of other shells
   subroutine test_points_affine(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_layout_type) :: grid, targets(2)
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g(:, :), y(:, :), exact(:, :)
      logical, allocatable :: covered(:)
      integer :: it

      call water(mol)
      call water_grid(mol, 16, 11, grid, err)
      if (.not. allocated(err)) call water_box(mol, 12, 0.8_wp, targets(1), err)
      if (.not. allocated(err)) call water_grid(mol, 21, 7, targets(2), err)
      call check(error, .not. allocated(err), "grids not built")
      if (allocated(error)) return
      f = affine(grid%xyz)
      do it = 1, 2
         call targets(it)%points(y)
         exact = affine(y)
         call transfer%prepare(grid, targets(it), err)
         if (allocated(g)) deallocate (g)
         allocate (g(size(y, 2), 2), source=0.0_wp)
         if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
         call check(error, .not. allocated(err), "affine transfer failed")
         if (allocated(error)) return
         call check(error, report%ninterpolated > size(y, 2)/2, "too few interpolated targets")
         if (allocated(error)) return
         call check(error, maxval(abs(pack(g - exact, spread(covered, 2, 2)))) <= 1.0e-10_wp*maxval(abs(exact)), &
            & "affine field not reproduced")
         if (allocated(error)) return
      end do
   contains
      !> Two affine columns
      !>
      !> @param[in] x  points (3, n)
      pure function affine(x) result(v)
         !> Points (3, n)
         real(wp), intent(in) :: x(:, :)
         real(wp), allocatable :: v(:, :)
         allocate (v(size(x, 2), 2))
         v(:, 1) = 1.0_wp + 0.3_wp*x(1, :) - 0.2_wp*x(2, :) + 0.5_wp*x(3, :)
         v(:, 2) = 2.0_wp - x(1, :)
      end function affine
   end subroutine test_points_affine

   !> The error on Gaussians at the same targets falls when the molecular grid is refined
   subroutine test_points_refinement(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_layout_type) :: coarse, fine, targets
      real(wp), allocatable :: y(:, :), g_coarse(:), g_fine(:), exact(:)
      logical, allocatable :: covered_coarse(:), covered_fine(:), both(:)
      real(wp) :: err_coarse, err_fine
      character(len=96) :: text

      call water(mol)
      call water_grid(mol, 12, 7, coarse, err)
      if (.not. allocated(err)) call water_grid(mol, 24, 17, fine, err)
      if (.not. allocated(err)) call water_box(mol, 10, 0.6_wp, targets, err)
      call check(error, .not. allocated(err), "grids not built")
      if (allocated(error)) return
      call targets%points(y)
      exact = gaussians(mol, y)
      call gaussians_onto(coarse, g_coarse, covered_coarse)
      if (allocated(error)) return
      call gaussians_onto(fine, g_fine, covered_fine)
      if (allocated(error)) return
      both = covered_coarse .and. covered_fine
      call check(error, count(both) > size(both)/2, "few targets covered by both grids")
      if (allocated(error)) return
      err_coarse = relrms(pack(g_coarse, both), pack(exact, both))
      err_fine = relrms(pack(g_fine, both), pack(exact, both))
      write (text, "(a,es9.2,a,es9.2,a,i0)") "coarse ", err_coarse, ", fine ", err_fine, ", targets ", count(both)
      call check(error, err_fine < 0.5_wp*err_coarse, trim(text))
   contains
      !> The Gaussians of a source grid at the targets
      !>
      !> @param[in]  source   source layout
      !> @param[out] g        values at the targets
      !> @param[out] covered  covered targets
      subroutine gaussians_onto(source, g, covered)
         !> Source layout
         type(moist_math_grid_layout_type), intent(in) :: source
         !> Values at the targets
         real(wp), allocatable, intent(out) :: g(:)
         !> Covered targets
         logical, allocatable, intent(out) :: covered(:)
         type(grid_transfer_type) :: transfer
         real(wp), allocatable :: f(:, :), values(:, :)
         f = reshape(gaussians(mol, source%xyz), [size(source%xyz, 2), 1])
         allocate (values(size(y, 2), 1), source=0.0_wp)
         call transfer%prepare(source, targets, err)
         if (.not. allocated(err)) call transfer%apply(f, values, err, covered=covered)
         call check(error, .not. allocated(err), "transfer of the Gaussians failed")
         g = values(:, 1)
      end subroutine gaussians_onto
   end subroutine test_points_refinement

   !> Neighbourhoods a stencil cannot resolve leave their targets uncovered: coplanar and
   !> collinear source points, a tight cluster, a target far above a two-layer slab, and
   !> targets beyond the outer points
   subroutine test_points_uncovered(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g(:, :)
      logical, allocatable :: covered(:)
      integer :: i, j, l, n

      ! Coplanar: a square lattice in z = 0; targets above it and in the plane between points
      old%kind = "points"
      allocate (old%xyz(3, 121))
      n = 0
      do j = 1, 11
         do i = 1, 11
            n = n + 1
            old%xyz(:, n) = [0.5_wp*real(i - 6, wp), 0.5_wp*real(j - 6, wp), 0.0_wp]
         end do
      end do
      new%kind = "points"
      new%xyz = reshape([0.1_wp, 0.2_wp, 0.3_wp, 0.25_wp, 0.25_wp, 0.0_wp], [3, 2])
      call expect_uncovered("coplanar source points", "source neighbourhood flat or too small")
      if (allocated(error)) return
      ! Collinear: points along x
      old%xyz = reshape([(real(i, wp)*0.25_wp, 0.0_wp, 0.0_wp, i=1, 40)], [3, 40])
      new%xyz = reshape([2.1_wp, 0.1_wp, 0.0_wp], [3, 1])
      call expect_uncovered("collinear source points", "source neighbourhood flat or too small")
      if (allocated(error)) return
      ! A tight cluster: its points reach 4e-4 bohr at most
      old%xyz = reshape([(1.0e-3_wp*cos(real(i, wp)), 1.0e-3_wp*sin(real(i, wp)), 1.0e-4_wp*real(i, wp), i=1, 30)], [3, 30])
      new%xyz = reshape([0.5_wp, 0.0_wp, 0.0_wp], [3, 1])
      call expect_uncovered("a tight cluster", "covers no point of this grid")
      if (allocated(error)) return
      ! A slab of two layers 0.5 bohr apart: 1 bohr above the top layer the stencil
      ! extrapolates with weights near 3 and -2
      deallocate (old%xyz)
      allocate (old%xyz(3, 2*121))
      n = 0
      do l = 0, 1
         do j = 1, 11
            do i = 1, 11
               n = n + 1
               old%xyz(:, n) = [0.5_wp*real(i - 6, wp), 0.5_wp*real(j - 6, wp), 0.5_wp*real(l, wp)]
            end do
         end do
      end do
      new%xyz = reshape([0.05_wp, 0.05_wp, 1.5_wp], [3, 1])
      call expect_uncovered("a target far above a slab", "stencil amplification above its bound")
      if (allocated(error)) return
      ! Beyond the outer points of the slab
      new%xyz = reshape([0.0_wp, 0.0_wp, 4.0_wp, 9.0_wp, 0.0_wp, 0.25_wp], [3, 2])
      call expect_uncovered("targets beyond the outer points", "covers no point of this grid")
   contains
      !> Transfer onto `new` and expect no covered target and a reason
      !>
      !> @param[in] label   test label
      !> @param[in] reason  expected part of the reason
      subroutine expect_uncovered(label, reason)
         !> Test label and expected part of the reason
         character(len=*), intent(in) :: label, reason
         if (allocated(f)) deallocate (f, g)
         allocate (f(size(old%xyz, 2), 1), source=1.0_wp)
         allocate (g(size(new%xyz, 2), 1), source=0.0_wp)
         call transfer%prepare(old, new, err)
         if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
         call check(error, .not. allocated(err), label//": transfer failed")
         if (allocated(error)) return
         call check(error, .not. any(covered) .and. report%reason == reason, label//": "//report%reason)
      end subroutine expect_uncovered
   end subroutine test_points_uncovered

   !> Values and coverage do not depend on the thread count, and stencils fitted again
   !> when applied (no stencil memory) give the stored values
   subroutine test_points_threads(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_layout_type) :: grid, box
      type(grid_transfer_type) :: transfer
      type(grid_transfer_options_type) :: options
      real(wp), allocatable :: f(:, :), g(:, :, :)
      logical, allocatable :: covered(:, :), mask(:)
      integer :: run

      call water(mol)
      call water_grid(mol, 16, 11, grid, err)
      if (.not. allocated(err)) call water_box(mol, 14, 0.7_wp, box, err)
      call check(error, .not. allocated(err), "grids not built")
      if (allocated(error)) return
      f = reshape([gaussians(mol, grid%xyz), -gaussians(mol, grid%xyz)**2], [size(grid%xyz, 2), 2])
      allocate (g(product(box%shape), 2, 3), source=0.0_wp)
      allocate (covered(product(box%shape), 3))
      do run = 1, 3
         options = grid_transfer_options_type()
         options%nthreads = merge(1, 4, run == 1)
         if (run == 3) options%stencil_memory = 0
         call transfer%prepare(grid, box, err, options=options)
         if (.not. allocated(err)) call transfer%apply(f, g(:, :, run), err, covered=mask)
         call check(error, .not. allocated(err), "transfer failed")
         if (allocated(error)) return
         covered(:, run) = mask
      end do
      call check(error, all(g(:, :, 2) == g(:, :, 1)) .and. all(covered(:, 2) .eqv. covered(:, 1)), &
         & "values depend on the thread count")
      if (allocated(error)) return
      call check(error, all(g(:, :, 3) == g(:, :, 1)) .and. all(covered(:, 3) .eqv. covered(:, 1)), &
         & "stencils fitted again give other values")
   end subroutine test_points_threads

   !> The check of coverage: the copies and the first fitted target, nothing where no
   !> stencil is accepted although the targets are near; the transfer cannot be applied
   subroutine test_points_check_only(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_layout_type) :: grid, box, plane, above
      type(grid_transfer_type) :: transfer
      type(grid_transfer_options_type) :: options
      type(grid_transfer_report_type) :: fitted, checked
      real(wp), allocatable :: f(:, :), g(:, :)
      integer :: i, j, n

      call water(mol)
      call water_grid(mol, 16, 11, grid, err)
      if (.not. allocated(err)) call water_box(mol, 14, 0.7_wp, box, err)
      call check(error, .not. allocated(err), "grids not built")
      if (allocated(error)) return
      call transfer%prepare(grid, box, err)
      if (.not. allocated(err)) call transfer%get_report(fitted)
      options%check_only = .true.
      if (.not. allocated(err)) call transfer%prepare(grid, box, err, options=options)
      call check(error, .not. allocated(err), "prepare failed")
      if (allocated(error)) return
      call transfer%get_report(checked)
      call check(error, checked%method == fitted%method .and. fitted%ninterpolated > 0 .and. &
         & checked%ninterpolated == 1, "check of a covered grid: "//checked%method)
      if (allocated(error)) return
      allocate (f(size(grid%xyz, 2), 1), g(product(box%shape), 1), source=0.0_wp)
      call transfer%apply(f, g, err)
      call check(error, allocated(err), "a transfer prepared as a check was applied")
      if (allocated(error)) return
      ! Coplanar source points near the targets: within reach, no stencil
      plane%kind = "points"
      allocate (plane%xyz(3, 121))
      n = 0
      do j = 1, 11
         do i = 1, 11
            n = n + 1
            plane%xyz(:, n) = [0.5_wp*real(i - 6, wp), 0.5_wp*real(j - 6, wp), 0.0_wp]
         end do
      end do
      above%kind = "points"
      above%xyz = reshape([0.1_wp, 0.2_wp, 0.3_wp, -0.6_wp, 0.4_wp, 0.2_wp], [3, 2])
      call transfer%prepare(plane, above, err, options=options)
      call check(error, .not. allocated(err), "check of a flat source failed")
      if (allocated(error)) return
      call transfer%get_report(checked)
      call check(error, checked%ncopied + checked%ninterpolated == 0 .and. &
         & checked%reason == "source neighbourhood flat or too small", "a flat source passed the check: "//checked%reason)
   end subroutine test_points_check_only

   !> Coinciding source points: a target at a cluster of 20 coinciding points with other
   !> points around takes a finite fitted value; a target at a set of only coinciding points
   !> is not covered; points next to a coinciding pair keep the reach of their neighbourhood
   subroutine test_points_duplicates(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g(:, :)
      logical, allocatable :: covered(:)
      integer :: i, j, k, n

      ! A 5^3 lattice, spacing 0.5, with 20 extra points at its centre
      old%kind = "points"
      allocate (old%xyz(3, 125 + 20))
      n = 0
      do k = 1, 5
         do j = 1, 5
            do i = 1, 5
               n = n + 1
               old%xyz(:, n) = 0.5_wp*real([i, j, k] - 3, wp)
            end do
         end do
      end do
      old%xyz(:, 126:) = 0.0_wp
      allocate (f(145, 1))
      f(:, 1) = 1.0_wp + old%xyz(1, :) - 2.0_wp*old%xyz(3, :)
      new%kind = "points"
      new%xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.1_wp, 0.05_wp, 0.0_wp], [3, 2])
      allocate (g(2, 1), source=0.0_wp)
      call transfer%prepare(old, new, err)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err), "transfer at a cluster failed")
      if (allocated(error)) return
      call check(error, all(covered) .and. report%nduplicate == 1 .and. report%ncopied == 0, &
         & "targets at and next to a cluster: "//report%method)
      if (allocated(error)) return
      call check(error, maxval(abs(g(:, 1) - (1.0_wp + new%xyz(1, :) - 2.0_wp*new%xyz(3, :)))) < 1.0e-10_wp, &
         & "affine field at a cluster")
      if (allocated(error)) return
      ! Only coinciding points
      deallocate (old%xyz, f)
      allocate (old%xyz(3, 30), source=0.0_wp)
      allocate (f(30, 1), source=1.0_wp)
      new%xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1])
      deallocate (g)
      allocate (g(1, 1), source=-1.0_wp)
      call transfer%prepare(old, new, err)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered, report=report)
      call check(error, .not. allocated(err) .and. .not. any(covered) .and. g(1, 1) == -1.0_wp, &
         & "a target at only coinciding points was covered: "//report%reason)
   end subroutine test_points_duplicates

   !> Extrapolation is bounded: above a slab of two layers 0.5 bohr apart, a target 0.25 bohr
   !> above the top layer is extrapolated and reproduces an affine field; one 1 bohr above
   !> needs weights near 3 and -2 and is not covered
   subroutine test_points_extrapolation(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_layout_type) :: old, new
      type(grid_transfer_type) :: transfer
      real(wp), allocatable :: f(:, :), g(:, :)
      logical, allocatable :: covered(:)
      integer :: i, j, l, n

      old%kind = "points"
      allocate (old%xyz(3, 2*121))
      n = 0
      do l = 0, 1
         do j = 1, 11
            do i = 1, 11
               n = n + 1
               old%xyz(:, n) = [0.5_wp*real(i - 6, wp), 0.5_wp*real(j - 6, wp), 0.5_wp*real(l, wp)]
            end do
         end do
      end do
      f = reshape(2.0_wp - old%xyz(1, :) + 3.0_wp*old%xyz(3, :), [n, 1])
      new%kind = "points"
      new%xyz = reshape([0.05_wp, 0.05_wp, 0.75_wp, 0.05_wp, 0.05_wp, 1.5_wp], [3, 2])
      allocate (g(2, 1), source=0.0_wp)
      call transfer%prepare(old, new, err)
      if (.not. allocated(err)) call transfer%apply(f, g, err, covered=covered)
      call check(error, .not. allocated(err), "transfer above a slab failed")
      if (allocated(error)) return
      call check(error, covered(1) .and. .not. covered(2), "extrapolation not bounded as documented")
      if (allocated(error)) return
      call check(error, abs(g(1, 1) - (2.0_wp - 0.05_wp + 3.0_wp*0.75_wp)) < 1.0e-10_wp, "extrapolated affine value")
   end subroutine test_points_extrapolation

end module test_math_grid_transfer
