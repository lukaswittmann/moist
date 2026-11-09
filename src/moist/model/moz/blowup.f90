!> A field on one grid moved onto another ("blow-up"), as a basis-set projection moves
!> orbitals between basis sets
!>
!> - The quantity, its units and its columns stay the same; the grids may be given live
!>   (3D or radial) or by their layouts (`get_layout`)
!> - Each point of the new grid takes the first method that applies (see
!>   `moist_math_grid_transfer`): `copy` from a coinciding point of the old grid, else,
!>   from a Cartesian grid, `fourier` onto Cartesian grids and `trilinear` onto others,
!>   `points` (moving least squares) from atom-centered and other point grids, or
!>   `spline` on radial nodes; the points no method covers keep their values in
!>   `field_new`
!> - For the correlations of a solve, move the indirect correlation t and rebuild c and
!>   h by the closure on the new grid
!> - A restart from a checkpoint also fits the solute rigidly and moves the long-range
!>   tail separately: `moz_checkpoint_guess` in `moist_model_moz_checkpoint`
!>
!>    call moz_blowup(grid_old, t_old, grid_new, t_new, error, covered=covered)
module moist_model_moz_blowup
   use mctc_env, only: wp, error_type
   use moist_math_grid_layout, only: moist_math_grid_layout_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_transfer, only: grid_transfer_type, grid_transfer_options_type, grid_transfer_report_type, &
      & rigid_fit_type
   implicit none(type, external)
   private

   public :: moz_blowup

   !> Move a field between two grids, given live or by their layouts
   interface moz_blowup
      module procedure blowup_grid_3d
      module procedure blowup_grid_radial
      module procedure blowup_layout
   end interface moz_blowup

contains

   !> Move a field between two 3D grids
   !>
   !> @param[in]    grid_old   updated grid of the field
   !> @param[in]    field_old  field on the old grid (npoint_old, ncol)
   !> @param[in]    grid_new   updated grid to move the field to
   !> @param[inout] field_new  field on the new grid (npoint_new, ncol); covered rows are
   !>                          overwritten
   !> @param[out]   error      grids not described, an unknown Cartesian method or `check_only`
   !>                          in `options`, shapes that do not match, a failed interpolation,
   !>                          or a point not covered without `covered`
   !> @param[out]   covered    points of the new grid with a value
   !> @param[in]    map        rigid map of the new points into the frame of the old grid
   !> @param[in]    options    methods and threads
   !> @param[out]   report     coverage and methods
   subroutine blowup_grid_3d(grid_old, field_old, grid_new, field_new, error, covered, map, options, report)
      !> Old grid
      class(moist_math_grid_3d_type), intent(in) :: grid_old
      !> Field on the old grid
      real(wp), intent(in) :: field_old(:, :)
      !> New grid
      class(moist_math_grid_3d_type), intent(in) :: grid_new
      !> Field on the new grid
      real(wp), intent(inout) :: field_new(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Points with a value
      logical, allocatable, intent(out), optional :: covered(:)
      !> Map into the old frame
      type(rigid_fit_type), intent(in), optional :: map
      !> Methods and threads
      type(grid_transfer_options_type), intent(in), optional :: options
      !> Coverage and methods
      type(grid_transfer_report_type), intent(out), optional :: report
      type(moist_math_grid_layout_type) :: old, new

      call grid_old%get_layout(old)
      call grid_new%get_layout(new)
      call blowup_layout(old, field_old, new, field_new, error, covered, map, options, report)
   end subroutine blowup_grid_3d

   !> Move a field between two radial grids
   !>
   !> @param[in]    grid_old   built grid of the field
   !> @param[in]    field_old  field on the old nodes (npts_old, ncol)
   !> @param[in]    grid_new   built grid to move the field to
   !> @param[inout] field_new  field on the new nodes (npts_new, ncol); covered rows are
   !>                          overwritten
   !> @param[out]   error      grids not described, an unknown Cartesian method or `check_only`
   !>                          in `options`, shapes that do not match, or a node not covered
   !>                          without `covered`
   !> @param[out]   covered    nodes of the new grid with a value
   !> @param[in]    options    methods and threads
   !> @param[out]   report     coverage and methods
   subroutine blowup_grid_radial(grid_old, field_old, grid_new, field_new, error, covered, options, report)
      !> Old grid
      class(moist_math_grid_radial_type), intent(in) :: grid_old
      !> Field on the old nodes
      real(wp), intent(in) :: field_old(:, :)
      !> New grid
      class(moist_math_grid_radial_type), intent(in) :: grid_new
      !> Field on the new nodes
      real(wp), intent(inout) :: field_new(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Nodes with a value
      logical, allocatable, intent(out), optional :: covered(:)
      !> Methods and threads
      type(grid_transfer_options_type), intent(in), optional :: options
      !> Coverage and methods
      type(grid_transfer_report_type), intent(out), optional :: report
      type(moist_math_grid_layout_type) :: old, new

      call grid_old%get_layout(old)
      call grid_new%get_layout(new)
      call blowup_layout(old, field_old, new, field_new, error, covered, options=options, report=report)
   end subroutine blowup_grid_radial

   !> Move a field between two grids given by their layouts
   !>
   !> @param[in]    old        layout of the old grid
   !> @param[in]    field_old  field on the old grid (npoint_old, ncol)
   !> @param[in]    new        layout of the new grid
   !> @param[inout] field_new  field on the new grid (npoint_new, ncol); covered rows are
   !>                          overwritten
   !> @param[out]   error      layouts not described, radial and 3D mixed, an unknown Cartesian
   !>                          method or `check_only` in `options`, shapes that do not match, a
   !>                          failed interpolation, or a point not covered without `covered`
   !> @param[out]   covered    points of the new grid with a value
   !> @param[in]    map        rigid map of the new points into the frame of the old grid
   !> @param[in]    options    methods and threads
   !> @param[out]   report     coverage and methods
   subroutine blowup_layout(old, field_old, new, field_new, error, covered, map, options, report)
      !> Old layout
      type(moist_math_grid_layout_type), intent(in) :: old
      !> Field on the old grid
      real(wp), intent(in) :: field_old(:, :)
      !> New layout
      type(moist_math_grid_layout_type), intent(in) :: new
      !> Field on the new grid
      real(wp), intent(inout) :: field_new(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Points with a value
      logical, allocatable, intent(out), optional :: covered(:)
      !> Map into the old frame
      type(rigid_fit_type), intent(in), optional :: map
      !> Methods and threads
      type(grid_transfer_options_type), intent(in), optional :: options
      !> Coverage and methods
      type(grid_transfer_report_type), intent(out), optional :: report
      type(grid_transfer_type) :: transfer

      call transfer%prepare(old, new, error, options=options, map=map)
      if (allocated(error)) return
      call transfer%apply(field_old, field_new, error, covered=covered, report=report)
   end subroutine blowup_layout

end module moist_model_moz_blowup
