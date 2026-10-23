!> RISM UV grid settings and construction of Cartesian or molecular point grids
module moist_math_grid_3d_factory
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, default_molecular_recipe, &
      & moist_math_grid_atomic_shell_constant_type, new_constant_shell_policy
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_handymod_type, new_handymod_mapping
   implicit none(type, external)
   private
   public :: rism_grid_parameters_type, new_spatial_grid

   !> Grid kinds, resolved from the case-insensitive name by grid_kind_id
   integer, parameter :: grid_kind_unknown = 0, grid_kind_radial = 1, grid_kind_cartesian = 2, &
      & grid_kind_molecular = 3

   !> UV grid configuration; only keys used by the selected grid are registered
   type, extends(moist_model_parameters_type) :: rism_grid_parameters_type
      !> UV dimension: 1 (rism1d) or 3 (rism3d), default 3
      integer, private :: dimension = 3
      !> Grid kind: radial (1D only), cartesian or molecular, default cartesian
      character(len=16) :: kind = "cartesian"
      !> Grid spacing, bohr, positive, default 0.6 (0.015 for radial)
      real(wp) :: dr = 0.6_wp
      !> Radial npts, Cartesian points per axis; at least 2, defaults 8192 and 128
      integer :: npts = 8192, nx = 128, ny = 128, nz = 128
      !> Cartesian grid: box margin around the nuclei, bohr, nonnegative, default 10
      real(wp) :: margin = 10.0_wp
      !> Molecular grid: radial points (50), angular degree (17), Becke iterations (3)
      integer :: nrad = 50, angular_degree = 17, becke_k = 3
      !> Molecular grid: rmax (bohr, 10), kbuffer (bohr, 2), nufft tolerance (1e-10)
      real(wp) :: rmax = 10.0_wp, kbuffer = 2.0_wp, nufft_tolerance = 1.0e-10_wp
      !> HandyMod radial mapping exponent, dimensionless, positive, default 2
      real(wp) :: handymod_exponent = 2.0_wp
   contains
      procedure :: init_defaults => grid_defaults
      procedure :: register_entries => grid_entries
      procedure :: validate => grid_validate
      procedure :: set_dimension => grid_set_dimension
   end type rism_grid_parameters_type

contains

   !> Construct a Cartesian or bounded molecular point-grid template
   !>
   !> @param[in] params Spatial grid settings
   !> @param[out] grid Geometry-free grid, updated by the model
   !> @param[out] error Native grid or atomic-recipe construction failure
   subroutine new_spatial_grid(params, grid, error)
      !> Spatial grid settings
      type(rism_grid_parameters_type), intent(in) :: params
      !> Geometry-free grid template
      class(moist_math_grid_3d_type), allocatable, intent(out) :: grid
      !> Construction failure, allocated on error
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), allocatable :: cartesian
      type(moist_math_grid_3d_molecular_type), allocatable :: molecular
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_radial_mapping_handymod_type) :: mapping
      type(moist_math_grid_atomic_shell_constant_type) :: shells

      select case (grid_kind_id(params%kind))
      case (grid_kind_cartesian)
         allocate(cartesian)
         call new_cartesian_point_grid(cartesian, error, nx=params%nx, ny=params%ny, nz=params%nz, &
            & dr=params%dr, margin=params%margin)
         if (.not. allocated(error)) call move_alloc(cartesian, grid)
      case (grid_kind_molecular)
         call default_molecular_recipe(recipe, error)
         if (allocated(error)) return
         call new_handymod_mapping(mapping, 0.0_wp, params%rmax, params%handymod_exponent, error)
         if (.not. allocated(error)) call new_constant_shell_policy(shells, params%angular_degree, error)
         if (allocated(error)) return
         recipe%radial%npts = params%nrad
         deallocate(recipe%radial%mapping, recipe%shells)
         allocate(recipe%radial%mapping, source=mapping)
         allocate(recipe%shells, source=shells)
         allocate(molecular)
         call new_molecular_point_grid(molecular, error, recipe=recipe, dr=params%dr, &
            & kbuffer=params%kbuffer, nufft_tol=params%nufft_tolerance, becke_k=params%becke_k)
         if (.not. allocated(error)) call move_alloc(molecular, grid)
      case default
         call fatal_error(error, "Unknown spatial RISM grid: "//trim(params%kind))
      end select
   end subroutine new_spatial_grid

   !> Select the UV dimension before defaults are restored: 1 for rism1d, 3 for rism3d
   !>
   !> @param[inout] self Grid settings
   !> @param[in] dimension UV dimension
   subroutine grid_set_dimension(self, dimension)
      !> Grid settings
      class(rism_grid_parameters_type), intent(inout) :: self
      !> UV dimension: 1 or 3
      integer, intent(in) :: dimension
      self%dimension = dimension
   end subroutine grid_set_dimension

   !> Restore UV grid defaults for the selected dimension
   !>
   !> @param[inout] self Grid settings
   subroutine grid_defaults(self)
      !> Grid settings
      class(rism_grid_parameters_type), intent(inout) :: self
      type(rism_grid_parameters_type) :: defaults
      self%kind = "cartesian"
      self%dr = defaults%dr
      if (self%dimension == 1) then
         self%kind = "radial"
         self%dr = 0.015_wp
      end if
      self%npts = defaults%npts
      self%nx = defaults%nx
      self%ny = defaults%ny
      self%nz = defaults%nz
      self%margin = defaults%margin
      self%nrad = defaults%nrad
      self%angular_degree = defaults%angular_degree
      self%becke_k = defaults%becke_k
      self%rmax = defaults%rmax
      self%handymod_exponent = defaults%handymod_exponent
      self%kbuffer = defaults%kbuffer
      self%nufft_tolerance = defaults%nufft_tolerance
   end subroutine grid_defaults

   !> Register the selected grid's keys; unused grid keys fail strict input
   !>
   !> @param[inout] self Grid settings
   subroutine grid_entries(self)
      !> Grid settings
      class(rism_grid_parameters_type), intent(inout), target :: self
      call self%register_string("kind", self%kind, level=2)
      call self%register_real_scalar("dr", self%dr, level=2)
      select case (grid_kind_id(self%kind))
      case (grid_kind_radial)
         call self%register_int_scalar("npts", self%npts, level=2)
      case (grid_kind_cartesian)
         call self%register_int_scalar("nx", self%nx, level=2)
         call self%register_int_scalar("ny", self%ny, level=2)
         call self%register_int_scalar("nz", self%nz, level=2)
         call self%register_real_scalar("margin", self%margin, level=2)
      case (grid_kind_molecular)
         call self%register_int_scalar("nrad", self%nrad, level=2)
         call self%register_int_scalar("angular_degree", self%angular_degree, level=2)
         call self%register_real_scalar("rmax", self%rmax, level=2)
         call self%register_real_scalar("handymod_exponent", self%handymod_exponent)
         call self%register_int_scalar("becke_k", self%becke_k)
         call self%register_real_scalar("kbuffer", self%kbuffer)
         call self%register_real_scalar("nufft_tolerance", self%nufft_tolerance)
      case default
         ! Unknown grid kinds are rejected by grid_validate
         continue
      end select
   end subroutine grid_entries

   !> Check grid names and dimensions without realizing geometry
   !>
   !> @param[inout] self Grid settings
   !> @param[out] error Unsupported grid or invalid size
   subroutine grid_validate(self, error)
      !> Grid settings, kind normalized to lower case
      class(rism_grid_parameters_type), intent(inout) :: self
      !> Validation failure, allocated on error
      type(error_type), allocatable, intent(out) :: error
      self%kind = to_lower(trim(adjustl(self%kind)))
      if (.not. positive(self%dr)) then
         call fatal_error(error, "UV grid.dr must be positive (bohr)")
         return
      end if
      if (self%dimension == 1 .and. grid_kind_id(self%kind) /= grid_kind_radial) then
         call fatal_error(error, "rism1d needs grid.kind = 'radial'")
         return
      end if
      select case (grid_kind_id(self%kind))
      case (grid_kind_radial)
         if (self%dimension /= 1 .or. self%npts < 2) then
            call fatal_error(error, "Radial UV grids need rism1d and at least two points")
         end if
      case (grid_kind_cartesian)
         if (min(self%nx, self%ny, self%nz) < 2 .or. .not. nonnegative(self%margin)) then
            call fatal_error(error, "Cartesian UV grids need at least two points per axis and a nonnegative margin")
         end if
      case (grid_kind_molecular)
         if (self%nrad < 1 .or. self%angular_degree < 0 .or. self%becke_k < 1 .or. &
            & .not. positive(self%rmax) .or. .not. positive(self%handymod_exponent) &
            & .or. .not. nonnegative(self%kbuffer) &
            & .or. .not. positive(self%nufft_tolerance)) then
            call fatal_error(error, "Invalid molecular UV grid settings")
         end if
      case default
         call fatal_error(error, "Unknown UV grid.kind: "//trim(self%kind))
      end select
   end subroutine grid_validate

   !> Grid kind of a case-insensitive name; grid_kind_unknown for anything else
   !>
   !> @param[in] name Grid kind name
   pure integer function grid_kind_id(name) result(id)
      !> Grid kind name
      character(len=*), intent(in) :: name
      select case (to_lower(trim(adjustl(name))))
      case ("radial")
         id = grid_kind_radial
      case ("cartesian")
         id = grid_kind_cartesian
      case ("molecular")
         id = grid_kind_molecular
      case default
         id = grid_kind_unknown
      end select
   end function grid_kind_id

   !> Finite positive scalar
   !>
   !> @param[in] value Scalar to check
   pure logical function positive(value)
      !> Scalar to check
      real(wp), intent(in) :: value
      positive = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive

   !> Finite nonnegative scalar
   !>
   !> @param[in] value Scalar to check
   pure logical function nonnegative(value)
      !> Scalar to check
      real(wp), intent(in) :: value
      nonnegative = ieee_is_finite(value) .and. value >= 0.0_wp
   end function nonnegative
end module moist_math_grid_3d_factory
