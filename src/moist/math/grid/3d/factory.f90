!> RISM UV grid settings and construction of Cartesian or molecular point grids
module moist_math_grid_3d_factory
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_parameters, only: moist_model_parameters_type, print_hidden
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid, &
      & partition_becke, partition_ssf, partition_pvoronoi
   use moist_math_grid_3d_kernel_becke, only: default_stiffness
   use moist_math_grid_3d_kernel_ssf, only: default_ssf_a
   use moist_math_grid_3d_kernel_pvoronoi, only: default_power_width
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, default_molecular_recipe, &
      & moist_math_grid_atomic_shell_constant_type, new_constant_shell_policy
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_handymod_type, new_handymod_mapping
   implicit none(type, external)
   private
   public :: rism_grid_parameters_type, new_spatial_grid
   public :: rism_radial_grid_parameters_type, rism_cartesian_grid_parameters_type
   public :: rism_molecular_grid_parameters_type
   public :: rism_molecular_radial_parameters_type, rism_molecular_angular_parameters_type

   !> Grid kinds, resolved from the case-insensitive name by grid_kind_id
   integer, parameter :: grid_kind_unknown = 0, grid_kind_radial = 1, grid_kind_cartesian = 2, &
      & grid_kind_molecular = 3

   !> Uniform radial UV grid settings
   type, extends(moist_model_parameters_type) :: rism_radial_grid_parameters_type
      !> Radial point count, at least two
      integer :: npts = 8192
      !> Radial spacing, bohr, positive
      real(wp) :: dr = 0.015_wp
   contains
      procedure :: init_defaults => radial_defaults
      procedure :: register_entries => radial_entries
      procedure :: validate => radial_validate
   end type rism_radial_grid_parameters_type

   !> Cartesian UV grid settings
   type, extends(moist_model_parameters_type) :: rism_cartesian_grid_parameters_type
      !> Points per axis, at least two
      integer :: nx = 128, ny = 128, nz = 128
      !> Grid spacing, bohr, positive
      real(wp) :: dr = 0.6_wp
      !> Box margin around the nuclei, bohr, nonnegative
      real(wp) :: margin = 10.0_wp
   contains
      procedure :: init_defaults => cartesian_defaults
      procedure :: register_entries => cartesian_entries
      procedure :: validate => cartesian_validate
   end type rism_cartesian_grid_parameters_type

   !> Bounded radial shells of an atom-centered molecular grid
   type, extends(moist_model_parameters_type) :: rism_molecular_radial_parameters_type
      !> Radial shell count, at least one
      integer :: npts = 50
      !> HandyMod outer radius, bohr, and positive mapping exponent
      real(wp) :: rmax = 10.0_wp, handymod_exponent = 2.0_wp
   contains
      procedure :: init_defaults => molecular_radial_defaults
      procedure :: register_entries => molecular_radial_entries
      procedure :: validate => molecular_radial_validate
   end type rism_molecular_radial_parameters_type

   !> Angular quadrature on each molecular radial shell
   type, extends(moist_model_parameters_type) :: rism_molecular_angular_parameters_type
      !> Minimum angular exactness degree, nonnegative
      integer :: degree = 17
   contains
      procedure :: init_defaults => molecular_angular_defaults
      procedure :: register_entries => molecular_angular_entries
      procedure :: validate => molecular_angular_validate
   end type rism_molecular_angular_parameters_type

   !> Bounded molecular UV grid and partition settings
   type, extends(moist_model_parameters_type) :: rism_molecular_grid_parameters_type
      !> Radial shell count and placement
      type(rism_molecular_radial_parameters_type) :: radial
      !> Angular resolution of each radial shell
      type(rism_molecular_angular_parameters_type) :: angular
      !> Molecular partition: becke, ssf or pvoronoi
      character(len=16) :: partition = "becke"
      !> Becke polynomial iteration count, at least one
      integer :: becke_k = default_stiffness
      !> SSF switch half-width in (0, 1]
      real(wp) :: ssf_a = default_ssf_a
      !> Power-coordinate switch half-width, bohr**2, positive
      real(wp) :: power_width = default_power_width
      !> Absolute integration-weight cutoff, bohr**3, nonnegative
      real(wp) :: weight_pruning_cutoff = 1.0e-14_wp
      !> Bare owner-partition weight cutoff in [0, 1)
      real(wp) :: partition_pruning_cutoff = 0.0_wp
      !> Reciprocal-grid real-space spacing and period margin, bohr
      real(wp) :: dr = 0.6_wp, kbuffer = 2.0_wp
      !> Requested FINUFFT relative tolerance, positive
      real(wp) :: nufft_tolerance = 1.0e-10_wp
   contains
      procedure :: init_defaults => molecular_defaults
      procedure :: register_entries => molecular_entries
      procedure :: validate => molecular_validate
   end type rism_molecular_grid_parameters_type

   !> UV grid kind and independent settings for each kind
   !>
   !> Files keep every group; printed tables show only the selected grid
   type, extends(moist_model_parameters_type) :: rism_grid_parameters_type
      !> UV dimension: 1 (rism1d) or 3 (rism3d), default 3
      integer, private :: dimension = 3
      !> Grid kind: radial (1D only), cartesian or molecular, default cartesian
      character(len=16) :: kind = "cartesian"
      !> Uniform radial grid settings
      type(rism_radial_grid_parameters_type) :: radial
      !> Cartesian point-grid settings
      type(rism_cartesian_grid_parameters_type) :: cartesian
      !> Atom-centered point-grid settings
      type(rism_molecular_grid_parameters_type) :: molecular
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
      integer, allocatable :: becke_k
      real(wp), allocatable :: ssf_a, power_width
      integer :: scheme

      select case (grid_kind_id(params%kind))
      case (grid_kind_cartesian)
         allocate (cartesian)
         call new_cartesian_point_grid(cartesian, error, nx=params%cartesian%nx, ny=params%cartesian%ny, &
            & nz=params%cartesian%nz, dr=params%cartesian%dr, margin=params%cartesian%margin)
         if (.not. allocated(error)) call move_alloc(cartesian, grid)
      case (grid_kind_molecular)
         call default_molecular_recipe(recipe, error)
         if (allocated(error)) return
         call new_handymod_mapping(mapping, 0.0_wp, params%molecular%radial%rmax, &
            & params%molecular%radial%handymod_exponent, error)
         if (.not. allocated(error)) call new_constant_shell_policy(shells, params%molecular%angular%degree, error)
         if (allocated(error)) return
         recipe%radial%npts = params%molecular%radial%npts
         deallocate (recipe%radial%mapping, recipe%shells)
         allocate (recipe%radial%mapping, source=mapping)
         allocate (recipe%shells, source=shells)
         scheme = partition_id(params%molecular%partition)
         select case (scheme)
         case (partition_becke)
            becke_k = params%molecular%becke_k
         case (partition_ssf)
            ssf_a = params%molecular%ssf_a
         case (partition_pvoronoi)
            power_width = params%molecular%power_width
         case default
            call fatal_error(error, "Unknown grid.molecular.partition: "//trim(params%molecular%partition))
            return
         end select
         allocate (molecular)
         call new_molecular_point_grid(molecular, error, recipe=recipe, dr=params%molecular%dr, &
            & kbuffer=params%molecular%kbuffer, nufft_tol=params%molecular%nufft_tolerance, partition=scheme, &
            & becke_k=becke_k, ssf_a=ssf_a, power_width=power_width, &
            & weight_threshold=params%molecular%weight_pruning_cutoff, pruning_threshold=params%molecular%partition_pruning_cutoff)
         if (.not. allocated(error)) call move_alloc(molecular, grid)
      case default
         call fatal_error(error, "Unknown spatial RISM grid: "//trim(params%kind))
      end select
   end subroutine new_spatial_grid

   !> Select the UV dimension before defaults are restored: 1 for rism1d, 3 for rism3d
   !>
   !> @param[in,out] self Grid settings
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
   !> @param[in,out] self Grid settings
   subroutine grid_defaults(self)
      !> Grid settings
      class(rism_grid_parameters_type), intent(inout) :: self
      self%kind = "cartesian"
      if (self%dimension == 1) self%kind = "radial"
      call self%radial%init_defaults()
      call self%cartesian%init_defaults()
      call self%molecular%init_defaults()
   end subroutine grid_defaults

   !> Register every grid group and print only the selected kind
   !>
   !> @param[in,out] self Grid settings
   subroutine grid_entries(self)
      !> Grid settings
      class(rism_grid_parameters_type), intent(inout), target :: self
      call self%register_string("kind", self%kind, level=2)
      call self%register_group("radial", self%radial, &
         & level=merge(2, print_hidden, grid_kind_id(self%kind) == grid_kind_radial))
      call self%register_group("cartesian", self%cartesian, &
         & level=merge(2, print_hidden, grid_kind_id(self%kind) == grid_kind_cartesian))
      call self%register_group("molecular", self%molecular, &
         & level=merge(2, print_hidden, grid_kind_id(self%kind) == grid_kind_molecular))
   end subroutine grid_entries

   !> Check grid names and dimensions without realizing geometry
   !>
   !> @param[in,out] self Grid settings
   !> @param[out] error Unsupported grid or invalid size
   subroutine grid_validate(self, error)
      !> Grid settings, kind normalized to lower case
      class(rism_grid_parameters_type), intent(inout) :: self
      !> Validation failure, allocated on error
      type(error_type), allocatable, intent(out) :: error
      self%kind = to_lower(trim(adjustl(self%kind)))
      if (self%dimension == 1 .and. grid_kind_id(self%kind) /= grid_kind_radial) then
         call fatal_error(error, "rism1d needs grid.kind = 'radial'")
         return
      end if
      select case (grid_kind_id(self%kind))
      case (grid_kind_radial)
         if (self%dimension /= 1) call fatal_error(error, "Radial UV grids need rism1d")
      case (grid_kind_cartesian, grid_kind_molecular)
         continue
      case default
         call fatal_error(error, "Unknown UV grid.kind: "//trim(self%kind))
      end select
      if (allocated(error)) return
      call self%radial%validate(error)
      if (.not. allocated(error)) call self%cartesian%validate(error)
      if (.not. allocated(error)) call self%molecular%validate(error)
   end subroutine grid_validate

   !> Restore radial grid defaults
   !>
   !> @param[in,out] self Radial settings
   subroutine radial_defaults(self)
      !> Radial settings
      class(rism_radial_grid_parameters_type), intent(inout) :: self
      type(rism_radial_grid_parameters_type) :: defaults
      self%npts = defaults%npts
      self%dr = defaults%dr
   end subroutine radial_defaults

   !> Register radial grid settings
   !>
   !> @param[in,out] self Radial settings
   subroutine radial_entries(self)
      !> Radial settings
      class(rism_radial_grid_parameters_type), intent(inout), target :: self
      call self%register_int_scalar("npts", self%npts, level=2)
      call self%register_real_scalar("dr", self%dr, level=2)
   end subroutine radial_entries

   !> Check radial grid dimensions
   !>
   !> @param[in,out] self Radial settings
   !> @param[out] error Invalid point count or spacing
   subroutine radial_validate(self, error)
      !> Radial settings
      class(rism_radial_grid_parameters_type), intent(inout) :: self
      !> Invalid point count or spacing
      type(error_type), allocatable, intent(out) :: error
      if (self%npts < 2 .or. .not. positive(self%dr)) then
         call fatal_error(error, "grid.radial needs at least two points and positive dr (bohr)")
      end if
   end subroutine radial_validate

   !> Restore Cartesian grid defaults
   !>
   !> @param[in,out] self Cartesian settings
   subroutine cartesian_defaults(self)
      !> Cartesian settings
      class(rism_cartesian_grid_parameters_type), intent(inout) :: self
      type(rism_cartesian_grid_parameters_type) :: defaults
      self%nx = defaults%nx
      self%ny = defaults%ny
      self%nz = defaults%nz
      self%dr = defaults%dr
      self%margin = defaults%margin
   end subroutine cartesian_defaults

   !> Register Cartesian grid settings
   !>
   !> @param[in,out] self Cartesian settings
   subroutine cartesian_entries(self)
      !> Cartesian settings
      class(rism_cartesian_grid_parameters_type), intent(inout), target :: self
      call self%register_int_scalar("nx", self%nx, level=2)
      call self%register_int_scalar("ny", self%ny, level=2)
      call self%register_int_scalar("nz", self%nz, level=2)
      call self%register_real_scalar("dr", self%dr, level=2)
      call self%register_real_scalar("margin", self%margin, level=2)
   end subroutine cartesian_entries

   !> Check Cartesian grid dimensions
   !>
   !> @param[in,out] self Cartesian settings
   !> @param[out] error Invalid point count, spacing or margin
   subroutine cartesian_validate(self, error)
      !> Cartesian settings
      class(rism_cartesian_grid_parameters_type), intent(inout) :: self
      !> Invalid point count, spacing or margin
      type(error_type), allocatable, intent(out) :: error
      if (min(self%nx, self%ny, self%nz) < 2 .or. .not. positive(self%dr) .or. .not. nonnegative(self%margin)) then
         call fatal_error(error, "grid.cartesian needs at least two points per axis, positive dr and nonnegative margin (bohr)")
      end if
   end subroutine cartesian_validate

   !> Restore molecular grid defaults
   !>
   !> @param[in,out] self Molecular settings
   subroutine molecular_defaults(self)
      !> Molecular settings
      class(rism_molecular_grid_parameters_type), intent(inout) :: self
      type(rism_molecular_grid_parameters_type) :: defaults
      call self%radial%init_defaults()
      call self%angular%init_defaults()
      self%partition = defaults%partition
      self%becke_k = defaults%becke_k
      self%ssf_a = defaults%ssf_a
      self%power_width = defaults%power_width
      self%weight_pruning_cutoff = defaults%weight_pruning_cutoff
      self%partition_pruning_cutoff = defaults%partition_pruning_cutoff
      self%dr = defaults%dr
      self%kbuffer = defaults%kbuffer
      self%nufft_tolerance = defaults%nufft_tolerance
   end subroutine molecular_defaults

   !> Register molecular grid settings and print the selected partition controls
   !>
   !> @param[in,out] self Molecular settings
   subroutine molecular_entries(self)
      !> Molecular settings
      class(rism_molecular_grid_parameters_type), intent(inout), target :: self
      call self%register_group("radial", self%radial)
      call self%register_group("angular", self%angular)
      call self%register_string("partition", self%partition, level=2)
      call self%register_int_scalar("becke_k", self%becke_k, &
         & level=merge(3, print_hidden, partition_id(self%partition) == partition_becke))
      call self%register_real_scalar("ssf_a", self%ssf_a, &
         & level=merge(3, print_hidden, partition_id(self%partition) == partition_ssf))
      call self%register_real_scalar("power_width", self%power_width, &
         & level=merge(3, print_hidden, partition_id(self%partition) == partition_pvoronoi))
      call self%register_real_scalar("weight_pruning_cutoff", self%weight_pruning_cutoff)
      call self%register_real_scalar("partition_pruning_cutoff", self%partition_pruning_cutoff)
      call self%register_real_scalar("dr", self%dr, level=2)
      call self%register_real_scalar("kbuffer", self%kbuffer)
      call self%register_real_scalar("nufft_tolerance", self%nufft_tolerance)
   end subroutine molecular_entries

   !> Check molecular resolution, partition and pruning settings
   !>
   !> @param[in,out] self Molecular settings, partition name normalized to lower case
   !> @param[out] error Unknown partition or invalid numerical setting
   subroutine molecular_validate(self, error)
      !> Molecular settings
      class(rism_molecular_grid_parameters_type), intent(inout) :: self
      !> Unknown partition or invalid numerical setting
      type(error_type), allocatable, intent(out) :: error
      self%partition = to_lower(trim(adjustl(self%partition)))
      if (partition_id(self%partition) == 0) then
         call fatal_error(error, "Unknown grid.molecular.partition: "//trim(self%partition)//" (choose becke, ssf or pvoronoi)")
      else if (.not. positive(self%dr) .or. &
         & .not. nonnegative(self%kbuffer) .or. .not. positive(self%nufft_tolerance)) then
         call fatal_error(error, "Invalid grid.molecular reciprocal-grid settings")
      else if (self%becke_k < 1 .or. .not. positive(self%ssf_a) .or. self%ssf_a > 1.0_wp .or. &
         & .not. positive(self%power_width)) then
         call fatal_error(error, "grid.molecular needs becke_k >= 1, ssf_a in (0, 1] and positive power_width (bohr**2)")
      else if (.not. nonnegative(self%weight_pruning_cutoff) .or. .not. nonnegative(self%partition_pruning_cutoff) &
         & .or. self%partition_pruning_cutoff >= 1.0_wp) then
         call fatal_error(error, "grid.molecular needs weight_pruning_cutoff >= 0 (bohr**3) "// &
            & "and partition_pruning_cutoff in [0, 1)")
      end if
      if (allocated(error)) return
      call self%radial%validate(error)
      if (.not. allocated(error)) call self%angular%validate(error)
   end subroutine molecular_validate

   !> Restore molecular radial defaults
   !>
   !> @param[in,out] self Molecular radial settings
   subroutine molecular_radial_defaults(self)
      !> Molecular radial settings
      class(rism_molecular_radial_parameters_type), intent(inout) :: self
      type(rism_molecular_radial_parameters_type) :: defaults
      self%npts = defaults%npts
      self%rmax = defaults%rmax
      self%handymod_exponent = defaults%handymod_exponent
   end subroutine molecular_radial_defaults

   !> Register molecular radial settings
   !>
   !> @param[in,out] self Molecular radial settings
   subroutine molecular_radial_entries(self)
      !> Molecular radial settings
      class(rism_molecular_radial_parameters_type), intent(inout), target :: self
      call self%register_int_scalar("npts", self%npts, level=2)
      call self%register_real_scalar("rmax", self%rmax, level=2)
      call self%register_real_scalar("handymod_exponent", self%handymod_exponent)
   end subroutine molecular_radial_entries

   !> Check molecular radial count and mapping settings
   !>
   !> @param[in,out] self Molecular radial settings
   !> @param[out] error Invalid count, radius or exponent
   subroutine molecular_radial_validate(self, error)
      !> Molecular radial settings
      class(rism_molecular_radial_parameters_type), intent(inout) :: self
      !> Invalid count, radius or exponent
      type(error_type), allocatable, intent(out) :: error
      if (self%npts < 1 .or. .not. positive(self%rmax) .or. .not. positive(self%handymod_exponent)) then
         call fatal_error(error, "grid.molecular.radial needs npts >= 1, positive rmax (bohr) and positive handymod_exponent")
      end if
   end subroutine molecular_radial_validate

   !> Restore molecular angular defaults
   !>
   !> @param[in,out] self Molecular angular settings
   subroutine molecular_angular_defaults(self)
      !> Molecular angular settings
      class(rism_molecular_angular_parameters_type), intent(inout) :: self
      type(rism_molecular_angular_parameters_type) :: defaults
      self%degree = defaults%degree
   end subroutine molecular_angular_defaults

   !> Register molecular angular settings
   !>
   !> @param[in,out] self Molecular angular settings
   subroutine molecular_angular_entries(self)
      !> Molecular angular settings
      class(rism_molecular_angular_parameters_type), intent(inout), target :: self
      call self%register_int_scalar("degree", self%degree, level=2)
   end subroutine molecular_angular_entries

   !> Check molecular angular resolution
   !>
   !> @param[in,out] self Molecular angular settings
   !> @param[out] error Negative exactness degree
   subroutine molecular_angular_validate(self, error)
      !> Molecular angular settings
      class(rism_molecular_angular_parameters_type), intent(inout) :: self
      !> Negative exactness degree
      type(error_type), allocatable, intent(out) :: error
      if (self%degree < 0) call fatal_error(error, "grid.molecular.angular.degree must be nonnegative")
   end subroutine molecular_angular_validate

   !> Partition constant of a case-insensitive name; zero for anything else
   !>
   !> @param[in] name Partition name
   pure function partition_id(name) result(id)
      !> Partition name
      character(len=*), intent(in) :: name
      !> Partition constant or zero
      integer :: id
      select case (to_lower(trim(adjustl(name))))
      case ("becke")
         id = partition_becke
      case ("ssf")
         id = partition_ssf
      case ("pvoronoi")
         id = partition_pvoronoi
      case default
         id = 0
      end select
   end function partition_id

   !> Grid kind of a case-insensitive name; grid_kind_unknown for anything else
   !>
   !> @param[in] name Grid kind name
   pure function grid_kind_id(name) result(id)
      !> Grid kind name
      character(len=*), intent(in) :: name
      !> Grid kind constant or grid_kind_unknown
      integer :: id
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
   pure function positive(value) result(r)
      !> Scalar to check
      real(wp), intent(in) :: value
      !> Whether the value is finite and positive
      logical :: r
      r = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive

   !> Finite nonnegative scalar
   !>
   !> @param[in] value Scalar to check
   pure function nonnegative(value) result(r)
      !> Scalar to check
      real(wp), intent(in) :: value
      !> Whether the value is finite and nonnegative
      logical :: r
      r = ieee_is_finite(value) .and. value >= 0.0_wp
   end function nonnegative
end module moist_math_grid_3d_factory
