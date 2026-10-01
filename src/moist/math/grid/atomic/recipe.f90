!> Atomic grid recipe, shell policies, and per-element overrides
!>
!> - Recipe: radial recipe, angular generator, and shell policy
!> - Shell policy: turns the element and a shell radius into an angular
!>   request; the angular generator then chooses the rule. The per-shell
!>   `min_points` and `max_points` are hard bounds copied into every request
!>   - Constant degree: the same minimum angular exactness degree on every shell
!>   - Radial sector: the minimum angular degree of the sector containing the
!>     shell; sectors are bounded by ascending absolute radii
!>   - Target arc spacing: piecewise-constant target spacing h over radial
!>     bands; soft target of `4*pi*r**2/h**2` angular points, an average
!>     area-based spacing estimate rather than a bound on the largest gap
!>     between points; `new_arc_shell_policy` sets a hard floor and cap, 110
!>     and 5810 points unless given
!> - Overrides: element lists, each paired with its own recipe; the first
!>   override listing an element wins, otherwise the default recipe applies
!> - `default_molecular_recipe`: the bounded recipe a molecular grid uses
!>   when none is given; it fits a reciprocal grid
!> - `default_element_recipes`: the per-element table for integration
!>   accuracy; its mapping is unbounded
module moist_math_grid_atomic_recipe
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_chebyshev2_type, &
      & new_chebyshev2_rule, moist_math_grid_radial_rule_midpoint_type, new_midpoint_rule
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_becke_type, &
      & new_becke_mapping, moist_math_grid_radial_mapping_handymod_type, new_handymod_mapping
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_recipe_type
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_request_type, &
      & moist_math_grid_angular_generator_type, moist_math_grid_angular_generator_lebedev_type, &
      & new_lebedev_generator
   implicit none(type, external)
   private

   public :: moist_math_grid_atomic_shell_type
   public :: check_shell_bounds
   public :: check_sector_edges
   public :: find_sector
   public :: moist_math_grid_atomic_shell_constant_type
   public :: new_constant_shell_policy
   public :: moist_math_grid_atomic_shell_sector_type
   public :: new_sector_shell_policy
   public :: moist_math_grid_atomic_shell_arc_type
   public :: new_arc_shell_policy
   public :: default_arc_min_points
   public :: default_arc_max_points
   public :: moist_math_grid_atomic_recipe_type
   public :: moist_math_grid_atomic_recipe_override_type
   public :: element_override_index
   public :: get_element_recipe
   public :: default_molecular_recipe
   public :: default_element_recipes

   !> Abstract shell policy: (z, r) -> angular request
   type, abstract :: moist_math_grid_atomic_shell_type
      !> Hard per-shell floor on the angular point count
      integer :: min_points = 0
      !> Hard per-shell cap on the angular point count
      integer :: max_points = huge(0)
   contains
      !> Angular request for a shell of radius r around element z
      procedure(shell_request_i), deferred :: request
      !> Check the policy parameters; subtypes extend it with their own data
      procedure :: validate => shell_validate
   end type moist_math_grid_atomic_shell_type

   !> Deferred request shared by every shell policy
   abstract interface
      !> Angular request for a shell of radius r around element z
      !>
      !> @param[in] self  Shell policy
      !> @param[in] z     Atomic number
      !> @param[in] r     Shell radius (bohr)
      pure function shell_request_i(self, z, r) result(request)
         import :: moist_math_grid_atomic_shell_type, moist_math_grid_angular_request_type, wp
         implicit none(type, external)
         !> Shell policy
         class(moist_math_grid_atomic_shell_type), intent(in) :: self
         !> Atomic number
         integer, intent(in) :: z
         !> Shell radius (bohr)
         real(wp), intent(in) :: r
         !> Angular request carrying the policy's hard bounds
         type(moist_math_grid_angular_request_type) :: request
      end function shell_request_i
   end interface

   !> Constant-degree shell policy
   type, extends(moist_math_grid_atomic_shell_type) :: moist_math_grid_atomic_shell_constant_type
      !> Minimum algebraic exactness degree on every shell
      integer :: degree = 0
   contains
      !> Angular request for a shell of radius r around element z
      procedure :: request => constant_request
      !> Check the degree and the hard bounds
      procedure :: validate => constant_validate
   end type moist_math_grid_atomic_shell_constant_type

   !> Radial-sector shell policy
   !>
   !> Sector i covers (edges(i-1), edges(i)]; a shell exactly on an edge
   !> belongs to the inner sector
   type, extends(moist_math_grid_atomic_shell_type) :: moist_math_grid_atomic_shell_sector_type
      !> Ascending sector edges (bohr), shape (nedge)
      real(wp), allocatable :: edges(:)
      !> Minimum exactness degree per sector, shape (nedge + 1)
      integer, allocatable :: degrees(:)
   contains
      !> Angular request for a shell of radius r around element z
      procedure :: request => sector_request
      !> Check the sectors and the hard bounds
      procedure :: validate => sector_validate
   end type moist_math_grid_atomic_shell_sector_type

   !> Default per-shell floor of the arc policy
   !>
   !> Inner shells, where the spacing target asks for fewer points, are a
   !> small part of the grid
   integer, parameter :: default_arc_min_points = 110

   !> Default per-shell cap of the arc policy, the largest Lebedev rule
   integer, parameter :: default_arc_max_points = 5810

   !> Target arc-spacing shell policy
   !>
   !> Band i covers (edges(i-1), edges(i)]; a shell exactly on an edge
   !> belongs to the inner band
   type, extends(moist_math_grid_atomic_shell_type) :: moist_math_grid_atomic_shell_arc_type
      !> Ascending band edges (bohr), shape (nedge)
      real(wp), allocatable :: edges(:)
      !> Target arc spacing per band (bohr), shape (nedge + 1)
      real(wp), allocatable :: spacing(:)
   contains
      !> Angular request for a shell of radius r around element z
      procedure :: request => arc_request
      !> Check the bands and the hard bounds
      procedure :: validate => arc_validate
   end type moist_math_grid_atomic_shell_arc_type

   !> Atomic grid recipe: radial recipe, angular generator, and shell policy
   type :: moist_math_grid_atomic_recipe_type
      !> Radial recipe (rule, node count, mapping, optional cutoffs)
      type(moist_math_grid_radial_recipe_type) :: radial
      !> Angular generator; used exactly as given
      class(moist_math_grid_angular_generator_type), allocatable :: angular
      !> Shell policy turning (z, r) into an angular request
      class(moist_math_grid_atomic_shell_type), allocatable :: shells
   end type moist_math_grid_atomic_recipe_type

   !> Per-element override: the elements it applies to and their recipe
   type :: moist_math_grid_atomic_recipe_override_type
      !> Atomic numbers the override applies to, shape (nelem)
      integer, allocatable :: elements(:)
      !> Recipe used for those elements
      type(moist_math_grid_atomic_recipe_type) :: recipe
   end type moist_math_grid_atomic_recipe_override_type

contains

   !> Check the hard per-shell bounds
   !>
   !> @param[in]  self   Shell policy
   !> @param[out] error  Set for a negative floor or a cap below the floor
   subroutine shell_validate(self, error)
      !> Shell policy
      class(moist_math_grid_atomic_shell_type), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_shell_bounds(self%min_points, self%max_points, "Shell policy", error)
   end subroutine shell_validate

   !> Check a per-shell floor and cap
   !>
   !> @param[in]  min_points  Floor on the angular point count, >= 0
   !> @param[in]  max_points  Cap on the angular point count, >= min_points
   !> @param[in]  label       Error-message label
   !> @param[out] error       Set for a negative floor or a cap below the floor
   subroutine check_shell_bounds(min_points, max_points, label, error)
      !> Floor on the angular point count
      integer, intent(in) :: min_points
      !> Cap on the angular point count
      integer, intent(in) :: max_points
      !> Error-message label
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (min_points < 0) then
         call fatal_error(error, label//": min_points must be non-negative")
         return
      end if
      if (max_points < min_points) then
         call fatal_error(error, label//": max_points must not be below min_points")
         return
      end if
   end subroutine check_shell_bounds

   !> Check ascending sector edges against the per-sector values
   !>
   !> @param[in]  edges    Sector edges (bohr); finite, positive, strictly ascending
   !> @param[in]  nvalues  Number of per-sector values; must be size(edges) + 1
   !> @param[in]  label    Error-message label
   !> @param[out] error    Set when the edges do not define size(edges) + 1 sectors
   subroutine check_sector_edges(edges, nvalues, label, error)
      !> Sector edges (bohr)
      real(wp), intent(in) :: edges(:)
      !> Number of per-sector values
      integer, intent(in) :: nvalues
      !> Error-message label
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i

      if (nvalues /= size(edges) + 1) then
         call fatal_error(error, label//": needs exactly one more sector value than edges")
         return
      end if
      if (.not. all(ieee_is_finite(edges))) then
         call fatal_error(error, label//": edges must be finite")
         return
      end if
      if (any(edges <= 0.0_wp)) then
         call fatal_error(error, label//": edges must be positive")
         return
      end if
      do i = 2, size(edges)
         if (edges(i) <= edges(i - 1)) then
            call fatal_error(error, label//": edges must be strictly ascending")
            return
         end if
      end do
   end subroutine check_sector_edges

   !> Sector of radius r among ascending edges
   !>
   !> Sector i covers (edges(i-1), edges(i)]; a radius exactly on an edge
   !> belongs to the inner sector, and radii beyond the last edge to sector
   !> size(edges) + 1
   !>
   !> @param[in] edges  Ascending sector edges (bohr)
   !> @param[in] r      Shell radius (bohr)
   pure function find_sector(edges, r) result(isec)
      !> Ascending sector edges (bohr)
      real(wp), intent(in) :: edges(:)
      !> Shell radius (bohr)
      real(wp), intent(in) :: r
      !> Sector index, 1 to size(edges) + 1
      integer :: isec

      integer :: i

      isec = size(edges) + 1
      do i = 1, size(edges)
         if (r <= edges(i)) then
            isec = i
            return
         end if
      end do
   end function find_sector

   !> Create a constant-degree shell policy
   !>
   !> @param[out] self        Shell policy
   !> @param[in]  degree      Minimum exactness degree on every shell, >= 0
   !> @param[out] error       Set for a negative degree or invalid bounds
   !> @param[in]  min_points  Hard per-shell floor, default 0
   !> @param[in]  max_points  Hard per-shell cap, default huge(0)
   subroutine new_constant_shell_policy(self, degree, error, min_points, max_points)
      !> Shell policy
      type(moist_math_grid_atomic_shell_constant_type), intent(out) :: self
      !> Minimum exactness degree
      integer, intent(in) :: degree
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Hard per-shell floor
      integer, intent(in), optional :: min_points
      !> Hard per-shell cap
      integer, intent(in), optional :: max_points

      self%degree = degree
      if (present(min_points)) self%min_points = min_points
      if (present(max_points)) self%max_points = max_points
      call self%validate(error)
   end subroutine new_constant_shell_policy

   !> Check the degree and the hard bounds
   !>
   !> @param[in]  self   Shell policy
   !> @param[out] error  Set for a negative degree or invalid bounds
   subroutine constant_validate(self, error)
      !> Shell policy
      class(moist_math_grid_atomic_shell_constant_type), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_shell_bounds(self%min_points, self%max_points, "Constant shell policy", error)
      if (allocated(error)) return
      if (self%degree < 0) then
         call fatal_error(error, "Constant shell policy: degree must be non-negative")
         return
      end if
   end subroutine constant_validate

   !> Request `min_degree = degree` with the hard bounds, independent of r
   !>
   !> @param[in] self  Shell policy
   !> @param[in] z     Atomic number (unused)
   !> @param[in] r     Shell radius in bohr (unused)
   pure function constant_request(self, z, r) result(request)
      !> Shell policy
      class(moist_math_grid_atomic_shell_constant_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Shell radius (bohr)
      real(wp), intent(in) :: r
      !> Angular request
      type(moist_math_grid_angular_request_type) :: request

      request%min_degree = self%degree
      request%min_points = self%min_points
      request%max_points = self%max_points
   end function constant_request

   !> Create a radial-sector shell policy
   !>
   !> @param[out] self        Shell policy
   !> @param[in]  edges       Ascending sector edges (bohr), shape (nedge); may be empty
   !> @param[in]  degrees     Minimum exactness degree per sector, shape (nedge + 1), >= 0
   !> @param[out] error       Set for invalid sectors or bounds
   !> @param[in]  min_points  Hard per-shell floor, default 0
   !> @param[in]  max_points  Hard per-shell cap, default huge(0)
   subroutine new_sector_shell_policy(self, edges, degrees, error, min_points, max_points)
      !> Shell policy
      type(moist_math_grid_atomic_shell_sector_type), intent(out) :: self
      !> Ascending sector edges (bohr)
      real(wp), intent(in) :: edges(:)
      !> Minimum exactness degree per sector
      integer, intent(in) :: degrees(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Hard per-shell floor
      integer, intent(in), optional :: min_points
      !> Hard per-shell cap
      integer, intent(in), optional :: max_points

      self%edges = edges
      self%degrees = degrees
      if (present(min_points)) self%min_points = min_points
      if (present(max_points)) self%max_points = max_points
      call self%validate(error)
   end subroutine new_sector_shell_policy

   !> Check the sectors and the hard bounds
   !>
   !> @param[in]  self   Shell policy
   !> @param[out] error  Set for invalid sectors, a negative degree, or invalid bounds
   subroutine sector_validate(self, error)
      !> Shell policy
      class(moist_math_grid_atomic_shell_sector_type), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_shell_bounds(self%min_points, self%max_points, "Sector shell policy", error)
      if (allocated(error)) return
      if (.not. (allocated(self%edges) .and. allocated(self%degrees))) then
         call fatal_error(error, "Sector shell policy: edges and degrees must be set")
         return
      end if
      call check_sector_edges(self%edges, size(self%degrees), "Sector shell policy", error)
      if (allocated(error)) return
      if (any(self%degrees < 0)) then
         call fatal_error(error, "Sector shell policy: degrees must be non-negative")
         return
      end if
   end subroutine sector_validate

   !> Request the degree of the sector containing r, with the hard bounds
   !>
   !> @param[in] self  Shell policy
   !> @param[in] z     Atomic number (unused)
   !> @param[in] r     Shell radius (bohr)
   pure function sector_request(self, z, r) result(request)
      !> Shell policy
      class(moist_math_grid_atomic_shell_sector_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Shell radius (bohr)
      real(wp), intent(in) :: r
      !> Angular request
      type(moist_math_grid_angular_request_type) :: request

      request%min_degree = self%degrees(find_sector(self%edges, r))
      request%min_points = self%min_points
      request%max_points = self%max_points
   end function sector_request

   !> Create a target arc-spacing shell policy
   !>
   !> @param[out] self        Shell policy
   !> @param[in]  edges       Ascending band edges (bohr), shape (nedge); may be empty
   !> @param[in]  spacing     Target arc spacing per band (bohr), shape (nedge + 1), > 0
   !> @param[out] error       Set for invalid bands or bounds
   !> @param[in]  min_points  Hard per-shell floor, default 110
   !> @param[in]  max_points  Hard per-shell cap, default 5810
   subroutine new_arc_shell_policy(self, edges, spacing, error, min_points, max_points)
      !> Shell policy
      type(moist_math_grid_atomic_shell_arc_type), intent(out) :: self
      !> Ascending band edges (bohr)
      real(wp), intent(in) :: edges(:)
      !> Target arc spacing per band (bohr)
      real(wp), intent(in) :: spacing(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Hard per-shell floor
      integer, intent(in), optional :: min_points
      !> Hard per-shell cap
      integer, intent(in), optional :: max_points

      self%edges = edges
      self%spacing = spacing
      self%min_points = default_arc_min_points
      self%max_points = default_arc_max_points
      if (present(min_points)) self%min_points = min_points
      if (present(max_points)) self%max_points = max_points
      call self%validate(error)
   end subroutine new_arc_shell_policy

   !> Check the bands and the hard bounds
   !>
   !> Whether an angular rule exists within [min_points, max_points] depends
   !> on the generator and is reported when the grid is built
   !>
   !> @param[in]  self   Shell policy
   !> @param[out] error  Set for invalid bands, a non-positive spacing, or invalid bounds
   subroutine arc_validate(self, error)
      !> Shell policy
      class(moist_math_grid_atomic_shell_arc_type), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. (allocated(self%edges) .and. allocated(self%spacing))) then
         call fatal_error(error, "Arc shell policy: edges and spacing must be set")
         return
      end if
      call check_sector_edges(self%edges, size(self%spacing), "Arc shell policy", error)
      if (allocated(error)) return
      if (.not. all(ieee_is_finite(self%spacing))) then
         call fatal_error(error, "Arc shell policy: spacing must be finite")
         return
      end if
      if (any(self%spacing <= 0.0_wp)) then
         call fatal_error(error, "Arc shell policy: spacing must be positive")
         return
      end if
      call check_shell_bounds(self%min_points, self%max_points, "Arc shell policy", error)
      if (allocated(error)) return
   end subroutine arc_validate

   !> Request the soft target `4*pi*r**2/h**2` of the band containing r
   !>
   !> The target keeps the spelling `4.0_wp*pi*r*r/(h*h)`, evaluated left
   !> to right and never rounded; a rounding difference next to a rule
   !> size would select a different angular rule
   !>
   !> @param[in] self  Shell policy
   !> @param[in] z     Atomic number (unused)
   !> @param[in] r     Shell radius (bohr)
   pure function arc_request(self, z, r) result(request)
      !> Shell policy
      class(moist_math_grid_atomic_shell_arc_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Shell radius (bohr)
      real(wp), intent(in) :: r
      !> Angular request
      type(moist_math_grid_angular_request_type) :: request

      real(wp) :: h

      h = self%spacing(find_sector(self%edges, r))
      request%target_points = 4.0_wp*pi*r*r/(h*h)
      request%min_points = self%min_points
      request%max_points = self%max_points
   end function arc_request

   !> Index of the first override that lists element z
   !>
   !> An unallocated override array passed as actual argument counts as
   !> absent
   !>
   !> @param[in] z          Atomic number
   !> @param[in] overrides  Per-element overrides, shape (noverride)
   pure function element_override_index(z, overrides) result(idx)
      !> Atomic number
      integer, intent(in) :: z
      !> Per-element overrides
      type(moist_math_grid_atomic_recipe_override_type), intent(in), optional :: overrides(:)
      !> Override index, 0 if none lists z (the default recipe applies)
      integer :: idx

      integer :: i

      idx = 0
      if (.not. present(overrides)) return
      do i = 1, size(overrides)
         if (.not. allocated(overrides(i)%elements)) cycle
         if (any(overrides(i)%elements == z)) then
            idx = i
            return
         end if
      end do
   end function element_override_index

   !> Copy of the recipe that applies to element z
   !>
   !> The first override listing z, otherwise the default recipe
   !>
   !> @param[out] selected   Copy of the applicable recipe
   !> @param[in]  recipe     Default recipe
   !> @param[in]  z          Atomic number
   !> @param[in]  overrides  Per-element overrides, shape (noverride); unallocated counts as absent
   subroutine get_element_recipe(selected, recipe, z, overrides)
      !> Copy of the applicable recipe
      type(moist_math_grid_atomic_recipe_type), intent(out) :: selected
      !> Default recipe
      type(moist_math_grid_atomic_recipe_type), intent(in) :: recipe
      !> Atomic number
      integer, intent(in) :: z
      !> Per-element overrides
      type(moist_math_grid_atomic_recipe_override_type), intent(in), optional :: overrides(:)

      integer :: idx

      idx = element_override_index(z, overrides)
      if (idx == 0) then
         selected = recipe
      else
         selected = overrides(idx)%recipe
      end if
   end subroutine get_element_recipe

   !> Recipe a molecular grid uses when none is given
   !>
   !> Midpoint rule, 50 points, HandyMod(rmin = 0, rmax = 10 bohr, m = 2),
   !> Lebedev with positive weights only, constant minimum degree 17
   !> (110 points). Bounded at 10 bohr, so it fits a reciprocal grid
   !>
   !> @param[out] recipe  Default recipe
   !> @param[out] error   Set if a component cannot be built
   subroutine default_molecular_recipe(recipe, error)
      !> Default recipe
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_rule_midpoint_type) :: rule
      type(moist_math_grid_radial_mapping_handymod_type) :: mapping
      type(moist_math_grid_angular_generator_lebedev_type) :: generator
      type(moist_math_grid_atomic_shell_constant_type) :: shells

      call new_midpoint_rule(rule)
      call new_handymod_mapping(mapping, 0.0_wp, 10.0_wp, 2.0_wp, error)
      if (allocated(error)) return
      call new_lebedev_generator(generator, positive_weights_only=.true.)
      call new_constant_shell_policy(shells, 17, error)
      if (allocated(error)) return

      allocate (recipe%radial%rule, source=rule)
      recipe%radial%npts = 50
      allocate (recipe%radial%mapping, source=mapping)
      allocate (recipe%angular, source=generator)
      allocate (recipe%shells, source=shells)
   end subroutine default_molecular_recipe

   !> Per-element default recipes
   !>
   !> Chebyshev-II radial rule with a Becke mapping, a Lebedev generator
   !> with positive weights only, and a constant-degree shell policy. The
   !> mapping is unbounded: a molecular grid with a reciprocal grid needs
   !> an `rcut_upper` on these recipes:
   !>
   !>   H          50 shells, radius_factor 1.0, degree 29 (302 points)
   !>   He         50 shells, radius_factor 0.5, degree 29 (302 points)
   !>   Li..Ne     75 shells, radius_factor 0.5, degree 29 (302 points)
   !>   Na..Ar     75 shells, radius_factor 0.5, degree 35 (434 points)
   !>   K and up   99 shells, radius_factor 0.5, degree 41 (590 points), the default recipe
   !>
   !> @param[out] recipe     Default recipe (K and heavier elements)
   !> @param[out] overrides  Overrides for H, He, Li..Ne, and Na..Ar, shape (4)
   !> @param[out] error      Set if a component cannot be built
   subroutine default_element_recipes(recipe, overrides, error)
      !> Default recipe
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Per-element overrides
      type(moist_math_grid_atomic_recipe_override_type), allocatable, intent(out) :: overrides(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i

      call default_table_recipe(recipe, 99, 0.5_wp, 41, error)
      if (allocated(error)) return

      allocate (overrides(4))
      overrides(1)%elements = [1]
      call default_table_recipe(overrides(1)%recipe, 50, 1.0_wp, 29, error)
      if (allocated(error)) return
      overrides(2)%elements = [2]
      call default_table_recipe(overrides(2)%recipe, 50, 0.5_wp, 29, error)
      if (allocated(error)) return
      overrides(3)%elements = [(i, i=3, 10)]
      call default_table_recipe(overrides(3)%recipe, 75, 0.5_wp, 29, error)
      if (allocated(error)) return
      overrides(4)%elements = [(i, i=11, 18)]
      call default_table_recipe(overrides(4)%recipe, 75, 0.5_wp, 35, error)
      if (allocated(error)) return
   end subroutine default_element_recipes

   !> One row of the per-element default table
   !>
   !> @param[out] recipe         Recipe of that row
   !> @param[in]  nrad           Number of Chebyshev-II shells
   !> @param[in]  radius_factor  Becke scale as a multiple of the covalent radius
   !> @param[in]  degree         Constant minimum Lebedev degree
   !> @param[out] error          Set if a component cannot be built
   subroutine default_table_recipe(recipe, nrad, radius_factor, degree, error)
      !> Recipe of that row
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Number of Chebyshev-II shells
      integer, intent(in) :: nrad
      !> Becke scale as a multiple of the covalent radius
      real(wp), intent(in) :: radius_factor
      !> Constant minimum Lebedev degree
      integer, intent(in) :: degree
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_rule_chebyshev2_type) :: rule
      type(moist_math_grid_radial_mapping_becke_type) :: mapping
      type(moist_math_grid_angular_generator_lebedev_type) :: generator
      type(moist_math_grid_atomic_shell_constant_type) :: shells

      call new_chebyshev2_rule(rule)
      call new_becke_mapping(mapping, error, radius_factor=radius_factor)
      if (allocated(error)) return
      call new_lebedev_generator(generator, positive_weights_only=.true.)
      call new_constant_shell_policy(shells, degree, error)
      if (allocated(error)) return

      allocate (recipe%radial%rule, source=rule)
      recipe%radial%npts = nrad
      allocate (recipe%radial%mapping, source=mapping)
      allocate (recipe%angular, source=generator)
      allocate (recipe%shells, source=shells)
   end subroutine default_table_recipe

end module moist_math_grid_atomic_recipe
