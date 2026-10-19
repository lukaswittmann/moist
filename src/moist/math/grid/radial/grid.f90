!> Concrete radial grid, its recipe, and the tagged uniform radial pair
!>
!> A radial grid owns nodes, plain `dr` weights (no r^2, no 4*pi), the
!> transform tag that selects its radial transform, and the requested and
!> retained node counts
module moist_math_grid_radial_grid
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_type
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_type
   implicit none(type, external)
   private

   public :: moist_math_grid_radial_type
   public :: moist_math_grid_radial_recipe_type
   public :: new_radial_grid
   public :: new_uniform_radial_pair
   public :: integrand_radial
   public :: transform_quadrature
   public :: transform_dst4

   !> Transform tag: general quadrature (dense Fourier-Bessel kernel)
   integer, parameter :: transform_quadrature = 1
   !> Transform tag: DST-IV on a uniform pair with recorded spacing
   integer, parameter :: transform_dst4 = 2

   !> 4*pi, the solid angle applied by the volume integrators
   real(wp), parameter :: four_pi = 4.0_wp*pi

   !> Scalar-valued radial integrand f(r) used by `integrate`
   abstract interface
      !> Radial integrand f(r)
      !>
      !> @param[in] r  Radius (bohr)
      pure function integrand_radial(r) result(val)
         import :: wp
         implicit none(type, external)
         !> Radius (bohr)
         real(wp), intent(in) :: r
         !> Function value at r
         real(wp) :: val
      end function integrand_radial
   end interface

   !> Radial recipe: base rule, node count, mapping, and optional cutoffs
   type :: moist_math_grid_radial_recipe_type
      !> Reference rule on [-1, 1]
      class(moist_math_grid_radial_rule_type), allocatable :: rule
      !> Number of reference nodes requested from the rule
      integer :: npts = 0
      !> Mapping from [-1, 1] to radii
      class(moist_math_grid_radial_mapping_type), allocatable :: mapping
      !> Lower cutoff in bohr; nodes with r < rcut_lower are dropped (unallocated: none)
      real(wp), allocatable :: rcut_lower
      !> Upper cutoff in bohr; nodes with r > rcut_upper are dropped (unallocated: none)
      real(wp), allocatable :: rcut_upper
   end type moist_math_grid_radial_recipe_type

   !> Concrete radial grid: nodes, dr weights, and transform tag
   type :: moist_math_grid_radial_type
      !> Number of retained nodes
      integer :: npts = 0
      !> Number of nodes requested before cutoffs
      integer :: npts_requested = 0
      !> Nodes (bohr, or 1/bohr on a reciprocal grid), shape (npts)
      real(wp), allocatable :: r(:)
      !> Weights representing dr (or dk), shape (npts)
      real(wp), allocatable :: w(:)
      !> Transform tag, transform_quadrature or transform_dst4; 0 until built
      integer :: transform = 0
      !> Uniform node spacing; set only on transform_dst4 grids, 0 otherwise
      real(wp) :: spacing = 0.0_wp
   contains
      !> Volume quadrature of a field already tabulated on the nodes
      procedure :: integrate_field => radial_integrate_field
      !> Volume quadrature of an analytic radial integrand sampled at the nodes
      procedure :: integrate => radial_integrate
   end type moist_math_grid_radial_type

contains

   !> Build a radial grid for element z from a recipe
   !>
   !> Generates the rule on [-1, 1], applies the mapping, then drops nodes with
   !> r < rcut_lower or r > rcut_upper for each cutoff that is set; node order
   !> is the rule's order
   !>
   !> @param[out] grid    New radial grid, tagged transform_quadrature
   !> @param[in]  recipe  Radial recipe
   !> @param[in]  z       Atomic number, passed to the mapping
   !> @param[out] error   Set on an invalid recipe, a rule or mapping failure,
   !>                     or when the cutoffs remove every node
   subroutine new_radial_grid(grid, recipe, z, error)
      !> New radial grid
      type(moist_math_grid_radial_type), intent(out) :: grid
      !> Radial recipe
      type(moist_math_grid_radial_recipe_type), intent(in) :: recipe
      !> Atomic number
      integer, intent(in) :: z
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: x(:), wx(:), r(:), wr(:)
      logical, allocatable :: keep(:)
      integer :: nkeep

      if (.not. allocated(recipe%rule)) then
         call fatal_error(error, "Radial grid: recipe has no rule")
         return
      end if
      if (.not. allocated(recipe%mapping)) then
         call fatal_error(error, "Radial grid: recipe has no mapping")
         return
      end if
      if (recipe%npts < 1) then
         call fatal_error(error, "Radial grid: recipe needs npts >= 1")
         return
      end if
      if (allocated(recipe%rcut_lower)) then
         if (.not. ieee_is_finite(recipe%rcut_lower)) then
            call fatal_error(error, "Radial grid: rcut_lower must be finite")
            return
         end if
      end if
      if (allocated(recipe%rcut_upper)) then
         if (.not. ieee_is_finite(recipe%rcut_upper)) then
            call fatal_error(error, "Radial grid: rcut_upper must be finite")
            return
         end if
      end if
      if (allocated(recipe%rcut_lower) .and. allocated(recipe%rcut_upper)) then
         if (.not. recipe%rcut_lower < recipe%rcut_upper) then
            call fatal_error(error, "Radial grid: rcut_lower must be below rcut_upper")
            return
         end if
      end if

      call recipe%rule%generate(recipe%npts, x, wx, error)
      if (allocated(error)) return
      call recipe%mapping%transform(z, x, wx, r, wr, error)
      if (allocated(error)) return

      allocate (keep(size(r)))
      keep(:) = .true.
      if (allocated(recipe%rcut_lower)) keep = keep .and. .not. (r < recipe%rcut_lower)
      if (allocated(recipe%rcut_upper)) keep = keep .and. .not. (r > recipe%rcut_upper)
      nkeep = count(keep)
      if (nkeep == 0) then
         call fatal_error(error, "Radial grid: cutoffs removed every node")
         return
      end if

      grid%npts_requested = recipe%npts
      grid%npts = nkeep
      grid%r = pack(r, keep)
      grid%w = pack(wr, keep)
      grid%transform = transform_quadrature
      grid%spacing = 0.0_wp
   end subroutine new_radial_grid

   !> Build the uniform r- and k-grid pair for the DST-IV radial transform
   !>
   !> r_i = (i - 0.5)*dr, k_j = (j - 0.5)*dk with dk = acos(-1)/(npts*dr);
   !> weights dr and dk; both grids are tagged transform_dst4 and record their
   !> spacing
   !>
   !> @param[out] rgrid  r-space grid, nodes in bohr
   !> @param[out] kgrid  k-space grid, nodes in 1/bohr
   !> @param[in]  npts   Number of nodes on each grid (npts >= 1)
   !> @param[in]  dr     r-space spacing in bohr, finite, > 0
   !> @param[out] error  Set on invalid parameters or allocation failure
   subroutine new_uniform_radial_pair(rgrid, kgrid, npts, dr, error)
      !> r-space grid
      type(moist_math_grid_radial_type), intent(out) :: rgrid
      !> k-space grid
      type(moist_math_grid_radial_type), intent(out) :: kgrid
      !> Number of nodes on each grid
      integer, intent(in) :: npts
      !> r-space spacing in bohr
      real(wp), intent(in) :: dr
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i, stat
      real(wp) :: pi_dst, dk

      if (npts < 1) then
         call fatal_error(error, "Uniform radial pair needs at least one point")
         return
      end if
      if (.not. ieee_is_finite(dr)) then
         call fatal_error(error, "Uniform radial pair needs a finite spacing")
         return
      end if
      if (dr <= 0.0_wp) then
         call fatal_error(error, "Uniform radial pair needs a positive spacing")
         return
      end if

      ! The DST-IV trafo recomputes dk with this expression and compares with ==
      pi_dst = acos(-1.0_wp)
      dk = pi_dst/(real(npts, wp)*dr)

      allocate (rgrid%r(npts), stat=stat)
      if (stat == 0) allocate (rgrid%w(npts), stat=stat)
      if (stat == 0) allocate (kgrid%r(npts), stat=stat)
      if (stat == 0) allocate (kgrid%w(npts), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate uniform radial pair nodes")
         return
      end if
      do i = 1, npts
         rgrid%r(i) = (real(i, wp) - 0.5_wp)*dr
         kgrid%r(i) = (real(i, wp) - 0.5_wp)*dk
      end do
      rgrid%w(:) = dr
      kgrid%w(:) = dk

      rgrid%npts = npts
      rgrid%npts_requested = npts
      rgrid%transform = transform_dst4
      rgrid%spacing = dr
      kgrid%npts = npts
      kgrid%npts_requested = npts
      kgrid%transform = transform_dst4
      kgrid%spacing = dk
   end subroutine new_uniform_radial_pair

   !> Volume quadrature of tabulated values: sum_i 4*pi*r_i^2*w_i*f_i
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Per-node field values, shape (npts); any other length is an error
   !> @param[out] result  Quadrature result
   !> @param[out] error   Error handling
   subroutine radial_integrate_field(self, f, result, error)
      !> Grid instance
      class(moist_math_grid_radial_type), intent(in) :: self
      !> Per-node field values
      real(wp), intent(in) :: f(:)
      !> Quadrature result
      real(wp), intent(out) :: result
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i

      result = 0.0_wp
      if (size(f) /= self%npts) then
         call fatal_error(error, "radial grid: integrate_field needs one value per node")
         return
      end if
      do i = 1, self%npts
         result = result + four_pi*self%r(i)*self%r(i)*self%w(i)*f(i)
      end do
   end subroutine radial_integrate_field

   !> Volume quadrature of an analytic integrand: sum_i 4*pi*r_i^2*w_i*f(r_i)
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Radial integrand f(r)
   !> @param[out] result  Quadrature result
   subroutine radial_integrate(self, f, result)
      !> Grid instance
      class(moist_math_grid_radial_type), intent(in) :: self
      !> Radial integrand f(r)
      procedure(integrand_radial) :: f
      !> Quadrature result
      real(wp), intent(out) :: result

      integer :: i

      result = 0.0_wp
      do i = 1, self%npts
         result = result + four_pi*self%r(i)*self%r(i)*self%w(i)*f(self%r(i))
      end do
   end subroutine radial_integrate

end module moist_math_grid_radial_grid
